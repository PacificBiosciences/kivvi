use crate::bam_operation::{start_pos_on_read, ClippedReads};
use crate::realignment::realign::{force_call_d4z4, force_call_kiv2};
use crate::realignment::utilities::Variant;
use crate::repeat_unit::d4z4_variants::update_read_with_special_calls;
use crate::repeat_unit::fingerprint_utils::{
    clean_up_segment_raw_fps, get_good_variants, get_start_end_fps, select_fps, update_fps,
};
use crate::util::{invalid_data_error, missing_data_error, DError, FlankReads, RegionCoordinates};
use log::{debug, trace};
use rust_htslib::bam::ext::BamRecordExtensions;
use rust_htslib::{bam, bam::Read, faidx, htslib};
use std::collections::{BTreeMap, HashSet};
use std::path::PathBuf;

/// Read parameters
#[derive(Clone, Debug)]
pub struct ReadParameters {
    /// minimum base quality
    pub min_base_quality: u8,
    /// minimum variant support
    pub min_variant_support: i32,
    /// minimum fingerprint support
    pub min_fingerprint_support: i32,
    /// correct fingerprints with low support to a more common fingerprint
    pub max_read_count_to_correct: i32,
}

/// Read segment
#[derive(Clone, Debug)]
pub struct ReadSegment {
    /// read name
    name: String,
    /// starting position on read
    pos: i32,
    /// alignment length on read
    aln_len: i32,
    /// bases at variant sites
    fingerprint: Vec<u8>,
    /// fingerprint
    fp_index: Option<i32>,
    /// overlaps left flank
    is_start: bool,
    /// overlaps right flank
    is_end: bool,
}

/// Fingerprint information
#[derive(Clone, Debug)]
pub struct FingerprintInfo {
    /// full read -> vector of fingerprint names
    pub read_edges: BTreeMap<String, Vec<i32>>,
    /// read segment -> fingerprint name
    pub grouped_reads: BTreeMap<String, i32>,
    /// fingerprint seq -> count
    pub fp_count: BTreeMap<Vec<u8>, i32>,
    /// fingerprint name -> fingerprint seq
    pub good_name_to_seq: BTreeMap<i32, Vec<u8>>,
    /// full read -> vector of starting positions on read
    pub read_positions: BTreeMap<String, Vec<i32>>,
    /// read segment name -> (tid, pos) -> base
    pub read_bases: BTreeMap<String, BTreeMap<(i32, i64), Vec<u8>>>,
    /// fingerprint name -> tid
    pub fp_to_tid: BTreeMap<i32, i32>,
    /// retained variants grouped by position
    pub variants_by_position: BTreeMap<i64, Vec<Variant>>,
}

/// Get fingerprints from a realigned bam
/// Fingerprints are named as follows:
/// 0: unknown
/// -1: upstream flanking region
/// -10 < a < -1: starting fingerprints
/// -10: downstream flanking region
/// a < -10: ending fingerprints
/// other positive numbers: good fingerprints
/// # Arguments
/// * `realigned_bam` - path to the realigned bam
/// * `reference` - path to the reference file
/// * `region_coordinates` - information on regions to analyze
/// * `flanking_reads` - reads flanking starts or ends
/// * `read_length` - length of each read
/// * `read_parameters` - read parameters
/// * `read_whitelist` - read whitelist
/// * `is_d4z4` - whether the region is d4z4
/// * `sensitive` - whether to use sensitive mode
/// # Returns
/// * `FingerprintInfo` - fingerprint information
/// * `BTreeMap<String, String>` - base of each read segment at the pivot site
/// * `BTreeMap<String, BTreeMap<usize, usize>>` - cpg sites per read
pub fn get_fingerprint(
    realigned_bam: PathBuf,
    reference: &PathBuf,
    region_coordinates: RegionCoordinates,
    flanking_reads: FlankReads,
    clipped_reads: ClippedReads,
    read_length: &BTreeMap<String, usize>,
    read_parameters: ReadParameters,
    read_whitelist: Vec<String>,
    excluded_count_segments: Vec<String>,
    is_d4z4: bool,
    sensitive: bool,
) -> Result<
    (
        FingerprintInfo,
        BTreeMap<String, String>,
        BTreeMap<String, BTreeMap<usize, usize>>,
    ),
    DError,
> {
    // reference
    let ref_reader = faidx::Reader::from_path(reference)?;
    let ref_names = ref_reader.seq_names()?;
    let excluded_count_segments = excluded_count_segments.into_iter().collect::<HashSet<_>>();

    // read name -> pos -> base
    let mut read_info: BTreeMap<String, BTreeMap<(i32, i64), Vec<u8>>> = BTreeMap::new();
    // filtered variant sites
    let (_filtered_sites, unfiltered_sites, bases_at_pivot_site, cpg_sites_per_read, read_seq_full) =
        get_filtered_sites(
            realigned_bam.clone(),
            reference,
            region_coordinates.clone(),
            &mut read_info,
            read_parameters.clone(),
        )?;

    let mut unfiltered_sites_by_tid: BTreeMap<i32, BTreeMap<i64, Vec<Vec<u8>>>> = BTreeMap::new();
    for ((tid, pos), bases) in unfiltered_sites {
        unfiltered_sites_by_tid
            .entry(tid)
            .or_default()
            .insert(pos, bases);
    }

    let mut read_fps: BTreeMap<String, Vec<ReadSegment>> = BTreeMap::new();
    let mut fp_count_all: BTreeMap<Vec<u8>, i32> = BTreeMap::new();
    let mut good_name_to_seq_all = BTreeMap::new();
    //let mut fp_types = BTreeMap::new();
    let mut fp_to_tid = BTreeMap::new();
    let mut variants_by_position = BTreeMap::new();
    let mut starting_index = 1;
    let num_refs = ref_names.len();
    debug!("num_refs {num_refs}");
    for i in 0..num_refs {
        let ref_name = ref_reader.seq_name(i as i32)?;
        let ref_len = ref_reader.fetch_seq_len(&ref_name);
        let ref_seq = ref_reader.fetch_seq(&ref_name, 0, ref_len as usize)?;
        debug!("reference index {i} has length {ref_len}");
        if unfiltered_sites_by_tid.contains_key(&(i as i32)) {
            let unfiltered_sites = unfiltered_sites_by_tid.get(&(i as i32)).ok_or_else(|| {
                missing_data_error("unfiltered sites for reference tid", i.to_string())
            })?;

            let variant_calls =
                get_good_variants(&region_coordinates, unfiltered_sites, &ref_seq, is_d4z4)?;
            // read segment -> raw fps
            //let mut read_segment_raw_fp =
            //    get_raw_fps_from_alignment(&filtered_sites, &read_info, &region_coordinates)?;
            let mut read_segment_raw_fp: BTreeMap<String, Vec<u8>> = if is_d4z4 {
                force_call_d4z4(
                    &realigned_bam,
                    reference,
                    &variant_calls,
                    &region_coordinates,
                    &read_parameters,
                )?
            } else {
                force_call_kiv2(
                    &realigned_bam,
                    reference,
                    &variant_calls,
                    &region_coordinates,
                    &read_parameters,
                    i,
                )?
            };
            let (cleaned_read_segment_raw_fp, mut new_variants_by_position) =
                clean_up_segment_raw_fps(
                    &mut read_segment_raw_fp,
                    &variant_calls,
                    &flanking_reads,
                    &clipped_reads,
                    &region_coordinates,
                    read_parameters.min_variant_support as usize,
                )?;
            read_segment_raw_fp = cleaned_read_segment_raw_fp;
            if is_d4z4 {
                // for d4z4, genotype two special sites
                (read_segment_raw_fp, new_variants_by_position) = update_read_with_special_calls(
                    &read_segment_raw_fp,
                    &new_variants_by_position,
                    realigned_bam.clone(),
                    reference,
                    &region_coordinates,
                )?;
            }
            if let Some(raw_fp) = read_segment_raw_fp.values().next() {
                if new_variants_by_position.len() != raw_fp.len() {
                    return Err(invalid_data_error(format!(
                        "Variant-position count ({}) does not match fingerprint width ({}) after realignment cleanup for tid {i}",
                        new_variants_by_position.len(),
                        raw_fp.len()
                    )));
                }
            }
            for (pos, variants) in new_variants_by_position {
                variants_by_position
                    .entry(pos)
                    .or_insert_with(Vec::new)
                    .extend(variants);
            }
            // check any starting or ending fingerprints
            let start_end_fps = get_start_end_fps(is_d4z4, &read_segment_raw_fp, &flanking_reads)?;

            // get fingerprints
            let (fp_count, mut good_seq_to_name, mut good_name_to_seq, mut to_replace) =
                select_fps(
                    &read_segment_raw_fp,
                    &read_parameters,
                    &read_whitelist,
                    &excluded_count_segments,
                    &start_end_fps,
                    is_d4z4,
                    sensitive,
                    Some(starting_index),
                )?;

            if is_d4z4 {
                // for d4z4, extra step to rescue fps that are the only prev or next to an existing fp
                // populate read_fps
                // fp info on full reads
                let read_fps = populate_read_fps(
                    &read_segment_raw_fp,
                    &good_seq_to_name,
                    &to_replace,
                    &flanking_reads,
                    is_d4z4,
                )?;
                let mut next_nodes: BTreeMap<i32, Vec<i32>> = BTreeMap::new();
                let mut prev_nodes: BTreeMap<i32, Vec<i32>> = BTreeMap::new();
                let mut next_segs: BTreeMap<i32, HashSet<Vec<u8>>> = BTreeMap::new();
                let mut prev_segs: BTreeMap<i32, HashSet<Vec<u8>>> = BTreeMap::new();
                for (_each_read, mut each_read_info) in read_fps.into_iter() {
                    each_read_info.sort_by(|a, b| a.pos.cmp(&b.pos));
                    let n_segment = each_read_info.len();
                    for segment_index in 0..(n_segment - 1) {
                        let seg1 = &each_read_info[segment_index];
                        let seg2 = &each_read_info[segment_index + 1];
                        if !seg1.fingerprint.contains(&b'x') && !seg2.fingerprint.contains(&b'x') {
                            let seg1_fp = seg1.fp_index;
                            let seg2_fp = seg2.fp_index;
                            if let Some(seg1_fp) = seg1_fp {
                                next_segs
                                    .entry(seg1_fp)
                                    .or_default()
                                    .insert(seg2.fingerprint.clone());
                                if let Some(seg2_fp) = seg2_fp {
                                    next_nodes.entry(seg1_fp).or_default().push(seg2_fp);
                                }
                            }
                            if let Some(seg2_fp) = seg2_fp {
                                prev_segs
                                    .entry(seg2_fp)
                                    .or_default()
                                    .insert(seg1.fingerprint.clone());
                                if let Some(seg1_fp) = seg1_fp {
                                    prev_nodes.entry(seg2_fp).or_default().push(seg1_fp);
                                }
                            }
                        }
                    }
                }
                let mut fps_to_add = HashSet::new();
                for (node, segs_next) in next_segs {
                    if !next_nodes.contains_key(&node) {
                        if segs_next.len() == 1 {
                            let next_seg_vec = segs_next.into_iter().collect::<Vec<Vec<u8>>>();
                            let next_seg_seq = next_seg_vec.first().ok_or_else(|| {
                                format!("Expected one successor segment for fingerprint {node}")
                            })?;
                            fps_to_add.insert(next_seg_seq.to_vec());
                            debug!(
                                "for fp {node}, the only next segment is {:?}",
                                std::str::from_utf8(next_seg_seq)?
                            );
                        }
                    }
                }
                for (node, segs_prev) in prev_segs {
                    if !prev_nodes.contains_key(&node) {
                        if segs_prev.len() == 1 {
                            let prev_seg_vec = segs_prev.into_iter().collect::<Vec<Vec<u8>>>();
                            let prev_seg_seq = prev_seg_vec.first().ok_or_else(|| {
                                format!("Expected one predecessor segment for fingerprint {node}")
                            })?;
                            fps_to_add.insert(prev_seg_seq.to_vec());
                            debug!(
                                "for fp {node}, the only prev segment is {:?}",
                                std::str::from_utf8(prev_seg_seq)?
                            );
                        }
                    }
                }
                // update fps
                to_replace = update_fps(
                    &fp_count,
                    &mut good_seq_to_name,
                    &mut good_name_to_seq,
                    fps_to_add,
                )?;
            }

            starting_index += good_seq_to_name.len() as i32;
            for (k, v) in &fp_count {
                fp_count_all.insert(k.clone(), *v);
            }
            for (k, v) in &good_name_to_seq {
                good_name_to_seq_all.insert(*k, v.clone());
                fp_to_tid.insert(*k, i as i32);
            }

            // populate read_fps
            // fp info on full reads
            let this_tid_read_fps: BTreeMap<String, Vec<ReadSegment>> = populate_read_fps(
                &read_segment_raw_fp,
                &good_seq_to_name,
                &to_replace,
                &flanking_reads,
                is_d4z4,
            )?;
            for (read, segs) in this_tid_read_fps {
                for seg in segs {
                    read_fps.entry(read.clone()).or_default().push(seg);
                }
            }
        }
    }
    let unknown_full_fp_index = 0;

    // get read edges
    let mut read_edges = BTreeMap::new();
    let mut grouped_reads = BTreeMap::new();
    // full read names -> vector of starting positions of each segment
    let mut read_positions: BTreeMap<String, Vec<i32>> = BTreeMap::new();
    for (each_read, mut each_read_info) in read_fps.into_iter() {
        let this_read_length = *read_length
            .get(&each_read)
            .ok_or_else(|| missing_data_error("read length", &each_read))?
            as i32;
        each_read_info.sort_by(|a, b| a.pos.cmp(&b.pos));

        let mut this_read_edges: Vec<i32> = Vec::new();
        let mut this_read_positions: Vec<i32> = Vec::new();
        let first_seg = each_read_info.first().ok_or_else(|| {
            format!("Read '{each_read}' has no aligned segments after fingerprinting")
        })?;
        // starts
        if first_seg.is_start && first_seg.pos > 1000 {
            this_read_edges.push(-1);
            this_read_positions.push(first_seg.pos);
        }
        let mut prev_position = first_seg.pos;
        let mut prev_aln = first_seg.aln_len;
        let mut prev_fp_is_end = false;

        for (segment_index, each_read_segment) in each_read_info.clone().into_iter().enumerate() {
            let current_position = each_read_segment.pos;
            let segment_name = each_read_segment.name.clone();
            let name_fields = segment_name.split_terminator(':').collect::<Vec<_>>();
            let full_read_name = name_fields[0].to_string();
            let pos_on_read = name_fields[1].parse::<i32>()?;
            let segment_name_short = format!("{full_read_name}:{pos_on_read}");
            let fp_seq_string = std::str::from_utf8(&each_read_segment.fingerprint)?;
            debug!(
                "read {} {} {} {:?} {}",
                each_read,
                segment_name,
                each_read_segment.pos,
                each_read_segment.fp_index,
                fp_seq_string
            );
            if prev_position > 20 && !is_d4z4 {
                prev_aln = region_coordinates.repeat_len as i32;
            }
            if segment_index > 0 {
                let expected_current_position = prev_position + prev_aln;
                if !prev_fp_is_end
                    && region_coordinates.repeat_len > 800
                    && (current_position > expected_current_position + 1000
                        || current_position < expected_current_position - 500)
                {
                    // replace a suspicious unit with unknown
                    this_read_edges.pop();
                    this_read_edges.push(0);
                    debug!(
                        "read segment {:?} is suspicious: current position {} expected {}",
                        segment_name, current_position, expected_current_position
                    );
                }
                // if one copy is missing. This step is for repeats with shorter units, <2kb
                if !prev_fp_is_end && region_coordinates.repeat_len < 2000 {
                    let position_offset = current_position - expected_current_position;
                    if position_offset > region_coordinates.repeat_len as i32 - 100 {
                        // replace a suspicious unit with unknown
                        let number_of_unknown = ((position_offset + 100) as f64
                            / (region_coordinates.repeat_len as f64))
                            .round() as usize;
                        for k in 0..number_of_unknown {
                            this_read_edges.push(unknown_full_fp_index);
                            this_read_positions.push(
                                expected_current_position
                                    + (k * region_coordinates.repeat_len) as i32,
                            );

                            let n1 = expected_current_position as usize
                                + (k * region_coordinates.repeat_len);
                            let n2 = n1 + region_coordinates.repeat_len;
                            let unknown_fp_seq = std::str::from_utf8(
                                &read_seq_full
                                    .get(&segment_name_short)
                                    .ok_or_else(|| {
                                        format!(
                                            "Missing full read sequence for segment '{segment_name_short}' while inferring unknown fingerprints"
                                        )
                                    })?[n1..n2],
                            )?;
                            debug!(">{each_read}_unknown_fp_{n1} {unknown_fp_seq}");
                        }
                        debug!(
                            "read segment {:?} has unknown fp: current position {} expected {}. number_of_unknown {number_of_unknown}",
                            segment_name, current_position, expected_current_position
                        );
                    }
                }
            }
            this_read_positions.push(current_position);
            if each_read_segment.fp_index.is_none() {
                this_read_edges.push(0);
                grouped_reads.entry(segment_name_short).or_insert(0);
            } else {
                let fp_name = each_read_segment.fp_index.ok_or_else(|| {
                    format!(
                        "Read segment '{}' is missing a fingerprint index after classification",
                        each_read_segment.name
                    )
                })?;
                this_read_edges.push(fp_name);
                grouped_reads.entry(segment_name_short).or_insert(fp_name);
            }
            if each_read_segment.is_end && is_d4z4 {
                this_read_edges.push(-10);
                this_read_positions.push(current_position);
                prev_fp_is_end = true;
            } else {
                prev_fp_is_end = false;
            }
            prev_position = each_read_segment.pos;
            prev_aln = each_read_segment.aln_len;
        }
        if !is_d4z4 {
            let last_seg = each_read_info
                .last()
                .ok_or_else(|| missing_data_error("last read segment", &each_read))?;
            let last_fp_seq = std::str::from_utf8(&last_seg.fingerprint)?;
            if last_fp_seq.starts_with("xxxx") {
                // replace a suspicious unit with unknown
                // not for cases with D4Z4 qA ends
                this_read_edges.pop();
                this_read_edges.push(0);
                debug!(
                    "last read segment {:?} is suspicious {:?}, replace with unknown",
                    last_seg, last_fp_seq
                );
            }
            // ends
            if last_seg.is_end {
                if last_seg.pos + last_seg.aln_len < this_read_length - 1000 {
                    this_read_edges.push(-10);
                    this_read_positions.push(last_seg.pos);
                }
            }
        }
        let this_read_edges_string = this_read_edges
            .iter()
            .map(|x| x.to_string())
            .collect::<Vec<String>>()
            .join("-");
        debug!("read {each_read} length {this_read_length:?} {this_read_edges_string:?}");
        debug!("read {each_read} length {this_read_length:?} {this_read_positions:?}");
        if this_read_edges.len() != this_read_positions.len() {
            return Err(invalid_data_error(format!(
                "Read '{each_read}' produced {} fingerprint edges but {} fingerprint positions after segmentation",
                this_read_edges.len(),
                this_read_positions.len()
            )));
        }
        read_edges
            .entry(each_read.clone())
            .or_insert(this_read_edges);
        read_positions
            .entry(each_read.clone())
            .or_insert(this_read_positions);
    }

    // simplified read info
    // read name (:pos) -> pos -> base
    let mut read_info_simple: BTreeMap<String, BTreeMap<(i32, i64), Vec<u8>>> = BTreeMap::new();
    for (read, per_read_info) in read_info.iter() {
        let mut read_name_split = read
            .to_string()
            .split_terminator(':')
            .map(|x| x.to_string())
            .collect::<Vec<_>>();
        read_name_split.pop();
        let read_name_simple = read_name_split.join(":");
        read_info_simple
            .entry(read_name_simple)
            .or_insert(per_read_info.clone());
    }

    Ok((
        FingerprintInfo {
            read_edges,
            grouped_reads,
            fp_count: fp_count_all,
            good_name_to_seq: good_name_to_seq_all,
            read_positions,
            read_bases: read_info_simple,
            fp_to_tid,
            variants_by_position,
        },
        bases_at_pivot_site,
        cpg_sites_per_read,
    ))
}

/// Get filtered sites from bam.
/// Filtering is based on base counts and prior knowledge
/// # Arguments
/// * `realigned_bam` - realigned bam
/// * `reference` - reference file path
/// * `region_coordinates` - coordinates defined in this region
/// * `aln2seq` - store read sequences
/// * `read_info` - read name -> pos -> base; to be filled in
/// * `read_parameters` - read parameters
/// # Returns
/// * `BTreeMap<(i32, i64), Vec<Vec<u8>>>` - filtered sites
/// * `BTreeMap<(i32, i64), Vec<Vec<u8>>>` - unfiltered sites
/// * `BTreeMap<String, String>` - base of each read segment at the pivot site
/// * `BTreeMap<String, BTreeMap<usize, usize>>` - cpg sites per read
/// * `BTreeMap<String, Vec<u8>>` - full read sequences
fn get_filtered_sites(
    realigned_bam: PathBuf,
    reference: &PathBuf,
    region_coordinates: RegionCoordinates,
    read_info: &mut BTreeMap<String, BTreeMap<(i32, i64), Vec<u8>>>,
    read_parameters: ReadParameters,
) -> Result<
    (
        BTreeMap<(i32, i64), Vec<Vec<u8>>>,
        BTreeMap<(i32, i64), Vec<Vec<u8>>>,
        BTreeMap<String, String>,
        BTreeMap<String, BTreeMap<usize, usize>>,
        BTreeMap<String, Vec<u8>>,
    ),
    DError,
> {
    // store read sequences
    let mut aln2seq = BTreeMap::new();
    let mut aln2seq_full = BTreeMap::new();
    let mut bam_reader = bam::IndexedReader::from_path(realigned_bam.clone())?;
    let ref_reader = faidx::Reader::from_path(reference)?;
    let clip_variant_sites = region_coordinates.clip_variant_sites;
    // bases across each position, pos -> (base, count)
    let mut raw_piles = BTreeMap::new();
    for x in bam_reader.pileup() {
        let x = x?;
        let pos = i64::from(x.pos());
        let tid = x.tid() as i32;
        let ref_name = ref_reader.seq_name(tid)?;
        let ref_len = ref_reader.fetch_seq_len(&ref_name);
        let ref_seq = ref_reader.fetch_seq(&ref_name, 0, ref_len as usize)?;
        let base_counts = query_seq_counter(
            &x,
            &ref_seq,
            &mut aln2seq,
            &mut aln2seq_full,
            read_info,
            read_parameters.min_base_quality,
        )?;
        raw_piles.insert((tid, pos), base_counts);
    }
    let mut bam_reader = bam::IndexedReader::from_path(realigned_bam.clone())?;
    for x in bam_reader.pileup() {
        let x = x?;
        let pos = i64::from(x.pos());
        let adj_pos = pos + 2;
        let variant_pos = pos + 1;
        if clip_variant_sites.contains_key(&adj_pos) {
            //trace!("position to check clip variants: {}, {}", pos, adj_pos);
            let next_base_counts = query_seq_counter_at_clip_site(
                &x,
                aln2seq_full.clone(),
                read_info,
                clip_variant_sites.clone(),
            )?;

            if !raw_piles.contains_key(&(0, variant_pos)) {
                raw_piles.insert((0, variant_pos), next_base_counts);
            } else {
                for (base, count) in next_base_counts {
                    debug!(
                        "next_base_counts pos {} {} {}",
                        pos + 1,
                        std::str::from_utf8(&base)?,
                        count
                    );
                    let entry = raw_piles.entry((0, variant_pos)).or_insert(BTreeMap::new());
                    if !entry.contains_key(&base) {
                        entry.insert(base, count);
                    } else {
                        *entry.entry(base).or_insert(0) += count;
                    }
                }
            }
        }
    }

    // pivot site
    let mut bases_at_special_site: BTreeMap<String, String> = BTreeMap::new();
    if let Some(pivot_site) = region_coordinates.pivot_site {
        bases_at_special_site =
            check_pivot_site(realigned_bam.clone(), aln2seq.clone(), pivot_site)?;
    }
    // cpg sites
    let mut cpg_sites_per_read: BTreeMap<String, BTreeMap<usize, usize>> = BTreeMap::new();
    let methyl_sites = region_coordinates.methyl_sites;
    if !methyl_sites.is_empty() {
        cpg_sites_per_read = get_cpg_sites(
            realigned_bam.clone(),
            aln2seq.clone(),
            aln2seq_full.clone(),
            methyl_sites,
        )?;
    }

    let mut filtered_sites: BTreeMap<(i32, i64), Vec<Vec<u8>>> = BTreeMap::new();
    let mut unfiltered_sites: BTreeMap<(i32, i64), Vec<Vec<u8>>> = BTreeMap::new();
    for ((tid, pos), base_counts) in raw_piles.into_iter() {
        // if multiple reference, then use all positions
        // else, use only positions within the specified repeat length
        if !region_coordinates.exclude_sites.contains(&(pos + 1))
            && (pos + 1 < region_coordinates.repeat_len as i64
                || region_coordinates.genome_offset.len() > 1)
        {
            let mut filtered_bases_per_site = Vec::new();
            let mut all_including_indels_per_site = Vec::new();
            let mut keep_position = false;
            let mut indel_count = 0;
            let total_depth = base_counts.clone().into_values().sum::<i32>();
            let mut indel_thres = (total_depth as f64) * 0.025;

            if indel_thres < 9.0 {
                indel_thres = 9.0;
            }

            for (base, count) in base_counts.clone().into_iter() {
                if base.len() > 1 || base.contains(&b'*') {
                    indel_count += count;
                }
            }
            trace!("pos {pos:?} has bases {base_counts:?} indel_count {indel_count:?} indel_thres {indel_thres:?} total_depth {total_depth:?}");
            for (base, count) in base_counts.clone().into_iter() {
                let base_str = std::str::from_utf8(&base)?;
                trace!("pos {pos} has base1 {base_str} at count {count}");
                if !base.contains(&b'+')
                    && !base.contains(&b'-')
                    && !base.contains(&b'*')
                    && count >= read_parameters.min_variant_support
                {
                    filtered_bases_per_site.push(base.clone());
                    all_including_indels_per_site.push(base);
                } else if (base.contains(&b'+') || base.contains(&b'-'))
                    && count >= read_parameters.min_variant_support
                {
                    all_including_indels_per_site.push(base);
                }
            }

            if filtered_bases_per_site.len() > 1
                && indel_count <= 20
                && (indel_count as f64) < indel_thres
            {
                keep_position = true;
            }
            if keep_position {
                filtered_sites
                    .entry((tid, pos))
                    .or_insert(filtered_bases_per_site.clone());
            }
            unfiltered_sites
                .entry((tid, pos))
                .or_insert(all_including_indels_per_site.clone());
        }
    }

    Ok((
        filtered_sites,
        unfiltered_sites,
        bases_at_special_site,
        cpg_sites_per_read,
        aln2seq,
    ))
}

/// populate read_fps
/// # Arguments
/// * `read_segment_raw_fp` - read segment -> raw fps
/// * `good_seq_to_name` - fingerprint sequence to name
/// * `to_replace` - fingerprints to replace
/// * `flanking_reads` - flanking reads
/// * `is_d4z4` - whether the region is d4z4
/// # Returns
/// * `BTreeMap<String, Vec<ReadSegment>>` - read name -> read segments
fn populate_read_fps(
    read_segment_raw_fp: &BTreeMap<String, Vec<u8>>,
    good_seq_to_name: &BTreeMap<Vec<u8>, i32>,
    to_replace: &BTreeMap<Vec<u8>, i32>,
    flanking_reads: &FlankReads,
    is_d4z4: bool,
) -> Result<BTreeMap<String, Vec<ReadSegment>>, DError> {
    // fp info on full reads
    let mut read_fps: BTreeMap<String, Vec<ReadSegment>> = BTreeMap::new();
    for (read_segment, raw_fp) in read_segment_raw_fp {
        let name_fields = read_segment.split_terminator(':').collect::<Vec<_>>();
        let full_read_name = name_fields[0].to_string();
        let pos_on_read = name_fields[1].parse::<i32>()?;
        let seg_aln_len = name_fields[2].parse::<i32>()?;

        let mut final_fp: Option<i32> = good_seq_to_name.get(raw_fp).copied();
        if final_fp.is_none() {
            final_fp = to_replace.get(raw_fp).copied();
        }

        let read_is_start = flanking_reads.start.contains(&full_read_name);
        let read_is_end = if !is_d4z4 {
            flanking_reads.end.contains(&full_read_name)
        } else {
            flanking_reads.end.contains(read_segment)
        };
        let this_segment = ReadSegment {
            name: read_segment.to_string(),
            pos: pos_on_read,
            aln_len: seg_aln_len,
            fingerprint: raw_fp.to_vec(),
            fp_index: final_fp,
            is_start: read_is_start,
            is_end: read_is_end,
        };
        read_fps
            .entry(full_read_name)
            .or_default()
            .push(this_segment);
    }
    Ok(read_fps)
}

/// Return the raw position of the current alignment on the read
/// # Arguments
/// * `x` - reference to an alignment
/// # Returns
/// * `usize` - raw position of the current alignment on the read
fn raw_qpos<'a>(x: &'a bam::pileup::Alignment<'a>) -> usize {
    static_assertions::const_assert_eq!(
        std::mem::size_of::<&bam::pileup::Alignment<'_>>(),
        std::mem::size_of::<&htslib::bam_pileup1_t>()
    );
    let ptr = (x as *const bam::pileup::Alignment<'_>).cast::<&htslib::bam_pileup1_t>();
    unsafe { *ptr }.qpos as usize
}

/// Given a pileup at one position, return bases and their counts
/// By Daniel Baker
/// # Arguments
/// * `x` - a Pileup to work with.
/// * `ref_seq` - reference sequence.
/// * `aln2seq` - store read sequences.
/// * `read_info` - read name -> pos -> base; to be filled in
/// * `min_base_quality` - minimum base quality
/// # Returns
/// * `BTreeMap<Vec<u8>, i32>` - base counts
fn query_seq_counter(
    x: &bam::pileup::Pileup,
    ref_seq: &[u8],
    aln2seq: &mut BTreeMap<String, Vec<u8>>,
    _aln2seq_full: &mut BTreeMap<String, Vec<u8>>,
    read_info: &mut BTreeMap<String, BTreeMap<(i32, i64), Vec<u8>>>,
    min_base_quality: u8,
) -> Result<BTreeMap<Vec<u8>, i32>, DError> {
    // bases on reads
    let mut base_counts = BTreeMap::<Vec<u8>, i32>::new();
    let tid = x.tid();
    let pos = x.pos();
    for aln in x.alignments() {
        let query_pos_raw = raw_qpos(&aln);
        let record = aln.record();
        let qname = std::str::from_utf8(record.qname())?;
        let read_start_pos = start_pos_on_read(&record);
        let reference_start_pos = &record.pos();
        let reference_end_pos = &record.reference_end();
        let aln_len = reference_end_pos - reference_start_pos;
        let is_reverse = record.is_reverse();
        let bq = record.qual().get(query_pos_raw).copied().unwrap_or(0);
        // Cache query seq, as this is can be shared across thousands of sites.
        let entry = aln2seq
            .entry(format!("{qname}:{}", read_start_pos))
            .or_insert_with(|| record.seq().as_bytes());
        /*
        // no need now as all seqs have cigar S instead of H
        if !record.is_supplementary() && !record.is_secondary() && !aln2seq_full.contains_key(qname)
        {
            aln2seq_full
                .entry(qname.to_string())
                .or_insert(record.seq().as_bytes());
        }
        */
        let seq: &[u8] = &entry[..];
        let mut query_seq = Vec::new();

        let ambiguous_bases = vec![b'-', b'N', b'<', b'>', b'*'];
        let base = if bq < min_base_quality {
            b'-'
        } else if !aln.is_del() && !aln.is_refskip() {
            seq.get(query_pos_raw).copied().unwrap_or(b'N')
        } else if aln.is_refskip() {
            if is_reverse {
                b'<'
            } else {
                b'>'
            }
        } else {
            b'*'
        };
        query_seq.push(base);
        let pos = pos as usize;
        if !ambiguous_bases.contains(&base) {
            match aln.indel() {
                bam::pileup::Indel::Ins(x) => {
                    debug_assert!(x > 0);
                    query_seq.push(b'+');
                    query_seq.extend_from_slice(x.to_string().as_bytes());
                    // TODO: speed this up by using slice operations to copy out faster.
                    for j in 1..=(x as usize) {
                        query_seq.push(seq[j + query_pos_raw]);
                    }
                }
                bam::pileup::Indel::Del(x) => {
                    debug_assert!(x > 0);
                    query_seq.push(b'-');
                    query_seq.extend_from_slice(x.to_string().as_bytes());
                    for j in 1..=(x as usize) {
                        query_seq.push(ref_seq[j + pos]);
                    }
                }
                bam::pileup::Indel::None => {}
            }
        }
        //query_seq.make_ascii_uppercase();
        // use high base quality (13) for base counts
        if bq >= 13 {
            *base_counts.entry(query_seq.clone()).or_default() += 1;
        }

        let per_base: BTreeMap<(i32, i64), Vec<u8>> = BTreeMap::new();
        let query_seq_vec = query_seq.clone().to_vec();

        read_info
            .entry(format!("{qname}:{}:{}", read_start_pos, aln_len))
            .or_insert_with(|| per_base)
            .entry((tid as i32, pos as i64))
            .or_insert_with(|| query_seq_vec.clone());
    }
    Ok(base_counts)
}

/// Get the base of each read segment at the clip site
/// # Arguments
/// * `x` - a Pileup to work with.
/// * `_aln2seq_full` - full read sequence
/// * `read_info` - read name -> pos -> base; to be filled in
/// * `clip_variant_sites` - clip variant sites
/// # Returns
/// * `BTreeMap<Vec<u8>, i32>` - base counts
fn query_seq_counter_at_clip_site(
    x: &bam::pileup::Pileup,
    _aln2seq_full: BTreeMap<String, Vec<u8>>,
    read_info: &mut BTreeMap<String, BTreeMap<(i32, i64), Vec<u8>>>,
    clip_variant_sites: BTreeMap<i64, Vec<u8>>,
) -> Result<BTreeMap<Vec<u8>, i32>, DError> {
    let tid = x.tid() as i32;
    let mut base_counts_next = BTreeMap::<Vec<u8>, i32>::new();
    let pos = x.pos();
    let this_pos = pos as i64;
    let variant_pos = this_pos + 1;
    let adj_pos = this_pos + 2;
    let expected_base = clip_variant_sites
        .get(&adj_pos)
        .ok_or("adj_pos not in clip_variant_sites")?;
    for aln in x.alignments() {
        let query_pos_raw = raw_qpos(&aln);
        let record = aln.record();
        let qname = std::str::from_utf8(record.qname())?;
        let read_start_pos = start_pos_on_read(&record);
        let reference_start_pos = &record.pos();
        let reference_end_pos = &record.reference_end();
        let aln_len = reference_end_pos - reference_start_pos;
        // alignment ends here
        if *reference_end_pos == this_pos + 1 {
            let query_pos_raw_next = query_pos_raw + 1;
            let full_seq = record.seq().as_bytes();
            /*
            // only when using Hardclips
            if record.is_supplementary() {
                query_pos_raw_next += read_start_pos as usize;
            }
            let mut full_seq: Vec<u8> = Vec::new();
            if aln2seq_full.contains_key(qname) {
                let this_read_full_seq = aln2seq_full
                    .get(qname)
                    .ok_or("qname not found in aln2seq_full")?;
                full_seq = this_read_full_seq.to_vec();
            }
            */
            //trace!("{qname} full_seq {:?}", std::str::from_utf8(&full_seq)?);
            if query_pos_raw_next < full_seq.len() {
                let next_base = full_seq[query_pos_raw_next];
                if expected_base.contains(&next_base) {
                    let next_base_vec = vec![next_base];
                    trace!(
                        "{qname}:{}:{} ref position {} read position {} base {:?}",
                        read_start_pos,
                        aln_len,
                        this_pos,
                        query_pos_raw_next - 1,
                        std::str::from_utf8(&vec![full_seq[query_pos_raw_next - 1]])?,
                    );
                    trace!(
                        "{qname}:{}:{} ref position {} read position {} base {:?}",
                        read_start_pos,
                        aln_len,
                        variant_pos,
                        query_pos_raw_next,
                        std::str::from_utf8(&vec![full_seq[query_pos_raw_next]])?,
                    );
                    *base_counts_next.entry(next_base_vec.clone()).or_default() += 1;
                    read_info
                        .entry(format!("{qname}:{}:{}", read_start_pos, aln_len))
                        .or_insert_with(|| BTreeMap::new())
                        .entry((tid, variant_pos))
                        .or_insert_with(|| next_base_vec.clone());
                }
            }
        }
    }
    Ok(base_counts_next)
}

/// Get the base of each read segment at the pivot site
/// # Arguments
/// * `realigned_bam` - realigned bam
/// * `aln2seq` - read sequence
/// * `pivot_site` - pivot site
/// # Returns
/// * `BTreeMap<String, String>` - read segment -> base
pub fn check_pivot_site(
    realigned_bam: PathBuf,
    aln2seq: BTreeMap<String, Vec<u8>>,
    pivot_site: i64,
) -> Result<BTreeMap<String, String>, DError> {
    let mut bases_at_special_site = BTreeMap::new();
    // do not use reads whose previous base is an indel
    let mut read_blacklist = Vec::new();
    let mut bam_reader = bam::IndexedReader::from_path(realigned_bam.clone())?;
    for x in bam_reader.pileup() {
        let x = x?;
        let pos = i64::from(x.pos());
        if pos + 2 == pivot_site {
            for aln in x.alignments() {
                let record = aln.record();
                let read_start_pos = start_pos_on_read(&record);
                let qname = std::str::from_utf8(record.qname())?;
                let qname_new = format!("{qname}:{}", read_start_pos);
                if aln.indel() != rust_htslib::bam::pileup::Indel::None || aln.is_del() {
                    read_blacklist.push(qname_new);
                }
            }
        }
        if pos + 1 == pivot_site {
            for aln in x.alignments() {
                let query_pos_raw = raw_qpos(&aln);
                let record = aln.record();
                let read_start_pos = start_pos_on_read(&record);
                let qname = std::str::from_utf8(record.qname())?;
                let qname_new = format!("{qname}:{}", read_start_pos);
                let seq = aln2seq.get(&qname_new).ok_or("qname_new not in aln2seq")?;
                if seq.len() > query_pos_raw + 6 {
                    let base = std::str::from_utf8(&seq[query_pos_raw..(query_pos_raw + 6)])?;
                    bases_at_special_site.insert(qname_new.to_string(), base.to_string());
                }
            }
        }
    }
    Ok(bases_at_special_site)
}

/// Get all CpG sites in reads
/// Returns read segment -> position on ref  -> index of C
/// # Arguments
/// * `realigned_bam` - realigned bam
/// * `aln2seq` - read sequence
/// * `_aln2seq_full` - full read sequence
/// * `methyl_sites` - methyl sites
/// # Returns
/// * `BTreeMap<String, BTreeMap<usize, usize>>` - read segment -> position on ref -> index of C
pub fn get_cpg_sites(
    realigned_bam: PathBuf,
    aln2seq: BTreeMap<String, Vec<u8>>,
    _aln2seq_full: BTreeMap<String, Vec<u8>>,
    methyl_sites: Vec<usize>,
) -> Result<BTreeMap<String, BTreeMap<usize, usize>>, DError> {
    let mut cpg_sites_per_read: BTreeMap<String, BTreeMap<usize, usize>> = BTreeMap::new();
    let mut bam_reader = bam::IndexedReader::from_path(realigned_bam.clone())?;
    for x in bam_reader.pileup() {
        let x = x?;
        let pos = x.pos() as usize;
        let pos_adjust = pos + 1;
        if methyl_sites.contains(&pos_adjust) {
            for aln in x.alignments() {
                let query_pos_raw = raw_qpos(&aln);
                let record = aln.record();
                let read_start_pos = start_pos_on_read(&record);
                let qname = std::str::from_utf8(record.qname())?;
                let qname_new = format!("{qname}:{}", read_start_pos);
                if aln2seq.contains_key(&qname_new) {
                    // && aln2seq_full.contains_key(qname) {
                    let seq = aln2seq.get(&qname_new).ok_or_else(|| {
                        format!(
                            "Missing aligned read sequence for segment '{qname_new}' while collecting CpG sites"
                        )
                    })?;
                    // seq is full_seq
                    //let full_seq = aln2seq_full.get(qname).unwrap();
                    if !aln.is_del() && query_pos_raw + 1 < seq.len() {
                        let base = seq.get(query_pos_raw).copied().unwrap_or(b'N');
                        let base_next = seq.get(query_pos_raw + 1).copied().unwrap_or(b'N');
                        if base == b'C' && base_next == b'G' {
                            //trace!(
                            //"{qname} position {pos} read_start_pos {read_start_pos} supplementary {}",
                            //record.is_supplementary()
                            //);
                            // only when using Hardclips
                            //if record.is_supplementary() {
                            //    query_pos_raw += read_start_pos as usize;
                            //}
                            if !record.is_reverse() {
                                let current_c_site = seq[..query_pos_raw]
                                    .windows(2)
                                    .filter(|&x| x == "CG".as_bytes())
                                    .count();
                                //trace!("current_c_site {current_c_site}");
                                cpg_sites_per_read
                                    .entry(qname_new.to_string())
                                    .or_default()
                                    .insert(pos, current_c_site);
                            } else {
                                let current_c_site = seq[query_pos_raw..]
                                    .windows(2)
                                    .filter(|&x| x == "CG".as_bytes())
                                    .count()
                                    - 1;
                                //trace!("current_c_site {current_c_site}");
                                cpg_sites_per_read
                                    .entry(qname_new.to_string())
                                    .or_default()
                                    .insert(pos, current_c_site);
                            }
                        }
                    }
                }
            }
        }
    }
    Ok(cpg_sites_per_read)
}
