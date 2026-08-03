use crate::bam_operation::start_pos_on_read;
use crate::realignment::utilities::{AlleleType, Variant, VariantType};
use crate::realignment::wfa_graph::{NodeAlleleMap, WFAGraph, WFAResult};
use crate::repeat_unit::fingerprint::ReadParameters;
use crate::util::RegionCoordinates;
use log::{debug, trace};
use paraphase::detail::util::DError;
use rust_htslib::bam;
use rust_htslib::bam::ext::BamRecordExtensions;
use rust_htslib::bam::Read;
use rust_htslib::faidx;
use std::collections::{BTreeMap, HashMap};
use std::path::PathBuf;

/// Force call KIV2 variants
/// # Arguments
/// * `realigned_bam` - realigned bam
/// * `reference` - reference file path
/// * `variant_calls` - variant calls
/// * `region_coordinates` - region coordinates
/// * `read_parameters` - read parameters
/// * `ref_index` - reference index
/// # Returns
/// * `BTreeMap<String, Vec<u8>>` - read segment -> raw fps
pub fn force_call_kiv2(
    realigned_bam: &PathBuf,
    reference: &PathBuf,
    variant_calls: &[Variant],
    region_coordinates: &RegionCoordinates,
    read_parameters: &ReadParameters,
    ref_index: usize,
) -> Result<BTreeMap<String, Vec<u8>>, DError> {
    // read segment -> raw fps
    let mut read_segment_raw_fp: BTreeMap<String, Vec<u8>> = BTreeMap::new();
    // reference
    let ref_reader = faidx::Reader::from_path(reference)?;
    let ref_name = ref_reader.seq_name(ref_index as i32)?;
    let ref_len = ref_reader.fetch_seq_len(&ref_name);
    let ref_seq = ref_reader.fetch_seq(&ref_name, 0, ref_len as usize)?;
    let mut bam_reader = bam::IndexedReader::from_path(realigned_bam)?;
    bam_reader
        .fetch((&ref_name, 0, (region_coordinates.repeat_len as i64)))
        .map_err(|e| {
            std::io::Error::other(format!(
                "Failed to fetch KIV2 force-call region 0-{} on {ref_name}: {e}",
                region_coordinates.repeat_len
            ))
        })?;
    for read_entry in bam_reader.records() {
        let mut read = read_entry?;
        //build out the cigar info
        read.cache_cigar();
        let qname = std::str::from_utf8(read.qname())?;
        let read_start_pos = start_pos_on_read(&read);
        let reference_start_pos = &read.pos();
        let reference_end_pos = &read.reference_end();
        let aln_len = reference_end_pos - reference_start_pos;
        let segment_name = format!("{qname}:{}:{}", read_start_pos, aln_len);

        let read_seq = global_realignment(
            &read,
            &variant_calls,
            &vec![],
            &ref_seq,
            0,
            ref_len as i64,
            read_parameters.min_base_quality,
            &region_coordinates.clip_variant_sites,
        )?;
        debug!(
            "segment_name {segment_name} read_simplified {:?}",
            std::str::from_utf8(&read_seq)?
        );
        read_segment_raw_fp.insert(segment_name, read_seq);
    }

    Ok(read_segment_raw_fp)
}

/// Realign all reads to force call pre-selected variants for D4Z4
/// # Arguments
/// * `realigned_bam` - realigned bam
/// * `reference` - reference file
/// * `variant_calls` - variant calls
/// * `region_coordinates` - region coordinates
/// * `read_parameters` - read parameters
/// # Returns
/// * `Result<BTreeMap<String, Vec<u8>>, DError>` - read segment name -> calls at variant sites
pub fn force_call_d4z4(
    realigned_bam: &PathBuf,
    reference: &PathBuf,
    variant_calls: &[Variant],
    region_coordinates: &RegionCoordinates,
    read_parameters: &ReadParameters,
) -> Result<BTreeMap<String, Vec<u8>>, DError> {
    let realign_segments = &region_coordinates.realign_segments;
    // read segment -> raw fps
    let mut read_segment_raw_fp: BTreeMap<String, Vec<u8>> = BTreeMap::new();
    // reference
    let ref_reader = faidx::Reader::from_path(reference)?;
    let ref_name = ref_reader.seq_name(0)?;
    let ref_len = ref_reader.fetch_seq_len(&ref_name);
    let ref_seq = ref_reader.fetch_seq(&ref_name, 0, ref_len as usize)?;
    let mut bam_reader = bam::IndexedReader::from_path(realigned_bam)?;
    bam_reader
        .fetch((&ref_name, 0, (region_coordinates.repeat_len as i64)))
        .map_err(|e| {
            std::io::Error::other(format!(
                "Failed to fetch D4Z4 force-call region 0-{} on {ref_name}: {e}",
                region_coordinates.repeat_len
            ))
        })?;
    for read_entry in bam_reader.records() {
        let mut read = read_entry?;
        //build out the cigar info
        read.cache_cigar();
        let qname = std::str::from_utf8(read.qname())?;
        let read_start_pos = start_pos_on_read(&read);
        let reference_start_pos = &read.pos();
        let reference_end_pos = &read.reference_end();
        let aln_len = reference_end_pos - reference_start_pos;
        let segment_name = format!("{qname}:{}:{}", read_start_pos, aln_len);
        let variant_calls1 = variant_calls
            .iter()
            .filter(|x| x.position() + 1 < realign_segments[0])
            .map(|x| x.clone())
            .collect::<Vec<_>>();
        let mut read_seq = global_realignment(
            &read,
            &variant_calls1,
            &vec![],
            &ref_seq,
            0,
            realign_segments[0],
            read_parameters.min_base_quality,
            &BTreeMap::new(),
        )?;

        let variant_calls2 = variant_calls
            .iter()
            .filter(|x| {
                x.position() + 1 >= realign_segments[1] && x.position() + 1 < realign_segments[2]
            })
            .map(|x| x.clone())
            .collect::<Vec<_>>();
        let mut read_seq2 = global_realignment(
            &read,
            &variant_calls2,
            &vec![],
            &ref_seq,
            realign_segments[1],
            realign_segments[2],
            read_parameters.min_base_quality,
            &region_coordinates.clip_variant_sites,
        )?;
        read_seq.append(&mut read_seq2);

        let variant_calls3 = variant_calls
            .iter()
            .filter(|x| x.position() + 1 >= realign_segments[3])
            .map(|x| x.clone())
            .collect::<Vec<_>>();
        let mut read_seq3 = global_realignment(
            &read,
            &variant_calls3,
            &vec![],
            &ref_seq,
            realign_segments[3],
            region_coordinates.repeat_len as i64,
            read_parameters.min_base_quality,
            &region_coordinates.clip_variant_sites,
        )?;
        read_seq.append(&mut read_seq3);

        debug!(
            "segment_name {segment_name} read_simplified {:?}",
            std::str::from_utf8(&read_seq)?
        );
        read_segment_raw_fp.insert(segment_name, read_seq);
    }

    Ok(read_segment_raw_fp)
}

/// Global realignment of a read
/// # Arguments
/// * `read` - read
/// * `variant_calls` - variant calls
/// * `hom_calls` - homozygous variant calls
/// * `chrom_seq` - reference sequence
/// * `region_start` - start of the region
/// * `region_end` - end of the region
/// * `min_qual` - minimum quality score
/// * `clip_variant_sites` - clip variant sites
/// # Returns
/// * `Result<Vec<u8>, DError>` - read represented at variant sites
pub fn global_realignment(
    read: &bam::Record,
    variant_calls: &[Variant],
    hom_calls: &[Variant],
    chrom_seq: &[u8],
    region_start: i64,
    region_end: i64,
    min_qual: u8,
    clip_variant_sites: &BTreeMap<i64, Vec<u8>>,
) -> Result<Vec<u8>, DError> {
    let mut read_simplified = Vec::new();
    let mut pos_to_variant: BTreeMap<i64, usize> = BTreeMap::new();
    for variant in variant_calls {
        *pos_to_variant.entry(variant.position()).or_default() += 1;
    }
    let total_num_pos = pos_to_variant.len();

    let num_variants: usize = variant_calls.len();

    let qname = std::str::from_utf8(read.qname())?;
    let read_start_pos = start_pos_on_read(&read);
    let reference_start_pos = &read.pos();
    let reference_end_pos = &read.reference_end();
    let aln_len = reference_end_pos - reference_start_pos;
    let segment_name = format!("{qname}:{}:{}", read_start_pos, aln_len);
    //build a lookup from reference coordinate -> sequence coordinate
    let mut coordinate_lookup: HashMap<i64, i64> = Default::default();
    let mut min_position: i64 = i64::MAX;
    let mut max_position: i64 = i64::MIN;
    for bp in read.aligned_pairs() {
        let segment_index = bp[0];
        let ref_index = bp[1];
        if ref_index < region_end && ref_index >= region_start {
            coordinate_lookup.insert(ref_index, segment_index);
            min_position = min_position.min(ref_index);
            max_position = max_position.max(ref_index);
        }
    }

    if max_position < min_position {
        return Ok(vec![b'x'; total_num_pos]);
    }

    // consider soft clip site
    let adj_pos = *reference_end_pos + 1;
    // max_position is the last one that we found, so add +1 to include in the range
    let aligned_range = if !clip_variant_sites.contains_key(&adj_pos) {
        min_position..(max_position + 1)
    } else {
        min_position..(max_position + 4)
    };

    //we will populate these with the variant level info
    let mut num_overlaps: usize = 0;
    let mut first_overlap: Option<usize> = None;
    let mut last_overlap: usize = 0;
    for (i, variant) in variant_calls.iter().enumerate() {
        let variant_pos: i64 = variant.position();
        if aligned_range.contains(&variant_pos) {
            if first_overlap.is_none() {
                first_overlap = Some(i);
            }
            last_overlap = i + 1;
            num_overlaps += 1;
        }
    }

    // if this mapping overlaps no alleles, then there's no reason to look at it anymore
    if num_overlaps == 0 {
        return Ok(vec![b'x'; total_num_pos]);
    }

    // convert into a non-option
    let Some(first_overlap) = first_overlap else {
        return Ok(vec![b'x'; total_num_pos]);
    };
    assert_eq!(num_overlaps, last_overlap - first_overlap);

    // check for homozygous variants also
    let mut first_hom_overlap: Option<usize> = None;
    let mut last_hom_overlap: usize = 0;
    for (i, variant) in hom_calls.iter().enumerate() {
        let variant_pos: i64 = variant.position();
        if aligned_range.contains(&variant_pos) {
            if first_hom_overlap.is_none() {
                first_hom_overlap = Some(i);
            }
            last_hom_overlap = i + 1;
        }
    }
    let first_hom_overlap: usize = first_hom_overlap.unwrap_or(0);

    // .seq() returns Seq<'_> type, but we should just full decode
    let read_sequence: Vec<u8> = read.seq().as_bytes();
    let read_qualities: &[u8] = read.qual();
    assert_eq!(read_sequence.len(), read_qualities.len());

    // these should always exist based on how we set it up
    let Some(read_start) = coordinate_lookup.get(&min_position).copied() else {
        debug!("missing read start coordinate for segment {segment_name} at {min_position}");
        return Ok(vec![b'x'; total_num_pos]);
    };
    let Some(read_end) = coordinate_lookup.get(&max_position).copied() else {
        debug!("missing read end coordinate for segment {segment_name} at {max_position}");
        return Ok(vec![b'x'; total_num_pos]);
    };
    let read_start: usize = read_start as usize;
    let mut read_end: usize = read_end as usize;

    // consider soft clip site
    if clip_variant_sites.contains_key(&adj_pos) {
        let distance_to_end = read_sequence.len() - 1 - read_end;
        read_end += std::cmp::min(distance_to_end, 3);
    }

    // pull out the part of the read we're aligning against
    let read_align: &[u8] = &read_sequence[read_start..(read_end + 1)];

    /*
    Current state:
    - we have the reference genome
    - we have the part of the read that aligns in `read_align`, the full read sequence in `read_sequence`
    - we have the indices of the first and last variant overlaps in `first_overlap` and `last_overlap`

    We need to populate:
    - alleles
    - quals
    - read stats (see below)

    Game plan:
    - construct a graph representing just this reference location + relevant alleles
    - while constructing, assign alleles to each new branch (it may be reference allele)
    -- IF you have multiple alleles starting at the same coordinate (e.g. identical call), then do not create an in-between node; this should resolve in the tie-breaking as "identical"
    -- so each branch should get a variant index + an allele assignment (0/1); reference alleles may end up with multiple 0 alleles in the event of multi-start
    - align the read via POA
    - look at the traversed nodes and copy the allele assignments; if anything is unassigned at the end, it gets 2; any with conflicting assignments get 2 also
    - update stats according to the assignments, we can't really do exact right now (maybe we can look at score deltas from one node to the next?)
    */

    // we need to also provide any preset alleles
    let ref_end = if !clip_variant_sites.contains_key(&adj_pos) {
        max_position as usize + 1
    } else {
        max_position as usize + 4
    };
    let (wfa_graph, node_to_alleles): (WFAGraph, NodeAlleleMap) =
        WFAGraph::from_reference_variants_with_hom(
            chrom_seq,
            &variant_calls[first_overlap..last_overlap], // these are both range style indices
            &hom_calls[first_hom_overlap..last_hom_overlap],
            min_position as usize,
            ref_end,
            500,
        )
        .map_err(|e| {
            std::io::Error::other(format!(
                "Failed to build WFA graph for segment {segment_name} over reference range {min_position}-{max_position}: {e}"
            ))
        })?;

    // pass through for the WFA errors now
    let wfa_result: WFAResult = match wfa_graph.edit_distance_with_pruning(read_align, 500) {
        Ok(r) => r,
        Err(_e) => {
            debug!("realignment failed on read {segment_name}...");
            return Ok(vec![b'x'; total_num_pos]);
        }
    };

    trace!(
        "read {} WFAGraph result () => num_nodes: {}, read_len: {}, variant_overlaps: {}, edit_distance: {}", 
        segment_name, wfa_graph.get_num_nodes(), max_position-min_position+1, num_overlaps, wfa_result.score()
    );

    //we will populate these with the variant level info
    let mut alleles: Vec<AlleleType> = vec![AlleleType::NoOverlap; num_variants];
    for traversed_index in wfa_result.traversed_nodes().iter() {
        for &(var_index, allele_assignment) in node_to_alleles
            .get(traversed_index)
            .unwrap_or(&vec![])
            .iter()
        {
            let correct_index: usize = first_overlap + var_index;
            let this_allele_type =
                AlleleType::from_repr(allele_assignment).unwrap_or(AlleleType::NoOverlap);
            if alleles[correct_index] == AlleleType::NoOverlap {
                alleles[correct_index] = this_allele_type;
            } else if alleles[correct_index] != this_allele_type {
                alleles[correct_index] = AlleleType::Ambiguous;
            }
        }
    }

    let mut deleted_positions = Vec::new();
    let mut pos_index = 0;
    for (pos, num_var) in pos_to_variant {
        let mut this_pos_variants = Vec::new();
        let mut allele_calls = Vec::new();
        for i in pos_index..(pos_index + num_var) {
            let allele_call = alleles[i];
            let this_variant = &variant_calls[i];
            if this_variant.get_type() == VariantType::Deletion
                && allele_call == AlleleType::Alternate
            {
                let ref_len = this_variant.get_ref_len();
                let pos1 = pos as usize;
                for deleted_pos in (pos1 + 1)..(pos1 + ref_len) {
                    deleted_positions.push(deleted_pos as i64);
                }
            }
            allele_calls.push(allele_call);
            this_pos_variants.push((
                std::str::from_utf8(this_variant.get_allele0())?,
                std::str::from_utf8(this_variant.get_allele1())?,
            ));
        }
        debug!(
            "seg_name {segment_name} pos {} {:?} {:?}",
            pos, this_pos_variants, allele_calls
        );

        let mut low_qual = false;
        if coordinate_lookup.contains_key(&pos) {
            let Some(pos_on_read) = coordinate_lookup.get(&pos).copied() else {
                pos_index += num_var;
                continue;
            };
            let pos_on_read = pos_on_read as usize;
            let this_pos_qual = read_qualities[pos_on_read];
            if this_pos_qual < min_qual {
                low_qual = true;
                read_simplified.push(b'-');
                debug!("pos {pos} low_qual {this_pos_qual}",);
            }
        }
        if !low_qual {
            if allele_calls.len() == 1 {
                let allele_call = allele_calls[0];
                if allele_call == AlleleType::NoOverlap {
                    if deleted_positions.contains(&pos) {
                        read_simplified.push(b'0');
                    } else if pos > *reference_start_pos && pos < *reference_end_pos {
                        read_simplified.push(b'-');
                    } else {
                        read_simplified.push(b'x');
                    }
                } else {
                    if allele_call == AlleleType::Reference {
                        read_simplified.push(b'0');
                    } else if allele_call == AlleleType::Alternate {
                        read_simplified.push(b'1');
                    } else if allele_call == AlleleType::Ambiguous {
                        let check_alt = close_check_variant(
                            pos,
                            pos_index,
                            variant_calls,
                            hom_calls,
                            chrom_seq,
                            &read_sequence,
                            &coordinate_lookup,
                        );
                        if let Some(found_alt) = check_alt {
                            if found_alt {
                                read_simplified.push(b'1');
                                debug!("update pos {pos} to alternate");
                            } else {
                                read_simplified.push(b'0');
                                debug!("update pos {pos} to reference");
                            }
                        } else {
                            read_simplified.push(b'-');
                        }
                    } else {
                        read_simplified.push(b'-');
                    }
                }
            } else {
                if allele_calls == vec![AlleleType::NoOverlap; num_var] {
                    if deleted_positions.contains(&pos) {
                        read_simplified.push(b'0');
                    } else if pos > *reference_start_pos && pos < *reference_end_pos {
                        read_simplified.push(b'-');
                    } else {
                        read_simplified.push(b'x');
                    }
                } else if allele_calls == vec![AlleleType::Reference; num_var] {
                    read_simplified.push(b'0');
                } else {
                    let alt_count = allele_calls
                        .iter()
                        .filter(|x| **x == AlleleType::Alternate)
                        .count();
                    if alt_count == 1 {
                        if let Some(alt_index) = allele_calls
                            .iter()
                            .position(|x| *x == AlleleType::Alternate)
                        {
                            if alt_index == 0 {
                                read_simplified.push(b'1');
                            } else if alt_index <= 8 {
                                let allele_name = alt_index + 1;
                                read_simplified
                                    .extend_from_slice(allele_name.to_string().as_bytes());
                            } else {
                                read_simplified.push(b'-');
                            }
                        } else {
                            read_simplified.push(b'-');
                        }
                    } else {
                        read_simplified.push(b'-');
                    }
                }
            }
        }
        pos_index += num_var;
    }
    Ok(read_simplified)
}

/// Check read and reference sequence to update the call if the global realignment returns ambiguous
/// # Arguments
/// * `pos` - position of the variant
/// * `pos_index` - index of the variant
/// * `variant_calls` - variant calls
/// * `hom_calls` - homozygous variant calls
/// * `chrom_seq` - reference sequence
/// * `read_sequence` - read sequence
/// * `coordinate_lookup` - lookup from reference coordinate -> sequence coordinate
/// # Returns
/// * `Option<bool>` - true for found alternate, false for reference, None for ambiguous
pub fn close_check_variant(
    pos: i64,
    pos_index: usize,
    variant_calls: &[Variant],
    hom_calls: &[Variant],
    chrom_seq: &[u8],
    read_sequence: &Vec<u8>,
    coordinate_lookup: &HashMap<i64, i64>,
) -> Option<bool> {
    let pos_usize = pos as usize;
    let this_variant_ref = variant_calls[pos_index].get_allele0();
    let this_variant_allele1 = variant_calls[pos_index].get_allele1();
    if pos >= 2 && this_variant_allele1.len() == 1 && this_variant_ref.len() == 1 {
        // first define a region flanking the variant position
        let left_pos_on_ref = pos - 2;
        let right_pos_on_ref = pos + 2;
        let mut left_pos_on_read: Option<i64> = None;
        let mut right_pos_on_read: Option<i64> = None;
        if coordinate_lookup.contains_key(&left_pos_on_ref) {
            left_pos_on_read = coordinate_lookup.get(&left_pos_on_ref).copied();
        } else if coordinate_lookup.contains_key(&(left_pos_on_ref - 1)) {
            left_pos_on_read = coordinate_lookup
                .get(&(left_pos_on_ref - 1))
                .copied()
                .map(|value| value + 1);
        }
        if coordinate_lookup.contains_key(&right_pos_on_ref) {
            right_pos_on_read = coordinate_lookup.get(&right_pos_on_ref).copied();
        } else if coordinate_lookup.contains_key(&(right_pos_on_ref + 1)) {
            right_pos_on_read = coordinate_lookup
                .get(&(right_pos_on_ref + 1))
                .copied()
                .map(|value| value - 1);
        }

        if let (Some(left_pos_on_read), Some(right_pos_on_read)) =
            (left_pos_on_read, right_pos_on_read)
        {
            let region_length_on_ref = right_pos_on_ref + 1 - left_pos_on_ref;
            let this_range = (left_pos_on_ref as usize)..(right_pos_on_ref as usize + 1);
            let ref_seq_this_range = &chrom_seq[this_range.clone()];
            let ref_seq_left = &chrom_seq[(left_pos_on_ref as usize)..pos_usize];
            let ref_seq_right = &chrom_seq[pos_usize + 1..(right_pos_on_ref as usize + 1)];
            trace!(
                "reference sequence in this region is {:?}",
                String::from_utf8_lossy(ref_seq_this_range)
            );

            let variants_in_range = variant_calls
                .iter()
                .filter(|x| x.position() != pos && this_range.contains(&(x.position() as usize)))
                .map(|x| x.get_allele1())
                .collect::<Vec<_>>();
            let hom_variants_in_range = hom_calls
                .iter()
                .filter(|x| x.position() != pos && this_range.contains(&(x.position() as usize)))
                .map(|x| x.get_allele1())
                .collect::<Vec<_>>();
            let left_pos_on_read = left_pos_on_read as usize;
            let right_pos_on_read = right_pos_on_read as usize;
            let read_seq = &read_sequence[left_pos_on_read..(right_pos_on_read + 1)];
            let read_seq_len = read_seq.len() as i64;
            trace!(
                "read sequence in this region from pos {} to pos {} is {:?}",
                left_pos_on_read,
                right_pos_on_read,
                String::from_utf8_lossy(read_seq)
            );

            // read sequence does not contain expected variant and there is no deletion in this region on the read
            if !read_seq.contains(&this_variant_allele1[0]) && read_seq_len >= region_length_on_ref
            {
                return Some(false);
            }
            if !ref_seq_this_range.contains(&this_variant_allele1[0])
                && !variants_in_range.contains(&this_variant_allele1)
                && !hom_variants_in_range.contains(&this_variant_allele1)
                && read_seq.contains(&this_variant_allele1[0])
            {
                return Some(true);
            }

            // read sequence does not contain expected reference and there is no deletion in this region on the read
            // disable for now
            //if !read_seq.contains(&this_variant_ref[0]) && read_seq_len >= region_length_on_ref {
            //    return Some(true);
            //}
            if !ref_seq_left.contains(&this_variant_ref[0])
                && !ref_seq_right.contains(&this_variant_ref[0])
                && !variants_in_range.contains(&this_variant_ref)
                && !hom_variants_in_range.contains(&this_variant_ref)
                && read_seq.contains(&this_variant_ref[0])
            {
                return Some(false);
            }
        }

        // special case for homopolymer, where we allow deletions in the read
        let mut homopolymer_range_to_check = 0..0;
        // variant position is in the middle of a homopolymer
        if &chrom_seq[pos_usize - 1] == &this_variant_ref[0]
            && &chrom_seq[pos_usize + 1] == &this_variant_ref[0]
        {
            if &chrom_seq[pos_usize - 2] == &this_variant_ref[0] {
                homopolymer_range_to_check = (pos_usize - 2)..(pos_usize + 2);
            } else if &chrom_seq[pos_usize + 2] == &this_variant_ref[0] {
                homopolymer_range_to_check = (pos_usize - 1)..(pos_usize + 3);
            }
        } else if &chrom_seq[pos_usize - 1] != &this_variant_ref[0]
            && &chrom_seq[pos_usize + 1] == &this_variant_ref[0]
            && &chrom_seq[pos_usize + 2] == &this_variant_ref[0]
            && &chrom_seq[pos_usize + 3] == &this_variant_ref[0]
        // variant position is at the start of a homopolymer
        {
            homopolymer_range_to_check = (pos_usize - 1)..(pos_usize + 4);
        } else if pos_usize >= 3
            && &chrom_seq[pos_usize + 1] != &this_variant_ref[0]
            && &chrom_seq[pos_usize - 1] == &this_variant_ref[0]
            && &chrom_seq[pos_usize - 2] == &this_variant_ref[0]
            && &chrom_seq[pos_usize - 3] == &this_variant_ref[0]
        // variant position is at the end of a homopolymer
        {
            homopolymer_range_to_check = (pos_usize - 3)..(pos_usize + 2);
        }

        if !homopolymer_range_to_check.is_empty() {
            let ref_seq_this_range = &chrom_seq[homopolymer_range_to_check.clone()];
            trace!(
                "homopolymer_range_to_check {:?}",
                homopolymer_range_to_check
            );
            trace!(
                "reference sequence in this region is {:?}",
                String::from_utf8_lossy(ref_seq_this_range)
            );
            let bases_on_read = homopolymer_range_to_check
                .collect::<Vec<usize>>()
                .iter()
                .filter(|x| coordinate_lookup.contains_key(&(**x as i64)))
                .map(|x| coordinate_lookup[&(*x as i64)])
                .collect::<Vec<i64>>();
            if !bases_on_read.is_empty() {
                let Some(read_start) = bases_on_read.iter().min() else {
                    return None;
                };
                let Some(read_end) = bases_on_read.iter().max() else {
                    return None;
                };
                let bases_on_read = &read_sequence[*read_start as usize..(*read_end as usize + 1)];
                trace!("bases_on_read {:?}", String::from_utf8_lossy(bases_on_read));
                if bases_on_read.contains(&this_variant_ref[0])
                    && !bases_on_read.contains(&this_variant_allele1[0])
                {
                    return Some(false);
                }
            }
        }
    }
    return None;
}

#[cfg(test)]
mod tests {
    use super::close_check_variant;
    use crate::realignment::utilities::Variant;
    use std::collections::HashMap;

    fn create_test_variant(position: i64, ref_allele: u8, alt_allele: u8) -> Variant {
        Variant::new_snv(0, position, vec![ref_allele], vec![alt_allele], 0, 1).unwrap()
    }

    fn create_coordinate_lookup(positions: &[(i64, i64)]) -> HashMap<i64, i64> {
        positions.iter().cloned().collect()
    }

    #[test]
    fn test_close_check_variant_alternate_found() {
        // Test case: alternate allele found in read sequence
        let pos = 10;
        let variant = create_test_variant(pos, b'A', b'C');
        let variant_calls = vec![variant.clone()];
        let hom_calls = vec![];

        // Reference sequence: ATGCGATGTG...ATGCG...
        let chrom_seq = b"ATGCGATGTGATGCG";

        // Read sequence contains the alternate 'C' at position 10
        let read_sequence = b"ATGCGATGTGCTGCG".to_vec();

        // Coordinate lookup: ref pos -> read pos (1:1 mapping for simplicity)
        let mut coordinate_lookup = HashMap::new();
        for i in 0..15 {
            coordinate_lookup.insert(i, i);
        }

        let result = close_check_variant(
            pos,
            0,
            &variant_calls,
            &hom_calls,
            chrom_seq,
            &read_sequence,
            &coordinate_lookup,
        );

        assert_eq!(result, Some(true));
    }

    #[test]
    fn test_close_check_variant_alternate_in_reference_range() {
        // Test case: alternate allele exists in reference range (should not return true)
        let pos = 10;
        let variant = create_test_variant(pos, b'A', b'C');
        let variant_calls = vec![variant.clone()];
        let hom_calls = vec![];

        // Reference has 'C' in the range (at pos 11)
        let chrom_seq = b"ATGCGATGTGACGCG";

        // Read also has 'C'
        let read_sequence = b"ATGCGATGTGCCGCG".to_vec();

        let mut coordinate_lookup = HashMap::new();
        for i in 0..15 {
            coordinate_lookup.insert(i, i);
        }

        let result = close_check_variant(
            pos,
            0,
            &variant_calls,
            &hom_calls,
            chrom_seq,
            &read_sequence,
            &coordinate_lookup,
        );

        assert!(result.is_none());
    }

    #[test]
    fn test_close_check_variant_with_hom_calls() {
        // Test case: with homozygous variant calls in range
        let pos = 10;
        let variant = create_test_variant(pos, b'A', b'C');
        let hom_variant = create_test_variant(9, b'G', b'C');
        let variant_calls = vec![variant.clone()];
        let hom_calls = vec![hom_variant];

        // Reference sequence: ...ATGCG...
        let chrom_seq = b"ATGCGATGTGATGCG";

        // Read sequence contains the alternate 'C' at position 10
        let read_sequence = b"ATGCGATGTGCTGCG".to_vec();

        let mut coordinate_lookup = HashMap::new();
        for i in 0..15 {
            coordinate_lookup.insert(i, i);
        }

        let result = close_check_variant(
            pos,
            0,
            &variant_calls,
            &hom_calls,
            chrom_seq,
            &read_sequence,
            &coordinate_lookup,
        );

        assert!(result.is_none());
    }

    #[test]
    fn test_close_check_variant_with_variants_in_range() {
        // Test case: other variants in the range should be excluded
        let pos = 10;
        let variant = create_test_variant(pos, b'A', b'C');
        let other_variant = create_test_variant(9, b'G', b'C');
        let variant_calls = vec![other_variant.clone(), variant.clone()];
        let hom_calls = vec![];

        let chrom_seq = b"ATGCGATGTGATGCG";
        // Read has 'G' at pos 10, but also has 'C' at pos 9 (from other variant)
        let read_sequence = b"ATGCGATGTGCTGCG".to_vec();

        let mut coordinate_lookup = HashMap::new();
        for i in 0..15 {
            coordinate_lookup.insert(i, i);
        }

        let result = close_check_variant(
            pos,
            1, // pos_index is 1 because variant is second in the list
            &variant_calls,
            &hom_calls,
            chrom_seq,
            &read_sequence,
            &coordinate_lookup,
        );

        assert!(result.is_none());
    }

    #[test]
    fn test_close_check_variant_with_deletion() {
        // Test case: coordinate lookup uses adjacent positions when exact not found
        let pos = 10;
        let variant = create_test_variant(pos, b'A', b'C');
        let variant_calls = vec![variant.clone()];
        let hom_calls = vec![];

        // note base before deletion is C, the expected variant
        let chrom_seq = b"ATGCGATCTGATGCG";
        // deletion of T at pos 8
        // ATGCGATC...GCTGCG
        let read_sequence = b"ATGCGATGGCTGCG".to_vec();

        // Missing exact positions, but have adjacent ones
        let coordinate_lookup = create_coordinate_lookup(&[
            (7, 8),
            (9, 9),
            (10, 10), // pos
            (12, 12),
            (13, 13),
        ]);

        let result = close_check_variant(
            pos,
            0,
            &variant_calls,
            &hom_calls,
            chrom_seq,
            &read_sequence,
            &coordinate_lookup,
        );

        assert_eq!(result, Some(true));
    }

    #[test]
    fn test_close_check_variant_alternate_not_found() {
        // Test case: alternate allele found in read sequence
        let pos = 10;
        let variant = create_test_variant(pos, b'A', b'C');
        let variant_calls = vec![variant.clone()];
        let hom_calls = vec![];

        // Reference sequence: ...ATGCG...
        let chrom_seq = b"ATGCGATGTGATGCG";

        // AAAAA in the range in the read
        let read_sequence = b"ATGCGATGAAAAACG".to_vec();

        // Coordinate lookup: ref pos -> read pos (1:1 mapping for simplicity)
        let mut coordinate_lookup = HashMap::new();
        for i in 0..15 {
            coordinate_lookup.insert(i, i);
        }

        let result = close_check_variant(
            pos,
            0,
            &variant_calls,
            &hom_calls,
            chrom_seq,
            &read_sequence,
            &coordinate_lookup,
        );

        assert_eq!(result, Some(false));
    }

    #[test]
    fn test_close_check_variant_reference_found() {
        // Test case: reference allele found in read sequence
        let pos = 10;
        let variant = create_test_variant(pos, b'A', b'C');
        let variant_calls = vec![variant.clone()];
        let hom_calls = vec![];

        // Reference sequence: ...ATGCG... (has 'A' at pos 10)
        let chrom_seq = b"ATGCGATGTGATGCG";

        // Read sequence contains the reference 'A' at position 10
        let read_sequence = b"ATGCGATGTGATGCG".to_vec();

        // Coordinate lookup
        let mut coordinate_lookup = HashMap::new();
        for i in 0..15 {
            coordinate_lookup.insert(i, i);
        }

        let result = close_check_variant(
            pos,
            0,
            &variant_calls,
            &hom_calls,
            chrom_seq,
            &read_sequence,
            &coordinate_lookup,
        );

        assert_eq!(result, Some(false));
    }

    #[test]
    fn test_close_check_variant_homopolymer_middle() {
        // Test case: variant in middle of homopolymer (AAA -> AAG)
        let pos = 10;
        let variant = create_test_variant(pos, b'A', b'T');
        let variant_calls = vec![variant.clone()];
        let hom_calls = vec![];

        // Reference: ATGCATGCGAA...AAA... (homopolymer)
        // 0        9
        // ATGCATGCGAAAAATGCG
        // ||||||||||||||||||
        // ATGCATGCG-AAAATGCG
        let chrom_seq = b"ATGCATGCGAAAAATGCG";

        // Read sequence
        // delete one A in the homopolymer
        let read_sequence = b"ATGCATGCGAAAATGCG".to_vec();

        let mut coordinate_lookup = HashMap::new();
        for i in 0..9 {
            coordinate_lookup.insert(i, i);
        }
        for i in 10..18 {
            coordinate_lookup.insert(i, i - 1);
        }

        let result = close_check_variant(
            pos,
            0,
            &variant_calls,
            &hom_calls,
            chrom_seq,
            &read_sequence,
            &coordinate_lookup,
        );

        assert_eq!(result, Some(false));
    }

    #[test]
    fn test_close_check_variant_homopolymer_start() {
        // Test case: variant at start of homopolymer
        let pos = 9;
        let variant = create_test_variant(pos, b'A', b'T');
        let variant_calls = vec![variant.clone()];
        let hom_calls = vec![];

        // Reference: ...TAAAA... (homopolymer starts after T)
        let chrom_seq = b"ATGCATGCGAAAAATGCG";

        let read_sequence = b"ATGCATGCGAAAATGCG".to_vec();

        let mut coordinate_lookup = HashMap::new();
        for i in 0..9 {
            coordinate_lookup.insert(i, i);
        }
        for i in 10..18 {
            coordinate_lookup.insert(i, i - 1);
        }

        let result = close_check_variant(
            pos,
            0,
            &variant_calls,
            &hom_calls,
            chrom_seq,
            &read_sequence,
            &coordinate_lookup,
        );

        assert_eq!(result, Some(false));
    }

    #[test]
    fn test_close_check_variant_homopolymer_start_starting_base_same_as_variant() {
        // Test case: variant at start of homopolymer
        let pos = 9;
        let variant = create_test_variant(pos, b'A', b'G');
        let variant_calls = vec![variant.clone()];
        let hom_calls = vec![];

        // Reference: ...TAAAA... (homopolymer starts after T)
        let chrom_seq = b"ATGCATGCGAAAAATGCG";

        let read_sequence = b"ATGCATGCGAAAATGCG".to_vec();

        let mut coordinate_lookup = HashMap::new();
        for i in 0..9 {
            coordinate_lookup.insert(i, i);
        }
        for i in 10..18 {
            coordinate_lookup.insert(i, i - 1);
        }

        let result = close_check_variant(
            pos,
            0,
            &variant_calls,
            &hom_calls,
            chrom_seq,
            &read_sequence,
            &coordinate_lookup,
        );

        assert_eq!(result, None);
    }

    #[test]
    fn test_close_check_variant_homopolymer_end() {
        // Test case: variant at end of homopolymer
        let pos = 13;
        let variant = create_test_variant(pos, b'A', b'G');
        let variant_calls = vec![variant.clone()];
        let hom_calls = vec![];

        // Reference: ...AAAAT... (homopolymer ends before T)
        let chrom_seq = b"ATGCATGCGAAAAATGCG";

        let read_sequence = b"ATGCATGCGAAAATGCG".to_vec();

        let mut coordinate_lookup = HashMap::new();
        for i in 0..9 {
            coordinate_lookup.insert(i, i);
        }
        for i in 10..18 {
            coordinate_lookup.insert(i, i - 1);
        }

        let result = close_check_variant(
            pos,
            0,
            &variant_calls,
            &hom_calls,
            chrom_seq,
            &read_sequence,
            &coordinate_lookup,
        );

        assert_eq!(result, Some(false));
    }

    #[test]
    fn test_close_check_variant_homopolymer_end_ending_base_same_as_variant() {
        // Test case: variant at end of homopolymer
        let pos = 13;
        let variant = create_test_variant(pos, b'A', b'T');
        let variant_calls = vec![variant.clone()];
        let hom_calls = vec![];

        // Reference: ...AAAAT... (homopolymer ends before T)
        let chrom_seq = b"ATGCATGCGAAAAATGCG";

        let read_sequence = b"ATGCATGCGAAAATGCG".to_vec();

        let mut coordinate_lookup = HashMap::new();
        for i in 0..9 {
            coordinate_lookup.insert(i, i);
        }
        for i in 10..18 {
            coordinate_lookup.insert(i, i - 1);
        }

        let result = close_check_variant(
            pos,
            0,
            &variant_calls,
            &hom_calls,
            chrom_seq,
            &read_sequence,
            &coordinate_lookup,
        );

        assert_eq!(result, None);
    }
}
