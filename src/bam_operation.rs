use crate::methylation::{get_methyl_prob, get_methyl_tags};
use crate::util::{
    append_kivvi_pg_header, resolve_chrom_name_from_header, DError, DResult, FlankReads,
    RegionCoordinates,
};
use log::{debug, trace, warn};
use minimap2::Built;
use minimap2::{ffi, Aligner};
use rust_htslib::bam::header::HeaderRecord;
use rust_htslib::bam::record::Aux;
use rust_htslib::bam::record::CigarString;
use rust_htslib::bam::{
    self, ext::BamRecordExtensions, record::Cigar, Format, Header, HeaderView, Read, Writer,
};
use rust_htslib::faidx;
use std::collections::BTreeMap;
use std::collections::HashSet;
use std::path::PathBuf;

/// Realign reads to repeat unit and filter alignments
/// # Arguments
/// * `bam_name` - input bam, wgs
/// * `region_coordinates` - coordinates defined in the region
/// * `reference` - reference file
/// * `realigned_bam` - output bam
/// * `get_methyl` - whether to get methyl tags
/// # Returns
/// * `Vec<bam::Record>` - realigned records
/// * `BTreeMap<String, usize>` - read lengths lookup
/// * `Writer` - writer for realigned bam
/// * `BTreeMap<String, (String, Vec<u8>)>` - methyl tags
/// * `BTreeMap<String, Vec<u8>>` - methyl probabilities lookup
pub fn realign(
    bam_name: PathBuf,
    region_coordinates: RegionCoordinates,
    reference: &PathBuf,
    realigned_bam: PathBuf,
    get_methyl: bool,
) -> Result<
    (
        Vec<bam::Record>,
        BTreeMap<String, usize>,
        Writer,
        BTreeMap<String, (String, Vec<u8>)>,
        BTreeMap<String, Vec<u8>>,
    ),
    DError,
> {
    let mut methyl_tags = BTreeMap::new();
    let mut methyl_probs = BTreeMap::new();
    let mut read_length = BTreeMap::new();
    let mut bam_reader = bam::IndexedReader::from_path(bam_name)?;
    let ref_reader = faidx::Reader::from_path(reference)?;
    let ref_names = ref_reader.seq_names()?;
    let mut realn_records = Vec::new();
    let mut aligner = Aligner::builder()
        .with_cigar()
        .with_index(reference, None)?;
    aligner.mapopt.flag |= ffi::MM_F_SOFTCLIP as i64;
    aligner.idxopt.k = 19;
    aligner.idxopt.flag = 0; // Disable homopolymer compression

    let mut header = Header::new();
    aligner.populate_header(&mut header);
    append_kivvi_pg_header(&mut header);
    let header_view = HeaderView::from_header(&header);
    let writer = Writer::from_path(&realigned_bam, &header, Format::Bam)?;
    let regions = region_coordinates.extract_regions;
    for region in regions {
        let fields = region.split_terminator(':').collect::<Vec<_>>();
        let nchr = fields[0];
        let coords = fields[1].split_terminator('-').collect::<Vec<_>>();
        let start = coords[0].parse::<i64>()?;
        let stop = coords[1].parse::<i64>()?;
        let resolved_chr = resolve_chrom_name_from_header(bam_reader.header(), nchr)
            .ok_or_else(|| format!("Chromosome '{nchr}' not found in input BAM header"))?;
        debug!("fetching region: {resolved_chr:?} {start:?} {stop:?}");
        bam_reader
            .fetch((resolved_chr.as_str(), start, stop))
            .map_err(|e| format!("Failed to fetch region {resolved_chr}:{start}-{stop}: {e}"))?;

        for read in bam_reader.records() {
            let record = read?;
            if !record.is_secondary() && !record.is_supplementary() {
                if get_methyl {
                    if let Some((mm, ml)) = get_methyl_tags(&record)? {
                        let qname = std::str::from_utf8(record.qname())?.to_string();
                        methyl_tags.insert(qname, (mm, ml));
                    }
                    if let Some(meth_prob) = get_methyl_prob(&record)? {
                        let qname = std::str::from_utf8(record.qname())?.to_string();
                        methyl_probs.insert(qname, meth_prob);
                    }
                }
                let qname = std::str::from_utf8(record.qname())?;
                read_length
                    .entry(qname.to_string())
                    .or_insert(record.seq().len());
                let alignments = seq2seq(record, &header_view, &header, &ref_names, &aligner)?;
                for realn_record in alignments {
                    realn_records.push(realn_record);
                }
            }
        }
    }
    realn_records.sort_by(|a, b| a.tid().cmp(&b.tid()).then(a.pos().cmp(&b.pos())));

    Ok((
        realn_records,
        read_length,
        writer,
        methyl_tags,
        methyl_probs,
    ))
}

/// Align a `bam::Record` with an aligner and convert to a `bam::Record`.
/// Aligns sequence-orientation reads, adds sam tags + qualities.
/// Also corrects cigar representation to be SAM-compatible.
/// # Arguments
/// * `input_record` - input bam record
/// * `header_view` - header view
/// * `header` - header
/// * `ref_names` - reference sequence names
/// * `aligner` - minimap2 aligner
/// # Returns
/// * `Vec<bam::Record>` - realigned records
pub fn seq2seq(
    input_record: bam::Record,
    header_view: &bam::HeaderView,
    header: &bam::Header,
    ref_names: &Vec<String>,
    aligner: &minimap2::Aligner<Built>,
) -> Result<Vec<bam::Record>, DError> {
    // qual needs + 33 for use in map_to_sam, but it needs to be back to unshifted
    // for use in bam::Record.
    let mut qual = input_record
        .qual()
        .iter()
        .map(|x| *x + 33)
        .collect::<Vec<_>>();
    let original_qual = input_record.qual();
    let is_reverse = input_record.is_reverse();
    let mut seq = input_record.seq().as_bytes();
    seq.make_ascii_uppercase();
    let original_seq = seq.clone();
    if is_reverse {
        let seq1 = reverse_complement(&seq, &mut qual);
        seq = seq1;
    }
    let qname = input_record.qname();

    let mut mappings = aligner.map(
        &seq,
        /* output_cigar= */ true,
        /* output_md= */ true,
        /* max_frag_len= */ None,
        /* extra_flags= */ None,
        Some(&qname),
    )?;

    let sam_records =
        aligner.map_to_sam(&seq, Some(&qual), Some(&qname), header_view, None, None)?;
    qual.iter_mut().for_each(|x| *x -= 33);
    let records = mappings
        .iter_mut()
        .zip(sam_records)
        .map(|(mapping, sam_mapping)| {
            let mut record = minimap2::htslib::mapping_to_record(
                Some(mapping),
                &seq,
                header.clone(),
                Some(&qual),
                Some(&qname),
            );
            let query_name = std::str::from_utf8(qname).unwrap();
            mapping.query_name = Some(query_name.to_string().into());
            let this_ref_name = mapping.target_name.clone().unwrap();
            let ref_index = ref_names.iter().position(|x| *x == *this_ref_name).unwrap();
            record.set_tid(ref_index as i32);
            for aux in sam_mapping.aux_iter() {
                let (aux_name, aux_field) = aux.expect("Aux error");
                record
                    .push_aux(aux_name, aux_field)
                    .expect("push_aux error");
            }
            (mapping, record)
        })
        .collect::<Vec<_>>();
    let mut alignments = Vec::with_capacity(records.len());
    let original_orientation_tag = if is_reverse { b'R' } else { b'F' };
    for (mapping, mut record) in records {
        assert_eq!(
            mapping.strand == minimap2::Strand::Reverse,
            record.is_reverse()
        );
        record.unset_unmapped();
        record.unset_secondary();
        let mut cigar = mapping
            .alignment
            .as_ref()
            .and_then(|alignment| alignment.cigar.as_ref())
            .unwrap()
            .to_owned();
        let mapping_to_record_cig = record.cigar().to_string();
        let map_to_sam_cig = cigar_to_cigarstr(&cigar).to_string();
        if map_to_sam_cig != mapping_to_record_cig {
            warn!(
                "Cigar for mapping_to_record: {} is different from Cigar found in map_to_sam: {}",
                record.cigar(),
                cigar_to_cigarstr(&cigar)
            );
        }
        // Now add softclips
        let query_len = seq.len() as i32;
        let overhang = query_len - mapping.query_end;
        const SOFT_CLIP: u8 = 4;

        if mapping.query_start > 0 {
            if record.is_reverse() {
                cigar.push((mapping.query_start as u32, SOFT_CLIP)); // soft-clip
            } else {
                cigar.insert(0, (mapping.query_start as u32, SOFT_CLIP)); // soft-clip
            }
        }
        if overhang > 0 {
            if record.is_reverse() {
                cigar.insert(0, (overhang as u32, SOFT_CLIP));
            } else {
                cigar.push((overhang as u32, SOFT_CLIP)); // soft-clip
            }
        }
        let cigar_str = cigar_to_cigarstr(&cigar);
        record.set(&qname, Some(&cigar_str), &original_seq, original_qual);
        record.push_aux(b"or", bam::record::Aux::Char(original_orientation_tag))?;
        debug_assert_eq!(query_length_cigar(&cigar_str) as usize, record.seq_len(), "cigar qlen {} for cigar {cigar_str}/{cigar:?}, mapping {mapping:?} and record {record:?}", query_length_cigar(&cigar_str));
        alignments.push(record);
    }
    Ok(alignments)
}

/// Get the query length from a cigar string
/// # Arguments
/// * `x` - cigar string
/// # Returns
/// * `u32` - query length
#[must_use]
pub fn query_length_cigar(x: &[Cigar]) -> u32 {
    x.iter()
        .copied()
        .filter(|x| consumes_qry(*x))
        .map(rust_htslib::bam::record::Cigar::len)
        .sum::<u32>()
}

/// Check if a cigar operation consumes a query base
/// # Arguments
/// * `x` - cigar operation
/// # Returns
/// * `bool` - true if the cigar operation consumes a query base
#[must_use]
#[inline]
pub fn consumes_qry(x: bam::record::Cigar) -> bool {
    use bam::record::Cigar::{Diff, Equal, Ins, Match, SoftClip};
    matches!(x, Ins(_) | SoftClip(_) | Match(_) | Diff(_) | Equal(_))
}

/// From minimap2-rs
/// Convert minimap2-rs cigar to a cigar string.
/// # Arguments
/// * `cigar` - cigar vector
/// # Returns
/// * `CigarString` - cigar string
fn cigar_to_cigarstr(cigar: &Vec<(u32, u8)>) -> CigarString {
    let op_vec: Vec<Cigar> = cigar
        .to_owned()
        .iter()
        .map(|(len, op)| match op {
            0 => Cigar::Match(*len),
            1 => Cigar::Ins(*len),
            2 => Cigar::Del(*len),
            3 => Cigar::RefSkip(*len),
            4 => Cigar::SoftClip(*len),
            5 => Cigar::HardClip(*len),
            6 => Cigar::Pad(*len),
            7 => Cigar::Equal(*len),
            8 => Cigar::Diff(*len),
            _ => panic!("Unexpected cigar operation"),
        })
        .collect();
    CigarString(op_vec)
}

/// Get the starting position on the read
/// # Arguments
/// * `record` - bam record
/// # Returns
/// * `i64` - starting position on the read
pub fn start_pos_on_read(record: &bam::Record) -> i64 {
    let mut first_clip_len = 0;
    for x in record.cigar().iter() {
        match x {
            Cigar::HardClip(_len) | Cigar::SoftClip(_len) => first_clip_len += i64::from(x.len()),
            _ => break,
        }
    }
    first_clip_len
}

/// Filter alignments based on a more detailed mismatch calculation over a specific region
/// # Arguments
/// * `record` - bam record
/// * `genome_reference` - genome reference
/// * `start` - start position on reference
/// * `end` - end position on reference
/// # Returns
/// * `bool` - true if the record fulfils the interval mismatch criteria
fn interval_mismatch(
    record: &bam::Record,
    genome_reference: Option<&PathBuf>,
    start: i64,
    end: i64,
) -> Result<bool, DError> {
    if genome_reference.is_none() {
        return Ok(true);
    }
    let ref_reader = faidx::Reader::from_path(genome_reference.unwrap())?;
    //let qname = std::str::from_utf8(record.qname())?;
    let mut region_match = 0;
    let mut new_nm = 0;
    let seq = record.seq().as_bytes();
    let tid = record.tid();
    let ref_name = ref_reader.seq_name(tid)?;

    for [read_pos, ref_pos] in record.aligned_pairs() {
        if ref_pos + 1 > start && ref_pos + 1 < end {
            let read_pos = read_pos as usize;
            let ref_pos = ref_pos as usize;
            if let Some(read_base) = seq.get(read_pos) {
                let ref_base = ref_reader.fetch_seq(&ref_name, ref_pos, ref_pos)?;
                region_match += 1;
                if read_base != ref_base.first().unwrap() {
                    new_nm += 1;
                }
            }
        }
    }
    if (new_nm as f64) < (region_match as f64) * 0.005 {
        Ok(true)
    } else {
        Ok(false)
    }
}

/// Get reads overlapping start and end of a VNTR region
/// # Arguments
/// * `bam_name` - input bam, wgs
/// * `region_coordinates` - coordinates defined in the region
/// * `genome_reference` - genome reference
/// # Returns
/// * `FlankReads` - flanking reads
pub fn get_start_end_from_genome(
    bam_name: PathBuf,
    region_coordinates: RegionCoordinates,
    genome_reference: Option<&PathBuf>,
) -> Result<FlankReads, DError> {
    let mut starting_reads_flank = HashSet::new();
    let mut ending_reads_flank = HashSet::new();
    let mut bam_reader = bam::IndexedReader::from_path(bam_name)?;
    let regions = region_coordinates.flanking_regions.unwrap_or_default();
    for (i, region) in regions.iter().take(2).enumerate() {
        let fields = region.split_terminator(':').collect::<Vec<_>>();
        let nchr = fields[0];
        let coords = fields[1].split_terminator('-').collect::<Vec<_>>();
        let start = coords[0].parse::<i64>()?;
        let stop = coords[1].parse::<i64>()?;

        let resolved_chr = resolve_chrom_name_from_header(bam_reader.header(), nchr)
            .ok_or_else(|| format!("Chromosome '{nchr}' not found in input BAM header"))?;
        trace!("fetching region: {resolved_chr:?} {start:?} {stop:?}");
        bam_reader
            .fetch((resolved_chr.as_str(), start, stop))
            .map_err(|e| format!("Failed to fetch region {resolved_chr}:{start}-{stop}: {e}"))?;

        for read in bam_reader.records() {
            let record = read?;
            let keep_record = interval_mismatch(&record, genome_reference, start, stop)?;
            if !record.is_secondary() && record.mapq() >= 30 && keep_record {
                let qname = std::str::from_utf8(record.qname())?.to_string();
                if i == 0 {
                    starting_reads_flank.insert(qname);
                } else {
                    ending_reads_flank.insert(qname);
                }
            }
        }
    }

    let flank_reads = FlankReads {
        start: starting_reads_flank,
        end: ending_reads_flank,
        start_segment: HashSet::new(),
        end_segment: HashSet::new(),
        good_clips_p5: vec![],
        good_clips_p3: vec![],
    };
    Ok(flank_reads)
}

/// For storing clipped reads
#[derive(Clone, Debug, Default)]
pub struct ClippedReads {
    /// positions with 5p clips, (tid, position)
    pub good_clips_p5: Vec<(usize, i64)>,
    /// positions with 3p clips, (tid, position)
    pub good_clips_p3: Vec<(usize, i64)>,
    /// clipped reads, (segment_name, (clip_type, tid, position))
    pub clipped_reads: BTreeMap<String, Vec<(String, usize, i64)>>,
}

/// Get clipped reads
/// # Arguments
/// * `repeat_records` - realigned repeat records
/// * `reference` - repeat reference file
/// # Returns
/// * `ClippedReads` - clipped reads
pub fn get_clipped_reads(
    repeat_records: &Vec<bam::Record>,
    reference: &PathBuf,
) -> Result<ClippedReads, DError> {
    let mut clipped_reads: BTreeMap<String, Vec<(String, usize, i64)>> = BTreeMap::new();
    let ref_reader = faidx::Reader::from_path(reference)?;
    for record in repeat_records {
        let tid = record.tid();
        let ref_name = ref_reader.seq_name(tid as i32)?;
        let ref_len = ref_reader.fetch_seq_len(&ref_name);
        let record_cigar = record.cigar();
        let reference_start_pos = record.pos();
        let reference_end_pos = record.reference_end();
        let alignment_len = reference_end_pos - reference_start_pos;
        let read_start_pos = start_pos_on_read(&record);
        let qname = std::str::from_utf8(record.qname())?.to_string();
        let segment_name = format!("{qname}:{}:{}", read_start_pos, alignment_len);
        // clip on p5
        if reference_start_pos >= 5 {
            let first_cigar = record_cigar.first().unwrap();
            match first_cigar {
                Cigar::SoftClip(softclip_len) => {
                    if *softclip_len > 100 {
                        clipped_reads
                            .entry(segment_name.clone())
                            .or_default()
                            .push((String::from("p5"), tid as usize, reference_start_pos));
                    }
                }
                _ => {}
            }
        }
        // clip on p3
        if reference_end_pos < ref_len as i64 - 5 {
            let last_cigar = record_cigar.last().unwrap();
            match last_cigar {
                Cigar::SoftClip(softclip_len) => {
                    if *softclip_len > 100 {
                        clipped_reads
                            .entry(segment_name.clone())
                            .or_default()
                            .push((String::from("p3"), tid as usize, reference_end_pos));
                    }
                }
                _ => {}
            }
        }
    }

    let mut clips_p5 = Vec::new();
    let mut clips_p3 = Vec::new();
    for (_segment, segment_clips) in &clipped_reads {
        for v in segment_clips {
            if v.0 == String::from("p5") {
                clips_p5.push((v.1, v.2));
            }
            if v.0 == String::from("p3") {
                clips_p3.push((v.1, v.2));
            }
        }
    }
    let good_clips_p5: Vec<(usize, i64)> = clips_p5
        .into_iter()
        .collect::<counter::Counter<(usize, i64)>>()
        .most_common_ordered()
        .into_iter()
        .filter(|x| x.1 >= 5)
        .map(|x| x.0)
        .collect::<Vec<_>>();
    let good_clips_p3: Vec<(usize, i64)> = clips_p3
        .into_iter()
        .collect::<counter::Counter<(usize, i64)>>()
        .most_common_ordered()
        .into_iter()
        .filter(|x| x.1 >= 5)
        .map(|x| x.0)
        .collect::<Vec<_>>();
    debug!("good_clips {good_clips_p5:?} {good_clips_p3:?}");
    debug!("clipped_reads {clipped_reads:?}");
    Ok(ClippedReads {
        good_clips_p5,
        good_clips_p3,
        clipped_reads,
    })
}

/// Get reads overlapping start and end of D4Z4
/// # Arguments
/// * `bam_name` - realigned bam
/// * `reference` - repeat reference file
/// * `read_length` - read lengths lookup
/// * `region_coordinates` - coordinates defined in the region
/// # Returns
/// * `FlankReads` - flanking reads
/// * `ClippedReads` - clipped reads
pub fn get_start_end_d4z4(
    bam_name: PathBuf,
    reference: &PathBuf,
    read_length: BTreeMap<String, usize>,
    region_coordinates: RegionCoordinates,
) -> Result<(FlankReads, ClippedReads), DError> {
    let mut clipped_reads: BTreeMap<String, Vec<(String, usize, i64)>> = BTreeMap::new();
    let starting_reads_flank = HashSet::new();
    let mut ending_reads_flank = HashSet::new();
    let mut starting_segments_flank = HashSet::new();
    let mut ending_segments_flank = HashSet::new();
    let mut bam_reader = bam::IndexedReader::from_path(bam_name)?;
    let start_positions_flank = region_coordinates.start_positions_flank.unwrap();
    let end_positions_flank = region_coordinates.end_positions_flank.unwrap();
    let ref_reader = faidx::Reader::from_path(reference)?;
    let ref_name = ref_reader.seq_name(0)?;

    bam_reader
        .fetch((&ref_name, 1, (region_coordinates.repeat_len as i64)))
        .map_err(|e| format!("Failed to fetch region {ref_name}:1-{}: {e}", region_coordinates.repeat_len))?;
    for read in bam_reader.records() {
        let record = read?;
        let record_cigar = record.cigar();
        let reference_start_pos = record.pos();
        let reference_end_pos = record.reference_end();
        let alignment_len = reference_end_pos - reference_start_pos;
        let read_start_pos = start_pos_on_read(&record);
        let qname = std::str::from_utf8(record.qname())?.to_string();
        let _this_read_length = *read_length
            .get(&qname)
            .ok_or("key not found in read_length")? as i32;
        let segment_name = format!("{qname}:{}:{}", read_start_pos, alignment_len);
        if (reference_end_pos as usize) > region_coordinates.repeat_len + 50 {
            ending_reads_flank.insert(qname);
            ending_reads_flank.insert(segment_name.clone());
        }
        // clip on p5
        if reference_start_pos >= 8
            && reference_start_pos < region_coordinates.repeat_len as i64 - 50
        {
            let first_cigar = record_cigar.first().unwrap();
            match first_cigar {
                Cigar::SoftClip(softclip_len) => {
                    if *softclip_len > 100 {
                        clipped_reads
                            .entry(segment_name.clone())
                            .or_default()
                            .push((String::from("p5"), 0, reference_start_pos));
                        for start_position in &start_positions_flank {
                            if reference_start_pos - start_position <= 2
                                && reference_start_pos - start_position >= -2
                            {
                                starting_segments_flank.insert(segment_name.clone());
                            }
                        }
                    }
                }
                _ => {}
            }
        }
        // clip on p3
        if reference_end_pos > 50 && reference_end_pos < region_coordinates.repeat_len as i64 - 50 {
            let last_cigar = record_cigar.last().unwrap();
            match last_cigar {
                Cigar::SoftClip(softclip_len) => {
                    if *softclip_len > 100 {
                        clipped_reads
                            .entry(segment_name.clone())
                            .or_default()
                            .push((String::from("p3"), 0, reference_end_pos));
                        for end_position in &end_positions_flank {
                            if reference_end_pos - end_position <= 2
                                && reference_end_pos - end_position >= -2
                            {
                                ending_segments_flank.insert(segment_name.clone());
                            }
                        }
                    }
                }
                _ => {}
            }
        }
    }
    let mut clips_p5 = Vec::new();
    let mut clips_p3 = Vec::new();
    for (_segment, segment_clips) in &clipped_reads {
        for v in segment_clips {
            if v.0 == String::from("p5") {
                clips_p5.push((v.1, v.2));
            }
            if v.0 == String::from("p3") {
                clips_p3.push((v.1, v.2));
            }
        }
    }
    let good_clips_p5: Vec<(usize, i64)> = clips_p5
        .into_iter()
        .collect::<counter::Counter<(usize, i64)>>()
        .most_common_ordered()
        .into_iter()
        .filter(|x| {
            (x.0 .1 > 50 && x.1 >= 3)
                || (x.0 .1 <= 50 && x.1 >= 15)
                || start_positions_flank.contains(&x.0 .1)
        })
        .map(|x| x.0)
        .collect::<Vec<_>>();
    let good_clips_p3: Vec<(usize, i64)> = clips_p3
        .into_iter()
        .collect::<counter::Counter<(usize, i64)>>()
        .most_common_ordered()
        .into_iter()
        .filter(|x| x.1 >= 3 || end_positions_flank.contains(&x.0 .1))
        .map(|x| x.0)
        .collect::<Vec<_>>();
    debug!("good_clips {good_clips_p5:?} {good_clips_p3:?}");
    debug!("clipped_reads {clipped_reads:?}");
    let flank_reads = FlankReads {
        start: starting_reads_flank,
        end: ending_reads_flank,
        start_segment: starting_segments_flank,
        end_segment: ending_segments_flank,
        good_clips_p5: vec![],
        good_clips_p3: vec![],
    };
    Ok((
        flank_reads,
        ClippedReads {
            good_clips_p5,
            good_clips_p3,
            clipped_reads,
        },
    ))
}

/// Tag reads with HP tag to indicate fingerprint assigment in the final output bam
/// # Arguments
/// * `repeat_records` - all read realignment records
/// * `grouped_reads` - read to fingerprint map
/// * `realigned_bam` - output bam
/// * `region_coordinates` - coordinates defined in the region
/// * `methyl_tags` - methyl tags lookup
/// * `reference` - repeat reference file
pub fn tag_reads(
    repeat_records: Vec<bam::Record>,
    grouped_reads: BTreeMap<String, i32>,
    realigned_bam: PathBuf,
    region_coordinates: RegionCoordinates,
    methyl_tags: BTreeMap<String, (String, Vec<u8>)>,
    reference: &PathBuf,
) -> DResult {
    let ref_reader = faidx::Reader::from_path(reference)?;
    let mut header = bam::Header::new();
    let chr_name = region_coordinates.chromosome_output;
    let mut record = HeaderRecord::new(b"SQ");
    record.push_tag(b"SN", chr_name.clone());
    record.push_tag(b"LN", region_coordinates.chromosome_len);
    header.push_record(&record);
    append_kivvi_pg_header(&mut header);

    let mut writer = Writer::from_path(&realigned_bam, &header, Format::Bam)?;
    for mut record in repeat_records {
        let qname = std::str::from_utf8(record.qname())?;
        let read_name = std::str::from_utf8(record.clone().qname())?.to_string();
        let tid = record.tid();
        let ref_name = ref_reader.seq_name(tid as i32)?;
        let this_offset = region_coordinates.genome_offset.get(&ref_name).unwrap();
        let read_start_pos = start_pos_on_read(&record);
        let reference_start_pos = &record.pos();
        let reference_end_pos = &record.reference_end();
        let new_name = format!("{qname}:{}", read_start_pos);
        if grouped_reads.contains_key(&new_name) {
            let hp_tag = grouped_reads.get(&new_name).ok_or("key not found")?;
            if hp_tag == &0 {
                if *reference_start_pos < 10
                    && *reference_end_pos > (region_coordinates.repeat_len as i64) - 10
                {
                    record.push_aux(b"HP", bam::record::Aux::String("unknown_full"))?;
                } else {
                    record.push_aux(b"HP", bam::record::Aux::String("unknown_partial"))?;
                }
            } else {
                record.push_aux(b"HP", bam::record::Aux::String(&hp_tag.to_string()))?;
            }
        } else {
            debug!("read {new_name} has no fp mapping...");
        }
        // offset position
        let genome_position = *reference_start_pos + 1 + *this_offset as i64;
        record.set_tid(0);
        record.set_pos(genome_position);
        // change SA tag
        match record.aux(b"SA") {
            Ok(value) => {
                if let Aux::String(sa_value) = value {
                    trace!("SA values {:?}", sa_value);
                    let mut new_sa_records: Vec<String> = Vec::new();
                    let sa_records = sa_value.split(";");
                    for sa_record in sa_records {
                        let parts = sa_record.split(",").collect::<Vec<&str>>();
                        if parts.len() == 6 {
                            let mut new_sa_record: Vec<String> = Vec::new();
                            new_sa_record.push(chr_name.clone());
                            let this_sa_chr = parts[0];
                            let this_sa_offset =
                                region_coordinates.genome_offset.get(this_sa_chr).unwrap();
                            let this_sa_position =
                                parts[1].parse::<i64>()? + 1 + *this_sa_offset as i64;
                            new_sa_record.push(this_sa_position.to_string());
                            new_sa_record.push(parts[2].to_string());
                            new_sa_record.push(parts[3].to_string());
                            new_sa_record.push(parts[4].to_string());
                            new_sa_record.push(parts[5].to_string());
                            new_sa_records.push(new_sa_record.join(","));
                        }
                    }
                    let mut new_sa_joined = new_sa_records.join(";");
                    new_sa_joined.push(';');
                    trace!("new SA {:?}", new_sa_joined);
                    let _ = record.remove_aux(b"SA");
                    record.push_aux(b"SA", Aux::String(&new_sa_joined))?;
                }
            }
            Err(_e) => {}
        };
        // add ML/MM tags
        if !methyl_tags.is_empty() && methyl_tags.contains_key(&read_name)
        //&& !record.is_supplementary()
        //&& !record.is_secondary()
        {
            let (mm, ml) = methyl_tags
                .get(&read_name)
                .ok_or("qname not in methyl_tags")?;
            record.push_aux(b"Mm", Aux::String(mm))?;
            record.push_aux(b"Ml", Aux::ArrayU8(ml.into()))?;
        }
        writer.write(&record)?;
    }
    drop(writer);
    bam::index::build(&realigned_bam, None, bam::index::Type::Bai, 1)?;
    Ok(())
}

/// Reverse complement a sequence
/// # Arguments
/// * `seq` - sequence
/// * `qual` - quality scores
/// # Returns
/// * `Vec<u8>` - reverse complemented sequence
fn reverse_complement(seq: &[u8], qual: &mut [u8]) -> Vec<u8> {
    let seq1 = seq
        .iter()
        .rev()
        .map(
            |x| match x | 32 /* make lower-case so we have fewer values to switch on */ {
                        b'a' => b'T',
                        b'c' => b'G',
                        b'g' => b'C',
                        b't' => b'A',
                        _ => b'?',
                    },
        )
        .collect::<Vec<u8>>();
    qual.reverse();
    return seq1;
}
