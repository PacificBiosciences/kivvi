use crate::bam_operation::start_pos_on_read;
use crate::util::DError;
use log::{debug, trace};
use rust_htslib::bam::{self, ext::BamRecordExtensions, record::Cigar, Writer};
use rust_htslib::faidx;
use std::collections::HashSet;
use std::path::PathBuf;

// deletion start position, position padding, deletion length, length padding
type DeletionFilter = (i64, i64, i64, i64);

/// Filter alignments for KIV2
/// # Arguments
/// * `realn_records` - realigned records
/// * `writer` - writer
/// * `reference` - reference
/// * `realigned_bam` - realigned bam
/// # Returns
/// * `repeat_records` - filtered records
pub fn filter_realignments_kiv2(
    realn_records: Vec<bam::Record>,
    mut writer: Writer,
    reference: &PathBuf,
    realigned_bam: PathBuf,
) -> Result<Vec<bam::Record>, DError> {
    // reference
    let ref_reader = faidx::Reader::from_path(reference)?;
    let mut repeat_records = Vec::<bam::Record>::new();
    for realn_record in &realn_records {
        //writer.write(realn_record)?;

        let reference_start_pos = realn_record.pos();
        let reference_end_pos = realn_record.reference_end();
        let alignment_len = reference_end_pos - reference_start_pos;
        let keep_read = interval_mismatch_kiv2(realn_record, &ref_reader)?;
        if realn_record.pos() + 1 < 3100
            && !realn_record.is_secondary()
            && alignment_len > 1500
            && realn_record.mapq() > 0
            && keep_read
        {
            writer.write(realn_record)?;
            repeat_records.push(realn_record.clone());
        }
    }
    drop(writer);
    bam::index::build(&realigned_bam, None, bam::index::Type::Bai, 1)?;
    Ok(repeat_records)
}

/// Filter alignments for HRNR
/// # Arguments
/// * `realn_records` - realigned records
/// * `writer` - writer
/// * `reference` - reference
/// * `realigned_bam` - realigned bam
/// # Returns
/// * `repeat_records` - filtered records
pub fn filter_realignments_hrnr(
    realn_records: Vec<bam::Record>,
    mut writer: Writer,
    reference: &PathBuf,
    realigned_bam: PathBuf,
) -> Result<Vec<bam::Record>, DError> {
    // reference
    let ref_reader = faidx::Reader::from_path(reference)?;
    let mut repeat_records = Vec::<bam::Record>::new();
    for realn_record in &realn_records {
        //writer.write(realn_record)?;
        let reference_start_pos = realn_record.pos();
        let reference_end_pos = realn_record.reference_end();
        let alignment_len = reference_end_pos - reference_start_pos;
        let qname = std::str::from_utf8(realn_record.qname())?.to_string();
        let nm = realn_record.aux(b"NM");
        if let Err(_e) = nm {
            debug!("missing NM tag for read {qname}");
        } else {
            let nm = nm.expect("expect NM tag");
            let mut nm = i32::try_from(extract_int_tag(&nm).expect("Tag was not integral."))
                .expect("Could not store nm in i32");
            let (longest_insertion_length, longest_deletion_length) =
                get_longest_insertion_deletion(realn_record)?;
            if longest_insertion_length > 20 && nm > longest_insertion_length as i32 {
                nm -= longest_insertion_length as i32;
            }
            if longest_deletion_length > 20 && nm > longest_deletion_length as i32 {
                nm -= longest_deletion_length as i32;
            }
            let keep_record = interval_mismatch_p5_and_p3(realn_record, &ref_reader, 450, 350)?;
            if alignment_len > 650
                && realn_record.mapq() >= 10
                && !realn_record.is_secondary()
                //&& (nm as f64) < (alignment_len as f64) * 0.05
                && ((nm as f64) < (alignment_len as f64) * 0.05 || (keep_record && alignment_len > 680))
            {
                writer.write(realn_record)?;
                repeat_records.push(realn_record.clone());
            }
        }
    }
    drop(writer);
    bam::index::build(&realigned_bam, None, bam::index::Type::Bai, 1)?;
    Ok(repeat_records)
}

/// Filter alignments for NBPF
/// # Arguments
/// * `realn_records` - realigned records
/// * `writer` - writer
/// * `reference` - reference
/// * `realigned_bam` - realigned bam
/// # Returns
/// * `repeat_records` - filtered records
pub fn filter_realignments_nbpf(
    realn_records: Vec<bam::Record>,
    mut writer: Writer,
    reference: &PathBuf,
    realigned_bam: PathBuf,
) -> Result<Vec<bam::Record>, DError> {
    // reference
    let ref_reader = faidx::Reader::from_path(reference)?;
    let mut repeat_records = Vec::<bam::Record>::new();
    for realn_record in &realn_records {
        //writer.write(realn_record)?;
        let reference_start_pos = realn_record.pos();
        let reference_end_pos = realn_record.reference_end();
        let alignment_len = reference_end_pos - reference_start_pos;
        let qname = std::str::from_utf8(realn_record.qname())?.to_string();
        let nm = realn_record.aux(b"NM");
        if let Err(_e) = nm {
            debug!("missing NM tag for read {qname}");
        } else {
            let nm = nm.expect("expect NM tag");
            let mut nm = i32::try_from(extract_int_tag(&nm).expect("Tag was not integral."))
                .expect("Could not store nm in i32");
            let (longest_insertion_length, longest_deletion_length) =
                get_longest_insertion_deletion(realn_record)?;
            if longest_insertion_length > 120 && nm > longest_insertion_length as i32 {
                nm -= longest_insertion_length as i32;
            }
            if longest_deletion_length > 120 && nm > longest_deletion_length as i32 {
                nm -= longest_deletion_length as i32;
            }
            let keep_record = interval_mismatch_p5_and_p3(realn_record, &ref_reader, 1300, 300)?;
            if alignment_len > 1400
                && realn_record.mapq() > 0
                && !realn_record.is_secondary()
                //&& (nm as f64) < (alignment_len as f64) * 0.02
            && ((nm as f64) < (alignment_len as f64) * 0.02 || (keep_record && alignment_len > 1300))
            {
                writer.write(realn_record)?;
                repeat_records.push(realn_record.clone());
            }
        }
    }
    drop(writer);
    bam::index::build(&realigned_bam, None, bam::index::Type::Bai, 1)?;
    Ok(repeat_records)
}

/// Filter alignments for D4Z4
/// # Arguments
/// * `realn_records` - realigned records
/// * `writer` - writer
/// * `reference` - reference
/// * `realigned_bam` - realigned bam
/// # Returns
/// * `repeat_records` - filtered records
/// * `white_list_read_segments` - white list read segments
pub fn filter_realignments_d4z4(
    realn_records: Vec<bam::Record>,
    mut writer: Writer,
    _reference: &PathBuf,
    realigned_bam: PathBuf,
) -> Result<(Vec<bam::Record>, Vec<String>, Vec<String>), DError> {
    // repeat units with big insertions may sometimes be aligned as extra segment with deletions
    // we don't want these segments to become a different fingerprint
    const D4Z4_BLACKLIST_DELETIONS: [DeletionFilter; 2] = [(1574, 5, 1685, 20), (2822, 5, 324, 10)];

    let mut white_list_read_segments = Vec::new();
    let mut blacklist_segments = Vec::new();
    let mut records_to_keep = HashSet::new();
    let mut repeat_records = Vec::<bam::Record>::new();
    for realn_record in &realn_records {
        let reference_start_pos = realn_record.pos();
        let reference_end_pos = realn_record.reference_end();
        let alignment_len = reference_end_pos - reference_start_pos;
        let read_start_pos = start_pos_on_read(realn_record);
        let qname = std::str::from_utf8(realn_record.qname())?.to_string();
        let nm = realn_record.aux(b"NM");
        if let Err(_e) = nm {
            debug!("missing NM tag for read {qname}");
        } else {
            let nm = nm.expect("expect NM tag");
            let mut nm = i32::try_from(extract_int_tag(&nm).expect("Tag was not integral."))
                .expect("Could not store nm in i32");
            let (longest_insertion_length, longest_deletion_length) =
                get_longest_insertion_deletion(realn_record)?;
            if longest_insertion_length > 150 && nm > longest_insertion_length as i32 {
                nm -= longest_insertion_length as i32;
            }
            if longest_deletion_length > 150 && nm > longest_deletion_length as i32 {
                nm -= longest_deletion_length as i32;
            }
            let segment_name = format!("{qname}:{}:{}", read_start_pos, alignment_len);
            if longest_insertion_length > 1600 && longest_insertion_length < 1650 {
                white_list_read_segments.push(segment_name.clone());
            }
            //else if reference_end_pos > ref_len as i64 - 5 {
            //    white_list_read_segments.push(segment_name.clone());
            //}
            let mut keep_record = false;
            if !realn_record.is_secondary() {
                if alignment_len > 1000 && (nm as f64) < (alignment_len as f64) * 0.07 && nm < 200 {
                    keep_record = true;
                } else if alignment_len > 400 && (nm as f64) < (alignment_len as f64) * 0.03 {
                    keep_record = true;
                }
            }
            if keep_record
                && has_deletion_near_any_position(realn_record, &D4Z4_BLACKLIST_DELETIONS)
            {
                let cigar_string = realn_record.cigar().to_string();
                debug!(
                    "Found read segment with blacklist deletion: {segment_name} cigar {cigar_string}"
                );
                blacklist_segments.push(segment_name.clone());
            }
            if keep_record {
                //writer.write(realn_record)?;
                //repeat_records.push(realn_record.clone());
                records_to_keep.insert(qname.clone());
            } else {
                trace!(
                    "Filtering out read {qname}_{read_start_pos} alignment length {alignment_len} mismatch {nm}"
                );
            }
        }
    }
    for realn_record in &realn_records {
        let qname = std::str::from_utf8(realn_record.qname())?.to_string();
        if records_to_keep.contains(&qname) {
            writer.write(realn_record)?;
            repeat_records.push(realn_record.clone());
        }
    }
    drop(writer);
    bam::index::build(&realigned_bam, None, bam::index::Type::Bai, 1)?;
    Ok((repeat_records, white_list_read_segments, blacklist_segments))
}

fn has_deletion_near_any_position(
    record: &bam::Record,
    deletion_filters: &[DeletionFilter],
) -> bool {
    deletion_filters
        .iter()
        .any(|&(target_pos, pos_padding, target_len, len_padding)| {
            has_deletion_near_position(record, target_pos, pos_padding, target_len, len_padding)
        })
}

fn has_deletion_near_position(
    record: &bam::Record,
    target_pos: i64,
    pos_padding: i64,
    target_len: i64,
    len_padding: i64,
) -> bool {
    let mut ref_pos = record.pos();
    for cigar in &record.cigar() {
        match cigar {
            Cigar::Del(len) => {
                let deletion_start = ref_pos + 1;
                let deletion_len = i64::from(*len);
                if (deletion_start - target_pos).abs() <= pos_padding
                    && (deletion_len - target_len).abs() <= len_padding
                {
                    return true;
                }
                ref_pos += deletion_len;
            }
            Cigar::Match(len) | Cigar::Equal(len) | Cigar::Diff(len) | Cigar::RefSkip(len) => {
                ref_pos += i64::from(*len);
            }
            Cigar::Ins(_) | Cigar::SoftClip(_) | Cigar::HardClip(_) | Cigar::Pad(_) => {}
        }
    }
    false
}

/// Filter alignments based on a more detailed mismatch calculation, used for hrnr and nbpf
/// Allows a read with high identity at the 5' end or 3' end to be kept even if the overall identity is low
/// # Arguments
/// * `record` - bam record
/// * `ref_reader` - reference reader
/// * `p5_end` - 5 prime end coordinate
/// * `p3_start` - 3 prime start coordinate
/// # Returns
/// * `bool` - true if the record fulfils the interval mismatch criteria
fn interval_mismatch_p5_and_p3(
    record: &bam::Record,
    ref_reader: &faidx::Reader,
    p5_end: i64,
    p3_start: i64,
) -> Result<bool, DError> {
    //let qname = std::str::from_utf8(record.qname())?;
    let mut region_match_3p = 0;
    let mut new_nm_3p = 0;
    let mut region_match_5p = 0;
    let mut new_nm_5p = 0;
    let seq = record.seq().as_bytes();
    let tid = record.tid();
    let ref_name = ref_reader.seq_name(tid)?;

    for [read_pos, ref_pos] in record.aligned_pairs() {
        if ref_pos + 1 < p5_end {
            let read_pos = read_pos as usize;
            let ref_pos = ref_pos as usize;
            if let Some(read_base) = seq.get(read_pos) {
                let ref_base = ref_reader.fetch_seq(&ref_name, ref_pos, ref_pos)?;
                region_match_5p += 1;
                if read_base != ref_base.first().unwrap() {
                    new_nm_5p += 1;
                }
            }
        }
        if ref_pos + 1 > p3_start {
            let read_pos = read_pos as usize;
            let ref_pos = ref_pos as usize;
            if let Some(read_base) = seq.get(read_pos) {
                let ref_base = ref_reader.fetch_seq(&ref_name, ref_pos, ref_pos)?;
                region_match_3p += 1;
                if read_base != ref_base.first().unwrap() {
                    new_nm_3p += 1;
                }
            }
        }
    }
    if (new_nm_5p as f64) < (region_match_5p as f64) * 0.005
        || (new_nm_3p as f64) < (region_match_3p as f64) * 0.001
    {
        Ok(true)
    } else {
        Ok(false)
    }
}

/// Filter alignments based on a more detailed mismatch calculation for KIV2
/// # Arguments
/// * `record` - bam record
/// * `ref_reader` - reference reader
/// # Returns
/// * `bool` - true if the record fulfils the interval mismatch criteria
fn interval_mismatch_kiv2(
    record: &bam::Record,
    ref_reader: &faidx::Reader,
) -> Result<bool, DError> {
    let qname = std::str::from_utf8(record.qname())?;
    let mut region_match = 0;
    let mut new_nm = 0;
    let mut region_match_5p = 0;
    let mut new_nm_5p = 0;
    let seq = record.seq().as_bytes();
    let ref_name = ref_reader.seq_name(0)?;

    for [read_pos, ref_pos] in record.aligned_pairs() {
        if ref_pos + 1 < 3560 || ref_pos + 1 > 4760 {
            let read_pos = read_pos as usize;
            let ref_pos = ref_pos as usize;
            if let Some(read_base) = seq.get(read_pos) {
                let ref_base = ref_reader.fetch_seq(&ref_name, ref_pos, ref_pos)?;
                //debug!("read_pos {read_pos:?} ref_pos {ref_pos:?} read_base {read_base:?} ref_base {ref_base:?}");
                region_match += 1;
                region_match_5p += 1;
                if read_base != ref_base.first().unwrap() {
                    new_nm += 1;
                    if ref_pos + 1 < 3560 {
                        new_nm_5p += 1;
                    }
                }
            }
        }
    }
    trace!("{qname:?} region_match {region_match:?} new_nm {new_nm:?}");
    if (new_nm as f64) < (region_match as f64) * 0.02
        && (new_nm_5p as f64) < (region_match_5p as f64) * 0.02
    {
        Ok(true)
    } else {
        Ok(false)
    }
}

/// Get the longest insertion and deletion length from a bam record
/// # Arguments
/// * `record` - bam record
/// # Returns
/// * `insertion_length` - longest insertion length
/// * `deletion_length` - longest deletion length
#[must_use]
pub fn get_longest_insertion_deletion(record: &bam::Record) -> Result<(i64, i64), DError> {
    let mut insertion_lengths = Vec::new();
    let mut deletion_lengths = Vec::new();
    for x in record.cigar().iter() {
        match x {
            Cigar::Ins(_len) => insertion_lengths.push(i64::from(x.len())),
            Cigar::Del(_len) => deletion_lengths.push(i64::from(x.len())),
            _ => {}
        }
    }
    let ins_length = if insertion_lengths.is_empty() {
        0
    } else {
        *insertion_lengths
            .iter()
            .max()
            .ok_or("max not found in insertion_lengths")?
    };
    let del_length = if deletion_lengths.is_empty() {
        0
    } else {
        *deletion_lengths
            .iter()
            .max()
            .ok_or("max not found in deletion_lengths")?
    };
    Ok((ins_length, del_length))
}

/// Extract an integer tag from a bam record
/// # Arguments
/// * `tag` - tag
/// # Returns
/// * `int_tag` - integer tag
#[must_use]
fn extract_int_tag(tag: &bam::record::Aux) -> Option<i64> {
    match tag {
        rust_htslib::bam::record::Aux::I8(tag) => Some(i64::from(*tag)),
        rust_htslib::bam::record::Aux::I16(tag) => Some(i64::from(*tag)),
        rust_htslib::bam::record::Aux::I32(tag) => Some(i64::from(*tag)),
        rust_htslib::bam::record::Aux::U8(tag) => Some(i64::from(*tag)),
        rust_htslib::bam::record::Aux::U16(tag) => Some(i64::from(*tag)),
        rust_htslib::bam::record::Aux::U32(tag) => Some(i64::from(*tag)),
        _ => None,
    }
}
