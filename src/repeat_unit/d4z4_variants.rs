use crate::bam_operation::start_pos_on_read;
use crate::util::RegionCoordinates;
use crate::util::{invalid_data_error, DError};
use log::{debug, error, trace};
use rust_htslib::bam::ext::BamRecordExtensions;
use rust_htslib::{bam, bam::Read, faidx};
use std::collections::{BTreeMap, HashSet};
use std::path::PathBuf;

/// Update the read with special calls
/// # Arguments
/// * `read_segment_raw_fp` - read segment -> raw fps
/// * `new_variants_by_position` - retained variants grouped by position
/// * `realigned_bam` - realigned bam file
/// * `reference` - reference file
/// * `region_coordinates` - region coordinates
/// # Returns
/// * `(BTreeMap<String, Vec<u8>>, BTreeMap<i64, Vec<crate::realignment::utilities::Variant>>)` - read segment -> updated fps and retained variants grouped by position
pub fn update_read_with_special_calls(
    read_segment_raw_fp: &BTreeMap<String, Vec<u8>>,
    new_variants_by_position: &BTreeMap<i64, Vec<crate::realignment::utilities::Variant>>,
    realigned_bam: PathBuf,
    reference: &PathBuf,
    region_coordinates: &RegionCoordinates,
) -> Result<
    (
        BTreeMap<String, Vec<u8>>,
        BTreeMap<i64, Vec<crate::realignment::utilities::Variant>>,
    ),
    DError,
> {
    let (special_calls_homopolymer, success_homopolymer) =
        genotype_homopolymer(realigned_bam.clone(), reference, region_coordinates)?;
    let (special_calls_str, success_str) =
        genotype_str(realigned_bam.clone(), reference, region_coordinates)?;

    let mut read_segment_raw_fp_updated = BTreeMap::new();
    for (segment_name, fp) in read_segment_raw_fp {
        let this_call_homopolymer = special_calls_homopolymer
            .get(segment_name)
            .copied()
            .unwrap_or_else(|| {
                error!("Segment {segment_name} missing homopolymer special call");
                b'-'
            });
        let this_call_str = special_calls_str
            .get(segment_name)
            .copied()
            .unwrap_or_else(|| {
                error!("Segment {segment_name} missing str special call");
                b'-'
            });
        let mut new_fp = fp.clone();
        if success_str {
            if fp.starts_with(&[b'S']) {
                new_fp.insert(0, b'S');
            } else {
                new_fp.insert(0, this_call_str);
            }
        }
        if success_homopolymer {
            if fp.ends_with(&[b'S']) {
                new_fp.push(b'S');
            } else {
                new_fp.push(this_call_homopolymer);
            }
        }
        read_segment_raw_fp_updated.insert(segment_name.clone(), new_fp);
    }
    let mut new_variants_by_position_updated = new_variants_by_position.clone();
    if success_str {
        new_variants_by_position_updated.entry(0).or_insert(vec![]);
    }
    if success_homopolymer {
        new_variants_by_position_updated
            .entry(5000)
            .or_insert(vec![]);
    }
    Ok((
        read_segment_raw_fp_updated,
        new_variants_by_position_updated,
    ))
}

/// Genotype the homopolymer region
/// # Arguments
/// * `realigned_bam` - realigned bam file
/// * `reference` - reference file
/// * `region_coordinates` - region coordinates
/// # Returns
/// * `BTreeMap<String, u8>` - read segment -> call
/// * `bool` - success
fn genotype_homopolymer(
    realigned_bam: PathBuf,
    reference: &PathBuf,
    region_coordinates: &RegionCoordinates,
) -> Result<(BTreeMap<String, u8>, bool), DError> {
    debug!("Genotype the homopolymer region at position 3113...");
    let mut success = false;
    let mut expected_variant: BTreeMap<String, u8> = BTreeMap::new();
    expected_variant.insert(String::from("A"), b'0');
    expected_variant.insert(String::from("G"), b'1');
    expected_variant.insert(String::from("GT"), b'2');
    expected_variant.insert(String::from("GA"), b'3');
    expected_variant.insert(String::from("AT"), b'4');
    let mut map_segment_to_call: BTreeMap<String, u8> = BTreeMap::new();
    let ref_reader = faidx::Reader::from_path(reference)?;
    let ref_name = ref_reader.seq_name(0)?;
    let mut bam_reader = bam::IndexedReader::from_path(realigned_bam.clone())?;
    bam_reader
        .fetch((&ref_name, 0, region_coordinates.repeat_len as i64))
        .map_err(|e| {
            invalid_data_error(format!(
                "failed to fetch homopolymer genotyping region 0-{} on {ref_name}: {e}",
                region_coordinates.repeat_len
            ))
        })?;
    for read_entry in bam_reader.records() {
        let read = read_entry?;
        let qname = std::str::from_utf8(read.qname())?;
        let read_start_pos = start_pos_on_read(&read);
        let reference_start_pos = &read.pos();
        let reference_end_pos = &read.reference_end();
        let aln_len = reference_end_pos - reference_start_pos;
        let segment_name = format!("{qname}:{}:{}", read_start_pos, aln_len);
        if *reference_start_pos > 3110 || *reference_end_pos < 3115 {
            map_segment_to_call.insert(segment_name.clone(), b'x');
        } else {
            map_segment_to_call.insert(segment_name.clone(), b'-');
        }

        let mut read_start: Option<i64> = None;
        let mut read_end: Option<i64> = None;
        for bp in read.aligned_pairs() {
            let segment_index = bp[0];
            let ref_index = bp[1];
            if ref_index == 3110 {
                read_start = Some(segment_index);
            }
            if ref_index == 3115 {
                read_end = Some(segment_index);
            }
            if read_start.is_some() && read_end.is_some() {
                break;
            }
        }
        if let (Some(read_start), Some(read_end)) = (read_start, read_end) {
            let read_start = read_start as usize;
            let read_end = read_end as usize;
            let read_seq = read.seq().as_bytes();
            let read_seq = std::str::from_utf8(&read_seq[read_start..read_end])?;
            let read_seq_strip_c = read_seq.trim_start_matches("C").trim_end_matches("C");
            trace!(
                "{segment_name}, read_seq {:?}, read_seq_strip_c: {:?}",
                read_seq,
                read_seq_strip_c
            );
            let this_call = expected_variant.get(read_seq_strip_c).unwrap_or(&b'-');
            map_segment_to_call.insert(segment_name, *this_call);
        }
    }
    let mut calls_to_remove = HashSet::new();
    let mut base_counts = BTreeMap::<u8, i32>::new();
    for (_segment_name, call) in &map_segment_to_call {
        *base_counts.entry(*call).or_default() += 1;
    }
    let mut count_missing_bases = 0;
    let mut total_bases = 0;
    for (base, count) in base_counts {
        debug!("Base {:?} has count {count}", std::str::from_utf8(&[base])?);
        if base == b'-' {
            count_missing_bases += count;
        }
        if base != b'x' {
            total_bases += count;
        }
        if count <= 2 && base != b'-' && base != b'x' {
            debug!(
                "Removing call {:?} with count {count}",
                std::str::from_utf8(&[base])?
            );
            calls_to_remove.insert(base);
        }
    }
    let map_segment_to_call_clone = map_segment_to_call.clone();
    for (segment_name, call) in &map_segment_to_call_clone {
        if calls_to_remove.contains(&call) {
            debug!("Updating call {call} to unknown for segment {segment_name}");
            map_segment_to_call.insert(segment_name.clone(), b'-');
            count_missing_bases += 1;
        }
    }

    debug!("homopolymer region count_missing_bases {count_missing_bases}");
    if count_missing_bases as f64 <= 20.0_f64.max(total_bases as f64 * 0.025) {
        success = true;
    } else {
        debug!("Failed to genotype the homopolymer region");
    }

    Ok((map_segment_to_call, success))
}

/// Genotype the STR region
/// # Arguments
/// * `realigned_bam` - realigned bam file
/// * `reference` - reference file
/// * `region_coordinates` - region coordinates
/// # Returns
/// * `BTreeMap<String, u8>` - read segment -> call
/// * `bool` - success
fn genotype_str(
    realigned_bam: PathBuf,
    reference: &PathBuf,
    region_coordinates: &RegionCoordinates,
) -> Result<(BTreeMap<String, u8>, bool), DError> {
    debug!("Genotype the STR region around position 138...");
    let mut success = false;
    let mut expected_variant: BTreeMap<String, u8> = BTreeMap::new();
    expected_variant.insert(String::from("2_7_36"), b'0');
    expected_variant.insert(String::from("2_7_35"), b'0');
    expected_variant.insert(String::from("2_7_37"), b'0');

    expected_variant.insert(String::from("3_6_36"), b'1');
    expected_variant.insert(String::from("3_6_35"), b'1');
    expected_variant.insert(String::from("3_6_37"), b'1');

    expected_variant.insert(String::from("1_8_36"), b'2');
    expected_variant.insert(String::from("1_8_35"), b'2');
    expected_variant.insert(String::from("1_8_37"), b'2');

    expected_variant.insert(String::from("2_6_36"), b'3');
    expected_variant.insert(String::from("2_6_35"), b'3');
    expected_variant.insert(String::from("2_6_37"), b'3');

    expected_variant.insert(String::from("3_8_44"), b'4');
    expected_variant.insert(String::from("3_8_43"), b'4');
    expected_variant.insert(String::from("3_8_45"), b'4');
    expected_variant.insert(String::from("3_8_42"), b'4');
    expected_variant.insert(String::from("3_8_46"), b'4');

    expected_variant.insert(String::from("2_5_28"), b'5');
    expected_variant.insert(String::from("2_5_27"), b'5');
    expected_variant.insert(String::from("2_5_29"), b'5');

    expected_variant.insert(String::from("1_6_28"), b'6');
    expected_variant.insert(String::from("1_6_27"), b'6');
    expected_variant.insert(String::from("1_6_29"), b'6');

    expected_variant.insert(String::from("2_6_32"), b'7');
    expected_variant.insert(String::from("2_6_31"), b'7');
    expected_variant.insert(String::from("2_6_33"), b'7');

    expected_variant.insert(String::from("3_7_40"), b'8');
    expected_variant.insert(String::from("3_7_39"), b'8');

    expected_variant.insert(String::from("4_8_48"), b'9');
    expected_variant.insert(String::from("4_8_47"), b'9');
    expected_variant.insert(String::from("4_8_49"), b'9');

    expected_variant.insert(String::from("4_9_52"), b'a');
    expected_variant.insert(String::from("4_9_51"), b'a');
    expected_variant.insert(String::from("4_9_53"), b'a');

    expected_variant.insert(String::from("5_9_56"), b'b');
    expected_variant.insert(String::from("5_9_55"), b'b');
    expected_variant.insert(String::from("5_9_57"), b'b');

    let mut map_segment_to_call: BTreeMap<String, u8> = BTreeMap::new();
    let ref_reader = faidx::Reader::from_path(reference)?;
    let ref_name = ref_reader.seq_name(0)?;
    let mut bam_reader = bam::IndexedReader::from_path(realigned_bam.clone())?;
    bam_reader
        .fetch((&ref_name, 0, region_coordinates.repeat_len as i64))
        .map_err(|e| {
            invalid_data_error(format!(
                "failed to fetch STR genotyping region 0-{} on {ref_name}: {e}",
                region_coordinates.repeat_len
            ))
        })?;
    for read_entry in bam_reader.records() {
        let read = read_entry?;
        let qname = std::str::from_utf8(read.qname())?;
        let read_start_pos = start_pos_on_read(&read);
        let reference_start_pos = &read.pos();
        let reference_end_pos = &read.reference_end();
        let aln_len = reference_end_pos - reference_start_pos;
        let segment_name = format!("{qname}:{}:{}", read_start_pos, aln_len);
        if *reference_start_pos > 126 || *reference_end_pos < 162 {
            map_segment_to_call.insert(segment_name.clone(), b'x');
        } else {
            map_segment_to_call.insert(segment_name.clone(), b'-');
        }

        let mut read_start: Option<i64> = None;
        let mut read_end: Option<i64> = None;
        let mut downstream_end: Option<i64> = None;
        for bp in read.aligned_pairs() {
            let segment_index = bp[0];
            let ref_index = bp[1];
            if ref_index == 126 {
                read_start = Some(segment_index);
            }
            if ref_index == 162 {
                read_end = Some(segment_index);
            }
            if ref_index == 177 {
                downstream_end = Some(segment_index);
            }
            if read_start.is_some() && read_end.is_some() && downstream_end.is_some() {
                break;
            }
        }
        if let (Some(read_start), Some(read_end), Some(downstream_end)) =
            (read_start, read_end, downstream_end)
        {
            let read_start = read_start as usize;
            let read_end = read_end as usize;
            let downstream_end = downstream_end as usize;
            let read_seq = read.seq().as_bytes();
            // check region downstream for indels
            let downstream_len_read = downstream_end as i32 - read_end as i32;
            let downstream_len_ref = 177 - 162;
            if downstream_len_read - downstream_len_ref <= -2
                || downstream_len_read - downstream_len_ref >= 4
            {
                debug!("segment {segment_name} has indels in the downstream region. downstream_len_read {downstream_len_read} downstream_len_ref {downstream_len_ref}");
            } else {
                let read_seq = std::str::from_utf8(&read_seq[read_start..read_end])?;
                let count_ca = read_seq.matches("CA").count();
                let count_ga = read_seq.matches("GA").count();
                let seq_len = read_seq.len();
                let expression = format!("{count_ca}_{count_ga}_{seq_len}");
                let this_call = expected_variant.get(&expression).unwrap_or(&b'-');
                trace!(
                    "{segment_name}, read_seq {:?}, expression: {:?}",
                    read_seq,
                    expression
                );
                map_segment_to_call.insert(segment_name, *this_call);
            }
        }
    }
    let mut calls_to_remove = HashSet::new();
    let mut base_counts = BTreeMap::<u8, i32>::new();
    for (_segment_name, call) in &map_segment_to_call {
        *base_counts.entry(*call).or_default() += 1;
    }
    let mut count_missing_bases = 0;
    let mut total_bases = 0;
    for (base, count) in base_counts {
        debug!("Base {:?} has count {count}", std::str::from_utf8(&[base])?);
        if base == b'-' {
            count_missing_bases += count;
        }
        if base != b'x' {
            total_bases += count;
        }
        if count <= 3 && base != b'-' && base != b'x' {
            debug!(
                "Removing call {:?} with count {count}",
                std::str::from_utf8(&[base])?
            );
            calls_to_remove.insert(base);
        }
    }
    let map_segment_to_call_clone = map_segment_to_call.clone();
    for (segment_name, call) in &map_segment_to_call_clone {
        if calls_to_remove.contains(&call) {
            debug!("Updating call {call} to unknown for segment {segment_name}");
            map_segment_to_call.insert(segment_name.clone(), b'-');
            count_missing_bases += 1;
        }
    }

    debug!("str region count_missing_bases {count_missing_bases}");
    if count_missing_bases as f64 <= 20.0_f64.max(total_bases as f64 * 0.025) {
        success = true;
    } else {
        debug!("Failed to genotype the str region");
    }

    Ok((map_segment_to_call, success))
}
