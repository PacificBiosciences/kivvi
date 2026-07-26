use crate::bam_operation::ClippedReads;
use crate::realignment::utilities::{Variant, VariantType};
use crate::repeat_unit::fingerprint::{FingerprintInfo, ReadParameters};
use crate::util::{d4z4_coordinates, DError, FlankReads, RegionCoordinates};
use itertools::Itertools;
use log::{debug, trace};
use std::collections::BTreeMap;
use std::collections::HashSet;

/// Clean up segment raw fps
/// # Arguments
/// * `read_segment_raw_fp` - read segment raw fps
/// * `variant_calls` - variant calls
/// * `flanking_reads` - flanking reads
/// * `clipped_reads` - clipped reads
/// * `region_coordinates` - region coordinates
/// * `min_variant_support` - minimum variant support
/// # Returns
/// * `(BTreeMap<String, Vec<u8>>, BTreeMap<i64, Vec<Variant>>)` - cleaned up read segment raw fps and retained variants grouped by position
pub fn clean_up_segment_raw_fps(
    read_segment_raw_fp: &mut BTreeMap<String, Vec<u8>>,
    variant_calls: &Vec<Variant>,
    _flanking_reads: &FlankReads,
    clipped_reads: &ClippedReads,
    region_coordinates: &RegionCoordinates,
    min_variant_support: usize,
) -> Result<(BTreeMap<String, Vec<u8>>, BTreeMap<i64, Vec<Variant>>), DError> {
    let variant_positions = variant_calls
        .iter()
        .map(|x| x.position())
        .collect::<HashSet<i64>>()
        .iter()
        .sorted()
        .map(|x| *x)
        .collect::<Vec<i64>>();
    debug!("variant_positions {variant_positions:?}");

    // handle type2 sites
    let mut remove_type2_sites = false;
    let type2_sites = &region_coordinates.type2_sites;
    if !type2_sites.is_empty() {
        let s1: HashSet<i64> = variant_calls
            .iter()
            .map(|x| x.position())
            .collect::<HashSet<i64>>();
        let s2: HashSet<i64> = type2_sites.iter().cloned().map(|x| x - 1).collect();
        let s3: Vec<i64> = (&s1 - &s2).iter().cloned().collect::<Vec<_>>();
        if s1.len() - s3.len() > s2.len() - 10 {
            remove_type2_sites = true;
            debug!("sample has type2 sites...removing these sites...");
        }
    }

    let num_pos = variant_positions.len();
    let reads_clipped = &clipped_reads.clipped_reads;
    let good_clips_p5 = &clipped_reads.good_clips_p5;
    let good_clips_p3 = &clipped_reads.good_clips_p3;
    for (read_segment, raw_fp) in &mut *read_segment_raw_fp {
        assert_eq!(raw_fp.len(), num_pos);
        if reads_clipped.contains_key(read_segment) {
            let this_read_clips = reads_clipped.get(read_segment).unwrap();
            for (clip_side, tid, clip_pos) in this_read_clips {
                if clip_side == "p5" {
                    if good_clips_p5.contains(&(*tid, *clip_pos))
                        || good_clips_p5.contains(&(*tid, *clip_pos - 1))
                        || good_clips_p5.contains(&(*tid, *clip_pos + 1))
                    {
                        for (i, pos) in variant_positions.iter().enumerate() {
                            if *pos < *clip_pos {
                                raw_fp[i] = b'S';
                            }
                        }
                    }
                } else if clip_side == "p3" {
                    if good_clips_p3.contains(&(*tid, *clip_pos))
                        || good_clips_p3.contains(&(*tid, *clip_pos - 1))
                        || good_clips_p3.contains(&(*tid, *clip_pos + 1))
                    {
                        for (i, pos) in variant_positions.iter().enumerate() {
                            if *pos > *clip_pos {
                                raw_fp[i] = b'S';
                            }
                        }
                    }
                }
            }
        }
    }
    // filter out bad positions
    let mut new_variant_positions = Vec::new();
    for (i, pos) in variant_positions.iter().enumerate() {
        let mut is_good_pos = false;
        let this_position_bases = read_segment_raw_fp
            .values()
            .map(|x| x[i])
            .collect::<Vec<_>>();
        let count_missing_bases = this_position_bases.iter().filter(|x| **x == b'-').count();
        debug!("pos {pos} count_missing_bases {count_missing_bases}");
        if count_missing_bases as f64 <= 20.0_f64.max(this_position_bases.len() as f64 * 0.025) {
            let good_bases = this_position_bases
                .iter()
                .filter(|x| **x != b'-' && **x != b'x' && **x != b'S')
                .collect::<counter::Counter<_>>()
                .most_common_ordered();
            debug!("pos {pos} good_bases {good_bases:?}");
            if good_bases.len() > 1 {
                let num1 = good_bases[0].1;
                let num2 = good_bases[1].1;
                if num1 >= min_variant_support
                    && num2 >= min_variant_support
                    && (!remove_type2_sites
                        || !type2_sites.contains(&(pos + 1))
                        || (good_bases.len() > 2 && good_bases[2].1 >= min_variant_support))
                {
                    new_variant_positions.push(*pos);
                    is_good_pos = true;
                }
            }
        }
        if !is_good_pos {
            debug!("filter out bad position or type2 position {pos}");
        }
    }
    debug!("new_variant_positions {new_variant_positions:?}");
    let new_variant_position_set = new_variant_positions
        .iter()
        .copied()
        .collect::<HashSet<i64>>();
    let mut new_variants_by_position: BTreeMap<i64, Vec<Variant>> = BTreeMap::new();
    for variant in variant_calls
        .iter()
        .filter(|variant| new_variant_position_set.contains(&variant.position()))
    {
        new_variants_by_position
            .entry(variant.position())
            .or_default()
            .push(variant.clone());
    }
    let mut new_read_segment_raw_fp = BTreeMap::new();
    for (read_segment, raw_fp) in read_segment_raw_fp {
        let mut new_raw_fp = Vec::new();
        for (i, base) in raw_fp.iter().enumerate() {
            if new_variant_position_set.contains(&variant_positions[i]) {
                new_raw_fp.push(*base);
            }
        }
        trace!(
            "cleaned raw_fp {:?} to {:?}",
            std::str::from_utf8(&raw_fp.clone())?,
            std::str::from_utf8(&new_raw_fp.clone())?
        );
        new_read_segment_raw_fp.insert(read_segment.clone(), new_raw_fp);
    }
    Ok((new_read_segment_raw_fp, new_variants_by_position))
}

/// Get good variants from unfiltered sites
/// # Arguments
/// * `region_coordinates` - region coordinates
/// * `unfiltered_sites` - unfiltered sites
/// * `ref_seq` - reference sequence
/// * `is_d4z4` - whether the region is d4z4
/// # Returns
/// * `Vec<Variant>` - good variants
pub fn get_good_variants(
    region_coordinates: &RegionCoordinates,
    unfiltered_sites: &BTreeMap<i64, Vec<Vec<u8>>>,
    ref_seq: &Vec<u8>,
    _is_d4z4: bool,
) -> Result<Vec<Variant>, DError> {
    let mut variants_to_force_call = region_coordinates.variants_to_call.clone();
    let known_del_sites = variants_to_force_call
        .iter()
        .filter(|x| x.get_type() == VariantType::Deletion)
        .map(|x| (x.position(), x.get_ref_len() - 1))
        .collect::<Vec<_>>();
    let known_ins_sites = variants_to_force_call
        .iter()
        .filter(|x| x.get_type() == VariantType::Insertion)
        .map(|x| (x.position(), x.get_allele1().len() - 1))
        .collect::<Vec<_>>();
    let mut variant_calls = Vec::new();
    let mut indels = Vec::new();
    for (pos, bases) in unfiltered_sites.iter() {
        let ref_base = ref_seq[*pos as usize];
        for base in bases {
            if base.len() == 1 {
                if base[0] != ref_base {
                    debug!(
                        "filtered sites at pos {}, ref_base {}, alt base {:?}",
                        *pos + 1,
                        std::str::from_utf8(&[ref_base])?,
                        std::str::from_utf8(base)?
                    );
                    let variant =
                        Variant::new_snv(0, *pos, vec![ref_base], base.to_vec(), 0, 1).unwrap();
                    if region_coordinates.variants_to_exclude.contains(&variant) {
                        debug!(
                            "filtered site at pos {}, ref_base {}, alt base {:?}",
                            *pos + 1,
                            std::str::from_utf8(&[ref_base])?,
                            std::str::from_utf8(base)?
                        );
                        continue;
                    }
                    variant_calls.push(variant);
                }
            } else {
                trace!(
                    "indel sites at pos {}, ref_base {}, alt base {:?}",
                    *pos + 1,
                    std::str::from_utf8(&[ref_base])?,
                    std::str::from_utf8(base)?
                );
                if base.contains(&b'+') {
                    let mut insertion_allele1 = base[0..1].to_vec();
                    let mut inserted_base = base[2..(base.len())].to_vec();
                    let var_len = inserted_base.len();
                    if var_len >= 9 {
                        let to_include = match_existing_indels(&known_ins_sites, pos, var_len);
                        if to_include {
                            insertion_allele1.append(&mut inserted_base);
                            let variant = Variant::new_insertion(
                                0,
                                *pos,
                                vec![ref_base],
                                insertion_allele1.to_vec(),
                                0,
                                1,
                            )
                            .unwrap();
                            indels.push(variant);
                        }
                    }
                } else if base.contains(&b'-') {
                    let mut deletion_allele0 = base[0..1].to_vec();
                    let mut deleted_base = base[2..(base.len())].to_vec();
                    let var_len = deleted_base.len();
                    if var_len >= 9 {
                        let to_include = match_existing_indels(&known_del_sites, pos, var_len);
                        if to_include {
                            deletion_allele0.append(&mut deleted_base);
                            let variant = Variant::new_deletion(
                                0,
                                *pos,
                                base.len() - 1,
                                deletion_allele0.to_vec(),
                                vec![ref_base],
                                0,
                                1,
                            )
                            .unwrap();
                            indels.push(variant);
                        }
                    }
                }
            }
        }
    }
    debug!("indels {:?}", indels);
    let mut filtered_indels: Vec<Variant> = Vec::new();
    for var in indels {
        if var.get_type() == VariantType::Deletion {
            let existing_del_sites = filtered_indels
                .iter()
                .filter(|x| x.get_type() == VariantType::Deletion)
                .map(|x| (x.position(), x.get_ref_len() - 1))
                .collect::<Vec<_>>();
            let to_include =
                match_existing_indels(&existing_del_sites, &var.position(), var.get_ref_len() - 1);
            if to_include {
                filtered_indels.push(var.clone());
            }
        }
        if var.get_type() == VariantType::Insertion {
            let existing_ins_sites = filtered_indels
                .iter()
                .filter(|x| x.get_type() == VariantType::Insertion)
                .map(|x| (x.position(), x.get_allele1().len() - 1))
                .collect::<Vec<_>>();
            let to_include = match_existing_indels(
                &existing_ins_sites,
                &var.position(),
                var.get_allele1().len() - 1,
            );
            if to_include {
                filtered_indels.push(var.clone());
            }
        }
    }
    for var in &filtered_indels {
        debug!(
            "filtered_indel {}_{}_{}",
            var.position(),
            std::str::from_utf8(&var.get_allele0())?,
            std::str::from_utf8(&var.get_allele1())?,
        );
    }
    //debug!("filtered_indels {filtered_indels:?}");
    variant_calls.append(&mut variants_to_force_call);
    variant_calls.append(&mut filtered_indels);

    //debug!("filtered sites {filtered_sites:?}");
    variant_calls.sort_by(|a, b| a.position().cmp(&b.position()));
    Ok(variant_calls)
}

/// Check if an indel matches existing indels
/// # Arguments
/// * `known_sites` - known indel sites
/// * `pos` - position
/// * `var_len` - variant length
/// # Returns
/// * `bool` - true if a match is found
fn match_existing_indels(known_sites: &Vec<(i64, usize)>, pos: &i64, var_len: usize) -> bool {
    let mut found_match = Vec::new();
    let position_buffer = if var_len < 20 { var_len as i64 } else { 5 };
    let length_buffer = if var_len < 500 {
        2.max((var_len as f64 * 0.05).round() as usize)
    } else {
        (var_len as f64 * 0.02).round() as usize
    };
    let mut n1 = *pos - position_buffer;
    for _j in 0..(position_buffer * 2 + 1) {
        let mut n2 = var_len - length_buffer;
        for _k in 0..(length_buffer * 2 + 1) {
            //debug!("known_indel_sites {known_sites:?}, {n1} {n2}");
            if known_sites.contains(&(n1, n2)) {
                found_match.push(false);
            }
            n2 += 1;
        }
        n1 += 1;
    }
    if found_match.is_empty() {
        return true;
    }
    false
}

/// Given read segments -> raw fps, select good fps and map others to them
/// # Arguments
/// * `read_segment_raw_fp` - read segment -> raw fps
/// * `read_parameters` - read parameters
/// * `read_whitelist` - read whitelist
/// * `start_end_fps` - start and end fps
/// * `is_d4z4` - whether the region is d4z4
/// * `sensitive` - whether to use sensitive mode
/// * `starting_index` - starting index
/// # Returns
/// * `BTreeMap<Vec<u8>, i32>` - fingerprint - count lookup
/// * `BTreeMap<Vec<u8>, i32>` - fingerprint sequence to name
/// * `BTreeMap<i32, Vec<u8>>` - fingerprint name to sequence
/// * `BTreeMap<Vec<u8>, i32>` - fingerprints to replace
pub fn select_fps(
    read_segment_raw_fp: &BTreeMap<String, Vec<u8>>,
    read_parameters: &ReadParameters,
    read_whitelist: &Vec<String>,
    excluded_count_segments: &HashSet<String>,
    start_end_fps: &BTreeMap<Vec<u8>, i32>,
    _is_d4z4: bool,
    sensitive: bool,
    starting_index: Option<i32>,
) -> Result<
    (
        BTreeMap<Vec<u8>, i32>,
        BTreeMap<Vec<u8>, i32>,
        BTreeMap<i32, Vec<u8>>,
        BTreeMap<Vec<u8>, i32>,
    ),
    DError,
> {
    let mut fp_count: BTreeMap<Vec<u8>, i32> = BTreeMap::<Vec<u8>, i32>::new();
    let mut fp_whitelist = HashSet::new();

    for (read, read_seq) in read_segment_raw_fp {
        if !excluded_count_segments.contains(read) {
            *fp_count.entry(read_seq.clone()).or_default() += 1;
        } else {
            debug!(
                "excluded read segment {read} {:?}",
                std::str::from_utf8(read_seq)?
            );
        }
        // white list fps, always add
        let full_read_name = read
            .split(':')
            .next()
            .ok_or("next not in read")?
            .to_string();
        if read_whitelist.contains(read) || read_whitelist.contains(&full_read_name) {
            fp_whitelist.insert(read_seq.clone());
        }
    }
    trace!("fp_count:");
    for (fp, count) in &fp_count {
        trace!("{:?}, read count {count:?}", std::str::from_utf8(fp)?);
    }

    // fp analysis
    let mut fp_index = starting_index.unwrap_or(1);
    let mut partial_fps = HashSet::new();
    // complete fp name -> seq
    let mut good_name_to_seq: BTreeMap<i32, Vec<u8>> = BTreeMap::new();
    // complete fp seq -> name
    let mut good_seq_to_name: BTreeMap<Vec<u8>, i32> = BTreeMap::new();
    // first round, very loose criteria
    for (fp, count) in fp_count.clone().into_iter() {
        let fp_seq = std::str::from_utf8(&fp)?;
        let fp_seq_string = fp_seq.to_string();
        if !fp.contains(&b'x')
            && !fp.contains(&b'-')
            && (count >= read_parameters.min_fingerprint_support - 1 || fp_whitelist.contains(&fp))
        {
            debug!("{fp_seq_string:?}, read count {count:?}, index {fp_index:?}");
            good_name_to_seq
                .entry(fp_index)
                .or_insert_with(|| fp.clone());
            good_seq_to_name
                .entry(fp.clone())
                .or_insert_with(|| fp_index);
            fp_index += 1;
        } else {
            trace!("{fp_seq_string:?} is filtered, read count {count:?}");
            partial_fps.insert(fp.clone());
        }
    }
    trace!("good_name_to_seq:");
    for (name, seq) in &good_name_to_seq {
        trace!("{name:?}, {:?}", std::str::from_utf8(seq)?);
    }
    trace!("good_seq_to_name:");
    for (seq, name) in &good_seq_to_name {
        trace!("{:?}, {name:?}", std::str::from_utf8(seq)?);
    }
    // map others to good fps
    let (_to_replace, full_unknown, _full_unknown_no_match) =
        map_partial_to_full(&good_name_to_seq, &partial_fps, true)?;

    // second round, rescue complete fingerprints with missing info
    let rescued_fps = rescue_complete_fp(full_unknown)?;
    for new_fp in &rescued_fps {
        if !good_seq_to_name.contains_key(new_fp) {
            debug!(
                "rescued complete fp with - (merged from two fps) {:?}, at index {fp_index}",
                std::str::from_utf8(new_fp)?
            );
            good_name_to_seq
                .entry(fp_index)
                .or_insert_with(|| new_fp.clone());
            good_seq_to_name
                .entry(new_fp.clone())
                .or_insert_with(|| fp_index);
            let this_fp_count = rescued_fps
                .iter()
                .map(|x| x == new_fp)
                .collect::<Vec<_>>()
                .len();
            fp_count
                .entry(new_fp.to_vec())
                .or_insert(this_fp_count as i32);
            fp_index += 1;
        }
    }

    let (_to_replace, _full_unknown, full_unknown_no_match) =
        map_partial_to_full(&good_name_to_seq, &partial_fps, true)?;
    // third round, add full_unknown_no_match: selected complete fps with unknown sites, with no match to other fps
    let mut good_name_to_seq_updated = good_name_to_seq.clone();
    let mut full_unknown_no_match_considered = Vec::new();
    for fp in &full_unknown_no_match {
        let count_unknown = fp.iter().filter(|x| **x == b'-').count();
        if count_unknown == 1 || count_unknown == 2 || fp_whitelist.contains(fp) {
            debug!(
                "consider full fp with unknown sites {:?}, index {fp_index:?}",
                std::str::from_utf8(&fp)?
            );
            good_name_to_seq_updated
                .entry(fp_index)
                .or_insert_with(|| fp.clone());
            fp_index += 1;
            full_unknown_no_match_considered.push(fp.clone());
        }
    }
    // map others to good fps
    let (to_replace, _full_unknown, _full_unknown_no_match) =
        map_partial_to_full(&good_name_to_seq_updated, &partial_fps, true)?;
    let mut to_replace_reverse: BTreeMap<Vec<u8>, Vec<Vec<u8>>> = BTreeMap::new();
    for (k, v) in &to_replace {
        if good_name_to_seq_updated.contains_key(v) {
            let v_seq = good_name_to_seq_updated.get(v).ok_or("err")?;
            if k != v_seq {
                to_replace_reverse
                    .entry(v_seq.to_vec())
                    .or_default()
                    .push(k.to_vec());
            }
        }
    }
    // fourth round, more stringent criteria
    // use full count and partial count as count filter
    let mut fp_index = starting_index.unwrap_or(1);
    let mut partial_fps = HashSet::new();
    // complete fp name -> seq
    let mut good_name_to_seq = BTreeMap::new();
    // complete fp seq -> name
    let mut good_seq_to_name: BTreeMap<Vec<u8>, i32> = BTreeMap::new();
    for (fp, count) in fp_count.clone().into_iter() {
        let fp_seq = std::str::from_utf8(&fp)?;
        let fp_seq_string = fp_seq.to_string();
        let mut partial_count = 0;
        if to_replace_reverse.contains_key(&fp) {
            partial_count = to_replace_reverse.get(&fp).ok_or("err")?.len() as i32;
        }
        if !fp.contains(&b'x')
            && (!fp.contains(&b'-') || full_unknown_no_match_considered.contains(&fp))
            && (count >= read_parameters.min_fingerprint_support
                || fp_whitelist.contains(&fp)
                || (!sensitive && partial_count >= read_parameters.min_fingerprint_support)
                || (sensitive && partial_count + count >= read_parameters.min_fingerprint_support))
        {
            if !start_end_fps.contains_key(&fp) {
                debug!("{fp_seq_string:?}, read count {count:?}, partial_count {partial_count:?}, index {fp_index:?}");
                good_name_to_seq
                    .entry(fp_index)
                    .or_insert_with(|| fp.clone());
                good_seq_to_name
                    .entry(fp.clone())
                    .or_insert_with(|| fp_index);
                fp_index += 1;
            } else {
                let start_end_fp_index = start_end_fps.get(&fp).unwrap();
                debug!("{fp_seq_string:?}, read count {count:?}, partial_count {partial_count:?}, index {start_end_fp_index:?}");
                good_name_to_seq
                    .entry(*start_end_fp_index)
                    .or_insert_with(|| fp.clone());
                good_seq_to_name
                    .entry(fp.clone())
                    .or_insert_with(|| *start_end_fp_index);
            }
        } else {
            trace!("{fp_seq_string:?} is filtered, read count {count:?}");
            partial_fps.insert(fp.clone());
        }
    }
    trace!("good_name_to_seq:");
    for (name, seq) in &good_name_to_seq {
        trace!("{name:?}, {:?}", std::str::from_utf8(seq)?);
    }
    trace!("good_seq_to_name:");
    for (seq, name) in &good_seq_to_name {
        trace!("{:?}, {name:?}", std::str::from_utf8(seq)?);
    }

    let (to_replace, full_unknown, full_unknown_no_match) =
        map_partial_to_full(&good_name_to_seq, &partial_fps, false)?;
    debug!("full_unknown fingerprints are:");
    for fp in &full_unknown {
        let fp_str = std::str::from_utf8(fp)?;
        debug!("{fp_str}");
    }

    // final round, add unknown fps with good counts directly
    for (fp, count) in &fp_count {
        let fp_seq = std::str::from_utf8(fp)?;
        let fp_seq_string = fp_seq.to_string();
        if !fp.contains(&b'x')
            && fp.contains(&b'-')
            && full_unknown_no_match.contains(fp)
            && *count >= read_parameters.min_fingerprint_support + 2
        {
            debug!("finally add unknown fps {fp_seq_string:?}, read count {count:?}, index {fp_index:?}");
            good_name_to_seq
                .entry(fp_index)
                .or_insert_with(|| fp.clone());
            good_seq_to_name
                .entry(fp.clone())
                .or_insert_with(|| fp_index);
            fp_index += 1;
        }
    }

    Ok((fp_count, good_seq_to_name, good_name_to_seq, to_replace))
}

/// Rescue complete fingerprints with missing info by comparing them with each other
/// # Arguments
/// * `full_unknown` - vector of complete repeat units that do not have an assigned fingerprint yet
/// # Returns
/// * `HashSet<Vec<u8>>` - rescued fps
fn rescue_complete_fp(full_unknown: Vec<Vec<u8>>) -> Result<HashSet<Vec<u8>>, DError> {
    let mut rescued_fps = HashSet::new();
    let mut full_unknown_rescued: BTreeMap<Vec<u8>, Vec<Vec<u8>>> = BTreeMap::new();
    let full_unknown_len = full_unknown.len();
    // 1) match two fps with -, rescue if fully match and no more -
    for i in 0..full_unknown_len {
        for j in (i + 1)..full_unknown_len {
            let hap1 = &full_unknown[i];
            let hap2 = &full_unknown[j];
            let (num_diff, _diff_sites) = edit_dis(hap1, hap2, false);
            if num_diff == 0 {
                let mut new_hap = Vec::new();
                for (x, y) in hap1.iter().zip(hap2.iter()) {
                    if *x != b'-' {
                        new_hap.push(*x);
                    } else {
                        new_hap.push(*y);
                    }
                }
                if !new_hap.contains(&b'-') {
                    debug!(
                        "tentatively merge {:?} and {:?} into full fp {:?}",
                        std::str::from_utf8(hap1)?,
                        std::str::from_utf8(hap2)?,
                        std::str::from_utf8(&new_hap)?,
                    );
                    if hap1.iter().filter(|x| **x == b'-').count() <= 2
                        && hap2.iter().filter(|x| **x == b'-').count() <= 2
                    {
                        full_unknown_rescued
                            .entry(hap1.to_vec())
                            .or_default()
                            .push(new_hap.clone());
                        full_unknown_rescued
                            .entry(hap2.to_vec())
                            .or_default()
                            .push(new_hap.clone());
                    }
                }
            }
        }
    }
    for (hap, hap_candidates) in full_unknown_rescued {
        let hap_candidates_set = hap_candidates
            .clone()
            .into_iter()
            .collect::<HashSet<Vec<u8>>>();
        trace!(
            "hap {:?} has {} merge candidates",
            std::str::from_utf8(&hap)?,
            hap_candidates_set.len()
        );
        if hap_candidates_set.len() == 1 {
            let hap_candidate = hap_candidates.first().unwrap();
            rescued_fps.insert(hap_candidate.to_vec());
        }
    }
    Ok(rescued_fps)
}

/// Find edit distance between two fingerprints
/// When strict, 'x' counts as a difference
/// # Arguments
/// * `a`
/// * `b`
/// * `strict` - whether to use strict mode
/// # Returns
/// * `(i32, Vec<usize>)` - number of differences and difference sites
pub fn edit_dis(a: &[u8], b: &[u8], strict: bool) -> (i32, Vec<usize>) {
    let mut num_diff = 0;
    let mut diff_sites = Vec::new();
    for (i, (x, y)) in a.iter().zip(b.iter()).enumerate() {
        if !strict {
            if *x != b'x' && *x != b'-' && *y != b'x' && *y != b'-' && *x != *y {
                num_diff += 1;
                diff_sites.push(i + 1);
            }
        } else if *x != b'-' && *y != b'-' && *x != *y {
            num_diff += 1;
            diff_sites.push(i + 1);
        }
    }
    (num_diff, diff_sites)
}

/// Given read segments -> raw fps, find segments that correspond to starting and ends fps
/// # Arguments
/// * `is_d4z4` - whether the region is d4z4
/// * `read_segment_raw_fp` - read segment -> raw fps
/// * `flanking_reads` - flanking reads
/// # Returns
/// * `BTreeMap<Vec<u8>, i32>` - start and end fps
pub fn get_start_end_fps(
    is_d4z4: bool,
    read_segment_raw_fp: &BTreeMap<String, Vec<u8>>,
    flanking_reads: &FlankReads,
) -> Result<BTreeMap<Vec<u8>, i32>, DError> {
    let mut start_end_fps = BTreeMap::new();
    if !is_d4z4 {
        return Ok(start_end_fps);
    }
    let mut starting_fps = BTreeMap::<Vec<u8>, i32>::new();
    let mut ending_fps = BTreeMap::<Vec<u8>, i32>::new();

    for (read_segment, fp) in read_segment_raw_fp.clone().into_iter() {
        if flanking_reads.start_segment.contains(&read_segment) {
            *starting_fps.entry(fp.clone()).or_default() += 1;
        }
        if flanking_reads.end_segment.contains(&read_segment) {
            *ending_fps.entry(fp.clone()).or_default() += 1;
        }
    }
    let mut start_counter = -2;
    for (fp, count) in starting_fps.into_iter() {
        let unknown_bases = fp.iter().filter(|x| **x == b'-').count();
        if (count >= 2 && unknown_bases == 0) || (unknown_bases == 1 && count >= 5) {
            let fp_name = start_counter;
            start_end_fps.entry(fp.clone()).or_insert(fp_name);
            debug!(
                "starting copy name {} {:?}, count {}",
                fp_name.clone(),
                std::str::from_utf8(&fp)?,
                count
            );
            start_counter -= 1
        }
    }
    let mut end_counter = -11;
    for (fp, count) in ending_fps.into_iter() {
        let unknown_bases = fp.iter().filter(|x| **x == b'-').count();
        if (count >= 2 && unknown_bases == 0) || (unknown_bases == 1 && count >= 5) {
            let fp_name = end_counter;
            start_end_fps.entry(fp.clone()).or_insert(fp_name);
            debug!(
                "ending copy name {} {:?}, count {}",
                fp_name.clone(),
                std::str::from_utf8(&fp)?,
                count
            );
            end_counter -= 1
        }
    }
    Ok(start_end_fps)
}

pub fn handle_qal_units(fp_info: FingerprintInfo) -> Result<(FingerprintInfo, Vec<i32>), DError> {
    let mut new_read_edges = BTreeMap::new();
    let mut new_replace = BTreeMap::<i32, i32>::new();
    let mut qal_units = Vec::new();

    let d4z4_region_coordinates = d4z4_coordinates();
    let long_insertion_variants = d4z4_region_coordinates
        .variants_to_call
        .iter()
        .rev()
        .skip(1)
        .take(2)
        .cloned()
        .collect::<Vec<_>>();

    let long_insertion_variant_codes = fp_info
        .variants_by_position
        .iter()
        .enumerate()
        .filter_map(|(index, (_pos, variants_at_pos))| {
            for long_insertion_variant in &long_insertion_variants {
                if let Some(alt_index) = variants_at_pos
                    .iter()
                    .position(|variant| variant == long_insertion_variant)
                {
                    if alt_index <= 8 {
                        return Some((index, b'1' + alt_index as u8));
                    }
                }
            }
            None
        })
        .collect::<Vec<_>>();

    for (unit_name, unit_fp) in fp_info.good_name_to_seq.iter() {
        for (index, expected_code) in long_insertion_variant_codes.iter() {
            if unit_fp.get(*index) == Some(expected_code) {
                debug!(
                    "unit {unit_name:?} has long insertion variant at index {index}, unit fp: {:?}",
                    std::str::from_utf8(unit_fp)?
                );
                qal_units.push(*unit_name);
                let mut unit_seq_without_insertion = unit_fp.clone();
                unit_seq_without_insertion[*index] = b'0';
                for (other_unit_name, other_unit_fp) in fp_info.good_name_to_seq.iter() {
                    if other_unit_fp == &unit_seq_without_insertion {
                        debug!("Found matching unit {other_unit_name:?} for {unit_name:?}, unit fp: {:?}", std::str::from_utf8(other_unit_fp)?);
                        new_replace.entry(*unit_name).or_insert(*other_unit_name);
                        qal_units.push(*other_unit_name);
                        break;
                    }
                }
            }
        }
    }
    debug!("new_replace: {new_replace:?}");

    if new_replace.is_empty() {
        return Ok((fp_info.clone(), qal_units));
    }

    for (each_read, read_fps) in fp_info.read_edges.iter() {
        let mut new_fps = Vec::new();
        for fp in read_fps {
            if new_replace.contains_key(fp) {
                let to_replace = new_replace.get(fp).ok_or("key not found in new_replace")?;
                new_fps.push(*to_replace);
            } else {
                new_fps.push(*fp);
            }
        }
        if new_fps != read_fps.to_vec() {
            debug!("updated edges {each_read}: from {read_fps:?} to {new_fps:?}");
        }
        new_read_edges
            .entry(each_read.to_string())
            .or_insert(new_fps);
    }

    Ok((
        FingerprintInfo {
            read_edges: new_read_edges,
            grouped_reads: fp_info.grouped_reads,
            fp_count: fp_info.fp_count,
            good_name_to_seq: fp_info.good_name_to_seq,
            read_positions: fp_info.read_positions,
            read_bases: fp_info.read_bases,
            fp_to_tid: fp_info.fp_to_tid,
            variants_by_position: fp_info.variants_by_position,
        },
        qal_units,
    ))
}

pub fn handle_last_d4z4_long_insertion(
    fp_info: FingerprintInfo,
) -> Result<FingerprintInfo, DError> {
    let mut new_read_edges = BTreeMap::new();
    let mut new_replace = BTreeMap::<i32, i32>::new();

    let d4z4_region_coordinates = d4z4_coordinates();
    let long_insertion_variant = d4z4_region_coordinates
        .variants_to_call
        .last()
        .cloned()
        .ok_or("missing d4z4 forced-call variant")?;

    let long_insertion_variant_codes = fp_info
        .variants_by_position
        .iter()
        .enumerate()
        .filter_map(|(index, (_pos, variants_at_pos))| {
            variants_at_pos
                .iter()
                .position(|variant| variant == &long_insertion_variant)
                .and_then(|alt_index| {
                    if alt_index <= 8 {
                        Some((index, b'1' + alt_index as u8))
                    } else {
                        None
                    }
                })
        })
        .collect::<Vec<_>>();

    for (unit_name, unit_fp) in fp_info.good_name_to_seq.iter() {
        for (index, expected_code) in long_insertion_variant_codes.iter() {
            if unit_fp.get(*index) == Some(expected_code) {
                debug!(
                    "unit {unit_name:?} has final d4z4 long insertion at index {index}, unit fp: {:?}",
                    std::str::from_utf8(unit_fp)?
                );
                let mut unit_seq_without_insertion = unit_fp.clone();
                unit_seq_without_insertion[*index] = b'0';
                let mut unit_seq_without_insertion_replace_unknown =
                    unit_seq_without_insertion.clone();
                for x in &mut unit_seq_without_insertion_replace_unknown {
                    if *x == b'-' {
                        *x = b'0';
                    }
                }
                let mut matching_unit_fps = Vec::new();
                for (other_unit_name, other_unit_fp) in fp_info.good_name_to_seq.iter() {
                    if other_unit_fp == &unit_seq_without_insertion
                        || other_unit_fp == &unit_seq_without_insertion_replace_unknown
                    {
                        matching_unit_fps.push(*other_unit_name);
                    }
                }
                if matching_unit_fps.len() == 1 {
                    let matching_unit_fp = matching_unit_fps.first().unwrap();
                    debug!("Found matching unit {matching_unit_fp:?} for {unit_name:?}");
                    new_replace.entry(*unit_name).or_insert(*matching_unit_fp);
                }
            }
        }
    }
    debug!("new_replace for final d4z4 long insertion: {new_replace:?}");

    if new_replace.is_empty() {
        return Ok(fp_info);
    }

    for (each_read, read_fps) in fp_info.read_edges.iter() {
        let mut new_fps = Vec::new();
        for fp in read_fps {
            if new_replace.contains_key(fp) {
                let to_replace = new_replace.get(fp).ok_or("key not found in new_replace")?;
                new_fps.push(*to_replace);
            } else {
                new_fps.push(*fp);
            }
        }
        if new_fps != read_fps.to_vec() {
            debug!("updated edges {each_read}: from {read_fps:?} to {new_fps:?}");
        }
        new_read_edges
            .entry(each_read.to_string())
            .or_insert(new_fps);
    }

    Ok(FingerprintInfo {
        read_edges: new_read_edges,
        grouped_reads: fp_info.grouped_reads,
        fp_count: fp_info.fp_count,
        good_name_to_seq: fp_info.good_name_to_seq,
        read_positions: fp_info.read_positions,
        read_bases: fp_info.read_bases,
        fp_to_tid: fp_info.fp_to_tid,
        variants_by_position: fp_info.variants_by_position,
    })
}

/// Compare fingerprints and remove redundant ones
/// # Arguments
/// * `fp_info` - fingerprint information
/// * `max_read_count_to_correct` - maximum read count allowed. Do not correct if above this value.
/// * `include_all_links` - whether to include all links
/// # Returns
/// * `(FingerprintInfo, bool)` - fingerprint information and whether any changes were made
pub fn rm_redundant_finger_prints(
    fp_info: FingerprintInfo,
    max_read_count_to_correct: i32,
    include_all_links: bool,
) -> Result<(FingerprintInfo, bool), DError> {
    let mut new_good_name_to_seq = BTreeMap::new();
    let mut new_read_edges = BTreeMap::new();
    let mut new_grouped_reads = BTreeMap::new();

    let mut new_replace = BTreeMap::<i32, i32>::new();
    for (hap1_name, hap1) in fp_info.good_name_to_seq.iter() {
        let reads1 = fp_info
            .read_edges
            .clone()
            .into_values()
            .filter(|x| x.contains(hap1_name))
            .collect::<Vec<_>>();
        let (hap1_prev, hap1_next) = get_prev_next(reads1.clone(), *hap1_name, true);
        log::trace!(
            "Evaluating {hap1_name} previous nodes {hap1_prev:?}, next nodes {hap1_next:?}"
        );
        // For D4Z4, do not replace if a node has no previous or next nodes
        if include_all_links && (hap1_prev.is_empty() || hap1_next.is_empty()) {
            continue;
        }
        let mut hap1_candidate_list = Vec::new();
        // check # of difference between this hap and all other haps
        for (hap2_name, hap2) in fp_info.good_name_to_seq.iter() {
            if hap1 != hap2 {
                let (num_diff, diff_sites) = edit_dis(hap1, hap2, false);
                let count1 = fp_info.fp_count.get(hap1).ok_or("haplotype not found")?;
                let count2 = fp_info.fp_count.get(hap2).ok_or("haplotype not found")?;
                trace!("{hap1_name}: count {count1}, vs. {hap2_name}: count {count2}, {num_diff} mismatches at sites {diff_sites:?}");
                // require count difference
                if num_diff < 2 && *count1 <= max_read_count_to_correct && *count2 >= *count1 * 2 {
                    // check if hap1 can be replaced with hap2
                    let reads2 = fp_info
                        .read_edges
                        .clone()
                        .into_values()
                        .filter(|x| x.contains(hap2_name))
                        .collect::<Vec<_>>();
                    let (mut hap2_prev, mut hap2_next) =
                        get_prev_next(reads2.clone(), *hap2_name, include_all_links);
                    hap2_prev.sort();
                    hap2_next.sort();
                    let mut reads1_replaced: Vec<Vec<i32>> = Vec::new();
                    for each_read in &reads1 {
                        let mut new_read: Vec<i32> = Vec::new();
                        for fp in each_read {
                            if fp == hap1_name {
                                new_read.push(*hap2_name)
                            } else {
                                new_read.push(*fp);
                            }
                        }
                        reads1_replaced.push(new_read);
                    }
                    let mut reads2_add = reads2.clone();
                    reads2_add.append(&mut reads1_replaced);
                    let (mut hap2_add_prev, mut hap2_add_next) =
                        get_prev_next(reads2_add.clone(), *hap2_name, include_all_links);
                    hap2_add_prev.sort();
                    hap2_add_next.sort();
                    trace!("{hap2_name} previous nodes {hap2_prev:?}, now {hap2_add_prev:?}");
                    trace!("{hap2_name} next nodes {hap2_next:?}, now {hap2_add_next:?}");
                    if hap2_prev == hap2_add_prev && hap2_next == hap2_add_next {
                        hap1_candidate_list.push(*hap2_name);
                    }
                }
            }
        }
        if hap1_candidate_list.len() == 1 {
            let to_replace_hap1 = hap1_candidate_list
                .first()
                .ok_or("first not found in hap1_candidate_list")?;
            let reads1 = fp_info
                .read_edges
                .clone()
                .into_values()
                .filter(|x| x.contains(hap1_name))
                .collect::<Vec<_>>();
            let (hap1_prev, hap1_next) =
                get_prev_next(reads1.clone(), *hap1_name, include_all_links);
            if hap1_prev == hap1_next && hap1_prev.len() == 1 && hap1_prev[0] == *to_replace_hap1 {
                debug!("do not replace fingerprint {hap1_name} with {to_replace_hap1} because its prev and next nodes are both {to_replace_hap1}");
            } else {
                debug!("replace fingerprint {hap1_name} with {to_replace_hap1}");
                new_replace.entry(*hap1_name).or_insert(*to_replace_hap1);
                //if fp_info.good_name_to_seq.contains(hap1_name) {}
            }
        }
    }
    if new_replace.is_empty() {
        return Ok((fp_info.clone(), false));
    }
    for (fp_name, fp_seq) in fp_info.good_name_to_seq.iter() {
        if !new_replace.contains_key(fp_name) {
            new_good_name_to_seq
                .entry(*fp_name)
                .or_insert(fp_seq.to_vec());
        }
    }
    for (each_read, read_fps) in fp_info.read_edges.iter() {
        let mut new_fps = Vec::new();
        for fp in read_fps {
            if new_replace.contains_key(fp) {
                let to_replace = new_replace.get(fp).ok_or("key not found in new_replace")?;
                new_fps.push(*to_replace);
            } else {
                new_fps.push(*fp);
            }
        }
        if new_fps != read_fps.to_vec() {
            debug!("updated edges {each_read}: from {read_fps:?} to {new_fps:?}");
        }
        new_read_edges
            .entry(each_read.to_string())
            .or_insert(new_fps);
    }

    for (each_segment, fp) in fp_info.grouped_reads.iter() {
        if new_replace.contains_key(fp) {
            let to_replace = new_replace.get(fp).ok_or("key not found in new_replace")?;
            new_grouped_reads
                .entry(each_segment.to_string())
                .or_insert(*to_replace);
        } else {
            new_grouped_reads
                .entry(each_segment.to_string())
                .or_insert(*fp);
        }
    }

    Ok((
        FingerprintInfo {
            read_edges: infer_unknown_fingerprints(new_read_edges),
            grouped_reads: new_grouped_reads,
            fp_count: fp_info.fp_count,
            good_name_to_seq: new_good_name_to_seq,
            read_positions: fp_info.read_positions,
            read_bases: fp_info.read_bases,
            fp_to_tid: fp_info.fp_to_tid,
            variants_by_position: fp_info.variants_by_position,
        },
        true,
    ))
}

/// Get nodes before and after
/// # Arguments
/// * `read_fps` - reads represented as vectors of fingerprints
/// * `fp_to_check` - the query fingerprint
/// * `include_all_links` - whether to include all links
/// # Returns
/// * `(Vec<i32>, Vec<i32>)` - previous and next fingerprints
fn get_prev_next(
    read_fps: Vec<Vec<i32>>,
    fp_to_check: i32,
    include_all_links: bool,
) -> (Vec<i32>, Vec<i32>) {
    let mut prev_fps: BTreeMap<i32, i32> = BTreeMap::new();
    let mut next_fps: BTreeMap<i32, i32> = BTreeMap::new();
    for each_read in read_fps {
        for (i, fp) in each_read.iter().enumerate() {
            if *fp == fp_to_check {
                if i > 0 {
                    let prev_fp = each_read[i - 1];
                    if prev_fp != 0 {
                        *prev_fps.entry(prev_fp).or_default() += 1;
                    }
                }
                if i < each_read.len() - 1 {
                    let next_fp = each_read[i + 1];
                    if next_fp != 0 {
                        *next_fps.entry(next_fp).or_default() += 1;
                    }
                }
            }
        }
    }
    let mut threshold = 2;
    if include_all_links {
        threshold = 1;
    }
    let prev_fps_f = prev_fps
        .into_iter()
        .filter(|(_x, y)| *y >= threshold)
        .map(|(x, _y)| x)
        .collect::<Vec<i32>>();
    let next_fps_f = next_fps
        .into_iter()
        .filter(|(_x, y)| *y >= threshold)
        .map(|(x, _y)| x)
        .collect::<Vec<i32>>();
    (prev_fps_f, next_fps_f)
}

/// Update edges to infer unknown nodes between two known nodes on a read
/// # Arguments
/// * `read_edges` - vector of edges per read
/// # Returns
/// * `BTreeMap<String, Vec<i32>>` - updated edges
pub fn infer_unknown_fingerprints(
    read_edges: BTreeMap<String, Vec<i32>>,
) -> BTreeMap<String, Vec<i32>> {
    let mut new_read_edges = BTreeMap::new();
    let mut previous_nodes: BTreeMap<i32, Vec<i32>> = BTreeMap::new();
    let mut next_nodes: BTreeMap<i32, Vec<i32>> = BTreeMap::new();
    for (_read, edges) in read_edges.iter() {
        let edges_clone = edges.clone();
        let adjacent_pairs = edges.iter().zip(edges_clone.iter().skip(1));
        for (node1, node2) in adjacent_pairs {
            if *node1 != 0 && *node2 != 0 {
                next_nodes.entry(*node1).or_default().push(*node2);
                previous_nodes.entry(*node2).or_default().push(*node1);
            }
        }
    }
    for (read, edges) in read_edges.iter() {
        let mut to_update = BTreeMap::new();
        let edges_len = edges.len();
        if edges_len > 2 {
            for i in 1..(edges_len - 1) {
                let prev_node = edges[i - 1];
                let this_node = edges[i];
                let next_node = edges[i + 1];
                if prev_node != 0
                    && prev_node != -1
                    && this_node == 0
                    && next_node != 0
                    && next_node != -10
                    && previous_nodes.contains_key(&next_node)
                    && next_nodes.contains_key(&prev_node)
                {
                    let prev_node_next = next_nodes.get(&prev_node).unwrap();
                    let next_node_prev = previous_nodes.get(&next_node).unwrap();
                    if prev_node_next.len() > 1 && next_node_prev.len() > 1 {
                        let prev_node_next_set: HashSet<i32> =
                            prev_node_next.iter().cloned().collect();
                        let next_node_prev_set: HashSet<i32> =
                            next_node_prev.iter().cloned().collect();
                        if prev_node_next_set.len() == 1
                            && next_node_prev_set.len() == 1
                            && prev_node_next_set == next_node_prev_set
                        {
                            let prev_node_next_set_vec: Vec<i32> =
                                prev_node_next_set.iter().cloned().collect_vec();
                            let prev_node_next_set_vec_node =
                                prev_node_next_set_vec.first().unwrap();
                            to_update.entry(i).or_insert(*prev_node_next_set_vec_node);
                        }
                    }
                }
            }
        }
        if to_update.is_empty() {
            new_read_edges
                .entry(read.to_string())
                .or_insert(edges.to_vec());
        } else {
            let mut new_edge = Vec::new();
            for j in 0..edges_len {
                let j_node = if to_update.contains_key(&j) {
                    *to_update.get(&j).unwrap()
                } else {
                    edges[j]
                };
                new_edge.push(j_node);
            }
            debug!("infer unknown nodes and update edge {edges:?} to new_edge {new_edge:?}");
            new_read_edges
                .entry(read.to_string())
                .or_insert(new_edge.to_vec());
        }
    }
    new_read_edges
}

/// Update fps with a new set of fps
/// # Arguments
/// * `fp_count` - fingerprint - count lookup
/// * `good_seq_to_name` - fingerprint sequence to name
/// * `good_name_to_seq` - fingerprint name to sequence
/// * `fps_to_add` - fingerprints to add
/// # Returns
/// * `BTreeMap<Vec<u8>, i32>` - updated fps
pub fn update_fps(
    fp_count: &BTreeMap<Vec<u8>, i32>,
    good_seq_to_name: &mut BTreeMap<Vec<u8>, i32>,
    good_name_to_seq: &mut BTreeMap<i32, Vec<u8>>,
    fps_to_add: HashSet<Vec<u8>>,
) -> Result<BTreeMap<Vec<u8>, i32>, DError> {
    //BTreeMap<Vec<u8>, i32>,
    //BTreeMap<Vec<u8>, i32>,
    //BTreeMap<i32, Vec<u8>>,
    let mut partial_fps = HashSet::new();
    let mut fp_index = good_seq_to_name.values().max().unwrap() + 1;
    for (fp, _count) in fp_count {
        if !good_seq_to_name.contains_key(fp) {
            if fps_to_add.contains(fp) {
                debug!(
                    "rescued orphan fp {:?}, at index {fp_index}",
                    std::str::from_utf8(fp)?
                );
                good_name_to_seq
                    .entry(fp_index)
                    .or_insert_with(|| fp.clone());
                good_seq_to_name
                    .entry(fp.clone())
                    .or_insert_with(|| fp_index);
                fp_index += 1;
            } else {
                partial_fps.insert(fp.clone());
            }
        }
    }

    let (to_replace, _full_unknown, _full_unknown_no_match) =
        map_partial_to_full(&good_name_to_seq, &partial_fps, false)?;

    Ok(to_replace)
}

/// Map partial fps to known fingerprints
/// # Arguments
/// * `good_name_to_seq` - fingerprint name -> sequence
/// * `partial_fps` - a set of fingerprints (partial, not complete)
/// * `strict` - whether to use strict mode
/// # Returns
/// * `(BTreeMap<Vec<u8>, i32>, Vec<Vec<u8>>, Vec<Vec<u8>>)` - to replace, full unknown, full unknown no match
fn map_partial_to_full(
    good_name_to_seq: &BTreeMap<i32, Vec<u8>>,
    partial_fps: &HashSet<Vec<u8>>,
    strict: bool,
) -> Result<(BTreeMap<Vec<u8>, i32>, Vec<Vec<u8>>, Vec<Vec<u8>>), DError> {
    let mut full_unknown = Vec::new();
    let mut full_unknown_no_match = Vec::new();
    let mut to_replace: BTreeMap<Vec<u8>, i32> = BTreeMap::new();
    for partial_fp in partial_fps {
        let mut candidates_no_mismatch = Vec::new();
        let mut candidates_one_mismatch = Vec::new();
        for (fp_index, fp) in good_name_to_seq.iter() {
            let (num_diff, _diff_sites) = edit_dis(partial_fp, fp, false);
            if num_diff < 1 {
                candidates_no_mismatch.push(fp_index);
            }
            if num_diff < 2 {
                candidates_one_mismatch.push(fp_index);
            }
        }
        let partial_fp_string = std::str::from_utf8(partial_fp)?;
        trace!("map partial fp {partial_fp_string:?} to candidates_no_mismatch {candidates_no_mismatch:?}");
        trace!(
            "map partial fp {partial_fp_string:?} to candidates_one_mismatch {candidates_one_mismatch:?}"
        );
        // if segment is too short, only trust candidates_no_mismatch
        let no_info_sites = partial_fp.iter().filter(|x| **x == b'x').count();
        if candidates_no_mismatch.len() == 1 {
            let candidates_no_mismatch_first =
                candidates_no_mismatch.first().ok_or("first not found")?;
            to_replace
                .entry(partial_fp.to_vec())
                .or_insert(**candidates_no_mismatch_first);
        } else if candidates_one_mismatch.len() == 1
            && !strict
            && no_info_sites < partial_fp.len() / 2
        {
            let candidates_one_mismatch_first =
                candidates_one_mismatch.first().ok_or("first not found")?;
            to_replace
                .entry(partial_fp.to_vec())
                .or_insert(**candidates_one_mismatch_first);
        } else if !partial_fp.contains(&b'x') {
            if candidates_no_mismatch.is_empty() {
                full_unknown_no_match.push(partial_fp.to_vec());
            }
            full_unknown.push(partial_fp.to_vec());
        }
    }
    Ok((to_replace, full_unknown, full_unknown_no_match))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_edit_dis() {
        let v1: Vec<u8> = vec![120, 65, 65, 67, 71, 45];
        let v2: Vec<u8> = vec![120, 67, 65, 120, 71, 71];
        let (n, sites) = edit_dis(&v1, &v2, false);
        assert_eq!(n, 1);
        assert_eq!(sites, [2]);
        let (n2, sites2) = edit_dis(&v1, &v2, true);
        assert_eq!(n2, 2);
        assert_eq!(sites2, [2, 4]);
    }

    #[test]
    fn test_edit_dis_identical() {
        let v: Vec<u8> = vec![b'A', b'C', b'G', b'T'];
        let (n, sites) = edit_dis(&v, &v, false);
        assert_eq!(n, 0);
        assert!(sites.is_empty());
        let (n2, sites2) = edit_dis(&v, &v, true);
        assert_eq!(n2, 0);
        assert!(sites2.is_empty());
    }

    #[test]
    fn test_edit_dis_empty() {
        let v: Vec<u8> = vec![];
        let (n, sites) = edit_dis(&v, &v, false);
        assert_eq!(n, 0);
        assert!(sites.is_empty());
    }

    #[test]
    fn test_edit_dis_non_strict_ignores_x_and_dash() {
        // Non-strict: 'x' (120) and '-' (45) are ignored - don't count as differences
        let v1: Vec<u8> = vec![b'A', b'x', b'-', b'C'];
        let v2: Vec<u8> = vec![b'G', b'x', b'-', b'T'];
        let (n, sites) = edit_dis(&v1, &v2, false);
        // Position 0: A vs G - both valid, count. Position 1: x vs x - ignore. Position 2: - vs - ignore. Position 3: C vs T - both valid, count.
        assert_eq!(n, 2);
        assert_eq!(sites, [1, 4]);
    }

    #[test]
    fn test_edit_dis_non_strict_x_vs_different() {
        // x vs non-x in non-strict: x is ignored, so no count
        let v1: Vec<u8> = vec![b'x', b'A'];
        let v2: Vec<u8> = vec![b'C', b'A'];
        let (n, sites) = edit_dis(&v1, &v2, false);
        assert_eq!(n, 0);
        assert!(sites.is_empty());
    }

    #[test]
    fn test_edit_dis_strict_counts_x() {
        // Strict: 'x' is NOT ignored - x vs different counts
        let v1: Vec<u8> = vec![b'x', b'A'];
        let v2: Vec<u8> = vec![b'C', b'A'];
        let (n, sites) = edit_dis(&v1, &v2, true);
        assert_eq!(n, 1);
        assert_eq!(sites, [1]);
    }

    #[test]
    fn test_edit_dis_strict_ignores_dash() {
        // Strict: '-' is still ignored
        let v1: Vec<u8> = vec![b'-', b'A'];
        let v2: Vec<u8> = vec![b'C', b'A'];
        let (n, sites) = edit_dis(&v1, &v2, true);
        assert_eq!(n, 0);
        assert!(sites.is_empty());

        let (n, sites) = edit_dis(&v1, &v2, false);
        assert_eq!(n, 0);
        assert!(sites.is_empty());
    }

    #[test]
    fn test_edit_dis_one_based_sites() {
        // Verify diff_sites uses 1-based indexing
        let v1: Vec<u8> = vec![b'A', b'B', b'C'];
        let v2: Vec<u8> = vec![b'X', b'B', b'Z'];
        let (n, sites) = edit_dis(&v1, &v2, true);
        assert_eq!(n, 2);
        assert_eq!(sites, [1, 3]);
    }

    #[test]
    fn test_edit_dis_different_lengths() {
        // zip stops at shorter length - only first 2 positions compared
        let v1: Vec<u8> = vec![b'A', b'B'];
        let v2: Vec<u8> = vec![b'A', b'B', b'C', b'D'];
        let (n, sites) = edit_dis(&v1, &v2, true);
        assert_eq!(n, 0);
        assert!(sites.is_empty());

        let v1: Vec<u8> = vec![b'A', b'X', b'C'];
        let v2: Vec<u8> = vec![b'A', b'Y'];
        let (n, sites) = edit_dis(&v1, &v2, true);
        assert_eq!(n, 1);
        assert_eq!(sites, [2]);
    }

    #[test]
    fn test_edit_dis_all_differ() {
        let v1: Vec<u8> = vec![b'A', b'C', b'G'];
        let v2: Vec<u8> = vec![b'T', b'G', b'C'];
        let (n, sites) = edit_dis(&v1, &v2, true);
        assert_eq!(n, 3);
        assert_eq!(sites, [1, 2, 3]);
    }

    #[test]
    fn test_rescue_complete_fp() {
        let hap1 = vec![1, 2, b'-'];
        let hap2 = vec![b'-', 2, 3];
        let haps = vec![hap1, hap2];
        let rescued_fps = rescue_complete_fp(haps)
            .unwrap()
            .into_iter()
            .collect::<Vec<Vec<u8>>>();
        assert_eq!(rescued_fps, vec![vec![1, 2, 3]]);

        let hap1 = vec![1, 2, b'-', 4];
        let hap2 = vec![b'-', 2, 3, 4];
        let haps = vec![hap1, hap2];
        let rescued_fps = rescue_complete_fp(haps)
            .unwrap()
            .into_iter()
            .collect::<Vec<Vec<u8>>>();
        assert_eq!(rescued_fps, vec![vec![1, 2, 3, 4]]);

        let hap1 = vec![1, b'-', b'-', 4, b'-'];
        let hap2 = vec![b'-', 2, 3, 4, 5];
        let haps = vec![hap1, hap2];
        let rescued_fps = rescue_complete_fp(haps)
            .unwrap()
            .into_iter()
            .collect::<Vec<Vec<u8>>>();
        assert!(rescued_fps.is_empty());

        let hap1 = vec![1, 2, b'-', 4];
        let hap2 = vec![b'-', 2, 3, 5];
        let haps = vec![hap1, hap2];
        let rescued_fps = rescue_complete_fp(haps)
            .unwrap()
            .into_iter()
            .collect::<Vec<Vec<u8>>>();
        assert!(rescued_fps.is_empty());
    }

    #[test]
    fn test_get_prev_next() {
        let reads = vec![vec![1, 2, 3], vec![2, 3], vec![1, 2], vec![2, 4]];
        let (prev, next) = get_prev_next(reads.clone(), 2, false);
        assert_eq!(prev, vec![1]);
        assert_eq!(next, vec![3]);

        let (prev, next) = get_prev_next(reads.clone(), 2, true);
        assert_eq!(prev, vec![1]);
        assert_eq!(next, vec![3, 4]);

        let (prev, next) = get_prev_next(reads.clone(), 3, false);
        assert_eq!(prev, vec![2]);
        assert!(next.is_empty());

        let (prev, next) = get_prev_next(reads.clone(), 4, false);
        assert!(prev.is_empty());
        assert!(next.is_empty());

        let (prev, next) = get_prev_next(reads.clone(), 4, true);
        assert_eq!(prev, vec![2]);
        assert!(next.is_empty());
    }

    #[test]
    fn test_infer_unknown_fingerprints() {
        let mut read_edges = BTreeMap::new();
        read_edges
            .entry(String::from("read1"))
            .or_insert(vec![1, 0, 3]);
        read_edges
            .entry(String::from("read2"))
            .or_insert(vec![1, 2]);
        read_edges
            .entry(String::from("read3"))
            .or_insert(vec![2, 3]);
        let new_reads = infer_unknown_fingerprints(read_edges);
        assert_eq!(*new_reads.get("read1").unwrap(), vec![1, 0, 3]);
        // at least two reads supporting the edge
        let mut read_edges = BTreeMap::new();
        read_edges
            .entry(String::from("read1"))
            .or_insert(vec![1, 0, 3]);
        read_edges
            .entry(String::from("read2"))
            .or_insert(vec![1, 2]);
        read_edges
            .entry(String::from("read3"))
            .or_insert(vec![2, 3]);
        read_edges
            .entry(String::from("read4"))
            .or_insert(vec![1, 2]);
        read_edges
            .entry(String::from("read5"))
            .or_insert(vec![2, 3]);
        let new_reads = infer_unknown_fingerprints(read_edges);
        assert_eq!(*new_reads.get("read1").unwrap(), vec![1, 2, 3]);
        // 2 is not unique
        let mut read_edges = BTreeMap::new();
        read_edges
            .entry(String::from("read1"))
            .or_insert(vec![1, 0, 3]);
        read_edges
            .entry(String::from("read2"))
            .or_insert(vec![1, 2]);
        read_edges
            .entry(String::from("read3"))
            .or_insert(vec![2, 3]);
        read_edges
            .entry(String::from("read4"))
            .or_insert(vec![1, 2]);
        read_edges
            .entry(String::from("read5"))
            .or_insert(vec![2, 3]);
        read_edges
            .entry(String::from("read6"))
            .or_insert(vec![4, 3]);
        let new_reads = infer_unknown_fingerprints(read_edges);
        assert_eq!(*new_reads.get("read1").unwrap(), vec![1, 0, 3]);
        // 2 is not unique
        let mut read_edges = BTreeMap::new();
        read_edges
            .entry(String::from("read1"))
            .or_insert(vec![1, 0, 3]);
        read_edges
            .entry(String::from("read2"))
            .or_insert(vec![1, 2]);
        read_edges
            .entry(String::from("read3"))
            .or_insert(vec![2, 3]);
        read_edges
            .entry(String::from("read4"))
            .or_insert(vec![1, 2]);
        read_edges
            .entry(String::from("read5"))
            .or_insert(vec![2, 3]);
        read_edges
            .entry(String::from("read6"))
            .or_insert(vec![1, 5]);
        let new_reads = infer_unknown_fingerprints(read_edges);
        assert_eq!(*new_reads.get("read1").unwrap(), vec![1, 0, 3]);
    }

    #[test]
    fn test_handle_last_d4z4_long_insertion() {
        let last_variant = d4z4_coordinates().variants_to_call.last().cloned().unwrap();
        let fp_info = FingerprintInfo {
            read_edges: BTreeMap::from([("read1".to_string(), vec![7, -10])]),
            grouped_reads: BTreeMap::new(),
            fp_count: BTreeMap::new(),
            good_name_to_seq: BTreeMap::from([(7, vec![b'1']), (8, vec![b'0'])]),
            read_positions: BTreeMap::new(),
            read_bases: BTreeMap::new(),
            fp_to_tid: BTreeMap::new(),
            variants_by_position: BTreeMap::from([(last_variant.position(), vec![last_variant])]),
        };

        let updated = handle_last_d4z4_long_insertion(fp_info).unwrap();

        assert_eq!(updated.read_edges.get("read1"), Some(&vec![8, -10]));
    }
}
