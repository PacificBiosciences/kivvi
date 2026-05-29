use crate::caller::vec_to_string;
use crate::depth::median;
use crate::repeat_unit::fingerprint::FingerprintInfo;
use crate::util::DError;
use log::{debug, trace, warn};
use rust_htslib::bam::record::Aux;
use rust_htslib::bam::Record;
use std::collections::{BTreeMap, HashSet};

/// Represents methylation information extracted from a read.
///
/// # Attributes
/// * `poses` - Vector of positions where methylation is detected.
/// * `probs` - Vector of probabilities associated with the methylation calls.
#[derive(Debug)]
pub struct MethInfo {
    pub poses: Vec<usize>,
    pub probs: Vec<u8>,
}

/// Summary of methylation information for a sample
/// # Attributes
/// * `meth_per_pos_median` - median methylation per position in a sample
/// * `meth_per_fp_median` - median methylation per fingerprint in a sample
/// * `all_sites_methyl_median` - median methylation of all sites in a sample
#[derive(Clone, Debug, Default, serde::Serialize, serde::Deserialize)]
pub struct MethSummary {
    pub meth_per_pos_median: BTreeMap<usize, f32>,
    pub meth_per_fp_median: BTreeMap<i32, f32>,
    pub all_sites_methyl_median: Option<f32>,
}

/// Summary of methylation information for an assembled allele that consistes of fingerprints
/// # Attributes
/// * `median_methylation_per_unit` - median methylation per unit in an assembled allele
/// * `methylation_per_site` - methylation per site in an assembled allele
#[derive(Clone, Debug, Default, serde::Serialize, serde::Deserialize)]
pub struct MethOutput {
    /// median_methylation_per_unit: each repeat unit has one value. Size is repeat length.
    pub median_methylation_per_unit: BTreeMap<String, String>,
    /// methylation_per_site: each CpG site has one value. Size is repeat length x number of sites.
    pub methylation_per_site: BTreeMap<String, String>,
}

/// Full sample methylation information
/// # Arguments
/// * `methyl_probs` - read -> ml
/// * `cpg_sites_per_read` - read segment  -> position on ref -> index of C
/// * `fp_info` - fingerprint information
/// # Returns
/// * `MethSummary` - summary of methylation information for a sample
/// * `segment_methyl_prob` - pos -> read segment -> methylation probability
pub fn methyl_prob_by_position(
    methyl_probs: &BTreeMap<String, Vec<u8>>,
    cpg_sites_per_read: &BTreeMap<String, BTreeMap<usize, usize>>,
    fp_info: &FingerprintInfo,
) -> Result<(MethSummary, BTreeMap<usize, BTreeMap<String, u8>>), DError> {
    // read -> index of C -> methylation prob (unnormalized)
    let mut read_methyl_bases: BTreeMap<String, BTreeMap<usize, u8>> = BTreeMap::new();
    // pos -> segment -> methylation probability
    let mut segment_methyl_prob: BTreeMap<usize, BTreeMap<String, u8>> = BTreeMap::new();
    for (read, ml_value) in methyl_probs {
        for (site_index, ml_prob) in ml_value.iter().enumerate() {
            read_methyl_bases
                .entry(read.clone())
                .or_default()
                .insert(site_index, *ml_prob);
        }
    }
    //debug!("read_methyl_bases {:?}", read_methyl_bases);
    // get probablity at ref positions
    for (read_segment_name, site_to_pos) in cpg_sites_per_read.iter() {
        let read_name = read_segment_name
            .split_terminator(':')
            .collect::<Vec<_>>()
            .first()
            .unwrap()
            .to_string();
        if read_methyl_bases.contains_key(&read_name) {
            // index of C -> methylation prob (unnormalized)
            let read_methyl_mls = read_methyl_bases
                .get(&read_name)
                .ok_or("read_name not in read_methyl_bases")?;
            for (ref_pos, c_index) in site_to_pos {
                //debug!("read {read_segment_name}, ref_pos {ref_pos} c_index {c_index} c_index in read_methyl_mls {:?}", read_methyl_mls.contains_key(c_index));
                if read_methyl_mls.contains_key(c_index) {
                    let ml_this_site = read_methyl_mls
                        .get(c_index)
                        .ok_or("c_index not in read_methyl_mls")?;
                    let ml_prob = *ml_this_site;
                    segment_methyl_prob
                        .entry(*ref_pos + 1)
                        .or_default()
                        .insert(read_segment_name.clone(), ml_prob);
                    //debug!("read {read_segment_name}, ref_pos {ref_pos} c_index {c_index} ml_prob {ml_prob}");
                }
            }
        }
    }

    let mut all_sites_methyl = Vec::new();
    // median methylation per position
    let mut meth_per_pos: BTreeMap<usize, Vec<i32>> = BTreeMap::new();
    let mut meth_per_pos_median: BTreeMap<usize, f32> = BTreeMap::new();
    // median methylation per fp
    let mut meth_per_fp: BTreeMap<i32, Vec<i32>> = BTreeMap::new();
    let mut meth_per_fp_median: BTreeMap<i32, f32> = BTreeMap::new();
    for (pos, segment_to_ml) in &segment_methyl_prob {
        for (read_segment, ml_value) in segment_to_ml {
            meth_per_pos.entry(*pos).or_default().push(*ml_value as i32);
            all_sites_methyl.push(*ml_value as i32);
            if fp_info.grouped_reads.contains_key(read_segment) {
                let this_segment_fp = fp_info.grouped_reads.get(read_segment).unwrap();
                meth_per_fp
                    .entry(*this_segment_fp)
                    .or_default()
                    .push(*ml_value as i32);
            }
        }
    }
    for (pos, this_pos_meth) in &meth_per_pos {
        let this_pos_meth_median =
            (1000.0 * median(this_pos_meth).unwrap() / 255.0).round() / 1000.0;
        let this_pos_meth_len = this_pos_meth.len();
        meth_per_pos_median.insert(*pos, this_pos_meth_median);
        trace!("methyl positions {pos} nsize {this_pos_meth_len} median methyl value {this_pos_meth_median} methyl values {this_pos_meth:?} ");
    }
    for (fp, this_fp_meth) in &meth_per_fp {
        if *fp != 0 {
            let this_fp_meth_median =
                (1000.0 * median(this_fp_meth).unwrap() / 255.0).round() / 1000.0;
            let this_fp_meth_len = this_fp_meth.len();
            meth_per_fp_median.insert(*fp, this_fp_meth_median);
            trace!("fp {fp} nsize {this_fp_meth_len} median methyl value {this_fp_meth_median} methyl values {this_fp_meth:?} ");
        }
    }
    // median methylation of the sample
    let all_sites_methyl_median = median(&all_sites_methyl);
    let all_sites_methyl_median = if all_sites_methyl_median.is_none() {
        None
    } else {
        Some((1000.0 * all_sites_methyl_median.unwrap() / 255.0).round() / 1000.0)
    };

    Ok((
        MethSummary {
            meth_per_pos_median,
            meth_per_fp_median,
            all_sites_methyl_median,
        },
        segment_methyl_prob,
    ))
}

/// Get methylation values per read matching positions on ref
/// # Arguments
/// * `fp_info` - fingerprint information
/// * `segment_methyl_prob` - pos -> segment -> methylation probability
/// * `cpg_sites_per_read` - read segment  -> position on ref -> index of C
/// * `reads_match_allele_index` - allele index -> read name -> position on allele
/// * `methyl_sites` - methyl sites
/// # Returns
/// * `MethOutput` - summary of methylation information for an assembled allele
/// * `alleles_reads_methyl_value` - vector of allele -> read -> methylation value per site (for plotting)
pub fn get_methyl_info(
    fp_info: &FingerprintInfo,
    segment_methyl_prob: &BTreeMap<usize, BTreeMap<String, u8>>,
    cpg_sites_per_read: &BTreeMap<String, BTreeMap<usize, usize>>,
    reads_match_allele_index: &BTreeMap<Vec<i32>, Vec<(String, i32)>>,
    methyl_sites: &Vec<usize>,
) -> Result<(MethOutput, Vec<BTreeMap<String, Vec<usize>>>), DError> {
    let read_edges = &fp_info.read_edges;
    let read_positions = &fp_info.read_positions;
    let nsite = methyl_sites.len();
    // allele index -> read name -> full position -> ml
    let mut alleles_reads_methyl: BTreeMap<i32, BTreeMap<String, BTreeMap<usize, usize>>> =
        BTreeMap::new();
    // allele name -> methylation value, one per fp
    let mut median_methylation_per_unit: BTreeMap<String, Vec<f32>> = BTreeMap::new();
    // allele name -> methylation value, one per site
    let mut methylation_per_site: BTreeMap<String, Vec<f32>> = BTreeMap::new();
    // get probability on complete alleles
    let mut allele_index = 0;
    for allele in reads_match_allele_index.keys() {
        allele_index += 1;
        let mut allele_reads = HashSet::new();
        let mut fp_supporting_reads: BTreeMap<String, Vec<String>> = BTreeMap::new();
        for (read, read_index_on_allele) in reads_match_allele_index[allele].iter() {
            allele_reads.insert(read.clone());
            let this_read_nodes = read_edges[read].clone();
            let this_read_positions = read_positions[read].clone();
            let mut read_map_index = 0;
            for (read_node, read_position) in this_read_nodes.iter().zip(this_read_positions.iter())
            {
                read_map_index += 1;
                let fp_index_on_allele = *read_index_on_allele + read_map_index - 1;
                if *read_node >= 0 && fp_index_on_allele >= 0 {
                    let segment_name = format!("{}:{}", read, *read_position);
                    fp_supporting_reads
                        .entry(format!("{allele_index}.{fp_index_on_allele}"))
                        .or_default()
                        .push(segment_name);
                }
            }
        }
        trace!("fp_supporting_reads {:?}", fp_supporting_reads);
        let allele_is_complete = if allele[0] < 0 { true } else { false };
        let range = if allele_is_complete {
            1..(allele.len() - 1)
        } else {
            0..(allele.len() - 1)
        };
        for i in range {
            let uniq_fp_name = format!("{allele_index}.{i}");
            trace!("allele {allele:?} uniq_fp_name {uniq_fp_name}");
            if fp_supporting_reads.contains_key(&uniq_fp_name) {
                trace!("allele {allele:?} uniq_fp_name {uniq_fp_name} in fp_supporting_reads");
                let mut this_fp_methyl_probs = Vec::new();
                let support_segments = fp_supporting_reads
                    .get(&uniq_fp_name)
                    .ok_or("uniq_fp_name not in fp_supporting_reads")?;
                for k in 0..nsite {
                    let full_pos = if allele_is_complete {
                        (i - 1) * nsite + k
                    } else {
                        i * nsite + k
                    };
                    for each_read in &allele_reads {
                        alleles_reads_methyl
                            .entry(allele_index)
                            .or_default()
                            .entry(each_read.to_string())
                            .or_default()
                            .insert(full_pos, 400);
                    }
                    let this_methyl_site = methyl_sites[k];
                    for segment_name in support_segments {
                        if cpg_sites_per_read.contains_key(segment_name) {
                            let this_segment_pos = cpg_sites_per_read
                                .get(segment_name)
                                .ok_or("segment_name not in cpg_sites_per_read")?;
                            let read_name = segment_name
                                .split_terminator(':')
                                .collect::<Vec<_>>()
                                .first()
                                .unwrap()
                                .to_string();
                            if this_segment_pos.contains_key(&(this_methyl_site - 1)) {
                                alleles_reads_methyl
                                    .entry(allele_index)
                                    .or_default()
                                    .entry(read_name)
                                    .or_default()
                                    .entry(full_pos)
                                    .and_modify(|val| {
                                        *val = 500 as usize;
                                    })
                                    .or_insert(400);
                            }
                        }
                    }
                    if segment_methyl_prob.contains_key(&this_methyl_site) {
                        let mut this_fp_this_pos_methyl_probs = Vec::new();
                        let segment_to_prob = segment_methyl_prob
                            .get(&this_methyl_site)
                            .ok_or("this_methyl_site not in segment_methyl_prob")?;
                        for (segment_name, prob) in segment_to_prob {
                            if support_segments.contains(segment_name) {
                                //trace!("allele {allele:?} uniq_fp_name {uniq_fp_name} i {i} k {k} full_pos {full_pos} segment_name {} prob {}", segment_name.to_string(), *prob);
                                let read_name = segment_name
                                    .split_terminator(':')
                                    .collect::<Vec<_>>()
                                    .first()
                                    .unwrap()
                                    .to_string();
                                alleles_reads_methyl
                                    .entry(allele_index)
                                    .or_default()
                                    .entry(read_name)
                                    .or_default()
                                    .entry(full_pos)
                                    .and_modify(|val| {
                                        *val = *prob as usize;
                                    })
                                    .or_insert(400);
                                this_fp_this_pos_methyl_probs.push(*prob as i32);
                                trace!("uniq_fp_name {uniq_fp_name} segment_name {segment_name} this_methyl_site {this_methyl_site} prob {}", *prob);
                            }
                        }
                        if !this_fp_this_pos_methyl_probs.is_empty() {
                            this_fp_methyl_probs
                                .push(median(&this_fp_this_pos_methyl_probs).unwrap());
                        } else {
                            trace!(
                                "at methyl site {this_methyl_site} this_fp_this_pos_methyl_probs {this_fp_this_pos_methyl_probs:?} is empty, cannot get median"
                            );
                            this_fp_methyl_probs.push(f32::NAN);
                        }
                    }
                }
                let allele_name = vec_to_string(&vec![allele.clone()], "-");
                let mut this_fp_methyl_probs_f32 = this_fp_methyl_probs
                    .iter()
                    .map(|x| {
                        if x.is_nan() {
                            f32::NAN
                        } else {
                            (1000.0 * x / 255.0).round() / 1000.0
                        }
                    })
                    .collect::<Vec<_>>();
                methylation_per_site
                    .entry(allele_name[0].clone())
                    .or_default()
                    .append(&mut this_fp_methyl_probs_f32);
                let this_fp_methyl_probs_i32 = this_fp_methyl_probs
                    .iter()
                    .filter(|x| !x.is_nan())
                    .map(|x| x.round() as i32)
                    .collect::<Vec<_>>();
                if !this_fp_methyl_probs_i32.is_empty() {
                    let this_fp_methyl_probs_median =
                        median(&this_fp_methyl_probs_i32).unwrap() / 255.0;
                    median_methylation_per_unit
                        .entry(allele_name[0].clone())
                        .or_default()
                        .push((this_fp_methyl_probs_median * 1000.0).round() / 1000.0);
                } else {
                    median_methylation_per_unit
                        .entry(allele_name[0].clone())
                        .or_default()
                        .push(f32::NAN);
                }
            }
        }
    }
    // vec[read name -> ml]
    let mut alleles_reads_methyl_value = Vec::new();
    for (_allele, allele_methyl) in alleles_reads_methyl.clone().into_iter() {
        let mut allele_methyl_value = BTreeMap::new();
        for (read, read_methyl_value) in allele_methyl {
            let this_read_methyl = read_methyl_value.into_values().collect::<Vec<usize>>();
            allele_methyl_value.insert(read.clone(), this_read_methyl.clone());
        }
        alleles_reads_methyl_value.push(allele_methyl_value);
    }
    let mut allele_fps_methyl_value_reformat = BTreeMap::new();
    for (k, v) in methylation_per_site {
        allele_fps_methyl_value_reformat.insert(
            k,
            v.iter()
                .map(|x| x.to_string())
                .collect::<Vec<String>>()
                .join(","),
        );
    }
    let mut allele_methyl_median_reformat = BTreeMap::new();
    for (k, v) in median_methylation_per_unit {
        allele_methyl_median_reformat.insert(
            k,
            v.iter()
                .map(|x| x.to_string())
                .collect::<Vec<String>>()
                .join(","),
        );
    }
    Ok((
        MethOutput {
            median_methylation_per_unit: allele_methyl_median_reformat,
            methylation_per_site: allele_fps_methyl_value_reformat,
        },
        alleles_reads_methyl_value,
    ))
}

/// Get methylation tags from a BAM record
/// # Arguments
/// * `record` - BAM record
/// # Returns
/// * (MM tag, ML tag)
pub fn get_methyl_tags(
    record: &rust_htslib::bam::Record,
) -> Result<Option<(String, Vec<u8>)>, DError> {
    let qname = std::str::from_utf8(record.qname())?;
    let mm: String;
    let ml: Vec<u8>;
    if let Ok(Aux::String(mm_value)) = record.aux(b"Mm") {
        mm = mm_value.to_string();
    } else {
        if let Ok(Aux::String(mm_value)) = record.aux(b"MM") {
            mm = mm_value.to_string();
        } else {
            trace!("missing Mm tag for read {qname}");
            return Ok(None);
        }
    }
    //trace!("read {qname} has Mm tag {:?}", mm);

    if let Ok(Aux::ArrayU8(ml_array)) = record.aux(b"Ml") {
        ml = ml_array.iter().collect::<Vec<_>>();
    } else {
        if let Ok(Aux::ArrayU8(ml_array)) = record.aux(b"ML") {
            ml = ml_array.iter().collect::<Vec<_>>();
        } else {
            debug!("missing Ml tag for read {qname}");
            return Ok(None);
        }
    }
    //trace!("read {qname} has Ml tag {:?}", ml);

    Ok(Some((mm, ml)))
}

// Methylation extraction code below adapted from TRGT
// https://github.com/PacificBiosciences/trgt

/// Get methylation probabilities from a BAM record
/// # Arguments
/// * `record` - BAM record
/// # Returns
/// * `Option<Vec<u8>>` - methylation probabilities on this read
pub fn get_methyl_prob(record: &rust_htslib::bam::Record) -> Result<Option<Vec<u8>>, DError> {
    let qname = std::str::from_utf8(record.qname())?;
    let bases = record.seq().as_bytes();
    let meth = get_mm_tag(record).and_then(|mm_tag| {
        get_ml_tag(record)
            .and_then(|ml_tag| parse_meth_tags(mm_tag, ml_tag))
            .and_then(|tags| {
                if record.is_reverse() {
                    decode_on_minus(&bases, &tags)
                } else {
                    decode_on_plus(&bases, &tags)
                }
            })
    });
    if let Some(ref meth_values) = meth {
        trace!(
            "read {qname} methylation values {:?} CG count {}",
            meth_values,
            meth_values.len()
        );
    }
    Ok(meth)
}

/// Decode methylation probabilities on a plus strand
/// # Arguments
/// * `bases` - bases on the plus strand
/// * `meth` - methylation information
/// # Returns
/// * `Option<Vec<u8>>` - methylation probabilities on this read
pub fn decode_on_plus(bases: &[u8], meth: &MethInfo) -> Option<Vec<u8>> {
    let mut profile = Vec::new();

    let mut num_cs_skipped = 0;
    let mut cite_index = 0;
    let mut non_cpgs_called = 0;

    for (index, base) in bases.iter().enumerate() {
        if *base != b'C' || index + 1 == bases.len() {
            continue;
        }

        let dinuc = &bases[index..index + 2];

        if dinuc == b"CG" {
            profile.push(0);
        }

        // TODO: Check if the first condition needed
        if cite_index != meth.poses.len() && num_cs_skipped == meth.poses[cite_index] {
            if dinuc == b"CG" {
                *profile.last_mut().unwrap() = meth.probs[cite_index];
            } else {
                non_cpgs_called += 1;
            }

            num_cs_skipped = 0;
            cite_index += 1;
        } else {
            num_cs_skipped += 1;
        }
    }

    if non_cpgs_called > 0 {
        warn!("Warning: non_cpgs_called = {non_cpgs_called}");
    }

    Some(profile)
}

/// Decode methylation probabilities on a minus strand
/// # Arguments
/// * `bases` - bases on the minus strand
/// * `meth` - methylation information
/// # Returns
/// * `Option<Vec<u8>>` - methylation probabilities on this read
pub fn decode_on_minus(bases: &[u8], meth: &MethInfo) -> Option<Vec<u8>> {
    let mut profile = Vec::new();

    let mut num_cs_skipped = 0;
    let mut cite_index = 0;
    let mut non_cpgs_called = 0;

    for (index_rev, base) in bases.iter().rev().enumerate() {
        let index = bases.len() - 1 - index_rev;

        if *base != b'G' || index == 0 {
            continue;
        }

        let dinuc = &bases[index - 1..index + 1];

        if dinuc == b"CG" {
            profile.push(0);
        }

        if cite_index != meth.poses.len() && num_cs_skipped == meth.poses[cite_index] {
            if dinuc == b"CG" {
                *profile.last_mut().unwrap() = meth.probs[cite_index];
            } else {
                non_cpgs_called += 1;
            }

            num_cs_skipped = 0;
            cite_index += 1;
        } else {
            num_cs_skipped += 1;
        }
    }

    if non_cpgs_called > 0 {
        warn!("Warning: non_cpgs_called = {non_cpgs_called}");
    }

    Some(profile)
}

/// Parses methylation tags from a BAM record into a `MethInfo` struct.
///
/// # Arguments
/// * `mm_tag` - The MM tag from the BAM record.
/// * `ml_tag` - The ML tag from the BAM record.
///
/// # Returns
/// Returns an `Option<MethInfo>` which is `Some` if the tags could be parsed, otherwise `None`.
fn parse_meth_tags(mm_tag: Aux, ml_tag: Aux) -> Option<MethInfo> {
    let mm_tag = match mm_tag {
        Aux::String(tag) => tag,
        _ => panic!("Unexpected MM tag format: {:?}", mm_tag),
    };

    // consider other possible modifications in MM
    let mm_tag = mm_tag
        .split_terminator(';')
        .map(std::borrow::ToOwned::to_owned)
        .collect::<Vec<String>>();
    let mut counter = 0;
    let mut c_m_mod = None;
    for each_mod in mm_tag {
        if each_mod.contains("C+m") {
            c_m_mod = Some(each_mod);
            break;
        } else {
            let each_mod_count = each_mod.split(',').count();
            counter += if each_mod_count == 0 {
                0
            } else {
                each_mod_count - 1
            };
        }
    }
    if c_m_mod.is_none() {
        return None;
    }
    let c_m_mod = c_m_mod.unwrap();
    let c_m_mod = c_m_mod
        .strip_prefix("C+m?")
        .or_else(|| c_m_mod.strip_prefix("C+m"))?;
    let c_m_mod = c_m_mod.trim_matches(',');
    if c_m_mod == "" {
        return None;
    }
    let poses = c_m_mod
        .split(',')
        .map(|n| n.parse::<usize>().unwrap())
        .collect::<Vec<usize>>();

    let mut probs = match ml_tag {
        Aux::ArrayU8(tag) => tag.iter().collect::<Vec<_>>(),
        _ => panic!("Unexpected ML tag format: {:?}", ml_tag),
    };
    probs = probs[counter..(counter + poses.len())].to_vec();
    assert_eq!(poses.len(), probs.len());

    Some(MethInfo { poses, probs })
}

/// Retrieves the MM tag from a BAM record.
///
/// # Arguments
/// * `rec` - A reference to the BAM record.
///
/// # Returns
/// Returns an `Option<Aux<'_>>` which is `Some` if the MM tag is present, otherwise `None`.
fn get_mm_tag(rec: &Record) -> Option<Aux<'_>> {
    rec.aux(b"MM").or_else(|_| rec.aux(b"Mm")).ok()
}

/// Retrieves the ML tag from a BAM record.
///
/// # Arguments
/// * `rec` - A reference to the BAM record.
///
/// # Returns
/// Returns an `Option<Aux<'_>>` which is `Some` if the ML tag is present, otherwise `None`.
fn get_ml_tag(rec: &Record) -> Option<Aux<'_>> {
    rec.aux(b"ML").or_else(|_| rec.aux(b"Ml")).ok()
}
