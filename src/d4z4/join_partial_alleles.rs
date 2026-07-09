use crate::assembly::assembler::FpGraph;
use crate::assembly::assembler_utils::find_overlapping_alleles;
use crate::depth::median;
use crate::methylation::MethOutput;
use crate::repeat_unit::fingerprint::FingerprintInfo;
use crate::util::RegionCoordinates;
use crate::util::{DError, DResult};
use crate::variant::VariantReport;
use log::debug;
use std::cmp;
use std::collections::BTreeMap;

/// Summary information of a D4Z4 allele
#[derive(Clone, Debug, Default, serde::Serialize, serde::Deserialize)]
pub struct AlleleSummary {
    /// name of the allele
    pub allele_name: String,
    /// chromosome of the allele
    pub chromosome: String,
    /// distal haplotype of the allele
    pub distal_haplotype: String,
    /// assembly status of the allele
    pub allele_type: String,
    /// size of the allele
    pub allele_size: String,
    /// methylation level of the allele
    pub methylation: f32,
    /// whether the allele ends in a qAL unit
    pub ending_in_qAL: bool,
}

/// Classify a fingerprint
/// # Arguments
/// * `variants` - variants on the fingerprint
/// * `region_coordinates` - region coordinates
/// # Returns
/// * `(String, usize, usize)` - the allele type, the number of variants for qADisruptedPolyA and qB
pub(crate) fn classify_fingerprint(
    variants: &Vec<String>,
    region_coordinates: &RegionCoordinates,
) -> (String, usize, usize) {
    let variants_to_check = region_coordinates
        .variants_to_distinguish_allele_types
        .clone();
    let variants_qa_disrupted = variants_to_check.get("qADisruptedPolyA").unwrap();
    let variants_qb = variants_to_check.get("qB").unwrap();
    let count_qa_disrupted = variants
        .iter()
        .filter(|x| variants_qa_disrupted.contains(x))
        .count();
    let count_qb = variants.iter().filter(|x| variants_qb.contains(x)).count();
    if count_qa_disrupted >= 2 && count_qb == 0 {
        return (
            String::from("qADisruptedPolyA"),
            count_qa_disrupted,
            count_qb,
        );
    }
    if count_qa_disrupted == 0 && count_qb >= 1 {
        return (String::from("qB"), count_qa_disrupted, count_qb);
    }
    if count_qa_disrupted == 0 && count_qb == 0 {
        return (String::from("qAIntactPolyA"), count_qa_disrupted, count_qb);
    }
    if count_qa_disrupted == 1 && count_qb == 0 {
        return (String::from("type4"), count_qa_disrupted, count_qb);
    }
    return (String::from("unknown"), count_qa_disrupted, count_qb);
}

/// Classify an allele based on the fingerprints
/// # Arguments
/// * `this_allele_fps_classified` - classified fingerprints on the allele
/// # Returns
/// * `String` - the allele type
pub(crate) fn classify_allele(this_allele_fps_classified: &Vec<String>) -> String {
    let allele_cn = this_allele_fps_classified.len();
    if allele_cn >= 2 {
        let counter = this_allele_fps_classified
            .iter()
            .map(|x| x.clone())
            .collect::<counter::Counter<String, i64>>();
        let most_common = counter.most_common_ordered();
        let highest_count = most_common[0].1;
        if highest_count as f64 >= allele_cn as f64 * 0.8 {
            return most_common[0].0.clone();
        }
    }
    if allele_cn >= 4 {
        let counter = this_allele_fps_classified
            .iter()
            .skip(1)
            .map(|x| x.clone())
            .collect::<counter::Counter<String, i64>>();
        let most_common = counter.most_common_ordered();
        let highest_count = most_common[0].1;
        if highest_count == allele_cn as i64 - 1 {
            return most_common[0].0.clone();
        }
    }
    if allele_cn >= 10 {
        let counter = this_allele_fps_classified
            .iter()
            .map(|x| x.clone())
            .collect::<counter::Counter<String, i64>>();
        let most_common = counter.most_common_ordered();
        let highest_count = most_common[0].1;
        if highest_count as f64 >= allele_cn as f64 * 0.7 {
            return most_common[0].0.clone();
        }
    }
    return String::from("unknown");
}

pub(crate) fn is_cis_dup_by_read_start_offset(
    allele: &str,
    fp_info: &FingerprintInfo,
) -> Result<bool, DError> {
    debug!("checking cis dup by read start offset for allele {allele}");
    let allele_nodes = allele
        .split('-')
        .filter(|node| !node.contains("Flank"))
        .map(|node| node.parse::<i32>())
        .collect::<Result<Vec<_>, _>>()?;
    if allele_nodes.len() < 2 {
        return Ok(false);
    }

    let mut supporting_reads = 0;
    let mut delayed_start_reads = 0;
    for (read, read_nodes) in fp_info.read_edges.iter() {
        let Some(read_positions) = fp_info.read_positions.get(read) else {
            continue;
        };
        if read_nodes.len() != read_positions.len() {
            continue;
        }

        let start_idx = 0;
        if read_nodes[start_idx] == allele_nodes[0] {
            let overlap_len = cmp::min(read_nodes.len() - start_idx, allele_nodes.len());
            if overlap_len < 2 {
                continue;
            }

            let nodes_in_read = &read_nodes[start_idx..(start_idx + overlap_len)];
            let nodes_in_allele = &allele_nodes[..overlap_len];
            let mut match_count = 0;
            let mut has_mismatch = false;
            for (read_node, allele_node) in nodes_in_read.iter().zip(nodes_in_allele.iter()) {
                if *read_node == 0 {
                    continue;
                }
                if read_node == allele_node {
                    match_count += 1;
                } else {
                    has_mismatch = true;
                    break;
                }
            }

            if !has_mismatch && match_count > 1 {
                supporting_reads += 1;
                debug!("supporting read {read} edges {read_nodes:?} positions {read_positions:?}");
                if read_positions[start_idx] > 300 {
                    delayed_start_reads += 1;
                    debug!("delayed start read {read} edges {read_nodes:?} positions {read_positions:?}");
                }
            }
        }
    }
    debug!("supporting_reads {supporting_reads} delayed_start_reads {delayed_start_reads}");
    let delayed_start_threshold = (supporting_reads as f64 * 0.8).floor() as i32;
    Ok(supporting_reads >= 3
        && delayed_start_reads >= (supporting_reads - 1).min(delayed_start_threshold))
}

fn is_cis_dup(allele: &str, fp_info: &FingerprintInfo) -> Result<bool, DError> {
    let first_node = allele.split("-").next().unwrap_or("");
    if first_node == "LeftFlank" {
        return Ok(false);
    }
    if first_node == "RightFlank" {
        return Ok(true);
    }
    if first_node != "LeftFlank" {
        let first_node = first_node.parse::<i32>()?;
        if let Some(first_node_seq) = fp_info.good_name_to_seq.get(&first_node) {
            if first_node_seq[0] == b'S' {
                return Ok(true);
            }
        }
    }
    is_cis_dup_by_read_start_offset(allele, fp_info)
}

/// Check if the two partial alleles have a left flank for one and a right flank for the other
/// # Arguments
/// * `alleles` - two partial alleles to merge
/// # Returns
/// * `bool` - true if one allele has a left flank and the other has a right flank, false otherwise
fn check_flank_presence(alleles: &Vec<String>) -> bool {
    let mut has_left = false;
    let mut has_right = false;
    for allele in alleles.iter() {
        if allele.contains("LeftFlank") && !allele.contains("RightFlank") {
            has_left = true;
        }
        if allele.contains("RightFlank") && !allele.contains("LeftFlank") {
            has_right = true;
        }
    }
    has_left && has_right
}

/// Try to merge an unknown distal allele with a matching proximal allele
/// # Arguments
/// * `unknown_allele` - the unknown allele in the distal side
/// * `unknown_allele_distal` - the distal background of the unknown allele
/// * `allele_match` - map of allele types to their alleles
/// * `all_starts_hap_backgrounds` - background of all starts haplotypes
/// * `all_ends_hap_backgrounds` - background of all ends haplotypes
/// * `distal_alleles_handled` - mutable vector to track handled distal alleles
/// * `fp_graph` - fingerprint graph
/// * `fp_info` - fingerprint information
/// # Returns
/// * `Option<AlleleSummary>` - summary of merged allele if merge was successful
fn try_merge_unknown_distal_allele(
    unknown_allele: &String,
    unknown_allele_distal: &String,
    allele_match: &BTreeMap<String, Vec<String>>,
    pairs_of_alleles_to_merge: &mut Vec<Vec<String>>,
) -> DResult {
    // Try to merge with qAIntactPolyA
    if unknown_allele_distal == "qAIntactPolyA" {
        if let Some(qa_intact_polya_alleles) = allele_match.get("qAIntactPolyA") {
            if qa_intact_polya_alleles.len() == 1 {
                let qa_intact_polya_allele = &qa_intact_polya_alleles[0];
                if qa_intact_polya_allele.contains("LeftFlank") {
                    let merged_allele =
                        vec![qa_intact_polya_allele.clone(), unknown_allele.clone()];
                    debug!("merge qAIntactPolyA partial alleles {merged_allele:?}");
                    pairs_of_alleles_to_merge.push(merged_allele);
                }
            }
        }
    }
    // Try to merge with qB
    if unknown_allele_distal == "qB" {
        if let Some(qb_alleles) = allele_match.get("qB") {
            if qb_alleles.len() == 1 {
                let qb_allele = &qb_alleles[0];
                if qb_allele.contains("LeftFlank") {
                    let merged_allele = vec![qb_allele.clone(), unknown_allele.clone()];
                    debug!("merge qB partial alleles {merged_allele:?}");
                    pairs_of_alleles_to_merge.push(merged_allele);
                }
            }
        }
    }

    Ok(())
}

/// Try to merge an unknown proximal chr4 allele with a matching distal allele
/// # Arguments
/// * `unknown_allele` - the unknown allele in the proximal side
/// * `allele_match` - map of allele types to their alleles
/// * `assembled_chr4_allele` - number of assembled chr4 alleles
/// * `all_starts_hap_backgrounds` - background of all starts haplotypes
/// * `all_ends_hap_backgrounds` - background of all ends haplotypes
/// * `distal_alleles_handled` - mutable vector to track handled distal alleles
/// * `fp_graph` - fingerprint graph
/// * `fp_info` - fingerprint information
/// # Returns
/// * `Vec<AlleleSummary>` - summaries of merged alleles
fn try_merge_unknown_proximal_allele(
    unknown_allele: &String,
    allele_match: &BTreeMap<String, Vec<String>>,
    assembled_chr4_allele: usize,
    pairs_of_alleles_to_merge: &mut Vec<Vec<String>>,
) -> DResult {
    // Try to merge with qAIntactPolyA
    if let Some(qa_intact_polya_alleles) = allele_match.get("qAIntactPolyA") {
        if qa_intact_polya_alleles.len() == 1 {
            let qa_intact_polya_allele = &qa_intact_polya_alleles[0];
            if qa_intact_polya_allele.contains("RightFlank") {
                let should_merge = (!allele_match.contains_key("qB") && assembled_chr4_allele == 1)
                    || (allele_match.contains_key("qB")
                        && allele_match
                            .get("qB")
                            .map(|v| v.len() == 2)
                            .unwrap_or(false));

                if should_merge {
                    let merged_allele =
                        vec![unknown_allele.clone(), qa_intact_polya_allele.clone()];
                    debug!("merge qAIntactPolyA partial alleles {merged_allele:?}");
                    pairs_of_alleles_to_merge.push(merged_allele);

                    // Also merge qB if it has 2 alleles
                    if let Some(qb_alleles) = allele_match.get("qB") {
                        if qb_alleles.len() == 2 && check_flank_presence(qb_alleles) {
                            debug!("merge qB partial alleles {qb_alleles:?}");
                            if !pairs_of_alleles_to_merge.contains(&qb_alleles.clone()) {
                                pairs_of_alleles_to_merge.push(qb_alleles.clone());
                            }
                        }
                    }
                }
            }
        }
    }
    // Try to merge with qB
    else if let Some(qb_alleles) = allele_match.get("qB") {
        if qb_alleles.len() == 1 {
            let qb_allele = &qb_alleles[0];
            if qb_allele.contains("RightFlank") {
                let should_merge = (!allele_match.contains_key("qAIntactPolyA")
                    && assembled_chr4_allele == 1)
                    || (allele_match.contains_key("qAIntactPolyA")
                        && allele_match
                            .get("qAIntactPolyA")
                            .map(|v| v.len() == 2)
                            .unwrap_or(false));

                if should_merge {
                    let merged_allele = vec![unknown_allele.clone(), qb_allele.clone()];
                    debug!("merge qB partial alleles {merged_allele:?}");
                    pairs_of_alleles_to_merge.push(merged_allele);

                    // Also merge qAIntactPolyA if it has 2 alleles
                    if let Some(qa_intact_polya_alleles) = allele_match.get("qAIntactPolyA") {
                        if qa_intact_polya_alleles.len() == 2
                            && check_flank_presence(qa_intact_polya_alleles)
                        {
                            debug!(
                                "merge qAIntactPolyA partial alleles {qa_intact_polya_alleles:?}"
                            );
                            if !pairs_of_alleles_to_merge.contains(&qa_intact_polya_alleles.clone())
                            {
                                pairs_of_alleles_to_merge.push(qa_intact_polya_alleles.clone());
                            }
                        }
                    }
                }
            }
        }
    }

    Ok(())
}

/// Merge two partial alleles, consider overlap
/// # Arguments
/// * `alleles` - two partial alleles to merge
/// * `fp_graph` - fingerprint graph
/// * `fp_info` - fingerprint information
/// # Returns
/// * `(String, String)` - the merged allele name and size
fn merge_two_partial_alleles(
    alleles: &Vec<String>,
    fp_graph: &FpGraph,
    fp_info: &FingerprintInfo,
) -> (usize, String, String) {
    let cyclic_nodes = &fp_graph.cyclic_nodes;
    let grouped_reads = &fp_info.grouped_reads;

    let allele1 = alleles[0].clone();
    let allele2 = alleles[1].clone();
    let allele2_end = allele2.split("-").map(|x| x.to_string()).last().unwrap();
    let ct1 = allele1
        .split("-")
        .map(|x| x.to_string())
        .filter(|x| !x.contains("Flank"))
        .collect::<Vec<String>>();
    let mut ct1 = ct1
        .iter()
        .map(|x| x.parse::<i32>().unwrap())
        .collect::<Vec<i32>>();
    let ct2 = allele2
        .split("-")
        .map(|x| x.to_string())
        .filter(|x| !x.contains("Flank"))
        .collect::<Vec<String>>();
    let mut ct2 = ct2
        .iter()
        .map(|x| x.parse::<i32>().unwrap())
        .collect::<Vec<i32>>();
    let min_len = cmp::min(ct1.len(), ct2.len());
    // overlapping
    let mut ovl_len = 0;
    for j in 0..min_len {
        let part1 = &ct1[(ct1.len() - (min_len - j))..];
        let part2 = &ct2[..(min_len - j)];
        if part1 == part2 {
            ovl_len = min_len - j;
            break;
        }
    }
    if ovl_len == 0 {
        let allele_size = allele1.split("-").filter(|x| !x.contains("Flank")).count()
            + allele2.split("-").filter(|x| !x.contains("Flank")).count();
        return (ovl_len, alleles.join("..."), format!(">={allele_size}"));
    }
    let ovl_region = &ct2[..ovl_len];
    debug!("ovl_len {ovl_len} ovl_region {ovl_region:?}");
    let counter = ovl_region
        .iter()
        .map(|x| *x)
        .collect::<counter::Counter<i32, i64>>();
    let most_common = counter.most_common_ordered();
    // cyclic node
    if ovl_len >= 2 && most_common.len() == 1 && most_common[0].1 == ovl_len as i64 {
        let cyclic_node = most_common[0].0;
        if cyclic_nodes.contains_key(&cyclic_node) {
            let cyclic_node_depth = grouped_reads
                .iter()
                .filter(|(_, v)| **v == cyclic_node)
                .count();
            let cyclic_node_cn = cyclic_nodes
                .get(&cyclic_node)
                .unwrap()
                .iter()
                .max()
                .unwrap();
            debug!("cyclic_node {cyclic_node} cyclic_node_depth {cyclic_node_depth} cyclic_node_cn {cyclic_node_cn}");
            if cyclic_node_depth <= 40 && *cyclic_node_cn >= 3 {
                while let Some(&last_element) = ct1.last() {
                    if last_element == cyclic_node {
                        ct1.pop();
                    } else {
                        break;
                    }
                }
                while let Some(&first_element) = ct2.first() {
                    if first_element == cyclic_node {
                        ct2.remove(0);
                    } else {
                        break;
                    }
                }
                let merged_allele = [ct1, vec![cyclic_node; *cyclic_node_cn], ct2].concat();
                let allele_size = merged_allele.len();
                let mut new_allele = merged_allele
                    .iter()
                    .map(|x| x.to_string())
                    .collect::<Vec<String>>();
                new_allele.insert(0, String::from("LeftFlank"));
                new_allele.push(allele2_end);
                return (ovl_len, new_allele.join("-"), format!("{allele_size}"));
            }
        }
    }
    // non-cyclic scenario
    let new_ct = [&ct1[..], &ct2[ovl_len..]].concat();
    let allele_size = new_ct.len();
    return (ovl_len, alleles.join("..."), format!(">={allele_size}"));
}

/// Get the summary of an allele
/// # Arguments
/// * `alleles` - two partial alleles to merge
/// * `all_starts_hap_backgrounds` - background of all starts haplotypes
/// * `all_ends_hap_backgrounds` - background of all ends haplotypes
/// # Returns
/// * `AlleleSummary` - summary of the merged allele
fn get_allele_summary(
    alleles: &Vec<String>,
    all_starts_hap_backgrounds: &BTreeMap<String, String>,
    all_ends_hap_backgrounds: &BTreeMap<String, String>,
    distal_alleles_handled: &mut Vec<String>,
    proximal_alleles_handled: &mut Vec<String>,
    fp_graph: &FpGraph,
    fp_info: &FingerprintInfo,
    methyl_values: &BTreeMap<String, Vec<Vec<i32>>>,
    qal_alleles: &Vec<String>,
) -> Result<AlleleSummary, DError> {
    let mut sorted_alleles = alleles.clone();
    sorted_alleles.sort_by(|a, b| b.contains("LeftFlank").cmp(&a.contains("LeftFlank")));
    debug!("sorted_alleles {:?}", sorted_alleles);
    distal_alleles_handled.push(sorted_alleles[1].clone());
    proximal_alleles_handled.push(sorted_alleles[0].clone());
    let ending_in_qAL = qal_alleles.contains(&sorted_alleles[1]);
    let chr_info = if all_starts_hap_backgrounds.contains_key(&sorted_alleles[0]) {
        all_starts_hap_backgrounds.get(&sorted_alleles[0]).unwrap()
    } else {
        &String::from("unknown")
    };
    let polya_info = if all_ends_hap_backgrounds.contains_key(&sorted_alleles[1]) {
        all_ends_hap_backgrounds.get(&sorted_alleles[1]).unwrap()
    } else {
        &String::from("unknown")
    };
    let methylation_value = get_methylation_value(&sorted_alleles[1], methyl_values, None)?;

    let (_ovl_len, allele_name, allele_size) =
        merge_two_partial_alleles(&sorted_alleles, fp_graph, fp_info);
    Ok(AlleleSummary {
        allele_name: allele_name,
        chromosome: chr_info.clone().replace("chromosome_unknown", "unknown"),
        distal_haplotype: polya_info.clone(),
        allele_type: String::from("merged"),
        allele_size: allele_size,
        methylation: methylation_value,
        ending_in_qAL: ending_in_qAL,
    })
}

/// Collect pairs of partial alleles that should be merged
/// # Arguments
/// * `allele_match` - map of allele types to their alleles
/// * `assembled_chr4_allele` - number of assembled chr4 alleles
/// * `partial_allele_number_match` - whether the number of left and right partial alleles match
/// * `all_starts_hap_backgrounds` - background of all starts haplotypes
/// * `all_ends_hap_backgrounds` - background of all ends haplotypes
/// * `pairs_of_alleles_to_merge` - mutable vector to collect pairs of alleles to merge
/// # Returns
/// * `DResult` - result indicating success or failure
pub(crate) fn collect_partial_alleles_to_merge(
    allele_match: &BTreeMap<String, Vec<String>>,
    assembled_chr4_allele: usize,
    assembled_chr10_allele: usize,
    partial_allele_number_match: bool,
    all_starts_hap_backgrounds: &BTreeMap<String, String>,
    all_ends_hap_backgrounds: &BTreeMap<String, String>,
    pairs_of_alleles_to_merge: &mut Vec<Vec<String>>,
) -> DResult {
    // no unknown partial alleles
    if !allele_match.contains_key("unknown") {
        for (allele_type, alleles) in allele_match.iter() {
            if alleles.len() == 2 {
                if check_flank_presence(alleles) {
                    debug!("merge {allele_type} partial alleles {alleles:?}");
                    pairs_of_alleles_to_merge.push(alleles.clone());
                }
            }
        }
    } else if assembled_chr4_allele == 1 {
        // already assembled a complete chr4 allele
        for (allele_type, alleles) in allele_match.iter() {
            if alleles.len() == 2
                && (allele_type == &String::from("qAIntactPolyA")
                    || allele_type == &String::from("qB"))
            {
                if check_flank_presence(alleles) {
                    debug!("merge {allele_type} partial alleles {alleles:?}");
                    pairs_of_alleles_to_merge.push(alleles.clone());
                }
            }
        }
    } else if assembled_chr10_allele == 1 {
        // already assembled a complete chr10 allele
        for (allele_type, alleles) in allele_match.iter() {
            if alleles.len() == 2 && allele_type == &String::from("qADisruptedPolyA") {
                if check_flank_presence(alleles) {
                    debug!("merge {allele_type} partial alleles {alleles:?}");
                    pairs_of_alleles_to_merge.push(alleles.clone());
                }
            }
        }
    }
    if partial_allele_number_match {
        // the number of left and right partial alleles are the same and less than 4
        // there is an unknown partial allele that can find a match in the other end
        if let Some(unknown_alleles) = allele_match.get("unknown") {
            if unknown_alleles.len() == 1 {
                let unknown_allele = &unknown_alleles[0];

                // Try to merge unknown allele in the distal side
                if let Some(unknown_allele_distal) = all_ends_hap_backgrounds.get(unknown_allele) {
                    let _ = try_merge_unknown_distal_allele(
                        unknown_allele,
                        unknown_allele_distal,
                        allele_match,
                        pairs_of_alleles_to_merge,
                    )?;
                }

                // Try to merge unknown allele in the proximal side (must be chr4)
                if let Some(unknown_allele_proximal) =
                    all_starts_hap_backgrounds.get(unknown_allele)
                {
                    if unknown_allele_proximal.contains("chr4") {
                        let _ = try_merge_unknown_proximal_allele(
                            unknown_allele,
                            allele_match,
                            assembled_chr4_allele,
                            pairs_of_alleles_to_merge,
                        )?;
                    }
                }
            }
        }
    }
    Ok(())
}

/// Get the median methylation value of the last 505 sites of an allele
/// # Arguments
/// * `allele` - allele name
/// * `methyl_values` - methylation values of each allele
/// * `last_n_sites` - number of sites to consider for methylation value
/// # Returns
/// * `f32` - methylation value
fn get_methylation_value(
    allele: &String,
    methyl_values: &BTreeMap<String, Vec<Vec<i32>>>,
    last_n_sites: Option<usize>,
) -> Result<f32, DError> {
    if !methyl_values.contains_key(allele) {
        return Ok(f32::NAN);
    }
    let last_n_sites = 5;
    let methylation: Vec<Vec<i32>> = methyl_values
        .get(allele)
        .ok_or("allele not in methyl_values")?
        .clone();
    let start_site = methylation.len().saturating_sub(last_n_sites);
    let last_n = &methylation[start_site..].to_vec();
    let last_n_i32 = last_n.iter().flatten().cloned().collect::<Vec<i32>>();
    if last_n_i32.is_empty() {
        return Ok(f32::NAN);
    }
    let methylation_value = (1000.0 * last_n_i32.iter().filter(|x| **x >= 128).count() as f32
        / last_n_i32.len() as f32)
        .round()
        / 1000.0;
    Ok(methylation_value)
}

/// Remove redundant alleles
/// # Arguments
/// * `alleles_to_check` - alleles to check
/// * `num_turns` - number of turns to remove redundant alleles
/// # Returns
/// * `Vec<String>` - remaining alleles
fn remove_redundant_alleles(
    alleles_to_check: &Vec<String>,
    num_turns: usize,
) -> Result<Vec<String>, DError> {
    let mut alleles_to_remove = Vec::new();
    for _ in 0..num_turns {
        let alleles_to_check_haps: BTreeMap<Vec<i32>, String> = alleles_to_check
            .iter()
            .map(|x| {
                (
                    x.split("-")
                        .filter(|x| !x.contains("Flank"))
                        .map(|x| x.parse::<i32>().unwrap())
                        .collect::<Vec<i32>>(),
                    x.clone(),
                )
            })
            .collect::<BTreeMap<Vec<i32>, String>>();
        debug!("Find overlapping alleles in alleles to check: {alleles_to_check_haps:?}");
        let (_overlapping_alleles, overlapping_alleles_match) = find_overlapping_alleles(
            alleles_to_check_haps
                .keys()
                .cloned()
                .collect::<Vec<Vec<i32>>>(),
            None,
        )?;
        let mut removed_one_redundant = false;
        for (allele, allele_match_info) in overlapping_alleles_match.iter() {
            let allele_name = alleles_to_check_haps.get(allele).unwrap();
            let allele_size = allele.len();
            for (matching_allele, overlap_len) in allele_match_info.iter() {
                let matching_allele_name = alleles_to_check_haps.get(matching_allele).unwrap();
                let matching_allele_size = matching_allele.len();
                if allele_size < matching_allele_size
                    && *overlap_len >= 5
                    && *overlap_len > allele_size / 2
                {
                    alleles_to_remove.push(allele_name.clone());
                    debug!("remove redundant allele: {allele_name}, redundant with {matching_allele_name}");
                    removed_one_redundant = true;
                    break;
                }
            }
            if removed_one_redundant {
                break;
            }
        }
    }
    Ok(alleles_to_remove)
}
/// Join partial alleles
/// # Arguments
/// * `all_starts_hap_backgrounds` - background of all starts haplotypes
/// * `all_ends_hap_backgrounds` - background of all ends haplotypes
/// * `complete_hap_backgrounds` - background of complete haplotypes
/// * `variant_report` - variant report
/// * `region_coordinates` - region coordinates
/// # Returns
/// * `Vec<AlleleSummary>` - summary of all alleles
pub fn join_partial_alleles(
    all_starts_hap_backgrounds: &mut BTreeMap<String, String>,
    all_ends_hap_backgrounds: &mut BTreeMap<String, String>,
    complete_hap_backgrounds: &mut BTreeMap<String, String>,
    variant_report: &VariantReport,
    region_coordinates: &RegionCoordinates,
    fp_graph: &FpGraph,
    fp_info: &FingerprintInfo,
    all_ends_ml_per_allele: &Option<BTreeMap<String, Vec<Vec<i32>>>>,
    qal_alleles: &Vec<String>,
) -> Result<Vec<AlleleSummary>, DError> {
    let mut methyl_values: BTreeMap<String, Vec<Vec<i32>>> = BTreeMap::new();
    if all_ends_ml_per_allele.is_some() {
        methyl_values = all_ends_ml_per_allele.clone().unwrap().clone();
    }
    let mut merged_allele_summary = Vec::new();
    let mut distal_alleles_handled = Vec::new();
    let mut proximal_alleles_handled = Vec::new();
    let complete_allele_variants = variant_report.complete_allele_variants.clone();
    let mut variants = variant_report.fp_variants_on_incomplete_alleles.clone();
    for (_allele, allele_variants) in complete_allele_variants.iter() {
        for (fp_id, variant) in allele_variants.iter() {
            variants.insert(fp_id.1, variant.clone());
        }
    }

    // first remove some redundant alleles
    if all_starts_hap_backgrounds.len() == 5 && all_ends_hap_backgrounds.len() >= 4 {
        let proximal_to_remove = remove_redundant_alleles(
            &all_starts_hap_backgrounds
                .keys()
                .cloned()
                .collect::<Vec<String>>(),
            1,
        )?;
        if proximal_to_remove.len() == 1
            && !(all_ends_hap_backgrounds.len() == 4
                && all_ends_hap_backgrounds.contains_key(&proximal_to_remove[0]))
        {
            let all_starts_hap_backgrounds = all_starts_hap_backgrounds
                .remove(&proximal_to_remove[0])
                .unwrap();
            if complete_hap_backgrounds.contains_key(&proximal_to_remove[0]) {
                let _ = complete_hap_backgrounds.remove(&proximal_to_remove[0]);
            }
        }
    }
    if all_starts_hap_backgrounds.len() == 4 && all_ends_hap_backgrounds.len() == 5 {
        let distal_no_cis_dup = all_ends_hap_backgrounds
            .keys()
            .filter(|x| !is_cis_dup(x, fp_info).unwrap_or(false))
            .cloned()
            .collect::<Vec<String>>();
        if distal_no_cis_dup.len() == 5 {
            let distal_to_remove = remove_redundant_alleles(&distal_no_cis_dup, 1)?;
            if distal_to_remove.len() == 1
                && !(all_starts_hap_backgrounds.len() == 4
                    && all_starts_hap_backgrounds.contains_key(&distal_to_remove[0]))
            {
                let _ = all_ends_hap_backgrounds
                    .remove(&distal_to_remove[0])
                    .unwrap();
                if complete_hap_backgrounds.contains_key(&distal_to_remove[0]) {
                    let _ = complete_hap_backgrounds.remove(&distal_to_remove[0]);
                }
            }
        }
    }

    let mut assembled_chr4_allele = 0;
    let mut assembled_chr10_allele = 0;
    for (allele, background) in complete_hap_backgrounds.iter() {
        if background.contains("chr4") {
            assembled_chr4_allele += 1;
        }
        if background.contains("chr10") {
            assembled_chr10_allele += 1;
        }
        let allele_size = allele.split("-").count() - 2;
        let background_parts = background.split("-").collect::<Vec<&str>>();
        let distal_haplotype = background_parts[0].to_string();
        let chromosome = background_parts[1].to_string();
        let methylation_value = get_methylation_value(allele, &methyl_values, None)?;

        merged_allele_summary.push(AlleleSummary {
            allele_name: allele.clone(),
            chromosome: chromosome.replace("chromosome_unknown", "unknown"),
            distal_haplotype: distal_haplotype,
            allele_type: String::from("assembled"),
            allele_size: allele_size.to_string(),
            methylation: methylation_value,
            ending_in_qAL: qal_alleles.contains(allele),
        });
        distal_alleles_handled.push(allele.clone());
        proximal_alleles_handled.push(allele.clone());
    }

    let mut allele_match: BTreeMap<String, Vec<String>> = BTreeMap::new();
    let mut partial_allele_starts = 0;
    let mut partial_allele_ends = 0;
    for (allele, background) in all_starts_hap_backgrounds.iter() {
        if !complete_hap_backgrounds.contains_key(allele) {
            partial_allele_starts += 1;
            let mut this_allele_fps_classified = Vec::new();
            let nodes = allele.split("-").collect::<Vec<&str>>();
            for node in nodes {
                if !node.contains("Flank") {
                    let node_i32 = node.parse::<i32>()?;
                    if variants.contains_key(&node_i32) {
                        let this_node_variant = variants.get(&node_i32).unwrap();
                        let (allele_type, count_qa_disrupted, count_qb) =
                            classify_fingerprint(this_node_variant, region_coordinates);
                        debug!("{allele} node {node_i32} allele_type {allele_type} count_qb {count_qb} count_qa_disrupted {count_qa_disrupted}");
                        this_allele_fps_classified.push(allele_type);
                    }
                }
            }
            let allele_type = classify_allele(&this_allele_fps_classified);
            debug!("{allele} {background} allele_type {allele_type}");
            allele_match
                .entry(allele_type.clone())
                .or_default()
                .push(allele.clone());
        }
    }
    for (allele, background) in all_ends_hap_backgrounds.iter() {
        let is_cis_dup = is_cis_dup(allele, fp_info)?;
        if !complete_hap_backgrounds.contains_key(allele) && !is_cis_dup {
            partial_allele_ends += 1;
            let mut this_allele_fps_classified = Vec::new();
            let nodes = allele.split("-").collect::<Vec<&str>>();
            for node in nodes {
                if !node.contains("Flank") {
                    let node_i32 = node.parse::<i32>()?;
                    if variants.contains_key(&node_i32) {
                        let this_node_variant = variants.get(&node_i32).unwrap();
                        let (allele_type, count_qa_disrupted, count_qb) =
                            classify_fingerprint(this_node_variant, region_coordinates);
                        debug!("{allele} node {node_i32} allele_type {allele_type} count_qb {count_qb} count_qa_disrupted {count_qa_disrupted}");
                        this_allele_fps_classified.push(allele_type);
                    }
                }
            }
            let allele_type = classify_allele(&this_allele_fps_classified);
            debug!("{allele} {background} allele_type {allele_type}");
            allele_match
                .entry(allele_type.clone())
                .or_default()
                .push(allele.clone());
        }
    }

    let partial_allele_number_match =
        (partial_allele_starts == partial_allele_ends) && (partial_allele_starts <= 4);

    // merge partial alleles
    let mut pairs_of_alleles_to_merge = Vec::new();
    collect_partial_alleles_to_merge(
        &allele_match,
        assembled_chr4_allele,
        assembled_chr10_allele,
        partial_allele_number_match,
        all_starts_hap_backgrounds,
        all_ends_hap_backgrounds,
        &mut pairs_of_alleles_to_merge,
    )?;

    for alleles_to_merge in pairs_of_alleles_to_merge {
        let new_allele_summary = get_allele_summary(
            &alleles_to_merge,
            all_starts_hap_backgrounds,
            all_ends_hap_backgrounds,
            &mut distal_alleles_handled,
            &mut proximal_alleles_handled,
            fp_graph,
            fp_info,
            &methyl_values,
            qal_alleles,
        )?;
        merged_allele_summary.push(new_allele_summary);
    }

    // handle scenario where there are two pairs of proximal-distal partial alleles in each group
    // try to merge a pair by checking overlap
    // get chromosome information if both proximal alleles are on the same chromosome
    // add the size of the shortest proximal allele
    for allele_type in ["qB", "qAIntactPolyA", "qADisruptedPolyA"] {
        if let Some(this_type_alleles) = allele_match.get(allele_type) {
            if this_type_alleles.len() == 4 {
                let left_flanks = this_type_alleles
                    .iter()
                    .filter(|x| x.starts_with("LeftFlank") && !x.contains("RightFlank"))
                    .map(|x| x.to_string())
                    .collect::<Vec<String>>();
                let right_flanks = this_type_alleles
                    .iter()
                    .filter(|x| x.ends_with("RightFlank"))
                    .map(|x| x.to_string())
                    .collect::<Vec<String>>();
                if left_flanks.len() == 2 && right_flanks.len() == 2 {
                    debug!("evaluating four partial alleles of the same type {allele_type}: {left_flanks:?} and {right_flanks:?}");
                    // check if we have two pairs of overlapping partial alleles
                    let left_flank1 = left_flanks.first().unwrap().clone();
                    let left_flank2 = left_flanks.last().unwrap().clone();
                    let right_flank1 = right_flanks.first().unwrap().clone();
                    let right_flank2 = right_flanks.last().unwrap().clone();

                    let (ovl_len1, allele_name1, allele_size1) = merge_two_partial_alleles(
                        &vec![left_flank1.clone(), right_flank1.clone()],
                        fp_graph,
                        fp_info,
                    );
                    let (ovl_len2, allele_name2, allele_size2) = merge_two_partial_alleles(
                        &vec![left_flank2.clone(), right_flank2.clone()],
                        fp_graph,
                        fp_info,
                    );
                    let (ovl_len3, allele_name3, allele_size3) = merge_two_partial_alleles(
                        &vec![left_flank1.clone(), right_flank2.clone()],
                        fp_graph,
                        fp_info,
                    );
                    let (ovl_len4, allele_name4, allele_size4) = merge_two_partial_alleles(
                        &vec![left_flank2.clone(), right_flank1.clone()],
                        fp_graph,
                        fp_info,
                    );

                    if ovl_len1 >= 2 && ovl_len2 >= 2 && ovl_len3 == 0 && ovl_len4 == 0 {
                        // merge left_flank1 and right_flank1, merge left_flank2 and right_flank2
                        if !right_flank1.contains("RightFlank")
                            && !right_flank1.contains("LeftFlank")
                        {
                            let chr_info = all_starts_hap_backgrounds
                                .get(&left_flank1)
                                .unwrap()
                                .clone();
                            let background =
                                all_ends_hap_backgrounds.get(&right_flank1).unwrap().clone();
                            let methylation_value =
                                get_methylation_value(&right_flank1, &methyl_values, None)?;
                            let ending_in_qAL = qal_alleles.contains(&right_flank1);
                            debug!(
                                "adding merged allele from checking 4-way overlaps: {allele_name1}"
                            );
                            merged_allele_summary.push(AlleleSummary {
                                allele_name: allele_name1,
                                chromosome: chr_info
                                    .clone()
                                    .replace("chromosome_unknown", "unknown"),
                                distal_haplotype: background,
                                allele_type: String::from("merged"),
                                allele_size: allele_size1,
                                methylation: methylation_value,
                                ending_in_qAL: ending_in_qAL,
                            });
                            distal_alleles_handled.push(right_flank1.clone());
                            proximal_alleles_handled.push(left_flank1.clone());
                        }
                        if !right_flank2.contains("RightFlank")
                            && !right_flank2.contains("LeftFlank")
                        {
                            let chr_info = all_starts_hap_backgrounds
                                .get(&left_flank2)
                                .unwrap()
                                .clone();
                            let background =
                                all_ends_hap_backgrounds.get(&right_flank2).unwrap().clone();
                            let methylation_value =
                                get_methylation_value(&right_flank2, &methyl_values, None)?;
                            let ending_in_qAL = qal_alleles.contains(&right_flank2);
                            debug!(
                                "adding merged allele from checking 4-way overlaps: {allele_name2}"
                            );
                            merged_allele_summary.push(AlleleSummary {
                                allele_name: allele_name2,
                                chromosome: chr_info
                                    .clone()
                                    .replace("chromosome_unknown", "unknown"),
                                distal_haplotype: background,
                                allele_type: String::from("merged"),
                                allele_size: allele_size2,
                                methylation: methylation_value,
                                ending_in_qAL: ending_in_qAL,
                            });
                            distal_alleles_handled.push(right_flank2.clone());
                            proximal_alleles_handled.push(left_flank2.clone());
                        }
                    } else if ovl_len1 == 0 && ovl_len2 == 0 && ovl_len3 >= 2 && ovl_len4 >= 2 {
                        // merge left_flank1 and right_flank2, merge left_flank2 and right_flank1
                        if !right_flank2.contains("RightFlank")
                            && !right_flank2.contains("LeftFlank")
                        {
                            let chr_info = all_starts_hap_backgrounds
                                .get(&left_flank1)
                                .unwrap()
                                .clone();
                            let background =
                                all_ends_hap_backgrounds.get(&right_flank2).unwrap().clone();
                            let methylation_value =
                                get_methylation_value(&right_flank2, &methyl_values, None)?;
                            let ending_in_qAL = qal_alleles.contains(&right_flank2);
                            debug!(
                                "adding merged allele from checking 4-way overlaps: {allele_name3}"
                            );
                            merged_allele_summary.push(AlleleSummary {
                                allele_name: allele_name3,
                                chromosome: chr_info
                                    .clone()
                                    .replace("chromosome_unknown", "unknown"),
                                distal_haplotype: background,
                                allele_type: String::from("merged"),
                                allele_size: allele_size3,
                                methylation: methylation_value,
                                ending_in_qAL: ending_in_qAL,
                            });
                            distal_alleles_handled.push(right_flank2.clone());
                            proximal_alleles_handled.push(left_flank1.clone());
                        }
                        if !right_flank1.contains("RightFlank")
                            && !right_flank1.contains("LeftFlank")
                        {
                            let chr_info = all_starts_hap_backgrounds
                                .get(&left_flank2)
                                .unwrap()
                                .clone();
                            let background =
                                all_ends_hap_backgrounds.get(&right_flank1).unwrap().clone();
                            let methylation_value =
                                get_methylation_value(&right_flank1, &methyl_values, None)?;
                            let ending_in_qAL = qal_alleles.contains(&right_flank1);
                            debug!(
                                "adding merged allele from checking 4-way overlaps: {allele_name4}"
                            );
                            merged_allele_summary.push(AlleleSummary {
                                allele_name: allele_name4,
                                chromosome: chr_info
                                    .clone()
                                    .replace("chromosome_unknown", "unknown"),
                                distal_haplotype: background,
                                allele_type: String::from("merged"),
                                allele_size: allele_size4,
                                methylation: methylation_value,
                                ending_in_qAL: ending_in_qAL,
                            });
                            distal_alleles_handled.push(right_flank1.clone());
                            proximal_alleles_handled.push(left_flank2.clone());
                        }
                    }
                    // add the size of the shorter of the two left flanks
                    let left_flank_size_short = left_flanks
                        .iter()
                        .map(|x| x.split("-").filter(|x| !x.contains("Flank")).count())
                        .min()
                        .unwrap_or(0);
                    debug!(
                        "shorter of the two proximal alleles of this type {allele_type}: {left_flank_size_short}"
                    );
                    // get the last node of the shorter of the two left flanks
                    let mut left_flank_last_nodes = Vec::new();
                    for left_flank in &left_flanks {
                        let left_flank_size = left_flank
                            .split("-")
                            .filter(|x| !x.contains("Flank"))
                            .count();
                        if left_flank_size == left_flank_size_short {
                            let last_node = left_flank.split("-").last().unwrap();
                            left_flank_last_nodes.push(last_node);
                        }
                    }

                    // check which chromosome the two left flanks are on
                    let mut infered_chromosome = String::from("unknown");
                    let left_flank_chromosomes = left_flanks
                        .iter()
                        .map(|x| all_starts_hap_backgrounds.get(x))
                        .filter(|x| x.is_some())
                        .map(|x| x.unwrap().clone())
                        .collect::<Vec<String>>();
                    debug!("chromosomes of the two proximal alleles of this type {allele_type}: {left_flank_chromosomes:?}");
                    if left_flank_chromosomes.len() == 2 {
                        let left_flank_chromosome_1 = left_flank_chromosomes.first().unwrap();
                        let left_flank_chromosome_2 = left_flank_chromosomes.last().unwrap();
                        if left_flank_chromosome_1.contains("chr4")
                            && left_flank_chromosome_2.contains("chr4")
                        {
                            infered_chromosome =
                                if left_flank_chromosome_1 == left_flank_chromosome_2 {
                                    left_flank_chromosome_1.clone()
                                } else {
                                    String::from("chr4:upstream_group_unknown")
                                };
                        } else if left_flank_chromosome_1.contains("chr10")
                            && left_flank_chromosome_2.contains("chr10")
                        {
                            infered_chromosome =
                                if left_flank_chromosome_1 == left_flank_chromosome_2 {
                                    left_flank_chromosome_1.clone()
                                } else {
                                    String::from("chr10:upstream_group_unknown")
                                };
                        }
                    }

                    for allele in &right_flanks {
                        if !distal_alleles_handled.contains(allele) {
                            // we can add the shorter of the two left flanks
                            let allele_size =
                                allele.split("-").filter(|x| !x.contains("Flank")).count()
                                    + left_flank_size_short;
                            let this_allele_first_node = allele.split("-").next().unwrap();
                            let background = all_ends_hap_backgrounds.get(allele).unwrap().clone();
                            let methylation_value =
                                get_methylation_value(allele, &methyl_values, None)?;
                            let ending_in_qAL = qal_alleles.contains(allele);
                            if left_flank_last_nodes.contains(&this_allele_first_node) {
                                merged_allele_summary.push(AlleleSummary {
                                    allele_name: allele.clone(),
                                    chromosome: infered_chromosome.clone(),
                                    distal_haplotype: background,
                                    allele_type: String::from("partial"),
                                    allele_size: format!(">={allele_size}"),
                                    methylation: methylation_value,
                                    ending_in_qAL: ending_in_qAL,
                                });
                            } else {
                                merged_allele_summary.push(AlleleSummary {
                                    allele_name: allele.clone(),
                                    chromosome: infered_chromosome.clone(),
                                    distal_haplotype: background,
                                    allele_type: String::from("partial"),
                                    allele_size: format!(">{allele_size}"),
                                    methylation: methylation_value,
                                    ending_in_qAL: ending_in_qAL,
                                });
                            }
                            distal_alleles_handled.push(allele.clone());
                        }
                    }
                    for allele in &left_flanks {
                        proximal_alleles_handled.push(allele.clone());
                    }
                }
            }
        }
    }

    // here we can do some cleanup to remove redundant alleles
    let mut remaining_distal_alleles = all_ends_hap_backgrounds
        .iter()
        .filter(|(k, _v)| !distal_alleles_handled.contains(k))
        .map(|(k, v)| k.clone())
        .filter(|x| !is_cis_dup(x, fp_info).unwrap_or(false))
        .collect::<Vec<String>>();
    let mut remaining_proximal_alleles = all_starts_hap_backgrounds
        .iter()
        .filter(|(k, _v)| !proximal_alleles_handled.contains(k))
        .map(|(k, v)| k.clone())
        .collect::<Vec<String>>();
    /*
    if remaining_proximal_alleles.len() > remaining_distal_alleles.len()
        && all_starts_hap_backgrounds.len() >= 4
    {
        // remove redundant proximal alleles
        remaining_proximal_alleles = remove_redundant_alleles(&remaining_proximal_alleles, 2)?;
    } else if remaining_proximal_alleles.len() < remaining_distal_alleles.len()
        && all_ends_hap_backgrounds.len() >= 4
    {
        // remove redundant distal alleles
        remaining_distal_alleles = remove_redundant_alleles(&remaining_distal_alleles, 2)?;
    }
    */

    // for the remaining distal alleles, we want to see if there is any information on the chromosome
    // or the size that we can update

    // get the min size of all proximal alleles
    let start_allele_chr4 = all_starts_hap_backgrounds
        .iter()
        .filter(|(_, background)| background.contains("chr4"))
        .count();
    let start_allele_chr10 = all_starts_hap_backgrounds
        .iter()
        .filter(|(_, background)| background.contains("chr10"))
        .count();
    let all_start_min_size = if start_allele_chr4 >= 2 && start_allele_chr10 >= 2 {
        all_starts_hap_backgrounds
            .iter()
            .filter(|(k, _v)| !k.ends_with("RightFlank"))
            .map(|(k, _v)| k.split("-").filter(|x| !x.contains("Flank")).count())
            .min()
            .unwrap_or(0)
    } else {
        0
    };
    debug!("all_start_min_size: {all_start_min_size}");
    // get the last node of the shorter of all proximal alleles
    let mut all_start_alleles_last_nodes = Vec::new();
    for (allele, _background) in all_starts_hap_backgrounds.iter() {
        let allele_size = allele.split("-").filter(|x| !x.contains("Flank")).count();
        if allele_size == all_start_min_size && !allele.ends_with("RightFlank") {
            let allele_last_node = allele.split("-").last().unwrap();
            all_start_alleles_last_nodes.push(allele_last_node);
        }
    }
    // check if the remaining proximal alleles are on the same chromosome
    let mut remaining_proximal_allele_chromosome = String::from("unknown");
    if all_starts_hap_backgrounds.len() == 4
        && all_ends_hap_backgrounds.len() == 4
        && start_allele_chr4 == 2
        && start_allele_chr10 == 2
        && remaining_proximal_alleles.len() == remaining_distal_alleles.len()
        && remaining_distal_alleles.len() <= 2
    {
        let remaining_proximal_alleles_chromosomes = all_starts_hap_backgrounds
            .iter()
            .filter(|(k, _v)| remaining_proximal_alleles.contains(k))
            .map(|(k, v)| v.clone())
            .collect::<Vec<String>>();
        debug!("remaining proximal alleles are on the following chromosomes: {remaining_proximal_alleles_chromosomes:?}");
        if remaining_proximal_alleles_chromosomes.len() == 1 {
            remaining_proximal_allele_chromosome = remaining_proximal_alleles_chromosomes
                .first()
                .unwrap()
                .clone()
                .replace("chromosome_unknown", "unknown");
        } else if remaining_proximal_alleles_chromosomes.len() == 2 {
            let remaining_proximal_allele_chromosome_1 = remaining_proximal_alleles_chromosomes
                .first()
                .unwrap()
                .clone();
            let remaining_proximal_allele_chromosome_2 = remaining_proximal_alleles_chromosomes
                .last()
                .unwrap()
                .clone();
            if remaining_proximal_allele_chromosome_1.contains("chr4")
                && remaining_proximal_allele_chromosome_2.contains("chr4")
            {
                if remaining_proximal_allele_chromosome_1 == remaining_proximal_allele_chromosome_2
                {
                    remaining_proximal_allele_chromosome = remaining_proximal_allele_chromosome_1
                        .replace("chromosome_unknown", "unknown");
                } else {
                    remaining_proximal_allele_chromosome =
                        String::from("chr4:upstream_group_unknown");
                }
            } else if remaining_proximal_allele_chromosome_1.contains("chr10")
                && remaining_proximal_allele_chromosome_2.contains("chr10")
            {
                if remaining_proximal_allele_chromosome_1 == remaining_proximal_allele_chromosome_2
                {
                    remaining_proximal_allele_chromosome = remaining_proximal_allele_chromosome_1
                        .replace("chromosome_unknown", "unknown");
                } else {
                    remaining_proximal_allele_chromosome =
                        String::from("chr10:upstream_group_unknown");
                }
            }
        }
    }

    // add remaining distall alleles
    for (allele, background) in all_ends_hap_backgrounds.iter() {
        let ending_in_qAL = qal_alleles.contains(allele);
        if !distal_alleles_handled.contains(allele) {
            let methylation_value = get_methylation_value(allele, &methyl_values, None)?;
            let allele_size = allele.split("-").filter(|x| !x.contains("Flank")).count();
            if is_cis_dup(allele, fp_info)? {
                merged_allele_summary.push(AlleleSummary {
                    allele_name: allele.clone(),
                    chromosome: String::from("unknown"),
                    distal_haplotype: background.clone(),
                    allele_type: String::from("assembled_cis_duplication"),
                    allele_size: format!("{allele_size}"),
                    methylation: methylation_value,
                    ending_in_qAL: ending_in_qAL,
                });
            } else {
                if ((background == "qAIntactPolyA" && start_allele_chr4 >= 2)
                        || (background == "qB" && start_allele_chr4 >= 2)
                        || (background == "qADisruptedPolyA" && start_allele_chr10 >= 2))
                        // && allele_size <= 10   // TODO: evaluate if we want to do this for all partial alleles
                        && all_start_min_size > 0
                        && all_ends_hap_backgrounds.len() <= 4 // we could have cis dup alleles and size shouldn't be added
                        && !allele.starts_with("LeftFlank")
                        && !allele.starts_with("RightFlank")
                {
                    debug!("adding size {all_start_min_size} to partial allele {allele} {background} with size {allele_size}");
                    let allele_size = allele_size + all_start_min_size;
                    let this_allele_first_node = allele.split("-").next().unwrap();
                    if all_start_alleles_last_nodes.contains(&this_allele_first_node) {
                        merged_allele_summary.push(AlleleSummary {
                            allele_name: allele.clone(),
                            chromosome: remaining_proximal_allele_chromosome.clone(),
                            distal_haplotype: background.clone(),
                            allele_type: String::from("partial"),
                            allele_size: format!(">={allele_size}"),
                            methylation: methylation_value,
                            ending_in_qAL: ending_in_qAL,
                        });
                    } else {
                        merged_allele_summary.push(AlleleSummary {
                            allele_name: allele.clone(),
                            chromosome: remaining_proximal_allele_chromosome.clone(),
                            distal_haplotype: background.clone(),
                            allele_type: String::from("partial"),
                            allele_size: format!(">{allele_size}"),
                            methylation: methylation_value,
                            ending_in_qAL: ending_in_qAL,
                        });
                    }
                } else if allele.starts_with("LeftFlank") {
                    let chromosome_info = if all_starts_hap_backgrounds.contains_key(allele) {
                        all_starts_hap_backgrounds
                            .get(allele)
                            .unwrap()
                            .clone()
                            .replace("chromosome_unknown", "unknown")
                    } else {
                        remaining_proximal_allele_chromosome.clone()
                    };
                    merged_allele_summary.push(AlleleSummary {
                        allele_name: allele.clone(),
                        chromosome: chromosome_info.clone(),
                        distal_haplotype: background.clone(),
                        allele_type: String::from("partial"),
                        allele_size: format!(">={allele_size}"),
                        methylation: methylation_value,
                        ending_in_qAL: ending_in_qAL,
                    });
                } else {
                    merged_allele_summary.push(AlleleSummary {
                        allele_name: allele.clone(),
                        chromosome: remaining_proximal_allele_chromosome.clone(),
                        distal_haplotype: background.clone(),
                        allele_type: String::from("partial"),
                        allele_size: format!(">{allele_size}"),
                        methylation: methylation_value,
                        ending_in_qAL: ending_in_qAL,
                    });
                }
            }
        }
    }
    Ok(merged_allele_summary)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn create_empty_backgrounds() -> (BTreeMap<String, String>, BTreeMap<String, String>) {
        (BTreeMap::new(), BTreeMap::new())
    }

    #[test]
    fn test_get_methylation_value_numeric_values() {
        // Function uses last 505 sites; need at least 505 to avoid underflow. Use 600 values.
        let values: Vec<String> = (0..600).map(|_| "0.5".to_string()).collect();
        let mut methyl_values = BTreeMap::new();
        methyl_values.insert("allele-1".to_string(), values.join(","));
        let allele = "allele-1".to_string();
        let result = get_methylation_value(&allele, &methyl_values, None).unwrap();
        assert!(!result.is_nan());
        assert!(
            (result - 0.5).abs() < 0.001,
            "expected ~0.5, got {}",
            result
        );
    }

    #[test]
    fn test_get_methylation_value_all_nan_returns_nan() {
        // NaN is parsed as nan and excluded from median; no values remain -> return NaN.
        let values: Vec<String> = (0..505).map(|_| "NaN".to_string()).collect();
        let mut methyl_values = BTreeMap::new();
        methyl_values.insert("allele-nan".to_string(), values.join(","));
        let allele = "allele-nan".to_string();
        let result = get_methylation_value(&allele, &methyl_values, None).unwrap();
        assert!(
            result.is_nan(),
            "all NaN with no numeric values should return NaN"
        );
    }

    #[test]
    fn test_get_methylation_value_mixed_nan_and_numeric() {
        let parts = String::from("0.5,0.5,0.1,0.2,0.3,0.4,0.5");
        let mut methyl_values = BTreeMap::new();
        methyl_values.insert("allele-mixed".to_string(), parts);
        let allele = "allele-mixed".to_string();
        let result = get_methylation_value(&allele, &methyl_values, Some(5)).unwrap();
        assert!(!result.is_nan());
        assert!(
            (result - 0.3).abs() < 0.001,
            "expected ~0.3, got {}",
            result
        );

        let parts = String::from("0.5,0.5,0.1,0.2,0.3,NaN,0.4,0.5");
        let mut methyl_values = BTreeMap::new();
        methyl_values.insert("allele-mixed".to_string(), parts);
        let allele = "allele-mixed".to_string();
        let result = get_methylation_value(&allele, &methyl_values, Some(5)).unwrap();
        assert!(!result.is_nan());
        assert!(
            (result - 0.35).abs() < 0.001,
            "expected ~0.35, got {}",
            result
        );
    }

    #[test]
    fn test_collect_partial_alleles_no_unknown_with_flanks() {
        let mut allele_match = BTreeMap::new();
        allele_match.insert(
            "qAIntactPolyA".to_string(),
            vec![
                "LeftFlank-1-2-3".to_string(),
                "4-5-6-RightFlank".to_string(),
            ],
        );
        allele_match.insert(
            "qB".to_string(),
            vec!["LeftFlank-7-8".to_string(), "9-10-RightFlank".to_string()],
        );

        let (starts, ends) = create_empty_backgrounds();
        let mut pairs_to_merge = Vec::new();

        collect_partial_alleles_to_merge(
            &allele_match,
            0,
            0,
            false,
            &starts,
            &ends,
            &mut pairs_to_merge,
        )
        .unwrap();

        assert_eq!(pairs_to_merge.len(), 2);
        assert!(pairs_to_merge.contains(&vec![
            "LeftFlank-1-2-3".to_string(),
            "4-5-6-RightFlank".to_string()
        ]));
        assert!(pairs_to_merge.contains(&vec![
            "LeftFlank-7-8".to_string(),
            "9-10-RightFlank".to_string()
        ]));
    }

    #[test]
    fn test_collect_partial_alleles_no_unknown_without_flanks() {
        let mut allele_match = BTreeMap::new();
        allele_match.insert(
            "qAIntactPolyA".to_string(),
            vec!["1-2-3".to_string(), "4-5-6".to_string()],
        );

        let (starts, ends) = create_empty_backgrounds();
        let mut pairs_to_merge = Vec::new();

        collect_partial_alleles_to_merge(
            &allele_match,
            0,
            0,
            false,
            &starts,
            &ends,
            &mut pairs_to_merge,
        )
        .unwrap();

        assert_eq!(pairs_to_merge.len(), 0);
    }

    #[test]
    fn test_collect_partial_alleles_no_unknown_single_allele() {
        let mut allele_match = BTreeMap::new();
        allele_match.insert(
            "qAIntactPolyA".to_string(),
            vec!["LeftFlank-1-2-3".to_string()],
        );

        let (starts, ends) = create_empty_backgrounds();
        let mut pairs_to_merge = Vec::new();

        collect_partial_alleles_to_merge(
            &allele_match,
            0,
            0,
            false,
            &starts,
            &ends,
            &mut pairs_to_merge,
        )
        .unwrap();

        assert_eq!(pairs_to_merge.len(), 0);
    }

    #[test]
    fn test_collect_partial_alleles_assembled_chr4_merge_qb() {
        let mut allele_match = BTreeMap::new();
        allele_match.insert(
            "qB".to_string(),
            vec!["LeftFlank-1-2".to_string(), "3-4-RightFlank".to_string()],
        );
        allele_match.insert(
            "qADisruptedPolyA".to_string(),
            vec!["LeftFlank-5-6".to_string(), "7-8-RightFlank".to_string()],
        );
        allele_match.insert("unknown".to_string(), vec!["unknown-allele".to_string()]);

        let (starts, ends) = create_empty_backgrounds();
        let mut pairs_to_merge = Vec::new();

        collect_partial_alleles_to_merge(
            &allele_match,
            1, // assembled_chr4_allele == 1
            0,
            false,
            &starts,
            &ends,
            &mut pairs_to_merge,
        )
        .unwrap();

        assert_eq!(pairs_to_merge.len(), 1);
        assert!(pairs_to_merge.contains(&vec![
            "LeftFlank-1-2".to_string(),
            "3-4-RightFlank".to_string()
        ]));
        // qAIntactPolyA should not be merged when assembled_chr4_allele == 1
        assert!(!pairs_to_merge.contains(&vec![
            "LeftFlank-5-6".to_string(),
            "7-8-RightFlank".to_string()
        ]));
    }

    #[test]
    fn test_collect_partial_alleles_assembled_chr4_merge_qaintactpolya() {
        let mut allele_match = BTreeMap::new();
        allele_match.insert(
            "qAIntactPolyA".to_string(),
            vec!["LeftFlank-1-2".to_string(), "3-4-RightFlank".to_string()],
        );
        allele_match.insert("unknown".to_string(), vec!["unknown-allele".to_string()]);

        let (starts, ends) = create_empty_backgrounds();
        let mut pairs_to_merge = Vec::new();

        collect_partial_alleles_to_merge(
            &allele_match,
            1, // assembled_chr4_allele == 1
            0,
            false,
            &starts,
            &ends,
            &mut pairs_to_merge,
        )
        .unwrap();

        assert_eq!(pairs_to_merge.len(), 1);
        assert!(pairs_to_merge.contains(&vec![
            "LeftFlank-1-2".to_string(),
            "3-4-RightFlank".to_string()
        ]));
    }

    #[test]
    fn test_collect_partial_alleles_unknown_distal_qaintactpolya() {
        let mut allele_match = BTreeMap::new();
        let unknown_allele = "unknown-allele".to_string();
        allele_match.insert("unknown".to_string(), vec![unknown_allele.clone()]);
        allele_match.insert(
            "qAIntactPolyA".to_string(),
            vec!["LeftFlank-1-2-3".to_string()],
        );

        let starts = BTreeMap::new();
        let mut ends = BTreeMap::new();
        ends.insert(unknown_allele.clone(), "qAIntactPolyA".to_string());

        let mut pairs_to_merge = Vec::new();

        collect_partial_alleles_to_merge(
            &allele_match,
            0,
            0,
            true, // partial_allele_number_match
            &starts,
            &ends,
            &mut pairs_to_merge,
        )
        .unwrap();

        assert_eq!(pairs_to_merge.len(), 1);
        assert_eq!(
            pairs_to_merge[0],
            vec!["LeftFlank-1-2-3".to_string(), unknown_allele]
        );
    }

    #[test]
    fn test_collect_partial_alleles_unknown_distal_qb() {
        let mut allele_match = BTreeMap::new();
        let unknown_allele = "unknown-allele".to_string();
        allele_match.insert("unknown".to_string(), vec![unknown_allele.clone()]);
        allele_match.insert("qB".to_string(), vec!["LeftFlank-1-2".to_string()]);

        let starts = BTreeMap::new();
        let mut ends = BTreeMap::new();
        ends.insert(unknown_allele.clone(), "qB".to_string());

        let mut pairs_to_merge = Vec::new();

        collect_partial_alleles_to_merge(
            &allele_match,
            0,
            0,
            true, // partial_allele_number_match
            &starts,
            &ends,
            &mut pairs_to_merge,
        )
        .unwrap();

        assert_eq!(pairs_to_merge.len(), 1);
        assert_eq!(
            pairs_to_merge[0],
            vec!["LeftFlank-1-2".to_string(), unknown_allele]
        );
    }

    #[test]
    fn test_collect_partial_alleles_unknown_distal_no_match() {
        let mut allele_match = BTreeMap::new();
        let unknown_allele = "unknown-allele".to_string();
        allele_match.insert("unknown".to_string(), vec![unknown_allele.clone()]);
        allele_match.insert(
            "qAIntactPolyA".to_string(),
            vec!["LeftFlank-1-2-3".to_string()],
        );

        let starts = BTreeMap::new();
        let mut ends = BTreeMap::new();
        ends.insert(unknown_allele, "other".to_string());

        let mut pairs_to_merge = Vec::new();

        collect_partial_alleles_to_merge(
            &allele_match,
            0,
            0,
            true, // partial_allele_number_match
            &starts,
            &ends,
            &mut pairs_to_merge,
        )
        .unwrap();

        assert_eq!(pairs_to_merge.len(), 0);
    }

    #[test]
    fn test_collect_partial_alleles_unknown_proximal_chr4_qaintactpolya() {
        let mut allele_match = BTreeMap::new();
        let unknown_allele = "unknown-allele".to_string();
        allele_match.insert("unknown".to_string(), vec![unknown_allele.clone()]);
        allele_match.insert(
            "qAIntactPolyA".to_string(),
            vec!["1-2-3-RightFlank".to_string()],
        );

        let mut starts = BTreeMap::new();
        let ends = BTreeMap::new();
        starts.insert(unknown_allele.clone(), "chr4".to_string());

        let mut pairs_to_merge = Vec::new();

        collect_partial_alleles_to_merge(
            &allele_match,
            1, // assembled_chr4_allele == 1, no qB
            0,
            true, // partial_allele_number_match
            &starts,
            &ends,
            &mut pairs_to_merge,
        )
        .unwrap();

        assert_eq!(pairs_to_merge.len(), 1);
        assert_eq!(
            pairs_to_merge[0],
            vec![unknown_allele, "1-2-3-RightFlank".to_string()]
        );
    }

    #[test]
    fn test_collect_partial_alleles_unknown_proximal_chr4_qb() {
        let mut allele_match = BTreeMap::new();
        let unknown_allele = "unknown-allele".to_string();
        allele_match.insert("unknown".to_string(), vec![unknown_allele.clone()]);
        allele_match.insert("qB".to_string(), vec!["1-2-RightFlank".to_string()]);

        let mut starts = BTreeMap::new();
        let ends = BTreeMap::new();
        starts.insert(unknown_allele.clone(), "chr4".to_string());

        let mut pairs_to_merge = Vec::new();

        collect_partial_alleles_to_merge(
            &allele_match,
            1, // assembled_chr4_allele == 1, no qAIntactPolyA
            0,
            true, // partial_allele_number_match
            &starts,
            &ends,
            &mut pairs_to_merge,
        )
        .unwrap();

        assert_eq!(pairs_to_merge.len(), 1);
        assert_eq!(
            pairs_to_merge[0],
            vec![unknown_allele, "1-2-RightFlank".to_string()]
        );
    }

    #[test]
    fn test_collect_partial_alleles_unknown_proximal_not_chr4() {
        let mut allele_match = BTreeMap::new();
        let unknown_allele = "unknown-allele".to_string();
        allele_match.insert("unknown".to_string(), vec![unknown_allele.clone()]);
        allele_match.insert(
            "qAIntactPolyA".to_string(),
            vec!["1-2-3-RightFlank".to_string()],
        );

        let mut starts = BTreeMap::new();
        let ends = BTreeMap::new();
        starts.insert(unknown_allele, "not-chr4".to_string());

        let mut pairs_to_merge = Vec::new();

        collect_partial_alleles_to_merge(
            &allele_match,
            0,
            0,
            true, // partial_allele_number_match
            &starts,
            &ends,
            &mut pairs_to_merge,
        )
        .unwrap();

        assert_eq!(pairs_to_merge.len(), 0);
    }

    #[test]
    fn test_collect_partial_alleles_unknown_proximal_chr4_with_qb_pair() {
        let mut allele_match = BTreeMap::new();
        let unknown_allele = "unknown-allele".to_string();
        allele_match.insert("unknown".to_string(), vec![unknown_allele.clone()]);
        allele_match.insert(
            "qAIntactPolyA".to_string(),
            vec!["1-2-3-RightFlank".to_string()],
        );
        allele_match.insert(
            "qB".to_string(),
            vec!["LeftFlank-4-5".to_string(), "6-7-RightFlank".to_string()],
        );

        let mut starts = BTreeMap::new();
        let ends = BTreeMap::new();
        starts.insert(unknown_allele.clone(), "chr4".to_string());

        let mut pairs_to_merge = Vec::new();

        collect_partial_alleles_to_merge(
            &allele_match,
            0,
            0,
            true, // partial_allele_number_match
            &starts,
            &ends,
            &mut pairs_to_merge,
        )
        .unwrap();

        // Should merge unknown with qAIntactPolyA, and also merge the qB pair
        assert_eq!(pairs_to_merge.len(), 2);
        assert!(pairs_to_merge.contains(&vec![unknown_allele, "1-2-3-RightFlank".to_string()]));
        assert!(pairs_to_merge.contains(&vec![
            "LeftFlank-4-5".to_string(),
            "6-7-RightFlank".to_string()
        ]));
    }

    #[test]
    fn test_collect_partial_alleles_unknown_multiple_unknown_alleles() {
        let mut allele_match = BTreeMap::new();
        allele_match.insert(
            "unknown".to_string(),
            vec!["unknown1".to_string(), "unknown2".to_string()],
        );

        let (starts, ends) = create_empty_backgrounds();
        let mut pairs_to_merge = Vec::new();

        collect_partial_alleles_to_merge(
            &allele_match,
            0,
            0,
            true, // partial_allele_number_match
            &starts,
            &ends,
            &mut pairs_to_merge,
        )
        .unwrap();

        // Should not merge when there are multiple unknown alleles
        assert_eq!(pairs_to_merge.len(), 0);
    }

    #[test]
    fn test_collect_partial_alleles_unknown_distal_no_leftflank() {
        let mut allele_match = BTreeMap::new();
        let unknown_allele = "unknown-allele".to_string();
        allele_match.insert("unknown".to_string(), vec![unknown_allele.clone()]);
        // qAIntactPolyA allele without LeftFlank
        allele_match.insert("qAIntactPolyA".to_string(), vec!["1-2-3".to_string()]);

        let starts = BTreeMap::new();
        let mut ends = BTreeMap::new();
        ends.insert(unknown_allele, "qAIntactPolyA".to_string());

        let mut pairs_to_merge = Vec::new();

        collect_partial_alleles_to_merge(
            &allele_match,
            0,
            0,
            true, // partial_allele_number_match
            &starts,
            &ends,
            &mut pairs_to_merge,
        )
        .unwrap();

        // Should not merge because qAIntactPolyA doesn't have LeftFlank
        assert_eq!(pairs_to_merge.len(), 0);
    }

    #[test]
    fn test_collect_partial_alleles_unknown_proximal_no_rightflank() {
        let mut allele_match = BTreeMap::new();
        let unknown_allele = "unknown-allele".to_string();
        allele_match.insert("unknown".to_string(), vec![unknown_allele.clone()]);
        // qAIntactPolyA allele without RightFlank
        allele_match.insert("qAIntactPolyA".to_string(), vec!["1-2-3".to_string()]);

        let mut starts = BTreeMap::new();
        let ends = BTreeMap::new();
        starts.insert(unknown_allele, "chr4".to_string());

        let mut pairs_to_merge = Vec::new();

        collect_partial_alleles_to_merge(
            &allele_match,
            1,
            0,
            true, // partial_allele_number_match
            &starts,
            &ends,
            &mut pairs_to_merge,
        )
        .unwrap();

        // Should not merge because qAIntactPolyA doesn't have RightFlank
        assert_eq!(pairs_to_merge.len(), 0);
    }

    #[test]
    fn test_collect_partial_alleles_unknown_proximal_qaintactpolya_multiple_alleles() {
        let mut allele_match = BTreeMap::new();
        let unknown_allele = "unknown-allele".to_string();
        allele_match.insert("unknown".to_string(), vec![unknown_allele.clone()]);
        // qAIntactPolyA with multiple alleles (should not merge)
        allele_match.insert(
            "qAIntactPolyA".to_string(),
            vec![
                "1-2-3-RightFlank".to_string(),
                "4-5-6-RightFlank".to_string(),
            ],
        );

        let mut starts = BTreeMap::new();
        let ends = BTreeMap::new();
        starts.insert(unknown_allele, "chr4".to_string());

        let mut pairs_to_merge = Vec::new();

        collect_partial_alleles_to_merge(
            &allele_match,
            1,
            0,
            true, // partial_allele_number_match
            &starts,
            &ends,
            &mut pairs_to_merge,
        )
        .unwrap();

        // Should not merge because qAIntactPolyA has multiple alleles
        assert_eq!(pairs_to_merge.len(), 0);
    }

    #[test]
    fn test_collect_partial_alleles_empty_allele_match() {
        let allele_match = BTreeMap::new();
        let (starts, ends) = create_empty_backgrounds();
        let mut pairs_to_merge = Vec::new();

        collect_partial_alleles_to_merge(
            &allele_match,
            0,
            0,
            false,
            &starts,
            &ends,
            &mut pairs_to_merge,
        )
        .unwrap();

        assert_eq!(pairs_to_merge.len(), 0);
    }

    fn create_test_region_coordinates() -> RegionCoordinates {
        let mut variants_to_distinguish = BTreeMap::new();
        variants_to_distinguish.insert(
            "qADisruptedPolyA".to_string(),
            vec![
                "256:C>T".to_string(),
                "583:C>T".to_string(),
                "683:A>G".to_string(),
            ],
        );
        variants_to_distinguish.insert(
            "qB".to_string(),
            vec![
                "172:A>G".to_string(),
                "513:C>A".to_string(),
                "2891:T>C".to_string(),
            ],
        );

        RegionCoordinates {
            repeat_len: 3298,
            chromosome_output: String::from("d4z4_ref"),
            chromosome_len: 4203,
            genome_offset: BTreeMap::new(),
            reference_seq: String::new(),
            extract_regions: Vec::new(),
            flanking_regions: None,
            depth_regions: (1000, 3000),
            type2_sites: Vec::new(),
            exclude_sites: Vec::new(),
            exclude_sites_vcf: Vec::new(),
            genome_depth_sites: Vec::new(),
            clip_variant_sites: BTreeMap::new(),
            methyl_sites: Vec::new(),
            start_positions_flank: None,
            end_positions_flank: None,
            pivot_site: None,
            variants_to_call: Vec::new(),
            variants_to_exclude: Vec::new(),
            realign_segments: Vec::new(),
            variants_to_distinguish_allele_types: variants_to_distinguish,
        }
    }

    #[test]
    fn test_classify_fingerprint_qadisruptedpolya_two_variants() {
        let region_coords = create_test_region_coordinates();
        let variants = vec!["256:C>T".to_string(), "583:C>T".to_string()];

        let (allele_type, count_qa_disrupted, count_qb) =
            classify_fingerprint(&variants, &region_coords);

        assert_eq!(allele_type, "qADisruptedPolyA");
        assert_eq!(count_qa_disrupted, 2);
        assert_eq!(count_qb, 0);
    }

    #[test]
    fn test_classify_fingerprint_qadisruptedpolya_three_variants() {
        let region_coords = create_test_region_coordinates();
        let variants = vec![
            "256:C>T".to_string(),
            "583:C>T".to_string(),
            "683:A>G".to_string(),
        ];

        let (allele_type, count_qa_disrupted, count_qb) =
            classify_fingerprint(&variants, &region_coords);

        assert_eq!(allele_type, "qADisruptedPolyA");
        assert_eq!(count_qa_disrupted, 3);
        assert_eq!(count_qb, 0);
    }

    #[test]
    fn test_classify_fingerprint_qadisruptedpolya_with_other_variants() {
        let region_coords = create_test_region_coordinates();
        let variants = vec![
            "256:C>T".to_string(),
            "583:C>T".to_string(),
            "9999:A>T".to_string(), // not in either list
        ];

        let (allele_type, count_qa_disrupted, count_qb) =
            classify_fingerprint(&variants, &region_coords);

        assert_eq!(allele_type, "qADisruptedPolyA");
        assert_eq!(count_qa_disrupted, 2);
        assert_eq!(count_qb, 0);
    }

    #[test]
    fn test_classify_fingerprint_qb_one_variant() {
        let region_coords = create_test_region_coordinates();
        let variants = vec!["172:A>G".to_string()];

        let (allele_type, count_qa_disrupted, count_qb) =
            classify_fingerprint(&variants, &region_coords);

        assert_eq!(allele_type, "qB");
        assert_eq!(count_qb, 1);
        assert_eq!(count_qa_disrupted, 0);
    }

    #[test]
    fn test_classify_fingerprint_qb_multiple_variants() {
        let region_coords = create_test_region_coordinates();
        let variants = vec![
            "172:A>G".to_string(),
            "513:C>A".to_string(),
            "2891:T>C".to_string(),
        ];

        let (allele_type, count_qa_disrupted, count_qb) =
            classify_fingerprint(&variants, &region_coords);

        assert_eq!(allele_type, "qB");
        assert_eq!(count_qb, 3);
        assert_eq!(count_qa_disrupted, 0);
    }

    #[test]
    fn test_classify_fingerprint_qaintactpolya_no_variants() {
        let region_coords = create_test_region_coordinates();
        let variants = vec![];

        let (allele_type, count_qa_disrupted, count_qb) =
            classify_fingerprint(&variants, &region_coords);

        assert_eq!(allele_type, "qAIntactPolyA");
        assert_eq!(count_qb, 0);
        assert_eq!(count_qa_disrupted, 0);
    }

    #[test]
    fn test_classify_fingerprint_qaintactpolya_other_variants() {
        let region_coords = create_test_region_coordinates();
        let variants = vec!["9999:A>T".to_string(), "8888:G>C".to_string()];

        let (allele_type, count_qa_disrupted, count_qb) =
            classify_fingerprint(&variants, &region_coords);

        assert_eq!(allele_type, "qAIntactPolyA");
        assert_eq!(count_qb, 0);
        assert_eq!(count_qa_disrupted, 0);
    }

    #[test]
    fn test_classify_fingerprint_type4_one_qa_disrupted() {
        let region_coords = create_test_region_coordinates();
        let variants = vec!["256:C>T".to_string()];

        let (allele_type, count_qa_disrupted, count_qb) =
            classify_fingerprint(&variants, &region_coords);

        assert_eq!(allele_type, "type4");
        assert_eq!(count_qb, 0);
        assert_eq!(count_qa_disrupted, 1);
    }

    #[test]
    fn test_classify_fingerprint_type4_one_qa_disrupted_with_other() {
        let region_coords = create_test_region_coordinates();
        let variants = vec!["256:C>T".to_string(), "9999:A>T".to_string()];

        let (allele_type, count_qa_disrupted, count_qb) =
            classify_fingerprint(&variants, &region_coords);

        assert_eq!(allele_type, "type4");
        assert_eq!(count_qb, 0);
        assert_eq!(count_qa_disrupted, 1);
    }

    #[test]
    fn test_classify_fingerprint_unknown_both_present() {
        let region_coords = create_test_region_coordinates();
        let variants = vec!["256:C>T".to_string(), "172:A>G".to_string()];

        let (allele_type, count_qa_disrupted, count_qb) =
            classify_fingerprint(&variants, &region_coords);

        assert_eq!(allele_type, "unknown");
        assert_eq!(count_qb, 1);
        assert_eq!(count_qa_disrupted, 1);
    }

    #[test]
    fn test_classify_fingerprint_unknown_two_qa_disrupted_with_qb() {
        let region_coords = create_test_region_coordinates();
        let variants = vec![
            "256:C>T".to_string(),
            "583:C>T".to_string(),
            "172:A>G".to_string(),
        ];

        let (allele_type, count_qa_disrupted, count_qb) =
            classify_fingerprint(&variants, &region_coords);

        assert_eq!(allele_type, "unknown");
        assert_eq!(count_qb, 1);
        assert_eq!(count_qa_disrupted, 2);
    }

    #[test]
    fn test_classify_fingerprint_unknown_mixed_variants() {
        let region_coords = create_test_region_coordinates();
        let variants = vec![
            "256:C>T".to_string(),
            "172:A>G".to_string(),
            "9999:A>T".to_string(),
        ];

        let (allele_type, count_qa_disrupted, count_qb) =
            classify_fingerprint(&variants, &region_coords);

        assert_eq!(allele_type, "unknown");
        assert_eq!(count_qb, 1);
        assert_eq!(count_qa_disrupted, 1);
    }

    #[test]
    fn test_classify_fingerprint_qb_partial_match() {
        let region_coords = create_test_region_coordinates();
        let variants = vec!["172:A>G".to_string(), "513:C>A".to_string()];

        let (allele_type, count_qa_disrupted, count_qb) =
            classify_fingerprint(&variants, &region_coords);

        assert_eq!(allele_type, "qB");
        assert_eq!(count_qb, 2);
        assert_eq!(count_qa_disrupted, 0);
    }

    #[test]
    fn test_classify_fingerprint_qadisruptedpolya_exactly_two() {
        let region_coords = create_test_region_coordinates();
        let variants = vec!["256:C>T".to_string(), "583:C>T".to_string()];

        let (allele_type, count_qa_disrupted, count_qb) =
            classify_fingerprint(&variants, &region_coords);

        assert_eq!(allele_type, "qADisruptedPolyA");
        assert_eq!(count_qa_disrupted, 2);
        assert_eq!(count_qb, 0);
    }

    #[test]
    fn test_classify_allele_empty() {
        let fps_classified = vec![];
        let result = classify_allele(&fps_classified);
        assert_eq!(result, "unknown");
    }

    #[test]
    fn test_classify_allele_single_element() {
        let fps_classified = vec!["qAIntactPolyA".to_string()];
        let result = classify_allele(&fps_classified);
        assert_eq!(result, "unknown");
    }

    #[test]
    fn test_classify_allele_two_same() {
        let fps_classified = vec!["qAIntactPolyA".to_string(), "qAIntactPolyA".to_string()];
        let result = classify_allele(&fps_classified);
        // 2 elements, both same: 2 >= 2 * 0.8 = 1.6, so should return qAIntactPolyA
        assert_eq!(result, "qAIntactPolyA");
    }

    #[test]
    fn test_classify_allele_two_different() {
        let fps_classified = vec!["qAIntactPolyA".to_string(), "qB".to_string()];
        let result = classify_allele(&fps_classified);
        // 2 elements, different: 1 < 2 * 0.8 = 1.6, so should return unknown
        assert_eq!(result, "unknown");
    }

    #[test]
    fn test_classify_allele_three_all_same() {
        let fps_classified = vec![
            "qAIntactPolyA".to_string(),
            "qAIntactPolyA".to_string(),
            "qAIntactPolyA".to_string(),
        ];
        let result = classify_allele(&fps_classified);
        // 3 elements, all same: 3 >= 3 * 0.8 = 2.4, so should return qAIntactPolyA
        assert_eq!(result, "qAIntactPolyA");
    }

    #[test]
    fn test_classify_allele_three_two_same() {
        let fps_classified = vec![
            "qAIntactPolyA".to_string(),
            "qAIntactPolyA".to_string(),
            "qB".to_string(),
        ];
        let result = classify_allele(&fps_classified);
        // 3 elements, 2 same: 2 < 3 * 0.8 = 2.4, so should return unknown
        assert_eq!(result, "unknown");
    }

    #[test]
    fn test_classify_allele_four_all_same() {
        let fps_classified = vec![
            "qAIntactPolyA".to_string(),
            "qAIntactPolyA".to_string(),
            "qAIntactPolyA".to_string(),
            "qAIntactPolyA".to_string(),
        ];
        let result = classify_allele(&fps_classified);
        // 4 elements, all same: 4 >= 4 * 0.8 = 3.2, so should return qAIntactPolyA
        assert_eq!(result, "qAIntactPolyA");
    }

    #[test]
    fn test_classify_allele_four_first_different_rest_same() {
        let fps_classified = vec![
            "qB".to_string(),
            "qAIntactPolyA".to_string(),
            "qAIntactPolyA".to_string(),
            "qAIntactPolyA".to_string(),
        ];
        let result = classify_allele(&fps_classified);
        // 4 elements, first different: 3 < 4 * 0.8 = 3.2, so first condition fails
        // Skip first, check rest: 3 elements, all qAIntactPolyA: 3 == 4 - 1, so should return qAIntactPolyA
        assert_eq!(result, "qAIntactPolyA");
    }

    #[test]
    fn test_classify_allele_four_mixed() {
        let fps_classified = vec![
            "qAIntactPolyA".to_string(),
            "qB".to_string(),
            "qAIntactPolyA".to_string(),
            "qB".to_string(),
        ];
        let result = classify_allele(&fps_classified);
        // 4 elements, 2 of each: 2 < 4 * 0.8 = 3.2, so first condition fails
        // Skip first, check rest: 3 elements, 2 qB, 1 qAIntactPolyA: 2 != 4 - 1, so should return unknown
        assert_eq!(result, "unknown");
    }

    #[test]
    fn test_classify_allele_five_all_same() {
        let fps_classified = vec![
            "qAIntactPolyA".to_string(),
            "qAIntactPolyA".to_string(),
            "qAIntactPolyA".to_string(),
            "qAIntactPolyA".to_string(),
            "qAIntactPolyA".to_string(),
        ];
        let result = classify_allele(&fps_classified);
        // 5 elements, all same: 5 >= 5 * 0.8 = 4.0, so should return qAIntactPolyA
        assert_eq!(result, "qAIntactPolyA");
    }

    #[test]
    fn test_classify_allele_five_first_different_rest_same() {
        let fps_classified = vec![
            "qB".to_string(),
            "qAIntactPolyA".to_string(),
            "qAIntactPolyA".to_string(),
            "qAIntactPolyA".to_string(),
            "qAIntactPolyA".to_string(),
        ];
        let result = classify_allele(&fps_classified);
        // 5 elements, first different: 4 < 5 * 0.8 = 4.0, so first condition fails
        // Skip first, check rest: 4 elements, all qAIntactPolyA: 4 == 5 - 1, so should return qAIntactPolyA
        assert_eq!(result, "qAIntactPolyA");
    }

    #[test]
    fn test_classify_allele_five_mixed() {
        let fps_classified = vec![
            "qAIntactPolyA".to_string(),
            "qB".to_string(),
            "qAIntactPolyA".to_string(),
            "qB".to_string(),
            "qAIntactPolyA".to_string(),
        ];
        let result = classify_allele(&fps_classified);
        // 5 elements, 3 qAIntactPolyA, 2 qB: 3 < 5 * 0.8 = 4.0, so first condition fails
        // Skip first, check rest: 4 elements, 2 qAIntactPolyA, 2 qB: 2 != 5 - 1, so should return unknown
        assert_eq!(result, "unknown");
    }

    #[test]
    fn test_classify_allele_six_five_same() {
        let fps_classified = vec![
            "qB".to_string(),
            "qAIntactPolyA".to_string(),
            "qAIntactPolyA".to_string(),
            "qAIntactPolyA".to_string(),
            "qAIntactPolyA".to_string(),
            "qAIntactPolyA".to_string(),
        ];
        let result = classify_allele(&fps_classified);
        // 6 elements, first qB, rest qAIntactPolyA: 5 < 6 * 0.8 = 4.8, so first condition fails
        // Skip first, check rest: 5 elements, all qAIntactPolyA: 5 == 6 - 1, so should return qAIntactPolyA
        assert_eq!(result, "qAIntactPolyA");
    }

    #[test]
    fn test_classify_allele_qadisruptedpolya() {
        let fps_classified = vec![
            "qADisruptedPolyA".to_string(),
            "qADisruptedPolyA".to_string(),
            "qADisruptedPolyA".to_string(),
        ];
        let result = classify_allele(&fps_classified);
        assert_eq!(result, "qADisruptedPolyA");
    }

    #[test]
    fn test_classify_allele_qb() {
        let fps_classified = vec!["qB".to_string(), "qB".to_string()];
        let result = classify_allele(&fps_classified);
        assert_eq!(result, "qB");
    }

    #[test]
    fn test_classify_allele_type4() {
        let fps_classified = vec![
            "type4".to_string(),
            "type4".to_string(),
            "type4".to_string(),
        ];
        let result = classify_allele(&fps_classified);
        assert_eq!(result, "type4");
    }

    #[test]
    fn test_classify_allele_four_first_different_rest_mixed() {
        let fps_classified = vec![
            "qB".to_string(),
            "qAIntactPolyA".to_string(),
            "qAIntactPolyA".to_string(),
            "qB".to_string(),
        ];
        let result = classify_allele(&fps_classified);
        // 4 elements, first qB, rest: 2 qAIntactPolyA, 1 qB
        // First condition: 2 < 4 * 0.8 = 3.2, fails
        // Skip first, check rest: 3 elements, 2 qAIntactPolyA, 1 qB: 2 != 4 - 1, so should return unknown
        assert_eq!(result, "unknown");
    }

    #[test]
    fn test_classify_allele_ten_eight_same() {
        let fps_classified = vec![
            "qAIntactPolyA".to_string(),
            "qAIntactPolyA".to_string(),
            "qAIntactPolyA".to_string(),
            "qAIntactPolyA".to_string(),
            "qAIntactPolyA".to_string(),
            "qAIntactPolyA".to_string(),
            "qAIntactPolyA".to_string(),
            "qAIntactPolyA".to_string(),
            "qB".to_string(),
            "qB".to_string(),
        ];
        let result = classify_allele(&fps_classified);
        // 10 elements, 8 qAIntactPolyA: 8 >= 10 * 0.8 = 8.0, so should return qAIntactPolyA
        assert_eq!(result, "qAIntactPolyA");
    }

    #[test]
    fn test_is_cis_dup_by_read_start_offset_three_supporting_reads_two_delayed() {
        let mut read_edges = BTreeMap::new();
        read_edges.insert("read1".to_string(), vec![7, 8, 9]);
        read_edges.insert("read2".to_string(), vec![7, 8]);
        read_edges.insert("read3".to_string(), vec![7, 8]);

        let mut read_positions = BTreeMap::new();
        read_positions.insert("read1".to_string(), vec![650, 800, 950]);
        read_positions.insert("read2".to_string(), vec![620, 770]);
        read_positions.insert("read3".to_string(), vec![700, 850]);

        let fp_info = FingerprintInfo {
            read_edges,
            grouped_reads: BTreeMap::new(),
            fp_count: BTreeMap::new(),
            good_name_to_seq: BTreeMap::new(),
            read_positions,
            read_bases: BTreeMap::new(),
            fp_to_tid: BTreeMap::new(),
            variants_by_position: BTreeMap::new(),
        };

        assert!(is_cis_dup_by_read_start_offset("7-8-9-10", &fp_info).unwrap());
    }

    #[test]
    fn test_is_cis_dup_by_read_start_offset_requires_two_delayed_when_three_support() {
        let mut read_edges = BTreeMap::new();
        read_edges.insert("read1".to_string(), vec![7, 8, 9]);
        read_edges.insert("read2".to_string(), vec![7, 8]);
        read_edges.insert("read3".to_string(), vec![7, 8]);

        let mut read_positions = BTreeMap::new();
        read_positions.insert("read1".to_string(), vec![300, 450, 600]);
        read_positions.insert("read2".to_string(), vec![480, 630]);
        read_positions.insert("read3".to_string(), vec![700, 850]);

        let fp_info = FingerprintInfo {
            read_edges,
            grouped_reads: BTreeMap::new(),
            fp_count: BTreeMap::new(),
            good_name_to_seq: BTreeMap::new(),
            read_positions,
            read_bases: BTreeMap::new(),
            fp_to_tid: BTreeMap::new(),
            variants_by_position: BTreeMap::new(),
        };

        assert!(!is_cis_dup_by_read_start_offset("7-8-9-10", &fp_info).unwrap());
    }

    #[test]
    fn test_is_cis_dup_by_read_start_offset_requires_at_least_three_supporting_reads() {
        let mut read_edges = BTreeMap::new();
        read_edges.insert("read1".to_string(), vec![7, 8, 9]);
        read_edges.insert("read2".to_string(), vec![7, 8]);

        let mut read_positions = BTreeMap::new();
        read_positions.insert("read1".to_string(), vec![650, 800, 950]);
        read_positions.insert("read2".to_string(), vec![700, 850]);

        let fp_info = FingerprintInfo {
            read_edges,
            grouped_reads: BTreeMap::new(),
            fp_count: BTreeMap::new(),
            good_name_to_seq: BTreeMap::new(),
            read_positions,
            read_bases: BTreeMap::new(),
            fp_to_tid: BTreeMap::new(),
            variants_by_position: BTreeMap::new(),
        };

        assert!(!is_cis_dup_by_read_start_offset("7-8-9-10", &fp_info).unwrap());
    }

    #[test]
    fn test_is_cis_dup_by_read_start_offset_requires_matching_first_unit() {
        let mut read_edges = BTreeMap::new();
        read_edges.insert("read1".to_string(), vec![0, 7, 8, 9]);
        read_edges.insert("read2".to_string(), vec![7, 8, 4]);

        let mut read_positions = BTreeMap::new();
        read_positions.insert("read1".to_string(), vec![700, 820, 940, 1060]);
        read_positions.insert("read2".to_string(), vec![300, 450, 600]);

        let fp_info = FingerprintInfo {
            read_edges,
            grouped_reads: BTreeMap::new(),
            fp_count: BTreeMap::new(),
            good_name_to_seq: BTreeMap::new(),
            read_positions,
            read_bases: BTreeMap::new(),
            fp_to_tid: BTreeMap::new(),
            variants_by_position: BTreeMap::new(),
        };

        assert!(!is_cis_dup_by_read_start_offset("7-8-9-10", &fp_info).unwrap());
    }

    #[test]
    fn test_is_cis_dup_by_read_start_offset_requires_allele_first_unit_to_be_read_first_unit() {
        let mut read_edges = BTreeMap::new();
        read_edges.insert("read1".to_string(), vec![1, 7, 8, 9]);
        read_edges.insert("read2".to_string(), vec![7, 8, 4]);

        let mut read_positions = BTreeMap::new();
        read_positions.insert("read1".to_string(), vec![100, 820, 940, 1060]);
        read_positions.insert("read2".to_string(), vec![300, 450, 600]);

        let fp_info = FingerprintInfo {
            read_edges,
            grouped_reads: BTreeMap::new(),
            fp_count: BTreeMap::new(),
            good_name_to_seq: BTreeMap::new(),
            read_positions,
            read_bases: BTreeMap::new(),
            fp_to_tid: BTreeMap::new(),
            variants_by_position: BTreeMap::new(),
        };

        assert!(!is_cis_dup_by_read_start_offset("7-8-9-10", &fp_info).unwrap());
    }
}
