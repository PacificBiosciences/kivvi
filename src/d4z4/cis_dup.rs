use crate::assembly::assembler::{
    build_graph, match_reads_and_haplotypes, AssemblyResult, FpGraph,
};
use crate::assembly::assembler_utils::{compare_two_haps_same_length, find_overlapping_alleles};
use crate::caller::vec_to_string;
use crate::d4z4::join_partial_alleles::is_cis_dup_by_read_start_offset;
use crate::repeat_unit::fingerprint::FingerprintInfo;
use crate::util::{invalid_data_error, missing_data_error, DError};
use itertools::Itertools;
use log::{debug, trace};
use paraphase::io::json::GeneCall;
use std::cmp;
use std::collections::{BTreeMap, HashSet};

fn is_cis_dup_hap(hap: &[i32], fp_info: &FingerprintInfo) -> Result<bool, DError> {
    let Some(&first_node) = hap.first() else {
        return Ok(false);
    };
    if first_node < 0 && first_node > -10 {
        return Ok(false);
    }
    if first_node <= -10 {
        return Ok(true);
    }
    if let Some(first_node_seq) = fp_info.good_name_to_seq.get(&first_node) {
        if first_node_seq[0] == b'S' {
            return Ok(true);
        }
    }
    is_cis_dup_nodes(hap, fp_info)
}

pub(crate) fn is_cis_dup_nodes(nodes: &[i32], fp_info: &FingerprintInfo) -> Result<bool, DError> {
    if nodes.len() < 2 {
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
        if read_nodes[start_idx] == nodes[0] {
            let overlap_len = cmp::min(read_nodes.len() - start_idx, nodes.len());
            if overlap_len < 2 {
                continue;
            }

            let nodes_in_read = &read_nodes[start_idx..(start_idx + overlap_len)];
            let nodes_in_hap = &nodes[..overlap_len];
            let mut match_count = 0;
            let mut has_mismatch = false;
            for (read_node, hap_node) in nodes_in_read.iter().zip(nodes_in_hap.iter()) {
                if *read_node == 0 {
                    continue;
                }
                if read_node == hap_node {
                    match_count += 1;
                } else {
                    has_mismatch = true;
                    break;
                }
            }

            if !has_mismatch && match_count > 1 {
                supporting_reads += 1;
                if read_positions[start_idx] > 300 {
                    delayed_start_reads += 1;
                }
            }
        }
    }
    let delayed_start_threshold = (supporting_reads as f64 * 0.8).floor() as i32;
    Ok(supporting_reads >= 3
        && delayed_start_reads >= (supporting_reads - 1).min(delayed_start_threshold))
}

fn remove_redundant_haplotypes(
    haps_to_check: &[Vec<i32>],
    num_turns: usize,
    sensitive: bool,
) -> Result<BTreeMap<Vec<i32>, Vec<i32>>, DError> {
    let mut haps_to_remove = BTreeMap::<Vec<i32>, Vec<i32>>::new();
    for _turn_index in 0..num_turns {
        let haps_to_check = haps_to_check
            .iter()
            .filter(|hap| !haps_to_remove.contains_key(*hap))
            .cloned()
            .collect::<Vec<_>>();
        let (_overlapping_haps, overlapping_haps_match) =
            find_overlapping_alleles(haps_to_check.clone(), Some(2))?;
        let mut removed_one_redundant = false;
        for (hap, hap_match_info) in &overlapping_haps_match {
            let hap_size = hap.len();
            for (matching_hap, overlap_len) in hap_match_info {
                let matching_hap_size = matching_hap.len();
                let mut is_overlap = false;
                if *overlap_len >= 5 && (*overlap_len - 1) >= (hap_size - 1) / 2 {
                    is_overlap = true;
                } else if *overlap_len >= 2 {
                    let hap_unique_units = hap
                        .iter()
                        .filter(|x| !matching_hap.contains(x))
                        .collect::<HashSet<_>>();
                    let matching_hap_unique_units = matching_hap
                        .iter()
                        .filter(|x| !hap.contains(x))
                        .collect::<HashSet<_>>();
                    if *overlap_len >= 4
                        && (hap_unique_units.len() <= (hap_size as f64 * 0.2).floor() as usize
                            || matching_hap_unique_units.len()
                                <= (matching_hap_size as f64 * 0.2).floor() as usize)
                    {
                        is_overlap = true;
                    }
                    if sensitive
                        && *overlap_len == 3
                        && (hap_size == *overlap_len + 1 || matching_hap_size == *overlap_len + 1)
                    {
                        is_overlap = true;
                    }
                }
                if is_overlap {
                    if hap_size > matching_hap_size {
                        continue;
                    }
                    if hap_size == matching_hap_size {
                        let mut haps = vec![hap.clone(), matching_hap.clone()];
                        haps.sort();
                        let hap1 = haps[0].clone();
                        let hap2 = haps[1].clone();
                        if !haps_to_remove.contains_key(&hap1) {
                            haps_to_remove.insert(hap1, hap2);
                            removed_one_redundant = true;
                            break;
                        }
                    }
                    if hap_size < matching_hap_size && !haps_to_remove.contains_key(hap) {
                        haps_to_remove.insert(hap.clone(), matching_hap.clone());
                        removed_one_redundant = true;
                        break;
                    }
                }
            }
            if removed_one_redundant {
                break;
            }
        }
    }
    Ok(haps_to_remove)
}

pub fn process_alleles(
    assembly_result: &AssemblyResult,
    fp_info: &FingerprintInfo,
) -> Result<(Vec<Vec<i32>>, Vec<Vec<i32>>, Vec<Vec<i32>>), DError> {
    debug!("Identify proximal and distal alleles and remove redundant ones...");
    let mut all_starting_haps = HashSet::new();
    let mut all_ending_haps = HashSet::new();
    for hap in &assembly_result.complete {
        if let Some(hap_first) = hap.first() {
            if *hap_first < 0 && *hap_first > -10 {
                all_starting_haps.insert(hap.clone());
            }
        }
        if let Some(hap_end) = hap.last() {
            if *hap_end <= -10 {
                all_ending_haps.insert(hap.clone());
            }
        }
    }
    for hap in &assembly_result.incomplete {
        if let Some(hap_first) = hap.first() {
            if *hap_first < 0 && *hap_first > -10 {
                all_starting_haps.insert(hap.clone());
            }
        }
        if let Some(hap_end) = hap.last() {
            if *hap_end <= -10 {
                all_ending_haps.insert(hap.clone());
            }
        }
    }
    debug!("all_starting_haps before removing redundant {all_starting_haps:?}");
    debug!("all_ending_haps before removing redundant {all_ending_haps:?}");
    let mut kept_starting_haps: Vec<Vec<i32>> = all_starting_haps.into_iter().collect();
    let mut kept_ending_haps: Vec<Vec<i32>> = all_ending_haps.into_iter().collect();
    let mut kept_complete = assembly_result.complete.clone();
    let mut kept_complete_set = kept_complete.iter().cloned().collect::<HashSet<_>>();

    let distal_no_cis_dup = kept_ending_haps
        .iter()
        .filter(|hap| !is_cis_dup_hap(hap, fp_info).unwrap_or(false))
        .cloned()
        .collect::<Vec<Vec<i32>>>();
    let size1_allele = distal_no_cis_dup
        .iter()
        .filter(|hap| hap.len() == 2)
        .cloned()
        .collect::<Vec<Vec<i32>>>();
    if size1_allele.len() == 1 && distal_no_cis_dup.len() >= 5 {
        if let Some(size1_allele) = size1_allele.first() {
            debug!("removing size1_allele {size1_allele:?}");
            kept_ending_haps.retain(|hap| hap != size1_allele);
            kept_complete_set.remove(size1_allele);
        }
    }
    if kept_starting_haps.len() >= 5 {
        let num_turns = kept_starting_haps.len() - 4;
        debug!("removing redundant proximal alleles, num_turns {num_turns}");
        let proximal_to_remove =
            remove_redundant_haplotypes(&kept_starting_haps, num_turns, false)?;
        debug!("proximal_to_remove {proximal_to_remove:?}");
        if proximal_to_remove.len() <= num_turns {
            for (proximal_to_remove_allele, redundant_allele) in &proximal_to_remove {
                if !(kept_ending_haps.len() == 4
                    && kept_ending_haps.contains(proximal_to_remove_allele))
                    && kept_starting_haps.contains(proximal_to_remove_allele)
                {
                    kept_starting_haps.retain(|hap| hap != proximal_to_remove_allele);
                    kept_complete_set.remove(proximal_to_remove_allele);
                } else if !(kept_ending_haps.len() == 4
                    && kept_ending_haps.contains(redundant_allele))
                {
                    kept_starting_haps.retain(|hap| hap != redundant_allele);
                    kept_complete_set.remove(redundant_allele);
                }
            }
        }
    }
    if kept_starting_haps.len() >= 5 {
        let num_turns = kept_starting_haps.len() - 4;
        debug!("removing redundant proximal alleles, num_turns {num_turns}");
        let proximal_to_remove = remove_redundant_haplotypes(&kept_starting_haps, num_turns, true)?;
        debug!("proximal_to_remove {proximal_to_remove:?}");
        if proximal_to_remove.len() <= num_turns {
            for (proximal_to_remove_allele, redundant_allele) in &proximal_to_remove {
                if !(kept_ending_haps.len() == 4
                    && kept_ending_haps.contains(proximal_to_remove_allele))
                    && kept_starting_haps.contains(proximal_to_remove_allele)
                {
                    kept_starting_haps.retain(|hap| hap != proximal_to_remove_allele);
                    kept_complete_set.remove(proximal_to_remove_allele);
                } else if !(kept_ending_haps.len() == 4
                    && kept_ending_haps.contains(redundant_allele))
                {
                    kept_starting_haps.retain(|hap| hap != redundant_allele);
                    kept_complete_set.remove(redundant_allele);
                }
            }
        }
    }
    let distal_no_cis_dup = kept_ending_haps
        .iter()
        .filter(|hap| !is_cis_dup_hap(hap, fp_info).unwrap_or(false))
        .cloned()
        .collect::<Vec<Vec<i32>>>();
    if distal_no_cis_dup.len() >= 5 {
        let num_turns = distal_no_cis_dup.len() - 4;
        debug!("removing redundant distal alleles, num_turns {num_turns}");
        let distal_to_remove = remove_redundant_haplotypes(&distal_no_cis_dup, num_turns, false)?;
        debug!("distal_to_remove {distal_to_remove:?}");
        if distal_to_remove.len() <= num_turns {
            for (distal_to_remove_allele, redundant_allele) in &distal_to_remove {
                if !(kept_starting_haps.len() == 4
                    && kept_starting_haps.contains(distal_to_remove_allele))
                    && kept_ending_haps.contains(distal_to_remove_allele)
                {
                    debug!("removing distal_to_remove_allele {distal_to_remove_allele:?}");
                    kept_ending_haps.retain(|hap| hap != distal_to_remove_allele);
                    kept_complete_set.remove(distal_to_remove_allele);
                } else if !(kept_starting_haps.len() == 4
                    && kept_starting_haps.contains(redundant_allele))
                {
                    debug!("removing redundant_allele {redundant_allele:?}");
                    kept_ending_haps.retain(|hap| hap != redundant_allele);
                    kept_complete_set.remove(redundant_allele);
                }
            }
        }
    }
    let distal_no_cis_dup = kept_ending_haps
        .iter()
        .filter(|hap| !is_cis_dup_hap(hap, fp_info).unwrap_or(false))
        .cloned()
        .collect::<Vec<Vec<i32>>>();
    if distal_no_cis_dup.len() >= 5 {
        let num_turns = distal_no_cis_dup.len() - 4;
        debug!("removing redundant distal alleles, num_turns {num_turns}");
        let distal_to_remove = remove_redundant_haplotypes(&distal_no_cis_dup, num_turns, true)?;
        debug!("distal_to_remove {distal_to_remove:?}");
        if distal_to_remove.len() <= num_turns {
            for (distal_to_remove_allele, redundant_allele) in &distal_to_remove {
                if !(kept_starting_haps.len() == 4
                    && kept_starting_haps.contains(distal_to_remove_allele))
                    && kept_ending_haps.contains(distal_to_remove_allele)
                {
                    debug!("removing distal_to_remove_allele {distal_to_remove_allele:?}");
                    kept_ending_haps.retain(|hap| hap != distal_to_remove_allele);
                    kept_complete_set.remove(distal_to_remove_allele);
                } else if !(kept_starting_haps.len() == 4
                    && kept_starting_haps.contains(redundant_allele))
                {
                    debug!("removing redundant_allele {redundant_allele:?}");
                    kept_ending_haps.retain(|hap| hap != redundant_allele);
                    kept_complete_set.remove(redundant_allele);
                }
            }
        }
    }

    kept_complete.retain(|hap| kept_complete_set.contains(hap));

    debug!("all_starting_haps after removing redundant {kept_starting_haps:?}");
    debug!("all_ending_haps after removing redundant {kept_ending_haps:?}");
    debug!("complete_haps after removing redundant {kept_complete:?}");

    Ok((kept_starting_haps, kept_ending_haps, kept_complete))
}

pub fn find_qal_alleles(all_haps: &[Vec<i32>], qal_units: Vec<i32>) -> HashSet<Vec<i32>> {
    all_haps
        .iter()
        .filter(|hap| hap.last() == Some(&-10) && hap.len() >= 2)
        .filter_map(|hap| {
            let second_to_last_unit = hap[hap.len() - 2];
            if qal_units.contains(&second_to_last_unit) {
                Some(hap.clone())
            } else {
                None
            }
        })
        .collect::<HashSet<_>>()
}

/// Mutable bookkeeping used while assembling cis-dup paths from multiple evidence sources.
#[derive(Default)]
struct CisDupAssemblyState {
    haps_to_node_names: BTreeMap<Vec<i32>, i32>,
    read_edges_for_haps: BTreeMap<String, Vec<i32>>,
    next_node_name: i32,
}

impl CisDupAssemblyState {
    fn new() -> Self {
        Self {
            next_node_name: 1,
            ..Self::default()
        }
    }

    fn ensure_node_name(&mut self, hap: &[i32]) -> i32 {
        if let Some(node_name) = self.haps_to_node_names.get(hap) {
            *node_name
        } else {
            let node_name = self.next_node_name;
            self.haps_to_node_names.insert(hap.to_vec(), node_name);
            self.next_node_name += 1;
            node_name
        }
    }

    fn set_read_path(&mut self, read_name: String, haps: &[Vec<i32>]) {
        let node_path = haps
            .iter()
            .map(|hap| self.ensure_node_name(hap))
            .collect::<Vec<_>>();
        self.read_edges_for_haps
            .entry(read_name)
            .or_insert(node_path);
    }
}

/// Identify haplotypes that independently satisfy the cis-dup read-start-offset heuristic.
fn identify_cis_dup_alleles(
    all_haps: &[Vec<i32>],
    fp_info: &FingerprintInfo,
) -> Result<HashSet<Vec<i32>>, DError> {
    all_haps
        .iter()
        .filter_map(|hap| {
            let allele_name = vec_to_string(&vec![hap.clone()], "-")
                .into_iter()
                .next()
                .unwrap_or_default();
            match is_cis_dup_by_read_start_offset(&allele_name, fp_info) {
                Ok(true) => Some(Ok(hap.clone())),
                Ok(false) => None,
                Err(err) => Some(Err(err)),
            }
        })
        .collect()
}

/// Partition repeat haplotypes linked to a downstream paraphase haplotype into upstream and downstream sets.
fn linked_repeat_haps_for_downstream_haplotype(
    downstream_hap_reads: &[String],
    all_haps_support: &BTreeMap<String, Vec<Vec<i32>>>,
    fp_info: &FingerprintInfo,
) -> Result<(HashSet<Vec<i32>>, HashSet<Vec<i32>>), DError> {
    let mut upstream = Vec::new();
    let mut downstream = Vec::new();
    for downstream_hap_read in downstream_hap_reads {
        let fields = downstream_hap_read.split("_sup_").collect::<Vec<_>>();
        let downstream_hap_read_name = fields[0];
        let aln_pos = fields
            .get(1)
            .ok_or_else(|| {
                invalid_data_error(format!(
                    "Downstream haplotype read annotation is missing a '_sup_' position suffix: '{downstream_hap_read}'"
                ))
            })?
            .split('_')
            .next()
            .ok_or_else(|| {
                invalid_data_error(format!(
                    "Downstream haplotype read annotation has an empty '_sup_' position field: '{downstream_hap_read}'"
                ))
            })?
            .parse::<i32>()?;
        trace!("downstream_hap_read {downstream_hap_read} downstream_hap_read_name {downstream_hap_read_name} pos {aln_pos}");
        if all_haps_support.contains_key(downstream_hap_read_name) {
            let this_read_repeat_support = all_haps_support
                .get(downstream_hap_read_name)
                .ok_or_else(|| {
                    missing_data_error(
                        "repeat-support haplotypes for read",
                        downstream_hap_read_name.to_string(),
                    )
                })?;
            let this_read_repeat_edges = fp_info
                .read_edges
                .get(downstream_hap_read_name)
                .ok_or_else(|| {
                    missing_data_error(
                        "repeat-edge path for read",
                        downstream_hap_read_name.to_string(),
                    )
                })?;
            let this_read_repeat_positions = fp_info
                .read_positions
                .get(downstream_hap_read_name)
                .ok_or_else(|| {
                missing_data_error(
                    "repeat-position path for read",
                    downstream_hap_read_name.to_string(),
                )
            })?;
            trace!("this_read_repeat_edges {this_read_repeat_edges:?}");
            trace!("this_read_repeat_positions {this_read_repeat_positions:?}");
            for repeat_hap in this_read_repeat_support {
                let node_index =
                    match_read_allele_first_node_index(this_read_repeat_edges, repeat_hap);
                trace!("matching repeat_hap {repeat_hap:?} node_index {node_index:?}");
                if let Some(node_index_value) = node_index {
                    let matching_position_on_read = this_read_repeat_positions[node_index_value.0];
                    trace!("matching_position_on_read {matching_position_on_read}");
                    if matching_position_on_read < aln_pos {
                        upstream.push(repeat_hap.clone());
                    } else if matching_position_on_read > aln_pos {
                        downstream.push(repeat_hap.clone());
                    }
                }
            }
        }
    }
    Ok((
        upstream.into_iter().collect::<HashSet<_>>(),
        downstream.into_iter().collect::<HashSet<_>>(),
    ))
}

/// Add a simple linear haplotype path into the temporary cis-dup assembly state.
fn add_linear_path(
    allele_links: &mut BTreeMap<Vec<i32>, Vec<Vec<i32>>>,
    state: &mut CisDupAssemblyState,
    read_name: &str,
    haps: &[Vec<i32>],
) {
    for pair in haps.windows(2) {
        trace!("adding non-read links {:?} to {:?}", pair[0], pair[1]);
        allele_links
            .entry(pair[0].clone())
            .or_default()
            .push(pair[1].clone());
    }
    state.set_read_path(read_name.to_string(), haps);
}

/// Convert one downstream paraphase support pattern into haplotype links when it matches a supported cis-dup layout.
fn add_paraphase_supported_links(
    allele_links: &mut BTreeMap<Vec<i32>, Vec<Vec<i32>>>,
    state: &mut CisDupAssemblyState,
    downstream_hap: &str,
    upstream_haps: &HashSet<Vec<i32>>,
    downstream_haps: &HashSet<Vec<i32>>,
    cis_dup_alleles: &HashSet<Vec<i32>>,
    qal_alleles: &HashSet<Vec<i32>>,
) {
    debug!(
        "paraphase_hap {downstream_hap} linking repeat haps up {upstream_haps:?} down {downstream_haps:?}"
    );
    // one upstream array linked to one downstream array
    // array1 -> paraphase hap -> array2
    if upstream_haps.len() == 1 && downstream_haps.len() == 1 {
        if let (Some(upstream), Some(downstream)) =
            (upstream_haps.iter().next(), downstream_haps.iter().next())
        {
            let path = vec![upstream.to_vec(), downstream.to_vec()];
            add_linear_path(allele_links, state, downstream_hap, &path);
        }
        return;
    }
    // three arrays in a row, linked to the same paraphase haplotype
    // array1 -> paraphase hap -> array2 -> paraphase hap -> array3
    if upstream_haps.len() == 2 && downstream_haps.len() == 2 {
        let overlap = upstream_haps
            .iter()
            .filter(|hap| downstream_haps.contains(*hap))
            .cloned()
            .collect::<Vec<_>>();
        if overlap.len() == 1 {
            let middle = overlap[0].clone();
            let first = upstream_haps
                .iter()
                .filter(|hap| **hap != middle)
                .cloned()
                .next();
            let last = downstream_haps
                .iter()
                .filter(|hap| **hap != middle)
                .cloned()
                .next();
            if let (Some(first), Some(last)) = (first, last) {
                add_linear_path(allele_links, state, downstream_hap, &[first, middle, last]);
            }
        }
        return;
    }
    // one upstream array linked to one downstream array and they are linked to the same paraphase haplotype
    // array1 -> paraphase hap -> array2 -> paraphase hap
    if upstream_haps.len() == 2 && downstream_haps.len() == 1 {
        let overlap = upstream_haps
            .iter()
            .filter(|hap| downstream_haps.contains(*hap))
            .cloned()
            .collect::<Vec<_>>();
        if overlap.len() == 1 {
            let downstream = overlap[0].clone();
            let upstream = upstream_haps
                .iter()
                .filter(|hap| **hap != downstream)
                .cloned()
                .next();
            if let Some(upstream) = upstream {
                add_linear_path(allele_links, state, downstream_hap, &[upstream, downstream]);
            }
        }
        return;
    }
    if downstream_haps.len() > 2 {
        return;
    }
    // now we want to have a special logic around qAL arrays
    // our assumption is that qAL arrays are on the same allele
    // our scenario is when the paraphase haplotype is linked to many d4z4 arrays
    let downstream_not_cis_dup = downstream_haps
        .iter()
        .filter(|hap| !cis_dup_alleles.contains(*hap))
        .cloned()
        .collect::<HashSet<_>>();
    let upstream_not_cis_dup_and_qal = upstream_haps
        .iter()
        .filter(|hap| !cis_dup_alleles.contains(*hap))
        .filter(|hap| qal_alleles.contains(*hap))
        .cloned()
        .collect::<HashSet<_>>();
    let downstream_not_qal = downstream_haps
        .iter()
        .filter(|hap| !qal_alleles.contains(*hap))
        .cloned()
        .collect::<HashSet<_>>();

    debug!(
        "paraphase_hap_linking_repeat_haps_downstream_set_is_not_cis_dup {:?}",
        downstream_not_cis_dup
    );
    debug!(
        "paraphase_hap_linking_repeat_haps_downstream_set_is_not_qal {:?}",
        downstream_not_qal
    );
    debug!(
        "paraphase_hap_linking_repeat_haps_upstream_set_is_not_cis_dup_and_is_qal {:?}",
        upstream_not_cis_dup_and_qal
    );

    // we want a scenario where among the upstream arrays, there is only one that is qal but not cis-dup
    // and all the downstream arrays are cis-dup and qal
    if !downstream_not_cis_dup.is_empty()
        || !downstream_not_qal.is_empty()
        || upstream_not_cis_dup_and_qal.len() != 1
    {
        return;
    }

    // the logic is: there is one downstream array and it's cis-dup and qal.
    // There is more that one upstream array, but only one of them is qal.
    // Then we think we see a link between the qal upstream array and qal downstream array
    if downstream_haps.len() == 1 {
        if let (Some(upstream), Some(downstream)) = (
            upstream_not_cis_dup_and_qal.iter().next(),
            downstream_haps.iter().next(),
        ) {
            let path = vec![upstream.to_vec(), downstream.to_vec()];
            add_linear_path(allele_links, state, downstream_hap, &path);
        }
    } else if downstream_haps.len() == 2 {
        // if there are two downstream arrays, we sort them by size
        // note that here the order is not supported by reads. we should flag that - TODO
        let downstream_sorted = downstream_haps
            .iter()
            .cloned()
            .sorted_by_key(|hap| std::cmp::Reverse(hap.len()))
            .collect::<Vec<_>>();
        if let (Some(upstream), Some(first), Some(last)) = (
            upstream_not_cis_dup_and_qal.iter().next(),
            downstream_sorted.first(),
            downstream_sorted.last(),
        ) {
            let path = vec![upstream.to_vec(), first.to_vec(), last.to_vec()];
            add_linear_path(allele_links, state, downstream_hap, &path);
        }
    }
}

/// Split a read path into segments separated by `-10` terminal markers.
fn split_read_segments(read_nodes: &[i32]) -> Vec<Vec<i32>> {
    let mut segments = Vec::new();
    let mut starting_index = 0;
    for (index, node) in read_nodes.iter().enumerate() {
        if *node == -10 {
            segments.push(read_nodes[starting_index..(index + 1)].to_vec());
            starting_index = index + 1;
        }
    }
    let trailing_segment = &read_nodes[starting_index..];
    if !trailing_segment.is_empty() {
        segments.push(trailing_segment.to_vec());
    }
    segments
}

/// Resolve the single qualifying allele match for one read segment, if the segment is unambiguous.
fn matching_allele_for_segment(
    segment: &[i32],
    segment_index: usize,
    all_haps: &[Vec<i32>],
    special_incomplete_haps: &[Vec<i32>],
    first_read_position: i32,
) -> Vec<i32> {
    let mut dummy_read = BTreeMap::new();
    dummy_read.insert(String::from("read"), segment.to_vec());
    let segment_match = if segment_index == 0 {
        match_reads_and_haplotypes(&dummy_read, all_haps, None, false).by_read
    } else {
        match_reads_and_haplotypes(&dummy_read, special_incomplete_haps, None, false).by_read
    };
    let Some(dummy_read_matches) = segment_match.get("read") else {
        return vec![];
    };

    let qualifying_matches = dummy_read_matches
        .iter()
        .filter(|candidate| {
            if segment_index == 0 {
                if !special_incomplete_haps.contains(candidate) {
                    return true;
                }
                let node_index = match_read_allele_first_node_index(&segment.to_vec(), candidate);
                node_index.is_some_and(|node_index| {
                    node_index.0 == 0 && (first_read_position < 2000 || node_index.1 == 0)
                })
            } else {
                let node_index = match_read_allele_first_node_index(&segment.to_vec(), candidate);
                node_index.is_some_and(|node_index| node_index.0 == 0 && node_index.1 == 0)
            }
        })
        .cloned()
        .collect::<Vec<_>>();

    if qualifying_matches.len() == 1 {
        qualifying_matches[0].to_vec()
    } else {
        vec![]
    }
}

/// Normalize a matched allele and its read offset for downstream match-index bookkeeping.
fn normalized_match_index(
    read_name: &str,
    segment: &[i32],
    matched_hap: &[i32],
    prev_segments_len: usize,
    extra_offset: usize,
) -> Result<(Vec<i32>, (String, i32)), DError> {
    let node_index = match_read_allele_first_node_index(&segment.to_vec(), &matched_hap.to_vec())
        .ok_or_else(|| {
        invalid_data_error(format!(
            "Could not determine the first matching node index between segment {:?} and hap {:?}",
            segment, matched_hap
        ))
    })?;
    let mut index_on_read = if node_index.1 > 0 {
        node_index.1 as i32
    } else {
        0 - node_index.0 as i32 - prev_segments_len as i32 - extra_offset as i32
    };
    if matched_hap.starts_with(&[-10]) && node_index.1 == 0 {
        index_on_read += 1;
    }
    let mut match_seg = matched_hap.to_vec();
    let match_len = match_seg.len();
    if match_len > 6 {
        match_seg = vec![0, 0, 0, 0, 0, 0];
        match_seg.copy_from_slice(&matched_hap[(match_len - 6)..]);
        index_on_read -= match_len as i32 - 6;
    }
    Ok((match_seg, (read_name.to_string(), index_on_read)))
}

/// Add cis-dup links inferred directly from multi-segment reads and record their match offsets.
fn add_read_supported_links(
    allele_links: &mut BTreeMap<Vec<i32>, Vec<Vec<i32>>>,
    state: &mut CisDupAssemblyState,
    match_index_on_read: &mut BTreeMap<Vec<i32>, HashSet<(String, i32)>>,
    all_haps: &[Vec<i32>],
    fp_info: &FingerprintInfo,
    special_incomplete_haps: &[Vec<i32>],
) -> Result<(), DError> {
    for (read, read_nodes) in &fp_info.read_edges {
        if !read_nodes.contains(&-10) {
            continue;
        }
        let end_index = read_nodes.iter().position(|x| *x == -10).ok_or_else(|| {
            missing_data_error("end (-10) node in read", format!("{read_nodes:?}"))
        })?;
        if end_index == read_nodes.len() - 1 {
            continue;
        }
        let this_read_repeat_positions = fp_info.read_positions.get(read).ok_or_else(|| {
            missing_data_error(
                "repeat-position path while adding read-supported links",
                read.to_string(),
            )
        })?;
        let this_read_first_position = *this_read_repeat_positions.first().ok_or_else(|| {
            missing_data_error(
                "first repeat position while adding read-supported links",
                read.to_string(),
            )
        })?;
        let segments = split_read_segments(read_nodes);
        // we have separated the two or more segments of a read
        // now we want to identify which d4z4 array each segment unambiguously map to
        let matching_alleles = segments
            .iter()
            .enumerate()
            .map(|(index, segment)| {
                matching_allele_for_segment(
                    segment,
                    index,
                    all_haps,
                    special_incomplete_haps,
                    this_read_first_position,
                )
            })
            .collect::<Vec<_>>();

        trace!(
            "read_nodes {read_nodes:?} segments {segments:?} matching_alleles {matching_alleles:?}"
        );
        let mut prev_segments_len = 0;
        for segment_index in 0..matching_alleles.len().saturating_sub(1) {
            let match1 = &matching_alleles[segment_index];
            let match2 = &matching_alleles[segment_index + 1];
            if match1.is_empty() || match2.is_empty() {
                prev_segments_len += segments[segment_index].len();
                continue;
            }

            trace!("adding read links {match1:?} to {match2:?}");
            allele_links
                .entry(match1.to_vec())
                .or_default()
                .push(match2.to_vec());

            let (match1_seg, match1_entry) = normalized_match_index(
                read,
                &segments[segment_index],
                match1,
                prev_segments_len,
                0,
            )?;
            match_index_on_read
                .entry(match1_seg)
                .or_default()
                .insert(match1_entry);

            let (match2_seg, match2_entry) = normalized_match_index(
                read,
                &segments[segment_index + 1],
                match2,
                prev_segments_len,
                segments[segment_index].len(),
            )?;
            match_index_on_read
                .entry(match2_seg)
                .or_default()
                .insert(match2_entry);

            state.set_read_path(read.to_string(), &[match1.clone(), match2.clone()]);
            prev_segments_len += segments[segment_index].len();
        }
    }
    Ok(())
}

/// Assemble the temporary cis-dup graph into ordered haplotype paths.
fn assemble_cis_dup_paths(
    state: CisDupAssemblyState,
    special_incomplete_haps: &[Vec<i32>],
) -> Result<Vec<Vec<Vec<i32>>>, DError> {
    debug!("haps_to_node_names {:?}", state.haps_to_node_names);
    debug!("read_edges_for_haps {:?}", state.read_edges_for_haps);
    let mut node_names_to_haps = BTreeMap::new();
    for (hap, node_name) in &state.haps_to_node_names {
        node_names_to_haps.insert(*node_name, hap.to_vec());
    }
    debug!("special_incomplete {:?}", special_incomplete_haps);
    debug!("For cis-dups, assemble haplotypes into alleles...");
    let mut fp_graph = build_graph(state.read_edges_for_haps, 2);
    let allele_phase_result = fp_graph.run_simple()?;
    debug!("allele_phase_result {:?}", allele_phase_result.incomplete);
    Ok(allele_phase_result
        .incomplete
        .into_iter()
        .map(|cis_dup| {
            let cis_dup_assembled = cis_dup
                .iter()
                .map(|x| {
                    node_names_to_haps.get(x).cloned().ok_or_else(|| {
                        missing_data_error("haplotype for temporary cis-dup node", x.to_string())
                    })
                })
                .collect::<Result<Vec<_>, _>>()?;
            debug!("graph assembled allele: {cis_dup_assembled:?}");
            Ok(cis_dup_assembled)
        })
        .collect::<Result<Vec<_>, DError>>()?)
}

pub fn find_cis_dup(
    all_haps: &Vec<Vec<i32>>,
    fp_graph: &FpGraph,
    fp_info: &FingerprintInfo,
    phasing_result: &BTreeMap<String, GeneCall>,
    qal_units: Vec<i32>,
    special_incomplete_haps: &Vec<Vec<i32>>,
) -> Result<
    (
        Vec<Vec<Vec<i32>>>,
        BTreeMap<Vec<i32>, HashSet<(String, i32)>>,
    ),
    DError,
> {
    let mut state = CisDupAssemblyState::new();
    let qal_alleles = find_qal_alleles(all_haps, qal_units);
    debug!("qal_alleles identified by long insertion before -10 {qal_alleles:?}");
    let cis_dup_alleles = identify_cis_dup_alleles(all_haps, fp_info)?;
    debug!("cis_dup_alleles identified by read start offset {cis_dup_alleles:?}");
    let all_haps_support = fp_graph
        .process_complete_haps(&all_haps, Some(1), true, false)?
        .support_by_read;
    let mut allele_links: BTreeMap<Vec<i32>, Vec<Vec<i32>>> = BTreeMap::new();
    let downstream_phasing_result = phasing_result.get(&String::from("DUX4")).ok_or_else(|| {
        missing_data_error("DUX4 phasing result while finding cis-dup alleles", "DUX4")
    })?;
    let downstream_reads = &downstream_phasing_result.unique_supporting_reads;
    for (downstream_hap, downstream_hap_reads) in downstream_reads {
        trace!("checking reads for downstream_hap {downstream_hap}");
        // identify d4z4 alleles upstream of this paraphase haplotype
        // identify d4z4 alleles downstream of this paraphase haplotype
        let (upstream_haps, downstream_haps) = linked_repeat_haps_for_downstream_haplotype(
            downstream_hap_reads,
            &all_haps_support,
            fp_info,
        )?;
        // identify links between an upstream d4z4 allele and a downstream d4z4 allele
        add_paraphase_supported_links(
            &mut allele_links,
            &mut state,
            downstream_hap,
            &upstream_haps,
            &downstream_haps,
            &cis_dup_alleles,
            &qal_alleles,
        );
    }

    // some d4z4 reads might contain links between two d4z4 arrays
    let mut match_index_on_read: BTreeMap<Vec<i32>, HashSet<(String, i32)>> = BTreeMap::new();
    add_read_supported_links(
        &mut allele_links,
        &mut state,
        &mut match_index_on_read,
        all_haps,
        fp_info,
        special_incomplete_haps,
    )?;
    let cis_dups_assembled = assemble_cis_dup_paths(state, special_incomplete_haps)?;
    Ok((cis_dups_assembled, match_index_on_read))
}

pub fn match_read_allele_first_node_index(
    hap1: &Vec<i32>,
    hap2: &Vec<i32>,
) -> Option<(usize, usize)> {
    let mut hap2_mod = hap2.clone();
    if hap2_mod.starts_with(&[-10]) {
        hap2_mod.remove(0);
    }
    let hap2 = &hap2_mod;
    let hap1_len = hap1.len();
    let hap2_len = hap2.len();
    for i in 0..hap1_len {
        if i == 0 {
            for k in 0..hap2_len {
                let offset_index = cmp::min(hap2_len - k, hap1_len);
                if offset_index > 1 && k + hap1_len <= hap2_len {
                    let test_hap1 = &hap1[..offset_index];
                    let test_hap2 = &hap2[k..(k + offset_index)];
                    let (hap_match, mismatch) = compare_two_haps_same_length(test_hap1, test_hap2);
                    if hap_match.iter().sum::<i32>() >= 2 && mismatch == 0 {
                        return Some((i, k));
                    }
                }
            }
        } else {
            let offset_index = cmp::min(hap1_len - i, hap2_len);
            if offset_index > 1 && i + hap2_len >= hap1_len {
                let test_hap1 = &hap1[i..(i + offset_index)];
                let test_hap2 = &hap2[..offset_index];
                let (hap_match, mismatch) = compare_two_haps_same_length(test_hap1, test_hap2);
                if hap_match.iter().sum::<i32>() >= 2 && mismatch == 0 {
                    return Some((i, 0));
                }
            }
        }
    }
    None
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::BTreeMap;

    #[test]
    fn test_split_read_segments() {
        let read_nodes = vec![1, 2, -10, 3, 4, -10, 5, 6];

        let segments = split_read_segments(&read_nodes);

        assert_eq!(segments, vec![vec![1, 2, -10], vec![3, 4, -10], vec![5, 6]]);
    }

    #[test]
    fn test_match_read_allele_first_node_index() {
        let read = vec![4, 5, 6];
        let hap = vec![1, 2, 3, 0, 4, 5, 6];
        let res = match_read_allele_first_node_index(&read, &hap);
        assert!(!res.is_none());
        assert_eq!(res.unwrap(), (0, 4));

        let read = vec![7, 1, 2];
        let hap = vec![1, 2, 3, 0, 4, 5, 6];
        let res = match_read_allele_first_node_index(&read, &hap);
        assert!(!res.is_none());
        assert_eq!(res.unwrap(), (1, 0));

        let read = vec![4, 5, 6, 7];
        let hap = vec![1, 2, 3, 0, 4, 5, 6];
        let res = match_read_allele_first_node_index(&read, &hap);
        assert!(res.is_none());
    }

    #[test]
    fn test_linked_repeat_haps_for_downstream_haplotype_errors_on_missing_sup_suffix() {
        let error = linked_repeat_haps_for_downstream_haplotype(
            &[String::from("read1")],
            &BTreeMap::new(),
            &FingerprintInfo {
                read_edges: BTreeMap::new(),
                grouped_reads: BTreeMap::new(),
                fp_count: BTreeMap::new(),
                good_name_to_seq: BTreeMap::new(),
                read_positions: BTreeMap::new(),
                read_bases: BTreeMap::new(),
                fp_to_tid: BTreeMap::new(),
                variants_by_position: BTreeMap::new(),
            },
        )
        .expect_err("downstream read annotations without _sup_ suffix should error");

        assert!(
            error.to_string().contains(
                "invalid data: Downstream haplotype read annotation is missing a '_sup_' position suffix: 'read1'"
            ),
            "unexpected error: {error}"
        );
    }

    #[test]
    fn test_find_cis_dup_errors_on_missing_dux4_phasing_result() {
        let fp_graph = build_graph(BTreeMap::new(), 2);
        let fp_info = FingerprintInfo {
            read_edges: BTreeMap::new(),
            grouped_reads: BTreeMap::new(),
            fp_count: BTreeMap::new(),
            good_name_to_seq: BTreeMap::new(),
            read_positions: BTreeMap::new(),
            read_bases: BTreeMap::new(),
            fp_to_tid: BTreeMap::new(),
            variants_by_position: BTreeMap::new(),
        };

        let error = find_cis_dup(
            &vec![],
            &fp_graph,
            &fp_info,
            &BTreeMap::new(),
            vec![],
            &vec![],
        )
        .expect_err("missing DUX4 phasing results should error");

        assert!(
            error
                .to_string()
                .contains("missing DUX4 phasing result while finding cis-dup alleles: DUX4"),
            "unexpected error: {error}"
        );
    }
}
