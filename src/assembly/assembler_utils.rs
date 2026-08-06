use crate::util::{missing_data_error, DError};
use log::{debug, error, trace};
use std::cmp;
use std::collections::{BTreeMap, HashSet};

/// Find alleles from a group of alleles that are overlapping with each other.
/// Require at least 4 overlapping nodes.
/// # Arguments
/// * `haps_to_assess` - alleles to assess
/// * `min_overlap_len` - minimum overlap length
/// # Returns
/// * `Vec<Vec<i32>>` - overlapping alleles
/// * `BTreeMap<Vec<i32>, Vec<(Vec<i32>, usize)>>` - overlapping alleles and their overlaps
pub fn find_overlapping_alleles(
    haps_to_assess: Vec<Vec<i32>>,
    min_overlap_len: Option<usize>,
) -> Result<(Vec<Vec<i32>>, BTreeMap<Vec<i32>, Vec<(Vec<i32>, usize)>>), DError> {
    let min_overlap_len = min_overlap_len.unwrap_or(4);
    let mut overlapping_haps = Vec::new();
    let mut overlapping_haps_match = BTreeMap::<Vec<i32>, Vec<(Vec<i32>, usize)>>::new();
    for hap1 in &haps_to_assess {
        for hap2 in &haps_to_assess {
            if *hap1 != *hap2 {
                trace!("comparing haplotypes {hap1:?} and {hap2:?}");
                let mut is_overlapping = false;
                let mut forward_overlap_len: usize;
                let mut backward_overlap_len: usize;
                let hap1_len = hap1.len();
                let hap2_len = hap2.len();
                // forward compare
                let shorter_len = cmp::min(hap1_len, hap2_len);
                let mut i_index: usize = 2;
                for i in 2..(shorter_len + 1) {
                    i_index = i;
                    if hap1[..i] != hap2[..i] {
                        break;
                    }
                }
                forward_overlap_len = i_index - 1;
                // fully contained
                if hap1[..shorter_len] == hap2[..shorter_len] {
                    forward_overlap_len = shorter_len;
                }
                trace!("forward overlap_len {forward_overlap_len}");
                if forward_overlap_len >= min_overlap_len {
                    is_overlapping = true;
                }
                // backward compare
                let mut j_index: usize = 1;
                // originally in python: up to shorter_len
                for j in 1..(shorter_len + 1) {
                    j_index = j;
                    if hap1[(hap1_len - j)..] != hap2[(hap2_len - j)..] {
                        break;
                    }
                }
                backward_overlap_len = j_index - 1;
                // fully contained
                if hap1[(hap1_len - shorter_len)..] == hap2[(hap2_len - shorter_len)..] {
                    backward_overlap_len = shorter_len;
                }
                trace!("backward overlap_len {backward_overlap_len}");
                if backward_overlap_len >= min_overlap_len {
                    is_overlapping = true;
                }
                if is_overlapping {
                    let overlap_len = cmp::max(forward_overlap_len, backward_overlap_len);
                    overlapping_haps_match
                        .entry(hap1.clone())
                        .or_default()
                        .push((hap2.clone(), overlap_len));
                    let hap1_vec = hap1.to_vec();
                    if !overlapping_haps.contains(&hap1_vec) {
                        overlapping_haps.push(hap1_vec);
                    }
                    let hap2_vec = hap2.to_vec();
                    if !overlapping_haps.contains(&hap2_vec) {
                        overlapping_haps.push(hap2_vec);
                    }
                    debug!("found overlapping haplotypes {hap1:?} and {hap2:?} with overlap_len {overlap_len}");
                }
            }
        }
    }
    Ok((overlapping_haps, overlapping_haps_match))
}

/// Identify suspicious complete alleles
/// # Arguments
/// * `all_read_edges` - read name -> vec of read nodes
/// * `complete_alleles` - complete alleles
/// * `incomplete_alleles` - incomplete alleles
/// # Returns
/// * `Vec<Vec<i32>>` - suspicious complete alleles
pub fn filter_complete_alleles(
    all_read_edges: &BTreeMap<String, Vec<i32>>,
    complete_alleles: &Vec<Vec<i32>>,
    incomplete_alleles: &Vec<Vec<i32>>,
) -> Result<Vec<Vec<i32>>, DError> {
    debug!("filter short alleles with suspicious reads");
    let mut suspicious_complete_alleles = Vec::new();
    for allele in complete_alleles {
        if is_short_complete_allele(allele) {
            if is_short_complete_allele_suspicious(allele, all_read_edges) {
                suspicious_complete_alleles.push(allele.clone());
            }
        } else if is_long_complete_allele_suspicious(allele, complete_alleles, incomplete_alleles)?
            && !suspicious_complete_alleles.contains(allele)
        {
            suspicious_complete_alleles.push(allele.clone());
        }
    }
    debug!(
        "suspicious_complete_alleles {:?}",
        suspicious_complete_alleles
    );
    Ok(suspicious_complete_alleles)
}

#[derive(Debug, Default)]
struct ShortAlleleAnalysisState {
    repeat_pos_support: BTreeMap<usize, Vec<Vec<i32>>>,
    suspicious_forward: BTreeMap<usize, Vec<Vec<i32>>>,
    suspicious_reverse: BTreeMap<usize, Vec<Vec<i32>>>,
    suspicious_reads: HashSet<Vec<i32>>,
    sites_supported_by_three: HashSet<usize>,
    sites_supported_by_four: HashSet<usize>,
    is_spanning: bool,
}

/// Determine whether a short complete allele is suspicious using read-support heuristics.
/// # Arguments
/// * `allele` - short complete allele under review
/// * `all_read_edges` - read paths used to evaluate support for the allele
/// # Returns
/// * `bool` - true if the short complete allele should be flagged as suspicious
fn is_short_complete_allele_suspicious(
    allele: &[i32],
    all_read_edges: &BTreeMap<String, Vec<i32>>,
) -> bool {
    let repeat_pos = get_repeat(allele);
    debug!("allele {allele:?} has identical-unit stretches at {repeat_pos:?}");
    let mut short_state = ShortAlleleAnalysisState::default();
    for (_read_name, read_nodes) in all_read_edges {
        let mut read_matches_forward = scan_read_matches_forward(allele, read_nodes);
        if !read_matches_forward.is_empty() {
            read_matches_forward.sort_by(|a, b| b.1.cmp(&a.1).then(b.0.cmp(&a.0)));
            let best_match = read_matches_forward[0];
            record_forward_support(
                &mut short_state,
                allele,
                read_nodes,
                best_match,
                &repeat_pos,
            );
        }
        let mut read_matches_reverse = scan_read_matches_backward(allele, read_nodes);
        if !read_matches_reverse.is_empty() {
            read_matches_reverse.sort_by(|a, b| b.1.cmp(&a.1).then(b.0.cmp(&a.0)));
            let best_match = read_matches_reverse[0];
            record_backward_support(&mut short_state, allele, read_nodes, best_match);
        }
    }
    is_short_complete_allele_suspicious_from_state(allele, &repeat_pos, &short_state)
}

/// Determine whether a short complete allele is suspicious from precomputed support state.
/// # Arguments
/// * `allele` - short complete allele under review
/// * `repeat_pos` - repeated-unit positions for the allele
/// * `short_state` - accumulated read-support evidence for the allele
/// # Returns
/// * `bool` - true if the short complete allele should be flagged as suspicious
fn is_short_complete_allele_suspicious_from_state(
    allele: &[i32],
    repeat_pos: &BTreeMap<usize, usize>,
    short_state: &ShortAlleleAnalysisState,
) -> bool {
    let allele_len = allele.len();
    debug!("is_spanning {}", short_state.is_spanning);
    debug!("repeat_pos_support {:?}", short_state.repeat_pos_support);
    let highly_repetitive = short_state.repeat_pos_support.len() < repeat_pos.len();
    debug!("allele {allele:?} highly_repetitive {highly_repetitive:?}");
    let good_support =
        allele_len >= 2 && short_state.sites_supported_by_three.len() == allele_len - 2;
    let better_support =
        allele_len >= 3 && short_state.sites_supported_by_four.len() == allele_len - 3;
    debug!(
        "allele {allele:?} good_support {good_support} sites_supported_by_three {:?}",
        short_state.sites_supported_by_three
    );
    debug!(
        "allele {allele:?} better_support {better_support} sites_supported_by_four {:?}",
        short_state.sites_supported_by_four
    );
    // 1. good if fully spanned by a long read
    if short_state.is_spanning {
        return false;
    }

    let num_suspicious_reads = short_state.suspicious_reads.len();
    debug!("num_suspicious_reads {num_suspicious_reads}");
    debug!("suspicious_forward {:?}", short_state.suspicious_forward);
    debug!("suspicious_reverse {:?}", short_state.suspicious_reverse);
    // 2. suspicious if highly repetitive
    if highly_repetitive {
        debug!("allele {allele:?} is suspicious because it is highly repetitive");
        return true;
    }

    // 3. if only two units, should have spanning reads
    if allele_len <= 4 {
        if short_state.sites_supported_by_four.is_empty() {
            debug!("allele {allele:?} is suspicious because it is two units or shorter, and has no spanning reads.");
            return true;
        } else {
            return false;
        }
    }

    // 4. we want good support at both ends
    if !short_state.sites_supported_by_three.contains(&0)
        || !short_state
            .sites_supported_by_three
            .contains(&(allele_len - 3))
    {
        debug!("allele {allele:?} is suspicious because the left side or the right side is not supported by reads linking the next two sites.");
        return true;
    }

    // 5. two suspicious sites, one forward, one reverse
    if short_state.suspicious_forward.len() == 1
        && short_state.suspicious_reverse.len() == 1
        && num_suspicious_reads > 1
        && !better_support
    {
        if let (Some((forward_pos, _)), Some((reverse_pos, _))) = (
            short_state.suspicious_forward.first_key_value(),
            short_state.suspicious_reverse.first_key_value(),
        ) {
            if !good_support {
                if *reverse_pos <= *forward_pos - 1 {
                    // 5.1 reverse pos is before forward pos and not good_support
                    debug!("allele {allele:?} is suspicious because not every site is supported by reads linking the next two sites, and it has two suspicious sites, one with forward-matching suspicous reads and one with reverse-matching suspicious reads.");
                    return true;
                }
            } else if *reverse_pos == *forward_pos - 2
                && (*reverse_pos > allele_len - 4
                    || short_state.sites_supported_by_four.contains(reverse_pos))
                && (*forward_pos > allele_len - 4
                    || short_state.sites_supported_by_four.contains(forward_pos))
            {
                // 5.2 good if reverse pos is 2 positions before forward pos and both positions are supported by reads linking the next three sites (or both sites are at the end of the allele)
                return false;
            } else if *reverse_pos <= *forward_pos - 2 {
                // 5.3 suspicious if reverse pos is 2 positions before forward pos (good_support is true here)
                debug!("allele {allele:?} is suspicious because it has two suspicious sites, one with forward-matching suspicous reads and one with reverse-matching suspicious reads. And suspicious sites are not supported by reads linking the next three sites.");
                return true;
            }
        }
    } else if short_state.suspicious_forward.is_empty() && short_state.suspicious_reverse.len() == 1
    // 6. one suspicious site only, currently only looking at reverse-matching suspicious reads
    {
        if let Some((suspicious_site, reads)) = short_state.suspicious_reverse.first_key_value() {
            let suspicious_site_num_reads = reads.len();
            debug!(
                "only one suspicious_site at index {suspicious_site:?} with {suspicious_site_num_reads} suspicious reads"
            );
            // 6.1 good if the site is linked to the next three sites
            if *suspicious_site > allele_len - 4
                || short_state
                    .sites_supported_by_four
                    .contains(suspicious_site)
            {
                return false;
            }
            // 6.2 suspicious if at the end and both this site and the previous site are not supported by reads linking the next three sites
            if *suspicious_site == allele_len - 4 {
                let prev_site = *suspicious_site - 1;
                if !short_state.sites_supported_by_four.contains(&prev_site)
                    && !short_state
                        .sites_supported_by_four
                        .contains(suspicious_site)
                    && suspicious_site_num_reads >= 2
                {
                    debug!("allele {allele:?} is suspicious because the suspicious site at the end of the allele and both this site and the previous site are not supported by reads linking the next three sites.");
                    return true;
                }
            }
            // 6.3 good if the site is linked to the next two sites
            if *suspicious_site > allele_len - 3
                || short_state
                    .sites_supported_by_three
                    .contains(suspicious_site)
            {
                return false;
            }
            // 6.4 suspicious if the site is not supported by reads linking the next two or three sites
            if suspicious_site_num_reads >= 3 {
                debug!("allele {allele:?} is suspicious because the suspicious site is not supported by reads linking the next two or three sites.");
                return true;
            }
        }
    } else if num_suspicious_reads >= 5 && !good_support {
        // 7. too many suspicious reads and not good_support
        debug!("allele {allele:?} is suspicious because it has at least 5 suspicious reads, and not every site is supported by reads linking the next two sites.");
        return true;
    }

    false
}

/// Update short-allele support state from the best forward match for one read.
/// # Arguments
/// * `state` - accumulated short-allele evidence
/// * `allele` - allele being evaluated
/// * `read_nodes` - read nodes that were matched against the allele
/// * `best_match` - best forward match as `(start_index, matched_len)`
/// * `repeat_pos` - repeated-unit positions for the allele
fn record_forward_support(
    state: &mut ShortAlleleAnalysisState,
    allele: &[i32],
    read_nodes: &[i32],
    best_match: (usize, usize),
    repeat_pos: &BTreeMap<usize, usize>,
) {
    let read_nodes_len = read_nodes.len();
    debug!(
        "allele {allele:?} read {read_nodes:?} found {} matches starting at {}",
        best_match.1, best_match.0
    );
    if best_match.1 == read_nodes_len && check_spanning(read_nodes) && !state.is_spanning {
        state.is_spanning = true;
    }
    if best_match.1 < read_nodes_len && best_match.0 + best_match.1 < allele.len() {
        state.suspicious_reads.insert(read_nodes.to_vec());
        state
            .suspicious_forward
            .entry(best_match.0 + best_match.1)
            .or_default()
            .push(read_nodes.to_vec());
    }
    let i_index = best_match.0;
    let k_index = best_match.1;
    if k_index == read_nodes_len && k_index >= 3 {
        for q in 0..(k_index - 2) {
            if !read_nodes[q..(q + 3)].contains(&0) {
                state.sites_supported_by_three.insert(i_index + q);
            }
        }
        for q in 0..(k_index - 3) {
            if !read_nodes[q..(q + 4)].contains(&0) {
                state.sites_supported_by_four.insert(i_index + q);
            }
            let this_pos = i_index + q;
            if let Some(repeat_size) = repeat_pos.get(&this_pos) {
                let end_pos = q + *repeat_size + 2;
                if end_pos <= read_nodes_len && !read_nodes[q..end_pos].contains(&0) {
                    state
                        .repeat_pos_support
                        .entry(this_pos)
                        .or_default()
                        .push(read_nodes.to_vec());
                }
            }
        }
    }
}

/// Update short-allele support state from the best backward match for one read.
/// # Arguments
/// * `state` - accumulated short-allele evidence
/// * `allele` - allele being evaluated
/// * `read_nodes` - read nodes that were matched against the allele
/// * `best_match` - best backward match as `(offset_from_end, matched_len)`
fn record_backward_support(
    state: &mut ShortAlleleAnalysisState,
    allele: &[i32],
    read_nodes: &[i32],
    best_match: (usize, usize),
) {
    let read_nodes_len = read_nodes.len();
    debug!(
        "allele {allele:?} read {read_nodes:?} found {} matches ending at {}",
        best_match.1, best_match.0
    );
    if best_match.1 == read_nodes_len && check_spanning(read_nodes) && !state.is_spanning {
        state.is_spanning = true;
    }
    if best_match.1 < read_nodes_len && best_match.0 + best_match.1 < allele.len() {
        state.suspicious_reads.insert(read_nodes.to_vec());
        state
            .suspicious_reverse
            .entry(allele.len() - best_match.0 - best_match.1)
            .or_default()
            .push(read_nodes.to_vec());
    }
}

/// Scan a read against an allele from the left side and collect matching segments.
/// # Arguments
/// * `allele` - allele being evaluated
/// * `read_nodes` - read nodes to compare against the allele
/// # Returns
/// * `Vec<(usize, usize)>` - matching segments as `(start_index, matched_len)`
fn scan_read_matches_forward(allele: &[i32], read_nodes: &[i32]) -> Vec<(usize, usize)> {
    let allele_len = allele.len();
    let read_nodes_len = read_nodes.len();
    let mut read_matches_forward = Vec::new();
    for i in 0..allele_len {
        let i_index = i;
        let upper_bound = cmp::min(allele_len - i, read_nodes_len) + 1;
        for k in 2..upper_bound {
            let k_index = k;
            let nodes_in_allele = &allele[i..(i + k)];
            let nodes_in_reads = &read_nodes[0..k];
            let mut match_allele = Vec::new();
            for (j, read_node) in nodes_in_reads.iter().enumerate() {
                if *read_node == 0 {
                    match_allele.push(0)
                } else if *read_node == nodes_in_allele[j] {
                    match_allele.push(1);
                } else {
                    match_allele.push(-1);
                }
            }
            let match_count = match_allele
                .iter()
                .filter(|x| **x == 1)
                .collect::<Vec<_>>()
                .len();
            trace!(
                "allele {allele:?} read {read_nodes:?} {i_index} {k_index} {nodes_in_allele:?} {nodes_in_reads:?} {match_count}"
            );
            if !match_allele.contains(&(-1)) && match_count > 1 {
                read_matches_forward.push((i_index, k_index));
            }
            if match_allele.contains(&(-1)) {
                break;
            }
        }
    }
    read_matches_forward
}

/// Scan a read against an allele from the right side and collect matching segments.
/// # Arguments
/// * `allele` - allele being evaluated
/// * `read_nodes` - read nodes to compare against the allele
/// # Returns
/// * `Vec<(usize, usize)>` - matching segments as `(offset_from_end, matched_len)`
fn scan_read_matches_backward(allele: &[i32], read_nodes: &[i32]) -> Vec<(usize, usize)> {
    let allele_len = allele.len();
    let read_nodes_len = read_nodes.len();
    let mut read_matches_reverse = Vec::new();
    for i in 0..allele_len {
        let i_index = i;
        let upper_bound = cmp::min(allele_len - i, read_nodes_len) + 1;
        for k in 2..upper_bound {
            let k_index = k;
            let nodes_in_allele = &allele[(allele_len - i - k)..(allele_len - i)];
            let nodes_in_reads = &read_nodes[(read_nodes_len - k)..read_nodes_len];
            let mut match_allele = Vec::new();
            for (j, read_node) in nodes_in_reads.iter().enumerate() {
                if *read_node == 0 {
                    match_allele.push(0)
                } else if *read_node == nodes_in_allele[j] {
                    match_allele.push(1);
                } else {
                    match_allele.push(-1);
                }
            }
            let match_count = match_allele
                .iter()
                .filter(|x| **x == 1)
                .collect::<Vec<_>>()
                .len();
            trace!(
                "allele {allele:?} read {read_nodes:?} {i_index} {k_index} {nodes_in_allele:?} {nodes_in_reads:?} {match_count}"
            );
            if !match_allele.contains(&(-1)) && match_count > 1 {
                read_matches_reverse.push((i_index, k_index));
            }
            if match_allele.contains(&(-1)) {
                break;
            }
        }
    }
    read_matches_reverse
}

/// Return whether an allele should use the short-allele suspicion heuristics.
/// # Arguments
/// * `allele` - allele to classify by size
/// # Returns
/// * `bool` - true if the allele should use short-allele logic
fn is_short_complete_allele(allele: &[i32]) -> bool {
    allele.len() <= 16
}

/// Determine whether a long complete allele is suspicious based on overlap with
/// other complete or incomplete haplotypes.
/// # Arguments
/// * `allele` - complete allele under review
/// * `complete_alleles` - all complete alleles
/// * `incomplete_alleles` - all incomplete alleles
/// # Returns
/// * `bool` - true if the long complete allele is suspicious
fn is_long_complete_allele_suspicious(
    allele: &[i32],
    complete_alleles: &[Vec<i32>],
    incomplete_alleles: &[Vec<i32>],
) -> Result<bool, DError> {
    let mut haps_to_assess = complete_alleles.to_vec();
    haps_to_assess.extend_from_slice(incomplete_alleles);
    let overlapping_haps = find_overlapping_alleles(haps_to_assess.clone(), Some(10))?.1;
    let overlapping_haps_loose = find_overlapping_alleles(haps_to_assess, Some(7))?.1;
    let Some(hap1_overlaps) = overlapping_haps.get(allele) else {
        return Ok(false);
    };
    let Some(this_overlaps_loose) = overlapping_haps_loose.get(allele) else {
        return Ok(false);
    };

    let this_overlaps_loose_count = this_overlaps_loose.len();
    let mut overlap_count = 0;
    for (hap2, overlap_len) in this_overlaps_loose {
        if !redundant_haplotype_allowed(allele, hap2, overlap_len)? {
            overlap_count += 1;
        }
    }
    if overlap_count > 0 && this_overlaps_loose_count > 1 {
        debug!(
            "Complete haplotype {allele:?} is suspicious because it has {this_overlaps_loose_count} loose overlapping haplotypes, {overlap_count} of them are not redundant"
        );
        return Ok(true);
    }

    for (hap2, overlap_len) in hap1_overlaps {
        if !redundant_haplotype_allowed(allele, hap2, overlap_len)? {
            debug!(
                "Complete haplotype {allele:?} is suspicious because it is overlapping with another haplotype {hap2:?}"
            );
            return Ok(true);
        }
    }

    Ok(false)
}

/// Identify redundant haplotypes
/// # Arguments
/// * `hap1` - haplotype 1
/// * `hap2` - haplotype 2
/// * `overlap_len` - overlap length
/// # Returns
/// * `bool` - true if redundant, false otherwise
pub fn redundant_haplotype_allowed(
    hap1: &[i32],
    hap2: &[i32],
    overlap_len: &usize,
) -> Result<bool, DError> {
    let hap1_len = hap1.len() as i32;
    let hap2_len = hap2.len() as i32;
    let hap1_first = hap1
        .first()
        .ok_or_else(|| missing_data_error("first haplotype node", format!("{hap1:?}")))?;
    let hap1_last = hap1
        .last()
        .ok_or_else(|| missing_data_error("last haplotype node", format!("{hap1:?}")))?;
    let hap2_first = hap2
        .first()
        .ok_or_else(|| missing_data_error("first haplotype node", format!("{hap2:?}")))?;
    let hap2_last = hap2
        .last()
        .ok_or_else(|| missing_data_error("last haplotype node", format!("{hap2:?}")))?;
    // if the other haplotype is complete and is similar length, skip
    if *hap2_first < 0 && *hap2_first > -10 && *hap2_last <= -10 {
        if (hap1_len - hap2_len).abs() <= 2
            && *hap1_first < 0
            && *hap1_first > -10
            && *hap1_last <= -10
        {
            return Ok(true);
        }
    } else {
        // if the other haplotype is just one unit longer than the overlap, skip
        // the extra unit could be a redundant fingerprint
        if hap2_len - (*overlap_len as i32) < 2 {
            return Ok(true);
        }
    }
    // only differ by one fingerprint
    let hap1_len = hap1.len();
    let hap2_len = hap2.len();
    let min_len = cmp::min(hap1_len, hap2_len) as usize;
    let mut differ_by_one = false;
    for i in 1..(min_len - 1) {
        if &hap1[0..i] == &hap2[0..i] && &hap1[i + 1..min_len] == &hap2[i + 1..min_len] {
            differ_by_one = true;
            break;
        }
    }
    if differ_by_one {
        return Ok(true);
    }
    let mut differ_by_one = false;
    for i in 1..(min_len - 1) {
        if &hap1[(hap1_len - min_len)..(hap1_len - i - 1)]
            == &hap2[(hap2_len - min_len)..(hap2_len - i - 1)]
            && &hap1[(hap1_len - i)..] == &hap2[(hap2_len - i)..]
        {
            differ_by_one = true;
            break;
        }
    }
    if differ_by_one {
        return Ok(true);
    }
    Ok(false)
}

/// Identify a stretch of identical units in an allele
/// # Arguments
/// * `allele` - allele
/// # Returns
/// * `BTreeMap<usize, usize>` - starting index of the identical stretch, length of the stretch
pub fn get_repeat(allele: &[i32]) -> BTreeMap<usize, usize> {
    let mut repeat_pos = BTreeMap::new();
    let allele_len = allele.len();
    for i in 0..(allele_len - 1) {
        if allele[i] != allele[i + 1] {
            for k in 4..(allele_len - i) {
                if allele[i + k - 1] != allele[i + k]
                    && allele[(i + 1)..(i + k)]
                        .iter()
                        .collect::<HashSet<&i32>>()
                        .len()
                        == 1
                {
                    repeat_pos.insert(i, k - 1);
                }
            }
        }
    }
    repeat_pos
}

/// Check if a read is spanning the entire allele
/// # Arguments
/// * `read_nodes` - read nodes
/// # Returns
/// * `bool` - true if spanning, false otherwise
pub fn check_spanning(read_nodes: &[i32]) -> bool {
    let (Some(first_node), Some(last_node)) = (read_nodes.first(), read_nodes.last()) else {
        return false;
    };
    *first_node < 0 && *last_node < 0
}

/// Return whether two haplotypes (can be different length) are matching.
/// `hap1` is a read, `hap2` is a haplotype.
/// has to overlap `key_site` on `hap2` if provided.
/// # Arguments
/// * `hap1` - haplotype 1
/// * `hap2` - haplotype 2
/// * `key_site` - key site
/// # Returns
/// * `bool` - true if matching, false otherwise
pub fn two_haplotypes_are_matching(
    hap1: &[i32],
    hap2: &[i32],
    key_site: Option<(usize, usize)>,
) -> bool {
    let hap1_len = hap1.len();
    let hap2_len = hap2.len();
    if let Some((nstart, nend)) = key_site {
        for i in 0..hap1_len {
            // hap1/read      -------------
            // hap2/hapl ------------
            //          |--k--|offset|
            // hap1/read       -------
            // hap2/hapl ---------------
            //           |--k--|offset|
            if i == 0 {
                for k in 0..hap2_len {
                    let offset_index = cmp::min(hap2_len - k, hap1_len);
                    if offset_index > 1 {
                        // here we are requiring a match until the end of the haplotype (hap2)
                        if k <= nstart && k + offset_index >= nend {
                            let test_hap1 = &hap1[..offset_index];
                            let test_hap2 = &hap2[k..(k + offset_index)];
                            let (hap_match, mismatch) =
                                compare_two_haps_same_length(test_hap1, test_hap2);
                            let match_key_sites = &hap_match[(nstart - k)..(nend - k)];
                            trace!("hap1 {hap1:?} hap2 {hap2:?} key_site {key_site:?} {i} {k} offset_index {offset_index} {test_hap1:?} {test_hap2:?}");
                            if !match_key_sites.contains(&0) && mismatch == 0 {
                                return true;
                            }
                        }
                    }
                }
            }
            //           |-i-|
            // hap1/read  ----------
            // hap2/hapl     ---------------
            //               |offset|
            //           |-i-|
            // hap1/read  -------------
            // hap2/hapl     --------
            //               |offset|
            else {
                let offset_index = cmp::min(hap1_len - i, hap2_len);
                if offset_index > 1 {
                    // again requiring the match to extend all the way to the end of hap2
                    if offset_index >= nend {
                        // this was originally implemented in python
                        //let test_hap1 = &hap1[(hap1_len-offset_index)..];
                        let test_hap1 = &hap1[i..(i + offset_index)];
                        let test_hap2 = &hap2[..offset_index];
                        let (hap_match, mismatch) =
                            compare_two_haps_same_length(test_hap1, test_hap2);
                        let match_key_sites = &hap_match[nstart..nend];
                        trace!("hap1 {hap1:?} hap2 {hap2:?} key_site {key_site:?} {i} no k offset_index {offset_index} {test_hap1:?} {test_hap2:?}");
                        if !match_key_sites.contains(&0) && mismatch == 0 {
                            return true;
                        }
                    }
                }
            }
        }
    } else {
        for i in 0..hap1_len {
            if i == 0 {
                for k in 0..hap2_len {
                    let offset_index = cmp::min(hap2_len - k, hap1_len);
                    if offset_index > 1 {
                        // here we are requiring a match until the end of the haplotype (hap2)
                        let test_hap1 = &hap1[..offset_index];
                        let test_hap2 = &hap2[k..(k + offset_index)];
                        let (hap_match, mismatch) =
                            compare_two_haps_same_length(test_hap1, test_hap2);
                        if hap_match.iter().sum::<i32>() >= 2 && mismatch == 0 {
                            return true;
                        }
                    }
                    if offset_index >= 1 && k == 0 {
                        let test_hap1 = &hap1[..offset_index];
                        let test_hap2 = &hap2[k..(k + offset_index)];
                        // allow a single match of a starting unit for d4z4
                        if test_hap1[0] < 0 && test_hap1[0] > -10 {
                            let (_hap_match, mismatch) =
                                compare_two_haps_same_length(test_hap1, test_hap2);
                            if test_hap1[0] == test_hap2[0] && mismatch == 0 {
                                return true;
                            }
                        }
                    }
                }
            } else {
                let offset_index = cmp::min(hap1_len - i, hap2_len);
                if offset_index > 1 {
                    // this was originally implemented in python
                    //let test_hap1 = &hap1[(hap1_len-offset_index)..];
                    let test_hap1 = &hap1[i..(i + offset_index)];
                    let test_hap2 = &hap2[..offset_index];
                    let (hap_match, mismatch) = compare_two_haps_same_length(test_hap1, test_hap2);
                    if hap_match.iter().sum::<i32>() >= 2 && mismatch == 0 {
                        return true;
                    }
                }
            }
        }
    }
    false
}

/// Calculate number of matching/mismatching bases between two haplotypes of same length
/// # Arguments
/// * `hap1` - haplotype 1
/// * `hap2` - haplotype 2
/// # Returns
/// * `Vec<i32>` - matching bases
/// * `i32` - number of mismatching bases
pub fn compare_two_haps_same_length(hap1: &[i32], hap2: &[i32]) -> (Vec<i32>, i32) {
    let hap1_len = hap1.len();
    let hap2_len = hap2.len();
    if hap1_len != hap2_len {
        error!(
            "two haplotypes {:?} and {:?} are not the same length",
            hap1, hap2
        );
    }
    let mut hap_match = Vec::new();
    let mut mismatch = 0;
    for (j, base1) in hap1.iter().enumerate() {
        let base2 = hap2[j];
        hap_match.push(0);
        if *base1 != 0 && base2 != 0 {
            if *base1 == base2 {
                hap_match[j] = 1;
            } else {
                mismatch += 1;
            }
        }
    }
    (hap_match, mismatch)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_is_short_complete_allele_suspicious_from_state_spanning_is_not_suspicious() {
        // scenario 1
        let allele = vec![-2, 2, 3, -10];
        let repeat_pos = BTreeMap::new();
        let short_state = ShortAlleleAnalysisState {
            is_spanning: true,
            ..Default::default()
        };

        let is_suspicious =
            is_short_complete_allele_suspicious_from_state(&allele, &repeat_pos, &short_state);

        assert!(!is_suspicious);
    }

    #[test]
    fn test_is_short_complete_allele_suspicious_from_state_highly_repetitive_is_suspicious() {
        // scenario 2
        let allele = vec![-2, 2, 2, 2, -10];
        let repeat_pos = BTreeMap::from([(1, 3)]);
        let short_state = ShortAlleleAnalysisState::default();

        let is_suspicious =
            is_short_complete_allele_suspicious_from_state(&allele, &repeat_pos, &short_state);

        assert!(is_suspicious);
    }

    #[test]
    fn test_is_short_complete_allele_suspicious_from_state_short_allele_without_support_by_four_is_suspicious(
    ) {
        // scenario 3, short allele, supported_by_four is empty
        let allele = vec![-2, 2, 3, -10];
        let repeat_pos = BTreeMap::new();
        let short_state = ShortAlleleAnalysisState {
            sites_supported_by_three: HashSet::from([0, 1]),
            ..Default::default()
        };
        let is_suspicious =
            is_short_complete_allele_suspicious_from_state(&allele, &repeat_pos, &short_state);
        assert!(is_suspicious);

        let short_state = ShortAlleleAnalysisState {
            sites_supported_by_three: HashSet::from([0, 1]),
            sites_supported_by_four: HashSet::from([0]),
            ..Default::default()
        };
        let is_suspicious =
            is_short_complete_allele_suspicious_from_state(&allele, &repeat_pos, &short_state);
        assert!(!is_suspicious);
    }

    #[test]
    fn test_is_short_complete_allele_suspicious_from_state_missing_end_support_is_suspicious() {
        // scenario 4, requiring both ends
        let allele = vec![-2, 2, 3, 4, 5, -10];
        let repeat_pos = BTreeMap::new();
        let short_state = ShortAlleleAnalysisState {
            sites_supported_by_three: HashSet::from([0, 1, 2]),
            ..Default::default()
        };
        let is_suspicious =
            is_short_complete_allele_suspicious_from_state(&allele, &repeat_pos, &short_state);
        assert!(is_suspicious);

        let short_state = ShortAlleleAnalysisState {
            sites_supported_by_three: HashSet::from([0, 1, 3]),
            ..Default::default()
        };
        let is_suspicious =
            is_short_complete_allele_suspicious_from_state(&allele, &repeat_pos, &short_state);
        assert!(!is_suspicious);
    }

    #[test]
    fn test_is_short_complete_allele_suspicious_from_state_many_suspicious_reads_without_good_support_is_suspicious(
    ) {
        // scenario 7, many suspicious reads and not good_support
        let allele = vec![-2, 2, 3, 4, 5, -10];
        let repeat_pos = BTreeMap::new();
        let short_state = ShortAlleleAnalysisState {
            suspicious_reads: HashSet::from([
                vec![1, 2],
                vec![2, 3],
                vec![3, 4],
                vec![4, 5],
                vec![5, 6],
            ]),
            sites_supported_by_three: HashSet::from([0, 1, 2]),
            ..Default::default()
        };

        let is_suspicious =
            is_short_complete_allele_suspicious_from_state(&allele, &repeat_pos, &short_state);

        assert!(is_suspicious);
    }

    #[test]
    fn test_is_short_complete_allele_suspicious_from_state_two_suspicious_sites_without_good_support_is_suspicious(
    ) {
        // scenario 5, one forward suspicious site and one reverse suspicious site
        let allele = vec![-2, 1, 2, 3, 4, 5, 6, -10];
        let repeat_pos = BTreeMap::new();
        let short_state = ShortAlleleAnalysisState {
            suspicious_forward: BTreeMap::from([(4, vec![vec![3, 4, 7]])]),
            suspicious_reverse: BTreeMap::from([(3, vec![vec![7, 3, 4]])]),
            suspicious_reads: HashSet::from([vec![3, 4, 7], vec![7, 3, 4]]),
            sites_supported_by_three: HashSet::from([0, 1, 5]),
            sites_supported_by_four: HashSet::from([0, 1, 2, 3, 4]), // better_support
            ..Default::default()
        };
        let is_suspicious =
            is_short_complete_allele_suspicious_from_state(&allele, &repeat_pos, &short_state);
        assert!(!is_suspicious);

        let short_state = ShortAlleleAnalysisState {
            suspicious_forward: BTreeMap::from([(4, vec![vec![3, 4, 7]])]),
            suspicious_reverse: BTreeMap::from([(3, vec![vec![7, 3, 4]])]),
            suspicious_reads: HashSet::from([vec![3, 4, 7], vec![7, 3, 4]]),
            sites_supported_by_three: HashSet::from([0, 1]),
            sites_supported_by_four: HashSet::from([0]),
            ..Default::default()
        };
        let is_suspicious =
            is_short_complete_allele_suspicious_from_state(&allele, &repeat_pos, &short_state);
        assert!(is_suspicious);

        // good support, site index differ by one
        let short_state = ShortAlleleAnalysisState {
            suspicious_forward: BTreeMap::from([(4, vec![vec![3, 4, 7]])]),
            suspicious_reverse: BTreeMap::from([(3, vec![vec![7, 3, 4]])]),
            suspicious_reads: HashSet::from([vec![3, 4, 7], vec![7, 3, 4]]),
            sites_supported_by_three: HashSet::from([0, 1, 2, 3, 4, 5]), // good_support
            sites_supported_by_four: HashSet::from([0]),
            ..Default::default()
        };
        let is_suspicious =
            is_short_complete_allele_suspicious_from_state(&allele, &repeat_pos, &short_state);
        assert!(!is_suspicious);

        // good support, site index differ by more than one
        let short_state = ShortAlleleAnalysisState {
            suspicious_forward: BTreeMap::from([(4, vec![vec![3, 4, 7]])]),
            suspicious_reverse: BTreeMap::from([(2, vec![vec![7, 2, 3]])]),
            suspicious_reads: HashSet::from([vec![3, 4, 7], vec![7, 2, 3]]),
            sites_supported_by_three: HashSet::from([0, 1, 2, 3, 4, 5]), // good_support
            sites_supported_by_four: HashSet::from([0]),
            ..Default::default()
        };
        let is_suspicious =
            is_short_complete_allele_suspicious_from_state(&allele, &repeat_pos, &short_state);
        assert!(is_suspicious);

        // good support, site index differ by two, but both sites in suported_by_four
        let short_state = ShortAlleleAnalysisState {
            suspicious_forward: BTreeMap::from([(4, vec![vec![3, 4, 7]])]),
            suspicious_reverse: BTreeMap::from([(2, vec![vec![7, 2, 3]])]),
            suspicious_reads: HashSet::from([vec![3, 4, 7], vec![7, 2, 3]]),
            sites_supported_by_three: HashSet::from([0, 1, 2, 3, 4, 5]), // good_support
            sites_supported_by_four: HashSet::from([0, 2, 4]),
            ..Default::default()
        };
        let is_suspicious =
            is_short_complete_allele_suspicious_from_state(&allele, &repeat_pos, &short_state);
        assert!(!is_suspicious);

        // good support, site index differ by two, but both sites in suported_by_four or index towards the end of allele
        let short_state = ShortAlleleAnalysisState {
            suspicious_forward: BTreeMap::from([(5, vec![vec![4, 5, 7]])]), // 5 is towards the end and is in support_by_three
            suspicious_reverse: BTreeMap::from([(3, vec![vec![7, 3, 4]])]),
            suspicious_reads: HashSet::from([vec![4, 5, 7], vec![7, 3, 4]]),
            sites_supported_by_three: HashSet::from([0, 1, 2, 3, 4, 5]), // good_support
            sites_supported_by_four: HashSet::from([0, 3]),
            ..Default::default()
        };
        let is_suspicious =
            is_short_complete_allele_suspicious_from_state(&allele, &repeat_pos, &short_state);
        assert!(!is_suspicious);
    }

    #[test]
    fn test_get_repeat() {
        let allele = &vec![1, 2, 2, 2, 3];
        let repeat_pos = get_repeat(allele);
        assert!(repeat_pos.contains_key(&0));
        assert_eq!(repeat_pos[&0], 3);

        let allele = &vec![0, 1, 2, 2, 2, 2, 3];
        let repeat_pos = get_repeat(allele);
        assert!(repeat_pos.contains_key(&1));
        assert_eq!(repeat_pos[&1], 4);
    }

    #[test]
    fn test_find_overlapping_alleles() {
        let haps = vec![vec![1, 2, 3, 4, 5], vec![6, 7, 8], vec![9, 2, 3, 4, 5]];
        let (ol_haps, ovl_len) = find_overlapping_alleles(haps.clone(), None).unwrap();
        assert_eq!(ol_haps, vec![vec![1, 2, 3, 4, 5], vec![9, 2, 3, 4, 5]]);
        let ovl_len = ovl_len.get(&haps[0]).unwrap();
        assert!(ovl_len.contains(&(haps[2].clone(), 4)));

        // fully contained #1
        let haps = vec![vec![1, 2, 3, 4, 5], vec![6, 7, 8], vec![2, 3, 4, 5]];
        let (ol_haps, ovl_len) = find_overlapping_alleles(haps.clone(), None).unwrap();
        assert_eq!(ol_haps, vec![vec![1, 2, 3, 4, 5], vec![2, 3, 4, 5]]);
        let ovl_len = ovl_len.get(&haps[0]).unwrap();
        assert!(ovl_len.contains(&(haps[2].clone(), 4)));

        // fully contained #2
        let haps = vec![vec![1, 2, 3, 4, 5], vec![6, 7, 8], vec![1, 2, 3, 4]];
        let (ol_haps, ovl_len) = find_overlapping_alleles(haps.clone(), None).unwrap();
        assert_eq!(ol_haps, vec![vec![1, 2, 3, 4, 5], vec![1, 2, 3, 4]]);
        let ovl_len = ovl_len.get(&haps[0]).unwrap();
        assert!(ovl_len.contains(&(haps[2].clone(), 4)));

        // fully contained #3
        let haps = vec![vec![1, 2, 3, 4, 5, 6], vec![6, 7, 8], vec![1, 2, 3, 4, 5]];
        let (ol_haps, ovl_len) = find_overlapping_alleles(haps.clone(), Some(5)).unwrap();
        assert_eq!(ol_haps, vec![vec![1, 2, 3, 4, 5, 6], vec![1, 2, 3, 4, 5]]);
        let ovl_len = ovl_len.get(&haps[0]).unwrap();
        assert!(ovl_len.contains(&(haps[2].clone(), 5)));

        let haps = vec![vec![1, 2, 3, 4, 5], vec![6, 7, 8], vec![1, 3, 4, 5]];
        let (ol_haps, _ovl_len) = find_overlapping_alleles(haps.clone(), None).unwrap();
        assert!(ol_haps.is_empty());

        let haps = vec![vec![1, 2, 3, 4, 5], vec![6, 7, 8], vec![1, 3, 4, 5, 9, 10]];
        let (ol_haps, _ovl_len) = find_overlapping_alleles(haps.clone(), None).unwrap();
        assert!(ol_haps.is_empty());

        let haps = vec![vec![1, 2, 3, 4, 6], vec![6, 7, 8], vec![1, 2, 3, 4, 5]];
        let (ol_haps, ovl_len) = find_overlapping_alleles(haps.clone(), None).unwrap();
        assert_eq!(ol_haps, vec![vec![1, 2, 3, 4, 6], vec![1, 2, 3, 4, 5]]);
        let ovl_len = ovl_len.get(&haps[0]).unwrap();
        assert!(ovl_len.contains(&(haps[2].clone(), 4)));

        let haps = vec![vec![1, 2, 3, 4, 6], vec![6, 7, 8], vec![1, 2, 3, 4, 5]];
        let (ol_haps, _ovl_len) = find_overlapping_alleles(haps.clone(), Some(5)).unwrap();
        assert!(ol_haps.is_empty());

        let haps = vec![
            vec![1, 2, 3, 4, 6, 9],
            vec![6, 7, 8],
            vec![1, 2, 3, 4, 6, 5],
        ];
        let (ol_haps, ovl_len) = find_overlapping_alleles(haps.clone(), Some(5)).unwrap();
        assert_eq!(
            ol_haps,
            vec![vec![1, 2, 3, 4, 6, 9], vec![1, 2, 3, 4, 6, 5]]
        );
        let ovl_len = ovl_len.get(&haps[0]).unwrap();
        assert!(ovl_len.contains(&(haps[2].clone(), 5)));

        let haps = vec![
            vec![1, 2, 3, 4, 6, 9, 10],
            vec![6, 7, 8],
            vec![1, 2, 3, 4, 6, 5, 11],
        ];
        let (ol_haps, ovl_len) = find_overlapping_alleles(haps.clone(), Some(5)).unwrap();
        assert_eq!(
            ol_haps,
            vec![vec![1, 2, 3, 4, 6, 9, 10], vec![1, 2, 3, 4, 6, 5, 11]]
        );
        let ovl_len = ovl_len.get(&haps[0]).unwrap();
        assert!(ovl_len.contains(&(haps[2].clone(), 5)));

        let haps = vec![
            vec![10, 13, 2, 3, 4, 6, 9],
            vec![6, 7, 8],
            vec![11, 12, 2, 3, 4, 6, 9],
        ];
        let (ol_haps, ovl_len) = find_overlapping_alleles(haps.clone(), Some(5)).unwrap();
        assert_eq!(
            ol_haps,
            vec![vec![10, 13, 2, 3, 4, 6, 9], vec![11, 12, 2, 3, 4, 6, 9]]
        );
        let ovl_len = ovl_len.get(&haps[0]).unwrap();
        assert!(ovl_len.contains(&(haps[2].clone(), 5)));

        let haps = vec![
            vec![10, 2, 3, 4, 6, 9],
            vec![6, 7, 8],
            vec![11, 2, 3, 4, 6, 9],
        ];
        let (ol_haps, ovl_len) = find_overlapping_alleles(haps.clone(), Some(5)).unwrap();
        assert_eq!(
            ol_haps,
            vec![vec![10, 2, 3, 4, 6, 9], vec![11, 2, 3, 4, 6, 9]]
        );
        let ovl_len = ovl_len.get(&haps[0]).unwrap();
        assert!(ovl_len.contains(&(haps[2].clone(), 5)));
    }

    #[test]
    fn test_two_haplotypes_are_matching() {
        let hap1 = vec![4, 5, 6];
        let hap2 = vec![1, 2, 3, 0, 4, 5, 6];
        let is_matching = two_haplotypes_are_matching(&hap1, &hap2, Some((4, 7)));
        assert!(is_matching);

        let hap1 = vec![4, 5, 6];
        let hap2 = vec![1, 2, 3, 0, 4, 5, 6];
        let is_matching = two_haplotypes_are_matching(&hap1, &hap2, Some((5, 7)));
        assert!(is_matching);

        let hap1 = vec![4, 5, 6];
        let hap2 = vec![1, 2, 3, 0, 4, 5, 6];
        let is_matching = two_haplotypes_are_matching(&hap1, &hap2, Some((3, 7)));
        assert!(!is_matching);

        let hap1 = vec![3, 4, 5, 6];
        let hap2 = vec![1, 2, 3, 0, 4, 5, 6];
        let is_matching = two_haplotypes_are_matching(&hap1, &hap2, Some((3, 7)));
        assert!(!is_matching);

        let hap1 = vec![3, 4, 5, 6];
        let hap2 = vec![1, 2, 3, 2, 4, 5, 6];
        let is_matching = two_haplotypes_are_matching(&hap1, &hap2, Some((3, 7)));
        assert!(!is_matching);

        let hap1 = vec![4, 5, 6, 7, 8];
        let hap2 = vec![1, 2, 3, 0, 4, 5, 6];
        let is_matching = two_haplotypes_are_matching(&hap1, &hap2, Some((4, 7)));
        assert!(is_matching);

        let hap1 = vec![4, 5, 6, 7, 8];
        let hap2 = vec![5, 6];
        let is_matching = two_haplotypes_are_matching(&hap1, &hap2, Some((0, 2)));
        assert!(is_matching);

        let hap1 = vec![4, 5, 6, 7, 8];
        let hap2 = vec![7, 8, 1, 2];
        let is_matching = two_haplotypes_are_matching(&hap1, &hap2, Some((0, 4)));
        assert!(!is_matching);

        let hap1 = vec![4, 5, 6, 7, 8];
        let hap2 = vec![6, 7, 8, 1, 2, 3];
        let is_matching = two_haplotypes_are_matching(&hap1, &hap2, Some((0, 3)));
        assert!(is_matching);

        let hap1 = vec![4, 5, 6, 7, 8];
        let hap2 = vec![4, 7, 8, 1, 2, 3];
        let is_matching = two_haplotypes_are_matching(&hap1, &hap2, Some((0, 3)));
        assert!(!is_matching);

        let hap1 = vec![4, 7, 8];
        let hap2 = vec![4, 7, 8, 1, 2, 3];
        let is_matching = two_haplotypes_are_matching(&hap1, &hap2, None);
        assert!(is_matching);

        let hap1 = vec![4, 6, 8];
        let hap2 = vec![4, 7, 8, 1, 2, 3];
        let is_matching = two_haplotypes_are_matching(&hap1, &hap2, None);
        assert!(!is_matching);

        let hap1 = vec![4, 0, 8];
        let hap2 = vec![4, 7, 8, 1, 2, 3];
        let is_matching = two_haplotypes_are_matching(&hap1, &hap2, None);
        assert!(is_matching);

        let hap1 = vec![0, 8];
        let hap2 = vec![4, 7, 8, 1, 2, 3];
        let is_matching = two_haplotypes_are_matching(&hap1, &hap2, None);
        assert!(!is_matching);

        let hap1 = vec![4];
        let hap2 = vec![4, 7, 8, 1, 2, 3];
        let is_matching = two_haplotypes_are_matching(&hap1, &hap2, None);
        assert!(!is_matching);

        let hap1 = vec![-2];
        let hap2 = vec![-2, 1, 2, 3];
        let is_matching = two_haplotypes_are_matching(&hap1, &hap2, None);
        assert!(is_matching);

        let hap1 = vec![-2, 0];
        let hap2 = vec![-2, 1, 2, 3];
        let is_matching = two_haplotypes_are_matching(&hap1, &hap2, None);
        assert!(is_matching);

        let hap1 = vec![-10];
        let hap2 = vec![1, 2, 3, -10];
        let is_matching = two_haplotypes_are_matching(&hap1, &hap2, None);
        assert!(!is_matching);

        let hap1 = vec![-10];
        let hap2 = vec![-10, 1, 2, 3];
        let is_matching = two_haplotypes_are_matching(&hap1, &hap2, None);
        assert!(!is_matching);
    }

    #[test]
    fn test_compare_two_haps_same_length() {
        let hap1 = vec![1, 2, 3, 0];
        let hap2 = vec![1, 4, 0, 0];
        let (hap_match, mismatch) = compare_two_haps_same_length(&hap1, &hap2);
        assert_eq!(mismatch, 1);
        assert_eq!(hap_match, vec![1, 0, 0, 0]);
    }

    #[test]
    fn test_redundant_haplotype_allowed() {
        // if the other haplotype is complete and is similar length
        let hap1 = vec![-2, 1, 2, 3, 4, 5, 6, -10];
        let hap2 = vec![-2, 1, 2, 3, 4, 5, -10];
        let overlap_len = 6;
        let is_allowed = redundant_haplotype_allowed(&hap1, &hap2, &overlap_len).unwrap();
        assert!(is_allowed);

        let hap1 = vec![-2, 1, 2, 3, 4, 5, 6, 7, -10];
        let hap2 = vec![-2, 1, 2, 3, 4, 5, -10];
        let overlap_len = 6;
        let is_allowed = redundant_haplotype_allowed(&hap1, &hap2, &overlap_len).unwrap();
        assert!(is_allowed);

        let hap1 = vec![-2, 1, 2, 3, 4, 5, 6, 7, 8, -10];
        let hap2 = vec![-2, 1, 2, 3, 4, 5, -10];
        let overlap_len = 6;
        let is_allowed = redundant_haplotype_allowed(&hap1, &hap2, &overlap_len).unwrap();
        assert!(is_allowed == false);

        // cases where hap2 is not complete
        let hap1 = vec![-2, 1, 2, 3, 4, 5, 6, 7, 8, -10];
        let hap2 = vec![-2, 1, 2, 3, 4, 5, 9];
        let overlap_len = 6;
        let is_allowed = redundant_haplotype_allowed(&hap1, &hap2, &overlap_len).unwrap();
        assert!(is_allowed);

        let hap1 = vec![-2, 1, 2, 3, 4, 5, 6, 7, 8, -10];
        let hap2 = vec![-2, 1, 2, 3, 4, 5, 9, 10];
        let overlap_len = 6;
        let is_allowed = redundant_haplotype_allowed(&hap1, &hap2, &overlap_len).unwrap();
        assert!(is_allowed == false);

        let hap1 = vec![-2, 1, 2, 3, 4, 5, 6, 7, 8, -10];
        let hap2 = vec![-2, 1, 2, 3, 4, 9, 6];
        let overlap_len = 5;
        let is_allowed = redundant_haplotype_allowed(&hap1, &hap2, &overlap_len).unwrap();
        assert!(is_allowed);

        let hap1 = vec![-2, 1, 2, 3, 4, 5, 6, 7, 8, -10];
        let hap2 = vec![-2, 1, 2, 3, 9, 5, 6];
        let overlap_len = 4;
        let is_allowed = redundant_haplotype_allowed(&hap1, &hap2, &overlap_len).unwrap();
        assert!(is_allowed);

        let hap1 = vec![-2, 1, 2, 3, 4, 5, 6, 7, 8, -10];
        let hap2 = vec![-2, 1, 2, 3, 10, 9, 6];
        let overlap_len = 4;
        let is_allowed = redundant_haplotype_allowed(&hap1, &hap2, &overlap_len).unwrap();
        assert!(is_allowed == false);

        let hap1 = vec![-2, 1, 2, 3, 4, 5, 6, 7, 8, -10];
        let hap2 = vec![3, 9, 5, 6, 7, 8, -10];
        let overlap_len = 5;
        let is_allowed = redundant_haplotype_allowed(&hap1, &hap2, &overlap_len).unwrap();
        assert!(is_allowed);

        let hap1 = vec![-2, 1, 2, 3, 4, 5, 6, 7, 8, -10];
        let hap2 = vec![3, 9, 10, 6, 7, 8, -10];
        let overlap_len = 4;
        let is_allowed = redundant_haplotype_allowed(&hap1, &hap2, &overlap_len).unwrap();
        assert!(is_allowed == false);

        let hap1 = vec![
            -2, 26, 23, 9, 4, 20, 20, 2, 12, 20, 20, 22, 25, 13, 15, 10, 14, 5, -11,
        ];
        let hap2 = vec![-2, 26, 23, 9, 4, 20, 20, 2, 12, 20, 20, 22, 25, 13, 22, 25];
        let overlap_len = 14;
        let is_allowed = redundant_haplotype_allowed(&hap1, &hap2, &overlap_len).unwrap();
        assert!(is_allowed == false);
    }
}
