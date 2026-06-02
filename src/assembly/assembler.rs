use crate::assembly::assembler_utils::{
    filter_complete_alleles, find_overlapping_alleles, two_haplotypes_are_matching,
};
use crate::util::{DError, DResult};
use log::{debug, trace};
use rand::{prelude::SliceRandom, SeedableRng};

use std::cmp;
use std::collections::{BTreeMap, HashSet};

pub const MIN_EDGE_SUPPORT: i32 = 2;
pub const MIN_ALLELE_SUPPORT: usize = 3;

// Using the terms "haplotypes" and "alleles" interchangeably in this work

/// Graph parameters
#[derive(Clone, Debug)]
pub struct GraphParameters {
    /// expected number of alleles
    pub expect_n_allele: i32,
    /// minimum overlap length
    pub min_overlap: i32,
    /// whether to run twice
    pub run_twice: bool,
    /// whether to perform less filtering
    pub less_filtering: bool,
}

/// Haplotype and its supporting reads
#[derive(Clone, Debug)]
pub struct HaplotypeReads {
    /// unique supporting reads. haplotype -> reads
    pub unique: BTreeMap<Vec<i32>, Vec<Vec<i32>>>,
    /// all supporting reads. read name -> vec of haplotypes that the read possibly support
    pub by_read: BTreeMap<String, Vec<Vec<i32>>>,
}

/// Result of processing complete haplotypes
#[derive(Clone, Debug)]
pub struct ProcessCompleteHaplotypesResult {
    /// supporting reads. allele -> reads
    pub supporting_reads: BTreeMap<Vec<i32>, HashSet<String>>,
    /// alleles that would be removed from complete alleles
    haps_to_remove: HashSet<Vec<i32>>,
    /// all supporting reads. read name -> all alleles it is consistent with
    pub support_by_read: BTreeMap<String, Vec<Vec<i32>>>,
    /// need to perform a second assembly run
    need_for_2nd_run: bool,
}

fn pick_supported_hap_with_seed(
    rng: &mut rand::rngs::SmallRng,
    supported_haps: &[Vec<i32>],
) -> Option<Vec<i32>> {
    if supported_haps.is_empty() {
        return None;
    }
    let mut sorted_haps = supported_haps.to_vec();
    sorted_haps.sort();
    sorted_haps.choose(rng).cloned()
}

/// Result of assembly
#[derive(Clone, Debug)]
pub struct AssemblyResult {
    /// complete alleles
    pub complete: Vec<Vec<i32>>,
    /// incomplete alleles
    pub incomplete: Vec<Vec<i32>>,
    /// cis duplications
    pub special_incomplete: Vec<Vec<i32>>,
    /// supporting reads. allele -> reads
    pub supporting_reads: BTreeMap<Vec<i32>, HashSet<String>>,
    /// nonunique reads
    pub nonunique_reads: Vec<String>,
}

/// Fingerprint graph
#[derive(Clone, Debug)]
pub struct FpGraph {
    /// read name to nodes
    pub reads: BTreeMap<String, Vec<i32>>,
    /// clean version of reads after removing low-support edges
    pub clean_reads: BTreeMap<String, Vec<i32>>,
    /// nodes in the graph
    pub nodes: HashSet<i32>,
    /// (node1, node2) to reads supporting the edges
    pub edges: BTreeMap<(i32, i32), Vec<String>>,
    /// clean version of edges after removing low-support edges
    pub clean_edges: BTreeMap<(i32, i32), usize>,
    /// cyclic nodes and a vector of the number of times they are seen in each read
    pub cyclic_nodes: BTreeMap<i32, Vec<usize>>,
    /// node to its connecting next nodes
    pub next_per_node: BTreeMap<i32, Vec<i32>>,
    /// node to its previous connecting nodes
    pub previous_per_node: BTreeMap<i32, Vec<i32>>,
    /// node to its connecting next nodes
    pub next_per_node_raw: BTreeMap<i32, Vec<i32>>,
    /// node to its previous connecting nodes
    pub previous_per_node_raw: BTreeMap<i32, Vec<i32>>,
    /// from edges to read names
    pub back_to_reads: BTreeMap<Vec<i32>, Vec<String>>,
    /// minimum of fps to overlap before extending
    pub min_overlap: i32,
}

/// Build graph from reads represented as vectors of fingerprints
/// # Arguments
/// * `read_edges` - read name -> vec of fingerprints
/// * `min_overlap` - minimum overlap length
/// # Returns
/// * `FpGraph` - fingerprint graph
pub fn build_graph(read_edges: BTreeMap<String, Vec<i32>>, min_overlap: i32) -> FpGraph {
    FpGraph {
        reads: read_edges,
        clean_reads: BTreeMap::new(),
        nodes: HashSet::new(),
        edges: BTreeMap::new(),
        clean_edges: BTreeMap::new(),
        cyclic_nodes: BTreeMap::new(),
        next_per_node: BTreeMap::new(),
        previous_per_node: BTreeMap::new(),
        next_per_node_raw: BTreeMap::new(),
        previous_per_node_raw: BTreeMap::new(),
        back_to_reads: BTreeMap::new(),
        min_overlap,
    }
}

impl FpGraph {
    /// Main function, from node edges to assembled alleles
    /// # Arguments
    /// * `graph_parameters` - graph parameters
    /// # Returns
    /// * `AssemblyResult` - assembly result
    pub fn run(&mut self, graph_parameters: GraphParameters) -> Result<AssemblyResult, DError> {
        let diploid_mode = graph_parameters.expect_n_allele == 2;
        self.get_edges();
        self.analyze_nodes(false);
        self.next_per_node_raw = self.next_per_node.clone();
        self.previous_per_node_raw = self.previous_per_node.clone();
        self.clean_up_edges(!diploid_mode)?;
        for (node1, node2) in self.clean_edges.keys() {
            self.nodes.insert(*node1);
            self.nodes.insert(*node2);
        }
        self.analyze_nodes(true);
        let assembly_result = self.assemble(!diploid_mode)?;
        let mut complete_haps = assembly_result.complete;
        let mut incomplete_haps = assembly_result.incomplete;
        let mut check_complete_haps = self.process_complete_haps(
            complete_haps.clone(),
            None,
            false,
            graph_parameters.less_filtering,
        )?;
        debug!(
            "After first round, removed alleles are {:?} ",
            check_complete_haps.haps_to_remove
        );
        trace!("support_by_read {:?}", check_complete_haps.support_by_read);

        let clean_reads_without_assembled_haps = self.clean_reads.clone();
        if !check_complete_haps.need_for_2nd_run && !graph_parameters.run_twice {
            debug!("No need for second run");
        } else {
            // assemble a second time adding the assembled contigs
            // from the first run into reads
            for hap in &complete_haps {
                let hap_string = hap
                    .iter()
                    .map(|i| i.to_string())
                    .collect::<Vec<_>>()
                    .join(",");
                self.clean_reads.entry(hap_string).or_insert(hap.to_vec());
                debug!("adding complete allele to read");
            }
            for hap in &incomplete_haps {
                let hap_string = hap
                    .iter()
                    .map(|i| i.to_string())
                    .collect::<Vec<_>>()
                    .join(",");
                self.clean_reads.entry(hap_string).or_insert(hap.to_vec());
                debug!("adding incomplete allele to read");
            }
            debug!("new reads with newly assembled contigs incorporated.");
            debug!("second round of assembly");
            let assembly_result = self.assemble(!diploid_mode)?;
            complete_haps = assembly_result.complete;
            incomplete_haps = assembly_result.incomplete;
            check_complete_haps = self.process_complete_haps(
                complete_haps.clone(),
                None,
                false,
                graph_parameters.less_filtering,
            )?;
        }
        complete_haps = check_complete_haps
            .supporting_reads
            .keys()
            .cloned()
            .collect();
        for hap in check_complete_haps.haps_to_remove {
            debug!("removing redundant allele in haps_to_remove {hap:?}");
            incomplete_haps.push(hap);
        }
        let nodes_not_used = self.get_unused_nodes(complete_haps.clone());
        debug!("nodes not used: {nodes_not_used:?}");

        // remove suspicious complete allele with a very similar incomplete allele that's longer
        let complete_haps_overlapping_incomplete =
            complete_overlapping_incomplete(complete_haps.clone(), incomplete_haps.clone());
        for hap in complete_haps_overlapping_incomplete {
            if complete_haps.contains(&hap) {
                let index = complete_haps
                    .iter()
                    .position(|x| *x == hap)
                    .ok_or("item not found")?;
                complete_haps.remove(index);
            }
            if !incomplete_haps.contains(&hap) {
                incomplete_haps.push(hap);
            }
        }
        if graph_parameters.expect_n_allele != 2 {
            let complete_suspicious = filter_complete_alleles(
                &clean_reads_without_assembled_haps,
                &complete_haps,
                &incomplete_haps,
            )?;
            for hap in complete_suspicious {
                if complete_haps.contains(&hap) {
                    let index = complete_haps
                        .iter()
                        .position(|x| *x == hap)
                        .ok_or("item not found")?;
                    complete_haps.remove(index);
                }
                if !incomplete_haps.contains(&hap) {
                    incomplete_haps.push(hap);
                }
            }
        }

        // more filtering
        if graph_parameters.expect_n_allele == 2 && complete_haps.len() > 1 {
            // check among those reads that do not have allele matches
            // check how many of those contain start/end
            let num_start_or_end_reads_with_no_match = self.check_start_or_end_reads_without_match(
                check_complete_haps.support_by_read.clone(),
            );
            if num_start_or_end_reads_with_no_match >= 10 {
                debug!("{num_start_or_end_reads_with_no_match} starting or ending reads have no allele match");
                for hap in complete_haps {
                    incomplete_haps.push(hap);
                }
                complete_haps = Vec::new();
            }
            // also check nodes that are not used
            else if nodes_not_used.len() > 3 && !graph_parameters.less_filtering {
                debug!("Too many nodes are not used. Removing complete alleles.");
                for hap in complete_haps {
                    incomplete_haps.push(hap);
                }
                complete_haps = Vec::new();
            }
        }

        for hap in &complete_haps {
            let hap_cn = hap.len() - 2;
            debug!("complete haplotype {hap:?} copy number {hap_cn:?}");
        }

        check_complete_haps = self.process_complete_haps(
            complete_haps.clone(),
            None,
            false,
            graph_parameters.less_filtering,
        )?;
        let mut rng = rand::rngs::SmallRng::seed_from_u64(1);
        let mut nonunique_reads: Vec<String> = Vec::new();
        let mut supporting_including_nonuniq = check_complete_haps.supporting_reads.clone();
        for (read, supported_haps) in check_complete_haps.support_by_read.iter() {
            if supported_haps.len() > 1 {
                nonunique_reads.push(read.to_string());
                let support_hap_picked = pick_supported_hap_with_seed(&mut rng, supported_haps)
                    .ok_or("error with seeded choice")?;
                trace!(
                    "nonunique read {:?} haps {:?} picked {:?}",
                    read.to_string(),
                    supported_haps,
                    support_hap_picked
                );
                supporting_including_nonuniq
                    .entry(support_hap_picked)
                    .or_default()
                    .insert(read.to_string());
            }
        }

        Ok(AssemblyResult {
            complete: complete_haps,
            incomplete: incomplete_haps,
            special_incomplete: assembly_result.special_incomplete,
            supporting_reads: supporting_including_nonuniq,
            nonunique_reads,
        })
    }

    /// Check how many of those reads that do not have allele matches contain start/end
    /// # Arguments
    /// * `support_by_read` - read name -> all supporting haps
    /// # Returns
    /// * `usize` - number of reads that do not have allele matches and contain start/end
    pub fn check_start_or_end_reads_without_match(
        &self,
        support_by_read: BTreeMap<String, Vec<Vec<i32>>>,
    ) -> usize {
        let mut start_or_end_reads_with_no_match = Vec::new();
        for (read, read_match) in support_by_read.iter() {
            if read_match.is_empty() {
                let read_info = self.clean_reads.get(read);
                if let Some(read_nodes) = read_info {
                    if has_start(read_nodes) || has_end(read_nodes) {
                        let internal_copies =
                            read_nodes.iter().filter(|x| **x > 0).collect::<Vec<_>>();
                        if !internal_copies.is_empty() {
                            start_or_end_reads_with_no_match.push(read);
                        }
                    }
                }
            }
        }
        start_or_end_reads_with_no_match.len()
    }

    /// Get nodes that are not used given the a set of haplotypes
    /// # Arguments
    /// * `haps` - haplotypes as vectors of nodes
    /// # Returns
    /// * `Vec<i32>` - nodes that are not used
    pub fn get_unused_nodes(&self, haps: Vec<Vec<i32>>) -> Vec<i32> {
        let mut nodes_not_used = Vec::new();
        for node in &self.nodes {
            let mut node_used = false;
            for hap in &haps {
                if hap.contains(node) {
                    node_used = true;
                }
            }
            if !node_used {
                nodes_not_used.push(*node);
            }
        }
        nodes_not_used
    }

    /// Evaluate assembled haplotypes and perform filtering
    /// # Arguments
    /// * `haps_to_assess` - haplotypes to assess
    /// * `min_read_support` - minimum number of supporting reads
    /// * `get_supporting_reads_only` - whether to get only supporting reads
    /// * `less_filtering` - whether to perform less filtering
    /// # Returns
    /// * `ProcessCompleteHaplotypesResult` - result of processing complete haplotypes
    pub fn process_complete_haps(
        &self,
        haps_to_assess: Vec<Vec<i32>>,
        min_read_support: Option<usize>,
        get_supporting_reads_only: bool,
        less_filtering: bool,
    ) -> Result<ProcessCompleteHaplotypesResult, DError> {
        let min_support = min_read_support.unwrap_or(MIN_ALLELE_SUPPORT);
        let mut haps_to_remove = HashSet::new();
        let mut need_for_2nd_run = true;
        let mut remaining_haps = Vec::new();
        if !get_supporting_reads_only {
            let nodes_not_used = self.get_unused_nodes(haps_to_assess.clone());
            // find overlapping haplotypes
            let overlapping_haps = if less_filtering {
                vec![]
            } else {
                find_overlapping_alleles(haps_to_assess.clone(), None)?.0
            };
            for hap in &haps_to_assess {
                if !overlapping_haps.contains(hap) {
                    remaining_haps.push(hap.to_vec());
                }
            }

            if overlapping_haps.is_empty() {
                need_for_2nd_run = false;
            } else {
                if haps_to_assess.len() == 2 && nodes_not_used.len() <= 3 {
                } else if remaining_haps.len() == 1 {
                    let (picked_hap, has_support) =
                        self.pick_from_candidates(overlapping_haps.clone())?;
                    if has_support {
                        need_for_2nd_run = false;
                        for hap in overlapping_haps {
                            if hap != picked_hap {
                                haps_to_remove.insert(hap.to_vec());
                            }
                        }
                    } else {
                        haps_to_remove = overlapping_haps.iter().cloned().collect();
                    }
                } else {
                    haps_to_remove = overlapping_haps.iter().cloned().collect();
                }
                for hap in &haps_to_remove {
                    debug!("haps_to_remove: overlapping haplotype {hap:?}");
                }
            }
        }
        // get supporting reads
        remaining_haps = Vec::new();
        for hap in haps_to_assess {
            if !haps_to_remove.contains(&hap) {
                remaining_haps.push(hap);
            }
        }
        let read_support =
            match_reads_and_haplotypes(self.reads.clone(), remaining_haps, None, false);
        let good_reads = read_support.unique;
        let mut final_supporting_reads = BTreeMap::<Vec<i32>, HashSet<String>>::new();
        for (test_hap, test_hap_support) in good_reads.iter() {
            let mut test_hap_support_filtered = Vec::new();
            for read in test_hap_support {
                let read_num_known_nodes = read.iter().filter(|x| **x != 0).count();
                if get_supporting_reads_only || read_num_known_nodes > 1 {
                    test_hap_support_filtered.push(read.to_vec());
                }
            }
            let nread = test_hap_support_filtered.len();
            debug!("{test_hap:?} has {nread} supporting reads.");
            if nread >= min_support {
                for partial_hap in &test_hap_support_filtered {
                    for read in self.back_to_reads.get(partial_hap).ok_or("key not found")? {
                        final_supporting_reads
                            .entry(test_hap.to_vec())
                            .or_default()
                            .insert(read.to_string());
                    }
                }
            } else {
                debug!("removing low support haplotype {test_hap:?} with {nread} reads.");
            }
        }
        Ok(ProcessCompleteHaplotypesResult {
            supporting_reads: final_supporting_reads,
            haps_to_remove,
            support_by_read: read_support.by_read,
            need_for_2nd_run,
        })
    }

    /// Given a few candidates, pick the best one based on read support or length
    /// # Arguments
    /// * `candidates` - candidates
    /// # Returns
    /// * `(Vec<i32>, bool)` - best candidate and whether it has support
    pub fn pick_from_candidates(
        &self,
        candidates: Vec<Vec<i32>>,
    ) -> Result<(Vec<i32>, bool), DError> {
        let candidates_read_support =
            match_reads_and_haplotypes(self.reads.clone(), candidates.clone(), None, false);
        let mut has_support = false;
        let candidates_read_support_uniq = candidates_read_support.unique;
        if !candidates_read_support_uniq.is_empty() {
            // sort by number of supporting reads
            let mut hap_reads: Vec<(Vec<i32>, usize)> = candidates_read_support_uniq
                .iter()
                .map(|(x, y)| (x.to_vec(), y.len()))
                .collect::<Vec<(Vec<i32>, usize)>>();
            // reverse sort
            hap_reads.sort_by(|a, b| b.1.cmp(&a.1));
            let (most_supported_hap, most_supported_nread) =
                hap_reads.first().ok_or("first not found")?;
            if *most_supported_nread >= MIN_ALLELE_SUPPORT {
                has_support = true;
                return Ok((most_supported_hap.to_vec(), has_support));
            }
        }
        // sort by haplotype length
        let mut candidates_sort = candidates.clone();
        candidates_sort.sort_by_key(|a| a.len());
        // pick the shortest one
        let shortest = candidates_sort.first().ok_or("first not found")?;
        Ok((shortest.to_vec(), has_support))
    }

    /// Get edges from reads
    pub fn get_edges(&mut self) {
        for (read, read_nodes) in self.reads.iter() {
            self.back_to_reads
                .entry(read_nodes.clone())
                .or_default()
                .push(read.to_string());
            if read_nodes.len() > 1 {
                let read_nodes_clone = read_nodes.clone();
                let adjacent_pairs = read_nodes.iter().zip(read_nodes_clone.iter().skip(1));
                for (node1, node2) in adjacent_pairs {
                    // trace!("read {:?}, {}, {}", read, node1, node2);
                    // exclude unknown
                    if *node1 != 0 && *node2 != 0 {
                        self.edges
                            .entry((*node1, *node2))
                            .or_default()
                            .push(read.to_string());
                        // record cyclic nodes
                        if *node1 == *node2 {
                            let node_count = read_nodes.iter().filter(|&n| *n == *node1).count();
                            self.cyclic_nodes
                                .entry(*node1)
                                .or_default()
                                .push(node_count);
                        }
                    }
                }
            }
        }
    }

    /// Get next and previous nodes of all nodes
    /// # Arguments
    /// * `use_clean` - use clean edges or not
    pub fn analyze_nodes(&mut self, use_clean: bool) {
        self.next_per_node = BTreeMap::new();
        self.previous_per_node = BTreeMap::new();
        if use_clean {
            for (node1, node2) in self.clean_edges.keys() {
                self.next_per_node.entry(*node1).or_default().push(*node2);
                self.previous_per_node
                    .entry(*node2)
                    .or_default()
                    .push(*node1);
            }
        } else {
            for (node1, node2) in self.edges.keys() {
                self.next_per_node.entry(*node1).or_default().push(*node2);
                self.previous_per_node
                    .entry(*node2)
                    .or_default()
                    .push(*node1);
            }
        }

        for (node, nodes_next) in self.next_per_node.iter() {
            debug!("node {:?} next nodes: {:?}", node, nodes_next);
        }
        for (node, nodes_prev) in self.previous_per_node.iter() {
            debug!("node {:?} previous nodes: {:?}", node, nodes_prev);
        }
    }

    /// Remove low support edges
    /// # Arguments
    /// * `allow_low_support_edges_to_ends` - whether to allow low support edges to an ending node (this is for D4Z4)
    pub fn clean_up_edges(&mut self, allow_low_support_edges_to_ends: bool) -> DResult {
        let mut removed_edges = Vec::new();
        for ((node1, node2), edge_reads) in self.edges.iter() {
            let nread = edge_reads.len();
            let mut remove_edge = false;
            if nread < MIN_EDGE_SUPPORT as usize
                && self.next_per_node.contains_key(node1)
                && self.previous_per_node.contains_key(node2)
            {
                let nodes_next = self.next_per_node.get(node1).ok_or("key not found")?;
                let nodes_prev = self.previous_per_node.get(node2).ok_or("key not found")?;
                if nodes_next.len() > 1 && nodes_prev.len() > 1 {
                    let mut this_vec = vec![*node1, *node2];
                    this_vec.sort();
                    let mut nodes_next2 = nodes_next.to_vec();
                    nodes_next2.sort();
                    let mut nodes_prev2 = nodes_prev.to_vec();
                    nodes_prev2.sort();
                    if nodes_next2 != this_vec && nodes_prev2 != this_vec {
                        // keep edges that end with -10
                        if allow_low_support_edges_to_ends && *node2 == -10 && nodes_prev.len() <= 3
                        {
                            debug!("keep edge {:?} to {:?}, with {} reads", node1, node2, nread);
                        } else {
                            removed_edges.push((node1, node2));
                            remove_edge = true;
                            debug!(
                                "remove low support edge {:?} to {:?}, with {} reads",
                                node1, node2, nread
                            );
                        }
                    }
                }
            }
            if !remove_edge {
                self.clean_edges.entry((*node1, *node2)).or_insert(nread);
            }
        }
        for (read, read_nodes) in self.reads.iter() {
            let read_nodes_clone = read_nodes.clone();
            let adjacent_pairs = read_nodes.iter().zip(read_nodes_clone.iter().skip(1));
            let mut found_removed_edge = false;
            for (node1, node2) in adjacent_pairs {
                if removed_edges.contains(&(node1, node2)) {
                    found_removed_edge = true;
                }
            }
            if !found_removed_edge {
                self.clean_reads
                    .entry(read.to_string())
                    .or_insert(read_nodes.to_vec());
            }
        }
        Ok(())
    }

    /// Main assembler function: first forward, then backward
    /// # Arguments
    /// * `strict_for_cyclic_nodes` - whether to be strict when extending cyclic nodes
    /// # Returns
    /// * `AssemblyResult` - assembly result
    pub fn assemble(&mut self, strict_for_cyclic_nodes: bool) -> Result<AssemblyResult, DError> {
        let mut assembled_haps: Vec<Vec<i32>> = Vec::new();
        // forward
        let mut complete_haps_forward: Vec<Vec<i32>> = Vec::new();
        let mut incomplete_haps_forward: Vec<Vec<i32>> = Vec::new();
        let starts = self
            .next_per_node
            .keys()
            .filter(|x| **x < 0 && **x > -10)
            .collect::<Vec<_>>();
        for starting_node in starts {
            let starting_next_nodes = self
                .next_per_node
                .get(starting_node)
                .ok_or("key not found in next_per_node")?;
            for this_node in starting_next_nodes {
                assembled_haps.push(vec![*starting_node, *this_node]);
            }
            let mut nstep = 0;
            loop {
                nstep += 1;
                debug!(
                    "Step {} going forward: start with {:?}",
                    nstep, assembled_haps
                );
                let mut all_extended: Vec<Vec<i32>> = Vec::new();
                for this_hap in assembled_haps {
                    debug!("extend node {this_hap:?}");
                    let extended_haps =
                        self.assemble_next(this_hap.clone(), strict_for_cyclic_nodes)?;
                    debug!("extended {:?} to {:?}", this_hap, extended_haps);
                    if extended_haps.is_empty() {
                        if !incomplete_haps_forward.contains(&this_hap) {
                            incomplete_haps_forward.push(this_hap);
                        }
                    } else {
                        for each_extended_hap in &extended_haps {
                            let first_unit = each_extended_hap.first().ok_or("first not found")?;
                            let last_unit = each_extended_hap.last().ok_or("last not found")?;
                            if *first_unit > -10
                                && *first_unit < 0
                                && *last_unit <= -10
                                && !complete_haps_forward.contains(each_extended_hap)
                            {
                                if extended_haps.len() == 1 {
                                    complete_haps_forward.push(each_extended_hap.to_vec());
                                } else {
                                    incomplete_haps_forward.push(each_extended_hap.to_vec());
                                }
                            } else {
                                all_extended.push(each_extended_hap.to_vec());
                            }
                        }
                    }
                }
                assembled_haps = all_extended
                    .iter()
                    .filter(|x| {
                        !complete_haps_forward.contains(x) && !incomplete_haps_forward.contains(x)
                    })
                    .map(|x| x.to_vec())
                    .collect::<Vec<_>>();
                if assembled_haps.is_empty() {
                    break;
                }
                if nstep > 100 || assembled_haps.len() > 50 {
                    for a in &assembled_haps {
                        incomplete_haps_forward.push(a.to_vec());
                    }
                    debug!("stopping assembly at step {nstep}");
                    break;
                }
            }
        }
        debug!("forward incomplete: {:?}", incomplete_haps_forward);
        debug!("forward complete: {:?}", complete_haps_forward);

        // backward
        let mut assembled_haps: Vec<Vec<i32>> = Vec::new();
        let mut complete_haps_backward: Vec<Vec<i32>> = Vec::new();
        let mut incomplete_haps_backward: Vec<Vec<i32>> = Vec::new();
        let mut special_incomplete: Vec<Vec<i32>> = Vec::new();
        let ends = self
            .previous_per_node
            .keys()
            .filter(|x| **x <= -10)
            .collect::<Vec<_>>();
        for ending_node in ends {
            let ending_next_nodes = self
                .previous_per_node
                .get(ending_node)
                .ok_or("key not found in previous_per_node")?;
            for this_node in ending_next_nodes {
                assembled_haps.push(vec![*this_node, *ending_node]);
            }
            let mut nstep = 0;
            loop {
                nstep += 1;
                debug!(
                    "Step {} going backward: start with {:?}",
                    nstep, assembled_haps
                );
                let mut all_extended: Vec<Vec<i32>> = Vec::new();
                for this_hap in assembled_haps {
                    debug!("extend node {this_hap:?}");
                    let extended_haps =
                        self.assemble_prev(this_hap.clone(), strict_for_cyclic_nodes)?;
                    debug!("extended {:?} to {:?}", this_hap, extended_haps);
                    if extended_haps.is_empty() {
                        if !incomplete_haps_backward.contains(&this_hap) {
                            incomplete_haps_backward.push(this_hap);
                        }
                    } else {
                        for each_extended_hap in &extended_haps {
                            let first_unit = each_extended_hap.first().ok_or("first not found")?;
                            let last_unit = each_extended_hap.last().ok_or("last not found")?;
                            if *first_unit > -10
                                && *first_unit < 0
                                && *last_unit <= -10
                                && !complete_haps_backward.contains(each_extended_hap)
                            {
                                if extended_haps.len() == 1 {
                                    complete_haps_backward.push(each_extended_hap.to_vec());
                                } else {
                                    incomplete_haps_backward.push(each_extended_hap.to_vec());
                                }
                            } else {
                                if *first_unit == -10 {
                                    let mut each_extended_hap_remove_first_end =
                                        each_extended_hap.clone();
                                    each_extended_hap_remove_first_end.remove(0);
                                    if !incomplete_haps_backward
                                        .contains(&each_extended_hap_remove_first_end)
                                    {
                                        debug!(
                                            "add cis-dup hap {each_extended_hap_remove_first_end:?} to incomplete"
                                        );
                                        incomplete_haps_backward.push(each_extended_hap.clone());
                                        special_incomplete.push(each_extended_hap.clone());
                                    }
                                } else {
                                    all_extended.push(each_extended_hap.to_vec());
                                }
                            }
                        }
                    }
                }
                assembled_haps = all_extended
                    .iter()
                    .filter(|x| {
                        !complete_haps_backward.contains(x) && !incomplete_haps_backward.contains(x)
                    })
                    .map(|x| x.to_vec())
                    .collect::<Vec<_>>();
                if assembled_haps.is_empty() {
                    break;
                }
                if nstep > 100 || assembled_haps.len() > 50 {
                    for a in &assembled_haps {
                        incomplete_haps_backward.push(a.to_vec());
                    }
                    debug!("stopping assembly at step {nstep}");
                    break;
                }
            }
        }
        debug!("backward incomplete: {:?}", incomplete_haps_backward);
        debug!("backward complete: {:?}", complete_haps_backward);

        let s1: HashSet<Vec<i32>> = complete_haps_forward.iter().cloned().collect();
        let s2: HashSet<Vec<i32>> = complete_haps_backward.iter().cloned().collect();
        let complete_haps = s1.union(&s2).cloned().collect::<Vec<_>>();

        let s3: HashSet<Vec<i32>> = incomplete_haps_forward.iter().cloned().collect();
        let s4: HashSet<Vec<i32>> = incomplete_haps_backward.iter().cloned().collect();
        let mut incomplete_haps = Vec::new();
        for hap in s3.union(&s4) {
            if has_start(hap) || has_end(hap) {
                let mut found_in_complete: bool = false;
                for complete_hap in &complete_haps {
                    let is_contained = b_is_contained_in_a(complete_hap, hap);
                    trace!(
                        "comparing incomplete {hap:?} to complete {complete_hap:?} {is_contained}"
                    );
                    if is_contained {
                        found_in_complete = true;
                        break;
                    }
                }
                if !found_in_complete {
                    incomplete_haps.push(hap.to_vec());
                }
            }
        }
        let mut redundant_haps = Vec::new();
        for hap in &incomplete_haps {
            for hap1 in &incomplete_haps {
                if hap != hap1 {
                    let is_contained = b_is_contained_in_a(hap1, hap);
                    if is_contained {
                        redundant_haps.push(hap.to_vec());
                        break;
                    }
                }
            }
        }
        for hap in &redundant_haps {
            let index = incomplete_haps
                .iter()
                .position(|x| *x == *hap)
                .ok_or("item not found")?;
            log::debug!("removing redundant incomplete haplotype {hap:?}");
            incomplete_haps.remove(index);
        }
        let merged_result = self.merge_two_incomplete(complete_haps, incomplete_haps)?;
        Ok(AssemblyResult {
            complete: merged_result.complete,
            incomplete: merged_result.incomplete,
            special_incomplete,
            supporting_reads: BTreeMap::new(),
            nonunique_reads: Vec::new(),
        })
    }

    /// Merge two incomplete haplotypes that each has a start or end.
    /// They are either overlapping or non-overlapping but with long read support
    /// # Arguments
    /// * `complete_haps` - complete alleles
    /// * `incomplete_haps` - incomplete alleles
    /// # Returns
    /// * `AssemblyResult` - assembly result
    pub fn merge_two_incomplete(
        &self,
        complete_haps: Vec<Vec<i32>>,
        incomplete_haps: Vec<Vec<i32>>,
    ) -> Result<AssemblyResult, DError> {
        let mut complete_clone = complete_haps.clone();
        let mut incomplete_clone = incomplete_haps.clone();
        let mut all_haps = complete_haps.clone();
        let mut incomplete_clone2 = incomplete_haps.clone();
        all_haps.append(&mut incomplete_clone2);
        let nodes_not_used = self.get_unused_nodes(all_haps.clone());
        debug!(
            "complete: {}, incomplete: {}, unused nodes: {:?}",
            complete_haps.len(),
            incomplete_haps.len(),
            nodes_not_used
        );
        if incomplete_haps.len() == 2 && complete_haps.len() == 1 && nodes_not_used.is_empty() {
            debug!("trying to merge the only two incomplete haps {incomplete_haps:?}");
            let mut ct1 = incomplete_haps.first().ok_or("first not found")?;
            let mut ct2 = incomplete_haps.last().ok_or("last not found")?;
            if has_end(ct1) && has_start(ct2) {
                ct1 = incomplete_haps.last().ok_or("last not found")?;
                ct2 = incomplete_haps.first().ok_or("first not found")?;
            }
            if has_end(ct2) && has_start(ct1) {
                let mut new_ct_candidates = Vec::new();
                let min_len = cmp::min(ct1.len(), ct2.len());
                // overlapping
                for j in 0..min_len {
                    let part1 = &ct1[(ct1.len() - (min_len - j))..];
                    let part2 = &ct2[..(min_len - j)];
                    if part1 == part2 {
                        let new_ct = [&ct1[..], &ct2[(min_len - j)..]].concat();
                        new_ct_candidates.push(new_ct.to_vec());
                    }
                }
                // nonoverlapping
                let ct1_last = ct1.last().ok_or("last not found")?;
                let ct2_first = ct2.first().ok_or("first not found")?;
                if self.edges.contains_key(&(*ct1_last, *ct2_first)) {
                    // require read support over a larger region (3 nodes)
                    let mut segments_supported: Vec<bool> = vec![false, false];
                    let mut segment1 = ct1[(ct1.len() - 2)..].to_vec();
                    segment1.push(*ct2_first);
                    let mut segment2 = ct2[..2].to_vec();
                    segment2.insert(0, *ct1_last);
                    for (segment_index, segment) in [segment1, segment2].iter().enumerate() {
                        for read_nodes in self.reads.values() {
                            let read_nodes_len = read_nodes.len();
                            if read_nodes_len > 2 {
                                for i in 0..(read_nodes_len - 2) {
                                    let read_segment = &read_nodes[i..(i + 3)];
                                    if read_segment == segment && !read_segment.contains(&0) {
                                        segments_supported[segment_index] = true;
                                    }
                                }
                            }
                        }
                    }
                    if !segments_supported.contains(&false) {
                        let new_ct = [&ct1[..], &ct2[..]].concat();
                        new_ct_candidates.push(new_ct.to_vec());
                    }
                }
                if !new_ct_candidates.is_empty() {
                    let (picked_candidate, _has_support) =
                        self.pick_from_candidates(new_ct_candidates.clone())?;
                    debug!("merge candidates: {new_ct_candidates:?}");
                    debug!("picked candidate: {picked_candidate:?}");
                    complete_clone.push(picked_candidate.to_vec());
                    let index = incomplete_clone
                        .iter()
                        .position(|x| *x == *ct1)
                        .ok_or("item not found")?;
                    incomplete_clone.remove(index);
                    let index = incomplete_clone
                        .iter()
                        .position(|x| *x == *ct2)
                        .ok_or("item not found")?;
                    incomplete_clone.remove(index);
                }
            }
        }
        Ok(AssemblyResult {
            complete: complete_clone,
            incomplete: incomplete_clone,
            special_incomplete: vec![],
            supporting_reads: BTreeMap::new(),
            nonunique_reads: Vec::new(),
        })
    }

    /// Extend given haplotype by one node, forward
    /// # Arguments
    /// * `this_hap` - haplotype
    /// * `strict_for_cyclic_nodes` - whether to be strict when extending cyclic nodes
    /// # Returns
    /// * `Vec<Vec<i32>>` - extended haplotypes
    fn assemble_next(
        &self,
        this_hap: Vec<i32>,
        strict_for_cyclic_nodes: bool,
    ) -> Result<Vec<Vec<i32>>, DError> {
        let last_unit = this_hap.last().ok_or("last not found")?;
        if !self.next_per_node.contains_key(last_unit) {
            return Ok(vec![]);
        }
        // next node cannot be starting node
        let next_nodes = &self
            .next_per_node
            .get(last_unit)
            .ok_or("key not found in next_per_node")?
            .iter()
            .filter(|x| **x >= 0 || **x <= -10)
            .map(|x| *x)
            .collect::<Vec<_>>();
        let mut haps_candidates = Vec::new();
        for next_node in next_nodes {
            let mut extended_hap = this_hap.clone();
            extended_hap.push(*next_node);
            let to_include =
                self.include_cyclic_node_forward(*last_unit, *next_node, extended_hap.clone())?;
            if to_include {
                haps_candidates.push(extended_hap);
            }
        }
        debug!("candidates: {:?}", haps_candidates);
        if haps_candidates.is_empty() {
            return Ok(vec![]);
        }
        if haps_candidates.len() == 1 {
            if self.min_overlap == 2 {
                return Ok(haps_candidates);
            }
            let next_node = next_nodes.first().ok_or("err")?;
            let next_nodes_raw = self.next_per_node_raw.get(last_unit).ok_or("err")?;
            if next_nodes_raw.len() == 1 && *next_node == -10 {
                return Ok(haps_candidates);
            }
            let this_node = last_unit;
            if self.previous_per_node_raw.contains_key(next_node) && next_nodes_raw.len() == 1 {
                let next_node_prev = self.previous_per_node_raw.get(next_node).ok_or("err")?;
                if next_node_prev.len() == 1 && next_node_prev.first().ok_or("err")? == this_node {
                    let this_node_prev = self.previous_per_node_raw.get(this_node);
                    let next_node_next = self.next_per_node_raw.get(next_node);
                    if let Some(this_node_prev1) = this_node_prev {
                        if let Some(next_node_next1) = next_node_next {
                            if next_node_next1.len() == 1 && this_node_prev1.len() == 1 {
                                let mut nodes_set = HashSet::new();
                                nodes_set.insert(this_node);
                                nodes_set.insert(next_node);
                                nodes_set.insert(next_node_next1.first().ok_or("err")?);
                                nodes_set.insert(this_node_prev1.first().ok_or("err")?);
                                if nodes_set.len() == 4 {
                                    return Ok(haps_candidates);
                                }
                            }
                        }
                    }
                }
            }
        }
        let mut check_overlap_length: usize = (self
            .require_overlap_length_check_backward(this_hap.clone(), strict_for_cyclic_nodes)?
            as usize)
            + 2;
        if check_overlap_length < (self.min_overlap as usize) {
            check_overlap_length = self.min_overlap as usize;
        }
        if check_overlap_length > this_hap.len() + 1 {
            check_overlap_length = this_hap.len() + 1;
        }
        debug!("check positions n={check_overlap_length}");
        let read_support = match_reads_and_haplotypes(
            self.reads.clone(),
            haps_candidates.clone(),
            Some(check_overlap_length),
            false,
        );
        debug!("read_support {:?}", read_support.unique);
        // get support length on each read
        let mut support_length: BTreeMap<Vec<i32>, Vec<(Vec<i32>, usize)>> = BTreeMap::new();
        let candidate_len = haps_candidates
            .first()
            .ok_or("first not found in haps_candidates")?
            .len();
        let n_max = cmp::min(10, candidate_len - check_overlap_length + 1);
        for n in 0..n_max {
            let read_support_stringent = match_reads_and_haplotypes(
                self.reads.clone(),
                haps_candidates.clone(),
                Some(check_overlap_length + n),
                false,
            );
            let read_support_stringent_uniq = read_support_stringent.unique;
            if !read_support_stringent_uniq.is_empty() {
                for (test_hap, test_hap_reads) in read_support_stringent_uniq.iter() {
                    for test_hap_read in test_hap_reads {
                        support_length
                            .entry(test_hap.to_vec())
                            .or_default()
                            .push((test_hap_read.to_vec(), check_overlap_length + n));
                    }
                }
            } else {
                break;
            }
        }
        let filtered_candidates = filter_assembled_candidates(read_support.unique, support_length)?;
        Ok(filtered_candidates)
    }

    /// Extend given haplotype by one node, going backward
    /// # Arguments
    /// * `this_hap` - haplotype
    /// * `strict_for_cyclic_nodes` - whether to be strict when extending cyclic nodes
    /// # Returns
    /// * `Vec<Vec<i32>>` - extended haplotypes
    fn assemble_prev(
        &self,
        this_hap: Vec<i32>,
        strict_for_cyclic_nodes: bool,
    ) -> Result<Vec<Vec<i32>>, DError> {
        let first_unit = this_hap.first().ok_or("last not found")?;
        if !self.previous_per_node.contains_key(first_unit) {
            return Ok(vec![]);
        }
        // prev node can be ending node, filter afterwards
        let prev_nodes = &self
            .previous_per_node
            .get(first_unit)
            .ok_or("key not found in previous_per_node")?
            .iter()
            //.filter(|x| **x > -10)
            .map(|x| *x)
            .collect::<Vec<_>>();
        let mut haps_candidates = Vec::new();
        for prev_node in prev_nodes {
            let mut extended_hap = this_hap.clone();
            extended_hap.insert(0, *prev_node);
            let to_include =
                self.include_cyclic_node_backward(*first_unit, *prev_node, extended_hap.clone())?;
            if to_include {
                haps_candidates.push(extended_hap);
            }
        }
        debug!("candidates: {:?}", haps_candidates);
        if haps_candidates.is_empty() {
            return Ok(vec![]);
        }
        // this was not originally implemented below in python, 2 line below
        //if haps_candidates.len() == 1 && self.min_overlap == 2 {
        //    return Ok(haps_candidates);
        //}
        if haps_candidates.len() == 1 {
            if self.min_overlap == 2 {
                return Ok(haps_candidates);
            }
            // further check number of prev/next nodes
            let prev_nodes_raw = self.previous_per_node_raw.get(first_unit).ok_or("err")?;
            let prev_node = prev_nodes.first().ok_or("err")?;
            let this_node = first_unit;
            if self.next_per_node_raw.contains_key(prev_node) && prev_nodes_raw.len() == 1 {
                let prev_node_next = self.next_per_node_raw.get(prev_node).ok_or("err")?;
                if prev_node_next.len() == 1 && prev_node_next.first().ok_or("err")? == this_node {
                    let this_node_next = self.next_per_node_raw.get(this_node);
                    let prev_node_prev = self.previous_per_node_raw.get(prev_node);
                    if let Some(this_node_next1) = this_node_next {
                        if let Some(prev_node_prev1) = prev_node_prev {
                            if this_node_next1.len() == 1 && prev_node_prev1.len() == 1 {
                                let mut nodes_set = HashSet::new();
                                nodes_set.insert(this_node);
                                nodes_set.insert(prev_node);
                                nodes_set.insert(this_node_next1.first().ok_or("err")?);
                                nodes_set.insert(prev_node_prev1.first().ok_or("err")?);
                                if nodes_set.len() == 4 {
                                    return Ok(haps_candidates);
                                }
                            }
                        }
                    }
                }
            }
        }
        let mut check_overlap_length: usize = (self
            .require_overlap_length_check_forward(this_hap.clone(), strict_for_cyclic_nodes)?
            as usize)
            + 2;
        if check_overlap_length < (self.min_overlap as usize) {
            check_overlap_length = self.min_overlap as usize;
        }
        if check_overlap_length > this_hap.len() + 1 {
            check_overlap_length = this_hap.len() + 1;
        }
        debug!("check positions n={check_overlap_length}");
        let read_support = match_reads_and_haplotypes(
            self.reads.clone(),
            haps_candidates.clone(),
            Some(check_overlap_length),
            true,
        );
        debug!("read_support {:?}", read_support.unique);
        // get support length on each read
        let mut support_length: BTreeMap<Vec<i32>, Vec<(Vec<i32>, usize)>> = BTreeMap::new();
        let candidate_len = haps_candidates
            .first()
            .ok_or("first not found in haps_candidates")?
            .len();
        let n_max = cmp::min(10, candidate_len - check_overlap_length + 1);
        for n in 0..n_max {
            let read_support_stringent = match_reads_and_haplotypes(
                self.reads.clone(),
                haps_candidates.clone(),
                Some(check_overlap_length + n),
                true,
            );
            let read_support_stringent_uniq = read_support_stringent.unique;
            if !read_support_stringent_uniq.is_empty() {
                for (test_hap, test_hap_reads) in read_support_stringent_uniq.iter() {
                    for test_hap_read in test_hap_reads {
                        support_length
                            .entry(test_hap.to_vec())
                            .or_default()
                            .push((test_hap_read.to_vec(), check_overlap_length + n));
                    }
                }
            } else {
                break;
            }
        }
        let filtered_candidates = filter_assembled_candidates(read_support.unique, support_length)?;
        Ok(filtered_candidates)
    }

    /// Get the needed overlap length backward in order to assemble forward
    /// # Arguments
    /// * `this_hap` - haplotype to extend
    /// * `strict_for_cyclic_nodes` - whether to be strict when extending cyclic nodes
    /// # Returns
    /// * `i32` - needed overlap length
    pub fn require_overlap_length_check_backward(
        &self,
        this_hap: Vec<i32>,
        strict_for_cyclic_nodes: bool,
    ) -> Result<i32, DError> {
        let hap_len = this_hap.len();
        // if cyclic, ovl_len should be at least the number of cylic nodes
        let last_node = this_hap.last().ok_or("last node not found")?;
        let mut repeat_count = 0;
        for i in 1..(hap_len + 1) {
            if this_hap[hap_len - i] == *last_node {
                repeat_count += 1;
            } else {
                break;
            }
        }
        if repeat_count <= 2 {
            repeat_count = 0;
        }
        debug!("this hap ends with a repeat of {last_node} repeat_count {repeat_count}");

        let mut multi_mapping_segment = false;
        let mut i_index = 0;
        for i in 0..hap_len {
            i_index = i;
            let test_partial_hap = &this_hap[i..];
            let mut test_partial_hap_previous = HashSet::new();
            for read in self.clean_reads.values() {
                let read_node_len = read.len();
                for j in 0..read_node_len {
                    let j_end = j + test_partial_hap.len();
                    if j > 0 && j_end <= read_node_len && *test_partial_hap == read[j..j_end] {
                        let read_previous_node =
                            read.get(j - 1).ok_or("index not found in read")?;
                        if *read_previous_node != 0 {
                            test_partial_hap_previous.insert(*read_previous_node);
                        }
                    }
                }
            }
            if test_partial_hap_previous.len() > 1 {
                trace!("Found multi mapping segment (backward) {test_partial_hap_previous:?} at i_index {i_index}");
                multi_mapping_segment = true;
                break;
            }
        }
        let ovl_len: i32 = if !multi_mapping_segment {
            0
        } else if strict_for_cyclic_nodes {
            ((hap_len - i_index) as i32).max(repeat_count as i32)
        } else {
            (hap_len - i_index) as i32
        };
        debug!("require ovl_len {ovl_len}");
        Ok(ovl_len)
    }

    /// Get the needed overlap length forward in order to assemble backward
    /// # Arguments
    /// * `this_hap` - haplotype to extend
    /// * `strict_for_cyclic_nodes` - whether to be strict when extending cyclic nodes
    /// # Returns
    /// * `i32` - needed overlap length
    pub fn require_overlap_length_check_forward(
        &self,
        this_hap: Vec<i32>,
        strict_for_cyclic_nodes: bool,
    ) -> Result<i32, DError> {
        let hap_len = this_hap.len();
        // if cyclic, ovl_len should be at least the number of cylic nodes
        let first_node = this_hap.first().ok_or("first node not found")?;
        let mut repeat_count = 0;
        for i in 0..hap_len {
            if this_hap[i] == *first_node {
                repeat_count += 1;
            } else {
                break;
            }
        }
        if repeat_count <= 2 {
            repeat_count = 0;
        }
        debug!("this hap starts with a repeat of {first_node} repeat_count {repeat_count}");

        let mut multi_mapping_segment = false;
        let mut i_index = 0;
        for i in 0..hap_len {
            i_index = i;
            let test_partial_hap = &this_hap[..(hap_len - i)];
            let test_partial_hap_len = test_partial_hap.len();
            let mut test_partial_hap_next = HashSet::new();
            for read in self.clean_reads.values() {
                let read_node_len = read.len();
                for j in 0..read_node_len {
                    if j >= test_partial_hap_len && j <= read_node_len {
                        let read_nodes_equivalent = &read[(j - test_partial_hap_len)..j];
                        if *test_partial_hap == *read_nodes_equivalent {
                            let read_next_node = read.get(j).ok_or("index not found in read")?;
                            if *read_next_node != 0 {
                                test_partial_hap_next.insert(*read_next_node);
                            }
                        }
                    }
                }
            }
            if test_partial_hap_next.len() > 1 {
                trace!(
                    "Found multi mapping segment (forward) {test_partial_hap_next:?} at i_index {i_index}"
                );
                multi_mapping_segment = true;
                break;
            }
        }
        let ovl_len: i32 = if !multi_mapping_segment {
            0
        } else if strict_for_cyclic_nodes {
            ((hap_len - i_index) as i32).max(repeat_count as i32)
        } else {
            (hap_len - i_index) as i32
        };
        debug!("require ovl_len {ovl_len}");
        Ok(ovl_len)
    }

    /// Stop extending cyclic nodes when its count is larger than longest ever seen in reads
    /// # Arguments
    /// * `this_node` - ending of current hap
    /// * `next_node` - candidate node for extension
    /// * `hap_nodes` - the extended haplotype
    /// # Returns
    /// * `bool` - true if to include, false otherwise
    pub fn include_cyclic_node_forward(
        &self,
        this_node: i32,
        next_node: i32,
        hap_nodes: Vec<i32>,
    ) -> Result<bool, DError> {
        let mut to_include = true;
        if this_node == next_node {
            let seen_cn = self.cyclic_nodes.get(&this_node).ok_or("key not found")?;
            let highest_cn: &usize = seen_cn.iter().max().ok_or("max not found")?;
            let mut cyclic_node_cn: usize = 0;
            for (i, node) in hap_nodes.iter().rev().enumerate() {
                cyclic_node_cn = i;
                if *node != this_node {
                    break;
                }
            }
            trace!(
                "hap {hap_nodes:?} cyclic_node_cn {cyclic_node_cn} highest_cn {:?}",
                *highest_cn
            );
            if cyclic_node_cn > *highest_cn {
                to_include = false
            }
        }
        Ok(to_include)
    }

    /// Stop extending cyclic nodes when its count is larger than longest ever seen in reads
    /// # Arguments
    /// * `this_node` - ending of current hap
    /// * `prev_node` - candidate node for extension
    /// * `hap_nodes` - the extended haplotype
    /// # Returns
    /// * `bool` - true if to include, false otherwise
    pub fn include_cyclic_node_backward(
        &self,
        this_node: i32,
        prev_node: i32,
        hap_nodes: Vec<i32>,
    ) -> Result<bool, DError> {
        let mut to_include = true;
        if this_node == prev_node {
            let seen_cn = self.cyclic_nodes.get(&this_node).ok_or("key not found")?;
            let highest_cn: &usize = seen_cn.iter().max().ok_or("max not found")?;
            let mut cyclic_node_cn: usize = 0;
            for (i, node) in hap_nodes.iter().enumerate() {
                cyclic_node_cn = i;
                if *node != this_node {
                    break;
                }
            }
            trace!(
                "hap {hap_nodes:?} cyclic_node_cn {cyclic_node_cn} highest_cn {:?}",
                *highest_cn
            );
            // originally in python: if cyclic_node_cn+1 > *highest_cn
            if cyclic_node_cn > *highest_cn {
                to_include = false
            }
        }
        Ok(to_include)
    }
}

/// Check if a haplotype contains the starting flank
/// # Arguments
/// * `hap` - a vector representing the haplotype
/// # Returns
/// * `bool` - true if contains the starting flank, false otherwise
/// # Examples
///  ```rust
/// use kivvi::assembly::assembler::has_start;
/// let test_hap = vec![1, 2, 3];
/// let test_hap_has_start = has_start(&test_hap);
/// assert!(!test_hap_has_start);
/// let test_hap = vec![-1, 2, 3];
/// let test_hap_has_start = has_start(&test_hap);
/// assert!(test_hap_has_start);
/// ```
pub fn has_start(hap: &Vec<i32>) -> bool {
    for node in hap {
        if *node > -10 && *node < 0 {
            return true;
        }
    }
    false
}

/// Check if a haplotype contains the ending flank
/// # Arguments
/// * `hap` - a vector representing the haplotype
/// # Returns
/// * `bool` - true if contains the ending flank, false otherwise
/// # Examples
///  ```rust
/// use kivvi::assembly::assembler::has_end;
/// let test_hap = vec![1, 2, 3];
/// let test_hap_has_end = has_end(&test_hap);
/// assert!(!test_hap_has_end);
/// let test_hap = vec![2, 3, -10];
/// let test_hap_has_end = has_end(&test_hap);
/// assert!(test_hap_has_end);
/// ```
pub fn has_end(hap: &Vec<i32>) -> bool {
    for node in hap {
        if *node <= -10 {
            return true;
        }
    }
    false
}

/// Check if array b is contained in array a
/// # Arguments
/// * `a` - array a
/// * `b` - array b
/// # Returns
/// * `bool` - true if b is contained in a, false otherwise
pub fn b_is_contained_in_a(a: &[i32], b: &[i32]) -> bool {
    let a_len = a.len();
    let b_len = b.len();
    if b == a {
        return true;
    }
    if b_len > a_len {
        return false;
    };
    if b_len == 0 {
        return true;
    }
    if a_len == 0 {
        return false;
    }

    'outer: for i in 0..(a_len - b_len + 1) {
        for j in 0..b_len {
            if a[i + j] != b[j] {
                continue 'outer;
            }
        }
        return true;
    }
    false
}

/// Pick from assembled candidates based on number of supporting reads and overlapping length of supporting reads.
/// Compare number of reads, how many are longer in support length, also the diff of max support
/// # Arguments
/// * `unique_support` - haplotype -> uniquely supporting read haplotypes
/// * `support_length` - haplotype -> vec(supporting read haplotype, overlap length tried)
/// # Returns
/// * `Vec<Vec<i32>>` - filtered candidates
pub fn filter_assembled_candidates(
    unique_support: BTreeMap<Vec<i32>, Vec<Vec<i32>>>,
    support_length: BTreeMap<Vec<i32>, Vec<(Vec<i32>, usize)>>,
) -> Result<Vec<Vec<i32>>, DError> {
    // filtered candidates to return
    let mut filtered_candidates = Vec::new();
    // haplotype -> list of support length (each value is a read)
    let mut support_length_per_hap = BTreeMap::new();
    for (test_hap, hap_reads) in unique_support.iter() {
        // read -> all support values
        let mut all_support_per_read: BTreeMap<Vec<i32>, Vec<usize>> = BTreeMap::new();
        let hap_support_length = support_length
            .get(test_hap)
            .ok_or("key not found in support_length")?;

        let mut nread_stringent = Vec::new();

        for (read, sup_len) in hap_support_length {
            all_support_per_read
                .entry(read.to_vec())
                .or_default()
                .push(*sup_len);
        }
        for read in hap_reads {
            let read_support_all = all_support_per_read
                .get(read)
                .ok_or("key not found in all_support_per_read")?;
            let read_support_max = read_support_all
                .iter()
                .max()
                .ok_or("cannot get max of read_support_all")?;
            nread_stringent.push(*read_support_max);
        }

        support_length_per_hap
            .entry(test_hap)
            .or_insert(nread_stringent.clone());

        let nread = hap_reads.len();
        debug!("haplotype {test_hap:?} has {nread} supports and checking support length {nread_stringent:?}");

        if nread > 0 {
            filtered_candidates.push(test_hap.to_vec());
        }
    }
    // more filtering based on support length
    if filtered_candidates.len() > 1 && !support_length.is_empty() {
        // values of support_length_per_hap
        let mut all_support_length = Vec::new();
        let mut all_nread = Vec::new();
        for (_hap, lens) in support_length_per_hap.iter() {
            for len in lens {
                all_support_length.push(*len);
                all_nread.push(lens.len());
            }
        }
        let max_support_len = all_support_length
            .iter()
            .max()
            .ok_or("max not found in all_support_length")?;
        let mut longest_support_haps = Vec::new();
        for (_hap, lens) in support_length_per_hap.iter() {
            let lens_max = lens.iter().max().ok_or("max not found in lens")?;
            if *lens_max == *max_support_len {
                longest_support_haps.push(lens.to_vec());
            }
        }
        longest_support_haps.sort_by_key(|b| std::cmp::Reverse(b.iter().sum::<usize>()));
        let longest_support = longest_support_haps
            .first()
            .ok_or("first not found in longest_support_haps")?;
        let longest_support_sum = longest_support.iter().sum::<usize>();
        let longest_support_max = longest_support
            .iter()
            .max()
            .ok_or("max not found in longest_support")?;
        let max_nread = all_nread.iter().max().ok_or("max not found in all_nread")?;
        // one read vs 10 reads or more
        for test_hap in unique_support.keys() {
            let this_hap_support = support_length_per_hap
                .get(test_hap)
                .ok_or("key not found in support_length_per_hap")?;
            let this_hap_support_sum = this_hap_support.iter().sum::<usize>();
            if this_hap_support.len() == 1
                && *max_nread >= 10
                && this_hap_support_sum as i32 <= longest_support_sum as i32 - 20
            {
                if filtered_candidates.contains(test_hap) {
                    let index = filtered_candidates
                        .iter()
                        .position(|x| *x == *test_hap)
                        .ok_or("item not found")?;
                    filtered_candidates.remove(index);
                }
            }
        }
        // longest support appear at least twice, and a diff of one in longest support
        // and low support for the one to be filtered
        // this scenario is probably better covered in the third one next?
        let longest_support_count = longest_support
            .iter()
            .filter(|x| *x == max_support_len)
            .collect::<Vec<_>>()
            .len();
        if longest_support_count >= 2 {
            for (test_hap, hap_reads) in unique_support.iter() {
                let nread = hap_reads.len();
                let this_hap_support = support_length_per_hap
                    .get(test_hap)
                    .ok_or("key not found in support_length_per_hap")?;
                let this_hap_support_max = this_hap_support
                    .iter()
                    .max()
                    .ok_or("max not found in this_hap_support")?;
                // originally in python: no (max_nread > 2)
                if nread <= 2 && *this_hap_support_max < *longest_support_max && *max_nread > 2 {
                    if filtered_candidates.contains(test_hap) {
                        let index = filtered_candidates
                            .iter()
                            .position(|x| *x == *test_hap)
                            .ok_or("item not found")?;
                        filtered_candidates.remove(index);
                    }
                }
            }
        }
        // cases with support of one or two reads, require a big difference in support length
        // originally in python: else if
        if *max_nread == 1 {
            for test_hap in unique_support.keys() {
                let this_hap_support = support_length_per_hap
                    .get(test_hap)
                    .ok_or("key not found in support_length_per_hap")?;
                let this_hap_support_max = this_hap_support
                    .iter()
                    .max()
                    .ok_or("max not found in this_hap_support")?;
                let this_hap_support_sum = this_hap_support.iter().sum::<usize>();
                if *this_hap_support_max as i32 <= *longest_support_max as i32 - 2
                    && this_hap_support_sum as i32 <= longest_support_sum as i32 - 2
                {
                    if filtered_candidates.contains(test_hap) {
                        let index = filtered_candidates
                            .iter()
                            .position(|x| *x == *test_hap)
                            .ok_or("item not found")?;
                        filtered_candidates.remove(index);
                    }
                }
            }
        } else if *max_nread == 2 {
            for test_hap in unique_support.keys() {
                let this_hap_support = support_length_per_hap
                    .get(test_hap)
                    .ok_or("key not found in support_length_per_hap")?;
                let this_hap_support_max = this_hap_support
                    .iter()
                    .max()
                    .ok_or("max not found in this_hap_support")?;
                let this_hap_support_sum = this_hap_support.iter().sum::<usize>();
                if *this_hap_support_max as i32 <= *longest_support_max as i32 - 2
                    && this_hap_support_sum as i32 <= longest_support_sum as i32 - 3
                {
                    if filtered_candidates.contains(test_hap) {
                        let index = filtered_candidates
                            .iter()
                            .position(|x| *x == *test_hap)
                            .ok_or("item not found")?;
                        filtered_candidates.remove(index);
                    }
                }
            }
        }
        // a diff of two in longest support, and a diff of 4 in sum of support length
        // at least two reads with longer support when the best candidate is compared against this one
        else {
            for test_hap in unique_support.keys() {
                let this_hap_support = support_length_per_hap
                    .get(test_hap)
                    .ok_or("key not found in support_length_per_hap")?;
                let this_hap_support_max = this_hap_support
                    .iter()
                    .max()
                    .ok_or("max not found in this_hap_support")?;
                let this_hap_support_sum = this_hap_support.iter().sum::<usize>();
                if (*this_hap_support_max as i32 <= *longest_support_max as i32 - 2
                    && this_hap_support_sum as i32 <= longest_support_sum as i32 - 4)
                    || (*this_hap_support_max as i32 <= *longest_support_max as i32 - 1
                        && this_hap_support_sum as i32 <= longest_support_sum as i32 - 9)
                {
                    let longest_support_above_this_max = longest_support
                        .iter()
                        .filter(|x| **x > *this_hap_support_max)
                        .collect::<Vec<_>>();
                    if longest_support_above_this_max.len() >= 2 {
                        if filtered_candidates.contains(test_hap) {
                            let index = filtered_candidates
                                .iter()
                                .position(|x| *x == *test_hap)
                                .ok_or("item not found")?;
                            filtered_candidates.remove(index);
                        }
                    }
                }
            }
        }
    }
    Ok(filtered_candidates)
}

/// Get reads that support each haplotype and possible haplotypes of each read
/// # Arguments
/// * `haplotype_per_read` - read name -> read hap (vec of fingerprints)
/// * `hap_list` - list of haplotypes
/// * `req_overlap_len` - required length of overlap for supporting reads
/// * `reverse` - checking backward from the end if true
/// # Returns
/// * `HaplotypeReads` - haplotype reads
pub fn match_reads_and_haplotypes(
    haplotype_per_read: BTreeMap<String, Vec<i32>>,
    hap_list: Vec<Vec<i32>>,
    req_overlap_len: Option<usize>,
    reverse: bool,
) -> HaplotypeReads {
    let mut support_reads: BTreeMap<Vec<i32>, Vec<Vec<i32>>> = BTreeMap::new();
    let mut supporting_haps_per_read = BTreeMap::new();
    for (read_name, read_hap) in haplotype_per_read.iter() {
        let mut matching_haplotypes = Vec::new();
        for haplotype_to_extend in &hap_list {
            let is_match: bool = if let Some(overlap_len) = req_overlap_len {
                if haplotype_to_extend.len() < overlap_len {
                    false
                } else if !reverse {
                    let haplotype_to_extend_len = haplotype_to_extend.len();
                    trace!("haplotype_to_extend_len {haplotype_to_extend_len} overlap_len {overlap_len}");
                    two_haplotypes_are_matching(
                        read_hap.to_vec(),
                        haplotype_to_extend.clone(),
                        Some((
                            haplotype_to_extend_len - overlap_len,
                            haplotype_to_extend_len,
                        )),
                    )
                } else {
                    let haplotype_to_extend_len = haplotype_to_extend.len();
                    trace!("haplotype_to_extend_len {haplotype_to_extend_len} overlap_len {overlap_len}");
                    two_haplotypes_are_matching(
                        read_hap.to_vec(),
                        haplotype_to_extend.clone(),
                        Some((0, overlap_len)),
                    )
                }
            } else {
                two_haplotypes_are_matching(read_hap.to_vec(), haplotype_to_extend.clone(), None)
            };
            if is_match {
                matching_haplotypes.push(haplotype_to_extend.clone());
            }
        }
        supporting_haps_per_read
            .entry(read_name.to_string())
            .or_insert(matching_haplotypes.clone());
        if matching_haplotypes.len() == 1 {
            support_reads
                .entry(matching_haplotypes[0].clone())
                .or_default()
                .push(read_hap.to_vec());
        }
    }
    HaplotypeReads {
        unique: support_reads,
        by_read: supporting_haps_per_read,
    }
}

/// Remove a complete allele if there is an overlapping incomplete allele longer than it.
/// The complete allele has to be at least CN3.
/// The complete allele can have at most one unique node than the incomplete allele.
/// # Arguments
/// * `complete` - complete alleles
/// * `incomplete` - incomplete alleles
/// # Returns
/// * `HashSet<Vec<i32>>` - complete alleles that overlap with incomplete alleles
pub fn complete_overlapping_incomplete(
    complete: Vec<Vec<i32>>,
    incomplete: Vec<Vec<i32>>,
) -> HashSet<Vec<i32>> {
    let mut complete_hap_overlapping_incomplete = HashSet::new();
    for a in complete {
        for b in &incomplete {
            let a_nodes: Vec<i32> = a.iter().filter(|x| **x > 0).copied().collect::<Vec<_>>();
            let b_nodes: Vec<i32> = b.iter().filter(|x| **x > 0).copied().collect::<Vec<_>>();
            let a_nodes_len = a_nodes.len() as i32;
            let b_nodes_len = b_nodes.len() as i32;
            let a_set_len = a_nodes.iter().cloned().collect::<HashSet<i32>>().len() as i32;
            let b_set_len = b_nodes.iter().cloned().collect::<HashSet<i32>>().len() as i32;
            let a_uniq_node_count = a_set_len - b_set_len;
            if a_uniq_node_count <= 1
                && a_nodes_len >= 3
                && b_nodes_len > a_nodes_len
                && (b_is_contained_in_a(&b_nodes, &a_nodes[1..])
                    || b_is_contained_in_a(&b_nodes, &a_nodes[..(a_nodes_len as usize - 1)]))
            {
                debug!("complete hap {a:?} overlaps an incomplete hap {b:?}");
                complete_hap_overlapping_incomplete.insert(a.clone());
            }
        }
    }
    complete_hap_overlapping_incomplete
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_get_edges() {
        let mut read_edges = BTreeMap::new();
        read_edges
            .entry(String::from("read1"))
            .or_insert(vec![1, 2, 2]);
        read_edges
            .entry(String::from("read2"))
            .or_insert(vec![2, 2, 3]);
        read_edges
            .entry(String::from("read3"))
            .or_insert(vec![1, 2]);
        read_edges
            .entry(String::from("read4"))
            .or_insert(vec![1, 0]);
        let mut graph = build_graph(read_edges, 2);
        graph.get_edges();
        assert!(graph.edges.contains_key(&(2, 3)));
        assert!(graph.edges.contains_key(&(2, 2)));
        assert!(graph.edges.contains_key(&(1, 2)));
        assert!(!graph.edges.contains_key(&(1, 0)));
        assert!(graph.cyclic_nodes.contains_key(&2));
        assert_eq!(graph.cyclic_nodes[&2], vec![2, 2]);
    }

    #[test]
    fn test_clean_up_edges() {
        let mut read_edges = BTreeMap::new();
        read_edges
            .entry(String::from("read1"))
            .or_insert(vec![1, 2, 3]);
        read_edges
            .entry(String::from("read2"))
            .or_insert(vec![2, 3, 4]);
        read_edges
            .entry(String::from("read3"))
            .or_insert(vec![1, 2]);
        // this is a suspicious edge
        read_edges
            .entry(String::from("read4"))
            .or_insert(vec![3, 2]);
        read_edges
            .entry(String::from("read5"))
            .or_insert(vec![3, 4]);
        let mut graph = build_graph(read_edges, 2);
        graph.get_edges();
        graph.analyze_nodes(false);
        let _ = graph.clean_up_edges(false);
        assert!(!graph.clean_edges.contains_key(&(3, 2)));

        let mut read_edges = BTreeMap::new();
        read_edges
            .entry(String::from("read1"))
            .or_insert(vec![1, 2, 3]);
        read_edges
            .entry(String::from("read2"))
            .or_insert(vec![2, 3, 4]);
        read_edges
            .entry(String::from("read3"))
            .or_insert(vec![1, 2]);
        // now two reads support this edge
        read_edges
            .entry(String::from("read4"))
            .or_insert(vec![3, 2]);
        read_edges
            .entry(String::from("read6"))
            .or_insert(vec![3, 2]);
        read_edges
            .entry(String::from("read5"))
            .or_insert(vec![3, 4]);
        let mut graph = build_graph(read_edges, 2);
        graph.get_edges();
        graph.analyze_nodes(false);
        let _ = graph.clean_up_edges(false);
        assert!(graph.clean_edges.contains_key(&(3, 2)));

        let mut read_edges = BTreeMap::new();
        read_edges
            .entry(String::from("read1"))
            .or_insert(vec![1, 2, 3]);
        read_edges
            .entry(String::from("read2"))
            .or_insert(vec![2, 3, 4]);
        read_edges
            .entry(String::from("read3"))
            .or_insert(vec![1, 2]);
        // this is edge has one count but 4 is the only node leading to 5
        read_edges
            .entry(String::from("read4"))
            .or_insert(vec![4, 5]);
        read_edges
            .entry(String::from("read5"))
            .or_insert(vec![3, 4]);
        let mut graph = build_graph(read_edges, 2);
        graph.get_edges();
        graph.analyze_nodes(false);
        let _ = graph.clean_up_edges(false);
        assert!(graph.clean_edges.contains_key(&(4, 5)));

        // no removal. complex scenario,
        // the next nodes of 3 are [2, 3] and prev nodes of 2 are [2, 3]
        let mut read_edges = BTreeMap::new();
        read_edges
            .entry(String::from("read1"))
            .or_insert(vec![1, 2, 2]);
        read_edges
            .entry(String::from("read2"))
            .or_insert(vec![2, 3]);
        read_edges
            .entry(String::from("read4"))
            .or_insert(vec![3, 2]);
        read_edges
            .entry(String::from("read6"))
            .or_insert(vec![2, 2]);
        read_edges
            .entry(String::from("read5"))
            .or_insert(vec![2, 3]);
        read_edges
            .entry(String::from("read5"))
            .or_insert(vec![3, 3]);
        let mut graph = build_graph(read_edges, 2);
        graph.get_edges();
        graph.analyze_nodes(false);
        let _ = graph.clean_up_edges(false);
        assert!(graph.clean_edges.contains_key(&(3, 2)));

        // the next nodes of 3 are [2, 3] and prev nodes of 2 are [2, 3, 4]
        // now can be removed
        let mut read_edges = BTreeMap::new();
        read_edges
            .entry(String::from("read1"))
            .or_insert(vec![1, 2, 2]);
        read_edges
            .entry(String::from("read2"))
            .or_insert(vec![2, 3, 4]);
        read_edges
            .entry(String::from("read4"))
            .or_insert(vec![3, 2]);
        read_edges
            .entry(String::from("read6"))
            .or_insert(vec![2, 2]);
        read_edges
            .entry(String::from("read5"))
            .or_insert(vec![2, 3]);
        read_edges
            .entry(String::from("read5"))
            .or_insert(vec![3, 3]);
        let mut graph = build_graph(read_edges, 2);
        graph.get_edges();
        graph.analyze_nodes(false);
        let _ = graph.clean_up_edges(false);
        assert!(!graph.clean_edges.contains_key(&(3, 2)));
    }

    #[test]
    fn test_clean_up_edges_allow_low_support_edges_to_ends_false() {
        let mut read_edges = BTreeMap::new();
        read_edges
            .entry(String::from("read1"))
            .or_insert(vec![1, 2, 3]);
        read_edges
            .entry(String::from("read2"))
            .or_insert(vec![2, 3, 4]);
        read_edges
            .entry(String::from("read3"))
            .or_insert(vec![1, 2]);
        // this is a suspicious edge
        read_edges
            .entry(String::from("read4"))
            .or_insert(vec![3, -10]);
        read_edges
            .entry(String::from("read5"))
            .or_insert(vec![3, 4]);
        read_edges
            .entry(String::from("read6"))
            .or_insert(vec![4, -10]);
        read_edges
            .entry(String::from("read7"))
            .or_insert(vec![4, -10]);
        let mut graph = build_graph(read_edges, 2);
        graph.get_edges();
        graph.analyze_nodes(false);
        let _ = graph.clean_up_edges(false);
        assert!(!graph.clean_edges.contains_key(&(3, -10)));
    }

    #[test]
    fn test_clean_up_edges_allow_low_support_edges_to_ends_true() {
        let mut read_edges = BTreeMap::new();
        read_edges
            .entry(String::from("read1"))
            .or_insert(vec![1, 2, 3]);
        read_edges
            .entry(String::from("read2"))
            .or_insert(vec![2, 3, 4]);
        read_edges
            .entry(String::from("read3"))
            .or_insert(vec![1, 2]);
        // this is a suspicious edge
        read_edges
            .entry(String::from("read4"))
            .or_insert(vec![3, -10]);
        read_edges
            .entry(String::from("read5"))
            .or_insert(vec![3, 4]);
        read_edges
            .entry(String::from("read6"))
            .or_insert(vec![4, -10]);
        read_edges
            .entry(String::from("read7"))
            .or_insert(vec![4, -10]);
        let mut graph = build_graph(read_edges, 2);
        graph.get_edges();
        graph.analyze_nodes(false);
        let _ = graph.clean_up_edges(true);
        assert!(graph.clean_edges.contains_key(&(3, -10)));
    }

    #[test]
    fn test_clean_up_edges_allow_low_support_edges_to_ends_true_too_many_links_to_end() {
        let mut read_edges = BTreeMap::new();
        read_edges
            .entry(String::from("read1"))
            .or_insert(vec![1, 2, 3]);
        read_edges
            .entry(String::from("read2"))
            .or_insert(vec![2, 3, 4]);
        read_edges
            .entry(String::from("read3"))
            .or_insert(vec![1, 2]);
        // this is a suspicious edge
        read_edges
            .entry(String::from("read4"))
            .or_insert(vec![3, -10]);
        read_edges
            .entry(String::from("read5"))
            .or_insert(vec![3, 4]);
        read_edges
            .entry(String::from("read6"))
            .or_insert(vec![4, -10]);
        read_edges
            .entry(String::from("read7"))
            .or_insert(vec![4, -10]);
        read_edges
            .entry(String::from("read8"))
            .or_insert(vec![5, -10]);
        read_edges
            .entry(String::from("read9"))
            .or_insert(vec![5, -10]);
        read_edges
            .entry(String::from("read10"))
            .or_insert(vec![6, -10]);
        read_edges
            .entry(String::from("read11"))
            .or_insert(vec![6, -10]);
        let mut graph = build_graph(read_edges, 2);
        graph.get_edges();
        graph.analyze_nodes(false);
        let _ = graph.clean_up_edges(true);
        assert!(!graph.clean_edges.contains_key(&(3, -10)));
    }

    #[test]
    fn test_include_cyclic_node_forward() {
        let mut read_edges = BTreeMap::new();
        read_edges
            .entry(String::from("read1"))
            .or_insert(vec![1, 2, 3]);
        let mut graph = build_graph(read_edges, 2);
        graph.cyclic_nodes.entry(1).or_insert(vec![1, 2]);
        let to_include = graph
            .include_cyclic_node_forward(1, 1, vec![2, 1, 1])
            .unwrap();
        assert_eq!(to_include, true);
        let to_include = graph
            .include_cyclic_node_forward(1, 1, vec![2, 1, 1, 1])
            .unwrap();
        assert_eq!(to_include, false);
    }

    #[test]
    fn test_include_cyclic_node_backward() {
        let mut read_edges = BTreeMap::new();
        read_edges
            .entry(String::from("read1"))
            .or_insert(vec![1, 2, 3]);
        let mut graph = build_graph(read_edges, 2);
        graph.cyclic_nodes.entry(1).or_insert(vec![1, 2]);
        let to_include = graph
            .include_cyclic_node_backward(1, 1, vec![1, 1, 2])
            .unwrap();
        assert_eq!(to_include, true);
        let to_include = graph
            .include_cyclic_node_backward(1, 1, vec![1, 1, 1, 2])
            .unwrap();
        assert_eq!(to_include, false);
    }

    #[test]
    fn test_require_overlap_length_check_forward() {
        let mut read_edges = BTreeMap::new();
        read_edges
            .entry(String::from("read1"))
            .or_insert(vec![1, 2, 3, 4, 6]);
        read_edges
            .entry(String::from("read2"))
            .or_insert(vec![6, 2, 3, 4, 7]);
        let mut graph = build_graph(read_edges, 2);
        graph.clean_reads = graph.reads.clone();
        let this_hap = vec![2, 3, 4, 5, 6];
        let ovl_len = graph
            .require_overlap_length_check_forward(this_hap, false)
            .unwrap();
        assert_eq!(ovl_len, 3);

        let this_hap = vec![1, 2, 3, 5, 6];
        let ovl_len = graph
            .require_overlap_length_check_forward(this_hap, false)
            .unwrap();
        assert_eq!(ovl_len, 0);
    }

    #[test]
    fn test_require_overlap_length_check_forward_strict_low_cyclic_count() {
        let mut read_edges = BTreeMap::new();
        read_edges
            .entry(String::from("read2"))
            .or_insert(vec![1, 3, 3]);
        read_edges
            .entry(String::from("read3"))
            .or_insert(vec![3, 3, 2]);
        let mut graph = build_graph(read_edges, 2);
        graph.clean_reads = graph.reads.clone();
        let this_hap = vec![3, 3, 2];
        let ovl_len = graph
            .require_overlap_length_check_forward(this_hap, true)
            .unwrap();
        assert_eq!(ovl_len, 1);

        let this_hap = vec![3, 3, 2];
        let ovl_len = graph
            .require_overlap_length_check_forward(this_hap, false)
            .unwrap();
        assert_eq!(ovl_len, 1);
    }

    #[test]
    fn test_require_overlap_length_check_forward_strict_high_cyclic_count() {
        let mut read_edges = BTreeMap::new();
        read_edges
            .entry(String::from("read2"))
            .or_insert(vec![1, 3, 3]);
        read_edges
            .entry(String::from("read3"))
            .or_insert(vec![3, 3, 2]);
        let mut graph = build_graph(read_edges, 2);
        graph.clean_reads = graph.reads.clone();
        let this_hap = vec![3, 3, 3, 2];
        let ovl_len = graph
            .require_overlap_length_check_forward(this_hap, true)
            .unwrap();
        assert_eq!(ovl_len, 3);

        let this_hap = vec![3, 3, 3, 2];
        let ovl_len = graph
            .require_overlap_length_check_forward(this_hap, false)
            .unwrap();
        assert_eq!(ovl_len, 1);
    }

    #[test]
    fn test_require_overlap_length_check_backward() {
        let mut read_edges = BTreeMap::new();
        read_edges
            .entry(String::from("read1"))
            .or_insert(vec![1, 2, 3, 4, 6]);
        read_edges
            .entry(String::from("read2"))
            .or_insert(vec![1, 5, 3, 4, 7]);
        let mut graph = build_graph(read_edges, 2);
        graph.clean_reads = graph.reads.clone();
        let this_hap = vec![1, 2, 3, 4];
        let ovl_len = graph
            .require_overlap_length_check_backward(this_hap, false)
            .unwrap();
        assert_eq!(ovl_len, 2);

        let this_hap = vec![1, 2, 3, 4, 6];
        let ovl_len = graph
            .require_overlap_length_check_backward(this_hap, false)
            .unwrap();
        assert_eq!(ovl_len, 0);
    }

    #[test]
    fn test_require_overlap_length_check_backward_strict_low_cyclic_count() {
        let mut read_edges = BTreeMap::new();
        read_edges
            .entry(String::from("read2"))
            .or_insert(vec![1, 3, 3]);
        read_edges
            .entry(String::from("read3"))
            .or_insert(vec![3, 3, 2]);
        let mut graph = build_graph(read_edges, 2);
        graph.clean_reads = graph.reads.clone();
        let this_hap = vec![1, 3, 3];
        let ovl_len = graph
            .require_overlap_length_check_backward(this_hap, true)
            .unwrap();
        assert_eq!(ovl_len, 1);

        let this_hap = vec![1, 3, 3];
        let ovl_len = graph
            .require_overlap_length_check_backward(this_hap, false)
            .unwrap();
        assert_eq!(ovl_len, 1);
    }

    #[test]
    fn test_require_overlap_length_check_backward_strict_high_cyclic_count() {
        let mut read_edges = BTreeMap::new();
        read_edges
            .entry(String::from("read2"))
            .or_insert(vec![1, 3, 3]);
        read_edges
            .entry(String::from("read3"))
            .or_insert(vec![3, 3, 2]);
        let mut graph = build_graph(read_edges, 2);
        graph.clean_reads = graph.reads.clone();
        let this_hap = vec![1, 3, 3, 3];
        let ovl_len = graph
            .require_overlap_length_check_backward(this_hap, true)
            .unwrap();
        assert_eq!(ovl_len, 3);

        let this_hap = vec![1, 3, 3, 3];
        let ovl_len = graph
            .require_overlap_length_check_backward(this_hap, false)
            .unwrap();
        assert_eq!(ovl_len, 1);
    }

    #[test]
    fn test_get_unused_nodes() {
        let mut read_edges = BTreeMap::new();
        read_edges
            .entry(String::from("read1"))
            .or_insert(vec![1, 2, 3, 4, 6]);
        let mut graph = build_graph(read_edges, 2);
        graph.nodes.insert(1);
        graph.nodes.insert(2);
        graph.nodes.insert(3);
        graph.nodes.insert(4);
        let test_haps: Vec<Vec<i32>> = vec![vec![1, 2], vec![1, 2, 3]];
        let unused_nodes = graph.get_unused_nodes(test_haps);
        assert_eq!(unused_nodes, vec![4]);

        let test_haps: Vec<Vec<i32>> = vec![vec![1, 2], vec![1, 2, 3], vec![2, 4]];
        let unused_nodes = graph.get_unused_nodes(test_haps);
        assert!(unused_nodes.is_empty());
    }

    #[test]
    fn test_pick_from_candidates() {
        let mut read_edges = BTreeMap::new();
        read_edges
            .entry(String::from("read1"))
            .or_insert(vec![1, 2, 3, 4, 6]);
        read_edges
            .entry(String::from("read2"))
            .or_insert(vec![1, 2, 3, 4, 5]);
        read_edges
            .entry(String::from("read3"))
            .or_insert(vec![2, 3, 4, 5]);

        let graph = build_graph(read_edges.clone(), 2);
        // highest support is 2, too low
        // pick the shortest one, but same length here
        let candidates = vec![vec![3, 4, 5], vec![3, 4, 6]];
        let (picked, has_support) = graph.pick_from_candidates(candidates).unwrap();
        assert!(!has_support);
        assert_eq!(picked, vec![3, 4, 5]);

        read_edges
            .entry(String::from("read4"))
            .or_insert(vec![3, 4, 5]);
        let graph = build_graph(read_edges.clone(), 2);
        // pick the one with highest support
        let candidates = vec![vec![3, 4, 5], vec![3, 4, 6]];
        let (picked, has_support) = graph.pick_from_candidates(candidates).unwrap();
        assert!(has_support);
        assert_eq!(picked, vec![3, 4, 5]);

        // without support, pick the shortest one
        let candidates = vec![vec![2, 4, 5], vec![1, 2, 4, 6]];
        let (picked, has_support) = graph.pick_from_candidates(candidates).unwrap();
        assert!(!has_support);
        assert_eq!(picked, vec![2, 4, 5]);
    }

    #[test]
    fn test_merge_two_incomplete() {
        let mut read_edges = BTreeMap::new();
        let graph = build_graph(read_edges.clone(), 2);

        let complete = vec![vec![-1, 9, 2, 11, 12, 5, -10]];
        // overlap by 2
        let incomplete = vec![vec![2, 3, 4, 5, 6, -10], vec![-1, 2, 3]];
        let merge_result = graph.merge_two_incomplete(complete, incomplete).unwrap();
        assert_eq!(
            merge_result.complete,
            vec![vec![-1, 9, 2, 11, 12, 5, -10], vec![-1, 2, 3, 4, 5, 6, -10]]
        );
        assert!(merge_result.incomplete.is_empty());

        let complete = vec![vec![-1, 9, 2, 11, 12, 5, -10]];
        // overlap by 1
        let incomplete = vec![vec![2, 3, 4, 5, 6, -10], vec![-1, 2]];
        let merge_result = graph.merge_two_incomplete(complete, incomplete).unwrap();
        assert_eq!(
            merge_result.complete,
            vec![vec![-1, 9, 2, 11, 12, 5, -10], vec![-1, 2, 3, 4, 5, 6, -10]]
        );
        assert!(merge_result.incomplete.is_empty());

        let complete = vec![vec![-1, 9, 2, 11, 12, 5, -10]];
        // no overlap
        let incomplete = vec![vec![3, 4, 5, 6, -10], vec![-1, 2]];
        let merge_result = graph.merge_two_incomplete(complete, incomplete).unwrap();
        assert_eq!(merge_result.complete, vec![vec![-1, 9, 2, 11, 12, 5, -10]]);
        assert_eq!(
            merge_result.incomplete,
            vec![vec![3, 4, 5, 6, -10], vec![-1, 2]]
        );

        // no overlap, without right reads
        read_edges
            .entry(String::from("read1"))
            .or_insert(vec![-1, 2, 3]);
        let mut graph = build_graph(read_edges.clone(), 2);
        graph.get_edges();
        let complete = vec![vec![-1, 9, 2, 11, 12, 5, -10]];
        let incomplete = vec![vec![3, 4, 5, 6, -10], vec![-1, 2]];
        let merge_result = graph.merge_two_incomplete(complete, incomplete).unwrap();
        assert_eq!(merge_result.complete, vec![vec![-1, 9, 2, 11, 12, 5, -10]]);
        assert_eq!(
            merge_result.incomplete,
            vec![vec![3, 4, 5, 6, -10], vec![-1, 2]]
        );

        // no overlap, with reads
        read_edges
            .entry(String::from("read3"))
            .or_insert(vec![-1, 2, 3, 4]);
        let mut graph = build_graph(read_edges.clone(), 2);
        graph.get_edges();
        let complete = vec![vec![-1, 9, 2, 11, 12, 5, -10]];
        let incomplete = vec![vec![3, 4, 5, 6, -10], vec![-1, 2]];
        let merge_result = graph.merge_two_incomplete(complete, incomplete).unwrap();
        assert_eq!(
            merge_result.complete,
            vec![vec![-1, 9, 2, 11, 12, 5, -10], vec![-1, 2, 3, 4, 5, 6, -10]]
        );
        assert!(merge_result.incomplete.is_empty());
    }

    #[test]
    fn test_b_is_contained_in_a() {
        let a = vec![1, 2, 3, 4, 5];
        let b = vec![3, 4];
        assert!(b_is_contained_in_a(&a, &b));

        let a = vec![1, 2, 3, 4, 5];
        let b = vec![1, 2];
        assert!(b_is_contained_in_a(&a, &b));

        let a = vec![1, 2, 3, 4, 5];
        let b = vec![4, 5];
        assert!(b_is_contained_in_a(&a, &b));

        let a = vec![1, 2, 3, 4, 5];
        let b = vec![3, 5];
        assert!(!b_is_contained_in_a(&a, &b));

        let a = vec![1, 2];
        let b = vec![3, 5];
        assert!(!b_is_contained_in_a(&a, &b));
        // a and b are identical
        let a = vec![3, 5];
        let b = vec![3, 5];
        assert!(b_is_contained_in_a(&a, &b));
        // b is bigger
        let a = vec![3, 4];
        let b = vec![1, 2, 3, 4, 5];
        assert!(!b_is_contained_in_a(&a, &b));
    }

    #[test]
    fn test_complete_overlapping_incomplete() {
        let complete = vec![vec![-1, 1, 2, 3, 4, 5, -10], vec![-1, 6, 7, 8, 9, 10, -10]];
        let incomplete = vec![vec![1, 2, 3, 4, 5, 6, -10]];
        let overlap = complete_overlapping_incomplete(complete, incomplete);
        let mut expected_overlap = HashSet::new();
        expected_overlap.insert(vec![-1, 1, 2, 3, 4, 5, -10]);
        assert_eq!(overlap, expected_overlap);
        // the complete allele has at most one unique node
        let complete = vec![
            vec![-1, 1, 2, 3, 4, 5, 7, -10],
            vec![-1, 6, 7, 8, 9, 10, -10],
        ];
        let incomplete = vec![vec![1, 2, 3, 4, 5, 6, 6, -10]];
        let overlap = complete_overlapping_incomplete(complete, incomplete);
        let mut expected_overlap = HashSet::new();
        expected_overlap.insert(vec![-1, 1, 2, 3, 4, 5, 7, -10]);
        assert_eq!(overlap, expected_overlap);
        // the complete allele has to be at least CN3
        let complete = vec![vec![-1, 1, 2, -10], vec![-1, 6, 7, 8, 9, 10, -10]];
        let incomplete = vec![vec![1, 2, 3, -10]];
        let overlap = complete_overlapping_incomplete(complete, incomplete);
        let expected_overlap = HashSet::new();
        assert_eq!(overlap, expected_overlap);
        // the incomplete allele has to be longer than the complete allele
        let complete = vec![
            vec![-1, 1, 2, 3, 4, 5, 7, -10],
            vec![-1, 6, 7, 8, 9, 10, -10],
        ];
        let incomplete = vec![vec![1, 2, 3, 4, 5, 6, -10]];
        let overlap = complete_overlapping_incomplete(complete, incomplete);
        let expected_overlap = HashSet::new();
        assert_eq!(overlap, expected_overlap);
    }

    #[test]
    fn test_check_start_or_end_reads_without_match() {
        let mut read_edges = BTreeMap::new();
        read_edges
            .entry(String::from("read1"))
            .or_insert(vec![-1, 2, 3, 4, 6]);
        read_edges
            .entry(String::from("read2"))
            .or_insert(vec![1, 5, 3, 4, -10]);
        read_edges
            .entry(String::from("read3"))
            .or_insert(vec![0, -10]);
        let mut graph = build_graph(read_edges, 2);
        graph.clean_reads = graph.reads.clone();
        let mut support_by_read = BTreeMap::new();
        support_by_read.entry(String::from("read1")).or_default();
        support_by_read.entry(String::from("read2")).or_default();
        support_by_read.entry(String::from("read3")).or_default();
        let nread = graph.check_start_or_end_reads_without_match(support_by_read);
        assert_eq!(nread, 2);
    }

    #[test]
    fn test_filter_assembled_candidates() {
        let mut unique_support = BTreeMap::new();
        let mut support_length = BTreeMap::new();
        unique_support
            .entry(vec![1, 2, 3, 4, 5])
            .or_insert(vec![vec![3, 4, 5], vec![4, 5]]);
        unique_support
            .entry(vec![1, 2, 3, 4, 6])
            .or_insert(vec![vec![3, 4, 6], vec![4, 6]]);
        support_length.entry(vec![1, 2, 3, 4, 5]).or_insert(vec![
            (vec![3 as i32, 4 as i32, 5 as i32], 3 as usize),
            (vec![4 as i32, 5 as i32], 2 as usize),
        ]);
        support_length.entry(vec![1, 2, 3, 4, 6]).or_insert(vec![
            (vec![3 as i32, 4 as i32, 6 as i32], 3 as usize),
            (vec![4 as i32, 6 as i32], 2 as usize),
        ]);
        let candidates = filter_assembled_candidates(unique_support, support_length).unwrap();
        assert_eq!(candidates, vec![vec![1, 2, 3, 4, 5], vec![1, 2, 3, 4, 6]]);

        // one read vs 10 reads or more
        let mut unique_support = BTreeMap::new();
        let mut support_length = BTreeMap::new();
        unique_support.entry(vec![1, 2, 3, 4, 5]).or_insert(vec![
            vec![3, 4, 5],
            vec![3, 4, 5],
            vec![4, 5],
            vec![4, 5],
            vec![4, 5],
            vec![4, 5],
            vec![4, 5],
            vec![4, 5],
            vec![4, 5],
            vec![4, 5],
        ]);
        unique_support
            .entry(vec![1, 2, 3, 4, 6])
            .or_insert(vec![vec![4, 6]]);
        support_length.entry(vec![1, 2, 3, 4, 5]).or_insert(vec![
            (vec![3 as i32, 4 as i32, 5 as i32], 3 as usize),
            (vec![3 as i32, 4 as i32, 5 as i32], 3 as usize),
            (vec![4 as i32, 5 as i32], 2 as usize),
            (vec![4 as i32, 5 as i32], 2 as usize),
            (vec![4 as i32, 5 as i32], 2 as usize),
            (vec![4 as i32, 5 as i32], 2 as usize),
            (vec![4 as i32, 5 as i32], 2 as usize),
            (vec![4 as i32, 5 as i32], 2 as usize),
            (vec![4 as i32, 5 as i32], 2 as usize),
            (vec![4 as i32, 5 as i32], 2 as usize),
        ]);
        support_length
            .entry(vec![1, 2, 3, 4, 6])
            .or_insert(vec![(vec![4 as i32, 6 as i32], 2 as usize)]);
        let candidates = filter_assembled_candidates(unique_support, support_length).unwrap();
        assert_eq!(candidates, vec![vec![1, 2, 3, 4, 5]]);

        // better one has the max support seen twice and the worse one has only two read support
        let mut unique_support = BTreeMap::new();
        let mut support_length = BTreeMap::new();
        unique_support.entry(vec![1, 2, 3, 4, 5]).or_insert(vec![
            vec![3, 4, 5],
            vec![3, 4, 5],
            vec![4, 5],
        ]);
        unique_support
            .entry(vec![1, 2, 3, 4, 6])
            .or_insert(vec![vec![4, 6], vec![4, 6]]);
        support_length.entry(vec![1, 2, 3, 4, 5]).or_insert(vec![
            (vec![3 as i32, 4 as i32, 5 as i32], 3 as usize),
            (vec![3 as i32, 4 as i32, 5 as i32], 3 as usize),
            (vec![4 as i32, 5 as i32], 2 as usize),
        ]);
        support_length.entry(vec![1, 2, 3, 4, 6]).or_insert(vec![
            (vec![4 as i32, 6 as i32], 2 as usize),
            (vec![4 as i32, 6 as i32], 2 as usize),
        ]);
        let candidates = filter_assembled_candidates(unique_support, support_length).unwrap();
        assert_eq!(candidates, vec![vec![1, 2, 3, 4, 5]]);

        // better one has the max support seen twice and the worse one has only two read support
        let mut unique_support = BTreeMap::new();
        let mut support_length = BTreeMap::new();
        unique_support
            .entry(vec![1, 2, 3, 4, 5])
            .or_insert(vec![vec![3, 4, 5], vec![3, 4, 5]]);
        unique_support.entry(vec![1, 2, 3, 4, 6]).or_insert(vec![
            vec![4, 6],
            vec![4, 6],
            vec![4, 6],
        ]);
        support_length.entry(vec![1, 2, 3, 4, 5]).or_insert(vec![
            (vec![3 as i32, 4 as i32, 5 as i32], 3 as usize),
            (vec![3 as i32, 4 as i32, 5 as i32], 3 as usize),
        ]);
        support_length.entry(vec![1, 2, 3, 4, 6]).or_insert(vec![
            (vec![4 as i32, 6 as i32], 2 as usize),
            (vec![4 as i32, 6 as i32], 2 as usize),
            (vec![4 as i32, 6 as i32], 2 as usize),
        ]);
        let candidates = filter_assembled_candidates(unique_support, support_length).unwrap();
        assert_eq!(candidates, vec![vec![1, 2, 3, 4, 5], vec![1, 2, 3, 4, 6]]);

        // a diff of two in longest support, and a diff of 4 in sum of support length
        // at least two reads with longer support when the best candidate is compared against this one
        let mut unique_support = BTreeMap::new();
        let mut support_length = BTreeMap::new();
        unique_support
            .entry(vec![1, 2, 3, 4, 5])
            .or_insert(vec![vec![2, 3, 4, 5], vec![2, 3, 4, 5]]);
        unique_support
            .entry(vec![1, 2, 3, 4, 6])
            .or_insert(vec![vec![4, 6], vec![4, 6]]);
        // sum 8 vs. 4
        // longest 4 vs. 2
        support_length.entry(vec![1, 2, 3, 4, 5]).or_insert(vec![
            (vec![2 as i32, 3 as i32, 4 as i32, 5 as i32], 4 as usize),
            (vec![2 as i32, 3 as i32, 4 as i32, 5 as i32], 4 as usize),
        ]);
        support_length.entry(vec![1, 2, 3, 4, 6]).or_insert(vec![
            (vec![4 as i32, 6 as i32], 2 as usize),
            (vec![4 as i32, 6 as i32], 2 as usize),
        ]);
        let candidates = filter_assembled_candidates(unique_support, support_length).unwrap();
        assert_eq!(candidates, vec![vec![1, 2, 3, 4, 5]]);
        // a diff of 1 in longest support, and a diff of 8 in sum of support length
        // at least two reads with longer support when the best candidate is compared against this one
        let mut unique_support = BTreeMap::new();
        let mut support_length = BTreeMap::new();
        unique_support.entry(vec![1, 2, 3, 4, 5]).or_insert(vec![
            vec![3, 4, 5],
            vec![3, 4, 5],
            vec![4, 5],
            vec![4, 5],
            vec![4, 5],
            vec![4, 5],
            vec![4, 5],
        ]);
        unique_support.entry(vec![1, 2, 3, 4, 6]).or_insert(vec![
            vec![4, 6],
            vec![4, 6],
            vec![4, 6],
        ]);
        // sum 16 vs. 6
        // longest 3 vs. 2
        support_length.entry(vec![1, 2, 3, 4, 5]).or_insert(vec![
            (vec![3 as i32, 4 as i32, 5 as i32], 3 as usize),
            (vec![3 as i32, 4 as i32, 5 as i32], 3 as usize),
            (vec![4 as i32, 5 as i32], 2 as usize),
            (vec![4 as i32, 5 as i32], 2 as usize),
            (vec![4 as i32, 5 as i32], 2 as usize),
            (vec![4 as i32, 5 as i32], 2 as usize),
            (vec![4 as i32, 5 as i32], 2 as usize),
        ]);
        support_length.entry(vec![1, 2, 3, 4, 6]).or_insert(vec![
            (vec![4 as i32, 6 as i32], 2 as usize),
            (vec![4 as i32, 6 as i32], 2 as usize),
            (vec![4 as i32, 6 as i32], 2 as usize),
        ]);
        let candidates = filter_assembled_candidates(unique_support, support_length).unwrap();
        assert_eq!(candidates, vec![vec![1, 2, 3, 4, 5]]);
        // sum is not high enough
        let mut unique_support = BTreeMap::new();
        let mut support_length = BTreeMap::new();
        unique_support.entry(vec![1, 2, 3, 4, 5]).or_insert(vec![
            vec![3, 4, 5],
            vec![3, 4, 5],
            vec![4, 5],
            vec![4, 5],
            vec![4, 5],
            vec![4, 5],
        ]);
        unique_support.entry(vec![1, 2, 3, 4, 6]).or_insert(vec![
            vec![4, 6],
            vec![4, 6],
            vec![4, 6],
        ]);
        // sum 14 vs. 6
        // longest 3 vs. 2
        support_length.entry(vec![1, 2, 3, 4, 5]).or_insert(vec![
            (vec![3 as i32, 4 as i32, 5 as i32], 3 as usize),
            (vec![3 as i32, 4 as i32, 5 as i32], 3 as usize),
            (vec![4 as i32, 5 as i32], 2 as usize),
            (vec![4 as i32, 5 as i32], 2 as usize),
            (vec![4 as i32, 5 as i32], 2 as usize),
            (vec![4 as i32, 5 as i32], 2 as usize),
        ]);
        support_length.entry(vec![1, 2, 3, 4, 6]).or_insert(vec![
            (vec![4 as i32, 6 as i32], 2 as usize),
            (vec![4 as i32, 6 as i32], 2 as usize),
            (vec![4 as i32, 6 as i32], 2 as usize),
        ]);
        let candidates = filter_assembled_candidates(unique_support, support_length).unwrap();
        assert_eq!(candidates, vec![vec![1, 2, 3, 4, 5], vec![1, 2, 3, 4, 6]]);
        // only one read each
        let mut unique_support = BTreeMap::new();
        let mut support_length = BTreeMap::new();
        unique_support
            .entry(vec![1, 2, 3, 4, 5])
            .or_insert(vec![vec![2, 3, 4, 5]]);
        unique_support
            .entry(vec![1, 2, 3, 4, 6])
            .or_insert(vec![vec![4, 6]]);
        // longest 4 vs. 2
        support_length.entry(vec![1, 2, 3, 4, 5]).or_insert(vec![(
            vec![2 as i32, 3 as i32, 4 as i32, 5 as i32],
            4 as usize,
        )]);
        support_length
            .entry(vec![1, 2, 3, 4, 6])
            .or_insert(vec![(vec![4 as i32, 6 as i32], 2 as usize)]);
        let candidates = filter_assembled_candidates(unique_support, support_length).unwrap();
        assert_eq!(candidates, vec![vec![1, 2, 3, 4, 5]]);

        // only one or two reads each
        let mut unique_support = BTreeMap::new();
        let mut support_length = BTreeMap::new();
        unique_support
            .entry(vec![1, 2, 3, 4, 5])
            .or_insert(vec![vec![2, 3, 4, 5], vec![3, 4, 5]]);
        unique_support
            .entry(vec![1, 2, 3, 4, 6])
            .or_insert(vec![vec![4, 6], vec![4, 6]]);
        // longest 4 vs. 2
        // sum 7 vs. 4
        support_length.entry(vec![1, 2, 3, 4, 5]).or_insert(vec![
            (vec![2 as i32, 3 as i32, 4 as i32, 5 as i32], 4 as usize),
            (vec![3 as i32, 4 as i32, 5 as i32], 3 as usize),
        ]);
        support_length.entry(vec![1, 2, 3, 4, 6]).or_insert(vec![
            (vec![4 as i32, 6 as i32], 2 as usize),
            (vec![4 as i32, 6 as i32], 2 as usize),
        ]);
        let candidates = filter_assembled_candidates(unique_support, support_length).unwrap();
        assert_eq!(candidates, vec![vec![1, 2, 3, 4, 5]]);
        // sum is not high enough
        let mut unique_support = BTreeMap::new();
        let mut support_length = BTreeMap::new();
        unique_support
            .entry(vec![1, 2, 3, 4, 5])
            .or_insert(vec![vec![2, 3, 4, 5], vec![4, 5]]);
        unique_support
            .entry(vec![1, 2, 3, 4, 6])
            .or_insert(vec![vec![4, 6], vec![4, 6]]);
        // longest 4 vs. 2
        // sum 6 vs. 4
        support_length.entry(vec![1, 2, 3, 4, 5]).or_insert(vec![
            (vec![2 as i32, 3 as i32, 4 as i32, 5 as i32], 4 as usize),
            (vec![4 as i32, 5 as i32], 2 as usize),
        ]);
        support_length.entry(vec![1, 2, 3, 4, 6]).or_insert(vec![
            (vec![4 as i32, 6 as i32], 2 as usize),
            (vec![4 as i32, 6 as i32], 2 as usize),
        ]);
        let candidates = filter_assembled_candidates(unique_support, support_length).unwrap();
        assert_eq!(candidates, vec![vec![1, 2, 3, 4, 5], vec![1, 2, 3, 4, 6]]);
    }

    #[test]
    fn test_pick_supported_hap_with_seed() {
        let supported_haps = vec![vec![4, 5, 6], vec![1, 2, 3], vec![7, 8, 9]];
        let mut rng = rand::rngs::SmallRng::seed_from_u64(1);
        let first = pick_supported_hap_with_seed(&mut rng, &supported_haps).unwrap();

        let mut rng = rand::rngs::SmallRng::seed_from_u64(1);
        let second = pick_supported_hap_with_seed(&mut rng, &supported_haps).unwrap();
        assert_eq!(first, second);

        let reordered_haps = vec![vec![7, 8, 9], vec![4, 5, 6], vec![1, 2, 3]];
        let mut rng = rand::rngs::SmallRng::seed_from_u64(1);
        let reordered_pick = pick_supported_hap_with_seed(&mut rng, &reordered_haps).unwrap();
        assert_eq!(first, reordered_pick);
    }
}
