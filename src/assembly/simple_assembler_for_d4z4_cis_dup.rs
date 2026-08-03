use crate::assembly::assembler::{AssemblyResult, FpGraph};
use crate::util::DError;
use log::debug;
use std::collections::BTreeMap;

impl FpGraph {
    /// Simple run for D4Z4 cis duplications
    /// # Returns
    /// * `AssemblyResult` - assembly result
    pub fn run_simple(&mut self) -> Result<AssemblyResult, DError> {
        self.get_edges();
        self.analyze_nodes(false);
        self.clean_reads = self.reads.clone();
        for (node1, node2) in self.edges.keys() {
            self.nodes.insert(*node1);
            self.nodes.insert(*node2);
        }
        let mut starts = Vec::new();
        let mut ends = Vec::new();
        for node in &self.nodes {
            let this_node_prev = self.previous_per_node.get(node);
            let this_node_next = self.next_per_node.get(node);
            if this_node_prev.is_none() {
                starts.push(*node);
            }
            if this_node_next.is_none() {
                ends.push(*node);
            }
        }
        debug!("starting_nodes {starts:?} ending_nodes {ends:?}");
        let mut assembled_haps: Vec<Vec<i32>> = Vec::new();
        // forward
        let mut incomplete_haps_forward: Vec<Vec<i32>> = Vec::new();
        for starting_node in &starts {
            if let Some(starting_next_nodes) = self.next_per_node.get(starting_node) {
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
                        let extended_haps = self.assemble_next_simple(this_hap.clone())?;
                        debug!("extended {:?} to {:?}", this_hap, extended_haps);
                        if extended_haps.is_empty() {
                            if !incomplete_haps_forward.contains(&this_hap) {
                                incomplete_haps_forward.push(this_hap);
                            }
                        } else {
                            for each_extended_hap in &extended_haps {
                                all_extended.push(each_extended_hap.to_vec());
                            }
                        }
                    }
                    assembled_haps = all_extended
                        .iter()
                        .filter(|x| !incomplete_haps_forward.contains(x))
                        .map(|x| x.to_vec())
                        .collect::<Vec<_>>();
                    if assembled_haps.is_empty() {
                        break;
                    }
                    if nstep > 10 || assembled_haps.len() > 20 {
                        for a in &assembled_haps {
                            incomplete_haps_forward.push(a.to_vec());
                        }
                        debug!("stopping assembly at step {nstep}");
                        break;
                    }
                }
            }
        }
        debug!("forward incomplete: {:?}", incomplete_haps_forward);

        // backward
        let mut assembled_haps: Vec<Vec<i32>> = Vec::new();
        let mut incomplete_haps_backward: Vec<Vec<i32>> = Vec::new();
        for ending_node in &ends {
            if let Some(ending_next_nodes) = self.previous_per_node.get(ending_node) {
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
                        let extended_haps = self.assemble_prev_simple(this_hap.clone())?;
                        debug!("extended {:?} to {:?}", this_hap, extended_haps);
                        if extended_haps.is_empty() {
                            if !incomplete_haps_backward.contains(&this_hap) {
                                incomplete_haps_backward.push(this_hap);
                            }
                        } else {
                            for each_extended_hap in &extended_haps {
                                all_extended.push(each_extended_hap.to_vec());
                            }
                        }
                    }
                    assembled_haps = all_extended
                        .iter()
                        .filter(|x| !incomplete_haps_backward.contains(x))
                        .map(|x| x.to_vec())
                        .collect::<Vec<_>>();
                    if assembled_haps.is_empty() {
                        break;
                    }
                    if nstep > 10 || assembled_haps.len() > 20 {
                        for a in &assembled_haps {
                            incomplete_haps_backward.push(a.to_vec());
                        }
                        debug!("stopping assembly at step {nstep}");
                        break;
                    }
                }
            }
        }
        debug!("backward incomplete: {:?}", incomplete_haps_backward);

        let mut incomplete_haps = Vec::new();
        for hap in &incomplete_haps_forward {
            if !incomplete_haps.contains(hap) {
                incomplete_haps.push(hap.to_vec());
            }
        }
        for hap in &incomplete_haps_backward {
            if !incomplete_haps.contains(hap) {
                incomplete_haps.push(hap.to_vec());
            }
        }

        let mut nstep = 0;
        loop {
            nstep += 1;
            let nodes_not_used = self.get_unused_nodes(incomplete_haps.clone());
            debug!("nodes not used: {nodes_not_used:?}");

            if nodes_not_used.is_empty() || nstep > 5 {
                break;
            }
            let Some(node) = nodes_not_used.first() else {
                break;
            };
            debug!("checking unused node {node}");
            let extended_haps = self.assemble_next_simple(vec![*node])?;
            if let Some(this_node_extended) = extended_haps.first() {
                if !incomplete_haps.contains(this_node_extended) {
                    incomplete_haps.push(this_node_extended.to_vec());
                    debug!("adding {this_node_extended:?} to incomplete_haps");
                }
            } else {
                let extended_haps = self.assemble_prev_simple(vec![*node])?;
                if let Some(this_node_extended) = extended_haps.first() {
                    if !incomplete_haps.contains(this_node_extended) {
                        incomplete_haps.push(this_node_extended.to_vec());
                        debug!("adding {this_node_extended:?} to incomplete_haps");
                    }
                }
            }
        }

        Ok(AssemblyResult {
            complete: vec![],
            incomplete: incomplete_haps,
            special_incomplete: vec![],
            supporting_reads: BTreeMap::new(),
            nonunique_reads: vec![],
        })
    }

    /// Extend given haplotype by one node, forward, only allows one-to-one match
    /// # Arguments
    /// * `this_hap` - haplotype
    /// # Returns
    /// * `Vec<Vec<i32>>` - extended haplotypes
    fn assemble_next_simple(&self, this_hap: Vec<i32>) -> Result<Vec<Vec<i32>>, DError> {
        let last_unit = this_hap.last().ok_or("last not found")?;
        if !self.next_per_node.contains_key(last_unit) {
            return Ok(vec![]);
        }
        let next_nodes = &self
            .next_per_node
            .get(last_unit)
            .ok_or("key not found in next_per_node")?
            .iter()
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
            return Ok(haps_candidates);
        }
        return Ok(vec![]);
    }

    /// Extend given haplotype by one node, going backward, only allows one-to-one match
    /// # Arguments
    /// * `this_hap` - haplotype
    /// # Returns
    /// * `Vec<Vec<i32>>` - extended haplotypes
    fn assemble_prev_simple(&self, this_hap: Vec<i32>) -> Result<Vec<Vec<i32>>, DError> {
        let first_unit = this_hap.first().ok_or("last not found")?;
        if !self.previous_per_node.contains_key(first_unit) {
            return Ok(vec![]);
        }
        let prev_nodes = &self
            .previous_per_node
            .get(first_unit)
            .ok_or("key not found in previous_per_node")?
            .iter()
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
        if haps_candidates.len() == 1 {
            return Ok(haps_candidates);
        }
        return Ok(vec![]);
    }
}
