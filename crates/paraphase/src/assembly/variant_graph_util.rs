use crate::assembly::assembly_result::AssembledPaths;
use crate::assembly::node_datum::NodeDatum;
use crate::assembly::variant_graph::{Graph, IxType, ReadHapPair, SegmentationClass};
use crate::detail::hapcmp::HapCompare;
use crate::detail::util::{DError, HashMap};
use crate::io::json::ReadFingerprintMap;
use vstr::{VStr, VString};

use itertools::Itertools;
use petgraph::graph::EdgeIndex;

use std::collections::{BTreeMap, BTreeSet};

/// Compute the midpoints of an array.
/// Returns None if the array is empty.
/// If the array length is odd, returns two references to the same location.
#[must_use]
pub fn midpoints<T>(x: &[T]) -> Option<(&T, &T)> {
    let len = x.len();
    match len {
        0 => None,
        _ => {
            let midpoint = len / 2;
            Some((&x[midpoint], &x[midpoint - 1]))
        }
    }
}

/// Compute the median of a slice of integers.
/// Assigns a float, averages midpoints in case of an even-sized slice.
/// ```
/// use paraphase::assembly::variant_graph_util::median_i32;
/// let x = &[0, 3i32, 4];
/// assert_eq!(3., median_i32(x));
/// assert!(median_i32(&[] as &[i32]).is_nan());
/// let x = &[4i32, 2];
/// assert_eq!(median_i32(x), 3.);
/// let x = &[1i32, 2];
/// assert_eq!(median_i32(x), 1.5f32);
/// ```
#[must_use]
pub fn median_i32(x: &[i32]) -> f32 {
    if x.is_empty() {
        return f32::NAN;
    }
    let copied = x.iter().copied().sorted().collect::<Vec<_>>();
    let len = copied.len();
    let midpoint = len / 2;
    if len & 1 != 0 {
        copied[midpoint] as f32
    } else {
        (copied[midpoint] + copied[midpoint - 1]) as f32 * 0.5
    }
}

/// Compute percentile over an array.
/// Matches np.percentile with method = "linear"
///```
/// use paraphase::assembly::variant_graph_util::percentile_i32;
/// let a: &[i32] = &[10, 7, 4, 3, 2, 1];
/// assert!((percentile_i32(a, 50) - 3.5).abs() < 1e-10);
/// let a: &[i32] = &[
///     41, 20, 17, 48, 11, 16, 2, 21, 26, 40, 12, 7, 48, 14, 28, 44, 37, 23, 15, 41, 24, 18, 49,
///     39, 35, 18, 0, 22, 19, 34, 39, 29, 4, 10, 4, 11, 10, 30, 26, 33, 14, 17, 35, 9, 32, 47, 28,
///     11, 19, 30, 11, 35, 39, 42, 39, 37, 26, 47, 49, 1, 31, 35, 22, 11, 3, 25, 30, 9, 44, 24,
///     19, 14, 9, 35, 27, 36, 16, 10, 20, 44, 14, 6, 42, 13, 23, 36, 3, 43, 4, 8, 13, 31, 25, 10,
///     41, 20, 28, 47, 47, 9,
/// ];
/// assert!(
///     (percentile_i32(a, 80) - 39.0).abs() < 1e-10,
///     "{} vs expected 39.0",
///     percentile_i32(a, 80)
/// );
/// assert!(
///     (percentile_i32(a, 60) - 28.4).abs() < 1e-10,
///     "{} vs expected 28.4",
///     percentile_i32(a, 60)
/// );
/// ```
#[must_use]
pub fn percentile_i32(x: &[i32], percentile: i32) -> f64 {
    if x.is_empty() {
        return f64::NAN;
    }
    let copied = x.iter().copied().sorted().collect::<Vec<_>>();
    let len = copied.len();
    const ALPHA: i32 = 1;
    const BETA: i32 = 1;
    let percentile = f64::from(percentile) * 0.01;
    // Numpy virtual index: (q / 100) * ( n - ALPHA - BETA + 1 ) + ALPHA
    // let n_mul = (len as i32 - 1) as f64;
    let n_mul = f64::from(len as i32 - ALPHA - BETA + 1);

    #[cfg(feature = "nightly")]
    let virtual_index = fma::fma(percentile, n_mul, ALPHA as f64);

    #[cfg(not(feature = "nightly"))]
    let virtual_index = percentile * n_mul + f64::from(ALPHA);

    let accessed_index = (virtual_index as usize) - 1;
    let fract = virtual_index.fract();
    let accessed_val = f64::from(copied[accessed_index]);
    if fract != 0.0 {
        #[cfg(feature = "nightly")]
        {
            fma::fma(
                accessed_val,
                1. - fract,
                copied[accessed_index + 1] as f64 * fract,
            )
        }
        #[cfg(not(feature = "nightly"))]
        {
            accessed_val * (1. - fract) + f64::from(copied[accessed_index + 1]) * fract
        }
        // Convex combination
        // Equivalent to: accessed_val * (1. - fract) + copied[accessed_index + 1] as f64 * fract
        // but better precision
    } else {
        accessed_val
    }
}

/// Find last index that is not 'x', or 0 if no such base exists.
/// ```
/// paraphase::detail::util::init_log(log::LevelFilter::Debug);
/// use paraphase::assembly::variant_graph_util::count_suffix_x;
/// assert_eq!(count_suffix_x(&[0, 1, b'x']), 1, "0, 1, x");
/// assert_eq!(count_suffix_x(&[0, 1]), 1, "0, 1");
/// assert_eq!(count_suffix_x(&[0, 1, b'x', 1]), 3, "0x1x");
/// assert_eq!(count_suffix_x(b"xxx122"), 5, "xxx122");
/// assert_eq!(count_suffix_x(b"xxxxxx"), 0, "xxxxxx");
/// assert_eq!(count_suffix_x(b"1xxxxx"), 0, "1xxxxx");
/// assert_eq!(count_suffix_x(b"1xxxxx1"), 6, "1xxxxx1");
/// ```
#[must_use]
pub fn count_suffix_x(x: &[u8]) -> usize {
    x.iter().rposition(|item| *item != b'x').unwrap_or(0)
}

impl Graph {
    /// `paraphase_format_node`: format the paraph-rs format (stored in graph) to the string representation used in paraphase.
    /// This is mostly for readability and testing.
    #[must_use]
    pub fn paraphase_format_node(
        &self,
        node_index: impl TryInto<IxType>,
        node_id_lookup: Option<&HashMap<u32, &NodeDatum>>,
    ) -> Option<String> {
        let node_index = node_index.try_into().ok()?;
        // .unwrap_or_else(|_| panic!("Node type not convertible to integral IxType"));
        let node_id_local = if node_id_lookup.is_none() {
            Some(self.get_id_lookup())
        } else {
            None
        };
        let node_id_lookup = node_id_lookup.unwrap_or(node_id_local.as_ref().unwrap());
        node_id_lookup
            .get(&node_index)
            .map(|node| format!("{}-{}", node.hap, node.pos))
    }

    /// `paraphase_format_edge`: format the paraph-rs format (stored in graph) to the string representation used in paraphase.
    /// This is mostly for readability and testing.
    #[must_use]
    pub fn paraphase_format_edge(
        &self,
        edge: EdgeIndex<IxType>,
        node_id_lookup: Option<&HashMap<u32, &NodeDatum>>,
    ) -> Option<(String, String)> {
        self.graph
            .edge_endpoints(edge)
            .map(|(source, target)| {
                (
                    self.paraphase_format_node(source.index(), node_id_lookup),
                    self.paraphase_format_node(target.index(), node_id_lookup),
                )
            })
            .and_then(|x| match x {
                (Some(y), Some(z)) => Some((y, z)),
                _ => None,
            })
    }

    /// Given strings to the left and right and connections between sub-haplotypes,
    /// generate additional possible complete haplotypes.
    #[must_use]
    pub(crate) fn rescue_missing<'a>(
        left: &[&'a VString],
        right: &[&'a VString],
        dnext: &BTreeMap<VString, BTreeSet<VString>>,
        dbefore: &BTreeMap<VString, BTreeSet<VString>>,
    ) -> BTreeSet<VString> {
        let sub_rescue = |side: &[&VString],
                          other_side: &[&VString],
                          match_dict: &BTreeMap<VString, BTreeSet<VString>>,
                          is_left: bool|
         -> BTreeSet<VString> {
            let merge_slices_base = |x: &VString, y: &VString| -> VString {
                let mut ret = x.clone();
                ret.extend_from_slice(&y[..]);
                ret
            };
            let merge_slices = |x: &VString, y: &VString, is_left: bool| -> VString {
                if is_left {
                    merge_slices_base(x, y)
                } else {
                    merge_slices_base(y, x)
                }
            };
            side.iter()
                .filter_map(|x| match_dict.get(*x).map(|matches| (x, matches)))
                .filter(|x| bytecount::count(x.0, b'x') <= 1)
                .filter_map(|(side_missing, next_hits)| {
                    if next_hits.len() == 1
                        && (other_side.is_empty() || other_side[..] == [next_hits.iter().next().unwrap()])
                    {
                        Some(
                            merge_slices(side_missing, next_hits.iter().next().unwrap(), is_left)
                        )
                    } else {
                        let cand_missing = other_side
                            .iter()
                            .filter(|x| next_hits.contains(**x))
                            .collect::<Vec<_>>();
                        if cand_missing.len() == 1
                            && other_side.len() == 1
                            && &other_side[0] == cand_missing[0]
                            && side.len() == 1
                        {
                            log::trace!(
                        "Merging {side_missing:?}. {next_hits:?} and {cand_missing:?} from {}", if is_left {"left"} else {"right"}
                    );
                            Some(merge_slices(side_missing, cand_missing[0], is_left))
                        } else {
                            None
                        }
                }
                }).collect::<BTreeSet<_>>()
        };
        let mut ret = sub_rescue(left, right, dnext, /* is_left= */ true);
        for item in sub_rescue(right, left, dbefore, false) {
            ret.insert(item);
        }
        log::trace!("After rescuing: {ret:?}");
        ret
    }

    /// Find potential partial haplotypes that match haps at left/right
    /// to use for assembling.
    #[must_use]
    pub(crate) fn get_missing<'a>(
        sub_haps_assembled: &BTreeSet<VString>,
        left: &'a [VString],
        right: &'a [VString],
    ) -> (Vec<&'a VString>, Vec<&'a VString>) {
        let left = left
            .iter()
            .filter(|x| {
                let len = x.len();
                !sub_haps_assembled
                    .iter()
                    .any(move |sub| sub.len() >= len && sub[..len] == x[..])
            })
            .collect::<Vec<&VString>>();
        let right = right
            .iter()
            .filter(|x| {
                let len = x.len();
                !sub_haps_assembled
                    .iter()
                    .any(move |sub| sub.len() >= len && sub[sub.len() - len..] == x[..])
            })
            .collect::<Vec<&VString>>();
        log::trace!("Getting missing for left {left:?}/right {right:?} use sub haps {sub_haps_assembled:?}. Found: {left:?}, {right:?}");
        (left, right)
    }

    /// Merge blocks together in order of decreasing size.
    pub fn merge_blocks_by_size(&mut self) -> Result<(), DError> {
        let mut iternum = 0usize;
        log::debug!(
            "[merge_blocks_by_size] self.pos_edge_map {:?}",
            self.pos_edge_map
        );
        while !self.pos_edge_map.is_empty() {
            iternum += 1;
            let order = self.order_blocks_by_size();
            log::trace!("Order: {order:?} at iteration {iternum} for merge_blocks_by_size");
            if order.len() == 1 {
                log::trace!(
                    "Only one block remains. Break!. Order: {order:?}. pem: {:?}",
                    self.formatted_edges()
                );
                break;
            }
            let mut i = 0;
            let mut total_success = false;
            for region in &order {
                let pos1 = region[0];
                let pos2 = region[1];
                i += 1;
                log::trace!("About to assemble subregion {pos1}..={pos2}");
                let (success, haps, _x, _y) = self.subregion_assembly(pos1, pos2, false)?;
                log::trace!("Assembled subregion {pos1}..={pos2}");
                if success {
                    log::trace!(
                        "Success for {pos1}..={pos2} yielded haps {haps:?} at iternum {iternum}. Current state before rm_add_edges {pos1}..={pos2}: {:?}", self.formatted_edges()
                    );
                    self.rm_add_edges(pos1, pos2, haps.iter())?;
                    log::trace!(
                        "Success for {pos1}..={pos2} yielded haps {haps:?} at iternum {iternum}. Current state after rm_add_edges {pos1}..={pos2}: {:?}", self.formatted_edges()
                    );
                    total_success = true;
                    break;
                }
                log::trace!("At iternum {iternum}, subregion assembly for {pos1}..={pos2} did not have success.");
            }
            if i == order.len() && !total_success {
                break;
            }
            log::trace!(
                "Success: {total_success}. Remaining in pos edge map: {:?} at iternum {iternum}",
                self.pos_edge_map.len()
            );
        }
        Ok(())
    }

    /// Sort blocks by size.
    /// Returns sorted vector of (node1, node2, minlen, maxlen)
    /// ordered by increasing (minlen, maxlen)
    #[must_use]
    pub(crate) fn order_blocks_by_size(&self) -> Vec<[u32; 4]> {
        let lengths = self
            .node_iter()
            .map(|x| (x.pos, x.hap.len() as u32))
            .collect::<BTreeSet<_>>()
            .into_iter()
            .collect::<Vec<_>>();
        lengths
            .windows(2)
            .map(|window| {
                let len1 = window[0].1;
                let len2 = window[1].1;
                let node1 = window[0].0;
                let node2 = window[1].0;
                let minlen = std::cmp::min(len1, len2);
                let maxlen = std::cmp::max(len1, len2);
                [node1, node2, minlen, maxlen]
            })
            .sorted_by_key(|x| (-(x[2] as i32), -(x[3] as i32)))
            .collect::<Vec<_>>()
    }

    /// Return the highest possible copy number given the assembled haplotypes
    /// Returns (highest count, `candidate_hap_sets`)
    pub(crate) fn get_highest_cn<T>(
        &self,
        paths: impl IntoIterator<Item = T>,
    ) -> Result<(usize, Vec<Vec<VString>>), DError>
    where
        T: std::convert::Into<VString>,
    {
        let seqs = paths
            .into_iter()
            .map(std::convert::Into::into)
            .collect::<Vec<_>>();
        let mut cand = Vec::with_capacity(self.nvar / 2);
        for set in (0..self.nvar)
            .rev()
            .map(|idx| seqs.iter().filter(|x| x[idx] != b'x').collect::<Vec<_>>())
        {
            if !cand.contains(&set) {
                cand.push(set);
            }
        }
        let cand = cand
            .into_iter()
            .sorted_by_key(|x| {
                (
                    std::cmp::Reverse(x.len()),
                    x.into_iter()
                        .map(|hap| hap.iter().filter(|&n| *n == b'x').collect::<Vec<_>>().len())
                        .sum::<usize>(),
                )
            })
            .map(|x| x.into_iter().map(VString::from).collect::<Vec<_>>())
            .collect::<Vec<_>>();
        log::trace!("[get_highest_cn] candidate haplotype sets: {cand:?}");
        // let max_len = cand.first().map(|x| x.len());
        Ok(cand
            .first()
            .map(std::vec::Vec::len)
            .map(|x| (x, cand))
            .ok_or_else(|| {
                anyhow::anyhow!("get_highest_cn failed because there were no passing haps")
            })?)
    }

    /// Choose the best candidates from an input set.
    #[must_use]
    fn filter_candidates(
        &self,
        hap: VStr<'_>,
        candidates: &BTreeSet<VString>,
        len_x: usize,
        candidate_first: bool,
    ) -> Vec<VString> {
        let mut read_support = BTreeMap::<VStr<'_>, Vec<i32>>::new();
        let mut ret = vec![];
        for c in candidates {
            let len = c.len();
            let mut hap2 = VString::from(vec![
                b'x';
                if candidate_first {
                    len_x - len
                } else {
                    self.nvar - len_x
                }
            ]);
            hap2.extend_from_slice(&c[..]);
            hap2.resize(self.nvar, b'x');
            let hap2 = hap2; // make const
            for read_seq in self.reads.values() {
                let cmp1 = HapCompare::from_haps(hap, read_seq).unwrap();
                let cmp2 = HapCompare::from_haps(&hap2, read_seq).unwrap();
                if cmp1.mismatches == 0
                    && cmp2.mismatches == 0
                    && cmp1.matches > 0
                    && cmp2.matches > 0
                {
                    read_support
                        .entry(c.vstr())
                        .or_default()
                        .push(std::cmp::min(cmp1.matches, cmp2.matches));
                }
            }
        }
        log::trace!("read_support: {read_support:?}");
        match read_support.len() {
            0 => {
                return ret;
            }
            1 => {
                ret.push(read_support.iter().next().unwrap().0.into());
                return ret;
            }
            _ => {}
        }
        let medians_and_maxes = read_support
            .values()
            .map(|x| {
                (
                    median_i32(x),
                    x.iter().copied().max().map_or(f32::NAN, |x| x as f32),
                )
            })
            .collect::<Vec<(f32, f32)>>();
        for ((max_other, median, id), (hap, read_ids)) in medians_and_maxes
            .iter()
            .map(|x| x.0)
            .enumerate()
            .map(|(id, median)| {
                let max_other = medians_and_maxes
                    .iter()
                    .map(|(_median, max)| max)
                    .copied()
                    .enumerate()
                    .filter_map(
                        |(other_id, max)| {
                            if other_id == id {
                                None
                            } else {
                                Some(max)
                            }
                        },
                    )
                    .max_by(|x, y| x.partial_cmp(y).unwrap_or(std::cmp::Ordering::Equal))
                    .unwrap();
                (max_other, median, id)
            })
            .zip(read_support.iter())
        {
            log::trace!("median: {median}. max other: {max_other} for id {id}");
            if median > max_other && median >= 4. && read_ids.len() >= 4 {
                ret.push(hap.into());
                break;
            }
        }
        if ret.is_empty() {
            ret = candidates.iter().cloned().collect::<Vec<_>>();
        }
        ret
    }

    /// If x has no 'x' charactes, simply return it.
    /// Otherwise, extend.
    fn extended_hap(&mut self, hap: VStr<'_>) -> Result<VString, DError> {
        if !hap.contains(&b'x') {
            return Ok(VString::from(hap));
        }
        let prefix_x_count = hap.iter().position(|x| *x != b'x').unwrap_or(hap.len());
        log::trace!("Hap: {hap} has length {}", hap.len());
        let suffix_x_count = count_suffix_x(&hap);
        log::trace!("pref x count: {prefix_x_count}. suffix: {suffix_x_count} for hap {hap}");
        let mut identical_bases_right = VString::default();
        let mut identical_bases_left = VString::default();
        // extend right
        if suffix_x_count < (self.nvar - 2) {
            log::trace!(
                "x count {suffix_x_count} < nvar - 2 ({}) for {hap}. dnhap: {:?}",
                self.nvar - 2,
                self.dnhap
            );
            let next_pos = (suffix_x_count + 1) as u32;
            // check that this is not a site very distant from next site
            let dnhap_is_present_and_assigned = matches!(
                self.dnhap
                    .get(&(suffix_x_count as u32))
                    .map(SegmentationClass::is_assigned),
                Some(true)
            );
            let min_support = self
                .edge_info
                .get(suffix_x_count)
                .and_then(|x| x.iter().copied().min())
                .unwrap_or(0);
            log::trace!(
                "dnhap_is_present_and_assigned on right for {hap}. pos edge map has minimum {min_support}. Pos edge map at {suffix_x_count}", 
            );
            if self.pos_edge_map.contains_key(&suffix_x_count) {
                log::trace!(
                    "[extended_hap] self.pos_edge_map[&suffix_x_count] {:?}",
                    self.pos_edge_map[&suffix_x_count]
                );
            }
            if dnhap_is_present_and_assigned && min_support >= 2 {
                log::trace!("dnhap had the suffix, it was assigned, and pos_edge_map had at least min count of 2 for hap {hap}.");
                let previous_pos = self.get_previous_pos(next_pos).ok_or_else(|| anyhow::anyhow!("Failed to get previous pos for {next_pos} with suffix x count {suffix_x_count}"))?;
                log::trace!("Assembling subregion from {previous_pos}..={next_pos}");
                let (_, _, dnext, _dbefore) =
                    self.subregion_assembly(previous_pos, next_pos, false)?;
                log::trace!("dnext {dnext:?} and dbefore {_dbefore:?}");
                let block = VString::from(&hap[(previous_pos as usize)..=suffix_x_count]);
                if let Some(candidates) = dnext.get(&block).map(|candidates| {
                    log::trace!("candidates from dnext: {candidates:?}");
                    if candidates.len() <= 1 {
                        candidates.iter().cloned().collect::<Vec<_>>()
                    } else {
                        log::trace!("filtering!");
                        self.filter_candidates(
                            hap,
                            candidates,
                            self.nvar - 1 - suffix_x_count,
                            false,
                        )
                    }
                }) {
                    log::trace!("filtered candidates : {candidates:?}");
                    if let Some(next_hap_len) = candidates.first().map(|x| x.len()) {
                        for idx in 0..next_hap_len {
                            let bases_per_candidate =
                                candidates.iter().map(|x| x[idx]).collect::<BTreeSet<_>>();
                            if bases_per_candidate.len() != 1 {
                                break;
                            }
                            identical_bases_right.push(*bases_per_candidate.first().unwrap());
                        }
                    }
                }
            } else {
                log::trace!(
                    "Did not find the dnhap match for {suffix_x_count} as key and dnhap {:?}",
                    self.dnhap
                );
            }
        }
        // extend right
        else {
            log::trace!("No extend right for {hap}");
        }
        // extend left
        if prefix_x_count >= 2 {
            let next_pos = prefix_x_count;
            let nend = self.get_nstart_nend().expect("Missing nend").1;
            log::trace!("next pos: {next_pos}, nend = {nend}");
            let hap_block = VString::from(if next_pos < nend as usize {
                let next_next_pos = self
                    .get_next_pos(next_pos as u32)
                    .expect("Missing next next pos");
                &hap[prefix_x_count..next_next_pos as usize]
            } else {
                &hap[prefix_x_count..]
            });
            // # check that this is not a site very distant from previous site
            let dnhap_is_assigned = self
                .dnhap
                .get(&((prefix_x_count - 1) as u32))
                .map(SegmentationClass::is_assigned)
                == Some(true);
            let passing_edge_support = self
                .edge_info
                .get(prefix_x_count - 1)
                .and_then(|x| x.iter().copied().min())
                .unwrap_or(0);
            let extends_left = dnhap_is_assigned && passing_edge_support >= 2;
            log::trace!("Extending hap {hap:?} left? {extends_left}. dnhap: {dnhap_is_assigned}. Edge support: {passing_edge_support} when querying index {}. Edges: {:?}", (prefix_x_count - 1) as u32, self.edge_info.get(prefix_x_count - 1));
            if extends_left {
                let previous_pos = self
                    .get_previous_pos(next_pos as u32)
                    .ok_or_else(|| anyhow::anyhow!("prev pos missing in extend_left {next_pos}"))?;
                let dbefore = self
                    .subregion_assembly(previous_pos, next_pos as u32, /* allow_x */ false)?
                    .3;
                if let Some(candidates) = dbefore.get(&hap_block).map(|candidates| {
                    if candidates.len() <= 1 {
                        candidates.iter().cloned().collect::<Vec<_>>()
                    } else {
                        self.filter_candidates(hap, candidates, prefix_x_count, true)
                    }
                }) {
                    if !candidates.is_empty() {
                        let previous_hap_len = candidates.first().map(|x| x.len()).unwrap();
                        for idx in (0..previous_hap_len).rev() {
                            let bases_per_candidate =
                                candidates.iter().map(|x| x[idx]).collect::<BTreeSet<_>>();
                            if bases_per_candidate.len() != 1 {
                                break;
                            }
                            identical_bases_left.push(*bases_per_candidate.first().unwrap());
                        }
                        identical_bases_left.reverse();
                    }
                }
            }
        } else {
            log::trace!("No extend left");
        }
        log::trace!("from left: {identical_bases_left}. from right: {identical_bases_right}. pre count {prefix_x_count}. suf {suffix_x_count} for hap {hap}");
        let mut ret = VString::default();
        ret.append(&mut vec![b'x'; prefix_x_count - identical_bases_left.len()]);
        ret.extend_from_slice(&identical_bases_left);
        ret.extend_from_slice(&hap[prefix_x_count..=suffix_x_count]);
        ret.extend_from_slice(&identical_bases_right);
        ret.resize(self.nvar, b'x');
        log::trace!("Starting with {hap}. from left: {identical_bases_left}. from right: {identical_bases_right}. Final ret: {ret}");
        Ok(ret)
    }

    /// Extend haplotypes from pivot blocks.
    /// Takes a set of haplotypes and extends them with `extended_hap`.
    #[must_use]
    pub(crate) fn extend_pivot_blocks<'a, T: 'a>(
        &mut self,
        x: impl IntoIterator<Item = T>,
    ) -> AssembledPaths
    where
        T: std::convert::Into<VStr<'a>>,
    {
        AssembledPaths::from_seqs(x.into_iter().map(|x| {
            let hap = x.into();
            let res = self
                .extended_hap(hap)
                .unwrap_or_else(|e| panic!("Failed to extend hap {hap}. Error: {e:?}"));
            log::trace!("Original hap {hap:?} was extended to {res:?}");
            res
        }))
    }

    /// Filter low-support haplotypes.
    /// Takes an `AssembledPaths` object and returns a `Result<AssembledPaths, DError>`
    /// where `DError` is an alias for `Box<dyn std::error::Error>`.
    pub fn filter_low_support_haps(
        &self,
        init_haps: &AssembledPaths,
    ) -> Result<AssembledPaths, DError> {
        Self::filter_low_support_haps_detail(
            init_haps,
            self.settings.min_hap_support as usize,
            &self.reads_original,
        )
    }

    pub fn filter_low_support_haps_detail(
        init_haps: &AssembledPaths,
        min_count: usize,
        reads: &ReadFingerprintMap,
    ) -> Result<AssembledPaths, DError> {
        log::trace!("Beginning filter_low_support_haps!!! min_count: {min_count}. Haps: {init_haps:?}. reads: {reads:?}");
        let flat_init_haps = init_haps
            .iter()
            .map(|x| x.vstr())
            .collect::<Vec<VStr<'_>>>();
        let read_hap_support =
            Self::match_reads_and_haps(reads, &flat_init_haps[..], /* min_match= */ None)?;
        let good_reads = read_hap_support.support_reads;
        log::trace!("good reads, initial: {good_reads:?}");
        let mut filtered_ass_haps = good_reads
            .into_iter()
            .filter(|(_id, support)| !support.is_empty())
            .collect::<Vec<_>>();
        log::trace!(
            "{} passing haps after {} initial. Min count {min_count}",
            filtered_ass_haps.len(),
            flat_init_haps.len()
        );
        let mut iternum = 0;
        loop {
            let haps_to_assess = filtered_ass_haps
                .iter()
                .map(|x| flat_init_haps[x.0 as usize])
                .collect::<Vec<_>>();
            log::trace!("At iternum {iternum}, we have {haps_to_assess:?} as haps");
            let read_hap_support =
                Self::match_reads_and_haps(reads, &haps_to_assess[..], /* min_match= */ None)?;
            log::trace!(
                "Number of uniquely-matching reads: {}. Support: {:?}",
                read_hap_support.support_reads.len(),
                read_hap_support.support_reads
            );
            let filtered_indices = filtered_ass_haps
                .iter()
                .map(|x| x.0 as usize)
                .collect::<Vec<_>>();
            let human_readable_support = read_hap_support
                .support_reads
                .iter()
                .map(|(k, v)| (flat_init_haps[filtered_indices[*k as usize]], &v[..]))
                .collect::<BTreeMap<VStr<'_>, &[ReadHapPair<'_>]>>();
            log::trace!("human readable support at iternum {iternum}: {human_readable_support:?}");
            let local_passing_haps = read_hap_support.support_reads.iter().filter_map(|(k, v)| {
                if v.len() >= min_count {
                    Some(k)
                } else {
                    log::trace!(
                        "Removing read {k}/{} for support which is length {}/{v:?}",
                        v.len(),
                        flat_init_haps[filtered_indices[*k as usize]]
                    );
                    None
                }
            }); // These have indices into the current haps_to_assess vector.
            log::trace!("Remaining passing: {local_passing_haps:?} at iternum {iternum}");

            // Replace local offsets with offsets into good reads by extracting by position.
            let offset_passing_haps = local_passing_haps
                .map(|x| std::mem::take(&mut filtered_ass_haps[*x as usize]))
                .collect::<Vec<_>>();
            filtered_ass_haps = offset_passing_haps;
            log::trace!("Filtered haps after filtering for passing status: {filtered_ass_haps:?}");
            iternum += 1;
            log::trace!(
                "{} haps currently after {iternum} iterations",
                filtered_ass_haps.len()
            );
            if filtered_ass_haps.len() == haps_to_assess.len() {
                log::trace!("Filtered haps {filtered_ass_haps:?} have the same length as haps to assess: {haps_to_assess:?}, breaking loop.");
                break;
            }
        }
        let filtered_ass_haps = AssembledPaths::from_seqs(
            filtered_ass_haps
                .into_iter()
                .map(|x| flat_init_haps[x.0 as usize]),
        );
        log::debug!(
            "{}{filtered_ass_haps:?} haps remain after filtering from {}/{init_haps:?}",
            filtered_ass_haps.len(),
            init_haps.len(),
        );
        Ok(filtered_ass_haps)
    }

    /// Display all nodes from the graph as a vector of strings.
    #[must_use]
    pub fn formatted_nodes(&self) -> Vec<String> {
        self.node_iter()
            .map(|x| {
                self.paraphase_format_node(*self.node_id_map.get(x).unwrap(), None)
                    .unwrap()
            })
            .sorted()
            .collect::<Vec<_>>()
    }

    /// Display an edge as a string.
    #[must_use]
    pub fn format_edge(&self, x: &(u32, u32, u32)) -> Option<(String, String)> {
        let edge_index = self.find_edge(x.0.into(), x.1.into());
        edge_index.and_then(|edge_index| self.paraphase_format_edge(edge_index, None))
    }

    /// Display all nodes from the graph as a vector of string pairs in (from, to) format.
    #[must_use]
    pub fn formatted_edges(&self) -> Vec<Vec<(String, String)>> {
        self.pos_edge_map
            .iter()
            .map(|(_pos, edges)| {
                edges
                    .iter()
                    .map(|edge| {
                        self.format_edge(edge)
                            .unwrap_or_else(|| panic!("missing edge for {edge:?}"))
                    })
                    .sorted()
                    .collect::<Vec<(String, String)>>()
            })
            .collect::<Vec<Vec<(String, String)>>>()
    }
} // Graph

#[cfg(test)]
mod tests {
    use super::*;
    use crate::assembly::variant_graph;
    use crate::detail::range::Range;
    use crate::detail::util::{DError, DResult};
    use crate::io::json::{ReadAlignmentId, ReadFingerprintMap};

    #[test]
    fn test_get_next_pos() -> Result<(), std::boxed::Box<dyn std::error::Error>> {
        use vstr::VString;
        // Set up graph.
        let mut reads = ReadFingerprintMap::new();
        reads.insert(ReadAlignmentId::from_name("r1"), VString::from("11"));
        reads.insert(ReadAlignmentId::from_name("r2"), VString::from("22"));
        let mut graph = Graph::try_new(reads.clone(), None)?;
        let mut add_node = |x: VString, y: u32| {
            let _added = graph.fully_add_node(NodeDatum::new(x, y));
        };
        add_node(VString::from("121"), 0);
        add_node(VString::from("111"), 0);
        add_node(VString::from("222"), 3);
        add_node(VString::from("111"), 3);

        // Get next for 1
        let res = graph.get_next_pos(1).unwrap();
        assert_eq!(res, 3u32);
        let prev_haps = graph.get_hap_by_pos(res).collect::<Vec<_>>();
        assert_eq!(prev_haps[0], "111");
        assert_eq!(prev_haps[1], "222");

        // Get next for 0
        let res = graph.get_next_pos(0).unwrap();
        assert_eq!(res, 3u32);
        let prev_haps = graph.get_hap_by_pos(res).collect::<Vec<_>>();
        assert_eq!(prev_haps[0], "111");
        assert_eq!(prev_haps[1], "222");

        // Get next for 3. Should be 0.
        assert!(graph.get_next_pos(3).is_none());
        Ok(())
    }

    #[test]
    fn test_get_previous_pos() -> Result<(), std::boxed::Box<dyn std::error::Error>> {
        use vstr::VString;
        // Set up graph.
        let mut reads = ReadFingerprintMap::new();
        reads.insert(ReadAlignmentId::from_name("r1"), VString::from("11"));
        reads.insert(ReadAlignmentId::from_name("r2"), VString::from("22"));
        let mut graph = Graph::try_new(reads.clone(), None)?;
        let mut add_node = |x: VString, y: u32| {
            let _added = graph.fully_add_node(NodeDatum::new(x, y));
        };
        add_node(VString::from("121"), 0);
        add_node(VString::from("111"), 0);
        add_node(VString::from("222"), 3);
        add_node(VString::from("111"), 3);

        // Get previous for 1
        let res = graph.get_previous_pos(1).unwrap();
        assert_eq!(res, 0u32);
        let prev_haps = graph.get_hap_by_pos(res).collect::<Vec<_>>();
        assert_eq!(prev_haps[0], "111");
        assert_eq!(prev_haps[1], "121");

        // Get previous for 3
        let res = graph.get_previous_pos(3).unwrap();
        assert_eq!(res, 0u32);
        let prev_haps = graph.get_hap_by_pos(res).collect::<Vec<_>>();
        assert_eq!(prev_haps[0], "111");
        assert_eq!(prev_haps[1], "121");

        // Get previous for 0
        assert!(graph.get_previous_pos(0).is_none());
        Ok(())
    }

    fn make_graph(reads: &ReadFingerprintMap) -> crate::assembly::variant_graph::Graph {
        use petgraph::stable_graph::NodeIndex;

        const NODES: &[(u8, i32)] = &[
            (b'1', 0), // 0
            (b'1', 1), // 1
            (b'1', 2), // 2
            (b'1', 3), // 3
            (b'2', 0), // 4
            (b'2', 1), // 5
            (b'2', 2), // 6
            (b'2', 3), // 7
        ];

        let mut graph = Graph::try_new(reads.clone(), None).expect("Failed to create Graph");
        for node in NODES
            .iter()
            .map(|(chr, pos)| NodeDatum::from((*chr, *pos as usize)))
        {
            // graph.nodes = ["1-0", "1-1", "1-2", "1-3", "2-0", "2-1", "2-2", "2-3"]
            let (_id, was_new) = graph.fully_add_node(node);
            assert!(was_new);
        }
        let str2idx = |x: &&str| -> u32 {
            let (seq, pos) = x.split_terminator('-').next_tuple().unwrap();
            let seq = seq.as_bytes()[0];
            let pos = pos.parse::<i32>().unwrap();
            NODES
                .iter()
                .position(|x| x == &(seq, pos))
                .expect("Need to find node in set") as u32
        };
        let edge_from = ["1-0", "2-0", "1-1", "2-1", "1-2", "2-2"];
        let edge_to = ["1-1", "2-1", "1-2", "2-2", "1-3", "2-3"];
        for (from, to) in edge_from
            .iter()
            .map(str2idx)
            .zip(edge_to.iter().map(str2idx))
        {
            graph.add_edge(NodeIndex::from(from), NodeIndex::from(to), vec![]);
        }
        graph.pos_edge_map = BTreeMap::new();
        //graph.pos_edge_map = vec![vec![]; 3];
        graph.pos_edge_map.insert(
            0,
            vec![
                (str2idx(&"1-0"), str2idx(&"1-1"), 2),
                (str2idx(&"2-0"), str2idx(&"2-1"), 2),
            ],
        );
        graph.pos_edge_map.insert(
            1,
            vec![
                (str2idx(&"1-1"), str2idx(&"1-2"), 2),
                (str2idx(&"2-1"), str2idx(&"2-2"), 2),
            ],
        );
        graph.pos_edge_map.insert(
            2,
            vec![
                (str2idx(&"1-2"), str2idx(&"1-3"), 2),
                (str2idx(&"2-2"), str2idx(&"2-3"), 2),
            ],
        );
        graph
    }

    #[test]
    fn test_rm_add_edges_ok() -> Result<(), std::boxed::Box<dyn std::error::Error>> {
        use std::collections::BTreeSet;
        use vstr::{VStr, VString};
        if let Err(e) = env_logger::builder()
            .format_timestamp_millis()
            .filter_level(log::LevelFilter::Debug)
            .try_init()
        {
            log::trace!("Logger already activated. Error: {e:?}");
        }
        let mut reads = ReadFingerprintMap::new();
        for (read_id, read) in ["x11x", "11xx", "x11x", "x22x", "xx22"].iter().enumerate() {
            reads.insert(
                ReadAlignmentId::from_name(format!("r{}", read_id + 1)),
                VString::from(*read),
            );
        }

        let mut graph = make_graph(&reads);
        let paths = AssembledPaths::from_seqs(["11", "22"]);

        graph.rm_add_edges(1, 2, paths.iter())?;
        let nodes = graph
            .node_iter()
            .map(|x| (x.hap.vstr(), x.pos))
            .collect::<BTreeSet<_>>();
        let expected_nodes = [("1", 0), ("1", 3), ("2", 0), ("2", 3), ("11", 1), ("22", 1)]
            .into_iter()
            .map(|(x, y)| (VStr::from(x), y))
            .collect::<BTreeSet<_>>();
        assert_eq!(expected_nodes, nodes);

        let format_edge = |x: &(u32, u32, u32)| -> Option<(String, String)> {
            let edge_index = graph.find_edge(x.0.into(), x.1.into());
            edge_index.and_then(|edge_index| graph.paraphase_format_edge(edge_index, None))
        };
        let edges = graph
            .pos_edge_map
            .iter()
            .flat_map(|(num, iter)| iter.clone().into_iter().map(move |x| (num, x)))
            .filter_map(|(num, edge)| format_edge(&edge).map(|x| (num, x)))
            .collect::<BTreeSet<_>>();
        let edges = edges.into_iter().map(|x| x.1).collect::<BTreeSet<_>>();
        let expected = [
            ("1-0", "11-1"),
            ("11-1", "1-3"),
            ("2-0", "22-1"),
            ("22-1", "2-3"),
        ]
        .into_iter()
        .map(|(x, y)| (x.to_owned(), y.to_owned()))
        .collect::<BTreeSet<_>>();
        assert_eq!(edges, expected);
        let get_data = |idx1: usize, idx2: usize| -> ((String, u32), (String, u32)) {
            let hap_pos = |x: &NodeDatum| -> (String, u32) { (x.hap.to_string(), x.pos) };
            let (from, to, _weight) = graph.pos_edge_map[&idx1][idx2];
            (
                hap_pos(graph.node_weight(from.into()).expect("No from")),
                hap_pos(graph.node_weight(to.into()).expect("No to")),
            )
        };
        let test = |x: usize, y: usize, s1: &str, s2: &str, z: u32, w: u32| {
            assert_eq!(
                get_data(x, y),
                ((String::from(s1), z), (String::from(s2), w))
            );
        };
        test(0, 0, "1", "11", 0, 1);
        test(0, 1, "2", "22", 0, 1);
        test(1, 0, "11", "1", 1, 3);
        test(1, 1, "22", "2", 1, 3);

        // Now test on edge
        let mut graph = make_graph(&reads);
        graph.rm_add_edges(0, 1, paths.iter())?;

        let nodes = graph
            .node_iter()
            .map(|x| (x.hap.vstr(), x.pos))
            .collect::<BTreeSet<_>>();
        let expected_nodes = [("1", 2), ("1", 3), ("11", 0), ("2", 2), ("2", 3), ("22", 0)]
            .into_iter()
            .map(|(x, y)| (VStr::from(x), y))
            .collect::<BTreeSet<_>>();
        assert_eq!(expected_nodes, nodes);
        let get_data = |idx1: usize, idx2: usize| -> ((String, u32), (String, u32)) {
            let hap_pos = |x: &NodeDatum| -> (String, u32) { (x.hap.to_string(), x.pos) };
            let (from, to, _weight) = graph.pos_edge_map[&idx1][idx2];
            (
                hap_pos(graph.node_weight(from.into()).expect("No from")),
                hap_pos(graph.node_weight(to.into()).expect("No to")),
            )
        };
        let test = |x: usize, y: usize, s1: &str, s2: &str, z: u32, w: u32| {
            assert_eq!(
                get_data(x, y),
                ((String::from(s1), z), (String::from(s2), w))
            );
        };
        let format_edge = |x: &(u32, u32, u32)| -> Option<(String, String)> {
            let edge_index = graph.find_edge(x.0.into(), x.1.into());
            edge_index.and_then(|edge_index| graph.paraphase_format_edge(edge_index, None))
        };

        let edges = graph
            .pos_edge_map
            .iter()
            .flat_map(|(num, iter)| iter.clone().into_iter().map(move |x| (num, x)))
            .filter_map(|(num, edge)| format_edge(&edge).map(|x| (num, x)))
            .collect::<BTreeSet<_>>();

        let expected = [
            ("1-2", "1-3"),
            ("11-0", "1-2"),
            ("2-2", "2-3"),
            ("22-0", "2-2"),
        ]
        .into_iter()
        .map(|(x, y)| (x.to_owned(), y.to_owned()))
        .collect::<BTreeSet<_>>();
        assert_eq!(
            edges.into_iter().map(|x| x.1).collect::<BTreeSet<_>>(),
            expected
        );

        test(0, 0, "11", "1", 0, 2);
        test(0, 1, "22", "2", 0, 2);
        test(2, 0, "1", "1", 2, 3);
        test(2, 1, "2", "2", 2, 3);
        assert!(!graph.pos_edge_map.contains_key(&1));

        Ok(())
    }

    #[test]
    fn get_segments_ok() -> Result<(), std::boxed::Box<dyn std::error::Error>> {
        use crate::assembly::variant_graph::Graph;
        use crate::assembly::variant_graph::SegmentationClass::*;
        use crate::detail::range::Range;
        use crate::io::json::{ReadAlignmentId, ReadFingerprintMap};
        use vstr::VString;
        let mut reads = ReadFingerprintMap::new();
        reads.insert(ReadAlignmentId::from_name("r1"), VString::from("11"));
        reads.insert(ReadAlignmentId::from_name("r2"), VString::from("22"));
        let mut graph = Graph::try_new(reads.clone(), None)?;
        graph.dnhap.insert(0, Two);
        graph.dnhap.insert(1, Ten);
        graph.dnhap.insert(2, Ten);
        graph.dnhap.insert(3, Ten);
        let segments = graph.segment_graph()?;
        assert_eq!(
            *segments
                .get(&Range::<u32>::new(0, 0))
                .expect("Range must be present"),
            Two
        );
        assert_eq!(
            *segments
                .get(&Range::<u32>::new(1, 3))
                .expect("Range must be present"),
            Ten
        );

        let mut graph = Graph::try_new(reads.clone(), None)?;
        graph.dnhap.insert(0, Two);
        graph.dnhap.insert(1, Two);
        graph.dnhap.insert(2, Ten);
        graph.dnhap.insert(3, Ten);
        let segments = graph.segment_graph()?;
        assert_eq!(
            *segments
                .get(&Range::<u32>::new(0, 1))
                .expect("Range must be present"),
            Two
        );
        assert_eq!(
            *segments
                .get(&Range::<u32>::new(2, 3))
                .expect("Range must be present"),
            Ten
        );

        let mut graph = Graph::try_new(reads, None)?;
        graph.dnhap.insert(0, Ten);
        graph.dnhap.insert(1, Two);
        graph.dnhap.insert(2, Two);
        graph.dnhap.insert(3, Ten);
        let segments = graph.segment_graph()?;
        log::trace!("segments: {segments:?}");
        let ranges = segments.keys().cloned().collect::<Vec<_>>();
        let expected_ranges = vec![
            Range::<u32>::new(0, 0),
            Range::<u32>::new(1, 2),
            Range::<u32>::new(3, 3),
        ];
        assert_eq!(ranges, expected_ranges);
        for (ivl, val) in &segments {
            match ivl.start {
                0 => {
                    assert_eq!(ivl.end, 0);
                    assert_eq!(val, &Ten);
                }
                1 => {
                    assert_eq!(ivl.end, 2);
                    assert_eq!(val, &Two);
                }
                3 => {
                    assert_eq!(ivl.end, 3);
                    assert_eq!(val, &Ten);
                }
                _ => {
                    panic!("Not expected");
                }
            };
        }
        Ok(())
    }

    #[test]
    fn path_ok() -> Result<(), std::boxed::Box<dyn std::error::Error>> {
        use super::*;
        use vstr::VString;
        let mut reads = ReadFingerprintMap::new();
        reads.insert(ReadAlignmentId::from_name("r1"), VString::from("11"));
        reads.insert(ReadAlignmentId::from_name("r2"), VString::from("22"));
        let mut graph = Graph::try_new(reads.clone(), None)?;
        graph.init()?;
        let paths = graph.path(Range::<u32>::new(0, 1))?;
        assert_eq!(paths.len(), 2, "paths: {paths:?}");
        assert_eq!(
            paths
                .iter()
                .map(std::string::ToString::to_string)
                .sorted()
                .collect::<Vec<_>>(),
            vec![String::from("11"), String::from("22")]
        );
        Ok(())
    }

    #[test]
    fn match_reads_and_haps_ok() -> Result<(), std::boxed::Box<dyn std::error::Error>> {
        use vstr::VString;
        let haps = ["21xx", "xx12"]
            .into_iter()
            .map(VString::from)
            .collect::<Vec<_>>();
        let haplotype_per_read = [("r1", "x11x"), ("r2", "x11x")]
            .into_iter()
            .map(|(x, y)| (ReadAlignmentId::from_name(x), VString::from(y)))
            .collect::<ReadFingerprintMap>();
        let res = Graph::match_reads_and_haps(
            &haplotype_per_read,
            &haps[..],
            /* min_match= */ None,
        )?;
        let unique = &res.support_reads;
        log::trace!("unique: {unique:?}");
        let match0 = unique.get(&0).unwrap();
        let match1 = unique.get(&1).unwrap();
        assert_eq!(
            match0
                .iter()
                .map(|x| x.path.to_string())
                .collect::<Vec<_>>(),
            &["x11x", "x11x"]
        );
        assert_eq!(
            match1
                .iter()
                .map(|x| x.path.to_string())
                .collect::<Vec<_>>(),
            &["x11x", "x11x"]
        );
        Ok(())
    }

    #[test]
    fn test_get_matching_hap_ok() -> Result<(), std::boxed::Box<dyn std::error::Error>> {
        use crate::assembly::variant_graph::Graph;
        use std::collections::BTreeMap;
        use vstr::VStr;

        let nodes = [(0, ["11", "22"]), (2, ["11", "22"])]
            .into_iter()
            .map(|(key, vec)| (key, (vec.into_iter().map(VStr::from).collect::<Vec<_>>())))
            .collect::<BTreeMap<_, _>>();
        let (hap, matches) = Graph::get_matching_hap(&nodes, 0, VStr::from("x11x"))?;
        assert_eq!(hap, "x1");
        let (score, hits) = matches.expect("No matches");
        assert_eq!(score, 1);
        assert_eq!(hits, vec![VStr::from("11")]);

        let (hap, matches) = Graph::get_matching_hap(&nodes, 2, VStr::from("x11x"))?;
        assert_eq!(hap, "1x");
        let (score, hits) = matches.expect("No matches");
        assert_eq!(score, 1);
        assert_eq!(hits, vec![VStr::from("11")]);
        Ok(())
    }

    #[test]
    fn test_merge_two_pos_ok() -> Result<(), std::boxed::Box<dyn std::error::Error>> {
        use crate::assembly::variant_graph::Graph;
        use crate::io::json::{ReadAlignmentId, ReadFingerprintMap};
        use vstr::{VStr, VString};

        if let Err(e) = env_logger::builder()
            .format_timestamp_millis()
            .filter_level(log::LevelFilter::Debug)
            .try_init()
        {
            log::trace!("Logger already activated. Error: {e:?}");
        }

        let data = ["x11x", "11xx", "x11x", "x22x", "xx22"]
            .iter()
            .enumerate()
            .map(|(idx, hap)| {
                (
                    ReadAlignmentId::from_name(format!("r{}", idx + 1)),
                    VString::from(*hap),
                )
            })
            .collect::<ReadFingerprintMap>();
        let graph = Graph::try_new(data, None)?;
        let nodes = [(0, ["11", "22"]), (2, ["11", "22"])]
            .into_iter()
            .map(|(key, v)| {
                (
                    key as u32,
                    v.into_iter()
                        .map(VStr::from)
                        .collect::<std::collections::BTreeSet<_>>(),
                )
            })
            .collect::<std::collections::BTreeMap<_, _>>();
        let (sub_hap, dnext, dbefore) = graph.merge_two_pos(&nodes)?;
        assert_eq!(
            dnext[&VStr::from("11")].iter().next(),
            Some(&VStr::from("11"))
        );
        assert_eq!(
            dnext[&VStr::from("22")].iter().next(),
            Some(&VStr::from("22"))
        );
        assert_eq!(
            dbefore[&VStr::from("11")].iter().next(),
            Some(&VStr::from("11"))
        );
        assert_eq!(
            dbefore[&VStr::from("22")].iter().next(),
            Some(&VStr::from("22"))
        );
        assert_eq!(
            sub_hap[&"1111".into()].clone(),
            vec![
                &ReadAlignmentId::from_name("r1"),
                &ReadAlignmentId::from_name("r3")
            ]
        );
        assert_eq!(
            sub_hap[&"2222".into()].clone(),
            vec![&ReadAlignmentId::from_name("r4")]
        );
        assert_eq!(sub_hap.len(), 2);
        Ok(())
    }

    #[test]
    fn rescue_missing_ok() {
        use crate::assembly::variant_graph::Graph;
        use std::collections::*;
        use vstr::VString;
        let left_seq = VString::from("11");
        let right_seq = VString::from("22");
        let left = vec![&left_seq];
        let right = vec![&right_seq];
        let dnext: BTreeMap<VString, BTreeSet<VString>> = [("11", "22")]
            .into_iter()
            .map(|(x, y)| {
                (
                    VString::from(x),
                    [y].into_iter().map(VString::from).collect::<BTreeSet<_>>(),
                )
            })
            .collect::<_>();
        let dbefore: BTreeMap<VString, BTreeSet<VString>> = [("22", "11")]
            .into_iter()
            .map(|(x, y)| {
                (
                    VString::from(x),
                    [y].into_iter().map(VString::from).collect::<BTreeSet<_>>(),
                )
            })
            .collect::<_>();

        let rescued = Graph::rescue_missing(&left, &right, &dnext, &dbefore);
        assert_eq!(
            rescued,
            [VString::from("1122")].into_iter().collect::<BTreeSet<_>>()
        );
    }

    #[test]
    fn get_missing_ok() {
        use vstr::VString;
        // Set up graph.
        let haps_left = ["11", "22"]
            .into_iter()
            .map(VString::from)
            .collect::<Vec<_>>();
        let haps_right = ["12", "21"]
            .into_iter()
            .map(VString::from)
            .collect::<Vec<_>>();
        let test_haps = [VString::from("1112")]
            .into_iter()
            .collect::<std::collections::BTreeSet<_>>();
        let (left, right) = Graph::get_missing(&test_haps, &haps_left[..], &haps_right);
        assert_eq!(left, vec![&VString::from("22")]);
        assert_eq!(right, vec![&VString::from("21")]);
    }

    #[test]
    fn get_highest_cn_ok() -> Result<(), DError> {
        use crate::assembly::variant_graph::Graph;
        use crate::io::json::ReadAlignmentId;
        use vstr::VString;
        let reads = [("r1", "11111"), ("r2", "22222")]
            .into_iter()
            .map(|(name, fprint)| (ReadAlignmentId::from_name(name), VString::from(fprint)))
            .collect::<std::collections::BTreeMap<_, _>>();
        let graph = Graph::try_new(reads, None)?;
        let (highest_cn, _) = graph.get_highest_cn(["11111", "22222"].into_iter())?;
        assert_eq!(highest_cn, 2);
        let (highest_cn, _) = graph.get_highest_cn(["121xxx", "111111", "xxx111"].into_iter())?;
        assert_eq!(highest_cn, 2);
        let (highest_cn, _) = graph.get_highest_cn(["1211xx", "111111", "xxx111"].into_iter())?;
        assert_eq!(highest_cn, 3);
        Ok(())
    }

    #[test]
    fn order_blocks_by_size_ok() -> Result<(), DError> {
        use crate::assembly::node_datum::NodeDatum;
        use crate::assembly::variant_graph::Graph;
        use crate::io::json::ReadAlignmentId;
        use vstr::VString;
        let reads = ["111x", "211x", "x111", "x112", "xx12", "1111"]
            .into_iter()
            .enumerate()
            .map(|(idx, item)| {
                (
                    ReadAlignmentId::from_name(format!("r{}", idx + 1)),
                    VString::from(item),
                )
            })
            .collect::<std::collections::BTreeMap<_, _>>();
        let mut graph = Graph::try_new(reads, None)?;
        for node in ["11-0", "22-0", "11-2", "22-2", "111-4", "222-4"]
            .into_iter()
            .map(|x| NodeDatum::try_from(x).unwrap())
        {
            let _added = graph.fully_add_node(node);
        }
        let block_order = graph.order_blocks_by_size();
        assert_eq!(block_order, &[[2, 4, 2, 3u32], [0, 2, 2, 2u32]]);
        Ok(())
    }

    #[test]
    fn extend_pivot_blocks_ok() -> Result<(), DError> {
        use crate::assembly::node_datum::NodeDatum;
        use crate::assembly::variant_graph::Graph;
        use crate::assembly::variant_graph::SegmentationClass::*;
        use crate::io::json::ReadAlignmentId;
        let reads = [
            "x2111x", "x11xxx", "xxx22x", "x2111x", "xxx122", "1211xx", "xx1122", "xx11xx",
            "x111x1", "x1x12x",
        ]
        .into_iter()
        .enumerate()
        .map(|(idx, hap)| {
            (
                ReadAlignmentId::from_name(format!("r{}", (idx % 10))),
                vstr::VString::from(hap),
            )
        })
        .collect::<std::collections::BTreeMap<_, _>>();
        let mut graph = Graph::try_new(reads, None)?;
        for node in ["121-0", "111-0", "122-3", "111-3"]
            .into_iter()
            .map(|x| NodeDatum::try_from(x).unwrap())
        {
            let _added = graph.fully_add_node(node);
        }
        graph.add_edge(0.into(), 0.into(), vec![]);
        graph.dnhap.insert(2, Ten);
        //graph.pos_edge_map = vec![vec![], vec![], vec![(0, 0, 2), (0, 0, 3)], vec![]];
        graph.pos_edge_map = BTreeMap::new();
        graph.pos_edge_map.insert(0, vec![]);
        graph.pos_edge_map.insert(1, vec![]);
        graph.pos_edge_map.insert(2, vec![(0, 0, 2), (0, 0, 3)]);
        graph.pos_edge_map.insert(3, vec![]);
        graph.edge_info = vec![vec![], vec![], vec![2, 3], vec![]];
        let extended_haps = graph.extend_pivot_blocks(["111xxx"]);
        assert_eq!(extended_haps.len(), 1);
        assert_eq!(
            extended_haps
                .iter()
                .next()
                .map(std::string::ToString::to_string)
                .unwrap(),
            "1111xx"
        );
        Ok(())
    }

    #[test]
    fn filter_low_support_haps_ok() -> DResult {
        let reads = [
            "11x", "22x", "x12", "x22", "x2x", "x22", "x11", "xx1", "xx1",
        ]
        .into_iter()
        .enumerate()
        .map(|(idx, hap)| {
            (
                ReadAlignmentId::from_name(format!("r{}", idx + 1)),
                vstr::VString::from(hap),
            )
        })
        .collect::<std::collections::BTreeMap<_, _>>();
        let mut graph = Graph::try_new(reads, None)?;
        graph.settings.min_hap_support = 3;
        let filtered_haps = graph.filter_low_support_haps(&AssembledPaths::from_seqs(
            ["112", "111", "222"].into_iter(),
        ))?;
        assert_eq!(
            filtered_haps.iter().map(|x| x.vstr()).collect::<Vec<_>>(),
            vec![vstr::VStr::from("111"), vstr::VStr::from("222")]
        );
        graph.settings.min_hap_support = 1;
        let filtered_haps = graph.filter_low_support_haps(&AssembledPaths::from_seqs(
            ["112", "111", "222"].into_iter(),
        ))?;
        assert_eq!(
            filtered_haps.iter().map(|x| x.vstr()).collect::<Vec<_>>(),
            vec![
                vstr::VStr::from("111"),
                vstr::VStr::from("112"),
                vstr::VStr::from("222")
            ]
        );
        Ok(())
    }

    #[test]
    fn run_ok() -> DResult {
        use super::*;
        use vstr::VString;
        let reads = ["11", "22", "11", "22", "11", "22"]
            .into_iter()
            .enumerate()
            .map(|(idx, hap)| {
                (
                    ReadAlignmentId::from_name(format!("r{}", idx + 1)),
                    VString::from(hap),
                )
            })
            .collect::<std::collections::BTreeMap<_, _>>();
        let settings = variant_graph::Settings::from_pivot(Some(1));
        let mut graph = Graph::try_new(reads, Some(settings))?;
        graph.init()?;
        let res = graph.assemble_haps()?;
        assert_eq!(
            res.into_inner()
                .into_iter()
                .map(|x| x.to_string())
                .sorted()
                .collect::<Vec<_>>(),
            ["11", "22"]
                .into_iter()
                .map(String::from)
                .collect::<Vec<_>>()
        );
        Ok(())
    }

    #[test]
    fn filter_low_complex_ok() -> DResult {
        use itertools::Itertools;
        let haps = [
            "121xxxxxxxxxx",
            "212xxxxxxxxxx",
            "xxx1xxxxxxxxx",
            "xxx2122222222",
            "xxxx222222222",
            "xxxxx11111111",
        ]
        .into_iter()
        .map(vstr::VStr::from)
        .sorted()
        .collect::<Vec<_>>();
        let reads = [
            ("m54086U_220112_213432/103743586/ccs", "xxx1xxxxxxxxx"),
            ("m54086U_220112_213432/108529606/ccs", "xxxx1222222xx"),
            ("m54086U_220112_213432/110036868/ccs", "xxx21xxxxxxxx"),
            ("m54086U_220112_213432/11665889/ccs", "x12xxxxxxxxxx"),
            ("m54086U_220112_213432/119212720/ccs", "xxx21xxxxxxxx"),
            ("m54086U_220112_213432/123667055/ccs", "xxxxxx2222x22"),
            ("m54086U_220112_213432/129172024/ccs", "xxx12xxxxxxxx"),
            ("m54086U_220112_213432/132579760/ccs", "xxxxx22222222"),
            ("m54086U_220112_213432/133891184/ccs", "xxxxxxx222222"),
            ("m54086U_220112_213432/135923990/ccs", "xxxxx11111111"),
            ("m54086U_220112_213432/136381648/ccs", "x21xxxxxxxxxx"),
            ("m54086U_220112_213432/137366144/ccs", "xxx2xxxxxxxxx"),
            ("m54086U_220112_213432/146474465/ccs", "xxxx2xxxxxxxx"),
            ("m54086U_220112_213432/147720103/ccs", "212xxxxxxxxxx"),
            ("m54086U_220112_213432/150339693/ccs", "x21xxxxxxxxxx"),
            ("m54086U_220112_213432/156894020/ccs", "xxxxxx1111111"),
            ("m54086U_220112_213432/157878867/ccs", "121xxxxxxxxxx"),
            ("m54086U_220112_213432/159712931/ccs", "xxx2xxxxxxxxx"),
            ("m54086U_220112_213432/165676282/ccs", "121xxxxxxxxxx"),
            ("m54086U_220112_213432/16909085/ccs", "xxx1xxxxxxxxx"),
            ("m54086U_220112_213432/171443538/ccs", "xxx1xxxxxxxxx"),
            ("m54086U_220112_213432/178193776/ccs", "xxxxx22222222"),
            ("m54086U_220112_213432/179963762/ccs", "121xxxxxxxxxx"),
            ("m54086U_220112_213432/26346554/ccs", "xxx2xxxxxxxxx"),
            ("m54086U_220112_213432/29624360/ccs", "xxx1xxxxxxxxx"),
            ("m54086U_220112_213432/30540767/ccs", "xxxx122222222"),
            ("m54086U_220112_213432/43581827/ccs", "xxxx1xxxxxxxx"),
            ("m54086U_220112_213432/525604/ccs", "xxxxx11111111"),
            ("m54086U_220112_213432/56231844/ccs", "212xxxxxxxxxx"),
            ("m54086U_220112_213432/57672451/ccs", "121xxxxxxxxxx"),
            ("m54086U_220112_213432/60623188/ccs", "212xxxxxxxxxx"),
            ("m54086U_220112_213432/61015098/ccs", "xxxxxx1111111"),
            ("m54086U_220112_213432/62915752/ccs", "xxxxxxxxx1111"),
            ("m54086U_220112_213432/66456232/ccs", "212xxxxxxxxxx"),
            ("m54086U_220112_213432/67699442/ccs", "xxxx22xxxxxxx"),
            ("m54086U_220112_213432/71567283/ccs", "1xxxxxxxxxxxx"),
            ("m54086U_220112_213432/73205791/ccs", "2xxxxxxxxxxxx"),
            ("m54086U_220112_213432/74252490/ccs", "xxx1xxxxxxxxx"),
            ("m54086U_220112_213432/78777158/ccs", "xxxxx222xx2x2"),
            ("m54086U_220112_213432/84085293/ccs", "21xxxxxxxxxxx"),
            ("m54086U_220112_213432/86442681/ccs", "xxx1xxxxxxxxx"),
            ("m54086U_220112_213432/91161760/ccs", "xxxxx11111111"),
            ("m54086U_220112_213432/93062608/ccs", "x21xxxxxxxxxx"),
            ("m54086U_220112_213432/93128452/ccs", "121xxxxxxxxxx"),
            ("m54086U_220112_213432/9962031/ccs", "xxx1xxxxxxxxx"),
        ]
        .into_iter()
        .map(|(k, v)| (ReadAlignmentId::from_name(k), vstr::VString::from(v)))
        .collect::<ReadFingerprintMap>();
        let res = Graph::filter_low_support_haps_detail(
            &AssembledPaths::from_seqs(haps.iter()),
            4,
            &reads,
        )?;
        let expected = AssembledPaths::from_seqs([
            "212xxxxxxxxxx",
            "121xxxxxxxxxx",
            "xxx2122222222",
            "xxx1xxxxxxxxx",
            "xxxxx11111111",
        ]);
        assert_eq!(expected, res);
        Ok(())
    }
}
