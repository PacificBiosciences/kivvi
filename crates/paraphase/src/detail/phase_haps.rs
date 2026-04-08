use crate::assembly::assembly_result::{AssembledPaths, AssemblyResult, AssemblyResultForJson};
use crate::assembly::variant_graph::{self, Graph as VariantGraph, ReadHapPair, ReadHapSupport};
use crate::detail::phaser_util::{base_qual, defaults, get_start_end};
use crate::detail::range;
use crate::detail::site_selection::{self, CandidateSite};
use crate::detail::util::{raw_qpos, DError, HashMap};
use crate::io::json::{ReadAlignmentId, ReadFingerprintMap};
use crate::phaser;

use vstr::{VStr, VString};

use itertools::{enumerate, iproduct, Itertools};
use rust_htslib::bam::{self, pileup::Indel, Read};

use std::cmp::Ordering;
use std::collections::{BTreeMap, BTreeSet};

type HapToReads = BTreeMap<VString, Vec<String>>;
type ReadToPossibleHaps = BTreeMap<String, Vec<VString>>;

/// Tuple containing assemblies, read assignments, raw read haps, and counts for each haplotype.
#[derive(Debug, Default, Clone)]
pub struct PhasedResult {
    pub assemblies: AssemblyResult, // final haps, main haps, hcn
    pub uniquely_supporting_reads: HapToReads,
    pub nonuniquely_supporting_reads: ReadToPossibleHaps,
    pub raw_read_haps: ReadFingerprintMap,
    pub read_counts: BTreeMap<VString, i32>,
}

/// `PhasedResult`, but converted to easily-displayed types for JSON output.
#[derive(Debug, Default, Clone, serde::Deserialize, serde::Serialize)]
pub struct PhasedResultForJson {
    pub assemblies: AssemblyResultForJson, // final haps, main haps, hcn
    pub uniquely_supporting_reads: BTreeMap<String, Vec<String>>,
    pub nonuniquely_supporting_reads: BTreeMap<String, Vec<String>>,
    pub raw_read_haps: BTreeMap<String, String>,
    pub read_counts: BTreeMap<String, i32>,
}

impl PhasedResultForJson {
    /// Generate a struct which can be printed via `serde` from a `PhasedResult` object.
    #[must_use]
    pub fn new(x: &PhasedResult) -> Self {
        let assemblies = AssemblyResultForJson::new(&x.assemblies);
        let uniquely_supporting_reads = x
            .uniquely_supporting_reads
            .iter()
            .map(|(hap, reads)| (hap.to_string(), reads.clone()))
            .collect::<BTreeMap<String, _>>();
        let nonuniquely_supporting_reads = x
            .nonuniquely_supporting_reads
            .iter()
            .map(|(read, haps)| {
                (
                    read.clone(),
                    haps.iter()
                        .map(std::string::ToString::to_string)
                        .collect::<Vec<_>>(),
                )
            })
            .collect::<BTreeMap<String, _>>();
        let raw_read_haps = x
            .raw_read_haps
            .iter()
            .map(|(k, v)| (k.to_string(), v.to_string()))
            .collect::<BTreeMap<String, _>>();
        let read_counts = x
            .read_counts
            .iter()
            .map(|x| (x.0.to_string(), *x.1))
            .collect::<BTreeMap<String, i32>>();
        Self {
            assemblies,
            uniquely_supporting_reads,
            nonuniquely_supporting_reads,
            raw_read_haps,
            read_counts,
        }
    }
}

pub type ReadSupport<'a> = (
    BTreeMap<VStr<'a>, Vec<VStr<'a>>>,
    BTreeMap<VStr<'a>, &'a [u32]>,
    BTreeMap<VStr<'a>, i32>,
);

/// Assignment for token - reference, alt, or other (missing or special)>
#[derive(Debug, PartialEq, Copy, Clone, Eq, PartialOrd, Ord, Hash)]
pub enum Assignment {
    Ref,
    Alt,
    Dot,
}

#[derive(Debug, PartialEq, Copy, Clone, Eq)]
pub enum Genotype {
    One,
    Zero,
    Dot,
}

/// Haplotype summary information.
/// Set of variants, reference coordinate range, and whether or not it is truncated.
/// `boundary_gene2` marks gene2 locations.
#[derive(Clone, Debug)]
pub struct HapInfo {
    pub variants: Vec<CandidateSite>,
    pub boundary: range::I64,
    pub boundary_gene2: Option<range::I64>,
    pub is_truncated: Vec<String>,
}

#[derive(Clone, Debug, serde::Deserialize, serde::Serialize)]
pub struct HapInfoForJson {
    pub variants: Vec<String>,
    pub boundary: String,
    pub boundary_gene2: Option<String>,
    pub is_truncated: Vec<String>,
}

impl std::convert::From<&HapInfo> for HapInfoForJson {
    fn from(info: &HapInfo) -> Self {
        let variants = info
            .variants
            .iter()
            .map(std::string::ToString::to_string)
            .collect::<Vec<_>>();
        let bound_to_string =
            |x: &range::I64| -> String { format!("{}-{}", x.start + 1, x.end + 1) };
        let boundary = bound_to_string(&info.boundary);
        let boundary_gene2 = info.boundary_gene2.as_ref().map(bound_to_string);
        let is_truncated = info.is_truncated.clone();
        Self {
            variants,
            boundary,
            boundary_gene2,
            is_truncated,
        }
    }
}

impl phaser::Phaser {
    ///
    /// Generate hap: read map from read->hap map.
    ///
    ///```
    /// use paraphase::io::json::{ReadAlignmentId, ReadFingerprintMap};
    /// use paraphase::phaser::Phaser;
    /// use vstr::{VStr, VString};
    ///
    /// let haps = [("r1", "12"), ("r2", "21")]
    ///     .into_iter()
    ///     .map(|(k, v)| (ReadAlignmentId::from_name(k), VString::from(v)))
    ///     .collect::<ReadFingerprintMap>();
    /// let (haplotypes_to_reads, reads_to_haplotypes) = Phaser::simplify_read_haps(&haps);
    /// assert_eq!(
    ///     haplotypes_to_reads[&VStr::from("12")],
    ///     vec![&ReadAlignmentId::from_name("r1")]
    /// );
    /// assert_eq!(
    ///     haplotypes_to_reads[&VStr::from("21")],
    ///     vec![&ReadAlignmentId::from_name("r2")]
    /// );
    /// assert_eq!(haplotypes_to_reads.len(), 2);
    /// assert_eq!(
    ///     reads_to_haplotypes[&ReadAlignmentId::from_name("r1")],
    ///     VString::from("12")
    /// );
    /// assert_eq!(
    ///     reads_to_haplotypes[&ReadAlignmentId::from_name("r2")],
    ///     VString::from("21")
    /// );
    /// assert_eq!(reads_to_haplotypes.len(), 2);
    /// ```
    #[must_use]
    pub fn simplify_read_haps(
        reads: &ReadFingerprintMap,
    ) -> (
        BTreeMap<VStr<'_>, Vec<&ReadAlignmentId>>,
        BTreeMap<&ReadAlignmentId, VStr<'_>>,
    ) {
        let mut haps_to_reads = BTreeMap::<VStr<'_>, Vec<&ReadAlignmentId>>::new();
        let mut reads_to_haps = BTreeMap::<&ReadAlignmentId, VStr<'_>>::new();
        for (read, hap) in reads {
            haps_to_reads.entry(hap.vstr()).or_default().push(read);
            reads_to_haps.insert(read, hap.vstr());
        }
        (haps_to_reads, reads_to_haps)
    }

    /// Identify problematic haplotypes and eliminate them.
    fn adjust_spurious_haplotypes<'a>(
        &self,
        current_hap_asm: &'a BTreeMap<VStr<'a>, Vec<VStr<'a>>>,
        flanking_bp: Option<i32>,
        min_base_quality: Option<u8>,
    ) -> Result<AssembledPaths, DError> {
        let min_base_quality = min_base_quality.unwrap_or(defaults::MIN_BASE_QUALITY);
        let flanking_bp = i64::from(flanking_bp.unwrap_or(defaults::FLANKING_SPURIOUS_BP));
        let mut passing = AssembledPaths::from_seqs(current_hap_asm.keys());
        let mut suspicion = Vec::new();
        // TODO: consider a single iteration for the N-choose-2 instead of all-pairs.
        for (hap1, hap2) in
            iproduct!(current_hap_asm.keys(), current_hap_asm.keys()).filter(|x| x.0 != x.1)
        {
            let (matches, mismatches, sites) = hap1
                .iter()
                .copied()
                .zip(hap2.iter().copied())
                .enumerate()
                .filter(|(_id, (x, y))| *x != b'x' && *y != b'x')
                .fold(
                    (0i32, 0i32, vec![]),
                    |(mut matches, mut mismatches, mut sites), (idx, (x, y))| {
                        let is_1_or_2 = |x: u8| -> bool { matches!(x, b'1' | b'2') };
                        if x == y {
                            matches += 1;
                        } else if is_1_or_2(x) && is_1_or_2(y) {
                            mismatches += 1;
                            sites.push(idx as i32);
                        }
                        (matches, mismatches, sites)
                    },
                );
            log::trace!("[adjust_spurious_haplotypes] hap1 {:?} vs. hap2 {:?}, matches {matches} mismatches {mismatches} sites {sites:?}", hap1, hap2);
            if matches >= 5 && mismatches == 1 && sites.len() == 1 {
                let mismatch_site = &self.het_sites[sites[0] as usize];
                let mismatch_pos = mismatch_site.pos;
                let hap1_reads = current_hap_asm.get(hap1).unwrap();
                let hap2_reads = current_hap_asm.get(hap2).unwrap();
                if let Some(pair) = if hap1_reads.len() <= 5 && hap2_reads.len() >= 6 {
                    Some((hap2, hap1))
                } else if hap2_reads.len() <= 5 && hap1_reads.len() >= 6 {
                    Some((hap1, hap2))
                } else {
                    None
                } {
                    suspicion.push((pair, mismatch_pos));
                }
            }
        }
        log::trace!("[adjust_spurious_haplotypes] suspicion {:?}", suspicion);

        // Now check the reads for these troublesome sites.
        suspicion.sort();
        suspicion.dedup();
        for ((hap1, hap2), mismatch_pos) in &suspicion {
            log::trace!(
                "[adjust_spurious_haplotypes] checking reads: hap1 {:?} vs. hap2 {:?}",
                hap1,
                hap2
            );
            let hap1_reads = current_hap_asm.get(*hap1).unwrap();
            let hap2_reads = current_hap_asm.get(*hap2).unwrap();
            let mut hap1_at_pos = BTreeSet::<VString>::new();
            let mut hap2_at_pos = BTreeSet::<VString>::new();
            let mut bam = self.realigned_bam();
            let tid = self.genome_tid().map(|x| x as i32).expect("No chr tid");
            bam.fetch((tid, *mismatch_pos - 1, *mismatch_pos))?; // off by one from Paraphase because already stored as 0-based.
            for x in bam.pileup() {
                let x = x.expect("Error in pileup");
                match i64::from(x.pos()).cmp(mismatch_pos) {
                    Ordering::Less => {
                        continue;
                    }
                    Ordering::Greater => {
                        break;
                    }
                    Ordering::Equal => {}
                };
                for aln in x.alignments().filter(|x| !x.is_del() && !x.is_refskip()) {
                    let record = aln.record();
                    if base_qual(&aln) < min_base_quality {
                        continue;
                    }
                    let qpos = raw_qpos(&aln) as i64;
                    let read_name = self.get_read_name(&record);
                    let has_hap1 = hap1_reads.contains(&VStr::from(&read_name[..]));
                    let has_hap2 = hap2_reads.contains(&VStr::from(&read_name[..]));
                    log::trace!("[adjust_spurious_haplotypes] read_name {read_name} has_hap1 {has_hap1} has_hap2 {has_hap2} qpos {qpos}");
                    if (has_hap1 || has_hap2)
                        && qpos >= flanking_bp
                        && (qpos + flanking_bp < record.seq_len() as i64)
                    {
                        let seq = record.seq().as_bytes();
                        let start = qpos - flanking_bp;
                        let end = qpos + flanking_bp;
                        let slice = &seq[start as usize..end as usize];
                        if has_hap1 {
                            hap1_at_pos.insert(slice.into());
                        }
                        if has_hap2 {
                            hap2_at_pos.insert(slice.into());
                        }
                    }
                }
            }
            log::trace!("[adjust_spurious_haplotypes] hap1 {:?} vs. hap2 {:?} reads {hap1_at_pos:?} {hap2_at_pos:?}", hap1, hap2);
            if hap1_at_pos.intersection(&hap2_at_pos).next().is_some() {
                passing.remove(&VString::from(*hap2));
            }
        } // end for each suspicious range
        Ok(passing)
    }

    /// Take read fingerprint maps and generate phased haplotypes from them.
    pub(crate) fn phase_haps(&self, reads: &ReadFingerprintMap) -> Result<PhasedResult, DError> {
        let mut min_hap_support = self.settings.min_hap_support as f32;
        if self.settings.targeted {
            let total_depth = self.region_avg_depth[0].0;
            min_hap_support =
                min_hap_support.max(total_depth * self.settings.min_haplotype_frequency as f32);
        }
        log::debug!("min_hap_support is {min_hap_support}");

        let het_sites = self.het_sites.clone();
        let (haps_to_reads, _raw_read_haps) = Self::simplify_read_haps(reads);
        assert!(_raw_read_haps
            .iter()
            .all(|(k, v)| reads.get(k) == Some(&VString::from(v))));
        let mut ret = PhasedResult {
            raw_read_haps: reads.clone(),
            ..Default::default()
        };

        let nvar = het_sites.len();
        if nvar == 0 {
            return Ok(ret);
        }
        ret.assemblies = if nvar == 1 {
            let final_haps = AssembledPaths::from_seqs(["1", "2"].into_iter().map(VString::from));
            let main_haps = final_haps.clone();
            AssemblyResult {
                final_haps,
                main_haps,
                highest_cn: 2,
            }
        } else {
            let pivot_index = self.get_pivot_index();
            let mut settings = variant_graph::Settings::from_pivot(pivot_index);
            settings.min_hap_support = min_hap_support.ceil() as i32;
            let mut graph = variant_graph::Graph::new(reads.clone(), Some(settings));
            graph.construct()?
        };
        log::debug!("Assembled before filtering: {ret:?}");
        if ret.assemblies.main_haps.is_empty() {
            return Ok(ret);
        }

        // And resolve current haps.
        // First, assign:
        let mut flat_haps = ret.assemblies.main_haps.iter().cloned().collect::<Vec<_>>();
        let mut support = VariantGraph::match_reads_and_haps(reads, &flat_haps[..], None)?;
        let mut read_support = Self::get_read_support(&haps_to_reads, &flat_haps[..], &support)?;
        let mut assembled_haps = self.adjust_spurious_haplotypes(&read_support.0, None, None)?;
        log::debug!(
            "[phase_haps] Num asm after adjust: {}. Before: {}. After: {assembled_haps:?}. Before: {flat_haps:?}",
            assembled_haps.len(),
            flat_haps.len()
        );

        flat_haps = assembled_haps.iter().cloned().collect::<Vec<_>>();
        support = VariantGraph::match_reads_and_haps(reads, &flat_haps[..], None)?;
        read_support = Self::get_read_support(&haps_to_reads, &flat_haps[..], &support)?;
        let uniquely_supporting_reads = &read_support.0;
        log::trace!("uniquely_supporting_reads: {uniquely_supporting_reads:?}");
        if true {
            // Filter low-support haps.
            let mut flat_read_counts = uniquely_supporting_reads
                .values()
                .map(std::vec::Vec::len)
                .sorted();
            log::trace!("flat_read_counts: {flat_read_counts:?}");
            let min_hap_support = if min_hap_support == 4.0
                && flat_read_counts.len() > 2
                && flat_read_counts.next().map_or(false, |x| x <= 4)
                && flat_read_counts.next().map_or(false, |x| x >= 12)
                && !assembled_haps.iter().any(|hap| hap.contains(&b'x'))
            {
                5.0
            } else {
                min_hap_support
            };
            log::debug!("min_hap_support {min_hap_support}");
            // Now re-assign
            assembled_haps =
            AssembledPaths::from_seqs(uniquely_supporting_reads.iter().filter_map(|x| {
                if x.1.len() as f32 >= min_hap_support {
                    Some(x.0)
                } else {
                    log::debug!("Eliminating hap with support {} as failing compared to {min_hap_support}. hap: {}", x.1.len(), x.0);
                    None
                }
            }));
            log::debug!(
            "[phase_haps]Num asm after low support filtering: {}. Before: {}. After: {assembled_haps:?}. Before: {flat_haps:?}",
            assembled_haps.len(),
            flat_haps.len()
        );
        }
        flat_haps = assembled_haps.iter().cloned().collect::<Vec<_>>();
        let support = VariantGraph::match_reads_and_haps(reads, &flat_haps[..], None)?;
        let read_support = Self::get_read_support(&haps_to_reads, &flat_haps[..], &support)?;
        let (ref uniquely_supporting_reads, ref nonuniquely_supporting_reads, ref read_counts) =
            read_support;

        // Assign outputs.
        ret.assemblies.main_haps = assembled_haps;
        ret.uniquely_supporting_reads = uniquely_supporting_reads
            .iter()
            .map(|(k, v)| {
                (
                    VString::from(k),
                    v.iter()
                        .map(std::string::ToString::to_string)
                        .collect::<Vec<_>>(),
                )
            })
            .collect::<BTreeMap<VString, Vec<String>>>();
        ret.nonuniquely_supporting_reads = nonuniquely_supporting_reads
            .iter()
            .map(|(name, ids)| {
                (
                    name.to_string(),
                    ids.iter()
                        .map(|id| VString::from(&flat_haps[*id as usize]))
                        .collect::<Vec<_>>(),
                )
            })
            .collect::<BTreeMap<String, Vec<VString>>>();
        ret.read_counts = read_counts
            .iter()
            .map(|(k, v)| (VString::from(k), *v))
            .collect::<BTreeMap<VString, i32>>();
        Ok(ret)
    }

    /// In Paraphase, `None` is returned if there is a failure.
    /// Instead, we return an empty dictionary.
    fn get_read_counts<'a>(
        hap2read: &'a HashMap<u32, Vec<ReadHapPair<'a>>>,
        paths: &'a [VString],
    ) -> BTreeMap<VStr<'a>, i32> {
        log::trace!(
            "[get_read_counts] hap2read of size {}/{} is {hap2read:?} for {paths:?}",
            hap2read.len(),
            paths.len()
        );
        let mut ret = BTreeMap::<VStr<'_>, i32>::new();
        if hap2read.is_empty() {
            return ret;
        }
        let nvar: VStr<'_> = paths.first().unwrap().vstr();
        let nvar = nvar.len();
        let nhap = hap2read.len();
        log::trace!("[get_read_counts] nvar {nvar}, nhap {nhap}");
        let mut hap_base_counts = vec![Vec::<i32>::new(); nhap];
        // assert_eq!(hap_base_counts.len(), paths.len());
        for (path_idx, hap_counts) in hap_base_counts.iter_mut().enumerate() {
            let matches = hap2read.get(&(path_idx as u32)).unwrap_or_else(|| {
                panic!("Missing path idx {path_idx} from hap2read {hap2read:?}");
            });
            // TODO: consider transposing loop and putting outside of path loop.
            for i in 0..nvar {
                let counted_bases = matches
                    .iter()
                    .filter(|x| !matches!(x.path[i], b'x' | b'0'))
                    .count();
                hap_counts.push(counted_bases as i32);
            }
        }
        // Identify low-support ranges.
        let mut ranges = Vec::new();
        assert!(hap_base_counts.iter().all(|x| x.len() == nvar));
        for fingerprint_idx in 0..nvar {
            if hap_base_counts
                .iter()
                .map(|x| x[fingerprint_idx])
                .min()
                .map_or(false, |x| x >= 5)
            {
                //  Loop i + 1 -> nvar instead of 0..nvar and skip loops where j_idx <= fingerprint_idx
                for j_idx in (fingerprint_idx + 1)..nvar {
                    if j_idx == nvar - 1
                        || hap_base_counts
                            .iter()
                            .map(|x| x[j_idx])
                            .min()
                            .map_or(true, |x| x < 5)
                    {
                        ranges.push(range::I64::new(fingerprint_idx as i64, j_idx as i64));
                        break;
                    }
                }
            }
        }
        let Some(longest_range) = ranges.iter().sorted_by(|x, y| y.len().cmp(&x.len())).next()
        else {
            return ret;
        };
        log::debug!("[get_read_counts] longest_range {:?}", longest_range);

        let mid = (longest_range.end + longest_range.start) / 2;
        let nstart = std::cmp::max(mid - 1, longest_range.start);
        let nend = std::cmp::min(mid + 1, longest_range.end);
        log::debug!("[get_read_counts] nstart {nstart} nend {nend}");
        if nend != nstart {
            for (path_idx, path) in paths.iter().enumerate() {
                let matches = hap2read.get(&(path_idx as u32)).unwrap();
                let l_reads = matches
                    .iter()
                    .filter(|x| {
                        x.path[nstart as usize..nend as usize]
                            .iter()
                            .any(|b| *b != b'x')
                    })
                    .count();
                ret.insert(path.vstr(), l_reads as i32);
            }
        }
        ret
    }

    /// Marks which reads match against which haplotypes.
    /// To reduce copying, these are integers corresponding into offsets into `assembled_haps`.
    fn get_read_support<'a>(
        haps_to_reads: &'a BTreeMap<VStr<'a>, Vec<&ReadAlignmentId>>,
        assembled_haps: &'a [VString],
        support: &'a ReadHapSupport<'a>,
    ) -> Result<ReadSupport<'a>, DError> {
        let uniquely_supporting_haps = &support.support_reads;
        let mut uniquely_supporting_reads = assembled_haps
            .iter()
            .map(|hap| (VStr::from(hap), vec![]))
            .collect::<BTreeMap<_, _>>();

        // 1. get_read_counts
        let read_counts = Self::get_read_counts(uniquely_supporting_haps, assembled_haps);
        for (haplotype_index, reads) in uniquely_supporting_haps.iter() {
            // XC - The function get_read_support doesn't seem to be returning assembled haplotypes and their supporting reads.
            // This line here seems to be adding partial haplotypes (those represented by single reads) to uniquely_supporting_reads.
            let ass_hap = VStr::from(&assembled_haps[*haplotype_index as usize]);
            let read_support = uniquely_supporting_reads.entry(ass_hap).or_default();
            for read in reads.iter() {
                let read_hap = read.path.vstr();
                for name in &haps_to_reads[&read_hap] {
                    read_support.push(VStr::from(&name.read_name[..]));
                }
            }
            read_support.sort();
            read_support.dedup();
        }

        // 2. Update uniquely_supporting_reads with haplotypes_to_reads
        // Handle non-unique support
        let mut nonuniquely_supporting_reads = BTreeMap::<VStr<'a>, &'a [u32]>::new();
        for (read, hap_ids) in support
            .supporting_haps_per_read
            .iter()
            .filter(|(_read, haps)| haps.len() > 1)
        {
            nonuniquely_supporting_reads.insert(VStr::from(&read.name), &hap_ids[..]);
        }
        Ok((
            uniquely_supporting_reads,
            nonuniquely_supporting_reads,
            read_counts,
        ))
    }

    /// Label genotype assignments based on genotypes.
    fn get_genotype_in_hap(
        &self,
        var_reads: &BTreeMap<String, Assignment>,
        hap_reads: &[String],
        hap_reads_nonunique: &BTreeMap<&String, &Vec<VString>>,
    ) -> Result<Genotype, DError> {
        let mut hap_reads_contain_var = hap_reads
            .iter()
            .filter(|x| var_reads.contains_key(&x[..]))
            .collect::<Vec<_>>();
        if hap_reads_contain_var.len() < 3 {
            for name in hap_reads_nonunique
                .keys()
                .filter(|x| var_reads.contains_key(&x[..]))
            {
                hap_reads_contain_var.push(name);
            }
        }
        // Convert back to status
        let hap_reads_contain_var = hap_reads_contain_var
            .iter()
            .map(|x| var_reads.get(&x[..]).copied().unwrap())
            .collect::<Vec<Assignment>>();

        if hap_reads_contain_var.len() >= 3 {
            let ctr = hap_reads_contain_var
                .into_iter()
                .collect::<counter::Counter<_, i32>>();
            let most_common = ctr.k_most_common_ordered(2);
            let top = most_common[0];
            let top_threshold = top.1 as f32 * 0.15;
            let threshold = if top_threshold > 2. {
                2.
            } else {
                top_threshold
            };
            if most_common.len() == 1 || most_common[1].1 as f32 <= threshold {
                match top.0 {
                    Assignment::Alt => {
                        return Ok(Genotype::One);
                    }
                    Assignment::Ref => {
                        return Ok(Genotype::Zero);
                    }
                    _ => {}
                }
            }
        }
        Ok(Genotype::Dot)
    }

    /// Check variants in haps
    pub fn check_variants_in_haplotypes(
        &self,
        site: &CandidateSite,
        ref_seq: &[u8],
        min_bq: Option<u8>,
    ) -> Result<BTreeMap<String, Assignment>, DError> {
        let mut dreads = BTreeMap::<String, Assignment>::new();
        let var_size = site.var_seq.len() as i64 - site.ref_seq.len() as i64;
        let indel_base_in_read = match var_size.cmp(&0) {
            Ordering::Less => {
                format!(
                    "{}{}{}",
                    std::str::from_utf8(&[site.ref_seq[0]])?.to_string(),
                    var_size,
                    VStr::from(&site.ref_seq[1..])
                )
            }
            Ordering::Greater => {
                format!(
                    "{}+{}{}",
                    std::str::from_utf8(&[site.ref_seq[0]])?.to_string(),
                    var_size,
                    VStr::from(&site.var_seq[1..])
                )
            }
            _ => String::default(),
        };
        //log::trace!("indel_base_in_read: {indel_base_in_read}");
        let ref_base = site.ref_seq[0];
        let alt_base = if var_size == 0 {
            VStr::from(&site.var_seq[..1])
        } else {
            VStr::from(&indel_base_in_read[..])
        };
        let mut reader = self.realigned_bam();
        let tid = self.genome_tid().map(|x| x as i32).expect("No chr tid");
        reader.fetch((tid, site.pos - 1, site.pos + 1))?;
        for pile in reader.pileup() {
            let pile = pile?;
            // If before position, skip.
            // If after position, break.
            // Otherwise, handle.
            match pile.pos().cmp(&(site.pos as u32)) {
                Ordering::Less => {
                    continue;
                }
                Ordering::Equal => {}
                Ordering::Greater => {
                    break;
                }
            }
            let min_base_quality = if min_bq.is_none() {
                self.settings
                    .site_selection_settings
                    .min_candidate_base_quality
            } else {
                min_bq.unwrap()
            };
            let offset = self.offset() as usize;
            for aln in pile.alignments() {
                // let qpos = aln.qpos();
                let query_pos_raw = raw_qpos(&aln);
                /*
                if aln.is_refskip() {
                    log::trace!("query_pos: {qpos:?}. Raw: {query_pos_raw:?}");
                    continue;
                }
                */
                let record = aln.record();
                let is_reverse = record.is_reverse();
                let refskip_char = if is_reverse { b'<' } else { b'>' };
                let _qual = record.qual();
                let bq = base_qual(&aln);
                if bq < min_base_quality {
                    //log::trace!(
                    //    "Skipping position for base quality {bq} < threshold {min_base_quality} at {query_pos_raw}. Quality: {qual:?}",
                    //);
                    continue;
                }
                //log::trace!("Using position for bq {bq}");
                let record = aln.record();
                let seq = record.seq().as_bytes();
                let mut query_seq = VString::default();
                let base = if !aln.is_del() {
                    let base = seq.get(query_pos_raw).copied().unwrap_or(b'N');
                    site_selection::maybe_strand_mark_char(base, is_reverse, false)
                } else if aln.is_refskip() {
                    refskip_char
                } else {
                    b'*'
                };
                query_seq.push(base);
                let pos = pile.pos() as usize;
                /*
                log::trace!(
                    "pos: {pos}. qpos: {query_pos_raw}. query_seq: {} bases, {query_seq} seq",
                    query_seq.len(),
                );
                */
                match aln.indel() {
                    Indel::Ins(x) => {
                        query_seq.push(b'+');
                        query_seq.extend_from_slice(x.to_string().as_bytes());
                        for j in 1..=(x as usize) {
                            query_seq.push(site_selection::maybe_strand_mark_char(
                                seq[j + query_pos_raw],
                                is_reverse,
                                false,
                            ));
                        }
                    }
                    Indel::Del(x) => {
                        query_seq.push(b'-');
                        query_seq.extend_from_slice(x.to_string().as_bytes());
                        for j in 1..=(x as usize) {
                            query_seq.push(site_selection::maybe_strand_mark_char(
                                ref_seq[j + pos - offset],
                                is_reverse,
                                false,
                            ));
                        }
                    }
                    bam::pileup::Indel::None => {}
                }
                let asn = match (query_seq == ref_base, query_seq == alt_base) {
                    (true, _) => Assignment::Ref,
                    (_, true) => Assignment::Alt,
                    _ => Assignment::Dot,
                };
                /*
                log::trace!(
                    "query_seq {} ref_base {} alt_base {}",
                    query_seq.to_string(),
                    std::str::from_utf8(&[ref_base])?.to_string(),
                    VString::from(alt_base).to_string(),
                );
                */
                let this_read_names = self.get_read_names(&record, None);
                let this_read_name = this_read_names
                    .first()
                    .ok_or("this_read_names is empty")?
                    .to_string();
                dreads.insert(this_read_name, asn);
            }
        }

        Ok(dreads)
    }

    ///get boundaries of (partial) haplotypes
    fn get_hap_variant_ranges(&self, hap: VStr<'_>) -> range::I64 {
        let (start, end) = get_start_end(&hap);
        let range = range::I64::new(start.into(), end.into());
        let nstart_previous_pos: i64 = if range.start == 0 {
            self.left_boundary_0based()
        } else {
            // get site before
            self.het_sites[(range.start - 1) as usize].pos + 1
        };
        let nend_next_pos = if end == (hap.len() - 1) as i32 {
            self.right_boundary_0based()
        } else {
            // get site after
            self.het_sites[(end + 1) as usize].pos - 1
        };
        range::I64::new(nstart_previous_pos, nend_next_pos)
    }

    /// Find corresponding coordinates in the secondary region.
    pub fn get_range_in_other_gene(&self, pos: i64, search_range: Option<i64>) -> Option<i64> {
        let search_range = search_range.unwrap_or(200);
        self.matches.get(&pos).copied().or_else(|| {
            (pos..pos + search_range)
                .filter_map(|x| self.matches.get(&x))
                .copied()
                .next()
        })
    }

    /// Given a haplotype, get its 5p clip position
    pub fn get_5pclip_from_hap(&self, hap: &VStr) -> Result<Option<i64>, DError> {
        let het_sites = &self.het_sites;
        let hap_len = hap.iter().len();
        assert_eq!(het_sites.len(), hap_len);
        let mut clips_not_present = Vec::new();
        for (index, base) in enumerate(hap.iter()) {
            if index < hap_len - 1 {
                let next_base = hap.get(index + 1).ok_or("index error")?;
                let site_before = het_sites.get(index).ok_or("index error")?.pos;
                let site_after = het_sites.get(index + 1).ok_or("index error")?.pos;
                for clip_position in &self.clip_5p_positions {
                    if *clip_position > site_before && *clip_position < site_after {
                        if *base == b'0' && *next_base != b'0' && *next_base != b'x' {
                            return Ok(Some(*clip_position));
                        }
                        if *next_base != b'0' && *next_base != b'x' {
                            if (*base != b'0' && *base != b'x')
                                || (*base == b'x' && site_after - *clip_position < 5000)
                            {
                                clips_not_present.push(*clip_position);
                            }
                        }
                    }
                }
            }
        }
        if clips_not_present == self.clip_5p_positions {
            return Ok(Some(0));
        }
        return Ok(None);
    }

    /// Given a haplotype, get its 3p clip position
    pub fn get_3pclip_from_hap(&self, hap: &VStr) -> Result<Option<i64>, DError> {
        let het_sites = &self.het_sites;
        let hap_len = hap.iter().len();
        assert_eq!(het_sites.len(), hap_len);
        let mut clips_not_present = Vec::new();
        for (index, base) in enumerate(hap.iter()) {
            if index < hap_len - 1 {
                let next_base = hap.get(index + 1).ok_or("index error")?;
                let site_before = het_sites.get(index).ok_or("index error")?.pos;
                let site_after = het_sites.get(index + 1).ok_or("index error")?.pos;
                for clip_position in &self.clip_3p_positions {
                    if *clip_position > site_before && *clip_position < site_after {
                        if *next_base == b'0' && *base != b'0' && *base != b'x' {
                            return Ok(Some(*clip_position));
                        }
                        if *base != b'0' && *base != b'x' {
                            if (*next_base != b'0' && *next_base != b'x')
                                || (*next_base == b'x' && *clip_position - site_before < 5000)
                            {
                                clips_not_present.push(*clip_position);
                            }
                        }
                    }
                }
            }
        }
        if clips_not_present == self.clip_3p_positions {
            return Ok(Some(0));
        }
        return Ok(None);
    }

    ///
    /// Summarize variants per hap.
    /// Output variants + genotypes.
    /// Haps may vary in length, range is reported.
    pub fn output_variants_in_haps(
        &mut self,
        result: &PhasedResult,
        known_del: &BTreeMap<char, String>,
        assembled_haps: BTreeMap<VStr<'_>, String>,
    ) -> Result<BTreeMap<String, HapInfo>, DError> {
        use std::str::FromStr;

        if result.assemblies.main_haps.is_empty() {
            return Ok(BTreeMap::new());
        }
        let het_sites = self.het_sites.clone(); // TODO: consider using a smart pointer instead of copying out for lifetime management.
        let no_phasing_sites = self.het_sites_no_phasing.clone();
        log::debug!("self.het_sites_no_phasing {:?}", self.het_sites_no_phasing);
        let mut hap_info = BTreeMap::<String, HapInfo>::new();
        let mut hap_variants = assembled_haps
            .values()
            .cloned()
            .map(|x| (x, BTreeSet::new()))
            .collect::<BTreeMap<String, BTreeSet<String>>>();

        let ref_seq = {
            let faidx = self.make_faidx()?;
            let (chrom, start, stop) = self
                .parsed_nchr_0based()
                .expect("Malformatted region string");

            faidx
                .fetch_seq(chrom, start as usize, stop as usize)?
                .to_owned()
        };
        if !result.uniquely_supporting_reads.is_empty() {
            log::trace!("There are unique reads");
            for var in &no_phasing_sites {
                let mut genotypes = vec![];
                let var_reads = self.check_variants_in_haplotypes(var, &ref_seq, None)?;
                log::trace!(
                    "Checked vars in haps for {var:?}. Current read haps: {:?}",
                    result
                        .uniquely_supporting_reads
                        .keys()
                        .collect::<counter::Counter<_, i32>>()
                );
                let mut haps_with_variant = Vec::new();
                for (hap, hap_name) in &assembled_haps {
                    let hap_vstring: VString = hap.into();
                    let hap_reads = &result
                        .uniquely_supporting_reads
                        .get(&hap_vstring)
                        .unwrap_or_else(|| {
                            panic!(
                                "hap missing {hap:?} in set {:?}",
                                result.uniquely_supporting_reads
                            )
                        })[..];
                    let hap_reads_nonunique = result
                        .nonuniquely_supporting_reads
                        .iter()
                        .filter(|(_read, hap_set)| hap_set.contains(&hap.into()))
                        .collect::<BTreeMap<_, _>>();
                    let genotype =
                        self.get_genotype_in_hap(&var_reads, hap_reads, &hap_reads_nonunique)?;
                    genotypes.push(genotype);
                    if genotype == Genotype::One {
                        haps_with_variant.push(hap_name);
                    }
                }
                if haps_with_variant.is_empty() {
                    self.het_sites_no_phasing.retain(|x| x != var); // Remove var
                } else {
                    for hap_name in haps_with_variant {
                        hap_variants
                            .get_mut(hap_name)
                            .expect("hap name missing from hap_variants")
                            .insert(var.to_string());
                    }
                }
            }
        }
        log::trace!("[output_variants_in_haps] finished unique reads. Now het sites not used.");

        // And now het sites
        for (hap, hap_name) in &assembled_haps {
            let mut hap_bounds = self.get_hap_variant_ranges(*hap);
            let mut is_truncated = Vec::new();
            // Get support rates.
            for (base, het_site) in hap.iter().zip(het_sites.iter()) {
                if *base == b'2' {
                    hap_variants
                        .entry(hap_name.clone())
                        .or_default()
                        .insert(het_site.to_string());
                } else if let Some(del_name) = known_del.get(&(*base as char)) {
                    hap_variants.entry(hap_name.clone()).or_default();
                    if !hap_variants[hap_name].contains(del_name) {
                        hap_variants
                            .entry(hap_name.clone())
                            .or_default()
                            .insert(del_name.clone());
                    }
                }
            }
            // Handle dels.
            let mut filtered_hom = self.hom_sites.clone();
            let del_names = known_del
                .values()
                .take(self.settings.max_number_deletions as usize)
                .collect::<Vec<_>>();
            for (del_name, del_data) in del_names.iter().zip(self.del_data.iter()) {
                if hap_variants[hap_name].contains(*del_name) {
                    filtered_hom.retain(|x| {
                        x.pos < del_data.threep().start || x.pos > del_data.fivep().end
                    });
                }
            }
            // Handle clips.
            let clip_position_5p = self.get_5pclip_from_hap(hap)?;
            if let Some(clip_position_5p_value) = clip_position_5p {
                if clip_position_5p_value != 0 {
                    filtered_hom.retain(|x| x.pos > clip_position_5p_value);
                    hap_bounds.start = std::cmp::max(hap_bounds.start, clip_position_5p_value);
                    hap_variants
                        .entry(hap_name.clone())
                        .or_default()
                        .insert(format!("{}_clip_5p", clip_position_5p_value + 1));
                    if clip_position_5p_value + 1 > self.gene_start() {
                        is_truncated.push(String::from("5p"));
                    }
                }
            }
            let clip_position_3p = self.get_3pclip_from_hap(hap)?;
            if let Some(clip_position_3p_value) = clip_position_3p {
                if clip_position_3p_value != 0 {
                    filtered_hom.retain(|x| x.pos < clip_position_3p_value);
                    hap_bounds.end = std::cmp::min(hap_bounds.end, clip_position_3p_value);
                    hap_variants
                        .entry(hap_name.clone())
                        .or_default()
                        .insert(format!("{}_clip_3p", clip_position_3p_value + 1));
                    if clip_position_3p_value + 1 < self.gene_end() {
                        is_truncated.push(String::from("3p"));
                    }
                }
            }

            log::trace!("check hom sites");
            if !result.uniquely_supporting_reads.is_empty() {
                for var in &filtered_hom {
                    let mut genotypes = vec![];
                    let var_reads = self.check_variants_in_haplotypes(var, &ref_seq, None)?;
                    log::trace!(
                        "Checked vars in haps for {var:?}. Current read haps: {:?}",
                        result
                            .uniquely_supporting_reads
                            .keys()
                            .collect::<counter::Counter<_, i32>>()
                    );
                    let mut haps_with_variant = Vec::new();
                    for (hap, hap_name) in &assembled_haps {
                        let hap_vstring: VString = hap.into();
                        let hap_reads = &result
                            .uniquely_supporting_reads
                            .get(&hap_vstring)
                            .unwrap_or_else(|| {
                                panic!(
                                    "hap missing {hap:?} in set {:?}",
                                    result.uniquely_supporting_reads
                                )
                            })[..];
                        let hap_reads_nonunique = result
                            .nonuniquely_supporting_reads
                            .iter()
                            .filter(|(_read, hap_set)| hap_set.contains(&hap.into()))
                            .collect::<BTreeMap<_, _>>();
                        let genotype =
                            self.get_genotype_in_hap(&var_reads, hap_reads, &hap_reads_nonunique)?;
                        genotypes.push(genotype);
                        if genotype == Genotype::One {
                            haps_with_variant.push(hap_name);
                        }
                    }
                    for hap_name in haps_with_variant {
                        hap_variants
                            .get_mut(hap_name)
                            .expect("hap name missing from hap_variants")
                            .insert(var.to_string());
                    }
                }
            } else {
                for site in &filtered_hom {
                    hap_variants
                        .entry(hap_name.clone())
                        .or_default()
                        .insert(site.to_string());
                }
            }
            hap_bounds.start = std::cmp::max(hap_bounds.start, self.left_boundary_0based());
            hap_bounds.end = std::cmp::min(hap_bounds.end, self.right_boundary_0based());
            // Handle gene boundaries
            let boundary_gene2 = if self
                .locus_config()
                .gene2_region(self.settings.genome == "37")
                .is_some()
            {
                let (start, end) = [hap_bounds.start, hap_bounds.end]
                    .into_iter()
                    .map(|x| self.get_range_in_other_gene(x, None))
                    .next_tuple()
                    .unwrap();
                match (start, end) {
                    (Some(start), Some(end)) => Some(range::I64::new(
                        std::cmp::min(start, end),
                        std::cmp::max(start, end),
                    )),
                    (None, Some(end)) => Some(range::I64::new(-2, end)),
                    (Some(start), None) => Some(range::I64::new(start, -2)),
                    (None, None) => None,
                }
            } else {
                None
            };
            let var_tmp = hap_variants.get(hap_name).unwrap();
            let getpos = |x: &str| {
                x.split_terminator('_')
                    .next()
                    .and_then(|x| x.parse::<i64>().ok())
            };
            // variant position so far is 1-based. hap bounds are 0-based.
            let var_tmp = var_tmp
                .iter()
                .filter(|x| {
                    let pos = getpos(x).unwrap();
                    pos >= hap_bounds.start + 1 && pos <= hap_bounds.end + 1
                })
                .sorted_by_cached_key(|x| getpos(&x[..]))
                .collect::<Vec<_>>();
            let mut variants = Vec::with_capacity(var_tmp.len());
            // variant position now is 0-based: CandidateSite::from_str(x)
            for var in var_tmp.into_iter().map(|x| CandidateSite::from_str(x)) {
                variants.push(var?);
            }
            let boundary = hap_bounds;
            let info = HapInfo {
                variants,
                boundary,
                boundary_gene2,
                is_truncated,
            };
            hap_info.insert(hap_name.into(), info);
        }

        Ok(hap_info)
    }
}
