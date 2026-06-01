use crate::config::Locus as LocusConfig;
use crate::detail::deletion::Datum as DeletionDatum;
use crate::detail::range;
use crate::detail::util::{self, DError, HashMap};
use crate::phaser::Exception;
use crate::phaser::Phaser;
use itertools::Itertools;
use rust_htslib::bam::{self, IndexedReader, Read as BamRead};
use std::fmt;
use vstr::{VStr, VString};

use std::collections::{BTreeMap, BTreeSet};
use std::hash::BuildHasherDefault;

/// `CandidateSite`
/// Holds position (`i64`), `ref_seq` and `var_seq` (both `VString`).
///
///
/// Use 0-based indexing internally but displays as 1-based to be consistent with SAM indexing.
///
/// For deletion entries, stored as `del` for `ref_seq` and `num_bases_deleted` as `var_seq`.
#[derive(Clone, Hash, PartialEq, PartialOrd, Ord, Eq, serde::Serialize, serde::Deserialize)]
pub struct CandidateSite {
    pub pos: i64,
    pub ref_seq: VString,
    pub var_seq: VString,
}

impl fmt::Display for CandidateSite {
    fn fmt(&self, format: &mut fmt::Formatter) -> fmt::Result {
        write!(format, "{}_{}_{}", self.pos + 1, self.ref_seq, self.var_seq)
    }
}
impl fmt::Debug for CandidateSite {
    fn fmt(&self, format: &mut fmt::Formatter) -> fmt::Result {
        write!(format, "{self}")
    }
}

impl CandidateSite {
    /// Create a `CandidateSite` from an integer, a reference sequence, and a variant sequence.
    /// # Panics
    /// Panics if `pos: impl TryInto<i64>` fails to convert to `i64`
    #[must_use]
    pub fn new(
        pos: impl TryInto<i64>,
        ref_seq: impl Into<VString>,
        var_seq: impl Into<VString>,
    ) -> Self {
        let ref_seq = ref_seq.into();
        let var_seq = var_seq.into();
        let pos = pos
            .try_into()
            .unwrap_or_else(|_| panic!("Failed to convert pos to i64"));
        Self {
            pos,
            ref_seq,
            var_seq,
        }
    }
    #[must_use]
    pub fn is_unit_length(&self) -> bool {
        self.ref_seq.len() == 1 && self.var_seq.len() == 1
    }

    #[must_use]
    pub fn reference_length(&self) -> usize {
        self.ref_seq.len()
    }
}

impl std::convert::From<&DeletionDatum> for CandidateSite {
    fn from(del: &DeletionDatum) -> Self {
        let pos = del.raw.start;
        let ref_seq = VString::from("del");
        let var_seq = VString::from(del.raw.len().to_string());
        CandidateSite {
            pos,
            ref_seq,
            var_seq,
        }
    }
}

impl std::str::FromStr for CandidateSite {
    type Err = Exception;
    fn from_str(x: &str) -> Result<Self, Self::Err> {
        let mut it = x.split_terminator('_');
        let (pos, ref_seq, var_seq) = it.next_tuple().ok_or_else(|| {
            Exception::new(String::from(
                "Malformatted CandidateSite string. Expected {pos}_{ref_seq}_{var_seq}",
            ))
        })?;
        let pos = pos.parse::<i64>().map_err(|e| {
            Exception::new(format!(
                "Failed to parse position (\"{pos}\") as i64. Parse error: {e:?}"
            ))
        })? - 1; // Correct for one-based indexing.
        Ok(Self {
            pos,
            ref_seq: ref_seq.into(),
            var_seq: var_seq.into(),
        })
    }
}

/// Settings for generating candidate sites for inclusion in the wfa graph.
///
/// The only major difference with paraphase at this point is that `max_indel_size` is stored as `25 - 1` and we check for `size <= threshold`.
/// The behavior will match, but we use <= instead of < because we chose a named threshold instead of a hard-coded number.
#[derive(Clone, Debug)]
pub struct Settings {
    pub min_vaf: f64,
    pub min_read_support: i32,
    // pub regions_to_check: Vec<range::I64>,
    pub min_base_quality: u8, // Min base quality for use in fingerprinting.
    pub min_candidate_base_quality: u8, // Min base quality for use in raw pileups.
    pub min_mean_base_quality: u8, // Min mean base quality for windows for inclustion.
    pub max_indel_size: i32,
    pub max_candidate_seqs: i32,
    pub trusted_read_support: i32,
    pub permit_list: BTreeMap<i64, VString>, /* a list of positions which are permitted even if filters fail. */
    pub targeted: bool,
}

impl std::default::Default for Settings {
    fn default() -> Self {
        Self {
            min_vaf: 0.11,
            min_read_support: 5,
            // regions_to_check: vec![],
            min_base_quality: 25, // Match the default in paraphase.phaser.Phaser
            min_mean_base_quality: 25,
            min_candidate_base_quality: 13, // Match the default in pysam.Pileup
            max_indel_size: 24,
            max_candidate_seqs: 3,
            trusted_read_support: 20,
            permit_list: BTreeMap::new(),
            targeted: false,
        }
    }
}

impl Settings {
    /// Determines if an indel is small enough to use.
    #[must_use]
    pub fn indel_size_passes(&self, size: usize) -> bool {
        size as i32 <= self.max_indel_size
    }

    /// Update `site_selection::Settings` from `LocusConfig` from yaml.
    ///
    /// # Panics
    /// 1. `white_list` or `permit_list` field is not a mapping as expected.
    /// 2. key from permit list is not integral as expected.
    /// 3. value from permit list is not a string as expected.
    pub fn update_from_settings(&mut self, config: &LocusConfig) {
        if let Some(list) = config.get("white_list").or(config.get("permit_list")) {
            let list = list.as_mapping().expect("permit_list should be a mapping");
            self.permit_list = list
                .into_iter()
                .map(|(k, v)| {
                    (
                        k.as_i64().expect("Key should be integral") - 1, // Subtract 1 to account for 1-based inputs.
                        VString::from(v.as_str().expect("Value should be string")),
                    )
                })
                .collect::<BTreeMap<_, _>>();
        }
    }

    /// Build `site_selection::Settings` from `LocusConfig` from yaml.
    ///
    /// # Panics
    /// 1. `white_list` or `permit_list` field is not a mapping as expected.
    /// 2. key from permit list is not integral as expected.
    /// 3. value from permit list is not a string as expected.
    #[must_use]
    pub fn new_from_settings(config: &LocusConfig) -> Self {
        let mut ret = Self::default();
        ret.update_from_settings(config);
        ret
    }
}

/*
cdef inline uint8_t strand_mark_char(uint8_t ch, bam1_t *b):
    if ch == b'=':
        if bam_is_rev(b):
            return b','
        else:
            return b'.'
    else:
        if bam_is_rev(b):
            return tolower(ch)
        else:
            return toupper(ch)
*/

#[inline]
#[must_use]
pub fn maybe_strand_mark_char(x: u8, is_rev: bool, mark_strand: bool) -> u8 {
    if mark_strand {
        strand_mark_char(x, is_rev)
    } else {
        x.to_ascii_uppercase()
    }
}

#[inline]
#[must_use]
pub fn strand_mark_char(x: u8, is_rev: bool) -> u8 {
    if is_rev {
        x.to_ascii_lowercase()
    } else {
        x.to_ascii_uppercase()
    }
}

/// Calculate sequences in a pileup, including indels.
/// Reproduces functionality in `pysam.libcalignedsegment.Pileup.get_query_sequences(add_indels=True)`.
///
/// Runtime improvement options: vectorize, remove array access checks.
///
///
/// # Arguments
///
/// * `x` - a Pileup to work with.
/// * `ref_seq` - reference sequence.
/// * `settings` - Settings, which specifies min base quality. paraphase's default is 13, but it is hidden inside pysam.
/// * `mark_strand` - whether or not to mark inserted sequences with strand. paraphase maps the data to upper-case.
///                   if false, all bases are upper-cased.
#[must_use]
pub fn query_seq_counter(
    x: &bam::pileup::Pileup,
    ref_seq: &[u8],
    settings: &Settings,
    offset: i64,
    aln2seq: &mut HashMap<String, VString>,
) -> BTreeMap<VString, i32> {
    let mut ret = BTreeMap::<VString, i32>::new();
    let min_base_quality = settings.min_candidate_base_quality;
    let pos = x.pos();
    let offset = offset as usize;
    for aln in x.alignments() {
        let query_pos_raw = util::raw_qpos(&aln);

        /*
        // pysam isn't filtering these, so we are leaving this in.
        if aln.is_refskip() {
            log::warn!("Refskip found in pileup. Genomic data expected. query_pos: {query_pos:?}. Raw: {query_pos_raw:?}. Position: {pos}/{}", pos + offset);
            continue;
        }
        let refskip_char = if is_reverse { b'<' } else { b'>' };
        */

        let record = aln.record();
        let qname = VStr::from(record.qname());
        let is_reverse = record.is_reverse();
        let bq = record.qual().get(query_pos_raw).copied().unwrap_or(0);
        if bq < min_base_quality {
            //log::trace!("Skipping read {qname} at pos {pos} for base quality {bq} < threshold {min_base_quality}",);
            continue;
        }
        // Cache query seq, as this is can be shared across thousands of sites.
        let entry = aln2seq
            .entry(format!("{qname}:{}:{}", record.pos(), record.flags()))
            .or_insert_with(|| record.seq().as_bytes().into());
        let seq = &entry[..];
        let mut query_seq = VString::default();
        let base = if !aln.is_del() && !aln.is_refskip() {
            seq.get(query_pos_raw).copied().unwrap_or(b'N')
        } else if aln.is_refskip() {
            if is_reverse {
                b'<'
            } else {
                b'>'
            }
        } else {
            b'*'
        };
        // is_del
        /*
            if p.is_refskip:
                if bam_is_rev(p.b):
                    kputc(b'<', buf)
                else:
                    kputc(b'>', buf)
            else:
                kputc(b'*', buf)

        */
        query_seq.push(base);
        let pos = pos as usize;
        //log::trace!(
        //    "{qname}@pos={pos}. qpos: {query_pos_raw}. query_seq: {} bases, {query_seq} seq",
        //    query_seq.len(),
        //);
        match aln.indel() {
            bam::pileup::Indel::Ins(x) => {
                debug_assert!(x > 0);
                query_seq.push(b'+');
                query_seq.extend_from_slice(x.to_string().as_bytes());
                // TODO: speed this up by using slice operations to copy out faster.
                for j in 1..=(x as usize) {
                    query_seq.push(seq[j + query_pos_raw]);
                }
            }
            bam::pileup::Indel::Del(x) => {
                debug_assert!(x > 0);
                query_seq.push(b'-');
                query_seq.extend_from_slice(x.to_string().as_bytes());
                for j in 1..=(x as usize) {
                    query_seq.push(ref_seq[j + pos - offset]);
                }
            }
            bam::pileup::Indel::None => {}
        }
        query_seq.make_ascii_uppercase();
        //log::trace!(
        //    "[query_seq_counter] Inserting query {query_seq} at pos {pos} for read name {qname} with flag {}",
        //    record.flags()
        //);
        *ret.entry(query_seq).or_default() += 1;
    }
    ret
}

pub type RawVariantCounts = HashMap<i64, BTreeMap<VString, i32>>;
pub type RawVariantCountsString = HashMap<i64, BTreeMap<String, i32>>;

#[must_use]
pub fn raw_variants_to_string(
    input: &HashMap<i64, BTreeMap<VString, i32>>,
) -> HashMap<i64, BTreeMap<String, i32>> {
    input
        .iter()
        .map(|(k, v)| {
            (
                *k,
                (v.iter()
                    .map(|(hap, count)| (hap.to_string(), *count))
                    .collect::<BTreeMap<_, _>>()),
            )
        })
        .collect::<HashMap<_, _>>()
}

/// Generates a multiset of query sequences at each position in the provided range.
/// # Errors
/// 1. Failure to seek `chr:start-stop` in `IndexedReader`e
fn position_seq_counter(
    x: &mut IndexedReader,
    chr: &str,
    positions: &range::I64,
    ref_seq: &[u8],
    settings: &Settings,
    offset: i64,
    aln2seq: &mut HashMap<String, VString>,
) -> Result<RawVariantCounts, DError> {
    let mut res: HashMap<i64, BTreeMap<VString, i32>> =
        HashMap::with_capacity_and_hasher(positions.len(), BuildHasherDefault::default());
    log::trace!(
        "Fetching {chr}:{}-{}. chrnames: {:?}",
        positions.start,
        positions.end,
        x.header()
            .target_names()
            .into_iter()
            .map(VStr::from)
            .collect::<Vec<_>>()
    );
    let target_names = x
        .header()
        .target_names()
        .into_iter()
        .map(VString::from)
        .collect::<Vec<_>>();
    let chr_bytes = &chr.as_bytes();
    let tid = target_names
        .iter()
        .position(|x| &&x[..] == chr_bytes)
        .or(target_names
            .iter()
            .position(|x| &x[..] == chr.split_terminator('_').next().unwrap().as_bytes()))
        .unwrap_or_else(|| {
            panic!("Failed to find chromosome for inputs {target_names:?} and query {chr}")
        }) as i32;
    let found = target_names
        .iter()
        .find(|x| &&x[..] == chr_bytes)
        .or(target_names
            .iter()
            .find(|x| &x[..] == chr.split_terminator('_').next().unwrap().as_bytes()))
        .unwrap_or_else(|| {
            panic!("Failed to find chromosome for inputs {target_names:?} and query {chr}")
        });
    log::trace!(
        "Piling up with tid = {tid}/{found} from chr {chr}:{}-{} and targets {target_names:?}",
        positions.start,
        positions.end
    );
    /*
    assert_eq!(
        x.header().target_names().len(),
        1,
        "Making sure there is only one allows us to just use tid = 0"
    );
    */
    let mut num_used = 0usize;
    let start = positions.start + 1;
    let end = positions.end + 1;
    //("Paraphase uses 1-based coordinates but queries 0-based in pysam. We add an offset (1) to use the same coordinates and get the same positions. This may be in error.");
    x.fetch((tid, start, end))?;
    for x in x.pileup() {
        let x = x?;
        let pos = i64::from(x.pos());
        if pos < start {
            continue;
        }
        if pos >= end {
            break;
        }
        num_used += 1;
        let seqs = query_seq_counter(&x, ref_seq, settings, offset, aln2seq);
        res.insert(pos, seqs);
    }
    log::debug!("Used {num_used} sites");
    Ok(res)
}

pub type VariantMap = HashMap<i64, Vec<(VString, VString)>>;
pub type FilteredSitesForJson = (
    HashMap<i64, Vec<(String, String)>>,
    HashMap<i64, (String, String)>,
);

#[derive(Clone, Debug, Default, serde::Serialize, serde::Deserialize)]
pub struct FilteredSites {
    pub variants: VariantMap,
    pub variants_no_phasing: HashMap<i64, (VString, VString)>,
}

impl FilteredSites {
    #[must_use]
    pub fn to_json_friendly(&self) -> FilteredSitesForJson {
        let variants = self
            .variants
            .iter()
            .map(|(key, vec)| {
                (
                    *key,
                    vec.iter()
                        .map(|(x, y)| (x.to_string(), y.to_string()))
                        .collect::<Vec<_>>(),
                )
            })
            .collect::<HashMap<_, _>>();
        let variants_no_phasing = self
            .variants_no_phasing
            .iter()
            .map(|(key, (x, y))| (*key, (x.to_string(), y.to_string())))
            .collect::<HashMap<_, _>>();
        (variants, variants_no_phasing)
    }
}

///```
/// use paraphase::detail::site_selection::raw_pile_to_string;
/// use vstr::VString;
/// let mut pile =
///     paraphase::detail::util::HashMap::<i64, std::collections::BTreeMap<VString, i32>>::default(
///     );
/// assert_eq!(raw_pile_to_string(&pile).unwrap(), "{}");
/// pile.insert(
///     0,
///     [(VString::from("ACGT"), 2)]
///         .into_iter()
///         .collect::<std::collections::BTreeMap<VString, i32>>(),
/// );
/// assert_eq!(raw_pile_to_string(&pile).unwrap(), "{0: {\"ACGT\": 2}\n}");
/// pile.entry(1).or_default().insert(VString::from("AGTG"), 1);
/// assert_eq!(
///     raw_pile_to_string(&pile).unwrap(),
///     "{0: {\"ACGT\": 2},\n1: {\"AGTG\": 1}\n}"
/// );
/// ```
pub fn raw_pile_to_string(
    x: &HashMap<i64, BTreeMap<VString, i32>>,
) -> Result<String, std::fmt::Error> {
    let x = x
        .iter()
        .map(|(k, v)| {
            (
                k,
                v.iter()
                    .map(|(vkey, vval)| (vkey.vstr(), vval))
                    .collect::<BTreeMap<_, _>>(),
            )
        })
        .collect::<BTreeMap<_, _>>();
    use std::fmt::Write;
    let mut ret = String::with_capacity(x.len() * 10);
    ret.push('{');
    let num_x = x.len();
    for (x_id, (pos, items)) in x.iter().enumerate() {
        write!(ret, "{pos}:")?;
        ret.push_str(" {");
        ret.push_str(
            &Itertools::intersperse(
                items
                    .iter()
                    .map(|(seq, count)| format!("\"{seq}\": {count}")),
                String::from(","),
            )
            .collect::<String>(),
        );
        ret.push_str(if x_id == num_x - 1 { "}\n" } else { "},\n" });
    }
    ret.push('}');
    Ok(ret)
}

///```
/// assert!(paraphase::detail::site_selection::seq_is_indel(b"BG-"));
/// assert!(!paraphase::detail::site_selection::seq_is_indel(b"BG"));
/// assert!(!paraphase::detail::site_selection::seq_is_indel(
///     b"ACGTG_1_3"
/// ));
/// ```
#[must_use]
pub fn seq_is_indel(x: &[u8]) -> bool {
    x.iter().any(|x| matches!(x, b'+' | b'-'))
}

///
/// Filter raw sites + select variant sites.
///
/// # Errors
/// 1. Faidx fetch failure.
/// 2. Failure to make local faidx.
/// 3. Error processing indel. (Should not happen.)
fn filtered_sites(
    settings: &Settings,
    phaser: &mut Phaser,
    raw_piles: &HashMap<i64, BTreeMap<VString, i32>>,
    regions_to_check: &[range::I64],
) -> Result<FilteredSites, DError> {
    let del_str = VString::from(&[b'*'][..]);
    let offset = phaser.offset();
    log::debug!("min read support: {}", settings.min_read_support);
    //log::trace!("Initial raw pileup: {}", raw_pile_to_string(raw_piles)?);
    log::trace!(
        "hpol sites at start of get_candidate_pos: {:?}",
        phaser.low_complexity_sites
    );
    let full_ref_seq = {
        log::debug!("Getting full ref seq");
        let faidx = phaser.make_local_faidx()?;
        log::debug!(
            "Loading seq {:?} for faidx at {:?}",
            phaser.local_chr(),
            &phaser.local_reference()
        );
        for id in 0..faidx.n_seqs() {
            log::debug!("Seq {id} is {}", faidx.seq_name(id as i32)?);
        }
        faidx
            .fetch_seq(
                phaser.local_chr().expect("No local chr"),
                0,
                i32::MAX as usize,
            )?
            .to_ascii_uppercase()
    };
    log::debug!("Fetched seq of size {}", full_ref_seq.len());
    let mut variants: HashMap<i64, Vec<(VString, VString)>> =
        HashMap::with_capacity_and_hasher(8, BuildHasherDefault::default());
    let mut variants_no_phasing: HashMap<i64, (VString, VString)> =
        HashMap::with_capacity_and_hasher(8, BuildHasherDefault::default());
    let faidx = phaser.make_faidx()?; // Notice that this is global/genomic, not local.
                                      // Homozygous
                                      // May need to save more/less than 32768.
    let chr = phaser.chr().expect("Impossible");
    log::trace!(
        "Fetching seq {chr} from faidx with names {:?}",
        util::faidx_names(&faidx)
    );
    let cached_faidx_seq = faidx
        .fetch_seq(chr, 0, i64::MAX as usize)?
        .to_ascii_uppercase();
    for (pos, pileup) in raw_piles {
        let pos = *pos;
        let depth = pileup.values().sum::<i32>();
        let del_count = pileup.get(&del_str).copied().unwrap_or(0);
        log::trace!("{pos} has {del_count} deletions and total depth {depth:?}.");
        assert!(
            depth >= del_count,
            "depth {depth} should be >= del_count {del_count}"
        );
        let depth_without_dels = depth - del_count;
        log::trace!("Depth without dels at {pos} is {depth_without_dels}");
        let offset_pos = (pos - offset) as usize;
        /*
            if total_depth >= min_read_support and (
                del_bases_count < min_read_support or self.allow_del_bases(pos)
            ):
        */
        let pass_for_dels = del_count < settings.min_read_support || phaser.allow_del_bases(pos);
        if depth < settings.min_read_support || !pass_for_dels {
            log::trace!(
                "variant filtered out due to low depth or high deletion count. depth {depth}, deletion count {del_count}, min_read_support {}",
                settings.min_read_support,
            );
            continue;
        }
        let ref_seq = VString::from(full_ref_seq[offset_pos..=offset_pos].to_ascii_uppercase());
        assert!(
            ref_seq.len() == 1,
            "ref seq {ref_seq:?} should be of size 1 because of the fetch."
        );
        // We are avoiding building a counter object.
        // Technically this could be faster using a heap, but we don't expect many variants at a site and so this is fair.
        let most_common = pileup
            .iter()
            .sorted_by_key(|(_item, c)| std::cmp::Reverse(*c))
            .filter(|x| &x.0[..] != b"*")
            .map(|(x, c)| (x.vstr(), c))
            .take(settings.max_candidate_seqs as usize)
            .collect::<Vec<_>>();
        assert!(
            most_common.windows(2).all(|x| x[1].1 <= x[0].1),
            "Most common should be sorted in descending order: {most_common:?}"
        );
        let counter_len = pileup.len() - usize::from(pileup.contains_key(&del_str));
        let total_depth_above_minimum = f64::from(depth_without_dels - settings.min_read_support);
        // TODO: eventually replace with a cdf.
        // use `paraphase::util::site_probably_het`

        let depth_without_dels_threshold = f64::from(depth_without_dels) * 0.85;
        let most_common_threshold = if total_depth_above_minimum
            .partial_cmp(&depth_without_dels_threshold)
            .unwrap_or(std::cmp::Ordering::Equal)
            == std::cmp::Ordering::Greater
        {
            total_depth_above_minimum
        } else {
            depth_without_dels_threshold
        };
        log::trace!("At pos {pos}, most common threshold is {most_common_threshold} and counts are {most_common:?}. Most common: {:?}", most_common[0]);
        let is_homozygous = counter_len == 1
            || (counter_len >= 2 && f64::from(*most_common[0].1) > most_common_threshold);

        // Warning:
        // We make low complexity sites match exactly
        let hp_site_data = phaser.low_complexity_sites.get_0based(pos);
        let is_hp_site = hp_site_data.is_some();
        let allow_del = phaser.allow_del_bases(pos);
        if is_homozygous {
            log::trace!("hom at site {offset_pos}/{}", pos);
            let var_seq: VStr<'_> = most_common[0].0;
            if var_seq != ref_seq[..] {
                if !seq_is_indel(&var_seq) {
                    // Process SNV
                    // Homozygous sites in deletions are added as heterozygous sites.
                    log::trace!(
                        "Process substitution at {pos}. Var {var_seq} Ref {}. Allow del: {}",
                        VStr::from(&ref_seq[..]),
                        allow_del
                    );
                    if allow_del && del_count >= settings.min_read_support && !is_hp_site {
                        log::trace!(
                            "is not a hp site. del count: {del_count}. Allow del bases: {}",
                            allow_del
                        );
                        variants
                            .entry(pos)
                            .or_default()
                            .push((ref_seq, VString::from(var_seq)));
                    } else if hp_site_data
                        .as_ref()
                        .map_or(true, |x| !(var_seq.len() == 1 && x.contains(&var_seq[0])))
                    {
                        phaser
                            .hom_sites
                            .push(CandidateSite::new(pos, ref_seq, var_seq));
                    }
                } else if !is_hp_site {
                    // Process indel
                    log::trace!("Process indel at {pos}.");
                    let (processed_refseq, processed_varseq, indel_len) =
                        phaser.process_indel(pos, &ref_seq, &var_seq, &cached_faidx_seq)?;
                    if indel_len <= settings.max_indel_size {
                        phaser.hom_sites.push(CandidateSite::new(
                            pos,
                            processed_refseq,
                            processed_varseq,
                        ));
                    }
                }
            }
        } else if counter_len >= 2 {
            log::trace!("potential het at site {offset_pos}/{}", pos);
            let found_ref = pileup
                .iter()
                .map(|(x, c)| (x.vstr(), c))
                .any(|x| x.0[..] == ref_seq[..]);
            //let found_ref = most_common.iter().any(|x| x.0[..] == ref_seq[..]);
            if found_ref || settings.permit_list.contains_key(&pos) {
                log::trace!("found ref in most common at {pos}.");
                for (var_seq, _var_count) in most_common.iter().filter(|(seq, count)| {
                    let count = **count;
                    let is_refseq = seq[..] == ref_seq[..];
                    let sufficient_read_count = count >= settings.min_read_support;
                    let sufficient_vaf = f64::from(count) >= settings.min_vaf * f64::from(depth);
                    let count_above_trusted = count >= settings.trusted_read_support;
                    log::trace!("At pos {pos} with var = {seq} and ref = {ref_seq}. Sufficient count: {sufficient_read_count}. Sufficient vaf: {sufficient_vaf}. Count {count} vs trusted {} vs min support {}", settings.trusted_read_support, settings.min_read_support);
                    !is_refseq && ((sufficient_read_count && sufficient_vaf) || (!settings.targeted && count_above_trusted))
                }) {
                    // Substitution
                    if !seq_is_indel(var_seq) {
                        log::trace!("At pos {pos}, {var_seq} is not an indel. Is it a hp site: {is_hp_site}");
                        debug_assert!(!var_seq.iter().any(|x| matches!(x, b'+' | b'-')));
                        if is_hp_site {
                            let var_seq_prohibited = var_seq.len() == 1
                                && hp_site_data
                                    .as_ref()
                                    .map_or(false, |forbid| forbid.contains(&var_seq[0]));
                            log::trace!("At pos {pos}, {var_seq} SNV is in a homopolymer site. Is it prohibited_bases? {var_seq_prohibited}. Forbidden: {:?}", hp_site_data.as_ref().map_or(vec![], |x| x.iter().copied().collect::<Vec<_>>()));
                            log::trace!("Forbidden at neighbors: {:?}, {:?}", if pos > 0 {phaser.low_complexity_sites.get(&(pos - 1))} else {None}, phaser.low_complexity_sites.get(&(pos + 1)));
                            if !var_seq_prohibited {
                                if hp_site_data.as_ref().map_or(false, |x| x.contains(&b'1')) {
                                    log::trace!("At pos {pos}, {var_seq} SNV is in a homopolymer site. var seq is not prohibited.");
                                    variants
                                        .entry(pos)
                                        .or_default()
                                        .push((ref_seq.clone(), (*var_seq).into()));
                                } else {
                                    log::trace!("At pos {pos}, {var_seq} SNV hp site data is {hp_site_data:?}, so we add to variants_no_phasing instead of variants.");
                                    variants_no_phasing.entry(pos).or_insert_with(|| {
                                        ((&ref_seq[..]).into(), (&var_seq[..]).into())
                                    });
                                }
                            }
                        } else {
                            log::trace!("At pos {pos}, {var_seq} SNV is not in a homopolymer site, so we are adding it.");
                            variants
                                .entry(pos)
                                .or_default()
                                .push(((&ref_seq[..]).into(), (&var_seq[..]).into()));
                        }
                    } else if !is_hp_site {
                        assert!(seq_is_indel(var_seq));
                        // Indel
                        let (ref_seq, var_seq, indel_len) =
                            phaser.process_indel(pos, &ref_seq, var_seq, &cached_faidx_seq)?;
                        if indel_len <= settings.max_indel_size {
                            variants_no_phasing
                                .entry(pos)
                                .or_insert_with(|| (ref_seq, var_seq));
                        }
                    } else {
                        log::trace!("is hp site, skipping indels.");
                    }
                }
            }
        }
        // Now handle exclusion logic.
    }

    let excluded_variants = regions_to_check
        .iter()
        .flat_map(|x| {
            variants
                .keys()
                .filter(|pos| x.contains(*pos) && **pos > x.start)
        })
        .collect::<BTreeSet<_>>();
    log::trace!(
        "Excluding {}/{} variant sites. Variants: {variants:?}",
        excluded_variants.len(),
        variants.len()
    );

    for (pos, variants) in variants.iter().filter(|x| !excluded_variants.contains(x.0)) {
        if variants.len() == 1 {
            log::trace!(
                "At {pos}, we have only one variant. Candidate site. Variants: {variants:?}"
            );
            let var = &variants[0];
            phaser
                .candidate_sites
                .insert(CandidateSite::new(*pos, &var.0, &var.1));
        } else if let Some(permitted_variant) = settings.permit_list.get(pos) {
            log::trace!("At {pos}, we have a permitted variant. Checking it....");
            if let Some((rseq, vseq)) = variants
                .iter()
                .find(|(_rseq, vseq)| vseq != permitted_variant)
            {
                phaser
                    .candidate_sites
                    .insert(CandidateSite::new(*pos, rseq, vseq));
            }
        } else {
            log::trace!("At {pos}, we have more than one variant and we don't accept anything at this position. Variants: {variants:?}");
        }
    }

    let excluded_variants = regions_to_check
        .iter()
        .flat_map(|x| {
            variants_no_phasing
                .keys()
                .filter(|pos| x.contains(*pos) && **pos > x.start)
        })
        .sorted()
        .unique()
        .collect::<BTreeSet<_>>();

    log::debug!(
        "Excluding {}/{} variants_no_phasing sites.",
        excluded_variants.len(),
        variants_no_phasing.len()
    );

    for site in variants_no_phasing
        .iter()
        .filter(|x| !excluded_variants.contains(&x.0))
        .map(|(pos, variant)| CandidateSite::new(*pos, &variant.0, &variant.1))
    {
        phaser.het_sites_no_phasing.push(site);
    }

    log::debug!(
        "{} het sites no phasing from {} initial variants_no_phasing positions",
        phaser.het_sites_no_phasing.len(),
        variants_no_phasing.len()
    );
    log::debug!("{} hom_sites ", phaser.hom_sites.len());
    phaser.het_sites = phaser
        .candidate_sites
        .iter()
        .sorted()
        .cloned()
        .collect::<Vec<_>>();
    log::debug!(
        "{} candidate sites. Het sites for phasing {:?}",
        phaser.candidate_sites.len(),
        phaser.candidate_sites
    );
    Ok(FilteredSites {
        variants,
        variants_no_phasing,
    })
}

/// Performs a pileup.
/// WARNING: positions are 0-based here to match htslib.
/// We are consistently using 0-based, which does not match paraphase. We are doing book-keeping to keep it in sync outside of this function.
/// # Errors
/// 1. Failure to fetch positions in `position_seq_counter`
/// 2. Errors in `filtered_sites`.
///   a. Faidx fetch failure.
///   b. Failure to make local faidx.
///   c. Error processing indel. (Should not happen.)
pub fn pileup(
    x: &mut IndexedReader,
    positions: &range::I64,
    ref_seq: &[u8],
    phaser: &mut Phaser,
    settings: Option<&Settings>,
    regions_to_check: &[range::I64],
    chr: &str,
) -> Result<(FilteredSites, RawVariantCounts), DError> {
    let default = Settings::default();
    let settings = settings.unwrap_or(&default);
    let offset = phaser.offset();
    let mut aln2seq =
        HashMap::<String, VString>::with_capacity_and_hasher(256, BuildHasherDefault::default());
    let raw_piles =
        position_seq_counter(x, chr, positions, ref_seq, settings, offset, &mut aln2seq)?;
    log::debug!(
        "Piling up with region {chr}. local chr {:?}, global chr {:?} for bam at {:?}",
        phaser.local_chr(),
        phaser.chr(),
        phaser.realigned_bam_path(),
    );
    //log::trace!("Raw piles: {raw_piles:?}");
    let f = filtered_sites(settings, phaser, &raw_piles, regions_to_check)?;
    Ok((f, raw_piles))
}

#[cfg(test)]
mod tests {
    use super::strand_mark_char;
    #[test]
    fn strand_mark_char_ok() {
        assert_eq!(b'A', strand_mark_char(b'A', false));
        assert_eq!(b'a', strand_mark_char(b'a', true));
        assert_eq!(b'A', strand_mark_char(b'A', false));
        assert_eq!(b'a', strand_mark_char(b'a', true));
        assert_eq!(b'G', strand_mark_char(b'G', false));
        assert_eq!(b'g', strand_mark_char(b'g', true));
        assert_eq!(b'G', strand_mark_char(b'G', false));
        assert_eq!(b'g', strand_mark_char(b'g', true));
        for base in b"1234567890QRS".iter().copied() {
            assert_eq!(base.to_ascii_lowercase(), strand_mark_char(base, true));
            assert_eq!(base.to_ascii_uppercase(), strand_mark_char(base, false));
        }
    }
}
