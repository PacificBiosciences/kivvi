use crate::assembly::variant_graph_util::{median_i32, percentile_i32};
use crate::config::{
    self, region::MatchMap, Gene as GeneConfig, Locus as LocusConfig, Region as RegionConfig,
};
use crate::depth::{self, count_pos, Result as GenomeDepthResult, Sex};
use crate::detail::deletion::{BigDeletionSettings, Datum as DeletionDatum};
use crate::detail::low_complexity::LowConfidenceSites;
use crate::detail::math::depth_prob;
use crate::detail::phase_haps::{HapInfo, HapInfoForJson, PhasedResult, PhasedResultForJson};
use crate::detail::phaser_util::{fiveprime_clip_length, threeprime_clip_length};
use crate::detail::range;
use crate::detail::site_selection::{
    maybe_strand_mark_char, CandidateSite, Settings as SiteSelectionSettings,
};
use crate::detail::util::{self, DError, DResult};
use crate::io::json::{GeneCall, ReadFingerprintMap};
use crate::realign::{align_mm2_intrinsic, RealignSettings};

use itertools::{intersperse, Itertools};
use rust_htslib::{bam, bam::ext::BamRecordExtensions, bam::Read, faidx};
use vstr::{VStr, VString};

use std::cmp::Ordering;
use std::collections::{BTreeMap, BTreeSet};
use std::fmt::Debug;
use std::path::PathBuf;

pub(crate) mod defaults {
    pub const MIN_BASE_QUALITY: u8 = 25u8;
    pub const NUM_SAMPLED_DEPTH_POSITIONS: usize = 100;
    pub const PERCENTILE: i32 = 80;
    pub const CLIP_OFFSET_HAPS_FROM_READS: i64 = 30;
}

/// Convert a read fingerprint map into a form that can be displayed easily with `serde`.
#[must_use]
pub fn to_string_map(map: &ReadFingerprintMap) -> BTreeMap<String, String> {
    map.clone()
        .into_iter()
        .map(|(id, path)| (id.unique_name(), path.to_string()))
        .collect::<_>()
}

/// Sample-level settings
///
/// # Members
/// `genome_reference` - `PathBuf`. Path to genome ref
/// `genome_bam` - `PathBuf`. Path to genome-aligned, indexed bam.
/// `outdir` - `PathBuf`. Path to output directory. Must be writeable to store output bams/vcfs.
/// `sample_id`: `String`,
/// `min_base_quality`: `u8`,
/// `depth`: `Option<GenomeDepthResult>`,
/// `sample_sex`: `Sex`,
/// `homopolymer_window_size`: `Option<usize>`,
/// `big_deletion_settings`: `BigDeletionSettings`,
/// `max_number_deletions`: `u32`,
/// `site_selection_settings`: `site_selection::Settings`,

#[derive(Debug, Clone)]
pub struct Settings {
    pub genome_reference: PathBuf,
    pub genome_bam: PathBuf, // Path to sorted, indexed genome-aligned bam
    pub gene_name: String,
    pub outdir: PathBuf,
    pub sample_id: String,
    pub min_base_quality: u8,
    pub depth: Option<GenomeDepthResult>,
    pub sample_sex: Sex,
    pub homopolymer_window_size: Option<usize>,
    pub big_deletion_settings: BigDeletionSettings,
    pub site_selection_settings: SiteSelectionSettings,
    pub region_config: RegionConfig,
    pub max_number_deletions: u32,
    pub allow_low_coverage: bool,
    pub min_hap_support: i32,
    pub genome: String,
    pub min_variant_frequency: Option<f64>,
    pub min_haplotype_frequency: f64,
    pub targeted: bool,
}

impl Settings {
    ///
    /// Construct settings from sample id, gene name, and paths to aligned bam + reference + output directory.
    ///
    /// Also requires a `paraphase::config::region::Config`, from which site selection parameters are extracted.
    ///
    /// # Arguments
    /// `sample_id`: String representation for sample.
    /// `genome_data`: Tuple (`genome_ref`, `genome_bam`) for genomic reference + reads aligned to this reference.
    /// `outdir`: Directory for storing outputs + intermediate results.
    /// `gene_name`: gene name as a string. Must be a key in the `region_config` file.
    /// `depth`: genomic depth.
    /// `sex`: sample sex. `Sex::{Male, Female, Other}`.
    ///
    /// # Panics
    /// If region config does not contain the gene in question.
    pub fn new(
        sample_id: impl Into<String>,
        genome_data: (impl Into<PathBuf>, impl Into<PathBuf>),
        outdir: impl Into<PathBuf>,
        gene_name: impl Into<String>,
        region_config: &config::region::Config,
        depth: Option<GenomeDepthResult>,
        sample_sex: Option<Sex>,
        genome: impl Into<String>,
        min_variant_frequency: Option<f64>,
        min_haplotype_frequency: f64,
        targeted: bool,
    ) -> Self {
        let sample_id = sample_id.into();
        let genome = genome.into();
        let (genome_reference, genome_bam) = genome_data;
        let genome_reference = genome_reference.into();
        let genome_bam = genome_bam.into();
        let outdir = outdir.into();
        let gene_name = gene_name.into();
        let sample_sex = sample_sex.unwrap_or(Sex::Other);
        let site_selection_settings = SiteSelectionSettings::new_from_settings(
            region_config
                .get(&gene_name)
                .unwrap_or_else(|| panic!("gene name {gene_name} not in region config.")),
        );
        let region_config = region_config.clone();
        Self {
            genome_reference,
            genome_bam,
            gene_name,
            outdir,
            sample_id,
            depth,
            sample_sex,
            site_selection_settings,
            region_config,
            genome,
            min_variant_frequency,
            min_haplotype_frequency,
            targeted,
            ..Default::default()
        }
    }
}

impl std::default::Default for Settings {
    fn default() -> Self {
        Self {
            sample_sex: Sex::Other,
            min_base_quality: defaults::MIN_BASE_QUALITY,
            sample_id: String::from("Unknown"),
            genome_reference: PathBuf::default(),
            genome_bam: PathBuf::default(),
            gene_name: String::default(),
            outdir: PathBuf::default(),
            depth: None,
            homopolymer_window_size: None,
            max_number_deletions: 2,
            big_deletion_settings: BigDeletionSettings::default(),
            site_selection_settings: SiteSelectionSettings::default(),
            allow_low_coverage: false,
            min_hap_support: 4,
            region_config: config::Region::default(),
            genome: String::from("38"),
            min_variant_frequency: None,
            min_haplotype_frequency: 0.03,
            targeted: false,
        }
    }
}

/// Bit flags for Phaser, which lets us pack 4 bools into one u8.
#[repr(u8)]
pub enum FlagBits {
    ToPhase = 1,
    UseSupplementary = 2,
    IsReverse = 4,
    ExpectCN2 = 8,
}

/// Merged config for both locus and genes.
#[derive(Debug, Clone)]
pub struct MergedConfig {
    pub locus: LocusConfig,
    pub gene: GeneConfig,
}

/// Core Phaser
#[derive(Debug)]
pub struct Phaser {
    // Main settings: reference + aligned bam paths, output folder.
    pub settings: Settings,
    pub realign_settings: RealignSettings,

    // Flags
    pub flag: u8,

    // Analysis parameters for specific gene
    pub config: MergedConfig,

    pub realign_region: String,

    pub left_boundary: Option<i64>,
    pub right_boundary: Option<i64>,

    pub gene_start: Option<i64>,
    pub gene_end: Option<i64>,

    pub add_sites: Vec<CandidateSite>,
    pub clip_3p_positions: Vec<i64>,
    pub clip_5p_positions: Vec<i64>,
    pub noisy_regions: Vec<range::I64>,
    pub pivot_site: Option<i64>,

    // Runtime details
    pub low_complexity_sites: LowConfidenceSites,
    pub het_sites: Vec<CandidateSite>,
    pub init_het_sites: Vec<CandidateSite>,
    pub het_sites_no_phasing: Vec<CandidateSite>,
    pub hom_sites: Vec<CandidateSite>,
    pub candidate_sites: BTreeSet<CandidateSite>,
    pub matches: MatchMap, // match in paraphase - maps between coordinates which match between gene and pseudogene.
    pub region_avg_depth: Vec<(f32, f32)>,

    // Deletion regions: usually 0-2 // Could be parameterized later, but using defaults for now.
    pub del_data: Vec<DeletionDatum>,
}

impl Phaser {
    /// Convert `chr:pos1-pos2` format to `chr_pos1_pos2`.
    /// This way we don't have to store both forms.
    /// Warning: these coordinates are one-based.
    #[must_use]
    pub fn realign_region_old(&self) -> Option<String> {
        let (chr, start, stop) = self.parsed_nchr()?;
        Some(format!("{chr}_{start}_{stop}"))
    }

    ///
    /// Format secondary region in chr:start-stop format.
    /// # Panics
    /// Panics if the gene2 region is the wrong type.
    #[must_use]
    pub fn secondary_region_old(&self) -> Option<String> {
        self.locus_config()
            .gene2_region(self.settings.genome == "37")
            .and_then(|region| {
                let (chr, coords) = region.split_terminator(':').next_tuple()?;
                let (start, stop) = coords.split_terminator('-').next_tuple()?;
                Some(format!("{chr}_{start}_{stop}"))
            })
    }

    /// Yields a Some(tuple) of 3 strings for chr, start, stop if valid, None otherwise.
    /// Warning: these coordinates are one-based.
    #[must_use]
    pub fn parsed_nchr(&self) -> Option<(&str, &str, &str)> {
        let (chr, coords) = self.realign_region.split_terminator(':').next_tuple()?;
        let (start, stop) = coords.split_terminator('-').next_tuple()?;
        Some((chr, start, stop))
    }

    /// Yields a `Some(tuple)` of 3 strings for chr, start, stop if valid, None otherwise.
    /// Warning: these coordinates are one-based.
    ///
    /// # Panics
    /// Malformatted region: requires `chr:pos1-pos2` format.
    #[must_use]
    pub fn parsed_nchr_0based(&self) -> Option<(&str, i64, i64)> {
        let (chr, coords) = self.realign_region.split_terminator(':').next_tuple()?;
        let (start, stop) = coords
            .split_terminator('-')
            .filter_map(|x| x.parse::<i64>().ok())
            .next_tuple()?;
        Some((chr, start - 1, stop - 1))
    }
    /// Yields a `Some(tuple)` of 3 strings for chr, start, stop if valid, None otherwise.
    /// Warning: these coordinates are one-based.
    ///
    /// # Panics
    /// Malformatted region: requires `chr:pos1-pos2` format.
    #[must_use]
    pub fn parsed_nchr_secondary_0based(&self) -> Option<(&str, i64, i64)> {
        self.locus_config()
            .gene2_region(self.settings.genome == "37")
            .and_then(|region| {
                let (chr, coords) = region.split_terminator(':').next_tuple()?;
                let (start, stop) = coords
                    .split_terminator('-')
                    .filter_map(|x| x.parse::<i64>().ok())
                    .next_tuple()?;
                Some((chr, start - 1, stop - 1))
            })
    }

    ///
    /// Get genomic target name. For example, `chr5` or `5`.
    ///
    /// # Returns
    /// None if `realign_region` is malformatted.
    /// `Some(&str)` for target name otherwise.
    #[must_use]
    pub fn chr(&self) -> Option<&str> {
        self.parsed_nchr().map(|x| x.0)
    }

    #[must_use]
    pub fn genome_tid(&self) -> Option<u32> {
        let bam = self.genome_bam();
        self.chr().and_then(|chr| bam.header().tid(chr.as_bytes()))
    }

    /// Build fasta reader from the path of `Phaser::reference`.
    /// We use this function instead of storing the reader for ease of borrowing.
    ///
    /// # Errors
    /// `rust_htslib::errors::Error` - thrown if fasta file not present, index not present, or index malformed.
    pub fn make_faidx(&self) -> Result<faidx::Reader, rust_htslib::errors::Error> {
        faidx::Reader::from_path(&self.settings.genome_reference)
    }

    /// Build fasta reader from the path of `Phaser::local_reference`.
    /// We use this function instead of storing the reader for ease of borrowing.
    /// # Errors
    /// `rust_htslib::errors::Error` - thrown if fasta file for local reference is not present, index not present, or index malformed.
    pub fn make_local_faidx(&self) -> Result<faidx::Reader, Exception> {
        log::debug!("[make_local_faidx] making local reference");
        let res = faidx::Reader::from_path(
            self.local_reference()
                .map_err(|e| {
                    Exception::new(format!(
                        "Error in make_local_faidx local reference path: {e:?}. Path: {:?}",
                        self.local_reference()
                    ))
                })?
                .0,
        )
        .map_err(|e| {
            Exception::new(format!(
                "Error in make_local_faidx from_path: {e:?}. Path: {:?}",
                self.local_reference()
            ))
        });
        log::debug!("local faidx ret: {res:?}");
        res
    }

    /// Build fasta reader from the path of `Phaser::secondary_reference`.
    /// We use this function instead of storing the reader for ease of borrowing.
    /// # Errors
    /// `rust_htslib::errors::Error` - thrown if fasta file for local reference is not present, index not present, or index malformed.
    pub fn make_local_faidx_gene2(&self) -> Result<faidx::Reader, Exception> {
        log::debug!("[make_local_faidx_gene2] making local reference");
        let res = faidx::Reader::from_path(
            self.secondary_reference()
                .map_err(|e| {
                    Exception::new(format!(
                        "Error in make_local_faidx_gene2 local reference path: {e:?}. Path: {:?}",
                        self.secondary_reference()
                    ))
                })?
                .0,
        )
        .map_err(|e| {
            Exception::new(format!(
                "Error in make_local_faidx_gene2 from_path: {e:?}. Path: {:?}",
                self.secondary_reference()
            ))
        });
        log::debug!("gene2 local faidx ret: {res:?}");
        res
    }

    /// Provides the index to access in the bam for the start of the relevant region.
    /// Note the subtraction of 1, as this is the 0-based bam index built from a samtools region string, which is 1-based.
    ///
    /// # Panics
    /// 1. If region string is malformatted.
    /// 2. If start coordinate cannot be parsed as an `i64`.
    #[must_use]
    pub fn offset(&self) -> i64 {
        let (_chr, start, _stop) = self.parsed_nchr().expect("Failed to parse region string");
        start.parse::<i64>().unwrap() - 1
    }

    /// Provides the index to access in the bam for the start of the relevant region.
    /// Note the subtraction of 1, as this is the 0-based bam index built from a samtools region string, which is 1-based.
    ///
    /// # Panics
    /// 1. If region string is malformatted.
    /// 2. If start coordinate cannot be parsed as an `i64`.
    #[must_use]
    pub fn secondary_offset(&self) -> Option<i64> {
        self.parsed_nchr_secondary_0based().map(|x| x.1)
    }

    /// Right boundary for gene. Usually 1-based.
    /// # Panics
    /// 1. malformatted region string: if not in expected format or `str::parse::<i64>()` returns an error.
    #[must_use]
    pub fn right_boundary(&self) -> i64 {
        if let Some(right) = self.right_boundary {
            right
        } else {
            self.parsed_nchr()
                .expect("Failed to parse region string")
                .2
                .parse::<i64>()
                .unwrap()
        }
    }

    #[must_use]
    /// # Panics
    /// 1. malformatted region string: if not in expected format or `str::parse::<i64>()` returns an error.
    pub fn right_boundary_0based(&self) -> i64 {
        self.right_boundary() - 1
    }

    /// Left boundary for gene. Usually 1-based.
    /// # Panics
    /// 1. malformatted region string: if not in expected format or `str::parse::<i64>()` returns an error.
    #[inline]
    #[must_use]
    pub fn left_boundary(&self) -> i64 {
        self.left_boundary.unwrap_or_else(|| {
            self.parsed_nchr()
                .expect("Failed to parse region string")
                .1
                .parse::<i64>()
                .unwrap()
        })
    }

    /// # Panics
    /// 1. malformatted region string: if not in expected format or `str::parse::<i64>()` returns an error.
    #[inline]
    #[must_use]
    pub fn left_boundary_0based(&self) -> i64 {
        self.left_boundary() - 1
    }

    #[inline]
    #[must_use]
    pub fn pivot_site_0based(&self) -> Option<i64> {
        self.pivot_site.map(|x| x - 1)
    }

    /// Coordinate for `gene_start`.
    /// 1-based.
    /// # Panics
    /// 1. malformatted region string: if not in expected format or `str::parse::<i64>()` returns an error.
    #[must_use]
    pub fn gene_start(&self) -> i64 {
        self.gene_start.unwrap_or_else(|| self.left_boundary())
    }

    /// Coordinate for `gene_end`.
    /// 1-based.
    /// # Panics
    /// 1. malformatted region string: if not in expected format or `str::parse::<i64>()` returns an error.
    #[must_use]
    pub fn gene_end(&self) -> i64 {
        self.gene_end.unwrap_or_else(|| self.right_boundary())
    }

    #[must_use]
    /// Open a `bam::IndexedReader` from `self.realigned_bam_path()`.
    ///
    /// # Panics
    /// 1. bam is not present.
    /// 2. bam is malformed.
    /// 3. index file is not present.
    pub fn realigned_bam(&self) -> bam::IndexedReader {
        util::read_indexed_bam(self.realigned_bam_path().display().to_string())
            .expect("Failed to open realigned bam file.")
    }

    /// Open a `bam::IndexedReader` from `self.genome_bam_path()`.
    ///
    /// # Panics
    /// 1. bam is not present.
    /// 2. bam is malformed.
    /// 3. index file is not present.
    #[must_use]
    pub fn genome_bam(&self) -> bam::IndexedReader {
        util::read_indexed_bam(self.genome_bam_path().display().to_string())
            .expect("Failed to open genome bam file.")
    }

    /// Computes median depth for a given region using a given stride for sampling.
    ///
    /// `left_boundary` and `right_boundary` are 1-based, so using one-based is usually correct.
    ///
    /// Arguments:
    /// * bam - `&mut bam::IndexedReader`. indexed bam to read from.
    /// * query - `&[range::I64]`, a set of intervals to query.
    /// * `num_intervals`: `Option<usize>` - stride for depth sampling. At least one position is sampled. Defaults to `100` if `None`.
    /// * `exclude_flag`: `Option<u16>`. bam flags to exclude. Defaults to `0x704u16` if None.
    /// * `one_based`: `Option<bool>`. Whether given positions are one-based. Defaults to true if None.
    ///
    /// Returns: `(Vec<(f32, f32)>)`, where the first is median depth per region, and the second is 80th percentile.
    #[must_use]
    pub fn regional_depth(
        bam: &mut bam::IndexedReader,
        tid: i32,
        query: &[range::I64],
        num_intervals: Option<usize>,
        exclude_flag: Option<u16>,
        one_based: Option<bool>,
        percentile: Option<i32>,
    ) -> Vec<(f32, f32)> {
        log::trace!(
            "tid = {tid}/{} for bam targets {:?}",
            VStr::from(bam.header().target_names()[tid as usize]),
            bam.header()
                .target_names()
                .into_iter()
                .map(VString::from)
                .collect::<Vec<_>>(),
        );
        // Defaults to one-based, like Paraphase currently.
        let one_based = usize::from(one_based.unwrap_or(true));
        let exclude_flag = exclude_flag.unwrap_or(0x704u16);
        let percentile = percentile.unwrap_or(defaults::PERCENTILE);
        log::debug!(
            "Getting depth of coverage for query {query:?} and num intervals {num_intervals:?}"
        );
        let num_intervals = num_intervals.unwrap_or(defaults::NUM_SAMPLED_DEPTH_POSITIONS);
        let depths = query
            .iter()
            .map(|q| {
                let step_by = std::cmp::max(1usize, q.len() / num_intervals);
                log::debug!("Step by {step_by} for {num_intervals}");
                let sampled_positions = (q.start..q.end)
                    .step_by(std::cmp::max(1usize, q.len() / num_intervals))
                    .collect::<Vec<_>>();
                log::debug!("Sampling depth from {sampled_positions:?}");
                let depths = (q.start..q.end)
                    .step_by(std::cmp::max(1usize, q.len() / num_intervals))
                    .map(|pos| count_pos(pos - one_based as i64, tid, bam, exclude_flag))
                    .collect::<Vec<_>>();
                log::debug!("Depths: {depths:?}");
                let median = median_i32(&depths);
                let percentile = percentile_i32(&depths, percentile) as f32;
                (median, percentile)
            })
            .collect::<Vec<_>>();
        depths
    }

    /// Determines if the sample should be failed for this region for coverage reasons.
    /// # Returns
    /// `bool`: whether or not coverage passes.
    ///         Usually a true means coverage is sufficient, but if `settings.allow_low_coverage` is set, it can be overridden.
    /// # Panics
    /// 1. If genome chr is missing from bam file.
    #[must_use]
    pub fn coverage_passes(&mut self) -> bool {
        let mut bam = self.realigned_bam(); // Not genome bam - credit XC.
        let left = self.left_boundary();
        let left_bound = self
            .clip_5p_positions
            .iter()
            .max()
            .map_or(left, |x| *x.max(&left));
        let right = self.right_boundary();
        let right_bound = self
            .clip_3p_positions
            .iter()
            .min()
            .map_or(right, |x| *x.min(&right));
        // left_boundary, clip_5p_positions, clip_3p_positions are all 1-based.
        let depth1 = Self::regional_depth(
            &mut bam,
            self.genome_tid().map(|x| x as i32).expect("No genome tid"),
            &[range::I64::new(left, right)],
            /* num_intervals (step) */ None,
            /* exclude_flag */ None,
            /* one_based */ Some(true),
            /* percentile */ Some(defaults::PERCENTILE),
        );
        self.region_avg_depth = depth1.clone();
        log::debug!("left_bound {left_bound} right_bound {right_bound}");
        if left_bound < right_bound && (left_bound != left || right_bound != right) {
            let depth2 = Self::regional_depth(
                &mut bam,
                self.genome_tid().map(|x| x as i32).expect("No genome tid"),
                &[range::I64::new(left_bound, right_bound)],
                /* num_intervals (step) */ None,
                /* exclude_flag */ None,
                /* one_based */ Some(false),
                /* percentile */ Some(defaults::PERCENTILE),
            );
            log::debug!("depth1 {:?} depth2 {:?}", depth1[0], depth2[0]);
            if depth2[0].0 > depth1[0].0 {
                self.region_avg_depth = depth2;
            }
        }
        // If either is NAN, it does not return true for these conditions, so we only need to check the parameters.
        let (median, percentile) = self.region_avg_depth[0];
        log::debug!("Coverage stats - median={median}, percentile={percentile}. Required median 10, {}th percentile {}", defaults::PERCENTILE, 50.);
        (median > 8. || percentile >= 50.) || self.settings.allow_low_coverage
    }

    /// We will make this behavior polymorphic when necessary.
    /// Currently only allows deleted bases within large deletions.
    #[must_use]
    pub fn allow_del_bases(&self, pos: i64) -> bool {
        self.del_data.iter().any(|deletion| {
            /*
            log::debug!(
                "{} reads partial for fivep {}..={}",
                deletion.del_reads_partial.len(),
                deletion.fivep().start,
                deletion.threep().end
            );
            */
            !deletion.del_reads_partial.is_empty()
                && deletion.threep().start <= pos
                && pos <= deletion.fivep().end
        })
    }

    /// Compute ref seq, variant seq, and indel length for indels.
    ///
    /// # Errors
    /// 1. insertion length > `i32::MAX`.
    /// 2. No deletion found in minus indel branch. (Should never happen.)
    /// 3. faidx fetch failure.
    /// 4. Malformatted region string (accessing target name).
    pub fn process_indel(
        &self,
        pos: i64,
        ref_seq: &[u8],
        var_seq: &[u8],
        cached_faidx: &[u8],
    ) -> Result<(VString, VString, i32), DError> {
        Self::free_process_indel(pos, ref_seq, var_seq, cached_faidx)
    }

    /// Compute ref seq, variant seq, and indel length for indels.
    /// Does not require borrowing the `Phaser` struct.
    ///
    /// # Errors
    /// 1. insertion length > `i32::MAX`.
    /// 2. No deletion found in minus indel branch. (Should never happen.)
    /// 3. faidx fetch failure.
    pub fn free_process_indel(
        pos: i64,
        ref_seq: &[u8],
        var_seq: &[u8],
        cached_faidx: &[u8],
    ) -> Result<(VString, VString, i32), DError> {
        let indel_size: i32;
        let mut var_ret = VString::default();
        let mut ref_ret = VString::default();
        if let Some(plus_index) = var_seq.iter().position(|&x| x == b'+') {
            let insertion = var_seq.split_at(plus_index + 1).1;
            let insertion = skip_digits(insertion);
            indel_size = i32::try_from(insertion.len())?;
            // var_seq = ref_seq + ins_base
            var_ret.extend_from_slice(ref_seq);
            var_ret.extend_from_slice(insertion);
            // ref_ret needs to be assigned
            ref_ret.extend_from_slice(ref_seq);
        } else {
            var_ret.extend_from_slice(ref_seq);
            let minus_index = if let Some(index) = var_seq.iter().position(|&x| x == b'-') {
                Result::<usize, simple_error::SimpleError>::Ok(index)
            } else {
                simple_error::bail!(
                    "- and + not found in indel seqs. var seq: {}",
                    vstr::VStr::from(var_seq)
                )
            }?;
            let deletion_len = extract_del_length(var_seq.split_at(minus_index + 1).1);
            indel_size = deletion_len;
            let deletion_len = deletion_len as usize;
            //let offset_pos = (pos - offset) as usize;
            let pos = pos as usize;
            // pysam's fetch subtracts 1 from the end.
            // paraphase adds one to the query.
            // We do both and let compile-time constant folding take care of this.

            // Subtract 1 from the coordinate to get the ref base at this position.
            let cached_seq = &cached_faidx[pos..pos + deletion_len + 1 as usize]; // Get deletion length + first base.
            debug_assert_eq!(cached_seq, cached_seq.to_ascii_uppercase());
            debug_assert_eq!(cached_seq.len(), deletion_len as usize + 1);
            ref_ret.extend_from_slice(cached_seq);
        }
        ref_ret.make_ascii_uppercase();
        var_ret.make_ascii_uppercase();
        Ok((ref_ret, var_ret, indel_size))
    }
    /// Get records overlapping pos1..pos2 and store as a Vec<bam::Record>.
    /// Returns a Result which only fails if the fetch does.
    fn get_records_for_range(&mut self, range: &range::I64) -> Result<Vec<bam::Record>, DError> {
        self.get_records(range.start, range.end)
    }

    /// Get records overlapping pos1..pos2 and store as a Vec<bam::Record>.
    /// Returns a Result which only fails if the fetch does.
    fn get_records(&mut self, pos1: i64, pos2: i64) -> Result<Vec<bam::Record>, DError> {
        let mut realigned_bam = self.realigned_bam();
        let tid = self.genome_tid().map(|x| x as i32).expect("No chr tid");
        realigned_bam.fetch((tid, pos1, pos2))?;
        Ok(realigned_bam
            .records()
            .filter_map(Result::ok)
            .collect::<Vec<_>>())
    }

    /// Find big deletions in data, use them for variant identification.
    /// # Panics
    /// If sanity check fails. (If read with deletion does not have deletion cigar operation, something went wrong.)
    fn discover_big_dels(&mut self) -> DResult {
        let mut realigned_bam = self.realigned_bam();
        let tid = self.genome_tid().map(|x| x as i32).expect("No chr tid");
        realigned_bam.fetch((
            tid,
            self.left_boundary_0based(),
            self.right_boundary_0based(),
        ))?;
        let del_reads = realigned_bam
            .rc_records()
            .filter_map(Result::ok)
            .map(|x| {
                let max_len = i64::from(
                    x.cigar()
                        .iter()
                        .map(|x| match x {
                            bam::record::Cigar::Del(x) => *x,
                            _ => 0,
                        })
                        .max()
                        .unwrap_or(0),
                );
                log::trace!(
                    "[discover_big_dels] raw deletion reads candidates {} deletion length {}",
                    std::str::from_utf8(&x.qname()).unwrap(),
                    max_len
                );
                (x, max_len)
            })
            .filter(|(_x, max_len)| *max_len >= self.settings.big_deletion_settings.min_size)
            .map(|(x, max_len)| {
                let del_position = x
                    .cigar()
                    .iter()
                    .position(|&x| x == bam::record::Cigar::Del(max_len as u32))
                    .expect("Sanity check");
                let del_pos = x.pos() - 1
                    + x.cigar()
                        .iter()
                        .take(del_position)
                        .map(|x| {
                            if let bam::record::Cigar::Match(_len) = x {
                                i64::from(x.len())
                            } else if let bam::record::Cigar::Del(_len) = x {
                                i64::from(x.len())
                            } else if let bam::record::Cigar::Equal(_len) = x {
                                i64::from(x.len())
                            } else if let bam::record::Cigar::Diff(_len) = x {
                                i64::from(x.len())
                            } else {
                                0
                            }
                        })
                        .sum::<i64>();
                log::trace!(
                    "[discover_big_dels] reads with long deletions {} start {} end {}",
                    std::str::from_utf8(&x.qname()).unwrap(),
                    del_pos,
                    del_pos + max_len
                );
                (del_pos, del_pos + max_len)
            })
            .collect::<Vec<(i64, i64)>>();
        log::debug!("[discover_big_dels] del_reads {del_reads:?}");
        let common = del_reads
            .iter()
            .copied()
            .collect::<counter::Counter<(i64, i64), i64>>()
            .k_most_common_ordered(self.settings.max_number_deletions as usize);

        let padding = self.settings.big_deletion_settings.padding;
        self.del_data = common
            .into_iter()
            .filter_map(|x| {
                if x.1 >= self.settings.big_deletion_settings.min_count {
                    Some(x.0)
                } else {
                    None
                }
            })
            .map(|item| {
                let raw = item.into();
                DeletionDatum::new(raw, Some(padding), None, None, None, None)
            })
            .collect::<Vec<_>>();
        log::debug!("Deletion data: {:?}", self.del_data);
        Ok(())
    }

    /// Call recurrent clip sites from data
    pub fn find_clip_site(
        &mut self,
        min_clip_length: Option<u32>,
        min_count: Option<usize>,
        padding: Option<i64>,
    ) -> DResult {
        let min_clip_length = min_clip_length.unwrap_or(800);
        let min_count = min_count.unwrap_or(6);
        let padding = padding.unwrap_or(1000);
        let mut clip_reads = Vec::new();
        let mut realigned_bam = self.realigned_bam();
        let tid = self.genome_tid().map(|x| x as i32).expect("No chr tid");
        realigned_bam.fetch((
            tid,
            self.left_boundary_0based(),
            self.right_boundary_0based(),
        ))?;
        let records = realigned_bam
            .records()
            .filter_map(Result::ok)
            .collect::<Vec<_>>();
        for record in records {
            let cigar = record.cigar();
            let clip_len_5p = fiveprime_clip_length(&cigar);
            if clip_len_5p >= min_clip_length {
                let clip_pos = record.reference_start();
                if clip_pos > self.left_boundary_0based() + padding
                    && clip_pos < self.right_boundary_0based() - padding
                {
                    clip_reads.push((clip_pos, "5p"));
                }
            }
            let clip_len_3p = threeprime_clip_length(&cigar);
            if clip_len_3p >= min_clip_length {
                let clip_pos = record.reference_end();
                if clip_pos > self.left_boundary_0based() + padding
                    && clip_pos < self.right_boundary_0based() - padding
                {
                    clip_reads.push((clip_pos, "3p"));
                }
            }
        }
        let clip_counter = clip_reads.iter().copied().collect::<counter::Counter<_>>();
        for ((pos, clip_direction), this_count) in clip_counter {
            if this_count >= min_count {
                if clip_direction == "5p" {
                    // Check if position is too close to existing clip positions
                    let too_close = self
                        .clip_5p_positions
                        .iter()
                        .any(|&existing_pos| (pos - 100..pos + 100).contains(&existing_pos));
                    if !too_close {
                        let mut in_deletion = false;
                        for each_known_deletion in &self.del_data {
                            if ranges_overlap(
                                pos - 100,
                                pos + 100,
                                each_known_deletion.fivep().start,
                                each_known_deletion.fivep().end,
                            ) {
                                in_deletion = true;
                            }
                        }
                        if !in_deletion {
                            self.clip_5p_positions.push(pos);
                        }
                    }
                }
                if clip_direction == "3p" {
                    // Check if position is too close to existing clip positions
                    let too_close = self
                        .clip_3p_positions
                        .iter()
                        .any(|&existing_pos| (pos - 100..pos + 100).contains(&existing_pos));
                    if !too_close {
                        let mut in_deletion = false;
                        for each_known_deletion in &self.del_data {
                            if ranges_overlap(
                                pos - 100,
                                pos + 100,
                                each_known_deletion.threep().start,
                                each_known_deletion.threep().end,
                            ) {
                                in_deletion = true;
                            }
                        }
                        if !in_deletion {
                            self.clip_3p_positions.push(pos);
                        }
                    }
                }
            }
        }
        // Sort positions
        self.clip_5p_positions.sort();
        self.clip_3p_positions.sort();

        Ok(())
    }

    ///
    /// Updates `self.del_data` vector with relevant read names per deletion.
    /// Because we use a vector instead of fixed members for del1, del2, we just do a loop.
    /// Also, we return nothing besides a `Result<(), dyn std::error::Error>` because we update the member.
    ///
    /// # Errors
    /// Failure to get records for deletion sites.
    ///
    /// # Future improvements
    /// This is the function in which we would use graph wfa for variant adjudication.
    /// After using `discover_big_dels` and `get_candidate_pos` to find candidate variants, we would insert them into this graph before realigning and reassigning.
    fn label_big_dels(&mut self) -> DResult {
        // Copy out coords to avoid conflicting mutable + immutable borrows.
        log::debug!("Label big deletions");
        //log::warn!("We fetch reads inefficiently. We can store them in memory for the whole processs and filter at some point for more speed.");
        // Get relevant reads.
        let mut fetched_reads = Vec::with_capacity(self.del_data.len());
        for (threep, fivep) in self.del_data.clone().iter().map(|x| {
            let threep = self.get_records_for_range(&x.threep());
            let fivep = self.get_records_for_range(&x.fivep());
            (threep, fivep)
        }) {
            fetched_reads.push((threep?, fivep?));
        }
        let min_clip_len = self.settings.big_deletion_settings.min_clip_len;
        let min_extend = self.settings.big_deletion_settings.min_extend;
        let padding_negative_reads = self.settings.big_deletion_settings.padding_negative_reads;
        let use_supplementary = self.use_supplementary();
        for ((threep, fivep), del_data) in fetched_reads.iter_mut().zip(self.del_data.iter_mut()) {
            let mut p3_reads = BTreeSet::<String>::new();
            let mut p5_reads = BTreeSet::<String>::new();
            let mut del_reads = BTreeSet::<String>::new();
            // Handle 3'
            for read in threep {
                let read_name = Self::get_read_name_free(read, use_supplementary);
                read.cache_cigar();
                if read.reference_start() < del_data.threep().start - padding_negative_reads
                    && read.reference_end() > del_data.threep().end + padding_negative_reads
                {
                    del_data.del_negative_reads.insert(read_name.clone());
                }

                let reference_start_cutoff = del_data.threep().start - min_extend;
                let threep_length = threeprime_clip_length(&read.cigar());
                let end = read.reference_end();

                if (i64::from(threep_length) >= min_clip_len)
                    && (del_data.threep().start < end)
                    && (end < del_data.threep().end)
                    && read.pos() < reference_start_cutoff
                {
                    p3_reads.insert(read_name.clone());
                }
                let has_deletion_in_cigar =
                    crate::detail::phaser_util::check_del(read, del_data.size(), del_data.threep());
                if has_deletion_in_cigar {
                    del_reads.insert(read_name.clone());
                }
            }

            // Handle 5'
            for read in fivep {
                let read_name = Self::get_read_name_free(read, use_supplementary);
                read.cache_cigar();
                if read.reference_start() < del_data.fivep().start - padding_negative_reads
                    && read.reference_end() > del_data.fivep().end + padding_negative_reads
                {
                    del_data.del_negative_reads.insert(read_name.clone());
                }

                let reference_end_cutoff = del_data.fivep().end + min_extend;
                let fivep_length = fiveprime_clip_length(&read.cigar());
                let pos = read.pos();

                if (i64::from(fivep_length) >= min_clip_len)
                    && (del_data.fivep().start < pos)
                    && (pos < del_data.fivep().end)
                    && read.reference_end() > reference_end_cutoff
                {
                    p5_reads.insert(read_name.clone());
                }
            }
            let del_reads = del_reads;
            log::debug!("del: {del_reads:?}. p3: {p3_reads:?}. p5: {p5_reads:?}");
            let p3_p5_del_union;
            let full_union;
            if !del_reads.is_empty() || (!p3_reads.is_empty() && !p5_reads.is_empty()) {
                let p3_p5_intersect = p3_reads
                    .intersection(&p5_reads)
                    .cloned()
                    .collect::<BTreeSet<_>>();
                p3_p5_del_union = p3_p5_intersect
                    .union(&del_reads)
                    .cloned()
                    .collect::<BTreeSet<_>>();
                full_union = del_reads
                    .union(&p3_reads)
                    .chain(p5_reads.iter())
                    .cloned()
                    .collect::<BTreeSet<_>>();
            } else {
                p3_p5_del_union = BTreeSet::new();
                full_union = BTreeSet::new();
            }
            del_data.del_reads = p3_p5_del_union;
            del_data.del_reads_partial = full_union;
            log::debug!("DeletionDatum: {del_data:?} after labeling");
        }
        Ok(())
    }

    ///
    /// # Errors
    /// 1. Failure to get path to local reference.
    /// 2. Failure to align.
    /// 3. Failure to postprocess bam.
    /// 4. Failure to remove temporary file.
    /// 5. Failure to build bam index.
    ///
    /// # Panics
    /// Panics if `self.chr()` is `None`.
    pub fn align(&mut self) -> DResult {
        let local_bam = self.realigned_bam_path();
        log::debug!("[align] making local reference");
        let (local_ref, _was_new) = self.local_reference()?;
        let regions_to_extract = self
            .locus_config()
            .extract_regions(self.settings.genome == "37");
        let regions_to_extract_view = regions_to_extract
            .iter()
            .map(|x| &x[..])
            .collect::<Vec<_>>();
        log::debug!("Aligning {local_bam:?} to ref {local_ref:?}");

        let opts = (
            1,
            self.locus_config().chain_bandwidth(),
            self.realign_settings,
            self.offset(),
        );
        align_mm2_intrinsic(
            &self.genome_bam_path(),
            &local_bam,
            &local_ref,
            &regions_to_extract_view,
            opts,
        )?;

        bam::index::build(&local_bam, None, bam::index::Type::Bai, 1)?;
        if self
            .locus_config()
            .gene2_region(self.settings.genome == "37")
            .is_some()
        {
            let (secondary_ref, _was_new) = self.secondary_reference()?;
            let opts = (
                1,
                self.locus_config().chain_bandwidth(),
                self.realign_settings,
                self.secondary_offset()
                    .expect("Required: secondary offset for secondary reference region alignment"),
            );
            let secondary_bam = self.realigned_gene2_bam_path();
            align_mm2_intrinsic(
                &self.genome_bam_path(),
                &secondary_bam,
                &secondary_ref,
                &regions_to_extract_view,
                opts,
            )?;
            bam::index::build(&secondary_bam, None, bam::index::Type::Bai, 1)?;
        }
        Ok(())
    }

    /// Identify which haplotypes may have two copies based on depth.
    pub(crate) fn compare_depth(
        &self,
        haps: &BTreeMap<String, HapInfo>,
        assembled_haps: &BTreeMap<VStr<'_>, String>,
        loose: bool,
        stringent: bool,
    ) -> Result<Vec<String>, DError> {
        if haps.len() <= 1 {
            return Ok(vec![]);
        }
        let mut two_cp_haps = Vec::<String>::new();

        let mut hap_name_to_seq = BTreeMap::new();
        for (hap_seq, hap_name) in assembled_haps {
            hap_name_to_seq.insert(hap_name, hap_seq);
        }
        let mut bound_start = Vec::new();
        let mut bound_end = Vec::new();
        for (hap, hap_info) in haps {
            let bound = &hap_info.boundary;
            let hap_seq = hap_name_to_seq
                .get(hap)
                .ok_or("hap not in hap_name_to_seq")?;
            let hap_clip_5p = self.get_5pclip_from_hap(*hap_seq)?;
            let hap_clip_3p = self.get_3pclip_from_hap(*hap_seq)?;
            let mut true_start = bound.start;
            if let Some(hap_clip_5p_value) = hap_clip_5p {
                if hap_clip_5p_value != 0 {
                    true_start = hap_clip_5p_value;
                }
            }
            bound_start.push(true_start);
            let mut true_end = bound.end;
            if let Some(hap_clip_3p_value) = hap_clip_3p {
                if hap_clip_3p_value != 0 {
                    true_end = hap_clip_3p_value;
                }
            }
            bound_end.push(true_end);
        }
        let start = bound_start.iter().max().unwrap();
        let end = bound_end.iter().min().unwrap();
        let range = range::I64::new(*start, *end);
        let vars = haps
            .values()
            .flat_map(|info| {
                // TODO: check for a variant type without 3 components; possibly ins/del
                info.variants.iter().filter(|x| {
                    range.contains(&x.pos) && self.het_sites.contains(x) && x.is_unit_length()
                })
            })
            .collect::<BTreeSet<_>>();
        log::trace!("[compare_depth] vars {:?}", vars);
        let faidx = self.make_faidx()?;
        let (chrom, start, stop) = self
            .parsed_nchr_0based()
            .expect("Malformatted region string");
        let ref_seq = faidx.fetch_seq(chrom, start as usize, stop as usize)?;
        let offset = self.offset() as usize;
        let threshold: f64 = match (stringent, loose) {
            (true, _) => 0.8,
            (_, true) => 0.5,
            _ => 0.6,
        };
        let mut reader = self.realigned_bam();
        for (hap_name, info) in haps {
            let mut sites = BTreeMap::new();
            let other_haps = haps.iter().filter(|x| x.0 != hap_name).collect::<Vec<_>>();
            let other_cn = other_haps.len();
            let this_hap_var = &info.variants[..];
            let other_haps_var = other_haps
                .iter()
                .flat_map(|(_name, info)| info.variants.iter())
                .collect::<Vec<_>>();
            for var in &vars {
                let in_this = this_hap_var.contains(var);
                let in_other = other_haps_var.contains(var);
                if in_this && !in_other {
                    sites.insert(var.pos, &var.var_seq);
                } else if !in_this
                    && other_haps_var.iter().filter(|x| *x == var).count() == other_cn
                {
                    sites.insert(var.pos, &var.ref_seq);
                }
            }
            // avoid making this vector, perform prob test during sites loop.
            // let mut base_counts = Vec::<(i32, i32)>::new(); // base count, non-base count
            log::trace!("[compare_depth] hap_name {hap_name} sites {:?}", sites);
            let mut double_prob_count = 0usize;
            let tid = self.genome_tid().map(|x| x as i32).expect("No chr tid");
            for (site, base) in &sites {
                reader
                    .fetch((tid, *site, site + 1))
                    .expect("Failed to fetch");
                for pileup in reader.pileup() {
                    let pileup = pileup.expect("Failed pileup");
                    let pos = i64::from(pileup.pos());
                    match pos.cmp(site) {
                        Ordering::Less => {
                            continue;
                        }
                        Ordering::Greater => {
                            break;
                        }
                        Ordering::Equal => {}
                    }
                    assert_eq!(self.settings.min_base_quality, 25, "base quality here should be 25. In the future we can modify this as a runtime parameter.");
                    // TODO: consider just counting, no dictionary for speed.
                    let got_counts =
                        counts(&pileup, &ref_seq, self.settings.min_base_quality, offset);
                    log::trace!("got_counts: {got_counts:?}");
                    debug_assert!(
                        got_counts.keys().all(|x| **x == x.to_ascii_uppercase()),
                        "counts had non uppercase keys: {got_counts:?}"
                    );
                    let mut base_count = 0;
                    let mut total = 0;
                    for (each_base, count) in &got_counts {
                        //if *each_base == b'*' {
                        //    continue;
                        //}
                        total += count;
                        let each_base_first = *each_base.first().unwrap();
                        if **base == each_base_first {
                            base_count += count;
                        }
                    }
                    let other_count = total - base_count;
                    let double_prob =
                        depth_prob(base_count, f64::from(other_count) / other_cn as f64);
                    log::trace!(
                        "[compare_depth] hap_name {hap_name} site {site} total {total} base_count {base_count} double_prob {:?}", double_prob
                    );
                    if let Some(double_prob_value) = double_prob {
                        if double_prob_value[0] < 0.25 {
                            double_prob_count += 1;
                        }
                    }
                }
            }
            log::trace!(
                "[compare_depth] hap_name {hap_name} nsite {} double_prob_count {}",
                sites.len(),
                double_prob_count
            );
            if double_prob_count > 0
                && sites.len() >= 5
                && double_prob_count as f64 >= sites.len() as f64 * threshold
            {
                two_cp_haps.push(hap_name.clone());
            }
        }

        if two_cp_haps.len() > 1 {
            // There can be only one.
            two_cp_haps.clear();
        }

        log::debug!("Found {} two cp haps. {two_cp_haps:?}", two_cp_haps.len());
        Ok(two_cp_haps)
    }

    pub fn compare_depth_by_read_count(
        &mut self,
        assembled_haps: &BTreeMap<VStr<'_>, String>,
        phase_results: &PhasedResult,
        prob_threshold: f32,
    ) -> Vec<String> {
        let min_read_count = if self.settings.targeted { 15 } else { 10 };
        let mut two_cp = Vec::new();
        if assembled_haps.len() < 2 || phase_results.read_counts.len() < 2 {
            return two_cp;
        }
        let read_counts = &phase_results.read_counts;
        let (cp2_hap, max_count) = read_counts.iter().max_by_key(|x| x.1).unwrap();
        let others_max = read_counts
            .iter()
            .sorted_by_key(|x| std::cmp::Reverse(x.1))
            .nth(1)
            .unwrap();
        let probs = depth_prob(*max_count, f64::from(*others_max.1));
        log::debug!(
            "read_counts {:?} cp2_hap {:?} max_count {:?} probs {:?}",
            read_counts,
            cp2_hap,
            max_count,
            probs
        );
        if let Some(probs_value) = probs {
            if probs_value[0] < prob_threshold && *others_max.1 >= min_read_count {
                two_cp.push(assembled_haps[&cp2_hap.vstr()].clone());
            }
        }
        two_cp
    }

    pub fn realign(&mut self) -> Result<Vec<u8>, DError> {
        log::debug!("Getting local reference sequence.");
        let local_chr = self.local_chr().ok_or(Exception::new(format!(
            "Missing target name for region {}",
            self.realign_region
        )))?;
        let faidx = self.make_local_faidx()?;
        log::debug!(
            "Made local faidx from file at {:?}. Querying with {local_chr}:0-{}",
            self.local_reference()?,
            i32::MAX
        );
        // Copy out ref seq
        let seq = faidx
            .fetch_seq(&local_chr, 0, i32::MAX as usize)
            .map_err(|e| {
                Exception::new(format!(
                    "Failed to query region {local_chr} for faidx at {:?}. Error: {e:?}",
                    self.local_reference_path()
                ))
            })?
            .to_ascii_uppercase();
        log::debug!("Got local reference sequence.");
        log::debug!(
            "Reference seq of size {}. Sampled first 100bp: {}",
            seq.len(),
            VStr::from(&seq[..std::cmp::min(seq.len(), 100)])
        );

        log::debug!("Aligning reads to local reference");
        self.align()?;
        assert!(
            self.realigned_bam_path().exists(),
            "{:?} was not created as expected.",
            self.realigned_bam_path()
        );
        Ok(seq)
    }

    pub fn get_sites(
        &mut self,
        seq: &Vec<u8>,
        min_no_var_region_size: Option<i64>,
        min_vaf: Option<f64>,
    ) -> Result<(Vec<CandidateSite>, Vec<CandidateSite>), DError> {
        log::debug!("Coverage passes. Finding homopolymer regions");
        // 1. Get homopolymers.
        self.low_complexity_sites =
            LowConfidenceSites::new(seq, self.offset(), self.settings.homopolymer_window_size);

        // get predefined deletion
        self.parse_deletions_from_config()?;

        // 2. Discover + label big deletions
        log::debug!("Discover big deletions");
        if self.del_data.is_empty() {
            self.discover_big_dels()?;
        }
        self.label_big_dels()?;

        // 3. Handle regions to check
        let regions_to_check = self
            .del_data
            .iter()
            .flat_map(|x| {
                let mut ret = vec![];
                if !x.del_reads_partial.is_empty() {
                    ret.push(x.threep());
                    ret.push(x.fivep());
                }
                ret.into_iter()
            })
            .collect::<Vec<range::I64>>();
        log::debug!("Regions to check for deletions: {regions_to_check:?}");

        // 4. Get candidate sites + assign het sites
        let (_filtered_variants, _raw_variant_counts) =
            self.get_candidate_pos(&regions_to_check, seq, min_vaf)?;
        //call.raw_variant_counts = site_selection::raw_variants_to_string(&raw_variant_counts); // Consider removing, it is bulky.
        //call.init_filtered_sites = filtered_variants.to_json_friendly();

        // 5. Remove noisy sites + add homozygous
        self.remove_noisy_sites();
        self.init_het_sites = self.het_sites.clone();
        let hom_sites_to_add = self.add_hom_sites(min_no_var_region_size, None, seq); // Could be parameterized later, but using defaults for now.
        self.remove_noisy_sites(); // outside of add_hom_sites for ownership management.

        log::debug!("[get_sites] hom_sites_to_add: {hom_sites_to_add:?}");
        // 6. Get haps from reads

        let mut add_sites = self.add_sites.clone();
        // add pivot site
        if let Some(pivot_site) = self.pivot_site_0based() {
            let het_sites_all_pos = self.het_sites.iter().map(|x| x.pos).collect::<Vec<_>>();
            let add_sites_all_pos = add_sites.iter().map(|x| x.pos).collect::<Vec<_>>();
            if !het_sites_all_pos.contains(&pivot_site) && !add_sites_all_pos.contains(&pivot_site)
            {
                let pos_on_ref = pivot_site - self.offset();
                let ref_base_u8 = seq[pos_on_ref as usize];
                let non_ref_bases = vec![b'A', b'C', b'G', b'T']
                    .iter()
                    .filter(|x| **x != ref_base_u8)
                    .map(|x| *x)
                    .collect::<Vec<_>>();
                let var_base_u8 = non_ref_bases.first().unwrap();
                let ref_base = std::str::from_utf8(&[ref_base_u8]).unwrap().to_string();
                let var_base = std::str::from_utf8(&[*var_base_u8]).unwrap().to_string();
                let new_variant = CandidateSite::new(pivot_site, ref_base, var_base);
                add_sites.push(new_variant);
            }
        }

        log::debug!(
            "[get_sites] Haps from reads. Adding {} sites from config.",
            add_sites.len()
        );

        log::debug!(
            "[get_sites] Assigning fingerprints using site set {:?} ({})",
            self.het_sites,
            intersperse(
                self.het_sites.iter().map(ToString::to_string),
                String::from(",")
            )
            .collect::<String>()
        );
        Ok((hom_sites_to_add, add_sites))
    }

    pub fn update_indel_and_phase(
        &mut self,
        init_read_hap_map: BTreeMap<crate::io::json::ReadAlignmentId, VString>,
        call: &mut GeneCall,
    ) -> Result<(PhasedResult, BTreeMap<char, String>), DError> {
        log::debug!("Read hap map (init): {init_read_hap_map:?}");

        let mut read_hap_map = init_read_hap_map;
        // 7. Handle known deletions
        let known_del = self.update_for_deletions(&mut read_hap_map)?;
        call.read_details = to_string_map(&read_hap_map);

        // 8. Phase haps
        let phase_results = self.phase_haps(&read_hap_map)?;
        Ok((phase_results, known_del))
    }

    pub fn adjust_depth(
        &mut self,
        assembled_haps: BTreeMap<VStr<'_>, String>,
        haps: BTreeMap<String, HapInfo>,
        phase_results: PhasedResult,
        expect_cn4: bool,
        stringent: bool,
        loose: bool,
        prob_threshold: f32,
    ) -> Result<(Vec<String>, usize), DError> {
        // 10. Check haplotype depths.
        let mut two_cp_haps = Vec::<String>::new();
        if assembled_haps.len() == 1 && self.init_het_sites.is_empty() {
            two_cp_haps = assembled_haps.values().cloned().collect::<Vec<_>>();
        } else {
            if self.settings.targeted {
                two_cp_haps = self.compare_depth_by_read_count(
                    &assembled_haps,
                    &phase_results,
                    prob_threshold,
                );
            } else if (assembled_haps.len() == 3
                && !self.expect_cn2()
                && self.gene_name() != "BPY2")
                || (self.gene_name() == "BPY2" && assembled_haps.len() < 3)
            {
                two_cp_haps = self.compare_depth(&haps, &assembled_haps, loose, stringent)?;
                if two_cp_haps.is_empty() && !phase_results.read_counts.is_empty() {
                    two_cp_haps = self.compare_depth_by_read_count(
                        &assembled_haps,
                        &phase_results,
                        prob_threshold,
                    );
                }
            }
        }

        // check gene1 haplotypes and update to cn2 if assume gene1 is never cn1
        // only for targeted mode
        if self.settings.targeted && two_cp_haps.is_empty() {
            let gene1_cn2 = self.locus_config().get("gene1_cn2");
            if !gene1_cn2.is_none() {
                let region_length = self.right_boundary() - self.left_boundary();
                let snp_count = self.het_sites.len();
                log::debug!("region_length {region_length} snp_count {snp_count}");
                if snp_count as f64 > region_length as f64 * 0.008 {
                    let mut gene1_haps = Vec::new();
                    let mut gene2_haps = Vec::new();
                    for (hap_seq, hap_name) in &assembled_haps {
                        // this is assuming gene2 is very different from gene1
                        let count2 = hap_seq.iter().filter(|x| **x == b'2').count();
                        if count2 as f64 > hap_seq.len() as f64 * 0.7 {
                            gene2_haps.push(hap_name);
                        } else {
                            gene1_haps.push(hap_name);
                        }
                    }
                    log::debug!("gene1_haps {gene1_haps:?} gene2_haps {gene2_haps:?}");
                    if gene1_haps.len() == 1 {
                        two_cp_haps.push(gene1_haps.first().unwrap().to_string())
                    }
                }
            }
        }

        let mut total_cn = assembled_haps.len() + two_cp_haps.len();
        // Fully homozygous
        if assembled_haps.is_empty() && self.init_het_sites.is_empty() {
            total_cn = 2;
        }
        // two identical haplotypes
        if total_cn == 2 && !self.expect_cn2() && self.gene_name() != "BPY2" {
            if expect_cn4 {
                two_cp_haps = assembled_haps.values().cloned().collect::<Vec<_>>();
                total_cn = 4;
            } else if let Some(depth) = self.settings.depth.as_ref() {
                if self.gene_config().uses_depth(self.gene_name())
                    && depth.status() == depth::MedianDepthStatus::Passing
                {
                    let copy_number_probs =
                        depth_prob(self.region_avg_depth[0].0 as i32, depth.median as f32);
                    if let Some(copy_number_probs_value) = copy_number_probs {
                        let same_cn_prob = copy_number_probs_value[0];
                        if same_cn_prob < 0.75 {
                            total_cn = 4;
                            if two_cp_haps.is_empty() && !assembled_haps.is_empty() {
                                two_cp_haps = assembled_haps.values().cloned().collect::<Vec<_>>();
                            }
                        }
                    }
                }
            }
        }
        Ok((two_cp_haps, total_cn))
    }

    /// Fill in the fields for the final report
    pub fn fill_in_call(&mut self, phase_results: PhasedResult, call: &mut GeneCall) {
        let phase_result = PhasedResultForJson::new(&phase_results);
        call.sites_for_phasing = self
            .het_sites
            .iter()
            .map(ToString::to_string)
            .collect::<Vec<_>>();
        self.init_het_sites.sort_by(|a, b| a.pos.cmp(&b.pos));
        call.heterozygous_sites = self
            .init_het_sites
            .iter()
            .map(ToString::to_string)
            .collect::<Vec<_>>();
        call.homozygous_sites = self
            .hom_sites
            .iter()
            .map(ToString::to_string)
            .collect::<Vec<_>>();
        self.het_sites_no_phasing.sort_by(|a, b| a.pos.cmp(&b.pos));
        call.het_sites_not_used_in_phasing = self
            .het_sites_no_phasing
            .iter()
            .map(ToString::to_string)
            .collect::<Vec<_>>();
        call.assembled_haplotypes = phase_result
            .assemblies
            .final_haps
            .iter()
            .map(std::string::ToString::to_string)
            .collect::<Vec<_>>();
        call.highest_total_cn = Some(phase_result.assemblies.highest_cn as i32);
        call.unique_supporting_reads = phase_result
            .uniquely_supporting_reads
            .iter()
            .map(|(k, v)| (k.to_string(), v.clone()))
            .collect::<BTreeMap<_, _>>();
        call.nonunique_supporting_reads = phase_result
            .nonuniquely_supporting_reads
            .iter()
            .map(|(k, v)| {
                (
                    k.clone(),
                    v.iter()
                        .map(std::borrow::ToOwned::to_owned)
                        .collect::<Vec<_>>(),
                )
            })
            .collect::<BTreeMap<_, _>>();
    }

    /// Go through reads and get bases at sites of interest.
    /// Two rounds, with variant site filtering in between.
    ///
    /// # Errors
    /// 1. Error in `haplotypes_from_reads_step`
    /// 2. Failed to get fasta index for local reference.
    /// 3. Failure to fetch sequence from reference, usualy invalid ranges.
    pub fn haplotypes_from_reads(
        &mut self,
        exclude_reads: Option<&BTreeSet<String>>,
        kept_sites: &[CandidateSite],
        add_sites: Option<&[CandidateSite]>,
        partial_deletion_reads: Option<&BTreeSet<String>>,
        options: (u8, bool, Option<u32>),
        tid: i32,
        clip_buffer: Option<i32>,
        hom_sites: &[CandidateSite],
    ) -> Result<ReadFingerprintMap, DError> {
        const CLIP_OFFSET: i64 = defaults::CLIP_OFFSET_HAPS_FROM_READS;
        const ACGT: &[u8] = b"ACGT";

        let (min_mapq, check_clip, min_clip_len) = options;
        let add_sites = add_sites.unwrap_or(&self.add_sites[..]).to_vec();
        let min_clip_len = min_clip_len.unwrap_or(50);
        log::debug!("[haplotypes_from_reads] with min mapq {min_mapq}, check clip {check_clip} and min clip len {min_clip_len}");
        let mut raw_read_haps = self.haplotypes_from_reads_step(
            exclude_reads,
            min_clip_len,
            check_clip,
            partial_deletion_reads,
            min_mapq,
            tid,
            clip_buffer,
            None,
        )?;
        log::debug!(
            "[haplotypes_from_reads] {} initially labeled reads from first haplotypes_from_reads_step",
            raw_read_haps.len()
        );
        log::debug!(
            "Removing variants from hap map of fingerprint length {} and kept_sites of length {}",
            raw_read_haps.values().next().map_or(0, |x| x.len()),
            kept_sites.len()
        );
        let absent_base_per_site = self.remove_var(&raw_read_haps, kept_sites, hom_sites);

        if !self.het_sites.is_empty() {
            self.het_sites.extend_from_slice(&add_sites);
            // Eliminate dups
            self.het_sites.sort();
            self.het_sites.dedup();
        }

        let faidx = self.make_local_faidx()?; // TODO: use in-memory sequence instead of faidx
        let chr = self.local_chr().ok_or(Exception::new("chr must be set"))?;
        log::debug!("chr (local): {chr}");
        let seq = faidx
            .fetch_seq(&chr, 0, i64::MAX as usize)?
            .to_ascii_uppercase();
        debug_assert!((0..seq.len())
            .all(|idx| seq[idx]
                == faidx.fetch_seq(&chr, idx, idx + 1).unwrap()[0].to_ascii_uppercase()));
        drop(faidx);

        // Handle clips
        // add variants before 5' clip sites or after 3' clip sites
        // Handle 5'
        let num_clip = self.clip_5p_positions.len();
        for i in 0..num_clip {
            let clip_pos = self.clip_5p_positions[i];
            let has_var_before_clip = if i == 0 {
                self.het_sites.iter().any(|x| x.pos < clip_pos)
            } else {
                self.het_sites
                    .iter()
                    .any(|x| x.pos > self.clip_5p_positions[i - 1] && x.pos < clip_pos)
            };
            if !has_var_before_clip {
                let var_pos = clip_pos - CLIP_OFFSET;
                let ref_base = seq[(var_pos - self.offset()) as usize];
                let var_base =
                    ACGT.iter()
                        .copied()
                        .find(|&x| x != ref_base)
                        .ok_or(Exception::new(format!(
                    "Expected a variant base at 5' {clip_pos} position with var_pos = {var_pos}"
                )))?;
                let new_var = CandidateSite::new(var_pos, vec![ref_base], vec![var_base]);
                self.het_sites.push(new_var);
            }
        }

        // Handle 3'
        let num_clip = self.clip_3p_positions.len();
        for i in (0..num_clip).rev() {
            let clip_pos = self.clip_3p_positions[i];
            let has_var_after_clip = if i == num_clip - 1 {
                self.het_sites.iter().any(|x| x.pos > clip_pos)
            } else {
                self.het_sites
                    .iter()
                    .any(|x| x.pos < self.clip_3p_positions[i + 1] && x.pos > clip_pos)
            };
            if !has_var_after_clip {
                let var_pos = clip_pos + CLIP_OFFSET;
                let ref_base = seq[(var_pos - self.offset()) as usize];
                let var_base =
                    ACGT.iter()
                        .copied()
                        .find(|&x| x != ref_base)
                        .ok_or(Exception::new(format!(
                    "Expected a variant base at 3' {clip_pos} position with var_pos = {var_pos}"
                )))?;
                let new_var = CandidateSite::new(var_pos, vec![ref_base], vec![var_base]);
                self.het_sites.push(new_var);
            }
        }

        self.het_sites.sort(); // Sort + dedup
        self.het_sites.dedup();

        raw_read_haps = self.haplotypes_from_reads_step(
            exclude_reads,
            min_clip_len,
            check_clip,
            partial_deletion_reads,
            min_mapq,
            tid,
            clip_buffer,
            Some(&absent_base_per_site),
        )?;

        Ok(raw_read_haps)
    }

    /// Add homozygous sites to fingerprinting.
    #[must_use]
    pub(crate) fn add_hom_sites(
        &mut self,
        min_no_var_region_size: Option<i64>,
        max_hom_var_to_add: Option<usize>,
        ref_seq: &[u8],
    ) -> Vec<CandidateSite> {
        let min_no_var_region_size = min_no_var_region_size.unwrap_or(10000);
        let max_hom_var_to_add = max_hom_var_to_add.unwrap_or(10);
        let mut ret = Vec::<CandidateSite>::new();
        self.het_sites.sort_by(|a, b| a.pos.cmp(&b.pos));
        self.hom_sites.sort_by(|a, b| a.pos.cmp(&b.pos));
        let het_pos = self
            .het_sites
            .iter()
            .filter(|x| x.pos < self.right_boundary_0based() && x.pos > self.left_boundary_0based())
            .map(|x| x.pos)
            .collect::<Vec<_>>();
        if het_pos.is_empty() {
            self.hom_sites
                .iter()
                .filter(|site| site.is_unit_length())
                .for_each(|x| ret.push(x.clone()));
            if self.hom_sites.is_empty() {
                // randomly pick a few non-variant sites
                let full_range = self.right_boundary_0based() - self.left_boundary_0based();
                let interval_size = full_range / 4;
                let positions = vec![
                    self.left_boundary_0based() + interval_size,
                    self.left_boundary_0based() + interval_size * 2,
                    self.left_boundary_0based() + interval_size * 3,
                ];
                for var_pos in positions {
                    let pos_on_ref = var_pos - self.offset();
                    let ref_base_u8 = ref_seq[pos_on_ref as usize];
                    let non_ref_bases = vec![b'A', b'C', b'G', b'T']
                        .iter()
                        .filter(|x| **x != ref_base_u8)
                        .map(|x| *x)
                        .collect::<Vec<_>>();
                    let var_base_u8 = non_ref_bases.first().unwrap();
                    let ref_base = std::str::from_utf8(&[ref_base_u8]).unwrap().to_string();
                    let var_base = std::str::from_utf8(&[*var_base_u8]).unwrap().to_string();
                    let new_variant = CandidateSite::new(var_pos, ref_base, var_base);
                    ret.push(new_variant);
                }
            }
        } else {
            let min_pos = het_pos.iter().min().unwrap();
            let max_pos = het_pos.iter().max().unwrap();
            // left side
            if *min_pos - self.left_boundary_0based() > min_no_var_region_size.into() {
                self.hom_sites
                    .iter()
                    .filter(|site| site.pos < *min_pos && site.is_unit_length())
                    .for_each(|x| ret.push(x.clone()));
            }
            // right side
            if self.right_boundary_0based() - *max_pos > min_no_var_region_size.into() {
                self.hom_sites
                    .iter()
                    .filter(|site| site.pos > *max_pos && site.is_unit_length())
                    .for_each(|x| ret.push(x.clone()));
            }
            // between two het site
            let het_sites_no_del = self
                .het_sites
                .iter()
                .filter(|x| x.is_unit_length())
                .collect::<Vec<_>>();
            let het_site_num = het_sites_no_del.len();
            for i in 0..(het_site_num - 1) {
                let interval_start = het_sites_no_del[i].pos;
                let interval_end = het_sites_no_del[i + 1].pos;
                if interval_end - interval_start > min_no_var_region_size {
                    self.hom_sites
                        .iter()
                        .filter(|site| {
                            site.pos < interval_end
                                && site.pos > interval_start
                                && site.is_unit_length()
                        })
                        .for_each(|x| ret.push(x.clone()));
                }
            }
        }
        // select variants evenly. TODO: evenly based on coordinates, rather than number of sites
        if !ret.is_empty() {
            ret.sort();
            let num_sites = ret.len();
            for hom_site in ret
                .iter()
                .step_by((num_sites + max_hom_var_to_add - 1) / max_hom_var_to_add)
            {
                log::debug!("Adding hom site: {hom_site:?}");
                self.het_sites.push((*hom_site).clone());
            }
            self.het_sites.sort();
        }
        ret.into_iter().map(|x| x.clone()).collect::<Vec<_>>()
    }

    /// removes variant siters within any regions marked as `noisy_regions`.
    pub(crate) fn remove_noisy_sites(&mut self) {
        self.het_sites.retain(|site| {
            // Keep only if it overlaps with none of the noisy regions.
            !self
                .noisy_regions
                .iter()
                .any(|region| region.contains(&site.pos))
        });
    }

    /// Yields a reference to the `LocusConfig`
    #[must_use]
    pub fn locus_config(&self) -> &LocusConfig {
        &self.config.locus
    }

    /// Yields a reference to the `GeneConfig`
    #[must_use]
    pub fn gene_config(&self) -> &GeneConfig {
        &self.config.gene
    }

    pub fn get_default_call(&self) -> GeneCall {
        let region_median_depth = self
            .region_avg_depth
            .iter()
            .map(|x| x.0)
            .collect::<Vec<_>>();
        let region_80percentile_depth = self
            .region_avg_depth
            .iter()
            .map(|x| x.1)
            .collect::<Vec<_>>();
        let mut region_depth = BTreeMap::new();
        region_depth.insert(
            String::from("median"),
            region_median_depth.first().unwrap().clone(),
        );
        region_depth.insert(
            String::from("percentile80"),
            region_80percentile_depth.first().unwrap().clone(),
        );
        GeneCall {
            gene_name: self.settings.gene_name.clone(),
            sample_sex: self.settings.sample_sex.to_string(),
            genome_depth: self.settings.depth.map(|x| x.median as f32),
            region_depth,
            phase_region: format!(
                "{}:{}:{}-{}",
                self.settings.genome,
                self.chr().unwrap(),
                self.left_boundary(),
                self.right_boundary()
            ),
            ..Default::default()
        }
    }

    /// Parse two predefined deletions
    pub fn parse_deletions_from_config(&mut self) -> DResult {
        // get predefined deletion
        let deletion1_name = self.locus_config().get("deletion1_name");
        if !deletion1_name.is_none() {
            let deletion1_name = deletion1_name
                .and_then(|x| x.as_str())
                .ok_or("Missing deletion1_name for gene config")?
                .to_string();
            let deletion1_start = deletion1_name
                .split("_")
                .map(std::borrow::ToOwned::to_owned)
                .collect::<Vec<_>>()
                .first()
                .ok_or("first not in split variant name")?
                .parse::<i64>()?
                - 1;
            let deletion1_size = self
                .locus_config()
                .get("deletion1_size")
                .and_then(|x| {
                    x.as_str()
                        .and_then(|s| s.parse::<i64>().ok())
                        .or(x.as_i64())
                })
                .ok_or("Missing deletion1_size for gene config")?;
            let del1_3p_pos1 = self
                .locus_config()
                .get("del1_3p_pos1")
                .and_then(|x| {
                    x.as_str()
                        .and_then(|s| s.parse::<i64>().ok())
                        .or(x.as_i64())
                })
                .ok_or("Missing del1_3p_pos1 for gene config")?;
            let del1_3p_pos2 = self
                .locus_config()
                .get("del1_3p_pos2")
                .and_then(|x| {
                    x.as_str()
                        .and_then(|s| s.parse::<i64>().ok())
                        .or(x.as_i64())
                })
                .ok_or("Missing del1_3p_pos2 for gene config")?;
            let del1_5p_pos1 = self
                .locus_config()
                .get("del1_5p_pos1")
                .and_then(|x| {
                    x.as_str()
                        .and_then(|s| s.parse::<i64>().ok())
                        .or(x.as_i64())
                })
                .ok_or("Missing del1_5p_pos1 for gene config")?;
            let del1_5p_pos2 = self
                .locus_config()
                .get("del1_5p_pos2")
                .and_then(|x| {
                    x.as_str()
                        .and_then(|s| s.parse::<i64>().ok())
                        .or(x.as_i64())
                })
                .ok_or("Missing del1_5p_pos2 for gene config")?;

            let deletion_end = deletion1_start + deletion1_size;
            let new_del_data = DeletionDatum::new(
                range::I64::new(deletion1_start, deletion_end),
                None,
                Some(del1_3p_pos1),
                Some(del1_3p_pos2),
                Some(del1_5p_pos1),
                Some(del1_5p_pos2),
            );
            self.del_data.push(new_del_data);
        }
        let deletion2_name = self.locus_config().get("deletion2_name");
        if !deletion2_name.is_none() {
            let deletion2_name = deletion2_name
                .and_then(|x| x.as_str())
                .ok_or("Missing deletion2_name for gene config")?
                .to_string();
            let deletion2_start = deletion2_name
                .split("_")
                .map(std::borrow::ToOwned::to_owned)
                .collect::<Vec<_>>()
                .first()
                .ok_or("first not in split variant name")?
                .parse::<i64>()?
                - 1;
            let deletion2_size = self
                .locus_config()
                .get("deletion2_size")
                .and_then(|x| {
                    x.as_str()
                        .and_then(|s| s.parse::<i64>().ok())
                        .or(x.as_i64())
                })
                .ok_or("Missing deletion2_size for gene config")?;
            let del2_3p_pos1 = self
                .locus_config()
                .get("del2_3p_pos1")
                .and_then(|x| {
                    x.as_str()
                        .and_then(|s| s.parse::<i64>().ok())
                        .or(x.as_i64())
                })
                .ok_or("Missing del2_3p_pos1 for gene config")?;
            let del2_3p_pos2 = self
                .locus_config()
                .get("del2_3p_pos2")
                .and_then(|x| {
                    x.as_str()
                        .and_then(|s| s.parse::<i64>().ok())
                        .or(x.as_i64())
                })
                .ok_or("Missing del2_3p_pos2 for gene config")?;
            let del2_5p_pos1 = self
                .locus_config()
                .get("del2_5p_pos1")
                .and_then(|x| {
                    x.as_str()
                        .and_then(|s| s.parse::<i64>().ok())
                        .or(x.as_i64())
                })
                .ok_or("Missing del2_5p_pos1 for gene config")?;
            let del2_5p_pos2 = self
                .locus_config()
                .get("del2_5p_pos2")
                .and_then(|x| {
                    x.as_str()
                        .and_then(|s| s.parse::<i64>().ok())
                        .or(x.as_i64())
                })
                .ok_or("Missing del2_5p_pos2 for gene config")?;

            let deletion_end = deletion2_start + deletion2_size;
            let new_del_data = DeletionDatum::new(
                range::I64::new(deletion2_start, deletion_end),
                None,
                Some(del2_3p_pos1),
                Some(del2_3p_pos2),
                Some(del2_5p_pos1),
                Some(del2_5p_pos2),
            );
            self.del_data.push(new_del_data);
        }
        log::debug!("del_data {:?}", self.del_data);
        Ok(())
    }

    ///
    /// Run Complete Phaser workflow.
    ///
    /// # Panics
    /// 1. Assert failure: if genome bam does not exist.
    /// 2. Downstream panics.
    pub fn run(&mut self) -> Result<GeneCall, DError> {
        log::trace!("Running caller using settings {:?}", self.settings);
        // Initial setup:
        // S1: Get local region.
        let seq = self.realign()?;
        // Check coverage after aligning to local reference.
        let coverage_passes = self.coverage_passes();
        let mut call = self.get_default_call();
        if !coverage_passes {
            log::debug!("Call failed for coverage");
            call.failed_for_coverage = true;
            return Ok(call);
        }

        let (hom_sites_to_add, add_sites) = self.get_sites(&seq, None, None)?;
        let tid = self.genome_tid().map(|x| x as i32).expect("No chr tid");
        let init_read_hap_map = self.haplotypes_from_reads(
            None,
            /* kept_sites */ &hom_sites_to_add,
            Some(&add_sites),
            /* partial_deletion_reads */ None,
            (
                /* min_mapq= */ 5,
                /* check_clip= */ true,
                /* min_clip_len */ Some(4500u32),
            ),
            tid,
            None,
            &hom_sites_to_add,
        )?;
        let mut phase_results: PhasedResult;
        let known_del: BTreeMap<char, String>;
        match self.update_indel_and_phase(init_read_hap_map.clone(), &mut call) {
            Ok(results) => {
                phase_results = results.0;
                known_del = results.1;
            }
            Err(e) => {
                log::debug!(
                    "Failed to update indel and phase with the following error: {:?}",
                    e
                );
                call.sites_for_phasing = self
                    .het_sites
                    .iter()
                    .map(ToString::to_string)
                    .collect::<Vec<_>>();
                self.init_het_sites.sort_by(|a, b| a.pos.cmp(&b.pos));
                call.heterozygous_sites = self
                    .init_het_sites
                    .iter()
                    .map(ToString::to_string)
                    .collect::<Vec<_>>();
                call.homozygous_sites = self
                    .hom_sites
                    .iter()
                    .map(ToString::to_string)
                    .collect::<Vec<_>>();
                self.het_sites_no_phasing.sort_by(|a, b| a.pos.cmp(&b.pos));
                call.het_sites_not_used_in_phasing = self
                    .het_sites_no_phasing
                    .iter()
                    .map(ToString::to_string)
                    .collect::<Vec<_>>();
                return Ok(call);
            }
        }

        // rename haplotypes
        let mut assembled_haps = BTreeMap::new();
        let main_haps_clone = phase_results.assemblies.main_haps.clone();
        let mod_gene_name =
            intersperse(self.gene_name().split_terminator('-'), ",").collect::<String>();
        for (idx, hap) in main_haps_clone.iter().enumerate() {
            assembled_haps.insert(hap.vstr(), format!("{mod_gene_name}_hap{}", idx + 1));
        }
        call.final_haplotypes = assembled_haps
            .clone()
            .into_iter()
            .map(|(k, v)| (k.to_string(), v))
            .collect::<BTreeMap<_, _>>();

        // Output variants
        let haps =
            self.output_variants_in_haps(&phase_results, &known_del, assembled_haps.clone())?;
        call.haplotype_details = haps
            .iter()
            .map(|(key, val)| (key.clone(), HapInfoForJson::from(val)))
            .collect::<BTreeMap<_, _>>();

        // report
        self.fill_in_call(phase_results, &mut call);
        Ok(call)
    }
}

pub type Exception = simple_error::SimpleError;

/// Takes a deletion length + count field, extracts the digit-only portion, and then converts it to an integer.
/// "1N"
#[must_use]
fn extract_del_length(x: &[u8]) -> i32 {
    let non_digit_idx = x
        .iter()
        .position(|&x| !x.is_ascii_digit())
        .unwrap_or(x.len());
    debug_assert!(!x.is_empty());
    // Technically, if it's guaranteed to be nonzero in lengthm, you could fold from the first digit.
    // But it shouldn't really matter.
    x[..non_digit_idx]
        .iter()
        .fold(0i32, |acc: i32, x: &u8| acc * 10 + i32::from(x - b'0'))
}

/// Takes a slice and returns a slice for the remainder of the string following all numeric characters.
/// ```
/// use paraphase::phaser::skip_digits;
/// assert_eq!(skip_digits(b"123ABC"), b"ABC");
/// assert_eq!(skip_digits(b"ABC"), b"ABC");
/// assert_eq!(skip_digits(b""), b"");
/// assert_eq!(skip_digits(b"1232322"), b"");
/// ```
#[must_use]
pub fn skip_digits(x: &[u8]) -> &[u8] {
    let seq_start = x
        .iter()
        .position(|&x| !x.is_ascii_digit())
        .unwrap_or(x.len());
    &x[seq_start..]
}

fn counts(
    x: &bam::pileup::Pileup,
    ref_seq: &[u8],
    min_base_quality: u8,
    offset: usize,
) -> counter::Counter<VString, i32> {
    use crate::detail::phaser_util::base_qual;
    let mut ret = counter::Counter::new();
    let mark_strand = false;
    for aln in x.alignments() {
        let record = aln.record();
        let query_pos = util::raw_qpos(&aln);
        let seq = record.seq().as_bytes();
        if base_qual(&aln) < min_base_quality {
            continue;
        }
        let is_reverse = record.is_reverse();
        let refskip_char = if is_reverse { b'<' } else { b'>' };
        let base = if !aln.is_del() && !aln.is_refskip() {
            seq.get(query_pos).copied().unwrap_or(b'N')
        } else if aln.is_refskip() {
            refskip_char
        } else {
            b'*'
        };
        let pos = record.pos() as usize;
        let mut query_seq = VString::from(vec![base]);
        match aln.indel() {
            bam::pileup::Indel::Ins(x) => {
                query_seq.push(b'+');
                query_seq.extend_from_slice(x.to_string().as_bytes());
                for j in 1..=(x as usize) {
                    query_seq.push(maybe_strand_mark_char(
                        seq[j + query_pos],
                        is_reverse,
                        mark_strand,
                    ));
                }
            }
            bam::pileup::Indel::Del(x) => {
                query_seq.push(b'-');
                query_seq.extend_from_slice(x.to_string().as_bytes());
                for j in 1..=(x as usize) {
                    query_seq.push(maybe_strand_mark_char(
                        ref_seq[j + pos - offset],
                        is_reverse,
                        mark_strand,
                    ));
                }
            }
            bam::pileup::Indel::None => {}
        }
        *ret.entry(query_seq).or_default() += 1;
    }
    ret
}

/// Helper function to check if two ranges overlap
pub fn ranges_overlap(start1: i64, end1: i64, start2: i64, end2: i64) -> bool {
    start1 < end2 && start2 < end1
}

#[cfg(test)]
mod tests {

    #[test]
    fn test_ranges_overlap() {
        assert!(crate::phaser::ranges_overlap(0, 100, 50, 150));
        assert!(!crate::phaser::ranges_overlap(0, 100, 100, 200));
        assert!(crate::phaser::ranges_overlap(0, 100, 50, 75));
    }

    #[test]
    fn skip_digits_ok() {
        assert_eq!(crate::phaser::skip_digits(b"ACGT"), b"ACGT");
        assert_eq!(crate::phaser::skip_digits(b"123ACGT"), b"ACGT");
        assert_eq!(crate::phaser::skip_digits(b"123"), &[] as &[u8]);
    }

    #[test]
    fn extract_del_len_ok() {
        use crate::phaser::extract_del_length;
        assert_eq!(extract_del_length(b"1N"), 1);
        assert_eq!(extract_del_length(b"100N"), 100);
        assert_eq!(extract_del_length(b"37FD34"), 37);
    }

    #[test]
    fn process_indel_ok() {
        use crate::detail::util;
        use crate::phaser::Phaser;
        let faidx = rust_htslib::faidx::Reader::from_path(util::test_file("smn1_ref.fa")).unwrap();
        let all_ref_seq = util::load_all_seqs(&faidx);
        let offset = 70_890_000;
        // both offset and pos are 1-based, which is okay.

        let (ref_seq, var_seq, indel_len) = Phaser::free_process_indel(
            70_940_935 - offset,
            b"A",
            b"A+3CCC",
            //70_889_999,88
            &all_ref_seq[0],
        )
        .unwrap();
        assert_eq!(*ref_seq, &b"A"[..]);
        assert_eq!(*var_seq, &b"ACCC"[..]);
        assert_eq!(indel_len, 3, "Insertion has wrong length.");

        let (ref_seq, var_seq, indel_len) = Phaser::free_process_indel(
            70_940_935 - offset,
            b"A",
            b"A-2NN",
            //70_889_999,
            &all_ref_seq[0],
        )
        .unwrap();
        assert_eq!(*ref_seq, &b"ACT"[..]);
        assert_eq!(*var_seq, &b"A"[..]);
        assert_eq!(indel_len, 2, "Deletion has wrong length");
    }
}
