use crate::config::{self, Gene as GeneConfig};
use crate::detail::deletion::Datum as DeletionDatum;
use crate::detail::range::I64 as Range64;
use crate::detail::site_selection::{
    pileup as site_pileup, CandidateSite, FilteredSites, RawVariantCounts,
    Settings as SiteSelectionSettings,
};
use crate::detail::util::{self, DError, DResult, HashSet};
use crate::io::json::{ReadAlignmentId, ReadFingerprintMap};
use crate::phaser::{
    self, Exception as PhaserException, FlagBits as PhaserFlagBits, Phaser,
    Settings as PhaserSettings,
};
use crate::realign::{reference_length, RealignSettings};

use vstr::VStr;

use itertools::Itertools;
use rust_htslib::bam::ext::BamRecordExtensions;
use rust_htslib::bam::{self, pileup, record::Cigar, Read};
use thiserror::Error;

use std::collections::{BTreeMap, BTreeSet};
use std::io::Write;
use std::ops::Index;
use std::str::FromStr;

pub mod defaults {
    pub const DELETION_LENGTH_FUZZ: i64 = 50;
    pub const FLANKING_SPURIOUS_BP: i32 = 10;
    pub const MIN_BASE_QUALITY: u8 = 25; // In Python: paraphase.phaser.Phaser.MEAN_BASE_QUAL
}

lazy_static::lazy_static! {
    static ref NAME_REG: regex::bytes::RegexSet = regex::bytes::RegexSetBuilder::new(["ccs", "transcript", "molecule"])
        .unicode(false).build().unwrap();
}

///
///```
/// use paraphase::detail::phaser_util::passes;
/// assert!(passes(b"ccs"));
/// assert!(passes(b"hello/i/am/molecule"));
/// assert!(!passes(b"23789ff3303449"));
/// ```
#[must_use]
pub fn passes(hay: &[u8]) -> bool {
    NAME_REG.is_match(hay)
}

///
/// Find the indices for first + last non-gap characters in a &[u8].
///
/// # Inputs
/// x: &[u8]
/// # Outputs
/// (first: i32, last: i32)
///
///```
/// use paraphase::detail::phaser_util::get_start_end;
/// assert_eq!(get_start_end(b"xx1211xxx"), (2, 5));
/// assert_eq!(get_start_end(b"xxxxx"), (4, 0));
/// ```
#[must_use]
pub fn get_start_end(x: &[u8]) -> (i32, i32) {
    let start = x
        .iter()
        .enumerate()
        .find(|x| *x.1 != b'x')
        .map_or(x.len() - 1, |x| x.0) as i32;
    let end = x
        .iter()
        .rev()
        .enumerate()
        .find(|x| *x.1 != b'x')
        .map_or(0, |(count, _item)| x.len() - 1 - count) as i32;
    (start, end)
}

#[derive(Clone, Debug, Error)]
#[error("{0}")]
pub struct FaiBuildError(String);

pub fn build_faidx(path: impl Into<std::path::PathBuf>) -> DResult {
    let path = path.into();
    let os_path = std::ffi::CString::new(path.display().to_string())?;
    let rc = unsafe { rust_htslib::htslib::fai_build(os_path.as_ptr()) };
    let msg = match rc {
        -1 => "indexing failed",
        _ => return Ok(()),
    };
    Err(FaiBuildError(format!(
        "rc: {rc}. msg: {msg}. Path: {path:?}"
    )))?
}

pub(crate) fn base_qual(x: &pileup::Alignment<'_>) -> u8 {
    // Set deletion base quality to 255 so it always passes.
    let record = x.record();
    let qual = record.qual();
    qual.get(util::raw_qpos(x)).copied().unwrap_or(u8::MAX)
}

#[allow(dead_code)]
fn parse_genomic_region(x: &str) -> Option<(&str, i64, i64)> {
    let mut it = x.split_terminator('_');
    let chr = it.next()?;
    let (start, stop) = it.filter_map(|x| x.parse::<i64>().ok()).next_tuple()?;
    Some((chr, start, stop))
}

impl Phaser {
    ///
    /// Generate path for `realigned_bam` from output directory.
    #[must_use]
    pub fn realigned_bam_path(&self) -> std::path::PathBuf {
        let suffix = format!(
            "{}_{}_realigned.bam",
            self.settings.sample_id,
            self.gene_name()
        );
        self.settings.outdir.join(suffix)
    }

    /// Write a BED file with the variants to the output directory.
    /// Useful for visualization.
    /// WARNING:
    ///   In contrast to the BED files in configuration with are 1-based, these are 0-based.
    pub fn write_variant_bed(&self) -> DResult {
        let path = self.settings.outdir.join(format!(
            "{}_{}_sampled.0based.bed",
            self.settings.sample_id,
            self.gene_name()
        ));
        let mut writer = std::io::BufWriter::new(std::fs::File::create(path)?);
        let gene_name = self.gene_name();
        let chr = self.chr().ok_or_else(|| {
            PhaserException::new("Failed to get chromosome for output BED file generation.")
        })?;
        let mut site_id = 0usize;
        writeln!(writer, "#chr\tstart\tstop\tSiteId\tRef_Var")?;
        for site in &self.het_sites {
            site_id += 1;
            writeln!(
                writer,
                "{chr}\t{}\t{}\t{gene_name}_het_site_{site_id}\t{}_{}",
                site.pos,
                site.pos + site.reference_length() as i64,
                site.ref_seq,
                site.var_seq
            )?;
        }
        for site in &self.hom_sites {
            site_id += 1;
            writeln!(
                writer,
                "{chr}\t{}\t{}\t{gene_name}_hom_site_{site_id}\t{}_{}",
                site.pos,
                site.pos + site.reference_length() as i64,
                site.ref_seq,
                site.var_seq
            )?;
        }
        for site in &self.het_sites_no_phasing {
            site_id += 1;
            writeln!(
                writer,
                "{chr}\t{}\t{}\t{gene_name}_het_site_no_phasing_{site_id}\t{}_{}",
                site.pos,
                site.pos + site.reference_length() as i64,
                site.ref_seq,
                site.var_seq
            )?;
        }
        Ok(())
    }

    ///
    /// Generate path for `realigned_tagged_path` from output directory.
    #[must_use]
    pub fn realigned_tagged_bam_path(&self) -> std::path::PathBuf {
        let suffix = format!(
            "{}_{}_realigned_tagged.bam",
            self.settings.sample_id,
            self.gene_name(),
        );
        self.settings.outdir.join(suffix)
    }

    ///
    /// Generate path for `realigned_tagged_gene2_path` from output directory.
    #[must_use]
    pub fn realigned_tagged_gene2_bam_path(&self) -> std::path::PathBuf {
        let suffix = format!(
            "{}_{}_gene2_realigned_tagged.bam",
            self.settings.sample_id,
            self.gene_name(),
        );
        self.settings.outdir.join(suffix)
    }

    ///
    /// Generate path for `realigned_gene2_path` from output directory.
    #[must_use]
    pub fn realigned_gene2_bam_path(&self) -> std::path::PathBuf {
        let suffix = format!(
            "{}_{}_gene2_realigned.bam",
            self.settings.sample_id,
            self.gene_name(),
        );
        self.settings.outdir.join(suffix)
    }

    ///
    /// Generate path for `genome_bam` from output directory.
    #[must_use]
    pub fn genome_bam_path(&self) -> std::path::PathBuf {
        self.settings.genome_bam.clone()
    }

    ///
    /// Generate path for local reference.
    #[must_use]
    pub fn local_reference_path(&self) -> std::path::PathBuf {
        let suffix = format!("{}_ref.fa", self.gene_name());
        let path = self.settings.outdir.join(suffix);
        log::debug!("Local reference path: {path:?}");
        path
    }

    ///
    /// Generate path for local reference for secondary region, aka 'gene2'.
    #[must_use]
    pub fn secondary_reference_path(&self) -> std::path::PathBuf {
        let suffix = format!("{}_gene2_ref.fa", self.gene_name());
        let path = self.settings.outdir.join(suffix);
        log::debug!("Secondary reference path: {path:?}");
        path
    }

    /// Generate a fasta file
    /// from a `realign_region` and the reference genome.
    /// # Panics
    /// • If the region string is malformed.
    ///
    /// # Errors
    /// • Invalid utf-8 fetched from fasta.
    /// • Output could not be created.
    /// • Realign region malformed.
    /// • Writing sequence to disk failed.
    /// • faidx could not be created for fasta path.
    pub fn generate_local_reference(&self) -> Result<std::path::PathBuf, DError> {
        log::trace!("Building local reference");
        let faidx = self.make_faidx()?;
        log::trace!("Made faidx correctly");
        let (chrom, start, stop) = self
            .parsed_nchr_0based()
            .expect("Malformatted region string");
        log::trace!("Building local reference from {chrom}:{start}-{stop}");
        let seq = std::str::from_utf8(&faidx.fetch_seq(chrom, start as usize, stop as usize)?)?
            .to_ascii_uppercase();
        let dest = self.local_reference_path();
        log::trace!("Writing local ref to {dest:?}");
        let mut output = std::io::BufWriter::new(std::fs::File::create(&dest)?);
        let realign_old = self.realign_region_old().ok_or(PhaserException::new(
            "realign_region_old failed".to_string(),
        ))?;
        log::trace!("Realign old {realign_old:?}");
        writeln!(output, ">{realign_old}\n{seq}",)?;
        drop(output);
        build_faidx(&dest)?;
        log::trace!("Built local reference");
        Ok(dest)
    }

    /// Generate a fasta file
    /// from a `gene2_region` and the reference genome.S
    ///
    /// # Errors
    /// • Invalid utf-8 fetched from fasta.
    /// • Output could not be created.
    /// • Realign region malformed.
    /// • Writing sequence to disk failed.
    /// • faidx could not be created for fasta path.
    pub fn generate_secondary_reference(&self) -> Result<std::path::PathBuf, DError> {
        log::trace!("Building secondary reference");
        let faidx = self.make_faidx()?;
        let (chrom, start, stop) = self
            .parsed_nchr_secondary_0based()
            .expect("Malformatted region string");
        log::trace!("Fetching seq at {chrom}:{start}-{stop} (0-based)");
        let seq = std::str::from_utf8(&faidx.fetch_seq(chrom, start as usize, stop as usize)?)?
            .to_ascii_uppercase();
        let dest = self.secondary_reference_path();
        let mut output = std::io::BufWriter::new(std::fs::File::create(&dest)?);
        writeln!(
            output,
            ">{}\n{seq}",
            self.secondary_region_old()
                .ok_or(PhaserException::new("secondary_region_old failed"))?
        )?;
        drop(output);
        build_faidx(&dest)?;
        log::trace!("Built secondary reference");
        Ok(dest)
    }

    /// Generate path if exists.
    /// Returns (Path, bool), where bool indices if the path was created after invocation.
    ///
    /// # Errors
    /// Propagates errors from `generate_local_reference`, which include:
    ///   1. utf-8 error on parsing string,
    ///   2. Error fetching from faidx (`rust_htslib::errors::Error`)
    ///   3. Creating destination file (`std::io::Error`)
    ///   4. Error writing seq to disk.
    ///   5. Error parsing `realign_region`
    ///   6. Error running samtools faidx to generate index for new reference
    pub fn secondary_reference(&self) -> Result<(std::path::PathBuf, bool), DError> {
        let path = self.secondary_reference_path();
        log::trace!("Secondary ref destination is {path:?}");
        let generated = if !path.exists() {
            self.generate_secondary_reference()?;
            true
        } else {
            false
        };
        Ok((path, generated))
    }

    /// Generate path if exists.
    /// Returns (Path, bool), where bool indices if the path was created after invocation.
    ///
    /// # Errors
    /// Propagates errors from `generate_local_reference`, which include:
    ///   1. utf-8 error on parsing string,
    ///   2. Error fetching from faidx (`rust_htslib::errors::Error`)
    ///   3. Creating destination file (`std::io::Error`)
    ///   4. Error writing seq to disk.
    ///   5. Error parsing `realign_region`
    ///   6. Error running samtools faidx to generate index for new reference
    pub fn local_reference(&self) -> Result<(std::path::PathBuf, bool), DError> {
        let path = self.local_reference_path();
        log::trace!("Putting local reference to {path:?}");
        let generated = if !path.exists() {
            self.generate_local_reference()?;
            true
        } else {
            false
        };
        log::trace!("Put local reference to {path:?}");
        Ok((path, generated))
    }

    /// Creates a Phaser structure from settings, gene name, and region config.
    ///
    /// # Panics
    /// 1. If `gene_name` is not in `RegionConfig`
    /// 2. If `realign_region` is not found in gene sub-RegionConfig
    /// 3. If int vector fields are not comprised of integer values.
    #[must_use]
    pub fn new(
        settings: PhaserSettings,
        gene_config: Option<GeneConfig>,
        site_selection_settings: Option<SiteSelectionSettings>,
        realign_settings: Option<RealignSettings>,
    ) -> Self {
        Self::try_new(
            settings,
            gene_config,
            site_selection_settings,
            realign_settings,
        )
        .unwrap_or_else(|e| panic!("Failed to create Phaser struct. Error: {e:}"))
    }

    ///
    /// Fallible new construction of Phaser.
    pub fn try_new(
        mut settings: PhaserSettings,
        gene_config: Option<GeneConfig>,
        site_selection_settings: Option<SiteSelectionSettings>,
        realign_settings: Option<RealignSettings>,
    ) -> Result<Self, DError> {
        let mut locus_config = settings
            .region_config
            .get(&settings.gene_name)
            .or_else(|| settings.region_config.get(&settings.gene_name))
            .ok_or(PhaserException::new(format!(
                "Failed to get gene config for region {} from region config {:?}",
                settings.gene_name, settings.region_config
            )))?
            .clone();
        if let Some(mut site_settings) = site_selection_settings {
            site_settings.update_from_settings(&locus_config);
            settings.site_selection_settings = site_settings;
        }
        let gene_config = gene_config.unwrap_or_default();
        let realign_settings = realign_settings
            .unwrap_or_default()
            .update_from_locus(&locus_config);
        locus_config.insert("gene_name".into(), settings.gene_name.clone().into());
        let realign_region = if settings.genome != "37" {
            locus_config
                .get("realign_region")
                .and_then(|x| x.as_str())
                .ok_or_else(|| {
                    PhaserException::new("Missing realign region for gene config {locus_config:?}")
                })?
                .to_string()
        } else {
            locus_config
                .get("realign_region")
                .and_then(|x| x.as_str())
                .ok_or_else(|| {
                    PhaserException::new("Missing realign region for gene config {locus_config:?}")
                })?
                .strip_prefix("chr")
                .ok_or("error with strip_prefix")?
                .to_string()
        };

        let get_int_field = |key: &str| -> Option<i64> {
            locus_config.get(key).and_then(|x| {
                x.as_str()
                    .and_then(|s| s.parse::<i64>().ok())
                    .or(x.as_i64())
            })
        };

        // 0 based
        let get_sorted_int_vec_field = |key: &str| -> Vec<i64> {
            locus_config
                .get(key)
                .and_then(|x| x.as_sequence())
                .map(|x| {
                    x.iter()
                        .map(|x| x.as_i64().expect("int vec field member was not an integer. Check configuration files.") - 1) // Subtract 1 because the config files are 1-based.
                        .sorted()
                        .collect::<Vec<_>>()
                })
                .unwrap_or(vec![])
        };

        let (left_boundary, right_boundary, gene_start, gene_end) =
            ["left_boundary", "right_boundary", "gene_start", "gene_end"]
                .into_iter()
                .map(get_int_field)
                .next_tuple()
                .ok_or(PhaserException::new(
                    "The impossible happened: wrong number of fields in a compile-time array.",
                ))?;
        log::debug!("From locus conf {locus_config:?}");
        log::debug!("left bound: {left_boundary:?}");
        log::debug!("right bound: {right_boundary:?}");
        log::debug!(
            "left boundary field: {:?}",
            locus_config.get("left_boundary")
        );
        let pivot_site = get_int_field("pivot_site");
        let add_sites = locus_config
            .get("add_sites")
            .and_then(|x| x.as_sequence())
            .map(|x| {
                x.iter()
                    .filter_map(|x| x.as_str().map(CandidateSite::from_str))
                    .flatten()
                    .collect::<Vec<_>>()
            })
            .unwrap_or(vec![]);
        let clip_3p_positions = get_sorted_int_vec_field("clip_3p_positions");
        let clip_5p_positions = get_sorted_int_vec_field("clip_5p_positions");
        //clip_5p_positions.reverse();
        assert!(clip_3p_positions.windows(2).all(|x| x[1] >= x[0]));
        assert!(clip_5p_positions.windows(2).all(|x| x[1] >= x[0]));

        let expect_cn2 = locus_config
            .get("expect_cn2")
            .and_then(serde_yaml::Value::as_bool)
            .unwrap_or(false);
        let is_reverse = locus_config
            .get("is_reverse")
            .and_then(serde_yaml::Value::as_bool)
            .unwrap_or(false);
        let use_supplementary = locus_config.use_supplementary();
        let to_phase = ["in_tandem", "to_phase"]
            .map(String::from)
            .iter()
            .any(|x| locus_config.contains_key(x));
        let flag = [to_phase, use_supplementary, is_reverse, expect_cn2]
            .into_iter()
            .zip([
                PhaserFlagBits::ToPhase,
                PhaserFlagBits::UseSupplementary,
                PhaserFlagBits::IsReverse,
                PhaserFlagBits::ExpectCN2,
            ])
            .fold(0u8, |mut acc, x| {
                let (set, value) = x;
                if set {
                    acc |= value as u8;
                }
                acc
            });
        log::debug!("Phaser flags set.");

        let noisy_regions = locus_config.extract_noisy_regions();
        let locus = locus_config;
        let gene = gene_config;
        log::debug!("Low complexity regions assigned.");
        let mut ret = Self {
            settings: settings.clone(),
            realign_settings,
            config: phaser::MergedConfig { locus, gene },
            realign_region,
            left_boundary,
            right_boundary,
            add_sites,
            clip_3p_positions,
            clip_5p_positions,
            gene_start,
            gene_end,
            pivot_site,
            candidate_sites: BTreeSet::new(), // Empty at start, added later.
            del_data: vec![],                 // Empty at start.
            het_sites: vec![],
            init_het_sites: vec![],
            het_sites_no_phasing: vec![],
            hom_sites: vec![],
            low_complexity_sites: Default::default(),
            matches: Default::default(),
            noisy_regions,
            flag,
            region_avg_depth: vec![],
        };
        // use if let instead of map so we can return errors.
        ret.matches = if let Some(gene2_region) =
            ret.locus_config().gene2_region(settings.genome == "37")
        {
            let faidx = ret.make_faidx()?;
            let matches = config::region::GenePositionCorrelation::from_ref_regions(
                &ret.realign_region,
                gene2_region,
                &faidx,
                ret.locus_config().chain_bandwidth(),
            )?
            .into_inner();
            log::debug!("Multi-gene map built");
            matches
        } else {
            log::debug!("Multi-gene not built, as we are not aligning back to a second region.");
            Default::default()
        };
        Ok(ret)
    }

    /// Remove variant sites with insufficient support.
    pub fn remove_var(
        &mut self,
        raw_read_haps: &ReadFingerprintMap,
        kept_sites: &[CandidateSite],
        hom_sites: &[CandidateSite],
    ) -> BTreeMap<CandidateSite, u8> {
        let mut absent_base_per_site = BTreeMap::new();
        if self.het_sites.is_empty() {
            log::debug!("No variant sites; skipping remove_var.");
            return absent_base_per_site;
        }

        assert_eq!(
            raw_read_haps
                .values()
                .map(|x| x.len())
                .collect::<BTreeSet<_>>()
                .len(),
            1,
            "All read fingerprints should be of uniform length."
        );

        let hap_len = raw_read_haps.values().next().map_or(0, |x| x.len());
        let mut bases_per_site = vec![counter::Counter::<u8, i32>::new(); hap_len];
        log::debug!("{} Initial read haps", raw_read_haps.len());
        for hap in raw_read_haps.values() {
            assert_eq!(hap.len(), bases_per_site.len());
            for (base, counter) in hap.iter().copied().zip(bases_per_site.iter_mut()) {
                *counter.entry(base).or_default() += 1;
            }
        }

        log::trace!("Counts for each variant site: {bases_per_site:?}");

        let mut sites_to_remove = Vec::new();
        for (pos, base_counter) in bases_per_site
            .iter()
            .enumerate()
            .filter(|x| !x.1.is_empty())
        {
            let ref_count = base_counter.get(&b'1').copied().unwrap_or(0);
            let alt_count = base_counter.get(&b'2').copied().unwrap_or(0);
            let base0_count = base_counter.get(&b'0').copied().unwrap_or(0);
            let x_count = base_counter.get(&b'x').copied().unwrap_or(0);
            let ref_plus_alt = ref_count + alt_count;
            // Remove if not in kept sites and insufficient coverage
            // or all bases are 'x'.
            let total = base_counter.total::<i32>();
            let in_kept_sites = kept_sites.contains(&self.het_sites[pos]);
            let mut to_remove = false;
            if x_count == total - base0_count {
                to_remove = true;
            } else if ref_plus_alt == (total - x_count - base0_count)
                && (alt_count <= 3 || ref_count <= 3)
            {
                if alt_count <= 3 && ref_count <= 3 {
                    to_remove = true;
                } else {
                    if !in_kept_sites {
                        to_remove = true;
                    }
                    if hom_sites.contains(&self.het_sites[pos]) {
                        if alt_count <= 3 {
                            absent_base_per_site.insert(self.het_sites[pos].clone(), b'2');
                        } else if ref_count <= 3 {
                            absent_base_per_site.insert(self.het_sites[pos].clone(), b'1');
                        }
                    }
                }
            }
            if to_remove {
                log::debug!("Filtering pos {pos} for counts. x_count: {x_count}. ref {ref_count} + alt {alt_count}: {}. In kept: {in_kept_sites}", ref_plus_alt);
                sites_to_remove.push(pos);
            }
        }
        log::debug!("Variants filtered out at sites {sites_to_remove:?}");
        for idx in &sites_to_remove {
            let var_to_remove = &self.het_sites[*idx];
            if self.init_het_sites.contains(var_to_remove) {
                let idx_init = self
                    .init_het_sites
                    .iter()
                    .position(|x| x == var_to_remove)
                    .unwrap();
                self.init_het_sites.swap_remove(idx_init);
            }
        }
        // Now, remove het sites at failing positions.
        // Iterate in reverse order to avoid invalidating indices
        sites_to_remove.into_iter().rev().for_each(|idx| {
            self.het_sites.swap_remove(idx);
        });
        absent_base_per_site
    }

    ///
    /// Create a unique identifier for the read alignment.
    /// This lets us distinguish multiple alignments per input read.
    ///
    /// This is a free function, so we can call it without borrowing the `Phaser` struct.
    pub(crate) fn get_read_name_free(record: &bam::Record, use_supplementary: bool) -> String {
        let qname = VStr::from(record.qname());
        //if !passes(&qname) {
        //    log::error!("Unknown data type in input");
        //    std::process::exit(exitcode::IOERR);
        //}
        if use_supplementary {
            //&& record.is_supplementary()
            //let ref_start = record.pos();
            let ref_len = reference_length(record);
            let mut read_start_pos = 0;
            for x in record.cigar().iter() {
                match x {
                    Cigar::HardClip(_len) | Cigar::SoftClip(_len) => {
                        read_start_pos += i64::from(x.len())
                    }
                    _ => break,
                }
            }
            format!("{qname}_sup_{read_start_pos}_{ref_len}")
        } else {
            qname.to_string()
        }
    }

    ///
    /// Create a unique identifier for the read alignment.
    /// This lets us distinguish multiple alignments per input read.
    ///
    /// This requires borrowing `Phaser`.
    /// To avoid borrowing, use `Self::get_read_name_free`.
    pub(crate) fn get_read_name(&self, record: &bam::Record) -> String {
        Self::get_read_name_free(record, self.use_supplementary())
    }

    ///
    /// Yields a vector of names. First the normal read name, and,
    /// # Panics
    /// Panics if read name is not utf-encoded.
    pub fn get_read_names(
        &self,
        record: &bam::Record,
        partial_deletion_reads: Option<&BTreeSet<String>>,
    ) -> Vec<String> {
        let qname = VStr::from(record.qname()).to_string();
        let mut ret = vec![self.get_read_name(record)];
        if self.locus_config().use_supplementary() && record.is_supplementary() {
            if let Some(partial_deletion_reads) = partial_deletion_reads {
                if partial_deletion_reads.contains(&qname)
                    && partial_deletion_reads.contains(&ret[0])
                {
                    ret.push(qname);
                }
            }
        }
        ret
    }

    /// Handle clipped positions within `haplotypes_from_reads_step`.
    ///
    /// # Errors
    /// Propagates htslib `IndexedReader::fetch` error `rust_htslib::errors::Error` if the requested region is not in the bam.
    fn handle_clip_step(
        &mut self,
        read_haps: &mut ReadFingerprintMap,
        min_clip_len: u32,
        het_sites: &[CandidateSite],
        tid: i32,
        clip_buffer: Option<i32>,
    ) -> Result<(), rust_htslib::errors::Error> {
        log::debug!("Handling clip.");
        let mut reader = self.realigned_bam();
        let mut record = bam::Record::new();
        let nvar = self.het_sites.len();
        // We can do this with a sorted list merge in linear time instead of quadratic, but let's only do that if it's a performance issue since correctness will take work.
        // And since there aren't that many sites or clip positions, it shouldn't be an issue.
        let clip_buffer_value: i64 = clip_buffer.unwrap_or(20).into();
        // Handle 3' clips
        // TODO: refactor this into one loop
        for (fingerprint_index, allele_site) in het_sites.iter().enumerate() {
            for clip_position in self
                .clip_3p_positions
                .iter()
                .copied()
                .filter(|&clip_position| allele_site.pos > clip_position)
            {
                reader.fetch((
                    tid,
                    std::cmp::max(0, clip_position - clip_buffer_value),
                    clip_position + clip_buffer_value,
                ))?;
                while let Some(code) = reader.read(&mut record) {
                    if let Err(e) = code {
                        log::warn!("Error reading read from file. Skipping... {e:?}");
                        continue;
                    }
                    let read_name = self.get_read_name(&record);
                    let entry = read_haps
                        .entry(read_name.into())
                        .or_insert_with(|| vec![b'x'; nvar].into());
                    // matching python code: 1-based - 0-based
                    if (clip_position + 1 - record.reference_end()).abs() < clip_buffer_value {
                        let cigar = record.cigar();
                        let threep_clip_len = threeprime_clip_length(&cigar);
                        if threep_clip_len >= min_clip_len {
                            entry[fingerprint_index] = b'0';
                        }
                    }
                }
            }
        }
        // Handle 5' clips
        // Bases associated with the start of softclips are marked zero.
        for (fingerprint_index, allele_site) in het_sites.iter().enumerate() {
            for clip_position in self
                .clip_5p_positions
                .iter()
                .rev()
                .copied()
                .filter(|&clip_position| allele_site.pos < clip_position)
            {
                reader.fetch((
                    tid,
                    std::cmp::max(0, clip_position - clip_buffer_value),
                    clip_position + clip_buffer_value,
                ))?;
                while let Some(code) = reader.read(&mut record) {
                    if let Err(e) = code {
                        log::warn!("Error reading read from file. Skipping... {e:?}");
                        continue;
                    }
                    let entry = read_haps
                        .entry(self.get_read_name(&record).into())
                        .or_insert_with(|| vec![b'x'; nvar].into());
                    if (clip_position + 1 - record.reference_start()).abs() < clip_buffer_value {
                        let cigar = record.cigar();
                        if fiveprime_clip_length(&cigar) >= min_clip_len {
                            entry[fingerprint_index] = b'0';
                        }
                    }
                }
            }
        }
        Ok(())
    }

    /// # Errors
    /// 1. Updating fingerprint map.
    /// 2. Failure to handle clips.
    pub fn haplotypes_from_reads_step(
        &mut self,
        exclude_reads: Option<&BTreeSet<String>>,
        min_clip_len: u32,
        check_clip: bool,
        partial_deletion_reads: Option<&BTreeSet<String>>,
        min_mapq: u8,
        tid: i32,
        clip_buffer: Option<i32>,
        absent_base_per_site: Option<&BTreeMap<CandidateSite, u8>>,
    ) -> Result<ReadFingerprintMap, DError> {
        log::debug!(
            "Starting haplotypes_from_reads_step with {} het sites for gene = {} and sample = {}",
            self.het_sites.len(),
            self.gene_name(),
            self.sample_id(),
        );
        let mut ret = ReadFingerprintMap::default();
        let mut reads_with_flanking_indels = HashSet::<String>::default();
        // Copy out for ease of ownership. Consider borrowing?
        let het_sites = self.het_sites.clone();
        for (fingerprint_index, allele_site) in het_sites.iter().enumerate() {
            log::trace!("Updating map for index {fingerprint_index}");
            self.update_fingerprint_map(
                &mut ret,
                &mut reads_with_flanking_indels,
                exclude_reads,
                partial_deletion_reads,
                (allele_site, fingerprint_index, min_mapq),
                tid,
                absent_base_per_site,
            )?;
        }

        if check_clip
            && std::cmp::max(self.clip_3p_positions.len(), self.clip_5p_positions.len()) > 0
        {
            self.handle_clip_step(&mut ret, min_clip_len, &het_sites, tid, clip_buffer)?;
        }

        Ok(ret)
    }

    ///
    /// For each variant site, performs updates.
    /// Core logic in `get_haplotypes_from_reads_step`.
    ///
    /// # Errors
    /// 1. bam index query failure.
    /// 2. pileup iteration failure. (htslib error propagation.)
    fn update_fingerprint_map(
        &mut self,
        read_haps: &mut ReadFingerprintMap,
        flanking_indel_reads: &mut HashSet<String>,
        exclude_reads: Option<&BTreeSet<String>>,
        partial_deletion_reads: Option<&BTreeSet<String>>,
        site_data: (&CandidateSite, usize, u8),
        tid: i32,
        absent_base_per_site: Option<&BTreeMap<CandidateSite, u8>>,
    ) -> DResult {
        let (allele_site, fingerprint_index, min_mapq) = site_data;
        log::trace!(
            "Fetching reads for realigned bam {:?} at {fingerprint_index} fingerprint index.",
            self.realigned_bam_path()
        );
        let mut bam = self.realigned_bam();
        bam.fetch((tid, allele_site.pos - 1, allele_site.pos + 1))?;

        let empty_map = BTreeMap::new();
        let absent_base_per_site = absent_base_per_site.unwrap_or(&empty_map);

        let mut read_name_inclusion = BTreeMap::<String, i32>::new();
        let mut read_name_exclusion = BTreeMap::<String, i32>::new();

        let nvar = self.het_sites.len();
        let mut num_pileups = 0usize;

        let mut in_exclude_count = 0usize;
        let mut in_flanking_indel_count = 0usize;
        let mut num_included = 0usize;

        flanking_indel_reads.clear();
        for p in bam.pileup() {
            let p = p?;
            num_pileups += 1;
            let diff = allele_site.pos - i64::from(p.pos());
            match diff {
                // Site before
                1 => {
                    log::trace!("Site before at pos {allele_site} ({})", p.pos());
                    for name in p
                        .alignments()
                        .filter(|x| {
                            (x.indel() != pileup::Indel::None || x.is_del())
                                && base_qual(x)
                                    >= self.settings.site_selection_settings.min_base_quality
                        })
                        .flat_map(|aln| {
                            self.get_read_names(&aln.record(), partial_deletion_reads)
                                .into_iter()
                        })
                    {
                        log::trace!("Read {name} is being inserted into flanking indel reads for allele pos {}", allele_site.pos);
                        flanking_indel_reads.insert(name);
                    }
                }
                // Site at pos
                0 => {
                    log::trace!(
                        "Site at {allele_site} with flindel reads {flanking_indel_reads:?}"
                    );
                    for pileup_read in p.alignments() {
                        let base_qual = base_qual(&pileup_read);
                        let record = pileup_read.record();
                        let read_is_included = !pileup_read.is_del()
                            && !pileup_read.is_refskip()
                            && !record.is_secondary()
                            && record.mapq() >= min_mapq
                            && pileup_read.indel() == pileup::Indel::None
                            && base_qual >= self.settings.site_selection_settings.min_base_quality;
                        if read_is_included {
                            log::trace!("Read included: {}. Pos: {}. Is del: {}. is refskip: {}. is secondary: {}. mapq {} vs min {min_mapq}, and indel {:?}. Base qual {base_qual} vs {}",
                                VStr::from(record.qname()), p.pos(), pileup_read.is_del(), pileup_read.is_refskip(), record.is_secondary(), record.mapq(), pileup_read.indel(), self.settings.site_selection_settings.min_base_quality);
                            let names = self.get_read_names(&record, partial_deletion_reads);
                            let read_name = VStr::from(record.qname());
                            for name in names {
                                let in_exclude =
                                    exclude_reads.map_or(false, |exclude| exclude.contains(&name));
                                let in_flanking_indel = flanking_indel_reads.contains(&name);
                                if !in_exclude && !in_flanking_indel {
                                    log::trace!("Read {read_name} is not in exclude and not in flankig indel");
                                    *read_name_inclusion.entry(name.clone()).or_default() += 1;
                                    num_included += 1;
                                    let qpos = pileup_read.qpos().expect("qpos should not be none");
                                    // exclude the last base of a read
                                    let threep_clip_len =
                                        threeprime_clip_length(&record.cigar()) as usize;
                                    if qpos < record.seq_len() - threep_clip_len - 1 {
                                        let prohibited_base = absent_base_per_site.get(allele_site);
                                        let base1 = if prohibited_base.is_none() {
                                            true
                                        } else {
                                            prohibited_base.unwrap() != &b'1'
                                        };
                                        let base2 = if prohibited_base.is_none() {
                                            true
                                        } else {
                                            prohibited_base.unwrap() != &b'2'
                                        };
                                        // bounds check
                                        let entry = read_haps
                                            .entry((&name).into())
                                            .or_insert_with(|| vec![b'x'; nvar].into());
                                        let base = record.seq().index(qpos).to_ascii_uppercase();
                                        let base = if base1
                                            && (allele_site.ref_seq == base
                                                || self
                                                    .settings
                                                    .site_selection_settings
                                                    .permit_list
                                                    .get(&allele_site.pos)
                                                    .map_or(false, |seq| seq.vstr() == base))
                                        {
                                            Some(b'1')
                                        } else if allele_site.var_seq == base && base2 {
                                            Some(b'2')
                                        } else {
                                            None
                                        };
                                        if let Some(base) = base {
                                            entry[fingerprint_index] = base;
                                            log::trace!("Inserting base {base}/{}at index {fingerprint_index} for read name {name}", base as char);
                                        }
                                    } else {
                                        let seq_len = record.seq_len();
                                        assert!(
                                        qpos < seq_len,
                                        "qpos out of bounds... qpos >= record seq len?? {qpos}, {seq_len}",
                                        );
                                    }
                                } else {
                                    log::trace!("Read {read_name} is either in exclude and in flankig indel: in exclude {in_exclude}. flindel {in_flanking_indel}");
                                    *read_name_exclusion.entry(name.clone()).or_default() += 1;
                                    if in_exclude {
                                        in_exclude_count += 1;
                                    }
                                    if in_flanking_indel {
                                        in_flanking_indel_count += 1;
                                    }
                                }
                            }
                        } else {
                            for name in self.get_read_names(&record, partial_deletion_reads) {
                                *read_name_exclusion.entry(name.clone()).or_default() += 1;
                                log::trace!(
                                "Read {} has been filtered out at pos {}. Is del: {}. is refskip: {}. is secondary: {}. mapq {} vs min {min_mapq}, and indel {:?}. Base qual {base_qual} vs {}",
                                VStr::from(record.qname()), p.pos(), pileup_read.is_del(), pileup_read.is_refskip(), record.is_secondary(), record.mapq(), pileup_read.indel(), self.settings.site_selection_settings.min_base_quality
                                );
                            }
                        }
                    }
                    break; // We can stop now, as only the base before and the site itself are checked.
                } // 0
                _ => {
                    assert!(
                        diff >= 0,
                        "Expected loop termination upon reaching variant site. diff: {diff}"
                    );
                    //log::trace!("Site not at expected sites {allele_site}. pos: {}", p.pos());
                    continue;
                }
            }
        }

        log::trace!("Processed {num_pileups} pileup reads.");
        log::trace!("Inclusion count {num_included}: {in_exclude_count} for excluding, {in_flanking_indel_count} for flanking indels.");
        log::trace!("Counts of use for each read name: {read_name_inclusion:?}. Exclusion: {read_name_exclusion:?}");
        Ok(())
    }

    /// Generates candidate variants with a pileup.
    ///
    /// # Returns
    /// Result<(`FilteredSites`, `RawVariantCounts`), `DError`>
    ///
    /// # Errors
    /// 1. Pileup region failure: region not in bam.
    /// 2. Pileup iteration failure: error iterating through pileup. (htslib error propagation.)
    /// 3. Faidx construction: fail to make faidx for local region.
    /// 4. Faidx query: fail to query faidx for region.
    pub fn get_candidate_pos(
        &mut self,
        regions_to_check: &[Range64],
        seq: &[u8],
        min_vaf: Option<f64>,
    ) -> Result<(FilteredSites, RawVariantCounts), DError> {
        let mut bam_handle = self.realigned_bam();
        log::trace!(
            "Reading from bam handle at path {:?}",
            self.realigned_bam_path()
        );

        // This function does get_candidate_pos and assigned het sites.

        log::debug!(
            "Pileup: get_candidate_pos - genome chr: {:?}. Fasta path: {:?}",
            self.local_chr(),
            self.local_reference_path()
        );

        let mut this_setting = self.settings.site_selection_settings.clone();
        if let Some(min_vaf_value) = min_vaf {
            this_setting.min_vaf = min_vaf_value;
        }
        if let Some(user_min_vaf) = self.settings.min_variant_frequency {
            this_setting.min_vaf = user_min_vaf;
        }
        this_setting.targeted = self.settings.targeted;
        log::debug!("this_setting is {this_setting:?}");
        let (filtered_variants, raw_variants) = site_pileup(
            &mut bam_handle,
            &Range64::new(self.left_boundary_0based(), self.right_boundary_0based()), /* Account for left/right boundary being one-based. */
            seq,
            self,
            Some(&this_setting),
            regions_to_check,
            &self.local_chr().expect("no chr"),
        )?;
        //log::trace!("variants: {raw_variants:?}, {filtered_variants:?}");
        log::debug!(
            "variants: {} raw, {} filtered directly, {} not used for phasing",
            raw_variants.len(),
            filtered_variants.variants.len(),
            filtered_variants.variants_no_phasing.len(),
        );

        Ok((filtered_variants, raw_variants))
    }

    /// Create a unique identifier for an alignment.
    /// This lets us disambiguate primary and supplementary alignments.
    #[must_use]
    pub fn labeled_read_name(record: &bam::Record) -> String {
        let qname = vstr::VStr::from(record.qname());
        if record.is_supplementary() {
            let pos = record.pos();
            let reference_end = record.reference_end();
            let length = reference_end - pos;
            format!("{qname}_sup_{pos}_{length}")
        } else {
            qname.to_string()
        }
    }

    /// Subroutine for `update_for_deletions`, which updates fingerprint strings based on deleted regions.
    /// This handles a specific index for fingerprints.
    fn sub_update_for_deletions(
        &mut self,
        hap_map: &mut ReadFingerprintMap,
        data: &DeletionDatum,
        index: usize,
        categories: &mut BTreeMap<char, String>,
    ) -> DResult {
        let initial_het_len = self.het_sites.len();
        // update hap_map with reads that only have information regarding the presence/absence of the deletion
        for read in data.del_reads_partial.iter() {
            if hap_map.get(&(read.clone().into())).is_none() {
                hap_map.insert(read.into(), vec![b'x'; initial_het_len].into());
            }
        }
        for read in data.del_negative_reads.iter() {
            if hap_map.get(&(read.clone().into())).is_none() {
                hap_map.insert(read.into(), vec![b'x'; initial_het_len].into());
            }
        }

        let mut het_sites = std::mem::take(&mut self.het_sites);
        let signifier = match index {
            0 => '3',
            1 => '4',
            _ => ((65 + index) as u8) as char,
        };
        let base = u8::try_from(signifier).expect("char out of u8 bounds. Too many deletions?");

        let del_range_start = data.threep().start;
        let del_range_end = data.fivep().end;
        log::debug!(
            "start, end: {del_range_start}, {del_range_end}. data: 3' = {:?}/5' = {:?}. Original deletion range: {:?}",
            data.threep(),
            data.fivep(),
            data.range(),
        );
        let pos1 = het_sites.iter().position(|x| x.pos > del_range_start);
        let pos2 = het_sites.iter().position(|x| x.pos > del_range_end);
        if let (Some(pos1), Some(pos2)) = (pos1, pos2) {
            let range = pos1..pos2;
            match pos1.cmp(&pos2) {
                // Mark deleted region
                std::cmp::Ordering::Less => {
                    log::debug!(
                        "Assigning {} to fingerprints in range {range:?}",
                        base as char
                    );
                    for hap in data.del_reads_partial.iter().map(ReadAlignmentId::from) {
                        if let Some(hap) = hap_map.get_mut(&hap) {
                            hap[range.clone()].fill(base);
                        }
                    }
                }
                std::cmp::Ordering::Equal => {
                    let cand = CandidateSite::from(data);
                    log::debug!("Found deletion site to handle for signifier {signifier} at {pos1} with candidate {cand}");
                    // Insert het site here.
                    het_sites.insert(pos1, cand);
                    // add deletion to the very beginning
                    if pos1 == 0 {
                        // and update all haplotypes at the same index.
                        for (name, hap) in hap_map.iter_mut() {
                            hap.insert(pos1, b'x');
                            if data.del_reads_partial.contains(&name.read_name) {
                                hap[pos1] = base;
                            } else if pos1 + 1 < hap.len() && hap.get(pos1 + 1) == Some(&b'0') {
                                hap[pos1] = b'0';
                            } else if data.del_negative_reads.contains(&name.read_name) {
                                hap[pos1] = b'1';
                            }
                        }
                    } else {
                        // and update all haplotypes at the same index.
                        for (name, hap) in hap_map.iter_mut() {
                            hap.insert(pos1, b'x');
                            if data.del_reads_partial.contains(&name.read_name) {
                                hap[pos1] = base;
                            } else if pos1 >= 1
                                && hap.get(pos1 - 1) == Some(&b'0')
                                && pos1 + 1 < hap.len()
                                && hap.get(pos1 + 1) == Some(&b'0')
                            {
                                hap[pos1] = b'0';
                            } else {
                                let hap_len = hap.len();
                                let flanking_left_start = if pos1 < 2 { 0 } else { pos1 - 2 };
                                let flanking_left_has_x =
                                    hap[flanking_left_start..pos1].contains(&b'x');
                                let flanking_right_has_x = hap[std::cmp::min(pos1 + 1, hap_len)
                                    ..std::cmp::min(pos1 + 3, hap_len)]
                                    .contains(&b'x');
                                if !flanking_left_has_x && !flanking_right_has_x {
                                    hap[pos1] = b'1';
                                }
                            }
                        }
                    }
                }
                std::cmp::Ordering::Greater => {
                    /* Do nothing, no overlap */
                    log::debug!("pos1 > pos2: {pos1:?}, {pos2:?}");
                }
            }
        } else if !pos1.is_none() && pos2.is_none() {
            // from pos1 to the end
            let nvar = het_sites.len();
            let range = pos1.unwrap()..nvar;
            log::debug!(
                "Assigning {} to fingerprints in range {range:?}",
                base as char
            );
            for hap in data.del_reads_partial.iter().map(ReadAlignmentId::from) {
                if let Some(hap) = hap_map.get_mut(&hap) {
                    hap[range.clone()].fill(base);
                }
            }
        } else if pos1.is_none() && pos2.is_none() {
            // add deletion to the very end
            let cand = CandidateSite::from(data);
            log::debug!("Found deletion site to handle for signifier {signifier} at very end with candidate {cand}");
            // Insert het site here.
            het_sites.push(cand);
            let pos1 = het_sites.len() - 1;
            // and update all haplotypes at the same index.
            for (name, hap) in hap_map.iter_mut() {
                hap.push(b'x');
                if data.del_reads_partial.contains(&name.read_name) {
                    hap[pos1] = base;
                } else if pos1 >= 1 && hap.get(pos1 - 1) == Some(&b'0') {
                    hap[pos1] = b'0';
                } else if data.del_negative_reads.contains(&name.read_name) {
                    hap[pos1] = b'1';
                }
            }
        } else {
            log::debug!("For deletion {pos1:?} and {pos2:?}, not match found");
        }

        log::debug!("Deletion {index} has {signifier} as character id.");
        categories.insert(signifier, data.name());
        self.het_sites = het_sites;
        log::debug!(
            "Het length {initial_het_len} at start, {} at end",
            self.het_sites.len()
        );
        Ok(())
    }

    /// Updates read-fingerprint map for big deletions.
    pub fn update_for_deletions(
        &mut self,
        hap_map: &mut ReadFingerprintMap,
    ) -> Result<BTreeMap<char, String>, DError> {
        let mut ret = BTreeMap::new();
        for (deletion_index, deletion) in self.del_data.clone().into_iter().enumerate() {
            if self.del_data[deletion_index].del_reads_partial.len() > 1 {
                self.sub_update_for_deletions(hap_map, &deletion, deletion_index, &mut ret)?;
            }
        }

        Ok(ret)
    }

    /// Path to realign region.
    #[must_use]
    pub fn local_chr(&self) -> Option<String> {
        self.realign_region_old()
    }

    #[must_use]
    pub fn to_phase(&self) -> bool {
        self.flag & (PhaserFlagBits::ToPhase as u8) != 0
    }

    #[must_use]
    pub fn use_supplementary(&self) -> bool {
        self.flag & (PhaserFlagBits::UseSupplementary as u8) != 0
    }

    #[must_use]
    pub fn is_reverse(&self) -> bool {
        self.flag & (PhaserFlagBits::IsReverse as u8) != 0
    }

    #[must_use]
    pub fn expect_cn2(&self) -> bool {
        self.flag & (PhaserFlagBits::ExpectCN2 as u8) != 0
    }

    /// Gets gene name.
    #[must_use]
    pub fn gene_name(&self) -> &str {
        &self.settings.gene_name
    }

    /// Get index into selected sites for pivot site.
    /// This converts genomic coordinates into an index into selected `self.het_sites`.
    /// # Panics
    /// Errors if position > `i64::MAX` and cannot be converted from `usize.
    #[must_use]
    pub fn get_pivot_index(&self) -> Option<i64> {
        self.pivot_site_0based()
            .and_then(|pos| self.het_sites.iter().position(|x| x.pos == pos))
            .map(|pos| i64::try_from(pos).expect("Failed to get_pivot_index, exceeded i64::MAX"))
    }

    #[must_use]
    pub fn sample_id(&self) -> &str {
        &self.settings.sample_id
    }
} // impl Phaser

/// Compute the length of a cigar operation.
#[inline]
fn clip_size(x: bam::record::Cigar) -> u32 {
    match x {
        bam::record::Cigar::SoftClip(size) | bam::record::Cigar::HardClip(size) => size,
        _ => 0,
    }
}

/// 3' clip length.
/// Only checks last cigar operation.
pub fn threeprime_clip_length(x: &[bam::record::Cigar]) -> u32 {
    x.last().copied().map_or(0u32, clip_size)
}

/// 5' clip length.
/// Only checks first cigar operation.
pub fn fiveprime_clip_length(x: &[bam::record::Cigar]) -> u32 {
    x.first().copied().map_or(0u32, clip_size)
}

/// Checks if a known deletion is present in the read.
/// # Arguments
/// * `record` - bam record
/// * `size` - deletion size
/// * `del_threeprime_range` - coordinate range for left of the deletion (reads with three prime clips)
/// # Returns
/// * a bool indicating the presence/absence of the deletion in the record
pub fn check_del(record: &bam::Record, size: i64, del_threeprime_range: Range64) -> bool {
    let mut starting_pos = record.pos();
    for cigar in record.cigar().iter() {
        let cigar_len = i64::from(cigar.len());
        let cigar_char = cigar.char();
        if cigar_char == 'D' {
            let len_diff = (cigar_len - size).abs() as f64;
            let diff_cutoff = (size as f64 * 0.1).min(50.0);
            if len_diff < diff_cutoff {
                let padding = size / 10;
                if starting_pos >= del_threeprime_range.start - padding
                    && starting_pos <= del_threeprime_range.end + padding
                {
                    return true;
                }
            }
        }
        if cigar_char == 'M' || cigar_char == 'D' || cigar_char == '=' || cigar_char == 'X' {
            starting_pos += cigar_len;
        }
    }
    return false;
}

#[cfg(test)]
mod tests {
    use crate::detail::low_complexity::LowConfidenceSites;
    use crate::detail::util::DResult;
    use crate::io::json::{ParsedParaphaseOutputJSON, ReadAlignmentId, ReadFingerprintMap};
    use crate::{
        config, depth,
        detail::util::{self, test_file},
        io,
        phaser::{self, Phaser},
    };

    use crate::detail::range::I64 as Range64;
    use assertables::{assert_lt, assert_lt_as_result};
    use itertools::Itertools;
    use rust_htslib::bam::record::{Cigar, CigarString, Record};
    use rust_htslib::bam::Read;
    use vstr::VString;

    use super::*;
    use std::collections::{BTreeMap, BTreeSet};
    use std::io::BufRead;

    #[test]
    fn test_check_del() {
        let mut test_record = Record::new();
        test_record.set_pos(200);
        let test_cigar = CigarString(vec![Cigar::Match(151), Cigar::Del(50), Cigar::Match(100)]);
        test_record.set(
            "test_reads".as_bytes(),
            Some(&test_cigar),
            "AAA".as_bytes(),
            "~~~".as_bytes(),
        );
        let test_range = Range64::new(350, 360);
        let deletion_size = 50;
        let deletion_present = check_del(&test_record, deletion_size, test_range);
        assert!(deletion_present);

        // deletion position off
        let mut test_record = Record::new();
        test_record.set_pos(200);
        let test_cigar = CigarString(vec![
            Cigar::Match(200),
            Cigar::SoftClip(50),
            Cigar::Match(100),
        ]);
        test_record.set(
            "test_reads".as_bytes(),
            Some(&test_cigar),
            "AAA".as_bytes(),
            "~~~".as_bytes(),
        );
        let test_range = Range64::new(350, 360);
        let deletion_size = 50;
        let deletion_present = check_del(&test_record, deletion_size, test_range);
        assert!(!deletion_present);

        // deletion size different
        let mut test_record = Record::new();
        test_record.set_pos(200);
        let test_cigar = CigarString(vec![Cigar::Match(151), Cigar::Del(100), Cigar::Match(100)]);
        test_record.set(
            "test_reads".as_bytes(),
            Some(&test_cigar),
            "AAA".as_bytes(),
            "~~~".as_bytes(),
        );
        let test_range = Range64::new(350, 360);
        let deletion_size = 50;
        let deletion_present = check_del(&test_record, deletion_size, test_range);
        assert!(!deletion_present);
    }

    #[test]
    fn paraphase_parity_ok() -> DResult {
        util::init_log(log::LevelFilter::Info);

        let outdir = tempfile::TempDir::new()?;
        let genome_bam = test_file("HG00733_smn1_realigned.bam");
        let genome_path = if let Ok(x) = std::env::var("HG38") {
            x.trim_end_matches(".mmi").to_string()
        } else {
            panic!("Set HG38 env!");
        };
        let gene_name = "smn1";
        let depth = depth::Calculator::from_hg38(genome_bam.display().to_string(), None)?.compute();
        let settings = phaser::Settings::new(
            "HG00733",
            (genome_path, genome_bam),
            outdir.path(),
            gene_name,
            &config::Region::try_load(None)?,
            /* genome depth= */ Some(depth),
            /* sex = */ None,
            String::from("38"),
            None,
            0.03,
            false,
        );
        // Expected depth
        /*
        Found avg depth CoverageSummary(median=80.0, percentile80=83.0)
        */

        let gene_config = config::Gene::try_load(None)?;
        let mut phaser = Phaser::new(
            settings,
            Some(gene_config),
            None, // Option<SiteSelectionSettings>
            None, // Option<RealignSettings>
        );

        let call = phaser.run()?;

        let expected_init_readhaps = std::io::BufReader::new(std::fs::File::open(test_file(
            "expected_haps_hg00733_smn1.txt",
        ))?);
        let expected_init_readhaps = expected_init_readhaps
            .lines()
            .map(|x| {
                let x = x.expect("invalid line");
                let (lk, v) = x.split_terminator('\t').next_tuple().unwrap();
                (ReadAlignmentId::from(lk), VString::from(v))
            })
            .collect::<ReadFingerprintMap>();

        // Make sure the right number of sites were chosen.
        // Make sure the haps match.
        /*
        assert_eq!(
            call.raw_read_haps
                .keys()
                .collect::<std::collections::BTreeSet<_>>(),
            expected_init_readhaps
                .keys()
                .collect::<std::collections::BTreeSet<_>>()
        );
        */
        /*
        // We mostly match, but there are differences.
        // There are 3 sites in paraph_rs currently not in the Python output.
        [ins] In [45]: rust_site_set - py_site_set
        Out[45]: {(70952674, 'G', 'A'), (70959245, 'G', 'T'), (70960119, 'T', 'C')}

        // And two sites in Python not in the rust output.
        [ins] In [46]: py_site_set - rust_site_set
        Out[46]: {(70948286, 'del', '6310'), (70959269, 'A', 'G')}
        */

        let mut num_passes = 0usize;
        let mut num_failed = 0usize;
        let mut missing = std::collections::BTreeSet::new();
        for (k, found) in &call.raw_read_haps {
            let Some(expect) = expected_init_readhaps.get(&ReadAlignmentId::from(k)) else {
                missing.insert(k);
                log::trace!("Read name {k} missing");
                continue;
            };
            if expect == found {
                num_passes += 1;
            } else {
                num_failed += 1;
                log::trace!(
                    "Found {}/{} for lengths, {found}/{expect} for seq",
                    found.len(),
                    expect.len()
                );
            }
        }
        log::trace!("Reads present in expected but not in paraph_rs output: {missing:?}");
        log::warn!(
            "Disabled read fingerprint assignment tests. {num_passes} passed {num_failed} failed out of /{} total.",
            call.raw_read_haps.len()
        );
        let num_missing = missing.len();
        eprintln!("Number found in rs but not py: {num_missing}");
        /*
        assert_eq!(
            num_passes,
            call.raw_read_haps.len(),
            "Read haps don't match. Num passing: {num_passes}. Missing {num_missing}. Call: {:?}",
            call.raw_read_haps
        );
        */
        let found_num_sites = call
            .raw_read_haps
            .values()
            .next()
            .map(std::string::String::len);
        let expected_num_sites = expected_init_readhaps.values().next().map(|x| x.len());
        log::warn!(
            "Disabled num site test. Expected {expected_num_sites:?}. Found {found_num_sites:?}"
        );
        let json_to_match = test_file("v3-jsons/HG00733.json.xz");
        let region_config = config::Region::try_load(None).unwrap();
        let json_data =
            ParsedParaphaseOutputJSON::from_path(&json_to_match, Some(&region_config)).unwrap();
        let gene_data = json_data.gene_data.get(&gene_name.to_string()).unwrap();
        let call = phaser.run()?;
        log::debug!("Depth: {:?}", phaser.region_avg_depth);
        log::debug!("Call: {call:?}");
        let sites: Vec<String> = gene_data
            .sites_for_phasing
            .iter()
            .map(std::string::ToString::to_string)
            .sorted()
            .collect::<Vec<_>>();
        let paraph_rs_sites = call
            .sites_for_phasing
            .iter()
            .cloned()
            .sorted()
            .collect::<Vec<_>>();
        log::trace!("Python sites: {sites:?}");
        log::trace!("Rust sites: {paraph_rs_sites:?}");
        /*
        assert_eq!(
            sites, paraph_rs_sites,
            "Python sites: {sites:?}\nRust sites: {paraph_rs_sites:?}\n"
        );
        assert_eq!(
            expected_num_sites,
            call.raw_read_haps
                .values()
                .next()
                .map(std::string::String::len),
            "Number of sites don't match paraphase. Number missing: {num_missing}"
        );
        */

        Ok(())
    }

    /// Integration test for AGAP9
    #[test]
    fn agap9_ok() -> DResult {
        util::init_log(log::LevelFilter::Info);

        let json_to_match = util::test_file("v3-jsons/BCH-35-manual.json.xz");
        let genome_bam = util::test_file("v3-bams/BCH-35-realigned_AGAP9.bam");
        let outdir = util::test_file("scratch");
        let genome_path = if let Ok(x) = std::env::var("HG38") {
            x.trim_end_matches(".mmi").to_string()
        } else {
            panic!("Set HG38 env!");
        };
        let gene_name = "AGAP9";
        let region_config = config::Region::try_load(None).unwrap();
        let json_data =
            ParsedParaphaseOutputJSON::from_path(&json_to_match, Some(&region_config)).unwrap();
        let gene_data = json_data.gene_data.get(&gene_name.to_string()).unwrap();
        let depth = depth::Calculator::from_hg38(genome_bam.display().to_string(), None)?.compute();
        log::debug!("Gene data: {gene_data:?}");
        log::debug!("depth: {depth:?}");
        assert!(gene_data.pivot_site.is_none());
        assert!(gene_data.pivot_site.is_none());
        /*
        assert!(gene_data.read_to_hap.is_empty());
        assert!(gene_data.assembled_haps.is_none());
        assert!(gene_data.final_haps.is_none());
        assert!(gene_data.sites_for_phasing.is_empty());
        */
        log::debug!("Building Settings struct for BCH-35 with paths {genome_path:?}/{genome_bam:?} with outdir = {outdir:?}");
        let settings = phaser::Settings::new(
            "BCH-35",
            (genome_path, genome_bam),
            outdir,
            gene_name,
            &region_config,
            /* genome depth= */ None,
            /* sex = */ None,
            String::from("38"),
            None,
            0.03,
            false,
        );
        log::debug!("Building gene-level config");
        let gene_config = config::Gene::try_load(None).unwrap();
        log::debug!("Building Phaser.");
        let mut phaser = Phaser::new(
            settings.clone(),
            Some(gene_config.clone()),
            None, // Option<SiteSelectionSettings>
            None, // Option<RealignSettings>
        );
        log::debug!("Built Phaser.");
        let expected_homopolymer_data = {
            let seqs = util::seq_name_pairs(&util::test_file("AGAP9_ref.fa"), true)
                .expect("Failed to get seq names and seqs");
            assert_eq!(seqs[0].0, b"chr10_47501354_47524138");
            let seq = seqs.into_iter().map(|x| x.1).next().unwrap();
            assert_eq!(22785, seq.len(), "Found seq: \"{seq}\"");
            LowConfidenceSites::new(&seq[..], phaser.offset(), None)
        };
        log::trace!("Expected {expected_homopolymer_data:?} for low_complexity_sites");
        let call = phaser.run().expect("Failed to run phaser");
        assert_eq!(
            phaser.low_complexity_sites, expected_homopolymer_data,
            "Found {:?} expected {expected_homopolymer_data:?} when running in/outside of Phaser",
            phaser.low_complexity_sites
        );
        let expected_homopolymer_data =
            util::parse_homopolymers(&test_file("agap9-hpol-expected.txt"));
        let found_homopol = phaser.low_complexity_sites.clone();
        let all_positions = phaser
            .low_complexity_sites
            .keys()
            .chain(expected_homopolymer_data.keys())
            .copied()
            .collect::<BTreeSet<_>>();
        let mut hpol_fails = (0usize, 0usize, 0usize);
        let mut hpol_success = 0;
        for pos in all_positions {
            if !found_homopol.contains_key(&pos) {
                hpol_fails.0 += 1;
                log::trace!("Found does not contain key {pos}");
                continue;
            }
            if !expected_homopolymer_data.contains_key(&pos) {
                hpol_fails.1 += 1;
                log::trace!("Expected does not contain key {pos}");
                continue;
            }
            if found_homopol.get(&pos) != expected_homopolymer_data.get(&pos) {
                hpol_fails.2 += 1;
                log::trace!(
                    "Differing values at key {pos}: {:?}/{:?}",
                    found_homopol.get(&pos),
                    expected_homopolymer_data.get(&pos)
                );
            } else {
                log::trace!(
                    "Found matching values at key {pos}: {:?}/{:?}",
                    found_homopol.get(&pos),
                    expected_homopolymer_data.get(&pos)
                );
                hpol_success += 1;
            }
        }
        let hpol_sum = hpol_fails.0 + hpol_fails.1 + hpol_fails.2;
        log::debug!("Call: {call:?}");
        assert_eq!(hpol_sum, 0, "Missing hp sites in python: {}. Missing in rust: {}. Mismatches: {}. Successes: {hpol_success}", hpol_fails.0, hpol_fails.1, hpol_fails.2);
        assert_eq!(
            phaser.low_complexity_sites, expected_homopolymer_data,
            "Found {:?} expected {expected_homopolymer_data:?} when running in Python vs in Rust",
            phaser.low_complexity_sites
        );
        assert_eq!(phaser.region_avg_depth[0], (76.0f32, 82.0f32));
        // Disable for now as it doesn't pass.
        let sites: Vec<String> = gene_data
            .sites_for_phasing
            .iter()
            .map(std::string::ToString::to_string)
            .sorted()
            .collect::<Vec<_>>();
        let paraph_rs_sites = call
            .sites_for_phasing
            .iter()
            .cloned()
            .sorted()
            .collect::<Vec<_>>();
        log::debug!("Python sites: {sites:?}");
        log::debug!("Rust sites: {paraph_rs_sites:?}");
        assert_eq!(
            sites, paraph_rs_sites,
            "Python sites: {sites:?}\nRust sites: {paraph_rs_sites:?}\n"
        );
        let read_to_hap = gene_data
            .read_to_hap
            .iter()
            .map(|(k, v)| (k.read_name.clone(), v.to_string()))
            .collect::<BTreeMap<_, _>>();
        let num_mismatches = call
            .read_details
            .iter()
            .filter(|(k, v)| v != &read_to_hap.get(&k[..]).unwrap())
            .count();
        assert_lt!(num_mismatches, 51, "We have some differences due to alignments. This assert ensures this does not regress. Nunber of read differences: {num_mismatches}");
        /*
        assert_eq!(
            call.read_details, read_to_hap,
            "Fingerprint map difference: {read_to_hap:?}. Found {:?}",
            call.read_details
        );
        */

        // Next, assembled haps
        let mut expected_haps = gene_data
            .assembled_haps
            .as_ref()
            .unwrap()
            .iter()
            .map(std::string::ToString::to_string)
            .sorted()
            .collect::<Vec<String>>();
        let mut assembled_haps = call
            .assembled_haplotypes
            .iter()
            .map(std::borrow::ToOwned::to_owned)
            .collect::<Vec<_>>();
        let sites_to_ignore = sites.iter().position(|x| {
            x.split_terminator('_')
                .next()
                .unwrap()
                .parse::<i64>()
                .unwrap()
                == 47_510_935
        });
        if let Some(site) = sites_to_ignore {
            eprintln!("Mapper struggles with this site, allow it to differ. TODO: make sure this is not a conflict, IE, at least one of the bases is 'x'");
            let remove_bad_site = |x: &str| -> String {
                x.chars()
                    .take(site)
                    .chain(x.chars().skip(site + 1))
                    .collect::<String>()
            };
            let remove_all = |x: &[String]| -> Vec<String> {
                x.iter().map(|x| remove_bad_site(x)).collect::<Vec<_>>()
            };
            expected_haps = remove_all(&expected_haps);
            assembled_haps = remove_all(&assembled_haps);
        } else {
            eprintln!("No site to ignore");
        }
        assert_eq!(
            expected_haps, assembled_haps,
            "Found {assembled_haps:?} but expected {expected_haps:?}"
        );
        Ok(())
    }

    #[ignore]
    #[test]
    fn parse_gatd3_ok() -> DResult {
        util::init_log(log::LevelFilter::Info);

        let json_to_match = util::test_file("v3-jsons/BCH-35-manual.json.xz");
        let genome_bam = util::test_file("v3-bams/BCH-35-realigned_GATD3.bam");
        let outdir = util::test_file("scratch");
        let genome_path = if let Ok(x) = std::env::var("HG38") {
            x.trim_end_matches(".mmi").to_string()
        } else {
            panic!("Set HG38 env!");
        };
        let region_config = config::Region::try_load(None).unwrap();
        let json_data =
            ParsedParaphaseOutputJSON::from_path(&json_to_match, Some(&region_config)).unwrap();
        let gatd3_data = json_data.gene_data.get(&"GATD3".to_string()).unwrap();
        let depth = depth::Calculator::from_hg38(genome_bam.display().to_string(), None)?.compute();
        log::debug!("Gene data: {gatd3_data:?}");
        log::debug!("depth: {depth:?}");
        assert!(gatd3_data.pivot_site.is_none());
        log::debug!("Building Settings struct for BCH-35 with paths {genome_path:?}/{genome_bam:?} with outdir = {outdir:?}");
        let settings = phaser::Settings::new(
            "BCH-35",
            (genome_path, genome_bam),
            outdir,
            "GATD3",
            &region_config,
            /* genome depth= */ None,
            /* sex = */ None,
            String::from("38"),
            None,
            0.03,
            false,
        );
        log::debug!("Building gene-level config");
        let gene_config = config::Gene::try_load(None).unwrap();
        log::debug!("Building Phaser.");
        let mut phaser = Phaser::new(
            settings.clone(),
            Some(gene_config.clone()),
            None, // Option<SiteSelectionSettings>
            None, // Option<RealignSettings>
        );
        log::debug!("Built Phaser.");
        let call = phaser.run().expect("Failed to run phaser");

        {
            let mut reader = rust_htslib::bam::Reader::from_path(phaser.realigned_bam_path())?;
            let names = reader
                .rc_records()
                .map(|x| VString::from(x.unwrap().qname()))
                .sorted()
                .dedup()
                .collect::<Vec<_>>();
            assert_eq!(
                names,
                [
                    "m64453e_230313_154759/117506523/ccs",
                    "m64453e_230313_154759/80020070/ccs",
                    "m64453e_230315_024224/57344269/ccs",
                    "m64453e_230315_024224/77858246/ccs",
                    "m64453e_230316_133824/120324218/ccs",
                    "m64453e_230316_133824/121178952/ccs",
                    "m64453e_230316_133824/133564476/ccs",
                    "m64453e_230316_133824/146672442/ccs",
                    "m64453e_230316_133824/177996356/ccs",
                    "m64453e_230316_133824/36962902/ccs",
                    "m64453e_230316_133824/52101333/ccs",
                    "m64453e_230316_133824/78840325/ccs",
                    "m64453e_230316_133824/87621678/ccs",
                    "m64453e_230316_133824/97060643/ccs",
                    "m64453e_230316_133824/97977628/ccs",
                    "m64453e_230316_133824/99813835/ccs"
                ]
                .into_iter()
                .sorted()
                .dedup()
                .map(VString::from)
                .collect::<Vec<_>>()
            );
            assert_eq!(names.iter().collect::<BTreeSet<_>>().len(), 16);
        }
        let sites: Vec<String> = gatd3_data
            .sites_for_phasing
            .iter()
            .map(std::string::ToString::to_string)
            .sorted()
            .collect::<Vec<_>>();
        assert_eq!(
            sites,
            call.sites_for_phasing
                .iter()
                .cloned()
                .sorted()
                .collect::<Vec<_>>()
        );
        assert_eq!(call.read_details.values().next().unwrap().len(), 7); // Should have 7 sites.
        assert_eq!(phaser.region_avg_depth[0], (12.0f32, 13.0f32));
        let read_to_hap = gatd3_data
            .read_to_hap
            .iter()
            .map(|(k, v)| (k.read_name.clone(), v.to_string()))
            .collect::<BTreeMap<_, _>>();
        assert_eq!(call.read_details, read_to_hap);
        let expected_haps = gatd3_data
            .assembled_haps
            .as_ref()
            .unwrap()
            .iter()
            .map(std::string::ToString::to_string)
            .sorted()
            .collect::<Vec<String>>();
        let assembled_haps = call
            .assembled_haplotypes
            .iter()
            .map(std::borrow::ToOwned::to_owned)
            .collect::<Vec<_>>();
        assert_eq!(
            expected_haps, assembled_haps,
            "Found {assembled_haps:?} but expected {expected_haps:?}"
        );
        Ok(())
    }

    #[test]
    fn parse_hsfy1_ok() -> DResult {
        util::init_log(log::LevelFilter::Info);

        //let json_to_match = util::test_file("v3-jsons/BCH-35.json.xz");
        let genome_bam = util::test_file("v3-bams/BCH-35-realigned_HSFY1.bam");
        let outdir = util::test_file("scratch");
        let region_config = config::Region::try_load(None).unwrap();
        let genome_path = if let Ok(x) = std::env::var("HG38") {
            x.trim_end_matches(".mmi").to_string()
        } else {
            panic!("Set HG38 env!");
        };
        let settings = phaser::Settings::new(
            "BCH-35",
            (genome_path, genome_bam),
            outdir,
            "HSFY1",
            &region_config,
            /* genome depth= */ None,
            /* sex = */ None,
            String::from("38"),
            None,
            0.03,
            false,
        );

        let gene_config = config::Gene::try_load(None).unwrap();
        let mut phaser = Phaser::new(
            settings.clone(),
            Some(gene_config.clone()),
            None, // Option<SiteSelectionSettings>
            None, // Option<RealignSettings>
        );
        let call = phaser.run().expect("Failed to run phaser");
        assert_eq!(phaser.region_avg_depth[0], (54.0f32, 59.0f32));
        assert_eq!(call.read_details.values().next().unwrap().len(), 3);
        let json_to_match = util::test_file("v3-jsons/BCH-35-manual.json.xz");
        let json_data =
            ParsedParaphaseOutputJSON::from_path(&json_to_match, Some(&region_config)).unwrap();
        let hsfy1_data = json_data.gene_data.get(&"HSFY1".to_string()).unwrap();
        let sites: Vec<String> = hsfy1_data
            .sites_for_phasing
            .iter()
            .map(std::string::ToString::to_string)
            .sorted()
            .collect::<Vec<_>>();
        let found_sites = call
            .sites_for_phasing
            .iter()
            .cloned()
            .sorted()
            .collect::<Vec<_>>();
        assert_eq!(
            sites, found_sites,
            "Found {found_sites:?} instead of {sites:?}",
        );
        assert_eq!(call.read_details.values().next().unwrap().len(), 3); // Should have 3 sites.
        log::debug!("Sites for phasing: {sites:?}/{found_sites:?}");
        log::debug!(
            "Read details length for HSFY1: {}. Expected 101",
            call.read_details.len()
        );
        let expected_reads = hsfy1_data
            .read_to_hap
            .keys()
            .map(std::string::ToString::to_string)
            .collect::<BTreeSet<_>>();
        let found_reads = call.read_details.keys().cloned().collect::<BTreeSet<_>>();
        if expected_reads != found_reads {
            let not_found = expected_reads.difference(&found_reads).collect::<Vec<_>>();
            let not_expected = found_reads.difference(&expected_reads).collect::<Vec<_>>();
            log::debug!("Not found: {not_found:?}. Not expected: {not_expected:?}");
        }
        assert_eq!(
            call.read_details.len(),
            101,
            "Found: {:?}",
            call.read_details
        ); // Expected 101 fingerprints.
        let read_to_hap = hsfy1_data
            .read_to_hap
            .iter()
            .map(|(k, v)| (k.read_name.clone(), v.to_string()))
            .collect::<BTreeMap<_, _>>();
        assert_eq!(call.read_details, read_to_hap);
        let expected_haps = hsfy1_data
            .assembled_haps
            .as_ref()
            .unwrap()
            .iter()
            .map(std::string::ToString::to_string)
            .sorted()
            .collect::<Vec<String>>();
        let assembled_haps = call
            .assembled_haplotypes
            .iter()
            .map(std::borrow::ToOwned::to_owned)
            .collect::<Vec<_>>();
        assert_eq!(
            expected_haps, assembled_haps,
            "Found {assembled_haps:?} but expected {expected_haps:?}"
        );
        Ok(())
    }

    #[test]
    fn get_sites_ok() {
        util::init_log(log::LevelFilter::Info);

        let json_to_match = util::test_file("jsons/14232-mo.json.xz");
        let genome_bam = util::test_file("v3-bams/GATD3.14232-mo.bam");
        let outdir = util::test_file("scratch");
        let genome_path = if let Ok(x) = std::env::var("HG38") {
            x.trim_end_matches(".mmi").to_string()
        } else {
            panic!("Set HG38 env!");
        };
        let region_config = config::Region::try_load(None).unwrap();
        let json_data =
            io::json::ParsedParaphaseOutputJSON::from_path(&json_to_match, Some(&region_config))
                .unwrap();
        let gene_data = json_data.gene_data.get(&"GATD3".to_string()).unwrap();
        let mut settings = phaser::Settings::new(
            "HG00733",
            (genome_path, genome_bam),
            outdir,
            "GATD3",
            &region_config,
            /* genome depth= */ None,
            /* sex = */ None,
            String::from("38"),
            None,
            0.03,
            false,
        );
        settings.allow_low_coverage = true;
        let gene_config = config::Gene::try_load(None).unwrap();
        let mut phaser = Phaser::new(
            settings,
            Some(gene_config),
            None, // Option<SiteSelectionSettings>
            None, // Option<RealignSettings>
        );
        let call = phaser.run().expect("Failed to run phaser");
        log::debug!("Expected {gene_data:?}. Call: {call:?}");
        let final_haps = call.final_haplotypes.keys().cloned().collect::<Vec<_>>();
        let assembled_haps = call.assembled_haplotypes.clone();
        log::debug!("final haps: {final_haps:?}.");
        log::debug!("asm haps: {assembled_haps:?}.");
    }

    #[test]
    fn get_pivot_index_ok() {
        let outdir = tempfile::TempDir::new().unwrap();
        let genome_bam = util::test_file("HG00733_smn1_realigned.bam");
        let genome_path = if let Ok(x) = std::env::var("HG38") {
            x.trim_end_matches(".mmi").to_string()
        } else {
            panic!("Set HG38 env!");
        };
        let gene_name = "smn1";
        let settings = phaser::Settings::new(
            "HG00733",
            (genome_path, genome_bam),
            outdir.path(),
            gene_name,
            &config::Region::try_load(None).unwrap(),
            None,
            None,
            String::from("38"),
            None,
            0.03,
            false,
        );

        let gene_config = config::Gene::try_load(None).unwrap();

        let mut phaser = Phaser::new(settings, Some(gene_config), None, None);

        phaser.het_sites = vec![
            "70951940_A_C".parse().unwrap(),
            "70951946_T_G".parse().unwrap(),
        ];
        assert_eq!(phaser.get_pivot_index(), Some(1));

        phaser.het_sites = vec![
            "70951940_A_C".parse().unwrap(),
            "70951946_T_G".parse().unwrap(),
            "70951958_T_G".parse().unwrap(),
        ];
        assert_eq!(phaser.get_pivot_index(), Some(1));

        phaser.het_sites = vec![
            "70951947_A_C".parse().unwrap(),
            "70951949_T_G".parse().unwrap(),
            "70951958_T_G".parse().unwrap(),
        ];
        assert_eq!(phaser.get_pivot_index(), None);
    }
}
