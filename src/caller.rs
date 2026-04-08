use crate::assembly::assembler::{build_graph, GraphParameters};
use crate::bam_operation::{
    get_start_end_d4z4, get_start_end_from_genome, realign, tag_reads, ClippedReads,
};
use crate::cli::Settings;
use crate::d4z4::d4z4_phasing::{
    find_cis_dup, get_background_for_allele_ends, get_background_for_allele_starts,
    haplotype_background, phase_flanking,
};
use crate::d4z4::join_partial_alleles::{join_partial_alleles, AlleleSummary};
use crate::depth::median;
use crate::depth::{depth_based_cn, DepthSummary};
use crate::methylation::{get_methyl_info, methyl_prob_by_position, MethOutput};
use crate::plot::plot_alleles::plot_alleles_and_reads;
use crate::read_filtering::{filter_realignments_d4z4, filter_realignments_kiv2};
use crate::repeat_unit::fingerprint::{get_fingerprint, ReadParameters};
use crate::repeat_unit::fingerprint_utils::rm_redundant_finger_prints;
use crate::util::{d4z4_coordinates, kiv2_coordinates, DError, DResult};
use crate::variant::report_variants;
use crate::vcf::write_vcf;
use log::{debug, info};
use paraphase::io::json::GeneCall;
//use std::cmp;
use std::collections::BTreeMap;
use std::fmt::Display;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::str;

/// Sample call report
#[derive(Clone, Debug, Default, serde::Serialize, serde::Deserialize)]
pub struct SampleCall {
    /// copy number of assembled alleles
    pub allele_cn: String,
    /// summary of depth information
    pub depth_summary: DepthSummary,
    /// completely assembled alleles
    pub complete_alleles: Vec<String>,
    /// partially assembled alleles
    pub partial_alleles: Vec<String>,
    /// supporting reads for complete alleles
    pub supporting_reads: BTreeMap<String, Vec<String>>,
    /// reporting summary information for D4Z4 alleles
    pub allele_info: Vec<AlleleSummary>,
    /// variants on fingerprints on complete alleles
    pub complete_allele_variants: BTreeMap<String, Vec<String>>,
    /// variants on fingerprints not found on complete alleles
    pub other_unit_variants: Vec<String>,
    /// summary of methylation information, including complete alleles and all allele ends
    pub methylation: BTreeMap<String, Option<MethOutput>>,
    /// additional key-value metadata
    pub additional: BTreeMap<String, serde_json::Value>,
}

/// Summary of paraphase phasing information
#[derive(Clone, Debug, Default, serde::Serialize, serde::Deserialize)]
pub struct AlleleFlankingPhasing {
    /// final haplotypes
    pub final_haplotypes: BTreeMap<String, String>,
    /// unique supporting reads
    pub unique_supporting_reads: BTreeMap<String, Vec<String>>,
    /// sites for phasing
    pub sites_for_phasing: Vec<String>,
    /// haplotype details, contains boundaries and variants
    pub haplotype_details: BTreeMap<String, paraphase::detail::phase_haps::HapInfoForJson>,
}

/// D4Z4 QC metrics
#[derive(Clone, Debug, Default, serde::Serialize, serde::Deserialize)]
pub struct D4Z4QCMetrics {
    /// median read length
    pub median_read_length: f32,
    /// per allele depth
    pub per_allele_depth: f32,
}

/// Convert alleles from vectors to strings
/// # Arguments
/// * `haps` - alleles represented as vectors of fingerprints
/// * `separater` - symbol to join fingerprints
/// # Returns
/// * `Vec<String>` - alleles represented as strings
pub fn vec_to_string<T: Display>(haps: &Vec<Vec<T>>, separater: &str) -> Vec<String> {
    let mut haps_string = Vec::new();
    for hap in haps {
        let mut renamed_hap = Vec::new();
        for fp in hap {
            let fp_string = fp.to_string();
            if fp_string.contains(&String::from("-")) {
                let fp_int = fp_string.parse::<i64>().unwrap();
                if fp_int < 0 {
                    if fp_int > -10 {
                        renamed_hap.push(String::from("LeftFlank"));
                    } else if fp_int == -10 {
                        renamed_hap.push(String::from("RightFlank"));
                    } else {
                        renamed_hap.push(String::from("RightFlankB"));
                    }
                }
            } else {
                renamed_hap.push(fp_string);
            }
        }
        let hap_string = renamed_hap.join(separater);
        haps_string.push(hap_string);
    }
    haps_string
}

/// Remove a file if it exists
/// # Arguments
/// * `path` - path to the file
fn remove_if_exists(path: impl AsRef<Path>) -> DResult {
    let path = path.as_ref();
    if path.exists() {
        std::fs::remove_file(path)?;
    }
    Ok(())
}

/// Parse user-provided variant list
/// # Arguments
/// * `variant_file_path` - path to the variant file
/// # Returns
/// * `Option<Vec<String>>` - predefined variants
pub fn get_predefined_variants(
    variant_file_path: Option<PathBuf>,
) -> Result<Option<Vec<String>>, DError> {
    let mut predefined_variants: Option<Vec<String>> = None;
    if let Some(variant_file) = variant_file_path.as_ref() {
        let variants = std::fs::read(variant_file)?;
        let variants = std::str::from_utf8(&variants)?;
        let variants = variants
            .split_terminator('\n')
            .map(std::borrow::ToOwned::to_owned)
            .collect::<Vec<_>>();
        debug!("predefined variants are {:?}", variants.clone());
        predefined_variants = Some(variants);
    }
    Ok(predefined_variants)
}

/// Main function to genotype KIV2
/// # Arguments
/// * `cli_settings` - CLI settings
pub fn call_kiv(cli_settings: Settings) -> DResult {
    let kiv2_read_parameters = ReadParameters {
        min_base_quality: 10,
        min_variant_support: 6,
        min_fingerprint_support: 2,
        max_read_count_to_correct: 4,
    };
    let kiv2_graph_parameters = GraphParameters {
        expect_n_allele: 2,
        min_overlap: 2,
        run_twice: false,
        less_filtering: false,
    };
    let sample_id = &cli_settings.prefix;
    let output_path = cli_settings.outdir.as_path();
    std::fs::create_dir_all(output_path)?;
    // output files
    let realigned_bam = output_path.join(format!("{sample_id}.kivvi.kiv2.bam"));
    let output_json = output_path.join(format!("{sample_id}.kivvi.kiv2.json"));
    let output_vcf = output_path.join(format!("{sample_id}.kivvi.kiv2.vcf"));
    let output_svg = output_path.join(format!("{sample_id}.kivvi.kiv2.svg"));
    // other region specific resources
    let region_coordinates = kiv2_coordinates();
    // create temporary reference file
    let reference = output_path.join(format!("{sample_id}.kiv2.ref.fa"));
    std::fs::write(&reference, region_coordinates.clone().reference_seq)
        .expect("Unable to write temporary reference file");
    // predefined variant list
    let predefined_variants = get_predefined_variants(cli_settings.variant_list)?;

    // realign reads and filter alignments
    debug!("Realign reads to repeat unit...");
    let (realn_records_unfiltered, read_length, writer, _methyl_tags, _methyl_probs) = realign(
        cli_settings.bam_filename.clone(),
        region_coordinates.clone(),
        &reference,
        realigned_bam.clone(),
        false,
    )?;
    let repeat_records = filter_realignments_kiv2(
        realn_records_unfiltered,
        writer,
        &reference,
        realigned_bam.clone(),
    )?;

    // get flanking reads
    debug!("Get flanking reads...");
    let flanking_reads = get_start_end_from_genome(
        cli_settings.bam_filename.clone(),
        region_coordinates.clone(),
        None,
    )?;
    debug!("starting_reads_flank: {:?}", flanking_reads.start);
    debug!("ending_reads_flank: {:?}", flanking_reads.end);

    // get genome depth and repeat depth
    debug!("Get read depth...");
    let depth_summary = depth_based_cn(
        cli_settings.bam_filename.clone(),
        realigned_bam.clone(),
        region_coordinates.clone(),
    )?;

    // get fingerprints
    debug!("Get fingerprints...");
    let (mut fp_info, _bases_at_pivot_site, _cpg_sites_per_read) = get_fingerprint(
        realigned_bam.clone(),
        &reference,
        region_coordinates.clone(),
        flanking_reads.clone(),
        ClippedReads::default(),
        &read_length,
        kiv2_read_parameters.clone(),
        vec![],
        false,
        cli_settings.sensitive,
    )?;

    // remove redundant fingerprints
    debug!("Remove redundant fingerprints...");
    loop {
        let (new_fp_info, changed) = rm_redundant_finger_prints(
            fp_info,
            kiv2_read_parameters.max_read_count_to_correct,
            false,
        )?;
        fp_info = new_fp_info.clone();
        if !changed {
            break;
        }
    }
    let read_info = fp_info.clone().read_bases;

    // tag reads in bam by fingerprints
    debug!("Tag reads with fingerprints...");
    let _ = tag_reads(
        repeat_records,
        fp_info.clone().grouped_reads,
        realigned_bam.clone(),
        region_coordinates.clone(),
        BTreeMap::new(),
        &reference,
    )?;

    // graph assembler
    debug!("Assemble alleles...");
    let mut fp_graph = build_graph(
        fp_info.clone().read_edges,
        kiv2_graph_parameters.min_overlap,
    );
    let assembly_result = fp_graph.run(kiv2_graph_parameters.clone())?;
    debug!(
        "complete_haps {:?} incomplete_haps {:?}",
        assembly_result.complete, assembly_result.incomplete
    );

    // call variants on fingerprints
    debug!("Call variants...");
    let variant_report = report_variants(
        fp_info.clone(),
        assembly_result.clone(),
        read_info.clone(),
        &reference,
        region_coordinates.clone(),
        predefined_variants,
    )?;

    // write to json
    debug!("Write to json...");
    // convert supporting reads to renamed alleles
    let mut support: BTreeMap<String, Vec<String>> = BTreeMap::new();
    for (hap, reads) in assembly_result.supporting_reads.iter() {
        let allele_name = vec_to_string(&vec![hap.to_vec()], "-");
        let hap_string = &allele_name[0];
        support
            .entry(hap_string.to_string())
            .or_insert(reads.iter().cloned().collect::<Vec<String>>());
    }
    let allele_cn = &assembly_result
        .complete
        .iter()
        .map(|x| x.len() - 2)
        .collect::<Vec<usize>>()
        .iter()
        .map(|i| i.to_string())
        .collect::<Vec<_>>()
        .join("/");
    let final_complete = vec_to_string(&assembly_result.complete, "-");
    // pack variants for reporting
    let mut complete_allele_variants = BTreeMap::new();
    for (allele, allele_variants) in variant_report.complete_allele_variants {
        let mut allele_variants_reformat = allele_variants.into_iter().collect::<Vec<_>>();
        allele_variants_reformat.sort_by(|a, b| a.0 .0.cmp(&b.0 .0));
        let allele_variants_reformat = allele_variants_reformat
            .iter()
            .map(|x| format!("{};{};{}", x.0 .0, x.0 .1, x.1.join(",")))
            .collect::<Vec<_>>();
        complete_allele_variants.insert(
            vec_to_string(&vec![allele], "-")
                .first()
                .unwrap()
                .to_string(),
            allele_variants_reformat,
        );
    }
    let sample_call = SampleCall {
        allele_cn: allele_cn.to_string(),
        depth_summary,
        complete_alleles: final_complete,
        partial_alleles: vec_to_string(&assembly_result.incomplete, "-"),
        supporting_reads: support,
        complete_allele_variants,
        methylation: BTreeMap::new(),
        other_unit_variants: variant_report
            .fp_variants_on_incomplete_alleles
            .into_iter()
            .map(|(a, b)| format!("{a};{}", b.join(",")))
            .collect::<Vec<String>>(),
        allele_info: Vec::new(),
        ..Default::default()
    };
    let mut writer = std::io::BufWriter::new(std::fs::File::create(output_json)?);
    writeln!(writer, "{}", serde_json::to_string_pretty(&sample_call)?)?;

    // write to vcf
    debug!("Write to VCF...");
    std::fs::File::create(output_vcf.clone())?;
    write_vcf(
        &output_vcf,
        &sample_id,
        variant_report.variant_summary,
        region_coordinates.clone(),
        assembly_result.complete.len(),
    )?;

    // plot
    debug!("Plot alleles...");
    let alleles_for_plot = variant_report.alleles_for_plot;
    if let Some(to_plot) = alleles_for_plot {
        plot_alleles_and_reads(&output_svg, to_plot)?;
    }

    // remove temporary reference file
    remove_if_exists(reference)?;
    let fai_file = output_path.join(format!("{sample_id}.kiv2.ref.fa.fai"));
    remove_if_exists(fai_file)?;

    info!("Completed kivvi analysis on KIV2...");
    Ok(())
}

/// Main function to genotype D4Z4
/// # Arguments
/// * `cli_settings` - CLI settings
pub fn call_d4z4(cli_settings: Settings) -> DResult {
    let d4z4_read_parameters = ReadParameters {
        min_base_quality: 3,
        min_variant_support: 4,
        min_fingerprint_support: 2,
        max_read_count_to_correct: 3,
    };
    let d4z4_graph_parameters = GraphParameters {
        expect_n_allele: 4,
        min_overlap: 3,
        run_twice: true,
        less_filtering: false,
    };
    let sample_id = &cli_settings.prefix;
    let output_path = cli_settings.outdir.as_path();
    std::fs::create_dir_all(output_path)?;
    // output files
    let realigned_bam = output_path.join(format!("{sample_id}.kivvi.d4z4.bam"));
    let output_json = output_path.join(format!("{sample_id}.kivvi.d4z4.json"));
    let output_vcf = output_path.join(format!("{sample_id}.kivvi.d4z4.vcf"));
    let output_svg = output_path.join(format!("{sample_id}.kivvi.d4z4.svg"));
    //let output_methyl_svg = output_path.join(format!("{sample_id}.kivvi.d4z4.methyl.svg"));
    // other region specific resources
    let region_coordinates = d4z4_coordinates();
    debug!("region_coordinates: {:?}", region_coordinates);
    // create temporary reference file
    let reference = output_path.join(format!("{sample_id}.d4z4.ref.fa"));
    std::fs::write(&reference, region_coordinates.clone().reference_seq)
        .expect("Unable to write temporary reference file");
    // create tempory file for the modified genome reference, specially made for d4z4
    lazy_static::lazy_static! {
        pub static ref GENOME_REFERENCE: &'static [u8] = std::include_bytes!(concat!(env!("CARGO_MANIFEST_DIR"), "/data/d4z4/chr4_mod.fa"));
    }
    let d4z4_genome_reference_seq = str::from_utf8(&GENOME_REFERENCE).unwrap().to_string();
    let genome_reference = output_path.join(format!("{sample_id}.d4z4.genome.fa"));
    std::fs::write(&genome_reference, d4z4_genome_reference_seq)
        .expect("Unable to write temporary genome reference file");

    // predefined variant list
    let predefined_variants = get_predefined_variants(cli_settings.variant_list)?;

    // phase flanking region with paraphase
    let mut phasing_result: BTreeMap<String, GeneCall> = BTreeMap::new();
    if !cli_settings.nopp {
        phasing_result = phase_flanking(
            sample_id,
            output_path,
            &cli_settings.bam_filename,
            &genome_reference,
            cli_settings.verbosity > 0,
        )?;
    } else {
        phasing_result.insert(String::from("DUX4p5"), GeneCall::default());
        phasing_result.insert(String::from("DUX4"), GeneCall::default());
    }
    //debug!("phasing_result {:?}", phasing_result);

    // realign reads and filter alignments
    debug!("Realign reads to repeat unit...");
    let (realn_records_unfiltered, read_length, writer, methyl_tags, methyl_probs) = realign(
        cli_settings.bam_filename.clone(),
        region_coordinates.clone(),
        &reference,
        realigned_bam.clone(),
        true,
    )?;
    let (repeat_records, mut white_list_read_segments) = filter_realignments_d4z4(
        realn_records_unfiltered,
        writer,
        &reference,
        realigned_bam.clone(),
    )?;
    //debug!("white_list_read_segments {:?}", white_list_read_segments);

    // get flanking reads
    debug!("Get flanking reads...");
    let (flanking_reads, clipped_reads) = get_start_end_d4z4(
        realigned_bam.clone(),
        &reference,
        read_length.clone(),
        region_coordinates.clone(),
    )?;
    debug!("flanking reads: {:?}", flanking_reads);
    // give green light to fully spanning reads
    for read_segment in &flanking_reads.start_segment {
        let read_name = read_segment
            .split(':')
            .next()
            .ok_or("next not in read_segment")?
            .to_string();
        if flanking_reads.end.contains(&read_name) {
            white_list_read_segments.push(read_name.to_string());
        }
    }

    debug!("white_list_read_segments {:?}", white_list_read_segments);

    // get genome depth and repeat depth
    debug!("Get read depth...");
    let depth_summary = depth_based_cn(
        cli_settings.bam_filename.clone(),
        realigned_bam.clone(),
        region_coordinates.clone(),
    )?;

    // get fingerprints
    debug!("Get fingerprints...");
    let (mut fp_info, bases_at_pivot_site, cpg_sites_per_read) = get_fingerprint(
        realigned_bam.clone(),
        &reference,
        region_coordinates.clone(),
        flanking_reads.clone(),
        clipped_reads,
        &read_length,
        d4z4_read_parameters.clone(),
        white_list_read_segments,
        true,
        cli_settings.sensitive,
    )?;

    // remove redundant fingerprints
    debug!("Remove redundant fingerprints...");
    if true {
        loop {
            let (new_fp_info, changed) = rm_redundant_finger_prints(
                fp_info,
                d4z4_read_parameters.max_read_count_to_correct,
                true,
            )?;
            fp_info = new_fp_info.clone();
            if !changed {
                break;
            }
        }
    }

    // qc metrics
    let all_read_length = read_length
        .iter()
        .filter(|(x, _y)| fp_info.read_edges.contains_key(*x))
        .map(|(_x, y)| *y as i32)
        .collect::<Vec<i32>>();
    let median_read_length = median(&all_read_length).unwrap();
    let mut starting_reads_count = 0.0;
    for (_read, read_edges) in &fp_info.read_edges {
        if read_edges.iter().any(|x| *x < 0 && *x > -10) {
            starting_reads_count += 1.0;
        }
    }
    let haploid_depth = starting_reads_count / 4.0;

    let read_info = fp_info.clone().read_bases;
    // tag reads in bam by fingerprints
    debug!("Tag reads with fingerprints...");
    let _tag_success = tag_reads(
        repeat_records,
        fp_info.clone().grouped_reads,
        realigned_bam.clone(),
        region_coordinates.clone(),
        methyl_tags,
        &reference,
    )?;

    // graph assembler
    debug!("Assemble alleles...");
    let mut fp_graph = build_graph(
        fp_info.clone().read_edges,
        d4z4_graph_parameters.min_overlap,
    );
    let assembly_result = fp_graph.run(d4z4_graph_parameters.clone())?;
    debug!(
        "complete_haps {:?} incomplete_haps {:?}",
        assembly_result.complete, assembly_result.incomplete
    );

    // methylation info
    let (meth_summary, segment_methyl_prob) =
        methyl_prob_by_position(&methyl_probs, &cpg_sites_per_read, &fp_info)?;

    // find in-cis duplications
    debug!("Find in-cis duplications...");
    let (cis_dups, cis_dups_match_index) =
        find_cis_dup(&assembly_result, &fp_graph, &fp_info, &phasing_result)?;
    debug!("cis_dups_match_index {cis_dups_match_index:?}");

    // get all starting haps
    let (all_starts_hap_backgrounds, all_starts_upstream_haplotypes) =
        get_background_for_allele_starts(
            &assembly_result,
            &fp_graph,
            &phasing_result,
            &bases_at_pivot_site,
            &fp_info,
        )?;

    // get all ending haps
    let (all_ends_hap_backgrounds, all_ends_reads_match_allele_index_renamed_to_string) =
        get_background_for_allele_ends(
            &assembly_result,
            &fp_graph,
            &phasing_result,
            &bases_at_pivot_site,
            &fp_info,
            &cis_dups_match_index,
        )?;
    // all ends methylation
    let mut all_ends_allele_methyl: Option<MethOutput> = None;
    if !methyl_probs.is_empty() {
        let (all_ends_meth_out, _all_ends_reads_methyl_value) = get_methyl_info(
            &fp_info,
            &segment_methyl_prob,
            &cpg_sites_per_read,
            &all_ends_reads_match_allele_index_renamed_to_string,
            &region_coordinates.methyl_sites,
        )?;
        all_ends_allele_methyl = Some(all_ends_meth_out);
    }

    // call variants on fingerprints
    debug!("Call variants...");
    let variant_report = report_variants(
        fp_info.clone(),
        assembly_result.clone(),
        read_info.clone(),
        &reference,
        region_coordinates.clone(),
        predefined_variants,
    )?;

    // methylation on fully assembled alleles
    let mut allele_methyl: Option<MethOutput> = None;
    if !methyl_probs.is_empty() {
        if variant_report.fp_suppporting_reads.is_some() {
            let (meth_out, _alleles_reads_methyl_value) = get_methyl_info(
                &fp_info,
                &segment_methyl_prob,
                &cpg_sites_per_read,
                &variant_report.reads_match_allele_index,
                &region_coordinates.methyl_sites,
            )?;
            //plot_methyl(&output_methyl_svg, alleles_reads_methyl_value)?;
            allele_methyl = Some(meth_out);
        }
    }

    // haplotype backgrounds
    let complete_haps = assembly_result.complete.clone();
    let complete_read_support = fp_graph
        .process_complete_haps(complete_haps, None, false, false)?
        .supporting_reads;
    let (mut hap_backgrounds, _upstream_haplotypes) = haplotype_background(
        &assembly_result.complete,
        &phasing_result,
        &complete_read_support,
        Some(&fp_info),
        &variant_report.reads_match_allele_index,
        &bases_at_pivot_site,
        true,
        true,
    )?;
    let mut updated_hap_backgrounds = BTreeMap::new();
    for (hap, hap_background) in &hap_backgrounds {
        let hap_background_parts = hap_background.split("-").collect::<Vec<&str>>();
        let mut distal_hap = hap_background_parts[0].to_string();
        let mut chr = hap_background_parts[1].to_string();
        if distal_hap.contains("unknown") {
            if all_ends_hap_backgrounds.contains_key(hap) {
                distal_hap = all_ends_hap_backgrounds.get(hap).unwrap().to_string();
            }
        }
        if chr.contains("unknown") {
            if all_starts_hap_backgrounds.contains_key(hap) {
                chr = all_starts_hap_backgrounds.get(hap).unwrap().to_string();
            }
        }
        updated_hap_backgrounds.insert(hap.clone(), format!("{distal_hap}-{chr}"));
    }
    hap_backgrounds = updated_hap_backgrounds;

    // join partial alleles
    debug!("Join partial alleles...");
    let allele_summaries = join_partial_alleles(
        &all_starts_hap_backgrounds,
        &all_ends_hap_backgrounds,
        &hap_backgrounds,
        &variant_report,
        &region_coordinates,
        &fp_graph,
        &fp_info,
        &all_ends_allele_methyl,
    )?;

    // write to json
    debug!("Write to json...");
    // convert supporting reads to renamed alleles
    let mut support: BTreeMap<String, Vec<String>> = BTreeMap::new();
    for allele in &assembly_result.complete {
        if assembly_result.supporting_reads.contains_key(allele) {
            let reads = assembly_result
                .supporting_reads
                .get(allele)
                .ok_or("hap not in assembly_result.supporting_reads.")?;
            let allele_name = vec_to_string(&vec![allele.clone()], "-");
            let hap_string = &allele_name[0];
            support
                .entry(hap_string.to_string())
                .or_insert(reads.iter().cloned().collect::<Vec<String>>());
        }
    }

    let allele_cn = &assembly_result
        .complete
        .iter()
        .map(|x| x.len() - 2)
        .collect::<Vec<usize>>()
        .iter()
        .map(|i| i.to_string())
        .collect::<Vec<_>>()
        .join(",");
    let final_complete = vec_to_string(&assembly_result.complete, "-");

    // summarize upstream/downstream haplotype backgrounds
    let upstream_phasing_result = phasing_result.get(&String::from("DUX4p5")).unwrap();
    let downstream_phasing_result = phasing_result.get(&String::from("DUX4")).unwrap();
    let mut flanking_phasing_info: BTreeMap<String, AlleleFlankingPhasing> = BTreeMap::new();
    let downstream_haps = &downstream_phasing_result.final_haplotypes;
    let downstream_reads = &downstream_phasing_result.unique_supporting_reads;
    let downstream_sites = &downstream_phasing_result.sites_for_phasing;
    let downstream_details = &downstream_phasing_result.haplotype_details;
    let upstream_haps = &upstream_phasing_result.final_haplotypes;
    let upstream_reads = &upstream_phasing_result.unique_supporting_reads;
    let upstream_sites = &upstream_phasing_result.sites_for_phasing;
    let upstream_details = &upstream_phasing_result.haplotype_details;
    flanking_phasing_info.insert(
        String::from("upstream"),
        AlleleFlankingPhasing {
            final_haplotypes: upstream_haps.clone(),
            unique_supporting_reads: upstream_reads.clone(),
            sites_for_phasing: upstream_sites.clone(),
            haplotype_details: upstream_details.clone(),
        },
    );
    flanking_phasing_info.insert(
        String::from("downstream"),
        AlleleFlankingPhasing {
            final_haplotypes: downstream_haps.clone(),
            unique_supporting_reads: downstream_reads.clone(),
            sites_for_phasing: downstream_sites.clone(),
            haplotype_details: downstream_details.clone(),
        },
    );

    let mut allele_methyl_info = BTreeMap::new();
    allele_methyl_info.insert(String::from("complete_alleles"), allele_methyl);
    allele_methyl_info.insert(String::from("all_allele_ends"), all_ends_allele_methyl);
    // pack variants for reporting
    let mut complete_allele_variants = BTreeMap::new();
    for (allele, allele_variants) in variant_report.complete_allele_variants {
        let mut allele_variants_reformat = allele_variants.into_iter().collect::<Vec<_>>();
        allele_variants_reformat.sort_by(|a, b| a.0 .0.cmp(&b.0 .0));
        let allele_variants_reformat = allele_variants_reformat
            .iter()
            .map(|x| format!("{};{};{}", x.0 .0, x.0 .1, x.1.join(",")))
            .collect::<Vec<_>>();
        complete_allele_variants.insert(
            vec_to_string(&vec![allele], "-")
                .first()
                .unwrap()
                .to_string(),
            allele_variants_reformat,
        );
    }

    let mut allele_background_assignment: BTreeMap<String, BTreeMap<String, String>> =
        BTreeMap::new();
    allele_background_assignment.insert(String::from("complete_alleles"), hap_backgrounds);
    allele_background_assignment.insert(
        String::from("all_allele_starts"),
        all_starts_hap_backgrounds,
    );
    allele_background_assignment.insert(String::from("all_allele_ends"), all_ends_hap_backgrounds);

    let mut sample_call = SampleCall {
        allele_cn: allele_cn.to_string(),
        depth_summary,
        complete_alleles: final_complete,
        partial_alleles: vec_to_string(&assembly_result.incomplete, "-"),
        supporting_reads: support,
        complete_allele_variants,
        methylation: allele_methyl_info,
        other_unit_variants: variant_report
            .fp_variants_on_incomplete_alleles
            .into_iter()
            .map(|(a, b)| format!("{a};{}", b.join(",")))
            .collect::<Vec<String>>(),
        allele_info: allele_summaries,
        ..Default::default()
    };

    sample_call.additional.insert(
        String::from("allele_background"),
        serde_json::to_value(&allele_background_assignment)?.into(),
    );
    sample_call.additional.insert(
        String::from("allele_starts_mapped_to_upstream_haplotypes"),
        serde_json::to_value(&all_starts_upstream_haplotypes)?.into(),
    );

    let raw_complete = &assembly_result
        .complete
        .iter()
        .map(|x| {
            let y = x.iter().map(|a| a.to_string()).collect::<Vec<_>>();
            y.join("-")
        })
        .collect::<Vec<_>>();
    sample_call.additional.insert(
        String::from("raw_complete_alleles"),
        serde_json::to_value(&raw_complete)?.into(),
    );
    let raw_incomplete = &assembly_result
        .incomplete
        .iter()
        .map(|x| {
            let y = x.iter().map(|a| a.to_string()).collect::<Vec<_>>();
            y.join("-")
        })
        .collect::<Vec<_>>();
    sample_call.additional.insert(
        String::from("raw_partial_alleles"),
        serde_json::to_value(&raw_incomplete)?.into(),
    );
    let mut d4z4_arrays_in_cis = Vec::new();
    for phased_haps in &cis_dups {
        let phased_haps_string = vec_to_string(phased_haps, "-");
        d4z4_arrays_in_cis.push(phased_haps_string);
    }
    sample_call.additional.insert(
        String::from("d4z4_arrays_in_cis"),
        serde_json::to_value(&d4z4_arrays_in_cis)?.into(),
    );

    let mut meth_per_fp_median_reformat = meth_summary
        .meth_per_fp_median
        .into_iter()
        .collect::<Vec<_>>();
    meth_per_fp_median_reformat.sort_by(|a, b| a.0.cmp(&b.0));
    let meth_per_fp_median_reformat = meth_per_fp_median_reformat
        .iter()
        .map(|x| format!("{}:{}", x.0, x.1))
        .collect::<Vec<_>>();
    let mut meth_per_pos_median_reformat = meth_summary
        .meth_per_pos_median
        .into_iter()
        .collect::<Vec<_>>();
    meth_per_pos_median_reformat.sort_by(|a, b| a.0.cmp(&b.0));
    let meth_per_pos_median_reformat = meth_per_pos_median_reformat
        .iter()
        .map(|x| format!("{}:{}", x.0, x.1))
        .collect::<Vec<_>>();

    sample_call.additional.insert(
        String::from("median_methylation_all_sites"),
        meth_summary.all_sites_methyl_median.into(),
    );
    sample_call.additional.insert(
        String::from("median_methylation_per_unit"),
        serde_json::to_value(&meth_per_fp_median_reformat)?.into(),
    );
    sample_call.additional.insert(
        String::from("median_methylation_per_site"),
        serde_json::to_value(&meth_per_pos_median_reformat)?.into(),
    );
    sample_call.additional.insert(
        String::from("phasing_of_flanking_regions"),
        serde_json::to_value(&flanking_phasing_info)?.into(),
    );
    let d4z4_qc_metrics = D4Z4QCMetrics {
        median_read_length,
        per_allele_depth: haploid_depth,
    };
    sample_call.additional.insert(
        String::from("qc_metrics"),
        serde_json::to_value(&d4z4_qc_metrics)?.into(),
    );

    let mut writer = std::io::BufWriter::new(std::fs::File::create(output_json)?);
    writeln!(writer, "{}", serde_json::to_string_pretty(&sample_call)?)?;

    // write to vcf
    debug!("Write to VCF...");
    std::fs::File::create(output_vcf.clone())?;
    write_vcf(
        &output_vcf,
        &sample_id,
        variant_report.variant_summary,
        region_coordinates.clone(),
        assembly_result.complete.len(),
    )?;

    // plot
    debug!("Plot alleles...");
    let alleles_for_plot = variant_report.alleles_for_plot;
    if let Some(to_plot) = alleles_for_plot {
        plot_alleles_and_reads(&output_svg, to_plot)?;
    }

    // remove temporary reference file
    remove_if_exists(reference)?;
    let fai_file = output_path.join(format!("{sample_id}.d4z4.ref.fa.fai"));
    remove_if_exists(fai_file)?;
    remove_if_exists(genome_reference)?;
    let fai_file = output_path.join(format!("{sample_id}.d4z4.genome.fa.fai"));
    remove_if_exists(fai_file)?;

    info!("Completed kivvi analysis on D4Z4...");
    Ok(())
}
