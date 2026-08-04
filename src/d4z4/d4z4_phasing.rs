use crate::assembly::assembler::FpGraph;
use crate::caller::vec_to_string;
use crate::repeat_unit::fingerprint::FingerprintInfo;
use crate::util::{create_kivvi_temp_dir, invalid_data_error, missing_data_error, DError};
use crate::variant::get_read_position_in_allele;
use log::{debug, trace};
use paraphase::config::region::try_load;
use paraphase::io::bam::BamWriter;
use paraphase::io::json::GeneCall;
use paraphase::{config, phaser};
use rust_htslib::bam::{self, Read, Record};
use std::cmp;
use std::collections::{BTreeMap, HashSet};
use std::path::{Path, PathBuf};

const PARAPHASE_D4Z4_CONFIG: &[u8] = std::include_bytes!(concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/data/d4z4/paraphase_d4z4_config.yaml"
));

fn strip_chr_region_token(token: &str) -> String {
    if let Some((chr, rest)) = token.split_once(':') {
        format!("{}:{rest}", chr.strip_prefix("chr").unwrap_or(chr))
    } else {
        token.to_string()
    }
}

fn strip_chr_in_region_config_yaml(input: &[u8]) -> Result<Vec<u8>, DError> {
    let text = std::str::from_utf8(input)?;
    let mut out = Vec::<String>::new();
    for line in text.lines() {
        let trimmed = line.trim_start();
        if let Some(value) = trimmed.strip_prefix("realign_region: ") {
            let updated = strip_chr_region_token(value.trim());
            out.push(format!("  realign_region: {updated}"));
        } else if let Some(value) = trimmed.strip_prefix("extract_regions: ") {
            let updated = value
                .split_whitespace()
                .map(strip_chr_region_token)
                .collect::<Vec<_>>()
                .join(" ");
            out.push(format!("  extract_regions: {updated}"));
        } else {
            out.push(line.to_string());
        }
    }
    Ok(out.join("\n").into_bytes())
}

fn load_region_config_for_bam(wgs_bam: &PathBuf) -> Result<paraphase::config::Region, DError> {
    let reader = bam::Reader::from_path(wgs_bam)?;
    let bam_uses_chr = reader
        .header()
        .target_names()
        .into_iter()
        .any(|name| name.starts_with(b"chr"));
    let loaded_config_bytes = if bam_uses_chr {
        PARAPHASE_D4Z4_CONFIG.to_vec()
    } else {
        strip_chr_in_region_config_yaml(PARAPHASE_D4Z4_CONFIG)?
    };
    let region_config = try_load(Some(&loaded_config_bytes))?;
    debug!("paraphase region config {:?}", region_config);
    Ok(region_config)
}

fn paraphase_gene_bam_path(sample: &str, output_path: &Path, gene: &str) -> PathBuf {
    output_path.join(format!("{sample}.kivvi.paraphase.{gene}.bam"))
}

fn combined_paraphase_bam_path(sample: &str, output_path: &Path) -> PathBuf {
    output_path.join(format!("{sample}.kivvi.paraphase.bam"))
}

/// Fetch a named Paraphase gene call from the phasing result map.
/// # Arguments
/// * `phasing_result` - Paraphase gene calls keyed by gene name
/// * `gene` - gene name to look up
/// # Returns
/// * `&GeneCall` - gene call for the requested gene
fn get_paraphase_gene_call<'a>(
    phasing_result: &'a BTreeMap<String, GeneCall>,
    gene: &str,
) -> Result<&'a GeneCall, DError> {
    phasing_result
        .get(gene)
        .ok_or_else(|| missing_data_error("Paraphase result for gene", gene))
}

pub fn phase_flanking_gene(
    sample: &str,
    output_path: &Path,
    wgs_bam: &PathBuf,
    genome_reference: &PathBuf,
    gene: &str,
    write_bam: bool,
) -> Result<(GeneCall, Option<PathBuf>), DError> {
    debug!("Running Paraphase for flanking region {gene}");
    let region_config = load_region_config_for_bam(wgs_bam)?;
    let tmp_dir = create_kivvi_temp_dir(output_path)?;
    let reader = bam::Reader::from_path(wgs_bam)?;
    let settings = phaser::Settings::new(
        sample,
        (genome_reference, wgs_bam),
        tmp_dir.path(),
        gene.to_string(),
        &region_config,
        /* genome depth= */ None,
        /* sex = */ None,
        String::from("38"),
        None,
        0.03,
        false,
    );

    let mut phaser = phaser::Phaser::new(
        settings,
        Some(config::Gene::default()),
        None, // Option<SiteSelectionSettings>
        None, // Option<RealignSettings>
    );
    let res = phaser.run()?;

    let bam_path = if write_bam {
        let output_bam = paraphase_gene_bam_path(sample, output_path, gene);
        let mut writer = bam::Writer::from_path(
            &output_bam,
            &bam::Header::from_template(reader.header()),
            bam::Format::Bam,
        )?;
        let bam_writer = BamWriter::new(&phaser, &res);
        for item in bam_writer.write_bams()? {
            writer.write(&item)?;
        }
        Some(output_bam)
    } else {
        None
    };
    tmp_dir.close()?;
    Ok((res, bam_path))
}

pub fn merge_phasing_bams(
    sample: &str,
    output_path: &Path,
    wgs_bam: &PathBuf,
    bam_paths: &[PathBuf],
) -> Result<(), DError> {
    let output_bam = combined_paraphase_bam_path(sample, output_path);
    let header_reader = bam::Reader::from_path(wgs_bam)?;
    let mut writer = bam::Writer::from_path(
        &output_bam,
        &bam::Header::from_template(header_reader.header()),
        bam::Format::Bam,
    )?;
    let mut records = Vec::<Record>::new();
    for bam_path in bam_paths {
        let mut reader = bam::Reader::from_path(bam_path)?;
        for record in reader.records() {
            records.push(record?);
        }
    }
    records.sort_by(|a, b| a.tid().cmp(&b.tid()).then(a.pos().cmp(&b.pos())));
    for record in &records {
        writer.write(record)?;
    }
    for bam_path in bam_paths {
        if bam_path.exists() {
            std::fs::remove_file(bam_path)?;
        }
    }
    Ok(())
}

pub fn remove_phasing_bam(sample: &str, output_path: &Path) -> Result<(), DError> {
    let output_bam = combined_paraphase_bam_path(sample, output_path);
    if output_bam.exists() {
        std::fs::remove_file(output_bam)?;
    }
    Ok(())
}

/// Use paraphase to phase upstream regions
/// # Arguments
/// * `sample` - sample name
/// * `output_path` - output path path
/// * `wgs_bam` - wgs bam path
/// * `genome_reference` - genome reference path
/// * `write_bam` - whether to write bam
/// # Returns
/// * `BTreeMap<String, GeneCall>` - region name -> paraphase gene calls
pub fn phase_flanking(
    sample: &String,
    output_path: &Path,
    wgs_bam: &PathBuf,
    genome_reference: &PathBuf,
    write_bam: bool,
) -> Result<BTreeMap<String, GeneCall>, DError> {
    debug!("Running Paraphase for flanking region");
    let (dux4p5, dux4p5_bam) = phase_flanking_gene(
        sample,
        output_path,
        wgs_bam,
        genome_reference,
        "DUX4p5",
        write_bam,
    )?;
    let (dux4, dux4_bam) = phase_flanking_gene(
        sample,
        output_path,
        wgs_bam,
        genome_reference,
        "DUX4",
        write_bam,
    )?;
    let mut ret = BTreeMap::<String, _>::new();
    ret.insert(String::from("DUX4p5"), dux4p5);
    ret.insert(String::from("DUX4"), dux4);
    if write_bam {
        let dux4p5_bam = dux4p5_bam
            .ok_or_else(|| missing_data_error("DUX4p5 Paraphase BAM path", "BAM output request"))?;
        let dux4_bam = dux4_bam
            .ok_or_else(|| missing_data_error("DUX4 Paraphase BAM path", "BAM output request"))?;
        merge_phasing_bams(sample, output_path, wgs_bam, &[dux4p5_bam, dux4_bam])?;
    } else {
        remove_phasing_bam(sample, output_path)?;
    }
    Ok(ret)
}

/// Determine the chromosome backgounds upstream/downstream of d4z4
/// # Arguments
/// * `complete_alleles` - complete alleles
/// * `phasing_result` - paraphase gene calls
/// * `supporting_reads` - supporting reads
/// * `fp_info` - fingerprint information
/// * `all_ends_reads_match_allele_index` - all ends reads match allele index
/// * `bases_at_pivot_site` - bases at pivot site
/// * `check_chromosome` - whether to check chromosome
/// * `check_polya` - whether to check polya
/// # Returns
/// * `BTreeMap<String, String>` - allele -> background
/// * `BTreeMap<String, Vec<String>>` - map an allele to its upstream paraphase haplotypes
pub fn haplotype_background(
    complete_alleles: &Vec<Vec<i32>>,
    phasing_result: &BTreeMap<String, GeneCall>,
    supporting_reads: &BTreeMap<Vec<i32>, HashSet<String>>,
    fp_info: Option<&FingerprintInfo>,
    all_ends_reads_match_allele_index: &BTreeMap<Vec<i32>, Vec<(String, i32)>>,
    bases_at_pivot_site: &BTreeMap<String, String>,
    check_chromosome: bool,
    check_polya: bool,
) -> Result<(BTreeMap<String, String>, BTreeMap<String, Vec<String>>), DError> {
    let mut hap_backgrounds = BTreeMap::new();
    let mut upstream_haplotypes = BTreeMap::new();
    // phasing upstream
    let upstream_phasing = get_paraphase_gene_call(phasing_result, "DUX4p5")?;
    let paraphase_reads = &upstream_phasing.unique_supporting_reads;
    let paraphase_reads_nonunique = &upstream_phasing.nonunique_supporting_reads;
    let paraphase_phasing_sites = &upstream_phasing.sites_for_phasing;
    let important_sites = vec![
        "54271_T_A",
        "54812_C_G",
        "54902_T_C",
        "54922_T_C",
        "54933_A_G",
        "54937_C_T",
        "54995_A_G",
        "55046_T_G",
        "55069_T_A",
        "55077_G_A",
    ]
    .iter()
    .map(|x| x.to_string())
    .collect::<Vec<String>>();
    let important_sites_index = important_sites
        .iter()
        .filter(|x| paraphase_phasing_sites.contains(x))
        .filter_map(|x| paraphase_phasing_sites.iter().position(|y| y == x))
        .collect::<Vec<usize>>();
    let mut var_to_group = BTreeMap::new();
    var_to_group.insert("1112111111".to_string(), String::from("Group1.1"));
    var_to_group.insert("1111111111".to_string(), String::from("Group1.2"));
    var_to_group.insert("1221222222".to_string(), String::from("Group2.1"));
    var_to_group.insert("2221222222".to_string(), String::from("Group2.2"));
    let paraphase_haplotypes_assignment = assign_paraphase_haplotypes_to_chromsome(phasing_result)?;
    debug!(
        "paraphase_haplotypes_assignment {:?}",
        paraphase_haplotypes_assignment
    );
    let mut paraphase_read_to_hap = BTreeMap::new();
    for (hap, reads) in paraphase_reads {
        for read in reads {
            paraphase_read_to_hap.insert(read.to_string(), hap.to_string());
        }
    }
    for allele in complete_alleles {
        debug!("checking background for allele {allele:?}");
        let mut flanking = Vec::new();
        let mut chromosome = String::from("chromosome_unknown");
        let mut upstream_group = String::from("upstream_group_unknown");
        let mut polya = String::from("unknown");
        if check_chromosome {
            let mut reads = Vec::new();
            let mut this_allele_paraphase_hap = Vec::new();
            if supporting_reads.contains_key(allele) {
                reads = supporting_reads
                    .get(allele)
                    .ok_or_else(|| {
                        missing_data_error("supporting reads for allele", format!("{allele:?}"))
                    })?
                    .iter()
                    .map(|x| x.to_string())
                    .collect::<Vec<_>>();
                let reads_found_by_paraphase = reads
                    .iter()
                    .filter(|x| paraphase_read_to_hap.contains_key(*x))
                    .collect::<Vec<_>>();
                debug!("reads_found_by_paraphase {reads_found_by_paraphase:?}");
                this_allele_paraphase_hap = reads_found_by_paraphase
                    .iter()
                    .map(|x| {
                        paraphase_read_to_hap.get(*x).cloned().ok_or_else(|| {
                            missing_data_error("Paraphase read-to-haplotype entry", x.to_string())
                        })
                    })
                    .collect::<Result<Vec<String>, _>>()?;
            }
            debug!("reads first attempt {reads:?}");
            debug!("this_allele_paraphase_hap first attempt {this_allele_paraphase_hap:?}");
            if this_allele_paraphase_hap.is_empty() {
                if let Some(fp_info) = fp_info {
                    let read_nodes = &fp_info.read_edges;
                    for (read, this_read_nodes) in read_nodes.iter() {
                        let Some(this_read_nodes_first) = this_read_nodes.first() else {
                            continue;
                        };
                        let this_read_len = this_read_nodes.len();
                        if *this_read_nodes_first < 0 && *this_read_nodes_first > -10 {
                            let check_size = cmp::min(3, allele.len());
                            if this_read_len >= check_size
                                && this_read_nodes[0..check_size] == allele[0..check_size]
                            {
                                let read_name = read
                                    .split_terminator(':')
                                    .next()
                                    .ok_or_else(|| {
                                        invalid_data_error(format!(
                                            "Read segment key is missing a read name prefix: '{read}'"
                                        ))
                                    })?
                                    .to_string();
                                reads.push(read_name);
                            }
                        }
                    }
                    this_allele_paraphase_hap = reads
                        .iter()
                        .filter(|x| paraphase_read_to_hap.contains_key(*x))
                        .map(|x| {
                            paraphase_read_to_hap.get(x).cloned().ok_or_else(|| {
                                missing_data_error(
                                    "Paraphase read-to-haplotype entry",
                                    x.to_string(),
                                )
                            })
                        })
                        .collect::<Result<Vec<String>, _>>()?;
                    debug!("reads looser check{reads:?}");
                    debug!("this_allele_paraphase_hap looser check {this_allele_paraphase_hap:?}");
                }
            }

            let mut matching_paraphase_haplotype_segments = Vec::new();
            for hap in &this_allele_paraphase_hap {
                let this_hap_chars = hap.clone().chars().collect::<Vec<char>>();
                let this_hap_segment = important_sites_index
                    .iter()
                    .map(|x| this_hap_chars[*x].to_string())
                    .collect::<Vec<String>>()
                    .join("");
                matching_paraphase_haplotype_segments.push(this_hap_segment);
            }
            debug!(
                "matching_paraphase_haplotype_segments {matching_paraphase_haplotype_segments:?}"
            );
            if !matching_paraphase_haplotype_segments.is_empty() {
                let all_same = matching_paraphase_haplotype_segments
                    .iter()
                    .all(|x| *x == matching_paraphase_haplotype_segments[0]);
                if all_same {
                    let variant_segment = &matching_paraphase_haplotype_segments[0];
                    if let Some(group) = var_to_group.get(variant_segment) {
                        upstream_group = group.to_string();
                    }
                }
            }

            let this_allele_paraphase_hap_assignment = this_allele_paraphase_hap
                .iter()
                .map(|x| {
                    paraphase_haplotypes_assignment
                        .get(x)
                        .cloned()
                        .ok_or_else(|| {
                            missing_data_error(
                                "chromosome assignment for Paraphase haplotype",
                                x.to_string(),
                            )
                        })
                })
                .collect::<Result<Vec<String>, _>>()?
                .into_iter()
                .collect::<counter::Counter<String, i64>>()
                .most_common_ordered();

            for hap in &this_allele_paraphase_hap {
                if !flanking.contains(hap) {
                    flanking.push(hap.to_string());
                }
            }
            if !this_allele_paraphase_hap_assignment.is_empty() {
                if this_allele_paraphase_hap_assignment.len() == 1 {
                    chromosome = this_allele_paraphase_hap_assignment[0].0.clone();
                } else if this_allele_paraphase_hap_assignment[0].1 > 1
                    && this_allele_paraphase_hap_assignment[1].1 <= 1
                {
                    chromosome = this_allele_paraphase_hap_assignment[0].0.clone();
                }
            } else {
                // no read is uniquely assigned haplotypes
                // use nonunique read assignments
                debug!("Using nonunique reads to determine chromosome and upstream group...");
                let mut this_allele_paraphase_hap_chrs = Vec::new();
                let mut this_allele_paraphase_hap_groups = Vec::new();
                let mut this_allele_paraphase_haps = HashSet::new();
                for read in &reads {
                    if paraphase_reads_nonunique.contains_key(read) {
                        let mut this_read_assignment_nonunique = Vec::new();
                        let nonunique = paraphase_reads_nonunique.get(read).ok_or_else(|| {
                            missing_data_error(
                                "Paraphase nonunique supporting reads for read",
                                read.to_string(),
                            )
                        })?;
                        for hap in nonunique {
                            let hap_assignment = paraphase_haplotypes_assignment
                                .get(hap)
                                .cloned()
                                .ok_or_else(|| {
                                    missing_data_error(
                                        "chromosome assignment for Paraphase haplotype",
                                        hap.to_string(),
                                    )
                                })?;
                            this_allele_paraphase_haps.insert(hap.clone());
                            this_read_assignment_nonunique.push(hap_assignment);

                            let this_hap_chars = hap.clone().chars().collect::<Vec<char>>();
                            let this_hap_found_segment = important_sites_index
                                .iter()
                                .map(|x| this_hap_chars[*x].to_string())
                                .collect::<Vec<String>>()
                                .join("");
                            this_allele_paraphase_hap_groups.push(this_hap_found_segment);
                        }
                        let this_read_assignment_nonunique_counter = this_read_assignment_nonunique
                            .into_iter()
                            .collect::<counter::Counter<String, i64>>()
                            .most_common_ordered();
                        if !this_read_assignment_nonunique_counter.is_empty()
                            && this_read_assignment_nonunique_counter.len() == 1
                        {
                            let this_read_chr = this_read_assignment_nonunique_counter[0].0.clone();
                            if this_read_chr != String::from("chromosome_unknown") {
                                this_allele_paraphase_hap_chrs.push(this_read_chr);
                            }
                        }
                    }
                }
                flanking = this_allele_paraphase_haps
                    .into_iter()
                    .collect::<Vec<String>>();
                let this_allele_paraphase_hap_chrs_counter = this_allele_paraphase_hap_chrs
                    .into_iter()
                    .collect::<counter::Counter<String, i64>>()
                    .most_common_ordered();
                if !this_allele_paraphase_hap_chrs_counter.is_empty()
                    && this_allele_paraphase_hap_chrs_counter.len() == 1
                {
                    chromosome = this_allele_paraphase_hap_chrs_counter[0].0.clone();
                }
                // check upstream groups
                debug!("Using nonunique reads: this_allele_paraphase_hap_groups {this_allele_paraphase_hap_groups:?}");
                if !this_allele_paraphase_hap_groups.is_empty() {
                    let all_same = this_allele_paraphase_hap_groups
                        .iter()
                        .all(|x| *x == this_allele_paraphase_hap_groups[0]);
                    if all_same {
                        let variant_segment = &this_allele_paraphase_hap_groups[0];
                        if let Some(group) = var_to_group.get(variant_segment) {
                            upstream_group = group.to_string();
                        }
                    }
                }
            }
        }

        if all_ends_reads_match_allele_index.contains_key(allele) {
            let reads = all_ends_reads_match_allele_index
                .get(allele)
                .ok_or_else(|| {
                    missing_data_error(
                        "reads matching allele index for allele",
                        format!("{allele:?}"),
                    )
                })?;
            if let Some(fp_info) = fp_info {
                polya = get_polya(allele, reads, bases_at_pivot_site, fp_info)?;
            }
        }
        let allele_name = vec_to_string(&vec![allele.clone()], "-");
        let hap_string = &allele_name[0];
        upstream_haplotypes.insert(hap_string.clone(), flanking.clone());
        if check_chromosome && check_polya {
            hap_backgrounds.insert(
                hap_string.to_string(),
                format!("{polya}-{chromosome}:{upstream_group}"),
            );
        } else if check_chromosome {
            hap_backgrounds.insert(
                hap_string.to_string(),
                format!("{chromosome}:{upstream_group}"),
            );
        } else if check_polya {
            hap_backgrounds.insert(hap_string.to_string(), polya);
        }
    }
    Ok((hap_backgrounds, upstream_haplotypes))
}

/// polyA site
/// # Arguments
/// * `allele_old` - allele
/// * `reads` - reads
/// * `bases_at_pivot_site` - bases at pivot site
/// * `fp_info` - fingerprint information
/// # Returns
/// * `String` - polya site
fn get_polya(
    allele_old: &Vec<i32>,
    reads: &Vec<(String, i32)>,
    bases_at_pivot_site: &BTreeMap<String, String>,
    fp_info: &FingerprintInfo,
) -> Result<String, DError> {
    let hap_last = allele_old
        .last()
        .ok_or_else(|| missing_data_error("last node of allele", format!("{allele_old:?}")))?;
    if *hap_last < -10 {
        return Ok(String::from("qB"));
    }

    let mut polya = String::from("unknown");
    let mut polya_site_this_allele = Vec::new();
    let read_positions = &fp_info.read_positions;
    let read_nodes = &fp_info.read_edges;
    for (read, index) in reads {
        let this_read_nodes = read_nodes.get(read).ok_or_else(|| {
            missing_data_error(
                "read-node path while inferring polyA status",
                read.to_string(),
            )
        })?;
        let this_read_positions = read_positions.get(read).ok_or_else(|| {
            missing_data_error(
                "read-position path while inferring polyA status",
                read.to_string(),
            )
        })?;
        let end_index_on_read = allele_old.len() as i32 - 2 - *index;
        trace!(
            "{allele_old:?} {read} {this_read_nodes:?} {this_read_positions:?} {} {end_index_on_read}",
            *index
        );
        if end_index_on_read < this_read_nodes.len() as i32 {
            let end_position_on_read = this_read_positions[end_index_on_read as usize];
            let segment_name = format!("{read}:{end_position_on_read}");
            trace!("{allele_old:?} {read} {this_read_nodes:?} {this_read_positions:?} {} {end_index_on_read} {end_position_on_read}", *index);
            if bases_at_pivot_site.contains_key(&segment_name) {
                let read_base = bases_at_pivot_site.get(&segment_name).ok_or_else(|| {
                    missing_data_error("pivot-site base for read segment", segment_name.to_string())
                })?;
                polya_site_this_allele.push(read_base);
            }
        }
    }
    // if no unique reads supporting the allele end, just match the second last unit
    if polya_site_this_allele.is_empty() && *hap_last == -10 {
        let allele_size = allele_old.len();
        let second_last_node = &allele_old[allele_size - 2];
        debug!(
            "checking ends loosely for allele {allele_old:?} second_last_node {second_last_node}"
        );
        for (read, this_read_nodes) in read_nodes.iter() {
            let read_size = this_read_nodes.len();
            let mut found_ending_node = false;
            let mut ending_node_index = 0;
            for (i, node) in this_read_nodes.iter().enumerate() {
                if i < read_size - 1 && node == second_last_node && this_read_nodes[i + 1] == -10 {
                    ending_node_index = i;
                    found_ending_node = true;
                    break;
                }
            }
            if found_ending_node {
                let this_read_positions = read_positions.get(read).ok_or_else(|| {
                    missing_data_error(
                        "read-position path while inferring polyA status",
                        read.to_string(),
                    )
                })?;
                let end_position_on_read = this_read_positions[ending_node_index as usize];
                let segment_name = format!("{read}:{end_position_on_read}");
                trace!("{allele_old:?} {read} {this_read_nodes:?} {this_read_positions:?} {ending_node_index} {end_position_on_read}");
                if bases_at_pivot_site.contains_key(&segment_name) {
                    let read_base = bases_at_pivot_site.get(&segment_name).ok_or_else(|| {
                        missing_data_error(
                            "pivot-site base for read segment",
                            segment_name.to_string(),
                        )
                    })?;
                    polya_site_this_allele.push(read_base);
                }
            }
        }
    }
    if !polya_site_this_allele.is_empty() {
        let all_count = polya_site_this_allele.len();
        let count_t = polya_site_this_allele
            .iter()
            .filter(|x| **x == "ATTAAA")
            .count();
        let count_c = polya_site_this_allele
            .iter()
            .filter(|x| **x == "ATCAAA")
            .count();
        let count_alt = polya_site_this_allele
            .iter()
            .filter(|x| **x == "ATTTAA")
            .count();
        let base_count = polya_site_this_allele
            .into_iter()
            .collect::<counter::Counter<&String, i64>>()
            .most_common_ordered();
        if !base_count.is_empty() {
            let common_base = base_count[0].0.clone();
            if common_base == "ATTAAA" && count_c <= 1 && count_alt <= 1 {
                polya = String::from("qAIntactPolyA");
            } else if count_t <= 1 {
                polya = String::from("qADisruptedPolyA");
            }
            debug!(
                "allele {allele_old:?} polyA site {common_base:?} all_count {all_count} count_t {count_t}"
            );
        }
    }
    Ok(polya)
}

/// Given a haplotype, determine whether it's chr4 or chr10
/// # Arguments
/// * `phasing_result` - paraphase gene calls
/// # Returns
/// * `BTreeMap<String, String>` - haplotype -> chromosome
fn assign_paraphase_haplotypes_to_chromsome(
    phasing_result: &BTreeMap<String, GeneCall>,
) -> Result<BTreeMap<String, String>, DError> {
    let data = std::include_bytes!(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/data/d4z4/chr4_chr10_diff_sites.txt"
    ));
    let diff_sites = std::str::from_utf8(data)?
        .split_terminator('\n')
        .map(std::borrow::ToOwned::to_owned)
        .collect::<Vec<_>>();
    let phasing_result = get_paraphase_gene_call(phasing_result, "DUX4p5")?;
    let mut chromosome_assignment = BTreeMap::new();
    for (hap_seq, hap_name) in &phasing_result.final_haplotypes {
        let mut assignment = String::from("chromosome_unknown");
        let first_site = hap_seq.as_bytes().first().ok_or_else(|| {
            invalid_data_error(format!(
                "Paraphase haplotype sequence is empty for hap '{hap_name}'"
            ))
        })?;
        if *first_site != b'x' {
            if *first_site == b'0' {
                assignment = String::from("chr10");
            } else {
                assignment = String::from("chr4");
            }
        } else {
            let hap_detail = &phasing_result
                .haplotype_details
                .get(hap_name)
                .ok_or_else(|| {
                    missing_data_error("Paraphase haplotype details", hap_name.to_string())
                })?;
            let hap_boundary = &hap_detail.boundary;
            let bounds = hap_boundary
                .split_terminator('-')
                .map(|x| {
                    x.parse::<i64>().map_err(|e| {
                        format!(
                            "Failed to parse Paraphase boundary '{hap_boundary}' for hap '{hap_name}': {e}"
                        )
                    })
                })
                .collect::<Result<Vec<_>, _>>()?;
            let hap_variants = &hap_detail.variants;
            let nsites_covered = diff_sites
                .iter()
                .filter(|x| {
                    let pos = x
                        .split_terminator('_')
                        .next()
                        .and_then(|value| value.parse::<i64>().ok());
                    pos.is_some_and(|pos| pos > bounds[0] && pos < bounds[1])
                })
                .count();
            let nvariants_overlap = diff_sites
                .iter()
                .filter(|x| hap_variants.contains(*x))
                .count();
            debug!("paraphase hap {hap_seq} {hap_name} nsites_covered {nsites_covered} nvariants_overlap {nvariants_overlap}");
            if nsites_covered < 10 && nsites_covered >= 3 {
                if nvariants_overlap == nsites_covered {
                    assignment = String::from("chr10");
                } else if nvariants_overlap == 0 {
                    assignment = String::from("chr4");
                }
            } else if nsites_covered >= 10 {
                if nvariants_overlap >= nsites_covered - 1 {
                    assignment = String::from("chr10");
                } else if nvariants_overlap <= 1 {
                    assignment = String::from("chr4");
                }
            }
        }
        chromosome_assignment.insert(hap_seq.clone(), assignment);
    }
    // if only one missing, infer it
    let count_unknown = chromosome_assignment
        .values()
        .filter(|x| *x == "chromosome_unknown")
        .count();
    debug!("number of alleles missing chromosome assignment {count_unknown}");
    if count_unknown == 1 {
        let count_chr10 = chromosome_assignment
            .values()
            .filter(|x| *x == "chr10")
            .count();
        let count_chr4 = chromosome_assignment
            .values()
            .filter(|x| *x == "chr4")
            .count();
        debug!("count_chr10 {count_chr10} count_chr4 {count_chr4}");
        let unknown_hap = chromosome_assignment
            .iter()
            .filter(|(_x, y)| *y == "chromosome_unknown")
            .next()
            .ok_or_else(|| {
                missing_data_error(
                    "unknown chromosome assignment entry",
                    "exactly one unresolved Paraphase haplotype",
                )
            })?
            .0
            .clone();
        debug!("unknown_hap {unknown_hap}");

        if count_chr10 == 2 && count_chr4 == 1 {
            debug!("assigning chr4 to unknown_hap {unknown_hap}");
            chromosome_assignment.insert(unknown_hap, String::from("chr4"));
        } else if count_chr10 == 1 && count_chr4 == 2 {
            debug!("assigning chr10 to unknown_hap {unknown_hap}");
            chromosome_assignment.insert(unknown_hap, String::from("chr10"));
        }
    }
    Ok(chromosome_assignment)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_get_paraphase_gene_call_errors_on_missing_gene() {
        let error = get_paraphase_gene_call(&BTreeMap::new(), "DUX4p5")
            .expect_err("missing Paraphase genes should error");

        assert!(
            error
                .to_string()
                .contains("missing Paraphase result for gene: DUX4p5"),
            "unexpected error: {error}"
        );
    }

    #[test]
    fn test_get_polya_errors_on_missing_read_nodes() {
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

        let error = get_polya(
            &vec![1, 2, 3],
            &vec![(String::from("read1"), 0)],
            &BTreeMap::new(),
            &fp_info,
        )
        .expect_err("missing read paths should error while inferring polyA");

        assert!(
            error
                .to_string()
                .contains("missing read-node path while inferring polyA status: read1"),
            "unexpected error: {error}"
        );
    }
}

/// Remove redundant haps from a list of haps
/// # Arguments
/// * `haps_to_assess` - list of haps to assess
/// * `assembly_result` - assembly result
/// # Returns
/// * `Vec<Vec<i32>>` - list of haps after removing redundant ones
/// get the background of all starting haps
/// # Arguments
/// * `assembly_result` - assembly result
/// * `fp_graph` - fingerprint graph
/// * `phasing_result` - paraphase gene calls
/// * `bases_at_pivot_site` - bases at pivot site
/// * `fp_info` - fingerprint information
/// # Returns
/// * `BTreeMap<String, String>` - allele -> background
/// * `BTreeMap<String, Vec<String>>` - map an allele to its upstream paraphase haplotypes
pub fn get_background_for_allele_starts(
    kept_starting_haps: &[Vec<i32>],
    fp_graph: &FpGraph,
    phasing_result: &BTreeMap<String, GeneCall>,
    bases_at_pivot_site: &BTreeMap<String, String>,
    fp_info: &FingerprintInfo,
) -> Result<(BTreeMap<String, String>, BTreeMap<String, Vec<String>>), DError> {
    let all_starting_read_support = fp_graph
        .process_complete_haps(kept_starting_haps.to_vec(), Some(1), true, false)?
        .supporting_reads;
    debug!("all_starting_read_support {:?}", all_starting_read_support);

    // haplotype backgrounds
    let (all_starts_hap_backgrounds, upstream_haplotypes) = haplotype_background(
        &kept_starting_haps.to_vec(),
        phasing_result,
        &all_starting_read_support,
        Some(fp_info),
        &BTreeMap::new(),
        bases_at_pivot_site,
        true,
        false,
    )?;
    debug!(
        "all_starts_hap_backgrounds {:?}",
        all_starts_hap_backgrounds
    );
    Ok((all_starts_hap_backgrounds, upstream_haplotypes))
}

/// get the background of all ending haps
/// # Arguments
/// * `assembly_result` - assembly result
/// * `fp_graph` - fingerprint graph
/// * `phasing_result` - paraphase gene calls
/// * `bases_at_pivot_site` - bases at pivot site
/// * `fp_info` - fingerprint information
/// * `cis_dups_match_index` - cis duplicates match index
/// # Returns
/// * `BTreeMap<String, String>` - allele -> background
/// * `BTreeMap<Vec<i32>, Vec<(String, i32)>>` - all_ends_reads_match_allele_index
pub fn get_background_for_allele_ends(
    kept_ending_haps: &[Vec<i32>],
    fp_graph: &FpGraph,
    phasing_result: &BTreeMap<String, GeneCall>,
    bases_at_pivot_site: &BTreeMap<String, String>,
    fp_info: &FingerprintInfo,
    cis_dups_match_index: &BTreeMap<Vec<i32>, HashSet<(String, i32)>>,
) -> Result<
    (
        BTreeMap<String, String>,
        BTreeMap<Vec<i32>, Vec<(String, i32)>>,
    ),
    DError,
> {
    let all_ending_read_support = fp_graph
        .process_complete_haps(kept_ending_haps.to_vec(), Some(1), true, false)?
        .supporting_reads;
    debug!("all_ending_read_support {:?}", all_ending_read_support);
    // find index on reads
    let mut all_ends_reads_match_allele_index = get_read_position_in_allele(
        &fp_info.read_edges,
        kept_ending_haps,
        &all_ending_read_support,
        true,
    )?;
    for (allele, reads) in cis_dups_match_index {
        if all_ends_reads_match_allele_index.contains_key(allele) {
            if let Some(val) = all_ends_reads_match_allele_index.get_mut(allele) {
                for read in reads {
                    if !val.contains(read) {
                        val.push(read.clone());
                    }
                }
            }
        } else {
            for read in reads {
                all_ends_reads_match_allele_index
                    .entry(allele.clone())
                    .or_default()
                    .push(read.clone());
            }
        }
    }
    debug!(
        "all_ends_reads_match_allele_index {:?}",
        all_ends_reads_match_allele_index
    );

    // haplotype backgrounds
    let (all_ends_hap_backgrounds, _upstream_haplotypes) = haplotype_background(
        &kept_ending_haps.to_vec(),
        phasing_result,
        &all_ending_read_support,
        Some(fp_info),
        &all_ends_reads_match_allele_index,
        bases_at_pivot_site,
        false,
        true,
    )?;
    debug!("all_ends_hap_backgrounds {:?}", all_ends_hap_backgrounds);
    Ok((all_ends_hap_backgrounds, all_ends_reads_match_allele_index))
}
