use crate::assembly::assembler::{
    build_graph, match_reads_and_haplotypes, AssemblyResult, FpGraph,
};
use crate::assembly::assembler_utils::{
    compare_two_haps_same_length, find_overlapping_alleles, redundant_haplotype_allowed,
};
use crate::caller::vec_to_string;
use crate::d4z4::join_partial_alleles::is_cis_dup_by_read_start_offset;
use crate::repeat_unit::fingerprint::FingerprintInfo;
use crate::util::DError;
use crate::variant::get_read_position_in_allele;
use itertools::Itertools;
use log::{debug, trace};
use paraphase::config::region::try_load;
use paraphase::io::bam::BamWriter;
use paraphase::io::json::GeneCall;
use paraphase::{config, phaser};
use rust_htslib::bam::{self, Read, Record};
use std::cmp;
use std::collections::{BTreeMap, HashSet};
use std::path::{Path, PathBuf};

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
    const DATA: &[u8] = std::include_bytes!(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/data/d4z4/paraphase_d4z4_config.yaml"
    ));
    let mut ret = BTreeMap::<String, _>::new();
    let mut bam_ret = std::collections::BTreeMap::<u64, Vec<bam::Record>>::new();

    // temp dir
    let tmp_dir = tempfile::TempDir::new()?;
    let reader = bam::Reader::from_path(wgs_bam)?;
    let bam_uses_chr = reader
        .header()
        .target_names()
        .into_iter()
        .any(|name| name.starts_with(b"chr"));
    let loaded_config_bytes = if bam_uses_chr {
        DATA.to_vec()
    } else {
        strip_chr_in_region_config_yaml(DATA)?
    };
    let region_config = try_load(Some(&loaded_config_bytes))?;
    let genes = region_config.keys().cloned().rev().collect::<Vec<_>>();
    debug!("paraphase region config {:?}", region_config);
    let output_bam = output_path.join(format!("{sample}.kivvi.paraphase.bam"));
    let mut writer = bam::Writer::from_path(
        &output_bam,
        &bam::Header::from_template(reader.header()),
        bam::Format::Bam,
    )?;

    // Compute for each gene
    for gene in genes {
        let settings = phaser::Settings::new(
            sample,
            (genome_reference, wgs_bam),
            //&args.outdir,
            tmp_dir.path(),
            gene.clone(),
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

        // write to bam
        if write_bam {
            let bam_writer = BamWriter::new(&phaser, &res);
            let bam_records_out: Vec<Record> = bam_writer.write_bams()?;
            for item in bam_records_out {
                let tid_pos = ((item.tid() as u64) << 32) | item.pos() as u64;
                bam_ret.entry(tid_pos).or_default().push(item);
            }
        }
        ret.insert(gene, res);
    }
    let alignments = bam_ret
        .into_values()
        .flat_map(std::iter::IntoIterator::into_iter)
        .collect::<Vec<_>>();
    for align in &alignments {
        writer.write(align)?;
    }
    //bam::index::build(&output_bam, None, bam::index::Type::Bai, 1)?;
    tmp_dir.close()?;
    if !write_bam && output_bam.exists() {
        std::fs::remove_file(output_bam)?;
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
    let paraphase_reads = &phasing_result
        .get(&String::from("DUX4p5"))
        .unwrap()
        .unique_supporting_reads
        .clone();
    let paraphase_reads_nonunique = &phasing_result
        .get(&String::from("DUX4p5"))
        .unwrap()
        .nonunique_supporting_reads
        .clone();
    let paraphase_phasing_sites = &phasing_result
        .get(&String::from("DUX4p5"))
        .unwrap()
        .sites_for_phasing
        .clone();
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
        .map(|x| paraphase_phasing_sites.iter().position(|y| y == x).unwrap())
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
                    .ok_or("hap not in assembly_result.supporting_reads.")?
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
                    .map(|x| paraphase_read_to_hap.get(*x).unwrap().to_string())
                    .collect::<Vec<String>>();
            }
            debug!("reads first attempt {reads:?}");
            debug!("this_allele_paraphase_hap first attempt {this_allele_paraphase_hap:?}");
            if this_allele_paraphase_hap.is_empty() {
                if fp_info.is_some() {
                    let fp_info = fp_info.unwrap();
                    let read_nodes = &fp_info.read_edges;
                    for (read, this_read_nodes) in read_nodes.iter() {
                        let this_read_nodes_first = this_read_nodes.first().unwrap();
                        let this_read_len = this_read_nodes.len();
                        if *this_read_nodes_first < 0 && *this_read_nodes_first > -10 {
                            let check_size = cmp::min(3, allele.len());
                            if this_read_len >= check_size
                                && this_read_nodes[0..check_size] == allele[0..check_size]
                            {
                                let read_name =
                                    read.split_terminator(':').collect::<Vec<_>>()[0].to_string();
                                reads.push(read_name);
                            }
                        }
                    }
                    this_allele_paraphase_hap = reads
                        .iter()
                        .filter(|x| paraphase_read_to_hap.contains_key(*x))
                        .map(|x| paraphase_read_to_hap.get(x).unwrap().to_string())
                        .collect::<Vec<String>>();
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
                    if var_to_group.contains_key(variant_segment) {
                        upstream_group = var_to_group.get(variant_segment).unwrap().to_string();
                    }
                }
            }

            let this_allele_paraphase_hap_assignment = this_allele_paraphase_hap
                .iter()
                .map(|x| paraphase_haplotypes_assignment.get(x).unwrap().to_string())
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
                        let nonunique = paraphase_reads_nonunique.get(read).unwrap();
                        for hap in nonunique {
                            let hap_assignment = paraphase_haplotypes_assignment
                                .get(hap)
                                .unwrap()
                                .to_string();
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
                        if var_to_group.contains_key(variant_segment) {
                            upstream_group = var_to_group.get(variant_segment).unwrap().to_string();
                        }
                    }
                }
            }
        }

        if all_ends_reads_match_allele_index.contains_key(allele) {
            let reads = all_ends_reads_match_allele_index
                .get(allele)
                .ok_or("hap not in assembly_result.supporting_reads.")?;
            if !fp_info.is_none() {
                polya = get_polya(allele, reads, bases_at_pivot_site, fp_info.unwrap())?;
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
    let hap_last = allele_old.last().ok_or("last not in allele_old")?;
    if *hap_last < -10 {
        return Ok(String::from("qB"));
    }

    let mut polya = String::from("unknown");
    let mut polya_site_this_allele = Vec::new();
    let read_positions = &fp_info.read_positions;
    let read_nodes = &fp_info.read_edges;
    for (read, index) in reads {
        let this_read_nodes = read_nodes.get(read).unwrap();
        let this_read_positions = read_positions.get(read).unwrap();
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
                let read_base = bases_at_pivot_site
                    .get(&segment_name)
                    .ok_or("segment_name not in bases_at_pivot_site")?;
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
                let this_read_positions = read_positions.get(read).unwrap();
                let end_position_on_read = this_read_positions[ending_node_index as usize];
                let segment_name = format!("{read}:{end_position_on_read}");
                trace!("{allele_old:?} {read} {this_read_nodes:?} {this_read_positions:?} {ending_node_index} {end_position_on_read}");
                if bases_at_pivot_site.contains_key(&segment_name) {
                    let read_base = bases_at_pivot_site
                        .get(&segment_name)
                        .ok_or("segment_name not in bases_at_pivot_site")?;
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
            log::debug!("allele {allele_old:?} polyA site {common_base:?} all_count {all_count} count_t {count_t}");
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
    let diff_sites = std::str::from_utf8(data)
        .unwrap()
        .split_terminator('\n')
        .map(std::borrow::ToOwned::to_owned)
        .collect::<Vec<_>>();
    let phasing_result = phasing_result.get(&String::from("DUX4p5")).unwrap();
    let mut chromosome_assignment = BTreeMap::new();
    for (hap_seq, hap_name) in &phasing_result.final_haplotypes {
        let mut assignment = String::from("chromosome_unknown");
        let first_site = hap_seq.as_bytes().first().unwrap();
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
                .ok_or("hap_name not in haplotype_details")?;
            let hap_boundary = &hap_detail.boundary;
            let bounds = hap_boundary
                .split_terminator('-')
                .map(std::borrow::ToOwned::to_owned)
                .map(|x| x.parse::<i64>().unwrap())
                .collect::<Vec<_>>();
            let hap_variants = &hap_detail.variants;
            let nsites_covered = diff_sites
                .iter()
                .filter(|x| {
                    let pos = x.split_terminator('_').collect::<Vec<_>>()[0]
                        .parse::<i64>()
                        .unwrap();
                    pos > bounds[0] && pos < bounds[1]
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
            .unwrap()
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

/// Remove redundant haps from a list of haps
/// # Arguments
/// * `haps_to_assess` - list of haps to assess
/// * `assembly_result` - assembly result
/// # Returns
/// * `Vec<Vec<i32>>` - list of haps after removing redundant ones
#[allow(dead_code)]
fn remove_redundant_haps(
    haps_to_assess: &Vec<Vec<i32>>,
    assembly_result: &AssemblyResult,
) -> Result<Vec<Vec<i32>>, DError> {
    let complete = &assembly_result.complete;
    // if there are less than 4 haps, return the original list
    if haps_to_assess.len() <= 4 {
        return Ok(haps_to_assess.clone());
    }
    let overlapping_haps = find_overlapping_alleles(haps_to_assess.clone(), Some(5))?.1;
    let mut redundant_haps = Vec::new();
    for (hap1, hap1_overlaps) in overlapping_haps.iter() {
        let hap1_len = hap1.len();
        for (hap2, overlap_len) in hap1_overlaps.iter() {
            let hap2_len = hap2.len();
            if hap1_len == hap2_len && redundant_haps.contains(hap1) {
                continue;
            }
            if complete.contains(hap1) {
                if hap2_len == *overlap_len {
                    redundant_haps.push(hap2.clone());
                } else if redundant_haplotype_allowed(hap1, hap2, overlap_len)? {
                    redundant_haps.push(hap2.clone());
                }
            } else if hap1_len >= hap2_len {
                if hap2_len == *overlap_len {
                    redundant_haps.push(hap2.clone());
                } else if redundant_haplotype_allowed(hap1, hap2, overlap_len)? {
                    redundant_haps.push(hap2.clone());
                }
            }
        }
    }
    for hap1 in haps_to_assess {
        for hap2 in haps_to_assess {
            if hap1 != hap2 {
                let hap1_len = hap1.len();
                let hap2_len = hap2.len();
                if hap1_len == hap2_len && redundant_haps.contains(hap1) {
                    continue;
                }
                if hap1_len >= hap2_len {
                    if hap2_len > 5 && hap1[1..hap2_len] == hap2[1..hap2_len] {
                        if !redundant_haps.contains(hap2) {
                            redundant_haps.push(hap2.clone());
                            debug!("{hap2:?} is redundant with {hap1:?}");
                        }
                    }
                }
            }
        }
    }
    let redundant_haps = redundant_haps
        .iter()
        .filter(|x| !complete.contains(x))
        .cloned()
        .collect::<Vec<_>>();
    debug!("redundant_haps {redundant_haps:?}");
    let haps_to_return: Vec<Vec<i32>> = haps_to_assess
        .iter()
        .filter(|x| !redundant_haps.contains(x))
        .cloned()
        .collect();
    if haps_to_return.len() < 4 {
        return Ok(haps_to_assess.clone());
    }
    Ok(haps_to_return)
}

fn is_cis_dup_hap(hap: &[i32], fp_info: &FingerprintInfo) -> Result<bool, DError> {
    let Some(&first_node) = hap.first() else {
        return Ok(false);
    };
    if first_node < 0 && first_node > -10 {
        return Ok(false);
    }
    if first_node <= -10 {
        return Ok(true);
    }
    if let Some(first_node_seq) = fp_info.good_name_to_seq.get(&first_node) {
        if first_node_seq[0] == b'S' {
            return Ok(true);
        }
    }

    let mut supporting_reads = 0;
    let mut delayed_start_reads = 0;
    for (read, read_nodes) in &fp_info.read_edges {
        let Some(read_positions) = fp_info.read_positions.get(read) else {
            continue;
        };
        if read_nodes.len() != read_positions.len() {
            continue;
        }

        let start_idx = 0;
        if read_nodes[start_idx] == first_node {
            let overlap_len = cmp::min(read_nodes.len() - start_idx, hap.len());
            if overlap_len < 2 {
                continue;
            }

            let nodes_in_read = &read_nodes[start_idx..(start_idx + overlap_len)];
            let nodes_in_hap = &hap[..overlap_len];
            let mut match_count = 0;
            let mut has_mismatch = false;
            for (read_node, hap_node) in nodes_in_read.iter().zip(nodes_in_hap.iter()) {
                if *read_node == 0 {
                    continue;
                }
                if read_node == hap_node {
                    match_count += 1;
                } else {
                    has_mismatch = true;
                    break;
                }
            }

            if !has_mismatch && match_count > 1 {
                supporting_reads += 1;
                if read_positions[start_idx] > 300 {
                    delayed_start_reads += 1;
                }
            }
        }
    }
    let delayed_start_threshold = (supporting_reads as f64 * 0.8).floor() as i32;
    Ok(supporting_reads >= 3
        && delayed_start_reads >= (supporting_reads - 1).min(delayed_start_threshold))
}

fn remove_redundant_haplotypes(
    haps_to_check: &[Vec<i32>],
    num_turns: usize,
) -> Result<BTreeMap<Vec<i32>, Vec<i32>>, DError> {
    let mut haps_to_remove = BTreeMap::<Vec<i32>, Vec<i32>>::new();
    for _turn_index in 0..num_turns {
        let haps_to_check = haps_to_check
            .iter()
            .filter(|hap| !haps_to_remove.contains_key(*hap))
            .cloned()
            .collect::<Vec<_>>();
        let (_overlapping_haps, overlapping_haps_match) =
            find_overlapping_alleles(haps_to_check.clone(), None)?;
        let mut removed_one_redundant = false;
        for (hap, hap_match_info) in &overlapping_haps_match {
            let hap_size = hap.len();
            for (matching_hap, overlap_len) in hap_match_info {
                let matching_hap_size = matching_hap.len();
                if *overlap_len >= 4 && *overlap_len >= hap_size / 2 {
                    if hap_size > matching_hap_size {
                        continue;
                    }
                    if hap_size == matching_hap_size {
                        let mut haps = vec![hap.clone(), matching_hap.clone()];
                        haps.sort();
                        let hap1 = haps[0].clone();
                        let hap2 = haps[1].clone();
                        if !haps_to_remove.contains_key(&hap1) {
                            haps_to_remove.insert(hap1, hap2);
                            removed_one_redundant = true;
                            break;
                        }
                    }
                    if hap_size < matching_hap_size {
                        if !haps_to_remove.contains_key(hap) {
                            haps_to_remove.insert(hap.clone(), matching_hap.clone());
                            removed_one_redundant = true;
                            break;
                        }
                    }
                }
            }
            if removed_one_redundant {
                break;
            }
        }
    }
    Ok(haps_to_remove)
}

pub fn process_alleles(
    assembly_result: &AssemblyResult,
    fp_info: &FingerprintInfo,
) -> Result<(Vec<Vec<i32>>, Vec<Vec<i32>>, Vec<Vec<i32>>), DError> {
    let mut all_starting_haps = HashSet::new();
    let mut all_ending_haps = HashSet::new();
    for hap in &assembly_result.complete {
        let hap_first = hap.first().unwrap();
        if *hap_first < 0 && *hap_first > -10 {
            all_starting_haps.insert(hap.clone());
        }
        let hap_end = hap.last().unwrap();
        if *hap_end <= -10 {
            all_ending_haps.insert(hap.clone());
        }
    }
    for hap in &assembly_result.incomplete {
        let hap_first = hap.first().unwrap();
        if *hap_first < 0 && *hap_first > -10 {
            all_starting_haps.insert(hap.clone());
        }
        let hap_end = hap.last().unwrap();
        if *hap_end <= -10 {
            all_ending_haps.insert(hap.clone());
        }
    }
    debug!("all_starting_haps before removing redundant {all_starting_haps:?}");
    debug!("all_ending_haps before removing redundant {all_ending_haps:?}");
    let mut kept_starting_haps: Vec<Vec<i32>> = all_starting_haps.into_iter().collect();
    let mut kept_ending_haps: Vec<Vec<i32>> = all_ending_haps.into_iter().collect();
    let mut kept_complete = assembly_result.complete.clone();
    // let kept_starting_haps = remove_redundant_haps(&kept_starting_haps, assembly_result)?;
    // let kept_ending_haps = remove_redundant_haps(&kept_ending_haps, assembly_result)?;
    let mut kept_complete_set = kept_complete.iter().cloned().collect::<HashSet<_>>();

    let distal_no_cis_dup = kept_ending_haps
        .iter()
        .filter(|hap| !is_cis_dup_hap(hap, fp_info).unwrap_or(false))
        .cloned()
        .collect::<Vec<Vec<i32>>>();
    let size1_allele = distal_no_cis_dup
        .iter()
        .filter(|hap| hap.len() == 2)
        .cloned()
        .collect::<Vec<Vec<i32>>>();
    if size1_allele.len() == 1 && distal_no_cis_dup.len() >= 5 {
        if let Some(size1_allele) = size1_allele.first() {
            kept_ending_haps.retain(|hap| hap != size1_allele);
            kept_complete_set.remove(size1_allele);
        }
    }
    if kept_starting_haps.len() >= 5 && kept_ending_haps.len() >= 4 {
        let num_turns = kept_starting_haps.len() - 4;
        let proximal_to_remove = remove_redundant_haplotypes(&kept_starting_haps, num_turns)?;
        if proximal_to_remove.len() <= num_turns {
            for (proximal_to_remove_allele, redundant_allele) in &proximal_to_remove {
                if !(kept_ending_haps.len() == 4
                    && kept_ending_haps.contains(proximal_to_remove_allele))
                {
                    kept_starting_haps.retain(|hap| hap != proximal_to_remove_allele);
                    kept_complete_set.remove(proximal_to_remove_allele);
                } else if !(kept_ending_haps.len() == 4
                    && kept_ending_haps.contains(redundant_allele))
                {
                    kept_starting_haps.retain(|hap| hap != redundant_allele);
                    kept_complete_set.remove(redundant_allele);
                }
            }
        }
    }
    let distal_no_cis_dup = kept_ending_haps
        .iter()
        .filter(|hap| !is_cis_dup_hap(hap, fp_info).unwrap_or(false))
        .cloned()
        .collect::<Vec<Vec<i32>>>();
    if kept_starting_haps.len() == 4 && distal_no_cis_dup.len() >= 5 {
        let num_turns = distal_no_cis_dup.len() - 4;
        let distal_to_remove = remove_redundant_haplotypes(&distal_no_cis_dup, num_turns)?;
        if distal_to_remove.len() <= num_turns {
            for (distal_to_remove_allele, redundant_allele) in &distal_to_remove {
                if !(kept_starting_haps.len() == 4
                    && kept_starting_haps.contains(distal_to_remove_allele))
                {
                    kept_ending_haps.retain(|hap| hap != distal_to_remove_allele);
                    kept_complete_set.remove(distal_to_remove_allele);
                } else if !(kept_starting_haps.len() == 4
                    && kept_starting_haps.contains(redundant_allele))
                {
                    kept_ending_haps.retain(|hap| hap != redundant_allele);
                    kept_complete_set.remove(redundant_allele);
                }
            }
        }
    }

    kept_complete.retain(|hap| kept_complete_set.contains(hap));

    debug!("all_starting_haps after removing redundant {kept_starting_haps:?}");
    debug!("all_ending_haps after removing redundant {kept_ending_haps:?}");
    debug!("complete_haps after removing redundant {kept_complete:?}");

    Ok((kept_starting_haps, kept_ending_haps, kept_complete))
}

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
        fp_info.read_edges.clone(),
        kept_ending_haps.to_vec(),
        all_ending_read_support.clone(),
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

pub fn find_qal_alleles(all_haps: &[Vec<i32>], qal_units: Vec<i32>) -> HashSet<Vec<i32>> {
    all_haps
        .iter()
        .filter(|hap| hap.last() == Some(&-10) && hap.len() >= 2)
        .filter_map(|hap| {
            let second_to_last_unit = hap[hap.len() - 2];
            if qal_units.contains(&second_to_last_unit) {
                Some(hap.clone())
            } else {
                None
            }
        })
        .collect::<HashSet<_>>()
}

/// find in-cis duplications
/// # Arguments
/// * `assembly_result` - assembly result
/// * `fp_graph` - fingerprint graph
/// * `fp_info` - fingerprint information
/// * `phasing_result` - paraphase gene calls
/// # Returns
/// * `Vec<Vec<Vec<i32>>>` - cis duplications identified (allele of alleles)
/// * `BTreeMap<Vec<i32>, HashSet<(String, i32)>>` - allele -> (read, match_index_on_allele)
pub fn find_cis_dup(
    all_haps: &Vec<Vec<i32>>,
    fp_graph: &FpGraph,
    fp_info: &FingerprintInfo,
    phasing_result: &BTreeMap<String, GeneCall>,
    qal_units: Vec<i32>,
    special_incomplete_haps: &Vec<Vec<i32>>,
) -> Result<
    (
        Vec<Vec<Vec<i32>>>,
        BTreeMap<Vec<i32>, HashSet<(String, i32)>>,
    ),
    DError,
> {
    let mut haps_to_node_names: BTreeMap<Vec<i32>, i32> = BTreeMap::new();
    let mut read_edges_for_haps: BTreeMap<String, Vec<i32>> = BTreeMap::new();
    let mut node_name = 1;

    let qal_alleles = find_qal_alleles(&all_haps, qal_units);
    debug!("qal_alleles identified by long insertion before -10 {qal_alleles:?}");
    let cis_dup_alleles = all_haps
        .iter()
        .filter_map(|hap| {
            let allele_name = vec_to_string(&vec![hap.clone()], "-")
                .into_iter()
                .next()
                .unwrap_or_default();
            match is_cis_dup_by_read_start_offset(&allele_name, fp_info) {
                Ok(true) => Some(Ok(hap.clone())),
                Ok(false) => None,
                Err(err) => Some(Err(err)),
            }
        })
        .collect::<Result<HashSet<_>, _>>()?;
    debug!("cis_dup_alleles identified by read start offset {cis_dup_alleles:?}");
    let all_haps_support = fp_graph
        .process_complete_haps(all_haps.clone(), Some(1), true, false)?
        .support_by_read;
    let mut allele_links: BTreeMap<Vec<i32>, Vec<Vec<i32>>> = BTreeMap::new();
    let downstream_phasing_result = phasing_result.get(&String::from("DUX4")).unwrap();
    let downstream_reads = &downstream_phasing_result.unique_supporting_reads;
    // check links supported by paraphase haplotypes
    for (downstream_hap, downstream_hap_reads) in downstream_reads {
        trace!("checking reads for downstream_hap {downstream_hap}");
        let mut paraphase_hap_linking_repeat_haps_upstream = Vec::new();
        let mut paraphase_hap_linking_repeat_haps_downstream = Vec::new();
        for downstream_hap_read in downstream_hap_reads {
            let fields = downstream_hap_read.split("_sup_").collect::<Vec<_>>();
            let downstream_hap_read_name = fields[0];
            let aln_pos = fields[1]
                .split('_')
                .next()
                .ok_or("next not found")?
                .parse::<i32>()?;
            trace!("downstream_hap_read {downstream_hap_read} downstream_hap_read_name {downstream_hap_read_name} pos {aln_pos}");
            if all_haps_support.contains_key(downstream_hap_read_name) {
                let this_read_repeat_support =
                    all_haps_support.get(downstream_hap_read_name).unwrap();
                let this_read_repeat_edges =
                    fp_info.read_edges.get(downstream_hap_read_name).unwrap();
                let this_read_repeat_positions = fp_info
                    .read_positions
                    .get(downstream_hap_read_name)
                    .unwrap();
                trace!("this_read_repeat_edges {this_read_repeat_edges:?}");
                trace!("this_read_repeat_positions {this_read_repeat_positions:?}");
                for repeat_hap in this_read_repeat_support {
                    let node_index =
                        match_read_allele_first_node_index(this_read_repeat_edges, repeat_hap);
                    trace!("matching repeat_hap {repeat_hap:?} node_index {node_index:?}");
                    if let Some(node_index_value) = node_index {
                        let matching_position_on_read =
                            this_read_repeat_positions[node_index_value.0];
                        trace!("matching_position_on_read {matching_position_on_read}");
                        if matching_position_on_read < aln_pos {
                            paraphase_hap_linking_repeat_haps_upstream.push(repeat_hap.clone());
                        } else if matching_position_on_read > aln_pos {
                            paraphase_hap_linking_repeat_haps_downstream.push(repeat_hap.clone());
                        }
                    }
                }
            }
        }
        let paraphase_hap_linking_repeat_haps_upstream_set =
            paraphase_hap_linking_repeat_haps_upstream
                .iter()
                .map(|x| x.clone())
                .collect::<HashSet<Vec<i32>>>();
        let paraphase_hap_linking_repeat_haps_downstream_set =
            paraphase_hap_linking_repeat_haps_downstream
                .iter()
                .map(|x| x.clone())
                .collect::<HashSet<Vec<i32>>>();
        debug!("paraphase_hap {downstream_hap} linking repeat haps up {paraphase_hap_linking_repeat_haps_upstream_set:?} down {paraphase_hap_linking_repeat_haps_downstream_set:?}");
        if paraphase_hap_linking_repeat_haps_upstream_set.len() == 1
            && paraphase_hap_linking_repeat_haps_downstream_set.len() == 1
        {
            let a = paraphase_hap_linking_repeat_haps_upstream_set
                .iter()
                .next()
                .unwrap();
            let b = paraphase_hap_linking_repeat_haps_downstream_set
                .iter()
                .next()
                .unwrap();
            trace!("adding non-read links {a:?} to {b:?}");
            allele_links.entry(a.to_vec()).or_default().push(b.to_vec());
            if !haps_to_node_names.contains_key(a) {
                let a_name = node_name;
                haps_to_node_names.insert(a.to_vec(), a_name);
                node_name += 1;
            }
            if !haps_to_node_names.contains_key(b) {
                let b_name = node_name;
                haps_to_node_names.insert(b.to_vec(), b_name);
                node_name += 1;
            }
            let a_name = haps_to_node_names.get(a).unwrap();
            let b_name = haps_to_node_names.get(b).unwrap();
            read_edges_for_haps.insert(downstream_hap.to_string(), vec![*a_name, *b_name]);
            //allele_links.entry(b.to_vec()).or_default().push(a.to_vec());
        } else if paraphase_hap_linking_repeat_haps_upstream_set.len() == 2
            && paraphase_hap_linking_repeat_haps_downstream_set.len() == 2
        {
            let ovl = paraphase_hap_linking_repeat_haps_upstream_set
                .iter()
                .filter(|x| paraphase_hap_linking_repeat_haps_downstream_set.contains(*x))
                .map(|x| x.clone())
                .collect::<Vec<_>>();
            if ovl.len() == 1 {
                let middle = ovl.into_iter().next().unwrap();
                let first = paraphase_hap_linking_repeat_haps_upstream_set
                    .iter()
                    .filter(|x| **x != middle)
                    .map(|x| x.clone())
                    .collect::<Vec<_>>()
                    .first()
                    .unwrap()
                    .to_vec();
                let last = paraphase_hap_linking_repeat_haps_downstream_set
                    .iter()
                    .filter(|x| **x != middle)
                    .map(|x| x.clone())
                    .collect::<Vec<_>>()
                    .first()
                    .unwrap()
                    .to_vec();
                allele_links
                    .entry(first.to_vec())
                    .or_default()
                    .push(middle.to_vec());
                allele_links
                    .entry(middle.to_vec())
                    .or_default()
                    .push(last.to_vec());
                if !haps_to_node_names.contains_key(&first) {
                    let first_name = node_name;
                    haps_to_node_names.insert(first.to_vec(), first_name);
                    node_name += 1;
                }
                if !haps_to_node_names.contains_key(&middle) {
                    let middle_name = node_name;
                    haps_to_node_names.insert(middle.to_vec(), middle_name);
                    node_name += 1;
                }
                if !haps_to_node_names.contains_key(&last) {
                    let last_name = node_name;
                    haps_to_node_names.insert(last.to_vec(), last_name);
                    node_name += 1;
                }
                let first_name = haps_to_node_names.get(&first).unwrap();
                let middle_name = haps_to_node_names.get(&middle).unwrap();
                let last_name = haps_to_node_names.get(&last).unwrap();
                read_edges_for_haps.insert(
                    downstream_hap.to_string(),
                    vec![*first_name, *middle_name, *last_name],
                );
                trace!("adding non-read links {first:?} to {middle:?}");
                trace!("adding non-read links {middle:?} to {last:?}");
            }
        } else if paraphase_hap_linking_repeat_haps_upstream_set.len() == 2
            && paraphase_hap_linking_repeat_haps_downstream_set.len() == 1
        {
            let ovl = paraphase_hap_linking_repeat_haps_upstream_set
                .iter()
                .filter(|x| paraphase_hap_linking_repeat_haps_downstream_set.contains(*x))
                .map(|x| x.clone())
                .collect::<Vec<_>>();
            if ovl.len() == 1 {
                let b = ovl.into_iter().next().unwrap();
                let a = paraphase_hap_linking_repeat_haps_upstream_set
                    .iter()
                    .filter(|x| **x != b)
                    .map(|x| x.clone())
                    .collect::<Vec<_>>()
                    .first()
                    .unwrap()
                    .to_vec();
                allele_links.entry(a.to_vec()).or_default().push(b.to_vec());
                if !haps_to_node_names.contains_key(&a) {
                    let a_name = node_name;
                    haps_to_node_names.insert(a.to_vec(), a_name);
                    node_name += 1;
                }
                if !haps_to_node_names.contains_key(&b) {
                    let b_name = node_name;
                    haps_to_node_names.insert(b.to_vec(), b_name);
                    node_name += 1;
                }
                let a_name = haps_to_node_names.get(&a).unwrap();
                let b_name = haps_to_node_names.get(&b).unwrap();
                read_edges_for_haps.insert(downstream_hap.to_string(), vec![*a_name, *b_name]);
                trace!("adding non-read links {a:?} to {b:?}");
            }
        } else if paraphase_hap_linking_repeat_haps_downstream_set.len() <= 2 {
            let paraphase_hap_linking_repeat_haps_downstream_set_is_not_cis_dup =
                paraphase_hap_linking_repeat_haps_downstream_set
                    .iter()
                    .filter(|x| !cis_dup_alleles.contains(*x))
                    .map(|x| x.clone())
                    .collect::<HashSet<_>>();
            let paraphase_hap_linking_repeat_haps_upstream_set_is_not_cis_dup_and_is_qal =
                paraphase_hap_linking_repeat_haps_upstream_set
                    .iter()
                    .filter(|x| !cis_dup_alleles.contains(*x))
                    .filter(|x| qal_alleles.contains(*x))
                    .map(|x| x.clone())
                    .collect::<HashSet<_>>();
            let paraphase_hap_linking_repeat_haps_downstream_set_is_not_qal =
                paraphase_hap_linking_repeat_haps_downstream_set
                    .iter()
                    .filter(|x| !qal_alleles.contains(*x))
                    .map(|x| x.clone())
                    .collect::<HashSet<_>>();

            debug!(
                "paraphase_hap_linking_repeat_haps_downstream_set_is_not_cis_dup {:?}",
                paraphase_hap_linking_repeat_haps_downstream_set_is_not_cis_dup
            );
            debug!(
                "paraphase_hap_linking_repeat_haps_downstream_set_is_not_qal {:?}",
                paraphase_hap_linking_repeat_haps_downstream_set_is_not_qal
            );
            debug!(
                "paraphase_hap_linking_repeat_haps_upstream_set_is_not_cis_dup_and_is_qal {:?}",
                paraphase_hap_linking_repeat_haps_upstream_set_is_not_cis_dup_and_is_qal
            );

            if paraphase_hap_linking_repeat_haps_downstream_set_is_not_cis_dup.is_empty()
                && paraphase_hap_linking_repeat_haps_downstream_set_is_not_qal.is_empty()
                && paraphase_hap_linking_repeat_haps_upstream_set_is_not_cis_dup_and_is_qal.len()
                    == 1
            {
                if paraphase_hap_linking_repeat_haps_downstream_set.len() == 1 {
                    let b = paraphase_hap_linking_repeat_haps_downstream_set
                        .into_iter()
                        .next()
                        .unwrap();
                    let a =
                        paraphase_hap_linking_repeat_haps_upstream_set_is_not_cis_dup_and_is_qal
                            .into_iter()
                            .next()
                            .unwrap();
                    allele_links.entry(a.to_vec()).or_default().push(b.to_vec());
                    if !haps_to_node_names.contains_key(&a) {
                        let a_name = node_name;
                        haps_to_node_names.insert(a.to_vec(), a_name);
                        node_name += 1;
                    }
                    if !haps_to_node_names.contains_key(&b) {
                        let b_name = node_name;
                        haps_to_node_names.insert(b.to_vec(), b_name);
                        node_name += 1;
                    }
                    let a_name = haps_to_node_names.get(&a).unwrap();
                    let b_name = haps_to_node_names.get(&b).unwrap();
                    read_edges_for_haps.insert(downstream_hap.to_string(), vec![*a_name, *b_name]);
                    trace!("adding non-read links {a:?} to {b:?}");
                } else if paraphase_hap_linking_repeat_haps_downstream_set.len() == 2 {
                    let paraphase_hap_linking_repeat_haps_downstream_set =
                        paraphase_hap_linking_repeat_haps_downstream_set
                            .into_iter()
                            .sorted_by_key(|inner| std::cmp::Reverse(inner.len()))
                            .collect::<Vec<_>>();

                    let middle = paraphase_hap_linking_repeat_haps_downstream_set
                        .iter()
                        .next()
                        .unwrap()
                        .to_vec();
                    let first =
                        paraphase_hap_linking_repeat_haps_upstream_set_is_not_cis_dup_and_is_qal
                            .iter()
                            .next()
                            .unwrap()
                            .to_vec();
                    let last = paraphase_hap_linking_repeat_haps_downstream_set
                        .iter()
                        .last()
                        .unwrap()
                        .to_vec();
                    allele_links
                        .entry(first.to_vec())
                        .or_default()
                        .push(middle.to_vec());
                    allele_links
                        .entry(middle.to_vec())
                        .or_default()
                        .push(last.to_vec());
                    if !haps_to_node_names.contains_key(&first) {
                        let first_name = node_name;
                        haps_to_node_names.insert(first.to_vec(), first_name);
                        node_name += 1;
                    }
                    if !haps_to_node_names.contains_key(&middle) {
                        let middle_name = node_name;
                        haps_to_node_names.insert(middle.to_vec(), middle_name);
                        node_name += 1;
                    }
                    if !haps_to_node_names.contains_key(&last) {
                        let last_name = node_name;
                        haps_to_node_names.insert(last.to_vec(), last_name);
                        node_name += 1;
                    }
                    let first_name = haps_to_node_names.get(&first).unwrap();
                    let middle_name = haps_to_node_names.get(&middle).unwrap();
                    let last_name = haps_to_node_names.get(&last).unwrap();
                    read_edges_for_haps.insert(
                        downstream_hap.to_string(),
                        vec![*first_name, *middle_name, *last_name],
                    );
                    trace!("adding non-read links {first:?} to {middle:?}");
                    trace!("adding non-read links {middle:?} to {last:?}");
                }
            }
        }
    }
    // check reads now
    let mut match_index_on_read: BTreeMap<Vec<i32>, HashSet<(String, i32)>> = BTreeMap::new();
    for (read, read_nodes) in &fp_info.read_edges {
        if read_nodes.contains(&(-10)) {
            let end_index = read_nodes
                .iter()
                .position(|x| *x == -10)
                .ok_or("end (-10) not found")?;
            if end_index != read_nodes.len() - 1 {
                let this_read_repeat_positions = fp_info.read_positions.get(read).unwrap();
                let this_read_first_position = this_read_repeat_positions.first().unwrap();
                let mut segments: Vec<Vec<i32>> = Vec::new();
                let mut matching_alleles = Vec::new();
                //cis_read_edges.insert(read.clone(), read_nodes.to_vec());
                let mut starting_index = 0;
                let read_len = read_nodes.len();
                for i in 0..read_len {
                    if read_nodes[i] == -10 {
                        let this_seg = &read_nodes[starting_index..(i + 1)];
                        segments.push(this_seg.into());
                        starting_index = i + 1;
                    }
                }
                let this_seg = &read_nodes[starting_index..];
                if !this_seg.is_empty() {
                    segments.push(this_seg.into());
                }
                for k in 0..(segments.len()) {
                    let segment = &segments[k];
                    let mut dummy_read = BTreeMap::new();
                    dummy_read.insert(String::from("read"), segment.to_vec());
                    let segment_match = if k == 0 {
                        match_reads_and_haplotypes(dummy_read, all_haps.clone(), None, false)
                            .by_read
                    } else {
                        match_reads_and_haplotypes(
                            dummy_read,
                            special_incomplete_haps.clone(),
                            None,
                            false,
                        )
                        .by_read
                    };
                    if segment_match.contains_key("read") {
                        let dummy_read_matches = segment_match.get("read").unwrap();
                        if k == 0 {
                            let mut qualifying_matches = Vec::new();
                            for dummy_read_match in dummy_read_matches {
                                // only match against regular incomplete haps when the first position on read is small enough
                                if !special_incomplete_haps.contains(dummy_read_match) {
                                    qualifying_matches.push(dummy_read_match.to_vec());
                                } else {
                                    let node_index = match_read_allele_first_node_index(
                                        segment,
                                        dummy_read_match,
                                    );
                                    if !node_index.is_none()
                                        && node_index.unwrap().0 == 0
                                        && (*this_read_first_position < 2000
                                            || node_index.unwrap().1 == 0)
                                    {
                                        qualifying_matches.push(dummy_read_match.to_vec());
                                    }
                                }
                            }
                            if qualifying_matches.len() == 1 {
                                let qualifying_match = qualifying_matches.first().unwrap();
                                matching_alleles.push(qualifying_match.to_vec());
                            } else {
                                matching_alleles.push(vec![]);
                            }
                        } else {
                            let mut qualifying_matches = Vec::new();
                            for dummy_read_match in dummy_read_matches {
                                let node_index =
                                    match_read_allele_first_node_index(segment, dummy_read_match);
                                if !node_index.is_none()
                                    && node_index.unwrap().0 == 0
                                    && node_index.unwrap().1 == 0
                                {
                                    qualifying_matches.push(dummy_read_match.to_vec());
                                }
                            }
                            if qualifying_matches.len() == 1 {
                                let qualifying_match = qualifying_matches.first().unwrap();
                                matching_alleles.push(qualifying_match.to_vec());
                            } else {
                                matching_alleles.push(vec![]);
                            }
                        }
                    } else {
                        matching_alleles.push(vec![]);
                    }
                }
                trace!("read_nodes {read_nodes:?} segments {segments:?} matching_alleles {matching_alleles:?}");
                let n = segments.len();
                let mut prev_segments_len = 0;
                for j in 0..(n - 1) {
                    let match1 = &matching_alleles[j];
                    let match2 = &matching_alleles[j + 1];
                    if !match1.is_empty() && !match2.is_empty() {
                        trace!("adding read links {match1:?} to {match2:?}");
                        allele_links
                            .entry(match1.to_vec())
                            .or_default()
                            .push(match2.to_vec());

                        let node_index = match_read_allele_first_node_index(&segments[j], match1);
                        let mut match1_index_on_read: i32;
                        if node_index.unwrap().1 > 0 {
                            match1_index_on_read = node_index.unwrap().1 as i32;
                        } else {
                            match1_index_on_read =
                                0 - node_index.unwrap().0 as i32 - prev_segments_len as i32;
                            if match1.starts_with(&[-10]) {
                                match1_index_on_read += 1;
                            }
                        }
                        let mut match1_seg = match1.clone();
                        let match1_len = match1.len();
                        if match1_len > 6 {
                            match1_seg = vec![0, 0, 0, 0, 0, 0];
                            match1_seg.copy_from_slice(&match1[(match1_len - 6)..]);
                            match1_index_on_read -= match1_len as i32 - 6;
                        }
                        match_index_on_read
                            .entry(match1_seg.clone())
                            .or_default()
                            .insert((read.clone(), match1_index_on_read as i32));
                        let node_index =
                            match_read_allele_first_node_index(&segments[j + 1], match2);
                        let mut match2_index_on_read = 0
                            - node_index.unwrap().0 as i32
                            - prev_segments_len as i32
                            - segments[j].len() as i32;
                        if match2.starts_with(&[-10]) {
                            match2_index_on_read += 1;
                        }
                        let mut match2_seg = match2.clone();
                        let match2_len = match2.len();
                        if match2_len > 6 {
                            match2_seg = vec![0, 0, 0, 0, 0, 0];
                            match2_seg.copy_from_slice(&match2[(match2_len - 6)..]);
                            match2_index_on_read -= match2_len as i32 - 6;
                        }
                        match_index_on_read
                            .entry(match2_seg.clone())
                            .or_default()
                            .insert((read.clone(), match2_index_on_read));

                        if !haps_to_node_names.contains_key(match1) {
                            let match1_name = node_name;
                            haps_to_node_names.insert(match1.to_vec(), match1_name);
                            node_name += 1;
                        }
                        if !haps_to_node_names.contains_key(match2) {
                            let match2_name = node_name;
                            haps_to_node_names.insert(match2.to_vec(), match2_name);
                            node_name += 1;
                        }
                        let match1_name = haps_to_node_names.get(match1).unwrap();
                        let match2_name = haps_to_node_names.get(match2).unwrap();
                        if !read_edges_for_haps.contains_key(read) {
                            read_edges_for_haps
                                .insert(read.to_string(), vec![*match1_name, *match2_name]);
                        } else {
                            let new_read_name = format!("{read}_2");
                            if !read_edges_for_haps.contains_key(read) {
                                read_edges_for_haps
                                    .insert(new_read_name, vec![*match1_name, *match2_name]);
                            } else {
                                let new_read_name = format!("{read}_3");
                                if !read_edges_for_haps.contains_key(read) {
                                    read_edges_for_haps
                                        .insert(new_read_name, vec![*match1_name, *match2_name]);
                                }
                            }
                        }
                    }
                    prev_segments_len += segments[j].len();
                }
            }
        }
    }
    debug!("haps_to_node_names {haps_to_node_names:?}");
    debug!("read_edges_for_haps {read_edges_for_haps:?}");
    let mut node_names_to_haps = BTreeMap::new();
    for (a, b) in &haps_to_node_names {
        node_names_to_haps.insert(*b, a.to_vec());
    }
    debug!("special_incomplete {:?}", special_incomplete_haps);
    // simpler graph assembler
    debug!("For cis-dups, assemble haplotypes into alleles...");
    let mut fp_graph = build_graph(read_edges_for_haps, 2);
    let allele_phase_result = fp_graph.run_simple()?;
    debug!("allele_phase_result {:?}", allele_phase_result.incomplete);
    let mut cis_dups_assembled = Vec::new();
    for cis_dup in allele_phase_result.incomplete {
        let cis_dup_assembled = cis_dup
            .iter()
            .map(|x| node_names_to_haps.get(x).unwrap().to_vec())
            .collect::<Vec<_>>();
        debug!("graph assembled allele: {cis_dup_assembled:?}");
        cis_dups_assembled.push(cis_dup_assembled);
    }
    Ok((cis_dups_assembled, match_index_on_read))
}

/// Find the index of the first matching node on a read
/// and the haplotype
/// Note that this does not allow read to surpass the end of the haplotype
/// i.e. assuming haplotype ends with the ending node
/// `hap1` is a read, `hap2` is a haplotype.
/// Returns (match_index_on_read, match_index_on_hap)
/// # Arguments
/// * `hap1` - read
/// * `hap2` - haplotype
/// # Returns
/// * `Option<(usize, usize)>` - match index on read, match index on haplotype
pub fn match_read_allele_first_node_index(
    hap1: &Vec<i32>,
    hap2: &Vec<i32>,
) -> Option<(usize, usize)> {
    let mut hap2_mod = hap2.clone();
    if hap2_mod.starts_with(&[-10]) {
        hap2_mod.remove(0);
    }
    let hap2 = &hap2_mod;
    let hap1_len = hap1.len();
    let hap2_len = hap2.len();
    for i in 0..hap1_len {
        // below is not allowed
        // hap1/read      -------------
        // hap2/hapl ------------
        //          |--k--|offset|
        //
        // below is allowed
        // hap1/read       -------
        // hap2/hapl ---------------
        //           |--k--|offset|
        if i == 0 {
            for k in 0..hap2_len {
                let offset_index = cmp::min(hap2_len - k, hap1_len);
                if offset_index > 1 && k + hap1_len <= hap2_len {
                    // here we are requiring a match until the end of the haplotype (hap2)
                    let test_hap1 = &hap1[..offset_index];
                    let test_hap2 = &hap2[k..(k + offset_index)];
                    let (hap_match, mismatch) = compare_two_haps_same_length(test_hap1, test_hap2);
                    if hap_match.iter().sum::<i32>() >= 2 && mismatch == 0 {
                        return Some((i, k));
                    }
                }
            }
        } else {
            // below is allowed
            //           |-i-|
            // hap1/read  ----------
            // hap2/hapl     ---------------
            //               |offset|
            //
            // below is not allowed
            //           |-i-|
            // hap1/read  -------------
            // hap2/hapl     --------
            //               |offset|
            let offset_index = cmp::min(hap1_len - i, hap2_len);
            if offset_index > 1 && i + hap2_len >= hap1_len {
                // this was originally implemented in python
                //let test_hap1 = &hap1[(hap1_len-offset_index)..];
                let test_hap1 = &hap1[i..(i + offset_index)];
                let test_hap2 = &hap2[..offset_index];
                let (hap_match, mismatch) = compare_two_haps_same_length(test_hap1, test_hap2);
                if hap_match.iter().sum::<i32>() >= 2 && mismatch == 0 {
                    return Some((i, 0));
                }
            }
        }
    }
    None
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_match_read_allele_first_node_index() {
        let read = vec![4, 5, 6];
        let hap = vec![1, 2, 3, 0, 4, 5, 6];
        let res = match_read_allele_first_node_index(&read, &hap);
        assert!(!res.is_none());
        assert_eq!(res.unwrap(), (0, 4));

        let read = vec![7, 1, 2];
        let hap = vec![1, 2, 3, 0, 4, 5, 6];
        let res = match_read_allele_first_node_index(&read, &hap);
        assert!(!res.is_none());
        assert_eq!(res.unwrap(), (1, 0));

        let read = vec![4, 5, 6, 7];
        let hap = vec![1, 2, 3, 0, 4, 5, 6];
        let res = match_read_allele_first_node_index(&read, &hap);
        assert!(res.is_none());
    }
}
