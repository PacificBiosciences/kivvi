use crate::assembly::assembler::AssemblyResult;
use crate::repeat_unit::fingerprint::FingerprintInfo;
use crate::util::{DError, RegionCoordinates};
use itertools::Itertools;
use log::{debug, trace};
use regex::Regex;
use rust_htslib::faidx;
//use std::cmp;
use std::collections::{BTreeMap, HashMap, HashSet};
use std::path::PathBuf;

/// Represent a read in the data structure for plotting
#[derive(Clone, Debug, PartialEq)]
pub struct ReadInfoForPlotting {
    /// start position on the allele (first position without missing info)
    pub start_position: i64,
    /// bases at variant sites following start_position
    pub bases: Vec<usize>,
    /// whether read is nonunique
    pub is_nonuniq: bool,
}

/// Represent alleles for plotting
#[derive(Clone, Debug, PartialEq)]
pub struct AlleleInfoForPlotting {
    /// Each allele group is a `Vec<ReadInfoForPlotting>`. The first item of each allele group is the allele itself. Rest are reads.
    pub reads: Vec<Vec<ReadInfoForPlotting>>,
    /// number of variant sites used for plotting each copy
    pub variant_count_per_copy: i64,
}

/// Represent the fingerprints on an allele and their variants
#[derive(Clone, Debug)]
pub struct AlleleFingerprintVariantInfo {
    /// unique fp name (position on allele) -> original fp name (as index)
    pub fp_names: BTreeMap<String, i32>,
    /// unique fp name -> supporting reads
    pub suppporting_reads: BTreeMap<String, HashSet<String>>,
    /// unique fp name -> pos -> bases, unique reads only
    pub bases: BTreeMap<String, BTreeMap<(i32, i64), Vec<Vec<u8>>>>,
    /// unique fp name -> pos -> bases, including nonuniq reads
    pub bases_all: BTreeMap<String, BTreeMap<(i32, i64), Vec<Vec<u8>>>>,
}

/// Represent variant information on alleles
#[derive(Clone, Debug)]
pub struct VariantReport {
    /// variants on fingerprints on complete alleles. Alleles -> (unit index, fingerprint name) -> variants
    pub complete_allele_variants: BTreeMap<Vec<i32>, BTreeMap<(usize, i32), Vec<String>>>,
    /// variants on fingerprints not found on complete alleles
    pub fp_variants_on_incomplete_alleles: BTreeMap<i32, Vec<String>>,
    /// variant name -> variant information
    pub variant_summary: BTreeMap<String, Vec<VariantInfoByVariant>>,
    /// alleles -> (read name, starting site index on allele)
    pub reads_match_allele_index: BTreeMap<Vec<i32>, Vec<(String, i32)>>,
    /// data for plotting
    pub alleles_for_plot: Option<AlleleInfoForPlotting>,
    /// unique fp name -> supporting reads
    pub fp_suppporting_reads: Option<BTreeMap<String, HashSet<String>>>,
}

/// Variant at a position from a fingerprint (a set of reads)
#[derive(Clone, Debug)]
pub struct VariantInfoByFP {
    /// variant bases
    pub base: Option<String>,
    /// reference bases
    pub ref_base: String,
    /// total depth
    pub depth: usize,
    /// number of reads supporting the variant
    pub nread: usize,
    /// original consensus call format (with + or - for indels)
    pub original_base: Option<String>,
}

/// Report informatin about a variant
#[derive(Clone, Debug)]
pub struct VariantInfoByVariant {
    /// fingerprint that has the variant
    pub fingerprint: i32,
    /// fingerprint's position on allele
    pub fingerprint_in_allele: Option<String>,
    /// total depth
    pub depth: usize,
    /// supporting reads for the variant
    pub nread: usize,
}

/// Report variants on each KIV2 unit of each allele
/// # Arguments
/// * `fp_info` - fingerprint information
/// * `assembly_result` - assembly result from assembler
/// * `read_info` - read -> pos -> bases
/// * `reference` - reference file
/// * `region_coordinates` - coordinates defined for this region
/// * `predefined_variant_list` - variant list
/// # Returns
/// * `VariantReport` - variant report
pub fn report_variants(
    fp_info: FingerprintInfo,
    assembly_result: AssemblyResult,
    read_info: BTreeMap<String, BTreeMap<(i32, i64), Vec<u8>>>,
    reference: &PathBuf,
    region_coordinates: RegionCoordinates,
    predefined_variant_list: Option<Vec<String>>,
) -> Result<VariantReport, DError> {
    let alleles = assembly_result.complete;
    let supporting_reads = assembly_result.supporting_reads;
    let nonunique_reads = assembly_result.nonunique_reads.clone();
    let fp_info_clone = fp_info.clone();
    let read_edges = fp_info_clone.read_edges;
    let read_positions = fp_info_clone.read_positions;
    let grouped_reads = fp_info_clone.grouped_reads;
    let good_name_to_seq = fp_info_clone.good_name_to_seq;
    // ref
    let ref_reader = faidx::Reader::from_path(reference)?;

    let mut fps_on_incomplete_alleles: BTreeMap<i32, Vec<String>> = BTreeMap::new();
    let mut fps_on_complete_alleles = Vec::new();
    let mut complete_allele_variants: BTreeMap<Vec<i32>, BTreeMap<(usize, i32), Vec<String>>> =
        BTreeMap::new();
    // alleles -> (read name, starting site index on allele)
    let mut variant_summary: BTreeMap<String, Vec<VariantInfoByVariant>> = BTreeMap::new();
    let mut variant_name_old_format: BTreeMap<String, String> = BTreeMap::new();
    // position reads onto alleles
    let reads_match_allele_index =
        get_read_position_in_allele(read_edges.clone(), alleles, supporting_reads, false)?;
    debug!("reads_match_allele_index {reads_match_allele_index:?}");
    let mut fp_to_read: Option<BTreeMap<String, HashSet<String>>> = None;
    if !reads_match_allele_index.is_empty() {
        // for each fingerprint, get all reads and all bases
        let fp_bases = get_fp_bases(
            read_edges.clone(),
            read_positions,
            reads_match_allele_index.clone(),
            read_info.clone(),
            nonunique_reads,
        )?;
        fp_to_read = Some(fp_bases.suppporting_reads);
        let fp_names = fp_bases.fp_names;
        let fp_pileup = fp_bases.bases;
        let fp_pileup_all = fp_bases.bases_all;

        let mut allele_index = 0;
        // walk through each allele and get consensus bases at every position
        for allele in reads_match_allele_index.keys() {
            allele_index += 1;
            let allele_len = allele.len();
            for q in 1..(allele_len - 1) {
                let fp_name_on_allele = format!("{}.{}", allele_index, q);
                if fp_names.contains_key(&fp_name_on_allele) {
                    let original_fp_name = fp_names
                        .get(&fp_name_on_allele)
                        .ok_or("key not found: fp_name_on_allele in fp_names")?;
                    let this_fp_bases = fp_pileup.get(&fp_name_on_allele);
                    let this_fp_bases_all = fp_pileup_all.get(&fp_name_on_allele);
                    let mut fp_var = Vec::new();
                    if !this_fp_bases_all.is_none() {
                        let this_fp_bases_all = this_fp_bases_all.unwrap();
                        for (pos, site_all_bases) in this_fp_bases_all.iter() {
                            let tid = pos.0;
                            let ref_name = ref_reader.seq_name(tid as i32)?;
                            let ref_len = ref_reader.fetch_seq_len(&ref_name);
                            let ref_seq = ref_reader.fetch_seq(&ref_name, 0, ref_len as usize)?;
                            let this_offset =
                                region_coordinates.genome_offset.get(&ref_name).unwrap();
                            let this_offset = *this_offset as i64;
                            if !region_coordinates
                                .exclude_sites_vcf
                                .contains(&((*pos).1 + 1))
                                && (*pos).1 < ref_seq.len() as i64
                                && ((*pos).1 < region_coordinates.repeat_len as i64
                                    || region_coordinates.genome_offset.len() > 1)
                            // if multiple reference, then use all positions
                            // else, use only positions within the specified repeat length
                            {
                                let mut fp_base_consensus: VariantInfoByFP =
                                    get_consensus_var(site_all_bases.clone(), (*pos).1, &ref_seq)?;
                                // use unique reads if there are enough of them
                                if !this_fp_bases.is_none() {
                                    let this_fp_bases = this_fp_bases.unwrap();
                                    if this_fp_bases.contains_key(pos) {
                                        let site_all_bases_uniq = this_fp_bases
                                            .get(pos)
                                            .ok_or("position not found in this_fp_bases")?;
                                        if site_all_bases_uniq.len() >= 3 {
                                            fp_base_consensus = get_consensus_var(
                                                site_all_bases_uniq.clone(),
                                                (*pos).1,
                                                &ref_seq,
                                            )?;
                                        }
                                    }
                                }
                                if let Some(consensus) = fp_base_consensus.base {
                                    let variant_name = format!(
                                        "{}:{}>{}",
                                        (*pos).1 + 2 + this_offset,
                                        fp_base_consensus.ref_base,
                                        consensus
                                    );
                                    let original_base = fp_base_consensus.original_base.unwrap();
                                    let ref_base = vec![ref_seq[(*pos).1 as usize]];
                                    let ref_base_string =
                                        std::str::from_utf8(&ref_base)?.to_string();
                                    let variant_name_old = format!(
                                        "{}-{}:{}>{}",
                                        tid,
                                        (*pos).1,
                                        ref_base_string,
                                        original_base
                                    );
                                    variant_name_old_format
                                        .entry(variant_name_old.clone())
                                        .or_insert(variant_name.clone());

                                    fp_var.push(variant_name.clone());
                                    let this_var = VariantInfoByVariant {
                                        fingerprint: *original_fp_name,
                                        fingerprint_in_allele: Some(fp_name_on_allele.clone()),
                                        depth: fp_base_consensus.depth,
                                        nread: fp_base_consensus.nread,
                                    };
                                    variant_summary
                                        .entry(variant_name.clone())
                                        .or_default()
                                        .push(this_var);
                                }
                            }
                        }
                    } else {
                        debug!(
                        "this_fp_bases_all not found for fp {original_fp_name} fp_name_on_allele {}",
                        fp_name_on_allele.clone()
                    );
                    }
                    fp_var.sort_by(|a, b| {
                        let pos1 = a.split(":").nth(0).unwrap().parse::<i64>().unwrap();
                        let pos2 = b.split(":").nth(0).unwrap().parse::<i64>().unwrap();
                        pos1.cmp(&pos2)
                    });
                    fps_on_complete_alleles.push(*original_fp_name);
                    complete_allele_variants
                        .entry(allele.to_vec())
                        .or_default()
                        .insert((q, *original_fp_name), fp_var.clone());
                }
            }
        }
    }
    // for fingerprints that do not have unique supporting reads to place onto an allele
    // use all supporting reads for variant calling
    // this is particularly for fingerprints not in complete alleles
    for fp_name in good_name_to_seq.keys() {
        if !fps_on_complete_alleles.contains(fp_name) {
            let mut fp_var = Vec::new();
            let mut this_fp_pileup: BTreeMap<(i32, i64), Vec<Vec<u8>>> = BTreeMap::new();
            for (read, read_bases) in read_info.iter() {
                if grouped_reads.contains_key(read) {
                    if grouped_reads[read] == *fp_name {
                        for (pos, bases) in read_bases.iter() {
                            this_fp_pileup.entry(*pos).or_default().push(bases.to_vec());
                        }
                    }
                }
            }
            for (pos, fp_bases) in this_fp_pileup.iter() {
                let tid = pos.0;
                let ref_name = ref_reader.seq_name(tid as i32)?;
                let ref_len = ref_reader.fetch_seq_len(&ref_name);
                let ref_seq = ref_reader.fetch_seq(&ref_name, 0, ref_len as usize)?;
                let this_offset = region_coordinates.genome_offset.get(&ref_name).unwrap();
                let this_offset = *this_offset as i64;
                if !region_coordinates
                    .exclude_sites_vcf
                    .contains(&((*pos).1 + 1))
                    && (*pos).1 < ref_seq.len() as i64
                    && ((*pos).1 < region_coordinates.repeat_len as i64
                        || region_coordinates.genome_offset.len() > 1)
                {
                    let fp_base_consensus =
                        get_consensus_var(fp_bases.clone(), (*pos).1, &ref_seq)?;
                    if let Some(consensus) = fp_base_consensus.base {
                        let variant_name = format!(
                            "{}:{}>{}",
                            (*pos).1 + 2 + this_offset,
                            fp_base_consensus.ref_base,
                            consensus
                        );
                        let original_base = fp_base_consensus.original_base.unwrap();
                        let ref_base = vec![ref_seq[(*pos).1 as usize]];
                        let ref_base_string = std::str::from_utf8(&ref_base)?.to_string();
                        let variant_name_old =
                            format!("{}-{}:{}>{}", tid, (*pos).1, ref_base_string, original_base);
                        variant_name_old_format
                            .entry(variant_name_old.clone())
                            .or_insert(variant_name.clone());
                        fp_var.push(variant_name.clone());
                        let this_var = VariantInfoByVariant {
                            fingerprint: *fp_name,
                            fingerprint_in_allele: None,
                            depth: fp_base_consensus.depth,
                            nread: fp_base_consensus.nread,
                        };
                        variant_summary
                            .entry(variant_name.clone())
                            .or_default()
                            .push(this_var);
                    }
                }
            }
            fp_var.sort_by(|a, b| {
                let pos1 = a.split(":").nth(0).unwrap().parse::<i64>().unwrap();
                let pos2 = b.split(":").nth(0).unwrap().parse::<i64>().unwrap();
                pos1.cmp(&pos2)
            });
            fps_on_incomplete_alleles.entry(*fp_name).or_insert(fp_var);
        }
    }
    debug!("variant_name_old_format {:?}", variant_name_old_format);
    if reads_match_allele_index.is_empty() {
        return Ok(VariantReport {
            complete_allele_variants,
            fp_variants_on_incomplete_alleles: fps_on_incomplete_alleles,
            variant_summary,
            reads_match_allele_index,
            alleles_for_plot: None,
            fp_suppporting_reads: fp_to_read,
        });
    }
    // prepare for plotting
    let alleles_for_plot = make_data_for_alleles(
        &fp_info,
        &complete_allele_variants,
        &fps_on_incomplete_alleles,
        &reads_match_allele_index,
        &read_info,
        &variant_name_old_format,
        &assembly_result.nonunique_reads,
        &ref_reader,
        predefined_variant_list,
    )?;
    Ok(VariantReport {
        complete_allele_variants,
        fp_variants_on_incomplete_alleles: fps_on_incomplete_alleles,
        variant_summary,
        reads_match_allele_index,
        alleles_for_plot: Some(alleles_for_plot),
        fp_suppporting_reads: fp_to_read,
    })
}

/// Get consensus among a set of bases
/// # Arguments
/// * `bases` - all bases at this position, each read at this position is a vec<u8>
/// * `pos` - position on repeat unit
/// * `ref_seq` - reference sequence
/// # Returns
/// * `VariantInfoByFP` - variant at this position from this fingerprint (this set of reads)
fn get_consensus_var(
    bases: Vec<Vec<u8>>,
    pos: i64,
    ref_seq: &[u8],
) -> Result<VariantInfoByFP, DError> {
    let depth = bases.len();
    let ref_base = vec![ref_seq[pos as usize]];
    let ref_base_string = std::str::from_utf8(&ref_base)?.to_string();
    let none_var = VariantInfoByFP {
        base: None,
        ref_base: ref_base_string.clone(),
        depth,
        nread: 0,
        original_base: None,
    };

    // count number of unique bases
    let mut bases_count: HashMap<Vec<u8>, usize> = HashMap::new();
    for base in bases {
        *bases_count.entry(base).or_default() += 1;
    }
    let mut bases_count2 = bases_count
        .iter()
        .map(|(x, y)| (x.clone(), *y))
        .collect::<Vec<_>>();
    // Make equal-depth consensus selection deterministic by falling back to a
    // stable lexical ordering of the observed allele strings.
    bases_count2.sort_by(|a, b| b.1.cmp(&a.1).then_with(|| a.0.cmp(&b.0)));
    if bases_count2.is_empty() {
        return Ok(none_var.clone());
    }
    if bases_count2.len() > 1 {
        let count1 = bases_count2[0].1;
        let count2 = bases_count2[1].1;
        if count1 == count2 && (bases_count2[0].0 == ref_base || bases_count2[1].0 == ref_base) {
            return Ok(none_var.clone());
        }
    }
    let consensus_count = bases_count2[0].1;
    let fp_base_consensus = bases_count2[0].0.clone();
    if fp_base_consensus != ref_base && !fp_base_consensus.contains(&b'*') {
        if consensus_count < 2 {
            return Ok(none_var.clone());
        }
        let fp_base_consensus_string = std::str::from_utf8(&fp_base_consensus)?.to_string();
        // check neighboring bases on ref for indels
        if fp_base_consensus.len() > 1 && fp_base_consensus.contains(&b'-') {
            let re = Regex::new(r"\-\d+")?;
            let del_seq = re
                .split(&fp_base_consensus_string)
                .collect::<Vec<&str>>()
                .last()
                .ok_or("cannot find last")?
                .to_string();
            /*
            let ref_len = ref_seq.len() as i64;
            let del_len: i64 = del_seq.len() as i64;
            let max_pos = cmp::min(pos + del_len + 3, ref_len) as usize;
            let nstart = (pos + del_len + 1) as usize;
            let ref_bases_after_del = std::str::from_utf8(&ref_seq[nstart..max_pos])?.to_string();
            //debug!("pos {pos} {fp_base_consensus_string:?} del_seq {del_seq:?}, ref_bases_after_del {ref_bases_after_del:?}");
            let s1: HashSet<char> = del_seq.chars().collect();
            let s2: HashSet<char> = ref_bases_after_del.to_string().chars().collect();
            if s1.len() == 1 && s2.len() == 1 && s1 == s2 {
                return Ok(none_var.clone());
            }
            */
            let mut new_ref_base = ref_base.clone();
            new_ref_base.extend_from_slice(del_seq.as_bytes());
            return Ok(VariantInfoByFP {
                base: Some(ref_base_string.clone()),
                ref_base: std::str::from_utf8(&new_ref_base)?.to_string(),
                depth,
                nread: consensus_count,
                original_base: Some(fp_base_consensus_string),
            });
        }
        if fp_base_consensus.len() > 1 && fp_base_consensus.contains(&b'+') {
            let re = Regex::new(r"\+\d+")?;
            let ins_seq = re
                .split(&fp_base_consensus_string)
                .collect::<Vec<&str>>()
                .last()
                .ok_or("cannot find last")?
                .to_string();

            /*
            let max_pos = cmp::min(pos + 3, ref_len) as usize;
            let nstart = (pos + 1) as usize;
            let ref_bases_after_ins = std::str::from_utf8(&ref_seq[nstart..max_pos])?.to_string();
            //debug!("pos {pos} {fp_base_consensus_string:?} ins_seq {ins_seq:?}, ref_bases_after_ins {ref_bases_after_ins:?}");
            let s1: HashSet<char> = ins_seq.chars().sorted().collect();
            let s2: HashSet<char> = ref_bases_after_ins.to_string().chars().sorted().collect();
            if s1.len() == 1 && s2.len() == 1 && s1 == s2 {
                return Ok(none_var.clone());
            };
            */
            let mut variant_base = ref_base.clone();
            variant_base.extend_from_slice(ins_seq.as_bytes());
            return Ok(VariantInfoByFP {
                base: Some(std::str::from_utf8(&variant_base)?.to_string()),
                ref_base: ref_base_string.clone(),
                depth,
                nread: consensus_count,
                original_base: Some(fp_base_consensus_string),
            });
        }
        if fp_base_consensus_string != String::from("x")
            && fp_base_consensus_string != String::from("-")
        {
            return Ok(VariantInfoByFP {
                base: Some(fp_base_consensus_string.clone()),
                ref_base: ref_base_string.clone(),
                depth,
                nread: consensus_count,
                original_base: Some(fp_base_consensus_string.clone()),
            });
        }
    }
    Ok(none_var.clone())
}

/// Get the supporting reads and all bases of each fingerprint on each allele
/// # Arguments
/// * `read_edges` - read name -> vec of fingerprints
/// * `read_positions` - read name -> vec of starting positions of each fingerprint
/// * `reads_match_allele_index` - alleles -> (read name, starting site index on allele)
/// * `read_info` - read -> pos -> bases
/// * `nonunique_reads` - vector of nonunique read names
/// # Returns
/// * `AlleleFingerprintVariantInfo` - fingerprints on each allele and their variants
fn get_fp_bases(
    read_edges: BTreeMap<String, Vec<i32>>,
    read_positions: BTreeMap<String, Vec<i32>>,
    reads_match_allele_index: BTreeMap<Vec<i32>, Vec<(String, i32)>>,
    read_info: BTreeMap<String, BTreeMap<(i32, i64), Vec<u8>>>,
    nonunique_reads: Vec<String>,
) -> Result<AlleleFingerprintVariantInfo, DError> {
    let mut fp_names: BTreeMap<String, i32> = BTreeMap::new();
    let mut suppporting_reads: BTreeMap<String, HashSet<String>> = BTreeMap::new();
    let mut bases: BTreeMap<String, BTreeMap<(i32, i64), Vec<Vec<u8>>>> = BTreeMap::new();
    let mut bases_all: BTreeMap<String, BTreeMap<(i32, i64), Vec<Vec<u8>>>> = BTreeMap::new();
    let mut allele_index = 0;
    for (allele, read_fp_matches) in reads_match_allele_index.iter() {
        allele_index += 1;
        for (read, fp_index) in read_fp_matches {
            if read_positions.contains_key(read) {
                let this_read_positions = read_positions[read].clone();
                let read_nodes = read_edges
                    .get(read)
                    .ok_or("key not found: read in read_edges")?;
                assert!(this_read_positions.len() == read_nodes.len());
                for (j, this_fp) in read_nodes.iter().enumerate() {
                    let fp_index_on_allele = j + *fp_index as usize;
                    if *this_fp > 0
                        && fp_index_on_allele > 0
                        && fp_index_on_allele < allele.len() - 1
                    {
                        let uniq_fp_name = format!("{}.{}", allele_index, fp_index_on_allele);
                        fp_names.entry(uniq_fp_name.clone()).or_insert(*this_fp);

                        let this_fp_start_pos = this_read_positions[j];
                        let read_new_name = format!("{read}:{}", this_fp_start_pos);
                        suppporting_reads
                            .entry(uniq_fp_name.clone())
                            .or_default()
                            .insert(read_new_name.clone());
                        if read_info.contains_key(&read_new_name) {
                            let read_bases = read_info.get(&read_new_name).ok_or(format!(
                                "key not found: read_new_name {read_new_name:?} in read_info"
                            ))?;
                            if !nonunique_reads.contains(read) {
                                for (pos, base) in read_bases.iter() {
                                    bases
                                        .entry(uniq_fp_name.clone())
                                        .or_default()
                                        .entry(*pos)
                                        .or_default()
                                        .push(base.to_vec());
                                }
                            }
                            for (pos, base) in read_bases.iter() {
                                bases_all
                                    .entry(uniq_fp_name.clone())
                                    .or_default()
                                    .entry(*pos)
                                    .or_default()
                                    .push(base.to_vec());
                            }
                        }
                    }
                }
            }
        }
    }
    Ok(AlleleFingerprintVariantInfo {
        fp_names,
        suppporting_reads,
        bases,
        bases_all,
    })
}

/// Find the index of each read on each allele
/// # Arguments
/// * `read_edges` - read name -> vec of fingerprints
/// * `alleles` - complete alleles
/// * `supporting_reads` - allele -> reads
/// * `allow_incomplete_allele` - whether to allow incomplete alleles
/// # Returns
/// * `reads_match_allele_index` - alleles -> (read name, starting index on allele)
pub fn get_read_position_in_allele(
    read_edges: BTreeMap<String, Vec<i32>>,
    alleles: Vec<Vec<i32>>,
    supporting_reads: BTreeMap<Vec<i32>, HashSet<String>>,
    allow_incomplete_allele: bool,
) -> Result<BTreeMap<Vec<i32>, Vec<(String, i32)>>, DError> {
    let mut reads_match_allele_index: BTreeMap<Vec<i32>, Vec<(String, i32)>> = BTreeMap::new();
    for allele in alleles {
        let allele_len = allele.len();
        if supporting_reads.contains_key(&allele) {
            let allele_reads = supporting_reads
                .get(&allele)
                .ok_or("key not found: allele in supporting_reads")?;
            for read in allele_reads {
                if read_edges.contains_key(read) {
                    let read_nodes = read_edges
                        .get(read)
                        .ok_or("key not found: read in read_edges")?;
                    let read_nodes_len = read_nodes.len();
                    let mut found_match = false;
                    // allele incomplete. read can start earlier than allele
                    if allow_incomplete_allele {
                        let mut i_index: i32 = 0;
                        for i in 0..read_nodes_len {
                            i_index = i as i32;
                            if read_nodes_len - i >= 2 && read_nodes_len - i <= allele_len {
                                let nodes_in_read = &read_nodes[i..];
                                let nodes_in_allele = &allele[..(read_nodes_len - i)];
                                let mut match_allele = Vec::new();
                                for (j, read_node) in nodes_in_read.iter().enumerate() {
                                    if *read_node == 0 {
                                        match_allele.push(0)
                                    } else if *read_node == nodes_in_allele[j] {
                                        match_allele.push(1);
                                    } else {
                                        match_allele.push(-1);
                                    }
                                }
                                let match_count = match_allele
                                    .iter()
                                    .filter(|x| **x == 1)
                                    .collect::<Vec<_>>()
                                    .len();
                                if !match_allele.contains(&(-1)) && match_count > 1 {
                                    found_match = true;
                                    break;
                                }
                            }
                        }
                        if found_match {
                            reads_match_allele_index
                                .entry(allele.clone())
                                .or_default()
                                .push((read.to_string(), 0 - i_index));
                        }
                    }
                    if !found_match && allele_len >= read_nodes_len {
                        let mut i_index: i32 = 0;
                        for i in 0..(allele_len - read_nodes_len + 1) {
                            i_index = i as i32;
                            let nodes_in_allele = &allele[i..(i + read_nodes_len)];
                            let mut match_allele = Vec::new();
                            for (j, read_node) in read_nodes.iter().enumerate() {
                                if *read_node == 0 {
                                    match_allele.push(0)
                                } else if *read_node == nodes_in_allele[j] {
                                    match_allele.push(1);
                                } else {
                                    match_allele.push(-1);
                                }
                            }
                            let match_count = match_allele
                                .iter()
                                .filter(|x| **x == 1)
                                .collect::<Vec<_>>()
                                .len();
                            if !match_allele.contains(&(-1)) && match_count > 1 {
                                found_match = true;
                                break;
                            }
                        }
                        if found_match {
                            reads_match_allele_index
                                .entry(allele.clone())
                                .or_default()
                                .push((read.to_string(), i_index));
                        }
                    }
                }
            }
        }
    }
    Ok(reads_match_allele_index)
}

/// Prepare data for plotting
/// # Arguments
/// * `fp_info` - fingerprint information
/// * `complete_allele_variants` - variants on complete alleles. allele -> (index on allele, fp id) -> vector of variants
/// * `fps_on_incomplete_alleles` - variants on incomplete alleles. fp id -> vector of variants
/// * `reads_match_allele_index` - alleles -> (read name, starting index on allele)
/// * `read_info` - read -> pos -> bases
/// * `variant_name_old_format` - variant name: tid-pos:A>T (old)-> pos:A>T (new)
/// * `nonunique_reads` - vector of nonuniq read names
/// * `ref_reader` - reference reader
/// * `variant_list` - variant list
/// # Colors
/// * 0 -> reference, yellow
/// * 1 -> variant, black
/// * 2 -> missing info, pink
/// * 3 -> no data (reads not overlapping), white. Note: not used at this step yet, introduced later when plotting
/// * 4 -> flank, orange
/// # Returns
/// * `AlleleInfoForPlotting` - data for plotting
pub fn make_data_for_alleles(
    fp_info: &FingerprintInfo,
    complete_allele_variants: &BTreeMap<Vec<i32>, BTreeMap<(usize, i32), Vec<String>>>,
    fps_on_incomplete_alleles: &BTreeMap<i32, Vec<String>>,
    reads_match_allele_index: &BTreeMap<Vec<i32>, Vec<(String, i32)>>,
    read_info: &BTreeMap<String, BTreeMap<(i32, i64), Vec<u8>>>,
    variant_name_old_format: &BTreeMap<String, String>,
    nonunique_reads: &Vec<String>,
    ref_reader: &faidx::Reader,
    variant_list: Option<Vec<String>>,
) -> Result<AlleleInfoForPlotting, DError> {
    let mut out_vec = Vec::new();
    let read_edges = &fp_info.read_edges;
    let read_positions = &fp_info.read_positions;
    let alleles = reads_match_allele_index.keys();
    let mut all_var_original = Vec::new();
    if !variant_list.is_none() {
        all_var_original = variant_list.unwrap();
    } else {
        for variant in variant_name_old_format.keys() {
            if !all_var_original.contains(variant) {
                all_var_original.push(variant.to_string());
            }
        }
    }
    let mut all_var_pos_original = Vec::new();
    let mut indel_pos = Vec::new();
    for var in &all_var_original {
        let fields = var.split_terminator(':').collect::<Vec<_>>();
        //let var_pos = fields.first().unwrap().parse::<i64>().unwrap();
        let var_pos = fields
            .first()
            .unwrap()
            .split_terminator('-')
            .collect::<Vec<_>>()
            .last()
            .unwrap()
            .parse::<i64>()?;
        let fields2 = fields
            .last()
            .unwrap()
            .split_terminator('>')
            .collect::<Vec<_>>();
        let var_alt = fields2.last().unwrap();
        if var_alt.contains('-') {
            // handle positions in a deletion
            let var_len = var_alt
                .split('-')
                .nth(1)
                .unwrap()
                .chars()
                .take_while(|c| c.is_ascii_digit())
                .collect::<String>()
                .parse::<i64>()
                .unwrap();
            for pos in var_pos..(var_pos + var_len + 1) {
                all_var_pos_original.push(pos);
            }
        } else {
            all_var_pos_original.push(var_pos);
        }
        if var_alt.contains('-') || var_alt.contains('+') {
            indel_pos.push(var_pos);
        }
    }
    let mut all_var_pos = Vec::new();
    for var in &all_var_original {
        let fields = var.split_terminator(':').collect::<Vec<_>>();
        let pos = fields
            .first()
            .unwrap()
            .split_terminator('-')
            .collect::<Vec<_>>()
            .last()
            .unwrap()
            .parse::<i64>()?;
        let pos_count = all_var_pos_original.iter().filter(|x| **x == pos).count();
        // take only positions with only one variant
        if pos_count == 1 && !indel_pos.contains(&pos) {
            all_var_pos.push(pos);
        }
    }
    all_var_pos.sort();
    let all_var_sorted = all_var_original
        .iter()
        .filter(|x| {
            let fields = x.split_terminator(':').collect::<Vec<_>>();
            let pos = fields
                .first()
                .unwrap()
                .split_terminator('-')
                .collect::<Vec<_>>()
                .last()
                .unwrap()
                .parse::<i64>()
                .unwrap();
            all_var_pos.contains(&pos)
        })
        .map(|x| {
            let fields = x.split_terminator(':').collect::<Vec<_>>();
            let pos = fields
                .first()
                .unwrap()
                .split_terminator('-')
                .collect::<Vec<_>>()
                .last()
                .unwrap()
                .parse::<i64>()
                .unwrap();
            (x, pos)
        })
        .sorted_by(|a, b| a.1.cmp(&b.1))
        .map(|(x, _y)| x.clone())
        .collect::<Vec<String>>();

    debug!("all_var_sorted {all_var_sorted:?}");

    let nvar = all_var_sorted.len() as i64;
    let mut allele_lens = Vec::new();
    for allele in alleles.clone() {
        let allele_len = allele.len() - 2;
        allele_lens.push(allele_len);
    }
    let allele_len_max = *allele_lens
        .iter()
        .max()
        .ok_or("cannot find max of allele_lens")? as i64;
    // go through each allele
    for allele in alleles.clone() {
        let this_allele_variants = complete_allele_variants
            .get(allele)
            .ok_or("key not found: allele in complete_allele_variants")?;
        let mut out_vec_allele = Vec::new();
        let allele_len = allele.len() - 2;
        let allele2 = allele.to_vec();
        // first the allele itself
        let mut allele_line: Vec<usize> = Vec::new();
        // add left flank
        allele_line.push(4);
        for repeat_index in 0..allele_len_max {
            if repeat_index < allele_len as i64 {
                let mut fp_var = &vec![];
                let node_name = allele2[repeat_index as usize + 1];
                let node_on_allele = (repeat_index as usize + 1, node_name);
                if this_allele_variants.contains_key(&node_on_allele) {
                    fp_var = this_allele_variants.get(&node_on_allele).unwrap();
                } else {
                    if fps_on_incomplete_alleles.contains_key(&node_name) {
                        fp_var = fps_on_incomplete_alleles.get(&node_name).unwrap();
                    }
                }
                for var in &all_var_sorted {
                    let new_var_name = variant_name_old_format.get(var).unwrap();
                    if fp_var.contains(new_var_name) {
                        allele_line.push(1);
                    } else {
                        allele_line.push(0);
                    }
                }
            }
        }
        // add right flank
        allele_line.push(4);
        out_vec_allele.push(ReadInfoForPlotting {
            start_position: 0,
            bases: allele_line,
            is_nonuniq: false,
        });
        // go through each read to get bases at each position for this allele
        let mut read_lines = Vec::new();
        for (read, read_index_on_allele) in reads_match_allele_index[allele].iter() {
            let mut read_line = Vec::new();
            // before read start
            let n_copy_before: i64 = if *read_index_on_allele < 1 {
                0
            } else {
                *read_index_on_allele as i64 - 1
            };

            let this_read_nodes = read_edges[read].clone();
            let this_read_positions = read_positions[read].clone();
            let mut read_map_index = 0;
            for (read_node, read_position) in this_read_nodes.iter().zip(this_read_positions.iter())
            {
                read_map_index += 1;
                if *read_index_on_allele + read_map_index - 1 > 0
                    && *read_index_on_allele + read_map_index - 1 < allele_len as i32 + 1
                    && *read_node >= 0
                {
                    let new_read_name = format!("{}:{}", read, *read_position);
                    if read_info.contains_key(&new_read_name) {
                        let this_read_bases = read_info
                            .get(&new_read_name)
                            .ok_or("key not found: new_read_name in read_info")?;
                        for (variant_index, pos) in all_var_pos.iter().enumerate() {
                            let expected_variant_old = &all_var_sorted[variant_index];
                            let expected_variant = expected_variant_old
                                .split_terminator('>')
                                .collect::<Vec<_>>()
                                .last()
                                .ok_or("last not found in expected_variant_new")?
                                .parse::<String>()?;
                            let expected_ref = expected_variant_old
                                .split_terminator(':')
                                .collect::<Vec<_>>()
                                .last()
                                .unwrap()
                                .split_terminator('>')
                                .collect::<Vec<_>>()
                                .first()
                                .unwrap()
                                .parse::<String>()?;
                            let fields = expected_variant_old
                                .split_terminator(':')
                                .collect::<Vec<_>>();
                            let variant_tid = fields
                                .first()
                                .unwrap()
                                .split_terminator('-')
                                .collect::<Vec<_>>()
                                .first()
                                .unwrap()
                                .parse::<i32>()?;
                            let tid = this_read_bases.first_key_value().unwrap().0 .0;
                            let ref_name = ref_reader.seq_name(tid as i32)?;
                            let ref_len = ref_reader.fetch_seq_len(&ref_name);
                            if *pos >= ref_len as i64 {
                                debug!("read {new_read_name:?} node {read_node} {expected_variant_old:?} pos {pos} is out of range for ref {ref_name}");
                                read_line.push(0);
                            } else if !this_read_bases.contains_key(&(tid, *pos)) {
                                trace!("read {new_read_name:?} node {read_node} {expected_variant_old:?}, ({tid}, {pos}) not in this_read_bases");
                                read_line.push(2);
                            } else if variant_tid == tid {
                                let this_read_base = this_read_bases
                                    .get(&(tid, *pos))
                                    .ok_or("index not found in this_read_bases")?;
                                let this_read_base_string =
                                    std::str::from_utf8(this_read_base)?.to_string();
                                if this_read_base_string == expected_variant {
                                    read_line.push(1);
                                } else if this_read_base_string == expected_ref {
                                    read_line.push(0);
                                } else {
                                    trace!("read {new_read_name:?} node {read_node} pos {pos} {expected_variant_old:?} expected_variant {expected_variant:?} this_read_base_string {this_read_base_string:?}");
                                    read_line.push(2);
                                }
                            } else {
                                // variant is on a different tid
                                read_line.push(0);
                            }
                        }
                    } else if *read_node == 0 {
                        debug!("read {new_read_name:?} node {read_node} not in read_info");
                        for (_variant_index, _pos) in all_var_pos.iter().enumerate() {
                            read_line.push(2);
                        }
                    }
                }
            }
            assert_eq!(read_map_index, this_read_nodes.len() as i32);
            let mut beginning_unknown: usize = 0;
            for a in &read_line {
                if *a == 2 {
                    beginning_unknown += 1;
                } else {
                    break;
                }
            }
            let start_position = nvar * n_copy_before + beginning_unknown as i64;
            let mut read_line_new = read_line[beginning_unknown..].to_vec();
            // add left and right flanks
            if start_position == 0 {
                let first_node = this_read_nodes
                    .first()
                    .ok_or("first not found in this_read_nodes")?;
                if *first_node < 0 && *first_node > -10 {
                    read_line_new.insert(0, 4);
                }
            }
            let last_node = this_read_nodes
                .last()
                .ok_or("last not found in this_read_nodes")?;
            if *last_node <= -10 {
                read_line_new.push(4);
            } else {
                while *read_line_new
                    .last()
                    .ok_or("last not found in read_line_new")?
                    == 2
                {
                    read_line_new.pop();
                }
                //let total_expect_len = nvar * (allele_len as i64);
                //let read_end = nvar * n_copy_before + read_line_new.len() as i64;
                //while read_end > total_expect_len {
                //    read_line.pop();
                //}
            }
            read_lines.push(ReadInfoForPlotting {
                start_position,
                bases: read_line_new,
                is_nonuniq: nonunique_reads.contains(read),
            });
        }
        read_lines.sort_by(|a, b| a.start_position.cmp(&b.start_position));
        for a in read_lines {
            if a.bases.contains(&4) || a.bases.len() > nvar as usize {
                // other reads not overlapping left flank need to shift right by one
                let new_start = if a.bases.starts_with(&[4]) {
                    a.start_position
                } else {
                    a.start_position + 1
                };
                out_vec_allele.push(ReadInfoForPlotting {
                    start_position: new_start,
                    bases: a.bases,
                    is_nonuniq: a.is_nonuniq,
                });
            }
        }
        out_vec.push(out_vec_allele);
    }
    Ok(AlleleInfoForPlotting {
        reads: out_vec,
        variant_count_per_copy: nvar,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_get_read_position_in_allele() {
        let alleles = vec![vec![1, 2, 3], vec![4, 5, 3]];
        let mut read_edges = BTreeMap::new();
        read_edges
            .entry(String::from("read1"))
            .or_insert(vec![1, 2]);
        read_edges
            .entry(String::from("read2"))
            .or_insert(vec![2, 3]);
        read_edges
            .entry(String::from("read3"))
            .or_insert(vec![1, 3]);
        let mut supporting_reads: BTreeMap<Vec<i32>, HashSet<String>> = BTreeMap::new();
        let mut set1: HashSet<String> = HashSet::new();
        set1.insert(String::from("read1"));
        set1.insert(String::from("read2"));
        supporting_reads.entry(vec![1, 2, 3]).or_insert(set1);
        let reads_match_allele_index =
            get_read_position_in_allele(read_edges, alleles, supporting_reads, false).unwrap();
        println!("reads_match_allele_index {:?}", reads_match_allele_index);
        assert!(reads_match_allele_index.contains_key(&vec![1, 2, 3]));
        let hap1_reads = reads_match_allele_index[&vec![1, 2, 3]].clone();
        assert!(hap1_reads.contains(&(String::from("read1"), 0)));
        assert!(hap1_reads.contains(&(String::from("read2"), 1)));
        assert!(!reads_match_allele_index.contains_key(&vec![4, 5, 3]));

        let alleles = vec![vec![1, 2, 3, 7], vec![4, 5, 3]];
        let mut read_edges = BTreeMap::new();
        read_edges
            .entry(String::from("read1"))
            .or_insert(vec![6, 1, 2]);
        read_edges
            .entry(String::from("read2"))
            .or_insert(vec![8, 6, 1, 2, 3]);
        read_edges
            .entry(String::from("read3"))
            .or_insert(vec![1, 2]);
        let mut supporting_reads: BTreeMap<Vec<i32>, HashSet<String>> = BTreeMap::new();
        let mut set1: HashSet<String> = HashSet::new();
        set1.insert(String::from("read1"));
        set1.insert(String::from("read2"));
        set1.insert(String::from("read3"));
        supporting_reads.entry(vec![1, 2, 3, 7]).or_insert(set1);
        let reads_match_allele_index =
            get_read_position_in_allele(read_edges, alleles, supporting_reads, true).unwrap();
        assert!(reads_match_allele_index.contains_key(&vec![1, 2, 3, 7]));
        let hap1_reads = reads_match_allele_index[&vec![1, 2, 3, 7]].clone();
        assert!(hap1_reads.contains(&(String::from("read1"), -1)));
        assert!(hap1_reads.contains(&(String::from("read2"), -2)));
        assert!(hap1_reads.contains(&(String::from("read3"), 0)));
        assert!(!reads_match_allele_index.contains_key(&vec![4, 5, 3]));
    }

    #[test]
    fn test_get_consensus_var() {
        let bases = vec![
            String::from("A").as_bytes().to_vec(),
            String::from("A").as_bytes().to_vec(),
            String::from("C").as_bytes().to_vec(),
        ];
        let refseq = &[b'T', b'A', b'T'];
        let consensus_var = get_consensus_var(bases, 0, refseq).unwrap();
        assert_eq!(consensus_var.base, Some(String::from("A")));

        let bases = vec![
            String::from("A").as_bytes().to_vec(),
            String::from("C").as_bytes().to_vec(),
        ];
        let refseq = &[b'T', b'A', b'T'];
        let consensus_var = get_consensus_var(bases, 0, refseq).unwrap();
        assert_eq!(consensus_var.base, None);
        assert_eq!(consensus_var.nread, 0);

        let bases = vec![
            String::from("A").as_bytes().to_vec(),
            String::from("A").as_bytes().to_vec(),
            String::from("C").as_bytes().to_vec(),
        ];
        let refseq = &[b'T', b'A', b'T'];
        let consensus_var = get_consensus_var(bases, 1, refseq).unwrap();
        assert_eq!(consensus_var.base, None);

        let bases = vec![
            String::from("A").as_bytes().to_vec(),
            String::from("A").as_bytes().to_vec(),
            String::from("C").as_bytes().to_vec(),
            String::from("C").as_bytes().to_vec(),
        ];
        let refseq = &[b'T', b'A', b'T'];
        let consensus_var = get_consensus_var(bases, 1, refseq).unwrap();
        assert_eq!(consensus_var.base, None);

        // indel, unfiltered
        let bases = vec![
            String::from("A+1T").as_bytes().to_vec(),
            String::from("A+1T").as_bytes().to_vec(),
            String::from("A+1T").as_bytes().to_vec(),
            String::from("C").as_bytes().to_vec(),
        ];
        let refseq = &[b'T', b'A', b'T', b'C'];
        let consensus_var = get_consensus_var(bases, 1, refseq).unwrap();
        assert_eq!(consensus_var.base, Some(String::from("AT")));
        assert_eq!(consensus_var.nread, 3);
        assert_eq!(consensus_var.depth, 4);

        // indel, unfiltered
        let bases = vec![
            String::from("A-1T").as_bytes().to_vec(),
            String::from("A-1T").as_bytes().to_vec(),
            String::from("A-1T").as_bytes().to_vec(),
            String::from("C").as_bytes().to_vec(),
        ];
        let refseq = &[b'T', b'A', b'T', b'C'];
        let consensus_var = get_consensus_var(bases, 1, refseq).unwrap();
        assert_eq!(consensus_var.base, Some(String::from("A")));
        assert_eq!(consensus_var.ref_base, String::from("AT"));

        /*
        // indel, filtered
        let bases = vec![
            String::from("A+1T").as_bytes().to_vec(),
            String::from("A+1T").as_bytes().to_vec(),
            String::from("A+1T").as_bytes().to_vec(),
            String::from("C").as_bytes().to_vec(),
        ];
        let refseq = &[b'T', b'A', b'T', b'T'];
        let consensus_var = get_consensus_var(bases, 1, refseq).unwrap();
        assert_eq!(consensus_var.base, None);

        // indel, filtered
        let bases = vec![
            String::from("A-1T").as_bytes().to_vec(),
            String::from("A-1T").as_bytes().to_vec(),
            String::from("A-1T").as_bytes().to_vec(),
            String::from("C").as_bytes().to_vec(),
        ];
        let refseq = &[b'T', b'A', b'T', b'T'];
        let consensus_var = get_consensus_var(bases, 1, refseq).unwrap();
        assert_eq!(consensus_var.base, None);
        */
    }
}
