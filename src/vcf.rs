use crate::util::{DError, DResult, RegionCoordinates};
use crate::variant::VariantInfoByVariant;
use itertools::Itertools;
use rust_htslib::bcf::{self, record::GenotypeAllele, Format};
use std::collections::{BTreeMap, HashSet};
use std::env;
use std::path::PathBuf;

/// Header lines defining the INFO and FORMAT fields for the VCF file.
const VCF_LINES: [&str; 3] = [
    r#"##FILTER=<ID=PASS,Description="All filters passed">"#,
    r#"##INFO=<ID=RU,Number=.,Type=String,Description="Repeat unit that the variant is in. The four values for each repeat unit are repeat unit ID, repeat unit position on allele, read depth and number of reads supporting the variant">"#,
    r#"##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">"#,
];

/// Write variants to VCF
/// # Arguments
/// * `output_path` - output VCF
/// * `_sample_name` - sample ID for column name (not used)
/// * `variant_summary` - variant name -> variant information on each fingerprint that has it
/// * `region_coordinates` - region coordinates
/// * `allele_len` - number of alleles
pub fn write_vcf(
    output_path: &PathBuf,
    _sample_name: &str,
    variant_summary: BTreeMap<String, Vec<VariantInfoByVariant>>,
    region_coordinates: RegionCoordinates,
    allele_len: usize,
) -> DResult {
    // sort variant_summary
    let variant_summary_pos = variant_summary
        .keys()
        .map(|x| {
            x.split_terminator(':')
                .collect::<Vec<_>>()
                .first()
                .unwrap()
                .parse::<i64>()
                .unwrap()
        })
        .collect::<Vec<i64>>();
    let all_var_sorted = variant_summary
        .keys()
        .zip(variant_summary_pos.iter())
        .sorted_by(|a, b| a.1.cmp(b.1))
        .map(|(x, _y)| x.clone())
        .collect::<Vec<String>>();

    let mut vcf_header = bcf::header::Header::new();
    // add header
    for line in VCF_LINES.iter() {
        vcf_header.push_record(line.as_bytes());
    }

    let contig_line = format!(
        r#"##contig=<ID={},length={}>"#,
        region_coordinates.chromosome_output, region_coordinates.chromosome_len
    );
    vcf_header.push_record(contig_line.as_bytes());

    let args: Vec<String> = env::args().collect();
    let command_line = args.join(" ");
    let line = format!("##{}Command={}", env!("CARGO_PKG_NAME"), command_line);
    vcf_header.push_record(line.as_bytes());
    //vcf_header.push_sample(sample_name.as_bytes());
    for i in 0..allele_len {
        let allele_name = format!("allele{}", i + 1);
        vcf_header.push_sample(allele_name.as_bytes());
    }

    let mut writer = bcf::Writer::from_path(output_path, &vcf_header, true, Format::Vcf)
        .map_err(|_| format!("Invalid VCF output path: {}", output_path.display()))?;

    for variant in all_var_sorted {
        let variant_info = variant_summary.get(&variant).ok_or("key not found")?;
        let mut record = writer.empty_record();

        let contig = region_coordinates.chromosome_output.as_bytes();
        let rid = writer.header().name2rid(contig)?;
        record.set_rid(Some(rid));

        let variant_pos = variant
            .clone()
            .split_terminator(':')
            .collect::<Vec<_>>()
            .first()
            .ok_or("first not found")?
            .parse::<i64>()?;
        record.set_pos(variant_pos - 1);

        // variant quality?
        record.set_qual(f32::from_bits(0x7F800001));
        //record.push_info_integer(b"RPOS", &[variant_pos as i32])?;

        let (data, alleles_per_variant) = encode_ru_field(variant_info.to_vec())?;
        record.push_info_string(b"RU", &[data.as_bytes()])?;

        let ref_base = variant
            .clone()
            .split_terminator(':')
            .collect::<Vec<_>>()
            .last()
            .ok_or("last not found")?
            .split_terminator('>')
            .collect::<Vec<_>>()
            .first()
            .ok_or("first not found")?
            .to_string();
        let alt_base = variant
            .clone()
            .split_terminator(':')
            .collect::<Vec<_>>()
            .last()
            .ok_or("last not found")?
            .split_terminator('>')
            .collect::<Vec<_>>()
            .last()
            .ok_or("last not found")?
            .to_string();

        let alleles: &[&[u8]] = &[ref_base.as_bytes(), alt_base.as_bytes()];
        record.set_alleles(alleles)?;
        record.set_filters(&["PASS".as_bytes()])?;

        // heterozygous
        let mut gts = Vec::new();
        for i in 0..allele_len {
            let allele_index = (i + 1) as i64;
            if alleles_per_variant.contains(&allele_index) {
                gts.push(GenotypeAllele::Unphased(1));
            } else {
                gts.push(GenotypeAllele::Unphased(0));
            }
        }
        record.push_genotypes(&gts)?;
        writer.write(&record)?;
    }

    Ok(())
}

/// Get RU INFO field
/// # Arguments
/// * `results` - variant information
/// # Returns
/// * `(String, HashSet<i64>)` - encoding and the value
fn encode_ru_field(results: Vec<VariantInfoByVariant>) -> Result<(String, HashSet<i64>), DError> {
    let mut alleles = HashSet::new();
    let mut encoding = String::new();
    for hap in results {
        if !encoding.is_empty() {
            encoding += ",";
        }
        encoding += &hap.fingerprint.to_string();
        encoding += ":";
        let allele_name = &hap.fingerprint_in_allele;
        if let Some(allele_name_string) = allele_name {
            encoding += allele_name_string;
            let allele_index = allele_name_string
                .to_string()
                .split_terminator('.')
                .collect::<Vec<_>>()
                .first()
                .ok_or("first not found")?
                .parse::<i64>()
                .unwrap();
            alleles.insert(allele_index);
        } else {
            encoding += "Unknown";
        }

        encoding += ":";
        encoding += &hap.depth.to_string();
        encoding += ":";
        encoding += &hap.nread.to_string();
    }
    Ok((encoding, alleles))
}
