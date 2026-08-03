use crate::util::{
    invalid_data_error, missing_data_error, DError, DResult, RegionCoordinates, FULL_VERSION,
};
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
        .map(|x| parse_variant_key(x).map(|parsed| parsed.position))
        .collect::<Result<Vec<i64>, DError>>()?;
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
    let version_line = format!("##{}Version={}", env!("CARGO_PKG_NAME"), &*FULL_VERSION);
    vcf_header.push_record(version_line.as_bytes());
    //vcf_header.push_sample(sample_name.as_bytes());
    for i in 0..allele_len {
        let allele_name = format!("allele{}", i + 1);
        vcf_header.push_sample(allele_name.as_bytes());
    }

    let mut writer = bcf::Writer::from_path(output_path, &vcf_header, true, Format::Vcf)
        .map_err(|_| format!("Invalid VCF output path: {}", output_path.display()))?;

    for variant in all_var_sorted {
        let variant_info = variant_summary
            .get(&variant)
            .ok_or_else(|| missing_data_error("variant summary entry", &variant))?;
        let mut record = writer.empty_record();
        let parsed_variant = parse_variant_key(&variant)?;

        let contig = region_coordinates.chromosome_output.as_bytes();
        let rid = writer.header().name2rid(contig)?;
        record.set_rid(Some(rid));

        let variant_pos = parsed_variant.position;
        record.set_pos(variant_pos - 1);

        // variant quality?
        record.set_qual(f32::from_bits(0x7F800001));
        //record.push_info_integer(b"RPOS", &[variant_pos as i32])?;

        let (data, alleles_per_variant) = encode_ru_field(variant_info.to_vec())?;
        record.push_info_string(b"RU", &[data.as_bytes()])?;

        let alleles: &[&[u8]] = &[
            parsed_variant.ref_base.as_bytes(),
            parsed_variant.alt_base.as_bytes(),
        ];
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
                .split_terminator('.')
                .next()
                .ok_or_else(|| {
                    invalid_data_error(format!(
                        "Unexpected allele name format in RU field: '{allele_name_string}'"
                    ))
                })?
                .parse::<i64>()
                .map_err(|e| {
                    format!(
                        "Failed to parse allele index from RU field '{allele_name_string}': {e}"
                    )
                })?;
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

/// Parsed components of the internal variant key format used while building VCF
/// records.
struct ParsedVariantKey<'a> {
    position: i64,
    ref_base: &'a str,
    alt_base: &'a str,
}

/// Parse a variant key of the form `POSITION:...:REF>ALT`.
/// # Arguments
/// * `variant` - internal variant identifier used by Kivvi
/// # Returns
/// * `ParsedVariantKey` - parsed position and allele components for VCF output
fn parse_variant_key(variant: &str) -> Result<ParsedVariantKey<'_>, DError> {
    let mut fields = variant.split_terminator(':');
    let position = fields
        .next()
        .ok_or_else(|| missing_data_error("variant key position", variant))?
        .parse::<i64>()
        .map_err(|e| format!("Failed to parse variant position from '{variant}': {e}"))?;
    let allele_field = fields
        .next_back()
        .ok_or_else(|| missing_data_error("variant key allele field", variant))?;
    let (ref_base, alt_base) = allele_field.split_once('>').ok_or_else(|| {
        invalid_data_error(format!(
            "Variant allele field is not REF>ALT in '{variant}'"
        ))
    })?;
    Ok(ParsedVariantKey {
        position,
        ref_base,
        alt_base,
    })
}
