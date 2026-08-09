use crate::realignment::utilities::Variant;
use log::debug;
use rust_htslib::bam;
use rust_htslib::bam::header::HeaderRecord;
use std::collections::{BTreeMap, HashSet};
use std::fmt::Display;
use std::path::{Path, PathBuf};
use std::str;
use std::str::FromStr;

pub type DError = std::boxed::Box<dyn std::error::Error>;
pub type DResult = Result<(), DError>;
pub type Exception = simple_error::SimpleError;

/// Build a standardized invalid-data error for runtime state or input that
/// violates expectations but should not panic the process.
pub fn invalid_data_error(message: impl Into<String>) -> DError {
    std::io::Error::other(format!("invalid data: {}", message.into())).into()
}

/// Build a standardized missing-data error for required values that could not
/// be found in the current execution context.
pub fn missing_data_error(subject: &str, context: impl Into<String>) -> DError {
    std::io::Error::other(format!("missing {subject}: {}", context.into())).into()
}

lazy_static::lazy_static! {
    pub static ref GIT_DESCRIBE: String = option_env!("VERGEN_GIT_DESCRIBE")
        .unwrap_or("unknown")
        .to_string();
    pub static ref FULL_VERSION: String = env!("CARGO_PKG_VERSION").to_string();
    pub static ref FULL_VERSION_PROGRAM: String =
        format!("{}-{}", env!("CARGO_PKG_NAME"), env!("CARGO_PKG_VERSION"));
    pub static ref CLI_COMMAND: String = std::env::args().collect::<Vec<_>>().join(" ");
}

/// Append this Kivvi invocation metadata as a `@PG` record.
pub fn append_kivvi_pg_header(header: &mut bam::Header) {
    let mut pg = HeaderRecord::new(b"PG");
    pg.push_tag(b"ID", env!("CARGO_PKG_NAME"));
    pg.push_tag(b"PN", env!("CARGO_PKG_NAME"));
    pg.push_tag(b"VN", &FULL_VERSION[..]);
    pg.push_tag(b"CL", &CLI_COMMAND[..]);
    header.push_record(&pg);
}

#[must_use]
/// Build a BAM header with `@PG` metadata for this Kivvi invocation.
pub fn output_bam_header(template: &bam::HeaderView) -> bam::Header {
    let mut header = bam::Header::from_template(template);
    append_kivvi_pg_header(&mut header);
    header
}

pub fn kivvi_temp_root(output_dir: &Path) -> PathBuf {
    output_dir.join(".kivvi-tmp")
}

pub fn create_kivvi_temp_dir(output_dir: &Path) -> Result<tempfile::TempDir, DError> {
    let temp_root = kivvi_temp_root(output_dir);
    std::fs::create_dir_all(&temp_root)?;
    Ok(tempfile::Builder::new()
        .prefix("kivvi-")
        .tempdir_in(temp_root)?)
}

#[must_use]
pub fn normalize_chrom_name(name: &str) -> &str {
    name.strip_prefix("chr").unwrap_or(name)
}

#[must_use]
pub fn resolve_chrom_name_from_header(header: &bam::HeaderView, requested: &str) -> Option<String> {
    let targets = header
        .target_names()
        .into_iter()
        .map(|x| String::from_utf8_lossy(x).to_string())
        .collect::<Vec<_>>();

    if targets.iter().any(|x| x == requested) {
        return Some(requested.to_string());
    }

    let requested_normalized = normalize_chrom_name(requested);
    targets
        .into_iter()
        .find(|x| normalize_chrom_name(x) == requested_normalized)
}

/// Decode embedded UTF-8 bytes into a string, falling back to lossy decoding
/// when bundled data contains invalid UTF-8.
fn embedded_utf8_string(bytes: &[u8], label: &str) -> String {
    match str::from_utf8(bytes) {
        Ok(value) => value.to_string(),
        Err(error) => {
            debug!("Failed to decode embedded {label} as UTF-8: {error}");
            String::from_utf8_lossy(bytes).to_string()
        }
    }
}

/// Parse newline-delimited embedded values, skipping blank or malformed lines
/// instead of aborting region initialization.
fn parse_embedded_numeric_lines<T>(bytes: &[u8], label: &str) -> Vec<T>
where
    T: FromStr,
    T::Err: Display,
{
    embedded_utf8_string(bytes, label)
        .lines()
        .filter_map(|line| {
            let trimmed = line.trim();
            if trimmed.is_empty() {
                return None;
            }
            match trimmed.parse::<T>() {
                Ok(value) => Some(value),
                Err(error) => {
                    debug!("Skipping invalid {label} entry {trimmed:?}: {error}");
                    None
                }
            }
        })
        .collect()
}

/// Parse the last tab-delimited column from embedded BED-like data as numeric
/// positions, skipping malformed lines.
fn parse_embedded_last_column_i64(bytes: &[u8], label: &str) -> Vec<i64> {
    embedded_utf8_string(bytes, label)
        .lines()
        .filter_map(|line| {
            let trimmed = line.trim();
            if trimmed.is_empty() {
                return None;
            }
            let Some(last_column) = trimmed.split('\t').next_back() else {
                debug!("Skipping malformed {label} line without tab-delimited columns: {trimmed:?}");
                return None;
            };
            match last_column.parse::<i64>() {
                Ok(value) => Some(value),
                Err(error) => {
                    debug!(
                        "Skipping invalid {label} last-column value {last_column:?} from line {trimmed:?}: {error}"
                    );
                    None
                }
            }
        })
        .collect()
}

/// Append a statically configured variant when it can be constructed, logging
/// invalid embedded definitions instead of panicking.
fn push_configured_variant<E>(target: &mut Vec<Variant>, variant: Result<Variant, E>, label: &str)
where
    E: Display,
{
    match variant {
        Ok(variant) => target.push(variant),
        Err(error) => debug!("Skipping invalid configured variant {label}: {error}"),
    }
}

/// For storing starting and ending reads
#[derive(Clone, Debug)]
pub struct FlankReads {
    // full read overlapping upstream of repeat
    pub start: HashSet<String>,
    // full read overlapping downstream of repeat
    pub end: HashSet<String>,
    // read segment for the first copy of repeat
    pub start_segment: HashSet<String>,
    // read segment for the last copy of repeat
    pub end_segment: HashSet<String>,
    // positions with 5p clips
    pub good_clips_p5: Vec<i64>,
    // positions with 3p clips
    pub good_clips_p3: Vec<i64>,
}

/// For storing important coordinates of each VNTR region
#[derive(Clone, Debug)]
pub struct RegionCoordinates {
    /// length of the repeat
    pub repeat_len: usize,
    /// chromosome name for output
    pub chromosome_output: String,
    /// length of the chromosome
    pub chromosome_len: i64,
    /// genome offset for each reference
    pub genome_offset: BTreeMap<String, i32>,
    /// reference sequence
    pub reference_seq: String,
    /// regions to extract reads containing the repeat
    pub extract_regions: Vec<String>,
    /// regions to extract flanking reads
    pub flanking_regions: Option<Vec<String>>,
    /// regions used for depth estimation
    pub depth_regions: (i64, i64),
    /// sites containing SNPs carried by a second subtype of kiv2
    pub type2_sites: Vec<i64>,
    /// sites to exclude for fingerprinting
    pub exclude_sites: Vec<i64>,
    /// sites to exclude for VCF output
    pub exclude_sites_vcf: Vec<i64>,
    /// regions to use to estimate genome depth
    pub genome_depth_sites: Vec<i64>,
    /// variants which are called by checking clip position
    pub clip_variant_sites: BTreeMap<i64, Vec<u8>>,
    /// sites to check methylation
    pub methyl_sites: Vec<usize>,
    /// clipping positions for starting reads
    pub start_positions_flank: Option<Vec<i64>>,
    /// clipping positions for ending reads
    pub end_positions_flank: Option<Vec<i64>>,
    /// site for determining allele type
    pub pivot_site: Option<i64>,
    /// variants to force genotype
    pub variants_to_call: Vec<Variant>,
    /// variants to exclude from genotype
    pub variants_to_exclude: Vec<Variant>,
    /// segments of the reference for global realignment
    pub realign_segments: Vec<i64>,
    /// variants to distinguish allele types
    pub variants_to_distinguish_allele_types: BTreeMap<String, Vec<String>>,
}

/// Get coordinates defined for KIV2
pub fn kiv2_coordinates() -> RegionCoordinates {
    lazy_static::lazy_static! {
        pub static ref HOMOPOLYMER: &'static [u8] = std::include_bytes!(concat!(env!("CARGO_MANIFEST_DIR"), "/data/kiv2/kiv2_homopolymer.txt"));
        pub static ref TYPE2SITES: &'static [u8] = std::include_bytes!(concat!(env!("CARGO_MANIFEST_DIR"), "/data/kiv2/kiv2_type2.txt"));
        pub static ref GENOMESITES: &'static [u8] = std::include_bytes!(concat!(env!("CARGO_MANIFEST_DIR"), "/data/genome_region.bed"));
        pub static ref REFERENCE: &'static [u8] = std::include_bytes!(concat!(env!("CARGO_MANIFEST_DIR"), "/data/kiv2/kiv2_ref.fa"));
    }
    let reference_seq = embedded_utf8_string(&REFERENCE, "kiv2 reference");
    let type2_sites = parse_embedded_numeric_lines::<i64>(&TYPE2SITES, "kiv2 type2 sites");
    debug!("type2 sites are: {type2_sites:?}");

    let genome_depth_sites = parse_embedded_last_column_i64(&GENOMESITES, "genome depth sites");
    debug!("genome_depth_sites sites are: {genome_depth_sites:?}");

    // here are some presets for kiv2
    let extract_regions = vec![String::from("chr6:160611535-160646860")];
    debug!("extract_regions are: {extract_regions:?}");

    let flanking_regions = vec![
        String::from("chr6:160605000-160610000"),
        String::from("chr6:160648000-160655000"),
    ];
    debug!("flanking_regions are: {flanking_regions:?}");

    let depth_regions = (1505, 3095);
    debug!("depth_regions are: {depth_regions:?}");

    // 1-based
    let mut exclude_sites: Vec<i64> = (1731..1768).map(|x| x as i64).collect::<Vec<_>>();
    exclude_sites.append(&mut vec![1534, 1535, 2883, 2884]);
    exclude_sites.sort();

    let mut variants_to_call = Vec::new();
    // use 0 based coordinate here
    push_configured_variant(
        &mut variants_to_call,
        Variant::new_deletion(
            0,
            193,
            9,
            "GTCCTTCTC".as_bytes().to_vec(),
            "G".as_bytes().to_vec(),
            0,
            1,
        ),
        "kiv2 deletion 193",
    );
    push_configured_variant(
        &mut variants_to_call,
        Variant::new_deletion(
            0,
            1174,
            4,
            "TGAC".as_bytes().to_vec(),
            "T".as_bytes().to_vec(),
            0,
            1,
        ),
        "kiv2 deletion 1174",
    );
    push_configured_variant(
        &mut variants_to_call,
        Variant::new_snv(
            0,
            1533,
            "T".as_bytes().to_vec(),
            "C".as_bytes().to_vec(),
            0,
            1,
        ),
        "kiv2 snv 1533",
    );
    push_configured_variant(
        &mut variants_to_call,
        Variant::new_snv(
            0,
            2882,
            "A".as_bytes().to_vec(),
            "G".as_bytes().to_vec(),
            0,
            1,
        ),
        "kiv2 snv 2882",
    );
    let mut genome_offset = BTreeMap::new();
    genome_offset.insert(String::from("chr6:160613619-160619170"), 160613617);

    RegionCoordinates {
        repeat_len: 5552,
        chromosome_output: String::from("chr6"),
        chromosome_len: 170805979,
        // To get IGV coordinates: 1-based, add 160613618; 0-based, add 160613619
        genome_offset,
        reference_seq,
        extract_regions,
        flanking_regions: Some(flanking_regions),
        depth_regions,
        type2_sites,
        exclude_sites,
        exclude_sites_vcf: vec![],
        genome_depth_sites,
        clip_variant_sites: BTreeMap::new(),
        methyl_sites: vec![],
        start_positions_flank: None,
        end_positions_flank: None,
        pivot_site: None,
        variants_to_call,
        variants_to_exclude: vec![],
        realign_segments: vec![],
        variants_to_distinguish_allele_types: BTreeMap::new(),
    }
}

/// Get coordinates defined for D4Z4
pub fn d4z4_coordinates() -> RegionCoordinates {
    lazy_static::lazy_static! {
        pub static ref HOMOPOLYMER: &'static [u8] = std::include_bytes!(concat!(env!("CARGO_MANIFEST_DIR"), "/data/d4z4/d4z4_homopolymer.txt"));
        pub static ref METHYLSITES: &'static [u8] = std::include_bytes!(concat!(env!("CARGO_MANIFEST_DIR"), "/data/d4z4/d4z4_methyl_sites.txt"));
        pub static ref GENOMESITES: &'static [u8] = std::include_bytes!(concat!(env!("CARGO_MANIFEST_DIR"), "/data/genome_region.bed"));
        pub static ref REFERENCE: &'static [u8] = std::include_bytes!(concat!(env!("CARGO_MANIFEST_DIR"), "/data/d4z4/d4z4_ref.fa"));
    }
    let reference_seq = embedded_utf8_string(&REFERENCE, "d4z4 reference");
    let methyl_sites = parse_embedded_numeric_lines::<usize>(&METHYLSITES, "d4z4 methyl sites");
    //debug!("methyl_sites sites are: {methyl_sites:?}");

    let genome_depth_sites = parse_embedded_last_column_i64(&GENOMESITES, "genome depth sites");
    //debug!("genome_depth_sites sites are: {genome_depth_sites:?}");

    // here are some presets for kiv2
    let extract_regions = vec![
        String::from("chr4:190065229-190093263"),
        String::from("chr4:190173122-190175903"),
        String::from("chr10:133664430-133685491"),
        String::from("chr10:133740609-133761980"),
    ];
    //debug!("extract_regions are: {extract_regions:?}");

    // 1-based
    let mut exclude_sites: Vec<i64> = (1136..1209)
        //.chain(1093..1098)
        .chain(127..173)
        .chain(3111..3116)
        .map(|x| x as i64)
        .collect::<Vec<_>>();
    exclude_sites.append(&mut vec![3001]);
    exclude_sites.sort();
    let exclude_sites_vcf: Vec<i64> = (1136..1209).map(|x| x as i64).collect::<Vec<_>>();

    let mut variants_to_exclude = Vec::new();
    // use 0 based coordinate here
    push_configured_variant(
        &mut variants_to_exclude,
        Variant::new_snv(
            0,
            683,
            "A".as_bytes().to_vec(),
            "G".as_bytes().to_vec(),
            0,
            1,
        ),
        "d4z4 exclude snv 683 A>G",
    );
    push_configured_variant(
        &mut variants_to_exclude,
        Variant::new_snv(
            0,
            2115,
            "C".as_bytes().to_vec(),
            "A".as_bytes().to_vec(),
            0,
            1,
        ),
        "d4z4 exclude snv 2115 C>A",
    );
    push_configured_variant(
        &mut variants_to_exclude,
        Variant::new_snv(
            0,
            2115,
            "C".as_bytes().to_vec(),
            "G".as_bytes().to_vec(),
            0,
            1,
        ),
        "d4z4 exclude snv 2115 C>G",
    );
    push_configured_variant(
        &mut variants_to_exclude,
        Variant::new_snv(
            0,
            3105,
            "A".as_bytes().to_vec(),
            "C".as_bytes().to_vec(),
            0,
            1,
        ),
        "d4z4 exclude snv 3105 A>C",
    );

    let start_positions_flank = vec![2058];
    let end_positions_flank = vec![429];
    let mut clip_variant_sites = BTreeMap::new();
    clip_variant_sites.entry(3296).or_insert(vec![b'A']);

    let mut variants_to_call = Vec::new();
    // use 0 based coordinate here
    push_configured_variant(
        &mut variants_to_call,
        Variant::new_deletion(
            0,
            215,
            5,
            "GGGAC".as_bytes().to_vec(),
            "G".as_bytes().to_vec(),
            0,
            1,
        ),
        "d4z4 deletion 215",
    );
    push_configured_variant(
        &mut variants_to_call,
        Variant::new_insertion(
            0,
            259,
            "G".as_bytes().to_vec(),
            "GCCTCCGGGAGTAGCGGGACCCCCGC".as_bytes().to_vec(),
            0,
            1,
        ),
        "d4z4 insertion 259",
    );
    push_configured_variant(
        &mut variants_to_call,
        Variant::new_deletion(
            0,
            306,
            42,
            "TCCGGCGCGGGCTGAGGGCTGGGCCCACAGCCGCCGCGCCGG"
                .as_bytes()
                .to_vec(),
            "T".as_bytes().to_vec(),
            0,
            1,
        ),
        "d4z4 deletion 306",
    );
    push_configured_variant(
        &mut variants_to_call,
        Variant::new_deletion(
            0,
            315,
            45,
            "GGCTGAGGGCTGGGCCCACAGCCGCCGCGCCGGCCGGCGCGGCAC"
                .as_bytes()
                .to_vec(),
            "G".as_bytes().to_vec(),
            0,
            1,
        ),
        "d4z4 deletion 315",
    );
    push_configured_variant(
        &mut variants_to_call,
        Variant::new_deletion(
            0,
            356,
            4,
            "GCAC".as_bytes().to_vec(),
            "G".as_bytes().to_vec(),
            0,
            1,
        ),
        "d4z4 deletion 356",
    );
    push_configured_variant(
        &mut variants_to_call,
        Variant::new_insertion(
            0,
            1042,
            "C".as_bytes().to_vec(),
            "CCGCG".as_bytes().to_vec(),
            0,
            1,
        ),
        "d4z4 insertion 1042",
    );
    push_configured_variant(
        &mut variants_to_call,
        Variant::new_deletion(
            0,
            1087,
            8,
            "GCCCCCCT".as_bytes().to_vec(),
            "G".as_bytes().to_vec(),
            0,
            1,
        ),
        "d4z4 deletion 1087",
    );
    push_configured_variant(
        &mut variants_to_call,
        Variant::new_deletion(
            0,
            1094,
            10,
            "TCCCCCCTCC".as_bytes().to_vec(),
            "T".as_bytes().to_vec(),
            0,
            1,
        ),
        "d4z4 deletion 1094",
    );
    push_configured_variant(
        &mut variants_to_call,
        Variant::new_deletion(
            0,
            1216,
            18,
            "CCCAGGCCTCGACGCCCT".as_bytes().to_vec(),
            "C".as_bytes().to_vec(),
            0,
            1,
        ),
        "d4z4 deletion 1216",
    );
    push_configured_variant(
        &mut variants_to_call,
        Variant::new_insertion(
            0,
            1327,
            "C".as_bytes().to_vec(),
            "CGGCT".as_bytes().to_vec(),
            0,
            1,
        ),
        "d4z4 insertion 1327",
    );
    push_configured_variant(
        &mut variants_to_call,
        Variant::new_deletion(
            0,
            1753,
            29,
            "CGAAAGCGGACCGCCGTCACCGGATCCCA".as_bytes().to_vec(),
            "C".as_bytes().to_vec(),
            0,
            1,
        ),
        "d4z4 deletion 1753",
    );
    push_configured_variant(
        &mut variants_to_call,
        Variant::new_deletion(
            0,
            2447,
            326,
            "CCGGGGCAGCTCCACCTCCCCAGCCCGCGCCCCCGGACGCCTCCGCCTCCGCGCGGCAGGGGCAGATGCAAGGCATCCCGGCGCCCTCCCAGGCGCTCCAGGAGCCGGCGCCCTGGTCTGCACTCCCCTGCGGCCTGCTGCTGGATGAGCTCCTGGCGAGCCCGGAGTTTCTGCAGCAGGCGCAACCTCTCCTAGAAACGGAGGCCCCGGGGGAGCTGGAGGCCTCGGAAGAGGCCGCCTCGCTGGAAGCACCCCTCAGCGAGGAAGAATACCGGGCTCTGCTGGAGGAGCTTTAGGACGCGGGGTTGGGACGGGGTCGGGTGGTT".as_bytes().to_vec(),
            "C".as_bytes().to_vec(),
            0,
            1,
        ),
        "d4z4 deletion 2447",
    );
    push_configured_variant(
        &mut variants_to_call,
        Variant::new_deletion(
            0,
            2483,
            7,
            "ACGCCTC".as_bytes().to_vec(),
            "A".as_bytes().to_vec(),
            0,
            1,
        ),
        "d4z4 deletion 2483",
    );
    push_configured_variant(
        &mut variants_to_call,
        Variant::new_deletion(
            0,
            2821,
            325,
            "CGGAGGGGCGTGTCTCCGCCCCGCCCCCTCCACCGGGCTGACCGGCCTGGGATTCCTGCCTTCTAGGTCTAGGCCCGGTGAGAGACTCCACTCCGCGGAGAACTGCCTTTCTTTCCTGGGCATCCCGGGGATCCCAGAGCCGGCCCAGGTACCAGCAGGTGGGCCGCCTACTGCGCACGCGCGGGTTTGCGGGCAGCCGCCTGGGCTGTGGGAGCAGCCCGGGCAGAGCTCTCCTGCCTCTCCACCAGCCCACCCCGCCGCCTGACCGCCCCCTCCCCACCCCCACCCCCCACCCCCGGAAAACGCGTCGTCCCCTGGGCTGGGT".as_bytes().to_vec(),
            "C".as_bytes().to_vec(),
            0,
            1,
        ),
        "d4z4 deletion 2821",
    );
    push_configured_variant(
        &mut variants_to_call,
        Variant::new_deletion(
            0,
            3088,
            14,
            "GCCCCCTCCCCACC".as_bytes().to_vec(),
            "G".as_bytes().to_vec(),
            0,
            1,
        ),
        "d4z4 deletion 3088 len14",
    );
    push_configured_variant(
        &mut variants_to_call,
        Variant::new_deletion(
            0,
            3088,
            19,
            "GCCCCCTCCCCACCCCCAC".as_bytes().to_vec(),
            "G".as_bytes().to_vec(),
            0,
            1,
        ),
        "d4z4 deletion 3088 len19",
    );
    push_configured_variant(
        &mut variants_to_call,
        Variant::new_insertion(
            0,
            3289,
            "G".as_bytes().to_vec(),
            "GCCTGGCGGCGGAACGCAGACCCCAGGCCCGGCGCACACCGGGGACGCTGAGCGTTCCAGGCGGGAGGGAAGGCGGGCAGAGATGGAGAGAGGAACGGGAGACCTAGAGGGGCGGAAGGACGGGCGGAGGGACGTTAGGAGGGAGGGAGGGAGGCAGGGAGGCAGGGAGGAACGGAGGAAAGACAGAGCGACGCAGGGACTGGGGGCGGGCGGGAGGGAGCCGGGGACGGACGGGGGGAGGAAGGCAGGGAGGAAAAGCGGTCCTCGGCCTCCGGGAGTAGCGGGACCCCCGCCCTCCGGGAAAACGGTCAGCGTCCGGCGCGGGCTGAGGGCTGGGCCCACAGCCGCCGCGCCGGCCGGCGGGGCACCACCCATTCGCCCCGGTTCCGGGGCCCAGGGAGTGGGCGGTTTCCTCCGGGACAAAAGACCGGGACTCGGGTTGCCGTCGGGTCTTCACCCGCGCGGTTCACAGACCGCACATCCCCAGGCTGAGCCCTGCAACGCGGCGCGAGGCCGACAGCCCCGGCCACGGAGGAGCCACACGCAGGACGACGGAGGCGTGATTTTGGTTTCCGCGTGGCTTTGCCCTCCGCAAGGCGGCCTGTTGCTCACGTCTCTCCGGCCCCCGAAAGGCTGGCCATGCCGACTGTTTGCTCCCGGAGCTCTGCGGGCACCCGGAAACATGCAGGGAAGGGTGCAAGCCGGCACGGTGCCTTCGCTCTCCTTGCCAGGTTCCAAACCGGCCACACTGCAGACTCCCCACGTTGCCGCACGCGGGAATCCATCGTCAGGCCATCACGCCGGGGAGGCATCTCCTCTCTGGGGTCTCGCTCTGGTCTTCTACGTGGAAATGAACGAGAGCCACACGCCTGCGTGTGCGAGACCGTCCCGGCAACGGCGACGCCCACAGGCATTGCCTCCTTCACGGAGAGAGGGCCTGGCACACTCAAGACTCCCACGGAGGTTCAGTTCCACACTCCCCTCCACCCTCCCAGGCTGGTTTCTCCCTGCTGCCGACGCGTGGGAGCCCAGAGAGCGGCTTCCCGTTCCCGCGGGATCCCTGGAGAGGTCCGGAGAGCCGGCCCCCGAAACGCGCCCCCCTCCCCCCTCCCCCCCTCCTCCCCGTTCCTCTTCGTCTGTCCGGCCCCACCACCACCACCGCCACCACGCCCTCCCCCCCCACCCCCCCCCCCCACCACCACCACCACCACCACCCCGCCGGCCGGCCCCAGGCCTCGACGCCCTGGGTCCCTTCCGGGGTGGGGCGGGCTGTCCCAGGGGGGCTCACCGCCATTCATGAAGGGGTGGAGCCTGCCTGCCTGTGGGCCTTTACAAGGGCGGCTGGCTGGCTGGCTGGCTGTCCGGGCAGGCCCCCTGGCTGCACCTGCCGCAGTGCACAGTCCGGCTGAGGTGCACGGGAGCCCGCCGGCCTCTCTCTGCCCGCGTCCGTCCGTGAAATTCCGGCCGGGGCTCACCGCGATGGCCCTCCCGACACCCTCGGACAGCACCCTCCCCGCGGAAGCCCGGGGACGAGGACGGCGACGGAGACTCGTTTGGACCCCGAGCCAAAGCGAGGCCCTGCGAGCCTGCAGCCTCCCAGCTGCCAGCGCGGAGCT".as_bytes().to_vec(),
            0,
            1,
        ),
        "d4z4 insertion 3289",
    );
    push_configured_variant(
        &mut variants_to_call,
        Variant::new_insertion(
            0,
            3276,
            "G".as_bytes().to_vec(),
            "GCCAGCACGGAGCGCCTGGCGGCGGAACGCAGACCCCAGGCCCGGCGCACACCCGGGGGACGCTGAGCGTTCCAGGCGGGAGGGAAGGCGGGCAGAGATGGAGAGAGGAACGGGAGACCTAGAGGGGCGGAAGGATGGGCGGAGGGACGTTAGGAGGGAGGGAGGGAGGCAGGGAGGCAGGGAGGAACGGAGGGAAAGACAGAGCGACGCAGGGACTGGGGGCGGGCGGGAGGGAGCCGGGGACGGGGGGAGGAAGGCAGGGAGGAAAAGCGGTCCTCGGCCTCCGGGAGTAGCGGACCCCCGCCCTCCGGGAAAACGGTCAGCGTCCGGCGCGGGCTGAGGGCTGGGCCCACAGCCGCCGCGCCGGCCGGCGGGGCACCACCCATTCGCCCCGGTTCCGGGGCCCAGGGAGTGGGCGGTTTCCTCCGGGACAAAAGACCGGGACTCGGGTTGCCGTCGGGTTTTCACCCGCGCGGTTCACAGACCGCACATCCCCAGGCTGAGCCCTGCAACGGGGCGCGAGGCCGACAGCCCCGGCCACGGAGGAGCCACACGCAGGACGACGGAGGCGTGATTTTGGTTTCCGCGTGGCTTTGCCCTCTGCAAGGCGGCCTGTTGCTCACGTCTCTCCGGCCCCCGAAAGCTGGCCATGCCGACTGTTTGCTCCCGGAGCTCTGCGGGCACCCGGAAACATGCAGGGAAGGGTGCAAGGCCCGGCACGGTGCCTTCGCTCTCCTTGCCAGGTTCCAAACCGGCCACACTGCAGACTCCCCACGTTGCCGCACGCGGGAATCCATCGTCAGGCCATCACGCCGGGGAGGCATCTCCTCTCTGGGGTCTCGCTCTGGTCTTCTACGTGGAAATGAACGAGAGCCACACGCCTGCGTGTGCGAGACCGTCCCGGCAACGGCGACGCCCACAGGCATTGCCTCCTTCACGGAGAGAGGGCCTGGCACACTCAAGACTCCCACGGAGGTTCAGTTCCACACTCCCCTCCACCCTCCCAGGCTGGTTTCTCCCTGCTGCCGACGCGTGGGAGCCCAGAGAGCGGCTTCCCGTTCCCGCGGGATCCCTGGAGAGGTCCGGAGAGCCGGCCCCCGAAACGCGCCCCCCTCCCCCCTCCCCCCTCTCCCCCTTCCTCTTCGTCTCTCCGGCCCCACCACCACCACCGCCACCACGCCCTCCCCCACCACCCCCCCCCCCCACCACCACCACCACCCCGCCGGCCGGCCCCAGGCCTCGACGCCCTGGGTCCCTTCCGGGGTGGGGCGGGCTGTCCCAGGGGGGGCTCACCGCCATTCATGAAGGGGTGGAGCCTGCCTGCCTGTGGGCCTTTACAAGGGCGGCTGGCTGGCTGGGCTGGCTGTCCGGGCAGGCCTCCCTGGCTGCACCTGCCGCAGCGCACAGTCCGGCTGAGGTGCACGGGAGCCCGCCGGCCTCTCTCTGCCCGCGTCCGTCCGTGAAATTCCGGCCGGGGCTCACCGCGATGGCCCTCCCGACACCTTCGGACAGCACCCTCCCCGCGGAAGCCCGGGGACGAGGACGGCGACGGAGACTCGTTTGGACCCCGAGCCAAAGCGAGGCCCTGCGAGCCTGCTTTGAGCGGAACCCGTACCCGGGCATCGCCACCAGAGAACGGCTGGCCCAGGCCATCGGCATTCCGGAGCCCAGGGTCCAGATTTGGTTTCAGAATGAGAGGTCACGCCAGCTGAGGCAGCACCGGCGGGAATCTCGGCCCTGGCCCGGGAGACGCGGCCCGCCAGAAGGCCGGCGAAAGCGGACCGCCGTCACCGGATCCCAGACCGCCCTGCTCCCTCCGAGCCTTTGAGAAGGATCGCTTTCCAGGCATCGCCGCCCGGGAGGAGCTGGCCAGAGAGACGGGCCTCCCGGAGTCCAGGATTCAGATCTGGTTTCAGAATCGAAGGGCCAGGCACCCGGGACAGGGTGGCAGGGCGCCCGCGCAGGCAGGCGGCCTGTGCAACGCGGCCCCCGGCGGGGGTCACCCTGCTCCCTCCGTGGGTCGCCTTCGCCCACACCGGCGCGTGGGGAACGGGGCTTCCCGCACCCCACGTGCCCTGCGCGCCTGGGGCTCTCCCACAGGGGGCTTTCGTGAGCCAGGCAGCGAGGGCCGCCCCCGCGCTGCAGCCCAGCCAGGCCGCGCCGGCAGAGGGGATCTCCCAACCTGCCCCGGCGCGCGGGGATTTCGCCTACGCCGCCCCGGCTCCTCCGGACGGGGCGCTCTCCCACCCTCAGGCTCCTCGGTGGCCTCCGCACCCGGGCAAAAGCCGGGAGGACCGGGACCCAGCAGCGCGACGGCCTGCCGGGCCCTGCGCGGTGGCACAGCCTGGGCCCGCTCAAGCGGGGCCGCAGGGCCAAGGGGTGCTTGCGCCACCCACGTCCCAGGGGAGTCCGTGGTGGGGCTGGGGCCGGGGTCCCCAGGTCGCCGGGGCGGCGTGGAACCCCAAGCCGGGGCAGCTCCACCTCCCCAGCCC".as_bytes().to_vec(),
            0,
            1,
        ),
        "d4z4 insertion 3276",
    );
    push_configured_variant(
        &mut variants_to_call,
        Variant::new_insertion(
            0,
            3295,
            "C".as_bytes().to_vec(),
            "CGCGGAACGCAGACCCCAGGCCCGGCGCACACCGGGGACGCTGAGCGTTCCAGGCGGGAGGGAAGGCGGGCAGAGATGGAGAGAGGAACGGGAGACCTAGAGGGGCGGAAGGATGGGCGGAGGGACGTTAGGAGGGAGGGAGGCAGGGAGGCAGGGAGGCAGGGAGGAACGGAGGGAAAGACAGAGCGACGCAGGGACTGGGGGCGGGCGGGAGGGAGCCGGGGACGGACGGGGGGAGGAAGGCAGGGAGGAAAAGCGGTCTTCGGCCTCCGGGAGTAGCGGGACCCCCGCCCTCCGGGAAAACGGTCAGCGTCCGGCGCGGGCTGAGGGCTGGGCCCACAGCCGCCGCGCCGGCCGGCGGGGCACCACCCATTCGCCCCGGTTCCGGGGCCCAGGGAGTGGGCGGTTTCCTCCGGGACAAAAGACCGGGACTCGGGTTGCCGTCGGGTCTTCACCCGCGCGGTTCACAGACCGCACATCCCCAGGCTGAGCCCTGCAACGCGGCGCGAGGCCGACAGCCCCGGCCACGGAGGAGCCACACGCAGGACGACGGAGGCGTGATTTTGGTTTCCGCGTGGCTTTGCCCTCCGCAAGGCGGCCTGTTGCTCACGTCTCTCCGGCCCCCGAAAGGCCGGCCATGCCGACTGTTTGCTCCCGGAGCTCTGCCGGCACCCGGAAACATGCAGGGAAGGGTGCAAGCCCGGCACGGTGCCTTCGCTCTCCTTGCCAGGTTCCAAACGGCCACACTGCAGACTCCCCACGTTGCCGCACGCGGGAATCCATCGTCAGGCCATCACGCCGGGGAGGCATCTCCTCTCTGGGGTCTCGCTCTGGTCTTCTACGTGGAAATGAACGAGAGCCACACGCCTGCGTGTGCGAGACCGTCCCGGCAACGGCGACGCCCACAGGCATTGCCTCCTTCACGGAGAGAGGGCCTGGCACACTCGAGACTCCCACGGAGGTTCAGTTCCACACTCCCCTCCACCCTCCCAGGCTGGTTTCTCCCTGCTGCCGACGCGTGGGAGCCCAGAGAGCGGCTTCCCGTTCCCGCGGGATCCCTGGAGAGGTCCGGAGAGCCGGCCCCCGAAACGCGCCCCCCCTCCCCCCTCCCCCCTCTCCCCCTTCCTCTTCGTCTCTCCGGCCCCACCACCACCACCGCCACCACGCCCTCCCCCACCACCCCCCCCCCACCACCACCACCACCACCACCACCCCGCCGGCCGGCCCCAGGCCTCGACGCCCTGGGTCCCTTCCGGGGTGGGGCGGGCTGTCCCAGGGGGGCTCACCGCCATTCATGAAGGGGTGGAGCCTGCCTGCCTGTGGGCCTTTACAAGGGCGGCTGGCTGGCTGGCTGGCTGGCTGTCCGGGCAGGCCTCCTGGCTGCACCTGCCGCAGTGCACAGTCCGGCTGAGGTGCACGGGAGCCCGCCGGCCTCTCTCTGCCCGCGTCCGTCCGTGAAATTGCGGCCGGGGCTCACCGCGATGGCCCTCCCGACACCTTCGGACAGCACCCTCCCCGCGGAAGCCCGGGGACGAGGACGGCGACGGAGACTCGTTTGGACCCCGAGCCAAAGCGAGGCCCTGCGAGCCTGCTTTGAGCGGAACCCGTACCCGGGCATCGCCACCAGAGAACGGCTGGCCCAGGCCATCGGCATTCCGGAGCCCAGGGTCCAGATTTGGTTTCAGAATGAGAGGTCACGCCAGCTGAGGCAGCACCGGCGGGAATCTCGGCCCTGGCCCGGGAGACGCGGCCCGCCAGAAGGCCCGGCGAAAGCGGACCGCCGTCACCGGATCCCAGACCGCCCTGCTCCTCCGAGCCTTTGAGAAGGATCGCTTTCCAGGCATCGCCGCCCGGGAGGAGCTGGCCAGAGAGACGGGCCTCCCGGAGTCCAGGATTCAGATCTGGTTTCAGAATCGAAGGGCCAGGCACCCGGGACAGGGTGGCAGGGCGCCCGCGCAGGCAGGCGGCCTGTGCAGCGCGGCCCCCGGCGGGGGTCACCCTGCTCCCTCGTGGGTCGCCTTCGCCCACACCGGCGCGTGGGGAACGGGGCTTCCCGCACCCCACGTGCCCTGCGCGCCTGGGGCTCTCCCACAGGGGGCTTTCGTGAGCCAGGCAGCGAGGGCCGCCCCCGCGCTGCAGCCCAGCCAGGCCGCGCCGGCAGAGGGGGTCTCCCAACCTGCCCCGGCGCGCGGGGGATTTCGCCTACGCCGCCCCGGCTCCTCCGGACGGGGCGCTCTCCCACCCTCAGGCTCCTCGGTGGCCTCCGCACCCGGCAAAAGCCGGGAGGACCGGGACCCGCAGCGCGACGGCCTGCCGGGCCCCTGCGCGGTGGCACAGCCTGGGCCCGCTCAAGCGGGGCCGCAGGGCCAAGGGGTGCTTGCGCCACCCACGTCCCAGGGGAGTCCGTGGTGGGGCTGGGGCCGGGGTCCCCCAGGTCGCCGGGGCGGCGTGGGAACCCCAAGCCGGGGCAGCTCCACCTCCCCAGCCCGCGCCCCCAGGACGCCTCCGCCTCCGCGCGGCAGGGGCAGATGCAAGGCATCCCGGCGCCCTCCCAGGCGCTCCGGGAGCCGGCGCCCTGGTCTGCACTCCCCTGCGGCCTGCTGCTGGATGAGCTCCTGGCGAGCCCGGAGTTTCTGCAGCAGGCGCAACCTCTCCTAGAAACGGAGGCCCCGGGGGAGCTGGAGGCCTCGGAAGAGGCCGCCTCGCTGGAAGCACCCCTCAGCGAGGAAGAATACCGGGCTCTGCTGGAGGAGCTTTAGGACGCGGGGTTGGGACGGGGTCGGGTGGCTCGGGGCAGGGCGGTGGCCTCTCTTTCGCGGGGAACACCTGGCTGGCTACGGAGACCCCCGTCCCGCGAAACACCGGGCCCCGCGCAGCGTCCGGGCCTGACACCGCTCCGGCGGCTCGCCTCCTCTGCGCCCCCGCGCCACCGTCGCCCGCCCGCCCGGGCCCCTGCAGCCTCCCAGCTGCCAGCAGGGAGCGCCTGGCT".as_bytes().to_vec(),
            0,
            1,
        ),
        "d4z4 insertion 3295",
    );

    let mut genome_offset = BTreeMap::new();
    genome_offset.insert(String::from("d4z4_ref"), -1);

    let mut variants_to_distinguish_allele_types = BTreeMap::new();
    variants_to_distinguish_allele_types.insert(
        String::from("qADisruptedPolyA"),
        vec![
            String::from("256:C>T"),
            String::from("583:C>T"),
            String::from("683:A>G"),
        ],
    );
    variants_to_distinguish_allele_types.insert(
        String::from("qB"),
        vec![
            String::from("172:A>G"),
            String::from("513:C>A"),
            String::from("2891:T>C"),
        ],
    );

    RegionCoordinates {
        repeat_len: 3298,
        //chromosome_len: 190214555,
        //genome_offset: 190172603,
        chromosome_output: String::from("d4z4_ref"),
        chromosome_len: 4203,
        genome_offset,
        reference_seq,
        extract_regions,
        flanking_regions: None,
        depth_regions: (1000, 3000),
        type2_sites: vec![],
        exclude_sites,
        exclude_sites_vcf,
        genome_depth_sites,
        clip_variant_sites,
        methyl_sites,
        start_positions_flank: Some(start_positions_flank),
        end_positions_flank: Some(end_positions_flank),
        pivot_site: Some(3509),
        variants_to_call,
        variants_to_exclude,
        realign_segments: vec![127, 162, 1136, 1206],
        variants_to_distinguish_allele_types,
    }
}
