use crate::assembly::assembly_result::AssembledPaths;
use crate::config::{self, Region as RegionConfig};
use crate::detail::phase_haps::HapInfoForJson;
//use crate::detail::site_selection::{self, FilteredSitesForJson};
use crate::detail::util::DError;

use vstr::VString;

use itertools::Itertools;

use std::collections::BTreeMap;
use std::fmt;

pub type Error = simple_error::SimpleError;

#[allow(dead_code)]
fn parse_json(path: &std::path::Path) -> Result<serde_json::Value, DError> {
    Ok(serde_json::from_reader(std::fs::File::open(path)?)?)
}

///
/// Primary return type of Phaser.
/// Has fixed attributes for results from the pipeline.
///
/// Also has additional `BTreeMap<String, serde_json::Value>` member for holding flexible key-value metadata storage.
#[must_use]
#[derive(Clone, Debug, Default, serde::Serialize)]
pub struct GeneCall {
    pub gene_name: String,
    pub phase_region: String,
    pub sample_sex: String,
    pub genome_depth: Option<f32>,
    pub region_depth: BTreeMap<String, f32>,
    pub failed_for_coverage: bool,

    pub total_cn: Option<i32>,
    pub final_haplotypes: BTreeMap<String, String>,
    pub two_copy_haplotypes: Vec<String>,
    pub region_specific_info: BTreeMap<String, serde_json::Value>, // additional key-value metadata

    pub sites_for_phasing: Vec<String>,
    pub assembled_haplotypes: Vec<String>,
    pub unique_supporting_reads: BTreeMap<String, Vec<String>>,

    pub highest_total_cn: Option<i32>,
    pub heterozygous_sites: Vec<String>,
    pub het_sites_not_used_in_phasing: Vec<String>,
    pub homozygous_sites: Vec<String>,
    pub haplotype_details: BTreeMap<String, HapInfoForJson>,
    pub nonunique_supporting_reads: BTreeMap<String, Vec<String>>,
    pub read_details: BTreeMap<String, String>,
}

/// Represents a read and an optional alignment id.
/// `align_id` lets us distinguish between multiple matches per read.
/// For instance, multiple repeat units in one read.
#[derive(
    Hash, PartialEq, Clone, PartialOrd, Ord, Eq, Default, serde::Deserialize, serde::Serialize,
)]
pub struct ReadAlignmentId {
    pub read_name: String,
}

impl fmt::Display for ReadAlignmentId {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        /*
        if let Some(id) = self.align_id {
            write!(f, "ReadAlignmentId{{name: {}, id: {id}}}", self.read_name)
        } else {
            write!(f, "{}", self.read_name)
        }
        */
        write!(f, "{}", self.read_name)
    }
}

/*
impl<'de> serde::Deserializer<'de> for ReadAlignmentId {
    fn deserialize<D>(deserializer: D) -> Result<i32, D::Error>
    where
        D: Deserializer<'de>,
    {
        fn deserialize<'de, D>(deserializer: D) -> Result<ReadAlignmentId, D::Error>
        where
            D: Deserializer<'de>,
        {
            struct ReadAlignmentIdVisitor;
            impl<'de> serde::de::Visitor<'de> for ReadAlignmentIdVisitor {
                type Value = ReadAlignmentId;

                fn expecting(&self, formatter: &mut fmt::Formatter) -> fmt::Result {
                    formatter.write_str("`secs` or `nanos`")
                }

                fn visit_str<E>(self, value: &str) -> Result<ReadAlignmentId, E>
                where
                    E: serde::de::Error,
                {
                    Ok(ReadAlignmentId::from_str(value).unwrap())
                }
            } // Visitor
            deserializer.deserialize_str(ReadAlignmentIdVisitor)
        }
    }
}
*/

impl fmt::Debug for ReadAlignmentId {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{self}")
    }
}

impl ReadAlignmentId {
    /// Construct from a name, which implies that there is only one alignment for this read.
    pub fn from_name(read_name: impl Into<String>) -> Self {
        let read_name = read_name.into();
        Self {
            read_name,
            //align_id: None,
        }
    }
    /// Construct a `ReadAlignmentId` from a name and a count. This is used with multiple alignments per read.
    pub fn from_name_and_id(read_name: impl Into<String>, id: impl Into<i32>) -> Self {
        let read_name = read_name.into();
        let _align_id = Some(id.into());
        Self {
            read_name,
            //align_id,
        }
    }
    #[must_use]
    pub fn unique_name(&self) -> String {
        /*
        if let Some(id) = self.align_id {
            format!("{}|{id}", self.read_name)
        } else {
            format!("{}", self.read_name)
        }
        */
        self.read_name.to_string()
    }
}

impl<T: Into<String>> std::convert::From<T> for ReadAlignmentId {
    fn from(x: T) -> ReadAlignmentId {
        ReadAlignmentId::from_name(x.into())
    }
}

/// Struct for phasing site.
#[derive(Clone, Debug, Default)]
pub struct PhasingSite {
    pub pos: i64,
    pub ref_base: String,
    pub var_base: String,
    pub deletion: String,
}

impl fmt::Display for PhasingSite {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}_{}_{}", self.pos, self.ref_base, self.var_base)
    }
}

impl PhasingSite {
    /// Create `PhasingSite` from an input string.
    /// Panics if malformed. Expects "{int}_{ref}_{var}".
    #[must_use]
    pub fn new(x: &str) -> Self {
        if let Some((pos, ref_base, var_base)) = x.split_terminator('_').next_tuple() {
            let pos = pos.parse::<i64>().unwrap();
            let ref_base = ref_base.to_owned();
            let var_base = var_base.to_owned();
            Self {
                pos,
                ref_base,
                var_base,
                deletion: String::new(),
            }
        } else if let Some((pos, var)) = x.split_terminator('_').next_tuple() {
            let pos = pos.parse::<i64>().unwrap();
            if var.starts_with("del") {
                Self {
                    pos,
                    deletion: var.to_owned(),
                    ..Default::default()
                }
            } else {
                panic!("var {var}");
            }
        } else {
            panic!("input {x}");
        }
    }
}

pub type ReadFingerprintMap = BTreeMap<ReadAlignmentId, VString>;

/// This struct contains the inputs to `VariantGraph` (`read_to_hap`),
/// and, if present, the assembled  haps in result json being parsed.
#[derive(Debug, Clone, Default)]
pub struct ParaphaseOutput {
    pub pivot_site: Option<i64>,
    pub gene: String,
    pub read_to_hap: ReadFingerprintMap,
    pub assembled_haps: Option<AssembledPaths>,
    pub final_haps: Option<BTreeMap<VString, String>>,
    pub sites_for_phasing: Vec<PhasingSite>,
}

type JsonMap = serde_json::Map<String, serde_json::Value>;

/// Takes an input json and fallibly extracts the `ReadFingerprintMap`.
fn parse_read_to_hap(value: &JsonMap) -> Result<ReadFingerprintMap, DError> {
    // log::trace!("value to parse: {value:?}");
    let read_details = value
        .get("read_details")
        .ok_or(Error::new("No read_details field"))?;
    let mut ret = ReadFingerprintMap::new();
    if let Some(read_details) = read_details.as_object() {
        for (read_name, hap) in read_details {
            ret.insert(
                ReadAlignmentId::from_name(read_name),
                VString::from(
                    hap.as_str()
                        .ok_or(Error::new("read_details hap was not a string type"))?,
                ),
            );
        }
    }
    Ok(ret)
}

fn parse_final_haps(value: &JsonMap) -> Option<BTreeMap<VString, String>> {
    value.get("final_haplotypes").and_then(|x| {
        // x is now an array type
        let mut ret = BTreeMap::<VString, String>::new();
        if let serde_json::Value::Object(dict) = x {
            for (key, val) in dict {
                let key = VString::from(key);
                let val = val.as_str()?;
                ret.insert(key, val.into());
            }
            Some(ret)
        } else {
            None
        }
    })
}

fn parse_hap_asm(value: &JsonMap) -> Option<AssembledPaths> {
    value.get("assembled_haplotypes").and_then(|x| {
        // x is now an array type
        if let serde_json::Value::Array(arr) = x {
            Some(AssembledPaths::from_seqs(arr.iter().map(|x| {
                x.as_str()
                    .expect("hap json did not have a string haplotype")
            })))
        } else {
            None
        }
    })
}

impl ParaphaseOutput {
    #[must_use]
    pub fn pivot_index(&self, pivot_site: Option<i64>) -> i64 {
        let pivot_site = pivot_site.or(self.pivot_site);
        pivot_site
            .and_then(|pos| {
                self.sites_for_phasing
                    .iter()
                    .position(|x| x.pos == pos)
                    .map(|x| x as i64)
            })
            .unwrap_or(-1)
    }
    #[must_use]
    pub fn new(read_to_hap: ReadFingerprintMap) -> Self {
        Self {
            read_to_hap,
            gene: String::from("NoGene"),
            ..Default::default()
        }
    }
    pub fn from_json(
        value: &serde_json::Value,
        gene: String,
        config: Option<&RegionConfig>,
    ) -> Result<Self, DError> {
        match value {
            serde_json::Value::Array(arr) => {
                Err(std::boxed::Box::new(Error::new(format!("Unexpected array type found. You should use a different function, as we expect a dictionary type. {value:?}, {arr:?}"))))
            }
            serde_json::Value::Object(dict) => {
                let read_to_hap = parse_read_to_hap(dict)?;
                let assembled_haps = parse_hap_asm(dict);
                let final_haps = parse_final_haps(dict);
                let sites_for_phasing: Vec<PhasingSite> = dict.get("sites_for_phasing").and_then(|x| x.as_array()).map(|x| x.iter().map(|x| PhasingSite::new(x.as_str().expect("Phasing site was not a string"))).collect::<Vec<_>>()).unwrap_or_default();

                let pivot_site = config.unwrap_or(&config::CONFIG).get(&gene[..]).and_then(|x| x.get("pivot_site").and_then(serde_yaml::Value::as_i64));
                // log::trace!("pivot site for gene {gene} is {pivot_site:?}");
                Ok(Self { pivot_site, gene, read_to_hap, assembled_haps, final_haps, sites_for_phasing })
            }
            x => {
                Err(std::boxed::Box::new(Error::new(format!("Unexpected value type: {x:?}"))))
            },
        }
    }

    /// Read a `ParaphaseOutput` struct from a `std::io::Read` object.
    /// # Errors
    /// 1. Malformatted json - could not parse.
    /// 2. Malformatted paraphase config - expected a json dictionary object.
    pub fn from_reader(
        reader: impl std::io::Read,
        config: Option<&RegionConfig>,
    ) -> Result<Self, DError> {
        Self::from_json(
            &serde_json::from_reader(reader)?,
            String::from("NoGene"),
            config,
        )
    }

    /// Read a `ParaphaseOutput` struct from a `&std::path::Path` object.
    /// Calls `from_reader` after opening.
    /// # Errors
    /// 1. Malformatted json - could not parse.
    /// 2. Malformatted paraphase config - expected a json dictionary object.
    pub fn from_path(
        path: &std::path::Path,
        config: Option<&RegionConfig>,
    ) -> Result<Self, DError> {
        Self::from_reader(std::io::BufReader::new(std::fs::File::open(path)?), config)
    }
}

/// Yields the json object from which it was parsed
/// as well as a map from the gene names to the assembly inputs.
/// The `ParaphaseOutput` can then be fed to `VariantGraph` for assembly.
#[derive(Debug, Clone)]
pub struct ParsedParaphaseOutputJSON {
    pub object: serde_json::Map<String, serde_json::Value>,
    pub gene_data: BTreeMap<String, ParaphaseOutput>,
}

impl ParsedParaphaseOutputJSON {
    /// # Errors
    /// 1. Mal-formatted json.
    pub fn from_reader(
        reader: impl std::io::Read,
        config: Option<&RegionConfig>,
    ) -> Result<Self, DError> {
        Self::from_json(serde_json::from_reader(reader)?, config)
    }

    /// # Errors
    /// 1. Fail to read from input path.
    /// 2. If ends with `.xz`, fail to decompress with `xz` executable.
    /// 3. Mal-formatted json.
    ///
    /// # Panics
    /// 1. If `xz` command fails to create `stdout`. This should never happen.
    pub fn from_path(
        path: impl Into<std::path::PathBuf>,
        config: Option<&RegionConfig>,
    ) -> Result<Self, DError> {
        let path = path.into();
        if path
            .extension()
            .map_or(false, |ext| ext.eq_ignore_ascii_case("xz"))
        {
            let res = std::process::Command::new("xz")
                .arg("-dc")
                .arg(path)
                .stdout(std::process::Stdio::piped())
                .spawn()?;
            Self::from_reader(std::io::BufReader::new(res.stdout.unwrap()), config)
        } else {
            Self::from_reader(std::io::BufReader::new(std::fs::File::open(path)?), config)
        }
    }

    /// # Errors
    /// 1. Unexpected json format - expected a dictionary type but found something else.
    pub fn from_json(
        object: serde_json::Value,
        config: Option<&RegionConfig>,
    ) -> Result<Self, DError> {
        let serde_json::Value::Object(object) = object else {
            return Err(Box::new(Error::new(
                "Not a dictionary for ParsedParaphaseOutputJSON",
            )));
        };
        log::debug!("Parsing from_json from json {object:?}");
        Self::from_object(object, config)
    }

    /// From `serde_json::Object`, which is a dictionary, build a specific config.
    /// # Errors
    /// 1. Unexpected json format.
    pub fn from_object(object: JsonMap, config: Option<&RegionConfig>) -> Result<Self, DError> {
        let mut gene_data = BTreeMap::new();
        for (key, value) in object.iter().map(|(key, value)| {
            (
                String::from(key),
                ParaphaseOutput::from_json(value, String::from(key), config),
            )
        }) {
            gene_data.insert(key, value?);
        }
        Ok(Self { object, gene_data })
    }
}

pub fn write_outputs(
    x: &BTreeMap<String, GeneCall>,
    writer: &mut impl std::io::Write,
) -> Result<(), DError> {
    writeln!(writer, "{}", serde_json::to_string_pretty(x)?)?;
    Ok(())
}
