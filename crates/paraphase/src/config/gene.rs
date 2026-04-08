use std::collections::BTreeSet;

use crate::detail::util::DError;

/// Configuration specifying logic by gene name.
#[derive(Clone, Debug, Default)]
pub struct Config {
    pub no_genome_depth_genes: BTreeSet<String>, // Don't check genomic depth when calling
    pub no_vcf_genes: BTreeSet<String>,          // Don't write a vcf
    pub genes_to_call: BTreeSet<String>,         // Specify which genes to analyze in this run
    pub check_sex_genes: BTreeSet<String>,       // Which genes require checking sample sex.
    pub two_reference_regions_genes: BTreeSet<String>, // Which genes have multiple reference regions to resolve.
}

impl Config {
    /// Builds `Config` from a file location.
    ///
    /// # Errors
    ///
    /// Returns `Err(config::Error)` on error.
    ///
    /// Malformatted yaml.
    /// Missing file at given path or not a file.
    pub fn try_from_path(x: impl Into<std::path::PathBuf>) -> Result<Self, DError> {
        Self::try_load(Some(&std::fs::read(x.into())?))
    }

    /// Parse a `gene::Config` from a path.
    /// Used to provide specific configurations.
    ///
    /// # Panics
    /// Malformatted yaml.
    /// Missing file at given path or not a file.
    #[must_use]
    pub fn from_path(x: impl Into<std::path::PathBuf>) -> Self {
        Self::try_from_path(x).expect("Failed to load gene config from path")
    }

    /// Try to yaml data from a byte slice.
    /// # Errors
    /// Errors if fails to parse.
    pub fn try_load(data: Option<&[u8]>) -> Result<Self, DError> {
        load(data)
    }

    /// Load yaml data from a byte slice.
    /// # Panics
    /// Errors if fails to parse.
    #[must_use]
    pub fn load(data: Option<&[u8]>) -> Self {
        Self::try_load(data).expect("Failed to parse gene config.")
    }

    ///
    /// Determines if a gene uses background depth calculation.
    #[must_use]
    pub fn uses_depth(&self, x: &str) -> bool {
        !self.no_genome_depth_genes.contains(x)
    }
}

#[must_use]
fn serde_value_to_stringset(x: &serde_yaml::Sequence) -> BTreeSet<String> {
    x.iter()
        .map(|x| x.as_str().expect("gene should be a string").to_string())
        .collect::<_>()
}

/// Parse gene confug from `&[u8]` slice.
/// To read from a file, load the file or mmap it.
/// If `None` is provided, use a genes yaml embedded in the executable.
///
/// # Errors
/// If the slice in invalid yaml.
///
/// # Panics
/// • If the parsed yaml is not a dictionary at the root level.
/// • If the keys are not string.
/// • If the values are not arrays of strings.
/// • If there is an unexpected key in the dictionary.
fn load(data: Option<&[u8]>) -> Result<Config, DError> {
    //const DATA: &[u8] =
    //    std::include_bytes!(concat!(env!("CARGO_MANIFEST_DIR"), "/data/genes.yaml"));
    //let data = data.unwrap_or(DATA);
    let data = data.unwrap();
    let data = serde_yaml::from_slice::<serde_yaml::Value>(data)?;
    let data = data
        .as_mapping()
        .expect("gene config yaml was not a mapping.");
    let mut ret = Config::default();
    let mut genome_depth_genes = BTreeSet::new();
    for (k, v) in data.into_iter().map(|(key, value)| {
        (
            key.as_str().expect("key was not a string"),
            value
                .as_sequence()
                .expect("Gene config should always be an array"),
        )
    }) {
        let stringset = || serde_value_to_stringset(v);
        match k {
            "genes_to_call" => {
                ret.genes_to_call = stringset();
            }
            "check_sex_genes" => {
                ret.check_sex_genes = stringset();
            }
            "no_vcf_genes" => {
                ret.no_vcf_genes = stringset();
            }
            "genome_depth_genes" => {
                genome_depth_genes = stringset();
            }
            "no_genome_depth_genes" => {
                ret.no_genome_depth_genes = stringset();
            }
            "two_reference_regions_genes" => {
                ret.two_reference_regions_genes = stringset();
            }
            _ => {
                panic!("Unexpected field {k} in gene config.");
            }
        }
    }

    let all_genes_to_call = if ret.genes_to_call.is_empty() {
        [
            &ret.check_sex_genes,
            &ret.no_vcf_genes,
            &ret.two_reference_regions_genes,
        ]
        .iter()
        .flat_map(|x| x.iter())
        .cloned()
        .collect::<BTreeSet<_>>()
    } else {
        ret.genes_to_call.clone()
    };

    // Account for old paraphase using genome_depth_genes and new paraphase using no_genome_depth_genes.
    // Use set subtraction to get parity.
    if !genome_depth_genes.is_empty() {
        ret.no_genome_depth_genes = all_genes_to_call
            .difference(&genome_depth_genes)
            .cloned()
            .collect::<_>();
    }

    Ok(ret)
}
