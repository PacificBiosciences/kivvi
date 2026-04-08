use crate::assembly::variant_graph::VGraph;
use crate::detail::low_complexity::LowConfidenceSites;

use vstr::VString;

use indexmap::IndexSet;
use itertools::Itertools;
use petgraph::{stable_graph::NodeIndex, Direction::Outgoing};
use rust_htslib::{bam, htslib};
use std::{
    collections::BTreeMap,
    error,
    path::{Path, PathBuf},
};

pub type DError = std::boxed::Box<dyn std::error::Error>;
pub type DResult = Result<(), DError>;

// NOTE: I think this reaches a nice middle ground between Result<(), DError>; and Result<T, String>;
// TODO: Replace DResult and DError use
pub type MResult<T> = std::result::Result<T, Box<dyn error::Error>>;

pub type HashMap<K, V> = fnv::FnvHashMap<K, V>;
pub type HashSet<K> = fnv::FnvHashSet<K>;

/// `NotImplementedError`
/// Thrown when code that has not been implemented is executed.
#[derive(Debug, Clone, Default)]
pub struct NotImplementedError {
    pub msg: String,
}

impl std::error::Error for NotImplementedError {}

impl std::fmt::Display for NotImplementedError {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "NotImplementedError{{msg: {}}}", self.msg)
    }
}

impl NotImplementedError {
    #[must_use]
    pub fn new(x: impl Into<String>) -> Self {
        Self { msg: x.into() }
    }
}

/// A port of `all_simple_paths` from `petgraph`, but skipping deleted nodes.
/// Necessary to support dynamic graph.
///
/// # Panics
/// 1. weight for child in graph is not present. (Should never happen.)
#[must_use]
pub fn all_simple_paths(
    graph: &VGraph,
    from: NodeIndex<u32>,
    to: NodeIndex<u32>,
    min_intermediate_nodes: usize,
    max_intermediate_nodes: Option<usize>,
) -> Vec<Vec<NodeIndex<u32>>> {
    type NodeId = NodeIndex<u32>;
    // how many nodes are allowed in simple path up to target node
    // it is min/max allowed path length minus one, because it is more appropriate when implementing lookahead
    // than constantly add 1 to length of current path
    let max_length = if let Some(l) = max_intermediate_nodes {
        l + 1
    } else {
        graph.node_count() - 1
    };

    let min_length = min_intermediate_nodes + 1;

    // list of visited nodes
    let mut visited: IndexSet<NodeId> = IndexSet::from_iter(Some(from));
    // list of childs of currently exploring path nodes,
    // last elem is list of childs of last visited node
    let mut stack = vec![graph.neighbors_directed(from, Outgoing)];

    std::iter::from_fn(move || {
        while let Some(children) = stack.last_mut() {
            // Here is where we skip the deleted nodes.
            let mut get_next = || -> Option<petgraph::stable_graph::NodeIndex<_>> {
                let mut ret = None;
                for res in children.by_ref() {
                    if !graph.node_weight(res).unwrap().is_del() {
                        ret = Some(res);
                        break;
                    }
                }
                ret
            };
            if let Some(child) = get_next() {
                if visited.len() < max_length {
                    if child == to {
                        if visited.len() >= min_length {
                            let path = visited.iter().copied().chain(Some(to)).collect::<Vec<_>>();
                            return Some(path);
                        }
                    } else if !visited.contains(&child) {
                        visited.insert(child);
                        stack.push(graph.neighbors_directed(child, Outgoing));
                    }
                } else {
                    if (child == to || children.any(|v| v == to)) && visited.len() >= min_length {
                        let path = visited.iter().copied().chain(Some(to)).collect::<Vec<_>>();
                        return Some(path);
                    }
                    stack.pop();
                    visited.pop();
                }
            } else {
                stack.pop();
                visited.pop();
            }
        }
        None
    })
    .collect::<Vec<_>>()
}

#[must_use]
pub fn test_file(x: &str) -> std::path::PathBuf {
    let mut path = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    for s in &["tests", "data", x] {
        path.push(s);
    }
    path
}

// from https://stackoverflow.com/questions/35045996/check-if-a-command-is-in-path-executable-as-process
#[must_use]
pub fn get_program_from_path(program: &str) -> Option<std::path::PathBuf> {
    std::env::var("PATH").ok().and_then(|path| {
        path.split_terminator(':')
            .map(std::path::PathBuf::from)
            .map(|p| p.join(program))
            .find(|p| std::fs::metadata(p).is_ok())
    })
}

pub fn init_log(level: log::LevelFilter) {
    if let Err(e) = env_logger::builder()
        .format_timestamp_millis()
        .filter_level(level)
        .try_init()
    {
        log::debug!("Logger already activated. Error: {e:?}");
    }
}

/// Read an indexed bam file. This can take a url or a local path.
/// If the index is not directly adjacent to the bam file, use `read_indexed_bam_with_index`.
/// # Errors
/// 1. `rust_htslib::errors::Error` if file not found or corrupted.
pub fn read_indexed_bam(
    input: impl Into<String>,
) -> Result<bam::IndexedReader, rust_htslib::errors::Error> {
    let input = input.into();
    if let Ok(url) = url::Url::parse(&input) {
        bam::IndexedReader::from_url(&url)
    } else {
        bam::IndexedReader::from_path(&input)
    }
}

///
/// Terribly unsafe hack.
/// `rust_htslib` doesn't expose the inner, so parity with Paraphase isn't possible without this.
///
/// This comes from pysam's base-quality filtering in the pileups, which uses the next base position after a deletion.
/// in `rust-htslib`, it returns `Option<usize>` and `None` if `is_del`.
/// This function lets us get `qpos` from `bam_pileup1_t` even if the read is a deletion.
///```
/// use paraphase::detail::util::{raw_qpos, test_file};
/// use rust_htslib::bam::Read;
/// let mut reader =
///     rust_htslib::bam::Reader::from_path(test_file("HG00733_smn1_realigned.bam")).unwrap();
/// for pileup in reader.pileup() {
///     let pileup = pileup.unwrap();
///     let qposes = pileup.alignments().map(|x| x.qpos()).collect::<Vec<_>>();
///     let qraws = pileup
///         .alignments()
///         .map(|x| raw_qpos(&x))
///         .collect::<Vec<_>>();
///     for (qp, qr) in qposes.iter().zip(qraws.iter()) {
///         assert!(
///             qp.is_none() || qp.unwrap() == *qr,
///             "qpos: {qposes:?} qr: {qraws:?}"
///         );
///     }
/// }
/// ```
#[must_use]
pub fn raw_qpos<'a>(x: &'a bam::pileup::Alignment<'a>) -> usize {
    static_assertions::const_assert_eq!(
        std::mem::size_of::<&bam::pileup::Alignment<'_>>(),
        std::mem::size_of::<&htslib::bam_pileup1_t>()
    );
    let ptr = (x as *const bam::pileup::Alignment<'_>).cast::<&htslib::bam_pileup1_t>();
    unsafe { *ptr }.qpos as usize
}

/*
Not possible currently in rust-htslib due to encapsulation.

/// Read an indexed bam file. This can take a url or a local path to either or both bam and index files.
/// # Errors
/// 1. `rust_htslib::errors::Error` if file not found or corrupted.
/// 2. `std::ffi::CString conversion` failure.
pub fn read_indexed_bam_with_index(
    input: impl Into<String>,
    index_path: impl Into<String>,
) -> Result<bam::IndexedReader, DError> {
    use rust_htslib::htslib::{self, hts_open};
    use std::{ffi, rc::Rc, str};

    let input = input.into();
    let index_path = index_path.into();
    let htsfile = hts_open(&input, b"r");
    let header = unsafe { htslib::sam_hdr_read(htsfile) };
    let c_str_path = ffi::CString::new(input)?;
    let c_str_index_path = ffi::CString::new(index_path)?;
    let idx =
        unsafe { htslib::sam_index_load2(htsfile, c_str_path.as_ptr(), c_str_index_path.as_ptr()) };
    if idx.is_null() {
        Err(Box::new(rust_htslib::errors::Error::BamInvalidIndex {
            target: str::from_utf8(path)?.to_owned(),
        }))
    } else {
        Ok(IndexedReader {
            htsfile,
            header: Rc::new(bam::HeaderView::new(header)),
            idx: Rc::new(bam::IndexView::new(idx)),
            itr: None,
            tpool: None,
        })
    }
}
*/

pub type GeneLevelConfig = BTreeMap<String, serde_yaml::Value>;

lazy_static::lazy_static! {
    pub static ref CRATE_NAME: String = std::env::var("CARGO_PKG_NAME").unwrap().to_string();
    //pub static ref FULL_VERSION: String = format!("{}", env!("CRATE_GIT_SHA"));
    //pub static ref FULL_VERSION_PROGRAM: String =
    //    format!("{}-{}",
    //    env!("CARGO_PKG_NAME"),
    //    env!("CRATE_GIT_SHA"));
}

/// Read a bam file. This can take a url or a local path.
/// # Errors
/// 1. `rust_htslib::errors::Error` if file not found or corrupted.
pub fn read_bam(input: impl AsRef<Path>) -> Result<bam::Reader, String> {
    let input = input.as_ref().to_string_lossy().into_owned();
    if let Ok(url) = url::Url::parse(&input) {
        bam::Reader::from_url(&url)
    } else if input == "-" || input == "/dev/stdin" {
        bam::Reader::from_stdin()
    } else {
        bam::Reader::from_path(&input)
    }
    .map_err(|x| format!("Error: {x}"))
}

/// Count records in a bam file.
pub fn count_records(x: &mut impl bam::Read) -> usize {
    x.rc_records().count()
}

/// Accesses all sample names from a bam header.
/// Since a given sample can have several RG tags, we return a map
/// from sample to a vector of matching RG tags.
///
///```
/// use paraphase::detail::util::{sample_names, test_file};
/// use rust_htslib::bam::{Header, Read, Reader};
/// let bam = Reader::from_path(test_file("header_only.bam")).unwrap();
/// let expected_names = [("UnnamedSample", "default")]
///     .into_iter()
///     .map(|(key, val)| (key.to_string(), vec![val.to_string()]))
///     .collect::<std::collections::BTreeMap<_, _>>();
/// assert_eq!(
///     sample_names(&Header::from_template(bam.header())),
///     expected_names,
/// );
/// ```
#[must_use]
pub fn sample_names(view: &bam::Header) -> BTreeMap<String, Vec<String>> {
    let mut ret = BTreeMap::<String, Vec<_>>::new();
    if let Some(map) = view.to_hashmap().remove("RG") {
        for (sample, id) in map.iter().filter_map(|x| {
            x.get("SM")
                .and_then(|sample| x.get("ID").map(|id| (sample.clone(), id.clone())))
        }) {
            ret.entry(sample).or_default().push(id);
        }
    }
    ret
}

#[must_use]
pub fn sample_names_from_input(x: &PathBuf) -> Option<BTreeMap<String, Vec<String>>> {
    use bam::Read;
    let reader = read_bam(x).ok()?;
    Some(sample_names(&bam::Header::from_template(reader.header())))
}

/*
 * BAM_CIGAR_TYPE  QUERY  REFERENCE
 * --------------------------------
 * BAM_CMATCH      1      1
 * BAM_CINS        1      0
 * BAM_CDEL        0      1
 * BAM_CREF_SKIP   0      1
 * BAM_CSOFT_CLIP  1      0
 * BAM_CHARD_CLIP  0      0
 * BAM_CPAD        0      0
 * BAM_CEQUAL      1      1
 * BAM_CDIFF       1      1
 * BAM_CBACK       0      0
 */
#[must_use]
#[inline]
pub fn consumes_ref(x: bam::record::Cigar) -> bool {
    use bam::record::Cigar::{Del, Diff, Equal, Match, RefSkip};
    matches!(x, Del(_) | RefSkip(_) | Match(_) | Diff(_) | Equal(_))
}

#[must_use]
#[inline]
pub fn consumes_qry(x: bam::record::Cigar) -> bool {
    use bam::record::Cigar::{Diff, Equal, Ins, Match, SoftClip};
    matches!(x, Ins(_) | SoftClip(_) | Match(_) | Diff(_) | Equal(_))
}

#[must_use]
pub fn crate_name() -> &'static str {
    &CRATE_NAME[..]
}

#[must_use]
pub fn load_all_seqs(index: &rust_htslib::faidx::Reader) -> Vec<vstr::VString> {
    load_all_seqs_view(index)
}

pub fn seq_name_pairs(
    path: &std::path::Path,
    uppercase: bool,
) -> Result<Vec<(VString, VString)>, DError> {
    use needletail::Sequence;
    let mut file = needletail::parse_fastx_file(path)?;
    let mut ret = Vec::with_capacity(8);
    while let Some(record) = file.next() {
        let record = record?;
        let mut sequence = VString::from(&record.sequence().strip_returns()[..]);
        if uppercase {
            sequence.make_ascii_uppercase();
        }
        let name = VString::from(record.id());
        ret.push((name, sequence));
    }
    Ok(ret)
}

#[must_use]
pub fn faidx_names(index: &rust_htslib::faidx::Reader) -> Vec<vstr::VString> {
    let mut ret = vec![];
    for v in 0..index.n_seqs() as i32 {
        ret.push(index.seq_name(v).expect("Missing seq name").into());
    }
    ret
}

#[must_use]
pub fn load_all_seqs_view(index: &rust_htslib::faidx::Reader) -> Vec<VString> {
    let seq_names =
        (0..index.n_seqs()).map(|x| index.seq_name(x as i32).expect("Missing seq name"));
    seq_names
        .map(|x| {
            index
                .fetch_seq(&x, 0, i64::MAX as usize)
                .expect("Failed to load seq from faidx")
        })
        .map(|x| VString::from(x))
        .collect::<Vec<VString>>()
}

#[must_use]
/// Parse homopolymer sites from a text file.
/// Format: "{pos}\t{characters}"
/// # Panics
/// If the position is not an integer.
pub fn parse_homopolymers(path: &std::path::Path) -> LowConfidenceSites {
    LowConfidenceSites::from_map(
        slurp::iterate_all_lines(path)
            .map(|x| x.expect("Line is not utf-8' in agap9-hpol-expected"))
            .map(|x| {
                let (k, v) = x.split_terminator('\t').next_tuple().unwrap();
                (
                    k.parse::<i64>().unwrap(),
                    v.split_terminator(',')
                        .map(|x| x.chars().next().unwrap() as u8)
                        .collect::<linear_map::set::LinearSet<_>>(),
                )
            })
            .collect::<BTreeMap<_, _>>(),
    )
}

pub trait DeletionInsensitiveCompare {
    /// Computes if two fingerprints are matches apart from gaps/deletions and soft-clips.
    /// Implemented for structures containing bytes.
    #[must_use]
    fn same_without_dels(&self, other: &Self) -> bool;
    #[must_use]
    fn is_del(x: u8) -> bool {
        matches!(x, b'x' | b'0')
    }
}

impl DeletionInsensitiveCompare for str {
    fn same_without_dels(&self, other: &Self) -> bool {
        self.as_bytes().same_without_dels(&other.as_bytes())
    }
}

impl DeletionInsensitiveCompare for &[u8] {
    /// This may need extension for leading xs followed by softclips,
    /// but let's hope it is rare e.g., " xxx0024".
    fn same_without_dels(&self, other: &Self) -> bool {
        let possible_match = self.len() == other.len()
            && self
                .iter()
                .zip(other.iter())
                .all(|(x, y)| x == y || Self::is_del(*x) || Self::is_del(*y));
        if !possible_match {
            return false;
        }
        let mut self_g = self
            .iter()
            .group_by(|x| *x)
            .into_iter()
            .map(|(k, v)| (*k, v.count() as i32))
            .collect::<Vec<_>>();
        let mut other_g = other
            .iter()
            .group_by(|x| *x)
            .into_iter()
            .map(|(k, v)| (*k, v.count() as i32))
            .collect::<Vec<_>>();
        let trim_ends = |x: &mut Vec<(u8, i32)>| {
            for idx in (0..x.len()).rev() {
                let reg = &x[idx];
                if !Self::is_del(reg.0) {
                    break;
                }
                x.pop();
            }
            let end_pos = x.iter().position(|x| !Self::is_del(x.0)).unwrap_or(x.len());
            x.drain(0..end_pos);
            /*
            while x.len() > 0 && Self::is_del(x.first().unwrap().0) {
                x.swap_remove(0);
            }
            */
        };
        trim_ends(&mut self_g);
        trim_ends(&mut other_g);
        let no_remaining_clips = |x: &[(u8, i32)]| -> bool {
            [x.last(), x.first()].iter().flatten().all(|x| x.0 != b'0')
        };
        [self_g, other_g].iter().all(|x| no_remaining_clips(x))
    }
}

impl<const N: usize> DeletionInsensitiveCompare for [u8; N] {
    fn same_without_dels(&self, other: &Self) -> bool {
        (&self[..]).same_without_dels(&&other[..])
    }
}

impl DeletionInsensitiveCompare for vstr::VStr<'_> {
    fn same_without_dels(&self, other: &Self) -> bool {
        (&self[..]).same_without_dels(&&other[..])
    }
}

impl DeletionInsensitiveCompare for vstr::VString {
    fn same_without_dels(&self, other: &Self) -> bool {
        self.vstr().same_without_dels(&other.vstr())
    }
}

#[cfg(test)]
mod tests {
    use crate::detail::util::{DResult, DeletionInsensitiveCompare};

    #[test]
    fn test_config_ok() {
        let res = crate::config::Region::load(None);
        assert_eq!(res["smn1"]["genes"].as_str(), Some("SMN1,SMN2"));
        assert_eq!(res["smn1"]["pivot_site"].as_i64(), Some(70_951_946));
        assert_eq!(res["CYP2D6"]["genes"].as_str().unwrap(), "CYP2D6");
    }

    #[test]
    fn test_gene_config_parser_ok() -> DResult {
        use std::collections::BTreeSet;
        let conf = crate::config::Gene::try_load(None)?;
        assert!(conf.genes_to_call.is_empty());
        assert_eq!(
            conf.no_vcf_genes,
            ["CFH", "CFHR3"]
                .into_iter()
                .map(String::from)
                .collect::<BTreeSet<_>>()
        );
        assert_eq!(
            conf.two_reference_regions_genes,
            ["smn1", "pms2", "strc", "ikbkg", "ncf1"]
                .into_iter()
                .map(String::from)
                .collect::<BTreeSet<_>>()
        );
        assert_eq!(
            conf.no_genome_depth_genes,
            ["pms2", "neb", "cfc1", "ikbkg", "opn1lw", "rccx"]
                .into_iter()
                .map(String::from)
                .collect::<BTreeSet<_>>()
        );
        Ok(())
    }

    #[test]
    fn same_without_dels_ok() {
        use DeletionInsensitiveCompare;
        // Simple case: matches
        assert!("12111".same_without_dels("12x11"));
        assert!(b"12111".same_without_dels(b"12x11"));
        let lhs = vstr::VString::from("121x1x11");
        let rhs = vstr::VString::from("121211x1");
        assert!(lhs.same_without_dels(&rhs));

        // Mismatched length: should be false
        let lhs = vstr::VString::from("121x1x11");
        let rhs = vstr::VString::from("121211x1a");
        assert!(!lhs.same_without_dels(&rhs));

        // Mismatched: mismatches
        //                                   *
        let lhs = vstr::VString::from("121x1x11");
        let rhs = vstr::VString::from("12121121");
        assert!(!lhs.same_without_dels(&rhs));

        // Allow softclip at end
        let lhs = vstr::VString::from("121x1x11");
        let rhs = vstr::VString::from("12121000");
        assert!(lhs.same_without_dels(&rhs));

        // Do not allow internal softclip.
        let lhs = vstr::VString::from("121x1x11");
        let rhs = vstr::VString::from("12121000");
        assert!(lhs.same_without_dels(&rhs));

        // Do allow internal gaps.
        let lhs = vstr::VString::from("121x1x11");
        let rhs = vstr::VString::from("12x21000");
        assert!(lhs.same_without_dels(&rhs));
    }
}
