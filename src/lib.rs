/// Allele assembly
pub mod assembly;
/// Working with BAM input/output
pub mod bam_operation;
/// Main callers for each locus
pub mod caller;
/// CLI functionality and checks
pub mod cli;
/// functions specific to D4Z4
pub mod d4z4;
/// Depth based calls
pub mod depth;
/// Check methylation
pub mod methylation;
/// Plot assembled alleles
pub mod plot;
/// functions specific to read filtering
pub mod read_filtering;
/// WFA graph realignment adapted from HiPhase
pub mod realignment;
/// Fingerprint analysis
pub mod repeat_unit;
/// Utilities and locus specific parameters
pub mod util;
/// Variant calling
pub mod variant;
/// Write to vcf
pub mod vcf;
