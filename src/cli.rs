use chrono::Datelike;
use clap::{ArgGroup, Parser, Subcommand};
use log::{error, info};
use std::path::PathBuf;

#[derive(Parser)]
#[clap(author,
    version = env!("CARGO_PKG_VERSION"),
    about,
    after_help = format!("Copyright (C) 2004-{}     Pacific Biosciences of California, Inc.
This program comes with ABSOLUTELY NO WARRANTY; it is intended for
Research Use Only and not for use in diagnostic procedures.", chrono::Utc::now().year()))]
pub struct Settings {
    #[command(subcommand)]
    pub command: Command,

    /// Input alignment file in BAM format
    #[clap(required = true)]
    #[clap(short = 'b')]
    #[clap(long = "bam")]
    #[clap(value_name = "BAM")]
    #[clap(help_heading = Some("Input/Output"))]
    pub bam_filename: PathBuf,

    /// Output directory
    #[clap(required = true)]
    #[clap(short = 'o')]
    #[clap(long = "out")]
    #[clap(value_name = "OUTDIR")]
    #[clap(help_heading = Some("Input/Output"))]
    pub outdir: PathBuf,

    /// Output prefix
    #[clap(required = true)]
    #[clap(short = 'p')]
    #[clap(long = "prefix")]
    #[clap(value_name = "PREFIX")]
    #[clap(help_heading = Some("Input/Output"))]
    pub prefix: String,

    /// Path to reference genome FASTA
    #[clap(required = true)]
    #[clap(short = 'r')]
    #[clap(long = "reference")]
    #[clap(value_name = "FASTA")]
    #[clap(help_heading = Some("Input/Output"))]
    pub reference: PathBuf,

    /// sensitive mode
    #[clap(long, action)]
    #[clap(help = "If specified, paraphase will use lower cutoff for fingerprint")]
    pub sensitive: bool,

    /// Optionally specify a list of variants (SNVs only) for allele plotting
    #[clap(short = 'l')]
    #[clap(long)]
    #[clap(hide = true)]
    pub variant_list: Option<PathBuf>,

    /// do not run paraphase
    #[clap(long, action)]
    #[clap(hide = true)]
    #[clap(help = "If specified, will not run paraphase to phase up/downstream regions")]
    pub nopp: bool,

    /// Enable verbose output
    #[clap(short = 'v')]
    #[clap(long = "verbose")]
    #[clap(action = clap::ArgAction::Count, help = "Specify multiple times to increase verbosity level (e.g., -vv for more verbosity)")]
    pub verbosity: u8,
}

#[derive(Subcommand)]
pub enum Command {
    #[clap(about = "LPA KIV2 repeat")]
    Kiv2(KivArgs),
    #[clap(about = "D4Z4 repeat")]
    D4z4(D4z4Args),
}

#[derive(Parser, Debug)]
#[command(group(ArgGroup::new("kiv2")))]
#[command(arg_required_else_help(false))]
pub struct KivArgs {}

#[derive(Parser, Debug)]
#[command(group(ArgGroup::new("d4z4")))]
#[command(arg_required_else_help(false))]
pub struct D4z4Args {}

/// Parse settings
pub fn get_raw_settings() -> Settings {
    Settings::parse()
}

/// Checks if required files exist
pub fn check_settings(settings: Settings) -> Settings {
    if !&settings.bam_filename.exists() {
        error!(
            "Alignment file does not exist: \"{}\"",
            &settings.bam_filename.display()
        );
        std::process::exit(1);
    } else {
        info!("Alignment file: \"{}\"", &settings.bam_filename.display());
    }
    settings
}
