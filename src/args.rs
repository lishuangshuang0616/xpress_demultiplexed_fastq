use clap::{Parser, Subcommand};
use std::ffi::OsString;
use std::path::PathBuf;

#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
pub struct Cli {
    #[command(subcommand)]
    pub command: Command,
}

#[derive(Subcommand, Debug)]
pub enum Command {
    /// Run FASTQ merge/transformation
    Run(RunArgs),
    /// Run and print per-stage timing metrics
    Benchmark(RunArgs),
}

/// Merge and modify FASTQ files based on a samplesheet.
#[derive(Parser, Debug, Clone)]
pub struct RunArgs {
    /// Input directory containing FASTQ files
    #[arg(short, long)]
    pub input_dir: PathBuf,

    /// Path to the samplesheet (CSV/TSV).
    /// Column 1: Barcode ID (matches filename).
    /// Column 2: Barcode Sequence (to append to R2).
    #[arg(short, long)]
    pub samplesheet: PathBuf,

    /// Gzip compression level (0-9).
    #[arg(long, default_value_t = 3)]
    pub gzip_level: u32,

    /// Output directory where merged files will be placed.
    #[arg(short, long, default_value = "output")]
    pub output_dir: PathBuf,

    /// Number of threads to use. Default is 8.
    #[arg(short, long, default_value = "8")]
    pub threads: usize,
}

pub fn parse_cli_with_legacy_run_compat() -> Cli {
    // Backward compatibility:
    // old style: fastq_merger --input-dir ...
    // new style: fastq_merger run --input-dir ...
    let mut argv: Vec<OsString> = std::env::args_os().collect();
    if argv.len() > 1 {
        let first = argv[1].to_string_lossy();
        if first.starts_with('-') {
            argv.insert(1, OsString::from("run"));
        }
    }
    Cli::parse_from(argv)
}
