mod args;
mod io_utils;
mod matcher;
mod pipeline;
mod samplesheet;

use anyhow::{Context, Result};
use mimalloc::MiMalloc;
use std::fs;
use std::time::Instant;

use crate::args::{parse_cli_with_legacy_run_compat, Command, RunArgs};
use crate::matcher::{build_file_pair_index, collect_candidates};
use crate::pipeline::{build_sample_jobs, process_jobs_to_final_outputs};
use crate::samplesheet::read_samplesheet;

#[global_allocator]
static GLOBAL: MiMalloc = MiMalloc;

fn main() -> Result<()> {
    let cli = parse_cli_with_legacy_run_compat();

    match cli.command {
        Command::Run(args) => run(args, false),
        Command::Benchmark(args) => run(args, true),
    }
}

fn run(args: RunArgs, benchmark: bool) -> Result<()> {
    if !args.output_dir.exists() {
        fs::create_dir_all(&args.output_dir).context("Failed to create output directory")?;
    }

    let ext = "fastq.gz";
    let out_r1_path = args.output_dir.join(format!("demultiplexed.R1.{}", ext));
    let out_r2_path = args.output_dir.join(format!("demultiplexed.R2.{}", ext));

    if out_r1_path.exists() || out_r2_path.exists() {
        return Err(anyhow::anyhow!(
            "Output files already exist in {:?}. Please remove them or use a different output directory.",
            args.output_dir
        ));
    }

    println!("Output R1: {:?}", out_r1_path);
    println!("Output R2: {:?}", out_r2_path);

    let t0 = Instant::now();
    let candidates = collect_candidates(&args.input_dir);
    let t_scan = t0.elapsed();

    let t1 = Instant::now();
    let parsed_sheet = read_samplesheet(&args.samplesheet)?;
    let t_sheet = t1.elapsed();
    println!(
        "Processing {} samples in parallel...",
        parsed_sheet.samples.len()
    );

    let t2 = Instant::now();
    let file_pair_index = build_file_pair_index(&candidates, &parsed_sheet.barcode_ids);
    let sample_jobs = build_sample_jobs(parsed_sheet.samples, &file_pair_index)?;
    let t_match = t2.elapsed();

    println!("\nResolved sample inputs:");
    for job in &sample_jobs {
        println!(
            "  - sample={} barcode={} R1={} R2={}",
            job.barcode_id,
            job.barcode_seq,
            job.r1_file.display(),
            job.r2_file.display()
        );
    }

    let t3 = Instant::now();
    let stats = process_jobs_to_final_outputs(
        &sample_jobs,
        &out_r1_path,
        &out_r2_path,
        args.gzip_level,
        args.threads,
    )?;
    let t_process = t3.elapsed();

    println!("Done. Merged files created.");

    if benchmark {
        let total = t_scan + t_sheet + t_match + t_process;
        println!("\nBenchmark:");
        println!("  - Scan input files:   {:.3}s", t_scan.as_secs_f64());
        println!("  - Parse samplesheet:  {:.3}s", t_sheet.as_secs_f64());
        println!("  - Match file pairs:   {:.3}s", t_match.as_secs_f64());
        println!("  - Process + write:    {:.3}s", t_process.as_secs_f64());
        println!("  - Total:              {:.3}s", total.as_secs_f64());
        println!("  - Read pairs written: {}", stats.read_pairs);
        if total.as_secs_f64() > 0.0 {
            println!(
                "  - Throughput:         {:.2} read pairs/s",
                stats.read_pairs as f64 / total.as_secs_f64()
            );
        }
    }

    Ok(())
}
