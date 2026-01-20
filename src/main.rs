use anyhow::{Context, Result};
use bio::io::fastq;
use bio::io::fastq::FastqRead;
use clap::Parser;
use flate2::read::MultiGzDecoder;
use flate2::write::GzEncoder;
use flate2::Compression;
use rayon::prelude::*;
use std::fs::File;
use std::io::{BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};
use walkdir::WalkDir;
use std::fs;

/// Merge and modify FASTQ files based on a samplesheet.
#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    /// Input directory containing FASTQ files
    #[arg(short, long)]
    input_dir: PathBuf,

    /// Path to the samplesheet (CSV/TSV). 
    /// Column 1: Barcode ID (matches filename). 
    /// Column 2: Barcode Sequence (to append to R2).
    #[arg(short, long)]
    samplesheet: PathBuf,

    /// Compress output with Gzip
    #[arg(short, long)]
    compress: bool,

    /// Output directory where merged files and temporary chunks will be placed.
    #[arg(short, long, default_value = "output")]
    output_dir: PathBuf,

    /// Number of threads to use. Default is 0 (auto-detect all available cores).
    #[arg(short, long, default_value = "0")]
    threads: usize,
}

// Enum to handle both plain and compressed writers uniformly
enum OutputWriter {
    Plain(BufWriter<File>),
    Compressed(Box<GzEncoder<BufWriter<File>>>),
}

impl Write for OutputWriter {
    fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
        match self {
            OutputWriter::Plain(w) => w.write(buf),
            OutputWriter::Compressed(w) => w.write(buf),
        }
    }

    fn flush(&mut self) -> std::io::Result<()> {
        match self {
            OutputWriter::Plain(w) => w.flush(),
            OutputWriter::Compressed(w) => w.flush(),
        }
    }
}

fn get_writer(path: &Path, compress: bool) -> Result<OutputWriter> {
    let file = File::create(path).with_context(|| format!("Failed to create output file {:?}", path))?;
    // Increased buffer size to 512KB
    let buf_writer = BufWriter::with_capacity(512 * 1024, file);

    if compress {
        // Use fast compression (level 1) for speed
        Ok(OutputWriter::Compressed(Box::new(GzEncoder::new(
            buf_writer,
            Compression::new(3),
        ))))
    } else {
        Ok(OutputWriter::Plain(buf_writer))
    }
}

// Helper to open a FASTQ reader that handles both .gz and plain text
fn get_reader(path: &Path) -> Result<fastq::Reader<BufReader<Box<dyn std::io::Read + Send>>>> {
    let file = File::open(path).with_context(|| format!("Failed to open input file {:?}", path))?;
    // Buffer the file input to improve decompression speed (512KB)
    let file_buf = BufReader::with_capacity(512 * 1024, file);
    
    let reader: Box<dyn std::io::Read + Send> = if path.extension().and_then(|s| s.to_str()) == Some("gz") {
        // flate2's MultiGzDecoder works best with a buffered reader
        Box::new(MultiGzDecoder::new(file_buf))
    } else {
        Box::new(file_buf)
    };
    
    // We pass the reader to fastq::Reader::new which wraps it in a BufReader.
    Ok(fastq::Reader::new(reader))
}

fn find_file_pair(candidates: &[PathBuf], barcode_id: &str) -> Result<(PathBuf, PathBuf)> {
    let mut r1_path: Option<PathBuf> = None;
    let mut r2_path: Option<PathBuf> = None;

    for path in candidates {
        let fname_str = path.file_name().unwrap().to_string_lossy();
        let fname_lower = fname_str.to_lowercase();

        // Strict extension check: Must be a FASTQ file
        if !fname_lower.ends_with(".fastq") 
            && !fname_lower.ends_with(".fq") 
            && !fname_lower.ends_with(".fastq.gz") 
            && !fname_lower.ends_with(".fq.gz") {
            continue;
        }
        
        // Barcode check
        if !fname_str.contains(&format!("_{}_", barcode_id)) && !fname_str.contains(&format!("-{}-", barcode_id)) {
            continue;
        }

        if fname_str.contains("_1.") || fname_str.contains("_R1.") || fname_str.contains(".R1.") {
            r1_path = Some(path.clone());
        } else if fname_str.contains("_2.") || fname_str.contains("_R2.") || fname_str.contains(".R2.") {
            r2_path = Some(path.clone());
        }
    }

    match (r1_path, r2_path) {
        (Some(r1), Some(r2)) => Ok((r1, r2)),
        _ => Err(anyhow::anyhow!("Could not find pair for barcode ID '{}' in directory.", barcode_id)),
    }
}

fn merge_files(output_path: &Path, parts: &[PathBuf]) -> Result<()> {
    let mut final_file = File::create(output_path)?;
    for part in parts {
        let mut part_file = File::open(part)?;
        std::io::copy(&mut part_file, &mut final_file)?;
    }
    Ok(())
}

fn main() -> Result<()> {
    let args = Args::parse();

    // Configure thread pool if user specified a limit
    if args.threads > 0 {
        rayon::ThreadPoolBuilder::new()
            .num_threads(args.threads)
            .build_global()
            .context("Failed to build thread pool")?;
    }

    // Create output directory if it doesn't exist
    if !args.output_dir.exists() {
        fs::create_dir_all(&args.output_dir).context("Failed to create output directory")?;
    }

    let ext = if args.compress { "fastq.gz" } else { "fastq" };
    let out_r1_path = args.output_dir.join(format!("demultiplexed.R1.{}", ext));
    let out_r2_path = args.output_dir.join(format!("demultiplexed.R2.{}", ext));

    // Check if output files already exist to prevent overwriting
    if out_r1_path.exists() || out_r2_path.exists() {
        return Err(anyhow::anyhow!("Output files already exist in {:?}. Please remove them or use a different output directory.", args.output_dir));
    }

    // Temporary files will be placed in the output directory
    let temp_dir = &args.output_dir;

    println!("Output R1: {:?}", out_r1_path);
    println!("Output R2: {:?}", out_r2_path);

    // Pre-scan directory
    let candidates: Vec<PathBuf> = WalkDir::new(&args.input_dir)
        .max_depth(1)
        .into_iter()
        .filter_map(|e| e.ok())
        .filter(|e| e.file_type().is_file())
        .map(|e| e.path().to_owned())
        .collect();

    let mut csv_rdr = csv::ReaderBuilder::new()
        .has_headers(false)
        .delimiter(if args.samplesheet.extension().map_or(false, |e| e == "tsv") { b'\t' } else { b',' })
        .from_path(&args.samplesheet)
        .context("Failed to open samplesheet")?;

    // Collect all records first to allow parallel processing
    let records: Vec<Result<csv::StringRecord, csv::Error>> = csv_rdr.records().collect();

    // Verify header and filter
    let valid_records: Vec<(usize, csv::StringRecord)> = records.into_iter().enumerate().filter_map(|(idx, res)| {
        match res {
            Ok(rec) => {
                 if rec.len() < 2 {
                    eprintln!("Warning: Skipping line {} (not enough columns)", idx + 1);
                    return None;
                }
                let barcode_id = rec[0].trim();
                let barcode_seq = rec[1].trim();
                // Skip header
                if idx == 0 && (barcode_id.to_lowercase().contains("id") || barcode_seq.to_lowercase().contains("seq")) {
                    println!("Detected header line, skipping: {:?}", rec);
                    return None;
                }
                Some((idx, rec))
            },
            Err(e) => {
                eprintln!("Error parsing CSV line {}: {}", idx + 1, e);
                None
            }
        }
    }).collect();

    println!("Processing {} samples in parallel...", valid_records.len());

    // Validation: Check for duplicate Barcode IDs and Sequences
    let mut seen_ids = std::collections::HashSet::new();
    let mut seen_seqs = std::collections::HashSet::new();
    let mut has_duplicates = false;

    for (_, rec) in &valid_records {
        let id = rec[0].trim();
        let seq = rec[1].trim();

        if !seen_ids.insert(id) {
            eprintln!("Error: Duplicate Barcode ID found in samplesheet: {}", id);
            has_duplicates = true;
        }
        if !seen_seqs.insert(seq) {
            eprintln!("Error: Duplicate Sequence found in samplesheet: {}", seq);
            has_duplicates = true;
        }
    }

    if has_duplicates {
        return Err(anyhow::anyhow!("Duplicate entries detected in samplesheet. Please fix them and retry."));
    }

    // Process in parallel
    // We return a list of created temporary files to merge later
    let chunk_files: Vec<Result<(PathBuf, PathBuf)>> = valid_records.par_iter().map(|(idx, record)| {
        let barcode_id = record[0].trim();
        let barcode_seq = record[1].trim();
        
        // Define temporary file paths unique to this job
        let temp_r1 = temp_dir.join(format!("tmp_{}_{}.R1.{}", idx, barcode_id, ext));
        let temp_r2 = temp_dir.join(format!("tmp_{}_{}.R2.{}", idx, barcode_id, ext));

        // Find inputs
        let (r1_file, r2_file) = match find_file_pair(&candidates, barcode_id) {
            Ok(pair) => pair,
            Err(e) => {
                eprintln!("Error: Barcode {} -> {}", barcode_id, e);
                return Err(e);
            }
        };
        
        println!("Barcode {} -> R1: {:?}, R2: {:?}", barcode_id, r1_file.file_name().unwrap_or_default(), r2_file.file_name().unwrap_or_default());

        let mut r1_reader = get_reader(&r1_file)?;
        let mut r2_reader = get_reader(&r2_file)?;
        
        let mut writer_r1 = fastq::Writer::new(get_writer(&temp_r1, args.compress)?);
        let mut writer_r2 = fastq::Writer::new(get_writer(&temp_r2, args.compress)?);

        let mut rec1 = fastq::Record::new();
        let mut rec2 = fastq::Record::new();
        
        // Local reusable buffers
        let mut new_seq = Vec::with_capacity(256);
        let mut new_qual = Vec::with_capacity(256);

        loop {
            let res1 = r1_reader.read(&mut rec1);
            let res2 = r2_reader.read(&mut rec2);

            match (res1, res2) {
                (Ok(()), Ok(())) => {
                    if rec1.is_empty() && rec2.is_empty() { break; }
                    if rec1.is_empty() || rec2.is_empty() {
                         let err = anyhow::anyhow!("Mismatch/Truncated ID {} (R1 empty: {}, R2 empty: {})", barcode_id, rec1.is_empty(), rec2.is_empty());
                         eprintln!("Error: {}", err);
                         return Err(err);
                    }
                    
                    let id1 = rec1.id();
                    let id2 = rec2.id();
                    
                    if id1 != id2 {
                        // Allow mismatch if it's just the /1 and /2 suffix (Casava < 1.8 style)
                        let is_suffix_mismatch = id1.ends_with("/1") 
                            && id2.ends_with("/2") 
                            && id1.len() > 2 
                            && id2.len() > 2 
                            && &id1[..id1.len()-2] == &id2[..id2.len()-2];

                        if !is_suffix_mismatch {
                            let err = anyhow::anyhow!("ID Mismatch for {}. R1: '{}', R2: '{}'", barcode_id, id1, id2);
                            eprintln!("Error: {}", err);
                            return Err(err);
                        }
                    }

                    writer_r1.write_record(&rec1)?;

                    new_seq.clear();
                    new_seq.extend_from_slice(rec2.seq());
                    new_seq.extend_from_slice(barcode_seq.as_bytes());

                    new_qual.clear();
                    new_qual.extend_from_slice(rec2.qual());
                    new_qual.extend(std::iter::repeat(b'I').take(barcode_seq.len()));

                    let new_rec2 = fastq::Record::with_attrs(rec2.id(), rec2.desc(), &new_seq, &new_qual);
                    writer_r2.write_record(&new_rec2)?;
                }
                (Err(e), _) | (_, Err(e)) => {
                    eprintln!("Error: Read error for barcode {}: {}", barcode_id, e);
                    return Err(anyhow::anyhow!("Read error: {}", e));
                }
            }
        }
        
        // Ensure writers are flushed/dropped
        drop(writer_r1);
        drop(writer_r2);

        println!("Successfully processed Barcode ID: {}", barcode_id);
        Ok((temp_r1, temp_r2))
    }).collect();

    // Check for errors and collect valid paths
    let mut temp_r1_paths = Vec::new();
    let mut temp_r2_paths = Vec::new();
    let mut failed_count = 0;

    for res in chunk_files {
        match res {
            Ok((p1, p2)) => {
                temp_r1_paths.push(p1);
                temp_r2_paths.push(p2);
            },
            Err(_) => {
                failed_count += 1;
            }
        }
    }

    println!("\nSummary:");
    println!("  - Succeeded: {}", temp_r1_paths.len());
    println!("  - Failed:    {}", failed_count);

    if failed_count > 0 {
        println!("\nCleaning up temporary files due to errors...");
        for p in temp_r1_paths.iter().chain(temp_r2_paths.iter()) {
            let _ = std::fs::remove_file(p);
        }
        return Err(anyhow::anyhow!("Processing failed for {} samples. Final merge cancelled.", failed_count));
    }

    if temp_r1_paths.is_empty() {
        println!("No files processed.");
        return Ok(());
    }

    println!("\nMerging {} temporary files...", temp_r1_paths.len());
    
    merge_files(&out_r1_path, &temp_r1_paths)?;
    merge_files(&out_r2_path, &temp_r2_paths)?;

    println!("Cleaning up temporary files...");
    for p in temp_r1_paths.iter().chain(temp_r2_paths.iter()) {
        let _ = std::fs::remove_file(p);
    }

    println!("Done. Merged files created.");

    Ok(())
}
