use crate::io_utils::{get_reader, merge_files};
use crate::samplesheet::SampleSpec;
use anyhow::Result;
use bio::io::fastq;
use bio::io::fastq::FastqRead;
use crossbeam_channel::bounded;
use flate2::write::GzEncoder;
use flate2::Compression;
use indicatif::{ProgressBar, ProgressStyle};
use rayon::prelude::*;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};

const BATCH_FLUSH_THRESHOLD: usize = 1024 * 1024;
const PROGRESS_READ_STEP: u64 = 100_000;
const SAMPLE_CHUNK_CHANNEL_CAPACITY: usize = 4;

#[derive(Clone)]
pub struct SampleJob {
    pub order: usize,
    pub barcode_id: String,
    pub barcode_seq: String,
    pub r1_file: PathBuf,
    pub r2_file: PathBuf,
}

pub struct ProcessingStats {
    pub read_pairs: u64,
}

struct SampleChunk {
    r1: Vec<u8>,
    r2: Vec<u8>,
}

pub fn build_sample_jobs(
    samples: Vec<SampleSpec>,
    file_pair_index: &HashMap<String, (Option<PathBuf>, Option<PathBuf>)>,
) -> Result<Vec<SampleJob>> {
    let mut sample_jobs = Vec::with_capacity(samples.len());
    let mut unresolved = 0usize;

    for (order, sample) in samples.into_iter().enumerate() {
        match file_pair_index.get(&sample.barcode_id) {
            Some((Some(r1_file), Some(r2_file))) => {
                sample_jobs.push(SampleJob {
                    order,
                    barcode_id: sample.barcode_id,
                    barcode_seq: sample.barcode_seq,
                    r1_file: r1_file.clone(),
                    r2_file: r2_file.clone(),
                });
            }
            _ => {
                eprintln!(
                    "Warning: Barcode {} -> missing R1 or R2, skipping this sample.",
                    sample.barcode_id
                );
                unresolved += 1;
            }
        }
    }

    if unresolved > 0 {
        eprintln!(
            "Warning: {} sample(s) skipped due to missing R1/R2 pairs.",
            unresolved
        );
    }

    if sample_jobs.is_empty() {
        return Err(anyhow::anyhow!(
            "No valid sample pairs found. Nothing to process."
        ));
    }

    Ok(sample_jobs)
}

pub fn process_jobs_to_final_outputs(
    sample_jobs: &[SampleJob],
    out_r1_path: &Path,
    out_r2_path: &Path,
    gzip_level: u32,
    threads: usize,
) -> Result<ProcessingStats> {
    let total_samples = sample_jobs.len();
    let total_threads = effective_threads(threads);

    let progress = ProgressBar::new_spinner();
    progress.set_style(
        ProgressStyle::default_spinner()
            .template("{spinner:.green} [{elapsed_precise}] {pos} reads ({per_sec}) {msg}")
            .unwrap(),
    );
    progress.enable_steady_tick(std::time::Duration::from_millis(100));
    progress.set_message(format!(
        "workers {} threads, samples 0/{}",
        total_threads, total_samples
    ));

    let tmp_dir = out_r1_path
        .parent()
        .unwrap_or_else(|| Path::new("."))
        .join(format!(".fastq_merger_tmp_{}", std::process::id()));
    std::fs::create_dir_all(&tmp_dir)?;

    let process_fn = || {
        sample_jobs
            .par_iter()
            .map(|job| {
                process_one_sample_to_temp(
                    job,
                    gzip_level,
                    &tmp_dir,
                    &progress,
                    total_samples,
                    total_threads,
                )
            })
            .collect::<Vec<Result<(usize, PathBuf, PathBuf, u64)>>>()
    };

    let results = if threads > 0 {
        rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()?
            .install(process_fn)
    } else {
        process_fn()
    };

    let mut ok_results = Vec::with_capacity(results.len());
    let mut failed = 0usize;
    let mut total_read_pairs = 0u64;
    for res in results {
        match res {
            Ok((order, r1, r2, reads)) => {
                ok_results.push((order, r1, r2));
                total_read_pairs += reads;
            }
            Err(e) => {
                failed += 1;
                eprintln!("Error: {}", e);
            }
        }
    }

    if failed > 0 {
        cleanup_tmp_dir(&tmp_dir);
        progress.finish_with_message(format!("Failed: {}/{} samples", failed, total_samples));
        return Err(anyhow::anyhow!("Processing failed for {} samples.", failed));
    }

    ok_results.sort_by_key(|(order, _, _)| *order);
    let r1_parts: Vec<PathBuf> = ok_results.iter().map(|(_, r1, _)| r1.clone()).collect();
    let r2_parts: Vec<PathBuf> = ok_results.iter().map(|(_, _, r2)| r2.clone()).collect();

    progress.set_message("Merging outputs (R1/R2)...".to_string());

    std::thread::scope(|scope| -> Result<()> {
        let r1_parts_ref = &r1_parts;
        let r2_parts_ref = &r2_parts;
        let t1 = scope.spawn(move || merge_files(out_r1_path, r1_parts_ref));
        let t2 = scope.spawn(move || merge_files(out_r2_path, r2_parts_ref));
        t1.join()
            .map_err(|_| anyhow::anyhow!("R1 merge thread panicked"))??;
        t2.join()
            .map_err(|_| anyhow::anyhow!("R2 merge thread panicked"))??;
        Ok(())
    })?;

    cleanup_tmp_dir(&tmp_dir);

    println!("\nSummary:");
    println!("  - Succeeded: {}", total_samples);
    println!("  - Failed:    0");

    progress.set_position(total_read_pairs);
    progress.finish_with_message(format!(
        "Done! processed {} reads, samples {}/{}",
        total_read_pairs, total_samples, total_samples
    ));

    Ok(ProcessingStats {
        read_pairs: total_read_pairs,
    })
}

fn process_one_sample_to_temp(
    job: &SampleJob,
    gzip_level: u32,
    tmp_dir: &Path,
    progress: &ProgressBar,
    total_samples: usize,
    total_threads: usize,
) -> Result<(usize, PathBuf, PathBuf, u64)> {
    let barcode_id = job.barcode_id.as_str();
    let barcode_seq = job.barcode_seq.as_bytes();
    let barcode_qual_tail = vec![b'I'; barcode_seq.len()];

    let temp_r1 = tmp_dir.join(format!("tmp_{:08}_{}.R1.fastq.gz", job.order, barcode_id));
    let temp_r2 = tmp_dir.join(format!("tmp_{:08}_{}.R2.fastq.gz", job.order, barcode_id));

    let mut r1_reader = get_reader(&job.r1_file)?;
    let mut r2_reader = get_reader(&job.r2_file)?;
    let mut rec1 = fastq::Record::new();
    let mut rec2 = fastq::Record::new();

    let mut batch_r1 = Vec::with_capacity(BATCH_FLUSH_THRESHOLD + 4096);
    let mut batch_r2 = Vec::with_capacity(BATCH_FLUSH_THRESHOLD + 4096);
    let mut read_pairs = 0u64;
    let mut progress_chunk = 0u64;

    let (chunk_tx, chunk_rx) = bounded::<SampleChunk>(SAMPLE_CHUNK_CHANNEL_CAPACITY);

    std::thread::scope(|scope| -> Result<()> {
        let t_r1 = temp_r1.clone();
        let t_r2 = temp_r2.clone();
        let writer = scope.spawn(move || -> Result<()> {
            let f1 = File::create(&t_r1)?;
            let f2 = File::create(&t_r2)?;
            let mut w1 = GzEncoder::new(
                BufWriter::with_capacity(8 * 1024 * 1024, f1),
                Compression::new(gzip_level.min(9)),
            );
            let mut w2 = GzEncoder::new(
                BufWriter::with_capacity(8 * 1024 * 1024, f2),
                Compression::new(gzip_level.min(9)),
            );
            while let Ok(chunk) = chunk_rx.recv() {
                w1.write_all(&chunk.r1)?;
                w2.write_all(&chunk.r2)?;
            }
            let mut buf_w1 = w1.finish()?;
            let mut buf_w2 = w2.finish()?;
            buf_w1.flush()?;
            buf_w2.flush()?;
            Ok(())
        });

        loop {
            let res1 = r1_reader.read(&mut rec1);
            let res2 = r2_reader.read(&mut rec2);
            match (res1, res2) {
                (Ok(()), Ok(())) => {
                    if rec1.is_empty() && rec2.is_empty() {
                        break;
                    }
                    if rec1.is_empty() || rec2.is_empty() {
                        return Err(anyhow::anyhow!(
                            "Barcode {} -> Mismatch/Truncated (R1 empty: {}, R2 empty: {})",
                            barcode_id,
                            rec1.is_empty(),
                            rec2.is_empty()
                        ));
                    }

                    let id1 = rec1.id();
                    let id2 = rec2.id();
                    if id1 != id2 {
                        let is_suffix_mismatch = id1.ends_with("/1")
                            && id2.ends_with("/2")
                            && id1.len() > 2
                            && id2.len() > 2
                            && id1[..id1.len() - 2] == id2[..id2.len() - 2];
                        if !is_suffix_mismatch {
                            return Err(anyhow::anyhow!(
                                "Barcode {} -> ID Mismatch R1='{}' R2='{}'",
                                barcode_id,
                                id1,
                                id2
                            ));
                        }
                    }

                    append_fastq_record(
                        &mut batch_r1,
                        rec1.id(),
                        rec1.desc(),
                        rec1.seq(),
                        rec1.qual(),
                    );

                    append_fastq_header(&mut batch_r2, rec2.id(), rec2.desc());
                    batch_r2.extend_from_slice(rec2.seq());
                    batch_r2.extend_from_slice(barcode_seq);
                    batch_r2.push(b'\n');
                    batch_r2.extend_from_slice(b"+\n");
                    batch_r2.extend_from_slice(rec2.qual());
                    batch_r2.extend_from_slice(&barcode_qual_tail);
                    batch_r2.push(b'\n');

                    read_pairs += 1;
                    progress_chunk += 1;
                    if progress_chunk >= PROGRESS_READ_STEP {
                        progress.inc(progress_chunk);
                        progress_chunk = 0;
                    }

                    if batch_r1.len() >= BATCH_FLUSH_THRESHOLD {
                        chunk_tx.send(SampleChunk {
                            r1: std::mem::replace(
                                &mut batch_r1,
                                Vec::with_capacity(BATCH_FLUSH_THRESHOLD + 4096),
                            ),
                            r2: std::mem::replace(
                                &mut batch_r2,
                                Vec::with_capacity(BATCH_FLUSH_THRESHOLD + 4096),
                            ),
                        })?;
                    }
                }
                (Err(e), _) | (_, Err(e)) => {
                    return Err(anyhow::anyhow!(
                        "Barcode {} -> Read error: {}",
                        barcode_id,
                        e
                    ));
                }
            }
        }

        if !batch_r1.is_empty() {
            chunk_tx.send(SampleChunk {
                r1: batch_r1,
                r2: batch_r2,
            })?;
        }
        drop(chunk_tx);

        match writer.join() {
            Ok(r) => r?,
            Err(_) => return Err(anyhow::anyhow!("Writer thread panicked for {}", barcode_id)),
        }

        Ok(())
    })?;

    if progress_chunk > 0 {
        progress.inc(progress_chunk);
    }
    progress.set_message(format!(
        "workers {} threads, samples processing.../{}, latest={}",
        total_threads, total_samples, barcode_id
    ));

    Ok((job.order, temp_r1, temp_r2, read_pairs))
}

fn append_fastq_record(buf: &mut Vec<u8>, id: &str, desc: Option<&str>, seq: &[u8], qual: &[u8]) {
    append_fastq_header(buf, id, desc);
    buf.extend_from_slice(seq);
    buf.push(b'\n');
    buf.extend_from_slice(b"+\n");
    buf.extend_from_slice(qual);
    buf.push(b'\n');
}

fn append_fastq_header(buf: &mut Vec<u8>, id: &str, desc: Option<&str>) {
    buf.push(b'@');
    buf.extend_from_slice(id.as_bytes());
    if let Some(d) = desc {
        buf.push(b' ');
        buf.extend_from_slice(d.as_bytes());
    }
    buf.push(b'\n');
}

fn effective_threads(threads: usize) -> usize {
    if threads > 0 {
        return threads;
    }
    std::thread::available_parallelism()
        .map(|n| n.get())
        .unwrap_or(8)
}

fn cleanup_tmp_dir(tmp_dir: &Path) {
    let _ = std::fs::remove_dir_all(tmp_dir);
}
