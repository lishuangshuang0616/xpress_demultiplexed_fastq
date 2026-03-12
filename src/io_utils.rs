use anyhow::{Context, Result};
use bio::io::fastq;
use flate2::read::MultiGzDecoder;
use rayon::prelude::*;
use std::fs::File;
use std::io::{BufReader, Read};
use std::os::unix::fs::FileExt;
use std::path::Path;

const IO_BUFFER_SIZE: usize = 1024 * 1024;

// Helper to open a FASTQ reader that handles both .gz and plain text
pub fn get_reader(path: &Path) -> Result<fastq::Reader<BufReader<Box<dyn std::io::Read + Send>>>> {
    let file = File::open(path).with_context(|| format!("Failed to open input file {:?}", path))?;
    let file_buf = BufReader::with_capacity(IO_BUFFER_SIZE, file);

    let reader: Box<dyn std::io::Read + Send> =
        if path.extension().and_then(|s| s.to_str()) == Some("gz") {
            Box::new(MultiGzDecoder::new(file_buf))
        } else {
            Box::new(file_buf)
        };

    Ok(fastq::Reader::new(reader))
}

pub fn merge_files(output_path: &Path, parts: &[std::path::PathBuf]) -> Result<()> {
    let out = File::create(output_path)
        .with_context(|| format!("Failed to create merged output {:?}", output_path))?;

    // 1. Assess file capacities and compute block offsets
    let mut total_len = 0u64;
    let mut offsets = Vec::with_capacity(parts.len());
    for part in parts {
        let size = std::fs::metadata(part)?.len();
        offsets.push(total_len);
        total_len += size;
    }

    // Predict block alignments and prevent APFS/Ext4 file fragmentation
    out.set_len(total_len)
        .context("Failed to pre-allocate output file length")?;

    // 2. Concurrently push distinct payload blocks into the unified target boundary
    parts.par_iter().zip(offsets.into_par_iter()).try_for_each(
        |(part, start_offset)| -> Result<()> {
            let mut input = File::open(part)
                .with_context(|| format!("Failed to open temp chunk {:?}", part))?;

            let mut buffer = vec![0u8; 4 * 1024 * 1024]; // 4MB Thread-local payload buffer
            let mut current_offset = start_offset;

            loop {
                let n = input.read(&mut buffer)?;
                if n == 0 {
                    break;
                }
                out.write_at(&buffer[..n], current_offset)?;
                current_offset += n as u64;
            }

            // Instantly recover filesystem footprint immediately upon discrete chunk ingestion completion
            let _ = std::fs::remove_file(part);

            Ok(())
        },
    )?;

    Ok(())
}
