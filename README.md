# fastq_merger

`fastq_merger` merges paired-end FASTQ samples into one final R1/R2 output pair based on a samplesheet.

Processing flow:
- Read samplesheet (`Barcode_ID`, `Sequence`)
- Match input R1/R2 files by `Barcode_ID`
- Append the barcode sequence to the end of each R2 read
- Append `I` quality characters with the same barcode length
- Process samples in parallel, then merge in samplesheet order into:
  - `demultiplexed.R1.fastq.gz`
  - `demultiplexed.R2.fastq.gz`

## 1. Commands

```bash
fastq_merger <COMMAND>
```

Available subcommands:
- `run`: execute processing
- `benchmark`: execute processing and print stage timing

Legacy compatibility: if no subcommand is provided, the tool defaults to `run`.

## 2. CLI Options

`run` and `benchmark` share the same options:

- `-i, --input-dir <INPUT_DIR>`: input directory (non-recursive scan)
- `-s, --samplesheet <SAMPLESHEET>`: samplesheet path (CSV/TSV)
- `--gzip-level <GZIP_LEVEL>`: gzip level `0-9`, default `3`
- `-o, --output-dir <OUTPUT_DIR>`: output directory, default `output`
- `-t, --threads <THREADS>`: thread count, default `8`

Notes:
- This version always writes gzip output (`.fastq.gz`).
- If output files already exist, the run aborts to prevent overwrite.

## 3. Samplesheet Format

At least 2 columns:
- Column 1: `Barcode_ID` (used for filename matching)
- Column 2: `Sequence` (appended to R2)

Example:

```csv
Barcode_ID,Sequence
1,ATGCATGC
2,CGTACGTA
```

Rules:
- First line is auto-skipped as header if it contains `id`/`seq` (case-insensitive)
- Duplicate `Barcode_ID` or `Sequence` is rejected

## 4. File Matching Rules

Supported extensions:
- `.fastq` `.fq` `.fastq.gz` `.fq.gz`

`Barcode_ID` matching patterns:
- `_ID_` (for example `Sample_1_L01_R1_001.fastq.gz`)
- `-ID-` (for example `Lib-1-Rep1_R1.fq.gz`)
- Prefix `ID<delimiter>` (for example `1_R1_001.fastq.gz`)

R1/R2 detection patterns:
- R1: `_R1.` `_R1_` `.R1.` `.R1_` `_1.` `_1_`
- R2: `_R2.` `_R2_` `.R2.` `.R2_` `_2.` `_2_`

Also:
- Symlinks in `input-dir` are included if they point to files.

## 5. Examples

### Standard run

```bash
fastq_merger run \
  --input-dir ./raw_data \
  --samplesheet ./sample.sheet \
  --output-dir ./out \
  --threads 8
```

### Run with timing metrics

```bash
fastq_merger benchmark \
  --input-dir ./raw_data \
  --samplesheet ./sample.sheet \
  --output-dir ./out \
  --threads 8
```

## 6. Output and Progress

Before processing, the tool prints:
- Output file paths
- Resolved sample mapping (`sample`, `barcode`, `R1`, `R2`)

Progress bar shows:
- Processed read count
- Current throughput
- Sample processing status
- Merge-stage message: `Merging outputs (R1/R2)...`

## 7. Common Errors

- `Duplicate Barcode ID found` / `Duplicate Sequence found`
  - Duplicates in samplesheet
- `Could not find pair for barcode ID 'xxx'`
  - R1/R2 files not found or filename does not match rules
- `ID Mismatch for xxx`
  - R1/R2 read IDs differ (except `/1` vs `/2` suffix compatibility)
- `Mismatch/Truncated`
  - R1/R2 length mismatch or truncated/corrupted input
- `Output files already exist`
  - Target output files already exist in output directory

## 8. Performance Tips

- Prefer `--gzip-level 1` for speed
- Start from `--threads 8`, then tune per machine/storage
- On network storage (NFS), too many threads can reduce throughput
