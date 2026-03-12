use anyhow::{Context, Result};
use std::collections::HashSet;
use std::path::Path;

#[derive(Clone, Debug)]
pub struct SampleSpec {
    pub barcode_id: String,
    pub barcode_seq: String,
}

pub struct ParsedSamplesheet {
    pub samples: Vec<SampleSpec>,
    pub barcode_ids: HashSet<String>,
}

pub fn read_samplesheet(samplesheet_path: &Path) -> Result<ParsedSamplesheet> {
    let delimiter = if samplesheet_path.extension().map_or(false, |e| e == "tsv") {
        b'\t'
    } else {
        b','
    };

    let mut csv_rdr = csv::ReaderBuilder::new()
        .has_headers(false)
        .delimiter(delimiter)
        .from_path(samplesheet_path)
        .context("Failed to open samplesheet")?;

    let mut samples = Vec::new();
    let mut seen_ids: HashSet<String> = HashSet::new();
    let mut seen_seqs: HashSet<String> = HashSet::new();
    let mut has_duplicates = false;

    for (idx, rec_res) in csv_rdr.records().enumerate() {
        let rec = match rec_res {
            Ok(rec) => rec,
            Err(e) => {
                eprintln!("Error parsing CSV line {}: {}", idx + 1, e);
                continue;
            }
        };

        if rec.len() < 2 {
            eprintln!("Warning: Skipping line {} (not enough columns)", idx + 1);
            continue;
        }

        let barcode_id = rec[0].trim();
        let barcode_seq = rec[1].trim();

        // Skip header
        if idx == 0
            && (barcode_id.to_lowercase().contains("id")
                || barcode_seq.to_lowercase().contains("seq"))
        {
            println!("Detected header line, skipping: {:?}", rec);
            continue;
        }

        let barcode_id = barcode_id.to_string();
        let barcode_seq = barcode_seq.to_string();

        if !seen_ids.insert(barcode_id.clone()) {
            eprintln!(
                "Error: Duplicate Barcode ID found in samplesheet: {}",
                barcode_id
            );
            has_duplicates = true;
        }

        if !seen_seqs.insert(barcode_seq.clone()) {
            eprintln!(
                "Error: Duplicate Sequence found in samplesheet: {}",
                barcode_seq
            );
            has_duplicates = true;
        }

        samples.push(SampleSpec {
            barcode_id,
            barcode_seq,
        });
    }

    if has_duplicates {
        return Err(anyhow::anyhow!(
            "Duplicate entries detected in samplesheet. Please fix them and retry."
        ));
    }

    Ok(ParsedSamplesheet {
        samples,
        barcode_ids: seen_ids,
    })
}
