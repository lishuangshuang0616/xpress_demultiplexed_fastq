use std::collections::{HashMap, HashSet};
use std::path::{Path, PathBuf};
use walkdir::WalkDir;

pub fn collect_candidates(input_dir: &Path) -> Vec<PathBuf> {
    WalkDir::new(input_dir)
        .max_depth(1)
        .into_iter()
        .filter_map(|e| e.ok())
        .filter(|e| {
            let ft = e.file_type();
            if ft.is_file() {
                return true;
            }
            if ft.is_symlink() {
                return std::fs::metadata(e.path())
                    .map(|m| m.is_file())
                    .unwrap_or(false);
            }
            false
        })
        .map(|e| e.path().to_owned())
        .collect()
}

pub fn build_file_pair_index(
    candidates: &[PathBuf],
    barcode_ids: &HashSet<String>,
) -> HashMap<String, (Option<PathBuf>, Option<PathBuf>)> {
    let mut index: HashMap<String, (Option<PathBuf>, Option<PathBuf>)> = HashMap::new();

    for path in candidates {
        let fname = path.file_name().unwrap().to_string_lossy();
        if !is_fastq_filename(&fname) {
            continue;
        }

        let is_r1 = is_r1_filename(&fname);
        let is_r2 = is_r2_filename(&fname);
        if !is_r1 && !is_r2 {
            continue;
        }

        // Prefix pattern: <id><delimiter>..., e.g. C18A1_R1.fq.gz
        if let Some(prefix_id) = extract_prefix_id_candidate(&fname) {
            if barcode_ids.contains(prefix_id) {
                update_index_pair(&mut index, prefix_id, path, is_r1, is_r2);
            }
        }

        // Wrapped pattern: _<id>_ and -<id>-
        match_wrapped_ids(&fname, '_', barcode_ids, |barcode_id| {
            update_index_pair(&mut index, barcode_id, path, is_r1, is_r2);
        });

        match_wrapped_ids(&fname, '-', barcode_ids, |barcode_id| {
            update_index_pair(&mut index, barcode_id, path, is_r1, is_r2);
        });
    }

    index
}

fn is_fastq_filename(fname: &str) -> bool {
    let fname_lower = fname.to_lowercase();
    fname_lower.ends_with(".fastq")
        || fname_lower.ends_with(".fq")
        || fname_lower.ends_with(".fastq.gz")
        || fname_lower.ends_with(".fq.gz")
}

fn extract_prefix_id_candidate(fname: &str) -> Option<&str> {
    fname
        .find(['_', '-', '.'])
        .filter(|pos| *pos > 0)
        .map(|pos| &fname[..pos])
}

fn match_wrapped_ids<F>(
    fname: &str,
    delimiter: char,
    barcode_ids: &HashSet<String>,
    mut on_match: F,
) where
    F: FnMut(&str),
{
    // For pattern _<id>_ (or -<id>-), valid IDs are only middle segments.
    let mut parts = fname.split(delimiter);
    let _first = parts.next();
    let mut prev = parts.next();

    for next in parts {
        if let Some(candidate) = prev {
            if !candidate.is_empty() && barcode_ids.contains(candidate) {
                on_match(candidate);
            }
        }
        prev = Some(next);
    }
}

fn update_index_pair(
    index: &mut HashMap<String, (Option<PathBuf>, Option<PathBuf>)>,
    barcode_id: &str,
    path: &Path,
    is_r1: bool,
    is_r2: bool,
) {
    let entry = index.entry(barcode_id.to_string()).or_insert((None, None));

    if is_r1 {
        entry.0 = Some(path.to_path_buf());
    }
    if is_r2 {
        entry.1 = Some(path.to_path_buf());
    }
}

fn is_r1_filename(fname: &str) -> bool {
    fname.contains("_R1.")
        || fname.contains("_R1_")
        || fname.contains(".R1.")
        || fname.contains(".R1_")
        || fname.contains("_1.")
        || fname.contains("_1_")
}

fn is_r2_filename(fname: &str) -> bool {
    fname.contains("_R2.")
        || fname.contains("_R2_")
        || fname.contains(".R2.")
        || fname.contains(".R2_")
        || fname.contains("_2.")
        || fname.contains("_2_")
}
