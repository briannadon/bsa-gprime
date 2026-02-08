use crate::types::{GStatisticResult, SignificanceResult};

/// Metadata for a single chromosome panel.
pub struct ChromPanel {
    pub chrom: String,
    /// Indices into the original results slice for this chromosome, sorted by position.
    pub indices: Vec<usize>,
}

/// Natural chromosome sort: chr1 < chr2 < ... < chr10 < chrX
fn natural_chrom_cmp(a: &str, b: &str) -> std::cmp::Ordering {
    let a_num = extract_chrom_number(a);
    let b_num = extract_chrom_number(b);
    match (a_num, b_num) {
        (Some(an), Some(bn)) => an.cmp(&bn),
        (Some(_), None) => std::cmp::Ordering::Less,
        (None, Some(_)) => std::cmp::Ordering::Greater,
        (None, None) => a.cmp(b),
    }
}

fn extract_chrom_number(chrom: &str) -> Option<u64> {
    let stripped = chrom
        .strip_prefix("chr")
        .or_else(|| chrom.strip_prefix("Chr"))
        .or_else(|| chrom.strip_prefix("CHR"))
        .unwrap_or(chrom);
    stripped.parse::<u64>().ok()
}

/// Group significance results by chromosome, returning at most `max_chroms` panels
/// in natural sort order.
pub fn prepare_chrom_panels_sig(
    results: &[SignificanceResult],
    max_chroms: usize,
) -> Vec<ChromPanel> {
    let mut chrom_map: std::collections::HashMap<String, Vec<usize>> =
        std::collections::HashMap::new();

    for (i, r) in results.iter().enumerate() {
        chrom_map
            .entry(r.raw.variant.chrom.clone())
            .or_default()
            .push(i);
    }

    let mut chroms: Vec<String> = chrom_map.keys().cloned().collect();
    chroms.sort_by(|a, b| natural_chrom_cmp(a, b));
    chroms.truncate(max_chroms);

    chroms
        .into_iter()
        .map(|chrom| {
            let mut indices = chrom_map.remove(&chrom).unwrap();
            indices.sort_by_key(|&i| results[i].raw.variant.pos);
            ChromPanel { chrom, indices }
        })
        .collect()
}

/// Group basic results by chromosome, returning at most `max_chroms` panels
/// in natural sort order.
pub fn prepare_chrom_panels_basic(
    results: &[GStatisticResult],
    max_chroms: usize,
) -> Vec<ChromPanel> {
    let mut chrom_map: std::collections::HashMap<String, Vec<usize>> =
        std::collections::HashMap::new();

    for (i, r) in results.iter().enumerate() {
        chrom_map
            .entry(r.variant.chrom.clone())
            .or_default()
            .push(i);
    }

    let mut chroms: Vec<String> = chrom_map.keys().cloned().collect();
    chroms.sort_by(|a, b| natural_chrom_cmp(a, b));
    chroms.truncate(max_chroms);

    chroms
        .into_iter()
        .map(|chrom| {
            let mut indices = chrom_map.remove(&chrom).unwrap();
            indices.sort_by_key(|&i| results[i].variant.pos);
            ChromPanel { chrom, indices }
        })
        .collect()
}
