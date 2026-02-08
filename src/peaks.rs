use crate::types::{GStatisticResult, SignificanceResult};
use anyhow::Result;
use csv::Writer;
use std::path::Path;

/// A detected QTL peak region.
#[derive(Debug, Clone)]
pub struct QtlPeak {
    pub chrom: String,
    pub start: u64,
    pub end: u64,
    pub width_bp: u64,
    pub n_snps: u64,
    pub max_neg_log10_p: Option<f64>,
    pub max_g_prime: Option<f64>,
    pub mean_delta_snp_index: f64,
}

/// Identify QTL peaks from significance results.
///
/// Groups contiguous runs of significant variants on each chromosome,
/// merging clusters within `merge_distance` bp.
pub fn identify_peaks(results: &[SignificanceResult], merge_distance: u64) -> Vec<QtlPeak> {
    // Group by chromosome, preserving sort order
    let mut chroms: Vec<String> = Vec::new();
    for r in results {
        let c = &r.raw.variant.chrom;
        if chroms.last().map_or(true, |last| last != c) {
            if !chroms.contains(c) {
                chroms.push(c.clone());
            }
        }
    }
    chroms.sort_by(|a, b| natural_chrom_cmp(a, b));

    let mut peaks = Vec::new();

    for chrom in &chroms {
        let mut chrom_sig: Vec<&SignificanceResult> = results
            .iter()
            .filter(|r| r.raw.variant.chrom == *chrom && r.significant)
            .collect();
        chrom_sig.sort_by_key(|r| r.raw.variant.pos);

        if chrom_sig.is_empty() {
            continue;
        }

        // Group into clusters by merge_distance
        let mut clusters: Vec<Vec<&SignificanceResult>> = Vec::new();
        let mut current_cluster: Vec<&SignificanceResult> = vec![chrom_sig[0]];

        for r in &chrom_sig[1..] {
            let last_pos = current_cluster.last().unwrap().raw.variant.pos;
            if r.raw.variant.pos <= last_pos + merge_distance {
                current_cluster.push(r);
            } else {
                clusters.push(std::mem::take(&mut current_cluster));
                current_cluster.push(r);
            }
        }
        clusters.push(current_cluster);

        for cluster in clusters {
            let start = cluster.first().unwrap().raw.variant.pos;
            let end = cluster.last().unwrap().raw.variant.pos;
            let n_snps = cluster.len() as u64;

            let min_p = cluster.iter().map(|r| r.p_value).fold(f64::INFINITY, f64::min);
            let max_neg_log10_p = -min_p.max(1e-300).log10();

            let max_gp = cluster
                .iter()
                .map(|r| r.g_prime)
                .fold(f64::NEG_INFINITY, f64::max);

            let mean_delta: f64 = cluster.iter().map(|r| r.raw.delta_snp_index).sum::<f64>()
                / cluster.len() as f64;

            peaks.push(QtlPeak {
                chrom: chrom.clone(),
                start,
                end,
                width_bp: end - start,
                n_snps,
                max_neg_log10_p: Some(max_neg_log10_p),
                max_g_prime: Some(max_gp),
                mean_delta_snp_index: mean_delta,
            });
        }
    }

    // Sort by max_neg_log10_p descending
    peaks.sort_by(|a, b| {
        b.max_neg_log10_p
            .unwrap_or(0.0)
            .partial_cmp(&a.max_neg_log10_p.unwrap_or(0.0))
            .unwrap()
    });

    peaks
}

/// Identify QTL peaks by percentile threshold on G-statistic or G'.
///
/// Fallback for when significance testing is not available.
pub fn identify_peaks_by_percentile(
    results: &[GStatisticResult],
    percentile: f64,
    merge_distance: u64,
) -> Vec<QtlPeak> {
    if results.is_empty() {
        return Vec::new();
    }

    // Compute percentile threshold
    let mut g_values: Vec<f64> = results.iter().map(|r| r.g_statistic).collect();
    g_values.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let idx = ((percentile / 100.0) * g_values.len() as f64) as usize;
    let threshold = g_values[idx.min(g_values.len() - 1)];

    // Group by chromosome
    let mut chroms: Vec<String> = Vec::new();
    for r in results {
        let c = &r.variant.chrom;
        if !chroms.contains(c) {
            chroms.push(c.clone());
        }
    }
    chroms.sort_by(|a, b| natural_chrom_cmp(a, b));

    let mut peaks = Vec::new();

    for chrom in &chroms {
        let mut chrom_above: Vec<&GStatisticResult> = results
            .iter()
            .filter(|r| r.variant.chrom == *chrom && r.g_statistic > threshold)
            .collect();
        chrom_above.sort_by_key(|r| r.variant.pos);

        if chrom_above.is_empty() {
            continue;
        }

        let mut clusters: Vec<Vec<&GStatisticResult>> = Vec::new();
        let mut current_cluster: Vec<&GStatisticResult> = vec![chrom_above[0]];

        for r in &chrom_above[1..] {
            let last_pos = current_cluster.last().unwrap().variant.pos;
            if r.variant.pos <= last_pos + merge_distance {
                current_cluster.push(r);
            } else {
                clusters.push(std::mem::take(&mut current_cluster));
                current_cluster.push(r);
            }
        }
        clusters.push(current_cluster);

        for cluster in clusters {
            let start = cluster.first().unwrap().variant.pos;
            let end = cluster.last().unwrap().variant.pos;
            let n_snps = cluster.len() as u64;

            let max_g = cluster
                .iter()
                .map(|r| r.g_statistic)
                .fold(f64::NEG_INFINITY, f64::max);

            let mean_delta: f64 =
                cluster.iter().map(|r| r.delta_snp_index).sum::<f64>() / cluster.len() as f64;

            peaks.push(QtlPeak {
                chrom: chrom.clone(),
                start,
                end,
                width_bp: end - start,
                n_snps,
                max_neg_log10_p: None,
                max_g_prime: Some(max_g),
                mean_delta_snp_index: mean_delta,
            });
        }
    }

    // Sort by max_g_prime descending
    peaks.sort_by(|a, b| {
        b.max_g_prime
            .unwrap_or(0.0)
            .partial_cmp(&a.max_g_prime.unwrap_or(0.0))
            .unwrap()
    });

    peaks
}

/// Write QTL peaks to a CSV file.
pub fn write_peaks_csv(peaks: &[QtlPeak], path: &Path) -> Result<()> {
    let mut wtr = Writer::from_path(path)?;

    wtr.write_record([
        "chrom",
        "start",
        "end",
        "width_bp",
        "n_snps",
        "max_neg_log10_p",
        "max_g_prime",
        "mean_delta_snp_index",
    ])?;

    for peak in peaks {
        wtr.write_record([
            &peak.chrom,
            &peak.start.to_string(),
            &peak.end.to_string(),
            &peak.width_bp.to_string(),
            &peak.n_snps.to_string(),
            &peak
                .max_neg_log10_p
                .map_or("NA".to_string(), |v| format!("{:.4}", v)),
            &peak
                .max_g_prime
                .map_or("NA".to_string(), |v| format!("{:.4}", v)),
            &format!("{:.6}", peak.mean_delta_snp_index),
        ])?;
    }

    wtr.flush()?;
    Ok(())
}

/// Print a formatted summary table of QTL peaks to stderr.
pub fn report_peaks(peaks: &[QtlPeak]) {
    eprintln!();
    eprintln!("{}", "=".repeat(90));

    if peaks.iter().any(|p| p.max_neg_log10_p.is_some()) {
        eprintln!("QTL Peak Detection (FDR-corrected significance)");
    } else {
        eprintln!("QTL Peak Detection (percentile-based threshold)");
    }

    eprintln!("{}", "=".repeat(90));
    eprintln!();

    if peaks.is_empty() {
        eprintln!("No QTL peaks found.");
        eprintln!();
        eprintln!("{}", "=".repeat(90));
        return;
    }

    eprintln!("Found {} QTL peak(s):", peaks.len());
    eprintln!();

    // Header
    eprintln!(
        "{:<12} {:>12} {:>12} {:>12} {:>8} {:>16} {:>12} {:>12}",
        "chrom", "start", "end", "width_bp", "n_snps", "-log10(p)_max", "g_prime_max", "mean_dSNP"
    );
    eprintln!("{}", "-".repeat(90));

    for peak in peaks {
        eprintln!(
            "{:<12} {:>12} {:>12} {:>12} {:>8} {:>16} {:>12} {:>12}",
            peak.chrom,
            peak.start,
            peak.end,
            peak.width_bp,
            peak.n_snps,
            peak.max_neg_log10_p
                .map_or("NA".to_string(), |v| format!("{:.2}", v)),
            peak.max_g_prime
                .map_or("NA".to_string(), |v| format!("{:.2}", v)),
            format!("{:.4}", peak.mean_delta_snp_index),
        );
    }

    eprintln!();
    eprintln!("{}", "=".repeat(90));
}

/// Natural chromosome sort: chr1 < chr2 < ... < chr10 < chr11 < chrX
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::types::Variant;

    fn make_sig_result(chrom: &str, pos: u64, significant: bool, p_value: f64) -> SignificanceResult {
        SignificanceResult {
            raw: GStatisticResult {
                variant: Variant {
                    chrom: chrom.to_string(),
                    pos,
                    ref_allele: "A".to_string(),
                    alt_allele: "T".to_string(),
                    resistant_ref_depth: 20,
                    resistant_alt_depth: 10,
                    resistant_dp: 30,
                    resistant_gq: 99,
                    susceptible_ref_depth: 25,
                    susceptible_alt_depth: 5,
                    susceptible_dp: 30,
                    susceptible_gq: 99,
                },
                g_statistic: 5.0,
                snp_index_resistant: 0.33,
                snp_index_susceptible: 0.17,
                delta_snp_index: 0.17,
            },
            g_prime: 4.5,
            p_value,
            q_value: p_value * 2.0,
            significant,
        }
    }

    #[test]
    fn test_identify_peaks_basic() {
        let results = vec![
            make_sig_result("chr1", 1000, true, 1e-10),
            make_sig_result("chr1", 2000, true, 1e-12),
            make_sig_result("chr1", 3000, true, 1e-8),
            make_sig_result("chr1", 100000, false, 0.5),
            make_sig_result("chr2", 5000, true, 1e-6),
        ];

        let peaks = identify_peaks(&results, 500_000);
        assert_eq!(peaks.len(), 2);
        // Should be sorted by -log10(p) descending
        assert!(peaks[0].max_neg_log10_p.unwrap() >= peaks[1].max_neg_log10_p.unwrap());
    }

    #[test]
    fn test_identify_peaks_merging() {
        let results = vec![
            make_sig_result("chr1", 1000, true, 1e-10),
            make_sig_result("chr1", 600_000, true, 1e-8),
        ];

        // With large merge distance, should be one peak
        let peaks = identify_peaks(&results, 1_000_000);
        assert_eq!(peaks.len(), 1);

        // With small merge distance, should be two peaks
        let peaks = identify_peaks(&results, 100_000);
        assert_eq!(peaks.len(), 2);
    }

    #[test]
    fn test_identify_peaks_empty() {
        let peaks = identify_peaks(&[], 500_000);
        assert!(peaks.is_empty());
    }

    #[test]
    fn test_natural_chrom_cmp() {
        assert_eq!(natural_chrom_cmp("chr1", "chr2"), std::cmp::Ordering::Less);
        assert_eq!(natural_chrom_cmp("chr2", "chr10"), std::cmp::Ordering::Less);
        assert_eq!(natural_chrom_cmp("chr10", "chr9"), std::cmp::Ordering::Greater);
        assert_eq!(natural_chrom_cmp("chrX", "chr1"), std::cmp::Ordering::Greater);
    }
}
