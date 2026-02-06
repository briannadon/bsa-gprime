use crate::types::GStatisticResult;
use rayon::prelude::*;

/// Tricube kernel: w(u) = (1 - |u|^3)^3 for |u| <= 1, else 0
fn tricube(u: f64) -> f64 {
    let abs_u = u.abs();
    if abs_u >= 1.0 {
        0.0
    } else {
        let t = 1.0 - abs_u.powi(3);
        t.powi(3)
    }
}

/// Compute Σkⱼ² (sum of squared normalized kernel weights) for variance estimation.
///
/// This is needed for equation 12 from Magwene et al. (2011):
/// Var[G'] ≈ Var[G] × Σkⱼ²
///
/// The tricube kernel weights are: kⱼ = (1 - Dⱼ³)³ / S_W
/// where S_W = Σ(1 - Dⱼ³)³ is the normalization factor
///
/// Returns the average Σkⱼ² across all SNPs, which can be used to scale
/// the parametric null distribution for smoothed G' values.
pub fn compute_smoothing_factor(results: &[GStatisticResult], bandwidth: u64) -> f64 {
    let n = results.len();
    if n == 0 {
        return 1.0;
    }

    // Collect chromosome boundary ranges
    let mut ranges: Vec<(usize, usize)> = Vec::new();
    let mut chrom_start = 0;
    while chrom_start < n {
        let chrom = &results[chrom_start].variant.chrom;
        let mut chrom_end = chrom_start;
        while chrom_end < n && results[chrom_end].variant.chrom == *chrom {
            chrom_end += 1;
        }
        ranges.push((chrom_start, chrom_end));
        chrom_start = chrom_end;
    }

    let bw = bandwidth as f64;
    let mut sum_k_squared = 0.0;
    let mut count = 0.0;

    for (start, end) in ranges.iter() {
        // Two-pointer window bounds
        let mut left = *start;
        let mut right = *start;

        for i in *start..*end {
            let center_pos = results[i].variant.pos as f64;

            // Advance left pointer
            while left < *end && (results[left].variant.pos as f64) < center_pos - bw {
                left += 1;
            }

            // Advance right pointer
            while right < *end && (results[right].variant.pos as f64) <= center_pos + bw {
                right += 1;
            }

            // Compute Σwⱼ and Σwⱼ² for this window
            let mut sum_w = 0.0;
            let mut sum_w_sq = 0.0;
            for j in left..right {
                let dist = (results[j].variant.pos as f64 - center_pos) / bw;
                let w = tricube(dist);
                sum_w += w;
                sum_w_sq += w * w;
            }

            // Normalized weights: kⱼ = wⱼ / Σw
            // Σkⱼ² = Σ(wⱼ²) / (Σw)²
            if sum_w > 0.0 {
                let sum_k_sq = sum_w_sq / (sum_w * sum_w);
                sum_k_squared += sum_k_sq;
            }
            count += 1.0;
        }
    }

    if count > 0.0 {
        sum_k_squared / count
    } else {
        1.0
    }
}

/// Smooth G-statistics using a tricube kernel, operating per-chromosome.
///
/// Assumes `results` are sorted by (chrom, pos) as they come from a sorted VCF.
/// Uses a two-pointer sliding window for O(n) per chromosome.
///
/// Returns a Vec<f64> of smoothed G' values, one per input result.
pub fn smooth_g_statistics(results: &[GStatisticResult], bandwidth: u64) -> Vec<f64> {
    let n = results.len();
    let mut g_prime = vec![0.0f64; n];

    // Collect chromosome boundary ranges
    let mut ranges: Vec<(usize, usize)> = Vec::new();
    let mut chrom_start = 0;
    while chrom_start < n {
        let chrom = &results[chrom_start].variant.chrom;
        let mut chrom_end = chrom_start;
        while chrom_end < n && results[chrom_end].variant.chrom == *chrom {
            chrom_end += 1;
        }
        ranges.push((chrom_start, chrom_end));
        chrom_start = chrom_end;
    }

    // Process each chromosome in parallel, each producing its own chunk
    let chunks: Vec<(usize, Vec<f64>)> = ranges
        .into_par_iter()
        .map(|(start, end)| {
            let mut chunk = vec![0.0f64; end - start];
            smooth_chromosome_into(results, &mut chunk, start, end, bandwidth);
            (start, chunk)
        })
        .collect();

    for (start, chunk) in chunks {
        g_prime[start..start + chunk.len()].copy_from_slice(&chunk);
    }

    g_prime
}

/// Smooth a single chromosome's worth of results using sliding window.
/// Writes into `out_slice` which corresponds to results[start..end].
fn smooth_chromosome_into(
    results: &[GStatisticResult],
    out_slice: &mut [f64],
    start: usize,
    end: usize,
    bandwidth: u64,
) {
    let bw = bandwidth as f64;

    // Two-pointer window bounds
    let mut left = start;
    let mut right = start;

    for i in start..end {
        let center_pos = results[i].variant.pos as f64;

        // Advance left pointer: skip positions too far left
        while left < end && (results[left].variant.pos as f64) < center_pos - bw {
            left += 1;
        }

        // Advance right pointer: include positions within bandwidth
        while right < end && (results[right].variant.pos as f64) <= center_pos + bw {
            right += 1;
        }

        // Compute weighted average over window [left, right)
        let mut sum_wg = 0.0;
        let mut sum_w = 0.0;

        for j in left..right {
            let dist = (results[j].variant.pos as f64 - center_pos) / bw;
            let w = tricube(dist);
            sum_wg += w * results[j].g_statistic;
            sum_w += w;
        }

        out_slice[i - start] = if sum_w > 0.0 { sum_wg / sum_w } else { results[i].g_statistic };
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::types::Variant;
    use approx::assert_relative_eq;

    fn make_result(chrom: &str, pos: u64, g: f64) -> GStatisticResult {
        GStatisticResult {
            variant: Variant {
                chrom: chrom.to_string(),
                pos,
                ref_allele: "A".to_string(),
                alt_allele: "T".to_string(),
                resistant_ref_depth: 50,
                resistant_alt_depth: 50,
                resistant_dp: 100,
                resistant_gq: 99,
                susceptible_ref_depth: 50,
                susceptible_alt_depth: 50,
                susceptible_dp: 100,
                susceptible_gq: 99,
            },
            g_statistic: g,
            snp_index_resistant: 0.5,
            snp_index_susceptible: 0.5,
            delta_snp_index: 0.0,
        }
    }

    #[test]
    fn test_tricube_at_zero() {
        assert_relative_eq!(tricube(0.0), 1.0, epsilon = 1e-10);
    }

    #[test]
    fn test_tricube_at_boundary() {
        assert_relative_eq!(tricube(1.0), 0.0, epsilon = 1e-10);
        assert_relative_eq!(tricube(-1.0), 0.0, epsilon = 1e-10);
    }

    #[test]
    fn test_tricube_outside() {
        assert_eq!(tricube(1.5), 0.0);
        assert_eq!(tricube(-2.0), 0.0);
    }

    #[test]
    fn test_single_snp() {
        let results = vec![make_result("chr1", 100, 5.0)];
        let g_prime = smooth_g_statistics(&results, 1_000_000);
        assert_relative_eq!(g_prime[0], 5.0, epsilon = 1e-10);
    }

    #[test]
    fn test_uniform_values_unchanged() {
        let results = vec![
            make_result("chr1", 100, 3.0),
            make_result("chr1", 200, 3.0),
            make_result("chr1", 300, 3.0),
        ];
        let g_prime = smooth_g_statistics(&results, 1_000_000);
        for gp in &g_prime {
            assert_relative_eq!(*gp, 3.0, epsilon = 1e-10);
        }
    }

    #[test]
    fn test_separate_chromosomes() {
        // Two chromosomes with very different G values
        let results = vec![
            make_result("chr1", 100, 1.0),
            make_result("chr1", 200, 1.0),
            make_result("chr2", 100, 100.0),
            make_result("chr2", 200, 100.0),
        ];
        let g_prime = smooth_g_statistics(&results, 1_000_000);
        // chr1 should stay near 1.0, chr2 near 100.0 - no cross-contamination
        assert_relative_eq!(g_prime[0], 1.0, epsilon = 0.1);
        assert_relative_eq!(g_prime[2], 100.0, epsilon = 0.1);
    }

    #[test]
    fn test_smoothing_reduces_spike() {
        // A spike surrounded by low values should be reduced
        let results = vec![
            make_result("chr1", 100, 1.0),
            make_result("chr1", 200, 1.0),
            make_result("chr1", 300, 100.0),
            make_result("chr1", 400, 1.0),
            make_result("chr1", 500, 1.0),
        ];
        let g_prime = smooth_g_statistics(&results, 1_000_000);
        // The spike at position 300 should be smoothed down
        assert!(g_prime[2] < 100.0);
        assert!(g_prime[2] > 1.0);
    }

    #[test]
    fn test_compute_smoothing_factor() {
        // Test that Σkⱼ² is computed correctly
        // Using 5 SNPs all at very similar positions with 1MB bandwidth
        // All SNPs are within the window, weights are equal
        let results = vec![
            make_result("chr1", 100, 1.0),
            make_result("chr1", 200, 1.0),
            make_result("chr1", 300, 1.0),
            make_result("chr1", 400, 1.0),
            make_result("chr1", 500, 1.0),
        ];
        let factor = compute_smoothing_factor(&results, 1_000_000);
        // With 5 SNPs all in the window, each gets weight ~1, normalized kⱼ = 1/5
        // Σkⱼ² = 5 × (1/5)² = 0.2
        assert_relative_eq!(factor, 0.2, epsilon = 0.01);
    }

    #[test]
    fn test_compute_smoothing_factor_sparse() {
        // Test with sparse SNPs - larger Σkⱼ² (fewer SNPs in window)
        let results = vec![
            make_result("chr1", 100, 1.0),
            make_result("chr1", 200, 1.0),
            make_result("chr1", 1000000, 1.0),  // Far apart
            make_result("chr1", 1000100, 1.0),
            make_result("chr1", 1000200, 1.0),
        ];
        let factor = compute_smoothing_factor(&results, 1000);  // 1kb bandwidth
        // With only 2 SNPs in each local window, Σkⱼ² ≈ 2 × (1/2)² = 0.5
        assert!(factor > 0.3 && factor < 0.7);
    }

    #[test]
    fn test_compute_smoothing_factor_empty() {
        let factor = compute_smoothing_factor(&[], 1_000_000);
        assert_relative_eq!(factor, 1.0, epsilon = 0.001);
    }
}
