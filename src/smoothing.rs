use crate::types::GStatisticResult;
use indicatif::ProgressBar;
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

/// Smoothing factors needed for equation 12 (Var[G'] estimation).
pub struct SmoothingFactors {
    /// Average Σk_j² across all SNPs (diagonal variance term)
    pub sum_k_squared: f64,
    /// Average Σ_{i≠j} (1-2r_ij)² k_i k_j across all SNPs (covariance term from linkage)
    pub covariance_weight: f64,
}

/// Compute smoothing factors for equation 12 from Magwene et al. (2011):
///
/// Var[G'] = Var[G] × Σk_j² + C²(4n_s-1)/(8n_s³) × Σ_{i≠j}(1-2r_ij)²k_ik_j
///
/// The tricube kernel weights are: k_j = (1 - D_j³)³ / S_W
/// where S_W = Σ(1 - D_j³)³ is the normalization factor.
///
/// * `recomb_rate_cm_per_mb` - recombination rate in cM/Mb, used to compute
///   pairwise recombination fractions r_ij via the Haldane mapping function.
///   If None, only the diagonal term (Σk_j²) is computed.
///
/// Returns `SmoothingFactors` with both the diagonal and covariance components,
/// averaged across all SNPs.
pub fn compute_smoothing_factors(
    results: &[GStatisticResult],
    bandwidth: u64,
    recomb_rate_cm_per_mb: Option<f64>,
) -> SmoothingFactors {
    let n = results.len();
    if n == 0 {
        return SmoothingFactors {
            sum_k_squared: 1.0,
            covariance_weight: 0.0,
        };
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

    // Pre-compute Haldane factor for (1-2r)² = exp(-4 * d_Morgans)
    // d_Morgans = d_bp * recomb_rate / 1e8  (cM/Mb × bp → cM → Morgans)
    let haldane_factor = recomb_rate_cm_per_mb.map(|rate| 4.0 * rate / 1e8);

    let mut total_sum_k_sq = 0.0;
    let mut total_cov_weight = 0.0;
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

            // Compute unnormalized weights for this window
            let window_size = right - left;
            let mut weights: Vec<f64> = Vec::with_capacity(window_size);
            let mut sum_w = 0.0;
            for j in left..right {
                let dist = (results[j].variant.pos as f64 - center_pos) / bw;
                let w = tricube(dist);
                weights.push(w);
                sum_w += w;
            }

            if sum_w > 0.0 {
                let inv_sum_w = 1.0 / sum_w;

                // Diagonal term: Σk_j² = Σ(w_j/Σw)²
                let sum_k_sq: f64 = weights.iter().map(|&w| {
                    let k = w * inv_sum_w;
                    k * k
                }).sum();
                total_sum_k_sq += sum_k_sq;

                // Covariance term: Σ_{j≠l} (1-2r_jl)² k_j k_l
                // Using Haldane mapping: (1-2r)² = exp(-4 * d_Morgans)
                if let Some(hf) = haldane_factor {
                    let mut cov = 0.0;
                    for j_idx in 0..window_size {
                        let k_j = weights[j_idx] * inv_sum_w;
                        let pos_j = results[left + j_idx].variant.pos as f64;
                        for l_idx in (j_idx + 1)..window_size {
                            let k_l = weights[l_idx] * inv_sum_w;
                            let d_bp = (pos_j - results[left + l_idx].variant.pos as f64).abs();
                            let linkage_sq = (-hf * d_bp).exp();
                            cov += linkage_sq * k_j * k_l;
                        }
                    }
                    // Factor of 2 for symmetry (j<l counted once, but sum is over all i≠j)
                    total_cov_weight += 2.0 * cov;
                }
            }
            count += 1.0;
        }
    }

    SmoothingFactors {
        sum_k_squared: if count > 0.0 { total_sum_k_sq / count } else { 1.0 },
        covariance_weight: if count > 0.0 { total_cov_weight / count } else { 0.0 },
    }
}

/// Smooth G-statistics using a tricube kernel, operating per-chromosome.
///
/// Assumes `results` are sorted by (chrom, pos) as they come from a sorted VCF.
/// Uses a two-pointer sliding window for O(n) per chromosome.
///
/// Returns a Vec<f64> of smoothed G' values, one per input result.
pub fn smooth_g_statistics(
    results: &[GStatisticResult],
    bandwidth: u64,
    pb: Option<&ProgressBar>,
) -> Vec<f64> {
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
            smooth_chromosome_into(results, &mut chunk, start, end, bandwidth, pb);
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
    pb: Option<&ProgressBar>,
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

        if let Some(pb) = pb {
            pb.inc(1);
        }
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
        let g_prime = smooth_g_statistics(&results, 1_000_000, None);
        assert_relative_eq!(g_prime[0], 5.0, epsilon = 1e-10);
    }

    #[test]
    fn test_uniform_values_unchanged() {
        let results = vec![
            make_result("chr1", 100, 3.0),
            make_result("chr1", 200, 3.0),
            make_result("chr1", 300, 3.0),
        ];
        let g_prime = smooth_g_statistics(&results, 1_000_000, None);
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
        let g_prime = smooth_g_statistics(&results, 1_000_000, None);
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
        let g_prime = smooth_g_statistics(&results, 1_000_000, None);
        // The spike at position 300 should be smoothed down
        assert!(g_prime[2] < 100.0);
        assert!(g_prime[2] > 1.0);
    }

    #[test]
    fn test_compute_smoothing_factors() {
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
        let factors = compute_smoothing_factors(&results, 1_000_000, None);
        // With 5 SNPs all in the window, each gets weight ~1, normalized kⱼ = 1/5
        // Σkⱼ² = 5 × (1/5)² = 0.2
        assert_relative_eq!(factors.sum_k_squared, 0.2, epsilon = 0.01);
        // No recombination rate provided, covariance_weight should be 0
        assert_relative_eq!(factors.covariance_weight, 0.0, epsilon = 1e-10);
    }

    #[test]
    fn test_compute_smoothing_factors_sparse() {
        // Test with sparse SNPs - larger Σkⱼ² (fewer SNPs in window)
        let results = vec![
            make_result("chr1", 100, 1.0),
            make_result("chr1", 200, 1.0),
            make_result("chr1", 1000000, 1.0),  // Far apart
            make_result("chr1", 1000100, 1.0),
            make_result("chr1", 1000200, 1.0),
        ];
        let factors = compute_smoothing_factors(&results, 1000, None);  // 1kb bandwidth
        // With only 2 SNPs in each local window, Σkⱼ² ≈ 2 × (1/2)² = 0.5
        assert!(factors.sum_k_squared > 0.3 && factors.sum_k_squared < 0.7);
    }

    #[test]
    fn test_compute_smoothing_factors_empty() {
        let factors = compute_smoothing_factors(&[], 1_000_000, None);
        assert_relative_eq!(factors.sum_k_squared, 1.0, epsilon = 0.001);
        assert_relative_eq!(factors.covariance_weight, 0.0, epsilon = 1e-10);
    }

    #[test]
    fn test_compute_smoothing_factors_with_recombination() {
        // Test covariance term with recombination rate
        // 5 tightly linked SNPs (100bp apart, high recomb rate won't matter much)
        let results = vec![
            make_result("chr1", 100, 1.0),
            make_result("chr1", 200, 1.0),
            make_result("chr1", 300, 1.0),
            make_result("chr1", 400, 1.0),
            make_result("chr1", 500, 1.0),
        ];
        // Very low recombination rate: SNPs are tightly linked
        // (1-2r)² ≈ 1 for all pairs, so covariance_weight ≈ Σ_{i≠j} k_i k_j
        // = (Σk_i)² - Σk_i² = 1 - 0.2 = 0.8
        let factors_linked = compute_smoothing_factors(&results, 1_000_000, Some(0.001));
        assert!(factors_linked.covariance_weight > 0.0);
        assert_relative_eq!(factors_linked.covariance_weight, 0.8, epsilon = 0.01);
        // sum_k_squared unchanged
        assert_relative_eq!(factors_linked.sum_k_squared, 0.2, epsilon = 0.01);

        // Very high recombination rate: SNPs are unlinked
        // (1-2r)² ≈ 0 for all pairs, covariance_weight ≈ 0
        let factors_unlinked = compute_smoothing_factors(&results, 1_000_000, Some(1e6));
        assert_relative_eq!(factors_unlinked.covariance_weight, 0.0, epsilon = 0.01);
    }
}
