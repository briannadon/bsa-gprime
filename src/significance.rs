use crate::types::{GStatisticResult, NullDistributionParams, SignificanceResult};
use rayon::prelude::*;
use statrs::distribution::{ContinuousCDF, LogNormal};

/// Estimate null distribution parameters using the parametric method (Magwene et al. 2011).
///
/// * `n_s` - effective population size per bulk (individuals * ploidy)
/// * `avg_coverage` - average read depth across all SNPs (both bulks combined)
pub fn estimate_null_parametric(n_s: f64, avg_coverage: f64) -> NullDistributionParams {
    let mean_raw = 1.0 + 1.0 / n_s + 1.0 / avg_coverage;
    let var_raw = 2.0 + 4.0 / n_s + 2.0 / avg_coverage;

    let variance_ratio = var_raw / (mean_raw * mean_raw);
    let mu = variance_ratio.ln_1p() * -0.5 + mean_raw.ln();
    let sigma = variance_ratio.ln_1p().sqrt();

    NullDistributionParams { mu, sigma }
}

/// Estimate null distribution parameters using the non-parametric method.
///
/// Uses Hampel's outlier rule with left-MAD to robustly estimate the
/// null distribution from observed smoothed G' values.
pub fn estimate_null_nonparametric(g_prime_values: &[f64]) -> NullDistributionParams {
    // Collect log(G') for G' > 0
    let mut log_g: Vec<f64> = g_prime_values
        .iter()
        .filter(|&&g| g > 0.0)
        .map(|&g| g.ln())
        .collect();

    if log_g.is_empty() {
        return NullDistributionParams {
            mu: 0.0,
            sigma: 1.0,
        };
    }

    log_g.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let median = compute_median(&log_g);

    // Left-MAD: median of |log_g_i - median| for values <= median
    let left_mad = compute_left_mad(&log_g, median);

    // Hampel's rule: outlier if (log_g_i - median) > 5.2 * left_MAD
    let threshold = 5.2 * left_mad;
    let trimmed: Vec<f64> = log_g
        .iter()
        .filter(|&&x| (x - median) <= threshold)
        .copied()
        .collect();

    if trimmed.is_empty() {
        return NullDistributionParams {
            mu: median,
            sigma: 1.4826 * left_mad,
        };
    }

    let trimmed_median = compute_median(&trimmed);
    let trimmed_left_mad = compute_left_mad(&trimmed, trimmed_median);

    NullDistributionParams {
        mu: trimmed_median,
        sigma: 1.4826 * trimmed_left_mad,
    }
}

fn compute_median(sorted: &[f64]) -> f64 {
    let n = sorted.len();
    if n == 0 {
        return 0.0;
    }
    if n % 2 == 0 {
        (sorted[n / 2 - 1] + sorted[n / 2]) / 2.0
    } else {
        sorted[n / 2]
    }
}

fn compute_left_mad(sorted: &[f64], median: f64) -> f64 {
    let mut deviations: Vec<f64> = sorted
        .iter()
        .filter(|&&x| x <= median)
        .map(|&x| (x - median).abs())
        .collect();
    deviations.sort_by(|a, b| a.partial_cmp(b).unwrap());
    if deviations.is_empty() {
        return f64::MIN_POSITIVE;
    }
    let mad = compute_median(&deviations);
    if mad == 0.0 {
        f64::MIN_POSITIVE
    } else {
        mad
    }
}

/// Compute p-values from smoothed G' values and null distribution parameters.
///
/// p_value = 1 - CDF(G') under LogNormal(mu, sigma)
pub fn compute_p_values(g_prime_values: &[f64], null: &NullDistributionParams) -> Vec<f64> {
    let dist = LogNormal::new(null.mu, null.sigma).expect("Invalid log-normal parameters");

    g_prime_values
        .par_iter()
        .map(|&gp| {
            if gp <= 0.0 {
                1.0
            } else {
                1.0 - dist.cdf(gp)
            }
        })
        .collect()
}

/// Benjamini-Hochberg FDR correction.
///
/// Returns q-values (adjusted p-values) in the same order as input.
pub fn benjamini_hochberg(p_values: &[f64]) -> Vec<f64> {
    let n = p_values.len();
    if n == 0 {
        return vec![];
    }

    // Create indexed list and sort by p-value
    let mut indexed: Vec<(usize, f64)> = p_values.iter().copied().enumerate().collect();
    indexed.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap());

    let mut q_values = vec![0.0; n];
    let n_f64 = n as f64;

    // Compute adjusted p-values
    let mut cummin = f64::INFINITY;
    for i in (0..n).rev() {
        let (orig_idx, p) = indexed[i];
        let rank = (i + 1) as f64;
        let adjusted = (p * n_f64 / rank).min(1.0);
        cummin = cummin.min(adjusted);
        q_values[orig_idx] = cummin;
    }

    q_values
}

/// Compute average coverage across all SNPs (sum of both bulks' DP / n).
pub fn compute_avg_coverage(results: &[GStatisticResult]) -> f64 {
    if results.is_empty() {
        return 0.0;
    }
    let total: f64 = results
        .iter()
        .map(|r| (r.variant.resistant_dp + r.variant.susceptible_dp) as f64)
        .sum();
    total / results.len() as f64
}

/// Full significance testing pipeline.
///
/// Takes raw G-statistic results and smoothed G' values, estimates the null
/// distribution, computes p-values and q-values, and returns SignificanceResults.
pub fn run_significance_pipeline(
    results: Vec<GStatisticResult>,
    g_prime_values: &[f64],
    null_params: &NullDistributionParams,
    fdr_threshold: f64,
) -> Vec<SignificanceResult> {
    let p_values = compute_p_values(g_prime_values, null_params);
    let q_values = benjamini_hochberg(&p_values);

    results
        .into_iter()
        .enumerate()
        .map(|(i, raw)| SignificanceResult {
            raw,
            g_prime: g_prime_values[i],
            p_value: p_values[i],
            q_value: q_values[i],
            significant: q_values[i] < fdr_threshold,
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_parametric_null_basic() {
        let params = estimate_null_parametric(100.0, 100.0);
        // mu is negative because the log-normal location parameter accounts for variance
        assert!(params.mu.is_finite());
        assert_relative_eq!(params.mu, -0.5262, epsilon = 0.01);
        assert!(params.sigma > 0.0);
        assert_relative_eq!(params.sigma, 1.045, epsilon = 0.01);
    }

    #[test]
    fn test_nonparametric_null_uniform() {
        // Uniform log-normal-ish values
        let values: Vec<f64> = (1..=100).map(|i| i as f64 * 0.1).collect();
        let params = estimate_null_nonparametric(&values);
        assert!(params.mu.is_finite());
        assert!(params.sigma > 0.0);
    }

    #[test]
    fn test_nonparametric_null_with_outliers() {
        // Most values around 1.0, a few extreme outliers
        let mut values: Vec<f64> = vec![1.0; 100];
        // Add some outliers
        values.extend(vec![1000.0; 5]);
        let params = estimate_null_nonparametric(&values);
        // mu should be close to ln(1.0) = 0
        assert_relative_eq!(params.mu, 0.0, epsilon = 0.1);
    }

    #[test]
    fn test_p_values_monotonic_with_g() {
        let null = NullDistributionParams {
            mu: 0.0,
            sigma: 1.0,
        };
        let g_values = vec![0.1, 1.0, 5.0, 10.0, 50.0];
        let p_vals = compute_p_values(&g_values, &null);

        // Higher G' should give smaller p-values
        for i in 1..p_vals.len() {
            assert!(p_vals[i] <= p_vals[i - 1], "p-values should decrease with increasing G'");
        }
    }

    #[test]
    fn test_p_values_in_range() {
        let null = NullDistributionParams {
            mu: 0.0,
            sigma: 1.0,
        };
        let g_values = vec![0.0, 0.5, 1.0, 5.0, 100.0];
        let p_vals = compute_p_values(&g_values, &null);
        for p in &p_vals {
            assert!(*p >= 0.0 && *p <= 1.0, "p-value out of range: {}", p);
        }
    }

    #[test]
    fn test_bh_fdr_basic() {
        let p_values = vec![0.01, 0.04, 0.03, 0.50, 0.90];
        let q_values = benjamini_hochberg(&p_values);

        // q-values should be >= corresponding p-values
        for (p, q) in p_values.iter().zip(q_values.iter()) {
            assert!(*q >= *p, "q-value {} < p-value {}", q, p);
        }

        // q-values should be in [0, 1]
        for q in &q_values {
            assert!(*q >= 0.0 && *q <= 1.0, "q-value out of range: {}", q);
        }
    }

    #[test]
    fn test_bh_fdr_preserves_order() {
        let p_values = vec![0.001, 0.01, 0.05, 0.10, 0.50];
        let q_values = benjamini_hochberg(&p_values);

        // The smallest p-value should yield the smallest q-value
        let min_p_idx = p_values
            .iter()
            .enumerate()
            .min_by(|a, b| a.1.partial_cmp(b.1).unwrap())
            .unwrap()
            .0;
        let min_q_idx = q_values
            .iter()
            .enumerate()
            .min_by(|a, b| a.1.partial_cmp(b.1).unwrap())
            .unwrap()
            .0;
        assert_eq!(min_p_idx, min_q_idx);
    }

    #[test]
    fn test_bh_fdr_empty() {
        assert!(benjamini_hochberg(&[]).is_empty());
    }

    #[test]
    fn test_compute_median() {
        assert_relative_eq!(compute_median(&[1.0, 2.0, 3.0]), 2.0);
        assert_relative_eq!(compute_median(&[1.0, 2.0, 3.0, 4.0]), 2.5);
        assert_relative_eq!(compute_median(&[5.0]), 5.0);
    }
}
