use crate::types::{GStatisticResult, NullDistributionParams, SignificanceResult};
use rayon::prelude::*;
use statrs::distribution::{ContinuousCDF, LogNormal};

/// Estimate null distribution parameters using the parametric method (Magwene et al. 2011).
///
/// This function estimates the null distribution for the **raw G-statistic** (unsmoothed).
/// For smoothed G' values, use [`estimate_null_parametric_gprime`] instead.
///
/// * `n_s` - effective population size per bulk (individuals * ploidy)
/// * `avg_coverage` - average read depth per bulk (not combined across both bulks)
///
/// Equations 8-9 from Magwene et al. (2011):
/// - E[G] ≈ 1 + C/(2n_s)
/// - Var[G] ≈ 2 + 1/(2C) + (1+2C)/n_s + C(4n_s-1)/(8n_s³)
pub fn estimate_null_parametric(n_s: f64, avg_coverage: f64) -> NullDistributionParams {
    // Equation 8: E[G] = 1 + C/(2n_s)
    let mean_g = 1.0 + avg_coverage / (2.0 * n_s);
    
    // Equation 9: Var[G] = 2 + 1/(2C) + (1+2C)/n_s + C(4n_s-1)/(8n_s³)
    let var_g = 2.0
        + 1.0 / (2.0 * avg_coverage)
        + (1.0 + 2.0 * avg_coverage) / n_s
        + avg_coverage * (4.0 * n_s - 1.0) / (8.0 * n_s.powi(3));
    
    // Convert to log-normal parameters
    let variance_ratio = var_g / (mean_g * mean_g);
    let mu = variance_ratio.ln_1p() * -0.5 + mean_g.ln();
    let sigma = variance_ratio.ln_1p().sqrt();

    NullDistributionParams { mu, sigma }
}

/// Estimate null distribution parameters for smoothed G' values using the parametric method.
///
/// This function accounts for the variance reduction from smoothing by applying
/// a scaling factor based on the smoothing kernel.
///
/// * `n_s` - effective population size per bulk (individuals * ploidy)
/// * `avg_coverage` - average read depth per bulk (not combined across both bulks)
/// * `smoothing_factor` - scaling factor for variance reduction from smoothing.
///   For a tricube kernel, this is approximately Σkⱼ² ≈ 0.5-0.7 depending on SNP density.
///   If None, uses 1.0 (no smoothing correction).
///
/// Equation 12 from Magwene et al. (2011):
/// - Var[G'] ≈ Var[G] × Σkⱼ²
pub fn estimate_null_parametric_gprime(
    n_s: f64,
    avg_coverage: f64,
    smoothing_factor: Option<f64>,
) -> NullDistributionParams {
    // Equation 8: E[G] = 1 + C/(2n_s)
    // Mean is approximately unchanged by smoothing
    let mean_g = 1.0 + avg_coverage / (2.0 * n_s);
    
    // Equation 9: Var[G] = 2 + 1/(2C) + (1+2C)/n_s + C(4n_s-1)/(8n_s³)
    let var_g = 2.0
        + 1.0 / (2.0 * avg_coverage)
        + (1.0 + 2.0 * avg_coverage) / n_s
        + avg_coverage * (4.0 * n_s - 1.0) / (8.0 * n_s.powi(3));
    
    // Apply variance reduction from smoothing (equation 12)
    // smoothing_factor = Σkⱼ² ≈ 0.5-0.7 for typical tricube kernels
    let scale = smoothing_factor.unwrap_or(1.0);
    let var_gprime = var_g * scale;
    
    // Convert to log-normal parameters
    let variance_ratio = var_gprime / (mean_g * mean_g);
    let mu = variance_ratio.ln_1p() * -0.5 + mean_g.ln();
    let sigma = variance_ratio.ln_1p().sqrt();

    NullDistributionParams { mu, sigma }
}

/// Estimate null distribution parameters using the robust non-parametric method.
///
/// Implements the "robust empirical estimator" from Magwene et al. (2011) step 3b.
/// This method infers the null distribution from observed G' data by:
/// 1. Log-transforming G' values
/// 2. Computing median and left-MAD
/// 3. Applying Hampel's rule to identify and remove outliers (QTL regions)
/// 4. Estimating mode from trimmed data
/// 5. Computing log-normal parameters from median and mode
///
/// The observed G' data is a mixture of:
/// - Null distribution (non-QTL regions) - log-normal
/// - Contaminating distributions (QTL regions) - higher means
///
/// By removing outliers, we isolate the null distribution and estimate its
/// parameters robustly without requiring prior knowledge of n_s or C.
///
/// References:
/// - Magwene et al. (2011) PLoS Comput Biol 7(11): e1002255
/// - Bickel DR, Fruehwirth R (2006) Comput Stat Data An 50: 3500-3530
pub fn estimate_null_robust(g_prime_values: &[f64]) -> NullDistributionParams {
    // Step 1: Log-transform the G' values
    // W_G' = ln(X_G')
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

    // Keep a copy of original (non-log) values for mode estimation
    let original_values: Vec<f64> = g_prime_values
        .iter()
        .filter(|&&g| g > 0.0)
        .copied()
        .collect();

    log_g.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let mut orig_sorted = original_values.clone();
    orig_sorted.sort_by(|a, b| a.partial_cmp(b).unwrap());

    // Step 2: Compute median and left-MAD
    // s_W = MAD_l(W_G') = Median(|w_i - Median(W_G')|) for w_i <= Median
    let median_log = compute_median(&log_g);
    let left_mad = compute_left_mad(&log_g, median_log);

    // Step 3: Apply Hampel's rule to identify outliers
    // Outliers are w_i such that w_i > Median + g(N,α) * MAD_l
    // where g(N,α) ≈ 5.2 for normally distributed data
    let hampel_threshold = 5.2 * left_mad;
    
    // Find the cutoff: values above (median + 5.2 * MAD) are outliers
    let log_cutoff = median_log + hampel_threshold;
    let orig_cutoff = log_cutoff.exp();

    // Step 4: Trim outliers from both log-transformed and original values
    let trimmed_log: Vec<f64> = log_g
        .iter()
        .filter(|&&w| w <= log_cutoff)
        .copied()
        .collect();

    let trimmed_orig: Vec<f64> = orig_sorted
        .iter()
        .filter(|&&g| g <= orig_cutoff)
        .copied()
        .collect();

    if trimmed_log.len() < 3 || trimmed_orig.len() < 3 {
        // Not enough data after trimming, fall back to MAD-based estimation
        return NullDistributionParams {
            mu: median_log,
            sigma: 1.4826 * left_mad,
        };
    }

    // Compute statistics from trimmed data
    let trimmed_log_median = compute_median(&trimmed_log);
    let trimmed_log_mad = compute_left_mad(&trimmed_log, trimmed_log_median);

    // Estimate mode from trimmed original values using kernel density estimation
    // Mode_r(X_T) - robust mode estimator from Bickel & Fruehwirth (2006)
    let mode = estimate_mode_robust(&trimmed_orig);

    // Step 5: Compute log-normal parameters
    // For a log-normal distribution:
    // - Median = exp(μ) => μ = ln(Median)
    // - Mode = exp(μ - σ²) => σ² = μ - ln(Mode)
    //
    // From the trimmed log data:
    // μ = trimmed_log_median
    //
    // From mode estimation on original trimmed data:
    // σ² = μ - ln(Mode)
    //
    // Convert to μ and σ for log-normal parameterization
    let mu = trimmed_log_median;
    let sigma_sq = (mu - mode.ln()).abs();
    let sigma = sigma_sq.sqrt();

    // Ensure sigma is valid
    if sigma.is_nan() || sigma == 0.0 {
        NullDistributionParams {
            mu: trimmed_log_median,
            sigma: 1.4826 * trimmed_log_mad,
        }
    } else {
        NullDistributionParams { mu, sigma }
    }
}

/// Robust mode estimator using kernel density estimation.
///
/// Implements a simplified version of the mode estimator from:
/// Bickel DR, Fruehwirth R (2006) "On a fast, robust estimator of the mode"
/// Computational Statistics & Data Analysis 50: 3500-3530
///
/// Uses a Gaussian kernel with Silverman's rule of thumb for bandwidth selection.
fn estimate_mode_robust(values: &[f64]) -> f64 {
    let n = values.len();
    if n == 0 {
        return 1.0; // Default mode
    }

    // Compute mean and standard deviation
    let mean: f64 = values.iter().sum::<f64>() / n as f64;
    let variance: f64 = values
        .iter()
        .map(|v| (v - mean).powi(2))
        .sum::<f64>()
        / n as f64;
    let std = variance.sqrt();

    if std == 0.0 {
        return mean;
    }

    // Silverman's rule of thumb for bandwidth
    // h = 1.06 * σ * n^(-1/5)
    let h = 1.06 * std * (n as f64).powf(-0.2);

    // Find the value that maximizes the kernel density
    // Use the data range as search bounds
    let min_val = values[0];
    let max_val = values[n - 1];

    // Grid search for mode
    let num_points = 100.min(n);
    let step = (max_val - min_val) / (num_points as f64 - 1.0);

    let mut max_density = f64::MIN_POSITIVE;
    let mut mode_estimate = mean;

    for i in 0..num_points {
        let x = min_val + (i as f64) * step;
        // Gaussian kernel
        let density: f64 = values
            .iter()
            .map(|v| {
                let z = (v - x) / h;
                (-0.5 * z * z).exp() / (h * (2.0 * std::f64::consts::PI).sqrt())
            })
            .sum();

        if density > max_density {
            max_density = density;
            mode_estimate = x;
        }
    }

    mode_estimate
}

/// Legacy non-parametric estimator (uses MAD instead of mode).
///
/// This method uses the left-MAD to estimate sigma rather than mode estimation.
/// For the full robust estimator with mode estimation, use [`estimate_null_robust`].
#[deprecated(since = "0.1.0", note = "Use estimate_null_robust instead for full paper compliance")]
pub fn estimate_null_nonparametric(g_prime_values: &[f64]) -> NullDistributionParams {
    // Step 1: Log-transform the G' values
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

    // Step 2: Left-MAD computation
    let left_mad = compute_left_mad(&log_g, median);

    // Step 3: Hampel's rule
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
        // Test with n_s=100, avg_coverage=100 (per bulk)
        // Equation 8: E[G] = 1 + C/(2n_s) = 1 + 100/(2*100) = 1.5
        // Equation 9: Var[G] = 2 + 1/(2C) + (1+2C)/n_s + C(4n_s-1)/(8n_s³)
        //             = 2 + 1/200 + 201/100 + 100*399/(8*1000000)
        //             = 2 + 0.005 + 2.01 + 0.0049875 = 4.0199875
        let params = estimate_null_parametric(100.0, 100.0);
        // mu is negative because the log-normal location parameter accounts for variance
        assert!(params.mu.is_finite());
        assert_relative_eq!(params.mu, -0.107, epsilon = 0.01);
        assert!(params.sigma > 0.0);
        assert_relative_eq!(params.sigma, 1.012, epsilon = 0.01);
    }

    #[test]
    fn test_parametric_null_gprime() {
        // Test G' parametric estimation with smoothing factor
        // With smoothing_factor=0.6, variance should be reduced
        let params_raw = estimate_null_parametric(100.0, 100.0);
        let params_gprime = estimate_null_parametric_gprime(100.0, 100.0, Some(0.6));
        
        // G' should have smaller sigma due to smoothing
        assert!(params_gprime.sigma < params_raw.sigma);
        assert!(params_gprime.sigma > 0.0);
    }

    #[test]
    fn test_parametric_null_gprime_default() {
        // Test G' parametric estimation with default smoothing factor (1.0)
        // Should be same as raw G
        let params_raw = estimate_null_parametric(100.0, 100.0);
        let params_gprime = estimate_null_parametric_gprime(100.0, 100.0, None);
        
        assert_relative_eq!(params_raw.mu, params_gprime.mu, epsilon = 0.001);
        assert_relative_eq!(params_raw.sigma, params_gprime.sigma, epsilon = 0.001);
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

    #[test]
    fn test_estimate_null_robust_basic() {
        // Generate values from a known log-normal distribution
        // LogNormal(mu=0.0, sigma=0.5) => median=1.0, mode=exp(-0.25)≈0.78
        // Create synthetic log-normal data using Box-Muller transform
        let mut values: Vec<f64> = Vec::with_capacity(1000);
        let mut seed = 12345u64;
        for _ in 0..1000 {
            // Simple pseudo-random for reproducibility
            seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1);
            let x = (seed as f64 / u64::MAX as f64) * 2.0 * std::f64::consts::PI;
            seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1);
            let y = (seed as f64 / u64::MAX as f64);
            let z = (-2.0 * y.ln()).sqrt() * x.cos();
            values.push((0.0 + 0.5 * z).exp());
        }
        
        let params = estimate_null_robust(&values);
        // mu should be close to 0.0 (ln of median for LogNormal(0, 0.5))
        assert_relative_eq!(params.mu, 0.0, epsilon = 0.2);
        // sigma should be close to 0.5
        assert_relative_eq!(params.sigma, 0.5, epsilon = 0.2);
    }

    #[test]
    fn test_estimate_null_robust_with_outliers() {
        // 95% null values (log-normal around G'=1-2)
        let mut values: Vec<f64> = Vec::with_capacity(960);
        let mut seed = 54321u64;
        for _ in 0..950 {
            seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1);
            let x = (seed as f64 / u64::MAX as f64) * 2.0 * std::f64::consts::PI;
            seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1);
            let y = (seed as f64 / u64::MAX as f64);
            let z = (-2.0 * y.ln()).sqrt() * x.cos();
            values.push((0.0 + 0.3 * z).exp());
        }
        
        // 5% QTL-like outliers (very high G' values)
        values.extend(vec![10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0, 45.0, 50.0]);
        
        let params = estimate_null_robust(&values);
        // mu should be estimated from the null bulk (around 0)
        // not pulled up by the outliers
        assert!(params.mu < 2.0, "mu should not be inflated by outliers: mu={}", params.mu);
        assert!(params.sigma > 0.0);
    }

    #[test]
    fn test_estimate_null_robust_empty() {
        let params = estimate_null_robust(&[]);
        assert_eq!(params.mu, 0.0);
        assert_eq!(params.sigma, 1.0);
    }

    #[test]
    fn test_estimate_null_robust_all_zeros() {
        let params = estimate_null_robust(&[0.0, 0.0, 0.0]);
        // All zeros should return defaults
        assert!(params.mu.is_finite());
        assert!(params.sigma > 0.0);
    }

    #[test]
    fn test_estimate_mode_robust_basic() {
        // Mode of uniform-ish distribution
        let values: Vec<f64> = vec![1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0];
        let mode = estimate_mode_robust(&values);
        // Mode should be somewhere in the middle
        assert!(mode >= 1.0 && mode <= 4.0);
    }

    #[test]
    fn test_estimate_mode_robust_constant() {
        // All same values
        let values: Vec<f64> = vec![5.0; 100];
        let mode = estimate_mode_robust(&values);
        assert_relative_eq!(mode, 5.0, epsilon = 0.01);
    }

    #[test]
    fn test_robust_vs_legacy() {
        // Both methods should give similar results for clean data
        let values: Vec<f64> = (1..=100).map(|i| i as f64 * 0.1).collect();
        
        let robust = estimate_null_robust(&values);
        let legacy = estimate_null_nonparametric(&values);
        
        // Should be reasonably close for clean data
        assert_relative_eq!(robust.mu, legacy.mu, epsilon = 0.3);
    }

    #[test]
    fn test_hampel_threshold_applied() {
        // Create data with clear outliers
        let mut values: Vec<f64> = vec![];
        // Null region: values around 1.0-2.0
        for i in 0..100 {
            let val = 1.0 + (i as f64 / 10.0).floor() * 0.1;
            values.push(val);
        }
        // Outliers: very high values
        values.extend(vec![100.0; 20]);
        
        let params = estimate_null_robust(&values);
        // After trimming outliers, mu should be close to ln(1-2) ≈ 0.3
        assert!(params.mu < 1.0, "Outliers should be trimmed: mu={}", params.mu);
    }
}
