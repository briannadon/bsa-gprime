/// Represents a single variant with allele counts
#[derive(Debug, Clone)]
pub struct Variant {
    pub chrom: String,
    pub pos: u64,
    pub ref_allele: String,
    pub alt_allele: String,

    // Resistant bulk
    pub resistant_ref_depth: u32,
    pub resistant_alt_depth: u32,
    pub resistant_dp: u32,
    pub resistant_gq: u32,

    // Susceptible bulk
    pub susceptible_ref_depth: u32,
    pub susceptible_alt_depth: u32,
    pub susceptible_dp: u32,
    pub susceptible_gq: u32,
}

/// Results for a single variant
#[derive(Debug, Clone)]
pub struct GStatisticResult {
    pub variant: Variant,
    pub g_statistic: f64,
    pub snp_index_resistant: f64,
    pub snp_index_susceptible: f64,
    pub delta_snp_index: f64,
}

/// Parameters of the estimated log-normal null distribution
#[derive(Debug, Clone)]
pub struct NullDistributionParams {
    pub mu: f64,    // log-normal location (mean of underlying normal)
    pub sigma: f64, // log-normal scale (std dev of underlying normal)
}

/// Extended result with significance testing
#[derive(Debug, Clone)]
pub struct SignificanceResult {
    pub raw: GStatisticResult,
    pub g_prime: f64,      // smoothed G-statistic
    pub p_value: f64,      // from log-normal null
    pub q_value: f64,      // BH-adjusted p-value
    pub significant: bool, // q_value < fdr_threshold
}
