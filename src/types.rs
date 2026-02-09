/// Parental genotype information for a variant site
#[derive(Debug, Clone)]
pub struct ParentalInfo {
    pub parent_high_gt: String,  // e.g., "0/0" or "1/1"
    pub parent_low_gt: String,
    pub alleles_swapped: bool,   // Whether REF/ALT depths were swapped
}

/// Represents a single variant with allele counts
#[derive(Debug, Clone)]
pub struct Variant {
    pub chrom: String,
    pub pos: u64,
    pub ref_allele: String,
    pub alt_allele: String,

    // Resistant bulk
    pub high_ref_depth: u32,
    pub high_alt_depth: u32,
    pub high_dp: u32,
    pub high_gq: u32,

    // Susceptible bulk
    pub low_ref_depth: u32,
    pub low_alt_depth: u32,
    pub low_dp: u32,
    pub low_gq: u32,

    // Parental info (None when running without parental samples)
    pub parental_info: Option<ParentalInfo>,
}

/// Results for a single variant
#[derive(Debug, Clone)]
pub struct GStatisticResult {
    pub variant: Variant,
    pub g_statistic: f64,
    pub snp_index_high: f64,
    pub snp_index_low: f64,
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
