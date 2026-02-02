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
