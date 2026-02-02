use crate::types::GStatisticResult;
use anyhow::Result;
use csv::Writer;
use std::path::Path;

pub fn write_results(results: &[GStatisticResult], path: &Path) -> Result<()> {
    let mut wtr = Writer::from_path(path)?;

    // Write header
    wtr.write_record(&[
        "chrom",
        "pos",
        "ref",
        "alt",
        "resistant_ref_depth",
        "resistant_alt_depth",
        "resistant_dp",
        "resistant_gq",
        "susceptible_ref_depth",
        "susceptible_alt_depth",
        "susceptible_dp",
        "susceptible_gq",
        "g_statistic",
        "snp_index_resistant",
        "snp_index_susceptible",
        "delta_snp_index",
    ])?;

    // Write data
    for result in results {
        wtr.write_record(&[
            &result.variant.chrom,
            &result.variant.pos.to_string(),
            &result.variant.ref_allele,
            &result.variant.alt_allele,
            &result.variant.resistant_ref_depth.to_string(),
            &result.variant.resistant_alt_depth.to_string(),
            &result.variant.resistant_dp.to_string(),
            &result.variant.resistant_gq.to_string(),
            &result.variant.susceptible_ref_depth.to_string(),
            &result.variant.susceptible_alt_depth.to_string(),
            &result.variant.susceptible_dp.to_string(),
            &result.variant.susceptible_gq.to_string(),
            &format!("{:.6}", result.g_statistic),
            &format!("{:.6}", result.snp_index_resistant),
            &format!("{:.6}", result.snp_index_susceptible),
            &format!("{:.6}", result.delta_snp_index),
        ])?;
    }

    wtr.flush()?;
    Ok(())
}
