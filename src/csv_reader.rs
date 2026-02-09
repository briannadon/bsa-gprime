use crate::types::{GStatisticResult, ParentalInfo, SignificanceResult, Variant};
use anyhow::{Context, Result};
use std::path::Path;

/// Results loaded from a CSV file, auto-detecting whether significance columns are present.
pub enum LoadedResults {
    /// Basic results without significance testing (16 columns, optionally +3 parental)
    Basic(Vec<GStatisticResult>),
    /// Results with significance testing (20 columns, optionally +3 parental)
    WithSignificance(Vec<SignificanceResult>),
}

/// Read a results CSV file, auto-detecting whether significance columns are present.
///
/// Detection logic: if headers include `p_value`, `q_value`, and `significant`,
/// parse significance columns. If headers include `parent_high_gt`, also parse
/// parental columns. Backward compatible with old 16/20-column CSVs.
pub fn load_results_csv(path: &Path) -> Result<LoadedResults> {
    let mut rdr = csv::Reader::from_path(path)
        .with_context(|| format!("Failed to open CSV file: {}", path.display()))?;

    let headers = rdr.headers()?.clone();
    let has_significance = headers.iter().any(|h| h == "p_value")
        && headers.iter().any(|h| h == "q_value")
        && headers.iter().any(|h| h == "significant");
    let has_parental = headers.iter().any(|h| h == "parent_high_gt");

    if has_significance {
        let mut results = Vec::new();
        for (i, record) in rdr.records().enumerate() {
            let record = record.with_context(|| format!("Failed to parse CSV row {}", i + 1))?;
            let sig = parse_significance_record(&record, i + 1, has_parental)?;
            results.push(sig);
        }
        eprintln!("Loaded {} variants with significance data from {}", results.len(), path.display());
        Ok(LoadedResults::WithSignificance(results))
    } else {
        let mut results = Vec::new();
        for (i, record) in rdr.records().enumerate() {
            let record = record.with_context(|| format!("Failed to parse CSV row {}", i + 1))?;
            let g = parse_basic_record(&record, i + 1, has_parental)?;
            results.push(g);
        }
        eprintln!("Loaded {} variants (no significance data) from {}", results.len(), path.display());
        Ok(LoadedResults::Basic(results))
    }
}

fn parse_parental_columns(record: &csv::StringRecord, start_col: usize, row: usize) -> Result<Option<ParentalInfo>> {
    let ctx = || format!("row {}", row);

    let parent_high_gt = match record.get(start_col) {
        Some(v) => v.to_string(),
        None => return Ok(None),
    };
    let parent_low_gt = record
        .get(start_col + 1)
        .with_context(ctx)?
        .to_string();
    let alleles_swapped = record
        .get(start_col + 2)
        .with_context(ctx)?
        .parse::<bool>()
        .with_context(ctx)?;

    Ok(Some(ParentalInfo {
        parent_high_gt,
        parent_low_gt,
        alleles_swapped,
    }))
}

fn parse_basic_record(record: &csv::StringRecord, row: usize, has_parental: bool) -> Result<GStatisticResult> {
    let ctx = || format!("row {}", row);

    let parental_info = if has_parental {
        // Parental columns start at index 16 for basic records
        parse_parental_columns(record, 16, row)?
    } else {
        None
    };

    let variant = Variant {
        chrom: record.get(0).with_context(ctx)?.to_string(),
        pos: record.get(1).with_context(ctx)?.parse().with_context(ctx)?,
        ref_allele: record.get(2).with_context(ctx)?.to_string(),
        alt_allele: record.get(3).with_context(ctx)?.to_string(),
        high_ref_depth: record.get(4).with_context(ctx)?.parse().with_context(ctx)?,
        high_alt_depth: record.get(5).with_context(ctx)?.parse().with_context(ctx)?,
        high_dp: record.get(6).with_context(ctx)?.parse().with_context(ctx)?,
        high_gq: record.get(7).with_context(ctx)?.parse().with_context(ctx)?,
        low_ref_depth: record.get(8).with_context(ctx)?.parse().with_context(ctx)?,
        low_alt_depth: record.get(9).with_context(ctx)?.parse().with_context(ctx)?,
        low_dp: record.get(10).with_context(ctx)?.parse().with_context(ctx)?,
        low_gq: record.get(11).with_context(ctx)?.parse().with_context(ctx)?,
        parental_info,
    };

    Ok(GStatisticResult {
        variant,
        g_statistic: record.get(12).with_context(ctx)?.parse().with_context(ctx)?,
        snp_index_high: record.get(13).with_context(ctx)?.parse().with_context(ctx)?,
        snp_index_low: record.get(14).with_context(ctx)?.parse().with_context(ctx)?,
        delta_snp_index: record.get(15).with_context(ctx)?.parse().with_context(ctx)?,
    })
}

fn parse_significance_record(record: &csv::StringRecord, row: usize, has_parental: bool) -> Result<SignificanceResult> {
    let ctx = || format!("row {}", row);

    let parental_info = if has_parental {
        // Parental columns start at index 20 for significance records
        parse_parental_columns(record, 20, row)?
    } else {
        None
    };

    let variant = Variant {
        chrom: record.get(0).with_context(ctx)?.to_string(),
        pos: record.get(1).with_context(ctx)?.parse().with_context(ctx)?,
        ref_allele: record.get(2).with_context(ctx)?.to_string(),
        alt_allele: record.get(3).with_context(ctx)?.to_string(),
        high_ref_depth: record.get(4).with_context(ctx)?.parse().with_context(ctx)?,
        high_alt_depth: record.get(5).with_context(ctx)?.parse().with_context(ctx)?,
        high_dp: record.get(6).with_context(ctx)?.parse().with_context(ctx)?,
        high_gq: record.get(7).with_context(ctx)?.parse().with_context(ctx)?,
        low_ref_depth: record.get(8).with_context(ctx)?.parse().with_context(ctx)?,
        low_alt_depth: record.get(9).with_context(ctx)?.parse().with_context(ctx)?,
        low_dp: record.get(10).with_context(ctx)?.parse().with_context(ctx)?,
        low_gq: record.get(11).with_context(ctx)?.parse().with_context(ctx)?,
        parental_info,
    };

    let raw = GStatisticResult {
        variant,
        g_statistic: record.get(12).with_context(ctx)?.parse().with_context(ctx)?,
        snp_index_high: record.get(13).with_context(ctx)?.parse().with_context(ctx)?,
        snp_index_low: record.get(14).with_context(ctx)?.parse().with_context(ctx)?,
        delta_snp_index: record.get(15).with_context(ctx)?.parse().with_context(ctx)?,
    };

    Ok(SignificanceResult {
        raw,
        g_prime: record.get(16).with_context(ctx)?.parse().with_context(ctx)?,
        p_value: record.get(17).with_context(ctx)?.parse().with_context(ctx)?,
        q_value: record.get(18).with_context(ctx)?.parse().with_context(ctx)?,
        significant: record.get(19).with_context(ctx)?.parse::<bool>().with_context(ctx)?,
    })
}
