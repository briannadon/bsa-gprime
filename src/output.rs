use crate::types::{GStatisticResult, SignificanceResult};
use anyhow::Result;
use csv::Writer;
use std::io::Write;
use std::path::Path;

/// Check whether any result has parental info attached.
fn has_parental_info_basic(results: &[GStatisticResult]) -> bool {
    results.first().map_or(false, |r| r.variant.parental_info.is_some())
}

fn has_parental_info_sig(results: &[SignificanceResult]) -> bool {
    results.first().map_or(false, |r| r.raw.variant.parental_info.is_some())
}

pub fn write_results(results: &[GStatisticResult], path: &Path) -> Result<()> {
    let mut wtr = Writer::from_path(path)?;
    let parental = has_parental_info_basic(results);

    // Write header
    let mut headers = vec![
        "chrom",
        "pos",
        "ref",
        "alt",
        "high_ref_depth",
        "high_alt_depth",
        "high_dp",
        "high_gq",
        "low_ref_depth",
        "low_alt_depth",
        "low_dp",
        "low_gq",
        "g_statistic",
        "snp_index_high",
        "snp_index_low",
        "delta_snp_index",
    ];
    if parental {
        headers.extend_from_slice(&["parent_high_gt", "parent_low_gt", "alleles_swapped"]);
    }
    wtr.write_record(&headers)?;

    // Write data
    for result in results {
        let mut row = vec![
            result.variant.chrom.clone(),
            result.variant.pos.to_string(),
            result.variant.ref_allele.clone(),
            result.variant.alt_allele.clone(),
            result.variant.high_ref_depth.to_string(),
            result.variant.high_alt_depth.to_string(),
            result.variant.high_dp.to_string(),
            result.variant.high_gq.to_string(),
            result.variant.low_ref_depth.to_string(),
            result.variant.low_alt_depth.to_string(),
            result.variant.low_dp.to_string(),
            result.variant.low_gq.to_string(),
            format!("{:.6}", result.g_statistic),
            format!("{:.6}", result.snp_index_high),
            format!("{:.6}", result.snp_index_low),
            format!("{:.6}", result.delta_snp_index),
        ];
        if parental {
            if let Some(ref pi) = result.variant.parental_info {
                row.push(pi.parent_high_gt.clone());
                row.push(pi.parent_low_gt.clone());
                row.push(pi.alleles_swapped.to_string());
            }
        }
        wtr.write_record(&row)?;
    }

    wtr.flush()?;
    Ok(())
}

pub fn write_significance_results(results: &[SignificanceResult], path: &Path) -> Result<()> {
    let mut wtr = Writer::from_path(path)?;
    let parental = has_parental_info_sig(results);

    let mut headers = vec![
        "chrom",
        "pos",
        "ref",
        "alt",
        "high_ref_depth",
        "high_alt_depth",
        "high_dp",
        "high_gq",
        "low_ref_depth",
        "low_alt_depth",
        "low_dp",
        "low_gq",
        "g_statistic",
        "snp_index_high",
        "snp_index_low",
        "delta_snp_index",
        "g_prime",
        "p_value",
        "q_value",
        "significant",
    ];
    if parental {
        headers.extend_from_slice(&["parent_high_gt", "parent_low_gt", "alleles_swapped"]);
    }
    wtr.write_record(&headers)?;

    for result in results {
        let mut row = vec![
            result.raw.variant.chrom.clone(),
            result.raw.variant.pos.to_string(),
            result.raw.variant.ref_allele.clone(),
            result.raw.variant.alt_allele.clone(),
            result.raw.variant.high_ref_depth.to_string(),
            result.raw.variant.high_alt_depth.to_string(),
            result.raw.variant.high_dp.to_string(),
            result.raw.variant.high_gq.to_string(),
            result.raw.variant.low_ref_depth.to_string(),
            result.raw.variant.low_alt_depth.to_string(),
            result.raw.variant.low_dp.to_string(),
            result.raw.variant.low_gq.to_string(),
            format!("{:.6}", result.raw.g_statistic),
            format!("{:.6}", result.raw.snp_index_high),
            format!("{:.6}", result.raw.snp_index_low),
            format!("{:.6}", result.raw.delta_snp_index),
            format!("{:.6}", result.g_prime),
            format!("{:.6e}", result.p_value),
            format!("{:.6e}", result.q_value),
            result.significant.to_string(),
        ];
        if parental {
            if let Some(ref pi) = result.raw.variant.parental_info {
                row.push(pi.parent_high_gt.clone());
                row.push(pi.parent_low_gt.clone());
                row.push(pi.alleles_swapped.to_string());
            }
        }
        wtr.write_record(&row)?;
    }

    wtr.flush()?;
    Ok(())
}

/// Write significant regions as a BED file with merged intervals.
///
/// Adjacent significant SNPs on the same chromosome within `merge_distance` bp
/// are merged into contiguous regions.
///
/// BED format: chrom, start (0-based), end (exclusive), name, score, strand
pub fn write_bed(
    results: &[SignificanceResult],
    path: &Path,
    merge_distance: u64,
) -> Result<()> {
    let mut file = std::fs::File::create(path)?;

    let mut region_num = 0u64;
    // (chrom, start_pos, end_pos, max_gprime, snp_count)
    let mut current: Option<(String, u64, u64, f64, u64)> = None;

    let emit_region =
        |file: &mut std::fs::File, chrom: &str, start: u64, end: u64, max_gp: f64, num: u64| -> Result<()> {
            let score = (max_gp as u64).min(1000);
            writeln!(
                file,
                "{}\t{}\t{}\tQTL_region_{}\t{}\t.",
                chrom,
                start.saturating_sub(1), // convert 1-based pos to 0-based BED start
                end,                      // 1-based pos works as 0-based exclusive end
                num,
                score,
            )?;
            Ok(())
        };

    for r in results {
        if !r.significant {
            continue;
        }

        let pos = r.raw.variant.pos;
        let chrom = &r.raw.variant.chrom;

        match current.take() {
            Some((cur_chrom, start, end, max_gp, count))
                if *chrom == cur_chrom && pos <= end + merge_distance =>
            {
                // Extend current region
                current = Some((cur_chrom, start, pos, max_gp.max(r.g_prime), count + 1));
            }
            Some((cur_chrom, start, end, max_gp, _count)) => {
                // Emit previous region, start new one
                region_num += 1;
                emit_region(&mut file, &cur_chrom, start, end, max_gp, region_num)?;
                current = Some((chrom.clone(), pos, pos, r.g_prime, 1));
            }
            None => {
                current = Some((chrom.clone(), pos, pos, r.g_prime, 1));
            }
        }
    }

    // Emit final region
    if let Some((chrom, start, end, max_gp, _count)) = current {
        region_num += 1;
        emit_region(&mut file, &chrom, start, end, max_gp, region_num)?;
    }

    Ok(())
}
