use anyhow::Result;
use bsa_gprime::{output, statistics, types::GStatisticResult, vcf_parser};
use clap::Parser;
use std::path::Path;

#[derive(Parser)]
#[command(name = "bsa-gprime")]
#[command(about = "Calculate G-statistic for BSA-Seq data", long_about = None)]
struct Args {
    /// Input VCF file (can be gzipped)
    #[arg(short, long)]
    input: String,

    /// Output CSV file
    #[arg(short, long)]
    output: String,

    /// Minimum read depth per sample
    #[arg(long, default_value = "10")]
    min_depth: u32,

    /// Minimum genotype quality
    #[arg(long, default_value = "20")]
    min_gq: u32,

    /// SNPs only (exclude INDELs)
    #[arg(long)]
    snps_only: bool,
}

fn main() -> Result<()> {
    let args = Args::parse();

    eprintln!("BSA G-Statistic Calculator (MVP version)");
    eprintln!("=========================================");
    eprintln!("Input VCF: {}", args.input);
    eprintln!("Output CSV: {}", args.output);
    eprintln!("Min depth: {}", args.min_depth);
    eprintln!("Min GQ: {}", args.min_gq);
    eprintln!("SNPs only: {}", args.snps_only);
    eprintln!();

    eprintln!("Step 1: Parsing VCF and applying filters...");
    let variants = vcf_parser::parse_vcf(
        Path::new(&args.input),
        args.min_depth,
        args.min_gq,
        args.snps_only,
    )?;

    if variants.is_empty() {
        eprintln!("ERROR: No variants passed filters!");
        return Ok(());
    }

    eprintln!();
    eprintln!("Step 2: Calculating G-statistics...");
    let results: Vec<_> = variants
        .iter()
        .enumerate()
        .map(|(i, v)| {
            if (i + 1) % 100_000 == 0 {
                eprintln!("  Calculated G-statistic for {} variants...", i + 1);
            }

            let g = statistics::calculate_g_statistic(
                v.resistant_ref_depth,
                v.resistant_alt_depth,
                v.susceptible_ref_depth,
                v.susceptible_alt_depth,
            );

            let snp_index_resistant =
                statistics::snp_index(v.resistant_ref_depth, v.resistant_alt_depth);
            let snp_index_susceptible =
                statistics::snp_index(v.susceptible_ref_depth, v.susceptible_alt_depth);
            let delta_snp_index = statistics::delta_snp_index(
                v.resistant_ref_depth,
                v.resistant_alt_depth,
                v.susceptible_ref_depth,
                v.susceptible_alt_depth,
            );

            GStatisticResult {
                variant: v.clone(),
                g_statistic: g,
                snp_index_resistant,
                snp_index_susceptible,
                delta_snp_index,
            }
        })
        .collect();

    eprintln!();
    eprintln!("Step 3: Calculating summary statistics...");
    let mut g_values: Vec<f64> = results.iter().map(|r| r.g_statistic).collect();
    g_values.sort_by(|a, b| a.partial_cmp(b).unwrap());

    let percentile = |p: f64| -> f64 {
        let idx = (p * g_values.len() as f64) as usize;
        g_values[idx.min(g_values.len() - 1)]
    };

    eprintln!("  Mean G: {:.3}", g_values.iter().sum::<f64>() / g_values.len() as f64);
    eprintln!("  Median G: {:.3}", percentile(0.50));
    eprintln!("  90th percentile: {:.3}", percentile(0.90));
    eprintln!("  95th percentile: {:.3}", percentile(0.95));
    eprintln!("  99th percentile: {:.3}", percentile(0.99));
    eprintln!("  99.9th percentile: {:.3}", percentile(0.999));
    eprintln!("  Max G: {:.3}", g_values.last().unwrap());

    eprintln!();
    eprintln!("Step 4: Writing results to CSV...");
    output::write_results(&results, Path::new(&args.output))?;

    eprintln!();
    eprintln!("Done! Results written to: {}", args.output);
    eprintln!();
    eprintln!("Next steps:");
    eprintln!("  - Visualize results using R/Python");
    eprintln!("  - Look for QTL peaks (high G-statistic regions)");
    eprintln!("  - Consider re-running with smoothing for clearer peaks");

    Ok(())
}
