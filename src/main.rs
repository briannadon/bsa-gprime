use anyhow::Result;
use bsa_gprime::{output, significance, smoothing, statistics, types::GStatisticResult, vcf_parser};
use clap::Parser;
use indicatif::{ProgressBar, ProgressStyle};
use rayon::prelude::*;
use std::path::Path;

#[derive(Parser)]
#[command(name = "bsa-gprime")]
#[command(version)]
#[command(about = "Calculate G-statistic for BSA-Seq data", long_about = None)]
struct Args {
    /// Input VCF file (can be gzipped)
    #[arg(short, long)]
    input: String,

    /// Output CSV file path
    #[arg(short, long)]
    output: Option<String>,

    /// Output directory (file auto-named from input)
    #[arg(long)]
    output_dir: Option<String>,

    /// Minimum read depth per sample
    #[arg(long, default_value = "10")]
    min_depth: u32,

    /// Minimum genotype quality
    #[arg(long, default_value = "20")]
    min_gq: u32,

    /// SNPs only (exclude INDELs)
    #[arg(long)]
    snps_only: bool,

    /// Suppress progress output
    #[arg(short, long)]
    quiet: bool,

    /// Disable significance testing pipeline (smoothing, p-values, FDR)
    #[arg(long)]
    no_significance: bool,

    /// Tricube kernel bandwidth in base pairs
    #[arg(long, default_value = "1000000")]
    bandwidth: u64,

    /// Null distribution method: "parametric" or "nonparametric"
    #[arg(long, default_value = "nonparametric")]
    null_method: String,

    /// Number of individuals per bulk (required for parametric method)
    #[arg(long)]
    bulk_size: Option<f64>,

    /// Ploidy of the organism (used with --bulk-size for parametric method)
    #[arg(long, default_value = "2")]
    ploidy: u32,

    /// Recombination rate in cM/Mb (required for full equation 12 covariance term
    /// in parametric mode). If omitted, only the diagonal Σk_j² term is used.
    #[arg(long)]
    recombination_rate: Option<f64>,

    /// FDR threshold for significance calls
    #[arg(long, default_value = "0.05")]
    fdr_threshold: f64,

    /// Number of threads for parallel processing
    #[arg(long, default_value_t = num_cpus())]
    threads: usize,

    /// Output BED file of significant regions (requires significance testing)
    #[arg(long)]
    bed: Option<String>,
}

fn num_cpus() -> usize {
    std::thread::available_parallelism()
        .map(|n| n.get())
        .unwrap_or(1)
}

macro_rules! progress {
    ($quiet:expr) => {
        if !$quiet {
            eprintln!();
        }
    };
    ($quiet:expr, $($arg:tt)*) => {
        if !$quiet {
            eprintln!($($arg)*);
        }
    };
}

fn make_progress_bar(quiet: bool, len: u64) -> ProgressBar {
    if quiet {
        return ProgressBar::hidden();
    }
    let pb = ProgressBar::new(len);
    pb.set_style(
        ProgressStyle::with_template("  [{elapsed_precise}/{eta_precise}] {bar:40} {pos}/{len} {msg}")
            .unwrap()
            .progress_chars("=> "),
    );
    pb
}

fn make_spinner(quiet: bool) -> ProgressBar {
    if quiet {
        return ProgressBar::hidden();
    }
    let pb = ProgressBar::new_spinner();
    pb.set_style(
        ProgressStyle::with_template("  {spinner} [{elapsed_precise}] {pos} {msg}")
            .unwrap(),
    );
    pb
}

fn main() -> Result<()> {
    let args = Args::parse();

    // Configure rayon thread pool
    rayon::ThreadPoolBuilder::new()
        .num_threads(args.threads)
        .build_global()
        .unwrap();

    // Validate input file exists
    if !Path::new(&args.input).exists() {
        anyhow::bail!("Input file not found: {}", args.input);
    }

    // Validate --bed requires significance testing
    if args.bed.is_some() && args.no_significance {
        anyhow::bail!("--bed requires significance testing; cannot be used with --no-significance");
    }

    // Validate significance args
    if !args.no_significance && args.null_method == "parametric" && args.bulk_size.is_none() {
        anyhow::bail!("--bulk-size is required when using --null-method parametric");
    }

    if !["parametric", "nonparametric"].contains(&args.null_method.as_str()) {
        anyhow::bail!(
            "Invalid --null-method '{}'. Must be 'parametric' or 'nonparametric'",
            args.null_method
        );
    }

    // Resolve output path
    let output_path = match (&args.output, &args.output_dir) {
        (Some(output), _) => output.clone(),
        (None, Some(dir)) => {
            let input_path = Path::new(&args.input);
            let stem = input_path
                .file_stem()
                .and_then(|s| {
                    // Strip .vcf from .vcf.gz
                    let s_str = s.to_string_lossy();
                    if s_str.ends_with(".vcf") {
                        Some(s_str.trim_end_matches(".vcf").to_string())
                    } else {
                        Some(s_str.to_string())
                    }
                })
                .unwrap_or_else(|| "output".to_string());
            let dir_path = Path::new(dir);
            std::fs::create_dir_all(dir_path)?;
            dir_path
                .join(format!("{}_bsa_results.csv", stem))
                .to_string_lossy()
                .to_string()
        }
        (None, None) => {
            anyhow::bail!("Either --output or --output-dir must be specified");
        }
    };

    progress!(args.quiet, "BSA G-Statistic Calculator");
    progress!(args.quiet, "=========================================");
    progress!(args.quiet, "Input VCF: {}", args.input);
    progress!(args.quiet, "Output CSV: {}", output_path);
    progress!(args.quiet, "Min depth: {}", args.min_depth);
    progress!(args.quiet, "Min GQ: {}", args.min_gq);
    progress!(args.quiet, "SNPs only: {}", args.snps_only);
    progress!(args.quiet, "Threads: {}", args.threads);
    if args.no_significance {
        progress!(args.quiet, "Significance testing: disabled");
    } else {
        progress!(args.quiet, "Significance testing: enabled");
        progress!(args.quiet, "  Bandwidth: {} bp", args.bandwidth);
        progress!(args.quiet, "  Null method: {}", args.null_method);
        if args.null_method == "parametric" {
            let bulk_size = args.bulk_size.unwrap();
            let n_s = bulk_size * args.ploidy as f64;
            progress!(args.quiet, "  Bulk size: {} individuals", bulk_size);
            progress!(args.quiet, "  Ploidy: {}", args.ploidy);
            progress!(args.quiet, "  Effective n_s: {}", n_s);
        }
        progress!(args.quiet, "  FDR threshold: {}", args.fdr_threshold);
    }
    if let Some(ref bed_path) = args.bed {
        progress!(args.quiet, "BED output: {}", bed_path);
    }
    progress!(args.quiet);

    // Step 1: Parse VCF
    progress!(args.quiet, "Step 1: Parsing VCF and applying filters...");
    let pb_vcf = make_spinner(args.quiet);
    pb_vcf.set_message("variants parsed");
    let variants = vcf_parser::parse_vcf(
        Path::new(&args.input),
        args.min_depth,
        args.min_gq,
        args.snps_only,
        args.quiet,
    )?;
    pb_vcf.finish_and_clear();

    if variants.is_empty() {
        anyhow::bail!("No variants passed filters!");
    }

    // Step 2: Calculate G-statistics (parallelized)
    progress!(args.quiet);
    progress!(args.quiet, "Step 2: Calculating G-statistics...");
    let pb_gstat = make_progress_bar(args.quiet, variants.len() as u64);
    let results: Vec<_> = variants
        .par_iter()
        .map(|v| {
            pb_gstat.inc(1);

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
    pb_gstat.finish_and_clear();

    // Step 3: Summary statistics
    progress!(args.quiet);
    progress!(args.quiet, "Step 3: Calculating summary statistics...");
    let mut g_values: Vec<f64> = results.iter().map(|r| r.g_statistic).collect();
    g_values.sort_by(|a, b| a.partial_cmp(b).unwrap());

    let percentile = |p: f64| -> f64 {
        let idx = (p * g_values.len() as f64) as usize;
        g_values[idx.min(g_values.len() - 1)]
    };

    progress!(args.quiet, "  Mean G: {:.3}", g_values.iter().sum::<f64>() / g_values.len() as f64);
    progress!(args.quiet, "  Median G: {:.3}", percentile(0.50));
    progress!(args.quiet, "  90th percentile: {:.3}", percentile(0.90));
    progress!(args.quiet, "  95th percentile: {:.3}", percentile(0.95));
    progress!(args.quiet, "  99th percentile: {:.3}", percentile(0.99));
    progress!(args.quiet, "  99.9th percentile: {:.3}", percentile(0.999));
    progress!(args.quiet, "  Max G: {:.3}", g_values.last().unwrap());

    if !args.no_significance {
        // Significance testing pipeline
        progress!(args.quiet);
        progress!(args.quiet, "Step 4: Smoothing G-statistics (tricube kernel, bandwidth={} bp)...", args.bandwidth);
        let pb_smooth = make_progress_bar(args.quiet, results.len() as u64);
        let g_prime_values = smoothing::smooth_g_statistics(&results, args.bandwidth);
        pb_smooth.finish_and_clear();

        progress!(args.quiet, "Step 5: Estimating null distribution ({})...", args.null_method);
        let null_params = match args.null_method.as_str() {
            "parametric" => {
                let n_s = args.bulk_size.unwrap() * args.ploidy as f64;
                let avg_coverage = significance::compute_avg_coverage(&results);
                progress!(args.quiet, "  Average per-bulk coverage (C): {:.1}", avg_coverage);
                
                // Compute smoothing factors for equation 12 (Var[G'] estimation)
                let factors = smoothing::compute_smoothing_factors(
                    &results, args.bandwidth, args.recombination_rate,
                );
                progress!(args.quiet, "  Smoothing factor (Σkⱼ²): {:.4}", factors.sum_k_squared);
                if args.recombination_rate.is_some() {
                    progress!(args.quiet, "  Recombination rate: {:.2} cM/Mb", args.recombination_rate.unwrap());
                    progress!(args.quiet, "  Covariance weight: {:.4}", factors.covariance_weight);
                } else {
                    progress!(args.quiet, "  Covariance term: omitted (no --recombination-rate provided)");
                }
                
                significance::estimate_null_parametric_gprime(
                    n_s,
                    avg_coverage,
                    Some(factors.sum_k_squared),
                    Some(factors.covariance_weight),
                )
            }
            "nonparametric" => significance::estimate_null_robust(&g_prime_values),
            _ => unreachable!(),
        };
        progress!(args.quiet, "  Null distribution: mu={:.4}, sigma={:.4}", null_params.mu, null_params.sigma);

        progress!(args.quiet, "Step 6: Computing p-values and FDR correction...");
        let pb_pval = make_spinner(args.quiet);
        pb_pval.set_message("computing p-values");
        let sig_results = significance::run_significance_pipeline(
            results,
            &g_prime_values,
            &null_params,
            args.fdr_threshold,
        );
        pb_pval.finish_and_clear();

        let n_significant = sig_results.iter().filter(|r| r.significant).count();
        progress!(args.quiet, "  Significant SNPs (q < {}): {} / {}",
            args.fdr_threshold, n_significant, sig_results.len());

        progress!(args.quiet);
        progress!(args.quiet, "Step 7: Writing results to CSV...");
        output::write_significance_results(&sig_results, Path::new(&output_path))?;

        // Write BED file if requested
        if let Some(ref bed_path) = args.bed {
            progress!(args.quiet, "Step 8: Writing BED file of significant regions...");
            output::write_bed(&sig_results, Path::new(bed_path), args.bandwidth)?;
            let bed_regions = std::fs::read_to_string(bed_path)?
                .lines()
                .count();
            progress!(args.quiet, "  Wrote {} regions to {}", bed_regions, bed_path);
        }
    } else {
        // Original pipeline (no significance testing)
        progress!(args.quiet);
        progress!(args.quiet, "Step 4: Writing results to CSV...");
        output::write_results(&results, Path::new(&output_path))?;
    }

    progress!(args.quiet);
    progress!(args.quiet, "Done! Results written to: {}", output_path);
    progress!(args.quiet);
    progress!(args.quiet, "Next steps:");
    progress!(args.quiet, "  - Visualize results using R/Python");
    progress!(args.quiet, "  - Look for QTL peaks (high G-statistic regions)");
    if args.no_significance {
        progress!(args.quiet, "  - Re-run without --no-significance for p-values and FDR correction");
    }

    Ok(())
}
