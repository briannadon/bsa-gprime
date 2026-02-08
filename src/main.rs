use anyhow::Result;
use bsa_gprime::{csv_reader, output, peaks, significance, smoothing, statistics, types::GStatisticResult, vcf_parser};
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
    #[arg(short, long, required_unless_present = "plot_from")]
    input: Option<String>,

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

    /// Generate plots after analysis
    #[arg(long)]
    plot: bool,

    /// Generate plots from an existing results CSV (skips VCF analysis)
    #[arg(long, conflicts_with = "input")]
    plot_from: Option<String>,

    /// Output directory for plots (defaults to output file's directory)
    #[arg(long)]
    plot_dir: Option<String>,

    /// Plot output format: "png" (default) or "svg"
    #[arg(long, default_value = "png")]
    plot_format: String,

    /// Also generate delta SNP-index plot
    #[arg(long)]
    delta_snp_index: bool,

    /// Maximum number of chromosomes to plot
    #[arg(long, default_value = "12")]
    max_plot_chroms: usize,

    /// Write QTL peak summary to CSV file
    #[arg(long)]
    peaks_csv: Option<String>,
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

#[cfg(feature = "plotting")]
fn parse_plot_format(s: &str) -> Result<bsa_gprime::plotting::PlotFormat> {
    match s.to_lowercase().as_str() {
        "png" => Ok(bsa_gprime::plotting::PlotFormat::Png),
        "svg" => Ok(bsa_gprime::plotting::PlotFormat::Svg),
        other => anyhow::bail!("Invalid --plot-format '{}'. Must be 'png' or 'svg'", other),
    }
}

fn main() -> Result<()> {
    let args = Args::parse();

    // Configure rayon thread pool
    rayon::ThreadPoolBuilder::new()
        .num_threads(args.threads)
        .build_global()
        .unwrap();

    // ─── Path B: Plot-from-CSV mode ───
    if let Some(ref csv_path) = args.plot_from {
        return run_plot_from_csv(&args, csv_path);
    }

    // ─── Path A: Normal analysis run ───
    let input = args.input.as_ref().unwrap();

    // Validate input file exists
    if !Path::new(input).exists() {
        anyhow::bail!("Input file not found: {}", input);
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
            let input_path = Path::new(input);
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
            anyhow::bail!("Either --output or --output-dir must be specified (unless using --plot-from)");
        }
    };

    progress!(args.quiet, "BSA G-Statistic Calculator");
    progress!(args.quiet, "=========================================");
    progress!(args.quiet, "Input VCF: {}", input);
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
        Path::new(input),
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
        let g_prime_values = smoothing::smooth_g_statistics(&results, args.bandwidth, Some(&pb_smooth));
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

        // Peak detection
        let detected_peaks = peaks::identify_peaks(&sig_results, args.bandwidth);
        peaks::report_peaks(&detected_peaks);

        if let Some(ref peaks_path) = args.peaks_csv {
            peaks::write_peaks_csv(&detected_peaks, Path::new(peaks_path))?;
            progress!(args.quiet, "  Peak summary written to: {}", peaks_path);
        }

        // Plotting (gated on feature)
        #[cfg(feature = "plotting")]
        if args.plot {
            run_plots_sig(&args, &sig_results, &output_path)?;
        }
    } else {
        // Original pipeline (no significance testing)
        progress!(args.quiet);
        progress!(args.quiet, "Step 4: Writing results to CSV...");
        output::write_results(&results, Path::new(&output_path))?;

        // Peak detection by percentile
        let detected_peaks = peaks::identify_peaks_by_percentile(&results, 99.0, args.bandwidth);
        peaks::report_peaks(&detected_peaks);

        if let Some(ref peaks_path) = args.peaks_csv {
            peaks::write_peaks_csv(&detected_peaks, Path::new(peaks_path))?;
            progress!(args.quiet, "  Peak summary written to: {}", peaks_path);
        }

        // Plotting (gated on feature)
        #[cfg(feature = "plotting")]
        if args.plot {
            run_plots_basic(&args, &results, &output_path)?;
        }
    }

    progress!(args.quiet);
    progress!(args.quiet, "Done! Results written to: {}", output_path);

    Ok(())
}

/// Run plot-from-CSV mode: read existing CSV, detect peaks, generate plots.
fn run_plot_from_csv(args: &Args, csv_path: &str) -> Result<()> {
    progress!(args.quiet, "BSA G-Statistic Plotter (from CSV)");
    progress!(args.quiet, "=========================================");
    progress!(args.quiet, "Input CSV: {}", csv_path);
    progress!(args.quiet);

    let loaded = csv_reader::load_results_csv(Path::new(csv_path))?;

    match loaded {
        csv_reader::LoadedResults::WithSignificance(sig_results) => {
            // Peak detection
            let detected_peaks = peaks::identify_peaks(&sig_results, args.bandwidth);
            peaks::report_peaks(&detected_peaks);

            if let Some(ref peaks_path) = args.peaks_csv {
                peaks::write_peaks_csv(&detected_peaks, Path::new(peaks_path))?;
                progress!(args.quiet, "  Peak summary written to: {}", peaks_path);
            }

            #[cfg(feature = "plotting")]
            run_plots_sig(args, &sig_results, csv_path)?;

            #[cfg(not(feature = "plotting"))]
            if args.plot || args.delta_snp_index {
                eprintln!("Warning: plotting feature not enabled. Rebuild with default features to enable plots.");
            }
        }
        csv_reader::LoadedResults::Basic(results) => {
            // Peak detection by percentile
            let detected_peaks = peaks::identify_peaks_by_percentile(&results, 99.0, args.bandwidth);
            peaks::report_peaks(&detected_peaks);

            if let Some(ref peaks_path) = args.peaks_csv {
                peaks::write_peaks_csv(&detected_peaks, Path::new(peaks_path))?;
                progress!(args.quiet, "  Peak summary written to: {}", peaks_path);
            }

            #[cfg(feature = "plotting")]
            run_plots_basic(args, &results, csv_path)?;

            #[cfg(not(feature = "plotting"))]
            if args.plot || args.delta_snp_index {
                eprintln!("Warning: plotting feature not enabled. Rebuild with default features to enable plots.");
            }
        }
    }

    progress!(args.quiet);
    progress!(args.quiet, "Done!");

    Ok(())
}

/// Generate plots for significance results.
#[cfg(feature = "plotting")]
fn run_plots_sig(
    args: &Args,
    sig_results: &[bsa_gprime::types::SignificanceResult],
    reference_path: &str,
) -> Result<()> {
    use bsa_gprime::plotting;

    let format = parse_plot_format(&args.plot_format)?;
    let config = plotting::PlotConfig {
        width: 1800,
        row_height: 300,
        format,
        max_chromosomes: args.max_plot_chroms,
    };

    let (plot_dir, stem) = resolve_plot_paths(args, reference_path)?;
    std::fs::create_dir_all(&plot_dir)?;

    let ext = config.format.extension();

    progress!(args.quiet, "Generating plots...");

    // -log10(p-value) Manhattan plot
    let pval_path = plot_dir.join(format!("{}_neglog10p.{}", stem, ext));
    plotting::plot_neg_log10_pvalue(sig_results, &pval_path, &config)?;

    // G-prime plot
    let gprime_path = plot_dir.join(format!("{}_gprime.{}", stem, ext));
    plotting::plot_gprime_sig(sig_results, &gprime_path, &config, 99.0)?;

    // Delta SNP-index plot (optional)
    if args.delta_snp_index {
        let delta_path = plot_dir.join(format!("{}_delta_snp.{}", stem, ext));
        plotting::plot_delta_snp_index_sig(sig_results, &delta_path, &config)?;
    }

    Ok(())
}

/// Generate plots for basic results (no significance).
#[cfg(feature = "plotting")]
fn run_plots_basic(
    args: &Args,
    results: &[GStatisticResult],
    reference_path: &str,
) -> Result<()> {
    use bsa_gprime::plotting;

    let format = parse_plot_format(&args.plot_format)?;
    let config = plotting::PlotConfig {
        width: 1800,
        row_height: 300,
        format,
        max_chromosomes: args.max_plot_chroms,
    };

    let (plot_dir, stem) = resolve_plot_paths(args, reference_path)?;
    std::fs::create_dir_all(&plot_dir)?;

    let ext = config.format.extension();

    progress!(args.quiet, "Generating plots...");

    // G-stat plot
    let gstat_path = plot_dir.join(format!("{}_gstat.{}", stem, ext));
    plotting::plot_gstat(results, &gstat_path, &config, 99.0)?;

    // Delta SNP-index plot (optional)
    if args.delta_snp_index {
        let delta_path = plot_dir.join(format!("{}_delta_snp.{}", stem, ext));
        plotting::plot_delta_snp_index_basic(results, &delta_path, &config)?;
    }

    Ok(())
}

/// Resolve plot output directory and file stem from args and reference path.
#[cfg(feature = "plotting")]
fn resolve_plot_paths(
    args: &Args,
    reference_path: &str,
) -> Result<(std::path::PathBuf, String)> {
    let ref_path = Path::new(reference_path);

    let plot_dir = if let Some(ref dir) = args.plot_dir {
        std::path::PathBuf::from(dir)
    } else {
        ref_path
            .parent()
            .unwrap_or_else(|| Path::new("."))
            .to_path_buf()
    };

    let stem = ref_path
        .file_stem()
        .map(|s| s.to_string_lossy().to_string())
        .unwrap_or_else(|| "bsa_plot".to_string());

    Ok((plot_dir, stem))
}
