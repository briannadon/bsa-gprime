#!/usr/bin/env python3
"""
Visualize BSA G-prime results to identify QTL peaks.

When significance testing columns are present (p_value, q_value, significant),
the default plot shows -log10(p-value) Manhattan-style peaks. Falls back to
G-statistic / G' when significance columns are absent (--no-significance runs).
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import argparse
from pathlib import Path


def load_results(csv_path):
    """Load BSA results from CSV file."""
    df = pd.read_csv(csv_path)
    has_sig = "p_value" in df.columns
    print(f"Loaded {len(df):,} variants from {csv_path}")
    print(f"Chromosomes: {df['chrom'].nunique()}")
    if has_sig:
        n_sig = df["significant"].sum() if "significant" in df.columns else 0
        print(f"Significant variants: {n_sig:,}")
    else:
        print("No significance columns found (--no-significance run)")
    return df


def apply_smoothing(df, window_bp=2_000_000):
    """
    Apply simple moving average smoothing to G-statistics.

    This is a basic smoothing - the Rust implementation uses a tricube kernel.
    Only useful for --no-significance runs that lack a g_prime column.
    """
    print(f"Applying rolling window smoothing (window={window_bp:,} bp)...")

    smoothed = []
    for chrom in df["chrom"].unique():
        chrom_df = df[df["chrom"] == chrom].copy().sort_values("pos")

        chrom_df["pos_mb"] = chrom_df["pos"] / 1_000_000

        avg_density = len(chrom_df) / max(
            chrom_df["pos"].max() - chrom_df["pos"].min(), 1
        )
        window_snps = max(int(window_bp * avg_density), 1)

        chrom_df["g_prime"] = (
            chrom_df["g_statistic"]
            .rolling(window=window_snps, center=True, min_periods=1)
            .mean()
        )

        smoothed.append(chrom_df)

    return pd.concat(smoothed, ignore_index=True)


def plot_log_pvalue_manhattan(
    df, output_path=None, fdr_threshold=0.05, max_chromosomes=12
):
    """
    Create Manhattan-style plot of -log10(p-value) with facets per chromosome.

    Significant variants (q_value < fdr_threshold) are highlighted.
    A horizontal line is drawn at the -log10(p) value corresponding to the
    most liberal significant call (i.e. the largest p-value still called
    significant after FDR correction).
    """
    print("Plotting -log10(p-value) Manhattan plot...")

    # Compute -log10(p_value), clamping p=0 to a large finite value
    neg_log_p = -np.log10(df["p_value"].clip(lower=1e-300))

    chromosomes = sorted(df["chrom"].unique())[:max_chromosomes]
    n_chroms = len(chromosomes)

    print(f"Chromosomes to plot ({n_chroms}): {chromosomes}")

    n_cols = 3
    n_rows = int(np.ceil(n_chroms / n_cols))

    fig, axes = plt.subplots(
        n_rows, n_cols, figsize=(18, 3 * n_rows), sharex=False, sharey=True
    )
    axes = np.atleast_1d(axes).flatten()

    # Determine significance threshold line from FDR-corrected calls
    has_sig_col = "significant" in df.columns
    if has_sig_col:
        sig_mask = df["significant"].astype(bool)
    else:
        sig_mask = df["q_value"] < fdr_threshold

    # Threshold line: -log10(p) at the boundary of significance
    if sig_mask.any():
        # The largest p-value among significant calls gives the visual threshold
        boundary_p = df.loc[sig_mask, "p_value"].max()
        threshold_line = -np.log10(max(boundary_p, 1e-300))
        print(f"Significance threshold line at -log10(p) = {threshold_line:.2f}")
    else:
        threshold_line = None
        print("No significant variants found; no threshold line drawn")

    for idx, chrom in enumerate(chromosomes):
        ax = axes[idx]
        mask = df["chrom"] == chrom
        chrom_df = df[mask].sort_values("pos")
        chrom_neg_log_p = neg_log_p[mask].loc[chrom_df.index]

        pos_mb = chrom_df["pos"] / 1_000_000

        # Background points
        ax.scatter(
            pos_mb, chrom_neg_log_p, s=2, alpha=0.4, color="#2E86AB", linewidths=0
        )

        # Highlight significant points
        chrom_sig = sig_mask[mask].loc[chrom_df.index]
        if chrom_sig.any():
            ax.scatter(
                pos_mb[chrom_sig],
                chrom_neg_log_p[chrom_sig],
                s=8,
                alpha=0.8,
                color="#A23B72",
                linewidths=0,
                zorder=5,
                label="Significant",
            )

        # Threshold line
        if threshold_line is not None:
            ax.axhline(
                y=threshold_line,
                color="red",
                linestyle="--",
                linewidth=1,
                alpha=0.5,
            )

        ax.set_title(f"{chrom}", fontsize=10, fontweight="bold")
        ax.set_xlabel("Position (Mb)", fontsize=9)
        ax.grid(True, alpha=0.3, linestyle=":")

        if idx % n_cols == 0:
            ax.set_ylabel("$-\\log_{10}(p)$", fontsize=9)

    for idx in range(n_chroms, len(axes)):
        axes[idx].set_visible(False)

    fig.suptitle(
        "BSA QTL Mapping - $-\\log_{10}(p\\text{-value})$ across Genome",
        fontsize=14,
        fontweight="bold",
    )

    if n_chroms > 0 and sig_mask.any():
        axes[0].legend(fontsize=8, loc="upper right")

    plt.tight_layout(rect=[0, 0, 1, 0.98])

    if output_path:
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
        print(f"Plot saved to: {output_path}")
    else:
        plt.show()

    return fig


def plot_qtl_manhattan(
    df,
    output_path=None,
    threshold_percentile=99,
    use_smoothed=False,
    max_chromosomes=12,
):
    """
    Create Manhattan-style plot of G-statistic / G' with facets per chromosome.

    Used as fallback when significance columns are not available.
    """
    stat_col = (
        "g_prime" if use_smoothed and "g_prime" in df.columns else "g_statistic"
    )
    stat_label = "G' (smoothed)" if stat_col == "g_prime" else "G-statistic"

    print(f"Plotting {stat_label} (no p-values available)...")

    threshold = df[stat_col].quantile(threshold_percentile / 100)
    print(f"{threshold_percentile}th percentile threshold: {threshold:.3f}")

    chromosomes = sorted(df["chrom"].unique())[:max_chromosomes]
    n_chroms = len(chromosomes)

    print(f"Chromosomes to plot ({n_chroms}): {chromosomes}")

    n_cols = 3
    n_rows = int(np.ceil(n_chroms / n_cols))

    fig, axes = plt.subplots(
        n_rows, n_cols, figsize=(18, 3 * n_rows), sharex=False, sharey=True
    )
    axes = np.atleast_1d(axes).flatten()

    for idx, chrom in enumerate(chromosomes):
        ax = axes[idx]
        chrom_df = df[df["chrom"] == chrom].sort_values("pos")

        pos_mb = chrom_df["pos"] / 1_000_000
        g_values = chrom_df[stat_col]

        ax.plot(pos_mb, g_values, linewidth=0.5, alpha=0.7, color="#2E86AB")

        significant = g_values > threshold
        if significant.any():
            ax.scatter(
                pos_mb[significant],
                g_values[significant],
                color="#A23B72",
                s=10,
                alpha=0.8,
                zorder=5,
                label=f">{threshold_percentile}th percentile",
            )

        ax.axhline(
            y=threshold,
            color="red",
            linestyle="--",
            linewidth=1,
            alpha=0.5,
            label=f"{threshold_percentile}th percentile",
        )

        ax.set_title(f"{chrom}", fontsize=10, fontweight="bold")
        ax.set_xlabel("Position (Mb)", fontsize=9)
        ax.grid(True, alpha=0.3, linestyle=":")

        if idx % n_cols == 0:
            ax.set_ylabel(stat_label, fontsize=9)

    for idx in range(n_chroms, len(axes)):
        axes[idx].set_visible(False)

    fig.suptitle(
        f"BSA QTL Mapping - {stat_label} across Genome",
        fontsize=14,
        fontweight="bold",
    )

    if n_chroms > 0:
        axes[0].legend(fontsize=8, loc="upper right")

    plt.tight_layout(rect=[0, 0, 1, 0.98])

    if output_path:
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
        print(f"Plot saved to: {output_path}")
    else:
        plt.show()

    return fig


def plot_delta_snp_index(df, output_path=None, max_chromosomes=12):
    """
    Plot delta SNP index across chromosomes.

    Delta SNP index shows the difference in allele frequencies between bulks.
    """
    print("Plotting delta SNP index...")

    chromosomes = sorted(df["chrom"].unique())[:max_chromosomes]
    n_chroms = len(chromosomes)

    print(f"Chromosomes to plot ({n_chroms}): {chromosomes}")

    n_cols = 3
    n_rows = int(np.ceil(n_chroms / n_cols))

    fig, axes = plt.subplots(
        n_rows, n_cols, figsize=(18, 3 * n_rows), sharex=False, sharey=True
    )
    axes = np.atleast_1d(axes).flatten()

    for idx, chrom in enumerate(chromosomes):
        ax = axes[idx]
        chrom_df = df[df["chrom"] == chrom].sort_values("pos")

        pos_mb = chrom_df["pos"] / 1_000_000
        delta = chrom_df["delta_snp_index"]

        ax.plot(pos_mb, delta, linewidth=0.5, alpha=0.7, color="#F18F01")
        ax.axhline(y=0, color="black", linestyle="-", linewidth=0.5, alpha=0.3)

        ax.set_title(f"{chrom}", fontsize=10, fontweight="bold")
        ax.set_xlabel("Position (Mb)", fontsize=9)
        ax.set_ylim(-1.1, 1.1)
        ax.grid(True, alpha=0.3, linestyle=":")

        if idx % n_cols == 0:
            ax.set_ylabel("$\\Delta$(SNP-index)", fontsize=9)

    for idx in range(n_chroms, len(axes)):
        axes[idx].set_visible(False)

    fig.suptitle(
        "BSA - Delta SNP Index across Genome", fontsize=14, fontweight="bold"
    )

    plt.tight_layout(rect=[0, 0, 1, 0.98])

    if output_path:
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
        print(f"Plot saved to: {output_path}")
    else:
        plt.show()

    return fig


def identify_qtl_peaks(df, fdr_threshold=0.05, min_peak_width_bp=500_000):
    """
    Identify and report QTL peak regions.

    When p_value column is present, peaks are defined as contiguous runs of
    significant variants. Otherwise falls back to percentile-based thresholding
    on G-statistic / G'.
    """
    has_pval = "p_value" in df.columns

    if has_pval:
        # Use FDR-corrected significance calls
        if "significant" in df.columns:
            above = df["significant"].astype(bool)
        else:
            above = df["q_value"] < fdr_threshold
        score_col = "p_value"
        score_label = "-log10(p)"
    else:
        stat_col = "g_prime" if "g_prime" in df.columns else "g_statistic"
        threshold = df[stat_col].quantile(0.99)
        above = df[stat_col] > threshold
        score_col = stat_col
        score_label = stat_col

    print(f"\n{'='*70}")
    if has_pval:
        print(f"QTL Peak Detection (FDR < {fdr_threshold})")
    else:
        print(f"QTL Peak Detection ({score_label} > 99th percentile)")
    print(f"{'='*70}")

    peaks = []
    for chrom in sorted(df["chrom"].unique()):
        chrom_df = df[df["chrom"] == chrom].sort_values("pos")
        significant = chrom_df[above[chrom_df.index]]

        if len(significant) == 0:
            continue

        # Group nearby significant SNPs into peaks
        significant = significant.copy()
        significant["gap"] = significant["pos"].diff()
        significant["new_peak"] = significant["gap"] > min_peak_width_bp
        significant["peak_id"] = significant["new_peak"].cumsum()

        for peak_id in significant["peak_id"].unique():
            peak_snps = significant[significant["peak_id"] == peak_id]

            if has_pval:
                min_p = peak_snps["p_value"].min()
                peak_info = {
                    "chrom": chrom,
                    "start": peak_snps["pos"].min(),
                    "end": peak_snps["pos"].max(),
                    "width_mb": (peak_snps["pos"].max() - peak_snps["pos"].min())
                    / 1_000_000,
                    "n_snps": len(peak_snps),
                    "max_neg_log10_p": -np.log10(max(min_p, 1e-300)),
                    "max_g_prime": peak_snps["g_prime"].max()
                    if "g_prime" in peak_snps.columns
                    else np.nan,
                }
            else:
                peak_info = {
                    "chrom": chrom,
                    "start": peak_snps["pos"].min(),
                    "end": peak_snps["pos"].max(),
                    "width_mb": (peak_snps["pos"].max() - peak_snps["pos"].min())
                    / 1_000_000,
                    "n_snps": len(peak_snps),
                    "max_g": peak_snps[score_col].max(),
                    "mean_g": peak_snps[score_col].mean(),
                }
            peaks.append(peak_info)

    sort_col = "max_neg_log10_p" if has_pval else "max_g"
    peaks_df = pd.DataFrame(peaks)
    if not peaks_df.empty:
        peaks_df = peaks_df.sort_values(sort_col, ascending=False)

    print(f"\nFound {len(peaks_df)} QTL peaks:\n")
    if not peaks_df.empty:
        print(peaks_df.to_string(index=False))
    print(f"\n{'='*70}\n")

    return peaks_df


def main():
    parser = argparse.ArgumentParser(
        description="Visualize BSA G-prime results and identify QTL peaks"
    )
    parser.add_argument("input", help="Input CSV file from bsa-gprime")
    parser.add_argument("-o", "--output", help="Output plot file (PNG/PDF)")
    parser.add_argument(
        "--smooth",
        action="store_true",
        help="Apply smoothing to raw G-statistics (for --no-significance runs)",
    )
    parser.add_argument(
        "--window",
        type=int,
        default=2_000_000,
        help="Smoothing window size in bp (default: 2,000,000)",
    )
    parser.add_argument(
        "--fdr-threshold",
        type=float,
        default=0.05,
        help="FDR threshold for significance (default: 0.05)",
    )
    parser.add_argument(
        "--threshold",
        type=float,
        default=99,
        help="Percentile threshold (used only when p-values are absent; default: 99)",
    )
    parser.add_argument(
        "--delta-snp-index",
        action="store_true",
        help="Also plot delta SNP index",
    )
    parser.add_argument(
        "--max-chroms",
        type=int,
        default=12,
        help="Maximum number of chromosomes to plot (default: 12)",
    )

    args = parser.parse_args()

    # Load data
    df = load_results(args.input)

    has_pval = "p_value" in df.columns

    # Apply smoothing if requested (only useful without significance columns)
    if args.smooth and not has_pval:
        df = apply_smoothing(df, window_bp=args.window)
    elif args.smooth and has_pval:
        print(
            "Note: --smooth ignored; data already has significance columns from Rust pipeline"
        )

    # Determine output paths
    if args.output:
        output_base = Path(args.output).stem
        output_dir = Path(args.output).parent
        output_ext = Path(args.output).suffix
    else:
        output_base = "qtl_plot"
        output_dir = Path(".")
        output_ext = ".png"

    # Create main plot
    if has_pval:
        # -log10(p-value) Manhattan plot
        main_output = (
            output_dir / f"{output_base}_neglog10p{output_ext}"
            if args.output
            else None
        )
        plot_log_pvalue_manhattan(
            df,
            output_path=main_output,
            fdr_threshold=args.fdr_threshold,
            max_chromosomes=args.max_chroms,
        )
    else:
        # Fall back to G-statistic / G' plot
        use_smoothed = "g_prime" in df.columns
        g_output = (
            output_dir / f"{output_base}_gstat{output_ext}" if args.output else None
        )
        plot_qtl_manhattan(
            df,
            output_path=g_output,
            threshold_percentile=args.threshold,
            use_smoothed=use_smoothed,
            max_chromosomes=args.max_chroms,
        )

    # Create delta SNP index plot if requested
    if args.delta_snp_index:
        delta_output = (
            output_dir / f"{output_base}_delta{output_ext}" if args.output else None
        )
        plot_delta_snp_index(
            df, output_path=delta_output, max_chromosomes=args.max_chroms
        )

    # Identify and report peaks
    peaks = identify_qtl_peaks(
        df,
        fdr_threshold=args.fdr_threshold,
        min_peak_width_bp=args.window,
    )

    # Save peaks to CSV
    if args.output and not peaks.empty:
        peaks_csv = output_dir / f"{output_base}_peaks.csv"
        peaks.to_csv(peaks_csv, index=False)
        print(f"Peak summary saved to: {peaks_csv}")

    if not args.output:
        plt.show()


if __name__ == "__main__":
    main()
