#!/usr/bin/env python3
"""
Visualize BSA G-statistic results to identify QTL peaks.

This script creates faceted plots showing G-statistics across all chromosomes,
making it easy to identify potential QTL regions.
"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import argparse
import numpy as np
from pathlib import Path


def load_results(csv_path):
    """Load BSA results from CSV file."""
    df = pd.read_csv(csv_path)
    print(f"Loaded {len(df):,} variants from {csv_path}")
    print(f"Chromosomes: {df['chrom'].nunique()}")
    return df


def apply_smoothing(df, window_bp=2_000_000):
    """
    Apply simple moving average smoothing to G-statistics.

    This is a basic smoothing - the Rust implementation would use tricube kernel.
    """
    print(f"Applying rolling window smoothing (window={window_bp:,} bp)...")

    smoothed = []
    for chrom in df['chrom'].unique():
        chrom_df = df[df['chrom'] == chrom].copy().sort_values('pos')

        # For simplicity, use a position-based window
        # Convert to Mb for easier handling
        chrom_df['pos_mb'] = chrom_df['pos'] / 1_000_000

        # Use rolling window based on number of SNPs (approximate)
        # Estimate window size as number of SNPs in the window
        avg_density = len(chrom_df) / (chrom_df['pos'].max() - chrom_df['pos'].min())
        window_snps = int(window_bp * avg_density)
        window_snps = max(window_snps, 1)  # At least 1

        chrom_df['g_prime'] = chrom_df['g_statistic'].rolling(
            window=window_snps, center=True, min_periods=1
        ).mean()

        smoothed.append(chrom_df)

    return pd.concat(smoothed, ignore_index=True)


def plot_qtl_manhattan(df, output_path=None, threshold_percentile=99,
                       use_smoothed=False, max_chromosomes=12):
    """
    Create Manhattan-style plot with facets for each chromosome.

    Parameters:
    -----------
    df : DataFrame
        Results dataframe with columns: chrom, pos, g_statistic, g_prime (optional)
    output_path : str, optional
        Path to save the plot
    threshold_percentile : float
        Percentile for drawing significance threshold line
    use_smoothed : bool
        Use g_prime instead of g_statistic if available
    max_chromosomes : int
        Maximum number of chromosomes to display
    """
    # Determine which statistic to plot
    stat_col = 'g_prime' if use_smoothed and 'g_prime' in df.columns else 'g_statistic'
    stat_label = "G' (smoothed)" if stat_col == 'g_prime' else "G-statistic"

    print(f"Plotting {stat_label}...")

    # Calculate threshold
    threshold = df[stat_col].quantile(threshold_percentile / 100)
    print(f"{threshold_percentile}th percentile threshold: {threshold:.3f}")

    # Get chromosome list (sorted)
    chromosomes = sorted(df['chrom'].unique())[:max_chromosomes]
    n_chroms = len(chromosomes)

    print(f"Chromosomes to plot ({n_chroms}): {chromosomes}")

    # Create grid layout
    n_cols = 3
    n_rows = int(np.ceil(n_chroms / n_cols))

    print(f"Grid layout: {n_rows} rows × {n_cols} columns = {n_rows * n_cols} panels")

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(18, 3 * n_rows), sharex=False, sharey=True)

    # Handle axes array properly
    if n_rows == 1 and n_cols == 1:
        axes = [axes]
    elif n_rows == 1 or n_cols == 1:
        axes = axes.flatten()
    else:
        axes = axes.flatten()

    # Plot each chromosome
    for idx, chrom in enumerate(chromosomes):
        ax = axes[idx]
        chrom_df = df[df['chrom'] == chrom].sort_values('pos')

        # Convert position to Mb for better readability
        pos_mb = chrom_df['pos'] / 1_000_000
        g_values = chrom_df[stat_col]

        # Plot line
        ax.plot(pos_mb, g_values, linewidth=0.5, alpha=0.7, color='#2E86AB')

        # Highlight significant peaks
        significant = g_values > threshold
        if significant.any():
            ax.scatter(pos_mb[significant], g_values[significant],
                      color='#A23B72', s=10, alpha=0.8, zorder=5,
                      label=f'>{threshold_percentile}th percentile')

        # Add threshold line
        ax.axhline(y=threshold, color='red', linestyle='--',
                  linewidth=1, alpha=0.5, label=f'{threshold_percentile}th percentile')

        # Formatting
        ax.set_title(f'{chrom}', fontsize=10, fontweight='bold')
        ax.set_xlabel('Position (Mb)', fontsize=9)
        ax.grid(True, alpha=0.3, linestyle=':')

        # Only show ylabel on leftmost plots
        if idx % n_cols == 0:
            ax.set_ylabel(stat_label, fontsize=9)

    # Hide unused subplots
    for idx in range(n_chroms, len(axes)):
        axes[idx].set_visible(False)

    # Add overall title
    fig.suptitle(f'BSA QTL Mapping - {stat_label} across Genome',
                 fontsize=14, fontweight='bold')

    # Add legend to first plot
    if n_chroms > 0:
        axes[0].legend(fontsize=8, loc='upper right')

    plt.tight_layout(rect=[0, 0, 1, 0.98])  # Leave room for suptitle

    if output_path:
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
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

    chromosomes = sorted(df['chrom'].unique())[:max_chromosomes]
    n_chroms = len(chromosomes)

    print(f"Chromosomes to plot ({n_chroms}): {chromosomes}")

    n_cols = 3
    n_rows = int(np.ceil(n_chroms / n_cols))

    print(f"Grid layout: {n_rows} rows × {n_cols} columns = {n_rows * n_cols} panels")

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(18, 3 * n_rows), sharex=False, sharey=True)

    # Handle axes array properly
    if n_rows == 1 and n_cols == 1:
        axes = [axes]
    elif n_rows == 1 or n_cols == 1:
        axes = axes.flatten()
    else:
        axes = axes.flatten()

    for idx, chrom in enumerate(chromosomes):
        ax = axes[idx]
        chrom_df = df[df['chrom'] == chrom].sort_values('pos')

        pos_mb = chrom_df['pos'] / 1_000_000
        delta = chrom_df['delta_snp_index']

        # Plot line
        ax.plot(pos_mb, delta, linewidth=0.5, alpha=0.7, color='#F18F01')

        # Add zero line
        ax.axhline(y=0, color='black', linestyle='-', linewidth=0.5, alpha=0.3)

        # Formatting
        ax.set_title(f'{chrom}', fontsize=10, fontweight='bold')
        ax.set_xlabel('Position (Mb)', fontsize=9)
        ax.set_ylim(-1.1, 1.1)
        ax.grid(True, alpha=0.3, linestyle=':')

        if idx % n_cols == 0:
            ax.set_ylabel('Δ(SNP-index)', fontsize=9)

    # Hide unused subplots
    for idx in range(n_chroms, len(axes)):
        axes[idx].set_visible(False)

    fig.suptitle('BSA - Delta SNP Index across Genome',
                 fontsize=14, fontweight='bold')

    plt.tight_layout(rect=[0, 0, 1, 0.98])  # Leave room for suptitle

    if output_path:
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        print(f"Plot saved to: {output_path}")
    else:
        plt.show()

    return fig


def identify_qtl_peaks(df, threshold_percentile=99, min_peak_width_bp=500_000):
    """
    Identify and report QTL peak regions.
    """
    stat_col = 'g_prime' if 'g_prime' in df.columns else 'g_statistic'
    threshold = df[stat_col].quantile(threshold_percentile / 100)

    print(f"\n{'='*70}")
    print(f"QTL Peak Detection (threshold: {threshold:.3f}, {threshold_percentile}th percentile)")
    print(f"{'='*70}")

    peaks = []
    for chrom in sorted(df['chrom'].unique()):
        chrom_df = df[df['chrom'] == chrom].sort_values('pos')
        significant = chrom_df[chrom_df[stat_col] > threshold]

        if len(significant) > 0:
            # Group nearby significant SNPs into peaks
            significant = significant.copy()
            significant['gap'] = significant['pos'].diff()
            significant['new_peak'] = significant['gap'] > min_peak_width_bp
            significant['peak_id'] = significant['new_peak'].cumsum()

            for peak_id in significant['peak_id'].unique():
                peak_snps = significant[significant['peak_id'] == peak_id]

                peak_info = {
                    'chrom': chrom,
                    'start': peak_snps['pos'].min(),
                    'end': peak_snps['pos'].max(),
                    'width_mb': (peak_snps['pos'].max() - peak_snps['pos'].min()) / 1_000_000,
                    'n_snps': len(peak_snps),
                    'max_g': peak_snps[stat_col].max(),
                    'mean_g': peak_snps[stat_col].mean(),
                }
                peaks.append(peak_info)

    # Sort by max G value
    peaks_df = pd.DataFrame(peaks).sort_values('max_g', ascending=False)

    print(f"\nFound {len(peaks_df)} QTL peaks:\n")
    print(peaks_df.to_string(index=False))
    print(f"\n{'='*70}\n")

    return peaks_df


def main():
    parser = argparse.ArgumentParser(description='Visualize BSA G-statistic results')
    parser.add_argument('input', help='Input CSV file from bsa-gprime')
    parser.add_argument('-o', '--output', help='Output plot file (PNG/PDF)')
    parser.add_argument('--smooth', action='store_true',
                       help='Apply smoothing to raw G-statistics')
    parser.add_argument('--window', type=int, default=2_000_000,
                       help='Smoothing window size in bp (default: 2,000,000)')
    parser.add_argument('--threshold', type=float, default=99,
                       help='Percentile threshold for significance (default: 99)')
    parser.add_argument('--delta-snp-index', action='store_true',
                       help='Also plot delta SNP index')
    parser.add_argument('--max-chroms', type=int, default=12,
                       help='Maximum number of chromosomes to plot (default: 12)')

    args = parser.parse_args()

    # Load data
    df = load_results(args.input)

    # Apply smoothing if requested
    if args.smooth:
        df = apply_smoothing(df, window_bp=args.window)
        use_smoothed = True
    else:
        use_smoothed = 'g_prime' in df.columns

    # Determine output path
    if args.output:
        output_base = Path(args.output).stem
        output_dir = Path(args.output).parent
        output_ext = Path(args.output).suffix
    else:
        output_base = 'qtl_plot'
        output_dir = Path('.')
        output_ext = '.png'

    # Create main G-statistic plot
    g_output = output_dir / f"{output_base}_gstat{output_ext}" if args.output else None
    plot_qtl_manhattan(df, output_path=g_output,
                      threshold_percentile=args.threshold,
                      use_smoothed=use_smoothed,
                      max_chromosomes=args.max_chroms)

    # Create delta SNP index plot if requested
    if args.delta_snp_index:
        delta_output = output_dir / f"{output_base}_delta{output_ext}" if args.output else None
        plot_delta_snp_index(df, output_path=delta_output,
                           max_chromosomes=args.max_chroms)

    # Identify and report peaks
    peaks = identify_qtl_peaks(df, threshold_percentile=args.threshold)

    # Save peaks to CSV
    if args.output:
        peaks_csv = output_dir / f"{output_base}_peaks.csv"
        peaks.to_csv(peaks_csv, index=False)
        print(f"Peak summary saved to: {peaks_csv}")

    if not args.output:
        plt.show()


if __name__ == '__main__':
    main()
