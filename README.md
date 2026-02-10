# BSA G-Prime: Rust Implementation for Bulk Segregant Analysis

A fast, efficient implementation of the G' (G-prime) statistic for Bulk Segregant Analysis using Next Generation Sequencing data.  Currently implemented for diploid organisms.

## Overview

This tool calculates the G-statistic for identifying QTL (Quantitative Trait Loci) from BSA-Seq data. It takes a VCF file with two bulk samples (high and low) and computes G-statistics, smoothed G' values, p-values, and FDR-corrected q-values to identify genomic regions associated with the trait of interest.

## Features

- Takes a VCF (or .vcf.gz) as input, with either 2 bulks (high + low) or 2 bulks + parents, and outputs a range of Bulk Segregant Analysis statistics as a CSV
- Supports G, G', ΔSNP index outputs
- Built-in plotting (with `--plot`) as well as optional python script for visualization
- Fast and multi-threaded (written in Rust)


## Installation

### Pre-built binaries

Go to the [releases](https://github.com/briannadon/bsa-gprime/releases) page and download the binary that matches your computing system (Linux, MacOS supported - Windows **not** supported).

If you want to use it from anywhere on the command line, place the binary in your $PATH. 

### Build from source

#### Prerequisites

- Rust 1.70+ (install from https://rustup.rs/)
- (Optional) Python3 with pandas, matplotlib, and numpy installed (for python plotting script)


```bash
# Clone the repository
cd bsa-gprime-rs

# Build in release mode (optimized)
cargo build --release

# The binary will be at: ./target/release/bsa-gprime
```

Plotting is enabled by default. To build without plotting support:

```bash
cargo build --release --no-default-features
```

## Usage

### Basic Usage

By default the tool runs the full significance testing pipeline (smoothing, p-values, FDR correction).

#### Bulks only

If your VCF contains only the two bulk samples, you can run the analysis directly. The first sample in the VCF is assumed to be the high bulk and the second the low bulk (or specify them explicitly with `--high-bulk` and `--low-bulk`):

```bash
bsa-gprime \
    --input data/bulks.vcf.gz \
    --output results/bsa_results.csv \
    --min-depth 10 \
    --min-gq 20 \
    --snps-only
```

#### Bulks + parental lines

If your VCF also contains the parental lines, provide their sample names with `--parent-high` and `--parent-low`. Parental genotypes are used to refine allele depth estimates in each bulk. You can optionally set separate depth/quality filters for the parents with `--min-parent-depth` and `--min-parent-gq` (they default to the bulk filter values):

```bash
bsa-gprime \
    --input data/cross.vcf.gz \
    --output results/bsa_results.csv \
    --high-bulk high_bulk \
    --low-bulk low_bulk \
    --parent-high parent_a \
    --parent-low parent_b \
    --min-depth 10 \
    --min-gq 20 \
    --snps-only
```

### Parametric null distribution

If you know the number of individuals per bulk, you can use the parametric method for the null distribution calculations (Magwene et al. 2011). Specify the number of individuals with `--bulk-size` and the organism ploidy with `--ploidy` (default: 2). The effective population parameter n_s is computed as `bulk_size * ploidy`.

```bash
bsa-gprime \
    --input data/bulks.vcf.gz \
    --output results/bsa_results.csv \
    --null-method parametric \
    --bulk-size 25 \
    --ploidy 2
```

### Raw G-statistics only (no significance testing)

To skip smoothing and significance testing entirely:

```bash
bsa-gprime \
    --input data/bulks.vcf.gz \
    --output results/bsa_results.csv \
    --no-significance
```

### Using --output-dir (auto-names the output file)

```bash
bsa-gprime \
    --input data/bulks.vcf.gz \
    --output-dir results/ \
    --snps-only
```

This creates `results/bulks_bsa_results.csv` automatically.

### Command-line Options

```
Options:
  -i, --input <INPUT>                  Input VCF file (can be gzipped)
  -o, --output <OUTPUT>                Output CSV file path
      --output-dir <OUTPUT_DIR>        Output directory (file auto-named from input)
      --min-depth <MIN_DEPTH>          Minimum read depth per sample [default: 10]
      --min-gq <MIN_GQ>               Minimum genotype quality [default: 20]
      --snps-only                      SNPs only (exclude INDELs)
  -q, --quiet                          Suppress progress output
      --no-significance                Disable significance testing pipeline
      --bandwidth <BANDWIDTH>          Tricube kernel bandwidth in bp [default: 1000000]
      --null-method <NULL_METHOD>      Null distribution method [default: nonparametric]
      --bulk-size <BULK_SIZE>          Number of individuals per bulk (required for parametric)
      --ploidy <PLOIDY>                Ploidy of the organism [default: 2]
      --recombination-rate <RATE>      Recombination rate in cM/Mb (parametric covariance term)
      --fdr-threshold <FDR_THRESHOLD>  FDR threshold for significance calls [default: 0.05]
      --threads <THREADS>              Number of threads for parallel processing [default: auto]
      --bed <BED>                      Output BED file of significant regions
      --peaks-csv <PEAKS_CSV>          Write QTL peak summary to CSV file
      --plot                           Generate plots after analysis
      --plot-from <CSV>                Generate plots from existing results CSV (skips VCF analysis)
      --plot-dir <PLOT_DIR>            Output directory for plots (defaults to output file's directory)
      --plot-format <FORMAT>           Plot output format: "png" or "svg" [default: png]
      --delta-snp-index                Also generate delta SNP-index plot
      --max-plot-chroms <N>            Maximum number of chromosomes to plot [default: 12]
  -V, --version                        Print version
  -h, --help                           Print help
```

Either `--output` or `--output-dir` must be provided when running analysis. When using `--output-dir`, the output file is automatically named `{input_stem}_bsa_results.csv`.

### Expected Input Format

Your VCF file should have:
- Exactly 2 samples (high/low, resistant/susceptible, etc.)
- The FIRST sample must be your "high" bulk, second is "low"
- FORMAT fields: GT, AD, DP, GQ
- Biallelic variants (single ALT allele)

Example VCF line:
```
Chr1  100543  .  A  G  .  .  .  GT:AD:DP:GQ  0/0:45,55:100:30  1/1:15,85:100:40
```

### Output Format

When significance testing is enabled (default), the output CSV contains 20 columns:

- `chrom`: Chromosome/contig name
- `pos`: Position (1-based)
- `ref`: Reference allele
- `alt`: Alternate allele
- `high_ref_depth`: REF allele depth in high bulk
- `high_alt_depth`: ALT allele depth in high bulk
- `high_dp`: Total depth in high bulk
- `high_gq`: Genotype quality in high bulk
- `low_ref_depth`: REF allele depth in low bulk
- `low_alt_depth`: ALT allele depth in low bulk
- `low_dp`: Total depth in low bulk
- `low_gq`: Genotype quality in low bulk
- `g_statistic`: Raw G-statistic value
- `snp_index_high`: SNP index for high bulk (ALT/(REF+ALT))
- `snp_index_low`: SNP index for low bulk
- `delta_snp_index`: Difference in SNP indices (low - high)
- `g_prime`: Smoothed G-statistic (tricube kernel)
- `p_value`: P-value from log-normal null distribution
- `q_value`: Benjamini-Hochberg adjusted p-value
- `significant`: Whether q_value < FDR threshold (`true`/`false`)

With `--no-significance`, the last 4 columns are omitted and only the first 16 are written.

### BED Output

When `--bed <PATH>` is provided (requires significance testing), a BED file of significant QTL regions is written for loading into genome browsers (IGV, UCSC, JBrowse, etc.). Adjacent significant SNPs on the same chromosome within `--bandwidth` bp are merged into contiguous regions.

```bash
./target/release/bsa-gprime \
    --input data/bulks.vcf.gz \
    --output results/bsa_results.csv \
    --bed results/qtl_regions.bed
```

Each line contains: `chrom  start  end  QTL_region_N  score  .` where coordinates are 0-based half-open (BED standard) and score is the maximum G' in the region (capped at 1000).

## Visualization

Built-in plotting generates faceted Manhattan plots directly from the Rust binary. No external dependencies are required.

You may optionally use the python `visualize_qtl.py` script in `scripts/` to plot the results from the output CSV. Requires Python3, Pandas, Numpy, Matplotlib.

### Plotting with analysis

Generate plots alongside the analysis:

```bash
./target/release/bsa-gprime \
    --input data/bulks.vcf.gz \
    --output results/bsa_results.csv \
    --plot \
    --delta-snp-index
```

This produces:
- `*_neglog10p.png` — -log10(p-value) Manhattan plot
- `*_gprime.png` — G' (smoothed G-statistic) Manhattan plot
- `*_delta_snp.png` — Delta SNP-index plot (with `--delta-snp-index`)

### Plotting from existing results

Re-generate or customize plots from a previously computed results CSV without re-running the full analysis:

```bash
./target/release/bsa-gprime \
    --plot-from results/bulks_bsa_results.csv \
    --plot-dir results/plots/ \
    --delta-snp-index \
    --plot-format png
```

### Plot options

| Flag | Description |
|---|---|
| `--plot` | Generate plots after analysis |
| `--plot-from <CSV>` | Plot from existing results CSV (no VCF needed) |
| `--plot-dir <DIR>` | Output directory for plots |
| `--plot-format <FMT>` | `png` (default) or `svg` |
| `--delta-snp-index` | Also generate delta SNP-index plot |
| `--max-plot-chroms <N>` | Max chromosomes to plot [default: 12] |

Data is automatically downsampled to the output resolution using min-max bucketing, so even datasets with millions of variants render in seconds.

### Legacy Python visualization

A Python script (`scripts/visualize_qtl.py`) is also included for alternative visualization. It requires Python 3.9+ with matplotlib and related dependencies (see `scripts/environment.yml` or `scripts/requirements.txt`). The built-in Rust plotting is recommended for most use cases.

## Statistical Details

All statistical calculations and theoretical underpinnings follow [Magwene, Willis, and Kelly 2011](https://pmc.ncbi.nlm.nih.gov/articles/PMC3207950/). Refer to equations 8-12.

1. **Tricube kernel smoothing**: G-statistics are smoothed per-chromosome using a tricube kernel with a configurable bandwidth (default 1 Mb). This produces the G' statistic from the raw G statistic values.

2. **Null distribution estimation**: The null distribution of G' is modeled as log-normal. Two methods are available:
   - **Non-parametric** (default): Robust estimation from the observed G' distribution using the median and left-MAD, with Hampel's outlier rule to exclude QTL signals. This method is recommended for most use cases as it directly estimates the null distribution from the data, and requires less input from you. It is also 
   - **Parametric**: Requires `--bulk-size` (and optionally `--ploidy`). Computes expected mean and variance using equations 8-9.
     - For smoothed G' values, variance is scaled by the sum of squared normalized kernel weights per equation 12. The implementation computes this automatically from the SNP density and smoothing window.

3. **P-value computation**: Survival function (1 - CDF) of the fitted log-normal distribution evaluated at each G' value.

4. **FDR correction**: Benjamini-Hochberg procedure to control the false discovery rate. SNPs with q-value below the threshold (default 0.05) are marked as significant.

### Statistics implementation

The implementation provides three methods for estimating the null distribution. Refer to 

1. **`estimate_null_parametric()`**: For raw G-statistics using equations 8-9
2. **`estimate_null_parametric_gprime()`**: For smoothed G' values, accounts for variance reduction via the sum of squared kernel weights
3. **`estimate_null_robust()`**: **Recommended** - Robust empirical estimator using Hampel's rule and mode estimation (step 3b of the paper's pipeline)

#### Robust Empirical Estimator (Recommended)

This method does the following:

1. **Log-transform**: W_G' = ln(G'_observed)
2. **Compute median & left-MAD**: MAD_l(W_G') = Median(|w_i - Median|) for w_i <= Median
3. **Apply Hampel's rule**: Remove outliers where w_i > Median + 5.2 x MAD_l
4. **Estimate mode**: Using kernel density estimation on trimmed data
5. **Compute log-normal parameters**:
   - mu = ln(Median of trimmed original values)
   - sigma^2 = mu - ln(Mode)

This approach:
- Infers the null distribution from observed G' data **without requiring n_s or C**
- Automatically identifies and removes QTL regions (as outliers)
- Is robust to model violations in the hierarchical sampling
- Matches Figure S1 validation from the original paper

#### Smoothing Factor

The sum of squared normalized kernel weights is computed as:
- k_j = w_j / sum(w_j) (normalized tricube weights)
- sum(k_j^2) = sum(w_j^2) / (sum(w_j))^2

For typical SNP densities (1 SNP per 1-10 kb), this ranges from 0.1 to 0.5, meaning smoothing reduces variance by a factor of 2-10.


## Performance

- Processes ~10M variants in 5-15 minutes on a modern laptop - even faster (2 min or so) with 4+ cores.
- Parallel G-statistic calculation, smoothing, and p-value computation via rayon (`--threads`)
- Handles gzipped VCF files natively (no decompression needed)
- Memory-efficient streaming parser
- Plotting uses min-max downsampling to handle millions of data points efficiently

## Interpreting Results

### G-statistic Values

- **G ~ 0**: No association between locus and trait (null hypothesis)
- **G > 0**: Evidence of association (higher values = stronger association)
- **Typical QTL regions**: G > 20-50
- **Strong QTL**: G > 100

### G' (Smoothed G-statistic)

G' is the tricube-smoothed version of the raw G-statistic. It reduces noise from individual SNPs and produces clearer QTL peaks. The bandwidth parameter controls the degree of smoothing.

### P-values and Q-values

- **p_value**: Probability of observing this G' or higher under the null distribution. Smaller = stronger evidence of a QTL.
- **q_value**: FDR-adjusted p-value (Benjamini-Hochberg). Controls the expected proportion of false positives among significant calls.
- **significant**: `true` if q_value < FDR threshold (default 0.05).

### SNP Index

- **SNP index = 0**: All reference alleles
- **SNP index = 1**: All alternate alleles
- **SNP index = 0.5**: Equal allele frequencies

### Delta SNP Index

- **Delta(SNP-index) ~ 0**: Similar allele frequencies in both bulks (no QTL)
- **Delta(SNP-index) > 0**: Susceptible bulk enriched for ALT allele
- **Delta(SNP-index) < 0**: Resistant bulk enriched for ALT allele

## Testing

Run the test suite:
```bash
cargo test
```

## License

This project is licensed under the GNU General Public License v3.0 - see the [LICENSE](LICENSE) file for details.

## References

```
Magwene PM, Willis JH, Kelly JK (2011) The Statistics of Bulk Segregant Analysis
Using Next Generation Sequencing. PLOS Computational Biology 7(11): e1002255.
```
## Support

For issues or questions, please open an issue on the GitHub repository.
