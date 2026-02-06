# BSA G-Prime: Rust Implementation for Bulk Segregant Analysis

A fast, efficient implementation of the G' (G-prime) statistic for Bulk Segregant Analysis using Next Generation Sequencing data.

## Overview

This tool calculates the G-statistic for identifying QTL (Quantitative Trait Loci) from BSA-Seq data. It takes a VCF file with two bulk samples (resistant and susceptible) and computes G-statistics, smoothed G' values, p-values, and FDR-corrected q-values to identify genomic regions associated with the trait of interest.

## Features

- Fast VCF parsing with native support for gzipped files
- G-statistic calculation for BSA-Seq data
- Tricube kernel smoothing of G-statistics (G' computation)
- Statistical significance testing via log-normal null distribution
- Parametric and non-parametric null distribution estimation
- Benjamini-Hochberg FDR correction
- SNP index and delta SNP index calculation
- Quality filtering (depth, genotype quality)
- Biallelic SNP filtering
- Parallel processing with configurable thread count (rayon)
- Progress bars with ETA for long-running steps
- CSV output for downstream analysis
- BED output of significant QTL regions for genome browsers

## Installation

### Prerequisites

- Rust 1.70+ (install from https://rustup.rs/)
- For visualization: Python 3.9+ with conda/mamba or pip

### Build from source

```bash
# Clone the repository
cd bsa-gprime-rs

# Build in release mode (optimized)
cargo build --release

# The binary will be at: ./target/release/bsa-gprime
```

## Usage

### Basic Usage

By default the tool runs the full significance testing pipeline (smoothing, p-values, FDR correction):

```bash
./target/release/bsa-gprime \
    --input data/bulks.vcf.gz \
    --output results/bsa_results.csv \
    --min-depth 10 \
    --min-gq 20 \
    --snps-only
```

### Parametric null distribution

If you know the number of individuals per bulk, you can use the parametric method (Magwene et al. 2011). Specify the number of individuals with `--bulk-size` and the organism ploidy with `--ploidy` (default: 2). The effective population parameter n_s is computed as `bulk_size * ploidy`.

```bash
./target/release/bsa-gprime \
    --input data/bulks.vcf.gz \
    --output results/bsa_results.csv \
    --null-method parametric \
    --bulk-size 25 \
    --ploidy 2
```

### Raw G-statistics only (no significance testing)

To skip smoothing and significance testing entirely:

```bash
./target/release/bsa-gprime \
    --input data/bulks.vcf.gz \
    --output results/bsa_results.csv \
    --no-significance
```

### Using --output-dir (auto-names the output file)

```bash
./target/release/bsa-gprime \
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
      --fdr-threshold <FDR_THRESHOLD>  FDR threshold for significance calls [default: 0.05]
      --threads <THREADS>              Number of threads for parallel processing [default: auto]
      --bed <BED>                      Output BED file of significant regions
  -V, --version                        Print version
  -h, --help                           Print help
```

Either `--output` or `--output-dir` must be provided. When using `--output-dir`, the output file is automatically named `{input_stem}_bsa_results.csv`.

### Expected Input Format

Your VCF file should have:
- Exactly 2 samples (resistant bulk and susceptible bulk)
- FORMAT fields: GT, AD, DP, GQ
- Biallelic variants (single ALT allele)

Example VCF line:
```
CP124720.1  100543  .  A  G  .  .  .  GT:AD:DP:GQ  0/0:45,55:100:30  1/1:15,85:100:40
```

### Output Format

When significance testing is enabled (default), the output CSV contains 20 columns:

- `chrom`: Chromosome/contig name
- `pos`: Position (1-based)
- `ref`: Reference allele
- `alt`: Alternate allele
- `resistant_ref_depth`: REF allele depth in resistant bulk
- `resistant_alt_depth`: ALT allele depth in resistant bulk
- `resistant_dp`: Total depth in resistant bulk
- `resistant_gq`: Genotype quality in resistant bulk
- `susceptible_ref_depth`: REF allele depth in susceptible bulk
- `susceptible_alt_depth`: ALT allele depth in susceptible bulk
- `susceptible_dp`: Total depth in susceptible bulk
- `susceptible_gq`: Genotype quality in susceptible bulk
- `g_statistic`: Raw G-statistic value
- `snp_index_resistant`: SNP index for resistant bulk (ALT/(REF+ALT))
- `snp_index_susceptible`: SNP index for susceptible bulk
- `delta_snp_index`: Difference in SNP indices (susceptible - resistant)
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

## Significance Testing

The significance testing pipeline follows Magwene et al. (2011):

1. **Tricube kernel smoothing**: G-statistics are smoothed per-chromosome using a tricube kernel with a configurable bandwidth (default 1 Mb). This produces the G' statistic.

2. **Null distribution estimation**: The null distribution of G' is modeled as log-normal. Two methods are available:
   - **Non-parametric** (default): Robust estimation from the observed G' distribution using the median and left-MAD, with Hampel's outlier rule to exclude QTL signals.
   - **Parametric**: Requires `--bulk-size` (and optionally `--ploidy`). Computes expected mean and variance of G under the null from the effective population size and average coverage.

3. **P-value computation**: Survival function (1 - CDF) of the fitted log-normal distribution evaluated at each G' value.

4. **FDR correction**: Benjamini-Hochberg procedure to control the false discovery rate. SNPs with q-value below the threshold (default 0.05) are marked as significant.

## Visualization

A Python script is provided in `scripts/` to visualize the results and identify QTL peaks.

### Setup Python Environment

**Using conda/mamba:**
```bash
mamba env create -f scripts/environment.yml
mamba activate bsa-gprime-viz
```

**Using pip:**
```bash
pip install -r scripts/requirements.txt
```

### Visualize Results

**Basic usage:**
```bash
python scripts/visualize_qtl.py bsa_results.csv
```

**With smoothing:**
```bash
python scripts/visualize_qtl.py bsa_results.csv \
    --smooth \
    --window 2000000 \
    --output qtl_peaks.png
```

**With delta SNP index plot:**
```bash
python scripts/visualize_qtl.py bsa_results.csv \
    --smooth \
    --delta-snp-index \
    --output qtl_analysis.png
```

### Visualization Options

```
Options:
  input                        Input CSV file from bsa-gprime
  -o, --output OUTPUT          Output plot file (PNG/PDF)
  --smooth                     Apply smoothing to raw G-statistics
  --window WINDOW              Smoothing window size in bp (default: 2,000,000)
  --threshold THRESHOLD        Percentile threshold (default: 99)
  --delta-snp-index            Also plot delta SNP index
  --max-chroms MAX_CHROMS      Maximum chromosomes to plot (default: 12)
```

The visualization script will:
1. Create faceted plots showing G-statistics across all chromosomes
2. Highlight peaks above the significance threshold
3. Identify and report QTL peak regions
4. Save a summary CSV of detected peaks

## Example Workflow

```bash
# 1. Run analysis with significance testing (default)
./target/release/bsa-gprime \
    --input data/bulks.vcf.gz \
    --output-dir results/ \
    --min-depth 10 \
    --min-gq 20 \
    --snps-only

# 2. Visualize results
python scripts/visualize_qtl.py results/bulks_bsa_results.csv \
    --smooth \
    --window 2000000 \
    --threshold 99 \
    --delta-snp-index \
    --output qtl_analysis.png

# 3. Check the outputs:
#    - qtl_analysis_gstat.png: G-statistic Manhattan plot
#    - qtl_analysis_delta.png: Delta SNP index plot
#    - qtl_analysis_peaks.csv: Summary of detected QTL peaks
```

## Performance

- Processes ~10M variants in 5-15 minutes on a modern laptop
- Parallel G-statistic calculation, smoothing, and p-value computation via rayon (`--threads`)
- Handles gzipped VCF files natively (no decompression needed)
- Memory-efficient streaming parser

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

Magwene PM, Willis JH, Kelly JK (2011) The Statistics of Bulk Segregant Analysis Using Next Generation Sequencing. *PLOS Computational Biology* 7(11): e1002255. https://doi.org/10.1371/journal.pcbi.1002255

## Citation

If you use this tool in your research, please cite the original G-statistic paper:
```
Magwene PM, Willis JH, Kelly JK (2011) The Statistics of Bulk Segregant Analysis
Using Next Generation Sequencing. PLOS Computational Biology 7(11): e1002255.
```

## Support

For issues or questions, please open an issue on the GitHub repository.
