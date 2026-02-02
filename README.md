# BSA G-Prime: Rust Implementation for Bulk Segregant Analysis

A fast, efficient implementation of the G' (G-prime) statistic for Bulk Segregant Analysis using Next Generation Sequencing data.

## Overview

This tool calculates the G-statistic for identifying QTL (Quantitative Trait Loci) from BSA-Seq data. It takes a VCF file with two bulk samples (resistant and susceptible) and computes G-statistics to identify genomic regions associated with the trait of interest.

**Reference**: Magwene et al. (2011) "The Statistics of Bulk Segregant Analysis Using Next Generation Sequencing" *PLOS Computational Biology* 7(11): e1002255

## Features

- ✅ Fast VCF parsing with native support for gzipped files
- ✅ G-statistic calculation for BSA-Seq data
- ✅ SNP index and delta SNP index calculation
- ✅ Quality filtering (depth, genotype quality)
- ✅ Biallelic SNP filtering
- ✅ Progress reporting for large datasets
- ✅ CSV output for downstream analysis
- 🚧 Tricube kernel smoothing (planned)
- 🚧 Significance testing (planned)

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

```bash
./target/release/bsa-gprime \
    --input results/gvcf/bulks.vcf.gz \
    --output bsa_results.csv \
    --min-depth 10 \
    --min-gq 20 \
    --snps-only
```

### Command-line Options

```
Options:
  -i, --input <INPUT>          Input VCF file (can be gzipped)
  -o, --output <OUTPUT>        Output CSV file
      --min-depth <MIN_DEPTH>  Minimum read depth per sample [default: 10]
      --min-gq <MIN_GQ>        Minimum genotype quality [default: 20]
      --snps-only              SNPs only (exclude INDELs)
  -h, --help                   Print help
```

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

The tool generates a CSV file with the following columns:

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
- `g_statistic`: G-statistic value
- `snp_index_resistant`: SNP index for resistant bulk (ALT/(REF+ALT))
- `snp_index_susceptible`: SNP index for susceptible bulk
- `delta_snp_index`: Difference in SNP indices (susceptible - resistant)

## Visualization

A Python script is provided to visualize the results and identify QTL peaks.

### Setup Python Environment

**Using conda/mamba:**
```bash
mamba env create -f environment.yml
conda activate bsa-gprime-viz
```

**Using pip:**
```bash
pip install -r requirements.txt
```

### Visualize Results

**Basic usage:**
```bash
python visualize_qtl.py bsa_results.csv
```

**With smoothing:**
```bash
python visualize_qtl.py bsa_results.csv \
    --smooth \
    --window 2000000 \
    --output qtl_peaks.png
```

**With delta SNP index plot:**
```bash
python visualize_qtl.py bsa_results.csv \
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
# 1. Run the G-statistic analysis
./target/release/bsa-gprime \
    --input my_bulks.vcf.gz \
    --output results.csv \
    --min-depth 10 \
    --min-gq 20 \
    --snps-only

# 2. Visualize results with smoothing
python visualize_qtl.py results.csv \
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
- Handles gzipped VCF files natively (no decompression needed)
- Memory-efficient streaming parser
- Expected to handle your 271MB gzipped VCF without issues

## Interpreting Results

### G-statistic Values

- **G ≈ 0**: No association between locus and trait (null hypothesis)
- **G > 0**: Evidence of association (higher values = stronger association)
- **G follows χ²(1) distribution** under null hypothesis
- **Typical QTL regions**: G > 20-50
- **Strong QTL**: G > 100

### SNP Index

- **SNP index = 0**: All reference alleles
- **SNP index = 1**: All alternate alleles
- **SNP index = 0.5**: Equal allele frequencies

### Delta SNP Index

- **Δ(SNP-index) ≈ 0**: Similar allele frequencies in both bulks (no QTL)
- **Δ(SNP-index) > 0**: Susceptible bulk enriched for ALT allele
- **Δ(SNP-index) < 0**: Resistant bulk enriched for ALT allele

## Current Limitations (MVP Version)

This is the MVP (Minimum Viable Product) version. Planned enhancements:

- [ ] Tricube kernel smoothing (currently done in Python visualization)
- [ ] Statistical significance testing
- [ ] Parallel processing by chromosome
- [ ] Progress bar with ETA
- [ ] Multiple testing correction
- [ ] BED format output for genome browsers

## Testing

Run the test suite:
```bash
cargo test
```

## License

This project implements algorithms described in:
Magwene PM, Willis JH, Kelly JK (2011) The Statistics of Bulk Segregant Analysis Using Next Generation Sequencing. PLOS Computational Biology 7(11): e1002255. https://doi.org/10.1371/journal.pcbi.1002255

## Citation

If you use this tool in your research, please cite the original G-statistic paper:
```
Magwene PM, Willis JH, Kelly JK (2011) The Statistics of Bulk Segregant Analysis
Using Next Generation Sequencing. PLOS Computational Biology 7(11): e1002255.
```

## Support

For issues or questions, please open an issue on the GitHub repository.
