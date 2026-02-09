use crate::types::{ParentalInfo, Variant};
use anyhow::{Context, Result};
use rust_htslib::bcf::record::GenotypeAllele;
use rust_htslib::{bcf, bcf::Read};
use std::path::Path;

macro_rules! progress {
    ($quiet:expr, $($arg:tt)*) => {
        if !$quiet {
            eprintln!($($arg)*);
        }
    };
}

/// Configuration for parental sample filtering
pub struct ParentalConfig {
    pub parent_high_name: String, // parent of the high/resistant bulk
    pub parent_low_name: String,  // parent of the low/susceptible bulk
    pub min_depth: u32,
    pub min_gq: u32,
}

/// Configuration for sample name resolution
pub struct SampleConfig {
    pub high_bulk_name: Option<String>,
    pub low_bulk_name: Option<String>,
    pub parental: Option<ParentalConfig>,
}

/// Statistics from parental filtering
#[derive(Default)]
struct ParentalFilterStats {
    missing_gt: u64,
    heterozygous: u64,
    monomorphic: u64,
    low_depth: u64,
    low_gq: u64,
}

impl ParentalFilterStats {
    fn total(&self) -> u64 {
        self.missing_gt + self.heterozygous + self.monomorphic + self.low_depth + self.low_gq
    }
}

/// Resolve a sample name to its column index in the VCF header.
fn resolve_sample_index(header: &bcf::header::HeaderView, name: &str) -> Result<usize> {
    let samples: Vec<String> = header
        .samples()
        .iter()
        .map(|s| String::from_utf8_lossy(s).to_string())
        .collect();
    samples
        .iter()
        .position(|s| s == name)
        .with_context(|| {
            format!(
                "Sample '{}' not found in VCF header. Available samples: {:?}",
                name, samples
            )
        })
}

/// Check whether a genotype is homozygous-ref (0/0).
fn is_hom_ref(gt: &rust_htslib::bcf::record::Genotype) -> bool {
    gt.iter().all(|a| matches!(a, GenotypeAllele::Unphased(0) | GenotypeAllele::Phased(0)))
}

/// Check whether a genotype is homozygous-alt (1/1).
fn is_hom_alt(gt: &rust_htslib::bcf::record::Genotype) -> bool {
    gt.iter().all(|a| matches!(a, GenotypeAllele::Unphased(1) | GenotypeAllele::Phased(1)))
}

/// Check whether a genotype has any missing allele.
fn has_missing_allele(gt: &rust_htslib::bcf::record::Genotype) -> bool {
    gt.iter().any(|a| matches!(a, GenotypeAllele::UnphasedMissing | GenotypeAllele::PhasedMissing))
}

/// Format a genotype as a string like "0/0" or "1/1".
fn gt_to_string(gt: &rust_htslib::bcf::record::Genotype) -> String {
    let alleles: Vec<String> = gt
        .iter()
        .map(|a| match a {
            GenotypeAllele::Unphased(n) | GenotypeAllele::Phased(n) => n.to_string(),
            _ => ".".to_string(),
        })
        .collect();
    alleles.join("/")
}

pub fn parse_vcf(
    path: &Path,
    min_depth: u32,
    min_gq: u32,
    snps_only: bool,
    quiet: bool,
    sample_config: Option<&SampleConfig>,
) -> Result<Vec<Variant>> {
    let mut reader = bcf::Reader::from_path(path).context("Failed to open VCF file")?;

    let header = reader.header().clone();
    let sample_count = header.sample_count() as usize;

    let sample_names: Vec<String> = header
        .samples()
        .iter()
        .map(|s| String::from_utf8_lossy(s).to_string())
        .collect();

    // Resolve bulk sample indices
    let (high_idx, low_idx) = match sample_config {
        Some(cfg) => {
            let hi = match &cfg.high_bulk_name {
                Some(name) => resolve_sample_index(&header, name)?,
                None => {
                    if sample_count < 2 {
                        anyhow::bail!("VCF has {} sample(s); need at least 2", sample_count);
                    }
                    0
                }
            };
            let lo = match &cfg.low_bulk_name {
                Some(name) => resolve_sample_index(&header, name)?,
                None => {
                    if sample_count < 2 {
                        anyhow::bail!("VCF has {} sample(s); need at least 2", sample_count);
                    }
                    1
                }
            };
            (hi, lo)
        }
        None => {
            // Legacy: exactly 2 samples, positional
            if sample_count != 2 {
                anyhow::bail!(
                    "VCF has {} samples. Use --high-bulk and --low-bulk to specify sample names \
                     when the VCF contains more than 2 samples.",
                    sample_count
                );
            }
            (0, 1)
        }
    };

    progress!(quiet, "High bulk (sample {}): {}", high_idx, sample_names[high_idx]);
    progress!(quiet, "Low bulk  (sample {}): {}", low_idx, sample_names[low_idx]);

    // Resolve parental sample indices (if provided)
    let parental_indices = if let Some(cfg) = sample_config {
        if let Some(ref pcfg) = cfg.parental {
            let phi = resolve_sample_index(&header, &pcfg.parent_high_name)?;
            let plo = resolve_sample_index(&header, &pcfg.parent_low_name)?;
            progress!(quiet, "Parent-high (sample {}): {}", phi, sample_names[phi]);
            progress!(quiet, "Parent-low  (sample {}): {}", plo, sample_names[plo]);
            Some((phi, plo, pcfg.min_depth, pcfg.min_gq))
        } else {
            None
        }
    } else {
        None
    };

    let mut variants = Vec::new();
    let mut total_count: u64 = 0;
    let mut filtered_count: u64 = 0;
    let mut parental_stats = ParentalFilterStats::default();

    for result in reader.records() {
        total_count += 1;
        if total_count % 100_000 == 0 {
            progress!(quiet, "Processed {} variants...", total_count);
        }

        let record = result.context("Failed to read VCF record")?;

        // Get basic variant info
        let rid = record.rid().context("No reference ID")?;
        let chrom = String::from_utf8_lossy(header.rid2name(rid)?).to_string();
        let pos = record.pos() as u64 + 1; // BCF is 0-based, we want 1-based

        // Get alleles
        let alleles = record.alleles();
        if alleles.len() != 2 {
            // Not biallelic
            filtered_count += 1;
            continue;
        }

        let ref_allele = String::from_utf8_lossy(alleles[0]).to_string();
        let alt_allele = String::from_utf8_lossy(alleles[1]).to_string();

        // Skip INDELs if requested
        if snps_only && (ref_allele.len() != 1 || alt_allele.len() != 1) {
            filtered_count += 1;
            continue;
        }

        // Helper to extract sample data: (ref_depth, alt_depth, dp, gq)
        let extract_sample_data = |sample_idx: usize| -> Option<(u32, u32, u32, u32)> {
            let genotypes = record.genotypes().context("No genotypes").ok()?;
            let gt = genotypes.get(sample_idx);

            if has_missing_allele(&gt) {
                return None;
            }

            let ad = record.format(b"AD").integer().ok()?;
            let ad_sample = ad.get(sample_idx)?;
            if ad_sample.len() < 2 {
                return None;
            }
            let ref_depth = ad_sample[0] as u32;
            let alt_depth = ad_sample[1] as u32;

            let dp = record.format(b"DP").integer().ok()?;
            let dp_sample = dp.get(sample_idx)?;
            if dp_sample.is_empty() {
                return None;
            }
            let dp_val = dp_sample[0] as u32;

            let gq = record.format(b"GQ").integer().ok()?;
            let gq_sample = gq.get(sample_idx)?;
            if gq_sample.is_empty() {
                return None;
            }
            let gq_val = gq_sample[0] as u32;

            Some((ref_depth, alt_depth, dp_val, gq_val))
        };

        // ── Parental filtering ──
        let mut parental_info: Option<ParentalInfo> = None;
        let mut should_swap = false;

        if let Some((phi, plo, p_min_depth, p_min_gq)) = parental_indices {
            let genotypes = match record.genotypes() {
                Ok(g) => g,
                Err(_) => {
                    parental_stats.missing_gt += 1;
                    continue;
                }
            };

            let gt_high = genotypes.get(phi);
            let gt_low = genotypes.get(plo);

            // 1. Skip if either parent has missing genotype
            if has_missing_allele(&gt_high) || has_missing_allele(&gt_low) {
                parental_stats.missing_gt += 1;
                continue;
            }

            // 2. Skip if either parent is heterozygous (both must be homozygous)
            let high_hom_ref = is_hom_ref(&gt_high);
            let high_hom_alt = is_hom_alt(&gt_high);
            let low_hom_ref = is_hom_ref(&gt_low);
            let low_hom_alt = is_hom_alt(&gt_low);

            if !(high_hom_ref || high_hom_alt) || !(low_hom_ref || low_hom_alt) {
                parental_stats.heterozygous += 1;
                continue;
            }

            // 3. Skip if both parents have same genotype (uninformative)
            if (high_hom_ref && low_hom_ref) || (high_hom_alt && low_hom_alt) {
                parental_stats.monomorphic += 1;
                continue;
            }

            // 4. Apply parental min_depth / min_gq thresholds
            let check_parental_quality = |idx: usize| -> Option<bool> {
                let dp = record.format(b"DP").integer().ok()?;
                let dp_sample = dp.get(idx)?;
                if dp_sample.is_empty() {
                    return Some(false);
                }
                let dp_val = dp_sample[0] as u32;

                let gq = record.format(b"GQ").integer().ok()?;
                let gq_sample = gq.get(idx)?;
                if gq_sample.is_empty() {
                    return Some(false);
                }
                let gq_val = gq_sample[0] as u32;

                Some(dp_val >= p_min_depth && gq_val >= p_min_gq)
            };

            let ph_ok = check_parental_quality(phi).unwrap_or(false);
            let pl_ok = check_parental_quality(plo).unwrap_or(false);

            if !ph_ok {
                if !pl_ok {
                    // Both fail - count once
                    parental_stats.low_depth += 1;
                } else {
                    parental_stats.low_depth += 1;
                }
                continue;
            }
            if !pl_ok {
                parental_stats.low_depth += 1;
                continue;
            }

            // AD swapping: if parent-high is ALT-homozygous (1/1), swap REF/ALT
            should_swap = high_hom_alt;

            parental_info = Some(ParentalInfo {
                parent_high_gt: gt_to_string(&gt_high),
                parent_low_gt: gt_to_string(&gt_low),
                alleles_swapped: should_swap,
            });
        }

        // Extract bulk sample data
        let high_data = match extract_sample_data(high_idx) {
            Some(data) => data,
            None => {
                filtered_count += 1;
                continue;
            }
        };

        let low_data = match extract_sample_data(low_idx) {
            Some(data) => data,
            None => {
                filtered_count += 1;
                continue;
            }
        };

        let (mut high_ref, mut high_alt, high_dp, high_gq) = high_data;
        let (mut low_ref, mut low_alt, low_dp, low_gq) = low_data;

        // Apply bulk filters
        if high_dp < min_depth || low_dp < min_depth {
            filtered_count += 1;
            continue;
        }
        if high_gq < min_gq || low_gq < min_gq {
            filtered_count += 1;
            continue;
        }

        // Apply AD swap if needed (parent-high is 1/1 -> swap so REF = parent-high allele)
        let (final_ref_allele, final_alt_allele) = if should_swap {
            std::mem::swap(&mut high_ref, &mut high_alt);
            std::mem::swap(&mut low_ref, &mut low_alt);
            (alt_allele, ref_allele)
        } else {
            (ref_allele, alt_allele)
        };

        variants.push(Variant {
            chrom,
            pos,
            ref_allele: final_ref_allele,
            alt_allele: final_alt_allele,
            high_ref_depth: high_ref,
            high_alt_depth: high_alt,
            high_dp,
            high_gq,
            low_ref_depth: low_ref,
            low_alt_depth: low_alt,
            low_dp,
            low_gq,
            parental_info,
        });
    }

    progress!(quiet, "Total variants: {}", total_count);
    progress!(quiet, "Filtered variants (bulk QC): {}", filtered_count);

    if parental_indices.is_some() {
        progress!(quiet, "Parental filters:");
        progress!(quiet, "  Missing genotype: {}", parental_stats.missing_gt);
        progress!(quiet, "  Heterozygous parent: {}", parental_stats.heterozygous);
        progress!(quiet, "  Monomorphic (same GT): {}", parental_stats.monomorphic);
        progress!(quiet, "  Low depth/GQ: {}", parental_stats.low_depth);
        progress!(quiet, "  Total parental filtered: {}", parental_stats.total());

        let swapped = variants.iter().filter(|v| {
            v.parental_info.as_ref().map_or(false, |p| p.alleles_swapped)
        }).count();
        progress!(quiet, "  AD-swapped sites: {} ({:.1}%)",
            swapped,
            100.0 * swapped as f64 / variants.len().max(1) as f64
        );
    }

    progress!(
        quiet,
        "Retained variants: {} ({:.1}%)",
        variants.len(),
        100.0 * variants.len() as f64 / total_count.max(1) as f64
    );

    Ok(variants)
}
