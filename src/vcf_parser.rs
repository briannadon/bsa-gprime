use crate::types::Variant;
use anyhow::{Context, Result};
use rust_htslib::bcf::record::GenotypeAllele;
use rust_htslib::{bcf, bcf::Read};
use std::path::Path;

pub fn parse_vcf(
    path: &Path,
    min_depth: u32,
    min_gq: u32,
    snps_only: bool,
) -> Result<Vec<Variant>> {
    let mut reader = bcf::Reader::from_path(path).context("Failed to open VCF file")?;

    let header = reader.header().clone();

    // Get sample names - should be 2 samples
    let sample_count = header.sample_count();
    if sample_count != 2 {
        anyhow::bail!("Expected exactly 2 samples, found {}", sample_count);
    }

    let sample_names: Vec<_> = header.samples().iter().map(|s| String::from_utf8_lossy(s).to_string()).collect();
    eprintln!("Sample 0: {}", sample_names[0]);
    eprintln!("Sample 1: {}", sample_names[1]);

    let mut variants = Vec::new();
    let mut total_count = 0;
    let mut filtered_count = 0;

    for result in reader.records() {
        total_count += 1;
        if total_count % 100_000 == 0 {
            eprintln!("Processed {} variants...", total_count);
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
        if snps_only {
            if ref_allele.len() != 1 || alt_allele.len() != 1 {
                filtered_count += 1;
                continue;
            }
        }

        // Helper function to extract sample data
        let extract_sample_data = |sample_idx: usize| -> Option<(u32, u32, u32, u32)> {
            // Get GT
            let genotypes = record.genotypes().context("No genotypes").ok()?;
            let gt = genotypes.get(sample_idx);

            // Check if genotype is called (not missing)
            let mut has_missing = false;
            for allele in gt.iter() {
                if let GenotypeAllele::UnphasedMissing | GenotypeAllele::PhasedMissing = allele {
                    has_missing = true;
                    break;
                }
            }
            if has_missing {
                return None;
            }

            // Get AD (allelic depths)
            let ad = record.format(b"AD").integer().ok()?;
            let ad_sample = ad.get(sample_idx)?;
            if ad_sample.len() < 2 {
                return None;
            }
            let ref_depth = ad_sample[0] as u32;
            let alt_depth = ad_sample[1] as u32;

            // Get DP
            let dp = record.format(b"DP").integer().ok()?;
            let dp_sample = dp.get(sample_idx)?;
            if dp_sample.is_empty() {
                return None;
            }
            let dp_val = dp_sample[0] as u32;

            // Get GQ
            let gq = record.format(b"GQ").integer().ok()?;
            let gq_sample = gq.get(sample_idx)?;
            if gq_sample.is_empty() {
                return None;
            }
            let gq_val = gq_sample[0] as u32;

            Some((ref_depth, alt_depth, dp_val, gq_val))
        };

        // Extract data for both samples
        let resistant_data = match extract_sample_data(0) {
            Some(data) => data,
            None => {
                filtered_count += 1;
                continue;
            }
        };

        let susceptible_data = match extract_sample_data(1) {
            Some(data) => data,
            None => {
                filtered_count += 1;
                continue;
            }
        };

        let (resistant_ref, resistant_alt, resistant_dp, resistant_gq) = resistant_data;
        let (susceptible_ref, susceptible_alt, susceptible_dp, susceptible_gq) =
            susceptible_data;

        // Apply filters
        if resistant_dp < min_depth || susceptible_dp < min_depth {
            filtered_count += 1;
            continue;
        }
        if resistant_gq < min_gq || susceptible_gq < min_gq {
            filtered_count += 1;
            continue;
        }

        // Create variant
        variants.push(Variant {
            chrom,
            pos: pos as u64,
            ref_allele,
            alt_allele,
            resistant_ref_depth: resistant_ref,
            resistant_alt_depth: resistant_alt,
            resistant_dp,
            resistant_gq,
            susceptible_ref_depth: susceptible_ref,
            susceptible_alt_depth: susceptible_alt,
            susceptible_dp,
            susceptible_gq,
        });
    }

    eprintln!("Total variants: {}", total_count);
    eprintln!("Filtered variants: {}", filtered_count);
    eprintln!(
        "Retained variants: {} ({:.1}%)",
        variants.len(),
        100.0 * variants.len() as f64 / total_count as f64
    );

    Ok(variants)
}
