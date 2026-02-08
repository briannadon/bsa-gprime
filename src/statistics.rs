/// Calculate G-statistic from 2x2 contingency table
///
/// Contingency table:
///                 Resistant  |  Susceptible
/// REF allele:        o1      |      o3
/// ALT allele:        o2      |      o4
///
/// Returns: G-statistic value (>= 0)
pub fn calculate_g_statistic(o1: u32, o2: u32, o3: u32, o4: u32) -> f64 {
    // Convert to f64 for calculations
    let o1 = o1 as f64;
    let o2 = o2 as f64;
    let o3 = o3 as f64;
    let o4 = o4 as f64;

    let total = o1 + o2 + o3 + o4;

    // Handle edge case: no data
    if total == 0.0 {
        return 0.0;
    }

    // Calculate expected values under independence
    let e1 = (o1 + o2) * (o1 + o3) / total;
    let e2 = (o1 + o2) * (o2 + o4) / total;
    let e3 = (o3 + o4) * (o1 + o3) / total;
    let e4 = (o3 + o4) * (o2 + o4) / total;

    // Log-likelihood ratio components
    // Note: lim(x→0) x*ln(x) = 0
    let llr = |observed: f64, expected: f64| -> f64 {
        if observed == 0.0 || expected == 0.0 {
            0.0
        } else {
            2.0 * observed * (observed / expected).ln()
        }
    };

    let g = llr(o1, e1) + llr(o2, e2) + llr(o3, e3) + llr(o4, e4);

    // G-statistic should be non-negative
    g.max(0.0)
}

/// Calculate SNP index for a bulk
/// SNP_index = ALT / (REF + ALT)
pub fn snp_index(ref_depth: u32, alt_depth: u32) -> f64 {
    let total = ref_depth + alt_depth;
    if total == 0 {
        f64::NAN
    } else {
        alt_depth as f64 / total as f64
    }
}

/// Calculate delta SNP index
/// Δ(SNP-index) = SNP_index_low - SNP_index_high
pub fn delta_snp_index(
    high_ref: u32,
    high_alt: u32,
    low_ref: u32,
    low_alt: u32,
) -> f64 {
    let si_high = snp_index(high_ref, high_alt);
    let si_low = snp_index(low_ref, low_alt);

    if si_high.is_nan() || si_low.is_nan() {
        f64::NAN
    } else {
        si_low - si_high
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_perfect_association() {
        // Resistant: all REF, Susceptible: all ALT
        let g = calculate_g_statistic(100, 0, 0, 100);
        assert!(g > 100.0); // Should be very high
    }

    #[test]
    fn test_no_association() {
        // Equal allele frequencies in both bulks
        let g = calculate_g_statistic(50, 50, 50, 50);
        assert!(g < 0.1); // Should be near zero
    }

    #[test]
    fn test_all_zeros() {
        let g = calculate_g_statistic(0, 0, 0, 0);
        assert_eq!(g, 0.0);
    }

    #[test]
    fn test_snp_index() {
        assert_relative_eq!(snp_index(50, 50), 0.5, epsilon = 0.001);
        assert_relative_eq!(snp_index(100, 0), 0.0, epsilon = 0.001);
        assert_relative_eq!(snp_index(0, 100), 1.0, epsilon = 0.001);
        assert!(snp_index(0, 0).is_nan());
    }
}
