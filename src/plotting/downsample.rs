/// Estimate the drawable pixel width per panel given total image width and column count.
///
/// Accounts for axis labels, margins, and inter-panel spacing (~60px overhead per panel).
pub fn panel_pixel_width(total_width: u32, n_cols: usize) -> usize {
    let per_panel = total_width as usize / n_cols.max(1);
    per_panel.saturating_sub(60)
}

/// Min-max downsample a slice of (x, y) points into `n_buckets` equal-width x-bins.
///
/// For each bucket, retains only the points with the minimum and maximum y-value,
/// preserving all visual peaks and valleys at the target resolution. Returns `None`
/// if the input already fits (fewer than `2 * n_buckets` points), so the caller can
/// use the original data unchanged. First and last points are always preserved.
///
/// The output is sorted by x (within each bucket, min-y point comes before max-y point
/// unless they are the same point).
pub fn minmax_downsample(points: &[(f64, f64)], n_buckets: usize) -> Option<Vec<(f64, f64)>> {
    let n = points.len();
    if n_buckets == 0 || n <= 2 * n_buckets {
        return None;
    }

    // Determine x range
    let (mut x_min, mut x_max) = (f64::INFINITY, f64::NEG_INFINITY);
    for &(x, _) in points {
        if x < x_min {
            x_min = x;
        }
        if x > x_max {
            x_max = x;
        }
    }

    let x_range = x_max - x_min;
    if x_range <= 0.0 {
        // All points share the same x; just return first and last
        return Some(vec![points[0], points[n - 1]]);
    }

    let bucket_width = x_range / n_buckets as f64;

    // Each bucket tracks the point with min-y and max-y
    struct Bucket {
        min_pt: (f64, f64),
        max_pt: (f64, f64),
        has_data: bool,
    }

    let mut buckets: Vec<Bucket> = (0..n_buckets)
        .map(|_| Bucket {
            min_pt: (0.0, f64::INFINITY),
            max_pt: (0.0, f64::NEG_INFINITY),
            has_data: false,
        })
        .collect();

    // Single pass over all points
    for &(x, y) in points {
        let bi = ((x - x_min) / bucket_width) as usize;
        let bi = bi.min(n_buckets - 1); // clamp last edge
        let b = &mut buckets[bi];
        b.has_data = true;
        if y < b.min_pt.1 {
            b.min_pt = (x, y);
        }
        if y > b.max_pt.1 {
            b.max_pt = (x, y);
        }
    }

    // Collect results; for each bucket emit min then max (skip duplicate if same point)
    let mut out = Vec::with_capacity(2 * n_buckets + 2);

    // Always include the very first point
    out.push(points[0]);

    for b in &buckets {
        if !b.has_data {
            continue;
        }
        // Emit in x-order within the bucket
        if b.min_pt.0 <= b.max_pt.0 {
            out.push(b.min_pt);
            if b.min_pt.0 != b.max_pt.0 || b.min_pt.1 != b.max_pt.1 {
                out.push(b.max_pt);
            }
        } else {
            out.push(b.max_pt);
            if b.min_pt.0 != b.max_pt.0 || b.min_pt.1 != b.max_pt.1 {
                out.push(b.min_pt);
            }
        }
    }

    // Always include the very last point
    let last = points[n - 1];
    if let Some(&prev) = out.last() {
        if prev.0 != last.0 || prev.1 != last.1 {
            out.push(last);
        }
    }

    Some(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn small_input_returns_none() {
        let pts: Vec<(f64, f64)> = (0..10).map(|i| (i as f64, (i as f64).sin())).collect();
        assert!(minmax_downsample(&pts, 10).is_none());
    }

    #[test]
    fn large_input_downsamples() {
        let pts: Vec<(f64, f64)> = (0..10_000)
            .map(|i| (i as f64, (i as f64 * 0.01).sin()))
            .collect();
        let ds = minmax_downsample(&pts, 100).unwrap();
        // Should have at most ~202 points (2 per bucket + first/last)
        assert!(ds.len() <= 210);
        assert!(ds.len() >= 100); // at least one per bucket
    }

    #[test]
    fn preserves_extremes() {
        // Create data with a clear spike
        let mut pts: Vec<(f64, f64)> = (0..1000).map(|i| (i as f64, 1.0)).collect();
        pts[500] = (500.0, 100.0); // spike
        let ds = minmax_downsample(&pts, 50).unwrap();
        // The spike should be preserved
        assert!(ds.iter().any(|&(_, y)| y == 100.0));
    }

    #[test]
    fn panel_pixel_width_basic() {
        assert_eq!(panel_pixel_width(1800, 3), 540);
    }
}
