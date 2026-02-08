use anyhow::Result;
use plotters::prelude::*;

use super::downsample;
use super::faceted::ChromPanel;
use super::{COLOR_GRID, COLOR_MAGENTA, COLOR_RED, COLOR_STEEL_BLUE, N_COLS};
use crate::types::SignificanceResult;

/// Draw a -log10(p-value) Manhattan plot on the given drawing area.
pub fn draw_neg_log10_pvalue<DB: DrawingBackend>(
    root: &DrawingArea<DB, plotters::coord::Shift>,
    results: &[SignificanceResult],
    panels: &[ChromPanel],
) -> Result<()>
where
    DB::ErrorType: 'static,
{
    root.fill(&WHITE)?;

    // Title area + chart area
    let (title_area, chart_area) = root.split_vertically(60);
    title_area.titled(
        "BSA QTL Mapping - -log10(p-value) across Genome",
        ("sans-serif", 22).into_font().color(&BLACK),
    )?;

    let n_rows = (panels.len() + N_COLS - 1) / N_COLS;
    let panel_areas = chart_area.split_evenly((n_rows, N_COLS));

    // Compute -log10(p) for all results
    let neg_log_p: Vec<f64> = results
        .iter()
        .map(|r| -r.p_value.max(1e-300).log10())
        .collect();

    // Compute shared y-axis max
    let y_max = neg_log_p
        .iter()
        .copied()
        .fold(f64::NEG_INFINITY, f64::max)
        * 1.05;
    let y_max = if y_max.is_finite() && y_max > 0.0 {
        y_max
    } else {
        10.0
    };

    // Compute threshold line: -log10(p) at the boundary of significance
    // = the largest p-value among significant calls
    let threshold_line: Option<f64> = {
        let max_sig_p = results
            .iter()
            .filter(|r| r.significant)
            .map(|r| r.p_value)
            .fold(f64::NEG_INFINITY, f64::max);
        if max_sig_p > 0.0 && max_sig_p.is_finite() {
            Some(-max_sig_p.max(1e-300).log10())
        } else {
            None
        }
    };

    let pixel_width = downsample::panel_pixel_width(root.dim_in_pixel().0, N_COLS);

    for (idx, panel) in panels.iter().enumerate() {
        if idx >= panel_areas.len() {
            break;
        }

        let positions_mb: Vec<f64> = panel
            .indices
            .iter()
            .map(|&i| results[i].raw.variant.pos as f64 / 1_000_000.0)
            .collect();

        let x_min = positions_mb.first().copied().unwrap_or(0.0);
        let x_max = positions_mb.last().copied().unwrap_or(1.0);
        let x_margin = (x_max - x_min).max(0.1) * 0.02;

        let y_label = if idx % N_COLS == 0 {
            "-log10(p)"
        } else {
            ""
        };

        let mut chart = ChartBuilder::on(&panel_areas[idx])
            .caption(&panel.chrom, ("sans-serif", 14).into_font().color(&BLACK))
            .margin(5)
            .x_label_area_size(25)
            .y_label_area_size(if idx % N_COLS == 0 { 50 } else { 30 })
            .build_cartesian_2d(
                (x_min - x_margin)..(x_max + x_margin),
                0.0..y_max,
            )?;

        chart
            .configure_mesh()
            .x_desc("Position (Mb)")
            .y_desc(y_label)
            .x_label_style(("sans-serif", 10))
            .y_label_style(("sans-serif", 10))
            .light_line_style(COLOR_GRID.mix(0.3))
            .draw()?;

        // Background points (all variants) - small blue dots
        let bg_points: Vec<(f64, f64)> = panel
            .indices
            .iter()
            .enumerate()
            .map(|(j, &i)| (positions_mb[j], neg_log_p[i]))
            .collect();

        let bg_ds = downsample::minmax_downsample(&bg_points, pixel_width);
        let bg_draw = bg_ds.as_deref().unwrap_or(&bg_points);

        chart.draw_series(bg_draw.iter().map(|&(x, y)| {
            Circle::new((x, y), 1, COLOR_STEEL_BLUE.mix(0.4).filled())
        }))?;

        // Significant points - larger magenta dots
        let sig_points: Vec<(f64, f64)> = panel
            .indices
            .iter()
            .enumerate()
            .filter(|(_, &i)| results[i].significant)
            .map(|(j, &i)| (positions_mb[j], neg_log_p[i]))
            .collect();

        if !sig_points.is_empty() {
            chart.draw_series(sig_points.iter().map(|&(x, y)| {
                Circle::new((x, y), 2, COLOR_MAGENTA.mix(0.8).filled())
            }))?;
        }

        // Threshold dashed line
        if let Some(thresh) = threshold_line {
            chart.draw_series(DashedLineSeries::new(
                vec![
                    (x_min - x_margin, thresh),
                    (x_max + x_margin, thresh),
                ],
                5,
                3,
                COLOR_RED.mix(0.5).into(),
            ))?;
        }
    }

    Ok(())
}
