use anyhow::Result;
use plotters::prelude::*;

use super::faceted::ChromPanel;
use super::{COLOR_GRID, COLOR_MAGENTA, COLOR_RED, COLOR_STEEL_BLUE, N_COLS};
use crate::types::{GStatisticResult, SignificanceResult};

/// Draw a G-prime plot from significance results.
pub fn draw_gprime_sig<DB: DrawingBackend>(
    root: &DrawingArea<DB, plotters::coord::Shift>,
    results: &[SignificanceResult],
    panels: &[ChromPanel],
    threshold_percentile: f64,
) -> Result<()>
where
    DB::ErrorType: 'static,
{
    root.fill(&WHITE)?;

    let (title_area, chart_area) = root.split_vertically(60);
    title_area.titled(
        "BSA QTL Mapping - G' (smoothed) across Genome",
        ("sans-serif", 22).into_font().color(&BLACK),
    )?;

    let n_rows = (panels.len() + N_COLS - 1) / N_COLS;
    let panel_areas = chart_area.split_evenly((n_rows, N_COLS));

    // Compute percentile threshold on g_prime
    let mut g_values: Vec<f64> = results.iter().map(|r| r.g_prime).collect();
    g_values.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let threshold = if g_values.is_empty() {
        0.0
    } else {
        let idx = ((threshold_percentile / 100.0) * g_values.len() as f64) as usize;
        g_values[idx.min(g_values.len() - 1)]
    };

    // Shared y-axis max
    let y_max = g_values.last().copied().unwrap_or(10.0) * 1.05;
    let y_max = if y_max.is_finite() && y_max > 0.0 {
        y_max
    } else {
        10.0
    };

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

        let y_label = if idx % N_COLS == 0 { "G'" } else { "" };

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

        // Line plot of G' values
        let line_data: Vec<(f64, f64)> = panel
            .indices
            .iter()
            .enumerate()
            .map(|(j, &i)| (positions_mb[j], results[i].g_prime))
            .collect();

        chart.draw_series(LineSeries::new(
            line_data.iter().copied(),
            COLOR_STEEL_BLUE.mix(0.7).stroke_width(1),
        ))?;

        // Scatter points above threshold
        let sig_points: Vec<(f64, f64)> = line_data
            .iter()
            .filter(|(_, y)| *y > threshold)
            .copied()
            .collect();

        if !sig_points.is_empty() {
            chart.draw_series(sig_points.iter().map(|&(x, y)| {
                Circle::new((x, y), 2, COLOR_MAGENTA.mix(0.8).filled())
            }))?;
        }

        // Threshold dashed line
        chart.draw_series(DashedLineSeries::new(
            vec![
                (x_min - x_margin, threshold),
                (x_max + x_margin, threshold),
            ],
            5,
            3,
            COLOR_RED.mix(0.5).into(),
        ))?;
    }

    Ok(())
}

/// Draw a G-statistic plot from basic results (no significance testing).
pub fn draw_gstat<DB: DrawingBackend>(
    root: &DrawingArea<DB, plotters::coord::Shift>,
    results: &[GStatisticResult],
    panels: &[ChromPanel],
    threshold_percentile: f64,
) -> Result<()>
where
    DB::ErrorType: 'static,
{
    root.fill(&WHITE)?;

    let (title_area, chart_area) = root.split_vertically(60);
    title_area.titled(
        "BSA QTL Mapping - G-statistic across Genome",
        ("sans-serif", 22).into_font().color(&BLACK),
    )?;

    let n_rows = (panels.len() + N_COLS - 1) / N_COLS;
    let panel_areas = chart_area.split_evenly((n_rows, N_COLS));

    // Compute percentile threshold
    let mut g_values: Vec<f64> = results.iter().map(|r| r.g_statistic).collect();
    g_values.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let threshold = if g_values.is_empty() {
        0.0
    } else {
        let idx = ((threshold_percentile / 100.0) * g_values.len() as f64) as usize;
        g_values[idx.min(g_values.len() - 1)]
    };

    let y_max = g_values.last().copied().unwrap_or(10.0) * 1.05;
    let y_max = if y_max.is_finite() && y_max > 0.0 {
        y_max
    } else {
        10.0
    };

    for (idx, panel) in panels.iter().enumerate() {
        if idx >= panel_areas.len() {
            break;
        }

        let positions_mb: Vec<f64> = panel
            .indices
            .iter()
            .map(|&i| results[i].variant.pos as f64 / 1_000_000.0)
            .collect();

        let x_min = positions_mb.first().copied().unwrap_or(0.0);
        let x_max = positions_mb.last().copied().unwrap_or(1.0);
        let x_margin = (x_max - x_min).max(0.1) * 0.02;

        let y_label = if idx % N_COLS == 0 {
            "G-statistic"
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

        // Line plot
        let line_data: Vec<(f64, f64)> = panel
            .indices
            .iter()
            .enumerate()
            .map(|(j, &i)| (positions_mb[j], results[i].g_statistic))
            .collect();

        chart.draw_series(LineSeries::new(
            line_data.iter().copied(),
            COLOR_STEEL_BLUE.mix(0.7).stroke_width(1),
        ))?;

        // Points above threshold
        let sig_points: Vec<(f64, f64)> = line_data
            .iter()
            .filter(|(_, y)| *y > threshold)
            .copied()
            .collect();

        if !sig_points.is_empty() {
            chart.draw_series(sig_points.iter().map(|&(x, y)| {
                Circle::new((x, y), 2, COLOR_MAGENTA.mix(0.8).filled())
            }))?;
        }

        // Threshold dashed line
        chart.draw_series(DashedLineSeries::new(
            vec![
                (x_min - x_margin, threshold),
                (x_max + x_margin, threshold),
            ],
            5,
            3,
            COLOR_RED.mix(0.5).into(),
        ))?;
    }

    Ok(())
}
