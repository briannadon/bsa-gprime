use anyhow::Result;
use plotters::prelude::*;

use super::downsample;
use super::faceted::ChromPanel;
use super::{COLOR_BLACK, COLOR_GRID, COLOR_ORANGE, N_COLS};
use crate::types::{GStatisticResult, SignificanceResult};

/// Draw delta SNP-index plot from significance results.
pub fn draw_delta_snp_sig<DB: DrawingBackend>(
    root: &DrawingArea<DB, plotters::coord::Shift>,
    results: &[SignificanceResult],
    panels: &[ChromPanel],
) -> Result<()>
where
    DB::ErrorType: 'static,
{
    root.fill(&WHITE)?;

    let (title_area, chart_area) = root.split_vertically(60);
    title_area.titled(
        "BSA - Delta SNP Index across Genome",
        ("sans-serif", 22).into_font().color(&BLACK),
    )?;

    let n_rows = (panels.len() + N_COLS - 1) / N_COLS;
    let panel_areas = chart_area.split_evenly((n_rows, N_COLS));

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
            "Delta(SNP-index)"
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
                -1.1f64..1.1f64,
            )?;

        chart
            .configure_mesh()
            .x_desc("Position (Mb)")
            .y_desc(y_label)
            .x_label_style(("sans-serif", 10))
            .y_label_style(("sans-serif", 10))
            .light_line_style(COLOR_GRID.mix(0.3))
            .draw()?;

        // Line plot of delta SNP index
        let line_data: Vec<(f64, f64)> = panel
            .indices
            .iter()
            .enumerate()
            .map(|(j, &i)| (positions_mb[j], results[i].raw.delta_snp_index))
            .collect();

        let line_ds = downsample::minmax_downsample(&line_data, pixel_width);
        let line_draw = line_ds.as_deref().unwrap_or(&line_data);

        chart.draw_series(LineSeries::new(
            line_draw.iter().copied(),
            COLOR_ORANGE.mix(0.7).stroke_width(1),
        ))?;

        // Zero line
        chart.draw_series(LineSeries::new(
            vec![(x_min - x_margin, 0.0), (x_max + x_margin, 0.0)],
            COLOR_BLACK.mix(0.3).stroke_width(1),
        ))?;
    }

    Ok(())
}

/// Draw delta SNP-index plot from basic results.
pub fn draw_delta_snp_basic<DB: DrawingBackend>(
    root: &DrawingArea<DB, plotters::coord::Shift>,
    results: &[GStatisticResult],
    panels: &[ChromPanel],
) -> Result<()>
where
    DB::ErrorType: 'static,
{
    root.fill(&WHITE)?;

    let (title_area, chart_area) = root.split_vertically(60);
    title_area.titled(
        "BSA - Delta SNP Index across Genome",
        ("sans-serif", 22).into_font().color(&BLACK),
    )?;

    let n_rows = (panels.len() + N_COLS - 1) / N_COLS;
    let panel_areas = chart_area.split_evenly((n_rows, N_COLS));

    let pixel_width = downsample::panel_pixel_width(root.dim_in_pixel().0, N_COLS);

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
            "Delta(SNP-index)"
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
                -1.1f64..1.1f64,
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
            .map(|(j, &i)| (positions_mb[j], results[i].delta_snp_index))
            .collect();

        let line_ds = downsample::minmax_downsample(&line_data, pixel_width);
        let line_draw = line_ds.as_deref().unwrap_or(&line_data);

        chart.draw_series(LineSeries::new(
            line_draw.iter().copied(),
            COLOR_ORANGE.mix(0.7).stroke_width(1),
        ))?;

        // Zero line
        chart.draw_series(LineSeries::new(
            vec![(x_min - x_margin, 0.0), (x_max + x_margin, 0.0)],
            COLOR_BLACK.mix(0.3).stroke_width(1),
        ))?;
    }

    Ok(())
}
