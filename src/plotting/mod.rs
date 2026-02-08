mod delta_snp_index;
mod downsample;
mod faceted;
mod gprime;
mod manhattan;

use anyhow::Result;
use plotters::prelude::*;
use std::path::Path;
use std::sync::Once;

use crate::types::{GStatisticResult, SignificanceResult};

static FONT_INIT: Once = Once::new();

/// Register an embedded font for the ab_glyph backend (no-op after first call).
fn ensure_fonts() {
    FONT_INIT.call_once(|| {
        // DejaVu Sans is available on virtually all Linux systems
        let font_data: &'static [u8] =
            include_bytes!("/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf");
        plotters::style::register_font("sans-serif", FontStyle::Normal, font_data)
            .unwrap_or_else(|_| panic!("failed to register sans-serif Normal font"));
        plotters::style::register_font("sans-serif", FontStyle::Bold, font_data)
            .unwrap_or_else(|_| panic!("failed to register sans-serif Bold font"));
        plotters::style::register_font("sans-serif", FontStyle::Italic, font_data)
            .unwrap_or_else(|_| panic!("failed to register sans-serif Italic font"));
    });
}

/// Output format for plots.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PlotFormat {
    Png,
    Svg,
}

impl PlotFormat {
    pub fn extension(&self) -> &'static str {
        match self {
            PlotFormat::Png => "png",
            PlotFormat::Svg => "svg",
        }
    }
}

/// Configuration for plot generation.
#[derive(Debug, Clone)]
pub struct PlotConfig {
    pub width: u32,
    pub row_height: u32,
    pub format: PlotFormat,
    pub max_chromosomes: usize,
}

impl Default for PlotConfig {
    fn default() -> Self {
        Self {
            width: 1800,
            row_height: 300,
            format: PlotFormat::Png,
            max_chromosomes: 12,
        }
    }
}

// Color constants matching the Python visualization script
pub const COLOR_STEEL_BLUE: RGBColor = RGBColor(46, 134, 171);   // #2E86AB
pub const COLOR_MAGENTA: RGBColor = RGBColor(162, 59, 114);      // #A23B72
pub const COLOR_ORANGE: RGBColor = RGBColor(241, 143, 1);        // #F18F01
pub const COLOR_RED: RGBColor = RGBColor(220, 50, 50);           // threshold lines
pub const COLOR_BLACK: RGBColor = RGBColor(0, 0, 0);
pub const COLOR_GRID: RGBColor = RGBColor(200, 200, 200);

const N_COLS: usize = 3;

/// Generate a -log10(p-value) Manhattan plot from significance results.
pub fn plot_neg_log10_pvalue(
    results: &[SignificanceResult],
    path: &Path,
    config: &PlotConfig,
) -> Result<()> {
    ensure_fonts();
    let panels = faceted::prepare_chrom_panels_sig(results, config.max_chromosomes);
    if panels.is_empty() {
        anyhow::bail!("No data to plot");
    }

    let n_rows = (panels.len() + N_COLS - 1) / N_COLS;
    let height = config.row_height * n_rows as u32 + 80; // 80px for title

    match config.format {
        PlotFormat::Png => {
            let root = BitMapBackend::new(path, (config.width, height)).into_drawing_area();
            manhattan::draw_neg_log10_pvalue(&root, results, &panels)?;
            root.present()?;
        }
        PlotFormat::Svg => {
            let root = SVGBackend::new(path, (config.width, height)).into_drawing_area();
            manhattan::draw_neg_log10_pvalue(&root, results, &panels)?;
            root.present()?;
        }
    }

    eprintln!("  Plot saved to: {}", path.display());
    Ok(())
}

/// Generate a G-prime/G-stat Manhattan plot from significance results.
pub fn plot_gprime_sig(
    results: &[SignificanceResult],
    path: &Path,
    config: &PlotConfig,
    threshold_percentile: f64,
) -> Result<()> {
    ensure_fonts();
    let panels = faceted::prepare_chrom_panels_sig(results, config.max_chromosomes);
    if panels.is_empty() {
        anyhow::bail!("No data to plot");
    }

    let n_rows = (panels.len() + N_COLS - 1) / N_COLS;
    let height = config.row_height * n_rows as u32 + 80;

    match config.format {
        PlotFormat::Png => {
            let root = BitMapBackend::new(path, (config.width, height)).into_drawing_area();
            gprime::draw_gprime_sig(&root, results, &panels, threshold_percentile)?;
            root.present()?;
        }
        PlotFormat::Svg => {
            let root = SVGBackend::new(path, (config.width, height)).into_drawing_area();
            gprime::draw_gprime_sig(&root, results, &panels, threshold_percentile)?;
            root.present()?;
        }
    }

    eprintln!("  Plot saved to: {}", path.display());
    Ok(())
}

/// Generate a G-stat Manhattan plot from basic results (no significance).
pub fn plot_gstat(
    results: &[GStatisticResult],
    path: &Path,
    config: &PlotConfig,
    threshold_percentile: f64,
) -> Result<()> {
    ensure_fonts();
    let panels = faceted::prepare_chrom_panels_basic(results, config.max_chromosomes);
    if panels.is_empty() {
        anyhow::bail!("No data to plot");
    }

    let n_rows = (panels.len() + N_COLS - 1) / N_COLS;
    let height = config.row_height * n_rows as u32 + 80;

    match config.format {
        PlotFormat::Png => {
            let root = BitMapBackend::new(path, (config.width, height)).into_drawing_area();
            gprime::draw_gstat(&root, results, &panels, threshold_percentile)?;
            root.present()?;
        }
        PlotFormat::Svg => {
            let root = SVGBackend::new(path, (config.width, height)).into_drawing_area();
            gprime::draw_gstat(&root, results, &panels, threshold_percentile)?;
            root.present()?;
        }
    }

    eprintln!("  Plot saved to: {}", path.display());
    Ok(())
}

/// Generate a delta SNP-index plot from significance results.
pub fn plot_delta_snp_index_sig(
    results: &[SignificanceResult],
    path: &Path,
    config: &PlotConfig,
) -> Result<()> {
    ensure_fonts();
    let panels = faceted::prepare_chrom_panels_sig(results, config.max_chromosomes);
    if panels.is_empty() {
        anyhow::bail!("No data to plot");
    }

    let n_rows = (panels.len() + N_COLS - 1) / N_COLS;
    let height = config.row_height * n_rows as u32 + 80;

    match config.format {
        PlotFormat::Png => {
            let root = BitMapBackend::new(path, (config.width, height)).into_drawing_area();
            delta_snp_index::draw_delta_snp_sig(&root, results, &panels)?;
            root.present()?;
        }
        PlotFormat::Svg => {
            let root = SVGBackend::new(path, (config.width, height)).into_drawing_area();
            delta_snp_index::draw_delta_snp_sig(&root, results, &panels)?;
            root.present()?;
        }
    }

    eprintln!("  Plot saved to: {}", path.display());
    Ok(())
}

/// Generate a delta SNP-index plot from basic results.
pub fn plot_delta_snp_index_basic(
    results: &[GStatisticResult],
    path: &Path,
    config: &PlotConfig,
) -> Result<()> {
    ensure_fonts();
    let panels = faceted::prepare_chrom_panels_basic(results, config.max_chromosomes);
    if panels.is_empty() {
        anyhow::bail!("No data to plot");
    }

    let n_rows = (panels.len() + N_COLS - 1) / N_COLS;
    let height = config.row_height * n_rows as u32 + 80;

    match config.format {
        PlotFormat::Png => {
            let root = BitMapBackend::new(path, (config.width, height)).into_drawing_area();
            delta_snp_index::draw_delta_snp_basic(&root, results, &panels)?;
            root.present()?;
        }
        PlotFormat::Svg => {
            let root = SVGBackend::new(path, (config.width, height)).into_drawing_area();
            delta_snp_index::draw_delta_snp_basic(&root, results, &panels)?;
            root.present()?;
        }
    }

    eprintln!("  Plot saved to: {}", path.display());
    Ok(())
}
