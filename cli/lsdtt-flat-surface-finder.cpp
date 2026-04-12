// lsdtt-flat-surface-finder.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html
//
// -----------------------------------------------------------------------------
// lsdtt-flat-surface-finder
//
// Identifies broadly flat surfaces (floodplains, alluvial fans, fluvial and
// agricultural terraces) in a DEM using surface texture metrics.
//
// Metrics:
//   MADk2         — second-order roughness; slope-independent; recommended
//   Residual IQR  — IQR of elevation residuals to a fitted reference surface;
//                   slope-independent by construction; three fitting methods:
//                     fftw_highpass      (DEM minus FFT lowpass surface)
//                     wavelet_quadratic  (DEM minus local wavelet polynomial fit)
//                     moving_window_mean (DEM minus local moving-window mean)
//   NVD s3        — smallest eigenvalue of surface normal orientation tensor;
//                   requires FFTW; disabled by default
//
// Outputs (always written):
//   <stem>_labelled_flat_mask.tif      — UInt32 labelled flat components
//   <stem>_madk2_qq.csv                — QQ data for MADk2 (if enabled)
//   <stem>_iqr_qq.csv                  — QQ data for residual IQR (if enabled)
//   <stem>_nvd_qq.csv                  — QQ data for NVD (if enabled)
//
// Optional outputs:
//   <stem>_madk2_roughness.tif
//   <stem>_iqr_residuals.tif
//   <stem>_iqr_texture.tif
//   <stem>_nvd_s3.tif
//   <stem>_madk2_flat_mask.tif
//   <stem>_iqr_flat_mask.tif
//   <stem>_nvd_flat_mask.tif
//   <stem>_combined_flat_mask.tif
//
// Usage:
//   lsdtt-flat-surface-finder <path> <parameter_file>
//   lsdtt-flat-surface-finder <parameter_file>
//   lsdtt-flat-surface-finder -h
// -----------------------------------------------------------------------------

#include <cmath>
#include <cstdlib>
#include <initializer_list>
#include <iostream>
#include <limits>
#include <map>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "../src/core/lsdtt_data_loader.hpp"
#include "../src/core/lsdtt_param_parser.hpp"
#include "../src/cpu/lsdtt_flat_surface_finder.hpp"
#include "../src/cpu/lsdtt_raster_cpu.hpp"
#include "../src/utils/lsdtt_distributions.hpp"
#include "../src/utils/lsdtt_histogram.hpp"

static const std::string k_version  = "0.1.0";
static const std::string k_citation = "https://github.com/LSDtopotools/lsdtt3";

static void print_banner()
{
    std::cout
        << "=========================================================\n"
        << "|| Welcome to lsdtt-flat-surface-finder                 ||\n"
        << "|| I will identify flat surfaces in a DEM for you       ||\n"
        << "|| Part of lsdtt3 by Simon M. Mudd                      ||\n"
        << "=========================================================\n"
        << "|| Version : " << k_version  << "\n"
        << "|| Citation: " << k_citation << "\n"
        << "=========================================================\n\n";
}

static std::pair<std::string, std::string> parse_args(int argc, char* argv[])
{
    auto is_flag = [](const std::string& s, std::initializer_list<std::string> flags)
    {
        for (const auto& f : flags) if (s == f) return true;
        return false;
    };
    for (int i = 1; i < argc; ++i)
    {
        std::string a = argv[i];
        if (is_flag(a, {"-h","--help","help","-H"}))
        {
            std::cout << "Usage: lsdtt-flat-surface-finder <path> <parameter_file>\n"
                      << "       lsdtt-flat-surface-finder <parameter_file>\n"
                      << "       lsdtt-flat-surface-finder -h\n";
            std::exit(0);
        }
        if (is_flag(a, {"-v","--version","version"}))
        {
            std::cout << "lsdtt-flat-surface-finder version " << k_version << "\n";
            std::exit(0);
        }
    }
    std::string path = "./";
    std::string fname;
    if (argc == 1) { std::cerr << "No parameter file. Use -h for help.\n"; std::exit(1); }
    else if (argc == 2) { fname = argv[1]; }
    else { path = argv[1]; fname = argv[2]; }
    return {path, fname};
}

static void write_qq_csv(
    lsdtt::lsdtt_data_loader&         loader,
    const std::string&                out_path,
    const std::vector<lsdtt::qq_bin>& qq)
{
    std::vector<std::vector<std::string>> rows;
    rows.reserve(qq.size());
    for (const auto& b : qq)
        rows.push_back({ std::to_string(b.standard_normal_quantile),
                         std::to_string(b.data_quantile),
                         std::to_string(b.normal_reference) });
    loader.write_csv(out_path,
                     {"standard_normal_variate", "data_quantile", "normal_reference"},
                     rows);
}

// Write one CSV summarising all four distribution fits for a metric.
// Columns: family, param1, param2, ks_statistic, max_divergence_x, n_values,
//          ks_critical_value_05, valid
//
// ks_critical_value_05 = 1.36 / sqrt(n_values) is the Kolmogorov-Smirnov
// critical value at the 5 % significance level for a *known-parameter* test.
// Because parameters are estimated from the data (MLE), the actual critical
// value is somewhat smaller (Lilliefors correction), so this column gives
// a conservative upper bound: if ks_statistic > ks_critical_value_05 the
// distribution is definitively rejected; if smaller, it may still be a poor
// fit due to the MLE bias.
static void write_dist_fits_csv(
    lsdtt::lsdtt_data_loader&                          loader,
    const std::string&                                 out_path,
    const std::vector<lsdtt::lsdtt_distribution_fit>&  fits)
{
    std::vector<std::vector<std::string>> rows;
    rows.reserve(fits.size());
    for (const auto& f : fits)
    {
        const double ks_crit = (f.n_values > 0)
            ? 1.36 / std::sqrt(static_cast<double>(f.n_values))
            : std::numeric_limits<double>::quiet_NaN();

        rows.push_back({
            lsdtt::distribution_family_name(f.family),
            std::to_string(f.param1),
            std::to_string(f.param2),
            std::to_string(f.ks_statistic),
            std::to_string(f.max_divergence_x),
            std::to_string(f.n_values),
            std::to_string(ks_crit),
            f.valid ? "true" : "false"
        });
    }
    loader.write_csv(out_path,
        {"family","param1","param2","ks_statistic","max_divergence_x",
         "n_values","ks_critical_value_05","valid"},
        rows);
}

// Write a smoothed PDF CSV with three columns:
//   bin_centre, smoothed_pdf, fitted_pdf
// where fitted_pdf is the best-fitting distribution's PDF evaluated at each
// bin centre.  The best fit is the valid entry in `fits` with the lowest
// KS statistic.
static void write_smoothed_pdf_csv(
    lsdtt::lsdtt_data_loader&                          loader,
    const std::string&                                 out_path,
    const lsdtt::lsdtt_smoothed_pdf&                   spdf,
    const std::vector<lsdtt::lsdtt_distribution_fit>&  fits)
{
    if (spdf.empty()) return;

    // Find the best valid fit.
    const lsdtt::lsdtt_distribution_fit* best = nullptr;
    for (const auto& f : fits)
    {
        if (!f.valid || !std::isfinite(f.ks_statistic)) continue;
        if (!best || f.ks_statistic < best->ks_statistic)
            best = &f;
    }

    std::vector<std::vector<std::string>> rows;
    rows.reserve(spdf.bin_centres.size());
    for (std::size_t i = 0; i < spdf.bin_centres.size(); ++i)
    {
        const double x   = spdf.bin_centres[i];
        const double pdf = spdf.smoothed_pdf[i];
        const double fit_pdf = best ? lsdtt::eval_pdf(*best, x)
                                    : std::numeric_limits<double>::quiet_NaN();
        rows.push_back({ std::to_string(x),
                         std::to_string(pdf),
                         std::to_string(fit_pdf) });
    }

    // Include the best-fit family name in the header comment via the column
    // name so the user knows which distribution is plotted.
    const std::string fit_col = best
        ? std::string("fitted_pdf_") + lsdtt::distribution_family_name(best->family)
        : "fitted_pdf";

    loader.write_csv(out_path,
                     {"bin_centre", "smoothed_pdf", fit_col},
                     rows);
}

// =============================================================================
// main
// =============================================================================
int main(int argc, char* argv[])
{
    print_banner();
    auto [path, fname] = parse_args(argc, argv);

    lsdtt::lsdtt_param_parser pp;
    try { pp = lsdtt::lsdtt_param_parser(path, fname); }
    catch (const std::exception& e)
    { std::cerr << "Error loading parameter file: " << e.what() << "\n"; return 1; }

    // -------------------------------------------------------------------------
    // Parameter declarations
    // -------------------------------------------------------------------------
    std::map<std::string, bool>                     bool_default_map;
    std::map<std::string, float>                    float_default_map;
    std::map<std::string, int>                      int_default_map;
    std::map<std::string, std::string>              string_default_map;
    std::map<std::string, std::vector<std::string>> help_map;

    int_default_map["n_threads"] = 0;
    help_map["n_threads"] = {"int","0","OpenMP thread count.","0 = all hardware threads."};

    string_default_map["write_output_format"] = "cog";
    help_map["write_output_format"] = {"string","cog","Raster output format.",
        "cog, gtiff, or same_as_input."};

    // --- MADk2 ---
    bool_default_map["use_madk2"] = true;
    help_map["use_madk2"] = {"bool","true",
        "Enable MADk2 second-order roughness metric.",
        "Recommended for all landscape types including sloping flat surfaces. "
        "Second-order differences cancel the linear slope trend."};

    float_default_map["madk2_window_radius_metres"] = 100.0f;
    help_map["madk2_window_radius_metres"] = {"float","100",
        "Moving-window radius (m) for MADk2.",
        "Choose smaller than the width of the target flat surfaces."};

    float_default_map["madk2_threshold_percentile"] = 25.0f;
    help_map["madk2_threshold_percentile"] = {"float","25",
        "Percentile of the MADk2 distribution below which a pixel is classed flat.",
        "Inspect the QQ CSV to identify the distribution break-point.  Range: 0–100."};

    // --- Residual IQR ---
    bool_default_map["use_residual_iqr"] = true;
    help_map["use_residual_iqr"] = {"bool","true",
        "Enable residual IQR metric (IQR of elevation residuals to a fitted surface).",
        "Slope-independent by construction.  Fitting method is set by iqr_fit_method."};

    string_default_map["iqr_fit_method"] = "moving_window_mean";
    help_map["iqr_fit_method"] = {"string","moving_window_mean",
        "Method used to compute the reference surface for residual IQR.",
        "fftw_highpass     — residuals = fftw2D_highpass_filter(dem, iqr_fit_lambda_metres). "
        "  Global spectral separation; cleanest frequency response; requires FFTW. "
        "wavelet_quadratic — residuals = dem − wavelet_quadratic_fit(dem, iqr_fit_lambda_metres).f. "
        "  Local polynomial detrending; handles spatially varying background; requires FFTW. "
        "moving_window_mean — residuals = dem − moving_window_mean(dem, iqr_detrend_radius_metres). "
        "  No FFTW dependency; simpler but less frequency-selective."};

    float_default_map["iqr_fit_lambda_metres"] = 500.0f;
    help_map["iqr_fit_lambda_metres"] = {"float","500",
        "Reference-surface wavelength (m) for fftw_highpass and wavelet_quadratic methods.",
        "Set to a value longer than the widest flat surface to be detected so that "
        "the fitted surface captures the regional trend.  Ignored for moving_window_mean."};

    float_default_map["iqr_detrend_radius_metres"] = 500.0f;
    help_map["iqr_detrend_radius_metres"] = {"float","500",
        "Moving-window mean radius (m) for the moving_window_mean detrending method.",
        "Ignored for fftw_highpass and wavelet_quadratic."};

    float_default_map["iqr_window_radius_metres"] = 50.0f;
    help_map["iqr_window_radius_metres"] = {"float","50",
        "Moving-window radius (m) for computing IQR of residuals.",
        "Independent of the detrending scale.  Should be smaller than the "
        "width of the flat surfaces to be detected."};

    string_default_map["iqr_window_shape"] = "circular";
    help_map["iqr_window_shape"] = {"string","circular",
        "Window shape for the residual IQR computation.","circular or square."};

    float_default_map["iqr_threshold_percentile"] = 25.0f;
    help_map["iqr_threshold_percentile"] = {"float","25",
        "Percentile of the residual IQR distribution below which a pixel is classed flat.",
        "Range: 0–100."};

    // --- NVD ---
    bool_default_map["use_nvd"] = false;
    help_map["use_nvd"] = {"bool","false",
        "Enable NVD (Normal Vector Dispersion s3) metric.",
        "Requires FFTW.  Low s3 = co-planar normals = flat surface."};

    float_default_map["nvd_polyfit_lambda"] = 10.0f;
    help_map["nvd_polyfit_lambda"] = {"float","10",
        "Wavelet wavelength (m) for deriving surface normals (d/e coefficients).",
        "Only used when use_nvd is true."};

    float_default_map["nvd_window_radius_metres"] = 100.0f;
    help_map["nvd_window_radius_metres"] = {"float","100",
        "Moving-window radius (m) for NVD orientation tensor accumulation.",
        "Only used when use_nvd is true."};

    string_default_map["nvd_window_shape"] = "circular";
    help_map["nvd_window_shape"] = {"string","circular",
        "Window shape for NVD computation.","circular or square."};

    float_default_map["nvd_threshold_percentile"] = 25.0f;
    help_map["nvd_threshold_percentile"] = {"float","25",
        "Percentile of the NVD s3 distribution below which a pixel is classed flat.",
        "Range: 0–100."};

    // --- Combination ---
    bool_default_map["require_all_metrics_pass"] = true;
    help_map["require_all_metrics_pass"] = {"bool","true",
        "Require a pixel to pass ALL enabled metric thresholds.",
        "true (AND, conservative) or false (OR, permissive)."};

    int_default_map["min_component_pixels"] = 100;
    help_map["min_component_pixels"] = {"int","100",
        "Minimum 8-connected component size for inclusion in the labelled flat mask.",
        "Smaller components are excluded.  Set to 1 to include all."};

    int_default_map["qq_n_bins"] = 100;
    help_map["qq_n_bins"] = {"int","100",
        "Number of quantile bins in the QQ output CSVs.","100 = one row per percentile."};

    int_default_map["pdf_n_bins"] = 100;
    help_map["pdf_n_bins"] = {"int","100",
        "Number of histogram bins for the smoothed PDF and peak detection.",
        "Increase for finer resolution; decrease to suppress noise."};

    int_default_map["pdf_smooth_half_window"] = 3;
    help_map["pdf_smooth_half_window"] = {"int","3",
        "Half-width of the box-smoothing kernel (bins) used for PDF peak detection.",
        "Total kernel width = 2*N+1.  Increase to merge spurious peaks."};

    int_default_map["curvature_window_pixels"] = 3;
    help_map["curvature_window_pixels"] = {"int","3",
        "Square window size (pixels) for total curvature (Laplacian) computation.",
        "Must be odd and >= 3.  3 = 3x3 neighbourhood (standard 5-point stencil)."};

    int_default_map["relief_window_pixels"] = 5;
    help_map["relief_window_pixels"] = {"int","5",
        "Square window size (pixels) for local relief (max - min) computation.",
        "Must be odd and >= 3.  5 = 5x5 neighbourhood."};

    bool_default_map["use_relief"] = true;
    help_map["use_relief"] = {"bool","true",
        "Use local relief as a flat-surface masking metric.",
        "Pixels below relief_threshold_percentile are classed flat."};

    float_default_map["relief_threshold_percentile"] = 10.0f;
    help_map["relief_threshold_percentile"] = {"float","10",
        "Percentile of the local relief distribution below which a pixel is classed flat.",
        "Range: 0-100.  Lower values = stricter (only lowest-relief pixels pass)."};


    // --- Output toggles ---
    bool_default_map["write_texture_rasters"] = false;
    help_map["write_texture_rasters"] = {"bool","false",
        "Write per-metric texture rasters (MADk2 roughness, IQR residuals, "
        "IQR texture, NVD s3).",
        "Useful for threshold selection and visual inspection."};

    bool_default_map["write_metric_flat_masks"] = false;
    help_map["write_metric_flat_masks"] = {"bool","false",
        "Write binary flat masks for each individual enabled metric.",
        "Useful for diagnosing which metric drives the combined mask."};

    bool_default_map["write_combined_flat_mask"] = true;
    help_map["write_combined_flat_mask"] = {"bool","true",
        "Write the combined binary flat mask.",
        "1 = flat, 0 = not flat."};

    // -------------------------------------------------------------------------
    // Parse
    // -------------------------------------------------------------------------
    pp.parse_all_parameters(float_default_map, int_default_map,
                            bool_default_map, string_default_map);
    auto bf = pp.get_bool_parameters();
    auto ff = pp.get_float_parameters();
    auto ip = pp.get_int_parameters();
    auto sf = pp.get_string_parameters();

    // -------------------------------------------------------------------------
    // Cross-parameter resolution
    // -------------------------------------------------------------------------
    lsdtt::raster_output_format out_format = lsdtt::raster_output_format::cog;
    {
        const std::string& fmt = sf["write_output_format"];
        if (fmt == "gtiff" || fmt == "geotiff")
            out_format = lsdtt::raster_output_format::gtiff;
        else if (fmt == "same_as_input")
            out_format = lsdtt::raster_output_format::same_as_input;
    }

    auto resolve_window_shape = [](const std::string& s)
    {
        return (s == "square") ? lsdtt::window_shape::square
                               : lsdtt::window_shape::circular;
    };

    lsdtt::iqr_fit_method iqr_method = lsdtt::iqr_fit_method::fftw_highpass;
    {
        const std::string& m = sf["iqr_fit_method"];
        if      (m == "wavelet_quadratic")  iqr_method = lsdtt::iqr_fit_method::wavelet_quadratic;
        else if (m == "moving_window_mean") iqr_method = lsdtt::iqr_fit_method::moving_window_mean;
        else if (m == "fftw_highpass")      iqr_method = lsdtt::iqr_fit_method::fftw_highpass;
        else
        {
            std::cerr << "Error: unknown iqr_fit_method '" << m
                      << "'. Allowed: fftw_highpass, wavelet_quadratic, moving_window_mean\n";
            return 1;
        }
    }

    if (!bf["use_madk2"] && !bf["use_residual_iqr"] && !bf["use_nvd"])
    {
        std::cerr << "Error: at least one metric must be enabled.\n";
        return 1;
    }

    // -------------------------------------------------------------------------
    // Load DEM
    // -------------------------------------------------------------------------
    lsdtt::lsdtt_data_loader loader;
    lsdtt::lsdtt_raster_cpu  dem;
    {
        const std::string in_path = pp.get_full_read_path();
        std::cout << "Loading DEM: " << in_path << "\n";
        try { loader.open_raster(in_path); }
        catch (const std::exception& e)
        { std::cerr << "Error loading DEM: " << e.what() << "\n"; return 1; }
        dem = lsdtt::lsdtt_raster_cpu(loader);
        try { dem.get_meta().check_projected_metres(); }
        catch (const std::exception& e)
        { std::cerr << "Projection error: " << e.what() << "\n"; return 1; }
        std::cout << "  rows=" << dem.get_meta().n_rows
                  << "  cols=" << dem.get_meta().n_cols
                  << "  cell=" << dem.get_meta().cell_size << " m\n\n";
    }

    // -------------------------------------------------------------------------
    // Build request
    // -------------------------------------------------------------------------
    lsdtt::flat_surface_request req;
    req.use_madk2                  = bf["use_madk2"];
    req.madk2_window_radius_metres = static_cast<double>(ff["madk2_window_radius_metres"]);
    req.madk2_threshold_percentile = static_cast<double>(ff["madk2_threshold_percentile"]);
    req.use_residual_iqr           = bf["use_residual_iqr"];
    req.iqr_method                 = iqr_method;
    req.iqr_fit_lambda_metres      = static_cast<double>(ff["iqr_fit_lambda_metres"]);
    req.iqr_detrend_radius_metres  = static_cast<double>(ff["iqr_detrend_radius_metres"]);
    req.iqr_window_radius_metres   = static_cast<double>(ff["iqr_window_radius_metres"]);
    req.iqr_window_shape           = resolve_window_shape(sf["iqr_window_shape"]);
    req.iqr_threshold_percentile   = static_cast<double>(ff["iqr_threshold_percentile"]);
    req.use_nvd                    = bf["use_nvd"];
    req.nvd_polyfit_lambda         = static_cast<double>(ff["nvd_polyfit_lambda"]);
    req.nvd_window_radius_metres   = static_cast<double>(ff["nvd_window_radius_metres"]);
    req.nvd_window_shape           = resolve_window_shape(sf["nvd_window_shape"]);
    req.nvd_threshold_percentile   = static_cast<double>(ff["nvd_threshold_percentile"]);
    req.require_all_metrics_pass   = bf["require_all_metrics_pass"];
    req.min_component_pixels       = ip["min_component_pixels"];
    req.qq_n_bins                  = ip["qq_n_bins"];
    req.pdf_n_bins                 = ip["pdf_n_bins"];
    req.pdf_smooth_half_window     = ip["pdf_smooth_half_window"];
    req.curvature_window_pixels      = ip["curvature_window_pixels"];
    req.relief_window_pixels         = ip["relief_window_pixels"];
    req.use_relief                   = bf["use_relief"];
    req.relief_threshold_percentile  = static_cast<double>(ff["relief_threshold_percentile"]);
    req.n_threads                    = ip["n_threads"];

    // -------------------------------------------------------------------------
    // Run
    // -------------------------------------------------------------------------
    std::cout << "Running flat surface finder...\n";
    if (req.use_madk2)
        std::cout << "  MADk2        radius=" << req.madk2_window_radius_metres
                  << " m  threshold=" << req.madk2_threshold_percentile << "th pct\n";
    if (req.use_residual_iqr)
    {
        std::cout << "  Residual IQR method=" << sf["iqr_fit_method"];
        if (iqr_method == lsdtt::iqr_fit_method::moving_window_mean)
            std::cout << "  detrend_radius=" << req.iqr_detrend_radius_metres << " m";
        else
            std::cout << "  fit_lambda=" << req.iqr_fit_lambda_metres << " m";
        std::cout << "  iqr_radius=" << req.iqr_window_radius_metres
                  << " m  threshold=" << req.iqr_threshold_percentile << "th pct\n";
    }
    if (req.use_nvd)
        std::cout << "  NVD s3       radius=" << req.nvd_window_radius_metres
                  << " m  threshold=" << req.nvd_threshold_percentile << "th pct\n";
    if (req.use_relief)
        std::cout << "  Local relief  window=" << req.relief_window_pixels
                  << " px  threshold=" << req.relief_threshold_percentile << "th pct\n";
    std::cout << "  combination: " << (req.require_all_metrics_pass ? "AND" : "OR")
              << "  min_component=" << req.min_component_pixels << " px\n\n";

    lsdtt::flat_surface_result result;
    try { result = lsdtt::find_flat_surfaces(dem, req); }
    catch (const std::exception& e)
    { std::cerr << "Error: " << e.what() << "\n"; return 1; }

    // -------------------------------------------------------------------------
    // Write outputs
    // -------------------------------------------------------------------------
    const std::string stem = pp.get_full_write_path();

    if (bf["write_texture_rasters"])
    {
        if (req.use_madk2)
        {
            const std::string p = stem + "_madk2_roughness.tif";
            std::cout << "Writing MADk2 roughness: " << p << "\n";
            loader.write_full(result.madk2_roughness.get_data(),
                              result.madk2_roughness.get_meta(), p, out_format);
        }
        if (req.use_residual_iqr)
        {
            {
                const std::string p = stem + "_iqr_residuals.tif";
                std::cout << "Writing IQR residuals: " << p << "\n";
                loader.write_full(result.iqr_residuals.get_data(),
                                  result.iqr_residuals.get_meta(), p, out_format);
            }
            {
                const std::string p = stem + "_iqr_texture.tif";
                std::cout << "Writing IQR texture: " << p << "\n";
                loader.write_full(result.iqr_texture.get_data(),
                                  result.iqr_texture.get_meta(), p, out_format);
            }
        }
        if (req.use_nvd)
        {
            const std::string p = stem + "_nvd_s3.tif";
            std::cout << "Writing NVD s3: " << p << "\n";
            loader.write_full(result.nvd_s3.get_data(),
                              result.nvd_s3.get_meta(), p, out_format);
        }

        // Diagnostic rasters — always computed, written here
        {
            const std::string p = stem + "_curvature.tif";
            std::cout << "Writing total curvature: " << p << "\n";
            loader.write_full(result.curvature.get_data(),
                              result.curvature.get_meta(), p, out_format);
        }
        {
            const std::string p = stem + "_relief.tif";
            std::cout << "Writing local relief: " << p << "\n";
            loader.write_full(result.relief.get_data(),
                              result.relief.get_meta(), p, out_format);
        }
    }

    if (bf["write_metric_flat_masks"])
    {
        if (req.use_madk2)
        {
            const std::string p = stem + "_madk2_flat_mask.tif";
            std::cout << "Writing MADk2 flat mask: " << p << "\n";
            loader.write_full(result.madk2_flat_mask.get_data(),
                              result.madk2_flat_mask.get_meta(), p, out_format);
        }
        if (req.use_residual_iqr)
        {
            const std::string p = stem + "_iqr_flat_mask.tif";
            std::cout << "Writing residual IQR flat mask: " << p << "\n";
            loader.write_full(result.iqr_flat_mask.get_data(),
                              result.iqr_flat_mask.get_meta(), p, out_format);
        }
        if (req.use_nvd)
        {
            const std::string p = stem + "_nvd_flat_mask.tif";
            std::cout << "Writing NVD flat mask: " << p << "\n";
            loader.write_full(result.nvd_flat_mask.get_data(),
                              result.nvd_flat_mask.get_meta(), p, out_format);
        }
        if (req.use_relief)
        {
            const std::string p = stem + "_relief_flat_mask.tif";
            std::cout << "Writing relief flat mask: " << p << "\n";
            loader.write_full(result.relief_flat_mask.get_data(),
                              result.relief_flat_mask.get_meta(), p, out_format);
        }
    }

    if (bf["write_combined_flat_mask"])
    {
        const std::string p = stem + "_combined_flat_mask.tif";
        std::cout << "Writing combined flat mask: " << p << "\n";
        loader.write_full(result.combined_flat_mask.get_data(),
                          result.combined_flat_mask.get_meta(), p, out_format);
    }

    // Labelled mask — always written
    {
        const std::string p = stem + "_labelled_flat_mask.tif";
        std::cout << "Writing labelled flat mask (UInt32): " << p << "\n";
        loader.write_full(result.labelled_flat_mask.get_data(),
                          result.labelled_flat_mask.get_meta(), p, out_format,
                          GDT_UInt32);
    }

    // QQ CSVs — always written for each enabled metric
    if (req.use_madk2 && !result.madk2_qq.empty())
    {
        const std::string p = stem + "_madk2_qq.csv";
        std::cout << "Writing MADk2 QQ: " << p << "\n";
        write_qq_csv(loader, p, result.madk2_qq);
    }
    if (req.use_residual_iqr && !result.iqr_qq.empty())
    {
        const std::string p = stem + "_iqr_qq.csv";
        std::cout << "Writing residual IQR QQ: " << p << "\n";
        write_qq_csv(loader, p, result.iqr_qq);
    }
    if (req.use_nvd && !result.nvd_qq.empty())
    {
        const std::string p = stem + "_nvd_qq.csv";
        std::cout << "Writing NVD QQ: " << p << "\n";
        write_qq_csv(loader, p, result.nvd_qq);
    }
    if (!result.curvature_qq.empty())
    {
        const std::string p = stem + "_curvature_qq.csv";
        std::cout << "Writing curvature QQ: " << p << "\n";
        write_qq_csv(loader, p, result.curvature_qq);
    }
    if (!result.relief_qq.empty())
    {
        const std::string p = stem + "_relief_qq.csv";
        std::cout << "Writing relief QQ: " << p << "\n";
        write_qq_csv(loader, p, result.relief_qq);
    }

    // Smoothed PDF CSVs — written for each metric
    // Three columns: bin_centre, smoothed_pdf, fitted_pdf_<family>
    if (req.use_madk2 && !result.madk2_smoothed_pdf.empty())
    {
        const std::string p = stem + "_madk2_smoothed_pdf.csv";
        std::cout << "Writing MADk2 smoothed PDF: " << p << "\n";
        write_smoothed_pdf_csv(loader, p,
                               result.madk2_smoothed_pdf, result.madk2_dist_fits);
    }
    if (req.use_residual_iqr && !result.iqr_smoothed_pdf.empty())
    {
        const std::string p = stem + "_iqr_smoothed_pdf.csv";
        std::cout << "Writing residual IQR smoothed PDF: " << p << "\n";
        write_smoothed_pdf_csv(loader, p,
                               result.iqr_smoothed_pdf, result.iqr_dist_fits);
    }
    if (req.use_nvd && !result.nvd_smoothed_pdf.empty())
    {
        const std::string p = stem + "_nvd_smoothed_pdf.csv";
        std::cout << "Writing NVD smoothed PDF: " << p << "\n";
        write_smoothed_pdf_csv(loader, p,
                               result.nvd_smoothed_pdf, result.nvd_dist_fits);
    }
    if (!result.curvature_smoothed_pdf.empty())
    {
        const std::string p = stem + "_curvature_smoothed_pdf.csv";
        std::cout << "Writing curvature smoothed PDF: " << p << "\n";
        write_smoothed_pdf_csv(loader, p,
                               result.curvature_smoothed_pdf, result.curvature_dist_fits);
    }
    if (!result.relief_smoothed_pdf.empty())
    {
        const std::string p = stem + "_relief_smoothed_pdf.csv";
        std::cout << "Writing relief smoothed PDF: " << p << "\n";
        write_smoothed_pdf_csv(loader, p,
                               result.relief_smoothed_pdf, result.relief_dist_fits);
    }

    // Distribution fit CSVs — always written for each metric
    auto write_dist = [&](const std::string& label,
                          const std::vector<lsdtt::lsdtt_distribution_fit>& fits)
    {
        if (fits.empty()) return;
        const std::string p = stem + "_" + label + "_dist_fits.csv";
        std::cout << "Writing " << label << " distribution fits: " << p << "\n";
        write_dist_fits_csv(loader, p, fits);
    };

    if (req.use_madk2)        write_dist("madk2",     result.madk2_dist_fits);
    if (req.use_residual_iqr) write_dist("iqr",       result.iqr_dist_fits);
    if (req.use_nvd)          write_dist("nvd",       result.nvd_dist_fits);
    write_dist("curvature", result.curvature_dist_fits);
    write_dist("relief",    result.relief_dist_fits);

    std::cout << "\n";
    pp.print_parameters();
    std::cout << "\nlsdtt-flat-surface-finder complete.\n";
    return 0;
}
