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
// agricultural terraces) in a DEM using surface texture metrics.  "Flat" means
// low surface roughness/texture — not necessarily level in terms of regional
// slope.
//
// Metrics:
//   MADk2   — second-order roughness (effective for sloping smooth surfaces)
//   IQR     — moving-window interquartile range of raw elevations
//   NVD s3  — smallest eigenvalue of the surface normal orientation tensor
//             (requires FFTW; disabled by default)
//
// Threshold selection:
//   Each metric is thresholded at a user-specified percentile.  The
//   quantile-quantile (QQ) CSV outputs compare empirical metric quantiles
//   against a standard normal variate and help the user identify the
//   appropriate threshold for their landscape.
//
// Outputs:
//   <stem>_madk2_roughness.tif         — MADk2 omnidirectional roughness
//   <stem>_window_iqr.tif              — moving-window IQR
//   <stem>_nvd_s3.tif                  — NVD smallest eigenvalue
//   <stem>_madk2_flat_mask.tif         — binary: 1 = flat by MADk2
//   <stem>_iqr_flat_mask.tif           — binary: 1 = flat by IQR
//   <stem>_nvd_flat_mask.tif           — binary: 1 = flat by NVD
//   <stem>_combined_flat_mask.tif      — combined binary flat mask
//   <stem>_labelled_flat_mask.tif      — UInt32 labelled components
//   <stem>_madk2_qq.csv                — QQ data for MADk2
//   <stem>_iqr_qq.csv                  — QQ data for IQR
//   <stem>_nvd_qq.csv                  — QQ data for NVD
//
// Usage:
//   lsdtt-flat-surface-finder <path> <parameter_file>
//   lsdtt-flat-surface-finder <parameter_file>   (path assumed to be ./)
//   lsdtt-flat-surface-finder -h                 (print help and exit)
// -----------------------------------------------------------------------------

#include <cstdlib>
#include <initializer_list>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "../src/core/lsdtt_data_loader.hpp"
#include "../src/core/lsdtt_param_parser.hpp"
#include "../src/cpu/lsdtt_flat_surface_finder.hpp"
#include "../src/cpu/lsdtt_raster_cpu.hpp"

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
    auto is_flag = [](const std::string& s,
                      std::initializer_list<std::string> flags)
    {
        for (const auto& f : flags) if (s == f) return true;
        return false;
    };

    for (int i = 1; i < argc; ++i)
    {
        std::string a = argv[i];
        if (is_flag(a, {"-h", "--help", "help", "-H"}))
        {
            std::cout
                << "Usage: lsdtt-flat-surface-finder <path> <parameter_file>\n"
                << "       lsdtt-flat-surface-finder <parameter_file>\n"
                << "       lsdtt-flat-surface-finder -h   (this help)\n"
                << "       lsdtt-flat-surface-finder -v   (version)\n";
            std::exit(0);
        }
        if (is_flag(a, {"-v", "--version", "version"}))
        {
            std::cout << "lsdtt-flat-surface-finder version " << k_version << "\n";
            std::exit(0);
        }
    }

    std::string path = "./";
    std::string fname;

    if (argc == 1)
    {
        std::cerr << "No parameter file specified. Use -h for help.\n";
        std::exit(1);
    }
    else if (argc == 2) { fname = argv[1]; }
    else                { path = argv[1]; fname = argv[2]; }

    return {path, fname};
}

// Helper: write a QQ dataset to CSV
static void write_qq_csv(
    lsdtt::lsdtt_data_loader&              loader,
    const std::string&                     out_path,
    const std::vector<lsdtt::qq_bin>&      qq)
{
    std::vector<std::vector<std::string>> rows;
    rows.reserve(qq.size());
    for (const auto& b : qq)
    {
        rows.push_back({
            std::to_string(b.standard_normal_quantile),
            std::to_string(b.data_quantile)
        });
    }
    loader.write_csv(out_path, {"standard_normal_quantile", "data_quantile"}, rows);
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
    {
        std::cerr << "Error loading parameter file: " << e.what() << "\n";
        return 1;
    }

    // -------------------------------------------------------------------------
    // Parameter declarations
    // -------------------------------------------------------------------------
    std::map<std::string, bool>                     bool_default_map;
    std::map<std::string, float>                    float_default_map;
    std::map<std::string, int>                      int_default_map;
    std::map<std::string, std::string>              string_default_map;
    std::map<std::string, std::vector<std::string>> help_map;

    // Threading / format
    int_default_map["n_threads"] = 0;
    help_map["n_threads"] = {"int","0","OpenMP thread count.","0 = all hardware threads."};

    string_default_map["write_output_format"] = "cog";
    help_map["write_output_format"] = {"string","cog","Raster output format.",
        "cog, gtiff, or same_as_input."};

    // -------------------------------------------------------------------------
    // MADk2
    // -------------------------------------------------------------------------
    bool_default_map["use_madk2"] = true;
    help_map["use_madk2"] = {"bool","true",
        "Enable MADk2 second-order roughness metric.",
        "Recommended for all landscape types, including sloping flat surfaces "
        "(alluvial fans, old terraces).  Second-order differences remove the "
        "linear slope trend so regional slope does not inflate roughness."};

    float_default_map["madk2_window_radius_metres"] = 100.0f;
    help_map["madk2_window_radius_metres"] = {"float","100",
        "Moving-window radius (m) for MADk2.",
        "Choose a radius somewhat smaller than the width of the flat surfaces "
        "to be detected so the window stays within a single flat unit."};

    float_default_map["madk2_threshold_percentile"] = 25.0f;
    help_map["madk2_threshold_percentile"] = {"float","25",
        "Percentile of the MADk2 distribution below which a pixel is classed flat.",
        "Inspect the MADk2 QQ CSV to identify the break-point in the distribution "
        "and adjust this value accordingly.  Range: 0–100."};

    // -------------------------------------------------------------------------
    // Moving-window IQR
    // -------------------------------------------------------------------------
    bool_default_map["use_window_iqr"] = true;
    help_map["use_window_iqr"] = {"bool","true",
        "Enable moving-window IQR metric.",
        "IQR of raw elevations within the window.  Effective for nearly level "
        "flats (floodplains).  For sloping surfaces, use a small window radius "
        "relative to the expected flat width, or prefer MADk2."};

    float_default_map["iqr_window_radius_metres"] = 50.0f;
    help_map["iqr_window_radius_metres"] = {"float","50",
        "Moving-window radius (m) for IQR.",
        "Smaller than the MADk2 radius by default because IQR is more sensitive "
        "to regional slope at large windows."};

    string_default_map["iqr_window_shape"] = "circular";
    help_map["iqr_window_shape"] = {"string","circular",
        "Window shape for IQR computation.",
        "circular or square."};

    float_default_map["iqr_threshold_percentile"] = 25.0f;
    help_map["iqr_threshold_percentile"] = {"float","25",
        "Percentile of the IQR distribution below which a pixel is classed flat.",
        "Range: 0–100."};

    // -------------------------------------------------------------------------
    // NVD
    // -------------------------------------------------------------------------
    bool_default_map["use_nvd"] = false;
    help_map["use_nvd"] = {"bool","false",
        "Enable NVD (Normal Vector Dispersion s3) metric.",
        "Requires FFTW.  Low s3 = co-planar normals = flat surface.  "
        "Slower than MADk2 and IQR due to the wavelet quadratic fit."};

    float_default_map["nvd_polyfit_lambda"] = 10.0f;
    help_map["nvd_polyfit_lambda"] = {"float","10",
        "Wavelet wavelength (m) for the quadratic fit used to derive d/e coefficients.",
        "Should be roughly equal to the cell size × a few cells to capture "
        "local slope direction while suppressing noise.  Only used when use_nvd is true."};

    float_default_map["nvd_window_radius_metres"] = 100.0f;
    help_map["nvd_window_radius_metres"] = {"float","100",
        "Moving-window radius (m) for NVD orientation tensor accumulation.",
        "Only used when use_nvd is true."};

    string_default_map["nvd_window_shape"] = "circular";
    help_map["nvd_window_shape"] = {"string","circular",
        "Window shape for NVD computation.",
        "circular or square."};

    float_default_map["nvd_threshold_percentile"] = 25.0f;
    help_map["nvd_threshold_percentile"] = {"float","25",
        "Percentile of the NVD s3 distribution below which a pixel is classed flat.",
        "Range: 0–100."};

    // -------------------------------------------------------------------------
    // Combination and component filtering
    // -------------------------------------------------------------------------
    bool_default_map["require_all_metrics_pass"] = true;
    help_map["require_all_metrics_pass"] = {"bool","true",
        "Require a pixel to pass ALL enabled metric thresholds to be classed flat.",
        "true (AND, default): conservative; use when false positive rate must be low. "
        "false (OR): permissive; use when detecting subtle or degraded flat surfaces."};

    int_default_map["min_component_pixels"] = 100;
    help_map["min_component_pixels"] = {"int","100",
        "Minimum 8-connected component size for inclusion in the labelled flat mask.",
        "Components smaller than this are excluded.  Set to 1 to include all components."};

    // -------------------------------------------------------------------------
    // QQ
    // -------------------------------------------------------------------------
    int_default_map["qq_n_bins"] = 100;
    help_map["qq_n_bins"] = {"int","100",
        "Number of quantile bins in the QQ output CSVs.",
        "100 gives one row per percentile point."};

    // -------------------------------------------------------------------------
    // Output toggles
    // -------------------------------------------------------------------------
    bool_default_map["write_texture_rasters"] = false;
    help_map["write_texture_rasters"] = {"bool","false",
        "Write per-metric texture rasters (MADk2, IQR, NVD s3).",
        "Useful for threshold selection and visual inspection."};

    bool_default_map["write_metric_flat_masks"] = false;
    help_map["write_metric_flat_masks"] = {"bool","false",
        "Write binary flat masks for each individual enabled metric.",
        "Useful for diagnosing which metric contributes most to the combined mask."};

    bool_default_map["write_combined_flat_mask"] = true;
    help_map["write_combined_flat_mask"] = {"bool","true",
        "Write the combined binary flat mask.",
        "1 = flat by the combined criterion, 0 = not flat."};

    // (labelled mask and QQ CSVs are always written)

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
    // Cross-parameter validation and resolution
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
        if (s == "square") return lsdtt::window_shape::square;
        return lsdtt::window_shape::circular;
    };

    if (!bf["use_madk2"] && !bf["use_window_iqr"] && !bf["use_nvd"])
    {
        std::cerr << "Error: at least one metric must be enabled "
                     "(use_madk2, use_window_iqr, or use_nvd).\n";
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
        {
            std::cerr << "Error loading DEM: " << e.what() << "\n";
            return 1;
        }
        dem = lsdtt::lsdtt_raster_cpu(loader);
        try { dem.get_meta().check_projected_metres(); }
        catch (const std::exception& e)
        {
            std::cerr << "Projection error: " << e.what() << "\n";
            return 1;
        }
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
    req.use_window_iqr             = bf["use_window_iqr"];
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
    req.n_threads                  = ip["n_threads"];

    // -------------------------------------------------------------------------
    // Run
    // -------------------------------------------------------------------------
    std::cout << "Running flat surface finder...\n";
    if (req.use_madk2)
        std::cout << "  MADk2    radius=" << req.madk2_window_radius_metres
                  << " m  threshold=" << req.madk2_threshold_percentile << "th pct\n";
    if (req.use_window_iqr)
        std::cout << "  IQR      radius=" << req.iqr_window_radius_metres
                  << " m  threshold=" << req.iqr_threshold_percentile << "th pct\n";
    if (req.use_nvd)
        std::cout << "  NVD s3   radius=" << req.nvd_window_radius_metres
                  << " m  threshold=" << req.nvd_threshold_percentile << "th pct\n";
    std::cout << "  combination: " << (req.require_all_metrics_pass ? "AND" : "OR") << "\n";
    std::cout << "  min component pixels: " << req.min_component_pixels << "\n\n";

    lsdtt::flat_surface_result result;
    try { result = lsdtt::find_flat_surfaces(dem, req); }
    catch (const std::exception& e)
    {
        std::cerr << "Error running find_flat_surfaces: " << e.what() << "\n";
        return 1;
    }

    // -------------------------------------------------------------------------
    // Write outputs
    // -------------------------------------------------------------------------
    const std::string stem = pp.get_full_write_path();

    // Texture rasters (optional)
    if (bf["write_texture_rasters"])
    {
        if (req.use_madk2)
        {
            const std::string p = stem + "_madk2_roughness.tif";
            std::cout << "Writing MADk2 roughness: " << p << "\n";
            loader.write_full(result.madk2_roughness.get_data(),
                              result.madk2_roughness.get_meta(), p, out_format);
        }
        if (req.use_window_iqr)
        {
            const std::string p = stem + "_window_iqr.tif";
            std::cout << "Writing window IQR: " << p << "\n";
            loader.write_full(result.window_iqr.get_data(),
                              result.window_iqr.get_meta(), p, out_format);
        }
        if (req.use_nvd)
        {
            const std::string p = stem + "_nvd_s3.tif";
            std::cout << "Writing NVD s3: " << p << "\n";
            loader.write_full(result.nvd_s3.get_data(),
                              result.nvd_s3.get_meta(), p, out_format);
        }
    }

    // Per-metric flat masks (optional)
    if (bf["write_metric_flat_masks"])
    {
        if (req.use_madk2)
        {
            const std::string p = stem + "_madk2_flat_mask.tif";
            std::cout << "Writing MADk2 flat mask: " << p << "\n";
            loader.write_full(result.madk2_flat_mask.get_data(),
                              result.madk2_flat_mask.get_meta(), p, out_format);
        }
        if (req.use_window_iqr)
        {
            const std::string p = stem + "_iqr_flat_mask.tif";
            std::cout << "Writing IQR flat mask: " << p << "\n";
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
    }

    // Combined flat mask
    if (bf["write_combined_flat_mask"])
    {
        const std::string p = stem + "_combined_flat_mask.tif";
        std::cout << "Writing combined flat mask: " << p << "\n";
        loader.write_full(result.combined_flat_mask.get_data(),
                          result.combined_flat_mask.get_meta(), p, out_format);
    }

    // Labelled flat mask (always written — it is the primary result)
    {
        const std::string p = stem + "_labelled_flat_mask.tif";
        std::cout << "Writing labelled flat mask (UInt32): " << p << "\n";
        loader.write_full(result.labelled_flat_mask.get_data(),
                          result.labelled_flat_mask.get_meta(), p, out_format,
                          GDT_UInt32);
    }

    // QQ CSVs (always written for each enabled metric)
    if (req.use_madk2 && !result.madk2_qq.empty())
    {
        const std::string p = stem + "_madk2_qq.csv";
        std::cout << "Writing MADk2 QQ: " << p << "\n";
        write_qq_csv(loader, p, result.madk2_qq);
    }
    if (req.use_window_iqr && !result.iqr_qq.empty())
    {
        const std::string p = stem + "_iqr_qq.csv";
        std::cout << "Writing IQR QQ: " << p << "\n";
        write_qq_csv(loader, p, result.iqr_qq);
    }
    if (req.use_nvd && !result.nvd_qq.empty())
    {
        const std::string p = stem + "_nvd_qq.csv";
        std::cout << "Writing NVD QQ: " << p << "\n";
        write_qq_csv(loader, p, result.nvd_qq);
    }

    std::cout << "\n";
    pp.print_parameters();
    std::cout << "\nlsdtt-flat-surface-finder complete.\n";
    return 0;
}
