// lsdtt-mountain-front.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html
//
// -----------------------------------------------------------------------------
// lsdtt-mountain-front
//
// Identifies mountain masses in a DEM via adaptive relief thresholding.
//
// Pipeline:
//   Stage 1 — Local relief (square window, relief_window_pixels)
//   Stage 2 — Adaptive threshold: 500-bin empirical CDF vs best-fit parametric
//             distribution; the x-value of maximum positive divergence (where
//             flat pixels inflate the CDF above the fitted model) defines the
//             relief threshold.  Prints chosen percentile and relief value.
//             Writes a CDF CSV.
//   Stage 3 — Flat mask (relief <= threshold).
//   Stage 4 — Remove small flat components (< min_flat_component_pixels).
//   Stage 5 — Dilate flat mask by 1 pixel (8-connected square SE) to bridge
//             riser pixels between adjacent flat patches.
//   Stage 6 — Fill holes in flat mask: enclosed non-flat pixels become flat.
//   Stage 7 — Invert flat mask → initial mountain mask.
//   Stage 8 — Opening by reconstruction (erosion_radius_metres): erode to sever
//             thin necks, then BFS-reconstruct thick bodies back to original extent.
//   Stage 9 — Remove small mountain components (< min_component_pixels).
//   Stage 10 — Fill holes in mountain mask (enclosed background → mountain).
//   Stage 11 — Valley infilling: morphological closing (infill_valley_radius_metres).
//              Disabled when infill_valley_radius_metres <= 0.
//   Stage 12 — Outlets: mountain edge pixels with upstream area >= min_contributing_pixels.
//
// Nodata semantics:
//   Pixels that are nodata in the original DEM are nodata in all output rasters.
//
// Outputs:
//   <stem>_relief_cdf.csv            — empirical + fitted CDF (always written)
//   <stem>_mountain_mask.tif         — binary mountain raster
//   <stem>_relief.tif                — local relief raster (optional)
//   <stem>_flat_mask.tif             — binary flat mask after processing (optional)
//   <stem>_mountain_labels.tif       — labelled mountain components (optional)
//   <stem>_mountain_front_outlets.*  — outlet points (optional)
//
// Usage:
//   lsdtt-mountain-front <path> <parameter_file>
//   lsdtt-mountain-front <parameter_file>   (path assumed to be ./)
//   lsdtt-mountain-front -h
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
#include "../src/core/lsdtt_flow_base.hpp"
#include "../src/core/lsdtt_param_parser.hpp"
#include "../src/core/lsdtt_stat_csv.hpp"
#include "../src/core/lsdtt_raster_meta.hpp"
#include "../src/core/lsdtt_vector_utils.hpp"
#include "../src/cpu/lsdtt_mountain_front_cpu.hpp"
#include "../src/cpu/lsdtt_raster_cpu.hpp"
#include "../src/cpu/lsdtt_raster_feature_finder_cpu.hpp"
#include "../src/cpu/lsdtt_raster_flow_preprocessing_cpu.hpp"
#include "../src/cpu/lsdtt_raster_surface_metrics_cpu.hpp"
#include "../src/utils/lsdtt_distributions.hpp"
#include "../src/utils/lsdtt_histogram.hpp"

static const std::string k_version  = "0.2.0";
static const std::string k_citation = "https://github.com/LSDtopotools/lsdtt3";

static void print_banner()
{
    std::cout
        << "=========================================================\n"
        << "|| Welcome to lsdtt-mountain-front                      ||\n"
        << "|| I will identify mountain masses in a DEM for you     ||\n"
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
        if (is_flag(a, {"-h","--help","help","-H","--H"}))
        {
            std::cout
                << "Usage: lsdtt-mountain-front <path> <parameter_file>\n"
                << "       lsdtt-mountain-front <parameter_file>\n"
                << "       lsdtt-mountain-front -h\n";
            std::exit(0);
        }
        if (is_flag(a, {"-v","--version","version","-V"}))
        { std::cout << "lsdtt-mountain-front version " << k_version << "\n"; std::exit(0); }
        if (is_flag(a, {"-c","--cite","citation","-citation","--citation"}))
        { std::cout << "Citation: " << k_citation << "\n"; std::exit(0); }
    }
    std::string path = "./";
    std::string fname;
    if (argc == 1) { std::cerr << "No parameter file. Use -h for help.\n"; std::exit(1); }
    else if (argc == 2) { fname = argv[1]; }
    else { path = argv[1]; fname = argv[2]; }
    return {path, fname};
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

    // --- Stage 1: relief ---
    int_default_map["relief_window_pixels"] = 5;
    help_map["relief_window_pixels"] = {"int","5",
        "Square window size (pixels) for local relief (max - min).",
        "Must be odd and >= 3."};

    // --- Stages 4-6: flat mask post-processing ---
    int_default_map["min_flat_component_pixels"] = 500;
    help_map["min_flat_component_pixels"] = {"int","500",
        "Minimum 8-connected component size (pixels) to keep in the flat mask.",
        "Small isolated flat patches (e.g. individual floodplain pixels) are "
        "removed before the flat mask is used to define mountain pixels."};

    // --- Stages 8-9: mountain mask post-processing ---
    float_default_map["erosion_radius_metres"] = 500.0f;
    help_map["erosion_radius_metres"] = {"float","500",
        "Erosion radius (m) for opening-by-reconstruction (Stage 8).",
        "Severs necks narrower than 2x this radius; thick mountain bodies are "
        "restored to their exact original extent.  Set to roughly half the neck "
        "width to remove.  0 = skip."};

    int_default_map["min_component_pixels"] = 1000;
    help_map["min_component_pixels"] = {"int","1000",
        "Minimum 8-connected component size (pixels) to keep in the mountain mask.",
        "Applied to the eroded mountain mask after neck severing."};

    // --- Stage 11: valley infilling ---
    float_default_map["infill_valley_radius_metres"] = 0.0f;
    help_map["infill_valley_radius_metres"] = {"float","0",
        "Radius (m) of the circular SE for morphological closing (valley infilling).",
        "Seals valley notches cut by incised rivers into the mountain front. "
        "Set to ~half the valley width to be filled.  0 = disabled."};

    // --- Stage 12: outlets ---
    string_default_map["contributing_pixels_fname"] = "";
    help_map["contributing_pixels_fname"] = {"string","",
        "Filename of a pre-computed n_contributing_pixels raster.",
        "Leave empty to auto-compute via fill + D8 flow routing."};

    int_default_map["min_contributing_pixels"] = 1000;
    help_map["min_contributing_pixels"] = {"int","1000",
        "Minimum upstream contributing pixels for a mountain edge pixel to be an outlet.",
        "Only used when write_mountain_front_outlets is true."};

    // --- Outputs ---
    bool_default_map["write_relief_raster"] = false;
    help_map["write_relief_raster"] = {"bool","false",
        "Write the local relief raster.","Useful for threshold inspection."};

    bool_default_map["write_flat_mask_raster"] = false;
    help_map["write_flat_mask_raster"] = {"bool","false",
        "Write the processed flat mask (after component filter, dilation, hole fill).",
        "The mountain mask is the inverse of this."};

    bool_default_map["write_mountain_mask_raster"] = true;
    help_map["write_mountain_mask_raster"] = {"bool","true",
        "Write the final mountain mask raster.","1 = mountain, 0 = flat / piedmont."};

    bool_default_map["write_mountain_labels_raster"] = false;
    help_map["write_mountain_labels_raster"] = {"bool","false",
        "Write a labelled mountain components raster (UInt32).",
        "Each 8-connected mountain component receives a unique integer label."};

    bool_default_map["write_mountain_front_outlets"] = false;
    help_map["write_mountain_front_outlets"] = {"bool","false",
        "Write mountain front outlet points (FlatGeobuf + CSV).",
        "Mountain mask edge pixels whose upstream area meets min_contributing_pixels."};

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
    // Cross-parameter logic
    // -------------------------------------------------------------------------
    lsdtt::raster_output_format out_format = lsdtt::raster_output_format::cog;
    {
        const std::string& fmt = sf["write_output_format"];
        if      (fmt == "gtiff" || fmt == "geotiff") out_format = lsdtt::raster_output_format::gtiff;
        else if (fmt == "same_as_input")             out_format = lsdtt::raster_output_format::same_as_input;
    }
    const int n_threads = ip["n_threads"];

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

    const std::string stem = pp.get_full_write_path();

    // -------------------------------------------------------------------------
    // Stage 1: local relief (square window)
    // -------------------------------------------------------------------------
    const int wpix = ip["relief_window_pixels"];
    if (wpix < 3 || wpix % 2 == 0)
    { std::cerr << "Error: relief_window_pixels must be odd and >= 3\n"; return 1; }

    std::cout << "Stage 1: computing local relief (square window, "
              << wpix << " pixels)...\n";

    const double relief_radius =
        static_cast<double>(std::max(1, (wpix - 1) / 2)) * dem.get_meta().cell_size;
    const lsdtt::lsdtt_raster_cpu relief =
        lsdtt::moving_window_local_relief(dem, relief_radius, lsdtt::window_shape::square);

    if (bf["write_relief_raster"])
    {
        const std::string p = stem + "_relief.tif";
        std::cout << "  Writing relief raster: " << p << "\n";
        loader.write_full(relief.get_data(), relief.get_meta(), p, out_format);
    }

    // -------------------------------------------------------------------------
    // Stage 2: adaptive relief threshold via 500-bin CDF divergence
    // -------------------------------------------------------------------------
    std::cout << "Stage 2: computing adaptive relief threshold...\n";

    const lsdtt::lsdtt_raster_meta& rmeta = relief.get_meta();
    const std::vector<float>&        rdata = relief.get_data();
    const std::vector<uint8_t>&      rmask = relief.get_mask();

    std::vector<double> relief_vals;
    relief_vals.reserve(rdata.size() / 4);
    for (std::size_t i = 0; i < rdata.size(); ++i)
    {
        if (rmask[i] == 0) continue;
        const float v = rdata[i];
        if (!std::isfinite(v) || rmeta.is_nodata(v)) continue;
        relief_vals.push_back(static_cast<double>(v));
    }

    if (relief_vals.size() < 10)
    { std::cerr << "Error: too few valid relief values to compute adaptive threshold.\n"; return 1; }

    const lsdtt::lsdtt_smoothed_cdf               scdf = lsdtt::compute_smoothed_cdf(relief_vals, 500);
    const std::vector<lsdtt::lsdtt_distribution_fit> fits = lsdtt::fit_distributions(relief_vals);

    const lsdtt::lsdtt_distribution_fit* best_fit = nullptr;
    for (const auto& f : fits)
    {
        if (!f.valid || !std::isfinite(f.ks_statistic)) continue;
        if (!best_fit || f.ks_statistic < best_fit->ks_statistic) best_fit = &f;
    }
    if (!best_fit)
    { std::cerr << "Error: no valid distribution fit for relief CDF.\n"; return 1; }

    double relief_threshold  = std::numeric_limits<double>::quiet_NaN();
    double relief_percentile = std::numeric_limits<double>::quiet_NaN();
    {
        double      max_diff = -std::numeric_limits<double>::infinity();
        std::size_t max_idx  = 0;
        for (std::size_t i = 0; i < scdf.bin_centres.size(); ++i)
        {
            const double fitted = lsdtt::eval_cdf(*best_fit, scdf.bin_centres[i]);
            if (!std::isfinite(fitted)) continue;
            const double diff = scdf.cdf[i] - fitted;
            if (diff > max_diff) { max_diff = diff; max_idx = i; }
        }
        relief_threshold  = scdf.bin_centres[max_idx];
        relief_percentile = scdf.cdf[max_idx] * 100.0;
    }

    std::cout << "  Distribution fit  : "
              << lsdtt::distribution_family_name(best_fit->family) << "\n"
              << "  Relief threshold  : " << relief_threshold << " m"
              << "  (=" << relief_percentile << "th percentile)\n\n";

    {
        const std::string p = stem + "_relief_cdf.csv";
        std::cout << "Writing relief CDF: " << p << "\n";
        write_smoothed_cdf_csv(loader, p, scdf, fits);
    }

    // -------------------------------------------------------------------------
    // Stage 3: build initial flat mask
    //
    // Nodata semantics: pixels nodata in the DEM are nodata in all masks.
    // We use the DEM mask as the authoritative nodata source.
    // -------------------------------------------------------------------------
    std::cout << "Stage 3: building flat mask (relief <= "
              << relief_threshold << " m)...\n";

    const std::vector<uint8_t>& dem_mask = dem.get_mask();
    const std::size_t npix =
        static_cast<std::size_t>(rmeta.n_rows) * static_cast<std::size_t>(rmeta.n_cols);

    const float nd_val = static_cast<float>(rmeta.no_data_value);
    std::vector<float>   flat_data(npix, nd_val);
    std::vector<uint8_t> flat_msk(npix, uint8_t{0});

    for (std::size_t i = 0; i < npix; ++i)
    {
        if (dem_mask[i] == 0) continue;          // DEM nodata → nodata in output
        const float v = rdata[i];
        if (rmask[i] == 0 || !std::isfinite(v) || rmeta.is_nodata(v)) continue;
        flat_msk[i]  = uint8_t{1};
        flat_data[i] = (static_cast<double>(v) <= relief_threshold) ? 1.0f : 0.0f;
    }

    lsdtt::lsdtt_raster_cpu flat_mask(
        std::move(flat_data), std::move(flat_msk), rmeta);

    // -------------------------------------------------------------------------
    // Stage 4: remove small flat components
    // -------------------------------------------------------------------------
    const int min_flat = ip["min_flat_component_pixels"];
    std::cout << "Stage 4: removing flat components < " << min_flat << " pixels...\n";
    flat_mask = lsdtt::filter_binary_components_by_size(flat_mask, min_flat);

    // -------------------------------------------------------------------------
    // Stage 5: dilate flat mask by 1 pixel (square SE) to bridge riser gaps
    // -------------------------------------------------------------------------
    std::cout << "Stage 5: dilating flat mask by 1 pixel to bridge risers...\n";
    flat_mask = lsdtt::binary_morphological_dilation(
        flat_mask,
        dem.get_meta().cell_size,       // radius = 1 cell → 3×3 square SE
        lsdtt::window_shape::square);

    // -------------------------------------------------------------------------
    // Stage 6: fill holes in flat mask
    // -------------------------------------------------------------------------
    std::cout << "Stage 6: filling holes in flat mask...\n";
    flat_mask = lsdtt::binary_fill_holes(flat_mask);

    // Enforce DEM nodata: any pixel that is nodata in the DEM must remain
    // nodata in the flat mask regardless of what intermediate stages produced.
    {
        std::vector<float>   fd = flat_mask.get_data();
        std::vector<uint8_t> fm = flat_mask.get_mask();
        for (std::size_t i = 0; i < npix; ++i)
            if (dem_mask[i] == 0) { fm[i] = uint8_t{0}; fd[i] = nd_val; }
        flat_mask = lsdtt::lsdtt_raster_cpu(std::move(fd), std::move(fm), rmeta);
    }

    if (bf["write_flat_mask_raster"])
    {
        const std::string p = stem + "_flat_mask.tif";
        std::cout << "  Writing flat mask: " << p << "\n";
        loader.write_full(flat_mask.get_data(), flat_mask.get_meta(), p, out_format);
    }

    // -------------------------------------------------------------------------
    // Stage 7: invert flat mask → mountain mask
    // -------------------------------------------------------------------------
    std::cout << "Stage 7: inverting flat mask to mountain mask...\n";
    {
        const std::vector<float>&   fd = flat_mask.get_data();
        const std::vector<uint8_t>& fm = flat_mask.get_mask();

        std::vector<float>   mtn_data(npix, 0.0f);
        std::vector<uint8_t> mtn_msk(npix, uint8_t{0});
        for (std::size_t i = 0; i < npix; ++i)
        {
            if (fm[i] == 0) continue;
            mtn_msk[i]  = uint8_t{1};
            mtn_data[i] = (fd[i] > 0.5f) ? 0.0f : 1.0f;
        }
        flat_mask = lsdtt::lsdtt_raster_cpu(
            std::move(mtn_data), std::move(mtn_msk), rmeta);
    }
    lsdtt::lsdtt_raster_cpu mountain_mask = std::move(flat_mask);

    // -------------------------------------------------------------------------
    // Stage 8: opening by reconstruction — sever thin necks, restore thick bodies
    // -------------------------------------------------------------------------
    const double erosion_r = static_cast<double>(ff["erosion_radius_metres"]);
    if (erosion_r > 0.0)
    {
        std::cout << "Stage 8: opening by reconstruction (radius="
                  << erosion_r << " m)...\n";
        mountain_mask = lsdtt::binary_opening_by_reconstruction(mountain_mask, erosion_r);
    }
    else
    {
        std::cout << "Stage 8: opening by reconstruction skipped (erosion_radius_metres <= 0).\n";
    }

    // -------------------------------------------------------------------------
    // Stage 9: remove small mountain components
    // -------------------------------------------------------------------------
    const int min_comp = ip["min_component_pixels"];
    std::cout << "Stage 9: removing mountain components < " << min_comp << " pixels...\n";
    mountain_mask = lsdtt::filter_binary_components_by_size(mountain_mask, min_comp);

    // -------------------------------------------------------------------------
    // Stage 10: fill holes in mountain mask
    // -------------------------------------------------------------------------
    std::cout << "Stage 10: filling holes in mountain mask...\n";
    mountain_mask = lsdtt::binary_fill_holes(mountain_mask);

    // -------------------------------------------------------------------------
    // Stage 11: valley infilling (morphological closing)
    // -------------------------------------------------------------------------
    const double infill_r = static_cast<double>(ff["infill_valley_radius_metres"]);
    if (infill_r > 0.0)
    {
        std::cout << "Stage 11: valley infilling (closing radius="
                  << infill_r << " m)...\n";
        mountain_mask = lsdtt::binary_morphological_closing(mountain_mask, infill_r);
    }
    else
    {
        std::cout << "Stage 11: valley infilling skipped.\n";
    }

    // Enforce DEM nodata on final mountain mask.
    {
        std::vector<float>   md = mountain_mask.get_data();
        std::vector<uint8_t> mm = mountain_mask.get_mask();
        for (std::size_t i = 0; i < npix; ++i)
            if (dem_mask[i] == 0) { mm[i] = uint8_t{0}; md[i] = nd_val; }
        mountain_mask = lsdtt::lsdtt_raster_cpu(std::move(md), std::move(mm), rmeta);
    }

    // -------------------------------------------------------------------------
    // Write mountain mask
    // -------------------------------------------------------------------------
    if (bf["write_mountain_mask_raster"])
    {
        const std::string p = stem + "_mountain_mask.tif";
        std::cout << "Writing mountain mask: " << p << "\n";
        loader.write_full(mountain_mask.get_data(),
                          mountain_mask.get_meta(), p, out_format);
    }

    if (bf["write_mountain_labels_raster"])
    {
        const std::string p = stem + "_mountain_labels.tif";
        std::cout << "Writing mountain component labels: " << p << "\n";
        const lsdtt::lsdtt_raster_cpu labels =
            lsdtt::label_binary_components_by_size(mountain_mask, std::max(1, min_comp));
        loader.write_full(labels.get_data(), labels.get_meta(), p, out_format, GDT_UInt32);
    }

    // -------------------------------------------------------------------------
    // Stage 12: outlets
    // -------------------------------------------------------------------------
    if (bf["write_mountain_front_outlets"])
    {
        std::cout << "\nStage 12: computing mountain front outlets...\n";

        lsdtt::lsdtt_raster_cpu ca_raster;
        bool ca_ok = false;

        const std::string ca_fname = sf["contributing_pixels_fname"];
        if (ca_fname.empty())
        {
            std::cout << "  contributing_pixels_fname not set — computing from DEM.\n";
            try
            {
                std::cout << "  Filling pits...\n";
                const lsdtt::lsdtt_raster_cpu filled =
                    lsdtt::fill_pits(dem, 0.0001f, n_threads);

                std::cout << "  Building D8 flow routing...\n";
                lsdtt::lsdtt_flowinfo_options fi_opts;
                fi_opts.interior_nodata_is_outlet = true;
                fi_opts.require_acyclic_graph     = true;
                fi_opts.build_topological_order   = true;
                fi_opts.build_s_vector            = true;
                fi_opts.n_threads                 = n_threads;

                const std::unique_ptr<lsdtt::lsdtt_flow_base> flow =
                    lsdtt::make_flow_base(
                        filled.get_data_array(), filled.get_mask_array(),
                        filled.get_meta(), fi_opts);

                std::cout << "  Extracting contributing pixels raster...\n";
                ca_raster = lsdtt::lsdtt_raster_cpu(
                    flow->get_n_contributing_pixels_raster());
                ca_ok = true;
                std::cout << "  Done.\n";
            }
            catch (const std::exception& e)
            {
                std::cerr << "Warning: auto-compute of n_contributing_pixels failed: "
                          << e.what() << " — skipping outlets.\n";
            }
        }
        else
        {
            const std::string ca_path = pp.get_read_path() + "/" + ca_fname;
            std::cout << "  Loading contributing pixels raster: " << ca_path << "\n";
            try
            {
                lsdtt::lsdtt_data_loader ca_loader;
                ca_loader.open_raster(ca_path);
                ca_raster = lsdtt::lsdtt_raster_cpu(ca_loader);
                ca_ok = true;
            }
            catch (const std::exception& e)
            {
                std::cerr << "Warning: failed to load contributing pixels raster: "
                          << e.what() << " — skipping outlets.\n";
            }
        }

        if (ca_ok)
        {
            try
            {
                const lsdtt::lsdtt_raster_cpu outlets =
                    lsdtt::find_mountain_front_outlets(
                        mountain_mask, ca_raster, ip["min_contributing_pixels"]);

                const lsdtt::lsdtt_raster_meta& mta  = outlets.get_meta();
                const std::vector<float>&        odat = outlets.get_data();
                const std::vector<uint8_t>&      omsk = outlets.get_mask();
                const std::vector<float>&        cdat = ca_raster.get_data();
                const std::vector<uint8_t>&      cmsk = ca_raster.get_mask();
                const int onrows = mta.n_rows;
                const int oncols = mta.n_cols;

                // FlatGeobuf
                {
                    const std::string fgb_path = stem + "_mountain_front_outlets.fgb";
                    std::cout << "Writing mountain front outlets (FlatGeobuf): "
                              << fgb_path << "\n";
                    lsdtt::lsdtt_vector_schema schema;
                    schema.field_names = {"row","col","contributing_pixels"};
                    schema.field_types  = {
                        lsdtt::lsdtt_vector_field_type::int64_type,
                        lsdtt::lsdtt_vector_field_type::int64_type,
                        lsdtt::lsdtt_vector_field_type::int64_type };
                    loader.begin_vector_write(
                        fgb_path, "mountain_front_outlets",
                        schema, mta.crs_wkt,
                        lsdtt::vector_output_format::flatgeobuf, false);
                    for (int row = 0; row < onrows; ++row)
                    {
                        for (int col = 0; col < oncols; ++col)
                        {
                            const std::size_t idx =
                                static_cast<std::size_t>(row) *
                                static_cast<std::size_t>(oncols) + col;
                            if (omsk[idx] == 0 || odat[idx] <= 0.5f) continue;
                            const lsdtt::lsdtt_xy xy =
                                lsdtt::cell_center_from_row_col(mta, row, col);
                            const std::int64_t ca_px =
                                (cmsk[idx] != 0 && std::isfinite(cdat[idx]))
                                    ? static_cast<std::int64_t>(cdat[idx]) : 0LL;
                            loader.write_point(xy.x, xy.y, {
                                static_cast<std::int64_t>(row),
                                static_cast<std::int64_t>(col), ca_px });
                        }
                    }
                    loader.finish_vector_write();
                }

                // CSV (lat/lon)
                {
                    const std::string csv_path = stem + "_mountain_front_outlets.csv";
                    std::cout << "Writing mountain front outlets (CSV): "
                              << csv_path << "\n";
                    lsdtt::lsdtt_coordinate_transform to_wgs84 =
                        lsdtt::lsdtt_coordinate_transform::from_epsg(mta.crs_wkt, 4326);
                    std::vector<std::vector<std::string>> rows;
                    for (int row = 0; row < onrows; ++row)
                    {
                        for (int col = 0; col < oncols; ++col)
                        {
                            const std::size_t idx =
                                static_cast<std::size_t>(row) *
                                static_cast<std::size_t>(oncols) + col;
                            if (omsk[idx] == 0 || odat[idx] <= 0.5f) continue;
                            const lsdtt::lsdtt_lonlat ll =
                                lsdtt::cell_center_to_lonlat(mta, row, col, to_wgs84);
                            rows.push_back({
                                std::to_string(ll.latitude),
                                std::to_string(ll.longitude) });
                        }
                    }
                    loader.write_csv(csv_path, {"latitude","longitude"}, rows);
                }
            }
            catch (const std::exception& e)
            {
                std::cerr << "Error computing mountain front outlets: "
                          << e.what() << "\n";
            }
        }
    }

    std::cout << "\n";
    pp.print_parameters();
    std::cout << "\nlsdtt-mountain-front complete.\n";
    return 0;
}
