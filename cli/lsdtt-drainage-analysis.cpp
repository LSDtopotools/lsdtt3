// lsdtt-drainage-analysis.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html
//
// -----------------------------------------------------------------------------
// lsdtt-drainage-analysis
//
// Compute drainage-related rasters from a DEM:
//
//   1. Optionally fill pits if the raster is not already filled
//   2. Compute drainage area using D8 (lsdtt_flowinfo), D-infinity, or MFD
//      (lsdtt_flowgraph_cpu)
//   3. Compute the Topographic Wetness Index (TWI)
//   4. Placeholders for connectivity and height_above_channel (coming soon)
//
// Usage:
//   lsdtt-drainage-analysis <path> <parameter_file>
//   lsdtt-drainage-analysis <parameter_file>   (path assumed to be ./)
//   lsdtt-drainage-analysis -h                 (print help and exit)
//   lsdtt-drainage-analysis -v                 (print version and exit)
//   lsdtt-drainage-analysis -c                 (print citation and exit)
// -----------------------------------------------------------------------------

#include <algorithm>
#include <cctype>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <functional>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "../src/core/lsdtt_data_loader.hpp"
#include "../src/core/lsdtt_flow_base.hpp"
#include "../src/core/lsdtt_raster_reproject.hpp"
#include "../src/core/lsdtt_types.hpp"
#include "../src/core/lsdtt_common_parameters.hpp"
#include "../src/core/lsdtt_param_parser.hpp"
#include "../src/cpu/lsdtt_connectivity_cpu.hpp"
#include "../src/cpu/lsdtt_flowgraph_cpu.hpp"
#include "../src/cpu/lsdtt_raster_algebra.hpp"
#include "../src/cpu/lsdtt_raster_cpu.hpp"
#include "../src/cpu/lsdtt_raster_flow_preprocessing_cpu.hpp"
#include "../src/cpu/lsdtt_raster_surface_metrics_cpu.hpp"
#include "../src/cpu/lsdtt_sti_cpu.hpp"
#include "../src/cpu/lsdtt_surface_texture_cpu.hpp"
#include "../src/utils/lsdtt_array2d.hpp"
#include "../src/utils/lsdtt_string_utils.hpp"

static const std::string k_version  = "0.1.0";
static const std::string k_citation = "https://github.com/LSDtopotools/lsdtt3";
static const std::string k_tool_name = "lsdtt-drainage-analysis";

// -----------------------------------------------------------------------------
// print_banner
// -----------------------------------------------------------------------------
static void print_banner()
{
    std::cout
        << "=========================================================\n"
        << "|| Welcome to lsdtt-drainage-analysis                  ||\n"
        << "|| I am going to compute drainage rasters for you      ||\n"
        << "|| Part of lsdtt3 by Simon M. Mudd                     ||\n"
        << "|| This CLI is for looking at various runoff metrics   ||\n"
        << "|| such as TWI, SWI and connectivity                   ||\n"
        << "=========================================================\n"
        << "|| CLI version: " << k_version  << "\n"
        << "|| Citation   : " << k_citation << "\n"
        << "|| -v or --version : full lsdtt3 release version\n"
        << "|| -h              : help explaining the parameters\n"
        << "=========================================================\n\n";
}

// -----------------------------------------------------------------------------
// keep_top_percentile_pixels
//
// Keep only valid pixels whose value is at or above the requested percentile
// threshold; all other valid pixels are set to nodata and masked out.
//
// Percentile uses type-7 linear interpolation (project decision D4),
// matching lsdtt::raster_percentile. It was nearest-rank before.
// -----------------------------------------------------------------------------
static lsdtt::lsdtt_raster_cpu keep_top_percentile_pixels(
    const lsdtt::lsdtt_raster_cpu& input,
    const double                   percentile,
    float&                         threshold,
    std::int64_t&                  retained_pixels)
{
    const auto& src_data = input.get_data();
    const auto& src_mask = input.get_mask();
    const auto& meta     = input.get_meta();

    if (src_data.size() != src_mask.size())
    {
        throw std::runtime_error(
            "keep_top_percentile_pixels: data/mask size mismatch");
    }

    std::vector<float> valid_values;
    valid_values.reserve(src_data.size());

    for (std::size_t i = 0; i < src_data.size(); ++i)
    {
        if (src_mask[i] == 0) continue;
        const float v = src_data[i];
        if (!std::isfinite(v) || meta.is_nodata(v)) continue;
        valid_values.push_back(v);
    }

    std::vector<float> out_data = src_data;
    std::vector<std::uint8_t> out_mask = src_mask;
    retained_pixels = 0;

    if (valid_values.empty())
    {
        threshold = std::numeric_limits<float>::quiet_NaN();
        return lsdtt::lsdtt_raster_cpu(std::move(out_data), std::move(out_mask), meta);
    }

    // TYPE-7 (linear interpolation), matching lsdtt::raster_percentile
    // (src/cpu/lsdtt_raster_algebra.cpp) and lsdtt::compute_percentile. This
    // replaced a NEAREST-RANK selection, ceil(p/100 * n), which returned an
    // observed data value; type-7 interpolates between the two values
    // straddling a fractional rank and is what R's quantile(), numpy and
    // MATLAB return by default. Project decision D4.
    //
    // WHY THE FORMULA IS REPEATED HERE RATHER THAN CALLING raster_percentile.
    // The library's collect_valid() admits a pixel on its MASK BYTE ALONE,
    // while the loop above additionally rejects non-finite values and the
    // nodata sentinel. On a raster obeying the usual invariant (mask 0
    // wherever nodata) the two agree, but the extra guards here are
    // deliberate, and a NaN admitted into the vector would corrupt the
    // nth_element ordering rather than fail loudly. Calling the library
    // directly would therefore change validity semantics as well as the
    // percentile definition, which is more than D4 decided. Resolve the mask
    // invariant first, then this whole function can go.
    const double p = std::clamp(percentile, 0.0, 100.0);
    const std::size_t n = valid_values.size();

    if (n == 1)
    {
        threshold = valid_values[0];
    }
    else
    {
        const double rank = (p / 100.0) * static_cast<double>(n - 1);
        const std::size_t lower_index = static_cast<std::size_t>(std::floor(rank));
        const std::size_t upper_index = static_cast<std::size_t>(std::ceil(rank));

        std::nth_element(valid_values.begin(),
                         valid_values.begin() + static_cast<std::ptrdiff_t>(upper_index),
                         valid_values.end());
        const float upper_value = valid_values[upper_index];

        if (lower_index == upper_index)
        {
            threshold = upper_value;
        }
        else
        {
            std::nth_element(valid_values.begin(),
                             valid_values.begin() + static_cast<std::ptrdiff_t>(lower_index),
                             valid_values.begin() + static_cast<std::ptrdiff_t>(upper_index) + 1);
            const float lower_value = valid_values[lower_index];
            const double frac = rank - std::floor(rank);
            threshold = static_cast<float>(
                static_cast<double>(lower_value) +
                frac * (static_cast<double>(upper_value) - static_cast<double>(lower_value)));
        }
    }

    const float nodata_f = static_cast<float>(meta.no_data_value);
    for (std::size_t i = 0; i < src_data.size(); ++i)
    {
        if (src_mask[i] == 0) continue;
        const float v = src_data[i];
        if (!std::isfinite(v) || meta.is_nodata(v)) continue;

        if (v >= threshold)
        {
            ++retained_pixels;
        }
        else
        {
            out_data[i] = nodata_f;
            out_mask[i] = 0;
        }
    }

    return lsdtt::lsdtt_raster_cpu(std::move(out_data), std::move(out_mask), meta);
}

// =============================================================================
// main
// =============================================================================
int main(int argc, char* argv[])
{
    print_banner();

    lsdtt::lsdtt_cli_args cli_args = lsdtt::lsdtt_param_parser::parse_args_deferred(
        argc, argv, k_tool_name, k_version, k_citation);

    // -------------------------------------------------------------------------
    // Parameter file
    // -------------------------------------------------------------------------
    
    // -------------------------------------------------------------------------
    // Parameter declarations
    // -------------------------------------------------------------------------
    std::map<std::string, bool>                        bool_default_map;
    std::map<std::string, float>                       float_default_map;
    std::map<std::string, int>                         int_default_map;
    std::map<std::string, std::string>                 string_default_map;
    lsdtt::lsdtt_help_registry    help_map;

    // --- Threading / routing ---
    int_default_map["n_threads"] = 0;
    help_map["n_threads"] = {"int","0",
        "Number of OpenMP threads.",
        "0 = use all available hardware threads. "
        "Applies to pit filling, flow routing, and polyfit.",
        "Threading / routing", "io", ""};

    string_default_map["boundary_conditions"] = "oooo";
    help_map["boundary_conditions"] = {"string","oooo",
        "Boundary conditions for flow routing.",
        "Four-character string for north, east, south, west. "
        "o = open, n = no-flux, p = periodic, b = forced outlet.",
        "Threading / routing", "param", ""};

    // --- DEM preparation ---
    lsdtt::declare_raster_is_filled_parameter(bool_default_map, help_map);

    float_default_map["fill_min_slope"] = 0.0001f;
    help_map["fill_min_slope"] = {"float","0.0001",
        "Minimum slope imposed by pit filling.",
        "Used only when raster_is_filled is false. "
        "Passed to fill_pits() to gently enforce downstream descent.",
        "DEM preparation", "param", ""};

    // --- MFD options ---
    float_default_map["mfd_exponent"] = 1.0f;
    help_map["mfd_exponent"] = {"float","1.0",
        "Slope-weighting exponent for MFD flow partitioning.",
        "Quinn et al. (1991) use 1.0; Freeman (1991) suggest 1.1. "
        "Only used when write_mfd_drainage_area_raster is true or "
        "twi_drainage_area_method is MFD.",
        "MFD options", "param", ""};

    // --- Drainage area outputs ---
    bool_default_map["build_flow_s_vector"] = true;
    help_map["build_flow_s_vector"] = {"bool","true",
        "Build the D8 s-vector during flow routing.",
        "The s-vector stores each cell's upslope set as a contiguous interval, "
        "which makes is_upslope_of an O(1) test and is what chi analysis, basin "
        "summaries, source finding and the cosmogenic routines rely on. It also "
        "costs about 57 percent of the D8 build (~106 ms of ~186 ms on a 1000x1000 "
        "tile) and roughly 24 bytes per cell of the flowinfo object's ~43. "
        "Drainage area does NOT need it: accumulating over the topological order "
        "produces a bit-identical raster for less work. Set false when drainage "
        "area is all you want, or to measure lsdtt3 against tools that build no "
        "equivalent structure. Ignored, with a warning, when another requested "
        "output needs it.",
        "Flow routing", "param", ""};

    bool_default_map["write_d8_drainage_area_raster"] = false;
    help_map["write_d8_drainage_area_raster"] = {"bool","false",
        "Write a D8 drainage area raster.",
        "Uses lsdtt_flowinfo (D8 single-receiver routing). "
        "Area is in square map units (m^2 for projected CRS).",
        "Drainage area outputs", "switch", ""};

    bool_default_map["write_dinf_drainage_area_raster"] = false;
    help_map["write_dinf_drainage_area_raster"] = {"bool","false",
        "Write a D-infinity drainage area raster.",
        "Uses Tarboton (1997) D-infinity routing via lsdtt_flowgraph_cpu. "
        "Distributed flow gives smoother hillslope area estimates than D8.",
        "Drainage area outputs", "switch", ""};

    bool_default_map["write_mfd_drainage_area_raster"] = false;
    help_map["write_mfd_drainage_area_raster"] = {"bool","false",
        "Write a multiple-flow-direction (MFD) drainage area raster.",
        "Uses Quinn et al. (1991) / Freeman (1991) MFD routing via "
        "lsdtt_flowgraph_cpu. Area is partitioned among all downslope "
        "neighbours weighted by slope^mfd_exponent.",
        "Drainage area outputs", "switch", ""};

    // --- Topographic Wetness Index ---
    bool_default_map["write_twi_raster"] = false;
    help_map["write_twi_raster"] = {"bool","false",
        "Write a Topographic Wetness Index (TWI) raster.",
        "TWI = ln( (A / pixel_area) / slope ), where A is the upstream "
        "drainage area (m^2), pixel_area = cell_size^2, and slope is the "
        "gradient magnitude (m/m) from a local quadratic fit. "
        "Only valid where area > 0 and slope > 0.",
        "Topographic Wetness Index", "switch", ""};

    string_default_map["twi_drainage_area_method"] = "Dinf";
    help_map["twi_drainage_area_method"] = {"string","Dinf",
        "Flow routing method used to compute drainage area for TWI.",
        "Allowed values: Dinf (D-infinity) or MFD. "
        "D-infinity (default) is generally preferred for hillslope TWI.",
        "Topographic Wetness Index", "param", "write_twi_raster"};

    float_default_map["twi_window_radius_metres"] = 100.0f;
    help_map["twi_window_radius_metres"] = {"float","100.0",
        "Polyfit window radius (metres) for slope computation in TWI.",
        "Passed to direct_polyfit(). Larger values smooth out local "
        "roughness; smaller values capture fine-scale slope variation.",
        "Topographic Wetness Index", "param", "write_twi_raster"};

    // --- Soil Topographic Index (STI) ---
    bool_default_map["write_sti_raster"] = false;
    help_map["write_sti_raster"] = {"bool","false",
        "Write STI raster using the ln(Ks*D) denominator term.",
        "Output file: <write_fname>_sti.tif. Requires clipped soil and land-use "
        "vector overlays and the SSKIB Ksat horizon CSV.",
        "Soil Topographic Index (STI)", "switch", ""};

    bool_default_map["write_sti_proxy_raster"] = false;
    help_map["write_sti_proxy_raster"] = {"bool","false",
        "Write STI proxy raster using a runoff-based SPR denominator term.",
        "Output file: <write_fname>_sti_proxy.tif. Uses SERIES_SPR.csv "
        "joined by normalised series code.",
        "Soil Topographic Index (STI)", "switch", ""};

    bool_default_map["write_sti_intermediate_rasters"] = false;
    help_map["write_sti_intermediate_rasters"] = {"bool","false",
        "Write STI debugging rasters (coded series, coded land use, term, valid mask).",
        "Useful for checking joins and nodata coverage.",
        "Soil Topographic Index (STI)", "switch", ""};

    bool_default_map["write_sti_lookup_csv"] = false;
    help_map["write_sti_lookup_csv"] = {"bool","false",
        "Write per-series/land-use lookup diagnostics CSV.",
        "CSV fields: series_code, land_use, n_profiles, ln_kd_median, "
        "spr_prc, spr_term, source.",
        "Soil Topographic Index (STI)", "switch", ""};

    bool_default_map["sti_resample_to_1m_before_calculation"] = false;
    help_map["sti_resample_to_1m_before_calculation"] = {"bool","false",
        "Resample DEM to 1 m (cubic) before all calculations in this run.",
        "When true, the input raster is GDAL-warped to 1 m spacing using cubic "
        "resampling before pit filling/flow/STI. Uses lsdtt reproject/warp "
        "machinery and writes the resampled raster to write_path.",
        "Soil Topographic Index (STI)", "switch", ""};

    bool_default_map["write_sti_95th_percentile_rasters"] = false;
    help_map["write_sti_95th_percentile_rasters"] = {"bool","false",
        "Write STI threshold rasters keeping only the top 5% of values.",
        "After each requested STI output is computed, evaluates the 95th "
        "percentile over valid pixels and writes a companion raster where "
        "only pixels >= P95 are retained; all others are set to nodata. "
        "Outputs: <write_fname>_sti_p95.tif and/or <write_fname>_sti_proxy_p95.tif.",
        "Soil Topographic Index (STI)", "switch", ""};

    bool_default_map["write_sti_p95_jenks_rasters"] = false;
    help_map["write_sti_p95_jenks_rasters"] = {"bool","false",
        "Write Jenks-classified STI category rasters from STI P95 outputs.",
        "Runs Jenks natural breaks on valid values in each P95 raster and writes "
        "class rasters with IDs 1..K. Requires write_sti_95th_percentile_rasters=true.",
        "Soil Topographic Index (STI)", "switch", ""};

    int_default_map["sti_p95_jenks_n_classes"] = 5;
    help_map["sti_p95_jenks_n_classes"] = {"int","5",
        "Number of Jenks classes for STI P95 categorisation.",
        "Used when write_sti_p95_jenks_rasters is true. If there are fewer "
        "distinct occupied histogram bins than requested classes, the effective "
        "number of classes is reduced automatically.",
        "Soil Topographic Index (STI)", "param", "write_sti_p95_jenks_rasters"};

    int_default_map["sti_p95_jenks_max_histogram_bins"] = 2048;
    help_map["sti_p95_jenks_max_histogram_bins"] = {"int","2048",
        "Maximum histogram bins used by Jenks fitting on P95 rasters.",
        "Higher values better preserve fine value structure but increase fit time. "
        "Used when write_sti_p95_jenks_rasters is true.",
        "Soil Topographic Index (STI)", "param", "write_sti_p95_jenks_rasters"};

    string_default_map["sti_soil_vector_path"] =
        "data/scotland_soils/hutton_series_code_min.fgb";
    string_default_map["sti_soil_layer_name"] = "hutton_series_code_min";
    string_default_map["sti_landuse_vector_path"] =
        "data/scotland_soils/sg_land_use_min.fgb";
    string_default_map["sti_landuse_layer_name"] = "sg_land_use_min";
    string_default_map["sti_ksat_csv_path"] =
        "data/scotland_soils/SSKIB_C_MED_KSAT.csv";
    string_default_map["sti_spr_csv_path"] =
        "data/scotland_soils/SERIES_SPR.csv";
    string_default_map["sti_xlsx_path"] =
        "data/scotland_soils/SSKIB_OCT18_V1.2.xlsx";
    string_default_map["sti_drainage_area_method"] = "Dinf";
    string_default_map["sti_term_mode"] = "ksd";
    string_default_map["sti_spr_transform"] = "ln_1_minus_spr";

    float_default_map["sti_window_radius_metres"] = 100.0f;
    float_default_map["sti_depth_cap_cm"] = 100.0f;
    float_default_map["sti_min_valid_depth_cm"] = 30.0f;
    float_default_map["sti_ks_unit_scale_to_m_per_day"] = 1.0f;

    // --- Connectivity index (Cavalli et al. 2013) ---
    bool_default_map["write_connectivity_raster"] = false;
    help_map["write_connectivity_raster"] = {"bool","false",
        "Write the sediment connectivity index IC = log10(D_up / D_down).",
        "Follows Cavalli et al. (2013). D8 routing is always used for "
        "D_up and D_down path calculations. The area raster method is "
        "controlled by connectivity_drainage_area_method.",
        "Connectivity index (Cavalli et al. 2013)", "switch", ""};

    string_default_map["connectivity_drainage_area_method"] = "d8";
    help_map["connectivity_drainage_area_method"] = {"string","d8",
        "Flow routing method for the drainage area fed into D_up.",
        "Allowed values: d8, dinf, mfd. D8 is the simplest; dinf and mfd "
        "give smoother hillslope area estimates. D8 routing is always used "
        "for the upslope mean and downstream path computations.",
        "Connectivity index (Cavalli et al. 2013)", "param", ""};

    float_default_map["connectivity_area_threshold_m2"] = 1.0e6f;
    help_map["connectivity_area_threshold_m2"] = {"float","1e6",
        "Map drainage-area threshold in square metres defining channel / sink pixels.",
        "This remains total map area in m^2 and is not a contributing-pixel or "
        "specific-area threshold. Pixels at or above it are treated as channels. "
        "D_down accumulation stops when the flow path reaches a channel.",
        "Connectivity index (Cavalli et al. 2013)", "param", ""};

    float_default_map["connectivity_window_radius_metres"] = 100.0f;
    help_map["connectivity_window_radius_metres"] = {"float","100.0",
        "Moving-window radius (metres) for slope and roughness computation.",
        "Used for the polyfit slope and for the roughness metric (MADk2 or "
        "IQR of residuals). Square window always used.",
        "Connectivity index (Cavalli et al. 2013)", "param", ""};

    string_default_map["connectivity_slope_method"] = "polyfit_direct";
    help_map["connectivity_slope_method"] = {"string","polyfit_direct",
        "Method for computing slope used in IC.",
        "Allowed values: polyfit_direct (fast, square window, no FFTW required), "
        "wavelet (FFT-based, requires FFTW). Both fit a local quadratic surface; "
        "slope magnitude = sqrt(dz/dx^2 + dz/dy^2).",
        "Connectivity index (Cavalli et al. 2013)", "param", ""};

    string_default_map["connectivity_roughness_method"] = "madk2";
    help_map["connectivity_roughness_method"] = {"string","madk2",
        "Roughness metric used to compute the W weighting factor.",
        "Allowed values: madk2 (Median Absolute Deviation of k2 differences, "
        "Trevisani & Guth 2025 - default), iqr_residual (IQR of polyfit "
        "residuals, square window). W = max(1 - R/R_max, 0.001).",
        "Connectivity index (Cavalli et al. 2013)", "param", ""};

    bool_default_map["write_height_above_channel_raster"] = false;
    help_map["write_height_above_channel_raster"] = {"bool","false",
        "Write a height-above-channel raster.",
        "For each pixel: elevation difference to the nearest downstream channel "
        "pixel along the D8 flow path. Channel pixels are those whose D8 "
        "contributing-pixel count >= threshold_contributing_pixels.",
        "Connectivity index (Cavalli et al. 2013)", "switch", ""};

    lsdtt::declare_channel_initiation_parameters(int_default_map, help_map);

    // --- Output format ---
    string_default_map["write_output_format"] = "cog";
    help_map["write_output_format"] = {"string","cog",
        "Output format for all raster outputs.",
        "Allowed values: cog, gtiff, same_as_input.",
        "Output format", "io", ""};

    // -------------------------------------------------------------------------
    // Parse
    // -------------------------------------------------------------------------
    // -------------------------------------------------------------------------
    // Help - rendered from the declarations above; no parameter file needed
    // -------------------------------------------------------------------------
    if (cli_args.help_requested)
    {
        lsdtt::lsdtt_param_parser::print_terminal_help(
            k_tool_name, k_version, k_citation, help_map, cli_args);
        return 0;
    }

    // -------------------------------------------------------------------------
    // Load parameter file
    // -------------------------------------------------------------------------
    lsdtt::lsdtt_param_parser pp;
    try
    {
        pp = lsdtt::lsdtt_param_parser(cli_args.path, cli_args.fname);
    }
    catch (const std::exception& e)
    {
        std::cerr << "Error loading parameter file: " << e.what() << "\n";
        return 1;
    }

    pp.set_deprecated_aliases({
        {"channel_threshold_contributing_pixels",
         "threshold_contributing_pixels"},
        {"twi_window_radius", "twi_window_radius_metres"},
        {"sti_window_radius", "sti_window_radius_metres"},
        {"connectivity_window_radius", "connectivity_window_radius_metres"} });

    pp.parse_all_parameters(float_default_map, int_default_map,
                            bool_default_map, string_default_map);

    auto bf  = pp.get_bool_parameters();
    auto ff  = pp.get_float_parameters();
    auto ifs = pp.get_int_parameters();
    auto sf  = pp.get_string_parameters();

    // -------------------------------------------------------------------------
    // Cross-parameter logic
    // -------------------------------------------------------------------------
    lsdtt::raster_output_format out_fmt =
        lsdtt::resolve_raster_output_format(sf["write_output_format"]);

    // Normalise method strings to lowercase for comparison.
    std::string twi_method = sf["twi_drainage_area_method"];
    twi_method = lsdtt::to_lower_ascii(twi_method);

    std::string conn_area_method = sf["connectivity_drainage_area_method"];
    conn_area_method = lsdtt::to_lower_ascii(conn_area_method);

    std::string conn_slope_method = sf["connectivity_slope_method"];
    conn_slope_method = lsdtt::to_lower_ascii(conn_slope_method);

    std::string conn_roughness_method = sf["connectivity_roughness_method"];
    conn_roughness_method = lsdtt::to_lower_ascii(conn_roughness_method);

    std::string sti_method = sf["sti_drainage_area_method"];
    sti_method = lsdtt::to_lower_ascii(sti_method);

    std::string sti_term_mode_str = sf["sti_term_mode"];
    sti_term_mode_str = lsdtt::to_lower_ascii(sti_term_mode_str);

    std::string sti_spr_transform_str = sf["sti_spr_transform"];
    sti_spr_transform_str = lsdtt::to_lower_ascii(sti_spr_transform_str);

    const bool write_sti = bf["write_sti_raster"];
    const bool write_sti_proxy = bf["write_sti_proxy_raster"];
    const bool write_sti_p95 = bf["write_sti_95th_percentile_rasters"];
    const bool write_sti_p95_jenks = bf["write_sti_p95_jenks_rasters"];
    const bool sti_resample_to_1m = bf["sti_resample_to_1m_before_calculation"];
    const bool need_sti = write_sti || write_sti_proxy;

    if (write_sti_p95 && !need_sti)
    {
        std::cout << "WARNING: write_sti_95th_percentile_rasters=true but no base STI output "
                  << "requested. Enable write_sti_raster and/or write_sti_proxy_raster.\n";
    }

    if (write_sti_p95_jenks && !write_sti_p95)
    {
        std::cerr << "FATAL: write_sti_p95_jenks_rasters=true requires "
                  << "write_sti_95th_percentile_rasters=true.\n";
        return 1;
    }

    if (write_sti_p95_jenks)
    {
        if (ifs["sti_p95_jenks_n_classes"] <= 1)
        {
            std::cerr << "FATAL: sti_p95_jenks_n_classes must be > 1. Got: "
                      << ifs["sti_p95_jenks_n_classes"] << "\n";
            return 1;
        }
        if (ifs["sti_p95_jenks_max_histogram_bins"] <= 0)
        {
            std::cerr << "FATAL: sti_p95_jenks_max_histogram_bins must be > 0. Got: "
                      << ifs["sti_p95_jenks_max_histogram_bins"] << "\n";
            return 1;
        }
    }

    if (sti_resample_to_1m && !need_sti)
    {
        std::cout << "WARNING: sti_resample_to_1m_before_calculation=true but no STI output "
                  << "requested. Resampling will still be applied to this run.\n";
    }

    if (bf["write_twi_raster"] && twi_method != "dinf" && twi_method != "mfd")
    {
        std::cerr << "FATAL: twi_drainage_area_method must be 'Dinf' or 'MFD'. "
                  << "Got: " << sf["twi_drainage_area_method"] << "\n";
        return 1;
    }

    if (bf["write_connectivity_raster"])
    {
        if (conn_area_method != "d8" && conn_area_method != "dinf" &&
            conn_area_method != "mfd")
        {
            std::cerr << "FATAL: connectivity_drainage_area_method must be "
                         "'d8', 'dinf', or 'mfd'. Got: "
                      << sf["connectivity_drainage_area_method"] << "\n";
            return 1;
        }
        if (conn_slope_method != "polyfit_direct" && conn_slope_method != "wavelet")
        {
            std::cerr << "FATAL: connectivity_slope_method must be "
                         "'polyfit_direct' or 'wavelet'. Got: "
                      << sf["connectivity_slope_method"] << "\n";
            return 1;
        }
        if (conn_roughness_method != "madk2" && conn_roughness_method != "iqr_residual")
        {
            std::cerr << "FATAL: connectivity_roughness_method must be "
                         "'madk2' or 'iqr_residual'. Got: "
                      << sf["connectivity_roughness_method"] << "\n";
            return 1;
        }
    }

    if (need_sti)
    {
        if (sti_method != "dinf" && sti_method != "mfd")
        {
            std::cerr << "FATAL: sti_drainage_area_method must be 'Dinf' or 'MFD'. "
                      << "Got: " << sf["sti_drainage_area_method"] << "\n";
            return 1;
        }

        if (sti_term_mode_str != "ksd" && sti_term_mode_str != "spr")
        {
            std::cerr << "FATAL: sti_term_mode must be 'ksd' or 'spr'. "
                      << "Got: " << sf["sti_term_mode"] << "\n";
            return 1;
        }

        if (sti_spr_transform_str != "ln_1_minus_spr" &&
            sti_spr_transform_str != "ln_spr")
        {
            std::cerr << "FATAL: sti_spr_transform must be 'ln_1_minus_spr' "
                      << "or 'ln_spr'. Got: " << sf["sti_spr_transform"] << "\n";
            return 1;
        }

        if (ff["sti_depth_cap_cm"] <= 0.0f)
        {
            std::cerr << "FATAL: sti_depth_cap_cm must be > 0. "
                      << "Got: " << ff["sti_depth_cap_cm"] << "\n";
            return 1;
        }

        if (ff["sti_min_valid_depth_cm"] <= 0.0f)
        {
            std::cerr << "FATAL: sti_min_valid_depth_cm must be > 0. "
                      << "Got: " << ff["sti_min_valid_depth_cm"] << "\n";
            return 1;
        }

        if (ff["sti_min_valid_depth_cm"] > ff["sti_depth_cap_cm"])
        {
            std::cerr << "FATAL: sti_min_valid_depth_cm cannot exceed "
                      << "sti_depth_cap_cm. Got min_valid_depth_cm="
                      << ff["sti_min_valid_depth_cm"] << ", depth_cap_cm="
                      << ff["sti_depth_cap_cm"] << "\n";
            return 1;
        }

        if (ff["sti_window_radius_metres"] <= 0.0f)
        {
            std::cerr << "FATAL: sti_window_radius_metres must be > 0. "
                      << "Got: " << ff["sti_window_radius_metres"] << "\n";
            return 1;
        }

        if (ff["sti_ks_unit_scale_to_m_per_day"] <= 0.0f)
        {
            std::cerr << "FATAL: sti_ks_unit_scale_to_m_per_day must be > 0. "
                      << "Got: " << ff["sti_ks_unit_scale_to_m_per_day"] << "\n";
            return 1;
        }

        if (write_sti)
        {
            if (sf["sti_soil_vector_path"].empty() || sf["sti_soil_layer_name"].empty())
            {
                std::cerr << "FATAL: STI requires non-empty sti_soil_vector_path "
                          << "and sti_soil_layer_name when write_sti_raster is true.\n";
                return 1;
            }
            if (sf["sti_landuse_vector_path"].empty() || sf["sti_landuse_layer_name"].empty())
            {
                std::cerr << "FATAL: STI requires non-empty sti_landuse_vector_path "
                          << "and sti_landuse_layer_name when write_sti_raster is true.\n";
                return 1;
            }
            if (sf["sti_ksat_csv_path"].empty())
            {
                std::cerr << "FATAL: STI ksd mode requires sti_ksat_csv_path.\n";
                return 1;
            }
        }

        if (write_sti_proxy)
        {
            if (sf["sti_soil_vector_path"].empty() || sf["sti_soil_layer_name"].empty())
            {
                std::cerr << "FATAL: STI proxy requires non-empty sti_soil_vector_path "
                          << "and sti_soil_layer_name when write_sti_proxy_raster is true.\n";
                return 1;
            }
            if (sf["sti_landuse_vector_path"].empty() || sf["sti_landuse_layer_name"].empty())
            {
                std::cerr << "FATAL: STI proxy requires non-empty sti_landuse_vector_path "
                          << "and sti_landuse_layer_name when write_sti_proxy_raster is true.\n";
                return 1;
            }
            if (sf["sti_spr_csv_path"].empty())
            {
                std::cerr << "FATAL: STI proxy mode requires sti_spr_csv_path.\n";
                return 1;
            }
        }
    }

    // -------------------------------------------------------------------------
    // Warn about unimplemented placeholders
    // -------------------------------------------------------------------------


    // -------------------------------------------------------------------------
    // Determine which flow objects are needed
    // -------------------------------------------------------------------------
    const bool need_d8   = bf["write_d8_drainage_area_raster"] ||
                           bf["write_connectivity_raster"]      ||   // D8 always needed for IC routing
                           bf["write_height_above_channel_raster"];
    const bool need_dinf = bf["write_dinf_drainage_area_raster"] ||
                            (bf["write_twi_raster"] && twi_method == "dinf") ||
                            (need_sti && sti_method == "dinf") ||
                            (bf["write_connectivity_raster"] && conn_area_method == "dinf");
    const bool need_mfd  = bf["write_mfd_drainage_area_raster"] ||
                            (bf["write_twi_raster"] && twi_method == "mfd") ||
                            (need_sti && sti_method == "mfd") ||
                            (bf["write_connectivity_raster"] && conn_area_method == "mfd");
    const bool need_flow = need_d8 || need_dinf || need_mfd;

    if (!need_flow)
    {
        std::cout << "No outputs requested. "
                  << "Enable at least one write_*_raster parameter.\n";
        pp.print_parameters();
        return 0;
    }

    // -------------------------------------------------------------------------
    // Open DEM
    // -------------------------------------------------------------------------
    std::string input_path = pp.get_full_read_path();
    std::cout << "\nOpening DEM: " << input_path << "\n";

    lsdtt::lsdtt_data_loader loader;
    try
    {
        loader.open_raster(input_path);
    }
    catch (const std::exception& e)
    {
        std::cerr << "Error opening DEM: " << e.what() << "\n";
        return 1;
    }

    loader.print_raster_info();

    if (sti_resample_to_1m)
    {
        std::cout << "\n--- Optional pre-calculation resampling ---\n";
        std::cout << "Resampling DEM to 1 m using cubic interpolation before analysis.\n";
        try
        {
            lsdtt::lsdtt_reproject_options ropt;
            ropt.enforce_projection = true;
            ropt.force_reproject = true;
            ropt.resampling_method = "cubic";
            ropt.force_resampling = true;
            ropt.forced_grid_spacing = 1.0;
            ropt.grid_spacing = -1.0;
            ropt.allow_lambert_projection = false;
            ropt.n_threads = ifs["n_threads"];
            ropt.warp_mem_mb = 512;

            const lsdtt::lsdtt_reproject_result res = lsdtt::lsdtt_reproject_if_needed(
                input_path,
                pp.get_write_path(),
                pp.get_write_fname() + "_sti_precalc_1m",
                loader.get_meta(),
                ropt
            );

            if (!res.reprojected || res.out_path.empty())
            {
                throw std::runtime_error(
                    "sti_resample_to_1m_before_calculation was requested but no output raster "
                    "was produced by the reprojection/warp step.");
            }

            input_path = res.out_path;
            std::cout << "Using resampled DEM: " << input_path << "\n";
            loader.open_raster(input_path);
            loader.print_raster_info();
        }
        catch (const std::exception& e)
        {
            std::cerr << "FATAL: 1 m pre-calculation resampling failed: "
                      << e.what() << "\n";
            return 1;
        }
    }

    try
    {
        loader.get_meta().check_projected_metres();
    }
    catch (const std::exception& e)
    {
        std::cerr << "\nFATAL: " << e.what() << "\n";
        return 1;
    }

    // -------------------------------------------------------------------------
    // Load into CPU raster
    // -------------------------------------------------------------------------
    lsdtt::lsdtt_raster_cpu dem(loader);
    lsdtt::lsdtt_raster_cpu routed_dem = dem;

    // -------------------------------------------------------------------------
    // Step 1 - Ensure the DEM is ready for routing
    // -------------------------------------------------------------------------
    if (!bf["raster_is_filled"])
    {
        std::cout << "\n--- Pit filling ---\n";
        std::cout << "Raster is not marked as filled. Running fill_pits().\n";
        std::cout << "Minimum slope: " << ff["fill_min_slope"] << "\n";

        try
        {
            routed_dem = lsdtt::fill_pits(
                dem,
                ff["fill_min_slope"],
                sf["boundary_conditions"],
                ifs["n_threads"]
            );
        }
        catch (const std::exception& e)
        {
            std::cerr << "FATAL: pit filling failed: " << e.what() << "\n";
            return 1;
        }
    }
    else
    {
        std::cout << "\nRaster is marked as already filled. "
                  << "Skipping fill step.\n";
    }

    // Convenience references to the filled DEM arrays used for flow routing.
    const lsdtt::lsdtt_array2d_f&  dem_data = routed_dem.get_data_array();
    const lsdtt::lsdtt_array2d_u8& dem_mask = routed_dem.get_mask_array();
    const lsdtt::lsdtt_raster_meta dem_meta = routed_dem.get_meta();

    // -------------------------------------------------------------------------
    // Step 2a - D8 flow routing (if needed)
    // -------------------------------------------------------------------------
    std::unique_ptr<lsdtt::lsdtt_flow_base> d8_flow;

    if (need_d8)
    {
        std::cout << "\n--- D8 flow routing ---\n";
        std::cout << "Boundary conditions: " << sf["boundary_conditions"] << "\n";
        std::cout << "Threads            : " << ifs["n_threads"] << "\n";

        lsdtt::lsdtt_flowinfo_options fi_opts;
        fi_opts.boundary_conditions       = sf["boundary_conditions"];
        fi_opts.interior_nodata_is_outlet = true;
        fi_opts.require_acyclic_graph     = true;
        fi_opts.build_topological_order   = true;

        // The s-vector is only optional when drainage area is the ONLY D8
        // product asked for. Connectivity routing and height-above-channel
        // both traverse upslope sets, so dropping it there would trade a
        // measurable saving for a runtime failure -- refuse rather than let
        // that happen, and say why.
        bool build_s_vector = bf["build_flow_s_vector"];
        if (!build_s_vector)
        {
            const bool needs_s_vector = bf["write_connectivity_raster"] ||
                                        bf["write_height_above_channel_raster"];
            if (needs_s_vector)
            {
                std::cout << "  NOTE: build_flow_s_vector: false was ignored - "
                             "connectivity and height-above-channel need the "
                             "s-vector. Building it.\n";
                build_s_vector = true;
            }
            else
            {
                std::cout << "  s-vector           : NOT built "
                             "(drainage area accumulates over the topological "
                             "order instead; identical result)\n";
            }
        }
        fi_opts.build_s_vector            = build_s_vector;
        fi_opts.n_threads                 = ifs["n_threads"];

        // Timed so the routing build can be separated from process startup
        // and GDAL I/O. The "Compute time:" wording matches
        // lsdtt-surface-metrics and is what
        // tests/scripts/run_benchmark_suite.py parses; the build and the
        // accumulation are reported separately because they are separately
        // meaningful and the reference implementations split them the same
        // way.
        const auto d8_build_t0 = std::chrono::steady_clock::now();
        try
        {
            d8_flow = lsdtt::make_flow_base(
                dem_data,
                dem_mask,
                dem_meta,
                fi_opts
            );
        }
        catch (const std::exception& e)
        {
            std::cerr << "FATAL: D8 flow routing failed: " << e.what() << "\n";
            return 1;
        }
        const auto d8_build_t1 = std::chrono::steady_clock::now();
        std::cout << "  Compute time: "
                  << std::chrono::duration<double>(
                         d8_build_t1 - d8_build_t0).count()
                  << " s\n";
    }

    // -------------------------------------------------------------------------
    // Step 2b - D-infinity flow routing (if needed)
    // -------------------------------------------------------------------------
    std::unique_ptr<lsdtt::lsdtt_flowgraph_cpu> dinf_graph;

    if (need_dinf)
    {
        std::cout << "\n--- D-infinity flow routing (Tarboton 1997) ---\n";
        std::cout << "Boundary conditions: " << sf["boundary_conditions"] << "\n";
        std::cout << "Threads            : " << ifs["n_threads"] << "\n";

        lsdtt::lsdtt_flowgraph_options fg_opts;
        fg_opts.boundary_conditions       = sf["boundary_conditions"];
        fg_opts.interior_nodata_is_outlet = true;
        fg_opts.method                    = lsdtt::lsdtt_flowgraph_method::d_infinity;
        fg_opts.n_threads                 = ifs["n_threads"];

        const auto dinf_build_t0 = std::chrono::steady_clock::now();
        try
        {
            dinf_graph = std::make_unique<lsdtt::lsdtt_flowgraph_cpu>(
                dem_data, dem_mask, dem_meta, fg_opts
            );
        }
        catch (const std::exception& e)
        {
            std::cerr << "FATAL: D-infinity flow routing failed: "
                      << e.what() << "\n";
            return 1;
        }
        const auto dinf_build_t1 = std::chrono::steady_clock::now();
        std::cout << "  Compute time: "
                  << std::chrono::duration<double>(
                         dinf_build_t1 - dinf_build_t0).count()
                  << " s\n";
    }

    // -------------------------------------------------------------------------
    // Step 2c - MFD flow routing (if needed)
    // -------------------------------------------------------------------------
    std::unique_ptr<lsdtt::lsdtt_flowgraph_cpu> mfd_graph;

    if (need_mfd)
    {
        std::cout << "\n--- MFD flow routing (Quinn et al. 1991 / Freeman 1991) ---\n";
        std::cout << "Boundary conditions: " << sf["boundary_conditions"] << "\n";
        std::cout << "MFD exponent       : " << ff["mfd_exponent"] << "\n";
        std::cout << "Threads            : " << ifs["n_threads"] << "\n";

        lsdtt::lsdtt_flowgraph_options fg_opts;
        fg_opts.boundary_conditions       = sf["boundary_conditions"];
        fg_opts.interior_nodata_is_outlet = true;
        fg_opts.method                    = lsdtt::lsdtt_flowgraph_method::mfd;
        fg_opts.mfd_exponent              = ff["mfd_exponent"];
        fg_opts.n_threads                 = ifs["n_threads"];

        const auto mfd_build_t0 = std::chrono::steady_clock::now();
        try
        {
            mfd_graph = std::make_unique<lsdtt::lsdtt_flowgraph_cpu>(
                dem_data, dem_mask, dem_meta, fg_opts
            );
        }
        catch (const std::exception& e)
        {
            std::cerr << "FATAL: MFD flow routing failed: "
                      << e.what() << "\n";
            return 1;
        }
        const auto mfd_build_t1 = std::chrono::steady_clock::now();
        std::cout << "  Compute time: "
                  << std::chrono::duration<double>(
                         mfd_build_t1 - mfd_build_t0).count()
                  << " s\n";
    }

    // -------------------------------------------------------------------------
    // Shared helper for TWI/STI paths
    // -------------------------------------------------------------------------
    auto compute_twi_raster =
        [&](const std::string& area_method,
            float              window_radius,
            const std::string& context_label) -> lsdtt::lsdtt_raster_cpu
    {
        std::cout << "\n--- " << context_label << " TWI computation ---\n";
        std::cout << "Drainage area method : " << area_method << "\n";
        std::cout << "Slope window radius  : " << window_radius << " m\n";

        lsdtt::poly_fit_coeffs coeffs;
        try
        {
            std::cout << "Computing slope via direct polyfit...\n";
            coeffs = lsdtt::direct_polyfit(
                routed_dem,
                window_radius,
                ifs["n_threads"]
            );
        }
        catch (const std::exception& e)
        {
            throw std::runtime_error(
                context_label + " polyfit for slope failed: " + std::string(e.what()));
        }

        lsdtt::polyfit_metric_request slope_req;
        slope_req.metrics   = { lsdtt::polyfit_metric::slope };
        slope_req.n_threads = ifs["n_threads"];

        lsdtt::polyfit_metric_result slope_res;
        try
        {
            slope_res = lsdtt::compute_polyfit_metrics(coeffs, slope_req);
        }
        catch (const std::exception& e)
        {
            throw std::runtime_error(
                context_label + " slope metric computation failed: " + std::string(e.what()));
        }

        const lsdtt::lsdtt_raster_cpu& slope_raster =
            slope_res.at(lsdtt::polyfit_metric::slope);

        std::string method_lc = area_method;
        method_lc = lsdtt::to_lower_ascii(method_lc);

        lsdtt::lsdtt_raster_result area_result;
        try
        {
            if (method_lc == "dinf")
            {
                area_result = dinf_graph->accumulate_raster();
            }
            else
            {
                area_result = mfd_graph->accumulate_raster();
            }
        }
        catch (const std::exception& e)
        {
            throw std::runtime_error(
                context_label + " drainage area accumulation failed: " +
                std::string(e.what()));
        }

        lsdtt::lsdtt_raster_cpu area_raster(std::move(area_result));

        const int    n_rows     = dem_meta.n_rows;
        const int    n_cols     = dem_meta.n_cols;
        const double pixel_area = dem_meta.cell_size * dem_meta.cell_size;
        const float  nodata_val = static_cast<float>(dem_meta.no_data_value);

        lsdtt::lsdtt_array2d_f  twi_data(n_rows, n_cols, nodata_val);
        lsdtt::lsdtt_array2d_u8 twi_mask(n_rows, n_cols, uint8_t(0));

        const lsdtt::lsdtt_array2d_f&  area_arr  = area_raster.get_data_array();
        const lsdtt::lsdtt_array2d_u8& area_msk  = area_raster.get_mask_array();
        const lsdtt::lsdtt_array2d_f&  slope_arr = slope_raster.get_data_array();
        const lsdtt::lsdtt_array2d_u8& slope_msk = slope_raster.get_mask_array();

        for (int row = 0; row < n_rows; ++row)
        {
            for (int col = 0; col < n_cols; ++col)
            {
                if (area_msk(row, col) == 0) continue;
                if (slope_msk(row, col) == 0) continue;

                const float area_val  = area_arr(row, col);
                const float slope_val = slope_arr(row, col);

                if (area_val > 0.0f && slope_val > 0.0f)
                {
                    twi_data(row, col) = std::log(
                        (area_val / static_cast<float>(pixel_area)) / slope_val
                    );
                    twi_mask(row, col) = uint8_t(1);
                }
            }
        }

        return lsdtt::lsdtt_raster_cpu(std::move(twi_data), std::move(twi_mask), dem_meta);
    };

    // -------------------------------------------------------------------------
    // Step 3a - Write D8 drainage area
    // -------------------------------------------------------------------------
    if (bf["write_d8_drainage_area_raster"])
    {
        std::cout << "\n--- Writing D8 drainage area ---\n";
        const std::string out_file = pp.make_output_path( "_d8_drainage_area", ".tif");
        std::cout << "Writing: " << out_file << "\n";

        try
        {
            const auto acc_t0 = std::chrono::steady_clock::now();
            lsdtt::lsdtt_raster_result area = d8_flow->get_drainage_area_raster();
            const auto acc_t1 = std::chrono::steady_clock::now();
            std::cout << "  Compute time: "
                      << std::chrono::duration<double>(acc_t1 - acc_t0).count()
                      << " s\n";

            loader.write_full(area, out_file, out_fmt);
            const auto acc_t2 = std::chrono::steady_clock::now();
            std::cout << "  Write time: "
                      << std::chrono::duration<double>(acc_t2 - acc_t1).count()
                      << " s\n";
        }
        catch (const std::exception& e)
        {
            std::cerr << "FATAL: writing D8 drainage area failed: "
                      << e.what() << "\n";
            return 1;
        }
    }

    // -------------------------------------------------------------------------
    // Step 3b - Write D-infinity drainage area
    // -------------------------------------------------------------------------
    if (bf["write_dinf_drainage_area_raster"])
    {
        std::cout << "\n--- Writing D-infinity drainage area ---\n";
        const std::string out_file = pp.make_output_path( "_dinf_drainage_area", ".tif");
        std::cout << "Writing: " << out_file << "\n";

        try
        {
            const auto acc_t0 = std::chrono::steady_clock::now();
            lsdtt::lsdtt_raster_result area = dinf_graph->accumulate_raster();
            const auto acc_t1 = std::chrono::steady_clock::now();
            std::cout << "  Compute time: "
                      << std::chrono::duration<double>(acc_t1 - acc_t0).count()
                      << " s\n";

            loader.write_full(area, out_file, out_fmt);
            const auto acc_t2 = std::chrono::steady_clock::now();
            std::cout << "  Write time: "
                      << std::chrono::duration<double>(acc_t2 - acc_t1).count()
                      << " s\n";
        }
        catch (const std::exception& e)
        {
            std::cerr << "FATAL: writing D-infinity drainage area failed: "
                      << e.what() << "\n";
            return 1;
        }
    }

    // -------------------------------------------------------------------------
    // Step 3c - Write MFD drainage area
    // -------------------------------------------------------------------------
    if (bf["write_mfd_drainage_area_raster"])
    {
        std::cout << "\n--- Writing MFD drainage area ---\n";
        const std::string out_file = pp.make_output_path( "_mfd_drainage_area", ".tif");
        std::cout << "Writing: " << out_file << "\n";

        try
        {
            const auto acc_t0 = std::chrono::steady_clock::now();
            lsdtt::lsdtt_raster_result area = mfd_graph->accumulate_raster();
            const auto acc_t1 = std::chrono::steady_clock::now();
            std::cout << "  Compute time: "
                      << std::chrono::duration<double>(acc_t1 - acc_t0).count()
                      << " s\n";

            loader.write_full(area, out_file, out_fmt);
            const auto acc_t2 = std::chrono::steady_clock::now();
            std::cout << "  Write time: "
                      << std::chrono::duration<double>(acc_t2 - acc_t1).count()
                      << " s\n";
        }
        catch (const std::exception& e)
        {
            std::cerr << "FATAL: writing MFD drainage area failed: "
                      << e.what() << "\n";
            return 1;
        }
    }

    // -------------------------------------------------------------------------
    // Step 4 - Topographic Wetness Index
    //
    // TWI = ln( (A / pixel_area) / slope )
    //
    //   A          - upstream drainage area in m^2 from D-infinity or MFD
    //   pixel_area - cell_size^2 (m^2)
    //   slope      - gradient magnitude (m/m) from a local quadratic fit
    //
    // Note: (A / pixel_area) gives the number of contributing pixels rather
    // than the standard specific contributing area (A / cell_size). This
    // matches the original lsdtt formula for continuity with earlier code.
    // -------------------------------------------------------------------------
    if (bf["write_twi_raster"])
    {
        lsdtt::lsdtt_raster_cpu twi_raster;
        try
        {
            twi_raster = compute_twi_raster(
                sf["twi_drainage_area_method"],
                ff["twi_window_radius_metres"],
                "TWI"
            );
        }
        catch (const std::exception& e)
        {
            std::cerr << "FATAL: " << e.what() << "\n";
            return 1;
        }

        const std::string twi_file = pp.make_output_path( "_twi", ".tif");
        std::cout << "Writing: " << twi_file << "\n";

        try
        {
            loader.write_full(
                twi_raster.get_data(),
                twi_raster.get_meta(),
                twi_file,
                out_fmt
            );
        }
        catch (const std::exception& e)
        {
            std::cerr << "FATAL: writing TWI raster failed: "
                      << e.what() << "\n";
            return 1;
        }
    }

    // -------------------------------------------------------------------------
    // Soil Topographic Index (STI)
    // -------------------------------------------------------------------------
    if (need_sti)
    {
        std::cout << "\n--- Soil Topographic Index (STI) ---\n";
        std::cout << "Requested outputs       : "
                  << (write_sti ? "_sti.tif " : "")
                  << (write_sti_proxy ? "_sti_proxy.tif " : "")
                  << (write_sti_p95 ? "(+ _sti_p95.tif / _sti_proxy_p95.tif)" : "")
                  << (write_sti_p95_jenks ? " (+ Jenks class rasters)" : "")
                  << "\n";
        std::cout << "STI drainage method     : " << sf["sti_drainage_area_method"] << "\n";
        std::cout << "STI slope window radius : " << ff["sti_window_radius_metres"] << " m\n";
        std::cout << "STI term mode           : " << sf["sti_term_mode"] << "\n";
        std::cout << "STI SPR transform       : " << sf["sti_spr_transform"] << "\n";
        std::cout << "\nWARNING: STI unit assumptions are not yet validated automatically.\n";
        std::cout << "WARNING: assuming HORZ_TOP/HORZ_BOTTOM are in cm and using "
                  << "sti_depth_cap_cm=" << ff["sti_depth_cap_cm"] << ".\n";
        std::cout << "WARNING: assuming COND_0_HYPRES is compatible with Ks proxy use.\n";
        std::cout << "WARNING: STI magnitudes should be treated as relative/proxy until units are verified.\n";
        std::cout << "WARNING: applying sti_ks_unit_scale_to_m_per_day="
                  << ff["sti_ks_unit_scale_to_m_per_day"] << " to COND_0_HYPRES before ln(Ks*D).\n";

        const std::string soil_clip_path = pp.make_output_path("_sti_soil_clip", ".fgb");
        const std::string land_clip_path = pp.make_output_path("_sti_landuse_clip", ".fgb");

        try
        {
            std::cout << "Clipping soil vector to DEM extent...\n";
            lsdtt::lsdtt_data_loader::clip_vector_to_raster_bbox(
                sf["sti_soil_vector_path"],
                soil_clip_path,
                dem_meta,
                lsdtt::vector_output_format::flatgeobuf
            );

            std::cout << "Clipping land-use vector to DEM extent...\n";
            lsdtt::lsdtt_data_loader::clip_vector_to_raster_bbox(
                sf["sti_landuse_vector_path"],
                land_clip_path,
                dem_meta,
                lsdtt::vector_output_format::flatgeobuf
            );
        }
        catch (const std::exception& e)
        {
            std::cerr << "FATAL: STI vector clipping failed: " << e.what() << "\n";
            return 1;
        }

        int soil_clip_count = -1;
        int land_clip_count = -1;
        try
        {
            soil_clip_count = lsdtt::lsdtt_data_loader::count_features_in_first_layer(soil_clip_path);
            land_clip_count = lsdtt::lsdtt_data_loader::count_features_in_first_layer(land_clip_path);
        }
        catch (const std::exception& e)
        {
            std::cerr << "FATAL: STI clipped feature counting failed: " << e.what() << "\n";
            return 1;
        }

        if (soil_clip_count <= 0)
        {
            std::cerr << "FATAL: clipped STI soil vector has zero features over DEM extent.\n"
                      << "Check sti_soil_vector_path/layer and DEM overlap.\n";
            return 1;
        }
        if (land_clip_count <= 0)
        {
            std::cerr << "FATAL: clipped STI land-use vector has zero features over DEM extent.\n"
                      << "Check sti_landuse_vector_path/layer and DEM overlap.\n";
            return 1;
        }

        std::cout << "Clipped soil features    : " << soil_clip_count << "\n";
        std::cout << "Clipped land-use features: " << land_clip_count << "\n";

        lsdtt::sti_polygon_string_rasterise_result series_codes;
        lsdtt::sti_polygon_string_rasterise_result land_use_codes;
        try
        {
            series_codes = lsdtt::rasterise_polygon_string_field(
                routed_dem,
                soil_clip_path,
                sf["sti_soil_layer_name"],
                "SERIES_CODE",
                8
            );
        }
        catch (const std::exception&)
        {
            try
            {
                series_codes = lsdtt::rasterise_polygon_string_field(
                    routed_dem,
                    soil_clip_path,
                    sf["sti_soil_layer_name"],
                    "SERCDE1",
                    8
                );
                std::cout << "STI note: soil field 'SERIES_CODE' not found; used fallback field 'SERCDE1'.\n";
            }
            catch (const std::exception& e)
            {
                std::cerr << "FATAL: STI soil rasterization failed for fields SERIES_CODE and SERCDE1: "
                          << e.what() << "\n";
                return 1;
            }
        }

        try
        {
            land_use_codes = lsdtt::rasterise_polygon_string_field(
                routed_dem,
                land_clip_path,
                sf["sti_landuse_layer_name"],
                "LAND_USE",
                8
            );
        }
        catch (const std::exception& e)
        {
            std::cerr << "FATAL: STI land-use rasterization failed for field LAND_USE: "
                      << e.what() << "\n";
            return 1;
        }

        std::cout << "STI soil categories rasterized   : " << series_codes.lookup.size() << "\n";
        std::cout << "STI land-use categories rasterized: " << land_use_codes.lookup.size() << "\n";

        const bool xlsx_ok = lsdtt::sti_check_xlsx_available(sf["sti_xlsx_path"]);
        if (!xlsx_ok)
        {
            std::cout << "WARNING: sti_xlsx_path not found/readable: " << sf["sti_xlsx_path"]
                      << " (continuing; xlsx consistency checks are lightweight in this first pass).\n";
        }
        else
        {
            std::cout << "STI xlsx consistency source found: " << sf["sti_xlsx_path"] << "\n";
        }

        lsdtt::sti_lnkd_table lnkd_table;
        lsdtt::sti_spr_table spr_table;

        if (write_sti || bf["write_sti_lookup_csv"])
        {
            try
            {
                lsdtt::sti_lnkd_options lnkd_opts;
                lnkd_opts.depth_cap_cm = ff["sti_depth_cap_cm"];
                lnkd_opts.min_valid_depth_cm = ff["sti_min_valid_depth_cm"];
                lnkd_opts.ks_unit_scale_to_m_per_day = ff["sti_ks_unit_scale_to_m_per_day"];

                const lsdtt::sti_lnkd_build_result lnkd =
                    lsdtt::build_sti_lnkd_table_from_csv(sf["sti_ksat_csv_path"], lnkd_opts);
                lnkd_table = lnkd.table;

                std::cout << "STI ksat rows read/used        : " << lnkd.rows_read
                          << " / " << lnkd.rows_used << "\n";
                std::cout << "STI profiles total/valid       : " << lnkd.profiles_total
                          << " / " << lnkd.profiles_valid << "\n";
                std::cout << "STI series-landuse pairs built : " << lnkd.pairs_built << "\n";
            }
            catch (const std::exception& e)
            {
                std::cerr << "FATAL: STI ln(Ks*D) table build failed: " << e.what() << "\n";
                return 1;
            }
        }

        if (write_sti_proxy || bf["write_sti_lookup_csv"])
        {
            try
            {
                const lsdtt::sti_spr_parse_result spr =
                    lsdtt::load_sti_spr_table_from_csv(sf["sti_spr_csv_path"]);
                spr_table = spr.table;
                std::cout << "STI SPR rows read/used         : " << spr.rows_read
                          << " / " << spr.rows_used << "\n";
                std::cout << "STI SPR clamped rows           : " << spr.clamped_rows << "\n";
            }
            catch (const std::exception& e)
            {
                std::cerr << "FATAL: STI SPR table load failed: " << e.what() << "\n";
                return 1;
            }
        }

        lsdtt::lsdtt_raster_cpu sti_twi;
        try
        {
            sti_twi = compute_twi_raster(
                sf["sti_drainage_area_method"],
                ff["sti_window_radius_metres"],
                "STI"
            );
        }
        catch (const std::exception& e)
        {
            std::cerr << "FATAL: " << e.what() << "\n";
            return 1;
        }

        lsdtt::sti_spr_transform spr_transform = lsdtt::sti_spr_transform::ln_1_minus_spr;
        if (sti_spr_transform_str == "ln_spr")
        {
            spr_transform = lsdtt::sti_spr_transform::ln_spr;
        }

        if (bf["write_sti_intermediate_rasters"])
        {
            try
            {
                loader.write_full(
                    sti_twi.get_data(),
                    sti_twi.get_meta(),
                    pp.make_output_path("_sti_twi", ".tif"),
                    out_fmt
                );
            }
            catch (const std::exception& e)
            {
                std::cerr << "FATAL: writing STI TWI intermediate raster failed: "
                          << e.what() << "\n";
                return 1;
            }
        }

        std::unordered_map<std::string, std::string> lookup_sources;

        if (write_sti)
        {
            lsdtt::sti_compute_options compute_opts;
            compute_opts.mode = lsdtt::sti_term_mode::ksd;
            compute_opts.spr_transform = spr_transform;

            lsdtt::sti_compute_result ksd_result;
            try
            {
                ksd_result = lsdtt::compute_sti_from_rasters(
                    sti_twi,
                    series_codes,
                    land_use_codes,
                    lnkd_table,
                    spr_table,
                    compute_opts
                );
            }
            catch (const std::exception& e)
            {
                std::cerr << "FATAL: STI ksd assembly failed: " << e.what() << "\n";
                return 1;
            }

            const std::string sti_file = pp.make_output_path("_sti", ".tif");
            std::cout << "Writing: " << sti_file << "\n";
            try
            {
                loader.write_full(
                    ksd_result.sti_raster.get_data(),
                    ksd_result.sti_raster.get_meta(),
                    sti_file,
                    out_fmt
                );
            }
            catch (const std::exception& e)
            {
                std::cerr << "FATAL: writing STI raster failed: " << e.what() << "\n";
                return 1;
            }

            std::cout << "STI ksd valid output pixels     : " << ksd_result.stats.valid_output_pixels
                      << " / " << ksd_result.stats.valid_twi_pixels << "\n";
            std::cout << "STI ksd missing lookup pixels   : " << ksd_result.stats.missing_lookup_pixels << "\n";
            std::cout << "STI ksd fallback pixels         : " << ksd_result.stats.fallback_pixels << "\n";

            if (write_sti_p95)
            {
                try
                {
                    float p95_threshold = std::numeric_limits<float>::quiet_NaN();
                    std::int64_t kept_pixels = 0;
                    lsdtt::lsdtt_raster_cpu sti_p95 = keep_top_percentile_pixels(
                        ksd_result.sti_raster, 95.0, p95_threshold, kept_pixels);
                    const std::string sti_p95_file = pp.make_output_path("_sti_p95", ".tif");
                    std::cout << "Writing: " << sti_p95_file << "\n";
                    loader.write_full(
                        sti_p95.get_data(),
                        sti_p95.get_meta(),
                        sti_p95_file,
                        out_fmt
                    );
                    std::cout << "STI ksd p95 threshold value     : " << p95_threshold << "\n";
                    std::cout << "STI ksd p95 retained pixels     : " << kept_pixels
                              << " / " << ksd_result.stats.valid_output_pixels << "\n";

                    if (write_sti_p95_jenks)
                    {
                        lsdtt::algebra::raster_jenks_classification_result classification =
                            lsdtt::algebra::classify_jenks_categories(
                            sti_p95,
                            ifs["sti_p95_jenks_n_classes"],
                            ifs["sti_p95_jenks_max_histogram_bins"]
                        );
                        const std::string sti_p95_jenks_file = pp.make_output_path(
                            "_sti_p95_jenks" + std::to_string(ifs["sti_p95_jenks_n_classes"]),
                            ".tif"
                        );
                        std::cout << "Writing: " << sti_p95_jenks_file << "\n";
                        loader.write_full(
                            classification.classified_raster.get_data(),
                            classification.classified_raster.get_meta(),
                            sti_p95_jenks_file,
                            out_fmt
                        );

                        std::cout << "STI ksd Jenks classes requested/effective: "
                                  << classification.jenks.requested_classes << " / "
                                  << classification.jenks.effective_classes << "\n";
                        std::cout << "STI ksd Jenks classified pixels          : "
                                  << classification.classified_pixels << " / " << kept_pixels << "\n";
                        std::cout << "STI ksd Jenks GVF                        : "
                                  << classification.jenks.goodness_of_variance_fit << "\n";
                        std::cout << "STI ksd Jenks upper breaks               : ";
                        for (std::size_t bi = 0; bi < classification.jenks.upper_breaks.size(); ++bi)
                        {
                            if (bi > 0) std::cout << ", ";
                            std::cout << classification.jenks.upper_breaks[bi];
                        }
                        std::cout << "\n";
                    }
                }
                catch (const std::exception& e)
                {
                    std::cerr << "FATAL: writing STI p95 raster failed: " << e.what() << "\n";
                    return 1;
                }
            }

            lookup_sources.insert(ksd_result.source_by_series_landuse.begin(),
                                  ksd_result.source_by_series_landuse.end());

            if (bf["write_sti_intermediate_rasters"])
            {
                try
                {
                    loader.write_full(
                        series_codes.raster_codes.get_data(),
                        series_codes.raster_codes.get_meta(),
                        pp.make_output_path("_sti_series_codes", ".tif"),
                        out_fmt
                    );
                    loader.write_full(
                        land_use_codes.raster_codes.get_data(),
                        land_use_codes.raster_codes.get_meta(),
                        pp.make_output_path("_sti_landuse_codes", ".tif"),
                        out_fmt
                    );
                    loader.write_full(
                        ksd_result.term_raster.get_data(),
                        ksd_result.term_raster.get_meta(),
                        pp.make_output_path("_sti_term_ksd", ".tif"),
                        out_fmt
                    );
                    loader.write_full(
                        ksd_result.valid_mask_raster.get_data(),
                        ksd_result.valid_mask_raster.get_meta(),
                        pp.make_output_path("_sti_valid_mask_ksd", ".tif"),
                        out_fmt
                    );
                }
                catch (const std::exception& e)
                {
                    std::cerr << "FATAL: writing STI intermediate rasters failed: "
                              << e.what() << "\n";
                    return 1;
                }
            }
        }

        if (write_sti_proxy)
        {
            lsdtt::sti_compute_options compute_opts;
            compute_opts.mode = lsdtt::sti_term_mode::spr;
            compute_opts.spr_transform = spr_transform;

            lsdtt::sti_compute_result spr_result;
            try
            {
                spr_result = lsdtt::compute_sti_from_rasters(
                    sti_twi,
                    series_codes,
                    land_use_codes,
                    lnkd_table,
                    spr_table,
                    compute_opts
                );
            }
            catch (const std::exception& e)
            {
                std::cerr << "FATAL: STI proxy assembly failed: " << e.what() << "\n";
                return 1;
            }

            const std::string sti_proxy_file = pp.make_output_path("_sti_proxy", ".tif");
            std::cout << "Writing: " << sti_proxy_file << "\n";
            try
            {
                loader.write_full(
                    spr_result.sti_raster.get_data(),
                    spr_result.sti_raster.get_meta(),
                    sti_proxy_file,
                    out_fmt
                );
            }
            catch (const std::exception& e)
            {
                std::cerr << "FATAL: writing STI proxy raster failed: " << e.what() << "\n";
                return 1;
            }

            std::cout << "STI proxy valid output pixels   : " << spr_result.stats.valid_output_pixels
                      << " / " << spr_result.stats.valid_twi_pixels << "\n";
            std::cout << "STI proxy missing lookup pixels : " << spr_result.stats.missing_lookup_pixels << "\n";
            std::cout << "STI proxy clamped pixels        : " << spr_result.stats.clamped_pixels << "\n";

            if (write_sti_p95)
            {
                try
                {
                    float p95_threshold = std::numeric_limits<float>::quiet_NaN();
                    std::int64_t kept_pixels = 0;
                    lsdtt::lsdtt_raster_cpu sti_proxy_p95 = keep_top_percentile_pixels(
                        spr_result.sti_raster, 95.0, p95_threshold, kept_pixels);
                    const std::string sti_proxy_p95_file = pp.make_output_path("_sti_proxy_p95", ".tif");
                    std::cout << "Writing: " << sti_proxy_p95_file << "\n";
                    loader.write_full(
                        sti_proxy_p95.get_data(),
                        sti_proxy_p95.get_meta(),
                        sti_proxy_p95_file,
                        out_fmt
                    );
                    std::cout << "STI proxy p95 threshold value   : " << p95_threshold << "\n";
                    std::cout << "STI proxy p95 retained pixels   : " << kept_pixels
                              << " / " << spr_result.stats.valid_output_pixels << "\n";

                    if (write_sti_p95_jenks)
                    {
                        lsdtt::algebra::raster_jenks_classification_result classification =
                            lsdtt::algebra::classify_jenks_categories(
                            sti_proxy_p95,
                            ifs["sti_p95_jenks_n_classes"],
                            ifs["sti_p95_jenks_max_histogram_bins"]
                        );
                        const std::string sti_proxy_p95_jenks_file = pp.make_output_path(
                            "_sti_proxy_p95_jenks" + std::to_string(ifs["sti_p95_jenks_n_classes"]),
                            ".tif"
                        );
                        std::cout << "Writing: " << sti_proxy_p95_jenks_file << "\n";
                        loader.write_full(
                            classification.classified_raster.get_data(),
                            classification.classified_raster.get_meta(),
                            sti_proxy_p95_jenks_file,
                            out_fmt
                        );

                        std::cout << "STI proxy Jenks classes requested/effective: "
                                  << classification.jenks.requested_classes << " / "
                                  << classification.jenks.effective_classes << "\n";
                        std::cout << "STI proxy Jenks classified pixels          : "
                                  << classification.classified_pixels << " / " << kept_pixels << "\n";
                        std::cout << "STI proxy Jenks GVF                        : "
                                  << classification.jenks.goodness_of_variance_fit << "\n";
                        std::cout << "STI proxy Jenks upper breaks               : ";
                        for (std::size_t bi = 0; bi < classification.jenks.upper_breaks.size(); ++bi)
                        {
                            if (bi > 0) std::cout << ", ";
                            std::cout << classification.jenks.upper_breaks[bi];
                        }
                        std::cout << "\n";
                    }
                }
                catch (const std::exception& e)
                {
                    std::cerr << "FATAL: writing STI proxy p95 raster failed: " << e.what() << "\n";
                    return 1;
                }
            }

            lookup_sources.insert(spr_result.source_by_series_landuse.begin(),
                                  spr_result.source_by_series_landuse.end());

            if (bf["write_sti_intermediate_rasters"])
            {
                try
                {
                    loader.write_full(
                        series_codes.raster_codes.get_data(),
                        series_codes.raster_codes.get_meta(),
                        pp.make_output_path("_sti_series_codes", ".tif"),
                        out_fmt
                    );
                    loader.write_full(
                        land_use_codes.raster_codes.get_data(),
                        land_use_codes.raster_codes.get_meta(),
                        pp.make_output_path("_sti_landuse_codes", ".tif"),
                        out_fmt
                    );
                    loader.write_full(
                        spr_result.term_raster.get_data(),
                        spr_result.term_raster.get_meta(),
                        pp.make_output_path("_sti_term_spr", ".tif"),
                        out_fmt
                    );
                    loader.write_full(
                        spr_result.valid_mask_raster.get_data(),
                        spr_result.valid_mask_raster.get_meta(),
                        pp.make_output_path("_sti_valid_mask_spr", ".tif"),
                        out_fmt
                    );
                }
                catch (const std::exception& e)
                {
                    std::cerr << "FATAL: writing STI intermediate rasters failed: "
                              << e.what() << "\n";
                    return 1;
                }
            }
        }

        if (bf["write_sti_lookup_csv"])
        {
            try
            {
                const std::vector<lsdtt::sti_lookup_row> rows =
                    lsdtt::build_sti_lookup_rows(
                        lookup_sources,
                        lnkd_table,
                        spr_table,
                        spr_transform
                    );
                const std::string lookup_csv = pp.make_output_path("_sti_lookup", ".csv");
                std::cout << "Writing: " << lookup_csv << "\n";
                lsdtt::write_sti_lookup_csv(rows, lookup_csv);
                std::cout << "STI lookup rows written         : " << rows.size() << "\n";
            }
            catch (const std::exception& e)
            {
                std::cerr << "FATAL: writing STI lookup CSV failed: " << e.what() << "\n";
                return 1;
            }
        }
    }

    // -------------------------------------------------------------------------
    // Connectivity index (Cavalli et al. 2013)
    //   IC = log10(D_up / D_down)
    //   D_up  = mean_W_upslope * mean_S_upslope * sqrt(A)
    //   D_down = sum d_i/(W_i*S_i) along D8 path to nearest channel
    // -------------------------------------------------------------------------
    if (bf["write_connectivity_raster"])
    {
        std::cout << "\n--- Connectivity index (Cavalli et al. 2013) ---\n";
        std::cout << "Area method      : " << sf["connectivity_drainage_area_method"] << "\n";
        std::cout << "Slope method     : " << sf["connectivity_slope_method"] << "\n";
        std::cout << "Roughness method : " << sf["connectivity_roughness_method"] << "\n";
        std::cout << "Window radius    : " << ff["connectivity_window_radius_metres"] << " m\n";
        std::cout << "Area threshold   : " << ff["connectivity_area_threshold_m2"] << " m^2\n";

        const double   conn_radius = static_cast<double>(ff["connectivity_window_radius_metres"]);
        const float    conn_thresh = ff["connectivity_area_threshold_m2"];

        try
        {
            // -- Drainage area for IC --
            lsdtt::lsdtt_raster_result conn_area;
            if (conn_area_method == "d8")
            {
                conn_area = d8_flow->get_drainage_area_raster();
            }
            else if (conn_area_method == "dinf")
            {
                conn_area = dinf_graph->accumulate_raster();
            }
            else  // mfd
            {
                conn_area = mfd_graph->accumulate_raster();
            }

            // -- Polyfit coefficients (needed for slope and optionally residual) --
            // Computed once and reused if both slope and roughness need it.
            const bool need_polyfit = conn_slope_method == "polyfit_direct" ||
                                      conn_roughness_method == "iqr_residual";
            const bool need_wavelet = conn_slope_method == "wavelet";

            lsdtt::poly_fit_coeffs conn_coeffs;
            if (need_polyfit)
            {
                std::cout << "Computing direct polyfit (square window)...\n";
                conn_coeffs = lsdtt::direct_polyfit(dem, conn_radius, ifs["n_threads"]);
            }

            lsdtt::poly_fit_coeffs conn_wav_coeffs;
            if (need_wavelet)
            {
                std::cout << "Computing wavelet quadratic fit...\n";
                conn_wav_coeffs = lsdtt::wavelet_quadratic_fit(
                    dem, conn_radius, ifs["n_threads"]);
            }

            // -- Slope --
            lsdtt::polyfit_metric_request slope_req;
            slope_req.metrics   = { lsdtt::polyfit_metric::slope };
            slope_req.n_threads = ifs["n_threads"];

            lsdtt::polyfit_metric_result slope_res =
                lsdtt::compute_polyfit_metrics(
                    conn_slope_method == "wavelet" ? conn_wav_coeffs : conn_coeffs,
                    slope_req);

            const lsdtt::lsdtt_raster_cpu& conn_slope =
                slope_res.at(lsdtt::polyfit_metric::slope);

            // -- Roughness for W --
            lsdtt::lsdtt_raster_cpu conn_roughness;
            if (conn_roughness_method == "madk2")
            {
                std::cout << "Computing MADk2 roughness...\n";
                lsdtt::madk2_result mk2 = lsdtt::calculate_madk2(
                    dem, conn_radius, 1, ifs["n_threads"]);
                conn_roughness = std::move(mk2.omnidirectional_roughness);
            }
            else  // iqr_residual
            {
                std::cout << "Computing IQR of polyfit residuals...\n";
                lsdtt::lsdtt_raster_cpu residual =
                    lsdtt::algebra::raster_subtract(dem, conn_coeffs.f);
                conn_roughness = lsdtt::moving_window_iqr(
                    residual, conn_radius,
                    lsdtt::window_shape::square,
                    ifs["n_threads"]);
            }

            // -- Compute IC --
            std::cout << "Computing connectivity index...\n";
            lsdtt::lsdtt_raster_result ic = lsdtt::compute_connectivity_index(
                *d8_flow, conn_roughness, conn_slope, conn_area, conn_thresh);

            const std::string ic_file = pp.make_output_path( "_connectivity_index", ".tif");
            std::cout << "Writing: " << ic_file << "\n";
            loader.write_full(ic, ic_file, out_fmt);

            // Optionally write D_up and D_down for debugging
            lsdtt::lsdtt_raster_result W_result =
                lsdtt::compute_W_raster(conn_roughness);

            const std::string dup_file = pp.make_output_path( "_connectivity_d_up", ".tif");
            std::cout << "Writing: " << dup_file << "\n";
            loader.write_full(
                lsdtt::compute_D_up_raster(*d8_flow, W_result, conn_slope, conn_area),
                dup_file, out_fmt);

            const std::string ddn_file = pp.make_output_path( "_connectivity_d_down", ".tif");
            std::cout << "Writing: " << ddn_file << "\n";
            loader.write_full(
                lsdtt::compute_D_down_raster(*d8_flow, W_result, conn_slope,
                                             conn_area, conn_thresh),
                ddn_file, out_fmt);
        }
        catch (const std::exception& e)
        {
            std::cerr << "FATAL: connectivity index computation failed: "
                      << e.what() << "\n";
            return 1;
        }
    }

    // -------------------------------------------------------------------------
    // Height above channel
    // -------------------------------------------------------------------------
    if (bf["write_height_above_channel_raster"])
    {
        std::cout << "\n--- Height above channel ---\n";
        std::cout << "Channel threshold: "
                  << ifs["threshold_contributing_pixels"]
                  << " contributing pixels\n";

        try
        {
            lsdtt::lsdtt_raster_result hac = d8_flow->get_height_above_channel_raster(
                routed_dem.get_data(),
                routed_dem.get_mask(),
                ifs["threshold_contributing_pixels"]
            );

            const std::string hac_file =
                pp.make_output_path( "_height_above_channel", ".tif");
            std::cout << "Writing: " << hac_file << "\n";
            loader.write_full(hac, hac_file, out_fmt);
        }
        catch (const std::exception& e)
        {
            std::cerr << "FATAL: height above channel failed: "
                      << e.what() << "\n";
            return 1;
        }
    }

    // -------------------------------------------------------------------------
    // Parameter echo
    // -------------------------------------------------------------------------
    pp.print_help(help_map, pp.get_write_fname() + "_README",
                  k_version, k_citation);
    pp.print_parameters();

    std::cout << "\nlsdtt-drainage-analysis: done.\n";
    return 0;
}
