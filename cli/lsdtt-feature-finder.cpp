// lsdtt-feature-finder.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2025 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html
//
// -----------------------------------------------------------------------------
// lsdtt-feature-finder
//
// Experimental feature-detection CLI for LiDAR-derived DEM analysis.
//
// Current first-draft functionality:
//   1. channel_bottoms
//
// The intent of the first feature mode is to build a binary raster of likely
// natural channel bottoms so that later workflows can preferentially connect
// those features across obstructions such as roads and culverts.
//
// Current channel-bottom workflow:
//   1. Compute polyfit curvature and slope at a small wavelength
//   2. Compute polyfit curvature at a larger wavelength
//   3. Convert curvature to a positive channelness signal using a configurable
//      sign convention
//   4. Normalise both signals to [0, 1]
//   5. Compute a multiscale persistence raster using pixelwise min()
//   6. Compute a local-depth signal using moving_window_mean(dem) - dem
//   7. Optionally include a weak slope term
//   8. Form a weighted score raster
//   9. Threshold to a binary channel-bottom mask
//
// Usage:
//   lsdtt-feature-finder <path> <parameter_file>
//   lsdtt-feature-finder <parameter_file>   (path assumed to be ./)
//   lsdtt-feature-finder -h                 (print help and exit)
//   lsdtt-feature-finder -v                 (print version and exit)
//   lsdtt-feature-finder -c                 (print citation and exit)
// -----------------------------------------------------------------------------

#include <algorithm>
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
#include "../src/cpu/lsdtt_raster_cpu.hpp"
#include "../src/cpu/lsdtt_raster_spectral_cpu.hpp"
#include "../src/cpu/lsdtt_raster_surface_metrics_cpu.hpp"

static const std::string k_version  = "0.2.0";
static const std::string k_citation = "https://github.com/LSDtopotools/lsdtt3";

// -----------------------------------------------------------------------------
// print_banner
// -----------------------------------------------------------------------------
static void print_banner()
{
    std::cout
        << "=========================================================\n"
        << "|| Welcome to lsdtt-feature-finder                     ||\n"
        << "|| I am going to look for terrain features for you     ||\n"
        << "|| Part of lsdtt3 by Simon M. Mudd                     ||\n"
        << "=========================================================\n"
        << "|| Version : " << k_version  << "\n"
        << "|| Citation: " << k_citation << "\n"
        << "=========================================================\n\n";
}

// -----------------------------------------------------------------------------
// parse_args
// -----------------------------------------------------------------------------
static std::pair<std::string, std::string> parse_args(int argc, char* argv[])
{
    auto is_flag = [](const std::string& s,
                      std::initializer_list<std::string> flags)
    {
        for (const auto& f : flags)
        {
            if (s == f) return true;
        }
        return false;
    };

    for (int i = 1; i < argc; ++i)
    {
        const std::string a = argv[i];

        if (is_flag(a, {"-h", "--help", "help", "-H", "--H"}))
        {
            std::cout
                << "Usage: lsdtt-feature-finder <path> <parameter_file>\n"
                << "       lsdtt-feature-finder <parameter_file>\n"
                << "       lsdtt-feature-finder -h   (this help)\n"
                << "       lsdtt-feature-finder -v   (version)\n"
                << "       lsdtt-feature-finder -c   (citation)\n";
            std::exit(0);
        }

        if (is_flag(a, {"-v", "--version", "version", "-V"}))
        {
            std::cout << "lsdtt-feature-finder version " << k_version << "\n";
            std::exit(0);
        }

        if (is_flag(a, {"-c", "--cite", "citation", "-citation", "--citation"}))
        {
            std::cout << "Citation: " << k_citation << "\n";
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
    else if (argc == 2)
    {
        fname = argv[1];
    }
    else
    {
        path  = argv[1];
        fname = argv[2];
    }

    return {path, fname};
}

// -----------------------------------------------------------------------------
// out_path
// -----------------------------------------------------------------------------
static std::string out_path(
    const lsdtt::lsdtt_param_parser& pp,
    const std::string&               suffix)
{
    return pp.get_full_write_path() + suffix;
}

// -----------------------------------------------------------------------------
// resolve_output_format
// -----------------------------------------------------------------------------
static lsdtt::raster_output_format resolve_output_format(
    const std::string& fmt_str)
{
    if (fmt_str == "gtiff")
    {
        return lsdtt::raster_output_format::gtiff;
    }
    if (fmt_str == "same_as_input")
    {
        return lsdtt::raster_output_format::same_as_input;
    }
    if (fmt_str != "cog")
    {
        std::cerr << "Warning: unknown write_output_format '"
                  << fmt_str << "' — defaulting to cog.\n";
    }
    return lsdtt::raster_output_format::cog;
}

// -----------------------------------------------------------------------------
// resolve_window_shape
// -----------------------------------------------------------------------------
static lsdtt::window_shape resolve_window_shape(const std::string& s)
{
    if (s == "square")
    {
        return lsdtt::window_shape::square;
    }
    if (s == "circular")
    {
        return lsdtt::window_shape::circular;
    }

    throw std::runtime_error(
        "Unknown depth_window_shape: " + s +
        ". Expected square or circular."
    );
}

// -----------------------------------------------------------------------------
// resolve_curvature_metric
// -----------------------------------------------------------------------------
static lsdtt::polyfit_metric resolve_curvature_metric(const std::string& s)
{
    if (s == "plan_curvature")
    {
        return lsdtt::polyfit_metric::plan_curvature;
    }
    if (s == "tangential_curvature")
    {
        return lsdtt::polyfit_metric::tangential_curvature;
    }
    if (s == "profile_curvature")
    {
        return lsdtt::polyfit_metric::profile_curvature;
    }
    if (s == "mean_curvature")
    {
        return lsdtt::polyfit_metric::mean_curvature;
    }
    if (s == "maximum_curvature")
    {
        return lsdtt::polyfit_metric::maximum_curvature;
    }
    if (s == "minimum_curvature")
    {
        return lsdtt::polyfit_metric::minimum_curvature;
    }

    throw std::runtime_error(
        "Unknown polyfit_curvature_metric: " + s
    );
}

// -----------------------------------------------------------------------------
// resolve_curvature_sign
// -----------------------------------------------------------------------------
static bool resolve_curvature_sign_positive_is_good(const std::string& s)
{
    if (s == "positive") return true;
    if (s == "negative") return false;

    throw std::runtime_error(
        "Unknown polyfit_curvature_sign: " + s +
        ". Expected positive or negative."
    );
}

// -----------------------------------------------------------------------------
// assert_same_meta
// -----------------------------------------------------------------------------
static void assert_same_meta(
    const lsdtt::lsdtt_raster_meta& a,
    const lsdtt::lsdtt_raster_meta& b,
    const std::string&              a_name,
    const std::string&              b_name)
{
    if (a.n_rows != b.n_rows ||
        a.n_cols != b.n_cols ||
        a.cell_size != b.cell_size ||
        a.x_min != b.x_min ||
        a.y_min != b.y_min ||
        a.no_data_value != b.no_data_value)
    {
        throw std::runtime_error(
            "lsdtt-feature-finder: incompatible rasters: "
            + a_name + " and " + b_name
        );
    }
}

// -----------------------------------------------------------------------------
// raster_cell_is_valid
// -----------------------------------------------------------------------------
static bool raster_cell_is_valid(
    const lsdtt::lsdtt_raster_cpu& raster,
    std::size_t                    idx)
{
    const lsdtt::lsdtt_raster_meta& meta = raster.get_meta();
    const std::vector<float>&       data = raster.get_data();
    const std::vector<uint8_t>&     mask = raster.get_mask();

    return mask[idx] != 0 &&
           !meta.is_nodata(data[idx]) &&
           std::isfinite(data[idx]);
}

// -----------------------------------------------------------------------------
// make_raster
// -----------------------------------------------------------------------------
static lsdtt::lsdtt_raster_cpu make_raster(
    std::vector<float>&&   data,
    std::vector<uint8_t>&& mask,
    const lsdtt::lsdtt_raster_meta& meta)
{
    return lsdtt::lsdtt_raster_cpu(
        std::move(data),
        std::move(mask),
        meta
    );
}

// -----------------------------------------------------------------------------
// write_raster
// -----------------------------------------------------------------------------
static void write_raster(
    lsdtt::lsdtt_data_loader&        loader,
    const lsdtt::lsdtt_raster_cpu&   raster,
    const std::string&               path,
    lsdtt::raster_output_format      fmt)
{
    loader.write_full(
        raster.get_data(),
        raster.get_meta(),
        path,
        fmt
    );
}

// -----------------------------------------------------------------------------
// curvature_to_positive_signal
//
// Converts a curvature raster to a positive-only signal:
//   positive-is-good:  signal = max(curvature, 0)
//   negative-is-good:  signal = max(-curvature, 0)
// -----------------------------------------------------------------------------
static lsdtt::lsdtt_raster_cpu curvature_to_positive_signal(
    const lsdtt::lsdtt_raster_cpu& curvature,
    bool                           positive_is_good)
{
    const lsdtt::lsdtt_raster_meta& meta = curvature.get_meta();
    const std::vector<float>&       in   = curvature.get_data();
    const std::vector<uint8_t>&     mask = curvature.get_mask();

    const std::size_t npix =
        static_cast<std::size_t>(meta.n_rows) *
        static_cast<std::size_t>(meta.n_cols);

    const float nodata = static_cast<float>(meta.no_data_value);

    std::vector<float>   out_data(npix, nodata);
    std::vector<uint8_t> out_mask(npix, uint8_t{0});

    for (std::size_t i = 0; i < npix; ++i)
    {
        if (mask[i] == 0) continue;
        if (meta.is_nodata(in[i])) continue;
        if (!std::isfinite(in[i])) continue;

        double v = positive_is_good
                 ? static_cast<double>(in[i])
                 : -static_cast<double>(in[i]);

        if (v < 0.0) v = 0.0;

        out_data[i] = static_cast<float>(v);
        out_mask[i] = uint8_t{1};
    }

    return make_raster(std::move(out_data), std::move(out_mask), meta);
}

// -----------------------------------------------------------------------------
// clamp_negative_to_zero
//
// Generic helper for signals that should remain non-negative.
// -----------------------------------------------------------------------------
static lsdtt::lsdtt_raster_cpu clamp_negative_to_zero(
    const lsdtt::lsdtt_raster_cpu& raster)
{
    const lsdtt::lsdtt_raster_meta& meta = raster.get_meta();
    const std::vector<float>&       in   = raster.get_data();
    const std::vector<uint8_t>&     mask = raster.get_mask();

    const std::size_t npix =
        static_cast<std::size_t>(meta.n_rows) *
        static_cast<std::size_t>(meta.n_cols);

    const float nodata = static_cast<float>(meta.no_data_value);

    std::vector<float>   out_data(npix, nodata);
    std::vector<uint8_t> out_mask(npix, uint8_t{0});

    for (std::size_t i = 0; i < npix; ++i)
    {
        if (mask[i] == 0) continue;
        if (meta.is_nodata(in[i])) continue;
        if (!std::isfinite(in[i])) continue;

        double v = static_cast<double>(in[i]);
        if (v < 0.0) v = 0.0;

        out_data[i] = static_cast<float>(v);
        out_mask[i] = uint8_t{1};
    }

    return make_raster(std::move(out_data), std::move(out_mask), meta);
}

// -----------------------------------------------------------------------------
// normalise_unit_interval
//
// Min-max normalisation over valid pixels only.
// If the raster is constant over valid pixels, all valid outputs become 0.
// -----------------------------------------------------------------------------
static lsdtt::lsdtt_raster_cpu normalise_unit_interval(
    const lsdtt::lsdtt_raster_cpu& raster)
{
    const lsdtt::lsdtt_raster_meta& meta = raster.get_meta();
    const std::vector<float>&       in   = raster.get_data();
    const std::vector<uint8_t>&     mask = raster.get_mask();

    const std::size_t npix =
        static_cast<std::size_t>(meta.n_rows) *
        static_cast<std::size_t>(meta.n_cols);

    const float nodata = static_cast<float>(meta.no_data_value);

    double vmin =  std::numeric_limits<double>::infinity();
    double vmax = -std::numeric_limits<double>::infinity();
    bool   any_valid = false;

    for (std::size_t i = 0; i < npix; ++i)
    {
        if (mask[i] == 0) continue;
        if (meta.is_nodata(in[i])) continue;
        if (!std::isfinite(in[i])) continue;

        const double v = static_cast<double>(in[i]);
        vmin = std::min(vmin, v);
        vmax = std::max(vmax, v);
        any_valid = true;
    }

    std::vector<float>   out_data(npix, nodata);
    std::vector<uint8_t> out_mask(npix, uint8_t{0});

    if (!any_valid)
    {
        return make_raster(std::move(out_data), std::move(out_mask), meta);
    }

    const double range = vmax - vmin;

    for (std::size_t i = 0; i < npix; ++i)
    {
        if (mask[i] == 0) continue;
        if (meta.is_nodata(in[i])) continue;
        if (!std::isfinite(in[i])) continue;

        double out = 0.0;
        if (range > 0.0)
        {
            out = (static_cast<double>(in[i]) - vmin) / range;
        }

        if (out < 0.0) out = 0.0;
        if (out > 1.0) out = 1.0;

        out_data[i] = static_cast<float>(out);
        out_mask[i] = uint8_t{1};
    }

    return make_raster(std::move(out_data), std::move(out_mask), meta);
}

// -----------------------------------------------------------------------------
// pixelwise_min
// -----------------------------------------------------------------------------
static lsdtt::lsdtt_raster_cpu pixelwise_min(
    const lsdtt::lsdtt_raster_cpu& a,
    const lsdtt::lsdtt_raster_cpu& b)
{
    const lsdtt::lsdtt_raster_meta& meta = a.get_meta();
    assert_same_meta(meta, b.get_meta(), "a", "b");

    const std::vector<float>&   a_data = a.get_data();
    const std::vector<float>&   b_data = b.get_data();
    const std::vector<uint8_t>& a_mask = a.get_mask();
    const std::vector<uint8_t>& b_mask = b.get_mask();

    const std::size_t npix =
        static_cast<std::size_t>(meta.n_rows) *
        static_cast<std::size_t>(meta.n_cols);

    const float nodata = static_cast<float>(meta.no_data_value);

    std::vector<float>   out_data(npix, nodata);
    std::vector<uint8_t> out_mask(npix, uint8_t{0});

    for (std::size_t i = 0; i < npix; ++i)
    {
        const bool valid =
            a_mask[i] != 0 && b_mask[i] != 0 &&
            !meta.is_nodata(a_data[i]) && !meta.is_nodata(b_data[i]) &&
            std::isfinite(a_data[i]) && std::isfinite(b_data[i]);

        if (!valid) continue;

        out_data[i] = std::min(a_data[i], b_data[i]);
        out_mask[i] = uint8_t{1};
    }

    return make_raster(std::move(out_data), std::move(out_mask), meta);
}

// -----------------------------------------------------------------------------
// positive_difference
//
// Returns max(a - b, 0).
// Useful here for local depth:
//   depth = moving_window_mean(dem) - dem
// -----------------------------------------------------------------------------
static lsdtt::lsdtt_raster_cpu positive_difference(
    const lsdtt::lsdtt_raster_cpu& a,
    const lsdtt::lsdtt_raster_cpu& b)
{
    const lsdtt::lsdtt_raster_meta& meta = a.get_meta();
    assert_same_meta(meta, b.get_meta(), "a", "b");

    const std::vector<float>&   a_data = a.get_data();
    const std::vector<float>&   b_data = b.get_data();
    const std::vector<uint8_t>& a_mask = a.get_mask();
    const std::vector<uint8_t>& b_mask = b.get_mask();

    const std::size_t npix =
        static_cast<std::size_t>(meta.n_rows) *
        static_cast<std::size_t>(meta.n_cols);

    const float nodata = static_cast<float>(meta.no_data_value);

    std::vector<float>   out_data(npix, nodata);
    std::vector<uint8_t> out_mask(npix, uint8_t{0});

    for (std::size_t i = 0; i < npix; ++i)
    {
        const bool valid =
            a_mask[i] != 0 && b_mask[i] != 0 &&
            !meta.is_nodata(a_data[i]) && !meta.is_nodata(b_data[i]) &&
            std::isfinite(a_data[i]) && std::isfinite(b_data[i]);

        if (!valid) continue;

        double v = static_cast<double>(a_data[i]) - static_cast<double>(b_data[i]);
        if (v < 0.0) v = 0.0;

        out_data[i] = static_cast<float>(v);
        out_mask[i] = uint8_t{1};
    }

    return make_raster(std::move(out_data), std::move(out_mask), meta);
}

// -----------------------------------------------------------------------------
// weighted_score
//
// Forms a weighted mean of 0..1 signals. Any signal with zero weight is ignored.
// -----------------------------------------------------------------------------
static lsdtt::lsdtt_raster_cpu weighted_score(
    const lsdtt::lsdtt_raster_cpu& persistent,
    const lsdtt::lsdtt_raster_cpu& depth,
    const lsdtt::lsdtt_raster_cpu& slope,
    double                         weight_persistent,
    double                         weight_depth,
    double                         weight_slope)
{
    const lsdtt::lsdtt_raster_meta& meta = persistent.get_meta();

    assert_same_meta(meta, depth.get_meta(),      "persistent", "depth");
    assert_same_meta(meta, slope.get_meta(),      "persistent", "slope");

    if (weight_persistent < 0.0 || weight_depth < 0.0 || weight_slope < 0.0)
    {
        throw std::runtime_error(
            "lsdtt-feature-finder: score weights must be >= 0"
        );
    }

    const double weight_sum = weight_persistent + weight_depth + weight_slope;
    if (weight_sum <= 0.0)
    {
        throw std::runtime_error(
            "lsdtt-feature-finder: at least one score weight must be > 0"
        );
    }

    const std::vector<float>&   p_data = persistent.get_data();
    const std::vector<float>&   d_data = depth.get_data();
    const std::vector<float>&   s_data = slope.get_data();

    const std::vector<uint8_t>& p_mask = persistent.get_mask();
    const std::vector<uint8_t>& d_mask = depth.get_mask();
    const std::vector<uint8_t>& s_mask = slope.get_mask();

    const std::size_t npix =
        static_cast<std::size_t>(meta.n_rows) *
        static_cast<std::size_t>(meta.n_cols);

    const float nodata = static_cast<float>(meta.no_data_value);

    std::vector<float>   out_data(npix, nodata);
    std::vector<uint8_t> out_mask(npix, uint8_t{0});

    for (std::size_t i = 0; i < npix; ++i)
    {
        bool valid = true;

        if (weight_persistent > 0.0)
        {
            valid = valid &&
                    p_mask[i] != 0 &&
                    !meta.is_nodata(p_data[i]) &&
                    std::isfinite(p_data[i]);
        }

        if (weight_depth > 0.0)
        {
            valid = valid &&
                    d_mask[i] != 0 &&
                    !meta.is_nodata(d_data[i]) &&
                    std::isfinite(d_data[i]);
        }

        if (weight_slope > 0.0)
        {
            valid = valid &&
                    s_mask[i] != 0 &&
                    !meta.is_nodata(s_data[i]) &&
                    std::isfinite(s_data[i]);
        }

        if (!valid) continue;

        double sum = 0.0;
        double used = 0.0;

        if (weight_persistent > 0.0)
        {
            sum  += weight_persistent * static_cast<double>(p_data[i]);
            used += weight_persistent;
        }

        if (weight_depth > 0.0)
        {
            sum  += weight_depth * static_cast<double>(d_data[i]);
            used += weight_depth;
        }

        if (weight_slope > 0.0)
        {
            sum  += weight_slope * static_cast<double>(s_data[i]);
            used += weight_slope;
        }

        if (used <= 0.0) continue;

        double v = sum / used;
        if (v < 0.0) v = 0.0;
        if (v > 1.0) v = 1.0;

        out_data[i] = static_cast<float>(v);
        out_mask[i] = uint8_t{1};
    }

    return make_raster(std::move(out_data), std::move(out_mask), meta);
}

// -----------------------------------------------------------------------------
// threshold_binary
// -----------------------------------------------------------------------------
static lsdtt::lsdtt_raster_cpu threshold_binary(
    const lsdtt::lsdtt_raster_cpu& raster,
    double                         threshold)
{
    const lsdtt::lsdtt_raster_meta& meta = raster.get_meta();
    const std::vector<float>&       in   = raster.get_data();
    const std::vector<uint8_t>&     mask = raster.get_mask();

    const std::size_t npix =
        static_cast<std::size_t>(meta.n_rows) *
        static_cast<std::size_t>(meta.n_cols);

    std::vector<float>   out_data(npix, 0.0f);
    std::vector<uint8_t> out_mask(npix, uint8_t{0});

    for (std::size_t i = 0; i < npix; ++i)
    {
        if (mask[i] == 0) continue;
        if (meta.is_nodata(in[i])) continue;
        if (!std::isfinite(in[i])) continue;

        out_data[i] = (static_cast<double>(in[i]) >= threshold) ? 1.0f : 0.0f;
        out_mask[i] = uint8_t{1};
    }

    return make_raster(std::move(out_data), std::move(out_mask), meta);
}

// =============================================================================
// main
// =============================================================================
int main(int argc, char* argv[])
{
    print_banner();

    auto [path, fname] = parse_args(argc, argv);

    // -------------------------------------------------------------------------
    // Parameter file
    // -------------------------------------------------------------------------
    lsdtt::lsdtt_param_parser pp;
    try
    {
        pp = lsdtt::lsdtt_param_parser(path, fname);
    }
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

    // --- General ---
    string_default_map["feature_type"] = "channel_bottoms";
    help_map["feature_type"] = {"string","channel_bottoms",
        "Feature mode to run.",
        "Currently only channel_bottoms is implemented in this first draft."};

    int_default_map["n_threads"] = 0;
    help_map["n_threads"] = {"int","0",
        "Number of OpenMP threads.",
        "0 = use all hardware threads. Applies to wavelet fitting and "
        "polyfit metric computation."};

    string_default_map["write_output_format"] = "cog";
    help_map["write_output_format"] = {"string","cog",
        "Output raster format.",
        "cog: Cloud Optimised GeoTIFF (default). "
        "gtiff: plain compressed GeoTIFF. "
        "same_as_input: match the input raster format."};

    // --- Channel-bottom detector ---
    float_default_map["polyfit_lambda_small"] = 4.0f;
    help_map["polyfit_lambda_small"] = {"float","4.0",
        "Small wavelength for the local quadratic fit, in metres.",
        "Use a scale near the width of the narrow channel bottom."};

    float_default_map["polyfit_lambda_large"] = 10.0f;
    help_map["polyfit_lambda_large"] = {"float","10.0",
        "Larger wavelength for the local quadratic fit, in metres.",
        "Use a broader scale that captures the surrounding hollow or valley form."};

    string_default_map["polyfit_curvature_metric"] = "plan_curvature";
    help_map["polyfit_curvature_metric"] = {"string","plan_curvature",
        "Polyfit metric used as the main channelness curvature signal.",
        "Suggested starting values are plan_curvature or tangential_curvature."};

    string_default_map["polyfit_curvature_sign"] = "positive";
    help_map["polyfit_curvature_sign"] = {"string","positive",
        "Whether channel bottoms are expressed by positive or negative curvature.",
        "Allowed values: positive, negative. This is configurable because sign "
        "conventions can vary between formulations and interpretations."};

    float_default_map["polyfit_z_factor"] = 1.0f;
    help_map["polyfit_z_factor"] = {"float","1.0",
        "Vertical scaling applied before interpreting polyfit derivatives.",
        "Useful if you want to amplify or suppress relief in the polyfit metrics."};

    float_default_map["polyfit_flat_slope_tolerance"] = 1e-12f;
    help_map["polyfit_flat_slope_tolerance"] = {"float","1e-12",
        "Flat-slope tolerance passed to the polyfit metric engine.",
        "Very small default so only nearly flat terrain is treated as directionless."};

    int_default_map["depth_window_radius"] = 12;
    help_map["depth_window_radius"] = {"int","12",
        "Radius of the moving window used to estimate local depth, in pixels.",
        "The depth signal is moving_window_mean(dem) - dem, clamped at zero."};

    string_default_map["depth_window_shape"] = "circular";
    help_map["depth_window_shape"] = {"string","circular",
        "Shape of the moving window used for local depth.",
        "Allowed values: square, circular."};

    float_default_map["channel_weight_persistent"] = 0.7f;
    help_map["channel_weight_persistent"] = {"float","0.7",
        "Weight for the multiscale curvature-persistence signal.",
        "This is usually the dominant term in the first-draft detector."};

    float_default_map["channel_weight_depth"] = 0.3f;
    help_map["channel_weight_depth"] = {"float","0.3",
        "Weight for the local-depth signal.",
        "Helps favour locally low, channel-bottom-like pixels."};

    float_default_map["channel_weight_slope"] = 0.0f;
    help_map["channel_weight_slope"] = {"float","0.0",
        "Optional weight for the slope signal.",
        "Usually kept small. Set to zero to ignore slope entirely."};

    float_default_map["channel_score_threshold"] = 0.5f;
    help_map["channel_score_threshold"] = {"float","0.5",
        "Threshold applied to the final channel score to produce a binary mask.",
        "Lower values produce a larger mask; higher values are more conservative."};

    bool_default_map["write_channel_score_raster"] = true;
    help_map["write_channel_score_raster"] = {"bool","true",
        "Write the final channel-bottom score raster.",
        "The score is a weighted combination of multiscale curvature persistence, "
        "local depth, and optional slope."};

    bool_default_map["write_channel_bottom_mask_raster"] = true;
    help_map["write_channel_bottom_mask_raster"] = {"bool","true",
        "Write the thresholded binary channel-bottom mask.",
        "Output pixels are 1.0 for likely channel bottoms and 0.0 otherwise."};

    bool_default_map["write_channel_diagnostics"] = true;
    help_map["write_channel_diagnostics"] = {"bool","true",
        "Write intermediate diagnostic rasters.",
        "Includes raw curvature rasters, derived channelness signals, curvature "
        "persistence, depth signal, and slope signal."};

    // -------------------------------------------------------------------------
    // Parse parameters
    // -------------------------------------------------------------------------
    pp.parse_all_parameters(
        float_default_map,
        int_default_map,
        bool_default_map,
        string_default_map
    );

    const std::map<std::string, bool>        bf  = pp.get_bool_parameters();
    const std::map<std::string, float>       ff  = pp.get_float_parameters();
    const std::map<std::string, int>         ifs = pp.get_int_parameters();
    const std::map<std::string, std::string> sf  = pp.get_string_parameters();

    const std::string feature_type = sf.at("feature_type");
    const int         n_threads    = ifs.at("n_threads");

    lsdtt::raster_output_format out_fmt;
    lsdtt::polyfit_metric       curvature_metric;
    lsdtt::window_shape         depth_window_shape;
    bool                        positive_curvature_is_good = true;

    try
    {
        out_fmt = resolve_output_format(sf.at("write_output_format"));
        curvature_metric =
            resolve_curvature_metric(sf.at("polyfit_curvature_metric"));
        depth_window_shape =
            resolve_window_shape(sf.at("depth_window_shape"));
        positive_curvature_is_good =
            resolve_curvature_sign_positive_is_good(sf.at("polyfit_curvature_sign"));
    }
    catch (const std::exception& e)
    {
        std::cerr << "Parameter error: " << e.what() << "\n";
        return 1;
    }

    if (feature_type != "channel_bottoms")
    {
        std::cerr << "Error: unsupported feature_type '" << feature_type
                  << "'. This first draft only supports channel_bottoms.\n";
        return 1;
    }

    if (ff.at("polyfit_lambda_small") <= 0.0f ||
        ff.at("polyfit_lambda_large") <= 0.0f)
    {
        std::cerr << "Error: polyfit wavelengths must be > 0.\n";
        return 1;
    }

    if (ff.at("polyfit_lambda_large") < ff.at("polyfit_lambda_small"))
    {
        std::cerr << "Warning: polyfit_lambda_large is smaller than "
                  << "polyfit_lambda_small. Proceeding anyway, but this is "
                  << "usually not what you want.\n";
    }

    if (ifs.at("depth_window_radius") < 1)
    {
        std::cerr << "Error: depth_window_radius must be >= 1.\n";
        return 1;
    }

    if (ff.at("channel_score_threshold") < 0.0f ||
        ff.at("channel_score_threshold") > 1.0f)
    {
        std::cerr << "Error: channel_score_threshold must be in [0, 1].\n";
        return 1;
    }

    // -------------------------------------------------------------------------
    // Open DEM
    // -------------------------------------------------------------------------
    const std::string input_path = pp.get_full_read_path();
    std::cout << "Output format: " << sf.at("write_output_format") << "\n";
    std::cout << "Feature type : " << feature_type << "\n";
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

    try
    {
        loader.get_meta().check_projected_metres();
    }
    catch (const std::exception& e)
    {
        std::cerr << "\nFATAL: " << e.what() << "\n";
        return 1;
    }

    lsdtt::lsdtt_raster_cpu dem(loader);

    // -------------------------------------------------------------------------
    // Feature mode: channel_bottoms
    // -------------------------------------------------------------------------
    try
    {
        std::cout << "\n--- Channel-bottom detector ---\n";
        std::cout << "Small lambda         : " << ff.at("polyfit_lambda_small") << " m\n";
        std::cout << "Large lambda         : " << ff.at("polyfit_lambda_large") << " m\n";
        std::cout << "Curvature metric     : " << sf.at("polyfit_curvature_metric") << "\n";
        std::cout << "Curvature sign       : " << sf.at("polyfit_curvature_sign") << "\n";
        std::cout << "Depth radius         : " << ifs.at("depth_window_radius") << " px\n";
        std::cout << "Depth window shape   : " << sf.at("depth_window_shape") << "\n";
        std::cout << "Weight persistent    : " << ff.at("channel_weight_persistent") << "\n";
        std::cout << "Weight depth         : " << ff.at("channel_weight_depth") << "\n";
        std::cout << "Weight slope         : " << ff.at("channel_weight_slope") << "\n";
        std::cout << "Score threshold      : " << ff.at("channel_score_threshold") << "\n";

        // ---------------------------------------------------------------------
        // Small-scale polyfit
        // ---------------------------------------------------------------------
        std::cout << "\nComputing small-scale quadratic fit...\n";
        lsdtt::poly_fit_coeffs poly_small =
            lsdtt::wavelet_quadratic_fit(
                dem,
                ff.at("polyfit_lambda_small"),
                n_threads
            );

        lsdtt::polyfit_metric_request request_small;
        request_small.metrics =
        {
            curvature_metric,
            lsdtt::polyfit_metric::slope
        };
        request_small.z_factor =
            static_cast<double>(ff.at("polyfit_z_factor"));
        request_small.n_threads = n_threads;
        request_small.flat_slope_tolerance =
            static_cast<double>(ff.at("polyfit_flat_slope_tolerance"));

        std::cout << "Computing small-scale polyfit metrics...\n";
        lsdtt::polyfit_metric_result metrics_small =
            lsdtt::compute_polyfit_metrics(poly_small, request_small);

        const lsdtt::lsdtt_raster_cpu& curvature_small_raw =
            metrics_small.at(curvature_metric);
        const lsdtt::lsdtt_raster_cpu& slope_small_raw =
            metrics_small.at(lsdtt::polyfit_metric::slope);

        // ---------------------------------------------------------------------
        // Large-scale polyfit
        // ---------------------------------------------------------------------
        std::cout << "Computing large-scale quadratic fit...\n";
        lsdtt::poly_fit_coeffs poly_large =
            lsdtt::wavelet_quadratic_fit(
                dem,
                ff.at("polyfit_lambda_large"),
                n_threads
            );

        lsdtt::polyfit_metric_request request_large;
        request_large.metrics = { curvature_metric };
        request_large.z_factor =
            static_cast<double>(ff.at("polyfit_z_factor"));
        request_large.n_threads = n_threads;
        request_large.flat_slope_tolerance =
            static_cast<double>(ff.at("polyfit_flat_slope_tolerance"));

        std::cout << "Computing large-scale polyfit metrics...\n";
        lsdtt::polyfit_metric_result metrics_large =
            lsdtt::compute_polyfit_metrics(poly_large, request_large);

        const lsdtt::lsdtt_raster_cpu& curvature_large_raw =
            metrics_large.at(curvature_metric);

        // ---------------------------------------------------------------------
        // Convert curvature to positive channel signals
        // ---------------------------------------------------------------------
        std::cout << "Building curvature-derived channelness signals...\n";
        lsdtt::lsdtt_raster_cpu curvature_small_signal =
            curvature_to_positive_signal(
                curvature_small_raw,
                positive_curvature_is_good
            );

        lsdtt::lsdtt_raster_cpu curvature_large_signal =
            curvature_to_positive_signal(
                curvature_large_raw,
                positive_curvature_is_good
            );

        lsdtt::lsdtt_raster_cpu curvature_small_norm =
            normalise_unit_interval(curvature_small_signal);

        lsdtt::lsdtt_raster_cpu curvature_large_norm =
            normalise_unit_interval(curvature_large_signal);

        lsdtt::lsdtt_raster_cpu curvature_persistent =
            pixelwise_min(curvature_small_norm, curvature_large_norm);

        // ---------------------------------------------------------------------
        // Local-depth signal
        // ---------------------------------------------------------------------
        std::cout << "Computing local-depth signal...\n";
        lsdtt::lsdtt_raster_cpu local_mean =
            lsdtt::moving_window_mean(
                dem,
                ifs.at("depth_window_radius"),
                depth_window_shape
            );

        lsdtt::lsdtt_raster_cpu depth_signal_raw =
            positive_difference(local_mean, dem);

        lsdtt::lsdtt_raster_cpu depth_signal =
            normalise_unit_interval(depth_signal_raw);

        // ---------------------------------------------------------------------
        // Optional slope signal
        // ---------------------------------------------------------------------
        lsdtt::lsdtt_raster_cpu slope_signal_raw =
            clamp_negative_to_zero(slope_small_raw);

        lsdtt::lsdtt_raster_cpu slope_signal =
            normalise_unit_interval(slope_signal_raw);

        // ---------------------------------------------------------------------
        // Final weighted score
        // ---------------------------------------------------------------------
        std::cout << "Computing final channel-bottom score...\n";
        lsdtt::lsdtt_raster_cpu channel_score =
            weighted_score(
                curvature_persistent,
                depth_signal,
                slope_signal,
                static_cast<double>(ff.at("channel_weight_persistent")),
                static_cast<double>(ff.at("channel_weight_depth")),
                static_cast<double>(ff.at("channel_weight_slope"))
            );

        lsdtt::lsdtt_raster_cpu channel_mask =
            threshold_binary(
                channel_score,
                static_cast<double>(ff.at("channel_score_threshold"))
            );

        // ---------------------------------------------------------------------
        // Write outputs
        // ---------------------------------------------------------------------
        if (bf.at("write_channel_score_raster"))
        {
            const std::string o = out_path(pp, "_channel_score.tif");
            std::cout << "Writing: " << o << "\n";
            write_raster(loader, channel_score, o, out_fmt);
        }

        if (bf.at("write_channel_bottom_mask_raster"))
        {
            const std::string o = out_path(pp, "_channel_bottom_mask.tif");
            std::cout << "Writing: " << o << "\n";
            write_raster(loader, channel_mask, o, out_fmt);
        }

        if (bf.at("write_channel_diagnostics"))
        {
            const std::string o1 = out_path(pp, "_channel_curvature_small.tif");
            const std::string o2 = out_path(pp, "_channel_curvature_large.tif");
            const std::string o3 = out_path(pp, "_channel_signal_small.tif");
            const std::string o4 = out_path(pp, "_channel_signal_large.tif");
            const std::string o5 = out_path(pp, "_channel_persistent.tif");
            const std::string o6 = out_path(pp, "_channel_depth_signal.tif");
            const std::string o7 = out_path(pp, "_channel_slope_signal.tif");

            std::cout << "Writing: " << o1 << "\n";
            write_raster(loader, curvature_small_raw, o1, out_fmt);

            std::cout << "Writing: " << o2 << "\n";
            write_raster(loader, curvature_large_raw, o2, out_fmt);

            std::cout << "Writing: " << o3 << "\n";
            write_raster(loader, curvature_small_norm, o3, out_fmt);

            std::cout << "Writing: " << o4 << "\n";
            write_raster(loader, curvature_large_norm, o4, out_fmt);

            std::cout << "Writing: " << o5 << "\n";
            write_raster(loader, curvature_persistent, o5, out_fmt);

            std::cout << "Writing: " << o6 << "\n";
            write_raster(loader, depth_signal, o6, out_fmt);

            std::cout << "Writing: " << o7 << "\n";
            write_raster(loader, slope_signal, o7, out_fmt);
        }
    }
    catch (const std::exception& e)
    {
        std::cerr << "\nFATAL: feature finding failed: " << e.what() << "\n";
        return 1;
    }

    // -------------------------------------------------------------------------
    // Parameter echo
    // -------------------------------------------------------------------------
    pp.print_parameters();

    std::cout << "\nlsdtt-feature-finder: done.\n";
    return 0;
}
