// lsdtt_raster_synthetic_tests.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2025 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html
//
// -----------------------------------------------------------------------------
// lsdtt_raster_synthetic_tests
//
// Generates synthetic rasters for testing analysis routines and optionally
// tests wavelet_quadratic_fit() against analytical reference surfaces.
//
// Usage:
//   lsdtt_raster_synthetic_tests <path> <parameter_file>
//   lsdtt_raster_synthetic_tests <parameter_file>
//   lsdtt_raster_synthetic_tests -h
//   lsdtt_raster_synthetic_tests -v
//   lsdtt_raster_synthetic_tests -c
//
// The output uses write_path / write_fname from the parameter file.
// -----------------------------------------------------------------------------

#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "../src/core/lsdtt_param_parser.hpp"
#include "../src/core/lsdtt_data_loader.hpp"
#include "../src/core/lsdtt_common_parameters.hpp"
#include "../src/cpu/lsdtt_raster_spectral_cpu.hpp"
#include "../src/cpu/lsdtt_raster_synthetic_surfaces_cpu.hpp"
#include "../src/core/lsdtt_stat_csv.hpp"
#include "../src/cpu/lsdtt_raster_algebra.hpp"

static const std::string k_version  = "0.1.0";
static const std::string k_citation = "https://github.com/LSDtopotools/lsdtt3";
static const std::string k_tool_name = "lsdtt-raster-synthetic-tests";

// -----------------------------------------------------------------------------
// print_banner
// -----------------------------------------------------------------------------
static void print_banner()
{
    std::cout
        << "=========================================================\n"
        << "|| Welcome to lsdtt-raster-synthetic-tests             ||\n"
        << "|| Synthetic raster generation and testing             ||\n"
        << "|| Part of lsdtt3 by Simon M. Mudd                     ||\n"
        << "=========================================================\n"
        << "|| CLI version: " << k_version  << "\n"
        << "|| Citation   : " << k_citation << "\n"
        << "|| -v or --version : full lsdtt3 release version\n"
        << "|| -h              : help explaining the parameters\n"
        << "=========================================================\n\n";
}

// -----------------------------------------------------------------------------
// coeff_name / coeff_raster
// -----------------------------------------------------------------------------
static const char* coeff_name(int i)
{
    switch (i)
    {
        case 0: return "a";
        case 1: return "b";
        case 2: return "c";
        case 3: return "d";
        case 4: return "e";
        case 5: return "f";
        default: throw std::runtime_error("Invalid coefficient index");
    }
}

static const lsdtt::lsdtt_raster_cpu& coeff_raster(
    const lsdtt::poly_fit_coeffs& coeffs,
    int i)
{
    switch (i)
    {
        case 0: return coeffs.a;
        case 1: return coeffs.b;
        case 2: return coeffs.c;
        case 3: return coeffs.d;
        case 4: return coeffs.e;
        case 5: return coeffs.f;
        default: throw std::runtime_error("Invalid coefficient index");
    }
}

// -----------------------------------------------------------------------------
// write_coeff_set
// -----------------------------------------------------------------------------
static void write_coeff_set(
    lsdtt::lsdtt_data_loader&     loader,
    const lsdtt::poly_fit_coeffs& coeffs,
    const std::string&            prefix,
    lsdtt::raster_output_format   out_fmt)
{
    for (int k = 0; k < 6; ++k)
    {
        const std::string out = prefix + "_" + coeff_name(k) + ".tif";
        loader.write_full(coeff_raster(coeffs, k).get_data(), coeff_raster(coeffs, k).get_meta(), out, out_fmt);
    }
}

// -----------------------------------------------------------------------------
// write_periodic_components_csv
// -----------------------------------------------------------------------------
static void write_periodic_components_csv(
    lsdtt::lsdtt_data_loader&                              loader,
    const std::vector<lsdtt::lsdtt_periodic_wave_component>& comps,
    const std::string&                                     out_csv)
{
    std::vector<std::string> headers =
    {
        "mode_x", "mode_y", "amplitude", "phase"
    };

    std::vector<std::vector<std::string>> rows;
    rows.reserve(comps.size());

    for (const auto& c : comps)
    {
        rows.push_back(
        {
            std::to_string(c.mode_x),
            std::to_string(c.mode_y),
            lsdtt::csv_number(c.amplitude),
            lsdtt::csv_number(c.phase)
        });
    }

    loader.write_csv(out_csv, headers, rows);
}

// -----------------------------------------------------------------------------
// write_global_quadratic_coeffs_csv
// -----------------------------------------------------------------------------
static void write_global_quadratic_coeffs_csv(
    lsdtt::lsdtt_data_loader&                           loader,
    const lsdtt::lsdtt_global_quadratic_coefficients&   coeffs,
    const std::string&                                  out_csv)
{
    std::vector<std::string> headers = { "coefficient", "value" };

    std::vector<std::vector<std::string>> rows =
    {
        {"a", lsdtt::csv_number(coeffs.a)},
        {"b", lsdtt::csv_number(coeffs.b)},
        {"c", lsdtt::csv_number(coeffs.c)},
        {"d", lsdtt::csv_number(coeffs.d)},
        {"e", lsdtt::csv_number(coeffs.e)},
        {"f", lsdtt::csv_number(coeffs.f)}
    };

    loader.write_csv(out_csv, headers, rows);
}

// Error statistics against the analytic reference come from the library:
// lsdtt::algebra::raster_difference_stats (src/cpu/lsdtt_raster_algebra.hpp).
// There was never anything synthetic-specific about masked-intersection
// RMSE/MAE/bias, and the benchmarking work wants the same function.
// -----------------------------------------------------------------------------
// run_wavelet_test
//
// Runs wavelet_quadratic_fit() on dem, computes error stats against analytic
// reference, optionally writes fitted and residual rasters, and appends a
// summary row per coefficient to wavelet_summary_rows.
//
// is_exact_oracle:
//   true  - the analytic surface is exact for wavelet_quadratic_fit away from
//            invalid-support edges (global quadratic surface). Prints a warning
//            if any coefficient RMSE exceeds a tight threshold.
//   false - approximate comparison only (periodic surface). No threshold check.
// -----------------------------------------------------------------------------
static void run_wavelet_test(
    lsdtt::lsdtt_data_loader&                 loader,
    const std::string&                        surface_name,
    bool                                      is_exact_oracle,
    const lsdtt::lsdtt_raster_cpu&            dem,
    const lsdtt::poly_fit_coeffs&             analytic,
    double                                    lambda,
    int                                       n_threads,
    bool                                      write_fitted,
    bool                                      write_residuals,
    const std::string&                        write_prefix,
    lsdtt::raster_output_format               out_fmt,
    std::vector<std::vector<std::string>>&    wavelet_summary_rows)
{
    std::cout << "Running wavelet_quadratic_fit on " << surface_name
              << "  lambda=" << lambda << "m\n";

    lsdtt::poly_fit_coeffs fitted = lsdtt::wavelet_quadratic_fit(
        dem, lambda, n_threads
    );

    for (int k = 0; k < 6; ++k)
    {
        const lsdtt::lsdtt_raster_cpu& fit_coeff = coeff_raster(fitted,   k);
        const lsdtt::lsdtt_raster_cpu& ana_coeff = coeff_raster(analytic, k);

        lsdtt::algebra::raster_difference_result stats = lsdtt::algebra::raster_difference_stats(fit_coeff, ana_coeff);

        const char* cn = coeff_name(k);

        std::cout << "  coeff " << cn
                  << "  n=" << stats.count
                  << "  RMSE="    << lsdtt::csv_number(stats.rmse)
                  << "  MAE="     << lsdtt::csv_number(stats.mae)
                  << "  maxAbs="  << lsdtt::csv_number(stats.max_abs)
                  << "  bias="    << lsdtt::csv_number(stats.bias)
                  << "\n";

        if (is_exact_oracle && stats.count > 0 && stats.rmse > 1e-3)
        {
            std::cerr << "  WARNING: coeff " << cn
                      << " RMSE=" << lsdtt::csv_number(stats.rmse)
                      << " exceeds 1e-3 on exact oracle surface.\n";
        }

        wavelet_summary_rows.push_back(
        {
            surface_name,
            is_exact_oracle ? "true" : "false",
            lsdtt::csv_number(lambda),
            std::string(cn),
            std::to_string(stats.count),
            lsdtt::csv_number(stats.rmse),
            lsdtt::csv_number(stats.mae),
            lsdtt::csv_number(stats.max_abs),
            lsdtt::csv_number(stats.bias)
        });

        if (write_fitted)
        {
            const std::string path = write_prefix
                + "_" + surface_name
                + "_wavelet_fitted_" + cn + ".tif";
            loader.write_full(fit_coeff.get_data(), fit_coeff.get_meta(), path, out_fmt);
        }

        if (write_residuals)
        {
            // Compute residual raster: fitted - analytic
            const lsdtt::lsdtt_raster_meta& meta = fit_coeff.get_meta();
            const std::size_t n = static_cast<std::size_t>(meta.n_rows)
                                * static_cast<std::size_t>(meta.n_cols);

            std::vector<float>   res_data(n, static_cast<float>(meta.no_data_value));
            std::vector<uint8_t> res_mask(n, uint8_t{0});

            const std::vector<float>&   fd = fit_coeff.get_data();
            const std::vector<uint8_t>& fm = fit_coeff.get_mask();
            const std::vector<float>&   ad = ana_coeff.get_data();
            const std::vector<uint8_t>& am = ana_coeff.get_mask();

            for (std::size_t i = 0; i < n; ++i)
            {
                if (fm[i] == 0 || am[i] == 0) continue;
                if (meta.is_nodata(fd[i]) || meta.is_nodata(ad[i])) continue;

                res_data[i] = fd[i] - ad[i];
                res_mask[i] = uint8_t{1};
            }

            lsdtt::lsdtt_raster_cpu res_raster(
                std::move(res_data),
                std::move(res_mask),
                meta
            );

            const std::string path = write_prefix
                + "_" + surface_name
                + "_wavelet_residual_" + cn + ".tif";
            loader.write_full(res_raster.get_data(), res_raster.get_meta(), path, out_fmt);
        }
    }
}

// -----------------------------------------------------------------------------
// main
// -----------------------------------------------------------------------------
int main(int argc, char* argv[])
{
    print_banner();

    lsdtt::lsdtt_cli_args cli_args = lsdtt::lsdtt_param_parser::parse_args_deferred(
        argc, argv, k_tool_name, k_version, k_citation);

    
    // -------------------------------------------------------------------------
    // Default parameter maps
    // -------------------------------------------------------------------------
    std::map<std::string, bool> bool_defaults =
    {
        // Stochastic surfaces
        {"write_uniform_random_surface",            false},
        {"write_red_noise_surface",                 false},
        {"write_blue_noise_surface",                false},
        {"write_diamond_square_surface",            false},

        // Periodic wave surface
        {"write_periodic_wave_surface",             false},
        {"periodic_random_sign",                    true},

        // Global quadratic surface
        {"write_global_quadratic_surface",          false},
        {"quadratic_random_sign",                   true},

        // Terrace staircase and alluvial fan (fluvial-surface test domains)
        {"write_terrace_staircase_surface",         false},
        {"write_terrace_staircase_tread_id_raster",        true},
        {"write_terrace_staircase_channel_mask",    true},
        {"write_alluvial_fan_surface",              false},
        {"write_alluvial_fan_lobe_id_raster",              true},

        // Hillslope strip (transient-solver test domain)
        {"write_hillslope_strip_surface",           false},
        {"write_hillslope_strip_channel_mask",      true},
        {"strip_cross_axis_is_x",                   true},
        {"strip_nonlinear_profile",                 false},

        // Wavelet tests
        {"test_wavelet_on_periodic_wave",           false},
        {"test_wavelet_on_global_quadratic",        false},

        // Output options
        {"write_true_coefficient_rasters",          false},
        {"write_wavelet_fitted_coefficient_rasters",false},
        {"write_wavelet_residual_rasters",          false},
    };

    std::map<std::string, float> float_defaults =
    {
        {"cell_size_metres",            1.0f},

        // Uniform noise
        {"uniform_min_value",          -1.0f},
        {"uniform_max_value",           1.0f},

        // Red noise
        {"red_noise_mean",              0.0f},
        {"red_noise_stddev",            1.0f},

        // Blue noise
        {"blue_noise_mean",             0.0f},
        {"blue_noise_stddev",           1.0f},

        // Diamond-square
        {"diamond_square_roughness",    0.55f},
        {"diamond_square_initial_scale",1.0f},
        {"diamond_square_mean",         0.0f},
        {"diamond_square_stddev",       1.0f},
        {"diamond_square_parabola_relief", 0.0f},
        {"diamond_square_noise_relief",    0.0f},

        // Periodic wave
        {"periodic_amplitude_min",      1.0f},
        {"periodic_amplitude_max",     10.0f},

        // Quadratic surface
        {"quadratic_abs_max",           1e-3f},
        {"quadratic_linear_abs_max",    1.0f},
        {"quadratic_constant_abs_max",100.0f},

        // Wavelet fit
        {"polyfit_lambda_metres",     100.0f},

        // Terrace staircase
        {"terrace_tread_width_metres",             300.0f},
        {"terrace_riser_height_metres",              8.0f},
        {"terrace_valley_gradient",            0.002f},
        {"terrace_riser_height_downstream_taper", 0.0f},
        {"terrace_tread_cross_gradient",       0.002f},
        {"terrace_riser_gap_width_metres",         120.0f},
        {"terrace_noise_amplitude_metres",           0.0f},
        {"terrace_channel_elevation",         100.0f},
        {"terrace_valley_wall_gradient",        0.5f},

        // Alluvial fan
        {"fan_apex_x_fraction",                 0.5f},
        {"fan_apex_y_fraction",                0.95f},
        {"fan_apex_elevation",                500.0f},
        {"fan_gradient",                       0.03f},
        {"fan_concavity",                       0.0f},
        {"fan_lobe_inset_metres",                    4.0f},
        {"fan_lobe_half_angle_degrees",            25.0f},
        {"fan_noise_amplitude_metres",               0.0f},

        // Hillslope strip
        {"strip_baselevel_elevation",   0.0f},
        {"strip_half_width_metres",     0.0f},   // <=0 = auto from geometry
        {"strip_uplift",               1e-3f},
        {"strip_diffusivity",          1e-2f},
        {"strip_critical_slope",        0.8f},
        {"strip_init_amplitude",       50.0f},
    };

    std::map<std::string, int> int_defaults =
    {
        {"n_rows",                    128},
        {"n_cols",                    128},
        {"random_seed",                33},
        {"n_threads",                   0},

        // Red/blue noise
        {"red_noise_blur_passes",       6},
        {"blue_noise_blur_passes",      6},

        // Periodic wave
        {"periodic_n_components",       8},
        {"periodic_mode_min",           1},
        {"periodic_mode_max",           6},

        // Terrace staircase
        {"terrace_n_terraces",          4},
        {"terrace_n_riser_gaps",        0},

        // Alluvial fan
        {"fan_n_inset_lobes",           2},
    };

    std::map<std::string, std::string> string_defaults =
    {
        {"write_output_format", "cog"},

        // Hillslope strip
        {"strip_channel_config", "both_ends"},   // both_ends | single_end_low | single_end_high
        {"strip_init_profile",   "flat"},        // flat | analytic_parabola | raised_plateau | perturbed | linear_ramp
    };

    // -------------------------------------------------------------------------
    // Help map
    // -------------------------------------------------------------------------
    std::map<std::string, std::vector<std::string>> help_map =
    {
        {"n_rows",
            {"int","128","Number of rows in the synthetic raster.",""}},
        {"n_cols",
            {"int","128","Number of columns in the synthetic raster.",""}},
        {"cell_size_metres",
            {"float","1.0","Cell size in metres.",
             "Note: synthetic rasters have no CRS. wavelet_quadratic_fit() "
              "is called with the supplied cell_size_metres."}},
        {"random_seed",
            {"int","33","Random seed for all stochastic generators.",""}},
        {"n_threads",
            {"int","0","Number of CPU threads (0 = all available).",""}},
        {"write_output_format",
            {"string","cog","Output format: cog, gtiff, same_as_input.",
             "cog falls back to gtiff for rasters below ~262k pixels."}},
        {"write_uniform_random_surface",
            {"bool","false","Write a uniform random noise surface.",""}},
        {"write_red_noise_surface",
            {"bool","false","Write a red (low-frequency) noise surface.",""}},
        {"write_blue_noise_surface",
            {"bool","false","Write a blue (high-frequency) noise surface.",""}},
        {"write_diamond_square_surface",
            {"bool","false","Write a fractal diamond-square terrain surface.",""}},
        {"write_periodic_wave_surface",
            {"bool","false","Write a periodic cosine-mode test surface.",""}},
        {"write_global_quadratic_surface",
            {"bool","false","Write a globally quadratic test surface.",""}},
        {"test_wavelet_on_periodic_wave",
            {"bool","false",
             "Run wavelet_quadratic_fit() on the periodic surface and report errors.",
             "The periodic surface is an approximate (not exact) oracle."}},
        {"test_wavelet_on_global_quadratic",
            {"bool","false",
             "Run wavelet_quadratic_fit() on the global quadratic surface and report errors.",
             "The global quadratic is the exact oracle for wavelet_quadratic_fit()."}},
        {"write_true_coefficient_rasters",
            {"bool","false","Write analytic coefficient rasters alongside the DEM.",""}},
        {"write_wavelet_fitted_coefficient_rasters",
            {"bool","false","Write wavelet-fitted coefficient rasters.",""}},
        {"write_wavelet_residual_rasters",
            {"bool","false",
             "Write residual (fitted - analytic) rasters for each coefficient.",""}},
        {"polyfit_lambda_metres",
            {"float","100","Wavelength lambda in metres for wavelet_quadratic_fit().",
              "Should be >= cell_size_metres. Larger values give smoother fits."}},
        {"periodic_n_components",
            {"int","8","Number of cosine modes in the periodic test surface.",""}},
        {"periodic_mode_min",
            {"int","1","Minimum mode number for the periodic surface.",""}},
        {"periodic_mode_max",
            {"int","6","Maximum mode number for the periodic surface.",""}},
        {"periodic_amplitude_min",
            {"float","1.0","Minimum mode amplitude.",""}},
        {"periodic_amplitude_max",
            {"float","10.0","Maximum mode amplitude.",""}},
        {"periodic_random_sign",
            {"bool","true","Allow random sign on mode amplitudes.",""}},
        {"quadratic_abs_max",
            {"float","1e-3","Max absolute value for quadratic coefficients a, b, c.",""}},
        {"quadratic_linear_abs_max",
            {"float","1.0","Max absolute value for linear coefficients d, e.",""}},
        {"quadratic_constant_abs_max",
            {"float","100.0","Max absolute value for constant coefficient f.",""}},
        {"quadratic_random_sign",
            {"bool","true","Allow random sign on quadratic coefficients.",""}},
        {"diamond_square_roughness",
            {"float","0.55","Roughness decay per level for diamond-square.",
             "0.5 gives Brownian noise; lower = smoother; higher = rougher."}},
        {"diamond_square_initial_scale",
            {"float","1.0","Initial perturbation amplitude for diamond-square.",""}},
        {"diamond_square_parabola_relief",
            {"float","0.0",
             "Peak relief (m) of an N-S parabola superimposed on the "
             "diamond-square surface. 0 disables it.",
             "A raw diamond-square surface has no regional gradient, so it is "
             "full of closed depressions; filling them creates flats, and "
             "routing across a flat gives feathered parallel drainage that is "
             "an artefact of the fill. The parabola is zero at the north and "
             "south boundaries and peaks at the centre row, giving a regional "
             "slope that drains to both boundaries so most depressions never "
             "form. Pair with diamond_square_noise_relief. "
             "CALIBRATION: measured on a 400x400 diamond-square of 537 m "
             "relief, cells lying in depressions went 16.2% (no parabola) -> "
             "12.0% (300 m) -> 2.7% (1075 m), and the MEAN FILL DEPTH, which is "
             "what actually governs feathering, went 10.3 m -> 3.2 m -> 1.3 m. "
             "But 1075 m (2x the fractal relief) over-corrects: the parabola "
             "then dominates and the channels run near-parallel down the "
             "regional slope, losing the dendritic structure the fractal was "
             "there to provide. About 0.5x the fractal relief was the useful "
             "setting - enough regional gradient to organise the drainage, not "
             "enough to impose its direction."}},
        {"diamond_square_noise_relief",
            {"float","0.0",
             "Amplitude (m) of white noise added to the diamond-square "
             "surface. 0 disables it.",
             "Roughens whatever depressions survive the parabola, so that when "
             "they are filled the result still has a gradient to route on "
             "rather than a flat. Applied everywhere, not only inside "
             "depressions. SIZE IT AGAINST THE MEAN FILL DEPTH, not against the "
             "relief: the noise has to be able to break up a filled pool, so it "
             "must be comparable to or larger than the depth the fill adds. In "
             "the calibration above, 5 m of noise cleared the feathering when "
             "the parabola had brought the mean fill to 3.2 m, while 2 m did "
             "not when the mean fill was 10.3 m. Small enough not to disturb "
             "the fractal structure otherwise."}},

        // Terrace staircase (fluvial-surface stratification test domain)
        {"write_terrace_staircase_surface",
            {"bool","false",
             "Write a symmetrical flight of fluvial terraces either side of a "
             "trunk channel, for height-stratification tests.",
             "Terrace j sits j riser heights above the local channel and "
             "occupies the band [(j-1)*tread_width, j*tread_width) from the "
             "channel on both flanks."}},
        {"write_terrace_staircase_tread_id_raster",
            {"bool","true",
             "Write the ground-truth tread id raster (1..n_terraces on treads, "
             "0 on the channel, on risers, and inside a riser gap).",""}},
        {"write_terrace_staircase_channel_mask",
            {"bool","true",
             "Write the trunk channel mask (1 = channel).",""}},
        {"terrace_n_terraces",
            {"int","4","Number of terrace treads on each flank.",""}},
        {"terrace_tread_width_metres",
            {"float","300","Cross-valley width of each tread in metres.",""}},
        {"terrace_riser_height_metres",
            {"float","8",
             "Height of each riser in metres.",
             "Set this BELOW the flat threshold to make every riser pass the "
             "flat test, which reproduces wholesale under-segmentation with no "
             "riser gaps at all."}},
        {"terrace_valley_gradient",
            {"float","0.002","Downstream fall of the channel in m/m.",""}},
        {"terrace_riser_height_downstream_taper",
            {"float","0.0",
             "Fraction of the riser height lost from the upstream to the "
             "downstream end; 0 gives a flight of constant height.",
             "0.667 with an 12 m riser runs 12 m upstream to 4 m downstream, "
             "so a reach-relative ladder can be compared with a global one."}},
        {"terrace_tread_cross_gradient",
            {"float","0.002",
             "Gradient (m/m) at which each tread falls toward the channel.",
             "NOT cosmetic: at 0 the treads drain downstream off the grid "
             "instead of to the channel, and height above channel is then "
             "undefined over most of the DEM. Must be less than "
             "riser_height_m / tread_width_m."}},
        {"terrace_n_riser_gaps",
            {"int","0",
             "Number of road- or spill-like cuts that break riser continuity.",
             "At each gap the stepped cross-section is replaced by a "
             "continuous ramp, so the gap reads as flat and connects treads "
             "that are metres apart in height."}},
        {"terrace_riser_gap_width_metres",
            {"float","120",
             "Along-valley length of each riser gap in metres.",
             "Must comfortably exceed the roughness window of whatever "
             "detector is being tested, or a window centred in the gap still "
             "reaches the stepped rows either side and never reads as flat. "
             "At a 40 m window a 30 m gap returns 4.5 m of local relief; a "
             "120 m gap returns about 1 m."}},
        {"terrace_noise_amplitude_metres",
            {"float","0","Symmetric per-pixel white noise amplitude, metres.",
             "0 gives an exact surface whose tread elevations can be asserted "
             "to the metre."}},
        {"terrace_channel_elevation",
            {"float","100","Channel elevation at the upstream edge, metres.",""}},
        {"terrace_valley_wall_gradient",
            {"float","0.5",
             "Gradient of the valley wall beyond the outermost tread, m/m.",""}},

        // Alluvial fan (fitted-surface reference-field test domain)
        {"write_alluvial_fan_surface",
            {"bool","false",
             "Write a radial alluvial fan with younger lobes inset below the "
             "older surface.",
             "A fan climbs continuously in height above its channel from toe "
             "to apex, so a height-above-channel ladder cuts it into arcs; the "
             "lobes sit a constant depth below the older surface, so a "
             "residual to a fitted fan separates them cleanly."}},
        {"write_alluvial_fan_lobe_id_raster",
            {"bool","true",
             "Write the ground-truth lobe id raster (0 = older surface, "
             "1..n_inset_lobes = inset lobes).",""}},
        {"fan_apex_x_fraction",
            {"float","0.5","Apex easting as a fraction of the grid width.",""}},
        {"fan_apex_y_fraction",
            {"float","0.95",
             "Apex northing as a fraction of the grid height, increasing "
             "southward.",
             "The default puts the apex near the south edge so the fan spreads "
             "north across the grid."}},
        {"fan_apex_elevation",
            {"float","500","Elevation at the apex, metres.",""}},
        {"fan_gradient",
            {"float","0.03","Radial fall at the apex in m/m.",""}},
        {"fan_concavity",
            {"float","0.0",
             "Reduction of the radial gradient with distance; 0 gives a "
             "straight cone.",
             "Must be less than fan_gradient. A positive value gives the "
             "concave-up profile a real fan has, which a degree-2 polynomial "
             "fits and a plane does not."}},
        {"fan_n_inset_lobes",
            {"int","2","Number of inset lobes radiating from the apex.",""}},
        {"fan_lobe_inset_metres",
            {"float","4","Depth of each lobe below the older surface, metres.",""}},
        {"fan_lobe_half_angle_degrees",
            {"float","25","Half-angle of each lobe wedge in degrees.",""}},
        {"fan_noise_amplitude_metres",
            {"float","0","Symmetric per-pixel white noise amplitude, metres.",""}},

        // Hillslope strip (transient Roering hillslope solver test domain)
        {"write_hillslope_strip_surface",
            {"bool","false",
             "Write a 1-D hillslope-strip DEM (and channel mask) for the "
             "transient hillslope solver validation tests.",""}},
        {"write_hillslope_strip_channel_mask",
            {"bool","true",
             "Write the geometry-defined channel mask (1 = channel/Dirichlet "
             "sink, 0 = hillslope) alongside the strip DEM.",""}},
        {"strip_channel_config",
            {"string","both_ends",
             "Channel placement: both_ends (divide at centre), single_end_low, "
             "or single_end_high (no-flux divide at the opposite end).",""}},
        {"strip_init_profile",
            {"string","flat",
             "Initial profile: flat, analytic_parabola, raised_plateau, "
             "perturbed, or linear_ramp.",
             "analytic_parabola is the steady fixed point; raised_plateau / "
             "perturbed seed transient-relaxation tests."}},
        {"strip_cross_axis_is_x",
            {"bool","true",
             "true: elevation varies E-W (channels are N-S end columns); "
             "false: varies N-S (channels are E-W end rows).",""}},
        {"strip_nonlinear_profile",
            {"bool","false",
             "For analytic_parabola: true uses the nonlinear Roering (2007) "
             "steady profile; false uses the linear parabola.",""}},
        {"strip_baselevel_elevation",
            {"float","0.0","Channel (Dirichlet) baselevel elevation in metres.",""}},
        {"strip_half_width_metres",
            {"float","0.0",
             "Hillslope half-width L (divide-to-channel distance, m) for the "
             "analytic profile. <= 0 auto-derives L from the grid (recommended; "
             "makes analytic_parabola a discrete fixed point).",""}},
        {"strip_uplift",
            {"float","1e-3","Uplift rate U (m/yr) used to build analytic_parabola.",""}},
        {"strip_diffusivity",
            {"float","1e-2","Diffusivity D (m^2/yr) used to build analytic_parabola.",""}},
        {"strip_critical_slope",
            {"float","0.8","Critical slope Sc (m/m); used only for the nonlinear profile.",""}},
        {"strip_init_amplitude",
            {"float","50.0",
             "Height scale (m): plateau height (raised_plateau), noise amplitude "
             "(perturbed), or crest height (linear_ramp).",""}},
    };

    // -------------------------------------------------------------------------
    // Help - rendered from the declarations above; no parameter file needed
    // -------------------------------------------------------------------------
    // This CLI declares its help map with brace initialisation, so copy it
    // into a registry for the help renderer (alphabetical order; the
    // entries carry no group metadata yet).
    lsdtt::lsdtt_help_registry help_registry;
    for (const auto& [help_key, help_entry] : help_map)
    {
        help_registry[help_key] = help_entry;
    }

    if (cli_args.help_requested)
    {
        lsdtt::lsdtt_param_parser::print_terminal_help(
            k_tool_name, k_version, k_citation, help_registry, cli_args);
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

    // Parameter renames onto the project conventions (docs/cli_parameter_conventions.md):
    // units are spelled out, and a bool that writes a raster ends in _raster.
    // The old keys keep working and print a warning naming the replacement.
    pp.set_deprecated_aliases({
        {"cell_size",                         "cell_size_metres"},
        {"seed",                              "random_seed"},
        {"strip_half_width",                  "strip_half_width_metres"},
        {"fan_lobe_half_angle_deg",          "fan_lobe_half_angle_degrees"},
        {"fan_lobe_inset_m",                 "fan_lobe_inset_metres"},
        {"fan_noise_amplitude_m",            "fan_noise_amplitude_metres"},
        {"terrace_noise_amplitude_m",        "terrace_noise_amplitude_metres"},
        {"terrace_riser_gap_width_m",        "terrace_riser_gap_width_metres"},
        {"terrace_riser_height_m",           "terrace_riser_height_metres"},
        {"terrace_tread_width_m",            "terrace_tread_width_metres"},
        {"wavelet_lambda",                   "polyfit_lambda_metres"},
        {"write_alluvial_fan_lobe_id",       "write_alluvial_fan_lobe_id_raster"},
        {"write_terrace_staircase_tread_id", "write_terrace_staircase_tread_id_raster"},
    });

    pp.parse_all_parameters(float_defaults, int_defaults, bool_defaults, string_defaults);

    auto bf  = pp.get_bool_parameters();
    auto ff  = pp.get_float_parameters();
    auto ifs = pp.get_int_parameters();
    auto sf  = pp.get_string_parameters();

    // Resolve output format
    lsdtt::raster_output_format out_fmt = lsdtt::raster_output_format::cog;
    {
        const std::string& fmt_str = sf["write_output_format"];

        if (fmt_str == "gtiff")
        {
            out_fmt = lsdtt::raster_output_format::gtiff;
        }
        else if (fmt_str == "same_as_input")
        {
            out_fmt = lsdtt::raster_output_format::same_as_input;
        }
        else if (fmt_str != "cog")
        {
            std::cerr << "Warning: unknown write_output_format '"
                      << fmt_str << "' - defaulting to cog.\n";
        }

        std::cout << "Output format: " << fmt_str << "\n";
    }

    const int   n_rows    = ifs["n_rows"];
    const int   n_cols    = ifs["n_cols"];
    const int   seed      = ifs["random_seed"];
    const int   n_threads = ifs["n_threads"];
    const float cell_size = ff["cell_size_metres"];

    lsdtt::lsdtt_data_loader loader;

    std::vector<std::vector<std::string>> wavelet_summary_rows;
    const std::string summary_csv = pp.make_output_path( "_wavelet_test_summary.csv");

    try
    {
        // ---------------------------------------------------------------------
        // Random / stochastic surfaces
        // ---------------------------------------------------------------------
        if (bf["write_uniform_random_surface"])
        {
            std::cout << "\n--- Uniform random noise ---\n";

            lsdtt::lsdtt_raster_cpu r =
                lsdtt::generate_uniform_random_noise_surface(
                    n_rows, n_cols, cell_size,
                    seed,
                    ff["uniform_min_value"],
                    ff["uniform_max_value"],
                    n_threads
                );

            loader.write_full(r.get_data(), r.get_meta(), pp.make_output_path( "_uniform_random.tif"), out_fmt);
        }

        if (bf["write_red_noise_surface"])
        {
            std::cout << "\n--- Red noise ---\n";

            lsdtt::lsdtt_raster_cpu r =
                lsdtt::generate_red_noise_surface(
                    n_rows, n_cols, cell_size,
                    seed,
                    ff["red_noise_mean"],
                    ff["red_noise_stddev"],
                    ifs["red_noise_blur_passes"],
                    n_threads
                );

            loader.write_full(r.get_data(), r.get_meta(), pp.make_output_path( "_red_noise.tif"), out_fmt);
        }

        if (bf["write_blue_noise_surface"])
        {
            std::cout << "\n--- Blue noise ---\n";

            lsdtt::lsdtt_raster_cpu r =
                lsdtt::generate_blue_noise_surface(
                    n_rows, n_cols, cell_size,
                    seed,
                    ff["blue_noise_mean"],
                    ff["blue_noise_stddev"],
                    ifs["blue_noise_blur_passes"],
                    n_threads
                );

            loader.write_full(r.get_data(), r.get_meta(), pp.make_output_path( "_blue_noise.tif"), out_fmt);
        }

        if (bf["write_diamond_square_surface"])
        {
            std::cout << "\n--- Diamond-square terrain ---\n";

            lsdtt::lsdtt_raster_cpu r =
                lsdtt::generate_diamond_square_surface(
                    n_rows, n_cols, cell_size,
                    seed,
                    ff["diamond_square_roughness"],
                    ff["diamond_square_initial_scale"],
                    ff["diamond_square_mean"],
                    ff["diamond_square_stddev"],
                    n_threads
                );

            // -----------------------------------------------------------------
            // Drainage conditioning: parabolic template, then roughening.
            //
            // A raw diamond-square surface has no regional gradient, so it is
            // full of closed depressions. Pit filling then flattens each one
            // into a level pool, and routing across a flat produces the
            // characteristic "feathered" parallel drainage that is an artefact
            // of the fill, not of the fractal.
            //
            // The fix (after LSDRasterMaker::create_diamond_square_surface in
            // LSDTT_Development) is two-part:
            //
            //   1. superimpose a low-relief N-S parabola, zero at the north and
            //      south boundaries and peaking at the centre row, so the
            //      surface has a regional slope draining to both boundaries and
            //      most depressions never form;
            //   2. add white noise, which roughens whatever depressions remain
            //      so that when they are filled the resulting surface still has
            //      a gradient to route on, instead of a flat.
            //
            // Both default to 0 (off), which reproduces the previous raw
            // behaviour exactly. Note the roughening is applied everywhere, not
            // only inside depressions -- as in the reference -- which is
            // simpler and has the same effect on the pools.
            if (ff["diamond_square_parabola_relief"] > 0.0f)
            {
                std::cout << "  Superimposing N-S parabola, peak "
                          << ff["diamond_square_parabola_relief"] << " m\n";
                r = lsdtt::superimpose_ns_parabola(
                        r, ff["diamond_square_parabola_relief"]);
            }

            if (ff["diamond_square_noise_relief"] > 0.0f)
            {
                std::cout << "  Roughening with white noise, amplitude "
                          << ff["diamond_square_noise_relief"] << " m\n";

                // Uniform, pixelwise independent: the point is to break flats,
                // so a correlated noise would defeat it.
                const lsdtt::lsdtt_raster_cpu noise =
                    lsdtt::generate_uniform_random_noise_surface(
                        n_rows, n_cols, cell_size,
                        seed + 977,          // decorrelated from the DS seed
                        0.0,
                        static_cast<double>(ff["diamond_square_noise_relief"]),
                        n_threads
                    );

                float* z = r.get_data().data();
                const float* nz = noise.get_data().data();
                const std::size_t n_cells =
                    static_cast<std::size_t>(n_rows) *
                    static_cast<std::size_t>(n_cols);
                for (std::size_t i = 0; i < n_cells; ++i) z[i] += nz[i];
            }

            loader.write_full(r.get_data(), r.get_meta(), pp.make_output_path( "_diamond_square.tif"), out_fmt);
        }

        // ---------------------------------------------------------------------
        // Periodic wave surface and optional wavelet test
        // ---------------------------------------------------------------------
        if (bf["write_periodic_wave_surface"] || bf["test_wavelet_on_periodic_wave"])
        {
            std::cout << "\n--- Periodic wave surface ---\n";

            lsdtt::lsdtt_periodic_wave_surface periodic =
                lsdtt::generate_periodic_wave_test_surface(
                    n_rows, n_cols, cell_size,
                    ifs["periodic_n_components"],
                    seed,
                    ff["periodic_amplitude_min"],
                    ff["periodic_amplitude_max"],
                    ifs["periodic_mode_min"],
                    ifs["periodic_mode_max"],
                    bf["periodic_random_sign"],
                    n_threads
                );

            if (bf["write_periodic_wave_surface"])
            {
                loader.write_full(periodic.dem.get_data(), periodic.dem.get_meta(), pp.make_output_path( "_periodic_wave.tif"), out_fmt);

                write_periodic_components_csv(
                    loader,
                    periodic.components,
                    pp.make_output_path( "_periodic_wave_components.csv")
                );

                if (bf["write_true_coefficient_rasters"])
                {
                    write_coeff_set(
                        loader,
                        periodic.analytic_taylor,
                        pp.make_output_path( "_periodic_wave_true"),
                        out_fmt
                    );
                }
            }

            if (bf["test_wavelet_on_periodic_wave"])
            {
                run_wavelet_test(
                    loader,
                    "periodic_wave",
                    false,  // not an exact oracle for finite lambda
                    periodic.dem,
                    periodic.analytic_taylor,
                    ff["polyfit_lambda_metres"],
                    n_threads,
                    bf["write_wavelet_fitted_coefficient_rasters"],
                    bf["write_wavelet_residual_rasters"],
                    pp.get_full_write_path(),
                    out_fmt,
                    wavelet_summary_rows
                );
            }
        }

        // ---------------------------------------------------------------------
        // Global quadratic surface and optional wavelet test
        // ---------------------------------------------------------------------
        if (bf["write_global_quadratic_surface"] || bf["test_wavelet_on_global_quadratic"])
        {
            std::cout << "\n--- Global quadratic surface ---\n";

            lsdtt::lsdtt_global_quadratic_surface quad =
                lsdtt::generate_random_quadratic_test_surface(
                    n_rows, n_cols, cell_size,
                    seed,
                    ff["quadratic_abs_max"],
                    ff["quadratic_linear_abs_max"],
                    ff["quadratic_constant_abs_max"],
                    bf["quadratic_random_sign"],
                    n_threads
                );

            if (bf["write_global_quadratic_surface"])
            {
                loader.write_full(quad.dem.get_data(), quad.dem.get_meta(), pp.make_output_path( "_global_quadratic.tif"), out_fmt);

                write_global_quadratic_coeffs_csv(
                    loader,
                    quad.global_coeffs,
                    pp.make_output_path( "_global_quadratic_coefficients.csv")
                );

                if (bf["write_true_coefficient_rasters"])
                {
                    write_coeff_set(
                        loader,
                        quad.analytic_taylor,
                        pp.make_output_path( "_global_quadratic_true"),
                        out_fmt
                    );
                }
            }

            if (bf["test_wavelet_on_global_quadratic"])
            {
                run_wavelet_test(
                    loader,
                    "global_quadratic",
                    true,   // exact oracle away from invalid support edges
                    quad.dem,
                    quad.analytic_taylor,
                    ff["polyfit_lambda_metres"],
                    n_threads,
                    bf["write_wavelet_fitted_coefficient_rasters"],
                    bf["write_wavelet_residual_rasters"],
                    pp.get_full_write_path(),
                    out_fmt,
                    wavelet_summary_rows
                );
            }
        }

        // ---------------------------------------------------------------------
        // Hillslope strip (transient-solver test domain)
        // ---------------------------------------------------------------------
        // ---------------------------------------------------------------------
        // Terrace staircase
        // ---------------------------------------------------------------------
        if (bf["write_terrace_staircase_surface"])
        {
            lsdtt::lsdtt_terrace_staircase_surface st =
                lsdtt::generate_terrace_staircase_surface(
                    n_rows, n_cols, cell_size,
                    ifs["terrace_n_terraces"],
                    ff["terrace_tread_width_metres"],
                    ff["terrace_riser_height_metres"],
                    ff["terrace_valley_gradient"],
                    ff["terrace_riser_height_downstream_taper"],
                    ff["terrace_tread_cross_gradient"],
                    ifs["terrace_n_riser_gaps"],
                    ff["terrace_riser_gap_width_metres"],
                    ff["terrace_noise_amplitude_metres"],
                    ff["terrace_channel_elevation"],
                    ff["terrace_valley_wall_gradient"],
                    seed,
                    n_threads
                );

            loader.write_full(st.dem.get_data(), st.dem.get_meta(),
                pp.make_output_path("_terrace_staircase.tif"), out_fmt);

            if (bf["write_terrace_staircase_tread_id_raster"])
            {
                loader.write_full(
                    st.tread_id.get_data(), st.tread_id.get_meta(),
                    pp.make_output_path("_terrace_staircase_tread_id.tif"),
                    out_fmt);
            }
            if (bf["write_terrace_staircase_channel_mask"])
            {
                loader.write_full(
                    st.channel_mask.get_data(), st.channel_mask.get_meta(),
                    pp.make_output_path("_terrace_staircase_channel_mask.tif"),
                    out_fmt);
            }
        }

        // ---------------------------------------------------------------------
        // Alluvial fan
        // ---------------------------------------------------------------------
        if (bf["write_alluvial_fan_surface"])
        {
            lsdtt::lsdtt_alluvial_fan_surface fan =
                lsdtt::generate_alluvial_fan_surface(
                    n_rows, n_cols, cell_size,
                    ff["fan_apex_x_fraction"],
                    ff["fan_apex_y_fraction"],
                    ff["fan_apex_elevation"],
                    ff["fan_gradient"],
                    ff["fan_concavity"],
                    ifs["fan_n_inset_lobes"],
                    ff["fan_lobe_inset_metres"],
                    ff["fan_lobe_half_angle_degrees"],
                    ff["fan_noise_amplitude_metres"],
                    seed,
                    n_threads
                );

            loader.write_full(fan.dem.get_data(), fan.dem.get_meta(),
                pp.make_output_path("_alluvial_fan.tif"), out_fmt);

            if (bf["write_alluvial_fan_lobe_id_raster"])
            {
                loader.write_full(
                    fan.lobe_id.get_data(), fan.lobe_id.get_meta(),
                    pp.make_output_path("_alluvial_fan_lobe_id.tif"), out_fmt);
            }
        }

        if (bf["write_hillslope_strip_surface"])
        {
            std::cout << "\n--- Hillslope strip ---\n";

            // Resolve channel placement
            lsdtt::lsdtt_strip_channel_config cfg =
                lsdtt::lsdtt_strip_channel_config::both_ends;
            {
                const std::string& s = sf["strip_channel_config"];
                if (s == "single_end_low")
                    cfg = lsdtt::lsdtt_strip_channel_config::single_end_low;
                else if (s == "single_end_high")
                    cfg = lsdtt::lsdtt_strip_channel_config::single_end_high;
                else if (s != "both_ends")
                    std::cerr << "Warning: unknown strip_channel_config '" << s
                              << "' - defaulting to both_ends.\n";
            }

            // Resolve init profile
            lsdtt::lsdtt_strip_init_profile prof =
                lsdtt::lsdtt_strip_init_profile::flat;
            {
                const std::string& s = sf["strip_init_profile"];
                if (s == "analytic_parabola")
                    prof = lsdtt::lsdtt_strip_init_profile::analytic_parabola;
                else if (s == "raised_plateau")
                    prof = lsdtt::lsdtt_strip_init_profile::raised_plateau;
                else if (s == "perturbed")
                    prof = lsdtt::lsdtt_strip_init_profile::perturbed;
                else if (s == "linear_ramp")
                    prof = lsdtt::lsdtt_strip_init_profile::linear_ramp;
                else if (s != "flat")
                    std::cerr << "Warning: unknown strip_init_profile '" << s
                              << "' - defaulting to flat.\n";
            }

            lsdtt::lsdtt_hillslope_strip_surface strip =
                lsdtt::generate_hillslope_strip_surface(
                    n_rows, n_cols, cell_size,
                    cfg, prof,
                    bf["strip_cross_axis_is_x"],
                    ff["strip_baselevel_elevation"],
                    ff["strip_half_width_metres"],
                    ff["strip_uplift"],
                    ff["strip_diffusivity"],
                    ff["strip_critical_slope"],
                    bf["strip_nonlinear_profile"],
                    ff["strip_init_amplitude"],
                    seed,
                    n_threads
                );

            loader.write_full(strip.dem.get_data(), strip.dem.get_meta(),
                pp.make_output_path("_hillslope_strip.tif"), out_fmt);

            if (bf["write_hillslope_strip_channel_mask"])
            {
                loader.write_full(
                    strip.channel_mask.get_data(), strip.channel_mask.get_meta(),
                    pp.make_output_path("_hillslope_strip_channel_mask.tif"), out_fmt);
            }
        }

        // ---------------------------------------------------------------------
        // Summary CSV for wavelet tests
        // ---------------------------------------------------------------------
        if (!wavelet_summary_rows.empty())
        {
            std::vector<std::string> headers =
            {
                "surface",
                "exact_oracle",
                "lambda",
                "coefficient",
                "n_valid",
                "rmse",
                "mae",
                "max_abs",
                "bias"
            };

            loader.write_csv(summary_csv, headers, wavelet_summary_rows);
            std::cout << "Wavelet summary written to: " << summary_csv << "\n";
        }

        // ---------------------------------------------------------------------
        // Help file + parameter echo
        // ---------------------------------------------------------------------
        pp.print_help(
            help_map,
            pp.get_write_fname() + "_README",
            k_version,
            k_citation
        );

        pp.print_parameters();

        std::cout << "\nlsdtt_raster_synthetic_tests: done.\n";
    }
    catch (const std::exception& e)
    {
        std::cerr << "FATAL: " << e.what() << "\n";
        return 1;
    }

    return 0;
}
