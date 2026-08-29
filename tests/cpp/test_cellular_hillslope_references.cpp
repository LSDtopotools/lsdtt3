// test_cellular_hillslope_references.cpp
//
// Small Newton-Roering reference fixtures for the cellular hillslope plan.
// These are intentionally narrow or small: they freeze converged reference
// fields and coupling behaviour without putting the 250x250 OCN benchmark in
// the normal development loop.

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <utility>
#include <vector>

#include "../../src/cpu/lsdtt_flowinfo_cpu.hpp"
#include "../../src/cpu/lsdtt_hillslope_newton_cpu.hpp"
#include "../../src/cpu/lsdtt_raster_synthetic_surfaces_cpu.hpp"
#include "../../src/cpu/lsdtt_stream_power_cpu.hpp"
#include "../../src/cpu/lsdtt_synthetic_network_cpu.hpp"

namespace
{

using clock_type = std::chrono::steady_clock;

constexpr double PI = 3.14159265358979323846;
constexpr float  DX = 5.0f;
constexpr float  D = 0.1f;
constexpr float  U = 4.0e-4f;
constexpr float  SC = 0.8f;

int g_checks = 0;
int g_failures = 0;
int g_timing_repeats = 1;

void check(bool ok, const std::string& what)
{
    ++g_checks;
    if (ok) std::cout << "  ok  : " << what << "\n";
    else
    {
        ++g_failures;
        std::cout << "  FAIL: " << what << "\n";
    }
}

double elapsed_ms(clock_type::time_point start)
{
    return std::chrono::duration<double, std::milli>(clock_type::now() - start).count();
}

double median(std::vector<double> values)
{
    std::sort(values.begin(), values.end());
    return values[values.size() / 2];
}

lsdtt::lsdtt_raster_meta make_meta(int nr, int nc, double dx = DX)
{
    lsdtt::lsdtt_raster_meta m;
    m.n_rows = nr;
    m.n_cols = nc;
    m.cell_size = dx;
    return m;
}

std::vector<std::uint8_t> channel_mask_from_raster(
    const lsdtt::lsdtt_raster_cpu& mask_raster)
{
    std::vector<std::uint8_t> out(mask_raster.get_data().size(), 0u);
    for (std::size_t i = 0; i < out.size(); ++i)
        out[i] = mask_raster.get_data()[i] != 0.0f ? 1u : 0u;
    return out;
}

std::vector<std::uint8_t> north_south_channels(int nr, int nc)
{
    std::vector<std::uint8_t> out(static_cast<std::size_t>(nr) * nc, 0u);
    for (int c = 0; c < nc; ++c)
    {
        out[static_cast<std::size_t>(c)] = 1u;
        out[static_cast<std::size_t>(nr - 1) * nc + c] = 1u;
    }
    return out;
}

std::vector<std::uint8_t> edge_channels(int nr, int nc)
{
    std::vector<std::uint8_t> out(static_cast<std::size_t>(nr) * nc, 0u);
    for (int r = 0; r < nr; ++r)
        for (int c = 0; c < nc; ++c)
            if (r == 0 || c == 0 || r == nr - 1 || c == nc - 1)
                out[static_cast<std::size_t>(r) * nc + c] = 1u;
    return out;
}

void perturb_strip(lsdtt::lsdtt_raster_cpu& dem,
                   const std::vector<std::uint8_t>& channels,
                   double amplitude)
{
    const int nc = dem.get_meta().n_cols;
    std::vector<float>& z = dem.get_data();
    for (std::size_t i = 0; i < z.size(); ++i)
    {
        if (channels[i]) continue;
        const int c = static_cast<int>(i % static_cast<std::size_t>(nc));
        const double x = static_cast<double>(c) / static_cast<double>(nc - 1);
        z[i] += static_cast<float>(amplitude * std::sin(2.0 * PI * x));
    }
}

struct reference_row
{
    std::string name;
    std::size_t cells = 0;
    double setup_ms = 0.0;
    double fluvial_ms = 0.0;
    double newton_ms = 0.0;
    int iterations = 0;
    int line_search_steps = 0;
    double residual = 0.0;
    double max_slope_ratio = 0.0;
    double elevation_sum = 0.0;
    double weighted_sum = 0.0;
    double fluvial_removed_m3 = 0.0;
    double hillslope_exported_m3 = 0.0;
    double uplift_added_m3 = 0.0;
    bool converged = false;
};

double elevation_sum(const lsdtt::lsdtt_raster_cpu& dem)
{
    double sum = 0.0;
    for (std::size_t i = 0; i < dem.get_data().size(); ++i)
        if (dem.get_mask()[i]) sum += static_cast<double>(dem.get_data()[i]);
    return sum;
}

void fingerprint(const lsdtt::lsdtt_raster_cpu& dem, reference_row& row)
{
    const std::vector<float>& z = dem.get_data();
    for (std::size_t i = 0; i < z.size(); ++i)
    {
        if (!dem.get_mask()[i]) continue;
        row.elevation_sum += static_cast<double>(z[i]);
        row.weighted_sum += static_cast<double>(i + 1) * static_cast<double>(z[i]);
    }
}

lsdtt::lsdtt_hillslope_newton_options newton_options(
    double diffusivity = D,
    double uplift = U,
    const std::string& boundary_conditions = "nnnn")
{
    lsdtt::lsdtt_hillslope_newton_options o;
    o.D = static_cast<float>(diffusivity);
    o.U = static_cast<float>(uplift);
    o.Sc = SC;
    o.nonlinear = true;
    o.boundary_conditions = boundary_conditions;
    o.residual_tol = 1e-10;
    o.max_iterations = 60;
    o.solver_choice = 0;
    o.n_threads = 1;
    return o;
}

reference_row solve_reference(
    const std::string& name,
    const lsdtt::lsdtt_raster_cpu& dem,
    const std::vector<std::uint8_t>& channels,
    const lsdtt::lsdtt_hillslope_newton_options& options,
    const std::vector<float>& old_elevation = {},
    lsdtt::lsdtt_raster_cpu* elevation_out = nullptr)
{
    reference_row row;
    row.name = name;
    row.cells = dem.get_data().size();

    std::vector<double> timings;
    timings.reserve(static_cast<std::size_t>(g_timing_repeats));
    lsdtt::lsdtt_hillslope_newton_result result;
    for (int repeat = 0; repeat < g_timing_repeats; ++repeat)
    {
        const auto start = clock_type::now();
        result = lsdtt::compute_hillslope_newton(
            dem, channels, options, old_elevation);
        timings.push_back(elapsed_ms(start));
    }
    row.newton_ms = median(std::move(timings));
    row.iterations = result.n_iterations;
    row.line_search_steps = result.n_line_search_steps;
    row.residual = result.final_residual;
    row.max_slope_ratio = result.max_slope_ratio;
    row.converged = result.converged && result.had_eigen;
    fingerprint(result.elevation, row);

    check(result.had_eigen, name + " used Eigen/Newton");
    check(result.converged, name + " Newton solve converged: " + result.message);
    check(std::isfinite(row.elevation_sum) && std::isfinite(row.weighted_sum),
          name + " produced finite reference fingerprints");
    if (elevation_out) *elevation_out = result.elevation;
    return row;
}

lsdtt::lsdtt_raster_cpu make_cross_slope_dem(int nr, int nc)
{
    std::vector<float> z(static_cast<std::size_t>(nr) * nc, 0.0f);
    std::vector<std::uint8_t> mask(z.size(), 1u);
    for (int r = 0; r < nr; ++r)
        for (int c = 0; c < nc; ++c)
        {
            const double y = static_cast<double>(r) / static_cast<double>(nr - 1);
            const double x = static_cast<double>(c) / static_cast<double>(nc - 1);
            z[static_cast<std::size_t>(r) * nc + c] = static_cast<float>(
                6.0 * std::sin(PI * y) * (0.8 + 0.2 * x * x)
                + 0.5 * std::sin(2.0 * PI * x));
        }
    return lsdtt::lsdtt_raster_cpu(std::move(z), std::move(mask), make_meta(nr, nc));
}

reference_row solve_periodic_reference()
{
    constexpr int nr = 17;
    constexpr int nc = 33;
    constexpr int shift = 7;
    const lsdtt::lsdtt_raster_cpu dem = make_cross_slope_dem(nr, nc);
    const std::vector<std::uint8_t> channels = north_south_channels(nr, nc);
    const auto options = newton_options(D, U, "npnp");

    reference_row row = solve_reference("periodic_17x33", dem, channels, options);

    std::vector<float> rolled_z(dem.get_data().size(), 0.0f);
    for (int r = 0; r < nr; ++r)
        for (int c = 0; c < nc; ++c)
            rolled_z[static_cast<std::size_t>(r) * nc + ((c + shift) % nc)] =
                dem.get_data()[static_cast<std::size_t>(r) * nc + c];
    lsdtt::lsdtt_raster_cpu rolled(
        std::move(rolled_z), std::vector<std::uint8_t>(dem.get_mask()), make_meta(nr, nc));

    const auto rolled_result =
        lsdtt::compute_hillslope_newton(rolled, channels, options);
    check(rolled_result.converged, "periodic rolled Newton solve converged");

    const auto reference_result =
        lsdtt::compute_hillslope_newton(dem, channels, options);
    double worst = 0.0;
    for (int r = 0; r < nr; ++r)
        for (int c = 0; c < nc; ++c)
            worst = std::max(worst, std::fabs(
                static_cast<double>(reference_result.elevation.get_data()[static_cast<std::size_t>(r) * nc + c])
                - static_cast<double>(rolled_result.elevation.get_data()[static_cast<std::size_t>(r) * nc + ((c + shift) % nc)])));
    check(worst < 1e-6, "periodic reference is translation invariant");
    return row;
}

struct ocn_fixture
{
    lsdtt::lsdtt_raster_cpu dem;
    std::vector<double> areas;
    std::int64_t outlet = -1;
    double setup_ms = 0.0;
};

ocn_fixture make_small_ocn()
{
    lsdtt::lsdtt_ocn_options o;
    o.n_rows = 32;
    o.n_cols = 32;
    o.cell_size = 100.0;
    o.outlet_mode = lsdtt::ocn_outlet_mode::single_south;
    o.gamma = 0.5;
    o.n_iter = 20000;
    o.temperature_scale = 0.0;
    o.seed = 7;
    o.energy_length_weighted = false;
    o.relief_mode = true;
    o.target_relief = 20.0f;

    const auto start = clock_type::now();
    lsdtt::lsdtt_synthetic_network_cpu engine(o);
    engine.run();
    lsdtt::lsdtt_raster_result raster = engine.get_elevation_raster();

    ocn_fixture fixture{
        lsdtt::lsdtt_raster_cpu(std::move(raster)),
        engine.areas(),
        engine.outlet_flat(),
        elapsed_ms(start)};
    return fixture;
}

reference_row solve_ocn_coupling(
    const std::string& name,
    const ocn_fixture& fixture,
    bool imposed_channels)
{
    constexpr float dt = 1000.0f;
    constexpr float dx = 100.0f;
    constexpr double channel_pixels = 16.0;

    lsdtt::lsdtt_flowinfo_options flow_options;
    flow_options.boundary_conditions = "nnnn";
    flow_options.allow_pits_as_outlets = true;

    const float threshold = imposed_channels
        ? static_cast<float>(channel_pixels * dx * dx)
        : 0.0f;
    std::vector<double> fluvial_timings;
    fluvial_timings.reserve(static_cast<std::size_t>(g_timing_repeats));
    lsdtt::lsdtt_raster_cpu post_fluvial;
    for (int repeat = 0; repeat < g_timing_repeats; ++repeat)
    {
        const auto fluvial_start = clock_type::now();
        lsdtt::lsdtt_flowinfo flow(
            fixture.dem.get_data_array(), fixture.dem.get_mask_array(),
            fixture.dem.get_meta(), flow_options);
        post_fluvial = lsdtt::compute_d8_stream_power_incision(
            fixture.dem, flow, 1.0e-6f, 0.45f, 1.0f, dt,
            threshold, false, 1.0e-5f);
        fluvial_timings.push_back(elapsed_ms(fluvial_start));
    }
    const double fluvial_ms = median(std::move(fluvial_timings));

    std::vector<std::uint8_t> channels(fixture.areas.size(), 0u);
    if (imposed_channels)
        for (std::size_t i = 0; i < channels.size(); ++i)
            channels[i] = fixture.areas[i] >= channel_pixels ? 1u : 0u;

    auto options = newton_options(0.1, 0.0, "nnnn");
    options.inv_dt = 1.0 / static_cast<double>(dt);
    lsdtt::lsdtt_raster_cpu post_hillslope;
    reference_row row = solve_reference(
        name, post_fluvial, channels, options, post_fluvial.get_data(),
        &post_hillslope);
    row.setup_ms = fixture.setup_ms;
    row.fluvial_ms = fluvial_ms;

    const double cell_area = dx * dx;
    row.fluvial_removed_m3 =
        (elevation_sum(fixture.dem) - elevation_sum(post_fluvial)) * cell_area;
    row.hillslope_exported_m3 =
        (elevation_sum(post_fluvial) - elevation_sum(post_hillslope)) * cell_area;

    // Complete the short coupling sequence with uplift, then restore the one
    // fixed OCN outlet exactly as run_lem_step restores boundary flats.
    lsdtt::lsdtt_raster_cpu post_uplift = post_hillslope;
    const float outlet_z = post_uplift.get_data()[static_cast<std::size_t>(fixture.outlet)];
    for (std::size_t i = 0; i < post_uplift.get_data().size(); ++i)
        if (post_uplift.get_mask()[i]) post_uplift.get_data()[i] += U * dt;
    post_uplift.get_data()[static_cast<std::size_t>(fixture.outlet)] = outlet_z;
    row.uplift_added_m3 =
        (elevation_sum(post_uplift) - elevation_sum(post_hillslope)) * cell_area;

    const std::size_t n_valid = static_cast<std::size_t>(std::count(
        fixture.dem.get_mask().begin(), fixture.dem.get_mask().end(), 1u));
    const double expected_uplift_m3 =
        static_cast<double>(U) * dt * cell_area * static_cast<double>(n_valid - 1);

    const std::size_t n_channels =
        static_cast<std::size_t>(std::count(channels.begin(), channels.end(), 1u));
    if (imposed_channels)
        check(n_channels > 0 && n_channels < channels.size(),
              "small OCN imposed-channel mask is non-empty and not full-domain");
    else
        check(n_channels == 0, "small OCN coupled variant has an empty channel set");
    check(row.fluvial_removed_m3 > 0.0,
          name + " records positive fluvial removal");
    if (imposed_channels)
        check(row.hillslope_exported_m3 >= 0.0,
              "small OCN imposed channels receive non-negative hillslope export");
    else
        check(std::fabs(row.hillslope_exported_m3) <= 1.0,
              "small OCN closed hillslope step conserves volume");
    check(std::fabs(row.uplift_added_m3 - expected_uplift_m3)
              <= 1e-5 * expected_uplift_m3,
          name + " uplift budget excludes the restored outlet and closes");
    check(row.max_slope_ratio < 0.85,
          name + " remains below r_match for direct Roering comparison");
    return row;
}

void write_rows(const std::filesystem::path& path,
                const std::vector<reference_row>& rows)
{
    std::filesystem::create_directories(path.parent_path());
    std::ofstream out(path);
    if (!out) throw std::runtime_error("cannot write reference CSV: " + path.string());
    out << "fixture,cells,setup_ms,fluvial_ms,newton_ms,newton_iterations,"
           "line_search_steps,residual_m_per_yr,max_slope_ratio,elevation_sum,"
           "weighted_elevation_sum,fluvial_removed_m3,hillslope_exported_m3,"
           "uplift_added_m3,converged\n";
    out << std::setprecision(17);
    for (const reference_row& row : rows)
        out << row.name << ',' << row.cells << ',' << row.setup_ms << ','
            << row.fluvial_ms << ',' << row.newton_ms << ',' << row.iterations << ','
            << row.line_search_steps << ',' << row.residual << ','
            << row.max_slope_ratio << ',' << row.elevation_sum << ','
            << row.weighted_sum << ',' << row.fluvial_removed_m3 << ','
            << row.hillslope_exported_m3 << ',' << row.uplift_added_m3 << ','
            << (row.converged ? 1 : 0) << '\n';
}

} // namespace

int main(int argc, char** argv)
{
    using namespace lsdtt;

    const auto suite_start = clock_type::now();
    std::cout << "\n=== Small Newton-Roering cellular reference fixtures ===\n";
    if (argc > 2)
    {
        g_timing_repeats = std::max(1, std::stoi(argv[2]));
        std::cout << "Timing repeats per measured operation: "
                  << g_timing_repeats << " (CSV reports medians)\n";
    }

    // Skip cleanly in the supported Eigen-free build.
    {
        const lsdtt_raster_cpu probe_dem = make_cross_slope_dem(5, 5);
        const auto probe = compute_hillslope_newton(
            probe_dem, edge_channels(5, 5), newton_options());
        if (!probe.had_eigen)
        {
            std::cout << "SKIP: Newton reference requires Eigen\n";
            return 0;
        }
    }

    std::vector<reference_row> rows;

    {
        lsdtt_hillslope_strip_surface strip = generate_hillslope_strip_surface(
            5, 65, DX, lsdtt_strip_channel_config::both_ends,
            lsdtt_strip_init_profile::analytic_parabola, true,
            0.0, 0.0, U, D, SC, true, 0.0, 1, 1);
        const auto channels = channel_mask_from_raster(strip.channel_mask);
        perturb_strip(strip.dem, channels, 0.02);
        rows.push_back(solve_reference(
            "strip_low_5x65", strip.dem, channels, newton_options(D, U, "popo")));
        check(rows.back().max_slope_ratio < 0.75,
              "low strip remains comfortably below r_match");
    }

    {
        constexpr double near_D = 0.128;
        lsdtt_hillslope_strip_surface strip = generate_hillslope_strip_surface(
            7, 257, DX, lsdtt_strip_channel_config::both_ends,
            lsdtt_strip_init_profile::analytic_parabola, true,
            0.0, 0.0, U, near_D, SC, true, 0.0, 1, 1);
        const auto channels = channel_mask_from_raster(strip.channel_mask);
        perturb_strip(strip.dem, channels, 0.01);
        rows.push_back(solve_reference(
            "strip_near_7x257", strip.dem, channels,
            newton_options(near_D, U, "popo")));
        check(rows.back().max_slope_ratio > 0.75
              && rows.back().max_slope_ratio < 0.85,
              "long strip reaches the intended near-critical comparison band");
    }

    {
        const lsdtt_raster_cpu dem = make_cross_slope_dem(17, 33);
        rows.push_back(solve_reference(
            "cross_slope_17x33", dem, north_south_channels(17, 33),
            newton_options()));
        check(rows.back().max_slope_ratio < 0.85,
              "cross-slope fixture remains below r_match");
    }

    rows.push_back(solve_periodic_reference());

    {
        constexpr int nr = 17;
        constexpr int nc = 33;
        std::vector<float> z(static_cast<std::size_t>(nr) * nc, 0.0f);
        std::vector<std::uint8_t> mask(z.size(), 1u);
        z[static_cast<std::size_t>(nr / 2) * nc + nc / 2] = 6.0f;
        lsdtt_raster_cpu dem(std::move(z), std::move(mask), make_meta(nr, nc));
        auto options = newton_options(D, 0.0);
        options.inv_dt = 1.0 / 100.0;
        rows.push_back(solve_reference(
            "supercritical_patch_17x33", dem, edge_channels(nr, nc), options,
            dem.get_data()));
        check(rows.back().max_slope_ratio < 1.0,
              "Newton reference restores the isolated supercritical patch");
    }

    {
        lsdtt_hillslope_strip_surface strip = generate_hillslope_strip_surface(
            7, 129, DX, lsdtt_strip_channel_config::both_ends,
            lsdtt_strip_init_profile::analytic_parabola, true,
            0.0, 0.0, U, D, SC, true, 0.0, 1, 1);
        const auto channels = channel_mask_from_raster(strip.channel_mask);
        const std::vector<float> old_z = strip.dem.get_data();
        for (std::size_t i = 0; i < strip.dem.get_data().size(); ++i)
            if (channels[i]) strip.dem.get_data()[i] -= 0.25f;
        auto options = newton_options(D, 0.0, "popo");
        options.inv_dt = 1.0 / 500.0;
        rows.push_back(solve_reference(
            "incision_pulse_7x129", strip.dem, channels, options, old_z));
        check(rows.back().max_slope_ratio < 0.85,
              "incision-pulse fixture remains below r_match");
    }

    const ocn_fixture ocn = make_small_ocn();
    rows.push_back(solve_ocn_coupling("ocn_imposed_32x32", ocn, true));
    rows.push_back(solve_ocn_coupling("ocn_coupled_32x32", ocn, false));

    const std::filesystem::path output = argc > 1
        ? std::filesystem::path(argv[1]) / "cellular_hillslope_newton_references.csv"
        : std::filesystem::path("cellular_hillslope_newton_references.csv");
    write_rows(output, rows);

    const double suite_ms = elapsed_ms(suite_start);
    std::cout << "\nReference CSV: " << output.string() << "\n";
    std::cout << "Small-fixture suite wall time: " << suite_ms << " ms"
              << (suite_ms <= 20000.0 ? " (within 20 s P0 budget)" : " (OVER P0 budget)")
              << "\n";
    std::cout << (g_checks - g_failures) << " / " << g_checks << " checks passed\n";
    return g_failures == 0 ? 0 : 1;
}
