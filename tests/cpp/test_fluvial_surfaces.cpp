// test_fluvial_surfaces.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
//
// -----------------------------------------------------------------------------
// Phase 0a of docs/plan_docs/plan_fluvial_surfaces_cli_split.md.
//
// lsdtt-flat-surface-finder and lsdtt-valley-centreline shipped with NO test
// coverage of any kind (plan finding C4): no smoke test, no parameter template,
// no fixture. Every later phase of that plan changes their behaviour, so this
// driver exists to lay down a baseline BEFORE anything moves.
//
// What it writes (--write-cli-fixture MODE):
//
//   fluvsurf_dem.tif       A projected, deterministic DEM: a Gaussian ridge
//                          running N-S with flat piedmont on both flanks
//                          (generate_linear_range_surface). The contrast that
//                          matters here is roughness, not height - the mountain
//                          zone carries 150 m of fractal noise and the piedmont
//                          0.5 m, so a local-relief threshold separates them
//                          cleanly and the flat-surface finder has something
//                          real to find. The ridge also sheds a channel network
//                          onto the piedmont, which is what the valley
//                          centreline needs.
//
//                          The ridge is NARROW on purpose
//                          (k_ridge_half_width_fraction = 0.10). A wider one
//                          leaves the outer piedmont sitting on the Gaussian
//                          tail rather than on flat ground, and the regional
//                          gradient there dominates local relief: at
//                          half_width 0.20 and a 1200 m ridge the "flat"
//                          piedmont measured 30.6 m of relief against the
//                          mountains' 49.8 m, which is no contrast at all.
//                          The self-checks below assert the separation for
//                          exactly that reason.
//
//   fluvsurf_sources.fgb   One channel source, in the raster CRS.
//                          lsdtt-valley-centreline requires source points and
//                          cannot run without them.
//
// The source is chosen by a stated deterministic rule rather than a hand-picked
// coordinate, so it survives a change to the DEM dimensions: among all cells
// whose contributing-pixel count is at least k_source_min_contributing_pixels,
// take the one of greatest elevation, breaking ties by row-major order. That is
// a channel head high on the range flank.
//
// The DEM is generated from the raw-dimension overload, which stamps
// meta.set_default_synthetic_georeference() - a projected UTM CRS in metres.
// Both CLIs call check_projected_metres() and refuse a geographic raster, which
// is why tests/data/tiny30_dem.tif cannot be used here.
//
// Usage:
//   lsdtt-test-fluvial-surfaces --write-cli-fixture <output_dir>
//   lsdtt-test-fluvial-surfaces                      (self-check only)
//
// Prints one PASS/FAIL line per check and a final "RESULT: PASS|FAIL".
// -----------------------------------------------------------------------------

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "lsdtt_data_loader.hpp"
#include "lsdtt_flow_base.hpp"
#include "lsdtt_raster_cpu.hpp"
#include "lsdtt_raster_flow_preprocessing_cpu.hpp"
#include "lsdtt_raster_geometry.hpp"
#include "lsdtt_raster_surface_metrics_cpu.hpp"
#include "lsdtt_raster_synthetic_surfaces_cpu.hpp"

namespace
{

// =============================================================================
// Tiny check harness (same shape as tests/cpp/test_internal_drainage.cpp)
// =============================================================================
int g_n_checks = 0;
int g_n_failed = 0;

void check(bool condition, const std::string& what)
{
    ++g_n_checks;
    if (!condition) ++g_n_failed;
    std::cout << (condition ? "PASS  " : "FAIL  ") << what << "\n";
}

// =============================================================================
// Fixture geometry - one place, so the tests and the writer cannot disagree
// =============================================================================
constexpr int    k_n_rows    = 200;
constexpr int    k_n_cols    = 200;
constexpr double k_cell_size = 10.0;   // metres -> a 2 km square

// Roughness contrast is the point of this surface: 150 m of noise in the
// mountains against 0.5 m on the piedmont.
constexpr double k_ridge_elevation          = 600.0;
constexpr double k_piedmont_elevation       = 0.0;
constexpr double k_ridge_position_fraction  = 0.5;
constexpr double k_ridge_half_width_fraction = 0.10;
constexpr double k_mountain_noise_amplitude = 150.0;
constexpr double k_piedmont_noise_amplitude = 0.5;
constexpr double k_noise_roughness          = 0.55;
constexpr int    k_noise_seed               = 1;

// Deterministic source-selection rule (see the header comment).
constexpr std::uint64_t k_source_min_contributing_pixels = 200;

lsdtt::lsdtt_mountain_range_surface make_range_surface()
{
    return lsdtt::generate_linear_range_surface(
        k_n_rows, k_n_cols, k_cell_size,
        k_ridge_elevation,
        k_piedmont_elevation,
        k_ridge_position_fraction,
        k_ridge_half_width_fraction,
        k_mountain_noise_amplitude,
        k_piedmont_noise_amplitude,
        /*ridge_runs_north_south=*/true,
        k_noise_roughness,
        k_noise_seed,
        /*n_threads=*/1);
}

// -----------------------------------------------------------------------------
// find_source_flat
//
// The deterministic rule: highest cell among those draining at least
// k_source_min_contributing_pixels, ties broken by row-major order. Returns -1
// if no cell qualifies, which would mean the fixture geometry has drifted.
// -----------------------------------------------------------------------------
std::int64_t find_source_flat(const lsdtt::lsdtt_flow_base& flow,
                              const lsdtt::lsdtt_raster_cpu& dem)
{
    const std::vector<float>&   z    = dem.get_data();
    const std::vector<uint8_t>& mask = dem.get_mask();

    std::int64_t best_flat = -1;
    float        best_z    = -std::numeric_limits<float>::infinity();

    for (std::size_t i = 0; i < z.size(); ++i)
    {
        if (mask[i] == 0) continue;
        const std::int64_t flat = static_cast<std::int64_t>(i);
        if (flow.s_vector_size(flat) < k_source_min_contributing_pixels) continue;
        if (z[i] > best_z)
        {
            best_z    = z[i];
            best_flat = flat;
        }
    }
    return best_flat;
}

// =============================================================================
// write_cli_fixture
// =============================================================================
void write_cli_fixture(const std::string& output_dir)
{
    std::filesystem::create_directories(output_dir);

    const lsdtt::lsdtt_mountain_range_surface range = make_range_surface();
    const lsdtt::lsdtt_raster_cpu&            dem   = range.dem;

    lsdtt::lsdtt_data_loader loader;
    loader.write_full(
        dem.get_data(), dem.get_meta(),
        output_dir + "/fluvsurf_dem.tif",
        lsdtt::raster_output_format::gtiff);

    // Flow routing on the filled surface, only to locate the source point.
    const lsdtt::lsdtt_raster_cpu filled =
        lsdtt::fill_pits(dem, 0.0001f, /*n_threads=*/1);

    lsdtt::lsdtt_flowinfo_options fi_opts;
    fi_opts.boundary_conditions     = "oooo";
    fi_opts.build_topological_order = true;
    fi_opts.build_s_vector          = true;
    fi_opts.n_threads               = 1;

    const std::unique_ptr<lsdtt::lsdtt_flow_base> flow = lsdtt::make_flow_base(
        filled.get_data_array(), filled.get_mask_array(),
        filled.get_meta(), fi_opts);

    const std::int64_t source_flat = find_source_flat(*flow, filled);
    if (source_flat < 0)
    {
        throw std::runtime_error(
            "test_fluvial_surfaces: no cell reached "
            "k_source_min_contributing_pixels; the fixture geometry has drifted");
    }

    const int row = flow->row_from_flat_unchecked(source_flat);
    const int col = flow->col_from_flat_unchecked(source_flat);
    const lsdtt::lsdtt_xy xy =
        lsdtt::cell_center_from_row_col(dem.get_meta(), row, col);

    lsdtt::lsdtt_vector_schema schema;
    schema.field_names = {"source_id"};
    schema.field_types = {lsdtt::lsdtt_vector_field_type::int64_type};

    loader.begin_vector_write(
        output_dir + "/fluvsurf_sources.fgb", "fluvsurf_sources",
        schema, dem.get_meta().crs_wkt,
        lsdtt::vector_output_format::flatgeobuf, /*spatial_index=*/false);
    loader.write_point(xy.x, xy.y, {static_cast<std::int64_t>(0)});
    loader.finish_vector_write();

    std::cout << "fixture written to " << output_dir << "\n"
              << "  dem     : " << k_n_rows << "x" << k_n_cols
              << " @ " << k_cell_size << " m\n"
              << "  source  : row " << row << ", col " << col
              << "  (x=" << std::fixed << std::setprecision(1) << xy.x
              << ", y=" << xy.y << ")\n";
}

// =============================================================================
// Self-checks
//
// These guard the FIXTURE, not the CLIs. If the generated surface stops having
// a roughness contrast, or stops shedding a channel, the CLI smoke tests below
// would fail for a reason that has nothing to do with the code under test - so
// the properties they depend on are asserted here, where the diagnosis is
// immediate.
// =============================================================================
void run_self_checks()
{
    const lsdtt::lsdtt_mountain_range_surface range = make_range_surface();
    const lsdtt::lsdtt_raster_cpu&            dem   = range.dem;

    check(dem.get_meta().n_rows == k_n_rows && dem.get_meta().n_cols == k_n_cols,
          "fixture DEM has the expected dimensions");

    // Projected metres: both CLIs refuse anything else.
    bool projected = true;
    try { dem.get_meta().check_projected_metres(); }
    catch (const std::exception&) { projected = false; }
    check(projected, "fixture DEM is projected in metres");

    // The roughness contrast that the relief threshold separates. A 5-cell
    // window at 10 m is a 20 m radius, matching relief_window_pixels: 5.
    const lsdtt::lsdtt_raster_cpu relief =
        lsdtt::moving_window_local_relief(dem, 2.0 * k_cell_size);

    const std::vector<float>&   rv   = relief.get_data();
    const std::vector<uint8_t>& rm   = relief.get_mask();
    const lsdtt::lsdtt_raster_meta& m = dem.get_meta();

    // Piedmont = the outer eighth of the cross-axis (x), mountains = the
    // central tenth. Both well inside their zones, so the check does not
    // depend on exactly where the ridge flank ends.
    double piedmont_sum = 0.0, mountain_sum = 0.0;
    std::int64_t piedmont_n = 0, mountain_n = 0;

    for (int r = 1; r < m.n_rows - 1; ++r)
    {
        for (int c = 1; c < m.n_cols - 1; ++c)
        {
            const std::size_t i = static_cast<std::size_t>(r) * m.n_cols + c;
            if (rm[i] == 0) continue;
            const double f = static_cast<double>(c) / (m.n_cols - 1);
            if (f < 0.125)
            {
                piedmont_sum += rv[i];
                ++piedmont_n;
            }
            else if (f > 0.45 && f < 0.55)
            {
                mountain_sum += rv[i];
                ++mountain_n;
            }
        }
    }

    check(piedmont_n > 0 && mountain_n > 0,
          "fixture has cells in both the piedmont and the mountain zone");

    const double piedmont_relief = piedmont_n ? piedmont_sum / piedmont_n : 0.0;
    const double mountain_relief = mountain_n ? mountain_sum / mountain_n : 0.0;

    std::cout << "      piedmont mean relief = " << std::fixed
              << std::setprecision(3) << piedmont_relief << " m\n"
              << "      mountain mean relief = "
              << mountain_relief << " m\n";

    check(piedmont_relief < 2.0,
          "piedmont mean local relief is below the 2 m flat threshold");
    check(mountain_relief > 10.0,
          "mountain mean local relief is well above the flat threshold");
    check(mountain_relief > 5.0 * piedmont_relief,
          "the roughness contrast is at least fivefold");

    // The channel the valley centreline traces has to exist.
    const lsdtt::lsdtt_raster_cpu filled =
        lsdtt::fill_pits(dem, 0.0001f, /*n_threads=*/1);

    lsdtt::lsdtt_flowinfo_options fi_opts;
    fi_opts.boundary_conditions     = "oooo";
    fi_opts.build_topological_order = true;
    fi_opts.build_s_vector          = true;
    fi_opts.n_threads               = 1;

    const std::unique_ptr<lsdtt::lsdtt_flow_base> flow = lsdtt::make_flow_base(
        filled.get_data_array(), filled.get_mask_array(),
        filled.get_meta(), fi_opts);

    const std::int64_t source_flat = find_source_flat(*flow, filled);
    check(source_flat >= 0,
          "a channel source meeting the contributing-pixel rule exists");

    if (source_flat >= 0)
    {
        const int row = flow->row_from_flat_unchecked(source_flat);
        const int col = flow->col_from_flat_unchecked(source_flat);
        std::cout << "      source at row " << row << ", col " << col
                  << ", contributing pixels = "
                  << flow->s_vector_size(source_flat) << "\n";
        check(row > 0 && row < k_n_rows - 1 && col > 0 && col < k_n_cols - 1,
              "the source is interior, not on a boundary cell");
    }
}

}  // namespace

int main(int argc, char* argv[])
{
    try
    {
        if (argc == 3 && std::string(argv[1]) == "--write-cli-fixture")
        {
            write_cli_fixture(argv[2]);
            return 0;
        }

        std::cout << "lsdtt-test-fluvial-surfaces\n"
                  << "Phase 0a fixture self-checks\n\n";
        run_self_checks();
    }
    catch (const std::exception& e)
    {
        std::cerr << "EXCEPTION: " << e.what() << "\n";
        return 1;
    }

    std::cout << "\n" << (g_n_checks - g_n_failed) << " / " << g_n_checks
              << " checks passed\n"
              << "RESULT: " << (g_n_failed == 0 ? "PASS" : "FAIL") << "\n";
    return g_n_failed == 0 ? 0 : 1;
}
