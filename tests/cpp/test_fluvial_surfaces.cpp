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
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

#include "lsdtt_data_loader.hpp"
#include "lsdtt_flow_base.hpp"
#include "lsdtt_raster_cpu.hpp"
#include "lsdtt_raster_flow_preprocessing_cpu.hpp"
#include "lsdtt_raster_linear_features_cpu.hpp"
#include "lsdtt_raster_geometry.hpp"
#include "lsdtt_raster_surface_metrics_cpu.hpp"
#include "lsdtt_raster_synthetic_surfaces_cpu.hpp"
#include "lsdtt_reference_field_ladder_cpu.hpp"
#include "lsdtt_valley_centreline_cpu.hpp"
#include "lsdtt_channel_scope_cpu.hpp"
#include "lsdtt_flat_surface_channel_classify_cpu.hpp"

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

// =============================================================================
// Phase 0b - the terrace staircase and alluvial fan generators
//
// These assert the GEOMETRY the later tests depend on. T1 (dissected treads),
// T3 (fan not shredded) and T6 (reach-relative beats global) are all written
// against properties claimed in the generator's header comment, so if a claim
// is wrong the diagnosis belongs here and not four phases downstream.
// =============================================================================

// Staircase fixture geometry. 4 terraces x 300 m either side of the channel is
// 2400 m of flight, so a 300-column grid at 10 m spans it with valley wall to
// spare on both flanks.
constexpr int    k_st_n_rows      = 200;
constexpr int    k_st_n_cols      = 300;
constexpr double k_st_cell_size   = 10.0;
constexpr int    k_st_n_terraces  = 4;
constexpr double k_st_tread_width = 300.0;
constexpr double k_st_riser       = 8.0;
constexpr double k_st_valley_grad = 0.002;
constexpr double k_st_cross_grad  = 0.002;
// Must exceed the 40 m relief window several times over, or the window straddles
// the gap and reads the risers either side of it. See the generator header.
constexpr double k_st_gap_width   = 120.0;

void run_terrace_staircase_checks()
{
    std::cout << "\n-- terrace staircase --\n";

    // Exact surface, no noise, no cross-fall: tread elevations are then single
    // values that can be counted.
    const lsdtt::lsdtt_terrace_staircase_surface st =
        lsdtt::generate_terrace_staircase_surface(
            k_st_n_rows, k_st_n_cols, k_st_cell_size,
            k_st_n_terraces, k_st_tread_width, k_st_riser, k_st_valley_grad,
            /*taper=*/0.0, /*cross_gradient=*/0.0,
            /*n_riser_gaps=*/0, /*riser_gap_width_m=*/30.0,
            /*noise=*/0.0, /*channel_elevation=*/100.0,
            /*valley_wall_gradient=*/0.5, /*seed=*/1, /*n_threads=*/1);

    const lsdtt::lsdtt_raster_meta& m = st.dem.get_meta();
    const std::vector<float>& z  = st.dem.get_data();
    const std::vector<float>& tid = st.tread_id.get_data();
    const std::vector<float>& ch  = st.channel_mask.get_data();

    // The plan's exit criterion: exactly n_terraces distinct tread elevations
    // at a fixed downstream station.
    const int probe_row = k_st_n_rows / 2;
    std::vector<double> tread_levels;
    for (int c = 0; c < m.n_cols; ++c)
    {
        const std::size_t i = static_cast<std::size_t>(probe_row) * m.n_cols + c;
        if (tid[i] <= 0.0f) continue;
        const double zv = z[i];
        bool seen = false;
        for (const double t : tread_levels)
        {
            if (std::abs(t - zv) < 1e-3) { seen = true; break; }
        }
        if (!seen) tread_levels.push_back(zv);
    }
    std::sort(tread_levels.begin(), tread_levels.end());

    std::cout << "      distinct tread elevations at row " << probe_row << ":";
    for (const double t : tread_levels) std::cout << " " << t;
    std::cout << "\n";

    check(static_cast<int>(tread_levels.size()) == k_st_n_terraces,
          "staircase has exactly n_terraces distinct tread elevations");

    // Consecutive treads are exactly one riser apart.
    bool risers_ok = tread_levels.size() >= 2;
    for (std::size_t i = 1; i < tread_levels.size(); ++i)
    {
        if (std::abs((tread_levels[i] - tread_levels[i - 1]) - k_st_riser) > 1e-3)
        {
            risers_ok = false;
        }
    }
    check(risers_ok, "consecutive tread elevations differ by exactly riser_height_m");

    // The channel is one cell wide and every tread id is in range.
    std::int64_t n_channel = 0;
    bool tread_ids_in_range = true;
    for (std::size_t i = 0; i < z.size(); ++i)
    {
        if (ch[i] > 0.5f) ++n_channel;
        const int t = static_cast<int>(tid[i]);
        if (t < 0 || t > k_st_n_terraces) tread_ids_in_range = false;
    }
    check(n_channel == k_st_n_rows,
          "the channel is exactly one cell wide down the whole valley");
    check(tread_ids_in_range, "every tread id is within 0..n_terraces");

    // --- the taper, which is what test T6 rests on ---
    const lsdtt::lsdtt_terrace_staircase_surface tp =
        lsdtt::generate_terrace_staircase_surface(
            k_st_n_rows, k_st_n_cols, k_st_cell_size,
            k_st_n_terraces, k_st_tread_width, /*riser=*/12.0, k_st_valley_grad,
            /*taper=*/0.667, /*cross_gradient=*/0.0,
            0, 30.0, 0.0, 100.0, 0.5, 1, 1);

    auto first_tread_height_above_channel = [&](const lsdtt::lsdtt_raster_cpu& dem,
                                                const lsdtt::lsdtt_raster_cpu& tids,
                                                int row) -> double
    {
        const lsdtt::lsdtt_raster_meta& mm = dem.get_meta();
        const std::vector<float>& zz = dem.get_data();
        const std::vector<float>& tt = tids.get_data();
        const std::size_t centre =
            static_cast<std::size_t>(row) * mm.n_cols + (mm.n_cols - 1) / 2;
        double z_ch = zz[centre];
        for (int c = (mm.n_cols - 1) / 2; c < mm.n_cols; ++c)
        {
            const std::size_t i = static_cast<std::size_t>(row) * mm.n_cols + c;
            if (static_cast<int>(tt[i]) == 1) return zz[i] - z_ch;
        }
        return 0.0;
    };

    const double h_up   = first_tread_height_above_channel(tp.dem, tp.tread_id, 2);
    const double h_down = first_tread_height_above_channel(
        tp.dem, tp.tread_id, k_st_n_rows - 3);

    std::cout << "      T1 height above channel: upstream " << h_up
              << " m, downstream " << h_down << " m\n";
    check(h_up > 11.5 && h_up < 12.1,
          "tapered staircase starts near riser_height_m upstream");
    check(h_down > 3.9 && h_down < 4.6,
          "tapered staircase falls to about a third of it downstream");

    // --- the riser gap, which is what test T1 rests on ---
    //
    // The claim in the generator header is that the gap ramp reads as FLAT to a
    // roughness threshold while the risers around it do not. That is the whole
    // mechanism by which a dissected riser leaks, so it is asserted rather than
    // assumed.
    const lsdtt::lsdtt_terrace_staircase_surface gp =
        lsdtt::generate_terrace_staircase_surface(
            k_st_n_rows, k_st_n_cols, k_st_cell_size,
            k_st_n_terraces, k_st_tread_width, k_st_riser, k_st_valley_grad,
            0.0, k_st_cross_grad,
            /*n_riser_gaps=*/3, /*riser_gap_width_m=*/k_st_gap_width,
            0.0, 100.0, 0.5, 1, 1);

    const lsdtt::lsdtt_raster_cpu gap_relief =
        lsdtt::moving_window_local_relief(gp.dem, 2.0 * k_st_cell_size);

    const std::vector<float>& gr  = gap_relief.get_data();
    const std::vector<float>& gid = gp.tread_id.get_data();

    // Gap rows carry tread_id 0 across the whole flight; a riser row carries
    // tread_id 0 only in the one-cell step between treads. Count zero-tread
    // cells per row inside the flight to tell them apart.
    std::vector<int> n_zero_per_row(static_cast<std::size_t>(k_st_n_rows), 0);
    for (int r = 0; r < k_st_n_rows; ++r)
    {
        int n_zero = 0;
        for (int c = 0; c < k_st_n_cols; ++c)
        {
            const std::size_t i = static_cast<std::size_t>(r) * k_st_n_cols + c;
            const double x = std::abs(c - 0.5 * (k_st_n_cols - 1)) * k_st_cell_size;
            if (x < k_st_n_terraces * k_st_tread_width && gid[i] == 0.0f) ++n_zero;
        }
        n_zero_per_row[static_cast<std::size_t>(r)] = n_zero;
    }

    // Take a row well INSIDE a gap, not at its leading edge. A relief window
    // centred one row inside the gap still reaches the stepped rows outside it
    // and returns their riser relief -- which is the same reason a gap narrower
    // than the window never reads as flat at all. Require the whole window
    // depth either side to be gap as well.
    const int half_window_rows =
        static_cast<int>(std::ceil(2.0 * k_st_cell_size / k_st_cell_size));
    int gap_row = -1, plain_row = -1;
    for (int r = half_window_rows; r < k_st_n_rows - half_window_rows; ++r)
    {
        if (plain_row < 0 && n_zero_per_row[static_cast<std::size_t>(r)] < 20)
        {
            plain_row = r;
        }
        if (gap_row >= 0) continue;

        bool window_all_gap = true;
        for (int d = -half_window_rows; d <= half_window_rows; ++d)
        {
            if (n_zero_per_row[static_cast<std::size_t>(r + d)] <= 100)
            {
                window_all_gap = false;
                break;
            }
        }
        if (window_all_gap) gap_row = r;
    }

    check(gap_row >= 0 && plain_row >= 0,
          "the staircase has both gap rows and ordinary stepped rows");

    if (gap_row >= 0 && plain_row >= 0)
    {
        // Sample a column that sits inside the flight and away from the channel.
        const int probe_col = (k_st_n_cols - 1) / 2 +
            static_cast<int>(1.5 * k_st_tread_width / k_st_cell_size);

        const double gap_relief_v =
            gr[static_cast<std::size_t>(gap_row) * k_st_n_cols + probe_col];

        // The steepest thing on a plain row: the riser itself.
        double max_plain_relief = 0.0;
        for (int c = 0; c < k_st_n_cols; ++c)
        {
            const std::size_t i =
                static_cast<std::size_t>(plain_row) * k_st_n_cols + c;
            const double x = std::abs(c - 0.5 * (k_st_n_cols - 1)) * k_st_cell_size;
            if (x < k_st_n_terraces * k_st_tread_width)
            {
                max_plain_relief = std::max(max_plain_relief,
                                            static_cast<double>(gr[i]));
            }
        }

        std::cout << "      gap-row relief at the probe column = "
                  << gap_relief_v << " m\n"
                  << "      max riser relief on a plain row     = "
                  << max_plain_relief << " m\n";

        check(gap_relief_v < 2.0,
              "the gap ramp reads as flat at a 2 m relief threshold");
        check(max_plain_relief > 4.0,
              "the risers do not, so only the gap can leak");
    }

    // --- drainage, which every HAND-based test rests on ---
    //
    // The header warns that a zero cross-fall leaves treads draining downstream
    // off the grid rather than to the channel. Assert that the warning is real
    // AND that the documented remedy works, so neither can rot unnoticed.
    auto fraction_reaching_channel = [](const lsdtt::lsdtt_terrace_staircase_surface& s)
        -> double
    {
        const lsdtt::lsdtt_raster_cpu filled =
            lsdtt::fill_pits(s.dem, 0.0001f, /*n_threads=*/1);

        lsdtt::lsdtt_flowinfo_options o;
        o.boundary_conditions     = "oooo";
        o.build_topological_order = true;
        o.build_s_vector          = true;
        o.n_threads               = 1;

        const std::unique_ptr<lsdtt::lsdtt_flow_base> f = lsdtt::make_flow_base(
            filled.get_data_array(), filled.get_mask_array(),
            filled.get_meta(), o);

        const std::vector<float>& cm = s.channel_mask.get_data();
        std::vector<uint8_t> chan_mask(cm.size(), 0);
        for (std::size_t i = 0; i < cm.size(); ++i)
            chan_mask[i] = cm[i] > 0.5f ? 1 : 0;

        const lsdtt::height_above_channel_donor_result hd =
            f->get_height_above_channel_mask_donor_raster(
                s.dem.get_data(), s.dem.get_mask(), chan_mask);

        const std::vector<uint8_t>& hm = hd.height.mask;
        std::int64_t n_def = 0, n_tot = 0;
        for (std::size_t i = 0; i < hm.size(); ++i)
        {
            ++n_tot;
            if (hm[i] != 0) ++n_def;
        }
        return n_tot ? static_cast<double>(n_def) / n_tot : 0.0;
    };

    const lsdtt::lsdtt_terrace_staircase_surface flat_treads =
        lsdtt::generate_terrace_staircase_surface(
            k_st_n_rows, k_st_n_cols, k_st_cell_size,
            k_st_n_terraces, k_st_tread_width, k_st_riser, k_st_valley_grad,
            0.0, /*cross_gradient=*/0.0, 0, 30.0, 0.0, 100.0, 0.5, 1, 1);

    const double f_flat  = fraction_reaching_channel(flat_treads);
    const double f_tilted = fraction_reaching_channel(gp);

    std::cout << "      HAND defined, zero cross-fall    = "
              << std::setprecision(3) << 100.0 * f_flat << " %\n"
              << "      HAND defined, 0.002 cross-fall   = "
              << 100.0 * f_tilted << " %\n";

    check(f_tilted > 0.90,
          "a positive tread_cross_gradient drains the treads to the channel");
    check(f_tilted > f_flat,
          "and it drains strictly more of the grid than a zero cross-fall");
}

// Fan fixture geometry. A 2 km grid with a 0.03 radial gradient gives roughly
// 60 m of climb from toe to apex, which is what makes a HAND ladder shred it.
constexpr int    k_fan_n_rows    = 200;
constexpr int    k_fan_n_cols    = 200;
constexpr double k_fan_cell_size = 10.0;
constexpr double k_fan_gradient  = 0.03;
constexpr double k_fan_inset     = 4.0;
constexpr int    k_fan_n_lobes   = 2;

void run_alluvial_fan_checks()
{
    std::cout << "\n-- alluvial fan --\n";

    const lsdtt::lsdtt_alluvial_fan_surface fan =
        lsdtt::generate_alluvial_fan_surface(
            k_fan_n_rows, k_fan_n_cols, k_fan_cell_size,
            /*apex_x_fraction=*/0.5, /*apex_y_fraction=*/0.95,
            /*apex_elevation=*/500.0, k_fan_gradient, /*concavity=*/0.0,
            k_fan_n_lobes, k_fan_inset, /*lobe_half_angle_deg=*/25.0,
            /*noise=*/0.0, /*seed=*/1, /*n_threads=*/1);

    const lsdtt::lsdtt_raster_meta& m = fan.dem.get_meta();
    const std::vector<float>& z  = fan.dem.get_data();
    const std::vector<float>& lid = fan.lobe_id.get_data();

    std::int64_t n_old = 0;
    std::vector<std::int64_t> n_lobe(
        static_cast<std::size_t>(k_fan_n_lobes) + 1, 0);
    for (std::size_t i = 0; i < lid.size(); ++i)
    {
        const int l = static_cast<int>(lid[i]);
        if (l == 0) ++n_old;
        else if (l >= 1 && l <= k_fan_n_lobes) ++n_lobe[static_cast<std::size_t>(l)];
    }

    std::cout << "      older surface cells = " << n_old << ", lobes =";
    for (int l = 1; l <= k_fan_n_lobes; ++l)
        std::cout << " " << n_lobe[static_cast<std::size_t>(l)];
    std::cout << "\n";

    bool all_lobes_present = true;
    for (int l = 1; l <= k_fan_n_lobes; ++l)
        if (n_lobe[static_cast<std::size_t>(l)] < 100) all_lobes_present = false;
    check(all_lobes_present, "every inset lobe covers a usable area");
    check(n_old > 0, "the older fan surface survives between the lobes");

    // The surface falls away from the apex.
    const double apex_col = 0.5 * (m.n_cols - 1);
    const double apex_row = 0.95 * (m.n_rows - 1);
    const std::size_t near_i =
        static_cast<std::size_t>(static_cast<int>(apex_row) - 10) * m.n_cols
        + static_cast<int>(apex_col);
    const std::size_t far_i =
        static_cast<std::size_t>(5) * m.n_cols + static_cast<int>(apex_col);
    check(z[near_i] > z[far_i], "the fan falls from apex to toe");

    // THE property that makes a second reference field necessary: the radial
    // climb across the fan is many times the inset that separates the lobes,
    // so a ladder banding height-above-channel cuts arcs across lobe
    // boundaries instead of following them.
    const double radial_climb = z[near_i] - z[far_i];
    std::cout << "      radial climb toe->apex = " << radial_climb
              << " m, lobe inset = " << k_fan_inset << " m\n";
    check(radial_climb > 6.0 * k_fan_inset,
          "the radial climb dwarfs the lobe inset, which is why a HAND ladder "
          "shreds a fan");

    // And the property that makes the fitted-surface field work: a lobe sits a
    // CONSTANT depth below the older surface, so a residual separates them by a
    // single offset rather than by position on the fan.
    //
    // Compare each lobe cell with the older-surface elevation at the same
    // radius, reconstructed analytically.
    double max_inset_error = 0.0;
    std::int64_t n_checked = 0;
    for (int r = 5; r < m.n_rows - 5; r += 3)
    {
        for (int c = 5; c < m.n_cols - 5; c += 3)
        {
            const std::size_t i = static_cast<std::size_t>(r) * m.n_cols + c;
            if (static_cast<int>(lid[i]) == 0) continue;
            const double dx = (c - apex_col) * k_fan_cell_size;
            const double dy = (r - apex_row) * k_fan_cell_size;
            const double rr = std::sqrt(dx * dx + dy * dy);
            const double z_old = 500.0 - k_fan_gradient * rr;
            const double inset = z_old - z[i];
            max_inset_error = std::max(max_inset_error,
                                       std::abs(inset - k_fan_inset));
            ++n_checked;
        }
    }
    std::cout << "      max |lobe inset - " << k_fan_inset << " m| over "
              << n_checked << " cells = " << max_inset_error << " m\n";
    check(n_checked > 0 && max_inset_error < 1e-2,
          "every lobe cell sits exactly lobe_inset_m below the older surface");

    // Concavity must produce a genuinely non-planar profile, or the poly2
    // branch of the fitted-surface field has nothing to distinguish it.
    const lsdtt::lsdtt_alluvial_fan_surface concave =
        lsdtt::generate_alluvial_fan_surface(
            k_fan_n_rows, k_fan_n_cols, k_fan_cell_size,
            0.5, 0.95, 500.0, k_fan_gradient, /*concavity=*/0.015,
            0, k_fan_inset, 25.0, 0.0, 1, 1);

    const std::vector<float>& zc = concave.dem.get_data();
    // Three points down the fan axis. A cone is linear in r and its midpoint
    // lands exactly on the chord. A geomorphically CONCAVE fan profile - steep
    // at the apex, slackening toward the toe - is a mathematically CONVEX
    // function of radius, so its midpoint sits BELOW the chord. The naming
    // clash is a standing trap: "concave profile" is a statement about the
    // shape as drawn, not about the sign of the second derivative.
    const int col_axis = static_cast<int>(apex_col);
    const std::size_t p0 = static_cast<std::size_t>(180) * m.n_cols + col_axis;
    const std::size_t p1 = static_cast<std::size_t>(95)  * m.n_cols + col_axis;
    const std::size_t p2 = static_cast<std::size_t>(10)  * m.n_cols + col_axis;
    const double chord_mid = 0.5 * (zc[p0] + zc[p2]);
    std::cout << "      concave profile: midpoint " << zc[p1]
              << " m vs chord midpoint " << chord_mid << " m\n";
    check(zc[p1] < chord_mid - 0.5,
          "fan_concavity bows the profile below its chord (not a cone)");
}

// =============================================================================
// write_phase0b_fixtures
// =============================================================================
void write_phase0b_fixtures(const std::string& output_dir)
{
    std::filesystem::create_directories(output_dir);
    lsdtt::lsdtt_data_loader loader;

    // The dissected staircase: 4 treads, 3 road-like gaps, treads tilted so the
    // grid drains to the channel and HAND is defined.
    const lsdtt::lsdtt_terrace_staircase_surface st =
        lsdtt::generate_terrace_staircase_surface(
            k_st_n_rows, k_st_n_cols, k_st_cell_size,
            k_st_n_terraces, k_st_tread_width, k_st_riser, k_st_valley_grad,
            /*taper=*/0.0, k_st_cross_grad,
            /*n_riser_gaps=*/3, /*riser_gap_width_m=*/k_st_gap_width,
            /*noise=*/0.0, 100.0, 0.5, 1, 1);

    loader.write_full(st.dem.get_data(), st.dem.get_meta(),
                      output_dir + "/staircase_dem.tif",
                      lsdtt::raster_output_format::gtiff);
    loader.write_full(st.tread_id.get_data(), st.tread_id.get_meta(),
                      output_dir + "/staircase_tread_id.tif",
                      lsdtt::raster_output_format::gtiff);
    loader.write_full(st.channel_mask.get_data(), st.channel_mask.get_meta(),
                      output_dir + "/staircase_channel_mask.tif",
                      lsdtt::raster_output_format::gtiff);

    const lsdtt::lsdtt_alluvial_fan_surface fan =
        lsdtt::generate_alluvial_fan_surface(
            k_fan_n_rows, k_fan_n_cols, k_fan_cell_size,
            0.5, 0.95, 500.0, k_fan_gradient, 0.0,
            k_fan_n_lobes, k_fan_inset, 25.0, 0.0, 1, 1);

    loader.write_full(fan.dem.get_data(), fan.dem.get_meta(),
                      output_dir + "/fan_dem.tif",
                      lsdtt::raster_output_format::gtiff);
    loader.write_full(fan.lobe_id.get_data(), fan.lobe_id.get_meta(),
                      output_dir + "/fan_lobe_id.tif",
                      lsdtt::raster_output_format::gtiff);

    std::cout << "phase 0b fixtures written to " << output_dir << "\n";
}

// =============================================================================
// Phase 0c - the reference-field ladder engine
//
// T2, T4, T5 and T7 of section 7 of the plan.
// =============================================================================

// HAND above the staircase's own trunk channel, plus the flow-distance raster
// and per-pixel source id the reach-relative path needs.
struct staircase_hand
{
    lsdtt::lsdtt_raster_cpu          hand;
    lsdtt::lsdtt_raster_cpu          flow_distance;
    lsdtt::lsdtt_raster_cpu          admitted;
    std::vector<std::int32_t>        source_id;
    std::vector<std::int32_t>        channel_source_id;
    std::vector<std::int64_t>        donor_flat;
    std::vector<std::uint8_t>        channel_mask;
};

staircase_hand hand_for_staircase(
    const lsdtt::lsdtt_terrace_staircase_surface& s)
{
    const lsdtt::lsdtt_raster_cpu filled =
        lsdtt::fill_pits(s.dem, 0.0001f, /*n_threads=*/1);

    lsdtt::lsdtt_flowinfo_options o;
    o.boundary_conditions     = "oooo";
    o.build_topological_order = true;
    o.build_s_vector          = true;
    o.n_threads               = 1;

    const std::unique_ptr<lsdtt::lsdtt_flow_base> f = lsdtt::make_flow_base(
        filled.get_data_array(), filled.get_mask_array(), filled.get_meta(), o);

    const std::vector<float>& cm = s.channel_mask.get_data();
    std::vector<uint8_t> chan(cm.size(), 0);
    for (std::size_t i = 0; i < cm.size(); ++i) chan[i] = cm[i] > 0.5f ? 1 : 0;

    lsdtt::height_above_channel_donor_result hd =
        f->get_height_above_channel_mask_donor_raster(
            s.dem.get_data(), s.dem.get_mask(), chan);

    staircase_hand out;
    out.hand          = lsdtt::lsdtt_raster_cpu(std::move(hd.height));
    out.flow_distance = lsdtt::lsdtt_raster_cpu(
        f->get_distance_from_outlet_raster());

    // Admit everything with a defined HAND that is not the channel itself.
    const std::vector<float>&   hv = out.hand.get_data();
    const std::vector<uint8_t>& hm = out.hand.get_mask();
    std::vector<float>   adm(hv.size(), 0.0f);
    std::vector<uint8_t> adm_mask(hv.size(), 1);
    for (std::size_t i = 0; i < hv.size(); ++i)
    {
        if (hm[i] && std::isfinite(hv[i]) && chan[i] == 0) adm[i] = 1.0f;
    }
    out.admitted = lsdtt::lsdtt_raster_cpu(
        std::move(adm), std::move(adm_mask), out.hand.get_meta());

    // One channel: every pixel with a donor belongs to source 0.
    out.source_id.assign(hv.size(), -1);
    out.channel_source_id.assign(hv.size(), -1);
    out.channel_mask = chan;
    out.donor_flat = hd.donor_flat;
    for (std::size_t i = 0; i < hv.size(); ++i)
    {
        if (hd.donor_flat[i] >= 0) out.source_id[i] = 0;
        if (chan[i]) out.channel_source_id[i] = 0;
    }
    return out;
}


// =============================================================================
// Phase 1 - T1: band-then-label keeps dissected treads separate
//
// The standing assertions are on the NEW path only. The old path's leak is
// measured once and reported, but is deliberately not asserted: it disappears
// with the label_within_bands gate, and a test that pins a defect breaks when
// the defect is deleted.
// =============================================================================
struct t1_measurement
{
    std::int64_t n_labels                = 0;
    std::int64_t largest_label_pixels    = 0;
    int          largest_label_treads    = 0;
    int          max_treads_in_any_label = 0;
    double       purity                  = 0.0;   // fraction of labelled tread
                                                  // pixels in single-tread labels
};

// How many distinct true treads does each label touch, and how pure is the
// labelling overall?
t1_measurement measure_labels(const lsdtt::lsdtt_raster_cpu& labels,
                              const lsdtt::lsdtt_raster_cpu& tread_id)
{
    const std::vector<float>& lv = labels.get_data();
    const std::vector<float>& tv = tread_id.get_data();

    std::unordered_map<std::int64_t, std::unordered_map<int, std::int64_t>> hits;
    std::unordered_map<std::int64_t, std::int64_t> sizes;

    for (std::size_t i = 0; i < lv.size(); ++i)
    {
        if (lv[i] <= 0.0f) continue;
        const std::int64_t l = static_cast<std::int64_t>(lv[i]);
        ++sizes[l];
        const int t = static_cast<int>(tv[i]);
        if (t > 0) ++hits[l][t];
    }

    t1_measurement m;
    m.n_labels = static_cast<std::int64_t>(sizes.size());

    std::int64_t pure_px = 0, total_px = 0;
    for (const auto& [l, per_tread] : hits)
    {
        const int n_treads = static_cast<int>(per_tread.size());
        m.max_treads_in_any_label = std::max(m.max_treads_in_any_label, n_treads);
        std::int64_t sum = 0;
        for (const auto& [t, c] : per_tread) sum += c;
        total_px += sum;
        if (n_treads == 1) pure_px += sum;
    }
    m.purity = total_px ? static_cast<double>(pure_px) / total_px : 0.0;

    for (const auto& [l, sz] : sizes)
    {
        if (sz > m.largest_label_pixels)
        {
            m.largest_label_pixels = sz;
            auto it = hits.find(l);
            m.largest_label_treads =
                (it == hits.end()) ? 0 : static_cast<int>(it->second.size());
        }
    }
    return m;
}

void run_band_then_label_checks()
{
    std::cout << "lsdtt-test-fluvial-surfaces\n"
              << "Phase 1 band-then-label (T1)\n\n";

    // The dissected staircase: 4 treads, 3 gaps wide enough to read as flat.
    // A gentle valley so the cross-fall wins and HAND is a tread height (see
    // the generator header).
    const lsdtt::lsdtt_terrace_staircase_surface st =
        lsdtt::generate_terrace_staircase_surface(
            k_st_n_rows, k_st_n_cols, k_st_cell_size,
            k_st_n_terraces, k_st_tread_width, k_st_riser,
            /*valley_gradient=*/0.0002, 0.0,
            0.6 * k_st_riser / k_st_tread_width,
            /*n_riser_gaps=*/3, k_st_gap_width,
            0.0, 100.0, 0.5, 1, 1);

    staircase_hand sh = hand_for_staircase(st);

    // Flat ground = everything off the risers and off the wall. Emulating what
    // the flat-surface finder admits: the treads AND the gap ramps, which is
    // precisely the set in which the gaps can connect two treads.
    const lsdtt::lsdtt_raster_cpu relief =
        lsdtt::moving_window_local_relief(st.dem, 2.0 * k_st_cell_size);
    {
        const std::vector<float>& rv = relief.get_data();
        std::vector<float>&       ad = sh.admitted.get_data();
        for (std::size_t i = 0; i < ad.size(); ++i)
            if (!(rv[i] < 2.0f)) ad[i] = 0.0f;
    }

    // One label over the whole admitted set: what the OLD ordering produces
    // before any height is consulted.
    const lsdtt::banded_label_result whole =
        lsdtt::label_components_within_bands(
            /*band_index=*/sh.admitted,     // a single band: every admitted
                                            // pixel carries band 0
            sh.admitted, sh.hand, 50, lsdtt::pixel_connectivity::d4);

    // Careful: passing `admitted` as the band index gives band 1 for admitted
    // pixels, which is still ONE band, so this is exactly a whole-raster
    // labelling.
    const t1_measurement old_m = measure_labels(whole.labels, st.tread_id);

    std::cout << "-- measurement: the OLD ordering on the dissected flight --\n"
              << "      labels = " << old_m.n_labels
              << ", largest spans " << old_m.largest_label_treads
              << " true treads (" << old_m.largest_label_pixels << " px)\n"
              << "      worst label spans " << old_m.max_treads_in_any_label
              << " treads, tread purity = "
              << std::setprecision(4) << old_m.purity << "\n"
              << "      (recorded, not asserted: this path is deleted in phase 6)\n\n";

    // ---- the NEW path ----
    lsdtt::flat_surface_classify_request request;
    request.label_within_bands = true;
    request.ladder_break_algorithm = lsdtt::ladder_break_choice::jenks;
    request.ladder_level_count = k_st_n_terraces;
    request.band_reach_window_metres = -1.0;
    request.vertical_accuracy_metres = 0.5;
    request.band_min_component_pixels = 50;
    request.floodplain_max_height_metres = 5.0;
    request.terrace_min_height_metres = 5.0;

    const lsdtt::flat_surface_classify_result classified =
        lsdtt::classify_flat_surfaces_by_channel(
            whole.labels, sh.hand, sh.channel_mask, sh.channel_source_id,
            request, &sh.donor_flat, &sh.flow_distance);

    const t1_measurement new_m =
        measure_labels(classified.banded_labels_raster, st.tread_id);

    std::set<int> bands;
    const std::vector<float>& labels = classified.banded_labels_raster.get_data();
    const std::vector<float>& band_index = classified.hand_band_raster.get_data();
    for (std::size_t i = 0; i < labels.size(); ++i)
        if (labels[i] > 0.0f && band_index[i] >= 0.0f)
            bands.insert(static_cast<int>(band_index[i]));

    std::map<std::pair<int, int>, std::set<std::int64_t>> true_patch_labels;
    std::unordered_map<std::int64_t, std::set<std::pair<int, int>>> label_patches;
    const std::vector<float>& truth = st.tread_id.get_data();
    for (std::size_t i = 0; i < labels.size(); ++i)
    {
        if (labels[i] <= 0.0f || truth[i] <= 0.0f) continue;
        const int tread = static_cast<int>(truth[i]);
        const int col = static_cast<int>(i % static_cast<std::size_t>(k_st_n_cols));
        const std::pair<int, int> key{tread, col < k_st_n_cols / 2 ? 0 : 1};
        const std::int64_t label = static_cast<std::int64_t>(labels[i]);
        true_patch_labels[key].insert(label);
        label_patches[label].insert(key);
    }

    bool one_label_per_true_patch = true;
    for (const auto& [key, patch_labels] : true_patch_labels)
        one_label_per_true_patch = one_label_per_true_patch && patch_labels.size() == 1;
    bool one_true_patch_per_label = true;
    for (const auto& [label, patches] : label_patches)
        one_true_patch_per_label = one_true_patch_per_label && patches.size() == 1;

    std::cout << "-- T1: the NEW ordering --\n"
              << "      bands = " << bands.size()
              << ", labels = " << new_m.n_labels
              << ", worst label spans " << new_m.max_treads_in_any_label
              << " treads, tread purity = " << new_m.purity << "\n";

    check(classified.used_band_then_label,
          "the classifier reports that the band-then-label path ran");
    check(bands == std::set<int>({0, 1, 2, 3}),
          "the ladder recovers exactly n_terraces bands");
    check(new_m.max_treads_in_any_label == 1,
          "no label spans more than one true tread");
    check(new_m.purity >= 0.999,
          "tread purity is 1.0 under band-then-label");
    check(true_patch_labels.size() == 2 * k_st_n_terraces,
          "the truth raster contains two spatial patches per tread");
    check(one_label_per_true_patch && one_true_patch_per_label,
          "component count per band equals the spatially separate tread patches");
    check(classified.patches.size() == classified.labels.size(),
          "every replacement patch has one classification record");
    bool metadata_agrees = true;
    for (const lsdtt::flat_surface_label_class& rec : classified.labels)
    {
        if (rec.label <= 0 ||
            static_cast<std::size_t>(rec.label) > classified.patches.size())
        {
            metadata_agrees = false;
            continue;
        }
        const lsdtt::patch_record& patch =
            classified.patches[static_cast<std::size_t>(rec.label - 1)];
        metadata_agrees = metadata_agrees && rec.hand_band == patch.band &&
            std::abs(rec.median_height - patch.median_field) < 1e-6;
    }
    check(metadata_agrees,
          "patch and classification records agree on band and median");

    // Window integration regressions found during the Phase 1 review.
    lsdtt::ladder_request windowed;
    windowed.break_algorithm = lsdtt::ladder_break_choice::jenks;
    windowed.level_count = k_st_n_terraces;
    windowed.reach_window_metres = 500.0;
    windowed.vertical_accuracy_metres = 0.5;
    windowed.min_window_pixels = 1000000;  // force global fallback
    std::vector<std::int32_t> partial_source = sh.source_id;
    std::size_t missing_donor_pixel = partial_source.size();
    for (std::size_t i = 0; i < partial_source.size(); ++i)
    {
        if (sh.admitted.get_data()[i] > 0.5f && sh.hand.get_mask()[i])
        {
            partial_source[i] = -1;
            missing_donor_pixel = i;
            break;
        }
    }
    const lsdtt::ladder_result sparse =
        lsdtt::build_reference_field_ladder(
            sh.hand, sh.admitted, windowed, &partial_source, &sh.flow_distance);
    bool sparse_metadata_agrees = !sparse.windows.empty();
    bool windows_sorted = true;
    for (std::size_t i = 0; i < sparse.windows.size(); ++i)
    {
        const lsdtt::ladder_window& w = sparse.windows[i];
        sparse_metadata_agrees = sparse_metadata_agrees &&
            w.status == lsdtt::ladder_status::degenerate &&
            w.used_global_breaks &&
            w.k_selected == static_cast<int>(w.breaks_m.size()) + 1;
        if (i > 0)
        {
            const lsdtt::ladder_window& prev = sparse.windows[i - 1];
            windows_sorted = windows_sorted &&
                (prev.source_id < w.source_id ||
                 (prev.source_id == w.source_id &&
                  prev.window_centre_m <= w.window_centre_m));
        }
    }
    check(sparse_metadata_agrees,
          "sparse windows report the global partition they actually apply");
    check(windows_sorted,
          "ladder windows have deterministic source and centre ordering");
    check(missing_donor_pixel < partial_source.size() &&
          sparse.band_index.get_data()[missing_donor_pixel] >= 0.0f,
          "an admitted pixel without a donor receives the global fallback band");

    lsdtt::ladder_request unresolved = windowed;
    unresolved.level_count = 0;
    unresolved.max_level_count = k_st_n_terraces;
    unresolved.min_window_pixels = 1;
    unresolved.vertical_accuracy_metres = 1000000.0;
    const lsdtt::ladder_result unresolved_result =
        lsdtt::build_reference_field_ladder(
            sh.hand, sh.admitted, unresolved, &sh.source_id, &sh.flow_distance);
    bool unresolved_is_one_band = !unresolved_result.windows.empty();
    for (const lsdtt::ladder_window& w : unresolved_result.windows)
        unresolved_is_one_band = unresolved_is_one_band &&
            w.status == lsdtt::ladder_status::unresolved &&
            w.k_selected == 1 && w.breaks_m.empty();
    check(unresolved_is_one_band && unresolved_result.max_band_index == 0,
          "auto-unresolved windows retain their own honest one-band partition");
}

void run_ladder_checks()
{
    std::cout << "lsdtt-test-fluvial-surfaces\n"
              << "Phase 0c reference-field ladder\n\n";

    // -------------------------------------------------------------------------
    // FINDING, and it decides which algorithm the rest of these tests use.
    //
    // find_population_breaks bins on log10. A terrace flight is spaced LINEARLY
    // in height above the channel, so log binning compresses it: treads at
    // 8 / 16 / 24 / 32 m land in bins 39, 45, 48 and 50 of a 64-bin histogram -
    // spacings of 6, 3 and 2 bins against a smoothing kernel about 5 bins wide.
    // The upper treads merge into one mode and no antimode survives between
    // them, so the finder returns NO breaks at all on a regular staircase.
    //
    // That is not a defect in the finder. It was built to separate a bimodal
    // valley-floor-versus-hillslope population, which is exactly a two-band cut
    // and which it does well. It is a statement about which reference
    // populations it can ladder, and the staircase is not one of them.
    // Fisher-Jenks, working on the linear field, resolves the same flight.
    //
    // Both facts are asserted below so neither can rot, and T2 accordingly runs
    // its resolvability arms under jenks.
    // -------------------------------------------------------------------------
    std::cout << "-- the antimode finder on a linearly spaced flight --\n";

    const lsdtt::lsdtt_terrace_staircase_surface st =
        lsdtt::generate_terrace_staircase_surface(
            k_st_n_rows, k_st_n_cols, k_st_cell_size,
            k_st_n_terraces, k_st_tread_width, k_st_riser, k_st_valley_grad,
            0.0, k_st_cross_grad, 0, 30.0, 0.0, 100.0, 0.5, 1, 1);

    const staircase_hand sh = hand_for_staircase(st);

    lsdtt::ladder_request base_req;
    base_req.reference_field     = lsdtt::reference_field_choice::hand;
    base_req.break_algorithm     = lsdtt::ladder_break_choice::antimode;
    base_req.reach_window_metres = -1.0;   // global
    base_req.hist_bins           = 64;
    base_req.min_field_metres    = 0.05;
    base_req.max_valley_fraction = 0.6;

    {
        const lsdtt::ladder_result anti =
            lsdtt::build_reference_field_ladder(sh.hand, sh.admitted, base_req);

        lsdtt::ladder_request jr = base_req;
        jr.break_algorithm = lsdtt::ladder_break_choice::jenks;
        const lsdtt::ladder_result jen =
            lsdtt::build_reference_field_ladder(sh.hand, sh.admitted, jr);

        std::cout << "      antimode candidates on the staircase = "
                  << anti.global_breaks_m.size()
                  << ", jenks k = " << jen.windows.at(0).k_selected << "\n";

        check(anti.global_breaks_m.empty(),
              "the antimode finder returns no breaks on a linearly spaced "
              "flight (log binning merges the upper treads)");
        check(jen.windows.at(0).k_selected >= 4,
              "jenks on the linear field resolves the same flight");
    }

    // -------------------------------------------------------------------------
    // T5 - level_count = 2 reproduces the legacy valley cut EXACTLY.
    //
    // This is the pin that lets the valley engine be re-plumbed onto the ladder
    // without changing its answer. find_population_breaks is shared with
    // lsdtt-valley-centreline, so the ladder must WRAP it, not replace it.
    //
    // Run on a BIMODAL population, because that is what the legacy cut is for
    // and what it can actually separate (see the finding above): a tight
    // low-lying valley-floor mode and a broad elevated hillslope mode.
    // -------------------------------------------------------------------------
    std::cout << "\n-- T5: k = 2 reproduces the legacy valley cut --\n";
    {
        const int nr = 200, nc = 200;
        const std::size_t np = static_cast<std::size_t>(nr) * nc;
        lsdtt::lsdtt_raster_meta bm = st.dem.get_meta();
        bm.n_rows = nr; bm.n_cols = nc;

        // Two well-separated log-space modes plus a sparse high tail.
        //
        // The tail is not decoration. find_population_breaks scans for local
        // maxima over bins 1..nb-2, so a mode occupying the TOPMOST histogram
        // bin is invisible to it and only one maximum is found - whereupon it
        // returns nothing at all. Real HAND populations always carry a valley-
        // wall tail above the upper mode; a synthetic one has to supply it.
        std::vector<float>   bv(np, 0.0f);
        std::vector<uint8_t> bmk(np, 1);
        std::vector<float>   ba(np, 1.0f);
        std::vector<uint8_t> bam(np, 1);
        for (std::size_t i = 0; i < np; ++i)
        {
            const double u = static_cast<double>((i * 2654435761u) % 1000) / 1000.0;
            if (i % 97 == 0)     bv[i] = static_cast<float>(150.0 + 100.0 * u);
            else if (i % 3 != 0) bv[i] = static_cast<float>(1.5 + 1.0 * u);
            else                 bv[i] = static_cast<float>(30.0 + 15.0 * u);
        }
        const lsdtt::lsdtt_raster_cpu field(
            std::move(bv), std::move(bmk), bm);
        const lsdtt::lsdtt_raster_cpu adm(
            std::move(ba), std::move(bam), bm);

        // The legacy call, verbatim, on the same population.
        const int nbins = 64;
        const double floor_m = 0.05;
        const double lo_log = std::log10(floor_m);
        double hi_log = lo_log;
        const std::vector<float>& fvv = field.get_data();
        for (std::size_t i = 0; i < np; ++i)
            hi_log = std::max(hi_log,
                std::log10(std::max<double>(fvv[i], floor_m)));
        hi_log += 1e-6;

        std::vector<std::int64_t> hist(nbins, 0);
        for (std::size_t i = 0; i < np; ++i)
        {
            const double l = std::log10(std::max<double>(fvv[i], floor_m));
            int b = static_cast<int>((l - lo_log) / (hi_log - lo_log) * nbins);
            ++hist[static_cast<std::size_t>(std::clamp(b, 0, nbins - 1))];
        }
        const std::vector<double> legacy =
            lsdtt::find_population_breaks(hist, lo_log, hi_log, 0.6);

        lsdtt::ladder_request r2 = base_req;
        r2.level_count = 2;
        const lsdtt::ladder_result lr =
            lsdtt::build_reference_field_ladder(field, adm, r2);

        check(!legacy.empty(),
              "the legacy break finder separates a bimodal population");
        check(lr.windows.size() == 1 && lr.windows.at(0).breaks_m.size() == 1,
              "the k = 2 ladder returns exactly one break");
        if (!legacy.empty() && !lr.windows.at(0).breaks_m.empty())
        {
            std::cout << "      legacy break = " << std::setprecision(12)
                      << legacy[0] << " m, ladder break = "
                      << lr.windows.at(0).breaks_m[0] << " m\n";
            check(lr.windows.at(0).breaks_m[0] == legacy[0],
                  "the k = 2 ladder break is BIT-EQUAL to the legacy cut");
        }
    }

    // -------------------------------------------------------------------------
    // T2 - UNRESOLVED is reported, never fabricated.
    //
    // Both arms matter. A test that only checks the unresolvable case can pass
    // by always saying unresolved, so the resolvable case is asserted beside it.
    // -------------------------------------------------------------------------
    std::cout << "\n-- T2: unresolved is reported, not fabricated --\n";

    auto ladder_on = [&](double riser_height, double ev, int level_count)
        -> lsdtt::ladder_result
    {
        const lsdtt::lsdtt_terrace_staircase_surface s =
            lsdtt::generate_terrace_staircase_surface(
                k_st_n_rows, k_st_n_cols, k_st_cell_size,
                k_st_n_terraces, k_st_tread_width, riser_height,
                // A gentle valley so the cross-fall is the steeper gradient at
                // BOTH riser heights: the monotonicity cap forces a small
                // cross-fall when the risers are small, and if the valley
                // gradient then wins, the treads drain downstream and HAND
                // stops being a tread height at all. The two arms must differ
                // only in riser height.
                /*valley_gradient=*/0.0002, 0.0,
                0.6 * riser_height / k_st_tread_width,
                0, 30.0, 0.0, 100.0, 0.5, 1, 1);
        staircase_hand h = hand_for_staircase(s);

        // Admit ONLY tread cells. With the valley wall included the population
        // spans ~180 m and Fisher-Jenks simply spends its bands up the wall,
        // where the gaps are large - so a flight of 0.4 m risers would come
        // back "resolved" on structure that has nothing to do with the
        // terraces. Restricting to the treads makes the test about the flight.
        {
            const std::vector<float>& tid = s.tread_id.get_data();
            std::vector<float>&       adm = h.admitted.get_data();
            for (std::size_t i = 0; i < adm.size(); ++i)
                if (tid[i] <= 0.0f) adm[i] = 0.0f;
        }

        lsdtt::ladder_request r = base_req;
        // Jenks, not antimode: the finding above shows antimode cannot ladder
        // a linearly spaced flight at all, so under it every arm of this test
        // would return k = 1 and the resolvability floor would never be
        // exercised.
        r.break_algorithm          = lsdtt::ladder_break_choice::jenks;
        r.vertical_accuracy_metres = ev;
        r.resolvability_kappa      = 3.0;
        r.level_count              = level_count;
        return lsdtt::build_reference_field_ladder(h.hand, h.admitted, r);
    };

    {
        // Risers of 0.4 m against a required gap of 1.5 m: nothing the DEM can
        // separate.
        const lsdtt::ladder_result lo_r = ladder_on(0.4, 0.5, 0);
        const lsdtt::ladder_window& w = lo_r.windows.at(0);
        std::cout << "      riser 0.4 m, required gap "
                  << w.required_gap_m << " m -> k = " << w.k_selected
                  << ", min centre gap = " << w.min_band_centre_gap_m << " m\n";
        check(w.status == lsdtt::ladder_status::unresolved,
              "a sub-resolution staircase reports UNRESOLVED");
        check(w.k_selected == 1,
              "and emits one band rather than inventing levels");
        check(lo_r.n_windows_unresolved == 1 && lo_r.n_windows_resolved == 0,
              "the unresolved count is reported");

        // The same flight with 8 m risers is comfortably resolvable.
        const lsdtt::ladder_result hi_r = ladder_on(8.0, 0.5, 0);
        const lsdtt::ladder_window& w2 = hi_r.windows.at(0);
        std::cout << "      riser 8.0 m, required gap "
                  << w2.required_gap_m << " m -> k = " << w2.k_selected
                  << ", min centre gap = " << w2.min_band_centre_gap_m << " m\n";
        check(w2.status == lsdtt::ladder_status::resolved,
              "an 8 m staircase resolves, so the test cannot pass by always "
              "saying unresolved");
        check(w2.k_selected >= 2,
              "and recovers more than one level");
        check(w2.min_band_centre_gap_m >= w2.required_gap_m,
              "the selected partition really does clear the floor");

        // An EXPLICIT k that cannot be resolved is honoured and flagged, never
        // silently replaced with a different number.
        const lsdtt::ladder_result forced = ladder_on(0.4, 0.5, 4);
        const lsdtt::ladder_window& w3 = forced.windows.at(0);
        std::cout << "      explicit k = 4 on the 0.4 m flight -> k = "
                  << w3.k_selected << ", status "
                  << static_cast<int>(w3.status) << "\n";
        check(w3.status == lsdtt::ladder_status::unresolved,
              "an inadmissible explicit level_count reports unresolved");
    }

    // -------------------------------------------------------------------------
    // T4 - labels are globally unique and each maps to exactly one band.
    //
    // The band is NOT packed into the label (see the header): it travels in the
    // band_index raster and in patch_record, so what has to hold is that the
    // label -> band map is single-valued and that every label round-trips
    // through the float32 backing store exactly.
    // -------------------------------------------------------------------------
    std::cout << "\n-- T4: labels globally unique, band map single-valued --\n";
    {
        lsdtt::ladder_request r = base_req;
        r.break_algorithm          = lsdtt::ladder_break_choice::jenks;
        r.vertical_accuracy_metres = 0.5;
        const lsdtt::ladder_result lr =
            lsdtt::build_reference_field_ladder(sh.hand, sh.admitted, r);

        const lsdtt::banded_label_result bl =
            lsdtt::label_components_within_bands(
                lr.band_index, sh.admitted, sh.hand,
                /*min_component_pixels=*/50, lsdtt::pixel_connectivity::d4);

        const std::vector<float>& lab = bl.labels.get_data();
        const std::vector<float>& bnd = lr.band_index.get_data();

        std::unordered_map<std::int64_t, int> label_band;
        bool single_valued = true;
        bool exact_float   = true;
        std::int64_t max_label = 0;

        for (std::size_t i = 0; i < lab.size(); ++i)
        {
            if (lab[i] <= 0.0f) continue;
            const std::int64_t l = static_cast<std::int64_t>(lab[i]);
            max_label = std::max(max_label, l);
            if (static_cast<float>(l) != lab[i]) exact_float = false;
            const int b = static_cast<int>(bnd[i]);
            auto it = label_band.find(l);
            if (it == label_band.end()) label_band[l] = b;
            else if (it->second != b) single_valued = false;
        }

        std::cout << "      " << bl.patches.size() << " patches, max label "
                  << max_label << ", " << (lr.max_band_index + 1) << " bands\n";

        check(!bl.patches.empty(), "band-then-label produces patches");
        check(single_valued,
              "every label maps to exactly one band");
        check(exact_float,
              "every label round-trips through the float32 store exactly");
        check(max_label == static_cast<std::int64_t>(bl.patches.size()),
              "labels are a dense global running counter 1..N");

        // The patch table agrees with the raster.
        bool table_ok = true;
        for (std::size_t p = 0; p < bl.patches.size(); ++p)
        {
            const lsdtt::patch_record& rec = bl.patches[p];
            if (rec.label != static_cast<std::int64_t>(p) + 1) table_ok = false;
            auto it = label_band.find(rec.label);
            if (it == label_band.end() || it->second != rec.band) table_ok = false;
            if (rec.n_pixels < 50) table_ok = false;
        }
        check(table_ok,
              "patch_record index, label, band and size all agree with the raster");

        // No patch spans two bands - the property band-then-label exists for.
        std::cout << "      (the 2^24 patch-count guard is not exercised: "
                  << "constructing 16.7M patches needs a ~34M-cell grid)\n";
    }

    // -------------------------------------------------------------------------
    // T7 - selection is a CONSTRAINT with an extremum, not a trade-off.
    //
    // Arm 1 runs under jenks. The antimode finder bins on log10 with a floor,
    // so multiplying the field by 100 re-bins the histogram and moves the
    // candidate breaks for reasons that have nothing to do with the selection
    // rule; jenks works on the linear field and takes k directly, which is what
    // isolates the property under test.
    //
    // Arm 2 is the behaviour that distinguishes a constraint from a trade-off:
    // raising the required gap can only ever reduce the admissible k, and past
    // some point admits nothing. Under an additive score a large gap
    // requirement could be paid for by a good variance fit, and k would not
    // fall monotonically.
    // -------------------------------------------------------------------------
    std::cout << "\n-- T7: selection is a constraint, not a trade-off --\n";
    {
        auto scaled_ladder = [&](double scale) -> int
        {
            const std::vector<float>&   hv = sh.hand.get_data();
            const std::vector<uint8_t>& hm = sh.hand.get_mask();
            std::vector<float>   sv(hv.size(), 0.0f);
            std::vector<uint8_t> sm(hm.begin(), hm.end());
            for (std::size_t i = 0; i < hv.size(); ++i)
                sv[i] = static_cast<float>(hv[i] * scale);
            const lsdtt::lsdtt_raster_cpu scaled(
                std::move(sv), std::move(sm), sh.hand.get_meta());

            lsdtt::ladder_request r = base_req;
            r.break_algorithm          = lsdtt::ladder_break_choice::jenks;
            r.vertical_accuracy_metres = 0.5 * scale;
            r.resolvability_kappa      = 3.0;
            const lsdtt::ladder_result lr =
                lsdtt::build_reference_field_ladder(scaled, sh.admitted, r);
            return lr.windows.at(0).k_selected;
        };

        const int k_1   = scaled_ladder(1.0);
        const int k_100 = scaled_ladder(100.0);
        std::cout << "      jenks k at scale 1 = " << k_1
                  << ", at scale 100 = " << k_100 << "\n";
        check(k_1 == k_100,
              "scaling the field and E_v together leaves k unchanged");

        // Monotone response to a tightening constraint.
        std::vector<int> ks;
        for (const double ev : {0.1, 1.0, 5.0, 20.0, 100.0})
        {
            lsdtt::ladder_request r = base_req;
            r.break_algorithm          = lsdtt::ladder_break_choice::jenks;
            r.vertical_accuracy_metres = ev;
            r.resolvability_kappa      = 3.0;
            const lsdtt::ladder_result lr =
                lsdtt::build_reference_field_ladder(sh.hand, sh.admitted, r);
            ks.push_back(lr.windows.at(0).k_selected);
        }
        std::cout << "      k as E_v rises (0.1, 1, 5, 20, 100 m):";
        for (const int k : ks) std::cout << " " << k;
        std::cout << "\n";

        bool monotone = true;
        for (std::size_t i = 1; i < ks.size(); ++i)
            if (ks[i] > ks[i - 1]) monotone = false;
        check(monotone,
              "k never RISES as the resolvability floor is raised");
        check(ks.front() > ks.back(),
              "and it genuinely falls, so the floor is actually binding");
        check(ks.back() == 1,
              "a floor no partition can clear collapses to one band");
    }
}

// =============================================================================
// T3 - the fitted-surface reference field does not shred a fan
//
// Phase 2 of docs/plan_docs/plan_fluvial_surfaces_cli_split.md.
//
// The claim under test is about WHICH FIELD, not which algorithm: the same
// ladder, the same break finder and the same k selection are run twice, and the
// only thing that changes is the raster being banded. HAND on a fan climbs
// continuously from toe to apex, so its bands are annuli of constant height
// above the channel and they cut clean across the lobe boundaries. The residual
// to a fitted surface has no radial trend left in it, so its bands ARE the
// depositional surfaces.
//
// The comparison uses level_count = 0 (auto) for both arms deliberately. Fixing
// k would let the test pass by having been told the answer; auto makes each
// field argue for its own level count against the same resolvability floor, so
// what is measured is the field's structure rather than the caller's.
// =============================================================================
struct fan_hand
{
    lsdtt::lsdtt_raster_cpu   hand;
    lsdtt::lsdtt_raster_cpu   admitted;
    std::vector<std::int64_t> donor_flat;
};

fan_hand hand_for_fan(const lsdtt::lsdtt_alluvial_fan_surface& fan,
                      int channel_threshold_pixels)
{
    const lsdtt::lsdtt_raster_cpu filled =
        lsdtt::fill_pits(fan.dem, 0.0001f, /*n_threads=*/1);

    lsdtt::lsdtt_flowinfo_options o;
    o.boundary_conditions     = "oooo";
    o.build_topological_order = true;
    o.build_s_vector          = true;
    o.n_threads               = 1;

    const std::unique_ptr<lsdtt::lsdtt_flow_base> f = lsdtt::make_flow_base(
        filled.get_data_array(), filled.get_mask_array(), filled.get_meta(), o);

    // No channel is supplied with the fan, so take one the way a user would:
    // a contributing-pixel threshold on the fan's own routing. The trunk runs
    // down the fan axis, which is exactly the geometry that makes HAND climb
    // with distance from the apex.
    const lsdtt::lsdtt_raster_result npx = f->get_n_contributing_pixels_raster();
    std::vector<std::uint8_t> chan(npx.data.size(), 0);
    for (std::size_t i = 0; i < npx.data.size(); ++i)
    {
        if (npx.mask[i] && npx.data[i] >= static_cast<float>(channel_threshold_pixels))
            chan[i] = 1;
    }

    lsdtt::height_above_channel_donor_result hd =
        f->get_height_above_channel_mask_donor_raster(
            fan.dem.get_data(), fan.dem.get_mask(), chan);

    fan_hand out;
    out.hand       = lsdtt::lsdtt_raster_cpu(std::move(hd.height));
    out.donor_flat = hd.donor_flat;

    // Admit every off-channel cell with a defined HAND. Both arms band exactly
    // this set, so a component-count difference cannot come from one arm having
    // been handed a different domain.
    const std::vector<float>&   hv = out.hand.get_data();
    const std::vector<uint8_t>& hm = out.hand.get_mask();
    std::vector<float>   adm(hv.size(), 0.0f);
    std::vector<uint8_t> adm_mask(hv.size(), 1);
    for (std::size_t i = 0; i < hv.size(); ++i)
    {
        if (hm[i] && std::isfinite(hv[i]) && chan[i] == 0) adm[i] = 1.0f;
    }
    out.admitted = lsdtt::lsdtt_raster_cpu(
        std::move(adm), std::move(adm_mask), out.hand.get_meta());
    return out;
}

// How many connected components does a labelling place on the cells of one
// ground-truth region, and how pure are they?
struct region_measurement
{
    std::int64_t n_components = 0;   // labels touching this region at all
    double       purity       = 0.0; // fraction of the region's labelled cells
                                     // sitting in labels whose majority IS this
                                     // region
};

region_measurement measure_region(
    const lsdtt::lsdtt_raster_cpu& labels,
    const lsdtt::lsdtt_raster_cpu& region_id,
    int                            want_region)
{
    const std::vector<float>& lv = labels.get_data();
    const std::vector<float>& rv = region_id.get_data();

    // per label: cells in the wanted region, and cells in total
    std::map<std::int64_t, std::int64_t> in_region;
    std::map<std::int64_t, std::int64_t> total;
    for (std::size_t i = 0; i < lv.size(); ++i)
    {
        if (!(lv[i] > 0.0f)) continue;
        const std::int64_t lab = static_cast<std::int64_t>(std::llround(lv[i]));
        ++total[lab];
        if (static_cast<int>(std::llround(rv[i])) == want_region) ++in_region[lab];
    }

    region_measurement out;
    std::int64_t n_region_cells = 0;
    std::int64_t n_pure_cells   = 0;
    for (const auto& kv : in_region)
    {
        if (kv.second == 0) continue;
        ++out.n_components;
        n_region_cells += kv.second;
        // "Pure" means the label is mostly this region, so a label that merely
        // clips a corner of it does not count as covering it.
        if (2 * kv.second >= total[kv.first]) n_pure_cells += kv.second;
    }
    out.purity = n_region_cells > 0
        ? static_cast<double>(n_pure_cells) / static_cast<double>(n_region_cells)
        : 0.0;
    return out;
}

void run_fitted_surface_checks()
{
    std::cout << "lsdtt-test-fluvial-surfaces\n"
              << "Phase 2 fitted-surface reference field (T3)\n\n";

    // 200 x 200 at 10 m = a 2 km fan, so HAND climbs tens of metres from toe to
    // apex while the lobes sit a constant 4 m below the older surface.
    const lsdtt::lsdtt_alluvial_fan_surface fan =
        lsdtt::generate_alluvial_fan_surface(
            k_fan_n_rows, k_fan_n_cols, k_fan_cell_size,
            /*apex_x_fraction=*/0.5, /*apex_y_fraction=*/0.95,
            /*apex_elevation=*/500.0, k_fan_gradient, /*concavity=*/0.0,
            k_fan_n_lobes, k_fan_inset, /*lobe_half_angle_deg=*/25.0,
            /*noise=*/0.0, /*seed=*/1, /*n_threads=*/1);

    const fan_hand fh = hand_for_fan(fan, /*channel_threshold_pixels=*/500);

    // The same request for both arms. Only reference_field and the banded
    // raster differ.
    lsdtt::ladder_request lr;
    lr.break_algorithm          = lsdtt::ladder_break_choice::jenks;
    lr.level_count              = 0;      // auto: let each field argue its own k
    lr.max_level_count          = 12;
    lr.vertical_accuracy_metres = 0.5;
    lr.resolvability_kappa      = 3.0;    // required band-centre gap 1.5 m
    lr.reach_window_metres      = -1.0;   // global, so windowing cannot confound
    lr.n_threads                = 1;

    // ---- arm 1: HAND ----------------------------------------------------
    lsdtt::ladder_request hr = lr;
    hr.reference_field = lsdtt::reference_field_choice::hand;

    const lsdtt::ladder_result hand_ladder =
        lsdtt::build_reference_field_ladder(fh.hand, fh.admitted, hr);
    const lsdtt::banded_label_result hand_labels =
        lsdtt::label_components_within_bands(
            hand_ladder.band_index, fh.admitted, fh.hand,
            /*min_component_pixels=*/50, lsdtt::pixel_connectivity::d4);

    const region_measurement hand_old =
        measure_region(hand_labels.labels, fan.lobe_id, 0);

    std::cout << "-- HAND field --\n"
              << "      k = " << hand_ladder.windows.front().k_selected
              << ", bands = " << (hand_ladder.max_band_index + 1)
              << ", patches = " << hand_labels.patches.size() << "\n"
              << "      older surface (lobe_id 0) is cut into "
              << hand_old.n_components << " component(s)\n";

    // ---- arm 2: the residual to a fitted plane ---------------------------
    lsdtt::fitted_surface_request fs;
    fs.algorithm  = lsdtt::fitted_surface_choice::cone;
    fs.n_threads  = 1;
    // The apex is found from the DEM, not supplied: the auto path is what a
    // caller gets by default, so it is the one worth testing.

    const lsdtt::fitted_surface_field field =
        lsdtt::build_fitted_surface_reference_field(fan.dem, fh.admitted, fs);
    check(field.valid, "the cone fit to the fan succeeds");
    if (!field.valid) return;

    std::cout << "\n-- fitted-surface field (cone) --\n"
              << "      apex (" << field.cone_coeffs.apex_x << ", "
              << field.cone_coeffs.apex_y << "), z_apex = "
              << field.cone_coeffs.apex_elevation << " m, gradient = "
              << field.cone_coeffs.gradient << "\n"
              << "      rmse = " << field.rmse << " m over "
              << field.n_pixels_used << " cell(s)\n";

    // The auto apex must land on the real one, or every radius is measured from
    // the wrong origin and the residual carries a spurious radial trend.
    check(std::abs(field.cone_coeffs.gradient + k_fan_gradient) < 0.1 * k_fan_gradient,
          "the fitted cone recovers the fan's radial gradient to within 10%");

    lsdtt::ladder_request fr = lr;
    fr.reference_field = lsdtt::reference_field_choice::fitted_surface;

    const lsdtt::ladder_result fit_ladder =
        lsdtt::build_reference_field_ladder(field.field, fh.admitted, fr);
    const lsdtt::banded_label_result fit_labels =
        lsdtt::label_components_within_bands(
            fit_ladder.band_index, fh.admitted, field.field,
            /*min_component_pixels=*/50, lsdtt::pixel_connectivity::d4);

    const region_measurement fit_old =
        measure_region(fit_labels.labels, fan.lobe_id, 0);

    std::cout << "      k = " << fit_ladder.windows.front().k_selected
              << ", bands = " << (fit_ladder.max_band_index + 1)
              << ", patches = " << fit_labels.patches.size() << "\n"
              << "      older surface (lobe_id 0) is cut into "
              << fit_old.n_components << " component(s)\n";

    for (int l = 1; l <= k_fan_n_lobes; ++l)
    {
        const region_measurement m =
            measure_region(fit_labels.labels, fan.lobe_id, l);
        std::cout << "      lobe " << l << ": " << m.n_components
                  << " component(s), purity " << std::setprecision(4)
                  << m.purity << "\n";
    }

    // ---- the assertions --------------------------------------------------
    std::cout << "\n";
    check(hand_old.n_components > 6,
          "HAND shreds the older fan surface into more than 6 components");
    check(fit_old.n_components <= 2,
          "the fitted-surface field keeps the older fan surface whole "
          "(<= 2 components)");
    check(fit_old.n_components < hand_old.n_components,
          "and it is strictly less fragmented than the HAND arm");

    bool lobes_clean = true;
    for (int l = 1; l <= k_fan_n_lobes; ++l)
    {
        const region_measurement m =
            measure_region(fit_labels.labels, fan.lobe_id, l);
        if (m.n_components != 1 || m.purity < 0.95) lobes_clean = false;
    }
    check(lobes_clean,
          "every inset lobe is a single component with purity >= 0.95");

    // The C6 note is a behaviour, not a comment: a fan has no along-channel
    // coordinate, so the fitted-surface field must be banded globally even when
    // window arguments are supplied.
    {
        std::vector<std::int32_t> src(fh.hand.get_data().size(), 0);
        lsdtt::ladder_request wr = fr;
        wr.reach_window_metres = 0.0;    // auto windows: would window on HAND
        const lsdtt::ladder_result w = lsdtt::build_reference_field_ladder(
            field.field, fh.admitted, wr, &src, &fh.hand);
        check(w.windows.size() == 1,
              "the fitted-surface field is banded globally even when window "
              "arguments are supplied (C6)");
    }

    // The antimode finder bins on log10, so a signed field must be refused
    // rather than silently folded into the lowest bin.
    {
        lsdtt::ladder_request ar = fr;
        ar.break_algorithm = lsdtt::ladder_break_choice::antimode;
        bool threw = false;
        try
        {
            (void)lsdtt::build_reference_field_ladder(
                field.field, fh.admitted, ar);
        }
        catch (const std::invalid_argument&) { threw = true; }
        check(threw,
              "the antimode finder refuses a field containing negative values "
              "instead of banding a destroyed histogram");
    }
}

// =============================================================================
// The trough's upstream-of-head exclusion
//
// Terrain beyond the channel head is still valley ground and stays in the
// valley mask, but it must not be carved into the trough. A broad pocket behind
// the head forms a deep distance-to-wall bowl whose lowest spill is the
// upstream tile edge, so the least-cost walk escapes that way and the centreline
// runs off the top of the DEM instead of down the valley.
//
// The guard is a SIGNED along-channel coordinate: each corridor pixel takes the
// along-channel distance of its nearest channel node plus its own offset from
// that node projected onto the local channel direction. Pixels beyond the head
// project backwards and go negative.
//
// This went inert at phase 5d, when the corridor HAND break search it happened
// to live inside moved to lsdtt-fluvial-surfaces, and it was restored as its own
// step afterwards. The fixtures in the smoke tests all have their channel head
// at a DEM edge, where there is no ground beyond the head and the exclusion
// count is legitimately zero - so none of them can tell a working guard from a
// dead one. This test builds the geometry that can.
// =============================================================================
void run_trough_upstream_checks()
{
    std::cout << "lsdtt-test-fluvial-surfaces\n"
              << "Trough upstream-of-head exclusion\n\n";

    // A straight valley running down the grid: a flat floor incised into
    // shoulders, falling to the south. Simple on purpose - the property under
    // test is geometric, and a realistic surface would only add ways to fail
    // for unrelated reasons.
    constexpr int    n_rows = 120;
    constexpr int    n_cols = 60;
    constexpr double cell   = 10.0;
    constexpr int    axis   = n_cols / 2;
    constexpr int    half_w = 6;      // valley half-width in cells

    lsdtt::lsdtt_raster_meta meta;
    meta.n_rows = n_rows;
    meta.n_cols = n_cols;
    meta.cell_size = cell;
    meta.x_min = 0.0;
    meta.y_min = 0.0;
    meta.no_data_value = -9999.0;
    meta.set_default_synthetic_georeference();

    const std::size_t n = static_cast<std::size_t>(n_rows) * n_cols;
    std::vector<float>   z(n, 0.0f);
    std::vector<uint8_t> zm(n, 1u);
    // Above the head the valley opens into a POCKET - wider floor, and the
    // ground stops falling southward. That is the geometry the guard exists
    // for: a broad flat pocket makes a deeper distance-to-wall bowl than the
    // narrow corridor below it, so the least-cost walk is drawn INTO the pocket
    // and out over the northern tile edge. A uniformly-falling valley of
    // constant width does not provoke it, and a test built on one would pass
    // whether the guard worked or not.
    constexpr int head_row  = n_rows / 3;
    constexpr int pocket_w  = 22;     // pocket half-width in cells
    for (int r = 0; r < n_rows; ++r)
    {
        const bool  in_pocket = r < head_row;
        const int   w         = in_pocket ? pocket_w : half_w;
        // The pocket does not drain: it sits at the head's elevation, so its
        // floor is a genuine flat basin behind the channel.
        const double fall = in_pocket ? (0.01 * head_row * cell)
                                      : (0.01 * r * cell);
        for (int c = 0; c < n_cols; ++c)
        {
            const int    d    = std::abs(c - axis);
            const double wall = (d <= w) ? 0.0 : 20.0 * (d - w);
            z[static_cast<std::size_t>(r) * n_cols + c] =
                static_cast<float>(200.0 - fall + wall);
        }
    }
    const lsdtt::lsdtt_raster_cpu dem(std::move(z), std::move(zm), meta);

    // The valley mask is the whole floor, top to bottom - INCLUDING the stretch
    // above where the channel starts. That is the point: the mask keeps it, the
    // trough must not.
    std::vector<float>   vm(n, 0.0f);
    std::vector<uint8_t> vmm(n, 1u);
    for (int r = 0; r < n_rows; ++r)
    {
        const int w = (r < head_row) ? pocket_w : half_w;
        for (int c = axis - w; c <= axis + w; ++c)
            vm[static_cast<std::size_t>(r) * n_cols + c] = 1.0f;
    }
    const lsdtt::lsdtt_raster_cpu valley_in(
        std::move(vm), std::move(vmm), meta);

    // The channel starts at the pocket's downstream lip. On every smoke fixture
    // the head sits at a DEM edge and this region does not exist at all.
    lsdtt::channel_scope_result scope;
    scope.channel_mask.assign(n, 0u);
    scope.channel_source_id.assign(n, -1);
    std::vector<std::int64_t> nodes;
    for (int r = head_row; r < n_rows; ++r)
    {
        const std::size_t f = static_cast<std::size_t>(r) * n_cols + axis;
        nodes.push_back(static_cast<std::int64_t>(f));
        scope.channel_mask[f] = 1u;
        scope.channel_source_id[f] = 0;
    }
    scope.channel_node_lists.push_back(nodes);
    scope.resolved_sources.push_back(nodes.front());

    lsdtt::valley_centreline_request vreq;
    vreq.supplied_valley_mask = &valley_in;
    vreq.corridor_width_metres = 400.0;
    vreq.nodata_buffer_pixels = 0;
    vreq.min_centreline_nodes = 5;
    vreq.n_threads = 1;

    lsdtt::valley_centreline_result vres =
        lsdtt::extract_valley_centrelines(dem, scope, vreq);

    check(!vres.channels.empty(), "the engine returns a channel record");
    if (vres.channels.empty()) return;

    const lsdtt::valley_channel_centreline& ch = vres.channels.front();
    std::cout << "      valley px = " << ch.n_valley_pixels
              << ", excluded from trough (upstream of head) = "
              << ch.n_upstream_of_head_excluded
              << ", centreline = " << ch.centreline_flats.size() << " node(s)\n";

    // The head is 40 rows up a 13-cell-wide floor, so there are on the order of
    // 500 valley cells beyond it. Requiring a substantial count rather than
    // "> 0" is what makes this fail if the signed coordinate degenerates to a
    // handful of cells around the head itself.
    check(ch.n_upstream_of_head_excluded > 200,
          "valley ground beyond the channel head is excluded from the trough");

    // ... and the exclusion is the TROUGH's, not the mask's. If this ever fails
    // the guard has started deleting valley rather than protecting the carve.
    const std::vector<float>& out_vm = vres.valley_mask.get_data();
    std::int64_t n_above_head_in_mask = 0;
    for (int r = 0; r < head_row - 2; ++r)
        for (int c = axis - pocket_w; c <= axis + pocket_w; ++c)
            if (out_vm[static_cast<std::size_t>(r) * n_cols + c] > 0.0f)
                ++n_above_head_in_mask;
    std::cout << "      valley cells above the head retained in the MASK = "
              << n_above_head_in_mask << "\n";
    check(n_above_head_in_mask > 200,
          "and those cells are still in the valley mask");

    check(ch.routing_succeeded && ch.centreline_flats.size() >= 5,
          "the centreline still routes");

    // A cheap invariant, NOT a demonstration. The centreline should never sit
    // north of the head, and this catches it if it ever does.
    //
    // Stated plainly: it does NOT discriminate here. Disabling the guard leaves
    // this assertion passing, because the walk starts at the head node and
    // proceeds downstream, so it cannot reach the pocket on this geometry. The
    // downstream escape the guard was written for needs the filled trough to
    // spill north, which this fixture does not reproduce - the exclusion count
    // above is what actually fails when the guard is broken (verified by
    // disabling it: 1571 -> 0). If someone later builds a fixture where the
    // routed path really does escape, this assertion is where it belongs.
    int min_row = n_rows;
    for (const std::int64_t f : ch.centreline_flats)
        min_row = std::min(min_row, static_cast<int>(f / n_cols));
    std::cout << "      northernmost centreline row = " << min_row
              << " (channel head at row " << head_row << ")\n";
    check(min_row >= head_row - 2,
          "and it does not escape north past the channel head");
}

}  // namespace


// =============================================================================
// Perona-Malik conditioning (Phase 7b)
//
// THE CLAIM the filter is adopted on: it removes within-facet noise while
// leaving the riser standing. Both halves have to be tested, because either one
// alone is trivially satisfiable - a no-op preserves every riser, and a heavy
// blur removes every bit of noise.
//
// So the fixture is a two-tread staircase: a lower tread at 100 m, a one-cell
// riser, an upper tread at 130 m, with +/- 2 m of deterministic noise on both
// treads. The mean filter is computed alongside as the CONTROL, because the
// argument for anisotropic diffusion is not that it smooths - anything smooths -
// but that it does so without moving the break. If a future change quietly
// turned the conductance off, P-M and the mean filter would agree and the
// step-preservation check would catch it.
// =============================================================================
void run_perona_malik_checks()
{
    std::cout << "\n--- Perona-Malik conditioning (7b) ---\n";

    const int n_rows = 80, n_cols = 80;
    const double cell = 10.0;
    const int riser_col = 40;
    const double lower = 100.0, upper = 130.0;

    lsdtt::lsdtt_raster_meta meta;
    meta.n_rows = n_rows;
    meta.n_cols = n_cols;
    meta.cell_size = cell;
    meta.x_min = 0.0;
    meta.y_min = 0.0;
    meta.no_data_value = -9999.0;

    std::vector<float>   d(static_cast<std::size_t>(n_rows) * n_cols, 0.0f);
    std::vector<uint8_t> m(d.size(), uint8_t{1});

    // Deterministic zig-zag noise: no RNG, so the numbers below are pinned.
    for (int r = 0; r < n_rows; ++r)
    {
        for (int c = 0; c < n_cols; ++c)
        {
            const std::size_t i = static_cast<std::size_t>(r) * n_cols + c;
            const double base = (c < riser_col) ? lower : upper;
            const double noise = 2.0 * (((r + c) % 2 == 0) ? 1.0 : -1.0);
            d[i] = static_cast<float>(base + noise);
        }
    }
    const lsdtt::lsdtt_raster_cpu dem(std::move(d), std::move(m), meta);

    lsdtt::perona_malik_options opt;
    opt.n_iterations        = 10;
    opt.lambda              = 0.20;
    opt.conductance         = lsdtt::pm_conductance_choice::exponential;
    opt.k_gradient_quantile = 0.90;

    const lsdtt::perona_malik_result pmr =
        lsdtt::perona_malik_diffusion(dem, opt);
    const std::vector<float>& sd = pmr.smoothed.get_data();
    const std::vector<float>& od = dem.get_data();

    check(pmr.k_metres_per_cell > 0.0,
          "PM1: K was estimated from the DEM's own gradients (K = "
          + std::to_string(pmr.k_metres_per_cell) + " m/cell)");

    // 3x3 mean filter as the control.
    std::vector<float> mean(od.size(), 0.0f);
    for (int r = 0; r < n_rows; ++r)
    {
        for (int c = 0; c < n_cols; ++c)
        {
            double acc = 0.0; int cnt = 0;
            for (int dr = -1; dr <= 1; ++dr)
                for (int dc = -1; dc <= 1; ++dc)
                {
                    const int rr = r + dr, cc = c + dc;
                    if (rr < 0 || rr >= n_rows || cc < 0 || cc >= n_cols) continue;
                    acc += od[static_cast<std::size_t>(rr) * n_cols + cc];
                    ++cnt;
                }
            mean[static_cast<std::size_t>(r) * n_cols + c] =
                static_cast<float>(acc / cnt);
        }
    }

    // Tread roughness: mean |z - tread base|, well away from the riser.
    auto tread_roughness = [&](const std::vector<float>& v) {
        double acc = 0.0; int cnt = 0;
        for (int r = 2; r < n_rows - 2; ++r)
            for (int c = 5; c < riser_col - 5; ++c)
            {
                acc += std::fabs(v[static_cast<std::size_t>(r) * n_cols + c] - lower);
                ++cnt;
            }
        return acc / cnt;
    };

    // Step height measured across the riser, one cell either side.
    auto step_height = [&](const std::vector<float>& v) {
        double acc = 0.0; int cnt = 0;
        for (int r = 5; r < n_rows - 5; ++r)
        {
            const std::size_t l = static_cast<std::size_t>(r) * n_cols + (riser_col - 1);
            const std::size_t h = static_cast<std::size_t>(r) * n_cols + riser_col;
            acc += static_cast<double>(v[h]) - v[l];
            ++cnt;
        }
        return acc / cnt;
    };

    const double rough_raw  = tread_roughness(od);
    const double rough_pm   = tread_roughness(sd);
    const double rough_mean = tread_roughness(mean);
    const double step_raw   = step_height(od);
    const double step_pm    = step_height(sd);
    const double step_mean  = step_height(mean);

    std::cout << "  tread roughness  raw " << rough_raw
              << "  PM " << rough_pm << "  mean-filter " << rough_mean << "\n";
    std::cout << "  step height      raw " << step_raw
              << "  PM " << step_pm << "  mean-filter " << step_mean << "\n";

    check(rough_pm < 0.2 * rough_raw,
          "PM2: within-tread roughness falls to under a fifth of the raw DEM");
    check(step_pm > 0.95 * step_raw,
          "PM3: the riser keeps more than 95% of its height");
    check(step_pm > step_mean + 1.0,
          "PM4: the riser survives P-M by more than a metre over the mean "
          "filter, so the conductance is doing something");

    // Nodata must not be diffused into, and must not leak a value outward.
    {
        std::vector<float>   d2 = dem.get_data();
        std::vector<uint8_t> m2 = dem.get_mask();
        const std::size_t hole = static_cast<std::size_t>(40) * n_cols + 20;
        d2[hole] = static_cast<float>(meta.no_data_value);
        m2[hole] = uint8_t{0};
        const lsdtt::lsdtt_raster_cpu holed(std::move(d2), std::move(m2), meta);
        const lsdtt::perona_malik_result hr =
            lsdtt::perona_malik_diffusion(holed, opt);
        check(hr.smoothed.get_mask()[hole] == 0,
              "PM5: a nodata cell stays nodata");
        const double neigh = hr.smoothed.get_data()[hole - 1];
        check(neigh > lower - 5.0 && neigh < lower + 5.0,
              "PM6: the cell beside the hole is not pulled toward nodata");
    }

    // Stability guard is a refusal, not a warning.
    {
        bool threw = false;
        lsdtt::perona_malik_options bad = opt;
        bad.lambda = 0.30;
        try { lsdtt::perona_malik_diffusion(dem, bad); }
        catch (const std::exception&) { threw = true; }
        check(threw, "PM7: lambda above the 0.25 stability limit is rejected");
    }
}


// =============================================================================
// Coherence-enhancing diffusion (Phase 7c)
//
// WHY A SECOND FILTER, AND WHAT SEPARATES IT FROM THE FIRST. Perona-Malik keys
// on gradient MAGNITUDE, so it keeps a riser because the riser is steep - and
// keeps a gully cut through that riser for exactly the same reason. Weickert
// keys on the structure tensor: it smooths ALONG the dominant local orientation
// and barely at all across it, so it keeps the riser because the riser is
// LINEAR, and it will close a break in the riser rather than preserve it.
//
// That difference is the whole justification for shipping both, so the fixture
// is a riser BROKEN BY NOTCHES: a two-tread staircase whose riser is cut back
// to the lower tread every eighth row. The discriminating statistic is the
// MINIMUM step height along the riser. Both filters should leave the intact
// riser standing; only the orientation-aware one should raise the floor in the
// notches. If a future change made the two behave alike, CED4 fails.
// =============================================================================
void run_coherence_diffusion_checks()
{
    std::cout << "\n--- Coherence-enhancing diffusion (7c) ---\n";

    const int n_rows = 80, n_cols = 80;
    const double cell = 10.0;
    const int riser_col = 40;
    const double lower = 100.0, upper = 130.0;

    lsdtt::lsdtt_raster_meta meta;
    meta.n_rows = n_rows;
    meta.n_cols = n_cols;
    meta.cell_size = cell;
    meta.x_min = 0.0;
    meta.y_min = 0.0;
    meta.no_data_value = -9999.0;

    std::vector<float>   d(static_cast<std::size_t>(n_rows) * n_cols, 0.0f);
    std::vector<uint8_t> m(d.size(), uint8_t{1});
    std::vector<int>     notched(n_rows, 0);

    for (int r = 0; r < n_rows; ++r)
    {
        // Every 8th row the upper tread is cut back past the riser, so the
        // riser line is interrupted. This is what gullying does to a flight.
        notched[r] = (r % 8 == 0) ? 1 : 0;
        for (int c = 0; c < n_cols; ++c)
        {
            const std::size_t i = static_cast<std::size_t>(r) * n_cols + c;
            double z = (c < riser_col) ? lower : upper;
            if (notched[r] && c >= riser_col && c < riser_col + 6) z = lower;
            const double noise = 1.0 * (((r + c) % 2 == 0) ? 1.0 : -1.0);
            d[i] = static_cast<float>(z + noise);
        }
    }
    const lsdtt::lsdtt_raster_cpu dem(std::move(d), std::move(m), meta);

    // Step height across the riser at one row. Averaged over four cells each
    // side, not read off single cells: the checkerboard noise alternates with
    // column parity, so a two-cell reading returns the noise amplitude rather
    // than the step and a notched row scores 2 m instead of 0.
    auto step_at = [&](const std::vector<float>& v, int r) {
        double lo_sum = 0.0, hi_sum = 0.0;
        for (int k = 1; k <= 4; ++k)
        {
            lo_sum += v[static_cast<std::size_t>(r) * n_cols + (riser_col - k)];
            hi_sum += v[static_cast<std::size_t>(r) * n_cols + (riser_col + k)];
        }
        return (hi_sum - lo_sum) / 4.0;
    };
    auto step_stats = [&](const std::vector<float>& v, bool want_notched) {
        double acc = 0.0; int cnt = 0;
        for (int r = 5; r < n_rows - 5; ++r)
        {
            if ((notched[r] != 0) != want_notched) continue;
            acc += step_at(v, r); ++cnt;
        }
        return cnt ? acc / cnt : 0.0;
    };

    lsdtt::perona_malik_options pmo;
    pmo.n_iterations        = 20;
    pmo.lambda              = 0.20;
    pmo.k_gradient_quantile = 0.50;
    // Hold the result: binding a reference straight to .smoothed.get_data() on
    // the temporary dangles, and the freed buffer then reads back as whatever
    // is written into it next - which made both filters report identical
    // numbers and very nearly passed as agreement.
    const lsdtt::perona_malik_result pm_res =
        lsdtt::perona_malik_diffusion(dem, pmo);
    const std::vector<float>& pm = pm_res.smoothed.get_data();

    const lsdtt::lsdtt_raster_cpu ced_r =
        lsdtt::coherence_enhancing_diffusion(
            dem, /*n_iterations=*/20, /*time_step=*/0.15,
            /*derivative_radius_metres=*/cell,
            /*integration_radius_metres=*/3.0 * cell,
            /*contrast=*/0.2, /*min_diffusivity=*/0.001, /*n_threads=*/1);
    const std::vector<float>& ce = ced_r.get_data();
    const std::vector<float>& raw = dem.get_data();

    const double intact_raw = step_stats(raw, false);
    const double intact_pm  = step_stats(pm,  false);
    const double intact_ce  = step_stats(ce,  false);
    const double notch_raw  = step_stats(raw, true);
    const double notch_pm   = step_stats(pm,  true);
    const double notch_ce   = step_stats(ce,  true);

    std::cout << "  step at intact rows   raw " << intact_raw
              << "  PM " << intact_pm << "  CED " << intact_ce << "\n";
    std::cout << "  step at notched rows  raw " << notch_raw
              << "  PM " << notch_pm << "  CED " << notch_ce << "\n";

    check(intact_ce > 0.90 * intact_raw,
          "CED1: the intact riser keeps more than 90% of its height");
    check(notch_raw < 1.0,
          "CED2: the fixture's notches really do cut the riser out (raw step "
          "there = " + std::to_string(notch_raw) + " m)");
    check(notch_ce > notch_pm,
          "CED3: the orientation-aware filter closes the notch further than "
          "the gradient-magnitude one");
    check(notch_ce > notch_raw + 1.0,
          "CED4: the notch is bridged by more than a metre, so smoothing ALONG "
          "the riser is doing something a magnitude filter does not");

    // Nodata and the stability guard, as for Perona-Malik.
    {
        std::vector<float>   d2 = dem.get_data();
        std::vector<uint8_t> m2 = dem.get_mask();
        const std::size_t hole = static_cast<std::size_t>(40) * n_cols + 20;
        d2[hole] = static_cast<float>(meta.no_data_value);
        m2[hole] = uint8_t{0};
        const lsdtt::lsdtt_raster_cpu holed(std::move(d2), std::move(m2), meta);
        const lsdtt::lsdtt_raster_cpu hr =
            lsdtt::coherence_enhancing_diffusion(
                holed, 5, 0.15, cell, 3.0 * cell, 0.2, 0.001, 1);
        check(hr.get_mask()[hole] == 0, "CED5: a nodata cell stays nodata");
    }
    {
        bool threw = false;
        try
        {
            lsdtt::coherence_enhancing_diffusion(
                dem, 5, 0.30, cell, 3.0 * cell, 0.2, 0.001, 1);
        }
        catch (const std::exception&) { threw = true; }
        check(threw, "CED6: a time step above the 0.25 stability limit is "
                     "rejected");
    }
}

int main(int argc, char* argv[])
{
    try
    {
        if (argc == 2 && std::string(argv[1]) == "--band-then-label")
        {
            run_band_then_label_checks();
            std::cout << "\n" << (g_n_checks - g_n_failed) << " / " << g_n_checks
                      << " checks passed\n"
                      << "RESULT: " << (g_n_failed == 0 ? "PASS" : "FAIL") << "\n";
            return g_n_failed == 0 ? 0 : 1;
        }

        if (argc == 2 && std::string(argv[1]) == "--fitted-surface")
        {
            run_fitted_surface_checks();
            std::cout << "\n" << (g_n_checks - g_n_failed) << " / " << g_n_checks
                      << " checks passed\n"
                      << "RESULT: " << (g_n_failed == 0 ? "PASS" : "FAIL") << "\n";
            return g_n_failed == 0 ? 0 : 1;
        }

        if (argc == 2 && std::string(argv[1]) == "--trough-upstream")
        {
            run_trough_upstream_checks();
            std::cout << "\n" << (g_n_checks - g_n_failed) << " / " << g_n_checks
                      << " checks passed\n"
                      << "RESULT: " << (g_n_failed == 0 ? "PASS" : "FAIL") << "\n";
            return g_n_failed == 0 ? 0 : 1;
        }

        if (argc == 2 && std::string(argv[1]) == "--coherence-diffusion")
        {
            run_coherence_diffusion_checks();
            std::cout << "\n" << (g_n_checks - g_n_failed) << " / " << g_n_checks
                      << " checks passed\n"
                      << "RESULT: " << (g_n_failed == 0 ? "PASS" : "FAIL") << "\n";
            return g_n_failed == 0 ? 0 : 1;
        }

        if (argc == 2 && std::string(argv[1]) == "--perona-malik")
        {
            run_perona_malik_checks();
            std::cout << "\n" << (g_n_checks - g_n_failed) << " / " << g_n_checks
                      << " checks passed\n"
                      << "RESULT: " << (g_n_failed == 0 ? "PASS" : "FAIL") << "\n";
            return g_n_failed == 0 ? 0 : 1;
        }

        if (argc == 2 && std::string(argv[1]) == "--ladder")
        {
            run_ladder_checks();
            std::cout << "\n" << (g_n_checks - g_n_failed) << " / " << g_n_checks
                      << " checks passed\n"
                      << "RESULT: " << (g_n_failed == 0 ? "PASS" : "FAIL") << "\n";
            return g_n_failed == 0 ? 0 : 1;
        }

        if (argc == 3 && std::string(argv[1]) == "--write-cli-fixture")
        {
            write_cli_fixture(argv[2]);
            write_phase0b_fixtures(argv[2]);
            return 0;
        }

        std::cout << "lsdtt-test-fluvial-surfaces\n"
                  << "Phase 0a/0b fixture self-checks\n\n";
        run_self_checks();
        run_terrace_staircase_checks();
        run_alluvial_fan_checks();
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
