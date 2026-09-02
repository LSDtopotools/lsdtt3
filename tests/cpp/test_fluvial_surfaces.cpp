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

}  // namespace

int main(int argc, char* argv[])
{
    try
    {
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
