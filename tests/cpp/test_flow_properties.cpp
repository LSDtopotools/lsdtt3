// test_flow_properties.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
//
// -----------------------------------------------------------------------------
// Part of the flow-routing benchmark (docs/testing_and_benchmarking.md): property tests for depression
// handling and flow routing.
//
// These are *property* tests, not comparisons: each check asserts an invariant
// that must hold for any correct implementation, so nothing here needs RichDEM,
// network access, or a downloaded DEM.  All fixtures are deterministic
// synthetic surfaces generated in-process.  That makes this driver part of the
// default CTest run; the RichDEM comparisons of Phases 3-4 sit behind the
// opt-in `benchmark` label.
//
// Properties checked:
//   P1  filled DEM contains no interior local minimum (depression-free)
//   P2  fill never lowers the surface            (filled  >= original)
//   P3  breach never raises the surface          (breached <= original)
//   P4  total D8 accumulation reaching outlets equals the valid-cell area
//   P5  a receiver's accumulation is >= each donor's
//   P6  filling an already-filled DEM is a no-op (idempotence)
//   P7  D-inf / MFD receiver weights sum to 1.0 per cell that has receivers
//   P8  D8 accumulation is exact on a planar ramp (analytic oracle; no reference)
//
// Section 8 additionally *characterises* (reports, does not assert) the
// alternative fill algorithm exposed by lsdtt-raster-preprocessing, per Simon's
// request: relative speed and how far its result diverges from Barnes 2014.
//
// Usage: lsdtt-test-flow-properties [output_dir]
//   output_dir (default ".") receives fill_algorithm_characterisation.csv.
// Prints one PASS/FAIL line per check and a final "RESULT: PASS|FAIL".
// -----------------------------------------------------------------------------

#include <chrono>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "lsdtt_raster_cpu.hpp"
#include "lsdtt_raster_flow_preprocessing_cpu.hpp"
#include "lsdtt_raster_synthetic_surfaces_cpu.hpp"
#include "lsdtt_flowinfo_cpu.hpp"
#include "lsdtt_flowgraph_cpu.hpp"

namespace
{

// =============================================================================
// Tiny check harness (same shape as tests/cpp/test_cosmo_physics.cpp)
// =============================================================================
int g_n_checks = 0;
int g_n_failed = 0;

void check(const std::string& name, bool ok, const std::string& detail = "")
{
    ++g_n_checks;
    if (!ok) ++g_n_failed;
    std::cout << (ok ? "PASS" : "FAIL") << "  " << name;
    if (!ok && !detail.empty()) std::cout << "   [" << detail << "]";
    std::cout << "\n";
}


// =============================================================================
// Fixtures
//
// Deterministic synthetic surfaces.  Two relief regimes are covered because
// they stress depression handling very differently: a high-relief fractal
// behaves like the Alpine benchmark tile, while a low-relief surface produces
// many shallow, broad depressions and large flats -- the case where filling
// degenerates and epsilon gradients matter most.
// =============================================================================
struct fixture
{
    std::string      name;
    lsdtt::lsdtt_raster_cpu dem;
};


lsdtt::lsdtt_raster_meta make_meta(int n_rows, int n_cols, double cell_size)
{
    lsdtt::lsdtt_raster_meta meta;
    meta.n_rows        = n_rows;
    meta.n_cols        = n_cols;
    meta.cell_size     = cell_size;
    meta.no_data_value = -9999.0;
    meta.set_default_synthetic_georeference();
    return meta;
}


// Hand-built cone with a carved interior pit.  Small enough that the expected
// behaviour is obvious by inspection, which makes it the fixture that localises
// a failure when the fractal surfaces also fail.
lsdtt::lsdtt_raster_cpu make_cone_with_pit(int n = 41, double cell_size = 10.0)
{
    lsdtt::lsdtt_raster_cpu dem;
    lsdtt::lsdtt_raster_meta meta = make_meta(n, n, cell_size);

    lsdtt::lsdtt_array2d_f  data(n, n, 0.0f);
    lsdtt::lsdtt_array2d_u8 mask(n, n, 1);

    const double centre = 0.5 * (n - 1);
    for (int r = 0; r < n; ++r)
    {
        for (int c = 0; c < n; ++c)
        {
            const double dr = r - centre;
            const double dc = c - centre;
            // Cone falling away from the centre: unambiguous radial drainage.
            data(r, c) = static_cast<float>(100.0 - 1.5 * std::sqrt(dr * dr + dc * dc));
        }
    }

    // Carve two interior pits well away from the edges.
    for (int r = 10; r <= 13; ++r)
        for (int c = 10; c <= 13; ++c)
            data(r, c) -= 12.0f;

    for (int r = 26; r <= 30; ++r)
        for (int c = 24; c <= 28; ++c)
            data(r, c) -= 7.0f;

    dem.set_from_vectors(std::move(data.get_vector()),
                         std::move(mask.get_vector()),
                         meta);
    return dem;
}


// Punch an interior nodata hole so the nodata-handling paths are exercised.
lsdtt::lsdtt_raster_cpu with_interior_nodata(const lsdtt::lsdtt_raster_cpu& in)
{
    lsdtt::lsdtt_raster_cpu out = in;
    const lsdtt::lsdtt_raster_meta& meta = out.meta();

    auto& data = out.get_data_array();
    auto& mask = const_cast<lsdtt::lsdtt_array2d_u8&>(out.get_mask_array());

    const int r0 = meta.n_rows / 3;
    const int c0 = meta.n_cols / 3;
    for (int r = r0; r < r0 + 5; ++r)
    {
        for (int c = c0; c < c0 + 5; ++c)
        {
            data(r, c) = static_cast<float>(meta.no_data_value);
            mask(r, c) = 0;
        }
    }
    return out;
}


std::vector<fixture> build_fixtures()
{
    std::vector<fixture> f;

    f.push_back({"cone_with_pits", make_cone_with_pit()});

    // High relief: stands in for the Alpine benchmark tile.
    f.push_back({"fractal_high_relief",
                 lsdtt::generate_diamond_square_surface(
                     make_meta(129, 129, 30.0),
                     /*seed=*/17, /*roughness=*/0.55, /*initial_scale=*/1.0,
                     /*target_mean=*/1200.0, /*target_stddev=*/300.0)});

    // Low relief: broad shallow depressions and near-flats.  Standard deviation
    // is small relative to the 30 m cell size, so gradients are gentle
    // everywhere and filled flats dominate.
    f.push_back({"fractal_low_relief",
                 lsdtt::generate_power_law_fractal_surface_fft(
                     make_meta(129, 129, 30.0),
                     /*seed=*/23, /*beta=*/2.0,
                     /*target_mean=*/50.0, /*target_stddev=*/1.5)});

    f.push_back({"fractal_high_relief_nodata",
                 with_interior_nodata(f[1].dem)});

    return f;
}


// =============================================================================
// Helpers
// =============================================================================

// True if (r,c) is an interior cell whose full 8-neighbourhood is in bounds and
// valid.  Edge cells and cells touching nodata are legitimate drainage exits
// under open boundary conditions, so they can never be "spurious" minima.
bool is_interior_valid(const lsdtt::lsdtt_raster_cpu& dem, int r, int c)
{
    const lsdtt::lsdtt_raster_meta& meta = dem.meta();
    const auto& mask = dem.get_mask_array();

    if (r <= 0 || c <= 0 || r >= meta.n_rows - 1 || c >= meta.n_cols - 1) return false;
    if (mask(r, c) == 0) return false;

    for (int dr = -1; dr <= 1; ++dr)
        for (int dc = -1; dc <= 1; ++dc)
            if (mask(r + dr, c + dc) == 0) return false;

    return true;
}


// Counts interior cells strictly lower than all eight neighbours: the raster
// level definition of a depression, computed without consulting the flow
// router.  Keeping this independent of lsdtt_flowinfo matters -- a shared bug
// between filling and routing would otherwise cancel out and hide itself.
std::int64_t count_interior_minima(const lsdtt::lsdtt_raster_cpu& dem)
{
    const lsdtt::lsdtt_raster_meta& meta = dem.meta();
    const auto& data = dem.get_data_array();

    std::int64_t n_minima = 0;
    for (int r = 1; r < meta.n_rows - 1; ++r)
    {
        for (int c = 1; c < meta.n_cols - 1; ++c)
        {
            if (!is_interior_valid(dem, r, c)) continue;

            const float z = data(r, c);
            bool has_non_higher_neighbour = false;

            for (int dr = -1; dr <= 1 && !has_non_higher_neighbour; ++dr)
            {
                for (int dc = -1; dc <= 1; ++dc)
                {
                    if (dr == 0 && dc == 0) continue;
                    if (data(r + dr, c + dc) <= z)
                    {
                        has_non_higher_neighbour = true;
                        break;
                    }
                }
            }

            if (!has_non_higher_neighbour) ++n_minima;
        }
    }
    return n_minima;
}


struct diff_stats
{
    double       max_abs        = 0.0;
    double       mean_abs       = 0.0;
    std::int64_t n_differing    = 0;
    std::int64_t n_compared     = 0;
    double       signed_min     = 0.0;   // most negative (a - b)
    double       signed_max     = 0.0;   // most positive (a - b)
};


diff_stats compare(const lsdtt::lsdtt_raster_cpu& a, const lsdtt::lsdtt_raster_cpu& b, double tol)
{
    diff_stats s;
    const lsdtt::lsdtt_raster_meta& meta = a.meta();
    const auto& da = a.get_data_array();
    const auto& db = b.get_data_array();
    const auto& ma = a.get_mask_array();
    const auto& mb = b.get_mask_array();

    double sum_abs = 0.0;
    for (int r = 0; r < meta.n_rows; ++r)
    {
        for (int c = 0; c < meta.n_cols; ++c)
        {
            if (ma(r, c) == 0 || mb(r, c) == 0) continue;

            const double d = static_cast<double>(da(r, c)) - static_cast<double>(db(r, c));
            ++s.n_compared;
            sum_abs += std::abs(d);
            if (std::abs(d) > s.max_abs) s.max_abs = std::abs(d);
            if (std::abs(d) > tol)       ++s.n_differing;
            if (d < s.signed_min)        s.signed_min = d;
            if (d > s.signed_max)        s.signed_max = d;
        }
    }
    s.mean_abs = (s.n_compared > 0) ? sum_abs / static_cast<double>(s.n_compared) : 0.0;
    return s;
}


double seconds_since(std::chrono::steady_clock::time_point t0)
{
    const auto t1 = std::chrono::steady_clock::now();
    return std::chrono::duration<double>(t1 - t0).count();
}


// =============================================================================
// P1 / P2 / P6 -- fill properties
// =============================================================================
void test_fill_properties(const fixture& fx)
{
    // Two min_slope regimes.  0.0 is the pure Barnes 2014 result and is what
    // the Phase 3 tight RichDEM comparison will use; 0.0001 is the CLI default
    // and imposes an epsilon gradient across filled flats.  Both must be
    // depression-free -- that is the whole point of the routine.
    for (const float min_slope : {0.0f, 0.0001f})
    {
        const std::string tag =
            fx.name + " min_slope=" + (min_slope == 0.0f ? "0" : "1e-4");

        const lsdtt::lsdtt_raster_cpu filled = lsdtt::fill_pits(fx.dem, min_slope);

        // -- P1: depression-free -------------------------------------------
        const std::int64_t n_minima_before = count_interior_minima(fx.dem);
        const std::int64_t n_minima_after  = count_interior_minima(filled);

        check("P1 depression-free: " + tag,
              n_minima_after == 0,
              "interior minima before=" + std::to_string(n_minima_before) +
              " after=" + std::to_string(n_minima_after));

        // A fixture that had no depressions to begin with would make P1
        // vacuous, so assert the fixture is actually doing work.
        check("P1 fixture is non-trivial: " + tag,
              n_minima_before > 0,
              "fixture contains no depressions; P1 proves nothing");

        // -- P2: fill never lowers the surface ------------------------------
        const diff_stats d = compare(filled, fx.dem, 0.0);
        check("P2 fill is non-lowering: " + tag,
              d.signed_min >= -1e-6,
              "most negative (filled - original) = " + std::to_string(d.signed_min));

        // -- P6: idempotence -------------------------------------------------
        const lsdtt::lsdtt_raster_cpu refilled = lsdtt::fill_pits(filled, min_slope);
        const diff_stats d2 = compare(refilled, filled, 1e-6);
        check("P6 fill is idempotent: " + tag,
              d2.max_abs <= 1e-4,
              "max |refill - fill| = " + std::to_string(d2.max_abs));
    }
}


// =============================================================================
// P3 -- breach properties
// =============================================================================
void test_breach_properties(const fixture& fx)
{
    {
        // NOTE: lsdtt3's breach_depressions is NOT a pure breach.  Before
        // carving, init_breach_arrays_owned() (lsdtt_raster_flow_preprocessing_cpu.cpp
        // ~line 655) raises every pit/flat cell to its lowest neighbour:
        //
        //     if (ba.carved(r, c) <= lowest) { ba.carved(r, c) = lowest; ... }
        //
        // so a naive "breaching never raises the surface" assertion fails by
        // design -- measured up to ~80 m on the high-relief fractal.  The
        // invariant that DOES hold, and that still catches real bugs, is the
        // bounded form: only pit/flat cells may rise, and never above the
        // elevation of their lowest original neighbour.
        const lsdtt::lsdtt_raster_cpu breached = lsdtt::breach_depressions(fx.dem);

        const lsdtt::lsdtt_raster_meta& meta = fx.dem.meta();
        const auto& z0 = fx.dem.get_data_array();
        const auto& z1 = breached.get_data_array();

        std::int64_t n_non_pit_raised = 0;   // a non-pit cell rose: real bug
        std::int64_t n_over_lowest    = 0;   // rose above lowest neighbour: real bug
        std::int64_t n_pits_raised    = 0;   // expected, reported for context
        double       worst_over       = 0.0;

        for (int r = 0; r < meta.n_rows; ++r)
        {
            for (int c = 0; c < meta.n_cols; ++c)
            {
                if (!is_interior_valid(fx.dem, r, c)) continue;

                const double rise = static_cast<double>(z1(r, c)) - static_cast<double>(z0(r, c));
                if (rise <= 1e-6) continue;

                // Lowest neighbour in the ORIGINAL DEM.
                float lowest = std::numeric_limits<float>::max();
                for (int dr = -1; dr <= 1; ++dr)
                    for (int dc = -1; dc <= 1; ++dc)
                        if (dr != 0 || dc != 0)
                            lowest = std::min(lowest, z0(r + dr, c + dc));

                const bool was_pit = (z0(r, c) <= lowest);
                if (!was_pit) ++n_non_pit_raised;
                else          ++n_pits_raised;

                const double over = static_cast<double>(z1(r, c)) - static_cast<double>(lowest);
                if (over > 1e-4)
                {
                    ++n_over_lowest;
                    if (over > worst_over) worst_over = over;
                }
            }
        }

        check("P3 breach (lindsay2016) raises only pit/flat cells: " + fx.name,
              n_non_pit_raised == 0,
              "non-pit cells raised = " + std::to_string(n_non_pit_raised) +
              " (pit cells raised = " + std::to_string(n_pits_raised) + ")");

        check("P3 breach (lindsay2016) never raises above lowest neighbour: " + fx.name,
              n_over_lowest == 0,
              "cells above lowest neighbour = " + std::to_string(n_over_lowest) +
              ", worst excess = " + std::to_string(worst_over) + " m");
    }

    {
        // Least-cost breaching has no RichDEM counterpart, so
        // property tests are the only correctness check it will ever get.
        // fill_remainder=true means the result should also be depression-free.
        const lsdtt::lsdtt_raster_cpu breached = lsdtt::breach_depressions_least_cost(
            fx.dem, /*max_dist=*/100, /*max_cost=*/std::numeric_limits<float>::max(),
            /*min_dist=*/true, /*flat_increment=*/0.0f, /*fill_remainder=*/true);

        check("P3 breach (least_cost, fill_remainder) is depression-free: " + fx.name,
              count_interior_minima(breached) == 0,
              "interior minima = " + std::to_string(count_interior_minima(breached)));
    }
}


// =============================================================================
// P4 / P5 -- D8 accumulation properties
// =============================================================================
void test_d8_accumulation_properties(const fixture& fx)
{
    const lsdtt::lsdtt_raster_cpu filled = lsdtt::fill_pits(fx.dem, 0.0001f);

    lsdtt::lsdtt_flowinfo_options opts;
    opts.boundary_conditions = "oooo";

    lsdtt::lsdtt_flowinfo fi(filled.get_data_array(),
                             filled.get_mask_array(),
                             filled.meta(),
                             opts);

    const std::vector<double> accum = fi.accumulate();
    const lsdtt::lsdtt_raster_meta&  meta  = filled.meta();
    const double cell_area = meta.cell_size * meta.cell_size;

    // -- P4: mass conservation ---------------------------------------------
    // Every valid cell contributes one cell area and drains to exactly one
    // outlet, so the outlet total must equal the whole valid area.
    double outlet_total = 0.0;
    std::int64_t n_valid = 0;
    for (std::int64_t flat = 0; flat < meta.n_pixels(); ++flat)
    {
        if (!fi.is_valid_cell(flat)) continue;
        ++n_valid;
        if (fi.is_outlet(flat)) outlet_total += accum[flat];
    }

    const double expected = static_cast<double>(n_valid) * cell_area;
    const double rel_err  = (expected > 0.0)
                          ? std::abs(outlet_total - expected) / expected
                          : 0.0;

    check("P4 D8 mass conservation: " + fx.name,
          rel_err <= 1e-9,
          "outlet total = " + std::to_string(outlet_total) +
          ", expected = " + std::to_string(expected) +
          ", rel err = " + std::to_string(rel_err));

    // -- P5: accumulation increases downstream ------------------------------
    std::int64_t n_violations = 0;
    double worst = 0.0;
    for (std::int64_t flat = 0; flat < meta.n_pixels(); ++flat)
    {
        if (!fi.is_valid_cell(flat) || fi.is_outlet(flat)) continue;

        const std::int64_t rec = fi.receiver_of_unchecked(flat);
        if (rec < 0) continue;

        // Strictly, accum[rec] >= accum[flat] + cell_area, since the receiver
        // also contributes its own cell.  Test the weaker >= to stay robust to
        // float accumulation order, and record the worst shortfall.
        const double deficit = accum[flat] - accum[rec];
        if (deficit > 1e-6 * cell_area)
        {
            ++n_violations;
            if (deficit > worst) worst = deficit;
        }
    }

    check("P5 D8 accumulation is monotonic downstream: " + fx.name,
          n_violations == 0,
          "violations = " + std::to_string(n_violations) +
          ", worst deficit = " + std::to_string(worst));
}


// =============================================================================
// P8 -- exact D8 accumulation on a planar ramp (analytic oracle)
//
// A surface tilted purely north-south has an exactly known answer: every cell
// drains straight down its own column, so row r accumulates r+1 cells.
//
// This is the check that distinguishes true steepest descent (drop / distance)
// from "lowest neighbour by elevation". On a planar ramp the S, SE and SW
// neighbours all sit at the same elevation, but the diagonals are sqrt(2)
// further away, so only the cardinal S is the steepest DESCENT. An
// implementation that compares raw elevation instead ties, breaks toward a
// diagonal, and produces diagonal stripes rather than columns.
//
// RichDEM does exactly that (richdem/flowmet/d8_flowdirs.hpp:63-71 selects the
// minimum-elevation neighbour with a diagonal-preferring tie-break), which is
// why it is NOT a usable reference for D8 accumulation -- see
// docs/testing_and_benchmarking.md Sec7d. Hence this analytic oracle rather
// than a comparison.
// =============================================================================
void test_d8_accumulation_on_ramp()
{
    const int    n         = 50;
    const double cell_size = 30.0;

    lsdtt::lsdtt_raster_meta meta = make_meta(n, n, cell_size);
    lsdtt::lsdtt_array2d_f   data(n, n, 0.0f);
    lsdtt::lsdtt_array2d_u8  mask(n, n, 1);

    for (int r = 0; r < n; ++r)
        for (int c = 0; c < n; ++c)
            data(r, c) = static_cast<float>(100.0 - r);   // tilts north -> south

    lsdtt::lsdtt_raster_cpu ramp;
    ramp.set_from_vectors(std::move(data.get_vector()),
                          std::move(mask.get_vector()), meta);

    lsdtt::lsdtt_flowinfo_options opts;
    opts.boundary_conditions = "oooo";
    lsdtt::lsdtt_flowinfo fi(ramp.get_data_array(), ramp.get_mask_array(),
                             ramp.meta(), opts);

    const std::vector<double> accum = fi.accumulate();
    const double cell_area = cell_size * cell_size;

    std::int64_t n_wrong = 0;
    double worst = 0.0;
    int    worst_row = -1;

    // Skip the outermost columns: with open boundaries their upslope set is
    // clipped by the edge, so the column-count oracle does not apply there.
    for (int r = 0; r < n; ++r)
    {
        for (int c = 1; c < n - 1; ++c)
        {
            const double expected = static_cast<double>(r + 1) * cell_area;
            const double got      = accum[static_cast<std::size_t>(r) * n + c];
            if (std::fabs(got - expected) > 1e-6 * cell_area)
            {
                ++n_wrong;
                if (std::fabs(got - expected) > worst)
                {
                    worst     = std::fabs(got - expected);
                    worst_row = r;
                }
            }
        }
    }

    check("P8 D8 accumulation exact on a planar ramp",
          n_wrong == 0,
          "cells wrong = " + std::to_string(n_wrong) +
          ", worst error = " + std::to_string(worst / cell_area) +
          " cells at row " + std::to_string(worst_row));

    // State the discriminating property directly: flow must go S, not SE/SW.
    std::int64_t n_diagonal = 0;
    for (int r = 0; r < n - 1; ++r)
    {
        for (int c = 1; c < n - 1; ++c)
        {
            const std::int64_t flat = static_cast<std::int64_t>(r) * n + c;
            if (!fi.is_valid_cell(flat) || fi.is_outlet(flat)) continue;
            // Direction code 4 == due south (see lsdtt_flowinfo_cpu.hpp).
            if (fi.receiver_dir_of_unchecked(flat) != 4) ++n_diagonal;
        }
    }

    check("P8 ramp routes due south, not diagonally",
          n_diagonal == 0,
          "cells not routing south = " + std::to_string(n_diagonal) +
          " (steepest DESCENT must divide the drop by distance)");
}


// =============================================================================
// P7 -- D-inf / MFD weight normalisation
// =============================================================================
void test_flowgraph_weight_properties(const fixture& fx)
{
    const lsdtt::lsdtt_raster_cpu filled = lsdtt::fill_pits(fx.dem, 0.0001f);

    struct method_case
    {
        std::string                  name;
        lsdtt::lsdtt_flowgraph_method method;
    };

    const std::vector<method_case> cases = {
        {"d_infinity", lsdtt::lsdtt_flowgraph_method::d_infinity},
        {"mfd",        lsdtt::lsdtt_flowgraph_method::mfd},
    };

    for (const auto& mc : cases)
    {
        lsdtt::lsdtt_flowgraph_options opts;
        opts.boundary_conditions = "oooo";
        opts.method              = mc.method;

        lsdtt::lsdtt_flowgraph_cpu fg(filled.get_data_array(),
                                      filled.get_mask_array(),
                                      filled.meta(),
                                      opts);

        const auto& offsets = fg.receiver_offsets();
        const auto& weights = fg.receiver_weights();

        std::int64_t n_bad = 0;
        std::int64_t n_with_receivers = 0;
        double worst = 0.0;

        for (std::size_t flat = 0; flat + 1 < offsets.size(); ++flat)
        {
            const std::int64_t begin = offsets[flat];
            const std::int64_t end   = offsets[flat + 1];
            if (end == begin) continue;   // sink: no receivers, nothing to sum

            ++n_with_receivers;

            double sum = 0.0;
            for (std::int64_t k = begin; k < end; ++k) sum += weights[k];

            const double err = std::abs(sum - 1.0);
            if (err > 1e-6)
            {
                ++n_bad;
                if (err > worst) worst = err;
            }
        }

        check("P7 " + mc.name + " weights sum to 1: " + fx.name,
              n_bad == 0,
              "cells with receivers = " + std::to_string(n_with_receivers) +
              ", bad = " + std::to_string(n_bad) +
              ", worst |sum-1| = " + std::to_string(worst));

        // D-inf routes to at most two receivers (Tarboton 1997); a cell with
        // three would mean the facet logic has broken down.
        if (mc.method == lsdtt::lsdtt_flowgraph_method::d_infinity)
        {
            std::int64_t n_over = 0;
            for (std::size_t flat = 0; flat + 1 < offsets.size(); ++flat)
                if (offsets[flat + 1] - offsets[flat] > 2) ++n_over;

            check("P7 d_infinity has <= 2 receivers per cell: " + fx.name,
                  n_over == 0,
                  "cells with >2 receivers = " + std::to_string(n_over));
        }
    }
}


// =============================================================================
// Section 8 -- fill algorithm characterisation (reported, not asserted)
//
// How the alternative depression-removal routine compares to Barnes 2014 on
// speed and on how far the answer diverges.  The other fill exposed by the CLI
// is Fill-Spill-Merge (Barnes et al. 2021), which is a hydrologically different
// operation -- it distributes a finite rainfall volume, so depressions that
// receive less water than their capacity stay unfilled BY DESIGN.  It is
// therefore characterised separately in Phase 5, not held to P1.
// =============================================================================
void characterise_fill_algorithms(const std::vector<fixture>& fixtures,
                                  const std::string& output_dir)
{
    std::cout << "\n--- Fill algorithm characterisation (reported, not asserted) ---\n";

    const std::string csv_path = output_dir + "/fill_algorithm_characterisation.csv";
    std::ofstream csv(csv_path);
    if (!csv)
    {
        std::cout << "  WARNING: could not open " << csv_path << " for writing\n";
        return;
    }

    csv << "fixture,n_cells,algorithm,seconds,max_abs_diff_vs_barnes2014,"
           "mean_abs_diff_vs_barnes2014,frac_cells_differing,interior_minima_remaining\n";

    for (const auto& fx : fixtures)
    {
        const std::int64_t n_cells = fx.dem.meta().n_pixels();

        // Reference: Barnes 2014, no epsilon.
        const auto t0 = std::chrono::steady_clock::now();
        const lsdtt::lsdtt_raster_cpu barnes = lsdtt::fill_pits(fx.dem, 0.0f);
        const double barnes_s = seconds_since(t0);

        csv << fx.name << ',' << n_cells << ",barnes2014," << barnes_s
            << ",0,0,0," << count_interior_minima(barnes) << '\n';

        // Alternative currently reachable from the CLI: least-cost breaching
        // with fill remainder.  Fill-Spill-Merge needs a rainfall rate and a
        // depression hierarchy, so it is characterised separately in Phase 5
        // rather than bolted on here.
        const auto t1 = std::chrono::steady_clock::now();
        const lsdtt::lsdtt_raster_cpu lc = lsdtt::breach_depressions_least_cost(
            fx.dem, /*max_dist=*/100, /*max_cost=*/std::numeric_limits<float>::max(),
            /*min_dist=*/true, /*flat_increment=*/0.0f, /*fill_remainder=*/true);
        const double lc_s = seconds_since(t1);

        const diff_stats d = compare(lc, barnes, 0.01);
        const double frac = (d.n_compared > 0)
                          ? static_cast<double>(d.n_differing) / static_cast<double>(d.n_compared)
                          : 0.0;

        csv << fx.name << ',' << n_cells << ",breach_least_cost," << lc_s << ','
            << d.max_abs << ',' << d.mean_abs << ',' << frac << ','
            << count_interior_minima(lc) << '\n';

        std::cout << "  " << std::left << std::setw(28) << fx.name
                  << " barnes2014 " << std::fixed << std::setprecision(4) << barnes_s << " s"
                  << " | least_cost " << lc_s << " s"
                  << " | max diff " << std::setprecision(3) << d.max_abs << " m"
                  << " | cells differing " << std::setprecision(4) << (100.0 * frac) << " %\n";
    }

    csv.close();
    std::cout << "  wrote " << csv_path << "\n";
}

} // namespace


int main(int argc, char** argv)
{
    const std::string output_dir = (argc > 1) ? argv[1] : ".";

    std::cout << "lsdtt3 flow-routing property tests "
                 " (docs/testing_and_benchmarking.md)\n\n";

    const std::vector<fixture> fixtures = build_fixtures();

    for (const auto& fx : fixtures)
    {
        std::cout << "--- fixture: " << fx.name << " ("
                  << fx.dem.meta().n_rows << "x" << fx.dem.meta().n_cols
                  << ", cell " << fx.dem.meta().cell_size << " m) ---\n";

        test_fill_properties(fx);
        test_breach_properties(fx);
        test_d8_accumulation_properties(fx);
        test_flowgraph_weight_properties(fx);
        std::cout << "\n";
    }

    // Analytic oracle: independent of any fixture, so run it once.
    std::cout << "--- analytic: planar ramp ---\n";
    test_d8_accumulation_on_ramp();
    std::cout << "\n";

    characterise_fill_algorithms(fixtures, output_dir);

    std::cout << "\n" << g_n_checks - g_n_failed << " / " << g_n_checks
              << " checks passed\n";
    std::cout << "RESULT: " << (g_n_failed == 0 ? "PASS" : "FAIL") << "\n";

    return (g_n_failed == 0) ? 0 : 1;
}
