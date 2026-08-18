// test_internal_drainage.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
//
// -----------------------------------------------------------------------------
// Unit tests for phases 0, 1 and 3 of the internal-drainage work
// (docs/plan_docs/internal_drainage_plan.md, section 6.1).
//
// Detection and the vector round-trip remain in the separate phase 2b test;
// the CLI belongs to phase 4.
//
//   Phase 0
//     generate_crater_field_surface()  the analytic oracle: a tilted plane with
//                                      closed depressions cut into it, whose
//                                      pit, spill, area and volume are known in
//                                      closed form
//     fill_pits_with_sinks()           priority-flood filling in which the
//                                      designated sinks are drainage exits and
//                                      are never raised
//
//   Phase 1
//     lsdtt_flowinfo_options::forced_outlet_flats
//                                      caller-nominated flow terminations,
//                                      honoured by build_receiver_graph_d8()
//
//   Phase 3
//     punch_sink_holes()               immutable nodata-hole preparation
//     tag_upslope_of_sinks()           forced-outlet and nodata-hole basin tags,
//                                      cutoffs, conservation and validation
//
// Everything is generated in process, so there is no fixture on disk and the
// driver is fast enough for the inner loop.
//
// Usage: lsdtt-test-internal-drainage
// Prints one PASS/FAIL line per check and a final "RESULT: PASS|FAIL".
// -----------------------------------------------------------------------------

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include "lsdtt_raster_cpu.hpp"
#include "lsdtt_raster_flow_preprocessing_cpu.hpp"
#include "lsdtt_raster_synthetic_surfaces_cpu.hpp"
#include "lsdtt_flowinfo_cpu.hpp"
#include "lsdtt_internal_drainage_cpu.hpp"

namespace
{

// =============================================================================
// Tiny check harness (same shape as tests/cpp/test_flow_properties.cpp)
// =============================================================================
int g_n_checks = 0;
int g_n_failed = 0;

void check(const std::string& name, bool ok, const std::string& detail = "")
{
    ++g_n_checks;
    if (!ok) ++g_n_failed;
    std::cout << (ok ? "PASS" : "FAIL") << "  " << name;
    if (!detail.empty()) std::cout << "   [" << detail << "]";
    std::cout << std::endl;
}

// Runs f and reports whether it threw std::invalid_argument whose message
// contains `needle`. The message is checked because these validators exist to
// tell the caller which entry was wrong; a bare throw would not do that.
template <typename F>
void check_throws(const std::string& name, const F& f, const std::string& needle)
{
    bool        threw = false;
    bool        matched = false;
    std::string what;

    try
    {
        f();
    }
    catch (const std::invalid_argument& e)
    {
        threw   = true;
        what    = e.what();
        matched = what.find(needle) != std::string::npos;
    }
    catch (const std::exception& e)
    {
        what = std::string("wrong exception type: ") + e.what();
    }

    check(name, threw && matched,
          threw ? ("message: " + what) : ("did not throw; " + what));
}


// =============================================================================
// Geometry helpers
// =============================================================================

double pixel_centre_x(const lsdtt::lsdtt_raster_meta& meta, int col)
{
    return meta.x_min + (static_cast<double>(col) + 0.5) * meta.cell_size;
}

double pixel_centre_y(const lsdtt::lsdtt_raster_meta& meta, int row)
{
    return meta.y_max() - (static_cast<double>(row) + 0.5) * meta.cell_size;
}

std::int64_t flat_of(const lsdtt::lsdtt_raster_meta& meta, int row, int col)
{
    return static_cast<std::int64_t>(row) *
           static_cast<std::int64_t>(meta.n_cols) +
           static_cast<std::int64_t>(col);
}


// D8 neighbour offsets in the lsdtt direction-code order:
//   7 0 1
//   6 x 2
//   5 4 3
const int k_d_row[8] = { -1, -1,  0,  1,  1,  1,  0, -1 };
const int k_d_col[8] = {  0,  1,  1,  1,  0, -1, -1, -1 };


// -----------------------------------------------------------------------------
// steepest_descent_terminus
//
// Walks the steepest descent from `start` over the supplied surface and returns
// the flat index of the cell where the walk STOPS -- the first cell with no
// strictly lower neighbour, or a cell on the raster edge from which the walk
// leaves the grid (returned as -1).
//
// This is deliberately independent of lsdtt_flowinfo: the point of the
// fill_pits_with_sinks tests is the SURFACE, and routing the surface with the
// object whose behaviour is being characterised elsewhere would make the check
// circular.
// -----------------------------------------------------------------------------
std::int64_t steepest_descent_terminus(const lsdtt::lsdtt_raster_cpu& dem,
                                       std::int64_t start)
{
    const lsdtt::lsdtt_raster_meta& meta = dem.get_meta();
    const auto& z    = dem.get_data_array();
    const auto& mask = dem.get_mask_array();

    const std::int64_t n_cells =
        static_cast<std::int64_t>(meta.n_rows) *
        static_cast<std::int64_t>(meta.n_cols);

    std::int64_t here  = start;
    std::int64_t steps = 0;

    while (steps <= n_cells)
    {
        const int row = static_cast<int>(here / meta.n_cols);
        const int col = static_cast<int>(here % meta.n_cols);

        if (row == 0 || col == 0 ||
            row == meta.n_rows - 1 || col == meta.n_cols - 1)
        {
            return -1;   // left the grid
        }

        const double z_here    = static_cast<double>(z(row, col));
        double       best_drop = 0.0;
        std::int64_t best      = -1;

        for (int dir = 0; dir < 8; ++dir)
        {
            const int nr = row + k_d_row[dir];
            const int nc = col + k_d_col[dir];
            if (mask(nr, nc) == 0) continue;

            const double dz = z_here - static_cast<double>(z(nr, nc));
            if (dz <= 0.0) continue;

            const double dist = (dir % 2 == 0) ? meta.cell_size
                                               : meta.cell_size * std::sqrt(2.0);
            const double slope = dz / dist;
            if (slope > best_drop)
            {
                best_drop = slope;
                best      = flat_of(meta, nr, nc);
            }
        }

        if (best < 0) return here;   // local minimum: the walk terminates here
        here = best;
        ++steps;
    }

    return -2;   // cycled: cannot happen on a filled surface
}


// -----------------------------------------------------------------------------
// steepest_descent_passes_through
//
// True if the steepest-descent walk from `start` ever steps onto `target`.
// This is REACHABILITY, and the plan (section 6.1) records that it does not
// discriminate: after a plain fill a handful of cells still step down onto the
// old pit incidentally. It is computed here only so the driver can report the
// contrast alongside the termination counts that are actually asserted.
// -----------------------------------------------------------------------------
bool steepest_descent_passes_through(const lsdtt::lsdtt_raster_cpu& dem,
                                     std::int64_t start,
                                     std::int64_t target)
{
    const lsdtt::lsdtt_raster_meta& meta = dem.get_meta();
    const auto& z    = dem.get_data_array();
    const auto& mask = dem.get_mask_array();

    const std::int64_t n_cells =
        static_cast<std::int64_t>(meta.n_rows) *
        static_cast<std::int64_t>(meta.n_cols);

    std::int64_t here = start;

    for (std::int64_t steps = 0; steps <= n_cells; ++steps)
    {
        if (here == target) return true;

        const int row = static_cast<int>(here / meta.n_cols);
        const int col = static_cast<int>(here % meta.n_cols);

        if (row == 0 || col == 0 ||
            row == meta.n_rows - 1 || col == meta.n_cols - 1) return false;

        const double z_here    = static_cast<double>(z(row, col));
        double       best_drop = 0.0;
        std::int64_t best      = -1;

        for (int dir = 0; dir < 8; ++dir)
        {
            const int nr = row + k_d_row[dir];
            const int nc = col + k_d_col[dir];
            if (mask(nr, nc) == 0) continue;

            const double dz = z_here - static_cast<double>(z(nr, nc));
            if (dz <= 0.0) continue;

            const double dist = (dir % 2 == 0) ? meta.cell_size
                                               : meta.cell_size * std::sqrt(2.0);
            const double slope = dz / dist;
            if (slope > best_drop)
            {
                best_drop = slope;
                best      = flat_of(meta, nr, nc);
            }
        }

        if (best < 0) return false;
        here = best;
    }

    return false;
}


// Follows the D8 receiver chain to its terminus and returns that flat index.
std::int64_t receiver_terminus(const lsdtt::lsdtt_flowinfo& flow,
                               std::int64_t start,
                               std::int64_t max_steps)
{
    std::int64_t here  = start;
    std::int64_t steps = 0;

    while (steps <= max_steps)
    {
        const std::int64_t rec = flow.receiver_of_unchecked(here);
        if (rec == lsdtt::lsdtt_flowinfo::outlet_receiver) return here;
        if (rec == lsdtt::lsdtt_flowinfo::nodata_receiver) return -1;
        if (rec == here) return here;
        here = rec;
        ++steps;
    }

    return -2;
}


// =============================================================================
// Fixtures
// =============================================================================

// A single cone crater on a tilted plane. Centre placed exactly on a pixel
// centre, which the generator's header requires for the analytic pit to be the
// raster minimum (see the phase-0 note in section 6.1 of the plan).
struct single_crater_fixture
{
    lsdtt::lsdtt_crater_field_surface gt;
    int    pit_row          = 100;
    int    pit_col          = 100;
    double radius_m         = 300.0;
    double depth_m          = 100.0;
    double background_slope = 0.02;
};

single_crater_fixture make_single_crater()
{
    single_crater_fixture fx;

    const lsdtt::lsdtt_raster_meta meta =
        [] {
            lsdtt::lsdtt_raster_meta m;
            m.n_rows    = 300;
            m.n_cols    = 300;
            m.cell_size = 10.0;
            m.x_min     = 0.0;
            m.y_min     = 0.0;
            m.set_default_synthetic_georeference();
            return m;
        }();

    lsdtt::lsdtt_crater_spec cs;
    cs.centre_x       = pixel_centre_x(meta, fx.pit_col);
    cs.centre_y       = pixel_centre_y(meta, fx.pit_row);
    cs.radius_m       = fx.radius_m;
    cs.depth_m        = fx.depth_m;
    cs.floor_fraction = 0.0;

    fx.gt = lsdtt::generate_crater_field_surface(
        meta, { cs }, fx.background_slope, 1000.0);

    return fx;
}


// Three non-overlapping cone craters on the same tilted plane. This is the
// fixture the routing tests use: several independent internal basins plus a
// large externally drained remainder.
lsdtt::lsdtt_crater_field_surface make_crater_field()
{
    lsdtt::lsdtt_raster_meta meta;
    meta.n_rows    = 300;
    meta.n_cols    = 300;
    meta.cell_size = 10.0;
    meta.x_min     = 0.0;
    meta.y_min     = 0.0;
    meta.set_default_synthetic_georeference();

    const int rows[3] = { 100, 200,  60 };
    const int cols[3] = { 100, 210, 230 };
    const double radii[3]  = { 300.0, 200.0, 150.0 };
    const double depths[3] = { 100.0,  80.0,  60.0 };

    std::vector<lsdtt::lsdtt_crater_spec> specs;
    for (int i = 0; i < 3; ++i)
    {
        lsdtt::lsdtt_crater_spec cs;
        cs.centre_x       = pixel_centre_x(meta, cols[i]);
        cs.centre_y       = pixel_centre_y(meta, rows[i]);
        cs.radius_m       = radii[i];
        cs.depth_m        = depths[i];
        cs.floor_fraction = 0.0;
        specs.push_back(cs);
    }

    return lsdtt::generate_crater_field_surface(meta, specs, 0.02, 1000.0);
}


// A plain south-draining tilted ramp with no depressions at all. Used by the
// short-circuit-ordering test, where the point is that the forced outlet has a
// perfectly good downhill neighbour.
lsdtt::lsdtt_raster_cpu make_tilted_ramp(int n_rows, int n_cols, double cell_size)
{
    lsdtt::lsdtt_raster_meta meta;
    meta.n_rows    = n_rows;
    meta.n_cols    = n_cols;
    meta.cell_size = cell_size;
    meta.x_min     = 0.0;
    meta.y_min     = 0.0;
    meta.set_default_synthetic_georeference();

    // A crater field with no craters is exactly the background plane, so this
    // reuses the same generator rather than introducing a second surface
    // definition that could drift away from it.
    const lsdtt::lsdtt_crater_field_surface s =
        lsdtt::generate_crater_field_surface(
            meta, std::vector<lsdtt::lsdtt_crater_spec>{}, 0.02, 1000.0);

    return s.dem;
}


// =============================================================================
// Section 1: the crater generator as an oracle
// =============================================================================

void test_single_crater_pit_location()
{
    const single_crater_fixture fx = make_single_crater();
    const lsdtt::lsdtt_raster_cpu& dem = fx.gt.dem;
    const lsdtt::lsdtt_raster_meta& meta = dem.get_meta();
    const auto& z = dem.get_data_array();

    // Scan for the raster minimum, and count how many cells attain it, so a
    // near-tie is reported rather than being hidden by whichever cell the scan
    // happened to reach first.
    float        z_min = z.data()[0];
    std::int64_t argmin = 0;
    std::int64_t n_at_min = 1;

    const std::int64_t n_cells =
        static_cast<std::int64_t>(meta.n_rows) *
        static_cast<std::int64_t>(meta.n_cols);

    for (std::int64_t i = 1; i < n_cells; ++i)
    {
        const float v = z.data()[i];
        if (v < z_min) { z_min = v; argmin = i; n_at_min = 1; }
        else if (v == z_min) { ++n_at_min; }
    }

    check("single_crater_pit_location: raster minimum is the analytic pit",
          argmin == fx.gt.pit_flats[0],
          "raster argmin " + std::to_string(argmin) +
          ", analytic pit " + std::to_string(fx.gt.pit_flats[0]));

    check("single_crater_pit_location: cone pit is a unique minimum",
          n_at_min == 1,
          std::to_string(n_at_min) + " cells at the minimum");

    check("single_crater_pit_location: sink_mask marks exactly the pit",
          fx.gt.sink_mask.get_data_array().data()[fx.gt.pit_flats[0]] == 1.0f);
}


void test_single_crater_geometry()
{
    const single_crater_fixture fx = make_single_crater();
    const lsdtt::lsdtt_raster_cpu& dem  = fx.gt.dem;
    const lsdtt::lsdtt_raster_meta& meta = dem.get_meta();
    const auto& z = dem.get_data_array();

    const double cell_area = meta.cell_size * meta.cell_size;

    // --- pit elevation -------------------------------------------------------
    // The centre sits on a pixel centre, so d == 0 there, w == 0, and the
    // generator's own formula collapses to z_pit exactly. Only float rounding
    // separates the two.
    const double z_pit_raster =
        static_cast<double>(z.data()[fx.gt.pit_flats[0]]);
    check("single_crater_geometry: pit elevation matches the oracle",
          std::fabs(z_pit_raster - fx.gt.pit_elevation[0]) < 1.0e-3,
          "raster " + std::to_string(z_pit_raster) +
          " vs analytic " + std::to_string(fx.gt.pit_elevation[0]));

    // --- spill elevation, area, volume from the fill difference --------------
    // Filling with no imposed gradient raises the depression to exactly the
    // level of its lowest rim cell, so the filled value at the pit IS the
    // measured spill elevation and {filled > dem} IS the depression footprint.
    const lsdtt::lsdtt_raster_cpu filled = lsdtt::fill_pits(dem, 0.0f);
    const auto& zf = filled.get_data_array();

    const double z_spill_raster =
        static_cast<double>(zf.data()[fx.gt.pit_flats[0]]);

    // The grid resolves the rim only to the nearest cell, so the measured
    // spill sits above the analytic minimum of the continuous rim by at most
    // one cell diagonal of rim relief.
    const double rim_gradient = fx.depth_m / fx.radius_m + fx.background_slope;
    const double spill_tol    = rim_gradient * meta.cell_size * std::sqrt(2.0);

    check("single_crater_geometry: spill elevation matches the oracle",
          z_spill_raster >= fx.gt.spill_elevation[0] - 1.0e-3 &&
          z_spill_raster <= fx.gt.spill_elevation[0] + spill_tol,
          "raster " + std::to_string(z_spill_raster) +
          " vs analytic " + std::to_string(fx.gt.spill_elevation[0]) +
          " (tol +" + std::to_string(spill_tol) + ")");

    std::int64_t n_fill  = 0;
    double       vol_sum = 0.0;

    const std::int64_t n_cells =
        static_cast<std::int64_t>(meta.n_rows) *
        static_cast<std::int64_t>(meta.n_cols);

    for (std::int64_t i = 0; i < n_cells; ++i)
    {
        const double d = static_cast<double>(zf.data()[i]) -
                         static_cast<double>(z.data()[i]);
        if (d > 0.0)
        {
            ++n_fill;
            vol_sum += d * cell_area;
        }
    }

    // The pit cell itself is filled by zero (it IS the spill level only after
    // filling), so add it back: the depression footprint includes its pit.
    const double area_raster = static_cast<double>(n_fill + 1) * cell_area;

    // Discretisation error is one cell along the whole depression boundary,
    // i.e. of order 2 * pi * R * cell_size in area.
    const double area_tol = 2.0 * M_PI * fx.radius_m * meta.cell_size;

    check("single_crater_geometry: depression area matches the oracle",
          std::fabs(area_raster - fx.gt.area_m2[0]) < area_tol,
          "raster " + std::to_string(area_raster) +
          " vs analytic " + std::to_string(fx.gt.area_m2[0]) +
          " (tol " + std::to_string(area_tol) + ")");

    // ... and the same boundary strip in volume, weighted by the mean depth
    // there, which for a cone is a small fraction of the full depth. One cell
    // of slop on every linear dimension is the header's stated expectation, so
    // compare relatively at the cell_size / radius scale.
    const double vol_tol = 3.0 * (meta.cell_size / fx.radius_m) *
                           fx.gt.volume_m3[0];

    check("single_crater_geometry: fill volume matches the oracle",
          std::fabs(vol_sum - fx.gt.volume_m3[0]) < vol_tol,
          "raster " + std::to_string(vol_sum) +
          " vs analytic " + std::to_string(fx.gt.volume_m3[0]) +
          " (tol " + std::to_string(vol_tol) + ")");
}


void test_flat_floor_tie_break()
{
    lsdtt::lsdtt_raster_meta meta;
    meta.n_rows    = 200;
    meta.n_cols    = 200;
    meta.cell_size = 10.0;
    meta.x_min     = 0.0;
    meta.y_min     = 0.0;
    meta.set_default_synthetic_georeference();

    lsdtt::lsdtt_crater_spec cs;
    cs.centre_x       = pixel_centre_x(meta, 100);
    cs.centre_y       = pixel_centre_y(meta, 100);
    cs.radius_m       = 400.0;
    cs.depth_m        = 100.0;
    cs.floor_fraction = 0.4;      // a 160 m radius floor: 16 cells across

    const lsdtt::lsdtt_crater_field_surface a =
        lsdtt::generate_crater_field_surface(meta, { cs }, 0.02, 1000.0);
    const lsdtt::lsdtt_crater_field_surface b =
        lsdtt::generate_crater_field_surface(meta, { cs }, 0.02, 1000.0);

    const auto& za = a.dem.get_data_array();
    const auto& zb = b.dem.get_data_array();

    const std::int64_t n_cells =
        static_cast<std::int64_t>(meta.n_rows) *
        static_cast<std::int64_t>(meta.n_cols);

    // The generator claims to be deterministic with no RNG anywhere.
    check("flat_floor_tie_break: generator is bit-identical across runs",
          std::memcmp(za.data(), zb.data(),
                      static_cast<std::size_t>(n_cells) * sizeof(float)) == 0);

    // Collect the floor: the set of cells at the exact minimum elevation. The
    // generator fades the regional tilt out towards the floor precisely so that
    // this is an exact tie set, not merely a nearly-flat one.
    float z_min = za.data()[0];
    for (std::int64_t i = 1; i < n_cells; ++i)
    {
        if (za.data()[i] < z_min) z_min = za.data()[i];
    }

    std::vector<std::int64_t> floor_set;
    for (std::int64_t i = 0; i < n_cells; ++i)
    {
        if (za.data()[i] == z_min) floor_set.push_back(i);
    }

    check("flat_floor_tie_break: the floor is an exact tie set of many cells",
          floor_set.size() > 50,
          std::to_string(floor_set.size()) + " cells at the exact minimum");

    check("flat_floor_tie_break: floor elevation is the analytic pit elevation",
          std::fabs(static_cast<double>(z_min) - a.pit_elevation[0]) < 1.0e-3,
          "z_min " + std::to_string(z_min) +
          " vs analytic " + std::to_string(a.pit_elevation[0]));

    // A lowest-flat-index tie-break is the rule the plan fixes (section 5.5).
    // It must be reproducible, and it must land ON the floor -- but it will NOT
    // agree with the generator's reported pit_rc, because the generator reports
    // the cell containing the centre while this rule picks the north-westmost
    // floor cell. Both are correct answers to "which cell is the pit" for a
    // floor that is an exact tie; asserting equality here would be asserting a
    // coincidence.
    const std::int64_t tie_break_pit = floor_set.front();

    check("flat_floor_tie_break: lowest-flat-index pit is a floor member",
          std::find(floor_set.begin(), floor_set.end(), tie_break_pit) !=
              floor_set.end());

    check("flat_floor_tie_break: the generator's pit_rc is a floor member",
          std::find(floor_set.begin(), floor_set.end(), a.pit_flats[0]) !=
              floor_set.end(),
          "pit_flat " + std::to_string(a.pit_flats[0]));

    check("flat_floor_tie_break: tie-break is stable across a repeat run",
          [&] {
              float m = zb.data()[0];
              for (std::int64_t i = 1; i < n_cells; ++i)
                  if (zb.data()[i] < m) m = zb.data()[i];
              for (std::int64_t i = 0; i < n_cells; ++i)
                  if (zb.data()[i] == m) return i == tie_break_pit;
              return false;
          }());
}


void test_generator_validation()
{
    lsdtt::lsdtt_raster_meta meta;
    meta.n_rows    = 200;
    meta.n_cols    = 200;
    meta.cell_size = 10.0;
    meta.x_min     = 0.0;
    meta.y_min     = 0.0;
    meta.set_default_synthetic_georeference();

    // Overlapping rims: two depressions that merge into one whose pit and spill
    // neither spec describes.
    check_throws(
        "generator_validation: overlapping rims are refused",
        [&] {
            lsdtt::lsdtt_crater_spec a;
            a.centre_x = pixel_centre_x(meta, 90);
            a.centre_y = pixel_centre_y(meta, 100);
            a.radius_m = 300.0;
            a.depth_m  = 100.0;

            lsdtt::lsdtt_crater_spec b = a;
            b.centre_x = pixel_centre_x(meta, 110);

            (void)lsdtt::generate_crater_field_surface(meta, { a, b }, 0.0, 1000.0);
        },
        "overlap");

    // Too shallow to close against the regional tilt: the wall does not outrun
    // the tilt on the down-tilt side, so this is a notch in a hillslope, not a
    // depression.
    check_throws(
        "generator_validation: a crater too shallow for the tilt is refused",
        [&] {
            lsdtt::lsdtt_crater_spec a;
            a.centre_x = pixel_centre_x(meta, 100);
            a.centre_y = pixel_centre_y(meta, 100);
            a.radius_m = 300.0;
            a.depth_m  = 5.0;      // budget is 0.02 * 300 * 2 == 12 m

            (void)lsdtt::generate_crater_field_surface(meta, { a }, 0.02, 1000.0);
        },
        "must exceed background_slope");
}


// =============================================================================
// Section 2: fill_pits_with_sinks
// =============================================================================

void test_fill_with_sinks_retains_descent()
{
    const lsdtt::lsdtt_crater_field_surface gt = make_crater_field();
    const lsdtt::lsdtt_raster_cpu& dem  = gt.dem;
    const lsdtt::lsdtt_raster_meta& meta = dem.get_meta();

    const float min_slope = 0.0001f;

    const lsdtt::lsdtt_raster_cpu filled_sinks =
        lsdtt::fill_pits_with_sinks(dem, gt.pit_flats, min_slope);
    const lsdtt::lsdtt_raster_cpu filled_plain =
        lsdtt::fill_pits_priority_flood(dem, min_slope);

    const auto& id_data = gt.crater_id.get_data_array();
    const auto& id_mask = gt.crater_id.get_mask_array();

    const std::int64_t n_cells =
        static_cast<std::int64_t>(meta.n_rows) *
        static_cast<std::int64_t>(meta.n_cols);

    std::int64_t n_crater      = 0;
    std::int64_t n_term_sinks  = 0;   // terminates at its own analytic pit
    std::int64_t n_term_plain  = 0;
    std::int64_t n_reach_plain = 0;   // merely STEPS ON the pit at some point

    for (std::int64_t i = 0; i < n_cells; ++i)
    {
        if (id_mask.data()[i] == 0) continue;
        ++n_crater;

        const std::size_t which =
            static_cast<std::size_t>(std::lround(id_data.data()[i]));
        const std::int64_t own_pit = gt.pit_flats[which];

        if (steepest_descent_terminus(filled_sinks, i) == own_pit) ++n_term_sinks;
        if (steepest_descent_terminus(filled_plain, i) == own_pit) ++n_term_plain;

        // Reachability, for contrast: this is the assertion the plan warns off,
        // because plain filling leaves a handful of cells whose imposed
        // min_slope gradient happens to step down onto the old pit in passing.
        if (steepest_descent_passes_through(filled_plain, i, own_pit))
            ++n_reach_plain;
    }

    std::cout << "       with sinks: " << n_term_sinks << " / " << n_crater
              << " terminate at own pit; plain fill: " << n_term_plain
              << " terminate, " << n_reach_plain << " merely reach\n";

    check("fill_with_sinks_retains_descent: every crater cell terminates at its own pit",
          n_crater > 0 && n_term_sinks == n_crater,
          std::to_string(n_term_sinks) + " / " + std::to_string(n_crater) +
          " crater cells");

    check("fill_with_sinks_retains_descent: plain fill terminates NO crater cell at a pit",
          n_term_plain == 0,
          std::to_string(n_term_plain) + " / " + std::to_string(n_crater) +
          " crater cells");
}


void test_sinks_never_raised()
{
    const lsdtt::lsdtt_crater_field_surface gt = make_crater_field();
    const lsdtt::lsdtt_raster_cpu& dem = gt.dem;

    const lsdtt::lsdtt_raster_cpu filled_sinks =
        lsdtt::fill_pits_with_sinks(dem, gt.pit_flats, 0.0001f);
    const lsdtt::lsdtt_raster_cpu filled_plain =
        lsdtt::fill_pits_priority_flood(dem, 0.0001f);

    const auto& z0 = dem.get_data_array();
    const auto& zs = filled_sinks.get_data_array();
    const auto& zp = filled_plain.get_data_array();

    bool   all_identical = true;
    double max_plain_rise = 0.0;

    for (const std::int64_t f : gt.pit_flats)
    {
        if (zs.data()[f] != z0.data()[f]) all_identical = false;
        max_plain_rise = std::max(
            max_plain_rise,
            static_cast<double>(zp.data()[f]) - static_cast<double>(z0.data()[f]));
    }

    check("sinks_never_raised: pit elevations are bit-identical to the raw DEM",
          all_identical);

    // The contrast that makes the check above meaningful: plain filling raises
    // those same cells to the spill level.
    check("sinks_never_raised: plain fill raises the same cells substantially",
          max_plain_rise > 10.0,
          "max rise " + std::to_string(max_plain_rise) + " m");
}


void test_empty_sink_set_equivalence()
{
    const lsdtt::lsdtt_crater_field_surface gt = make_crater_field();
    const lsdtt::lsdtt_raster_cpu& dem = gt.dem;

    const std::vector<std::int64_t> no_sinks;

    for (const float min_slope : { 0.0f, 0.0001f })
    {
        const lsdtt::lsdtt_raster_cpu a =
            lsdtt::fill_pits_with_sinks(dem, no_sinks, min_slope);
        const lsdtt::lsdtt_raster_cpu b =
            lsdtt::fill_pits_priority_flood(dem, min_slope);

        const std::size_t n_bytes =
            static_cast<std::size_t>(dem.get_meta().n_rows) *
            static_cast<std::size_t>(dem.get_meta().n_cols) * sizeof(float);

        check("empty_sink_set_equivalence: bit-identical to priority flood"
              " (min_slope " + std::to_string(min_slope) + ")",
              std::memcmp(a.get_data_array().data(),
                          b.get_data_array().data(), n_bytes) == 0);
    }
}


void test_fill_with_sinks_validation()
{
    const lsdtt::lsdtt_crater_field_surface gt = make_crater_field();
    lsdtt::lsdtt_raster_cpu dem = gt.dem;

    const lsdtt::lsdtt_raster_meta& meta = dem.get_meta();
    const std::int64_t n_cells =
        static_cast<std::int64_t>(meta.n_rows) *
        static_cast<std::int64_t>(meta.n_cols);

    check_throws(
        "fill_with_sinks_validation: index past the end is refused",
        [&] { (void)lsdtt::fill_pits_with_sinks(dem, { n_cells }, 0.0001f); },
        "outside the raster");

    check_throws(
        "fill_with_sinks_validation: negative index is refused",
        [&] { (void)lsdtt::fill_pits_with_sinks(dem, { -1 }, 0.0001f); },
        "outside the raster");

    // Punch a nodata hole and aim a sink at it. A sink that is not a real cell
    // cannot act as an exit, and dropping it silently would leave the
    // depression filled -- which is the whole failure this function exists to
    // prevent.
    lsdtt::lsdtt_raster_cpu holed = dem;
    auto& hdata = holed.get_data_array();
    auto& hmask = const_cast<lsdtt::lsdtt_array2d_u8&>(holed.get_mask_array());
    const int hr = 20;
    const int hc = 20;
    hdata(hr, hc) = static_cast<float>(meta.no_data_value);
    hmask(hr, hc) = 0;

    const std::int64_t hole_flat = flat_of(meta, hr, hc);

    check_throws(
        "fill_with_sinks_validation: a sink on a nodata cell is refused",
        [&] { (void)lsdtt::fill_pits_with_sinks(holed, { hole_flat }, 0.0001f); },
        "nodata or masked cell");
}


// =============================================================================
// Section 3: forced_outlet_flats
// =============================================================================

// Everything in this section shares one prepared surface: the crater field
// filled with its own pits as sinks, which is the ONLY correct preparation for
// forced outlets (see the trap pinned by
// test_forced_outlet_without_sink_fill below).
struct routed_field
{
    lsdtt::lsdtt_crater_field_surface gt;
    lsdtt::lsdtt_raster_cpu           filled;
};

routed_field make_routed_field()
{
    routed_field rf;
    rf.gt     = make_crater_field();
    rf.filled = lsdtt::fill_pits_with_sinks(rf.gt.dem, rf.gt.pit_flats, 0.0001f);
    return rf;
}


std::vector<lsdtt::sink_candidate> make_routing_sinks(
    const std::vector<std::int64_t>& pits,
    const lsdtt::lsdtt_raster_meta& meta)
{
    static const int durable_ids[] = { 30, 10, 20, 40, 50, 60, 70, 80 };
    std::vector<lsdtt::sink_candidate> sinks;
    sinks.reserve(pits.size());

    for (std::size_t i = 0; i < pits.size(); ++i)
    {
        lsdtt::sink_candidate sink;
        sink.sink_id = durable_ids[i % 8] + static_cast<int>(i / 8) * 100;
        sink.pit_flat = pits[i];
        sink.pit_row = static_cast<int>(pits[i] / meta.n_cols);
        sink.pit_col = static_cast<int>(pits[i] % meta.n_cols);
        sink.accepted = true;
        sinks.push_back(sink);
    }

    return lsdtt::accepted_sinks(sinks);
}


void test_forced_outlet_terminates()
{
    const routed_field rf = make_routed_field();
    const lsdtt::lsdtt_raster_meta& meta = rf.filled.get_meta();

    lsdtt::lsdtt_flowinfo_options opts;
    opts.forced_outlet_flats = rf.gt.pit_flats;

    const lsdtt::lsdtt_flowinfo flow(rf.filled.get_data_array(),
                                     rf.filled.get_mask_array(),
                                     meta, opts);

    const std::int64_t n_cells =
        static_cast<std::int64_t>(meta.n_rows) *
        static_cast<std::int64_t>(meta.n_cols);

    bool all_outlets = true;
    for (const std::int64_t f : rf.gt.pit_flats)
    {
        if (!flow.is_outlet(f)) all_outlets = false;
        if (flow.receiver_of(f) != lsdtt::lsdtt_flowinfo::outlet_receiver)
            all_outlets = false;
        if (flow.receiver_dir_of(f) != lsdtt::lsdtt_flowinfo::no_direction)
            all_outlets = false;

        const auto& o = flow.outlets();
        if (std::find(o.begin(), o.end(), f) == o.end()) all_outlets = false;
    }

    check("forced_outlet_terminates: every pit is an outlet with no direction",
          all_outlets);

    // Every crater cell's path must end at ITS OWN pit -- not at another
    // crater's pit and not at the raster edge.
    const auto& id_data = rf.gt.crater_id.get_data_array();
    const auto& id_mask = rf.gt.crater_id.get_mask_array();

    std::int64_t n_crater   = 0;
    std::int64_t n_own_pit  = 0;
    std::int64_t n_off_grid = 0;

    for (std::int64_t i = 0; i < n_cells; ++i)
    {
        if (id_mask.data()[i] == 0) continue;
        ++n_crater;

        const std::size_t which =
            static_cast<std::size_t>(std::lround(id_data.data()[i]));
        const std::int64_t term = receiver_terminus(flow, i, n_cells);

        if (term == rf.gt.pit_flats[which]) ++n_own_pit;

        // A path that reaches the raster edge has escaped the depression, which
        // is exactly the failure this whole design prevents.
        if (term >= 0)
        {
            const int tr = static_cast<int>(term / meta.n_cols);
            const int tc = static_cast<int>(term % meta.n_cols);
            if (tr == 0 || tc == 0 ||
                tr == meta.n_rows - 1 || tc == meta.n_cols - 1) ++n_off_grid;
        }
    }

    check("forced_outlet_terminates: every crater cell ends at its own pit",
          n_crater > 0 && n_own_pit == n_crater,
          std::to_string(n_own_pit) + " / " + std::to_string(n_crater));

    check("forced_outlet_terminates: no crater cell reaches the raster edge",
          n_off_grid == 0,
          std::to_string(n_off_grid) + " cells");

    // Plan section 5.3b: the internal basin is STRICTLY LARGER than the crater
    // footprint on a tilted background, because the depression intercepts the
    // up-slope strip draining into it. Subset, never equality.
    const lsdtt::lsdtt_raster_result ncp =
        flow.get_n_contributing_pixels_raster();

    bool subset_ok      = true;
    bool strictly_bigger = true;
    bool ncp_matches    = true;

    std::vector<std::uint8_t> in_upslope(static_cast<std::size_t>(n_cells), 0);

    for (std::size_t k = 0; k < rf.gt.pit_flats.size(); ++k)
    {
        const std::int64_t pit = rf.gt.pit_flats[k];
        const lsdtt::lsdtt_flat_index_slice up = flow.upslope_nodes(pit);

        for (const std::int64_t f : up)
            in_upslope[static_cast<std::size_t>(f)] = 1;

        std::int64_t n_this_crater = 0;
        for (std::int64_t i = 0; i < n_cells; ++i)
        {
            if (id_mask.data()[i] == 0) continue;
            if (static_cast<std::size_t>(std::lround(id_data.data()[i])) != k)
                continue;
            ++n_this_crater;
            if (in_upslope[static_cast<std::size_t>(i)] == 0) subset_ok = false;
        }

        if (static_cast<std::int64_t>(up.size()) <= n_this_crater)
            strictly_bigger = false;

        if (static_cast<double>(ncp.data[static_cast<std::size_t>(pit)]) !=
            static_cast<double>(up.size()))
            ncp_matches = false;

        std::cout << "       crater " << k << ": " << n_this_crater
                  << " crater cells, " << up.size() << " upslope cells\n";

        std::fill(in_upslope.begin(), in_upslope.end(), std::uint8_t{0});
    }

    check("forced_outlet_terminates: crater cells are a SUBSET of the pit's upslope set",
          subset_ok);
    check("forced_outlet_terminates: the internal basin is strictly larger than the crater",
          strictly_bigger);
    check("forced_outlet_terminates: n_contributing_pixels equals the upslope set size",
          ncp_matches);
}


void test_forced_outlet_beats_steepest_descent()
{
    // A plain south-draining ramp: no depressions anywhere, so the nominated
    // cell unambiguously HAS a strictly downhill neighbour. If the short circuit
    // were placed after the steepest-descent search this cell would keep that
    // neighbour as its receiver and the check would fail.
    const lsdtt::lsdtt_raster_cpu ramp = make_tilted_ramp(120, 120, 10.0);
    const lsdtt::lsdtt_raster_meta& meta = ramp.get_meta();

    const int row = 60;
    const int col = 60;
    const std::int64_t f = flat_of(meta, row, col);

    // The surface really does descend away from this cell.
    const auto& z = ramp.get_data_array();
    bool has_downhill = false;
    for (int dir = 0; dir < 8; ++dir)
    {
        if (z(row + k_d_row[dir], col + k_d_col[dir]) < z(row, col))
            has_downhill = true;
    }
    check("forced_outlet_beats_steepest_descent: the test cell has a downhill neighbour",
          has_downhill);

    lsdtt::lsdtt_flowinfo_options plain;
    const lsdtt::lsdtt_flowinfo flow_plain(ramp.get_data_array(),
                                           ramp.get_mask_array(), meta, plain);

    check("forced_outlet_beats_steepest_descent: by default the cell routes normally",
          !flow_plain.is_outlet(f) &&
          flow_plain.receiver_of(f) != lsdtt::lsdtt_flowinfo::outlet_receiver,
          "receiver " + std::to_string(flow_plain.receiver_of(f)));

    lsdtt::lsdtt_flowinfo_options forced;
    forced.forced_outlet_flats = { f };
    const lsdtt::lsdtt_flowinfo flow_forced(ramp.get_data_array(),
                                            ramp.get_mask_array(), meta, forced);

    check("forced_outlet_beats_steepest_descent: the forced outlet wins over the neighbour",
          flow_forced.is_outlet(f) &&
          flow_forced.receiver_of(f) == lsdtt::lsdtt_flowinfo::outlet_receiver &&
          flow_forced.receiver_dir_of(f) == lsdtt::lsdtt_flowinfo::no_direction,
          "receiver " + std::to_string(flow_forced.receiver_of(f)) +
          ", dir " + std::to_string(
              static_cast<int>(flow_forced.receiver_dir_of(f))));

    // ... and it really does terminate flow: the cells that used to route
    // through it now stop there.
    check("forced_outlet_beats_steepest_descent: the forced outlet gains a basin",
          flow_forced.upslope_nodes(f).size() > 1,
          std::to_string(flow_forced.upslope_nodes(f).size()) + " upslope cells");
}


void test_empty_forced_outlet_set_equivalence()
{
    // A PLAINLY filled crater field, deliberately: with no forced outlets there
    // must be no sinks left in the surface, and the sink-seeded fill leaves
    // three of them on purpose. This test is about the option's default being
    // inert, so the surface it runs on has to be one the default router can
    // actually route.
    const lsdtt::lsdtt_crater_field_surface gt = make_crater_field();
    const lsdtt::lsdtt_raster_cpu filled =
        lsdtt::fill_pits_priority_flood(gt.dem, 0.0001f);
    const lsdtt::lsdtt_raster_meta& meta = filled.get_meta();

    lsdtt::lsdtt_flowinfo_options defaulted;

    lsdtt::lsdtt_flowinfo_options explicit_empty;
    explicit_empty.forced_outlet_flats.clear();

    const lsdtt::lsdtt_flowinfo a(filled.get_data_array(),
                                  filled.get_mask_array(), meta, defaulted);
    const lsdtt::lsdtt_flowinfo b(filled.get_data_array(),
                                  filled.get_mask_array(), meta, explicit_empty);

    const std::int64_t n_cells =
        static_cast<std::int64_t>(meta.n_rows) *
        static_cast<std::int64_t>(meta.n_cols);

    bool rec_same = true;
    bool dir_same = true;
    for (std::int64_t i = 0; i < n_cells; ++i)
    {
        if (a.receiver_array().data()[i] != b.receiver_array().data()[i])
            rec_same = false;
        if (a.receiver_dir_array().data()[i] != b.receiver_dir_array().data()[i])
            dir_same = false;
    }

    check("empty_forced_outlet_set_equivalence: receiver_ is unchanged", rec_same);
    check("empty_forced_outlet_set_equivalence: receiver_dir_ is unchanged", dir_same);
}


void test_tag_area_conservation_and_disjointness()
{
    const routed_field rf = make_routed_field();
    const lsdtt::lsdtt_raster_meta& meta = rf.filled.get_meta();

    lsdtt::lsdtt_flowinfo_options opts;
    opts.forced_outlet_flats = rf.gt.pit_flats;

    const lsdtt::lsdtt_flowinfo flow(rf.filled.get_data_array(),
                                     rf.filled.get_mask_array(), meta, opts);

    const std::int64_t n_cells =
        static_cast<std::int64_t>(meta.n_rows) *
        static_cast<std::int64_t>(meta.n_cols);

    std::int64_t n_valid = 0;
    for (std::int64_t i = 0; i < n_cells; ++i)
    {
        if (rf.filled.get_mask_array().data()[i] != 0) ++n_valid;
    }

    std::vector<std::int64_t> tag(static_cast<std::size_t>(n_cells), -1);

    std::int64_t n_internal = 0;
    bool         disjoint   = true;

    for (std::size_t k = 0; k < rf.gt.pit_flats.size(); ++k)
    {
        const lsdtt::lsdtt_flat_index_slice up =
            flow.upslope_nodes(rf.gt.pit_flats[k]);
        n_internal += static_cast<std::int64_t>(up.size());

        for (const std::int64_t f : up)
        {
            if (tag[static_cast<std::size_t>(f)] >= 0) disjoint = false;
            tag[static_cast<std::size_t>(f)] = static_cast<std::int64_t>(k);
        }
    }

    // Everything that is not tagged must reach an outlet that is not a sink,
    // i.e. it drains off-grid.
    std::int64_t n_external = 0;
    for (const std::int64_t o : flow.outlets())
    {
        if (std::find(rf.gt.pit_flats.begin(), rf.gt.pit_flats.end(), o) !=
            rf.gt.pit_flats.end()) continue;
        n_external += static_cast<std::int64_t>(flow.upslope_nodes(o).size());
    }

    std::cout << "       " << n_internal << " internal + " << n_external
              << " external == " << n_valid << " valid\n";

    check("tag_area_conservation: internal + external == total valid, exactly",
          n_internal + n_external == n_valid,
          std::to_string(n_internal) + " + " + std::to_string(n_external) +
          " vs " + std::to_string(n_valid));

    check("tag_disjoint: no cell carries two sink ids", disjoint);
}


void test_forced_outlet_on_b_boundary()
{
    // 'b' is a forced baselevel edge, and is_forced_outlet_cell() already makes
    // every cell on it an outlet. A caller-supplied forced outlet landing on the
    // same cell must be a no-op: both paths write outlet_receiver /
    // no_direction, so the graph cannot depend on which test fired.
    const lsdtt::lsdtt_raster_cpu ramp = make_tilted_ramp(100, 100, 10.0);
    const lsdtt::lsdtt_raster_meta& meta = ramp.get_meta();

    const std::int64_t edge_flat = flat_of(meta, 0, 50);   // north edge

    lsdtt::lsdtt_flowinfo_options a;
    a.boundary_conditions = "booo";

    lsdtt::lsdtt_flowinfo_options b = a;
    b.forced_outlet_flats = { edge_flat };

    const lsdtt::lsdtt_flowinfo fa(ramp.get_data_array(), ramp.get_mask_array(),
                                   meta, a);
    const lsdtt::lsdtt_flowinfo fb(ramp.get_data_array(), ramp.get_mask_array(),
                                   meta, b);

    const std::int64_t n_cells =
        static_cast<std::int64_t>(meta.n_rows) *
        static_cast<std::int64_t>(meta.n_cols);

    bool same = true;
    for (std::int64_t i = 0; i < n_cells; ++i)
    {
        if (fa.receiver_array().data()[i] != fb.receiver_array().data()[i])
            same = false;
        if (fa.receiver_dir_array().data()[i] != fb.receiver_dir_array().data()[i])
            same = false;
    }

    check("forced_outlet_on_b_boundary: the cell is an outlet either way",
          fa.is_outlet(edge_flat) && fb.is_outlet(edge_flat));
    check("forced_outlet_on_b_boundary: the receiver graph is identical", same);
}


void test_forced_outlet_without_sink_fill()
{
    // THIS IS A TRAP FOR USERS, NOT A BUG.
    //
    // Forced outlets only do useful work on a DEM prepared with
    // fill_pits_with_sinks() using the SAME cell set. A plain fill_pits() has
    // already raised the depression to its spill level, and the min_slope
    // gradient it imposes rises INWARDS from the spill point -- so the old pit
    // is now the highest cell of the flattened depression, and nothing at all
    // drains to it. Nominating it as a forced outlet is then harmless but
    // useless: no error, and an essentially empty basin.
    //
    // The test pins that behaviour so it stays documented rather than becoming
    // a surprise. It is the documented failure mode of the option
    // (lsdtt_flow_base.hpp, forced_outlet_flats).
    const lsdtt::lsdtt_crater_field_surface gt = make_crater_field();

    const lsdtt::lsdtt_raster_cpu filled_plain =
        lsdtt::fill_pits(gt.dem, 0.0001f);

    lsdtt::lsdtt_flowinfo_options opts;
    opts.forced_outlet_flats = gt.pit_flats;

    bool        constructed = true;
    std::string what;
    std::vector<std::size_t> basin_sizes;

    try
    {
        const lsdtt::lsdtt_flowinfo flow(filled_plain.get_data_array(),
                                         filled_plain.get_mask_array(),
                                         filled_plain.get_meta(), opts);

        for (const std::int64_t f : gt.pit_flats)
            basin_sizes.push_back(flow.upslope_nodes(f).size());
    }
    catch (const std::exception& e)
    {
        constructed = false;
        what        = e.what();
    }

    check("forced_outlet_without_sink_fill: construction succeeds, no error",
          constructed, what);

    if (!constructed) return;

    // The reference: the SAME sinks on a correctly prepared surface.
    const lsdtt::lsdtt_raster_cpu filled_sinks =
        lsdtt::fill_pits_with_sinks(gt.dem, gt.pit_flats, 0.0001f);

    const lsdtt::lsdtt_flowinfo good(filled_sinks.get_data_array(),
                                     filled_sinks.get_mask_array(),
                                     filled_sinks.get_meta(), opts);

    bool        all_tiny = true;
    std::string detail;

    for (std::size_t k = 0; k < gt.pit_flats.size(); ++k)
    {
        const std::size_t got  = basin_sizes[k];
        const std::size_t want = good.upslope_nodes(gt.pit_flats[k]).size();

        const double ratio = static_cast<double>(got) /
                             static_cast<double>(want);
        if (!(ratio < 0.05)) all_tiny = false;

        detail += "sink " + std::to_string(k) + ": " + std::to_string(got) +
                  " of " + std::to_string(want) + "; ";
    }

    // Not "exactly zero": epsilon-filled flats are not monotone, so a few cells
    // adjacent to the drowned pit still happen to drain to it. The point is the
    // order of magnitude -- a couple of percent of the basin the user asked for.
    check("forced_outlet_without_sink_fill: the basins are a negligible fraction"
          " of the intended ones",
          all_tiny, detail);
}


void test_forced_outlet_validation()
{
    const routed_field rf = make_routed_field();
    const lsdtt::lsdtt_raster_meta& meta = rf.filled.get_meta();

    const std::int64_t n_cells =
        static_cast<std::int64_t>(meta.n_rows) *
        static_cast<std::int64_t>(meta.n_cols);

    check_throws(
        "forced_outlet_validation: index past the end is refused",
        [&] {
            lsdtt::lsdtt_flowinfo_options o;
            o.forced_outlet_flats = { n_cells };
            const lsdtt::lsdtt_flowinfo f(rf.filled.get_data_array(),
                                          rf.filled.get_mask_array(), meta, o);
            (void)f;
        },
        "outside the raster");

    check_throws(
        "forced_outlet_validation: negative index is refused",
        [&] {
            lsdtt::lsdtt_flowinfo_options o;
            o.forced_outlet_flats = { -1 };
            const lsdtt::lsdtt_flowinfo f(rf.filled.get_data_array(),
                                          rf.filled.get_mask_array(), meta, o);
            (void)f;
        },
        "outside the raster");

    // A forced outlet on a nodata cell has no node, so it cannot terminate
    // anything.
    lsdtt::lsdtt_raster_cpu holed = rf.filled;
    auto& hdata = holed.get_data_array();
    auto& hmask = const_cast<lsdtt::lsdtt_array2d_u8&>(holed.get_mask_array());
    hdata(20, 20) = static_cast<float>(meta.no_data_value);
    hmask(20, 20) = 0;

    check_throws(
        "forced_outlet_validation: an outlet on a nodata cell is refused",
        [&] {
            lsdtt::lsdtt_flowinfo_options o;
            o.forced_outlet_flats = { flat_of(meta, 20, 20) };
            const lsdtt::lsdtt_flowinfo f(holed.get_data_array(),
                                          holed.get_mask_array(), meta, o);
            (void)f;
        },
        "nodata or masked cell");

    // Duplicates are documented as absorbed silently: stamping the plane twice
    // is idempotent.
    lsdtt::lsdtt_flowinfo_options once;
    once.forced_outlet_flats = rf.gt.pit_flats;

    lsdtt::lsdtt_flowinfo_options twice;
    twice.forced_outlet_flats = rf.gt.pit_flats;
    twice.forced_outlet_flats.insert(twice.forced_outlet_flats.end(),
                                     rf.gt.pit_flats.begin(),
                                     rf.gt.pit_flats.end());

    const lsdtt::lsdtt_flowinfo fa(rf.filled.get_data_array(),
                                   rf.filled.get_mask_array(), meta, once);
    const lsdtt::lsdtt_flowinfo fb(rf.filled.get_data_array(),
                                   rf.filled.get_mask_array(), meta, twice);

    bool same = true;
    for (std::int64_t i = 0; i < n_cells; ++i)
    {
        if (fa.receiver_array().data()[i] != fb.receiver_array().data()[i])
            same = false;
    }

    check("forced_outlet_validation: duplicate entries are absorbed silently",
          same && fa.outlets().size() == fb.outlets().size(),
          std::to_string(fa.outlets().size()) + " vs " +
          std::to_string(fb.outlets().size()) + " outlets");
}


// =============================================================================
// Section 4: Phase 3 hole punching and sink tagging
// =============================================================================

void test_punch_sink_holes_contract()
{
    const lsdtt::lsdtt_crater_field_surface gt = make_crater_field();
    const lsdtt::lsdtt_raster_meta& meta = gt.dem.get_meta();
    const std::int64_t existing_hole = flat_of(meta, 10, 10);
    const std::int64_t edge_sink = flat_of(meta, 0, 150);

    std::vector<float> data = gt.dem.get_data();
    std::vector<std::uint8_t> mask = gt.dem.get_mask();
    data[static_cast<std::size_t>(existing_hole)] =
        static_cast<float>(meta.no_data_value);
    mask[static_cast<std::size_t>(existing_hole)] = 0u;

    const lsdtt::lsdtt_raster_cpu input(
        std::move(data), std::move(mask), meta);
    const std::vector<float> input_data = input.get_data();
    const std::vector<std::uint8_t> input_mask = input.get_mask();

    std::vector<std::int64_t> requested = gt.pit_flats;
    requested.push_back(edge_sink);
    const std::vector<lsdtt::sink_candidate> sinks =
        make_routing_sinks(requested, meta);

    const lsdtt::lsdtt_raster_cpu punched =
        lsdtt::punch_sink_holes(input, sinks);

    bool requested_are_holes = true;
    bool others_unchanged = true;
    for (std::size_t i = 0; i < input_data.size(); ++i)
    {
        const bool requested_here =
            std::find(requested.begin(), requested.end(),
                      static_cast<std::int64_t>(i)) != requested.end();
        if (requested_here)
        {
            if (punched.get_mask()[i] != 0u ||
                !meta.is_nodata(punched.get_data()[i]))
                requested_are_holes = false;
        }
        else if (punched.get_mask()[i] != input_mask[i] ||
                 punched.get_data()[i] != input_data[i])
        {
            others_unchanged = false;
        }
    }

    check("punch_sink_holes: every requested pit, including an edge cell, is nodata",
          requested_are_holes);
    check("punch_sink_holes: every other cell, including existing nodata, is unchanged",
          others_unchanged);
    check("punch_sink_holes: the input raster is not modified",
          input.get_data() == input_data && input.get_mask() == input_mask);

    const lsdtt::lsdtt_raster_cpu empty =
        lsdtt::punch_sink_holes(input, {});
    check("punch_sink_holes: empty input is an exact clone",
          empty.get_data() == input_data && empty.get_mask() == input_mask);

    lsdtt::sink_candidate rejected = sinks.front();
    rejected.accepted = false;
    check_throws("punch_sink_holes: rejected sinks are refused",
                 [&] { (void)lsdtt::punch_sink_holes(input, { rejected }); },
                 "accepted_sinks");

    lsdtt::sink_candidate masked = sinks.front();
    masked.pit_flat = existing_hole;
    masked.pit_row = 10;
    masked.pit_col = 10;
    masked.routing_sink_id = 0;
    check_throws("punch_sink_holes: a requested pre-existing hole is refused",
                 [&] { (void)lsdtt::punch_sink_holes(input, { masked }); },
                 "nodata or a masked cell");
}


void test_tag_forced_outlet_and_cutoffs()
{
    const routed_field rf = make_routed_field();
    const lsdtt::lsdtt_raster_meta& meta = rf.filled.get_meta();
    const std::vector<lsdtt::sink_candidate> sinks =
        make_routing_sinks(rf.gt.pit_flats, meta);

    lsdtt::lsdtt_flowinfo_options flow_opts;
    flow_opts.forced_outlet_flats = rf.gt.pit_flats;
    const lsdtt::lsdtt_flowinfo concrete(
        rf.filled.get_data_array(), rf.filled.get_mask_array(), meta, flow_opts);
    const lsdtt::lsdtt_flow_base& flow = concrete;

    const lsdtt::sink_tagging_result unlimited =
        lsdtt::tag_upslope_of_sinks(
            flow, sinks, lsdtt::sink_mode::forced_outlet);

    bool labels_exact = true;
    std::int64_t expected_internal = 0;
    for (const lsdtt::sink_candidate& sink : sinks)
    {
        const lsdtt::lsdtt_flat_index_slice upslope =
            flow.upslope_nodes(sink.pit_flat);
        expected_internal += static_cast<std::int64_t>(upslope.size());

        if (unlimited.n_pixels[static_cast<std::size_t>(sink.routing_sink_id)] !=
            static_cast<std::int64_t>(upslope.size()))
            labels_exact = false;
        if (unlimited.area_m2[static_cast<std::size_t>(sink.routing_sink_id)] !=
            static_cast<double>(upslope.size()) * meta.pixel_area())
            labels_exact = false;

        for (const std::int64_t flat : upslope)
        {
            const std::size_t index = static_cast<std::size_t>(flat);
            if (unlimited.sink_id_raster.get_mask()[index] == 0u ||
                unlimited.sink_id_raster.get_data()[index] !=
                    static_cast<float>(sink.routing_sink_id))
                labels_exact = false;
        }
    }

    check("tag_upslope_of_sinks: forced mode labels each exact S-vector slice",
          labels_exact);
    check("tag_upslope_of_sinks: durable ids are not used as raster labels",
          unlimited.sink_id_raster.get_data()[
              static_cast<std::size_t>(sinks.front().pit_flat)] == 0.0f &&
          sinks.front().sink_id != 0 &&
          unlimited.sink_id_raster.get_meta().no_data_value == -1.0);

    const double pixel_area = meta.pixel_area();
    const double unlimited_total =
        static_cast<double>(expected_internal) * pixel_area +
        unlimited.non_sink_outlet_area_m2 +
        unlimited.untagged_internal_area_m2;
    check("tag_upslope_of_sinks: unlimited area is conserved exactly",
          unlimited_total ==
              static_cast<double>(flow.n_valid_cells()) * pixel_area &&
          unlimited.untagged_internal_area_m2 == 0.0);

    lsdtt::sink_tagging_options finite_opts;
    finite_opts.max_upslope_distance_m = 100.0f;
    const lsdtt::sink_tagging_result finite =
        lsdtt::tag_upslope_of_sinks(
            flow, sinks, lsdtt::sink_mode::forced_outlet, finite_opts);
    const lsdtt::lsdtt_raster_result distance =
        flow.get_distance_from_outlet_raster();

    std::int64_t expected_finite = 0;
    for (const lsdtt::sink_candidate& sink : sinks)
    {
        for (const std::int64_t flat : flow.upslope_nodes(sink.pit_flat))
        {
            if (distance.data[static_cast<std::size_t>(flat)] < 100.0f)
                ++expected_finite;
        }
    }

    std::int64_t measured_finite = 0;
    for (const std::int64_t count : finite.n_pixels) measured_finite += count;
    check("tag_upslope_of_sinks: finite cutoff uses strict flow distance",
          measured_finite == expected_finite &&
          finite.untagged_internal_area_m2 ==
              static_cast<double>(expected_internal - expected_finite) * pixel_area);
    check("tag_upslope_of_sinks: cutoff does not change non-sink outlet area",
          finite.non_sink_outlet_area_m2 == unlimited.non_sink_outlet_area_m2);

    lsdtt::sink_tagging_options zero_opts;
    zero_opts.max_upslope_distance_m = 0.0f;
    const lsdtt::sink_tagging_result zero =
        lsdtt::tag_upslope_of_sinks(
            flow, sinks, lsdtt::sink_mode::forced_outlet, zero_opts);
    std::int64_t zero_count = 0;
    for (const std::int64_t count : zero.n_pixels) zero_count += count;
    check("tag_upslope_of_sinks: zero cutoff tags nothing",
          zero_count == 0 &&
          zero.untagged_internal_area_m2 ==
              static_cast<double>(expected_internal) * pixel_area);
}


void test_nodata_hole_equivalence()
{
    const lsdtt::lsdtt_crater_field_surface gt = make_crater_field();
    const lsdtt::lsdtt_raster_meta& meta = gt.dem.get_meta();
    const std::vector<lsdtt::sink_candidate> sinks =
        make_routing_sinks(gt.pit_flats, meta);

    const lsdtt::lsdtt_raster_cpu forced_filled =
        lsdtt::fill_pits_with_sinks(gt.dem, gt.pit_flats, 0.0001f);
    lsdtt::lsdtt_flowinfo_options forced_opts;
    forced_opts.forced_outlet_flats = gt.pit_flats;
    const lsdtt::lsdtt_flowinfo forced_flow(
        forced_filled.get_data_array(), forced_filled.get_mask_array(),
        meta, forced_opts);
    const lsdtt::sink_tagging_result forced =
        lsdtt::tag_upslope_of_sinks(
            forced_flow, sinks, lsdtt::sink_mode::forced_outlet);

    const lsdtt::lsdtt_raster_cpu punched =
        lsdtt::punch_sink_holes(gt.dem, sinks);
    const lsdtt::lsdtt_raster_cpu hole_filled =
        lsdtt::fill_pits_priority_flood(punched, 0.0001f);
    const lsdtt::lsdtt_flowinfo hole_flow(
        hole_filled.get_data_array(), hole_filled.get_mask_array(), meta, {});
    const lsdtt::sink_tagging_result holes =
        lsdtt::tag_upslope_of_sinks(
            hole_flow, sinks, lsdtt::sink_mode::nodata_hole);

    bool common_labels_equal = true;
    for (std::size_t i = 0; i < punched.full_grid_size(); ++i)
    {
        if (punched.get_mask()[i] == 0u) continue;
        if (forced.sink_id_raster.get_mask()[i] !=
                holes.sink_id_raster.get_mask()[i] ||
            (forced.sink_id_raster.get_mask()[i] != 0u &&
             forced.sink_id_raster.get_data()[i] !=
                holes.sink_id_raster.get_data()[i]))
            common_labels_equal = false;
    }

    check("nodata_hole_equivalence: labels match on the common valid domain",
          common_labels_equal);

    bool pit_policy_correct = true;
    for (const lsdtt::sink_candidate& sink : sinks)
    {
        const std::size_t index = static_cast<std::size_t>(sink.pit_flat);
        if (forced.sink_id_raster.get_mask()[index] == 0u ||
            holes.sink_id_raster.get_mask()[index] != 0u)
            pit_policy_correct = false;
    }
    check("nodata_hole_equivalence: forced pits are labelled and holes remain nodata",
          pit_policy_correct);

    double hole_internal = holes.untagged_internal_area_m2;
    for (const double area : holes.area_m2) hole_internal += area;
    check("nodata_hole_equivalence: hole-mode area is conserved",
          hole_internal + holes.non_sink_outlet_area_m2 ==
              static_cast<double>(hole_flow.n_valid_cells()) * meta.pixel_area());

    lsdtt::sink_tagging_options cutoff;
    cutoff.max_upslope_distance_m = 100.0f;
    const lsdtt::sink_tagging_result holes_finite =
        lsdtt::tag_upslope_of_sinks(
            hole_flow, sinks, lsdtt::sink_mode::nodata_hole, cutoff);

    const lsdtt::lsdtt_raster_result hole_distance =
        hole_flow.get_distance_from_outlet_raster();
    std::int64_t expected_hole_finite = 0;
    for (const lsdtt::sink_candidate& sink : sinks)
    {
        for (int direction = 0; direction < 8; ++direction)
        {
            const int row = sink.pit_row + k_d_row[direction];
            const int col = sink.pit_col + k_d_col[direction];
            if (!hole_flow.in_bounds(row, col) ||
                !hole_flow.is_valid_cell(row, col)) continue;

            const std::int64_t seed = hole_flow.flat_index(row, col);
            if (!hole_flow.is_outlet(seed)) continue;

            const double final_step = meta.cell_size *
                ((k_d_row[direction] != 0 && k_d_col[direction] != 0)
                     ? std::sqrt(2.0) : 1.0);
            for (const std::int64_t flat : hole_flow.upslope_nodes(seed))
            {
                if (static_cast<double>(
                        hole_distance.data[static_cast<std::size_t>(flat)]) +
                        final_step < 100.0)
                    ++expected_hole_finite;
            }
        }
    }
    std::int64_t measured_hole_finite = 0;
    for (const std::int64_t count : holes_finite.n_pixels)
        measured_hole_finite += count;
    check("nodata_hole_equivalence: finite cutoff includes the final step into the hole",
          measured_hole_finite == expected_hole_finite);
}


void test_phase3_validation()
{
    const routed_field rf = make_routed_field();
    const lsdtt::lsdtt_raster_meta& meta = rf.filled.get_meta();
    const std::vector<lsdtt::sink_candidate> sinks =
        make_routing_sinks(rf.gt.pit_flats, meta);

    lsdtt::lsdtt_flowinfo_options opts;
    opts.forced_outlet_flats = rf.gt.pit_flats;
    const lsdtt::lsdtt_flowinfo flow(
        rf.filled.get_data_array(), rf.filled.get_mask_array(), meta, opts);

    const lsdtt::sink_tagging_result empty =
        lsdtt::tag_upslope_of_sinks(
            flow, {}, lsdtt::sink_mode::forced_outlet);
    bool empty_mask = true;
    for (const std::uint8_t value : empty.sink_id_raster.get_mask())
    {
        if (value != 0u) empty_mask = false;
    }
    check("phase3_validation: empty sinks classify all valid area as non-sink",
          empty.non_sink_outlet_area_m2 ==
              static_cast<double>(flow.n_valid_cells()) * meta.pixel_area() &&
          empty.untagged_internal_area_m2 == 0.0 && empty_mask);

    lsdtt::sink_tagging_options bad_cutoff;
    bad_cutoff.max_upslope_distance_m = -1.0f;
    check_throws("phase3_validation: negative cutoff is refused",
                 [&] {
                     (void)lsdtt::tag_upslope_of_sinks(
                         flow, sinks, lsdtt::sink_mode::forced_outlet, bad_cutoff);
                 },
                 "non-negative");
    bad_cutoff.max_upslope_distance_m =
        std::numeric_limits<float>::quiet_NaN();
    check_throws("phase3_validation: NaN cutoff is refused",
                 [&] {
                     (void)lsdtt::tag_upslope_of_sinks(
                         flow, sinks, lsdtt::sink_mode::forced_outlet, bad_cutoff);
                 },
                 "non-negative");

    std::vector<lsdtt::sink_candidate> duplicate = sinks;
    duplicate[1].pit_flat = duplicate[0].pit_flat;
    duplicate[1].pit_row = duplicate[0].pit_row;
    duplicate[1].pit_col = duplicate[0].pit_col;
    check_throws("phase3_validation: duplicate pit cells are refused",
                 [&] {
                     (void)lsdtt::tag_upslope_of_sinks(
                         flow, duplicate, lsdtt::sink_mode::forced_outlet);
                 },
                 "duplicate pit_flat");

    std::vector<lsdtt::sink_candidate> sparse = sinks;
    sparse[0].routing_sink_id = 99;
    check_throws("phase3_validation: sparse routing ids are refused",
                 [&] {
                     (void)lsdtt::tag_upslope_of_sinks(
                         flow, sparse, lsdtt::sink_mode::forced_outlet);
                 },
                 "contiguous from zero");

    std::vector<lsdtt::sink_candidate> rejected = sinks;
    rejected[0].accepted = false;
    check_throws("phase3_validation: rejected sinks are refused",
                 [&] {
                     (void)lsdtt::tag_upslope_of_sinks(
                         flow, rejected, lsdtt::sink_mode::forced_outlet);
                 },
                 "accepted_sinks");

    check_throws("phase3_validation: nodata mode refuses a valid pit",
                 [&] {
                     (void)lsdtt::tag_upslope_of_sinks(
                         flow, sinks, lsdtt::sink_mode::nodata_hole);
                 },
                 "still valid");

    lsdtt::lsdtt_flowinfo_options no_s_opts = opts;
    no_s_opts.build_s_vector = false;
    const lsdtt::lsdtt_flowinfo no_s_vector(
        rf.filled.get_data_array(), rf.filled.get_mask_array(), meta, no_s_opts);
    check_throws("phase3_validation: an absent S-vector is refused",
                 [&] {
                     (void)lsdtt::tag_upslope_of_sinks(
                         no_s_vector, sinks, lsdtt::sink_mode::forced_outlet);
                 },
                 "no s_vector");

    std::int64_t non_outlet = -1;
    for (std::int64_t flat = 0;
         flat < static_cast<std::int64_t>(flow.full_grid_size()); ++flat)
    {
        if (flow.is_valid_cell(flat) && !flow.is_outlet(flat))
        {
            non_outlet = flat;
            break;
        }
    }
    std::vector<lsdtt::sink_candidate> wrong_outlet =
        make_routing_sinks({ non_outlet }, meta);
    check_throws("phase3_validation: forced mode requires an outlet",
                 [&] {
                     (void)lsdtt::tag_upslope_of_sinks(
                         flow, wrong_outlet, lsdtt::sink_mode::forced_outlet);
                 },
                 "not a valid outlet");

    // Two adjacent holes can share an outlet neighbour. There is no flow-graph
    // information that says which hole owns that outlet, so iteration order
    // must not decide the label.
    lsdtt::lsdtt_raster_meta small_meta;
    small_meta.n_rows = 20;
    small_meta.n_cols = 20;
    small_meta.cell_size = 10.0;
    small_meta.set_default_synthetic_georeference();
    const lsdtt::lsdtt_crater_field_surface flat_surface =
        lsdtt::generate_crater_field_surface(
            small_meta, std::vector<lsdtt::lsdtt_crater_spec>{}, 0.0, 1000.0);
    const std::vector<std::int64_t> adjacent_flats = {
        flat_of(small_meta, 10, 9), flat_of(small_meta, 10, 10)
    };
    const std::vector<lsdtt::sink_candidate> adjacent_sinks =
        make_routing_sinks(adjacent_flats, small_meta);
    const lsdtt::lsdtt_raster_cpu adjacent_punched =
        lsdtt::punch_sink_holes(flat_surface.dem, adjacent_sinks);
    const lsdtt::lsdtt_raster_cpu adjacent_filled =
        lsdtt::fill_pits_priority_flood(adjacent_punched, 0.0001f);
    const lsdtt::lsdtt_flowinfo adjacent_flow(
        adjacent_filled.get_data_array(), adjacent_filled.get_mask_array(),
        small_meta, {});
    check_throws("phase3_validation: a shared nodata-hole outlet is refused",
                 [&] {
                     (void)lsdtt::tag_upslope_of_sinks(
                         adjacent_flow, adjacent_sinks,
                         lsdtt::sink_mode::nodata_hole);
                 },
                 "attribution is ambiguous");

    const std::vector<lsdtt::sink_candidate> edge_hole =
        make_routing_sinks({ flat_of(small_meta, 0, 10) }, small_meta);
    const lsdtt::lsdtt_raster_cpu edge_punched =
        lsdtt::punch_sink_holes(flat_surface.dem, edge_hole);
    const lsdtt::lsdtt_raster_cpu edge_filled =
        lsdtt::fill_pits_priority_flood(edge_punched, 0.0001f);
    const lsdtt::lsdtt_flowinfo edge_flow(
        edge_filled.get_data_array(), edge_filled.get_mask_array(),
        small_meta, {});
    check_throws("phase3_validation: an edge hole is attribution-ambiguous",
                 [&] {
                     (void)lsdtt::tag_upslope_of_sinks(
                         edge_flow, edge_hole, lsdtt::sink_mode::nodata_hole);
                 },
                 "raster edge");

    std::vector<float> void_data = flat_surface.dem.get_data();
    std::vector<std::uint8_t> void_mask = flat_surface.dem.get_mask();
    const std::int64_t unrelated_void = flat_of(small_meta, 10, 12);
    void_data[static_cast<std::size_t>(unrelated_void)] =
        static_cast<float>(small_meta.no_data_value);
    void_mask[static_cast<std::size_t>(unrelated_void)] = 0u;
    const lsdtt::lsdtt_raster_cpu with_void(
        std::move(void_data), std::move(void_mask), small_meta);
    const std::vector<lsdtt::sink_candidate> nearby_hole =
        make_routing_sinks({ flat_of(small_meta, 10, 10) }, small_meta);
    const lsdtt::lsdtt_raster_cpu nearby_punched =
        lsdtt::punch_sink_holes(with_void, nearby_hole);
    const lsdtt::lsdtt_raster_cpu nearby_filled =
        lsdtt::fill_pits_priority_flood(nearby_punched, 0.0001f);
    const lsdtt::lsdtt_flowinfo nearby_flow(
        nearby_filled.get_data_array(), nearby_filled.get_mask_array(),
        small_meta, {});
    check_throws("phase3_validation: a hole beside unrelated nodata is ambiguous",
                 [&] {
                     (void)lsdtt::tag_upslope_of_sinks(
                         nearby_flow, nearby_hole,
                         lsdtt::sink_mode::nodata_hole);
                 },
                 "another nodata cell");
}

} // namespace


int main()
{
    std::cout << "lsdtt3 internal-drainage unit tests "
                 "(phases 0, 1 and 3)\n\n";

    std::cout << "--- crater generator: analytic oracle ---\n";
    test_single_crater_pit_location();
    test_single_crater_geometry();
    test_flat_floor_tie_break();
    test_generator_validation();

    std::cout << "\n--- fill_pits_with_sinks ---\n";
    test_fill_with_sinks_retains_descent();
    test_sinks_never_raised();
    test_empty_sink_set_equivalence();
    test_fill_with_sinks_validation();

    std::cout << "\n--- forced_outlet_flats ---\n";
    test_forced_outlet_terminates();
    test_forced_outlet_beats_steepest_descent();
    test_empty_forced_outlet_set_equivalence();
    test_tag_area_conservation_and_disjointness();
    test_forced_outlet_on_b_boundary();
    test_forced_outlet_without_sink_fill();
    test_forced_outlet_validation();

    std::cout << "\n--- Phase 3: hole punching and sink tagging ---\n";
    test_punch_sink_holes_contract();
    test_tag_forced_outlet_and_cutoffs();
    test_nodata_hole_equivalence();
    test_phase3_validation();

    std::cout << "\n" << g_n_checks - g_n_failed << " / " << g_n_checks
              << " checks passed\n";
    std::cout << "RESULT: " << (g_n_failed == 0 ? "PASS" : "FAIL") << "\n";

    return (g_n_failed == 0) ? 0 : 1;
}
