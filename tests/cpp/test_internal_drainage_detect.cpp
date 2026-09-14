// test_internal_drainage_detect.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
//
// -----------------------------------------------------------------------------
// Unit tests for phase 2b of the internal-drainage work
// (docs/plan_docs/internal_drainage_plan.md, sections 4.3, 5 and 6.1): the
// detection module src/cpu/lsdtt_internal_drainage_cpu.
//
// Phases 0 and 1 -- the crater generator, fill_pits_with_sinks and
// forced_outlet_flats -- are covered by tests/cpp/test_internal_drainage.cpp
// and are not repeated here.  Nothing in this driver routes anything; routing
// and tagging are phase 3.
//
// What is under test
// ------------------
//   find_internal_drainage(), method depression_hierarchy
//     pits against the generator's analytic oracle; depth / spill / area /
//     volume against the same oracle to one boundary cell; the deterministic
//     flat-floor tie-break; the BINARY-TREE node count of a nested group
//     (leaves + (leaves - 1) metas, so three nested hollows give five sinks,
//     not four); the fact that a meta-depression shares its pit CELL with its
//     lowest child, so a pit cell is not a unique key; the four size filters;
//     and the reattachment of a filtered leaf's cells to its nearest surviving
//     ancestor
//
//   find_internal_drainage(), method fill_difference (the legacy path)
//     agreement with the hierarchy on a simple crater field, and the
//     DOCUMENTED disagreement in COUNT on a nested one (plan section 5.4)
//
//   the review round-trip
//     write -> read -> write in CSV, FlatGeobuf and GeoPackage; .csv and .fgb
//     from one detection run parsing back to identical records; attribute
//     typing; the four position tiers; a hand-added feature; a candidate layer
//     supplied in EPSG:4326; format resolution errors; and -- plan 5.7 -- a
//     FlatGeobuf round trip with the SPATIAL INDEX ON, which reorders features
//     on write and must therefore be keyed on sink_id and never on file
//     position
//
//   sink_depressions_to_polygons()
//     the 1:1 join to the point layer in the leaves-only case, and the fact
//     that with leaves_only == false an intermediate node whose cells are all
//     claimed by its children has an EMPTY footprint and therefore no polygon
//
// A trap worth restating (plan section 5.3b): sink_candidate::area_m2 is the
// DEPRESSION FOOTPRINT, not the drainage basin of the sink.  Nothing here
// asserts that the two are equal, because they are not.
//
// All fixtures are generated in process; the vector and CSV files are written
// at run time into the directory given as argv[1] (default
// ./internal_drainage_detect_output), so nothing is committed to tests/data.
//
// Usage: lsdtt-test-internal-drainage-detect [output_dir]
// Prints one PASS/FAIL line per check and a final "RESULT: PASS|FAIL".
// -----------------------------------------------------------------------------

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <gdal.h>
#include <gdal_priv.h>
#include <ogr_geometry.h>
#include <ogr_spatialref.h>
#include <ogrsf_frmts.h>

#include "lsdtt_data_loader.hpp"
#include "lsdtt_internal_drainage_cpu.hpp"
#include "lsdtt_raster_cpu.hpp"
#include "lsdtt_raster_synthetic_surfaces_cpu.hpp"

namespace
{

// =============================================================================
// Tiny check harness (same shape as tests/cpp/test_internal_drainage.cpp)
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

// Reports a measurement that is deliberately NOT asserted -- a documented
// divergence or an environment-dependent observation.
void report(const std::string& name, const std::string& detail)
{
    std::cout << "NOTE  " << name << "   [" << detail << "]" << std::endl;
}

template <typename F>
void check_throws(const std::string& name, const F& f, const std::string& needle)
{
    bool        threw   = false;
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

// Signed shoelace area of a ring.
double ring_area(const std::vector<lsdtt::lsdtt_xy>& ring)
{
    if (ring.size() < 3) return 0.0;

    double twice = 0.0;
    for (std::size_t i = 0; i < ring.size(); ++i)
    {
        const lsdtt::lsdtt_xy& a = ring[i];
        const lsdtt::lsdtt_xy& b = ring[(i + 1) % ring.size()];
        twice += a.x * b.y - b.x * a.y;
    }

    return 0.5 * twice;
}

// Net planform area of a polygon: exterior less its holes.
double polygon_area(const lsdtt::lsdtt_polygon_geom& polygon)
{
    double area = std::fabs(ring_area(polygon.exterior));

    for (const std::vector<lsdtt::lsdtt_xy>& hole : polygon.holes)
    {
        area -= std::fabs(ring_area(hole));
    }

    return area;
}

bool ring_contains(const std::vector<lsdtt::lsdtt_xy>& ring, double x, double y)
{
    bool inside = false;

    for (std::size_t i = 0, j = ring.size() - 1; i < ring.size(); j = i++)
    {
        const double xi = ring[i].x, yi = ring[i].y;
        const double xj = ring[j].x, yj = ring[j].y;

        if (((yi > y) != (yj > y)) &&
            (x < (xj - xi) * (y - yi) / (yj - yi) + xi))
        {
            inside = !inside;
        }
    }

    return inside;
}

// Point in polygon, holes honoured.
bool polygon_contains(const lsdtt::lsdtt_polygon_geom& polygon,
                      double x,
                      double y)
{
    if (polygon.exterior.size() < 3)              return false;
    if (!ring_contains(polygon.exterior, x, y))   return false;

    for (const std::vector<lsdtt::lsdtt_xy>& hole : polygon.holes)
    {
        if (hole.size() >= 3 && ring_contains(hole, x, y)) return false;
    }

    return true;
}

// =============================================================================
// Fixtures
//
// simple crater field  three non-overlapping cones on a tilted plane: three
//                      independent leaf depressions, no nesting anywhere
//
// nested crater field  a flat-floored bowl with three cones cut into its floor.
//                      The generator refuses overlapping rims -- deliberately,
//                      because two merged craters have a pit and a spill that
//                      neither spec describes -- so the nested surface is built
//                      by SUPERPOSITION of two generated surfaces: the bowl,
//                      plus a zero-background field of three small cones whose
//                      values are <= 0.  Both halves are still the generator's
//                      own analytic surface, and the sum is exact.
// =============================================================================
struct crater_field_fixture
{
    lsdtt::lsdtt_crater_field_surface gt;
    std::vector<int>                  rows;
    std::vector<int>                  cols;
    std::vector<double>               radii;
    std::vector<double>               depths;
    double                            background_slope = 0.02;
};

lsdtt::lsdtt_raster_meta make_meta(int n_rows, int n_cols, double cell_size)
{
    lsdtt::lsdtt_raster_meta meta;
    meta.n_rows    = n_rows;
    meta.n_cols    = n_cols;
    meta.cell_size = cell_size;
    meta.x_min     = 0.0;
    meta.y_min     = 0.0;
    meta.set_default_synthetic_georeference();
    return meta;
}

crater_field_fixture make_simple_crater_field()
{
    crater_field_fixture fx;

    const lsdtt::lsdtt_raster_meta meta = make_meta(300, 300, 10.0);

    fx.rows   = { 100, 200,  60 };
    fx.cols   = { 100, 210, 230 };
    fx.radii  = { 300.0, 200.0, 150.0 };
    fx.depths = { 100.0,  80.0,  60.0 };

    std::vector<lsdtt::lsdtt_crater_spec> specs;
    for (std::size_t i = 0; i < fx.rows.size(); ++i)
    {
        lsdtt::lsdtt_crater_spec cs;
        cs.centre_x       = pixel_centre_x(meta, fx.cols[i]);
        cs.centre_y       = pixel_centre_y(meta, fx.rows[i]);
        cs.radius_m       = fx.radii[i];
        cs.depth_m        = fx.depths[i];
        cs.floor_fraction = 0.0;
        specs.push_back(cs);
    }

    fx.gt = lsdtt::generate_crater_field_surface(
        meta, specs, fx.background_slope, 1000.0);

    return fx;
}

// The nested fixture: one flat-floored bowl containing three cones.
struct nested_field_fixture
{
    lsdtt::lsdtt_raster_cpu   dem;
    std::vector<std::int64_t> cone_pits;   // the three analytic leaf pits
    std::int64_t              bowl_pit = -1;
};

nested_field_fixture make_nested_crater_field()
{
    nested_field_fixture fx;

    const lsdtt::lsdtt_raster_meta meta = make_meta(300, 300, 10.0);

    // The bowl: a wide crater with a large exactly flat floor, cut into a
    // tilted plane.  Its floor is where the hollows will go.
    lsdtt::lsdtt_crater_spec bowl;
    bowl.centre_x       = pixel_centre_x(meta, 150);
    bowl.centre_y       = pixel_centre_y(meta, 150);
    bowl.radius_m       = 1200.0;
    bowl.depth_m        = 200.0;
    bowl.floor_fraction = 0.6;            // a 720 m radius flat floor

    const lsdtt::lsdtt_crater_field_surface bowl_surface =
        lsdtt::generate_crater_field_surface(meta, { bowl }, 0.02, 1000.0);

    fx.bowl_pit = bowl_surface.pit_flats[0];

    // Three small cones on a ZERO background: outside their rims the surface is
    // exactly zero, inside it dips to -depth.  Adding this to the bowl cuts
    // three hollows into the bowl's flat floor and changes nothing outside
    // them.
    // Distinct radii as well as distinct depths, so that a size filter can
    // single out one hollow.
    const int    cone_rows[3]   = { 120, 175, 130 };
    const int    cone_cols[3]   = { 120, 140, 185 };
    const double cone_radius[3] = { 200.0, 170.0, 140.0 };
    const double cone_depth[3]  = {  30.0,  20.0,  10.0 };

    std::vector<lsdtt::lsdtt_crater_spec> cone_specs;
    for (int i = 0; i < 3; ++i)
    {
        lsdtt::lsdtt_crater_spec cs;
        cs.centre_x       = pixel_centre_x(meta, cone_cols[i]);
        cs.centre_y       = pixel_centre_y(meta, cone_rows[i]);
        cs.radius_m       = cone_radius[i];
        cs.depth_m        = cone_depth[i];
        cs.floor_fraction = 0.0;
        cone_specs.push_back(cs);
    }

    const lsdtt::lsdtt_crater_field_surface cone_surface =
        lsdtt::generate_crater_field_surface(meta, cone_specs, 0.0, 0.0);

    fx.cone_pits = cone_surface.pit_flats;

    fx.dem = bowl_surface.dem;

    auto&       z    = fx.dem.get_data_array();
    const auto& cone = cone_surface.dem.get_data_array();

    const std::int64_t n_cells =
        static_cast<std::int64_t>(meta.n_rows) *
        static_cast<std::int64_t>(meta.n_cols);

    for (std::int64_t i = 0; i < n_cells; ++i)
    {
        z.data()[i] = z.data()[i] + cone.data()[i];
    }

    return fx;
}

// =============================================================================
// Small utilities over a sink list
// =============================================================================
std::vector<std::int64_t> pit_flats_of(
    const std::vector<lsdtt::sink_candidate>& sinks)
{
    std::vector<std::int64_t> flats;
    flats.reserve(sinks.size());
    for (const lsdtt::sink_candidate& s : sinks) flats.push_back(s.pit_flat);
    std::sort(flats.begin(), flats.end());
    return flats;
}

const lsdtt::sink_candidate* sink_at_pit(
    const std::vector<lsdtt::sink_candidate>& sinks,
    std::int64_t                              pit_flat)
{
    for (const lsdtt::sink_candidate& s : sinks)
    {
        if (s.pit_flat == pit_flat) return &s;
    }
    return nullptr;
}

// Field-for-field equality of two candidate records.  Every field is compared:
// this is what "identical sink_candidate vectors" in plan 6.1 means, and it is
// what makes the cross-format tests worth anything.
bool records_equal(const lsdtt::sink_candidate& a, const lsdtt::sink_candidate& b)
{
    return a.sink_id         == b.sink_id         &&
           a.routing_sink_id == b.routing_sink_id &&
           a.pit_flat        == b.pit_flat        &&
           a.pit_row         == b.pit_row         &&
           a.pit_col         == b.pit_col         &&
           a.pit_x           == b.pit_x           &&
           a.pit_y           == b.pit_y           &&
           a.latitude        == b.latitude        &&
           a.longitude       == b.longitude       &&
           a.pit_elevation   == b.pit_elevation   &&
           a.spill_elevation == b.spill_elevation &&
           a.depth_m         == b.depth_m         &&
           a.area_m2         == b.area_m2         &&
           a.volume_m3       == b.volume_m3       &&
           a.n_pixels        == b.n_pixels        &&
           a.tagged_drainage_area_m2 == b.tagged_drainage_area_m2 &&
           a.tagged_n_pixels == b.tagged_n_pixels &&
           a.hierarchy_level == b.hierarchy_level &&
           a.parent_sink_id  == b.parent_sink_id  &&
           a.accepted        == b.accepted;
}

std::string describe(const lsdtt::sink_candidate& s)
{
    std::ostringstream out;
    // Full precision: a cross-format difference is often in the last digit, and
    // a failure message that rounds it away is useless.
    out.precision(17);
    out << "id " << s.sink_id
        << " routing_id " << s.routing_sink_id
        << " flat " << s.pit_flat
        << " z " << s.pit_elevation
        << " spill " << s.spill_elevation
        << " depth " << s.depth_m
        << " area " << s.area_m2
        << " vol " << s.volume_m3
        << " n " << s.n_pixels
        << " level " << s.hierarchy_level
        << " parent " << s.parent_sink_id
        << " accepted " << (s.accepted ? 1 : 0);
    return out.str();
}

bool vectors_equal(const std::vector<lsdtt::sink_candidate>& a,
                   const std::vector<lsdtt::sink_candidate>& b,
                   std::string&                              first_difference)
{
    if (a.size() != b.size())
    {
        first_difference = "sizes " + std::to_string(a.size()) + " vs " +
                           std::to_string(b.size());
        return false;
    }

    for (std::size_t i = 0; i < a.size(); ++i)
    {
        if (!records_equal(a[i], b[i]))
        {
            first_difference = "record " + std::to_string(i) + ": " +
                               describe(a[i]) + " vs " + describe(b[i]);
            return false;
        }
    }

    return true;
}

std::string read_whole_file(const std::string& path)
{
    std::ifstream in(path, std::ios::binary);
    std::ostringstream buffer;
    buffer << in.rdbuf();
    return buffer.str();
}

// =============================================================================
// Section A: detection with the depression hierarchy
// =============================================================================

void test_detected_pits_match_the_oracle()
{
    const crater_field_fixture fx = make_simple_crater_field();

    lsdtt::sink_detection_options opts;   // defaults: hierarchy, leaves only

    const std::vector<lsdtt::sink_candidate> sinks =
        lsdtt::find_internal_drainage(fx.gt.dem, opts);

    check("detected_pits_match_the_oracle: one sink per crater",
          sinks.size() == fx.gt.pit_flats.size(),
          std::to_string(sinks.size()) + " sinks, " +
          std::to_string(fx.gt.pit_flats.size()) + " craters");

    std::vector<std::int64_t> expected = fx.gt.pit_flats;
    std::sort(expected.begin(), expected.end());

    check("detected_pits_match_the_oracle: pit cells are exactly the analytic pits",
          pit_flats_of(sinks) == expected);

    // pit_row / pit_col / pit_flat must agree with each other and with the
    // elevation the DEM actually holds there: the record is read by a human in
    // a GIS, and a row/col that disagrees with the geometry is silently wrong.
    const lsdtt::lsdtt_raster_meta& meta = fx.gt.dem.get_meta();
    const auto& z = fx.gt.dem.get_data_array();

    bool consistent = true;
    for (const lsdtt::sink_candidate& s : sinks)
    {
        if (flat_of(meta, s.pit_row, s.pit_col) != s.pit_flat) consistent = false;
        if (z(s.pit_row, s.pit_col) != s.pit_elevation)        consistent = false;
        if (s.pit_x != pixel_centre_x(meta, s.pit_col))        consistent = false;
        if (s.pit_y != pixel_centre_y(meta, s.pit_row))        consistent = false;
    }

    check("detected_pits_match_the_oracle: row/col/flat/x/y/elevation agree",
          consistent);

    // Every leaf on a field of three separate craters, so no nesting is
    // reported and nothing is a parent of anything.  Three craters 1.5 km
    // apart cannot contain one another, and a parent_sink_id that says
    // otherwise would make the emitted tree -- which the polygons, the stats
    // and any downstream area accounting read -- wrong.
    //
    // MEASURED 2026-08-18, and reported as a defect rather than accommodated
    // here: the third crater comes back with parent_sink_id 0.
    // build_depression_hierarchy sets depression_node::parent for an
    // OCEAN LINK as well as for containment (see ocean_linked in
    // lsdtt_depression_hierarchy_cpu.hpp: "depressions that link to the ocean
    // through this one but are NOT sub-depressions"), and
    // nearest_selected_ancestor() in lsdtt_internal_drainage_cpu.cpp walks
    // .parent without testing .ocean_parent.  The probe reads
    //   d1 parent 0 ocean_parent 1 cells 1151
    //   d3 parent 1 ocean_parent 1 cells 657
    // so d3 is reported as nested inside d1.  sink_depressions_to_polygons()
    // is unaffected -- its cell walk does test ocean_parent -- so the symptom
    // is confined to parent_sink_id.
    bool all_leaves = true;
    for (const lsdtt::sink_candidate& s : sinks)
    {
        if (s.hierarchy_level != 0 || s.parent_sink_id != -1) all_leaves = false;
    }

    std::ostringstream tree_detail;
    for (const lsdtt::sink_candidate& s : sinks)
    {
        tree_detail << "level " << s.hierarchy_level
                    << " parent " << s.parent_sink_id << "; ";
    }

    check("detected_pits_match_the_oracle: separate craters are all top-level leaves",
          all_leaves, tree_detail.str());

    // Latitude and longitude are filled from the raster CRS through the OGR
    // transform, not left at zero.  The synthetic georeference is UTM 58S.
    bool lonlat_set = true;
    for (const lsdtt::sink_candidate& s : sinks)
    {
        if (s.latitude == 0.0 || s.longitude == 0.0) lonlat_set = false;
    }

    check("detected_pits_match_the_oracle: latitude/longitude are populated",
          lonlat_set,
          "first sink lat " + std::to_string(sinks[0].latitude) +
          " lon " + std::to_string(sinks[0].longitude));
}


void test_detected_geometry_matches_the_oracle()
{
    const crater_field_fixture fx = make_simple_crater_field();
    const lsdtt::lsdtt_raster_meta& meta = fx.gt.dem.get_meta();

    const std::vector<lsdtt::sink_candidate> sinks =
        lsdtt::find_internal_drainage(fx.gt.dem, lsdtt::sink_detection_options{});

    for (std::size_t c = 0; c < fx.gt.pit_flats.size(); ++c)
    {
        const std::string tag = "crater " + std::to_string(c);

        const lsdtt::sink_candidate* s = sink_at_pit(sinks, fx.gt.pit_flats[c]);

        if (s == nullptr)
        {
            check("detected_geometry_matches_the_oracle: " + tag + " found",
                  false);
            continue;
        }

        // The grid resolves the rim only to the nearest cell, so the measured
        // spill sits above the analytic minimum of the continuous rim by at
        // most one cell diagonal of rim relief.
        const double rim_gradient = fx.depths[c] / fx.radii[c] + fx.background_slope;
        const double spill_tol    = rim_gradient * meta.cell_size * std::sqrt(2.0);

        check("detected_geometry_matches_the_oracle: " + tag + " spill elevation",
              static_cast<double>(s->spill_elevation) >=
                  fx.gt.spill_elevation[c] - 1.0e-3 &&
              static_cast<double>(s->spill_elevation) <=
                  fx.gt.spill_elevation[c] + spill_tol,
              "detected " + std::to_string(s->spill_elevation) +
              " vs analytic " + std::to_string(fx.gt.spill_elevation[c]) +
              " (tol +" + std::to_string(spill_tol) + ")");

        const double analytic_depth =
            fx.gt.spill_elevation[c] - fx.gt.pit_elevation[c];

        check("detected_geometry_matches_the_oracle: " + tag + " depth",
              std::fabs(static_cast<double>(s->depth_m) - analytic_depth) <
                  spill_tol + 1.0e-3,
              "detected " + std::to_string(s->depth_m) +
              " vs analytic " + std::to_string(analytic_depth));

        // depth_m is defined as spill minus pit and must be exactly that, not
        // an independently measured relief.
        check("detected_geometry_matches_the_oracle: " + tag +
              " depth is spill minus pit",
              s->depth_m == s->spill_elevation - s->pit_elevation);

        // Discretisation error is one cell along the whole depression
        // boundary, i.e. of order 2 * pi * R * cell_size in area.
        const double area_tol = 2.0 * M_PI * fx.radii[c] * meta.cell_size;

        check("detected_geometry_matches_the_oracle: " + tag + " area",
              std::fabs(s->area_m2 - fx.gt.area_m2[c]) < area_tol,
              "detected " + std::to_string(s->area_m2) +
              " vs analytic " + std::to_string(fx.gt.area_m2[c]) +
              " (tol " + std::to_string(area_tol) + ")");

        check("detected_geometry_matches_the_oracle: " + tag +
              " area is n_pixels times the cell area",
              std::fabs(s->area_m2 -
                        static_cast<double>(s->n_pixels) * meta.pixel_area()) <
                  1.0e-6);

        // The same boundary strip in volume, weighted by the mean depth there.
        const double vol_tol = 3.0 * (meta.cell_size / fx.radii[c]) *
                               fx.gt.volume_m3[c];

        check("detected_geometry_matches_the_oracle: " + tag + " volume",
              std::fabs(s->volume_m3 - fx.gt.volume_m3[c]) < vol_tol,
              "detected " + std::to_string(s->volume_m3) +
              " vs analytic " + std::to_string(fx.gt.volume_m3[c]) +
              " (tol " + std::to_string(vol_tol) + ")");
    }
}


void test_flat_floor_tie_break_for_detection()
{
    const lsdtt::lsdtt_raster_meta meta = make_meta(200, 200, 10.0);

    lsdtt::lsdtt_crater_spec cs;
    cs.centre_x       = pixel_centre_x(meta, 100);
    cs.centre_y       = pixel_centre_y(meta, 100);
    cs.radius_m       = 400.0;
    cs.depth_m        = 100.0;
    cs.floor_fraction = 0.4;      // a 160 m radius floor: 16 cells across

    const lsdtt::lsdtt_crater_field_surface surface =
        lsdtt::generate_crater_field_surface(meta, { cs }, 0.02, 1000.0);

    const std::vector<lsdtt::sink_candidate> a =
        lsdtt::find_internal_drainage(surface.dem, lsdtt::sink_detection_options{});
    const std::vector<lsdtt::sink_candidate> b =
        lsdtt::find_internal_drainage(surface.dem, lsdtt::sink_detection_options{});

    check("flat_floor_tie_break: a flat-floored crater is one sink",
          a.size() == 1, std::to_string(a.size()) + " sinks");

    if (a.size() != 1 || b.size() != 1) return;

    // The floor: every cell at the exact minimum elevation.  The generator
    // fades the regional tilt out towards the floor precisely so that this is
    // an exact tie set.
    const auto& z = surface.dem.get_data_array();

    const std::int64_t n_cells =
        static_cast<std::int64_t>(meta.n_rows) *
        static_cast<std::int64_t>(meta.n_cols);

    float z_min = z.data()[0];
    for (std::int64_t i = 1; i < n_cells; ++i)
    {
        if (z.data()[i] < z_min) z_min = z.data()[i];
    }

    std::vector<std::int64_t> floor_cells;
    for (std::int64_t i = 0; i < n_cells; ++i)
    {
        if (z.data()[i] == z_min) floor_cells.push_back(i);
    }

    check("flat_floor_tie_break: the floor really is a tie set",
          floor_cells.size() > 1,
          std::to_string(floor_cells.size()) + " cells at the minimum");

    // Deliberately NOT compared with the generator's pit_rc: a floor is an
    // exact tie, the generator reports the centre cell, and a lowest-flat-index
    // tie-break picks the north-westmost floor cell.  Both are correct.
    const bool member =
        std::find(floor_cells.begin(), floor_cells.end(), a[0].pit_flat) !=
        floor_cells.end();

    check("flat_floor_tie_break: the pit is a member of the floor set",
          member,
          "pit " + std::to_string(a[0].pit_flat) +
          ", generator pit " + std::to_string(surface.pit_flats[0]));

    check("flat_floor_tie_break: the pit is repeatable",
          a[0].pit_flat == b[0].pit_flat && records_equal(a[0], b[0]));

    // The documented rule is lowest elevation, then LOWEST FLAT INDEX.
    check("flat_floor_tie_break: the tie is broken on the lowest flat index",
          a[0].pit_flat == floor_cells.front(),
          "pit " + std::to_string(a[0].pit_flat) +
          ", lowest floor cell " + std::to_string(floor_cells.front()));
}


void test_size_filters()
{
    const crater_field_fixture fx = make_simple_crater_field();

    const std::vector<lsdtt::sink_candidate> all =
        lsdtt::find_internal_drainage(fx.gt.dem, lsdtt::sink_detection_options{});

    if (all.size() != 3)
    {
        check("size_filters: fixture gives three sinks", false,
              std::to_string(all.size()) + " sinks");
        return;
    }

    // The thresholds are derived from the unfiltered run rather than written
    // down as magic numbers, so the test asserts the FILTER SEMANTICS -- drop
    // everything strictly below the threshold, keep the rest -- rather than
    // re-stating the fixture's dimensions.
    const auto smallest_two_dropped =
        [&](auto value_of, const std::string& what, auto set_threshold)
    {
        std::vector<double> values;
        for (const lsdtt::sink_candidate& s : all) values.push_back(value_of(s));

        std::vector<double> sorted = values;
        std::sort(sorted.begin(), sorted.end());

        // A threshold strictly between the smallest and the middle value drops
        // exactly one sink; between the middle and the largest, exactly two.
        const double drop_one = 0.5 * (sorted[0] + sorted[1]);
        const double drop_two = 0.5 * (sorted[1] + sorted[2]);

        for (int n_dropped = 1; n_dropped <= 2; ++n_dropped)
        {
            const double threshold = (n_dropped == 1) ? drop_one : drop_two;

            lsdtt::sink_detection_options opts;
            set_threshold(opts, threshold);

            const std::vector<lsdtt::sink_candidate> kept =
                lsdtt::find_internal_drainage(fx.gt.dem, opts);

            std::vector<std::int64_t> expected;
            for (const lsdtt::sink_candidate& s : all)
            {
                if (value_of(s) >= threshold) expected.push_back(s.pit_flat);
            }
            std::sort(expected.begin(), expected.end());

            check("size_filters: " + what + " drops the smallest " +
                  std::to_string(n_dropped),
                  pit_flats_of(kept) == expected &&
                  kept.size() == static_cast<std::size_t>(3 - n_dropped),
                  "threshold " + std::to_string(threshold) + ", kept " +
                  std::to_string(kept.size()));

            // Survivors are renumbered contiguously from zero.
            bool renumbered = true;
            for (std::size_t i = 0; i < kept.size(); ++i)
            {
                if (kept[i].sink_id != static_cast<int>(i)) renumbered = false;
            }

            check("size_filters: " + what + " renumbers the survivors from zero",
                  renumbered);
        }
    };

    smallest_two_dropped(
        [](const lsdtt::sink_candidate& s) { return static_cast<double>(s.n_pixels); },
        "min_pixels_in_patch",
        [](lsdtt::sink_detection_options& o, double t)
        { o.min_pixels_in_patch = static_cast<int>(std::ceil(t)); });

    smallest_two_dropped(
        [](const lsdtt::sink_candidate& s) { return static_cast<double>(s.depth_m); },
        "min_depth_m",
        [](lsdtt::sink_detection_options& o, double t)
        { o.min_depth_m = static_cast<float>(t); });

    smallest_two_dropped(
        [](const lsdtt::sink_candidate& s) { return s.area_m2; },
        "min_area_m2",
        [](lsdtt::sink_detection_options& o, double t) { o.min_area_m2 = t; });

    smallest_two_dropped(
        [](const lsdtt::sink_candidate& s) { return s.volume_m3; },
        "min_volume_m3",
        [](lsdtt::sink_detection_options& o, double t) { o.min_volume_m3 = t; });

    // A filter that rejects everything returns nothing rather than throwing.
    lsdtt::sink_detection_options brutal;
    brutal.min_depth_m = 1.0e6f;

    check("size_filters: a filter that rejects everything returns an empty list",
          lsdtt::find_internal_drainage(fx.gt.dem, brutal).empty());
}


// =============================================================================
// Section A (continued): nesting
// =============================================================================

void test_nested_craters_hierarchy()
{
    const nested_field_fixture fx = make_nested_crater_field();

    lsdtt::sink_detection_options leaves;
    leaves.leaves_only = true;

    lsdtt::sink_detection_options whole_tree;
    whole_tree.leaves_only = false;

    const std::vector<lsdtt::sink_candidate> leaf_sinks =
        lsdtt::find_internal_drainage(fx.dem, leaves);
    const std::vector<lsdtt::sink_candidate> tree_sinks =
        lsdtt::find_internal_drainage(fx.dem, whole_tree);

    check("nested_craters_hierarchy: leaves_only true gives one sink per hollow",
          leaf_sinks.size() == 3, std::to_string(leaf_sinks.size()) + " sinks");

    std::vector<std::int64_t> expected_leaf_pits = fx.cone_pits;
    std::sort(expected_leaf_pits.begin(), expected_leaf_pits.end());

    check("nested_craters_hierarchy: the leaves are the three hollows",
          pit_flats_of(leaf_sinks) == expected_leaf_pits);

    // THE BINARY-TREE COUNT.  build_depression_hierarchy merges depressions
    // PAIRWISE, so joining three leaves takes two merges and creates TWO
    // meta-nodes, not one enclosing lake: leaves + (leaves - 1) == 5.  Four
    // would be the intuitive answer and is wrong.
    check("nested_craters_hierarchy: leaves_only false gives leaves + (leaves - 1)",
          tree_sinks.size() == 5, std::to_string(tree_sinks.size()) + " sinks");

    if (tree_sinks.size() != 5) return;

    // Ids run in ascending depression label, so children precede their parents:
    // 0, 1, 2 are the leaves and 3, 4 the metas.
    bool leaves_first = true;
    for (int i = 0; i < 3; ++i)
    {
        if (tree_sinks[static_cast<std::size_t>(i)].hierarchy_level != 0)
            leaves_first = false;
    }

    check("nested_craters_hierarchy: the first three sinks are the leaves",
          leaves_first && tree_sinks[3].hierarchy_level > 0 &&
          tree_sinks[4].hierarchy_level > 0,
          "levels " +
          std::to_string(tree_sinks[0].hierarchy_level) + " " +
          std::to_string(tree_sinks[1].hierarchy_level) + " " +
          std::to_string(tree_sinks[2].hierarchy_level) + " " +
          std::to_string(tree_sinks[3].hierarchy_level) + " " +
          std::to_string(tree_sinks[4].hierarchy_level));

    // The measured chain of a three-leaf binary group: two leaves merge into
    // meta 3, meta 3 and the third leaf merge into meta 4, and meta 4 is the
    // top of the tree.
    std::ostringstream chain;
    for (const lsdtt::sink_candidate& s : tree_sinks)
    {
        chain << s.sink_id << "->" << s.parent_sink_id << " ";
    }

    // The SHAPE of the chain is the contract, not which particular hollow
    // merges first: two of the three leaves join into the first meta, the
    // remaining leaf joins that meta in the second, and the second meta is the
    // top of the tree.  (The plan records 0->3, 1->3, 2->4, 3->4, 4->-1 from
    // its own fixture; which two leaves pair up first depends on the geometry,
    // so the measured chain is reported alongside.)
    int n_children_of_3 = 0;
    int n_children_of_4 = 0;

    for (int i = 0; i < 3; ++i)
    {
        const int parent = tree_sinks[static_cast<std::size_t>(i)].parent_sink_id;
        if (parent == 3) ++n_children_of_3;
        if (parent == 4) ++n_children_of_4;
    }

    const bool expected_chain =
        n_children_of_3 == 2 &&
        n_children_of_4 == 1 &&
        tree_sinks[3].parent_sink_id == 4 &&
        tree_sinks[4].parent_sink_id == -1;

    check("nested_craters_hierarchy: two leaves merge into meta 3, the third and "
          "meta 3 into meta 4, and meta 4 is the top",
          expected_chain, chain.str());

    // A parent must be at least as deep as its children and must contain them.
    bool nesting_sane = true;
    for (const lsdtt::sink_candidate& s : tree_sinks)
    {
        if (s.parent_sink_id < 0) continue;

        const lsdtt::sink_candidate& parent =
            tree_sinks[static_cast<std::size_t>(s.parent_sink_id)];

        if (parent.hierarchy_level <= s.hierarchy_level) nesting_sane = false;
        if (parent.n_pixels        <  s.n_pixels)        nesting_sane = false;
        if (parent.spill_elevation <  s.spill_elevation) nesting_sane = false;
    }

    check("nested_craters_hierarchy: every parent is deeper, larger and spills higher",
          nesting_sane);

    // With leaves_only the reported set is exactly the leaves of the same tree.
    std::vector<std::int64_t> tree_leaf_pits;
    for (const lsdtt::sink_candidate& s : tree_sinks)
    {
        if (s.hierarchy_level == 0) tree_leaf_pits.push_back(s.pit_flat);
    }
    std::sort(tree_leaf_pits.begin(), tree_leaf_pits.end());

    check("nested_craters_hierarchy: leaves_only is the leaf subset of the tree",
          tree_leaf_pits == pit_flats_of(leaf_sinks));
}


void test_shared_pit_cell_is_not_a_key()
{
    const nested_field_fixture fx = make_nested_crater_field();

    lsdtt::sink_detection_options opts;
    opts.leaves_only = false;

    const std::vector<lsdtt::sink_candidate> sinks =
        lsdtt::find_internal_drainage(fx.dem, opts);

    // A meta-depression's pit IS its lowest child's pit, so the lake and the
    // hollow at the bottom of it share one cell.
    std::map<std::int64_t, std::vector<std::size_t>> by_pit;
    for (std::size_t i = 0; i < sinks.size(); ++i)
    {
        by_pit[sinks[i].pit_flat].push_back(i);
    }

    check("shared_pit_cell_is_not_a_key: fewer distinct pit cells than sinks",
          by_pit.size() < sinks.size(),
          std::to_string(by_pit.size()) + " distinct pit cells for " +
          std::to_string(sinks.size()) + " sinks");

    std::size_t n_shared = 0;
    bool        distinguishable = true;

    for (const auto& entry : by_pit)
    {
        if (entry.second.size() < 2) continue;

        ++n_shared;

        std::set<int> levels;
        std::set<int> ids;

        for (const std::size_t i : entry.second)
        {
            levels.insert(sinks[i].hierarchy_level);
            ids.insert(sinks[i].sink_id);
        }

        // Nothing may collapse two sinks into one: they keep distinct ids AND
        // are separable by hierarchy_level, which is the documented
        // disambiguator.
        if (levels.size() != entry.second.size()) distinguishable = false;
        if (ids.size()    != entry.second.size()) distinguishable = false;
    }

    check("shared_pit_cell_is_not_a_key: at least one pit cell is shared",
          n_shared >= 1, std::to_string(n_shared) + " shared pit cells");

    check("shared_pit_cell_is_not_a_key: sinks sharing a cell stay distinct, "
          "separated by hierarchy_level",
          distinguishable);

    // The deepest child and its parent chain really do share the cell: the
    // lowest-elevation leaf's pit must also be the pit of the top meta.
    const lsdtt::sink_candidate* lowest_leaf = nullptr;
    const lsdtt::sink_candidate* top_meta    = nullptr;

    for (const lsdtt::sink_candidate& s : sinks)
    {
        if (s.hierarchy_level == 0 &&
            (lowest_leaf == nullptr || s.pit_elevation < lowest_leaf->pit_elevation))
        {
            lowest_leaf = &s;
        }
        if (s.parent_sink_id == -1 &&
            (top_meta == nullptr || s.hierarchy_level > top_meta->hierarchy_level))
        {
            top_meta = &s;
        }
    }

    check("shared_pit_cell_is_not_a_key: the top meta's pit is the lowest leaf's pit",
          lowest_leaf != nullptr && top_meta != nullptr &&
          lowest_leaf->pit_flat == top_meta->pit_flat &&
          lowest_leaf->sink_id  != top_meta->sink_id);
}


void test_min_pixels_filters_a_leaf()
{
    const nested_field_fixture fx = make_nested_crater_field();
    const lsdtt::lsdtt_raster_meta& meta = fx.dem.get_meta();

    lsdtt::sink_detection_options all_opts;
    all_opts.leaves_only = false;

    const std::vector<lsdtt::sink_candidate> all =
        lsdtt::find_internal_drainage(fx.dem, all_opts);

    // The smallest leaf, and a threshold that removes it and nothing else.
    const lsdtt::sink_candidate* smallest = nullptr;
    int second_smallest = 0;

    std::vector<int> leaf_sizes;
    for (const lsdtt::sink_candidate& s : all)
    {
        if (s.hierarchy_level != 0) continue;
        leaf_sizes.push_back(s.n_pixels);
        if (smallest == nullptr || s.n_pixels < smallest->n_pixels) smallest = &s;
    }

    std::sort(leaf_sizes.begin(), leaf_sizes.end());
    if (leaf_sizes.size() >= 2) second_smallest = leaf_sizes[1];

    if (smallest == nullptr || second_smallest <= smallest->n_pixels)
    {
        check("min_pixels_filters_a_leaf: fixture has distinct leaf sizes", false);
        return;
    }

    const std::int64_t filtered_pit = smallest->pit_flat;

    lsdtt::sink_detection_options filtered = all_opts;
    filtered.min_pixels_in_patch = smallest->n_pixels + 1;

    const std::vector<lsdtt::sink_candidate> kept =
        lsdtt::find_internal_drainage(fx.dem, filtered);

    // Exactly one node leaves the reported tree.  A META whose pit happens to
    // be the filtered leaf's pit may legitimately remain -- a pit cell is not a
    // unique key -- so the leaf itself is checked separately below.
    check("min_pixels_filters_a_leaf: exactly one sink is dropped",
          kept.size() == all.size() - 1,
          std::to_string(all.size()) + " sinks before, " +
          std::to_string(kept.size()) + " after");

    bool leaf_gone = true;
    for (const lsdtt::sink_candidate& s : kept)
    {
        if (s.hierarchy_level == 0 && s.pit_flat == filtered_pit) leaf_gone = false;
    }

    check("min_pixels_filters_a_leaf: no leaf remains at the filtered pit",
          leaf_gone);

    // THE POINT OF THE TEST: the filtered leaf's cells do not vanish, they
    // reattach to the nearest surviving ancestor.  So the pit cell of the
    // dropped hollow now lies inside some surviving sink's footprint.
    const lsdtt::depression_polygons_result polygons =
        lsdtt::sink_depressions_to_polygons(fx.dem, kept, filtered);

    const int pit_row = static_cast<int>(filtered_pit / meta.n_cols);
    const int pit_col = static_cast<int>(filtered_pit % meta.n_cols);
    const double px   = pixel_centre_x(meta, pit_col);
    const double py   = pixel_centre_y(meta, pit_row);

    int owner = -1;
    for (const lsdtt::lsdtt_vector_record& record : polygons.records)
    {
        if (!std::holds_alternative<lsdtt::lsdtt_polygon_geom>(record.geometry))
            continue;

        const lsdtt::lsdtt_polygon_geom& polygon =
            std::get<lsdtt::lsdtt_polygon_geom>(record.geometry);

        if (polygon_contains(polygon, px, py))
        {
            owner = static_cast<int>(
                std::get<std::int64_t>(record.field_values[0]));
        }
    }

    check("min_pixels_filters_a_leaf: the dropped leaf's cells reattach to a "
          "surviving ancestor",
          owner >= 0,
          owner >= 0 ? ("claimed by sink " + std::to_string(owner))
                     : "no surviving footprint covers the dropped pit");

    // ... and that owner really is an ancestor, i.e. a meta-depression rather
    // than one of the sibling hollows.
    const lsdtt::sink_candidate* owning_sink = nullptr;
    for (const lsdtt::sink_candidate& s : kept)
    {
        if (s.sink_id == owner) owning_sink = &s;
    }

    check("min_pixels_filters_a_leaf: the new owner is an ancestor, not a sibling",
          owning_sink != nullptr && owning_sink->hierarchy_level > 0,
          owning_sink != nullptr
              ? ("level " + std::to_string(owning_sink->hierarchy_level))
              : "owner not in the kept list");
}


// =============================================================================
// min_pixels_in_patch must CONSERVE TERRAIN
//
// The size filter decides which SINKS are reported. It must never make DEM cells
// disappear from every footprint: a filtered leaf is merged into its nearest
// surviving ancestor, so the ground it covered is still inside some reported
// footprint, just attributed to a coarser sink.
//
// This is the invariant that phase 5b of
// docs/plan_docs/depression_hierarchy_fix_plan.md restored. It was violated
// because a sink whose cells are disconnected -- exactly what happens when a
// dropped leaf's patch reattaches to an ancestor across a sibling's footprint --
// had all but one of its parts silently discarded by the polygon tracer.
//
// It is pinned here on total covered AREA rather than on any one pit, so it does
// not depend on the crater geometry that happened to expose the bug. Footprints
// of different sinks do not overlap (asserted separately), so summing
// exterior-minus-holes over every record gives the union.
// =============================================================================
double total_traced_area(const lsdtt::depression_polygons_result& polygons)
{
    double total = 0.0;
    for (const lsdtt::lsdtt_vector_record& record : polygons.records)
    {
        if (!std::holds_alternative<lsdtt::lsdtt_polygon_geom>(record.geometry))
            continue;
        total += polygon_area(
            std::get<lsdtt::lsdtt_polygon_geom>(record.geometry));
    }
    return total;
}

void test_min_pixels_conserves_terrain()
{
    const nested_field_fixture fx = make_nested_crater_field();

    lsdtt::sink_detection_options all_opts;
    all_opts.leaves_only = false;

    const std::vector<lsdtt::sink_candidate> all =
        lsdtt::find_internal_drainage(fx.dem, all_opts);

    // The smallest leaf, and a threshold that removes it and nothing else.
    const lsdtt::sink_candidate* smallest = nullptr;
    for (const lsdtt::sink_candidate& s : all)
    {
        if (s.hierarchy_level != 0) continue;
        if (smallest == nullptr || s.n_pixels < smallest->n_pixels) smallest = &s;
    }

    check("min_pixels_conserves_terrain: fixture has a leaf to filter",
          smallest != nullptr);
    if (smallest == nullptr) return;

    lsdtt::sink_detection_options filtered = all_opts;
    filtered.min_pixels_in_patch = smallest->n_pixels + 1;

    const std::vector<lsdtt::sink_candidate> kept =
        lsdtt::find_internal_drainage(fx.dem, filtered);

    const double area_all =
        total_traced_area(sink_depressions_to_polygons(fx.dem, all, all_opts));
    const double area_kept =
        total_traced_area(sink_depressions_to_polygons(fx.dem, kept, filtered));

    // Areas are traced along cell boundaries, so equality is exact up to the
    // accumulation of many shoelace terms; a relative bound well below one cell
    // is the right tolerance.
    const double tol = 1.0e-9 * std::max(area_all, 1.0);

    check("min_pixels_conserves_terrain: filtering a leaf covers the same ground",
          std::fabs(area_all - area_kept) <= tol,
          "unfiltered " + std::to_string(area_all) +
              " m2 vs filtered " + std::to_string(area_kept) + " m2");
}

// =============================================================================
// Section B: the legacy fill-difference path
// =============================================================================

void test_fill_difference_agrees_on_a_simple_field()
{
    const crater_field_fixture fx = make_simple_crater_field();

    lsdtt::sink_detection_options hierarchy;

    lsdtt::sink_detection_options legacy;
    legacy.method = lsdtt::sink_detection_method::fill_difference;

    const std::vector<lsdtt::sink_candidate> from_hierarchy =
        lsdtt::find_internal_drainage(fx.gt.dem, hierarchy);
    const std::vector<lsdtt::sink_candidate> from_fill =
        lsdtt::find_internal_drainage(fx.gt.dem, legacy);

    check("fill_difference_agrees: same number of sinks on separate craters",
          from_fill.size() == from_hierarchy.size(),
          std::to_string(from_fill.size()) + " vs " +
          std::to_string(from_hierarchy.size()));

    // The two methods answer DIFFERENT questions -- the lowest cell, and the
    // cell the fill had to raise furthest -- so agreement is not automatic.  It
    // is structural HERE because build_fill_difference_state fills with
    // min_slope = 0: the fill is then flat across a component, so
    // (spill - z) is maximised at exactly the lowest cell and the two
    // definitions coincide.  A POSITIVE min_slope lays an epsilon gradient
    // across the flat and IS EXPECTED to move the deepest-fill cell by a pixel
    // or two (plan section 5.4).  If this test ever fails after a min_slope
    // change, that is the reason, and the Greenland cross-validation in plan
    // 6.2 is where the legacy value must be checked.
    check("fill_difference_agrees: the same pit cells, exactly",
          pit_flats_of(from_fill) == pit_flats_of(from_hierarchy));

    // Every fill-difference sink is a flat, parentless component: the method
    // has no tree.
    bool flat_tree = true;
    for (const lsdtt::sink_candidate& s : from_fill)
    {
        if (s.hierarchy_level != 0 || s.parent_sink_id != -1) flat_tree = false;
    }

    check("fill_difference_agrees: fill-difference sinks carry no tree",
          flat_tree);

    // Footprint size and volume measure the same depression, but not by the
    // same definition: a fill-difference component holds only the cells the
    // fill RAISED, whereas a hierarchy depression holds every cell it contains.
    // The cells sitting exactly at the spill elevation, and the pit itself, are
    // raised by zero and so are missing from the legacy footprint.  That is a
    // rim-length effect, of order one cell around the boundary, so the two
    // agree to a few per cent and not to a cell.
    bool        sizes_agree = true;
    std::string worst;
    std::ostringstream measured;

    for (const lsdtt::sink_candidate& s : from_fill)
    {
        const lsdtt::sink_candidate* h = sink_at_pit(from_hierarchy, s.pit_flat);
        if (h == nullptr) { sizes_agree = false; continue; }

        measured << s.n_pixels << " vs " << h->n_pixels << " px; ";

        // The legacy footprint is a SUBSET of the hierarchy footprint, never
        // larger: that direction is structural and is asserted, the magnitude
        // is bounded loosely.
        if (s.n_pixels > h->n_pixels)
        {
            sizes_agree = false;
            worst = "legacy footprint larger: " + std::to_string(s.n_pixels) +
                    " vs " + std::to_string(h->n_pixels);
        }

        if (h->n_pixels - s.n_pixels > h->n_pixels / 20)
        {
            sizes_agree = false;
            worst = "n_pixels " + std::to_string(s.n_pixels) + " vs " +
                    std::to_string(h->n_pixels);
        }

        if (std::fabs(s.volume_m3 - h->volume_m3) > 0.02 * h->volume_m3)
        {
            sizes_agree = false;
            worst = "volume " + std::to_string(s.volume_m3) + " vs " +
                    std::to_string(h->volume_m3);
        }
    }

    check("fill_difference_agrees: the legacy footprint is the hierarchy's, less "
          "its unraised rim",
          sizes_agree, worst.empty() ? measured.str() : worst);
}


void test_fill_difference_diverges_on_a_nested_field()
{
    const nested_field_fixture fx = make_nested_crater_field();

    lsdtt::sink_detection_options hierarchy;   // leaves only

    lsdtt::sink_detection_options legacy;
    legacy.method = lsdtt::sink_detection_method::fill_difference;

    const std::vector<lsdtt::sink_candidate> leaves =
        lsdtt::find_internal_drainage(fx.dem, hierarchy);
    const std::vector<lsdtt::sink_candidate> components =
        lsdtt::find_internal_drainage(fx.dem, legacy);

    // A DOCUMENTED DIVERGENCE, not a defect (plan section 5.4).  Connected
    // components of the fill difference merge a lake and every hollow inside it
    // into ONE patch, so the legacy method reports one sink where the hierarchy
    // reports three.
    check("fill_difference_diverges: one component swallows the nested lake",
          components.size() == 1,
          std::to_string(components.size()) + " components vs " +
          std::to_string(leaves.size()) + " hierarchy leaves");

    check("fill_difference_diverges: the hierarchy finds strictly more sinks",
          leaves.size() > components.size());

    if (components.size() != 1) return;

    // Its single pit sits exactly on one of the hierarchy's leaves: the pit
    // POSITIONS agree, it is the COUNT that differs.
    const bool on_a_leaf =
        sink_at_pit(leaves, components[0].pit_flat) != nullptr;

    check("fill_difference_diverges: its pit coincides with a hierarchy leaf",
          on_a_leaf, "pit " + std::to_string(components[0].pit_flat));

    // The swallowing component is larger than any single hollow.
    int largest_leaf = 0;
    for (const lsdtt::sink_candidate& s : leaves)
    {
        largest_leaf = std::max(largest_leaf, s.n_pixels);
    }

    check("fill_difference_diverges: the component is larger than any hollow",
          components[0].n_pixels > largest_leaf,
          std::to_string(components[0].n_pixels) + " cells vs largest hollow " +
          std::to_string(largest_leaf));

    report("fill_difference_diverges: measured component size",
           std::to_string(components[0].n_pixels) + " px, depth " +
           std::to_string(components[0].depth_m) + " m");
}


void test_fill_difference_min_pixels_is_legacy()
{
    const crater_field_fixture fx = make_simple_crater_field();

    lsdtt::sink_detection_options legacy;
    legacy.method = lsdtt::sink_detection_method::fill_difference;

    const std::vector<lsdtt::sink_candidate> all =
        lsdtt::find_internal_drainage(fx.gt.dem, legacy);

    if (all.size() != 3)
    {
        check("fill_difference_min_pixels: fixture gives three components", false,
              std::to_string(all.size()) + " components");
        return;
    }

    std::vector<int> sizes;
    for (const lsdtt::sink_candidate& s : all) sizes.push_back(s.n_pixels);
    std::sort(sizes.begin(), sizes.end());

    // The legacy meaning: minimum_pixels_in_patch is a filter on the SIZE OF
    // THE CONNECTED COMPONENT, applied when the components are labelled, so a
    // component smaller than the threshold is not merely unreported -- it is
    // not a patch at all.
    lsdtt::sink_detection_options filtered = legacy;
    filtered.min_pixels_in_patch = sizes[0] + 1;

    const std::vector<lsdtt::sink_candidate> kept =
        lsdtt::find_internal_drainage(fx.gt.dem, filtered);

    std::vector<std::int64_t> expected;
    for (const lsdtt::sink_candidate& s : all)
    {
        if (s.n_pixels >= sizes[0] + 1) expected.push_back(s.pit_flat);
    }
    std::sort(expected.begin(), expected.end());

    check("fill_difference_min_pixels: the smallest patch is dropped",
          kept.size() == 2 && pit_flats_of(kept) == expected,
          std::to_string(kept.size()) + " kept, threshold " +
          std::to_string(sizes[0] + 1));

    filtered.min_pixels_in_patch = sizes[2] + 1;

    check("fill_difference_min_pixels: a threshold above every patch keeps none",
          lsdtt::find_internal_drainage(fx.gt.dem, filtered).empty());

    // Zero is not a legal patch size, and the legacy code silently accepted it.
    lsdtt::sink_detection_options zero = legacy;
    zero.min_pixels_in_patch = 0;

    check_throws("fill_difference_min_pixels: a zero threshold is refused",
                 [&] { (void)lsdtt::find_internal_drainage(fx.gt.dem, zero); },
                 "min_pixels_in_patch must be >= 1");
}


// =============================================================================
// Section D: depression footprints as polygons
// =============================================================================

void test_depression_polygon_join()
{
    const crater_field_fixture fx = make_simple_crater_field();

    lsdtt::sink_detection_options opts;   // leaves only

    const std::vector<lsdtt::sink_candidate> sinks =
        lsdtt::find_internal_drainage(fx.gt.dem, opts);

    const lsdtt::depression_polygons_result polygons =
        lsdtt::sink_depressions_to_polygons(fx.gt.dem, sinks, opts);

    check("depression_polygon_join: one polygon per sink in the leaves-only case",
          polygons.records.size() == sinks.size(),
          std::to_string(polygons.records.size()) + " polygons for " +
          std::to_string(sinks.size()) + " sinks");

    check("depression_polygon_join: the schema leads with sink_id",
          !polygons.schema.field_names.empty() &&
          polygons.schema.field_names[0] == "sink_id" &&
          polygons.schema.field_types[0] == lsdtt::lsdtt_vector_field_type::int64_type);

    std::set<std::int64_t> polygon_ids;
    bool areas_match  = true;
    bool attributes_match = true;
    std::string detail;

    for (const lsdtt::lsdtt_vector_record& record : polygons.records)
    {
        const std::int64_t sink_id =
            std::get<std::int64_t>(record.field_values[0]);
        polygon_ids.insert(sink_id);

        const lsdtt::sink_candidate* s = nullptr;
        for (const lsdtt::sink_candidate& candidate : sinks)
        {
            if (candidate.sink_id == static_cast<int>(sink_id)) s = &candidate;
        }

        if (s == nullptr) { attributes_match = false; continue; }

        // The polygon's attributes are copied from the sink, so the point and
        // polygon layers can never disagree.
        if (std::get<double>(record.field_values[4]) != s->area_m2 ||
            std::get<double>(record.field_values[5]) != s->volume_m3 ||
            std::get<std::int64_t>(record.field_values[6]) !=
                static_cast<std::int64_t>(s->n_pixels))
        {
            attributes_match = false;
        }

        if (!std::holds_alternative<lsdtt::lsdtt_polygon_geom>(record.geometry))
        {
            areas_match = false;
            continue;
        }

        const double traced =
            polygon_area(std::get<lsdtt::lsdtt_polygon_geom>(record.geometry));

        // The footprint is traced along cell boundaries, so its area is the
        // pixel count times the cell area -- the same number area_m2 carries.
        if (std::fabs(traced - s->area_m2) > 1.0e-6 * s->area_m2)
        {
            areas_match = false;
            detail = "traced " + std::to_string(traced) +
                     " vs area_m2 " + std::to_string(s->area_m2);
        }
    }

    check("depression_polygon_join: sink_id joins 1:1 to the point layer",
          polygon_ids.size() == sinks.size());

    check("depression_polygon_join: polygon attributes come from the sink",
          attributes_match);

    check("depression_polygon_join: traced area matches area_m2",
          areas_match, detail);

    // A subset in gives exactly that subset's polygons back, which is what the
    // route stage needs for the accepted-only export.
    std::vector<lsdtt::sink_candidate> subset = { sinks[1] };

    const lsdtt::depression_polygons_result one =
        lsdtt::sink_depressions_to_polygons(fx.gt.dem, subset, opts);

    check("depression_polygon_join: an accepted-only subset gets its own polygons",
          one.records.size() == 1 &&
          std::get<std::int64_t>(one.records[0].field_values[0]) ==
              static_cast<std::int64_t>(sinks[1].sink_id),
          std::to_string(one.records.size()) + " polygons");

    std::vector<lsdtt::sink_candidate> sparse_ids = sinks;
    sparse_ids[0].sink_id = std::numeric_limits<int>::max();
    sparse_ids[1].sink_id = 1000000000;
    sparse_ids[2].sink_id = 42;
    const lsdtt::depression_polygons_result sparse_polygons =
        lsdtt::sink_depressions_to_polygons(fx.gt.dem, sparse_ids, opts);
    std::set<std::int64_t> sparse_polygon_ids;
    for (const lsdtt::lsdtt_vector_record& record : sparse_polygons.records)
    {
        sparse_polygon_ids.insert(
            std::get<std::int64_t>(record.field_values[0]));
    }
    check("depression_polygon_join: sparse durable ids are attributes, not dense indices",
          sparse_polygons.records.size() == sparse_ids.size() &&
          sparse_polygon_ids.count(std::numeric_limits<int>::max()) == 1 &&
          sparse_polygon_ids.count(1000000000) == 1 &&
          sparse_polygon_ids.count(42) == 1);

    check("depression_polygon_join: an empty sink list gives no polygons",
          lsdtt::sink_depressions_to_polygons(
              fx.gt.dem, std::vector<lsdtt::sink_candidate>{}, opts)
              .records.empty());
}


void test_nested_sink_may_have_no_polygon()
{
    const nested_field_fixture fx = make_nested_crater_field();

    lsdtt::sink_detection_options opts;
    opts.leaves_only = false;

    const std::vector<lsdtt::sink_candidate> sinks =
        lsdtt::find_internal_drainage(fx.dem, opts);

    const lsdtt::depression_polygons_result polygons =
        lsdtt::sink_depressions_to_polygons(fx.dem, sinks, opts);

    // Each cell is drawn ONCE, into the innermost depression that owns it, so
    // an intermediate node entirely covered by its children has an EMPTY
    // footprint and therefore no polygon at all.  "One polygon per sink" is
    // safe only in the leaves-only case.
    check("nested_sink_may_have_no_polygon: 5 sinks give 4 polygons",
          sinks.size() == 5 && polygons.records.size() == 4,
          std::to_string(sinks.size()) + " sinks, " +
          std::to_string(polygons.records.size()) + " polygons");

    std::set<int> with_polygon;
    for (const lsdtt::lsdtt_vector_record& record : polygons.records)
    {
        with_polygon.insert(
            static_cast<int>(std::get<std::int64_t>(record.field_values[0])));
    }

    std::vector<const lsdtt::sink_candidate*> missing;
    for (const lsdtt::sink_candidate& s : sinks)
    {
        if (with_polygon.find(s.sink_id) == with_polygon.end()) missing.push_back(&s);
    }

    check("nested_sink_may_have_no_polygon: exactly one sink has no footprint",
          missing.size() == 1, std::to_string(missing.size()) + " missing");

    check("nested_sink_may_have_no_polygon: the sink with no footprint is a meta",
          missing.size() == 1 && missing[0]->hierarchy_level > 0,
          missing.empty() ? "" : describe(*missing[0]));

    // It is still in the point layer with a non-zero area_m2: the depression
    // exists, it is only its EXCLUSIVE footprint that is empty.
    check("nested_sink_may_have_no_polygon: it still reports an area",
          missing.size() == 1 && missing[0]->area_m2 > 0.0);

    // The polygons partition the cells they cover: with each cell drawn once,
    // no two footprints may overlap, so the traced areas sum to at most the
    // total valid area.
    double traced_total = 0.0;
    for (const lsdtt::lsdtt_vector_record& record : polygons.records)
    {
        if (!std::holds_alternative<lsdtt::lsdtt_polygon_geom>(record.geometry))
            continue;
        traced_total +=
            polygon_area(std::get<lsdtt::lsdtt_polygon_geom>(record.geometry));
    }

    const lsdtt::lsdtt_raster_meta& meta = fx.dem.get_meta();
    const double total_area =
        static_cast<double>(meta.n_rows) * meta.n_cols * meta.pixel_area();

    // A sum of area_m2 over the same sinks double-counts, because a parent's
    // area includes its children.  That contrast is the point of plan 5.3b.
    double attribute_total = 0.0;
    for (const lsdtt::sink_candidate& s : sinks) attribute_total += s.area_m2;

    check("nested_sink_may_have_no_polygon: footprints do not overlap",
          traced_total < total_area && traced_total < attribute_total,
          "traced " + std::to_string(traced_total) +
          ", summed area_m2 " + std::to_string(attribute_total));
}


// =============================================================================
// Section C: the review round-trip
// =============================================================================

struct lon_lat
{
    double lon = 0.0;
    double lat = 0.0;
};

// Raster CRS -> EPSG:4326, so a fixture can be specified in the pixel whose
// row/col the test asserts.
lon_lat raster_xy_to_lon_lat(const lsdtt::lsdtt_raster_meta& meta,
                             double x,
                             double y)
{
    OGRSpatialReference src_srs;
    OGRSpatialReference dst_srs;

    src_srs.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);
    dst_srs.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);

    if (src_srs.importFromWkt(meta.crs_wkt.c_str()) != OGRERR_NONE ||
        dst_srs.importFromEPSG(4326) != OGRERR_NONE)
    {
        throw std::runtime_error(
            "test_internal_drainage_detect: could not build the CRS pair");
    }

    OGRCoordinateTransformation* tx =
        OGRCreateCoordinateTransformation(&src_srs, &dst_srs);

    if (tx == nullptr)
    {
        throw std::runtime_error(
            "test_internal_drainage_detect: could not build the transform");
    }

    double xx = x;
    double yy = y;
    const int ok = tx->Transform(1, &xx, &yy);
    OCTDestroyCoordinateTransformation(tx);

    if (ok == 0)
    {
        throw std::runtime_error(
            "test_internal_drainage_detect: transform to EPSG:4326 failed");
    }

    return lon_lat{xx, yy};
}

std::string epsg_4326_wkt()
{
    OGRSpatialReference srs;
    srs.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);

    if (srs.importFromEPSG(4326) != OGRERR_NONE)
    {
        throw std::runtime_error(
            "test_internal_drainage_detect: could not build EPSG:4326");
    }

    char* wkt = nullptr;
    srs.exportToWkt(&wkt);
    const std::string out = (wkt != nullptr) ? std::string(wkt) : std::string();
    CPLFree(wkt);

    return out;
}

// Writes a small CSV by hand, so the position-precedence fixtures can carry
// deliberately conflicting columns.
void write_csv_by_hand(const std::string&                           path,
                       const std::vector<std::string>&              headers,
                       const std::vector<std::vector<std::string>>& rows)
{
    std::ofstream out(path, std::ios::binary);

    for (std::size_t i = 0; i < headers.size(); ++i)
    {
        if (i != 0) out << ",";
        out << headers[i];
    }
    out << "\n";

    for (const std::vector<std::string>& row : rows)
    {
        for (std::size_t i = 0; i < row.size(); ++i)
        {
            if (i != 0) out << ",";
            out << row[i];
        }
        out << "\n";
    }
}

std::string number_to_string(double value)
{
    std::ostringstream out;
    out.precision(17);
    out << value;
    return out.str();
}


void test_candidate_round_trip_csv(const std::string&                        dir,
                                   const std::vector<lsdtt::sink_candidate>& sinks,
                                   const lsdtt::lsdtt_raster_meta&           meta)
{
    lsdtt::lsdtt_data_loader loader;

    const std::string first  = dir + "/round_trip_a.csv";
    const std::string second = dir + "/round_trip_b.csv";

    lsdtt::write_sink_candidates(loader, first, sinks, meta,
                                 lsdtt::sink_io_format::csv);

    const std::vector<lsdtt::sink_candidate> read_back =
        lsdtt::read_sink_candidates(first, meta, lsdtt::sink_io_format::csv);

    std::string difference;
    check("candidate_round_trip_csv: the records survive the round trip",
          vectors_equal(sinks, read_back, difference), difference);

    lsdtt::write_sink_candidates(loader, second, read_back, meta,
                                 lsdtt::sink_io_format::csv);

    check("candidate_round_trip_csv: write -> read -> write is byte-identical",
          read_whole_file(first) == read_whole_file(second));

    // accepted: 0 rows survive the read carrying the flag, and are dropped by
    // accepted_sinks() -- which is what the route stage consumes.
    std::vector<lsdtt::sink_candidate> edited = sinks;
    edited[1].accepted = false;

    const std::string edited_path = dir + "/round_trip_rejected.csv";
    lsdtt::write_sink_candidates(loader, edited_path, edited, meta,
                                 lsdtt::sink_io_format::csv);

    const std::vector<lsdtt::sink_candidate> edited_back =
        lsdtt::read_sink_candidates(edited_path, meta, lsdtt::sink_io_format::csv);

    check("candidate_round_trip_csv: the accepted flag survives the round trip",
          edited_back.size() == edited.size() &&
          !edited_back[1].accepted &&
          edited_back[0].accepted && edited_back[2].accepted);

    const std::vector<lsdtt::sink_candidate> kept =
        lsdtt::accepted_sinks(edited_back);

    check("candidate_round_trip_csv: accepted rows keep durable ids and get routing ids",
          kept.size() == edited.size() - 1 &&
          kept[0].sink_id == 0 && kept[1].sink_id == 2 &&
          kept[0].routing_sink_id == 0 && kept[1].routing_sink_id == 1 &&
          kept[0].pit_flat == sinks[0].pit_flat &&
          kept[1].pit_flat == sinks[2].pit_flat,
          std::to_string(kept.size()) + " accepted");

    // A file with NO accepted column keeps every row: the reviewer deleted the
    // rows they rejected, which is the legacy onlyvalid.csv semantics.
    const std::string no_column = dir + "/round_trip_no_accepted.csv";
    {
        std::istringstream in(read_whole_file(edited_path));
        std::ofstream      out(no_column, std::ios::binary);
        std::string        line;

        while (std::getline(in, line))
        {
            if (!line.empty() && line.back() == '\r') line.pop_back();
            const std::size_t last = line.find_last_of(',');
            out << ((last == std::string::npos) ? line : line.substr(0, last))
                << "\n";
        }
    }

    const std::vector<lsdtt::sink_candidate> without_column =
        lsdtt::read_sink_candidates(no_column, meta, lsdtt::sink_io_format::csv);

    bool all_accepted = true;
    for (const lsdtt::sink_candidate& s : without_column)
    {
        if (!s.accepted) all_accepted = false;
    }

    check("candidate_round_trip_csv: a file with no accepted column keeps every row",
          without_column.size() == edited.size() && all_accepted,
          std::to_string(without_column.size()) + " rows");
}


void test_candidate_round_trip_fgb(const std::string&                        dir,
                                   const std::vector<lsdtt::sink_candidate>& sinks,
                                   const lsdtt::lsdtt_raster_meta&           meta)
{
    lsdtt::lsdtt_data_loader loader;

    const std::string first = dir + "/round_trip_a.fgb";

    lsdtt::write_sink_candidates(loader, first, sinks, meta,
                                 lsdtt::sink_io_format::flatgeobuf);

    const std::vector<lsdtt::sink_candidate> read_back =
        lsdtt::read_sink_candidates(first, meta, lsdtt::sink_io_format::flatgeobuf);

    std::string difference;
    check("candidate_round_trip_fgb: the records survive the round trip",
          vectors_equal(sinks, read_back, difference), difference);

    // Geometry survives to full double precision: the coordinate written is the
    // coordinate read, bit for bit, not merely close.
    bool exact_geometry = true;
    for (std::size_t i = 0; i < read_back.size() && i < sinks.size(); ++i)
    {
        if (read_back[i].pit_x != sinks[i].pit_x ||
            read_back[i].pit_y != sinks[i].pit_y) exact_geometry = false;
    }

    check("candidate_round_trip_fgb: geometry survives to full double precision",
          exact_geometry);

    std::vector<lsdtt::sink_candidate> edited = sinks;
    edited[1].accepted = false;

    const std::string edited_path = dir + "/round_trip_rejected.fgb";
    lsdtt::write_sink_candidates(loader, edited_path, edited, meta,
                                 lsdtt::sink_io_format::flatgeobuf);

    const std::vector<lsdtt::sink_candidate> edited_back =
        lsdtt::read_sink_candidates(edited_path, meta,
                                    lsdtt::sink_io_format::flatgeobuf);

    check("candidate_round_trip_fgb: the accepted flag survives the round trip",
          edited_back.size() == edited.size() && !edited_back[1].accepted);

    check("candidate_round_trip_fgb: accepted: 0 features are dropped",
          lsdtt::accepted_sinks(edited_back).size() == edited.size() - 1);

    // A layer with no accepted field at all: every feature present is accepted.
    lsdtt::lsdtt_vector_schema schema;
    schema.field_names = { "sink_id", "easting", "northing" };
    schema.field_types = { lsdtt::lsdtt_vector_field_type::int64_type,
                           lsdtt::lsdtt_vector_field_type::double_type,
                           lsdtt::lsdtt_vector_field_type::double_type };

    const std::string no_field = dir + "/round_trip_no_accepted.fgb";

    loader.begin_vector_write(no_field, "sink_candidates", schema, meta.crs_wkt,
                              lsdtt::vector_output_format::flatgeobuf, false);

    for (const lsdtt::sink_candidate& s : sinks)
    {
        loader.write_point(s.pit_x, s.pit_y,
                           { static_cast<std::int64_t>(s.sink_id),
                             s.pit_x, s.pit_y });
    }

    loader.finish_vector_write();

    const std::vector<lsdtt::sink_candidate> without_field =
        lsdtt::read_sink_candidates(no_field, meta,
                                    lsdtt::sink_io_format::flatgeobuf);

    bool all_accepted = true;
    for (const lsdtt::sink_candidate& s : without_field)
    {
        if (!s.accepted) all_accepted = false;
    }

    check("candidate_round_trip_fgb: a layer with no accepted field keeps every "
          "feature",
          without_field.size() == sinks.size() && all_accepted);
}


void test_csv_fgb_record_equality(const std::string&                        dir,
                                  const std::vector<lsdtt::sink_candidate>& sinks,
                                  const lsdtt::lsdtt_raster_meta&           meta)
{
    lsdtt::lsdtt_data_loader loader;

    const std::string csv_path = dir + "/equality.csv";
    const std::string fgb_path = dir + "/equality.fgb";

    lsdtt::write_sink_candidates(loader, csv_path, sinks, meta,
                                 lsdtt::sink_io_format::csv);
    lsdtt::write_sink_candidates(loader, fgb_path, sinks, meta,
                                 lsdtt::sink_io_format::flatgeobuf);

    // Compared through the sink_candidate accessors, NEVER through the
    // attribute kind: a value written as 0 comes back kind::integer from CSV
    // (typed by inspection) and kind::real from OGR (typed by the schema), and
    // plan 6.4(4) records that consumers must not compare type across formats.
    const std::vector<lsdtt::sink_candidate> from_csv =
        lsdtt::read_sink_candidates(csv_path, meta, lsdtt::sink_io_format::csv);
    const std::vector<lsdtt::sink_candidate> from_fgb =
        lsdtt::read_sink_candidates(fgb_path, meta,
                                    lsdtt::sink_io_format::flatgeobuf);

    std::string difference;
    check("csv_fgb_record_equality: .csv and .fgb parse back to identical records",
          vectors_equal(from_csv, from_fgb, difference), difference);

    check("csv_fgb_record_equality: and both equal what was written",
          vectors_equal(sinks, from_csv, difference), difference);
}


void test_cross_format_round_trip(const std::string&                        dir,
                                  const std::vector<lsdtt::sink_candidate>& sinks,
                                  const lsdtt::lsdtt_raster_meta&           meta)
{
    lsdtt::lsdtt_data_loader loader;

    struct format_case
    {
        lsdtt::sink_io_format format;
        const char*           extension;
        const char*           name;
    };

    const format_case formats[3] = {
        { lsdtt::sink_io_format::csv,        ".csv",  "csv"  },
        { lsdtt::sink_io_format::flatgeobuf, ".fgb",  "fgb"  },
        { lsdtt::sink_io_format::gpkg,       ".gpkg", "gpkg" }
    };

    for (const format_case& from : formats)
    {
        for (const format_case& to : formats)
        {
            if (std::string(from.name) == to.name) continue;

            const std::string first =
                dir + "/cross_" + from.name + "_" + to.name + from.extension;
            const std::string second =
                dir + "/cross_" + from.name + "_" + to.name + to.extension;

            std::filesystem::remove(first);
            std::filesystem::remove(second);

            lsdtt::write_sink_candidates(loader, first, sinks, meta, from.format);

            const std::vector<lsdtt::sink_candidate> stage_one =
                lsdtt::read_sink_candidates(first, meta, from.format);

            lsdtt::write_sink_candidates(loader, second, stage_one, meta,
                                         to.format);

            const std::vector<lsdtt::sink_candidate> stage_two =
                lsdtt::read_sink_candidates(second, meta, to.format);

            std::string difference;
            check(std::string("cross_format_round_trip: ") + from.name +
                  " -> read -> " + to.name + " -> read is unchanged",
                  vectors_equal(sinks, stage_two, difference), difference);
        }
    }

    // The format may also be resolved from the extension, which is what the CLI
    // does, and must give the same answer.
    const std::string by_extension = dir + "/by_extension.fgb";
    lsdtt::write_sink_candidates(
        loader, by_extension, sinks, meta,
        lsdtt::resolve_sink_io_format(by_extension, "same_as_extension"));

    std::string difference;
    check("cross_format_round_trip: resolving the format from the extension "
          "gives the same records",
          vectors_equal(
              sinks,
              lsdtt::read_sink_candidates(
                  by_extension, meta,
                  lsdtt::resolve_sink_io_format(by_extension, "")),
              difference),
          difference);
}


void test_candidate_attribute_typing(const std::string&                        dir,
                                     const std::vector<lsdtt::sink_candidate>& sinks,
                                     const lsdtt::lsdtt_raster_meta&           meta)
{
    lsdtt::lsdtt_data_loader loader;

    const std::string path = dir + "/typing.fgb";
    lsdtt::write_sink_candidates(loader, path, sinks, meta,
                                 lsdtt::sink_io_format::flatgeobuf);

    lsdtt::lsdtt_data_loader reader;
    reader.open_vector(path);
    const std::vector<std::string> layer_names = reader.get_layer_names();
    const lsdtt::lsdtt_vector_layer layer =
        reader.read_vector_layer(layer_names.empty() ? std::string()
                                                     : layer_names.front());
    reader.close_vector();

    // THE RULE THAT MATTERS: no field is a string.  The attribute model offers
    // int64, double and string only, so ids, counts and flags are int64 and
    // everything else -- including elevations, which are float in the struct --
    // is double.  float/int are not options and must not be asserted.
    bool no_strings = true;
    for (const lsdtt::lsdtt_vector_field_type type : layer.schema.field_types)
    {
        if (type == lsdtt::lsdtt_vector_field_type::string_type) no_strings = false;
    }

    check("candidate_attribute_typing: no field is a string", no_strings);

    const std::vector<std::string> int64_fields = {
        "sink_id", "routing_sink_id", "pit_row", "pit_col", "pit_flat",
        "n_pixels", "tagged_n_pixels", "hierarchy_level", "parent_sink_id",
        "accepted"
    };
    const std::vector<std::string> double_fields = {
        "easting", "northing", "latitude", "longitude",
        "pit_elevation", "spill_elevation", "depth_m", "area_m2", "volume_m3",
        "tagged_drainage_area_m2"
    };

    bool typing_ok = true;
    std::string wrong;

    for (const std::string& name : int64_fields)
    {
        const int index = layer.schema.index_of(name);
        if (index < 0 ||
            layer.schema.field_types[static_cast<std::size_t>(index)] !=
                lsdtt::lsdtt_vector_field_type::int64_type)
        {
            typing_ok = false;
            wrong += name + " ";
        }
    }

    for (const std::string& name : double_fields)
    {
        const int index = layer.schema.index_of(name);
        if (index < 0 ||
            layer.schema.field_types[static_cast<std::size_t>(index)] !=
                lsdtt::lsdtt_vector_field_type::double_type)
        {
            typing_ok = false;
            wrong += name + " ";
        }
    }

    check("candidate_attribute_typing: ids, counts and flags are int64; "
          "coordinates, elevations and sizes are double",
          typing_ok, wrong.empty() ? "" : ("wrong or missing: " + wrong));

    check("candidate_attribute_typing: every schema field is written",
          layer.schema.size() == 20 &&
          layer.records.size() == sinks.size(),
          std::to_string(layer.schema.size()) + " fields, " +
          std::to_string(layer.records.size()) + " features");
}


void test_candidate_position_precedence(const std::string&              dir,
                                        const lsdtt::lsdtt_raster_meta& meta)
{
    lsdtt::lsdtt_data_loader loader;

    // Four cells, one per tier, so that a feature carrying conflicting
    // positions can be resolved unambiguously.
    const int rows[4] = {  10,  40,  70, 120 };
    const int cols[4] = {  20,  50,  80, 130 };

    double x[4];
    double y[4];
    lon_lat ll[4];

    for (int i = 0; i < 4; ++i)
    {
        x[i]  = pixel_centre_x(meta, cols[i]);
        y[i]  = pixel_centre_y(meta, rows[i]);
        ll[i] = raster_xy_to_lon_lat(meta, x[i], y[i]);
    }

    // --- tier 1 beats tier 2: geometry beats easting / northing --------------
    lsdtt::lsdtt_vector_schema schema;
    schema.field_names = { "sink_id", "easting", "northing",
                           "latitude", "longitude", "pit_row", "pit_col" };
    schema.field_types = { lsdtt::lsdtt_vector_field_type::int64_type,
                           lsdtt::lsdtt_vector_field_type::double_type,
                           lsdtt::lsdtt_vector_field_type::double_type,
                           lsdtt::lsdtt_vector_field_type::double_type,
                           lsdtt::lsdtt_vector_field_type::double_type,
                           lsdtt::lsdtt_vector_field_type::int64_type,
                           lsdtt::lsdtt_vector_field_type::int64_type };

    const std::string geometry_path = dir + "/precedence_geometry.fgb";
    std::filesystem::remove(geometry_path);

    loader.begin_vector_write(geometry_path, "sink_candidates", schema,
                              meta.crs_wkt,
                              lsdtt::vector_output_format::flatgeobuf, false);

    // Geometry at cell 0; every attribute says somewhere else.
    loader.write_point(x[0], y[0],
                       { static_cast<std::int64_t>(0),
                         x[1], y[1], ll[2].lat, ll[2].lon,
                         static_cast<std::int64_t>(rows[3]),
                         static_cast<std::int64_t>(cols[3]) });

    loader.finish_vector_write();

    const std::vector<lsdtt::sink_candidate> by_geometry =
        lsdtt::read_sink_candidates(geometry_path, meta,
                                    lsdtt::sink_io_format::flatgeobuf);

    check("candidate_position_precedence: geometry beats easting/northing",
          by_geometry.size() == 1 &&
          by_geometry[0].pit_row == rows[0] && by_geometry[0].pit_col == cols[0],
          by_geometry.empty()
              ? "no record"
              : ("row " + std::to_string(by_geometry[0].pit_row) +
                 " col " + std::to_string(by_geometry[0].pit_col)));

    // A moved point moves the sink: the same layer with the geometry dragged to
    // cell 3 resolves to cell 3, whatever the attributes still say.
    const std::string moved_path = dir + "/precedence_moved.fgb";
    std::filesystem::remove(moved_path);

    loader.begin_vector_write(moved_path, "sink_candidates", schema,
                              meta.crs_wkt,
                              lsdtt::vector_output_format::flatgeobuf, false);
    loader.write_point(x[3], y[3],
                       { static_cast<std::int64_t>(0),
                         x[1], y[1], ll[2].lat, ll[2].lon,
                         static_cast<std::int64_t>(rows[0]),
                         static_cast<std::int64_t>(cols[0]) });
    loader.finish_vector_write();

    const std::vector<lsdtt::sink_candidate> moved =
        lsdtt::read_sink_candidates(moved_path, meta,
                                    lsdtt::sink_io_format::flatgeobuf);

    check("candidate_position_precedence: a point moved in the layer moves the sink",
          moved.size() == 1 &&
          moved[0].pit_row == rows[3] && moved[0].pit_col == cols[3]);

    // --- tier 2 beats tier 3: easting / northing beat latitude / longitude ---
    // A CSV has no geometry of its own, so tier 1 does not apply and the
    // easting / northing pair -- which is already in the raster CRS and is
    // therefore exact -- wins over a projected lat / lon.
    const std::string en_path = dir + "/precedence_easting.csv";
    write_csv_by_hand(
        en_path,
        { "sink_id", "easting", "northing", "latitude", "longitude",
          "pit_row", "pit_col" },
        { { "0", number_to_string(x[1]), number_to_string(y[1]),
            number_to_string(ll[2].lat), number_to_string(ll[2].lon),
            std::to_string(rows[3]), std::to_string(cols[3]) } });

    const std::vector<lsdtt::sink_candidate> by_easting =
        lsdtt::read_sink_candidates(en_path, meta, lsdtt::sink_io_format::csv);

    check("candidate_position_precedence: easting/northing beat latitude/longitude",
          by_easting.size() == 1 &&
          by_easting[0].pit_row == rows[1] && by_easting[0].pit_col == cols[1],
          by_easting.empty()
              ? "no record"
              : ("row " + std::to_string(by_easting[0].pit_row) +
                 " col " + std::to_string(by_easting[0].pit_col)));

    // --- tier 3 beats tier 4: latitude / longitude beat row / col -----------
    const std::string ll_path = dir + "/precedence_latlon.csv";
    write_csv_by_hand(
        ll_path,
        { "sink_id", "latitude", "longitude", "pit_row", "pit_col" },
        { { "0", number_to_string(ll[2].lat), number_to_string(ll[2].lon),
            std::to_string(rows[3]), std::to_string(cols[3]) } });

    const std::vector<lsdtt::sink_candidate> by_lonlat =
        lsdtt::read_sink_candidates(ll_path, meta, lsdtt::sink_io_format::csv);

    check("candidate_position_precedence: latitude/longitude beat row/col",
          by_lonlat.size() == 1 &&
          by_lonlat[0].pit_row == rows[2] && by_lonlat[0].pit_col == cols[2],
          by_lonlat.empty()
              ? "no record"
              : ("row " + std::to_string(by_lonlat[0].pit_row) +
                 " col " + std::to_string(by_lonlat[0].pit_col)));

    // --- tier 4 resolves on its own -----------------------------------------
    // Blank coordinates, so nothing above tier 4 can place the feature and the
    // row / col columns are all that is left.
    const std::string rc_path = dir + "/precedence_rowcol.csv";
    write_csv_by_hand(
        rc_path,
        { "sink_id", "latitude", "longitude", "pit_row", "pit_col" },
        { { "0", "", "", std::to_string(rows[3]), std::to_string(cols[3]) } });

    const std::vector<lsdtt::sink_candidate> by_rowcol =
        lsdtt::read_sink_candidates(rc_path, meta, lsdtt::sink_io_format::csv);

    check("candidate_position_precedence: row/col resolve when nothing else can",
          by_rowcol.size() == 1 &&
          by_rowcol[0].pit_row == rows[3] && by_rowcol[0].pit_col == cols[3],
          by_rowcol.empty()
              ? "no record"
              : ("row " + std::to_string(by_rowcol[0].pit_row) +
                 " col " + std::to_string(by_rowcol[0].pit_col)));

    // Whichever tier placed it, the position is snapped to the cell centre and
    // latitude / longitude are recomputed from that cell, so all four agree
    // with each other on everything positional.
    const bool snapped =
        by_rowcol.size() == 1 && moved.size() == 1 &&
        by_rowcol[0].pit_x == moved[0].pit_x &&
        by_rowcol[0].pit_y == moved[0].pit_y &&
        by_rowcol[0].latitude  == moved[0].latitude &&
        by_rowcol[0].longitude == moved[0].longitude &&
        by_rowcol[0].pit_flat  == moved[0].pit_flat;

    check("candidate_position_precedence: every tier snaps to the same cell centre",
          snapped);

    // A feature that says nothing about where it is cannot be placed and is
    // skipped rather than invented.
    const std::string nowhere_path = dir + "/precedence_nowhere.csv";
    write_csv_by_hand(nowhere_path,
                      { "sink_id", "latitude", "longitude" },
                      { { "0", "", "" } });

    check("candidate_position_precedence: an unplaceable feature is skipped",
          lsdtt::read_sink_candidates(
              nowhere_path, meta, lsdtt::sink_io_format::csv).empty());

    // ... and so is one that lands outside the raster.  read_sink_candidates
    // documents that "features that cannot be placed on the raster are
    // skipped", and its own bounds test
    //
    //     if (sink.pit_row < 0 || sink.pit_row >= meta.n_rows || ...) continue;
    //
    // says the same.  A reviewer who drags one point off the DEM must not
    // abort the read of the other 147.  The call is wrapped because a throw
    // here would otherwise take the rest of the driver with it.
    const std::string outside_path = dir + "/precedence_outside.csv";
    write_csv_by_hand(
        outside_path,
        { "sink_id", "easting", "northing", "latitude", "longitude" },
        { { "0", number_to_string(meta.x_min - 100000.0),
            number_to_string(meta.y_min - 100000.0), "", "" } });

    // MEASURED 2026-08-18, and reported as a defect rather than accommodated
    // here: read_sink_candidates resolves the position and then calls
    // col_from_x / row_from_y, which THROW std::out_of_range for a coordinate
    // outside the extent, so control never reaches the bounds test two lines
    // later.  One point dragged off the DEM by a reviewer therefore aborts the
    // read of the whole candidate file.  A point placed by GEOMETRY does not
    // hit this, because the loader marks it outside_raster_extent first; it is
    // the easting/northing and latitude/longitude tiers that throw.
    bool        skipped = false;
    std::string outside_detail;

    try
    {
        skipped = lsdtt::read_sink_candidates(
                      outside_path, meta, lsdtt::sink_io_format::csv).empty();
    }
    catch (const std::exception& e)
    {
        outside_detail = std::string("threw: ") + e.what();
    }

    check("candidate_position_precedence: a feature outside the raster is skipped",
          skipped, outside_detail);

    // The same thing one cell beyond each edge, so the report says which edge.
    const double just_west  = meta.x_min - 0.5 * meta.cell_size;
    const double inside_y   = pixel_centre_y(meta, 10);

    const std::string edge_path = dir + "/precedence_just_outside.csv";
    write_csv_by_hand(
        edge_path,
        { "sink_id", "easting", "northing", "latitude", "longitude" },
        { { "0", number_to_string(just_west), number_to_string(inside_y),
            "", "" } });

    bool        edge_skipped = false;
    std::string edge_detail;

    try
    {
        edge_skipped = lsdtt::read_sink_candidates(
                           edge_path, meta, lsdtt::sink_io_format::csv).empty();
    }
    catch (const std::exception& e)
    {
        edge_detail = std::string("threw: ") + e.what();
    }

    check("candidate_position_precedence: a feature half a cell west of the grid "
          "is skipped",
          edge_skipped, edge_detail);
}


// Writes a point layer with GDAL directly, so that a feature can carry NO
// attributes at all -- which is what a point hand-digitised in QGIS looks like.
void write_bare_geometry_layer(const lsdtt::lsdtt_raster_meta& meta,
                               const std::string&              path,
                               const std::vector<lsdtt::lsdtt_xy>& points,
                               bool                            with_sink_id,
                               const std::vector<int>&         sink_ids)
{
    std::filesystem::remove_all(path);

    GDALDriver* driver =
        GetGDALDriverManager()->GetDriverByName("FlatGeobuf");

    if (driver == nullptr)
    {
        throw std::runtime_error(
            "test_internal_drainage_detect: no FlatGeobuf driver");
    }

    GDALDataset* dataset =
        driver->Create(path.c_str(), 0, 0, 0, GDT_Unknown, nullptr);

    if (dataset == nullptr)
    {
        throw std::runtime_error(
            "test_internal_drainage_detect: could not create " + path);
    }

    OGRSpatialReference srs;
    srs.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);
    srs.importFromWkt(meta.crs_wkt.c_str());

    OGRLayer* layer =
        dataset->CreateLayer("sink_candidates", &srs, wkbPoint, nullptr);

    if (layer == nullptr)
    {
        GDALClose(dataset);
        throw std::runtime_error(
            "test_internal_drainage_detect: could not create the layer");
    }

    if (with_sink_id)
    {
        OGRFieldDefn field("sink_id", OFTInteger64);
        layer->CreateField(&field);
    }

    for (std::size_t i = 0; i < points.size(); ++i)
    {
        OGRFeature* feature = OGRFeature::CreateFeature(layer->GetLayerDefn());

        // Left UNSET when the caller supplies no id for this feature, which is
        // an absent attribute rather than a zero.
        if (with_sink_id && i < sink_ids.size())
        {
            feature->SetField("sink_id",
                              static_cast<GIntBig>(sink_ids[i]));
        }

        OGRPoint point(points[i].x, points[i].y);
        feature->SetGeometry(&point);

        if (layer->CreateFeature(feature) != OGRERR_NONE)
        {
            OGRFeature::DestroyFeature(feature);
            GDALClose(dataset);
            throw std::runtime_error(
                "test_internal_drainage_detect: could not write a feature");
        }

        OGRFeature::DestroyFeature(feature);
    }

    GDALClose(dataset);
}


void test_hand_added_feature(const std::string&              dir,
                             const lsdtt::lsdtt_raster_meta& meta)
{
    // A layer with no attribute table whatsoever: pure geometry.
    const std::string bare = dir + "/hand_added_bare.fgb";

    const std::vector<lsdtt::lsdtt_xy> points = {
        { pixel_centre_x(meta, 33), pixel_centre_y(meta, 22) }
    };

    write_bare_geometry_layer(meta, bare, points, false, {});

    const std::vector<lsdtt::sink_candidate> bare_back =
        lsdtt::read_sink_candidates(bare, meta, lsdtt::sink_io_format::flatgeobuf);

    check("hand_added_feature: a feature carrying only geometry is read",
          bare_back.size() == 1 &&
          bare_back[0].pit_row == 22 && bare_back[0].pit_col == 33,
          bare_back.empty()
              ? "no record"
              : ("row " + std::to_string(bare_back[0].pit_row) +
                 " col " + std::to_string(bare_back[0].pit_col)));

    check("hand_added_feature: it is accepted and given a sink_id",
          bare_back.size() == 1 && bare_back[0].accepted &&
          bare_back[0].sink_id == 0,
          bare_back.empty() ? "" : describe(bare_back[0]));

    // Mixed with reviewed candidates: the hand-added feature's id continues
    // from the largest one already in the file rather than colliding with it.
    const std::string mixed = dir + "/hand_added_mixed.fgb";

    const std::vector<lsdtt::lsdtt_xy> mixed_points = {
        { pixel_centre_x(meta, 40), pixel_centre_y(meta, 40) },
        { pixel_centre_x(meta, 50), pixel_centre_y(meta, 50) },
        { pixel_centre_x(meta, 60), pixel_centre_y(meta, 60) }
    };

    // Ids 7 and 9 for the first two; the third is left unset.
    write_bare_geometry_layer(meta, mixed, mixed_points, true, { 7, 9 });

    const std::vector<lsdtt::sink_candidate> mixed_back =
        lsdtt::read_sink_candidates(mixed, meta, lsdtt::sink_io_format::flatgeobuf);

    std::set<int> ids;
    for (const lsdtt::sink_candidate& s : mixed_back) ids.insert(s.sink_id);

    std::ostringstream got;
    for (const lsdtt::sink_candidate& s : mixed_back) got << s.sink_id << " ";

    check("hand_added_feature: a fresh id continues from the largest in the file",
          mixed_back.size() == 3 && ids.size() == 3 &&
          ids.find(7) != ids.end() && ids.find(9) != ids.end() &&
          ids.find(10) != ids.end(),
          "ids " + got.str());

    const std::string several_bare = dir + "/hand_added_several_bare.fgb";
    const std::vector<lsdtt::lsdtt_xy> several_points = {
        { pixel_centre_x(meta, 220), pixel_centre_y(meta, 200) },
        { pixel_centre_x(meta, 20), pixel_centre_y(meta, 30) }
    };
    write_bare_geometry_layer(meta, several_bare, several_points, false, {});
    const std::vector<lsdtt::sink_candidate> several_back =
        lsdtt::read_sink_candidates(
            several_bare, meta, lsdtt::sink_io_format::flatgeobuf);
    bool ids_follow_cells = several_back.size() == 2;
    for (const lsdtt::sink_candidate& sink : several_back)
    {
        const int expected = (sink.pit_row == 30 && sink.pit_col == 20) ? 0 : 1;
        if (sink.sink_id != expected) ids_follow_cells = false;
    }
    check("hand_added_feature: multiple missing ids are assigned by cell, not feature order",
          ids_follow_cells);
}


void test_candidate_crs_transform(const std::string&                        dir,
                                  const std::vector<lsdtt::sink_candidate>& sinks,
                                  const lsdtt::lsdtt_raster_meta&           meta)
{
    lsdtt::lsdtt_data_loader loader;

    // The same candidates, supplied in EPSG:4326 instead of the raster CRS --
    // which is what comes back from a reviewer who reprojected the layer, and
    // is the Greenland case in miniature (EPSG:3413 there, UTM 58S here).
    const std::string wgs84_path = dir + "/candidates_4326.fgb";
    std::filesystem::remove(wgs84_path);

    lsdtt::lsdtt_vector_schema schema;
    schema.field_names = { "sink_id" };
    schema.field_types = { lsdtt::lsdtt_vector_field_type::int64_type };

    loader.begin_vector_write(wgs84_path, "sink_candidates", schema,
                              epsg_4326_wkt(),
                              lsdtt::vector_output_format::flatgeobuf, false);

    for (const lsdtt::sink_candidate& s : sinks)
    {
        const lon_lat ll = raster_xy_to_lon_lat(meta, s.pit_x, s.pit_y);
        loader.write_point(ll.lon, ll.lat,
                           { static_cast<std::int64_t>(s.sink_id) });
    }

    loader.finish_vector_write();

    const std::vector<lsdtt::sink_candidate> transformed =
        lsdtt::read_sink_candidates(wgs84_path, meta,
                                    lsdtt::sink_io_format::flatgeobuf);

    bool same_cells = transformed.size() == sinks.size();
    std::string detail;

    for (std::size_t i = 0; i < transformed.size() && i < sinks.size(); ++i)
    {
        if (transformed[i].pit_flat != sinks[i].pit_flat)
        {
            same_cells = false;
            detail = "sink " + std::to_string(i) + ": flat " +
                     std::to_string(transformed[i].pit_flat) + " vs " +
                     std::to_string(sinks[i].pit_flat);
        }
    }

    check("candidate_crs_transform: an EPSG:4326 layer lands on the same pixels",
          same_cells, detail);

    // ... and the positional fields are then identical to the originals,
    // because every tier snaps to the cell centre.
    bool positions_identical = true;
    for (std::size_t i = 0; i < transformed.size() && i < sinks.size(); ++i)
    {
        if (transformed[i].pit_x     != sinks[i].pit_x ||
            transformed[i].pit_y     != sinks[i].pit_y ||
            transformed[i].latitude  != sinks[i].latitude ||
            transformed[i].longitude != sinks[i].longitude)
        {
            positions_identical = false;
        }
    }

    check("candidate_crs_transform: the snapped positions are bit-identical",
          positions_identical);
}


void test_format_resolution(const std::string&                        dir,
                            const std::vector<lsdtt::sink_candidate>& sinks,
                            const lsdtt::lsdtt_raster_meta&           meta)
{
    check("format_resolution: extensions map to formats",
          lsdtt::resolve_sink_io_format("a.csv", "") == lsdtt::sink_io_format::csv &&
          lsdtt::resolve_sink_io_format("a.fgb", "") ==
              lsdtt::sink_io_format::flatgeobuf &&
          lsdtt::resolve_sink_io_format("a.gpkg", "") ==
              lsdtt::sink_io_format::gpkg &&
          lsdtt::resolve_sink_io_format("a.FGB", "same_as_extension") ==
              lsdtt::sink_io_format::flatgeobuf);

    check("format_resolution: an explicit request overrides the extension",
          lsdtt::resolve_sink_io_format("a.csv", "flatgeobuf") ==
              lsdtt::sink_io_format::flatgeobuf &&
          lsdtt::resolve_sink_io_format("a.fgb", "csv") ==
              lsdtt::sink_io_format::csv &&
          lsdtt::resolve_sink_io_format("a.fgb", "geopackage") ==
              lsdtt::sink_io_format::gpkg);

    check_throws("format_resolution: an unknown extension is an error",
                 [] { (void)lsdtt::resolve_sink_io_format("sinks.shp", ""); },
                 "cannot resolve a sink file format");

    check_throws("format_resolution: a file with no extension is an error",
                 [] { (void)lsdtt::resolve_sink_io_format("sinks", ""); },
                 "cannot resolve a sink file format");

    check_throws("format_resolution: an unknown requested format is an error",
                 [] { (void)lsdtt::resolve_sink_io_format("sinks.csv", "geojson"); },
                 "unknown format");

    // The error names the offending string, because that is what the user has
    // to correct.
    check_throws("format_resolution: the message names the bad format",
                 [] { (void)lsdtt::resolve_sink_io_format("sinks.csv", "parquet"); },
                 "'parquet'");

    // Reading something that is not there is an error too, not an empty list.
    bool threw = false;
    try
    {
        (void)lsdtt::read_sink_candidates(dir + "/does_not_exist.fgb", meta,
                                          lsdtt::sink_io_format::flatgeobuf);
    }
    catch (const std::exception&)
    {
        threw = true;
    }

    check("format_resolution: reading a missing file throws", threw);

    (void)sinks;
}


void test_fgb_spatial_index_round_trip(const std::string&              dir,
                                       const lsdtt::lsdtt_raster_meta& meta)
{
    lsdtt::lsdtt_data_loader loader;

    // FlatGeobuf's packed Hilbert R-tree REORDERS features on write (plan
    // section 5.7: written 7, 8, 9 came back 9, 7, 8).  The reader is faithful;
    // it is the writer that reorders.  So the sink layer must be keyed on
    // sink_id and nothing may key on file position.
    //
    // The ids here are deliberately NON-SEQUENTIAL and NOT in position order,
    // so a reader that assigned ids by file position would be caught whether or
    // not the driver happened to reorder this particular layer.
    const int         ids[9]  = { 105, 3, 47, 12, 900, 6, 71, 28, 54 };
    const int         rows[9] = {  20, 140,  60, 200,  35, 175,  90, 250, 115 };
    const int         cols[9] = { 240,  30, 160,  75, 210,  95, 130, 265,  55 };

    std::vector<lsdtt::sink_candidate> sinks;

    for (int i = 0; i < 9; ++i)
    {
        lsdtt::sink_candidate s;
        s.sink_id         = ids[i];
        s.pit_row         = rows[i];
        s.pit_col         = cols[i];
        s.pit_flat        = flat_of(meta, rows[i], cols[i]);
        s.pit_x           = pixel_centre_x(meta, cols[i]);
        s.pit_y           = pixel_centre_y(meta, rows[i]);
        s.pit_elevation   = static_cast<float>(100 + i);
        s.spill_elevation = static_cast<float>(150 + i);
        s.depth_m         = s.spill_elevation - s.pit_elevation;
        s.area_m2         = 1000.0 * (i + 1);
        s.volume_m3       = 250.0 * (i + 1);
        s.n_pixels        = 10 * (i + 1);
        s.hierarchy_level = i % 2;
        s.parent_sink_id  = (i == 0) ? -1 : ids[0];
        s.accepted        = (i % 3) != 0;

        const lon_lat ll = raster_xy_to_lon_lat(meta, s.pit_x, s.pit_y);
        s.latitude  = ll.lat;
        s.longitude = ll.lon;

        sinks.push_back(s);
    }

    const std::string indexed = dir + "/spatial_index.fgb";
    std::filesystem::remove(indexed);

    lsdtt::write_sink_candidates(loader, indexed, sinks, meta,
                                 lsdtt::sink_io_format::flatgeobuf,
                                 "sink_candidates",
                                 true);   // SPATIAL INDEX ON

    const std::vector<lsdtt::sink_candidate> read_back =
        lsdtt::read_sink_candidates(indexed, meta,
                                    lsdtt::sink_io_format::flatgeobuf);

    check("fgb_spatial_index_round_trip: every feature comes back",
          read_back.size() == sinks.size(),
          std::to_string(read_back.size()) + " of " +
          std::to_string(sinks.size()));

    // Whether the driver reordered THIS layer is an observation about GDAL, not
    // a contract, so it is reported rather than asserted.
    std::ostringstream written_order;
    std::ostringstream read_order;
    bool reordered = false;

    for (std::size_t i = 0; i < sinks.size(); ++i)
    {
        written_order << sinks[i].sink_id << " ";
    }
    for (std::size_t i = 0; i < read_back.size(); ++i)
    {
        read_order << read_back[i].sink_id << " ";
        if (i < sinks.size() && read_back[i].sink_id != sinks[i].sink_id)
        {
            reordered = true;
        }
    }

    report("fgb_spatial_index_round_trip: feature order",
           "written " + written_order.str() + "-> read " + read_order.str() +
           (reordered ? "(REORDERED by the writer)" : "(order preserved here)"));

    // THE ASSERTION THAT MATTERS: keyed on sink_id, every record is recovered
    // field for field.  Nothing here indexes by file position.
    std::map<int, const lsdtt::sink_candidate*> by_id;
    for (const lsdtt::sink_candidate& s : read_back) by_id[s.sink_id] = &s;

    bool        all_recovered = true;
    std::string first_bad;

    for (const lsdtt::sink_candidate& original : sinks)
    {
        const auto found = by_id.find(original.sink_id);

        if (found == by_id.end())
        {
            all_recovered = false;
            first_bad = "sink_id " + std::to_string(original.sink_id) + " missing";
            break;
        }

        if (!records_equal(original, *found->second))
        {
            all_recovered = false;
            first_bad = describe(original) + " vs " + describe(*found->second);
            break;
        }
    }

    check("fgb_spatial_index_round_trip: every record is recovered, keyed on sink_id",
          all_recovered && by_id.size() == sinks.size(), first_bad);

    // The same layer written WITHOUT the index carries the same records, so the
    // index is a storage detail and not a semantic one.
    const std::string streamed = dir + "/no_spatial_index.fgb";
    std::filesystem::remove(streamed);

    lsdtt::write_sink_candidates(loader, streamed, sinks, meta,
                                 lsdtt::sink_io_format::flatgeobuf,
                                 "sink_candidates",
                                 false);

    const std::vector<lsdtt::sink_candidate> streamed_back =
        lsdtt::read_sink_candidates(streamed, meta,
                                    lsdtt::sink_io_format::flatgeobuf);

    bool same_set = streamed_back.size() == read_back.size();

    for (const lsdtt::sink_candidate& s : streamed_back)
    {
        const auto found = by_id.find(s.sink_id);
        if (found == by_id.end() || !records_equal(s, *found->second))
        {
            same_set = false;
        }
    }

    check("fgb_spatial_index_round_trip: indexed and streamed layers hold the "
          "same records",
          same_set);

    const std::vector<lsdtt::sink_candidate> indexed_accepted =
        lsdtt::accepted_sinks(read_back);
    const std::vector<lsdtt::sink_candidate> streamed_accepted =
        lsdtt::accepted_sinks(streamed_back);

    std::string accepted_difference;
    bool ids_stable = indexed_accepted.size() == 6;
    for (std::size_t i = 0; i < indexed_accepted.size(); ++i)
    {
        if (indexed_accepted[i].routing_sink_id != static_cast<int>(i))
            ids_stable = false;
        if (i > 0 &&
            indexed_accepted[i - 1].sink_id >= indexed_accepted[i].sink_id)
            ids_stable = false;
    }

    check("fgb_spatial_index_round_trip: routing ids are deterministic after reordering",
          ids_stable && vectors_equal(indexed_accepted, streamed_accepted,
                                      accepted_difference),
          accepted_difference);
}

} // namespace


int main(int argc, char** argv)
{
    const std::string output_dir =
        (argc > 1) ? std::string(argv[1])
                   : std::string("internal_drainage_detect_output");

    std::cout << "lsdtt3 internal-drainage unit tests "
                 "(phase 2b: detection and the review round-trip)\n\n";

    try
    {
        GDALAllRegister();
        std::filesystem::create_directories(output_dir);

        std::cout << "--- detection: depression hierarchy ---\n";
        test_detected_pits_match_the_oracle();
        test_detected_geometry_matches_the_oracle();
        test_flat_floor_tie_break_for_detection();
        test_size_filters();

        std::cout << "\n--- detection: nesting ---\n";
        test_nested_craters_hierarchy();
        test_shared_pit_cell_is_not_a_key();
        test_min_pixels_filters_a_leaf();

        std::cout << "\n--- detection: the legacy fill-difference path ---\n";
        test_fill_difference_agrees_on_a_simple_field();
        test_fill_difference_diverges_on_a_nested_field();
        test_fill_difference_min_pixels_is_legacy();

        std::cout << "\n--- depression footprints ---\n";
        test_depression_polygon_join();
        test_nested_sink_may_have_no_polygon();
        test_min_pixels_conserves_terrain();

        // One detection run feeds every round-trip test, which is the point:
        // the .csv and the .fgb written from THE SAME run must parse back to
        // identical records.
        const crater_field_fixture fx = make_simple_crater_field();
        const lsdtt::lsdtt_raster_meta& meta = fx.gt.dem.get_meta();

        const std::vector<lsdtt::sink_candidate> sinks =
            lsdtt::find_internal_drainage(fx.gt.dem,
                                          lsdtt::sink_detection_options{});

        std::cout << "\n--- the review round-trip ---\n";
        test_candidate_round_trip_csv(output_dir, sinks, meta);
        test_candidate_round_trip_fgb(output_dir, sinks, meta);
        test_csv_fgb_record_equality(output_dir, sinks, meta);
        test_cross_format_round_trip(output_dir, sinks, meta);
        test_candidate_attribute_typing(output_dir, sinks, meta);
        test_candidate_position_precedence(output_dir, meta);
        test_hand_added_feature(output_dir, meta);
        test_candidate_crs_transform(output_dir, sinks, meta);
        test_format_resolution(output_dir, sinks, meta);
        test_fgb_spatial_index_round_trip(output_dir, meta);
    }
    catch (const std::exception& e)
    {
        std::cout << "FAIL  unexpected exception: " << e.what() << std::endl;
        ++g_n_checks;
        ++g_n_failed;
    }

    std::cout << "\n" << g_n_checks - g_n_failed << " / " << g_n_checks
              << " checks passed\n";
    std::cout << "RESULT: " << (g_n_failed == 0 ? "PASS" : "FAIL") << "\n";

    return (g_n_failed == 0) ? 0 : 1;
}
