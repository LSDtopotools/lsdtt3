// test_divide_migration.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
//
// -----------------------------------------------------------------------------
// Unit driver for the divide-migration compute layer
// (src/cpu/lsdtt_divide_migration_cpu.cpp), plus the generator for the
// synthetic fixture used by the CLI end-to-end tests.
//
//   argv[1] (optional) - directory into which the CLI fixture (a three-valley
//                        DEM and the outlet point layer) is written. Without it
//                        only the in-memory unit checks run.
//
// The checks target the contract stated in lsdtt_divide_migration_cpu.hpp:
//
//   N1  the nearest divide point is the closest point on a SEGMENT, not the
//       closest vertex -- every fixture here puts the answer strictly between
//       two vertices, because a projection bug is invisible when it does not
//   N2  a head sitting on its divide yields a NaN gradient rather than an
//       infinity, and does not crash
//   N3  negative relief is a legitimate result and is passed through
//   N4  a head whose unit borders no divide keeps NaN fields and an empty
//       pair_id_nearest
//   N5  a nearest point on nodata leaves the elevation-derived fields NaN but
//       still reports distance_to_divide_m and pair_id_nearest
//   N6  the per-unit index: a head is NEVER matched to a divide between two
//       OTHER units, even when that divide is far closer in space. This is the
//       property the per-unit index design exists to guarantee
//   B1  apply_common_base_level_to_chi truncates below the base level, re-zeroes
//       chi at the root of EACH retained fragment, and leaves outlet_flat alone
//   M1  merge_divide_fragments joins two contiguous fragments into one polyline
//       with no duplicated join vertex and a RECOMPUTED length
//   M2  two DISJOINT fragments sharing a pair_id stay two features -- the
//       over-merge guard, which is the whole risk of this routine
//   M3  fragments supplied reversed and out of order still chain correctly
//   M4  a chain that closes on itself is a ring: it terminates and is closed
//   M5  fragments with different pair_ids never merge, even end to end
//   M6  the contiguity criterion is one cell diagonal: a gap inside it joins, a
//       gap outside it does not
//   P1  compute_pair_summaries: side means, per-metric counts and A-minus-B
//       deltas on a hand-computable fixture
//   P2  the A-minus-B sign convention: side A is ALWAYS the lower-numbered
//       unit, so swapping which unit carries which values flips every delta
//   P3  each of the four uncertainty methods against an independently computed
//       value (std_dev / std_err / ttest against scipy; bootstrap against an
//       exactly-zero degenerate case and against 1.96 * std_err asymptotically)
//   P4  NaN metrics are SKIPPED, not propagated: a side with a NaN relief gives
//       the same mean as the same side with that head removed
//   P5  a single-head side and an empty side give NaN uncertainty and do not
//       crash; a pair with an empty side still emits a row
//   P6  bootstrap determinism: same seed twice is identical, a different seed
//       differs, and an unrelated extra pair changes nothing
//   P7  heads with an empty pair_id_nearest take no part in any row
//   P8  a head is never counted on the wrong side of its pair
//   A1  adjacency-driven pair generation finds exactly the pairs the old
//       all-vs-all loop did, on a grid holding a NON-adjacent pair, and the
//       pair it skips is one that yielded no divide anyway
//   P9  the pair set comes from the DIVIDES: a divide with no matched heads
//       still emits a row (all statistics NaN), two segments sharing a pair_id
//       still give one row, and a head naming a pair no segment carries invents
//       no row
// -----------------------------------------------------------------------------

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include <ogr_spatialref.h>

#include "lsdtt_chi_extraction_cpu.hpp"
#include "lsdtt_data_loader.hpp"
#include "lsdtt_divide_migration_cpu.hpp"
#include "lsdtt_flow_base.hpp"
#include "lsdtt_junction_basin_manager_cpu.hpp"
#include "lsdtt_junction_network_cpu.hpp"
#include "lsdtt_outlet_clusters.hpp"
#include "lsdtt_raster_flow_preprocessing_cpu.hpp"
#include "lsdtt_raster_cpu.hpp"
#include "lsdtt_vector_utils.hpp"

namespace
{

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

bool close_to(double got, double expected, double tol)
{
    return std::isfinite(got) && std::fabs(got - expected) <= tol;
}

std::string num(double v)
{
    return std::to_string(v);
}

// -----------------------------------------------------------------------------
// A raster whose lower-left corner is at (0, 0) in projected units, so that
// test geometry can be written in plain metres.
// -----------------------------------------------------------------------------
lsdtt::lsdtt_raster_cpu make_uniform_elevation(
    int n_rows, int n_cols, double cell_size, float z)
{
    lsdtt::lsdtt_raster_meta meta;
    meta.n_rows = n_rows;
    meta.n_cols = n_cols;
    meta.cell_size = cell_size;
    meta.no_data_value = -9999.0;
    meta.set_default_synthetic_georeference();
    meta.x_min = 0.0;
    meta.y_min = 0.0;

    std::vector<float> data(static_cast<std::size_t>(n_rows) * n_cols, z);
    std::vector<std::uint8_t> mask(data.size(), 1);
    return lsdtt::lsdtt_raster_cpu(std::move(data), std::move(mask), meta);
}

lsdtt::lsdtt_divide_segment_record make_segment(
    const std::string& pair_id, int unit_a, int unit_b,
    std::vector<lsdtt::lsdtt_xy> points)
{
    lsdtt::lsdtt_divide_segment_record rec;
    rec.pair_id = pair_id;
    rec.unit_a = unit_a;
    rec.unit_b = unit_b;
    rec.segment_id = 1;
    rec.points = std::move(points);
    return rec;
}

// A divide fragment as the CLI builds one: a pair_id, the two source-basin
// junction numbers, and a bogus length_m so that any output length copied from
// the input rather than recomputed is immediately visible.
lsdtt::lsdtt_divide_segment_record make_fragment(
    const std::string& pair_id, int unit_a, int unit_b,
    std::int32_t junc_a, std::int32_t junc_b,
    std::vector<lsdtt::lsdtt_xy> points)
{
    lsdtt::lsdtt_divide_segment_record rec;
    rec.pair_id = pair_id;
    rec.unit_a = unit_a;
    rec.unit_b = unit_b;
    rec.basin_junction_a = junc_a;
    rec.basin_junction_b = junc_b;
    rec.segment_id = 99;
    rec.length_m = -12345.0;
    rec.points = std::move(points);
    return rec;
}

bool points_equal(const lsdtt::lsdtt_xy& a, double x, double y)
{
    return close_to(a.x, x, 1e-9) && close_to(a.y, y, 1e-9);
}

// True when no two consecutive vertices of a polyline coincide, i.e. no join
// vertex has been written twice.
bool has_no_repeated_vertex(const std::vector<lsdtt::lsdtt_xy>& points)
{
    for (std::size_t i = 1; i < points.size(); ++i)
    {
        if (points_equal(points[i], points[i - 1].x, points[i - 1].y))
        {
            return false;
        }
    }
    return true;
}

lsdtt::lsdtt_channel_head_metric_record make_head(
    int head_id, int unit_id, double x, double y, double z)
{
    lsdtt::lsdtt_channel_head_metric_record rec;
    rec.head_id = head_id;
    rec.outlet_id = unit_id;
    rec.unit_id = unit_id;
    rec.head_x = x;
    rec.head_y = y;
    rec.head_elevation_m = z;
    return rec;
}

// =============================================================================
// N1, N2, N3 -- projection, zero distance, negative relief
// =============================================================================
void test_nearest_point_geometry()
{
    // Elevation is a uniform 40 m plate, so every sampled divide elevation is
    // known without reference to the grid geometry.
    const lsdtt::lsdtt_raster_cpu elevation =
        make_uniform_elevation(40, 40, 10.0, 40.0f);

    // One horizontal divide, y = 100, running from x = 0 to x = 300. The two
    // vertices are deliberately far from every head below.
    std::vector<lsdtt::lsdtt_divide_segment_record> segments;
    segments.push_back(make_segment("1|2", 1, 2,
        {lsdtt::lsdtt_xy{0.0, 100.0}, lsdtt::lsdtt_xy{300.0, 100.0}}));

    std::vector<lsdtt::lsdtt_channel_head_metric_record> heads;
    heads.push_back(make_head(1, 1, 150.0, 60.0, 10.0));   // N1: mid-segment
    heads.push_back(make_head(2, 1, 150.0, 100.0, 10.0));  // N2: on the divide
    heads.push_back(make_head(3, 2, 150.0, 60.0, 100.0));  // N3: head above divide

    const std::size_t n_matched =
        lsdtt::compute_nearest_divide_metrics(segments, elevation, heads);

    check("N0 every head matched", n_matched == 3,
          "n_matched=" + std::to_string(n_matched));

    // N1: the answer is the foot of the perpendicular at x = 150, 40 m away.
    // Nearest-VERTEX logic would return (0,100) or (300,100) at 156.5 m.
    check("N1 nearest point is mid-segment, not a vertex",
          close_to(heads[0].nearest_divide_x, 150.0, 1e-9) &&
          close_to(heads[0].nearest_divide_y, 100.0, 1e-9),
          "got (" + num(heads[0].nearest_divide_x) + ", " +
                    num(heads[0].nearest_divide_y) + ")");
    check("N1 distance is the perpendicular distance",
          close_to(heads[0].distance_to_divide_m, 40.0, 1e-9),
          num(heads[0].distance_to_divide_m));
    check("N1 pair_id is the owning segment", heads[0].pair_id_nearest == "1|2",
          heads[0].pair_id_nearest);
    check("N1 relief and gradient from the sampled divide elevation",
          close_to(heads[0].relief_m, 30.0, 1e-9) &&
          close_to(heads[0].gradient_to_divide, 0.75, 1e-9),
          "relief=" + num(heads[0].relief_m) +
          " gradient=" + num(heads[0].gradient_to_divide));

    // N2: zero distance must not divide by zero.
    check("N2 zero distance is reported as zero",
          close_to(heads[1].distance_to_divide_m, 0.0, 1e-12),
          num(heads[1].distance_to_divide_m));
    check("N2 zero distance yields a NaN gradient, not an infinity",
          std::isnan(heads[1].gradient_to_divide),
          num(heads[1].gradient_to_divide));
    check("N2 relief is still reported at zero distance",
          close_to(heads[1].relief_m, 30.0, 1e-9), num(heads[1].relief_m));

    // N3: a head above its divide is a real result, not an error to clamp.
    check("N3 negative relief is preserved, not clamped or absolutised",
          close_to(heads[2].relief_m, -60.0, 1e-9), num(heads[2].relief_m));
    check("N3 negative gradient is preserved",
          close_to(heads[2].gradient_to_divide, -1.5, 1e-9),
          num(heads[2].gradient_to_divide));
}

// =============================================================================
// N4 -- a unit that borders no divide
// =============================================================================
void test_unit_without_divide()
{
    const lsdtt::lsdtt_raster_cpu elevation =
        make_uniform_elevation(40, 40, 10.0, 40.0f);

    std::vector<lsdtt::lsdtt_divide_segment_record> segments;
    segments.push_back(make_segment("1|2", 1, 2,
        {lsdtt::lsdtt_xy{0.0, 100.0}, lsdtt::lsdtt_xy{300.0, 100.0}}));

    std::vector<lsdtt::lsdtt_channel_head_metric_record> heads;
    heads.push_back(make_head(1, 7, 150.0, 60.0, 10.0));   // unit 7 has no divide
    heads.push_back(make_head(2, 2, 150.0, 60.0, 10.0));

    const std::size_t n_matched =
        lsdtt::compute_nearest_divide_metrics(segments, elevation, heads);

    check("N4 only the head whose unit borders a divide is matched",
          n_matched == 1, "n_matched=" + std::to_string(n_matched));
    check("N4 unmatched head keeps an empty pair_id_nearest",
          heads[0].pair_id_nearest.empty(), heads[0].pair_id_nearest);
    check("N4 unmatched head keeps NaN nearest-divide fields",
          std::isnan(heads[0].nearest_divide_x) &&
          std::isnan(heads[0].nearest_divide_y) &&
          std::isnan(heads[0].divide_elevation_m) &&
          std::isnan(heads[0].distance_to_divide_m) &&
          std::isnan(heads[0].relief_m) &&
          std::isnan(heads[0].gradient_to_divide));
    check("N4 the neighbouring unit is unaffected",
          heads[1].pair_id_nearest == "1|2" &&
          close_to(heads[1].distance_to_divide_m, 40.0, 1e-9));
}

// =============================================================================
// N5 -- nearest point on nodata
// =============================================================================
void test_nearest_point_on_nodata()
{
    // 40 x 40 at 10 m from (0,0), so y_max = 400. The divide is put at y = 105
    // rather than on a cell boundary, so the cell containing the nearest point
    // (155, 105) is unambiguously col = 15, row = floor((400-105)/10) = 29.
    lsdtt::lsdtt_raster_meta meta;
    meta.n_rows = 40;
    meta.n_cols = 40;
    meta.cell_size = 10.0;
    meta.no_data_value = -9999.0;
    meta.set_default_synthetic_georeference();
    meta.x_min = 0.0;
    meta.y_min = 0.0;

    std::vector<float> data(1600, 40.0f);
    std::vector<std::uint8_t> mask(1600, 1);
    const std::size_t hole = static_cast<std::size_t>(29) * 40 + 15;
    data[hole] = -9999.0f;
    mask[hole] = 0;
    lsdtt::lsdtt_raster_cpu elevation(std::move(data), std::move(mask), meta);

    std::vector<lsdtt::lsdtt_divide_segment_record> segments;
    segments.push_back(make_segment("1|2", 1, 2,
        {lsdtt::lsdtt_xy{0.0, 105.0}, lsdtt::lsdtt_xy{300.0, 105.0}}));

    std::vector<lsdtt::lsdtt_channel_head_metric_record> heads;
    heads.push_back(make_head(1, 1, 155.0, 60.0, 10.0));

    const std::size_t n_matched =
        lsdtt::compute_nearest_divide_metrics(segments, elevation, heads);

    check("N5 a nodata nearest cell still counts as matched", n_matched == 1,
          "n_matched=" + std::to_string(n_matched));
    check("N5 geometry is still reported over nodata",
          heads[0].pair_id_nearest == "1|2" &&
          close_to(heads[0].distance_to_divide_m, 45.0, 1e-9) &&
          close_to(heads[0].nearest_divide_x, 155.0, 1e-9),
          "pair=" + heads[0].pair_id_nearest +
          " d=" + num(heads[0].distance_to_divide_m));
    check("N5 elevation-derived fields are NaN over nodata",
          std::isnan(heads[0].divide_elevation_m) &&
          std::isnan(heads[0].relief_m) &&
          std::isnan(heads[0].gradient_to_divide),
          "z=" + num(heads[0].divide_elevation_m) +
          " relief=" + num(heads[0].relief_m));
}

// =============================================================================
// N6 -- three units: never match across a divide the head's unit does not share
// =============================================================================
void test_three_unit_isolation()
{
    const lsdtt::lsdtt_raster_cpu elevation =
        make_uniform_elevation(60, 120, 10.0, 40.0f);

    // Units 1 | 2 | 3 side by side. The 1|2 divide runs up x = 100 and the 2|3
    // divide up x = 900. Unit 1 and unit 3 share no divide at all.
    std::vector<lsdtt::lsdtt_divide_segment_record> segments;
    segments.push_back(make_segment("1|2", 1, 2,
        {lsdtt::lsdtt_xy{100.0, 0.0}, lsdtt::lsdtt_xy{100.0, 500.0}}));
    segments.push_back(make_segment("2|3", 2, 3,
        {lsdtt::lsdtt_xy{900.0, 0.0}, lsdtt::lsdtt_xy{900.0, 500.0}}));

    std::vector<lsdtt::lsdtt_channel_head_metric_record> heads;
    // A unit-3 head parked 10 m from the 1|2 divide. The globally nearest
    // geometry is the WRONG answer: unit 3 does not border that divide, so the
    // 2|3 divide 790 m away is the only admissible match.
    heads.push_back(make_head(1, 3, 110.0, 250.0, 10.0));
    // The mirror case for unit 1 against the 2|3 divide.
    heads.push_back(make_head(2, 1, 890.0, 250.0, 10.0));
    // A unit-2 head, which borders both, must take the genuinely nearest one.
    heads.push_back(make_head(3, 2, 200.0, 250.0, 10.0));

    const std::size_t n_matched =
        lsdtt::compute_nearest_divide_metrics(segments, elevation, heads);
    check("N6 all three heads matched", n_matched == 3,
          "n_matched=" + std::to_string(n_matched));

    check("N6 a unit-3 head is never matched to the 1|2 divide",
          heads[0].pair_id_nearest == "2|3" &&
          close_to(heads[0].distance_to_divide_m, 790.0, 1e-9),
          "pair=" + heads[0].pair_id_nearest +
          " d=" + num(heads[0].distance_to_divide_m));
    check("N6 a unit-1 head is never matched to the 2|3 divide",
          heads[1].pair_id_nearest == "1|2" &&
          close_to(heads[1].distance_to_divide_m, 790.0, 1e-9),
          "pair=" + heads[1].pair_id_nearest +
          " d=" + num(heads[1].distance_to_divide_m));
    check("N6 a shared unit takes the genuinely nearest of its own divides",
          heads[2].pair_id_nearest == "1|2" &&
          close_to(heads[2].distance_to_divide_m, 100.0, 1e-9),
          "pair=" + heads[2].pair_id_nearest +
          " d=" + num(heads[2].distance_to_divide_m));
}

// =============================================================================
// B1 -- common base level applied to chi
// =============================================================================
void test_common_base_level()
{
    // A 1-column, 8-row grid at 10 m. Flat index i is row i, so elevation is
    // set per node directly.
    lsdtt::lsdtt_raster_meta meta;
    meta.n_rows = 8;
    meta.n_cols = 1;
    meta.cell_size = 10.0;
    meta.no_data_value = -9999.0;
    meta.set_default_synthetic_georeference();
    meta.x_min = 0.0;
    meta.y_min = 0.0;

    std::vector<float> data = {0.0f, 10.0f, 20.0f, 30.0f,
                               40.0f, 50.0f, 60.0f, 70.0f};
    std::vector<std::uint8_t> mask(8, 1);
    lsdtt::lsdtt_raster_cpu elevation(std::move(data), std::move(mask), meta);

    lsdtt::lsdtt_chi_collection_result chi;

    // Basin A: a single chain 5 -> 4 -> 3 -> 2 -> 1 -> 0, chi = node index.
    lsdtt::lsdtt_chi_basin_result basin_a;
    basin_a.outlet_flat = 0;
    basin_a.outlet_chi = 0.0f;
    basin_a.flats = {0, 1, 2, 3, 4, 5};
    basin_a.chi = {0.0f, 1.0f, 2.0f, 3.0f, 4.0f, 5.0f};
    basin_a.receiver_flats = {-1, 0, 1, 2, 3, 4};
    chi.basins.push_back(basin_a);

    // Basin B: node 0 with TWO tributaries, so truncating just above node 0
    // leaves two disjoint retained fragments, each of which must be re-zeroed
    // at its own root. A single global offset would fail this.
    //   0(z=0) <- 4(z=40) <- 6(z=60)
    //   0(z=0) <- 5(z=50) <- 7(z=70)
    lsdtt::lsdtt_chi_basin_result basin_b;
    basin_b.outlet_flat = 0;
    basin_b.outlet_chi = 0.0f;
    basin_b.flats = {0, 4, 5, 6, 7};
    basin_b.chi = {0.0f, 1.0f, 2.0f, 3.0f, 5.0f};
    basin_b.receiver_flats = {-1, 0, 0, 4, 5};
    chi.basins.push_back(basin_b);

    const std::size_t n_removed =
        lsdtt::apply_common_base_level_to_chi(chi, elevation, 35.0);

    // Basin A keeps z >= 35, i.e. nodes 4 and 5 (three removed);
    // basin B keeps 4, 5, 6, 7 (one removed).
    check("B1 removed-node count is the sum over basins", n_removed == 5,
          "n_removed=" + std::to_string(n_removed));

    const lsdtt::lsdtt_chi_basin_result& a = chi.basins[0];
    check("B1 basin A retains exactly the nodes at or above base level",
          a.flats.size() == 2 && a.flats[0] == 4 && a.flats[1] == 5,
          "n=" + std::to_string(a.flats.size()));
    check("B1 basin A chi is re-zeroed at the lowest retained node",
          a.chi.size() == 2 && close_to(a.chi[0], 0.0, 1e-6) &&
          close_to(a.chi[1], 1.0, 1e-6));
    check("B1 basin A outlet_flat is left untouched", a.outlet_flat == 0,
          std::to_string(a.outlet_flat));
    check("B1 basin A outlet_chi is reset to zero",
          close_to(a.outlet_chi, 0.0, 1e-9), num(a.outlet_chi));

    const lsdtt::lsdtt_chi_basin_result& b = chi.basins[1];
    check("B1 basin B retains all four upslope nodes",
          b.flats.size() == 4 && b.flats[0] == 4 && b.flats[1] == 5 &&
          b.flats[2] == 6 && b.flats[3] == 7,
          "n=" + std::to_string(b.flats.size()));
    // Fragment roots are nodes 4 (chi 1) and 5 (chi 2); each subtree is offset
    // by its OWN root, never by a single collection-wide value.
    check("B1 chi is re-zeroed per retained fragment, not globally",
          b.chi.size() == 4 &&
          close_to(b.chi[0], 0.0, 1e-6) && close_to(b.chi[1], 0.0, 1e-6) &&
          close_to(b.chi[2], 2.0, 1e-6) && close_to(b.chi[3], 3.0, 1e-6),
          "chi=" + num(b.chi[0]) + "," + num(b.chi[1]) + "," +
                   num(b.chi[2]) + "," + num(b.chi[3]));
    check("B1 basin B receivers are kept in step with the retained nodes",
          b.receiver_flats.size() == 4 && b.receiver_flats[0] == 0 &&
          b.receiver_flats[1] == 0 && b.receiver_flats[2] == 4 &&
          b.receiver_flats[3] == 5);
    check("B1 basin B outlet_flat is left untouched", b.outlet_flat == 0,
          std::to_string(b.outlet_flat));
}

// =============================================================================
// M1 -- two contiguous fragments become one polyline
// =============================================================================
void test_merge_two_contiguous_fragments()
{
    // 10 m cells, so the join tolerance is one cell diagonal, 14.142 m.
    // The two fragments share the vertex (100, 50) exactly.
    std::vector<lsdtt::lsdtt_divide_segment_record> fragments;
    fragments.push_back(make_fragment("1|2", 1, 2, 11, 22,
        {lsdtt::lsdtt_xy{0.0, 0.0}, lsdtt::lsdtt_xy{100.0, 0.0},
         lsdtt::lsdtt_xy{100.0, 50.0}}));
    fragments.push_back(make_fragment("1|2", 1, 2, 11, 33,
        {lsdtt::lsdtt_xy{100.0, 50.0}, lsdtt::lsdtt_xy{200.0, 50.0}}));

    const std::vector<lsdtt::lsdtt_divide_segment_record> merged =
        lsdtt::merge_divide_fragments(fragments, 10.0);

    check("M1 two contiguous fragments become one feature", merged.size() == 1,
          "n=" + std::to_string(merged.size()));
    if (merged.size() != 1) return;

    const lsdtt::lsdtt_divide_segment_record& m = merged[0];
    // 3 + 2 vertices minus the ONE shared vertex.
    check("M1 the join vertex is written once, not twice",
          m.points.size() == 4 && has_no_repeated_vertex(m.points),
          "n_vertices=" + std::to_string(m.points.size()));
    check("M1 the merged polyline runs end to end",
          m.points.size() == 4 &&
          points_equal(m.points.front(), 0.0, 0.0) &&
          points_equal(m.points.back(), 200.0, 50.0));
    // 100 + 50 + 100; the inputs carried a nonsense length_m, so summing them
    // or copying one of them cannot produce this.
    check("M1 length_m is recomputed from the merged geometry",
          close_to(m.length_m, 250.0, 1e-9), num(m.length_m));
    check("M1 n_fragments records how many fragments contributed",
          m.n_fragments == 2, std::to_string(m.n_fragments));
    check("M1 a merged record does not report one basin's junctions",
          m.basin_junction_a == -1 && m.basin_junction_b == -1,
          std::to_string(m.basin_junction_a) + "," +
          std::to_string(m.basin_junction_b));
    check("M1 segment_id is reassigned over the merged set", m.segment_id == 1,
          std::to_string(m.segment_id));
    check("M1 pair and unit ids survive the merge",
          m.pair_id == "1|2" && m.unit_a == 1 && m.unit_b == 2);
}

// =============================================================================
// M2 -- the over-merge guard: sharing a pair_id is not contiguity
// =============================================================================
void test_merge_leaves_disjoint_fragments_alone()
{
    // Two units can meet along two genuinely separate stretches of divide (they
    // wrap around a third unit). These two share a pair_id but are 400 m apart,
    // so joining them would invent a divide that does not exist.
    std::vector<lsdtt::lsdtt_divide_segment_record> fragments;
    fragments.push_back(make_fragment("1|2", 1, 2, 11, 22,
        {lsdtt::lsdtt_xy{0.0, 0.0}, lsdtt::lsdtt_xy{100.0, 0.0}}));
    fragments.push_back(make_fragment("1|2", 1, 2, 44, 55,
        {lsdtt::lsdtt_xy{500.0, 0.0}, lsdtt::lsdtt_xy{600.0, 0.0}}));

    const std::vector<lsdtt::lsdtt_divide_segment_record> merged =
        lsdtt::merge_divide_fragments(fragments, 10.0);

    check("M2 disjoint fragments sharing a pair_id stay two features",
          merged.size() == 2, "n=" + std::to_string(merged.size()));
    if (merged.size() != 2) return;

    check("M2 neither feature grew any geometry",
          merged[0].points.size() == 2 && merged[1].points.size() == 2 &&
          close_to(merged[0].length_m, 100.0, 1e-9) &&
          close_to(merged[1].length_m, 100.0, 1e-9),
          "len=" + num(merged[0].length_m) + "," + num(merged[1].length_m));
    check("M2 an unmerged feature keeps its source junction numbers",
          merged[0].n_fragments == 1 && merged[0].basin_junction_a == 11 &&
          merged[0].basin_junction_b == 22 && merged[1].basin_junction_a == 44 &&
          merged[1].basin_junction_b == 55);
    check("M2 segment_ids stay unique and dense",
          merged[0].segment_id == 1 && merged[1].segment_id == 2,
          std::to_string(merged[0].segment_id) + "," +
          std::to_string(merged[1].segment_id));
}

// =============================================================================
// M3 -- any order, either orientation
// =============================================================================
void test_merge_reversed_and_out_of_order()
{
    // The intended line is (0,0) -> (0,100) -> (100,100) -> (100,300). It is
    // supplied last-piece-first, with the first piece reversed, so a merge that
    // assumed input order or a consistent orientation cannot rebuild it.
    std::vector<lsdtt::lsdtt_divide_segment_record> fragments;
    fragments.push_back(make_fragment("1|2", 1, 2, 11, 22,
        {lsdtt::lsdtt_xy{100.0, 100.0}, lsdtt::lsdtt_xy{100.0, 300.0}}));
    fragments.push_back(make_fragment("1|2", 1, 2, 11, 33,
        {lsdtt::lsdtt_xy{0.0, 100.0}, lsdtt::lsdtt_xy{0.0, 0.0}}));
    fragments.push_back(make_fragment("1|2", 1, 2, 11, 44,
        {lsdtt::lsdtt_xy{0.0, 100.0}, lsdtt::lsdtt_xy{100.0, 100.0}}));

    const std::vector<lsdtt::lsdtt_divide_segment_record> merged =
        lsdtt::merge_divide_fragments(fragments, 10.0);

    check("M3 all three fragments chain into one feature", merged.size() == 1,
          "n=" + std::to_string(merged.size()));
    if (merged.size() != 1) return;

    const lsdtt::lsdtt_divide_segment_record& m = merged[0];
    check("M3 three fragments give 4 vertices, not 6",
          m.points.size() == 4 && has_no_repeated_vertex(m.points),
          "n_vertices=" + std::to_string(m.points.size()));
    check("M3 the chain runs from one free end to the other",
          m.points.size() == 4 &&
          ((points_equal(m.points.front(), 0.0, 0.0) &&
            points_equal(m.points.back(), 100.0, 300.0)) ||
           (points_equal(m.points.front(), 100.0, 300.0) &&
            points_equal(m.points.back(), 0.0, 0.0))));
    check("M3 the chained length is the true along-divide length",
          close_to(m.length_m, 400.0, 1e-9), num(m.length_m));
    check("M3 n_fragments counts all three", m.n_fragments == 3,
          std::to_string(m.n_fragments));
}

// =============================================================================
// M4 -- a chain that closes on itself
// =============================================================================
void test_merge_closed_ring()
{
    // The four sides of a 100 m square, shuffled and half of them reversed. A
    // ring is the case where a naive "keep extending until nothing matches"
    // loop can run round for ever, so this test is about termination as much as
    // about geometry.
    std::vector<lsdtt::lsdtt_divide_segment_record> fragments;
    fragments.push_back(make_fragment("1|2", 1, 2, 11, 22,
        {lsdtt::lsdtt_xy{0.0, 0.0}, lsdtt::lsdtt_xy{100.0, 0.0}}));
    fragments.push_back(make_fragment("1|2", 1, 2, 11, 33,
        {lsdtt::lsdtt_xy{100.0, 100.0}, lsdtt::lsdtt_xy{100.0, 0.0}}));
    fragments.push_back(make_fragment("1|2", 1, 2, 11, 44,
        {lsdtt::lsdtt_xy{100.0, 100.0}, lsdtt::lsdtt_xy{0.0, 100.0}}));
    fragments.push_back(make_fragment("1|2", 1, 2, 11, 55,
        {lsdtt::lsdtt_xy{0.0, 0.0}, lsdtt::lsdtt_xy{0.0, 100.0}}));

    const std::vector<lsdtt::lsdtt_divide_segment_record> merged =
        lsdtt::merge_divide_fragments(fragments, 10.0);

    check("M4 a ring is one feature", merged.size() == 1,
          "n=" + std::to_string(merged.size()));
    if (merged.size() != 1) return;

    const lsdtt::lsdtt_divide_segment_record& m = merged[0];
    check("M4 the ring closes: 5 vertices with the first repeated at the end",
          m.points.size() == 5 &&
          points_equal(m.points.front(), m.points.back().x, m.points.back().y),
          "n_vertices=" + std::to_string(m.points.size()));
    check("M4 no interior vertex is duplicated",
          m.points.size() == 5 &&
          has_no_repeated_vertex(m.points));
    check("M4 the ring length is the full perimeter",
          close_to(m.length_m, 400.0, 1e-9), num(m.length_m));
    check("M4 every fragment was consumed exactly once", m.n_fragments == 4,
          std::to_string(m.n_fragments));
}

// =============================================================================
// M5 -- pair_id is a hard boundary
// =============================================================================
void test_merge_never_crosses_pair_id()
{
    // These two fragments meet at (100, 0) exactly, but they are divides
    // between different unit pairs -- a triple point, not one divide.
    std::vector<lsdtt::lsdtt_divide_segment_record> fragments;
    fragments.push_back(make_fragment("1|2", 1, 2, 11, 22,
        {lsdtt::lsdtt_xy{0.0, 0.0}, lsdtt::lsdtt_xy{100.0, 0.0}}));
    fragments.push_back(make_fragment("1|3", 1, 3, 11, 33,
        {lsdtt::lsdtt_xy{100.0, 0.0}, lsdtt::lsdtt_xy{200.0, 0.0}}));
    fragments.push_back(make_fragment("2|3", 2, 3, 22, 33,
        {lsdtt::lsdtt_xy{100.0, 0.0}, lsdtt::lsdtt_xy{100.0, 100.0}}));

    const std::vector<lsdtt::lsdtt_divide_segment_record> merged =
        lsdtt::merge_divide_fragments(fragments, 10.0);

    check("M5 fragments of different pairs meeting at a point never merge",
          merged.size() == 3, "n=" + std::to_string(merged.size()));
    if (merged.size() != 3) return;

    check("M5 each pair keeps its own geometry and junctions",
          merged[0].pair_id == "1|2" && merged[1].pair_id == "1|3" &&
          merged[2].pair_id == "2|3" &&
          merged[0].n_fragments == 1 && merged[1].n_fragments == 1 &&
          merged[2].n_fragments == 1 &&
          merged[2].basin_junction_a == 22 && merged[2].basin_junction_b == 33);
}

// =============================================================================
// M6 -- the contiguity criterion is exactly one cell diagonal
// =============================================================================
void test_merge_tolerance_is_one_cell_diagonal()
{
    // 10 m cells: the tolerance is 14.142 m. A 10 m gap is inside it (two
    // outlines stopping on slightly different nodes at a triple point); a 20 m
    // gap has a whole cell between the ends and is real separation.
    const auto build = [](double gap_x) {
        std::vector<lsdtt::lsdtt_divide_segment_record> fragments;
        fragments.push_back(make_fragment("1|2", 1, 2, 11, 22,
            {lsdtt::lsdtt_xy{0.0, 0.0}, lsdtt::lsdtt_xy{100.0, 0.0}}));
        fragments.push_back(make_fragment("1|2", 1, 2, 11, 33,
            {lsdtt::lsdtt_xy{gap_x, 0.0}, lsdtt::lsdtt_xy{200.0, 0.0}}));
        return fragments;
    };

    const std::vector<lsdtt::lsdtt_divide_segment_record> inside =
        lsdtt::merge_divide_fragments(build(110.0), 10.0);
    check("M6 a gap inside one cell diagonal is contiguous",
          inside.size() == 1 && inside[0].points.size() == 3,
          "n=" + std::to_string(inside.size()));

    const std::vector<lsdtt::lsdtt_divide_segment_record> outside =
        lsdtt::merge_divide_fragments(build(120.0), 10.0);
    check("M6 a gap wider than one cell diagonal is not contiguous",
          outside.size() == 2, "n=" + std::to_string(outside.size()));

    // The same 20 m gap on a 100 m grid is well inside one cell diagonal, so
    // the criterion really is expressed in cells and not in fixed metres.
    const std::vector<lsdtt::lsdtt_divide_segment_record> coarse =
        lsdtt::merge_divide_fragments(build(120.0), 100.0);
    check("M6 the criterion scales with cell size, not fixed metres",
          coarse.size() == 1, "n=" + std::to_string(coarse.size()));
}

// =============================================================================
// Pair-summary helpers
//
// The expected values below are computed either by hand or with an INDEPENDENT
// implementation (scipy 1.15.3 / numpy 2.2.6), never read back out of the code
// under test. The local mean / sd helpers here exist so that the asymptotic
// bootstrap check does not lean on the library either.
// =============================================================================
constexpr double k_test_nan = std::numeric_limits<double>::quiet_NaN();

lsdtt::lsdtt_channel_head_metric_record make_summary_head(
    int head_id, int unit_id, const std::string& pair_id,
    double elevation, double gradient, double relief, double chi)
{
    lsdtt::lsdtt_channel_head_metric_record rec;
    rec.head_id = head_id;
    rec.outlet_id = unit_id;
    rec.unit_id = unit_id;
    rec.pair_id_nearest = pair_id;
    rec.head_elevation_m = elevation;
    rec.gradient_to_divide = gradient;
    rec.relief_m = relief;
    rec.chi = chi;
    return rec;
}

// The divide layer a head fixture implies: one segment per distinct non-empty
// pair_id_nearest. The pair set of a summary now comes from the segments, so a
// fixture written as heads alone has to be handed the divides those heads were
// matched to -- exactly what the CLI passes, where every pair_id_nearest was
// read off one of these same segments.
std::vector<lsdtt::lsdtt_divide_segment_record> segments_for(
    const std::vector<lsdtt::lsdtt_channel_head_metric_record>& heads)
{
    std::vector<lsdtt::lsdtt_divide_segment_record> segments;
    std::vector<std::string> seen;
    for (const auto& h : heads)
    {
        if (h.pair_id_nearest.empty()) continue;
        if (std::find(seen.begin(), seen.end(), h.pair_id_nearest) != seen.end())
        {
            continue;
        }
        seen.push_back(h.pair_id_nearest);

        const std::size_t bar = h.pair_id_nearest.find('|');
        const int unit_a = std::stoi(h.pair_id_nearest.substr(0, bar));
        const int unit_b = std::stoi(h.pair_id_nearest.substr(bar + 1));
        segments.push_back(make_segment(
            h.pair_id_nearest, unit_a, unit_b,
            {lsdtt::lsdtt_xy{0.0, 0.0}, lsdtt::lsdtt_xy{100.0, 0.0}}));
    }
    return segments;
}

const lsdtt::lsdtt_pair_summary_record* find_pair(
    const std::vector<lsdtt::lsdtt_pair_summary_record>& rows,
    const std::string& pair_id)
{
    for (const auto& r : rows)
    {
        if (r.pair_id == pair_id) return &r;
    }
    return nullptr;
}

double local_mean(const std::vector<double>& v)
{
    double s = 0.0;
    for (double x : v) s += x;
    return s / static_cast<double>(v.size());
}

double local_sd(const std::vector<double>& v)
{
    const double m = local_mean(v);
    double ss = 0.0;
    for (double x : v) ss += (x - m) * (x - m);
    return std::sqrt(ss / static_cast<double>(v.size() - 1));
}

// The P1/P2/P3 fixture, as a pair "1|2" with side A on unit 1.
//
//   side A (unit 1): elevation 100, 200, 300   gradient 0.1, 0.2, 0.3
//                    relief    10,  20,  30    chi      1,   2,   3
//   side B (unit 2): elevation 400, 700        gradient 0.5, 0.7
//                    relief     5,  15         chi      10,  20
//
// Every mean is a whole number and every delta is exact in binary except the
// gradients, which are checked to 1e-12.
std::vector<lsdtt::lsdtt_channel_head_metric_record> make_summary_fixture(
    int unit_for_side_a, int unit_for_side_b)
{
    std::vector<lsdtt::lsdtt_channel_head_metric_record> heads;
    const std::string pair =
        lsdtt::make_unit_pair_id(unit_for_side_a, unit_for_side_b);
    heads.push_back(make_summary_head(1, unit_for_side_a, pair, 100.0, 0.1, 10.0,  1.0));
    heads.push_back(make_summary_head(2, unit_for_side_a, pair, 200.0, 0.2, 20.0,  2.0));
    heads.push_back(make_summary_head(3, unit_for_side_a, pair, 300.0, 0.3, 30.0,  3.0));
    heads.push_back(make_summary_head(4, unit_for_side_b, pair, 400.0, 0.5,  5.0, 10.0));
    heads.push_back(make_summary_head(5, unit_for_side_b, pair, 700.0, 0.7, 15.0, 20.0));
    return heads;
}

// =============================================================================
// P1 -- side means, per-metric counts and A-minus-B deltas
// =============================================================================
void test_pair_summary_means_and_deltas()
{
    const std::vector<lsdtt::lsdtt_channel_head_metric_record> heads =
        make_summary_fixture(1, 2);

    lsdtt::lsdtt_pair_summary_options options;
    options.uncertainty_method = lsdtt::lsdtt_pair_uncertainty_method::std_err;

    const std::vector<lsdtt::lsdtt_pair_summary_record> rows =
        lsdtt::compute_pair_summaries(segments_for(heads), heads, options);

    check("P1 one pair gives one row", rows.size() == 1,
          "n=" + std::to_string(rows.size()));
    if (rows.size() != 1) return;

    const lsdtt::lsdtt_pair_summary_record& r = rows[0];
    check("P1 the row is tagged with its pair and both units",
          r.pair_id == "1|2" && r.unit_a == 1 && r.unit_b == 2,
          r.pair_id + " " + std::to_string(r.unit_a) + "," +
          std::to_string(r.unit_b));
    check("P1 heads are counted on the side their unit_id names",
          r.n_heads_a == 3 && r.n_heads_b == 2,
          std::to_string(r.n_heads_a) + "," + std::to_string(r.n_heads_b));

    // (100+200+300)/3 = 200; (400+700)/2 = 550; 200 - 550 = -350.
    check("P1 head elevation means and the A-B delta",
          close_to(r.mean_head_elevation_a, 200.0, 1e-12) &&
          close_to(r.mean_head_elevation_b, 550.0, 1e-12) &&
          close_to(r.delta_head_elevation, -350.0, 1e-12),
          num(r.mean_head_elevation_a) + " " + num(r.mean_head_elevation_b) +
          " " + num(r.delta_head_elevation));
    // (0.1+0.2+0.3)/3 = 0.2; (0.5+0.7)/2 = 0.6; 0.2 - 0.6 = -0.4.
    check("P1 gradient means and the A-B delta",
          close_to(r.mean_gradient_a, 0.2, 1e-12) &&
          close_to(r.mean_gradient_b, 0.6, 1e-12) &&
          close_to(r.delta_gradient, -0.4, 1e-12),
          num(r.mean_gradient_a) + " " + num(r.mean_gradient_b) + " " +
          num(r.delta_gradient));
    // (10+20+30)/3 = 20; (5+15)/2 = 10; 20 - 10 = +10. The only POSITIVE delta
    // in the fixture, so a routine that returned B-A everywhere fails here.
    check("P1 relief means and the A-B delta",
          close_to(r.mean_relief_a, 20.0, 1e-12) &&
          close_to(r.mean_relief_b, 10.0, 1e-12) &&
          close_to(r.delta_relief, 10.0, 1e-12),
          num(r.mean_relief_a) + " " + num(r.mean_relief_b) + " " +
          num(r.delta_relief));
    // (1+2+3)/3 = 2; (10+20)/2 = 15; 2 - 15 = -13.
    check("P1 chi means and the A-B delta",
          close_to(r.mean_chi_a, 2.0, 1e-12) &&
          close_to(r.mean_chi_b, 15.0, 1e-12) &&
          close_to(r.delta_chi, -13.0, 1e-12),
          num(r.mean_chi_a) + " " + num(r.mean_chi_b) + " " + num(r.delta_chi));

    check("P1 every metric carries its own valid-value count",
          r.n_head_elevation_a == 3 && r.n_head_elevation_b == 2 &&
          r.n_gradient_a == 3 && r.n_gradient_b == 2 &&
          r.n_relief_a == 3 && r.n_relief_b == 2 &&
          r.n_chi_a == 3 && r.n_chi_b == 2);
    check("P1 the method used is recorded in the row",
          r.uncertainty_method == "std_err", r.uncertainty_method);
}

// =============================================================================
// P2 -- side A is ALWAYS the lower-numbered unit
// =============================================================================
void test_pair_summary_sign_convention()
{
    lsdtt::lsdtt_pair_summary_options options;

    // The same five heads, but now the LOW-elevation triple sits on unit 2 and
    // the high pair on unit 1. Side A is still unit 1, so every delta must come
    // out with the OPPOSITE sign to P1. A routine that called whichever side it
    // met first "A" would return the P1 signs here.
    const std::vector<lsdtt::lsdtt_channel_head_metric_record> swapped =
        make_summary_fixture(2, 1);

    const std::vector<lsdtt::lsdtt_pair_summary_record> rows =
        lsdtt::compute_pair_summaries(segments_for(swapped), swapped, options);

    check("P2 swapping the sides keeps one row on the same pair_id",
          rows.size() == 1 && rows[0].pair_id == "1|2",
          "n=" + std::to_string(rows.size()));
    if (rows.size() != 1) return;

    const lsdtt::lsdtt_pair_summary_record& r = rows[0];
    check("P2 unit_a is the lower-numbered unit whichever side carries what",
          r.unit_a == 1 && r.unit_b == 2 && r.n_heads_a == 2 && r.n_heads_b == 3,
          std::to_string(r.unit_a) + "," + std::to_string(r.unit_b) + " n=" +
          std::to_string(r.n_heads_a) + "," + std::to_string(r.n_heads_b));
    check("P2 the means follow the units, not the input order",
          close_to(r.mean_head_elevation_a, 550.0, 1e-12) &&
          close_to(r.mean_head_elevation_b, 200.0, 1e-12),
          num(r.mean_head_elevation_a) + " " + num(r.mean_head_elevation_b));
    check("P2 every delta flips sign relative to P1",
          close_to(r.delta_head_elevation, 350.0, 1e-12) &&
          close_to(r.delta_gradient, 0.4, 1e-12) &&
          close_to(r.delta_relief, -10.0, 1e-12) &&
          close_to(r.delta_chi, 13.0, 1e-12),
          num(r.delta_head_elevation) + " " + num(r.delta_gradient) + " " +
          num(r.delta_relief) + " " + num(r.delta_chi));

    // The convention must hold for units that are not 1 and 2, and for a pair
    // supplied with the higher unit first.
    std::vector<lsdtt::lsdtt_channel_head_metric_record> high_units;
    high_units.push_back(make_summary_head(1, 7, "2|7", 900.0, 0.9, 90.0, 30.0));
    high_units.push_back(make_summary_head(2, 7, "2|7", 950.0, 0.9, 90.0, 30.0));
    high_units.push_back(make_summary_head(3, 2, "2|7", 100.0, 0.1, 10.0,  1.0));
    high_units.push_back(make_summary_head(4, 2, "2|7", 150.0, 0.1, 10.0,  1.0));

    const std::vector<lsdtt::lsdtt_pair_summary_record> high_rows =
        lsdtt::compute_pair_summaries(
            segments_for(high_units), high_units, options);
    check("P2 side A is unit 2 in pair 2|7, so the delta is negative",
          high_rows.size() == 1 && high_rows[0].unit_a == 2 &&
          high_rows[0].unit_b == 7 &&
          close_to(high_rows[0].mean_head_elevation_a, 125.0, 1e-12) &&
          close_to(high_rows[0].delta_head_elevation, -800.0, 1e-12),
          high_rows.empty() ? "no rows" : num(high_rows[0].delta_head_elevation));
}

// =============================================================================
// P3 -- the four uncertainty methods
// =============================================================================
void test_pair_summary_uncertainty_methods()
{
    const std::vector<lsdtt::lsdtt_channel_head_metric_record> heads =
        make_summary_fixture(1, 2);

    // Elevations: A = {100, 200, 300} (sd 100), B = {400, 700} (sd 212.132034).
    // Expected values from scipy 1.15.3 / numpy 2.2.6:
    //   pooled sd  sqrt((2*100^2 + 1*212.132034^2)/3)      = 147.19601443879745
    //   se(delta)  sqrt(100^2/3 + 212.132034^2/2)          = 160.72751268321590
    //   Welch p    stats.ttest_ind(A, B, equal_var=False)  =   0.22510998427378
    lsdtt::lsdtt_pair_summary_options options;

    options.uncertainty_method = lsdtt::lsdtt_pair_uncertainty_method::std_dev;
    const auto sd_rows =
        lsdtt::compute_pair_summaries(segments_for(heads), heads, options);
    check("P3 std_dev is the pooled standard deviation of the two sides",
          sd_rows.size() == 1 &&
          close_to(sd_rows[0].uncertainty_head_elevation,
                   147.19601443879745, 1e-9) &&
          sd_rows[0].uncertainty_method == "std_dev",
          sd_rows.empty() ? "no rows"
                          : num(sd_rows[0].uncertainty_head_elevation));

    options.uncertainty_method = lsdtt::lsdtt_pair_uncertainty_method::std_err;
    const auto se_rows =
        lsdtt::compute_pair_summaries(segments_for(heads), heads, options);
    check("P3 std_err is the standard error of the DIFFERENCE of means",
          se_rows.size() == 1 &&
          close_to(se_rows[0].uncertainty_head_elevation,
                   160.72751268321590, 1e-9) &&
          se_rows[0].uncertainty_method == "std_err",
          se_rows.empty() ? "no rows"
                          : num(se_rows[0].uncertainty_head_elevation));
    // The two must not be the same number, or neither assertion means anything.
    check("P3 std_dev and std_err are genuinely different statistics",
          sd_rows.size() == 1 && se_rows.size() == 1 &&
          std::fabs(sd_rows[0].uncertainty_head_elevation -
                    se_rows[0].uncertainty_head_elevation) > 1.0);

    options.uncertainty_method = lsdtt::lsdtt_pair_uncertainty_method::ttest;
    const auto tt_rows =
        lsdtt::compute_pair_summaries(segments_for(heads), heads, options);
    check("P3 ttest reports the two-sided Welch p-value",
          tt_rows.size() == 1 &&
          close_to(tt_rows[0].uncertainty_head_elevation,
                   0.22510998427378373, 1e-9) &&
          tt_rows[0].uncertainty_method == "ttest",
          tt_rows.empty() ? "no rows"
                          : num(tt_rows[0].uncertainty_head_elevation));

    // A second, clearly separated t-test case, so the p-value is pinned at the
    // far end of the scale too and not just near 0.2. scipy:
    //   stats.ttest_ind([1..5], [101..105], equal_var=False)
    //   -> t = -100, df = 8, p = 1.1167803048255944e-13
    // Checked on a RELATIVE tolerance, since an absolute one is meaningless at
    // 1e-13.
    std::vector<lsdtt::lsdtt_channel_head_metric_record> separated;
    for (int i = 0; i < 5; ++i)
    {
        separated.push_back(make_summary_head(
            i + 1, 1, "1|2", 1.0 + i, 0.1, 1.0, 1.0));
        separated.push_back(make_summary_head(
            i + 6, 2, "1|2", 101.0 + i, 0.1, 1.0, 1.0));
    }
    const auto sep_rows =
        lsdtt::compute_pair_summaries(
            segments_for(separated), separated, options);
    const double sep_p = sep_rows.empty()
        ? k_test_nan : sep_rows[0].uncertainty_head_elevation;
    check("P3 a cleanly separated pair gives a tiny p-value",
          sep_rows.size() == 1 && std::isfinite(sep_p) &&
          std::fabs(sep_p / 1.1167803048255944e-13 - 1.0) < 1e-6,
          num(sep_p * 1e13) + "e-13");

    // Bootstrap, exact case: both sides constant, so EVERY resample gives the
    // same delta and the confidence interval has zero width. No tolerance
    // needed and no reliance on the RNG.
    options.uncertainty_method = lsdtt::lsdtt_pair_uncertainty_method::bootstrap;
    options.bootstrap_replicates = 200;
    std::vector<lsdtt::lsdtt_channel_head_metric_record> constant_sides;
    for (int i = 0; i < 4; ++i)
    {
        constant_sides.push_back(make_summary_head(
            i + 1, 1, "1|2", 5.0, 0.1, 1.0, 1.0));
        constant_sides.push_back(make_summary_head(
            i + 5, 2, "1|2", 2.0, 0.1, 1.0, 1.0));
    }
    const auto const_rows =
        lsdtt::compute_pair_summaries(
            segments_for(constant_sides), constant_sides, options);
    check("P3 bootstrap of two constant sides has exactly zero width",
          const_rows.size() == 1 &&
          const_rows[0].uncertainty_head_elevation == 0.0 &&
          close_to(const_rows[0].delta_head_elevation, 3.0, 1e-12),
          const_rows.empty() ? "no rows"
                             : num(const_rows[0].uncertainty_head_elevation));

    // Bootstrap, asymptotic case: with 100 values a side the 95% CI half-width
    // of the delta must land near 1.96 * the analytic standard error of the
    // difference. The expected value is computed here from the sample, not read
    // back out of the library.
    std::vector<lsdtt::lsdtt_channel_head_metric_record> large;
    std::vector<double> side_a;
    std::vector<double> side_b;
    for (int i = 0; i < 100; ++i)
    {
        const double a = 100.0 + static_cast<double>(i % 10) * 2.0;
        const double b = 130.0 + static_cast<double>(i % 10) * 2.0;
        side_a.push_back(a);
        side_b.push_back(b);
        large.push_back(make_summary_head(i + 1,   1, "1|2", a, 0.1, 1.0, 1.0));
        large.push_back(make_summary_head(i + 201, 2, "1|2", b, 0.1, 1.0, 1.0));
    }
    const double sa = local_sd(side_a);
    const double sb = local_sd(side_b);
    const double se_delta = std::sqrt(sa * sa / 100.0 + sb * sb / 100.0);
    const double expected_half_width = 1.959963985 * se_delta;

    options.bootstrap_replicates = 2000;
    const auto large_rows =
        lsdtt::compute_pair_summaries(segments_for(large), large, options);
    const double got = large_rows.empty()
        ? k_test_nan : large_rows[0].uncertainty_head_elevation;
    check("P3 bootstrap converges on the 95% interval of the delta",
          large_rows.size() == 1 && std::isfinite(got) &&
          std::fabs(got / expected_half_width - 1.0) < 0.15,
          "got=" + num(got) + " expected~" + num(expected_half_width));
    check("P3 bootstrap records the method it used",
          large_rows.size() == 1 && large_rows[0].uncertainty_method == "bootstrap",
          large_rows.empty() ? "no rows" : large_rows[0].uncertainty_method);
}

// =============================================================================
// P4 -- NaN metrics are skipped, never propagated
// =============================================================================
void test_pair_summary_skips_nan()
{
    lsdtt::lsdtt_pair_summary_options options;
    options.uncertainty_method = lsdtt::lsdtt_pair_uncertainty_method::std_err;

    // Side A has three heads, one of which has a NaN relief and a NaN gradient
    // (the real case: its nearest divide point landed on nodata). Its elevation
    // and chi are perfectly good and must still count. Its elevation (260) and
    // chi (5) are deliberately off the arithmetic mean of the other two, so
    // that a run which dropped the whole head gives a DIFFERENT answer for
    // those two metrics -- otherwise the last check below proves nothing.
    std::vector<lsdtt::lsdtt_channel_head_metric_record> with_nan;
    with_nan.push_back(make_summary_head(1, 1, "1|2", 100.0, 0.1, 10.0, 1.0));
    with_nan.push_back(make_summary_head(
        2, 1, "1|2", 260.0, k_test_nan, k_test_nan, 5.0));
    with_nan.push_back(make_summary_head(3, 1, "1|2", 300.0, 0.3, 30.0, 3.0));
    with_nan.push_back(make_summary_head(4, 2, "1|2", 400.0, 0.5,  5.0, 10.0));
    with_nan.push_back(make_summary_head(5, 2, "1|2", 700.0, 0.7, 15.0, 20.0));

    // The same data with the offending head simply removed.
    std::vector<lsdtt::lsdtt_channel_head_metric_record> without;
    without.push_back(make_summary_head(1, 1, "1|2", 100.0, 0.1, 10.0, 1.0));
    without.push_back(make_summary_head(3, 1, "1|2", 300.0, 0.3, 30.0, 3.0));
    without.push_back(make_summary_head(4, 2, "1|2", 400.0, 0.5,  5.0, 10.0));
    without.push_back(make_summary_head(5, 2, "1|2", 700.0, 0.7, 15.0, 20.0));

    const auto nan_rows =
        lsdtt::compute_pair_summaries(
            segments_for(with_nan), with_nan, options);
    const auto ref_rows =
        lsdtt::compute_pair_summaries(segments_for(without), without, options);
    check("P4 both fixtures give one row",
          nan_rows.size() == 1 && ref_rows.size() == 1);
    if (nan_rows.size() != 1 || ref_rows.size() != 1) return;

    const lsdtt::lsdtt_pair_summary_record& n = nan_rows[0];
    const lsdtt::lsdtt_pair_summary_record& r = ref_rows[0];

    check("P4 a NaN relief does not poison the mean: it equals the mean "
          "without that head",
          std::isfinite(n.mean_relief_a) &&
          close_to(n.mean_relief_a, r.mean_relief_a, 1e-12) &&
          close_to(n.mean_relief_a, 20.0, 1e-12),
          num(n.mean_relief_a) + " vs " + num(r.mean_relief_a));
    check("P4 the same holds for gradient and for the delta",
          close_to(n.mean_gradient_a, 0.2, 1e-12) &&
          close_to(n.delta_relief, r.delta_relief, 1e-12),
          num(n.mean_gradient_a) + " " + num(n.delta_relief));
    check("P4 the uncertainty is computed from the surviving values only",
          close_to(n.uncertainty_relief, r.uncertainty_relief, 1e-12),
          num(n.uncertainty_relief) + " vs " + num(r.uncertainty_relief));

    // The whole reason the schema needed extending: n_heads_a no longer
    // describes the relief sample.
    check("P4 n_heads_a counts heads, the per-metric n counts valid values",
          n.n_heads_a == 3 && n.n_relief_a == 2 && n.n_gradient_a == 2 &&
          n.n_head_elevation_a == 3 && n.n_chi_a == 3,
          std::to_string(n.n_heads_a) + " relief=" +
          std::to_string(n.n_relief_a) + " elev=" +
          std::to_string(n.n_head_elevation_a));
    // (100 + 260 + 300)/3 = 220 and (1 + 5 + 3)/3 = 3, against 200 and 2 for the
    // reference: dropping the whole head because one of its metrics was NaN
    // would give the reference numbers here.
    check("P4 metrics that were never NaN keep the full sample",
          close_to(n.mean_head_elevation_a, 220.0, 1e-12) &&
          close_to(n.mean_chi_a, 3.0, 1e-12) &&
          close_to(r.mean_head_elevation_a, 200.0, 1e-12) &&
          !close_to(n.mean_head_elevation_a, r.mean_head_elevation_a, 1e-12),
          num(n.mean_head_elevation_a) + " vs " + num(r.mean_head_elevation_a));
}

// =============================================================================
// P5 -- single-head and empty sides
// =============================================================================
void test_pair_summary_small_samples()
{
    lsdtt::lsdtt_pair_summary_options options;

    // One head on side A. The mean is defined; no method has a within-side
    // spread to work with, so every uncertainty must be NaN rather than 0.
    std::vector<lsdtt::lsdtt_channel_head_metric_record> one_a;
    one_a.push_back(make_summary_head(1, 1, "1|2", 100.0, 0.1, 10.0, 1.0));
    one_a.push_back(make_summary_head(2, 2, "1|2", 400.0, 0.5,  5.0, 10.0));
    one_a.push_back(make_summary_head(3, 2, "1|2", 700.0, 0.7, 15.0, 20.0));

    const lsdtt::lsdtt_pair_uncertainty_method methods[4] = {
        lsdtt::lsdtt_pair_uncertainty_method::std_dev,
        lsdtt::lsdtt_pair_uncertainty_method::std_err,
        lsdtt::lsdtt_pair_uncertainty_method::bootstrap,
        lsdtt::lsdtt_pair_uncertainty_method::ttest};

    bool all_nan = true;
    bool mean_ok = true;
    for (const lsdtt::lsdtt_pair_uncertainty_method m : methods)
    {
        options.uncertainty_method = m;
        const auto rows =
            lsdtt::compute_pair_summaries(segments_for(one_a), one_a, options);
        if (rows.size() != 1) { all_nan = false; break; }
        all_nan = all_nan &&
                  std::isnan(rows[0].uncertainty_head_elevation) &&
                  std::isnan(rows[0].uncertainty_gradient) &&
                  std::isnan(rows[0].uncertainty_relief) &&
                  std::isnan(rows[0].uncertainty_chi);
        mean_ok = mean_ok &&
                  close_to(rows[0].mean_head_elevation_a, 100.0, 1e-12) &&
                  close_to(rows[0].delta_head_elevation, -450.0, 1e-12);
    }
    check("P5 a single-head side gives NaN uncertainty under every method, "
          "never a meaningless zero", all_nan);
    check("P5 a single-head side still reports its mean and delta", mean_ok);

    // No heads at all on side A. The row must still be emitted -- the divide
    // exists in the geometry output, and a row that says one side is empty is
    // more use than a divide that silently vanishes from the table.
    std::vector<lsdtt::lsdtt_channel_head_metric_record> empty_a;
    empty_a.push_back(make_summary_head(1, 2, "1|2", 400.0, 0.5,  5.0, 10.0));
    empty_a.push_back(make_summary_head(2, 2, "1|2", 700.0, 0.7, 15.0, 20.0));

    options.uncertainty_method = lsdtt::lsdtt_pair_uncertainty_method::std_err;
    const auto empty_rows =
        lsdtt::compute_pair_summaries(segments_for(empty_a), empty_a, options);
    check("P5 a pair with an empty side still emits a row",
          empty_rows.size() == 1, "n=" + std::to_string(empty_rows.size()));
    if (empty_rows.size() != 1) return;

    const lsdtt::lsdtt_pair_summary_record& e = empty_rows[0];
    // unit_a comes from the pair_id, so it is known even with no head on it.
    check("P5 the empty side's unit id is still recovered from the pair_id",
          e.unit_a == 1 && e.unit_b == 2 && e.n_heads_a == 0 && e.n_heads_b == 2,
          std::to_string(e.unit_a) + " n=" + std::to_string(e.n_heads_a));
    check("P5 an empty side gives NaN mean, delta and uncertainty",
          std::isnan(e.mean_head_elevation_a) &&
          std::isnan(e.delta_head_elevation) &&
          std::isnan(e.uncertainty_head_elevation) &&
          std::isnan(e.mean_chi_a) && std::isnan(e.delta_chi));
    check("P5 the populated side is still summarised",
          close_to(e.mean_head_elevation_b, 550.0, 1e-12),
          num(e.mean_head_elevation_b));

    // An input with no matched heads at all must give no rows, not a crash.
    const auto no_rows = lsdtt::compute_pair_summaries({}, {}, options);
    check("P5 an empty head list gives no rows and does not crash",
          no_rows.empty(), "n=" + std::to_string(no_rows.size()));
}

// =============================================================================
// P6 -- bootstrap determinism
// =============================================================================
void test_pair_summary_bootstrap_determinism()
{
    // A five-head fixture is too coarse for the seed sensitivity check below:
    // with 3 and 2 values a side the delta distribution has only a handful of
    // attainable values and the 2.5 / 97.5 percentiles hit the extremes under
    // almost any seed. This fixture has 14 and 11 irregular values a side, and
    // each metric gets a different spread so the per-metric streams are
    // distinguishable.
    std::vector<lsdtt::lsdtt_channel_head_metric_record> heads;
    for (int i = 0; i < 14; ++i)
    {
        const double t = static_cast<double>((i * 7) % 13);
        heads.push_back(make_summary_head(
            i + 1, 1, "1|2", 100.0 + 11.0 * t, 0.1 + 0.03 * t,
            10.0 + 4.0 * t, 1.0 + 0.7 * t));
    }
    for (int i = 0; i < 11; ++i)
    {
        const double t = static_cast<double>((i * 5) % 11);
        heads.push_back(make_summary_head(
            i + 101, 2, "1|2", 400.0 + 9.0 * t, 0.5 + 0.05 * t,
            5.0 + 6.0 * t, 10.0 + 1.3 * t));
    }

    lsdtt::lsdtt_pair_summary_options options;
    options.uncertainty_method = lsdtt::lsdtt_pair_uncertainty_method::bootstrap;
    options.bootstrap_replicates = 500;
    options.bootstrap_seed = 5489;

    const auto first  =
        lsdtt::compute_pair_summaries(segments_for(heads), heads, options);
    const auto second =
        lsdtt::compute_pair_summaries(segments_for(heads), heads, options);
    check("P6 the same seed twice gives bit-identical uncertainties",
          first.size() == 1 && second.size() == 1 &&
          first[0].uncertainty_head_elevation ==
              second[0].uncertainty_head_elevation &&
          first[0].uncertainty_gradient == second[0].uncertainty_gradient &&
          first[0].uncertainty_relief == second[0].uncertainty_relief &&
          first[0].uncertainty_chi == second[0].uncertainty_chi,
          first.empty() ? "no rows" : num(first[0].uncertainty_head_elevation));

    options.bootstrap_seed = 12345;
    const auto other_seed =
        lsdtt::compute_pair_summaries(segments_for(heads), heads, options);
    check("P6 a different seed gives a different answer",
          other_seed.size() == 1 && first.size() == 1 &&
          other_seed[0].uncertainty_head_elevation !=
              first[0].uncertainty_head_elevation,
          other_seed.empty() ? "no rows"
                             : num(other_seed[0].uncertainty_head_elevation));

    // Each metric must draw its own stream, or all four columns of a pair whose
    // metrics happen to share a scale would move together.
    options.bootstrap_seed = 5489;
    check("P6 each metric draws from its own stream, not one shared one",
          first.size() == 1 &&
          first[0].uncertainty_head_elevation != first[0].uncertainty_relief);

    // Seeding must not depend on iteration position. Prepending an unrelated
    // pair shifts 1|2 to the second row; its bootstrap value must not move.
    std::vector<lsdtt::lsdtt_channel_head_metric_record> with_extra;
    with_extra.push_back(make_summary_head(90, 3, "3|4", 11.0, 0.9, 90.0, 30.0));
    with_extra.push_back(make_summary_head(91, 3, "3|4", 13.0, 0.8, 80.0, 20.0));
    with_extra.push_back(make_summary_head(92, 4, "3|4", 21.0, 0.4, 40.0, 10.0));
    with_extra.push_back(make_summary_head(93, 4, "3|4", 27.0, 0.2, 20.0,  5.0));
    with_extra.insert(with_extra.end(), heads.begin(), heads.end());

    const auto extra_rows =
        lsdtt::compute_pair_summaries(
            segments_for(with_extra), with_extra, options);
    const lsdtt::lsdtt_pair_summary_record* one_two = find_pair(extra_rows, "1|2");
    check("P6 an unrelated extra pair does not perturb 1|2's bootstrap",
          extra_rows.size() == 2 && one_two != nullptr && first.size() == 1 &&
          one_two->uncertainty_head_elevation ==
              first[0].uncertainty_head_elevation,
          one_two == nullptr ? "1|2 missing"
                             : num(one_two->uncertainty_head_elevation));
    check("P6 rows come out ordered by pair_id",
          extra_rows.size() == 2 && extra_rows[0].pair_id == "1|2" &&
          extra_rows[1].pair_id == "3|4",
          extra_rows.size() == 2 ? extra_rows[0].pair_id + "," +
                                       extra_rows[1].pair_id
                                 : "n=" + std::to_string(extra_rows.size()));
}

// =============================================================================
// P7 -- unmatched heads take no part
// =============================================================================
void test_pair_summary_excludes_unmatched_heads()
{
    lsdtt::lsdtt_pair_summary_options options;

    std::vector<lsdtt::lsdtt_channel_head_metric_record> heads =
        make_summary_fixture(1, 2);
    const auto reference =
        lsdtt::compute_pair_summaries(segments_for(heads), heads, options);

    // Heads whose unit borders no extracted divide carry an EMPTY
    // pair_id_nearest. Their values here are wild enough that including them in
    // either side would be obvious.
    heads.push_back(make_summary_head(10, 1, "", 1.0e6, 99.0, 9999.0, 999.0));
    heads.push_back(make_summary_head(11, 2, "", 1.0e6, 99.0, 9999.0, 999.0));
    const auto rows =
        lsdtt::compute_pair_summaries(segments_for(heads), heads, options);

    check("P7 unmatched heads create no row of their own",
          rows.size() == 1, "n=" + std::to_string(rows.size()));
    if (rows.size() != 1 || reference.size() != 1) return;

    check("P7 unmatched heads are not counted on either side",
          rows[0].n_heads_a == 3 && rows[0].n_heads_b == 2,
          std::to_string(rows[0].n_heads_a) + "," +
          std::to_string(rows[0].n_heads_b));
    check("P7 every mean and delta is unchanged by their presence",
          rows[0].mean_head_elevation_a == reference[0].mean_head_elevation_a &&
          rows[0].mean_head_elevation_b == reference[0].mean_head_elevation_b &&
          rows[0].delta_relief == reference[0].delta_relief &&
          rows[0].delta_chi == reference[0].delta_chi,
          num(rows[0].mean_head_elevation_a));
}

// =============================================================================
// P8 -- a head is never counted on the wrong side of its pair
// =============================================================================
void test_pair_summary_side_assignment()
{
    lsdtt::lsdtt_pair_summary_options options;

    // Three units in a row, so unit 2 borders BOTH pairs. Its heads must split
    // by their own pair_id_nearest: the 1|2 head onto side B of 1|2, and the
    // 2|3 head onto side A of 2|3. Putting either on the other side, or on the
    // other pair, changes a mean.
    std::vector<lsdtt::lsdtt_channel_head_metric_record> heads;
    heads.push_back(make_summary_head(1, 1, "1|2", 100.0, 0.1, 10.0, 1.0));
    heads.push_back(make_summary_head(2, 1, "1|2", 200.0, 0.2, 20.0, 2.0));
    heads.push_back(make_summary_head(3, 2, "1|2", 500.0, 0.5, 50.0, 5.0));
    heads.push_back(make_summary_head(4, 2, "1|2", 700.0, 0.7, 70.0, 7.0));
    heads.push_back(make_summary_head(5, 2, "2|3", 900.0, 0.9, 90.0, 9.0));
    heads.push_back(make_summary_head(6, 2, "2|3", 1100.0, 1.1, 110.0, 11.0));
    heads.push_back(make_summary_head(7, 3, "2|3", 1300.0, 1.3, 130.0, 13.0));
    heads.push_back(make_summary_head(8, 3, "2|3", 1500.0, 1.5, 150.0, 15.0));

    // A head tagged with a pair its own unit does not belong to cannot come out
    // of compute_nearest_divide_metrics, which only matches a head to divides
    // bounding its own unit. It must not be assigned to a side by elimination.
    heads.push_back(make_summary_head(9, 9, "1|2", 1.0e6, 99.0, 9999.0, 999.0));

    const auto rows =
        lsdtt::compute_pair_summaries(segments_for(heads), heads, options);
    check("P8 two pairs give two rows", rows.size() == 2,
          "n=" + std::to_string(rows.size()));

    const lsdtt::lsdtt_pair_summary_record* p12 = find_pair(rows, "1|2");
    const lsdtt::lsdtt_pair_summary_record* p23 = find_pair(rows, "2|3");
    check("P8 both pairs are present", p12 != nullptr && p23 != nullptr);
    if (p12 == nullptr || p23 == nullptr) return;

    check("P8 a unit-2 head lands on side B of pair 1|2",
          p12->n_heads_a == 2 && p12->n_heads_b == 2 &&
          close_to(p12->mean_head_elevation_a, 150.0, 1e-12) &&
          close_to(p12->mean_head_elevation_b, 600.0, 1e-12),
          num(p12->mean_head_elevation_a) + "," +
          num(p12->mean_head_elevation_b));
    check("P8 a unit-2 head lands on side A of pair 2|3",
          p23->n_heads_a == 2 && p23->n_heads_b == 2 &&
          close_to(p23->mean_head_elevation_a, 1000.0, 1e-12) &&
          close_to(p23->mean_head_elevation_b, 1400.0, 1e-12),
          num(p23->mean_head_elevation_a) + "," +
          num(p23->mean_head_elevation_b));
    check("P8 unit 2's heads are split by pair, never pooled across pairs",
          p12->n_heads_b + p23->n_heads_a == 4 &&
          !close_to(p12->mean_head_elevation_b,
                    p23->mean_head_elevation_a, 1e-9));
    check("P8 a head whose unit is on neither side of its pair is dropped",
          p12->n_heads_a == 2 && p12->n_heads_b == 2 &&
          close_to(p12->mean_head_elevation_b, 600.0, 1e-12),
          std::to_string(p12->n_heads_a + p12->n_heads_b));
}

// =============================================================================
// P9 -- the pair set is taken from the DIVIDES, not from the heads
// =============================================================================
void test_pair_summary_pairs_come_from_divides()
{
    lsdtt::lsdtt_pair_summary_options options;
    options.uncertainty_method = lsdtt::lsdtt_pair_uncertainty_method::std_err;

    // Two extracted divides, but every head sits on 1|2. The 3|4 divide is in
    // the geometry output, so it must be in the table too -- summarised from
    // nothing, which is the empty-side convention taken one step further.
    std::vector<lsdtt::lsdtt_divide_segment_record> segments;
    segments.push_back(make_segment("1|2", 1, 2,
        {lsdtt::lsdtt_xy{0.0, 0.0}, lsdtt::lsdtt_xy{100.0, 0.0}}));
    segments.push_back(make_segment("3|4", 3, 4,
        {lsdtt::lsdtt_xy{0.0, 500.0}, lsdtt::lsdtt_xy{100.0, 500.0}}));

    const std::vector<lsdtt::lsdtt_channel_head_metric_record> heads =
        make_summary_fixture(1, 2);

    const auto rows = lsdtt::compute_pair_summaries(segments, heads, options);
    check("P9 a divide with no matched heads still gets a row",
          rows.size() == 2, "n=" + std::to_string(rows.size()));

    const lsdtt::lsdtt_pair_summary_record* empty_pair = find_pair(rows, "3|4");
    check("P9 the headless pair carries its units from the divide",
          empty_pair != nullptr && empty_pair->unit_a == 3 &&
          empty_pair->unit_b == 4,
          empty_pair == nullptr ? "3|4 missing"
                                : std::to_string(empty_pair->unit_a));
    if (empty_pair == nullptr) return;

    check("P9 a headless pair reports zero heads on both sides",
          empty_pair->n_heads_a == 0 && empty_pair->n_heads_b == 0 &&
          empty_pair->n_head_elevation_a == 0 &&
          empty_pair->n_head_elevation_b == 0 && empty_pair->n_relief_a == 0 &&
          empty_pair->n_chi_b == 0,
          std::to_string(empty_pair->n_heads_a) + "," +
          std::to_string(empty_pair->n_heads_b));
    check("P9 every statistic of a headless pair is NaN, never zero",
          std::isnan(empty_pair->mean_head_elevation_a) &&
          std::isnan(empty_pair->mean_head_elevation_b) &&
          std::isnan(empty_pair->delta_head_elevation) &&
          std::isnan(empty_pair->mean_gradient_a) &&
          std::isnan(empty_pair->delta_gradient) &&
          std::isnan(empty_pair->mean_relief_b) &&
          std::isnan(empty_pair->delta_relief) &&
          std::isnan(empty_pair->mean_chi_a) &&
          std::isnan(empty_pair->delta_chi) &&
          std::isnan(empty_pair->uncertainty_head_elevation) &&
          std::isnan(empty_pair->uncertainty_gradient) &&
          std::isnan(empty_pair->uncertainty_relief) &&
          std::isnan(empty_pair->uncertainty_chi),
          num(empty_pair->delta_head_elevation));
    check("P9 the headless row still records the method used",
          empty_pair->uncertainty_method == "std_err",
          empty_pair->uncertainty_method);

    // The populated pair must be untouched by the presence of the empty one.
    const lsdtt::lsdtt_pair_summary_record* full_pair = find_pair(rows, "1|2");
    check("P9 the populated pair is summarised exactly as before",
          full_pair != nullptr && full_pair->n_heads_a == 3 &&
          full_pair->n_heads_b == 2 &&
          close_to(full_pair->mean_head_elevation_a, 200.0, 1e-12) &&
          close_to(full_pair->delta_head_elevation, -350.0, 1e-12),
          full_pair == nullptr ? "1|2 missing"
                               : num(full_pair->delta_head_elevation));

    // Two segments sharing a pair_id are the same unit pair meeting along two
    // separate stretches of divide. That is ONE divide pair and must give ONE
    // row, with every head still counted once.
    std::vector<lsdtt::lsdtt_divide_segment_record> two_stretches;
    two_stretches.push_back(make_segment("1|2", 1, 2,
        {lsdtt::lsdtt_xy{0.0, 0.0}, lsdtt::lsdtt_xy{100.0, 0.0}}));
    two_stretches.push_back(make_segment("1|2", 1, 2,
        {lsdtt::lsdtt_xy{900.0, 0.0}, lsdtt::lsdtt_xy{1000.0, 0.0}}));

    const auto two_rows =
        lsdtt::compute_pair_summaries(two_stretches, heads, options);
    check("P9 two segments sharing a pair_id give one row, not two",
          two_rows.size() == 1, "n=" + std::to_string(two_rows.size()));
    check("P9 the duplicated pair counts each head exactly once",
          two_rows.size() == 1 && two_rows[0].n_heads_a == 3 &&
          two_rows[0].n_heads_b == 2 &&
          close_to(two_rows[0].mean_head_elevation_a, 200.0, 1e-12),
          two_rows.empty() ? "no rows"
                           : std::to_string(two_rows[0].n_heads_a));

    // A head naming a pair that no segment carries must not invent a row. This
    // cannot come out of compute_nearest_divide_metrics -- every
    // pair_id_nearest is read off one of the segments handed in here -- so the
    // invariant is asserted rather than quietly accommodated.
    std::vector<lsdtt::lsdtt_channel_head_metric_record> stray = heads;
    stray.push_back(make_summary_head(20, 5, "5|6", 1.0e6, 99.0, 9999.0, 999.0));
    stray.push_back(make_summary_head(21, 6, "5|6", 1.0e6, 99.0, 9999.0, 999.0));

    std::vector<lsdtt::lsdtt_divide_segment_record> only_one_two;
    only_one_two.push_back(make_segment("1|2", 1, 2,
        {lsdtt::lsdtt_xy{0.0, 0.0}, lsdtt::lsdtt_xy{100.0, 0.0}}));

    const auto stray_rows =
        lsdtt::compute_pair_summaries(only_one_two, stray, options);
    check("P9 a head naming a pair absent from the divides invents no row",
          stray_rows.size() == 1 && stray_rows[0].pair_id == "1|2",
          stray_rows.size() == 1 ? stray_rows[0].pair_id
                                 : "n=" + std::to_string(stray_rows.size()));
    check("P9 that head is not folded into any other pair either",
          stray_rows.size() == 1 && stray_rows[0].n_heads_a == 3 &&
          stray_rows[0].n_heads_b == 2 &&
          close_to(stray_rows[0].mean_head_elevation_a, 200.0, 1e-12),
          stray_rows.empty() ? "no rows"
                             : num(stray_rows[0].mean_head_elevation_a));

    // Divides but no heads at all: still one row per divide pair.
    const auto headless = lsdtt::compute_pair_summaries(segments, {}, options);
    check("P9 divides with no heads at all still give one row each",
          headless.size() == 2 && headless[0].pair_id == "1|2" &&
          headless[1].pair_id == "3|4",
          "n=" + std::to_string(headless.size()));
}

// =============================================================================
// CLI fixture: a deterministic three-valley DEM plus its outlet layer
//
// Elevation rises to the north and is corrugated across x with troughs at
// x = 2000, 6000 and 10000 m, so the grid holds three parallel south-draining
// valleys separated by two ridges. Units 1 and 3 are therefore NOT adjacent,
// which is exactly the arrangement N6 pins in memory.
//
// The outlet layer carries two cluster fields so that one fixture drives all
// three CLI cases:
//   cluster_id  -> 1, 1, 2  (two units: the middle divide is internal)
//   cluster_all -> 7, 7, 7  (one unit: must be rejected)
// =============================================================================
// The three-valley elevation surface itself, shared by the fixture writer and
// by the adjacency check (A1) so that both work on exactly the same grid.
const int    k_valley_n_rows = 80;
const int    k_valley_n_cols = 120;
const double k_valley_cell   = 100.0;
// One outlet at the mouth of each valley: the bottom row, on the trough.
const int    k_valley_outlet_cols[3] = {19, 59, 99};

lsdtt::lsdtt_raster_cpu make_three_valley_dem()
{
    const int    n_rows = k_valley_n_rows;
    const int    n_cols = k_valley_n_cols;
    const double cell   = k_valley_cell;

    lsdtt::lsdtt_raster_meta meta;
    meta.n_rows = n_rows;
    meta.n_cols = n_cols;
    meta.cell_size = cell;
    meta.no_data_value = -9999.0;
    meta.set_default_synthetic_georeference();

    constexpr double k_pi = 3.14159265358979323846;
    // The troughs sit on cell CENTRES (col 19, 59, 99), not on cell edges. A
    // trough on an edge is a perfect elevation tie between two columns, which
    // splits each valley into two half-basins and destroys the adjacency the
    // divide extraction depends on.
    const double period = 4000.0;     // valley spacing (m)
    const double x_trough = 1950.0;   // easting of the first trough

    std::vector<float> data(static_cast<std::size_t>(n_rows) * n_cols, 0.0f);
    for (int row = 0; row < n_rows; ++row)
    {
        for (int col = 0; col < n_cols; ++col)
        {
            const double x = (col + 0.5) * cell;               // local easting
            const double y = (n_rows - row - 0.5) * cell;      // local northing

            // Cross-valley corrugation: 0 at every trough, 120 m at every crest.
            const double cross =
                60.0 * (1.0 - std::cos(2.0 * k_pi * (x - x_trough) / period));
            // Down-valley fall to the south, plus a hillslope corrugation in y
            // that seeds tributaries (and therefore junctions) without
            // disturbing the trough, where cross is zero.
            const double down = 0.05 * y;
            const double trib = 1.0 + 0.12 * std::cos(2.0 * k_pi * y / 1600.0);

            data[static_cast<std::size_t>(row) * n_cols + col] =
                static_cast<float>(down + cross * trib);
        }
    }

    std::vector<std::uint8_t> mask(data.size(), 1);
    return lsdtt::lsdtt_raster_cpu(std::move(data), std::move(mask), meta);
}

void write_cli_fixture(const std::string& out_dir)
{
    const int    n_rows = k_valley_n_rows;
    const double cell   = k_valley_cell;

    const lsdtt::lsdtt_raster_cpu dem = make_three_valley_dem();
    const lsdtt::lsdtt_raster_meta meta = dem.get_meta();
    const std::vector<float>& data = dem.get_data();

    const double x0 = meta.x_min;
    const double y0 = meta.y_min;

    std::filesystem::create_directories(out_dir);
    const std::string dem_path = out_dir + "/divmig_valleys.tif";
    lsdtt::lsdtt_data_loader writer;
    writer.write_full(data, meta, dem_path,
                      lsdtt::raster_output_format::gtiff);
    std::cout << "Wrote CLI fixture DEM: " << dem_path << "\n";

    const int outlet_row = n_rows - 1;
    const std::vector<int> outlet_cols(std::begin(k_valley_outlet_cols),
                                       std::end(k_valley_outlet_cols));
    const std::vector<std::string> cluster_id  = {"1", "1", "2"};
    const std::vector<std::string> cluster_all = {"7", "7", "7"};

    lsdtt::lsdtt_vector_schema schema;
    schema.field_names = {"point_id", "cluster_id", "cluster_all"};
    schema.field_types = {lsdtt::lsdtt_vector_field_type::int64_type,
                          lsdtt::lsdtt_vector_field_type::string_type,
                          lsdtt::lsdtt_vector_field_type::string_type};

    const std::string outlet_path = out_dir + "/divmig_outlets.fgb";
    lsdtt::lsdtt_data_loader vwriter;
    vwriter.begin_vector_write(outlet_path, "outlets", schema, meta.crs_wkt,
                               lsdtt::vector_output_format::flatgeobuf, false);
    for (std::size_t i = 0; i < outlet_cols.size(); ++i)
    {
        const double x = x0 + (outlet_cols[i] + 0.5) * cell;
        const double y = y0 + (n_rows - outlet_row - 0.5) * cell;
        vwriter.write_point(x, y, {static_cast<long int>(i + 1),
                                   cluster_id[i], cluster_all[i]});
    }
    vwriter.finish_vector_write();
    std::cout << "Wrote CLI fixture outlets: " << outlet_path << "\n";

    // The same three outlets as a WGS84 CSV. Both forms are needed: the CLI
    // reads cluster fields from a CSV by column name and from a vector layer
    // through OGR, and only the CSV path is usable today (see the CLI defect
    // reported alongside these tests).
    OGRSpatialReference utm;
    OGRSpatialReference wgs84;
    if (utm.importFromWkt(meta.crs_wkt.c_str()) != OGRERR_NONE ||
        wgs84.importFromEPSG(4326) != OGRERR_NONE)
    {
        throw std::runtime_error("could not build the fixture CRS transform");
    }
    utm.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);
    wgs84.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);
    OGRCoordinateTransformation* ct =
        OGRCreateCoordinateTransformation(&utm, &wgs84);
    if (ct == nullptr)
    {
        throw std::runtime_error("could not create the fixture CRS transform");
    }

    const std::string csv_path = out_dir + "/divmig_outlets.csv";
    std::ofstream csv(csv_path);
    if (!csv) throw std::runtime_error("could not write " + csv_path);
    csv << "latitude,longitude,point_id,cluster_id,cluster_all\n";
    csv << std::setprecision(12);
    for (std::size_t i = 0; i < outlet_cols.size(); ++i)
    {
        double x = x0 + (outlet_cols[i] + 0.5) * cell;
        double y = y0 + (n_rows - outlet_row - 0.5) * cell;
        if (!ct->Transform(1, &x, &y))
        {
            throw std::runtime_error("fixture outlet transform failed");
        }
        csv << y << "," << x << "," << (i + 1) << ","
            << cluster_id[i] << "," << cluster_all[i] << "\n";
    }
    csv.close();
    OCTDestroyCoordinateTransformation(ct);
    std::cout << "Wrote CLI fixture outlet CSV: " << csv_path << "\n";
}

// =============================================================================
// A1 -- adjacency-driven pair generation finds exactly the pairs the old
//       all-vs-all loop did, on a grid that contains a NON-adjacent pair
//
// The three-valley surface is the right fixture for this: the outer valleys
// (basins 0 and 2) are separated by the middle one, so pair (0, 2) is the case
// the old loop only discovered was hopeless after paying for a full ridgeline
// extraction. The check runs both strategies over the same basins and compares
// the divides they produce, so the point is not merely that the adjacency pass
// is fast but that it loses nothing.
// =============================================================================
void test_adjacency_matches_all_vs_all()
{
    const lsdtt::lsdtt_raster_cpu dem = make_three_valley_dem();
    const lsdtt::lsdtt_raster_cpu filled = lsdtt::fill_pits(dem, 0.0001f);

    lsdtt::lsdtt_flowinfo_options fi_opts;
    fi_opts.boundary_conditions       = "oooo";
    fi_opts.interior_nodata_is_outlet = true;
    fi_opts.require_acyclic_graph     = true;
    fi_opts.build_topological_order   = true;
    fi_opts.build_s_vector            = true;

    std::unique_ptr<lsdtt::lsdtt_flow_base> flow = lsdtt::make_flow_base(
        filled.get_data_array(), filled.get_mask_array(), filled.get_meta(),
        fi_opts);

    const int n_rows = filled.get_meta().n_rows;
    const int n_cols = filled.get_meta().n_cols;

    std::vector<std::int64_t> outlet_flats;
    for (const int col : k_valley_outlet_cols)
    {
        outlet_flats.push_back(
            static_cast<std::int64_t>(n_rows - 1) * n_cols + col);
    }

    lsdtt::lsdtt_junction_network_options jn_opts;
    jn_opts.build_topological_order = true;
    jn_opts.build_s_vector          = true;
    lsdtt::lsdtt_junction_network_cpu junction_network;
    junction_network.create(*flow, flow->get_sources_index_threshold(1000),
                            jn_opts);

    lsdtt::lsdtt_junction_basin_manager_cpu basins(*flow, junction_network);
    basins.select_junction_basins_from_basin_outlets(outlet_flats);
    basins.filter_duplicates();

    check("A1 the fixture gives the three basins the check needs",
          basins.size() == 3, "n=" + std::to_string(basins.size()));
    if (basins.size() != 3) return;

    // One extracted divide, reduced to the things the outputs are built from.
    struct extracted
    {
        int i = -1;
        int j = -1;
        std::size_t n_points = 0;
        double length_m = 0.0;
    };

    auto extract = [&](std::size_t i, std::size_t j) -> extracted {
        extracted out;
        const auto& basin_i = basins.basin_at(i);
        const auto& basin_j = basins.basin_at(j);
        // The CLI's nesting rule, kept here so the two strategies are compared
        // under identical exclusions.
        if (flow->is_upslope_of(basin_i.get_basin_outlet_flat_index(),
                                basin_j.get_basin_outlet_flat_index()) ||
            flow->is_upslope_of(basin_j.get_basin_outlet_flat_index(),
                                basin_i.get_basin_outlet_flat_index()))
        {
            return out;
        }
        lsdtt::lsdtt_ridgeline_cpu ridge;
        try
        {
            ridge = basin_i.get_basin_outline().get_interfluvial_divide_as_ridgeline(
                basin_i.get_junction_number(), basin_j.get_junction_number(),
                junction_network);
        }
        catch (const std::exception&)
        {
            return out;
        }
        if (ridge.empty()) return out;
        out.i = static_cast<int>(i);
        out.j = static_cast<int>(j);
        out.n_points = ridge.positions().size();
        out.length_m = ridge.total_length();
        return out;
    };

    // The old strategy: every pair, non-adjacency discovered after the fact.
    std::vector<extracted> all_vs_all;
    for (std::size_t i = 0; i < basins.size(); ++i)
    {
        for (std::size_t j = i + 1; j < basins.size(); ++j)
        {
            const extracted e = extract(i, j);
            if (e.i >= 0) all_vs_all.push_back(e);
        }
    }

    // The new strategy: only pairs whose cells touch on the grid.
    std::vector<lsdtt::lsdtt_flat_index_slice> basin_cells;
    for (std::size_t i = 0; i < basins.size(); ++i)
    {
        basin_cells.push_back(basins.basin_at(i).get_slice());
    }
    const std::vector<std::pair<int, int>> adjacent =
        lsdtt::find_adjacent_basin_pairs(basin_cells, n_rows, n_cols);

    std::string adjacent_text;
    for (const auto& p : adjacent)
    {
        adjacent_text += std::to_string(p.first) + "-" +
                         std::to_string(p.second) + " ";
    }
    check("A1 the adjacency pass finds the two touching pairs and not the "
          "outer one",
          adjacent.size() == 2 && adjacent[0] == std::make_pair(0, 1) &&
          adjacent[1] == std::make_pair(1, 2),
          adjacent_text);

    // The pair the adjacency pass refuses to look at must be one the old loop
    // got nothing from, or the saving would be a loss.
    const extracted outer = extract(0, 2);
    check("A1 the pair the adjacency pass skips yields no divide anyway",
          outer.i < 0, "n_points=" + std::to_string(outer.n_points));

    std::vector<extracted> adjacency_driven;
    for (const auto& p : adjacent)
    {
        const extracted e = extract(static_cast<std::size_t>(p.first),
                                    static_cast<std::size_t>(p.second));
        if (e.i >= 0) adjacency_driven.push_back(e);
    }

    check("A1 both strategies produce the same number of divides",
          adjacency_driven.size() == all_vs_all.size(),
          std::to_string(adjacency_driven.size()) + " vs " +
          std::to_string(all_vs_all.size()));
    if (adjacency_driven.size() != all_vs_all.size()) return;

    bool same = true;
    for (std::size_t k = 0; k < all_vs_all.size(); ++k)
    {
        same = same && adjacency_driven[k].i == all_vs_all[k].i &&
               adjacency_driven[k].j == all_vs_all[k].j &&
               adjacency_driven[k].n_points == all_vs_all[k].n_points &&
               adjacency_driven[k].length_m == all_vs_all[k].length_m;
    }
    check("A1 the divides come out in the same order with the same geometry",
          same && all_vs_all.size() == 2);
}

// =============================================================================
// CLI fixture 2: a main ridge with one valley on one side and two on the other
//
// The three-valley fixture above cannot exercise fragment merging: every unit
// pair there is adjacent along exactly one stretch of divide, so there is
// nothing to join. This one is the smallest domain that does.
//
// A single E-W ridge crosses the middle of the grid. The northern half holds one
// valley draining north; the southern half holds two valleys draining south,
// with their troughs half a period out of phase with the northern one, so the
// northern basin faces BOTH southern basins across the ridge. Grouping the two
// southern outlets into one unit therefore leaves that unit bordering the
// northern unit along two fragments -- one either side of the point where the
// south-south divide meets the main ridge -- and those two fragments are
// contiguous at that triple point. The run must report ONE divide segment.
//
//   dy    = |y - y_ridge|
//   ramp  = 0.15 + 0.85 min(1, dy / 4000)
//   trib  = 1 + 0.12 cos(2 pi dy / 1600)
//   cross = 90 (1 - cos(2 pi (x - x_trough) / 6000)) ramp trib
//   z     = 400 - 0.06 dy + cross
//
// Design constraints behind those numbers:
//   - troughs and the ridge line sit on cell CENTRES (x = 2950, 5950, 8950;
//     y = 5950). A trough or crest on a cell edge is a perfect elevation tie
//     that splits a valley in two and destroys the adjacency the divide
//     extraction depends on -- the same trap noted on the fixture above,
//   - the ramp never reaches zero, so the ridge row is not a 120-cell flat,
//   - the ramp gains 180 * 0.85/4000 = 0.038 m/m on the crests, less than the
//     0.06 m/m fall away from the ridge, so crest lines still drain outwards,
//   - trib modulates ridge height along the valley, seeding the tributaries
//     that give each valley a junction to hang a basin on.
// =============================================================================
void write_merge_cli_fixture(const std::string& out_dir)
{
    const int    n_rows = 120;
    const int    n_cols = 120;
    const double cell   = 100.0;

    lsdtt::lsdtt_raster_meta meta;
    meta.n_rows = n_rows;
    meta.n_cols = n_cols;
    meta.cell_size = cell;
    meta.no_data_value = -9999.0;
    meta.set_default_synthetic_georeference();

    constexpr double k_pi = 3.14159265358979323846;
    const double period      = 6000.0;
    const double y_ridge     = 5950.0;   // row 60, a cell centre
    const double x_trough_n  = 5950.0;   // col 59: the single northern valley
    const double x_trough_s  = 2950.0;   // col 29: the western southern valley
                                         // (col 89 is its half-period partner)

    std::vector<float> data(static_cast<std::size_t>(n_rows) * n_cols, 0.0f);
    for (int row = 0; row < n_rows; ++row)
    {
        for (int col = 0; col < n_cols; ++col)
        {
            const double x = (col + 0.5) * cell;
            const double y = (n_rows - row - 0.5) * cell;
            const double dy = std::fabs(y - y_ridge);

            const double x_trough = (y > y_ridge) ? x_trough_n : x_trough_s;
            const double ramp = 0.15 + 0.85 * std::min(1.0, dy / 4000.0);
            const double trib = 1.0 + 0.12 * std::cos(2.0 * k_pi * dy / 1600.0);
            const double cross =
                90.0 * (1.0 - std::cos(2.0 * k_pi * (x - x_trough) / period)) *
                ramp * trib;

            data[static_cast<std::size_t>(row) * n_cols + col] =
                static_cast<float>(400.0 - 0.06 * dy + cross);
        }
    }

    std::filesystem::create_directories(out_dir);
    const std::string dem_path = out_dir + "/divmig_ridge.tif";
    lsdtt::lsdtt_data_loader writer;
    writer.write_full(data, meta, dem_path,
                      lsdtt::raster_output_format::gtiff);
    std::cout << "Wrote merge fixture DEM: " << dem_path << "\n";

    // One outlet at the mouth of each of the three valleys: the northern one on
    // the top row, the two southern ones on the bottom row. cluster_id puts the
    // two southern valleys in one unit, so that unit meets the northern unit
    // along two contiguous fragments.
    struct outlet_spec { int row; int col; const char* cluster; };
    const std::vector<outlet_spec> outlets = {
        {0,           59, "1"},   // north valley
        {n_rows - 1,  29, "2"},   // south-west valley
        {n_rows - 1,  89, "2"}    // south-east valley
    };

    lsdtt::lsdtt_vector_schema schema;
    schema.field_names = {"point_id", "cluster_id"};
    schema.field_types = {lsdtt::lsdtt_vector_field_type::int64_type,
                          lsdtt::lsdtt_vector_field_type::string_type};

    const std::string outlet_path = out_dir + "/divmig_ridge_outlets.fgb";
    lsdtt::lsdtt_data_loader vwriter;
    vwriter.begin_vector_write(outlet_path, "outlets", schema, meta.crs_wkt,
                               lsdtt::vector_output_format::flatgeobuf, false);
    for (std::size_t i = 0; i < outlets.size(); ++i)
    {
        const double x = meta.x_min + (outlets[i].col + 0.5) * cell;
        const double y = meta.y_min + (n_rows - outlets[i].row - 0.5) * cell;
        vwriter.write_point(x, y, {static_cast<long int>(i + 1),
                                   std::string(outlets[i].cluster)});
    }
    vwriter.finish_vector_write();
    std::cout << "Wrote merge fixture outlets: " << outlet_path << "\n";
}

} // namespace

int main(int argc, char* argv[])
{
    std::cout << "=== lsdtt3 divide-migration unit tests ===\n";

    test_nearest_point_geometry();
    test_unit_without_divide();
    test_nearest_point_on_nodata();
    test_three_unit_isolation();
    test_common_base_level();
    test_merge_two_contiguous_fragments();
    test_merge_leaves_disjoint_fragments_alone();
    test_merge_reversed_and_out_of_order();
    test_merge_closed_ring();
    test_merge_never_crosses_pair_id();
    test_merge_tolerance_is_one_cell_diagonal();
    test_pair_summary_means_and_deltas();
    test_pair_summary_sign_convention();
    test_pair_summary_uncertainty_methods();
    test_pair_summary_skips_nan();
    test_pair_summary_small_samples();
    test_pair_summary_bootstrap_determinism();
    test_pair_summary_excludes_unmatched_heads();
    test_pair_summary_side_assignment();
    test_pair_summary_pairs_come_from_divides();
    test_adjacency_matches_all_vs_all();

    if (argc > 1)
    {
        try
        {
            write_cli_fixture(argv[1]);
            write_merge_cli_fixture(argv[1]);
        }
        catch (const std::exception& e)
        {
            check("F1 CLI fixture written", false, e.what());
        }
    }

    std::cout << "\nChecks: " << g_n_checks
              << "   failures: " << g_n_failed << "\n";
    std::cout << (g_n_failed == 0 ? "RESULT: PASS\n" : "RESULT: FAIL\n");
    return g_n_failed == 0 ? 0 : 1;
}
