// test_hilltop_flow_routing.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
//
// -----------------------------------------------------------------------------
// Phase 1 of the hillslope-channel coupling work: ridge extraction.
//
// The whole layer rests on one indexing claim -- that divide node (dr, dc) is
// the corner shared by DEM pixels (dr-1,dc-1), (dr-1,dc), (dr,dc-1), (dr,dc),
// so that a cardinal EDGE between two adjacent corners separates exactly two
// DEM pixels.  That claim is verified against
// lsdtt_basin_outline_cpu.cpp:608-624, but a reading of the source is not a
// test, so the fixtures below are drawn by hand and the expected pixels written
// out explicitly.
//
// FIXTURE GEOMETRY (5x5 DEM, non-periodic, so a 6x6 divide grid)
//
//   DEM pixel (r,c) occupies the square whose NORTH-WEST corner is divide
//   node (r,c).  Corner indices therefore run 0..5 and pixel indices 0..4:
//
//        c=0     c=1     c=2     c=3     c=4
//     0-------1-------2-------3-------4-------5     <- corner row 0
//     |       |       |       |       |       |
//     | (0,0) | (0,1) | (0,2) | (0,3) | (0,4) |
//     |       |       |       |       |       |
//     0-------1-------2-------3-------4-------5     <- corner row 1
//     |       |       |       |       |       |
//     | (1,0) | (1,1) | (1,2) | (1,3) | (1,4) |
//     ...
//
//   A horizontal ridgeline along corner row 2 from corner (2,1) to (2,3)
//   therefore separates DEM row 1 (north) from DEM row 2 (south), over
//   columns 1 and 2.
// -----------------------------------------------------------------------------

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "lsdtt_hilltop_flow_routing_cpu.hpp"
#include "lsdtt_ridgeline_cpu.hpp"

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

lsdtt::lsdtt_raster_meta make_dem_meta(int n_rows, int n_cols, double cell)
{
    lsdtt::lsdtt_raster_meta m;
    m.n_rows        = n_rows;
    m.n_cols        = n_cols;
    m.cell_size     = cell;
    m.no_data_value = -9999.0;
    m.set_default_synthetic_georeference();
    return m;
}

// The divide grid: half a cell shifted, n+1 on a non-periodic axis and n on a
// periodic one (lsdtt_basin_outline_cpu.cpp:164-165).
lsdtt::lsdtt_raster_meta make_divide_meta(const lsdtt::lsdtt_raster_meta& dem,
                                          bool periodic_ns, bool periodic_ew)
{
    lsdtt::lsdtt_raster_meta d = dem;
    d.x_min  = dem.x_min - 0.5 * dem.cell_size;
    d.y_min  = dem.y_min - 0.5 * dem.cell_size;
    d.n_rows = dem.n_rows + (periodic_ns ? 0 : 1);
    d.n_cols = dem.n_cols + (periodic_ew ? 0 : 1);
    return d;
}

lsdtt::lsdtt_ridgeline_cpu make_ridgeline(
    const std::vector<std::pair<int,int>>& rc,
    const lsdtt::lsdtt_raster_meta& divide_meta)
{
    std::vector<lsdtt::lsdtt_pixel_rc> nodes;
    nodes.reserve(rc.size());
    for (const auto& p : rc) nodes.push_back(lsdtt::lsdtt_pixel_rc{p.first, p.second});
    return lsdtt::lsdtt_ridgeline_cpu(std::move(nodes), divide_meta);
}

std::set<std::pair<int,int>> pixels_of(
    const lsdtt::lsdtt_ridge_extraction_result& res)
{
    std::set<std::pair<int,int>> out;
    for (const auto& p : res.pixels) out.insert({p.row, p.col});
    return out;
}

std::set<std::pair<int,int>> pixels_of_side(
    const lsdtt::lsdtt_ridge_extraction_result& res, int side)
{
    std::set<std::pair<int,int>> out;
    for (const auto& p : res.pixels)
        if (p.ridge_side == side) out.insert({p.row, p.col});
    return out;
}

std::string show(const std::set<std::pair<int,int>>& s)
{
    std::string out = "{";
    for (const auto& p : s)
        out += "(" + std::to_string(p.first) + "," + std::to_string(p.second) + ")";
    return out + "}";
}

// =============================================================================
// 1. The indexing claim, drawn by hand
// =============================================================================

void test_horizontal_edge_pixels()
{
    const auto dem = make_dem_meta(5, 5, 10.0);
    const auto dgm = make_divide_meta(dem, false, false);

    // corner row 2, travelling EAST from (2,1) to (2,3): two cardinal edges
    auto rl = make_ridgeline({{2,1},{2,2},{2,3}}, dgm);

    const auto res = lsdtt::extract_hilltop_pixels({rl}, dem);

    const std::set<std::pair<int,int>> want_all =
        {{1,1},{2,1},{1,2},{2,2}};
    const std::set<std::pair<int,int>> want_left  = {{1,1},{1,2}};  // north
    const std::set<std::pair<int,int>> want_right = {{2,1},{2,2}};  // south

    check("horizontal edge: separates DEM row 1 from row 2, cols 1-2",
          pixels_of(res) == want_all,
          "got " + show(pixels_of(res)) + " want " + show(want_all));

    check("horizontal edge, travelling east: LEFT is north",
          pixels_of_side(res, 0) == want_left,
          "got " + show(pixels_of_side(res, 0)));

    check("horizontal edge, travelling east: RIGHT is south",
          pixels_of_side(res, 1) == want_right,
          "got " + show(pixels_of_side(res, 1)));

    check("horizontal edge: two cardinal edges counted",
          res.n_edges_cardinal == 2 && res.n_edges_total == 2,
          std::to_string(res.n_edges_cardinal));
}

void test_vertical_edge_pixels()
{
    const auto dem = make_dem_meta(5, 5, 10.0);
    const auto dgm = make_divide_meta(dem, false, false);

    // corner col 3, travelling SOUTH from (1,3) to (3,3)
    auto rl = make_ridgeline({{1,3},{2,3},{3,3}}, dgm);

    const auto res = lsdtt::extract_hilltop_pixels({rl}, dem);

    const std::set<std::pair<int,int>> want_all =
        {{1,2},{1,3},{2,2},{2,3}};
    const std::set<std::pair<int,int>> want_left  = {{1,3},{2,3}};  // east
    const std::set<std::pair<int,int>> want_right = {{1,2},{2,2}};  // west

    check("vertical edge: separates DEM col 2 from col 3, rows 1-2",
          pixels_of(res) == want_all,
          "got " + show(pixels_of(res)) + " want " + show(want_all));

    check("vertical edge, travelling south: LEFT is east",
          pixels_of_side(res, 0) == want_left,
          "got " + show(pixels_of_side(res, 0)));

    check("vertical edge, travelling south: RIGHT is west",
          pixels_of_side(res, 1) == want_right,
          "got " + show(pixels_of_side(res, 1)));
}

// Reversing a ridgeline must swap the sides and change nothing else.
void test_reversal_swaps_sides()
{
    const auto dem = make_dem_meta(5, 5, 10.0);
    const auto dgm = make_divide_meta(dem, false, false);

    auto fwd = make_ridgeline({{2,1},{2,2},{2,3}}, dgm);
    auto rev = make_ridgeline({{2,3},{2,2},{2,1}}, dgm);

    const auto rf = lsdtt::extract_hilltop_pixels({fwd}, dem);
    const auto rr = lsdtt::extract_hilltop_pixels({rev}, dem);

    check("reversal: same pixel set",
          pixels_of(rf) == pixels_of(rr));

    check("reversal: left and right swap",
          pixels_of_side(rf, 0) == pixels_of_side(rr, 1) &&
          pixels_of_side(rf, 1) == pixels_of_side(rr, 0),
          "fwd left " + show(pixels_of_side(rf, 0)) +
          " rev right " + show(pixels_of_side(rr, 1)));
}

// An L-shaped divide: the side assignment must stay consistent around the
// corner, which is the case a per-node (4-pixel) derivation gets wrong.
void test_corner_side_consistency()
{
    const auto dem = make_dem_meta(5, 5, 10.0);
    const auto dgm = make_divide_meta(dem, false, false);

    // east along corner row 2 to (2,3), then south to (4,3)
    auto rl = make_ridgeline({{2,1},{2,2},{2,3},{3,3},{4,3}}, dgm);

    const auto res = lsdtt::extract_hilltop_pixels({rl}, dem);

    // east legs: left = north = (1,1),(1,2); right = south = (2,1),(2,2)
    // south legs: left = east = (2,3),(3,3); right = west = (2,2),(3,2)
    // (2,2) is emitted first as a RIGHT pixel on the east leg, so dedup keeps
    // that record.
    const std::set<std::pair<int,int>> want_left  = {{1,1},{1,2},{2,3},{3,3}};
    const std::set<std::pair<int,int>> want_right = {{2,1},{2,2},{3,2}};

    check("L-shaped divide: left side is consistent around the corner",
          pixels_of_side(res, 0) == want_left,
          "got " + show(pixels_of_side(res, 0)) +
          " want " + show(want_left));

    check("L-shaped divide: right side is consistent around the corner",
          pixels_of_side(res, 1) == want_right,
          "got " + show(pixels_of_side(res, 1)) +
          " want " + show(want_right));

    check("L-shaped divide: no pixel appears on both sides",
          [&]{
              auto l = pixels_of_side(res, 0);
              for (const auto& p : pixels_of_side(res, 1))
                  if (l.count(p)) return false;
              return true;
          }());
}

// =============================================================================
// 2. Boundaries and periodicity
// =============================================================================

void test_divide_on_dem_boundary_drops_outside_pixel()
{
    const auto dem = make_dem_meta(5, 5, 10.0);
    const auto dgm = make_divide_meta(dem, false, false);

    // corner row 0 is the northern DEM boundary: nothing lies north of it
    auto rl = make_ridgeline({{0,1},{0,2},{0,3}}, dgm);

    const auto res = lsdtt::extract_hilltop_pixels({rl}, dem);

    const std::set<std::pair<int,int>> want = {{0,1},{0,2}};

    check("boundary divide: only the inside pixels are emitted",
          pixels_of(res) == want,
          "got " + show(pixels_of(res)) + " want " + show(want));

    check("boundary divide: the off-grid side is counted, not silently lost",
          res.n_off_grid_dropped == 2,
          "got " + std::to_string(res.n_off_grid_dropped));
}

void test_periodic_ew_wraps_at_the_seam()
{
    const auto dem = make_dem_meta(5, 5, 10.0);
    const auto dgm = make_divide_meta(dem, false, true);   // periodic east-west

    check("periodic EW: divide grid keeps n columns",
          dgm.n_cols == 5 && dgm.n_rows == 6,
          std::to_string(dgm.n_rows) + "x" + std::to_string(dgm.n_cols));

    bool pns = false, pew = false;
    lsdtt::infer_divide_grid_periodicity(dem, dgm, pns, pew);
    check("periodic EW: inferred from the meta pair", !pns && pew);

    // Travelling WEST across the seam: corner (2,0) -> (2,4) is a single step
    // through the seam, and the pixel pair is column 4.
    auto rl = make_ridgeline({{2,0},{2,4}}, dgm);
    const auto res = lsdtt::extract_hilltop_pixels({rl}, dem);

    check("periodic EW: seam step is cardinal, not disjoint",
          res.n_edges_cardinal == 1 && res.n_edges_disjoint == 0,
          "cardinal=" + std::to_string(res.n_edges_cardinal) +
          " disjoint=" + std::to_string(res.n_edges_disjoint));

    const std::set<std::pair<int,int>> want = {{1,4},{2,4}};
    check("periodic EW: seam step separates the wrapped column",
          pixels_of(res) == want,
          "got " + show(pixels_of(res)) + " want " + show(want));
}

void test_inconsistent_metas_throw()
{
    const auto dem = make_dem_meta(5, 5, 10.0);
    lsdtt::lsdtt_raster_meta bad = dem;
    bad.n_rows = 9;                       // neither n nor n+1

    bool threw = false;
    try
    {
        bool a = false, b = false;
        lsdtt::infer_divide_grid_periodicity(dem, bad, a, b);
    }
    catch (const std::invalid_argument&) { threw = true; }

    check("meta mismatch: infer_divide_grid_periodicity throws", threw);
}

// =============================================================================
// 3. Step classification and the diagonal policy
// =============================================================================

void test_step_classification()
{
    const auto dem = make_dem_meta(5, 5, 10.0);
    const auto dgm = make_divide_meta(dem, false, false);

    using K = lsdtt::ridgeline_step_kind;
    auto cls = [&](int r1,int c1,int r2,int c2) {
        return lsdtt::classify_ridgeline_step(
            {r1,c1}, {r2,c2}, dgm, false, false);
    };

    check("classify: east step is cardinal",  cls(2,1,2,2) == K::cardinal);
    check("classify: south step is cardinal", cls(1,3,2,3) == K::cardinal);
    check("classify: diagonal step",          cls(2,2,3,3) == K::diagonal);
    check("classify: repeated node",          cls(2,2,2,2) == K::repeated);
    check("classify: two-cell jump is disjoint", cls(2,2,2,4) == K::disjoint);
}

void test_diagonal_policies()
{
    const auto dem = make_dem_meta(5, 5, 10.0);
    const auto dgm = make_divide_meta(dem, false, false);

    auto rl = make_ridgeline({{2,1},{3,2}}, dgm);   // one diagonal step

    // skip (default)
    {
        const auto res = lsdtt::extract_hilltop_pixels({rl}, dem);
        check("diagonal skip: no pixels emitted",
              res.pixels.empty(), std::to_string(res.pixels.size()));
        check("diagonal skip: counted",
              res.n_edges_diagonal == 1,
              std::to_string(res.n_edges_diagonal));
    }

    // bridge
    {
        lsdtt::lsdtt_ridge_extraction_options o;
        o.diagonal_policy = lsdtt::diagonal_step_policy::bridge;
        const auto res = lsdtt::extract_hilltop_pixels({rl}, dem, o);

        check("diagonal bridge: emits pixels via the intervening corner",
              !res.pixels.empty(), std::to_string(res.pixels.size()));
        check("diagonal bridge: still counted as a diagonal step",
              res.n_edges_diagonal == 1);
    }

    // fail
    {
        lsdtt::lsdtt_ridge_extraction_options o;
        o.diagonal_policy = lsdtt::diagonal_step_policy::fail;
        bool threw = false;
        try { (void)lsdtt::extract_hilltop_pixels({rl}, dem, o); }
        catch (const std::runtime_error&) { threw = true; }
        check("diagonal fail: throws with the node named", threw);
    }
}

// =============================================================================
// 4. Ordering attributes and deduplication
// =============================================================================

void test_ordering_attributes()
{
    const auto dem = make_dem_meta(9, 9, 10.0);
    const auto dgm = make_divide_meta(dem, false, false);

    auto rl = make_ridgeline({{4,1},{4,2},{4,3},{4,4},{4,5}}, dgm);
    const auto res = lsdtt::extract_hilltop_pixels({rl}, dem);

    check("ordering: ridgeline_id set",
          !res.pixels.empty() && res.pixels.front().ridgeline_id == 0);

    // distance must increase with node index, and total length is 4 cells
    bool monotone = true;
    for (const auto& p : res.pixels)
    {
        const double want = static_cast<double>(p.ridgeline_node_index) * 10.0;
        if (std::abs(p.distance_along_ridgeline_m - want) > 1e-9) monotone = false;
    }
    check("ordering: distance_along_ridgeline_m equals node index * cell",
          monotone);

    check("ordering: total length is 4 cells",
          !res.pixels.empty() &&
          std::abs(res.pixels.front().ridgeline_total_length_m - 40.0) < 1e-9,
          res.pixels.empty() ? "" :
          std::to_string(res.pixels.front().ridgeline_total_length_m));

    check("ordering: node index spans the edges",
          res.pixels.back().ridgeline_node_index == 3,
          std::to_string(res.pixels.back().ridgeline_node_index));
}

// A divide shared by two basins is traced once per basin.  Without dedup the
// interfluve is double-weighted in every statistic.
void test_shared_divide_deduplication()
{
    const auto dem = make_dem_meta(5, 5, 10.0);
    const auto dgm = make_divide_meta(dem, false, false);

    auto a = make_ridgeline({{2,1},{2,2},{2,3}}, dgm);
    auto b = make_ridgeline({{2,3},{2,2},{2,1}}, dgm);   // same divide, reversed

    const auto res = lsdtt::extract_hilltop_pixels({a, b}, dem);

    check("shared divide: dedup keeps four unique pixels",
          res.pixels.size() == 4, std::to_string(res.pixels.size()));
    check("shared divide: the duplicates are counted",
          res.n_duplicates_dropped == 4,
          std::to_string(res.n_duplicates_dropped));

    lsdtt::lsdtt_ridge_extraction_options o;
    o.deduplicate_pixels = false;
    const auto raw = lsdtt::extract_hilltop_pixels({a, b}, dem, o);
    check("shared divide: dedup off gives every per-edge record",
          raw.pixels.size() == 8, std::to_string(raw.pixels.size()));
}

void test_basin_outlet_passthrough()
{
    const auto dem = make_dem_meta(5, 5, 10.0);
    const auto dgm = make_divide_meta(dem, false, false);
    auto rl = make_ridgeline({{2,1},{2,2}}, dgm);

    const auto res = lsdtt::extract_hilltop_pixels({rl}, dem, {}, {77});
    check("basin outlet: carried onto every pixel",
          !res.pixels.empty() &&
          std::all_of(res.pixels.begin(), res.pixels.end(),
                      [](const lsdtt::lsdtt_hilltop_pixel& p)
                      { return p.basin_outlet_flat == 77; }));

    bool threw = false;
    try { (void)lsdtt::extract_hilltop_pixels({rl}, dem, {}, {1, 2}); }
    catch (const std::invalid_argument&) { threw = true; }
    check("basin outlet: length mismatch throws", threw);
}

void test_degenerate_inputs()
{
    const auto dem = make_dem_meta(5, 5, 10.0);
    const auto dgm = make_divide_meta(dem, false, false);

    const auto empty = lsdtt::extract_hilltop_pixels({}, dem);
    check("degenerate: no ridgelines gives no pixels and no edges",
          empty.pixels.empty() && empty.n_edges_total == 0);

    auto single = make_ridgeline({{2,2}}, dgm);
    const auto one = lsdtt::extract_hilltop_pixels({single}, dem);
    check("degenerate: a one-node ridgeline has no edge and no pixels",
          one.pixels.empty() && one.n_edges_total == 0);

    auto repeated = make_ridgeline({{2,2},{2,2},{2,3}}, dgm);
    const auto rep = lsdtt::extract_hilltop_pixels({repeated}, dem);
    check("degenerate: repeated node counted and skipped",
          rep.n_edges_repeated == 1 && rep.n_edges_cardinal == 1,
          "rep=" + std::to_string(rep.n_edges_repeated) +
          " card=" + std::to_string(rep.n_edges_cardinal));
}

} // namespace

int main()
{
    std::cout << "=== hilltop flow routing: Phase 1 ridge extraction ===\n\n";

    std::cout << "-- 1. divide-grid indexing, drawn by hand --------------\n";
    test_horizontal_edge_pixels();
    test_vertical_edge_pixels();
    test_reversal_swaps_sides();
    test_corner_side_consistency();

    std::cout << "\n-- 2. boundaries and periodicity -----------------------\n";
    test_divide_on_dem_boundary_drops_outside_pixel();
    test_periodic_ew_wraps_at_the_seam();
    test_inconsistent_metas_throw();

    std::cout << "\n-- 3. step classification and diagonal policy ----------\n";
    test_step_classification();
    test_diagonal_policies();

    std::cout << "\n-- 4. ordering, dedup, degenerate inputs ---------------\n";
    test_ordering_attributes();
    test_shared_divide_deduplication();
    test_basin_outlet_passthrough();
    test_degenerate_inputs();

    std::cout << "\n=======================================================\n";
    std::cout << "checks: " << g_n_checks << "   failed: " << g_n_failed << "\n";
    return g_n_failed == 0 ? 0 : 1;
}
