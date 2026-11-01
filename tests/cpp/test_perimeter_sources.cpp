// test_perimeter_sources.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
//
// -----------------------------------------------------------------------------
// find_perimeter_sources (lsdtt_perimeter_sources_cpu) on synthetic DEMs.
//
// P1 trace: rectangle length and start, padding equivalence inside a nodata
//    frame, a rotated square (closed, diagonal steps), a one-cell neck
//    (terminates, repeated visits), two components and a hole.
// P2 runs: a flat-floored valley gives one crossing per edge it cuts; a
//    channel in the floor splits it without gap closing; a bench on a tilted
//    plane is not a low spot; a plateau on the trace start is one run across
//    the seam.
// P3/P4 outlets and sources: valley through (1 source, 1 outlet on the floor
//    exit, area summed over parallel D8 outlets); trunk plus tributary (2
//    sources); a perched flat notch (failed_inward); a basin-clipped DEM
//    (0 sources); the two-network shape with n = 1 and 2 (1 source both times);
//    a coast with a river mouth.
//
// Terrain: a 400 m cap cut by straight valleys, each with a flat floor of
// half-width hw (cells) and sides at gradient 0.3; z = min over valleys.
// 10 m cells on the synthetic georeference.
// -----------------------------------------------------------------------------

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "lsdtt_raster_cpu.hpp"
#include "lsdtt_perimeter_sources_cpu.hpp"
#include "lsdtt_raster_flow_preprocessing_cpu.hpp"
#include "../../src/core/lsdtt_data_loader.hpp"
#include "../../src/core/lsdtt_flow_base.hpp"

namespace
{
int g_failures = 0;
void check(bool ok, const std::string& what)
{
    std::cout << (ok ? "PASS  " : "FAIL  ") << what << std::endl;
    if (!ok) ++g_failures;
}

constexpr double k_cell = 10.0;

struct valley { double r0, c0, r1, c1, z0, z1, hw; };

struct grid
{
    int rows, cols;
    std::vector<float> z;
    std::vector<std::uint8_t> m;
    grid(int r, int c, float z0 = 400.0f)
        : rows(r), cols(c), z(static_cast<std::size_t>(r) * c, z0),
          m(static_cast<std::size_t>(r) * c, 1) {}
    float& at(int r, int c) { return z[static_cast<std::size_t>(r) * cols + c]; }
    std::uint8_t& mask(int r, int c) { return m[static_cast<std::size_t>(r) * cols + c]; }
    void cut(const valley& v, double side_gradient = 0.3)
    {
        const double dr = v.r1 - v.r0, dc = v.c1 - v.c0, len2 = dr * dr + dc * dc;
        for (int r = 0; r < rows; ++r)
            for (int c = 0; c < cols; ++c)
            {
                double t = ((r - v.r0) * dr + (c - v.c0) * dc) / len2;
                t = std::clamp(t, 0.0, 1.0);
                const double d = std::hypot(r - (v.r0 + t * dr), c - (v.c0 + t * dc));
                const double zv = v.z0 + t * (v.z1 - v.z0) +
                                  side_gradient * k_cell * std::max(0.0, d - v.hw);
                at(r, c) = std::min(at(r, c), static_cast<float>(zv));
            }
    }
};

lsdtt::lsdtt_raster_meta meta_for(int rows, int cols)
{
    lsdtt::lsdtt_raster_meta m;
    m.n_rows = rows;
    m.n_cols = cols;
    m.cell_size = k_cell;
    m.no_data_value = -9999.0;
    m.set_default_synthetic_georeference();
    return m;
}

lsdtt::perimeter_sources_result run(const grid& g, lsdtt::perimeter_sources_request rq)
{
    const lsdtt::lsdtt_raster_meta meta = meta_for(g.rows, g.cols);
    std::vector<float> z = g.z;
    for (std::size_t i = 0; i < z.size(); ++i)
        if (!g.m[i]) z[i] = static_cast<float>(meta.no_data_value);
    lsdtt::lsdtt_raster_cpu dem(std::vector<float>(z), std::vector<std::uint8_t>(g.m), meta);
    const lsdtt::lsdtt_raster_cpu filled = lsdtt::fill_pits(dem, 0.0001f, "oooo", 1);
    lsdtt::lsdtt_flowinfo_options fo;
    fo.boundary_conditions = "oooo";
    fo.n_threads = 1;
    const auto flow = lsdtt::make_flow_base(filled.get_data_array(), filled.get_mask_array(),
                                            meta, fo);
    return lsdtt::find_perimeter_sources(dem, filled, *flow, rq);
}

int count_verdict(const lsdtt::perimeter_sources_result& r, lsdtt::perimeter_run_verdict v)
{
    return static_cast<int>(std::count_if(r.runs.begin(), r.runs.end(),
        [&](const lsdtt::perimeter_run& p) { return p.verdict == v; }));
}

bool closed_and_adjacent(const lsdtt::perimeter_sources_result& r, int& n_diagonal)
{
    n_diagonal = 0;
    const std::size_t N = r.profile.size();
    for (std::size_t k = 0; k < N; ++k)
    {
        const auto& a = r.profile[k];
        const auto& b = r.profile[(k + 1) % N];
        const int dr = std::abs(a.row - b.row), dc = std::abs(a.col - b.col);
        if (dr > 1 || dc > 1 || (dr == 0 && dc == 0)) return false;
        if (dr && dc) ++n_diagonal;
    }
    return true;
}

void describe(const std::string& name, const lsdtt::perimeter_sources_result& r)
{
    std::cout << "  [" << name << "] visits " << r.profile.size() << ", runs " << r.runs.size()
              << ", outlet runs " << r.outlet_runs.size() << ", sources " << r.sources.size()
              << "\n";
    for (const auto& p : r.runs)
        if (p.width_metres >= 100.0)
            std::cout << "    run " << p.id << " (" << r.profile[p.lowest_visit].row << ","
                      << r.profile[p.lowest_visit].col << ") width " << p.width_metres
                      << " relief " << p.relief_metres << " "
                      << lsdtt::perimeter_run_verdict_name(p.verdict)
                      << " inward_max " << p.inward_max_gradient << "\n";
    for (std::size_t k = 0; k < std::min<std::size_t>(3, r.outlet_runs.size()); ++k)
    {
        const auto& o = r.outlet_runs[k];
        std::cout << "    outlet rank " << o.rank << " area " << o.summed_drainage_area_m2
                  << " cells " << o.n_cells << " on_crossing " << o.on_crossing
                  << " outlet flat " << o.outlet_flat << "\n";
    }
}

grid make_valley()
{
    grid g(120, 160);
    g.cut({-5, 80, 125, 80, 150, 130, 25});
    return g;
}

// Trunk-and-tributary shape: trunk entering north, tributary entering west.
grid make_tributary()
{
    grid g = make_valley();
    g.cut({60, -5, 60, 55, 160, 141, 14});
    return g;
}

// Two-network shape: main valley south -> north, a second network starting
// inside and leaving east.
grid make_two_network()
{
    grid g(140, 160);
    g.cut({145, 60, -5, 60, 150, 120, 25});
    g.cut({70, 120, 70, 165, 160, 150, 17});
    return g;
}

// The CLI pins (fixture_perimeter_sources_cli) read these.
void write_cli_fixture(const std::string& dir)
{
    std::filesystem::create_directories(dir);
    lsdtt::lsdtt_data_loader loader;
    for (const auto& [name, g] : {std::make_pair(std::string("ps_tributary_dem.tif"), make_tributary()),
                                  std::make_pair(std::string("ps_two_network_dem.tif"), make_two_network())})
        loader.write_full(g.z, meta_for(g.rows, g.cols), dir + "/" + name,
                          lsdtt::raster_output_format::gtiff);
}

} // namespace

int main(int argc, char* argv[])
{
    if (argc == 3 && std::string(argv[1]) == "--write-cli-fixture")
    {
        write_cli_fixture(argv[2]);
        return 0;
    }
    using V = lsdtt::perimeter_run_verdict;

    // =========================================================================
    // P1: trace
    // =========================================================================
    {
        grid g(20, 30);
        for (int r = 0; r < 20; ++r)
            for (int c = 0; c < 30; ++c) g.at(r, c) = static_cast<float>(100 + r + 2 * c);
        const auto res = run(g, {});
        const std::size_t N = res.profile.size();
        check(N == 2 * (20 + 30) - 4, "rectangle: profile has 2(nx+ny)-4 visits (" +
              std::to_string(N) + ")");
        check(std::abs(res.report.profile_length_metres - 96 * k_cell) < 1e-6,
              "rectangle: length = 96 cells");
        float zmax = 0;
        for (const auto& v : res.profile) zmax = std::max(zmax, v.elevation);
        check(res.profile[0].elevation == zmax, "rectangle: starts at the highest cell");
        check(res.report.n_repeated_visits == 0, "rectangle: no repeated visits");

        // Same rectangle inside a 3-cell nodata frame.
        grid f(26, 36);
        for (int r = 0; r < 26; ++r)
            for (int c = 0; c < 36; ++c)
            {
                const bool in = r >= 3 && r < 23 && c >= 3 && c < 33;
                f.mask(r, c) = in ? 1 : 0;
                if (in) f.at(r, c) = g.at(r - 3, c - 3);
            }
        const auto rf = run(f, {});
        bool same = rf.profile.size() == N;
        for (std::size_t k = 0; same && k < N; ++k)
            same = rf.profile[k].elevation == res.profile[k].elevation &&
                   rf.profile[k].row == res.profile[k].row + 3 &&
                   rf.profile[k].col == res.profile[k].col + 3;
        check(same, "nodata frame: identical profile (padding equivalence)");
    }
    {
        // Square rotated 45 degrees inside nodata: staircase edges.
        grid g(80, 80);
        const double a = 45.0 * M_PI / 180.0;
        for (int r = 0; r < 80; ++r)
            for (int c = 0; c < 80; ++c)
            {
                const double x = c - 39.5, y = r - 39.5;
                const double u = x * std::cos(a) + y * std::sin(a);
                const double w = -x * std::sin(a) + y * std::cos(a);
                g.mask(r, c) = (std::abs(u) < 28 && std::abs(w) < 28) ? 1 : 0;
                g.at(r, c) = static_cast<float>(100 + r);
            }
        const auto res = run(g, {});
        int nd = 0;
        check(closed_and_adjacent(res, nd) && nd > 0,
              "rotated square: closed, 8-adjacent steps, " + std::to_string(nd) + " diagonal");
        check(res.report.n_components_dropped == 0, "rotated square: one component");
    }
    {
        // Two 12x12 blocks joined by a one-cell-wide corridor 4 cells long.
        grid g(14, 30);
        for (int r = 0; r < 14; ++r)
            for (int c = 0; c < 30; ++c)
            {
                const bool block = r >= 1 && r <= 12 && ((c >= 1 && c <= 12) || (c >= 17 && c <= 28));
                const bool neck = r == 6 && c >= 13 && c <= 16;
                g.mask(r, c) = (block || neck) ? 1 : 0;
                g.at(r, c) = static_cast<float>(100 + c);
            }
        const auto res = run(g, {});
        int nd = 0;
        check(closed_and_adjacent(res, nd), "neck: trace closed");
        check(res.report.n_repeated_visits >= 4, "neck: corridor cells visited twice (" +
              std::to_string(res.report.n_repeated_visits) + " repeats)");
    }
    {
        // Two components and a hole.
        grid g(40, 60);
        for (int r = 0; r < 40; ++r)
            for (int c = 0; c < 60; ++c)
            {
                const bool big = r >= 2 && r < 38 && c >= 2 && c < 40;
                const bool hole = r >= 15 && r < 20 && c >= 15 && c < 20;
                const bool small = r >= 5 && r < 10 && c >= 50 && c < 55;
                g.mask(r, c) = ((big && !hole) || small) ? 1 : 0;
                g.at(r, c) = static_cast<float>(100 + r + c);
            }
        const auto res = run(g, {});
        check(res.report.n_components_dropped == 1 && res.report.n_cells_dropped == 25,
              "two components: larger traced, 25-cell one dropped");
        check(res.report.n_holes == 1 && res.report.n_hole_cells == 25, "one 25-cell hole");
        bool on_big = true;
        for (const auto& v : res.profile) on_big = on_big && v.col < 40;
        check(on_big, "two components: profile on the larger only");
    }

    // =========================================================================
    // P2-P4: valleys
    // =========================================================================
    // A: valley entering north, leaving south.
    const grid A = make_valley();
    {
        const auto res = run(A, {});
        describe("valley", res);
        check(res.sources.size() == 1 && res.sources[0].row == 0,
              "valley: one source, on the north edge");
        check(!res.outlet_runs.empty() && res.outlet_runs[0].on_crossing &&
              res.outlet_runs[0].n_cells > 1,
              "valley: top outlet is the floor exit, summed over parallel outlet cells (" +
              (res.outlet_runs.empty() ? std::string("none") :
               std::to_string(res.outlet_runs[0].n_cells)) + " cells)");
        check(count_verdict(res, V::outlet_crossing) == 1, "valley: one outlet crossing");
        // Floor area is 51 x 120 cells; the exit carries at least that.
        check(!res.outlet_runs.empty() &&
              res.outlet_runs[0].summed_drainage_area_m2 >= 51.0 * 120 * k_cell * k_cell,
              "valley: outlet run area >= floor area");
    }
    {
        // Edge noise: the entering floor's lowest edge cell is a 0.3 m low
        // that drains straight out. The source moves to
        // the lowest floor cell that drains inward.
        grid g = A;
        g.at(0, 80) -= 0.3f;
        const auto res = run(g, {});
        check(res.sources.size() == 1 && res.sources[0].row == 0 && res.sources[0].col != 80,
              "edge noise: still one north source, off the draining-out cell");
    }
    {
        // Channel cut in the floor: splits each crossing without gap closing.
        grid g = A;
        for (int r = 0; r < g.rows; ++r)
            for (int c = 80; c <= 81; ++c) g.at(r, c) -= 3.0f;
        lsdtt::perimeter_sources_request rq;
        rq.gradient_radius_metres = 10.0;
        rq.gap_closing_metres = 0.0;
        auto wide = [](const lsdtt::perimeter_sources_result& r) {
            return std::count_if(r.runs.begin(), r.runs.end(),
                                 [](const lsdtt::perimeter_run& p) { return p.width_metres >= 100.0; }); };
        const auto open = run(g, rq);
        rq.gap_closing_metres = 60.0;
        const auto closed = run(g, rq);
        check(wide(open) == 4, "channel: 4 floor halves without gap closing (" +
              std::to_string(wide(open)) + ")");
        check(wide(closed) == 2 && closed.sources.size() == 1,
              "channel: 2 crossings and 1 source with gap closing (" +
              std::to_string(wide(closed)) + ", " + std::to_string(closed.sources.size()) + ")");
    }
    {
        // Bench on a tilted plane.
        grid g(100, 160);
        for (int r = 0; r < 100; ++r)
            for (int c = 0; c < 160; ++c)
                g.at(r, c) = static_cast<float>(c < 60 ? 100 + 2 * c : c <= 100 ? 220 : 220 + 2 * (c - 100));
        const auto res = run(g, {});
        describe("bench", res);
        int benches = 0;
        for (const auto& p : res.runs)
            if (p.width_metres >= 300.0) ++benches;
        check(benches == 2 && count_verdict(res, V::not_low_spot) == 2 && res.sources.empty(),
              "bench: two wide runs, both not_low_spot, no sources");
    }
    {
        // Flat plateau over the top-left corner, where the trace starts.
        grid g(60, 60);
        for (int r = 0; r < 60; ++r)
            for (int c = 0; c < 60; ++c)
                g.at(r, c) = static_cast<float>(200 - 3 * std::max({r - 19, c - 19, 0}));
        lsdtt::perimeter_sources_request rq;
        rq.gradient_radius_metres = 10.0;
        const auto res = run(g, rq);
        const std::size_t N = res.profile.size();
        const std::int32_t r0 = res.profile[0].run_id;
        const bool straddles = r0 >= 0 && res.profile[N - 1].run_id == r0;
        check(straddles && res.runs[static_cast<std::size_t>(r0)].n_visits >= 30,
              "plateau: one run across the profile seam (" +
              (r0 >= 0 ? std::to_string(res.runs[static_cast<std::size_t>(r0)].n_visits) : std::string("-")) +
              " visits)");
    }
    {
        // B: trunk plus a tributary entering from the west edge.
        const grid g = make_tributary();
        lsdtt::perimeter_sources_request rq;
        rq.min_valley_width_metres = 200.0;
        const auto res = run(g, rq);
        describe("tributary", res);
        bool west = false, north = false;
        for (const auto& s : res.sources) { west |= s.col == 0; north |= s.row == 0; }
        check(res.sources.size() == 2 && west && north,
              "tributary: two sources, north and west edges");
        check(res.sources.size() == 2 && res.sources[0].outlet_run_id == res.sources[1].outlet_run_id,
              "tributary: both drain to the same outlet");
    }
    {
        // C: perched flat notch on the east edge, dropping steeply inward.
        grid g(120, 140);
        g.cut({-5, 70, 125, 70, 150, 130, 25});
        g.cut({40, 145, 40, 120, 262, 258, 20});   // floor tilts inward
        const auto res = run(g, {});
        describe("perched notch", res);
        bool notch_failed = false;
        for (const auto& p : res.runs)
            if (p.verdict == V::failed_inward && res.profile[p.lowest_visit].col == 139)
                notch_failed = p.inward_failure_flat >= 0 && p.inward_checked_metres <= 500.0;
        check(notch_failed, "perched notch: failed_inward, failure inside 500 m");
        check(res.sources.size() == 1, "perched notch: valley source still found");
    }
    {
        // D: DEM clipped to a basin (ellipse), valley leaving south.
        grid g(120, 100);
        g.cut({10, 50, 125, 50, 200, 150, 8});
        for (int r = 0; r < 120; ++r)
            for (int c = 0; c < 100; ++c)
            {
                const double u = (r - 60.0) / 58.0, w = (c - 50.0) / 35.0;
                g.mask(r, c) = (u * u + w * w <= 1.0) ? 1 : 0;
            }
        const auto res = run(g, {});
        describe("basin clipped", res);
        check(res.sources.empty(), "basin clipped: no sources");
        check(!res.outlet_runs.empty() &&
              res.outlet_runs[0].outlet_flat / g.cols > 100,
              "basin clipped: top outlet at the south notch");
    }
    {
        // E: two-network shape: main valley south -> north, a second network
        // starting inside and leaving east.
        const grid g = make_two_network();
        lsdtt::perimeter_sources_request rq;
        rq.min_valley_width_metres = 200.0;
        const auto r1 = run(g, rq);
        describe("two_network n=1", r1);
        check(r1.sources.size() == 1 && r1.sources[0].row == g.rows - 1,
              "two_network n=1: one source, south");
        check(count_verdict(r1, V::outlet_crossing) >= 2,
              "two_network n=1: north exit and east network are outlet crossings");
        check(r1.report.outlet_area_ratio > 1.0, "two_network n=1: n/n+1 area ratio reported (" +
              std::to_string(r1.report.outlet_area_ratio) + ")");
        rq.n_outlets = 2;
        const auto r2 = run(g, rq);
        check(r2.sources.size() == 1 && r2.sources[0].row == g.rows - 1,
              "two_network n=2: still one source, south");
    }
    {
        // F: coast. Sea (nodata) south of row 85; a valley reaches it. Steep
        // coastal slope (0.08), then a gentle coastal plain (0.02).
        for (const double coast_gradient : {0.08, 0.02})
        {
            grid g(100, 160);
            for (int r = 0; r < 100; ++r)
                for (int c = 0; c < 160; ++c)
                {
                    g.mask(r, c) = r < 85 ? 1 : 0;
                    g.at(r, c) = static_cast<float>(2.0 + coast_gradient * k_cell * (84 - r) +
                                                    0.5 * std::abs(c - 80) * 0.0);
                }
            // Coastal ground rises inland; hills behind.
            for (int r = 0; r < 85; ++r)
                for (int c = 0; c < 160; ++c)
                    g.at(r, c) = std::min(400.0f, g.at(r, c) + static_cast<float>(
                                          r < 40 ? 3.0 * (40 - r) : 0.0));
            g.cut({-5, 80, 84, 80, 60, 2, 15});
            lsdtt::perimeter_sources_request rq;
            rq.min_valley_width_metres = 200.0;
            const auto res = run(g, rq);
            const std::string tag = "coast " + std::to_string(coast_gradient).substr(0, 4);
            describe(tag, res);
            const bool river_first = !res.outlet_runs.empty() &&
                std::abs(res.outlet_runs[0].outlet_flat % g.cols - 80) <= 16;
            if (coast_gradient > 0.05)
                check(river_first, tag + ": the river mouth ranks first");
            else
                std::cout << "  MEASURE " << tag << ": river mouth ranks first: "
                          << (river_first ? "yes" : "no") << "\n";
        }
    }

    std::cout << (g_failures ? "FAILURES: " + std::to_string(g_failures) : std::string("ALL PASS"))
              << std::endl;
    return g_failures ? 1 : 0;
}
