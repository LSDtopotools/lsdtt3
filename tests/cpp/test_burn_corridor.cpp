// test_burn_corridor.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
//
// -----------------------------------------------------------------------------
// burn_to_corridor, phase VT2 (river-lines plan 12.10-12.15): the burned
// profile comes from the minima of the ORIGINAL DEM under the line, and the
// rest of the DEM is carved and filled to drain to it.
//
// Synthetic DEMs, 80 x 40 cells of 10 m. The standard valley falls 0.05 m/m
// east and rises 0.1 m/m away from a floor along row 24.
//
// Every burned case is also checked for the conditioned-DEM properties
// (plan tests 3 and 8): each path falls strictly downstream; every valid
// cell not on the DEM edge has a strictly lower neighbour; the steepest
// descent from every line cell (but a path's last) is a line cell; filling
// the output, or carving then filling it, changes no cell; and the D8
// channel network at a threshold contains the line downstream of where the
// line first reaches it.
//
// Cases:
//   1-D fit: a falling floor is all minima; a 10 m hole is used at the
//     default cap and rejected (every cell of it) at 5 m; chi with theta = 0
//     equals linear; chi interpolation matches the analytic value; the
//     fit never rises above the cleaned DEM.
//   Small-pit cleaning: a 1-cell noise pit on the side of a large bowl is
//     filled and the bowl is not; a low cell beside a NoData hole drains into
//     it and is not filled.
//   Embankment (test 1): an 8 m embankment across the valley, the line along
//     the floor, 1-cell noise pits on the floor upstream. Off the pits, the
//     line upstream keeps the original floor; the pit cells are at their
//     spill level; the embankment is cut to the interpolation between the
//     floor minima either side; no off-line cell upstream is raised > 0.1 m.
//   Noise pit pulling the profile (test 2): a 2 m single-cell pit cuts the
//     line downstream of it at N = 0, not at N = 4 (the pit and its spill
//     cell, at the same level, make a 2-cell depression, so N must exceed 2).
//   Side valley (test 4): an unmapped side valley dammed 2 m drains into the
//     line by carving through the dam.
//   Parallel trench (R1): the line one cell north of the floor, a dam on the
//     floor: no bank-zone cell is carved (the pit beside the line is filled).
//   Confluence (test 5): a tributary with a scour below the main stem at the
//     junction ends exactly on the main-stem value (raised); the main stem is
//     unchanged by the tributary; rank 0 is the main stem; the tributary's
//     flow distance at the join equals the main stem's there.
//   Distributary: starts at or below the main stem's value at the node.
//   Idempotence (test 6): a second pass changes no line cell.
//   Offset valley (tests 7, 10): the line 30 m north of the floor. X = 50 m:
//     the channel is in the floor, no banks; X = 10 m: within 15 m of the
//     line; X = 0: on the line, banks raised.
//   Interior end: a line ending in a pond behind a dam is given an exit.
//   Cut cap in 2-D (12.14): a 10 m deep, 30-cell hole on the line (N = 16):
//     used at the default cap, rejected at 5 m (cut <= 5 m), linear and chi.
//   NoData on the path (R4): a 3-cell NoData hole on the line: X = 20 m goes
//     round it; X = 0 skips the reach (no_corridor_path).
//
// Usage: lsdtt-test-burn-corridor
// -----------------------------------------------------------------------------

#include <algorithm>
#include <cmath>
#include <functional>
#include <iostream>
#include <set>
#include <string>
#include <vector>

#include "../../src/cpu/lsdtt_burn_corridor_cpu.hpp"
#include "../../src/cpu/lsdtt_raster_breach_least_cost_cpu.hpp"
#include "../../src/cpu/lsdtt_raster_cpu.hpp"
#include "../../src/cpu/lsdtt_raster_flow_preprocessing_cpu.hpp"

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

std::string num(double v) { return std::to_string(v); }

constexpr double k_x0 = 500000.0, k_y0 = 6000000.0, k_cs = 10.0;
constexpr int k_nc = 80, k_nr = 40;
constexpr int k_floor_row = 24;

double cx(int c) { return k_x0 + (c + 0.5) * k_cs; }
double cy(int r) { return k_y0 + (k_nr - r - 0.5) * k_cs; }
int fl(int r, int c) { return r * k_nc + c; }

using height_fn = std::function<double(int r, int c)>;

lsdtt::lsdtt_raster_cpu make_dem(const height_fn& height, const std::function<bool(int, int)>& nodata = nullptr)
{
    lsdtt::lsdtt_raster_meta meta;
    meta.n_rows = k_nr;
    meta.n_cols = k_nc;
    meta.x_min = k_x0;
    meta.y_min = k_y0;
    meta.cell_size = k_cs;
    meta.no_data_value = -9999.0f;
    meta.set_crs_from_epsg(32630);
    std::vector<float> data(static_cast<size_t>(k_nr * k_nc));
    std::vector<uint8_t> mask(data.size(), 1);
    for (int r = 0; r < k_nr; ++r)
        for (int c = 0; c < k_nc; ++c)
        {
            const size_t i = static_cast<size_t>(fl(r, c));
            if (nodata && nodata(r, c)) { data[i] = meta.no_data_value; mask[i] = 0; }
            else data[i] = static_cast<float>(height(r, c));
        }
    return lsdtt::lsdtt_raster_cpu(std::move(data), std::move(mask), meta);
}

// The standard valley: floor along `floor_row`, falling east.
double valley(int r, int c, int floor_row = k_floor_row)
{
    return 100.0 - 0.05 * k_cs * c + 0.1 * k_cs * std::abs(r - floor_row);
}

lsdtt::lsdtt_xy pt(int r, int c) { return {cx(c), cy(r)}; }
std::vector<lsdtt::lsdtt_xy> row_line(int r, int c0, int c1) { return {pt(r, c0), pt(r, c1)}; }

struct run_out
{
    lsdtt::burn_corridor_result r;
    lsdtt::lsdtt_raster_cpu out;
};

run_out run(const lsdtt::lsdtt_raster_cpu& dem, const std::vector<std::vector<lsdtt::lsdtt_xy>>& lines,
            const lsdtt::burn_corridor_options& o)
{
    run_out ro;
    ro.r = lsdtt::burn_vector_network_to_corridor(dem, ro.out, lines, {}, o);
    return ro;
}

float zat(const lsdtt::lsdtt_raster_cpu& d, std::int64_t f) { return d.get_data()[static_cast<size_t>(f)]; }
bool ok_cell(const lsdtt::lsdtt_raster_cpu& d, int r, int c)
{
    return r >= 0 && r < k_nr && c >= 0 && c < k_nc && d.get_mask()[static_cast<size_t>(fl(r, c))];
}
bool on_edge(const lsdtt::lsdtt_raster_cpu& d, int r, int c)
{
    for (int dr = -1; dr <= 1; ++dr)
        for (int dc = -1; dc <= 1; ++dc)
            if ((dr || dc) && !ok_cell(d, r + dr, c + dc)) return true;
    return false;
}
// Steepest-descent receiver (distance weighted, as the engine), -1 if none lower.
std::int64_t receiver(const lsdtt::lsdtt_raster_cpu& d, std::int64_t f)
{
    const int r = static_cast<int>(f / k_nc), c = static_cast<int>(f % k_nc);
    double best = 0.0;
    std::int64_t rec = -1;
    for (int dr = -1; dr <= 1; ++dr)
        for (int dc = -1; dc <= 1; ++dc)
        {
            if (!(dr || dc) || !ok_cell(d, r + dr, c + dc)) continue;
            const std::int64_t v = fl(r + dr, c + dc);
            const double s = (static_cast<double>(zat(d, f)) - zat(d, v)) / ((dr && dc) ? std::sqrt(2.0) : 1.0);
            if (s > best) { best = s; rec = v; }
        }
    return rec;
}

std::set<std::int64_t> line_cells(const lsdtt::burn_corridor_result& r)
{
    std::set<std::int64_t> s;
    for (const auto& p : r.paths) s.insert(p.cells.begin(), p.cells.end());
    return s;
}

// The conditioned-DEM properties (plan tests 3 and 8).
void check_conditioned(const std::string& tag, const lsdtt::lsdtt_raster_cpu& dem, const run_out& ro)
{
    const auto& out = ro.out;
    const auto lines = line_cells(ro.r);
    bool rises = false, profile_gradient = true;
    for (const auto& p : ro.r.paths)
        for (size_t i = 1; i < p.cells.size(); ++i)
        {
            if (!(zat(out, p.cells[i]) < zat(out, p.cells[i - 1]))) rises = true;
            const bool pinned_end = (i + 1 == p.cells.size() && p.joins_path_id >= 0);
            if (!p.pin_conflict && !pinned_end &&
                p.z_profile[i - 1] - p.z_profile[i] < 1.0e-4 * (p.distance_m[i] - p.distance_m[i - 1]) - 1e-9)
                profile_gradient = false;
        }
    check(tag + ": every path falls strictly downstream", !rises);
    check(tag + ": profile falls at least min_gradient per metre", profile_gradient);
    long n_sink = 0;
    for (int r = 0; r < k_nr; ++r)
        for (int c = 0; c < k_nc; ++c)
            if (ok_cell(out, r, c) && !on_edge(out, r, c) && receiver(out, fl(r, c)) < 0) ++n_sink;
    check(tag + ": every interior cell has a lower neighbour", n_sink == 0, num(n_sink) + " sinks");
    long n_off = 0;
    for (const auto& p : ro.r.paths)
        for (size_t i = 0; i + 1 < p.cells.size(); ++i)
        {
            const int r = static_cast<int>(p.cells[i] / k_nc), c = static_cast<int>(p.cells[i] % k_nc);
            if (on_edge(out, r, c)) continue;
            if (!lines.count(receiver(out, p.cells[i]))) ++n_off;
        }
    check(tag + ": steepest descent from every line cell stays on the line", n_off == 0, num(n_off));
    {
        const auto filled = lsdtt::fill_pits(out, 1.0e-4f, 1);
        lsdtt::breach_least_cost_options bo;
        const auto carved = lsdtt::fill_pits(lsdtt::breach_depressions_least_cost(out, bo), 1.0e-4f, 1);
        long nf = 0, nc = 0;
        for (size_t i = 0; i < out.get_data().size(); ++i)
        {
            if (!out.get_mask()[i]) continue;
            if (filled.get_data()[i] != out.get_data()[i]) ++nf;
            if (carved.get_data()[i] != out.get_data()[i]) ++nc;
        }
        check(tag + ": filling the output changes no cell", nf == 0, num(nf));
        check(tag + ": carving then filling the output changes no cell", nc == 0, num(nc));
    }
    {
        // D8 accumulation (cells), highest first; then the channel network at
        // 20 cells must contain each path from where it first reaches 20.
        std::vector<std::int64_t> cells;
        std::vector<double> a(out.get_data().size(), 0.0);
        for (std::int64_t f = 0; f < static_cast<std::int64_t>(a.size()); ++f)
            if (out.get_mask()[static_cast<size_t>(f)]) { cells.push_back(f); a[static_cast<size_t>(f)] = 1.0; }
        std::stable_sort(cells.begin(), cells.end(), [&](std::int64_t x, std::int64_t y) { return zat(out, x) > zat(out, y); });
        for (const auto f : cells)
        {
            const auto v = receiver(out, f);
            if (v >= 0) a[static_cast<size_t>(v)] += a[static_cast<size_t>(f)];
        }
        bool ok = true;
        for (const auto& p : ro.r.paths)
        {
            bool in = false;
            for (const auto f : p.cells)
            {
                const bool big = a[static_cast<size_t>(f)] >= 20.0;
                if (in && !big) ok = false;
                in = in || big;
            }
        }
        check(tag + ": the D8 network (20 cells) contains the line downstream of where it joins it", ok);
    }
    (void)dem;
}

} // namespace

int main()
{
    // ===== 1-D fit =================================================================
    {
        lsdtt::burn_profile_input in;
        for (int i = 0; i < 50; ++i)
        {
            in.s.push_back(10.0 * i);
            in.zc.push_back(100.0 - 0.5 * i);
        }
        in.x = in.s;
        auto o = lsdtt::fit_burn_profile(in);
        double dmax = 0.0;
        for (size_t i = 0; i < in.zc.size(); ++i) dmax = std::max(dmax, std::fabs(o.z[i] - in.zc[i]));
        check("fit: falling floor, every cell a minimum, profile = floor", dmax < 1e-9 && o.n_record_lows == 50, num(dmax));

        for (int i = 20; i <= 22; ++i) in.zc[static_cast<size_t>(i)] -= 10.0;
        o = lsdtt::fit_burn_profile(in);
        check("fit: 10 m hole used at the default cap, river below it held down",
              o.record_low[22] && o.z[30] <= in.zc[22] + 1e-9 && o.max_cut_m > 9.0, num(o.z[30]));
        in.max_cut = 5.0;
        o = lsdtt::fit_burn_profile(in);
        check("fit: every cell of the hole rejected at a 5 m cap, profile passes over it",
              o.rejected_low[20] && o.rejected_low[21] && o.rejected_low[22] && o.n_rejected_lows == 3 &&
                  o.z[21] > in.zc[21] + 9.0 && o.max_cut_m < 5.0 && std::fabs(o.z[30] - in.zc[30]) < 1e-9,
              "cut " + num(o.max_cut_m));
        bool above = false;
        for (size_t i = 0; i < in.zc.size(); ++i)
            if (!o.rejected_low[i] && o.z[i] > in.zc[i] + 1e-9) above = true;
        check("fit: never above the cleaned DEM except at rejected minima", !above);
        in.max_cut = 2000.0;
    }
    {
        // chi: anchors at 0 and 20 only (everything between stands higher)
        lsdtt::burn_profile_input in;
        const double theta = 0.5;
        double chi = 0.0;
        for (int i = 0; i <= 20; ++i)
        {
            in.s.push_back(10.0 * i);
            in.zc.push_back(i == 0 ? 100.0 : (i == 20 ? 90.0 : 120.0));
            if (i > 0) chi += std::pow(1.0 / (1000.0 * (i + 1)), theta) * 10.0;
            in.x.push_back(chi);
        }
        in.min_gradient = 1e-6;
        const auto o = lsdtt::fit_burn_profile(in);
        double dmax = 0.0;
        for (int i = 0; i <= 20; ++i)
            dmax = std::max(dmax, std::fabs(o.z[static_cast<size_t>(i)] - (100.0 - 10.0 * in.x[static_cast<size_t>(i)] / chi)));
        check("fit: chi interpolation matches the analytic value (1e-4 m)", dmax < 1e-4, num(dmax));
        bool concave = true;   // slopes fall downstream when A grows
        for (int i = 2; i <= 20; ++i)
            if (o.z[static_cast<size_t>(i - 2)] - o.z[static_cast<size_t>(i - 1)] <
                o.z[static_cast<size_t>(i - 1)] - o.z[static_cast<size_t>(i)] - 1e-12) concave = false;
        check("fit: chi profile concave up where area grows", concave);
        auto lin = in;
        lin.x = lin.s;
        auto zero = in;
        double c0 = 0.0;
        for (int i = 0; i <= 20; ++i) { if (i) c0 += std::pow(1.0 / (1000.0 * (i + 1)), 0.0) * 10.0; zero.x[static_cast<size_t>(i)] = c0; }
        const auto ol = lsdtt::fit_burn_profile(lin), oz = lsdtt::fit_burn_profile(zero);
        double dd = 0.0;
        for (size_t i = 0; i < ol.z.size(); ++i) dd = std::max(dd, std::fabs(ol.z[i] - oz.z[i]));
        check("fit: chi with theta = 0 equals linear", dd < 1e-12, num(dd));
    }

    // ===== small-pit cleaning ==========================================================
    {
        // A bowl (paraboloid, 6 m deep, centre (20, 40)) with a 1-cell noise pit
        // on its side at (20, 34); a NoData hole at (5, 70) with a low cell
        // beside it at (5, 69).
        auto bowl = [](int r, int c)
        {
            const double d2 = (r - 20.0) * (r - 20.0) + (c - 40.0) * (c - 40.0);
            double z = 100.0 + std::min(6.0, 6.0 * d2 / 100.0);
            if (r == 20 && c == 34) z -= 1.0;
            if (r == 5 && c == 69) z -= 3.0;
            return z;
        };
        const auto dem = make_dem(bowl, [](int r, int c) { return r == 5 && c == 70; });
        const auto f = lsdtt::fill_small_depressions(dem, 16);
        check("pits: the noise pit on the bowl side is filled",
              f.pit_id[static_cast<size_t>(fl(20, 34))] >= 0 && zat(f.cleaned, fl(20, 34)) > zat(dem, fl(20, 34)));
        check("pits: the bowl (>= 16 cells) is not filled", zat(f.cleaned, fl(20, 40)) == zat(dem, fl(20, 40)));
        check("pits: a low cell beside a NoData hole drains into it (not filled)",
              zat(f.cleaned, fl(5, 69)) == zat(dem, fl(5, 69)));
        check("pits: exactly one depression filled", f.n_pits == 1, num(static_cast<double>(f.n_pits)));
        const auto f0 = lsdtt::fill_small_depressions(dem, 0);
        check("pits: N = 0 fills nothing", f0.n_cells == 0);
    }

    lsdtt::burn_corridor_options base;

    // ===== embankment, with noise pits (test 1) =========================================
    {
        const std::set<int> pits{10, 20, 30};
        auto h = [&](int r, int c)
        {
            double z = valley(r, c);
            if (c == 40 || c == 41) z = std::max(z, valley(k_floor_row, 40) + 8.0);
            if (r == k_floor_row && pits.count(c)) z -= 1.5;
            return z;
        };
        const auto dem = make_dem(h);
        const auto ro = run(dem, {row_line(k_floor_row, 1, 78)}, base);
        check("embankment: one reach fitted, one path", ro.r.n_fitted == 1 && ro.r.paths.size() == 1);
        const auto& p = ro.r.paths.at(0);
        double up_dev = 0.0, pit_dev = 0.0, cut_dev = 0.0, lake = 0.0;
        bool pits_flagged = true;
        for (size_t i = 0; i < p.cells.size(); ++i)
        {
            const int c = static_cast<int>(p.cells[i] % k_nc);
            const double zo = zat(dem, p.cells[i]), zf = zat(ro.out, p.cells[i]);
            if (c < 40 && pits.count(c)) { pit_dev = std::max(pit_dev, std::fabs(zf - p.z_cleaned[i])); pits_flagged &= p.z_cleaned[i] > zo; }
            else if (c < 39 && !pits.count(c - 1)) up_dev = std::max(up_dev, std::fabs(zf - zo));
            if (c == 40 || c == 41)
            {
                const double z39 = zat(dem, fl(k_floor_row, 39)), z42 = zat(dem, fl(k_floor_row, 42));
                cut_dev = std::max(cut_dev, std::fabs(zf - (z39 + (z42 - z39) * (c - 39) / 3.0)));
            }
        }
        for (int r = 18; r <= 30; ++r)
            for (int c = 1; c < 40; ++c)
                if (r != k_floor_row) lake = std::max(lake, static_cast<double>(zat(ro.out, fl(r, c)) - zat(dem, fl(r, c))));
        check("embankment: line upstream keeps the original floor off the pits (1 cm)", up_dev < 0.01, num(up_dev));
        check("embankment: pit cells at their spill level, flagged as filled", pit_dev < 0.01 && pits_flagged, num(pit_dev));
        check("embankment: cut to the interpolation between the floor minima", cut_dev < 0.01, num(cut_dev));
        check("embankment: no lake (no off-line cell upstream raised > 0.1 m)", lake < 0.1, num(lake));
        check_conditioned("embankment", dem, ro);

        // idempotence (test 6)
        const auto ro2 = run(ro.out, {row_line(k_floor_row, 1, 78)}, base);
        double dz = 0.0;
        for (const auto f : line_cells(ro.r)) dz = std::max(dz, std::fabs(static_cast<double>(zat(ro2.out, f)) - zat(ro.out, f)));
        check("embankment: a second pass changes no line cell (1e-4 m)", dz < 1e-4, num(dz));
    }

    // ===== a noise pit pulls the profile at N = 0, not at N = 4 (test 2) ====================
    {
        auto h = [](int r, int c) { return valley(r, c) - ((r == k_floor_row && c == 20) ? 2.0 : 0.0); };
        const auto dem = make_dem(h);
        auto cut_below = [&](int N)
        {
            auto o = base;
            o.max_pit_pixels = N;
            const auto ro = run(dem, {row_line(k_floor_row, 1, 78)}, o);
            double cut = 0.0;
            for (const auto f : ro.r.paths.at(0).cells)
            {
                const int c = static_cast<int>(f % k_nc);
                if (c >= 21 && c <= 25) cut = std::max(cut, static_cast<double>(zat(dem, f) - zat(ro.out, f)));
            }
            return cut;
        };
        const double c0 = cut_below(0), c2 = cut_below(4);
        check("noise pit: N = 0 holds the line downstream down (> 1 m)", c0 > 1.0, num(c0));
        check("noise pit: N = 4 does not (< 0.3 m)", c2 < 0.3, num(c2));
    }

    // ===== an unmapped side valley drains to the line by carving (test 4) ==================
    {
        // A side valley down column 30 from the north, dammed 2 m at row 12.
        auto h = [](int r, int c)
        {
            double z = valley(r, c);
            if (r < k_floor_row)
            {
                const double side = 100.0 - 0.05 * k_cs * 30 + 0.05 * k_cs * (k_floor_row - r) + 0.2 * k_cs * std::abs(c - 30);
                z = std::min(z, side);
                if (r == 12 && std::abs(c - 30) <= 1) z += 2.0;
            }
            return z;
        };
        const auto dem = make_dem(h);
        const auto ro = run(dem, {row_line(k_floor_row, 1, 78)}, base);
        const auto lines = line_cells(ro.r);
        std::int64_t f = fl(2, 30);
        for (int step = 0; step < 200 && f >= 0 && !lines.count(f); ++step) f = receiver(ro.out, f);
        check("side valley: steepest descent from its head reaches the line", f >= 0 && lines.count(f));
        check("side valley: the dam is carved, not filled", zat(ro.out, fl(12, 30)) < zat(dem, fl(12, 30)),
              num(zat(ro.out, fl(12, 30))) + " vs " + num(zat(dem, fl(12, 30))));
        check_conditioned("side valley", dem, ro);
    }

    // ===== no breach path along the line (R1) ===========================================
    {
        // The line on row 23, one cell north of the floor (row 24); a 3 m dam
        // across the floor and the south side at column 50 makes a pit beside
        // the line at (24, 49).
        auto h = [](int r, int c)
        {
            double z = valley(r, c);
            if (c == 50 && r >= k_floor_row) z += 3.0;
            return z;
        };
        const auto dem = make_dem(h);
        const auto ro = run(dem, {row_line(23, 1, 78)}, base);
        const auto lines = line_cells(ro.r);
        long carved_bank = 0;
        for (int r = 0; r < k_nr; ++r)
            for (int c = 0; c < k_nc; ++c)
            {
                if (lines.count(fl(r, c))) continue;
                bool bank = false;
                for (int dr = -1; dr <= 1; ++dr)
                    for (int dc = -1; dc <= 1; ++dc)
                        if (ok_cell(dem, r + dr, c + dc) && lines.count(fl(r + dr, c + dc))) bank = true;
                if (bank && zat(ro.out, fl(r, c)) < zat(dem, fl(r, c)) - 1e-6f) ++carved_bank;
            }
        check("parallel trench: no cell beside the line is carved", carved_bank == 0, num(carved_bank));
        check("parallel trench: the pit beside the line is filled", zat(ro.out, fl(24, 49)) > zat(dem, fl(24, 49)));
        check_conditioned("parallel trench", dem, ro);
    }

    // ===== confluence (test 5), flow distance and rank ==========================================
    {
        // Main stem along the floor in two reaches meeting at column 40; a
        // tributary down column 40 from the north, with a 3 m scour in its last
        // four cells (rows 20-23), below the main stem at the junction.
        auto h = [](int r, int c)
        {
            double z = valley(r, c);
            if (r < k_floor_row)
            {
                const double side = valley(k_floor_row, 40) + 0.05 * k_cs * (k_floor_row - r) + 0.2 * k_cs * std::abs(c - 40);
                z = std::min(z, side);
                if (c == 40 && r >= 20) z -= 3.0;
            }
            return z;
        };
        const auto dem = make_dem(h);
        auto o = base;
        o.max_pit_pixels = 0;
        const std::vector<std::vector<lsdtt::lsdtt_xy>> net = {
            row_line(k_floor_row, 1, 40), row_line(k_floor_row, 40, 78), {pt(1, 40), pt(k_floor_row, 40)}};
        const auto ro = run(dem, net, o);
        const auto main_only = run(dem, {net[0], net[1]}, o);
        check("confluence: two paths", ro.r.paths.size() == 2, num(static_cast<double>(ro.r.paths.size())));
        const lsdtt::burn_corridor_path* mainp = nullptr;
        const lsdtt::burn_corridor_path* trib = nullptr;
        for (const auto& p : ro.r.paths) (p.joins_path_id >= 0 ? trib : mainp) = &p;
        if (mainp && trib)
        {
            check("confluence: rank 0 is the main stem", mainp->path_rank == 0 && trib->path_rank == 1);
            check("confluence: the tributary joins the main stem", trib->joins_path_id == mainp->path_id);
            const auto join = trib->cells.back();
            size_t jm = 0;
            while (jm < mainp->cells.size() && mainp->cells[jm] != join) ++jm;
            check("confluence: the tributary ends on a main-stem cell", jm < mainp->cells.size());
            if (jm < mainp->cells.size())
            {
                check("confluence: the tributary ends exactly at the main-stem value",
                      std::fabs(trib->z_profile.back() - mainp->z_profile[jm]) < 1e-9);
                check("confluence: flow distance at the join equals the main stem's",
                      std::fabs(trib->flow_distance_m.back() - mainp->flow_distance_m[jm]) < 1e-6 &&
                          std::fabs(trib->join_distance_m - mainp->distance_m[jm]) < 1e-6);
            }
            check("confluence: the tributary is raised near the join (scour below the main stem)",
                  trib->max_raise_m > 1.0 && ro.r.n_line_cells_raised > 0, num(trib->max_raise_m));
            double dz = 0.0;
            for (const auto f : main_only.r.paths.at(0).cells)
                dz = std::max(dz, std::fabs(static_cast<double>(zat(ro.out, f)) - zat(main_only.out, f)));
            check("confluence: the main stem is unchanged by the tributary", dz < 1e-6, num(dz));
        }
        check_conditioned("confluence", dem, ro);
    }

    // ===== distributary =============================================================
    {
        // The main stem along the floor; a second line leaving it at column 40
        // south-east to the bottom edge.
        auto h = [](int r, int c)
        {
            double z = valley(r, c);
            if (r > k_floor_row)
            {
                const double arm = valley(k_floor_row, 40) + 0.2 * k_cs * std::abs((c - 40) - (r - k_floor_row)) -
                                   0.05 * k_cs * (r - k_floor_row);
                z = std::min(z, arm);
            }
            return z;
        };
        const auto dem = make_dem(h);
        const std::vector<std::vector<lsdtt::lsdtt_xy>> net = {
            row_line(k_floor_row, 1, 40), row_line(k_floor_row, 40, 78), {pt(k_floor_row, 40), pt(39, 55)}};
        const auto ro = run(dem, net, base);
        const lsdtt::burn_corridor_path* d = nullptr;
        for (const auto& p : ro.r.paths) if (p.starts_on_path_id >= 0) d = &p;
        check("distributary: a path starts on the main stem", d != nullptr);
        if (d)
        {
            const auto& m = ro.r.paths.at(static_cast<size_t>(d->starts_on_path_id));
            size_t j = 0;
            while (j < m.cells.size() && m.cells[j] != d->cells.front()) ++j;
            bool below = j < m.cells.size();
            for (size_t i = 0; below && i < d->cells.size(); ++i)
                if (d->z_profile[i] > m.z_profile[j] + 1e-9) below = false;
            check("distributary: at or below the main stem's value at the node", below);
        }
        check_conditioned("distributary", dem, ro);
    }

    // ===== offset valley at X = 50 / 10 / 0 (tests 7 and 10) ==============================
    {
        const auto dem = make_dem([](int r, int c) { return valley(r, c); });
        const std::vector<std::vector<lsdtt::lsdtt_xy>> line = {row_line(21, 1, 78)};   // 30 m north of the floor
        for (const double X : {50.0, 10.0, 0.0})
        {
            auto o = base;
            o.corridor_half_width_metres = X;
            const auto ro = run(dem, line, o);
            const std::string tag = "offset valley X = " + num(X).substr(0, 4);
            const auto& p = ro.r.paths.at(0);
            double max_sep = 0.0, mid_sep = 0.0;
            long in_floor = 0, mid = 0;
            for (size_t i = 0; i < p.cells.size(); ++i)
            {
                const int r = static_cast<int>(p.cells[i] / k_nc), c = static_cast<int>(p.cells[i] % k_nc);
                const double sep = std::fabs(cy(r) - cy(21));
                max_sep = std::max(max_sep, sep);
                if (c > 6 && c < 73) { ++mid; mid_sep = std::max(mid_sep, sep); if (r == k_floor_row) ++in_floor; }
            }
            if (X == 50.0)
            {
                check(tag + ": the channel follows the valley floor", in_floor == mid && mid > 0, num(in_floor) + "/" + num(mid));
                check(tag + ": no banks (the corridor reaches the valley)", ro.r.n_bank_cells == 0, num(ro.r.n_bank_cells));
            }
            else if (X == 10.0)
            {
                check(tag + ": channel within 15 m of the line", max_sep <= 15.0 + 1e-6, num(max_sep));
            }
            else
            {
                check(tag + ": channel on the line away from the ends", mid_sep < 1e-6, num(mid_sep));
                check(tag + ": banks raised", ro.r.n_bank_cells > 0, num(ro.r.n_bank_cells));
            }
            check_conditioned(tag, dem, ro);
        }
    }

    // ===== interior end ==============================================================
    {
        // A 5 m dam across the valley at column 50; the line ends at column 45
        // in the pond behind it.
        auto h = [](int r, int c) { return valley(r, c) + (c == 50 ? 5.0 : 0.0); };
        const auto dem = make_dem(h);
        auto o = base;
        o.max_pit_pixels = 0;
        const auto ro = run(dem, {row_line(k_floor_row, 1, 45)}, o);
        check("interior end: drained", ro.r.n_ends_drained == 1 && ro.r.n_cells_daylighted > 0,
              num(ro.r.n_ends_drained) + " ends, " + num(ro.r.n_cells_daylighted) + " cells");
        check_conditioned("interior end", dem, ro);
    }

    // ===== cut cap in 2-D (12.14) =========================================================
    {
        // A 10 m deep hole of 30 cells (rows 23-25, columns 30-39) on the floor.
        auto h = [](int r, int c) { return valley(r, c) - ((r >= 23 && r <= 25 && c >= 30 && c <= 39) ? 10.0 : 0.0); };
        const auto dem = make_dem(h);
        for (const auto fit : {lsdtt::burn_profile_fit::linear, lsdtt::burn_profile_fit::chi})
        {
            const std::string tag = std::string("cut cap (") + lsdtt::burn_profile_fit_name(fit) + ")";
            auto o = base;
            o.fit = fit;
            const auto r0 = run(dem, {row_line(k_floor_row, 1, 78)}, o);
            const auto& p0 = r0.r.paths.at(0);
            check(tag + ": default cap uses the hole, river downstream cut > 5 m", p0.max_cut_m > 5.0 && p0.n_rejected_lows == 0,
                  num(p0.max_cut_m));
            o.max_cut_metres = 5.0;
            const auto r5 = run(dem, {row_line(k_floor_row, 1, 78)}, o);
            const auto& p5 = r5.r.paths.at(0);
            check(tag + ": 5 m cap rejects the hole, cut <= 5 m", p5.n_rejected_lows > 0 && p5.max_cut_m <= 5.0 + 1e-6 &&
                  !p5.cut_cap_exceeded, num(p5.n_rejected_lows) + " rejected, cut " + num(p5.max_cut_m));
            check_conditioned(tag + ", 5 m cap", dem, r5);
        }
    }

    // ===== NoData on the path (R4) ==================================================
    {
        const auto dem = make_dem([](int r, int c) { return valley(r, c); },
                                  [](int r, int c) { return r == k_floor_row && c >= 30 && c <= 32; });
        auto o = base;
        o.corridor_half_width_metres = 20.0;
        const auto r20 = run(dem, {row_line(k_floor_row, 1, 78)}, o);
        bool clean = !r20.r.paths.empty();
        for (const auto& p : r20.r.paths)
            for (const auto f : p.cells) if (!dem.get_mask()[static_cast<size_t>(f)]) clean = false;
        check("nodata hole: X = 20 m goes round it", r20.r.n_fitted == 1 && clean);
        check_conditioned("nodata hole X = 20", dem, r20);
        const auto r0 = run(dem, {row_line(k_floor_row, 1, 78)}, base);
        check("nodata hole: X = 0 skips the reach (no_corridor_path)",
              r0.r.n_skipped == 1 && r0.r.reaches.at(0).skip_reason == "no_corridor_path",
              r0.r.reaches.at(0).skip_reason);
    }

    std::cout << "\n" << (g_n_checks - g_n_failed) << " / " << g_n_checks << " checks passed\n";
    return g_n_failed == 0 ? 0 : 1;
}
