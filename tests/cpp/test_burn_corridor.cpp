// test_burn_corridor.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
//
// -----------------------------------------------------------------------------
// burn_vector_network_to_corridor (stream_burn_mode burn_to_corridor,
// river-lines plan phase VT), on synthetic DEMs (80 x 40 cells, 10 m).
//
// Offset valley: the surface falls 0.05 m/m east and rises 0.1 m/m away from a
// valley floor along row 24 (y = y0 + 155); the vector line runs east along
// row 21 (y = y0 + 185), 30 m north of the floor. Expected, from the corridor
// rule (DEM path accepted while within X of the line):
//   X = 50 m: steepest descent from the line's first cell goes down to the
//     floor and along it, never more than 30 m from the line: "followed",
//     nothing lowered, largest separation 30 m;
//   X = 10 m: the descent leaves the corridor at once: "burned", and every
//     channel cell is within 15 m (the 1.5-cell end-node disc) of the line;
//   X = 0: "burned" on the line's own cells: every channel cell more than
//     20 m from both ends lies on the line (separation 0).
// Properties for every burned case: the channel never rises downstream, and a
// second pass over the burned DEM burns nothing (the capturing burn holds the
// flow on the path).
//
// Confluence: a main stem along the valley floor in two reaches, and a
// tributary straight down the hillside into the junction: both main reaches
// followed, no channel rises, steepest descent from the tributary head reaches
// the outlet end, and a second pass burns nothing.
//
// Interior end: a line across a filled lake, against the lake surface's tilt,
// ends inside the DEM below the lake surface (a pit); it is given an exit, so
// a fill of the burned DEM raises no channel cell.
//
// Tunnel: the valley runs along the line, with a ridge 8 m high across it at
// columns 40-41 and a notch in the north side at column 30. After filling, the
// lake behind the ridge spills out through the notch, so steepest descent
// leaves a 10 m corridor: "burned", the channel crosses the ridge columns, the
// deepest cut is the ridge top above the filled lake (about 3 m, worked below),
// and a second pass follows the path.
//
// Usage: lsdtt-test-burn-corridor
// -----------------------------------------------------------------------------

#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include "../../src/cpu/lsdtt_burn_corridor_cpu.hpp"
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

constexpr double k_x0 = 500000.0, k_y0 = 6000000.0, k_cs = 10.0;
constexpr int k_nc = 80, k_nr = 40;

double cx(int c) { return k_x0 + (c + 0.5) * k_cs; }
double cy(int r) { return k_y0 + (k_nr - r - 0.5) * k_cs; }

template <typename F>
lsdtt::lsdtt_raster_cpu make_dem(F height)
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
            data[static_cast<size_t>(r * k_nc + c)] = static_cast<float>(height(cx(c), cy(r), r, c));
    return lsdtt::lsdtt_raster_cpu(std::move(data), std::move(mask), meta);
}

double z_at(const lsdtt::lsdtt_raster_cpu& d, const lsdtt::lsdtt_xy& p)
{
    const int c = static_cast<int>(std::floor((p.x - k_x0) / k_cs));
    const int r = static_cast<int>(std::floor((k_y0 + k_nr * k_cs - p.y) / k_cs));
    return d.get_data_array()(r, c);
}

bool never_rises(const lsdtt::lsdtt_raster_cpu& d, const std::vector<lsdtt::lsdtt_xy>& ch)
{
    for (size_t i = 1; i < ch.size(); ++i)
        if (z_at(d, ch[i]) > z_at(d, ch[i - 1])) return false;
    return true;
}

} // namespace

int main()
{
    const double y_floor = k_y0 + 155.0;
    const double y_line = k_y0 + 185.0;
    const std::vector<std::vector<lsdtt::lsdtt_xy>> line = {{{k_x0 + 15.0, y_line}, {k_x0 + 785.0, y_line}}};
    auto valley = [&](double x, double y, int, int) { return 100.0 - 0.05 * (x - k_x0) + 0.1 * std::fabs(y - y_floor); };

    // ---- offset valley ----------------------------------------------------------
    {
        auto dem = make_dem(valley);
        lsdtt::burn_corridor_options o;
        o.corridor_half_width_metres = 50.0;
        const auto r = lsdtt::burn_vector_network_to_corridor(dem, line, o);
        check("X = 50: one reach, followed", r.n_followed == 1 && r.n_burned == 0 && r.n_skipped == 0);
        check("X = 50: nothing lowered", r.n_cells_lowered == 0, std::to_string(r.n_cells_lowered));
        check("X = 50: channel reaches the floor, 30 m from the line",
              !r.reaches.empty() && std::fabs(r.reaches[0].max_separation_m - 30.0) < 1e-6,
              r.reaches.empty() ? "" : std::to_string(r.reaches[0].max_separation_m));
    }
    for (const double X : {10.0, 0.0})
    {
        const std::string tag = "X = " + std::to_string(static_cast<int>(X)) + ": ";
        auto dem = make_dem(valley);
        lsdtt::burn_corridor_options o;
        o.corridor_half_width_metres = X;
        const auto r = lsdtt::burn_vector_network_to_corridor(dem, line, o);
        check(tag + "burned", r.n_burned == 1 && r.n_followed == 0, "burned " + std::to_string(r.n_burned));
        const auto& ch = r.reaches.empty() ? std::vector<lsdtt::lsdtt_xy>() : r.reaches[0].channel;
        check(tag + "every channel cell within 15 m of the line",
              !r.reaches.empty() && r.reaches[0].max_separation_m <= 15.0 + 1e-9,
              r.reaches.empty() ? "" : std::to_string(r.reaches[0].max_separation_m));
        if (X == 0.0)
        {
            bool on_line = true;
            for (const auto& p : ch)
                if (p.x - k_x0 > 35.0 && p.x - k_x0 < 765.0) on_line &= std::fabs(p.y - y_line) < 1e-6;
            check(tag + "away from the ends the channel is the line's own cells", on_line);
        }
        check(tag + "channel ends within 15 m of the outlet end",
              !ch.empty() && std::hypot(ch.back().x - (k_x0 + 785.0), ch.back().y - y_line) <= 15.0 + 1e-9);
        check(tag + "channel never rises downstream", never_rises(dem, ch));
        const auto r2 = lsdtt::burn_vector_network_to_corridor(dem, line, o);
        check(tag + "second pass burns nothing (the burn captured the flow)",
              r2.n_burned == 0 && r2.n_cells_lowered == 0 && r2.n_followed == 1,
              "burned " + std::to_string(r2.n_burned) + ", lowered " + std::to_string(r2.n_cells_lowered));
    }

    // ---- tunnel -------------------------------------------------------------------
    {
        auto tunnel = [&](double x, double y, int, int c)
        {
            double h = 100.0 - 0.05 * (x - k_x0) + 0.1 * std::fabs(y - y_line);
            if (c == 40 || c == 41) h = std::max(h, 100.0 - 0.05 * (x - k_x0) + 8.0);
            if (c == 30 && y > y_line) h = 100.0 - 0.05 * (x - k_x0) - 0.002 * (y - y_line);  // notch, falls north
            return h;
        };
        auto raw = make_dem(tunnel);
        auto dem = lsdtt::fill_pits(raw, 1.0e-4f, 1);
        lsdtt::burn_corridor_options o;
        o.corridor_half_width_metres = 10.0;
        const auto r = lsdtt::burn_vector_network_to_corridor(dem, line, o);
        check("tunnel: burned (the filled lake spills out through the notch)", r.n_burned == 1,
              "followed " + std::to_string(r.n_followed));
        bool crosses = false;
        for (const auto& p : r.reaches[0].channel)
            crosses |= (p.x > cx(40) - 1.0 && p.x < cx(41) + 1.0 && std::fabs(p.y - y_line) <= 10.0 + 1e-9);
        check("tunnel: the channel crosses the ridge inside the corridor", crosses);
        // On the FILLED DEM the lake behind the ridge stands at the notch's spill
        // level (column 30). The ridge top at column 40 is 8 m above a floor that
        // is 5 m (100 m x 0.05) lower than the floor at the notch, so the cut is
        // about 8 - 5 = 3 m, plus the small fill and burn gradients.
        check("tunnel: deepest cut 3.0-3.1 m (ridge top above the filled lake)",
              r.reaches[0].max_burn_depth_m >= 3.0 && r.reaches[0].max_burn_depth_m <= 3.1,
              std::to_string(r.reaches[0].max_burn_depth_m));
        check("tunnel: channel never rises downstream", never_rises(dem, r.reaches[0].channel));
        const auto r2 = lsdtt::burn_vector_network_to_corridor(dem, line, o);
        check("tunnel: second pass follows the channel", r2.n_followed == 1 && r2.n_cells_lowered == 0,
              "lowered " + std::to_string(r2.n_cells_lowered));
    }

    // ---- confluence ---------------------------------------------------------------
    // Main stem along the valley floor (row 24), in two reaches meeting at
    // x0 + 505; a tributary from (x0 + 505, y0 + 345) straight down the
    // hillside into that junction, where no valley guides it. X = 10 m.
    {
        auto dem = make_dem(valley);
        const lsdtt::lsdtt_xy J{k_x0 + 505.0, y_floor};
        const std::vector<std::vector<lsdtt::lsdtt_xy>> net = {
            {{k_x0 + 15.0, y_floor}, J}, {J, {k_x0 + 785.0, y_floor}}, {{k_x0 + 505.0, k_y0 + 345.0}, J}};
        lsdtt::burn_corridor_options o;
        o.corridor_half_width_metres = 10.0;
        const auto r = lsdtt::burn_vector_network_to_corridor(dem, net, o);
        check("confluence: 3 reaches, 4 nodes", r.reaches.size() == 3 && r.n_nodes == 4);
        check("confluence: both main-stem reaches followed",
              r.reaches[0].status == "followed" && r.reaches[1].status == "followed",
              r.reaches[0].status + "/" + r.reaches[1].status);
        bool all_down = true;
        for (const auto& rc : r.reaches) all_down &= never_rises(dem, rc.channel);
        check("confluence: no channel rises downstream", all_down);
        // Steepest descent on the final DEM from the tributary head reaches the
        // outlet end (within the 15 m end disc).
        int rr = static_cast<int>(std::floor((k_y0 + k_nr * k_cs - (k_y0 + 345.0)) / k_cs));
        int cc = static_cast<int>(std::floor(505.0 / k_cs));
        const auto& a = dem.get_data_array();
        for (int step = 0; step < 400; ++step)
        {
            int br = -1, bc = -1;
            double best = 0.0;
            for (int dr = -1; dr <= 1; ++dr)
                for (int dc = -1; dc <= 1; ++dc)
                {
                    const int r2 = rr + dr, c2 = cc + dc;
                    if ((dr == 0 && dc == 0) || r2 < 0 || r2 >= k_nr || c2 < 0 || c2 >= k_nc) continue;
                    const double s = (a(rr, cc) - a(r2, c2)) / ((dr && dc) ? std::sqrt(2.0) : 1.0);
                    if (s > best) { best = s; br = r2; bc = c2; }
                }
            if (br < 0) break;
            rr = br;
            cc = bc;
        }
        check("confluence: steepest descent from the tributary head reaches the outlet end",
              std::hypot(cx(cc) - (k_x0 + 785.0), cy(rr) - y_floor) <= 15.0 + 1e-9,
              "ended at col " + std::to_string(cc) + ", row " + std::to_string(rr));
        const auto r2 = lsdtt::burn_vector_network_to_corridor(dem, net, o);
        check("confluence: second pass burns nothing", r2.n_burned == 0 && r2.n_cells_lowered == 0,
              "burned " + std::to_string(r2.n_burned));
    }

    // ---- channel end inside the DEM -----------------------------------------------
    // A plane falling 0.05 m/m east with a 5 m hollow (radius 150 m) centred at
    // (x0 + 400, y0 + 200), filled: a lake whose surface is tilted east towards
    // its spill point. A line digitised WEST across the lake, from x0 + 505 to
    // x0 + 405, against that tilt, ending inside the DEM: every burned cell has
    // to be lowered a little more, so the end finishes below the lake surface,
    // a pit. It must be given an exit, so that filling the burned DEM raises no
    // channel cell and the channel still runs downhill.
    {
        auto bowl = [&](double x, double y, int, int)
        {
            const double rr = std::hypot(x - (k_x0 + 400.0), y - (k_y0 + 200.0));
            return 100.0 - 0.05 * (x - k_x0) - 5.0 * std::max(0.0, 1.0 - rr / 150.0);
        };
        auto raw = make_dem(bowl);
        auto dem = lsdtt::fill_pits(raw, 1.0e-4f, 1);
        const std::vector<std::vector<lsdtt::lsdtt_xy>> west = {{{k_x0 + 505.0, k_y0 + 205.0}, {k_x0 + 405.0, k_y0 + 205.0}}};
        lsdtt::burn_corridor_options o;
        o.corridor_half_width_metres = 0.0;
        const auto r = lsdtt::burn_vector_network_to_corridor(dem, west, o);
        check("interior end: burned", r.n_burned == 1);
        check("interior end: the end was a pit and was given an exit", r.n_ends_drained == 1,
              std::to_string(r.n_ends_drained) + " drained");
        const auto filled = lsdtt::fill_pits(dem, 1.0e-4f, 1);
        int raised_on_channel = 0;
        for (const auto& p : r.reaches[0].channel)
            if (z_at(filled, p) > z_at(dem, p) + 1e-6) ++raised_on_channel;
        check("interior end: a fill raises no channel cell", raised_on_channel == 0,
              std::to_string(raised_on_channel) + " raised");
        check("interior end: channel never rises downstream", never_rises(dem, r.reaches[0].channel));
    }

    std::cout << "\n" << (g_n_checks - g_n_failed) << " / " << g_n_checks << " checks passed\n";
    return g_n_failed == 0 ? 0 : 1;
}
