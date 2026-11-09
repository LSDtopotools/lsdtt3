// test_sea_mask.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
//
// -----------------------------------------------------------------------------
// lsdtt_sea_mask_cpu: the sea from polygons or an open coastline.
//
// A 40 x 40 grid of 10 m cells. A straight coast runs from the north edge to
// the south edge, digitised north to south, so the sea (east of it) is on its
// LEFT. Land west of it stands 10 m, with a polder of -3 m behind the coast
// (below the sea, which stands -1 m). The line cells are found by sampling
// the line every 0.1 cell and marking each cell crossed together with its
// 4-neighbour joins (a chain a 4-connected flood cannot cross, as
// rasterise ALL_TOUCHED gives).
//
// Expected:
//   - lines, auto: the sea is exactly the cells whose centres lie east of the
//     line; the polder stays; evidence median -1 against 10 m;
//   - the same line digitised south to north: the same sea;
//   - sea_left on the north-to-south line: the same; sea_right: the west;
//   - a 6-cell gap in the line, or a line stopping half way: a leak, the
//     fallback within the distance band, reported;
//   - a flat DEM: auto refuses;
//   - polygons covering the sea: auto picks sea_inside; covering the land:
//     land_inside; both give the east; explicit sides as given;
//   - the sea already NoData and the polygon over it: auto decides by NoData;
//   - a side for the other geometry refused;
//   - an enclosed lagoon (a sea polygon not touching the edge) in a bowl:
//     after apply_sea_mask, fill_pits leaves the land beside it unraised, so
//     the masked sea is an outlet.
//
// Usage: lsdtt-test-sea-mask
// -----------------------------------------------------------------------------

#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "../../src/cpu/lsdtt_raster_cpu.hpp"
#include "../../src/cpu/lsdtt_raster_flow_preprocessing_cpu.hpp"
#include "../../src/cpu/lsdtt_sea_mask_cpu.hpp"

namespace
{
int g_n = 0, g_fail = 0;
void check(const std::string& name, bool ok, const std::string& detail = "")
{
    ++g_n;
    if (!ok) ++g_fail;
    std::cout << (ok ? "PASS" : "FAIL") << "  " << name;
    if (!ok && !detail.empty()) std::cout << "   [" << detail << "]";
    std::cout << "\n";
}
std::string s(double v) { return std::to_string(v); }

constexpr int N = 40;
constexpr double CS = 10.0, X0 = 1000.0, Y0 = 5000.0;
const double YMAX = Y0 + N * CS;
using lsdtt::lsdtt_xy;

lsdtt_xy centre(int r, int c) { return {X0 + (c + 0.5) * CS, YMAX - (r + 0.5) * CS}; }
// the coast: from (col 12.3, north edge) to (col 27.7, south edge)
const lsdtt_xy A{X0 + 12.3 * CS, YMAX}, B{X0 + 27.7 * CS, Y0};
// > 0 east of the line A -> B (that is its left)
double east(const lsdtt_xy& p) { return (B.x - A.x) * (p.y - A.y) - (B.y - A.y) * (p.x - A.x); }
bool in_polder(int r, int c) { return r >= 10 && r <= 20 && c >= 3 && c <= 8; }

lsdtt::lsdtt_raster_cpu make(const std::vector<float>& z, const std::vector<std::uint8_t>& m)
{
    lsdtt::lsdtt_raster_meta meta;
    meta.n_rows = N; meta.n_cols = N; meta.x_min = X0; meta.y_min = Y0; meta.cell_size = CS;
    meta.no_data_value = -9999.0f;
    return lsdtt::lsdtt_raster_cpu(std::vector<float>(z), std::vector<std::uint8_t>(m), meta);
}

std::vector<std::uint8_t> line_cells(const std::vector<std::vector<lsdtt_xy>>& lines)
{
    std::vector<std::uint8_t> b(N * N, 0);
    for (const auto& l : lines)
        for (size_t i = 0; i + 1 < l.size(); ++i)
        {
            int pr = -1, pc = -1;
            const int k = 4000;
            for (int s2 = 0; s2 <= k; ++s2)
            {
                const double t = static_cast<double>(s2) / k;
                const double x = l[i].x + t * (l[i + 1].x - l[i].x), y = l[i].y + t * (l[i + 1].y - l[i].y);
                const int c = std::min(N - 1, std::max(0, static_cast<int>((x - X0) / CS)));
                const int r = std::min(N - 1, std::max(0, static_cast<int>((YMAX - y) / CS)));
                b[static_cast<size_t>(r * N + c)] = 1;
                if (pr >= 0 && pr != r && pc != c) b[static_cast<size_t>(pr * N + c)] = 1;   // 4-connect a diagonal step
                pr = r; pc = c;
            }
        }
    return b;
}

std::vector<std::uint8_t> expected_east()
{
    std::vector<std::uint8_t> e(N * N, 0);
    for (int r = 0; r < N; ++r)
        for (int c = 0; c < N; ++c) e[static_cast<size_t>(r * N + c)] = east(centre(r, c)) > 0 ? 1 : 0;
    return e;
}
long count(const std::vector<std::uint8_t>& v) { long n = 0; for (auto x : v) n += x; return n; }
long diff(const std::vector<std::uint8_t>& a, const std::vector<std::uint8_t>& b)
{ long n = 0; for (size_t i = 0; i < a.size(); ++i) n += (a[i] != 0) != (b[i] != 0); return n; }
} // namespace

int main()
{
    using namespace lsdtt;
    std::vector<float> z(N * N);
    std::vector<std::uint8_t> m(N * N, 1);
    for (int r = 0; r < N; ++r)
        for (int c = 0; c < N; ++c)
            z[static_cast<size_t>(r * N + c)] = east(centre(r, c)) > 0 ? -1.0f : (in_polder(r, c) ? -3.0f : 10.0f);
    const auto dem = make(z, m);
    const auto east_cells = expected_east();
    const std::vector<std::vector<lsdtt_xy>> line_ns = {{A, B}}, line_sn = {{B, A}};
    const auto bar = line_cells(line_ns);
    sea_mask_options o;

    // ---- lines ----
    auto r1 = sea_mask_from_lines(dem, bar, line_ns, o);
    check("lines auto: the sea is exactly the cells east of the coast", diff(r1.sea, east_cells) == 0 && !r1.leaked,
          s(diff(r1.sea, east_cells)) + " cells differ; " + r1.evidence);
    bool polder_kept = true;
    for (int r = 0; r < N; ++r) for (int c = 0; c < N; ++c) if (in_polder(r, c) && r1.sea[static_cast<size_t>(r * N + c)]) polder_kept = false;
    check("lines auto: the polder (-3 m, below the sea) stays land", polder_kept);
    check("lines auto: sea on the left of the north-to-south line, by median elevation",
          r1.side_used == "sea_left" && r1.evidence.find("median") != std::string::npos, r1.side_used + "; " + r1.evidence);
    auto r2 = sea_mask_from_lines(dem, bar, line_sn, o);
    check("lines auto, digitised south to north: the same sea", diff(r2.sea, east_cells) == 0 && r2.side_used == "sea_right",
          r2.side_used);
    o.side = sea_side_choice::sea_left;
    check("sea_left on the north-to-south line: the east", diff(sea_mask_from_lines(dem, bar, line_ns, o).sea, east_cells) == 0);
    o.side = sea_side_choice::sea_right;
    const auto rw = sea_mask_from_lines(dem, bar, line_ns, o);
    long west_ok = 0;
    for (int r = 0; r < N; ++r) for (int c = 0; c < N; ++c)
    {
        const size_t i = static_cast<size_t>(r * N + c);
        if ((east(centre(r, c)) < 0) == (rw.sea[i] != 0)) ++west_ok;
    }
    check("sea_right on the north-to-south line: the west", west_ok == N * N, s(west_ok));
    o.side = sea_side_choice::automatic;

    // gap: the line broken between 40% and 55% of its length
    auto at = [](double t) { return lsdtt_xy{A.x + t * (B.x - A.x), A.y + t * (B.y - A.y)}; };
    const std::vector<std::vector<lsdtt_xy>> gapped = {{A, at(0.40)}, {at(0.55), B}};
    o.max_distance_metres = 50.0;
    const auto rg = sea_mask_from_lines(dem, line_cells(gapped), gapped, o);
    bool band_ok = true;
    for (int r = 0; r < N; ++r) for (int c = 0; c < N; ++c)
    {
        const size_t i = static_cast<size_t>(r * N + c);
        if (rg.sea[i] && !east_cells[i]) band_ok = false;   // never on the land side
    }
    check("a gap in the line: leak reported, fallback band used, never on the land side",
          rg.leaked && rg.fell_back && band_ok && !rg.notes.empty() && count(rg.sea) > 0 && count(rg.sea) < count(east_cells),
          s(count(rg.sea)));
    const std::vector<std::vector<lsdtt_xy>> half = {{A, at(0.5)}};
    const auto rh = sea_mask_from_lines(dem, line_cells(half), half, o);
    check("a line stopping half way: leak and fallback", rh.leaked && rh.fell_back);
    o.max_distance_metres = 500.0;

    // flat: no contrast
    const auto flat = make(std::vector<float>(N * N, 0.0f), m);
    bool threw = false;
    try { (void)sea_mask_from_lines(flat, bar, line_ns, o); } catch (const std::invalid_argument&) { threw = true; }
    check("a flat DEM: auto refuses", threw);

    // ---- polygons ----
    auto rp = sea_mask_from_polygons(dem, east_cells, o);
    check("polygons over the sea, auto: sea_inside, the east", rp.side_used == "sea_inside" && diff(rp.sea, east_cells) == 0,
          rp.side_used + "; " + rp.evidence);
    std::vector<std::uint8_t> west_cover(N * N);
    for (size_t i = 0; i < west_cover.size(); ++i) west_cover[i] = east_cells[i] ? 0 : 1;
    auto rl = sea_mask_from_polygons(dem, west_cover, o);
    check("polygons over the land, auto: land_inside, the east", rl.side_used == "land_inside" && diff(rl.sea, east_cells) == 0,
          rl.side_used);
    o.side = sea_side_choice::land_inside;
    check("land_inside given: the east", diff(sea_mask_from_polygons(dem, west_cover, o).sea, east_cells) == 0);
    o.side = sea_side_choice::automatic;
    // the sea already NoData
    std::vector<std::uint8_t> m2(m);
    for (size_t i = 0; i < m2.size(); ++i) if (east_cells[i]) m2[i] = 0;
    const auto clipped = make(z, m2);
    auto rn = sea_mask_from_polygons(clipped, east_cells, o);
    check("sea already NoData: auto decides by NoData, sea_inside, nothing valid removed",
          rn.side_used == "sea_inside" && rn.evidence.find("NoData") != std::string::npos && rn.n_sea_valid_cells == 0,
          rn.evidence);
    // wrong side family
    threw = false;
    o.side = sea_side_choice::sea_inside;
    try { (void)sea_mask_from_lines(dem, bar, line_ns, o); } catch (const std::invalid_argument&) { threw = true; }
    check("sea_inside with lines: refused", threw);
    threw = false;
    o.side = sea_side_choice::sea_left;
    try { (void)sea_mask_from_polygons(dem, east_cells, o); } catch (const std::invalid_argument&) { threw = true; }
    check("sea_left with polygons: refused", threw);
    threw = false;
    o.side = sea_side_choice::land_inside;
    try { (void)sea_mask_from_polygons(dem, std::vector<std::uint8_t>(N * N, 0), o); } catch (const std::invalid_argument&) { threw = true; }
    check("a side that would remove every cell: refused", threw);
    o.side = sea_side_choice::automatic;

    // ---- an enclosed lagoon is an outlet after masking ----
    {
        std::vector<float> bowl(N * N);
        std::vector<std::uint8_t> lagoon(N * N, 0);
        for (int r = 0; r < N; ++r)
            for (int c = 0; c < N; ++c)
            {
                const double d = std::hypot(r - 19.5, c - 19.5);
                bowl[static_cast<size_t>(r * N + c)] = static_cast<float>(5.0 + 0.5 * d);   // rises to the edges
                if (std::abs(r - 19.5) < 3 && std::abs(c - 19.5) < 3) { lagoon[static_cast<size_t>(r * N + c)] = 1; bowl[static_cast<size_t>(r * N + c)] = -0.5f; }
            }
        auto b = make(bowl, m);
        o.side = sea_side_choice::sea_inside;
        const auto rs = sea_mask_from_polygons(b, lagoon, o);
        const std::size_t n_removed = apply_sea_mask(b, rs.sea);
        const auto filled = fill_pits(b, 0.0001f, 1);
        long raised = 0;
        for (int r = 0; r < N; ++r)
            for (int c = 0; c < N; ++c)
            {
                const size_t i = static_cast<size_t>(r * N + c);
                if (!b.get_mask()[i]) continue;
                if (filled.get_data()[i] > b.get_data()[i] + 1e-3f) ++raised;
            }
        check("lagoon: 36 cells removed, and after fill_pits no land cell is raised (the masked sea is an outlet)",
              n_removed == 36 && raised == 0, s(static_cast<double>(n_removed)) + " removed, " + s(raised) + " raised");
    }

    std::cout << "\n" << (g_n - g_fail) << " / " << g_n << " checks passed\n";
    return g_fail == 0 ? 0 : 1;
}
