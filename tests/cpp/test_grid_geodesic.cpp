// test_grid_geodesic.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
//
// -----------------------------------------------------------------------------
// Characterisation test for grid_dijkstra_d8 (src/utils/lsdtt_grid_geodesic.hpp,
// consolidation C11).
//
// Three Dijkstra loops were replaced. Each is kept below verbatim, reduced to the
// inputs it reads, and compared with grid_dijkstra_d8 called the way the
// migrated site now calls it:
//
//   detour     lsdtt_raster_binary_mask_cpu 2b: float geodesic distance from
//              seeds through a mask
//   nearest    lsdtt_raster_binary_mask_cpu 2c: float distance plus the source
//              seed carried to every cell
//   centreline lsdtt_valley_centreline_cpu: double cost with a per-cell weight,
//              predecessors, and an early stop at the target
//
// Every distance is compared by bit pattern; sources and predecessors exactly.
// Grids are built with many exact ties (unit cell sizes, integer weights,
// several seeds), because tie order is what a careless rewrite changes.
//
// Usage: lsdtt-test-grid-geodesic
// Prints one line per check and a final "RESULT: PASS|FAIL".
// -----------------------------------------------------------------------------

#include <cmath>
#include <cstdint>
#include <cstring>
#include <functional>
#include <iostream>
#include <limits>
#include <queue>
#include <random>
#include <string>
#include <utility>
#include <vector>

#include "../../src/utils/lsdtt_grid_geodesic.hpp"
#include "../../src/utils/lsdtt_grid_neighbours.hpp"

using namespace lsdtt;

namespace
{

static constexpr const auto& k_dr = grid::d8_rowmajor_row;
static constexpr const auto& k_dc = grid::d8_rowmajor_col;

// ============================================================================
// Removed copies (references)
// ============================================================================

// src/cpu/lsdtt_raster_binary_mask_cpu.cpp:1451-1488 (as of 2026-10-07), 2b
std::vector<float> old_detour(const std::vector<std::uint8_t>& mask, const std::vector<std::int64_t>& seed_flats,
                              int n_rows, int n_cols, double cell)
{
    const std::size_t n = mask.size();
    constexpr float k_unreached = std::numeric_limits<float>::infinity();
    std::vector<float> geo(n, k_unreached);
    using qnode = std::pair<float, std::int64_t>;
    std::priority_queue<qnode, std::vector<qnode>,
                        std::greater<qnode>> pq;
    for (const std::int64_t cf : seed_flats)
    {
        if (cf < 0 || cf >= static_cast<std::int64_t>(n)) continue;
        const std::size_t fi = static_cast<std::size_t>(cf);
        if (!mask[fi]) continue;
        geo[fi] = 0.0f;
        pq.emplace(0.0f, cf);
    }
    const float straight = static_cast<float>(cell);
    const float diag = static_cast<float>(cell * grid::k_sqrt2_d);
    while (!pq.empty())
    {
        const auto [d, flat] = pq.top();
        pq.pop();
        const std::size_t fi = static_cast<std::size_t>(flat);
        if (d > geo[fi]) continue;
        const int r = static_cast<int>(flat / n_cols);
        const int c = static_cast<int>(flat % n_cols);
        for (int k = 0; k < 8; ++k)
        {
            const int rr = r + k_dr[k], cc = c + k_dc[k];
            if (rr < 0 || rr >= n_rows || cc < 0 || cc >= n_cols)
                continue;
            const std::size_t nf =
                static_cast<std::size_t>(rr) * n_cols + cc;
            if (!mask[nf]) continue;
            const float nd =
                d + ((k_dr[k] != 0 && k_dc[k] != 0) ? diag : straight);
            if (nd < geo[nf])
            {
                geo[nf] = nd;
                pq.emplace(nd, static_cast<std::int64_t>(nf));
            }
        }
    }
    return geo;
}

// src/cpu/lsdtt_raster_binary_mask_cpu.cpp:1519-1557 (as of 2026-10-07), 2c
void old_nearest(const std::vector<std::uint8_t>& distance_mask, const std::vector<std::int64_t>& seed_flats,
                 int n_rows, int n_cols, double cell,
                 std::vector<float>& gdist, std::vector<std::int64_t>& gsrc)
{
    const std::size_t n = distance_mask.size();
    gdist.assign(n, std::numeric_limits<float>::infinity());
    gsrc.assign(n, -1);
    using qnode = std::pair<float, std::int64_t>;
    std::priority_queue<qnode, std::vector<qnode>,
                        std::greater<qnode>> pq;
    for (const std::int64_t cf : seed_flats)
    {
        if (cf < 0 || cf >= static_cast<std::int64_t>(n)) continue;
        gdist[static_cast<std::size_t>(cf)] = 0.0f;
        gsrc[static_cast<std::size_t>(cf)]  = cf;
        pq.emplace(0.0f, cf);
    }
    const float straight = static_cast<float>(cell);
    const float diag = static_cast<float>(cell * grid::k_sqrt2_d);
    while (!pq.empty())
    {
        const auto [d, flat] = pq.top();
        pq.pop();
        const std::size_t fi = static_cast<std::size_t>(flat);
        if (d > gdist[fi]) continue;
        const int r = static_cast<int>(flat / n_cols);
        const int c = static_cast<int>(flat % n_cols);
        for (int k = 0; k < 8; ++k)
        {
            const int rr = r + k_dr[k], cc = c + k_dc[k];
            if (rr < 0 || rr >= n_rows || cc < 0 || cc >= n_cols)
                continue;
            const std::size_t nf =
                static_cast<std::size_t>(rr) * n_cols + cc;
            if (!distance_mask[nf]) continue;
            const float nd =
                d + ((k_dr[k] != 0 && k_dc[k] != 0) ? diag : straight);
            if (nd < gdist[nf])
            {
                gdist[nf] = nd;
                gsrc[nf]  = gsrc[fi];
                pq.emplace(nd, static_cast<std::int64_t>(nf));
            }
        }
    }
}

// src/cpu/lsdtt_valley_centreline_cpu.cpp:479-518 (as of 2026-10-07)
void old_centreline(const std::vector<std::uint8_t>& hole, const std::vector<float>& wd, float wmax,
                    std::int64_t start_f, std::int64_t target_f, int n_rows, int n_cols, double cell,
                    std::vector<double>& cost, std::vector<std::int64_t>& prev)
{
    const std::size_t n = hole.size();
    const double eps = 0.5;
    cost.assign(n, std::numeric_limits<double>::infinity());
    prev.assign(n, -1);
    using qnode = std::pair<double, std::int64_t>;
    std::priority_queue<qnode, std::vector<qnode>,
                        std::greater<qnode>> pq;
    cost[static_cast<std::size_t>(start_f)] = 0.0;
    pq.emplace(0.0, start_f);
    const double straight = cell;
    const double diag = cell * grid::k_sqrt2_d;
    while (!pq.empty())
    {
        const auto [d, flat] = pq.top();
        pq.pop();
        const std::size_t fi = static_cast<std::size_t>(flat);
        if (d > cost[fi]) continue;
        if (flat == target_f) break;
        const int r = static_cast<int>(flat / n_cols);
        const int c = static_cast<int>(flat % n_cols);
        for (int k = 0; k < 8; ++k)
        {
            const int rr = r + k_dr[k], cc = c + k_dc[k];
            if (rr < 0 || rr >= n_rows || cc < 0 || cc >= n_cols)
                continue;
            const std::size_t nf =
                static_cast<std::size_t>(rr) * n_cols + cc;
            if (!hole[nf]) continue;
            const double step =
                (k_dr[k] != 0 && k_dc[k] != 0) ? diag : straight;
            const double w = step *
                (static_cast<double>(wmax - wd[nf]) + eps);
            const double nd = d + w;
            if (nd < cost[nf])
            {
                cost[nf] = nd;
                prev[nf] = flat;
                pq.emplace(nd, static_cast<std::int64_t>(nf));
            }
        }
    }
}

// ============================================================================
// The migrated sites, written exactly as they now appear in the sources
// ============================================================================

std::vector<float> new_detour(const std::vector<std::uint8_t>& mask, const std::vector<std::int64_t>& seed_flats,
                              int n_rows, int n_cols, double cell)
{
    const std::size_t n = mask.size();
    std::vector<float> geo(n, std::numeric_limits<float>::infinity());
    std::vector<std::int64_t> sources;
    for (const std::int64_t cf : seed_flats)
    {
        if (cf < 0 || cf >= static_cast<std::int64_t>(n)) continue;
        const std::size_t fi = static_cast<std::size_t>(cf);
        if (!mask[fi]) continue;
        geo[fi] = 0.0f;
        sources.push_back(cf);
    }
    const float straight = static_cast<float>(cell);
    const float diag = static_cast<float>(cell * grid::k_sqrt2_d);
    grid_dijkstra_d8(n_rows, n_cols, geo, sources,
        [&](std::size_t nf) { return mask[nf] != 0; },
        [&](std::size_t, std::size_t, int k) {
            return (k_dr[k] != 0 && k_dc[k] != 0) ? diag : straight; },
        [](std::size_t, std::size_t) {},
        [](std::int64_t) { return false; });
    return geo;
}

void new_nearest(const std::vector<std::uint8_t>& distance_mask, const std::vector<std::int64_t>& seed_flats,
                 int n_rows, int n_cols, double cell,
                 std::vector<float>& gdist, std::vector<std::int64_t>& gsrc)
{
    const std::size_t n = distance_mask.size();
    gdist.assign(n, std::numeric_limits<float>::infinity());
    gsrc.assign(n, -1);
    std::vector<std::int64_t> sources;
    for (const std::int64_t cf : seed_flats)
    {
        if (cf < 0 || cf >= static_cast<std::int64_t>(n)) continue;
        gdist[static_cast<std::size_t>(cf)] = 0.0f;
        gsrc[static_cast<std::size_t>(cf)]  = cf;
        sources.push_back(cf);
    }
    const float straight = static_cast<float>(cell);
    const float diag = static_cast<float>(cell * grid::k_sqrt2_d);
    grid_dijkstra_d8(n_rows, n_cols, gdist, sources,
        [&](std::size_t nf) { return distance_mask[nf] != 0; },
        [&](std::size_t, std::size_t, int k) {
            return (k_dr[k] != 0 && k_dc[k] != 0) ? diag : straight; },
        [&](std::size_t from, std::size_t to) { gsrc[to] = gsrc[from]; },
        [](std::int64_t) { return false; });
}

void new_centreline(const std::vector<std::uint8_t>& hole, const std::vector<float>& wd, float wmax,
                    std::int64_t start_f, std::int64_t target_f, int n_rows, int n_cols, double cell,
                    std::vector<double>& cost, std::vector<std::int64_t>& prev)
{
    const std::size_t n = hole.size();
    const double eps = 0.5;
    cost.assign(n, std::numeric_limits<double>::infinity());
    prev.assign(n, -1);
    cost[static_cast<std::size_t>(start_f)] = 0.0;
    const double straight = cell;
    const double diag = cell * grid::k_sqrt2_d;
    grid_dijkstra_d8(n_rows, n_cols, cost, std::vector<std::int64_t>{start_f},
        [&](std::size_t nf) { return hole[nf] != 0; },
        [&](std::size_t, std::size_t nf, int k) {
            const double step = (k_dr[k] != 0 && k_dc[k] != 0) ? diag : straight;
            const double w = step * (static_cast<double>(wmax - wd[nf]) + eps);
            return w; },
        [&](std::size_t from, std::size_t to) { prev[to] = static_cast<std::int64_t>(from); },
        [&](std::int64_t flat) { return flat == target_f; });
}

// ============================================================================
// Harness
// ============================================================================

int g_failures = 0;

void check(const std::string& what, bool ok, long long n_cases, const std::string& detail = "")
{
    std::cout << (ok ? "PASS " : "FAIL ") << what << " (" << n_cases << " cases)";
    if (!ok && !detail.empty()) std::cout << " -- " << detail;
    std::cout << "\n";
    if (!ok) ++g_failures;
}

template <class T>
bool same_bits(const std::vector<T>& a, const std::vector<T>& b)
{
    return a.size() == b.size() && (a.empty() || std::memcmp(a.data(), b.data(), a.size() * sizeof(T)) == 0);
}

} // namespace

int main()
{
    std::cout << "== C11 grid Dijkstra vs the removed loops ==\n";
    std::mt19937 rng(42);
    std::uniform_real_distribution<double> u(0.0, 1.0);
    long long n_cases = 0, bad_detour = 0, bad_nearest = 0, bad_centre = 0;

    const std::pair<int, int> dims[] = {{1, 1}, {1, 12}, {9, 1}, {6, 6}, {20, 33}, {47, 29}};
    for (auto [nr, nc] : dims)
        for (double density : {0.4, 0.7, 1.0})
            for (double cell : {1.0, 2.0, 30.0, 0.3})
                for (int rep = 0; rep < 4; ++rep)
                {
                    const std::size_t n = static_cast<std::size_t>(nr) * nc;
                    std::vector<std::uint8_t> mask(n);
                    std::vector<float> wd(n);
                    for (std::size_t i = 0; i < n; ++i)
                    {
                        mask[i] = (u(rng) < density) ? 1 : 0;
                        wd[i] = static_cast<float>(std::floor(u(rng) * 4.0));   // integer weights: ties
                    }
                    std::vector<std::int64_t> seeds;
                    const int n_seeds = 1 + rep;
                    for (int s = 0; s < n_seeds; ++s)
                        seeds.push_back(static_cast<std::int64_t>(u(rng) * n));
                    seeds.push_back(-1);                                   // out of range: skipped
                    seeds.push_back(static_cast<std::int64_t>(n));
                    if (!seeds.empty()) seeds.push_back(seeds.front());   // duplicate seed
                    ++n_cases;

                    if (!same_bits(old_detour(mask, seeds, nr, nc, cell), new_detour(mask, seeds, nr, nc, cell)))
                        ++bad_detour;

                    std::vector<float> d1, d2;
                    std::vector<std::int64_t> s1, s2;
                    old_nearest(mask, seeds, nr, nc, cell, d1, s1);
                    new_nearest(mask, seeds, nr, nc, cell, d2, s2);
                    if (!same_bits(d1, d2) || s1 != s2) ++bad_nearest;

                    // Centreline: start and target inside the hole when possible.
                    std::vector<std::uint8_t> hole = mask;
                    const std::int64_t start_f = static_cast<std::int64_t>(u(rng) * n);
                    const std::int64_t target_f = static_cast<std::int64_t>(u(rng) * n);
                    hole[static_cast<std::size_t>(start_f)] = 1;
                    hole[static_cast<std::size_t>(target_f)] = 1;
                    float wmax = 0.0f;
                    for (std::size_t i = 0; i < n; ++i) if (hole[i]) wmax = std::max(wmax, wd[i]);
                    std::vector<double> c1, c2;
                    std::vector<std::int64_t> p1, p2;
                    old_centreline(hole, wd, wmax, start_f, target_f, nr, nc, cell, c1, p1);
                    new_centreline(hole, wd, wmax, start_f, target_f, nr, nc, cell, c2, p2);
                    if (!same_bits(c1, c2) || p1 != p2) ++bad_centre;
                }

    check("binary_mask detour distances (float bits)", bad_detour == 0, n_cases);
    check("binary_mask nearest-seed distances (float bits) and sources", bad_nearest == 0, n_cases);
    check("valley_centreline costs (double bits), predecessors, early stop", bad_centre == 0, n_cases);

    std::cout << "RESULT: " << (g_failures == 0 ? "PASS" : "FAIL") << "\n";
    return g_failures == 0 ? 0 : 1;
}
