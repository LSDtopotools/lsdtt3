// test_valley_relative_elevation.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
//
// -----------------------------------------------------------------------------
// compute_centreline_baseline and compute_valley_relative_elevation
// (lsdtt_valley_relative_elevation_cpu) on hand-checkable geometry.
//
// Baselines: a straight centreline along row 20 and a channel that zigzags
// across it, crossing at known columns with known elevations, one of which
// rises downstream. channel_crossings must pass through every crossing's
// elevation (rise included), be linear between and flat beyond;
// projected_monotonic must never rise and never sit above the channel.
//
// Relative elevation: two parallel centrelines in one valley strip, the
// "larger" (higher priority) one must take every pixel its discs cover.
// -----------------------------------------------------------------------------

#include <cmath>
#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

#include "lsdtt_raster_cpu.hpp"
#include "lsdtt_valley_relative_elevation_cpu.hpp"

namespace
{
int g_failures = 0;
void check(bool ok, const std::string& what)
{
    std::cout << (ok ? "PASS  " : "FAIL  ") << what << std::endl;
    if (!ok) ++g_failures;
}

lsdtt::lsdtt_raster_meta grid(int rows, int cols)
{
    lsdtt::lsdtt_raster_meta m;
    m.n_rows = rows;
    m.n_cols = cols;
    m.cell_size = 1.0;
    m.no_data_value = -9999.0;
    m.set_default_synthetic_georeference();
    return m;
}
} // namespace

int main()
{
    using lsdtt::centreline_baseline_choice;

    // -------------------------------------------------------------------------
    // Baselines
    // -------------------------------------------------------------------------
    {
        const int R = 41, C = 61;
        const lsdtt::lsdtt_raster_meta meta = grid(R, C);
        std::vector<float> z(static_cast<std::size_t>(R) * C, 0.0f);

        // Centreline: row 20, cols 0..60, west -> east (downstream).
        std::vector<std::int64_t> cl;
        for (int c = 0; c <= 60; ++c) cl.push_back(20 * C + c);

        // Channel: vertical runs joined along rows 15 and 25, crossing row 20
        // at cols 10, 20, 30, 40 and 50. Each crossing pixel is given its own
        // elevation; the channel elsewhere is set so it descends to the next.
        const int xc[5] = {10, 20, 30, 40, 50};
        const float zc[5] = {10.0f, 8.0f, 9.0f, 6.0f, 5.0f};   // rises at col 30
        std::vector<std::int64_t> ch;
        int row = 15;
        ch.push_back(15 * C + 5);
        for (int c = 6; c <= xc[0]; ++c) ch.push_back(15 * C + c);
        for (int i = 0; i < 5; ++i)
        {
            const int to = (row == 15) ? 25 : 15;
            const int step = (to > row) ? 1 : -1;
            for (int r = row + step; r != to + step; r += step) ch.push_back(r * C + xc[i]);
            row = to;
            const int next = (i < 4) ? xc[i + 1] : 55;
            for (int c = xc[i] + 1; c <= next; ++c) ch.push_back(row * C + c);
        }
        // Elevations: each crossing pixel gets zc; others the crossing value of
        // the run they belong to (only the crossing pixels matter here).
        int k = -1;
        for (const std::int64_t f : ch)
        {
            const int c = static_cast<int>(f % C);
            for (int i = 0; i < 5; ++i) if (c == xc[i]) k = i;
            z[static_cast<std::size_t>(f)] = k < 0 ? 11.0f : zc[k];
        }

        const auto bc = lsdtt::compute_centreline_baseline(
            cl, ch, z, meta, centreline_baseline_choice::channel_crossings);
        bool at_crossings = true;
        for (int i = 0; i < 5; ++i)
            if (std::abs(bc[static_cast<std::size_t>(xc[i])] - zc[i]) > 1e-4f) at_crossings = false;
        check(at_crossings, "channel_crossings passes through every crossing's elevation");
        check(std::abs(bc[30] - 9.0f) < 1e-4f && bc[30] > bc[20],
              "channel_crossings follows the channel where it rises downstream");
        check(std::abs(bc[15] - 9.0f) < 1e-4f && std::abs(bc[35] - 7.5f) < 1e-4f,
              "channel_crossings is linear between crossings");
        check(std::abs(bc[0] - 10.0f) < 1e-4f && std::abs(bc[60] - 5.0f) < 1e-4f,
              "channel_crossings is flat beyond the first and last crossing");

        const auto bp = lsdtt::compute_centreline_baseline(
            cl, ch, z, meta, centreline_baseline_choice::projected_monotonic);
        bool never_rises = true, below_pins = true;
        for (std::size_t i = 1; i < bp.size(); ++i)
            if (bp[i] > bp[i - 1] + 1e-5f) never_rises = false;
        for (int i = 0; i < 5; ++i)
            if (bp[static_cast<std::size_t>(xc[i])] > zc[i] + 1e-4f) below_pins = false;
        check(never_rises, "projected_monotonic never rises downstream");
        check(below_pins, "projected_monotonic is never above a crossing's channel elevation");
        check(bp[30] <= 8.0f + 1e-4f,
              "projected_monotonic holds the baseline down through the rise");

        // A centreline running along the channel's own pixels (collinear
        // segments): every shared node pins, so the baseline is the channel.
        {
            std::vector<std::int64_t> same;
            std::vector<float> zz(z.size(), 0.0f);
            for (int c = 0; c <= 60; ++c)
            {
                same.push_back(20 * C + c);
                zz[static_cast<std::size_t>(20 * C + c)] =
                    static_cast<float>(50.0 - 0.5 * c + ((c % 7 == 3) ? 0.8 : 0.0));
            }
            const auto bs = lsdtt::compute_centreline_baseline(
                same, same, zz, meta, centreline_baseline_choice::channel_crossings);
            bool equal = true;
            for (int c = 0; c <= 60; ++c)
                if (std::abs(bs[static_cast<std::size_t>(c)] -
                             zz[static_cast<std::size_t>(20 * C + c)]) > 1e-4f) equal = false;
            check(equal, "channel_crossings on a centreline sharing the channel's pixels is the channel profile");
        }

        // No crossing at all: channel_crossings has nothing to pin.
        std::vector<std::int64_t> away;
        for (int c = 0; c <= 60; ++c) away.push_back(5 * C + c);
        const auto bn = lsdtt::compute_centreline_baseline(
            cl, away, z, meta, centreline_baseline_choice::channel_crossings);
        check(!std::isfinite(bn[0]), "channel_crossings with no crossing is NaN");
    }

    // -------------------------------------------------------------------------
    // Relative elevation: the larger channel wins
    // -------------------------------------------------------------------------
    {
        const int R = 31, C = 80;
        const lsdtt::lsdtt_raster_meta meta = grid(R, C);
        const std::size_t n = static_cast<std::size_t>(R) * C;
        std::vector<float> z(n), v(n, 0.0f);
        for (int r = 0; r < R; ++r)
            for (int c = 0; c < C; ++c)
            {
                z[static_cast<std::size_t>(r) * C + c] = static_cast<float>(100.0 + r * 0.1);
                if (r >= 5 && r <= 25) v[static_cast<std::size_t>(r) * C + c] = 1.0f;
            }
        const lsdtt::lsdtt_raster_cpu dem(std::vector<float>(z), std::vector<std::uint8_t>(n, 1), meta);
        const lsdtt::lsdtt_raster_cpu vm(std::move(v), std::vector<std::uint8_t>(n, 1), meta);

        // Two centrelines across the strip: big along row 15 (the valley's
        // middle, discs radius ~11), small along row 8 near the wall (discs
        // radius ~4), both baseline-constant but different.
        lsdtt::valley_rer_centreline big, small;
        big.id = 1;
        small.id = 2;
        for (int c = 0; c < C; ++c)
        {
            big.flats.push_back(15 * C + c);
            small.flats.push_back(8 * C + c);
        }
        big.baseline.assign(big.flats.size(), 90.0f);
        small.baseline.assign(small.flats.size(), 95.0f);
        big.priority.assign(big.flats.size(), 1e6);
        small.priority.assign(small.flats.size(), 1e3);

        const auto rr = lsdtt::compute_valley_relative_elevation(dem, vm, {small, big});
        const auto& own = rr.owner.get_data();
        const auto& rer = rr.relative_elevation.get_data();
        bool big_wins = true;
        for (int r = 5; r <= 25; ++r)
            for (int c = 0; c < C; ++c)
                if (own[static_cast<std::size_t>(r) * C + c] != 1.0f) big_wins = false;
        check(big_wins, "every pixel in the larger centreline's valley is measured against it");
        check(rr.n_unassigned == 0, "every valley pixel is assigned");
        check(std::abs(rer[8 * C + 40] - (100.8f - 90.0f)) < 1e-4f,
              "relative elevation = elevation - the larger centreline's baseline");
        check(rr.n_nodes_dropped == static_cast<std::int64_t>(small.flats.size()),
              "the smaller centreline's nodes inside the larger's valley are dropped");
        check(own[0] == static_cast<float>(meta.no_data_value),
              "nothing outside the valley mask");

        // Without the big one, the small one takes the valley.
        const auto r2 = lsdtt::compute_valley_relative_elevation(dem, vm, {small});
        check(r2.owner.get_data()[20 * C + 40] == 2.0f && r2.n_unassigned == 0,
              "a lone centreline takes the whole connected valley");
    }

    std::cout << (g_failures == 0 ? "ALL PASSED" : "FAILURES: " + std::to_string(g_failures))
              << std::endl;
    return g_failures == 0 ? 0 : 1;
}
