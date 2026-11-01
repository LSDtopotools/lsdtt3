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

#include <algorithm>
#include <cmath>
#include <limits>
#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

#include "lsdtt_raster_cpu.hpp"
#include "lsdtt_valley_relative_elevation_cpu.hpp"
#include "lsdtt_valley_terrace_profile_cpu.hpp"
#include "../../src/utils/lsdtt_point_bucket_index.hpp"

#include <random>

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

        const auto rr = lsdtt::compute_valley_relative_elevation(dem, vm, vm, {small, big});
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
        const auto r2 = lsdtt::compute_valley_relative_elevation(dem, vm, vm, {small});
        check(r2.owner.get_data()[20 * C + 40] == 2.0f && r2.n_unassigned == 0,
              "a lone centreline takes the whole connected valley");
    }

    // -------------------------------------------------------------------------
    // A domain wider than the floor: terraces either side of a narrow floor,
    // a gully stream with no floor of its own crossing the terrace. The main
    // centreline's reach is scaled by its floor half-width, the gully's by one
    // cell, so the gully claims only the ground beside it.
    // -------------------------------------------------------------------------
    {
        const int R = 31, C = 80;
        const lsdtt::lsdtt_raster_meta meta = grid(R, C);
        const std::size_t n = static_cast<std::size_t>(R) * C;
        std::vector<float> z(n), fl(n, 0.0f), dm(n, 0.0f);
        for (int r = 0; r < R; ++r)
            for (int c = 0; c < C; ++c)
            {
                const std::size_t i = static_cast<std::size_t>(r) * C + c;
                z[i] = static_cast<float>(100.0 + std::abs(r - 15) * 1.0);
                if (r >= 11 && r <= 19) fl[i] = 1.0f;       // floor, half-width ~5
                if (r >= 2 && r <= 28) dm[i] = 1.0f;        // terraces
            }
        const lsdtt::lsdtt_raster_cpu dem(std::move(z), std::vector<std::uint8_t>(n, 1), meta);
        const lsdtt::lsdtt_raster_cpu floor(std::move(fl), std::vector<std::uint8_t>(n, 1), meta);
        const lsdtt::lsdtt_raster_cpu dom(std::move(dm), std::vector<std::uint8_t>(n, 1), meta);
        lsdtt::valley_rer_centreline main_c, gully;
        main_c.id = 1;
        gully.id = 2;
        for (int c = 0; c < C; ++c) main_c.flats.push_back(15 * C + c);
        for (int r = 2; r <= 10; ++r) gully.flats.push_back(r * C + 60);
        main_c.baseline.assign(main_c.flats.size(), 100.0f);
        gully.baseline.assign(gully.flats.size(), 105.0f);
        main_c.priority.assign(main_c.flats.size(), 1e6);
        gully.priority.assign(gully.flats.size(), 1e3);

        const auto rr = lsdtt::compute_valley_relative_elevation(dem, floor, dom, {main_c, gully});
        const auto& own = rr.owner.get_data();
        std::int64_t n_gully = 0, n_dom = 0;
        for (int r = 2; r <= 28; ++r)
            for (int c = 0; c < C; ++c)
            {
                ++n_dom;
                if (own[static_cast<std::size_t>(r) * C + c] == 2.0f) ++n_gully;
            }
        check(rr.n_unassigned == 0, "wide domain: every domain pixel is assigned");
        check(own[5 * C + 30] == 1.0f && own[25 * C + 60] == 1.0f,
              "wide domain: terrace away from the gully is measured against the main centreline");
        check(own[5 * C + 61] == 2.0f, "wide domain: the cell beside the gully is the gully's");
        check(n_gully < n_dom / 20,
              "wide domain: the floorless gully claims only the ground beside it (" +
                  std::to_string(n_gully) + " of " + std::to_string(n_dom) + ")");
        check(std::abs(rr.relative_elevation.get_data()[5 * C + 30] - 10.0f) < 1e-4f,
              "wide domain: terrace relative elevation is against the main baseline");
        check(own[0] == static_cast<float>(meta.no_data_value), "wide domain: nothing outside it");

        // A 4 m cap: terrace cells more than 4 m (ground) from every node are out.
        lsdtt::valley_rer_request req;
        req.max_distance_metres = 4.0;
        const auto rc = lsdtt::compute_valley_relative_elevation(dem, floor, dom, {main_c, gully}, req);
        check(rc.owner.get_data()[25 * C + 30] == static_cast<float>(meta.no_data_value) &&
                  rc.owner.get_data()[18 * C + 30] == 1.0f,
              "distance cap: 10 m away is left out, 3 m away is kept");

        // ---------------------------------------------------------------------
        // Terrace profile: patch 1 north (rows 3-6, cols 10-29), patch 2 south
        // (rows 24-27, cols 40-69), patch 3 straddling the gully (rows 3-6,
        // cols 50-69: the gully claims the few cells beside it, the main
        // centreline the rest). Relative elevation = |row - 15| against the main
        // centreline; 5 m bins along it (cell = 1 m).
        // ---------------------------------------------------------------------
        std::vector<float> lab(n, 0.0f);
        for (int r = 3; r <= 6; ++r)
            for (int c = 10; c <= 29; ++c) lab[static_cast<std::size_t>(r) * C + c] = 1.0f;
        for (int r = 24; r <= 27; ++r)
            for (int c = 40; c <= 69; ++c) lab[static_cast<std::size_t>(r) * C + c] = 2.0f;
        for (int r = 3; r <= 6; ++r)
            for (int c = 50; c <= 69; ++c) lab[static_cast<std::size_t>(r) * C + c] = 3.0f;
        const lsdtt::lsdtt_raster_cpu labels(std::move(lab), std::vector<std::uint8_t>(n, 1), meta);
        lsdtt::terrace_profile_request treq;
        treq.bin_metres = 5.0;
        const auto tp = lsdtt::compute_terrace_profile(dem, labels, {main_c, gully}, rr, treq);
        check(tp.patches.size() == 3, "terrace profile: three patches");
        int n1 = 0, n2 = 0;
        bool sizes = true, med1 = true;
        for (const auto& b : tp.bins)
        {
            if (b.patch_label == 1)
            {
                ++n1;
                if (b.n_pixels != 20) sizes = false;
                if (std::abs(b.relative_elevation[2] - 10.5f) > 1e-4f ||
                    std::abs(b.elevation[2] - 110.5f) > 1e-4f ||
                    std::abs(b.relative_elevation[0] - 9.0f) > 1e-4f ||
                    std::abs(b.relative_elevation[4] - 12.0f) > 1e-4f) med1 = false;
            }
            if (b.patch_label == 2) ++n2;
        }
        check(n1 == 4 && n2 == 6, "terrace profile: patches cut into 5 m bins along the centreline (" +
                                      std::to_string(n1) + ", " + std::to_string(n2) + ")");
        check(sizes, "terrace profile: 20 pixels in each of patch 1's bins");
        check(med1, "terrace profile: patch 1 bins span relative elevation 9-12 m, median 10.5 m");
        bool p3 = false;
        for (const auto& q : tp.patches)
            if (q.patch_label == 3)
                p3 = q.centreline_id == 1 && q.n_pixels_reassigned > 0 &&
                     q.n_pixels_reassigned < q.n_pixels / 2;
        check(p3, "terrace profile: a patch straddling the gully goes to the centreline owning "
                  "most of it, its gully pixels tied to that centreline");
        bool along = true;
        for (const auto& q : tp.patches)
            if (q.patch_label == 1 && (q.along_min_metres != 10.0 || q.along_max_metres != 29.0))
                along = false;
        check(along, "terrace profile: along-centreline extent of patch 1 is 10-29 m");
    }

    // -------------------------------------------------------------------------
    // point_bucket_index against brute force: nearest (including queries far
    // outside the points' extent) and within-radius, on a random cloud and on
    // a polyline-like run of points.
    // -------------------------------------------------------------------------
    {
        std::mt19937 rng(7);
        std::uniform_real_distribution<double> u(0.0, 1000.0), q(-500.0, 1500.0);
        std::vector<double> px, py;
        for (int i = 0; i < 2000; ++i) { px.push_back(u(rng)); py.push_back(u(rng)); }
        for (int i = 0; i < 300; ++i) { px.push_back(200.0 + i * 1.5); py.push_back(500.0); }
        const lsdtt::point_bucket_index idx(px, py, 20.0);
        bool nearest_ok = true, within_ok = true;
        for (int t = 0; t < 3000; ++t)
        {
            const double x = q(rng), y = q(rng);
            std::size_t bb = 0;
            double bd = std::numeric_limits<double>::infinity();
            for (std::size_t i = 0; i < px.size(); ++i)
            {
                const double d = (px[i] - x) * (px[i] - x) + (py[i] - y) * (py[i] - y);
                if (d < bd) { bd = d; bb = i; }
            }
            if (idx.nearest(x, y) != bb) nearest_ok = false;
            if (t < 300)
            {
                std::vector<std::size_t> got;
                idx.within(x, y, 37.0, got);
                std::sort(got.begin(), got.end());
                std::vector<std::size_t> want;
                for (std::size_t i = 0; i < px.size(); ++i)
                    if ((px[i] - x) * (px[i] - x) + (py[i] - y) * (py[i] - y) <= 37.0 * 37.0)
                        want.push_back(i);
                if (got != want) within_ok = false;
            }
        }
        check(nearest_ok, "point_bucket_index nearest matches brute force (3000 queries)");
        check(within_ok, "point_bucket_index within matches brute force (300 queries)");
    }

    std::cout << (g_failures == 0 ? "ALL PASSED" : "FAILURES: " + std::to_string(g_failures))
              << std::endl;
    return g_failures == 0 ? 0 : 1;
}
