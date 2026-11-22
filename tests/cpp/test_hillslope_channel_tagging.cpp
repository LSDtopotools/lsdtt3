// test_hillslope_channel_tagging.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
//
// -----------------------------------------------------------------------------
// Unit test for src/cpu/lsdtt_hillslope_channel_tagging_cpu:
//
//   * planar V-valleys, symmetric and asymmetric: exact HAND, ridge relief and
//     linear relative height, for the channel_pixel and entry_pixel groups;
//   * hand-built D8 paths: the farthest and the highest pixel, ties, a relative
//     height above 1, a channel fragment stopping the search, pixels reaching no
//     channel, zero relief;
//   * a rough routed surface: agreement with the existing HAND and
//     farthest-ridge routines, and the bounds and monotonicity of each method
//     against a brute-force search of each group;
//   * refusals.
//
// Usage: lsdtt-test-hillslope-channel-tagging
// Prints one line per check and a final "RESULT: PASS|FAIL".
// -----------------------------------------------------------------------------

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <iostream>
#include <map>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

#include "../../src/core/lsdtt_raster_meta.hpp"
#include "../../src/cpu/lsdtt_flowinfo_cpu.hpp"
#include "../../src/cpu/lsdtt_hillslope_channel_tagging_cpu.hpp"
#include "../../src/cpu/lsdtt_raster_cpu.hpp"
#include "../../src/cpu/lsdtt_raster_flow_preprocessing_cpu.hpp"

using namespace lsdtt;

namespace
{

int g_failures = 0;
int g_checks = 0;

void check(const std::string& name, bool ok, const std::string& detail = "")
{
    ++g_checks;
    if (!ok) ++g_failures;
    std::cout << (ok ? "  PASS  " : "  FAIL  ") << name;
    if (!detail.empty()) std::cout << "  (" << detail << ")";
    std::cout << "\n";
}

std::string num(double v)
{
    char buf[64];
    std::snprintf(buf, sizeof(buf), "%.6g", v);
    return buf;
}

lsdtt_raster_meta make_meta(int n_rows, int n_cols, double cs)
{
    lsdtt_raster_meta meta;
    meta.n_rows = n_rows;
    meta.n_cols = n_cols;
    meta.cell_size = cs;
    meta.x_min = 500000.0;
    meta.y_min = 4000000.0;
    meta.no_data_value = -9999.0;
    return meta;
}

lsdtt_flowinfo make_flow(const std::vector<float>& z, const std::vector<std::uint8_t>& mask,
                         const lsdtt_raster_meta& meta, const std::string& bc,
                         const std::vector<std::int64_t>& forced = {})
{
    const lsdtt_raster_cpu raster(std::vector<float>(z), std::vector<std::uint8_t>(mask), meta);
    lsdtt_flowinfo_options options;
    options.boundary_conditions = bc;
    options.build_topological_order = true;
    options.build_s_vector = true;
    options.interior_nodata_is_outlet = forced.empty();
    options.forced_outlet_flats = forced;
    return lsdtt_flowinfo(raster.get_data_array(), raster.get_mask_array(), meta, options);
}

// A valley running south, channel in column c0, planar sides of gradient S.
struct valley
{
    lsdtt_raster_meta meta;
    std::vector<float> z;
    std::vector<std::uint8_t> mask, channel;
};

valley make_valley(int nr, int nc, int c0, double S, double dx)
{
    valley v{make_meta(nr, nc, dx), std::vector<float>(nr * nc),
             std::vector<std::uint8_t>(nr * nc, 1u), std::vector<std::uint8_t>(nr * nc, 0u)};
    for (int r = 0; r < nr; ++r)
        for (int c = 0; c < nc; ++c)
        {
            v.z[r * nc + c] = static_cast<float>(1000.0 + S * dx * std::abs(c - c0) +
                                                 0.01 * dx * (nr - 1 - r));
            if (c == c0) v.channel[r * nc + c] = 1u;
        }
    return v;
}

// A single valid row between two nodata rows; only the forced outlets drain.
struct line_grid
{
    lsdtt_raster_meta meta;
    std::vector<float> z;
    std::vector<std::uint8_t> mask, channel;
    int nc;
    std::int64_t flat(int c) const { return nc + c; }
    std::size_t at(int c) const { return static_cast<std::size_t>(nc + c); }
};

line_grid make_line(const std::vector<float>& row, const std::vector<int>& channels, double dx)
{
    const int nc = static_cast<int>(row.size());
    line_grid g{make_meta(3, nc, dx), std::vector<float>(3 * nc, -9999.f),
                std::vector<std::uint8_t>(3 * nc, 0u), std::vector<std::uint8_t>(3 * nc, 0u), nc};
    for (int c = 0; c < nc; ++c) { g.z[nc + c] = row[c]; g.mask[nc + c] = 1u; }
    for (const int c : channels) g.channel[nc + c] = 1u;
    return g;
}

hillslope_relief_result relief(const lsdtt_flowinfo& flow, const std::vector<float>& z,
                               const std::vector<std::uint8_t>& mask,
                               const std::vector<std::uint8_t>& channel,
                               const std::string& ridge, const std::string& grouping)
{
    hillslope_relief_options o;
    o.ridge = ridge;
    o.grouping = grouping;
    return compute_hillslope_relief(flow, z, mask, channel, o);
}

void test_valleys()
{
    std::cout << "Planar V-valleys\n";
    const int nr = 20, nc = 11;
    const double S = 0.5, dx = 10.0;

    {
        const valley v = make_valley(nr, nc, 5, S, dx);
        const lsdtt_flowinfo flow = make_flow(v.z, v.mask, v.meta, "nnon");
        const auto res = relief(flow, v.z, v.mask, v.channel, "farthest", "channel_pixel");
        double err_h = 0, err_r = 0, err_f = 0, err_l = 0;
        for (int r = 0; r < nr - 1; ++r)
            for (int c = 0; c < nc; ++c)
            {
                if (c == 5) continue;
                const std::size_t i = r * nc + c;
                const int d = std::abs(c - 5);
                err_h = std::max(err_h, std::abs(res.hand.data[i] - S * dx * d));
                err_r = std::max(err_r, std::abs(res.ridge_relief.data[i] - S * dx * 5));
                err_f = std::max(err_f, std::abs(res.relative_height.data[i] - d / 5.0));
                err_l = std::max(err_l, std::abs(res.channel_distance.data[i] - dx * d));
            }
        check("symmetric valley: HAND", err_h < 1e-3, "max error " + num(err_h));
        check("symmetric valley: ridge relief", err_r < 1e-3, "max error " + num(err_r));
        check("symmetric valley: relative height is linear", err_f < 1e-4, "max error " + num(err_f));
        check("symmetric valley: D8 distance to channel", err_l < 1e-4, "max error " + num(err_l));
        check("symmetric valley: equal ridges resolve to the lower flat index",
              res.ridge_flat[3 * nc + 2] == 3 * nc + 0 && res.ridge_flat[3 * nc + 8] == 3 * nc + 0);
    }

    {
        const valley v = make_valley(nr, nc, 3, S, dx);
        const lsdtt_flowinfo flow = make_flow(v.z, v.mask, v.meta, "nnon");
        const auto shared = relief(flow, v.z, v.mask, v.channel, "farthest", "channel_pixel");
        const auto banks  = relief(flow, v.z, v.mask, v.channel, "farthest", "entry_pixel");
        double err_shared = 0, err_banks = 0;
        for (int r = 0; r < nr - 1; ++r)
            for (int c = 0; c < nc; ++c)
            {
                if (c == 3) continue;
                const std::size_t i = r * nc + c;
                const int d = std::abs(c - 3);
                err_shared = std::max(err_shared, std::abs(shared.relative_height.data[i] - d / 7.0));
                err_banks = std::max(err_banks,
                                     std::abs(banks.relative_height.data[i] - d / (c < 3 ? 3.0 : 7.0)));
            }
        check("asymmetric valley: channel_pixel shares the longer bank's ridge", err_shared < 1e-4,
              "max error " + num(err_shared));
        check("asymmetric valley: entry_pixel keeps each bank's ridge", err_banks < 1e-4,
              "max error " + num(err_banks));
    }
}

void test_paths()
{
    std::cout << "Hand-built D8 paths\n";
    const double dx = 10.0;

    {
        // A long gentle branch west of the channel and a short steep one east.
        const line_grid g = make_line({6, 5, 4, 3, 2, 1, 0, 50, 100}, {6}, dx);
        const lsdtt_flowinfo flow = make_flow(g.z, g.mask, g.meta, "nnnn", {g.flat(6)});

        const auto far = relief(flow, g.z, g.mask, g.channel, "farthest", "channel_pixel");
        check("farthest is the longest path, not the highest pixel",
              far.ridge_flat[g.at(3)] == g.flat(0) && far.ridge_flat[g.at(8)] == g.flat(0));
        check("farthest: ridge distance", std::abs(far.ridge_distance.data[g.at(8)] - 60.0) < 1e-4);
        check("farthest: a pixel above the ridge keeps its ratio above 1",
              std::abs(far.relative_height.data[g.at(8)] - 100.0 / 6.0) < 1e-4);
        check("farthest: pixels above the ridge are counted", far.n_above_ridge == 2,
              "n_above_ridge " + std::to_string(far.n_above_ridge));

        const auto high = relief(flow, g.z, g.mask, g.channel, "highest", "channel_pixel");
        check("highest is the highest pixel, not the longest path",
              high.ridge_flat[g.at(3)] == g.flat(8) && high.ridge_flat[g.at(0)] == g.flat(8));
        check("highest: nothing above the ridge", high.n_above_ridge == 0);
        check("highest: relative height", std::abs(high.relative_height.data[g.at(7)] - 0.5) < 1e-6);

        const auto own = relief(flow, g.z, g.mask, g.channel, "farthest", "pixel");
        check("pixel grouping: own upslope path", own.ridge_flat[g.at(3)] == g.flat(0) &&
              own.ridge_flat[g.at(7)] == g.flat(8));
        check("pixel grouping: a headwater is its own ridge",
              own.ridge_flat[g.at(8)] == g.flat(8) &&
              std::abs(own.relative_height.data[g.at(8)] - 1.0) < 1e-6);

        const auto entry = relief(flow, g.z, g.mask, g.channel, "farthest", "entry_pixel");
        check("entry_pixel: each branch has its own ridge",
              entry.ridge_flat[g.at(3)] == g.flat(0) && entry.ridge_flat[g.at(7)] == g.flat(8) &&
              entry.entry_flat[g.at(0)] == g.flat(5) && entry.entry_flat[g.at(8)] == g.flat(7));
    }

    {
        const line_grid higher = make_line({3, 2, 1, 0, 1, 2, 3.5f}, {3}, dx);
        const auto flow = make_flow(higher.z, higher.mask, higher.meta, "nnnn", {higher.flat(3)});
        const auto r = relief(flow, higher.z, higher.mask, higher.channel, "farthest", "channel_pixel");
        check("farthest tie: the higher pixel", r.ridge_flat[higher.at(1)] == higher.flat(6));

        const line_grid equal = make_line({3, 2, 1, 0, 1, 2, 3}, {3}, dx);
        const auto flow2 = make_flow(equal.z, equal.mask, equal.meta, "nnnn", {equal.flat(3)});
        const auto r2 = relief(flow2, equal.z, equal.mask, equal.channel, "farthest", "channel_pixel");
        check("farthest tie: then the lower flat index", r2.ridge_flat[equal.at(5)] == equal.flat(0));

        const line_grid high_tie = make_line({9, 2, 1, 0, 1, 2, 9, 3}, {3}, dx);
        const auto flow3 = make_flow(high_tie.z, high_tie.mask, high_tie.meta, "nnnn",
                                     {high_tie.flat(3), high_tie.flat(7)});
        const auto r3 = relief(flow3, high_tie.z, high_tie.mask, high_tie.channel, "highest",
                               "channel_pixel");
        check("highest tie: the lower flat index", r3.ridge_flat[high_tie.at(5)] == high_tie.flat(0));
    }

    {
        // A channel fragment at column 3, upslope of the channel at column 6.
        const line_grid g = make_line({6, 5, 4, 3, 2, 1, 0}, {3, 6}, dx);
        const auto flow = make_flow(g.z, g.mask, g.meta, "nnnn", {g.flat(6)});
        for (const std::string grouping : {"channel_pixel", "entry_pixel", "pixel"})
        {
            const auto r = relief(flow, g.z, g.mask, g.channel, "farthest", grouping);
            check("channel fragment stops the search (" + grouping + ")",
                  r.ridge_flat[g.at(5)] == g.flat(4) && r.channel_flat[g.at(5)] == g.flat(6) &&
                  r.channel_flat[g.at(1)] == g.flat(3) && r.ridge_flat[g.at(1)] == g.flat(0));
        }
    }

    {
        // Column 2 is a channel outlet; column 6 is an outlet that is not a channel.
        const line_grid g = make_line({2, 1, 0, 1, 2, 1, 0}, {2}, dx);
        const auto flow = make_flow(g.z, g.mask, g.meta, "nnnn", {g.flat(2), g.flat(6)});
        const auto r = relief(flow, g.z, g.mask, g.channel, "farthest", "channel_pixel");
        const int above = r.channel_flat[g.at(4)] >= 0 ? 3 : 2;
        const std::int64_t expected_none = 7 - 1 - (2 + above - 1);
        check("pixels reaching no channel are counted", r.n_no_channel == expected_none,
              "n_no_channel " + std::to_string(r.n_no_channel));
        check("pixels reaching no channel are nodata", r.hand.mask[g.at(6)] == 0u &&
              r.relative_height.mask[g.at(5)] == 0u && r.channel_flat[g.at(6)] == -1 &&
              r.ridge_status[g.at(6)] ==
                  static_cast<std::uint8_t>(hillslope_ridge_status::no_channel));
        check("a channel pixel is not a hillslope pixel", r.hand.mask[g.at(2)] == 0u &&
              r.channel_flat[g.at(2)] == g.flat(2) && r.ridge_flat[g.at(2)] == -1);
    }

    {
        // Routed on a descending row, with heights taken from a flat surface.
        const line_grid g = make_line({3, 2, 1, 0}, {3}, dx);
        const auto flow = make_flow(g.z, g.mask, g.meta, "nnnn", {g.flat(3)});
        std::vector<float> flat = g.z;
        for (int c = 0; c < 4; ++c) flat[g.at(c)] = 5.f;
        const auto r = relief(flow, flat, g.mask, g.channel, "farthest", "channel_pixel");
        check("zero relief: relative height is nodata and counted",
              r.n_zero_relief == 3 && r.relative_height.mask[g.at(1)] == 0u &&
              r.hand.mask[g.at(1)] == 1u && r.hand.data[g.at(1)] == 0.f);
    }
}

void test_surface()
{
    std::cout << "Rough routed surface\n";
    const int nr = 150, nc = 160;
    const double dx = 30.0;
    const lsdtt_raster_meta meta = make_meta(nr, nc, dx);

    std::mt19937 rng(20261122u);
    std::normal_distribution<float> noise(0.f, 6.f);
    std::vector<float> z0(nr * nc);
    for (int r = 0; r < nr; ++r)
        for (int c = 0; c < nc; ++c)
            z0[r * nc + c] = static_cast<float>(
                1200.0 + 2.0 * r + 120.0 * std::sin(0.07 * c) * std::cos(0.05 * r) +
                45.0 * std::sin(0.23 * c + 0.11 * r)) + noise(rng);
    const lsdtt_raster_cpu filled = fill_pits(
        lsdtt_raster_cpu(std::vector<float>(z0), std::vector<std::uint8_t>(nr * nc, 1u), meta));
    const std::vector<float> z = filled.get_data_array().get_vector();
    const std::vector<std::uint8_t> mask = filled.get_mask_array().get_vector();

    lsdtt_flowinfo_options options;
    options.boundary_conditions = "oooo";
    const lsdtt_flowinfo flow(filled.get_data_array(), filled.get_mask_array(), meta, options);

    const lsdtt_raster_result contributing = flow.get_n_contributing_pixels_raster();
    std::vector<std::uint8_t> channel(nr * nc, 0u);
    for (std::size_t i = 0; i < channel.size(); ++i)
        if (contributing.mask[i] && contributing.data[i] >= 60.f) channel[i] = 1u;

    const auto existing = flow.get_height_above_channel_mask_donor_raster(z, mask, channel);
    const auto far_own = relief(flow, z, mask, channel, "farthest", "pixel");
    {
        std::int64_t diff_link = 0, diff_hand = 0;
        for (std::size_t i = 0; i < z.size(); ++i)
        {
            if (channel[i]) continue;
            if (far_own.channel_flat[i] != existing.donor_flat[i]) ++diff_link;
            if (far_own.hand.mask[i] != existing.height.mask[i] ||
                (far_own.hand.mask[i] && far_own.hand.data[i] != existing.height.data[i]))
                ++diff_hand;
        }
        check("channel link equals the existing HAND donor", diff_link == 0,
              std::to_string(diff_link) + " differ");
        check("HAND is bit-identical to the existing routine", diff_hand == 0,
              std::to_string(diff_hand) + " differ");
        check("the surface has hillslope and channel pixels",
              far_own.n_hillslope > 10000 && far_own.n_channel > 200,
              std::to_string(far_own.n_hillslope) + " hillslope, " +
              std::to_string(far_own.n_channel) + " channel");
    }

    {
        // On a closed network, literal per-pixel grouping is the existing
        // farthest-ridge distance.
        const lsdtt_raster_result old_distance = flow.get_distance_to_farthest_ridge_raster();
        double worst = 0.0;
        for (std::size_t i = 0; i < z.size(); ++i)
            if (far_own.hand.mask[i])
                worst = std::max(worst, std::abs(static_cast<double>(far_own.ridge_distance.data[i]) -
                                                 far_own.channel_distance.data[i] -
                                                 old_distance.data[i]));
        check("pixel grouping matches the existing farthest-ridge distance", worst < 1e-2,
              "max difference " + num(worst) + " m");
    }

    for (const std::string ridge : {"farthest", "highest"})
        for (const std::string grouping : {"channel_pixel", "entry_pixel"})
        {
            const auto res = relief(flow, z, mask, channel, ridge, grouping);
            const std::string tag = ridge + "/" + grouping;

            // Brute-force ridge of every group.
            std::map<std::int64_t, std::int64_t> brute;
            const auto better = [&](std::int64_t a, std::int64_t b) {
                const double la = res.channel_distance.data[a], lb = res.channel_distance.data[b];
                if (ridge == "farthest")
                {
                    if (std::abs(la - lb) > 1e-3) return la > lb;
                    if (z[a] != z[b]) return z[a] > z[b];
                }
                else
                {
                    if (z[a] != z[b]) return z[a] > z[b];
                    if (std::abs(la - lb) > 1e-3) return la > lb;
                }
                return a < b;
            };
            for (std::int64_t i = 0; i < static_cast<std::int64_t>(z.size()); ++i)
            {
                if (!res.hand.mask[i]) continue;
                const std::int64_t key = grouping == "channel_pixel" ? res.channel_flat[i]
                                                                     : res.entry_flat[i];
                auto it = brute.find(key);
                if (it == brute.end()) brute[key] = i;
                else if (better(i, it->second)) it->second = i;
            }

            std::int64_t wrong_ridge = 0, rises = 0, out_of_bounds = 0, above = 0, wrong_above = 0;
            for (std::int64_t i = 0; i < static_cast<std::int64_t>(z.size()); ++i)
            {
                if (!res.hand.mask[i]) continue;
                const std::int64_t key = grouping == "channel_pixel" ? res.channel_flat[i]
                                                                     : res.entry_flat[i];
                if (res.ridge_flat[i] != brute[key]) ++wrong_ridge;
                if (!res.relative_height.mask[i]) continue;
                const float value = res.relative_height.data[i];
                if (value < 0.f) ++out_of_bounds;
                const double exact = (static_cast<double>(z[i]) - z[res.channel_flat[i]]) /
                                     (static_cast<double>(z[res.ridge_flat[i]]) - z[res.channel_flat[i]]);
                if (exact > 1.0 + 1e-6)
                {
                    ++above;
                    if (z[i] <= z[res.ridge_flat[i]]) ++wrong_above;
                }
                const std::int64_t rec = flow.receiver_of_unchecked(i);
                if (rec >= 0 && res.hand.mask[rec] && res.relative_height.mask[rec] &&
                    res.relative_height.data[rec] > value + 1e-5f)
                    ++rises;
            }
            check(tag + ": ridge equals a brute-force search of the group", wrong_ridge == 0,
                  std::to_string(wrong_ridge) + " differ");
            check(tag + ": relative height never rises downslope", rises == 0,
                  std::to_string(rises) + " links");
            check(tag + ": relative height is not negative", out_of_bounds == 0);
            check(tag + ": a ratio above 1 means the pixel is above its ridge", wrong_above == 0,
                  std::to_string(above) + " above 1");
            if (ridge == "highest")
                check(tag + ": relative height is at most 1", above == 0 && res.n_above_ridge == 0);
        }
}

void test_refusals()
{
    std::cout << "Refusals\n";
    const line_grid g = make_line({3, 2, 1, 0}, {3}, 10.0);
    const auto flow = make_flow(g.z, g.mask, g.meta, "nnnn", {g.flat(3)});
    const auto throws = [&](const std::string& ridge, const std::string& grouping) {
        try { (void)relief(flow, g.z, g.mask, g.channel, ridge, grouping); }
        catch (const std::invalid_argument&) { return true; }
        return false;
    };
    check("unknown ridge is refused", throws("nearest", "channel_pixel"));
    check("unknown grouping is refused", throws("farthest", "basin"));
    bool size_refused = false;
    try
    {
        (void)compute_hillslope_relief(flow, std::vector<float>(3, 0.f), g.mask, g.channel);
    }
    catch (const std::invalid_argument&) { size_refused = true; }
    check("wrong array size is refused", size_refused);
}

} // namespace

int main()
{
    try
    {
        test_valleys();
        test_paths();
        test_surface();
        test_refusals();
    }
    catch (const std::exception& e)
    {
        std::cout << "  FAIL  unexpected exception: " << e.what() << "\n";
        ++g_failures;
    }
    std::cout << g_checks << " checks, " << g_failures << " failed\n";
    std::cout << "RESULT: " << (g_failures == 0 ? "PASS" : "FAIL") << "\n";
    return g_failures == 0 ? 0 : 1;
}
