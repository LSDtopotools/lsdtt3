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
//   * divide edges: exact on hand-built networks, including a junction, a basin
//     with no channel and a diagonal-only contact;
//   * aspect rays: analytic valleys, every octant, exact corners, the crest and
//     at-or-above rules, each stopping condition, agreement with a fine-step
//     reference ray, polyfit and D-infinity aspect on planes, threads;
//   * tags: exact values through the D8 link, partial coverage and exclusion
//     counts, the choice of valid basins for chi, chi and k_sn against the
//     hilltop join;
//   * imported attributes: CSV and FlatGeobuf, exact, coordinate and line joins,
//     reprojection, medians, gap filling along a link, refusals;
//   * refusals.
//
// Usage: lsdtt-test-hillslope-channel-tagging
// Prints one line per check and a final "RESULT: PASS|FAIL".
// -----------------------------------------------------------------------------

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

#include "../../src/core/lsdtt_channel_attribute_input.hpp"
#include "../../src/core/lsdtt_data_loader.hpp"
#include "../../src/core/lsdtt_raster_meta.hpp"
#include "../../src/core/lsdtt_vector_utils.hpp"
#include "../../src/cpu/lsdtt_flowinfo_cpu.hpp"
#include "../../src/cpu/lsdtt_chi_extraction_cpu.hpp"
#include "../../src/cpu/lsdtt_chi_network_analysis_cpu.hpp"
#include "../../src/cpu/lsdtt_hillslope_channel_tagging_cpu.hpp"
#include "../../src/cpu/lsdtt_hilltop_flow_routing_cpu.hpp"
#include "../../src/cpu/lsdtt_raster_cpu.hpp"
#include "../../src/cpu/lsdtt_raster_flow_preprocessing_cpu.hpp"
#include "../../src/utils/lsdtt_raster_geometry.hpp"

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

// Aspect is downslope, compass degrees from grid north (up the grid).
lsdtt_raster_result constant_aspect(const lsdtt_raster_meta& meta, float degrees)
{
    const std::size_t n = static_cast<std::size_t>(meta.n_rows) * meta.n_cols;
    lsdtt_raster_result a;
    a.meta = meta;
    a.data.assign(n, degrees);
    a.mask.assign(n, 1u);
    return a;
}

hillslope_divide_edges no_divides(const lsdtt_raster_meta& meta)
{
    const std::size_t n = static_cast<std::size_t>(meta.n_rows) * meta.n_cols;
    hillslope_divide_edges d;
    d.n_rows = meta.n_rows;
    d.n_cols = meta.n_cols;
    d.east.assign(n, 0u);
    d.south.assign(n, 0u);
    return d;
}

// A surface used only for heights: routed on a tilted plane to one channel row,
// so that a ray's behaviour is set entirely by the heights and edges supplied.
struct ray_grid
{
    lsdtt_raster_meta meta;
    std::vector<float> route, z;
    std::vector<std::uint8_t> mask, channel;
    int nr, nc;
    std::size_t at(int r, int c) const { return static_cast<std::size_t>(r) * nc + c; }
};

ray_grid make_ray_grid(int nr, int nc, double dx = 10.0)
{
    ray_grid g{make_meta(nr, nc, dx), std::vector<float>(nr * nc), std::vector<float>(nr * nc, 100.f),
               std::vector<std::uint8_t>(nr * nc, 1u), std::vector<std::uint8_t>(nr * nc, 0u), nr, nc};
    for (int r = 0; r < nr; ++r)
        for (int c = 0; c < nc; ++c)
        {
            g.route[r * nc + c] = static_cast<float>(1000.0 + 5.0 * (nr - 1 - r) + 0.001 * c);
            if (r == nr - 1) g.channel[r * nc + c] = 1u;
        }
    return g;
}

hillslope_relief_result rays(const ray_grid& g, const lsdtt_flowinfo& flow,
                             const lsdtt_raster_result& aspect,
                             const hillslope_divide_edges& divides, int threads = 1)
{
    hillslope_relief_options o;
    o.ridge = "aspect_ray";
    o.ray_aspect = &aspect;
    o.divides = &divides;
    o.n_threads = threads;
    return compute_hillslope_relief(flow, g.z, g.mask, g.channel, o);
}

void test_divides()
{
    std::cout << "Divide edges\n";
    const double dx = 10.0;

    {
        // Two hillslopes on one channel row share no divide but meet at a peak.
        const line_grid g = make_line({0, 1, 2, 3, 2, 1, 0}, {0, 6}, dx);
        const auto flow = make_flow(g.z, g.mask, g.meta, "nnnn", {g.flat(0), g.flat(6)});
        const auto d = build_hillslope_divide_edges(flow, g.mask, g.channel);
        // Column 3 falls to one side, so the single divide is next to it.
        int found = -1, count = 0;
        for (int c = 0; c + 1 < g.nc; ++c)
            if (d.east[g.at(c)]) { found = c; ++count; }
        const bool to_west = flow.receiver_of_unchecked(g.flat(3)) == g.flat(2);
        check("line: one divide, at the peak", count == 1 && found == (to_west ? 3 : 2),
              "at column " + std::to_string(found));
        check("line: two channel links", d.n_channel_links == 2);
    }

    {
        // A Y network: two tributaries join at a junction and drain south.
        const int nr = 9, nc = 9;
        std::vector<float> z(nr * nc);
        std::vector<std::uint8_t> mask(nr * nc, 1u), channel(nr * nc, 0u);
        for (int r = 0; r < nr; ++r)
            for (int c = 0; c < nc; ++c)
            {
                const double valley = r >= 4 ? std::abs(c - 4) : std::min(std::abs(c - 2), std::abs(c - 6));
                z[r * nc + c] = static_cast<float>(100.0 + 3.0 * (nr - 1 - r) + 40.0 * valley);
            }
        for (int r = 4; r < nr; ++r) channel[r * nc + 4] = 1u;
        for (int r = 0; r < 3; ++r) { channel[r * nc + 2] = 1u; channel[r * nc + 6] = 1u; }
        channel[3 * nc + 3] = 1u;
        channel[3 * nc + 5] = 1u;
        z[3 * nc + 3] = 100.f + 3.f * 5 - 1.f;
        z[3 * nc + 5] = 100.f + 3.f * 5 - 1.f;
        z[3 * nc + 2] = 200.f;
        z[3 * nc + 6] = 200.f;
        z[3 * nc + 4] = 200.f;
        z[2 * nc + 2] = 100.f + 3.f * 6;
        z[2 * nc + 6] = 100.f + 3.f * 6;
        const auto flow = make_flow(z, mask, make_meta(nr, nc, dx), "nnon");
        const auto d = build_hillslope_divide_edges(flow, mask, channel);
        check("Y network: three links", d.n_channel_links == 3,
              std::to_string(d.n_channel_links) + " links");
        std::int64_t channel_edges = 0;
        for (int i = 0; i < nr * nc; ++i)
        {
            if (d.east[i] && (channel[i] || channel[i + 1])) ++channel_edges;
            if (d.south[i] && (channel[i] || channel[i + nc])) ++channel_edges;
        }
        check("Y network: no divide edge touches a channel", channel_edges == 0);
        check("Y network: the tributaries are separated", d.n_edges > 0,
              std::to_string(d.n_edges) + " edges");
    }

    {
        // A pixel reaching no channel is its own catchment; nodata is not a divide.
        const line_grid g = make_line({0, 1, 2, 1, 0}, {0}, dx);
        const auto flow = make_flow(g.z, g.mask, g.meta, "nnnn", {g.flat(0), g.flat(4)});
        const auto d = build_hillslope_divide_edges(flow, g.mask, g.channel);
        int count = 0;
        for (std::size_t i = 0; i < d.east.size(); ++i) count += d.east[i] + d.south[i];
        check("catchment with no channel: one divide, none against nodata", count == 1);
    }

    {
        // Two cells of one basin touching only at a corner (a jagged coast): the
        // labels, not an outline, define the edges, so this cannot fail.
        const int nr = 4, nc = 4;
        std::vector<float> z = {9,  -9999, -9999, -9999,
                                -9999, 5, -9999, -9999,
                                -9999, 4, 3, -9999,
                                -9999, -9999, -9999, 1};
        std::vector<std::uint8_t> mask(nr * nc, 0u), channel(nr * nc, 0u);
        for (int i = 0; i < nr * nc; ++i) mask[i] = z[i] > -9000.f ? 1u : 0u;
        channel[15] = 1u;
        const auto flow = make_flow(z, mask, make_meta(nr, nc, dx), "nnnn", {15});
        const auto d = build_hillslope_divide_edges(flow, mask, channel);
        check("diagonal-only contact: one catchment, no divide", d.n_edges == 0 && d.n_channel_links == 1);
    }
}

void test_rays()
{
    std::cout << "Aspect rays\n";

    {
        // Planar valley, ridge at both grid ends as a marked divide column.
        const int nr = 12, nc = 21, c0 = 10;
        const double S = 0.5, dx = 10.0;
        const valley v = make_valley(nr, nc, c0, S, dx);
        const lsdtt_flowinfo flow = make_flow(v.z, v.mask, v.meta, "nnon");
        lsdtt_raster_result aspect = constant_aspect(v.meta, 0.f);
        for (int r = 0; r < nr; ++r)
            for (int c = 0; c < nc; ++c)
                aspect.data[r * nc + c] = c < c0 ? 90.f : 270.f;
        hillslope_divide_edges d = no_divides(v.meta);
        for (int r = 0; r < nr; ++r) { d.east[r * nc + 3] = 1u; d.east[r * nc + 16] = 1u; }
        // A ridge along the divide lines: columns 3-4 and 16-17 are the crest.
        std::vector<float> z = v.z;
        for (int r = 0; r < nr; ++r)
            for (int c = 0; c < nc; ++c)
            {
                const int d_c = std::abs(c - c0);
                const int height = d_c <= 6 ? d_c : 13 - d_c;
                z[r * nc + c] = static_cast<float>(1000.0 + S * dx * height + 0.01 * dx * (nr - 1 - r));
            }
        hillslope_relief_options o;
        o.ridge = "aspect_ray";
        o.ray_aspect = &aspect;
        o.divides = &d;
        o.n_threads = 1;
        const auto res = compute_hillslope_relief(flow, z, v.mask, v.channel, o);
        double err = 0.0;
        for (int r = 0; r < nr - 1; ++r)
            for (const int c : {5, 6, 7, 8, 9, 11, 12, 13, 14, 15})
            {
                const std::size_t i = r * nc + c;
                const int d_c = std::abs(c - c0);
                // The two pixels on each divide are equally high: the lower flat index.
                const int ridge_c = c < c0 ? 3 : 16;
                err = std::max(err, std::abs(res.relative_height.data[i] - d_c / 6.0));
                err = std::max(err, std::abs(res.ridge_distance.data[i] - dx * std::abs(c - ridge_c)));
                if (res.ridge_flat[i] != static_cast<std::int64_t>(r) * nc + ridge_c) err = 1e9;
                if (res.ridge_status[i] != static_cast<std::uint8_t>(hillslope_ridge_status::crest_divide))
                    err = 1e9;
            }
        check("planar valley: ray finds the divide crest; ratio and distance", err < 1e-3,
              "max error " + num(err));
    }

    // Each compass octant from the centre of an open grid with a ring of divides.
    {
        ray_grid g = make_ray_grid(41, 41);
        const lsdtt_flowinfo flow = make_flow(g.route, g.mask, g.meta, "nnon");
        const int r0 = 20, c0 = 20, R = 6;
        struct oct { float aspect; int dr, dc; };  // the upslope step is opposite the aspect
        const std::vector<oct> octants = {{180.f, -1, 0}, {225.f, -1, 1}, {270.f, 0, 1}, {315.f, 1, 1},
                                          {0.f, 1, 0}, {45.f, 1, -1}, {90.f, 0, -1}, {135.f, -1, -1}};
        for (const oct& o : octants)
        {
            ray_grid h = g;
            hillslope_divide_edges d = no_divides(g.meta);
            const int rt = r0 + o.dr * R, ct = c0 + o.dc * R;
            const int rb = rt - o.dr, cb = ct - o.dc;
            // Mark the edges of the step into the target cell.
            if (o.dr != 0 && o.dc != 0)
            {
                d.south[h.at(std::min(rb, rt), cb)] = 1u;
                d.east[h.at(rt, std::min(cb, ct))] = 1u;
            }
            else if (o.dr != 0) d.south[h.at(std::min(rb, rt), ct)] = 1u;
            else d.east[h.at(rt, std::min(cb, ct))] = 1u;
            h.z[h.at(rt, ct)] = 150.f;
            const auto a = constant_aspect(g.meta, o.aspect);
            const auto res = rays(h, flow, a, d);
            const std::size_t i = h.at(r0, c0);
            const double expected = 10.0 * std::hypot(o.dr * R, o.dc * R);
            check("octant, aspect " + num(o.aspect) + ": stops at the marked divide",
                  res.ridge_flat[i] == static_cast<std::int64_t>(h.at(rt, ct)) &&
                  std::abs(res.ridge_distance.data[i] - expected) < 1e-3 &&
                  res.ridge_status[i] == static_cast<std::uint8_t>(hillslope_ridge_status::crest_divide));
        }
    }

    {
        // The next cells along an eastward upslope ray from column 2 of row 5.
        ray_grid g = make_ray_grid(11, 20);
        const lsdtt_flowinfo flow = make_flow(g.route, g.mask, g.meta, "nnon");
        const auto a = constant_aspect(g.meta, 270.f);  // faces west, so upslope is east
        const std::size_t start = g.at(5, 2);

        {
            ray_grid h = g;
            hillslope_divide_edges d = no_divides(g.meta);
            d.east[h.at(5, 5)] = 1u;
            d.east[h.at(5, 10)] = 1u;
            h.z[h.at(5, 5)] = 90.f;
            h.z[h.at(5, 6)] = 95.f;
            h.z[h.at(5, 10)] = 120.f;
            const auto r = rays(h, flow, a, d);
            check("a divide below the start is passed", r.ridge_flat[start] == h.at(5, 10));
        }
        {
            ray_grid h = g;
            hillslope_divide_edges d = no_divides(g.meta);
            d.east[h.at(5, 5)] = 1u;
            d.east[h.at(5, 12)] = 1u;
            h.z[h.at(5, 5)] = 110.f;
            h.z[h.at(5, 6)] = 115.f;
            h.z[h.at(5, 7)] = 130.f;
            h.z[h.at(5, 12)] = 140.f;
            const auto r = rays(h, flow, a, d);
            check("a non-crest divide with rising ground beyond is passed",
                  r.ridge_flat[start] == h.at(5, 12));
        }
        {
            ray_grid h = g;
            hillslope_divide_edges d = no_divides(g.meta);
            d.east[h.at(5, 5)] = 1u;
            h.z[h.at(5, 5)] = 110.f;
            h.z[h.at(5, 6)] = 120.f;
            h.z[h.at(5, 7)] = 115.f;
            const auto r = rays(h, flow, a, d);
            check("higher entered cell with lower ground beyond is the crest",
                  r.ridge_flat[start] == h.at(5, 6));
        }
        {
            ray_grid h = g;
            hillslope_divide_edges d = no_divides(g.meta);
            d.east[h.at(5, 5)] = 1u;
            h.z[h.at(5, 5)] = 110.f;
            h.z[h.at(5, 6)] = 105.f;
            h.z[h.at(5, 7)] = 200.f;
            const auto r = rays(h, flow, a, d);
            check("a crest stops the ray before higher ground further on",
                  r.ridge_flat[start] == h.at(5, 5));
        }
        {
            ray_grid h = g;
            h.z[h.at(5, 9)] = 170.f;
            h.z[h.at(5, 14)] = 130.f;
            const auto r = rays(h, flow, a, no_divides(g.meta));
            check("no divide: the highest cell on the ray", r.ridge_flat[start] == h.at(5, 9));
            check("no divide, off the grid: status", r.ridge_status[start] ==
                  static_cast<std::uint8_t>(hillslope_ridge_status::ray_boundary));
        }
        {
            ray_grid h = g;
            h.channel[h.at(5, 8)] = 1u;
            h.z[h.at(5, 6)] = 150.f;
            h.z[h.at(5, 12)] = 300.f;
            const auto r = rays(h, flow, a, no_divides(g.meta));
            check("a channel stops the ray: highest cell before it",
                  r.ridge_flat[start] == h.at(5, 6) && r.ridge_status[start] ==
                  static_cast<std::uint8_t>(hillslope_ridge_status::ray_channel));
        }
        {
            ray_grid h = g;
            h.z[start] = 500.f;
            const auto r = rays(h, flow, a, no_divides(g.meta));
            check("start pixel highest: its own ridge, ratio 1",
                  r.ridge_flat[start] == static_cast<std::int64_t>(start) &&
                  std::abs(r.relative_height.data[start] - 1.f) < 1e-6f);
        }
        {
            auto missing = a;
            missing.mask[start] = 0u;
            const auto r = rays(g, flow, missing, no_divides(g.meta));
            check("no aspect: HAND kept, no ridge, counted",
                  r.ridge_flat[start] == -1 && r.ridge_relief.mask[start] == 0u &&
                  r.relative_height.mask[start] == 0u && r.hand.mask[start] == 1u &&
                  r.n_no_aspect == 1 && r.ridge_status[start] ==
                  static_cast<std::uint8_t>(hillslope_ridge_status::no_aspect));
        }
    }

    {
        // Exact diagonal between two channel pixels: must not pass.
        ray_grid g = make_ray_grid(21, 21);
        const lsdtt_flowinfo flow = make_flow(g.route, g.mask, g.meta, "nnon");
        const auto a = constant_aspect(g.meta, 225.f);  // upslope is north-east
        ray_grid h = g;
        h.channel[h.at(9, 10)] = 1u;
        h.channel[h.at(10, 11)] = 1u;
        h.z[h.at(5, 15)] = 300.f;
        const auto r = rays(h, flow, a, no_divides(g.meta));
        check("an exact diagonal cannot pass between two channel cells",
              r.ridge_status[h.at(10, 10)] ==
              static_cast<std::uint8_t>(hillslope_ridge_status::ray_channel) &&
              r.ridge_flat[h.at(10, 10)] == static_cast<std::int64_t>(h.at(10, 10)));
    }

    {
        // Every touched cell against a very fine-step reference ray.
        const int nr = 61, nc = 61;
        ray_grid g = make_ray_grid(nr, nc);
        const lsdtt_flowinfo flow = make_flow(g.route, g.mask, g.meta, "nnon");
        std::mt19937 rng(7u);
        std::uniform_real_distribution<float> height(0.f, 100.f), bearing(0.f, 360.f);
        for (int i = 0; i < nr * nc; ++i) g.z[i] = 1000.f + height(rng);
        lsdtt_raster_result aspect = constant_aspect(g.meta, 0.f);
        for (float& b : aspect.data) b = bearing(rng);
        const auto res = rays(g, flow, aspect, no_divides(g.meta));

        std::int64_t wrong = 0, n = 0;
        for (int r = 5; r < nr - 2; r += 3)
            for (int c = 2; c < nc - 2; c += 3)
            {
                const std::size_t i = g.at(r, c);
                if (!res.hand.mask[i]) continue;
                const double b = (aspect.data[i] + 180.0) * 3.14159265358979323846 / 180.0;
                const double dr = -std::cos(b), dc = std::sin(b);
                std::int64_t best = static_cast<std::int64_t>(i);
                for (double t = 0.0; t < 200.0; t += 0.0005)
                {
                    const int rr = static_cast<int>(std::floor(r + 0.5 + t * dr));
                    const int cc = static_cast<int>(std::floor(c + 0.5 + t * dc));
                    if (rr < 0 || rr >= nr || cc < 0 || cc >= nc || g.channel[g.at(rr, cc)]) break;
                    const std::int64_t k = g.at(rr, cc);
                    if (g.z[k] > g.z[best] || (g.z[k] == g.z[best] && k < best)) best = k;
                }
                ++n;
                if (best != res.ridge_flat[i]) ++wrong;
            }
        check("highest cell on the ray equals a fine-step reference", wrong == 0 && n > 100,
              std::to_string(wrong) + " of " + std::to_string(n) + " differ");
    }

    {
        // Aspect sources on planes that rise east, south and north-east.
        const int nr = 41, nc = 41;
        const lsdtt_raster_meta meta = make_meta(nr, nc, 10.0);
        struct plane { double east, south, expected; };
        for (const plane& p : {plane{1, 0, 270}, plane{0, 1, 0}, plane{1, -1, 225}})
        {
            std::vector<float> z(nr * nc);
            for (int r = 0; r < nr; ++r)
                for (int c = 0; c < nc; ++c)
                    z[r * nc + c] = static_cast<float>(500.0 + 3.0 * (p.east * c + p.south * r));
            const lsdtt_raster_cpu surface(std::vector<float>(z),
                                           std::vector<std::uint8_t>(nr * nc, 1u), meta);
            const auto poly = compute_hillslope_ray_aspect(surface, "polyfit", 30.0);
            const auto dinf = compute_hillslope_ray_aspect(surface, "dinf", 0.0);
            const std::size_t i = 20 * nc + 20;
            const auto off = [&](float a) {
                const double e = std::fmod(std::abs(a - p.expected), 360.0);
                return std::min(e, 360.0 - e);
            };
            check("plane, expected aspect " + num(p.expected) + ": polyfit and D-infinity",
                  off(poly.data[i]) < 1e-3 && off(dinf.data[i]) < 1e-3,
                  "polyfit " + num(poly.data[i]) + ", dinf " + num(dinf.data[i]));
        }
    }

    {
        // Real pipeline on a rough surface: bounds, thread independence, identity.
        const int nr = 150, nc = 160;
        const lsdtt_raster_meta meta = make_meta(nr, nc, 30.0);
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
        lsdtt_flowinfo_options fo;
        fo.boundary_conditions = "oooo";
        const lsdtt_flowinfo flow(filled.get_data_array(), filled.get_mask_array(), meta, fo);
        const lsdtt_raster_result cp = flow.get_n_contributing_pixels_raster();
        std::vector<std::uint8_t> channel(nr * nc, 0u);
        for (std::size_t i = 0; i < channel.size(); ++i)
            if (cp.mask[i] && cp.data[i] >= 60.f) channel[i] = 1u;

        const auto aspect = compute_hillslope_ray_aspect(filled, "polyfit", 90.0);
        const auto d = build_hillslope_divide_edges(flow, mask, channel);
        const auto run = [&](int threads) {
            hillslope_relief_options o;
            o.ridge = "aspect_ray";
            o.ray_aspect = &aspect;
            o.divides = &d;
            o.n_threads = threads;
            return compute_hillslope_relief(flow, z, mask, channel, o);
        };
        const auto one = run(1);
        const auto two = run(2);
        const auto eight = run(8);
        check("rays: identical at 1, 2 and 8 threads",
              one.ridge_flat == two.ridge_flat && one.ridge_flat == eight.ridge_flat &&
              one.relative_height.data == eight.relative_height.data &&
              one.ridge_status == eight.ridge_status);

        std::int64_t out_of_bounds = 0;
        for (std::size_t i = 0; i < z.size(); ++i)
            if (one.relative_height.mask[i] &&
                (one.relative_height.data[i] < 0.f || one.relative_height.data[i] > 1.f + 1e-6f))
                ++out_of_bounds;
        check("rays: relative height lies in [0, 1]", out_of_bounds == 0,
              std::to_string(out_of_bounds) + " outside");
        check("rays: every hillslope pixel has a status",
              one.n_crest_divide + one.n_ray_channel + one.n_ray_boundary + one.n_no_aspect ==
              one.n_hillslope);

        const auto far = relief(flow, z, mask, channel, "farthest", "channel_pixel");
        check("the channel link and HAND do not depend on the ridge method",
              one.channel_flat == far.channel_flat && one.hand.data == far.hand.data &&
              one.channel_distance.data == far.channel_distance.data);

        std::int64_t touching = 0, same_label_side = 0;
        for (int r = 0; r < nr; ++r)
            for (int c = 0; c < nc; ++c)
            {
                const std::size_t i = static_cast<std::size_t>(r) * nc + c;
                if (d.east[i] && (channel[i] || channel[i + 1])) ++touching;
                if (d.south[i] && (channel[i] || channel[i + nc])) ++touching;
                if (d.east[i] && far.channel_flat[i] >= 0 &&
                    far.channel_flat[i] == far.channel_flat[i + 1]) ++same_label_side;
                if (d.south[i] && far.channel_flat[i] >= 0 &&
                    far.channel_flat[i] == far.channel_flat[i + nc]) ++same_label_side;
            }
        check("divides: none touches a channel or lies inside one channel pixel's hillslope",
              touching == 0 && same_label_side == 0 && d.n_edges > 500,
              std::to_string(d.n_edges) + " edges");
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

    const auto aspect = constant_aspect(g.meta, 0.f);
    const auto ray_throws = [&](hillslope_relief_options o) {
        try { (void)compute_hillslope_relief(flow, g.z, g.mask, g.channel, o); }
        catch (const std::invalid_argument&) { return true; }
        return false;
    };
    hillslope_relief_options o;
    o.ridge = "aspect_ray";
    check("aspect_ray without an aspect is refused", ray_throws(o));
    o.ray_aspect = &aspect;
    o.grouping = "entry_pixel";
    check("a grouping with aspect_ray is refused", ray_throws(o));

    bool bad_method = false;
    try
    {
        const lsdtt_raster_cpu surface(std::vector<float>(g.z), std::vector<std::uint8_t>(g.mask), g.meta);
        (void)compute_hillslope_ray_aspect(surface, "sobel", 30.0);
    }
    catch (const std::invalid_argument&) { bad_method = true; }
    check("unknown aspect method is refused", bad_method);

    {
        const int nr = 8, nc = 8;
        const lsdtt_raster_meta meta = make_meta(nr, nc, 10.0);
        std::vector<float> z(nr * nc);
        for (int r = 0; r < nr; ++r)
            for (int c = 0; c < nc; ++c) z[r * nc + c] = static_cast<float>(100.0 + 3.0 * (nr - r) + 0.1 * c);
        const std::vector<std::uint8_t> mask(nr * nc, 1u), channel(nr * nc, 0u);
        const auto periodic = make_flow(z, mask, meta, "popo");
        bool refused = false;
        try { (void)build_hillslope_divide_edges(periodic, mask, channel); }
        catch (const std::invalid_argument&) { refused = true; }
        check("periodic boundaries are refused", refused);
    }
}

// Two parallel valleys inside a nodata collar. The right valley stands higher,
// so its outlet is above the left one and only an elevation limit separates them.
struct twin_valleys
{
    lsdtt_raster_meta meta;
    std::vector<float> z;
    std::vector<std::uint8_t> mask, channel;
    std::int64_t left_outlet, right_outlet;
    int nr, nc;
};

twin_valleys make_twin_valleys()
{
    const int nr = 30, nc = 21, left = 5, right = 15;
    twin_valleys v{make_meta(nr, nc, 10.0), std::vector<float>(nr * nc, -9999.f),
                   std::vector<std::uint8_t>(nr * nc, 0u), std::vector<std::uint8_t>(nr * nc, 0u),
                   (nr - 2) * nc + left, (nr - 2) * nc + right, nr, nc};
    for (int r = 1; r < nr - 1; ++r)
        for (int c = 1; c < nc - 1; ++c)
        {
            const bool is_left = c <= 10;
            const int c0 = is_left ? left : right;
            v.z[r * nc + c] = static_cast<float>((is_left ? 10.0 : 200.0) + 5.0 * std::abs(c - c0) +
                                                 2.0 * (nr - 2 - r) + 0.0005 * r * r);
            v.mask[r * nc + c] = 1u;
            if (c == c0) v.channel[r * nc + c] = 1u;
        }
    return v;
}

void test_tags()
{
    std::cout << "\n--- tags ---\n";
    const twin_valleys v = make_twin_valleys();
    const std::size_t n = v.z.size();
    const lsdtt_flowinfo flow = make_flow(v.z, v.mask, v.meta, "oooo");
    const hillslope_relief_result relief = compute_hillslope_relief(flow, v.z, v.mask, v.channel);

    // A ramp indexed by flat recovers each pixel's channel link exactly.
    std::vector<float> ramp(n);
    for (std::size_t i = 0; i < n; ++i) ramp[i] = static_cast<float>(i);
    const auto all = tag_hillslopes_from_channel(v.meta, relief.channel_flat, v.mask, ramp, v.channel);
    bool exact = true;
    std::int64_t n_linked = 0, n_unlinked = 0, n_on_channel = 0;
    for (std::size_t i = 0; i < n; ++i)
    {
        if (v.mask[i] == 0u) { exact = exact && all.tag.mask[i] == 0u; continue; }
        const std::int64_t link = relief.channel_flat[i];
        if (link < 0) { ++n_unlinked; exact = exact && all.tag.mask[i] == 0u; continue; }
        if (link == static_cast<std::int64_t>(i)) ++n_on_channel; else ++n_linked;
        exact = exact && all.tag.mask[i] == 1u && all.tag.data[i] == ramp[link];
    }
    check("tag of a flat-index ramp recovers the channel link", exact);
    check("tag counts", all.n_hillslope_tagged == n_linked && all.n_channel_tagged == n_on_channel &&
                            all.n_no_channel == n_unlinked && all.n_no_value == 0 && all.n_excluded == 0,
          std::to_string(all.n_hillslope_tagged) + " hillslope");
    check("the fixture has hillslope pixels", n_linked > 100);

    // Partial coverage: alternate channel rows have no value.
    std::vector<std::uint8_t> partial(n, 0u);
    for (std::size_t i = 0; i < n; ++i)
        partial[i] = (v.channel[i] != 0u && (i / v.nc) % 2 == 0) ? 1u : 0u;
    const auto part = tag_hillslopes_from_channel(v.meta, relief.channel_flat, v.mask, ramp, partial);
    std::int64_t want_tagged = 0, want_missing = 0;
    for (std::size_t i = 0; i < n; ++i)
    {
        const std::int64_t link = relief.channel_flat[i];
        if (v.mask[i] == 0u || link < 0 || link == static_cast<std::int64_t>(i)) continue;
        if (partial[link]) ++want_tagged; else ++want_missing;
    }
    check("partial coverage is counted", part.n_hillslope_tagged == want_tagged &&
                                             part.n_no_value == want_missing && want_missing > 0);

    // A non-finite channel value is missing, never copied.
    std::vector<float> with_nan = ramp;
    with_nan[v.left_outlet] = std::numeric_limits<float>::quiet_NaN();
    const auto nan_tag = tag_hillslopes_from_channel(v.meta, relief.channel_flat, v.mask, with_nan, v.channel);
    check("NaN at a channel pixel gives nodata", nan_tag.tag.mask[v.left_outlet] == 0u);

    bool bad_size = false, off_grid = false;
    try { (void)tag_hillslopes_from_channel(v.meta, relief.channel_flat, v.mask, std::vector<float>(3), v.channel); }
    catch (const std::invalid_argument&) { bad_size = true; }
    std::vector<std::int64_t> wrong = relief.channel_flat;
    wrong[v.left_outlet] = static_cast<std::int64_t>(n);
    try { (void)tag_hillslopes_from_channel(v.meta, wrong, v.mask, ramp, v.channel); }
    catch (const std::invalid_argument&) { off_grid = true; }
    check("a wrong array size is refused", bad_size);
    check("a channel link off the grid is refused", off_grid);

    // Basin selection.
    hillslope_tag_basin_options by_height;
    by_height.max_outlet_elevation_m = 50.0;
    by_height.minimum_contributing_pixels = 5;
    const auto low = select_hillslope_tag_basins(flow, v.z, by_height);
    check("elevation limit keeps only the low basin",
          low.outlets == std::vector<std::int64_t>{v.left_outlet} && low.n_excluded_elevation > 0,
          std::to_string(low.outlets.size()) + " outlets");
    check("elevation limit: the high basin is outside the mask",
          low.in_selected_basin[v.left_outlet] == 1u && low.in_selected_basin[v.right_outlet] == 0u);
    hillslope_tag_basin_options too_big;
    too_big.outlet_flats = {v.left_outlet};
    too_big.minimum_contributing_pixels = static_cast<int>(n);
    const auto small = select_hillslope_tag_basins(flow, v.z, too_big);
    check("a basin below the channel threshold is excluded and counted",
          small.outlets.empty() && small.n_excluded_too_small == 1);

    hillslope_tag_basin_options by_list;
    by_list.outlet_flats = {v.right_outlet, v.right_outlet};
    by_list.minimum_contributing_pixels = 5;
    const auto listed = select_hillslope_tag_basins(flow, v.z, by_list);
    check("outlet list selects that basin once", listed.outlets == std::vector<std::int64_t>{v.right_outlet} &&
                                                     listed.n_candidates == 1);

    hillslope_tag_basin_options nested;
    const std::int64_t inner = (v.nr - 8) * v.nc + 5;
    nested.outlet_flats = {inner, v.left_outlet};
    const auto nest = select_hillslope_tag_basins(flow, v.z, nested);
    check("a nested outlet gives way to the outer one",
          nest.outlets == std::vector<std::int64_t>{v.left_outlet} && nest.n_excluded_nested == 1);

    bool neither = false, both = false, bad_outlet = false;
    try { (void)select_hillslope_tag_basins(flow, v.z, {}); }
    catch (const std::invalid_argument&) { neither = true; }
    hillslope_tag_basin_options two = by_list;
    two.max_outlet_elevation_m = 50.0;
    try { (void)select_hillslope_tag_basins(flow, v.z, two); }
    catch (const std::invalid_argument&) { both = true; }
    hillslope_tag_basin_options nodata_outlet;
    nodata_outlet.outlet_flats = {0};
    try { (void)select_hillslope_tag_basins(flow, v.z, nodata_outlet); }
    catch (const std::invalid_argument&) { bad_outlet = true; }
    check("neither basin setting is refused", neither);
    check("both basin settings are refused", both);
    check("an outlet on nodata is refused", bad_outlet);

    // A basin touching the grid edge is excluded, but a complete basin nested in
    // it survives: the edge test comes before the nesting test.
    {
        const int nr = 20, nc = 9, c0 = 4;
        const lsdtt_raster_meta meta = make_meta(nr, nc, 10.0);
        std::vector<float> z(nr * nc);
        const std::vector<std::uint8_t> full(nr * nc, 1u);
        for (int r = 0; r < nr; ++r)
            for (int c = 0; c < nc; ++c)
                z[r * nc + c] = static_cast<float>(5.0 * std::abs(c - c0) + 2.0 * (nr - 1 - r) + 10.0);
        const lsdtt_flowinfo open = make_flow(z, full, meta, "oooo");
        hillslope_tag_basin_options edge;
        edge.max_outlet_elevation_m = 1000.0;
        const auto none = select_hillslope_tag_basins(open, z, edge);
        check("basins on the grid edge are excluded",
              none.outlets.empty() && none.n_excluded_grid_edge > 0,
              std::to_string(none.n_excluded_grid_edge) + " excluded");

        std::vector<float> framed = z;
        std::vector<std::uint8_t> frame = full;
        for (int r = 0; r < nr; ++r)
            for (int c = 0; c < nc; ++c)
                if (r == 0 && c != c0)
                { framed[r * nc + c] = -9999.f; frame[r * nc + c] = 0u; }
        const lsdtt_flowinfo partly = make_flow(framed, frame, meta, "oooo");
        hillslope_tag_basin_options pair;
        const std::int64_t outer = (nr - 1) * nc + c0;   // basin reaches row 0 of the grid
        const std::int64_t middle = 10 * nc + c0;
        pair.outlet_flats = {outer, middle};
        const auto mixed = select_hillslope_tag_basins(partly, framed, pair);
        check("edge basin dropped before nesting is tested",
              mixed.n_excluded_grid_edge == 2 && mixed.n_excluded_nested == 0 && mixed.outlets.empty(),
              "both reach the top row");
    }

    // Chi and k_sn: excluded hillslopes, and agreement with the hilltop join.
    const lsdtt_raster_cpu surface{std::vector<float>(v.z), std::vector<std::uint8_t>(v.mask), v.meta};
    hillslope_chi_tag_options chi_opt;
    chi_opt.minimum_contributing_pixels = 5;
    const auto tags = compute_hillslope_chi_tags(flow, surface, relief.channel_flat, low, chi_opt);
    check("chi was fitted", tags.n_fitted_nodes > 10 && tags.n_segments >= 1,
          std::to_string(tags.n_fitted_nodes) + " nodes");

    std::int64_t want_excluded = 0;
    bool high_is_nodata = true, chi_matches = true;
    for (std::size_t i = 0; i < n; ++i)
    {
        const std::int64_t link = relief.channel_flat[i];
        if (v.mask[i] == 0u || link < 0) continue;
        if (!low.in_selected_basin[link])
        {
            if (link != static_cast<std::int64_t>(i)) ++want_excluded;
            high_is_nodata = high_is_nodata && tags.chi.tag.mask[i] == 0u && tags.k_sn.tag.mask[i] == 0u;
        }
        else if (tags.chi.tag.mask[i])
            chi_matches = chi_matches && tags.chi.tag.mask[link] &&
                          tags.chi.tag.data[i] == tags.chi.tag.data[link];
    }
    check("excluded hillslopes are nodata", high_is_nodata);
    check("excluded hillslopes are counted",
          tags.chi.n_excluded == want_excluded && tags.k_sn.n_excluded == want_excluded && want_excluded > 0,
          std::to_string(tags.chi.n_excluded));
    check("a hillslope's chi equals its channel pixel's", chi_matches);
    check("chi is zero at the outlet", tags.chi.tag.mask[v.left_outlet] && tags.chi.tag.data[v.left_outlet] == 0.f);

    const auto networks = get_chi_from_outlets_zero_base_threshold(flow, low.outlets, 5, chi_opt.chi);
    const auto fits = fit_chi_networks(networks, surface, {});
    const auto map = build_channel_property_map(export_chi_network_fit_nodes(fits),
                                                export_chi_network_fit_segments(fits), nullptr, nullptr);
    lsdtt_hillslope_metrics_result rows;
    std::vector<std::size_t> pixels;
    for (std::size_t i = 0; i < n; ++i)
    {
        const std::int64_t link = relief.channel_flat[i];
        if (v.mask[i] == 0u || link < 0 || link == static_cast<std::int64_t>(i) || !map.count(link)) continue;
        lsdtt_hilltop_metrics_row row;
        row.status = hilltop_trace_status::reached_channel;
        row.channel_row = static_cast<int>(link / v.nc);
        row.channel_col = static_cast<int>(link % v.nc);
        rows.rows.push_back(row);
        pixels.push_back(i);
    }
    (void)attach_channel_properties(rows, map, v.nc);
    bool join_matches = !pixels.empty();
    for (std::size_t k = 0; k < pixels.size(); ++k)
        join_matches = join_matches && tags.k_sn.tag.mask[pixels[k]] &&
                       tags.k_sn.tag.data[pixels[k]] == static_cast<float>(rows.rows[k].channel_k_sn);
    check("k_sn tag equals the hilltop join for the same channel pixel", join_matches,
          std::to_string(pixels.size()) + " pixels");

    hillslope_tag_basins empty_selection;
    empty_selection.in_selected_basin.assign(n, 0u);
    const auto no_tags = compute_hillslope_chi_tags(flow, surface, relief.channel_flat, empty_selection, chi_opt);
    check("no selected basins gives all-nodata tags",
          no_tags.chi.n_hillslope_tagged == 0 && no_tags.chi.n_channel_tagged == 0 &&
              no_tags.chi.n_excluded == all.n_hillslope_tagged);
}

// A Y-shaped network: a main stem in column 8 from the junction at row 10, and
// two diagonal tributaries above it.
struct y_network
{
    lsdtt_raster_meta meta;
    std::vector<float> z;
    std::vector<std::uint8_t> mask, channel;
    int nr, nc;
    std::int64_t flat(int r, int c) const { return static_cast<std::int64_t>(r) * nc + c; }
    lsdtt_xy centre(std::int64_t f) const { return cell_center_from_row_col(meta, static_cast<int>(f / nc), static_cast<int>(f % nc)); }
};

y_network make_y_network()
{
    const int nr = 22, nc = 17;
    y_network y{make_meta(nr, nc, 10.0), std::vector<float>(nr * nc),
                std::vector<std::uint8_t>(nr * nc, 1u), std::vector<std::uint8_t>(nr * nc, 0u), nr, nc};
    y.meta.set_crs_from_epsg(32633);
    for (int r = 0; r < nr; ++r)
        for (int c = 0; c < nc; ++c)
        {
            int offset;
            if (r >= 10) offset = std::abs(c - 8);
            else
            {
                const int spread = 10 - r;
                offset = std::min(std::abs(c - (8 - spread)), std::abs(c - (8 + spread)));
            }
            y.z[r * nc + c] = static_cast<float>(100.0 + 5.0 * offset + 1.0 * (nr - 1 - r));
            if (offset == 0 && r >= 2) y.channel[r * nc + c] = 1u;
        }
    return y;
}

channel_attribute_table exact_table(const std::vector<std::pair<std::int64_t, double>>& records)
{
    channel_attribute_table t;
    t.columns = {"v"};
    t.values.resize(1);
    t.n_records = static_cast<std::int64_t>(records.size());
    for (std::size_t k = 0; k < records.size(); ++k)
    {
        t.values[0].push_back(records[k].second);
        t.sample_record.push_back(static_cast<std::int64_t>(k));
        t.sample_flat.push_back(records[k].first);
        t.sample_x.push_back(0.0);
        t.sample_y.push_back(0.0);
    }
    return t;
}

channel_attribute_table point_table(double x, double y, double v)
{
    channel_attribute_table t;
    t.columns = {"v"};
    t.values = {{v}};
    t.n_records = 1;
    t.sample_record = {0};
    t.sample_flat = {-1};
    t.sample_x = {x};
    t.sample_y = {y};
    return t;
}

void write_text(const std::string& path, const std::string& text)
{
    std::ofstream f(path);
    f << text;
}

lsdtt_vector_schema value_schema(bool with_pixel)
{
    lsdtt_vector_schema schema;
    if (with_pixel)
    {
        schema.field_names = {"row", "col", "flat"};
        schema.field_types = {lsdtt_vector_field_type::int64_type, lsdtt_vector_field_type::int64_type,
                              lsdtt_vector_field_type::int64_type};
    }
    schema.field_names.push_back("v");
    schema.field_types.push_back(lsdtt_vector_field_type::double_type);
    return schema;
}

bool same_column(const channel_attribute_column& a, const channel_attribute_column& b)
{
    return a.mask == b.mask && a.values == b.values;
}

template <typename F>
bool refuses(F&& f)
{
    try { f(); }
    catch (const std::invalid_argument&) { return true; }
    return false;
}

void test_attributes()
{
    std::cout << "\n--- imported attributes ---\n";
    const y_network y = make_y_network();
    const std::size_t n = y.z.size();
    lsdtt_raster_meta flow_meta = y.meta;
    const lsdtt_flowinfo flow = make_flow(y.z, y.mask, flow_meta, "oooo");
    const auto join = [&](const channel_attribute_table& t, double limit) {
        return join_channel_attributes(flow, y.mask, y.channel, t, limit);
    };

    std::vector<std::int64_t> channel_pixels;
    for (std::size_t i = 0; i < n; ++i)
        if (y.channel[i]) channel_pixels.push_back(static_cast<std::int64_t>(i));
    const auto value_of = [](std::int64_t f) { return 0.5 * static_cast<double>(f) + 1.0; };

    const std::string dir = (std::filesystem::temp_directory_path() /
                             ("lsdtt_test_channel_attributes_" + std::to_string(std::random_device{}()))).string();
    std::filesystem::create_directories(dir);

    // Exact CSV, as written by write_channel_network_points.
    {
        std::string text = "row,col,flat,latitude,longitude,v\n";
        std::string coords_text = "latitude,longitude,v\n";
        lsdtt_raster_meta wgs84;
        wgs84.set_crs_from_epsg(4326);
        const lsdtt_coordinate_transform to_ll(y.meta.crs_wkt, wgs84.crs_wkt);
        for (const auto f : channel_pixels)
        {
            const lsdtt_xy c = y.centre(f);
            const lsdtt_xy ll = to_ll.transform_xy(c.x, c.y);
            char buf[160];
            std::snprintf(buf, sizeof(buf), "%d,%d,%lld,%.10f,%.10f,%.3f\n", static_cast<int>(f / y.nc),
                          static_cast<int>(f % y.nc), static_cast<long long>(f), ll.y, ll.x, value_of(f));
            text += buf;
            std::snprintf(buf, sizeof(buf), "%.10f,%.10f,%.3f\n", ll.y, ll.x, value_of(f));
            coords_text += buf;
        }
        write_text(dir + "/exact.csv", text);
        write_text(dir + "/coords.csv", coords_text);

        const auto exact = join(read_channel_attribute_table(dir + "/exact.csv", y.meta, {"v"}), 0.0);
        bool ok = exact.columns[0].n_direct == static_cast<std::int64_t>(channel_pixels.size());
        for (const auto f : channel_pixels)
            ok = ok && exact.columns[0].mask[f] && exact.columns[0].values[f] == static_cast<float>(value_of(f));
        check("exact CSV join gives every channel pixel its value", ok);
        check("exact join counts", exact.n_records_matched == exact.n_records && exact.n_records_dropped == 0);

        const auto coords = join(read_channel_attribute_table(dir + "/coords.csv", y.meta, {"v"}), 5.0);
        check("latitude and longitude give the same tags as the exact join", same_column(coords.columns[0], exact.columns[0]));

        // FlatGeobuf with pixel columns, in the DEM CRS.
        {
            lsdtt_data_loader w;
            w.begin_vector_write(dir + "/exact.fgb", "a", value_schema(true), y.meta.crs_wkt, vector_output_format::flatgeobuf);
            for (const auto f : channel_pixels)
            {
                const lsdtt_xy c = y.centre(f);
                w.write_point(c.x, c.y, {std::int64_t(f / y.nc), std::int64_t(f % y.nc), std::int64_t(f), value_of(f)});
            }
            w.finish_vector_write();
        }
        const auto exact_fgb = join(read_channel_attribute_table(dir + "/exact.fgb", y.meta, {"v"}), 0.0);
        check("exact FlatGeobuf equals exact CSV", same_column(exact_fgb.columns[0], exact.columns[0]));

        // FlatGeobuf with coordinates only, in another CRS.
        {
            lsdtt_raster_meta wgs;
            wgs.set_crs_from_epsg(4326);
            const lsdtt_coordinate_transform to_wgs(y.meta.crs_wkt, wgs.crs_wkt);
            lsdtt_data_loader w;
            w.begin_vector_write(dir + "/wgs.fgb", "a", value_schema(false), wgs.crs_wkt, vector_output_format::flatgeobuf);
            for (const auto f : channel_pixels)
            {
                const lsdtt_xy c = y.centre(f);
                const lsdtt_xy p = to_wgs.transform_xy(c.x, c.y);
                w.write_point(p.x, p.y, {value_of(f)});
            }
            w.finish_vector_write();
        }
        const auto reprojected = join(read_channel_attribute_table(dir + "/wgs.fgb", y.meta, {"v"}), 5.0);
        check("a file in another CRS gives the same tags", same_column(reprojected.columns[0], exact.columns[0]));

        // Hillslopes take the imported values through the D8 link.
        const auto relief = compute_hillslope_relief(flow, y.z, y.mask, y.channel);
        const auto tagged = tag_hillslopes_from_channel(y.meta, relief.channel_flat, y.mask,
                                                        exact.columns[0].values, exact.columns[0].mask);
        bool through_link = tagged.n_hillslope_tagged > 100;
        for (std::size_t i = 0; i < n; ++i)
            if (relief.channel_flat[i] >= 0 && !y.channel[i])
                through_link = through_link && tagged.tag.mask[i] &&
                               tagged.tag.data[i] == static_cast<float>(value_of(relief.channel_flat[i]));
        check("hillslopes take the imported value of their channel pixel", through_link);
    }

    // Snap distance: 30 m east of a main stem pixel centre, on a row with no
    // nearer channel pixel.
    {
        const std::int64_t f = y.flat(15, 8);
        const lsdtt_xy c = y.centre(f);
        const auto inside = join(point_table(c.x + 30.0, c.y, 4.0), 30.0);
        const auto outside = join(point_table(c.x + 30.0, c.y, 4.0), 29.99);
        check("a point exactly at the snap distance joins", inside.columns[0].n_direct == 1 && inside.columns[0].mask[f]);
        check("a snap of three cells is reported", inside.n_records_snapped_beyond_one_cell == 1);
        check("a point just beyond it is dropped and counted", outside.columns[0].n_direct == 0 && outside.n_records_dropped == 1);
        const auto nearby = join(point_table(c.x + 5.0, c.y, 4.0), 30.0);
        check("a snap within one cell is not reported as moved", nearby.n_records_snapped_beyond_one_cell == 0);
        const auto off_grid = join(point_table(c.x + 95.0, c.y, 4.0), 95.0);
        check("a point off the grid can still snap", off_grid.columns[0].mask[f] == 1u);
    }

    // Medians.
    {
        const std::int64_t f = y.flat(15, 8);
        const double nan = std::numeric_limits<double>::quiet_NaN();
        const auto odd = join(exact_table({{f, 9.0}, {f, 1.0}, {f, 4.0}}), 0.0);
        const auto even = join(exact_table({{f, 9.0}, {f, 1.0}, {f, 4.0}, {f, 6.0}}), 0.0);
        const auto with_nan = join(exact_table({{f, 9.0}, {f, nan}, {f, 1.0}}), 0.0);
        const auto all_nan = join(exact_table({{f, nan}}), 0.0);
        check("odd records give the median", odd.columns[0].values[f] == 4.f && odd.columns[0].n_multiple_records == 1);
        check("even records give the mean of the middle two", even.columns[0].values[f] == 5.f);
        check("NaN is ignored in the median", with_nan.columns[0].values[f] == 5.f);
        check("a NaN value stays nodata", all_nan.columns[0].mask[f] == 0u && all_nan.columns[0].n_direct == 0);
        const auto off = join(exact_table({{y.flat(15, 2), 3.0}}), 100.0);
        check("an exact record off the network is dropped and counted",
              off.n_exact_samples_off_channel == 1 && off.n_records_dropped == 1);
    }

    // Lines.
    {
        const lsdtt_xy top = y.centre(y.flat(11, 8)), bottom = y.centre(y.flat(21, 8));
        lsdtt_data_loader w;
        w.begin_vector_write(dir + "/lines.fgb", "a", value_schema(false), y.meta.crs_wkt, vector_output_format::flatgeobuf);
        w.write_line({top, bottom}, {3.0});
        // Runs over the same pixels many times: it must still count once per pixel.
        w.write_line({top, bottom, top, bottom, top}, {9.0});
        w.finish_vector_write();
        const auto both = join(read_channel_attribute_table(dir + "/lines.fgb", y.meta, {"v"}), 1.0);
        bool ok = true;
        for (int r = 11; r <= 21; ++r)
            ok = ok && both.columns[0].mask[y.flat(r, 8)] && both.columns[0].values[y.flat(r, 8)] == 6.f;
        check("overlapping lines give the median, one vote per line per pixel", ok);

        lsdtt_data_loader w1;
        w1.begin_vector_write(dir + "/line.fgb", "a", value_schema(false), y.meta.crs_wkt, vector_output_format::flatgeobuf);
        w1.write_line({top, bottom}, {3.0});
        w1.finish_vector_write();
        const auto single = join(read_channel_attribute_table(dir + "/line.fgb", y.meta, {"v"}), 1.0);
        bool ok1 = single.n_records == 1 && single.n_records_matched == 1;
        for (int r = 11; r <= 21; ++r)
            ok1 = ok1 && single.columns[0].values[y.flat(r, 8)] == 3.f;
        check("every channel pixel under a line takes its value", ok1);

        lsdtt_data_loader w2;
        w2.begin_vector_write(dir + "/multi.fgb", "a", value_schema(false), y.meta.crs_wkt, vector_output_format::flatgeobuf);
        w2.write_multiline({{y.centre(y.flat(11, 8)), y.centre(y.flat(13, 8))},
                            {y.centre(y.flat(17, 8)), y.centre(y.flat(19, 8))}}, {5.0});
        w2.finish_vector_write();
        const auto multi = join(read_channel_attribute_table(dir + "/multi.fgb", y.meta, {"v"}), 1.0);
        check("a multipart line is one record", multi.n_records == 1 && multi.columns[0].values[y.flat(12, 8)] == 5.f &&
                                                   multi.columns[0].values[y.flat(18, 8)] == 5.f);
    }

    // Empty and NaN attributes in a file.
    {
        const std::int64_t a = y.flat(15, 8), b = y.flat(16, 8);
        write_text(dir + "/blank.csv", "row,col,flat,v\n15,8," + std::to_string(a) + ",\n16,8," +
                                           std::to_string(b) + ",nan\n");
        const auto blank = join(read_channel_attribute_table(dir + "/blank.csv", y.meta, {"v"}), 0.0);
        check("empty and NaN values stay nodata, not 0",
              blank.columns[0].mask[a] == 0u && blank.columns[0].mask[b] == 0u && blank.columns[0].n_direct == 0);
    }

    // Gap filling.
    {
        const auto gap = join(exact_table({{y.flat(12, 8), 2.0}, {y.flat(18, 8), 8.0}}), 100.0);
        const auto& g = gap.columns[0];
        check("a gap takes the nearest record on its link",
              g.values[y.flat(13, 8)] == 2.f && g.values[y.flat(14, 8)] == 2.f &&
                  g.values[y.flat(16, 8)] == 8.f && g.values[y.flat(17, 8)] == 8.f);
        check("equidistant records give their median", g.values[y.flat(15, 8)] == 5.f);
        check("filled pixels are counted", g.n_direct == 2 && g.n_filled == 10, std::to_string(g.n_filled));

        const auto limited = join(exact_table({{y.flat(12, 8), 2.0}}), 20.0);
        const auto& l = limited.columns[0];
        check("a fill reaches the limit", l.mask[y.flat(14, 8)] == 1u && l.mask[y.flat(10, 8)] == 1u);
        check("a gap beyond the limit stays nodata", l.mask[y.flat(15, 8)] == 0u && l.mask[y.flat(21, 8)] == 0u);

        const auto main_only = join(exact_table({{y.flat(10, 8), 7.0}}), 1000.0);
        bool no_crossing = true;
        for (const auto f : channel_pixels)
            if (f / y.nc < 10) no_crossing = no_crossing && main_only.columns[0].mask[f] == 0u;
        check("a main-stem record does not fill a tributary", no_crossing);
        check("the main stem is filled", main_only.columns[0].mask[y.flat(21, 8)] == 1u);

        const auto trib = join(exact_table({{y.flat(5, 3), 6.0}}), 1000.0);
        check("a tributary record fills its own tributary", trib.columns[0].mask[y.flat(9, 7)] == 1u);
        check("a tributary record stops at the junction",
              trib.columns[0].mask[y.flat(10, 8)] == 0u && trib.columns[0].mask[y.flat(9, 9)] == 0u);

        const auto diag = join(exact_table({{y.flat(5, 3), 6.0}}), 15.0);
        check("diagonal steps count as sqrt(2) cells",
              diag.columns[0].mask[y.flat(6, 4)] == 1u && diag.columns[0].mask[y.flat(7, 5)] == 0u);
    }

    // Refusals.
    {
        const std::string a = std::to_string(y.flat(15, 8));
        write_text(dir + "/bad_row.csv", "row,col,flat,v\n99,8," + a + ",1\n");
        write_text(dir + "/bad_flat.csv", "row,col,flat,v\n15,8,3,1\n");
        write_text(dir + "/text.csv", "row,col,flat,v\n15,8," + a + ",wide\n");
        write_text(dir + "/no_location.csv", "v\n1\n");
        const auto read = [&](const std::string& name, std::vector<std::string> columns = {"v"}) {
            (void)read_channel_attribute_table(dir + "/" + name, y.meta, columns);
        };
        check("a row off this grid is refused", refuses([&] { read("bad_row.csv"); }));
        check("a flat disagreeing with row and col is refused", refuses([&] { read("bad_flat.csv"); }));
        check("a text value is refused", refuses([&] { read("text.csv"); }));
        check("a table with no location is refused", refuses([&] { read("no_location.csv"); }));
        check("a missing column is refused (CSV)", refuses([&] { read("exact.csv", {"w"}); }));
        check("a missing column is refused (vector)", refuses([&] { read("exact.fgb", {"w"}); }));
        check("no columns is refused", refuses([&] { read("exact.csv", {}); }));
        check("a repeated column is refused", refuses([&] { read("exact.csv", {"v", "v"}); }));

        lsdtt_raster_meta shifted = y.meta;
        shifted.x_min += 500.0;
        check("a table on a shifted grid with the same size is refused",
              refuses([&] { (void)read_channel_attribute_table(dir + "/exact.fgb", shifted, {"v"}); }));

        lsdtt_vector_schema text_schema;
        text_schema.field_names = {"v"};
        text_schema.field_types = {lsdtt_vector_field_type::string_type};
        lsdtt_data_loader w;
        w.begin_vector_write(dir + "/text.fgb", "a", text_schema, y.meta.crs_wkt, vector_output_format::flatgeobuf);
        const lsdtt_xy c = y.centre(y.flat(15, 8));
        w.write_point(c.x, c.y, {std::string("wide")});
        w.finish_vector_write();
        check("a text field is refused", refuses([&] { read("text.fgb"); }));

        lsdtt_data_loader wp;
        wp.begin_vector_write(dir + "/polygon.fgb", "a", value_schema(false), y.meta.crs_wkt, vector_output_format::flatgeobuf);
        wp.write_polygon({{c.x, c.y}, {c.x + 10, c.y}, {c.x + 10, c.y + 10}, {c.x, c.y}}, {1.0});
        wp.finish_vector_write();
        check("polygon geometry is refused", refuses([&] { read("polygon.fgb"); }));

        check("a negative snap distance is refused", refuses([&] { (void)join(exact_table({}), -1.0); }));
    }

    std::error_code ignored;
    std::filesystem::remove_all(dir, ignored);
}

} // namespace

int main()
{
    try
    {
        test_valleys();
        test_paths();
        test_surface();
        test_divides();
        test_rays();
        test_tags();
        test_attributes();
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
