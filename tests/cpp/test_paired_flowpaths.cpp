// test_paired_flowpaths.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
//
// -----------------------------------------------------------------------------
// Tests for lsdtt_paired_flowpaths_cpu -- the Goren & Shelef (2024) lengthwise
// asymmetry dL of paired flow paths.
//
// A NOTE ON WHAT IS BEING TESTED. These tests drive the computation with
// SYNTHETIC ridgelines: a straight line of divide-grid nodes laid over a
// convergent V-valley DEM. Those lines are not real drainage divides, and they
// are not meant to be -- pairing the module with real interfluvial divides from
// lsdtt_basin_outline_cpu is a separate concern. What a straight synthetic line
// buys is a set of flanking pixel pairs whose flow paths are exactly known, so
// the arithmetic, the datum handling, the rejection bookkeeping and the
// invariances can be pinned down without the divide extraction in the way.
//
// The load-bearing checks are the two that no amount of plausible-looking
// output would satisfy by accident:
//
//   * scale invariance -- dL is a ratio of lengths, so tripling the cell size
//     and the relief must reproduce every dL bit-for-bit;
//   * the junction datum -- L must be measured from the junction where the two
//     paths rejoin, not from the outlet. Two paths that rejoin share all of
//     their downstream length, and that shared part cancels in the numerator of
//     dL but not the denominator, so an outlet datum systematically drags dL
//     towards zero and makes every network look simple. The test asserts the
//     junction-datum median is the larger of the two.
// -----------------------------------------------------------------------------

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

#include "lsdtt_chi_extraction_cpu.hpp"
#include "lsdtt_flowinfo_cpu.hpp"
#include "lsdtt_paired_flowpaths_cpu.hpp"
#include "lsdtt_raster_cpu.hpp"

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

constexpr int k_n_rows = 41;
constexpr int k_n_cols = 41;
constexpr int k_valley_col = 20;

// -----------------------------------------------------------------------------
// A V-shaped valley running south along k_valley_col, tilted south.
//
// The cross-valley gradient (0.6 per cell) beats the down-valley one (0.4), so
// every cell steps diagonally towards the valley axis until it reaches it, then
// runs south. The whole grid therefore drains through the single cell
// (n_rows-1, k_valley_col): every pair of pixels rejoins, which is what makes
// the junction datum testable.
//
// `scale` multiplies both the cell size and the relief, so the surface is
// geometrically similar -- all lengths scale, all ratios do not.
// -----------------------------------------------------------------------------
lsdtt::lsdtt_raster_cpu make_v_valley_dem(double scale)
{
    lsdtt::lsdtt_raster_meta meta;
    meta.n_rows = k_n_rows;
    meta.n_cols = k_n_cols;
    meta.cell_size = 10.0 * scale;
    meta.no_data_value = -9999.0;
    meta.set_default_synthetic_georeference();

    const std::size_t n =
        static_cast<std::size_t>(k_n_rows) * static_cast<std::size_t>(k_n_cols);
    std::vector<float> data(n, 0.0f);
    std::vector<std::uint8_t> mask(n, 1);

    for (int row = 0; row < k_n_rows; ++row)
    {
        for (int col = 0; col < k_n_cols; ++col)
        {
            const double z =
                1000.0 - 0.4 * static_cast<double>(row) +
                0.6 * std::abs(static_cast<double>(col - k_valley_col));

            data[static_cast<std::size_t>(row) * k_n_cols + col] =
                static_cast<float>(z * scale);
        }
    }

    return lsdtt::lsdtt_raster_cpu(std::move(data), std::move(mask), meta);
}

// Two basins back to back: a ridge along k_valley_col, west flank draining to
// the west edge and east flank to the east edge. Nothing crosses.
lsdtt::lsdtt_raster_cpu make_two_basin_dem()
{
    lsdtt::lsdtt_raster_meta meta;
    meta.n_rows = k_n_rows;
    meta.n_cols = k_n_cols;
    meta.cell_size = 10.0;
    meta.no_data_value = -9999.0;
    meta.set_default_synthetic_georeference();

    const std::size_t n =
        static_cast<std::size_t>(k_n_rows) * static_cast<std::size_t>(k_n_cols);
    std::vector<float> data(n, 0.0f);
    std::vector<std::uint8_t> mask(n, 1);

    for (int row = 0; row < k_n_rows; ++row)
    {
        for (int col = 0; col < k_n_cols; ++col)
        {
            const double z =
                1000.0 - 0.6 * std::abs(static_cast<double>(col - k_valley_col));
            data[static_cast<std::size_t>(row) * k_n_cols + col] =
                static_cast<float>(z);
        }
    }

    return lsdtt::lsdtt_raster_cpu(std::move(data), std::move(mask), meta);
}

lsdtt::lsdtt_flowinfo make_flow(const lsdtt::lsdtt_raster_cpu& dem)
{
    lsdtt::lsdtt_flowinfo_options fopts;
    fopts.boundary_conditions = "oooo";
    fopts.build_topological_order = true;
    fopts.build_s_vector = true;

    return lsdtt::lsdtt_flowinfo(
        dem.get_data_array(), dem.get_mask_array(), dem.get_meta(), fopts);
}

// A straight run of divide-grid nodes down column `col`, rows 0..n-1. Divide
// node (r, c) is the top-left corner of DEM pixel (r, c), so this line separates
// DEM columns col-1 and col.
lsdtt::lsdtt_ridgeline_cpu make_straight_ridgeline(
    int col, int n_nodes, const lsdtt::lsdtt_raster_meta& dem_meta)
{
    std::vector<lsdtt::lsdtt_pixel_rc> nodes;
    nodes.reserve(static_cast<std::size_t>(n_nodes));
    for (int r = 0; r < n_nodes; ++r) nodes.push_back({r, col});

    // The divide grid sits on pixel corners: one more row and column than the
    // DEM, shifted half a cell.
    lsdtt::lsdtt_raster_meta dmeta = dem_meta;
    dmeta.n_rows = dem_meta.n_rows + 1;
    dmeta.n_cols = dem_meta.n_cols + 1;

    return lsdtt::lsdtt_ridgeline_cpu(std::move(nodes), dmeta);
}

struct scene
{
    lsdtt::lsdtt_raster_cpu dem;
    std::vector<float> fd_at_channel;
    std::vector<float> chi_at_channel;
    std::vector<float> fd_from_outlet;
    std::vector<float> chi_from_outlet;
    std::vector<std::uint8_t> channel_mask;
};

scene build_scene(const lsdtt::lsdtt_flowinfo& flow,
                  lsdtt::lsdtt_raster_cpu dem)
{
    scene s{std::move(dem), {}, {}, {}, {}, {}};

    const lsdtt::lsdtt_raster_result fd = flow.get_distance_from_outlet_raster();
    const lsdtt::lsdtt_raster_cpu chi = lsdtt::get_chi_raster(flow);

    s.fd_from_outlet = fd.data;
    {
        const lsdtt::lsdtt_array2d_f& c = chi.get_data_array();
        s.chi_from_outlet.assign(c.data(), c.data() + fd.data.size());
    }

    // Channels are cells with at least 15 contributing pixels.
    const lsdtt::lsdtt_raster_result ncp =
        flow.get_n_contributing_pixels_raster();

    s.channel_mask.assign(ncp.data.size(), 0);
    for (std::size_t i = 0; i < ncp.data.size(); ++i)
    {
        if (std::isfinite(ncp.data[i]) && ncp.data[i] >= 15.0f)
            s.channel_mask[i] = 1;
    }

    const lsdtt::lsdtt_array2d_u8& m = s.dem.get_mask_array();
    const std::vector<std::uint8_t> valid(m.data(), m.data() + ncp.data.size());

    s.fd_at_channel =
        lsdtt::sample_at_nearest_channel(flow, s.fd_from_outlet, valid,
                                         s.channel_mask).data;
    s.chi_at_channel =
        lsdtt::sample_at_nearest_channel(flow, s.chi_from_outlet, valid,
                                         s.channel_mask).data;

    return s;
}

lsdtt::lsdtt_paired_flowpath_result run(
    const lsdtt::lsdtt_flowinfo& flow,
    const scene& s,
    const lsdtt::lsdtt_paired_flowpath_options& opts)
{
    std::vector<lsdtt::lsdtt_ridgeline_cpu> rls;
    rls.push_back(make_straight_ridgeline(31, 21, s.dem.get_meta()));
    rls.push_back(make_straight_ridgeline(11, 21, s.dem.get_meta()));

    return lsdtt::compute_paired_flowpath_statistics(
        flow, rls, s.fd_at_channel, s.chi_at_channel, s.fd_from_outlet,
        s.chi_from_outlet, s.channel_mask, opts);
}

double median_of_copy(std::vector<double> v)
{
    if (v.empty()) return std::numeric_limits<double>::quiet_NaN();
    std::sort(v.begin(), v.end());
    const std::size_t m = v.size() / 2;
    return (v.size() % 2 == 1) ? v[m] : 0.5 * (v[m - 1] + v[m]);
}

} // namespace


int main()
{
    std::cout << "=== paired flow path (lengthwise asymmetry) tests ===\n\n";

    // -------------------------------------------------------------------------
    // PF1 - the common downstream node
    // -------------------------------------------------------------------------
    const lsdtt::lsdtt_raster_cpu dem = make_v_valley_dem(1.0);
    const lsdtt::lsdtt_flowinfo flow = make_flow(dem);

    const std::int64_t a = 5 * k_n_cols + 5;
    const std::int64_t b = 5 * k_n_cols + 35;

    const std::int64_t cj = lsdtt::find_common_downstream_node(flow, a, b);
    check("PF1a common downstream node exists in a single-outlet DEM", cj >= 0);
    check("PF1b common downstream node is symmetric in its arguments",
          cj == lsdtt::find_common_downstream_node(flow, b, a));
    check("PF1c a node is its own common downstream node",
          lsdtt::find_common_downstream_node(flow, a, a) == a);

    // The junction must lie on both flow paths, and be the FIRST node of b's
    // path that does so. Verified here by walking the receivers independently.
    {
        auto on_path = [&](std::int64_t start, std::int64_t target) -> bool
        {
            std::int64_t x = start;
            for (int g = 0; g < k_n_rows * k_n_cols + 1; ++g)
            {
                if (x == target) return true;
                const std::int64_t r = flow.receiver_of(x);
                if (r < 0 || r == x) return false;
                x = r;
            }
            return false;
        };
        check("PF1d junction lies on both flow paths",
              on_path(a, cj) && on_path(b, cj));

        std::int64_t first = -1;
        std::int64_t y = b;
        for (int g = 0; g < k_n_rows * k_n_cols + 1; ++g)
        {
            if (on_path(a, y)) { first = y; break; }
            const std::int64_t r = flow.receiver_of(y);
            if (r < 0 || r == y) break;
            y = r;
        }
        check("PF1e junction is the FIRST shared node, not merely a shared one",
              first == cj, "expected " + std::to_string(first) +
                           " got " + std::to_string(cj));
    }

    // Two basins that never meet.
    {
        const lsdtt::lsdtt_raster_cpu two = make_two_basin_dem();
        const lsdtt::lsdtt_flowinfo tflow = make_flow(two);
        const std::int64_t west = 20 * k_n_cols + 5;
        const std::int64_t east = 20 * k_n_cols + 35;
        check("PF1f paths in separate basins report no common node",
              lsdtt::find_common_downstream_node(tflow, west, east) < 0);
    }

    // -------------------------------------------------------------------------
    // PF2 - the statistics themselves
    // -------------------------------------------------------------------------
    const scene s = build_scene(flow, make_v_valley_dem(1.0));

    lsdtt::lsdtt_paired_flowpath_options opts;
    opts.threshold_divide_distance_m = 0.0;

    const lsdtt::lsdtt_paired_flowpath_result base = run(flow, s, opts);

    check("PF2a pairs were produced", base.n_pairs > 0,
          "n_pairs=" + std::to_string(base.n_pairs));
    check("PF2b every pair rejoins (single-outlet DEM)",
          base.n_rejected_no_junction == 0,
          "n_rejected_no_junction=" +
              std::to_string(base.n_rejected_no_junction));
    check("PF2c pairs never exceed edges examined",
          base.n_pairs <= base.n_edges_examined);

    {
        bool in_range = true;
        bool chi_non_negative = true;
        for (const lsdtt::lsdtt_paired_flowpath& p : base.pairs)
        {
            if (!std::isfinite(p.delta_L) || p.delta_L < 0.0 || p.delta_L > 2.0)
                in_range = false;
            // Both channels are upstream of the junction, so junction-relative
            // chi cannot be negative.
            if (p.chi_a < -1e-6 || p.chi_b < -1e-6) chi_non_negative = false;
        }
        check("PF2d every delta_L lies in [0, 2]", in_range);
        check("PF2e junction-relative chi is non-negative", chi_non_negative);
        check("PF2f median delta_L is finite",
              std::isfinite(base.median_delta_L));
    }

    // -------------------------------------------------------------------------
    // PF3 - scale invariance. dL is a ratio of lengths, so a geometrically
    // similar landscape must give identical values.
    // -------------------------------------------------------------------------
    {
        const lsdtt::lsdtt_raster_cpu big_dem = make_v_valley_dem(3.0);
        const lsdtt::lsdtt_flowinfo big_flow = make_flow(big_dem);
        const scene big = build_scene(big_flow, make_v_valley_dem(3.0));
        const lsdtt::lsdtt_paired_flowpath_result r = run(big_flow, big, opts);

        check("PF3a same pair count at triple scale",
              r.n_pairs == base.n_pairs,
              std::to_string(r.n_pairs) + " vs " + std::to_string(base.n_pairs));

        bool same = (r.n_pairs == base.n_pairs);
        double worst = 0.0;
        if (same)
        {
            for (std::size_t i = 0; i < r.pairs.size(); ++i)
            {
                const double d =
                    std::abs(r.pairs[i].delta_L - base.pairs[i].delta_L);
                worst = std::max(worst, d);
            }
            // The distance and chi rasters are float, so agreement is bounded by
            // single precision, not by the algorithm.
            same = worst < 1e-5;
        }
        check("PF3b every delta_L is unchanged at triple scale", same,
              "max difference " + std::to_string(worst));

        // Lengths themselves must scale by exactly three -- proof the invariance
        // above is a real cancellation and not both runs being degenerate.
        bool scaled = base.n_pairs > 0;
        for (std::size_t i = 0; i < r.pairs.size() && i < base.pairs.size(); ++i)
        {
            if (base.pairs[i].L_a > 1.0 &&
                std::abs(r.pairs[i].L_a - 3.0 * base.pairs[i].L_a) > 1e-3)
                scaled = false;
        }
        check("PF3c the underlying path lengths did scale by three", scaled);
    }

    // -------------------------------------------------------------------------
    // PF4 - the junction datum. Measuring L from the outlet instead leaves the
    // shared downstream length in both terms, which cancels in the numerator
    // but not the denominator and pushes dL towards zero.
    // -------------------------------------------------------------------------
    {
        std::vector<double> junction_datum, outlet_datum;
        for (const lsdtt::lsdtt_paired_flowpath& p : base.pairs)
        {
            const double fa =
                static_cast<double>(s.fd_at_channel[
                    static_cast<std::size_t>(p.pixel_a)]);
            const double fb =
                static_cast<double>(s.fd_at_channel[
                    static_cast<std::size_t>(p.pixel_b)]);
            if (fa + fb <= 0.0) continue;

            junction_datum.push_back(p.delta_L);
            outlet_datum.push_back(2.0 * std::abs(fa - fb) / (fa + fb));
        }

        const double m_j = median_of_copy(junction_datum);
        const double m_o = median_of_copy(outlet_datum);

        check("PF4a junction datum gives a larger median delta_L than the "
              "outlet datum",
              std::isfinite(m_j) && std::isfinite(m_o) && m_j > m_o,
              "junction=" + std::to_string(m_j) +
                  " outlet=" + std::to_string(m_o));
    }

    // -------------------------------------------------------------------------
    // PF5 - the filters
    // -------------------------------------------------------------------------
    {
        lsdtt::lsdtt_paired_flowpath_options o = opts;
        o.only_take_channel_heads = true;
        const lsdtt::lsdtt_paired_flowpath_result r = run(flow, s, o);

        check("PF5a the channel-head filter never adds pairs",
              r.n_pairs <= base.n_pairs,
              std::to_string(r.n_pairs) + " vs " + std::to_string(base.n_pairs));
        check("PF5b pairs dropped by the head filter are counted",
              r.n_pairs + r.n_rejected_not_head == base.n_pairs,
              std::to_string(r.n_pairs) + " + " +
                  std::to_string(r.n_rejected_not_head) + " vs " +
                  std::to_string(base.n_pairs));
    }

    {
        lsdtt::lsdtt_paired_flowpath_options o = opts;
        o.unique_pairs = true;
        const lsdtt::lsdtt_paired_flowpath_result r = run(flow, s, o);
        check("PF5c de-duplication never adds pairs", r.n_pairs <= base.n_pairs);
        check("PF5d duplicates dropped are counted",
              r.n_pairs + r.n_rejected_duplicate == base.n_pairs);
    }

    {
        lsdtt::lsdtt_paired_flowpath_options o = opts;
        o.threshold_divide_distance_m = 1.0e9;   // longer than any ridgeline
        const lsdtt::lsdtt_paired_flowpath_result r = run(flow, s, o);
        check("PF5e an unreachable divide-distance threshold rejects everything",
              r.n_pairs == 0 && r.n_rejected_distance > 0,
              "n_pairs=" + std::to_string(r.n_pairs));
    }

    std::cout << "\n" << (g_n_checks - g_n_failed) << " / " << g_n_checks
              << " checks passed\n";

    return (g_n_failed == 0) ? 0 : 1;
}
