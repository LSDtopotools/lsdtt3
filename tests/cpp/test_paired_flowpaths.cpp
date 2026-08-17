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
#include <filesystem>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "../../src/core/lsdtt_data_loader.hpp"
#include "lsdtt_chi_extraction_cpu.hpp"
#include "lsdtt_flowinfo_cpu.hpp"
#include "lsdtt_junction_basin_manager_cpu.hpp"
#include "lsdtt_junction_network_cpu.hpp"
#include "lsdtt_paired_flowpaths_cpu.hpp"
#include "lsdtt_raster_cpu.hpp"
#include "lsdtt_raster_flow_preprocessing_cpu.hpp"
#include "lsdtt_raster_synthetic_surfaces_cpu.hpp"

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
    std::vector<std::uint8_t> junction_mask;
    std::vector<float> sl_at_channel;
};

// A junction mask built straight from the flow graph and the channel mask:
// a channel pixel is a junction if it is a channel source (no channel donor),
// a confluence (two or more channel donors) or an outlet. That is the same rule
// lsdtt_junction_network_cpu::identify_junctions uses, reproduced here so the
// straight-link tests do not need a junction network object.
std::vector<std::uint8_t> make_junction_mask(
    const lsdtt::lsdtt_flowinfo& flow,
    const std::vector<std::uint8_t>& channel_mask)
{
    const std::size_t n = channel_mask.size();
    std::vector<int> n_channel_donors(n, 0);

    for (std::size_t i = 0; i < n; ++i)
    {
        if (!channel_mask[i]) continue;
        const std::int64_t flat = static_cast<std::int64_t>(i);
        const std::int64_t rec = flow.receiver_of(flat);
        if (rec >= 0 && rec != flat && channel_mask[static_cast<std::size_t>(rec)])
            ++n_channel_donors[static_cast<std::size_t>(rec)];
    }

    std::vector<std::uint8_t> jmask(n, 0);
    for (std::size_t i = 0; i < n; ++i)
    {
        if (!channel_mask[i]) continue;
        const std::int64_t flat = static_cast<std::int64_t>(i);
        const std::int64_t rec = flow.receiver_of(flat);
        const bool is_outlet = (rec < 0 || rec == flat);
        if (n_channel_donors[i] != 1 || is_outlet) jmask[i] = 1;
    }
    return jmask;
}

scene build_scene(const lsdtt::lsdtt_flowinfo& flow,
                  lsdtt::lsdtt_raster_cpu dem)
{
    scene s{std::move(dem), {}, {}, {}, {}, {}, {}, {}};

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

    s.junction_mask = make_junction_mask(flow, s.channel_mask);
    s.sl_at_channel =
        lsdtt::get_straight_link_distance_at_nearest_channel(
            flow, s.channel_mask, s.junction_mask).raster.data;

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
        flow, rls, s.fd_at_channel, s.chi_at_channel, s.sl_at_channel,
        s.fd_from_outlet, s.chi_from_outlet, s.channel_mask, opts);
}

// -----------------------------------------------------------------------------
// CLI fixture. Written when argv[1] is a directory, following the convention
// lsdtt-test-divide-migration uses, so a CLI driver test has a PROJECTED DEM to
// run lsdtt-ridgeline-analysis on. tests/data/tiny30_dem.tif is geographic and
// the CLI refuses it outright.
//
// A fractal surface rather than an analytic one: the pair source and the
// straight-link switch both need a real branching network with sibling groups
// in it, which a V-valley does not have.
// -----------------------------------------------------------------------------
void write_cli_fixture(const std::string& out_dir)
{
    const lsdtt::lsdtt_raster_cpu dem =
        lsdtt::generate_diamond_square_surface(
            120, 120, 30.0, /*seed*/ 1, /*roughness*/ 0.6,
            /*initial_scale*/ 1.0, /*target_mean*/ 1000.0,
            /*target_stddev*/ 40.0);

    std::filesystem::create_directories(out_dir);
    const std::string dem_path = out_dir + "/pfcli_dem.tif";

    lsdtt::lsdtt_data_loader writer;
    writer.write_full(dem.get_data(), dem.get_meta(), dem_path,
                      lsdtt::raster_output_format::gtiff);
    std::cout << "Wrote CLI fixture DEM: " << dem_path << "\n";
}

double median_of_copy(std::vector<double> v)
{
    if (v.empty()) return std::numeric_limits<double>::quiet_NaN();
    std::sort(v.begin(), v.end());
    const std::size_t m = v.size() / 2;
    return (v.size() % 2 == 1) ? v[m] : 0.5 * (v[m - 1] + v[m]);
}

} // namespace


int main(int argc, char* argv[])
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

    // -------------------------------------------------------------------------
    // SL - the straight-link variant
    //
    // The point of the variant is that it does NOT count channel sinuosity, so
    // the tests have to separate it from plain flow distance. On the V-valley
    // every channel link is straight (flanks run a straight diagonal to the
    // valley, the trunk runs straight south), so there SL and FD must AGREE --
    // which pins the accumulation but cannot tell the two quantities apart.
    // SL4 uses a deliberately sinuous channel to do that.
    // -------------------------------------------------------------------------
    {
        const lsdtt::straight_link_distance_result sl =
            lsdtt::get_straight_link_distance_at_nearest_channel(
                flow, s.channel_mask, make_junction_mask(flow, s.channel_mask));

        check("SL1 no channel node terminates at baselevel without a junction",
              sl.n_channel_nodes_terminated_at_baselevel == 0,
              "count=" +
                  std::to_string(sl.n_channel_nodes_terminated_at_baselevel));

        // Straight-line between two points can never exceed the path length
        // between them, so this must hold on ANY surface. It is the defining
        // inequality and catches a sign or accumulation error immediately.
        bool sl_le_fd = true;
        double worst_excess = 0.0;
        for (std::size_t i = 0; i < s.channel_mask.size(); ++i)
        {
            if (!s.channel_mask[i]) continue;
            const double a = static_cast<double>(sl.raster.data[i]);
            const double b = static_cast<double>(s.fd_at_channel[i]);
            if (!std::isfinite(a) || !std::isfinite(b)) continue;
            if (a > b + 1e-3) { sl_le_fd = false; worst_excess =
                std::max(worst_excess, a - b); }
        }
        check("SL2 straight-link distance never exceeds flow distance",
              sl_le_fd, "worst excess " + std::to_string(worst_excess) + " m");

        // Every link on this DEM is straight, so the two must coincide.
        bool agree = true;
        double worst_gap = 0.0;
        for (std::size_t i = 0; i < s.channel_mask.size(); ++i)
        {
            if (!s.channel_mask[i]) continue;
            const double a = static_cast<double>(sl.raster.data[i]);
            const double b = static_cast<double>(s.fd_at_channel[i]);
            if (!std::isfinite(a) || !std::isfinite(b)) continue;
            worst_gap = std::max(worst_gap, std::abs(a - b));
            if (std::abs(a - b) > 1e-2) agree = false;
        }
        check("SL3 on a DEM whose channel links are all straight, "
              "straight-link equals flow distance",
              agree, "worst gap " + std::to_string(worst_gap) + " m");
    }

    // SL4 - a sinuous channel, where the two quantities MUST differ.
    //
    // A zig-zag carved into a high plateau: D8 has nowhere to go but along it.
    // The path alternates south and east, so its flow distance is the
    // Manhattan length while the straight-link distance across the single
    // source-to-outlet link is the Euclidean diagonal. Both are known exactly.
    {
        // Each iteration steps south once then east once, so the path makes
        // n_step+1 south steps and n_step east steps and must end on the last
        // row -- otherwise its final cell is an interior pit and the router
        // rejects the DEM.
        constexpr int n_step = 6;
        constexpr int nr = n_step + 2;
        constexpr int nc = n_step + 2;
        constexpr double cs = 10.0;

        lsdtt::lsdtt_raster_meta meta;
        meta.n_rows = nr;
        meta.n_cols = nc;
        meta.cell_size = cs;
        meta.no_data_value = -9999.0;
        meta.set_default_synthetic_georeference();

        const std::size_t n = static_cast<std::size_t>(nr) * nc;
        std::vector<float> z(n, 1000.0f);
        std::vector<std::uint8_t> mask(n, 1);
        std::vector<std::uint8_t> chan(n, 0);

        // The plateau has to drain too, or the router refuses the DEM for
        // having a flat. A gentle south tilt does it, and stays far above the
        // carved channel so the channel still captures its own path.
        for (int rr = 0; rr < nr; ++rr)
            for (int cc = 0; cc < nc; ++cc)
                z[static_cast<std::size_t>(rr) * nc + cc] =
                    static_cast<float>(1000.0 - 0.5 * rr);

        // Carve the zig-zag, descending one metre per step.
        int r = 0, c = 0, step = 0;
        std::vector<std::size_t> path;
        auto put = [&](int rr, int cc)
        {
            const std::size_t idx = static_cast<std::size_t>(rr) * nc + cc;
            z[idx] = static_cast<float>(900.0 - step);
            chan[idx] = 1;
            path.push_back(idx);
            ++step;
        };
        put(r, c);
        for (int k = 0; k < n_step; ++k)
        {
            put(++r, c);          // south
            put(r, ++c);          // east
        }
        put(++r, c);              // final south step onto the outlet row

        lsdtt::lsdtt_raster_cpu dem_z(std::move(z), std::move(mask), meta);
        const lsdtt::lsdtt_flowinfo zflow = make_flow(dem_z);

        const std::vector<std::uint8_t> jmask = make_junction_mask(zflow, chan);
        const lsdtt::straight_link_distance_result sl =
            lsdtt::get_straight_link_distance_at_nearest_channel(
                zflow, chan, jmask);

        const std::size_t src = path.front();
        const std::size_t out = path.back();
        const int src_r = static_cast<int>(src) / nc;
        const int src_c = static_cast<int>(src) % nc;
        const int out_r = static_cast<int>(out) / nc;
        const int out_c = static_cast<int>(out) % nc;

        const double expect_sl =
            std::hypot((out_r - src_r) * cs, (out_c - src_c) * cs);
        const double got_sl = static_cast<double>(sl.raster.data[src]) -
                              static_cast<double>(sl.raster.data[out]);

        // Every carved step is cardinal, so the flow distance is the Manhattan
        // length: (rows travelled + cols travelled) * cell size.
        const double expect_fd =
            (std::abs(out_r - src_r) + std::abs(out_c - src_c)) * cs;

        check("SL4a straight-link across a sinuous link is the Euclidean "
              "distance, not the path length",
              std::abs(got_sl - expect_sl) < 1e-2,
              "expected " + std::to_string(expect_sl) +
                  " got " + std::to_string(got_sl));

        check("SL4b and it is strictly SHORTER than the flow distance, "
              "which is what removing sinuosity means",
              expect_sl < expect_fd - 1e-6,
              "SL=" + std::to_string(expect_sl) +
                  " FD=" + std::to_string(expect_fd));
    }

    // -------------------------------------------------------------------------
    // SIB - the sibling pair source, and THE GUARANTEE
    //
    // dL exists only for a divide whose two flanking basins share a parent. An
    // area or stream-order window over all junctions cannot promise that: it
    // selects on size, and sharing a parent is a topological property. On one
    // snapping-sweep surface a 400-4000 pixel window selected 106 basins and
    // every one of 10007 divide edges was rejected as "paths never rejoin".
    //
    // These checks are on a multi-outlet fractal surface, which is the case that
    // discriminates: the single-outlet V-valley above rejoins everything, so it
    // could not tell a guaranteed pair source from a lucky one. The contrast
    // arms deliberately show the window-style selection producing non-rejoining
    // pairs on the SAME surface where the sibling source produces none.
    // -------------------------------------------------------------------------
    {
        constexpr int  sib_nr = 100;
        constexpr int  sib_nc = 100;
        constexpr double sib_cs = 30.0;

        const lsdtt::lsdtt_raster_cpu ds_raw =
            lsdtt::generate_diamond_square_surface(
                sib_nr, sib_nc, sib_cs, /*seed*/ 1, /*roughness*/ 0.6,
                /*initial_scale*/ 1.0, /*target_mean*/ 1000.0,
                /*target_stddev*/ 40.0);
        const lsdtt::lsdtt_raster_cpu ds = lsdtt::fill_pits(ds_raw, 0.0001f);

        lsdtt::lsdtt_flowinfo_options fopts;
        fopts.boundary_conditions     = "oooo";
        fopts.build_topological_order = true;
        fopts.build_s_vector          = true;
        const lsdtt::lsdtt_flowinfo dflow(
            ds.get_data_array(), ds.get_mask_array(), ds.get_meta(), fopts);

        lsdtt::lsdtt_junction_network_options jopts;
        jopts.build_topological_order = true;
        jopts.build_s_vector          = true;
        lsdtt::lsdtt_junction_network_cpu jn;
        jn.create(dflow, dflow.get_sources_index_threshold(40), jopts);

        const std::vector<lsdtt::lsdtt_sibling_basin_group> groups =
            lsdtt::select_sibling_basin_groups(dflow, jn);
        const std::vector<lsdtt::lsdtt_sibling_basin_pair> sib_pairs =
            lsdtt::get_sibling_basin_pairs(groups);

        check("SIB1a the sibling source finds groups on a fractal surface",
              !groups.empty() && !sib_pairs.empty(),
              "groups=" + std::to_string(groups.size()) +
                  " pairs=" + std::to_string(sib_pairs.size()));

        // The structural claim: a donor of a parent has that parent as its
        // receiver. Everything else rests on this.
        {
            bool ok = !groups.empty();
            for (const auto& g : groups)
            {
                if (g.donor_junctions.size() < 2) ok = false;
                for (std::int32_t d : g.donor_junctions)
                {
                    if (jn.receiver_of_junction(d) != g.parent_junction)
                        ok = false;
                }
            }
            check("SIB1b every donor's receiver junction IS its group's parent",
                  ok);
        }

        // Siblings are disjoint: neither basin can be upslope of the other, so
        // no nesting test is needed downstream.
        {
            bool disjoint = !sib_pairs.empty();
            for (const auto& p : sib_pairs)
            {
                const std::int64_t fa = jn.basin_outlet_flat_of_junction(p.junction_a);
                const std::int64_t fb = jn.basin_outlet_flat_of_junction(p.junction_b);
                if (dflow.is_upslope_of(fa, fb) || dflow.is_upslope_of(fb, fa))
                    disjoint = false;
            }
            check("SIB1c sibling basins are never nested in one another",
                  disjoint);
        }

        // THE GUARANTEE, at the level of the flow graph: both donor outlets
        // drain into the parent junction pixel, so their paths must meet.
        {
            std::size_t n_no_common = 0;
            for (const auto& p : sib_pairs)
            {
                const std::int64_t fa = jn.basin_outlet_flat_of_junction(p.junction_a);
                const std::int64_t fb = jn.basin_outlet_flat_of_junction(p.junction_b);
                if (lsdtt::find_common_downstream_node(dflow, fa, fb) < 0)
                    ++n_no_common;
            }
            check("SIB1d every sibling pair has a common downstream node",
                  n_no_common == 0,
                  "pairs without one: " + std::to_string(n_no_common));
        }

        // CONTRAST. The same surface, selected the old way -- an area window
        // over all junctions -- admits pairs whose paths never meet. Without
        // this arm the check above would pass on any single-outlet DEM and prove
        // nothing about the pair source.
        {
            std::vector<std::int64_t> window_outlets;
            for (std::int32_t j = 0; j < jn.n_junctions(); ++j)
            {
                if (!jn.is_valid_junction(j)) continue;
                const std::int64_t f = jn.basin_outlet_flat_of_junction(j);
                if (f < 0 || !dflow.is_valid_cell(f)) continue;
                const std::int64_t np = dflow.s_vector_size(f);
                if (np >= 60 && np <= 600) window_outlets.push_back(f);
                if (window_outlets.size() >= 60) break;   // keep the O(n^2) arm cheap
            }

            std::size_t n_window_pairs = 0, n_window_no_common = 0;
            for (std::size_t i = 0; i + 1 < window_outlets.size(); ++i)
            {
                for (std::size_t j = i + 1; j < window_outlets.size(); ++j)
                {
                    // The same non-nesting test the CLI's all-pairs scan uses.
                    if (dflow.is_upslope_of(window_outlets[i], window_outlets[j]) ||
                        dflow.is_upslope_of(window_outlets[j], window_outlets[i]))
                        continue;
                    ++n_window_pairs;
                    if (lsdtt::find_common_downstream_node(
                            dflow, window_outlets[i], window_outlets[j]) < 0)
                        ++n_window_no_common;
                }
            }
            check("SIB1e an area window over all junctions DOES admit pairs "
                  "whose paths never meet",
                  n_window_no_common > 0,
                  std::to_string(n_window_no_common) + " of " +
                      std::to_string(n_window_pairs) + " window pairs");
        }

        // ---------------------------------------------------------------------
        // SIB2 - end to end, through the real interfluvial divide extraction and
        // the statistic itself. This is the check that matters: the "paths never
        // rejoin" rejection count must be ZERO.
        // ---------------------------------------------------------------------
        {
            lsdtt::lsdtt_sibling_basin_options sopts;
            sopts.max_groups = 40;          // bounded, so the test stays cheap
            const std::vector<lsdtt::lsdtt_sibling_basin_group> capped =
                lsdtt::select_sibling_basin_groups(dflow, jn, sopts);

            check("SIB2a the group cap is honoured and reproducible",
                  capped.size() == std::min<std::size_t>(40, groups.size()),
                  "capped=" + std::to_string(capped.size()));

            lsdtt::lsdtt_junction_basin_manager_cpu bm(dflow, jn);
            bm.select_junction_basins_from_basin_outlets(
                lsdtt::get_sibling_basin_outlet_flats(jn, capped));

            std::map<std::int32_t, std::size_t> index_of_junction;
            for (std::size_t i = 0; i < bm.size(); ++i)
                index_of_junction[bm.basin_at(i).get_junction_number()] = i;

            std::vector<lsdtt::lsdtt_ridgeline_cpu> sib_divides;
            for (const auto& g : capped)
            {
                for (std::size_t i = 0; i + 1 < g.donor_junctions.size(); ++i)
                {
                    auto it_i = index_of_junction.find(g.donor_junctions[i]);
                    if (it_i == index_of_junction.end()) continue;
                    auto BO = bm.basin_at(it_i->second).get_basin_outline();

                    for (std::size_t j = i + 1; j < g.donor_junctions.size(); ++j)
                    {
                        if (index_of_junction.find(g.donor_junctions[j]) ==
                            index_of_junction.end()) continue;

                        auto rl = BO.get_interfluvial_divide_as_ridgeline(
                            g.donor_junctions[i], g.donor_junctions[j], jn);
                        if (!rl.empty()) sib_divides.push_back(std::move(rl));
                    }
                }
            }

            check("SIB2b sibling pairs yield real interfluvial divides",
                  !sib_divides.empty(),
                  "divides=" + std::to_string(sib_divides.size()));

            const scene ss = build_scene(dflow, lsdtt::fill_pits(ds_raw, 0.0001f));

            lsdtt::lsdtt_paired_flowpath_options po;
            po.threshold_divide_distance_m = 0.0;

            const lsdtt::lsdtt_paired_flowpath_result sr =
                lsdtt::compute_paired_flowpath_statistics(
                    dflow, sib_divides, ss.fd_at_channel, ss.chi_at_channel,
                    ss.sl_at_channel, ss.fd_from_outlet, ss.chi_from_outlet,
                    ss.channel_mask, po);

            // Printed rather than only asserted: a guarantee over an empty or
            // trivial sample is worth nothing, so the sample size has to be
            // visible in the test log.
            std::cout << "   [sibling source: " << capped.size() << " groups, "
                      << sib_divides.size() << " divides, "
                      << sr.n_edges_examined << " edges examined, "
                      << sr.n_pairs << " pairs, "
                      << sr.n_rejected_no_junction << " non-rejoining, "
                      << sr.n_rejected_no_channel << " no-channel]\n";

            check("SIB2c pairs survive from the sibling source", sr.n_pairs > 0,
                  "n_pairs=" + std::to_string(sr.n_pairs) + " of " +
                      std::to_string(sr.n_edges_examined) + " edges");

            check("SIB2d THE GUARANTEE: not one pair is rejected as "
                  "'paths never rejoin'",
                  sr.n_rejected_no_junction == 0,
                  "n_rejected_no_junction=" +
                      std::to_string(sr.n_rejected_no_junction) +
                      " of " + std::to_string(sr.n_edges_examined) + " edges");

            check("SIB2e and the statistic itself is finite",
                  std::isfinite(sr.median_delta_L),
                  "median=" + std::to_string(sr.median_delta_L));
        }
    }

    // The CLI driver fixture, written only when a directory is supplied.
    if (argc > 1)
    {
        try
        {
            write_cli_fixture(argv[1]);
        }
        catch (const std::exception& e)
        {
            check("F1 CLI fixture written", false, e.what());
        }
    }

    std::cout << "\n" << (g_n_checks - g_n_failed) << " / " << g_n_checks
              << " checks passed\n";

    return (g_n_failed == 0) ? 0 : 1;
}
