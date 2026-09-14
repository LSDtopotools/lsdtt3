// test_hillslope_length_statistic.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
//
// -----------------------------------------------------------------------------
// "Mean hillslope length" is ambiguous by a factor of two.
//
// `estimate_hillslope_length_stats` averages the flow distance to the nearest
// channel over EVERY non-channel pixel. On a uniform planar hillslope of length
// L that mean is L/2. But "the average D8 path distance from the ridgetops" --
// the quantity anyone targeting a hillslope length actually means -- is L.
//
// A tool that hits "100 m" against the wrong statistic delivers 200 m
// hillslopes, and nothing about the number would look wrong. So the two are
// computed and reported side by side, and the factor of two between them is
// pinned here on a surface where both answers are known exactly rather than
// approximately.
//
//   H1  the analytic plane: ridge_mean = L and all_pixel_mean = (L+cell)/2,
//       EXACTLY, so the ratio is ~2 -- the whole point of the distinction
//   H2  the ridge set is the D8 no-donor set
//   H3  p90 tracks L, not L/2 (it is the statistic the stiffness prediction
//       uses, so it must not be confused with the mean)
//   H4  a two-sided tent, where the D8 ridge set is TWO cells wide -- the crest
//       drains one way and orphans its neighbour on the other side. Real, and
//       not what you would guess
//   H5  THE DOMAIN EDGE ACTS AS A CHANNEL: a path that leaves the grid is
//       recorded at distance 0. Pre-existing behaviour, load-bearing for
//       anything that searches the channel threshold
//   H6  the selector: names round-trip, and a typo is FATAL rather than
//       quietly picking one of three values that differ by 2x
//
// Usage: lsdtt-test-hillslope-length-statistic
// -----------------------------------------------------------------------------

#include <cmath>
#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

#include "../../src/cpu/lsdtt_landscape_equilibration_cpu.hpp"

namespace
{

int g_checks   = 0;
int g_failures = 0;

void check(bool ok, const std::string& what)
{
    ++g_checks;
    if (!ok) { ++g_failures; std::cout << "  FAIL: " << what << "\n"; }
    else     { std::cout << "  ok  : " << what << "\n"; }
}

void check_close(double got, double expected, double tol, const std::string& what)
{
    ++g_checks;
    if (!(std::fabs(got - expected) <= tol))
    {
        ++g_failures;
        std::cout << "  FAIL: " << what << "\n        got " << got
                  << ", expected " << expected << "\n";
    }
    else std::cout << "  ok  : " << what << "  (" << got << ")\n";
}

lsdtt::lsdtt_raster_meta make_meta(int nr, int nc, double cell)
{
    lsdtt::lsdtt_raster_meta meta;
    meta.n_rows = nr; meta.n_cols = nc; meta.cell_size = cell;
    meta.x_min = 0.0; meta.y_min = 0.0; meta.no_data_value = -9999.0;
    return meta;
}

// A plane tilted purely along +row, so every D8 receiver is due south and every
// flow path is a straight column. That is what makes L and L/2 exact rather
// than approximate: no diagonal steps, no ties.
struct plane_result
{
    lsdtt::hillslope_length_stats stats;
    double L;              // ridge-to-channel distance, metres
    int    n_ridge_rows;
};

plane_result uniform_plane(int nr, int nc, double cell, int channel_row)
{
    auto meta = make_meta(nr, nc, cell);

    std::vector<float>        z(static_cast<std::size_t>(nr) * nc);
    std::vector<std::uint8_t> m(z.size(), 1u);
    for (int r = 0; r < nr; ++r)
        for (int c = 0; c < nc; ++c)
            z[static_cast<std::size_t>(r) * nc + c] =
                static_cast<float>(1000.0 - 2.0 * r);

    const lsdtt::lsdtt_raster_cpu dem(std::move(z), std::move(m), meta);

    lsdtt::lsdtt_flowinfo_options fi;
    fi.boundary_conditions = "oooo";
    lsdtt::lsdtt_flowinfo flow(dem.get_data_array(), dem.get_mask_array(),
                               dem.get_meta(), fi);

    std::vector<std::uint8_t> chan(static_cast<std::size_t>(nr) * nc, 0u);
    for (int c = 0; c < nc; ++c)
        chan[static_cast<std::size_t>(channel_row) * nc + c] = 1u;

    plane_result out;
    out.stats        = lsdtt::estimate_hillslope_length_stats(flow, chan);
    out.L            = channel_row * cell;   // row 0 is the divide
    out.n_ridge_rows = 1;
    return out;
}

}  // namespace


int main()
{
    using namespace lsdtt;

    std::cout << "=====================================================\n"
              << " Hillslope length statistics\n"
              << "=====================================================\n";

    // -------------------------------------------------------------------------
    std::cout << "\nH1: the analytic plane -- ridge_mean = L, all_pixel_mean = (L+cell)/2\n";
    {
        // The channel is the LAST row, so the domain IS the hillslope and no
        // path leaves the grid. (Put the channel mid-domain and the rows below
        // it drain off the edge, which H5 shows is counted as distance 0 -- that
        // would contaminate the very statistic under test.)
        const int    nr = 11, nc = 11, channel_row = 10;
        const double cell = 10.0;
        const auto   p = uniform_plane(nr, nc, cell, channel_row);

        // Rows 0..9 drain south to the channel at row 10. Row 0 is the ridge, at
        // 10 cells = 100 m.
        check_close(p.stats.ridge_mean_m, p.L, 1e-9,
                    "ridge_mean_m == L (distance from the ridgetop)");

        // Distances 100,90,...,10 over rows 0..9: mean 55 = (L + cell)/2. It is
        // exactly L/2 only in the continuum; on a grid the pixel centres are
        // offset by half a cell, and stating that is more useful than tuning the
        // fixture to hide it.
        const double expected_all = (p.L + cell) / 2.0;
        check_close(p.stats.mean_m, expected_all, 1e-9,
                    "mean_m == (L + cell)/2, the discrete form of L/2");

        check_close(p.stats.ridge_mean_m / p.stats.mean_m, 100.0 / 55.0, 1e-9,
                    "the two differ by ~2x -- THE REASON THIS EXISTS");
        check(p.stats.ridge_mean_m > 1.7 * p.stats.mean_m,
              "and the factor is close enough to 2 to mislead in practice");
    }

    // -------------------------------------------------------------------------
    std::cout << "\nH2: the ridge set is the D8 no-donor set\n";
    {
        const auto p = uniform_plane(11, 11, 10.0, 10);
        // Only row 0 has no donor: every other hillslope row receives from the
        // row above. 11 columns.
        check(p.stats.n_ridge_pixels == 11,
              "one ridge row of 11 cells (got " +
              std::to_string(p.stats.n_ridge_pixels) + ")");
        check(p.stats.n_hillslope_pixels == 110,
              "10 hillslope rows x 11 cells (got " +
              std::to_string(p.stats.n_hillslope_pixels) + ")");
        check(p.stats.n_ridge_pixels < p.stats.n_hillslope_pixels,
              "the ridge set is a strict subset");
    }

    // -------------------------------------------------------------------------
    std::cout << "\nH3: p90 tracks L, not L/2\n";
    {
        const auto p = uniform_plane(11, 11, 10.0, 10);
        check(p.stats.p90_m > 0.8 * p.L,
              "p90_m is within 20% of L, not of L/2");
        check(p.stats.p90_m > 1.5 * p.stats.mean_m,
              "and is far from the all-pixel mean");
    }

    // -------------------------------------------------------------------------
    std::cout << "\nH4: two-sided tent -- the D8 ridge band is TWO cells wide\n";
    {
        // A tent: elevation falls away from row 10 in both directions, with
        // channels at rows 0 and 20. The ridge is now INTERIOR, which is the
        // real geometry -- H1's ridge sat on the domain edge.
        const int nr = 21, nc = 11;
        const double cell = 10.0;
        auto meta = make_meta(nr, nc, cell);

        std::vector<float>        z(static_cast<std::size_t>(nr) * nc);
        std::vector<std::uint8_t> m(z.size(), 1u);
        for (int r = 0; r < nr; ++r)
            for (int c = 0; c < nc; ++c)
                z[static_cast<std::size_t>(r) * nc + c] =
                    static_cast<float>(1000.0 - 2.0 * std::abs(r - 10));

        const lsdtt_raster_cpu dem(std::move(z), std::move(m), meta);
        lsdtt_flowinfo_options fi;
        fi.boundary_conditions = "oooo";
        lsdtt_flowinfo flow(dem.get_data_array(), dem.get_mask_array(),
                            dem.get_meta(), fi);

        std::vector<std::uint8_t> chan(static_cast<std::size_t>(nr) * nc, 0u);
        for (int c = 0; c < nc; ++c)
        {
            chan[static_cast<std::size_t>(0) * nc + c]  = 1u;
            chan[static_cast<std::size_t>(20) * nc + c] = 1u;
        }

        const auto st = estimate_hillslope_length_stats(flow, chan);

        // THE D8 RIDGE SET HERE IS TWO ROWS WIDE, NOT ONE, and that is correct
        // rather than a defect. The crest (row 10) has no donors and drains to
        // ONE side; the row on the other side then receives nothing from the
        // crest and nothing from further out (its outward neighbour drains away
        // from it), so it has no donors either. A symmetric divide cannot be
        // one cell wide in D8 -- the crest has to pick a side.
        //
        // The consequence for anyone targeting a hillslope length: ridge_mean is
        // an average over a divide band, so it sits slightly BELOW L, by about
        // half a cell. Here: rows 10 and 11 at 100 m and 90 m -> 95 m.
        check(st.n_ridge_pixels == 22,
              "the D8 ridge band is 2 rows x 11 cells (got " +
              std::to_string(st.n_ridge_pixels) + ")");
        check_close(st.ridge_mean_m, 95.0, 1e-9,
                    "ridge_mean_m = 95 = (L + (L - cell))/2, the divide band");
        check(std::fabs(st.ridge_mean_m - 100.0) <= 10.0,
              "which is within one cell of L -- the bias is bounded by the grid");
        check(st.mean_m < st.ridge_mean_m,
              "the all-pixel mean is still the shorter statistic");
    }

    // -------------------------------------------------------------------------
    std::cout << "\nH5: the domain edge acts as a channel\n";
    {
        // PRE-EXISTING BEHAVIOUR, pinned because it is easy to be caught by and
        // it matters for anything that searches the channel threshold: when a
        // flow path leaves the grid, `estimate_hillslope_length_stats` anchors
        // it at distance 0 (`known_tail = 0.0` on `!is_valid_cell`). The edge is
        // therefore indistinguishable from a channel.
        //
        // Consequence: on a domain where flow exits the edge, the statistics MIX
        // true hillslope distances with distance-to-edge. And if a threshold
        // search ever empties the channel set, this function does not fail --
        // it silently returns distances to the domain edge. Any such search must
        // reject an empty channel set itself.
        const int nr = 21, nc = 11, channel_row = 10;
        const double cell = 10.0;
        const auto p = uniform_plane(nr, nc, cell, channel_row);

        // Rows 11..20 have NO channel below them -- they drain off the bottom
        // edge -- yet they are counted: 20 hillslope rows, not 10.
        check(p.stats.n_hillslope_pixels == 220,
              "rows draining off the edge ARE counted (got " +
              std::to_string(p.stats.n_hillslope_pixels) + " of 220)");

        // And they are counted at distance-to-edge, which happens to mirror the
        // real hillslope here, pulling the mean to exactly 50 rather than 55.
        check_close(p.stats.mean_m, 50.0, 1e-9,
                    "mean_m is pulled to 50 by the off-edge rows (55 without them)");

        // The ridge statistic is immune in this fixture: the off-edge rows all
        // have donors, so none of them joins the ridge set.
        check(p.stats.n_ridge_pixels == 11,
              "the ridge set is unaffected here (got " +
              std::to_string(p.stats.n_ridge_pixels) + ")");
        check_close(p.stats.ridge_mean_m, 100.0, 1e-9,
                    "and ridge_mean_m still reports the true L");
    }

    // -------------------------------------------------------------------------
    std::cout << "\nH6: the selector, and a fatal typo\n";
    {
        check(parse_hillslope_length_statistic("ridge_mean") ==
              hillslope_length_statistic::ridge_mean, "ridge_mean parses");
        check(parse_hillslope_length_statistic("P90") ==
              hillslope_length_statistic::p90, "p90 parses, case-insensitively");
        check(std::string(name_of(hillslope_length_statistic::ridge_mean)) ==
              "ridge_mean", "name round-trips");

        for (const char* bad : {"ridge", "mean", "", "ridge-mean"})
        {
            bool threw = false;
            try { (void)parse_hillslope_length_statistic(bad); }
            catch (const std::exception&) { threw = true; }
            check(threw, std::string("rejects '") + bad +
                  "' rather than picking one of three 2x-apart values");
        }

        const auto p = uniform_plane(11, 11, 10.0, 10);
        check_close(value_of(p.stats, hillslope_length_statistic::ridge_mean),
                    p.stats.ridge_mean_m, 1e-12, "value_of selects ridge_mean");
        check_close(value_of(p.stats, hillslope_length_statistic::all_pixel_mean),
                    p.stats.mean_m, 1e-12, "value_of selects all_pixel_mean");
        check_close(value_of(p.stats, hillslope_length_statistic::p90),
                    p.stats.p90_m, 1e-12, "value_of selects p90");
    }

    std::cout << "\n=====================================================\n"
              << " " << (g_checks - g_failures) << " / " << g_checks
              << " checks passed\n"
              << "=====================================================\n";
    return g_failures == 0 ? 0 : 1;
}
