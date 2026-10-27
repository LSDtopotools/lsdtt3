// test_riser_finder.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
//
// -----------------------------------------------------------------------------
// find_risers (lsdtt_riser_finder_cpu), fluvial-surface methods plan phase 6.
//
// Fixtures, all 10 m cells on the synthetic georeference:
//
//   staircase   generate_terrace_staircase_surface: four 300 m treads, 8 m
//               risers, tread cross-fall 0.002, so the TRUE step across a
//               riser is 8 - 0.002 * 300 = 7.4 m. Six risers (three a side).
//
//   tanh risers two straight risers (8 m each, at x = 505 and 1005 m),
//               z = H/2 (1 + tanh((x - x_j) / 5 m)), searched inside a mask
//               with a road corridor (rows 76-84) cut out: the risers stop
//               and restart, with the DEM intact either side - what gap
//               closing is for. Gaps made in the DEM itself are no use here:
//               the staircase generator's gaps are graded ramps with vertical
//               edges up- and down-valley, so the riser runs into a real
//               across-valley scarp and has no free end to close; widening
//               the tanh to 100 m does not break the line at all; and fading
//               the riser height to zero puts a real across-valley step at
//               each end of the fade, which is correctly found as a riser.
//               Spur pruning trims ~10 px off every free end, so the opening
//               the gap closer sees is ~20 rows wider than the corridor; it
//               stays inside the 300 m search radius.
//
//   ramp        a plane falling 0.05 down-valley and 0.02 across: no riser.
//   noisy flat  a flat with 0.3 m deterministic noise: no riser.
//
// Usage: lsdtt-test-riser-finder        (prints PASS/FAIL per check and
//                                         RESULT: PASS|FAIL)
// -----------------------------------------------------------------------------

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "lsdtt_raster_cpu.hpp"
#include "lsdtt_raster_synthetic_surfaces_cpu.hpp"
#include "lsdtt_riser_finder_cpu.hpp"

namespace
{

int g_failures = 0;

void check(bool ok, const std::string& what, const std::string& detail = "")
{
    std::cout << (ok ? "PASS  " : "FAIL  ") << what;
    if (!detail.empty()) std::cout << "  [" << detail << "]";
    std::cout << "\n";
    if (!ok) ++g_failures;
}

lsdtt::lsdtt_raster_cpu make_raster(std::vector<float> z, int nr, int nc)
{
    lsdtt::lsdtt_raster_meta meta;
    meta.n_rows        = nr;
    meta.n_cols        = nc;
    meta.cell_size     = 10.0;
    meta.no_data_value = -9999.0;
    meta.set_default_synthetic_georeference();
    lsdtt::lsdtt_raster_cpu out;
    std::vector<std::uint8_t> valid(z.size(), 1);
    out.set_from_vectors(std::move(z), std::move(valid), meta);
    return out;
}

// Deterministic noise in [-1, 1] per cell.
double hash_noise(int r, int c)
{
    std::uint32_t h = static_cast<std::uint32_t>(r) * 73856093u ^
                      static_cast<std::uint32_t>(c) * 19349663u;
    h ^= h >> 13; h *= 0x5bd1e995u; h ^= h >> 15;
    return (static_cast<double>(h % 20001u) / 10000.0) - 1.0;
}

constexpr int    k_deg_rows = 160, k_deg_cols = 150;
constexpr int    k_deg_riser_col[2] = {50, 100};
constexpr double k_deg_riser_h = 8.0;
constexpr int    k_gap_row_lo = 76, k_gap_row_hi = 84;   // masked out

lsdtt::lsdtt_raster_cpu make_degraded()
{
    std::vector<float> z(static_cast<std::size_t>(k_deg_rows) * k_deg_cols);
    for (int r = 0; r < k_deg_rows; ++r)
        for (int c = 0; c < k_deg_cols; ++c)
        {
            const double x = c * 10.0;
            double v = 100.0 - 0.002 * r * 10.0 + 0.002 * x;
            for (int j = 0; j < 2; ++j)
                v += 0.5 * k_deg_riser_h *
                     (1.0 + std::tanh((x - k_deg_riser_col[j] * 10.0 - 5.0) / 5.0));
            z[static_cast<std::size_t>(r) * k_deg_cols + c] = static_cast<float>(v);
        }
    return make_raster(std::move(z), k_deg_rows, k_deg_cols);
}

// The search mask: everything but a road corridor across the valley.
lsdtt::lsdtt_raster_cpu make_road_mask()
{
    std::vector<float> m(static_cast<std::size_t>(k_deg_rows) * k_deg_cols, 1.0f);
    for (int r = k_gap_row_lo; r <= k_gap_row_hi; ++r)
        for (int c = 0; c < k_deg_cols; ++c)
            m[static_cast<std::size_t>(r) * k_deg_cols + c] = 0.0f;
    return make_raster(std::move(m), k_deg_rows, k_deg_cols);
}

lsdtt::riser_request default_request(lsdtt::riser_gap_closing_choice gc)
{
    lsdtt::riser_request rq;
    rq.gap_closing = gc;
    rq.skeleton.elevation_connect.search_radius_metres = 300.0;
    return rq;
}

// Rows in which the kept centreline has a pixel within 3 columns of col.
int rows_covered(const lsdtt::riser_result& res, int col)
{
    const auto& m  = res.centreline.get_meta();
    const auto& cl = res.centreline.get_data();
    int n = 0;
    for (int r = 0; r < m.n_rows; ++r)
        for (int c = std::max(0, col - 3); c <= std::min(m.n_cols - 1, col + 3); ++c)
            if (cl[static_cast<std::size_t>(r) * m.n_cols + c] > 0.5f) { ++n; break; }
    return n;
}

// Gap-closing pixels more than 4 columns from every riser: a bridge between
// two risers at different heights would put pixels there.
int gap_pixels_off_risers(const lsdtt::riser_result& res)
{
    const auto& m  = res.gap_closing_mask.get_meta();
    const auto& gm = res.gap_closing_mask.get_data();
    int n = 0;
    for (std::size_t i = 0; i < gm.size(); ++i)
    {
        if (!(gm[i] > 0.5f)) continue;
        const int c = static_cast<int>(i % m.n_cols);
        bool near = false;
        for (int rc : k_deg_riser_col) if (std::abs(c - rc) <= 4) near = true;
        if (!near) ++n;
    }
    return n;
}

void report(const std::string& tag, const lsdtt::riser_result& res)
{
    std::cout << "  " << tag << ": found " << res.n_segments_found
              << ", dropped short " << res.n_dropped_short
              << ", dropped low step " << res.n_dropped_low_step
              << ", kept " << res.segments.size() << "\n";
    for (const auto& s : res.segments)
        std::cout << "    segment " << s.id << ": " << s.n_pixels << " px, "
                  << s.length_m << " m, step " << s.step_height_m
                  << " m, face gradient " << s.face_gradient << "\n";
}

} // namespace

int main()
{
    using lsdtt::riser_gap_closing_choice;

    // ---- staircase: six full-length risers, step 7.4 m ---------------------
    {
        const auto st = lsdtt::generate_terrace_staircase_surface(
            200, 300, 10.0, 4, 300.0, 8.0, 0.002, 0.0, 0.002, 0, 120.0, 0.0,
            100.0, 0.5, 1, 1);
        const auto res = lsdtt::find_risers(st.dem, lsdtt::lsdtt_raster_cpu(),
                                            default_request(riser_gap_closing_choice::none));
        report("staircase", res);
        int n_full = 0, n_step_ok = 0;
        for (const auto& s : res.segments)
        {
            if (s.n_pixels >= 180) ++n_full;
            if (std::abs(s.step_height_m - 7.4) <= 0.74) ++n_step_ok;
        }
        check(res.segments.size() == 6 && n_full == 6,
              "staircase: the six risers, each one full-length segment",
              std::to_string(res.segments.size()) + " kept, " +
              std::to_string(n_full) + " >= 180 px");
        check(n_step_ok == static_cast<int>(res.segments.size()) && n_step_ok > 0,
              "staircase: every step height within 10 % of the true 7.4 m");

        auto rq = default_request(riser_gap_closing_choice::none);
        rq.n_threads = 1;
        const auto r1 = lsdtt::find_risers(st.dem, lsdtt::lsdtt_raster_cpu(), rq);
        rq.n_threads = 8;
        const auto r8 = lsdtt::find_risers(st.dem, lsdtt::lsdtt_raster_cpu(), rq);
        check(r1.centreline.get_data() == r8.centreline.get_data() &&
              r1.segments.size() == r8.segments.size(),
              "staircase: identical at 1 and 8 threads");
    }

    // ---- degraded risers: gap closing ----------------------------------------
    {
        const auto dem = make_degraded();
        const auto res0 = lsdtt::find_risers(dem, lsdtt::lsdtt_raster_cpu(),
                                             default_request(riser_gap_closing_choice::none));
        report("tanh risers, no mask", res0);
        check(res0.segments.size() == 2 &&
              rows_covered(res0, k_deg_riser_col[0]) >= k_deg_rows - 20 &&
              rows_covered(res0, k_deg_riser_col[1]) >= k_deg_rows - 20,
              "tanh risers, no mask: two segments, each near full length");

        const auto mask = make_road_mask();
        const int gap_rows = k_gap_row_hi - k_gap_row_lo + 1;
        for (const auto gc : {riser_gap_closing_choice::none,
                              riser_gap_closing_choice::astar})
        {
            const std::string name = gc == riser_gap_closing_choice::none ? "none" : "astar";
            const auto res = lsdtt::find_risers(dem, mask, default_request(gc));
            report("road mask, " + name, res);
            const int cov0 = rows_covered(res, k_deg_riser_col[0]);
            const int cov1 = rows_covered(res, k_deg_riser_col[1]);
            const std::string detail = "rows covered " + std::to_string(cov0) +
                                       ", " + std::to_string(cov1);
            if (gc == riser_gap_closing_choice::none)
            {
                check(cov0 <= k_deg_rows - gap_rows && cov1 <= k_deg_rows - gap_rows,
                      "gap closing none: the road gap stays open", detail);
            }
            else
            {
                check(res.segments.size() == 2 &&
                      cov0 >= k_deg_rows - 20 && cov1 >= k_deg_rows - 20,
                      "gap closing " + name + ": each riser one segment across the gap",
                      std::to_string(res.segments.size()) + " segments, " + detail);
                check(gap_pixels_off_risers(res) == 0,
                      "gap closing " + name + ": never bridges between the two risers",
                      std::to_string(gap_pixels_off_risers(res)) + " pixels off-riser");
            }
            for (const auto& s : res.segments)
                if (std::abs(s.step_height_m - 7.4) > 1.5 && s.n_pixels > 20)
                {
                    check(false, "road mask " + name + ": kept segment step near 7-8 m",
                          std::to_string(s.step_height_m));
                    break;
                }
        }
    }

    // ---- riser-free surfaces --------------------------------------------------
    {
        std::vector<float> ramp(120 * 120), noisy(120 * 120);
        for (int r = 0; r < 120; ++r)
            for (int c = 0; c < 120; ++c)
            {
                ramp[r * 120 + c]  = static_cast<float>(500.0 - 0.05 * r * 10.0 + 0.02 * c * 10.0);
                noisy[r * 120 + c] = static_cast<float>(200.0 + 0.3 * hash_noise(r, c));
            }
        const auto rr = lsdtt::find_risers(make_raster(ramp, 120, 120),
                                           lsdtt::lsdtt_raster_cpu(),
                                           default_request(riser_gap_closing_choice::none));
        report("ramp", rr);
        check(rr.segments.empty(), "planar ramp: no riser segments");
        const auto rn = lsdtt::find_risers(make_raster(noisy, 120, 120),
                                           lsdtt::lsdtt_raster_cpu(),
                                           default_request(riser_gap_closing_choice::none));
        report("noisy flat", rn);
        check(rn.segments.empty(), "noisy flat (0.3 m): no riser segments");
    }

    // ---- argument checks ------------------------------------------------------
    {
        bool threw = false;
        try { (void)lsdtt::parse_riser_gap_closing_choice("contour_trace"); }
        catch (const std::invalid_argument&) { threw = true; }
        check(threw, "unknown gap closing name (contour_trace, not offered) is refused");
        threw = false;
        lsdtt::riser_request rq = default_request(riser_gap_closing_choice::astar);
        rq.skeleton.elevation_connect.search_radius_metres = 0.0;
        try { (void)lsdtt::find_risers(make_degraded(), lsdtt::lsdtt_raster_cpu(), rq); }
        catch (const std::invalid_argument&) { threw = true; }
        check(threw, "astar with no search radius is refused");
    }

    std::cout << "RESULT: " << (g_failures == 0 ? "PASS" : "FAIL") << "\n";
    return g_failures == 0 ? 0 : 1;
}
