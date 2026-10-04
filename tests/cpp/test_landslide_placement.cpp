// test_landslide_placement.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
//
// -----------------------------------------------------------------------------
// Phase 1 of the triggered landslide module: placing a scar of a known size on
// a real flow network. Tests T6-T10d:
//
//   T6   no scar pixel has contributing area >= the channel threshold
//   T7   max hostable size stays below the channel threshold
//   T8   annealing blocks re-excavation, and releases at anneal_time
//   T9   scars are contiguous and connected to their initiation pixel
//   T10  the active mask is honoured
//   T10b every scar reaches its target area EXACTLY -- no stalls, no shrinkage
//   T10c every scar is a subset of upslope_nodes(initiation)
//   T10d eligibility rejects hosts with n_up < n_target, and an impossible
//        request is reported rather than silently producing nothing
//
// Usage: test-landslide-placement [output_dir]
// -----------------------------------------------------------------------------

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <vector>

#include "lsdtt_flowinfo_cpu.hpp"
#include "lsdtt_landslide_placement_cpu.hpp"
#include "lsdtt_raster_cpu.hpp"
#include "lsdtt_raster_flow_preprocessing_cpu.hpp"

namespace
{

int g_failures = 0;
int g_checks   = 0;

void check(bool ok, const std::string& what)
{
    ++g_checks;
    if (!ok) { ++g_failures; std::cout << "  FAIL: " << what << "\n"; }
    else      { std::cout << "  ok  : " << what << "\n"; }
}

// A ridge-and-valley surface: a regional tilt with superimposed cross-slope
// sinusoids, so flow converges into parallel valleys and there are real
// hillslopes with a range of contributing areas. Deterministic, so a failure
// is reproducible.
lsdtt::lsdtt_raster_cpu make_ridge_valley(int n_rows, int n_cols,
                                          double cell_size)
{
    lsdtt::lsdtt_raster_meta meta;
    meta.n_rows        = n_rows;
    meta.n_cols        = n_cols;
    meta.cell_size     = cell_size;
    meta.no_data_value = -9999.0;
    meta.set_default_synthetic_georeference();

    lsdtt::lsdtt_array2d_f data(n_rows, n_cols, 0.0f);
    lsdtt::lsdtt_array2d_u8 mask(n_rows, n_cols, 1);

    for (int r = 0; r < n_rows; ++r)
    {
        for (int c = 0; c < n_cols; ++c)
        {
            const double y = static_cast<double>(r);
            const double x = static_cast<double>(c);
            // Tilt down-row, plus valleys every ~16 columns, plus a gentle
            // long-wavelength term so valley heads are not all identical.
            const double z = 220.0 - 1.6 * y
                           + 9.0 * std::cos(2.0 * M_PI * x / 16.0)
                           + 2.5 * std::sin(2.0 * M_PI * (x + 2.0 * y) / 55.0);
            data(r, c) = static_cast<float>(z);
        }
    }
    return lsdtt::lsdtt_raster_cpu(std::move(data), std::move(mask), meta);
}

// 8-connectivity check on a set of flat indices, plus reachability from the
// initiation pixel. A scar that is merely "connected" but detached from where
// we said it started would still be wrong.
bool scar_is_connected(const lsdtt::LandslideScar& scar, int n_cols)
{
    std::set<std::int64_t> remaining(scar.pixels.begin(), scar.pixels.end());
    if (remaining.find(scar.initiation_flat) == remaining.end()) return false;

    std::vector<std::int64_t> stack{scar.initiation_flat};
    remaining.erase(scar.initiation_flat);

    while (!stack.empty())
    {
        const std::int64_t f = stack.back();
        stack.pop_back();
        const int r = static_cast<int>(f / n_cols);
        const int c = static_cast<int>(f % n_cols);
        for (int dr = -1; dr <= 1; ++dr)
        {
            for (int dc = -1; dc <= 1; ++dc)
            {
                if (dr == 0 && dc == 0) continue;
                const std::int64_t nb =
                    static_cast<std::int64_t>(r + dr) * n_cols + (c + dc);
                auto it = remaining.find(nb);
                if (it != remaining.end())
                {
                    remaining.erase(it);
                    stack.push_back(nb);
                }
            }
        }
    }
    return remaining.empty();
}

}  // namespace


int main(int argc, char** argv)
{
    const std::string outdir = (argc > 1) ? argv[1] : ".";
    std::cout << "=== Landslide placement tests (Phase 1) ===\n\n";

    const int    n_rows = 160, n_cols = 160;
    const double cell   = 10.0;

    lsdtt::lsdtt_raster_cpu raw    = make_ridge_valley(n_rows, n_cols, cell);
    lsdtt::lsdtt_raster_cpu filled = lsdtt::fill_pits(raw, 0.0001f);

    lsdtt::lsdtt_flowinfo_options opts;
    opts.boundary_conditions = "oooo";
    lsdtt::lsdtt_flowinfo fi(filled.get_data_array(), filled.get_mask_array(),
                             filled.meta(), opts);

    const lsdtt::lsdtt_raster_meta& meta = filled.meta();

    lsdtt::LandslidePlacementConfig cfg;
    cfg.channel_threshold_pixels = 300;
    cfg.max_host_attempts        = 64;

    lsdtt::LandslidePlacer placer(fi, meta, cfg);
    std::cout << "Landscape: " << n_rows << "x" << n_cols << " at " << cell
              << " m, channel threshold " << cfg.channel_threshold_pixels
              << " px\n"
              << "Largest hostable landslide: " << placer.max_hostable_pixels()
              << " px\n\n";

    // -------------------------------------------------------------------------
    // T7 -- the host list can never offer a landslide at or above the threshold
    // -------------------------------------------------------------------------
    std::cout << "T7: max hostable size below the channel threshold\n";
    check(placer.max_hostable_pixels() < cfg.channel_threshold_pixels,
          "max hostable < channel_threshold_pixels");
    check(placer.max_hostable_pixels() > 20,
          "landscape offers usefully large hosts (sanity on the fixture)");

    // -------------------------------------------------------------------------
    // T10d -- eligibility is monotone, and an impossible request fails cleanly
    // -------------------------------------------------------------------------
    std::cout << "\nT10d: eligible-host counts\n";
    {
        const std::size_t e1  = placer.n_eligible_hosts(1);
        const std::size_t e10 = placer.n_eligible_hosts(10);
        const std::size_t e50 = placer.n_eligible_hosts(50);
        std::cout << "        hosts for 1 px: " << e1 << ", 10 px: " << e10
                  << ", 50 px: " << e50 << "\n";
        check(e1 >= e10 && e10 >= e50, "eligible-host count is non-increasing in size");
        check(e1 > 0, "at least one host for a 1-pixel landslide");

        const int impossible =
            static_cast<int>(placer.max_hostable_pixels()) + 1;
        check(placer.n_eligible_hosts(impossible) == 0,
              "no hosts for a landslide larger than the landscape can hold");

        placer.begin_event();
        lsdtt::LandslideScar scar;
        std::uint64_t rng = 42;
        const bool placed = placer.place_scar_seeded(impossible, rng, scar);
        check(!placed, "impossible request returns false, not a silent empty scar");
        check(placer.diagnostics().n_failed == 1, "the failure is COUNTED");
    }

    // -------------------------------------------------------------------------
    // T6 / T9 / T10b / T10c -- the core invariants, over many placements
    // -------------------------------------------------------------------------
    std::cout << "\nT6/T9/T10b/T10c: invariants over 400 placements\n";
    {
        placer.begin_event();
        std::uint64_t rng = 20260803;

        int n_placed = 0, n_failed = 0;
        bool all_exact = true, all_connected = true, all_upslope = true;
        bool all_hillslope = true, no_reuse = true;
        std::set<std::int64_t> used;

        const std::vector<int> targets = {1, 2, 3, 5, 8, 13, 21, 34, 55, 89};

        for (int i = 0; i < 400; ++i)
        {
            const int n_target = targets[static_cast<std::size_t>(i) % targets.size()];
            lsdtt::LandslideScar scar;
            if (!placer.place_scar_seeded(n_target, rng, scar)) { ++n_failed; continue; }
            ++n_placed;

            // T10b: exactly the requested area, always.
            if (static_cast<int>(scar.pixels.size()) != n_target) all_exact = false;
            if (scar.n_pixels != n_target) all_exact = false;

            // T9: contiguous AND anchored at the initiation pixel.
            if (!scar_is_connected(scar, meta.n_cols)) all_connected = false;

            // T10c: subset of the host's upslope set.
            for (std::int64_t f : scar.pixels)
                if (!fi.is_upslope_of(f, scar.initiation_flat)) all_upslope = false;

            // T6: never a channel pixel. Structural, not enforced in the code.
            for (std::int64_t f : scar.pixels)
                if (static_cast<std::int64_t>(fi.s_vector_size_unchecked(f))
                        >= cfg.channel_threshold_pixels)
                    all_hillslope = false;

            // No pixel excavated twice in one event (forbid_overlap).
            for (std::int64_t f : scar.pixels)
                if (!used.insert(f).second) no_reuse = false;
        }

        std::cout << "        placed " << n_placed << ", failed " << n_failed
                  << ", distinct pixels " << used.size() << "\n";
        check(n_placed > 300, "most placements succeed on an empty landscape");
        check(all_exact,     "T10b: every scar has EXACTLY its target pixel count");
        check(all_connected, "T9: every scar is 8-connected and anchored at its initiation pixel");
        check(all_upslope,   "T10c: every scar pixel is upslope of its initiation pixel");
        check(all_hillslope, "T6: no scar pixel reaches the channel threshold");
        check(no_reuse,      "no pixel is excavated twice within one event");

        placer.finalise_diagnostics();
        std::cout << "        median host fill ratio = "
                  << placer.diagnostics().median_host_fill_ratio << "\n";
    }

    // -------------------------------------------------------------------------
    // Channel-head clustering: large landslides should fill their hosts more
    // completely than small ones. Reported and asserted, because it
    // is an imposed spatial correlation the user should know about.
    // -------------------------------------------------------------------------
    std::cout << "\nHost fill ratio vs landslide size\n";
    {
        const std::string csv = outdir + "/landslide_host_fill.csv";
        FILE* f = std::fopen(csv.c_str(), "w");
        if (f) std::fprintf(f, "n_target_px,n_placed,median_fill_ratio,mean_fill_ratio\n");

        std::cout << "  n_target   placed   median fill ratio\n";
        double small_fill = 0.0, large_fill = 0.0;
        for (int n_target : {2, 5, 20, 60, 120})
        {
            placer.begin_event();
            std::uint64_t rng = 7 + static_cast<std::uint64_t>(n_target);
            std::vector<double> fills;
            for (int i = 0; i < 60; ++i)
            {
                lsdtt::LandslideScar scar;
                if (placer.place_scar_seeded(n_target, rng, scar))
                    fills.push_back(scar.host_fill_ratio);
            }
            if (fills.empty()) { std::cout << "  " << n_target << "  (none)\n"; continue; }
            std::sort(fills.begin(), fills.end());
            const double med = fills[fills.size() / 2];
            double mean = 0.0; for (double v : fills) mean += v; mean /= fills.size();
            std::printf("  %8d %8zu %19.3f\n", n_target, fills.size(), med);
            if (f) std::fprintf(f, "%d,%zu,%.6f,%.6f\n", n_target, fills.size(), med, mean);
            if (n_target == 2)   small_fill = med;
            if (n_target == 120) large_fill = med;
        }
        if (f) std::fclose(f);
        if (large_fill > 0.0)
            check(large_fill > small_fill,
                  "large landslides fill more of their host basin than small ones");
    }

    // -------------------------------------------------------------------------
    // T8 -- annealing
    // -------------------------------------------------------------------------
    std::cout << "\nT8: annealing blocks and then releases\n";
    {
        // time_of_last_failure, all nodata = never failed.
        lsdtt::lsdtt_array2d_f tdata(n_rows, n_cols,
                                     static_cast<float>(meta.no_data_value));
        lsdtt::lsdtt_array2d_u8 tmask(n_rows, n_cols, 1);
        lsdtt::lsdtt_raster_cpu t_fail(std::move(tdata), std::move(tmask), meta);

        // Place one scar and record it as having failed at t = 1000.
        lsdtt::LandslidePlacer p2(fi, meta, cfg);
        p2.begin_event();
        std::uint64_t rng = 555;
        lsdtt::LandslideScar scar;
        const bool ok = p2.place_scar_seeded(40, rng, scar);
        check(ok, "reference scar placed");

        auto& td = t_fail.get_data_array();
        for (std::int64_t fl : scar.pixels)
            td(static_cast<int>(fl / n_cols), static_cast<int>(fl % n_cols)) = 1000.0f;

        const double anneal = 500.0;

        // Still annealing at t = 1200 (elapsed 200 < 500).
        {
            lsdtt::LandslidePlacer p3(fi, meta, cfg);
            p3.apply_annealing(t_fail, 1200.0, anneal);
            p3.begin_event();
            bool any_blocked_available = false;
            for (std::int64_t fl : scar.pixels)
                if (p3.is_available(fl)) any_blocked_available = true;
            check(!any_blocked_available,
                  "every scar pixel is unavailable 200 yr into a 500 yr anneal");
        }

        // Released at t = 1500 + eps (elapsed 500 >= 500).
        {
            lsdtt::LandslidePlacer p4(fi, meta, cfg);
            p4.apply_annealing(t_fail, 1500.1, anneal);
            p4.begin_event();
            bool all_available = true;
            for (std::int64_t fl : scar.pixels)
                if (!p4.is_available(fl)) all_available = false;
            check(all_available, "every scar pixel is available once anneal_time has elapsed");
        }

        // A blocked pixel must not appear in any new scar.
        {
            lsdtt::LandslidePlacer p5(fi, meta, cfg);
            p5.apply_annealing(t_fail, 1200.0, anneal);
            p5.begin_event();
            std::set<std::int64_t> blocked(scar.pixels.begin(), scar.pixels.end());
            std::uint64_t r2 = 9090;
            bool clean = true;
            for (int i = 0; i < 200; ++i)
            {
                lsdtt::LandslideScar s2;
                if (!p5.place_scar_seeded(15, r2, s2)) continue;
                for (std::int64_t fl : s2.pixels)
                    if (blocked.count(fl)) clean = false;
                // Connectivity must survive the presence of blocked ground --
                // this is the case that breaks a naive nearest-n selection.
                if (!scar_is_connected(s2, meta.n_cols)) clean = false;
                if (static_cast<int>(s2.pixels.size()) != 15) clean = false;
            }
            check(clean,
                  "T8: new scars avoid annealing ground, stay connected and stay exact");
        }
    }

    // -------------------------------------------------------------------------
    // T10 -- mask
    // -------------------------------------------------------------------------
    std::cout << "\nT10: active mask honoured\n";
    {
        lsdtt::lsdtt_array2d_f mdata(n_rows, n_cols, 0.0f);
        lsdtt::lsdtt_array2d_u8 mmask(n_rows, n_cols, 1);
        for (int r = 0; r < n_rows; ++r)
            for (int c = 0; c < n_cols / 2; ++c)
                mdata(r, c) = 1.0f;            // left half active only
        lsdtt::lsdtt_raster_cpu mask_ras(std::move(mdata), std::move(mmask), meta);

        lsdtt::LandslidePlacer pm(fi, meta, cfg, &mask_ras);
        pm.begin_event();
        std::uint64_t rng = 31337;
        bool all_init_in_mask = true;
        int placed = 0;
        for (int i = 0; i < 200; ++i)
        {
            lsdtt::LandslideScar scar;
            if (!pm.place_scar_seeded(10, rng, scar)) continue;
            ++placed;
            const int c = static_cast<int>(scar.initiation_flat % n_cols);
            if (c >= n_cols / 2) all_init_in_mask = false;
        }
        std::cout << "        placed " << placed << " scars in the masked half\n";
        check(placed > 0, "placement succeeds inside the mask");
        check(all_init_in_mask, "every initiation pixel lies inside the mask");
    }

    // -------------------------------------------------------------------------
    // The CLI bridge: build routing through make_flow_base() and drive the
    // placer from the ABSTRACT interface. This is the pattern a CLI must use,
    // because lsdtt_flowinfo_cpu.hpp may not be included at the CLI layer, and
    // it is worth exercising rather than merely assuming it compiles.
    // -------------------------------------------------------------------------
    std::cout << "\nCLI bridge: placement through lsdtt_flow_base\n";
    {
        std::unique_ptr<lsdtt::lsdtt_flow_base> flow_base =
            lsdtt::make_flow_base(filled.get_data_array(),
                                  filled.get_mask_array(), meta, opts);
        check(flow_base != nullptr, "make_flow_base returns a routing object");

        lsdtt::LandslidePlacer pb(*flow_base, meta, cfg);
        check(pb.max_hostable_pixels() == placer.max_hostable_pixels(),
              "abstract-interface placer sees the same hostable maximum");

        pb.begin_event();
        std::uint64_t rng = 20260803;
        lsdtt::LandslideScar scar;
        const bool ok = pb.place_scar_seeded(25, rng, scar);
        check(ok, "a scar places through the abstract interface");
        check(static_cast<int>(scar.pixels.size()) == 25,
              "and reaches its exact target area");
        check(scar_is_connected(scar, meta.n_cols),
              "and is connected");
    }

    // -------------------------------------------------------------------------
    // T13 -- a placement mask must contain the WHOLE scar, not just its seed
    //
    // The mask is a promise about where landslides may occur. Testing it only
    // at the initiation pixel breaks that promise, because a scar occupies its
    // host's upslope set and can spill across the boundary. Measured on a real
    // step-change tag mask before this was fixed: 40 of 1941 excavated pixels
    // lay outside the mask, 35 of them on ground upstream of the knickpoint --
    // exactly the population the mask existed to exclude.
    // -------------------------------------------------------------------------
    {
        std::cout << "\nT13: a mask contains the whole scar, not just its seed\n";

        // Left half active. The boundary runs down the middle of a landscape
        // whose flow crosses it, so an unmasked-upslope bug leaks immediately.
        lsdtt::lsdtt_raster_cpu mask_raster = filled;
        auto& md = mask_raster.get_data_array();
        for (int r = 0; r < meta.n_rows; ++r)
            for (int c = 0; c < meta.n_cols; ++c)
                md(r, c) = (c < meta.n_cols / 2) ? 1.0f : 0.0f;

        auto inside_mask = [&](std::int64_t f) -> bool
        {
            const int rr = static_cast<int>(f / meta.n_cols);
            const int cc = static_cast<int>(f % meta.n_cols);
            return mask_raster.get_mask_array()(rr, cc) &&
                   mask_raster.get_data_array()(rr, cc) != 0.0f;
        };

        lsdtt::LandslidePlacer pm(fi, meta, cfg, &mask_raster);
        pm.begin_event();

        std::uint64_t rng = 987654321;
        int  n_placed = 0;
        long n_outside = 0;
        for (int i = 0; i < 200; ++i)
        {
            lsdtt::LandslideScar scar;
            if (!pm.place_scar_seeded(20, rng, scar)) continue;
            ++n_placed;
            for (std::int64_t f : scar.pixels)
                if (!inside_mask(f)) ++n_outside;
        }

        std::cout << "  placed " << n_placed << " scars, "
                  << n_outside << " pixel(s) outside the mask\n";
        check(n_placed > 0, "scars still place inside a half-domain mask");
        check(n_outside == 0, "no scar pixel falls outside the mask");
    }

    // -------------------------------------------------------------------------
    std::cout << "\n=== " << (g_checks - g_failures) << " / " << g_checks
              << " checks passed ===\n";
    if (g_failures > 0)
        std::cout << "*** " << g_failures << " FAILURES ***\n";
    return g_failures == 0 ? 0 : 1;
}
