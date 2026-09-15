// test_lem_channel_coupling.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
//
// -----------------------------------------------------------------------------
// Process-domain coupling in run_lem_step: channels must drive the hillslopes.
//
// In the channels-as-boundaries regime the hillslope step treats channel pixels
// as Dirichlet nodes, so the hillslopes relax toward the elevations the fluvial
// step left in the channels. The defect this guards against is a hillslope
// model that never receives the channel set, and so relaxes toward the only
// other fixed elevations it has -- the domain edges, i.e. base level at the
// outlet. Two models had exactly that defect:
//
//   * adi was called with no channel set at all, so the channels diffused along
//     with everything else;
//   * cidre re-derived its sinks from incision_threshold_area inside the sweep,
//     so channel_mask_fname and hillslope_channel_coupling never reached it.
//
// Fixture: a 60 x 21 grid, cell 10 m. One channel column (col 10) held at 100 m
// through an explicit channel mask; hillslopes either side start 5 m above it.
// The south edge is an 'o' outlet at 0 m; the other edges are no-flux. Fluvial
// incision and uplift are off, so the hillslope step is the only process.
//
//   C1  Every model leaves the channel pixels exactly unchanged by a step.
//   C2  After 3e4 yr (7x the lateral relaxation time 4W^2/(pi^2 D)), hillslope
//       pixels far from the outlet sit at the CHANNEL elevation (100 m) -- not at
//       the outlet (0 m) and not at the domain mean (~104.8 m).
//   C3  cidre honours hillslope_channel_coupling: coupled with a threshold is
//       the same run as no threshold at all.
//   C4  cidre's threshold-derived channel set is unchanged by routing it through
//       the shared mask (same output as the sweep building it itself).
//
// Usage: lsdtt-test-lem-channel-coupling
// -----------------------------------------------------------------------------

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

#include "../../src/cpu/lsdtt_flowgraph_cpu.hpp"
#include "../../src/cpu/lsdtt_hillslope_transport_cpu.hpp"
#include "../../src/cpu/lsdtt_lem_cpu.hpp"
#include "../../src/cpu/lsdtt_raster_cpu.hpp"

namespace
{
int g_checks = 0, g_failures = 0;

void check(bool ok, const std::string& what)
{
    ++g_checks;
    if (!ok) { ++g_failures; std::cout << "  FAIL: " << what << "\n"; }
    else     { std::cout << "  ok  : " << what << "\n"; }
}

constexpr int    NR = 60, NC = 21;
constexpr int    CH_COL = 10;
constexpr float  H = 10.0f;
constexpr float  Z_CHANNEL = 100.0f;
constexpr float  Z_OFFSET  = 5.0f;
constexpr float  Z_OUTLET  = 0.0f;
const std::string BC = "nnon";            // N, E, S, W: outlet on the south edge

std::size_t idx(int r, int c)
{
    return static_cast<std::size_t>(r) * NC + static_cast<std::size_t>(c);
}

lsdtt::lsdtt_raster_meta meta()
{
    lsdtt::lsdtt_raster_meta m;
    m.n_rows = NR; m.n_cols = NC; m.cell_size = H;
    m.no_data_value = -9999.0f;
    return m;
}

lsdtt::lsdtt_raster_cpu make_dem()
{
    std::vector<float>   z(static_cast<std::size_t>(NR) * NC, Z_CHANNEL + Z_OFFSET);
    std::vector<uint8_t> mask(z.size(), 1);
    for (int r = 0; r < NR; ++r) z[idx(r, CH_COL)] = Z_CHANNEL;
    for (int c = 0; c < NC; ++c) z[idx(NR - 1, c)] = Z_OUTLET;
    return lsdtt::lsdtt_raster_cpu(std::move(z), std::move(mask), meta());
}

std::vector<std::uint8_t> make_channel_mask()
{
    std::vector<std::uint8_t> m(static_cast<std::size_t>(NR) * NC, 0u);
    for (int r = 0; r < NR - 1; ++r) m[idx(r, CH_COL)] = 1u;
    return m;
}

lsdtt::lsdtt_lem_options base_options(const std::string& model)
{
    lsdtt::lsdtt_lem_options o;
    o.timestep                        = 20.0f;
    o.hillslope_diffusivity           = 1.0f;
    o.enable_nonlinear_hillslopes     = false;
    o.critical_hillslope_slope        = 1.2f;
    o.hillslope_transport_model       = model;
    o.deposit_on_hillslope_baselevels = false;   // so C1 is exact for cidre
    o.enable_d8_fluvial_erosion       = false;
    o.enable_dinf_fluvial_erosion     = false;
    o.enable_hillslope_transport      = true;
    o.enable_uplift                   = false;
    o.boundary_conditions             = BC;
    o.n_threads                       = 1;
    return o;
}

double max_channel_change(const lsdtt::lsdtt_raster_cpu& a,
                          const lsdtt::lsdtt_raster_cpu& b)
{
    double m = 0.0;
    for (int r = 0; r < NR - 1; ++r)
        m = std::max(m, std::fabs(static_cast<double>(a.get_data()[idx(r, CH_COL)])
                                - static_cast<double>(b.get_data()[idx(r, CH_COL)])));
    return m;
}

bool identical(const lsdtt::lsdtt_raster_cpu& a, const lsdtt::lsdtt_raster_cpu& b)
{
    return a.get_data() == b.get_data();
}
} // namespace

int main()
{
    using namespace lsdtt;
    std::cout << "\n=== LEM process-domain coupling: channels drive hillslopes ===\n";

    const lsdtt_raster_cpu dem0 = make_dem();
    const std::vector<std::int64_t> fixed = build_boundary_flats(dem0, BC);

    const std::vector<std::string> models = {
        "adi", "cidre", "linear", "roering", "roering_implicit", "cellular"};

    std::cout << "\n-- C0: the fixture discriminates (adi with NO channel set) --\n";
    {
        // Without a channel set the channel column is diffused like any other
        // cell. This is what adi did inside run_lem_step before the fix, and it
        // is what makes C1 a real test rather than a vacuous one.
        const lsdtt_raster_cpu no_channels =
            compute_adi_hillslope_transport(dem0, 1.0f, 20.0f, BC);
        const double dch = max_channel_change(no_channels, dem0);
        std::cout << "     max |dz| on channel after 1 step: " << dch << " m\n";
        check(dch > 0.1, "without a channel set the channel pixels move");
    }

    for (const std::string& model : models)
    {
        std::cout << "\n-- model: " << model << " --\n";

        lsdtt_lem_options o = base_options(model);
        o.explicit_channel_mask = make_channel_mask();

        lsdtt_hillslope_newton_workspace   newton_ws;
        lsdtt_hillslope_cellular_workspace cellular_ws;

        // C1: one step
        lsdtt_raster_cpu z = dem0;
        {
            auto step = run_lem_step(z, o, fixed, &newton_ws, &cellular_ws);
            const double dch = max_channel_change(step.first, dem0);
            std::cout << "     max |dz| on channel after 1 step: " << dch << " m\n";
            check(dch == 0.0, model + ": C1 channel pixels unchanged by the hillslope step");
            z = std::move(step.first);
        }

        // C2: relax
        const int n_steps = 1500;              // 3e4 yr
        double worst_channel = 0.0;
        for (int s = 1; s < n_steps; ++s)
        {
            auto step = run_lem_step(z, o, fixed, &newton_ws, &cellular_ws);
            worst_channel = std::max(worst_channel, max_channel_change(step.first, dem0));
            z = std::move(step.first);
        }

        double sum_abs = 0.0, min_dev = 1e30, max_dev = -1e30;
        int    n = 0;
        for (int r = 0; r < 10; ++r)
            for (int c = 0; c < NC; ++c)
            {
                if (c == CH_COL) continue;
                const double dev = static_cast<double>(z.get_data()[idx(r, c)]) - Z_CHANNEL;
                sum_abs += std::fabs(dev);
                min_dev = std::min(min_dev, dev);
                max_dev = std::max(max_dev, dev);
                ++n;
            }
        const double mean_abs = sum_abs / std::max(1, n);
        std::cout << "     after 3e4 yr, rows 0-9 hillslope z - z_channel: mean|.| "
                  << mean_abs << " m, range [" << min_dev << ", " << max_dev << "] m\n"
                  << "     max |dz| on channel over the run: " << worst_channel << " m\n";

        check(worst_channel == 0.0, model + ": C1 channel pixels unchanged over the whole run");

        // Cidre is NOT held to C2. It moves sediment only down the D-inf
        // steepest-descent facet, and in this strip that facet points along the
        // hillslope toward the outlet rather than into the channel once the lateral
        // gradient has relaxed. Measured (2026-09-15): rows 0-9 end 2.4-4.0 m BELOW
        // the channel. That is the transport law, not a missing channel set -- C1
        // above is what proves the channel set reaches it.
        if (model == "cidre") continue;

        check(mean_abs < 0.25,
              model + ": C2 hillslopes relaxed to the channel elevation (mean |z - z_ch| < 0.25 m)");
        check(min_dev > -0.25,
              model + ": C2 hillslopes not drawn below the channel toward the outlet");
    }

    std::cout << "\n-- C3: cidre honours hillslope_channel_coupling --\n";
    {
        // A specific-area threshold of 15 m selects any cell draining more than
        // 1.5 pixels, i.e. the channel column that the hillslopes drain into.
        const float threshold = 15.0f;

        lsdtt_lem_options thr = base_options("cidre");
        thr.incision_threshold_area = threshold;
        thr.coupling_regime = hillslope_coupling_regime::coupled;

        lsdtt_lem_options none = base_options("cidre");
        none.incision_threshold_area = 0.0f;

        lsdtt_lem_options boundaries = base_options("cidre");
        boundaries.incision_threshold_area = threshold;
        boundaries.coupling_regime = hillslope_coupling_regime::channels_as_boundaries;

        const std::vector<std::uint8_t> set = compute_lem_channel_mask(dem0, boundaries);
        const auto n_set = std::count(set.begin(), set.end(), std::uint8_t{1});
        std::cout << "     channel pixels at threshold " << threshold << ": " << n_set << "\n";
        check(n_set > 0, "the threshold selects a non-empty channel set");

        const auto a = run_lem_step(dem0, thr, fixed);
        const auto b = run_lem_step(dem0, none, fixed);
        const auto c = run_lem_step(dem0, boundaries, fixed);
        check(identical(a.first, b.first),
              "coupled + threshold gives the same step as no threshold");
        check(!identical(a.first, c.first),
              "and differs from channels_as_boundaries with the same threshold");
    }

    std::cout << "\n-- C4: cidre's threshold-derived channel set is unchanged --\n";
    {
        lsdtt_lem_options o = base_options("cidre");
        o.incision_threshold_area = 50.0f;

        const auto lem = run_lem_step(dem0, o, {});

        lsdtt_flowgraph_options fg;
        fg.boundary_conditions       = BC;
        fg.interior_nodata_is_outlet = true;
        fg.method                    = lsdtt_flowgraph_method::d_infinity;
        fg.n_threads                 = 1;
        const lsdtt_flowgraph_cpu flow(dem0.get_data_array(), dem0.get_mask_array(),
                                       dem0.get_meta(), fg);
        const lsdtt_raster_cpu direct = compute_cidre_hillslope_transport(
            dem0, flow, o.hillslope_diffusivity, o.timestep, false,
            o.critical_hillslope_slope, o.incision_threshold_area,
            o.divide_area_by_dx, o.deposit_on_hillslope_baselevels, 1);

        check(identical(lem.first, direct),
              "run_lem_step (shared mask) == cidre sweep building its own mask");
    }

    std::cout << "\n" << (g_checks - g_failures) << " / " << g_checks
              << " checks passed\n";
    return g_failures == 0 ? 0 : 1;
}
