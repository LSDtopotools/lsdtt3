// lsdtt_lem_cpu.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html

#include "lsdtt_lem_cpu.hpp"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <vector>

#include "lsdtt_flowgraph_cpu.hpp"
#include "lsdtt_hillslope_implicit_cpu.hpp"
#include "lsdtt_hillslope_transport_cpu.hpp"
#include "lsdtt_stream_power_cpu.hpp"
#include <iostream>

namespace lsdtt
{

std::pair<lsdtt_raster_cpu, lsdtt_lem_step_stats>
run_lem_step(
    const lsdtt_raster_cpu&          dem,
    const lsdtt_lem_options&         opts,
    const std::vector<std::int64_t>& boundary_flats)
{
    const lsdtt_raster_meta& meta = dem.get_meta();
    const std::size_t n_cells = static_cast<std::size_t>(meta.n_rows) * meta.n_cols;

    // Save pre-step elevations of boundary pixels
    const std::vector<float>& z_pre = dem.get_data();
    std::vector<float> boundary_elev(boundary_flats.size());
    for (std::size_t i = 0; i < boundary_flats.size(); ++i)
    {
        const auto flat = static_cast<std::size_t>(boundary_flats[i]);
        if (flat < n_cells)
            boundary_elev[i] = z_pre[flat];
    }

    // Snapshot of pre-step elevations for statistics (valid cells only)
    const std::vector<float> z_before = z_pre;

    lsdtt_raster_cpu current = dem;

    // -------------------------------------------------------------------------
    // 1. Fluvial incision
    // -------------------------------------------------------------------------
    if (opts.enable_d8_fluvial_erosion)
    {
        // Build D8 flow graph (only when needed)
        lsdtt_flowinfo_options fi_opts;
        fi_opts.boundary_conditions       = opts.boundary_conditions;
        fi_opts.interior_nodata_is_outlet = true;
        fi_opts.allow_pits_as_outlets     = true; // LEM surface may have transient pits

        lsdtt_flowinfo flow_d8(
            dem.get_data_array(),
            dem.get_mask_array(),
            dem.get_meta(),
            fi_opts
        );

        current = compute_d8_stream_power_incision(
            current, flow_d8,
            opts.erodibility,
            opts.m_exponent,
            opts.n_exponent,
            opts.timestep,
            opts.incision_threshold_area,
            opts.divide_area_by_dx,
            opts.slope_clamp
        );
    }
    else if (opts.enable_dinf_fluvial_erosion)
    {
        // Build D-infinity flow graph for fluvial incision
        lsdtt_flowgraph_options fg_opts;
        fg_opts.boundary_conditions       = opts.boundary_conditions;
        fg_opts.interior_nodata_is_outlet = true;
        fg_opts.method                    = lsdtt_flowgraph_method::d_infinity;
        fg_opts.n_threads                 = opts.n_threads;

        lsdtt_flowgraph_cpu flow_dinf(
            dem.get_data_array(),
            dem.get_mask_array(),
            dem.get_meta(),
            fg_opts
        );

        current = compute_dinf_stream_power_incision(
            current, flow_dinf,
            opts.erodibility,
            opts.m_exponent,
            opts.n_exponent,
            opts.timestep,
            opts.incision_threshold_area,
            opts.divide_area_by_dx,
            opts.slope_clamp,
            opts.n_threads
        );
    }

    // -------------------------------------------------------------------------
    // 2. Hillslope transport
    //    Built from the post-incision topography.
    //
    //    Cidre (default): builds a D-inf flow graph and uses the topological
    //      forward sweep for both routing and channel identification.
    //    Roering: builds a D-inf flow graph only for channel identification,
    //      then applies an explicit FD stencil (no topological ordering needed).
    // -------------------------------------------------------------------------
    if (opts.enable_hillslope_transport)
    {

        if (opts.hillslope_transport_model == "adi") 
        {
            // Alternating directions implicit does not need to build a flow graph
            current = compute_adi_hillslope_transport(
                current,
                opts.hillslope_diffusivity,
                opts.timestep,
                opts.boundary_conditions
            );
        }
        else 
        {
            // The linear/roering and cidre hillslope options require a flow
            // graph to be made.
            lsdtt_flowgraph_options fg_opts;
            fg_opts.boundary_conditions       = opts.boundary_conditions;
            fg_opts.interior_nodata_is_outlet = true;
            fg_opts.method                    = lsdtt_flowgraph_method::d_infinity;
            fg_opts.n_threads                 = opts.n_threads;

            lsdtt_flowgraph_cpu flow_hillslope(
                current.get_data_array(),
                current.get_mask_array(),
                current.get_meta(),
                fg_opts
            );

            // Build the Dirichlet channel mask once for the channel-using explicit
            // FD models: prefer an explicit channel mask when supplied (§8), else
            // derive it from the drainage-area threshold.
            const std::int64_t n_hs =
                static_cast<std::int64_t>(meta.n_rows) * meta.n_cols;
            auto build_is_channel = [&]() -> std::vector<uint8_t>
            {
                const float dx = static_cast<float>(meta.cell_size);
                std::vector<uint8_t> is_channel(static_cast<std::size_t>(n_hs), 0);
                if (opts.explicit_channel_mask.size()
                        == static_cast<std::size_t>(n_hs))
                {
                    for (std::int64_t flat = 0; flat < n_hs; ++flat)
                        is_channel[static_cast<std::size_t>(flat)] =
                            opts.explicit_channel_mask[static_cast<std::size_t>(flat)] ? 1 : 0;
                }
                else if (opts.incision_threshold_area > 0.0f)
                {
                    std::vector<double> area = flow_hillslope.accumulate();
                    for (std::int64_t flat = 0; flat < n_hs; ++flat)
                    {
                        if (!flow_hillslope.is_valid_flat(flat)) continue;
                        double a = area[static_cast<std::size_t>(flat)];
                        if (opts.divide_area_by_dx) a /= static_cast<double>(dx);
                        if (a >= static_cast<double>(opts.incision_threshold_area))
                            is_channel[static_cast<std::size_t>(flat)] = 1;
                    }
                }
                return is_channel;
            };

            if (opts.hillslope_transport_model == "linear" ||
                opts.hillslope_transport_model == "roering")
            {
                // Explicit FD on the regular grid — flow graph used only for channel ID.
                std::vector<uint8_t> is_channel = build_is_channel();

                const bool nl = (opts.hillslope_transport_model == "roering")
                                && opts.enable_nonlinear_hillslopes;

                current = compute_roering_hillslope_transport(
                    current,
                    is_channel,
                    opts.hillslope_diffusivity,
                    opts.timestep,
                    nl,
                    opts.critical_hillslope_slope,
                    opts.n_threads
                );
            }
            else if (opts.hillslope_transport_model == "roering_implicit")
            {
                // Unconditionally-stable backward-Euler Roering solver (§8).
                std::vector<uint8_t> is_channel = build_is_channel();

                // Transient-coupling rule 2 (§8): with no Dirichlet channel sinks the
                // solve degenerates to sink-free diffusion toward the edges, which is
                // meaningless — skip the hillslope step entirely.
                const bool any_channel =
                    std::any_of(is_channel.begin(), is_channel.end(),
                                [](uint8_t v) { return v != 0; });
                if (!any_channel)
                {
                    std::cout << "[lem] roering_implicit: no channel pixels "
                                 "(empty is_channel) — skipping the hillslope step "
                                 "(needs incision_threshold_area > 0 or a non-empty "
                                 "channel_mask_fname).\n";
                }
                else
                {
                    lsdtt_roering_implicit_options ri = opts.roering_implicit_opts;
                    ri.D  = opts.hillslope_diffusivity;
                    ri.U  = 0.0f;                       // uplift applied in step 3
                    ri.Sc = opts.critical_hillslope_slope;
                    ri.nonlinear = opts.enable_nonlinear_hillslopes;
                    ri.boundary_conditions = opts.boundary_conditions;
                    ri.n_threads = opts.n_threads;

                    // Transient-coupling rule 1 (§8): integrate exactly one LEM step
                    // (dt_total = opts.timestep), so channels are re-incised each step.
                    lsdtt_roering_implicit_result ri_res =
                        compute_roering_hillslope_transport_implicit(
                            current, is_channel, opts.timestep, ri);
                    current = std::move(ri_res.elevation);
                }
            }
            else  // "cidre"
            {
                current = compute_cidre_hillslope_transport(
                    current, flow_hillslope,
                    opts.hillslope_diffusivity,
                    opts.timestep,
                    opts.enable_nonlinear_hillslopes,
                    opts.critical_hillslope_slope,
                    opts.incision_threshold_area,
                    opts.divide_area_by_dx,
                    opts.deposit_on_hillslope_baselevels,
                    opts.n_threads
                );
            }
        }
    }

    // -------------------------------------------------------------------------
    // 3. Uplift (uniform)
    // -------------------------------------------------------------------------
    const std::vector<uint8_t>& mask = current.get_mask();
    std::vector<float>& z = current.get_data();

    if (opts.enable_uplift)
    {
        const float uplift_dz = opts.uplift_rate * opts.timestep;
        for (std::size_t i = 0; i < n_cells; ++i)
        {
            if (mask[i])
                z[i] += uplift_dz;
        }
    }

    // -------------------------------------------------------------------------
    // 4. Statistics (computed before boundary restoration so they reflect the
    //    full-domain change)
    // -------------------------------------------------------------------------
    lsdtt_lem_step_stats stats{};
    {
        double sum_before = 0.0, sum_after = 0.0;
        double max_before = -1e30, max_after = -1e30;
        double max_local  = 0.0;
        double sum_erosion = 0.0;
        double max_erosion = 0.0;
        std::int64_t n_valid = 0;

        for (std::size_t i = 0; i < n_cells; ++i)
        {
            if (!mask[i])
                continue;
            ++n_valid;
            const double zb = static_cast<double>(z_before[i]);
            const double za = static_cast<double>(z[i]);

            sum_before += zb;
            sum_after  += za;
            max_before  = std::max(max_before, zb);
            max_after   = std::max(max_after,  za);

            const double local_change = std::abs(za - zb);
            max_local = std::max(max_local, local_change);

            // Erosion = how much lower the surface is now (before uplift, but
            // we already applied uplift, so we back it out for the rate)
            const double uplift_dz = opts.enable_uplift
                ? static_cast<double>(opts.uplift_rate * opts.timestep)
                : 0.0;
            const double erosion = zb - (za - uplift_dz);  // positive = erosion
            if (erosion > 0.0)
            {
                sum_erosion += erosion;
                max_erosion  = std::max(max_erosion, erosion);
            }
        }

        const double inv_dt = 1.0 / static_cast<double>(opts.timestep);
        const double inv_n  = n_valid > 0 ? 1.0 / static_cast<double>(n_valid) : 0.0;

        stats.delta_max_elevation        = static_cast<float>(std::abs(max_after - max_before));
        stats.delta_mean_elevation       = static_cast<float>(
            std::abs((sum_after - sum_before) * inv_n));
        stats.max_local_elevation_change = static_cast<float>(max_local);
        stats.mean_erosion_rate          = static_cast<float>(sum_erosion * inv_n * inv_dt);
        stats.max_erosion_rate           = static_cast<float>(max_erosion * inv_dt);
    }

    // -------------------------------------------------------------------------
    // 5. Restore fixed-elevation boundary pixels
    // -------------------------------------------------------------------------
    for (std::size_t i = 0; i < boundary_flats.size(); ++i)
    {
        const auto flat = static_cast<std::size_t>(boundary_flats[i]);
        if (flat < n_cells)
            z[flat] = boundary_elev[i];
    }

    return {std::move(current), stats};
}

} // namespace lsdtt
