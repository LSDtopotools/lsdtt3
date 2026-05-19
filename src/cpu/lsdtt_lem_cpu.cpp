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
#include "lsdtt_hillslope_transport_cpu.hpp"
#include "lsdtt_stream_power_cpu.hpp"

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

    // -------------------------------------------------------------------------
    // 1. Build D-infinity flow graph
    // -------------------------------------------------------------------------
    lsdtt_flowgraph_options fg_opts;
    fg_opts.boundary_conditions       = opts.boundary_conditions;
    fg_opts.interior_nodata_is_outlet = true;
    fg_opts.method                    = lsdtt_flowgraph_method::d_infinity;
    fg_opts.n_threads                 = opts.n_threads;

    lsdtt_flowgraph_cpu flow(
        dem.get_data_array(),
        dem.get_mask_array(),
        dem.get_meta(),
        fg_opts
    );

    lsdtt_raster_cpu current = dem;

    // -------------------------------------------------------------------------
    // 2. Fluvial incision
    // -------------------------------------------------------------------------
    if (opts.enable_fluvial_erosion)
    {
        current = compute_stream_power_incision(
            current, flow,
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
    // 3. Hillslope transport (uses the same flow graph built from the DEM at
    //    the start of the step — consistent with a first-order operator split)
    // -------------------------------------------------------------------------
    if (opts.enable_hillslope_transport)
    {
        current = compute_cidre_hillslope_transport(
            current, flow,
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

    // -------------------------------------------------------------------------
    // 4. Uplift (uniform)
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
    // 5. Statistics (computed before boundary restoration so they reflect the
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
    // 6. Restore fixed-elevation boundary pixels
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
