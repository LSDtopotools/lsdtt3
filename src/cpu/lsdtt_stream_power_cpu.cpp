// lsdtt_stream_power_cpu.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html

#include "lsdtt_stream_power_cpu.hpp"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <vector>
#include <iostream>

#ifdef _OPENMP
#  include <omp.h>
#endif

namespace lsdtt
{


// -----------------------------------------------------------------------------
// compute_d8_stream_power_incision
// -----------------------------------------------------------------------------
lsdtt_raster_cpu compute_d8_stream_power_incision(
    const lsdtt_raster_cpu&  dem,
    const lsdtt_flowinfo& flow,
    float K,
    float m,
    float n,
    float dt,
    float incision_threshold_area,
    bool  divide_area_by_dx,
    float slope_clamp)
{
    const lsdtt_raster_meta& meta       = dem.get_meta();
    const float              dx         = static_cast<float>(meta.cell_size);
    const float              dx_diag    = dx * std::sqrt(2.0f);
    const std::int64_t       n_cells    = static_cast<std::int64_t>(meta.n_rows)
                                          * meta.n_cols;

    // Working elevation array updated in-place
    std::vector<float>       z(dem.get_data());
    const std::vector<float>& z_old = dem.get_data();

    // Drainage area via D8 accumulation
    std::vector<double> area_d = flow.accumulate();

    // Topological order is headwater→outlet; we need outlet→headwater,
    // so iterate in reverse.
    const std::vector<std::int64_t>& topo = flow.topological_order();

    // NR parameters
    const int    NR_MAX_ITER = 100;
    const double NR_TOL      = 1e-6;
    const bool   n_is_one    = (std::abs(n - 1.0f) < 1e-6f);

    // Diagonal direction codes per the header convention:
    //   7  0  1
    //   6  x  2
    //   5  4  3
    // Odd codes (1, 3, 5, 7) are diagonal.
    auto is_diagonal_dir = [](std::uint8_t dir) -> bool {
        return (dir & 1u) == 1u;
    };

    // Iterate outlet → headwater (reverse of topological order)
    for (std::int64_t ti = static_cast<std::int64_t>(topo.size()) - 1;
         ti >= 0; --ti)
    {
        const std::int64_t flat = topo[static_cast<std::size_t>(ti)];

        // Skip outlets — they are baselevel, not eroded
        if (flow.is_outlet_unchecked(flat))
            continue;

        // Skip nodata
        if (!flow.is_valid_cell(flat))
            continue;

        // Receiver index and grid spacing
        const std::int64_t  recv = flow.receiver_of_unchecked(flat);
        const std::uint8_t  dir  = flow.receiver_dir_of_unchecked(flat);

        if (recv < 0 || recv >= n_cells)
            continue;

        const float dx_eff = is_diagonal_dir(dir) ? dx_diag : dx;

        // Drainage area
        double area = area_d[static_cast<std::size_t>(flat)];
        if (divide_area_by_dx)
            area /= dx;

        if (incision_threshold_area > 0.0f &&
            area < static_cast<double>(incision_threshold_area))
            continue;

        const float  z_i_old = z_old[static_cast<std::size_t>(flat)];
        const float  z_r     = z[static_cast<std::size_t>(recv)];   // already updated
        float        z_new   = z_i_old;

        // Node already at or below its receiver — no incision possible
        if (z_r >= z_i_old)
            continue;

        const double elev_tolerance = static_cast<double>(slope_clamp) * dx_eff;
        const double kam = static_cast<double>(K) *
                           std::pow(area, static_cast<double>(m));

        if (n_is_one)
        {
            // Closed-form implicit update (Braun & Willett 2013, eq. 2):
            //   z_new = (z_old + K*A^m*(dt/dx)*z_r) / (1 + K*A^m*(dt/dx))
            const double spf = kam * (dt / dx_eff);
            z_new = static_cast<float>((z_i_old + z_r * spf) / (1.0 + spf));

            // Guard: elevation must not increase
            if ((z_new - z_i_old) > static_cast<float>(elev_tolerance) ||
                std::isnan(z_new))
            {
                z_new = z_i_old;
            }
        }
        else
        {
            // Newton-Raphson for n ≠ 1:
            //   f(z)  = z - z_old + K*A^m * dt * ((z - z_r)/dx)^n = 0
            //   f'(z) = 1 + K*A^m * dt * (n/dx) * ((z - z_r)/dx)^(n-1)
            double epsilon = elev_tolerance + 1.0;
            int    iter_count = 0;

            do {
                const double slope = (z_new - z_r) / dx_eff;

                if (slope <= 0.0) {
                    // Slope gone negative: stall; overexcavation clamp below
                    epsilon = 0.0;
                } else {
                    const double f_z  = z_new - z_i_old
                                      + kam * dt * std::pow(slope, n);
                    const double f_dz = 1.0
                                      + kam * dt * (n / dx_eff)
                                        * std::pow(slope, n - 1.0);
                    epsilon = f_z / f_dz;
                }
                z_new -= static_cast<float>(epsilon);

                ++iter_count;
                if (iter_count > NR_MAX_ITER)
                    epsilon = elev_tolerance / 2.0;

            } while (std::abs(epsilon) > elev_tolerance);

            // Guard: elevation must not increase
            if ((z_new - z_i_old) > static_cast<float>(elev_tolerance) ||
                std::isnan(z_new))
            {
                z_new = z_i_old;
            }

            // Overexcavation: clamp to minimum slope above receiver
            if (z_new < z_r)
                z_new = static_cast<float>(z_r + slope_clamp * dx_eff);
        }

        z[static_cast<std::size_t>(flat)] = z_new;
    }

    return lsdtt_raster_cpu(
        std::move(z),
        std::vector<uint8_t>(dem.get_mask()),
        meta
    );
}
lsdtt_raster_cpu compute_dinf_stream_power_incision(
    const lsdtt_raster_cpu&    dem,
    const lsdtt_flowgraph_cpu& flow,
    float K,
    float m,
    float n,
    float dt,
    float incision_threshold_area,
    bool  divide_area_by_dx,
    float slope_clamp,
    int   n_threads)
{
#ifdef _OPENMP
    if (n_threads > 0)
        omp_set_num_threads(n_threads);
#else
    (void)n_threads;
#endif

    const lsdtt_raster_meta& meta = dem.get_meta();
    const float dx_cardinal  = static_cast<float>(meta.cell_size);
    const float dx_diagonal  = static_cast<float>(meta.cell_size) * std::sqrt(2.0f);
    const std::int64_t n_cells = static_cast<std::int64_t>(meta.n_rows) * meta.n_cols;

    // Working elevation array, updated in-place during the sweep
    std::vector<float> z(dem.get_data());
    const std::vector<float>& z_old = dem.get_data();

    // Drainage area from the flow graph (forward accumulation)
    std::vector<double> area_d = flow.accumulate();

    // Build levels for the reverse sweep (outlet → headwater)
    auto levels = build_reverse_levels(flow);

    // NR parameters
    const int   NR_MAX_ITER = 100;
    const float NR_TOL      = 1e-6f;
    const bool  n_is_one    = (std::abs(n - 1.0f) < 1e-6f);

    for (auto& level : levels)
    {
        const std::int64_t level_size = static_cast<std::int64_t>(level.size());

#ifdef _OPENMP
        #pragma omp parallel for schedule(static)
#endif
        for (std::int64_t li = 0; li < level_size; ++li)
        {
            const std::int64_t flat = level[static_cast<std::size_t>(li)];

            if (flow.is_outlet_flat(flat)) continue;
            if (!flow.is_valid_flat(flat)) continue;

            double area = area_d[static_cast<std::size_t>(flat)];
            if (divide_area_by_dx)
                area /= dx_cardinal;

            if (incision_threshold_area > 0.0f &&
                area < static_cast<double>(incision_threshold_area))
                continue;

            const float z_i_old = z_old[static_cast<std::size_t>(flat)];
            float       z_new   = z_i_old;

            const std::int64_t nr = flow.n_receivers(flat);

            // ----------------------------------------------------------------
            // Classify receiver configuration, mirroring fluvial_incision_Dinf:
            //   receiver 0 = cardinal (weight==1 → only cardinal)
            //   receiver 1 = diagonal (weight==1 → only diagonal)
            //   both weights non-zero → two-receiver case
            // ----------------------------------------------------------------
            float w0 = (nr > 0) ? flow.receiver_weight_at(flat, 0) : 0.0f;
            float w1 = (nr > 1) ? flow.receiver_weight_at(flat, 1) : 0.0f;

            // Snap near-unity weights to exactly 1.0 to match reference logic
            bool only_cardinal = (nr == 1) || (w0 >= 1.0f - 1e-6f);
            bool only_diagonal = (nr > 1)  && (w1 >= 1.0f - 1e-6f);
            bool two_receivers = (nr > 1)  && !only_cardinal && !only_diagonal;

            // Retrieve receiver elevations (already updated, earlier in sweep)
            auto recv_elev = [&](int k) -> float {
                if (k >= nr) return z_i_old;
                const std::int64_t idx = flow.receiver_index_at(flat, k);
                if (idx < 0 || idx >= n_cells) return z_i_old;
                return z[static_cast<std::size_t>(idx)];
            };

            const float z_e1 = recv_elev(0);   // cardinal receiver elevation
            const float z_e2 = recv_elev(1);   // diagonal receiver elevation (if present)

            const double kam = static_cast<double>(K) *
                               std::pow(area, static_cast<double>(m));
            const double elev_tolerance = static_cast<double>(slope_clamp) * dx_cardinal;

            // ----------------------------------------------------------------
            // Two-receiver case  (Anand et al. 2020, eq. 5)
            // f(z)  = z - z_old + K*A^m * [(z_e1-z)^2/dx^2 + (z_e2-z_e1)^2/dx^2]^(n/2) * dt
            // f'(z) = 1 - K*A^m * n * (z_e1-z)/dx^2 * sbt^(n/2 - 1) * dt
            // Both distances use dx_cardinal per the reference implementation.
            // ----------------------------------------------------------------
            if (two_receivers)
            {
                double epsilon = elev_tolerance + 1.0;
                int    iter_count = 0;

                do {
                    double sbt = std::pow((z_e1 - z_new) / dx_cardinal, 2.0)
                               + std::pow((z_e2 - z_e1) / dx_cardinal, 2.0);

                    if (sbt <= 0.0) {
                        epsilon = 0.0;
                    } else {
                        const double f_z  = z_new - z_i_old
                                          + kam * std::pow(sbt, 0.5 * n) * dt;
                        const double f_dz = 1.0
                                          - kam * n * (z_e1 - z_new)
                                            / (dx_cardinal * dx_cardinal)
                                            * std::pow(sbt, 0.5 * n - 1.0) * dt;
                        epsilon = f_z / f_dz;
                    }
                    z_new -= static_cast<float>(epsilon);

                    ++iter_count;
                    if (iter_count > NR_MAX_ITER)
                        epsilon = elev_tolerance / 2.0;

                } while (std::abs(epsilon) > elev_tolerance);

                // Instability: elevation increased — revert
                if ((z_new - z_i_old) > static_cast<float>(elev_tolerance) ||
                    std::isnan(z_new))
                {
                    z_new = z_i_old;
                }

                // Overexcavation: clamp above lowest receiver
                const float z_min_recv = std::min(z_e1, z_e2);
                if (z_new < z_min_recv)
                    z_new = z_min_recv + slope_clamp * dx_cardinal;
            }

            // ----------------------------------------------------------------
            // Single-receiver cases
            // ----------------------------------------------------------------
            else
            {
                // Pick the correct receiver elevation and grid spacing
                float z_r;
                float dx_eff;
                if (only_diagonal) {
                    z_r    = z_e2;
                    dx_eff = dx_diagonal;
                } else {
                    z_r    = z_e1;
                    dx_eff = dx_cardinal;
                }

                if (z_r >= z_i_old)
                {
                    // Node already at or below receiver — no incision
                    z_new = z_i_old;
                }
                else if (n_is_one)
                {
                    // Closed-form implicit update (Braun & Willett 2013, eq. 2)
                    const double spf = kam * (dt / dx_eff);
                    z_new = static_cast<float>(
                        (z_i_old + z_r * spf) / (1.0 + spf));

                    if ((z_new - z_i_old) > static_cast<float>(elev_tolerance) ||
                        std::isnan(z_new))
                    {
                        z_new = z_i_old;
                    }
                }
                else
                {
                    // Newton-Raphson for n ≠ 1
                    double epsilon = elev_tolerance + 1.0;
                    int    iter_count = 0;

                    do {
                        const double slope = (z_new - z_r) / dx_eff;

                        if (slope <= 0.0) {
                            epsilon = 0.0;
                        } else {
                            const double f_z  = z_new - z_i_old
                                              + kam * dt * std::pow(slope, n);
                            const double f_dz = 1.0
                                              + kam * dt * (n / dx_eff)
                                                * std::pow(slope, n - 1.0);
                            epsilon = f_z / f_dz;
                        }
                        z_new -= static_cast<float>(epsilon);

                        ++iter_count;
                        if (iter_count > NR_MAX_ITER)
                            epsilon = elev_tolerance / 2.0;

                    } while (std::abs(epsilon) > elev_tolerance);

                    // Instability: elevation increased — revert
                    if ((z_new - z_i_old) > static_cast<float>(elev_tolerance) ||
                        std::isnan(z_new))
                    {
                        z_new = z_i_old;
                    }

                    // Overexcavation: clamp above receiver
                    if (z_new < z_r)
                        z_new = static_cast<float>(z_r + slope_clamp * dx_eff);
                }
            }

            z[static_cast<std::size_t>(flat)] = z_new;
        }
    }

    return lsdtt_raster_cpu(
        std::move(z),
        std::vector<uint8_t>(dem.get_mask()),
        meta
    );
}

} // namespace lsdtt
