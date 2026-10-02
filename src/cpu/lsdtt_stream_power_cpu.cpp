// lsdtt_stream_power_cpu.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html

#include "lsdtt_stream_power_cpu.hpp"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <vector>
#include <iostream>
#include <stdexcept>

#ifdef _OPENMP
#  include <omp.h>
#endif

namespace lsdtt
{

namespace
{

std::vector<double> smooth_drainage_area(
    const std::vector<double>& area,
    const std::vector<std::uint8_t>& mask,
    int n_rows,
    int n_cols,
    int kernel_size,
    const std::string& boundary_conditions)
{
    if (kernel_size < 1 || kernel_size % 2 == 0)
        throw std::invalid_argument(
            "drainage area moving-average kernel size must be a positive odd number");
    if (kernel_size == 1) return area;

    const int radius = kernel_size / 2;
    const bool periodic_north = boundary_conditions.size() == 4 &&
                                (boundary_conditions[0] == 'p' || boundary_conditions[0] == 'P');
    const bool periodic_east = boundary_conditions.size() == 4 &&
                               (boundary_conditions[1] == 'p' || boundary_conditions[1] == 'P');
    const bool periodic_south = boundary_conditions.size() == 4 &&
                                (boundary_conditions[2] == 'p' || boundary_conditions[2] == 'P');
    const bool periodic_west = boundary_conditions.size() == 4 &&
                               (boundary_conditions[3] == 'p' || boundary_conditions[3] == 'P');
    auto wrapped_index = [](int index, int size) {
        return (index % size + size) % size;
    };

    auto resolve_index = [&](int index, int size, bool periodic_negative,
                             bool periodic_positive) {
        if (index < 0 && !periodic_negative) return -1;
        if (index >= size && !periodic_positive) return -1;
        return wrapped_index(index, size);
    };

    const std::size_t n_cells = static_cast<std::size_t>(n_rows) * n_cols;
    std::vector<double> horizontal_sum(n_cells, 0.0);
    std::vector<std::size_t> horizontal_count(n_cells, 0);

    // First pass: rolling horizontal sums and valid-cell counts.
    for (int row = 0; row < n_rows; ++row)
    {
        double sum = 0.0;
        std::size_t count = 0;
        auto add_column = [&](int raw_col, double sign) {
            const int column = resolve_index(
                raw_col, n_cols, periodic_west, periodic_east);
            if (column < 0) return;
            const std::size_t neighbour =
                static_cast<std::size_t>(row) * n_cols + column;
            if (neighbour < area.size() && neighbour < mask.size() &&
                mask[neighbour] != 0)
            {
                sum += sign * area[neighbour];
                if (sign > 0.0) ++count;
                else --count;
            }
        };

        for (int offset = -radius; offset <= radius; ++offset)
            add_column(offset, 1.0);
        for (int col = 0; col < n_cols; ++col)
        {
            const std::size_t flat = static_cast<std::size_t>(row) * n_cols + col;
            horizontal_sum[flat] = sum;
            horizontal_count[flat] = count;
            if (col + 1 < n_cols)
            {
                add_column(col - radius, -1.0);
                add_column(col + radius + 1, 1.0);
            }
        }
    }

    std::vector<double> smoothed(n_cells, 0.0);

    // Second pass: rolling vertical sums over the horizontal-pass results.
    for (int col = 0; col < n_cols; ++col)
    {
        double sum = 0.0;
        std::size_t count = 0;
        auto add_row = [&](int raw_row, double sign) {
            const int row = resolve_index(
                raw_row, n_rows, periodic_north, periodic_south);
            if (row < 0) return;
            const std::size_t neighbour =
                static_cast<std::size_t>(row) * n_cols + col;
            sum += sign * horizontal_sum[neighbour];
            if (sign > 0.0) count += horizontal_count[neighbour];
            else count -= horizontal_count[neighbour];
        };

        for (int offset = -radius; offset <= radius; ++offset)
            add_row(offset, 1.0);
        for (int row = 0; row < n_rows; ++row)
        {
            const std::size_t flat = static_cast<std::size_t>(row) * n_cols + col;
            if (mask[flat] != 0 && count != 0)
                smoothed[flat] = sum / static_cast<double>(count);
            if (row + 1 < n_rows)
            {
                add_row(row - radius, -1.0);
                add_row(row + radius + 1, 1.0);
            }
        }
    }
    return smoothed;
}

} // namespace


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
    float slope_clamp,
    int   drainage_area_moving_average_kernel_size_pixels,
    const std::string& boundary_conditions)
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
    area_d = smooth_drainage_area(
        area_d, dem.get_mask(), meta.n_rows, meta.n_cols,
        drainage_area_moving_average_kernel_size_pixels, boundary_conditions);

    // Topological order is headwater->outlet; we need outlet->headwater,
    // so iterate in reverse.
    const std::vector<std::int64_t>& topo = flow.topological_order();

    // NR parameters
    const int    NR_MAX_ITER = 100;
    const bool   n_is_one    = (std::abs(n - 1.0f) < 1e-6f);

    // Diagonal direction codes per the header convention:
    //   7  0  1
    //   6  x  2
    //   5  4  3
    // Odd codes (1, 3, 5, 7) are diagonal.
    auto is_diagonal_dir = [](std::uint8_t dir) -> bool {
        return (dir & 1u) == 1u;
    };

    // Iterate outlet -> headwater (reverse of topological order)
    for (std::int64_t ti = static_cast<std::int64_t>(topo.size()) - 1;
         ti >= 0; --ti)
    {
        const std::int64_t flat = topo[static_cast<std::size_t>(ti)];

        // Skip outlets - they are baselevel, not eroded
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

        // Node already at or below its receiver - no incision possible
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
            // Newton-Raphson for n != 1:
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
    int   n_threads,
    int   drainage_area_moving_average_kernel_size_pixels,
    const std::string& boundary_conditions)
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
    area_d = smooth_drainage_area(
        area_d, dem.get_mask(), meta.n_rows, meta.n_cols,
        drainage_area_moving_average_kernel_size_pixels, boundary_conditions);

    // Build levels for the reverse sweep (outlet -> headwater)
    auto levels = build_reverse_levels(flow);

    // NR parameters
    const int   NR_MAX_ITER = 100;
    const bool  n_is_one    = (std::abs(n - 1.0f) < 1e-6f);

    // ==============================================================
    // 1. Initialize global counters to track metrics across the timestep
    // ==============================================================
    std::size_t global_slope_clamp_count = 0;
    std::size_t global_overexcavation_count = 0;
    std::size_t global_total_evaluated_nodes = 0;

    for (auto& level : levels)
    {
        const std::int64_t level_size = static_cast<std::int64_t>(level.size());

        // 2. Initialize thread-safe counters for the current topological level
        std::size_t level_slope_clamp_count = 0;
        std::size_t level_overexcavation_count = 0;
        std::size_t level_evaluated_nodes = 0;

#ifdef _OPENMP
        #pragma omp parallel for schedule(static) \
            reduction(+:level_slope_clamp_count, level_overexcavation_count, level_evaluated_nodes)
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

            // 3. Mark the node as actively evaluated for erosion
            level_evaluated_nodes++;

            // Flags to ensure we don't double-count bouncing NR loops
            bool node_clamped = false;
            bool node_overexcavated = false;

            const float z_i_old = z_old[static_cast<std::size_t>(flat)];
            float       z_new   = z_i_old;

            const std::int64_t nr = flow.n_receivers(flat);

            // ----------------------------------------------------------------
            // Classify receiver configuration, mirroring fluvial_incision_Dinf:
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

            const float z_e1 = recv_elev(0);   
            const float z_e2 = recv_elev(1);   

            const double kam = static_cast<double>(K) *
                               std::pow(area, static_cast<double>(m));
            const double elev_tolerance = static_cast<double>(slope_clamp) * dx_cardinal;

            // ----------------------------------------------------------------
            // Two-receiver case
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
                        node_clamped = true; // Track NR loop slope stall
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

                // Instability: elevation increased - revert
                if ((z_new - z_i_old) > static_cast<float>(elev_tolerance) ||
                    std::isnan(z_new))
                {
                    z_new = z_i_old;
                }

                // Overexcavation: clamp above lowest receiver
                const float z_min_recv = std::min(z_e1, z_e2);
                if (z_new < z_min_recv) {
                    z_new = z_min_recv + slope_clamp * dx_cardinal;
                    node_overexcavated = true; // Track overexcavation
                }
            }

            // ----------------------------------------------------------------
            // Single-receiver cases
            // ----------------------------------------------------------------
            else
            {
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
                    z_new = z_i_old;
                }
                else if (n_is_one)
                {
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
                    // Newton-Raphson for n != 1
                    double epsilon = elev_tolerance + 1.0;
                    int    iter_count = 0;

                    do {
                        const double slope = (z_new - z_r) / dx_eff;

                        if (slope <= 0.0) {
                            epsilon = 0.0;
                            node_clamped = true; // Track NR loop slope stall
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

                    // Instability: elevation increased - revert
                    if ((z_new - z_i_old) > static_cast<float>(elev_tolerance) ||
                        std::isnan(z_new))
                    {
                        z_new = z_i_old;
                    }

                    // Overexcavation: clamp above receiver
                    if (z_new < z_r) {
                        z_new = static_cast<float>(z_r + slope_clamp * dx_eff);
                        node_overexcavated = true; // Track overexcavation
                    }
                }
            }
            
            // 4. Update reductions
            if (node_clamped) level_slope_clamp_count++;
            if (node_overexcavated) level_overexcavation_count++;

            z[static_cast<std::size_t>(flat)] = z_new;
        }

        // 5. Commit thread-safe level values up to the globals
        global_slope_clamp_count += level_slope_clamp_count;
        global_overexcavation_count += level_overexcavation_count;
        global_total_evaluated_nodes += level_evaluated_nodes;
    }

    // 6. Print to terminal
    std::cout << global_slope_clamp_count << "," 
              << global_overexcavation_count << "," 
              << global_total_evaluated_nodes << "\n";

    return lsdtt_raster_cpu(
        std::move(z),
        std::vector<uint8_t>(dem.get_mask()),
        meta
    );
}


} // namespace lsdtt

// =============================================================================
// The channel threshold. See the header for the units trap this centralises.
// =============================================================================
namespace lsdtt
{

double channel_threshold_quantity(double upslope_area_m2,
                                  double cell_size_m,
                                  bool   divide_area_by_dx) noexcept
{
    if (!divide_area_by_dx) return upslope_area_m2;
    return upslope_area_m2 / std::max(1e-12, cell_size_m);
}

channel_threshold_report describe_channel_threshold(double threshold,
                                                    double cell_size_m,
                                                    bool   divide_area_by_dx) noexcept
{
    const double dx = std::max(1e-12, cell_size_m);
    channel_threshold_report r;

    if (divide_area_by_dx)
    {
        r.specific_area_m = threshold;
        r.map_area_m2     = threshold * dx;
    }
    else
    {
        r.map_area_m2     = threshold;
        r.specific_area_m = threshold / dx;
    }
    r.equivalent_pixels = r.map_area_m2 / (dx * dx);
    return r;
}

std::vector<std::uint8_t> channel_mask_from_area_threshold(
    const lsdtt_flowinfo&      flow,
    const std::vector<double>& upslope_area_m2,
    double                     threshold,
    double                     cell_size_m,
    bool                       divide_area_by_dx)
{
    const std::size_t n = flow.full_grid_size();
    std::vector<std::uint8_t> is_channel(n, 0u);
    if (!(threshold > 0.0)) return is_channel;

    for (std::size_t i = 0; i < n && i < upslope_area_m2.size(); ++i)
    {
        if (!flow.is_valid_cell(static_cast<std::int64_t>(i))) continue;
        const double q = channel_threshold_quantity(
            upslope_area_m2[i], cell_size_m, divide_area_by_dx);
        if (q >= threshold) is_channel[i] = 1u;
    }
    return is_channel;
}

} // namespace lsdtt
