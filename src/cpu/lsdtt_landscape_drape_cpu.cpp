// lsdtt_landscape_drape_cpu.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html

#include "lsdtt_landscape_drape_cpu.hpp"
#include "../utils/lsdtt_grid_neighbours.hpp"

#include "../utils/lsdtt_thread_utils.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <vector>
#include <iostream>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace lsdtt {

lsdtt_cloth_drape_result compute_cloth_drape(
    const lsdtt_raster_cpu&          dem,
    const lsdtt_cloth_drape_options& opts)
{
    const lsdtt_raster_meta& meta  = dem.get_meta();
    const int nrows  = meta.n_rows;
    const int ncols  = meta.n_cols;
    const int n      = nrows * ncols;
    const float nodata = static_cast<float>(meta.no_data_value);

    const std::vector<float>&   dem_data = dem.get_data();
    const std::vector<uint8_t>& dem_mask = dem.get_mask();

    // ------------------------------------------------------------------
    // Find the maximum valid DEM elevation - cloth initialisation height
    // ------------------------------------------------------------------
    float max_elev = -std::numeric_limits<float>::max();
    for (int i = 0; i < n; ++i)
        if (dem_mask[i])
            max_elev = std::max(max_elev, dem_data[i]);

    if (max_elev == -std::numeric_limits<float>::max())
        throw std::runtime_error(
            "compute_cloth_drape: DEM contains no valid pixels");

    // ------------------------------------------------------------------
    // Allocate cloth state: z_curr, z_prev, z_next
    // Valid pixels start at max_elev (flat sheet); nodata pixels stay at nodata.
    // ------------------------------------------------------------------
    std::vector<float>   z_curr(n, nodata);
    std::vector<float>   z_prev(n, nodata);
    std::vector<float>   z_next(n, nodata);
    std::vector<uint8_t> pinned(n, 0u);

    for (int i = 0; i < n; ++i)
    {
        if (dem_mask[i])
        {
            z_curr[i] = max_elev;
            z_prev[i] = max_elev;
        }
    }

    // ------------------------------------------------------------------
    // 8-connected neighbour offsets
    // ------------------------------------------------------------------
    static constexpr const auto& k_drow = grid::d8_rowmajor_row;
    static constexpr const auto& k_dcol = grid::d8_rowmajor_col;

    const float dt2         = opts.time_step * opts.time_step;
    const int   n_threads   = lsdtt::resolve_thread_count(opts.n_threads);

    int   iter     = 0;
    float max_disp = std::numeric_limits<float>::max();

    // ------------------------------------------------------------------
    // Verlet integration loop
    // ------------------------------------------------------------------
    for (iter = 0; iter < opts.max_iterations; ++iter)
    {
        // ---- parallel update ----
        #pragma omp parallel for schedule(static) num_threads(n_threads)
        for (int row = 0; row < nrows; ++row)
        {
            for (int col = 0; col < ncols; ++col)
            {
                const int i = row * ncols + col;

                if (!dem_mask[i])
                {
                    z_next[i] = nodata;
                    continue;
                }

                if (pinned[i])
                {
                    z_next[i] = z_curr[i];
                    continue;
                }

                // Spring force: mean displacement from valid neighbours
                float spring    = 0.0f;
                int   n_valid   = 0;

                for (int d = 0; d < 8; ++d)
                {
                    const int nr = row + k_drow[d];
                    const int nc = col + k_dcol[d];

                    if (nr < 0 || nr >= nrows || nc < 0 || nc >= ncols) continue;

                    const int j = nr * ncols + nc;
                    if (!dem_mask[j]) continue;

                    spring += z_curr[j] - z_curr[i];
                    ++n_valid;
                }

                if (n_valid > 0)
                    spring = opts.rigidity * spring / static_cast<float>(n_valid);

                // Verlet update
                float z_new = 2.0f * z_curr[i] - z_prev[i]
                              + (spring - opts.gravity) * dt2;

                // Contact constraint: pin when cloth reaches terrain
                if (z_new < dem_data[i] + opts.class_threshold)
                {
                    z_new      = dem_data[i];
                    pinned[i]  = 1u;
                }

                z_next[i] = z_new;
            }
        }

        // ---- convergence check (sequential, O(n)) ----
        // After swap, z_prev will hold the old z_curr and z_curr will hold z_next.
        // Displacement = |z_next - z_curr| before the swap.
        max_disp = 0.0f;
        for (int i = 0; i < n; ++i)
        {
            if (!dem_mask[i] || pinned[i]) continue;
            const float disp = std::abs(z_next[i] - z_curr[i]);
            if (disp > max_disp) max_disp = disp;
        }

        std::swap(z_prev, z_curr);
        std::swap(z_curr, z_next);

        if (max_disp < opts.convergence_tol)
        {
            ++iter;   // report the iteration that converged
            break;
        }
    }

    // ------------------------------------------------------------------
    // Build output rasters
    // ------------------------------------------------------------------
    std::vector<float>   env_data(n,    nodata);
    std::vector<uint8_t> env_mask(n,    0u);
    std::vector<float>   inc_data(n,    nodata);
    std::vector<uint8_t> inc_mask(n,    0u);
    std::vector<float>   crest_data(n,  nodata);
    std::vector<uint8_t> crest_msk(n,   0u);

    for (int i = 0; i < n; ++i)
    {
        if (!dem_mask[i]) continue;

        env_data[i]   = z_curr[i];
        env_mask[i]   = 1u;

        inc_data[i]   = z_curr[i] - dem_data[i];
        inc_mask[i]   = 1u;

        crest_data[i] = static_cast<float>(pinned[i]);
        crest_msk[i]  = 1u;
    }

    lsdtt_cloth_drape_result result;
    result.envelope   = lsdtt_raster_cpu(std::move(env_data),   std::move(env_mask),   meta);
    result.incision   = lsdtt_raster_cpu(std::move(inc_data),   std::move(inc_mask),   meta);
    result.crest_mask = lsdtt_raster_cpu(std::move(crest_data), std::move(crest_msk),  meta);
    result.iterations_taken       = iter;
    result.final_max_displacement = max_disp;

    std::cout << "compute_cloth_drape: converged in " << iter
              << " iterations (max displacement: " << max_disp << " m)\n";

    return result;
}

} // namespace lsdtt
