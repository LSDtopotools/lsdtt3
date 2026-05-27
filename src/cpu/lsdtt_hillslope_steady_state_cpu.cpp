// lsdtt_hillslope_steady_state_cpu.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html

#include "lsdtt_hillslope_steady_state_cpu.hpp"

#include <algorithm>
#include <cmath>
#include <vector>

#ifdef _OPENMP
#  include <omp.h>
#endif

namespace lsdtt
{

lsdtt_hillslope_steady_state_result compute_hillslope_steady_state(
    const lsdtt_raster_cpu&                     dem,
    const std::vector<uint8_t>&                 is_channel,
    const lsdtt_hillslope_steady_state_options& opts)
{
#ifdef _OPENMP
    if (opts.n_threads > 0) omp_set_num_threads(opts.n_threads);
#endif

    const lsdtt_raster_meta& meta   = dem.get_meta();
    const int     nrows   = meta.n_rows;
    const int     ncols   = meta.n_cols;
    const float   h       = static_cast<float>(meta.cell_size);
    const float   h2      = h * h;
    const int64_t n_cells = static_cast<int64_t>(nrows) * ncols;

    // Boundary condition flags: boundary_conditions string is NESW order.
    const std::string& bc = opts.boundary_conditions;
    const bool per_n = (bc.size() > 0 && bc[0] == 'p');
    const bool per_e = (bc.size() > 1 && bc[1] == 'p');
    const bool per_s = (bc.size() > 2 && bc[2] == 'p');
    const bool per_w = (bc.size() > 3 && bc[3] == 'p');

    // RHS: U*h²/D — the uplift contribution that hillslope sediment flux must balance.
    const float rhs = opts.U * h2 / opts.D;
    const float Sc2 = opts.Sc * opts.Sc;

    const std::vector<uint8_t>& mask = dem.get_mask();

    // Jacobi requires two buffers: z_cur is read, z_nxt is written, then swapped.
    std::vector<float> z_cur(dem.get_data());
    std::vector<float> z_nxt(z_cur);

    // -------------------------------------------------------------------------
    // nbr_z: elevation of the neighbour displaced by (dr, dc) from (r, c).
    //
    // Handles:
    //   - periodic wrap ('p' boundary)
    //   - Neumann reflection for 'o' boundaries: ghost cell mirrors inward cell
    //     using centered-difference zero-flux: z_ghost = z_inward_neighbour
    //   - masked interior neighbours: same Neumann reflection inward
    // -------------------------------------------------------------------------
    auto nbr_z = [&](const std::vector<float>& z, int r, int c,
                     int dr, int dc) -> float
    {
        int nr = r + dr;
        int nc = c + dc;

        // North/South boundary
        if (nr < 0)
        {
            if (per_n) nr = nrows - 1;
            else       nr = std::max(0, r - dr);     // Neumann: reflect inward
        }
        else if (nr >= nrows)
        {
            if (per_s) nr = 0;
            else       nr = std::min(nrows - 1, r - dr);
        }

        // East/West boundary
        if (nc < 0)
        {
            if (per_w) nc = ncols - 1;
            else       nc = std::max(0, c - dc);
        }
        else if (nc >= ncols)
        {
            if (per_e) nc = 0;
            else       nc = std::min(ncols - 1, c - dc);
        }

        const size_t nidx = static_cast<size_t>(nr * ncols + nc);

        // Masked interior neighbour: same Neumann reflection
        if (!mask[nidx])
        {
            const int ir = std::max(0, std::min(nrows - 1, r - dr));
            const int ic = std::max(0, std::min(ncols - 1, c - dc));
            return z[static_cast<size_t>(ir * ncols + ic)];
        }

        return z[nidx];
    };

    // -------------------------------------------------------------------------
    // alpha: face diffusivity coefficient.
    //   Linear:    alpha = 1.
    //   Nonlinear: alpha = 1 / (1 - (S/Sc)²),  S = (z_nbr - z_c) / h.
    //   Slope is clamped to (1-ε)*Sc to keep alpha finite during iteration.
    // -------------------------------------------------------------------------
    auto face_alpha = [&](float z_c, float z_nbr) -> float
    {
        if (!opts.nonlinear) return 1.0f;
        const float dz_h = (z_nbr - z_c) / h;
        const float s2   = std::min(dz_h * dz_h / Sc2, 0.9999f);
        return 1.0f / (1.0f - s2);
    };

    // -------------------------------------------------------------------------
    // Jacobi iteration.
    //
    // For each non-channel cell the update is:
    //
    //   z_new = (α_N*z_N + α_S*z_S + α_E*z_E + α_W*z_W + U*h²/D)
    //           / (α_N + α_S + α_E + α_W)
    //
    // α_k = 1 for linear; α_k = 1/(1-(S_k/Sc)²) for Roering.
    // All α are evaluated from z_cur (Picard/lagged).
    // Channel and nodata cells are copied unchanged.
    // Convergence: max|Δz| < tolerance.
    // -------------------------------------------------------------------------
    float final_residual = opts.tolerance + 1.0f;
    int   n_iter         = 0;
    bool  converged      = false;

    for (int it = 0; it < opts.max_iterations; ++it)
    {
        float max_delta = 0.0f;

#ifdef _OPENMP
        #pragma omp parallel for schedule(static) reduction(max : max_delta)
#endif
        for (int64_t idx = 0; idx < n_cells; ++idx)
        {
            const size_t sidx = static_cast<size_t>(idx);

            if (!mask[sidx] || is_channel[sidx])
            {
                z_nxt[sidx] = z_cur[sidx];
                continue;
            }

            const int   row = static_cast<int>(idx / ncols);
            const int   col = static_cast<int>(idx % ncols);
            const float zc  = z_cur[sidx];

            const float zN = nbr_z(z_cur, row, col, -1,  0);
            const float zS = nbr_z(z_cur, row, col, +1,  0);
            const float zE = nbr_z(z_cur, row, col,  0, +1);
            const float zW = nbr_z(z_cur, row, col,  0, -1);

            const float aN = face_alpha(zc, zN);
            const float aS = face_alpha(zc, zS);
            const float aE = face_alpha(zc, zE);
            const float aW = face_alpha(zc, zW);
            const float aT = aN + aS + aE + aW;

            const float z_new = (aN * zN + aS * zS + aE * zE + aW * zW + rhs) / aT;
            const float delta = std::abs(z_new - zc);
            if (delta > max_delta) max_delta = delta;

            z_nxt[sidx] = z_new;
        }

        std::swap(z_cur, z_nxt);
        ++n_iter;
        final_residual = max_delta;

        if (max_delta <= opts.tolerance)
        {
            converged = true;
            break;
        }
    }

    lsdtt_hillslope_steady_state_result result;
    result.elevation      = lsdtt_raster_cpu(
        std::move(z_cur),
        std::vector<uint8_t>(mask),
        meta);
    result.n_iterations   = n_iter;
    result.final_residual = final_residual;
    result.converged      = converged;
    return result;
}

} // namespace lsdtt
