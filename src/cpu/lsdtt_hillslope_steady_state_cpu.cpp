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

#include "lsdtt_flowgraph_cpu.hpp"

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

    const float rhs = opts.U * h2 / opts.D;
    const float Sc2 = opts.Sc * opts.Sc;

    // -------------------------------------------------------------------------
    // Determine SOR relaxation factor omega.
    //
    // For the 2D Poisson equation on an NxM grid with Dirichlet BCs at all
    // four boundaries the optimal omega is derived from the Jacobi spectral radius:
    //
    //   rho_J  = [cos(pi/(N+1)) + cos(pi/(M+1))] / 2
    //   omega_opt = 2 / (1 + sqrt(1 - rho_J^2))
    //
    // When channels are densely distributed through the domain (typical LEM
    // setup) the effective grid spacing is smaller than the full domain, so
    // the above formula overestimates omega.  We therefore cap it at 1.95.
    // The user may override via sor_omega > 0.
    // -------------------------------------------------------------------------
    float omega;
    if (opts.sor_omega > 0.0f)
    {
        omega = opts.sor_omega;
    }
    else
    {
        static constexpr float k_pi_f = 3.14159265f;
        const float rho_j = 0.5f * (std::cos(k_pi_f / static_cast<float>(nrows + 1))
                                  + std::cos(k_pi_f / static_cast<float>(ncols + 1)));
        omega = 2.0f / (1.0f + std::sqrt(1.0f - rho_j * rho_j));
        omega = std::min(omega, 1.95f);
    }

    const std::vector<uint8_t>& mask = dem.get_mask();

    // Single elevation array - SOR updates in-place.
    std::vector<float> z(dem.get_data());

    // -------------------------------------------------------------------------
    // nbr_z: elevation of the neighbour at (r+dr, c+dc) read from z.
    //
    // Handles periodic wrap and Neumann (zero-flux) reflection at open
    // boundaries and masked interior neighbours.
    // -------------------------------------------------------------------------
    auto nbr_z = [&](int r, int c, int dr, int dc) -> float
    {
        int nr = r + dr;
        int nc = c + dc;

        if (nr < 0)
        {
            if (per_n) nr = nrows - 1;
            else       nr = std::max(0, r - dr);   // Neumann: reflect inward
        }
        else if (nr >= nrows)
        {
            if (per_s) nr = 0;
            else       nr = std::min(nrows - 1, r - dr);
        }

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
    // face_alpha: diffusivity at the face between z_c and z_nbr.
    //
    // Linear: always 1.
    //
    // Nonlinear (Roering 1999): the flux law uses the FULL gradient magnitude:
    //   qs = -D * gradz / (1 - |gradz|^2/Sc^2),   |gradz|^2 = (dz/dx)^2 + (dz/dy)^2
    //
    // Each face contributes a face-normal slope and a face-tangent slope.
    // Using only the face-normal component (the error in the old code)
    // underestimates alpha on any cell with a non-zero cross-slope, allowing
    // slopes to grow steeper than Sc in the cross direction.
    //
    // cross2_over_Sc2 = (tangent_slope)^2 / Sc^2  - pre-computed per face pair:
    //   N/S faces share the E-W tangent  -> cross2 = (zE - zW)^2 / (4h^2 Sc^2)
    //   E/W faces share the N-S tangent  -> cross2 = (zN - zS)^2 / (4h^2 Sc^2)
    // -------------------------------------------------------------------------
    auto face_alpha = [&](float z_c, float z_nbr, float cross2_over_Sc2) -> float
    {
        if (!opts.nonlinear) return 1.0f;
        const float dz_h = (z_nbr - z_c) / h;
        const float s2   = std::min(dz_h * dz_h / Sc2 + cross2_over_Sc2, 0.9999f);
        return 1.0f / (1.0f - s2);
    };

    // -------------------------------------------------------------------------
    // Red-Black SOR.
    //
    // Two colour passes per iteration.  In colour pass `color`, only cells
    // where (row + col) % 2 == color are updated.  Because the 4-connected
    // neighbours of a red cell are always black (and vice versa), updates
    // within a single pass are independent -> safe to parallelise with OMP.
    //
    // Update rule:
    //   z_Jacobi = (alphaN*zN + alphaS*zS + alphaE*zE + alphaW*zW + U*h^2/D) / (alphaN+alphaS+alphaE+alphaW)
    //   z_new    = (1 - omega) * z_old  +  omega * z_Jacobi
    //
    // The face alphas are evaluated from the current z (Picard linearisation
    // for the nonlinear mode), using the freshly-updated values of already-
    // processed neighbours - i.e. Gauss-Seidel style within each pass.
    // -------------------------------------------------------------------------
    float final_residual = opts.tolerance + 1.0f;
    int   n_iter         = 0;
    bool  converged      = false;

    for (int it = 0; it < opts.max_iterations; ++it)
    {
        float max_delta = 0.0f;

        for (int color = 0; color <= 1; ++color)
        {
            float pass_delta = 0.0f;

#ifdef _OPENMP
            #pragma omp parallel for schedule(static) reduction(max : pass_delta)
#endif
            for (int64_t idx = 0; idx < n_cells; ++idx)
            {
                const size_t sidx = static_cast<size_t>(idx);
                const int    row  = static_cast<int>(idx / ncols);
                const int    col  = static_cast<int>(idx % ncols);

                if ((row + col) % 2 != color)       continue;
                if (!mask[sidx] || is_channel[sidx]) continue;

                const float zc = z[sidx];

                const float zN = nbr_z(row, col, -1,  0);
                const float zS = nbr_z(row, col, +1,  0);
                const float zE = nbr_z(row, col,  0, +1);
                const float zW = nbr_z(row, col,  0, -1);

                const float sx2_norm = (zE - zW) * (zE - zW) / (4.0f * h2 * Sc2);
                const float sy2_norm = (zN - zS) * (zN - zS) / (4.0f * h2 * Sc2);

                const float aN = face_alpha(zc, zN, sx2_norm);
                const float aS = face_alpha(zc, zS, sx2_norm);
                const float aE = face_alpha(zc, zE, sy2_norm);
                const float aW = face_alpha(zc, zW, sy2_norm);
                const float aT = aN + aS + aE + aW;

                const float z_jac = (aN * zN + aS * zS + aE * zE + aW * zW + rhs) / aT;
                const float z_new = (1.0f - omega) * zc + omega * z_jac;

                const float delta = std::abs(z_new - zc);
                if (delta > pass_delta) pass_delta = delta;

                z[sidx] = z_new;
            }

            if (pass_delta > max_delta) max_delta = pass_delta;
        }

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
        std::move(z),
        std::vector<uint8_t>(mask),
        meta);
    result.n_iterations   = n_iter;
    result.final_residual = final_residual;
    result.converged      = converged;
    return result;
}


// -----------------------------------------------------------------------------
// compute_cidre_hillslope_steady_state
// -----------------------------------------------------------------------------
lsdtt_hillslope_steady_state_result compute_cidre_hillslope_steady_state(
    const lsdtt_raster_cpu&                     dem,
    const std::vector<uint8_t>&                 is_channel,
    const lsdtt_hillslope_steady_state_options& opts)
{
#ifdef _OPENMP
    if (opts.n_threads > 0) omp_set_num_threads(opts.n_threads);
#endif

    const lsdtt_raster_meta& meta  = dem.get_meta();
    const float   dx      = static_cast<float>(meta.cell_size);
    const int64_t n_cells = static_cast<int64_t>(meta.n_rows) * meta.n_cols;
    const float   Sc2     = opts.Sc * opts.Sc;

    const std::vector<uint8_t>& mask = dem.get_mask();

    // Pseudo-time step: CFL bound for linear Cidre diffusion (kappa*dt/dx <= 0.25*dx).
    // The slope-restoration branch of nonlinear Cidre is unconditionally stable, so
    // this conservative choice works for both modes.
    const float dt_ss      = 0.25f * dx / opts.D;   // -> kappa*dt/dx = 0.25*dx
    const float kdt_dx     = opts.D * dt_ss / dx;   // = 0.25*dx
    const float uplift_dz  = opts.U * dt_ss;

    std::vector<float> z(dem.get_data());

    // Save channel elevations - restored after every iteration.
    std::vector<float> channel_z(static_cast<std::size_t>(n_cells), 0.0f);
    for (int64_t i = 0; i < n_cells; ++i)
        if (is_channel[static_cast<std::size_t>(i)])
            channel_z[static_cast<std::size_t>(i)] = z[static_cast<std::size_t>(i)];

    // D-inf flow graph options (rebuilt each iteration for current topography).
    lsdtt_flowgraph_options fg_opts;
    fg_opts.boundary_conditions = opts.boundary_conditions;
    fg_opts.method              = lsdtt_flowgraph_method::d_infinity;
    fg_opts.n_threads           = opts.n_threads;

    float final_residual = opts.tolerance + 1.0f;
    int   n_iter         = 0;
    bool  converged      = false;

    for (int it = 0; it < opts.max_iterations; ++it)
    {
        const std::vector<float> z_before(z);

        // -----------------------------------------------------------------
        // Rebuild the D-inf flow graph from the current elevation field
        // (Picard linearization: routing directions update each iteration).
        // -----------------------------------------------------------------
        lsdtt_raster_cpu cur(std::vector<float>(z), std::vector<uint8_t>(mask), meta);

        lsdtt_flowgraph_cpu flow(
            cur.get_data_array(), cur.get_mask_array(), meta, fg_opts);

        const std::vector<float>& slopes = flow.dinf_slopes_raw();

        auto levels = build_forward_levels(flow);

        std::vector<float> flux(static_cast<std::size_t>(n_cells), 0.0f);

        // -----------------------------------------------------------------
        // Forward level sweep - mirrors compute_cidre_hillslope_transport
        // but uses the precomputed is_channel mask for Dirichlet BCs.
        // -----------------------------------------------------------------
        for (auto& level : levels)
        {
            const int64_t lsz = static_cast<int64_t>(level.size());

#ifdef _OPENMP
            #pragma omp parallel for schedule(static)
#endif
            for (int64_t li = 0; li < lsz; ++li)
            {
                const int64_t  flat = level[static_cast<std::size_t>(li)];
                const std::size_t sidx = static_cast<std::size_t>(flat);

                if (!flow.is_valid_flat(flat)) continue;

                const bool is_ch   = static_cast<bool>(is_channel[sidx]);
                const bool is_sink = is_ch || flow.is_outlet_flat(flat);
                const float flux_in = flux[sidx];

                if (is_sink)
                {
                    // Channel/outlet acts as a drain: incoming flux deposits but
                    // channel elevations are restored at the end of the iteration.
                    if (!is_ch)
                        z[sidx] += flux_in;  // non-channel outlet (rare on conditioned DEMs)
                    continue;
                }

                float S = slopes[sidx];
                if (S < 0.0f) S = 0.0f;

                float deposit, erosion, outgoing;

                if (!opts.nonlinear)
                {
                    // Linear Cidre: all incoming flux deposits; only erosion travels on.
                    erosion  = kdt_dx * S;
                    deposit  = flux_in;
                    outgoing = erosion;
                }
                else
                {
                    if (S < opts.Sc)
                    {
                        const float S2 = S * S;
                        erosion  = kdt_dx * S;
                        deposit  = flux_in * (1.0f - S2 / Sc2);
                        outgoing = erosion + flux_in * (S2 / Sc2);
                    }
                    else
                    {
                        // Slope-restoration: remove excess above Sc; all material moves on.
                        erosion  = dx * (S - opts.Sc);
                        deposit  = 0.0f;
                        outgoing = erosion + flux_in;
                    }
                }

                z[sidx] += deposit - erosion;

                const int64_t nr = flow.n_receivers(flat);
                for (int64_t k = 0; k < nr; ++k)
                {
                    const int64_t recv = flow.receiver_index_at(flat, k);
                    const float   w    = flow.receiver_weight_at(flat, k);
                    if (recv >= 0 && recv < n_cells)
                    {
#ifdef _OPENMP
                        #pragma omp atomic
#endif
                        flux[static_cast<std::size_t>(recv)] += outgoing * w;
                    }
                }
            }
        }

        // -----------------------------------------------------------------
        // Apply uplift to hillslope pixels; restore channel Dirichlet values.
        // -----------------------------------------------------------------
        for (int64_t i = 0; i < n_cells; ++i)
        {
            const std::size_t sidx = static_cast<std::size_t>(i);
            if (!mask[sidx]) continue;
            if (is_channel[sidx]) { z[sidx] = channel_z[sidx]; continue; }
            z[sidx] += uplift_dz;
        }

        // -----------------------------------------------------------------
        // Convergence check (hillslope pixels only).
        // -----------------------------------------------------------------
        float max_delta = 0.0f;
#ifdef _OPENMP
        #pragma omp parallel for schedule(static) reduction(max : max_delta)
#endif
        for (int64_t i = 0; i < n_cells; ++i)
        {
            const std::size_t sidx = static_cast<std::size_t>(i);
            if (!mask[sidx] || is_channel[sidx]) continue;
            const float d = std::abs(z[sidx] - z_before[sidx]);
            if (d > max_delta) max_delta = d;
        }

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
        std::move(z), std::vector<uint8_t>(mask), meta);
    result.n_iterations   = n_iter;
    result.final_residual = final_residual;
    result.converged      = converged;
    return result;
}

} // namespace lsdtt
