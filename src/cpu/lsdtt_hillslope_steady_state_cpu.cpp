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
#include <utility>
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

    // 'o' is a DIRICHLET edge: the boundary is base level, so edge cells are held
    // at their input elevation and contribute to their neighbours' updates exactly
    // like channel cells.  'n' is the Neumann (zero-flux) edge.
    //
    // These solvers previously treated 'o' and 'n' alike, both as Neumann
    // reflection, which disagreed with lsdtt_hillslope_implicit_cpu.cpp (which has
    // always honoured 'o' as Dirichlet) and silently converted a fixed-elevation
    // boundary into a no-flux one.  A run declaring 'opop' -- periodic E-W, fixed
    // elevation N-S -- was getting no-flux N-S instead.
    auto bc_code = [&](std::size_t i) -> char
    {
        return (bc.size() > i) ? bc[i] : 'n';
    };
    const char cN = bc_code(0), cE = bc_code(1), cS = bc_code(2), cW = bc_code(3);

    auto is_o_edge = [&](int r, int c) -> bool
    {
        return (r == 0         && cN == 'o') || (r == nrows - 1 && cS == 'o')
            || (c == 0         && cW == 'o') || (c == ncols - 1 && cE == 'o');
    };

    // The SOR arithmetic is carried in double (see the note on z below), so the
    // constants it consumes are formed in double too.
    const double hd   = static_cast<double>(h);
    const double h2d  = static_cast<double>(h2);
    const double Sc2d = static_cast<double>(opts.Sc) * static_cast<double>(opts.Sc);
    const double rhs  = static_cast<double>(opts.U) * h2d / static_cast<double>(opts.D);

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
    //
    // DOUBLE, not float, and deliberately so.  The convergence test is
    // max|z_new - z_old|, an increment that shrinks toward zero while z itself
    // stays O(100 m).  In single precision the representable step near 100 m is
    // ~1e-5 m, and the update forms a sum of four neighbours (O(400)) before
    // differencing it, so the increment cannot be resolved below ~2e-4 m.  The
    // default tolerance of 1e-4 m therefore sat BELOW the arithmetic noise
    // floor: the residual bottomed out around 2e-4 and wandered there forever,
    // so every solve ran to max_iterations no matter how well converged the
    // field was.  In double the floor drops by ~9 orders of magnitude and the
    // tolerance means what it says.  The raster is converted back to float on
    // the way out; only the iteration is carried at higher precision.
    const std::vector<float>& z_in = dem.get_data();
    std::vector<double> z(z_in.begin(), z_in.end());

    // -------------------------------------------------------------------------
    // nbr_z: elevation of the neighbour at (r+dr, c+dc) read from z.
    //
    // Handles periodic wrap and Neumann (zero-flux) reflection at open
    // boundaries and masked interior neighbours.
    // -------------------------------------------------------------------------
    auto nbr_rc = [&](int r, int c, int dr, int dc) -> std::pair<int, int>
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

        // Masked interior neighbour: same Neumann reflection
        if (!mask[static_cast<size_t>(nr * ncols + nc)])
        {
            nr = std::max(0, std::min(nrows - 1, r - dr));
            nc = std::max(0, std::min(ncols - 1, c - dc));
        }

        return {nr, nc};
    };

    auto nbr_z = [&](int r, int c, int dr, int dc) -> double
    {
        const auto rc = nbr_rc(r, c, dr, dc);
        return z[static_cast<size_t>(rc.first * ncols + rc.second)];
    };

    // Centred tangential slopes AT AN ARBITRARY CELL.  Needed because the face
    // coefficient must be a property of the face, not of the cell that happens
    // to be asking (see the note on face symmetry below).
    auto slope_x = [&](int r, int c) -> double
    {
        return (nbr_z(r, c, 0, +1) - nbr_z(r, c, 0, -1)) / (2.0 * hd);
    };
    auto slope_y = [&](int r, int c) -> double
    {
        return (nbr_z(r, c, +1, 0) - nbr_z(r, c, -1, 0)) / (2.0 * hd);
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
    // cross2_over_Sc2 = (tangent_slope)^2 / Sc^2, where the tangent is the
    // E-W slope for N/S faces and the N-S slope for E/W faces.
    //
    // FACE SYMMETRY.  The tangent is averaged ACROSS the face, from the two
    // cells that share it, rather than taken from the cell being updated.  This
    // matters: alpha is a property of the face, so the N face of (r,c) and the
    // S face of (r-1,c) are the same face and must get the same coefficient.
    // Evaluating the tangent at the updating cell alone gave them DIFFERENT
    // values, which makes the discrete operator non-symmetric and the scheme
    // non-conservative -- flux leaving one cell through a face did not equal
    // the flux entering its neighbour through that same face.  The face-normal
    // term is already symmetric, being squared.
    // -------------------------------------------------------------------------
    auto face_alpha = [&](double z_c, double z_nbr, double cross2_over_Sc2) -> double
    {
        if (!opts.nonlinear) return 1.0;
        const double dz_h = (z_nbr - z_c) / hd;
        const double s2   = std::min(dz_h * dz_h / Sc2d + cross2_over_Sc2, 0.9999);
        return 1.0 / (1.0 - s2);
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
    const double omega_d = static_cast<double>(omega);

    double final_residual = static_cast<double>(opts.tolerance) + 1.0;
    int    n_iter         = 0;
    bool   converged      = false;

    for (int it = 0; it < opts.max_iterations; ++it)
    {
        double max_delta = 0.0;

        for (int color = 0; color <= 1; ++color)
        {
            double pass_delta = 0.0;

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
                if (is_o_edge(row, col))            continue;   // Dirichlet base level

                const double zc = z[sidx];

                const double zN = nbr_z(row, col, -1,  0);
                const double zS = nbr_z(row, col, +1,  0);
                const double zE = nbr_z(row, col,  0, +1);
                const double zW = nbr_z(row, col,  0, -1);

                // Face tangents, averaged across each face so that the two
                // cells sharing a face agree on its coefficient.  Only needed
                // in nonlinear mode; the linear path ignores them entirely.
                double cN = 0.0, cS = 0.0, cE = 0.0, cW = 0.0;
                if (opts.nonlinear)
                {
                    const double sxc = slope_x(row, col);
                    const double syc = slope_y(row, col);

                    const auto rcN = nbr_rc(row, col, -1,  0);
                    const auto rcS = nbr_rc(row, col, +1,  0);
                    const auto rcE = nbr_rc(row, col,  0, +1);
                    const auto rcW = nbr_rc(row, col,  0, -1);

                    const double tN = 0.5 * (sxc + slope_x(rcN.first, rcN.second));
                    const double tS = 0.5 * (sxc + slope_x(rcS.first, rcS.second));
                    const double tE = 0.5 * (syc + slope_y(rcE.first, rcE.second));
                    const double tW = 0.5 * (syc + slope_y(rcW.first, rcW.second));

                    cN = tN * tN / Sc2d;
                    cS = tS * tS / Sc2d;
                    cE = tE * tE / Sc2d;
                    cW = tW * tW / Sc2d;
                }

                const double aN = face_alpha(zc, zN, cN);
                const double aS = face_alpha(zc, zS, cS);
                const double aE = face_alpha(zc, zE, cE);
                const double aW = face_alpha(zc, zW, cW);
                const double aT = aN + aS + aE + aW;

                const double z_jac = (aN * zN + aS * zS + aE * zE + aW * zW + rhs) / aT;
                const double z_new = (1.0 - omega_d) * zc + omega_d * z_jac;

                const double delta = std::abs(z_new - zc);
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

    std::vector<float> z_out(z.begin(), z.end());

    lsdtt_hillslope_steady_state_result result;
    result.elevation      = lsdtt_raster_cpu(
        std::move(z_out),
        std::vector<uint8_t>(mask),
        meta);
    result.n_iterations   = n_iter;
    result.final_residual = static_cast<float>(final_residual);
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
