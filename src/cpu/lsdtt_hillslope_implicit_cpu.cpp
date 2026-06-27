// lsdtt_hillslope_implicit_cpu.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html
//
// -----------------------------------------------------------------------------
// Implicit transient Roering hillslope solver.
//
// Phase 0 scaffolding: the function is callable and returns a valid, filled
// lsdtt_roering_implicit_result regardless of whether Eigen is present.  The
// Eigen branch currently logs its intent, builds the hillslope DOF map, and then
// delegates to the explicit subcycling solver; the backward-Euler + Picard core
// is staged behind TODO tags matching the phased checklist (§12 of
// docs/roering_implicit_transient_plan.md).  When Eigen is unavailable the
// function falls back to the explicit solver and prints a one-line notice, so
// callers need no #ifdef.
// -----------------------------------------------------------------------------

#include "lsdtt_hillslope_implicit_cpu.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

#include "lsdtt_hillslope_transport_cpu.hpp"

#ifdef LSDTT_HAS_EIGEN
#  include <Eigen/Sparse>
#endif

namespace lsdtt
{

namespace
{

// Compute the maximum |∇z|/Sc over the hillslope (non-channel, valid) cells of a
// field, using the same centred-difference full-gradient estimate as the explicit
// solver.  Returned ratio is *not* capped, so a value > 1 flags a super-critical
// face that the limiter would act on.  Used to populate the result diagnostic
// max_slope_ratio in both branches.
float compute_max_slope_ratio(
    const lsdtt_raster_cpu&     dem,
    const std::vector<uint8_t>& is_channel,
    float                       Sc)
{
    const lsdtt_raster_meta&    meta = dem.get_meta();
    const int                   nrows = meta.n_rows;
    const int                   ncols = meta.n_cols;
    const float                 h     = static_cast<float>(meta.cell_size);
    const std::vector<float>&   z     = dem.get_data();
    const std::vector<uint8_t>& mask  = dem.get_mask();

    if (Sc <= 0.0f || nrows < 1 || ncols < 1)
    {
        return 0.0f;
    }

    auto nbr_z = [&](int r, int c, int dr, int dc) -> float
    {
        const int nr = std::max(0, std::min(nrows - 1, r + dr));
        const int nc = std::max(0, std::min(ncols - 1, c + dc));
        const std::size_t nidx = static_cast<std::size_t>(nr * ncols + nc);
        if (!mask[nidx])
        {
            return z[static_cast<std::size_t>(r * ncols + c)];
        }
        return z[nidx];
    };

    float max_ratio = 0.0f;
    const std::int64_t n_cells = static_cast<std::int64_t>(nrows) * ncols;
    for (std::int64_t idx = 0; idx < n_cells; ++idx)
    {
        const std::size_t sidx = static_cast<std::size_t>(idx);
        if (!mask[sidx] || (sidx < is_channel.size() && is_channel[sidx]))
        {
            continue;
        }
        const int row = static_cast<int>(idx / ncols);
        const int col = static_cast<int>(idx % ncols);
        const float gx = (nbr_z(row, col, 0, +1) - nbr_z(row, col, 0, -1)) / (2.0f * h);
        const float gy = (nbr_z(row, col, -1, 0) - nbr_z(row, col, +1, 0)) / (2.0f * h);
        const float ratio = std::sqrt(gx * gx + gy * gy) / Sc;
        if (ratio > max_ratio)
        {
            max_ratio = ratio;
        }
    }
    return max_ratio;
}

} // namespace

// -----------------------------------------------------------------------------
// compute_roering_hillslope_transport_implicit
// -----------------------------------------------------------------------------
lsdtt_roering_implicit_result compute_roering_hillslope_transport_implicit(
    const lsdtt_raster_cpu&                 dem,
    const std::vector<uint8_t>&             is_channel,
    float                                   dt_total,
    const lsdtt_roering_implicit_options&   opts)
{
    lsdtt_roering_implicit_result result;

#ifdef LSDTT_HAS_EIGEN
    // -------------------------------------------------------------------------
    // Eigen path — backward-Euler + damped Picard core (§4).
    //
    // Phases 1 & 2: conservative cell-centred finite-volume backward Euler over
    // fixed sub-steps, with a damped Picard (lagged-coefficient) iteration for
    // the nonlinear face diffusivity (§3).  Linear mode (nonlinear == false) uses
    // D_eff = D and converges in one Picard pass.  Dirichlet sinks come from
    // is_channel AND the 'o' (fixed-elevation) boundary code; 'n' = no-flux,
    // 'p' = periodic (§3).  Uplift U is folded into the RHS for hillslope DOFs.
    // The system is SPD (symmetric face-centred D_eff) and solved with
    // Eigen::SimplicialLDLT.
    //
    // TODO(Phase 3): adaptive sub-stepping + post-step projection / rejection.
    // TODO(Phase 4): analyzePattern() reuse, OpenMP triplet assembly, CG backend.
    // TODO(Phase 6): optional Newton / TR-BDF2 cores.
    // -------------------------------------------------------------------------
    const lsdtt_raster_meta&    meta    = dem.get_meta();
    const int                   n_rows  = meta.n_rows;
    const int                   n_cols  = meta.n_cols;
    const std::int64_t          n_cells = static_cast<std::int64_t>(n_rows) * n_cols;
    const std::vector<uint8_t>& mask    = dem.get_mask();
    const double                h       = static_cast<double>(meta.cell_size);
    const double                h2      = h * h;
    const double                D       = static_cast<double>(opts.D);
    const double                U       = static_cast<double>(opts.U);
    const double                Sc      = (opts.Sc > 0.0f) ? static_cast<double>(opts.Sc) : 1.0;
    const bool                  nonlinear = opts.nonlinear;
    const double                r_max   = static_cast<double>(opts.slope_ratio_cap);
    const double                omega   = (opts.relaxation_omega > 0.0f)
                                          ? static_cast<double>(opts.relaxation_omega) : 1.0;
    const double                picard_tol = static_cast<double>(opts.picard_tol);
    const int                   picard_max_iter = (opts.picard_max_iter > 0)
                                                  ? opts.picard_max_iter : 50;

    // Boundary codes in NESW order; pad with 'o' (fixed elevation) if short.
    auto bc_code = [&](int i) -> char
    {
        return (static_cast<int>(opts.boundary_conditions.size()) > i)
               ? opts.boundary_conditions[static_cast<std::size_t>(i)] : 'o';
    };
    const char cN = bc_code(0), cE = bc_code(1), cS = bc_code(2), cW = bc_code(3);

    auto idx = [&](int r, int c) -> std::int64_t
    {
        return static_cast<std::int64_t>(r) * n_cols + c;
    };

    // Working field (double) — Dirichlet cells keep their (fixed) value throughout.
    std::vector<double> z(static_cast<std::size_t>(n_cells), 0.0);
    {
        const std::vector<float>& zin = dem.get_data();
        for (std::int64_t i = 0; i < n_cells; ++i)
            z[static_cast<std::size_t>(i)] = static_cast<double>(zin[static_cast<std::size_t>(i)]);
    }

    auto is_o_edge = [&](int r, int c) -> bool
    {
        return (r == 0          && cN == 'o') || (r == n_rows - 1 && cS == 'o')
            || (c == 0          && cW == 'o') || (c == n_cols - 1 && cE == 'o');
    };
    auto is_ch_at = [&](std::int64_t flat) -> bool
    {
        const std::size_t s = static_cast<std::size_t>(flat);
        return (s < is_channel.size()) && (is_channel[s] != 0);   // robust to empty/short
    };

    // Classify: Dirichlet (valid AND (channel OR 'o'-edge)) vs DOF (valid, not Dirichlet).
    std::vector<std::int64_t> flat_to_dof(static_cast<std::size_t>(n_cells), -1);
    std::int64_t n_dof = 0;
    for (int r = 0; r < n_rows; ++r)
        for (int c = 0; c < n_cols; ++c)
        {
            const std::int64_t f = idx(r, c);
            const std::size_t  s = static_cast<std::size_t>(f);
            if (!mask[s]) continue;
            const bool dirichlet = is_ch_at(f) || is_o_edge(r, c);
            if (!dirichlet) flat_to_dof[s] = n_dof++;
        }

    if (n_dof == 0 || dt_total <= 0.0f)
    {
        // Nothing to solve — return the input field unchanged.
        result.elevation           = dem;
        result.n_substeps          = 0;
        result.n_rejected_substeps = 0;
        result.total_picard_iters  = 0;
        result.final_dt            = dt_total;
        result.max_slope_ratio     = compute_max_slope_ratio(dem, is_channel, opts.Sc);
        result.limiter_fired       = false;
        return result;
    }

    // Reflecting/periodic neighbour value (for cross-slope stencils on the current iterate).
    auto nbr_z = [&](int r, int c, int dr, int dc) -> double
    {
        int rr = r + dr, cc = c + dc;
        if (rr < 0)        rr = (cN == 'p') ? n_rows - 1 : r;
        if (rr >= n_rows)  rr = (cS == 'p') ? 0          : r;
        if (cc < 0)        cc = (cW == 'p') ? n_cols - 1 : c;
        if (cc >= n_cols)  cc = (cE == 'p') ? 0          : c;
        if (!mask[static_cast<std::size_t>(idx(rr, cc))]) return z[static_cast<std::size_t>(idx(r, c))];
        return z[static_cast<std::size_t>(idx(rr, cc))];
    };

    // Cardinal face: returns true and sets (nr,nc) of a real neighbour cell, or
    // false for a no-flux face (OOB on 'n'/'o', or masked neighbour). dir: 0=N,1=E,2=S,3=W.
    auto get_face = [&](int r, int c, int dir, int& nr, int& nc) -> bool
    {
        int dr = 0, dc = 0; char code = 'o';
        if      (dir == 0) { dr = -1; code = cN; }
        else if (dir == 1) { dc =  1; code = cE; }
        else if (dir == 2) { dr =  1; code = cS; }
        else               { dc = -1; code = cW; }
        int rr = r + dr, cc = c + dc;
        if (rr < 0 || rr >= n_rows || cc < 0 || cc >= n_cols)
        {
            if (code != 'p') return false;               // no-flux ('n'/'o' OOB)
            rr = (rr + n_rows) % n_rows;
            cc = (cc + n_cols) % n_cols;
        }
        if (!mask[static_cast<std::size_t>(idx(rr, cc))]) return false;  // Neumann at mask
        nr = rr; nc = cc; return true;
    };

    bool limiter_fired = false;
    // Face-centred effective diffusivity (§3); symmetric in the two adjacent cells.
    auto face_Deff = [&](int r, int c, int nr, int nc, int dir) -> double
    {
        if (!nonlinear) return D;
        const double Snorm = (z[static_cast<std::size_t>(idx(nr, nc))]
                            - z[static_cast<std::size_t>(idx(r, c))]) / h;
        double Scross;
        if (dir == 0 || dir == 2)   // N/S face -> cross is the E-W gradient
        {
            const double gP = (nbr_z(r, c, 0, 1)  - nbr_z(r, c, 0, -1))  / (2.0 * h);
            const double gQ = (nbr_z(nr, nc, 0, 1) - nbr_z(nr, nc, 0, -1)) / (2.0 * h);
            Scross = 0.5 * (gP + gQ);
        }
        else                        // E/W face -> cross is the N-S gradient
        {
            const double gP = (nbr_z(r, c, -1, 0)  - nbr_z(r, c, 1, 0))  / (2.0 * h);
            const double gQ = (nbr_z(nr, nc, -1, 0) - nbr_z(nr, nc, 1, 0)) / (2.0 * h);
            Scross = 0.5 * (gP + gQ);
        }
        double ratio2 = (Snorm * Snorm + Scross * Scross) / (Sc * Sc);
        const double cap = r_max * r_max;
        if (ratio2 > cap) { ratio2 = cap; limiter_fired = true; }
        return D / (1.0 - ratio2);
    };

    // Fixed sub-stepping over dt_total (adaptive controller is Phase 3).
    double dt_sub = (opts.dt_init > 0.0f) ? static_cast<double>(opts.dt_init)
                                          : static_cast<double>(dt_total);
    int n_sub = static_cast<int>(std::ceil(static_cast<double>(dt_total) / dt_sub));
    if (n_sub < 1) n_sub = 1;
    dt_sub = static_cast<double>(dt_total) / n_sub;

    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
    std::vector<double> z_old(static_cast<std::size_t>(n_cells));
    int total_picard = 0;

    for (int sub = 0; sub < n_sub; ++sub)
    {
        z_old = z;                                  // start of step
        for (int k = 0; k < picard_max_iter; ++k)
        {
            std::vector<Eigen::Triplet<double>> trips;
            trips.reserve(static_cast<std::size_t>(n_dof) * 5);
            Eigen::VectorXd b(n_dof);

            for (int r = 0; r < n_rows; ++r)
                for (int c = 0; c < n_cols; ++c)
                {
                    const std::int64_t f = idx(r, c);
                    const std::int64_t row = flat_to_dof[static_cast<std::size_t>(f)];
                    if (row < 0) continue;           // skip Dirichlet / masked

                    double diag = 1.0 / dt_sub;
                    double rhs  = z_old[static_cast<std::size_t>(f)] / dt_sub + U;

                    for (int dir = 0; dir < 4; ++dir)
                    {
                        int nr = 0, nc = 0;
                        if (!get_face(r, c, dir, nr, nc)) continue;   // no-flux face
                        const double Deff = face_Deff(r, c, nr, nc, dir);
                        const double w    = Deff / h2;
                        diag += w;
                        const std::int64_t ncol =
                            flat_to_dof[static_cast<std::size_t>(idx(nr, nc))];
                        if (ncol >= 0)
                            trips.emplace_back(static_cast<int>(row), static_cast<int>(ncol), -w);
                        else                                           // Dirichlet neighbour
                            rhs += w * z[static_cast<std::size_t>(idx(nr, nc))];
                    }
                    trips.emplace_back(static_cast<int>(row), static_cast<int>(row), diag);
                    b(static_cast<int>(row)) = rhs;
                }

            Eigen::SparseMatrix<double> A(static_cast<int>(n_dof), static_cast<int>(n_dof));
            A.setFromTriplets(trips.begin(), trips.end());
            solver.compute(A);
            const Eigen::VectorXd zstar = solver.solve(b);
            ++total_picard;

            // Damped update of the DOF cells; converge on max|Δz|.
            double max_dz = 0.0;
            for (std::int64_t f = 0; f < n_cells; ++f)
            {
                const std::int64_t row = flat_to_dof[static_cast<std::size_t>(f)];
                if (row < 0) continue;
                const double znew = z[static_cast<std::size_t>(f)]
                    + omega * (zstar(static_cast<int>(row)) - z[static_cast<std::size_t>(f)]);
                max_dz = std::max(max_dz, std::abs(znew - z[static_cast<std::size_t>(f)]));
                z[static_cast<std::size_t>(f)] = znew;
            }
            if (!nonlinear || max_dz < picard_tol) break;   // linear: one solve is exact
        }
    }

    // Pack the result.
    std::vector<float>   zf(static_cast<std::size_t>(n_cells));
    std::vector<uint8_t> mcopy = mask;
    for (std::int64_t i = 0; i < n_cells; ++i)
    {
        const std::size_t s = static_cast<std::size_t>(i);
        zf[s] = mask[s] ? static_cast<float>(z[s]) : meta.no_data_value;
    }
    result.elevation           = lsdtt_raster_cpu(std::move(zf), std::move(mcopy), meta);
    result.n_substeps          = n_sub;
    result.n_rejected_substeps = 0;     // Phase 3
    result.total_picard_iters  = total_picard;
    result.final_dt            = dt_sub;
    result.max_slope_ratio     = compute_max_slope_ratio(result.elevation, is_channel, opts.Sc);
    result.limiter_fired       = limiter_fired;

#else
    // -------------------------------------------------------------------------
    // No Eigen — fall back to the explicit subcycling solver so callers need no
    // #ifdef (§8 recommendation).
    // -------------------------------------------------------------------------
    std::cout << "[roering_implicit] Eigen unavailable; using the explicit "
                 "subcycling fallback (compute_roering_hillslope_transport)."
              << std::endl;

    lsdtt_raster_cpu elev = compute_roering_hillslope_transport(
        dem, is_channel, opts.D, dt_total, opts.nonlinear, opts.Sc, opts.n_threads);

    result.elevation           = std::move(elev);
    result.n_substeps          = 0;     // explicit fallback does not report its substep count
    result.n_rejected_substeps = 0;
    result.total_picard_iters  = 0;     // not applicable to the explicit scheme
    result.final_dt            = dt_total;
    result.max_slope_ratio     = compute_max_slope_ratio(result.elevation, is_channel, opts.Sc);
    result.limiter_fired       = false;
#endif

    return result;
}

} // namespace lsdtt
