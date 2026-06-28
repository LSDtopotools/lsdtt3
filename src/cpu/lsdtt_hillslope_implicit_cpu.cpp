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
#include <limits>
#include <vector>

#include "lsdtt_hillslope_transport_cpu.hpp"

#ifdef _OPENMP
#  include <omp.h>
#endif

#ifdef LSDTT_HAS_EIGEN
#  include <Eigen/Sparse>
#  include <Eigen/IterativeLinearSolvers>
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
    // Phase 3 (done): adaptive iteration-count-based sub-stepping controller (§7)
    // with post-step face-slope projection / step rejection (§5 layer 3) and
    // omega auto-tuning on Picard divergence.
    // Phase 4 (done): the sparsity pattern is fixed for the whole call (it depends
    // only on the mask / boundary codes / DOF map, not on z), so analyzePattern()
    // is called once and only factorize() runs each Picard iteration; the triplet
    // assembly is OpenMP-parallel (per-thread buffers, then concatenate); a
    // CG + IncompleteCholesky backend (solver_choice == 1) with warm-start is
    // available for large grids; and a semi-implicit fast path takes a single
    // lagged-coefficient solve when the field is far below Sc.
    // TODO(Phase 6): optional Newton / TR-BDF2 cores.
    // -------------------------------------------------------------------------
#ifdef _OPENMP
    if (opts.n_threads > 0) omp_set_num_threads(opts.n_threads);
#endif
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

    // Squared face slope ratio (|∇z|²/Sc², uncapped) using the face-centred
    // cross-slope of Perron (2011) (§3); symmetric in the two adjacent cells.
    auto face_ratio2 = [&](int r, int c, int nr, int nc, int dir) -> double
    {
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
        return (Snorm * Snorm + Scross * Scross) / (Sc * Sc);
    };

    // Face-centred effective diffusivity (§3) — the r_max cap (limiter layer 1)
    // keeps D_eff finite.  The cap firing here is benign (§5/§9) so it does NOT
    // set the result limiter_fired flag, which is reserved for the post-step
    // projection (layer 3).
    auto face_Deff = [&](int r, int c, int nr, int nc, int dir) -> double
    {
        if (!nonlinear) return D;
        double ratio2 = face_ratio2(r, c, nr, nc, dir);
        const double cap = r_max * r_max;
        if (ratio2 > cap) ratio2 = cap;
        return D / (1.0 - ratio2);
    };

    // Maximum (uncapped) face slope ratio max|∇z|/Sc over all faces of valid
    // cells, in the face-centred convention used by the solver (§3).  This is the
    // quantity the limiter and the result diagnostic are defined against — it
    // differs from the per-cell centred-difference compute_max_slope_ratio helper.
    auto max_face_ratio = [&]() -> double
    {
        double mx = 0.0;
        for (int r = 0; r < n_rows; ++r)
            for (int c = 0; c < n_cols; ++c)
            {
                if (!mask[static_cast<std::size_t>(idx(r, c))]) continue;
                for (int dir = 0; dir < 4; ++dir)
                {
                    int nr = 0, nc = 0;
                    if (!get_face(r, c, dir, nr, nc)) continue;
                    mx = std::max(mx, std::sqrt(face_ratio2(r, c, nr, nc, dir)));
                }
            }
        return mx;
    };

    // Linear algebra backends (§6).  The sparsity pattern is fixed for the whole
    // call, so the direct solver's symbolic analysis is done once (analyzePattern)
    // and reused; only factorize() runs per Picard iteration.  solver_choice == 1
    // selects CG + IncompleteCholesky (warm-started) for large grids where direct
    // factor fill-in is prohibitive.
    const int solver_choice = opts.solver_choice;
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> ldlt;
    Eigen::ConjugateGradient<Eigen::SparseMatrix<double>,
                             Eigen::Lower | Eigen::Upper,
                             Eigen::IncompleteCholesky<double>> cg;
    bool pattern_analyzed = false;
    Eigen::SparseMatrix<double> A(static_cast<int>(n_dof), static_cast<int>(n_dof));

#ifdef _OPENMP
    const int n_assembly_threads = std::max(1, omp_get_max_threads());
#else
    const int n_assembly_threads = 1;
#endif

    // One backward-Euler sub-step over dt with a damped Picard iteration from the
    // start-of-step field z_old (z is the iterate, advanced in place).  Returns
    // whether Picard converged within iter_cap; reports the iteration count and the
    // max |z − z_old| change accumulated over the step.  When semi_implicit is set
    // (gentle field, §6 fast path) a single lagged-coefficient solve is taken and
    // treated as converged.
    auto attempt_step = [&](double dt, const std::vector<double>& z_old,
                            double cur_omega, int& iters_out,
                            double& step_max_dz, int iter_cap,
                            bool semi_implicit) -> bool
    {
        iters_out   = 0;
        step_max_dz = 0.0;
        bool converged = false;
        for (int k = 0; k < iter_cap; ++k)
        {
            // Per-thread triplet buffers, concatenated after the parallel region
            // (the matrix pattern is identical every iteration; b(row) writes are to
            // unique rows so they need no synchronisation).
            std::vector<std::vector<Eigen::Triplet<double>>> tl_trips(
                static_cast<std::size_t>(n_assembly_threads));
            Eigen::VectorXd b(n_dof);

            #pragma omp parallel
            {
#ifdef _OPENMP
                const int tid = omp_get_thread_num();
#else
                const int tid = 0;
#endif
                std::vector<Eigen::Triplet<double>>& trips =
                    tl_trips[static_cast<std::size_t>(tid)];
                trips.reserve(static_cast<std::size_t>(n_dof) * 5
                              / static_cast<std::size_t>(n_assembly_threads) + 16);
                #pragma omp for schedule(static)
                for (int r = 0; r < n_rows; ++r)
                    for (int c = 0; c < n_cols; ++c)
                    {
                        const std::int64_t f = idx(r, c);
                        const std::int64_t row = flat_to_dof[static_cast<std::size_t>(f)];
                        if (row < 0) continue;           // skip Dirichlet / masked

                        double diag = 1.0 / dt;
                        double rhs  = z_old[static_cast<std::size_t>(f)] / dt + U;

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
                            else                                          // Dirichlet neighbour
                                rhs += w * z[static_cast<std::size_t>(idx(nr, nc))];
                        }
                        trips.emplace_back(static_cast<int>(row), static_cast<int>(row), diag);
                        b(static_cast<int>(row)) = rhs;
                    }
            }

            std::size_t n_trips = 0;
            for (const auto& v : tl_trips) n_trips += v.size();
            std::vector<Eigen::Triplet<double>> trips;
            trips.reserve(n_trips);
            for (auto& v : tl_trips) trips.insert(trips.end(), v.begin(), v.end());

            A.setFromTriplets(trips.begin(), trips.end());

            Eigen::VectorXd zstar;
            if (solver_choice == 1)
            {
                // CG + IncompleteCholesky, warm-started from the current iterate.
                cg.compute(A);
                Eigen::VectorXd guess(n_dof);
                for (std::int64_t g = 0; g < n_cells; ++g)
                {
                    const std::int64_t row = flat_to_dof[static_cast<std::size_t>(g)];
                    if (row >= 0) guess(static_cast<int>(row)) = z[static_cast<std::size_t>(g)];
                }
                zstar = cg.solveWithGuess(b, guess);
            }
            else
            {
                // Direct: reuse the symbolic factorisation; refactor numerically.
                if (!pattern_analyzed) { ldlt.analyzePattern(A); pattern_analyzed = true; }
                ldlt.factorize(A);
                zstar = ldlt.solve(b);
            }
            ++iters_out;

            // Damped update of the DOF cells; converge on max|Δz| between iterates.
            double max_dz = 0.0;
            for (std::int64_t f = 0; f < n_cells; ++f)
            {
                const std::int64_t row = flat_to_dof[static_cast<std::size_t>(f)];
                if (row < 0) continue;
                const double znew = z[static_cast<std::size_t>(f)]
                    + cur_omega * (zstar(static_cast<int>(row)) - z[static_cast<std::size_t>(f)]);
                max_dz = std::max(max_dz, std::abs(znew - z[static_cast<std::size_t>(f)]));
                z[static_cast<std::size_t>(f)] = znew;
            }
            if (!nonlinear)    { converged = true; break; }  // linear: one solve is exact
            if (semi_implicit) { converged = true; break; }  // gentle field fast path
            if (max_dz < picard_tol) { converged = true; break; }
        }
        for (std::int64_t f = 0; f < n_cells; ++f)
        {
            if (flat_to_dof[static_cast<std::size_t>(f)] < 0) continue;
            step_max_dz = std::max(step_max_dz,
                std::abs(z[static_cast<std::size_t>(f)] - z_old[static_cast<std::size_t>(f)]));
        }
        return converged;
    };

    // Last-resort slope-restoration projection (limiter layer 3, §5): clamp DOF
    // cells by Gauss-Seidel sweeps so no cardinal face-normal slope exceeds
    // r_max·Sc.  Only lowers cells (over-steepening from uplift makes ridge cells
    // too high relative to downslope/channel neighbours).
    auto apply_projection = [&]()
    {
        const double cap_slope = r_max * Sc;
        for (int sweep = 0; sweep < 25; ++sweep)
        {
            double max_change = 0.0;
            for (int r = 0; r < n_rows; ++r)
                for (int c = 0; c < n_cols; ++c)
                {
                    const std::int64_t f = idx(r, c);
                    if (flat_to_dof[static_cast<std::size_t>(f)] < 0) continue;
                    double z_allowed = std::numeric_limits<double>::infinity();
                    for (int dir = 0; dir < 4; ++dir)
                    {
                        int nr = 0, nc = 0;
                        if (!get_face(r, c, dir, nr, nc)) continue;
                        z_allowed = std::min(z_allowed,
                            z[static_cast<std::size_t>(idx(nr, nc))] + cap_slope * h);
                    }
                    if (z[static_cast<std::size_t>(f)] > z_allowed)
                    {
                        max_change = std::max(max_change,
                                              z[static_cast<std::size_t>(f)] - z_allowed);
                        z[static_cast<std::size_t>(f)] = z_allowed;
                    }
                }
            if (max_change < 1e-12) break;
        }
    };

    // -------------------------------------------------------------------------
    // Time integration over dt_total.
    // -------------------------------------------------------------------------
    std::vector<double> z_old(static_cast<std::size_t>(n_cells));
    bool   limiter_fired = false;
    int    n_sub         = 0;
    int    n_rejected    = 0;
    int    total_picard  = 0;
    double final_dt      = 0.0;

    if (!nonlinear)
    {
        // Linear backward-Euler: fixed sub-stepping (one exact solve per step), no
        // limiter / rejection (no Sc constraint in the linear law).
        double dt_lin = (opts.dt_init > 0.0f)
                        ? std::min(static_cast<double>(opts.dt_init), static_cast<double>(dt_total))
                        : static_cast<double>(dt_total);
        int n_lin = static_cast<int>(std::ceil(static_cast<double>(dt_total) / dt_lin));
        if (n_lin < 1) n_lin = 1;
        dt_lin = static_cast<double>(dt_total) / n_lin;
        for (int s = 0; s < n_lin; ++s)
        {
            z_old = z;
            int iters = 0; double step_dz = 0.0;
            attempt_step(dt_lin, z_old, 1.0, iters, step_dz, picard_max_iter, false);
            total_picard += iters;
            ++n_sub;
        }
        final_dt = dt_lin;
    }
    else
    {
        // Adaptive, iteration-count-based controller (§7) + post-step projection /
        // rejection (§5 layer 3).
        const double dt_max_eff = (opts.dt_max > 0.0f) ? static_cast<double>(opts.dt_max)
                                                       : static_cast<double>(dt_total);
        double dt_min_eff = (opts.dt_min > 0.0f) ? static_cast<double>(opts.dt_min)
                                                 : static_cast<double>(dt_total) * 1e-4;
        dt_min_eff = std::min(dt_min_eff, dt_max_eff);
        const double grow   = (opts.grow_factor   > 1.0f) ? static_cast<double>(opts.grow_factor)   : 1.5;
        const double shrink = (opts.shrink_factor > 0.0f && opts.shrink_factor < 1.0f)
                              ? static_cast<double>(opts.shrink_factor) : 0.5;
        const double dz_target = static_cast<double>(opts.dz_target);   // 0 disables
        const double reject_ratio = r_max;          // accept needs max face ratio < r_max
        const int    k_lo = 5;                       // "easy" Picard convergence
        const double omega_floor = 0.1;
        // Semi-implicit fast path (§6): when the field is well below Sc the
        // coefficient lag is harmless, so a single lagged solve suffices.  As
        // slopes steepen toward Sc the controller drops back to full Picard.
        const double semi_implicit_ratio = 0.5;

        double cur_omega = omega;
        double dt = (opts.dt_init > 0.0f) ? static_cast<double>(opts.dt_init)
                                          : static_cast<double>(dt_total);
        dt = std::min(dt, std::min(dt_max_eff, static_cast<double>(dt_total)));
        double t_remaining = static_cast<double>(dt_total);
        const double t_eps = static_cast<double>(dt_total) * 1e-12;
        const int    max_attempts = 2000000;
        int          attempts = 0;

        while (t_remaining > t_eps && attempts < max_attempts)
        {
            ++attempts;
            dt = std::min(dt, t_remaining);
            dt = std::min(dt, dt_max_eff);
            if (dt < dt_min_eff) dt = dt_min_eff;

            z_old = z;
            // Use the semi-implicit single-solve path only when the field is gentle
            // *and* damping is undisturbed (cur_omega == omega); a tightened omega
            // signals a recent divergence, so fall back to full Picard.
            const bool semi = (max_face_ratio() < semi_implicit_ratio)
                              && (cur_omega >= omega - 1e-12);
            const int  iter_cap = semi ? 1 : picard_max_iter;
            // A semi-implicit step is a single lagged-coefficient (IMEX) solve, not a
            // fixed-point iteration, so it must be applied IN FULL (omega = 1).  Using
            // the damped omega here would leave the field at only omega of the solve —
            // an under-relaxed, under-diffused step (which biased coupled LEM steady
            // states high).  Damping only makes sense across multiple Picard iterates.
            const double step_omega = semi ? 1.0 : cur_omega;
            int iters = 0; double step_dz = 0.0;
            const bool converged = attempt_step(dt, z_old, step_omega, iters, step_dz,
                                                iter_cap, semi);
            const double mratio = max_face_ratio();
            const bool slope_ok = (mratio < reject_ratio);
            const bool dz_ok    = (dz_target <= 0.0) || (step_dz <= dz_target);
            const bool accept   = converged && slope_ok && dz_ok;

            if (accept)
            {
                t_remaining -= dt; ++n_sub; total_picard += iters; final_dt = dt;
                const bool easy = (iters <= k_lo)
                    && ((dz_target <= 0.0) || (step_dz <= 0.5 * dz_target));
                if (easy) dt = std::min(dt * grow, dt_max_eff);
                cur_omega = std::min(omega, cur_omega * 1.1);   // relax damping back
            }
            else if (dt <= dt_min_eff * (1.0 + 1e-9))
            {
                // At the floor and still failing: force progress.  Keep the
                // attempt's field, shave any super-critical slopes (layer 3), and
                // flag that the projection fired.
                apply_projection();
                limiter_fired = true;
                t_remaining -= dt; ++n_sub; total_picard += iters; final_dt = dt;
            }
            else
            {
                // Reject: discard the attempt, shrink dt, retry.  Tighten damping
                // when the failure was a Picard divergence.
                z = z_old;
                ++n_rejected;
                dt = std::max(dt * shrink, dt_min_eff);
                if (!converged) cur_omega = std::max(omega_floor, cur_omega * 0.7);
            }
        }
        if (final_dt <= 0.0) final_dt = dt;
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
    result.n_rejected_substeps = n_rejected;
    result.total_picard_iters  = total_picard;
    result.final_dt            = static_cast<float>(final_dt);
    result.max_slope_ratio     = static_cast<float>(max_face_ratio());
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
