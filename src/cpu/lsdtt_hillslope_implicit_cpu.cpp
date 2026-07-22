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
// The function is callable and returns a valid, filled
// lsdtt_roering_implicit_result regardless of whether Eigen is present.  With
// Eigen it runs the backward-Euler + damped-Picard core (and the TR-BDF2
// variant) described below; without Eigen it falls back to the explicit
// subcycling solver and prints a one-line notice, so callers need no #ifdef.
//
// Design and derivation: docs/developer_notes/lsdtt3_architecture.md
// (hillslope transport section); validation: docs/testing_and_benchmarking.md.
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

// Compute the maximum |gradz|/Sc over the hillslope (non-channel, valid) cells of a
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
    // Eigen path - backward-Euler + damped Picard core (Sec4).
    //
    // Phases 1 & 2: conservative cell-centred finite-volume backward Euler over
    // fixed sub-steps, with a damped Picard (lagged-coefficient) iteration for
    // the nonlinear face diffusivity (Sec3).  Linear mode (nonlinear == false) uses
    // D_eff = D and converges in one Picard pass.  Dirichlet sinks come from
    // is_channel AND the 'o' (fixed-elevation) boundary code; 'n' = no-flux,
    // 'p' = periodic (Sec3).  Uplift U is folded into the RHS for hillslope DOFs.
    // The system is SPD (symmetric face-centred D_eff) and solved with
    // Eigen::SimplicialLDLT.
    //
    // Phase 3 (done): adaptive iteration-count-based sub-stepping controller (Sec7)
    // with post-step face-slope projection / step rejection (Sec5 layer 3) and
    // omega auto-tuning on Picard divergence.
    // Phase 4 (done): the sparsity pattern is fixed for the whole call (it depends
    // only on the mask / boundary codes / DOF map, not on z), so the pattern is
    // assembled and analyzePattern() run exactly once; each Picard iteration
    // rewrites A.valuePtr() in place from face-diffusivity arrays built off
    // precomputed cell gradients (each face evaluated once), and only factorize()
    // repeats.  A CG + IncompleteCholesky backend (solver_choice == 1) with
    // warm-start is available for large grids; a semi-implicit fast path takes a
    // single lagged-coefficient solve when the field is far below Sc; the linear
    // fixed-dt path reuses one factorisation across all sub-steps.
    // Phase 6 (done, TR-BDF2 only): time_integrator == 1 selects a 2nd-order,
    // L-stable TR-BDF2 composite (trapezoidal stage over gamma*dt + BDF2 stage over
    // (1-gamma)*dt, gamma = 2-sqrt2), reusing the SPD Picard stage solver, with an embedded
    // local-error estimate driving a PI step-size controller for time-accuracy.
    // The Newton core remains deferred (Picard converges in a few iterations, so
    // its 9-point nonsymmetric Jacobian buys little over the TR-BDF2 accuracy path).
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

    // Working field (double) - Dirichlet cells keep their (fixed) value throughout.
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
        // Nothing to solve - return the input field unchanged.
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

    // -------------------------------------------------------------------------
    // Per-call invariant structure.  The DOF list, face topology and sparsity
    // pattern depend only on the mask, boundary codes and Dirichlet set - never
    // on z - so they are built exactly once per call.  Each Picard iteration
    // then only refreshes the face coefficients and writes them straight into
    // the matrix's value array: no triplet lists, no setFromTriplets re-sort,
    // no per-iteration allocation.
    // -------------------------------------------------------------------------

    // Compact DOF list (dof -> flat) so the hot loops scan only unknowns.
    std::vector<std::int64_t> dof_flat(static_cast<std::size_t>(n_dof));
    for (std::int64_t f = 0; f < n_cells; ++f)
    {
        const std::int64_t d = flat_to_dof[static_cast<std::size_t>(f)];
        if (d >= 0) dof_flat[static_cast<std::size_t>(d)] = f;
    }

    // E/S face topology over ALL valid cells: each face is owned once, by the
    // cell on its N/W side (periodic wrap faces by the high-index cell).  -1 =
    // no-flux ('n'/'o' edge or masked neighbour).  Faces of Dirichlet cells are
    // included so a DOF cell reaches a shared face through its owner.
    std::vector<std::int64_t> faceE_nbr(static_cast<std::size_t>(n_cells), -1);
    std::vector<std::int64_t> faceS_nbr(static_cast<std::size_t>(n_cells), -1);
    for (int r = 0; r < n_rows; ++r)
        for (int c = 0; c < n_cols; ++c)
        {
            const std::size_t s = static_cast<std::size_t>(idx(r, c));
            if (!mask[s]) continue;
            int nr = 0, nc = 0;
            if (get_face(r, c, 1, nr, nc)) faceE_nbr[s] = idx(nr, nc);
            if (get_face(r, c, 2, nr, nc)) faceS_nbr[s] = idx(nr, nc);
        }

    // Per-DOF face tables in N,E,S,W order: neighbour flat index (-1 = no-flux
    // face) and the owning slot in the combined face-diffusivity array (E faces
    // at [f], S faces at [n_cells + f]).
    std::vector<std::int64_t> fnbr (static_cast<std::size_t>(4 * n_dof), -1);
    std::vector<std::int64_t> fdeff(static_cast<std::size_t>(4 * n_dof), -1);
    for (std::int64_t d = 0; d < n_dof; ++d)
    {
        const std::int64_t f = dof_flat[static_cast<std::size_t>(d)];
        const int r = static_cast<int>(f / n_cols);
        const int c = static_cast<int>(f % n_cols);
        for (int dir = 0; dir < 4; ++dir)
        {
            int nr = 0, nc = 0;
            if (!get_face(r, c, dir, nr, nc)) continue;
            const std::int64_t nf = idx(nr, nc);
            const std::size_t  t  = static_cast<std::size_t>(4 * d + dir);
            fnbr[t] = nf;
            if      (dir == 1) fdeff[t] = f;              // own E face
            else if (dir == 3) fdeff[t] = nf;             // W = neighbour's E face
            else if (dir == 2) fdeff[t] = n_cells + f;    // own S face
            else               fdeff[t] = n_cells + nf;   // N = neighbour's S face
        }
    }

    // Linear algebra backends (Sec6) and the fixed sparsity pattern.  The pattern
    // is assembled from triplets exactly once (and the direct solver's symbolic
    // analysis done on it once); per-iteration assembly rewrites A.valuePtr() in
    // place.  slot_face holds each DOF face's (nbr_dof, dof) slot in the value
    // array (-1 when the neighbour is Dirichlet, whose contribution goes to the
    // RHS instead).  solver_choice == 1 selects CG + IncompleteCholesky
    // (warm-started) for large grids where direct factor fill-in is prohibitive.
    const int solver_choice = opts.solver_choice;
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> ldlt;
    Eigen::ConjugateGradient<Eigen::SparseMatrix<double>,
                             Eigen::Lower | Eigen::Upper,
                             Eigen::IncompleteCholesky<double>> cg;
    Eigen::SparseMatrix<double> A(static_cast<int>(n_dof), static_cast<int>(n_dof));
    std::vector<std::int64_t> slot_diag(static_cast<std::size_t>(n_dof), -1);
    std::vector<std::int64_t> slot_face(static_cast<std::size_t>(4 * n_dof), -1);
    {
        std::vector<Eigen::Triplet<double>> pat;
        pat.reserve(static_cast<std::size_t>(5 * n_dof));
        for (std::int64_t d = 0; d < n_dof; ++d)
        {
            pat.emplace_back(static_cast<int>(d), static_cast<int>(d), 1.0);
            for (int dir = 0; dir < 4; ++dir)
            {
                const std::int64_t nf = fnbr[static_cast<std::size_t>(4 * d + dir)];
                if (nf < 0) continue;
                const std::int64_t nd = flat_to_dof[static_cast<std::size_t>(nf)];
                if (nd >= 0)
                    pat.emplace_back(static_cast<int>(nd), static_cast<int>(d), 1.0);
            }
        }
        // setFromTriplets sums duplicates, so a 2-cell periodic wrap (E and W
        // faces sharing one neighbour) yields a single entry; update_system
        // accumulates into the slots accordingly.
        A.setFromTriplets(pat.begin(), pat.end());
        A.makeCompressed();

        // Map each (dof, dir) to its valuePtr() slot by walking the column's
        // sorted inner indices.
        const int* outer = A.outerIndexPtr();
        const int* inner = A.innerIndexPtr();
        for (std::int64_t d = 0; d < n_dof; ++d)
            for (int sp = outer[d]; sp < outer[d + 1]; ++sp)
            {
                const int rowi = inner[sp];
                if (rowi == static_cast<int>(d))
                {
                    slot_diag[static_cast<std::size_t>(d)] = sp;
                    continue;
                }
                for (int dir = 0; dir < 4; ++dir)
                {
                    const std::size_t t = static_cast<std::size_t>(4 * d + dir);
                    const std::int64_t nf = fnbr[t];
                    if (nf >= 0
                        && flat_to_dof[static_cast<std::size_t>(nf)] == rowi)
                        slot_face[t] = sp;
                }
            }
        if (solver_choice != 1) ldlt.analyzePattern(A);
    }

    // Per-iteration work arrays, allocated once for the whole call.
    std::vector<double> gx(static_cast<std::size_t>(n_cells), 0.0);
    std::vector<double> gy(static_cast<std::size_t>(n_cells), 0.0);
    std::vector<double> face_deff(static_cast<std::size_t>(2 * n_cells), 0.0);
    std::vector<double> b_dir(static_cast<std::size_t>(n_dof), 0.0);
    std::vector<double> g_field(static_cast<std::size_t>(n_cells), 0.0);
    Eigen::VectorXd b(n_dof), sol(n_dof), guess(n_dof);

    // Cell-centred full gradients of the current iterate (same stencil the
    // per-face lambdas used); feeds the face cross-slope so each cell's
    // gradient is computed once instead of once per adjacent face.
    auto compute_gradients = [&]()
    {
        #pragma omp parallel for schedule(static)
        for (std::int64_t f = 0; f < n_cells; ++f)
        {
            const std::size_t s = static_cast<std::size_t>(f);
            if (!mask[s]) continue;
            const int r = static_cast<int>(f / n_cols);
            const int c = static_cast<int>(f % n_cols);
            gx[s] = (nbr_z(r, c, 0, 1) - nbr_z(r, c, 0, -1)) / (2.0 * h);
            gy[s] = (nbr_z(r, c, -1, 0) - nbr_z(r, c, 1, 0)) / (2.0 * h);
        }
    };

    // Face-centred effective diffusivity (Sec3), each face computed exactly once
    // (cross-slope = mean of the two adjacent cells' transverse gradient, as in
    // Perron 2011; requires compute_gradients() first when nonlinear).  The
    // r_max cap (limiter layer 1) keeps D_eff finite; the cap firing here is
    // benign (Sec5/Sec9) so it does NOT set limiter_fired, which is reserved for
    // the post-step projection (layer 3).
    auto compute_face_coeffs = [&]()
    {
        const double Sc2 = Sc * Sc;
        const double cap = r_max * r_max;
        #pragma omp parallel for schedule(static)
        for (std::int64_t f = 0; f < n_cells; ++f)
        {
            const std::size_t s = static_cast<std::size_t>(f);
            const std::int64_t ne = faceE_nbr[s];
            if (ne >= 0)
            {
                double Deff = D;
                if (nonlinear)
                {
                    const std::size_t sn = static_cast<std::size_t>(ne);
                    const double Snorm  = (z[sn] - z[s]) / h;
                    const double Scross = 0.5 * (gy[s] + gy[sn]);
                    double r2 = (Snorm * Snorm + Scross * Scross) / Sc2;
                    if (r2 > cap) r2 = cap;
                    Deff = D / (1.0 - r2);
                }
                face_deff[s] = Deff;
            }
            const std::int64_t ns_ = faceS_nbr[s];
            if (ns_ >= 0)
            {
                double Deff = D;
                if (nonlinear)
                {
                    const std::size_t sn = static_cast<std::size_t>(ns_);
                    const double Snorm  = (z[sn] - z[s]) / h;
                    const double Scross = 0.5 * (gx[s] + gx[sn]);
                    double r2 = (Snorm * Snorm + Scross * Scross) / Sc2;
                    if (r2 > cap) r2 = cap;
                    Deff = D / (1.0 - r2);
                }
                face_deff[static_cast<std::size_t>(n_cells + f)] = Deff;
            }
        }
    };

    // Maximum (uncapped) SQUARED face slope ratio max(|gradz|/Sc)^2 over all faces
    // of valid cells (face-centred convention, Sec3) - the quantity the limiter
    // accept test and the max_slope_ratio diagnostic are defined against.
    // Callers compare against thresholds squared and only sqrt for reporting.
    auto max_face_ratio2 = [&]() -> double
    {
        compute_gradients();
        const double Sc2 = Sc * Sc;
        double mx = 0.0;
        #pragma omp parallel for schedule(static) reduction(max:mx)
        for (std::int64_t f = 0; f < n_cells; ++f)
        {
            const std::size_t s = static_cast<std::size_t>(f);
            const std::int64_t ne = faceE_nbr[s];
            if (ne >= 0)
            {
                const std::size_t sn = static_cast<std::size_t>(ne);
                const double Snorm  = (z[sn] - z[s]) / h;
                const double Scross = 0.5 * (gy[s] + gy[sn]);
                const double r2 = (Snorm * Snorm + Scross * Scross) / Sc2;
                if (r2 > mx) mx = r2;
            }
            const std::int64_t ns_ = faceS_nbr[s];
            if (ns_ >= 0)
            {
                const std::size_t sn = static_cast<std::size_t>(ns_);
                const double Snorm  = (z[sn] - z[s]) / h;
                const double Scross = 0.5 * (gx[s] + gx[sn]);
                const double r2 = (Snorm * Snorm + Scross * Scross) / Sc2;
                if (r2 > mx) mx = r2;
            }
        }
        return mx;
    };

    // Write the stage matrix [I - beta*L(D_eff)] straight into A's value array
    // (fixed pattern; every slot of column d is written only by DOF d, so the
    // loop is thread-safe) and record the Dirichlet contribution to the RHS in
    // b_dir.  Requires compute_face_coeffs() first.
    auto update_system = [&](double beta)
    {
        double*    av    = A.valuePtr();
        const int* outer = A.outerIndexPtr();
        #pragma omp parallel for schedule(static)
        for (std::int64_t d = 0; d < n_dof; ++d)
        {
            for (int sp = outer[d]; sp < outer[d + 1]; ++sp) av[sp] = 0.0;
            double diag = 1.0;
            double bd   = 0.0;
            for (int dir = 0; dir < 4; ++dir)
            {
                const std::size_t t = static_cast<std::size_t>(4 * d + dir);
                const std::int64_t nf = fnbr[t];
                if (nf < 0) continue;
                const double w =
                    beta * face_deff[static_cast<std::size_t>(fdeff[t])] / h2;
                diag += w;
                const std::int64_t sl = slot_face[t];
                if (sl >= 0) av[sl] -= w;
                else         bd += w * z[static_cast<std::size_t>(nf)];
            }
            av[slot_diag[static_cast<std::size_t>(d)]] += diag;
            b_dir[static_cast<std::size_t>(d)] = bd;
        }
    };

    // Generic implicit stage solver, shared by backward Euler (beta = dt: the
    // [I/dt - L] system scaled through by dt) and both TR-BDF2 stages: a damped
    // Picard iteration on
    //   [ I - beta*L(D_eff) ] y = g   (+ Dirichlet terms)
    // with the nonlinear face diffusivity lagged at the current iterate z, which
    // is advanced in place (warm-started by the caller).  g is a per-cell RHS
    // field (used only at DOF rows).  When reuse_matrix is set (linear fixed-dt
    // sub-stepping, where the operator never changes) the previous factorisation
    // is kept and only the RHS is rebuilt.  Returns whether Picard converged.
    bool   factor_valid = false;
    double factor_beta  = 0.0;
    auto solve_stage = [&](double beta, const std::vector<double>& g,
                           double cur_omega, int iter_cap, int& iters_out,
                           bool reuse_matrix) -> bool
    {
        iters_out = 0;
        bool converged = false;
        for (int k = 0; k < iter_cap; ++k)
        {
            if (!(reuse_matrix && factor_valid && factor_beta == beta))
            {
                if (nonlinear) compute_gradients();
                compute_face_coeffs();
                update_system(beta);
                if (solver_choice == 1) cg.compute(A);
                else                    ldlt.factorize(A);
                factor_valid = true;
                factor_beta  = beta;
            }
            #pragma omp parallel for schedule(static)
            for (std::int64_t d = 0; d < n_dof; ++d)
                b(d) = g[static_cast<std::size_t>(
                             dof_flat[static_cast<std::size_t>(d)])]
                     + b_dir[static_cast<std::size_t>(d)];

            if (solver_choice == 1)
            {
                // CG + IncompleteCholesky, warm-started from the current iterate.
                for (std::int64_t d = 0; d < n_dof; ++d)
                    guess(d) = z[static_cast<std::size_t>(
                                     dof_flat[static_cast<std::size_t>(d)])];
                sol = cg.solveWithGuess(b, guess);
            }
            else
            {
                sol = ldlt.solve(b);
            }
            ++iters_out;

            // Damped update of the DOF cells; converge on max|Deltaz| between iterates.
            double max_dz = 0.0;
            #pragma omp parallel for schedule(static) reduction(max:max_dz)
            for (std::int64_t d = 0; d < n_dof; ++d)
            {
                const std::size_t s = static_cast<std::size_t>(
                    dof_flat[static_cast<std::size_t>(d)]);
                const double znew = z[s] + cur_omega * (sol(d) - z[s]);
                const double dz   = std::abs(znew - z[s]);
                if (dz > max_dz) max_dz = dz;
                z[s] = znew;
            }
            if (!nonlinear)          { converged = true; break; }  // linear: one solve is exact
            if (max_dz < picard_tol) { converged = true; break; }
        }
        return converged;
    };

    // One backward-Euler sub-step over dt from the start-of-step field z_old,
    // via the shared stage solver with beta = dt.  Reports the Picard iteration
    // count and the max |z - z_old| change over the step.  A semi-implicit step
    // (Sec6 fast path) is a single lagged-coefficient IMEX solve treated as
    // converged; the caller passes omega = 1 and iter_cap = 1 for it.
    auto attempt_step = [&](double dt, const std::vector<double>& z_old,
                            double cur_omega, int& iters_out,
                            double& step_max_dz, int iter_cap,
                            bool semi_implicit, bool reuse_matrix) -> bool
    {
        #pragma omp parallel for schedule(static)
        for (std::int64_t d = 0; d < n_dof; ++d)
        {
            const std::size_t s = static_cast<std::size_t>(
                dof_flat[static_cast<std::size_t>(d)]);
            g_field[s] = z_old[s] + dt * U;
        }
        bool converged = solve_stage(dt, g_field, cur_omega, iter_cap, iters_out,
                                     reuse_matrix);
        if (semi_implicit) converged = true;   // gentle field fast path
        double smdz = 0.0;
        #pragma omp parallel for schedule(static) reduction(max:smdz)
        for (std::int64_t d = 0; d < n_dof; ++d)
        {
            const std::size_t s = static_cast<std::size_t>(
                dof_flat[static_cast<std::size_t>(d)]);
            const double dz = std::abs(z[s] - z_old[s]);
            if (dz > smdz) smdz = dz;
        }
        step_max_dz = smdz;
        return converged;
    };

    // -------------------------------------------------------------------------
    // Phase 6: TR-BDF2 building blocks (used only when time_integrator == 1).
    // -------------------------------------------------------------------------

    // Apply the spatial operator L(z)_P = (1/h^2) Sigma_f D_eff,f (z_nb - z_P) at the
    // current field z, over DOF cells (Dirichlet / masked rows left 0).  Used to
    // form the explicit f_n = U + L(z^n) term in the TR-BDF2 stages and error
    // estimate; D_eff is evaluated at the field as it currently stands in z.
    auto compute_L = [&](std::vector<double>& Lout)
    {
        Lout.assign(static_cast<std::size_t>(n_cells), 0.0);
        if (nonlinear) compute_gradients();
        compute_face_coeffs();
        #pragma omp parallel for schedule(static)
        for (std::int64_t d = 0; d < n_dof; ++d)
        {
            const std::int64_t f  = dof_flat[static_cast<std::size_t>(d)];
            const double       zP = z[static_cast<std::size_t>(f)];
            double acc = 0.0;
            for (int dir = 0; dir < 4; ++dir)
            {
                const std::size_t t = static_cast<std::size_t>(4 * d + dir);
                const std::int64_t nf = fnbr[t];
                if (nf < 0) continue;
                const double w =
                    face_deff[static_cast<std::size_t>(fdeff[t])] / h2;
                acc += w * (z[static_cast<std::size_t>(nf)] - zP);
            }
            Lout[static_cast<std::size_t>(f)] = acc;
        }
    };

    // TR-BDF2 coefficients (gamma = 2 - sqrt2): a trapezoidal stage over gamma*dt to z_gamma, then
    // a BDF2 stage over (1-gamma)*dt to z^{n+1}.  Both stages are SPD [I - beta*L]
    // solves dispatched through solve_stage.  err_c scales the embedded local-error
    // estimate (Hosea & Shampine 1996, "Analysis and implementation of TR-BDF2").
    const double tr_gamma = 2.0 - std::sqrt(2.0);
    const double bdf_c1 = 1.0 / (tr_gamma * (2.0 - tr_gamma));
    const double bdf_c2 = (1.0 - tr_gamma) * (1.0 - tr_gamma)
                          / (tr_gamma * (2.0 - tr_gamma));
    const double bdf_c3 = (1.0 - tr_gamma) / (2.0 - tr_gamma);
    const double err_c  = (-3.0 * tr_gamma * tr_gamma + 4.0 * tr_gamma - 2.0)
                          / (12.0 * (2.0 - tr_gamma));

    // One TR-BDF2 step from z^n (which z must hold on entry; z holds z^{n+1} on
    // return).  Produces a 2nd-order, L-stable update plus an embedded local-error
    // estimate err_out (elevation units) from the three stage derivatives f_n, f_gamma,
    // f_{n+1}.  The combination a0*f_n + a1*f_gamma + a2*f_{n+1} is a second divided
    // difference of f (a0+a1+a2 = 0), so it vanishes for constant f and scales as
    // O(dt^2) -> err = O(dt^3), the correct order for a 2nd-order method.  f_gamma and
    // f_{n+1} are recovered from the converged stage relations, so the estimate
    // costs no extra solves or operator applies beyond L(z^n).  Returns whether both
    // stages' Picard iterations converged.
    auto tr_bdf2_step = [&](double dt, const std::vector<double>& z_n,
                            double cur_omega, int iter_cap,
                            int& iters_out, double& err_out,
                            double& step_max_dz) -> bool
    {
        std::vector<double> Lzn;            // L(z^n); z currently holds z^n
        compute_L(Lzn);

        // TR stage: [I - (gamma*dt/2) L] z_gamma = z^n + (gamma*dt/2) L(z^n) + gamma*dt*U
        const double beta_tr = 0.5 * tr_gamma * dt;
        std::vector<double> g(static_cast<std::size_t>(n_cells), 0.0);
        for (std::int64_t f = 0; f < n_cells; ++f)
        {
            const std::size_t s = static_cast<std::size_t>(f);
            if (flat_to_dof[s] < 0) continue;
            g[s] = z_n[s] + beta_tr * Lzn[s] + tr_gamma * dt * U;
        }
        int it1 = 0;
        const bool c_tr = solve_stage(beta_tr, g, cur_omega, iter_cap, it1, false);
        std::vector<double> z_gamma = z;    // stage-1 solution

        // BDF2 stage: [I - c3*dt*L] z^{n+1} = c1*z_gamma - c2*z^n + c3*dt*U
        const double beta_bdf = bdf_c3 * dt;
        for (std::int64_t f = 0; f < n_cells; ++f)
        {
            const std::size_t s = static_cast<std::size_t>(f);
            if (flat_to_dof[s] < 0) continue;
            g[s] = bdf_c1 * z_gamma[s] - bdf_c2 * z_n[s] + bdf_c3 * dt * U;
        }
        int it2 = 0;
        const bool c_bdf = solve_stage(beta_bdf, g, cur_omega, iter_cap, it2, false);
        iters_out = it1 + it2;              // z now holds z^{n+1}

        const double a0 =  1.0 / tr_gamma;
        const double a1 = -1.0 / (tr_gamma * (1.0 - tr_gamma));
        const double a2 =  1.0 / (1.0 - tr_gamma);
        double err = 0.0, smdz = 0.0;
        for (std::int64_t f = 0; f < n_cells; ++f)
        {
            const std::size_t s = static_cast<std::size_t>(f);
            if (flat_to_dof[s] < 0) continue;
            const double fn   = U + Lzn[s];
            const double fg   = 2.0 * (z_gamma[s] - z_n[s]) / (tr_gamma * dt) - fn;
            const double fnp1 = (z[s] - bdf_c1 * z_gamma[s] + bdf_c2 * z_n[s])
                                / (bdf_c3 * dt);
            const double e = err_c * dt * (a0 * fn + a1 * fg + a2 * fnp1);
            err  = std::max(err, std::abs(e));
            smdz = std::max(smdz, std::abs(z[s] - z_n[s]));
        }
        err_out     = err;
        step_max_dz = smdz;
        return c_tr && c_bdf;
    };

    // Last-resort slope-restoration projection (limiter layer 3, Sec5): clamp DOF
    // cells by Gauss-Seidel sweeps so no cardinal face-normal slope exceeds
    // r_max*Sc.  Only lowers cells (over-steepening from uplift makes ridge cells
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

    if (opts.time_integrator == 1)
    {
        // ---------------------------------------------------------------------
        // Phase 6: TR-BDF2 with embedded-error PI step control (Sec7).  Each step
        // is a 2nd-order, L-stable composite solve; the embedded local-error
        // estimate drives a PI controller for genuine time-accuracy.  The same
        // slope-limiter accept/reject and floor projection as the backward-Euler
        // adaptive path guard sub-criticality.  Handles linear and nonlinear
        // (each stage solve simply converges in one Picard pass when linear).
        // ---------------------------------------------------------------------
        const double dt_max_eff = (opts.dt_max > 0.0f) ? static_cast<double>(opts.dt_max)
                                                       : static_cast<double>(dt_total);
        double dt_min_eff = (opts.dt_min > 0.0f) ? static_cast<double>(opts.dt_min)
                                                 : static_cast<double>(dt_total) * 1e-6;
        dt_min_eff = std::min(dt_min_eff, dt_max_eff);
        const double grow   = (opts.grow_factor   > 1.0f) ? static_cast<double>(opts.grow_factor)   : 1.5;
        const double shrink = (opts.shrink_factor > 0.0f && opts.shrink_factor < 1.0f)
                              ? static_cast<double>(opts.shrink_factor) : 0.5;
        const double reject_ratio = r_max;
        const double omega_floor  = 0.1;

        // Error tolerance (elevation units).  Auto: a small fraction of the
        // current relief, floored, when the user leaves error_tol == 0.
        double tol = static_cast<double>(opts.error_tol);
        if (tol <= 0.0)
        {
            double zmin =  std::numeric_limits<double>::infinity();
            double zmax = -std::numeric_limits<double>::infinity();
            for (std::int64_t f = 0; f < n_cells; ++f)
            {
                if (flat_to_dof[static_cast<std::size_t>(f)] < 0) continue;
                zmin = std::min(zmin, z[static_cast<std::size_t>(f)]);
                zmax = std::max(zmax, z[static_cast<std::size_t>(f)]);
            }
            const double relief = (zmax > zmin) ? (zmax - zmin) : 1.0;
            tol = std::max(1e-4 * relief, 1e-9);
        }

        const double exponent = 1.0 / 3.0;     // 1/(order+1); TR-BDF2 is 2nd order
        const double safety   = 0.9;

        double cur_omega = omega;
        double dt = (opts.dt_init > 0.0f) ? static_cast<double>(opts.dt_init)
                                          : static_cast<double>(dt_total);
        dt = std::min(dt, std::min(dt_max_eff, static_cast<double>(dt_total)));
        double t_remaining = static_cast<double>(dt_total);
        const double t_eps = static_cast<double>(dt_total) * 1e-12;
        double err_prev = -1.0;                 // <0 -> no PI memory yet
        const int max_attempts = 2000000;
        int       attempts = 0;

        while (t_remaining > t_eps && attempts < max_attempts)
        {
            ++attempts;
            dt = std::min(dt, t_remaining);
            dt = std::min(dt, dt_max_eff);
            if (dt < dt_min_eff) dt = dt_min_eff;

            z_old = z;                          // z^n
            int iters = 0; double err = 0.0, step_dz = 0.0;
            const bool converged = tr_bdf2_step(dt, z_old, cur_omega, picard_max_iter,
                                                iters, err, step_dz);
            const double mratio2  = max_face_ratio2();
            const bool   slope_ok = (mratio2 < reject_ratio * reject_ratio);
            const bool   err_ok   = (err <= tol);
            const bool   accept   = converged && slope_ok && err_ok;

            // PI step-size factor (Gustafsson) from the embedded error estimate.
            double fac;
            if (err <= 0.0)
                fac = grow;
            else if (err_prev > 0.0)
                fac = safety * std::pow(tol / err, 0.7 * exponent)
                             * std::pow(err_prev / err, 0.4 * exponent);
            else
                fac = safety * std::pow(tol / err, exponent);
            fac = std::max(shrink, std::min(grow, fac));

            if (accept)
            {
                t_remaining -= dt; ++n_sub; total_picard += iters; final_dt = dt;
                err_prev  = err;
                cur_omega = std::min(omega, cur_omega * 1.1);
                dt = std::min(dt * fac, dt_max_eff);
            }
            else if (dt <= dt_min_eff * (1.0 + 1e-9))
            {
                // At the floor and still failing: force progress, shave any
                // super-critical slopes (layer 3), flag the projection.
                apply_projection();
                limiter_fired = true;
                t_remaining -= dt; ++n_sub; total_picard += iters; final_dt = dt;
                if (err > 0.0) err_prev = err;
            }
            else
            {
                // Reject: discard the attempt, shrink dt (always by at least
                // `shrink`), retry; tighten damping after a Picard divergence.
                z = z_old;
                ++n_rejected;
                dt = std::max(dt * std::min(fac, shrink), dt_min_eff);
                if (!converged) cur_omega = std::max(omega_floor, cur_omega * 0.7);
            }
        }
        if (final_dt <= 0.0) final_dt = dt;
    }
    else if (!nonlinear)
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
            // reuse_matrix: the linear operator and dt_lin never change across
            // sub-steps, so the factorisation from step 0 serves all of them.
            attempt_step(dt_lin, z_old, 1.0, iters, step_dz, picard_max_iter,
                         false, true);
            total_picard += iters;
            ++n_sub;
        }
        final_dt = dt_lin;
    }
    else
    {
        // Adaptive, iteration-count-based controller (Sec7) + post-step projection /
        // rejection (Sec5 layer 3).
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
        // Semi-implicit fast path (Sec6): when the field is well below Sc the
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

        // Cached max face ratio^2 for the CURRENT z: between an accepted step's
        // post-check and the next attempt's semi decision z does not change (and
        // a rejected attempt restores z), so the pre-step scan is never repeated.
        double mratio2_cache = max_face_ratio2();

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
            const bool semi = (mratio2_cache < semi_implicit_ratio * semi_implicit_ratio)
                              && (cur_omega >= omega - 1e-12);
            const int  iter_cap = semi ? 1 : picard_max_iter;
            // A semi-implicit step is a single lagged-coefficient (IMEX) solve, not a
            // fixed-point iteration, so it must be applied IN FULL (omega = 1).  Using
            // the damped omega here would leave the field at only omega of the solve -
            // an under-relaxed, under-diffused step (which biased coupled LEM steady
            // states high).  Damping only makes sense across multiple Picard iterates.
            const double step_omega = semi ? 1.0 : cur_omega;
            int iters = 0; double step_dz = 0.0;
            const bool converged = attempt_step(dt, z_old, step_omega, iters, step_dz,
                                                iter_cap, semi, false);
            const double mratio2 = max_face_ratio2();
            const bool slope_ok = (mratio2 < reject_ratio * reject_ratio);
            const bool dz_ok    = (dz_target <= 0.0) || (step_dz <= dz_target);
            const bool accept   = converged && slope_ok && dz_ok;

            if (accept)
            {
                t_remaining -= dt; ++n_sub; total_picard += iters; final_dt = dt;
                mratio2_cache = mratio2;
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
                mratio2_cache = max_face_ratio2();   // projection changed z
                t_remaining -= dt; ++n_sub; total_picard += iters; final_dt = dt;
            }
            else
            {
                // Reject: discard the attempt, shrink dt, retry.  Tighten damping
                // when the failure was a Picard divergence.  mratio2_cache still
                // describes the restored z.
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
    result.max_slope_ratio     = static_cast<float>(std::sqrt(max_face_ratio2()));
    result.limiter_fired       = limiter_fired;

#else
    // -------------------------------------------------------------------------
    // No Eigen - fall back to the explicit subcycling solver so callers need no
    // #ifdef (Sec8 recommendation).
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
