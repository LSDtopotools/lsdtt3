// =============================================================================
// lsdtt_hillslope_newton_cpu.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
//
// Newton solve for the Roering nonlinear hillslope. See the header for why this
// exists and how it relates to the Picard/SOR path.
// =============================================================================

#include "lsdtt_hillslope_newton_cpu.hpp"
#include "lsdtt_hillslope_operator_cpu.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <memory>
#include <vector>

#ifdef _OPENMP
#  include <omp.h>
#endif

#ifdef LSDTT_HAS_EIGEN
#  include <Eigen/Sparse>
#  include <Eigen/IterativeLinearSolvers>
#endif

namespace lsdtt
{

// -----------------------------------------------------------------------------
// lsdtt_hillslope_newton_workspace
//
// Eigen types live in `impl` so the public header stays Eigen-free.
// -----------------------------------------------------------------------------
struct lsdtt_hillslope_newton_workspace::impl
{
#ifdef LSDTT_HAS_EIGEN
    Eigen::SparseMatrix<double>                 J;
    Eigen::SparseLU<Eigen::SparseMatrix<double>> lu;
#endif

    // Validity key. Anything that changes the MATRIX must appear here, or a stale
    // factorisation would be reused and silently produce a wrong answer.
    bool          have_key         = false;
    int           nrows = 0, ncols = 0;
    char          bc[4]            = {0, 0, 0, 0};
    double        D = 0.0, Sc = 0.0, inv_dt = 0.0, h = 0.0;
    bool          nonlinear        = true;
    long int      n_dof            = 0;
    std::uint64_t dof_hash         = 0;   // covers mask + channels + Dirichlet edges

    bool pattern_valid = false;   // analyzePattern done for this key
    bool numeric_valid = false;   // factorize done AND the values cannot have moved

    long int n_pattern_analyses = 0;
    long int n_factorisations   = 0;
    long int n_solves           = 0;

    void invalidate() { pattern_valid = false; numeric_valid = false; }
};

lsdtt_hillslope_newton_workspace::lsdtt_hillslope_newton_workspace()
    : p_(std::make_unique<impl>()) {}
lsdtt_hillslope_newton_workspace::~lsdtt_hillslope_newton_workspace() = default;
lsdtt_hillslope_newton_workspace::lsdtt_hillslope_newton_workspace(
    lsdtt_hillslope_newton_workspace&&) noexcept = default;
lsdtt_hillslope_newton_workspace&
lsdtt_hillslope_newton_workspace::operator=(
    lsdtt_hillslope_newton_workspace&&) noexcept = default;

void lsdtt_hillslope_newton_workspace::reset()
{
    const long int a = p_->n_pattern_analyses, f = p_->n_factorisations,
                   s = p_->n_solves;
    p_ = std::make_unique<impl>();
    p_->n_pattern_analyses = a; p_->n_factorisations = f; p_->n_solves = s;
}

long int lsdtt_hillslope_newton_workspace::n_pattern_analyses() const noexcept
{ return p_->n_pattern_analyses; }
long int lsdtt_hillslope_newton_workspace::n_factorisations() const noexcept
{ return p_->n_factorisations; }
long int lsdtt_hillslope_newton_workspace::n_solves() const noexcept
{ return p_->n_solves; }


#ifdef LSDTT_HAS_EIGEN

namespace
{

// -----------------------------------------------------------------------------
// The residual and Jacobian of
//
//     F_P(z) = inv_dt*(z_P - z_P^n) - U - (D/h^2) sum_f alpha_f(z)(z_nb,f - z_P)
//
// Derivation of the Jacobian. With alpha = 1/(1-r), r = (S_n^2 + S_t^2)/Sc^2:
//
//     d(alpha)/dr        = 1/(1-r)^2 = alpha^2
//     d(alpha_f)/dz_Q    = alpha_f^2 * 2*(S_n dS_n/dz_Q + S_t dS_t/dz_Q) / Sc^2
//
// For the E face between P=(i,j) and E=(i,j+1):
//
//     S_n = (z_E - z_P)/h          dS_n/dz_P = -1/h,  dS_n/dz_E = +1/h
//     S_t = 1/2[(z_{i+1,j}-z_{i-1,j})/(2h) + (z_{i+1,j+1}-z_{i-1,j+1})/(2h)]
//                                  dS_t/dz_{i+1,j}   = +1/(4h)
//                                  dS_t/dz_{i-1,j}   = -1/(4h)
//                                  dS_t/dz_{i+1,j+1} = +1/(4h)
//                                  dS_t/dz_{i-1,j+1} = -1/(4h)
//
// so the tangential term reaches the DIAGONAL neighbours: the Jacobian is
// 9-point. And because d(alpha_f)/dz_Q multiplies (z_nb - z_P), which is
// antisymmetric in the pair, it is also NON-SYMMETRIC.
//
// Where the cap binds, r is clamped and d(alpha)/dz = 0. That is a valid
// subgradient, it bounds the Jacobian, and it is what lets the limiter coexist
// with Newton instead of injecting spurious curvature.
// -----------------------------------------------------------------------------

struct jac_entry
{
    int    r;
    int    c;
    double w;
};

// Accumulate d(alpha_f)/dz_Q for one face into `out`, returning the number of
// entries written (at most 6: the two face cells and four diagonal/tangential
// contributors). `dir`: 0=N, 1=E, 2=S, 3=W.
int dalpha_dz(const lsdtt_hillslope_stencil<double>& st,
              int r, int c, int dir, double Sc, double ratio_cap,
              double alpha, double s_n, jac_entry out[8], bool& capped)
{
    static constexpr int k_dr[4] = { -1,  0, +1,  0 };
    static constexpr int k_dc[4] = {  0, +1,  0, -1 };

    const double h  = st.h;
    const double Sc2 = Sc * Sc;

    int nr = 0, nc = 0;
    st.dom.neighbour(r, c, k_dr[dir], k_dc[dir], nr, nc);

    const bool tangent_is_x = (dir == 0 || dir == 2);
    const double s_t = tangent_is_x
        ? 0.5 * (st.slope_x(r, c) + st.slope_x(nr, nc))
        : 0.5 * (st.slope_y(r, c) + st.slope_y(nr, nc));

    // If the cap binds, alpha is locally constant: zero derivative.
    const double r_raw = (s_n * s_n + s_t * s_t) / Sc2;
    if (r_raw >= ratio_cap) { capped = true; return 0; }

    const double pref = alpha * alpha * 2.0 / Sc2;

    int n = 0;
    // Face-normal part: dS_n/dz_P = -1/h, dS_n/dz_nb = +1/h.
    out[n++] = { r,  c,  pref * s_n * (-1.0 / h) };
    out[n++] = { nr, nc, pref * s_n * (+1.0 / h) };

    // Tangential part, averaged across the face -> both cells contribute their
    // own centred difference, each weighted 1/2 * 1/(2h) = 1/(4h).
    const double wt = pref * s_t * (1.0 / (4.0 * h));

    int ar = 0, ac = 0;
    if (tangent_is_x)
    {
        st.dom.neighbour(r, c, 0, +1, ar, ac);   out[n++] = { ar, ac, +wt };
        st.dom.neighbour(r, c, 0, -1, ar, ac);   out[n++] = { ar, ac, -wt };
        st.dom.neighbour(nr, nc, 0, +1, ar, ac); out[n++] = { ar, ac, +wt };
        st.dom.neighbour(nr, nc, 0, -1, ar, ac); out[n++] = { ar, ac, -wt };
    }
    else
    {
        st.dom.neighbour(r, c, +1, 0, ar, ac);   out[n++] = { ar, ac, +wt };
        st.dom.neighbour(r, c, -1, 0, ar, ac);   out[n++] = { ar, ac, -wt };
        st.dom.neighbour(nr, nc, +1, 0, ar, ac); out[n++] = { ar, ac, +wt };
        st.dom.neighbour(nr, nc, -1, 0, ar, ac); out[n++] = { ar, ac, -wt };
    }
    return n;
}

} // namespace

#endif // LSDTT_HAS_EIGEN


// -----------------------------------------------------------------------------
// compute_hillslope_newton
// -----------------------------------------------------------------------------
lsdtt_hillslope_newton_result compute_hillslope_newton(
    const lsdtt_raster_cpu&               dem,
    const std::vector<std::uint8_t>&      is_channel,
    const lsdtt_hillslope_newton_options& opts,
    const std::vector<float>&             z_old,
    lsdtt_hillslope_newton_workspace*     ws)
{
    lsdtt_hillslope_newton_result result;

#ifndef LSDTT_HAS_EIGEN
    (void)is_channel; (void)opts; (void)z_old; (void)ws;
    result.elevation = dem;
    result.had_eigen = false;
    result.converged = false;
    result.message   = "compute_hillslope_newton requires Eigen; lsdtt3 was built "
                       "without it (-DLSDTT_EIGEN=OFF). No solve was performed.";
    return result;
#else

#ifdef _OPENMP
    if (opts.n_threads > 0) omp_set_num_threads(opts.n_threads);
#endif

    result.had_eigen = true;

    const lsdtt_raster_meta& meta = dem.get_meta();
    const int    nrows = meta.n_rows;
    const int    ncols = meta.n_cols;
    const double h     = static_cast<double>(meta.cell_size);
    const std::vector<std::uint8_t>& mask = dem.get_mask();

    const lsdtt_hillslope_domain dom =
        lsdtt_hillslope_domain::make(nrows, ncols, opts.boundary_conditions, mask);

    const double D      = static_cast<double>(opts.D);
    const double U      = static_cast<double>(opts.U);
    const double Sc     = static_cast<double>(opts.Sc);
    const double inv_dt = opts.inv_dt;
    const double DA     = D / (h * h);

    // Working field in double, as for the SOR solver: the residual is a small
    // quantity measured on a large field.
    const std::vector<float>& z_in = dem.get_data();
    std::vector<double> z(z_in.begin(), z_in.end());

    std::vector<double> zn;
    if (inv_dt > 0.0)
    {
        if (!z_old.empty()) zn.assign(z_old.begin(), z_old.end());
        else                zn = z;      // no z^n supplied: treat the input as z^n
    }

    // -------------------------------------------------------------------------
    // Degrees of freedom: everything valid that is not a channel and not on an
    // 'o' edge. Dirichlet cells keep their input value and feed the residual of
    // their neighbours.
    // -------------------------------------------------------------------------
    const std::size_t n_cells = static_cast<std::size_t>(nrows) * ncols;
    std::vector<int>  dof_of(n_cells, -1);
    std::vector<int>  rc_of_dof;
    rc_of_dof.reserve(n_cells * 2);

    int n_dof = 0;
    for (int r = 0; r < nrows; ++r)
        for (int c = 0; c < ncols; ++c)
        {
            const std::size_t s = dom.index(r, c);
            if (!mask[s]) continue;
            if (s < is_channel.size() && is_channel[s]) continue;
            if (dom.is_dirichlet_edge(r, c)) continue;
            dof_of[s] = n_dof++;
            rc_of_dof.push_back(r);
            rc_of_dof.push_back(c);
        }

    if (n_dof == 0)
    {
        result.elevation = dem;
        result.converged = true;
        result.message   = "no hillslope degrees of freedom (all cells are masked, "
                           "channels, or Dirichlet edges); field returned unchanged";
        return result;
    }

    // -------------------------------------------------------------------------
    // Residual F, and the infinity norm over the DOFs.
    // -------------------------------------------------------------------------
    auto residual = [&](const std::vector<double>& zz, std::vector<double>& F) -> double
    {
        lsdtt_hillslope_stencil<double> st(dom, zz.data(), h);
        double norm = 0.0;

        for (int d = 0; d < n_dof; ++d)
        {
            const int r = rc_of_dof[static_cast<std::size_t>(2 * d)];
            const int c = rc_of_dof[static_cast<std::size_t>(2 * d) + 1];

            double a[4], sn[4];
            st.face_alphas(r, c, Sc, opts.nonlinear, a, sn, opts.slope_ratio_cap);

            // sum_f alpha_f (z_nb - z_P) == h * sum_f alpha_f * S_n,f
            double flux = 0.0;
            for (int k = 0; k < 4; ++k) flux += a[k] * sn[k] * h;

            double f = -U - DA * flux;
            if (inv_dt > 0.0)
                f += inv_dt * (zz[dom.index(r, c)] - zn[dom.index(r, c)]);

            F[static_cast<std::size_t>(d)] = f;
            norm = std::max(norm, std::fabs(f));
        }
        return norm;
    };

    // Largest |grad z|/Sc over the DOFs -- both a diagnostic and the line-search
    // feasibility test.
    auto max_ratio = [&](const std::vector<double>& zz) -> double
    {
        lsdtt_hillslope_stencil<double> st(dom, zz.data(), h);
        double m = 0.0;
        for (int d = 0; d < n_dof; ++d)
        {
            const int r = rc_of_dof[static_cast<std::size_t>(2 * d)];
            const int c = rc_of_dof[static_cast<std::size_t>(2 * d) + 1];
            const double sx = st.slope_x(r, c);
            const double sy = st.slope_y(r, c);
            m = std::max(m, std::sqrt(sx * sx + sy * sy) / Sc);
        }
        return m;
    };

    // -------------------------------------------------------------------------
    // Cache validity.
    //
    // The key must cover everything that changes the MATRIX. The DOF hash folds in
    // the mask, the channel set and the Dirichlet edges in one pass, which is what
    // actually determines the structure -- cheaper than comparing arrays and it
    // cannot miss a rearrangement that leaves n_dof the same.
    std::uint64_t dof_hash = 1469598103934665603ull;
    for (int d = 0; d < n_dof; ++d)
    {
        const std::uint64_t v =
            static_cast<std::uint64_t>(rc_of_dof[static_cast<std::size_t>(2 * d)]) * 1000003ull
          + static_cast<std::uint64_t>(rc_of_dof[static_cast<std::size_t>(2 * d) + 1]);
        dof_hash = (dof_hash ^ v) * 1099511628211ull;
    }

    lsdtt_hillslope_newton_workspace::impl* wsp = ws ? ws->p_.get() : nullptr;

    bool cache_ok = false;
    if (wsp)
    {
        cache_ok = wsp->have_key
                && wsp->nrows == nrows && wsp->ncols == ncols
                && wsp->bc[0] == dom.bc_n && wsp->bc[1] == dom.bc_e
                && wsp->bc[2] == dom.bc_s && wsp->bc[3] == dom.bc_w
                && wsp->D == D && wsp->Sc == Sc && wsp->inv_dt == inv_dt
                && wsp->h == h
                && wsp->nonlinear == opts.nonlinear
                && wsp->n_dof == static_cast<long int>(n_dof)
                && wsp->dof_hash == dof_hash;

        if (!cache_ok)
        {
            wsp->invalidate();
            wsp->have_key = true;
            wsp->nrows = nrows; wsp->ncols = ncols;
            wsp->bc[0] = dom.bc_n; wsp->bc[1] = dom.bc_e;
            wsp->bc[2] = dom.bc_s; wsp->bc[3] = dom.bc_w;
            wsp->D = D; wsp->Sc = Sc; wsp->inv_dt = inv_dt; wsp->h = h;
            wsp->nonlinear = opts.nonlinear;
            wsp->n_dof = static_cast<long int>(n_dof);
            wsp->dof_hash = dof_hash;
        }
    }

    // A LINEAR problem has a constant matrix: alpha is 1 regardless of z, so once
    // the key matches, both the pattern and the numeric factorisation stay valid
    // across any number of calls. Only the right-hand side moves.
    bool reuse_numeric =
        wsp && cache_ok && !opts.nonlinear && wsp->numeric_valid;

    std::vector<double> F(static_cast<std::size_t>(n_dof), 0.0);
    double fnorm = residual(z, F);

    // The matrix and its factorisation live in the workspace when one was given,
    // so they outlive this call.
    Eigen::SparseMatrix<double>                  J_local(n_dof, n_dof);
    Eigen::SparseLU<Eigen::SparseMatrix<double>> lu_local;

    Eigen::SparseMatrix<double>&                  J  = wsp ? wsp->J  : J_local;
    Eigen::SparseLU<Eigen::SparseMatrix<double>>& lu = wsp ? wsp->lu : lu_local;

    // A workspace-held matrix is default-constructed 0x0 on the first call, and
    // setFromTriplets does not resize -- it would write outside the matrix. Size it
    // here rather than at construction, because n_dof is not known until now.
    if (J.rows() != n_dof || J.cols() != n_dof)
    {
        J.resize(n_dof, n_dof);
        if (wsp) wsp->invalidate();      // any cached factorisation is meaningless
        reuse_numeric = false;
    }

    std::vector<Eigen::Triplet<double>> trips;
    if (!reuse_numeric) trips.reserve(static_cast<std::size_t>(n_dof) * 12);

    const double feasible_ratio = std::sqrt(opts.slope_ratio_cap);
    bool capped_anywhere = false;

    long int pattern_analyses_local = 0; (void)pattern_analyses_local;

    int it = 0;
    for (; it < opts.max_iterations; ++it)
    {
        if (fnorm <= opts.residual_tol) { result.converged = true; break; }

        // ---------------------------------------------------------------------
        // Assemble J = dF/dz -- unless a valid constant matrix is already cached.
        // ---------------------------------------------------------------------
        if (!reuse_numeric)
        {
        trips.clear();
        lsdtt_hillslope_stencil<double> st(dom, z.data(), h);

        static constexpr int k_dr[4] = { -1,  0, +1,  0 };
        static constexpr int k_dc[4] = {  0, +1,  0, -1 };

        for (int d = 0; d < n_dof; ++d)
        {
            const int r = rc_of_dof[static_cast<std::size_t>(2 * d)];
            const int c = rc_of_dof[static_cast<std::size_t>(2 * d) + 1];

            double a[4], sn[4];
            st.face_alphas(r, c, Sc, opts.nonlinear, a, sn, opts.slope_ratio_cap);

            auto add = [&](int rr, int cc, double w)
            {
                if (w == 0.0) return;
                const int dd = dof_of[dom.index(rr, cc)];
                if (dd < 0) return;                    // Dirichlet: no unknown
                trips.emplace_back(d, dd, w);
            };

            if (inv_dt > 0.0) add(r, c, inv_dt);

            for (int dir = 0; dir < 4; ++dir)
            {
                int nr = 0, nc = 0;
                dom.neighbour(r, c, k_dr[dir], k_dc[dir], nr, nc);

                // -(D/h^2) * d/dz [ alpha_f * (z_nb - z_P) ]
                //   = -(D/h^2) * [ alpha_f * (d z_nb - d z_P)
                //                + (z_nb - z_P) * d alpha_f ]
                add(nr, nc, -DA * a[dir]);
                add(r,  c,  +DA * a[dir]);

                if (opts.nonlinear)
                {
                    jac_entry ent[8];
                    bool capped = false;
                    const int n = dalpha_dz(st, r, c, dir, Sc, opts.slope_ratio_cap,
                                            a[dir], sn[dir], ent, capped);
                    if (capped) capped_anywhere = true;

                    const double dz_face = sn[dir] * h;      // (z_nb - z_P)
                    for (int k = 0; k < n; ++k)
                        add(ent[k].r, ent[k].c, -DA * dz_face * ent[k].w);
                }
            }
        }

        // Duplicate (row,col) pairs are summed by setFromTriplets, which is what
        // a periodic wrap or a repeated tangential contributor needs.
        J.setFromTriplets(trips.begin(), trips.end());
        J.makeCompressed();
        }   // end assembly

        Eigen::VectorXd rhs(n_dof);
        for (int d = 0; d < n_dof; ++d) rhs[d] = -F[static_cast<std::size_t>(d)];

        Eigen::VectorXd step(n_dof);
        bool solve_ok = true;

        // Resolve AUTO once we know the problem size.
        //
        // A CACHED direct factorisation changes the arithmetic completely. The
        // threshold exists because fill-in makes a one-off factorisation
        // expensive at large N -- but when a workspace is caching it for a linear
        // problem, that cost is paid ONCE and then amortised over every timestep,
        // while an iterative solve pays its full price on every call. So the
        // direct path stays worthwhile far past the one-off threshold, and AUTO
        // raises it (4x) when the cache can actually be used.
        const bool cache_will_amortise = (wsp != nullptr) && !opts.nonlinear;
        const long int direct_limit =
            opts.direct_solve_max_unknowns * (cache_will_amortise ? 4 : 1);

        const bool use_iterative =
            (opts.solver_choice == 1) ||
            (opts.solver_choice < 0 &&
             static_cast<long int>(n_dof) > direct_limit);

        if (use_iterative)
        {
            Eigen::BiCGSTAB<Eigen::SparseMatrix<double>,
                            Eigen::IncompleteLUT<double>> solver;
            solver.compute(J);
            if (solver.info() != Eigen::Success) solve_ok = false;
            else
            {
                step = solver.solve(rhs);
                if (solver.info() != Eigen::Success) solve_ok = false;
            }
            if (wsp) ++wsp->n_solves;
        }
        else
        {
            // The sparsity PATTERN is fixed for a given mask / BC set, so the
            // symbolic analysis is hoisted whenever the key matches -- true for the
            // nonlinear case too, where only the values move. The NUMERIC
            // factorisation is additionally reused when the matrix itself cannot
            // have changed (linear, fixed dt).
            if (!reuse_numeric)
            {
                if (!wsp || !wsp->pattern_valid)
                {
                    lu.analyzePattern(J);
                    if (wsp) { wsp->pattern_valid = true; ++wsp->n_pattern_analyses; }
                    else     ++pattern_analyses_local;
                }
                lu.factorize(J);
                if (wsp) ++wsp->n_factorisations;
                if (lu.info() != Eigen::Success) solve_ok = false;
                else if (wsp && !opts.nonlinear) wsp->numeric_valid = true;
            }

            if (solve_ok)
            {
                step = lu.solve(rhs);
                if (lu.info() != Eigen::Success) solve_ok = false;
                if (wsp) ++wsp->n_solves;
            }
        }

        if (!solve_ok)
        {
            result.message = "linear solve failed at Newton iteration "
                           + std::to_string(it)
                           + " (the Jacobian is 9-point and non-symmetric; try "
                             "solver_choice = 1 for BiCGSTAB+ILUT)";
            break;
        }

        // ---------------------------------------------------------------------
        // Line search. Newton without a safeguard steps through Sc; accept the
        // largest fraction that BOTH reduces the residual norm and leaves every
        // face sub-critical.
        // ---------------------------------------------------------------------
        double lambda = 1.0;
        std::vector<double> z_try(z.size());
        std::vector<double> F_try(static_cast<std::size_t>(n_dof), 0.0);
        double fnorm_try = fnorm;
        double max_dz    = 0.0;
        bool   accepted  = false;

        const int max_ls = opts.enable_line_search ? opts.max_line_search : 1;

        // Is the CURRENT iterate already super-critical?
        //
        // This matters, and getting it wrong makes the solver useless on real
        // input. A strict "every trial step must be sub-critical" filter cannot
        // be satisfied when the starting field is itself super-critical: shrinking
        // lambda drives the trial toward that same infeasible field, so every
        // backtrack fails and Newton gives up at iteration 0 without moving. Users
        // supply arbitrary DEMs, and an initial guess above Sc is ordinary.
        //
        // So the filter is conditional. From a feasible iterate, refuse to leave
        // the feasible region (the usual guard). From an infeasible one, run a
        // RESTORATION step instead: accept whatever reduces the slope ratio, and
        // let the residual test resume once the field is back under Sc.
        const double ratio_now = opts.nonlinear ? max_ratio(z) : 0.0;
        const bool   restoring = opts.nonlinear && (ratio_now > feasible_ratio);

        for (int ls = 0; ls < max_ls; ++ls)
        {
            z_try = z;
            max_dz = 0.0;
            for (int d = 0; d < n_dof; ++d)
            {
                const int r = rc_of_dof[static_cast<std::size_t>(2 * d)];
                const int c = rc_of_dof[static_cast<std::size_t>(2 * d) + 1];
                const double dz = lambda * step[d];
                z_try[dom.index(r, c)] += dz;
                max_dz = std::max(max_dz, std::fabs(dz));
            }

            if (!opts.enable_line_search) { accepted = true; break; }

            if (opts.nonlinear)
            {
                const double ratio_try = max_ratio(z_try);

                if (restoring)
                {
                    // Progress toward feasibility is the acceptance test.
                    if (ratio_try >= ratio_now)
                    {
                        lambda *= 0.5;
                        ++result.n_line_search_steps;
                        continue;
                    }
                    fnorm_try = residual(z_try, F_try);
                    accepted = true;
                    break;
                }

                if (ratio_try > feasible_ratio)
                {
                    lambda *= 0.5;
                    ++result.n_line_search_steps;
                    continue;
                }
            }

            fnorm_try = residual(z_try, F_try);
            if (fnorm_try < fnorm) { accepted = true; break; }

            lambda *= 0.5;
            ++result.n_line_search_steps;
        }

        if (!accepted)
        {
            result.message = "line search failed to reduce the residual at Newton "
                             "iteration " + std::to_string(it)
                           + "; the field is likely at or beyond Sc. Fall back to "
                             "damped Picard, or reduce dt.";
            break;
        }

        z    = z_try;
        F    = F_try;
        fnorm = opts.enable_line_search ? fnorm_try : residual(z, F);
        result.final_max_dz = max_dz;

        if (fnorm <= opts.residual_tol) { result.converged = true; ++it; break; }
    }

    result.n_iterations   = it;
    result.final_residual = fnorm;
    result.max_slope_ratio = max_ratio(z);

    // Whether the cap binds in the CONVERGED FIELD -- not whether it ever bound
    // during the iteration.
    //
    // The distinction matters because of what this flag is for: it tells the user
    // that near-Sc behaviour in the answer is set by the limiter rather than by
    // the flux law. A solve that STARTS super-critical clamps freely during the
    // restoration phase and then converges to a field nowhere near the cap;
    // reporting that as "limiter active" is a false alarm, and a false alarm on a
    // correctness flag is worse than no flag. `capped_anywhere` (the
    // during-iteration signal) is kept for the message below but does not set the
    // flag on its own.
    result.limiter_active = (result.max_slope_ratio >= std::sqrt(opts.slope_ratio_cap));

    if (capped_anywhere && !result.limiter_active && result.converged)
        result.message = "converged (the slope cap bound during iteration, but not "
                         "in the converged field)";

    if (result.converged && result.message.empty())
        result.message = "converged";
    else if (!result.converged && result.message.empty())
        result.message = "maximum Newton iterations reached; ||F||inf="
                       + std::to_string(result.final_residual) + " m/yr";

    std::vector<float> z_out(z.begin(), z.end());
    result.elevation = lsdtt_raster_cpu(std::move(z_out),
                                        std::vector<std::uint8_t>(mask), meta);
    return result;
#endif // LSDTT_HAS_EIGEN
}


// -----------------------------------------------------------------------------
// compute_hillslope_newton_transient
//
// Adaptive backward-Euler sub-stepping around the Newton core.
//
// Note what is NOT here: a CFL condition. Backward Euler is unconditionally
// stable, so sub-step size is chosen for accuracy and for nonlinear solvability,
// never for stability. The controller keys on the Newton iteration count -- a step
// that converged in a couple of iterations could have been larger; one that needed
// many was straining -- with an optional max|dz| gate for trajectory fidelity,
// because a stable step can still be an inaccurate one.
// -----------------------------------------------------------------------------
lsdtt_hillslope_newton_transient_result compute_hillslope_newton_transient(
    const lsdtt_raster_cpu&                         dem,
    const std::vector<std::uint8_t>&                is_channel,
    float                                           dt_total,
    const lsdtt_hillslope_newton_transient_options& opts,
    lsdtt_hillslope_newton_workspace*               ws)
{
    lsdtt_hillslope_newton_transient_result result;
    result.elevation = dem;

    if (!(dt_total > 0.0f))
    {
        result.converged = true;
        result.had_eigen = true;
        result.message   = "dt_total <= 0: nothing to integrate";
        return result;
    }

    const double total = static_cast<double>(dt_total);

    double dt = (opts.dt_init > 0.0f) ? static_cast<double>(opts.dt_init) : total;
    const double dt_min = (opts.dt_min > 0.0f) ? static_cast<double>(opts.dt_min)
                                               : total * 1e-9;
    const double dt_max = (opts.dt_max > 0.0f) ? static_cast<double>(opts.dt_max)
                                               : total;
    dt = std::min(dt, dt_max);

    // z^n for the first sub-step is the input field. Dirichlet cells (channels and
    // 'o' edges) keep their dem value throughout -- already at t^{n+1} by contract.
    std::vector<float> zn = dem.get_data();
    lsdtt_raster_cpu   current = dem;

    double remaining = total;
    result.converged = true;

    while (remaining > 0.0)
    {
        const double step = std::min(dt, remaining);

        lsdtt_hillslope_newton_options nopts = opts.newton;
        nopts.inv_dt = 1.0 / step;

        const lsdtt_hillslope_newton_result r =
            compute_hillslope_newton(current, is_channel, nopts, zn, ws);

        result.had_eigen = r.had_eigen;
        if (!r.had_eigen)
        {
            result.converged = false;
            result.message   = r.message;
            return result;
        }

        result.total_newton_iters += r.n_iterations;

        const bool dz_ok = (opts.dz_target <= 0.0f) ||
                           (r.final_max_dz <= static_cast<double>(opts.dz_target));

        if (!r.converged || !dz_ok)
        {
            // Reject and retry smaller. Backward Euler cannot be unstable here, so
            // a rejection means the NONLINEAR solve was too hard at this size, or
            // the step moved the field further than the accuracy gate allows.
            ++result.n_rejected_substeps;

            if (result.n_rejected_substeps > opts.max_rejects || step <= dt_min)
            {
                result.converged = false;
                result.message   = r.converged
                    ? "sub-step rejected by the accuracy gate at the minimum step size"
                    : ("Newton failed at the minimum sub-step: " + r.message);
                break;
            }

            dt = std::max(dt_min, step * static_cast<double>(opts.shrink_factor));
            continue;
        }

        // Accepted.
        current   = r.elevation;
        zn        = current.get_data();
        remaining -= step;

        ++result.n_substeps;
        result.final_dt        = step;
        result.max_slope_ratio = r.max_slope_ratio;
        if (r.limiter_active) result.limiter_active = true;

        if (r.n_iterations <= opts.iters_lo)
            dt = std::min(dt_max, step * static_cast<double>(opts.grow_factor));
        else if (r.n_iterations >= opts.iters_hi)
            dt = std::max(dt_min, step * static_cast<double>(opts.shrink_factor));
    }

    result.elevation = current;
    if (result.converged && result.message.empty())
        result.message = "converged";
    return result;
}


// -----------------------------------------------------------------------------
// compute_hillslope_steady_state_auto
// -----------------------------------------------------------------------------
lsdtt_hillslope_steady_state_result compute_hillslope_steady_state_auto(
    const lsdtt_raster_cpu&                     dem,
    const std::vector<std::uint8_t>&            is_channel,
    const lsdtt_hillslope_steady_state_options& opts,
    bool                                        force_sor,
    lsdtt_hillslope_solve_report*               report,
    lsdtt_hillslope_newton_workspace*           ws)
{
    if (!force_sor)
    {
        lsdtt_hillslope_newton_workspace local_ws;
        lsdtt_hillslope_newton_workspace* active_ws = ws ? ws : &local_ws;
        lsdtt_hillslope_newton_options nopts;
        nopts.D  = opts.D;
        nopts.U  = opts.U;
        nopts.Sc = opts.Sc;
        nopts.nonlinear           = opts.nonlinear;
        nopts.boundary_conditions = opts.boundary_conditions;
        nopts.inv_dt              = 0.0;              // the steady problem
        nopts.n_threads           = opts.n_threads;
        nopts.max_iterations      = std::min(opts.max_iterations, 50);

        lsdtt_hillslope_newton_result nr =
            compute_hillslope_newton(dem, is_channel, nopts,
                                     std::vector<float>(), active_ws);

        int continuation_iterations = 0;
        if (!nr.converged && nr.had_eigen && opts.nonlinear)
        {
            // A nearly flat initial field can have a steady Newton direction far
            // outside the sub-critical region. Backtracking that one direction
            // makes negligible progress. Pseudo-time continuation supplies the
            // missing globalization: easy backward-Euler solves approach the
            // steady basin before one final exact solve.
            lsdtt_raster_cpu current = dem;
            const double h = dem.get_meta().cell_size;
            double dt = h * h / std::max(static_cast<double>(opts.D), 1e-12);
            bool continuation_ok = true;

            for (int stage = 0; stage < 12; ++stage)
            {
                nopts.inv_dt = 1.0 / dt;
                nopts.max_iterations = 50;
                const std::vector<float> z_old = current.get_data();
                const lsdtt_hillslope_newton_result step =
                    compute_hillslope_newton(current, is_channel, nopts,
                                             z_old, active_ws);
                continuation_iterations += step.n_iterations;
                if (!step.converged)
                {
                    continuation_ok = false;
                    break;
                }
                current = step.elevation;
                dt *= 4.0;
            }

            if (continuation_ok)
            {
                nopts.inv_dt = 0.0;
                nopts.max_iterations = opts.max_iterations;
                nr = compute_hillslope_newton(current, is_channel, nopts,
                                              std::vector<float>(), active_ws);
                nr.n_iterations += continuation_iterations;
                if (nr.converged)
                    nr.message = "converged after pseudo-time continuation";
            }
        }

        if (nr.converged)
        {
            if (report)
            {
                report->solver_used = nr.had_eigen ? "Newton" : "Newton (no Eigen)";
                report->detail      = "iterations=" + std::to_string(nr.n_iterations)
                                    + "  ||F||inf=" + std::to_string(nr.final_residual)
                                    + " m/yr  " + nr.message;
                if (nr.limiter_active)
                    report->warning =
                        "the sub-critical limiter binds in the converged field; "
                        "near-Sc behaviour there is set by the cap, not the flux law";
            }

            lsdtt_hillslope_steady_state_result out;
            out.elevation      = nr.elevation;
            out.n_iterations   = nr.n_iterations;
            out.final_residual = static_cast<float>(nr.final_max_dz);
            out.converged      = true;
            return out;
        }

        if (report)
            report->warning = "the Newton hillslope solve did not converge ("
                            + nr.message + "); falling back to SOR with Picard "
                              "linearisation";
    }

    const lsdtt_hillslope_steady_state_result sr =
        compute_hillslope_steady_state(dem, is_channel, opts);

    if (report)
    {
        report->solver_used = "SOR";
        report->detail      = "sweeps=" + std::to_string(sr.n_iterations)
                            + "  max|dz|=" + std::to_string(sr.final_residual) + " m";
    }
    return sr;
}

} // namespace lsdtt
