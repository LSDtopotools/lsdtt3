// =============================================================================
// lsdtt_hillslope_operator_cpu.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
//
// -----------------------------------------------------------------------------
// THE shared discrete hillslope operator.
//
// Every hillslope solver in lsdtt3 -- the SOR steady state, the explicit
// transient solvers, the implicit transient solver, and the Newton scheme --
// must discretise the same equation:
//
//     dz/dt = U + (1/h^2) sum_f  D_eff,f (z_nb,f - z_P)
//
// with the Roering (1999) face diffusivity
//
//     D_eff,f = D / (1 - (S_n,f^2 + S_t,f^2)/Sc^2)
//
// where S_n is the slope NORMAL to the face and S_t the slope TANGENTIAL to it.
// Steady state sets the left-hand side to zero. A transient run at dt -> infinity
// must therefore land on the SOR steady state, and that agreement is one of the
// cheapest validations available -- but ONLY if the two solvers really do share
// one operator.
//
// WHY THIS HEADER EXISTS.  They did not. Three separate copies of the neighbour
// accessor grew up in lsdtt_hillslope_steady_state_cpu.cpp,
// lsdtt_hillslope_transport_cpu.cpp and lsdtt_hillslope_implicit_cpu.cpp, and
// they drifted -- in three different ways, each of which was a real defect:
//
//   1. FACE SYMMETRY.  The SOR solver evaluated the tangential slope at the cell
//      being updated, so the N face of (r,c) and the S face of (r-1,c) -- one
//      physical face -- received DIFFERENT diffusivities. The operator was
//      non-symmetric and the scheme non-conservative: flux leaving one cell
//      through a face did not equal flux entering its neighbour through it. The
//      nonlinear solve limit-cycled and slopes over-steepened.
//   2. THE MEANING OF 'o'.  The implicit solver treated 'o' as Dirichlet (base
//      level); SOR and ADI treated it as Neumann, identical to 'n'. So they were
//      solving different boundary-value problems from the same input.
//   3. PADDING A SHORT BC STRING.  The implicit solver padded with 'o', the
//      others with 'n' -- the same input meant different things in each.
//
// All three are now settled here, in one place, and the contract is:
//
//   'p'  periodic: wrap to the opposite edge.
//   'o'  DIRICHLET at base level: the cell is held at its input elevation and is
//        NOT a degree of freedom. 'o' is never a Neumann condition.
//   'n'  Neumann zero-flux: reflect the donor cell back across the boundary.
//
// A short or empty string pads with 'o'. That is the safer default: an all-'n'
// domain with no channels has no Dirichlet data anywhere, so the problem is
// ill posed and the field grows without bound.
//
// The functions here are header-inline and take a plain field pointer rather
// than a raster, so a solver can apply them to whatever working array it holds
// (float or double) without copying.
// =============================================================================

#ifndef LSDTT_HILLSLOPE_OPERATOR_CPU_HPP
#define LSDTT_HILLSLOPE_OPERATOR_CPU_HPP

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace lsdtt
{

// -----------------------------------------------------------------------------
// lsdtt_hillslope_domain
//
// The grid, its mask, and its boundary codes: everything the operator needs that
// does not depend on the elevation field. Cheap to copy; build one per solve.
// -----------------------------------------------------------------------------
struct lsdtt_hillslope_domain
{
    int n_rows = 0;
    int n_cols = 0;

    // NESW, as everywhere else in lsdtt3.
    char bc_n = 'o';
    char bc_e = 'o';
    char bc_s = 'o';
    char bc_w = 'o';

    const std::vector<std::uint8_t>* mask = nullptr;

    // -------------------------------------------------------------------------
    // Build from a boundary-condition string, padding short strings with 'o'.
    // -------------------------------------------------------------------------
    static lsdtt_hillslope_domain make(int nrows, int ncols,
                                       const std::string& bc,
                                       const std::vector<std::uint8_t>& msk)
    {
        auto code = [&](std::size_t i) -> char
        {
            return (bc.size() > i) ? bc[i] : 'o';
        };

        lsdtt_hillslope_domain d;
        d.n_rows = nrows;
        d.n_cols = ncols;
        d.bc_n = code(0); d.bc_e = code(1); d.bc_s = code(2); d.bc_w = code(3);
        d.mask = &msk;
        return d;
    }

    bool periodic_ns() const { return bc_n == 'p' || bc_s == 'p'; }
    bool periodic_ew() const { return bc_e == 'p' || bc_w == 'p'; }

    std::size_t index(int r, int c) const
    {
        return static_cast<std::size_t>(r) * static_cast<std::size_t>(n_cols)
             + static_cast<std::size_t>(c);
    }

    bool valid(int r, int c) const
    {
        return (*mask)[index(r, c)] != 0;
    }

    // -------------------------------------------------------------------------
    // is_dirichlet_edge: true when (r,c) lies on an edge declared 'o'.
    //
    // 'o' is base level. Such a cell keeps its input elevation and contributes to
    // its neighbours' updates exactly as a channel cell does. It is NOT solved for.
    // -------------------------------------------------------------------------
    bool is_dirichlet_edge(int r, int c) const
    {
        return (r == 0          && bc_n == 'o')
            || (r == n_rows - 1 && bc_s == 'o')
            || (c == 0          && bc_w == 'o')
            || (c == n_cols - 1 && bc_e == 'o');
    }

    // -------------------------------------------------------------------------
    // neighbour: resolve (r+dr, c+dc) to a real cell.
    //
    // Periodic edges wrap. Every other edge, and every masked interior
    // neighbour, resolves to the CELL ITSELF -- which is what makes the face a
    // genuine no-flux boundary, because the face-normal slope
    // (z_ghost - z_P)/h is then exactly zero.
    //
    // THIS IS NOT WHAT THE SOR SOLVER USED TO DO, and the difference is a real
    // defect rather than a stylistic choice. It resolved a no-flux neighbour as
    // `max(0, r - dr)` -- the OPPOSITE neighbour, mirroring about the cell
    // centre. That gives a zero CENTRED slope at the boundary cell, which looks
    // right, but a NON-ZERO face flux: row 0 then receives flux from row 1
    // twice, once through its real S face and once through its fictitious N
    // face. Mass is not conserved on a closed domain, which is exactly what
    // `unit_hillslope_transport_bc` T4 measures.
    //
    // The other three solvers already had it right in their own ways -- the
    // implicit solver skips a no-flux face outright, ADI zeroes its conductance,
    // and the explicit transport solver clamped (which returns the cell itself).
    // The mirror was the outlier, and this is the rule they now all share.
    //
    // Returning the INDEX rather than the value is what makes face-symmetric
    // coefficients possible: the caller can evaluate the neighbour's own centred
    // slope with identical boundary handling.
    // -------------------------------------------------------------------------
    void neighbour(int r, int c, int dr, int dc, int& nr, int& nc) const
    {
        nr = r + dr;
        nc = c + dc;

        if (nr < 0)             nr = (bc_n == 'p') ? n_rows - 1 : r;
        else if (nr >= n_rows)  nr = (bc_s == 'p') ? 0          : r;

        if (nc < 0)             nc = (bc_w == 'p') ? n_cols - 1 : c;
        else if (nc >= n_cols)  nc = (bc_e == 'p') ? 0          : c;

        // Masked interior neighbour: the same no-flux treatment.
        if (!valid(nr, nc)) { nr = r; nc = c; }
    }
};

// -----------------------------------------------------------------------------
// lsdtt_hillslope_stencil
//
// Reads a field through the domain's boundary rules. Templated on the field's
// scalar type so a solver iterating in double and one iterating in float share
// the same logic; all arithmetic is done in double regardless.
// -----------------------------------------------------------------------------
template <typename T>
struct lsdtt_hillslope_stencil
{
    const lsdtt_hillslope_domain& dom;
    const T*                      z;
    double                        h;      // cell size

    lsdtt_hillslope_stencil(const lsdtt_hillslope_domain& d, const T* field, double cell_size)
        : dom(d), z(field), h(cell_size) {}

    double at(int r, int c, int dr, int dc) const
    {
        int nr = 0, nc = 0;
        dom.neighbour(r, c, dr, dc, nr, nc);
        return static_cast<double>(z[dom.index(nr, nc)]);
    }

    double at(int r, int c) const
    {
        return static_cast<double>(z[dom.index(r, c)]);
    }

    // Centred slopes AT AN ARBITRARY CELL. Needed at a cell other than the one
    // being updated, so that a face coefficient can average across the face.
    double slope_x(int r, int c) const
    {
        return (at(r, c, 0, +1) - at(r, c, 0, -1)) / (2.0 * h);
    }
    double slope_y(int r, int c) const
    {
        return (at(r, c, +1, 0) - at(r, c, -1, 0)) / (2.0 * h);
    }

    // -------------------------------------------------------------------------
    // face_alpha: the dimensionless diffusivity multiplier on ONE face.
    //
    // dir: 0 = N, 1 = E, 2 = S, 3 = W.
    //
    // Returns 1 in linear mode. In nonlinear mode returns
    //
    //     alpha = 1 / (1 - (S_n^2 + S_t^2)/Sc^2)
    //
    // with S_t averaged ACROSS the face from the two cells sharing it. That
    // averaging is the whole point: it makes alpha a property of the FACE, so
    // the N face of (r,c) and the S face of its northern neighbour necessarily
    // agree, and the scheme is conservative. Evaluating S_t at the updating cell
    // alone -- which is what the SOR solver did until 2026-08-26 -- gives the two
    // sides different values and breaks conservation.
    //
    // s2 is clamped strictly below 1 (`ratio_cap`, default 0.9999) so the
    // singularity at S = Sc becomes a strong-but-finite restoring diffusivity
    // rather than a division by zero. A run in which the cap is active has its
    // near-Sc behaviour set by the limiter rather than by the flux law, so
    // callers that care about that regime should report when it fires.
    // -------------------------------------------------------------------------
    double face_alpha(int r, int c, int dir, double Sc,
                      bool nonlinear, double ratio_cap = 0.9999) const
    {
        if (!nonlinear) return 1.0;

        int dr = 0, dc = 0;
        bool tangent_is_x = false;
        switch (dir)
        {
            case 0: dr = -1; tangent_is_x = true;  break;   // N face: E-W tangent
            case 1: dc = +1; tangent_is_x = false; break;   // E face: N-S tangent
            case 2: dr = +1; tangent_is_x = true;  break;   // S face: E-W tangent
            default: dc = -1; tangent_is_x = false; break;  // W face: N-S tangent
        }

        int nr = 0, nc = 0;
        dom.neighbour(r, c, dr, dc, nr, nc);

        const double s_n = (at(nr, nc) - at(r, c)) / h;

        const double s_t = tangent_is_x
            ? 0.5 * (slope_x(r, c) + slope_x(nr, nc))
            : 0.5 * (slope_y(r, c) + slope_y(nr, nc));

        const double s2 = std::min((s_n * s_n + s_t * s_t) / (Sc * Sc), ratio_cap);
        return 1.0 / (1.0 - s2);
    }

    // -------------------------------------------------------------------------
    // face_alphas: all four face coefficients at once, in N,E,S,W order.
    //
    // Prefer this to four face_alpha() calls. The centre cell's two centred
    // slopes are shared by all four faces -- slope_x by N and S, slope_y by E
    // and W -- so computing them once saves four redundant centred-difference
    // evaluations (eight neighbour resolutions) per cell. Every solver needs all
    // four coefficients together, so this is the natural interface; face_alpha()
    // remains for the rare single-face case and for tests.
    //
    // Optionally returns the face-NORMAL slopes too (out_sn, N,E,S,W order):
    // Newton needs them for the Jacobian, and a limiter needs them to detect a
    // face crossing Sc, so recomputing them outside would be wasteful.
    // -------------------------------------------------------------------------
    void face_alphas(int r, int c, double Sc, bool nonlinear,
                     double out_alpha[4], double* out_sn = nullptr,
                     double ratio_cap = 0.9999) const
    {
        static constexpr int k_dr[4] = { -1,  0, +1,  0 };
        static constexpr int k_dc[4] = {  0, +1,  0, -1 };

        const double zc = at(r, c);

        // N and S share the E-W tangent; E and W share the N-S tangent.
        const double sxc = nonlinear ? slope_x(r, c) : 0.0;
        const double syc = nonlinear ? slope_y(r, c) : 0.0;

        for (int dir = 0; dir < 4; ++dir)
        {
            int nr = 0, nc = 0;
            dom.neighbour(r, c, k_dr[dir], k_dc[dir], nr, nc);

            const double s_n = (at(nr, nc) - zc) / h;
            if (out_sn) out_sn[dir] = s_n;

            if (!nonlinear) { out_alpha[dir] = 1.0; continue; }

            const bool tangent_is_x = (dir == 0 || dir == 2);
            const double s_t = tangent_is_x
                ? 0.5 * (sxc + slope_x(nr, nc))
                : 0.5 * (syc + slope_y(nr, nc));

            const double s2 = std::min((s_n * s_n + s_t * s_t) / (Sc * Sc), ratio_cap);
            out_alpha[dir] = 1.0 / (1.0 - s2);
        }
    }
};

} // namespace lsdtt

#endif // LSDTT_HILLSLOPE_OPERATOR_CPU_HPP
