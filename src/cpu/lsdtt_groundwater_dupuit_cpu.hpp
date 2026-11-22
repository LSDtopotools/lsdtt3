// lsdtt_groundwater_dupuit_cpu.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// Dupuit - Boussinesq groundwater in an unconfined aquifer on a sloping,
// impermeable base, with seepage to the surface.
//
// After the Landlab GroundwaterDupuitPercolator: Litwin, D. G., Tucker, G. E.,
// Barnhart, K. R. & Harman, C. J. (2020), GroundwaterDupuitPercolator: a
// Landlab component for groundwater flow, Journal of Open Source Software
// 5(46), 1935, https://doi.org/10.21105/joss.01935 (MIT licence); seepage after
// Marcais, de Dreuzy & Erhel (2017), Adv. Water Resour. 109; sloping base after
// Childs (1971), Water Resour. Res. 7. Written from the published equations
// and the component's documented algorithm; not a translation.
//
// Named for the method, not for slope stability: lsdtt-slope-stability uses it
// for pore pressure, and distributed hydrology may later.
//
// EQUATIONS (h saturated thickness above the base z_b; d = z - z_b the
// regolith thickness; z_w = z_b + h the water table)
//
//   q   = -T(h) cos^2(alpha) grad z_w        per unit width, on 4-connected links
//   q_s = G(h/d) R(f - div q),  G(u) = exp(-(1 - u)/eps),  R(u) = max(u, 0)
//   n dh/dt = f - q_s - div q
//
//   alpha  slope of the base along the link; f recharge; n drainable porosity
//   T(h) = K h                                       uniform conductivity (Landlab)
//   T(h) = K0 z_m (exp(-(d - h)/z_m) - exp(-d/z_m))  K = K0 exp(-depth / z_m)
//
// The exponential profile is TOPMODEL's (deficit D = n (d - h), m = n z_m,
// T0 = K0 z_m as d -> infinity): with it, and with the water table parallel to
// the surface, the steady state is TOPMODEL's index relation, except that the
// slope term is sin(beta) cos(beta) rather than tan(beta).
//
// NUMERICS (as Landlab, except where noted)
//   * The link thickness is UPWIND (the cell with the higher water table), times
//     cos(alpha); the link gradient is times cos(alpha), giving cos^2.
//   * Thickness update: with div q held fixed over a substep, the regularised
//     ODE n dh/dt = (f - div q)(1 - G(h/d)) has a closed form (Landlab's
//     _update_thickness); where f <= div q it is linear. h is clipped at 0.
//   * Substeps: dt = min(C_c dx / |v|, C_vn n dx^2 / (4 T cos^2), remaining),
//     v = |q| / (n h_upwind); C_c = 0.5, C_vn = 0.8 by default.
//   * Steady state (run_to_steady): the march starts with a substep cap of the
//     time to fill a tenth of the thinnest regolith (a dry aquifer has no
//     stability limit). The cap grows by 1.5 per substep while the residual
//     falls and halves (not below its start) when it rises; the Courant and
//     von Neumann limits always apply. Only the end state is meaningful.
//   * Seepage is reported as the substep's MEAN, f - div q - n dh/dt from the
//     unclipped update, so the water budget closes to round-off; Landlab
//     reports the regularised rate at the end of each step instead.
//   * Fixed-head cells (channels in lsdtt-slope-stability) keep their
//     thickness; flow into them leaves the aquifer. Nodata cells, cells with no
//     regolith (d <= 0) and the raster edge are no-flux.
//
// PARALLELISM: links are stored per cell (east and south of each cell), so a
// substep is two passes over rows: fluxes and the stability limit (a min
// reduction), then the per-cell update, gathering the four link fluxes. Budget
// terms are summed per row and the rows added in order, so every result,
// budget included, is identical for any thread count. Grids under 16,384
// cells run serially: a substep there is shorter than the thread fork.
//
// Double precision throughout (a steady-state stopping test is below float
// resolution otherwise).
// -----------------------------------------------------------------------------

#pragma once

#include <cstdint>
#include <string>
#include <vector>

#include "../core/lsdtt_raster_meta.hpp"

namespace lsdtt
{

enum class groundwater_transmissivity_profile
{
    uniform,
    exponential
};

struct groundwater_dupuit_parameters
{
    double hydraulic_conductivity_m_per_s = 1.0e-3;   // K, or K0 at the surface (exponential)
    double drainable_porosity = 0.2;                   // n
    groundwater_transmissivity_profile profile = groundwater_transmissivity_profile::uniform;
    double conductivity_decay_depth_metres = 1.0;      // z_m (exponential only)
    double regularisation = 0.01;                      // eps
    double courant_coefficient = 0.5;
    double von_neumann_coefficient = 0.8;
};

// Throws std::invalid_argument naming the field.
void validate_groundwater_dupuit_parameters(const groundwater_dupuit_parameters& p);

// Landlab's semi-analytic thickness update over dt, for regolith thickness d,
// recharge f, divergence div_q and porosity n, with regularisation eps.
// Exposed for the unit test. Not clipped.
[[nodiscard]] double dupuit_update_thickness(
    double dt, double h0, double d, double f, double div_q, double n, double eps);

// Cumulative water budget, in m^3 over the run so far.
struct groundwater_budget
{
    double recharge_in = 0.0;        // recharge on active (not fixed-head) cells
    double storage_change = 0.0;     // n * dh * cell area, active cells
    double outflow_fixed_head = 0.0; // into fixed-head cells (channels)
    double seepage_out = 0.0;        // return flow to the surface
    double clipped_added = 0.0;      // water added by clipping h at 0
    // recharge_in + clipped_added = storage_change + outflow_fixed_head + seepage_out
};

class lsdtt_groundwater_dupuit
{
public:
    // surface, base: elevations (m), one per cell, row-major on meta's grid.
    // valid: nonzero = cell in the domain.
    // fixed_head: empty, or nonzero = fixed-head cell (thickness held).
    // fixed_thickness: empty (fixed cells held at h = d, the water table at the
    //   surface), or one thickness per cell, used at fixed-head cells.
    lsdtt_groundwater_dupuit(const lsdtt_raster_meta& meta,
                             const std::vector<double>& surface,
                             const std::vector<double>& base,
                             const std::vector<std::uint8_t>& valid,
                             const std::vector<std::uint8_t>& fixed_head,
                             const std::vector<double>& fixed_thickness,
                             const groundwater_dupuit_parameters& params,
                             int n_threads = 0);   // 0 = OpenMP default

    // Initial thickness of the active cells (fixed cells keep theirs);
    // clipped to [0, d].
    void set_thickness(const std::vector<double>& h);

    // Advance by duration seconds with recharge f (m/s, one per cell).
    // Returns the number of substeps taken.
    std::int64_t run(double duration, const std::vector<double>& recharge);

    struct steady_report
    {
        bool converged = false;
        std::int64_t substeps = 0;     // explicit: substeps; implicit: pseudo-time (Newton) steps
        double time = 0.0;             // pseudo-time marched (s)
        double residual = 0.0;         // explicit: max n |dh/dt| / max f on the last substep;
                                       // implicit: the steady measure of run_to_steady_implicit
        // Implicit only:
        std::int64_t linear_iterations = 0;   // BiCGSTAB iterations (0 with the direct solver)
        std::int64_t line_search_cuts = 0;    // step halvings in the line search
        std::int64_t dt_cuts = 0;             // pseudo-time steps retried with a smaller dt
        std::int64_t active_set_changes = 0;  // cells whose bound status changed, summed over steps
        std::int64_t n_saturated = 0;         // unknowns at h = d at the end
        std::int64_t n_dry = 0;               // unknowns at h = 0 at the end
        bool used_direct_solver = false;
        std::string message;                  // why it stopped, when not converged
    };
    // March to steady state under constant recharge: stops when
    // max |n dh/dt| <= relative_tolerance * max(f). Throws if max(f) <= 0.
    steady_report run_to_steady(const std::vector<double>& recharge,
                                double relative_tolerance = 1.0e-6,
                                std::int64_t max_substeps = 10000000);

    // ---- Implicit steady state -------------------------------------------------
    // Solves the steady problem directly: on every unknown, 0 <= h <= d and
    //   F = 0 where 0 < h < d,  F <= 0 where h = d (seepage -F),  F >= 0 where h = 0,
    // F the steady residual (implicit_residual with inv_dt = 0). This is a
    // box-constrained complementarity problem; the regularised seepage of the
    // explicit march converges to the same state.
    //
    // METHOD. Pseudo-transient continuation (Kelley & Keyes 1998): each pseudo-
    // time step is ONE semismooth Newton step (primal-dual active set;
    // Hintermuller, Ito & Kunisch 2002) on the backward-Euler problem
    //   Phi(h) = h - P_[0,d](h - c Ft(h)) = 0,  Ft = n (h - h_k) / dt + F(h),
    // c = min(dt / n, dx^2 / T_max). Cells where h - c Ft leaves [0, d] are
    // active: their row becomes h = bound; the rest take the Jacobian row. The
    // step is projected onto [0, d] and backtracked (Armijo on ||Phi||^2, at
    // most max_line_search halvings); if that fails, dt is quartered and the
    // step retried. dt then grows by switched evolution relaxation
    // (Mulder & van Leer 1985): dt *= min(dt_growth_max, m_old / m_new), with
    // m the steady measure below; once n / dt is negligible the iteration is
    // Newton on the steady problem.
    //
    // CONVERGENCE. m = max over unknowns of |F| (0 < h < d), max(F, 0) (h = d) and
    // max(-F, 0) (h = 0), divided by max f; converged when m <= tolerance.
    //
    // LINEAR SOLVES. Eigen SparseLU for up to direct_solve_max_unknowns
    // unknowns, else BiCGSTAB with IncompleteLUT (both serial, so the result
    // does not depend on the thread count). Without Eigen it reports failure.
    //
    // Starts from the current thickness (set_thickness, else 0). Afterwards
    // thickness() is the solution, seepage() is -F on saturated cells (m/s),
    // and steady_rates() the steady budget. budget() is not changed.
    struct implicit_options
    {
        double relative_tolerance = 1.0e-6;     // as run_to_steady
        int max_steps = 200;                    // pseudo-time (Newton) steps
        double dt_initial = 0.0;                // s; 0 = 0.1 n d_min / f_max, as the explicit march
        double dt_growth_max = 10.0;
        double dt_max = 1.0e30;                 // s
        int solver_choice = -1;                 // -1 auto, 0 direct, 1 iterative
        long direct_solve_max_unknowns = 10000;   // measured: BiCGSTAB + ILUT is faster from ~4e4 unknowns
        double linear_tolerance = 1.0e-10;      // BiCGSTAB relative residual
        int max_line_search = 8;
        // run_implicit (transient) only:
        double max_dt_s = 3600.0;               // longest backward-Euler step (s)
        double dt_tolerance_metres = 0.0;            // step-doubling error control (m); 0 = off
        int max_newton_per_step = 30;
        int max_splits = 8;                     // halvings of a step that fails to converge
    };
    steady_report run_to_steady_implicit(const std::vector<double>& recharge, const implicit_options& options);
    steady_report run_to_steady_implicit(const std::vector<double>& recharge)
    {
        return run_to_steady_implicit(recharge, implicit_options{});
    }

    // Steady budget of the last implicit steady solve, in m^3/s:
    //   recharge = outflow_fixed_head + seepage + imbalance (imbalance from the
    //   residual left on the unknowns; ~ tolerance).
    struct steady_budget
    {
        double recharge = 0.0, outflow_fixed_head = 0.0, seepage = 0.0, imbalance = 0.0;
    };
    [[nodiscard]] const steady_budget& steady_rates() const { return steady_rates_; }

    // ---- Implicit transient stepping -----------------------------------------
    // Advance by duration seconds with recharge f by backward Euler: each step
    // solves n (h - h_old) / dt + F(h) = 0 with 0 <= h <= d (the same semismooth
    // Newton iteration as run_to_steady_implicit, h_old fixed, active-set scale
    // c = dt / n, to max |Phi| <= 1e-12 max(h_old, 1 mm) in thickness, or a
    // round-off stall within 1e-9 of it). Steps are at most max_dt_s; a step
    // that fails to converge is halved, at most max_splits times. With
    // dt_tolerance_metres > 0, each step is also taken as two half steps; where the
    // two differ by more than dt_tolerance_metres the step is halved, otherwise the
    // half-step result is kept and the next step may double (up to max_dt_s).
    // Call once per constant-recharge segment so steps land on its ends.
    //
    // Per step the seepage is -Ft on saturated cells at the new state, so
    //   recharge dt A = storage change + outflow to fixed heads dt + seepage dt
    // to the Newton tolerance; budget() accumulates these (m^3) as run() does,
    // seepage() is the mean over this call. Backward Euler is first order in dt
    // and damps fast transients: max_dt_s sets the resolution.
    struct implicit_run_report
    {
        bool ok = true;
        std::int64_t steps = 0;
        std::int64_t newton_iterations = 0;
        std::int64_t linear_iterations = 0;
        std::int64_t splits = 0;              // halvings after a failed step
        std::int64_t rejected_by_error = 0;   // halvings by the step-doubling test
        double max_error_estimate_metres = 0.0;    // largest accepted step-doubling difference
        double time = 0.0;                    // s advanced
        std::string message;
    };
    implicit_run_report run_implicit(double duration, const std::vector<double>& recharge,
                                     const implicit_options& options);
    implicit_run_report run_implicit(double duration, const std::vector<double>& recharge)
    {
        return run_implicit(duration, recharge, implicit_options{});
    }

    [[nodiscard]] const std::vector<double>& thickness() const { return h_; }
    [[nodiscard]] std::vector<double> water_table_depth() const;   // z - z_w (m); NaN off the domain
    [[nodiscard]] const std::vector<double>& seepage() const { return seepage_; }   // m/s, mean over the last run
    [[nodiscard]] const std::vector<double>& regolith_thickness() const { return d_; }
    [[nodiscard]] const groundwater_budget& budget() const { return budget_; }
    [[nodiscard]] double storage() const;   // n h A over active cells (m^3)

    // ---- Building blocks of the implicit (Newton) solves ---------------------
    // The unknowns are the active cells (not fixed-head, regolith > 0), in
    // row-major order. For a full-grid thickness h (values at fixed-head and
    // inactive cells are ignored: fixed cells keep their held thickness), the
    // residual of unknown i is
    //
    //   F_i(h) = n (h_i - h_old_i) inv_dt + (div q(h))_i - f_i      (m/s)
    //
    // with div q from exactly the link fluxes of the explicit substep (upwind
    // thickness, cos^2 of the base, fixed heads, no-flux edges and nodata).
    // inv_dt = 0 gives the steady residual (h_old is then not read). At steady
    // state F_i = 0 where 0 < h_i < d_i, and F_i <= 0 (seepage -F_i) where h_i = d_i.
    //
    // The Jacobian dF/dh is 5-point and non-symmetric. For a link (a -> b),
    // q = T(h_up) c2 (z_w,a - z_w,b) / dx with up = a if z_w,a >= z_w,b, and
    //   dq/dh_a =  T c2 / dx + [up = a] T'(h_a) c2 (z_w,a - z_w,b) / dx
    //   dq/dh_b = -T c2 / dx + [up = b] T'(h_b) c2 (z_w,a - z_w,b) / dx,
    // T' = K (uniform) or K0 exp(-(d - h) / z_m) (exponential), the derivative
    // from above at h = 0. At a tie (z_w,a = z_w,b) the upwind side is a, as in
    // the explicit kernel: one element of the generalised (Clarke) Jacobian.
    struct jacobian_entry
    {
        std::int64_t row, col;   // unknown indices
        double value;
    };
    [[nodiscard]] std::size_t n_unknowns() const { return unknown_cell_.size(); }
    [[nodiscard]] const std::vector<std::int64_t>& unknown_cell() const { return unknown_cell_; }   // cell of each unknown
    [[nodiscard]] const std::vector<std::int64_t>& cell_unknown() const { return cell_unknown_; }   // unknown of each cell, -1
    // residual: resized to n_unknowns().
    void implicit_residual(const std::vector<double>& h, const std::vector<double>& h_old, double inv_dt,
                           const std::vector<double>& recharge, std::vector<double>& residual) const;
    // entries: cleared and filled row by row (row-major unknown order, then the
    // order of the links east, west, south, north); duplicates are summed by
    // the caller (Eigen setFromTriplets does). The diagonal is always present.
    void implicit_jacobian(const std::vector<double>& h, double inv_dt, std::vector<jacobian_entry>& entries) const;
    // Flow from active cells into fixed-head cells for thickness h (m^3/s).
    [[nodiscard]] double outflow_to_fixed_head(const std::vector<double>& h) const;

private:
    double transmissivity(std::size_t cell, double h) const;
    double transmissivity_derivative(std::size_t cell, double h) const;
    double active_set_scale_reference() const;   // dx^2 / T_max (s)
    double held(const std::vector<double>& h, std::size_t i) const { return fixed_[i] ? h_[i] : h[i]; }
    // Flux q (positive from a to b) per unit width and its derivatives, for
    // the link a -> b with cos^2 c2, at thickness h.
    void link_flux(std::size_t a, std::size_t b, double c2, const std::vector<double>& h,
                   double& q, double& dq_da, double& dq_db) const;
    // One substep of at most dt_max; returns the dt taken and max n|dh/dt|.
    // Seepage: accumulate seep * dt into seep_accum_ (run), or overwrite
    // seepage_ with this substep's mean rate (run_to_steady).
    enum class seep_mode { accumulate, overwrite };
    double substep(double dt_max, const std::vector<double>& recharge, double& max_rate, seep_mode mode);

    lsdtt_raster_meta meta_;
    std::size_t n_;
    std::vector<double> z_, zb_, d_, h_, seepage_;
    std::vector<std::uint8_t> active_, fixed_;
    // Per-cell links: c2_east_[i] / c2_south_[i] = cos^2(alpha) of the link from
    // i to its east / south neighbour, or 0 for no link. q_east_ / q_south_:
    // flux per unit width, positive from i to the neighbour (work arrays).
    std::vector<double> c2_east_, c2_south_, q_east_, q_south_, seep_accum_;
    std::vector<double> row_terms_;   // per-row budget partial sums (work)
    std::vector<std::int64_t> unknown_cell_, cell_unknown_;
    int n_threads_ = 0;
    bool parallel_ = false;   // threads only on grids large enough to repay the fork (k_parallel_min_cells)
    groundwater_dupuit_parameters p_;
    groundwater_budget budget_;
    steady_budget steady_rates_;
};

} // namespace lsdtt
