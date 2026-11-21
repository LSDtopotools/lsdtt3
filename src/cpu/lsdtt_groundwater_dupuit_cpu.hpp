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
        std::int64_t substeps = 0;
        double time = 0.0;             // pseudo-time marched (s)
        double residual = 0.0;         // max n |dh/dt| / max f on the last substep
    };
    // March to steady state under constant recharge: stops when
    // max |n dh/dt| <= relative_tolerance * max(f). Throws if max(f) <= 0.
    steady_report run_to_steady(const std::vector<double>& recharge,
                                double relative_tolerance = 1.0e-6,
                                std::int64_t max_substeps = 10000000);

    [[nodiscard]] const std::vector<double>& thickness() const { return h_; }
    [[nodiscard]] std::vector<double> water_table_depth() const;   // z - z_w (m); NaN off the domain
    [[nodiscard]] const std::vector<double>& seepage() const { return seepage_; }   // m/s, mean over the last run
    [[nodiscard]] const std::vector<double>& regolith_thickness() const { return d_; }
    [[nodiscard]] const groundwater_budget& budget() const { return budget_; }
    [[nodiscard]] double storage() const;   // n h A over active cells (m^3)

private:
    double transmissivity(std::size_t cell, double h) const;
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
    int n_threads_ = 0;
    bool parallel_ = false;   // threads only on grids large enough to repay the fork (k_parallel_min_cells)
    groundwater_dupuit_parameters p_;
    groundwater_budget budget_;
};

} // namespace lsdtt
