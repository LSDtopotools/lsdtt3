// =============================================================================
// lsdtt_lem_timestep_cpu.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
//
// See the header for the criterion and the measurements behind it.
// =============================================================================

#include "lsdtt_lem_timestep_cpu.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <sstream>
#include <vector>

namespace lsdtt
{

// -----------------------------------------------------------------------------
double median_interior_gradient(const lsdtt_raster_cpu& dem)
{
    const lsdtt_raster_meta& meta = dem.get_meta();
    const int nrows = meta.n_rows, ncols = meta.n_cols;
    const double dx = meta.cell_size;
    if (nrows < 3 || ncols < 3 || !(dx > 0.0)) return 0.0;

    const std::vector<float>&   z    = dem.get_data();
    const std::vector<uint8_t>& mask = dem.get_mask();

    std::vector<double> g;
    g.reserve(static_cast<std::size_t>(nrows) * ncols / 2);

    for (int r = 1; r < nrows - 1; ++r)
        for (int c = 1; c < ncols - 1; ++c)
        {
            const std::size_t f = static_cast<std::size_t>(r) * ncols + c;
            const std::size_t fn = f - ncols, fs = f + ncols, fw = f - 1, fe = f + 1;
            if (!mask[f] || !mask[fn] || !mask[fs] || !mask[fw] || !mask[fe]) continue;
            const double gx = (static_cast<double>(z[fe]) - static_cast<double>(z[fw]))
                            / (2.0 * dx);
            const double gy = (static_cast<double>(z[fs]) - static_cast<double>(z[fn]))
                            / (2.0 * dx);
            g.push_back(std::sqrt(gx * gx + gy * gy));
        }

    if (g.empty()) return 0.0;
    std::nth_element(g.begin(), g.begin() + g.size() / 2, g.end());
    return g[g.size() / 2];
}


// -----------------------------------------------------------------------------
lsdtt_lem_timestep_diagnosis diagnose_lem_timestep(
    const lsdtt_raster_cpu& dem,
    double                  uplift_rate,
    double                  timestep_yr,
    double                  target_ratio)
{
    lsdtt_lem_timestep_diagnosis out;

    const double dx = dem.get_meta().cell_size;
    out.reference_slope = median_interior_gradient(dem);

    // A flat or near-flat starting surface carries no gradient information, so
    // there is nothing to compare against. Say nothing rather than invent a
    // reference: a spurious warning on every diamond-square run would train the
    // user to ignore this one.
    if (!(out.reference_slope > 1e-6) || !(dx > 0.0) ||
        !(uplift_rate > 0.0) || !(timestep_yr > 0.0))
        return out;

    out.cell_relief_m     = dx * out.reference_slope;
    out.uplift_per_step_m = uplift_rate * timestep_yr;
    out.ratio             = out.uplift_per_step_m / out.cell_relief_m;
    out.suggested_dt_yr   = target_ratio * out.cell_relief_m / uplift_rate;
    out.too_large         = (out.ratio > target_ratio);

    if (out.too_large)
    {
        std::ostringstream m;
        m.setf(std::ios::fixed);
        m << std::setprecision(2)
          << "one timestep of uplift (" << out.uplift_per_step_m
          << " m) is " << out.ratio << "x the relief across a single cell ("
          << out.cell_relief_m << " m at median gradient "
          << out.reference_slope << "). Above ~" << target_ratio
          << " the operator split stops resolving incision against diffusion and "
             "manufactures grid-scale roughness -- the run will still converge and "
             "still report ~100% of the uplift rate, but hillslope length can be "
             "wrong by several fold. Suggested timestep: <= "
          << std::setprecision(0) << out.suggested_dt_yr
          << " yr, or enable adaptive_timestep.";
        out.warning = m.str();
    }
    return out;
}


// -----------------------------------------------------------------------------
lsdtt_lem_timestep_controller::lsdtt_lem_timestep_controller(
    double dt_initial_yr,
    const lsdtt_lem_timestep_controller_options& opts)
    : opts_(opts), dt_(dt_initial_yr), dt_raw_(dt_initial_yr),
      dt_initial_(dt_initial_yr)
{
    dt_min_ = (opts.dt_min_yr > 0.0) ? opts.dt_min_yr : dt_initial_yr / 100.0;
    dt_max_ = (opts.dt_max_yr > 0.0) ? opts.dt_max_yr : dt_initial_yr * 100.0;
}


double lsdtt_lem_timestep_controller::update(double max_local_dz_m,
                                             const lsdtt_raster_cpu& dem)
{
    if (n_steps_ % std::max(1, opts_.slope_refresh_interval) == 0)
        slope_ = median_interior_gradient(dem);
    ++n_steps_;

    const double dx = dem.get_meta().cell_size;
    const double cell_relief = dx * slope_;

    // Nothing to key on: no gradient yet (a flat start), or no movement. Leave the
    // step alone rather than react to a ratio built from noise.
    if (!(cell_relief > 1e-9) || !(max_local_dz_m > 0.0)) return dt_;

    last_ratio_ = max_local_dz_m / cell_relief;

    // Decide on a SMOOTHED ratio, not the raw one.
    //
    // max_local_elevation_change is a maximum over the grid, so it jumps around
    // from step to step even in a perfectly steady run -- one cell somewhere is
    // always having an unusual step. Reacting to that noise makes the controller
    // oscillate across a ladder rung and pay a refactorisation each time (13 over
    // 60 steps, measured). An exponential moving average keeps genuine drift and
    // drops the jitter. The RAW ratio is still what gets reported, because that is
    // what the criterion is actually about.
    ratio_ema_ = (ratio_ema_ > 0.0)
        ? (1.0 - opts_.ratio_smoothing) * ratio_ema_
              + opts_.ratio_smoothing * last_ratio_
        : last_ratio_;
    const double r_dec = ratio_ema_;

    const double dt_before = dt_;

    // The controller accumulates on a CONTINUOUS value and only the reported
    // timestep is quantised. Applying a 1.15x growth directly to a quantised dt
    // never reaches the next rung of a 2x ladder, so the step would be pinned
    // forever one rung below where it belongs -- which is what happened on the
    // first attempt: 54 consecutive grow decisions with dt stuck at 500 yr and the
    // ratio sitting at 0.197 against a target of 0.7.
    if (r_dec > opts_.target_ratio)
    {
        const double scale = std::max(opts_.shrink_factor,
                                      opts_.target_ratio / r_dec);
        dt_raw_ = std::max(dt_min_, dt_raw_ * scale);
        ++n_shrinks_;
    }
    else if (r_dec < 0.5 * opts_.target_ratio)
    {
        dt_raw_ = std::min(dt_max_, dt_raw_ * opts_.grow_factor);
        ++n_grows_;
    }

    dt_ = dt_raw_;

    // Snap to the ladder, with a FULL RUNG of hysteresis.
    //
    // Rounding to the nearest rung is not enough: the continuous value then sits
    // near a boundary and re-crosses it repeatedly, and every crossing costs a
    // refactorisation downstream. Requiring a whole rung of travel before moving
    // means the step changes only when it genuinely belongs at a different scale
    // (9 factorisations over 60 steps became 2 in the measured case).
    if (opts_.ladder_ratio > 1.0 && dt_initial_ > 0.0)
    {
        const double r = opts_.ladder_ratio;

        // ASYMMETRIC, for the same reason the growth factor is: an oversized step
        // corrupts the landscape, an undersized one only costs time.
        //
        // DOWN promptly -- as soon as the achieved ratio breaches the target, drop
        // a rung. Symmetric full-rung hysteresis was tried first and left the step
        // parked one rung too high, running at a ratio of 0.93 against a target of
        // 0.7; saving a refactorisation is not worth running outside the criterion
        // the controller exists to enforce.
        //
        // UP reluctantly -- a full rung of travel, so growth cannot oscillate
        // across a boundary and pay a refactorisation each time.
        if (r_dec > opts_.target_ratio && rung_ > lowest_rung_seen_ - 1)
            --rung_;
        while (dt_raw_ >= dt_initial_ * std::pow(r, rung_ + 1)) ++rung_;
        while (dt_raw_ <= dt_initial_ * std::pow(r, rung_ - 1)) --rung_;
        lowest_rung_seen_ = std::min(lowest_rung_seen_, rung_);
        dt_ = std::min(dt_max_,
                       std::max(dt_min_,
                                dt_initial_ * std::pow(r, rung_)));
    }

    if (dt_ != dt_before) ++n_changes_;
    return dt_;
}


// -----------------------------------------------------------------------------
// lsdtt_lem_convergence_monitor
// -----------------------------------------------------------------------------
lsdtt_lem_convergence_monitor::lsdtt_lem_convergence_monitor(
    double dt_initial_yr, const lsdtt_lem_convergence_options& opts)
    : opts_(opts), dt_initial_(dt_initial_yr)
{
    history_.reserve(static_cast<std::size_t>(std::max(8, opts.stall_window)) + 1);
}


lem_convergence_state lsdtt_lem_convergence_monitor::update(
    double max_local_dz_m, double dt_used_yr)
{
    ++n_steps_;
    last_dz_ = max_local_dz_m;

    history_.push_back(max_local_dz_m);
    const std::size_t w = static_cast<std::size_t>(std::max(8, opts_.stall_window));
    if (history_.size() > w) history_.erase(history_.begin());

    // ---- converged? -------------------------------------------------------
    //
    // The dt guard is the point: an absolute per-step tolerance is trivially met
    // by a timestep approaching zero, so a run whose controller is collapsing must
    // not be allowed to call that steady state.
    const bool dt_ok = (dt_initial_ <= 0.0) ||
                       (dt_used_yr >= opts_.min_dt_fraction * dt_initial_);

    // Both tolerances must be met when both are set. The rate is the one that
    // means the same thing at any timestep; the per-step one is kept because it is
    // what a user naturally specifies.
    const bool below_abs  = (opts_.dz_tolerance_m <= 0.0) ||
                            (max_local_dz_m < opts_.dz_tolerance_m);
    const bool below_rate = (opts_.dz_rate_tolerance_m_per_yr <= 0.0) ||
                            (dt_used_yr <= 0.0) ||
                            (max_local_dz_m / dt_used_yr
                                 < opts_.dz_rate_tolerance_m_per_yr);

    if (below_abs && below_rate && dt_ok) ++n_below_;
    else                                  n_below_ = 0;

    if (n_below_ >= opts_.patience_steps)
    {
        state_ = lem_convergence_state::converged;

        // Name the criterion that was ACTUALLY applied. Reporting a disabled
        // tolerance reads as "no cell moved more than 0.000000 m", which is both
        // false and alarming.
        std::ostringstream m;
        m << "steady state for " << opts_.patience_steps
          << " consecutive steps (p99.9 of |dz|";
        if (opts_.dz_tolerance_m > 0.0)
            m << " < " << opts_.dz_tolerance_m << " m/step";
        if (opts_.dz_rate_tolerance_m_per_yr > 0.0)
            m << (opts_.dz_tolerance_m > 0.0 ? " and" : "")
              << " < " << opts_.dz_rate_tolerance_m_per_yr << " m/yr";
        m << "; last " << last_dz_ << " m/step)";
        message_ = m.str();
        return state_;
    }

    // ---- stalled? ---------------------------------------------------------
    if (history_.size() < w) return state_;

    // Compare the best of the first fifth of the window with the best of the last
    // fifth. Using minima rather than means keeps a single noisy step from masking
    // real progress, which a max statistic is prone to.
    const std::size_t q = std::max<std::size_t>(1, w / 5);
    const double early = *std::min_element(history_.begin(), history_.begin() + q);
    const double late  = *std::min_element(history_.end() - q, history_.end());
    decay_ = (early > 0.0) ? late / early : 1.0;

    if (decay_ < opts_.stall_decay_threshold) return state_;   // still improving

    // Not improving. Alternating steps mean a period-2 cycle rather than slow
    // convergence, and that distinction changes what the user should do about it.
    long alt = 0;
    for (std::size_t i = 2; i < history_.size(); ++i)
    {
        const double d1 = history_[i - 1] - history_[i - 2];
        const double d2 = history_[i]     - history_[i - 1];
        if (d1 * d2 < 0.0) ++alt;
    }
    oscillating_ = (static_cast<double>(alt) /
                    static_cast<double>(history_.size() - 2))
                   > opts_.oscillation_alternation_fraction;

    state_ = lem_convergence_state::stalled;
    message_ = oscillating_
        ? "OSCILLATING: max|dz| has not fallen over the last "
          + std::to_string(w) + " steps and alternates step to step -- a period-2 "
          "cycle, which points at the numerics (timestep too large, or a limit "
          "cycle in the nonlinear solve) rather than at needing more model time"
        : "STALLED: max|dz| has not fallen over the last " + std::to_string(w)
          + " steps but is still above tolerance -- the run is not converging";
    return state_;
}

} // namespace lsdtt
