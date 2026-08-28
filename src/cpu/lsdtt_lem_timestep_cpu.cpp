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
    : opts_(opts), dt_(dt_initial_yr)
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

    if (last_ratio_ > opts_.target_ratio)
    {
        // Overshoot: shrink, and scale the shrink to how far out we are so a
        // badly-oversized step is corrected in one move rather than a dozen.
        const double scale = std::max(opts_.shrink_factor,
                                      opts_.target_ratio / last_ratio_);
        dt_ = std::max(dt_min_, dt_ * scale);
        ++n_shrinks_;
    }
    else if (last_ratio_ < 0.5 * opts_.target_ratio)
    {
        // Comfortably inside: grow slowly. Growth is deliberately much gentler
        // than shrinkage -- an oversized step corrupts the landscape, an
        // undersized one only costs time.
        dt_ = std::min(dt_max_, dt_ * opts_.grow_factor);
        ++n_grows_;
    }

    return dt_;
}

} // namespace lsdtt
