// lsdtt_landslide_crn_cpu.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
//
// See lsdtt_landslide_crn_cpu.hpp. No new CRN physics is defined here.

#include "lsdtt_landslide_crn_cpu.hpp"

#include <cmath>
#include <map>
#include <stdexcept>
#include <vector>

#include "lsdtt_crn_column.hpp"

namespace lsdtt
{

namespace
{

double pick_nuclide(const CrnColumnSample& s, const std::string& nuclide)
{
    if (nuclide == "Be10") return s.conc_10Be;
    if (nuclide == "Al26") return s.conc_26Al;
    if (nuclide == "C14")  return s.conc_14C;
    throw std::invalid_argument(
        "compute_landslide_surface_concentration: unknown nuclide \"" + nuclide +
        "\". Options: Be10, Al26, C14.");
}

// Effective erosion rate (g/cm2/yr) -> physical rate (m/yr).
// to_eff() in lsdtt_crn_column is depth_m * rho_kg_m3 * 0.1, so the inverse is
// eff / (rho * 0.1).
double eff_to_m_per_yr(double eff_g_cm2_yr, double rho_kg_m3)
{
    return eff_g_cm2_yr / (rho_kg_m3 * 0.1);
}

}  // namespace


double landslide_crn_background_concentration(const LandslideCrnConfig& config)
{
    lsdtt_crn_column col(config.rock_density_kg_m3);
    const CrnColumnSample s = col.steady_state(
        0.0, config.background_eff_erate_g_cm2_yr,
        config.production_scaling, config.total_shielding, config.muon_scheme);
    return pick_nuclide(s, config.nuclide);
}


lsdtt_raster_cpu compute_landslide_surface_concentration(
    const lsdtt_raster_cpu&   excavation_depth_m,
    const lsdtt_raster_cpu&   time_of_last_failure,
    const LandslideCrnConfig& config)
{
    const lsdtt_raster_meta& meta = excavation_depth_m.meta();
    const lsdtt_raster_meta& tmeta = time_of_last_failure.meta();
    if (meta.n_rows != tmeta.n_rows || meta.n_cols != tmeta.n_cols)
        throw std::invalid_argument(
            "compute_landslide_surface_concentration: excavation depth and "
            "time-of-last-failure rasters differ in size");

    const auto& depth = excavation_depth_m.get_data_array();
    const auto& dmask = excavation_depth_m.get_mask_array();
    const auto& tfail = time_of_last_failure.get_data_array();
    const auto& tmask = time_of_last_failure.get_mask_array();

    const float depth_ndv = static_cast<float>(meta.no_data_value);
    const float tfail_ndv = static_cast<float>(tmeta.no_data_value);

    lsdtt_crn_column col(config.rock_density_kg_m3);
    const double eps_m_yr =
        eff_to_m_per_yr(config.background_eff_erate_g_cm2_yr, config.rock_density_kg_m3);

    // Undisturbed surface: the value every unfailed pixel takes.
    const double background = landslide_crn_background_concentration(config);

    // One transient solve per DISTINCT (depth, time-since-failure). Every pixel
    // in a scar shares both, so this collapses thousands of pixels onto a
    // handful of solves. Keys are quantised to avoid float noise splitting what
    // is physically one scar.
    std::map<std::pair<long long, long long>, double> cache;

    lsdtt_array2d_f  out(meta.n_rows, meta.n_cols, depth_ndv);
    lsdtt_array2d_u8 omask(meta.n_rows, meta.n_cols, 0);

    for (int r = 0; r < meta.n_rows; ++r)
    {
        for (int c = 0; c < meta.n_cols; ++c)
        {
            if (!dmask(r, c)) continue;

            const bool failed =
                tmask(r, c) && tfail(r, c) != tfail_ndv &&
                dmask(r, c) && depth(r, c) != depth_ndv && depth(r, c) > 0.0f;

            if (!failed)
            {
                out(r, c)   = static_cast<float>(background);
                omask(r, c) = 1;
                continue;
            }

            const double d_m      = static_cast<double>(depth(r, c));
            const double t_since  = config.t_now_yr - static_cast<double>(tfail(r, c));

            if (!(t_since > 0.0))
            {
                // Failed at (or after) the reporting time: the surface is the
                // rock that was at depth d, with no re-accumulation yet. That
                // is the inherited concentration, straight from steady state.
                const CrnColumnSample s = col.steady_state(
                    d_m * config.rock_density_kg_m3 * 0.1,
                    config.background_eff_erate_g_cm2_yr,
                    config.production_scaling, config.total_shielding,
                    config.muon_scheme);
                out(r, c)   = static_cast<float>(pick_nuclide(s, config.nuclide));
                omask(r, c) = 1;
                continue;
            }

            const std::pair<long long, long long> key{
                static_cast<long long>(std::llround(d_m * 1000.0)),
                static_cast<long long>(std::llround(t_since))};

            auto it = cache.find(key);
            if (it == cache.end())
            {
                // Erosion history, row 0 MOST RECENT, yr BP increasing:
                //   [0] t_since  : background erosion from failure to now
                //   [1] t_since  : the landslide -- instantaneous removal of
                //                  d_m metres (rate = 0 marks it as a removal,
                //                  and it legally shares the time above it)
                //   [2] long     : background erosion before the failure
                //   [3] -9999    : background rate sentinel
                const std::vector<double> t_yr   {t_since, t_since, 1.0e6, -9999.0};
                const std::vector<double> rate   {eps_m_yr, 0.0, eps_m_yr, eps_m_yr};
                const std::vector<double> removal{0.0, d_m, 0.0, 0.0};

                lsdtt_crn_column c2(config.rock_density_kg_m3);
                c2.set_erosion_time_series(t_yr, rate, removal);

                // Initial concentrations left at the default (-1), which starts
                // the column at steady state under the background rate -- the
                // pre-failure profile, which is what inheritance requires.
                const std::vector<CrnColumnSample> prof = c2.run_transient(
                    std::vector<double>{0.0},
                    config.production_scaling * config.total_shielding,
                    config.muon_scheme);

                const double v = prof.empty() ? background
                                              : pick_nuclide(prof[0], config.nuclide);
                it = cache.emplace(key, v).first;
            }

            out(r, c)   = static_cast<float>(it->second);
            omask(r, c) = 1;
        }
    }

    return lsdtt_raster_cpu(std::move(out), std::move(omask), meta);
}

}  // namespace lsdtt
