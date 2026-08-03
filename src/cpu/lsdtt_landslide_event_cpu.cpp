// lsdtt_landslide_event_cpu.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
//
// See lsdtt_landslide_event_cpu.hpp for the design.

#include "lsdtt_landslide_event_cpu.hpp"

#include <algorithm>
#include <cmath>
#include <functional>
#include <iomanip>
#include <limits>
#include <sstream>

namespace lsdtt
{

namespace
{

inline double next_uniform(std::uint64_t& state)
{
    state ^= state >> 12;
    state ^= state << 25;
    state ^= state >> 27;
    const std::uint64_t v = state * 2685821657736338717ULL;
    return static_cast<double>(v >> 11) * (1.0 / 9007199254740992.0);
}

}  // namespace


std::string LandslideEventResult::summary() const
{
    std::ostringstream os;
    os << std::fixed << std::setprecision(2);
    os << "landslide event: " << n_landslides << " placed";
    if (n_failed > 0) os << " (" << n_failed << " unplaceable)";
    os << "\n  volume requested " << requested_volume_m3 << " m3"
       << ", achieved " << achieved_volume_m3 << " m3";
    if (requested_volume_m3 > 0.0)
        os << " (" << 100.0 * achieved_volume_m3 / requested_volume_m3 << "%)";
    os << "\n  budget solve: " << (budget_converged ? "converged" : "NOT CONVERGED")
       << " in " << budget_iterations << " iteration(s)"
       << ", E[V] = " << mean_volume_used_m3 << " m3"
       << ", A_max = " << area_max_used_m2 << " m2";
    os << "\n  distribution truncated at A_max: "
       << std::setprecision(4) << 100.0 * truncated_fraction_number
       << "% by number, " << std::setprecision(2)
       << 100.0 * truncated_fraction_volume << "% by volume";
    if (rollover_unresolved)
        os << "\n  WARNING: cell size cannot resolve the small-area rollover; "
              "the sampled distribution is effectively a truncated power law";
    if (n_failed > 0)
        os << "\n  WARNING: " << n_failed << " landslide(s) could not be placed; "
              "the achieved volume is short of the sampled volume";
    return os.str();
}


std::vector<int> solve_landslide_count_for_budget(
    const LandslideSizeSampler& sampler,
    double target_volume_m3,
    double tol,
    int    max_iter,
    std::uint64_t& rng_state,
    bool&  converged,
    int&   iterations_used,
    double& sampled_volume_m3)
{
    converged        = false;
    iterations_used  = 0;
    sampled_volume_m3 = 0.0;

    std::vector<int> best;
    double best_rel = std::numeric_limits<double>::max();
    double best_vol = 0.0;

    if (!(target_volume_m3 > 0.0)) return best;

    const double mean_v = sampler.mean_volume_m3();
    if (!(mean_v > 0.0)) return best;

    // First guess from the mean. E[V] is computed over the TRUNCATED pmf, so
    // this already accounts for the A_max cutoff.
    long n = static_cast<long>(std::llround(target_volume_m3 / mean_v));
    if (n < 1) n = 1;

    std::vector<int> draw;
    for (int iter = 1; iter <= max_iter; ++iter)
    {
        iterations_used = iter;

        draw.clear();
        draw.reserve(static_cast<std::size_t>(n));
        for (long i = 0; i < n; ++i)
            draw.push_back(sampler.sample_pixels(next_uniform(rng_state)));

        // Per landslide, through the sampler's own scaling. A mean-area
        // shortcut would be wrong: V is non-linear in A.
        double vol = 0.0;
        for (int px : draw)
            vol += sampler.volume_of_pixels(px);

        const double rel = std::fabs(vol - target_volume_m3) / target_volume_m3;
        if (rel < best_rel) { best_rel = rel; best = draw; best_vol = vol; }

        if (rel < tol) { converged = true; break; }

        // Rescale N by how far off we are, then redraw. Both levers move: a bad
        // N is corrected by the ratio, a bad DRAW by resampling.
        //
        // When the ratio rounds back to the same N, RESAMPLE AT THAT N rather
        // than nudging it: N is already the best integer for this target and
        // the residual is draw scatter, which a fresh draw can fix and a forced
        // +/-1 step cannot.
        //
        // WHAT ACTUALLY LIMITS CONVERGENCE is not the search but the physics.
        // Volume is quantised by individual large landslides: with the Malamud
        // tail, one slide at A_max carries alpha*A_max^gamma, which for a 400
        // px threshold on a 10 m grid is ~197,000 m3. Against a 200,000 m3
        // budget that is 98% of the target in a SINGLE landslide, so no choice
        // of N can land within 2%. Tolerance is reachable only when the budget
        // is large compared with the biggest permitted landslide -- roughly
        // V_target >> alpha*A_max^gamma. Callers who need tight budget fidelity
        // must either raise the budget or lower A_max. The alternative --
        // trimming the last landslide to fit -- would buy exactness by
        // destroying the distribution at the large end, which is the one place
        // it must not be touched.
        if (!(vol > 0.0)) { ++n; continue; }
        long n_new = static_cast<long>(std::llround(
            static_cast<double>(n) * target_volume_m3 / vol));
        if (n_new < 1) n_new = 1;
        n = n_new;
    }

    sampled_volume_m3 = best_vol;
    return best;
}


LandslideEventResult trigger_landslide_event(
    lsdtt_raster_cpu&           dem,
    lsdtt_raster_cpu&           time_of_last_failure,
    lsdtt_raster_cpu&           excavation_depth,
    LandslidePlacer&            placer,
    const LandslideSizeSampler& sampler,
    const LandslideEventConfig& config)
{
    LandslideEventResult res;
    res.area_max_used_m2 =
        static_cast<double>(sampler.n_max_pixels()) * sampler.cell_area_m2();
    res.mean_volume_used_m3       = sampler.mean_volume_m3();
    res.truncated_fraction_number = sampler.truncated_fraction_number();
    res.truncated_fraction_volume = sampler.truncated_fraction_volume();
    res.rollover_unresolved       = sampler.rollover_unresolved();

    const lsdtt_raster_meta& meta = dem.meta();

    // Resolve the budget. volume-per-area is over the ACTIVE area, which is the
    // set of pixels that could host a landslide -- not the whole raster, since
    // a masked or channel-dominated DEM would otherwise get a budget it has no
    // room to spend.
    double target = config.target_volume_m3;
    if (!(target > 0.0) && config.volume_per_area_m3_per_m2 > 0.0)
    {
        const double active_area =
            static_cast<double>(placer.n_eligible_hosts(1)) * sampler.cell_area_m2();
        target = config.volume_per_area_m3_per_m2 * active_area;
    }
    res.requested_volume_m3 = target;
    if (!(target > 0.0)) return res;

    // --- 1. Budget solve: statistics only ---
    std::uint64_t rng = config.seed ? config.seed : 1;
    std::vector<int> sizes = solve_landslide_count_for_budget(
        sampler, target, config.budget_tol, config.budget_max_iter,
        rng, res.budget_converged, res.budget_iterations, res.sampled_volume_m3);

    res.n_requested = static_cast<int>(sizes.size());
    if (sizes.empty()) return res;

    // --- 2. Place, largest first ---
    // Small scars fragment the hillslopes and can leave no upslope set big
    // enough for a large slide that would have fitted on an empty landscape.
    // Descending order removes that ordering artefact.
    std::sort(sizes.begin(), sizes.end(), std::greater<int>());

    placer.begin_event();

    // Data arrays only: lsdtt_raster_cpu exposes no mutable mask accessor, so
    // both output rasters must arrive with a mask already covering the domain.
    // That is the right contract anyway -- an event does not change which cells
    // exist, only what is recorded in them.
    auto& dem_data = dem.get_data_array();
    auto& t_data   = time_of_last_failure.get_data_array();
    auto& exc_data = excavation_depth.get_data_array();

    LandslideScar scar;
    for (int n_px : sizes)
    {
        if (!placer.place_scar_seeded(n_px, rng, scar))
        {
            ++res.n_failed;
            continue;
        }

        const double area_m2 = static_cast<double>(n_px) * sampler.cell_area_m2();
        const double vol_m3  = sampler.volume_of_pixels(n_px);
        // Flat slab: uniform depth over the scar.
        const double depth_m = vol_m3 / area_m2;

        for (std::int64_t flat : scar.pixels)
        {
            const int r = static_cast<int>(flat / meta.n_cols);
            const int c = static_cast<int>(flat % meta.n_cols);

            dem_data(r, c) = static_cast<float>(
                static_cast<double>(dem_data(r, c)) - depth_m);

            // Accumulates, so repeated events on one pixel sum correctly.
            const float exc_ndv =
                static_cast<float>(excavation_depth.meta().no_data_value);
            const double prev = (exc_data(r, c) == exc_ndv)
                ? 0.0 : static_cast<double>(exc_data(r, c));
            exc_data(r, c) = static_cast<float>(prev + depth_m);

            // The annealing clock and the cosmogenic time-since-change field
            // are the same number, written once.
            t_data(r, c) = static_cast<float>(config.event_time_yr);
        }

        ++res.n_landslides;
        res.achieved_volume_m3 += vol_m3;
        res.achieved_area_m2   += area_m2;
        res.realised_pixels.push_back(n_px);
    }

    placer.finalise_diagnostics();
    res.median_host_fill_ratio = placer.diagnostics().median_host_fill_ratio;

    return res;
}

}  // namespace lsdtt
