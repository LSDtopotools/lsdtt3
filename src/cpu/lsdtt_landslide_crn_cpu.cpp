// lsdtt_landslide_crn_cpu.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
//
// See lsdtt_landslide_crn_cpu.hpp. No new CRN physics is defined here.

#include "lsdtt_crn_units.hpp"
#include "lsdtt_landslide_crn_cpu.hpp"

#include "lsdtt_cosmo_raster_cpu.hpp"
#include <tuple>
#include <sstream>

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

lsdtt_raster_cpu select_nuclide_raster(MultiNuclideRasters&& m,
                                       const std::string& nuclide)
{
    if (nuclide == "Be10") return std::move(m.Be10);
    if (nuclide == "Al26") return std::move(m.Al26);
    if (nuclide == "C14")  return std::move(m.C14);
    throw std::invalid_argument(
        "compute_landslide_step_change_concentration: unknown nuclide \""
        + nuclide + "\". Options: Be10, Al26, C14.");
}

// Effective erosion rate (g/cm2/yr) -> physical rate (m/yr).
// The inverse of lsdtt_crn_column's to_eff(). Kept as a named local so the
// call sites read as intent; the arithmetic lives in lsdtt_crn_units.hpp.
double eff_to_m_per_yr(double eff_g_cm2_yr, double rho_kg_m3)
{
    return eff_erate_to_m_yr(eff_g_cm2_yr, rho_kg_m3);
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


MultiNuclideRasters compute_landslide_surface_concentration_multi(
    const lsdtt_raster_cpu&   excavation_depth_m,
    const lsdtt_raster_cpu&   time_of_last_failure,
    const LandslideCrnConfig& config)
{
    // ONE KERNEL. This was a second, independent implementation of
    // landslide -> concentration until 2026-08-05: ~120 lines that built the
    // same four-row erosion history, cached on the same key and called the same
    // column, 260 lines away from the raster version and sharing nothing with
    // it. Nothing asserted the two agreed, and the density split of 2026-08-04
    // had to be applied to both by hand.
    //
    // It is now a THIN WRAPPER: a steady landscape is the step-change kernel
    // with equal rates, so the scalar inputs are expressed as the uniform
    // fields that kernel takes and the work happens in one place. The mapping:
    //
    //   production_scaling  -> a uniform production raster
    //   total_shielding     -> a uniform shielding raster
    //     (kept SEPARATE, as they must be: the column forms
    //      S_tot = scaling x shielding and solves the Vermeesch redistribution
    //      on the product. Merging them before this point is the mistake that
    //      cost this codebase 11.4 % once.)
    //   background rate     -> BOTH eff_erate_old and eff_erate_new (no step)
    //   time since change   -> 1e6 yr, i.e. "the landscape has been steady for
    //                          far longer than any scar age"
    //   flush window        -> 0, so every scar takes the scar treatment
    //
    // Pinned by unit_landslide_crn T17, which fixes the absolute values across
    // 384 combinations so the collapse cannot have moved any of them.
    if (excavation_depth_m.meta().n_rows != time_of_last_failure.meta().n_rows ||
        excavation_depth_m.meta().n_cols != time_of_last_failure.meta().n_cols)
        throw std::invalid_argument(
            "compute_landslide_surface_concentration: excavation depth and "
            "time-of-last-failure rasters differ in size");

    const lsdtt_raster_cpu t_since =
        time_since_failure_from_time_of_failure(time_of_last_failure,
                                                config.t_now_yr);

    const lsdtt_raster_cpu t_since_change =
        make_uniform_raster(excavation_depth_m, 1.0e6f, MaskedFill::Value);
    const lsdtt_raster_cpu prod =
        make_uniform_raster(excavation_depth_m,
                            static_cast<float>(config.production_scaling),
                            MaskedFill::Value);
    const lsdtt_raster_cpu shld =
        make_uniform_raster(excavation_depth_m,
                            static_cast<float>(config.total_shielding),
                            MaskedFill::Value);

    LandslideStepChangeConfig cfg;
    cfg.eff_erate_old_g_cm2_yr = config.background_eff_erate_g_cm2_yr;
    cfg.eff_erate_new_g_cm2_yr = config.background_eff_erate_g_cm2_yr;
    cfg.rock_density_kg_m3     = config.rock_density_kg_m3;
    cfg.scar_density_kg_m3     = config.scar_density_kg_m3;
    cfg.flush_time_threshold_yr = 0.0;
    cfg.muon_scheme            = config.muon_scheme;

    return compute_landslide_step_change_concentration_multi(
               excavation_depth_m, t_since, t_since_change, prod, shld, cfg)
           .concentration;
}


lsdtt_raster_cpu compute_landslide_surface_concentration(
    const lsdtt_raster_cpu&   excavation_depth_m,
    const lsdtt_raster_cpu&   time_of_last_failure,
    const LandslideCrnConfig& config)
{
    return select_nuclide_raster(
        compute_landslide_surface_concentration_multi(
            excavation_depth_m, time_of_last_failure, config),
        config.nuclide);
}


// -----------------------------------------------------------------------------
std::string LandslideCouplingResult::summary() const
{
    std::ostringstream os;
    os << "Landslide coupling\n"
       << "  failed pixels            : " << n_failed_pixels << "\n"
       << "    of which flushing      : " << n_flushing_pixels << "\n"
       << "    of which scars         : "
       << (n_failed_pixels - n_flushing_pixels) << "\n"
       << "  discarded (pre-wave)     : " << n_discarded_prewave << "\n"
       << "  distinct column solves   : " << n_solves << "\n";
    if (n_flushing_pixels > 0)
    {
        os << "  delivered by flushing    : " << total_delivered_kg_yr
           << " kg/yr\n";
    }
    if (n_discarded_prewave > 0)
    {
        os << "  NOTE: those failures pre-date the erosion-rate change at their\n"
              "  own pixel, so they cannot have been caused by it. They keep the\n"
              "  plain step-change value.\n";
    }
    return os.str();
}


std::string LandslideCouplingResultMulti::summary() const
{
    std::ostringstream os;
    os << "Landslide coupling\n"
       << "  failed pixels            : " << n_failed_pixels << "\n"
       << "    of which flushing      : " << n_flushing_pixels << "\n"
       << "    of which scars         : "
       << (n_failed_pixels - n_flushing_pixels) << "\n"
       << "  discarded (pre-wave)     : " << n_discarded_prewave << "\n"
       << "  distinct column solves   : " << n_solves << "\n";
    if (n_flushing_pixels > 0)
    {
        os << "  delivered by flushing    : " << total_delivered_kg_yr
           << " kg/yr\n";
    }
    if (n_discarded_prewave > 0)
    {
        os << "  NOTE: those failures pre-date the erosion-rate change at their\n"
              "  own pixel, so they cannot have been caused by it. They keep the\n"
              "  plain step-change value.\n";
    }
    return os.str();
}


// -----------------------------------------------------------------------------
lsdtt_raster_cpu time_since_failure_from_time_of_failure(
    const lsdtt_raster_cpu& time_of_last_failure,
    double                  t_report_yr)
{
    const lsdtt_raster_meta& meta = time_of_last_failure.meta();
    const auto& tf = time_of_last_failure.get_data_array();
    const auto& tm = time_of_last_failure.get_mask_array();
    const float ndv = static_cast<float>(meta.no_data_value);

    lsdtt_array2d_f  out(meta.n_rows, meta.n_cols, ndv);
    lsdtt_array2d_u8 omask(meta.n_rows, meta.n_cols, 0);

    for (int r = 0; r < meta.n_rows; ++r)
    {
        for (int c = 0; c < meta.n_cols; ++c)
        {
            if (!tm(r, c)) continue;
            omask(r, c) = 1;

            if (tf(r, c) == ndv || !std::isfinite(tf(r, c)))
            {
                out(r, c) = static_cast<float>(k_landslide_never_failed);
                continue;
            }
            // elapsed == 0 means the scar failed AT the reporting instant. That
            // IS a failure -- the surface is now the rock that was at depth d --
            // so it is reported as an elapsed time of zero, not as "never".
            // Only a failure in the FUTURE of the reporting time is "never".
            // See the header for the convention this replaced and why.
            const double elapsed = t_report_yr - static_cast<double>(tf(r, c));
            out(r, c) = (elapsed >= 0.0)
                ? static_cast<float>(elapsed)
                : static_cast<float>(k_landslide_never_failed);
        }
    }
    return lsdtt_raster_cpu(std::move(out), std::move(omask), meta);
}




namespace
{

// Gauss-Legendre nodes/weights on [-1, 1], generated by Newton iteration on the
// Legendre polynomial. Cheap for the orders used here and avoids a table that
// could be mistyped.
void gauss_legendre(int n, std::vector<double>& x, std::vector<double>& w)
{
    x.assign(n, 0.0);
    w.assign(n, 0.0);
    const double eps = 1.0e-14;
    for (int i = 0; i < n; ++i)
    {
        // Chebyshev starting guess, then Newton on P_n.
        double z = std::cos(M_PI * (static_cast<double>(i) + 0.75)
                            / (static_cast<double>(n) + 0.5));
        double pp = 0.0;
        for (int it = 0; it < 100; ++it)
        {
            double p0 = 1.0, p1 = 0.0;
            for (int j = 0; j < n; ++j)
            {
                const double p2 = p1;
                p1 = p0;
                p0 = ((2.0 * j + 1.0) * z * p1 - j * p2) / (j + 1.0);
            }
            pp = n * (z * p0 - p1) / (z * z - 1.0);
            const double dz = -p0 / pp;
            z += dz;
            if (std::fabs(dz) < eps) break;
        }
        x[i] = z;
        w[i] = 2.0 / ((1.0 - z * z) * pp * pp);
    }
}

}  // namespace


namespace
{

CrnColumnSample landslide_flush_concentration_sample(
    double d_eff_g_cm2,
    double eff_erate_old_g_cm2_yr,
    double eff_erate_new_g_cm2_yr,
    double time_since_change_at_failure_yr,
    double production_scaling,
    double total_shielding,
    const std::string& muon_scheme,
    double rock_density_kg_m3,
    int    n_nodes)
{
    lsdtt_crn_column col(rock_density_kg_m3);

    auto point = [&](double d) -> CrnColumnSample
    {
        return col.step_change(
            d, eff_erate_old_g_cm2_yr, eff_erate_new_g_cm2_yr,
            std::max(0.0, time_since_change_at_failure_yr),
            production_scaling, total_shielding, muon_scheme);
    };

    // A slab of zero thickness delivers nothing but is asked for at the limit;
    // return the surface value so the flush -> scar limit is continuous.
    if (!(d_eff_g_cm2 > 0.0)) return point(0.0);

    if (n_nodes < 2) n_nodes = 2;
    std::vector<double> x, w;
    gauss_legendre(n_nodes, x, w);

    // Map [-1,1] -> [0, d_eff]; the average is the integral over the thickness.
    const double half = 0.5 * d_eff_g_cm2;
    CrnColumnSample acc;
    acc.conc_10Be = 0.0;
    acc.conc_26Al = 0.0;
    acc.conc_14C = 0.0;
    for (int i = 0; i < n_nodes; ++i)
    {
        const CrnColumnSample s = point(half * (x[i] + 1.0));
        acc.conc_10Be += w[i] * s.conc_10Be;
        acc.conc_26Al += w[i] * s.conc_26Al;
        acc.conc_14C += w[i] * s.conc_14C;
    }
    acc.conc_10Be *= 0.5;
    acc.conc_26Al *= 0.5;
    acc.conc_14C *= 0.5;
    return acc;   // (half/d_eff) * acc = acc/2
}

} // namespace


double landslide_flush_concentration(
    double d_eff_g_cm2,
    double eff_erate_old_g_cm2_yr,
    double eff_erate_new_g_cm2_yr,
    double time_since_change_at_failure_yr,
    double production_scaling,
    double total_shielding,
    const std::string& muon_scheme,
    const std::string& nuclide,
    double rock_density_kg_m3,
    int    n_nodes)
{
    return pick_nuclide(
        landslide_flush_concentration_sample(
            d_eff_g_cm2, eff_erate_old_g_cm2_yr, eff_erate_new_g_cm2_yr,
            time_since_change_at_failure_yr, production_scaling,
            total_shielding, muon_scheme, rock_density_kg_m3, n_nodes),
        nuclide);
}


// -----------------------------------------------------------------------------
LandslideCouplingResultMulti compute_landslide_step_change_concentration_multi(
    const lsdtt_raster_cpu& excavation_depth_m,
    const lsdtt_raster_cpu& time_since_failure_yr,
    const lsdtt_raster_cpu& time_since_change_yr,
    const lsdtt_raster_cpu& production_raster,
    const lsdtt_raster_cpu& topo_shield,
    const LandslideStepChangeConfig& config,
    const lsdtt_raster_cpu* erosion_rate_g_cm2_yr)
{
    const lsdtt_raster_meta& meta = production_raster.meta();
    auto same_grid = [&](const lsdtt_raster_cpu& other, const char* what)
    {
        if (other.meta().n_rows != meta.n_rows ||
            other.meta().n_cols != meta.n_cols)
            throw std::invalid_argument(
                std::string("compute_landslide_step_change_concentration: ")
                + what + " does not match the production raster grid");
    };
    same_grid(excavation_depth_m,   "excavation depth");
    same_grid(time_since_failure_yr, "time since failure");
    same_grid(time_since_change_yr,  "time since change");
    same_grid(topo_shield,           "topographic shielding");
    if (erosion_rate_g_cm2_yr) same_grid(*erosion_rate_g_cm2_yr, "erosion rate");

    const bool do_flush = config.flush_time_threshold_yr > 0.0;
    // See the config comment: the selector and the rate denominator are
    // different quantities and only coincide when the window is finite.
    const double delivery_window_yr =
        (config.flush_delivery_window_yr > 0.0) ? config.flush_delivery_window_yr
                                                : config.flush_time_threshold_yr;

    // Unfailed pixels come from the ORDINARY step-change kernel, called on the
    // whole raster. Doing it this way is what makes "no landslides anywhere"
    // bitwise identical to a plain step-change run: there is no second
    // implementation of the step change here that could drift from the first.
    LandslideCouplingResultMulti result;
    result.concentration = calculate_crn_concentration_raster_step_change_multi(
        config.muon_scheme,
        time_since_change_yr, production_raster, topo_shield,
        config.eff_erate_old_g_cm2_yr, config.eff_erate_new_g_cm2_yr);

    auto& out10 = result.concentration.Be10.get_data_array();
    auto& out26 = result.concentration.Al26.get_data_array();
    auto& out14 = result.concentration.C14.get_data_array();
    const auto& omask = result.concentration.Be10.get_mask_array();

    const auto& depth = excavation_depth_m.get_data_array();
    const auto& dmask = excavation_depth_m.get_mask_array();
    const auto& tsf   = time_since_failure_yr.get_data_array();
    const auto& tsfm  = time_since_failure_yr.get_mask_array();
    const auto& tsc   = time_since_change_yr.get_data_array();
    const auto& tscm  = time_since_change_yr.get_mask_array();
    const auto& prod  = production_raster.get_data_array();
    const auto& shld  = topo_shield.get_data_array();

    const float dndv = static_cast<float>(excavation_depth_m.meta().no_data_value);

    const double rho_rock = config.rock_density_kg_m3;
    const double eps_old_m_yr = eff_to_m_per_yr(config.eff_erate_old_g_cm2_yr, rho_rock);
    const double eps_new_m_yr = eff_to_m_per_yr(config.eff_erate_new_g_cm2_yr, rho_rock);

    // One solve per distinct (depth, t_since_failure, t_since_change, site
    // factor). Every pixel in a scar shares the first three and, over a few
    // hundred metres, very nearly the fourth -- so the site factor is quantised
    // rather than used raw, or the cache would never hit.
    std::map<std::tuple<long long, long long, long long, long long>, CrnColumnSample> cache;

    for (int r = 0; r < meta.n_rows; ++r)
    {
        for (int c = 0; c < meta.n_cols; ++c)
        {
            if (!omask(r, c)) continue;
            if (!dmask(r, c) || !tsfm(r, c) || !tscm(r, c)) continue;

            const double d_m = static_cast<double>(depth(r, c));
            const double t_f = static_cast<double>(tsf(r, c));
            const double t_c = static_cast<double>(tsc(r, c));

            // t_f == 0 is a scar that failed AT the reporting instant, and it
            // counts. "Never failed" is a NEGATIVE sentinel
            // (k_landslide_never_failed), so >= 0 excludes it correctly.
            const bool failed = (t_f >= 0.0) && (depth(r, c) != dndv) && (d_m > 0.0);
            if (!failed) continue;

            // A failure older than the incision wave at this pixel cannot have
            // been caused by the change. Keep the plain step-change value.
            if (t_f > t_c)
            {
                ++result.n_discarded_prewave;
                continue;
            }
            ++result.n_failed_pixels;

            const double site = static_cast<double>(prod(r, c))
                              * static_cast<double>(shld(r, c));

            // FLUSH: the scar is young enough that its debris is still in the
            // channel, so what the sample contains is the whole excavated slab
            // -- the depth average of the PRE-FAILURE profile -- not the
            // surface of the scar.
            if (do_flush && t_f <= config.flush_time_threshold_yr)
            {
                ++result.n_flushing_pixels;
                const double d_eff =
                    metres_to_g_cm2(d_m, config.scar_density_kg_m3);

                const std::tuple<long long, long long, long long, long long>
                    fkey{static_cast<long long>(std::llround(d_m * 1000.0)),
                         static_cast<long long>(std::llround(t_f)),
                         static_cast<long long>(std::llround(t_c)),
                         // negative site key keeps flush and scar solves apart
                         -static_cast<long long>(std::llround(site * 1.0e6))};

                auto fit = cache.find(fkey);
                if (fit == cache.end())
                {
                    const CrnColumnSample v = landslide_flush_concentration_sample(
                        d_eff, config.eff_erate_old_g_cm2_yr,
                        config.eff_erate_new_g_cm2_yr,
                        t_c - t_f, static_cast<double>(prod(r, c)),
                        static_cast<double>(shld(r, c)), config.muon_scheme,
                        config.rock_density_kg_m3,
                        config.flush_quadrature_nodes);
                    fit = cache.emplace(fkey, v).first;
                    ++result.n_solves;
                }
                out10(r, c) = static_cast<float>(fit->second.conc_10Be);
                out26(r, c) = static_cast<float>(fit->second.conc_26Al);
                out14(r, c) = static_cast<float>(fit->second.conc_14C);
                continue;
            }

            const std::tuple<long long, long long, long long, long long> key{
                static_cast<long long>(std::llround(d_m * 1000.0)),
                static_cast<long long>(std::llround(t_f)),
                static_cast<long long>(std::llround(t_c)),
                static_cast<long long>(std::llround(site * 1.0e6))};

            auto it = cache.find(key);
            if (it == cache.end())
            {
                const double d_rock_equiv =
                    d_m * config.scar_density_kg_m3 / rho_rock;

                // t_f == 0: the scar failed AT the reporting instant, so there
                // is no re-accumulation yet and the surface is simply the
                // PRE-FAILURE column sampled at the scar depth -- inheritance.
                //
                // This needs its own branch rather than falling out of the
                // four-row history below, because that history would have a
                // ZERO-LENGTH most-recent segment. Measured: the degenerate
                // form gives 362,047 atoms/g where the correct answer is
                // 237,058, a factor of 1.53. So the pre-failure column is
                // evaluated directly, at depth, with no removal row.
                if (t_f == 0.0)
                {
                    // The pre-failure profile at the scar depth, via the same
                    // steady_state() call compute_landslide_surface_concentration
                    // makes -- which is the reference for this case.
                    //
                    // scaling and shielding are handed over as (site, 1.0): the
                    // column forms S_tot = scaling * shielding and solves the
                    // Vermeesch redistribution on the product, so passing the
                    // already-multiplied site value is exact, not an
                    // approximation.
                    //
                    // APPROXIMATION THAT IS REAL: under a genuine step change
                    // this uses the NEW rate's steady state rather than the
                    // mid-transient profile, so it is exact only where the
                    // incision wave passed long ago (t_c >> Gamma/eps). Good
                    // enough at the instant of failure, and it is the same
                    // assumption the scalar path has always made.
                    lsdtt_crn_column pre(rho_rock);
                    const CrnColumnSample sample = pre.steady_state(
                        metres_to_g_cm2(d_m, config.scar_density_kg_m3),
                        config.eff_erate_new_g_cm2_yr,
                        site, 1.0, config.muon_scheme);
                    const CrnColumnSample v = sample;
                    it = cache.emplace(key, v).first;
                    ++result.n_solves;
                    out10(r, c) = static_cast<float>(it->second.conc_10Be);
                    out26(r, c) = static_cast<float>(it->second.conc_26Al);
                    out14(r, c) = static_cast<float>(it->second.conc_14C);
                    continue;
                }

                // Four rows, most recent first, yr BP increasing with index.
                const std::vector<double> t_yr   {t_f, t_f, t_c, -9999.0};
                const std::vector<double> rate   {eps_new_m_yr, 0.0,
                                                  eps_new_m_yr, eps_old_m_yr};
                const std::vector<double> removal{0.0, d_rock_equiv, 0.0, 0.0};

                lsdtt_crn_column col(rho_rock);
                col.set_erosion_time_series(t_yr, rate, removal);

                const std::vector<CrnColumnSample> prof = col.run_transient(
                    std::vector<double>{0.0}, site, config.muon_scheme);

                CrnColumnSample v;
                if (prof.empty())
                {
                    v.conc_10Be = static_cast<double>(out10(r, c));
                    v.conc_26Al = static_cast<double>(out26(r, c));
                    v.conc_14C = static_cast<double>(out14(r, c));
                }
                else
                {
                    v = prof[0];
                }
                it = cache.emplace(key, v).first;
                ++result.n_solves;
            }

            out10(r, c) = static_cast<float>(it->second.conc_10Be);
            out26(r, c) = static_cast<float>(it->second.conc_26Al);
            out14(r, c) = static_cast<float>(it->second.conc_14C);
        }
    }

    // -------- delivery mass ------------------------------------------------
    // The weight a downstream mixing calculation should use. A flushing pixel
    // delivers its whole excavated slab over the flush window; everything else
    // delivers its background erosion. Both in kg m^-2 yr^-1, so they are
    // commensurable -- which is the point of not reusing the erosion-rate
    // raster as a weight.
    if (erosion_rate_g_cm2_yr)
    {
        const auto& er   = erosion_rate_g_cm2_yr->get_data_array();
        const auto& erm  = erosion_rate_g_cm2_yr->get_mask_array();
        const double cell_area_m2 = meta.cell_size * meta.cell_size;

        lsdtt_array2d_f  dm(meta.n_rows, meta.n_cols,
                            static_cast<float>(meta.no_data_value));
        lsdtt_array2d_u8 dmk(meta.n_rows, meta.n_cols, 0);

        for (int r = 0; r < meta.n_rows; ++r)
        {
            for (int c = 0; c < meta.n_cols; ++c)
            {
                if (!erm(r, c)) continue;

                double mass =
                    eff_erate_to_kg_m2_yr(static_cast<double>(er(r, c)));

                if (do_flush && dmask(r, c) && tsfm(r, c) && tscm(r, c))
                {
                    const double d_m = static_cast<double>(depth(r, c));
                    const double t_f = static_cast<double>(tsf(r, c));
                    const double t_c = static_cast<double>(tsc(r, c));
                    if (t_f >= 0.0 && depth(r, c) != dndv && d_m > 0.0 &&
                        t_f <= t_c && t_f <= config.flush_time_threshold_yr)
                    {
                        mass = config.scar_density_kg_m3 * d_m
                             / delivery_window_yr;
                        result.total_delivered_kg_yr += mass * cell_area_m2;
                    }
                }
                dm(r, c)  = static_cast<float>(mass);
                dmk(r, c) = 1;
            }
        }
        result.delivery_mass_kg_m2_yr =
            lsdtt_raster_cpu(std::move(dm), std::move(dmk), meta);
    }

    return result;
}


LandslideCouplingResult compute_landslide_step_change_concentration(
    const lsdtt_raster_cpu& excavation_depth_m,
    const lsdtt_raster_cpu& time_since_failure_yr,
    const lsdtt_raster_cpu& time_since_change_yr,
    const lsdtt_raster_cpu& production_raster,
    const lsdtt_raster_cpu& topo_shield,
    const LandslideStepChangeConfig& config,
    const lsdtt_raster_cpu* erosion_rate_g_cm2_yr)
{
    LandslideCouplingResultMulti mres =
        compute_landslide_step_change_concentration_multi(
            excavation_depth_m, time_since_failure_yr, time_since_change_yr,
            production_raster, topo_shield, config, erosion_rate_g_cm2_yr);

    LandslideCouplingResult out;
    out.concentration = select_nuclide_raster(std::move(mres.concentration),
                                              config.nuclide);
    out.delivery_mass_kg_m2_yr = std::move(mres.delivery_mass_kg_m2_yr);
    out.n_failed_pixels = mres.n_failed_pixels;
    out.n_discarded_prewave = mres.n_discarded_prewave;
    out.n_flushing_pixels = mres.n_flushing_pixels;
    out.n_solves = mres.n_solves;
    out.total_delivered_kg_yr = mres.total_delivered_kg_yr;
    return out;
}

}  // namespace lsdtt
