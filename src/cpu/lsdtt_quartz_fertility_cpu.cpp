// lsdtt_quartz_fertility_cpu.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// See the header for why quartz fertility is a correctness issue rather than a
// refinement, and for why the ice mask is a helper rather than a module.
// -----------------------------------------------------------------------------

#include "lsdtt_crn_units.hpp"
#include "lsdtt_quartz_fertility_cpu.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <sstream>
#include <stdexcept>

namespace lsdtt
{

namespace
{

void require_same_grid(const lsdtt_raster_meta& a, const lsdtt_raster_meta& b,
                       const char* what)
{
    if (a.n_rows != b.n_rows || a.n_cols != b.n_cols)
        throw std::invalid_argument(
            std::string("lsdtt_quartz_fertility_cpu: ") + what
            + " does not match the reference grid");
}

}  // namespace


// -----------------------------------------------------------------------------
double ice_effective_depth_g_cm2(double thickness_m, double ice_density_kg_m3)
{
    return metres_to_g_cm2(thickness_m, ice_density_kg_m3);
}


// -----------------------------------------------------------------------------
std::string QuartzWeightSummary::summary() const
{
    std::ostringstream os;
    os << "Quartz mixing weights\n"
       << "  valid pixels             : " << n_valid << "\n"
       << "  mean fertility           : " << mean_fertility << "\n"
       << "  quartz-free pixels       : " << n_zero_fertility;
    if (n_valid > 0)
        os << "  (" << (100.0 * static_cast<double>(n_zero_fertility)
                        / static_cast<double>(n_valid)) << " %)";
    os << "\n"
       << "  pixels contributing none : " << n_zero_weight << "\n"
       << "  total quartz delivered   : " << total_quartz_kg_yr << " kg/yr\n";

    if (yields_no_quartz())
    {
        os << "\n  WARNING: NO QUARTZ AT ALL is delivered anywhere in this "
              "domain.\n"
              "  A basin average taken from it is undefined -- there is nothing "
              "to\n"
              "  separate and nothing to measure. This is a meaningful answer, "
              "not a\n"
              "  failure: it is what a wholly quartz-free catchment gives. Check "
              "the\n"
              "  lithology raster and the ice mask before reading anything "
              "further.\n";
    }
    else if (n_zero_fertility > 0)
    {
        os << "\n  NOTE: " << n_zero_fertility << " pixel(s) deliver sediment "
              "but no quartz.\n"
              "  They are excluded from the mixture however fast they erode, "
              "which is\n"
              "  the whole point of weighting by fertility rather than by mass.\n";
    }
    return os.str();
}


// -----------------------------------------------------------------------------
lsdtt_raster_cpu make_quartz_fertility_raster(
    const lsdtt_raster_meta&     meta,
    const QuartzFertilityConfig& config,
    const lsdtt_raster_cpu*      lithology_raster,
    std::int64_t*                n_clamped)
{
    lsdtt_array2d_f  out(meta.n_rows, meta.n_cols,
                         static_cast<float>(meta.no_data_value));
    lsdtt_array2d_u8 mask(meta.n_rows, meta.n_cols, 0);

    std::int64_t clamped = 0;

    if (lithology_raster)
    {
        require_same_grid(meta, lithology_raster->meta(), "lithology raster");
        const auto& ld = lithology_raster->get_data_array();
        const auto& lm = lithology_raster->get_mask_array();

        for (int r = 0; r < meta.n_rows; ++r)
        {
            for (int c = 0; c < meta.n_cols; ++c)
            {
                if (!lm(r, c)) continue;
                double f = static_cast<double>(ld(r, c));
                if (!std::isfinite(f)) continue;
                if (f < 0.0 || f > 1.0)
                {
                    ++clamped;
                    f = std::clamp(f, 0.0, 1.0);
                }
                out(r, c)  = static_cast<float>(f);
                mask(r, c) = 1;
            }
        }
    }
    else
    {
        double f = config.uniform_fertility;
        if (f < 0.0 || f > 1.0)
        {
            clamped = static_cast<std::int64_t>(meta.n_rows) * meta.n_cols;
            f = std::clamp(f, 0.0, 1.0);
        }
        for (int r = 0; r < meta.n_rows; ++r)
            for (int c = 0; c < meta.n_cols; ++c)
            {
                out(r, c)  = static_cast<float>(f);
                mask(r, c) = 1;
            }
    }

    if (n_clamped) *n_clamped = clamped;
    return lsdtt_raster_cpu(std::move(out), std::move(mask), meta);
}


// -----------------------------------------------------------------------------
double delivered_quartz_fertility(const QuartzFertilityConfig& config,
                                  double surface_fertility,
                                  double depth_m)
{
    if (!config.depth_dependent()) return surface_fertility;
    if (!(depth_m > 0.0))          return config.regolith_fertility;

    const double d_reg = std::min(depth_m, config.regolith_depth_m);
    const double d_bed = std::max(0.0, depth_m - config.regolith_depth_m);
    return (config.regolith_fertility * d_reg
          + config.bedrock_fertility  * d_bed) / depth_m;
}


// -----------------------------------------------------------------------------
lsdtt_raster_cpu make_quartz_weight_raster(
    const lsdtt_raster_cpu&      mass_flux_kg_m2_yr,
    const lsdtt_raster_cpu&      fertility,
    const QuartzFertilityConfig& config,
    const lsdtt_raster_cpu*      excavation_depth_m,
    QuartzWeightSummary*         summary)
{
    const lsdtt_raster_meta& meta = mass_flux_kg_m2_yr.meta();
    require_same_grid(meta, fertility.meta(), "fertility raster");
    if (excavation_depth_m)
        require_same_grid(meta, excavation_depth_m->meta(), "excavation raster");

    const auto& md = mass_flux_kg_m2_yr.get_data_array();
    const auto& mm = mass_flux_kg_m2_yr.get_mask_array();
    const auto& fd = fertility.get_data_array();
    const auto& fm = fertility.get_mask_array();

    const double cell_area_m2 = meta.cell_size * meta.cell_size;

    lsdtt_array2d_f  out(meta.n_rows, meta.n_cols,
                         static_cast<float>(meta.no_data_value));
    lsdtt_array2d_u8 omask(meta.n_rows, meta.n_cols, 0);

    QuartzWeightSummary s;
    double fert_sum = 0.0;

    for (int r = 0; r < meta.n_rows; ++r)
    {
        for (int c = 0; c < meta.n_cols; ++c)
        {
            if (!mm(r, c) || !fm(r, c)) continue;

            double f = static_cast<double>(fd(r, c));

            // A failed pixel delivers its excavated slab, not its surface, so
            // under the two-layer model it carries the slab's average
            // fertility. This is the only place the depth dependence acts.
            if (excavation_depth_m && config.depth_dependent())
            {
                const auto& ed = excavation_depth_m->get_data_array();
                const auto& em = excavation_depth_m->get_mask_array();
                const float endv =
                    static_cast<float>(excavation_depth_m->meta().no_data_value);
                if (em(r, c) && ed(r, c) != endv && ed(r, c) > 0.0f)
                    f = delivered_quartz_fertility(
                        config, f, static_cast<double>(ed(r, c)));
            }

            const double w = static_cast<double>(md(r, c)) * f;

            out(r, c)   = static_cast<float>(w);
            omask(r, c) = 1;

            ++s.n_valid;
            fert_sum += f;
            if (!(f > 0.0)) ++s.n_zero_fertility;
            if (!(w > 0.0)) ++s.n_zero_weight;
            s.total_quartz_kg_yr += w * cell_area_m2;
        }
    }

    if (s.n_valid > 0) s.mean_fertility = fert_sum / static_cast<double>(s.n_valid);
    if (summary) *summary = s;

    return lsdtt_raster_cpu(std::move(out), std::move(omask), meta);
}


// -----------------------------------------------------------------------------
std::int64_t apply_ice_mask(const lsdtt_raster_cpu&  ice_mask,
                            const lsdtt_raster_cpu*  ice_thickness_raster,
                            const IceMaskConfig&     config,
                            lsdtt_raster_cpu&        snow_shield_g_cm2,
                            lsdtt_raster_cpu&        erosion_rate_g_cm2_yr,
                            lsdtt_raster_cpu&        quartz_fertility)
{
    const lsdtt_raster_meta& meta = ice_mask.meta();
    const auto& id = ice_mask.get_data_array();
    const auto& im = ice_mask.get_mask_array();
    const float indv = static_cast<float>(meta.no_data_value);

    if (ice_thickness_raster)
        require_same_grid(meta, ice_thickness_raster->meta(), "ice thickness");

    auto usable = [&](lsdtt_raster_cpu& r) -> bool
    {
        if (r.meta().n_rows == 0) return false;
        require_same_grid(meta, r.meta(), "an output raster");
        return true;
    };

    const bool do_snow = usable(snow_shield_g_cm2);
    const bool do_ero  = usable(erosion_rate_g_cm2_yr);
    const bool do_fert = usable(quartz_fertility);

    std::int64_t n_ice = 0;

    for (int r = 0; r < meta.n_rows; ++r)
    {
        for (int c = 0; c < meta.n_cols; ++c)
        {
            if (!im(r, c) || id(r, c) == indv) continue;
            if (!(id(r, c) != 0.0f)) continue;      // zero means no ice
            ++n_ice;

            const double thickness = ice_thickness_raster
                ? static_cast<double>(
                      ice_thickness_raster->get_data_array()(r, c))
                : config.ice_thickness_m;

            // Values are written only where the target pixel is already
            // valid. lsdtt_raster_cpu exposes no mutable mask accessor, by
            // design: an ice mask says what the surface is made of, not which
            // cells exist, and inventing domain would be a different and much
            // larger claim. Ice falling on nodata is counted and reported
            // rather than silently creating data.
            if (do_snow && snow_shield_g_cm2.get_mask_array()(r, c))
            {
                snow_shield_g_cm2.get_data_array()(r, c) = static_cast<float>(
                    ice_effective_depth_g_cm2(thickness, config.ice_density_kg_m3));
            }
            if (do_ero && erosion_rate_g_cm2_yr.get_mask_array()(r, c))
            {
                erosion_rate_g_cm2_yr.get_data_array()(r, c) =
                    static_cast<float>(config.erate_under_ice_g_cm2_yr);
            }
            if (do_fert && quartz_fertility.get_mask_array()(r, c))
            {
                quartz_fertility.get_data_array()(r, c) =
                    static_cast<float>(config.quartz_fertility_under_ice);
            }
        }
    }
    return n_ice;
}

}  // namespace lsdtt
