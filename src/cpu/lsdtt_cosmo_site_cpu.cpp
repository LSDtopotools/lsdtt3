// lsdtt_cosmo_site_cpu.cpp -- see the header for the design and its reasons.

#include "lsdtt_cosmo_site_cpu.hpp"

#include <iostream>
#include <sstream>
#include <stdexcept>

#include "../core/lsdtt_data_loader.hpp"
#include "../core/lsdtt_exe_path.hpp"

namespace lsdtt
{

namespace
{

// The messages below are reproduced verbatim from the two CLIs they were
// extracted from, minus the "FATAL: " prefix the caller adds. They were already
// identical in both, which is what made the unification clean.

AtmSource resolve_atmosphere(const std::string& name)
{
    if (name == "NCEP")     return AtmSource::NCEP;
    if (name == "standard") return AtmSource::Standard;
    if (name == "ERA40")    return AtmSource::ERA40;
    throw std::runtime_error(
        "unknown atmosphere \"" + name +
        "\". Allowed values: ERA40, NCEP, standard.");
}

ScalingScheme resolve_scaling_scheme(const std::string& name)
{
    if (name == "lsdn")  return ScalingScheme::LSDn;
    if (name == "stone") return ScalingScheme::Stone;
    throw std::runtime_error(
        "unknown scaling_scheme \"" + name + "\". Allowed values: stone, lsdn.");
}

bool is_emulator_scheme(const std::string& muon_scheme)
{
    return muon_scheme == "cronusEmulatorSt" ||
           muon_scheme == "cronusEmulatorLSDn";
}

}  // namespace


// -----------------------------------------------------------------------------
void validate_cosmo_site_config(const CosmoSiteConfig& cfg, CosmoSiteCaps caps)
{
    if (!has_cap(caps, CosmoSiteCaps::Production)) return;

    (void)resolve_atmosphere(cfg.atmosphere);
    const ScalingScheme scheme = resolve_scaling_scheme(cfg.scaling_scheme);

    if (scheme == ScalingScheme::LSDn &&
        cfg.nuclide != "Be10" && cfg.nuclide != "Al26" && cfg.nuclide != "C14")
    {
        throw std::runtime_error(
            "scaling_scheme: lsdn supports nuclide Be10, Al26 or C14 only (got \"" +
            cfg.nuclide + "\").");
    }
}


// -----------------------------------------------------------------------------
lsdtt_raster_cpu load_raster_or_uniform(const lsdtt_raster_cpu& reference,
                                        const std::string&      prefix,
                                        const std::string&      read_path,
                                        double                  uniform_value)
{
    if (prefix != "NULL")
    {
        const std::string path = read_path + prefix + ".tif";
        std::cout << "  Loading raster: " << path << "\n";
        lsdtt_data_loader ld;
        ld.open_raster(path);
        return lsdtt_raster_cpu(ld);
    }
    return make_uniform_raster(reference, static_cast<float>(uniform_value),
                               MaskedFill::Value);
}


// -----------------------------------------------------------------------------
lsdtt_raster_cpu CosmoSite::recompute_production(
    const lsdtt_raster_cpu& dem,
    double                  age_yr,
    const std::vector<std::int64_t>* eval_only) const
{
    return calculate_production_raster(dem, crn_nc_, age_yr, atm_source_,
                                       scaling_scheme_, nuclide_, eval_only);
}


// -----------------------------------------------------------------------------
std::string CosmoSite::summary() const
{
    std::ostringstream os;
    os << "cosmogenic site: " << (uses_lsdn() ? "LSDn" : "Stone 2000")
       << " scaling, " << atm_report_ << " atmosphere";
    if (production_from_disk_) os << ", production LOADED from disk";
    if (needs_pressure_)       os << ", pressure raster built";
    return os.str();
}


// -----------------------------------------------------------------------------
CosmoSite build_cosmo_site(const lsdtt_raster_cpu& dem,
                           const CosmoSiteConfig&  cfg,
                           CosmoSiteCaps           caps)
{
    CosmoSite site;
    site.nuclide_ = cfg.nuclide;

    if (has_cap(caps, CosmoSiteCaps::Production))
    {
        // Validate again rather than trusting the caller to have called
        // validate_cosmo_site_config(). It is cheap, and a component that
        // depends on being used correctly is a component that will be used
        // incorrectly.
        validate_cosmo_site_config(cfg, caps);

        site.atm_source_     = resolve_atmosphere(cfg.atmosphere);
        site.scaling_scheme_ = resolve_scaling_scheme(cfg.scaling_scheme);

        // Empty path -> the nc bundled beside the executable. Still empty after
        // that -> standard atmosphere, which is a supported configuration.
        site.crn_nc_path_ = cfg.crn_nc_file.empty()
                          ? default_crn_nc_path()
                          : cfg.read_path + cfg.crn_nc_file;
        if (!site.crn_nc_path_.empty())
        {
            std::cout << "Loading CRN NC data: " << site.crn_nc_path_ << "\n";
            site.crn_nc_ = lsdtt_crn_parameters::load_nc_data(site.crn_nc_path_);
        }

        if (site.uses_lsdn() && !site.crn_nc_.lsdn_loaded)
        {
            throw std::runtime_error(
                "scaling_scheme: lsdn requires the LSDn scaling grids in "
                "crn_parameters.nc, but none were found in \"" +
                (site.crn_nc_path_.empty()
                     ? std::string("<no crn_parameters.nc located>")
                     : site.crn_nc_path_) +
                "\".\nSet crn_nc_fname to an nc that carries the LSDn grids, or "
                "use scaling_scheme: stone.");
        }

        // What the atmosphere ACTUALLY resolved to after the documented
        // absent-grid fallback. Reported rather than assumed, because a run
        // that silently used `standard` where ERA40 was asked for biases every
        // pressure in the DEM.
        site.atm_report_ =
            (site.atm_source_ == AtmSource::NCEP  && site.crn_nc_.ncep_loaded) ? "NCEP"
          : (site.atm_source_ == AtmSource::ERA40 && site.crn_nc_.loaded)      ? "ERA40"
                                                                               : "standard";

        // Pressure raster: mandatory for the pressure-dependent emulator muon
        // schemes, and built before production so production can reuse it
        // rather than deriving the same values a second time.
        site.needs_pressure_ = is_emulator_scheme(cfg.muon_scheme);
        const bool build_pressure = site.needs_pressure_
                                 || has_cap(caps, CosmoSiteCaps::PressureAlways);
        if (build_pressure)
        {
            std::cout << "Calculating pressure raster ("
                      << site.atm_report_ << " atmosphere)...\n";
            try
            {
                site.pressure_ = calculate_pressure_raster(dem, site.crn_nc_,
                                                           site.atm_source_);
                site.has_pressure_ = true;
            }
            catch (const std::exception& e)
            {
                // Fatal only where something DEPENDS on it. Where it is a
                // reported value it degrades to nan with a warning, which is
                // what lsdtt-cosmo-tool already did.
                if (site.needs_pressure_)
                    throw std::runtime_error(
                        std::string("pressure raster failed: ") + e.what());
                std::cerr << "WARNING: pressure raster failed: " << e.what()
                          << "\n         effective_pressure_hpa will be "
                             "reported as nan.\n";
            }
            if (site.needs_pressure_ && cfg.muon_scheme == "cronusEmulatorLSDn")
            {
                site.lsdn_slhl_reference_ =
                    lsdtt_crn_parameters::cronus_emulator_lsdn_slhl_reference(
                        site.crn_nc_, cfg.nuclide);
            }
        }

        const bool from_disk = has_cap(caps, CosmoSiteCaps::ProductionFromDisk)
                            && cfg.load_production_raster
                            && cfg.production_raster_prefix != "NULL";
        if (from_disk)
        {
            const std::string prod_path =
                cfg.read_path + cfg.production_raster_prefix + ".tif";
            std::cout << "Loading production raster: " << prod_path << "\n";
            lsdtt_data_loader prod_loader;
            prod_loader.open_raster(prod_path);
            site.production_ = lsdtt_raster_cpu(prod_loader);
            site.production_from_disk_ = true;
        }
        else
        {
            std::cout << "Computing " << (site.uses_lsdn() ? "LSDn" : "Stone 2000")
                      << " production raster...\n";
            std::cout << "  Atmosphere: " << site.atm_report_ << "\n";
            site.production_ = site.recompute_production(dem, cfg.sample_age_yr,
                                                         nullptr);
            std::cout << "  Production raster computed.\n";
        }
    }

    // Shielding. With ShieldingRasters a prefix may point at a file; without
    // it, only the scalars exist and the raster branch is unreachable -- which
    // is the capability statement, not an optimisation.
    const bool rasters_ok = has_cap(caps, CosmoSiteCaps::ShieldingRasters);
    const std::string none = "NULL";
    site.topo_shield_ = load_raster_or_uniform(
        dem, rasters_ok ? cfg.topographic_shielding_raster_prefix : none,
        cfg.read_path, cfg.topo_shield_value);
    site.snow_shield_ = load_raster_or_uniform(
        dem, rasters_ok ? cfg.snow_shielding_raster_prefix : none,
        cfg.read_path, cfg.snow_shield_eff_thickness);
    site.self_shield_ = load_raster_or_uniform(
        dem, rasters_ok ? cfg.self_shielding_raster_prefix : none,
        cfg.read_path, cfg.self_shield_eff_thickness);

    return site;
}

}  // namespace lsdtt
