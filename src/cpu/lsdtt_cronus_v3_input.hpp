// =============================================================================
// lsdtt_cronus_v3_input.hpp
//
// Emit a CRONUS v3 "erosion_input_v3" text block: the batch format the CRONUS
// v3 online calculator accepts by copy-and-paste into its input box.
//
// Why this lives in the library rather than in a script: a user comparing an
// lsdtt3 rate against the online calculator should not have to run a separate
// Python tool to build the request.  Both lsdtt-cosmo-tool and
// lsdtt-cosmo-column-tool write this file, so the format lives in exactly one
// place and the two tools cannot drift apart.
//
// Format, two lines per sample (the second repeats for each nuclide):
//
//   <id> <lat> <lon> <elev|pressure> <flag> <thick_cm> <rho> <shield> <eps> <yr> ;
//   <id> Be-10 quartz <conc_atoms_g> <sigma_atoms_g> <standard> ;
//
// <flag> is "std" when the fourth field is an elevation in metres and the
// calculator should derive pressure from the standard atmosphere, or "pre"
// when it is already a pressure in hPa.  Basin-averaged samples must use
// "pre": a catchment has no single site elevation, so lsdtt-cosmo-tool passes
// the basin-effective pressure with the basin-centroid latitude.
// =============================================================================

#ifndef LSDTT_CRONUS_V3_INPUT_HPP
#define LSDTT_CRONUS_V3_INPUT_HPP

#include <cstddef>
#include <string>
#include <vector>

namespace lsdtt
{

// One nuclide measurement belonging to a site.
struct CronusV3Nuclide
{
    // lsdtt3 spelling ("Be10", "Al26", "C14"); translated to the calculator's
    // own spelling ("Be-10", ...) on write.
    std::string nuclide;
    double      conc_atoms_g   = 0.0;
    // 1-sigma, atoms/g.  Written as-is: a non-positive value is emitted as 0
    // rather than replaced by a guess, so a file never carries an invented
    // uncertainty that could be mistaken for a measurement.  Callers should
    // warn when they have none.
    double      sigma_atoms_g  = 0.0;
    // AMS standardisation.  Empty selects the calculator's own default for the
    // nuclide (07KNSTD for Be-10, KNSTD for Al-26; C-14 takes no standard).
    std::string standardisation;
};

struct CronusV3Sample
{
    std::string name;                       // sanitised on write
    double      latitude_deg  = 0.0;
    double      longitude_deg = 0.0;

    // true  -> elevation_or_pressure is a pressure in hPa      (flag "pre")
    // false -> elevation_or_pressure is an elevation in metres (flag "std")
    bool        use_pressure  = true;
    double      elevation_or_pressure = 0.0;

    double      thickness_cm  = 0.0;        // surface sample
    double      density_g_cm3 = 2.65;
    double      shielding     = 1.0;        // 1 = none; see DiBiase (2018)
    double      erate_guess_g_cm2_yr = 0.0001;
    int         sample_year   = 2010;

    std::vector<CronusV3Nuclide> nuclides;
};

// Make a token safe for the calculator's parser: it splits on whitespace and
// terminates records on ';', so anything outside [A-Za-z0-9_-] becomes '_'.
// An empty or fully-stripped name yields the supplied fallback.
std::string cronus_v3_sanitise_id(const std::string& raw,
                                  const std::string& fallback);

// Render one sample's records.  Ends with a newline.  Samples carrying no
// nuclide produce an empty string, since a site record alone is not a valid
// request.
std::string format_cronus_v3_sample(const CronusV3Sample& sample);

// Write a whole batch.  Ids are sanitised and de-duplicated (a repeated name
// gains a _2, _3, ... suffix) because the calculator keys its results by id and
// silently merges collisions.  Returns the number of samples written; throws
// std::runtime_error if the file cannot be opened.
std::size_t write_cronus_v3_input_file(const std::string&                  path,
                                       const std::vector<CronusV3Sample>&  samples);

}  // namespace lsdtt

#endif  // LSDTT_CRONUS_V3_INPUT_HPP
