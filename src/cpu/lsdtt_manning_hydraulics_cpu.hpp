// lsdtt_manning_hydraulics_cpu.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html
//
// -----------------------------------------------------------------------------
// Manning's-equation bankfull width estimation 
//
// Given a reach-averaged Manning's n and a target bankfull discharge
// back-calculated at one reference cross section per reach (external,
// done in Python), solves for the water depth at EVERY cross section in
// the reach that reproduces that discharge, using that cross section's
// OWN slope (looked up from the channel-profile CSV, not averaged over
// the reach). Reports the bank offsets and true x/y locations where the
// water surface intersects the profile.
//
// IMPORTANT: reaches are bounded by DISTANCE FROM SOURCE
// (start_distance_m / end_distance_m), not by cross_section_id.
// cross_section_id is assigned during extraction by a row-major scan of
// the raster and has NO relationship to position along the channel.
// A sinuous channel can produce wildly non-monotonic cross_section_id
// values along its length. Only distance-from-source (computed via a
// single downstream walk in compute_channel_profile) is a valid basis
// for defining a contiguous stretch of channel.
//
// write_channel_profile_csv assumes a SINGLE mainstem channel (one
// source point) so distance-from-source can be computed as a plain
// downstream walk from that source, with no branch/tributary
// bookkeeping. Slope is sampled directly from the DEM via
// direct_polyfit at each cross section's centre pixel -- it does not
// depend on channel-network topology or cross-section point spacing.
// -----------------------------------------------------------------------------

#pragma once

#include <cstdint>
#include <string>
#include <unordered_map>
#include <vector>

#include "lsdtt_cross_sections_cpu.hpp"
#include "lsdtt_raster_cpu.hpp"
#include "../core/lsdtt_raster_meta.hpp"

namespace lsdtt
{

struct lsdtt_wetted_geometry
{
    double area_m2            = 0.0;
    double wetted_perimeter_m = 0.0;
    double hydraulic_radius_m = 0.0;
    bool   valid              = false;
};

struct lsdtt_bank_crossing
{
    bool   found = false;
    double s = 0.0;
    double x = 0.0;
    double y = 0.0;
};

// One row read from the user's reach-parameter CSV.
// NOTE: bounded by distance-from-source, NOT cross_id. Slope is not
// here either -- it is looked up per cross section from the
// channel-profile CSV.
struct lsdtt_manning_reach
{
    double start_distance_m     = 0.0;
    double end_distance_m       = 0.0;
    double manning_n            = 0.04;
    double target_discharge_m3s = 0.0;
};

// Per-cross-section slope, drainage area, and distance-from-source,
// loaded from the CSV written by write_channel_profile_csv. Keyed by
// cross_section_id. distance_by_cross_id is required for matching
// cross sections to reaches (see lsdtt_manning_reach); slope/area are
// carried through for the Manning solve and for QA output.
//
// A cross_id absent from slope_by_cross_id means slope was invalid at
// extraction time (slope_valid == 0 in the source CSV) -- the Manning
// solve flags such sections rather than guessing a value.
struct lsdtt_channel_profile_lookup
{
    std::unordered_map<long long, double> slope_by_cross_id;
    std::unordered_map<long long, double> drainage_area_by_cross_id;
    std::unordered_map<long long, double> distance_by_cross_id;
};

struct manning_solver_options
{
    double min_depth_m       = 0.01;
    double max_depth_m       = 15.0;
    double depth_tolerance_m = 0.005;
    int    max_iterations    = 60;
};

struct lsdtt_manning_width_record
{
    long long cross_section_id = -1;
    double centre_x = 0.0;
    double centre_y = 0.0;

    // Per-XS values looked up from the channel-profile CSV, carried
    // through to the output for traceability/QA.
    double distance_from_source_m = 0.0;
    double slope_m_per_m          = 0.0;
    double drainage_area_m2       = 0.0;

    double thalweg_s = 0.0;
    double thalweg_z = 0.0;

    double solved_depth_m            = 0.0;
    double water_surface_elevation_m = 0.0;

    double left_bank_s  = 0.0;
    double right_bank_s = 0.0;
    double left_bank_x  = 0.0;
    double left_bank_y  = 0.0;
    double right_bank_x = 0.0;
    double right_bank_y = 0.0;
    double width_m      = 0.0;

    double area_m2                   = 0.0;
    double wetted_perimeter_m        = 0.0;
    double hydraulic_radius_m        = 0.0;
    double discharge_at_solution_m3s = 0.0;

    bool        valid = false;
    std::string quality_flag = "invalid_profile";
};

// One row per cross section: distance from source, drainage area,
// thalweg elevation, and raster-sampled slope.
struct lsdtt_channel_profile_row
{
    long long cross_section_id       = -1;
    double    distance_from_source_m = 0.0;
    double    drainage_area_m2       = 0.0;
    double    thalweg_elevation_m    = 0.0;
    double    slope_m_per_m          = 0.0;
    bool      slope_valid            = false;
};

// Reads a reach-parameter CSV with header:
//   start_distance_m,end_distance_m,manning_n,target_discharge_m3s
// Bounds are along-channel distance from source (metres), NOT cross_id.
std::vector<lsdtt_manning_reach> read_manning_reaches_csv(
    const std::string& csv_path);

// Reads the channel-profile CSV written by write_channel_profile_csv
// (header: cross_id,distance_from_source_m,drainage_area_m2,
// thalweg_elevation_m,slope_m_per_m,slope_valid) into a lookup keyed by
// cross_id. cross_id values must match the cross_section_id numbering
// of the sections passed to compute_manning_widths_all_reaches.
lsdtt_channel_profile_lookup read_channel_profile_csv(
    const std::string& csv_path);

lsdtt_wetted_geometry compute_wetted_geometry(
    const std::vector<double>& s,
    const std::vector<double>& z,
    double water_surface_elevation);

double manning_discharge(
    double area_m2, double hydraulic_radius_m,
    double slope, double manning_n);

// Bisects for the depth above the thalweg that reproduces
// target_discharge_m3s at the given slope. Returns false (does not
// throw) if the target isn't bracketed within
// [opts.min_depth_m, opts.max_depth_m] -- this is an expected outcome
// at reach boundaries or unusual cross sections, not an error.
bool solve_bankfull_depth(
    const std::vector<double>& s,
    const std::vector<double>& z,
    double thalweg_elevation,
    double slope,
    double manning_n,
    double target_discharge_m3s,
    const manning_solver_options& opts,
    double& solved_depth_m);

// Runs the Manning solve over every cross section in `sections` whose
// distance-from-source (looked up from `profile`) falls within
// [reach.start_distance_m, reach.end_distance_m]. reach.manning_n and
// reach.target_discharge_m3s are held fixed across the reach; slope is
// looked up per cross section from `profile`, so solved depth varies
// with each cross section's own slope and geometry even though Q and n
// are constant. Cross sections absent from `profile` entirely (i.e. not
// present in the channel-profile CSV) are silently skipped -- they
// cannot be assigned to any reach without a distance value.
std::vector<lsdtt_manning_width_record> compute_manning_widths_for_reach(
    const std::vector<lsdtt_cross_section>& sections,
    const lsdtt_manning_reach&              reach,
    const lsdtt_channel_profile_lookup&     profile,
    const manning_solver_options&           opts = manning_solver_options{});

// Runs compute_manning_widths_for_reach over every reach and concatenates
// the results, sorted by distance_from_source_m.
std::vector<lsdtt_manning_width_record> compute_manning_widths_all_reaches(
    const std::vector<lsdtt_cross_section>&    sections,
    const std::vector<lsdtt_manning_reach>&    reaches,
    const lsdtt_channel_profile_lookup&        profile,
    const manning_solver_options&              opts = manning_solver_options{});

// One point per bank (left/right) per valid cross section, via
// lsdtt_data_loader. FlatGeobuf preferred; gpkg/csv also supported.
void write_manning_bank_points(
    const std::vector<lsdtt_manning_width_record>& results,
    const lsdtt_raster_meta&                        meta,
    const std::string&                              out_path,
    point_output_format                             format,
    bool                                             spatial_index,
    const std::string&                               layer_name);

// One row per cross section (valid and invalid), with quality_flag.
std::size_t write_manning_widths_csv(
    const std::vector<lsdtt_manning_width_record>& results,
    const std::string&                              out_path);

// -----------------------------------------------------------------------------
// compute_channel_profile / write_channel_profile_csv
//
// Distance from source, drainage area, thalweg elevation, and slope per
// cross section -- input to the Python reach-identification workflow AND
// to the Manning solve above (via read_channel_profile_csv). Independent
// of the Manning solve; can be run before reach parameters exist.
//
// REQUIRES the channel network `sections` was extracted from was built
// from a SINGLE source point (source_flat). Distance-from-source is
// computed as one downstream walk from source_flat; any cross section
// not on that path (i.e. the network had more than one source) will
// cause a std::runtime_error.
//
// slope_window_radius_metres is the direct_polyfit window radius used
// to sample slope at each cross section's centre pixel.
// -----------------------------------------------------------------------------
std::vector<lsdtt_channel_profile_row> compute_channel_profile(
    const std::vector<lsdtt_cross_section>& sections,
    const lsdtt_flow_base&                  flow,
    const lsdtt_raster_cpu&                 dem,
    float                                    slope_window_radius_metres,
    std::int64_t                             source_flat,
    int                                       n_threads = 0);

void write_channel_profile_csv(
    const std::vector<lsdtt_channel_profile_row>& rows,
    const std::string&                            out_path);

} // namespace lsdtt