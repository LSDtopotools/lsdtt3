// lsdtt_manning_hydraulics_cpu.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Martha Hawksworth
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// lsdtt_manning_hydraulics_cpu
//
// Manning's-equation bankfull width estimation.
//
// For each cross section in a reach:
//
//   1. Look up the reach's Manning's n and target discharge (reach-averaged,
//      back-calculated at one reference cross section, external to this file)
//   2. Look up that cross section's own slope from the channel profile
//   3. Solve for the water depth that reproduces the target discharge at
//      this slope, not a reach-averaged slope
//   4. Find where that water surface intersects the profile on each side
//   5. Report the bank offsets and true x/y crossing locations
//
// Reaches are bounded by distance from source (start_distance_m /
// end_distance_m), not by cross_section_id. cross_section_id is assigned
// during extraction by a row-major scan of the raster and has no
// relationship to position along the channel. Distance-from-source,
// computed via a single downstream walk in compute_channel_profile, is the
// only valid basis for a contiguous stretch of channel.
//
// Reach membership is HALF-OPEN: start_distance_m <= d < end_distance_m, so a
// cross section sitting exactly on a shared boundary belongs to exactly one
// reach. Reaches must not overlap.
//
// compute_channel_profile assumes a single mainstem channel (one source
// point), so distance-from-source is a plain downstream walk from that
// source with no branch/tributary bookkeeping. Slope is the along-channel
// least-squares slope of the pit-filled DEM over a window of the D8 channel
// path centred on each cross section (see compute_channel_profile).
//
// Units: distances and depths are metres, slope is m/m, discharge is m3/s,
// drainage area is m2 (CSV column names keep the short _m / _m2 suffixes
// because downstream Python workflows read them).
// -----------------------------------------------------------------------------

#pragma once

#include <cstddef>
#include <cstdint>
#include <string>
#include <unordered_map>
#include <vector>

#include "lsdtt_cross_sections_cpu.hpp"
#include "lsdtt_raster_cpu.hpp"
#include "../core/lsdtt_data_loader.hpp"
#include "../core/lsdtt_flow_base.hpp"
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
// Bounded by distance-from-source, NOT cross_id. Slope is not here either,
// it is looked up per cross section from the channel profile.
struct lsdtt_manning_reach
{
    double start_distance_m     = 0.0;
    double end_distance_m       = 0.0;
    double manning_n            = 0.04;
    double target_discharge_m3s = 0.0;
};

// Per-cross-section slope, drainage area, and distance-from-source, keyed by
// cross_section_id. distance_by_cross_id is required for matching cross
// sections to reaches.
//
// A cross_id absent from slope_by_cross_id means slope was invalid
// (slope_valid == 0), the Manning solve flags such sections rather than
// guessing a value.
struct lsdtt_channel_profile_lookup
{
    std::unordered_map<long long, double> slope_by_cross_id;
    std::unordered_map<long long, double> drainage_area_by_cross_id;
    std::unordered_map<long long, double> distance_by_cross_id;
};

struct lsdtt_manning_solver_options
{
    double min_depth_m       = 0.01;
    double max_depth_m       = 15.0;
    double depth_tolerance_m = 0.005;
    int    max_iterations    = 60;

    // The thalweg is the lowest valid sample within +/- this distance of the
    // section centre (offset 0), NOT the global minimum of the whole section.
    // A road ditch, neighbouring gully or second channel elsewhere in the
    // section must not be able to capture the solve. Set to roughly half the
    // channel width (the D8 centre pixel can sit that far from the thalweg).
    double thalweg_search_half_width_m = 25.0;

    // After solving, the discharge at the solved depth must match the target
    // to within this fraction, else the section is flagged
    // "spill_discontinuity". Connected wetted area jumps when the water
    // surface overtops a bank crest into a lower area, so the target Q can
    // fall inside that jump and bisection lands on the crest, not on Q.
    double max_discharge_relative_error = 0.05;
};

struct lsdtt_manning_width_record
{
    long long cross_section_id = -1;
    double centre_x = 0.0;
    double centre_y = 0.0;

    // Per-XS values looked up from the channel profile, carried through to
    // the output for traceability/QA.
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
// thalweg elevation, and along-channel slope.
struct lsdtt_channel_profile_row
{
    long long cross_section_id       = -1;
    double    distance_from_source_m = 0.0;
    double    drainage_area_m2       = 0.0;
    double    thalweg_elevation_m    = 0.0;
    double    slope_m_per_m          = 0.0;
    bool      slope_valid            = false;
};

// One row per reach in the template CSV. Coordinates are in the raster CRS
// (typically UTM, metres) for QGIS visualization. manning_n and
// target_discharge_m3s are left empty for the user to fill in.
struct lsdtt_reach_template_row
{
    int reach_index = 0;
    double start_distance_m = 0.0;
    double end_distance_m = 0.0;
    double start_easting_m = 0.0;
    double start_northing_m = 0.0;
    double end_easting_m = 0.0;
    double end_northing_m = 0.0;
};

// One row per reference cross-section candidate. Coordinates are in the
// raster CRS (typically UTM, metres) for QGIS visualization.
struct lsdtt_reference_xs_candidate
{
    int reach_index = 0;
    double reach_start_distance_m = 0.0;
    double reach_end_distance_m = 0.0;
    long long cross_id = -1;
    double distance_from_source_m = 0.0;
    double easting_m = 0.0;
    double northing_m = 0.0;
    double drainage_area_m2 = 0.0;
    double slope_m_per_m = 0.0;
    int slope_valid = 0;
};

// How many extracted sections the profile lookup can actually serve.
struct lsdtt_manning_profile_coverage
{
    std::size_t n_sections         = 0;
    std::size_t n_missing_distance = 0;  // silently unassignable to any reach
    std::size_t n_missing_slope    = 0;  // will be flagged, not solved
};

// Reads a reach-parameter CSV with header:
//   start_distance_m,end_distance_m,manning_n,target_discharge_m3s
// Bounds are along-channel distance from source (metres), NOT cross_id.
// Throws std::runtime_error (naming the line) on a short row or an
// unparseable number rather than skipping it.
[[nodiscard]] std::vector<lsdtt_manning_reach> read_manning_reaches_csv(
    const std::string& csv_path);

// Throws std::invalid_argument if the list is empty, any reach has
// start >= end, manning_n <= 0 or target discharge <= 0, or two reaches
// overlap. Gaps between reaches are allowed.
void validate_manning_reaches(const std::vector<lsdtt_manning_reach>& reaches);

// Reads the channel-profile CSV written by write_channel_profile_csv
// (header: cross_id,distance_from_source_m,drainage_area_m2,
// thalweg_elevation_m,slope_m_per_m,slope_valid). Throws on malformed rows.
[[nodiscard]] lsdtt_channel_profile_lookup read_channel_profile_csv(
    const std::string& csv_path);

// Builds the lookup directly from compute_channel_profile output, so a single
// invocation never has to round-trip through a CSV that could belong to a
// different run.
[[nodiscard]] lsdtt_channel_profile_lookup make_channel_profile_lookup(
    const std::vector<lsdtt_channel_profile_row>& rows);

[[nodiscard]] lsdtt_manning_profile_coverage compute_manning_profile_coverage(
    const std::vector<lsdtt_cross_section>& sections,
    const lsdtt_channel_profile_lookup&     profile);

[[nodiscard]] lsdtt_wetted_geometry compute_wetted_geometry(
    const std::vector<double>& s,
    const std::vector<double>& z,
    double water_surface_elevation);

// Wetted geometry of ONLY the water body connected to the thalweg at
// `thalweg_idx`: integrates outward each side until the first sample at or
// above the water surface, clipping the final segment at the interpolated
// bank crossing. Disconnected low spots (ditches, gullies, pools behind a
// bank crest) are never counted, so area, wetted perimeter and the bank
// crossings from the solve describe the same water. Invalid if the thalweg
// sample is not below the water surface.
[[nodiscard]] lsdtt_wetted_geometry compute_connected_wetted_geometry(
    const std::vector<double>& s,
    const std::vector<double>& z,
    std::size_t thalweg_idx,
    double water_surface_elevation);

[[nodiscard]] double manning_discharge(
    double area_m2, double hydraulic_radius_m,
    double slope, double manning_n);

// Bisects for the depth above the thalweg that reproduces
// target_discharge_m3s at the given slope, using only water connected to the
// channel. Returns false (does not throw) if the target isn't bracketed
// within [opts.min_depth_m, opts.max_depth_m] -- an expected outcome at reach
// boundaries or unusual cross sections, not an error.
[[nodiscard]] bool solve_connected_bankfull_depth(
    const std::vector<double>& s,
    const std::vector<double>& z,
    std::size_t thalweg_idx,
    double slope,
    double manning_n,
    double target_discharge_m3s,
    const lsdtt_manning_solver_options& opts,
    double& solved_depth_m);

// Runs the Manning solve over every cross section whose distance-from-source
// (looked up from `profile`) satisfies start <= d < end. reach.manning_n and
// reach.target_discharge_m3s are held fixed across the reach; slope is looked
// up per cross section. Cross sections absent from `profile` are skipped
// (use compute_manning_profile_coverage to count them).
[[nodiscard]] std::vector<lsdtt_manning_width_record>
compute_manning_widths_for_reach(
    const std::vector<lsdtt_cross_section>& sections,
    const lsdtt_manning_reach&              reach,
    const lsdtt_channel_profile_lookup&     profile,
    const lsdtt_manning_solver_options&     opts = lsdtt_manning_solver_options{});

// Validates the reaches, runs compute_manning_widths_for_reach over each and
// concatenates the results, sorted by distance_from_source_m.
[[nodiscard]] std::vector<lsdtt_manning_width_record>
compute_manning_widths_all_reaches(
    const std::vector<lsdtt_cross_section>&    sections,
    const std::vector<lsdtt_manning_reach>&    reaches,
    const lsdtt_channel_profile_lookup&        profile,
    const lsdtt_manning_solver_options&        opts = lsdtt_manning_solver_options{});

// One point per bank (left/right) per valid cross section, via
// lsdtt_data_loader. Vector formats only (FlatGeobuf preferred, gpkg also
// supported). format == csv throws std::invalid_argument: a point CSV must
// carry EPSG:4326 latitude/longitude, which is not implemented here; use
// write_manning_widths_csv for tabular output.
void write_manning_bank_points(
    const std::vector<lsdtt_manning_width_record>& results,
    const lsdtt_raster_meta&                        meta,
    const std::string&                              out_path,
    point_output_format                             format,
    bool                                             spatial_index,
    const std::string&                               layer_name);

// One row per cross section (valid and invalid), with quality_flag.
// Returns the number of valid rows.
std::size_t write_manning_widths_csv(
    const std::vector<lsdtt_manning_width_record>& results,
    const std::string&                              out_path);

// -----------------------------------------------------------------------------
// compute_channel_profile / write_channel_profile_csv
//
// Distance from source, drainage area, thalweg elevation, and slope per
// cross section -- input to the Python reach-identification workflow AND
// to the Manning solve above. Independent of the Manning solve; can be run
// before reach parameters exist.
//
// REQUIRES that the channel network `sections` was extracted from was built
// from a SINGLE source point (source_flat). Distance-from-source is
// computed as one downstream walk from source_flat; any cross section
// not on that path (i.e. the network had more than one source) will
// cause a std::runtime_error.
//
// SLOPE is the along-channel bed slope: the least-squares slope of
// elevation against along-channel distance, fitted over the D8 channel
// path within +/- slope_window_half_length_metres of each cross section's
// centre. It is measured on `filled_dem` (the pit-filled DEM used for
// routing). It is NOT the magnitude of the local surface gradient, which
// would include the cross-channel component and bias slope high.
//
// CAVEAT: pit filling imposes a minimum gradient (fill_min_slope) on flats.
// On a low-gradient reach that imposed gradient is a real fraction of the
// fitted slope, and Manning discharge scales as sqrt(slope).
//
// slope_valid is false when the window holds fewer than 3 valid path cells,
// the valid cells span less than slope_window_half_length_metres along the
// channel (near the source or outlet a one-sided window can still pass), or
// the fitted slope is not strictly positive.
//
// thalweg_elevation_m is the lowest FINITE valid sample within +/-
// thalweg_search_half_width_metres of each section's centre (offset 0),
// matching the thalweg used by the Manning solve. NaN if no valid sample
// lies in the window.
// -----------------------------------------------------------------------------
[[nodiscard]] std::vector<lsdtt_channel_profile_row> compute_channel_profile(
    const std::vector<lsdtt_cross_section>& sections,
    const lsdtt_flow_base&                  flow,
    const lsdtt_raster_cpu&                 filled_dem,
    double                                   slope_window_half_length_metres,
    std::int64_t                             source_flat,
    double                                   thalweg_search_half_width_metres);

void write_channel_profile_csv(
    const std::vector<lsdtt_channel_profile_row>& rows,
    const std::string&                            out_path);

// -----------------------------------------------------------------------------
// Reach template generation
//
// detect_reach_breaks: identifies reach breaks from drainage area jumps,
// enforces minimum spacing and reach length, and builds reach bounds with
// coordinates for QGIS visualization.
//
// sample_reference_xs: selects N evenly-spaced cross sections per reach as
// candidates for reference XS (for discharge calculation).
//
// write_reaches_template_csv: writes reach bounds with empty n/Q columns.
//
// write_reference_xs_candidates_csv: writes candidate XS per reach.
// -----------------------------------------------------------------------------
[[nodiscard]] std::vector<lsdtt_reach_template_row> detect_reach_breaks(
    const std::vector<lsdtt_channel_profile_row>& profile,
    const std::vector<lsdtt_cross_section>&       sections,
    double relative_jump_threshold,
    double abs_jump_fraction_of_max,
    double min_break_spacing_metres,
    double min_reach_length_metres);

[[nodiscard]] std::vector<lsdtt_reference_xs_candidate> sample_reference_xs(
    const std::vector<lsdtt_channel_profile_row>&      profile,
    const std::vector<lsdtt_reach_template_row>&       reaches,
    const std::vector<lsdtt_cross_section>&            sections,
    int n_per_reach);

void write_reaches_template_csv(
    const std::vector<lsdtt_reach_template_row>& reaches,
    const std::string&                           out_path);

void write_reference_xs_candidates_csv(
    const std::vector<lsdtt_reference_xs_candidate>& candidates,
    const std::string&                               out_path);

void write_reaches_template_lines(
    const std::vector<lsdtt_reach_template_row>& reaches,
    const std::string&                           out_path,
    const std::string&                           crs_wkt,
    lsdtt::vector_output_format                  format,
    bool                                         spatial_index);

void write_reference_xs_candidates_points(
    const std::vector<lsdtt_reference_xs_candidate>& candidates,
    const std::string&                               out_path,
    const std::string&                               crs_wkt,
    lsdtt::vector_output_format                      format,
    bool                                             spatial_index);

} // namespace lsdtt