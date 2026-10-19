// lsdtt_erodibility_calibration_cpu.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// Calibration of the stream-power erodibility K from basin-averaged erosion
// rates (e.g. cosmogenic 10Be): the sample input.
//
// A sample is a basin outlet with an erosion rate E and its 1-sigma
// uncertainty, in units the user must state (there is no default unit):
//
//   mm_per_kyr  = m_per_myr   1e-6 m/yr
//   cm_per_kyr                1e-5 m/yr
//   mm_per_yr                 1e-3 m/yr
//   m_per_yr                  1    m/yr
//
// Units are matched case-insensitively after trimming; m_per_myr is an alias of
// mm_per_kyr. Mass fluxes (g cm^-2 yr^-1, as OCTOPUS's EBE_GCMYR) are refused
// with a message to divide by the rock density. When the file has an
// erosion_rate_units column, every row must name the same unit as the
// parameter (both canonicalised), otherwise the read throws naming the first
// row that disagrees; this is how a file written by lsdtt-fetch-erosion-rates
// protects itself against a wrong parameter.
//
// The reader (read_erosion_rate_samples) goes through
// load_point_layer_with_attributes, so a sample file is either a CSV with
// latitude / longitude (EPSG:4326) or any OGR point layer in any CRS, and
// every input row comes back in input order with its point status (a point
// outside the DEM is kept, with status outside_raster_extent).
//
// Per-row checks made here:
//   * bad_rate             the rate is missing, not a number or <= 0
//   * missing_uncertainty  the uncertainty is missing, not a number or <= 0
// Both rows are kept; status is "" for a row that passed (later steps of the
// calibration fill it). The reported area (km^2 in the file) is stored in m^2,
// NaN when the column is absent or the value is blank or <= 0.
//
// Nothing here prints.
// -----------------------------------------------------------------------------

#pragma once

#include <limits>
#include <string>
#include <vector>

#include "../core/lsdtt_point_input.hpp"
#include "../core/lsdtt_raster_meta.hpp"

namespace lsdtt
{

// The canonical name of an erosion-rate unit (one of mm_per_kyr, cm_per_kyr,
// mm_per_yr, m_per_yr). Throws std::invalid_argument for an empty or unknown
// unit (listing the accepted ones) and for a mass flux.
[[nodiscard]] std::string canonical_erosion_rate_units(const std::string& units);

// The factor converting a rate in `units` to m/yr (exact powers of ten).
// Throws as canonical_erosion_rate_units.
[[nodiscard]] double erosion_rate_units_to_m_per_yr(const std::string& units);

// "mm_per_kyr, m_per_myr, cm_per_kyr, mm_per_yr, m_per_yr", for messages.
[[nodiscard]] std::string accepted_erosion_rate_units_list();

struct lsdtt_erosion_rate_input_options
{
    std::string path;
    std::string layer_name;                                  // vector input only
    std::string sample_id_column = "sample_id";
    std::string latitude_column = "latitude";                // CSV only
    std::string longitude_column = "longitude";              // CSV only
    std::string rate_column = "erosion_rate";
    std::string uncertainty_column = "erosion_rate_uncertainty";
    std::string units;                                       // required
    std::string reported_area_column = "reported_area_km2"; // "" = none
};

struct lsdtt_erosion_rate_sample
{
    static constexpr double nan = std::numeric_limits<double>::quiet_NaN();
    int input_row = 0;                    // 1-based data row (feature) number
    std::string sample_id;
    lsdtt_input_point point;              // in the raster CRS, with its status
    double rate_input = nan;              // as in the file
    double uncertainty_input = nan;       // as in the file
    double E_m_per_yr = nan;
    double sigma_E_m_per_yr = nan;
    double reported_area_m2 = nan;        // NaN = none
    std::string status;                   // "" | bad_rate | missing_uncertainty
    std::string status_detail;
};

struct lsdtt_erosion_rate_read_result
{
    std::vector<lsdtt_erosion_rate_sample> samples;   // input order
    std::string units;                                // canonical
    double units_to_m_per_yr = 1.0;
    bool has_reported_area = false;                   // the area column is in the file
    bool has_sample_id = false;                       // the id column is in the file
    bool has_units_column = false;                    // the file has erosion_rate_units
    std::vector<std::string> columns;                 // every attribute column found
};

// Reads every row of opts.path into the raster CRS of `meta`. Throws
// std::invalid_argument for bad units, a missing file ("erosion_rate_fname
// does not exist: <path>"), a missing rate or uncertainty column (naming the
// parameter and listing the columns found), or a units column that disagrees
// with opts.units ("... disagrees with erosion_rate_units ..."); dataset-level
// read errors propagate as std::runtime_error.
[[nodiscard]] lsdtt_erosion_rate_read_result read_erosion_rate_samples(
    const lsdtt_raster_meta& meta, const lsdtt_erosion_rate_input_options& opts);

} // namespace lsdtt
