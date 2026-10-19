// lsdtt_lem_parameter_rasters_cpu.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// Per-pixel LEM parameter rasters: what their values mean and whether they are
// acceptable. Where the pixels sit (grids, clipping, warping) is not this
// module's concern; that is lsdtt_lem_input_grid_cpu.
//
// Erodibility calibration. A K value is only meaningful for the m (and n) it
// was calibrated with, and, with specific area (divide_area_by_dx), for the
// cell size: K_specific(dx) = K_map * dx^m. check_lem_erodibility_calibration
// is the single place that compares what a K source states against what the
// run uses, and returns the factor that converts the source's K into the run's
// convention. It never prints: notes and warnings are returned.
// -----------------------------------------------------------------------------

#pragma once

#include <cstddef>
#include <cstdint>
#include <optional>
#include <string>
#include <vector>

#include "lsdtt_lem_input_grid_cpu.hpp"
#include "lsdtt_raster_cpu.hpp"
#include "lsdtt_stream_power_cpu.hpp"

namespace lsdtt
{

// What a K source states about its calibration. Negative / nullopt = unstated.
struct lsdtt_erodibility_calibration
{
    float  m_exponent = -1.0f;
    float  n_exponent = -1.0f;
    std::optional<erodibility_area_convention> area_convention;
    double cell_size  = -1.0;   // calibration cell size (metres), for specific_area
};

// What the run uses.
struct lsdtt_erodibility_run_convention
{
    float  m_exponent        = 0.45f;
    float  n_exponent        = 1.0f;
    bool   divide_area_by_dx = true;
    double cell_size         = 1.0;
};

struct lsdtt_erodibility_calibration_check
{
    double conversion_factor = 1.0;      // K_run = K_stated * factor; exactly 1 if nothing converts
    std::string run_units;               // describe_erodibility_units for the run
    std::vector<std::string> notes;      // conversions performed, with the factor
    std::vector<std::string> warnings;   // unstated calibration keys
};

// Throws std::invalid_argument (message contains source_label and the key) if:
//   the stated m or n differs from the run's (by more than 1e-6);
//   specific_area is stated without a cell size and the run uses map area.
// Converts, with a note:
//   map_area stated, specific-area run: factor dx_run^m;
//   specific_area stated at a cell size != the run's: factor (dx_run/dx_cal)^m;
//   specific_area stated (with a cell size), map-area run: factor dx_cal^(-m).
// Warns:
//   m, n or the convention unstated (the convention warning names the factor
//   that would apply if the K were in map-area units);
//   specific_area stated without a cell size on a specific-area run (assumed
//   calibrated at the run's cell size).
[[nodiscard]] lsdtt_erodibility_calibration_check check_lem_erodibility_calibration(
    const lsdtt_erodibility_calibration&    stated,
    const lsdtt_erodibility_run_convention& run,
    const std::string&                      source_label);

// -----------------------------------------------------------------------------
// Erodibility field values
//
// A field holds ABSOLUTE K in the run's convention, flat on the model grid. On
// every DEM-valid pixel it must be finite and > 0 (K = 0 would mean no fluvial
// erosion; switch fluvial erosion off instead). DEM-nodata pixels are stored as
// 0 and never read.
// -----------------------------------------------------------------------------
struct lsdtt_lem_erodibility_field_options
{
    float reference_erodibility = 1e-3f;      // the scalar `erodibility` (K_ref), > 0
    lsdtt_erodibility_calibration stated;     // what the source says about its K
    float model_m_exponent = 0.45f;
    float model_n_exponent = 1.0f;
    bool  divide_area_by_dx = true;           // run cell size comes from the DEM grid
    std::string source_label = "erodibility raster";
};

struct lsdtt_lem_erodibility_field
{
    std::vector<float> values;                // run convention; DEM-nodata pixels = 0
    std::size_t n_valid = 0;
    float  min_value = 0.0f, max_value = 0.0f;
    double mean_value = 0.0;
    double conversion_factor = 1.0;           // exactly 1 when nothing was converted
    std::string units;                        // describe_erodibility_units for the run
    std::vector<std::string> notes, warnings;
};

// Values already on the DEM grid (a lithology rewrite, a NumPy array). values_valid
// empty = every value is data. Calibration (check_lem_erodibility_calibration) ->
// conversion K_run = float(double(K) * factor) -> checks. Throws
// std::invalid_argument if a DEM-valid pixel is nodata in the source
//   ("<label> has nodata on <count> DEM-valid pixels; first at (r, c)")
// or has a non-finite or <= 0 value ("<label>: K must be > 0; pixel (r, c) has K = <v>").
// Never prints.
lsdtt_lem_erodibility_field validate_lem_erodibility_values(
    const std::vector<float>&                  values,
    const std::vector<std::uint8_t>&           values_valid,
    const lsdtt_raster_cpu&                    dem,
    const lsdtt_lem_erodibility_field_options& opts);

// A raster that must already lie on the DEM grid: check_lem_field_grid_alignment
// (strict), then validate_lem_erodibility_values. The alignment warnings are
// merged into .warnings.
lsdtt_lem_erodibility_field validate_lem_erodibility_field(
    const lsdtt_raster_cpu&                    erodibility_raster,
    const lsdtt_raster_cpu&                    dem,
    const lsdtt_lem_erodibility_field_options& opts);

// From disk: prepare_lem_input_raster(path, dem grid, placement, source_label),
// then validate_lem_erodibility_field. Preparation notes and warnings come
// first in .notes / .warnings.
lsdtt_lem_erodibility_field load_lem_erodibility_field(
    const std::string&                              path,
    const lsdtt_raster_cpu&                         dem,
    const lsdtt_lem_input_raster_placement_options& placement,
    const lsdtt_lem_erodibility_field_options&      opts);

// The per-step safety check run_lem_step applies to opts.erodibility_field:
// size n_rows*n_cols, and finite K > 0 on every DEM-valid pixel. Throws
// std::invalid_argument. Serial; call it outside any parallel region.
void check_lem_erodibility_field(const std::vector<float>& field, const lsdtt_raster_cpu& dem);

// out[i] = float((double(base[i]) * double(applied_K)) / double(reference_K)) for
// base[i] > 0; base[i] == 0 (nodata) stays 0. Identities: applied_K == reference_K
// gives out == base bitwise; base[i] == reference_K gives out[i] == applied_K
// bitwise. Throws std::invalid_argument if reference_K or applied_K is not > 0,
// and std::range_error naming the pixel if a result exceeds FLT_MAX or rounds to 0.
void scale_lem_erodibility_field(const std::vector<float>& base, float reference_K,
                                 float applied_K, std::vector<float>& out);

struct lsdtt_lem_erodibility_summary
{
    float  min_value  = 0.0f;
    float  max_value  = 0.0f;
    double mean_value = 0.0;
    std::size_t n_valid = 0;
};

// Min, max and mean of the field over pixels where dem_mask != 0.
[[nodiscard]] lsdtt_lem_erodibility_summary summarise_lem_erodibility_field(
    const std::vector<float>& field, const std::vector<std::uint8_t>& dem_mask);

} // namespace lsdtt
