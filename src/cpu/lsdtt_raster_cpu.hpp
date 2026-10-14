// lsdtt_raster_cpu.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd and Rutger Vervoordeldonk
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// lsdtt_raster_cpu
//
// Concrete CPU implementation of lsdtt_raster_base.
// Owns the data and mask arrays for a raster tile or full raster.
//
// Storage:
//   Data and mask are stored as lsdtt_array2d (which wraps std::vector).
//   This means get_data() / get_mask() returning std::vector<float>& /
//   std::vector<uint8_t>& are zero-cost - they return a reference into the
//   underlying vector of the array2d. get_data_array() / get_mask_array()
//   return the lsdtt_array2d directly for code that needs 2D indexing.
//
// Analysis free functions:
//   All CPU analysis free functions in src/cpu/ take lsdtt_raster_cpu const&
//   directly. They access data through get_data(), get_mask(), get_meta(),
//   get_data_array(), or get_mask_array() as appropriate.
// -----------------------------------------------------------------------------

#pragma once

#include <cstddef>
#include <cstdint>
#include <string>
#include <ostream>
#include <vector>
#include <memory>

#include "../core/lsdtt_raster_base.hpp"
#include "../core/lsdtt_raster_meta.hpp"
#include "../core/lsdtt_raster_result.hpp"
#include "../core/lsdtt_data_loader.hpp"
#include "../utils/lsdtt_array2d.hpp"

namespace lsdtt {

class lsdtt_raster_cpu : public lsdtt_raster_base
{
public:

    // -------------------------------------------------------------------------
    // Constructors
    // -------------------------------------------------------------------------

    // Default - produces an empty, invalid raster
    lsdtt_raster_cpu() = default;

    // Construct from a data loader - reads the full raster into memory
    explicit lsdtt_raster_cpu(const lsdtt_data_loader& loader);

    // Construct from a tile read by the loader
    explicit lsdtt_raster_cpu(const lsdtt_raster_tile& tile);

    // Construct from lsdtt_array2d (move) - preferred for analysis results
    lsdtt_raster_cpu(
        lsdtt_array2d_f&&        data,
        lsdtt_array2d_u8&&       mask,
        const lsdtt_raster_meta& meta
    );

    // Construct from std::vector (move) - for backwards compatibility
    // with code that builds result vectors directly
    lsdtt_raster_cpu(
        std::vector<float>&&     data,
        std::vector<uint8_t>&&   mask,
        const lsdtt_raster_meta& meta
    );

    // Construct from a core-layer raster result (move).
    // Use this to promote an lsdtt_raster_result returned by a base-class
    // engine method into a full CPU raster for further analysis.
    explicit lsdtt_raster_cpu(lsdtt_raster_result&& result);

    // -------------------------------------------------------------------------
    // Cloning
    // -------------------------------------------------------------------------
    [[nodiscard]] std::unique_ptr<lsdtt_raster_base> clone() const override
    {
        return std::make_unique<lsdtt_raster_cpu>(*this);
    }


    // -------------------------------------------------------------------------
    // lsdtt_array2d accessors - preferred for analysis code
    // Zero-copy: returns reference into internal storage.
    // -------------------------------------------------------------------------

    const lsdtt_array2d_f&  get_data_array() const { return data_; }
    lsdtt_array2d_f&        get_data_array()       { return data_; }
    const lsdtt_array2d_u8& get_mask_array() const { return mask_; }
    lsdtt_array2d_u8&       get_mask_array()       { return mask_; }

    // -------------------------------------------------------------------------
    // std::vector accessors - zero-copy into the array2d-owned vector
    // -------------------------------------------------------------------------

    const std::vector<float>&   get_data() const { return data_.get_vector(); }
    std::vector<float>&         get_data()       { return data_.get_vector(); }
    const std::vector<uint8_t>& get_mask() const { return mask_.get_vector(); }

    // -------------------------------------------------------------------------
    // Metadata accessor retained for existing CPU analysis code
    // -------------------------------------------------------------------------

    const lsdtt_raster_meta& get_meta() const { return meta_; }

    // -------------------------------------------------------------------------
    // lsdtt_raster_base interface
    // -------------------------------------------------------------------------

    [[nodiscard]] bool empty() const noexcept override
    {
        return data_.size() == 0;
    }

    [[nodiscard]] std::size_t full_grid_size() const noexcept override
    {
        return data_.size();
    }

    [[nodiscard]] const lsdtt_raster_meta& meta() const noexcept override
    {
        return meta_;
    }

    [[nodiscard]] int n_rows() const noexcept override
    {
        return meta_.n_rows;
    }

    [[nodiscard]] int n_cols() const noexcept override
    {
        return meta_.n_cols;
    }

    [[nodiscard]] double cell_size() const noexcept override
    {
        return meta_.cell_size;
    }

    [[nodiscard]] const float* data_ptr() const noexcept override
    {
        return data_.data();
    }

    [[nodiscard]] const std::uint8_t* mask_ptr() const noexcept override
    {
        return mask_.data();
    }

    void set_from_vectors(
        std::vector<float>&&        data,
        std::vector<std::uint8_t>&& mask,
        const lsdtt_raster_meta&    meta
    ) override;

    // -------------------------------------------------------------------------
    // Mask consistency
    //
    // If a pixel has nodata in the data array, force its mask value to 0.
    // Returns the number of mask entries changed.
    // -------------------------------------------------------------------------
    [[nodiscard]] std::int64_t assert_mask();

    // -------------------------------------------------------------------------
    // Whole-raster operations over valid pixels only (mask != 0)
    // -------------------------------------------------------------------------

    [[nodiscard]] std::int64_t count_valid_pixels() const;

private:

    void validate_storage_sizes(const char* function_name) const;
    [[nodiscard]] std::int64_t checked_cell_count(const char* function_name) const;
    lsdtt_array2d_f   data_;
    lsdtt_array2d_u8  mask_;
    lsdtt_raster_meta meta_;
};


// -----------------------------------------------------------------------------
// make_uniform_raster
//
// A raster of one value, shaped like `reference` and carrying its valid mask.
//
// This replaces FOUR local variants that had drifted into TWO different
// behaviours -- `make_constant_raster` in lsdtt-cosmo-prediction, `like` in
// lsdtt-cosmo-landslides, a `make_uniform` lambda in lsdtt-cosmo-tool, and a
// `uniform` helper repeated across the test files.
//
// THE DIFFERENCE IS EXPLICIT HERE BECAUSE IT WAS ACCIDENTAL THERE. Two of the
// four filled the WHOLE array with `value`, including pixels outside the mask;
// the third wrote nodata outside the mask. `get_data()` returns the raw array
// and the mask is not applied on write, so for any raster that reaches disk
// those two choices produce different files.
//
// Neither is wrong, so neither is imposed: `MaskedFill` makes the caller say
// which it means. Standardising on one would have been a silent behaviour
// change to a shipped output (`_LS_EXCAVATION.tif` is built this way).
enum class MaskedFill
{
    Value,   // fill outside the mask with `value` too -- the array is uniform
    NoData   // fill outside the mask with the metadata's no-data value
};

// -----------------------------------------------------------------------------
// Vertical units
//
// A DEM's vertical unit is NOT recorded in any format lsdtt3 reads. ENVI
// headers, ESRI .hdr/.prj and GeoTIFF all describe the HORIZONTAL CRS only, so
// a raster on a 1 m UTM grid whose elevations are in feet is indistinguishable,
// to every reader, from one whose elevations are in metres. Much US state lidar
// is delivered exactly that way (NAVD88 feet on a metre grid).
//
// Every slope, curvature, relief and gradient lsdtt3 computes is then wrong by
// 3.28, and wrong CONSISTENTLY, so nothing downstream looks broken -- the
// landscape just appears far steeper than it is.
//
// THIS IS NEVER APPLIED AUTOMATICALLY, and must not be. The unit cannot be
// inferred from the data with any certainty: a genuinely alpine catchment in
// metres and a gentle one in feet overlap in every statistic. Guessing would
// silently corrupt correct data, which is worse than the bug it fixes. The user
// must establish the vertical datum from the data provider's metadata and say
// so explicitly.
// -----------------------------------------------------------------------------
enum class vertical_unit
{
    metres,           // no conversion
    feet,             // international foot, exactly 0.3048 m
    us_survey_feet    // 1200/3937 m -- differs from the international foot by
                      // 2 ppm, i.e. ~0.6 mm over the 327 m relief of a
                      // catchment. Distinguished because US datums specify it,
                      // not because the difference is ever measurable here.
};

// Parse a vertical_unit from a parameter string. Accepts "metres"/"meters"/"m",
// "feet"/"foot"/"ft", "us_survey_feet"/"us_feet". Throws std::runtime_error on
// anything else rather than defaulting, so a typo cannot silently mean metres.
[[nodiscard]] vertical_unit parse_vertical_unit(const std::string& name);

// Multiply every valid pixel by the factor that converts `from` to metres.
// A no-op for vertical_unit::metres.
[[nodiscard]] lsdtt_raster_cpu convert_vertical_units_to_metres(
    const lsdtt_raster_cpu& dem,
    vertical_unit           from);

// What a plausibility check found. Slopes are gradient magnitudes converted to
// degrees, over valid interior pixels only.
struct vertical_plausibility_report
{
    double      median_slope_deg      = 0.0;
    double      p90_slope_deg         = 0.0;
    double      p99_slope_deg         = 0.0;
    double      relief_m              = 0.0;
    // The same median if the elevations were really feet. Reported so the user
    // can see both candidates side by side; NOT a recommendation.
    double      median_slope_deg_if_feet = 0.0;
    std::size_t n_sampled             = 0;
    bool        implausible           = false;
};

// Measure the slope distribution and flag it as physically implausible when the
// MEDIAN slope exceeds `median_slope_limit_deg`.
//
// The threshold is a statement about landscapes, not about a number: soil and
// regolith fail at around 30-35 degrees, so a landscape whose MEDIAN pixel is
// steeper than that cannot be soil-mantled terrain. Isolated cliffs are
// expected and do not move a median, which is why the median is tested rather
// than a maximum. A real bedrock-cliff landscape can trip this legitimately --
// hence a warning, never an error and never an automatic conversion.
[[nodiscard]] vertical_plausibility_report check_vertical_plausibility(
    const lsdtt_raster_cpu& dem,
    double                  median_slope_limit_deg = 30.0);

// Run the check and, if it fails, write a warning to `out` naming the feet
// hypothesis and the parameter that would fix it. Returns the report.
vertical_plausibility_report warn_if_vertically_implausible(
    const lsdtt_raster_cpu& dem,
    std::ostream&           out,
    const std::string&      vertical_units_parameter = "vertical_units");

[[nodiscard]] lsdtt_raster_cpu make_uniform_raster(
    const lsdtt_raster_cpu& reference,
    float                   value,
    MaskedFill              outside);

} // namespace lsdtt
