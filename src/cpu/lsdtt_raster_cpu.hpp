// lsdtt_raster_cpu.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html
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

[[nodiscard]] lsdtt_raster_cpu make_uniform_raster(
    const lsdtt_raster_cpu& reference,
    float                   value,
    MaskedFill              outside);

} // namespace lsdtt
