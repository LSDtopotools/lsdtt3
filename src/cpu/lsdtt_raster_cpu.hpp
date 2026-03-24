// lsdtt_raster_cpu.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2025 Simon M. Mudd
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
//   std::vector<uint8_t>& are zero-cost — they return a reference into the
//   underlying vector of the array2d. get_data_array() / get_mask_array()
//   return the lsdtt_array2d directly for code that needs 2D indexing.
//
// Analysis free functions:
//   All CPU analysis free functions in src/cpu/ take lsdtt_raster_cpu const&
//   directly. They access data through get_data(), get_mask(), get_meta(),
//   get_data_array(), or get_mask_array() as appropriate.
// -----------------------------------------------------------------------------

#pragma once

#include <cstdint>
#include <vector>

#include "../core/lsdtt_raster_base.hpp"
#include "../core/lsdtt_raster_meta.hpp"
#include "../core/lsdtt_data_loader.hpp"
#include "../utils/lsdtt_array2d.hpp"

namespace lsdtt {

class lsdtt_raster_cpu : public lsdtt_raster_base
{
public:

    // -------------------------------------------------------------------------
    // Constructors
    // -------------------------------------------------------------------------

    // Default — produces an empty, invalid raster
    lsdtt_raster_cpu() = default;

    // Construct from a data loader — reads the full raster into memory
    explicit lsdtt_raster_cpu(const lsdtt_data_loader& loader);

    // Construct from a tile read by the loader
    explicit lsdtt_raster_cpu(const lsdtt_raster_tile& tile);

    // Construct from lsdtt_array2d (move) — preferred for analysis results
    lsdtt_raster_cpu(
        lsdtt_array2d_f&&        data,
        lsdtt_array2d_u8&&       mask,
        const lsdtt_raster_meta& meta
    );

    // Construct from std::vector (move) — for backwards compatibility
    // with code that builds result vectors directly
    lsdtt_raster_cpu(
        std::vector<float>&&     data,
        std::vector<uint8_t>&&   mask,
        const lsdtt_raster_meta& meta
    );

    // -------------------------------------------------------------------------
    // lsdtt_array2d accessors — preferred for analysis code
    // Zero-copy: returns reference into internal storage.
    // -------------------------------------------------------------------------

    const lsdtt_array2d_f&  get_data_array() const { return data_; }
    lsdtt_array2d_f&        get_data_array()       { return data_; }
    const lsdtt_array2d_u8& get_mask_array() const { return mask_; }

    // -------------------------------------------------------------------------
    // std::vector accessors — implement lsdtt_raster_base interface
    // Zero-copy: returns reference into the vector owned by lsdtt_array2d.
    // -------------------------------------------------------------------------

    const std::vector<float>&   get_data() const { return data_.get_vector(); }
    std::vector<float>&         get_data()       { return data_.get_vector(); }
    const std::vector<uint8_t>& get_mask() const { return mask_.get_vector(); }

    // -------------------------------------------------------------------------
    // Metadata accessor — implement lsdtt_raster_base interface
    // -------------------------------------------------------------------------

    const lsdtt_raster_meta& get_meta() const override { return meta_; }

    // -------------------------------------------------------------------------
    // I/O — implement lsdtt_raster_base interface
    // -------------------------------------------------------------------------

    void write_to_loader(lsdtt_data_loader& loader)        override;
    void read_from_loader(const lsdtt_data_loader& loader) override;

    // -------------------------------------------------------------------------
    // Core filtering — implement lsdtt_raster_base interface
    // -------------------------------------------------------------------------

    void smooth_gaussian(double sigma)                             override;
    void smooth_bilateral(double sigma_space, double sigma_colour) override;


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

    [[nodiscard]] float min_value() const;
    [[nodiscard]] float max_value() const;
    [[nodiscard]] float mean_value() const;
    [[nodiscard]] float median_value() const;
    [[nodiscard]] float percentile(double p) const;   // p in [0, 100]
    [[nodiscard]] float stddev_value() const;         // population stddev
    [[nodiscard]] float range_value() const;

    void add_value(float value);
    void multiply_value(float factor);
    void clip(float lower, float upper);

private:

    void validate_storage_sizes(const char* function_name) const;
    [[nodiscard]] std::int64_t checked_cell_count(const char* function_name) const;
    [[nodiscard]] std::vector<float> collect_valid_values() const;

    lsdtt_array2d_f   data_;
    lsdtt_array2d_u8  mask_;
    lsdtt_raster_meta meta_;
};

} // namespace lsdtt
