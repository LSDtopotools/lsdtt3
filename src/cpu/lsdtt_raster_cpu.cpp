// lsdtt_raster_cpu.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html

#include "lsdtt_raster_cpu.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace lsdtt {

// -----------------------------------------------------------------------------
// Constructors
// -----------------------------------------------------------------------------

lsdtt_raster_cpu::lsdtt_raster_cpu(const lsdtt_data_loader& loader)
{
    lsdtt_raster_tile tile = loader.read_full();

    set_from_vectors(
        std::move(tile.data),
        std::move(tile.mask),
        tile.meta
    );
}

lsdtt_raster_cpu::lsdtt_raster_cpu(const lsdtt_raster_tile& tile)
    : meta_(tile.meta)
{
    data_ = lsdtt_array2d_f(
        tile.meta.n_rows,
        tile.meta.n_cols,
        lsdtt_array2d_f::from_ptr,
        tile.data.data()
    );

    mask_ = lsdtt_array2d_u8(
        tile.meta.n_rows,
        tile.meta.n_cols,
        lsdtt_array2d_u8::from_ptr,
        tile.mask.data()
    );
}

lsdtt_raster_cpu::lsdtt_raster_cpu(
    lsdtt_array2d_f&&        data,
    lsdtt_array2d_u8&&       mask,
    const lsdtt_raster_meta& meta
)
    : data_(std::move(data))
    , mask_(std::move(mask))
    , meta_(meta)
{
}

lsdtt_raster_cpu::lsdtt_raster_cpu(
    std::vector<float>&&     data,
    std::vector<uint8_t>&&   mask,
    const lsdtt_raster_meta& meta
)
    : data_(meta.n_rows, meta.n_cols, std::move(data))
    , mask_(meta.n_rows, meta.n_cols, std::move(mask))
    , meta_(meta)
{
}

lsdtt_raster_cpu::lsdtt_raster_cpu(lsdtt_raster_result&& result)
    : data_(result.meta.n_rows, result.meta.n_cols, std::move(result.data))
    , mask_(result.meta.n_rows, result.meta.n_cols, std::move(result.mask))
    , meta_(result.meta)
{
}

// -----------------------------------------------------------------------------
// Loader-to-backend transfer
// -----------------------------------------------------------------------------

void lsdtt_raster_cpu::set_from_vectors(
    std::vector<float>&&        data,
    std::vector<std::uint8_t>&& mask,
    const lsdtt_raster_meta&    meta
)
{
    meta_ = meta;

    data_ = lsdtt_array2d_f(
        meta.n_rows,
        meta.n_cols,
        std::move(data)
    );

    mask_ = lsdtt_array2d_u8(
        meta.n_rows,
        meta.n_cols,
        std::move(mask)
    );
}

// -----------------------------------------------------------------------------
// Private helpers
// -----------------------------------------------------------------------------

void lsdtt_raster_cpu::validate_storage_sizes(const char* function_name) const
{
    if (data_.size() != mask_.size())
    {
        throw std::runtime_error(
            std::string("lsdtt_raster_cpu::") + function_name +
            "(): data and mask sizes do not match"
        );
    }
}

std::int64_t lsdtt_raster_cpu::checked_cell_count(const char* function_name) const
{
    validate_storage_sizes(function_name);

    if (data_.size() >
        static_cast<std::size_t>(std::numeric_limits<std::int64_t>::max()))
    {
        throw std::runtime_error(
            std::string("lsdtt_raster_cpu::") + function_name +
            "(): raster contains more cells than fit in int64_t"
        );
    }

    return static_cast<std::int64_t>(data_.size());
}

// -----------------------------------------------------------------------------
// Whole-raster operations over valid pixels (mask != 0)
// -----------------------------------------------------------------------------

std::int64_t lsdtt_raster_cpu::count_valid_pixels() const
{
    const std::int64_t n_cells = checked_cell_count("count_valid_pixels");
    const std::uint8_t* mask_ptr = mask_.data();

    std::int64_t n_valid = 0;

#ifdef _OPENMP
#pragma omp parallel for schedule(static) reduction(+:n_valid)
#endif
    for (std::int64_t i = 0; i < n_cells; ++i)
    {
        if (mask_ptr[static_cast<std::size_t>(i)] != 0)
        {
            ++n_valid;
        }
    }

    return n_valid;
}

std::int64_t lsdtt_raster_cpu::assert_mask()
{
    const std::int64_t n_cells = checked_cell_count("assert_mask");

    const float*        data_ptr = data_.data();
    std::uint8_t*       mask_ptr = mask_.data();
    const float         nodata   = static_cast<float>(meta_.no_data_value);
    const bool          nodata_is_nan = std::isnan(static_cast<double>(nodata));

    std::int64_t n_changed = 0;

#ifdef _OPENMP
#pragma omp parallel for schedule(static) reduction(+:n_changed)
#endif
    for (std::int64_t i = 0; i < n_cells; ++i)
    {
        const std::size_t idx = static_cast<std::size_t>(i);
        const float z = data_ptr[idx];

        const bool is_nodata =
            nodata_is_nan
                ? std::isnan(static_cast<double>(z))
                : (z == nodata || std::isnan(static_cast<double>(z)));

        if (is_nodata && mask_ptr[idx] != 0)
        {
            mask_ptr[idx] = 0;
            ++n_changed;
        }
    }

    return n_changed;
}


lsdtt_raster_cpu make_uniform_raster(const lsdtt_raster_cpu& reference,
                                     float                   value,
                                     MaskedFill              outside)
{
    const lsdtt_raster_meta& meta = reference.get_meta();
    const float fill = (outside == MaskedFill::Value)
                     ? value
                     : static_cast<float>(meta.no_data_value);

    lsdtt_array2d<float>   d(meta.n_rows, meta.n_cols, fill);
    lsdtt_array2d<uint8_t> m(meta.n_rows, meta.n_cols, 0);

    const auto& ref_mask = reference.get_mask_array();
    for (int r = 0; r < meta.n_rows; ++r)
        for (int c = 0; c < meta.n_cols; ++c)
            if (ref_mask(r, c)) { d(r, c) = value; m(r, c) = 1; }

    return lsdtt_raster_cpu(std::move(d), std::move(m), meta);
}

} // namespace lsdtt
