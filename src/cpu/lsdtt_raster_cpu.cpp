// lsdtt_raster_cpu.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2025 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html

#include "lsdtt_raster_cpu.hpp"

#include <stdexcept>

namespace lsdtt {

// -----------------------------------------------------------------------------
// Constructors
// -----------------------------------------------------------------------------

lsdtt_raster_cpu::lsdtt_raster_cpu(const lsdtt_data_loader& loader)
{
    read_from_loader(loader);
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

// -----------------------------------------------------------------------------
// I/O
// -----------------------------------------------------------------------------

void lsdtt_raster_cpu::read_from_loader(const lsdtt_data_loader& loader)
{
    lsdtt_raster_tile tile = loader.read_full();
    meta_ = tile.meta;

    data_ = lsdtt_array2d_f(
        meta_.n_rows,
        meta_.n_cols,
        std::move(tile.data)
    );

    mask_ = lsdtt_array2d_u8(
        meta_.n_rows,
        meta_.n_cols,
        std::move(tile.mask)
    );
}

void lsdtt_raster_cpu::write_to_loader(lsdtt_data_loader& loader)
{
    loader.write_full(data_.get_vector(), meta_, "");
}

// -----------------------------------------------------------------------------
// Core filtering — placeholders
// -----------------------------------------------------------------------------

void lsdtt_raster_cpu::smooth_gaussian(double /*sigma*/)
{
    throw std::runtime_error(
        "lsdtt_raster_cpu: smooth_gaussian() not yet implemented"
    );
}

void lsdtt_raster_cpu::smooth_bilateral(double /*sigma_space*/,
                                         double /*sigma_colour*/)
{
    throw std::runtime_error(
        "lsdtt_raster_cpu: smooth_bilateral() not yet implemented"
    );
}

} // namespace lsdtt
