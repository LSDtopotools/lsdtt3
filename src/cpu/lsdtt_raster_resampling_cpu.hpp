// Lightweight, GDAL-free resampling of in-memory float rasters.

#pragma once

#include "lsdtt_raster_cpu.hpp"

namespace lsdtt {

enum class raster_resampling_method
{
    nearest,
    linear,
    mean,
    sum,
    minimum,
    maximum,
    cubic,
    mode,
    conservative_mean
};

enum class power_of_two_direction
{
    lower,
    higher,
    nearest
};

// Fill an already allocated destination raster. Source and destination must
// use the same CRS; their extents may differ.
void resample_raster(
    const lsdtt_raster_cpu& source,
    lsdtt_raster_cpu& destination,
    raster_resampling_method method);

// Factor-two paths. Both rasters must cover the same extent.
void double_raster(
    const lsdtt_raster_cpu& source,
    lsdtt_raster_cpu& destination,
    raster_resampling_method method);

void halve_raster(
    const lsdtt_raster_cpu& source,
    lsdtt_raster_cpu& destination,
    raster_resampling_method method);

// Pad north and east with nodata. Cell size and south-west origin are kept.
lsdtt_raster_cpu pad_to_nearest_power_of_two(
    const lsdtt_raster_cpu& source);

// Resize both dimensions to a power of two. This preserves the full extent,
// so the row and column scale factors must agree for square cells.
lsdtt_raster_cpu interpolate_to_nearest_power_of_two(
    const lsdtt_raster_cpu& source,
    power_of_two_direction direction,
    raster_resampling_method method);

} // namespace lsdtt
