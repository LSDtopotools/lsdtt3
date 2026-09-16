// lsdtt_pixel_type.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// lsdtt_pixel_type
//
// The on-disk sample type of a raster band, as lsdtt3 describes it.
//
// This exists so that lsdtt_raster_meta -- which every compute object in the
// tree carries around -- does not have to name GDALDataType, and therefore does
// not have to pull <gdal.h> into every translation unit that touches raster
// metadata. Compute objects work internally in float32 regardless; this type is
// retained so an output can be written in the same type as its input.
//
// The GDAL mapping lives in lsdtt_pixel_type_gdal.hpp, which is included ONLY by
// the .cpp files at the GDAL boundary. Nothing else should need it.
// -----------------------------------------------------------------------------

#ifndef LSDTT_PIXEL_TYPE_HPP
#define LSDTT_PIXEL_TYPE_HPP

#include <string>

namespace lsdtt
{

// Sample types lsdtt3 reads and writes. Anything GDAL supports that is not
// listed here -- the complex types, and the 64-bit integers added in GDAL 3.5 --
// maps to `unknown`, which callers treat as "write float32 instead".
enum class lsdtt_pixel_type
{
    unknown = 0,
    uint8,
    int8,
    uint16,
    int16,
    uint32,
    int32,
    float32,
    float64
};

// Human-readable name, e.g. "Float32". Matches GDAL's spelling so that log
// output does not change wording depending on which side of the boundary
// printed it.
std::string pixel_type_name(lsdtt_pixel_type t);

// True for the floating-point types.
bool pixel_type_is_float(lsdtt_pixel_type t);

// -----------------------------------------------------------------------------
// tiff_predictor_for
//
// The TIFF predictor to pair with DEFLATE/LZW/ZSTD for a given sample type:
// 3 (floating-point predictor) for float data, 2 (horizontal differencing) for
// integer data, 1 (none) for anything else.
//
// Every place that sets a predictor on a GTiff or COG must use this rather than
// a literal: a wrong predictor corrupts rather than merely failing to compress,
// so `none` is the only safe answer for a type we do not recognise.
// -----------------------------------------------------------------------------
int tiff_predictor_for(lsdtt_pixel_type t);

} // namespace lsdtt

#endif // LSDTT_PIXEL_TYPE_HPP
