// lsdtt_pixel_type.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// Both sides of lsdtt_pixel_type: the GDAL-free queries declared in
// lsdtt_pixel_type.hpp, and the GDAL mapping declared in
// lsdtt_pixel_type_gdal.hpp. They live in one translation unit so the mapping
// tables cannot drift apart.

#include "lsdtt_pixel_type.hpp"
#include "lsdtt_pixel_type_gdal.hpp"

namespace lsdtt
{

std::string pixel_type_name(lsdtt_pixel_type t)
{
    switch (t)
    {
        case lsdtt_pixel_type::uint8:   return "Byte";
        case lsdtt_pixel_type::int8:    return "Int8";
        case lsdtt_pixel_type::uint16:  return "UInt16";
        case lsdtt_pixel_type::int16:   return "Int16";
        case lsdtt_pixel_type::uint32:  return "UInt32";
        case lsdtt_pixel_type::int32:   return "Int32";
        case lsdtt_pixel_type::float32: return "Float32";
        case lsdtt_pixel_type::float64: return "Float64";
        case lsdtt_pixel_type::unknown: break;
    }
    return "Unknown";
}

bool pixel_type_is_float(lsdtt_pixel_type t)
{
    return t == lsdtt_pixel_type::float32 || t == lsdtt_pixel_type::float64;
}

int tiff_predictor_for(lsdtt_pixel_type t)
{
    switch (t)
    {
        case lsdtt_pixel_type::float32:
        case lsdtt_pixel_type::float64:
            return 3;
        case lsdtt_pixel_type::uint8:
        case lsdtt_pixel_type::int8:
        case lsdtt_pixel_type::uint16:
        case lsdtt_pixel_type::int16:
        case lsdtt_pixel_type::uint32:
        case lsdtt_pixel_type::int32:
            return 2;
        case lsdtt_pixel_type::unknown:
            break;
    }
    // Complex types, and anything a future GDAL adds: no predictor is the only
    // safe answer, because a wrong one corrupts rather than merely compressing
    // worse.
    return 1;
}

lsdtt_pixel_type pixel_type_from_gdal(GDALDataType t)
{
    switch (t)
    {
        case GDT_Byte:    return lsdtt_pixel_type::uint8;
        case GDT_UInt16:  return lsdtt_pixel_type::uint16;
        case GDT_Int16:   return lsdtt_pixel_type::int16;
        case GDT_UInt32:  return lsdtt_pixel_type::uint32;
        case GDT_Int32:   return lsdtt_pixel_type::int32;
        case GDT_Float32: return lsdtt_pixel_type::float32;
        case GDT_Float64: return lsdtt_pixel_type::float64;
        default:
            // GDT_Int8 is only defined from GDAL 3.7, and the 64-bit integers
            // only from 3.5, so they are matched by name here rather than with
            // case labels that would not compile against an older GDAL.
#if defined(GDAL_VERSION_NUM) && GDAL_VERSION_NUM >= GDAL_COMPUTE_VERSION(3,7,0)
            if (t == GDT_Int8) return lsdtt_pixel_type::int8;
#endif
            return lsdtt_pixel_type::unknown;
    }
}

GDALDataType pixel_type_to_gdal(lsdtt_pixel_type t)
{
    switch (t)
    {
        case lsdtt_pixel_type::uint8:   return GDT_Byte;
        case lsdtt_pixel_type::uint16:  return GDT_UInt16;
        case lsdtt_pixel_type::int16:   return GDT_Int16;
        case lsdtt_pixel_type::uint32:  return GDT_UInt32;
        case lsdtt_pixel_type::int32:   return GDT_Int32;
        case lsdtt_pixel_type::float32: return GDT_Float32;
        case lsdtt_pixel_type::float64: return GDT_Float64;
        case lsdtt_pixel_type::int8:
#if defined(GDAL_VERSION_NUM) && GDAL_VERSION_NUM >= GDAL_COMPUTE_VERSION(3,7,0)
            return GDT_Int8;
#else
            // No native Int8 before GDAL 3.7. Int16 holds every Int8 value
            // exactly, so widening is lossless -- unlike GDT_Byte, which would
            // wrap the negatives.
            return GDT_Int16;
#endif
        case lsdtt_pixel_type::unknown:
            break;
    }
    return GDT_Unknown;
}

} // namespace lsdtt
