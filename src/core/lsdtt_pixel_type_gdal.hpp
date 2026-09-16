// lsdtt_pixel_type_gdal.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// The GDAL side of lsdtt_pixel_type.
//
// INCLUDE THIS ONLY FROM A .cpp AT THE GDAL BOUNDARY. It names GDALDataType, so
// including it from a header would defeat the point of lsdtt_pixel_type.hpp,
// which is to keep <gdal.h> out of every translation unit that merely carries
// raster metadata around.
// -----------------------------------------------------------------------------

#ifndef LSDTT_PIXEL_TYPE_GDAL_HPP
#define LSDTT_PIXEL_TYPE_GDAL_HPP

#include "lsdtt_pixel_type.hpp"

#include <gdal.h>

namespace lsdtt
{

// GDAL sample type -> lsdtt_pixel_type. Types lsdtt3 does not model (complex,
// and the 64-bit integers GDAL 3.5 added) become lsdtt_pixel_type::unknown.
lsdtt_pixel_type pixel_type_from_gdal(GDALDataType t);

// lsdtt_pixel_type -> GDAL sample type. lsdtt_pixel_type::unknown becomes
// GDT_Unknown, which GDAL's Create() rejects -- callers that need a concrete
// band type must substitute GDT_Float32 themselves rather than have this
// function silently pick one for them.
GDALDataType pixel_type_to_gdal(lsdtt_pixel_type t);

} // namespace lsdtt

#endif // LSDTT_PIXEL_TYPE_GDAL_HPP
