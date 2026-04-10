#pragma once

// lsdtt_binary_array2d_tracing.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Binary (mask) boundary tracing utilities.
// Implements a Moore-neighbour boundary tracing algorithm based on
// Matlab code presented in:
//
// Gonzales, Woods, Eddins, Digitial Image Processing Using MATLAB, 
// pp. 557-560.
//
// This module lives in utils because:
//   - operates on lsdtt_array2d_u8 and not an lsdtt_raster_cpu, 
//     because it is likely you want to crop your raster data to
//     to image of interest first. 
//   - returns row/col coordinates in the input array coordinate system
//   - does not know about CRS, cell size, flow routing, etc so no
//     raster metadata is known to this object.
//
// Conventions:
//   - Any nonzero input pixel is treated as foreground (1).
//   - Zero is treated as background.
//   - Returned boundary is in the order produced by the Moore algorithm,
//     using a clockwise neighbour search table.
//
// If no foreground pixel exists, returns an empty vector.
//
// Copyright (C) 2026 Rutger Vervoordeldonk
// Released under the GNU General Public Licence v3 (GPL-3.0)

#include <vector>

#include "lsdtt_array2d.hpp"
#include "../core/lsdtt_types.hpp"

namespace lsdtt
{

// -----------------------------------------------------------------------------
// trace_boundary_pixels
//
// Runs a Moore-neighbour boundary tracing algorithm on a binary image.
// This is a direct analogue of LSDTopoTools2 trace_boundary_pixels.
//
// Parameters
//   bw           - binary image (0 = background, nonzero = foreground)
//   connectivity - four or eight connectivity
//
// Returns
//   Vector of boundary pixel coordinates (row/col) in the input bw grid.
//
// Notes
//   - The algorithm pads the input with a 1-cell background border internally.
//   - The start pixel is chosen by scanning the padded image row-by-row,
//     selecting the first foreground pixel that has a background pixel above it.
// -----------------------------------------------------------------------------
[[nodiscard]] std::vector<lsdtt_pixel_rc>
trace_boundary_pixels(
    const lsdtt_array2d_u8& bw,
    int connectivity // either 4 or 8
);

} // namespace lsdtt
