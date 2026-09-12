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

#include <cstdint>
#include <array>
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

// -----------------------------------------------------------------------------
// resolve_active_mask_corners
//
// The four corners of an active mask, as flat indices, in the order
// NW, NE, SW, SE.
//
// Found by scanning the FIRST and LAST COLUMNS only: the topmost and bottommost
// active pixel of each. That is what makes it work on a non-rectangular domain
// -- a sinuous LEM mask, say -- where the active area does not reach the
// rectangle's corners and the true corners sit somewhere down the side edges.
//
// Throws std::runtime_error if either edge column is entirely masked out, since
// there is then no corner to find on that side.
//
// Takes the MASK rather than a raster: nothing here needs elevations, and this
// header does not otherwise depend on the raster types.
// -----------------------------------------------------------------------------
[[nodiscard]] std::array<std::int64_t, 4> resolve_active_mask_corners(
    const lsdtt_array2d_u8& mask
);

// -----------------------------------------------------------------------------
// split_boundary_into_cardinal_segments
//
// Trace the mask boundary and cut it into four segments at `corners`, returned
// indexed by side: [0] north, [1] east, [2] south, [3] west.
//
// `corners` is in the order resolve_active_mask_corners returns (NW, NE, SW,
// SE). Each corner is matched to its NEAREST traced boundary pixel, because a
// corner taken from the edge columns need not itself be on the traced loop;
// the matched positions are then sorted into trace order, so the four arcs come
// out in the order the trace visits them whichever way round it runs.
//
// A segment that straddles the end of the trace array is stitched across the
// wrap. Both endpoints are inclusive, so consecutive segments share a corner
// pixel.
//
// The side each arc belongs to is decided by WHICH TWO CORNERS bound it, not by
// its position in the trace, so the result is labelled correctly regardless of
// the trace's start point or direction.
//
// Throws std::runtime_error if the mask has no boundary pixels.
// -----------------------------------------------------------------------------
[[nodiscard]] std::vector<std::vector<std::int64_t>>
split_boundary_into_cardinal_segments(
    const lsdtt_array2d_u8&            mask,
    const std::array<std::int64_t, 4>& corners
);

} // namespace lsdtt
