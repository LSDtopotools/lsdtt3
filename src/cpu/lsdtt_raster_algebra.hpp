// lsdtt_raster_algebra.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2025 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html
//
// -----------------------------------------------------------------------------
// Raster algebra free functions (CPU)
//
// All operations are element-wise and produce a new lsdtt_raster_cpu.
// The output raster inherits its metadata from the left-hand operand.
//
// Nodata semantics:
//   An output cell is nodata whenever either input cell is nodata.
//   For raster_log / raster_log10 the output is additionally nodata if the
//   input value is <= 0 (i.e. the logarithm is undefined or -inf).
//
// Bounds checking:
//   Every function asserts that both inputs share identical n_rows, n_cols,
//   cell_size, x_min, y_min, and no_data_value. A std::invalid_argument
//   exception is thrown if they differ.
//
// Parallelism:
//   Inner loops are parallelised with OpenMP (if available at compile time).
//   No synchronisation is needed because each output element is written
//   by exactly one thread.
// -----------------------------------------------------------------------------

#pragma once

#include "lsdtt_raster_cpu.hpp"

namespace lsdtt {
namespace algebra {

// -----------------------------------------------------------------------------
// Metadata compatibility check
//
// Throws std::invalid_argument if the two rasters are not compatible for
// element-wise algebra (different dimensions, cell size, origin, or nodata).
// Call this from your own code if you want to validate before passing to
// the algebra functions (they call it internally too).
// -----------------------------------------------------------------------------

void assert_compatible(
    const lsdtt_raster_cpu& a,
    const lsdtt_raster_cpu& b
);

// -----------------------------------------------------------------------------
// Arithmetic operations
//
// result(i) = a(i) op b(i)
// result(i) = nodata  if a(i) is nodata OR b(i) is nodata
// -----------------------------------------------------------------------------

/// Element-wise addition: result = a + b
lsdtt_raster_cpu raster_add(
    const lsdtt_raster_cpu& a,
    const lsdtt_raster_cpu& b
);

/// Element-wise subtraction: result = a - b
lsdtt_raster_cpu raster_subtract(
    const lsdtt_raster_cpu& a,
    const lsdtt_raster_cpu& b
);

/// Element-wise multiplication: result = a * b
lsdtt_raster_cpu raster_multiply(
    const lsdtt_raster_cpu& a,
    const lsdtt_raster_cpu& b
);

/// Element-wise division: result = a / b
/// result(i) is additionally nodata if b(i) == 0.
lsdtt_raster_cpu raster_divide(
    const lsdtt_raster_cpu& a,
    const lsdtt_raster_cpu& b
);

/// Element-wise minimum: result = min(a, b)
lsdtt_raster_cpu raster_min(
    const lsdtt_raster_cpu& a,
    const lsdtt_raster_cpu& b
);

/// Positive difference: result = max(a - b, 0)
lsdtt_raster_cpu raster_positive_difference(
    const lsdtt_raster_cpu& a,
    const lsdtt_raster_cpu& b
);

// -----------------------------------------------------------------------------
// Unary transforms
// -----------------------------------------------------------------------------

/// Clamp valid pixels from below: result = max(a, lower)
lsdtt_raster_cpu raster_clamp_min(
    const lsdtt_raster_cpu& a,
    float                   lower
);

// -----------------------------------------------------------------------------
// Logarithm operations
//
// result(i) = log(a(i))   (natural) / log10(a(i))
// result(i) = nodata  if a(i) is nodata OR a(i) <= 0
// -----------------------------------------------------------------------------

/// Natural logarithm: result = ln(a)
lsdtt_raster_cpu raster_log(const lsdtt_raster_cpu& a);

/// Base-10 logarithm: result = log10(a)
lsdtt_raster_cpu raster_log10(const lsdtt_raster_cpu& a);

} // namespace algebra
} // namespace lsdtt
