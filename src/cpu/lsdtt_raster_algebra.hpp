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
// operations that take one input raster and produce one output raster by transforming each pixel value.
// -----------------------------------------------------------------------------

/// Absolute value: result = abs(a)
lsdtt_raster_cpu raster_abs(
    const lsdtt_raster_cpu& a
);

/// Clamp valid pixels from below: result = max(a, lower)
// If any pixel is less than lower, it gets raised to the value of lower
lsdtt_raster_cpu raster_clamp_min(
    const lsdtt_raster_cpu& a,
    float                   lower
);

/// Clamp valid pixels from above: result = min(a, upper)
// If any pixel is greater than upper, it gets lowered to the value of upper
lsdtt_raster_cpu raster_clamp_max(
    const lsdtt_raster_cpu& a,
    float                   upper
);


/// Raise valid pixels by a fixed amount: result = a + amount
lsdtt_raster_cpu raster_raise(
    const lsdtt_raster_cpu& a,
    float                   amount
);

/// Lower valid pixels by a fixed amount: result = a - amount
lsdtt_raster_cpu raster_lower(
    const lsdtt_raster_cpu& a,
    float                   amount
);

/// Scale valid pixels by a scalar factor: result = a * factor
lsdtt_raster_cpu raster_scale(
    const lsdtt_raster_cpu& a,
    float                   factor
);

/// Clamp valid pixels to [lower, upper]: result = clamp(a, lower, upper)
/// Throws std::invalid_argument if lower > upper.
lsdtt_raster_cpu raster_clip(
    const lsdtt_raster_cpu& a,
    float                   lower,
    float                   upper
);



/// Positive part: result = max(a, 0)
lsdtt_raster_cpu raster_positive_part(
    const lsdtt_raster_cpu& a
);

/// Negative part: result = max(-a, 0)
lsdtt_raster_cpu raster_negative_part(
    const lsdtt_raster_cpu& a
);

/// Min-max normalise valid pixels to [0, 1].
/// If all valid pixels have the same value, all valid outputs become 0.
/// Nodata remains nodata.
lsdtt_raster_cpu raster_normalise_unit_interval(
    const lsdtt_raster_cpu& a
);


// -----------------------------------------------------------------------------
// Masking and conditional selection
// -----------------------------------------------------------------------------

/// Binary threshold from above:
/// result = 1 where a >= threshold, else 0.
/// Valid input pixels remain valid in the output; nodata remains nodata.
lsdtt_raster_cpu raster_binary_true_above_threshold(
    const lsdtt_raster_cpu& a,
    float                   threshold
);

/// Binary threshold from below:
/// result = 1 where a <= threshold, else 0.
/// Valid input pixels remain valid in the output; nodata remains nodata.
lsdtt_raster_cpu raster_binary_true_below_threshold(
    const lsdtt_raster_cpu& a,
    float                   threshold
);



/// Keep the values raster where keep_mask > 0.5; set nodata elsewhere.
/// Both rasters must be compatible. Output metadata follow values.
lsdtt_raster_cpu raster_mask_by_binary_keep(
    const lsdtt_raster_cpu& values,
    const lsdtt_raster_cpu& keep_mask
);


/// Keep the values raster where selector <= max_value;
/// set nodata elsewhere. Both rasters must be compatible.
lsdtt_raster_cpu raster_mask_keep_where_less_or_eq(
    const lsdtt_raster_cpu& values,
    const lsdtt_raster_cpu& selector,
    float                   max_value
);

/// Keep the values raster where selector >= min_value;
/// set nodata elsewhere. Both rasters must be compatible.
lsdtt_raster_cpu raster_mask_keep_where_greater_or_eq(
    const lsdtt_raster_cpu& values,
    const lsdtt_raster_cpu& selector,
    float                   min_value
);


/// Keep foreground pixels from a binary mask only where
/// selector <= max_value. Background remains 0.
/// Validity follows the intersection of the two rasters.
lsdtt_raster_cpu raster_binary_keep_true_where_less_or_eq(
    const lsdtt_raster_cpu& binary_mask,
    const lsdtt_raster_cpu& selector,
    float                   max_value
);

/// Keep foreground pixels from a binary mask only where
/// selector >= min_value. Background remains 0.
/// Validity follows the intersection of the two rasters.
lsdtt_raster_cpu raster_binary_keep_true_where_greater_or_eq(
    const lsdtt_raster_cpu& binary_mask,
    const lsdtt_raster_cpu& selector,
    float                   min_value
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


// -----------------------------------------------------------------------------
// Scalar summary statistics
//
// Each function scans all valid pixels (mask != 0) and returns a single float.
// Throws std::runtime_error if the raster contains no valid pixels.
// raster_percentile additionally throws std::invalid_argument if p is outside
// [0, 100].
// -----------------------------------------------------------------------------

/// Minimum value over all valid pixels.
float raster_min_value(const lsdtt_raster_cpu& a);

/// Maximum value over all valid pixels.
float raster_max_value(const lsdtt_raster_cpu& a);

/// Arithmetic mean over all valid pixels.
float raster_mean_value(const lsdtt_raster_cpu& a);

/// Median value (50th percentile) over all valid pixels.
float raster_median_value(const lsdtt_raster_cpu& a);

/// Percentile over all valid pixels. p must be in [0, 100].
float raster_percentile(const lsdtt_raster_cpu& a, double p);

/// Population standard deviation over all valid pixels.
float raster_stddev_value(const lsdtt_raster_cpu& a);

/// Range (max − min) over all valid pixels.
float raster_range_value(const lsdtt_raster_cpu& a);

} // namespace algebra
} // namespace lsdtt
