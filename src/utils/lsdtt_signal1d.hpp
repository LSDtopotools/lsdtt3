// lsdtt_signal1d.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Rutger Vervoordeldonk and Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html
//
// -----------------------------------------------------------------------------
// 1D Signal Processing tools for profiles (moving average, LOWESS, derivative)
// -----------------------------------------------------------------------------

#pragma once

#include <vector>
#include <cstddef>
#include <iterator>

namespace lsdtt
{

// -----------------------------------------------------------------------------
// Smoothing Algorithms
// -----------------------------------------------------------------------------

// Applies an unweighted moving average using a sliding window context.
// By default, the window_size will be forced to an odd number internally.
std::vector<double> compute_moving_average_smoothing(
    const std::vector<double>& y,
    int window_size,
    bool periodic
);

// LOWESS smoothing (Locally Weighted Scatterplot Smoothing)
// Bandwidth determines the tricube kernel width in X-units.
// NOTE: Automatically conducts a fixed number of robustness iterations (typical 3)
std::vector<double> compute_lowess_smoothing(
    const std::vector<double>& x,
    const std::vector<double>& y,
    double bandwidth,
    int iterations = 3
);

// -----------------------------------------------------------------------------
// Profile Signal Operators
// -----------------------------------------------------------------------------

// Computes the central difference derivative (dy/dx) at each node.
// Supports spatial periodic boundaries naturally.
std::vector<double> compute_central_difference(
    const std::vector<double>& x,
    const std::vector<double>& y,
    bool periodic
);

// Searches for signal knees (points of max curvature) using the exact Kneedle 
// topology implementation (Satopaa et al., 2011).
// Returns a vector mapping internally generated indices to knees positions.
std::vector<std::size_t> find_knees(
    const std::vector<double>& x,
    const std::vector<double>& y,
    double sensitivity = 1.0
);

// -----------------------------------------------------------------------------
// Iterator and Wrapper Overloads
// -----------------------------------------------------------------------------

template<typename InputItY>
std::vector<double> compute_moving_average_smoothing(
    InputItY y_begin, InputItY y_end, int window_size, bool periodic)
{
    std::vector<double> y_dbl(y_begin, y_end);
    return compute_moving_average_smoothing(y_dbl, window_size, periodic);
}

template<typename InputItX, typename InputItY>
std::vector<double> compute_lowess_smoothing(
    InputItX x_begin, InputItX x_end, InputItY y_begin, InputItY y_end, double bandwidth, int iterations = 3)
{
    std::vector<double> x_dbl(x_begin, x_end);
    std::vector<double> y_dbl(y_begin, y_end);
    return compute_lowess_smoothing(x_dbl, y_dbl, bandwidth, iterations);
}

template<typename InputItX, typename InputItY>
std::vector<double> compute_central_difference(
    InputItX x_begin, InputItX x_end, InputItY y_begin, InputItY y_end, bool periodic)
{
    std::vector<double> x_dbl(x_begin, x_end);
    std::vector<double> y_dbl(y_begin, y_end);
    return compute_central_difference(x_dbl, y_dbl, periodic);
}

template<typename InputItX, typename InputItY>
std::vector<std::size_t> find_knees(
    InputItX x_begin, InputItX x_end, InputItY y_begin, InputItY y_end, double sensitivity = 1.0)
{
    std::vector<double> x_dbl(x_begin, x_end);
    std::vector<double> y_dbl(y_begin, y_end);
    return find_knees(x_dbl, y_dbl, sensitivity);
}

} // namespace lsdtt
