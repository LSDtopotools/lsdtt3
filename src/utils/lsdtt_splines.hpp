// lsdtt_splines.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html
//
// -----------------------------------------------------------------------------
// Parametric spline fitting for VALID ordered planar data.
//
// This module is intended for ordered curve data represented as x(s), y(s),
// where s is a strictly increasing parameter. In river applications, s will
// typically be flow distance or along-channel distance.
//
// IMPORTANT CONTRACT
// ------------------
// This module assumes that the supplied x, y, s data are already valid.
// It does NOT perform nodata filtering, raster flattening, moving-window
// extraction, or basin extraction.
//
// What it provides:
//   - natural cubic parametric spline fit to x(s), y(s)
//   - C2 continuity (twice continuously differentiable)
//   - first and second derivatives with respect to s
//   - signed curvature
//   - optional pre-smoothing of x(s), y(s) with a Gaussian kernel in s
//
// Notes:
//   - input s must be strictly increasing
//   - evaluate() clamps query points to the spline domain [s_min, s_max]
//   - curvature returned here is signed planar curvature:
//
//       kappa = (x' y'' - y' x'') / (x'^2 + y'^2)^(3/2)
//
//     If you want curvature magnitude only, use std::abs(kappa).
//
//   - if s is true arc length, then sqrt(x'^2 + y'^2) should be close to 1
// -----------------------------------------------------------------------------

#pragma once

#include <cstddef>
#include <limits>
#include <vector>

namespace lsdtt
{

struct lsdtt_parametric_xy_series
{
    std::vector<double> s;
    std::vector<double> x;
    std::vector<double> y;

    bool empty() const noexcept
    {
        return s.empty();
    }

    std::size_t size() const noexcept
    {
        return s.size();
    }
};

struct lsdtt_parametric_spline_evaluation
{
    double s       = std::numeric_limits<double>::quiet_NaN();
    double x       = std::numeric_limits<double>::quiet_NaN();
    double y       = std::numeric_limits<double>::quiet_NaN();

    double dx_ds   = std::numeric_limits<double>::quiet_NaN();
    double dy_ds   = std::numeric_limits<double>::quiet_NaN();

    double d2x_ds2 = std::numeric_limits<double>::quiet_NaN();
    double d2y_ds2 = std::numeric_limits<double>::quiet_NaN();

    double speed   = std::numeric_limits<double>::quiet_NaN();
    double curvature = std::numeric_limits<double>::quiet_NaN();
};

class lsdtt_parametric_cubic_spline
{
public:
    lsdtt_parametric_cubic_spline() = default;

    lsdtt_parametric_cubic_spline(
        const std::vector<double>& s,
        const std::vector<double>& x,
        const std::vector<double>& y
    );

    void fit(
        const std::vector<double>& s,
        const std::vector<double>& x,
        const std::vector<double>& y
    );

    bool empty() const noexcept
    {
        return s_.empty();
    }

    std::size_t n_knots() const noexcept
    {
        return s_.size();
    }

    double s_min() const noexcept;
    double s_max() const noexcept;

    lsdtt_parametric_spline_evaluation evaluate(double s_query) const;

    std::vector<lsdtt_parametric_spline_evaluation> evaluate(
        const std::vector<double>& s_queries
    ) const;

    std::vector<lsdtt_parametric_spline_evaluation> sample_uniform(
        std::size_t n_points
    ) const;

private:
    std::size_t locate_interval(double s_query) const;
    static double quiet_nan();

    // Shared spline knot vector
    std::vector<double> s_;

    // x(s) coefficients per interval:
    //   x = a + b t + c t^2 + d t^3
    std::vector<double> ax_;
    std::vector<double> bx_;
    std::vector<double> cx_;
    std::vector<double> dx_;

    // y(s) coefficients per interval:
    //   y = a + b t + c t^2 + d t^3
    std::vector<double> ay_;
    std::vector<double> by_;
    std::vector<double> cy_;
    std::vector<double> dy_;
};

// -----------------------------------------------------------------------------
// Smoothing helpers
//
// These smooth x(s) and y(s) separately using a Gaussian kernel in the
// parameter s, then return smoothed samples on the original s locations.
//
// This is useful as a first-pass denoising / trend-extraction step before
// fitting a C2 spline, but note that any purely along-s smoother will tend to
// shrink tight loops before it truly cuts across them.
// -----------------------------------------------------------------------------

lsdtt_parametric_xy_series gaussian_smooth_parametric_series(
    const std::vector<double>& s,
    const std::vector<double>& x,
    const std::vector<double>& y,
    double sigma_s,
    double truncate_sigma = 3.0
);

lsdtt_parametric_cubic_spline fit_gaussian_smoothed_parametric_cubic_spline(
    const std::vector<double>& s,
    const std::vector<double>& x,
    const std::vector<double>& y,
    double sigma_s,
    double truncate_sigma = 3.0
);

} // namespace lsdtt
