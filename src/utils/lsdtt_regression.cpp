// lsdtt_regression.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html
//
// -----------------------------------------------------------------------------
// Linear and orthogonal linear regression for VALID paired numeric samples.
// -----------------------------------------------------------------------------

#include "lsdtt_regression.hpp"

#include <cmath>
#include <stdexcept>
#include <vector>

namespace lsdtt
{

double predict_y(
    const lsdtt_linear_regression_result& fit,
    double x
)
{
    if (!fit.has_finite_slope)
    {
        throw std::runtime_error(
            "predict_y: regression does not have a finite slope"
        );
    }

    return fit.slope * x + fit.intercept;
}

double signed_distance_to_line(
    const lsdtt_orthogonal_linear_regression_result& fit,
    double x,
    double y
)
{
    if (!fit.has_defined_direction)
    {
        throw std::runtime_error(
            "signed_distance_to_line: regression line direction is undefined"
        );
    }

    // Since a,b are the components of a unit normal vector,
    // a*x + b*y + c is the signed perpendicular distance.
    return fit.a * x + fit.b * y + fit.c;
}

lsdtt_linear_regression_result compute_linear_regression(
    const std::vector<double>& x,
    const std::vector<double>& y
)
{
    return compute_linear_regression(
        x.begin(), x.end(),
        y.begin(), y.end()
    );
}

lsdtt_orthogonal_linear_regression_result compute_orthogonal_linear_regression(
    const std::vector<double>& x,
    const std::vector<double>& y
)
{
    return compute_orthogonal_linear_regression(
        x.begin(), x.end(),
        y.begin(), y.end()
    );
}

} // namespace lsdtt
