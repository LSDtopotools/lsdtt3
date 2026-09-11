// lsdtt_numeric_utils.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
// Released under the GNU General Public Licence v3 (GPL-3.0)

#pragma once

#include <cmath>
#include <stdexcept>
#include <vector>

namespace lsdtt {

// Generates start, start + step, ... through the inclusive end, allowing the
// final accumulated value to lie up to half a step beyond end. This preserves
// the established sweep convention used by parameterised numeric analyses.
inline std::vector<double> inclusive_numeric_range(
    double start,
    double end,
    double step)
{
    if (!(std::isfinite(start) && std::isfinite(end) && std::isfinite(step)))
        throw std::invalid_argument("inclusive_numeric_range: values must be finite");
    if (step <= 0.0)
        throw std::invalid_argument("inclusive_numeric_range: step must be positive");
    if (end < start)
        throw std::invalid_argument("inclusive_numeric_range: end must not precede start");

    std::vector<double> values;
    for (double value = start; value <= end + 0.5 * step; value += step)
        values.push_back(value);
    return values;
}

} // namespace lsdtt
