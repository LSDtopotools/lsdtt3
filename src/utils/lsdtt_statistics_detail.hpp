// lsdtt_statistics_detail.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// Internal helpers shared by the statistics implementation files
// (lsdtt_statistics.cpp, lsdtt_statistical_tests.cpp), which each had their own
// identical copy. Not part of the public statistics API.
//
// Header-only: a constant and a sort, no arithmetic that rounds.
// -----------------------------------------------------------------------------

#pragma once

#include <algorithm>
#include <limits>
#include <vector>

namespace lsdtt
{
namespace statistics_detail
{

inline double nan_value()
{
    return std::numeric_limits<double>::quiet_NaN();
}

inline std::vector<double> make_sorted_copy(const std::vector<double>& values)
{
    std::vector<double> sorted = values;
    std::sort(sorted.begin(), sorted.end());
    return sorted;
}

} // namespace statistics_detail
} // namespace lsdtt
