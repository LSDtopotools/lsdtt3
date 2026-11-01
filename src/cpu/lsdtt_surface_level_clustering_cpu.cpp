// lsdtt_surface_level_clustering_cpu.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html

#include "lsdtt_surface_level_clustering_cpu.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace lsdtt
{

namespace
{

std::uint64_t splitmix64(std::uint64_t& state)
{
    std::uint64_t z = (state += 0x9E3779B97F4A7C15ull);
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ull;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBull;
    return z ^ (z >> 31);
}

// Uniform on [0, 1) from the top 53 bits.
double unit_draw(std::uint64_t& state)
{
    return static_cast<double>(splitmix64(state) >> 11) * 0x1.0p-53;
}

// Sorted values -> empty fraction.
double empty_fraction_sorted(const std::vector<double>& v, double gap)
{
    const double range = v.back() - v.front();
    if (!(range > 0.0)) return 0.0;
    double empty = 0.0;
    for (std::size_t i = 1; i < v.size(); ++i)
    {
        const double g = v[i] - v[i - 1];
        if (g >= gap) empty += g;
    }
    return empty / range;
}

} // namespace

level_clustering_result test_level_clustering(std::vector<double> values,
                                              double              required_gap,
                                              int                 n_null,
                                              std::uint64_t       seed)
{
    if (!(required_gap > 0.0))
        throw std::invalid_argument("test_level_clustering: required_gap must be > 0");
    values.erase(std::remove_if(values.begin(), values.end(),
                                [](double v) { return !std::isfinite(v); }),
                 values.end());
    std::sort(values.begin(), values.end());

    level_clustering_result r;
    r.n_values = static_cast<int>(values.size());
    if (values.empty()) return r;

    std::size_t start = 0;
    for (std::size_t i = 1; i <= values.size(); ++i)
    {
        if (i < values.size() && values[i] - values[i - 1] < required_gap) continue;
        level_cluster c;
        c.lo = values[start];
        c.hi = values[i - 1];
        c.n  = static_cast<int>(i - start);
        const std::size_t m = start + (i - start) / 2;
        c.centre = ((i - start) % 2 == 1) ? values[m]
                                          : 0.5 * (values[m - 1] + values[m]);
        r.clusters.push_back(c);
        start = i;
    }

    if (values.size() < 3) return r;
    r.empty_fraction = empty_fraction_sorted(values, required_gap);

    const double lo = values.front(), range = values.back() - values.front();
    std::uint64_t state = seed;
    std::vector<double> draw(values.size());
    int n_ge = 0;
    double sum = 0.0;
    for (int k = 0; k < n_null; ++k)
    {
        // The two ends are fixed at the observed extremes so the null spans
        // exactly the same range; the rest fall uniformly inside it.
        draw.front() = lo;
        draw.back()  = lo + range;
        for (std::size_t i = 1; i + 1 < draw.size(); ++i)
            draw[i] = lo + range * unit_draw(state);
        std::sort(draw.begin(), draw.end());
        const double e = empty_fraction_sorted(draw, required_gap);
        sum += e;
        if (e >= r.empty_fraction) ++n_ge;
    }
    r.n_null = n_null;
    if (n_null > 0)
    {
        r.null_mean_empty_fraction = sum / n_null;
        r.p_value = static_cast<double>(1 + n_ge) / static_cast<double>(1 + n_null);
    }
    return r;
}

} // namespace lsdtt
