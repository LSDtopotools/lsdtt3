// lsdtt_surface_level_clustering_cpu.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// Do the fluvial surfaces cluster in height? (lsdtt-fluvial-surfaces, after
// the ladder has made its patches.)
//
// The values are one height per surface - each patch's median HAND. Sorted,
// they split into clusters wherever consecutive values are at least
// `required_gap` apart (the resolvability floor, kappa x E_v): a cluster is a
// set of surfaces at one level the DEM cannot tell apart.
//
// Whether that structure is more than chance: the statistic is the EMPTY
// FRACTION, the share of the value range taken by gaps of at least
// required_gap,
//     E = sum over gaps g_i >= required_gap of g_i  /  (max - min).
// Surfaces sitting on a few levels leave most of the range empty (E near 1);
// surfaces spread through the range leave little (E near 0). The null is the
// same number of values drawn uniformly over the same range, n_null times;
// p = (1 + #{E_null >= E_obs}) / (1 + n_null). A small p says the heights
// cluster more than surfaces scattered at random heights would.
//
// The draws use a generator written out here (splitmix64), so a p-value
// reproduces across platforms and standard libraries.
// -----------------------------------------------------------------------------

#pragma once

#include <cstdint>
#include <vector>

namespace lsdtt
{

struct level_cluster
{
    double lo = 0.0, hi = 0.0, centre = 0.0;   // centre = median of members
    int    n  = 0;
};

struct level_clustering_result
{
    int    n_values         = 0;
    std::vector<level_cluster> clusters;       // ascending
    double empty_fraction   = 0.0;
    double null_mean_empty_fraction = 0.0;
    double p_value          = 1.0;             // 1 when fewer than 3 values
    int    n_null           = 0;
};

// Throws std::invalid_argument when required_gap <= 0.
[[nodiscard]] level_clustering_result test_level_clustering(
    std::vector<double> values,
    double              required_gap,
    int                 n_null = 999,
    std::uint64_t       seed   = 1);

} // namespace lsdtt
