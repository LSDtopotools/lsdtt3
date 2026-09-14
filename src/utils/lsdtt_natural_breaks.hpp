// lsdtt_natural_breaks.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html

#pragma once

#include <cstdint>
#include <vector>

namespace lsdtt {

struct lsdtt_jenks_options
{
    int n_classes = 5;
    int max_histogram_bins = 2048;
};

struct lsdtt_jenks_result
{
    std::vector<double> upper_breaks;
    int requested_classes = 0;
    int effective_classes = 0;
    std::int64_t n_values = 0;
    std::int64_t n_histogram_bins_used = 0;
    double goodness_of_variance_fit = 0.0;
};

// Compute Jenks natural-break upper thresholds from numeric values.
// Non-finite values are ignored.
//
// The implementation uses a fixed-width histogram (max_histogram_bins) to keep
// runtime predictable on large rasters while still producing deterministic
// class breaks suitable for thematic mapping.
lsdtt_jenks_result compute_jenks_breaks(
    const std::vector<double>& values,
    const lsdtt_jenks_options& opts = lsdtt_jenks_options{}
);

// Classify one value into 1..K using Jenks upper breaks.
// Returns 0 when the value is outside all breaks or when breaks are empty.
int classify_jenks_value(
    double value,
    const std::vector<double>& upper_breaks
);

} // namespace lsdtt
