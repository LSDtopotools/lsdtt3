// lsdtt_surface_method_comparison_cpu.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// Comparing the floodplain / terrace classes of several methods on one grid
// (fluvial-surface methods plan, phase 7b). lsdtt-fluvial-surfaces runs every
// method listed in fluvial_surface_algorithm_choice and hands their class
// rasters here; the comparison is a pure function of those rasters, so the
// later compare-earlier-runs mode (7e) and the bindings reuse it.
//
// Classes are surface_class values: 0 none, 1 floodplain, 2 terrace. Cells
// outside `valid` take no part.
//
//   per pair:    the 3 x 3 confusion matrix (rows method a, columns b),
//                floodplain and terrace intersection-over-union, and the
//                fraction of cells classified by either method on which the
//                two agree;
//   agreement:   per cell, bit k set where method k put the given class.
// -----------------------------------------------------------------------------

#pragma once

#include <array>
#include <cstdint>
#include <string>
#include <vector>

namespace lsdtt
{

struct method_class_map
{
    std::string               name;
    std::vector<std::uint8_t> cls;   // per cell: 0 none, 1 floodplain, 2 terrace
};

struct method_pair_comparison
{
    std::size_t a = 0, b = 0;                    // indices into the methods
    std::array<std::array<std::int64_t, 3>, 3> confusion{};   // [cls_a][cls_b]
    double floodplain_iou = 0.0;   // NaN when neither method has the class
    double terrace_iou    = 0.0;
    // Of the cells either method classed floodplain or terrace, the fraction
    // given the same class by both. NaN when neither classed any.
    double classified_agreement = 0.0;
};

// Every pair (a < b). Throws std::invalid_argument on size mismatches or a
// class value above 2.
[[nodiscard]] std::vector<method_pair_comparison> compare_surface_classes(
    const std::vector<method_class_map>& methods,
    const std::vector<std::uint8_t>&     valid);

// Bit k (method k, list order) set where that method gave class `cls`; 0
// elsewhere. At most 8 methods.
[[nodiscard]] std::vector<std::uint8_t> surface_class_agreement_bits(
    const std::vector<method_class_map>& methods,
    const std::vector<std::uint8_t>&     valid,
    std::uint8_t                         cls);

} // namespace lsdtt
