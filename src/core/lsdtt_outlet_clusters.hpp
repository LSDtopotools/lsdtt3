// lsdtt_outlet_clusters.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html

#pragma once

#include <string>
#include <vector>

namespace lsdtt
{

enum class lsdtt_outlet_cluster_mode
{
    auto_mode,
    integer_mode,
    string_mode
};

[[nodiscard]] lsdtt_outlet_cluster_mode parse_outlet_cluster_mode(
    const std::string& mode_text);

[[nodiscard]] bool outlet_cluster_value_is_integer(
    const std::string& value);

[[nodiscard]] std::string make_unit_pair_id(
    int unit_a,
    int unit_b);

[[nodiscard]] std::vector<int> assign_unit_ids_from_keys(
    const std::vector<std::string>& unit_keys);

} // namespace lsdtt
