// lsdtt_distributed_hydrology_forcing.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html

#include "lsdtt_distributed_hydrology_forcing.hpp"

#include <cmath>
#include <stdexcept>

#include "lsdtt_csv_table.hpp"

namespace lsdtt
{

namespace
{

// What a forcing value is allowed to be. Non-finite values (NaN, inf) are
// rejected under both policies.
enum class sign_policy
{
    non_negative,  // depths: v >= 0 and finite
    signed_finite  // temperatures: any finite v
};

// The shared reader behind both public loaders.
//
//   column_primary  the column read as-is (scale 1).
//   column_mm       an alternative column in millimetres, scaled by 0.001 on
//                   read. Empty means column_primary is the only accepted one.
//
// column_primary takes precedence when both are present. The messages for the
// non-negative policy are byte-identical to the CLI's former static
// load_depth_series.
std::vector<double> load_series(const std::string& path,
                                const std::string& column_primary,
                                const std::string& column_mm,
                                const sign_policy  policy)
{
    const lsdtt_csv_table table = load_csv_table(path);

    const bool have_primary = csv_has_column(table, column_primary);
    const bool have_mm      = !column_mm.empty() && csv_has_column(table, column_mm);

    if (!have_primary && !have_mm)
    {
        if (!column_mm.empty())
        {
            throw std::runtime_error(
                "forcing file '" + path + "' has no column '" + column_primary +
                "' and no column '" + column_mm + "'. One of the two is required, "
                "holding the depth for each timestep, one row per step.");
        }
        throw std::runtime_error(
            "forcing file '" + path + "' has no column '" + column_primary +
            "'. It is required, holding the value for each timestep, one row "
            "per step.");
    }

    const std::string key   = have_primary ? column_primary : column_mm;
    const double      scale = have_primary ? 1.0 : 0.001;

    std::vector<double> series;
    series.reserve(table.rows.size());

    for (std::size_t i = 0; i < table.rows.size(); ++i)
    {
        const double v = csv_row_double_required(
            table.rows[i], key, path, static_cast<int>(i) + 1);

        if (policy == sign_policy::non_negative)
        {
            if (!(v >= 0.0) || !std::isfinite(v))
            {
                throw std::runtime_error(
                    "forcing file '" + path + "' row " + std::to_string(i + 1) +
                    ": '" + key + "' is negative or non-finite (" +
                    std::to_string(v) + "). A depth per timestep cannot be either.");
            }
        }
        else if (!std::isfinite(v))
        {
            throw std::runtime_error(
                "forcing file '" + path + "' row " + std::to_string(i + 1) +
                ": '" + key + "' is non-finite (" + std::to_string(v) +
                "). A value per timestep must be finite.");
        }
        series.push_back(v * scale);
    }

    if (series.empty())
    {
        throw std::runtime_error(
            "forcing file '" + path + "' has a header but no data rows.");
    }

    return series;
}

} // namespace

std::vector<double> load_distributed_hydrology_depth_series(
    const std::string& path,
    const std::string& column_metres,
    const std::string& column_mm)
{
    return load_series(path, column_metres, column_mm, sign_policy::non_negative);
}

std::vector<double> load_distributed_hydrology_temperature_series(
    const std::string& path,
    const std::string& column)
{
    return load_series(path, column, "", sign_policy::signed_finite);
}

} // namespace lsdtt
