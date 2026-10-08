// lsdtt_labelled_raster_stats_cpu.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html

#include "lsdtt_labelled_raster_stats_cpu.hpp"
#include "../utils/lsdtt_string_utils.hpp"

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "../core/lsdtt_data_loader.hpp"
#include "../utils/lsdtt_raster_geometry.hpp"
#include "../utils/lsdtt_statistical_tests.hpp"

namespace lsdtt
{
namespace
{

// -----------------------------------------------------------------------------
// validate_request
// -----------------------------------------------------------------------------
void validate_request(const labelled_raster_stats_request& request)
{
    if (request.min_values_per_label < 1)
    {
        throw std::invalid_argument(
            "lsdtt::analyse_labelled_raster_stats: min_values_per_label must be >= 1"
        );
    }

    if (!(request.ks_alpha >= 0.0 && request.ks_alpha <= 1.0))
    {
        throw std::invalid_argument(
            "lsdtt::analyse_labelled_raster_stats: ks_alpha must lie in [0, 1]"
        );
    }
}

// -----------------------------------------------------------------------------
// number_to_string
// -----------------------------------------------------------------------------
// csv_number at this table's precision (17): the same text for every finite
// value, and nan / inf / -inf spelled the same on every platform (decision D7,
// 2026-10-08; a negative NaN used to stream as "-nan").
std::string number_to_string(double v)
{
    return csv_number(v, 17);
}

// -----------------------------------------------------------------------------
// matrix_dimensions_are_valid
// -----------------------------------------------------------------------------
bool matrix_dimensions_are_valid(const labelled_raster_stats_result& result)
{
    const std::size_t n = result.records.size();

    if (result.same_distribution_matrix.size() != n)
    {
        return false;
    }

    for (const auto& row : result.same_distribution_matrix)
    {
        if (row.size() != n)
        {
            return false;
        }
    }

    return true;
}

} // unnamed namespace

// -----------------------------------------------------------------------------
// analyse_labelled_raster_stats
// -----------------------------------------------------------------------------
labelled_raster_stats_result analyse_labelled_raster_stats(
    const lsdtt_raster_cpu&          labelled_raster,
    const lsdtt_raster_cpu&          value_raster,
    const labelled_raster_stats_request& request)
{
    validate_request(request);

    const lsdtt_raster_meta& lmeta = labelled_raster.get_meta();
    const lsdtt_raster_meta& vmeta = value_raster.get_meta();

    if (!lmeta.is_valid())
    {
        throw std::runtime_error(
            "lsdtt::analyse_labelled_raster_stats: labelled_raster metadata is invalid"
        );
    }

    if (!vmeta.is_valid())
    {
        throw std::runtime_error(
            "lsdtt::analyse_labelled_raster_stats: value_raster metadata is invalid"
        );
    }

    if (!same_grid(lmeta, vmeta))
    {
        throw std::runtime_error(
            "lsdtt::analyse_labelled_raster_stats: labelled_raster and value_raster are not co-registered"
        );
    }

    const std::vector<float>&   ld = labelled_raster.get_data();
    const std::vector<uint8_t>& lm = labelled_raster.get_mask();
    const std::vector<float>&   vd = value_raster.get_data();
    const std::vector<uint8_t>& vm = value_raster.get_mask();

    if (ld.size() != lm.size() || vd.size() != vm.size() || ld.size() != vd.size())
    {
        throw std::runtime_error(
            "lsdtt::analyse_labelled_raster_stats: raster data/mask size mismatch"
        );
    }

    // label -> collected values
    std::map<std::int64_t, std::vector<double>> values_by_label;

    for (std::size_t i = 0; i < ld.size(); ++i)
    {
        if (lm[i] == 0 || vm[i] == 0)
        {
            continue;
        }

        const float label_f = ld[i];
        const float value_f = vd[i];

        if (!std::isfinite(label_f) || !std::isfinite(value_f))
        {
            continue;
        }

        if (lmeta.is_nodata(label_f) || vmeta.is_nodata(value_f))
        {
            continue;
        }

        const std::int64_t label =
            static_cast<std::int64_t>(std::llround(static_cast<double>(label_f)));

        if (label < request.minimum_label)
        {
            continue;
        }

        values_by_label[label].push_back(static_cast<double>(value_f));
    }

    labelled_raster_stats_result result;

    // Filter by minimum sample count and compute descriptive stats.
    for (const auto& kv : values_by_label)
    {
        const std::int64_t label = kv.first;
        const std::vector<double>& vals = kv.second;

        if (static_cast<int>(vals.size()) < request.min_values_per_label)
        {
            continue;
        }

        labelled_raster_stats_record rec;
        rec.label = label;
        rec.stats = compute_descriptive_stats(vals);

        result.records.push_back(std::move(rec));
    }

    if (request.verbose)
    {
        std::cout << "analyse_labelled_raster_stats: retained "
                  << result.records.size() << " label(s)\n";
    }

    // Optional K-S same-distribution matrix.
    if (request.compute_same_distribution_matrix)
    {
        const std::size_t n = result.records.size();
        result.same_distribution_matrix.assign(
            n, std::vector<bool>(n, false));

        for (std::size_t i = 0; i < n; ++i)
        {
            result.same_distribution_matrix[i][i] = true;
        }

        for (std::size_t i = 0; i < n; ++i)
        {
            const std::int64_t label_i = result.records[i].label;
            const std::vector<double>& vals_i = values_by_label.at(label_i);

            for (std::size_t j = i + 1; j < n; ++j)
            {
                const std::int64_t label_j = result.records[j].label;
                const std::vector<double>& vals_j = values_by_label.at(label_j);

                const lsdtt_ks_test_result ks =
                    compute_kolmogorov_smirnov_two_sample_test(vals_i, vals_j);

                const bool same_dist =
                    std::isfinite(ks.p_value_approx) &&
                    ks.p_value_approx >= request.ks_alpha;

                result.same_distribution_matrix[i][j] = same_dist;
                result.same_distribution_matrix[j][i] = same_dist;
            }
        }
    }

    return result;
}

// -----------------------------------------------------------------------------
// write_labelled_raster_stats_csv
// -----------------------------------------------------------------------------
void write_labelled_raster_stats_csv(
    const labelled_raster_stats_result& result,
    const std::string&                  out_path)
{
    lsdtt_data_loader loader;

    std::vector<std::string> header =
    {
        "label",
        "n_values",
        "min",
        "max",
        "range",
        "sum",
        "mean",
        "median",
        "variance_population",
        "variance_sample",
        "stddev_population",
        "stddev_sample",
        "mean_absolute_deviation",
        "median_absolute_deviation",
        "mad_normal_scaled",
        "q25",
        "q75",
        "iqr"
    };

    std::vector<std::vector<std::string>> rows;
    rows.reserve(result.records.size());

    for (const auto& rec : result.records)
    {
        const lsdtt_descriptive_stats& s = rec.stats;

        rows.push_back(
        {
            std::to_string(rec.label),
            std::to_string(s.n_values),
            number_to_string(s.min),
            number_to_string(s.max),
            number_to_string(s.range),
            number_to_string(s.sum),
            number_to_string(s.mean),
            number_to_string(s.median),
            number_to_string(s.variance_population),
            number_to_string(s.variance_sample),
            number_to_string(s.stddev_population),
            number_to_string(s.stddev_sample),
            number_to_string(s.mean_absolute_deviation),
            number_to_string(s.median_absolute_deviation),
            number_to_string(s.mad_normal_scaled),
            number_to_string(s.q25),
            number_to_string(s.q75),
            number_to_string(s.iqr)
        });
    }

    loader.write_csv(out_path, header, rows);
}

// -----------------------------------------------------------------------------
// write_labelled_raster_same_distribution_csv
// -----------------------------------------------------------------------------
void write_labelled_raster_same_distribution_csv(
    const labelled_raster_stats_result& result,
    const std::string&                  out_path)
{
    if (result.records.empty())
    {
        throw std::runtime_error(
            "lsdtt::write_labelled_raster_same_distribution_csv: no label records to write"
        );
    }

    if (!matrix_dimensions_are_valid(result))
    {
        throw std::runtime_error(
            "lsdtt::write_labelled_raster_same_distribution_csv: same_distribution_matrix dimensions do not match records"
        );
    }

    lsdtt_data_loader loader;

    std::vector<std::string> header;
    header.reserve(result.records.size() + 1);
    header.push_back("label");

    for (const auto& rec : result.records)
    {
        header.push_back(std::to_string(rec.label));
    }

    std::vector<std::vector<std::string>> rows;
    rows.reserve(result.records.size());

    for (std::size_t i = 0; i < result.records.size(); ++i)
    {
        std::vector<std::string> row;
        row.reserve(result.records.size() + 1);

        row.push_back(std::to_string(result.records[i].label));

        for (std::size_t j = 0; j < result.records.size(); ++j)
        {
            row.push_back(result.same_distribution_matrix[i][j]
                          ? "true"
                          : "false");
        }

        rows.push_back(std::move(row));
    }

    loader.write_csv(out_path, header, rows);
}

} // namespace lsdtt
