// lsdtt_chi_extraction_cpu.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html
//
// -----------------------------------------------------------------------------
// lsdtt_chi_extraction_cpu
//
// Chi-coordinate extraction utilities.
//
// This file implements chi extraction for:
//
//   1. whole-raster chi from all flow outlets with chi(outlet) = 0
//   2. outlet-seeded chi with each supplied outlet reset to chi = 0
//   3. outlet-seeded chi with outlet chi inherited from a precomputed chi raster
//
// Implemented domains:
//   - full upslope raster domain
//   - threshold-defined channel domain
//
// Junction-network-backed extraction is currently left as an explicit
// not-yet-implemented runtime error until the minimal required query API is
// fixed on lsdtt_junction_network_base.
// -----------------------------------------------------------------------------

#include "lsdtt_chi_extraction_cpu.hpp"

#include "lsdtt_basin_cpu.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace lsdtt
{
namespace
{

constexpr double k_sqrt_2 = 1.4142135623730950488;

// -----------------------------------------------------------------------------
// Helpers
// -----------------------------------------------------------------------------
void validate_flow_for_chi(const lsdtt_flow_base& flow)
{
    if (flow.empty())
    {
        throw std::runtime_error(
            "lsdtt_chi_extraction_cpu: flow object is empty."
        );
    }

    if (!flow.has_s_vector())
    {
        throw std::runtime_error(
            "lsdtt_chi_extraction_cpu: chi extraction requires a flow object "
            "with s_vector built."
        );
    }

    flow.meta().check_projected_metres();
}

void validate_chi_options(const lsdtt_chi_options& options)
{
    if (!std::isfinite(options.m_over_n))
    {
        throw std::invalid_argument(
            "lsdtt_chi_extraction_cpu: m_over_n must be finite."
        );
    }

    if (!std::isfinite(options.A_0) || options.A_0 <= 0.0)
    {
        throw std::invalid_argument(
            "lsdtt_chi_extraction_cpu: A_0 must be finite and > 0."
        );
    }
}

bool nearly_equal(double a, double b)
{
    const double scale =
        std::max(1.0, std::max(std::fabs(a), std::fabs(b)));

    return std::fabs(a - b) <= 1e-9 * scale;
}

void assert_raster_matches_flow_grid(
    const lsdtt_flow_base& flow,
    const lsdtt_raster_cpu& raster,
    const std::string& raster_name)
{
    const lsdtt_raster_meta& fm = flow.meta();
    const lsdtt_raster_meta& rm = raster.get_meta();

    if (fm.n_rows != rm.n_rows || fm.n_cols != rm.n_cols)
    {
        throw std::invalid_argument(
            "lsdtt_chi_extraction_cpu: " + raster_name +
            " dimensions do not match the flow grid."
        );
    }

    if (!nearly_equal(fm.cell_size, rm.cell_size) ||
        !nearly_equal(fm.x_min,     rm.x_min)     ||
        !nearly_equal(fm.y_min,     rm.y_min))
    {
        throw std::invalid_argument(
            "lsdtt_chi_extraction_cpu: " + raster_name +
            " georeferencing does not match the flow grid."
        );
    }
}

double step_length_from_receiver_dir(
    std::uint8_t dir,
    double cell_size)
{
    if (dir == lsdtt_flow_base::no_direction)
    {
        throw std::runtime_error(
            "lsdtt_chi_extraction_cpu: encountered no_direction on a "
            "non-seed chi step."
        );
    }

    return (dir % 2u == 0u) ? cell_size : cell_size * k_sqrt_2;
}

std::vector<std::int64_t> collect_flow_outlets(const lsdtt_flow_base& flow)
{
    std::vector<std::int64_t> outlets;
    outlets.reserve(flow.n_outlets());

    for (int row = 0; row < flow.n_rows(); ++row)
    {
        for (int col = 0; col < flow.n_cols(); ++col)
        {
            if (!flow.is_valid_cell(row, col))
            {
                continue;
            }

            const std::int64_t flat = flow.flat_index_unchecked(row, col);

            if (flow.receiver_of_unchecked(flat) == k_outlet_receiver)
            {
                outlets.push_back(flat);
            }
        }
    }

    if (flow.n_valid_cells() > 0 && outlets.empty())
    {
        throw std::runtime_error(
            "lsdtt_chi_extraction_cpu: no outlets found in a non-empty flow object."
        );
    }

    return outlets;
}

std::vector<std::int64_t> unique_valid_outlets(
    const lsdtt_flow_base& flow,
    const std::vector<std::int64_t>& outlet_flats)
{
    std::vector<std::int64_t> unique_outlets;
    unique_outlets.reserve(outlet_flats.size());

    std::unordered_set<std::int64_t> seen;
    seen.reserve(outlet_flats.size());

    for (const std::int64_t flat : outlet_flats)
    {
        if (!flow.is_valid_cell(flat))
        {
            throw std::invalid_argument(
                "lsdtt_chi_extraction_cpu: supplied outlet flat index is "
                "invalid or out of bounds."
            );
        }

        if (seen.insert(flat).second)
        {
            unique_outlets.push_back(flat);
        }
    }

    return unique_outlets;
}

float sample_seed_chi(
    const lsdtt_flow_base& flow,
    const lsdtt_raster_cpu& seed_chi_raster,
    std::int64_t outlet_flat)
{
    assert_raster_matches_flow_grid(flow, seed_chi_raster, "seed_chi_raster");

    const std::vector<float>&   data = seed_chi_raster.get_data();
    const std::vector<uint8_t>& mask = seed_chi_raster.get_mask();

    const std::size_t idx = static_cast<std::size_t>(outlet_flat);

    if (mask[idx] == 0)
    {
        throw std::runtime_error(
            "lsdtt_chi_extraction_cpu: seed_chi_raster is nodata at a "
            "supplied outlet."
        );
    }

    const float chi0 = data[idx];
    if (!std::isfinite(chi0))
    {
        throw std::runtime_error(
            "lsdtt_chi_extraction_cpu: seed_chi_raster has non-finite value "
            "at a supplied outlet."
        );
    }

    return chi0;
}


struct basin_full_result
{
    std::vector<std::int64_t> flats;
    std::vector<float>        chi;
};

template <class DenominatorFn>
basin_full_result compute_full_basin_chi(
    const lsdtt_flow_base& flow,
    std::int64_t outlet_flat,
    float seed_chi,
    const lsdtt_chi_options& options,
    DenominatorFn&& denominator_fn)
{
    if (!flow.is_valid_cell(outlet_flat))
    {
        throw std::invalid_argument(
            "lsdtt_chi_extraction_cpu: outlet_flat is invalid."
        );
    }

    const lsdtt_flat_index_slice slice = flow.upslope_nodes(outlet_flat);
    if (slice.empty())
    {
        throw std::runtime_error(
            "lsdtt_chi_extraction_cpu: empty upslope slice for outlet."
        );
    }

    const std::int64_t basin_start = flow.s_vector_start(outlet_flat);

    basin_full_result result;
    result.flats.assign(slice.begin(), slice.end());
    result.chi.assign(result.flats.size(), 0.0f);

    result.chi[0] = seed_chi;

    for (std::size_t i = 1; i < result.flats.size(); ++i)
    {
        const std::int64_t flat = result.flats[i];
        const std::int64_t rec  = flow.receiver_of_unchecked(flat);

        if (rec < 0)
        {
            throw std::runtime_error(
                "lsdtt_chi_extraction_cpu: encountered a negative receiver "
                "inside a non-root upslope slice."
            );
        }

        const std::int64_t rec_local =
            flow.s_vector_start_unchecked(rec) - basin_start;

        if (rec_local < 0 ||
            rec_local >= static_cast<std::int64_t>(i))
        {
            throw std::runtime_error(
                "lsdtt_chi_extraction_cpu: receiver is not in the expected "
                "upslope ordering for chi extraction."
            );
        }

        const double denom = denominator_fn(flat);
        if (!std::isfinite(denom) || denom <= 0.0)
        {
            throw std::runtime_error(
                "lsdtt_chi_extraction_cpu: chi denominator must be finite and > 0."
            );
        }

        const double step = step_length_from_receiver_dir(
            flow.receiver_dir_of_unchecked(flat),
            flow.cell_size()
        );

        result.chi[i] = static_cast<float>(
            static_cast<double>(
                result.chi[static_cast<std::size_t>(rec_local)]
            ) +
            step * std::pow(options.A_0 / denom, options.m_over_n)
        );
    }

    return result;
}

template <class DenominatorFn>
lsdtt_chi_basin_result compute_threshold_basin_chi(
    const lsdtt_flow_base& flow,
    std::int64_t outlet_flat,
    float seed_chi,
    int minimum_contributing_pixels,
    const lsdtt_chi_options& options,
    DenominatorFn&& denominator_fn)
{
    if (minimum_contributing_pixels <= 0)
    {
        throw std::invalid_argument(
            "lsdtt_chi_extraction_cpu: minimum_contributing_pixels must be > 0."
        );
    }

    if (!flow.is_valid_cell(outlet_flat))
    {
        throw std::invalid_argument(
            "lsdtt_chi_extraction_cpu: outlet_flat is invalid."
        );
    }

    if (static_cast<std::int64_t>(flow.s_vector_size(outlet_flat)) <
        static_cast<std::int64_t>(minimum_contributing_pixels))
    {
        throw std::invalid_argument(
            "lsdtt_chi_extraction_cpu: supplied outlet does not meet the "
            "minimum_contributing_pixels threshold."
        );
    }

    const lsdtt_flat_index_slice slice = flow.upslope_nodes(outlet_flat);

    lsdtt_chi_basin_result result;
    result.outlet_flat = outlet_flat;
    result.outlet_chi  = seed_chi;

    result.flats.reserve(slice.size());

    for (const std::int64_t* it = slice.begin(); it != slice.end(); ++it)
    {
        const std::int64_t flat = *it;
        if (static_cast<std::int64_t>(flow.s_vector_size_unchecked(flat)) >=
            static_cast<std::int64_t>(minimum_contributing_pixels))
        {
            result.flats.push_back(flat);
        }
    }

    if (result.flats.empty())
    {
        throw std::runtime_error(
            "lsdtt_chi_extraction_cpu: threshold selection produced no channel pixels."
        );
    }

    std::unordered_map<std::int64_t, std::size_t> local_index;
    local_index.reserve(result.flats.size());

    for (std::size_t i = 0; i < result.flats.size(); ++i)
    {
        local_index[result.flats[i]] = i;
    }

    result.chi.assign(result.flats.size(), 0.0f);
    result.receiver_flats.assign(result.flats.size(), k_nodata_receiver);

    if (result.flats[0] != outlet_flat)
    {
        throw std::runtime_error(
            "lsdtt_chi_extraction_cpu: thresholded chi extraction violated "
            "the expected root-first ordering."
        );
    }

    result.chi[0]            = seed_chi;
    result.receiver_flats[0] = k_outlet_receiver;

    for (std::size_t i = 1; i < result.flats.size(); ++i)
    {
        const std::int64_t flat = result.flats[i];
        const std::int64_t rec  = flow.receiver_of_unchecked(flat);

        const auto found = local_index.find(rec);
        if (found == local_index.end())
        {
            throw std::runtime_error(
                "lsdtt_chi_extraction_cpu: thresholded active set is missing "
                "the downstream receiver of a selected channel pixel."
            );
        }

        const std::size_t rec_local = found->second;
        if (rec_local >= i)
        {
            throw std::runtime_error(
                "lsdtt_chi_extraction_cpu: thresholded active set is not in "
                "receiver-before-donor order."
            );
        }

        const double denom = denominator_fn(flat);
        if (!std::isfinite(denom) || denom <= 0.0)
        {
            throw std::runtime_error(
                "lsdtt_chi_extraction_cpu: chi denominator must be finite and > 0."
            );
        }

        const double step = step_length_from_receiver_dir(
            flow.receiver_dir_of_unchecked(flat),
            flow.cell_size()
        );

        result.receiver_flats[i] = rec;
        result.chi[i] = static_cast<float>(
            static_cast<double>(result.chi[rec_local]) +
            step * std::pow(options.A_0 / denom, options.m_over_n)
        );
    }

    return result;
}


template <class DenominatorFn>
lsdtt_chi_basin_result compute_junction_network_basin_chi(
    const lsdtt_flow_base& flow,
    const lsdtt_junction_network_base& channel_network,
    std::int64_t outlet_flat,
    float seed_chi,
    const lsdtt_chi_options& options,
    DenominatorFn&& denominator_fn)
{
    if (!flow.is_valid_cell(outlet_flat))
    {
        throw std::invalid_argument(
            "lsdtt_chi_extraction_cpu: outlet_flat is invalid."
        );
    }

    if (!channel_network.is_channel(outlet_flat))
    {
        throw std::invalid_argument(
            "lsdtt_chi_extraction_cpu: supplied outlet_flat is not on the "
            "supplied junction network."
        );
    }

    const lsdtt_flat_index_slice slice = flow.upslope_nodes(outlet_flat);

    lsdtt_chi_basin_result result;
    result.outlet_flat = outlet_flat;
    result.outlet_chi  = seed_chi;

    result.flats.reserve(slice.size());

    for (const std::int64_t* it = slice.begin(); it != slice.end(); ++it)
    {
        const std::int64_t flat = *it;
        if (channel_network.is_channel(flat))
        {
            result.flats.push_back(flat);
        }
    }

    if (result.flats.empty())
    {
        throw std::runtime_error(
            "lsdtt_chi_extraction_cpu: junction-network selection produced "
            "no channel pixels."
        );
    }

    if (result.flats.front() != outlet_flat)
    {
        throw std::runtime_error(
            "lsdtt_chi_extraction_cpu: junction-network chi extraction "
            "violated the expected root-first ordering."
        );
    }

    std::unordered_map<std::int64_t, std::size_t> local_index;
    local_index.reserve(result.flats.size());

    for (std::size_t i = 0; i < result.flats.size(); ++i)
    {
        local_index[result.flats[i]] = i;
    }

    result.chi.assign(result.flats.size(), 0.0f);
    result.receiver_flats.assign(result.flats.size(), k_nodata_receiver);

    result.chi[0]            = seed_chi;
    result.receiver_flats[0] = k_outlet_receiver;

    for (std::size_t i = 1; i < result.flats.size(); ++i)
    {
        const std::int64_t flat = result.flats[i];
        const std::int64_t rec  = flow.receiver_of_unchecked(flat);

        const auto found = local_index.find(rec);
        if (found == local_index.end())
        {
            throw std::runtime_error(
                "lsdtt_chi_extraction_cpu: junction-network active set is "
                "missing the downstream receiver of a selected channel pixel."
            );
        }

        const std::size_t rec_local = found->second;
        if (rec_local >= i)
        {
            throw std::runtime_error(
                "lsdtt_chi_extraction_cpu: junction-network active set is not "
                "in receiver-before-donor order."
            );
        }

        const double denom = denominator_fn(flat);
        if (!std::isfinite(denom) || denom <= 0.0)
        {
            throw std::runtime_error(
                "lsdtt_chi_extraction_cpu: chi denominator must be finite and > 0."
            );
        }

        const double step = step_length_from_receiver_dir(
            flow.receiver_dir_of_unchecked(flat),
            flow.cell_size()
        );

        result.receiver_flats[i] = rec;
        result.chi[i] = static_cast<float>(
            static_cast<double>(result.chi[rec_local]) +
            step * std::pow(options.A_0 / denom, options.m_over_n)
        );
    }

    return result;
}






bool chi_values_consistent(float a, float b)
{
    const double scale =
        std::max(1.0, std::max(std::fabs(static_cast<double>(a)),
                               std::fabs(static_cast<double>(b))));

    return std::fabs(static_cast<double>(a) - static_cast<double>(b))
           <= 1e-6 * scale;
}

template <class SeedFn, class DenominatorFn>
lsdtt_raster_cpu compute_raster_from_seed_outlets(
    const lsdtt_flow_base& flow,
    const std::vector<std::int64_t>& outlet_flats,
    const lsdtt_chi_options& options,
    SeedFn&& seed_fn,
    DenominatorFn&& denominator_fn)
{
    const float ndv = static_cast<float>(flow.meta().no_data_value);

    std::vector<float>   data(flow.full_grid_size(), ndv);
    std::vector<uint8_t> mask(flow.full_grid_size(), 0u);

    const std::vector<std::int64_t> unique_outlets =
        unique_valid_outlets(flow, outlet_flats);

    for (const std::int64_t outlet_flat : unique_outlets)
    {
        const float seed_chi = seed_fn(outlet_flat);

        if (!std::isfinite(seed_chi))
        {
            throw std::runtime_error(
                "lsdtt_chi_extraction_cpu: seed chi must be finite."
            );
        }

        const basin_full_result basin = compute_full_basin_chi(
            flow,
            outlet_flat,
            seed_chi,
            options,
            denominator_fn
        );

        for (std::size_t i = 0; i < basin.flats.size(); ++i)
        {
            const std::size_t idx = static_cast<std::size_t>(basin.flats[i]);

            if (mask[idx] == 0u)
            {
                data[idx] = basin.chi[i];
                mask[idx] = 1u;
            }
            else if (!chi_values_consistent(data[idx], basin.chi[i]))
            {
                throw std::runtime_error(
                    "lsdtt_chi_extraction_cpu: overlapping outlet domains "
                    "produce inconsistent chi values in raster output. This "
                    "usually means the supplied outlets are nested but do not "
                    "share a consistent baselevel assumption for a single raster."
                );
            }
        }
    }

    return lsdtt_raster_cpu(std::move(data), std::move(mask), flow.meta());
}

template <class DenominatorFn>
lsdtt_chi_collection_result compute_threshold_collection(
    const lsdtt_flow_base& flow,
    const std::vector<std::int64_t>& outlet_flats,
    float (*seed_fn)(const lsdtt_flow_base&, const lsdtt_raster_cpu&, std::int64_t),
    const lsdtt_raster_cpu* seed_chi_raster,
    int minimum_contributing_pixels,
    const lsdtt_chi_options& options,
    DenominatorFn&& denominator_fn)
{
    lsdtt_chi_collection_result result;

    const std::vector<std::int64_t> unique_outlets =
        unique_valid_outlets(flow, outlet_flats);

    result.basins.reserve(unique_outlets.size());

    for (const std::int64_t outlet_flat : unique_outlets)
    {
        float chi0 = 0.0f;
        if (seed_fn != nullptr && seed_chi_raster != nullptr)
        {
            chi0 = seed_fn(flow, *seed_chi_raster, outlet_flat);
        }

        result.basins.push_back(
            compute_threshold_basin_chi(
                flow,
                outlet_flat,
                chi0,
                minimum_contributing_pixels,
                options,
                denominator_fn
            )
        );
    }

    return result;
}

} // anonymous namespace

// =============================================================================
// 1) Full-raster chi
// =============================================================================
lsdtt_raster_cpu get_chi_raster(
    const lsdtt_flow_base& flow,
    const lsdtt_chi_options& options)
{
    validate_flow_for_chi(flow);
    validate_chi_options(options);

    const double pixel_area = flow.cell_size() * flow.cell_size();

    auto denominator_fn = [&](std::int64_t flat) -> double
    {
        return static_cast<double>(flow.s_vector_size_unchecked(flat)) * pixel_area;
    };

    auto seed_fn = [](std::int64_t) -> float
    {
        return 0.0f;
    };

    return compute_raster_from_seed_outlets(
        flow,
        collect_flow_outlets(flow),
        options,
        seed_fn,
        denominator_fn
    );
}

lsdtt_raster_cpu get_chi_raster(
    const lsdtt_flow_base& flow,
    const lsdtt_raster_cpu& denominator_raster,
    const lsdtt_chi_options& options)
{
    validate_flow_for_chi(flow);
    validate_chi_options(options);
    assert_raster_matches_flow_grid(flow, denominator_raster, "denominator_raster");

    const std::vector<float>&   data = denominator_raster.get_data();
    const std::vector<uint8_t>& mask = denominator_raster.get_mask();

    auto denominator_fn = [&](std::int64_t flat) -> double
    {
        const std::size_t idx = static_cast<std::size_t>(flat);

        if (mask[idx] == 0)
        {
            throw std::runtime_error(
                "lsdtt_chi_extraction_cpu: denominator_raster is nodata at a "
                "pixel where chi was requested."
            );
        }

        return static_cast<double>(data[idx]);
    };

    auto seed_fn = [](std::int64_t) -> float
    {
        return 0.0f;
    };

    return compute_raster_from_seed_outlets(
        flow,
        collect_flow_outlets(flow),
        options,
        seed_fn,
        denominator_fn
    );
}

// =============================================================================
// 2a) Outlet-seeded chi on full upslope raster, zero base
// =============================================================================
lsdtt_raster_cpu get_chi_raster_from_outlets_zero_base(
    const lsdtt_flow_base& flow,
    const std::vector<std::int64_t>& outlet_flats,
    const lsdtt_chi_options& options)
{
    validate_flow_for_chi(flow);
    validate_chi_options(options);

    const double pixel_area = flow.cell_size() * flow.cell_size();

    auto denominator_fn = [&](std::int64_t flat) -> double
    {
        return static_cast<double>(flow.s_vector_size_unchecked(flat)) * pixel_area;
    };

    auto seed_fn = [](std::int64_t) -> float
    {
        return 0.0f;
    };

    return compute_raster_from_seed_outlets(
        flow,
        outlet_flats,
        options,
        seed_fn,
        denominator_fn
    );
}

lsdtt_raster_cpu get_chi_raster_from_outlets_zero_base(
    const lsdtt_flow_base& flow,
    const std::vector<std::int64_t>& outlet_flats,
    const lsdtt_raster_cpu& denominator_raster,
    const lsdtt_chi_options& options)
{
    validate_flow_for_chi(flow);
    validate_chi_options(options);
    assert_raster_matches_flow_grid(flow, denominator_raster, "denominator_raster");

    const std::vector<float>&   data = denominator_raster.get_data();
    const std::vector<uint8_t>& mask = denominator_raster.get_mask();

    auto denominator_fn = [&](std::int64_t flat) -> double
    {
        const std::size_t idx = static_cast<std::size_t>(flat);

        if (mask[idx] == 0)
        {
            throw std::runtime_error(
                "lsdtt_chi_extraction_cpu: denominator_raster is nodata at a "
                "pixel where chi was requested."
            );
        }

        return static_cast<double>(data[idx]);
    };

    auto seed_fn = [](std::int64_t) -> float
    {
        return 0.0f;
    };

    return compute_raster_from_seed_outlets(
        flow,
        outlet_flats,
        options,
        seed_fn,
        denominator_fn
    );
}

// =============================================================================
// 2b) Threshold-defined channel domain, zero base
// =============================================================================
lsdtt_chi_basin_result get_chi_from_outlet_zero_base_threshold(
    const lsdtt_flow_base& flow,
    std::int64_t outlet_flat,
    int minimum_contributing_pixels,
    const lsdtt_chi_options& options)
{
    validate_flow_for_chi(flow);
    validate_chi_options(options);

    const double pixel_area = flow.cell_size() * flow.cell_size();

    auto denominator_fn = [&](std::int64_t flat) -> double
    {
        return static_cast<double>(flow.s_vector_size_unchecked(flat)) * pixel_area;
    };

    return compute_threshold_basin_chi(
        flow,
        outlet_flat,
        0.0f,
        minimum_contributing_pixels,
        options,
        denominator_fn
    );
}

lsdtt_chi_basin_result get_chi_from_outlet_zero_base_threshold(
    const lsdtt_flow_base& flow,
    std::int64_t outlet_flat,
    int minimum_contributing_pixels,
    const lsdtt_raster_cpu& denominator_raster,
    const lsdtt_chi_options& options)
{
    validate_flow_for_chi(flow);
    validate_chi_options(options);
    assert_raster_matches_flow_grid(flow, denominator_raster, "denominator_raster");

    const std::vector<float>&   data = denominator_raster.get_data();
    const std::vector<uint8_t>& mask = denominator_raster.get_mask();

    auto denominator_fn = [&](std::int64_t flat) -> double
    {
        const std::size_t idx = static_cast<std::size_t>(flat);

        if (mask[idx] == 0)
        {
            throw std::runtime_error(
                "lsdtt_chi_extraction_cpu: denominator_raster is nodata at a "
                "selected threshold-channel pixel."
            );
        }

        return static_cast<double>(data[idx]);
    };

    return compute_threshold_basin_chi(
        flow,
        outlet_flat,
        0.0f,
        minimum_contributing_pixels,
        options,
        denominator_fn
    );
}

lsdtt_chi_collection_result get_chi_from_outlets_zero_base_threshold(
    const lsdtt_flow_base& flow,
    const std::vector<std::int64_t>& outlet_flats,
    int minimum_contributing_pixels,
    const lsdtt_chi_options& options)
{
    validate_flow_for_chi(flow);
    validate_chi_options(options);

    const double pixel_area = flow.cell_size() * flow.cell_size();

    auto denominator_fn = [&](std::int64_t flat) -> double
    {
        return static_cast<double>(flow.s_vector_size_unchecked(flat)) * pixel_area;
    };

    return compute_threshold_collection(
        flow,
        outlet_flats,
        nullptr,
        nullptr,
        minimum_contributing_pixels,
        options,
        denominator_fn
    );
}

lsdtt_chi_collection_result get_chi_from_outlets_zero_base_threshold(
    const lsdtt_flow_base& flow,
    const std::vector<std::int64_t>& outlet_flats,
    int minimum_contributing_pixels,
    const lsdtt_raster_cpu& denominator_raster,
    const lsdtt_chi_options& options)
{
    validate_flow_for_chi(flow);
    validate_chi_options(options);
    assert_raster_matches_flow_grid(flow, denominator_raster, "denominator_raster");

    const std::vector<float>&   data = denominator_raster.get_data();
    const std::vector<uint8_t>& mask = denominator_raster.get_mask();

    auto denominator_fn = [&](std::int64_t flat) -> double
    {
        const std::size_t idx = static_cast<std::size_t>(flat);

        if (mask[idx] == 0)
        {
            throw std::runtime_error(
                "lsdtt_chi_extraction_cpu: denominator_raster is nodata at a "
                "selected threshold-channel pixel."
            );
        }

        return static_cast<double>(data[idx]);
    };

    return compute_threshold_collection(
        flow,
        outlet_flats,
        nullptr,
        nullptr,
        minimum_contributing_pixels,
        options,
        denominator_fn
    );
}

// =============================================================================
// 2c) Junction-network domain, zero base
// =============================================================================
lsdtt_chi_basin_result get_chi_from_outlet_zero_base_junction_network(
    const lsdtt_flow_base& flow,
    const lsdtt_junction_network_base& channel_network,
    std::int64_t outlet_flat,
    const lsdtt_chi_options& options)
{
    validate_flow_for_chi(flow);
    validate_chi_options(options);

    const double pixel_area = flow.cell_size() * flow.cell_size();

    auto denominator_fn = [&](std::int64_t flat) -> double
    {
        return static_cast<double>(flow.s_vector_size_unchecked(flat)) * pixel_area;
    };

    return compute_junction_network_basin_chi(
        flow,
        channel_network,
        outlet_flat,
        0.0f,
        options,
        denominator_fn
    );
}

lsdtt_chi_basin_result get_chi_from_outlet_zero_base_junction_network(
    const lsdtt_flow_base& flow,
    const lsdtt_junction_network_base& channel_network,
    std::int64_t outlet_flat,
    const lsdtt_raster_cpu& denominator_raster,
    const lsdtt_chi_options& options)
{
    validate_flow_for_chi(flow);
    validate_chi_options(options);
    assert_raster_matches_flow_grid(flow, denominator_raster, "denominator_raster");

    const std::vector<float>&   data = denominator_raster.get_data();
    const std::vector<uint8_t>& mask = denominator_raster.get_mask();

    auto denominator_fn = [&](std::int64_t flat) -> double
    {
        const std::size_t idx = static_cast<std::size_t>(flat);

        if (mask[idx] == 0)
        {
            throw std::runtime_error(
                "lsdtt_chi_extraction_cpu: denominator_raster is nodata at a "
                "selected junction-network channel pixel."
            );
        }

        return static_cast<double>(data[idx]);
    };

    return compute_junction_network_basin_chi(
        flow,
        channel_network,
        outlet_flat,
        0.0f,
        options,
        denominator_fn
    );
}

lsdtt_chi_collection_result get_chi_from_outlets_zero_base_junction_network(
    const lsdtt_flow_base& flow,
    const lsdtt_junction_network_base& channel_network,
    const std::vector<std::int64_t>& outlet_flats,
    const lsdtt_chi_options& options)
{
    validate_flow_for_chi(flow);
    validate_chi_options(options);

    const double pixel_area = flow.cell_size() * flow.cell_size();

    auto denominator_fn = [&](std::int64_t flat) -> double
    {
        return static_cast<double>(flow.s_vector_size_unchecked(flat)) * pixel_area;
    };

    lsdtt_chi_collection_result result;

    const std::vector<std::int64_t> unique_outlets =
        unique_valid_outlets(flow, outlet_flats);

    result.basins.reserve(unique_outlets.size());

    for (const std::int64_t outlet_flat : unique_outlets)
    {
        result.basins.push_back(
            compute_junction_network_basin_chi(
                flow,
                channel_network,
                outlet_flat,
                0.0f,
                options,
                denominator_fn
            )
        );
    }

    return result;
}

lsdtt_chi_collection_result get_chi_from_outlets_zero_base_junction_network(
    const lsdtt_flow_base& flow,
    const lsdtt_junction_network_base& channel_network,
    const std::vector<std::int64_t>& outlet_flats,
    const lsdtt_raster_cpu& denominator_raster,
    const lsdtt_chi_options& options)
{
    validate_flow_for_chi(flow);
    validate_chi_options(options);
    assert_raster_matches_flow_grid(flow, denominator_raster, "denominator_raster");

    const std::vector<float>&   data = denominator_raster.get_data();
    const std::vector<uint8_t>& mask = denominator_raster.get_mask();

    auto denominator_fn = [&](std::int64_t flat) -> double
    {
        const std::size_t idx = static_cast<std::size_t>(flat);

        if (mask[idx] == 0)
        {
            throw std::runtime_error(
                "lsdtt_chi_extraction_cpu: denominator_raster is nodata at a "
                "selected junction-network channel pixel."
            );
        }

        return static_cast<double>(data[idx]);
    };

    lsdtt_chi_collection_result result;

    const std::vector<std::int64_t> unique_outlets =
        unique_valid_outlets(flow, outlet_flats);

    result.basins.reserve(unique_outlets.size());

    for (const std::int64_t outlet_flat : unique_outlets)
    {
        result.basins.push_back(
            compute_junction_network_basin_chi(
                flow,
                channel_network,
                outlet_flat,
                0.0f,
                options,
                denominator_fn
            )
        );
    }

    return result;
}


// =============================================================================
// 3a) Outlet-seeded chi on full upslope raster, inherited base
// =============================================================================
lsdtt_raster_cpu get_chi_raster_from_outlets_inherited_base(
    const lsdtt_flow_base& flow,
    const std::vector<std::int64_t>& outlet_flats,
    const lsdtt_raster_cpu& seed_chi_raster,
    const lsdtt_chi_options& options)
{
    validate_flow_for_chi(flow);
    validate_chi_options(options);
    assert_raster_matches_flow_grid(flow, seed_chi_raster, "seed_chi_raster");

    const double pixel_area = flow.cell_size() * flow.cell_size();

    auto denominator_fn = [&](std::int64_t flat) -> double
    {
        return static_cast<double>(flow.s_vector_size_unchecked(flat)) * pixel_area;
    };

    auto seed_fn = [&](std::int64_t outlet_flat) -> float
    {
        return sample_seed_chi(flow, seed_chi_raster, outlet_flat);
    };

    return compute_raster_from_seed_outlets(
        flow,
        outlet_flats,
        options,
        seed_fn,
        denominator_fn
    );
}

lsdtt_raster_cpu get_chi_raster_from_outlets_inherited_base(
    const lsdtt_flow_base& flow,
    const std::vector<std::int64_t>& outlet_flats,
    const lsdtt_raster_cpu& seed_chi_raster,
    const lsdtt_raster_cpu& denominator_raster,
    const lsdtt_chi_options& options)
{
    validate_flow_for_chi(flow);
    validate_chi_options(options);
    assert_raster_matches_flow_grid(flow, seed_chi_raster, "seed_chi_raster");
    assert_raster_matches_flow_grid(flow, denominator_raster, "denominator_raster");

    const std::vector<float>&   data = denominator_raster.get_data();
    const std::vector<uint8_t>& mask = denominator_raster.get_mask();

    auto denominator_fn = [&](std::int64_t flat) -> double
    {
        const std::size_t idx = static_cast<std::size_t>(flat);

        if (mask[idx] == 0)
        {
            throw std::runtime_error(
                "lsdtt_chi_extraction_cpu: denominator_raster is nodata at a "
                "pixel where chi was requested."
            );
        }

        return static_cast<double>(data[idx]);
    };

    auto seed_fn = [&](std::int64_t outlet_flat) -> float
    {
        return sample_seed_chi(flow, seed_chi_raster, outlet_flat);
    };

    return compute_raster_from_seed_outlets(
        flow,
        outlet_flats,
        options,
        seed_fn,
        denominator_fn
    );
}

// =============================================================================
// 3b) Threshold-defined channel domain, inherited base
// =============================================================================
lsdtt_chi_basin_result get_chi_from_outlet_inherited_base_threshold(
    const lsdtt_flow_base& flow,
    std::int64_t outlet_flat,
    const lsdtt_raster_cpu& seed_chi_raster,
    int minimum_contributing_pixels,
    const lsdtt_chi_options& options)
{
    validate_flow_for_chi(flow);
    validate_chi_options(options);
    assert_raster_matches_flow_grid(flow, seed_chi_raster, "seed_chi_raster");

    const double pixel_area = flow.cell_size() * flow.cell_size();

    auto denominator_fn = [&](std::int64_t flat) -> double
    {
        return static_cast<double>(flow.s_vector_size_unchecked(flat)) * pixel_area;
    };

    return compute_threshold_basin_chi(
        flow,
        outlet_flat,
        sample_seed_chi(flow, seed_chi_raster, outlet_flat),
        minimum_contributing_pixels,
        options,
        denominator_fn
    );
}

lsdtt_chi_basin_result get_chi_from_outlet_inherited_base_threshold(
    const lsdtt_flow_base& flow,
    std::int64_t outlet_flat,
    const lsdtt_raster_cpu& seed_chi_raster,
    int minimum_contributing_pixels,
    const lsdtt_raster_cpu& denominator_raster,
    const lsdtt_chi_options& options)
{
    validate_flow_for_chi(flow);
    validate_chi_options(options);
    assert_raster_matches_flow_grid(flow, seed_chi_raster, "seed_chi_raster");
    assert_raster_matches_flow_grid(flow, denominator_raster, "denominator_raster");

    const std::vector<float>&   data = denominator_raster.get_data();
    const std::vector<uint8_t>& mask = denominator_raster.get_mask();

    auto denominator_fn = [&](std::int64_t flat) -> double
    {
        const std::size_t idx = static_cast<std::size_t>(flat);

        if (mask[idx] == 0)
        {
            throw std::runtime_error(
                "lsdtt_chi_extraction_cpu: denominator_raster is nodata at a "
                "selected threshold-channel pixel."
            );
        }

        return static_cast<double>(data[idx]);
    };

    return compute_threshold_basin_chi(
        flow,
        outlet_flat,
        sample_seed_chi(flow, seed_chi_raster, outlet_flat),
        minimum_contributing_pixels,
        options,
        denominator_fn
    );
}

lsdtt_chi_collection_result get_chi_from_outlets_inherited_base_threshold(
    const lsdtt_flow_base& flow,
    const std::vector<std::int64_t>& outlet_flats,
    const lsdtt_raster_cpu& seed_chi_raster,
    int minimum_contributing_pixels,
    const lsdtt_chi_options& options)
{
    validate_flow_for_chi(flow);
    validate_chi_options(options);
    assert_raster_matches_flow_grid(flow, seed_chi_raster, "seed_chi_raster");

    const double pixel_area = flow.cell_size() * flow.cell_size();

    auto denominator_fn = [&](std::int64_t flat) -> double
    {
        return static_cast<double>(flow.s_vector_size_unchecked(flat)) * pixel_area;
    };

    return compute_threshold_collection(
        flow,
        outlet_flats,
        sample_seed_chi,
        &seed_chi_raster,
        minimum_contributing_pixels,
        options,
        denominator_fn
    );
}

lsdtt_chi_collection_result get_chi_from_outlets_inherited_base_threshold(
    const lsdtt_flow_base& flow,
    const std::vector<std::int64_t>& outlet_flats,
    const lsdtt_raster_cpu& seed_chi_raster,
    int minimum_contributing_pixels,
    const lsdtt_raster_cpu& denominator_raster,
    const lsdtt_chi_options& options)
{
    validate_flow_for_chi(flow);
    validate_chi_options(options);
    assert_raster_matches_flow_grid(flow, seed_chi_raster, "seed_chi_raster");
    assert_raster_matches_flow_grid(flow, denominator_raster, "denominator_raster");

    const std::vector<float>&   data = denominator_raster.get_data();
    const std::vector<uint8_t>& mask = denominator_raster.get_mask();

    auto denominator_fn = [&](std::int64_t flat) -> double
    {
        const std::size_t idx = static_cast<std::size_t>(flat);

        if (mask[idx] == 0)
        {
            throw std::runtime_error(
                "lsdtt_chi_extraction_cpu: denominator_raster is nodata at a "
                "selected threshold-channel pixel."
            );
        }

        return static_cast<double>(data[idx]);
    };

    return compute_threshold_collection(
        flow,
        outlet_flats,
        sample_seed_chi,
        &seed_chi_raster,
        minimum_contributing_pixels,
        options,
        denominator_fn
    );
}

// =============================================================================
// 3c) Junction-network domain, inherited base
// =============================================================================
lsdtt_chi_basin_result get_chi_from_outlet_inherited_base_junction_network(
    const lsdtt_flow_base& flow,
    const lsdtt_junction_network_base& channel_network,
    std::int64_t outlet_flat,
    const lsdtt_raster_cpu& seed_chi_raster,
    const lsdtt_chi_options& options)
{
    validate_flow_for_chi(flow);
    validate_chi_options(options);
    assert_raster_matches_flow_grid(flow, seed_chi_raster, "seed_chi_raster");

    const double pixel_area = flow.cell_size() * flow.cell_size();

    auto denominator_fn = [&](std::int64_t flat) -> double
    {
        return static_cast<double>(flow.s_vector_size_unchecked(flat)) * pixel_area;
    };

    return compute_junction_network_basin_chi(
        flow,
        channel_network,
        outlet_flat,
        sample_seed_chi(flow, seed_chi_raster, outlet_flat),
        options,
        denominator_fn
    );
}

lsdtt_chi_basin_result get_chi_from_outlet_inherited_base_junction_network(
    const lsdtt_flow_base& flow,
    const lsdtt_junction_network_base& channel_network,
    std::int64_t outlet_flat,
    const lsdtt_raster_cpu& seed_chi_raster,
    const lsdtt_raster_cpu& denominator_raster,
    const lsdtt_chi_options& options)
{
    validate_flow_for_chi(flow);
    validate_chi_options(options);
    assert_raster_matches_flow_grid(flow, seed_chi_raster, "seed_chi_raster");
    assert_raster_matches_flow_grid(flow, denominator_raster, "denominator_raster");

    const std::vector<float>&   data = denominator_raster.get_data();
    const std::vector<uint8_t>& mask = denominator_raster.get_mask();

    auto denominator_fn = [&](std::int64_t flat) -> double
    {
        const std::size_t idx = static_cast<std::size_t>(flat);

        if (mask[idx] == 0)
        {
            throw std::runtime_error(
                "lsdtt_chi_extraction_cpu: denominator_raster is nodata at a "
                "selected junction-network channel pixel."
            );
        }

        return static_cast<double>(data[idx]);
    };

    return compute_junction_network_basin_chi(
        flow,
        channel_network,
        outlet_flat,
        sample_seed_chi(flow, seed_chi_raster, outlet_flat),
        options,
        denominator_fn
    );
}

lsdtt_chi_collection_result get_chi_from_outlets_inherited_base_junction_network(
    const lsdtt_flow_base& flow,
    const lsdtt_junction_network_base& channel_network,
    const std::vector<std::int64_t>& outlet_flats,
    const lsdtt_raster_cpu& seed_chi_raster,
    const lsdtt_chi_options& options)
{
    validate_flow_for_chi(flow);
    validate_chi_options(options);
    assert_raster_matches_flow_grid(flow, seed_chi_raster, "seed_chi_raster");

    const double pixel_area = flow.cell_size() * flow.cell_size();

    auto denominator_fn = [&](std::int64_t flat) -> double
    {
        return static_cast<double>(flow.s_vector_size_unchecked(flat)) * pixel_area;
    };

    lsdtt_chi_collection_result result;

    const std::vector<std::int64_t> unique_outlets =
        unique_valid_outlets(flow, outlet_flats);

    result.basins.reserve(unique_outlets.size());

    for (const std::int64_t outlet_flat : unique_outlets)
    {
        result.basins.push_back(
            compute_junction_network_basin_chi(
                flow,
                channel_network,
                outlet_flat,
                sample_seed_chi(flow, seed_chi_raster, outlet_flat),
                options,
                denominator_fn
            )
        );
    }

    return result;
}

lsdtt_chi_collection_result get_chi_from_outlets_inherited_base_junction_network(
    const lsdtt_flow_base& flow,
    const lsdtt_junction_network_base& channel_network,
    const std::vector<std::int64_t>& outlet_flats,
    const lsdtt_raster_cpu& seed_chi_raster,
    const lsdtt_raster_cpu& denominator_raster,
    const lsdtt_chi_options& options)
{
    validate_flow_for_chi(flow);
    validate_chi_options(options);
    assert_raster_matches_flow_grid(flow, seed_chi_raster, "seed_chi_raster");
    assert_raster_matches_flow_grid(flow, denominator_raster, "denominator_raster");

    const std::vector<float>&   data = denominator_raster.get_data();
    const std::vector<uint8_t>& mask = denominator_raster.get_mask();

    auto denominator_fn = [&](std::int64_t flat) -> double
    {
        const std::size_t idx = static_cast<std::size_t>(flat);

        if (mask[idx] == 0)
        {
            throw std::runtime_error(
                "lsdtt_chi_extraction_cpu: denominator_raster is nodata at a "
                "selected junction-network channel pixel."
            );
        }

        return static_cast<double>(data[idx]);
    };

    lsdtt_chi_collection_result result;

    const std::vector<std::int64_t> unique_outlets =
        unique_valid_outlets(flow, outlet_flats);

    result.basins.reserve(unique_outlets.size());

    for (const std::int64_t outlet_flat : unique_outlets)
    {
        result.basins.push_back(
            compute_junction_network_basin_chi(
                flow,
                channel_network,
                outlet_flat,
                sample_seed_chi(flow, seed_chi_raster, outlet_flat),
                options,
                denominator_fn
            )
        );
    }

    return result;
}

// =============================================================================
// reject_chi_basins_with_few_channel_pixels
// =============================================================================
std::size_t reject_chi_basins_with_few_channel_pixels(
    lsdtt_chi_collection_result& chi,
    std::size_t                  min_channel_pixels
)
{
    if (min_channel_pixels <= 1)
    {
        return 0;
    }

    const std::size_t n_before = chi.basins.size();

    chi.basins.erase(
        std::remove_if(
            chi.basins.begin(),
            chi.basins.end(),
            [min_channel_pixels](const lsdtt_chi_basin_result& basin)
            {
                return basin.flats.size() < min_channel_pixels;
            }
        ),
        chi.basins.end()
    );

    return n_before - chi.basins.size();
}

// =============================================================================
// keep_n_largest_complete_chi_basins
// =============================================================================
std::size_t keep_n_largest_complete_chi_basins(
    lsdtt_chi_collection_result& chi,
    const lsdtt_flow_base&       flow,
    std::size_t                  n_largest
)
{
    if (n_largest == 0)
    {
        return 0;  // 0 == keep all complete basins (no size cap).
    }

    const std::size_t n_before = chi.basins.size();

    // Rank complete basins by total basin pixels (drainage area proxy).
    struct ranked_basin
    {
        std::size_t  original_index;
        std::int32_t n_basin_pixels;
    };

    std::vector<ranked_basin> complete;
    complete.reserve(chi.basins.size());

    for (std::size_t bi = 0; bi < chi.basins.size(); ++bi)
    {
        const lsdtt_chi_basin_result& basin = chi.basins[bi];
        if (basin.flats.empty())
        {
            continue;
        }

        lsdtt_basin_cpu basin_obj(basin.outlet_flat, flow);
        if (!basin_obj.check_basin_completeness())
        {
            continue;
        }

        complete.push_back({bi, basin_obj.get_number_of_basin_pixels()});
    }

    // Largest first.
    std::sort(complete.begin(), complete.end(),
        [](const ranked_basin& a, const ranked_basin& b)
        {
            return a.n_basin_pixels > b.n_basin_pixels;
        });

    if (complete.size() > n_largest)
    {
        complete.resize(n_largest);
    }

    // Restore original ordering so basin_key / basin_index stay stable.
    std::sort(complete.begin(), complete.end(),
        [](const ranked_basin& a, const ranked_basin& b)
        {
            return a.original_index < b.original_index;
        });

    std::vector<lsdtt_chi_basin_result> kept;
    kept.reserve(complete.size());
    for (const ranked_basin& r : complete)
    {
        kept.push_back(std::move(chi.basins[r.original_index]));
    }

    chi.basins = std::move(kept);

    return n_before - chi.basins.size();
}


} // namespace lsdtt
