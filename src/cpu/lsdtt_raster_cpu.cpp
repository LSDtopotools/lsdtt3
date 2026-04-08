// lsdtt_raster_cpu.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2025 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html

#include "lsdtt_raster_cpu.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace lsdtt {

// -----------------------------------------------------------------------------
// Constructors
// -----------------------------------------------------------------------------

lsdtt_raster_cpu::lsdtt_raster_cpu(const lsdtt_data_loader& loader)
{
    lsdtt_raster_tile tile = loader.read_full();

    set_from_vectors(
        std::move(tile.data),
        std::move(tile.mask),
        tile.meta
    );
}

lsdtt_raster_cpu::lsdtt_raster_cpu(const lsdtt_raster_tile& tile)
    : meta_(tile.meta)
{
    data_ = lsdtt_array2d_f(
        tile.meta.n_rows,
        tile.meta.n_cols,
        lsdtt_array2d_f::from_ptr,
        tile.data.data()
    );

    mask_ = lsdtt_array2d_u8(
        tile.meta.n_rows,
        tile.meta.n_cols,
        lsdtt_array2d_u8::from_ptr,
        tile.mask.data()
    );
}

lsdtt_raster_cpu::lsdtt_raster_cpu(
    lsdtt_array2d_f&&        data,
    lsdtt_array2d_u8&&       mask,
    const lsdtt_raster_meta& meta
)
    : data_(std::move(data))
    , mask_(std::move(mask))
    , meta_(meta)
{
}

lsdtt_raster_cpu::lsdtt_raster_cpu(
    std::vector<float>&&     data,
    std::vector<uint8_t>&&   mask,
    const lsdtt_raster_meta& meta
)
    : data_(meta.n_rows, meta.n_cols, std::move(data))
    , mask_(meta.n_rows, meta.n_cols, std::move(mask))
    , meta_(meta)
{
}

// -----------------------------------------------------------------------------
// Loader-to-backend transfer
// -----------------------------------------------------------------------------

void lsdtt_raster_cpu::set_from_vectors(
    std::vector<float>&&        data,
    std::vector<std::uint8_t>&& mask,
    const lsdtt_raster_meta&    meta
)
{
    meta_ = meta;

    data_ = lsdtt_array2d_f(
        meta.n_rows,
        meta.n_cols,
        std::move(data)
    );

    mask_ = lsdtt_array2d_u8(
        meta.n_rows,
        meta.n_cols,
        std::move(mask)
    );
}

// -----------------------------------------------------------------------------
// Private helpers
// -----------------------------------------------------------------------------

void lsdtt_raster_cpu::validate_storage_sizes(const char* function_name) const
{
    if (data_.size() != mask_.size())
    {
        throw std::runtime_error(
            std::string("lsdtt_raster_cpu::") + function_name +
            "(): data and mask sizes do not match"
        );
    }
}

std::int64_t lsdtt_raster_cpu::checked_cell_count(const char* function_name) const
{
    validate_storage_sizes(function_name);

    if (data_.size() >
        static_cast<std::size_t>(std::numeric_limits<std::int64_t>::max()))
    {
        throw std::runtime_error(
            std::string("lsdtt_raster_cpu::") + function_name +
            "(): raster contains more cells than fit in int64_t"
        );
    }

    return static_cast<std::int64_t>(data_.size());
}

std::vector<float> lsdtt_raster_cpu::collect_valid_values() const
{
    const std::int64_t n_valid = count_valid_pixels();

    if (n_valid == 0)
    {
        throw std::runtime_error(
            "lsdtt_raster_cpu::collect_valid_values(): raster contains no valid pixels"
        );
    }

    const std::int64_t n_cells = checked_cell_count("collect_valid_values");

    const float*        data_ptr = data_.data();
    const std::uint8_t* mask_ptr = mask_.data();

    std::vector<float> values(static_cast<std::size_t>(n_valid));

#ifdef _OPENMP
    const int n_threads = omp_get_max_threads();
    std::vector<std::int64_t> thread_counts(static_cast<std::size_t>(n_threads), 0);

#pragma omp parallel
    {
        const int tid = omp_get_thread_num();

#pragma omp for schedule(static)
        for (std::int64_t i = 0; i < n_cells; ++i)
        {
            const std::size_t idx = static_cast<std::size_t>(i);
            if (mask_ptr[idx] != 0)
            {
                ++thread_counts[static_cast<std::size_t>(tid)];
            }
        }

#pragma omp barrier

        std::int64_t offset = 0;
        for (int t = 0; t < tid; ++t)
        {
            offset += thread_counts[static_cast<std::size_t>(t)];
        }

#pragma omp for schedule(static)
        for (std::int64_t i = 0; i < n_cells; ++i)
        {
            const std::size_t idx = static_cast<std::size_t>(i);
            if (mask_ptr[idx] != 0)
            {
                values[static_cast<std::size_t>(offset++)] = data_ptr[idx];
            }
        }
    }
#else
    std::int64_t write_index = 0;
    for (std::int64_t i = 0; i < n_cells; ++i)
    {
        const std::size_t idx = static_cast<std::size_t>(i);
        if (mask_ptr[idx] != 0)
        {
            values[static_cast<std::size_t>(write_index++)] = data_ptr[idx];
        }
    }
#endif

    return values;
}

// -----------------------------------------------------------------------------
// Whole-raster statistics over valid pixels (mask != 0)
// -----------------------------------------------------------------------------

std::int64_t lsdtt_raster_cpu::count_valid_pixels() const
{
    const std::int64_t n_cells = checked_cell_count("count_valid_pixels");
    const std::uint8_t* mask_ptr = mask_.data();

    std::int64_t n_valid = 0;

#ifdef _OPENMP
#pragma omp parallel for schedule(static) reduction(+:n_valid)
#endif
    for (std::int64_t i = 0; i < n_cells; ++i)
    {
        if (mask_ptr[static_cast<std::size_t>(i)] != 0)
        {
            ++n_valid;
        }
    }

    return n_valid;
}

float lsdtt_raster_cpu::min_value() const
{
    const std::int64_t n_cells = checked_cell_count("min_value");

    const float*        data_ptr = data_.data();
    const std::uint8_t* mask_ptr = mask_.data();

    std::int64_t n_valid = 0;
    double min_val = std::numeric_limits<double>::infinity();

#ifdef _OPENMP
#pragma omp parallel for schedule(static) reduction(+:n_valid) reduction(min:min_val)
#endif
    for (std::int64_t i = 0; i < n_cells; ++i)
    {
        const std::size_t idx = static_cast<std::size_t>(i);

        if (mask_ptr[idx] == 0) continue;

        const double v = static_cast<double>(data_ptr[idx]);
        if (v < min_val) min_val = v;
        ++n_valid;
    }

    if (n_valid == 0)
    {
        throw std::runtime_error(
            "lsdtt_raster_cpu::min_value(): raster contains no valid pixels"
        );
    }

    return static_cast<float>(min_val);
}

float lsdtt_raster_cpu::max_value() const
{
    const std::int64_t n_cells = checked_cell_count("max_value");

    const float*        data_ptr = data_.data();
    const std::uint8_t* mask_ptr = mask_.data();

    std::int64_t n_valid = 0;
    double max_val = -std::numeric_limits<double>::infinity();

#ifdef _OPENMP
#pragma omp parallel for schedule(static) reduction(+:n_valid) reduction(max:max_val)
#endif
    for (std::int64_t i = 0; i < n_cells; ++i)
    {
        const std::size_t idx = static_cast<std::size_t>(i);

        if (mask_ptr[idx] == 0) continue;

        const double v = static_cast<double>(data_ptr[idx]);
        if (v > max_val) max_val = v;
        ++n_valid;
    }

    if (n_valid == 0)
    {
        throw std::runtime_error(
            "lsdtt_raster_cpu::max_value(): raster contains no valid pixels"
        );
    }

    return static_cast<float>(max_val);
}

float lsdtt_raster_cpu::mean_value() const
{
    const std::int64_t n_cells = checked_cell_count("mean_value");

    const float*        data_ptr = data_.data();
    const std::uint8_t* mask_ptr = mask_.data();

    std::int64_t n_valid = 0;
    double sum = 0.0;

#ifdef _OPENMP
#pragma omp parallel for schedule(static) reduction(+:n_valid,sum)
#endif
    for (std::int64_t i = 0; i < n_cells; ++i)
    {
        const std::size_t idx = static_cast<std::size_t>(i);

        if (mask_ptr[idx] == 0) continue;

        sum += static_cast<double>(data_ptr[idx]);
        ++n_valid;
    }

    if (n_valid == 0)
    {
        throw std::runtime_error(
            "lsdtt_raster_cpu::mean_value(): raster contains no valid pixels"
        );
    }

    return static_cast<float>(sum / static_cast<double>(n_valid));
}

float lsdtt_raster_cpu::median_value() const
{
    return percentile(50.0);
}

float lsdtt_raster_cpu::percentile(double p) const
{
    if (p < 0.0 || p > 100.0)
    {
        throw std::invalid_argument(
            "lsdtt_raster_cpu::percentile(): percentile must be in [0, 100]"
        );
    }

    std::vector<float> values = collect_valid_values();

    const std::size_t n = values.size();
    if (n == 0)
    {
        throw std::runtime_error(
            "lsdtt_raster_cpu::percentile(): raster contains no valid pixels"
        );
    }

    if (n == 1)
    {
        return values[0];
    }

    const double rank =
        (p / 100.0) * static_cast<double>(n - 1);

    const std::size_t lower_index =
        static_cast<std::size_t>(std::floor(rank));
    const std::size_t upper_index =
        static_cast<std::size_t>(std::ceil(rank));

    std::nth_element(
        values.begin(),
        values.begin() + static_cast<std::ptrdiff_t>(upper_index),
        values.end()
    );
    const float upper_value = values[upper_index];

    if (lower_index == upper_index)
    {
        return upper_value;
    }

    std::nth_element(
        values.begin(),
        values.begin() + static_cast<std::ptrdiff_t>(lower_index),
        values.begin() + static_cast<std::ptrdiff_t>(upper_index) + 1
    );
    const float lower_value = values[lower_index];

    const double weight_upper = rank - static_cast<double>(lower_index);
    const double weight_lower = 1.0 - weight_upper;

    return static_cast<float>(
        weight_lower * static_cast<double>(lower_value) +
        weight_upper * static_cast<double>(upper_value)
    );
}

float lsdtt_raster_cpu::stddev_value() const
{
    const std::int64_t n_cells = checked_cell_count("stddev_value");

    const float*        data_ptr = data_.data();
    const std::uint8_t* mask_ptr = mask_.data();

    std::int64_t n_valid = 0;
    double sum = 0.0;

#ifdef _OPENMP
#pragma omp parallel for schedule(static) reduction(+:n_valid,sum)
#endif
    for (std::int64_t i = 0; i < n_cells; ++i)
    {
        const std::size_t idx = static_cast<std::size_t>(i);

        if (mask_ptr[idx] == 0) continue;

        sum += static_cast<double>(data_ptr[idx]);
        ++n_valid;
    }

    if (n_valid == 0)
    {
        throw std::runtime_error(
            "lsdtt_raster_cpu::stddev_value(): raster contains no valid pixels"
        );
    }

    const double mean = sum / static_cast<double>(n_valid);

    double squared_diff_sum = 0.0;

#ifdef _OPENMP
#pragma omp parallel for schedule(static) reduction(+:squared_diff_sum)
#endif
    for (std::int64_t i = 0; i < n_cells; ++i)
    {
        const std::size_t idx = static_cast<std::size_t>(i);

        if (mask_ptr[idx] == 0) continue;

        const double diff = static_cast<double>(data_ptr[idx]) - mean;
        squared_diff_sum += diff * diff;
    }

    const double variance =
        squared_diff_sum / static_cast<double>(n_valid);

    return static_cast<float>(std::sqrt(variance));
}

float lsdtt_raster_cpu::range_value() const
{
    const std::int64_t n_cells = checked_cell_count("range_value");

    const float*        data_ptr = data_.data();
    const std::uint8_t* mask_ptr = mask_.data();

    std::int64_t n_valid = 0;
    double min_val = std::numeric_limits<double>::infinity();
    double max_val = -std::numeric_limits<double>::infinity();

#ifdef _OPENMP
#pragma omp parallel for schedule(static) reduction(+:n_valid) reduction(min:min_val) reduction(max:max_val)
#endif
    for (std::int64_t i = 0; i < n_cells; ++i)
    {
        const std::size_t idx = static_cast<std::size_t>(i);

        if (mask_ptr[idx] == 0) continue;

        const double v = static_cast<double>(data_ptr[idx]);
        if (v < min_val) min_val = v;
        if (v > max_val) max_val = v;
        ++n_valid;
    }

    if (n_valid == 0)
    {
        throw std::runtime_error(
            "lsdtt_raster_cpu::range_value(): raster contains no valid pixels"
        );
    }

    return static_cast<float>(max_val - min_val);
}

std::int64_t lsdtt_raster_cpu::assert_mask()
{
    const std::int64_t n_cells = checked_cell_count("assert_mask");

    const float*        data_ptr = data_.data();
    std::uint8_t*       mask_ptr = mask_.data();
    const float         nodata   = static_cast<float>(meta_.no_data_value);
    const bool          nodata_is_nan = std::isnan(static_cast<double>(nodata));

    std::int64_t n_changed = 0;

#ifdef _OPENMP
#pragma omp parallel for schedule(static) reduction(+:n_changed)
#endif
    for (std::int64_t i = 0; i < n_cells; ++i)
    {
        const std::size_t idx = static_cast<std::size_t>(i);
        const float z = data_ptr[idx];

        const bool is_nodata =
            nodata_is_nan
                ? std::isnan(static_cast<double>(z))
                : (z == nodata || std::isnan(static_cast<double>(z)));

        if (is_nodata && mask_ptr[idx] != 0)
        {
            mask_ptr[idx] = 0;
            ++n_changed;
        }
    }

    return n_changed;
}

// -----------------------------------------------------------------------------
// Whole-raster arithmetic over valid pixels (mask != 0)
// -----------------------------------------------------------------------------

void lsdtt_raster_cpu::add_value(float value)
{
    const std::int64_t n_cells = checked_cell_count("add_value");

    float*              data_ptr = data_.data();
    const std::uint8_t* mask_ptr = mask_.data();

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (std::int64_t i = 0; i < n_cells; ++i)
    {
        const std::size_t idx = static_cast<std::size_t>(i);

        if (mask_ptr[idx] != 0)
        {
            data_ptr[idx] += value;
        }
    }
}

void lsdtt_raster_cpu::multiply_value(float factor)
{
    const std::int64_t n_cells = checked_cell_count("multiply_value");

    float*              data_ptr = data_.data();
    const std::uint8_t* mask_ptr = mask_.data();

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (std::int64_t i = 0; i < n_cells; ++i)
    {
        const std::size_t idx = static_cast<std::size_t>(i);

        if (mask_ptr[idx] != 0)
        {
            data_ptr[idx] *= factor;
        }
    }
}

void lsdtt_raster_cpu::clip(float lower, float upper)
{
    if (lower > upper)
    {
        throw std::invalid_argument(
            "lsdtt_raster_cpu::clip(): lower must be <= upper"
        );
    }

    const std::int64_t n_cells = checked_cell_count("clip");

    float*              data_ptr = data_.data();
    const std::uint8_t* mask_ptr = mask_.data();

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (std::int64_t i = 0; i < n_cells; ++i)
    {
        const std::size_t idx = static_cast<std::size_t>(i);

        if (mask_ptr[idx] == 0) continue;

        if (data_ptr[idx] < lower)      data_ptr[idx] = lower;
        else if (data_ptr[idx] > upper) data_ptr[idx] = upper;
    }
}

} // namespace lsdtt
