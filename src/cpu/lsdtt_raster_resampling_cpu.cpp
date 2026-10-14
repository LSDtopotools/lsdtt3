// lsdtt_raster_resampling_cpu.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Boris Gailleton and Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html

#include <cctype>
#include "lsdtt_raster_resampling_cpu.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <map>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace lsdtt {
namespace {

struct sample
{
    double value = 0.0;
    bool valid = false;
};

bool close(double a, double b)
{
    const double scale = std::max({1.0, std::abs(a), std::abs(b)});
    return std::abs(a - b) <= 1e-10 * scale;
}

bool same_extent(const lsdtt_raster_meta& a, const lsdtt_raster_meta& b)
{
    return close(a.x_min, b.x_min) && close(a.y_min, b.y_min) &&
           close(a.x_max(), b.x_max()) && close(a.y_max(), b.y_max());
}

void validate(const lsdtt_raster_cpu& raster, const char* name)
{
    const auto& meta = raster.get_meta();
    if (!meta.is_valid() ||
        raster.get_data().size() != static_cast<std::size_t>(meta.n_pixels()) ||
        raster.get_mask().size() != static_cast<std::size_t>(meta.n_pixels()))
    {
        throw std::invalid_argument(std::string(name) + " is not a valid raster");
    }
}

void validate_pair(const lsdtt_raster_cpu& source, const lsdtt_raster_cpu& destination)
{
    validate(source, "source");
    validate(destination, "destination");
    const auto& a = source.get_meta();
    const auto& b = destination.get_meta();
    const bool one_crs_missing = a.crs_wkt.empty() != b.crs_wkt.empty();
    const bool known_epsg_match = a.epsg_code > 0 && b.epsg_code > 0 &&
                                  a.epsg_code == b.epsg_code;
    if (one_crs_missing ||
        (!a.crs_wkt.empty() && !known_epsg_match && a.crs_wkt != b.crs_wkt))
    {
        throw std::invalid_argument("source and destination CRS must match");
    }
}

bool usable(const lsdtt_raster_cpu& raster, int row, int col)
{
    if (row < 0 || row >= raster.n_rows() || col < 0 || col >= raster.n_cols())
        return false;
    const std::size_t flat = static_cast<std::size_t>(row) * raster.n_cols() + col;
    return raster.get_mask()[flat] != 0 && std::isfinite(raster.get_data()[flat]);
}

sample nearest_sample(const lsdtt_raster_cpu& source, double x, double y)
{
    const auto& meta = source.get_meta();
    if (x < meta.x_min || x > meta.x_max() || y < meta.y_min || y > meta.y_max())
        return {};

    int col = static_cast<int>(std::floor((x - meta.x_min) / meta.cell_size));
    int row = static_cast<int>(std::floor((meta.y_max() - y) / meta.cell_size));
    col = std::clamp(col, 0, meta.n_cols - 1);
    row = std::clamp(row, 0, meta.n_rows - 1);
    if (!usable(source, row, col)) return {};
    return {source.get_data_array()(row, col), true};
}

sample linear_sample(const lsdtt_raster_cpu& source, double x, double y)
{
    const auto& meta = source.get_meta();
    if (x < meta.x_min || x > meta.x_max() || y < meta.y_min || y > meta.y_max())
        return {};

    const double col_f = std::clamp(
        (x - meta.x_min) / meta.cell_size - 0.5, 0.0,
        static_cast<double>(meta.n_cols - 1));
    const double row_f = std::clamp(
        (meta.y_max() - y) / meta.cell_size - 0.5, 0.0,
        static_cast<double>(meta.n_rows - 1));
    const int col0 = static_cast<int>(std::floor(col_f));
    const int row0 = static_cast<int>(std::floor(row_f));
    const int col1 = std::min(col0 + 1, meta.n_cols - 1);
    const int row1 = std::min(row0 + 1, meta.n_rows - 1);
    const double dc = col_f - col0;
    const double dr = row_f - row0;
    const int rows[4] = {row0, row0, row1, row1};
    const int cols[4] = {col0, col1, col0, col1};
    const double weights[4] = {
        (1.0 - dr) * (1.0 - dc), (1.0 - dr) * dc,
        dr * (1.0 - dc), dr * dc};

    double total = 0.0;
    double weight = 0.0;
    for (int i = 0; i < 4; ++i)
    {
        if (weights[i] > 0.0 && usable(source, rows[i], cols[i]))
        {
            total += weights[i] * source.get_data_array()(rows[i], cols[i]);
            weight += weights[i];
        }
    }
    return weight > 0.0 ? sample{total / weight, true} : sample{};
}

double cubic_weight(double x)
{
    x = std::abs(x);
    if (x <= 1.0) return ((1.5 * x - 2.5) * x) * x + 1.0;
    if (x < 2.0) return ((-0.5 * x + 2.5) * x - 4.0) * x + 2.0;
    return 0.0;
}

sample cubic_sample(const lsdtt_raster_cpu& source, double x, double y)
{
    const auto& meta = source.get_meta();
    if (x < meta.x_min || x > meta.x_max() || y < meta.y_min || y > meta.y_max())
        return {};

    const double col_f = std::clamp(
        (x - meta.x_min) / meta.cell_size - 0.5, 0.0,
        static_cast<double>(meta.n_cols - 1));
    const double row_f = std::clamp(
        (meta.y_max() - y) / meta.cell_size - 0.5, 0.0,
        static_cast<double>(meta.n_rows - 1));
    const int col0 = static_cast<int>(std::floor(col_f));
    const int row0 = static_cast<int>(std::floor(row_f));

    double total = 0.0;
    double weight = 0.0;
    for (int j = -1; j <= 2; ++j)
    {
        const int row = std::clamp(row0 + j, 0, meta.n_rows - 1);
        const double wy = cubic_weight(row_f - (row0 + j));
        for (int i = -1; i <= 2; ++i)
        {
            const int col = std::clamp(col0 + i, 0, meta.n_cols - 1);
            const double w = wy * cubic_weight(col_f - (col0 + i));
            if (std::abs(w) <= 1e-15) continue;
            // Cubic interpolation across nodata is poorly defined. Fall back
            // to the local linear stencil rather than ringing across a hole.
            if (!usable(source, row, col)) return linear_sample(source, x, y);
            total += w * source.get_data_array()(row, col);
            weight += w;
        }
    }
    return std::abs(weight) > 1e-12 ? sample{total / weight, true} : sample{};
}

sample point_sample(
    const lsdtt_raster_cpu& source,
    double x,
    double y,
    raster_resampling_method method)
{
    if (method == raster_resampling_method::nearest)
        return nearest_sample(source, x, y);
    if (method == raster_resampling_method::linear)
        return linear_sample(source, x, y);
    return cubic_sample(source, x, y);
}

bool is_point_method(raster_resampling_method method)
{
    return method == raster_resampling_method::nearest ||
           method == raster_resampling_method::linear ||
           method == raster_resampling_method::cubic;
}

sample aggregate_sample(
    const lsdtt_raster_cpu& source,
    double x0,
    double x1,
    double y0,
    double y1,
    raster_resampling_method method)
{
    const auto& meta = source.get_meta();
    const double sx0 = std::max(x0, meta.x_min);
    const double sx1 = std::min(x1, meta.x_max());
    const double sy0 = std::max(y0, meta.y_min);
    const double sy1 = std::min(y1, meta.y_max());
    if (sx1 <= sx0 || sy1 <= sy0) return {};

    int col0 = static_cast<int>(std::floor((sx0 - meta.x_min) / meta.cell_size));
    int col1 = static_cast<int>(std::ceil((sx1 - meta.x_min) / meta.cell_size)) - 1;
    int row0 = static_cast<int>(std::floor((meta.y_max() - sy1) / meta.cell_size));
    int row1 = static_cast<int>(std::ceil((meta.y_max() - sy0) / meta.cell_size)) - 1;
    col0 = std::clamp(col0, 0, meta.n_cols - 1);
    col1 = std::clamp(col1, 0, meta.n_cols - 1);
    row0 = std::clamp(row0, 0, meta.n_rows - 1);
    row1 = std::clamp(row1, 0, meta.n_rows - 1);

    double total = 0.0;
    double area_total = 0.0;
    std::size_t count = 0;
    double extremum = method == raster_resampling_method::minimum
        ? std::numeric_limits<double>::infinity()
        : -std::numeric_limits<double>::infinity();
    std::map<float, double> frequencies;

    for (int row = row0; row <= row1; ++row)
    {
        const double cell_y1 = meta.y_max() - row * meta.cell_size;
        const double cell_y0 = cell_y1 - meta.cell_size;
        const double overlap_y = std::max(0.0, std::min(sy1, cell_y1) - std::max(sy0, cell_y0));
        for (int col = col0; col <= col1; ++col)
        {
            if (!usable(source, row, col)) continue;
            const double cell_x0 = meta.x_min + col * meta.cell_size;
            const double cell_x1 = cell_x0 + meta.cell_size;
            const double overlap_x = std::max(0.0, std::min(sx1, cell_x1) - std::max(sx0, cell_x0));
            const double area = overlap_x * overlap_y;
            if (area <= 0.0) continue;

            const float value = source.get_data_array()(row, col);
            ++count;
            if (method == raster_resampling_method::mean) total += value;
            else if (method == raster_resampling_method::sum)
                total += value * area / meta.pixel_area();
            else if (method == raster_resampling_method::conservative_mean)
                total += value * area;
            else if (method == raster_resampling_method::minimum)
                extremum = std::min(extremum, static_cast<double>(value));
            else if (method == raster_resampling_method::maximum)
                extremum = std::max(extremum, static_cast<double>(value));
            else if (method == raster_resampling_method::mode)
                frequencies[value] += area;
            area_total += area;
        }
    }

    if (count == 0) return {};
    if (method == raster_resampling_method::mean) return {total / count, true};
    if (method == raster_resampling_method::sum) return {total, true};
    if (method == raster_resampling_method::conservative_mean)
        return {total / area_total, true};
    if (method == raster_resampling_method::minimum ||
        method == raster_resampling_method::maximum)
        return {extremum, true};

    auto best = frequencies.begin();
    for (auto it = frequencies.begin(); it != frequencies.end(); ++it)
        if (it->second > best->second) best = it;
    return {best->first, true};
}

void write_sample(lsdtt_raster_cpu& destination, int row, int col, sample value)
{
    const std::size_t flat = static_cast<std::size_t>(row) * destination.n_cols() + col;
    if (value.valid && std::isfinite(value.value))
    {
        destination.get_data()[flat] = static_cast<float>(value.value);
        destination.get_mask_array().get_vector()[flat] = 1;
    }
    else
    {
        destination.get_data()[flat] = static_cast<float>(destination.get_meta().no_data_value);
        destination.get_mask_array().get_vector()[flat] = 0;
    }
}

void generic_resample(
    const lsdtt_raster_cpu& source,
    lsdtt_raster_cpu& destination,
    raster_resampling_method method)
{
    const auto& target = destination.get_meta();

#ifdef _OPENMP
#pragma omp parallel for collapse(2) schedule(static)
#endif
    for (int row = 0; row < target.n_rows; ++row)
    {
        for (int col = 0; col < target.n_cols; ++col)
        {
            const double x0 = target.x_min + col * target.cell_size;
            const double x1 = x0 + target.cell_size;
            const double y1 = target.y_max() - row * target.cell_size;
            const double y0 = y1 - target.cell_size;
            const sample value = is_point_method(method)
                ? point_sample(source, 0.5 * (x0 + x1), 0.5 * (y0 + y1), method)
                : aggregate_sample(source, x0, x1, y0, y1, method);
            write_sample(destination, row, col, value);
        }
    }
}

void double_impl(
    const lsdtt_raster_cpu& source,
    lsdtt_raster_cpu& destination,
    raster_resampling_method method)
{
    const auto& target = destination.get_meta();
    const bool direct = !is_point_method(method) ||
                        method == raster_resampling_method::nearest;

#ifdef _OPENMP
#pragma omp parallel for collapse(2) schedule(static)
#endif
    for (int row = 0; row < target.n_rows; ++row)
    {
        for (int col = 0; col < target.n_cols; ++col)
        {
            sample value;
            if (direct)
            {
                const int source_row = row / 2;
                const int source_col = col / 2;
                if (usable(source, source_row, source_col))
                {
                    value = {source.get_data_array()(source_row, source_col), true};
                    if (method == raster_resampling_method::sum) value.value *= 0.25;
                }
            }
            else
            {
                const double x = target.x_min + (col + 0.5) * target.cell_size;
                const double y = target.y_max() - (row + 0.5) * target.cell_size;
                value = point_sample(source, x, y, method);
            }
            write_sample(destination, row, col, value);
        }
    }
}

void halve_impl(
    const lsdtt_raster_cpu& source,
    lsdtt_raster_cpu& destination,
    raster_resampling_method method)
{
    if (method == raster_resampling_method::cubic)
    {
        generic_resample(source, destination, method);
        return;
    }

    const auto& target = destination.get_meta();
#ifdef _OPENMP
#pragma omp parallel for collapse(2) schedule(static)
#endif
    for (int row = 0; row < target.n_rows; ++row)
    {
        for (int col = 0; col < target.n_cols; ++col)
        {
            const int row0 = 2 * row;
            const int col0 = 2 * col;
            if (method == raster_resampling_method::nearest)
            {
                sample value;
                if (usable(source, row0 + 1, col0 + 1))
                    value = {source.get_data_array()(row0 + 1, col0 + 1), true};
                write_sample(destination, row, col, value);
                continue;
            }

            double total = 0.0;
            double extremum = method == raster_resampling_method::minimum
                ? std::numeric_limits<double>::infinity()
                : -std::numeric_limits<double>::infinity();
            int count = 0;
            std::map<float, int> frequencies;
            for (int dr = 0; dr < 2; ++dr)
            {
                for (int dc = 0; dc < 2; ++dc)
                {
                    if (!usable(source, row0 + dr, col0 + dc)) continue;
                    const float v = source.get_data_array()(row0 + dr, col0 + dc);
                    ++count;
                    total += v;
                    extremum = method == raster_resampling_method::minimum
                        ? std::min(extremum, static_cast<double>(v))
                        : std::max(extremum, static_cast<double>(v));
                    ++frequencies[v];
                }
            }
            sample value;
            if (count > 0)
            {
                if (method == raster_resampling_method::sum) value = {total, true};
                else if (method == raster_resampling_method::minimum ||
                         method == raster_resampling_method::maximum)
                    value = {extremum, true};
                else if (method == raster_resampling_method::mode)
                {
                    auto best = frequencies.begin();
                    for (auto it = frequencies.begin(); it != frequencies.end(); ++it)
                        if (it->second > best->second) best = it;
                    value = {best->first, true};
                }
                else value = {total / count, true};
            }
            write_sample(destination, row, col, value);
        }
    }
}

int lower_power_of_two(int value)
{
    int power = 1;
    while (power <= value / 2) power *= 2;
    return power;
}

int higher_power_of_two(int value)
{
    const int lower = lower_power_of_two(value);
    return lower == value ? value : lower * 2;
}

int selected_power_of_two(int value, power_of_two_direction direction)
{
    const int lower = lower_power_of_two(value);
    const int higher = higher_power_of_two(value);
    if (direction == power_of_two_direction::lower) return lower;
    if (direction == power_of_two_direction::higher) return higher;
    return value - lower < higher - value ? lower : higher;
}

lsdtt_raster_cpu empty_raster(const lsdtt_raster_meta& meta)
{
    const std::size_t size = static_cast<std::size_t>(meta.n_pixels());
    return lsdtt_raster_cpu(
        std::vector<float>(size, static_cast<float>(meta.no_data_value)),
        std::vector<std::uint8_t>(size, 0), meta);
}

} // namespace

void resample_raster(
    const lsdtt_raster_cpu& source,
    lsdtt_raster_cpu& destination,
    raster_resampling_method method)
{
    validate_pair(source, destination);
    if (&source == &destination) return;

    const auto& a = source.get_meta();
    const auto& b = destination.get_meta();
    if (same_extent(a, b) && b.n_rows == 2 * a.n_rows && b.n_cols == 2 * a.n_cols &&
        close(2.0 * b.cell_size, a.cell_size))
    {
        double_impl(source, destination, method);
        return;
    }
    if (same_extent(a, b) && a.n_rows == 2 * b.n_rows && a.n_cols == 2 * b.n_cols &&
        close(2.0 * a.cell_size, b.cell_size))
    {
        halve_impl(source, destination, method);
        return;
    }
    generic_resample(source, destination, method);
}

void double_raster(
    const lsdtt_raster_cpu& source,
    lsdtt_raster_cpu& destination,
    raster_resampling_method method)
{
    validate_pair(source, destination);
    const auto& a = source.get_meta();
    const auto& b = destination.get_meta();
    if (!same_extent(a, b) || b.n_rows != 2 * a.n_rows || b.n_cols != 2 * a.n_cols ||
        !close(2.0 * b.cell_size, a.cell_size))
        throw std::invalid_argument("double_raster requires a same-extent 2x destination");
    double_impl(source, destination, method);
}

void halve_raster(
    const lsdtt_raster_cpu& source,
    lsdtt_raster_cpu& destination,
    raster_resampling_method method)
{
    validate_pair(source, destination);
    const auto& a = source.get_meta();
    const auto& b = destination.get_meta();
    if (!same_extent(a, b) || a.n_rows != 2 * b.n_rows || a.n_cols != 2 * b.n_cols ||
        !close(2.0 * a.cell_size, b.cell_size))
        throw std::invalid_argument("halve_raster requires a same-extent 1/2x destination");
    halve_impl(source, destination, method);
}

lsdtt_raster_cpu pad_to_nearest_power_of_two(const lsdtt_raster_cpu& source)
{
    validate(source, "source");
    const auto& source_meta = source.get_meta();
    auto target_meta = source_meta;
    target_meta.n_rows = higher_power_of_two(source_meta.n_rows);
    target_meta.n_cols = higher_power_of_two(source_meta.n_cols);
    lsdtt_raster_cpu destination = empty_raster(target_meta);

    const int row_offset = target_meta.n_rows - source_meta.n_rows;
    for (int row = 0; row < source_meta.n_rows; ++row)
    {
        for (int col = 0; col < source_meta.n_cols; ++col)
        {
            const std::size_t source_flat = static_cast<std::size_t>(row) * source_meta.n_cols + col;
            const std::size_t target_flat = static_cast<std::size_t>(row + row_offset) * target_meta.n_cols + col;
            destination.get_data()[target_flat] = source.get_data()[source_flat];
            destination.get_mask_array().get_vector()[target_flat] = source.get_mask()[source_flat];
        }
    }
    return destination;
}

lsdtt_raster_cpu interpolate_to_nearest_power_of_two(
    const lsdtt_raster_cpu& source,
    power_of_two_direction direction,
    raster_resampling_method method)
{
    validate(source, "source");
    const auto& source_meta = source.get_meta();
    auto target_meta = source_meta;
    target_meta.n_rows = selected_power_of_two(source_meta.n_rows, direction);
    target_meta.n_cols = selected_power_of_two(source_meta.n_cols, direction);

    const double row_scale = static_cast<double>(source_meta.n_rows) / target_meta.n_rows;
    const double col_scale = static_cast<double>(source_meta.n_cols) / target_meta.n_cols;
    if (!close(row_scale, col_scale))
    {
        throw std::invalid_argument(
            "power-of-two interpolation cannot preserve this raster extent with square cells; "
            "supply explicit compatible destination metadata to resample_raster");
    }
    target_meta.cell_size = source_meta.cell_size * col_scale;
    lsdtt_raster_cpu destination = empty_raster(target_meta);
    resample_raster(source, destination, method);
    return destination;
}


namespace
{
// Wrap an index onto [0, n) when the axis is periodic; otherwise report "off
// the edge" as -1.
inline int wrap_index(const int idx, const int n, const bool periodic)
{
    if (periodic) return (idx % n + n) % n;
    return (idx >= 0 && idx < n) ? idx : -1;
}

inline bool axis_is_periodic(const std::string& bc, const std::size_t a,
                             const std::size_t b)
{
    return bc.size() >= 4 &&
           std::tolower(static_cast<unsigned char>(bc[a])) == 'p' &&
           std::tolower(static_cast<unsigned char>(bc[b])) == 'p';
}
} // namespace

lsdtt_raster_cpu pad_raster_by_boundary_conditions(
    const lsdtt_raster_cpu& raster,
    const std::string&      boundary_conditions,
    const int               border_pixels
)
{
    if (border_pixels < 0)
    {
        throw std::invalid_argument(
            "pad_raster_by_boundary_conditions: border_pixels must be >= 0");
    }

    const lsdtt_raster_meta& meta = raster.get_meta();
    const int nrows = meta.n_rows;
    const int ncols = meta.n_cols;

    // N and S are indices 0 and 2; E and W are 1 and 3.
    const bool periodic_ns = axis_is_periodic(boundary_conditions, 0, 2);
    const bool periodic_ew = axis_is_periodic(boundary_conditions, 1, 3);

    const int p_rows = nrows + 2 * border_pixels;
    const int p_cols = ncols + 2 * border_pixels;

    std::vector<float> p_data(
        static_cast<std::size_t>(p_rows) * static_cast<std::size_t>(p_cols),
        static_cast<float>(meta.no_data_value));
    std::vector<std::uint8_t> p_mask(
        static_cast<std::size_t>(p_rows) * static_cast<std::size_t>(p_cols), 0u);

    const std::vector<float>&        d = raster.get_data();
    const std::vector<std::uint8_t>& m = raster.get_mask();

    for (int r = -border_pixels; r < nrows + border_pixels; ++r)
    {
        const int wr = wrap_index(r, nrows, periodic_ns);
        if (wr < 0) continue;

        for (int c = -border_pixels; c < ncols + border_pixels; ++c)
        {
            const int wc = wrap_index(c, ncols, periodic_ew);
            if (wc < 0) continue;

            const std::size_t src =
                static_cast<std::size_t>(wr) * static_cast<std::size_t>(ncols) +
                static_cast<std::size_t>(wc);
            const std::size_t dst =
                static_cast<std::size_t>(r + border_pixels) *
                static_cast<std::size_t>(p_cols) +
                static_cast<std::size_t>(c + border_pixels);

            p_data[dst] = d[src];
            p_mask[dst] = m[src];
        }
    }

    lsdtt_raster_meta p_meta = meta;
    p_meta.n_rows = p_rows;
    p_meta.n_cols = p_cols;

    return lsdtt_raster_cpu(std::move(p_data), std::move(p_mask), p_meta);
}


lsdtt_raster_cpu crop_raster_border(
    const lsdtt_raster_cpu&  padded,
    const lsdtt_raster_meta& original_meta,
    const int                border_pixels
)
{
    const int nrows = original_meta.n_rows;
    const int ncols = original_meta.n_cols;
    const int p_cols = padded.get_meta().n_cols;

    if (nrows <= 0 || ncols <= 0)
    {
        throw std::invalid_argument(
            "crop_raster_border: cropping would leave an empty raster");
    }

    std::vector<float> c_data(
        static_cast<std::size_t>(nrows) * static_cast<std::size_t>(ncols),
        static_cast<float>(original_meta.no_data_value));
    std::vector<std::uint8_t> c_mask(
        static_cast<std::size_t>(nrows) * static_cast<std::size_t>(ncols), 0u);

    const std::vector<float>&        pd = padded.get_data();
    const std::vector<std::uint8_t>& pm = padded.get_mask();

    for (int r = 0; r < nrows; ++r)
    {
        for (int c = 0; c < ncols; ++c)
        {
            const std::size_t src =
                static_cast<std::size_t>(r + border_pixels) *
                static_cast<std::size_t>(p_cols) +
                static_cast<std::size_t>(c + border_pixels);
            const std::size_t dst =
                static_cast<std::size_t>(r) * static_cast<std::size_t>(ncols) +
                static_cast<std::size_t>(c);

            c_data[dst] = pd[src];
            c_mask[dst] = pm[src];
        }
    }

    return lsdtt_raster_cpu(std::move(c_data), std::move(c_mask), original_meta);
}

} // namespace lsdtt
