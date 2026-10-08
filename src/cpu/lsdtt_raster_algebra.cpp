// lsdtt_raster_algebra.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html

#include "lsdtt_raster_algebra.hpp"
#include "../utils/lsdtt_statistics.hpp"

#include "../core/lsdtt_data_loader.hpp"
#include "../core/lsdtt_stat_csv.hpp"
#include "../utils/lsdtt_distributions.hpp"
#include "../utils/lsdtt_histogram.hpp"
#include "../utils/lsdtt_natural_breaks.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#include <iostream>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace lsdtt {
namespace algebra {

lsdtt_raster_cpu raster_transform(
    const lsdtt_raster_cpu&             source,
    const std::function<float(float)>& transform)
{
    const auto& data = source.get_data();
    const auto& mask = source.get_mask();
    const float nodata = static_cast<float>(source.get_meta().no_data_value);

    std::vector<float> out_data(data.size(), nodata);
    std::vector<uint8_t> out_mask(mask.size(), 0u);
    for (std::size_t i = 0; i < data.size(); ++i)
    {
        if (!mask[i]) continue;
        out_data[i] = transform(data[i]);
        out_mask[i] = 1u;
    }

    return lsdtt_raster_cpu(
        std::move(out_data), std::move(out_mask), source.get_meta());
}

raster_value_set_mask_result raster_mask_from_rounded_value_set(
    const lsdtt_raster_cpu& source,
    const std::vector<long>& values)
{
    const auto& src_data = source.get_data();
    const auto& src_mask = source.get_mask();
    if (src_data.size() != src_mask.size())
    {
        throw std::runtime_error(
            "raster_mask_from_rounded_value_set: data/mask size mismatch");
    }

    raster_value_set_mask_result result;
    std::vector<float> out_data(src_data.size(), 0.0f);
    std::vector<std::uint8_t> out_mask = src_mask;
    for (std::size_t i = 0; i < src_data.size(); ++i)
    {
        if (src_mask[i] == 0 || !std::isfinite(src_data[i])) continue;
        const long rounded = std::lround(static_cast<double>(src_data[i]));
        if (std::find(values.begin(), values.end(), rounded) != values.end())
        {
            out_data[i] = 1.0f;
            ++result.selected_pixels;
        }
    }
    result.raster = lsdtt_raster_cpu(
        std::move(out_data), std::move(out_mask), source.get_meta());
    return result;
}

raster_active_mask_result raster_active_mask_and_not_excluded(
    const lsdtt_raster_cpu& domain,
    const lsdtt_raster_cpu* active_mask,
    const lsdtt_raster_cpu& exclusion_mask)
{
    const auto& meta = domain.get_meta();
    const auto dimensions_match = [&meta](const lsdtt_raster_cpu& raster)
    {
        return raster.get_meta().n_rows == meta.n_rows &&
               raster.get_meta().n_cols == meta.n_cols;
    };
    if (!dimensions_match(exclusion_mask) ||
        (active_mask != nullptr && !dimensions_match(*active_mask)))
    {
        throw std::invalid_argument(
            "raster_active_mask_and_not_excluded: raster dimensions do not match domain");
    }

    const auto& ice_data = exclusion_mask.get_data();
    const auto& ice_mask = exclusion_mask.get_mask();
    const float ice_nodata =
        static_cast<float>(exclusion_mask.get_meta().no_data_value);
    const std::size_t n = static_cast<std::size_t>(meta.n_pixels());
    if (ice_data.size() != n || ice_mask.size() != n ||
        (active_mask != nullptr &&
         (active_mask->get_data().size() != n || active_mask->get_mask().size() != n)))
    {
        throw std::runtime_error(
            "raster_active_mask_and_not_excluded: data/mask size mismatch");
    }

    raster_active_mask_result result;
    std::vector<float> out_data(n, 0.0f);
    std::vector<std::uint8_t> out_mask(n, 1);
    for (std::size_t i = 0; i < n; ++i)
    {
        const bool in_active = active_mask == nullptr ||
            (active_mask->get_mask()[i] != 0 && active_mask->get_data()[i] != 0.0f);
        const bool excluded = ice_mask[i] != 0 && ice_data[i] != ice_nodata &&
                              ice_data[i] != 0.0f;
        if (excluded && in_active) ++result.excluded_pixels;
        if (in_active && !excluded)
        {
            out_data[i] = 1.0f;
            ++result.active_pixels;
        }
    }
    result.raster = lsdtt_raster_cpu(
        std::move(out_data), std::move(out_mask), meta);
    return result;
}

raster_jenks_classification_result classify_jenks_categories(
    const lsdtt_raster_cpu& input,
    const int               requested_classes,
    const int               max_histogram_bins)
{
    const auto& src_data = input.get_data();
    const auto& src_mask = input.get_mask();
    const auto& meta = input.get_meta();

    if (src_data.size() != src_mask.size())
    {
        throw std::runtime_error("classify_jenks_categories: data/mask size mismatch");
    }

    std::vector<double> values;
    values.reserve(src_data.size());
    for (std::size_t i = 0; i < src_data.size(); ++i)
    {
        if (src_mask[i] == 0) continue;
        const float v = src_data[i];
        if (!std::isfinite(v) || meta.is_nodata(v)) continue;
        values.push_back(static_cast<double>(v));
    }

    lsdtt_jenks_options opts;
    opts.n_classes = requested_classes;
    opts.max_histogram_bins = max_histogram_bins;

    raster_jenks_classification_result result;
    result.jenks = compute_jenks_breaks(values, opts);

    std::vector<float> out_data = src_data;
    std::vector<std::uint8_t> out_mask = src_mask;
    const float nodata_f = static_cast<float>(meta.no_data_value);

    for (std::size_t i = 0; i < src_data.size(); ++i)
    {
        if (src_mask[i] == 0)
        {
            out_data[i] = nodata_f;
            continue;
        }
        const float v = src_data[i];
        if (!std::isfinite(v) || meta.is_nodata(v))
        {
            out_data[i] = nodata_f;
            out_mask[i] = 0;
            continue;
        }

        const int cls = classify_jenks_value(
            static_cast<double>(v), result.jenks.upper_breaks);
        if (cls <= 0)
        {
            out_data[i] = nodata_f;
            out_mask[i] = 0;
            continue;
        }
        out_data[i] = static_cast<float>(cls);
        ++result.classified_pixels;
    }

    result.classified_raster = lsdtt_raster_cpu(
        std::move(out_data), std::move(out_mask), meta);
    return result;
}

// =============================================================================
// assert_compatible
// =============================================================================

void assert_compatible(
    const lsdtt_raster_cpu& a,
    const lsdtt_raster_cpu& b)
{
    const lsdtt_raster_meta& ma = a.get_meta();
    const lsdtt_raster_meta& mb = b.get_meta();

    std::ostringstream err;
    bool ok = true;

    if (ma.n_rows != mb.n_rows || ma.n_cols != mb.n_cols)
    {
        err << "  dimensions differ: ("
            << ma.n_rows << "x" << ma.n_cols << ") vs ("
            << mb.n_rows << "x" << mb.n_cols << ")\n";
        ok = false;
    }
    if (ma.cell_size != mb.cell_size)
    {
        err << "  cell_size differs: " << ma.cell_size
            << " vs " << mb.cell_size << "\n";
        ok = false;
    }
    if (ma.x_min != mb.x_min || ma.y_min != mb.y_min)
    {
        err << "  origin differs: ("
            << ma.x_min << ", " << ma.y_min << ") vs ("
            << mb.x_min << ", " << mb.y_min << ")\n";
        ok = false;
    }
    if (ma.no_data_value != mb.no_data_value)
    {
        err << "  no_data_value differs: "
            << ma.no_data_value << " vs " << mb.no_data_value << "\n";
        ok = false;
    }

    if (!ok)
    {
        throw std::invalid_argument(
            "lsdtt raster algebra: incompatible rasters -\n" + err.str()
        );
    }
}

// =============================================================================
// Arithmetic operations
// =============================================================================

lsdtt_raster_cpu raster_add(
    const lsdtt_raster_cpu& a,
    const lsdtt_raster_cpu& b)
{
    assert_compatible(a, b);

    const std::vector<float>&   da = a.get_data();
    const std::vector<uint8_t>& ma = a.get_mask();
    const std::vector<float>&   db = b.get_data();
    const std::vector<uint8_t>& mb = b.get_mask();

    const int   n      = static_cast<int>(da.size());
    const float nodata = static_cast<float>(a.get_meta().no_data_value);

    std::vector<float>   out_data(n);
    std::vector<uint8_t> out_mask(n);

    #pragma omp parallel for schedule(static)
    for (int i = 0; i < n; ++i)
    {
        if (ma[i] && mb[i])
        {
            out_data[i] = da[i] + db[i];
            out_mask[i] = 1u;
        }
        else
        {
            out_data[i] = nodata;
            out_mask[i] = 0u;
        }
    }

    return lsdtt_raster_cpu(
        std::move(out_data), std::move(out_mask), a.get_meta()
    );
}

lsdtt_raster_cpu raster_subtract(
    const lsdtt_raster_cpu& a,
    const lsdtt_raster_cpu& b)
{
    assert_compatible(a, b);

    const std::vector<float>&   da = a.get_data();
    const std::vector<uint8_t>& ma = a.get_mask();
    const std::vector<float>&   db = b.get_data();
    const std::vector<uint8_t>& mb = b.get_mask();

    const int   n      = static_cast<int>(da.size());
    const float nodata = static_cast<float>(a.get_meta().no_data_value);

    std::vector<float>   out_data(n);
    std::vector<uint8_t> out_mask(n);

    #pragma omp parallel for schedule(static)
    for (int i = 0; i < n; ++i)
    {
        if (ma[i] && mb[i])
        {
            out_data[i] = da[i] - db[i];
            out_mask[i] = 1u;
        }
        else
        {
            out_data[i] = nodata;
            out_mask[i] = 0u;
        }
    }

    return lsdtt_raster_cpu(
        std::move(out_data), std::move(out_mask), a.get_meta()
    );
}

lsdtt_raster_cpu raster_multiply(
    const lsdtt_raster_cpu& a,
    const lsdtt_raster_cpu& b)
{
    assert_compatible(a, b);

    const std::vector<float>&   da = a.get_data();
    const std::vector<uint8_t>& ma = a.get_mask();
    const std::vector<float>&   db = b.get_data();
    const std::vector<uint8_t>& mb = b.get_mask();

    const int   n      = static_cast<int>(da.size());
    const float nodata = static_cast<float>(a.get_meta().no_data_value);

    std::vector<float>   out_data(n);
    std::vector<uint8_t> out_mask(n);

    #pragma omp parallel for schedule(static)
    for (int i = 0; i < n; ++i)
    {
        if (ma[i] && mb[i])
        {
            out_data[i] = da[i] * db[i];
            out_mask[i] = 1u;
        }
        else
        {
            out_data[i] = nodata;
            out_mask[i] = 0u;
        }
    }

    return lsdtt_raster_cpu(
        std::move(out_data), std::move(out_mask), a.get_meta()
    );
}

lsdtt_raster_cpu raster_divide(
    const lsdtt_raster_cpu& a,
    const lsdtt_raster_cpu& b)
{
    assert_compatible(a, b);

    const std::vector<float>&   da = a.get_data();
    const std::vector<uint8_t>& ma = a.get_mask();
    const std::vector<float>&   db = b.get_data();
    const std::vector<uint8_t>& mb = b.get_mask();

    const int   n      = static_cast<int>(da.size());
    const float nodata = static_cast<float>(a.get_meta().no_data_value);

    std::vector<float>   out_data(n);
    std::vector<uint8_t> out_mask(n);

    #pragma omp parallel for schedule(static)
    for (int i = 0; i < n; ++i)
    {
        if (ma[i] && mb[i] && db[i] != 0.0f)
        {
            out_data[i] = da[i] / db[i];
            out_mask[i] = 1u;
        }
        else
        {
            out_data[i] = nodata;
            out_mask[i] = 0u;
        }
    }

    return lsdtt_raster_cpu(
        std::move(out_data), std::move(out_mask), a.get_meta()
    );
}

lsdtt_raster_cpu raster_min(
    const lsdtt_raster_cpu& a,
    const lsdtt_raster_cpu& b)
{
    assert_compatible(a, b);

    const std::vector<float>&   da = a.get_data();
    const std::vector<uint8_t>& ma = a.get_mask();
    const std::vector<float>&   db = b.get_data();
    const std::vector<uint8_t>& mb = b.get_mask();

    const int   n      = static_cast<int>(da.size());
    const float nodata = static_cast<float>(a.get_meta().no_data_value);

    std::vector<float>   out_data(n);
    std::vector<uint8_t> out_mask(n);

    #pragma omp parallel for schedule(static)
    for (int i = 0; i < n; ++i)
    {
        if (ma[i] && mb[i])
        {
            out_data[i] = std::min(da[i], db[i]);
            out_mask[i] = 1u;
        }
        else
        {
            out_data[i] = nodata;
            out_mask[i] = 0u;
        }
    }

    return lsdtt_raster_cpu(
        std::move(out_data), std::move(out_mask), a.get_meta()
    );
}

lsdtt_raster_cpu raster_max_union(
    const lsdtt_raster_cpu& a,
    const lsdtt_raster_cpu& b)
{
    assert_compatible(a, b);

    const std::vector<float>&   da = a.get_data();
    const std::vector<uint8_t>& ma = a.get_mask();
    const std::vector<float>&   db = b.get_data();
    const std::vector<uint8_t>& mb = b.get_mask();

    const int   n      = static_cast<int>(da.size());
    const float nodata = static_cast<float>(a.get_meta().no_data_value);

    std::vector<float>   out_data(n);
    std::vector<uint8_t> out_mask(n);

    #pragma omp parallel for schedule(static)
    for (int i = 0; i < n; ++i)
    {
        if (ma[i] && mb[i])
        {
            out_data[i] = std::max(da[i], db[i]);
            out_mask[i] = 1u;
        }
        else if (ma[i])
        {
            out_data[i] = da[i];
            out_mask[i] = 1u;
        }
        else if (mb[i])
        {
            out_data[i] = db[i];
            out_mask[i] = 1u;
        }
        else
        {
            out_data[i] = nodata;
            out_mask[i] = 0u;
        }
    }

    return lsdtt_raster_cpu(
        std::move(out_data), std::move(out_mask), a.get_meta()
    );
}

lsdtt_raster_cpu raster_positive_difference(
    const lsdtt_raster_cpu& a,
    const lsdtt_raster_cpu& b)
{
    assert_compatible(a, b);

    const std::vector<float>&   da = a.get_data();
    const std::vector<uint8_t>& ma = a.get_mask();
    const std::vector<float>&   db = b.get_data();
    const std::vector<uint8_t>& mb = b.get_mask();

    const int   n      = static_cast<int>(da.size());
    const float nodata = static_cast<float>(a.get_meta().no_data_value);

    std::vector<float>   out_data(n);
    std::vector<uint8_t> out_mask(n);

    #pragma omp parallel for schedule(static)
    for (int i = 0; i < n; ++i)
    {
        if (ma[i] && mb[i])
        {
            const float diff = da[i] - db[i];
            out_data[i] = (diff > 0.0f) ? diff : 0.0f;
            out_mask[i] = 1u;
        }
        else
        {
            out_data[i] = nodata;
            out_mask[i] = 0u;
        }
    }

    return lsdtt_raster_cpu(
        std::move(out_data), std::move(out_mask), a.get_meta()
    );
}




// =============================================================================
// Unary transforms
// operations that take one input raster and produce one output raster by transforming each pixel value.
// =============================================================================
// raster_clamp_min: If any pixel is less than lower, it gets raised to the value of lower
lsdtt_raster_cpu raster_clamp_min(
    const lsdtt_raster_cpu& a,
    float                   lower)
{
    const std::vector<float>&   da = a.get_data();
    const std::vector<uint8_t>& ma = a.get_mask();

    const int   n      = static_cast<int>(da.size());
    const float nodata = static_cast<float>(a.get_meta().no_data_value);

    std::vector<float>   out_data(n);
    std::vector<uint8_t> out_mask(n);

    #pragma omp parallel for schedule(static)
    for (int i = 0; i < n; ++i)
    {
        if (ma[i])
        {
            out_data[i] = std::max(da[i], lower);
            out_mask[i] = 1u;
        }
        else
        {
            out_data[i] = nodata;
            out_mask[i] = 0u;
        }
    }

    return lsdtt_raster_cpu(
        std::move(out_data), std::move(out_mask), a.get_meta()
    );
}


// If any pixel is greater than upper, it gets lowered to the value of upper
lsdtt_raster_cpu raster_clamp_max(
    const lsdtt_raster_cpu& a,
    float                   upper)
{
    const std::vector<float>&   da = a.get_data();
    const std::vector<uint8_t>& ma = a.get_mask();

    const int   n      = static_cast<int>(da.size());
    const float nodata = static_cast<float>(a.get_meta().no_data_value);

    std::vector<float>   out_data(n);
    std::vector<uint8_t> out_mask(n);

    #pragma omp parallel for schedule(static)
    for (int i = 0; i < n; ++i)
    {
        if (ma[i])
        {
            out_data[i] = std::min(da[i], upper);
            out_mask[i] = 1u;
        }
        else
        {
            out_data[i] = nodata;
            out_mask[i] = 0u;
        }
    }

    return lsdtt_raster_cpu(
        std::move(out_data), std::move(out_mask), a.get_meta()
    );
}

lsdtt_raster_cpu raster_raise(
    const lsdtt_raster_cpu& a,
    float                   amount)
{
    const std::vector<float>&   da = a.get_data();
    const std::vector<uint8_t>& ma = a.get_mask();

    const int   n      = static_cast<int>(da.size());
    const float nodata = static_cast<float>(a.get_meta().no_data_value);

    std::vector<float>   out_data(n);
    std::vector<uint8_t> out_mask(n);

    #pragma omp parallel for schedule(static)
    for (int i = 0; i < n; ++i)
    {
        if (ma[i])
        {
            out_data[i] = da[i] + amount;
            out_mask[i] = 1u;
        }
        else
        {
            out_data[i] = nodata;
            out_mask[i] = 0u;
        }
    }

    return lsdtt_raster_cpu(
        std::move(out_data), std::move(out_mask), a.get_meta()
    );
}

lsdtt_raster_cpu raster_lower(
    const lsdtt_raster_cpu& a,
    float                   amount)
{
    const std::vector<float>&   da = a.get_data();
    const std::vector<uint8_t>& ma = a.get_mask();

    const int   n      = static_cast<int>(da.size());
    const float nodata = static_cast<float>(a.get_meta().no_data_value);

    std::vector<float>   out_data(n);
    std::vector<uint8_t> out_mask(n);

    #pragma omp parallel for schedule(static)
    for (int i = 0; i < n; ++i)
    {
        if (ma[i])
        {
            out_data[i] = da[i] - amount;
            out_mask[i] = 1u;
        }
        else
        {
            out_data[i] = nodata;
            out_mask[i] = 0u;
        }
    }

    return lsdtt_raster_cpu(
        std::move(out_data), std::move(out_mask), a.get_meta()
    );
}


lsdtt_raster_cpu raster_scale(
    const lsdtt_raster_cpu& a,
    float                   factor)
{
    const std::vector<float>&   da = a.get_data();
    const std::vector<uint8_t>& ma = a.get_mask();

    const int   n      = static_cast<int>(da.size());
    const float nodata = static_cast<float>(a.get_meta().no_data_value);

    std::vector<float>   out_data(n);
    std::vector<uint8_t> out_mask(n);

    #pragma omp parallel for schedule(static)
    for (int i = 0; i < n; ++i)
    {
        if (ma[i])
        {
            out_data[i] = da[i] * factor;
            out_mask[i] = 1u;
        }
        else
        {
            out_data[i] = nodata;
            out_mask[i] = 0u;
        }
    }

    return lsdtt_raster_cpu(
        std::move(out_data), std::move(out_mask), a.get_meta()
    );
}

lsdtt_raster_cpu raster_negate(
    const lsdtt_raster_cpu& a)
{
    const std::vector<float>&   da = a.get_data();
    const std::vector<uint8_t>& ma = a.get_mask();

    const int   n      = static_cast<int>(da.size());
    const float nodata = static_cast<float>(a.get_meta().no_data_value);

    std::vector<float>   out_data(n);
    std::vector<uint8_t> out_mask(n);

    #pragma omp parallel for schedule(static)
    for (int i = 0; i < n; ++i)
    {
        if (ma[i])
        {
            out_data[i] = -da[i];
            out_mask[i] = 1u;
        }
        else
        {
            out_data[i] = nodata;
            out_mask[i] = 0u;
        }
    }

    return lsdtt_raster_cpu(
        std::move(out_data), std::move(out_mask), a.get_meta()
    );
}

lsdtt_raster_cpu raster_scalar_subtract(
    float                   c,
    const lsdtt_raster_cpu& a)
{
    const std::vector<float>&   da = a.get_data();
    const std::vector<uint8_t>& ma = a.get_mask();

    const int   n      = static_cast<int>(da.size());
    const float nodata = static_cast<float>(a.get_meta().no_data_value);

    std::vector<float>   out_data(n);
    std::vector<uint8_t> out_mask(n);

    #pragma omp parallel for schedule(static)
    for (int i = 0; i < n; ++i)
    {
        if (ma[i])
        {
            out_data[i] = c - da[i];
            out_mask[i] = 1u;
        }
        else
        {
            out_data[i] = nodata;
            out_mask[i] = 0u;
        }
    }

    return lsdtt_raster_cpu(
        std::move(out_data), std::move(out_mask), a.get_meta()
    );
}

lsdtt_raster_cpu raster_clip(
    const lsdtt_raster_cpu& a,
    float                   lower,
    float                   upper)
{
    if (lower > upper)
    {
        throw std::invalid_argument(
            "lsdtt::algebra::raster_clip(): lower must be <= upper"
        );
    }

    const std::vector<float>&   da = a.get_data();
    const std::vector<uint8_t>& ma = a.get_mask();

    const int   n      = static_cast<int>(da.size());
    const float nodata = static_cast<float>(a.get_meta().no_data_value);

    std::vector<float>   out_data(n);
    std::vector<uint8_t> out_mask(n);

    #pragma omp parallel for schedule(static)
    for (int i = 0; i < n; ++i)
    {
        if (ma[i])
        {
            float v = da[i];
            if (v < lower)      v = lower;
            else if (v > upper) v = upper;
            out_data[i] = v;
            out_mask[i] = 1u;
        }
        else
        {
            out_data[i] = nodata;
            out_mask[i] = 0u;
        }
    }

    return lsdtt_raster_cpu(
        std::move(out_data), std::move(out_mask), a.get_meta()
    );
}

// =============================================================================
// Modification of bounds
// =============================================================================
lsdtt_raster_cpu raster_positive_part(
    const lsdtt_raster_cpu& a)
{
    const std::vector<float>&   da = a.get_data();
    const std::vector<uint8_t>& ma = a.get_mask();

    const int   n      = static_cast<int>(da.size());
    const float nodata = static_cast<float>(a.get_meta().no_data_value);

    std::vector<float>   out_data(n);
    std::vector<uint8_t> out_mask(n);

    #pragma omp parallel for schedule(static)
    for (int i = 0; i < n; ++i)
    {
        if (ma[i])
        {
            out_data[i] = std::max(da[i], 0.0f);
            out_mask[i] = 1u;
        }
        else
        {
            out_data[i] = nodata;
            out_mask[i] = 0u;
        }
    }

    return lsdtt_raster_cpu(
        std::move(out_data), std::move(out_mask), a.get_meta()
    );
}

lsdtt_raster_cpu raster_negative_part(
    const lsdtt_raster_cpu& a)
{
    const std::vector<float>&   da = a.get_data();
    const std::vector<uint8_t>& ma = a.get_mask();

    const int   n      = static_cast<int>(da.size());
    const float nodata = static_cast<float>(a.get_meta().no_data_value);

    std::vector<float>   out_data(n);
    std::vector<uint8_t> out_mask(n);

    #pragma omp parallel for schedule(static)
    for (int i = 0; i < n; ++i)
    {
        if (ma[i])
        {
            out_data[i] = std::max(-da[i], 0.0f);
            out_mask[i] = 1u;
        }
        else
        {
            out_data[i] = nodata;
            out_mask[i] = 0u;
        }
    }

    return lsdtt_raster_cpu(
        std::move(out_data), std::move(out_mask), a.get_meta()
    );
}

lsdtt_raster_cpu raster_normalise_unit_interval(
    const lsdtt_raster_cpu& a)
{
    const lsdtt_raster_meta& meta = a.get_meta();
    const std::vector<float>&   da = a.get_data();
    const std::vector<uint8_t>& ma = a.get_mask();

    const int   n      = static_cast<int>(da.size());
    const float nodata = static_cast<float>(meta.no_data_value);

    double vmin =  std::numeric_limits<double>::infinity();
    double vmax = -std::numeric_limits<double>::infinity();
    bool any_valid = false;

    for (int i = 0; i < n; ++i)
    {
        if (!ma[i]) continue;
        if (meta.is_nodata(da[i])) continue;
        if (!std::isfinite(da[i])) continue;

        const double v = static_cast<double>(da[i]);
        if (v < vmin) vmin = v;
        if (v > vmax) vmax = v;
        any_valid = true;
    }

    std::vector<float>   out_data(n, nodata);
    std::vector<uint8_t> out_mask(n, 0u);

    if (!any_valid)
    {
        return lsdtt_raster_cpu(
            std::move(out_data), std::move(out_mask), meta
        );
    }

    const double range = vmax - vmin;

    #pragma omp parallel for schedule(static)
    for (int i = 0; i < n; ++i)
    {
        if (!ma[i]) continue;
        if (meta.is_nodata(da[i])) continue;
        if (!std::isfinite(da[i])) continue;

        double out = 0.0;
        if (range > 0.0)
        {
            out = (static_cast<double>(da[i]) - vmin) / range;
        }

        if (out < 0.0) out = 0.0;
        if (out > 1.0) out = 1.0;

        out_data[i] = static_cast<float>(out);
        out_mask[i] = 1u;
    }

    return lsdtt_raster_cpu(
        std::move(out_data), std::move(out_mask), meta
    );
}


lsdtt_raster_cpu raster_normalise_nonnegative(
    const lsdtt_raster_cpu& a)
{
    const lsdtt_raster_meta& meta = a.get_meta();
    const std::vector<float>&   da = a.get_data();
    const std::vector<uint8_t>& ma = a.get_mask();

    const int   n      = static_cast<int>(da.size());
    const float nodata = static_cast<float>(meta.no_data_value);

    // Find the maximum of clamped-to-zero values over valid pixels.
    double vmax = 0.0;
    bool any_valid = false;

    for (int i = 0; i < n; ++i)
    {
        if (!ma[i]) continue;
        if (meta.is_nodata(da[i])) continue;
        if (!std::isfinite(da[i])) continue;

        const double v = static_cast<double>(da[i]);
        if (v > vmax) vmax = v;
        any_valid = true;
    }

    std::vector<float>   out_data(n, nodata);
    std::vector<uint8_t> out_mask(n, 0u);

    if (!any_valid)
    {
        return lsdtt_raster_cpu(
            std::move(out_data), std::move(out_mask), meta
        );
    }

    #pragma omp parallel for schedule(static)
    for (int i = 0; i < n; ++i)
    {
        if (!ma[i]) continue;
        if (meta.is_nodata(da[i])) continue;
        if (!std::isfinite(da[i])) continue;

        const double v = static_cast<double>(da[i]);
        double out = (v > 0.0) ? v : 0.0;
        if (vmax > 0.0) out /= vmax;

        out_data[i] = static_cast<float>(out);
        out_mask[i] = 1u;
    }

    return lsdtt_raster_cpu(
        std::move(out_data), std::move(out_mask), meta
    );
}


// =============================================================================
// Masking and conditional selection
// =============================================================================
lsdtt_raster_cpu raster_binary_true_above_threshold(
    const lsdtt_raster_cpu& a,
    float                   threshold)
{
    const lsdtt_raster_meta& meta = a.get_meta();
    const std::vector<float>&   da = a.get_data();
    const std::vector<uint8_t>& ma = a.get_mask();

    const int   n      = static_cast<int>(da.size());
    const float nodata = static_cast<float>(meta.no_data_value);

    std::vector<float>   out_data(n, nodata);
    std::vector<uint8_t> out_mask(n, 0u);

    #pragma omp parallel for schedule(static)
    for (int i = 0; i < n; ++i)
    {
        if (!ma[i]) continue;
        if (meta.is_nodata(da[i])) continue;
        if (!std::isfinite(da[i])) continue;

        out_data[i] = (da[i] >= threshold) ? 1.0f : 0.0f;
        out_mask[i] = 1u;
    }

    return lsdtt_raster_cpu(
        std::move(out_data), std::move(out_mask), meta
    );
}

lsdtt_raster_cpu raster_binary_true_above_but_excluding_threshold(
    const lsdtt_raster_cpu& a,
    float                   threshold)
{
    const lsdtt_raster_meta& meta = a.get_meta();
    const std::vector<float>&   da = a.get_data();
    const std::vector<uint8_t>& ma = a.get_mask();

    const int   n      = static_cast<int>(da.size());
    const float nodata = static_cast<float>(meta.no_data_value);

    std::vector<float>   out_data(n, nodata);
    std::vector<uint8_t> out_mask(n, 0u);

    #pragma omp parallel for schedule(static)
    for (int i = 0; i < n; ++i)
    {
        if (!ma[i]) continue;
        if (meta.is_nodata(da[i])) continue;
        if (!std::isfinite(da[i])) continue;

        out_data[i] = (da[i] > threshold) ? 1.0f : 0.0f;
        out_mask[i] = 1u;
    }

    return lsdtt_raster_cpu(
        std::move(out_data), std::move(out_mask), meta
    );
}

lsdtt_raster_cpu raster_binary_true_below_threshold(
    const lsdtt_raster_cpu& a,
    float                   threshold)
{
    const lsdtt_raster_meta& meta = a.get_meta();
    const std::vector<float>&   da = a.get_data();
    const std::vector<uint8_t>& ma = a.get_mask();

    const int   n      = static_cast<int>(da.size());
    const float nodata = static_cast<float>(meta.no_data_value);

    std::vector<float>   out_data(n, nodata);
    std::vector<uint8_t> out_mask(n, 0u);

    #pragma omp parallel for schedule(static)
    for (int i = 0; i < n; ++i)
    {
        if (!ma[i]) continue;
        if (meta.is_nodata(da[i])) continue;
        if (!std::isfinite(da[i])) continue;

        out_data[i] = (da[i] <= threshold) ? 1.0f : 0.0f;
        out_mask[i] = 1u;
    }

    return lsdtt_raster_cpu(
        std::move(out_data), std::move(out_mask), meta
    );
}

lsdtt_raster_cpu invert_binary_mask(
    const lsdtt_raster_cpu& a)
{
    const lsdtt_raster_meta& meta = a.get_meta();
    const std::vector<float>&   da = a.get_data();
    const std::vector<uint8_t>& ma = a.get_mask();

    const int   n      = static_cast<int>(da.size());
    const float nodata = static_cast<float>(meta.no_data_value);

    std::vector<float>   out_data(n, nodata);
    std::vector<uint8_t> out_mask(n, 0u);

    #pragma omp parallel for schedule(static)
    for (int i = 0; i < n; ++i)
    {
        if (!ma[i]) continue;
        out_data[i] = (da[i] > 0.5f) ? 0.0f : 1.0f;
        out_mask[i] = 1u;
    }

    return lsdtt_raster_cpu(
        std::move(out_data), std::move(out_mask), meta
    );
}



// keep the values raster only where keep_mask is true; set everything else to nodata.
lsdtt_raster_cpu raster_mask_by_binary_keep(
    const lsdtt_raster_cpu& values,
    const lsdtt_raster_cpu& keep_mask)
{
    assert_compatible(values, keep_mask);

    const lsdtt_raster_meta& meta = values.get_meta();

    const std::vector<float>&   v_data = values.get_data();
    const std::vector<uint8_t>& v_mask = values.get_mask();
    const std::vector<float>&   k_data = keep_mask.get_data();
    const std::vector<uint8_t>& k_mask = keep_mask.get_mask();

    const int   n      = static_cast<int>(v_data.size());
    const float nodata = static_cast<float>(meta.no_data_value);

    std::vector<float>   out_data(n, nodata);
    std::vector<uint8_t> out_mask(n, 0u);

    #pragma omp parallel for schedule(static)
    for (int i = 0; i < n; ++i)
    {
        const bool valid =
            v_mask[i] != 0 &&
            k_mask[i] != 0 &&
            !meta.is_nodata(v_data[i]) &&
            !meta.is_nodata(k_data[i]) &&
            std::isfinite(v_data[i]) &&
            std::isfinite(k_data[i]);

        if (!valid) continue;

        if (k_data[i] > 0.5f)
        {
            out_data[i] = v_data[i];
            out_mask[i] = 1u;
        }
    }

    return lsdtt_raster_cpu(
        std::move(out_data), std::move(out_mask), meta
    );
}


// keep the values raster only where the threshold raster 
// is less than or equal to a cutoff; set everything else to nodata.
// Take this raster of values, but only keep pixels where some other raster is not too large.
lsdtt_raster_cpu raster_mask_keep_where_less_or_eq(
    const lsdtt_raster_cpu& values,
    const lsdtt_raster_cpu& selector,
    float                   max_value)
{
    assert_compatible(values, selector);

    const lsdtt_raster_meta& meta = values.get_meta();

    const std::vector<float>&   v_data = values.get_data();
    const std::vector<uint8_t>& v_mask = values.get_mask();
    const std::vector<float>&   s_data = selector.get_data();
    const std::vector<uint8_t>& s_mask = selector.get_mask();

    const int   n      = static_cast<int>(v_data.size());
    const float nodata = static_cast<float>(meta.no_data_value);

    std::vector<float>   out_data(n, nodata);
    std::vector<uint8_t> out_mask(n, 0u);

    #pragma omp parallel for schedule(static)
    for (int i = 0; i < n; ++i)
    {
        const bool valid =
            v_mask[i] != 0 &&
            s_mask[i] != 0 &&
            !meta.is_nodata(v_data[i]) &&
            !meta.is_nodata(s_data[i]) &&
            std::isfinite(v_data[i]) &&
            std::isfinite(s_data[i]);

        if (!valid) continue;

        if (s_data[i] <= max_value)
        {
            out_data[i] = v_data[i];
            out_mask[i] = 1u;
        }
    }

    return lsdtt_raster_cpu(
        std::move(out_data), std::move(out_mask), meta
    );
}

// Take this raster of values, but only keep pixels where some other raster is not too small.
lsdtt_raster_cpu raster_mask_keep_where_greater_or_eq(
    const lsdtt_raster_cpu& values,
    const lsdtt_raster_cpu& selector,
    float                   min_value)
{
    assert_compatible(values, selector);

    const lsdtt_raster_meta& meta = values.get_meta();

    const std::vector<float>&   v_data = values.get_data();
    const std::vector<uint8_t>& v_mask = values.get_mask();
    const std::vector<float>&   s_data = selector.get_data();
    const std::vector<uint8_t>& s_mask = selector.get_mask();

    const int   n      = static_cast<int>(v_data.size());
    const float nodata = static_cast<float>(meta.no_data_value);

    std::vector<float>   out_data(n, nodata);
    std::vector<uint8_t> out_mask(n, 0u);

    #pragma omp parallel for schedule(static)
    for (int i = 0; i < n; ++i)
    {
        const bool valid =
            v_mask[i] != 0 &&
            s_mask[i] != 0 &&
            !meta.is_nodata(v_data[i]) &&
            !meta.is_nodata(s_data[i]) &&
            std::isfinite(v_data[i]) &&
            std::isfinite(s_data[i]);

        if (!valid) continue;

        if (s_data[i] >= min_value)
        {
            out_data[i] = v_data[i];
            out_mask[i] = 1u;
        }
    }

    return lsdtt_raster_cpu(
        std::move(out_data), std::move(out_mask), meta
    );
}


lsdtt_raster_cpu raster_binary_keep_true_where_less_or_eq(
    const lsdtt_raster_cpu& binary_mask,
    const lsdtt_raster_cpu& selector,
    float                   max_value)
{
    assert_compatible(binary_mask, selector);

    const lsdtt_raster_meta& meta = binary_mask.get_meta();

    const std::vector<float>&   b_data = binary_mask.get_data();
    const std::vector<uint8_t>& b_mask = binary_mask.get_mask();
    const std::vector<float>&   s_data = selector.get_data();
    const std::vector<uint8_t>& s_mask = selector.get_mask();

    const int   n      = static_cast<int>(b_data.size());
    const float nodata = static_cast<float>(meta.no_data_value);

    std::vector<float>   out_data(n, nodata);
    std::vector<uint8_t> out_mask(n, 0u);

    #pragma omp parallel for schedule(static)
    for (int i = 0; i < n; ++i)
    {
        const bool valid =
            b_mask[i] != 0 &&
            s_mask[i] != 0 &&
            !meta.is_nodata(b_data[i]) &&
            !meta.is_nodata(s_data[i]) &&
            std::isfinite(b_data[i]) &&
            std::isfinite(s_data[i]);

        if (!valid) continue;

        out_mask[i] = 1u;
        out_data[i] = (b_data[i] > 0.5f && s_data[i] <= max_value) ? 1.0f : 0.0f;
    }

    return lsdtt_raster_cpu(
        std::move(out_data), std::move(out_mask), meta
    );
}

lsdtt_raster_cpu raster_binary_keep_true_where_greater_or_eq(
    const lsdtt_raster_cpu& binary_mask,
    const lsdtt_raster_cpu& selector,
    float                   min_value)
{
    assert_compatible(binary_mask, selector);

    const lsdtt_raster_meta& meta = binary_mask.get_meta();

    const std::vector<float>&   b_data = binary_mask.get_data();
    const std::vector<uint8_t>& b_mask = binary_mask.get_mask();
    const std::vector<float>&   s_data = selector.get_data();
    const std::vector<uint8_t>& s_mask = selector.get_mask();

    const int   n      = static_cast<int>(b_data.size());
    const float nodata = static_cast<float>(meta.no_data_value);

    std::vector<float>   out_data(n, nodata);
    std::vector<uint8_t> out_mask(n, 0u);

    #pragma omp parallel for schedule(static)
    for (int i = 0; i < n; ++i)
    {
        const bool valid =
            b_mask[i] != 0 &&
            s_mask[i] != 0 &&
            !meta.is_nodata(b_data[i]) &&
            !meta.is_nodata(s_data[i]) &&
            std::isfinite(b_data[i]) &&
            std::isfinite(s_data[i]);

        if (!valid) continue;

        out_mask[i] = 1u;
        out_data[i] = (b_data[i] > 0.5f && s_data[i] >= min_value) ? 1.0f : 0.0f;
    }

    return lsdtt_raster_cpu(
        std::move(out_data), std::move(out_mask), meta
    );
}




// =============================================================================
// Absolute values
// =============================================================================
lsdtt_raster_cpu raster_abs(
    const lsdtt_raster_cpu& a)
{
    const std::vector<float>&   da = a.get_data();
    const std::vector<uint8_t>& ma = a.get_mask();

    const int   n      = static_cast<int>(da.size());
    const float nodata = static_cast<float>(a.get_meta().no_data_value);

    std::vector<float>   out_data(n);
    std::vector<uint8_t> out_mask(n);

    #pragma omp parallel for schedule(static)
    for (int i = 0; i < n; ++i)
    {
        if (ma[i])
        {
            out_data[i] = std::abs(da[i]);
            out_mask[i] = 1u;
        }
        else
        {
            out_data[i] = nodata;
            out_mask[i] = 0u;
        }
    }

    return lsdtt_raster_cpu(
        std::move(out_data), std::move(out_mask), a.get_meta()
    );
}


// =============================================================================
// Logarithm operations
// =============================================================================

lsdtt_raster_cpu raster_log(const lsdtt_raster_cpu& a)
{
    const std::vector<float>&   da = a.get_data();
    const std::vector<uint8_t>& ma = a.get_mask();

    const int   n      = static_cast<int>(da.size());
    const float nodata = static_cast<float>(a.get_meta().no_data_value);

    std::vector<float>   out_data(n);
    std::vector<uint8_t> out_mask(n);

    #pragma omp parallel for schedule(static)
    for (int i = 0; i < n; ++i)
    {
        if (ma[i] && da[i] > 0.0f)
        {
            out_data[i] = std::log(da[i]);
            out_mask[i] = 1u;
        }
        else
        {
            out_data[i] = nodata;
            out_mask[i] = 0u;
        }
    }

    return lsdtt_raster_cpu(
        std::move(out_data), std::move(out_mask), a.get_meta()
    );
}

lsdtt_raster_cpu raster_log10(const lsdtt_raster_cpu& a)
{
    const std::vector<float>&   da = a.get_data();
    const std::vector<uint8_t>& ma = a.get_mask();

    const int   n      = static_cast<int>(da.size());
    const float nodata = static_cast<float>(a.get_meta().no_data_value);

    std::vector<float>   out_data(n);
    std::vector<uint8_t> out_mask(n);

    #pragma omp parallel for schedule(static)
    for (int i = 0; i < n; ++i)
    {
        if (ma[i] && da[i] > 0.0f)
        {
            out_data[i] = std::log10(da[i]);
            out_mask[i] = 1u;
        }
        else
        {
            out_data[i] = nodata;
            out_mask[i] = 0u;
        }
    }

    return lsdtt_raster_cpu(
        std::move(out_data), std::move(out_mask), a.get_meta()
    );
}

// =============================================================================
// Scalar summary statistics
// =============================================================================

namespace {

// Collects all valid pixel values into a vector. Throws if none exist.
std::vector<float> collect_valid(const lsdtt_raster_cpu& a)
{
    const std::vector<float>&    da = a.get_data();
    const std::vector<uint8_t>&  ma = a.get_mask();
    const int n = static_cast<int>(da.size());

    std::vector<float> out;
    out.reserve(static_cast<std::size_t>(n));
    for (int i = 0; i < n; ++i)
    {
        if (ma[i]) out.push_back(da[i]);
    }

    if (out.empty())
    {
        throw std::runtime_error(
            "lsdtt::algebra: raster contains no valid pixels"
        );
    }

    return out;
}

} // anonymous namespace

float raster_min_value(const lsdtt_raster_cpu& a)
{
    const std::vector<float>&   da = a.get_data();
    const std::vector<uint8_t>& ma = a.get_mask();
    const int n = static_cast<int>(da.size());

    double min_val = std::numeric_limits<double>::infinity();
    std::int64_t n_valid = 0;

    #pragma omp parallel for schedule(static) reduction(+:n_valid) reduction(min:min_val)
    for (int i = 0; i < n; ++i)
    {
        if (!ma[i]) continue;
        const double v = static_cast<double>(da[i]);
        if (v < min_val) min_val = v;
        ++n_valid;
    }

    if (n_valid == 0)
    {
        throw std::runtime_error(
            "lsdtt::algebra::raster_min_value(): raster contains no valid pixels"
        );
    }

    return static_cast<float>(min_val);
}

float raster_max_value(const lsdtt_raster_cpu& a)
{
    const std::vector<float>&   da = a.get_data();
    const std::vector<uint8_t>& ma = a.get_mask();
    const int n = static_cast<int>(da.size());

    double max_val = -std::numeric_limits<double>::infinity();
    std::int64_t n_valid = 0;

    #pragma omp parallel for schedule(static) reduction(+:n_valid) reduction(max:max_val)
    for (int i = 0; i < n; ++i)
    {
        if (!ma[i]) continue;
        const double v = static_cast<double>(da[i]);
        if (v > max_val) max_val = v;
        ++n_valid;
    }

    if (n_valid == 0)
    {
        throw std::runtime_error(
            "lsdtt::algebra::raster_max_value(): raster contains no valid pixels"
        );
    }

    return static_cast<float>(max_val);
}

float raster_mean_value(const lsdtt_raster_cpu& a)
{
    const std::vector<float>&   da = a.get_data();
    const std::vector<uint8_t>& ma = a.get_mask();
    const int n = static_cast<int>(da.size());

    double sum = 0.0;
    std::int64_t n_valid = 0;

    #pragma omp parallel for schedule(static) reduction(+:n_valid,sum)
    for (int i = 0; i < n; ++i)
    {
        if (!ma[i]) continue;
        sum += static_cast<double>(da[i]);
        ++n_valid;
    }

    if (n_valid == 0)
    {
        throw std::runtime_error(
            "lsdtt::algebra::raster_mean_value(): raster contains no valid pixels"
        );
    }

    return static_cast<float>(sum / static_cast<double>(n_valid));
}

float raster_median_value(const lsdtt_raster_cpu& a)
{
    return raster_percentile(a, 50.0);
}

float raster_percentile(const lsdtt_raster_cpu& a, double p)
{
    if (p < 0.0 || p > 100.0)
    {
        throw std::invalid_argument(
            "lsdtt::algebra::raster_percentile(): p must be in [0, 100]"
        );
    }

    std::vector<float> values = collect_valid(a);
    return static_cast<float>(compute_quantile_interpolated_in_place(values, p / 100.0));
}

float raster_stddev_value(const lsdtt_raster_cpu& a)
{
    const std::vector<float>&   da = a.get_data();
    const std::vector<uint8_t>& ma = a.get_mask();
    const int n = static_cast<int>(da.size());

    double sum = 0.0;
    std::int64_t n_valid = 0;

    #pragma omp parallel for schedule(static) reduction(+:n_valid,sum)
    for (int i = 0; i < n; ++i)
    {
        if (!ma[i]) continue;
        sum += static_cast<double>(da[i]);
        ++n_valid;
    }

    if (n_valid == 0)
    {
        throw std::runtime_error(
            "lsdtt::algebra::raster_stddev_value(): raster contains no valid pixels"
        );
    }

    const double mean = sum / static_cast<double>(n_valid);
    double sq_sum = 0.0;

    #pragma omp parallel for schedule(static) reduction(+:sq_sum)
    for (int i = 0; i < n; ++i)
    {
        if (!ma[i]) continue;
        const double d = static_cast<double>(da[i]) - mean;
        sq_sum += d * d;
    }

    return static_cast<float>(std::sqrt(sq_sum / static_cast<double>(n_valid)));
}

float raster_range_value(const lsdtt_raster_cpu& a)
{
    const std::vector<float>&   da = a.get_data();
    const std::vector<uint8_t>& ma = a.get_mask();
    const int n = static_cast<int>(da.size());

    double min_val =  std::numeric_limits<double>::infinity();
    double max_val = -std::numeric_limits<double>::infinity();
    std::int64_t n_valid = 0;

    #pragma omp parallel for schedule(static) reduction(+:n_valid) reduction(min:min_val) reduction(max:max_val)
    for (int i = 0; i < n; ++i)
    {
        if (!ma[i]) continue;
        const double v = static_cast<double>(da[i]);
        if (v < min_val) min_val = v;
        if (v > max_val) max_val = v;
        ++n_valid;
    }

    if (n_valid == 0)
    {
        throw std::runtime_error(
            "lsdtt::algebra::raster_range_value(): raster contains no valid pixels"
        );
    }

    return static_cast<float>(max_val - min_val);
}

// =============================================================================
// raster_cdf_divergence_threshold
// =============================================================================
double raster_cdf_divergence_threshold(
    const lsdtt_raster_cpu& raster,
    bool                    verbose,
    bool                    write_cdf,
    lsdtt_data_loader*      loader,
    const std::string&      cdf_path)
{
    const lsdtt_raster_meta&     meta = raster.get_meta();
    const std::vector<float>&   da   = raster.get_data();
    const std::vector<uint8_t>& ma   = raster.get_mask();

    std::vector<double> vals;
    vals.reserve(da.size() / 4);
    for (std::size_t i = 0; i < da.size(); ++i)
    {
        if (!ma[i]) continue;
        const float v = da[i];
        if (!std::isfinite(v) || meta.is_nodata(v)) continue;
        vals.push_back(static_cast<double>(v));
    }

    if (vals.size() < 10)
        throw std::runtime_error(
            "lsdtt::algebra::raster_cdf_divergence_threshold: too few valid pixels");

    const lsdtt_smoothed_cdf                  scdf = compute_smoothed_cdf(vals, 500);
    const std::vector<lsdtt_distribution_fit> fits = fit_distributions(vals);

    const lsdtt_distribution_fit* best_fit = nullptr;
    for (const auto& f : fits)
    {
        if (!f.valid || !std::isfinite(f.ks_statistic)) continue;
        if (!best_fit || f.ks_statistic < best_fit->ks_statistic) best_fit = &f;
    }

    if (!best_fit)
        throw std::runtime_error(
            "lsdtt::algebra::raster_cdf_divergence_threshold: no valid distribution fit");

    double      max_diff = -std::numeric_limits<double>::infinity();
    std::size_t max_idx  = 0;
    for (std::size_t i = 0; i < scdf.bin_centres.size(); ++i)
    {
        const double fitted = eval_cdf(*best_fit, scdf.bin_centres[i]);
        if (!std::isfinite(fitted)) continue;
        const double diff = scdf.cdf[i] - fitted;
        if (diff > max_diff) { max_diff = diff; max_idx = i; }
    }

    const double threshold  = scdf.bin_centres[max_idx];
    const double percentile = scdf.cdf[max_idx] * 100.0;

    if (verbose)
    {
        std::cout << "  Best-fit distribution: "
                  << distribution_family_name(best_fit->family) << "\n"
                  << "  Max CDF divergence   : " << max_diff << "\n"
                  << "  Threshold percentile : " << percentile << "th\n";
    }

    if (write_cdf && loader && !cdf_path.empty())
        write_smoothed_cdf_csv(*loader, cdf_path, scdf, fits);

    return threshold;
}

// =============================================================================
// raster_adaptive_relief_threshold
// =============================================================================
lsdtt_raster_cpu raster_adaptive_relief_threshold(
    const lsdtt_raster_cpu& raster,
    bool                    verbose,
    bool                    write_cdf,
    lsdtt_data_loader*      loader,
    const std::string&      cdf_path)
{
    const double threshold = raster_cdf_divergence_threshold(
        raster, verbose, write_cdf, loader, cdf_path);
    return raster_binary_true_below_threshold(raster, static_cast<float>(threshold));
}

raster_difference_result raster_difference_stats(
    const lsdtt_raster_cpu& a,
    const lsdtt_raster_cpu& b
)
{
    const lsdtt_raster_meta& ma = a.get_meta();
    const lsdtt_raster_meta& mb = b.get_meta();

    if (ma.n_rows != mb.n_rows || ma.n_cols != mb.n_cols)
    {
        throw std::invalid_argument(
            "lsdtt::algebra::raster_difference_stats(): rasters differ in size"
        );
    }

    const std::vector<float>&        da = a.get_data();
    const std::vector<std::uint8_t>& xa = a.get_mask();
    const std::vector<float>&        db = b.get_data();
    const std::vector<std::uint8_t>& xb = b.get_mask();

    const std::size_t n = da.size();

    double sum_sq  = 0.0;
    double sum_abs = 0.0;
    double sum_err = 0.0;
    double max_abs = 0.0;
    std::size_t count = 0;

    for (std::size_t i = 0; i < n; ++i)
    {
        if (xa[i] == 0 || xb[i] == 0) continue;
        if (ma.is_nodata(da[i]) || mb.is_nodata(db[i])) continue;
        if (!std::isfinite(da[i]) || !std::isfinite(db[i])) continue;

        const double err = static_cast<double>(da[i]) - static_cast<double>(db[i]);

        sum_sq  += err * err;
        sum_abs += std::abs(err);
        sum_err += err;
        max_abs  = std::max(max_abs, std::abs(err));
        ++count;
    }

    if (count == 0) return raster_difference_result{};

    const double n_d = static_cast<double>(count);

    raster_difference_result s;
    s.count   = count;
    s.rmse    = std::sqrt(sum_sq / n_d);
    s.mae     = sum_abs / n_d;
    s.max_abs = max_abs;
    s.bias    = sum_err / n_d;
    return s;
}

} // namespace algebra
} // namespace lsdtt
