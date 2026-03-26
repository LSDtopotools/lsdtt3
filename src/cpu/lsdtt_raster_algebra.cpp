// lsdtt_raster_algebra.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2025 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html

#include "lsdtt_raster_algebra.hpp"

#include <algorithm>
#include <cmath>
#include <sstream>
#include <stdexcept>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace lsdtt {
namespace algebra {

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
            "lsdtt raster algebra: incompatible rasters —\n" + err.str()
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
// =============================================================================

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

} // namespace algebra
} // namespace lsdtt
