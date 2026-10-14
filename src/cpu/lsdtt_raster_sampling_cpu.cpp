// lsdtt_raster_sampling_cpu.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd and Rutger Vervoordeldonk
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html

#include "lsdtt_raster_sampling_cpu.hpp"

#include <cmath>

namespace lsdtt
{

bool raster_sample_bilinear(const lsdtt_raster_cpu& raster,
                            double                  x,
                            double                  y,
                            float&                  value)
{
    const lsdtt_raster_meta& meta = raster.get_meta();
    const lsdtt_array2d_f& data = raster.get_data_array();
    const lsdtt_array2d_u8& mask = raster.get_mask_array();

    if (x < meta.x_min || x > meta.x_max() ||
        y < meta.y_min || y > meta.y_max())
    {
        return false;
    }

    const double col_f = (x - meta.x_min) / meta.cell_size - 0.5;
    const double row_f = (meta.y_max() - y) / meta.cell_size - 0.5;

    if (col_f < 0.0 || row_f < 0.0 ||
        col_f > static_cast<double>(meta.n_cols - 1) ||
        row_f > static_cast<double>(meta.n_rows - 1))
    {
        return false;
    }

    const int col0 = static_cast<int>(std::floor(col_f));
    const int row0 = static_cast<int>(std::floor(row_f));

    // Edge fallback: nearest-neighbour if a full 2x2 bilinear stencil
    // is not available.
    if (col0 >= meta.n_cols - 1 || row0 >= meta.n_rows - 1)
    {
        const int col_nn = static_cast<int>(std::llround(col_f));
        const int row_nn = static_cast<int>(std::llround(row_f));

        if (row_nn < 0 || row_nn >= meta.n_rows ||
            col_nn < 0 || col_nn >= meta.n_cols)
        {
            return false;
        }

        if (mask(row_nn, col_nn) == 0) return false;

        value = data(row_nn, col_nn);
        return true;
    }

    const int row1 = row0 + 1;
    const int col1 = col0 + 1;

    if (mask(row0, col0) == 0 ||
        mask(row0, col1) == 0 ||
        mask(row1, col0) == 0 ||
        mask(row1, col1) == 0)
    {
        return false;
    }

    const double dc = col_f - static_cast<double>(col0);
    const double dr = row_f - static_cast<double>(row0);

    const double z00 = static_cast<double>(data(row0, col0));
    const double z01 = static_cast<double>(data(row0, col1));
    const double z10 = static_cast<double>(data(row1, col0));
    const double z11 = static_cast<double>(data(row1, col1));

    const double z0 = (1.0 - dc) * z00 + dc * z01;
    const double z1 = (1.0 - dc) * z10 + dc * z11;
    const double z  = (1.0 - dr) * z0  + dr * z1;

    value = static_cast<float>(z);
    return true;
}

} // namespace lsdtt
