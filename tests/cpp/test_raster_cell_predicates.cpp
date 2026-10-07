// test_raster_cell_predicates.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
//
// -----------------------------------------------------------------------------
// Characterisation test for src/cpu/lsdtt_raster_cell_predicates_cpu.hpp.
//
// The shared predicates replaced private copies in five modules. Each removed
// copy is kept below VERBATIM, with the file and line it came from, as the
// reference: the shared function must give the same answer as every copy on
// every input, including nodata (finite and NaN nodata values), NaN, +-inf,
// -0.0, denormals, values either side of the 0.5 foreground threshold, and
// mask bytes 0, 1 and 255. Index helpers are compared over in- and out-of-range
// rows and columns, including grids large enough to overflow int.
//
// The private is_valid_scalar_ in lsdtt_tidal_channel_cpu.cpp:39-48 was never
// called; it was deleted, not migrated, so it has no reference here.
//
// Usage: lsdtt-test-raster-cell-predicates
// Prints one line per check and a final "RESULT: PASS|FAIL".
// -----------------------------------------------------------------------------

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <limits>
#include <random>
#include <string>
#include <vector>

#include "../../src/core/lsdtt_raster_meta.hpp"
#include "../../src/cpu/lsdtt_raster_cell_predicates_cpu.hpp"
#include "../../src/utils/lsdtt_array2d.hpp"

using namespace lsdtt;

namespace
{

// ---- removed copies (references) --------------------------------------------

// src/cpu/lsdtt_marsh_platform_cpu.cpp:43-73 (as of 2026-10-07)
namespace marsh
{
inline bool is_valid_scalar_(
    const lsdtt_raster_meta&       meta,
    const std::vector<float>&      data,
    const std::vector<uint8_t>&    mask,
    std::size_t                    idx)
{
    return mask[idx] != 0 &&
           !meta.is_nodata(data[idx]) &&
           std::isfinite(data[idx]);
}

inline bool is_binary_fg_(
    const lsdtt_raster_meta&       meta,
    const std::vector<float>&      data,
    const std::vector<uint8_t>&    mask,
    std::size_t                    idx)
{
    return is_valid_scalar_(meta, data, mask, idx) && data[idx] > 0.5f;
}

inline std::size_t idx_of_(int row, int col, int ncols)
{
    return static_cast<std::size_t>(row) * static_cast<std::size_t>(ncols) +
           static_cast<std::size_t>(col);
}

inline bool in_bounds_(int row, int col, int nrows, int ncols)
{
    return row >= 0 && row < nrows && col >= 0 && col < ncols;
}
} // namespace marsh

// src/cpu/lsdtt_tidal_skeleton_metrics_cpu.cpp:42-62 (as of 2026-10-07)
namespace tidal_skeleton
{
inline bool is_valid_scalar_(
    const lsdtt_raster_meta&    meta,
    const std::vector<float>&   data,
    const std::vector<uint8_t>& mask,
    std::size_t                 idx)
{
    return mask[idx] != 0 &&
           !meta.is_nodata(data[idx]) &&
           std::isfinite(data[idx]);
}

inline std::size_t idx_of_(int row, int col, int ncols)
{
    return static_cast<std::size_t>(row) * static_cast<std::size_t>(ncols) +
           static_cast<std::size_t>(col);
}

inline bool in_bounds_(int row, int col, int nrows, int ncols)
{
    return row >= 0 && row < nrows && col >= 0 && col < ncols;
}
} // namespace tidal_skeleton

// src/cpu/lsdtt_terrace_edge_tracer_cpu.cpp:39-48 (as of 2026-10-07)
namespace terrace
{
static bool is_valid_scalar_(
    const lsdtt_raster_meta& meta,
    const std::vector<float>& data,
    const std::vector<uint8_t>& mask,
    std::size_t idx)
{
    return mask[idx] != 0 &&
           !meta.is_nodata(data[idx]) &&
           std::isfinite(data[idx]);
}
} // namespace terrace

// src/cpu/lsdtt_raster_coregistration_cpu.cpp:155-168 (as of 2026-10-07)
namespace coregistration
{
static bool valid_cell_value(
    const lsdtt_raster_meta& meta,
    const lsdtt_array2d_f& data,
    const lsdtt_array2d_u8& mask,
    int row,
    int col)
{
    if (row < 0 || col < 0 || row >= meta.n_rows || col >= meta.n_cols) return false;
    if (mask(row, col) == 0u) return false;
    const float v = data(row, col);
    if (!std::isfinite(v)) return false;
    if (meta.is_nodata(v)) return false;
    return true;
}
} // namespace coregistration

// The body coregistration now uses (lsdtt_raster_coregistration_cpu.cpp).
static bool coregistration_new(
    const lsdtt_raster_meta& meta,
    const lsdtt_array2d_f& data,
    const lsdtt_array2d_u8& mask,
    int row,
    int col)
{
    return raster_in_bounds(row, col, meta.n_rows, meta.n_cols) &&
           raster_value_is_valid(meta, mask(row, col), data(row, col));
}

// ---- harness -----------------------------------------------------------------

int g_failures = 0;

void check(const std::string& what, bool ok, long long n_cases)
{
    std::cout << (ok ? "PASS " : "FAIL ") << what << " (" << n_cases << " cases)\n";
    if (!ok) ++g_failures;
}

// Values chosen to hit every branch of the predicates.
std::vector<float> special_values(float nodata)
{
    const float inf = std::numeric_limits<float>::infinity();
    return {
        nodata, std::numeric_limits<float>::quiet_NaN(), inf, -inf,
        0.0f, -0.0f, 0.5f, std::nextafter(0.5f, 1.0f), std::nextafter(0.5f, 0.0f),
        1.0f, -1.0f, std::numeric_limits<float>::denorm_min(),
        std::numeric_limits<float>::max(), std::numeric_limits<float>::lowest(),
        -9999.0f, 255.0f,
    };
}

struct raster_case
{
    lsdtt_raster_meta    meta;
    std::vector<float>   data;
    std::vector<uint8_t> mask;
};

raster_case make_case(double nodata, unsigned seed)
{
    raster_case rc;
    rc.meta.n_rows        = 37;
    rc.meta.n_cols        = 53;
    rc.meta.cell_size     = 2.0;
    rc.meta.no_data_value = nodata;
    const std::size_t n = static_cast<std::size_t>(rc.meta.n_rows) * rc.meta.n_cols;

    std::mt19937 rng(seed);
    std::uniform_real_distribution<float> uni(-2.0f, 2.0f);
    std::uniform_int_distribution<int>    pick(0, 3);
    const std::vector<float> specials = special_values(static_cast<float>(nodata));
    std::uniform_int_distribution<std::size_t> pick_special(0, specials.size() - 1);
    const uint8_t mask_bytes[3] = {0, 1, 255};
    std::uniform_int_distribution<int> pick_mask(0, 2);

    rc.data.resize(n);
    rc.mask.resize(n);
    for (std::size_t i = 0; i < n; ++i)
    {
        rc.data[i] = (pick(rng) == 0) ? specials[pick_special(rng)] : uni(rng);
        rc.mask[i] = mask_bytes[pick_mask(rng)];
    }
    // Every special value with every mask byte, so no combination is left to chance.
    std::size_t i = 0;
    for (float v : specials)
        for (uint8_t m : mask_bytes)
        {
            rc.data[i] = v;
            rc.mask[i] = m;
            ++i;
        }
    return rc;
}

} // namespace

int main()
{
    std::cout << "== raster cell predicates vs the removed private copies ==\n";

    const double nodatas[] = {-9999.0, std::numeric_limits<double>::quiet_NaN(), 0.0, 1.0};
    unsigned seed = 11;
    for (double nd : nodatas)
    {
        const raster_case rc = make_case(nd, seed++);
        const std::size_t n = rc.data.size();
        bool ok_valid = true, ok_fg = true, ok_coreg = true;

        lsdtt_array2d_f  data2d(rc.meta.n_rows, rc.meta.n_cols, std::vector<float>(rc.data));
        lsdtt_array2d_u8 mask2d(rc.meta.n_rows, rc.meta.n_cols, std::vector<uint8_t>(rc.mask));

        for (std::size_t i = 0; i < n; ++i)
        {
            const bool v_new = raster_cell_is_valid(rc.meta, rc.data, rc.mask, i);
            ok_valid = ok_valid &&
                v_new == marsh::is_valid_scalar_(rc.meta, rc.data, rc.mask, i) &&
                v_new == tidal_skeleton::is_valid_scalar_(rc.meta, rc.data, rc.mask, i) &&
                v_new == terrace::is_valid_scalar_(rc.meta, rc.data, rc.mask, i) &&
                v_new == raster_value_is_valid(rc.meta, rc.mask[i], rc.data[i]);
            ok_fg = ok_fg &&
                raster_cell_is_foreground(rc.meta, rc.data, rc.mask, i) ==
                marsh::is_binary_fg_(rc.meta, rc.data, rc.mask, i);
        }
        long long n_coreg = 0;
        for (int r = -2; r < rc.meta.n_rows + 2; ++r)
            for (int c = -2; c < rc.meta.n_cols + 2; ++c, ++n_coreg)
                ok_coreg = ok_coreg &&
                    coregistration_new(rc.meta, data2d, mask2d, r, c) ==
                    coregistration::valid_cell_value(rc.meta, data2d, mask2d, r, c);

        const std::string tag = "nodata=" + std::to_string(nd);
        check("raster_cell_is_valid == marsh/tidal_skeleton/terrace copies, " + tag,
              ok_valid, static_cast<long long>(n));
        check("raster_cell_is_foreground == marsh is_binary_fg_, " + tag,
              ok_fg, static_cast<long long>(n));
        check("coregistration body == valid_cell_value (incl. out of bounds), " + tag,
              ok_coreg, n_coreg);
    }

    // Index helpers: in and out of range, and a grid large enough that
    // row * n_cols overflows int.
    {
        bool ok_idx = true, ok_bounds = true;
        long long n_cases = 0;
        const int n_cols_list[] = {1, 7, 53, 70001};
        const int n_rows_list[] = {1, 5, 37, 70001};
        for (int nr : n_rows_list)
            for (int nc : n_cols_list)
            {
                const int rows[] = {-3, -1, 0, 1, nr / 2, nr - 1, nr, nr + 2};
                const int cols[] = {-3, -1, 0, 1, nc / 2, nc - 1, nc, nc + 2};
                for (int r : rows)
                    for (int c : cols)
                    {
                        ++n_cases;
                        const std::size_t idx = raster_flat_index(r, c, nc);
                        ok_idx = ok_idx &&
                            idx == marsh::idx_of_(r, c, nc) &&
                            idx == tidal_skeleton::idx_of_(r, c, nc);
                        const bool ib = raster_in_bounds(r, c, nr, nc);
                        ok_bounds = ok_bounds &&
                            ib == marsh::in_bounds_(r, c, nr, nc) &&
                            ib == tidal_skeleton::in_bounds_(r, c, nr, nc);
                    }
            }
        check("raster_flat_index == marsh/tidal_skeleton idx_of_", ok_idx, n_cases);
        check("raster_in_bounds == marsh/tidal_skeleton in_bounds_", ok_bounds, n_cases);
        check("raster_flat_index does not overflow int (70000 x 70001 grid)",
              raster_flat_index(70000, 70000, 70001) ==
                  std::size_t{70000} * std::size_t{70001} + std::size_t{70000}, 1);
    }

    std::cout << "RESULT: " << (g_failures == 0 ? "PASS" : "FAIL") << "\n";
    return g_failures == 0 ? 0 : 1;
}
