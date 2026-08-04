// lsdtt_raster_coregistration_cpu.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html

#include "lsdtt_raster_coregistration_cpu.hpp"

#include "../utils/lsdtt_linear_algebra.hpp"
#include "../utils/lsdtt_raster_geometry.hpp"
#include "../utils/lsdtt_statistics.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <random>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace lsdtt {

namespace {

struct sampled_value
{
    float value = 0.0f;
    bool valid = false;
};

static sampled_value bilinear_sample_raster(
    const lsdtt_raster_cpu& raster,
    double x,
    double y)
{
    const lsdtt_raster_meta& meta = raster.get_meta();
    const lsdtt_array2d_f& data = raster.get_data_array();
    const lsdtt_array2d_u8& mask = raster.get_mask_array();

    sampled_value out;

    if (x < meta.x_min || x > meta.x_max() || y < meta.y_min || y > meta.y_max())
    {
        return out;
    }

    const double col_f = (x - meta.x_min) / meta.cell_size - 0.5;
    const double row_f = (meta.y_max() - y) / meta.cell_size - 0.5;

    if (col_f < 0.0 || row_f < 0.0 ||
        col_f > static_cast<double>(meta.n_cols - 1) ||
        row_f > static_cast<double>(meta.n_rows - 1))
    {
        return out;
    }

    const int col0 = static_cast<int>(std::floor(col_f));
    const int row0 = static_cast<int>(std::floor(row_f));

    if (col0 >= meta.n_cols - 1 || row0 >= meta.n_rows - 1)
    {
        const int col_nn = static_cast<int>(std::llround(col_f));
        const int row_nn = static_cast<int>(std::llround(row_f));

        if (row_nn < 0 || row_nn >= meta.n_rows ||
            col_nn < 0 || col_nn >= meta.n_cols)
        {
            return out;
        }

        if (mask(row_nn, col_nn) == 0) return out;

        out.value = data(row_nn, col_nn);
        out.valid = std::isfinite(out.value);
        return out;
    }

    const int row1 = row0 + 1;
    const int col1 = col0 + 1;

    if (mask(row0, col0) == 0 ||
        mask(row0, col1) == 0 ||
        mask(row1, col0) == 0 ||
        mask(row1, col1) == 0)
    {
        return out;
    }

    const double dc = col_f - static_cast<double>(col0);
    const double dr = row_f - static_cast<double>(row0);

    const double z00 = static_cast<double>(data(row0, col0));
    const double z01 = static_cast<double>(data(row0, col1));
    const double z10 = static_cast<double>(data(row1, col0));
    const double z11 = static_cast<double>(data(row1, col1));

    const double z0 = (1.0 - dc) * z00 + dc * z01;
    const double z1 = (1.0 - dc) * z10 + dc * z11;
    const double z = (1.0 - dr) * z0 + dr * z1;

    out.value = static_cast<float>(z);
    out.valid = std::isfinite(out.value);
    return out;
}

static void validate_nuth_kaab_options(const nuth_kaab_coreg_options& opt)
{
    if (opt.nuth_kaab_coreg_max_iterations < 1)
    {
        throw std::invalid_argument("nuth_kaab_coreg_max_iterations must be >= 1");
    }
    if (opt.nuth_kaab_coreg_offset_threshold_pixels <= 0.0)
    {
        throw std::invalid_argument("nuth_kaab_coreg_offset_threshold_pixels must be > 0");
    }
    if (opt.nuth_kaab_coreg_slope_min_degrees < 0.0 ||
        opt.nuth_kaab_coreg_slope_min_degrees > 90.0)
    {
        throw std::invalid_argument("nuth_kaab_coreg_slope_min_degrees must be in [0,90]");
    }
    if (opt.nuth_kaab_coreg_slope_max_degrees <= opt.nuth_kaab_coreg_slope_min_degrees ||
        opt.nuth_kaab_coreg_slope_max_degrees > 90.0)
    {
        throw std::invalid_argument("nuth_kaab_coreg_slope_max_degrees must be in (slope_min,90]");
    }
    if (opt.nuth_kaab_coreg_aspect_bin_count < 8)
    {
        throw std::invalid_argument("nuth_kaab_coreg_aspect_bin_count must be >= 8");
    }
    if (opt.nuth_kaab_coreg_min_samples < 25)
    {
        throw std::invalid_argument("nuth_kaab_coreg_min_samples must be >= 25");
    }
    if (opt.nuth_kaab_coreg_subsample_count == 0)
    {
        throw std::invalid_argument("nuth_kaab_coreg_subsample_count must be -1 or >= 1");
    }
}

static void solve_weighted_3param_lsq(
    const std::vector<double>& c0,
    const std::vector<double>& c1,
    const std::vector<double>& c2,
    const std::vector<double>& y,
    double& p0,
    double& p1,
    double& p2);

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

static void compute_slope_tan_and_aspect(
    const lsdtt_raster_cpu& reference,
    std::vector<double>& slope_tan,
    std::vector<double>& aspect)
{
    const lsdtt_raster_meta& meta = reference.get_meta();
    const lsdtt_array2d_f& dem = reference.get_data_array();
    const lsdtt_array2d_u8& mask = reference.get_mask_array();

    const std::size_t n = static_cast<std::size_t>(meta.n_rows) *
                          static_cast<std::size_t>(meta.n_cols);
    slope_tan.assign(n, std::numeric_limits<double>::quiet_NaN());
    aspect.assign(n, std::numeric_limits<double>::quiet_NaN());

    const double two_pi = 2.0 * M_PI;

    for (int row = 0; row < meta.n_rows; ++row)
    {
        for (int col = 0; col < meta.n_cols; ++col)
        {
            if (!valid_cell_value(meta, dem, mask, row, col)) continue;

            double gx = std::numeric_limits<double>::quiet_NaN();
            double gy = std::numeric_limits<double>::quiet_NaN();

            if (col == 0)
            {
                if (!valid_cell_value(meta, dem, mask, row, col + 1)) continue;
                gx = static_cast<double>(dem(row, col + 1)) - static_cast<double>(dem(row, col));
            }
            else if (col == meta.n_cols - 1)
            {
                if (!valid_cell_value(meta, dem, mask, row, col - 1)) continue;
                gx = static_cast<double>(dem(row, col)) - static_cast<double>(dem(row, col - 1));
            }
            else
            {
                if (!valid_cell_value(meta, dem, mask, row, col + 1) ||
                    !valid_cell_value(meta, dem, mask, row, col - 1))
                {
                    continue;
                }
                gx = (static_cast<double>(dem(row, col + 1)) -
                      static_cast<double>(dem(row, col - 1))) * 0.5;
            }

            if (row == 0)
            {
                if (!valid_cell_value(meta, dem, mask, row + 1, col)) continue;
                gy = static_cast<double>(dem(row + 1, col)) - static_cast<double>(dem(row, col));
            }
            else if (row == meta.n_rows - 1)
            {
                if (!valid_cell_value(meta, dem, mask, row - 1, col)) continue;
                gy = static_cast<double>(dem(row, col)) - static_cast<double>(dem(row - 1, col));
            }
            else
            {
                if (!valid_cell_value(meta, dem, mask, row + 1, col) ||
                    !valid_cell_value(meta, dem, mask, row - 1, col))
                {
                    continue;
                }
                gy = (static_cast<double>(dem(row + 1, col)) -
                      static_cast<double>(dem(row - 1, col))) * 0.5;
            }

            const double s = std::sqrt(gx * gx + gy * gy);
            if (s <= 0.0 || !std::isfinite(s)) continue;

            double a = std::atan2(-gx, gy) + M_PI;
            if (a < 0.0) a += two_pi;
            if (a >= two_pi) a -= two_pi;

            const std::size_t idx = static_cast<std::size_t>(row) *
                                    static_cast<std::size_t>(meta.n_cols) +
                                    static_cast<std::size_t>(col);
            slope_tan[idx] = s;
            aspect[idx] = a;
        }
    }
}

static void fit_nuth_kaab_cosine_binned(
    const std::vector<double>& aspect,
    const std::vector<double>& y,
    int bin_count,
    double& easting_offset_pixels,
    double& northing_offset_pixels,
    double& vertical_term)
{
    std::vector<std::vector<double>> bins(static_cast<std::size_t>(bin_count));
    const double two_pi = 2.0 * M_PI;

    for (std::size_t i = 0; i < y.size(); ++i)
    {
        const double a = aspect[i];
        const double v = y[i];
        if (!std::isfinite(a) || !std::isfinite(v)) continue;
        int bi = static_cast<int>(std::floor((a / two_pi) * static_cast<double>(bin_count)));
        if (bi < 0) bi = 0;
        if (bi >= bin_count) bi = bin_count - 1;
        bins[static_cast<std::size_t>(bi)].push_back(v);
    }

    std::vector<double> x_fit;
    std::vector<double> y_fit;
    x_fit.reserve(static_cast<std::size_t>(bin_count));
    y_fit.reserve(static_cast<std::size_t>(bin_count));

    for (int bi = 0; bi < bin_count; ++bi)
    {
        const auto& bv = bins[static_cast<std::size_t>(bi)];
        if (bv.empty()) continue;
        const double x_center = (static_cast<double>(bi) + 0.5) * (two_pi / static_cast<double>(bin_count));
        x_fit.push_back(x_center);
        y_fit.push_back(compute_median(bv));
    }

    if (x_fit.size() < 3)
    {
        throw std::runtime_error("Nuth-Kaab coreg: insufficient populated aspect bins for cosine fit");
    }

    std::vector<double> c0(x_fit.size());
    std::vector<double> c1(x_fit.size());
    std::vector<double> c2(x_fit.size(), 1.0);
    for (std::size_t i = 0; i < x_fit.size(); ++i)
    {
        c0[i] = std::cos(x_fit[i]);
        c1[i] = std::sin(x_fit[i]);
    }

    double A = 0.0;
    double B = 0.0;
    double C = 0.0;
    solve_weighted_3param_lsq(c0, c1, c2, y_fit, A, B, C);

    northing_offset_pixels = A;
    easting_offset_pixels = B;
    vertical_term = C;
}

static void solve_weighted_3param_lsq(
    const std::vector<double>& c0,
    const std::vector<double>& c1,
    const std::vector<double>& c2,
    const std::vector<double>& y,
    double& p0,
    double& p1,
    double& p2)
{
    if (c0.size() != c1.size() || c0.size() != c2.size() || c0.size() != y.size())
    {
        throw std::invalid_argument("internal coreg solve: size mismatch");
    }

    double ata[3][3] = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
    double aty[3] = {0.0, 0.0, 0.0};

    for (std::size_t i = 0; i < y.size(); ++i)
    {
        const double a0 = c0[i];
        const double a1 = c1[i];
        const double a2 = c2[i];
        const double yi = y[i];

        ata[0][0] += a0 * a0;
        ata[0][1] += a0 * a1;
        ata[0][2] += a0 * a2;
        ata[1][1] += a1 * a1;
        ata[1][2] += a1 * a2;
        ata[2][2] += a2 * a2;

        aty[0] += a0 * yi;
        aty[1] += a1 * yi;
        aty[2] += a2 * yi;
    }

    ata[1][0] = ata[0][1];
    ata[2][0] = ata[0][2];
    ata[2][1] = ata[1][2];

    double inv[3][3] = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
    if (!invert3x3(ata, inv))
    {
        throw std::runtime_error("Nuth-Kaab coreg: singular normal matrix");
    }

    p0 = inv[0][0] * aty[0] + inv[0][1] * aty[1] + inv[0][2] * aty[2];
    p1 = inv[1][0] * aty[0] + inv[1][1] * aty[1] + inv[1][2] * aty[2];
    p2 = inv[2][0] * aty[0] + inv[2][1] * aty[1] + inv[2][2] * aty[2];
}

static lsdtt_raster_cpu apply_shift_to_moving(
    const lsdtt_raster_cpu& moving,
    double dx_m,
    double dy_m,
    double dz_m)
{
    const lsdtt_raster_meta& meta = moving.get_meta();
    const float nodata = static_cast<float>(meta.no_data_value);

    std::vector<float> out_data(static_cast<std::size_t>(meta.n_rows) *
                                static_cast<std::size_t>(meta.n_cols), nodata);
    std::vector<uint8_t> out_mask(out_data.size(), 0u);

    for (int row = 0; row < meta.n_rows; ++row)
    {
        for (int col = 0; col < meta.n_cols; ++col)
        {
            const lsdtt_xy xy = cell_center_from_row_col(meta, row, col);
            const sampled_value s = bilinear_sample_raster(moving, xy.x + dx_m, xy.y + dy_m);
            if (!s.valid) continue;

            const std::size_t idx = static_cast<std::size_t>(row) *
                                    static_cast<std::size_t>(meta.n_cols) +
                                    static_cast<std::size_t>(col);
            out_data[idx] = static_cast<float>(static_cast<double>(s.value) + dz_m);
            out_mask[idx] = 1u;
        }
    }

    return lsdtt_raster_cpu(std::move(out_data), std::move(out_mask), meta);
}

static dem_coreg_result run_nuth_kaab(
    const lsdtt_raster_cpu& reference,
    const lsdtt_raster_cpu& moving,
    const nuth_kaab_coreg_options& opt)
{
    validate_nuth_kaab_options(opt);

    const lsdtt_raster_meta& ref_meta = reference.get_meta();
    const lsdtt_raster_meta& mov_meta = moving.get_meta();

    if (!same_grid(ref_meta, mov_meta))
    {
        throw std::runtime_error(
            "Nuth-Kaab coreg requires reference and moving rasters on the same grid");
    }

    ref_meta.check_projected_metres();
    mov_meta.check_projected_metres();

    const std::vector<float>& ref_data = reference.get_data();
    const std::vector<uint8_t>& ref_mask = reference.get_mask();

    std::vector<double> slope_tan;
    std::vector<double> aspect;
    compute_slope_tan_and_aspect(reference, slope_tan, aspect);

    std::vector<std::size_t> candidate_indices;
    candidate_indices.reserve(ref_data.size());

    const double deg_per_rad = 180.0 / M_PI;
    for (std::size_t idx = 0; idx < ref_data.size(); ++idx)
    {
        if (ref_mask[idx] == 0u || !std::isfinite(ref_data[idx])) continue;
        if (!std::isfinite(slope_tan[idx]) || !std::isfinite(aspect[idx])) continue;

        const double slope_deg = std::atan(slope_tan[idx]) * deg_per_rad;
        if (slope_deg < opt.nuth_kaab_coreg_slope_min_degrees ||
            slope_deg > opt.nuth_kaab_coreg_slope_max_degrees)
        {
            continue;
        }

        candidate_indices.push_back(idx);
    }

    if (static_cast<int>(candidate_indices.size()) < opt.nuth_kaab_coreg_min_samples)
    {
        throw std::runtime_error("Nuth-Kaab coreg: too few candidate samples before iteration");
    }

    if (opt.nuth_kaab_coreg_subsample_count > 0 &&
        static_cast<int>(candidate_indices.size()) > opt.nuth_kaab_coreg_subsample_count)
    {
        std::mt19937 rng(42u);
        std::shuffle(candidate_indices.begin(), candidate_indices.end(), rng);
        candidate_indices.resize(static_cast<std::size_t>(opt.nuth_kaab_coreg_subsample_count));
    }

    dem_coreg_result result;

    double dx_m = opt.nuth_kaab_coreg_initial_dx_m;
    double dy_m = opt.nuth_kaab_coreg_initial_dy_m;
    double dz_m = opt.nuth_kaab_coreg_apply_vertical_shift
        ? opt.nuth_kaab_coreg_initial_dz_m
        : 0.0;

    for (int iter = 1; iter <= opt.nuth_kaab_coreg_max_iterations; ++iter)
    {
        std::vector<double> dh_step;
        std::vector<double> slope_step;
        std::vector<double> aspect_step;
        dh_step.reserve(candidate_indices.size());
        slope_step.reserve(candidate_indices.size());
        aspect_step.reserve(candidate_indices.size());

        for (std::size_t k = 0; k < candidate_indices.size(); ++k)
        {
            const std::size_t idx = candidate_indices[k];
            const int row = static_cast<int>(idx / static_cast<std::size_t>(ref_meta.n_cols));
            const int col = static_cast<int>(idx % static_cast<std::size_t>(ref_meta.n_cols));

            const lsdtt_xy xy = cell_center_from_row_col(ref_meta, row, col);
            const sampled_value smp = bilinear_sample_raster(moving, xy.x + dx_m, xy.y + dy_m);
            if (!smp.valid) continue;

            const double dh = static_cast<double>(ref_data[idx]) -
                              static_cast<double>(smp.value);
            if (!std::isfinite(dh)) continue;

            dh_step.push_back(dh);
            slope_step.push_back(slope_tan[idx]);
            aspect_step.push_back(aspect[idx]);
        }

        if (static_cast<int>(dh_step.size()) < opt.nuth_kaab_coreg_min_samples)
        {
            throw std::runtime_error(
                "Nuth-Kaab coreg: the subsample contains no more valid values "
                "(or too few) at this iteration");
        }

        const double vshift = compute_median(dh_step);
        std::vector<double> y;
        std::vector<double> x_aspect;
        y.reserve(dh_step.size());
        x_aspect.reserve(dh_step.size());

        for (std::size_t i = 0; i < dh_step.size(); ++i)
        {
            const double st = slope_step[i];
            if (!std::isfinite(st) || st <= 0.0) continue;
            const double yi = (dh_step[i] - vshift) / st;
            if (!std::isfinite(yi)) continue;
            y.push_back(yi);
            x_aspect.push_back(aspect_step[i]);
        }

        if (static_cast<int>(y.size()) < opt.nuth_kaab_coreg_min_samples)
        {
            throw std::runtime_error("Nuth-Kaab coreg: too few samples after slope-tangent filtering");
        }

        double easting_px = 0.0;
        double northing_px = 0.0;
        double fit_vertical = 0.0;
        fit_nuth_kaab_cosine_binned(
            x_aspect,
            y,
            opt.nuth_kaab_coreg_aspect_bin_count,
            easting_px,
            northing_px,
            fit_vertical
        );
        (void)fit_vertical;

        double ddx = easting_px * ref_meta.cell_size;
        double ddy = northing_px * ref_meta.cell_size;
        double ddz = opt.nuth_kaab_coreg_apply_vertical_shift ? vshift : 0.0;

        if (opt.nuth_kaab_coreg_max_horizontal_step_m > 0.0)
        {
            const double step_m = std::sqrt(ddx * ddx + ddy * ddy);
            const double max_step_m = opt.nuth_kaab_coreg_max_horizontal_step_m;
            if (step_m > max_step_m)
            {
                const double scale = max_step_m / step_m;
                ddx *= scale;
                ddy *= scale;
            }
        }

        dx_m += ddx;
        dy_m += ddy;
        if (opt.nuth_kaab_coreg_apply_vertical_shift)
        {
            dz_m = opt.nuth_kaab_coreg_initial_dz_m + ddz;
        }

        dem_coreg_iteration_diagnostics di;
        di.iteration = iter;
        di.n_samples = static_cast<int>(dh_step.size());
        di.median_residual = vshift;
        di.nmad_residual = 1.482602218505602 * compute_mad(dh_step);
        di.ddx_m = ddx;
        di.ddy_m = ddy;
        di.ddz_m = ddz;
        di.horizontal_step_pixels = std::sqrt(easting_px * easting_px + northing_px * northing_px);
        di.cumulative_dx_m = dx_m;
        di.cumulative_dy_m = dy_m;
        di.cumulative_dz_m = dz_m;
        result.diagnostics.push_back(di);

        result.iterations_run = iter;
        if (iter > 1 && di.horizontal_step_pixels <= opt.nuth_kaab_coreg_offset_threshold_pixels)
        {
            result.converged = true;
            break;
        }
    }

    result.dx_m = dx_m;
    result.dy_m = dy_m;
    result.dz_m = dz_m;
    result.corrected_moving = apply_shift_to_moving(moving, dx_m, dy_m, dz_m);
    return result;
}

} // anonymous namespace

dem_coreg_method parse_dem_coreg_method(const std::string& method_name)
{
    if (method_name == "nuth_kaab")
    {
        return dem_coreg_method::nuth_kaab;
    }
    throw std::invalid_argument(
        "Unknown dem_coregistration_method '" + method_name +
        "'. Allowed: nuth_kaab");
}

dem_coreg_result dem_coregister_to_reference(
    const lsdtt_raster_cpu& reference,
    const lsdtt_raster_cpu& moving,
    const dem_coreg_options& options)
{
    switch (options.method)
    {
        case dem_coreg_method::none:
        {
            dem_coreg_result out;
            out.corrected_moving = moving;
            out.converged = true;
            out.iterations_run = 0;
            return out;
        }
        case dem_coreg_method::nuth_kaab:
            return run_nuth_kaab(reference, moving, options.nuth_kaab);
    }

    throw std::invalid_argument("Unsupported dem_coreg_method enum value");
}

} // namespace lsdtt
