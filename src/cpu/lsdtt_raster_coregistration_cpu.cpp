// lsdtt_raster_coregistration_cpu.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html

#include "lsdtt_raster_coregistration_cpu.hpp"

#include "lsdtt_raster_surface_metrics_cpu.hpp"

#include "../utils/lsdtt_linear_algebra.hpp"
#include "../utils/lsdtt_raster_geometry.hpp"
#include "../utils/lsdtt_statistics.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
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
    if (opt.nuth_kaab_coreg_tolerance_horizontal_m <= 0.0)
    {
        throw std::invalid_argument("nuth_kaab_coreg_tolerance_horizontal_m must be > 0");
    }
    if (opt.nuth_kaab_coreg_slope_min_degrees < 0.0 ||
        opt.nuth_kaab_coreg_slope_min_degrees >= 90.0)
    {
        throw std::invalid_argument("nuth_kaab_coreg_slope_min_degrees must be in [0,90)");
    }
    if (opt.nuth_kaab_coreg_slope_max_degrees <= opt.nuth_kaab_coreg_slope_min_degrees ||
        opt.nuth_kaab_coreg_slope_max_degrees >= 90.0)
    {
        throw std::invalid_argument("nuth_kaab_coreg_slope_max_degrees must be in (slope_min,90)");
    }
    if (opt.nuth_kaab_coreg_residual_nmad_clip <= 0.0)
    {
        throw std::invalid_argument("nuth_kaab_coreg_residual_nmad_clip must be > 0");
    }
    if (opt.nuth_kaab_coreg_min_samples < 25)
    {
        throw std::invalid_argument("nuth_kaab_coreg_min_samples must be >= 25");
    }
    if (opt.nuth_kaab_coreg_subsample_step < 1)
    {
        throw std::invalid_argument("nuth_kaab_coreg_subsample_step must be >= 1");
    }
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

    polyfit_metric_request req;
    req.metrics = {
        polyfit_metric::dzdx,
        polyfit_metric::dzdy,
        polyfit_metric::slope_degrees
    };
    req.n_threads = opt.nuth_kaab_coreg_polyfit_threads;

    const poly_fit_coeffs poly = direct_polyfit(
        reference,
        2.0 * ref_meta.cell_size,
        opt.nuth_kaab_coreg_polyfit_threads
    );
    const polyfit_metric_result deriv = compute_polyfit_metrics(poly, req);

    const std::vector<float>& ref_data = reference.get_data();
    const std::vector<uint8_t>& ref_mask = reference.get_mask();

    const std::vector<float>& p_data = deriv.at(polyfit_metric::dzdx).get_data();
    const std::vector<uint8_t>& p_mask = deriv.at(polyfit_metric::dzdx).get_mask();

    const std::vector<float>& q_data = deriv.at(polyfit_metric::dzdy).get_data();
    const std::vector<uint8_t>& q_mask = deriv.at(polyfit_metric::dzdy).get_mask();

    const std::vector<float>& slope_data = deriv.at(polyfit_metric::slope_degrees).get_data();
    const std::vector<uint8_t>& slope_mask = deriv.at(polyfit_metric::slope_degrees).get_mask();

    dem_coreg_result result;

    double dx_m = 0.0;
    double dy_m = 0.0;
    double dz_m = 0.0;

    const int subsample_step = opt.nuth_kaab_coreg_subsample_step;
    const float slope_min = static_cast<float>(opt.nuth_kaab_coreg_slope_min_degrees);
    const float slope_max = static_cast<float>(opt.nuth_kaab_coreg_slope_max_degrees);

    for (int iter = 1; iter <= opt.nuth_kaab_coreg_max_iterations; ++iter)
    {
        std::vector<double> a_dx;
        std::vector<double> a_dy;
        std::vector<double> a_dz;
        std::vector<double> residuals;
        a_dx.reserve(ref_data.size() / 8);
        a_dy.reserve(ref_data.size() / 8);
        a_dz.reserve(ref_data.size() / 8);
        residuals.reserve(ref_data.size() / 8);

        for (int row = 0; row < ref_meta.n_rows; row += subsample_step)
        {
            for (int col = 0; col < ref_meta.n_cols; col += subsample_step)
            {
                const std::size_t idx = static_cast<std::size_t>(row) *
                                        static_cast<std::size_t>(ref_meta.n_cols) +
                                        static_cast<std::size_t>(col);

                if (ref_mask[idx] == 0u || p_mask[idx] == 0u || q_mask[idx] == 0u ||
                    slope_mask[idx] == 0u)
                {
                    continue;
                }

                const float slope_deg = slope_data[idx];
                if (!std::isfinite(slope_deg) || slope_deg < slope_min || slope_deg > slope_max)
                {
                    continue;
                }

                const float ref_z = ref_data[idx];
                const float p = p_data[idx];
                const float q = q_data[idx];
                if (!std::isfinite(ref_z) || !std::isfinite(p) || !std::isfinite(q))
                {
                    continue;
                }

                const lsdtt_xy xy = cell_center_from_row_col(ref_meta, row, col);
                const sampled_value smp = bilinear_sample_raster(moving, xy.x + dx_m, xy.y + dy_m);
                if (!smp.valid) continue;

                const double residual = static_cast<double>(ref_z) -
                                        (static_cast<double>(smp.value) + dz_m);

                // First-order Nuth-Kaab linearisation in x/y with local slopes:
                // residual ~= -(ddx * dz/dx + ddy * dz/dy + ddz)
                a_dx.push_back(-static_cast<double>(p));
                a_dy.push_back(-static_cast<double>(q));
                a_dz.push_back(opt.nuth_kaab_coreg_apply_vertical_shift ? -1.0 : 0.0);
                residuals.push_back(residual);
            }
        }

        if (static_cast<int>(residuals.size()) < opt.nuth_kaab_coreg_min_samples)
        {
            throw std::runtime_error("Nuth-Kaab coreg: too few samples after slope filtering");
        }

        const double median_resid = compute_median(residuals);
        const double nmad = 1.482602218505602 * compute_mad(residuals);

        std::vector<double> f_a_dx;
        std::vector<double> f_a_dy;
        std::vector<double> f_a_dz;
        std::vector<double> f_resid;
        f_a_dx.reserve(residuals.size());
        f_a_dy.reserve(residuals.size());
        f_a_dz.reserve(residuals.size());
        f_resid.reserve(residuals.size());

        const double clip = opt.nuth_kaab_coreg_residual_nmad_clip;
        const double tol = clip * (std::isfinite(nmad) ? nmad : 0.0);
        const bool use_clip = std::isfinite(tol) && tol > 0.0;

        for (std::size_t i = 0; i < residuals.size(); ++i)
        {
            const double r = residuals[i];
            if (use_clip && std::abs(r - median_resid) > tol) continue;
            f_a_dx.push_back(a_dx[i]);
            f_a_dy.push_back(a_dy[i]);
            f_a_dz.push_back(a_dz[i]);
            f_resid.push_back(r);
        }

        if (static_cast<int>(f_resid.size()) < opt.nuth_kaab_coreg_min_samples)
        {
            throw std::runtime_error("Nuth-Kaab coreg: too few samples after robust clipping");
        }

        double ddx = 0.0;
        double ddy = 0.0;
        double ddz = 0.0;
        solve_weighted_3param_lsq(f_a_dx, f_a_dy, f_a_dz, f_resid, ddx, ddy, ddz);

        if (!opt.nuth_kaab_coreg_apply_vertical_shift)
        {
            ddz = 0.0;
        }

        if (opt.nuth_kaab_coreg_max_horizontal_step_m > 0.0)
        {
            const double step = std::sqrt(ddx * ddx + ddy * ddy);
            const double max_step = opt.nuth_kaab_coreg_max_horizontal_step_m;
            if (step > max_step)
            {
                const double scale = max_step / step;
                ddx *= scale;
                ddy *= scale;
            }
        }

        dx_m += ddx;
        dy_m += ddy;
        dz_m += ddz;

        dem_coreg_iteration_diagnostics di;
        di.iteration = iter;
        di.n_samples = static_cast<int>(f_resid.size());
        di.median_residual = median_resid;
        di.nmad_residual = nmad;
        di.ddx_m = ddx;
        di.ddy_m = ddy;
        di.ddz_m = ddz;
        di.cumulative_dx_m = dx_m;
        di.cumulative_dy_m = dy_m;
        di.cumulative_dz_m = dz_m;
        result.diagnostics.push_back(di);

        result.iterations_run = iter;
        if (std::sqrt(ddx * ddx + ddy * ddy) <= opt.nuth_kaab_coreg_tolerance_horizontal_m)
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
    if (method_name.empty() || method_name == "none")
    {
        return dem_coreg_method::none;
    }
    if (method_name == "nuth_kaab")
    {
        return dem_coreg_method::nuth_kaab;
    }
    throw std::invalid_argument(
        "Unknown dem_coregistration_method '" + method_name +
        "'. Allowed: none, nuth_kaab");
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
