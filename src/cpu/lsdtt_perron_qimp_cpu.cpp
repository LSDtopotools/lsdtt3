// lsdtt_perron_qimp_cpu.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html

// =============================================================================
// lsdtt_perron_qimp_cpu.cpp
// =============================================================================
#include "lsdtt_perron_qimp_cpu.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

#ifdef LSDTT_HAS_EIGEN
#  include <Eigen/Sparse>
#  include <Eigen/SparseLU>
#endif

namespace lsdtt
{
namespace
{

bool valid_options(const lsdtt_perron_qimp_options& o)
{
    return o.K > 0.0 && o.Sc > 0.0 && o.dx > 0.0 && o.dy > 0.0 && o.dt > 0.0
        && std::isfinite(o.K) && std::isfinite(o.U) && std::isfinite(o.Sc)
        && std::isfinite(o.dx) && std::isfinite(o.dy) && std::isfinite(o.dt);
}

} // namespace

lsdtt_perron_qimp_row_1d evaluate_perron_qimp_row_1d(
    double zm, double zc, double zp, const lsdtt_perron_qimp_options& o)
{
    lsdtt_perron_qimp_row_1d out;
    if (!valid_options(o)) return out;
    const double p = (zp - zm) / (2.0 * o.dx);
    const double r = (zp - 2.0 * zc + zm) / (o.dx * o.dx);
    const double b = 1.0 / (o.Sc * o.Sc);
    const double bp2 = b * p * p;
    const double denominator = 1.0 - bp2;
    out.slope_ratio = std::fabs(p) / o.Sc;
    if (!(denominator > 0.0) || !std::isfinite(denominator)) return out;

    const double R = o.K * (1.0 + bp2) / (denominator * denominator);
    const double P = o.K * r * 2.0 * b * p * (3.0 + bp2)
                   / (denominator * denominator * denominator);
    out.rate = o.U + R * r;
    out.j_minus = R / (o.dx * o.dx) - P / (2.0 * o.dx);
    out.j_center = -2.0 * R / (o.dx * o.dx);
    out.j_plus = R / (o.dx * o.dx) + P / (2.0 * o.dx);
    out.valid = std::isfinite(out.rate) && std::isfinite(out.j_minus)
             && std::isfinite(out.j_center) && std::isfinite(out.j_plus);
    return out;
}

lsdtt_perron_qimp_row_2d evaluate_perron_qimp_row_2d(
    double zc, double zxm, double zxp, double zym, double zyp,
    double zxmym, double zxpym, double zxmyp, double zxpyp,
    const lsdtt_perron_qimp_options& o)
{
    lsdtt_perron_qimp_row_2d out;
    if (!valid_options(o)) return out;
    const double p = (zxp - zxm) / (2.0 * o.dx);
    const double q = (zyp - zym) / (2.0 * o.dy);
    const double r = (zxp - 2.0 * zc + zxm) / (o.dx * o.dx);
    const double s = (zyp - 2.0 * zc + zym) / (o.dy * o.dy);
    const double t = (zxpyp - zxmyp - zxpym + zxmym) / (4.0 * o.dx * o.dy);
    const double b = 1.0 / (o.Sc * o.Sc);
    const double denominator = 1.0 - b * (p * p + q * q);
    out.slope_ratio = std::sqrt(p * p + q * q) / o.Sc;
    if (!(denominator > 0.0) || !std::isfinite(denominator)) return out;
    const double H = 1.0 / denominator;
    const double L = r + s;
    const double M = p * p * r + q * q * s + 2.0 * p * q * t;
    const double P = o.K * (2.0 * b * p * L * H * H
        + 4.0 * b * (p * r + q * t) * H * H
        + 8.0 * b * b * p * M * H * H * H);
    const double Q = o.K * (2.0 * b * q * L * H * H
        + 4.0 * b * (q * s + p * t) * H * H
        + 8.0 * b * b * q * M * H * H * H);
    const double R = o.K * (H + 2.0 * b * p * p * H * H);
    const double S = o.K * (H + 2.0 * b * q * q * H * H);
    const double T = 4.0 * o.K * b * p * q * H * H;

    out.rate = o.U + o.K * (L * H + 2.0 * b * M * H * H);
    out.jacobian[0] = -2.0 * R / (o.dx * o.dx) - 2.0 * S / (o.dy * o.dy);
    out.jacobian[1] = -P / (2.0 * o.dx) + R / (o.dx * o.dx);
    out.jacobian[2] =  P / (2.0 * o.dx) + R / (o.dx * o.dx);
    out.jacobian[3] = -Q / (2.0 * o.dy) + S / (o.dy * o.dy);
    out.jacobian[4] =  Q / (2.0 * o.dy) + S / (o.dy * o.dy);
    out.jacobian[5] =  T / (4.0 * o.dx * o.dy);
    out.jacobian[6] = -T / (4.0 * o.dx * o.dy);
    out.jacobian[7] = -T / (4.0 * o.dx * o.dy);
    out.jacobian[8] =  T / (4.0 * o.dx * o.dy);
    out.valid = std::isfinite(out.rate);
    for (double value : out.jacobian) out.valid = out.valid && std::isfinite(value);
    return out;
}

lsdtt_perron_qimp_result compute_perron_qimp_step_1d(
    const std::vector<double>& z_old, const std::vector<std::uint8_t>& is_fixed,
    const std::vector<double>& fixed_values, const lsdtt_perron_qimp_options& o)
{
    lsdtt_perron_qimp_result out;
    out.elevation = z_old;
#ifndef LSDTT_HAS_EIGEN
    (void)is_fixed; (void)fixed_values; (void)o;
    out.message = "Perron Q-imp requires Eigen; build with -DLSDTT_EIGEN=ON";
    return out;
#else
    out.had_eigen = true;
    out.solver_used = "Eigen::SparseLU";
    const int n = static_cast<int>(z_old.size());
    if (n < 3 || is_fixed.size() != z_old.size()
        || (!fixed_values.empty() && fixed_values.size() != z_old.size())
        || !valid_options(o) || !is_fixed.front() || !is_fixed.back())
    {
        out.message = "invalid 1-D Q-imp inputs; endpoints must be fixed";
        return out;
    }
    using Triplet = Eigen::Triplet<double>;
    std::vector<Triplet> entries;
    entries.reserve(static_cast<std::size_t>(3 * n));
    Eigen::VectorXd rhs(n);
    double max_ratio = 0.0;
    for (int i = 0; i < n; ++i)
    {
        if (is_fixed[static_cast<std::size_t>(i)])
        {
            entries.emplace_back(i, i, 1.0);
            const double target = fixed_values.empty() ? z_old[static_cast<std::size_t>(i)]
                                                       : fixed_values[static_cast<std::size_t>(i)];
            rhs[i] = target - z_old[static_cast<std::size_t>(i)];
            continue;
        }
        const auto row = evaluate_perron_qimp_row_1d(
            z_old[static_cast<std::size_t>(i - 1)], z_old[static_cast<std::size_t>(i)],
            z_old[static_cast<std::size_t>(i + 1)], o);
        if (!row.valid)
        {
            out.message = "Q-imp coefficient is non-finite or at/beyond Sc";
            return out;
        }
        max_ratio = std::max(max_ratio, row.slope_ratio);
        entries.emplace_back(i, i - 1, -o.dt * row.j_minus);
        entries.emplace_back(i, i, 1.0 - o.dt * row.j_center);
        entries.emplace_back(i, i + 1, -o.dt * row.j_plus);
        rhs[i] = o.dt * row.rate;
    }
    Eigen::SparseMatrix<double> A(n, n);
    A.setFromTriplets(entries.begin(), entries.end());
    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
    solver.compute(A);
    if (solver.info() != Eigen::Success)
    {
        out.message = "SparseLU factorization failed";
        return out;
    }
    Eigen::VectorXd increment = solver.solve(rhs);
    if (solver.info() != Eigen::Success || !increment.allFinite())
    {
        out.message = "SparseLU solve failed";
        return out;
    }
    for (int i = 0; i < n; ++i)
        if (is_fixed[static_cast<std::size_t>(i)]) increment[i] = rhs[i];
    out.linear_residual_inf = (A * increment - rhs).lpNorm<Eigen::Infinity>();
    for (int i = 0; i < n; ++i) out.elevation[static_cast<std::size_t>(i)] += increment[i];
    out.max_slope_ratio = max_ratio;
    out.n_unknowns = n;
    out.completed = true;
    out.message = "completed";
    return out;
#endif
}

lsdtt_perron_qimp_result compute_perron_qimp_step_2d(
    const std::vector<double>& z_old, int n_rows, int n_cols,
    const std::vector<std::uint8_t>& is_fixed,
    const std::vector<double>& fixed_values, const lsdtt_perron_qimp_options& o)
{
    lsdtt_perron_qimp_result out;
    out.elevation = z_old;
#ifndef LSDTT_HAS_EIGEN
    (void)n_rows; (void)n_cols; (void)is_fixed; (void)fixed_values; (void)o;
    out.message = "Perron Q-imp requires Eigen; build with -DLSDTT_EIGEN=ON";
    return out;
#else
    out.had_eigen = true;
    out.solver_used = "Eigen::SparseLU";
    const int n = n_rows * n_cols;
    if (n_rows < 3 || n_cols < 3 || static_cast<int>(z_old.size()) != n
        || is_fixed.size() != z_old.size()
        || (!fixed_values.empty() && fixed_values.size() != z_old.size())
        || !valid_options(o))
    {
        out.message = "invalid 2-D Q-imp inputs";
        return out;
    }
    auto index = [n_cols](int r, int c) { return r * n_cols + c; };
    for (int r = 0; r < n_rows; ++r)
        for (int c = 0; c < n_cols; ++c)
            if ((r == 0 || c == 0 || r == n_rows - 1 || c == n_cols - 1)
                && !is_fixed[static_cast<std::size_t>(index(r, c))])
            {
                out.message = "all rectangular edge cells must be fixed";
                return out;
            }

    using Triplet = Eigen::Triplet<double>;
    std::vector<Triplet> entries;
    entries.reserve(static_cast<std::size_t>(9 * n));
    Eigen::VectorXd rhs(n);
    double max_ratio = 0.0;
    static constexpr int dr[9] = {0, 0, 0, -1, 1, -1, -1, 1, 1};
    static constexpr int dc[9] = {0, -1, 1, 0, 0, -1, 1, -1, 1};
    for (int r = 0; r < n_rows; ++r)
        for (int c = 0; c < n_cols; ++c)
        {
            const int i = index(r, c);
            if (is_fixed[static_cast<std::size_t>(i)])
            {
                entries.emplace_back(i, i, 1.0);
                const double target = fixed_values.empty() ? z_old[static_cast<std::size_t>(i)]
                                                           : fixed_values[static_cast<std::size_t>(i)];
                rhs[i] = target - z_old[static_cast<std::size_t>(i)];
                continue;
            }
            const auto row = evaluate_perron_qimp_row_2d(
                z_old[static_cast<std::size_t>(i)],
                z_old[static_cast<std::size_t>(index(r, c - 1))],
                z_old[static_cast<std::size_t>(index(r, c + 1))],
                z_old[static_cast<std::size_t>(index(r - 1, c))],
                z_old[static_cast<std::size_t>(index(r + 1, c))],
                z_old[static_cast<std::size_t>(index(r - 1, c - 1))],
                z_old[static_cast<std::size_t>(index(r - 1, c + 1))],
                z_old[static_cast<std::size_t>(index(r + 1, c - 1))],
                z_old[static_cast<std::size_t>(index(r + 1, c + 1))], o);
            if (!row.valid)
            {
                out.message = "Q-imp coefficient is non-finite or at/beyond Sc";
                return out;
            }
            max_ratio = std::max(max_ratio, row.slope_ratio);
            for (int k = 0; k < 9; ++k)
            {
                const int j = index(r + dr[k], c + dc[k]);
                entries.emplace_back(i, j, (k == 0 ? 1.0 : 0.0) - o.dt * row.jacobian[k]);
            }
            rhs[i] = o.dt * row.rate;
        }
    Eigen::SparseMatrix<double> A(n, n);
    A.setFromTriplets(entries.begin(), entries.end());
    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
    solver.compute(A);
    if (solver.info() != Eigen::Success)
    {
        out.message = "SparseLU factorization failed";
        return out;
    }
    Eigen::VectorXd increment = solver.solve(rhs);
    if (solver.info() != Eigen::Success || !increment.allFinite())
    {
        out.message = "SparseLU solve failed";
        return out;
    }
    for (int i = 0; i < n; ++i)
        if (is_fixed[static_cast<std::size_t>(i)]) increment[i] = rhs[i];
    out.linear_residual_inf = (A * increment - rhs).lpNorm<Eigen::Infinity>();
    for (int i = 0; i < n; ++i) out.elevation[static_cast<std::size_t>(i)] += increment[i];
    out.max_slope_ratio = max_ratio;
    out.n_unknowns = n;
    out.completed = true;
    out.message = "completed";
    return out;
#endif
}

} // namespace lsdtt
