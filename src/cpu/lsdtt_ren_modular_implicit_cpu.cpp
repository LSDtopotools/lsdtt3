// =============================================================================
// lsdtt_ren_modular_implicit_cpu.cpp
// =============================================================================
#include "lsdtt_ren_modular_implicit_cpu.hpp"

#include <algorithm>
#include <cmath>

#ifdef LSDTT_HAS_EIGEN
#  include <Eigen/Sparse>
#  include <Eigen/SparseLU>
#endif

namespace lsdtt
{
namespace
{
bool valid_options(const lsdtt_ren_modular_options& o)
{
    return o.D > 0.0 && o.Sc > 0.0 && o.cell_size > 0.0 && o.dt > 0.0
        && o.archive_slope_margin >= 0.0 && o.archive_slope_margin < o.Sc
        && std::isfinite(o.D) && std::isfinite(o.U) && std::isfinite(o.Sc)
        && std::isfinite(o.cell_size) && std::isfinite(o.dt);
}
} // namespace

lsdtt_ren_modular_row evaluate_ren_modular_matrix_row(
    double phi, double phi_x, double phi_y,
    const lsdtt_ren_modular_options& o)
{
    lsdtt_ren_modular_row out;
    if (!valid_options(o) || !(phi > 0.0) || !std::isfinite(phi)
        || !std::isfinite(phi_x) || !std::isfinite(phi_y)) return out;
    const double h = o.cell_size;
    const double axial = phi * 4.0 / (6.0 * h * h);
    const double diagonal = phi / (6.0 * h * h);
    out.matrix[0] = 1.0 / o.dt + phi * 20.0 / (6.0 * h * h);
    out.matrix[1] = -axial + phi_x / (2.0 * h);
    out.matrix[2] = -axial - phi_x / (2.0 * h);
    out.matrix[3] = -axial + phi_y / (2.0 * h);
    out.matrix[4] = -axial - phi_y / (2.0 * h);
    out.matrix[5] = -diagonal;
    out.matrix[6] = -diagonal;
    out.matrix[7] = -diagonal;
    out.matrix[8] = -diagonal;
    out.valid = true;
    for (double value : out.matrix) out.valid = out.valid && std::isfinite(value);
    return out;
}

lsdtt_ren_modular_result compute_ren_modular_step_2d(
    const std::vector<double>& z_old, int n_rows, int n_cols,
    const std::vector<std::uint8_t>& is_fixed,
    const std::vector<double>& fixed_values, const lsdtt_ren_modular_options& o)
{
    lsdtt_ren_modular_result out;
    out.elevation = z_old;
#ifndef LSDTT_HAS_EIGEN
    (void)n_rows; (void)n_cols; (void)is_fixed; (void)fixed_values; (void)o;
    out.message = "Ren modular implicit requires Eigen; build with -DLSDTT_EIGEN=ON";
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
        out.message = "invalid Ren modular implicit inputs";
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

    const double h = o.cell_size;
    std::vector<double> phi(static_cast<std::size_t>(n), o.D);
    double max_ratio = 0.0;
    bool clamp_active = false;
    auto derivative_x = [&](const std::vector<double>& field, int r, int c) {
        if (c == 0) return (field[static_cast<std::size_t>(index(r, 1))]
                          - field[static_cast<std::size_t>(index(r, 0))]) / h;
        if (c == n_cols - 1) return (field[static_cast<std::size_t>(index(r, c))]
                                   - field[static_cast<std::size_t>(index(r, c - 1))]) / h;
        return (field[static_cast<std::size_t>(index(r, c + 1))]
              - field[static_cast<std::size_t>(index(r, c - 1))]) / (2.0 * h);
    };
    auto derivative_y = [&](const std::vector<double>& field, int r, int c) {
        if (r == 0) return (field[static_cast<std::size_t>(index(1, c))]
                          - field[static_cast<std::size_t>(index(0, c))]) / h;
        if (r == n_rows - 1) return (field[static_cast<std::size_t>(index(r, c))]
                                   - field[static_cast<std::size_t>(index(r - 1, c))]) / h;
        return (field[static_cast<std::size_t>(index(r + 1, c))]
              - field[static_cast<std::size_t>(index(r - 1, c))]) / (2.0 * h);
    };
    const double clamp_slope = o.Sc - o.archive_slope_margin;
    for (int r = 0; r < n_rows; ++r)
        for (int c = 0; c < n_cols; ++c)
        {
            const int i = index(r, c);
            const double gx = derivative_x(z_old, r, c);
            const double gy = derivative_y(z_old, r, c);
            const double raw_slope = std::sqrt(gx * gx + gy * gy);
            max_ratio = std::max(max_ratio, raw_slope / o.Sc);
            double slope = raw_slope;
            if (o.use_archive_slope_clamp && slope > clamp_slope)
            {
                slope = clamp_slope;
                clamp_active = true;
            }
            const double denominator = 1.0 - (slope / o.Sc) * (slope / o.Sc);
            if (!(denominator > 0.0) || !std::isfinite(denominator))
            {
                out.message = "Ren coefficient is non-finite or at/beyond Sc";
                return out;
            }
            phi[static_cast<std::size_t>(i)] = o.D / denominator;
        }

    std::vector<double> phi_x(static_cast<std::size_t>(n));
    std::vector<double> phi_y(static_cast<std::size_t>(n));
    for (int r = 0; r < n_rows; ++r)
        for (int c = 0; c < n_cols; ++c)
        {
            const int i = index(r, c);
            phi_x[static_cast<std::size_t>(i)] = derivative_x(phi, r, c);
            phi_y[static_cast<std::size_t>(i)] = derivative_y(phi, r, c);
        }

    using Triplet = Eigen::Triplet<double>;
    std::vector<Triplet> entries;
    entries.reserve(static_cast<std::size_t>(9 * n));
    Eigen::VectorXd rhs(n);
    static constexpr int dr[9] = {0, 0, 0, -1, 1, -1, -1, 1, 1};
    static constexpr int dc[9] = {0, -1, 1, 0, 0, -1, 1, -1, 1};
    for (int r = 0; r < n_rows; ++r)
        for (int c = 0; c < n_cols; ++c)
        {
            const int i = index(r, c);
            if (is_fixed[static_cast<std::size_t>(i)])
            {
                entries.emplace_back(i, i, 1.0);
                rhs[i] = fixed_values.empty() ? z_old[static_cast<std::size_t>(i)]
                                              : fixed_values[static_cast<std::size_t>(i)];
                continue;
            }
            const auto row = evaluate_ren_modular_matrix_row(
                phi[static_cast<std::size_t>(i)], phi_x[static_cast<std::size_t>(i)],
                phi_y[static_cast<std::size_t>(i)], o);
            if (!row.valid)
            {
                out.message = "invalid Ren modular matrix row";
                return out;
            }
            for (int k = 0; k < 9; ++k)
                entries.emplace_back(i, index(r + dr[k], c + dc[k]), row.matrix[k]);
            rhs[i] = z_old[static_cast<std::size_t>(i)] / o.dt + o.U;
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
    Eigen::VectorXd solution = solver.solve(rhs);
    if (solver.info() != Eigen::Success || !solution.allFinite())
    {
        out.message = "SparseLU solve failed";
        return out;
    }
    for (int i = 0; i < n; ++i)
        if (is_fixed[static_cast<std::size_t>(i)]) solution[i] = rhs[i];
    out.linear_residual_inf = (A * solution - rhs).lpNorm<Eigen::Infinity>();
    out.elevation.assign(solution.data(), solution.data() + n);
    out.slope_clamp_active = clamp_active;
    out.max_raw_slope_ratio = max_ratio;
    out.n_unknowns = n;
    out.completed = true;
    out.message = "completed";
    return out;
#endif
}

} // namespace lsdtt
