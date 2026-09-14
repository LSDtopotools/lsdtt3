// lsdtt_chi_ksn_inversion_cpu.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// lsdtt_chi_ksn_inversion_cpu
//
// Spatial k_sn / u* inversion with m/n model selection.
//
// This file implements a first version of the spatial chi-based inversion:
// for each tested m/n value, chi is recomputed on a fixed channel network,
// a sparse design matrix A is assembled from downstream chi increments,
// and a coarse spatial field of u* is solved by regularized least squares.
//
// Physical interpretation:
// ------------------------
// For a fixed trial m/n, the chi coordinate is fixed. In chi-z space the local
// slope is u* (equivalent to k_sn when A_0 = 1 m^2). Elevation at a node is
// therefore the sum of downstream delta-chi contributions weighted by the u*
// values assigned to the coarse inversion cells crossed by the downstream path.
//
// Numerically, for each node i we form one row of a sparse matrix A such that:
//
//   z_rel = A * u
//
// where:
//   z_rel : node elevations relative to the outlet elevation
//   u     : one u* parameter per active coarse inversion cell
//
// Regularization:
// ---------------
// A discrete 4-neighbour Laplacian L is built on the active coarse cells.
// The version-1 inversion solves the regularized unconstrained system:
//
//   (A^TA + lambda L^TL + epsilon I) u = A^T z_rel
//
// where:
//   lambda : roughness_weight
//   epsilon : diagonal_jitter
//
// Matlab relation:
// ----------------
// This reproduces the core spatial k_sn / u* inversion logic discussed in the
// paper and in the Matlab code, but in a simplified first version:
//
//   - one spatial field of unknowns
//   - no box constraints yet
//   - no second parameter block
//   - no temporal uplift-history extension
//
// -----------------------------------------------------------------------------

#include "lsdtt_chi_ksn_inversion_cpu.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <fstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>
#include <Eigen/SparseQR>

namespace lsdtt
{
namespace
{

constexpr double k_sqrt_2 = 1.4142135623730950488;

// -----------------------------------------------------------------------------
// Helpers
// -----------------------------------------------------------------------------
float quiet_nan_f()
{
    return std::numeric_limits<float>::quiet_NaN();
}

double quiet_nan_d()
{
    return std::numeric_limits<double>::quiet_NaN();
}

bool is_finite(float x)
{
    return std::isfinite(static_cast<double>(x));
}

bool is_finite(double x)
{
    return std::isfinite(x);
}

bool nearly_equal(double a, double b)
{
    const double scale =
        std::max(1.0, std::max(std::fabs(a), std::fabs(b)));

    return std::fabs(a - b) <= 1e-9 * scale;
}

void validate_options(
    const lsdtt_chi_ksn_inversion_options& options,
    const std::vector<double>& tested_m_over_n)
{
    if (!is_finite(options.A_0) || options.A_0 <= 0.0)
    {
        throw std::invalid_argument(
            "lsdtt_chi_ksn_inversion_cpu: A_0 must be finite and > 0."
        );
    }

    if (!is_finite(options.coarse_cell_size) || options.coarse_cell_size <= 0.0)
    {
        throw std::invalid_argument(
            "lsdtt_chi_ksn_inversion_cpu: coarse_cell_size must be finite and > 0."
        );
    }

    if (!is_finite(options.roughness_weight) || options.roughness_weight < 0.0)
    {
        throw std::invalid_argument(
            "lsdtt_chi_ksn_inversion_cpu: roughness_weight must be finite and >= 0."
        );
    }

    if (!is_finite(options.diagonal_jitter) || options.diagonal_jitter < 0.0)
    {
        throw std::invalid_argument(
            "lsdtt_chi_ksn_inversion_cpu: diagonal_jitter must be finite and >= 0."
        );
    }

    if (tested_m_over_n.empty())
    {
        throw std::invalid_argument(
            "lsdtt_chi_ksn_inversion_cpu: tested_m_over_n is empty."
        );
    }

    for (const double mn : tested_m_over_n)
    {
        if (!is_finite(mn))
        {
            throw std::invalid_argument(
                "lsdtt_chi_ksn_inversion_cpu: tested_m_over_n contains a non-finite value."
            );
        }
    }
}

void validate_flow_for_inversion(const lsdtt_flow_base& flow)
{
    if (flow.empty())
    {
        throw std::runtime_error(
            "lsdtt_chi_ksn_inversion_cpu: flow object is empty."
        );
    }

    if (!flow.has_s_vector())
    {
        throw std::runtime_error(
            "lsdtt_chi_ksn_inversion_cpu: inversion requires a flow object "
            "with s_vector built."
        );
    }

    flow.meta().check_projected_metres();
}

void assert_raster_matches_flow_grid(
    const lsdtt_flow_base& flow,
    const lsdtt_raster_cpu& raster,
    const std::string& raster_name)
{
    const lsdtt_raster_meta& fm = flow.meta();
    const lsdtt_raster_meta& rm = raster.get_meta();

    if (fm.n_rows != rm.n_rows || fm.n_cols != rm.n_cols)
    {
        throw std::invalid_argument(
            "lsdtt_chi_ksn_inversion_cpu: " + raster_name +
            " dimensions do not match the flow grid."
        );
    }

    if (!nearly_equal(fm.cell_size, rm.cell_size) ||
        !nearly_equal(fm.x_min,     rm.x_min)     ||
        !nearly_equal(fm.y_min,     rm.y_min))
    {
        throw std::invalid_argument(
            "lsdtt_chi_ksn_inversion_cpu: " + raster_name +
            " georeferencing does not match the flow grid."
        );
    }
}

void validate_inversion_network(
    const lsdtt_flow_base& flow,
    const lsdtt_chi_ksn_inversion_network& network)
{
    const std::size_t n = network.flats.size();

    if (n == 0)
    {
        throw std::invalid_argument(
            "lsdtt_chi_ksn_inversion_cpu: inversion network is empty."
        );
    }

    if (network.receiver_flats.size() != n ||
        network.observed_elevation.size() != n)
    {
        throw std::invalid_argument(
            "lsdtt_chi_ksn_inversion_cpu: inversion network vectors have inconsistent lengths."
        );
    }

    if (network.outlet_flat == k_nodata_receiver)
    {
        throw std::invalid_argument(
            "lsdtt_chi_ksn_inversion_cpu: inversion network has invalid outlet_flat."
        );
    }

    if (network.flats[0] != network.outlet_flat)
    {
        throw std::invalid_argument(
            "lsdtt_chi_ksn_inversion_cpu: inversion network must begin with the outlet node."
        );
    }

    if (network.receiver_flats[0] != k_outlet_receiver)
    {
        throw std::invalid_argument(
            "lsdtt_chi_ksn_inversion_cpu: outlet node must have receiver_flats[0] == k_outlet_receiver."
        );
    }

    std::unordered_map<std::int64_t, int> local_index;
    local_index.reserve(n);

    for (std::size_t i = 0; i < n; ++i)
    {
        const std::int64_t flat = network.flats[i];

        if (!flow.is_valid_cell(flat))
        {
            throw std::invalid_argument(
                "lsdtt_chi_ksn_inversion_cpu: inversion network contains an invalid flat index."
            );
        }

        if (!is_finite(network.observed_elevation[i]))
        {
            throw std::invalid_argument(
                "lsdtt_chi_ksn_inversion_cpu: inversion network contains non-finite observed elevations."
            );
        }

        const auto inserted = local_index.emplace(flat, static_cast<int>(i));
        if (!inserted.second)
        {
            throw std::invalid_argument(
                "lsdtt_chi_ksn_inversion_cpu: inversion network contains duplicate flat indices."
            );
        }
    }

    for (std::size_t i = 1; i < n; ++i)
    {
        const std::int64_t flat = network.flats[i];
        const std::int64_t rec  = network.receiver_flats[i];

        if (rec == k_outlet_receiver)
        {
            throw std::invalid_argument(
                "lsdtt_chi_ksn_inversion_cpu: non-root node routes to k_outlet_receiver."
            );
        }

        const auto found = local_index.find(rec);
        if (found == local_index.end())
        {
            throw std::invalid_argument(
                "lsdtt_chi_ksn_inversion_cpu: receiver flat is not present in the inversion network."
            );
        }

        if (found->second >= static_cast<int>(i))
        {
            throw std::invalid_argument(
                "lsdtt_chi_ksn_inversion_cpu: inversion network is not in receiver-before-donor order."
            );
        }

        if (flow.receiver_of_unchecked(flat) != rec)
        {
            throw std::invalid_argument(
                "lsdtt_chi_ksn_inversion_cpu: inversion network receiver does not match the flow object."
            );
        }
    }
}

double step_length_from_receiver_dir(std::uint8_t dir, double cell_size)
{
    if (dir == lsdtt_flow_base::no_direction)
    {
        throw std::runtime_error(
            "lsdtt_chi_ksn_inversion_cpu: encountered no_direction on a non-root node."
        );
    }

    return (dir % 2u == 0u) ? cell_size : cell_size * k_sqrt_2;
}

float sample_dem_at_flat(
    const lsdtt_raster_cpu& dem,
    std::int64_t flat)
{
    const std::vector<float>&   data = dem.get_data();
    const std::vector<uint8_t>& mask = dem.get_mask();

    if (flat < 0 || flat >= static_cast<std::int64_t>(data.size()))
    {
        throw std::invalid_argument(
            "lsdtt_chi_ksn_inversion_cpu: flat index lies outside DEM bounds."
        );
    }

    const std::size_t idx = static_cast<std::size_t>(flat);

    if (mask[idx] == 0)
    {
        throw std::runtime_error(
            "lsdtt_chi_ksn_inversion_cpu: DEM is nodata at an inversion-network node."
        );
    }

    const float z = data[idx];
    if (!is_finite(z))
    {
        throw std::runtime_error(
            "lsdtt_chi_ksn_inversion_cpu: DEM contains non-finite elevation at an inversion-network node."
        );
    }

    return z;
}

// -----------------------------------------------------------------------------
// Internal inverter
// -----------------------------------------------------------------------------
class chi_ksn_inverter
{
public:
    chi_ksn_inverter(
        const lsdtt_flow_base& flow,
        const lsdtt_chi_ksn_inversion_network& network,
        const std::vector<double>& tested_m_over_n,
        const lsdtt_chi_ksn_inversion_options& options
    )
        : flow_(flow)
        , network_(network)
        , tested_m_over_n_(tested_m_over_n)
        , options_(options)
    {
    }

    lsdtt_chi_ksn_inversion_result run_area_based()
    {
        validate_flow_for_inversion(flow_);
        validate_inversion_network(flow_, network_);
        validate_options(options_, tested_m_over_n_);

        build_static_graph();
        build_coarse_grid_mapping();

        auto denominator_fn = [&](std::int64_t flat) -> double
        {
            const double pixel_area = flow_.cell_size() * flow_.cell_size();
            return static_cast<double>(flow_.s_vector_size_unchecked(flat)) * pixel_area;
        };

        return run_with_denominator(denominator_fn);
    }

    lsdtt_chi_ksn_inversion_result run_custom_denominator(
        const lsdtt_raster_cpu& denominator_raster)
    {
        validate_flow_for_inversion(flow_);
        validate_inversion_network(flow_, network_);
        validate_options(options_, tested_m_over_n_);
        assert_raster_matches_flow_grid(flow_, denominator_raster, "denominator_raster");

        build_static_graph();
        build_coarse_grid_mapping();

        const std::vector<float>&   data = denominator_raster.get_data();
        const std::vector<uint8_t>& mask = denominator_raster.get_mask();

        auto denominator_fn = [&](std::int64_t flat) -> double
        {
            const std::size_t idx = static_cast<std::size_t>(flat);
            if (mask[idx] == 0)
            {
                throw std::runtime_error(
                    "lsdtt_chi_ksn_inversion_cpu: denominator_raster is nodata at an inversion-network node."
                );
            }

            const double denom = static_cast<double>(data[idx]);
            if (!is_finite(denom) || denom <= 0.0)
            {
                throw std::runtime_error(
                    "lsdtt_chi_ksn_inversion_cpu: denominator_raster must be finite and > 0."
                );
            }

            return denom;
        };

        return run_with_denominator(denominator_fn);
    }

private:
    const lsdtt_flow_base& flow_;
    const lsdtt_chi_ksn_inversion_network& network_;
    const std::vector<double>& tested_m_over_n_;
    const lsdtt_chi_ksn_inversion_options& options_;

    // Static graph / indexing
    std::unordered_map<std::int64_t, int> local_index_of_flat_;
    std::vector<int> receiver_local_;
    std::vector<int> row_of_node_;
    std::vector<int> col_of_node_;
    std::vector<double> step_to_receiver_;

    // Coarse grid
    lsdtt_raster_meta coarse_meta_{};
    int coarse_n_rows_ = 0;
    int coarse_n_cols_ = 0;

    struct coarse_rc
    {
        int row = -1;
        int col = -1;
    };

    std::vector<coarse_rc> param_rc_of_id_;
    std::vector<int> coarse_cell_flat_of_id_;
    std::vector<std::int32_t> coarse_cell_id_of_node_;
    std::vector<int> support_count_of_param_;

private:
    void build_static_graph()
    {
        const std::size_t n = network_.flats.size();

        local_index_of_flat_.clear();
        local_index_of_flat_.reserve(n);

        receiver_local_.assign(n, -1);
        row_of_node_.assign(n, 0);
        col_of_node_.assign(n, 0);
        step_to_receiver_.assign(n, 0.0);

        const int n_cols = flow_.meta().n_cols;
        const double cell_size = flow_.cell_size();

        for (std::size_t i = 0; i < n; ++i)
        {
            local_index_of_flat_.emplace(network_.flats[i], static_cast<int>(i));

            const std::int64_t flat = network_.flats[i];
            row_of_node_[i] = static_cast<int>(flat / static_cast<std::int64_t>(n_cols));
            col_of_node_[i] = static_cast<int>(flat % static_cast<std::int64_t>(n_cols));
        }

        for (std::size_t i = 1; i < n; ++i)
        {
            const std::int64_t flat = network_.flats[i];
            const std::int64_t rec  = network_.receiver_flats[i];

            receiver_local_[i] = local_index_of_flat_.at(rec);
            step_to_receiver_[i] = step_length_from_receiver_dir(
                flow_.receiver_dir_of_unchecked(flat),
                cell_size
            );
        }
    }

    void build_coarse_grid_mapping()
    {
        const lsdtt_raster_meta& fm = flow_.meta();

        const double width  = static_cast<double>(fm.n_cols) * fm.cell_size;
        const double height = static_cast<double>(fm.n_rows) * fm.cell_size;

        coarse_n_cols_ = std::max(
            1,
            static_cast<int>(std::ceil(width / options_.coarse_cell_size))
        );
        coarse_n_rows_ = std::max(
            1,
            static_cast<int>(std::ceil(height / options_.coarse_cell_size))
        );

        coarse_meta_ = fm;
        coarse_meta_.n_cols = coarse_n_cols_;
        coarse_meta_.n_rows = coarse_n_rows_;
        coarse_meta_.cell_size = options_.coarse_cell_size;
        coarse_meta_.x_min = fm.x_min;
        coarse_meta_.y_min = fm.y_min;

        std::unordered_map<int, int> param_id_of_coarse_flat;
        param_id_of_coarse_flat.reserve(network_.flats.size());

        coarse_cell_id_of_node_.assign(network_.flats.size(), -1);
        param_rc_of_id_.clear();
        coarse_cell_flat_of_id_.clear();
        support_count_of_param_.clear();

        for (std::size_t i = 0; i < network_.flats.size(); ++i)
        {
            const int row = row_of_node_[i];
            const int col = col_of_node_[i];

            const double x_center =
                fm.x_min + (static_cast<double>(col) + 0.5) * fm.cell_size;

            const double y_center =
                fm.y_min +
                (static_cast<double>(fm.n_rows - row) - 0.5) * fm.cell_size;

            int coarse_col = static_cast<int>(
                std::floor((x_center - coarse_meta_.x_min) / options_.coarse_cell_size)
            );
            int coarse_row_from_south = static_cast<int>(
                std::floor((y_center - coarse_meta_.y_min) / options_.coarse_cell_size)
            );

            coarse_col = std::max(0, std::min(coarse_n_cols_ - 1, coarse_col));
            coarse_row_from_south =
                std::max(0, std::min(coarse_n_rows_ - 1, coarse_row_from_south));

            const int coarse_row = coarse_n_rows_ - 1 - coarse_row_from_south;
            const int coarse_flat = coarse_row * coarse_n_cols_ + coarse_col;

            auto found = param_id_of_coarse_flat.find(coarse_flat);
            if (found == param_id_of_coarse_flat.end())
            {
                const int new_id = static_cast<int>(param_rc_of_id_.size());
                param_id_of_coarse_flat.emplace(coarse_flat, new_id);

                coarse_rc rc;
                rc.row = coarse_row;
                rc.col = coarse_col;
                param_rc_of_id_.push_back(rc);
                coarse_cell_flat_of_id_.push_back(coarse_flat);
                support_count_of_param_.push_back(0);

                found = param_id_of_coarse_flat.find(coarse_flat);
            }

            coarse_cell_id_of_node_[i] = found->second;
            ++support_count_of_param_[static_cast<std::size_t>(found->second)];
        }

        if (param_rc_of_id_.empty())
        {
            throw std::runtime_error(
                "lsdtt_chi_ksn_inversion_cpu: no active coarse inversion cells were constructed."
            );
        }
    }

    template <class DenominatorFn>
    lsdtt_chi_ksn_inversion_result run_with_denominator(
        DenominatorFn&& denominator_fn)
    {
        const std::size_t n_nodes = network_.flats.size();
        const int n_params = static_cast<int>(param_rc_of_id_.size());

        lsdtt_chi_ksn_inversion_result result;
        result.network = network_;
        result.tested_m_over_n = tested_m_over_n_;
        result.trials.assign(
            tested_m_over_n_.size(),
            lsdtt_chi_ksn_trial_result{}
        );

        result.best_trial_index = -1;
        result.best_m_over_n = quiet_nan_d();
        result.best_mean_absolute_error = quiet_nan_f();
        result.best_root_mean_square_error = quiet_nan_f();

        result.chi.assign(n_nodes, quiet_nan_f());
        result.predicted_elevation.assign(n_nodes, quiet_nan_f());
        result.residual_elevation.assign(n_nodes, quiet_nan_f());
        result.node_u_star.assign(n_nodes, quiet_nan_f());
        result.coarse_cell_id_of_node = coarse_cell_id_of_node_;

        const float outlet_elevation = network_.observed_elevation[0];

        std::vector<float> chi_trial(n_nodes, 0.0f);
        std::vector<double> delta_chi(n_nodes, 0.0);
        std::vector<float> predicted_trial(n_nodes, quiet_nan_f());
        std::vector<float> residual_trial(n_nodes, quiet_nan_f());
        std::vector<float> node_u_star_trial(n_nodes, quiet_nan_f());

        std::vector<double> best_u_vector;
        best_u_vector.clear();

        for (std::size_t trial_idx = 0; trial_idx < tested_m_over_n_.size(); ++trial_idx)
        {
            const double mn = tested_m_over_n_[trial_idx];
            lsdtt_chi_ksn_trial_result& trial = result.trials[trial_idx];
            trial.m_over_n = mn;
            trial.n_nodes = static_cast<std::int32_t>(n_nodes);
            trial.n_parameters = n_params;
            trial.n_constrained_parameters = 0;
            trial.mean_absolute_error = quiet_nan_f();
            trial.root_mean_square_error = quiet_nan_f();

            bool trial_valid = true;

            try
            {
                // -------------------------------------------------------------
                // 1) Recompute chi on the fixed channel network
                // -------------------------------------------------------------
                chi_trial[0] = 0.0f;
                delta_chi[0] = 0.0;

                for (std::size_t i = 1; i < n_nodes; ++i)
                {
                    const int rec_local = receiver_local_[i];
                    const std::int64_t flat = network_.flats[i];

                    const double denom = denominator_fn(flat);
                    if (!is_finite(denom) || denom <= 0.0)
                    {
                        throw std::runtime_error(
                            "lsdtt_chi_ksn_inversion_cpu: denominator must be finite and > 0."
                        );
                    }

                    const double this_delta =
                        step_to_receiver_[i] *
                        std::pow(options_.A_0 / denom, mn);

                    if (!is_finite(this_delta))
                    {
                        throw std::runtime_error(
                            "lsdtt_chi_ksn_inversion_cpu: non-finite delta chi encountered."
                        );
                    }

                    delta_chi[i] = this_delta;
                    chi_trial[i] = static_cast<float>(
                        static_cast<double>(chi_trial[static_cast<std::size_t>(rec_local)]) +
                        this_delta
                    );
                }

                // -------------------------------------------------------------
                // 2) Assemble sparse design matrix A
                //
                // Each row corresponds to one node. Row i contains the summed
                // downstream delta-chi contributions grouped by coarse
                // inversion cell, so that z_rel = A * u.
                // -------------------------------------------------------------
                std::vector<Eigen::Triplet<double>> A_triplets;
                A_triplets.reserve(n_nodes * 8);

                std::vector<std::vector<std::pair<int,double>>> row_entries(n_nodes);
                std::vector<uint8_t> constrained_param(n_params, 0u);

                for (std::size_t i = 0; i < n_nodes; ++i)
                {
                    if (i == 0)
                    {
                        continue; // outlet row is zero in relative-elevation coordinates
                    }

                    const int rec_local = receiver_local_[i];
                    row_entries[i] = row_entries[static_cast<std::size_t>(rec_local)];

                    const int param_id = coarse_cell_id_of_node_[i];
                    if (param_id < 0 || param_id >= n_params)
                    {
                        throw std::runtime_error(
                            "lsdtt_chi_ksn_inversion_cpu: invalid coarse-cell parameter mapping."
                        );
                    }

                    bool found = false;
                    for (auto& item : row_entries[i])
                    {
                        if (item.first == param_id)
                        {
                            item.second += delta_chi[i];
                            found = true;
                            break;
                        }
                    }
                    if (!found)
                    {
                        row_entries[i].emplace_back(param_id, delta_chi[i]);
                    }

                    for (const auto& item : row_entries[i])
                    {
                        if (std::fabs(item.second) > 0.0)
                        {
                            A_triplets.emplace_back(
                                static_cast<int>(i),
                                item.first,
                                item.second
                            );
                            constrained_param[static_cast<std::size_t>(item.first)] = 1u;
                        }
                    }
                }

                for (const uint8_t flag : constrained_param)
                {
                    if (flag != 0u)
                    {
                        ++trial.n_constrained_parameters;
                    }
                }

                Eigen::SparseMatrix<double> A(
                    static_cast<int>(n_nodes),
                    n_params
                );
                A.setFromTriplets(A_triplets.begin(), A_triplets.end());
                A.makeCompressed();

                // -------------------------------------------------------------
                // 3) Assemble Laplacian roughness operator on active coarse cells
                // -------------------------------------------------------------
                std::unordered_map<int,int> param_id_of_coarse_flat;
                param_id_of_coarse_flat.reserve(static_cast<std::size_t>(n_params));
                for (int pid = 0; pid < n_params; ++pid)
                {
                    param_id_of_coarse_flat.emplace(coarse_cell_flat_of_id_[pid], pid);
                }

                std::vector<Eigen::Triplet<double>> L_triplets;
                L_triplets.reserve(static_cast<std::size_t>(n_params) * 5u);

                for (int pid = 0; pid < n_params; ++pid)
                {
                    const int coarse_flat = coarse_cell_flat_of_id_[pid];
                    const int row = coarse_flat / coarse_n_cols_;
                    const int col = coarse_flat % coarse_n_cols_;

                    int degree = 0;

                    const int nrow[4] = { row - 1, row + 1, row, row };
                    const int ncol[4] = { col, col, col - 1, col + 1 };

                    for (int k = 0; k < 4; ++k)
                    {
                        if (nrow[k] < 0 || nrow[k] >= coarse_n_rows_ ||
                            ncol[k] < 0 || ncol[k] >= coarse_n_cols_)
                        {
                            continue;
                        }

                        const int neigh_flat = nrow[k] * coarse_n_cols_ + ncol[k];
                        const auto found = param_id_of_coarse_flat.find(neigh_flat);
                        if (found == param_id_of_coarse_flat.end())
                        {
                            continue;
                        }

                        ++degree;
                        L_triplets.emplace_back(pid, found->second, -1.0);
                    }

                    L_triplets.emplace_back(pid, pid, static_cast<double>(degree));
                }

                Eigen::SparseMatrix<double> L(n_params, n_params);
                L.setFromTriplets(L_triplets.begin(), L_triplets.end());
                L.makeCompressed();

                // -------------------------------------------------------------
                // 4) Build z relative to outlet elevation
                // -------------------------------------------------------------
                Eigen::VectorXd z_rel(static_cast<int>(n_nodes));
                for (std::size_t i = 0; i < n_nodes; ++i)
                {
                    z_rel[static_cast<int>(i)] =
                        static_cast<double>(network_.observed_elevation[i] - outlet_elevation);
                }

                // -------------------------------------------------------------
                // 5) Solve regularized system
                //
                //   (A^TA + lambda L^TL + epsilon I) u = A^T z_rel
                //
                // Version 1 uses an unconstrained solve. If LDLT factorization
                // fails, a SparseQR solve on the stacked system is attempted.
                // -------------------------------------------------------------
                const Eigen::SparseMatrix<double> AtA = A.transpose() * A;
                const Eigen::SparseMatrix<double> LtL = L.transpose() * L;

                Eigen::SparseMatrix<double> H = AtA + options_.roughness_weight * LtL;

                if (options_.diagonal_jitter > 0.0)
                {
                    std::vector<Eigen::Triplet<double>> jitter_triplets;
                    jitter_triplets.reserve(static_cast<std::size_t>(n_params));
                    for (int pid = 0; pid < n_params; ++pid)
                    {
                        jitter_triplets.emplace_back(pid, pid, options_.diagonal_jitter);
                    }
                    Eigen::SparseMatrix<double> J(n_params, n_params);
                    J.setFromTriplets(jitter_triplets.begin(), jitter_triplets.end());
                    H = H + J;
                }

                const Eigen::VectorXd b = A.transpose() * z_rel;

                Eigen::VectorXd u(static_cast<int>(n_params));
                bool solved = false;

                {
                    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
                    solver.compute(H);
                    if (solver.info() == Eigen::Success)
                    {
                        u = solver.solve(b);
                        if (solver.info() == Eigen::Success)
                        {
                            solved = true;
                        }
                    }
                }

                if (!solved)
                {
                    const int n_aug_rows = static_cast<int>(n_nodes) + n_params + n_params;

                    std::vector<Eigen::Triplet<double>> B_triplets;
                    B_triplets.reserve(
                        static_cast<std::size_t>(A.nonZeros()) +
                        static_cast<std::size_t>(L.nonZeros()) +
                        static_cast<std::size_t>(n_params)
                    );

                    for (int k = 0; k < A.outerSize(); ++k)
                    {
                        for (Eigen::SparseMatrix<double>::InnerIterator it(A, k); it; ++it)
                        {
                            B_triplets.emplace_back(it.row(), it.col(), it.value());
                        }
                    }

                    const double sqrt_lambda = std::sqrt(std::max(0.0, options_.roughness_weight));
                    for (int k = 0; k < L.outerSize(); ++k)
                    {
                        for (Eigen::SparseMatrix<double>::InnerIterator it(L, k); it; ++it)
                        {
                            B_triplets.emplace_back(
                                static_cast<int>(n_nodes) + it.row(),
                                it.col(),
                                sqrt_lambda * it.value()
                            );
                        }
                    }

                    const double sqrt_jitter = std::sqrt(std::max(0.0, options_.diagonal_jitter));
                    if (sqrt_jitter > 0.0)
                    {
                        for (int pid = 0; pid < n_params; ++pid)
                        {
                            B_triplets.emplace_back(
                                static_cast<int>(n_nodes) + n_params + pid,
                                pid,
                                sqrt_jitter
                            );
                        }
                    }

                    Eigen::SparseMatrix<double> B(n_aug_rows, n_params);
                    B.setFromTriplets(B_triplets.begin(), B_triplets.end());
                    B.makeCompressed();

                    Eigen::VectorXd rhs = Eigen::VectorXd::Zero(n_aug_rows);
                    rhs.head(static_cast<int>(n_nodes)) = z_rel;

                    Eigen::SparseQR<
                        Eigen::SparseMatrix<double>,
                        Eigen::COLAMDOrdering<int>
                    > qr;

                    qr.compute(B);
                    if (qr.info() != Eigen::Success)
                    {
                        throw std::runtime_error(
                            "lsdtt_chi_ksn_inversion_cpu: sparse factorization failed for this trial."
                        );
                    }

                    u = qr.solve(rhs);
                    if (qr.info() != Eigen::Success)
                    {
                        throw std::runtime_error(
                            "lsdtt_chi_ksn_inversion_cpu: sparse solve failed for this trial."
                        );
                    }
                }

                for (int pid = 0; pid < n_params; ++pid)
                {
                    if (!std::isfinite(u[pid]))
                    {
                        throw std::runtime_error(
                            "lsdtt_chi_ksn_inversion_cpu: solver returned non-finite parameter values."
                        );
                    }
                }

                // -------------------------------------------------------------
                // 6) Predict, compute residuals, and trial misfit
                // -------------------------------------------------------------
                const Eigen::VectorXd pred_rel = A * u;

                double sum_abs = 0.0;
                double sum_sq  = 0.0;

                for (std::size_t i = 0; i < n_nodes; ++i)
                {
                    const float pred_abs = static_cast<float>(
                        pred_rel[static_cast<int>(i)] + static_cast<double>(outlet_elevation)
                    );

                    predicted_trial[i] = pred_abs;
                    residual_trial[i]  = network_.observed_elevation[i] - pred_abs;

                    const double abs_res = std::fabs(static_cast<double>(residual_trial[i]));
                    sum_abs += abs_res;
                    sum_sq  += abs_res * abs_res;
                }

                trial.mean_absolute_error = static_cast<float>(
                    sum_abs / static_cast<double>(n_nodes)
                );
                trial.root_mean_square_error = static_cast<float>(
                    std::sqrt(sum_sq / static_cast<double>(n_nodes))
                );

                for (std::size_t i = 0; i < n_nodes; ++i)
                {
                    const int pid = coarse_cell_id_of_node_[i];
                    if (pid >= 0 && pid < n_params)
                    {
                        node_u_star_trial[i] = static_cast<float>(u[pid]);
                    }
                }

                // -------------------------------------------------------------
                // 7) Retain best trial by selected metric
                // -------------------------------------------------------------
                const float this_score =
                    options_.use_mean_absolute_error_for_model_selection
                        ? trial.mean_absolute_error
                        : trial.root_mean_square_error;

                const bool take_this_trial =
                    (result.best_trial_index < 0) ||
                    (!std::isnan(this_score) &&
                     this_score <
                     (options_.use_mean_absolute_error_for_model_selection
                        ? result.best_mean_absolute_error
                        : result.best_root_mean_square_error));

                if (take_this_trial)
                {
                    result.best_trial_index = static_cast<std::int32_t>(trial_idx);
                    result.best_m_over_n = mn;
                    result.best_mean_absolute_error = trial.mean_absolute_error;
                    result.best_root_mean_square_error = trial.root_mean_square_error;

                    result.chi = chi_trial;
                    result.predicted_elevation = predicted_trial;
                    result.residual_elevation = residual_trial;
                    result.node_u_star = node_u_star_trial;

                    best_u_vector.assign(
                        u.data(),
                        u.data() + u.size()
                    );
                }
            }
            catch (...)
            {
                trial_valid = false;
            }

            if (!trial_valid)
            {
                trial.mean_absolute_error = quiet_nan_f();
                trial.root_mean_square_error = quiet_nan_f();
            }
        }

        // -------------------------------------------------------------
        // 8) Rasterize best coarse-grid solution and support counts
        // -------------------------------------------------------------
        {
            const std::size_t coarse_size =
                static_cast<std::size_t>(coarse_n_rows_) * coarse_n_cols_;

            std::vector<float> support_data(coarse_size, static_cast<float>(coarse_meta_.no_data_value));
            std::vector<uint8_t> support_mask(coarse_size, 0u);

            for (std::size_t pid = 0; pid < coarse_cell_flat_of_id_.size(); ++pid)
            {
                const std::size_t idx = static_cast<std::size_t>(coarse_cell_flat_of_id_[pid]);
                support_data[idx] = static_cast<float>(support_count_of_param_[pid]);
                support_mask[idx] = 1u;
            }

            result.support_count_field =
                lsdtt_raster_cpu(std::move(support_data), std::move(support_mask), coarse_meta_);
        }

        {
            const std::size_t coarse_size =
                static_cast<std::size_t>(coarse_n_rows_) * coarse_n_cols_;

            std::vector<float> u_data(coarse_size, static_cast<float>(coarse_meta_.no_data_value));
            std::vector<uint8_t> u_mask(coarse_size, 0u);

            if (result.best_trial_index >= 0 &&
                best_u_vector.size() == coarse_cell_flat_of_id_.size())
            {
                for (std::size_t pid = 0; pid < coarse_cell_flat_of_id_.size(); ++pid)
                {
                    const std::size_t idx = static_cast<std::size_t>(coarse_cell_flat_of_id_[pid]);
                    u_data[idx] = static_cast<float>(best_u_vector[pid]);
                    u_mask[idx] = 1u;
                }
            }

            result.u_star_field =
                lsdtt_raster_cpu(std::move(u_data), std::move(u_mask), coarse_meta_);
        }

        return result;
    }
};

} // anonymous namespace

// =============================================================================
// Helper builders
// =============================================================================
lsdtt_chi_ksn_inversion_network make_chi_ksn_inversion_network(
    const lsdtt_chi_basin_result& chi_basin,
    const lsdtt_raster_cpu& dem)
{
    const std::size_t n = chi_basin.flats.size();

    if (n == 0)
    {
        throw std::invalid_argument(
            "lsdtt_chi_ksn_inversion_cpu: chi_basin is empty."
        );
    }

    if (chi_basin.receiver_flats.size() != n)
    {
        throw std::invalid_argument(
            "lsdtt_chi_ksn_inversion_cpu: chi_basin has inconsistent vector lengths."
        );
    }

    lsdtt_chi_ksn_inversion_network network;
    network.outlet_flat = chi_basin.outlet_flat;
    network.flats = chi_basin.flats;
    network.receiver_flats = chi_basin.receiver_flats;
    network.observed_elevation.resize(n, 0.0f);

    for (std::size_t i = 0; i < n; ++i)
    {
        network.observed_elevation[i] = sample_dem_at_flat(dem, network.flats[i]);
    }

    return network;
}

lsdtt_chi_ksn_inversion_network make_chi_ksn_inversion_network(
    const lsdtt_chi_network_fit_result& fit_result)
{
    const std::size_t n = fit_result.flats.size();

    if (n == 0)
    {
        throw std::invalid_argument(
            "lsdtt_chi_ksn_inversion_cpu: fit_result is empty."
        );
    }

    if (fit_result.receiver_flats.size() != n ||
        fit_result.observed_elevation.size() != n)
    {
        throw std::invalid_argument(
            "lsdtt_chi_ksn_inversion_cpu: fit_result has inconsistent vector lengths."
        );
    }

    lsdtt_chi_ksn_inversion_network network;
    network.outlet_flat = fit_result.outlet_flat;
    network.flats = fit_result.flats;
    network.receiver_flats = fit_result.receiver_flats;
    network.observed_elevation = fit_result.observed_elevation;

    return network;
}

// =============================================================================
// Single-network inversion
// =============================================================================
lsdtt_chi_ksn_inversion_result invert_chi_ksn(
    const lsdtt_flow_base& flow,
    const lsdtt_chi_ksn_inversion_network& network,
    const std::vector<double>& tested_m_over_n,
    const lsdtt_chi_ksn_inversion_options& options)
{
    chi_ksn_inverter inverter(flow, network, tested_m_over_n, options);
    return inverter.run_area_based();
}

lsdtt_chi_ksn_inversion_result invert_chi_ksn(
    const lsdtt_flow_base& flow,
    const lsdtt_chi_ksn_inversion_network& network,
    const lsdtt_raster_cpu& denominator_raster,
    const std::vector<double>& tested_m_over_n,
    const lsdtt_chi_ksn_inversion_options& options)
{
    chi_ksn_inverter inverter(flow, network, tested_m_over_n, options);
    return inverter.run_custom_denominator(denominator_raster);
}

// =============================================================================
// Collection overloads
// =============================================================================
lsdtt_chi_ksn_inversion_collection_result invert_chi_ksn(
    const lsdtt_flow_base& flow,
    const std::vector<lsdtt_chi_ksn_inversion_network>& networks,
    const std::vector<double>& tested_m_over_n,
    const lsdtt_chi_ksn_inversion_options& options)
{
    lsdtt_chi_ksn_inversion_collection_result result;
    result.basins.reserve(networks.size());

    for (const auto& network : networks)
    {
        result.basins.push_back(
            invert_chi_ksn(flow, network, tested_m_over_n, options)
        );
    }

    return result;
}

lsdtt_chi_ksn_inversion_collection_result invert_chi_ksn(
    const lsdtt_flow_base& flow,
    const std::vector<lsdtt_chi_ksn_inversion_network>& networks,
    const lsdtt_raster_cpu& denominator_raster,
    const std::vector<double>& tested_m_over_n,
    const lsdtt_chi_ksn_inversion_options& options)
{
    lsdtt_chi_ksn_inversion_collection_result result;
    result.basins.reserve(networks.size());

    for (const auto& network : networks)
    {
        result.basins.push_back(
            invert_chi_ksn(flow, network, denominator_raster, tested_m_over_n, options)
        );
    }

    return result;
}

void write_chi_ksn_inversion_trials_csv(
    const std::string& filename,
    const lsdtt_chi_ksn_inversion_collection_result& results)
{
    std::ofstream out(filename);
    if (!out)
        throw std::runtime_error(
            "write_chi_ksn_inversion_trials_csv: cannot open " + filename);

    out << "basin_index,outlet_flat,trial_index,m_over_n,mean_absolute_error,"
           "root_mean_square_error,n_nodes,n_parameters,n_constrained_parameters\n";

    for (std::size_t b = 0; b < results.basins.size(); ++b)
    {
        const auto& basin = results.basins[b];
        for (std::size_t t = 0; t < basin.trials.size(); ++t)
        {
            const auto& tr = basin.trials[t];
            out << static_cast<std::int32_t>(b) << ','
                << basin.network.outlet_flat << ','
                << static_cast<std::int32_t>(t) << ','
                << tr.m_over_n << ','
                << tr.mean_absolute_error << ','
                << tr.root_mean_square_error << ','
                << tr.n_nodes << ','
                << tr.n_parameters << ','
                << tr.n_constrained_parameters << '\n';
        }
    }
}

void write_chi_ksn_inversion_best_summary_csv(
    const std::string& filename,
    const lsdtt_chi_ksn_inversion_collection_result& results)
{
    std::ofstream out(filename);
    if (!out)
        throw std::runtime_error(
            "write_chi_ksn_inversion_best_summary_csv: cannot open " + filename);

    out << "basin_index,outlet_flat,best_trial_index,best_m_over_n,"
           "best_mean_absolute_error,best_root_mean_square_error,"
           "n_nodes_best,n_parameters_best,n_constrained_parameters_best\n";

    for (std::size_t b = 0; b < results.basins.size(); ++b)
    {
        const auto& basin = results.basins[b];

        std::int32_t n_nodes_best  = 0;
        std::int32_t n_params_best = 0;
        std::int32_t n_const_best  = 0;

        if (basin.best_trial_index >= 0 &&
            static_cast<std::size_t>(basin.best_trial_index) < basin.trials.size())
        {
            const auto& tr =
                basin.trials[static_cast<std::size_t>(basin.best_trial_index)];
            n_nodes_best  = tr.n_nodes;
            n_params_best = tr.n_parameters;
            n_const_best  = tr.n_constrained_parameters;
        }

        out << static_cast<std::int32_t>(b) << ','
            << basin.network.outlet_flat << ','
            << basin.best_trial_index << ','
            << basin.best_m_over_n << ','
            << basin.best_mean_absolute_error << ','
            << basin.best_root_mean_square_error << ','
            << n_nodes_best << ','
            << n_params_best << ','
            << n_const_best << '\n';
    }
}

} // namespace lsdtt
