// lsdtt_ridgeline_cpu.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Rutger Vervoordeldonk
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html

#include "lsdtt_ridgeline_cpu.hpp"

#include "lsdtt_flowinfo_cpu.hpp" // For dynamic_cast accessing CPU-specific config
#include "../utils/lsdtt_regression.hpp"
#include "../utils/lsdtt_signal1d.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <unordered_set>
#include <set>
#include <sstream>
#include <iomanip>
#include <cctype>

namespace lsdtt
{

namespace
{
// Helper to format floats closely matching precision of LSDTopoTools2
std::string format_float(double value)
{
    if (std::isnan(value)) return "NaN";
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(9) << value;
    return oss.str();
}


struct EdgePixels {
    int left_r, left_c;
    int right_r, right_c;
    bool valid;
};

// Gets the row/col DEM pixels strictly to the Left and Right of an oriented 4-connected divide grid step
EdgePixels get_edge_pixels(const lsdtt_pixel_rc& p1, const lsdtt_pixel_rc& p2) {
    EdgePixels ep;
    ep.valid = true;
    int dr = p2.row - p1.row;
    int dc = p2.col - p1.col;

    if (dr == 0 && dc == 1) { // East
        ep.left_r = p1.row - 1; ep.left_c = p1.col;
        ep.right_r = p1.row;    ep.right_c = p1.col;
    } else if (dr == 0 && dc == -1) { // West
        ep.left_r = p1.row;     ep.left_c = p1.col - 1;
        ep.right_r = p1.row - 1; ep.right_c = p1.col - 1;
    } else if (dr == 1 && dc == 0) { // South
        ep.left_r = p1.row;     ep.left_c = p1.col;
        ep.right_r = p1.row;    ep.right_c = p1.col - 1;
    } else if (dr == -1 && dc == 0) { // North
        ep.left_r = p1.row - 1; ep.left_c = p1.col - 1;
        ep.right_r = p1.row - 1; ep.right_c = p1.col;
    } else {
        ep.valid = false; // Not a strict 4-connected step
    }
    return ep;
}
}

// =============================================================================
// Construction
// =============================================================================

lsdtt_ridgeline_cpu::lsdtt_ridgeline_cpu(
    std::vector<lsdtt_pixel_rc> nodes,
    const lsdtt_raster_meta&    divide_grid_meta)
    : divide_grid_meta_(divide_grid_meta)
    , nodes_(std::move(nodes))
{
    compute_positions_and_distances_();
}

void lsdtt_ridgeline_cpu::compute_positions_and_distances_()
{
    positions_.clear();
    distances_.clear();

    if (nodes_.empty()) return;

    positions_.reserve(nodes_.size());
    distances_.reserve(nodes_.size());

    for (const lsdtt_pixel_rc& rc : nodes_)
    {
        positions_.push_back(cell_center_from_row_col(divide_grid_meta_, rc.row, rc.col));
    }

    distances_.push_back(0.0);
    for (std::size_t i = 1; i < positions_.size(); ++i)
    {
        double dx = positions_[i].x - positions_[i - 1].x;
        double dy = positions_[i].y - positions_[i - 1].y;

        // Correct coordinate jumps when crossing periodic boundaries
        int d_row = nodes_[i].row - nodes_[i - 1].row;
        int d_col = nodes_[i].col - nodes_[i - 1].col;

        bool wrapped = false;
        int unwrapped_row = nodes_[i].row;
        int unwrapped_col = nodes_[i].col;

        if (divide_grid_meta_.n_rows > 0 && std::abs(d_row) > divide_grid_meta_.n_rows / 2)
        {
            int true_d_row = (d_row > 0) ? d_row - divide_grid_meta_.n_rows : d_row + divide_grid_meta_.n_rows;
            unwrapped_row = nodes_[i - 1].row + true_d_row;
            wrapped = true;
        }

        if (divide_grid_meta_.n_cols > 0 && std::abs(d_col) > divide_grid_meta_.n_cols / 2)
        {
            int true_d_col = (d_col > 0) ? d_col - divide_grid_meta_.n_cols : d_col + divide_grid_meta_.n_cols;
            unwrapped_col = nodes_[i - 1].col + true_d_col;
            wrapped = true;
        }

        // If a boundary was crossed, compute dx/dy using the locally un-wrapped coordinates
        if (wrapped)
        {
            lsdtt_xy unwrapped_pos = cell_center_from_row_col(divide_grid_meta_, unwrapped_row, unwrapped_col);
            dx = unwrapped_pos.x - positions_[i - 1].x;
            dy = unwrapped_pos.y - positions_[i - 1].y;
        }

        distances_.push_back(distances_.back() + std::sqrt(dx * dx + dy * dy));
    }
}

bool lsdtt_ridgeline_cpu::check_if_nodes_four_connected(bool periodic_ns, bool periodic_ew) const
{
    if (nodes_.size() < 2) return true;

    for (std::size_t i = 1; i < nodes_.size(); ++i)
    {
        int prev_row = nodes_[i - 1].row;
        int prev_col = nodes_[i - 1].col;
        int curr_row = nodes_[i].row;
        int curr_col = nodes_[i].col;

        int d_row = std::abs(curr_row - prev_row);
        int d_col = std::abs(curr_col - prev_col);

        bool row_connected = (d_row == 1 && d_col == 0);
        bool col_connected = (d_row == 0 && d_col == 1);

        if (periodic_ns && d_row == divide_grid_meta_.n_rows - 2 && d_col == 0) row_connected = true;
        if (periodic_ew && d_col == divide_grid_meta_.n_cols - 2 && d_row == 0) col_connected = true;

        if (!(row_connected || col_connected)) return false;
    }
    return true;
}

// =============================================================================
// Elevation filling
// =============================================================================

void lsdtt_ridgeline_cpu::fill_elevation(const lsdtt_raster_cpu& dem)
{
    elevations_.clear();
    elevations_.reserve(nodes_.size());

    const lsdtt_raster_meta& dem_meta = dem.get_meta();
    const auto& data = dem.get_data_array();
    const auto& mask = dem.get_mask_array();
    const float ndv  = dem_meta.no_data_value;

    const int dem_nrows = dem_meta.n_rows;
    const int dem_ncols = dem_meta.n_cols;

    for (const lsdtt_pixel_rc& dg_rc : nodes_)
    {
        float sum_elev = 0.0f;
        int   n_valid  = 0;

        const int offsets[4][2] = { {-1,-1}, {-1, 0}, {0,-1}, {0, 0} };

        for (const auto& off : offsets)
        {
            const int r = dg_rc.row + off[0];
            const int c = dg_rc.col + off[1];

            if (r < 0 || r >= dem_nrows || c < 0 || c >= dem_ncols) continue;
            if (mask(r, c) == 0) continue;

            sum_elev += data(r, c);
            ++n_valid;
        }

        elevations_.push_back(n_valid > 0 ? (sum_elev / static_cast<float>(n_valid)) : ndv);
    }
}

std::vector<float> lsdtt_ridgeline_cpu::sample_raster_values(const lsdtt_raster_cpu& raster) const
{
    std::vector<float> sampled_vals;
    sampled_vals.reserve(nodes_.size());

    const lsdtt_raster_meta& raster_meta = raster.get_meta();
    const auto& data = raster.get_data_array();
    const auto& mask = raster.get_mask_array();
    const float ndv  = raster_meta.no_data_value;

    const int nrows = raster_meta.n_rows;
    const int ncols = raster_meta.n_cols;

    for (const lsdtt_pixel_rc& dg_rc : nodes_)
    {
        float sum_val = 0.0f;
        int   n_valid = 0;

        const int offsets[4][2] = { {-1,-1}, {-1, 0}, {0,-1}, {0, 0} };

        for (const auto& off : offsets)
        {
            const int r = dg_rc.row + off[0];
            const int c = dg_rc.col + off[1];

            if (r < 0 || r >= nrows || c < 0 || c >= ncols) continue;
            if (mask(r, c) == 0) continue;

            sum_val += data(r, c);
            ++n_valid;
        }

        sampled_vals.push_back(n_valid > 0 ? (sum_val / static_cast<float>(n_valid)) : ndv);
    }
    return sampled_vals;
}

// =============================================================================
// Mutating operations
// =============================================================================

void lsdtt_ridgeline_cpu::reverse()
{
    std::reverse(nodes_.begin(), nodes_.end());
    std::reverse(positions_.begin(), positions_.end());
    std::reverse(elevations_.begin(), elevations_.end());
    
    if (!fitted_elevations_.empty())
    {
        std::reverse(fitted_elevations_.begin(), fitted_elevations_.end());
        std::reverse(detrended_elevations_.begin(), detrended_elevations_.end());
        std::reverse(smoothed_detrended_elevations_.begin(), smoothed_detrended_elevations_.end());
        std::reverse(smoothed_detrended_slopes_.begin(), smoothed_detrended_slopes_.end());
    }

    std::reverse(nearest_channel_elevations_left_.begin(), nearest_channel_elevations_left_.end());
    std::reverse(nearest_channel_elevations_right_.begin(), nearest_channel_elevations_right_.end());

    compute_positions_and_distances_();
}

void lsdtt_ridgeline_cpu::orient_periodic_start(bool columns)
{
    if (nodes_.empty()) return;
    auto rotate_to_minimum = [&]()
    {
        const auto index = std::min_element(
            nodes_.begin(), nodes_.end(),
            [columns](const lsdtt_pixel_rc& a, const lsdtt_pixel_rc& b)
            { return columns ? a.col < b.col : a.row < b.row; });
        const auto offset = static_cast<std::size_t>(index - nodes_.begin());
        if (offset == 0) return;
        std::rotate(nodes_.begin(), nodes_.begin() + offset, nodes_.end());
        std::rotate(positions_.begin(), positions_.begin() + offset, positions_.end());
        if (!elevations_.empty())
            std::rotate(elevations_.begin(), elevations_.begin() + offset, elevations_.end());
        if (!fitted_elevations_.empty())
            std::rotate(fitted_elevations_.begin(), fitted_elevations_.begin() + offset, fitted_elevations_.end());
        if (!detrended_elevations_.empty())
            std::rotate(detrended_elevations_.begin(), detrended_elevations_.begin() + offset, detrended_elevations_.end());
        if (!smoothed_detrended_elevations_.empty())
            std::rotate(smoothed_detrended_elevations_.begin(), smoothed_detrended_elevations_.begin() + offset, smoothed_detrended_elevations_.end());
        if (!smoothed_detrended_slopes_.empty())
            std::rotate(smoothed_detrended_slopes_.begin(), smoothed_detrended_slopes_.begin() + offset, smoothed_detrended_slopes_.end());
    };
    rotate_to_minimum();
    const auto coordinate = [columns](const lsdtt_pixel_rc& node)
    { return columns ? node.col : node.row; };
    const int coordinate_period = columns ? divide_grid_meta_.n_cols
                                          : divide_grid_meta_.n_rows;
    for (std::size_t i = 1; i < nodes_.size(); ++i)
    {
        const int delta = coordinate(nodes_[i]) - coordinate(nodes_[0]);
        if (delta != 0 && (delta < 0 || delta > coordinate_period / 2))
        {
            reverse();
            rotate_to_minimum();
            break;
        }
        if (delta != 0) break;
    }

    compute_positions_and_distances_();
}

void lsdtt_ridgeline_cpu::calculate_nearest_channel_elevations(
    const lsdtt_flow_base& flow,
    const lsdtt_raster_cpu& dem,
    const std::vector<std::uint8_t>& channel_mask)
{
    if (channel_mask.size() != dem.full_grid_size())
        throw std::invalid_argument("channel mask does not match DEM size");

    nearest_channel_elevations_left_.assign(
        nodes_.size(), std::numeric_limits<float>::quiet_NaN());
    nearest_channel_elevations_right_.assign(
        nodes_.size(), std::numeric_limits<float>::quiet_NaN());

    const auto& data = dem.get_data();
    const auto& dem_mask = dem.get_mask();

    auto trace_to_channel = [&](int row, int col) -> float {
        if (!flow.in_bounds(row, col)) return std::numeric_limits<float>::quiet_NaN();
        std::int64_t current = flow.flat_index_unchecked(row, col);
        std::size_t steps = 0;
        while (flow.is_valid_cell(current) && steps++ <= channel_mask.size())
        {
            if (channel_mask[static_cast<std::size_t>(current)] != 0 &&
                dem_mask[static_cast<std::size_t>(current)] != 0)
                return data[static_cast<std::size_t>(current)];

            const std::int64_t next = flow.receiver_of_unchecked(current);
            if (next < 0 || next == current) break;
            current = next;
        }
        return std::numeric_limits<float>::quiet_NaN();
    };

    for (std::size_t i = 0; i + 1 < nodes_.size(); ++i)
    {
        const EdgePixels edge = get_edge_pixels(nodes_[i], nodes_[i + 1]);
        if (!edge.valid) continue;
        nearest_channel_elevations_left_[i] =
            trace_to_channel(edge.left_r, edge.left_c);
        nearest_channel_elevations_right_[i] =
            trace_to_channel(edge.right_r, edge.right_c);
    }
}

void lsdtt_ridgeline_cpu::append_without_check(const lsdtt_ridgeline_cpu& other)
{
    if (other.empty()) return;

    nodes_.insert(nodes_.end(), other.nodes_.begin(), other.nodes_.end());

    if (!elevations_.empty() && !other.elevations_.empty())
    {
        elevations_.insert(elevations_.end(), other.elevations_.begin(), other.elevations_.end());
    }
    else
    {
        elevations_.clear();
    }

    compute_positions_and_distances_();
}

// =============================================================================
// MDD assembly helper
// =============================================================================

std::vector<lsdtt_ridgeline_cpu>
lsdtt_ridgeline_cpu::extract_contiguous_shared_ridgelines(const lsdtt_ridgeline_cpu& other) const
{
    if (empty() || other.empty()) return {};

    auto node_hash = [](const lsdtt_pixel_rc& rc) -> std::size_t {
        const std::size_t r = static_cast<std::size_t>(rc.row);
        const std::size_t c = static_cast<std::size_t>(rc.col);
        return ((r + c) * (r + c + 1)) / 2 + c;
    };

    std::unordered_set<std::size_t> other_set;
    other_set.reserve(other.nodes_.size() * 2);
    for (const lsdtt_pixel_rc& rc : other.nodes_)
    {
        other_set.insert(node_hash(rc));
    }

    std::vector<lsdtt_ridgeline_cpu> result;
    std::vector<lsdtt_pixel_rc> current_run;

    for (const lsdtt_pixel_rc& rc : nodes_)
    {
        if (other_set.count(node_hash(rc)) != 0)
        {
            current_run.push_back(rc);
        }
        else
        {
            if (!current_run.empty())
            {
                result.emplace_back(std::move(current_run), divide_grid_meta_);
                current_run.clear();
            }
        }
    }

    if (!current_run.empty())
    {
        result.emplace_back(std::move(current_run), divide_grid_meta_);
    }

    return result;
}

// =============================================================================
// Ridgeline analysis tools
// =============================================================================

void lsdtt_ridgeline_cpu::segment_ridgeline_snout_by_curvature(float smoothing_window_fraction)
{
    if (distances_.empty() || elevations_.empty()) return;

    double bandwidth = distances_.back() * smoothing_window_fraction;

    std::vector<double> smoothed_elev = compute_lowess_smoothing(
        distances_, 
        std::vector<double>(elevations_.begin(), elevations_.end()), 
        bandwidth, 
        3
    );

    std::vector<std::size_t> knees = find_knees(distances_, smoothed_elev, 1.0);
    snoutidx_ = knees.empty() ? 0 : static_cast<int>(knees.front());
}

void lsdtt_ridgeline_cpu::fit_ridgeline_profile(
    bool use_polynomial, bool use_constant, int moving_average_window_pixels,
    bool additional_output)
{
    if (static_cast<std::size_t>(snoutidx_) >= nodes_.size()) return;
    if (moving_average_window_pixels < 1 ||
        (moving_average_window_pixels > 1 && moving_average_window_pixels % 2 == 0))
    {
        throw std::invalid_argument(
            "ridgeline detrend moving-average window must be a positive odd number of pixels");
    }
    const int configured_moving_average_window = moving_average_window_pixels;

    std::vector<double> dist_main(distances_.begin() + snoutidx_, distances_.end());
    std::vector<double> elev_main(elevations_.begin() + snoutidx_, elevations_.end());

    if (!additional_output)
    {
        fitted_elevations_.assign(nodes_.size(), divide_grid_meta_.no_data_value);
        detrended_elevations_.assign(nodes_.size(), divide_grid_meta_.no_data_value);
    }

    // Always compute polynomial fit coefficients for the statistics
    lsdtt_polynomial2_regression_result poly_fit = compute_polynomial2_regression(dist_main, elev_main);
    poly2_A_ = static_cast<float>(poly_fit.A);
    poly2_B_ = static_cast<float>(poly_fit.B);
    poly2_C_ = static_cast<float>(poly_fit.C);
    poly2_r2_ = static_cast<float>(poly_fit.r_squared);

    // Always compute linear fit for the slope and r2 statistics
    lsdtt_linear_regression_result lin_fit = compute_linear_regression(dist_main, elev_main);
    slope_fit_     = static_cast<float>(lin_fit.slope);
    intercept_fit_ = static_cast<float>(lin_fit.intercept);
    r2_fit_        = static_cast<float>(lin_fit.r_squared);

    if (additional_output)
    {
        auto existing = std::find(moving_window_sizes_pixels_.begin(),
                                  moving_window_sizes_pixels_.end(),
                                  configured_moving_average_window);
        if (existing != moving_window_sizes_pixels_.end())
            throw std::invalid_argument("duplicate moving-window ridgeline fit size");
        std::vector<float> moving_fitted(nodes_.size(), divide_grid_meta_.no_data_value);
        std::vector<float> moving_detrended(nodes_.size(), divide_grid_meta_.no_data_value);
        const int radius = configured_moving_average_window / 2;
        const std::size_t profile_size = nodes_.size() - snoutidx_;
        for (std::size_t i = snoutidx_; i < nodes_.size(); ++i)
        {
            const std::size_t local_i = i - snoutidx_;
            std::size_t first = 0;
            std::size_t last = 0;
            if (profile_size >= static_cast<std::size_t>(configured_moving_average_window))
            {
                first = local_i > static_cast<std::size_t>(radius)
                    ? local_i - static_cast<std::size_t>(radius) : 0;
                if (first + static_cast<std::size_t>(configured_moving_average_window) > profile_size)
                    first = profile_size - static_cast<std::size_t>(configured_moving_average_window);
                last = first + static_cast<std::size_t>(configured_moving_average_window) - 1;
            }
            else
            {
                first = 0;
                last = profile_size - 1;
            }
            double sum = 0.0;
            std::size_t count = 0;
            for (std::size_t j = first; j <= last; ++j)
            {
                const float z = elevations_[snoutidx_ + j];
                if (std::isfinite(z) && !divide_grid_meta_.is_nodata(z))
                {
                    sum += z;
                    ++count;
                }
            }
            const float fit_z = count > 0
                ? static_cast<float>(sum / static_cast<double>(count))
                : elevations_[i];
            moving_fitted[i] = fit_z;
            moving_detrended[i] = elevations_[i] - fit_z;
        }
        moving_window_sizes_pixels_.push_back(configured_moving_average_window);
        moving_window_fitted_elevations_.push_back(std::move(moving_fitted));
        moving_window_detrended_elevations_.push_back(std::move(moving_detrended));
    }
    if (additional_output) return;

    if (use_constant)
    {
        float mean_z = static_cast<float>(
            std::accumulate(elev_main.begin(), elev_main.end(), 0.0) / elev_main.size()
        );
        for (std::size_t i = snoutidx_; i < nodes_.size(); ++i)
        {
            fitted_elevations_[i] = mean_z;
            detrended_elevations_[i] = elevations_[i] - mean_z;
        }
    }
    else if (use_polynomial)
    {
        for (std::size_t i = snoutidx_; i < nodes_.size(); ++i)
        {
            float d = static_cast<float>(distances_[i] - distances_[snoutidx_]);
            float fit_z = poly2_A_ * d * d + poly2_B_ * d + poly2_C_;
            fitted_elevations_[i] = fit_z;
            detrended_elevations_[i] = elevations_[i] - fit_z;
        }
    }
    else
    {
        for (std::size_t i = snoutidx_; i < nodes_.size(); ++i)
        {
            float d = static_cast<float>(distances_[i] - distances_[snoutidx_]);
            float fit_z = slope_fit_ * d + intercept_fit_;
            fitted_elevations_[i] = fit_z;
            detrended_elevations_[i] = elevations_[i] - fit_z;
        }
    }
}

void lsdtt_ridgeline_cpu::smoothen_and_derive_detrended_ridgeline(int smoothing_window_pixels)
{
    if (static_cast<std::size_t>(snoutidx_) >= nodes_.size() || detrended_elevations_.empty()) return;

    std::vector<double> rl_dist;
    std::vector<double> rl_elev;

    bool periodic = false;
    if (nodes_.size() > 1 && nodes_[snoutidx_].row == nodes_.back().row && nodes_[snoutidx_].col == nodes_.back().col)
    {
        periodic = true;
    }

    std::size_t target_nodes = nodes_.size() - snoutidx_ - (periodic ? 1 : 0);
    
    for (std::size_t i = snoutidx_; i < snoutidx_ + target_nodes; ++i)
    {
        rl_dist.push_back(distances_[i]);
        rl_elev.push_back(detrended_elevations_[i]);
    }

    std::vector<double> smoothed_elev;
    if (smoothing_window_pixels == 1)
    {
        smoothed_elev = rl_elev;
    }
    else
    {
        smoothed_elev = compute_moving_average_smoothing(rl_elev, smoothing_window_pixels, periodic);
    }

    std::vector<double> derivs = compute_central_difference(rl_dist, smoothed_elev, periodic);

    smoothed_detrended_elevations_.assign(nodes_.size(), divide_grid_meta_.no_data_value);
    smoothed_detrended_slopes_.assign(nodes_.size(), divide_grid_meta_.no_data_value);

    for (std::size_t i = 0; i < target_nodes; ++i)
    {
        smoothed_detrended_elevations_[snoutidx_ + i] = static_cast<float>(smoothed_elev[i]);
        smoothed_detrended_slopes_[snoutidx_ + i] = static_cast<float>(derivs[i]);
    }

    if (periodic)
    {
        smoothed_detrended_elevations_.back() = static_cast<float>(smoothed_elev.front());
        smoothed_detrended_slopes_.back()     = static_cast<float>(derivs.front());
    }
}

float lsdtt_ridgeline_cpu::calculate_ridgeline_steepness(float basin_drainage_area, float concavity) const
{
    if (std::isnan(slope_fit_)) return std::numeric_limits<float>::quiet_NaN();
    return slope_fit_ / std::pow(basin_drainage_area, -concavity);
}

std::vector<lsdtt_ridgeline_extremum> lsdtt_ridgeline_cpu::get_ridgeline_extrema(
    const lsdtt_flow_base&             flow, 
    const lsdtt_junction_network_base& junction_network
) const
{
    std::vector<lsdtt_ridgeline_extremum> extrema;
    if (static_cast<std::size_t>(snoutidx_) >= nodes_.size()) return extrema;

    bool periodic = (nodes_[snoutidx_].row == nodes_.back().row && nodes_[snoutidx_].col == nodes_.back().col);
    int N = static_cast<int>(nodes_.size() - snoutidx_ - (periodic ? 1 : 0));
    
    if (N < 2) return extrema;

    std::vector<int> extremum_indices;
    std::vector<std::string> extremum_types;

    for (int i = 0; i < N; ++i)
    {
        bool is_min = true;
        bool is_max = true;

        int prev = (i - 1 >= 0) ? i - 1 : (periodic ? N - 1 : -1);
        int next = (i + 1 < N) ? i + 1 : (periodic ? 0 : -1);

        float curr_elev = smoothed_detrended_elevations_[snoutidx_ + i];
        
        if (prev != -1)
        {
            if (curr_elev >= smoothed_detrended_elevations_[snoutidx_ + prev]) is_min = false;
            if (curr_elev <= smoothed_detrended_elevations_[snoutidx_ + prev]) is_max = false;
        }
        if (next != -1)
        {
            if (curr_elev >= smoothed_detrended_elevations_[snoutidx_ + next]) is_min = false;
            if (curr_elev <= smoothed_detrended_elevations_[snoutidx_ + next]) is_max = false;
        }

        if (is_min)
        {
            extremum_indices.push_back(i);
            extremum_types.push_back("minimum");
        }
        else if (is_max)
        {
            extremum_indices.push_back(i);
            extremum_types.push_back("maximum");
        }
    }

    for (std::size_t i = 0; i < extremum_indices.size(); ++i)
    {
        int global_idx = snoutidx_ + extremum_indices[i];
        
        double d_next = divide_grid_meta_.no_data_value;
        double r_next = divide_grid_meta_.no_data_value;

        if (i < extremum_indices.size() - 1)
        {
            d_next = distances_[snoutidx_ + extremum_indices[i + 1]] - distances_[global_idx];
            r_next = std::abs(
                smoothed_detrended_elevations_[snoutidx_ + extremum_indices[i + 1]] - 
                smoothed_detrended_elevations_[global_idx]
            );
        }
        else if (periodic && extremum_indices.size() > 1)
        {
            d_next = (distances_[snoutidx_ + N] - distances_[global_idx]) + 
                     (distances_[snoutidx_ + extremum_indices.front()] - distances_[snoutidx_]);
            
            r_next = std::abs(
                smoothed_detrended_elevations_[snoutidx_ + extremum_indices.front()] -
                smoothed_detrended_elevations_[global_idx]
            );
        }

        float tdj_dist = find_distance_to_closest_tdj(global_idx, flow, junction_network);

        extrema.push_back({ global_idx, extremum_types[i], d_next, r_next, tdj_dist });
    }

    return extrema;
}

float lsdtt_ridgeline_cpu::find_distance_to_closest_tdj(
    int                                node_index,
    const lsdtt_flow_base&             flow,
    const lsdtt_junction_network_base& junction_network
) const
{
    // Safe fallback tracking (no NaN == NaN issues)
    float distance_to_tdj = -9999.0f;
    bool  found_tdj       = false;

    // Bounds check on the ridgeline node index
    if (node_index < 0 || static_cast<std::size_t>(node_index) >= distances_.size()) 
    {
        return -9999.0f;
    }

    const float distance_at_node = static_cast<float>(distances_[node_index]);

    // Helper lambda: Safely retrieves the upstream junction ID for any given DEM pixel
    auto get_basin_junc = [&](int r, int c) -> std::int32_t 
    {
        // Must stay inside DEM bounds
        if (!flow.in_bounds(r, c)) return -1;
        
        std::int64_t current = flow.flat_index_unchecked(r, c);

        // Discard invalid / masked cells entirely
        if (!flow.is_valid_cell(current)) return -1;

        // Route strictly downstream until we hit the channel network (or exit the DEM)
        while (!flow.is_outlet_unchecked(current) && !junction_network.is_channel(current)) 
        {
            current = flow.receiver_of_unchecked(current);

            // Catch disconnected/faulty topology loops just in case
            if (current < 0) break; 
        }

        // If we securely landed on an active channel link, grab its upstream basin junction
        if (current >= 0 && junction_network.is_channel(current)) 
        {
            return junction_network.junction_upstream_of_flat(current); 
        }
        
        return -1;
    };

    // Helper lambda: Scans forward or backward along the ridgeline
    auto scan_tdj = [&](int start_idx, int end_idx, int step) 
    {
        for (int i = start_idx; i != end_idx; i += step) 
        {
            const lsdtt_pixel_rc& rc = nodes_[i];
            std::unordered_set<std::int32_t> basin_juncs;

            // Divide-grid node (R,C) sits geometrically between 4 DEM pixels: 
            // (R-1, C-1), (R-1, C), (R, C-1), (R, C)
            const int r_offsets[2] = {rc.row - 1, rc.row};
            const int c_offsets[2] = {rc.col - 1, rc.col};

            for (int dr : r_offsets) 
            {
                for (int dc : c_offsets) 
                {
                    std::int32_t bj = get_basin_junc(dr, dc);
                    if (bj >= 0) 
                    {
                        basin_juncs.insert(bj);
                    }
                }
            }

            // A Triple (or Quadruple) Divide Junction is detected if 3+ unique basins meet here
            if (basin_juncs.size() > 2) 
            {
                const float thisdist = std::abs(static_cast<float>(distances_[i]) - distance_at_node);
                if (!found_tdj || thisdist < distance_to_tdj) 
                {
                    distance_to_tdj = thisdist;
                    found_tdj = true;
                }
                break; // Stop looking further in this direction, we found the closest TDJ
            }
        }
    };

    // Scan sequentially Forwards
    scan_tdj(node_index, static_cast<int>(distances_.size()), 1);
    
    // Scan sequentially Backwards
    scan_tdj(node_index, -1, -1);

    // Return the minimum found distance, or default to -9999.0
    return found_tdj ? distance_to_tdj : -9999.0f;
}

std::vector<lsdtt_ridgeline_triple_divide> lsdtt_ridgeline_cpu::get_triple_divides(
    const lsdtt_flow_base&             flow,
    const lsdtt_junction_network_base& junction_network
) const
{
    std::vector<lsdtt_ridgeline_triple_divide> result;
    for (std::size_t i = 1; i + 1 < nodes_.size(); ++i)
    {
        std::set<std::int32_t> basins;
        const lsdtt_pixel_rc& rc = nodes_[i];
        lsdtt_ridgeline_triple_divide divide;
        divide.node_index = static_cast<int>(i);
        const std::array<std::pair<int, int>, 4> offsets{{
            {-1, -1}, {-1, 0}, {0, 0}, {0, -1}}}; // NW, NE, SE, SW
        for (std::size_t pixel = 0; pixel < offsets.size(); ++pixel)
        {
            const int r = rc.row + offsets[pixel].first;
            const int c = rc.col + offsets[pixel].second;
            if (!flow.in_bounds(r, c)) continue;
            std::int64_t flat = flow.flat_index_unchecked(r, c);
            if (!flow.is_valid_cell(flat)) continue;
            while (!flow.is_outlet_unchecked(flat) &&
                   !junction_network.is_channel(flat))
            {
                flat = flow.receiver_of_unchecked(flat);
                if (flat < 0) break;
            }
            if (flat >= 0 && junction_network.is_channel(flat))
            {
                const std::int32_t j =
                    junction_network.junction_upstream_of_flat(flat);
                if (j >= 0)
                {
                    divide.surrounding_junctions[pixel] = j;
                    basins.insert(j);
                }
            }
        }
        if (basins.size() > 2) result.push_back(divide);
    }
    return result;
}

std::vector<lsdtt_ridgeline_topological_triple_divide>
lsdtt_ridgeline_cpu::get_topological_distance_triple_divides(
    const lsdtt_flow_base& flow,
    float                  threshold
) const
{
    std::vector<lsdtt_ridgeline_topological_triple_divide> result;
    if (nodes_.size() < 3) return result;
    for (std::size_t i = 1; i + 1 < nodes_.size(); ++i)
    {
        const EdgePixels in_edge = get_edge_pixels(nodes_[i - 1], nodes_[i]);
        const EdgePixels out_edge = get_edge_pixels(nodes_[i], nodes_[i + 1]);
        if (!in_edge.valid || !out_edge.valid) continue;
        auto flat = [&](int r, int c) -> std::int64_t {
            return flow.in_bounds(r, c) ? flow.flat_index_unchecked(r, c) : -1;
        };
        const std::int64_t li = flat(in_edge.left_r, in_edge.left_c);
        const std::int64_t lo = flat(out_edge.left_r, out_edge.left_c);
        const std::int64_t ri = flat(in_edge.right_r, in_edge.right_c);
        const std::int64_t ro = flat(out_edge.right_r, out_edge.right_c);
        lsdtt_ridgeline_topological_triple_divide divide;
        divide.node_index = static_cast<int>(i);
        if (flow.is_valid_cell(li) && flow.is_valid_cell(lo))
        {
            const auto metrics = calculate_flow_topological_distance(li, lo, flow);
            divide.left = {metrics.distance, metrics.junction_flat};
        }
        if (flow.is_valid_cell(ri) && flow.is_valid_cell(ro))
        {
            const auto metrics = calculate_flow_topological_distance(ri, ro, flow);
            divide.right = {metrics.distance, metrics.junction_flat};
        }
        const bool left_qualifies = divide.left.junction_flat >= 0 &&
                                    divide.left.distance > threshold;
        const bool right_qualifies = divide.right.junction_flat >= 0 &&
                                     divide.right.distance > threshold;
        if (left_qualifies || right_qualifies) result.push_back(divide);
    }
    return result;
}

std::vector<lsdtt_ridgeline_topological_flowpath>
lsdtt_ridgeline_cpu::get_topological_distance_flowpaths(
    const lsdtt_flow_base& flow,
    int                    node_index
) const
{
    std::vector<lsdtt_ridgeline_topological_flowpath> result;
    if (node_index <= 0 || static_cast<std::size_t>(node_index + 1) >= nodes_.size())
        return result;

    const EdgePixels in_edge = get_edge_pixels(nodes_[node_index - 1], nodes_[node_index]);
    const EdgePixels out_edge = get_edge_pixels(nodes_[node_index], nodes_[node_index + 1]);
    if (!in_edge.valid || !out_edge.valid) return result;

    auto flat = [&](int r, int c) -> std::int64_t {
        return flow.in_bounds(r, c) ? flow.flat_index_unchecked(r, c) : -1;
    };
    const std::int64_t left_in = flat(in_edge.left_r, in_edge.left_c);
    const std::int64_t left_out = flat(out_edge.left_r, out_edge.left_c);
    const std::int64_t right_in = flat(in_edge.right_r, in_edge.right_c);
    const std::int64_t right_out = flat(out_edge.right_r, out_edge.right_c);

    auto trace_to_junction = [&](std::int64_t start, std::int64_t junction) {
        std::vector<lsdtt_xy> path;
        if (start < 0 || junction < 0 || !flow.is_valid_cell(start)) return path;
        std::int64_t current = start;
        for (std::size_t step = 0; step <= nodes_.size() * 4; ++step)
        {
            path.push_back(cell_center_from_row_col(
                flow.meta(), flow.row_from_flat_unchecked(current),
                flow.col_from_flat_unchecked(current)));
            if (current == junction) return path;
            if (flow.is_outlet_unchecked(current)) break;
            current = flow.receiver_of_unchecked(current);
            if (current < 0 || !flow.is_valid_cell(current)) break;
        }
        path.clear();
        return path;
    };

    auto add_side = [&](const std::string& side, std::int64_t first,
                        std::int64_t second) {
        if (first < 0 || second < 0) return;
        const lsdtt_flow_topology_metrics metrics =
            calculate_flow_topological_distance(first, second, flow);
        if (metrics.junction_flat < 0) return;
        const auto first_path = trace_to_junction(first, metrics.junction_flat);
        const auto second_path = trace_to_junction(second, metrics.junction_flat);
        if (first_path.empty() || second_path.empty()) return;
        result.push_back({side, "incoming", first_path, metrics.distance, metrics.junction_flat});
        result.push_back({side, "outgoing", second_path, metrics.distance, metrics.junction_flat});
    };

    add_side("left", left_in, left_out);
    add_side("right", right_in, right_out);
    return result;
}




// =============================================================================
// CSV serialization utilities
// =============================================================================

std::vector<std::vector<std::string>> lsdtt_ridgeline_cpu::get_ridgeline_csv_rows(
    const std::string&                id_string,
    const std::string&                ridgeline_type,
    const lsdtt_coordinate_transform& transform
) const
{
    std::vector<std::vector<std::string>> rows;
    rows.reserve(nodes_.size());

    for (std::size_t i = 0; i < nodes_.size(); ++i)
    {
        lsdtt_lonlat ll = xy_to_lonlat(positions_[i].x, positions_[i].y, transform);

        std::string fit_z   = (i < fitted_elevations_.size()) ? format_float(fitted_elevations_[i]) : "NaN";
        std::string det_z   = (i < detrended_elevations_.size()) ? format_float(detrended_elevations_[i]) : "NaN";
        std::string s_det_z = (i < smoothed_detrended_elevations_.size()) ? format_float(smoothed_detrended_elevations_[i]) : "NaN";
        std::string s_det_s = (i < smoothed_detrended_slopes_.size()) ? format_float(smoothed_detrended_slopes_[i]) : "NaN";

        std::vector<std::string> row = {
            id_string,
            std::to_string(nodes_[i].row),
            std::to_string(nodes_[i].col),
            format_float(ll.latitude),
            format_float(ll.longitude),
            ridgeline_type,
            format_float(distances_[i]),
            (i < elevations_.size()) ? format_float(elevations_[i]) : "NaN",
            fit_z,
            det_z,
            s_det_z,
            s_det_s
        };
        if (!nearest_channel_elevations_left_.empty() ||
            !nearest_channel_elevations_right_.empty())
        {
            row.push_back((i < nearest_channel_elevations_left_.size())
                ? format_float(nearest_channel_elevations_left_[i]) : "NaN");
            row.push_back((i < nearest_channel_elevations_right_.size())
                ? format_float(nearest_channel_elevations_right_[i]) : "NaN");
        }
        for (std::size_t window_index = 0;
             window_index < moving_window_sizes_pixels_.size(); ++window_index)
        {
            const auto& fitted = moving_window_fitted_elevations_[window_index];
            const auto& detrended = moving_window_detrended_elevations_[window_index];
            row.push_back(i < fitted.size() ? format_float(fitted[i]) : "NaN");
            row.push_back(i < detrended.size() ? format_float(detrended[i]) : "NaN");
        }
        rows.push_back(std::move(row));
    }

    return rows;
}

std::vector<std::string> lsdtt_ridgeline_cpu::get_stats_csv_row_single_basin(
    int                basin_junction_number,
    const std::string& ridgeline_type,
    float              basin_drainage_area,
    float              concavity
) const
{
    return {
        std::to_string(basin_junction_number),
        ridgeline_type,
        format_float(basin_drainage_area),
        format_float(slope_fit_),
        format_float(r2_fit_),
        format_float(calculate_ridgeline_steepness(basin_drainage_area, concavity)),
        format_float(poly2_r2_),
        format_float(poly2_A_),
        format_float(poly2_B_),
        format_float(poly2_C_)
    };
}

std::vector<std::string> lsdtt_ridgeline_cpu::get_stats_csv_row_dual_basin(
    int                ridgeline_number,
    const std::string& ridgeline_type,
    int                left_basin_junction_number,
    int                right_basin_junction_number,
    float              left_basin_drainage_area,
    float              right_basin_drainage_area,
    float              concavity
) const
{
    float NoDataValue = divide_grid_meta_.no_data_value;
    std::string s_left  = left_basin_junction_number >= 0 ? std::to_string(left_basin_junction_number) : format_float(NoDataValue);
    std::string s_right = right_basin_junction_number >= 0 ? std::to_string(right_basin_junction_number) : format_float(NoDataValue);

    float stp_left  = calculate_ridgeline_steepness(left_basin_drainage_area, concavity);
    float stp_right = calculate_ridgeline_steepness(right_basin_drainage_area, concavity);

    return {
        std::to_string(ridgeline_number),
        ridgeline_type,
        s_left,
        s_right,
        (left_basin_drainage_area >= 0.0) ? format_float(left_basin_drainage_area) : format_float(NoDataValue),
        (right_basin_drainage_area >= 0.0) ? format_float(right_basin_drainage_area) : format_float(NoDataValue),
        format_float(slope_fit_),
        format_float(r2_fit_),
        (left_basin_drainage_area >= 0.0) ? format_float(stp_left) : format_float(NoDataValue),
        (right_basin_drainage_area >= 0.0) ? format_float(stp_right) : format_float(NoDataValue),
        format_float(poly2_r2_),
        format_float(poly2_A_),
        format_float(poly2_B_),
        format_float(poly2_C_)
    };
}

std::vector<std::vector<std::string>> lsdtt_ridgeline_cpu::get_extrema_csv_rows(
    const std::string&                           id_string,
    const std::string&                           ridgeline_type,
    const std::vector<lsdtt_ridgeline_extremum>& extrema,
    const lsdtt_coordinate_transform&            transform
) const
{
    std::vector<std::vector<std::string>> rows;
    rows.reserve(extrema.size());

    for (const auto& ext : extrema)
    {
        lsdtt_lonlat ll = xy_to_lonlat(positions_[ext.node_index].x, positions_[ext.node_index].y, transform);
        
        std::vector<std::string> row = {
            id_string,
            ridgeline_type,
            format_float(ll.latitude),
            format_float(ll.longitude),
            format_float(distances_[ext.node_index]),
            ext.type,
            format_float(ext.distance_to_next),
            format_float(ext.relief_to_next),
            format_float(ext.distance_to_tdj)
        };
        if (!nearest_channel_elevations_left_.empty() ||
            !nearest_channel_elevations_right_.empty())
        {
            row.push_back(ext.node_index < static_cast<int>(nearest_channel_elevations_left_.size())
                ? format_float(nearest_channel_elevations_left_[ext.node_index]) : "NaN");
            row.push_back(ext.node_index < static_cast<int>(nearest_channel_elevations_right_.size())
                ? format_float(nearest_channel_elevations_right_[ext.node_index]) : "NaN");
        }
        rows.push_back(std::move(row));
    }

    return rows;
}

std::vector<std::vector<std::string>> lsdtt_ridgeline_cpu::get_flow_topological_distance_csv_rows(
    const std::string&                id_string,
    const std::string&                ridgeline_type,
    const lsdtt_flow_base&            flow,
    const lsdtt_raster_cpu&           dem,
    const lsdtt_coordinate_transform& transform
) const
{
    std::vector<std::vector<std::string>> rows;
    if (nodes_.empty()) return rows; 

    rows.reserve(nodes_.size());

    const auto& data = dem.get_data_array();
    const auto& mask = dem.get_mask_array();
    const float ndv  = dem.get_meta().no_data_value;

    for (std::size_t i = 0; i < nodes_.size(); ++i) {
        
        lsdtt_flow_topology_metrics left_met  = {-9999.0f, -1};
        lsdtt_flow_topology_metrics right_met = {-9999.0f, -1};
        
        // Exclude start and end node from getting distances (no fully formed incoming/outgoing edge limits)
        if (i > 0 && i < nodes_.size() - 1) {
            EdgePixels in_edge = get_edge_pixels(nodes_[i-1], nodes_[i]);
            EdgePixels out_edge = get_edge_pixels(nodes_[i], nodes_[i+1]);

            if (in_edge.valid && out_edge.valid) {
                std::int64_t left_in_flat = flow.in_bounds(in_edge.left_r, in_edge.left_c) ? 
                    flow.flat_index_unchecked(in_edge.left_r, in_edge.left_c) : -1;
                std::int64_t left_out_flat = flow.in_bounds(out_edge.left_r, out_edge.left_c) ? 
                    flow.flat_index_unchecked(out_edge.left_r, out_edge.left_c) : -1;

                std::int64_t right_in_flat = flow.in_bounds(in_edge.right_r, in_edge.right_c) ? 
                    flow.flat_index_unchecked(in_edge.right_r, in_edge.right_c) : -1;
                std::int64_t right_out_flat = flow.in_bounds(out_edge.right_r, out_edge.right_c) ? 
                    flow.flat_index_unchecked(out_edge.right_r, out_edge.right_c) : -1;

                if (left_in_flat >= 0 && left_out_flat >= 0) {
                    left_met = calculate_flow_topological_distance(left_in_flat, left_out_flat, flow);
                }
                if (right_in_flat >= 0 && right_out_flat >= 0) {
                    right_met = calculate_flow_topological_distance(right_in_flat, right_out_flat, flow);
                }
            }
        }

        // Establish the main node metrics
        lsdtt_lonlat ll = xy_to_lonlat(positions_[i].x, positions_[i].y, transform);
        float ridge_z = (i < elevations_.size() && elevations_[i] != ndv) ? elevations_[i] : -9999.0f;

        // Perform the elevation lookups and differences
        float left_junc_z = -9999.0f;
        float right_junc_z = -9999.0f;
        float left_diff_z = -9999.0f;
        float right_diff_z = -9999.0f;

        if (left_met.junction_flat >= 0) {
            int r = flow.row_from_flat_unchecked(left_met.junction_flat);
            int c = flow.col_from_flat_unchecked(left_met.junction_flat);
            if (mask(r,c) != 0 && data(r,c) != ndv) {
                left_junc_z = data(r,c);
                if (ridge_z != -9999.0f) left_diff_z = ridge_z - left_junc_z;
            }
        }
        
        if (right_met.junction_flat >= 0) {
            int r = flow.row_from_flat_unchecked(right_met.junction_flat);
            int c = flow.col_from_flat_unchecked(right_met.junction_flat);
            if (mask(r,c) != 0 && data(r,c) != ndv) {
                right_junc_z = data(r,c);
                if (ridge_z != -9999.0f) right_diff_z = ridge_z - right_junc_z;
            }
        }

        std::vector<std::string> row = {
            id_string,
            ridgeline_type,
            format_float(ll.latitude),
            format_float(ll.longitude),
            format_float(distances_[i]),
            format_float(ridge_z),
            format_float(left_met.distance),
            format_float(right_met.distance),
            format_float(left_junc_z),
            format_float(right_junc_z),
            format_float(left_diff_z),
            format_float(right_diff_z)
        };
        if (!nearest_channel_elevations_left_.empty() ||
            !nearest_channel_elevations_right_.empty())
        {
            row.push_back(i < nearest_channel_elevations_left_.size()
                ? format_float(nearest_channel_elevations_left_[i]) : "NaN");
            row.push_back(i < nearest_channel_elevations_right_.size()
                ? format_float(nearest_channel_elevations_right_[i]) : "NaN");
        }
        rows.push_back(std::move(row));
    }

    return rows;
}




} // namespace lsdtt
