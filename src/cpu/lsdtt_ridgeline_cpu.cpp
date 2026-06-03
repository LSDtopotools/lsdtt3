// lsdtt_ridgeline_cpu.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Rutger Vervoordeldonk and Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html

#include "lsdtt_ridgeline_cpu.hpp"
#include "lsdtt_flowinfo_cpu.hpp" // For dynamic_cast accessing CPU-specific config
#include "../utils/lsdtt_regression.hpp"
#include "../utils/lsdtt_signal1d.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <stdexcept>
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
        const double dx = positions_[i].x - positions_[i - 1].x;
        const double dy = positions_[i].y - positions_[i - 1].y;
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

    compute_positions_and_distances_();
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

void lsdtt_ridgeline_cpu::fit_ridgeline_profile(bool use_polynomial, bool use_constant)
{
    if (static_cast<std::size_t>(snoutidx_) >= nodes_.size()) return;

    std::vector<double> dist_main(distances_.begin() + snoutidx_, distances_.end());
    std::vector<double> elev_main(elevations_.begin() + snoutidx_, elevations_.end());

    fitted_elevations_.assign(nodes_.size(), divide_grid_meta_.no_data_value);
    detrended_elevations_.assign(nodes_.size(), divide_grid_meta_.no_data_value);

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
    if (nodes_[snoutidx_].row == nodes_.back().row && nodes_[snoutidx_].col == nodes_.back().col)
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
    float dist_to_tdj = divide_grid_meta_.no_data_value;

    if (node_index < 0 || static_cast<std::size_t>(node_index) >= nodes_.size()) return dist_to_tdj;

    std::string bc = "oooo"; 
    if (const auto* fi = dynamic_cast<const lsdtt_flowinfo*>(&flow))
    {
        bc = fi->boundary_conditions();
    }

    bool periodic_ns = (bc.size() == 4 && std::tolower(static_cast<unsigned char>(bc[0])) == 'p' && std::tolower(static_cast<unsigned char>(bc[2])) == 'p');
    bool periodic_ew = (bc.size() == 4 && std::tolower(static_cast<unsigned char>(bc[1])) == 'p' && std::tolower(static_cast<unsigned char>(bc[3])) == 'p');
    
    int n_rows = flow.n_rows();
    int n_cols = flow.n_cols();

    auto wrap_idx = [](int idx, int max_size, bool is_periodic) -> int {
        if (is_periodic) return (idx % max_size + max_size) % max_size;
        return (idx >= 0 && idx < max_size) ? idx : -1;
    };

    auto check_node_for_tdj = [&](int current_index) -> bool 
    {
        int divrow = nodes_[current_index].row;
        int divcol = nodes_[current_index].col;

        int raw_xs[2] = {divcol - 1, divcol};
        int raw_ys[2] = {divrow - 1, divrow};

        std::set<std::int32_t> basin_junctions;

        for (int x : raw_xs)
        {
            for (int y : raw_ys)
            {
                int r = wrap_idx(y, n_rows, periodic_ns);
                int c = wrap_idx(x, n_cols, periodic_ew);

                if (r != -1 && c != -1 && flow.is_valid_cell(r, c))
                {
                    std::int64_t flat = flow.flat_index_unchecked(r, c);
                    std::int64_t curr = flat;
                    
                    while (curr != lsdtt_flow_base::nodata_receiver && curr != lsdtt_flow_base::outlet_receiver)
                    {
                        if (junction_network.is_channel(curr))
                        {
                            std::int32_t j_idx = junction_network.junction_of_flat(curr);
                            if (j_idx >= 0) basin_junctions.insert(j_idx);
                            break;
                        }
                        curr = flow.receiver_of_unchecked(curr);
                    }
                }
            }
        }

        if (basin_junctions.size() > 2)
        {
            float current_tdj_dist = static_cast<float>(std::abs(distances_[current_index] - distances_[node_index]));
            if (std::isnan(dist_to_tdj) || current_tdj_dist < dist_to_tdj)
            {
                dist_to_tdj = current_tdj_dist;
            }
            return true;
        }
        return false;
    };

    for (std::size_t i = node_index; i < nodes_.size(); ++i)
    {
        if (check_node_for_tdj(static_cast<int>(i))) break;
    }

    for (int i = node_index; i >= 0; --i)
    {
        if (check_node_for_tdj(i)) break;
    }

    return dist_to_tdj;
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

        rows.push_back({
            id_string,
            format_float(ll.latitude),
            format_float(ll.longitude),
            ridgeline_type,
            format_float(distances_[i]),
            (i < elevations_.size()) ? format_float(elevations_[i]) : "NaN",
            fit_z,
            det_z,
            s_det_z,
            s_det_s
        });
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
        
        rows.push_back({
            id_string,
            ridgeline_type,
            format_float(ll.latitude),
            format_float(ll.longitude),
            format_float(distances_[ext.node_index]),
            ext.type,
            format_float(ext.distance_to_next),
            format_float(ext.relief_to_next),
            format_float(ext.distance_to_tdj)
        });
    }

    return rows;
}


} // namespace lsdtt
