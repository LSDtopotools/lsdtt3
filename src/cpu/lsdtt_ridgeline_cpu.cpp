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
