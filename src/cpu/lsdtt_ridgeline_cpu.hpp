// lsdtt_ridgeline_cpu.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Rutger Vervoordeldonk and Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html
//
// -----------------------------------------------------------------------------
// lsdtt_ridgeline_cpu
//
// Holds a sequence of divide-grid nodes that define a ridgeline (a segment of
// a drainage divide). Ridgelines are produced by the basin-outline splitting
// and divide-extraction methods on lsdtt_basin_outline_cpu and are consumed by
// the ridgeline analysis routines in lsdtt-ridgeline-analysis.
// -----------------------------------------------------------------------------

#pragma once

#include <cstdint>
#include <string>
#include <vector>
#include <optional>
#include <array>

#include "../core/lsdtt_flow_base.hpp"
#include "../core/lsdtt_junction_network_base.hpp"
#include "lsdtt_raster_cpu.hpp"
#include "../core/lsdtt_raster_meta.hpp"
#include "../core/lsdtt_types.hpp"
#include "../utils/lsdtt_raster_geometry.hpp"
#include "../core/lsdtt_vector_utils.hpp"

namespace lsdtt
{

// Structure to hold information about a ridgeline extremum (min/max)
struct lsdtt_ridgeline_extremum
{
    int         node_index;
    std::string type;
    double      distance_to_next;
    double      relief_to_next;
    double      distance_to_tdj;
};

struct lsdtt_ridgeline_topological_flowpath
{
    std::string side;
    std::string branch;
    std::vector<lsdtt_xy> positions;
    float pair_distance = -9999.0f;
    std::int64_t junction_flat = -1;
};

// =============================================================================
// lsdtt_ridgeline_cpu
// =============================================================================
class lsdtt_ridgeline_cpu
{
public:
    // -------------------------------------------------------------------------
    // Construction
    // -------------------------------------------------------------------------
    lsdtt_ridgeline_cpu() = default;

    lsdtt_ridgeline_cpu(
        std::vector<lsdtt_pixel_rc> nodes,
        const lsdtt_raster_meta&    divide_grid_meta
    );

    // -------------------------------------------------------------------------
    // Basic queries
    // -------------------------------------------------------------------------
    [[nodiscard]] bool        empty() const noexcept { return nodes_.empty(); }
    [[nodiscard]] std::size_t size() const noexcept { return nodes_.size(); }

    [[nodiscard]] const std::vector<lsdtt_pixel_rc>& nodes() const noexcept { return nodes_; }
    [[nodiscard]] const std::vector<lsdtt_xy>&       positions() const noexcept { return positions_; }
    [[nodiscard]] const std::vector<double>&         distances() const noexcept { return distances_; }
    [[nodiscard]] const std::vector<float>&          elevations() const noexcept { return elevations_; }

    [[nodiscard]] double total_length() const noexcept
    {
        return distances_.empty() ? 0.0 : distances_.back();
    }

    [[nodiscard]] const lsdtt_raster_meta& divide_grid_meta() const noexcept
    {
        return divide_grid_meta_;
    }
    
    // Fit Statistics
    [[nodiscard]] int   snout_index() const noexcept { return snoutidx_; }
    [[nodiscard]] float slope_fit() const noexcept { return slope_fit_; }
    [[nodiscard]] float r_squared_fit() const noexcept { return r2_fit_; }
    [[nodiscard]] float polynomial_A() const noexcept { return poly2_A_; }
    [[nodiscard]] float polynomial_B() const noexcept { return poly2_B_; }
    [[nodiscard]] float polynomial_C() const noexcept { return poly2_C_; }
    [[nodiscard]] float polynomial_r_squared() const noexcept { return poly2_r2_; }

    // Check if the current ridgeline nodes are 4-connected
    [[nodiscard]] bool check_if_nodes_four_connected(bool periodic_ns, bool periodic_ew) const;

    // -------------------------------------------------------------------------
    // Elevation filling
    // -------------------------------------------------------------------------
    
    // Sample DEM elevation at each divide-grid node and store in elevations_.
    void fill_elevation(const lsdtt_raster_cpu& dem);

    // Sample continuous data attributes from any input raster over divide-nodes.
    [[nodiscard]] std::vector<float> sample_raster_values(const lsdtt_raster_cpu& raster) const;

    // -------------------------------------------------------------------------
    // Mutating operations
    // -------------------------------------------------------------------------

    void reverse();
    void append_without_check(const lsdtt_ridgeline_cpu& other);

    // -------------------------------------------------------------------------
    // MDD assembly helper
    // -------------------------------------------------------------------------

    [[nodiscard]] std::vector<lsdtt_ridgeline_cpu> extract_contiguous_shared_ridgelines(
        const lsdtt_ridgeline_cpu& other
    ) const;

    // -------------------------------------------------------------------------
    // Ridgeline analysis tools
    // -------------------------------------------------------------------------

    // Segment the ridge snout off by max curvature following Robinson et al. (2025).
    void segment_ridgeline_snout_by_curvature(float smoothing_window_fraction);

    // Fit the ridgeline elevation profile with a linear and optionally polynomial trend.
    // Skips points in the segmented snout (if identified).
    void fit_ridgeline_profile(bool use_polynomial, bool use_constant);

    // Compute the detrended elevation profile and smooth it.
    void smoothen_and_derive_detrended_ridgeline(int smoothing_window_pixels);

    // Calculate normalized ridgeline steepness for basin-outline based ridgelines.
    [[nodiscard]] float calculate_ridgeline_steepness(
        float basin_drainage_area, 
        float concavity
    ) const;

    // Find local maxima and minima along the smoothened, detrended ridgeline.
    [[nodiscard]] std::vector<lsdtt_ridgeline_extremum> get_ridgeline_extrema(
        const lsdtt_flow_base&             flow,
        const lsdtt_junction_network_base& junction_network
    ) const;

    // Calculate nearest distance along the ridgeline to a Triple/Quadruple Divide Junction.
    [[nodiscard]] float find_distance_to_closest_tdj(
        int                                node_index,
        const lsdtt_flow_base&             flow,
        const lsdtt_junction_network_base& junction_network
    ) const;

    // -------------------------------------------------------------------------
    // CSV output serialization routines
    // -------------------------------------------------------------------------
    
    [[nodiscard]] std::vector<std::vector<std::string>> get_ridgeline_csv_rows(
        const std::string&                id_string,
        const std::string&                ridgeline_type,
        const lsdtt_coordinate_transform& transform
    ) const;

    [[nodiscard]] std::vector<std::string> get_stats_csv_row_single_basin(
        int                basin_junction_number,
        const std::string& ridgeline_type,
        float              basin_drainage_area,
        float              concavity
    ) const;

    [[nodiscard]] std::vector<std::string> get_stats_csv_row_dual_basin(
        int                ridgeline_number,
        const std::string& ridgeline_type,
        int                left_basin_junction_number,
        int                right_basin_junction_number,
        float              left_basin_drainage_area,
        float              right_basin_drainage_area,
        float              concavity
    ) const;

    [[nodiscard]] std::vector<std::vector<std::string>> get_extrema_csv_rows(
        const std::string&                           id_string,
        const std::string&                           ridgeline_type,
        const std::vector<lsdtt_ridgeline_extremum>& extrema,
        const lsdtt_coordinate_transform&            transform
    ) const;

    [[nodiscard]] std::vector<std::vector<std::string>> get_flow_topological_distance_csv_rows(
        const std::string&                id_string,
        const std::string&                ridgeline_type,
        const lsdtt_flow_base&            flow,
        const lsdtt_raster_cpu&           dem,
        const lsdtt_coordinate_transform& transform
    ) const;

    // Return internal divide nodes that are triple/quadruple junctions, either
    // by the surrounding-basin test or by the left/right flow-topological
    // distance threshold.
    [[nodiscard]] std::vector<int> get_triple_divide_node_indices(
        const lsdtt_flow_base&             flow,
        const lsdtt_junction_network_base& junction_network
    ) const;

    [[nodiscard]] std::vector<int> get_topological_distance_triple_divide_node_indices(
        const lsdtt_flow_base& flow,
        float                  threshold
    ) const;

    [[nodiscard]] std::vector<lsdtt_ridgeline_topological_flowpath>
    get_topological_distance_flowpaths(
        const lsdtt_flow_base& flow,
        int                    node_index
    ) const;




private:
    void compute_positions_and_distances_();

    lsdtt_raster_meta           divide_grid_meta_;
    std::vector<lsdtt_pixel_rc> nodes_;
    std::vector<lsdtt_xy>       positions_;
    std::vector<double>         distances_;
    
    // Core property vectors
    std::vector<float> elevations_;
    std::vector<float> fitted_elevations_;
    std::vector<float> detrended_elevations_;
    std::vector<float> smoothed_detrended_elevations_;
    std::vector<float> smoothed_detrended_slopes_;

    // Core state properties
    int   snoutidx_      = 0;
    float slope_fit_     = std::numeric_limits<float>::quiet_NaN();
    float intercept_fit_ = std::numeric_limits<float>::quiet_NaN();
    float r2_fit_        = std::numeric_limits<float>::quiet_NaN();
    
    float poly2_A_       = std::numeric_limits<float>::quiet_NaN();
    float poly2_B_       = std::numeric_limits<float>::quiet_NaN();
    float poly2_C_       = std::numeric_limits<float>::quiet_NaN();
    float poly2_r2_      = std::numeric_limits<float>::quiet_NaN();
};

} // namespace lsdtt
