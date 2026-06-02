// lsdtt_ridgeline_cpu.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Rutger Vervoordeldonk and Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
//
// See lsdtt_ridgeline_cpu.hpp for full documentation.

#include "lsdtt_ridgeline_cpu.hpp"
#include "lsdtt_raster_cpu.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <unordered_set>

namespace lsdtt
{

// =============================================================================
// Construction
// =============================================================================

lsdtt_ridgeline_cpu::lsdtt_ridgeline_cpu(
    std::vector<lsdtt_pixel_rc>  nodes,
    const lsdtt_raster_meta&     divide_grid_meta)
    : divide_grid_meta_(divide_grid_meta)
    , nodes_(std::move(nodes))
{
    compute_positions_and_distances_();
}

// =============================================================================
// Private helpers
// =============================================================================

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

// =============================================================================
// Elevation filling
// =============================================================================

void lsdtt_ridgeline_cpu::fill_elevation(const lsdtt_raster_cpu& dem)
{
    elevations_.clear();
    elevations_.reserve(nodes_.size());

    const lsdtt_raster_meta& dem_meta = dem.get_meta();
    const auto& data  = dem.get_data_array();
    const auto& mask  = dem.get_mask_array();
    const float ndv   = dem_meta.no_data_value;

    const int dem_nrows = dem_meta.n_rows;
    const int dem_ncols = dem_meta.n_cols;

    for (const lsdtt_pixel_rc& dg_rc : nodes_)
    {
        // The divide-grid is shifted by half a cell relative to the DEM.
        // A divide-grid node (dg_row, dg_col) sits at the corner between the
        // four DEM pixels:
        //   (dg_row - 1, dg_col - 1)   (dg_row - 1, dg_col)
        //   (dg_row,     dg_col - 1)   (dg_row,     dg_col)
        // (in non-periodic DEM-pixel coordinates).
        //
        // We average all valid surrounding pixels; fall back to the nearest
        // valid pixel if none of the four are valid (e.g. at DEM edges).

        float sum_elev  = 0.0f;
        int   n_valid   = 0;

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

        if (n_valid > 0)
        {
            elevations_.push_back(sum_elev / static_cast<float>(n_valid));
        }
        else
        {
            elevations_.push_back(ndv);
        }
    }
}

// =============================================================================
// Mutating operations
// =============================================================================

void lsdtt_ridgeline_cpu::reverse()
{
    std::reverse(nodes_.begin(),     nodes_.end());
    std::reverse(positions_.begin(), positions_.end());
    std::reverse(elevations_.begin(), elevations_.end());

    // Recompute distances (they are cumulative from the new start).
    compute_positions_and_distances_();
}

void lsdtt_ridgeline_cpu::append_without_check(const lsdtt_ridgeline_cpu& other)
{
    if (other.empty()) return;

    // Append nodes.
    nodes_.insert(nodes_.end(), other.nodes_.begin(), other.nodes_.end());

    // Append elevations if both have them; otherwise clear.
    if (!elevations_.empty() && !other.elevations_.empty())
    {
        elevations_.insert(elevations_.end(),
                           other.elevations_.begin(),
                           other.elevations_.end());
    }
    else
    {
        elevations_.clear();
    }

    // Recompute positions and distances from the updated node sequence.
    compute_positions_and_distances_();
}

// =============================================================================
// MDD assembly helper
// =============================================================================

std::vector<lsdtt_ridgeline_cpu>
lsdtt_ridgeline_cpu::extract_contiguous_shared_ridgelines(
    const lsdtt_ridgeline_cpu& other) const
{
    if (empty() || other.empty()) return {};

    // Build a lookup set of nodes in *other for O(1) membership queries.
    // We use a flat hash of (row * large_prime + col) as key.
    // Collisions are possible but exceedingly unlikely for realistic DEM sizes.
    // A full map<pair> would be cleaner but slower.
    auto node_hash = [](const lsdtt_pixel_rc& rc) -> std::size_t
    {
        // Cantor pairing — unique for non-negative integers.
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

    // Walk *this and split into contiguous shared runs.
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
// Placeholder stubs
// =============================================================================

void lsdtt_ridgeline_cpu::segment_ridgeline_snout_by_curvature(
    float /*smoothing_window_fraction*/)
{
    // TODO: port Robinson et al. (2025) / Satopaa et al. (2011) kneed algorithm.
    // This requires LOWESS smoothing followed by maximum-curvature detection.
    // For now this is a no-op so that the rest of the pipeline can run.
}

void lsdtt_ridgeline_cpu::fit_ridgeline_profile(
    bool /*use_polynomial*/,
    bool /*use_constant*/)
{
    // TODO: fit linear (and optionally quadratic) trend to
    // elevation vs. distance along ridgeline.
}

void lsdtt_ridgeline_cpu::smoothen_and_derive_detrended_ridgeline(
    int /*smoothing_window_pixels*/)
{
    // TODO: subtract fitted trend, apply moving-average smoothing, and compute
    // first-difference slope of the smoothed detrended profile.
}

} // namespace lsdtt