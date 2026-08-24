// lsdtt_hilltop_flow_routing_cpu.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
//
// See lsdtt_hilltop_flow_routing_cpu.hpp for the divide-grid convention and the
// side-assignment rule.

#include "lsdtt_hilltop_flow_routing_cpu.hpp"

#include <algorithm>
#include <stdexcept>
#include <string>
#include <unordered_set>

namespace lsdtt
{

namespace
{

// Signed one-step difference on an axis that may wrap.
// Returns 0 if the two indices are equal, +1/-1 for a single step (through the
// seam when periodic), and something with |d| > 1 otherwise.
int step_delta(int a, int b, int n, bool periodic)
{
    int d = b - a;
    if (periodic && n > 0)
    {
        if (d >  n / 2) d -= n;
        if (d < -n / 2) d += n;
    }
    return d;
}

int wrap_index(int idx, int n, bool periodic)
{
    if (periodic && n > 0) return ((idx % n) + n) % n;
    return idx;   // may be out of range; the caller checks
}

} // namespace

// =============================================================================
// infer_divide_grid_periodicity
// =============================================================================
void infer_divide_grid_periodicity(
    const lsdtt_raster_meta& dem_meta,
    const lsdtt_raster_meta& divide_grid_meta,
    bool& periodic_ns,
    bool& periodic_ew)
{
    const int dr = divide_grid_meta.n_rows - dem_meta.n_rows;
    const int dc = divide_grid_meta.n_cols - dem_meta.n_cols;

    if ((dr != 0 && dr != 1) || (dc != 0 && dc != 1))
    {
        throw std::invalid_argument(
            "lsdtt_hilltop_flow_routing_cpu: the divide grid must be n or n+1 "
            "on each axis relative to the DEM (got divide "
            + std::to_string(divide_grid_meta.n_rows) + "x"
            + std::to_string(divide_grid_meta.n_cols) + " against DEM "
            + std::to_string(dem_meta.n_rows) + "x"
            + std::to_string(dem_meta.n_cols)
            + "). These two metas are not a consistent pair.");
    }

    // n+1 corners on an axis means the axis is NOT periodic.
    periodic_ns = (dr == 0);
    periodic_ew = (dc == 0);
}

// =============================================================================
// classify_ridgeline_step
// =============================================================================
ridgeline_step_kind classify_ridgeline_step(
    const lsdtt_pixel_rc&    a,
    const lsdtt_pixel_rc&    b,
    const lsdtt_raster_meta& divide_grid_meta,
    bool                     periodic_ns,
    bool                     periodic_ew)
{
    const int dr = step_delta(a.row, b.row, divide_grid_meta.n_rows, periodic_ns);
    const int dc = step_delta(a.col, b.col, divide_grid_meta.n_cols, periodic_ew);

    const int adr = dr < 0 ? -dr : dr;
    const int adc = dc < 0 ? -dc : dc;

    if (adr == 0 && adc == 0) return ridgeline_step_kind::repeated;
    if (adr > 1 || adc > 1)   return ridgeline_step_kind::disjoint;
    if (adr == 1 && adc == 1) return ridgeline_step_kind::diagonal;
    return ridgeline_step_kind::cardinal;
}

// =============================================================================
// filter_hilltop_pixels
// =============================================================================
std::vector<lsdtt_hilltop_pixel> filter_hilltop_pixels(
    const std::vector<lsdtt_hilltop_pixel>& pixels,
    const lsdtt_raster_meta&                dem_meta,
    const lsdtt_hilltop_filter_options&     opts,
    const std::vector<float>*               channel_mask,
    const std::vector<float>*               gradient,
    const std::vector<float>*               curvature,
    lsdtt_hilltop_filter_report&            report)
{
    const int nr = dem_meta.n_rows;
    const int nc = dem_meta.n_cols;
    const std::size_t ncell = static_cast<std::size_t>(nr) * nc;

    auto check_size = [&](const std::vector<float>* v, const char* what) {
        if (v && v->size() != ncell)
            throw std::invalid_argument(
                std::string("lsdtt_hilltop_flow_routing_cpu: ") + what +
                " raster does not match the DEM grid (" +
                std::to_string(v->size()) + " against " +
                std::to_string(ncell) + " cells).");
    };
    check_size(channel_mask, "channel mask");
    check_size(gradient, "gradient");
    check_size(curvature, "curvature");

    report = lsdtt_hilltop_filter_report{};
    report.n_in = static_cast<std::int64_t>(pixels.size());

    if (channel_mask)
        for (const auto& p : pixels)
            if ((*channel_mask)[static_cast<std::size_t>(p.flat)] > 0.5f)
                ++report.n_on_channel_before;

    const int pad = opts.stream_network_padding_pixels;

    std::vector<lsdtt_hilltop_pixel> out;
    out.reserve(pixels.size());

    for (const auto& p : pixels)
    {
        const std::size_t k = static_cast<std::size_t>(p.flat);

        if (channel_mask && pad > 0)
        {
            // Chebyshev neighbourhood: pad = 1 tests the pixel itself only.
            const int reach = pad - 1;
            bool near = false;
            for (int dr = -reach; dr <= reach && !near; ++dr)
                for (int dc = -reach; dc <= reach && !near; ++dc)
                {
                    const int rr = p.row + dr, cc = p.col + dc;
                    if (rr < 0 || rr >= nr || cc < 0 || cc >= nc) continue;
                    if ((*channel_mask)[static_cast<std::size_t>(rr) * nc + cc] > 0.5f)
                        near = true;
                }
            if (near) { ++report.n_dropped_near_channel; continue; }
        }

        if (gradient && opts.threshold_hilltop_gradient > 0.0f &&
            (*gradient)[k] > opts.threshold_hilltop_gradient)
        { ++report.n_dropped_steep; continue; }

        if (curvature && opts.remove_positive_curvature_hilltops &&
            (*curvature)[k] > 0.0f)
        { ++report.n_dropped_convex_up; continue; }

        out.push_back(p);
    }

    report.n_out = static_cast<std::int64_t>(out.size());
    return out;
}

// =============================================================================
// extract_hilltop_pixels
// =============================================================================
lsdtt_ridge_extraction_result extract_hilltop_pixels(
    const std::vector<lsdtt_ridgeline_cpu>& ridgelines,
    const lsdtt_raster_meta&                dem_meta,
    const lsdtt_ridge_extraction_options&   opts,
    const std::vector<std::int64_t>&        basin_outlet_flats)
{
    lsdtt_ridge_extraction_result result;

    if (!basin_outlet_flats.empty() &&
        basin_outlet_flats.size() != ridgelines.size())
    {
        throw std::invalid_argument(
            "lsdtt_hilltop_flow_routing_cpu: basin_outlet_flats must be empty "
            "or the same length as ridgelines (got "
            + std::to_string(basin_outlet_flats.size()) + " against "
            + std::to_string(ridgelines.size()) + ").");
    }

    const int nr = dem_meta.n_rows;
    const int nc = dem_meta.n_cols;

    bool periodic_ns = false, periodic_ew = false;
    bool have_periodicity = false;

    std::unordered_set<std::int64_t> seen;

    // Emit one hilltop pixel for a given DEM row/col, honouring wrap and the
    // deduplication policy.
    auto emit = [&](int row, int col,
                    std::int32_t rid, std::int32_t node_index,
                    double dist_m, double total_m,
                    std::int8_t side, std::int64_t basin_flat)
    {
        const int r = wrap_index(row, nr, periodic_ns);
        const int c = wrap_index(col, nc, periodic_ew);

        if (r < 0 || r >= nr || c < 0 || c >= nc)
        {
            ++result.n_off_grid_dropped;
            return;
        }

        const std::int64_t flat = static_cast<std::int64_t>(r) * nc + c;

        if (opts.deduplicate_pixels && !seen.insert(flat).second)
        {
            ++result.n_duplicates_dropped;
            return;
        }

        lsdtt_hilltop_pixel px;
        px.flat                       = flat;
        px.row                        = r;
        px.col                        = c;
        px.ridgeline_id               = rid;
        px.ridgeline_node_index       = node_index;
        px.distance_along_ridgeline_m = dist_m;
        px.ridgeline_total_length_m   = total_m;
        px.ridge_side                 = side;
        px.basin_outlet_flat          = basin_flat;

        result.pixels.push_back(px);
    };

    // One cardinal edge: emit the two DEM pixels it separates.
    //
    // Travelling EAST  along corner row r: left  = (r-1, c_lo), right = (r, c_lo)
    // Travelling WEST                    : left  = (r,   c_lo), right = (r-1, c_lo)
    // Travelling SOUTH along corner col c: left  = (r_lo, c),   right = (r_lo, c-1)
    // Travelling NORTH                   : left  = (r_lo, c-1), right = (r_lo, c)
    //
    // where c_lo / r_lo is the lower of the two varying corner indices, taken
    // through the seam when the step wraps.
    auto emit_edge = [&](const lsdtt_pixel_rc& a, int dr, int dc,
                         std::int32_t rid, std::int32_t node_index,
                         double dist_m, double total_m,
                         std::int64_t basin_flat)
    {
        if (dr == 0)
        {
            // Horizontal edge. c_lo is the DEM column of the pixel pair: the
            // western of the two corners. Left unwrapped here -- emit() applies
            // the wrap (or drops it off-grid) against the DEM dimensions.
            const int c_lo = (dc > 0) ? a.col : a.col - 1;
            const int r = a.row;

            const int north_row = r - 1;
            const int south_row = r;

            if (dc > 0)   // travelling east: left is north
            {
                emit(north_row, c_lo, rid, node_index, dist_m, total_m, 0, basin_flat);
                emit(south_row, c_lo, rid, node_index, dist_m, total_m, 1, basin_flat);
            }
            else          // travelling west: left is south
            {
                emit(south_row, c_lo, rid, node_index, dist_m, total_m, 0, basin_flat);
                emit(north_row, c_lo, rid, node_index, dist_m, total_m, 1, basin_flat);
            }
        }
        else
        {
            // Vertical edge. r_lo is the DEM row of the pixel pair: the
            // northern of the two corners. Unwrapped; emit() handles it.
            const int r_lo = (dr > 0) ? a.row : a.row - 1;
            const int c = a.col;

            const int west_col = c - 1;
            const int east_col = c;

            if (dr > 0)   // travelling south: left is east
            {
                emit(r_lo, east_col, rid, node_index, dist_m, total_m, 0, basin_flat);
                emit(r_lo, west_col, rid, node_index, dist_m, total_m, 1, basin_flat);
            }
            else          // travelling north: left is west
            {
                emit(r_lo, west_col, rid, node_index, dist_m, total_m, 0, basin_flat);
                emit(r_lo, east_col, rid, node_index, dist_m, total_m, 1, basin_flat);
            }
        }
    };

    for (std::size_t ri = 0; ri < ridgelines.size(); ++ri)
    {
        const lsdtt_ridgeline_cpu& rl = ridgelines[ri];
        if (rl.size() < 2) continue;

        const lsdtt_raster_meta& dgm = rl.divide_grid_meta();

        if (!have_periodicity)
        {
            infer_divide_grid_periodicity(dem_meta, dgm, periodic_ns, periodic_ew);
            result.periodic_ns = periodic_ns;
            result.periodic_ew = periodic_ew;
            have_periodicity = true;
        }

        const auto&  nodes     = rl.nodes();
        const auto&  distances = rl.distances();
        const double total_m   = rl.total_length();

        const std::int32_t rid = static_cast<std::int32_t>(ri);
        const std::int64_t basin_flat =
            basin_outlet_flats.empty() ? -1 : basin_outlet_flats[ri];

        for (std::size_t i = 0; i + 1 < nodes.size(); ++i)
        {
            const lsdtt_pixel_rc& a = nodes[i];
            const lsdtt_pixel_rc& b = nodes[i + 1];

            ++result.n_edges_total;

            const ridgeline_step_kind kind =
                classify_ridgeline_step(a, b, dgm, periodic_ns, periodic_ew);

            const double dist_m =
                (i < distances.size()) ? distances[i] : 0.0;
            const std::int32_t node_index = static_cast<std::int32_t>(i);

            switch (kind)
            {
                case ridgeline_step_kind::repeated:
                    ++result.n_edges_repeated;
                    continue;

                case ridgeline_step_kind::disjoint:
                    ++result.n_edges_disjoint;
                    continue;

                case ridgeline_step_kind::diagonal:
                {
                    ++result.n_edges_diagonal;

                    if (opts.diagonal_policy == diagonal_step_policy::fail)
                    {
                        throw std::runtime_error(
                            "lsdtt_hilltop_flow_routing_cpu: ridgeline "
                            + std::to_string(ri) + " has a diagonal step at node "
                            + std::to_string(i)
                            + ". A diagonal corner-to-corner step does not "
                              "separate a well-defined pixel pair. Set "
                              "diagonal_policy to skip or bridge.");
                    }
                    if (opts.diagonal_policy == diagonal_step_policy::skip)
                        continue;

                    // bridge: route through the intervening corner
                    // (row of a, col of b), giving two cardinal edges.
                    const int dr = step_delta(a.row, b.row, dgm.n_rows, periodic_ns);
                    const int dc = step_delta(a.col, b.col, dgm.n_cols, periodic_ew);

                    // first leg: a -> (a.row, b.col), horizontal
                    emit_edge(a, 0, dc, rid, node_index, dist_m, total_m, basin_flat);

                    // second leg: (a.row, b.col) -> b, vertical
                    lsdtt_pixel_rc mid{a.row, b.col};
                    emit_edge(mid, dr, 0, rid, node_index, dist_m, total_m, basin_flat);
                    continue;
                }

                case ridgeline_step_kind::cardinal:
                default:
                {
                    ++result.n_edges_cardinal;
                    const int dr = step_delta(a.row, b.row, dgm.n_rows, periodic_ns);
                    const int dc = step_delta(a.col, b.col, dgm.n_cols, periodic_ew);
                    emit_edge(a, dr, dc, rid, node_index, dist_m, total_m, basin_flat);
                    continue;
                }
            }
        }
    }

    return result;
}

} // namespace lsdtt
