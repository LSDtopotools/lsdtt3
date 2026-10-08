// lsdtt_hilltop_flow_routing_cpu.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
//
// See lsdtt_hilltop_flow_routing_cpu.hpp for the divide-grid convention and the
// side-assignment rule.

#include "lsdtt_hilltop_flow_routing_cpu.hpp"
#include "../utils/lsdtt_output_file.hpp"

#include <algorithm>
#include <sstream>
#include <stdexcept>
#include <string>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <unordered_map>
#include <map>
#include <unordered_set>
#include <utility>

#include "lsdtt_streamlines_cpu.hpp"

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
    const lsdtt_raster_cpu*                 channel_mask,
    const lsdtt_raster_cpu*                 gradient,
    const lsdtt_raster_cpu*                 curvature,
    lsdtt_hilltop_filter_report&            report)
{
    const int nr = dem_meta.n_rows;
    const int nc = dem_meta.n_cols;

    auto check_size = [&](const lsdtt_raster_cpu* r, const char* what) {
        if (r && (r->get_meta().n_rows != nr || r->get_meta().n_cols != nc))
            throw std::invalid_argument(
                std::string("lsdtt_hilltop_flow_routing_cpu: ") + what +
                " raster does not match the DEM grid.");
    };

    // Read a value only where the raster says it is valid.
    auto valid_at = [](const lsdtt_raster_cpu* r, int row, int col) -> bool {
        return r && r->get_mask_array().at(row, col) != 0;
    };
    check_size(channel_mask, "channel mask");
    check_size(gradient, "gradient");
    check_size(curvature, "curvature");

    report = lsdtt_hilltop_filter_report{};
    report.n_in = static_cast<std::int64_t>(pixels.size());

    if (channel_mask)
        for (const auto& p : pixels)
            if (valid_at(channel_mask, p.row, p.col) &&
                channel_mask->get_data_array().at(p.row, p.col) > 0.5f)
                ++report.n_on_channel_before;

    const int pad = opts.stream_network_padding_pixels;

    std::vector<lsdtt_hilltop_pixel> out;
    out.reserve(pixels.size());

    for (const auto& p : pixels)
    {
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
                    if (valid_at(channel_mask, rr, cc) &&
                        channel_mask->get_data_array().at(rr, cc) > 0.5f)
                        near = true;
                }
            if (near) { ++report.n_dropped_near_channel; continue; }
        }

        if (gradient && opts.threshold_hilltop_gradient > 0.0f &&
            valid_at(gradient, p.row, p.col) &&
            gradient->get_data_array().at(p.row, p.col) >
                opts.threshold_hilltop_gradient)
        { ++report.n_dropped_steep; continue; }

        if (curvature && opts.remove_positive_curvature_hilltops &&
            valid_at(curvature, p.row, p.col) &&
            curvature->get_data_array().at(p.row, p.col) > 0.0f)
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


// =============================================================================
// PHASE 2 -- HILLSLOPE METRICS
// =============================================================================

namespace
{

constexpr double k_nan = std::numeric_limits<double>::quiet_NaN();

hilltop_trace_status status_from_outcome(streamline_outcome o)
{
    switch (o)
    {
        case streamline_outcome::reached_sink:  return hilltop_trace_status::reached_channel;
        case streamline_outcome::hit_edge:      return hilltop_trace_status::hit_edge;
        case streamline_outcome::loop_detected: return hilltop_trace_status::loop_detected;
        case streamline_outcome::stalled_on_flat: return hilltop_trace_status::stalled_on_flat;
        case streamline_outcome::no_aspect:     return hilltop_trace_status::no_aspect;
    }
    return hilltop_trace_status::hit_edge;
}

} // namespace

lsdtt_hillslope_metrics_result compute_hillslope_metrics(
    const std::vector<lsdtt_hilltop_pixel>& hilltop_pixels,
    const lsdtt_raster_cpu&                 dem,
    const lsdtt_raster_cpu&                 channel_mask,
    const lsdtt_raster_cpu*                 curvature,
    const lsdtt_raster_cpu*                 gradient,
    const lsdtt_raster_cpu*                 plan_curvature,
    const lsdtt_hillslope_metrics_options&  opts)
{
    if (!(opts.critical_slope_Sc > 0.0))
        throw std::invalid_argument(
            "lsdtt_hilltop_flow_routing_cpu: critical_slope_Sc must be > 0 (got "
            + std::to_string(opts.critical_slope_Sc) + "). E* and R* are "
              "undefined otherwise.");

    const lsdtt_raster_meta& meta = dem.get_meta();
    const int nr = meta.n_rows;
    const int nc = meta.n_cols;
    const std::size_t ncell = static_cast<std::size_t>(nr) * nc;

    auto check = [&](const lsdtt_raster_cpu* r, const char* what) {
        if (r && (r->get_meta().n_rows != nr || r->get_meta().n_cols != nc))
            throw std::invalid_argument(
                std::string("lsdtt_hilltop_flow_routing_cpu: ") + what +
                " raster does not match the DEM grid.");
    };

    // Read a value only where the raster's own mask says it is valid; a nodata
    // sentinel must never reach a metric.
    auto sample = [](const lsdtt_raster_cpu* r, int row, int col) -> double {
        if (!r) return k_nan;
        if (r->get_mask_array().at(row, col) == 0) return k_nan;
        const double v = r->get_data_array().at(row, col);
        return std::isnan(v) ? k_nan : v;
    };
    check(curvature, "curvature");
    check(gradient, "gradient");
    check(plan_curvature, "plan curvature");

    // ---- build the source mask from the hilltop pixels ----------------------
    std::vector<float>        src(ncell, 0.0f);
    std::vector<std::uint8_t> src_valid(ncell, 1);
    for (const auto& p : hilltop_pixels)
        if (p.flat >= 0 && static_cast<std::size_t>(p.flat) < ncell)
            src[static_cast<std::size_t>(p.flat)] = 1.0f;

    lsdtt_raster_cpu source_mask(std::move(src), std::move(src_valid), meta);

    // ---- trace ---------------------------------------------------------------
    // save_nodes is required, not optional: the divergent/planar counts are
    // sampled along the path and there is no other way to know which cells the
    // trace visited.
    lsdtt_streamline_options topts;
    topts.save_nodes = true;
    topts.max_steps  = opts.max_steps;
    topts.flat_terminus_gradient = opts.flat_terminus_gradient;
    topts.n_threads  = opts.n_threads;

    const lsdtt_streamline_result traces =
        trace_streamlines(dem, source_mask, channel_mask, topts);

    // ---- index the traces by source pixel ------------------------------------
    std::unordered_map<std::int64_t, const lsdtt_streamline*> by_source;
    by_source.reserve(traces.streamlines.size() * 2);
    for (const auto& t : traces.streamlines)
    {
        if (t.source_row < 0 || t.source_col < 0) continue;
        by_source.emplace(
            static_cast<std::int64_t>(t.source_row) * nc + t.source_col, &t);
    }

    const double cell = meta.cell_size;
    const double Sc   = opts.critical_slope_Sc;

    lsdtt_hillslope_metrics_result result;
    result.rows.reserve(hilltop_pixels.size());

    for (const auto& p : hilltop_pixels)
    {
        lsdtt_hilltop_metrics_row row;
        row.ridge = p;

        row.hilltop_curvature = sample(curvature, p.row, p.col);
        row.hilltop_gradient  = sample(gradient,  p.row, p.col);

        auto it = by_source.find(p.flat);
        if (it == by_source.end())
        {
            // Never offered to the tracer: it sits on the outer ring.
            row.status = hilltop_trace_status::not_traced_border;
            ++result.n_not_traced_border;
            result.rows.push_back(row);
            continue;
        }

        const lsdtt_streamline& t = *it->second;
        row.status = status_from_outcome(t.outcome);

        switch (row.status)
        {
            case hilltop_trace_status::reached_channel: ++result.n_reached_channel; break;
            case hilltop_trace_status::hit_edge:        ++result.n_hit_edge;        break;
            case hilltop_trace_status::loop_detected:   ++result.n_loop_detected;   break;
            case hilltop_trace_status::stalled_on_flat: ++result.n_stalled_on_flat; break;
            case hilltop_trace_status::no_aspect:       ++result.n_no_aspect;       break;
            default: break;
        }

        // Path diagnostics are meaningful for any trace that actually moved.
        if (plan_curvature && !t.nodes.empty())
        {
            const double thresh = opts.plan_curvature_threshold;
            for (const auto& nd : t.nodes)
            {
                // map the sub-pixel node back to its cell
                const int c = static_cast<int>(
                    std::floor((nd.x - meta.x_min) / cell));
                const int r = nr - 1 - static_cast<int>(
                    std::floor((nd.y - meta.y_min) / cell));
                if (r < 0 || r >= nr || c < 0 || c >= nc) continue;

                const double pc = sample(plan_curvature, r, c);
                if (std::isnan(pc)) continue;
                if (std::abs(pc) > thresh) ++row.divergent_count;
                else                       ++row.planar_count;
            }
        }

        if (row.status != hilltop_trace_status::reached_channel)
        {
            // Lh and R are meaningless for a trace that never found a channel;
            // leave them NaN rather than reporting the distance it happened to
            // cover before failing.
            result.rows.push_back(row);
            continue;
        }

        row.hillslope_length_m = t.length_m;
        row.relief_m           = t.relief_m;
        row.channel_row        = t.sink_row;
        row.channel_col        = t.sink_col;

        const double dr = static_cast<double>(t.sink_row - t.source_row) * cell;
        const double dc = static_cast<double>(t.sink_col - t.source_col) * cell;
        row.euclidean_distance_m = std::sqrt(dr * dr + dc * dc);

        if (t.length_m > 0.0)
        {
            // Legacy definition: S is relief over path length, not a mean of
            // sampled gradients (LSDFlowInfo.cpp:7028). R* is therefore S / Sc.
            row.mean_path_gradient = t.relief_m / t.length_m;
            row.R_star             = t.relief_m / (t.length_m * Sc);

            if (!std::isnan(row.hilltop_curvature))
                row.E_star = 2.0 * std::abs(row.hilltop_curvature)
                             * t.length_m / Sc;
        }
        else
        {
            ++result.n_undefined_nondimensional;
        }

        result.rows.push_back(row);
    }

    return result;
}

// =============================================================================
// write_hilltop_metrics_csv
// =============================================================================
void write_hilltop_metrics_csv(
    const lsdtt_hillslope_metrics_result& result,
    const lsdtt_raster_meta&              meta,
    const std::string&                    out_path)
{
    std::ofstream f = open_output_file(
        out_path, "lsdtt_hilltop_flow_routing_cpu: cannot open ");
    f << std::setprecision(10);
    f << "row,col,easting,northing,"
         "ridgeline_id,ridgeline_node_index,distance_along_ridgeline_m,"
         "ridgeline_total_length_m,ridge_side,basin_outlet_flat,"
         "Cht,hilltop_gradient,mean_path_gradient,relief_m,hillslope_length_m,"
         "euclidean_distance_m,channel_row,channel_col,E_star,R_star,"
         "divergent_count,planar_count,trace_status,"
         "channel_k_sn,channel_chi,channel_elevation,channel_flow_distance_m,"
         "channel_drainage_area_m2,channel_basin_index,channel_segment_id\n";

    auto num = [](double v) -> std::string {
        // NaN is written as an empty field so a reader sees a null rather than
        // the string "nan" or, worse, a magic number.
        if (std::isnan(v)) return std::string();
        std::ostringstream os;
        os << std::setprecision(10) << v;
        return os.str();
    };

    for (const auto& r : result.rows)
    {
        const double easting  = meta.x_min + (r.ridge.col + 0.5) * meta.cell_size;
        const double northing = meta.y_min +
            (meta.n_rows - r.ridge.row - 1 + 0.5) * meta.cell_size;

        f << r.ridge.row << ',' << r.ridge.col << ','
          << easting << ',' << northing << ','
          << r.ridge.ridgeline_id << ',' << r.ridge.ridgeline_node_index << ','
          << r.ridge.distance_along_ridgeline_m << ','
          << r.ridge.ridgeline_total_length_m << ','
          << static_cast<int>(r.ridge.ridge_side) << ','
          << r.ridge.basin_outlet_flat << ','
          << num(r.hilltop_curvature) << ',' << num(r.hilltop_gradient) << ','
          << num(r.mean_path_gradient) << ',' << num(r.relief_m) << ','
          << num(r.hillslope_length_m) << ',' << num(r.euclidean_distance_m) << ','
          << r.channel_row << ',' << r.channel_col << ','
          << num(r.E_star) << ',' << num(r.R_star) << ','
          << r.divergent_count << ',' << r.planar_count << ','
          << static_cast<int>(r.status) << ','
          << num(r.channel_k_sn) << ',' << num(r.channel_chi) << ','
          << num(r.channel_elevation) << ',' << num(r.channel_flow_distance_m) << ','
          << num(r.channel_drainage_area_m2) << ','
          << r.channel_basin_index << ',' << r.channel_segment_id << '\n';
    }
}


// =============================================================================
// PHASE 5 -- THE CHI JOIN
// =============================================================================

lsdtt_channel_property_map build_channel_property_map(
    const lsdtt_chi_network_fit_node_table&    nodes,
    const lsdtt_chi_network_fit_segment_table& segments,
    const std::vector<float>*                  flow_distance,
    const std::vector<float>*                  drainage_area)
{
    lsdtt_channel_property_map out;
    const std::size_t n_nodes = nodes.flat.size();
    if (n_nodes == 0) return out;
    out.reserve(n_nodes * 2);

    // Index the segments by (basin_index, branch_id) so the containment test is
    // done within a branch rather than across the whole network -- chi
    // intervals from different branches overlap and matching across them would
    // silently attach the wrong k_sn.
    struct seg_ref { double lo, hi, k_sn; std::int32_t segment_id; };
    std::map<std::pair<std::int32_t, std::int32_t>, std::vector<seg_ref>> by_group;

    const std::size_t n_seg = segments.segment_id.size();
    for (std::size_t s = 0; s < n_seg; ++s)
    {
        double lo = static_cast<double>(segments.downstream_chi[s]);
        double hi = static_cast<double>(segments.upstream_chi[s]);
        if (lo > hi) std::swap(lo, hi);   // endpoints are stored in either order
        by_group[{segments.basin_index[s], segments.branch_id[s]}].push_back(
            seg_ref{lo, hi, static_cast<double>(segments.k_sn[s]),
                    segments.segment_id[s]});
    }

    for (std::size_t i = 0; i < n_nodes; ++i)
    {
        lsdtt_channel_properties props;
        props.basin_index = nodes.basin_index[i];
        props.branch_id   = nodes.branch_id[i];
        props.chi         = static_cast<double>(nodes.chi[i]);
        props.channel_elevation =
            static_cast<double>(nodes.observed_elevation[i]);

        const std::int64_t flat = nodes.flat[i];
        if (flat >= 0)
        {
            if (flow_distance && static_cast<std::size_t>(flat) < flow_distance->size())
                props.flow_distance_m = (*flow_distance)[static_cast<std::size_t>(flat)];
            if (drainage_area && static_cast<std::size_t>(flat) < drainage_area->size())
                props.drainage_area_m2 = (*drainage_area)[static_cast<std::size_t>(flat)];
        }

        const auto git = by_group.find({nodes.basin_index[i], nodes.branch_id[i]});
        if (git != by_group.end() && !git->second.empty())
        {
            const double chi = props.chi;
            const seg_ref* chosen = nullptr;
            for (const auto& sr : git->second)
                if (chi >= sr.lo && chi <= sr.hi) { chosen = &sr; break; }

            if (!chosen)
            {
                // Just outside every interval: endpoint rounding. Take the
                // nearest segment in this branch rather than dropping the node,
                // which is what lsdviztools3's assign_node_ksn does.
                double best = std::numeric_limits<double>::max();
                for (const auto& sr : git->second)
                {
                    const double d = chi < sr.lo ? (sr.lo - chi)
                                   : chi > sr.hi ? (chi - sr.hi) : 0.0;
                    if (d < best) { best = d; chosen = &sr; }
                }
            }

            if (chosen)
            {
                props.k_sn       = chosen->k_sn;
                props.segment_id = chosen->segment_id;
            }
        }

        if (flat >= 0) out[flat] = props;
    }

    return out;
}

std::int64_t attach_channel_properties(
    lsdtt_hillslope_metrics_result&   metrics,
    const lsdtt_channel_property_map& channel_properties,
    int                               n_cols)
{
    std::int64_t n_joined = 0;

    for (auto& row : metrics.rows)
    {
        // Only a trace that reached a channel has a channel to join to.
        if (row.status != hilltop_trace_status::reached_channel) continue;
        if (row.channel_row < 0 || row.channel_col < 0) continue;

        const std::int64_t flat =
            static_cast<std::int64_t>(row.channel_row) * n_cols + row.channel_col;

        const auto it = channel_properties.find(flat);
        if (it == channel_properties.end()) continue;   // not a fitted node

        const lsdtt_channel_properties& p = it->second;
        row.channel_k_sn             = p.k_sn;
        row.channel_chi              = p.chi;
        row.channel_elevation        = p.channel_elevation;
        row.channel_flow_distance_m  = p.flow_distance_m;
        row.channel_drainage_area_m2 = p.drainage_area_m2;
        row.channel_basin_index      = p.basin_index;
        row.channel_segment_id       = p.segment_id;

        if (!std::isnan(p.k_sn)) ++n_joined;
    }

    return n_joined;
}

} // namespace lsdtt
