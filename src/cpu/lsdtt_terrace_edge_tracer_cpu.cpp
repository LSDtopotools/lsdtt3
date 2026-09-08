// lsdtt_terrace_edge_tracer_cpu.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html

#include "lsdtt_terrace_edge_tracer_cpu.hpp"

#include "lsdtt_raster_skeleton_graph_cpu.hpp"
#include "../utils/lsdtt_raster_geometry.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <numeric>
#include <queue>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <utility>
#include <vector>

namespace lsdtt
{
namespace
{

struct rc_pixel
{
    int row = -1;
    int col = -1;
};

static bool is_valid_scalar_(
    const lsdtt_raster_meta& meta,
    const std::vector<float>& data,
    const std::vector<uint8_t>& mask,
    std::size_t idx)
{
    return mask[idx] != 0 &&
           !meta.is_nodata(data[idx]) &&
           std::isfinite(data[idx]);
}

static double median_copy_(std::vector<double> vals)
{
    if (vals.empty())
        return std::numeric_limits<double>::quiet_NaN();

    const std::size_t n = vals.size();
    const std::size_t mid = n / 2;

    std::nth_element(vals.begin(), vals.begin() + mid, vals.end());
    double med = vals[mid];

    if (n % 2 == 0)
    {
        std::nth_element(vals.begin(), vals.begin() + mid - 1, vals.end());
        med = 0.5 * (med + vals[mid - 1]);
    }

    return med;
}

// -----------------------------------------------------------------------------
// label_centreline_components_d8_
//
// Labels 8-connected foreground components of the centreline.
// labels[idx] = 0 -> background
// labels[idx] = 1..n -> component ID
// comp_pixels[label] -> list of pixel indices in that component
// -----------------------------------------------------------------------------
static std::vector<int> label_centreline_components_d8_(
    const lsdtt_raster_cpu&                  centreline,
    std::vector<std::vector<std::size_t>>&   comp_pixels)
{
    const lsdtt_raster_meta& meta = centreline.get_meta();
    const int nrows = meta.n_rows;
    const int ncols = meta.n_cols;
    const std::size_t npix =
        static_cast<std::size_t>(nrows) * static_cast<std::size_t>(ncols);

    const std::vector<float>&   data = centreline.get_data();
    const std::vector<uint8_t>& mask = centreline.get_mask();

    std::vector<int> labels(npix, 0);
    comp_pixels.clear();
    comp_pixels.push_back({}); // label 0 unused

    auto idx_of = [ncols](int row, int col) -> std::size_t
    {
        return static_cast<std::size_t>(row) * static_cast<std::size_t>(ncols) +
               static_cast<std::size_t>(col);
    };

    constexpr int DR8[8] = {-1,-1,-1, 0, 0, 1, 1, 1};
    constexpr int DC8[8] = {-1, 0, 1,-1, 1,-1, 0, 1};

    int next_label = 1;

    for (int row = 0; row < nrows; ++row)
    {
        for (int col = 0; col < ncols; ++col)
        {
            const std::size_t seed = idx_of(row, col);

            if (!is_valid_scalar_(meta, data, mask, seed)) continue;
            if (data[seed] <= 0.5f) continue;
            if (labels[seed] != 0) continue;

            labels[seed] = next_label;
            comp_pixels.push_back({});
            std::queue<std::size_t> q;
            q.push(seed);

            while (!q.empty())
            {
                const std::size_t cur = q.front();
                q.pop();
                comp_pixels[next_label].push_back(cur);

                const int cr = static_cast<int>(cur / static_cast<std::size_t>(ncols));
                const int cc = static_cast<int>(cur % static_cast<std::size_t>(ncols));

                for (int d = 0; d < 8; ++d)
                {
                    const int nr = cr + DR8[d];
                    const int nc = cc + DC8[d];
                    if (nr < 0 || nr >= nrows || nc < 0 || nc >= ncols) continue;

                    const std::size_t nidx = idx_of(nr, nc);
                    if (!is_valid_scalar_(meta, data, mask, nidx)) continue;
                    if (data[nidx] <= 0.5f) continue;
                    if (labels[nidx] != 0) continue;

                    labels[nidx] = next_label;
                    q.push(nidx);
                }
            }

            ++next_label;
        }
    }

    return labels;
}

// -----------------------------------------------------------------------------
// gather_endpoint_pixels_
//
// Uses endpoint_mask if non-empty; otherwise derives endpoints internally.
// Returns endpoint coordinates grouped by centreline component label.
// -----------------------------------------------------------------------------
static std::vector<std::vector<rc_pixel>> gather_endpoint_pixels_(
    const lsdtt_raster_cpu&      centreline,
    const lsdtt_raster_cpu&      endpoint_mask,
    const std::vector<int>&      labels,
    int                          n_labels)
{
    const lsdtt_raster_meta& cmeta = centreline.get_meta();
    const int nrows = cmeta.n_rows;
    const int ncols = cmeta.n_cols;

    std::vector<std::vector<rc_pixel>> endpoints_by_label(
        static_cast<std::size_t>(n_labels + 1));

    auto idx_of = [ncols](int row, int col) -> std::size_t
    {
        return static_cast<std::size_t>(row) * static_cast<std::size_t>(ncols) +
               static_cast<std::size_t>(col);
    };

    std::vector<rc_pixel> endpoint_pixels;

    if (endpoint_mask.get_data().empty())
    {
        const skeleton_node_result snr = classify_skeleton_nodes(centreline);
        endpoint_pixels.reserve(snr.endpoints.size());
        for (const skeleton_node& nd : snr.endpoints)
        {
            endpoint_pixels.push_back({nd.row, nd.col});
        }
    }
    else
    {
        const lsdtt_raster_meta& emeta = endpoint_mask.get_meta();
        if (!same_grid(cmeta, emeta))
        {
            throw std::runtime_error(
                "lsdtt::trace_terrace_edges_along_contours: endpoint_mask is not co-registered with centreline"
            );
        }

        const std::vector<float>&   ed = endpoint_mask.get_data();
        const std::vector<uint8_t>& em = endpoint_mask.get_mask();

        for (int row = 0; row < nrows; ++row)
        {
            for (int col = 0; col < ncols; ++col)
            {
                const std::size_t idx = idx_of(row, col);
                if (!is_valid_scalar_(emeta, ed, em, idx)) continue;
                if (ed[idx] <= 0.5f) continue;

                endpoint_pixels.push_back({row, col});
            }
        }
    }

    for (const rc_pixel& ep : endpoint_pixels)
    {
        const std::size_t idx = idx_of(ep.row, ep.col);
        const int lbl = labels[idx];
        if (lbl > 0 && lbl <= n_labels)
            endpoints_by_label[lbl].push_back(ep);
    }

    return endpoints_by_label;
}

static terrace_edge_seed_info compute_seed_stats_(
    int                                 label,
    const std::vector<std::size_t>&     pixels,
    const lsdtt_raster_cpu&             dem,
    int                                 n_endpoints)
{
    terrace_edge_seed_info info;
    info.label = static_cast<std::int32_t>(label);
    info.n_pixels = static_cast<std::int32_t>(pixels.size());
    info.n_endpoints = static_cast<std::int32_t>(n_endpoints);

    const lsdtt_raster_meta& meta = dem.get_meta();
    const std::vector<float>&   data = dem.get_data();
    const std::vector<uint8_t>& mask = dem.get_mask();

    std::vector<double> z;
    z.reserve(pixels.size());

    for (const std::size_t idx : pixels)
    {
        if (!is_valid_scalar_(meta, data, mask, idx)) continue;
        z.push_back(static_cast<double>(data[idx]));
    }

    info.n_valid_elevation_pixels = static_cast<std::int32_t>(z.size());

    if (z.empty())
        return info;

    const double sum = std::accumulate(z.begin(), z.end(), 0.0);
    const double mean = sum / static_cast<double>(z.size());

    double var = 0.0;
    for (double v : z)
    {
        const double dv = v - mean;
        var += dv * dv;
    }
    var /= static_cast<double>(z.size());

    info.mean_elevation_metres   = mean;
    info.median_elevation_metres = median_copy_(z);
    info.stddev_elevation_metres = std::sqrt(std::max(0.0, var));

    return info;
}

static bool component_is_traceable_(
    const terrace_edge_seed_info&      info,
    const terrace_edge_trace_request&  request)
{
    if (info.n_pixels < request.min_seed_component_pixels)
        return false;

    if (info.n_valid_elevation_pixels < 2)
        return false;

    if (!(info.stddev_elevation_metres <=
          request.max_seed_elevation_stddev_metres))
        return false;

    if (request.require_exactly_two_endpoints)
        return info.n_endpoints == 2;

    return info.n_endpoints >= 1;
}

// -----------------------------------------------------------------------------
// trace_from_endpoint_
//
// Greedy contour-guided outward tracing from one endpoint of one accepted seed.
// -----------------------------------------------------------------------------
static terrace_edge_trace_record trace_from_endpoint_(
    int                                     seed_label,
    int                                     endpoint_index,
    const rc_pixel&                         endpoint,
    double                                  target_elevation,
    const std::vector<int>&                 labels,
    const lsdtt_raster_cpu&                 dem,
    const lsdtt_raster_cpu&                 centreline,
    const lsdtt_raster_cpu&                 flat_mask,
    const terrace_edge_trace_request&       request,
    std::vector<uint8_t>&                   extension_work)
{
    terrace_edge_trace_record rec;
    rec.seed_label = static_cast<std::int32_t>(seed_label);
    rec.endpoint_index = static_cast<std::int32_t>(endpoint_index);
    rec.start_row = endpoint.row;
    rec.start_col = endpoint.col;
    rec.target_elevation_metres = target_elevation;

    const lsdtt_raster_meta& meta = dem.get_meta();
    const int nrows = meta.n_rows;
    const int ncols = meta.n_cols;

    const std::vector<float>&   zd = dem.get_data();
    const std::vector<uint8_t>& zm = dem.get_mask();

    const std::vector<float>&   cd = centreline.get_data();
    const std::vector<uint8_t>& cm = centreline.get_mask();

    const lsdtt_raster_meta& fmeta = flat_mask.get_meta();
    const std::vector<float>&   fd = flat_mask.get_data();
    const std::vector<uint8_t>& fm = flat_mask.get_mask();

    auto idx_of = [ncols](int row, int col) -> std::size_t
    {
        return static_cast<std::size_t>(row) * static_cast<std::size_t>(ncols) +
               static_cast<std::size_t>(col);
    };

    auto in_bounds = [nrows, ncols](int row, int col) -> bool
    {
        return row >= 0 && row < nrows && col >= 0 && col < ncols;
    };

    auto flat_fg = [&](std::size_t idx) -> bool
    {
        return fm[idx] != 0 &&
               !fmeta.is_nodata(fd[idx]) &&
               std::isfinite(fd[idx]) &&
               fd[idx] > 0.5f;
    };

    auto dem_valid = [&](std::size_t idx) -> bool
    {
        return zm[idx] != 0 &&
               !meta.is_nodata(zd[idx]) &&
               std::isfinite(zd[idx]);
    };

    constexpr int DR8[8] = {-1,-1,-1, 0, 0, 1, 1, 1};
    constexpr int DC8[8] = {-1, 0, 1,-1, 1,-1, 0, 1};

    const std::size_t seed_idx = idx_of(endpoint.row, endpoint.col);

    // Find the unique inward neighbour on the same labelled seed line.
    std::vector<std::size_t> inward_same_label;
    inward_same_label.reserve(2);

    for (int d = 0; d < 8; ++d)
    {
        const int nr = endpoint.row + DR8[d];
        const int nc = endpoint.col + DC8[d];
        if (!in_bounds(nr, nc)) continue;

        const std::size_t nidx = idx_of(nr, nc);
        if (labels[nidx] == seed_label)
            inward_same_label.push_back(nidx);
    }

    if (inward_same_label.size() != 1)
    {
        rec.stop_reason = "bad_seed_endpoint_neighbourhood";
        return rec;
    }

    std::size_t prev = inward_same_label.front();
    std::size_t curr = seed_idx;

    int heading_dr =
        static_cast<int>(curr / static_cast<std::size_t>(ncols)) -
        static_cast<int>(prev / static_cast<std::size_t>(ncols));
    int heading_dc =
        static_cast<int>(curr % static_cast<std::size_t>(ncols)) -
        static_cast<int>(prev % static_cast<std::size_t>(ncols));

    std::unordered_set<std::size_t> local_visited;
    local_visited.reserve(static_cast<std::size_t>(
        std::max(8, request.max_trace_length_pixels * 2)));
    local_visited.insert(seed_idx);

    for (int step = 0; step < request.max_trace_length_pixels; ++step)
    {
        const int crow = static_cast<int>(curr / static_cast<std::size_t>(ncols));
        const int ccol = static_cast<int>(curr % static_cast<std::size_t>(ncols));

        // Stop if we are already adjacent to another labelled line.
        for (int d = 0; d < 8; ++d)
        {
            const int nr = crow + DR8[d];
            const int nc = ccol + DC8[d];
            if (!in_bounds(nr, nc)) continue;

            const std::size_t nidx = idx_of(nr, nc);
            const int other_label = labels[nidx];

            if (other_label > 0 && other_label != seed_label)
            {
                rec.hit_label = static_cast<std::int32_t>(other_label);
                rec.stop_reason = "hit_line";
                return rec;
            }

            if (extension_work[nidx] != 0 && nidx != curr)
            {
                rec.stop_reason = "hit_extension";
                return rec;
            }
        }

        double best_score = std::numeric_limits<double>::infinity();
        std::size_t best_idx = static_cast<std::size_t>(-1);

        for (int d = 0; d < 8; ++d)
        {
            const int nr = crow + DR8[d];
            const int nc = ccol + DC8[d];
            if (!in_bounds(nr, nc)) continue;

            if (!request.allow_diagonal_steps &&
                std::abs(DR8[d]) + std::abs(DC8[d]) == 2)
            {
                continue;
            }

            const std::size_t nidx = idx_of(nr, nc);
            if (nidx == prev) continue;
            if (local_visited.find(nidx) != local_visited.end()) continue;

            // Must be inside flat mask.
            if (!flat_fg(nidx)) continue;

            // Need valid elevation.
            if (!dem_valid(nidx)) continue;

            // Do not step back onto the original centreline of this same label.
            if (labels[nidx] == seed_label) continue;

            // Another line is handled as a stop condition above.
            if (labels[nidx] > 0 && labels[nidx] != seed_label) continue;

            const double zdiff =
                std::abs(static_cast<double>(zd[nidx]) - target_elevation);

            if (zdiff > request.contour_elevation_tolerance_metres)
                continue;

            const int step_dr = nr - crow;
            const int step_dc = nc - ccol;

            const double hnorm =
                std::sqrt(static_cast<double>(heading_dr * heading_dr +
                                              heading_dc * heading_dc));
            const double snorm =
                std::sqrt(static_cast<double>(step_dr * step_dr +
                                              step_dc * step_dc));

            double cosang = 1.0;
            if (hnorm > 0.0 && snorm > 0.0)
            {
                cosang =
                    (static_cast<double>(heading_dr * step_dr +
                                         heading_dc * step_dc)) / (hnorm * snorm);
                cosang = std::max(-1.0, std::min(1.0, cosang));
            }

            const double heading_penalty = 1.0 - cosang;

            const double score =
                request.elevation_weight * zdiff +
                request.heading_weight   * heading_penalty;

            if (score < best_score)
            {
                best_score = score;
                best_idx = nidx;
            }
        }

        if (best_idx == static_cast<std::size_t>(-1))
        {
            rec.stop_reason = "no_next_step";
            return rec;
        }

        // Burn the pixel into the extension mask unless already part of original centreline.
        const bool is_original_line =
            cm[best_idx] != 0 &&
            !centreline.get_meta().is_nodata(cd[best_idx]) &&
            std::isfinite(cd[best_idx]) &&
            cd[best_idx] > 0.5f;

        if (!is_original_line)
        {
            extension_work[best_idx] = uint8_t{1};
            ++rec.traced_pixels;
        }

        local_visited.insert(best_idx);

        const std::size_t new_prev = curr;
        const std::size_t new_curr = best_idx;

        heading_dr =
            static_cast<int>(new_curr / static_cast<std::size_t>(ncols)) -
            static_cast<int>(new_prev / static_cast<std::size_t>(ncols));
        heading_dc =
            static_cast<int>(new_curr % static_cast<std::size_t>(ncols)) -
            static_cast<int>(new_prev % static_cast<std::size_t>(ncols));

        prev = new_prev;
        curr = new_curr;
    }

    rec.stop_reason = "max_trace_length";
    return rec;
}

} // unnamed namespace

// -----------------------------------------------------------------------------
// trace_terrace_edges_along_contours
// -----------------------------------------------------------------------------
terrace_edge_trace_result trace_terrace_edges_along_contours(
    const lsdtt_raster_cpu&           dem,
    const lsdtt_raster_cpu&           centreline,
    const lsdtt_raster_cpu&           flat_mask,
    const lsdtt_raster_cpu&           endpoint_mask,
    const terrace_edge_trace_request& request)
{
    const lsdtt_raster_meta& dmeta = dem.get_meta();
    const lsdtt_raster_meta& cmeta = centreline.get_meta();
    const lsdtt_raster_meta& fmeta = flat_mask.get_meta();

    if (!dmeta.is_valid() || !cmeta.is_valid() || !fmeta.is_valid())
    {
        throw std::runtime_error(
            "lsdtt::trace_terrace_edges_along_contours: one or more input raster metadata objects are invalid"
        );
    }

    if (!same_grid(dmeta, cmeta) || !same_grid(dmeta, fmeta))
    {
        throw std::runtime_error(
            "lsdtt::trace_terrace_edges_along_contours: dem, centreline, and flat_mask are not co-registered"
        );
    }

    if (!endpoint_mask.get_data().empty() && !same_grid(dmeta, endpoint_mask.get_meta()))
    {
        throw std::runtime_error(
            "lsdtt::trace_terrace_edges_along_contours: endpoint_mask is not co-registered with dem"
        );
    }

    const int nrows = dmeta.n_rows;
    const int ncols = dmeta.n_cols;
    const std::size_t npix =
        static_cast<std::size_t>(nrows) * static_cast<std::size_t>(ncols);

    const std::vector<float>&   dd = dem.get_data();
    const std::vector<uint8_t>& dm = dem.get_mask();

    // -------------------------------------------------------------------------
    // Label centreline components.
    // -------------------------------------------------------------------------
    std::vector<std::vector<std::size_t>> comp_pixels;
    const std::vector<int> labels =
        label_centreline_components_d8_(centreline, comp_pixels);

    const int n_labels =
        static_cast<int>(comp_pixels.size()) - 1;

    // -------------------------------------------------------------------------
    // Endpoint pixels grouped by component label.
    // -------------------------------------------------------------------------
    const std::vector<std::vector<rc_pixel>> endpoints_by_label =
        gather_endpoint_pixels_(centreline, endpoint_mask, labels, n_labels);

    // -------------------------------------------------------------------------
    // Seed statistics and acceptance.
    // -------------------------------------------------------------------------
    std::vector<terrace_edge_seed_info> seed_info(
        static_cast<std::size_t>(n_labels + 1));

    std::vector<uint8_t> accepted_label(
        static_cast<std::size_t>(n_labels + 1), uint8_t{0});

    for (int lbl = 1; lbl <= n_labels; ++lbl)
    {
        terrace_edge_seed_info info =
            compute_seed_stats_(
                lbl,
                comp_pixels[static_cast<std::size_t>(lbl)],
                dem,
                static_cast<int>(endpoints_by_label[static_cast<std::size_t>(lbl)].size()));

        info.accepted_seed =
            info.n_pixels >= request.min_seed_component_pixels &&
            info.n_valid_elevation_pixels >= 2 &&
            info.stddev_elevation_metres <= request.max_seed_elevation_stddev_metres;

        info.traceable = component_is_traceable_(info, request);

        seed_info[static_cast<std::size_t>(lbl)] = info;

        if (info.traceable)
            accepted_label[static_cast<std::size_t>(lbl)] = uint8_t{1};
    }

    // -------------------------------------------------------------------------
    // Extension work mask.
    // -------------------------------------------------------------------------
    std::vector<uint8_t> extension_work(npix, uint8_t{0});
    std::vector<terrace_edge_trace_record> trace_records;

    for (int lbl = 1; lbl <= n_labels; ++lbl)
    {
        if (!accepted_label[static_cast<std::size_t>(lbl)])
            continue;

        const terrace_edge_seed_info& info =
            seed_info[static_cast<std::size_t>(lbl)];

        const auto& eps = endpoints_by_label[static_cast<std::size_t>(lbl)];

        for (std::size_t ei = 0; ei < eps.size(); ++ei)
        {
            terrace_edge_trace_record rec =
                trace_from_endpoint_(
                    lbl,
                    static_cast<int>(ei),
                    eps[ei],
                    info.median_elevation_metres,
                    labels,
                    dem,
                    centreline,
                    flat_mask,
                    request,
                    extension_work);

            trace_records.push_back(std::move(rec));
        }
    }

    // -------------------------------------------------------------------------
    // Build output rasters.
    // Use DEM-valid domain as the output valid mask.
    // -------------------------------------------------------------------------
    std::vector<uint8_t> out_valid(npix, uint8_t{0});
    for (std::size_t i = 0; i < npix; ++i)
    {
        if (dm[i] != 0 &&
            !dmeta.is_nodata(dd[i]) &&
            std::isfinite(dd[i]))
        {
            out_valid[i] = uint8_t{1};
        }
    }

    // labelled_line_components
    std::vector<float> labelled_data(npix, 0.0f);
    std::vector<uint8_t> labelled_mask = out_valid;
    for (std::size_t i = 0; i < npix; ++i)
    {
        if (!out_valid[i]) continue;
        labelled_data[i] = static_cast<float>(labels[i]);
    }

    // seed_mask
    std::vector<float> seed_data(npix, 0.0f);
    std::vector<uint8_t> seed_mask = out_valid;
    for (int lbl = 1; lbl <= n_labels; ++lbl)
    {
        if (!seed_info[static_cast<std::size_t>(lbl)].accepted_seed)
            continue;

        for (const std::size_t idx : comp_pixels[static_cast<std::size_t>(lbl)])
        {
            if (!out_valid[idx]) continue;
            seed_data[idx] = 1.0f;
        }
    }

    // trace_extension_mask
    std::vector<float> ext_data(npix, 0.0f);
    std::vector<uint8_t> ext_mask = out_valid;
    for (std::size_t i = 0; i < npix; ++i)
    {
        if (!out_valid[i]) continue;
        ext_data[i] = (extension_work[i] != 0) ? 1.0f : 0.0f;
    }

    // connected_mask = original centreline OR extension
    const std::vector<float>&   cd = centreline.get_data();
    const std::vector<uint8_t>& cm = centreline.get_mask();

    std::vector<float> connected_data(npix, 0.0f);
    std::vector<uint8_t> connected_mask = out_valid;

    for (std::size_t i = 0; i < npix; ++i)
    {
        if (!out_valid[i]) continue;

        const bool orig_fg =
            cm[i] != 0 &&
            !centreline.get_meta().is_nodata(cd[i]) &&
            std::isfinite(cd[i]) &&
            cd[i] > 0.5f;

        connected_data[i] =
            (orig_fg || extension_work[i] != 0) ? 1.0f : 0.0f;
    }

    terrace_edge_trace_result result;
    result.labelled_line_components =
        lsdtt_raster_cpu(std::move(labelled_data), std::move(labelled_mask), dmeta);
    result.seed_mask =
        lsdtt_raster_cpu(std::move(seed_data), std::move(seed_mask), dmeta);
    result.trace_extension_mask =
        lsdtt_raster_cpu(std::move(ext_data), std::move(ext_mask), dmeta);
    result.connected_mask =
        lsdtt_raster_cpu(std::move(connected_data), std::move(connected_mask), dmeta);

    result.seed_info = std::move(seed_info);
    result.trace_records = std::move(trace_records);

    return result;
}

} // namespace lsdtt
