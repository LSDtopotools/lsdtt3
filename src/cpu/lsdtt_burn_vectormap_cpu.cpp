// lsdtt_burn_vectormap_cpu.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html

#include "lsdtt_burn_vectormap_cpu.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <variant>

#include "../core/lsdtt_data_loader.hpp"
#include "lsdtt_least_cost_path.hpp"

namespace lsdtt
{

namespace
{

inline bool xy_to_rc(const lsdtt_raster_meta& m, double x, double y,
                     int& row, int& col)
{
    col = static_cast<int>(std::floor((x - m.x_min) / m.cell_size));
    row = static_cast<int>(std::floor((m.y_max() - y) / m.cell_size));
    return row >= 0 && row < m.n_rows && col >= 0 && col < m.n_cols;
}

inline lsdtt_xy rc_to_centre(const lsdtt_raster_meta& m, int row, int col)
{
    return { m.x_min + (static_cast<double>(col) + 0.5) * m.cell_size,
             m.y_max() - (static_cast<double>(row) + 0.5) * m.cell_size };
}

// Densely sample a polyline at `step` metres.
std::vector<lsdtt_xy> sample_polyline(const std::vector<lsdtt_xy>& xy, double step)
{
    std::vector<lsdtt_xy> out;
    if (xy.empty()) return out;
    out.push_back(xy.front());
    for (std::size_t k = 1; k < xy.size(); ++k)
    {
        const lsdtt_xy& a = xy[k - 1];
        const lsdtt_xy& b = xy[k];
        const double seg = std::hypot(b.x - a.x, b.y - a.y);
        const int n = std::max(1, static_cast<int>(seg / step));
        for (int s = 1; s <= n; ++s)
        {
            const double t = static_cast<double>(s) / n;
            out.push_back({ a.x + t * (b.x - a.x), a.y + t * (b.y - a.y) });
        }
    }
    return out;
}

} // namespace

// ---------------------------------------------------------------------------
// channel_heads_from_accumulation
// ---------------------------------------------------------------------------
std::vector<std::int64_t> channel_heads_from_accumulation(
    const lsdtt_flow_base& flow,
    std::int64_t           threshold)
{
    const int n_rows = flow.n_rows();
    const int n_cols = flow.n_cols();
    const std::size_t n = static_cast<std::size_t>(n_rows) * n_cols;

    const lsdtt_raster_result ncp = flow.get_n_contributing_pixels_raster();

    auto is_channel = [&](std::size_t f) -> bool
    {
        if (f >= ncp.data.size()) return false;
        if (f < ncp.mask.size() && ncp.mask[f] == 0) return false;
        return ncp.data[f] >= static_cast<float>(threshold);
    };

    // A cell is "fed" by an upstream channel cell when some valid cell that is a
    // channel drains into it.
    std::vector<std::uint8_t> has_channel_donor(n, 0);
    for (std::size_t f = 0; f < n; ++f)
    {
        const std::int64_t flat = static_cast<std::int64_t>(f);
        if (!flow.is_valid_cell(flat)) continue;
        if (!is_channel(f)) continue;
        if (flow.is_outlet_unchecked(flat)) continue;
        const std::int64_t r = flow.receiver_of_unchecked(flat);
        if (r >= 0 && r != flat)
            has_channel_donor[static_cast<std::size_t>(r)] = 1;
    }

    std::vector<std::int64_t> heads;
    for (std::size_t f = 0; f < n; ++f)
    {
        const std::int64_t flat = static_cast<std::int64_t>(f);
        if (!flow.is_valid_cell(flat)) continue;
        if (is_channel(f) && !has_channel_donor[f])
            heads.push_back(flat);
    }
    return heads;
}

// ---------------------------------------------------------------------------
// burn_vectormap_road_gaps
// ---------------------------------------------------------------------------
vectormap_burn_result burn_vectormap_road_gaps(
    lsdtt_raster_cpu&             dem,
    const lsdtt_flow_base&        flow,
    const std::string&            vector_path,
    const vectormap_burn_options& options)
{
    vectormap_burn_result result;

    const lsdtt_raster_meta& meta = dem.get_meta();
    const int    n_rows = meta.n_rows;
    const int    n_cols = meta.n_cols;
    const double cs     = meta.cell_size;

    const double step = (options.sample_step_m > 0.0)
                            ? options.sample_step_m : 0.5 * cs;
    const int snap_cells =
        std::max(1, static_cast<int>(std::ceil(options.snap_radius_m / cs)));

    // ---- D8 channel mask from the accumulation threshold ------------------
    std::vector<std::uint8_t> channel_mask(
        static_cast<std::size_t>(n_rows) * n_cols, 0);
    {
        const lsdtt_raster_result ncp = flow.get_n_contributing_pixels_raster();
        for (std::size_t i = 0; i < channel_mask.size() && i < ncp.data.size(); ++i)
        {
            const bool valid = (i < ncp.mask.size()) ? (ncp.mask[i] != 0) : true;
            if (valid && ncp.data[i] >=
                static_cast<float>(options.accumulation_threshold_pixels))
            {
                channel_mask[i] = 1;
                ++result.channel_mask_cells;
            }
        }
    }

    // Nearest channel cell to (row,col) within snap_cells. Returns flat index or
    // -1. Bounded Chebyshev-ring search, stops once no closer cell is possible.
    auto nearest_channel = [&](int row, int col) -> std::int64_t
    {
        std::int64_t best_flat = -1;
        double best = std::numeric_limits<double>::infinity();
        for (int r = 0; r <= snap_cells; ++r)
        {
            if (static_cast<double>(r) * cs > best) break;
            auto test = [&](int rr, int cc)
            {
                if (rr < 0 || rr >= n_rows || cc < 0 || cc >= n_cols) return;
                if (!channel_mask[static_cast<std::size_t>(rr) * n_cols + cc]) return;
                const double d = std::hypot(double(rr - row), double(cc - col));
                if (d < best)
                {
                    best = d;
                    best_flat = static_cast<std::int64_t>(rr) * n_cols + cc;
                }
            };
            if (r == 0) { test(row, col); continue; }
            for (int cc = col - r; cc <= col + r; ++cc)
            { test(row - r, cc); test(row + r, cc); }
            for (int rr = row - r + 1; rr <= row + r - 1; ++rr)
            { test(rr, col - r); test(rr, col + r); }
        }
        return best_flat;
    };

    // ---- Read vector lines ------------------------------------------------
    lsdtt_data_loader loader;
    loader.open_vector(vector_path);
    const auto layer_names = loader.get_layer_names();
    if (layer_names.empty())
    {
        loader.close_vector();
        throw std::runtime_error(
            "burn_vectormap_road_gaps: no layers in " + vector_path);
    }
    const lsdtt_vector_layer vlayer = loader.read_vector_layer(layer_names[0]);
    loader.close_vector();

    const lsdtt_array2d_f& zdata = dem.get_data_array();
    auto elev_at = [&](std::int64_t flat) -> double
    {
        const int r = static_cast<int>(flat / n_cols);
        const int c = static_cast<int>(flat % n_cols);
        return static_cast<double>(zdata(r, c));
    };

    // Margin (samples) added on each side of a gap so the trench overlaps the
    // real channel and reliably connects to it.
    const int margin = std::max(2, static_cast<int>(std::ceil(cs / step)));

    for (const lsdtt_vector_record& rec : vlayer.records)
    {
        const lsdtt_line_geom* line = std::get_if<lsdtt_line_geom>(&rec.geometry);
        if (!line || line->points.size() < 2) continue;
        ++result.n_lines;

        const std::vector<lsdtt_xy> samples = sample_polyline(line->points, step);
        const int N = static_cast<int>(samples.size());
        if (N < 2) continue;

        // Classify each sample: on-channel (channel within snap radius) or gap.
        std::vector<std::uint8_t> on_channel(N, 0);
        for (int i = 0; i < N; ++i)
        {
            int row, col;
            if (!xy_to_rc(meta, samples[i].x, samples[i].y, row, col)) continue;
            on_channel[i] = (nearest_channel(row, col) >= 0) ? 1 : 0;
        }

        // Find maximal gap runs (on_channel == 0) that are bracketed by
        // on-channel samples on BOTH sides (a true obstruction between two
        // channel reaches).
        int i = 0;
        while (i < N)
        {
            if (on_channel[i]) { ++i; continue; }
            int j = i;
            while (j < N && !on_channel[j]) ++j;
            // gap run [i, j-1]
            const bool anchored = (i > 0) && (j < N);
            const double gap_len = static_cast<double>(j - i) * step;

            vectormap_burn_gap gap;
            gap.length_m = gap_len;
            ++result.n_gaps_detected;

            if (!anchored)
            {
                gap.burned = false;
                gap.skip_reason = "unanchored (line end has no channel)";
                ++result.n_gaps_skipped;
                result.gaps.push_back(std::move(gap));
                i = j; continue;
            }
            if (gap_len < options.min_gap_length_m)
            {
                i = j; continue;   // noise, not a real gap
            }
            if (options.max_gap_length_m > 0.0 &&
                gap_len > options.max_gap_length_m)
            {
                gap.burned = false;
                gap.skip_reason = "gap longer than max_gap_length_m";
                ++result.n_gaps_skipped;
                result.gaps.push_back(std::move(gap));
                i = j; continue;
            }

            // Build the burn path from samples [i-margin .. j-1+margin], so the
            // trench starts and ends a few cells inside the real channel.
            const int a = std::max(0, i - margin);
            const int b = std::min(N - 1, (j - 1) + margin);

            std::vector<std::int64_t> path;
            for (int k = a; k <= b; ++k)
            {
                int row, col;
                if (!xy_to_rc(meta, samples[k].x, samples[k].y, row, col)) continue;
                const std::int64_t flat =
                    static_cast<std::int64_t>(row) * n_cols + col;
                if (path.empty() || path.back() != flat) path.push_back(flat);
            }
            if (path.size() < 2)
            {
                i = j; continue;
            }

            // Orient so the trench descends (burn_capturing_path descends from
            // front to back).
            if (elev_at(path.front()) < elev_at(path.back()))
                std::reverse(path.begin(), path.end());

            burn_capturing_path(dem, path, options.burn_min_slope);

            gap.burned = true;
            gap.path_xy.reserve(path.size());
            for (std::int64_t f : path)
                gap.path_xy.push_back(rc_to_centre(meta,
                    static_cast<int>(f / n_cols), static_cast<int>(f % n_cols)));
            ++result.n_gaps_burned;
            result.total_burn_length_m += gap_len;
            result.gaps.push_back(std::move(gap));

            i = j;
        }
    }

    return result;
}

} // namespace lsdtt
