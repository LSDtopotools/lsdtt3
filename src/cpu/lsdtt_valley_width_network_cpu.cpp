// lsdtt_valley_width_network_cpu.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// See lsdtt_valley_width_network_cpu.hpp.
// -----------------------------------------------------------------------------

#include "lsdtt_valley_width_network_cpu.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <vector>

#include "../utils/lsdtt_grid_geodesic.hpp"
#include "../utils/lsdtt_grid_neighbours.hpp"
#include "../utils/lsdtt_raster_geometry.hpp"
#include "../utils/lsdtt_splines.hpp"
#include "lsdtt_proximity_cpu.hpp"

namespace lsdtt
{

namespace
{

static constexpr const auto& k_dr = grid::d8_rowmajor_row;
static constexpr const auto& k_dc = grid::d8_rowmajor_col;

// A valid cell on the grid border or next to nodata.
bool is_edge_cell(const lsdtt_flow_base& flow, std::int64_t f,
                  int n_rows, int n_cols)
{
    const int r = static_cast<int>(f / n_cols);
    const int c = static_cast<int>(f % n_cols);
    if (r == 0 || c == 0 || r == n_rows - 1 || c == n_cols - 1) return true;
    return flow.is_influenced_by_invalid(f);
}

double percentile(std::vector<double> v, double p)
{
    if (v.empty()) return 0.0;
    std::sort(v.begin(), v.end());
    const double pos = p * static_cast<double>(v.size() - 1);
    const std::size_t lo = static_cast<std::size_t>(std::floor(pos));
    const std::size_t hi = std::min(lo + 1, v.size() - 1);
    const double t = pos - static_cast<double>(lo);
    return v[lo] * (1.0 - t) + v[hi] * t;
}

double median(std::vector<double> v)
{
    return percentile(std::move(v), 0.5);
}

struct xy { double x, y; };

// Width line of a section: left side, the samples between the sides, right side.
std::vector<xy> width_line(const valley_cross_section& xs)
{
    std::vector<xy> pts;
    pts.push_back({xs.left_x, xs.left_y});
    for (const valley_cross_section_sample& sp : xs.samples)
        if (sp.offset_metres > xs.left_offset_metres &&
            sp.offset_metres < xs.right_offset_metres)
            pts.push_back({sp.x, sp.y});
    pts.push_back({xs.right_x, xs.right_y});
    return pts;
}

// The polyline with its first and last segments lengthened by `ext`.
std::vector<xy> extended(std::vector<xy> pts, double ext)
{
    if (pts.size() < 2 || ext <= 0.0) return pts;
    auto push_out = [ext](xy& end, const xy& prev) {
        const double dx = end.x - prev.x, dy = end.y - prev.y;
        const double len = std::hypot(dx, dy);
        if (len > 0.0) { end.x += dx / len * ext; end.y += dy / len * ext; }
    };
    push_out(pts.front(), pts[1]);
    push_out(pts.back(), pts[pts.size() - 2]);
    return pts;
}

struct bbox
{
    double x0 =  std::numeric_limits<double>::infinity();
    double y0 =  std::numeric_limits<double>::infinity();
    double x1 = -std::numeric_limits<double>::infinity();
    double y1 = -std::numeric_limits<double>::infinity();
    void add(const xy& p)
    {
        x0 = std::min(x0, p.x); y0 = std::min(y0, p.y);
        x1 = std::max(x1, p.x); y1 = std::max(y1, p.y);
    }
    bool overlaps(const bbox& o) const
    {
        return x0 <= o.x1 && o.x0 <= x1 && y0 <= o.y1 && o.y0 <= y1;
    }
};

double orient(const xy& a, const xy& b, const xy& c)
{
    return (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x);
}

bool on_segment(const xy& a, const xy& b, const xy& p)
{
    return std::min(a.x, b.x) <= p.x && p.x <= std::max(a.x, b.x) &&
           std::min(a.y, b.y) <= p.y && p.y <= std::max(a.y, b.y);
}

// Closed-segment intersection (touching counts).
bool segments_intersect(const xy& p1, const xy& p2, const xy& q1, const xy& q2)
{
    const double d1 = orient(q1, q2, p1);
    const double d2 = orient(q1, q2, p2);
    const double d3 = orient(p1, p2, q1);
    const double d4 = orient(p1, p2, q2);
    if (((d1 > 0 && d2 < 0) || (d1 < 0 && d2 > 0)) &&
        ((d3 > 0 && d4 < 0) || (d3 < 0 && d4 > 0)))
        return true;
    if (d1 == 0 && on_segment(q1, q2, p1)) return true;
    if (d2 == 0 && on_segment(q1, q2, p2)) return true;
    if (d3 == 0 && on_segment(p1, p2, q1)) return true;
    if (d4 == 0 && on_segment(p1, p2, q2)) return true;
    return false;
}

bool polylines_intersect(const std::vector<xy>& a, const std::vector<xy>& b)
{
    for (std::size_t i = 1; i < a.size(); ++i)
        for (std::size_t j = 1; j < b.size(); ++j)
            if (segments_intersect(a[i - 1], a[i], b[j - 1], b[j]))
                return true;
    return false;
}

} // namespace

network_centreline_choice parse_network_centreline_choice(const std::string& name)
{
    if (name == "medial_path") return network_centreline_choice::medial_path;
    if (name == "smoothed_channel")
        return network_centreline_choice::smoothed_channel;
    throw std::invalid_argument(
        "unknown network centreline choice '" + name +
        "'. Valid: medial_path, smoothed_channel.");
}

const char* to_string(network_section_exclusion e) noexcept
{
    switch (e)
    {
        case network_section_exclusion::none:                return "none";
        case network_section_exclusion::not_in_valley:       return "not_in_valley";
        case network_section_exclusion::side_not_found:      return "side_not_found";
        case network_section_exclusion::near_parent_channel: return "near_parent_channel";
        case network_section_exclusion::tributary_crossing:  return "tributary_crossing";
    }
    return "none";
}

valley_width_network_result build_valley_width_network(
    const lsdtt_raster_cpu&             dem,
    const lsdtt_flow_base&              flow,
    const std::vector<std::uint8_t>&    channel_mask,
    const lsdtt_raster_cpu&             valley_mask,
    const valley_width_network_request& request)
{
    const lsdtt_raster_meta& meta = dem.get_meta();
    const int n_rows = meta.n_rows;
    const int n_cols = meta.n_cols;
    const double cell = meta.cell_size;
    const std::size_t n = static_cast<std::size_t>(n_rows) * n_cols;

    if (channel_mask.size() != n || flow.full_grid_size() != n)
        throw std::invalid_argument(
            "build_valley_width_network: channel mask / flow / DEM size "
            "mismatch.");
    if (valley_mask.get_data().size() != n)
        throw std::invalid_argument(
            "build_valley_width_network: valley mask / DEM size mismatch.");

    const std::vector<float>&        elev  = dem.get_data();
    const std::vector<std::uint8_t>& dmask = dem.get_mask();
    const std::vector<float>&        vdata = valley_mask.get_data();
    const std::vector<std::uint8_t>& vmsk  = valley_mask.get_mask();

    std::vector<std::uint8_t> in_valley(n, 0);
    for (std::size_t i = 0; i < n; ++i)
        in_valley[i] = (dmask[i] && vmsk[i] && vdata[i] > 0.5f) ? 1 : 0;

    // -------------------------------------------------------------------------
    // Distance to the valley wall, once for the tile (0 outside the valley).
    // -------------------------------------------------------------------------
    std::vector<float> wd(n, 0.0f);
    {
        std::vector<float> wall_src(n, 0.0f);
        for (std::size_t i = 0; i < n; ++i)
            if (dmask[i] && !in_valley[i]) wall_src[i] = 1.0f;
        lsdtt_raster_cpu dist = euclidean_distance_transform(
            lsdtt_raster_cpu(std::move(wall_src),
                             std::vector<std::uint8_t>(dmask), meta));
        const std::vector<float>&        dd = dist.get_data();
        const std::vector<std::uint8_t>& dm = dist.get_mask();
        for (std::size_t i = 0; i < n; ++i)
            if (in_valley[i] && dm[i] && std::isfinite(dd[i])) wd[i] = dd[i];
    }

    // -------------------------------------------------------------------------
    // 1. Channel tree: heads traced longest-first, each stopping at a claimed
    //    pixel.
    // -------------------------------------------------------------------------
    std::vector<std::uint8_t> has_channel_donor(n, 0);
    std::size_t n_channel = 0;
    for (std::size_t i = 0; i < n; ++i)
    {
        if (!channel_mask[i] || !dmask[i]) continue;
        ++n_channel;
        const std::int64_t f = static_cast<std::int64_t>(i);
        if (flow.is_outlet_unchecked(f)) continue;
        const std::int64_t r = flow.receiver_of_unchecked(f);
        if (r >= 0 && r != f && channel_mask[static_cast<std::size_t>(r)])
            has_channel_donor[static_cast<std::size_t>(r)] = 1;
    }
    if (n_channel == 0)
        throw std::invalid_argument(
            "build_valley_width_network: the channel mask is empty.");

    const std::vector<float> fd = flow.get_distance_from_outlet_raster().data;
    std::vector<std::int64_t> heads;
    for (std::size_t i = 0; i < n; ++i)
        if (channel_mask[i] && dmask[i] && !has_channel_donor[i])
            heads.push_back(static_cast<std::int64_t>(i));
    std::sort(heads.begin(), heads.end(),
              [&](std::int64_t a, std::int64_t b) {
                  const float da = fd[static_cast<std::size_t>(a)];
                  const float db = fd[static_cast<std::size_t>(b)];
                  if (da != db) return da > db;
                  return a < b;
              });

    valley_width_network_result result;
    std::vector<valley_network_segment>& segs = result.segments;
    segs.reserve(heads.size());
    std::vector<std::int32_t> owner(n, -1);

    for (const std::int64_t h : heads)
    {
        valley_network_segment seg;
        seg.segment_id = static_cast<std::int32_t>(segs.size());
        seg.head_flat  = h;
        std::int64_t f = h;
        while (true)
        {
            const std::size_t fi = static_cast<std::size_t>(f);
            if (owner[fi] >= 0)
            {
                seg.confluence_flat   = f;
                seg.parent_segment_id = owner[fi];
                break;
            }
            owner[fi] = seg.segment_id;
            seg.channel_flats.push_back(f);
            if (flow.is_outlet_unchecked(f)) break;
            const std::int64_t g = flow.receiver_of_unchecked(f);
            if (g < 0 || g == f || !dmask[static_cast<std::size_t>(g)]) break;
            f = g;
        }
        if (seg.parent_segment_id >= 0)
        {
            const valley_network_segment& par =
                segs[static_cast<std::size_t>(seg.parent_segment_id)];
            seg.nesting_depth  = par.nesting_depth + 1;
            seg.drains_to_edge = par.drains_to_edge;
        }
        else if (!seg.channel_flats.empty())
        {
            seg.drains_to_edge =
                is_edge_cell(flow, seg.channel_flats.back(), n_rows, n_cols);
        }
        segs.push_back(std::move(seg));
    }

    // Catchment touches the edge: edge cells upslope of the segment's last own
    // pixel, less those in the 5 x 5 block round it that drain to it.
    {
        std::vector<float> edge_val(n, 0.0f);
        std::vector<std::uint8_t> edge_valid(n, 0);
        for (std::size_t i = 0; i < n; ++i)
        {
            if (!dmask[i]) continue;
            edge_valid[i] = 1;
            if (is_edge_cell(flow, static_cast<std::int64_t>(i), n_rows, n_cols))
                edge_val[i] = 1.0f;
        }
        const std::vector<float> n_edge_upslope =
            flow.get_upslope_sum_raster(edge_val, edge_valid, true).data;

        auto drains_to = [&](std::int64_t u, std::int64_t target) -> bool
        {
            for (int step = 0; step < 64; ++step)
            {
                if (u == target) return true;
                if (flow.is_outlet_unchecked(u)) return false;
                const std::int64_t g = flow.receiver_of_unchecked(u);
                if (g < 0 || g == u) return false;
                u = g;
            }
            return false;
        };

        for (valley_network_segment& seg : segs)
        {
            if (seg.channel_flats.empty()) continue;
            const std::int64_t last = seg.channel_flats.back();
            double count = n_edge_upslope[static_cast<std::size_t>(last)];
            const int lr = static_cast<int>(last / n_cols);
            const int lc = static_cast<int>(last % n_cols);
            for (int dr = -2; dr <= 2; ++dr)
                for (int dc = -2; dc <= 2; ++dc)
                {
                    const int rr = lr + dr, cc = lc + dc;
                    if (rr < 0 || cc < 0 || rr >= n_rows || cc >= n_cols) continue;
                    const std::size_t u =
                        static_cast<std::size_t>(rr) * n_cols + cc;
                    if (edge_val[u] > 0.5f &&
                        drains_to(static_cast<std::int64_t>(u), last))
                        count -= 1.0;
                }
            seg.catchment_touches_edge = count > 0.5;
        }
    }

    // -------------------------------------------------------------------------
    // 2. One centreline per segment.
    // -------------------------------------------------------------------------
    result.centrelines.channels_share_valley_mask = true;
    {
        std::vector<float> vm(n, 0.0f);
        for (std::size_t i = 0; i < n; ++i) vm[i] = in_valley[i] ? 1.0f : 0.0f;
        result.centrelines.valley_mask =
            lsdtt_raster_cpu(std::move(vm), std::vector<std::uint8_t>(dmask), meta);
    }

    auto flat_xy = [&](std::int64_t f) -> xy
    {
        const lsdtt_xy p = cell_center_from_row_col(
            meta, static_cast<int>(f / n_cols), static_cast<int>(f % n_cols));
        return {p.x, p.y};
    };

    for (valley_network_segment& seg : segs)
    {
        valley_channel_centreline rec;
        rec.source_id = seg.segment_id;
        rec.valley_break_derived = false;
        rec.valley_hand_break_metres = std::numeric_limits<double>::quiet_NaN();

        {
            std::vector<double> w;
            w.reserve(seg.channel_flats.size());
            for (const std::int64_t f : seg.channel_flats)
                w.push_back(static_cast<double>(wd[static_cast<std::size_t>(f)]));
            seg.valley_width_estimate_metres = 2.0 * percentile(w, 0.9);
        }

        if (static_cast<int>(seg.channel_flats.size()) < request.min_segment_pixels)
        {
            result.centrelines.channels.push_back(std::move(rec));
            continue;
        }

        // The channel pixels the centreline runs along: own pixels, then the
        // confluence pixel so a tributary reaches its parent.
        std::vector<std::int64_t> nodes = seg.channel_flats;
        if (seg.confluence_flat >= 0) nodes.push_back(seg.confluence_flat);
        const std::int64_t start  = nodes.front();
        const std::int64_t target = nodes.back();

        std::vector<std::int64_t> path;

        if (request.centreline == network_centreline_choice::medial_path)
        {
            double wmax_ch = 0.0;
            int r0 = n_rows, r1 = -1, c0 = n_cols, c1 = -1;
            for (const std::int64_t f : nodes)
            {
                wmax_ch = std::max(wmax_ch,
                                   static_cast<double>(wd[static_cast<std::size_t>(f)]));
                const int r = static_cast<int>(f / n_cols);
                const int c = static_cast<int>(f % n_cols);
                r0 = std::min(r0, r); r1 = std::max(r1, r);
                c0 = std::min(c0, c); c1 = std::max(c1, c);
            }
            const int margin = static_cast<int>(std::ceil(
                request.window_margin_multiple * wmax_ch / cell)) + 10;
            r0 = std::max(0, r0 - margin); r1 = std::min(n_rows - 1, r1 + margin);
            c0 = std::max(0, c0 - margin); c1 = std::min(n_cols - 1, c1 + margin);
            const int wr = r1 - r0 + 1;
            const int wc = c1 - c0 + 1;
            const std::size_t wn = static_cast<std::size_t>(wr) * wc;

            auto to_local = [&](std::int64_t f) -> std::int64_t
            {
                const int r = static_cast<int>(f / n_cols) - r0;
                const int c = static_cast<int>(f % n_cols) - c0;
                return static_cast<std::int64_t>(r) * wc + c;
            };
            auto to_global = [&](std::size_t l) -> std::size_t
            {
                const int r = static_cast<int>(l / wc) + r0;
                const int c = static_cast<int>(l % wc) + c0;
                return static_cast<std::size_t>(r) * n_cols + c;
            };

            std::vector<std::uint8_t> passable(wn, 0);
            std::vector<float>        lwd(wn, 0.0f);
            float wmax = 0.0f;
            for (std::size_t l = 0; l < wn; ++l)
            {
                const std::size_t g = to_global(l);
                if (in_valley[g])
                {
                    passable[l] = 1;
                    lwd[l] = wd[g];
                    wmax = std::max(wmax, wd[g]);
                }
            }
            for (const std::int64_t f : nodes)
                passable[static_cast<std::size_t>(to_local(f))] = 1;

            const std::int64_t ls = to_local(start);
            const std::int64_t lt = to_local(target);
            // Small positive length penalty (m), as in the single-channel
            // engine: every step costs something, so no absurd detours.
            const double eps = 0.5;
            const double diag = cell * grid::k_sqrt2_d;
            std::vector<double> cost(wn, std::numeric_limits<double>::infinity());
            std::vector<std::int64_t> prev(wn, -1);
            cost[static_cast<std::size_t>(ls)] = 0.0;
            grid_dijkstra_d8(wr, wc, cost, std::vector<std::int64_t>{ls},
                [&](std::size_t nf) { return passable[nf] != 0; },
                [&](std::size_t, std::size_t nf, int k) {
                    const double step =
                        (k_dr[k] != 0 && k_dc[k] != 0) ? diag : cell;
                    return step * (static_cast<double>(wmax - lwd[nf]) + eps); },
                [&](std::size_t from, std::size_t to) {
                    prev[to] = static_cast<std::int64_t>(from); },
                [&](std::int64_t l) { return l == lt; });

            if (std::isfinite(cost[static_cast<std::size_t>(lt)]))
            {
                std::int64_t cur = lt;
                while (cur != -1)
                {
                    path.push_back(static_cast<std::int64_t>(
                        to_global(static_cast<std::size_t>(cur))));
                    if (cur == ls) break;
                    cur = prev[static_cast<std::size_t>(cur)];
                }
                std::reverse(path.begin(), path.end());
            }
            else
            {
                seg.medial_path_failed = true;
                path = nodes;
            }
        }
        else   // smoothed_channel
        {
            const std::size_t nn = nodes.size();
            std::vector<double> s(nn, 0.0), x(nn), y(nn);
            for (std::size_t i = 0; i < nn; ++i)
            {
                const xy p = flat_xy(nodes[i]);
                x[i] = p.x; y[i] = p.y;
                if (i > 0) s[i] = s[i - 1] + std::hypot(x[i] - x[i - 1], y[i] - y[i - 1]);
            }
            const double sigma = std::max(
                2.0 * cell,
                request.smoothing_width_multiple * seg.valley_width_estimate_metres);

            // A Gaussian over an open curve pulls its ends inward by about
            // sigma, which collapses a short segment. Pad each end with the
            // curve point-reflected through that end (x' = 2 x_end - x): the
            // smoothed end then stays on the end pixel and the line keeps its
            // full length. Trimmed back after smoothing.
            const double pad_len = 3.0 * sigma;
            std::vector<double> ps, px, py;
            for (std::size_t k = nn - 1; k >= 1; --k)
            {
                if (s[k] > pad_len) continue;
                ps.push_back(-s[k]);
                px.push_back(2.0 * x[0] - x[k]);
                py.push_back(2.0 * y[0] - y[k]);
            }
            const std::size_t lead = ps.size();
            ps.insert(ps.end(), s.begin(), s.end());
            px.insert(px.end(), x.begin(), x.end());
            py.insert(py.end(), y.begin(), y.end());
            for (std::size_t k = nn - 1; k-- > 0;)
            {
                const double back = s[nn - 1] - s[k];
                if (back > pad_len) break;
                ps.push_back(s[nn - 1] + back);
                px.push_back(2.0 * x[nn - 1] - x[k]);
                py.push_back(2.0 * y[nn - 1] - y[k]);
            }
            const lsdtt_parametric_xy_series sm =
                gaussian_smooth_parametric_series(ps, px, py, sigma);
            const std::vector<double> ss(sm.s.begin() + lead, sm.s.begin() + lead + nn);
            const std::vector<double> sx(sm.x.begin() + lead, sm.x.begin() + lead + nn);
            const std::vector<double> sy(sm.y.begin() + lead, sm.y.begin() + lead + nn);
            const lsdtt_parametric_cubic_spline spline(ss, sx, sy);

            std::size_t j = 0;   // channel node nearest the current s
            const double s_end = s.back();
            for (double sq = 0.0; ; sq += cell)
            {
                const double sk = std::min(sq, s_end);
                while (j + 1 < nn &&
                       std::abs(s[j + 1] - sk) <= std::abs(s[j] - sk)) ++j;
                const lsdtt_parametric_spline_evaluation ev = spline.evaluate(sk);

                const int col = static_cast<int>(std::floor((ev.x - meta.x_min) / cell));
                const int row = static_cast<int>(std::floor((meta.y_max() - ev.y) / cell));
                const bool inside = row >= 0 && col >= 0 && row < n_rows && col < n_cols;
                const std::int64_t sf = inside
                    ? static_cast<std::int64_t>(row) * n_cols + col : -1;
                const std::int64_t node = nodes[j];

                std::int64_t chosen;
                if (sf >= 0 && in_valley[static_cast<std::size_t>(sf)])
                    chosen = sf;
                else if (in_valley[static_cast<std::size_t>(node)])
                {
                    // Off the mask while the channel is on it: the nearest
                    // in-mask cell to the spline point, so the line runs
                    // along the mask edge. It is no farther than the channel
                    // pixel, which is itself in the mask.
                    const xy np = flat_xy(node);
                    const double reach = std::hypot(ev.x - np.x, ev.y - np.y);
                    const int rad = static_cast<int>(std::ceil(reach / cell)) + 1;
                    const int cr = std::clamp(row, 0, n_rows - 1);
                    const int cc = std::clamp(col, 0, n_cols - 1);
                    chosen = node;
                    double best = std::numeric_limits<double>::infinity();
                    for (int rr = std::max(0, cr - rad);
                         rr <= std::min(n_rows - 1, cr + rad); ++rr)
                        for (int c2 = std::max(0, cc - rad);
                             c2 <= std::min(n_cols - 1, cc + rad); ++c2)
                        {
                            const std::size_t g =
                                static_cast<std::size_t>(rr) * n_cols + c2;
                            if (!in_valley[g]) continue;
                            const xy gp = flat_xy(static_cast<std::int64_t>(g));
                            const double d2 = (gp.x - ev.x) * (gp.x - ev.x) +
                                              (gp.y - ev.y) * (gp.y - ev.y);
                            if (d2 < best)
                            {
                                best = d2;
                                chosen = static_cast<std::int64_t>(g);
                            }
                        }
                    ++seg.n_samples_snapped_to_mask_edge;
                }
                else
                    chosen = node;

                if (path.empty() || path.back() != chosen) path.push_back(chosen);
                if (sk >= s_end) break;
            }
        }

        rec.routing_succeeded = !seg.medial_path_failed && path.size() >= 2;

        // Trough-elevation field kept for the shared outputs: the running
        // minimum of the DEM downstream along the centreline.
        float run = std::numeric_limits<float>::infinity();
        for (const std::int64_t f : path)
        {
            const std::size_t fi = static_cast<std::size_t>(f);
            if (dmask[fi] && elev[fi] < run) run = elev[fi];
            rec.centreline_trough_elevation.push_back(
                std::isfinite(run) ? run : elev[fi]);
        }

        // Baseline: pinned to the segment's own channel (and its confluence)
        // by the requested method. is_pin marks centreline nodes that are
        // themselves channel pixels.
        {
            std::vector<std::int64_t> own_sorted = nodes;
            std::sort(own_sorted.begin(), own_sorted.end());
            rec.centreline_is_pin.assign(path.size(), 0);
            for (std::size_t q = 0; q < path.size(); ++q)
                if (std::binary_search(own_sorted.begin(), own_sorted.end(), path[q]))
                    rec.centreline_is_pin[q] = 1;
            rec.centreline_baseline_elevation = compute_centreline_baseline(
                path, nodes, elev, meta, request.baseline);
        }
        rec.centreline_flats = std::move(path);
        result.centrelines.channels.push_back(std::move(rec));
    }

    return result;
}

std::vector<std::int64_t> match_points_to_channel(
    const std::vector<lsdtt_xy>&     points,
    const std::vector<double>&       along,
    double                           total_length,
    const std::vector<std::int64_t>& channel_flats,
    const lsdtt_raster_meta&         meta)
{
    std::vector<std::int64_t> match(points.size(), -1);
    if (points.empty() || channel_flats.empty()) return match;
    const int n_cols = meta.n_cols;
    std::vector<xy> cxy(channel_flats.size());
    for (std::size_t k = 0; k < channel_flats.size(); ++k)
    {
        const lsdtt_xy p = cell_center_from_row_col(
            meta, static_cast<int>(channel_flats[k] / n_cols),
            static_cast<int>(channel_flats[k] % n_cols));
        cxy[k] = {p.x, p.y};
    }
    const std::int64_t n = static_cast<std::int64_t>(channel_flats.size());
    const std::int64_t half = std::max<std::int64_t>(
        10, static_cast<std::int64_t>(std::ceil(
                3.0 * static_cast<double>(n) / static_cast<double>(points.size()))));
    std::int64_t prev = 0;
    for (std::size_t i = 0; i < points.size(); ++i)
    {
        const double f = total_length > 0.0
            ? std::clamp(along[i] / total_length, 0.0, 1.0) : 0.0;
        const std::int64_t expect =
            static_cast<std::int64_t>(std::llround(f * static_cast<double>(n - 1)));
        const std::int64_t lo = std::max(prev, expect - half);
        const std::int64_t hi = std::min(n - 1, std::max(lo, expect + half));
        std::int64_t best_k = lo;
        double best = std::numeric_limits<double>::infinity();
        for (std::int64_t k = lo; k <= hi; ++k)
        {
            const double d2 = (cxy[k].x - points[i].x) * (cxy[k].x - points[i].x) +
                              (cxy[k].y - points[i].y) * (cxy[k].y - points[i].y);
            if (d2 < best) { best = d2; best_k = k; }
        }
        match[i] = best_k;
        prev = best_k;
    }
    return match;
}

std::vector<std::int64_t> match_sections_to_channel(
    const valley_width_network_result& network,
    const valley_cross_section_result& sections,
    const lsdtt_raster_meta&           meta)
{
    const std::vector<valley_cross_section>& xs = sections.sections;
    const std::size_t n_seg = network.segments.size();
    const int n_cols = meta.n_cols;
    std::vector<std::int64_t> match(xs.size(), -1);

    std::vector<std::vector<std::size_t>> of_seg(n_seg);
    for (std::size_t i = 0; i < xs.size(); ++i)
    {
        const std::int32_t s = xs[i].source_id;
        if (s >= 0 && static_cast<std::size_t>(s) < n_seg)
            of_seg[static_cast<std::size_t>(s)].push_back(i);
    }

    for (std::size_t s = 0; s < n_seg; ++s)
    {
        std::vector<std::size_t>& st = of_seg[s];
        if (st.empty()) continue;
        std::sort(st.begin(), st.end(), [&](std::size_t a, std::size_t b) {
            return xs[a].along_distance_metres < xs[b].along_distance_metres; });

        // Centreline length, for each station's fractional position.
        const std::vector<std::int64_t>& cl = network.centrelines.channels[s].centreline_flats;
        double cl_len = 0.0;
        for (std::size_t k = 1; k < cl.size(); ++k)
        {
            const lsdtt_xy a = cell_center_from_row_col(
                meta, static_cast<int>(cl[k - 1] / n_cols), static_cast<int>(cl[k - 1] % n_cols));
            const lsdtt_xy b = cell_center_from_row_col(
                meta, static_cast<int>(cl[k] / n_cols), static_cast<int>(cl[k] % n_cols));
            cl_len += std::hypot(b.x - a.x, b.y - a.y);
        }
        std::vector<lsdtt_xy> pts;
        std::vector<double> along;
        for (const std::size_t i : st)
        {
            pts.push_back({xs[i].centre_x, xs[i].centre_y});
            along.push_back(xs[i].along_distance_metres);
        }
        const std::vector<std::int64_t> m = match_points_to_channel(
            pts, along, cl_len, network.segments[s].channel_flats, meta);
        for (std::size_t q = 0; q < st.size(); ++q) match[st[q]] = m[q];
    }
    return match;
}

std::vector<network_section_exclusion> apply_network_section_rules(
    const valley_width_network_result&   network,
    const valley_cross_section_result&   sections,
    const lsdtt_raster_meta&             meta,
    const network_section_rules_request& request)
{
    const std::vector<valley_cross_section>& xs = sections.sections;
    const std::size_t n_sec = xs.size();
    const std::size_t n_seg = network.segments.size();
    const int n_cols = meta.n_cols;

    // A station whose centre is outside the valley mask has no valley to
    // measure (its width is 0 or a sliver): the channel runs where the mask
    // has nothing. Then stations with a side not found.
    const std::vector<float>& vm = network.centrelines.valley_mask.get_data();
    std::vector<network_section_exclusion> reason(n_sec, network_section_exclusion::none);
    for (std::size_t i = 0; i < n_sec; ++i)
    {
        const int col = static_cast<int>(
            std::floor((xs[i].centre_x - meta.x_min) / meta.cell_size));
        const int row = static_cast<int>(
            std::floor((meta.y_max() - xs[i].centre_y) / meta.cell_size));
        const bool in_grid = row >= 0 && col >= 0 &&
                             row < meta.n_rows && col < meta.n_cols;
        const bool in_valley = in_grid && vm.size() ==
            static_cast<std::size_t>(meta.n_rows) * meta.n_cols &&
            vm[static_cast<std::size_t>(row) * n_cols + col] > 0.5f;
        if (!in_valley)
            reason[i] = network_section_exclusion::not_in_valley;
        else if (!(xs[i].left_found && xs[i].right_found))
            reason[i] = network_section_exclusion::side_not_found;
    }

    // Sections of each segment, upstream -> downstream.
    std::vector<std::vector<std::size_t>> of_seg(n_seg);
    for (std::size_t i = 0; i < n_sec; ++i)
    {
        const std::int32_t s = xs[i].source_id;
        if (s >= 0 && static_cast<std::size_t>(s) < n_seg)
            of_seg[static_cast<std::size_t>(s)].push_back(i);
    }
    for (std::vector<std::size_t>& v : of_seg)
        std::sort(v.begin(), v.end(), [&](std::size_t a, std::size_t b) {
            return xs[a].along_distance_metres < xs[b].along_distance_metres; });

    auto flat_xy = [&](std::int64_t f) -> xy
    {
        const lsdtt_xy p = cell_center_from_row_col(
            meta, static_cast<int>(f / n_cols), static_cast<int>(f % n_cols));
        return {p.x, p.y};
    };

    // Segments by nesting depth, so a parent's own stops come first.
    std::vector<std::size_t> order(n_seg);
    std::iota(order.begin(), order.end(), std::size_t{0});
    std::stable_sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) {
        return network.segments[a].nesting_depth < network.segments[b].nesting_depth; });

    // -------------------------------------------------------------------------
    // (a) Tributary stations stop near the parent's centreline.
    // -------------------------------------------------------------------------
    if (request.stop_tributaries_at_parent)
    {
        for (const std::size_t c : order)
        {
            const valley_network_segment& seg = network.segments[c];
            if (seg.parent_segment_id < 0 || seg.confluence_flat < 0) continue;
            const std::size_t p = static_cast<std::size_t>(seg.parent_segment_id);
            const std::vector<std::int64_t>& pline =
                network.centrelines.channels[p].centreline_flats;
            if (pline.empty() || of_seg[c].empty()) continue;

            // Parent's kept sections.
            std::vector<std::size_t> kept;
            for (const std::size_t i : of_seg[p])
                if (reason[i] == network_section_exclusion::none) kept.push_back(i);
            if (kept.empty()) continue;

            // The parent's local half-width on the side a point lies: the
            // median offset on that side over the parent's kept stations
            // within +/- window of the station nearest the point. Judged per
            // tributary station, not once at the confluence - the confluence is
            // on the parent's own channel, close to its centreline, so its side
            // says little about where the tributary's valley is. The window is
            // a multiple of the TRIBUTARY's valley width so it is always longer
            // than the mouth, whose sections are widened by the tributary and
            // must stay a minority of the median at any station spacing.
            const double window = std::max(
                request.parent_width_window_multiple *
                    seg.valley_width_estimate_metres,
                0.0);
            auto parent_half_width = [&](double px, double py) -> double
            {
                std::size_t k_near = 0;
                double d_near = std::numeric_limits<double>::infinity();
                for (std::size_t k = 0; k < kept.size(); ++k)
                {
                    const double d = std::hypot(xs[kept[k]].centre_x - px,
                                                xs[kept[k]].centre_y - py);
                    if (d < d_near) { d_near = d; k_near = k; }
                }
                const valley_cross_section& near = xs[kept[k_near]];
                const bool right_side =
                    (px - near.centre_x) * near.normal_dx +
                    (py - near.centre_y) * near.normal_dy > 0.0;
                const double a0 = near.along_distance_metres;
                std::vector<double> half;
                for (std::size_t k = 0; k < kept.size(); ++k)
                {
                    const valley_cross_section& s = xs[kept[k]];
                    const bool neighbour = k + 1 >= k_near && k <= k_near + 1;
                    if (!neighbour &&
                        std::abs(s.along_distance_metres - a0) > window)
                        continue;
                    half.push_back(right_side ? s.right_offset_metres
                                              : -s.left_offset_metres);
                }
                return median(half);
            };

            std::vector<xy> ppts;
            ppts.reserve(pline.size());
            for (const std::int64_t f : pline) ppts.push_back(flat_xy(f));

            bool stopped = false;
            for (const std::size_t i : of_seg[c])
            {
                if (!stopped)
                {
                    double dmin = std::numeric_limits<double>::infinity();
                    for (const xy& q : ppts)
                        dmin = std::min(dmin, std::hypot(xs[i].centre_x - q.x,
                                                         xs[i].centre_y - q.y));
                    const double h = parent_half_width(xs[i].centre_x,
                                                       xs[i].centre_y);
                    stopped = h > 0.0 &&
                              dmin <= request.tributary_stop_width_fraction * h;
                    // A width line that crosses the parent's centreline is
                    // measuring along the parent valley, wherever its centre is
                    // (an oblique section where the tributary's line bends
                    // into the confluence).
                    if (!stopped && reason[i] == network_section_exclusion::none)
                        stopped = polylines_intersect(width_line(xs[i]), ppts);
                }
                if (stopped && reason[i] != network_section_exclusion::not_in_valley)
                    reason[i] = network_section_exclusion::near_parent_channel;
            }
        }
    }

    // -------------------------------------------------------------------------
    // (b) Parent sections that cross a direct child's kept sections.
    // -------------------------------------------------------------------------
    if (request.remove_tributary_crossing_sections)
    {
        std::vector<std::vector<std::size_t>> children(n_seg);
        for (std::size_t c = 0; c < n_seg; ++c)
        {
            const std::int32_t p = network.segments[c].parent_segment_id;
            if (p >= 0) children[static_cast<std::size_t>(p)].push_back(c);
        }

        // Tributary width lines are lengthened by a cell at each end: a parent
        // section that has climbed into a tributary valley lies inside its
        // floor, so it crosses the tributary's width lines - except along the
        // floor's very edge, where it only touches their ends.
        std::vector<std::vector<xy>> line(n_sec), child_line(n_sec);
        std::vector<bbox> box(n_sec), child_box_of(n_sec);
        for (std::size_t i = 0; i < n_sec; ++i)
        {
            if (reason[i] == network_section_exclusion::side_not_found ||
                reason[i] == network_section_exclusion::not_in_valley) continue;
            line[i] = width_line(xs[i]);
            for (const xy& q : line[i]) box[i].add(q);
            child_line[i] = extended(line[i], meta.cell_size);
            for (const xy& q : child_line[i]) child_box_of[i].add(q);
        }

        // Decide from the reasons after (a) only, so the order parents are
        // visited in does not matter.
        const std::vector<network_section_exclusion> after_a = reason;
        for (std::size_t p = 0; p < n_seg; ++p)
        {
            if (children[p].empty()) continue;
            std::vector<std::size_t> child_kept;
            bbox child_box;
            for (const std::size_t c : children[p])
                for (const std::size_t t : of_seg[c])
                    if (after_a[t] == network_section_exclusion::none)
                    {
                        child_kept.push_back(t);
                        child_box.x0 = std::min(child_box.x0, child_box_of[t].x0);
                        child_box.y0 = std::min(child_box.y0, child_box_of[t].y0);
                        child_box.x1 = std::max(child_box.x1, child_box_of[t].x1);
                        child_box.y1 = std::max(child_box.y1, child_box_of[t].y1);
                    }
            if (child_kept.empty()) continue;

            for (const std::size_t s : of_seg[p])
            {
                if (after_a[s] != network_section_exclusion::none) continue;
                if (!box[s].overlaps(child_box)) continue;
                for (const std::size_t t : child_kept)
                {
                    if (!box[s].overlaps(child_box_of[t])) continue;
                    if (polylines_intersect(line[s], child_line[t]))
                    {
                        reason[s] = network_section_exclusion::tributary_crossing;
                        break;
                    }
                }
            }
        }
    }

    return reason;
}

} // namespace lsdtt
