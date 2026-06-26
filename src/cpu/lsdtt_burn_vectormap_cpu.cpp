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
#include <unordered_map>
#include <unordered_set>
#include <variant>

#include "../core/lsdtt_data_loader.hpp"
#include "../utils/lsdtt_polygon_geometry.hpp"
#include "lsdtt_least_cost_path.hpp"
#include "lsdtt_vector_to_junction_network.hpp"

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

// Cumulative arc-length at each vertex of a polyline (cum[0] = 0).
std::vector<double> cumulative_lengths(const std::vector<lsdtt_xy>& pts)
{
    std::vector<double> cum(pts.size(), 0.0);
    for (std::size_t i = 1; i < pts.size(); ++i)
        cum[i] = cum[i - 1] + std::hypot(pts[i].x - pts[i - 1].x,
                                         pts[i].y - pts[i - 1].y);
    return cum;
}

// Point on a polyline at arc-length s (clamped to [0, total]).
lsdtt_xy point_at_arclength(const std::vector<lsdtt_xy>& pts,
                            const std::vector<double>& cum, double s)
{
    if (pts.empty()) return {0.0, 0.0};
    if (s <= 0.0) return pts.front();
    if (s >= cum.back()) return pts.back();
    std::size_t i = 1;
    while (i < cum.size() && cum[i] < s) ++i;
    const double seg = cum[i] - cum[i - 1];
    const double t = seg > 0.0 ? (s - cum[i - 1]) / seg : 0.0;
    return { pts[i - 1].x + t * (pts[i].x - pts[i - 1].x),
             pts[i - 1].y + t * (pts[i].y - pts[i - 1].y) };
}

// Sub-polyline between arc-lengths s0 < s1.
std::vector<lsdtt_xy> subline(const std::vector<lsdtt_xy>& pts,
                              const std::vector<double>& cum,
                              double s0, double s1)
{
    std::vector<lsdtt_xy> out;
    out.push_back(point_at_arclength(pts, cum, s0));
    for (std::size_t i = 0; i < pts.size(); ++i)
        if (cum[i] > s0 + 1e-9 && cum[i] < s1 - 1e-9)
            out.push_back(pts[i]);
    out.push_back(point_at_arclength(pts, cum, s1));
    return out;
}

} // namespace

// ---------------------------------------------------------------------------
// burn_vectormap_trace_divergence
// ---------------------------------------------------------------------------
vectormap_trace_burn_result burn_vectormap_trace_divergence(
    lsdtt_raster_cpu&                   dem,
    const lsdtt_flow_base&              flow,
    const std::string&                  vector_path,
    const vectormap_trace_burn_options& options)
{
    vectormap_trace_burn_result result;

    const lsdtt_raster_meta& meta = dem.get_meta();
    const int    n_rows = meta.n_rows;
    const int    n_cols = meta.n_cols;
    const double cs     = meta.cell_size;
    const double seg_len = std::max(cs, options.segment_length_m);
    const double tol     = options.max_divergence_m;
    const double step    = (options.sample_step_m > 0.0)
                               ? options.sample_step_m : 0.5 * cs;
    const int max_steps = options.max_trace_steps > 0
                              ? options.max_trace_steps
                              : 10 * std::max(n_rows, n_cols);

    const lsdtt_array2d_f&  zdata = dem.get_data_array();
    const lsdtt_array2d_u8& mdata = dem.get_mask_array();

    // Elevation at a point, but only on a VALID (in-mask) cell; NaN otherwise.
    // Using the mask is essential for orientation: a vector vertex that falls on
    // a nodata cell returns the raster's nodata sentinel (a large finite value),
    // which silently fools an endpoint-elevation comparison.
    auto valid_elev_at_xy = [&](const lsdtt_xy& p, double& z) -> bool
    {
        int r, c;
        if (!xy_to_rc(meta, p.x, p.y, r, c)) return false;
        if (mdata(r, c) == 0) return false;
        z = static_cast<double>(zdata(r, c));
        return std::isfinite(z);
    };

    // ---- Read + clip vector lines to the DEM extent -----------------------
    // Clip to the in-DEM linework using the shared boundary handling (splits
    // boundary meanders into in-DEM pieces, drops short stubs and shallow
    // pass-in/pass-out meanders) so the burn never processes off-DEM lines.
    const std::vector<std::vector<lsdtt_xy>> clipped_lines =
        clip_vector_lines_to_dem(vector_path, dem);

    // ---- Build all vector lines, oriented downhill ------------------------
    // Each line's points are oriented so arc-length increases downstream (first
    // and last VALID, in-mask vertices compared: VML lines run into the DEM
    // nodata footprint, so a nodata endpoint would otherwise reverse the line).
    struct oriented_line
    {
        std::vector<lsdtt_xy> pts;
        std::vector<double>   cum;
        double                total = 0.0;
    };
    std::vector<oriented_line> all_lines;
    all_lines.reserve(clipped_lines.size());
    for (const std::vector<lsdtt_xy>& line_pts : clipped_lines)
    {
        if (line_pts.size() < 2) continue;
        oriented_line ol;
        ol.pts = line_pts;
        double z_first = 0.0, z_last = 0.0;
        bool ok_first = false, ok_last = false;
        for (std::size_t a = 0; a < ol.pts.size() && !ok_first; ++a)
            ok_first = valid_elev_at_xy(ol.pts[a], z_first);
        for (std::size_t a = ol.pts.size(); a-- > 0 && !ok_last; )
            ok_last = valid_elev_at_xy(ol.pts[a], z_last);
        if (ok_first && ok_last && z_first < z_last)
            std::reverse(ol.pts.begin(), ol.pts.end());
        ol.cum = cumulative_lengths(ol.pts);
        ol.total = ol.cum.back();
        if (ol.total < cs) continue;
        all_lines.push_back(std::move(ol));
    }

    // Nearest downstream continuation of line `li`: another line whose UPSTREAM
    // end (arc-length 0) lies within `connect_tol_m` of line li's downstream end.
    // Lets a burned run keep descending across the gap where one VML fragment
    // ends at an obstruction and the channel resumes as a separate fragment.
    const double connect_tol_m = (options.channel_anchor_radius_m > 0.0)
                                     ? options.channel_anchor_radius_m : 30.0;
    auto find_continuation = [&](int li, const lsdtt_xy& endp) -> int
    {
        int best = -1;
        double best_d = connect_tol_m;
        for (int j = 0; j < static_cast<int>(all_lines.size()); ++j)
        {
            if (j == li) continue;
            const lsdtt_xy& up = all_lines[j].pts.front();
            const double d = std::hypot(up.x - endp.x, up.y - endp.y);
            if (d < best_d) { best_d = d; best = j; }
        }
        return best;
    };

    // Rasterise a straight cell line a->b (Bresenham), inclusive.
    auto line_cells = [&](std::int64_t a, std::int64_t b)
        -> std::vector<std::int64_t>
    {
        int r0 = static_cast<int>(a / n_cols), c0 = static_cast<int>(a % n_cols);
        const int r1 = static_cast<int>(b / n_cols), c1 = static_cast<int>(b % n_cols);
        std::vector<std::int64_t> out;
        const int dr = std::abs(r1 - r0), dc = std::abs(c1 - c0);
        const int sr = r0 < r1 ? 1 : -1, sc = c0 < c1 ? 1 : -1;
        int err = dc - dr;
        while (true)
        {
            out.push_back(static_cast<std::int64_t>(r0) * n_cols + c0);
            if (r0 == r1 && c0 == c1) break;
            const int e2 = 2 * err;
            if (e2 > -dr) { err -= dr; c0 += sc; }
            if (e2 <  dc) { err += dc; r0 += sr; }
        }
        return out;
    };

    // Lateral fallback: nearest cell LOWER than z_ref within the anchor radius,
    // preferring the downstream half-plane (dot with the run heading > 0) so we
    // do not connect into a wrong-way ditch. Used when following the vector down
    // does not reach lower ground but a lower channel sits beside the run end.
    const int anchor_rad_cells =
        (options.channel_anchor_radius_m > 0.0)
            ? std::max(1, static_cast<int>(std::ceil(
                  options.channel_anchor_radius_m / cs)))
            : 0;
    auto find_lower_nearby = [&](int r0, int c0, double z_ref,
                                 double dirx, double diry) -> std::int64_t
    {
        if (anchor_rad_cells <= 0) return -1;
        const double dn = std::hypot(dirx, diry);
        const bool have_dir = dn > 1.0e-9;
        if (have_dir) { dirx /= dn; diry /= dn; }
        std::int64_t best = -1; double best_z = z_ref - 0.1;
        for (int dr = -anchor_rad_cells; dr <= anchor_rad_cells; ++dr)
            for (int dc = -anchor_rad_cells; dc <= anchor_rad_cells; ++dc)
            {
                if (dr == 0 && dc == 0) continue;
                const int r = r0 + dr, c = c0 + dc;
                if (r < 0 || r >= n_rows || c < 0 || c >= n_cols) continue;
                if (mdata(r, c) == 0) continue;
                if (have_dir &&
                    static_cast<double>(dc) * dirx +
                    static_cast<double>(-dr) * diry <= 0.0) continue;
                const double z = static_cast<double>(zdata(r, c));
                if (z < best_z) { best_z = z;
                    best = static_cast<std::int64_t>(r) * n_cols + c; }
            }
        return best;
    };

    // Build a burn path that descends from arc-length s0 on line `line0` and
    // keeps walking DOWN the vector — continuing along the line, and across into
    // the next connected fragment when one ends at an obstruction — until it
    // reaches a cell genuinely LOWER than the run start (so the trench daylights
    // below the obstruction instead of perching on its crest). `found_lower` is
    // set when such an endpoint is reached. The path stops at the first lower
    // cell, is capped at `max_walk_m`, and never revisits a line.
    const double max_walk_m = 8.0 * std::max(seg_len, 50.0);
    auto build_descending_path =
        [&](int line0, double s0, double s1, double z_start, bool& found_lower)
        -> std::vector<std::int64_t>
    {
        found_lower = false;
        std::vector<std::int64_t> path;
        auto push_cell = [&](const lsdtt_xy& p, double& z_out) -> bool
        {
            int r, c;
            if (!xy_to_rc(meta, p.x, p.y, r, c)) { z_out = 0.0; return false; }
            const std::int64_t f = static_cast<std::int64_t>(r) * n_cols + c;
            if (path.empty() || path.back() != f) path.push_back(f);
            z_out = static_cast<double>(zdata(r, c));
            return mdata(r, c) != 0;
        };

        std::vector<char> visited(all_lines.size(), 0);
        int    cur = line0;
        double cur_s = s0;
        double walked = 0.0;
        while (cur >= 0 && !visited[static_cast<std::size_t>(cur)]
               && walked < max_walk_m)
        {
            visited[static_cast<std::size_t>(cur)] = 1;
            const oriented_line& Ln = all_lines[static_cast<std::size_t>(cur)];
            for (double s = cur_s; s <= Ln.total + 1e-6; s += step)
            {
                const double ss = std::min(s, Ln.total);
                double z = 0.0;
                const bool valid =
                    push_cell(point_at_arclength(Ln.pts, Ln.cum, ss), z);
                walked += step;
                // Burn the whole divergent run first; only LOOK for a lower
                // daylight endpoint once we are past the run end s1 (on the
                // starting line). Otherwise a normally descending run would stop
                // a few cells in and never bridge the obstruction.
                const bool past_run = (cur != line0) || (ss >= s1 - 1e-6);
                const double need = std::max(0.1, options.burn_min_slope * walked);
                if (past_run && valid && z < z_start - need)
                {
                    found_lower = true;
                    return path;
                }
                if (walked >= max_walk_m) return path;
            }
            // Reached this line's downstream end; try to continue into the next
            // connected fragment, bridging the gap between them with a straight
            // cell line so the burned trench stays contiguous.
            const lsdtt_xy endp = Ln.pts.back();
            const int nxt = find_continuation(cur, endp);
            if (nxt < 0) break;
            int fr, fc;
            if (!path.empty() &&
                xy_to_rc(meta,
                         all_lines[static_cast<std::size_t>(nxt)].pts.front().x,
                         all_lines[static_cast<std::size_t>(nxt)].pts.front().y,
                         fr, fc))
            {
                const std::int64_t b = static_cast<std::int64_t>(fr) * n_cols + fc;
                for (std::int64_t f : line_cells(path.back(), b))
                    if (path.back() != f) path.push_back(f);
            }
            cur = nxt;
            cur_s = 0.0;
        }
        return path;
    };

    // Carve a divergent run starting at arc-length s0 on line `line_idx`, walking
    // down the vector network to a lower daylight endpoint. `anchored` is set when
    // a genuinely lower endpoint was reached. Returns true when a burnable path
    // was produced (skipped when dry_run).
    auto burn_run = [&](int line_idx, double s0, double s1, bool& anchored) -> bool
    {
        anchored = false;
        const oriented_line& Ln = all_lines[static_cast<std::size_t>(line_idx)];
        int sr, sc;
        if (!xy_to_rc(meta, point_at_arclength(Ln.pts, Ln.cum, s0).x,
                            point_at_arclength(Ln.pts, Ln.cum, s0).y, sr, sc))
            return false;
        const double z_start = static_cast<double>(zdata(sr, sc));

        std::vector<std::int64_t> path =
            build_descending_path(line_idx, s0, s1, z_start, anchored);
        if (path.size() < 2) return false;

        // Always pull the trench end to the lowest nearby cell (lower than the
        // current path end), so it daylights in the real channel whether that
        // resumes ALONG the vector (handled by the walk above) or BESIDE the run
        // end (a channel off to the side, not along a vector continuation).
        {
            const lsdtt_xy a = point_at_arclength(Ln.pts, Ln.cum, s0);
            const lsdtt_xy b = point_at_arclength(Ln.pts, Ln.cum,
                                                  std::min(s1, Ln.total));
            const std::int64_t end = path.back();
            const double z_end = static_cast<double>(
                zdata(static_cast<int>(end / n_cols),
                      static_cast<int>(end % n_cols)));
            const std::int64_t tgt = find_lower_nearby(
                static_cast<int>(end / n_cols), static_cast<int>(end % n_cols),
                z_end, b.x - a.x, b.y - a.y);
            if (tgt >= 0 && tgt != end)
            {
                for (std::int64_t f : line_cells(end, tgt))
                    if (path.back() != f) path.push_back(f);
                anchored = true;
            }
        }
        if (path.size() < 2) return false;
        if (!options.dry_run)
            burn_capturing_path(dem, path, options.burn_min_slope,
                                /*capture_margin=*/1.0e-3, options.max_burn_depth_m);
        return true;
    };

    for (int line_index = 0; line_index < static_cast<int>(all_lines.size());
         ++line_index)
    {
        const std::vector<lsdtt_xy>& pts = all_lines[line_index].pts;
        const std::vector<double>&   cum = all_lines[line_index].cum;
        const double                 total = all_lines[line_index].total;
        ++result.n_lines;

        const int m = std::max(1, static_cast<int>(std::ceil(total / seg_len)));
        result.n_segments += m;

        // Outcome of tracing the D8 flow from the top of segment k: how far
        // (segment index) it follows the vector within tol, the largest
        // trace-to-vector separation observed, and optionally the trace path.
        // Why a trace stopped (diagnostic): which terminating condition fired.
        // 0 advanced (reach>k), 1 diverged (>tol), 2 reached line end,
        // 3 hit DEM outlet, 4 no receiver, 5 max steps, 6 off-DEM top.
        enum trace_stop { STOP_ADVANCED=0, STOP_DIVERGED=1, STOP_END=2,
                          STOP_OUTLET=3, STOP_NORECV=4, STOP_MAXSTEPS=5,
                          STOP_OFFDEM=6 };

        struct trace_outcome
        {
            int    reach = 0;
            double max_divergence = 0.0;
            int    stop = STOP_MAXSTEPS;
            std::vector<lsdtt_xy> path;
        };

        auto do_trace = [&](int k, bool want_path) -> trace_outcome
        {
            trace_outcome oc;
            oc.reach = k;
            const double s_top = k * seg_len;
            const lsdtt_xy ptop = point_at_arclength(pts, cum, s_top);
            int row, col;
            if (!xy_to_rc(meta, ptop.x, ptop.y, row, col))
            {
                oc.reach = k + 1;   // off-DEM: treat as cleared
                oc.stop  = STOP_OFFDEM;
                return oc;
            }
            std::int64_t cur = flow.flat_index_unchecked(row, col);
            double max_s = s_top;
            for (int s = 0; s < max_steps; ++s)
            {
                const int rr = flow.row_from_flat_unchecked(cur);
                const int cc = flow.col_from_flat_unchecked(cur);
                const lsdtt_xy here = rc_to_centre(meta, rr, cc);
                if (want_path) oc.path.push_back(here);
                const lsdtt_nearest_polyline_point np =
                    nearest_point_on_polyline(here, pts);
                if (np.distance > oc.max_divergence) oc.max_divergence = np.distance;
                if (np.distance > tol) { oc.stop = STOP_DIVERGED; break; }
                const double s_here =
                    cum[np.segment_index] +
                    np.t * (cum[np.segment_index + 1] - cum[np.segment_index]);
                if (s_here > max_s) max_s = s_here;
                if (max_s >= total - 0.5 * cs) { oc.stop = STOP_END; break; }
                if (flow.is_outlet_unchecked(cur)) { oc.stop = STOP_OUTLET; break; }
                const std::int64_t nxt = flow.receiver_of_unchecked(cur);
                if (nxt == cur || nxt < 0) { oc.stop = STOP_NORECV; break; }
                cur = nxt;
            }
            oc.reach = static_cast<int>(max_s / seg_len + 1e-6);
            if (oc.reach > k) oc.stop = STOP_ADVANCED;
            return oc;
        };

        // Per-segment classification for this line (collect_diagnostics).
        std::vector<std::uint8_t> seg_status;   // 0 followed, 1 burned/diverging
        std::vector<double>       seg_div;
        std::vector<std::uint8_t> seg_reason;   // trace_stop code
        if (options.collect_diagnostics)
        {
            seg_status.assign(static_cast<std::size_t>(m), 0);
            seg_div.assign(static_cast<std::size_t>(m), 0.0);
            seg_reason.assign(static_cast<std::size_t>(m), STOP_ADVANCED);
        }

        int i = 0;
        while (i < m)
        {
            const trace_outcome oc = do_trace(i, false);

            // Only a genuine LATERAL divergence (trace strayed > tol from the
            // vector) is a problem to burn. A trace that merely fails to advance
            // a full segment because it reached the line end, ran off the valid
            // DEM (a base-level outlet), or hit a sink — all while staying within
            // tolerance — has FOLLOWED the channel correctly; there is nothing to
            // carve. Treating those terminations as divergence (the old reach<=k
            // test) painted whole valleys that drain to the data edge as failures
            // and carved spurious trenches under them.
            if (oc.stop != STOP_DIVERGED)
            {
                const int j = std::max(oc.reach, i + 1); // always make progress
                if (options.collect_diagnostics)
                    for (int t = i; t < j && t < m; ++t)
                    {
                        seg_status[static_cast<std::size_t>(t)] = 0;
                        seg_div[static_cast<std::size_t>(t)] =
                            std::max(seg_div[static_cast<std::size_t>(t)],
                                     oc.max_divergence);
                        seg_reason[static_cast<std::size_t>(t)] =
                            static_cast<std::uint8_t>(oc.stop);
                    }
                i = j;
                continue;
            }

            // Divergent run: extend over consecutive divergent segments so the
            // whole reach is burned as ONE continuous, draining channel ending
            // where flow rejoins the vector.
            int k = i;
            double run_div = oc.max_divergence;
            while (k < m)
            {
                const trace_outcome ock = do_trace(k, false);
                if (ock.stop != STOP_DIVERGED) break;
                run_div = std::max(run_div, ock.max_divergence);
                if (options.collect_diagnostics)
                    seg_reason[static_cast<std::size_t>(k)] =
                        static_cast<std::uint8_t>(ock.stop);
                ++k;
            }
            if (k == i) ++k;                             // safety: force progress

            const double s0 = i * seg_len;
            const double s1 = std::min(static_cast<double>(k) * seg_len, total);
            bool anchored = false;
            if (burn_run(line_index, s0, s1, anchored))
            {
                result.n_segments_burned += (k - i);
                result.total_burn_length_m += (s1 - s0);
                if (anchored) ++result.n_runs_anchored;
                if (!options.dry_run)
                {
                    vectormap_burn_gap g;
                    g.path_xy = subline(pts, cum, s0, s1);
                    g.length_m = s1 - s0;
                    g.burned = true;
                    result.burned.push_back(std::move(g));
                }
            }
            if (options.collect_diagnostics)
                for (int t = i; t < k && t < m; ++t)
                {
                    seg_status[static_cast<std::size_t>(t)] = 1;
                    seg_div[static_cast<std::size_t>(t)] =
                        std::max(seg_div[static_cast<std::size_t>(t)], run_div);
                }
            i = k;
        }

        // Emit per-segment diagnostic records for this line.
        if (options.collect_diagnostics)
        {
            for (int t = 0; t < m; ++t)
            {
                vectormap_trace_segment vs;
                vs.line_index    = line_index;
                vs.segment_index = t;
                vs.divergence_m  = seg_div[static_cast<std::size_t>(t)];
                vs.status        = seg_status[static_cast<std::size_t>(t)];
                vs.stop_reason   = seg_reason[static_cast<std::size_t>(t)];
                const double s0 = t * seg_len;
                const double s1 = std::min(static_cast<double>(t + 1) * seg_len,
                                           total);
                vs.xy = subline(pts, cum, s0, s1);
                result.segments.push_back(std::move(vs));
            }
        }

        // Emit sampled D8 traces (flow-trace overlay).
        if (options.collect_traces)
        {
            const double spacing = options.trace_sample_spacing_m > 0.0
                                       ? options.trace_sample_spacing_m : 200.0;
            const int step_seg =
                std::max(1, static_cast<int>(std::lround(spacing / seg_len)));
            for (int t = 0; t < m; t += step_seg)
            {
                const trace_outcome oc = do_trace(t, true);
                if (oc.path.size() < 2) continue;
                vectormap_trace_polyline tp;
                tp.line_index       = line_index;
                tp.segment_index    = t;
                tp.xy               = oc.path;
                tp.max_divergence_m = oc.max_divergence;
                result.traces.push_back(std::move(tp));
            }
        }
    }

    return result;
}

// ---------------------------------------------------------------------------
// burn_vectormap_source_paths
// ---------------------------------------------------------------------------
vectormap_source_burn_result burn_vectormap_source_paths(
    lsdtt_raster_cpu&                         dem,
    const lsdtt_flow_base&                    flow,
    const std::vector<std::vector<lsdtt_xy>>& source_paths,
    const vectormap_trace_burn_options&       options)
{
    vectormap_source_burn_result result;

    const lsdtt_raster_meta& meta   = dem.get_meta();
    const int    n_rows = meta.n_rows;
    const int    n_cols = meta.n_cols;
    const double cs     = meta.cell_size;
    const double tol    = options.max_divergence_m;
    const double step   = (options.sample_step_m > 0.0)
                              ? options.sample_step_m : 0.5 * cs;
    const int max_steps = options.max_trace_steps > 0
                              ? options.max_trace_steps
                              : 10 * std::max(n_rows, n_cols);

    const lsdtt_array2d_f&  zdata = dem.get_data_array();
    const lsdtt_array2d_u8& mdata = dem.get_mask_array();
    std::vector<float>&     elevw = dem.get_data();

    // Sloping-cut carve: lower every cell on `sub` to a SINGLE straight
    // descending line from the first cell's elevation to the last cell's
    // elevation. The caller chooses `sub` to span an obstruction (berm/road) from
    // a point above it to where the channel resumes below it (the lower-convex-
    // hull anchor), so by construction the natural surface lies above this line in
    // between: the cut slices cleanly through the obstruction with a constant
    // grade and daylights at both ends — no step, no flat. Cells are only lowered.
    auto burn_subpath = [&](const std::vector<lsdtt_xy>& sub) -> bool
    {
        std::vector<std::int64_t> path;
        for (const lsdtt_xy& p : sample_polyline(sub, step))
        {
            int r, c;
            if (!xy_to_rc(meta, p.x, p.y, r, c)) continue;
            const std::int64_t f = static_cast<std::int64_t>(r) * n_cols + c;
            if (path.empty() || path.back() != f) path.push_back(f);
        }
        if (path.size() < 2) return false;

        std::vector<double> d(path.size(), 0.0);
        for (std::size_t k = 1; k < path.size(); ++k)
        {
            const int r0 = static_cast<int>(path[k-1] / n_cols);
            const int c0 = static_cast<int>(path[k-1] % n_cols);
            const int r1 = static_cast<int>(path[k]   / n_cols);
            const int c1 = static_cast<int>(path[k]   % n_cols);
            d[k] = d[k-1] + std::hypot((r1 - r0) * cs, (c1 - c0) * cs);
        }
        const double tot = d.back();
        if (tot <= 0.0) return false;

        const double zt = static_cast<double>(elevw[path.front()]);
        const double zb = static_cast<double>(elevw[path.back()]);
        if (zt <= zb) return false;   // not descending: nothing to dig

        // Pure straight sloping cut: lower each cell to the line zt -> zb. The
        // caller passes a path SNAPPED to the local thalweg and truncated at the
        // lower-convex-hull outlet, so the natural surface lies above this line in
        // between -- a clean constant-grade cut through the obstruction, daylighting
        // at both ends, with no step and no flat. Cells are only lowered.
        for (std::size_t k = 0; k < path.size(); ++k)
        {
            const double zl = zt - (zt - zb) * d[k] / tot;
            if (static_cast<double>(elevw[path[k]]) > zl)
                elevw[path[k]] = static_cast<float>(zl);
        }
        return true;
    };

    // Deep forcing carve: lower the cells of `sub` to a strictly-descending
    // trench `depth` m below the natural surface, daylighting at the (low) outlet.
    // Used to FORCE flow onto the vector through a real obstruction (a road berm
    // with a culvert) where the LiDAR's competing route is a drainage ditch that
    // must be beaten. A FAR low outlet keeps the interior at natural-minus-depth
    // (a clean slope parallel to the surface) rather than a flat. (See
    // burn_relative_depth_path.)
    auto burn_deep = [&](const std::vector<lsdtt_xy>& sub, double depth) -> bool
    {
        std::vector<std::int64_t> path;
        for (const lsdtt_xy& p : sample_polyline(sub, step))
        {
            int r, c;
            if (!xy_to_rc(meta, p.x, p.y, r, c)) continue;
            const std::int64_t f = static_cast<std::int64_t>(r) * n_cols + c;
            if (path.empty() || path.back() != f) path.push_back(f);
        }
        if (path.size() < 2) return false;
        burn_relative_depth_path(dem, path, depth, options.burn_min_slope);
        return true;
    };

    // Snap a point to the lowest valid DEM cell within `radius` m -- the local
    // thalweg. So a vector drawn a few metres off the valley bottom is burned INTO
    // the thalweg, while a point over an obstruction (no low cell nearby) barely
    // moves, leaving the sloping cut to dig through it. Returns the cell centre and
    // sets z_out to its elevation.
    auto snap_to_thalweg = [&](const lsdtt_xy& p, double radius,
                               double& z_out) -> lsdtt_xy
    {
        int r, c;
        if (!xy_to_rc(meta, p.x, p.y, r, c))
        { z_out = std::numeric_limits<double>::max(); return p; }
        const int rad = std::max(0, static_cast<int>(std::ceil(radius / cs)));
        const double r2 = radius * radius;
        double best = std::numeric_limits<double>::max();
        int br = r, bc = c;
        for (int dr = -rad; dr <= rad; ++dr)
            for (int dc = -rad; dc <= rad; ++dc)
            {
                const int rr = r + dr, cc = c + dc;
                if (rr < 0 || rr >= n_rows || cc < 0 || cc >= n_cols) continue;
                if (mdata(rr, cc) == 0) continue;
                if (static_cast<double>(dr*dr + dc*dc) * cs * cs > r2) continue;
                const double z = static_cast<double>(zdata(rr, cc));
                if (z < best) { best = z; br = rr; bc = cc; }
            }
        z_out = best;
        return rc_to_centre(meta, br, bc);
    };

    const double seg_len = std::max(cs, options.segment_length_m);

    // Trace the D8 path from arc-length s_top and report how far along the path
    // it follows (within tol). On divergence, `escape_z` is set to the elevation
    // of the cell the flow strayed TO (lower than the path, the competing route)
    // — the carve must reach a path pixel below this to actually capture the
    // flow, not merely below the divergence point. `diverged` distinguishes a
    // genuine stray from reaching the path end / a base-level outlet.
    // `trace_out`/`maxdiv_out` (optional, for _algoviz): when non-null, record the
    // D8 cell-centre path walked and the maximum trace-to-path separation seen.
    auto trace_reach_s = [&](const std::vector<lsdtt_xy>& path,
                             const std::vector<double>& cum, double total,
                             double s_top, double& escape_z, bool& diverged,
                             std::vector<lsdtt_xy>* trace_out = nullptr,
                             double* maxdiv_out = nullptr)
        -> double
    {
        escape_z = std::numeric_limits<double>::max();
        diverged = false;
        if (maxdiv_out) *maxdiv_out = 0.0;
        if (trace_out)  trace_out->clear();
        int row, col;
        const lsdtt_xy ptop = point_at_arclength(path, cum, s_top);
        if (!xy_to_rc(meta, ptop.x, ptop.y, row, col)) return total; // off-DEM: clear
        std::int64_t cur = flow.flat_index_unchecked(row, col);
        double max_s = s_top;
        for (int s = 0; s < max_steps; ++s)
        {
            const int rr = flow.row_from_flat_unchecked(cur);
            const int cc = flow.col_from_flat_unchecked(cur);
            if (trace_out) trace_out->push_back(rc_to_centre(meta, rr, cc));
            const lsdtt_nearest_polyline_point np =
                nearest_point_on_polyline(rc_to_centre(meta, rr, cc), path);
            if (maxdiv_out && np.distance > *maxdiv_out) *maxdiv_out = np.distance;
            if (np.distance > tol)
            {
                escape_z = static_cast<double>(zdata(rr, cc));
                diverged = true;
                break;
            }
            const double s_here = cum[np.segment_index] +
                np.t * (cum[np.segment_index + 1] - cum[np.segment_index]);
            if (s_here > max_s) max_s = s_here;
            if (max_s >= total - 0.5 * cs) break;
            if (flow.is_outlet_unchecked(cur)) break;
            const std::int64_t nxt = flow.receiver_of_unchecked(cur);
            if (nxt == cur || nxt < 0) break;
            cur = nxt;
        }
        return max_s;
    };

    // Lowest-cell elevation along the path at arc-length s (sampled cell).
    auto path_elev = [&](const std::vector<lsdtt_xy>& path,
                         const std::vector<double>& cum, double s,
                         double& z) -> bool
    {
        int r, c;
        const lsdtt_xy p = point_at_arclength(path, cum, s);
        if (!xy_to_rc(meta, p.x, p.y, r, c)) return false;
        if (mdata(r, c) == 0) return false;
        z = static_cast<double>(zdata(r, c));
        return true;
    };

    for (std::size_t spi = 0; spi < source_paths.size(); ++spi)
    {
        const std::vector<lsdtt_xy>& path = source_paths[spi];
        const bool viz =
            std::find(options.algoviz_source_indices.begin(),
                      options.algoviz_source_indices.end(),
                      static_cast<int>(spi)) !=
            options.algoviz_source_indices.end();

        if (path.size() < 2) continue;
        ++result.n_sources;

        const std::vector<double> cum = cumulative_lengths(path);
        const double total = cum.back();
        if (total < cs) continue;

        // Fallback: force the WHOLE path into the DEM (last resort for sources no
        // sequence of reach-burns can resolve), but only if it actually diverges.
        if (options.force_full_path)
        {
            double esc; bool dv;
            if (trace_reach_s(path, cum, total, 0.0, esc, dv) >= total - 0.5 * cs)
                continue;                                   // already follows
            const std::vector<lsdtt_xy> sub = subline(path, cum, 0.0, total);
            if (burn_subpath(sub))
            {
                ++result.n_burned;
                result.total_burn_length_m += total;
                vectormap_burn_gap g;
                g.path_xy = sub; g.length_m = total; g.burned = true;
                result.burned.push_back(std::move(g));
            }
            continue;
        }

        // Walk the whole path in segments and burn EVERY divergent crossing, not
        // just the first one — so a stuck crossing upstream does not prevent the
        // (often easier) crossings downstream from being burned in the same pass.
        const int m = std::max(1, static_cast<int>(std::ceil(total / seg_len)));
        int i = 0;
        int viz_iter = 0;
        while (i < m)
        {
            const double s_top = i * seg_len;
            double escape_z = 0.0; bool diverged = false;
            std::vector<lsdtt_xy> viz_trace; double viz_maxdiv = 0.0;
            const double reach =
                trace_reach_s(path, cum, total, s_top, escape_z, diverged,
                              viz ? &viz_trace : nullptr,
                              viz ? &viz_maxdiv : nullptr);
            const int reach_i = static_cast<int>(reach / seg_len + 1e-6);

            // Pre-build the algoviz step (geometry/divergence common to all
            // outcomes); the outcome and burn geometry are filled in below.
            vectormap_algoviz_step vstep;
            if (viz)
            {
                vstep.source_index = static_cast<int>(spi);
                vstep.pass         = options.algoviz_pass;
                vstep.iteration    = viz_iter++;
                vstep.s_top        = s_top;
                vstep.d8_trace     = viz_trace;
                vstep.divergence_m = viz_maxdiv;
                vstep.diverged     = diverged;
                vstep.escape_z     = escape_z;
            }

            if (!diverged)
            {
                // The trace stayed within tolerance to the end of its walk (line
                // end / DEM outlet / base level): this span is already carried by
                // the LiDAR. ONE good trace clears many segments, so skip ahead.
                if (viz) { vstep.outcome = 0; vstep.s_low = reach;
                           result.algoviz_steps.push_back(std::move(vstep)); }
                i = (reach_i > i) ? reach_i : i + 1;
                continue;
            }

            // Diverged: the flow followed the path to arc-length `reach` (within
            // tolerance) and then left it for a LOWER competing cell (escape_z).
            // Burn the path from the divergence location downstream to the first
            // pixel below the escape, forcing flow back onto the vector.
            //
            // NOTE: we must act on the divergence even when the trace advanced
            // past segment i before straying. Skipping ahead whenever reach_i > i
            // (the previous behaviour) silently dropped every divergence that was
            // not immediate — so a channel that followed its first segment or two
            // and then got captured by a road downstream was never burned.
            const double s_div = (reach > s_top) ? reach : s_top;

            // Build the SNAPPED-THALWEG polyline downstream of the divergence: at
            // each sample, snap the vector to the lowest cell within `snap_radius`
            // (the local valley bottom). Where the vector is a few metres off the
            // thalweg the snapped path sits IN the thalweg (so the final flow
            // follows the real channel, not the drawn line); where the vector
            // crosses an obstruction the snapped point stays high (no low nearby).
            // FORCE the flow onto the vector (a culvert under a road berm; the
            // LiDAR's competing route is a drainage ditch it must beat). Carve the
            // RAW vector line (no thalweg snap -- the nearby low IS the ditch here),
            // deep enough to sit below the ditch, down to the local channel just
            // below the obstruction.
            const double force_window = 400.0;  // obstruction crossings are local
            const double s_far = std::min(total, s_div + force_window);
            double z_top;
            if (!path_elev(path, cum, s_div, z_top)) { ++i; continue; }

            // Outlet: the first point the vector drops below the escape (the ditch),
            // then extended downstream to the local minimum (the channel bottom)
            // so the deep trench daylights descending, not into a flat.
            double s_out = -1.0, z_prev = z_top;
            bool below = false;
            for (double s = s_div + step; s <= s_far; s += step)
            {
                double z;
                if (!path_elev(path, cum, s, z)) continue;
                if (!below && z < escape_z) below = true;
                if (below)
                {
                    if (z <= z_prev) { s_out = s; z_prev = z; }
                    else break;          // started rising: local min reached
                }
            }
            if (s_out < 0.0)
            {
                // Vector never drops below the ditch within the window -> leave it.
                if (viz) { vstep.outcome = 2; vstep.s_low = s_div;
                           result.algoviz_steps.push_back(std::move(vstep)); }
                i = (reach_i > i) ? reach_i + 1 : i + 1;
                continue;
            }

            // Depth: enough to beat the ditch (escape) plus head-room, at least the
            // configured force depth.
            const double cfg_depth =
                (options.relative_burn_depth_m > 0.0)
                    ? options.relative_burn_depth_m : 4.0;
            const double force_depth =
                std::min(10.0, std::max(cfg_depth, (z_top - escape_z) + 2.0));
            std::vector<lsdtt_xy> sub = subline(path, cum, s_div, s_out);
            const double s_anchor = s_out;
            if (burn_deep(sub, force_depth))
            {
                ++result.n_burned;
                result.total_burn_length_m += (s_anchor - s_div);
                vectormap_burn_gap g;
                g.path_xy = sub; g.length_m = s_anchor - s_div; g.burned = true;
                result.burned.push_back(std::move(g));
                if (viz) { vstep.outcome = 1; vstep.s_low = s_anchor;
                           vstep.burned_segment = sub; }
            }
            if (viz) result.algoviz_steps.push_back(std::move(vstep));
            // Continue downstream past the carved span.
            const int next_i = static_cast<int>(s_anchor / seg_len + 1e-6);
            i = (next_i > i) ? next_i : i + 1;
        }
    }

    return result;
}

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
// vectormap_line_source_cells
// ---------------------------------------------------------------------------
std::vector<std::int64_t> vectormap_line_source_cells(
    const lsdtt_raster_cpu& dem,
    const std::string&      vector_path,
    double                  source_exclusion_radius_m)
{
    const lsdtt_raster_meta& meta = dem.get_meta();
    const int    n_cols = meta.n_cols;
    const lsdtt_array2d_f&  zdata = dem.get_data_array();
    const lsdtt_array2d_u8& mdata = dem.get_mask_array();

    auto valid_elev = [&](const lsdtt_xy& p, double& z) -> bool
    {
        int r, c;
        if (!xy_to_rc(meta, p.x, p.y, r, c)) return false;
        if (mdata(r, c) == 0) return false;
        z = static_cast<double>(zdata(r, c));
        return std::isfinite(z);
    };

    // Walk a fragment from one end inward (densified at cell size) to the first
    // cell with a valid in-DEM elevation. Used where a channel ENTERS the DEM
    // across the NoData fringe: the fragment's source end is off-DEM, so seed at
    // the first valid (entry) cell instead of dropping the source altogether.
    auto first_valid_inward = [&](const std::vector<lsdtt_xy>& pts,
                                  bool from_front,
                                  lsdtt_xy& out, double& out_z) -> bool
    {
        const double step = meta.cell_size;
        auto test = [&](const lsdtt_xy& p) -> bool
        { double z; if (valid_elev(p, z)) { out = p; out_z = z; return true; }
          return false; };
        if (from_front)
            for (std::size_t k = 1; k < pts.size(); ++k)
            {
                const double L=std::hypot(pts[k].x-pts[k-1].x, pts[k].y-pts[k-1].y);
                for (double s=0; s<=L; s+=step){ const double t=(L>0)?s/L:0.0;
                    if (test({pts[k-1].x+(pts[k].x-pts[k-1].x)*t,
                              pts[k-1].y+(pts[k].y-pts[k-1].y)*t})) return true; }
            }
        else
            for (std::size_t k = pts.size(); k-- > 1; )
            {
                const double L=std::hypot(pts[k].x-pts[k-1].x, pts[k].y-pts[k-1].y);
                for (double s=0; s<=L; s+=step){ const double t=(L>0)?s/L:0.0;
                    if (test({pts[k].x+(pts[k-1].x-pts[k].x)*t,
                              pts[k].y+(pts[k-1].y-pts[k].y)*t})) return true; }
            }
        return false;
    };

    // Clip to the in-DEM linework (same boundary handling as the burn) so heads
    // are only seeded for genuine in-DEM channels, not off-DEM fragments.
    const std::vector<std::vector<lsdtt_xy>> clipped_lines =
        clip_vector_lines_to_dem(vector_path, dem);

    // Candidate head = the upstream (higher) end of each fragment.
    struct candidate { lsdtt_xy up; double z; };
    std::vector<candidate> cands;
    // All fragment vertices with a valid elevation, for the "is something higher
    // draining toward this head?" test below.
    std::vector<std::pair<lsdtt_xy, double>> verts;
    for (const std::vector<lsdtt_xy>& pts : clipped_lines)
    {
        if (pts.size() < 2) continue;
        double z_first = 0.0, z_last = 0.0;
        bool ok_first = false, ok_last = false;
        for (std::size_t a = 0; a < pts.size() && !ok_first; ++a)
            ok_first = valid_elev(pts[a], z_first);
        for (std::size_t a = pts.size(); a-- > 0 && !ok_last; )
            ok_last = valid_elev(pts[a], z_last);
        if (!ok_first && !ok_last) continue;
        lsdtt_xy up; double zu = 0.0;
        if (ok_first && ok_last)
        {
            // both ends in-DEM: head is the higher end (unchanged behaviour)
            up = (z_first >= z_last) ? pts.front() : pts.back();
            if (!valid_elev(up, zu)) continue;
        }
        else
        {
            // one end off-DEM: the channel enters the DEM from that (source)
            // side across the NoData fringe. Seed at the first valid cell when
            // walking inward from the off-DEM end (the channel's DEM entry).
            const bool source_is_front = !ok_first;
            if (!first_valid_inward(pts, source_is_front, up, zu)) continue;
        }
        cands.push_back({up, zu});
        for (const lsdtt_xy& p : pts)
        {
            double zp = 0.0;
            if (valid_elev(p, zp)) verts.push_back({p, zp});
        }
    }

    // A candidate is a GENUINE source only if no channel vertex within
    // source_exclusion_radius_m is meaningfully HIGHER — i.e. nothing drains
    // toward it from above. The OS VectorMap "repaired" network is not actually
    // joined at the repairs, so every fragment downstream of a road has an
    // upstream end sitting mid-network, just below the road, with the upstream
    // channel (higher) a short distance away. Without this test each such end
    // seeds a spurious channel that starts in the middle of a real one. A true
    // headwater tip is the local high point and has nothing higher nearby.
    const double R = (source_exclusion_radius_m > 0.0)
                         ? source_exclusion_radius_m : 75.0;
    const double R2 = R * R;
    const double higher_margin = 0.5;   // m

    std::vector<std::int64_t> heads;
    std::vector<std::uint8_t> seen;
    for (const candidate& cand : cands)
    {
        bool higher_nearby = false;
        for (const auto& v : verts)
        {
            if (v.second <= cand.z + higher_margin) continue;
            const double dx = v.first.x - cand.up.x;
            const double dy = v.first.y - cand.up.y;
            if (dx * dx + dy * dy <= R2) { higher_nearby = true; break; }
        }
        if (higher_nearby) continue;   // mid-network / below-road, not a source

        int r, c;
        if (!xy_to_rc(meta, cand.up.x, cand.up.y, r, c)) continue;
        if (mdata(r, c) == 0) continue;
        const std::int64_t flat = static_cast<std::int64_t>(r) * n_cols + c;
        if (seen.empty())
            seen.assign(static_cast<std::size_t>(meta.n_rows) * n_cols, 0);
        if (seen[static_cast<std::size_t>(flat)]) continue;
        seen[static_cast<std::size_t>(flat)] = 1;
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

// ---------------------------------------------------------------------------
// segmented_trace_test
// ---------------------------------------------------------------------------
std::vector<vectormap_segmented_trace> segmented_trace_test(
    const lsdtt_raster_cpu&                   dem,
    const lsdtt_flow_base&                    flow,
    const std::vector<std::vector<lsdtt_xy>>& lines,
    double                                    spacing_m,
    double                                    max_divergence_m,
    int                                       min_near_pixels,
    double                                    endpoint_tol_m,
    int                                       spur_min_pixels,
    const std::vector<int>*                   line_strahler,
    std::vector<std::vector<lsdtt_xy>>*       accepted_segments,
    std::vector<int>*                         accepted_seg_line,
    std::vector<std::vector<lsdtt_xy>>*       accepted_clean,
    std::vector<int>*                         accepted_clean_line,
    int                                       max_trace_steps)
{
    std::vector<vectormap_segmented_trace> out;

    const lsdtt_raster_meta& meta   = dem.get_meta();
    const int    n_rows = meta.n_rows;
    const int    n_cols = meta.n_cols;
    const double cs     = meta.cell_size;
    const double tol    = max_divergence_m;
    const double endpoint_tol =
        (endpoint_tol_m > 0.0) ? endpoint_tol_m : 0.5 * tol;
    const double spacing = std::max(cs, spacing_m);
    const int max_steps = max_trace_steps > 0
                              ? max_trace_steps
                              : 10 * std::max(n_rows, n_cols);

    const lsdtt_array2d_f&  zdata = dem.get_data_array();
    const lsdtt_array2d_u8& mdata = dem.get_mask_array();

    // Elevation of the nearest valid cell to a vertex (for downhill orientation).
    auto vertex_z = [&](const lsdtt_xy& p, double& z) -> bool
    {
        int r, c;
        if (!xy_to_rc(meta, p.x, p.y, r, c)) return false;
        if (mdata(r, c) == 0) return false;
        z = static_cast<double>(zdata(r, c));
        return true;
    };

    const bool want_accept = (accepted_segments != nullptr) ||
                             (accepted_clean    != nullptr);
    // Each accepted (and trimmed) trace: cell-centre xy, the source line (edge)
    // index, and the VML arc-length interval [s0, s1] it covers.
    std::vector<std::vector<lsdtt_xy>>     acc_xy;
    std::vector<int>                       acc_line;
    std::vector<double>                    acc_s0, acc_s1;

    for (std::size_t li = 0; li < lines.size(); ++li)
    {
        std::vector<lsdtt_xy> pts = lines[li];
        if (pts.size() < 2) continue;

        // Orient downhill: first valid vertex should be higher than the last.
        double z_first = 0.0, z_last = 0.0;
        bool have_f = false, have_l = false;
        for (std::size_t k = 0; k < pts.size() && !have_f; ++k)
            have_f = vertex_z(pts[k], z_first);
        for (std::size_t k = pts.size(); k-- > 0 && !have_l; )
            have_l = vertex_z(pts[k], z_last);
        if (have_f && have_l && z_last > z_first)
            std::reverse(pts.begin(), pts.end());

        const std::vector<double> cum = cumulative_lengths(pts);
        const double total = cum.back();
        if (total < cs) continue;

        int seg = 0;
        for (double s0 = 0.0; s0 <= total + 1e-6; s0 += spacing, ++seg)
        {
            vectormap_segmented_trace tr;
            tr.line_index = static_cast<int>(li);
            tr.seg_index  = seg;
            tr.start_s_m  = s0;
            tr.reach_s_m  = s0;

            const lsdtt_xy ptop = point_at_arclength(pts, cum, s0);
            int row, col;
            if (!xy_to_rc(meta, ptop.x, ptop.y, row, col))
            { tr.stop_reason = 5; out.push_back(std::move(tr)); continue; }

            std::vector<std::int64_t> flats;   // trace cells, parallel to tr.xy
            std::vector<double>       dists;    // each cell's distance to the line

            std::int64_t cur = flow.flat_index_unchecked(row, col);
            tr.stop_reason = 4;   // default: max steps
            for (int step = 0; step < max_steps; ++step)
            {
                const int rr = flow.row_from_flat_unchecked(cur);
                const int cc = flow.col_from_flat_unchecked(cur);
                const lsdtt_xy here = rc_to_centre(meta, rr, cc);
                tr.xy.push_back(here);
                flats.push_back(cur);

                const lsdtt_nearest_polyline_point np =
                    nearest_point_on_polyline(here, pts);
                dists.push_back(np.distance);
                if (np.distance > tr.max_divergence_m)
                    tr.max_divergence_m = np.distance;
                if (np.distance > tol) { tr.stop_reason = 0; break; }  // strayed

                ++tr.near_count;                            // this cell is near
                const double s_here = cum[np.segment_index] +
                    np.t * (cum[np.segment_index + 1] - cum[np.segment_index]);
                if (s_here > tr.reach_s_m) tr.reach_s_m = s_here;

                // Reached (near) the downstream end of the line, still on it.
                if (s_here >= total - 0.5 * cs) { tr.stop_reason = 1; break; }

                if (flow.is_outlet_unchecked(cur)) { tr.stop_reason = 2; break; }
                const std::int64_t nxt = flow.receiver_of_unchecked(cur);
                if (nxt == cur || nxt < 0)         { tr.stop_reason = 3; break; }
                cur = nxt;
            }

            // Accept traces that reach the line end or stay near it > N cells.
            if (want_accept &&
                (tr.stop_reason == 1 || tr.near_count > min_near_pixels))
            {
                // Trim back from the downstream end until the last kept cell is
                // within endpoint_tol of the line (so the endpoint sits close to
                // the channel even though the trace may wander mid-way).
                int e = static_cast<int>(flats.size()) - 1;
                while (e >= 0 && dists[e] > endpoint_tol) --e;
                if (e >= 1)
                {
                    acc_xy.emplace_back(tr.xy.begin(), tr.xy.begin() + e + 1);
                    acc_line.push_back(static_cast<int>(li));
                    acc_s0.push_back(tr.start_s_m);
                    acc_s1.push_back(tr.reach_s_m);
                }
            }

            if (tr.xy.size() >= 2) out.push_back(std::move(tr));
        }
    }

    // Initial accepted segments: every accepted, trimmed trace (these overlap).
    if (accepted_segments)
        for (std::size_t i = 0; i < acc_xy.size(); ++i)
            if (acc_xy[i].size() >= 2)
            {
                accepted_segments->push_back(acc_xy[i]);
                if (accepted_seg_line) accepted_seg_line->push_back(acc_line[i]);
            }

    // Cleaned accepted segments: drop contained reaches, then spurs.
    if (accepted_clean && !acc_xy.empty())
    {
        const std::size_t n = acc_xy.size();
        const double eps = 0.5 * cs;

        // Longest interval first.
        std::vector<std::size_t> order(n);
        for (std::size_t i = 0; i < n; ++i) order[i] = i;
        std::sort(order.begin(), order.end(),
            [&](std::size_t a, std::size_t b)
            { return (acc_s1[a] - acc_s0[a]) > (acc_s1[b] - acc_s0[b]); });

        // (1) CONTAINMENT (by arc-length interval, within the same edge): keep a
        // reach only if no already-kept reach on the SAME line starts upstream of
        // AND ends downstream of it. Different edges never contain each other, so
        // tributaries (distinct edges) are always kept separately.
        std::vector<std::size_t> kept;
        for (std::size_t idx : order)
        {
            bool contained = false;
            for (std::size_t k : kept)
            {
                if (acc_line[k] != acc_line[idx]) continue;
                if (acc_s0[k] <= acc_s0[idx] + eps &&
                    acc_s1[k] >= acc_s1[idx] - eps) { contained = true; break; }
            }
            if (!contained) kept.push_back(idx);
        }

        // (2) SPUR removal: drop a kept reach only if its line is a Strahler-1
        // (leaf) edge shorter than spur_min_pixels cells. Trunks (Strahler >= 2)
        // are never removed regardless of length.
        for (std::size_t k : kept)
        {
            if (acc_xy[k].size() < 2) continue;
            const int li = acc_line[k];
            const bool is_leaf = line_strahler &&
                li >= 0 && li < static_cast<int>(line_strahler->size()) &&
                (*line_strahler)[li] <= 1;
            if (is_leaf &&
                static_cast<int>(acc_xy[k].size()) < spur_min_pixels)
                continue;                                   // short leaf -> spur
            accepted_clean->push_back(acc_xy[k]);
            if (accepted_clean_line) accepted_clean_line->push_back(li);
        }
    }

    return out;
}

// ---------------------------------------------------------------------------
// identify_burn_segments
// ---------------------------------------------------------------------------
std::vector<burn_segment_def> identify_burn_segments(
    const lsdtt_raster_cpu&                   dem,
    const lsdtt_vector_to_junction_network&   net,
    const std::vector<std::vector<lsdtt_xy>>& accepted,
    const std::vector<int>&                   accepted_edge_id,
    double tol, double parallel_tol, double min_drop)
{
    std::vector<burn_segment_def> out;
    const lsdtt_raster_meta& meta = dem.get_meta();
    const int n_cols = meta.n_cols;
    // Geometric snap tolerance for "does the downstream edge pass through the
    // receiver junction" (CASE 3). This is topology bookkeeping, NOT the gap
    // tolerance: the downstream edge starts at the receiver junction by
    // construction, so a few cells is plenty. The GAP threshold (`tol`) is kept
    // separate and is 0 by default so that EVERY gap is burned.
    const double geom_tol = 5.0 * meta.cell_size;
    const lsdtt_array2d_f&  zdata = dem.get_data_array();
    const lsdtt_array2d_u8& mdata = dem.get_mask_array();
    const double PI = 3.14159265358979323846;

    auto to_cell = [&](const lsdtt_xy& p) -> std::int64_t
    { int r,c; if(!xy_to_rc(meta,p.x,p.y,r,c)) return -1;
      return static_cast<std::int64_t>(r)*n_cols + c; };
    auto zat = [&](const lsdtt_xy& p, double& z) -> bool
    { int r,c; if(!xy_to_rc(meta,p.x,p.y,r,c)) return false;
      if(mdata(r,c)==0) return false; z=static_cast<double>(zdata(r,c)); return true; };
    auto plen = [&](const std::vector<lsdtt_xy>& a) -> double
    { double L=0; for(std::size_t k=1;k<a.size();++k)
        L+=std::hypot(a[k].x-a[k-1].x,a[k].y-a[k-1].y); return L; };
    auto bearing = [&](const std::vector<lsdtt_xy>& a) -> double
    { return std::atan2(a.back().y-a.front().y, a.back().x-a.front().x); };

    const std::size_t R = accepted.size();

    // ---- 0. drop parallel duplicates (sample the shorter against the longer) -
    std::vector<char>   keep(R, 1);
    std::vector<double> len(R, 0.0);
    for (std::size_t i=0;i<R;++i) len[i]=plen(accepted[i]);
    for (std::size_t i=0;i<R;++i)
    {
        if(!keep[i]||accepted[i].size()<2) continue;
        for (std::size_t j=0;j<R;++j)
        {
            if(i==j||!keep[j]||accepted[j].size()<2) continue;
            if(len[j]>=len[i]) continue;                 // test shorter j vs longer i
            int nin=0,ntot=0;
            for(const lsdtt_xy& q : sample_polyline(accepted[j], 5.0))
            { ++ntot; if(nearest_point_on_polyline(q,accepted[i]).distance<=parallel_tol) ++nin; }
            if(ntot==0) continue;
            double db=std::fabs(bearing(accepted[i])-bearing(accepted[j]));
            if(db>PI) db=std::fabs(db-2*PI);
            if(static_cast<double>(nin)/ntot>0.7 && db<30.0*PI/180.0) keep[j]=0;
        }
    }

    // ---- 1. oriented reaches (up->down) grouped by edge ---------------------
    const auto& edges = net.edges();
    const auto& juncs = net.junctions();
    std::unordered_map<int,int> jidx, eidx;
    for(std::size_t i=0;i<juncs.size();++i) jidx[juncs[i].junction_id]=static_cast<int>(i);
    for(std::size_t i=0;i<edges.size();++i) eidx[edges[i].edge_id]=static_cast<int>(i);
    auto jpos = [&](int jid, lsdtt_xy& p) -> bool
    { auto it=jidx.find(jid); if(it==jidx.end()) return false;
      p.x=juncs[it->second].x; p.y=juncs[it->second].y; return true; };

    struct reach { std::vector<lsdtt_xy> xy; int edge_id=-1; int s_pos=0; };
    std::unordered_map<int,std::vector<reach>> by_edge;
    for(std::size_t i=0;i<R;++i)
    {
        if(!keep[i]||accepted[i].size()<2) continue;
        const int eid = (i<accepted_edge_id.size())? accepted_edge_id[i] : -1;
        auto it=eidx.find(eid); if(it==eidx.end()) continue;
        const auto& E=edges[it->second];
        reach rr; rr.edge_id=eid; rr.xy=accepted[i];
        lsdtt_xy up; if(jpos(E.upstream_junction,up))
        {
            const double d0=std::hypot(rr.xy.front().x-up.x, rr.xy.front().y-up.y);
            const double d1=std::hypot(rr.xy.back().x-up.x,  rr.xy.back().y-up.y);
            if(d1<d0) std::reverse(rr.xy.begin(), rr.xy.end());
        }
        rr.s_pos = nearest_point_on_polyline(rr.xy.front(), E.xy).segment_index;
        by_edge[eid].push_back(std::move(rr));
    }
    for(auto& kv:by_edge)
        std::sort(kv.second.begin(),kv.second.end(),
                  [](const reach&a,const reach&b){return a.s_pos<b.s_pos;});

    // downstream edge: junction -> the edge whose upstream_junction == junction
    std::unordered_map<int,int> out_edge;
    for(const auto& E:edges) if(!E.is_loop_edge) out_edge[E.upstream_junction]=E.edge_id;
    auto down_edge_of=[&](int eid)->int
    { auto it=eidx.find(eid); if(it==eidx.end()) return -1;
      auto o=out_edge.find(edges[it->second].receiver_junction);
      return o==out_edge.end()? -1 : o->second; };
    auto first_reach_down=[&](int eid,int& found)->const reach*
    { int cur=eid,g=0; while(cur>=0&&g++<2000){ auto it=by_edge.find(cur);
        if(it!=by_edge.end()&&!it->second.empty()){found=cur;return &it->second.front();}
        cur=down_edge_of(cur);} return nullptr; };

    // condition endpoints so z(start) > z(end)+min_drop by walking the reaches.
    // up_pl (start's reach, up->down) stepped from its END upward; dn_pl (end's
    // reach) stepped from its START downward. up_pl may be null (source).
    auto condition=[&](const std::vector<lsdtt_xy>* up_pl,
                       const std::vector<lsdtt_xy>* dn_pl,
                       lsdtt_xy s, lsdtt_xy e,
                       std::int64_t& sc, std::int64_t& ec)
    {
        std::size_t si = up_pl? up_pl->size()-1 : 0;  // index of current start vertex
        std::size_t ei = 0;                            // index of current end vertex
        bool up_turn=true;
        for(int it=0; it<200; ++it)
        {
            double zs,ze;
            if(!zat(s,zs)||!zat(e,ze)) break;
            if(zs > ze + min_drop) break;
            bool moved=false;
            if(up_turn && up_pl && si>0){ --si; s=(*up_pl)[si]; moved=true; }
            else if(dn_pl && ei+1<dn_pl->size()){ ++ei; e=(*dn_pl)[ei]; moved=true; }
            else if(up_pl && si>0){ --si; s=(*up_pl)[si]; moved=true; }
            up_turn=!up_turn;
            if(!moved) break;
        }
        sc=to_cell(s); ec=to_cell(e);
    };

    auto emit=[&](int bcase,int upE,int downE,
                  const std::vector<lsdtt_xy>* up_pl,const std::vector<lsdtt_xy>* dn_pl,
                  const lsdtt_xy& s,const lsdtt_xy& e)
    {
        burn_segment_def b; b.burn_case=bcase; b.up_edge=upE; b.down_edge=downE;
        condition(up_pl,dn_pl,s,e,b.start_cell,b.end_cell);
        if(b.start_cell>=0 && b.end_cell>=0 && b.start_cell!=b.end_cell)
            out.push_back(b);
    };

    // ---- CASE 1: source edges whose coverage doesn't start at the source ----
    for(const auto& E:edges)
    {
        if(E.is_loop_edge) continue;
        auto ju=jidx.find(E.upstream_junction); if(ju==jidx.end()) continue;
        const auto& J=juncs[ju->second];
        if(!J.donor_junctions.empty() || J.is_outlet) continue;   // not a source
        int fe=-1; const reach* fr=first_reach_down(E.edge_id, fe);
        if(!fr) continue;                                          // no coverage at all
        lsdtt_xy src{J.x,J.y};
        const lsdtt_xy up_end=fr->xy.front();
        // Snap the start to the first VALID (in-DEM) cell along the source edge:
        // VML sources commonly sit on the clip's NoData fringe (the channel
        // enters from off-DEM), so the raw source cell is NoData and the burn
        // can never descend from it. Walk the edge from the source end inward.
        lsdtt_xy start_pt = src;
        {
            std::vector<lsdtt_xy> eg = E.xy;
            if(eg.size()>=2 &&
               std::hypot(eg.back().x-src.x, eg.back().y-src.y) <
               std::hypot(eg.front().x-src.x, eg.front().y-src.y))
                std::reverse(eg.begin(), eg.end());            // orient source->down
            bool found=false; double zz;
            for(const lsdtt_xy& q : sample_polyline(eg, meta.cell_size))
                if(zat(q,zz)){ start_pt=q; found=true; break; }
            if(!found) continue;                               // edge entirely off-DEM
        }
        if(std::hypot(start_pt.x-up_end.x, start_pt.y-up_end.y) <= tol) continue;
        emit(1, E.edge_id, fe, nullptr, &fr->xy, start_pt, up_end);
    }

    // ---- CASE 2/3: gaps between consecutive reaches along the flow path -----
    for(auto& kv:by_edge)
    {
        const int eid=kv.first; auto& rs=kv.second;
        for(std::size_t k=0;k<rs.size();++k)
        {
            const reach& Rr=rs[k];
            const lsdtt_xy down_end=Rr.xy.back();
            const reach* nextr=nullptr; int nextE=eid;
            if(k+1<rs.size()) nextr=&rs[k+1];
            else { int de=down_edge_of(eid); if(de>=0) nextr=first_reach_down(de,nextE); }
            if(!nextr)
            {
                // channel end. CASE 3 only if receiver is not a DEM outlet.
                auto it=eidx.find(eid); if(it==eidx.end()) continue;
                auto rj=jidx.find(edges[it->second].receiver_junction);
                if(rj!=jidx.end() && !juncs[rj->second].is_outlet)
                {
                    // burn to the next channel down near the receiver junction
                    int de=down_edge_of(eid);
                    auto dit=(de>=0)? eidx.find(de):eidx.end();
                    if(dit!=eidx.end())
                    {
                        const lsdtt_xy& rjp = lsdtt_xy{juncs[rj->second].x,juncs[rj->second].y};
                        lsdtt_nearest_polyline_point np=
                            nearest_point_on_polyline(rjp, edges[dit->second].xy);
                        if(np.distance<=geom_tol)
                            emit(3, eid, de, &Rr.xy, &edges[dit->second].xy,
                                 down_end, np.point);
                    }
                }
                continue;
            }
            const lsdtt_xy up_start=nextr->xy.front();
            if(std::hypot(down_end.x-up_start.x, down_end.y-up_start.y) > tol)
                emit(2, eid, nextE, &Rr.xy, &nextr->xy, down_end, up_start);
        }
    }

    return out;
}

// ---------------------------------------------------------------------------
// burn_lcp_segment
// ---------------------------------------------------------------------------
lcp_burn_result burn_lcp_segment(
    lsdtt_raster_cpu&         dem,
    const std::vector<float>& da_score,
    std::int64_t             start_cell,
    std::int64_t             end_cell,
    double                   min_slope,
    double                   capture_margin,
    double                   uphill_mult,
    double                   channel_alpha,
    int                      window_pad)
{
    lcp_burn_result result;

    const lsdtt_raster_meta& meta = dem.get_meta();
    std::vector<float>&      elev = dem.get_data();
    const lsdtt_array2d_u8&  mask = dem.get_mask_array();
    const int    n_rows = meta.n_rows;
    const int    n_cols = meta.n_cols;
    const double cs      = meta.cell_size;
    const double cs_diag = cs * 1.4142135623730950488;

    auto rc = [&](std::int64_t f, int& r, int& c)
    { r = static_cast<int>(f / n_cols); c = static_cast<int>(f % n_cols); };
    auto steplen = [&](std::int64_t a, std::int64_t b) -> double
    {
        const int dr = std::abs(static_cast<int>(a / n_cols - b / n_cols));
        const int dc = std::abs(static_cast<int>(a % n_cols - b % n_cols));
        return (dr == 1 && dc == 1) ? cs_diag : cs;
    };
    static constexpr int DR[8] = {-1, -1, 0, 1, 1, 1, 0, -1};
    static constexpr int DC[8] = { 0,  1, 1, 1, 0, -1, -1, -1};

    int sr, sc, er, ec; rc(start_cell, sr, sc); rc(end_cell, er, ec);
    if (sr < 0 || sr >= n_rows || sc < 0 || sc >= n_cols ||
        er < 0 || er >= n_rows || ec < 0 || ec >= n_cols) return result;
    if (static_cast<double>(elev[start_cell]) <=
        static_cast<double>(elev[end_cell])) return result;  // caller conditions

    // -- Step 2: least-cost path, length vs drainage area --------------------
    lcp_params lp;
    lp.uphill_mult         = uphill_mult;
    lp.channel_score       = &da_score;
    lp.channel_score_alpha = channel_alpha;
    lp.goal_flats          = { end_cell };
    lp.r_min = std::max(0, std::min(sr, er) - window_pad);
    lp.r_max = std::min(n_rows - 1, std::max(sr, er) + window_pad);
    lp.c_min = std::max(0, std::min(sc, ec) - window_pad);
    lp.c_max = std::min(n_cols - 1, std::max(sc, ec) + window_pad);
    lcp_result lcp = raster_least_cost_path(dem, sr, sc, er, ec, lp);
    if (!lcp.found || lcp.path.size() < 2) return result;
    result.path = lcp.path;

    const std::vector<std::int64_t>& P = result.path;
    const std::size_t n = P.size();
    std::unordered_set<std::int64_t> path_set(P.begin(), P.end());

    // Natural elevations and cumulative flow distance along the path.
    std::vector<double> nat(n), d(n, 0.0);
    for (std::size_t k = 0; k < n; ++k) nat[k] = static_cast<double>(elev[P[k]]);
    for (std::size_t k = 1; k < n; ++k) d[k] = d[k-1] + steplen(P[k-1], P[k]);
    const double total = d.back();
    if (total <= 0.0) return result;
    const double z_start = nat.front();
    const double z_end   = nat.back();

    // -- Steps 3+4: straight-line base, concave-deepen capture, raise -------
    std::vector<double> z(n);
    z[0] = z_start;                                   // start kept
    elev[P[0]] = static_cast<float>(z[0]);
    for (std::size_t k = 1; k < n; ++k)
    {
        const std::int64_t pv = P[k-1], cur = P[k];
        const double L = steplen(pv, cur);

        if (k == n - 1) { z[k] = z_end; elev[cur] = static_cast<float>(z[k]); break; }

        // straight-line base, lower-only, strictly descending
        const double line = z_start - (z_start - z_end) * d[k] / total;
        double zk = std::min(nat[k], line);
        zk = std::min(zk, z[k-1] - min_slope * L);

        // best competing OFF-PATH neighbour slope from pv
        int pr, pc; rc(pv, pr, pc);
        double best_comp = 0.0;
        for (int dir = 0; dir < 8; ++dir)
        {
            const int nr = pr + DR[dir], nc = pc + DC[dir];
            if (nr < 0 || nr >= n_rows || nc < 0 || nc >= n_cols) continue;
            if (mask(nr, nc) == 0) continue;
            const std::int64_t nf = static_cast<std::int64_t>(nr) * n_cols + nc;
            if (nf == cur || path_set.count(nf)) continue;          // off-path only
            const double dz = z[k-1] - static_cast<double>(elev[nf]);
            if (dz <= 0.0) continue;
            const double s = dz / ((dir % 2 == 0) ? cs : cs_diag);
            if (s > best_comp) best_comp = s;
        }

        // floor: must still descend to the fixed target over remaining steps
        const double floor_z =
            z_end + static_cast<double>(n - 1 - k) * min_slope * cs_diag;

        const double needed = z[k-1] - (best_comp * L + capture_margin);
        if (needed < zk)                       // must deepen to capture
        {
            if (needed >= floor_z) zk = needed;          // concave deepen
            else zk = floor_z;                           // floored -> raise below
        }
        if (zk > nat[k])        zk = nat[k];             // never raise a cell
        if (zk >= z[k-1])       zk = z[k-1] - min_slope * L;
        if (zk < floor_z)       zk = floor_z;

        z[k] = zk;
        elev[cur] = static_cast<float>(zk);

        // GUARANTEE D8 follows: cur must be the steepest descent from pv. Concave
        // deepening above captures wherever it can without floating the target;
        // any off-path neighbour of pv that still competes (equal/steeper slope)
        // is raised to just below the slope toward cur. This raise is the spec's
        // LAST RESORT -- it only fires for competitors deepening could not beat
        // (it raises nothing where the carve already wins).
        {
            const double slope_cur = (z[k-1] - zk) / L;
            for (int dir = 0; dir < 8; ++dir)
            {
                const int nr = pr + DR[dir], nc = pc + DC[dir];
                if (nr < 0 || nr >= n_rows || nc < 0 || nc >= n_cols) continue;
                if (mask(nr, nc) == 0) continue;
                const std::int64_t nf = static_cast<std::int64_t>(nr)*n_cols + nc;
                if (nf == cur || path_set.count(nf)) continue;
                const double Lnf = (dir % 2 == 0) ? cs : cs_diag;
                const double slope_nf =
                    (z[k-1] - static_cast<double>(elev[nf])) / Lnf;
                if (slope_nf >= slope_cur)
                {
                    const double new_nf =
                        z[k-1] - Lnf * slope_cur + capture_margin;
                    if (new_nf > static_cast<double>(elev[nf]))
                    { elev[nf] = static_cast<float>(new_nf);
                      result.raised.push_back(nf); }
                }
            }
        }
    }

    // -- Step 5: verify D8 from start reaches end staying on the path corridor.
    // D8 may shortcut across a meander bend to a LATER path cell; that is still
    // connected. Only leaving the path set (an off-path diverter) is a failure.
    {
        std::unordered_map<std::int64_t,std::size_t> pidx;
        for (std::size_t k = 0; k < n; ++k) pidx[P[k]] = k;
        std::int64_t cur = P[0];
        std::size_t idx = 0;
        result.connected = false;
        std::size_t brk_idx = 0; int brk_reason = 0;
        for (int guard = 0; guard < static_cast<int>(n) + 8; ++guard)
        {
            if (cur == end_cell) { result.connected = true; break; }
            int r0, c0; rc(cur, r0, c0);
            std::int64_t best = -1; double best_s = 0.0;
            for (int dir = 0; dir < 8; ++dir)
            {
                const int nr = r0 + DR[dir], nc = c0 + DC[dir];
                if (nr < 0 || nr >= n_rows || nc < 0 || nc >= n_cols) continue;
                if (mask(nr, nc) == 0) continue;
                const std::int64_t nf = static_cast<std::int64_t>(nr)*n_cols + nc;
                const double dz = static_cast<double>(elev[cur]) -
                                  static_cast<double>(elev[nf]);
                if (dz <= 0.0) continue;
                const double s = dz / ((dir % 2 == 0) ? cs : cs_diag);
                if (s > best_s) { best_s = s; best = nf; }
            }
            if (best < 0) { brk_idx=idx; brk_reason=1; break; }   // no descent
            auto it = pidx.find(best);
            if (it == pidx.end()) { brk_idx=idx; brk_reason=2; break; }   // off-path
            if (it->second <= idx) { brk_idx=idx; brk_reason=3; break; }  // backward
            idx = it->second; cur = best;                       // forward (shortcut ok)
        }
        // Reaching the terminus's neighbour (idx >= n-2) is connection: the last
        // step daylights into the downstream channel (which we must not dam).
        if (!result.connected && idx + 2 >= n) result.connected = true;
        result.break_reason = brk_reason; result.break_idx = static_cast<int>(brk_idx);
    }

    return result;
}

} // namespace lsdtt
