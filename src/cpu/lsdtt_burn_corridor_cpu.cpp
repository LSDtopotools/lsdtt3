// lsdtt_burn_corridor_cpu.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
//
// -----------------------------------------------------------------------------
// See lsdtt_burn_corridor_cpu.hpp. GDAL-free.
// -----------------------------------------------------------------------------

#include "lsdtt_burn_corridor_cpu.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <queue>
#include <unordered_map>
#include <unordered_set>
#include <utility>

#include "lsdtt_vector_ops_cpu.hpp"

namespace lsdtt
{

namespace
{

constexpr int k_dr[8] = {-1, -1, -1, 0, 0, 1, 1, 1};
constexpr int k_dc[8] = {-1, 0, 1, -1, 1, -1, 0, 1};

double seg_dist(const lsdtt_xy& p, const lsdtt_xy& a, const lsdtt_xy& b)
{
    const double dx = b.x - a.x, dy = b.y - a.y;
    const double len2 = dx * dx + dy * dy;
    double t = len2 > 0.0 ? ((p.x - a.x) * dx + (p.y - a.y) * dy) / len2 : 0.0;
    t = std::min(1.0, std::max(0.0, t));
    return std::hypot(p.x - (a.x + t * dx), p.y - (a.y + t * dy));
}

double line_dist(const lsdtt_xy& p, const std::vector<lsdtt_xy>& line)
{
    double d = std::numeric_limits<double>::infinity();
    for (size_t i = 0; i + 1 < line.size(); ++i) d = std::min(d, seg_dist(p, line[i], line[i + 1]));
    return d;
}

double line_length(const std::vector<lsdtt_xy>& line)
{
    double s = 0.0;
    for (size_t i = 1; i < line.size(); ++i) s += std::hypot(line[i].x - line[i - 1].x, line[i].y - line[i - 1].y);
    return s;
}

struct grid
{
    int nr = 0, nc = 0;
    double x0 = 0.0, ymax = 0.0, cs = 1.0;

    std::int64_t flat(int r, int c) const { return static_cast<std::int64_t>(r) * nc + c; }
    int row(std::int64_t f) const { return static_cast<int>(f / nc); }
    int col(std::int64_t f) const { return static_cast<int>(f % nc); }
    bool rc_of(double x, double y, int& r, int& c) const
    {
        const double fc = std::floor((x - x0) / cs), fr = std::floor((ymax - y) / cs);
        if (!(fc >= 0 && fc < nc && fr >= 0 && fr < nr)) return false;
        r = static_cast<int>(fr);
        c = static_cast<int>(fc);
        return true;
    }
    lsdtt_xy centre(std::int64_t f) const
    {
        return {x0 + (col(f) + 0.5) * cs, ymax - (row(f) + 0.5) * cs};
    }
};

struct reach_data
{
    std::vector<lsdtt_xy> line;
    int from = -1, to = -1;
    double length = 0.0;
};

} // namespace

burn_corridor_result burn_vector_network_to_corridor(lsdtt_raster_cpu& dem,
                                                     const std::vector<std::vector<lsdtt_xy>>& lines,
                                                     const burn_corridor_options& opts)
{
    burn_corridor_result res;
    const lsdtt_raster_meta& meta = dem.get_meta();
    lsdtt_array2d_f& z = dem.get_data_array();
    const lsdtt_array2d_u8& mask = dem.get_mask_array();
    grid g;
    g.nr = meta.n_rows;
    g.nc = meta.n_cols;
    g.x0 = meta.x_min;
    g.ymax = meta.y_max();
    g.cs = meta.cell_size;
    const double cs = g.cs;
    const double X = std::max(0.0, opts.corridor_half_width_metres);
    const double R = std::max(X, 1.5 * cs);           // end-node disc radius
    auto valid = [&](int r, int c) { return r >= 0 && r < g.nr && c >= 0 && c < g.nc && mask(r, c) != 0; };
    auto zf = [&](std::int64_t f) -> float& { return z(g.row(f), g.col(f)); };

    // ---- network --------------------------------------------------------------
    lsdtt_point_snapper snapper(opts.snap_tolerance_metres >= 0.0 ? opts.snap_tolerance_metres : 0.5 * cs);
    std::vector<reach_data> rs;
    int n_self = 0;
    for (const auto& l : lines)
    {
        if (l.size() < 2 || line_length(l) <= 0.0) continue;
        reach_data rd;
        rd.line = l;
        rd.from = snapper.id_for(l.front());
        rd.to = snapper.id_for(l.back());
        rd.length = line_length(l);
        if (rd.from == rd.to) { ++n_self; continue; }
        rs.push_back(std::move(rd));
    }
    const int nn = static_cast<int>(snapper.size());
    res.n_nodes = nn;
    if (n_self > 0) res.notes.push_back(std::to_string(n_self) + " closed line(s) skipped");
    const int nreach = static_cast<int>(rs.size());
    std::vector<std::vector<int>> in_reaches(static_cast<size_t>(nn));
    for (int k = 0; k < nreach; ++k) in_reaches[static_cast<size_t>(rs[static_cast<size_t>(k)].to)].push_back(k);

    // Topological order (Kahn on reaches) and upstream length.
    std::vector<int> order;
    {
        std::vector<int> pending(static_cast<size_t>(nreach));
        std::vector<std::vector<int>> out_of(static_cast<size_t>(nn));
        for (int k = 0; k < nreach; ++k)
        {
            pending[static_cast<size_t>(k)] =
                static_cast<int>(in_reaches[static_cast<size_t>(rs[static_cast<size_t>(k)].from)].size());
            out_of[static_cast<size_t>(rs[static_cast<size_t>(k)].from)].push_back(k);
        }
        std::queue<int> q;
        for (int k = 0; k < nreach; ++k) if (pending[static_cast<size_t>(k)] == 0) q.push(k);
        std::vector<bool> done(static_cast<size_t>(nreach), false);
        while (!q.empty())
        {
            const int k = q.front();
            q.pop();
            order.push_back(k);
            done[static_cast<size_t>(k)] = true;
            for (const int d : out_of[static_cast<size_t>(rs[static_cast<size_t>(k)].to)])
                if (--pending[static_cast<size_t>(d)] == 0) q.push(d);
        }
        int n_cyc = 0;
        for (int k = 0; k < nreach; ++k) if (!done[static_cast<size_t>(k)]) { order.push_back(k); ++n_cyc; }
        if (n_cyc > 0) res.notes.push_back(std::to_string(n_cyc) + " reach(es) on directed cycles, processed last");
    }
    std::vector<double> up_len(static_cast<size_t>(nreach), 0.0);
    for (const int k : order)
    {
        double s = rs[static_cast<size_t>(k)].length;
        for (const int u : in_reaches[static_cast<size_t>(rs[static_cast<size_t>(k)].from)])
            s += up_len[static_cast<size_t>(u)];
        up_len[static_cast<size_t>(k)] = s;
    }
    std::vector<int> main_inflow(static_cast<size_t>(nn), -1);
    for (int n = 0; n < nn; ++n)
        for (const int k : in_reaches[static_cast<size_t>(n)])
            if (main_inflow[static_cast<size_t>(n)] < 0 ||
                up_len[static_cast<size_t>(k)] > up_len[static_cast<size_t>(main_inflow[static_cast<size_t>(n)])])
                main_inflow[static_cast<size_t>(n)] = k;

    // ---- per-reach processing ----------------------------------------------------
    res.reaches.resize(static_cast<size_t>(nreach));
    std::vector<std::int64_t> arrival(static_cast<size_t>(nn), -1);
    std::vector<std::vector<std::int64_t>> channel_cells(static_cast<size_t>(nreach));
    // cell -> (reach, index in its channel), for channels already processed.
    std::unordered_map<std::int64_t, std::pair<int, size_t>> channel_owner;
    std::vector<bool> processed(static_cast<size_t>(nreach), false);

    auto process = [&](int k)
    {
        const reach_data& rd = rs[static_cast<size_t>(k)];
        burn_corridor_reach& out = res.reaches[static_cast<size_t>(k)];
        out.reach_id = k;
        out.from_node = rd.from;
        out.to_node = rd.to;
        out.length_m = rd.length;
        out.upstream_length_m = up_len[static_cast<size_t>(k)];
        const lsdtt_xy pfrom = snapper.point(rd.from), pto = snapper.point(rd.to);

        // 1. corridor: cell -> distance to the line (m)
        std::unordered_map<std::int64_t, double> corr;
        auto mark = [&](int r, int c, double d)
        {
            if (!valid(r, c)) return;
            const std::int64_t f = g.flat(r, c);
            auto it = corr.find(f);
            if (it == corr.end() || d < it->second) corr[f] = d;
        };
        const int pad = static_cast<int>(std::ceil(X / cs)) + 1;
        for (size_t i = 0; i + 1 < rd.line.size(); ++i)
        {
            const lsdtt_xy a = rd.line[i], b = rd.line[i + 1];
            if (X > 0.0)
            {
                int r0, c0, r1, c1;
                const double xa = std::min(a.x, b.x), xb = std::max(a.x, b.x);
                const double ya = std::min(a.y, b.y), yb = std::max(a.y, b.y);
                r0 = static_cast<int>(std::floor((g.ymax - yb) / cs)) - pad;
                r1 = static_cast<int>(std::floor((g.ymax - ya) / cs)) + pad;
                c0 = static_cast<int>(std::floor((xa - g.x0) / cs)) - pad;
                c1 = static_cast<int>(std::floor((xb - g.x0) / cs)) + pad;
                for (int r = std::max(0, r0); r <= std::min(g.nr - 1, r1); ++r)
                    for (int c = std::max(0, c0); c <= std::min(g.nc - 1, c1); ++c)
                    {
                        const double d = seg_dist(g.centre(g.flat(r, c)), a, b);
                        if (d <= X) mark(r, c, d);
                    }
            }
            // the line's own cells
            const double len = std::hypot(b.x - a.x, b.y - a.y);
            const int n = std::max(1, static_cast<int>(std::ceil(len / (0.25 * cs))));
            for (int s = 0; s <= n; ++s)
            {
                const double t = static_cast<double>(s) / n;
                int r, c;
                if (g.rc_of(a.x + t * (b.x - a.x), a.y + t * (b.y - a.y), r, c) && valid(r, c))
                    mark(r, c, seg_dist(g.centre(g.flat(r, c)), a, b));
            }
        }
        for (const lsdtt_xy& p : {pfrom, pto})
        {
            int rc, cc;
            const int rad = static_cast<int>(std::ceil(R / cs)) + 1;
            const int r_ = static_cast<int>(std::floor((g.ymax - p.y) / cs));
            const int c_ = static_cast<int>(std::floor((p.x - g.x0) / cs));
            for (rc = r_ - rad; rc <= r_ + rad; ++rc)
                for (cc = c_ - rad; cc <= c_ + rad; ++cc)
                {
                    if (!valid(rc, cc)) continue;
                    const lsdtt_xy q = g.centre(g.flat(rc, cc));
                    if (std::hypot(q.x - p.x, q.y - p.y) <= R) mark(rc, cc, line_dist(q, rd.line));
                }
        }
        auto in_corr = [&](std::int64_t f) { return corr.count(f) > 0; };
        auto near_node = [&](std::int64_t f, const lsdtt_xy& p)
        {
            const lsdtt_xy q = g.centre(f);
            return std::hypot(q.x - p.x, q.y - p.y) <= R;
        };

        // 2. start
        std::int64_t start = arrival[static_cast<size_t>(rd.from)];
        if (start < 0)
        {
            for (size_t i = 0; i + 1 < rd.line.size() && start < 0; ++i)
            {
                const lsdtt_xy a = rd.line[i], b = rd.line[i + 1];
                const int n = std::max(1, static_cast<int>(std::ceil(std::hypot(b.x - a.x, b.y - a.y) / (0.5 * cs))));
                for (int s = 0; s <= n && start < 0; ++s)
                {
                    const double t = static_cast<double>(s) / n;
                    int r, c;
                    if (g.rc_of(a.x + t * (b.x - a.x), a.y + t * (b.y - a.y), r, c) && valid(r, c)) start = g.flat(r, c);
                }
            }
        }
        if (start < 0 || corr.empty())
        {
            out.status = "skipped";
            out.skip_reason = "no_valid_cell";
            ++res.n_skipped;
            return;
        }
        corr.emplace(start, line_dist(g.centre(start), rd.line));

        // 3. target
        std::unordered_set<std::int64_t> target;
        const std::int64_t to_arrival = arrival[static_cast<size_t>(rd.to)];
        if (to_arrival < 0)
        {
            for (const auto& kv : corr) if (near_node(kv.first, pto)) target.insert(kv.first);
        }
        else
        {
            // A later inflow joins a channel already laid down at its node: the
            // channel cells of reaches that start or end at the node, within
            // the node disc plus two cells.
            target.insert(to_arrival);
            corr.emplace(to_arrival, line_dist(g.centre(to_arrival), rd.line));
            const double join_r = R + 2.0 * cs;
            for (const auto& kv : channel_owner)
            {
                const reach_data& other = rs[static_cast<size_t>(kv.second.first)];
                if (other.from != rd.to && other.to != rd.to) continue;
                const lsdtt_xy q = g.centre(kv.first);
                if (std::hypot(q.x - pto.x, q.y - pto.y) > join_r) continue;
                target.insert(kv.first);
                corr.emplace(kv.first, line_dist(q, rd.line));
            }
        }
        // Cells of channels already laid down (other than the start and the
        // targets) cannot be crossed: a reach never cuts through another's
        // channel.
        auto blocked = [&](std::int64_t f) { return f != start && !target.count(f) && channel_owner.count(f); };

        if (target.empty())
        {
            // The end node is over NoData (a reach leaving the DEM across its
            // NoData margin): aim at the last valid corridor cell along the
            // line instead, walking back from the end.
            for (size_t i = rd.line.size() - 1; i > 0 && target.empty(); --i)
            {
                const lsdtt_xy a = rd.line[i], b = rd.line[i - 1];
                const int n = std::max(1, static_cast<int>(std::ceil(std::hypot(b.x - a.x, b.y - a.y) / (0.5 * cs))));
                for (int s2 = 0; s2 <= n && target.empty(); ++s2)
                {
                    const double t = static_cast<double>(s2) / n;
                    int r, c;
                    if (g.rc_of(a.x + t * (b.x - a.x), a.y + t * (b.y - a.y), r, c) && valid(r, c) &&
                        in_corr(g.flat(r, c)) && g.flat(r, c) != start)
                        target.insert(g.flat(r, c));
                }
            }
        }
        if (target.empty())
        {
            out.status = "skipped";
            out.skip_reason = "no_target_cell";
            ++res.n_skipped;
            return;
        }

        // 4. steepest descent on the current DEM
        std::vector<std::int64_t> path{start};
        bool followed = target.count(start) > 0;
        {
            const int cap = static_cast<int>(4.0 * rd.length / cs) + 50;
            std::int64_t cur = start;
            for (int step = 0; step < cap && !followed; ++step)
            {
                const int r = g.row(cur), c = g.col(cur);
                double best = 0.0;
                std::int64_t nxt = -1;
                for (int d = 0; d < 8; ++d)
                {
                    const int rr = r + k_dr[d], cc = c + k_dc[d];
                    if (!valid(rr, cc)) continue;
                    const double dist = (k_dr[d] != 0 && k_dc[d] != 0) ? std::sqrt(2.0) : 1.0;
                    const double s = (static_cast<double>(z(r, c)) - z(rr, cc)) / dist;
                    if (s > best) { best = s; nxt = g.flat(rr, cc); }
                }
                if (nxt < 0 || !in_corr(nxt) || blocked(nxt)) break;
                path.push_back(nxt);
                cur = nxt;
                if (target.count(cur)) followed = true;
            }
        }

        if (followed)
        {
            out.status = "followed";
            ++res.n_followed;
        }
        else
        {
            // 5. least-cost path inside the corridor
            const int rf = static_cast<int>(std::ceil((X + cs) / cs));
            std::unordered_map<std::int64_t, double> floor_z;
            for (const auto& kv : corr)
            {
                const int r = g.row(kv.first), c = g.col(kv.first);
                double m = z(r, c);
                for (int dr = -rf; dr <= rf; ++dr)
                    for (int dc = -rf; dc <= rf; ++dc)
                    {
                        const std::int64_t f = g.flat(std::min(g.nr - 1, std::max(0, r + dr)),
                                                      std::min(g.nc - 1, std::max(0, c + dc)));
                        if (corr.count(f)) m = std::min(m, static_cast<double>(zf(f)));
                    }
                floor_z[kv.first] = m;
            }
            std::unordered_map<std::int64_t, double> best;
            std::unordered_map<std::int64_t, std::int64_t> prev;
            using item = std::pair<double, std::int64_t>;
            std::priority_queue<item, std::vector<item>, std::greater<item>> pq;
            best[start] = 0.0;
            pq.push({0.0, start});
            std::int64_t reached = -1;
            while (!pq.empty())
            {
                const auto [d, u] = pq.top();
                pq.pop();
                if (d > best[u]) continue;
                if (target.count(u)) { reached = u; break; }
                const int r = g.row(u), c = g.col(u);
                for (int k8 = 0; k8 < 8; ++k8)
                {
                    const int rr = r + k_dr[k8], cc = c + k_dc[k8];
                    if (!valid(rr, cc)) continue;
                    const std::int64_t v = g.flat(rr, cc);
                    if (!in_corr(v) || blocked(v)) continue;
                    const double step = (k_dr[k8] != 0 && k_dc[k8] != 0) ? std::sqrt(2.0) : 1.0;
                    const double h = std::max(0.0, static_cast<double>(zf(v)) - floor_z[v]);
                    const double nd = d + step * (1.0 + h);
                    auto it = best.find(v);
                    if (it == best.end() || nd < it->second)
                    {
                        best[v] = nd;
                        prev[v] = u;
                        pq.push({nd, v});
                    }
                }
            }
            if (reached < 0)
            {
                out.status = "skipped";
                out.skip_reason = "no_corridor_path";
                ++res.n_skipped;
                return;
            }
            path.clear();
            for (std::int64_t v = reached;; v = prev[v])
            {
                path.push_back(v);
                if (v == start) break;
            }
            std::reverse(path.begin(), path.end());

            // Capturing burn: every path cell must drain to the next path cell
            // more steeply than to any neighbour off the path, or steepest
            // descent on the burned DEM would leave the path again (the next
            // cell is lowered to below z_i - step x (steepest off-path drop +
            // min_gradient)). Target cells (a channel being joined) are valid
            // receivers too, and are not lowered for capture.
            std::unordered_set<std::int64_t> on_path(path.begin(), path.end());
            std::int64_t regrade_from = -1;
            for (size_t i = 0; i + 1 < path.size(); ++i)
            {
                const std::int64_t fi = path[i], fj = path[i + 1];
                const int ri = g.row(fi), ci = g.col(fi);
                const double zi = zf(fi);
                const bool diag = g.row(fj) != ri && g.col(fj) != ci;
                const double step_m = (diag ? std::sqrt(2.0) : 1.0) * cs;
                float& zc = zf(fj);
                if (channel_owner.count(fj))
                {
                    // Joining a channel: only make sure it is below us; if we
                    // arrive below it, lower it and re-grade downstream.
                    const double allowed = zi - step_m * opts.min_gradient;
                    if (zc > allowed)
                    {
                        out.max_burn_depth_m = std::max(out.max_burn_depth_m, static_cast<double>(zc) - allowed);
                        zc = static_cast<float>(allowed);
                        ++res.n_cells_lowered;
                        regrade_from = fj;
                    }
                    continue;
                }
                double side = 0.0;
                for (int d = 0; d < 8; ++d)
                {
                    const int rr = ri + k_dr[d], cc = ci + k_dc[d];
                    if (!valid(rr, cc)) continue;
                    const std::int64_t fn = g.flat(rr, cc);
                    if (on_path.count(fn) || target.count(fn)) continue;
                    const double dist = ((k_dr[d] != 0 && k_dc[d] != 0) ? std::sqrt(2.0) : 1.0) * cs;
                    side = std::max(side, (zi - zf(fn)) / dist);
                }
                const double allowed = zi - step_m * (side + opts.min_gradient);
                if (zc > allowed)
                {
                    out.max_burn_depth_m = std::max(out.max_burn_depth_m, static_cast<double>(zc) - allowed);
                    out.burn_length_m += step_m;
                    zc = static_cast<float>(allowed);
                    ++res.n_cells_lowered;
                }
            }
            // Re-grade the joined channel downstream of the lowered cell, through
            // the reaches below it, lowering only.
            if (regrade_from >= 0)
            {
                auto it = channel_owner.find(regrade_from);
                int m = it->second.first;
                size_t j = it->second.second;
                for (int guard = 0; m >= 0 && guard < nreach; ++guard)
                {
                    const auto& ch = channel_cells[static_cast<size_t>(m)];
                    for (size_t q = j + 1; q < ch.size(); ++q)
                    {
                        const bool dg = g.row(ch[q]) != g.row(ch[q - 1]) && g.col(ch[q]) != g.col(ch[q - 1]);
                        const double allowed = zf(ch[q - 1]) - opts.min_gradient * (dg ? std::sqrt(2.0) : 1.0) * cs;
                        float& zq = zf(ch[q]);
                        if (zq > allowed) { zq = static_cast<float>(allowed); ++res.n_cells_lowered; }
                    }
                    int next = -1;
                    for (int k2 = 0; k2 < nreach; ++k2)
                        if (processed[static_cast<size_t>(k2)] && rs[static_cast<size_t>(k2)].from == rs[static_cast<size_t>(m)].to &&
                            !channel_cells[static_cast<size_t>(k2)].empty() &&
                            channel_cells[static_cast<size_t>(k2)].front() == ch.back())
                            next = k2;
                    m = next;
                    j = 0;
                }
            }
            out.status = "burned";
            ++res.n_burned;
        }

        if (arrival[static_cast<size_t>(rd.to)] < 0) arrival[static_cast<size_t>(rd.to)] = path.back();
        channel_cells[static_cast<size_t>(k)] = path;
        for (size_t j = 0; j < path.size(); ++j) channel_owner.emplace(path[j], std::make_pair(k, j));
        for (const std::int64_t f : path)
        {
            out.channel.push_back(g.centre(f));
            auto it = corr.find(f);
            out.max_separation_m = std::max(out.max_separation_m,
                                            it != corr.end() ? it->second : line_dist(g.centre(f), rd.line));
        }
    };

    // Longest path first: a reach is ready once the main inflow of its start
    // node is done (and, for a later inflow, once its end node's main inflow is
    // done); among ready reaches the one with the longest upstream length goes
    // next, so main stems are laid down before the tributaries that join them.
    for (int done = 0; done < nreach; ++done)
    {
        int pick = -1;
        for (int k = 0; k < nreach; ++k)
        {
            if (processed[static_cast<size_t>(k)]) continue;
            const int mf = main_inflow[static_cast<size_t>(rs[static_cast<size_t>(k)].from)];
            const int mt = main_inflow[static_cast<size_t>(rs[static_cast<size_t>(k)].to)];
            bool ready = (mf < 0 || processed[static_cast<size_t>(mf)]);
            if (mt >= 0 && mt != k) ready &= processed[static_cast<size_t>(mt)];
            if (!ready) continue;
            if (pick < 0 || up_len[static_cast<size_t>(k)] > up_len[static_cast<size_t>(pick)]) pick = k;
        }
        if (pick < 0)
        {
            // only on directed cycles: take the longest remaining
            for (int k = 0; k < nreach; ++k)
                if (!processed[static_cast<size_t>(k)] &&
                    (pick < 0 || up_len[static_cast<size_t>(k)] > up_len[static_cast<size_t>(pick)]))
                    pick = k;
        }
        process(pick);
        processed[static_cast<size_t>(pick)] = true;
    }
    // Final pass: re-grade every channel once more in network order, upstream
    // first (a reach starts on its upstream channel's last cell), so no channel
    // rises from where its upstream channel arrived. A later join or cut can
    // otherwise leave a step on a channel processed earlier, a pit that a
    // re-fill would raise back to the filled surface, tilted towards that
    // depression's own spill point.
    // Channels share cells where one joins another, so re-grading one channel
    // can lower a cell of another already re-graded: repeat until nothing
    // changes (only lowering, so it converges).
    auto final_regrade_once = [&]() -> int
    {
    int changed = 0;
    for (const int k : order)
    {
        const auto& ch = channel_cells[static_cast<size_t>(k)];
        for (size_t q = 1; q < ch.size(); ++q)
        {
            const bool dg = g.row(ch[q]) != g.row(ch[q - 1]) && g.col(ch[q]) != g.col(ch[q - 1]);
            const double allowed = zf(ch[q - 1]) - opts.min_gradient * (dg ? std::sqrt(2.0) : 1.0) * cs;
            float& zq = zf(ch[q]);
            if (zq > allowed)
            {
                zq = static_cast<float>(allowed);
                ++res.n_cells_regraded;
                ++changed;
            }
        }
    }
    return changed;
    };
    auto final_regrade = [&]()
    {
        for (int it = 0; it < 1000 && final_regrade_once() > 0; ++it) {}
    };

    // Drain channel ends inside the DEM. Where the vector network stops inside
    // the DEM (no reach leaves the node), the burned channel can end below its
    // surroundings: a pit, which a re-fill would raise back to a flat tilted
    // towards the depression's own spill point, i.e. back up the channel. From
    // the end cell, find the cheapest route (cost: height above the end cell,
    // plus a little length) to a lower cell or the DEM edge, and carve it so it
    // never rises.
    auto drain_ends = [&]()
    {
        // Nodes with no CHANNEL leaving them: no reach, or only skipped ones.
        std::vector<int> n_out(static_cast<size_t>(nn), 0);
        for (int k2 = 0; k2 < nreach; ++k2)
            if (!channel_cells[static_cast<size_t>(k2)].empty()) ++n_out[static_cast<size_t>(rs[static_cast<size_t>(k2)].from)];
        auto on_edge = [&](std::int64_t f)
        {
            const int r = g.row(f), c = g.col(f);
            for (int d = 0; d < 8; ++d)
                if (!valid(r + k_dr[d], c + k_dc[d])) return true;
            return false;
        };
        for (int n = 0; n < nn; ++n)
        {
            if (n_out[static_cast<size_t>(n)] > 0 || arrival[static_cast<size_t>(n)] < 0) continue;
            const std::int64_t end = arrival[static_cast<size_t>(n)];
            // Follow steepest descent from the end: drained if it reaches the
            // DEM edge; otherwise it stops in a pit, which is carved out to the
            // edge along the cheapest route (a lower neighbour alone is not
            // enough: it can lie in a closed hollow).
            std::int64_t pit = end;
            bool drained = false;
            for (int step = 0; step < g.nr * g.nc; ++step)
            {
                if (on_edge(pit)) { drained = true; break; }
                const int r = g.row(pit), c = g.col(pit);
                double best_s = 0.0;
                std::int64_t nxt = -1;
                for (int d = 0; d < 8; ++d)
                {
                    const int rr = r + k_dr[d], cc = c + k_dc[d];
                    if (!valid(rr, cc)) continue;
                    const double dist = (k_dr[d] != 0 && k_dc[d] != 0) ? std::sqrt(2.0) : 1.0;
                    const double sl = (static_cast<double>(z(r, c)) - z(rr, cc)) / dist;
                    if (sl > best_s) { best_s = sl; nxt = g.flat(rr, cc); }
                }
                if (nxt < 0) break;
                pit = nxt;
            }
            if (drained) continue;
            const double z0 = zf(pit);
            std::unordered_map<std::int64_t, double> best;
            std::unordered_map<std::int64_t, std::int64_t> prev;
            using item = std::pair<double, std::int64_t>;
            std::priority_queue<item, std::vector<item>, std::greater<item>> pq;
            best[pit] = 0.0;
            pq.push({0.0, pit});
            std::int64_t exit_cell = -1;
            while (!pq.empty())
            {
                const auto [d, u] = pq.top();
                pq.pop();
                if (d > best[u]) continue;
                // Exit at the DEM edge, or on a channel cell already lower than
                // the pit (joining a lower channel drains it).
                if (u != pit && (on_edge(u) || (channel_owner.count(u) && zf(u) < z0))) { exit_cell = u; break; }
                for (int k8 = 0; k8 < 8; ++k8)
                {
                    const int rr = g.row(u) + k_dr[k8], cc = g.col(u) + k_dc[k8];
                    if (!valid(rr, cc)) continue;
                    const std::int64_t v = g.flat(rr, cc);
                    // Never route up a channel at or above the pit: carving it
                    // would reverse a channel that drains into this pit.
                    if (channel_owner.count(v) && zf(v) >= z0) continue;
                    const double step = (k_dr[k8] != 0 && k_dc[k8] != 0) ? std::sqrt(2.0) : 1.0;
                    const double nd = d + std::max(0.0, static_cast<double>(zf(v)) - z0) + 1.0e-3 * step;
                    auto it = best.find(v);
                    if (it == best.end() || nd < it->second) { best[v] = nd; prev[v] = u; pq.push({nd, v}); }
                }
            }
            if (exit_cell < 0) continue;
            std::vector<std::int64_t> route;
            for (std::int64_t v = exit_cell;; v = prev[v]) { route.push_back(v); if (v == pit) break; }
            std::reverse(route.begin(), route.end());
            for (size_t q = 1; q < route.size(); ++q)
            {
                const bool dg = g.row(route[q]) != g.row(route[q - 1]) && g.col(route[q]) != g.col(route[q - 1]);
                const double allowed = zf(route[q - 1]) - opts.min_gradient * (dg ? std::sqrt(2.0) : 1.0) * cs;
                float& zq = zf(route[q]);
                if (zq > allowed) { zq = static_cast<float>(allowed); ++res.n_cells_daylighted; }
            }
            ++res.n_ends_drained;
        }
    };

    // Re-grade, then drain ends, until a drain pass carves nothing: the
    // re-grade can lower an end cell below the exit just carved for it, and a
    // carve can cross another channel, so neither alone is final.
    final_regrade();
    for (int pass = 0; pass < 50; ++pass)
    {
        const int before = res.n_cells_daylighted;
        drain_ends();
        if (res.n_cells_daylighted == before) break;
        final_regrade();
    }
    return res;
}

} // namespace lsdtt
