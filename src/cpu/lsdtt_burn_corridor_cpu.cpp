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
#include <numeric>
#include <queue>
#include <set>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>
#include <utility>

#include "lsdtt_depression_hierarchy_cpu.hpp"
#include "lsdtt_vector_ops_cpu.hpp"

namespace lsdtt
{

burn_profile_fit parse_burn_profile_fit(const std::string& text, const std::string& key)
{
    if (text == "linear") return burn_profile_fit::linear;
    if (text == "chi") return burn_profile_fit::chi;
    throw std::invalid_argument(key + " must be one of linear, chi; got '" + text + "'");
}

const char* burn_profile_fit_name(burn_profile_fit fit)
{
    return fit == burn_profile_fit::chi ? "chi" : "linear";
}

// -----------------------------------------------------------------------------
// fit_burn_profile
// -----------------------------------------------------------------------------
burn_profile_output fit_burn_profile(const burn_profile_input& in)
{
    burn_profile_output out;
    const size_t n = in.zc.size();
    out.z.assign(n, 0.0);
    out.record_low.assign(n, 0);
    out.rejected_low.assign(n, 0);
    if (n == 0) return out;
    if (in.s.size() != n || in.x.size() != n)
        throw std::invalid_argument("fit_burn_profile: s, x and zc must have the same length");
    const bool up = std::isfinite(in.upstream_pin);
    const bool dn = std::isfinite(in.downstream_pin);
    const double g = std::max(0.0, in.min_gradient);

    // Anchors: the accepted minima. A candidate k (lower than the last
    // accepted minimum a) is judged by the cut it CAUSES: the cut downstream
    // of it with it, minus the cut there without it.
    //   b: the end of the run falling from k (k and the cells it leads into);
    //   with k: a line from k to J, the first cell past b lower than k; the
    //     cut is how far the DEM in (b, J) stands above that line (held at
    //     z(k) to the end if there is no J);
    //   without k (and the run): a line from a to J2, the first cell past b
    //     lower than a; the cut is how far the DEM in (b, J2) stands above it.
    // A hole is rejected (its whole run: each cell of it is judged with the
    // cells it leads into); the cell just upstream of an obstruction is not
    // (the obstruction is cut either way). Lines are straight in x.
    std::vector<size_t> run_end(n);
    run_end[n - 1] = n - 1;
    for (size_t i = n - 1; i-- > 0;) run_end[i] = (in.zc[i + 1] <= in.zc[i]) ? run_end[i + 1] : i;
    auto cut_above = [&](size_t from, double z_from, size_t lo, size_t hi_excl, size_t to, bool has_to)
    {
        double cut = 0.0;
        for (size_t j = lo; j < hi_excl; ++j)
        {
            double L = z_from;
            if (has_to)
            {
                const double dx = in.x[to] - in.x[from];
                L = dx > 0.0 ? z_from + (in.zc[to] - z_from) * (in.x[j] - in.x[from]) / dx : in.zc[to];
            }
            cut = std::max(cut, in.zc[j] - L);
        }
        return cut;
    };
    std::vector<size_t> anchors{0};
    std::vector<double> value{up ? in.upstream_pin : in.zc[0]};
    if (!up) out.record_low[0] = 1;
    double za = value[0];
    size_t ia = 0;
    for (size_t i = 1; i < n; ++i)
    {
        if (!(in.zc[i] < za)) continue;
        const size_t bnd = run_end[i];
        size_t J = bnd + 1;
        while (J < n && !(in.zc[J] < in.zc[i])) ++J;
        size_t J2 = bnd + 1;
        while (J2 < n && !(in.zc[J2] < za)) ++J2;
        const double with_k = cut_above(i, in.zc[i], bnd + 1, std::min(J, n), J, J < n);
        const double without = cut_above(ia, za, bnd + 1, std::min(J2, n), J2, J2 < n);
        if (with_k - without <= in.max_cut)
        {
            anchors.push_back(i);
            value.push_back(in.zc[i]);
            za = in.zc[i];
            ia = i;
            out.record_low[i] = 1;
        }
        else
        {
            out.rejected_low[i] = 1;
        }
    }
    for (size_t q = 0; q + 1 < anchors.size(); ++q)
    {
        const size_t p = anchors[q], r = anchors[q + 1];
        const double dx = in.x[r] - in.x[p];
        for (size_t j = p; j <= r; ++j)
            out.z[j] = dx > 0.0 ? value[q] + (value[q + 1] - value[q]) * (in.x[j] - in.x[p]) / dx : value[q + 1];
    }
    for (size_t j = anchors.back(); j < n; ++j) out.z[j] = value.back();
    for (size_t i = 1; i < n; ++i) out.z[i] = std::min(out.z[i], out.z[i - 1] - g * (in.s[i] - in.s[i - 1]));

    if (dn)
    {
        const double L = in.s[n - 1];
        if (up && in.upstream_pin - g * L < in.downstream_pin)
        {
            out.pin_conflict = true;
            for (size_t i = 0; i + 1 < n; ++i) out.z[i] = in.upstream_pin - g * in.s[i];
        }
        else
        {
            for (size_t i = 0; i < n; ++i) out.z[i] = std::max(out.z[i], in.downstream_pin + g * (L - in.s[i]));
        }
        out.z[n - 1] = in.downstream_pin;
    }
    if (up) out.z[0] = in.upstream_pin;

    for (size_t i = 0; i < n; ++i)
    {
        if ((i == 0 && up) || (i + 1 == n && dn)) continue;   // pinned cells belong to another path
        out.max_cut_m = std::max(out.max_cut_m, in.zc[i] - out.z[i]);
        out.max_raise_m = std::max(out.max_raise_m, out.z[i] - in.zc[i]);
        out.n_record_lows += out.record_low[i];
        out.n_rejected_lows += out.rejected_low[i];
    }
    out.cut_cap_exceeded = out.max_cut_m > in.max_cut + 1e-9;
    return out;
}

// -----------------------------------------------------------------------------
// fill_small_depressions
// -----------------------------------------------------------------------------
small_pit_fill_result fill_small_depressions(const lsdtt_raster_cpu& dem, int max_pixels,
                                             const std::string& boundary_conditions)
{
    small_pit_fill_result res;
    res.cleaned = dem;
    const size_t n = dem.get_data().size();
    res.pit_id.assign(n, -1);
    if (max_pixels <= 0) return res;

    depression_hierarchy_options o;
    o.boundary_conditions = boundary_conditions;
    o.nodata_is_ocean = true;   // NoData is an outlet, as in the fill (plan 12.13 R3)
    std::vector<dep_label_t> label;
    std::vector<uint8_t> flowdir;
    const lsdtt_depression_hierarchy deps = build_depression_hierarchy(dem, label, flowdir, o);

    std::vector<float>& z = res.cleaned.get_data();
    const std::vector<uint8_t>& m = dem.get_mask();
    const uint32_t N = static_cast<uint32_t>(max_pixels);
    auto is_ocean = [&](dep_label_t d) { return d == k_ocean_label || d == k_no_parent || d >= deps.size(); };
    std::set<dep_label_t> used;
    for (size_t i = 0; i < n; ++i)
    {
        if (!m[i]) continue;
        dep_label_t d = label[i];
        if (is_ocean(d)) continue;
        // the smallest depression containing the cell
        const float zi = z[i];
        while (!is_ocean(d))
        {
            if (zi <= deps[d].out_elev) break;
            if (deps[d].ocean_parent) { d = k_ocean_label; break; }
            d = deps[d].parent;
        }
        if (is_ocean(d)) continue;
        // the largest depression containing it with fewer than N cells
        dep_label_t best = k_no_label;
        for (dep_label_t e = d; !is_ocean(e);)
        {
            if (deps[e].cell_count >= N) break;
            best = e;
            if (deps[e].ocean_parent) break;
            e = deps[e].parent;
        }
        if (best == k_no_label) continue;
        const float spill = deps[best].out_elev;
        if (std::isfinite(spill) && zi < spill)
        {
            z[i] = spill;
            res.pit_id[i] = static_cast<std::int64_t>(best);
            ++res.n_cells;
            used.insert(best);
        }
    }
    res.n_pits = static_cast<long>(used.size());
    return res;
}

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
    double step(std::int64_t a, std::int64_t b) const
    {
        return (row(a) != row(b) && col(a) != col(b)) ? std::sqrt(2.0) * cs : cs;
    }
};

struct reach_data
{
    std::vector<lsdtt_xy> line;
    std::string name;
    int line_index = -1;
    int from = -1, to = -1;
    double length = 0.0;
};

struct path_build
{
    std::vector<int> reaches;
    std::vector<std::int64_t> cells;
    std::vector<int> cell_reach;
};

} // namespace

// -----------------------------------------------------------------------------
// burn_vector_network_to_corridor
// -----------------------------------------------------------------------------
burn_corridor_result burn_vector_network_to_corridor(const lsdtt_raster_cpu& original,
                                                     lsdtt_raster_cpu& out_dem,
                                                     const std::vector<std::vector<lsdtt_xy>>& lines,
                                                     const std::vector<std::string>& names,
                                                     const burn_corridor_options& opts)
{
    burn_corridor_result res;
    const lsdtt_raster_meta& meta = original.get_meta();
    grid g;
    g.nr = meta.n_rows;
    g.nc = meta.n_cols;
    g.x0 = meta.x_min;
    g.ymax = meta.y_max();
    g.cs = meta.cell_size;
    const double cs = g.cs;
    const double X = std::max(0.0, opts.corridor_half_width_metres);
    const double R = std::max(X, 1.5 * cs);           // end-node disc radius
    const double gmin = std::max(1.0e-9, opts.min_gradient);
    const std::int64_t ncell = static_cast<std::int64_t>(g.nr) * g.nc;
    const std::vector<float>& zo = original.get_data();
    const std::vector<uint8_t>& mk = original.get_mask();
    auto valid = [&](int r, int c) { return r >= 0 && r < g.nr && c >= 0 && c < g.nc && mk[static_cast<size_t>(g.flat(r, c))] != 0; };
    auto on_edge = [&](std::int64_t f)
    {
        const int r = g.row(f), c = g.col(f);
        for (int d = 0; d < 8; ++d)
            if (!valid(r + k_dr[d], c + k_dc[d])) return true;
        return false;
    };

    // ---- 1. small-pit cleaning ---------------------------------------------------
    const small_pit_fill_result spf =
        fill_small_depressions(original, opts.max_pit_pixels, opts.boundary_conditions);
    const std::vector<float>& zcl = spf.cleaned.get_data();
    res.n_pits_cleaned = spf.n_pits;
    res.n_pit_cells = spf.n_cells;
    auto zc = [&](std::int64_t f) { return static_cast<double>(zcl[static_cast<size_t>(f)]); };

    // ---- network --------------------------------------------------------------
    lsdtt_point_snapper snapper(opts.snap_tolerance_metres >= 0.0 ? opts.snap_tolerance_metres : 0.5 * cs);
    std::vector<reach_data> rs;
    int n_self = 0;
    for (size_t li = 0; li < lines.size(); ++li)
    {
        const auto& l = lines[li];
        if (l.size() < 2 || line_length(l) <= 0.0) continue;
        reach_data rd;
        rd.line = l;
        if (li < names.size()) rd.name = names[li];
        rd.line_index = static_cast<int>(li);
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

    // ---- 2. the cells under each reach ---------------------------------------------
    res.reaches.resize(static_cast<size_t>(nreach));
    std::vector<std::int64_t> arrival(static_cast<size_t>(nn), -1);
    std::vector<std::vector<std::int64_t>> channel_cells(static_cast<size_t>(nreach));
    std::unordered_map<std::int64_t, int> channel_owner;   // cell -> reach, channels laid down so far
    std::vector<bool> processed(static_cast<size_t>(nreach), false);

    auto adjacent = [&](std::int64_t a, std::int64_t b)
    {
        return a != b && std::abs(g.row(a) - g.row(b)) <= 1 && std::abs(g.col(a) - g.col(b)) <= 1;
    };
    auto process = [&](int k)
    {
        const reach_data& rd = rs[static_cast<size_t>(k)];
        burn_corridor_reach& out = res.reaches[static_cast<size_t>(k)];
        out.reach_id = k;
        out.line_index = rd.line_index;
        out.from_node = rd.from;
        out.to_node = rd.to;
        out.name = rd.name;
        out.length_m = rd.length;
        out.upstream_length_m = up_len[static_cast<size_t>(k)];
        const lsdtt_xy pfrom = snapper.point(rd.from), pto = snapper.point(rd.to);

        // corridor: cell -> distance to the line (m)
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
                const double xa = std::min(a.x, b.x), xb = std::max(a.x, b.x);
                const double ya = std::min(a.y, b.y), yb = std::max(a.y, b.y);
                const int r0 = static_cast<int>(std::floor((g.ymax - yb) / cs)) - pad;
                const int r1 = static_cast<int>(std::floor((g.ymax - ya) / cs)) + pad;
                const int c0 = static_cast<int>(std::floor((xa - g.x0) / cs)) - pad;
                const int c1 = static_cast<int>(std::floor((xb - g.x0) / cs)) + pad;
                for (int r = std::max(0, r0); r <= std::min(g.nr - 1, r1); ++r)
                    for (int c = std::max(0, c0); c <= std::min(g.nc - 1, c1); ++c)
                    {
                        const double d = seg_dist(g.centre(g.flat(r, c)), a, b);
                        if (d <= X) mark(r, c, d);
                    }
            }
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
            const int rad = static_cast<int>(std::ceil(R / cs)) + 1;
            const int r_ = static_cast<int>(std::floor((g.ymax - p.y) / cs));
            const int c_ = static_cast<int>(std::floor((p.x - g.x0) / cs));
            for (int rc = r_ - rad; rc <= r_ + rad; ++rc)
                for (int cc = c_ - rad; cc <= c_ + rad; ++cc)
                {
                    if (!valid(rc, cc)) continue;
                    const lsdtt_xy q = g.centre(g.flat(rc, cc));
                    if (std::hypot(q.x - p.x, q.y - p.y) <= R) mark(rc, cc, line_dist(q, rd.line));
                }
        }
        auto in_corr = [&](std::int64_t f) { return corr.count(f) > 0; };

        // start
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
            return;
        }
        corr.emplace(start, line_dist(g.centre(start), rd.line));

        // target
        std::unordered_set<std::int64_t> target;
        const std::int64_t to_arrival = arrival[static_cast<size_t>(rd.to)];
        if (to_arrival < 0)
        {
            for (const auto& kv : corr)
            {
                const lsdtt_xy q = g.centre(kv.first);
                if (std::hypot(q.x - pto.x, q.y - pto.y) <= R) target.insert(kv.first);
            }
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
                const reach_data& other = rs[static_cast<size_t>(kv.second)];
                if (other.from != rd.to && other.to != rd.to) continue;
                const lsdtt_xy q = g.centre(kv.first);
                if (std::hypot(q.x - pto.x, q.y - pto.y) > join_r) continue;
                target.insert(kv.first);
                corr.emplace(kv.first, line_dist(q, rd.line));
            }
        }
        auto blocked = [&](std::int64_t f) { return f != start && !target.count(f) && channel_owner.count(f); };
        if (target.empty())
        {
            // End node over NoData: aim at the last valid corridor cell along the line.
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
            return;
        }

        // least-cost path inside the corridor, on the cleaned DEM. The local
        // floor is the lowest cell of the corridor BAND (within X of the line;
        // the line's own cells at X = 0) near each cell: the end-node discs
        // are left out, so a deep tributary mouth inside a disc does not make
        // the main valley floor look high.
        const int rf = static_cast<int>(std::ceil((X + cs) / cs));
        const double band = std::max(X, 0.75 * cs);
        std::unordered_map<std::int64_t, double> floor_z;
        for (const auto& kv : corr)
        {
            const int r = g.row(kv.first), c = g.col(kv.first);
            double m = zc(kv.first);
            for (int dr = -rf; dr <= rf; ++dr)
                for (int dc = -rf; dc <= rf; ++dc)
                {
                    const std::int64_t f = g.flat(std::min(g.nr - 1, std::max(0, r + dr)),
                                                  std::min(g.nc - 1, std::max(0, c + dc)));
                    auto it = corr.find(f);
                    if (it != corr.end() && it->second <= band) m = std::min(m, zc(f));
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
                const double h = std::max(0.0, zc(v) - floor_z[v]);
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
            return;
        }
        std::vector<std::int64_t> path;
        for (std::int64_t v = reached;; v = prev[v])
        {
            path.push_back(v);
            if (v == start) break;
        }
        std::reverse(path.begin(), path.end());
        // Remove corners: a cell whose neighbours on the path are themselves
        // neighbours carries no flow (steepest descent takes the diagonal), so
        // it would be a line cell with almost no drainage area.
        {
            std::vector<std::int64_t> pulled;
            pulled.reserve(path.size());
            for (const std::int64_t f : path)
            {
                while (pulled.size() >= 2 && adjacent(pulled[pulled.size() - 2], f)) pulled.pop_back();
                pulled.push_back(f);
            }
            path.swap(pulled);
        }

        out.status = "fitted";
        if (arrival[static_cast<size_t>(rd.to)] < 0) arrival[static_cast<size_t>(rd.to)] = path.back();
        channel_cells[static_cast<size_t>(k)] = path;
        for (const std::int64_t f : path) channel_owner.emplace(f, k);
        for (const std::int64_t f : path)
        {
            out.channel.push_back(g.centre(f));
            auto it = corr.find(f);
            out.max_separation_m = std::max(out.max_separation_m,
                                            it != corr.end() ? it->second : line_dist(g.centre(f), rd.line));
        }
    };

    // Longest path first (as in phase VT): a reach is ready once the main
    // inflow of its start node is done (and, for a later inflow, once its end
    // node's main inflow is done); among ready reaches the longest upstream
    // length goes next.
    std::vector<int> proc_order;
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
            for (int k = 0; k < nreach; ++k)
                if (!processed[static_cast<size_t>(k)] &&
                    (pick < 0 || up_len[static_cast<size_t>(k)] > up_len[static_cast<size_t>(pick)]))
                    pick = k;
        }
        process(pick);
        processed[static_cast<size_t>(pick)] = true;
        proc_order.push_back(pick);
    }

    // ---- 3. paths ------------------------------------------------------------------
    // Cells where a reach channel starts or ends, and how many do.
    std::unordered_map<std::int64_t, int> n_endpoints;
    for (int k = 0; k < nreach; ++k)
    {
        const auto& ch = channel_cells[static_cast<size_t>(k)];
        if (ch.empty()) continue;
        ++n_endpoints[ch.front()];
        ++n_endpoints[ch.back()];
    }
    std::vector<path_build> pb;
    std::vector<int> path_of(static_cast<size_t>(nreach), -1);
    std::vector<size_t> reach_begin(static_cast<size_t>(nreach), 0), reach_end(static_cast<size_t>(nreach), 0);
    std::vector<bool> node_continued(static_cast<size_t>(nn), false);
    for (const int k : proc_order)
    {
        const auto& ch = channel_cells[static_cast<size_t>(k)];
        if (ch.empty()) continue;
        const int from = rs[static_cast<size_t>(k)].from;
        const int u = main_inflow[static_cast<size_t>(from)];
        const bool cont = u >= 0 && path_of[static_cast<size_t>(u)] >= 0 &&
                          pb[static_cast<size_t>(path_of[static_cast<size_t>(u)])].reaches.back() == u &&
                          !node_continued[static_cast<size_t>(from)] &&
                          ch.front() == channel_cells[static_cast<size_t>(u)].back();
        int P;
        size_t first = 0;
        if (cont)
        {
            P = path_of[static_cast<size_t>(u)];
            node_continued[static_cast<size_t>(from)] = true;
            first = 1;
            // A corner at the junction cell (the cells either side of it are
            // neighbours) carries no flow: steepest descent cuts it. Drop it,
            // unless another reach starts or ends there.
            auto& pp = pb[static_cast<size_t>(P)];
            const std::int64_t T = ch.front();
            if (pp.cells.size() >= 2 && ch.size() >= 2 && pp.cells.back() == T &&
                adjacent(pp.cells[pp.cells.size() - 2], ch[1]) && n_endpoints[T] <= 2 &&
                reach_end[static_cast<size_t>(u)] - reach_begin[static_cast<size_t>(u)] >= 2)
            {
                pp.cells.pop_back();
                pp.cell_reach.pop_back();
                --reach_end[static_cast<size_t>(u)];
            }
        }
        else
        {
            P = static_cast<int>(pb.size());
            pb.emplace_back();
        }
        auto& p = pb[static_cast<size_t>(P)];
        p.reaches.push_back(k);
        reach_begin[static_cast<size_t>(k)] = p.cells.size();
        for (size_t i = first; i < ch.size(); ++i)
        {
            p.cells.push_back(ch[i]);
            p.cell_reach.push_back(k);
        }
        reach_end[static_cast<size_t>(k)] = p.cells.size();
        path_of[static_cast<size_t>(k)] = P;
    }
    const int np = static_cast<int>(pb.size());
    std::vector<int> owner_path(static_cast<size_t>(ncell), -1);
    std::vector<int> owner_idx(static_cast<size_t>(ncell), -1);
    // A cell belongs to the first reach that laid it down (processing order),
    // so a later reach ending on an earlier channel finds that channel's path.
    // A path's first cell is owned by the path it starts on, if any.
    for (const int k : proc_order)
    {
        const int P = path_of[static_cast<size_t>(k)];
        if (P < 0) continue;
        const auto& p = pb[static_cast<size_t>(P)];
        const auto& ch = channel_cells[static_cast<size_t>(k)];
        // the reach's own first cell (the start: on another path, or this path's head)
        const size_t f0 = static_cast<size_t>(ch.front());
        if (owner_path[f0] < 0 && reach_begin[static_cast<size_t>(k)] < p.cells.size() &&
            p.cells[reach_begin[static_cast<size_t>(k)]] == ch.front())
        {
            owner_path[f0] = P;
            owner_idx[f0] = static_cast<int>(reach_begin[static_cast<size_t>(k)]);
        }
        for (size_t i = reach_begin[static_cast<size_t>(k)]; i < reach_end[static_cast<size_t>(k)]; ++i)
        {
            const size_t f = static_cast<size_t>(p.cells[i]);
            if (owner_path[f] < 0)
            {
                owner_path[f] = P;
                owner_idx[f] = static_cast<int>(i);
            }
        }
    }

    res.paths.resize(static_cast<size_t>(np));
    std::vector<std::vector<double>> s_of(static_cast<size_t>(np));
    for (int P = 0; P < np; ++P)
    {
        auto& rp = res.paths[static_cast<size_t>(P)];
        const auto& p = pb[static_cast<size_t>(P)];
        rp.path_id = P;
        rp.reach_ids = p.reaches;
        rp.cells = p.cells;
        rp.cell_reach = p.cell_reach;
        auto& s = s_of[static_cast<size_t>(P)];
        s.assign(p.cells.size(), 0.0);
        for (size_t i = 1; i < p.cells.size(); ++i) s[i] = s[i - 1] + g.step(p.cells[i - 1], p.cells[i]);
        rp.length_m = s.empty() ? 0.0 : s.back();
        const size_t f0 = static_cast<size_t>(p.cells.front()), f1 = static_cast<size_t>(p.cells.back());
        if (owner_path[f0] != P) rp.starts_on_path_id = owner_path[f0];
        if (p.cells.size() > 1 && owner_path[f1] != P) rp.joins_path_id = owner_path[f1];
    }

    // Fitting order: a path after the paths it starts on or joins; among the
    // ready ones the longest first.
    std::vector<int> fit_order;
    {
        std::vector<int> pending(static_cast<size_t>(np), 0);
        std::vector<std::vector<int>> dependants(static_cast<size_t>(np));
        for (int P = 0; P < np; ++P)
            for (const int Q : {res.paths[static_cast<size_t>(P)].starts_on_path_id, res.paths[static_cast<size_t>(P)].joins_path_id})
                if (Q >= 0 && Q != P)
                {
                    ++pending[static_cast<size_t>(P)];
                    dependants[static_cast<size_t>(Q)].push_back(P);
                }
        std::vector<bool> done(static_cast<size_t>(np), false);
        using item = std::pair<double, int>;
        std::priority_queue<item> ready;
        for (int P = 0; P < np; ++P)
            if (pending[static_cast<size_t>(P)] == 0) ready.push({res.paths[static_cast<size_t>(P)].length_m, -P});
        for (;;)
        {
            if (ready.empty())
            {
                // dependency cycle: take the longest remaining
                int pick = -1;
                for (int P = 0; P < np; ++P)
                    if (!done[static_cast<size_t>(P)] &&
                        (pick < 0 || res.paths[static_cast<size_t>(P)].length_m > res.paths[static_cast<size_t>(pick)].length_m))
                        pick = P;
                if (pick < 0) break;
                res.notes.push_back("path " + std::to_string(pick) + " is on a dependency cycle; fitted without its unfitted pins");
                ready.push({0.0, -pick});
                pending[static_cast<size_t>(pick)] = 0;
            }
            const int P = -ready.top().second;
            ready.pop();
            if (done[static_cast<size_t>(P)]) continue;
            done[static_cast<size_t>(P)] = true;
            fit_order.push_back(P);
            for (const int D : dependants[static_cast<size_t>(P)])
                if (--pending[static_cast<size_t>(D)] == 0) ready.push({res.paths[static_cast<size_t>(D)].length_m, -D});
        }
    }
    for (size_t q = 0; q < fit_order.size(); ++q) res.paths[static_cast<size_t>(fit_order[q])].path_rank = static_cast<int>(q);

    // ---- 4. fit -----------------------------------------------------------------
    std::vector<bool> fitted(static_cast<size_t>(np), false);
    auto fit_all = [&](const std::vector<double>* area)
    {
        std::fill(fitted.begin(), fitted.end(), false);
        const double theta = opts.m_over_n, A0 = std::max(1e-12, opts.chi_A0);
        for (const int P : fit_order)
        {
            auto& rp = res.paths[static_cast<size_t>(P)];
            const auto& s = s_of[static_cast<size_t>(P)];
            const size_t n = rp.cells.size();
            burn_profile_input in;
            in.s = s;
            in.zc.resize(n);
            for (size_t i = 0; i < n; ++i) in.zc[i] = zc(rp.cells[i]);
            in.x = s;
            rp.drainage_area_m2.assign(n, std::numeric_limits<double>::quiet_NaN());
            rp.chi.assign(n, std::numeric_limits<double>::quiet_NaN());
            if (area)
            {
                double chi = 0.0;
                for (size_t i = 0; i < n; ++i)
                {
                    const double A = std::max(cs * cs, (*area)[static_cast<size_t>(rp.cells[i])]);
                    if (i > 0) chi += std::pow(A0 / A, theta) * (s[i] - s[i - 1]);
                    rp.drainage_area_m2[i] = A;
                    rp.chi[i] = chi;
                    in.x[i] = chi;
                }
            }
            in.min_gradient = gmin;
            in.max_cut = opts.max_cut_metres;
            auto pin_of = [&](int Q, std::int64_t cell) -> double
            {
                if (Q < 0 || !fitted[static_cast<size_t>(Q)]) return std::numeric_limits<double>::quiet_NaN();
                return res.paths[static_cast<size_t>(Q)].z_profile[static_cast<size_t>(owner_idx[static_cast<size_t>(cell)])];
            };
            in.upstream_pin = pin_of(rp.starts_on_path_id, rp.cells.front());
            if (n > 1) in.downstream_pin = pin_of(rp.joins_path_id, rp.cells.back());
            const burn_profile_output o = fit_burn_profile(in);
            rp.z_profile = o.z;
            rp.is_record_low = o.record_low;
            rp.is_rejected_low = o.rejected_low;
            rp.n_record_lows = o.n_record_lows;
            rp.n_rejected_lows = o.n_rejected_lows;
            rp.pin_conflict = o.pin_conflict;
            rp.cut_cap_exceeded = o.cut_cap_exceeded;
            rp.max_cut_m = o.max_cut_m;
            rp.max_raise_m = o.max_raise_m;
            fitted[static_cast<size_t>(P)] = true;
        }
    };

    // ---- 5-7. set the line, drain ends, carve and fill the rest, banks ---------------
    std::vector<float> work;
    std::vector<uint8_t> fixed;
    auto fill_with_fixed = [&]()
    {
        // Priority flood seeded at the edge / NoData and at every fixed cell;
        // fixed cells are never raised. Same rise rule as fill_pits.
        std::vector<uint8_t> seen(static_cast<size_t>(ncell), 0);
        using item = std::pair<float, std::int64_t>;
        std::priority_queue<item, std::vector<item>, std::greater<item>> pq;
        for (std::int64_t f = 0; f < ncell; ++f)
        {
            if (!mk[static_cast<size_t>(f)]) continue;
            if (fixed[static_cast<size_t>(f)] || on_edge(f))
            {
                seen[static_cast<size_t>(f)] = 1;
                pq.push({work[static_cast<size_t>(f)], f});
            }
        }
        long raised = 0;
        while (!pq.empty())
        {
            const auto [zu, u] = pq.top();
            pq.pop();
            const int r = g.row(u), c = g.col(u);
            for (int d = 0; d < 8; ++d)
            {
                const int rr = r + k_dr[d], cc = c + k_dc[d];
                if (!valid(rr, cc)) continue;
                const std::int64_t v = g.flat(rr, cc);
                if (seen[static_cast<size_t>(v)]) continue;
                seen[static_cast<size_t>(v)] = 1;
                float& zv = work[static_cast<size_t>(v)];
                if (zv <= zu)
                {
                    const float dist = static_cast<float>((k_dr[d] != 0 && k_dc[d] != 0) ? std::sqrt(2.0) * cs : cs);
                    const float req = zu + static_cast<float>(gmin) * dist;
                    zv = std::max(req, std::nextafter(zu, std::numeric_limits<float>::infinity()));
                    ++raised;
                }
                pq.push({zv, v});
            }
        }
        return raised;
    };
    auto drain_ends = [&]()
    {
        for (const auto& rp : res.paths)
        {
            if (rp.cells.empty() || rp.joins_path_id >= 0) continue;
            const std::int64_t e = rp.cells.back();
            if (on_edge(e)) continue;
            const double z0 = work[static_cast<size_t>(e)];
            std::unordered_map<std::int64_t, double> best;
            std::unordered_map<std::int64_t, std::int64_t> prev;
            using item = std::pair<double, std::int64_t>;
            std::priority_queue<item, std::vector<item>, std::greater<item>> pq;
            best[e] = 0.0;
            pq.push({0.0, e});
            std::int64_t exit_cell = -1;
            while (!pq.empty())
            {
                const auto [d, u] = pq.top();
                pq.pop();
                if (d > best[u]) continue;
                if (u != e && (on_edge(u) || (fixed[static_cast<size_t>(u)] && work[static_cast<size_t>(u)] < z0)))
                {
                    exit_cell = u;
                    break;
                }
                for (int k8 = 0; k8 < 8; ++k8)
                {
                    const int rr = g.row(u) + k_dr[k8], cc = g.col(u) + k_dc[k8];
                    if (!valid(rr, cc)) continue;
                    const std::int64_t v = g.flat(rr, cc);
                    // never route up a line at or above the end
                    if (fixed[static_cast<size_t>(v)] && work[static_cast<size_t>(v)] >= z0) continue;
                    const double step = (k_dr[k8] != 0 && k_dc[k8] != 0) ? std::sqrt(2.0) : 1.0;
                    const double nd = d + std::max(0.0, static_cast<double>(work[static_cast<size_t>(v)]) - z0) + 1.0e-3 * step;
                    auto it = best.find(v);
                    if (it == best.end() || nd < it->second) { best[v] = nd; prev[v] = u; pq.push({nd, v}); }
                }
            }
            if (exit_cell < 0) continue;
            std::vector<std::int64_t> route;
            for (std::int64_t v = exit_cell;; v = prev[v]) { route.push_back(v); if (v == e) break; }
            std::reverse(route.begin(), route.end());
            const bool to_line = fixed[static_cast<size_t>(exit_cell)] != 0;
            if (to_line)
            {
                // graded exactly from the end to the lower line cell
                const double zx = work[static_cast<size_t>(exit_cell)];
                double L = 0.0;
                for (size_t q = 1; q < route.size(); ++q) L += g.step(route[q - 1], route[q]);
                double s = 0.0;
                for (size_t q = 1; q + 1 < route.size(); ++q)
                {
                    s += g.step(route[q - 1], route[q]);
                    float& zq = work[static_cast<size_t>(route[q])];
                    const float v = static_cast<float>(z0 - (z0 - zx) * s / L);
                    if (zq != v) ++res.n_cells_daylighted;
                    zq = v;
                    fixed[static_cast<size_t>(route[q])] = 1;
                }
            }
            else
            {
                for (size_t q = 1; q < route.size(); ++q)
                {
                    const double allowed = work[static_cast<size_t>(route[q - 1])] - gmin * g.step(route[q - 1], route[q]);
                    float& zq = work[static_cast<size_t>(route[q])];
                    if (zq > allowed) { zq = static_cast<float>(allowed); ++res.n_cells_daylighted; }
                    fixed[static_cast<size_t>(route[q])] = 1;
                }
            }
            ++res.n_ends_drained;
        }
    };
    // where each line cell occurs: (path, index), for the cross-path raise
    std::unordered_map<std::int64_t, std::vector<std::pair<int, int>>> occurrences;
    for (const auto& rp : res.paths)
        for (size_t i = 0; i < rp.cells.size(); ++i)
            occurrences[rp.cells[i]].emplace_back(rp.path_id, static_cast<int>(i));
    auto banks = [&]()
    {
        // Each line cell must drain to its own next cell. A lower off-line
        // receiver is raised into a bank; a line cell of ANOTHER path that
        // would take the flow (a distributary at a bifurcation, a parallel
        // channel) is raised until it is less steep than the own-next step,
        // where its own path stays monotone (otherwise it is counted).
        long n = 0, n_raised_line = 0;
        for (const auto& rp : res.paths)
        {
            for (size_t i = 0; i + 1 < rp.cells.size(); ++i)
            {
                const std::int64_t f = rp.cells[i];
                if (owner_path[static_cast<size_t>(f)] != rp.path_id) continue;
                const std::int64_t own_next = rp.cells[i + 1];
                for (int guard = 0; guard < 16; ++guard)
                {
                    const int r = g.row(f), c = g.col(f);
                    double best_s = 0.0;
                    std::int64_t recv = -1;
                    double recv_dist = cs;
                    for (int d = 0; d < 8; ++d)
                    {
                        const int rr = r + k_dr[d], cc = c + k_dc[d];
                        if (!valid(rr, cc)) continue;
                        const std::int64_t v = g.flat(rr, cc);
                        const double dist = (k_dr[d] != 0 && k_dc[d] != 0) ? std::sqrt(2.0) * cs : cs;
                        const double sl = (static_cast<double>(work[static_cast<size_t>(f)]) - work[static_cast<size_t>(v)]) / dist;
                        // ties go against the own next cell, so it must be strictly steeper
                        if (sl > best_s || (sl == best_s && sl > 0.0 && v != own_next))
                        {
                            best_s = sl;
                            recv = v;
                            recv_dist = dist;
                        }
                    }
                    if (recv < 0 || recv == own_next) break;
                    const float zf = work[static_cast<size_t>(f)];
                    if (!fixed[static_cast<size_t>(recv)])
                    {
                        work[static_cast<size_t>(recv)] = std::max(zf + static_cast<float>(gmin * recv_dist),
                                                                   std::nextafter(zf, std::numeric_limits<float>::infinity()));
                        ++n;
                        if (rp.z_profile[i] > zc(f) + 1e-9) ++n_raised_line;
                        continue;
                    }
                    // a fixed cell of another path (or an exit route)
                    const int Q = owner_path[static_cast<size_t>(recv)];
                    int qi = Q >= 0 ? owner_idx[static_cast<size_t>(recv)] : -1;
                    const double s_own = (static_cast<double>(zf) - work[static_cast<size_t>(own_next)]) / g.step(f, own_next);
                    const double want = static_cast<double>(zf) - 0.5 * s_own * recv_dist;
                    // every path through the cell must stay monotone
                    double cap = std::numeric_limits<double>::infinity();
                    bool head_of_a_path = false;
                    auto oc = occurrences.find(recv);
                    if (oc != occurrences.end())
                        for (const auto& [P2, j2] : oc->second)
                        {
                            if (j2 == 0) { head_of_a_path = true; continue; }
                            const std::int64_t pred = res.paths[static_cast<size_t>(P2)].cells[static_cast<size_t>(j2 - 1)];
                            cap = std::min(cap, static_cast<double>(work[static_cast<size_t>(pred)]) - gmin * g.step(pred, recv));
                        }
                    if (head_of_a_path) qi = 0;
                    if (Q < 0 || qi <= 0 || want > cap)
                    {
                        ++res.n_cross_path_unresolved;
                        break;
                    }
                    work[static_cast<size_t>(recv)] = static_cast<float>(want);
                    ++res.n_cross_path_raised;
                }
            }
        }
        res.n_bank_cells += n;
        res.n_bank_cells_at_raised_line += n_raised_line;
        return n;
    };
    auto condition = [&]()
    {
        work = zo;
        fixed.assign(static_cast<size_t>(ncell), 0);
        res.n_ends_drained = 0;
        res.n_cells_daylighted = 0;
        res.n_bank_cells = 0;
        res.n_bank_cells_at_raised_line = 0;
        res.n_cross_path_raised = 0;
        res.n_cross_path_unresolved = 0;
        for (const auto& rp : res.paths)
            for (size_t i = 0; i < rp.cells.size(); ++i)
            {
                const size_t f = static_cast<size_t>(rp.cells[i]);
                if (owner_path[f] != rp.path_id) continue;
                work[f] = static_cast<float>(rp.z_profile[i]);
                fixed[f] = 1;
            }
        drain_ends();
        // carve: least-cost breach with the fixed cells as outlets
        {
            lsdtt_raster_cpu tmp(std::vector<float>(work), std::vector<uint8_t>(mk), meta);
            breach_least_cost_options bo = opts.breach;
            bo.fill_remainder = false;
            bo.fixed_cells = &fixed;
            res.breach = breach_least_cost_diagnostics{};
            const lsdtt_raster_cpu carved = breach_depressions_least_cost(tmp, bo, &res.breach);
            const std::vector<float>& cz = carved.get_data();
            for (std::int64_t f = 0; f < ncell; ++f)
                if (mk[static_cast<size_t>(f)] && !fixed[static_cast<size_t>(f)]) work[static_cast<size_t>(f)] = cz[static_cast<size_t>(f)];
        }
        fill_with_fixed();
        if (banks() > 0) fill_with_fixed();
    };
    auto d8_area = [&]()
    {
        std::vector<double> area(static_cast<size_t>(ncell), 0.0);
        std::vector<std::int64_t> cells;
        for (std::int64_t f = 0; f < ncell; ++f)
            if (mk[static_cast<size_t>(f)]) { cells.push_back(f); area[static_cast<size_t>(f)] = cs * cs; }
        std::stable_sort(cells.begin(), cells.end(), [&](std::int64_t a, std::int64_t b)
                         { return work[static_cast<size_t>(a)] > work[static_cast<size_t>(b)]; });
        for (const std::int64_t f : cells)
        {
            const int r = g.row(f), c = g.col(f);
            double best_s = 0.0;
            std::int64_t recv = -1;
            for (int d = 0; d < 8; ++d)
            {
                const int rr = r + k_dr[d], cc = c + k_dc[d];
                if (!valid(rr, cc)) continue;
                const std::int64_t v = g.flat(rr, cc);
                const double dist = (k_dr[d] != 0 && k_dc[d] != 0) ? std::sqrt(2.0) : 1.0;
                const double sl = (static_cast<double>(work[static_cast<size_t>(f)]) - work[static_cast<size_t>(v)]) / dist;
                if (sl > best_s) { best_s = sl; recv = v; }
            }
            if (recv >= 0) area[static_cast<size_t>(recv)] += area[static_cast<size_t>(f)];
        }
        return area;
    };

    fit_all(nullptr);
    condition();
    if (opts.fit == burn_profile_fit::chi)
    {
        const std::vector<double> a1 = d8_area();
        fit_all(&a1);
        condition();
        const std::vector<double> a2 = d8_area();
        std::vector<double> ch;
        for (const auto& rp : res.paths)
            for (const std::int64_t f : rp.cells)
            {
                const double A1 = a1[static_cast<size_t>(f)], A2 = a2[static_cast<size_t>(f)];
                if (A1 > 0.0) ch.push_back(std::fabs(A2 - A1) / A1);
            }
        if (!ch.empty())
        {
            std::sort(ch.begin(), ch.end());
            res.max_area_change = ch.back();
            res.median_area_change = ch[ch.size() / 2];
            res.frac_area_change_over_10pct =
                static_cast<double>(ch.end() - std::upper_bound(ch.begin(), ch.end(), 0.1)) / static_cast<double>(ch.size());
        }
    }

    // ---- outputs ---------------------------------------------------------------------
    for (const int P : fit_order)
    {
        auto& rp = res.paths[static_cast<size_t>(P)];
        const auto& s = s_of[static_cast<size_t>(P)];
        const size_t n = rp.cells.size();
        rp.xy.resize(n);
        rp.distance_m = s;
        rp.z_original.resize(n);
        rp.z_cleaned.resize(n);
        rp.pit_id.resize(n);
        rp.flow_distance_m.resize(n);
        double base = 0.0;
        if (rp.joins_path_id >= 0)
        {
            const auto& q = res.paths[static_cast<size_t>(rp.joins_path_id)];
            const int j = owner_idx[static_cast<size_t>(rp.cells.back())];
            if (!q.flow_distance_m.empty()) base = q.flow_distance_m[static_cast<size_t>(j)];
            rp.join_distance_m = s_of[static_cast<size_t>(rp.joins_path_id)][static_cast<size_t>(j)];
        }
        for (size_t i = 0; i < n; ++i)
        {
            const size_t f = static_cast<size_t>(rp.cells[i]);
            rp.xy[i] = g.centre(rp.cells[i]);
            rp.z_original[i] = zo[f];
            rp.z_cleaned[i] = zcl[f];
            rp.pit_id[i] = spf.pit_id[f];
            rp.flow_distance_m[i] = base + (s.back() - s[i]);
        }
        res.n_record_lows += rp.n_record_lows;
        res.n_rejected_lows += rp.n_rejected_lows;
        if (rp.pin_conflict) ++res.n_pin_conflicts;
        if (rp.cut_cap_exceeded) ++res.n_cut_cap_exceeded;
        // per reach
        for (size_t i = 0; i < n; ++i)
        {
            if ((i == 0 && rp.starts_on_path_id >= 0) || (i + 1 == n && rp.joins_path_id >= 0)) continue;
            auto& rc = res.reaches[static_cast<size_t>(rp.cell_reach[i])];
            rc.path_id = P;
            rc.max_cut_m = std::max(rc.max_cut_m, rp.z_cleaned[i] - rp.z_profile[i]);
            rc.max_raise_m = std::max(rc.max_raise_m, rp.z_profile[i] - rp.z_cleaned[i]);
            rc.n_record_lows += rp.is_record_low[i];
            rc.n_rejected_lows += rp.is_rejected_low[i];
        }
        for (const int k : rp.reach_ids) res.reaches[static_cast<size_t>(k)].path_id = P;
    }
    for (const auto& rc : res.reaches)
    {
        if (rc.status == "fitted") ++res.n_fitted;
        else ++res.n_skipped;
    }
    for (std::int64_t f = 0; f < ncell; ++f)
    {
        const size_t i = static_cast<size_t>(f);
        if (!mk[i]) continue;
        const bool on_line = owner_path[i] >= 0;
        if (on_line)
        {
            ++res.n_line_cells;
            if (work[i] < zo[i]) ++res.n_line_cells_cut;
            else if (work[i] > zo[i]) ++res.n_line_cells_raised;
        }
        else
        {
            if (work[i] < zo[i]) ++res.n_offline_cells_carved;
            else if (work[i] > zo[i]) ++res.n_offline_cells_raised;
        }
    }
    out_dem = lsdtt_raster_cpu(std::vector<float>(work), std::vector<uint8_t>(mk), meta);
    return res;
}

} // namespace lsdtt
