// lsdtt_valley_relative_elevation_cpu.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// See lsdtt_valley_relative_elevation_cpu.hpp.
// -----------------------------------------------------------------------------

#include "lsdtt_valley_relative_elevation_cpu.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <stdexcept>

#include "../utils/lsdtt_grid_geodesic.hpp"
#include "../utils/lsdtt_grid_neighbours.hpp"
#include "../utils/lsdtt_point_bucket_index.hpp"
#include "../utils/lsdtt_raster_geometry.hpp"
#include "lsdtt_proximity_cpu.hpp"
#include "lsdtt_valley_centreline_cpu.hpp"

namespace lsdtt
{

namespace
{

struct pt { double x, y; };

pt flat_pt(const lsdtt_raster_meta& meta, std::int64_t f)
{
    const lsdtt_xy p = cell_center_from_row_col(
        meta, static_cast<int>(f / meta.n_cols), static_cast<int>(f % meta.n_cols));
    return {p.x, p.y};
}

double cross2(double ax, double ay, double bx, double by) { return ax * by - ay * bx; }

// Segment p->p2 against q->q2. On intersection, u is the fraction along p->p2
// and v along q->q2 (LSDTopoTools2 segmentsIntersect2D, endpoints included).
bool segment_intersection(const pt& p, const pt& p2, const pt& q, const pt& q2,
                          double& u, double& v)
{
    const double rx = p2.x - p.x, ry = p2.y - p.y;
    const double sx = q2.x - q.x, sy = q2.y - q.y;
    const double rxs = cross2(rx, ry, sx, sy);
    if (std::abs(rxs) < 1e-12) return false;   // parallel or collinear
    const double qpx = q.x - p.x, qpy = q.y - p.y;
    u = cross2(qpx, qpy, sx, sy) / rxs;
    v = cross2(qpx, qpy, rx, ry) / rxs;
    const double e = 1e-9;
    return u >= -e && u <= 1.0 + e && v >= -e && v <= 1.0 + e;
}

std::vector<double> along_of(const std::vector<pt>& p)
{
    std::vector<double> a(p.size(), 0.0);
    for (std::size_t i = 1; i < p.size(); ++i)
        a[i] = a[i - 1] + std::hypot(p[i].x - p[i - 1].x, p[i].y - p[i - 1].y);
    return a;
}

} // namespace

centreline_baseline_choice parse_centreline_baseline_choice(const std::string& name)
{
    if (name == "projected_monotonic") return centreline_baseline_choice::projected_monotonic;
    if (name == "channel_crossings") return centreline_baseline_choice::channel_crossings;
    throw std::invalid_argument("unknown centreline baseline choice '" + name +
                                "'. Valid: projected_monotonic, channel_crossings.");
}

std::vector<float> compute_centreline_baseline(
    const std::vector<std::int64_t>& centreline_flats,
    const std::vector<std::int64_t>& channel_flats,
    const std::vector<float>&        elevation,
    const lsdtt_raster_meta&         meta,
    centreline_baseline_choice       choice)
{
    const std::size_t nc = centreline_flats.size();
    const float nan = std::numeric_limits<float>::quiet_NaN();
    std::vector<float> base(nc, nan);
    if (nc == 0 || channel_flats.empty()) return base;

    std::vector<pt> cl(nc), ch(channel_flats.size());
    for (std::size_t i = 0; i < nc; ++i) cl[i] = flat_pt(meta, centreline_flats[i]);
    for (std::size_t i = 0; i < ch.size(); ++i) ch[i] = flat_pt(meta, channel_flats[i]);
    const std::vector<double> t = along_of(cl);
    auto z_of = [&](std::size_t i) {
        return elevation[static_cast<std::size_t>(channel_flats[i])]; };

    // Bucket index over the centreline nodes: both methods ask, for each
    // channel pixel, which centreline nodes are near it.
    std::vector<double> clx(nc), cly(nc);
    for (std::size_t k = 0; k < nc; ++k) { clx[k] = cl[k].x; cly[k] = cl[k].y; }
    const point_bucket_index cl_index(clx, cly, 4.0 * meta.cell_size);

    if (choice == centreline_baseline_choice::projected_monotonic)
    {
        // Every channel pixel pins its nearest centreline node; lowest wins.
        std::vector<float> pins(nc, nan);
        for (std::size_t i = 0; i < ch.size(); ++i)
        {
            const std::size_t best = cl_index.nearest(ch[i].x, ch[i].y);
            const float z = z_of(i);
            if (std::isfinite(z) && (!std::isfinite(pins[best]) || z < pins[best]))
                pins[best] = z;
        }
        return fit_pinned_monotonic_baseline(t, pins);
    }

    // channel_crossings: pins where the channel crosses the centreline.
    struct crossing { double t; double z; };
    std::vector<crossing> cr;
    // A centreline node that IS a channel pixel is a crossing too. Where the
    // two lines run along the same pixels their segments are collinear, which
    // the intersection test rejects (it would only register where the lines
    // change direction), leaving the baseline interpolated across a reach
    // where the channel is right there. LSDTopoTools2's centreline was a
    // smoothed line that never coincided with the channel, so it never met
    // this case.
    {
        std::vector<std::pair<std::int64_t, std::size_t>> chs;
        chs.reserve(channel_flats.size());
        for (std::size_t i = 0; i < channel_flats.size(); ++i) chs.push_back({channel_flats[i], i});
        std::sort(chs.begin(), chs.end());
        for (std::size_t k = 0; k < nc; ++k)
        {
            const auto it = std::lower_bound(
                chs.begin(), chs.end(), std::make_pair(centreline_flats[k], std::size_t{0}));
            if (it != chs.end() && it->first == centreline_flats[k])
                cr.push_back({t[k], z_of(it->second)});
        }
    }
    // A centreline segment that crosses channel segment i has both ends within
    // (its own length + i's length) of channel node i, so only the segments
    // touching the nodes found within that radius are tested.
    double max_cl_seg = 0.0;
    for (std::size_t j = 0; j + 1 < nc; ++j)
        max_cl_seg = std::max(max_cl_seg, std::hypot(cl[j + 1].x - cl[j].x, cl[j + 1].y - cl[j].y));
    std::vector<std::size_t> near_nodes, cand;
    for (std::size_t i = 0; i + 1 < ch.size(); ++i)
    {
        const double cx0 = std::min(ch[i].x, ch[i + 1].x), cx1 = std::max(ch[i].x, ch[i + 1].x);
        const double cy0 = std::min(ch[i].y, ch[i + 1].y), cy1 = std::max(ch[i].y, ch[i + 1].y);
        const double reach = max_cl_seg +
            std::hypot(ch[i + 1].x - ch[i].x, ch[i + 1].y - ch[i].y) + 1e-6;
        near_nodes.clear();
        cand.clear();
        cl_index.within(ch[i].x, ch[i].y, reach, near_nodes);
        for (const std::size_t q : near_nodes)
        {
            if (q > 0) cand.push_back(q - 1);
            if (q + 1 < nc) cand.push_back(q);
        }
        std::sort(cand.begin(), cand.end());
        cand.erase(std::unique(cand.begin(), cand.end()), cand.end());
        for (const std::size_t j : cand)
        {
            if (std::max(cl[j].x, cl[j + 1].x) < cx0 || std::min(cl[j].x, cl[j + 1].x) > cx1 ||
                std::max(cl[j].y, cl[j + 1].y) < cy0 || std::min(cl[j].y, cl[j + 1].y) > cy1)
                continue;
            double u, v;
            if (!segment_intersection(ch[i], ch[i + 1], cl[j], cl[j + 1], u, v)) continue;
            u = std::clamp(u, 0.0, 1.0);
            v = std::clamp(v, 0.0, 1.0);
            const double z = z_of(i) + u * (static_cast<double>(z_of(i + 1)) - z_of(i));
            cr.push_back({t[j] + v * (t[j + 1] - t[j]), z});
        }
    }
    if (cr.empty()) return base;
    std::sort(cr.begin(), cr.end(), [](const crossing& a, const crossing& b) { return a.t < b.t; });
    // Crossings closer than 1 mm along the centreline: keep the later one
    // (LSDTopoTools2 behaviour; shared vertices are found twice).
    std::vector<crossing> uq{cr.front()};
    for (std::size_t i = 1; i < cr.size(); ++i)
    {
        if (std::abs(cr[i].t - uq.back().t) > 1e-3) uq.push_back(cr[i]);
        else uq.back() = cr[i];
    }
    if (uq.size() == 1)
    {
        std::fill(base.begin(), base.end(), static_cast<float>(uq[0].z));
        return base;
    }
    std::size_t c = 0;
    for (std::size_t k = 0; k < nc; ++k)
    {
        const double tk = t[k];
        if (tk <= uq.front().t) { base[k] = static_cast<float>(uq.front().z); continue; }
        if (tk >= uq.back().t)  { base[k] = static_cast<float>(uq.back().z);  continue; }
        while (c + 1 < uq.size() && uq[c + 1].t < tk) ++c;
        const double t0 = uq[c].t, t1 = uq[c + 1].t;
        const double a = t1 > t0 ? (tk - t0) / (t1 - t0) : 0.0;
        base[k] = static_cast<float>(uq[c].z + a * (uq[c + 1].z - uq[c].z));
    }
    return base;
}

valley_rer_result compute_valley_relative_elevation(
    const lsdtt_raster_cpu&                   dem,
    const lsdtt_raster_cpu&                   floor_mask,
    const lsdtt_raster_cpu&                   domain_mask,
    const std::vector<valley_rer_centreline>& centrelines,
    const valley_rer_request&                 request)
{
    const lsdtt_raster_meta& meta = dem.get_meta();
    const int n_rows = meta.n_rows, n_cols = meta.n_cols;
    const double cell = meta.cell_size;
    const std::size_t n = static_cast<std::size_t>(n_rows) * n_cols;
    if (floor_mask.get_data().size() != n || domain_mask.get_data().size() != n)
        throw std::invalid_argument(
            "compute_valley_relative_elevation: floor / domain mask / DEM size "
            "mismatch.");
    for (const valley_rer_centreline& c : centrelines)
        if (c.baseline.size() != c.flats.size() || c.priority.size() != c.flats.size())
            throw std::invalid_argument(
                "compute_valley_relative_elevation: baseline / priority length "
                "differs from the centreline's.");

    const std::vector<float>&        elev  = dem.get_data();
    const std::vector<std::uint8_t>& dmask = dem.get_mask();
    auto on = [](const lsdtt_raster_cpu& r, std::size_t i) {
        return r.get_mask()[i] && r.get_data()[i] > 0.5f; };

    // Floor and domain (the floor is always in the domain).
    std::vector<std::uint8_t> is_floor(n, 0), in_dom(n, 0);
    for (std::size_t i = 0; i < n; ++i)
    {
        if (!dmask[i]) continue;
        is_floor[i] = on(floor_mask, i) ? 1 : 0;
        in_dom[i] = (is_floor[i] || on(domain_mask, i)) ? 1 : 0;
    }

    // Compact index over floor pixels (the discs only ever touch the floor).
    std::vector<std::int32_t> fid(n, -1);
    std::vector<std::int64_t> fflat;
    for (std::size_t i = 0; i < n; ++i)
        if (is_floor[i])
        {
            fid[i] = static_cast<std::int32_t>(fflat.size());
            fflat.push_back(static_cast<std::int64_t>(i));
        }
    const std::size_t nf = fflat.size();

    valley_rer_result res;
    for (std::size_t i = 0; i < n; ++i) if (in_dom[i]) ++res.n_domain_pixels;

    // Floor half-width: distance from a floor pixel to the edge of the floor
    // (0 off the floor).
    std::vector<float> wd(n, 0.0f);
    {
        std::vector<float> src(n, 0.0f);
        for (std::size_t i = 0; i < n; ++i) if (dmask[i] && !is_floor[i]) src[i] = 1.0f;
        const lsdtt_raster_cpu dist = euclidean_distance_transform(
            lsdtt_raster_cpu(std::move(src), std::vector<std::uint8_t>(dmask), meta));
        const std::vector<float>& dd = dist.get_data();
        const std::vector<std::uint8_t>& dm = dist.get_mask();
        for (std::size_t i = 0; i < n; ++i)
            if (is_floor[i] && dm[i] && std::isfinite(dd[i])) wd[i] = dd[i];
    }

    // -------------------------------------------------------------------------
    // 1. Floor discs: larger channel wins on the floor.
    // -------------------------------------------------------------------------
    const float no_pri = -std::numeric_limits<float>::infinity();
    std::vector<std::int32_t> best_c(nf, -1), best_k(nf, -1);
    std::vector<float>        best_pri(nf, no_pri), best_d2(nf, 0.0f);
    std::vector<std::int32_t> tmp_k(nf, -1);
    std::vector<float>        tmp_d2(nf, 0.0f);
    std::vector<std::int32_t> touched;

    for (std::size_t ci = 0; ci < centrelines.size(); ++ci)
    {
        const valley_rer_centreline& c = centrelines[ci];
        const std::size_t m = c.flats.size();
        if (m == 0) continue;
        touched.clear();
        // A disc at every node. Spacing them out (a quarter of the local
        // radius) was tried and halved the time on a wide floodplain, but it
        // trims the floor's edges between centres and those cells then fall to
        // the ground-distance rule: 3-5 % of cells changed, by up to ~8 m where
        // they landed on another centreline (2026-10-30). Kept exact.
        for (std::size_t k = 0; k < m; ++k)
        {
            const std::int64_t f = c.flats[k];
            const double r = static_cast<double>(wd[static_cast<std::size_t>(f)]) + cell;
            const std::size_t k_here = k;
            if (!is_floor[static_cast<std::size_t>(f)]) continue;  // no floor, no disc
            const int rc = static_cast<int>(std::ceil(r / cell));
            const int r0 = static_cast<int>(f / n_cols), c0 = static_cast<int>(f % n_cols);
            for (int rr = std::max(0, r0 - rc); rr <= std::min(n_rows - 1, r0 + rc); ++rr)
                for (int cc = std::max(0, c0 - rc); cc <= std::min(n_cols - 1, c0 + rc); ++cc)
                {
                    const std::size_t g = static_cast<std::size_t>(rr) * n_cols + cc;
                    const std::int32_t v = fid[g];
                    if (v < 0) continue;
                    const double d2 = (static_cast<double>(rr - r0) * (rr - r0) +
                                       static_cast<double>(cc - c0) * (cc - c0)) * cell * cell;
                    if (d2 > r * r) continue;
                    if (tmp_k[v] < 0) touched.push_back(v);
                    if (tmp_k[v] < 0 || d2 < tmp_d2[v])
                    {
                        tmp_k[v] = static_cast<std::int32_t>(k_here);
                        tmp_d2[v] = static_cast<float>(d2);
                    }
                }
        }
        for (const std::int32_t v : touched)
        {
            const float pri = static_cast<float>(c.priority[static_cast<std::size_t>(tmp_k[v])]);
            if (pri > best_pri[v] || (pri == best_pri[v] && tmp_d2[v] < best_d2[v]))
            {
                best_pri[v] = pri;
                best_d2[v] = tmp_d2[v];
                best_c[v] = static_cast<std::int32_t>(ci);
                best_k[v] = tmp_k[v];
            }
            tmp_k[v] = -1;
        }
    }
    for (std::size_t v = 0; v < nf; ++v)
    {
        if (best_c[v] < 0) continue;
        const valley_rer_centreline& c = centrelines[static_cast<std::size_t>(best_c[v])];
        const int s = 1;   // tmp_k is the nearest covering node; check its neighbours
        const pt p = flat_pt(meta, fflat[v]);
        const std::int64_t k0 = std::max<std::int64_t>(0, best_k[v] - s);
        const std::int64_t k1 = std::min<std::int64_t>(static_cast<std::int64_t>(c.flats.size()) - 1,
                                                       best_k[v] + s);
        double bd = std::numeric_limits<double>::infinity();
        for (std::int64_t k = k0; k <= k1; ++k)
        {
            const pt q = flat_pt(meta, c.flats[static_cast<std::size_t>(k)]);
            const double d = (q.x - p.x) * (q.x - p.x) + (q.y - p.y) * (q.y - p.y);
            if (d < bd) { bd = d; best_k[v] = static_cast<std::int32_t>(k); }
        }
    }

    // -------------------------------------------------------------------------
    // 2. Scaled ground distance through the domain from the kept nodes: cost
    //    of a step = step length / the floor half-width of the node the front
    //    came from (at least one cell). Metres along the winning path are
    //    tracked alongside: the distance cap is applied AFTER the partition, so
    //    a cell too far from the node it belongs to is left out rather than
    //    handed to some other, nearer node.
    // -------------------------------------------------------------------------
    std::vector<float> gd(n, std::numeric_limits<float>::infinity());
    std::vector<float> gm(n, std::numeric_limits<float>::infinity());
    std::vector<std::int64_t> glab(n, -1);   // ci << 32 | k
    std::vector<float> src_pri(n, no_pri);
    std::vector<std::int64_t> sources;
    std::vector<std::vector<float>> scale(centrelines.size());
    for (std::size_t ci = 0; ci < centrelines.size(); ++ci)
    {
        const valley_rer_centreline& c = centrelines[ci];
        scale[ci].resize(c.flats.size());
        for (std::size_t k = 0; k < c.flats.size(); ++k)
        {
            const std::size_t f = static_cast<std::size_t>(c.flats[k]);
            scale[ci][k] = static_cast<float>(std::max(static_cast<double>(wd[f]), cell));
            if (!in_dom[f]) continue;
            const std::int32_t v = fid[f];
            if (v >= 0 && best_c[v] >= 0 && best_c[v] != static_cast<std::int32_t>(ci))
            {
                ++res.n_nodes_dropped;
                continue;
            }
            const float pri = static_cast<float>(c.priority[k]);
            if (glab[f] < 0) sources.push_back(static_cast<std::int64_t>(f));
            if (glab[f] < 0 || pri > src_pri[f])
            {
                glab[f] = (static_cast<std::int64_t>(ci) << 32) | static_cast<std::int64_t>(k);
                src_pri[f] = pri;
                gd[f] = 0.0f;
                gm[f] = 0.0f;
            }
        }
    }
    const double cap = request.max_distance_metres > 0.0
        ? request.max_distance_metres : std::numeric_limits<double>::infinity();
    const float diag = static_cast<float>(cell * grid::k_sqrt2_d);
    const float straight = static_cast<float>(cell);
    auto step_of = [&](int k) {
        return (grid::d8_rowmajor_row[k] != 0 && grid::d8_rowmajor_col[k] != 0)
            ? diag : straight; };
    grid_dijkstra_d8(n_rows, n_cols, gd, sources,
        [&](std::size_t to) { return in_dom[to] != 0; },
        [&](std::size_t from, std::size_t, int k) {
            const float st = step_of(k);
            const std::int64_t lab = glab[from];
            const float sc = scale[static_cast<std::size_t>(lab >> 32)]
                                  [static_cast<std::size_t>(lab & 0xffffffffLL)];
            return st / sc; },
        [&](std::size_t from, std::size_t to) {
            glab[to] = glab[from];
            const int dr = static_cast<int>(to / n_cols) - static_cast<int>(from / n_cols);
            const int dc = static_cast<int>(to % n_cols) - static_cast<int>(from % n_cols);
            gm[to] = gm[from] + ((dr != 0 && dc != 0) ? diag : straight); },
        [](std::int64_t) { return false; });

    // -------------------------------------------------------------------------
    // Relative elevation.
    // -------------------------------------------------------------------------
    const float ndv = static_cast<float>(meta.no_data_value);
    std::vector<float> rer(n, ndv), own(n, ndv);
    std::vector<std::uint8_t> rm(n, 0), om(n, 0);
    res.owner_centreline.assign(n, -1);
    res.owner_node.assign(n, -1);
    for (std::size_t f = 0; f < n; ++f)
    {
        if (!in_dom[f]) continue;
        std::int64_t ci = -1, k = -1;
        const std::int32_t v = fid[f];
        if (v >= 0 && best_c[v] >= 0) { ci = best_c[v]; k = best_k[v]; ++res.n_by_disc; }
        else if (glab[f] >= 0)
        {
            if (static_cast<double>(gm[f]) > cap) { ++res.n_unassigned; continue; }
            ci = glab[f] >> 32;
            k = glab[f] & 0xffffffffLL;
            ++res.n_by_geodesic;
        }
        if (ci < 0) { ++res.n_unassigned; continue; }
        res.owner_centreline[f] = static_cast<std::int32_t>(ci);
        res.owner_node[f] = static_cast<std::int32_t>(k);
        const valley_rer_centreline& c = centrelines[static_cast<std::size_t>(ci)];
        own[f] = static_cast<float>(c.id);
        om[f] = 1;
        const float b = c.baseline[static_cast<std::size_t>(k)];
        if (std::isfinite(b) && std::isfinite(elev[f]))
        {
            rer[f] = elev[f] - b;
            rm[f] = 1;
        }
    }
    res.relative_elevation = lsdtt_raster_cpu(std::move(rer), std::move(rm), meta);
    res.owner = lsdtt_raster_cpu(std::move(own), std::move(om), meta);
    return res;
}

} // namespace lsdtt
