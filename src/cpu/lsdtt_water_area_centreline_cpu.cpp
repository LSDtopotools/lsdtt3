// lsdtt_water_area_centreline_cpu.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html

#include "lsdtt_water_area_centreline_cpu.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <functional>
#include <map>
#include <queue>
#include <set>
#include <stdexcept>
#include <utility>
#include <variant>
#include <vector>

#include "lsdtt_raster_cpu.hpp"
#include "lsdtt_raster_skeleton_graph_cpu.hpp"
#include "lsdtt_raster_binary_mask_cpu.hpp"
#include "lsdtt_proximity_cpu.hpp"
#include "../core/lsdtt_raster_meta.hpp"
#include "../utils/lsdtt_raster_geometry.hpp"
#include "../utils/lsdtt_polygon_geometry.hpp"

namespace lsdtt
{
namespace
{

constexpr double k_inf = std::numeric_limits<double>::infinity();

// --- read provenance fields from an input record ---------------------------
std::string read_string_field(
    const lsdtt_vector_layer&  layer,
    const lsdtt_vector_record& rec,
    const std::string&         name)
{
    const int idx = layer.schema.index_of(name);
    if (idx < 0 || idx >= static_cast<int>(rec.field_values.size())) return "";
    if (const auto* s =
            std::get_if<std::string>(&rec.field_values[static_cast<std::size_t>(idx)]))
        return *s;
    return "";
}

std::int64_t read_int_field(
    const lsdtt_vector_layer&  layer,
    const lsdtt_vector_record& rec,
    const std::string&         name)
{
    const int idx = layer.schema.index_of(name);
    if (idx < 0 || idx >= static_cast<int>(rec.field_values.size())) return 0;
    if (const auto* v =
            std::get_if<std::int64_t>(&rec.field_values[static_cast<std::size_t>(idx)]))
        return *v;
    return 0;
}

struct bbox2d { double xmin, ymin, xmax, ymax; };

bbox2d bbox_of_ring(const std::vector<lsdtt_xy>& ring)
{
    bbox2d b{k_inf, k_inf, -k_inf, -k_inf};
    for (const lsdtt_xy& p : ring)
    {
        b.xmin = std::min(b.xmin, p.x);
        b.ymin = std::min(b.ymin, p.y);
        b.xmax = std::max(b.xmax, p.x);
        b.ymax = std::max(b.ymax, p.y);
    }
    return b;
}

// Build a synthesised reference grid covering bbox + padding at cell_size.
// Auto-coarsens cell_size if the grid would exceed a cell cap.
lsdtt_raster_meta make_grid_meta(
    const bbox2d& bbox, double cell_size, double padding, int epsg)
{
    double cs = (cell_size > 0.0) ? cell_size : 1.0;
    const double w = (bbox.xmax - bbox.xmin) + 2.0 * padding;
    const double h = (bbox.ymax - bbox.ymin) + 2.0 * padding;

    const long long cap = 50000000LL;   // ~50 M cells safety cap
    for (int guard = 0; guard < 32; ++guard)
    {
        const long long nc = static_cast<long long>(std::ceil(w / cs));
        const long long nr = static_cast<long long>(std::ceil(h / cs));
        if (nc * nr <= cap || nc < 1 || nr < 1) break;
        cs *= 1.5;   // coarsen and retry
    }

    lsdtt_raster_meta meta;
    meta.cell_size     = cs;
    meta.n_cols        = std::max(1, static_cast<int>(std::ceil(w / cs)));
    meta.n_rows        = std::max(1, static_cast<int>(std::ceil(h / cs)));
    meta.x_min         = bbox.xmin - padding;
    meta.y_min         = bbox.ymin - padding;
    meta.no_data_value = -9999.0;
    try { meta.set_crs_from_epsg(epsg); } catch (...) { /* leave CRS empty */ }
    return meta;
}

// Scanline-rasterise a polygon (exterior + holes) into a binary mask aligned to
// meta. data = 1 inside (even-odd rule over all rings), 0 outside; mask = 1
// everywhere (the synthesised grid is the valid domain).
lsdtt_raster_cpu rasterise_polygon_scanline(
    const lsdtt_raster_meta&  meta,
    const lsdtt_polygon_geom& poly)
{
    const int nr = meta.n_rows;
    const int nc = meta.n_cols;
    const double cs = meta.cell_size;

    lsdtt_array2d_f  data(nr, nc, 0.0f);
    lsdtt_array2d_u8 mask(nr, nc, std::uint8_t{1});

    // Gather all rings (exterior + holes).
    std::vector<const std::vector<lsdtt_xy>*> rings;
    rings.push_back(&poly.exterior);
    for (const auto& h : poly.holes) rings.push_back(&h);

    std::vector<double> xs;
    for (int row = 0; row < nr; ++row)
    {
        const double yc = cell_centre_y(meta, row);
        xs.clear();

        for (const std::vector<lsdtt_xy>* ring_ptr : rings)
        {
            const std::vector<lsdtt_xy>& ring = *ring_ptr;
            const std::size_t n = ring.size();
            if (n < 3) continue;
            for (std::size_t i = 0; i < n; ++i)
            {
                const lsdtt_xy& a = ring[i];
                const lsdtt_xy& b = ring[(i + 1) % n];
                // Half-open crossing test avoids double-counting shared vertices.
                const bool crosses =
                    (a.y <= yc && yc < b.y) || (b.y <= yc && yc < a.y);
                if (!crosses) continue;
                const double t = (yc - a.y) / (b.y - a.y);
                xs.push_back(a.x + t * (b.x - a.x));
            }
        }

        if (xs.size() < 2) continue;
        std::sort(xs.begin(), xs.end());

        for (std::size_t k = 0; k + 1 < xs.size(); k += 2)
        {
            const double x0 = xs[k];
            const double x1 = xs[k + 1];
            // cell centre x = x_min + (col + 0.5) * cs  ->  col in [lo, hi]
            int col_lo = static_cast<int>(std::ceil((x0 - meta.x_min) / cs - 0.5));
            int col_hi = static_cast<int>(std::floor((x1 - meta.x_min) / cs - 0.5));
            col_lo = std::max(0, col_lo);
            col_hi = std::min(nc - 1, col_hi);
            for (int col = col_lo; col <= col_hi; ++col)
                data(row, col) = 1.0f;
        }
    }

    return lsdtt_raster_cpu(std::move(data), std::move(mask), meta);
}

// Open-polyline Douglas-Peucker (does NOT close the ring, unlike the polygon
// helper). Keeps endpoints; recursively keeps the farthest point above tol.
void dp_open_recurse(
    const std::vector<lsdtt_xy>& pts, std::size_t i0, std::size_t i1,
    double tol, std::vector<char>& keep)
{
    if (i1 <= i0 + 1) return;
    double dmax = -1.0;
    std::size_t idx = i0;
    for (std::size_t i = i0 + 1; i < i1; ++i)
    {
        const double d = point_to_segment_distance(
            pts[i].x, pts[i].y, pts[i0].x, pts[i0].y, pts[i1].x, pts[i1].y);
        if (d > dmax) { dmax = d; idx = i; }
    }
    if (dmax > tol)
    {
        keep[idx] = 1;
        dp_open_recurse(pts, i0, idx, tol, keep);
        dp_open_recurse(pts, idx, i1, tol, keep);
    }
}

// Laplacian smoothing of an open polyline (endpoints fixed). Each interior
// vertex moves toward the average of its neighbours, removing the pixel-skeleton
// staircase wobble while preserving the overall path.
std::vector<lsdtt_xy> smooth_open_polyline(
    const std::vector<lsdtt_xy>& in, int iterations)
{
    if (iterations <= 0 || in.size() < 3) return in;
    std::vector<lsdtt_xy> p = in;
    for (int it = 0; it < iterations; ++it)
    {
        std::vector<lsdtt_xy> q = p;
        for (std::size_t i = 1; i + 1 < p.size(); ++i)
        {
            q[i].x = 0.5 * p[i].x + 0.25 * (p[i - 1].x + p[i + 1].x);
            q[i].y = 0.5 * p[i].y + 0.25 * (p[i - 1].y + p[i + 1].y);
        }
        p.swap(q);
    }
    return p;
}

std::vector<lsdtt_xy> simplify_open_polyline(
    const std::vector<lsdtt_xy>& pts, double tol)
{
    const std::size_t n = pts.size();
    if (tol <= 0.0 || n < 3) return pts;
    std::vector<char> keep(n, 0);
    keep[0]     = 1;
    keep[n - 1] = 1;
    dp_open_recurse(pts, 0, n - 1, tol, keep);
    std::vector<lsdtt_xy> out;
    out.reserve(n);
    for (std::size_t i = 0; i < n; ++i)
        if (keep[i]) out.push_back(pts[i]);
    return out;
}

// --- centreline chains (graph diameter + narrow channel spurs) --------------
// Returns the water body's centreline as ordered pixel chains: the skeleton's
// longest path (graph diameter = one clean long-axis line) plus, when
// detect_spurs is set, narrow elongated branches (river channels entering or
// leaving the body) identified by a small local half-width (EDT) at the branch
// tip. Wide bays are ignored. Returns the longest single segment if there are
// no edges, or empty if the skeleton is empty.
std::vector<std::vector<skeleton_pixel>> extract_centreline_chains(
    const std::vector<skeleton_segment>& segs,
    const std::vector<float>&            edt_half_width,   // metres, row-major
    int                                  n_cols,
    bool                                 detect_spurs,
    double                               channel_max_half_width_m,
    double                               min_spur_len_m,
    double                               cell_size_m)
{
    std::vector<std::vector<skeleton_pixel>> chains;

    std::map<std::pair<int,int>, int> nid;
    std::vector<std::pair<int,int>>   rc;
    auto node = [&](int r, int c) -> int
    {
        const std::pair<int,int> k{r, c};
        auto it = nid.find(k);
        if (it != nid.end()) return it->second;
        const int id = static_cast<int>(rc.size());
        nid.emplace(k, id);
        rc.push_back(k);
        return id;
    };

    struct E { int a; int b; std::size_t seg; double w; };
    std::vector<E> edges;
    for (std::size_t si = 0; si < segs.size(); ++si)
    {
        const skeleton_segment& s = segs[si];
        if (s.is_loop || s.pixels.size() < 2) continue;
        const int a = node(s.start_node.row, s.start_node.col);
        const int b = node(s.end_node.row,   s.end_node.col);
        if (a == b) continue;
        edges.push_back(E{a, b, si, static_cast<double>(s.pixels.size())});
    }

    if (edges.empty())
    {
        std::size_t best = segs.size(), bn = 0;
        for (std::size_t si = 0; si < segs.size(); ++si)
            if (segs[si].pixels.size() > bn) { bn = segs[si].pixels.size(); best = si; }
        if (best < segs.size()) chains.push_back(segs[best].pixels);
        return chains;
    }

    std::vector<std::vector<std::pair<int,int>>> adj(rc.size());
    for (int e = 0; e < static_cast<int>(edges.size()); ++e)
    {
        adj[static_cast<std::size_t>(edges[e].a)].push_back({edges[e].b, e});
        adj[static_cast<std::size_t>(edges[e].b)].push_back({edges[e].a, e});
    }

    auto dij = [&](int src, std::vector<double>& dist,
                   std::vector<int>& pn, std::vector<int>& pe)
    {
        dist.assign(rc.size(), k_inf);
        pn.assign(rc.size(), -1);
        pe.assign(rc.size(), -1);
        using qe = std::pair<double,int>;
        std::priority_queue<qe, std::vector<qe>, std::greater<qe>> pq;
        dist[static_cast<std::size_t>(src)] = 0.0;
        pq.push({0.0, src});
        while (!pq.empty())
        {
            const auto [d, u] = pq.top(); pq.pop();
            if (d > dist[static_cast<std::size_t>(u)]) continue;
            for (const auto& [v, e] : adj[static_cast<std::size_t>(u)])
            {
                const double nd = d + edges[static_cast<std::size_t>(e)].w;
                if (nd < dist[static_cast<std::size_t>(v)])
                {
                    dist[static_cast<std::size_t>(v)] = nd;
                    pn[static_cast<std::size_t>(v)] = u;
                    pe[static_cast<std::size_t>(v)] = e;
                    pq.push({nd, v});
                }
            }
        }
    };
    auto farthest = [&](const std::vector<double>& d) -> int
    {
        int bi = -1; double bd = -1.0;
        for (int i = 0; i < static_cast<int>(d.size()); ++i)
            if (std::isfinite(d[i]) && d[i] > bd) { bd = d[i]; bi = i; }
        return bi;
    };

    // Reconstruct the node/edge path src..dst into an ordered pixel chain,
    // recording visited nodes into `visited`.
    auto reconstruct = [&](int src, int dst, const std::vector<int>& pn,
                           const std::vector<int>& pe, std::set<int>& visited)
        -> std::vector<skeleton_pixel>
    {
        std::vector<int> seq, es;
        for (int cur = dst; cur != src; )
        {
            const int e = pe[static_cast<std::size_t>(cur)];
            const int p = pn[static_cast<std::size_t>(cur)];
            if (e < 0 || p < 0) return {};
            seq.push_back(cur); es.push_back(e); cur = p;
        }
        seq.push_back(src);
        std::reverse(seq.begin(), seq.end());
        std::reverse(es.begin(), es.end());

        std::vector<skeleton_pixel> chain;
        for (std::size_t i = 0; i < es.size(); ++i)
        {
            const skeleton_segment& s = segs[edges[static_cast<std::size_t>(es[i])].seg];
            const std::pair<int,int>& urc = rc[static_cast<std::size_t>(seq[i])];
            const bool fwd =
                (s.start_node.row == urc.first && s.start_node.col == urc.second);
            if (fwd)
                for (std::size_t k = 0; k < s.pixels.size(); ++k)
                { if (!chain.empty() && k == 0) continue; chain.push_back(s.pixels[k]); }
            else
                for (std::size_t k = s.pixels.size(); k-- > 0; )
                { if (!chain.empty() && k == s.pixels.size() - 1) continue;
                  chain.push_back(s.pixels[k]); }
        }
        for (int n : seq) visited.insert(n);
        return chain;
    };

    // Main stems = the diameter of EACH connected component. A large lake's
    // thinned skeleton can fragment into several components (narrow necks,
    // islands); a single Dijkstra from one node would cover only that node's
    // component and leave the rest of the lake with no centreline.
    std::set<int> retained;
    {
        std::vector<char> seen(rc.size(), 0);
        for (int s = 0; s < static_cast<int>(rc.size()); ++s)
        {
            if (seen[static_cast<std::size_t>(s)]) continue;
            std::vector<double> ds; std::vector<int> sn, se;
            dij(s, ds, sn, se);                       // reaches s's component
            for (int i = 0; i < static_cast<int>(rc.size()); ++i)
                if (std::isfinite(ds[static_cast<std::size_t>(i)]))
                    seen[static_cast<std::size_t>(i)] = 1;
            const int A = farthest(ds);               // farthest in component
            if (A < 0) { seen[static_cast<std::size_t>(s)] = 1; continue; }
            std::vector<double> da; std::vector<int> pan, pae;
            dij(A, da, pan, pae);
            const int B = farthest(da);               // component diameter A..B
            if (B < 0) continue;
            std::vector<skeleton_pixel> stem = reconstruct(A, B, pan, pae, retained);
            if (stem.size() >= 2) chains.push_back(std::move(stem));
        }
    }
    if (chains.empty()) return chains;

    // Narrow channel spurs: retain elongated branches whose tip half-width is
    // small (a river channel), ignoring wide bays. Iterate to a fixpoint so
    // nested arms connect once their target branch is retained.
    if (detect_spurs && !edt_half_width.empty())
    {
        auto hw_at = [&](int nodeid) -> double
        {
            const std::pair<int,int>& p = rc[static_cast<std::size_t>(nodeid)];
            const std::size_t idx =
                static_cast<std::size_t>(p.first) * static_cast<std::size_t>(n_cols)
                + static_cast<std::size_t>(p.second);
            return (idx < edt_half_width.size())
                ? static_cast<double>(edt_half_width[idx]) : k_inf;
        };

        bool added = true;
        int  guard = 0;
        while (added && guard++ < 16)
        {
            added = false;
            for (int v = 0; v < static_cast<int>(rc.size()); ++v)
            {
                if (adj[static_cast<std::size_t>(v)].size() != 1) continue; // tips
                if (retained.count(v)) continue;
                if (hw_at(v) > channel_max_half_width_m) continue;          // wide bay
                std::vector<double> dv; std::vector<int> pvn, pve;
                dij(v, dv, pvn, pve);
                int    tgt = -1; double tb = k_inf;
                for (int rn : retained)
                    if (std::isfinite(dv[static_cast<std::size_t>(rn)]) &&
                        dv[static_cast<std::size_t>(rn)] < tb)
                    { tb = dv[static_cast<std::size_t>(rn)]; tgt = rn; }
                if (tgt < 0) continue;
                if (tb * cell_size_m < min_spur_len_m) continue;            // too short
                std::vector<skeleton_pixel> ch =
                    reconstruct(v, tgt, pvn, pve, retained);
                if (ch.size() >= 2) { chains.push_back(std::move(ch)); added = true; }
            }
        }
    }

    return chains;
}

// --- clip water lines against water-area polygons ---------------------------
struct clipped_trib
{
    std::vector<lsdtt_xy> xy;
    std::string           fid;
    std::int64_t          feat = 0;
    bool                  front_cut = false;  // front end was cut at a water area
    bool                  back_cut  = false;  // back end was cut at a water area
};

// Drop Water_Line portions that are the outline of (or run inside) a water
// area, keeping maximal runs of segments whose midpoint is neither inside a
// polygon nor within tol of any polygon ring. Each kept run becomes a tributary.
std::vector<clipped_trib> clip_water_lines(
    const lsdtt_vector_layer&              water_lines,
    const std::vector<lsdtt_polygon_geom>& polys,
    double                                 tol,
    double                                 min_len)
{
    // Closed rings (exterior + holes) for boundary-distance tests.
    std::vector<std::vector<lsdtt_xy>> rings;
    for (const lsdtt_polygon_geom& poly : polys)
    {
        auto add_closed = [&](const std::vector<lsdtt_xy>& r)
        {
            if (r.size() < 2) return;
            std::vector<lsdtt_xy> c = r;
            if (!(c.front().x == c.back().x && c.front().y == c.back().y))
                c.push_back(c.front());
            rings.push_back(std::move(c));
        };
        add_closed(poly.exterior);
        for (const auto& h : poly.holes) add_closed(h);
    }

    auto in_or_near_water = [&](const lsdtt_xy& p) -> bool
    {
        for (const lsdtt_polygon_geom& poly : polys)
            if (point_in_polygon(p, poly.exterior, poly.holes)) return true;
        for (const std::vector<lsdtt_xy>& ring : rings)
            if (nearest_point_on_polyline(p, ring).distance <= tol) return true;
        return false;
    };

    std::vector<clipped_trib> out;
    for (const lsdtt_vector_record& rec : water_lines.records)
    {
        const auto* line = std::get_if<lsdtt_line_geom>(&rec.geometry);
        if (!line || line->points.size() < 2) continue;
        const std::vector<lsdtt_xy>& P = line->points;
        const std::string  fid  = read_string_field(water_lines, rec, "Fid");
        const std::int64_t feat = read_int_field(water_lines, rec, "FeatCode");

        const std::size_t nseg = P.size() - 1;
        std::size_t i = 0;
        while (i < nseg)
        {
            const lsdtt_xy mid{0.5 * (P[i].x + P[i + 1].x),
                               0.5 * (P[i].y + P[i + 1].y)};
            if (in_or_near_water(mid)) { ++i; continue; }   // outline/inside: drop

            const std::size_t start_i = i;
            clipped_trib t;
            t.fid = fid;
            t.feat = feat;
            t.front_cut = (start_i > 0);   // a dropped segment precedes this run
            t.xy.push_back(P[i]);
            while (i < nseg)
            {
                const lsdtt_xy m{0.5 * (P[i].x + P[i + 1].x),
                                 0.5 * (P[i].y + P[i + 1].y)};
                if (in_or_near_water(m)) break;
                t.xy.push_back(P[i + 1]);
                ++i;
            }
            t.back_cut = (i < nseg);       // stopped because the next segment is water
            if (t.xy.size() >= 2 && polyline_length(t.xy) >= min_len)
                out.push_back(std::move(t));
        }
    }
    return out;
}

} // anonymous namespace


// ===========================================================================
// extract_water_area_centrelines
// ===========================================================================
water_centreline_result extract_water_area_centrelines(
    const lsdtt_vector_layer&       water_areas,
    const lsdtt_vector_layer&       water_lines,
    const water_centreline_options& opt)
{
    if (opt.cell_size_m <= 0.0)
        throw std::runtime_error(
            "extract_water_area_centrelines: cell_size_m must be > 0.");

    water_centreline_result result;

    // Build polygon list (exterior + holes) once.
    std::vector<lsdtt_polygon_geom> polys;
    for (const lsdtt_vector_record& rec : water_areas.records)
    {
        const auto* poly = std::get_if<lsdtt_polygon_geom>(&rec.geometry);
        if (poly && poly->exterior.size() >= 3) polys.push_back(*poly);
    }

    // Clip Water_Line against Water_Area: drop the outline/inside portions,
    // keep the tributary runs outside the water bodies. These clipped tributaries
    // drive anchoring, connectors, and the pass-through output.
    std::vector<clipped_trib> tribs;
    if (opt.clip_lines_to_areas)
    {
        tribs = clip_water_lines(water_lines, polys,
                                 opt.line_overlap_tolerance_m,
                                 opt.min_tributary_length_m);
    }
    else
    {
        for (const lsdtt_vector_record& rec : water_lines.records)
        {
            const auto* line = std::get_if<lsdtt_line_geom>(&rec.geometry);
            if (!line || line->points.size() < 2) continue;
            clipped_trib t;
            t.xy   = line->points;
            t.fid  = read_string_field(water_lines, rec, "Fid");
            t.feat = read_int_field(water_lines, rec, "FeatCode");
            tribs.push_back(std::move(t));
        }
    }

    std::int64_t seg_id = 0;

    // Collected centreline geometries (+ source polygon Fid) for the connector
    // pass after all polygons are processed.
    std::vector<std::pair<std::vector<lsdtt_xy>, std::string>> all_centrelines;

    // Output schema.
    result.output_layer.schema.field_names = {
        "role", "source_fid", "feat_code", "length_m", "segment_id"
    };
    result.output_layer.schema.field_types = {
        lsdtt_vector_field_type::string_type,
        lsdtt_vector_field_type::string_type,
        lsdtt_vector_field_type::int64_type,
        lsdtt_vector_field_type::double_type,
        lsdtt_vector_field_type::int64_type
    };

    auto push_segment = [&](const water_centreline_segment& s, const char* role_str)
    {
        lsdtt_vector_record rec;
        rec.geometry = lsdtt_line_geom{s.xy};
        rec.field_values = {
            std::string(role_str),
            s.source_fid.empty() ? lsdtt_vector_field_value{std::monostate{}}
                                 : lsdtt_vector_field_value{s.source_fid},
            s.feat_code,
            s.length_m,
            s.segment_id
        };
        result.output_layer.records.push_back(std::move(rec));
    };

    // -----------------------------------------------------------------------
    // Per polygon
    // -----------------------------------------------------------------------
    for (const lsdtt_vector_record& area_rec : water_areas.records)
    {
        const auto* poly = std::get_if<lsdtt_polygon_geom>(&area_rec.geometry);
        if (!poly || poly->exterior.size() < 3) continue;

        ++result.n_polygons;

        water_centreline_polygon_result pr;
        pr.source_fid = read_string_field(water_areas, area_rec, "Fid");
        const std::int64_t feat_code =
            read_int_field(water_areas, area_rec, "FeatCode");

        // Loch detection from Feret aspect of the exterior hull.
        const std::vector<lsdtt_xy> hull = convex_hull(poly->exterior);
        const lsdtt_feret_diameters fer = feret_diameters(hull);
        const double aspect =
            (fer.d_min > 1e-9) ? (fer.d_max / fer.d_min) : k_inf;
        pr.treated_as_loch = (aspect < opt.max_centreline_aspect_for_loch);
        if (pr.treated_as_loch) ++result.n_lochs;

        std::vector<water_centreline_segment> centrelines;

        const bool make_centreline =
            !(pr.treated_as_loch && opt.skip_loch_centreline);

        if (make_centreline)
        {
            // 1-2. grid + binary mask.
            const bbox2d bb = bbox_of_ring(poly->exterior);
            const lsdtt_raster_meta meta =
                make_grid_meta(bb, opt.cell_size_m, opt.padding_m, opt.epsg);
            pr.n_grid_rows = meta.n_rows;
            pr.n_grid_cols = meta.n_cols;

            lsdtt_raster_cpu mask = rasterise_polygon_scanline(meta, *poly);
            if (opt.fill_holes) mask = binary_fill_holes(mask);

            // 3. skeletonise.
            lsdtt_raster_cpu skel =
                thin_distance_guided(mask, opt.preserve_distance_pixels);
            if (opt.max_spur_pixels > 0)
                skel = prune_skeleton_spurs(skel, opt.max_spur_pixels);
            if (opt.max_loop_pixels > 0)
                skel = prune_terminal_loops(skel, opt.max_loop_pixels);

            std::vector<skeleton_segment> segs = extract_skeleton_segments(skel);
            pr.n_skeleton_segments = static_cast<int>(segs.size());

            // 4. Centreline = the skeleton's longest path (graph diameter) plus
            //    narrow channel spurs. The latter use the EDT half-width
            //    (distance to bank): build an inverted mask (outside = source)
            //    so the EDT at each interior pixel is its local half-width.
            std::vector<float> edt_hw;
            if (opt.detect_channel_spurs)
            {
                const std::vector<float>& md = mask.get_data();
                std::vector<float>        inv(md.size());
                std::vector<std::uint8_t> invm(md.size(), std::uint8_t{1});
                for (std::size_t i = 0; i < md.size(); ++i)
                    inv[i] = (md[i] > 0.5f) ? 0.0f : 1.0f;   // outside = source
                lsdtt_raster_cpu inv_r(std::move(inv), std::move(invm), meta);
                edt_hw = euclidean_distance_transform(inv_r).get_data();
            }

            std::vector<std::vector<skeleton_pixel>> retained_chains =
                extract_centreline_chains(
                    segs, edt_hw, meta.n_cols,
                    opt.detect_channel_spurs,
                    opt.channel_max_half_width_m,
                    opt.min_channel_spur_len_m,
                    meta.cell_size);

            // 5. pixel chains -> world polylines -> simplify.
            for (const std::vector<skeleton_pixel>& chain : retained_chains)
            {
                if (chain.size() < 2) continue;
                skeleton_segment tmp;
                tmp.pixels = chain;
                std::vector<lsdtt_xy> poly_xy =
                    skeleton_segment_to_polyline(tmp, meta);
                poly_xy = smooth_open_polyline(poly_xy, opt.centreline_smooth_iterations);
                poly_xy = simplify_open_polyline(poly_xy, opt.dp_tolerance_m);
                if (poly_xy.size() < 2) continue;

                water_centreline_segment cs;
                cs.xy = std::move(poly_xy);
                cs.length_m = polyline_length(cs.xy);
                if (cs.length_m < opt.min_centreline_component_m) continue; // declutter
                cs.role = water_centreline_role::centreline;
                cs.source_fid = pr.source_fid;
                cs.feat_code = feat_code;
                cs.segment_id = seg_id++;
                centrelines.push_back(std::move(cs));
            }

            // Guarantee a line for EVERY polygon: if the skeleton produced
            // nothing (tiny / degenerate water body), fall back to the polygon's
            // long axis (min-area-rectangle major axis through the centroid).
            if (centrelines.empty())
            {
                const std::vector<lsdtt_xy> fallback_hull = convex_hull(poly->exterior);
                const lsdtt_min_area_rect rect = min_area_bounding_rect(fallback_hull);
                const lsdtt_xy c = polygon_centroid(poly->exterior);
                const double hl = 0.5 * rect.length;
                const double ct = std::cos(rect.theta_rad);
                const double st = std::sin(rect.theta_rad);
                water_centreline_segment cs;
                cs.xy = { lsdtt_xy{c.x - hl * ct, c.y - hl * st},
                          lsdtt_xy{c.x + hl * ct, c.y + hl * st} };
                cs.role = water_centreline_role::centreline;
                cs.source_fid = pr.source_fid;
                cs.feat_code = feat_code;
                cs.length_m = polyline_length(cs.xy);
                cs.segment_id = seg_id++;
                centrelines.push_back(std::move(cs));
            }
        }

        pr.centrelines = centrelines;
        for (const auto& c : pr.centrelines) push_segment(c, "centreline");
        result.n_centrelines += static_cast<int>(pr.centrelines.size());
        for (const auto& c : pr.centrelines)
            all_centrelines.push_back({c.xy, pr.source_fid});

        result.per_polygon.push_back(std::move(pr));
    }

    // -----------------------------------------------------------------------
    // 6. Tributary connectors (spurs). Every tributary mouth - an end that was
    //    cut by a water area, or an end lying within tributary_snap_tolerance of
    //    a water-area boundary - is joined to the nearest centreline point.
    // -----------------------------------------------------------------------
    {
        const double half_cell = 0.5 * opt.cell_size_m;

        std::vector<std::vector<lsdtt_xy>> ext_rings;
        for (const lsdtt_polygon_geom& poly : polys)
        {
            std::vector<lsdtt_xy> c = poly.exterior;
            if (c.size() >= 2 &&
                !(c.front().x == c.back().x && c.front().y == c.back().y))
                c.push_back(c.front());
            ext_rings.push_back(std::move(c));
        }
        auto near_any_boundary = [&](const lsdtt_xy& p) -> bool
        {
            for (const std::vector<lsdtt_xy>& r : ext_rings)
                if (nearest_point_on_polyline(p, r).distance
                        <= opt.tributary_snap_tolerance_m)
                    return true;
            return false;
        };
        auto connect_end = [&](const lsdtt_xy& ep)
        {
            if (all_centrelines.empty()) return;
            lsdtt_xy    best_pt;
            double      best_d = k_inf;
            std::string best_fid;
            for (const auto& cl : all_centrelines)
            {
                const lsdtt_nearest_polyline_point np =
                    nearest_point_on_polyline(ep, cl.first);
                if (np.distance < best_d)
                { best_d = np.distance; best_pt = np.point; best_fid = cl.second; }
            }
            if (best_d > opt.centreline_snap_max_m) return;
            if (best_d <= half_cell) return;   // already on a centreline
            water_centreline_segment conn;
            conn.xy = { ep, best_pt };
            conn.role = water_centreline_role::connector;
            conn.source_fid = best_fid;
            conn.length_m = polyline_length(conn.xy);
            conn.segment_id = seg_id++;
            push_segment(conn, "connector");
            ++result.n_connectors;
        };

        for (const clipped_trib& t : tribs)
        {
            if (t.xy.size() < 2) continue;
            const lsdtt_xy& f = t.xy.front();
            const lsdtt_xy& b = t.xy.back();
            if (t.front_cut || near_any_boundary(f)) connect_end(f);
            if (t.back_cut  || near_any_boundary(b)) connect_end(b);
        }
    }

    // -----------------------------------------------------------------------
    // Pass-through of the CLIPPED tributary lines (outline portions removed).
    // -----------------------------------------------------------------------
    if (opt.emit_unchanged_tributaries)
    {
        for (const clipped_trib& ct : tribs)
        {
            water_centreline_segment t;
            t.xy = ct.xy;
            t.role = water_centreline_role::tributary;
            t.source_fid = ct.fid;
            t.feat_code = ct.feat;
            t.length_m = polyline_length(t.xy);
            t.segment_id = seg_id++;
            push_segment(t, "tributary");
            ++result.n_tributaries;
        }
    }

    return result;
}

} // namespace lsdtt
