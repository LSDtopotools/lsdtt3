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

} // namespace

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

    // Candidate head = the upstream (higher) end of each fragment. `edge` flags a
    // head that is a genuine inflow across the DEM boundary / NoData fringe (the
    // chosen end's original vertex is off-DEM) - these bypass the higher_nearby
    // exclusion below so "channels coming in from the edge" are always seeded.
    struct candidate { lsdtt_xy up; double z; bool edge; };
    std::vector<candidate> cands;
    // All fragment vertices with a valid elevation, for the "is something higher
    // draining toward this head?" test below.
    std::vector<std::pair<lsdtt_xy, double>> verts;
    for (const std::vector<lsdtt_xy>& pts : clipped_lines)
    {
        if (pts.size() < 2) continue;
        // First in-DEM cell walking inward from each end. When an end is itself
        // in-DEM this is just that endpoint; when an end is off-DEM (boundary /
        // NoData fringe) this is the channel's DEM-entry cell.
        lsdtt_xy pf, pb; double zf = 0.0, zb = 0.0;
        const bool okf = first_valid_inward(pts, true,  pf, zf);
        const bool okb = first_valid_inward(pts, false, pb, zb);
        if (!okf && !okb) continue;
        // Was each original end off-DEM? (governs the edge-inflow bypass)
        double ztmp = 0.0;
        const bool front_off = !valid_elev(pts.front(), ztmp);
        const bool back_off  = !valid_elev(pts.back(),  ztmp);
        // Head = the higher of the two inward entries (the upstream end).
        lsdtt_xy up; double zu; bool from_front;
        if (okf && okb)
        { if (zf >= zb) { up = pf; zu = zf; from_front = true; }
          else          { up = pb; zu = zb; from_front = false; } }
        else if (okf)   { up = pf; zu = zf; from_front = true; }
        else            { up = pb; zu = zb; from_front = false; }
        const bool edge = from_front ? front_off : back_off;
        cands.push_back({up, zu, edge});
        for (const lsdtt_xy& p : pts)
        {
            double zp = 0.0;
            if (valid_elev(p, zp)) verts.push_back({p, zp});
        }
    }

    // A candidate is a GENUINE source only if no channel vertex within
    // source_exclusion_radius_m is meaningfully HIGHER - i.e. nothing drains
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
        if (!cand.edge)   // edge inflows are always genuine sources
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
