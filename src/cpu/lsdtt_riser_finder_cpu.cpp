// lsdtt_riser_finder_cpu.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html

#include "lsdtt_riser_finder_cpu.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>

#include "lsdtt_raster_skeleton_graph_cpu.hpp"

namespace lsdtt
{

riser_gap_closing_choice parse_riser_gap_closing_choice(const std::string& name)
{
    if (name == "none")          return riser_gap_closing_choice::none;
    if (name == "astar")         return riser_gap_closing_choice::astar;
    throw std::invalid_argument("unknown riser gap closing '" + name +
                                "'. Valid: none, astar.");
}

namespace
{

// DEM value at (row, col), NaN outside the grid or on nodata.
double sample(const lsdtt_raster_cpu& dem, int r, int c)
{
    const lsdtt_raster_meta& m = dem.get_meta();
    if (r < 0 || c < 0 || r >= m.n_rows || c >= m.n_cols)
        return std::numeric_limits<double>::quiet_NaN();
    const std::size_t i = static_cast<std::size_t>(r) * m.n_cols + c;
    const float v = dem.get_data()[i];
    if (!dem.get_mask()[i] || !std::isfinite(v) || m.is_nodata(v))
        return std::numeric_limits<double>::quiet_NaN();
    return v;
}

// Removes the staircase corners Zhang-Suen leaves on a diagonal line: a pixel
// with both a vertical and a horizontal 4-neighbour whose neighbours stay
// 8-connected without it. Such a pixel has three neighbours but crossing
// number 2, so extract_skeleton_segments never starts a walk from it, and the
// chain between two of them comes back in 2-pixel pieces. Sequential, so of a
// "##" pair only one goes.
void remove_staircase_corners(std::vector<float>& line, int nr, int nc)
{
    static constexpr int dr[8] = {-1, -1, -1, 0, 0, 1, 1, 1};
    static constexpr int dc[8] = {-1, 0, 1, -1, 1, -1, 0, 1};
    auto on = [&](int r, int c)
    {
        return r >= 0 && c >= 0 && r < nr && c < nc &&
               line[static_cast<std::size_t>(r) * nc + c] > 0.5f;
    };
    bool changed = true;
    while (changed)
    {
        changed = false;
        for (int r = 0; r < nr; ++r)
            for (int c = 0; c < nc; ++c)
            {
                if (!on(r, c)) continue;
                const bool vert = on(r - 1, c) || on(r + 1, c);
                const bool horz = on(r, c - 1) || on(r, c + 1);
                if (!vert || !horz) continue;
                // 8-components among the neighbours, by union-find on <= 8.
                int nb[8], n_nb = 0, parent[8];
                for (int k = 0; k < 8; ++k)
                    if (on(r + dr[k], c + dc[k])) nb[n_nb++] = k;
                if (n_nb < 2) continue;
                for (int a = 0; a < n_nb; ++a) parent[a] = a;
                auto find = [&](int a) { while (parent[a] != a) a = parent[a]; return a; };
                for (int a = 0; a < n_nb; ++a)
                    for (int b = a + 1; b < n_nb; ++b)
                        if (std::abs(dr[nb[a]] - dr[nb[b]]) <= 1 &&
                            std::abs(dc[nb[a]] - dc[nb[b]]) <= 1)
                            parent[find(a)] = find(b);
                int n_comp = 0;
                for (int a = 0; a < n_nb; ++a) if (find(a) == a) ++n_comp;
                if (n_comp != 1) continue;
                line[static_cast<std::size_t>(r) * nc + c] = 0.0f;
                changed = true;
            }
    }
}

double median_of(std::vector<double> v)
{
    if (v.empty()) return std::numeric_limits<double>::quiet_NaN();
    std::sort(v.begin(), v.end());
    const std::size_t h = v.size() / 2;
    return (v.size() % 2) ? v[h] : 0.5 * (v[h - 1] + v[h]);
}

} // namespace

riser_result find_risers(const lsdtt_raster_cpu& dem,
                         const lsdtt_raster_cpu& mask_in,
                         const riser_request&    request)
{
    const lsdtt_raster_meta& meta = dem.get_meta();
    const std::size_t n = dem.get_data().size();
    if (!mask_in.get_data().empty() && mask_in.get_data().size() != n)
        throw std::invalid_argument("find_risers: mask and DEM differ in size");
    if (request.gap_closing == riser_gap_closing_choice::astar &&
        !(request.skeleton.elevation_connect.search_radius_metres > 0.0))
        throw std::invalid_argument(
            "find_risers: gap closing astar needs a search radius > 0");

    // Every valid cell when no mask is given.
    const lsdtt_raster_cpu mask = mask_in.get_data().empty()
        ? lsdtt_raster_cpu(std::vector<float>(n, 1.0f),
                           std::vector<uint8_t>(dem.get_mask()), meta)
        : mask_in;

    riser_result out;
    vec_skeleton_request vsr = request.skeleton;
    vsr.n_threads = request.n_threads;
    vsr.elevation_connect.search_radius_metres = 0.0;   // connected below
    out.skeleton = find_vec_skeleton(dem, mask, vsr);

    // ---- gap closing -----------------------------------------------------
    // The skeleton is valid only inside the search mask (the VEC is masked to
    // it). The bridges are routed on a copy valid wherever the DEM is: a road
    // or a gully cut out of the mask is exactly the gap a riser needs closing
    // across. Detection stays inside the mask; a bridge is bounded by the
    // search radius and the elevation tolerance.
    lsdtt_raster_cpu line = out.skeleton.skeleton;
    if (request.gap_closing == riser_gap_closing_choice::astar)
    {
        std::vector<float> d(n, 0.0f);
        const std::vector<float>&   sd = out.skeleton.skeleton.get_data();
        const std::vector<uint8_t>& sm = out.skeleton.skeleton.get_mask();
        for (std::size_t i = 0; i < n; ++i)
            if (sm[i] && sd[i] > 0.5f) d[i] = 1.0f;
        const lsdtt_raster_cpu skel_on_dem(std::move(d),
                                           std::vector<uint8_t>(dem.get_mask()), meta);
        elevation_connect_result ec = connect_skeleton_by_elevation(
            skel_on_dem, dem, request.skeleton.elevation_connect);
        out.skeleton.break_bridge_mask  = ec.bridge_mask;
        out.skeleton.skeleton_connected = ec.skeleton_connected;
        line = std::move(ec.skeleton_connected);
    }
    {
        std::vector<float> added(n, 0.0f);
        const std::vector<float>& a = line.get_data();
        const std::vector<float>& b = out.skeleton.skeleton.get_data();
        for (std::size_t i = 0; i < n && i < a.size(); ++i)
            if (a[i] > 0.5f && !(i < b.size() && b[i] > 0.5f)) added[i] = 1.0f;
        out.gap_closing_mask = lsdtt_raster_cpu(
            std::move(added), std::vector<uint8_t>(dem.get_mask()), meta);
    }

    // ---- segments and their step heights ----------------------------------
    // Gap closing can leave the line more than one pixel wide; segments need a
    // one-pixel skeleton with no staircase corners.
    lsdtt_raster_cpu thin = thin_zhang_suen(line);
    {
        std::vector<float> d = thin.get_data();
        remove_staircase_corners(d, meta.n_rows, meta.n_cols);
        thin = lsdtt_raster_cpu(std::move(d), std::vector<uint8_t>(thin.get_mask()), meta);
    }
    const std::vector<skeleton_segment> segs = extract_skeleton_segments(thin);
    out.n_segments_found = static_cast<std::int64_t>(segs.size());

    const double cell = meta.cell_size;
    // The skeleton follows the curvature maximum, which the polyfit window
    // smears off the step itself (25 m off a sharp step at lambda 100 m), so
    // the default window straddles half the polyfit window.
    const double d_m  = (request.step_window_radius_metres > 0.0)
        ? request.step_window_radius_metres
        : std::max(3.0 * cell, 0.5 * request.skeleton.polyfit_lambda_metres);
    const int    hw   = std::max(1, request.tangent_half_window_pixels);

    std::vector<float> keep(n, 0.0f);
    std::int64_t next_id = 1;
    for (const skeleton_segment& sg : segs)
    {
        const int np = static_cast<int>(sg.pixels.size());
        if (np < std::max(1, request.min_segment_pixels))
        {
            ++out.n_dropped_short;
            continue;
        }
        std::vector<double> steps, grads, elev;
        double length = 0.0;
        for (int k = 0; k < np; ++k)
        {
            const skeleton_pixel& p = sg.pixels[static_cast<std::size_t>(k)];
            if (k > 0)
            {
                const skeleton_pixel& q = sg.pixels[static_cast<std::size_t>(k - 1)];
                length += std::hypot(p.row - q.row, p.col - q.col) * cell;
            }
            const double z0 = sample(dem, p.row, p.col);
            if (std::isfinite(z0)) elev.push_back(z0);

            // Local direction from the pixels either side along the segment.
            const int a = std::max(0, k - hw);
            const int b = std::min(np - 1, k + hw);
            const double tr = sg.pixels[static_cast<std::size_t>(b)].row -
                              sg.pixels[static_cast<std::size_t>(a)].row;
            const double tc = sg.pixels[static_cast<std::size_t>(b)].col -
                              sg.pixels[static_cast<std::size_t>(a)].col;
            const double tl = std::hypot(tr, tc);
            if (!(tl > 0.0)) continue;
            const double nr = -tc / tl, nc = tr / tl;   // unit normal, in cells
            const double ur = tr / tl, uc = tc / tl;   // unit tangent
            // Mean elevation over a short run along the tangent at offset
            // s cells along the normal: one-pixel samples let 0.3 m noise
            // pass a 1 m step floor.
            auto z_at = [&](double s)
            {
                double sum = 0.0;
                int    cnt = 0;
                for (int t = -hw; t <= hw; ++t)
                {
                    const double z = sample(dem,
                        static_cast<int>(std::lround(p.row + s * nr + t * ur)),
                        static_cast<int>(std::lround(p.col + s * nc + t * uc)));
                    if (std::isfinite(z)) { sum += z; ++cnt; }
                }
                return cnt ? sum / cnt : std::numeric_limits<double>::quiet_NaN();
            };
            auto dz = [&](double r_m)
            {
                const double s = r_m / cell;
                return z_at(s) - z_at(-s);
            };
            const double s1 = dz(d_m), s2 = dz(2.0 * d_m);
            if (std::isfinite(s1) && std::isfinite(s2))
                steps.push_back(std::fabs(2.0 * s1 - s2));
            // Steepest centred gradient along the normal within +-d: the
            // skeleton sits off the step, so the face is somewhere in there.
            double g_max = std::numeric_limits<double>::quiet_NaN();
            const int n_off = static_cast<int>(std::lround(d_m / cell));
            for (int o = -n_off; o <= n_off; ++o)
            {
                const double g = std::fabs(z_at(o + 1.0) - z_at(o - 1.0)) / (2.0 * cell);
                if (std::isfinite(g) && !(g <= g_max)) g_max = g;
            }
            if (std::isfinite(g_max)) grads.push_back(g_max);
        }
        const double step = median_of(steps);
        if (!(step >= request.min_step_height_metres))
        {
            ++out.n_dropped_low_step;
            continue;
        }
        riser_segment rs;
        rs.id            = next_id++;
        rs.n_pixels      = np;
        rs.length_m      = length;
        rs.step_height_m = step;
        rs.face_gradient = median_of(grads);
        rs.is_loop       = sg.is_loop;
        double se = 0.0;
        for (double z : elev) se += z;
        rs.mean_elevation_m = elev.empty() ? 0.0 : se / static_cast<double>(elev.size());
        rs.polyline = skeleton_segment_to_polyline(sg, meta);
        for (const skeleton_pixel& p : sg.pixels)
        {
            const std::size_t i = static_cast<std::size_t>(p.row) * meta.n_cols + p.col;
            rs.flats.push_back(static_cast<std::int64_t>(i));
            keep[i] = 1.0f;
        }
        out.segments.push_back(std::move(rs));
    }
    out.centreline = lsdtt_raster_cpu(std::move(keep),
                                      std::vector<uint8_t>(dem.get_mask()), meta);
    return out;
}

} // namespace lsdtt
