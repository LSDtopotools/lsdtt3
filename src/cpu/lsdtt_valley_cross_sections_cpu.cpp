// lsdtt_valley_cross_sections_cpu.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html
//
// -----------------------------------------------------------------------------
// lsdtt_valley_cross_sections_cpu
//
// Implementation. See the header for the algorithm outline.
// -----------------------------------------------------------------------------

#include "lsdtt_valley_cross_sections_cpu.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

#include "../utils/lsdtt_raster_geometry.hpp"

namespace lsdtt
{
namespace
{

// A point on the centreline polyline plus its cumulative along-distance.
struct poly_node
{
    double x = 0.0, y = 0.0, along = 0.0;
};

struct station_geom
{
    double cx = 0.0, cy = 0.0, nx = 0.0, ny = 0.0;
    std::size_t node = 0;
};

// Sample the DEM at a world coordinate (nearest cell). Returns validity.
bool sample_dem(const lsdtt_raster_meta& meta,
                const std::vector<float>& elev,
                const std::vector<std::uint8_t>& emask,
                double x, double y, float& out)
{
    if (!point_within_raster(meta, x, y)) return false;
    const double fcol = (x - meta.x_min) / meta.cell_size;
    const double frow = (meta.y_max() - y) / meta.cell_size;
    int col = static_cast<int>(std::floor(fcol));
    int row = static_cast<int>(std::floor(frow));
    if (col < 0 || col >= meta.n_cols || row < 0 || row >= meta.n_rows)
        return false;
    const std::size_t i =
        static_cast<std::size_t>(row) * meta.n_cols + col;
    if (i >= emask.size() || !emask[i]) return false;
    out = elev[i];
    return true;
}

// Sample an arbitrary raster at a world coordinate (nearest cell). Returns
// validity (in-bounds and the raster's own mask valid).
bool sample_raster(const lsdtt_raster_meta& meta,
                   const lsdtt_raster_cpu& r,
                   double x, double y, float& out)
{
    if (!point_within_raster(meta, x, y)) return false;
    const double fcol = (x - meta.x_min) / meta.cell_size;
    const double frow = (meta.y_max() - y) / meta.cell_size;
    int col = static_cast<int>(std::floor(fcol));
    int row = static_cast<int>(std::floor(frow));
    if (col < 0 || col >= meta.n_cols || row < 0 || row >= meta.n_rows)
        return false;
    const std::size_t i =
        static_cast<std::size_t>(row) * meta.n_cols + col;
    const std::vector<std::uint8_t>& m = r.get_mask();
    if (i >= m.size() || !m[i]) return false;
    out = r.get_data()[i];
    return true;
}

// Report prominent local maxima of a 1-D profile as sample INDICES. A value is
// a peak if it is a strict local maximum whose prominence (peak minus the
// higher of the two flanking valley minima, walking until a higher peak or an
// end) is at least min_prom_frac of the profile's valid value range. Invalid
// samples break runs. Edge-margin samples are ignored.
std::vector<int> prominent_peaks(
    const std::vector<float>& v, const std::vector<char>& valid,
    double min_prom_frac, int edge_margin)
{
    std::vector<int> out;
    const int n = static_cast<int>(v.size());
    if (n < 3) return out;
    float vmin = 0.0f, vmax = 0.0f; bool any = false;
    for (int i = 0; i < n; ++i)
        if (valid[i])
        {
            if (!any) { vmin = vmax = v[i]; any = true; }
            else { vmin = std::min(vmin, v[i]); vmax = std::max(vmax, v[i]); }
        }
    if (!any) return out;
    const double range = static_cast<double>(vmax - vmin);
    if (range <= 0.0) return out;
    const double min_prom = min_prom_frac * range;
    const int lo = std::max(1, edge_margin);
    const int hi = std::min(n - 2, n - 1 - edge_margin);
    for (int i = lo; i <= hi; ++i)
    {
        if (!valid[i] || !valid[i - 1] || !valid[i + 1]) continue;
        if (!(v[i] > v[i - 1] && v[i] >= v[i + 1])) continue;
        // Left key col: walk left to a higher sample or an invalid/end.
        double left_min = v[i];
        for (int j = i - 1; j >= 0 && valid[j]; --j)
        {
            if (v[j] > v[i]) break;
            left_min = std::min(left_min, static_cast<double>(v[j]));
        }
        double right_min = v[i];
        for (int j = i + 1; j < n && valid[j]; ++j)
        {
            if (v[j] > v[i]) break;
            right_min = std::min(right_min, static_cast<double>(v[j]));
        }
        const double prom = v[i] - std::max(left_min, right_min);
        if (prom >= min_prom) out.push_back(i);
    }
    return out;
}

struct hydraulic_curves
{
    std::vector<double> stage;
    std::vector<double> width;
    double ds = 0.0;
};

// Flood the profile and compute W(h) over the wetted run connected to the
// thalweg (matching the Python chart routine).
hydraulic_curves cross_section_hydraulics(
    const std::vector<double>& offset,
    const std::vector<float>& elevation,
    const std::vector<char>& valid,
    double stage_step)
{
    hydraulic_curves out;
    if (stage_step <= 0.0) return out;

    std::vector<std::pair<double,double>> p;
    p.reserve(offset.size());
    for (std::size_t i = 0; i < offset.size() && i < elevation.size(); ++i)
        if (i < valid.size() && valid[i])
            p.emplace_back(offset[i], static_cast<double>(elevation[i]));
    if (p.size() < 3) return out;
    std::sort(p.begin(), p.end(),
              [](const auto& a, const auto& b){ return a.first < b.first; });

    std::vector<double> off, z;
    off.reserve(p.size());
    z.reserve(p.size());
    for (const auto& q : p) { off.push_back(q.first); z.push_back(q.second); }

    std::vector<double> d;
    d.reserve(off.size() - 1);
    for (std::size_t i = 1; i < off.size(); ++i)
        d.push_back(off[i] - off[i - 1]);
    if (d.empty()) return out;
    std::sort(d.begin(), d.end());
    out.ds = d[d.size() / 2];
    if (!(out.ds > 0.0)) return out;

    const auto it_min = std::min_element(z.begin(), z.end());
    const auto it_max = std::max_element(z.begin(), z.end());
    if (it_min == z.end() || it_max == z.end()) return out;
    const double zmin = *it_min;
    const double zmax = *it_max;
    if (!(zmax > zmin)) return out;

    const int n = static_cast<int>(z.size());
    const int thalweg = static_cast<int>(std::distance(z.begin(), it_min));
    for (double h = zmin + stage_step; h < zmax; h += stage_step)
    {
        if (z[thalweg] >= h) continue;
        int lo = thalweg;
        while (lo - 1 >= 0 && z[lo - 1] < h) --lo;
        int hi = thalweg;
        while (hi + 1 < n && z[hi + 1] < h) ++hi;
        const int nseg = hi - lo + 1;
        if (nseg <= 0) continue;
        const double width = static_cast<double>(nseg) * out.ds;
        out.stage.push_back(h);
        out.width.push_back(width);
    }
    return out;
}

void combine_width_steps(
    const std::vector<double>& stage,
    const std::vector<double>& width,
    double ds,
    double min_step_height_metres,
    int min_step_width_pixels,
    std::vector<double>& boundaries,
    std::vector<double>& mids)
{
    boundaries.clear();
    mids.clear();
    if (stage.size() < 2 || width.size() != stage.size() || !(ds > 0.0))
        return;

    std::vector<std::pair<double,double>> sw;
    sw.reserve(stage.size());
    for (std::size_t i = 0; i < stage.size(); ++i)
        if (std::isfinite(stage[i]) && std::isfinite(width[i]))
            sw.emplace_back(stage[i], width[i]);
    if (sw.size() < 2) return;
    std::sort(sw.begin(), sw.end(),
              [](const auto& a, const auto& b){ return a.first < b.first; });

    std::vector<double> s, w;
    s.reserve(sw.size());
    w.reserve(sw.size());
    for (const auto& q : sw) { s.push_back(q.first); w.push_back(q.second); }

    std::vector<float> dwdh(s.size(), 0.0f);
    if (s.size() >= 3)
    {
        for (std::size_t i = 1; i + 1 < s.size(); ++i)
        {
            const double ds_local = s[i + 1] - s[i - 1];
            dwdh[i] = (std::abs(ds_local) > 0.0)
                ? static_cast<float>((w[i + 1] - w[i - 1]) / ds_local)
                : 0.0f;
        }
        dwdh.front() = dwdh[1];
        dwdh.back() = dwdh[dwdh.size() - 2];
    }
    std::vector<char> valid(dwdh.size(), 1);
    const std::vector<int> pk = (dwdh.size() >= 3)
        ? prominent_peaks(dwdh, valid, 0.10, 1)
        : std::vector<int>{};

    std::vector<int> cand;
    cand.reserve(pk.size() + 2);
    cand.push_back(0);
    for (int i : pk)
        if (i > 0 && i < static_cast<int>(s.size()) - 1)
            cand.push_back(i);
    cand.push_back(static_cast<int>(s.size()) - 1);
    std::sort(cand.begin(), cand.end());
    cand.erase(std::unique(cand.begin(), cand.end()), cand.end());
    if (cand.size() < 2) return;

    const double min_w = static_cast<double>(std::max(0, min_step_width_pixels))
                         * ds;
    std::vector<int> keep;
    keep.push_back(cand.front());
    int ref = cand.front();
    for (std::size_t j = 1; j + 1 < cand.size(); ++j)
    {
        const int i = cand[j];
        if ((w[i] - w[ref]) >= min_w &&
            (s[i] - s[ref]) >= min_step_height_metres)
        {
            keep.push_back(i);
            ref = i;
        }
    }
    if (keep.back() != cand.back()) keep.push_back(cand.back());

    for (int i : keep)
    {
        const double z = s[i];
        if (boundaries.empty() || (z - boundaries.back()) > 1e-9)
            boundaries.push_back(z);
    }
    if (boundaries.size() < 2) { boundaries.clear(); return; }

    mids.reserve(boundaries.size() - 1);
    for (std::size_t k = 0; k + 1 < boundaries.size(); ++k)
        mids.push_back(0.5 * (boundaries[k] + boundaries[k + 1]));
}

std::vector<double> water_line_crossings(
    const std::vector<double>& offset,
    const std::vector<float>& elevation,
    const std::vector<char>& valid,
    double level)
{
    std::vector<std::pair<double,double>> p;
    p.reserve(offset.size());
    for (std::size_t i = 0; i < offset.size() && i < elevation.size(); ++i)
        if (i < valid.size() && valid[i] && std::isfinite(offset[i]))
            p.emplace_back(offset[i], static_cast<double>(elevation[i]));
    if (p.size() < 2) return {};
    std::sort(p.begin(), p.end(),
              [](const auto& a, const auto& b){ return a.first < b.first; });

    std::vector<double> out;
    out.reserve(8);
    for (std::size_t i = 0; i + 1 < p.size(); ++i)
    {
        const double x0 = p[i].first;
        const double z0 = p[i].second;
        const double x1 = p[i + 1].first;
        const double z1 = p[i + 1].second;
        const double a = z0 - level;
        const double b = z1 - level;
        if (a == 0.0) out.push_back(x0);
        if ((a < 0.0 && b > 0.0) || (a > 0.0 && b < 0.0))
        {
            const double t = a / (a - b);
            out.push_back(x0 + t * (x1 - x0));
        }
    }
    if ((p.back().second - level) == 0.0) out.push_back(p.back().first);
    return out;
}

double optimize_step_midpoint(
    const std::vector<double>& offset,
    const std::vector<float>& elevation,
    const std::vector<char>& valid,
    double zmid,
    double zlo,
    double zhi,
    int n_candidates,
    int min_crossings)
{
    if (!std::isfinite(zlo) || !std::isfinite(zhi) || !(zhi > zlo))
        return zmid;

    const double eps = std::max(1e-9, 1e-6 * (zhi - zlo));
    const double lo = zlo + eps;
    const double hi = zhi - eps;
    if (!(hi > lo)) return 0.5 * (zlo + zhi);

    double z0 = std::isfinite(zmid) ? zmid : 0.5 * (lo + hi);
    z0 = std::max(lo, std::min(hi, z0));

    const int n = std::max(5, n_candidates);
    std::vector<double> cands;
    cands.reserve(static_cast<std::size_t>(n) + 1);
    for (int i = 0; i < n; ++i)
    {
        const double t = (n > 1) ? static_cast<double>(i) / (n - 1) : 0.0;
        cands.push_back(lo + t * (hi - lo));
    }
    cands.push_back(z0);
    std::sort(cands.begin(), cands.end());
    cands.erase(std::unique(cands.begin(), cands.end()), cands.end());

    auto crossings_count = [&](double z)
    {
        return static_cast<int>(water_line_crossings(offset, elevation, valid, z)
                                    .size());
    };

    int best_any_i = -1;
    int best_any_n = std::numeric_limits<int>::max();
    double best_any_dz = std::numeric_limits<double>::infinity();

    int best_pref_i = -1;
    int best_pref_n = std::numeric_limits<int>::max();
    double best_pref_dz = std::numeric_limits<double>::infinity();

    for (std::size_t i = 0; i < cands.size(); ++i)
    {
        const int ncr = crossings_count(cands[i]);
        const double dz = std::abs(cands[i] - z0);

        if (ncr < best_any_n || (ncr == best_any_n && dz < best_any_dz))
        {
            best_any_i = static_cast<int>(i);
            best_any_n = ncr;
            best_any_dz = dz;
        }

        if (ncr >= min_crossings)
        {
            if (ncr < best_pref_n ||
                (ncr == best_pref_n && dz < best_pref_dz))
            {
                best_pref_i = static_cast<int>(i);
                best_pref_n = ncr;
                best_pref_dz = dz;
            }
        }
    }

    if (best_pref_i >= 0) return cands[best_pref_i];
    return (best_any_i >= 0) ? cands[best_any_i] : z0;
}

std::vector<valley_cross_section_hydraulic_edge> hydraulic_edges_from_samples(
    const valley_cross_section& xs,
    const valley_cross_section_request& request)
{
    std::vector<valley_cross_section_hydraulic_edge> edges;
    if (xs.samples.size() < 3) return edges;

    std::vector<double> off;
    std::vector<float> elev;
    std::vector<char> valid;
    off.reserve(xs.samples.size());
    elev.reserve(xs.samples.size());
    valid.reserve(xs.samples.size());
    for (const valley_cross_section_sample& sp : xs.samples)
    {
        off.push_back(sp.offset_metres);
        elev.push_back(sp.elevation);
        valid.push_back(sp.elevation_valid ? 1 : 0);
    }

    hydraulic_curves hyd = cross_section_hydraulics(
        off, elev, valid, request.hydraulic_stage_step_metres);
    if (hyd.stage.empty()) return edges;

    std::vector<double> boundaries, mids;
    combine_width_steps(
        hyd.stage,
        hyd.width,
        hyd.ds,
        request.hydraulic_min_step_height_metres,
        request.hydraulic_min_step_width_pixels,
        boundaries,
        mids);
    if (boundaries.size() < 2 || mids.empty()) return edges;

    for (std::size_t k = 0; k + 1 < boundaries.size(); ++k)
    {
        const double zmid = (k < mids.size()) ? mids[k]
                                              : 0.5 * (boundaries[k] + boundaries[k + 1]);
        const double zopt = optimize_step_midpoint(
            off, elev, valid,
            zmid,
            boundaries[k], boundaries[k + 1],
            request.hydraulic_midpoint_search_samples,
            request.hydraulic_min_crossings_preferred);

        const std::vector<double> xoff = water_line_crossings(off, elev, valid, zopt);
        for (double o : xoff)
        {
            valley_cross_section_hydraulic_edge e;
            e.step = static_cast<std::int32_t>(k);
            e.stage_metres = zopt;
            e.offset_metres = o;
            edges.push_back(e);
        }
    }
    return edges;
}

// Is (x,y) inside the valley mask (per-channel raster, data 1 = valley)?
bool in_valley_at(const lsdtt_raster_meta& meta,
                  const std::vector<float>& vmask,
                  double x, double y)
{
    if (!point_within_raster(meta, x, y)) return false;
    const double fcol = (x - meta.x_min) / meta.cell_size;
    const double frow = (meta.y_max() - y) / meta.cell_size;
    int col = static_cast<int>(std::floor(fcol));
    int row = static_cast<int>(std::floor(frow));
    if (col < 0 || col >= meta.n_cols || row < 0 || row >= meta.n_rows)
        return false;
    const std::size_t i =
        static_cast<std::size_t>(row) * meta.n_cols + col;
    return i < vmask.size() && vmask[i] > 0.5f;
}

// March outward from the centre along +dir*normal (dir = +1 right, -1 left)
// and return the offset (>= 0) at which the contiguous in-valley run ends.
// Small holes up to gap_tol consecutive out-steps are bridged. If the march
// reaches march_limit still in the valley, found = false and the limit is
// returned.
double find_valley_side(const lsdtt_raster_meta& meta,
                        const std::vector<float>& vmask,
                        double cx, double cy,
                        double nx, double ny, int dir,
                        double step, double march_limit,
                        int gap_tol, bool& found)
{
    double last_in = 0.0;   // last offset confirmed in-valley
    int consecutive_out = 0;
    found = false;
    for (double off = step; off <= march_limit + 0.5 * step; off += step)
    {
        const double x = cx + dir * off * nx;
        const double y = cy + dir * off * ny;
        if (in_valley_at(meta, vmask, x, y))
        {
            last_in = off;
            consecutive_out = 0;
        }
        else
        {
            ++consecutive_out;
            if (consecutive_out > gap_tol)
            {
                // The side is the last in-valley offset (the first out step
                // of this run). This is the contiguous run containing centre.
                found = true;
                return last_in;
            }
        }
    }
    // Reached the limit still (mostly) inside the valley.
    return march_limit;
}

double neighbour_clip_limit(
    const std::vector<station_geom>& geom,
    std::size_t s,
    int dir,
    int k_window,
    double base_limit)
{
    const station_geom& gs = geom[s];
    double limit = base_limit;
    const int n = static_cast<int>(geom.size());
    const int is = static_cast<int>(s);
    const int k = std::max(1, k_window);
    for (int dt = -k; dt <= k; ++dt)
    {
        if (dt == 0) continue;
        const int t = is + dt;
        if (t < 0 || t >= n) continue;
        const station_geom& gt = geom[static_cast<std::size_t>(t)];
        const double dx = gt.cx - gs.cx;
        const double dy = gt.cy - gs.cy;
        const double d2 = dx * dx + dy * dy;
        if (!(d2 > 0.0)) continue;
        const double proj = gs.nx * dx + gs.ny * dy;
        const double denom = 2.0 * static_cast<double>(dir) * proj;
        if (denom <= 1e-9) continue;
        const double oeq = d2 / denom;
        if (oeq > 0.0) limit = std::min(limit, oeq);
    }
    return limit;
}

} // namespace

valley_cross_section_result extract_valley_cross_sections(
    const lsdtt_raster_cpu&             dem,
    const valley_centreline_result&     centrelines,
    const valley_cross_section_request& request,
    const lsdtt_raster_cpu*             curvature,
    const lsdtt_raster_cpu*             hillshade,
    const lsdtt_raster_cpu*             channel_mask)
{
    const lsdtt_raster_meta& meta = dem.get_meta();
    const std::vector<float>&        elev  = dem.get_data();
    const std::vector<std::uint8_t>& emask = dem.get_mask();
    const int    n_cols = meta.n_cols;
    const double cell    = meta.cell_size;

    const double step = (request.sample_step_metres > 0.0)
        ? request.sample_step_metres : cell;
    const double spacing = std::max(request.spacing_metres, cell);
    const int    win = std::max(1, request.tangent_window_nodes);
    const int    gap_tol = std::max(0, request.side_gap_tolerance_steps);

    valley_cross_section_result result;

    for (const valley_channel_centreline& ch : centrelines.channels)
    {
        if (ch.centreline_flats.size() < 2) continue;
        const std::vector<float>& vmask = ch.valley_mask.get_data();
        if (vmask.empty()) continue;
        if (vmask.size() != elev.size())
            throw std::invalid_argument(
                "extract_valley_cross_sections: valley mask / DEM size "
                "mismatch.");

        // Build the centreline polyline with cumulative along-distance.
        const std::size_t nn = ch.centreline_flats.size();
        std::vector<poly_node> poly(nn);
        for (std::size_t i = 0; i < nn; ++i)
        {
            const std::int64_t f = ch.centreline_flats[i];
            const int r = static_cast<int>(f / n_cols);
            const int c = static_cast<int>(f % n_cols);
            const lsdtt_xy xy = cell_center_from_row_col(meta, r, c);
            poly[i].x = xy.x;
            poly[i].y = xy.y;
            if (i > 0)
                poly[i].along = poly[i - 1].along +
                    std::hypot(poly[i].x - poly[i - 1].x,
                               poly[i].y - poly[i - 1].y);
        }
        const double total_along = poly.back().along;

        // Station along-distances (regular spacing, offset half a spacing in
        // so the head/outlet stubs are not clipped).
        std::vector<double> station_along;
        for (double d = 0.5 * spacing; d <= total_along; d += spacing)
            station_along.push_back(d);
        if (station_along.empty())
            station_along.push_back(0.5 * total_along);

        // For each station: interpolate the centre and compute the smoothed
        // local tangent/normal. Also record the nearest node index.
        std::vector<station_geom> geom(station_along.size());

        auto node_at_along = [&](double d) -> std::size_t
        {
            // First node whose along >= d (binary search).
            std::size_t lo = 0, hi = nn - 1;
            while (lo < hi)
            {
                const std::size_t mid = (lo + hi) / 2;
                if (poly[mid].along < d) lo = mid + 1; else hi = mid;
            }
            return lo;
        };

        int station_seq_out = 0;
        for (std::size_t s = 0; s < station_along.size(); ++s)
        {
            const double d = station_along[s];
            const std::size_t k = node_at_along(d);
            // Interpolate the centre between k-1 and k.
            double cx = poly[k].x, cy = poly[k].y;
            if (k > 0)
            {
                const double d0 = poly[k - 1].along, d1 = poly[k].along;
                const double t = (d1 > d0) ? (d - d0) / (d1 - d0) : 0.0;
                cx = poly[k - 1].x + t * (poly[k].x - poly[k - 1].x);
                cy = poly[k - 1].y + t * (poly[k].y - poly[k - 1].y);
            }
            // Smoothed tangent from a node window around k.
            const std::size_t ka =
                (k > static_cast<std::size_t>(win)) ? k - win : 0;
            const std::size_t kb =
                std::min(nn - 1, k + static_cast<std::size_t>(win));
            double tx = poly[kb].x - poly[ka].x;
            double ty = poly[kb].y - poly[ka].y;
            const double tlen = std::hypot(tx, ty);
            if (tlen > 0.0) { tx /= tlen; ty /= tlen; }
            else            { tx = 1.0; ty = 0.0; }
            // Right-of-downstream normal in world coords (x east, y north):
            // rotate tangent by -90 deg -> (ty, -tx).
            geom[s] = { cx, cy, ty, -tx, k };
        }

        // Pass 1: valley width at every station (march both sides), capped at
        // max_half_length so an unbounded valley (mask connected to a
        // reservoir / coastal flat / the sea) cannot run away. Only stations
        // where BOTH sides are found within the cap — a genuinely bounded
        // valley — contribute to the per-channel width used to fix the
        // transect length; a reservoir-crossing station reaches the cap
        // without finding its far side and is excluded.
        const double march_limit =
            std::max(request.max_half_length_metres, 4.0 * cell);
        std::vector<double> left_off(station_along.size(), 0.0);
        std::vector<double> right_off(station_along.size(), 0.0);
        std::vector<char>   left_ok(station_along.size(), 0);
        std::vector<char>   right_ok(station_along.size(), 0);
        double max_width = 0.0;
        for (std::size_t s = 0; s < station_along.size(); ++s)
        {
            bool lf = false, rf = false;
            const double lo = find_valley_side(
                meta, vmask, geom[s].cx, geom[s].cy,
                geom[s].nx, geom[s].ny, -1, step, march_limit, gap_tol, lf);
            const double ro = find_valley_side(
                meta, vmask, geom[s].cx, geom[s].cy,
                geom[s].nx, geom[s].ny, +1, step, march_limit, gap_tol, rf);
            left_off[s] = lo;  left_ok[s] = lf ? 1 : 0;
            right_off[s] = ro; right_ok[s] = rf ? 1 : 0;
            if (lf && rf) max_width = std::max(max_width, lo + ro);
        }
        // No bounded station (fully unbounded valley): fall back to the cap.
        if (max_width <= 0.0) max_width = 2.0 * march_limit;

        double half_len = 0.5 * request.half_length_multiplier * max_width;
        half_len = std::max(half_len, request.min_half_length_metres);
        half_len = std::min(half_len, march_limit);

        result.channel_source_ids.push_back(ch.source_id);
        result.channel_half_length_metres.push_back(half_len);
        result.channel_max_valley_width_metres.push_back(max_width);

        // Pass 2: sample each transect over clipped limits (Phase-1 non-overlap
        // clips concave-side rays at neighbour bisectors).
        std::vector<double> left_limit(station_along.size(), half_len);
        std::vector<double> right_limit(station_along.size(), half_len);
        std::vector<char> left_trunc(station_along.size(), 0);
        std::vector<char> right_trunc(station_along.size(), 0);
        if (request.clip_overlapping_sections)
        {
            for (std::size_t s = 0; s < station_along.size(); ++s)
            {
                const double l = neighbour_clip_limit(
                    geom, s, -1, request.neighbour_clip_window, half_len);
                const double r = neighbour_clip_limit(
                    geom, s, +1, request.neighbour_clip_window, half_len);
                left_limit[s] = std::max(step, std::min(half_len, l));
                right_limit[s] = std::max(step, std::min(half_len, r));
                left_trunc[s] = (left_limit[s] < (half_len - 1e-6)) ? 1 : 0;
                right_trunc[s] = (right_limit[s] < (half_len - 1e-6)) ? 1 : 0;
            }
        }

        for (std::size_t s = 0; s < station_along.size(); ++s)
        {
            const double full_len = 2.0 * half_len;
            const double clipped_len = left_limit[s] + right_limit[s];
            if (request.clip_overlapping_sections && full_len > 0.0)
            {
                const double frac = clipped_len / full_len;
                if (frac < request.min_section_reach_fraction) continue;
            }

            valley_cross_section xs;
            xs.source_id   = ch.source_id;
            xs.station_seq = static_cast<std::int32_t>(station_seq_out++);
            xs.along_distance_metres = station_along[s];
            xs.centre_x = geom[s].cx;
            xs.centre_y = geom[s].cy;
            xs.normal_dx = geom[s].nx;
            xs.normal_dy = geom[s].ny;
            xs.left_truncated_by_neighbour = left_trunc[s] ? true : false;
            xs.right_truncated_by_neighbour = right_trunc[s] ? true : false;

            const double march_left = left_limit[s];
            const double march_right = right_limit[s];

            // Re-evaluate valley sides against this station's clip limits.
            bool lf = false, rf = false;
            const double lo = find_valley_side(
                meta, vmask, geom[s].cx, geom[s].cy,
                geom[s].nx, geom[s].ny, -1, step, march_left, gap_tol, lf);
            const double ro = find_valley_side(
                meta, vmask, geom[s].cx, geom[s].cy,
                geom[s].nx, geom[s].ny, +1, step, march_right, gap_tol, rf);
            xs.left_found  = lf;
            xs.right_found = rf;
            xs.left_offset_metres  = -lo;
            xs.right_offset_metres =  ro;
            xs.left_x  = geom[s].cx - lo * geom[s].nx;
            xs.left_y  = geom[s].cy - lo * geom[s].ny;
            xs.right_x = geom[s].cx + ro * geom[s].nx;
            xs.right_y = geom[s].cy + ro * geom[s].ny;
            xs.valley_width_metres = lo + ro;

            const int n_left = static_cast<int>(
                std::floor(march_left / step + 0.5));
            const int n_right = static_cast<int>(
                std::floor(march_right / step + 0.5));
            xs.samples.reserve(n_left + n_right + 1);
            for (int j = -n_left; j <= n_right; ++j)
            {
                const double off = j * step;
                valley_cross_section_sample sp;
                sp.offset_metres = off;
                sp.x = geom[s].cx + off * geom[s].nx;
                sp.y = geom[s].cy + off * geom[s].ny;
                float e = 0.0f;
                sp.elevation_valid =
                    sample_dem(meta, elev, emask, sp.x, sp.y, e) ? 1 : 0;
                sp.elevation = e;
                sp.in_valley =
                    in_valley_at(meta, vmask, sp.x, sp.y) ? 1 : 0;
                float cv = 0.0f;
                if (curvature &&
                    sample_raster(meta, *curvature, sp.x, sp.y, cv))
                { sp.curvature = cv; sp.curvature_valid = 1; }
                float hs = 0.0f;
                if (hillshade &&
                    sample_raster(meta, *hillshade, sp.x, sp.y, hs))
                { sp.hillshade = hs; sp.hillshade_valid = 1; }
                float ch_v = 0.0f;
                if (channel_mask &&
                    sample_raster(meta, *channel_mask, sp.x, sp.y, ch_v))
                    sp.is_channel = (ch_v > 0.5f) ? 1 : 0;
                xs.samples.push_back(sp);
            }

            // --- Per-section feature detection ------------------------------
            const int ns = static_cast<int>(xs.samples.size());
            // (a) D8 channel crossings: the centre offset of each contiguous
            //     run of channel pixels.
            if (channel_mask)
            {
                int run_lo = -1;
                for (int j = 0; j <= ns; ++j)
                {
                    const bool ch = (j < ns) && xs.samples[j].is_channel;
                    if (ch && run_lo < 0) run_lo = j;
                    else if (!ch && run_lo >= 0)
                    {
                        const int mid = (run_lo + j - 1) / 2;
                        xs.channel_crossing_offsets.push_back(
                            xs.samples[mid].offset_metres);
                        run_lo = -1;
                    }
                }
            }
            // (c)/(d) prominent local maxima of curvature and hillshade.
            auto peak_offsets = [&](bool use_curv)
            {
                std::vector<float> vals(ns, 0.0f);
                std::vector<char>  ok(ns, 0);
                for (int j = 0; j < ns; ++j)
                {
                    const valley_cross_section_sample& sp = xs.samples[j];
                    vals[j] = use_curv ? sp.curvature : sp.hillshade;
                    ok[j]   = use_curv ? sp.curvature_valid : sp.hillshade_valid;
                }
                std::vector<double> offs;
                for (int idx : prominent_peaks(
                         vals, ok,
                         request.feature_peak_min_prominence_fraction,
                         request.feature_peak_edge_margin_samples))
                    offs.push_back(xs.samples[idx].offset_metres);
                return offs;
            };
            if (curvature) xs.curvature_peak_offsets = peak_offsets(true);
            if (hillshade) xs.hillshade_peak_offsets = peak_offsets(false);

            // (e) Hydraulic step-midpoint boundaries (terrace/floodplain
            // edges) from the sampled elevation profile.
            xs.hydraulic_edges = hydraulic_edges_from_samples(xs, request);

            result.sections.push_back(std::move(xs));
        }
    }

    return result;
}

} // namespace lsdtt
