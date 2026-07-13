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
#include "lsdtt_proximity_cpu.hpp"

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
    std::vector<double> stage;      // absolute water-surface elevation h
    std::vector<double> width;      // total flow width W
    std::vector<double> area;       // wetted cross-sectional area A
    std::vector<double> perimeter;  // wetted perimeter P
    std::vector<double> hyd_depth;  // hydraulic depth A / W
    double zmin = 0.0;              // thalweg elevation
    double ds = 0.0;
};

// Flood the profile and compute the stage-hydraulic-geometry curves over the
// wetted run connected to the thalweg (matching the Python chart routine).
//
// The HYDRAULIC DEPTH is D_h = A / W (wetted area over flow width). Its purpose
// is to mimic BankfullMapper (Delchiaro et al. 2025, Computers & Geosciences,
// doi:10.1016/j.cageo.2025.106001): that tool detects bankfull / inset
// flat-surface levels as breakpoints in the hydraulic-depth function -- the
// area-to-width ratio plotted against elevation above the thalweg. At a flat
// surface (floodplain / terrace tread) width W increases rapidly while D_h
// stays roughly constant (the flow spreads laterally without deepening), which
// is exactly the breakpoint BankfullMapper picks out.
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
    out.zmin = zmin;

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
        // Wetted area A = integral of flow depth over the run; wetted perimeter
        // P = bed length across the run; hydraulic depth D_h = A / W (the
        // BankfullMapper area-to-width metric, see the function comment above).
        double area = 0.0;
        for (int k = lo; k <= hi; ++k) area += (h - z[k]);
        area *= out.ds;
        // Bed length across the run: one hypot segment per adjacent bed-point
        // pair; a single wetted cell has perimeter = ds.
        double perim = (nseg > 1) ? 0.0 : out.ds;
        for (int k = lo + 1; k <= hi; ++k)
            perim += std::hypot(out.ds, z[k] - z[k - 1]);
        const double hyd_depth = (width > 0.0) ? area / width : 0.0;
        out.stage.push_back(h);
        out.width.push_back(width);
        out.area.push_back(area);
        out.perimeter.push_back(perim);
        out.hyd_depth.push_back(hyd_depth);
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
    const valley_cross_section_request& request,
    std::vector<valley_cross_section_hydraulic_sample>* out_curve = nullptr,
    std::vector<double>* out_boundaries = nullptr)
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

    // (f) hand back the full BankfullMapper hydraulic-depth curve if asked.
    if (out_curve)
    {
        out_curve->clear();
        out_curve->reserve(hyd.stage.size());
        for (std::size_t i = 0; i < hyd.stage.size(); ++i)
        {
            valley_cross_section_hydraulic_sample hs;
            hs.stage_metres = hyd.stage[i];
            hs.depth_above_thalweg_metres = hyd.stage[i] - hyd.zmin;
            hs.width_metres = hyd.width[i];
            hs.area_metres2 = hyd.area[i];
            hs.wetted_perimeter_metres = hyd.perimeter[i];
            hs.hydraulic_depth_metres = hyd.hyd_depth[i];
            out_curve->push_back(hs);
        }
    }

    std::vector<double> boundaries, mids;
    combine_width_steps(
        hyd.stage,
        hyd.width,
        hyd.ds,
        request.hydraulic_min_step_height_metres,
        request.hydraulic_min_step_width_pixels,
        boundaries,
        mids);
    // Hand back the flat-surface separator stages (thalweg, break1, break2, ...)
    // regardless of whether any step midpoints were found.
    if (out_boundaries) *out_boundaries = boundaries;
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

// =============================================================================
// Curved (boundary-fitted) cross sections
//
// Each section is traced as an isocontour of a smooth along-valley (streamwise)
// coordinate field D, so it is straight in straight reaches and curves through
// bends like the constant-streamwise grid lines of a boundary-fitted mesh. D is
// the swath along-distance field, Laplacian-smoothed within the valley corridor
// with the centreline pixels held fixed (Dirichlet anchor) to remove the
// projection creases that would otherwise make sections splay at bends.
// =============================================================================

// Bilinear sample of a raster at world (x,y). All four surrounding cell centres
// must be valid in `mask`, else returns false.
bool bilinear_raster(const lsdtt_raster_meta& meta,
                     const std::vector<float>& data,
                     const std::vector<std::uint8_t>& mask,
                     double x, double y, double& out)
{
    const double cs = meta.cell_size;
    const double col_f = (x - meta.x_min) / cs - 0.5;
    const double row_f = (meta.y_max() - y) / cs - 0.5;
    const int c0 = static_cast<int>(std::floor(col_f));
    const int r0 = static_cast<int>(std::floor(row_f));
    if (c0 < 0 || c0 + 1 >= meta.n_cols || r0 < 0 || r0 + 1 >= meta.n_rows)
        return false;
    const double tc = col_f - c0, tr = row_f - r0;
    const int rr[4] = { r0, r0, r0 + 1, r0 + 1 };
    const int cc[4] = { c0, c0 + 1, c0, c0 + 1 };
    double v[4];
    for (int k = 0; k < 4; ++k)
    {
        const std::size_t i =
            static_cast<std::size_t>(rr[k]) * meta.n_cols + cc[k];
        if (i >= mask.size() || !mask[i]) return false;
        v[k] = data[i];
    }
    const double top = v[0] * (1.0 - tc) + v[1] * tc;
    const double bot = v[2] * (1.0 - tc) + v[3] * tc;
    out = top * (1.0 - tr) + bot * tr;
    return true;
}

// Gradient of the streamwise field D via central differences (bilinear).
bool grad_D(const lsdtt_raster_meta& meta,
            const std::vector<float>& D,
            const std::vector<std::uint8_t>& Dm,
            double x, double y, double& gx, double& gy)
{
    const double h = meta.cell_size;
    double xp, xm, yp, ym;
    if (!bilinear_raster(meta, D, Dm, x + h, y, xp)) return false;
    if (!bilinear_raster(meta, D, Dm, x - h, y, xm)) return false;
    if (!bilinear_raster(meta, D, Dm, x, y + h, yp)) return false;
    if (!bilinear_raster(meta, D, Dm, x, y - h, ym)) return false;
    gx = (xp - xm) / (2.0 * h);
    gy = (yp - ym) / (2.0 * h);
    return true;
}

// One traced side of a curved section: the isocontour points, their cumulative
// arc-length from the centre, and the straight-normal distance to the valley
// wall on this side (used as the shorten-to-wall floor).
struct curved_side
{
    std::vector<lsdtt_xy> pts;
    std::vector<double>   arc;
    double wall = 0.0;
    bool   wall_found = false;
};

// Trace the D-isocontour through (cx,cy) outward on one side, up to arc-length
// max_arc. `nx,ny` is the straight right-of-downstream normal and `sign` (+1
// right, -1 left) picks which way to leave the centre. Marching follows the
// contour tangent (perpendicular to grad D) with a short Newton correction back
// onto the level set each step.
curved_side trace_iso_side(const lsdtt_raster_meta& meta,
                           const std::vector<float>& D,
                           const std::vector<std::uint8_t>& Dm,
                           double cx, double cy, double level,
                           double nx, double ny, int sign,
                           double stp, double max_arc)
{
    curved_side sd;
    sd.pts.push_back({ cx, cy });
    sd.arc.push_back(0.0);
    double px = cx, py = cy, arc = 0.0;
    while (arc < max_arc)
    {
        double gx, gy;
        if (!grad_D(meta, D, Dm, px, py, gx, gy)) break;
        const double gmag = std::hypot(gx, gy);
        if (gmag < 1e-9) break;
        double tx = -gy / gmag, ty = gx / gmag;   // along the contour
        if ((tx * nx + ty * ny) * sign < 0.0) { tx = -tx; ty = -ty; }
        double qx = px + tx * stp;
        double qy = py + ty * stp;
        bool bad = false;
        for (int it = 0; it < 3; ++it)             // project back onto level set
        {
            double d, ggx, ggy;
            if (!bilinear_raster(meta, D, Dm, qx, qy, d)) { bad = true; break; }
            if (!grad_D(meta, D, Dm, qx, qy, ggx, ggy)) break;
            const double gg = ggx * ggx + ggy * ggy;
            if (gg < 1e-12) break;
            const double f = (d - level) / gg;
            qx -= f * ggx; qy -= f * ggy;
        }
        if (bad || !std::isfinite(qx) || !std::isfinite(qy)) break;
        double dchk;
        if (!bilinear_raster(meta, D, Dm, qx, qy, dchk)) break;
        const double seg = std::hypot(qx - px, qy - py);
        if (seg < 0.25 * stp) break;               // stalled
        arc += seg;
        px = qx; py = qy;
        sd.pts.push_back({ px, py });
        sd.arc.push_back(arc);
    }
    return sd;
}

// First index at which any neighbour centre is strictly closer than own centre
// (the Voronoi bisector crossing). Returns pts.size() if never crossed.
std::size_t bisector_cross(const std::vector<lsdtt_xy>& pts,
                           const lsdtt_xy& own,
                           const std::vector<lsdtt_xy>& nbs)
{
    for (std::size_t i = 0; i < pts.size(); ++i)
    {
        const double dself = std::hypot(pts[i].x - own.x, pts[i].y - own.y);
        double dmin = std::numeric_limits<double>::max();
        for (const lsdtt_xy& q : nbs)
            dmin = std::min(dmin, std::hypot(pts[i].x - q.x, pts[i].y - q.y));
        if (dmin < dself - 1e-9) return i;
    }
    return pts.size();
}

// Do the two line segments p1p2 and p3p4 properly intersect?
bool segments_intersect(const lsdtt_xy& p1, const lsdtt_xy& p2,
                        const lsdtt_xy& p3, const lsdtt_xy& p4)
{
    auto cross = [](double ax, double ay, double bx, double by)
    { return ax * by - ay * bx; };
    const double d1 = cross(p4.x - p3.x, p4.y - p3.y, p1.x - p3.x, p1.y - p3.y);
    const double d2 = cross(p4.x - p3.x, p4.y - p3.y, p2.x - p3.x, p2.y - p3.y);
    const double d3 = cross(p2.x - p1.x, p2.y - p1.y, p3.x - p1.x, p3.y - p1.y);
    const double d4 = cross(p2.x - p1.x, p2.y - p1.y, p4.x - p1.x, p4.y - p1.y);
    return ((d1 > 0) != (d2 > 0)) && ((d3 > 0) != (d4 > 0));
}

// Do two curved-section polylines cross? (subsampled for speed.)
bool sections_cross(const std::vector<lsdtt_xy>& a,
                    const std::vector<lsdtt_xy>& b)
{
    for (std::size_t i = 1; i < a.size(); ++i)
        for (std::size_t j = 1; j < b.size(); ++j)
            if (segments_intersect(a[i - 1], a[i], b[j - 1], b[j]))
                return true;
    return false;
}

valley_cross_section_result extract_curved_impl(
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
    const int    n_rows = meta.n_rows;
    const double cell   = meta.cell_size;

    const double step = (request.sample_step_metres > 0.0)
        ? request.sample_step_metres : cell;
    const double spacing = std::max(request.spacing_metres, cell);
    const int    win = std::max(1, request.tangent_window_nodes);
    const int    gap_tol = std::max(0, request.side_gap_tolerance_steps);
    const double clear = std::max(0.0, request.curved_wall_clearance_fraction);
    const int    smooth_passes =
        std::max(0, request.curved_streamwise_smoothing_passes);

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

        // Centreline polyline (with along-distance) and axis-point list.
        const std::size_t nn = ch.centreline_flats.size();
        std::vector<poly_node>  poly(nn);
        std::vector<lsdtt_xy>   axis(nn);
        for (std::size_t i = 0; i < nn; ++i)
        {
            const std::int64_t f = ch.centreline_flats[i];
            const int r = static_cast<int>(f / n_cols);
            const int c = static_cast<int>(f % n_cols);
            const lsdtt_xy xy = cell_center_from_row_col(meta, r, c);
            poly[i].x = xy.x; poly[i].y = xy.y;
            axis[i] = xy;
            if (i > 0)
                poly[i].along = poly[i - 1].along +
                    std::hypot(poly[i].x - poly[i - 1].x,
                               poly[i].y - poly[i - 1].y);
        }
        const double total_along = poly.back().along;

        std::vector<double> station_along;
        for (double d = 0.5 * spacing; d <= total_along; d += spacing)
            station_along.push_back(d);
        if (station_along.empty())
            station_along.push_back(0.5 * total_along);

        auto node_at_along = [&](double d) -> std::size_t
        {
            std::size_t lo = 0, hi = nn - 1;
            while (lo < hi)
            {
                const std::size_t mid = (lo + hi) / 2;
                if (poly[mid].along < d) lo = mid + 1; else hi = mid;
            }
            return lo;
        };

        std::vector<station_geom> geom(station_along.size());
        for (std::size_t s = 0; s < station_along.size(); ++s)
        {
            const double d = station_along[s];
            const std::size_t k = node_at_along(d);
            double cx = poly[k].x, cy = poly[k].y;
            if (k > 0)
            {
                const double d0 = poly[k - 1].along, d1 = poly[k].along;
                const double t = (d1 > d0) ? (d - d0) / (d1 - d0) : 0.0;
                cx = poly[k - 1].x + t * (poly[k].x - poly[k - 1].x);
                cy = poly[k - 1].y + t * (poly[k].y - poly[k - 1].y);
            }
            const std::size_t ka =
                (k > static_cast<std::size_t>(win)) ? k - win : 0;
            const std::size_t kb =
                std::min(nn - 1, k + static_cast<std::size_t>(win));
            double tx = poly[kb].x - poly[ka].x;
            double ty = poly[kb].y - poly[ka].y;
            const double tlen = std::hypot(tx, ty);
            if (tlen > 0.0) { tx /= tlen; ty /= tlen; }
            else            { tx = 1.0; ty = 0.0; }
            geom[s] = { cx, cy, ty, -tx, k };
        }

        // --- Streamwise coordinate field D over the valley corridor ----------
        const swath_distances swath = compute_swath_distances(dem, axis, 0);
        const std::vector<float>&        along_src = swath.along_axis.get_data();
        const std::vector<std::uint8_t>& along_msk = swath.along_axis.get_mask();
        const std::vector<float>&        perp_src  = swath.perp_distance.get_data();
        const std::vector<std::uint8_t>& perp_msk  = swath.perp_distance.get_mask();

        // Fixed transect half-length L from the valley half-width (99th pct of
        // the perpendicular distance over the valley mask), matching the
        // straight path's half_length_multiplier * max_valley_width.
        std::vector<float> perp_in;
        perp_in.reserve(1024);
        for (std::size_t i = 0; i < vmask.size(); ++i)
            if (vmask[i] > 0.5f && i < perp_msk.size() && perp_msk[i])
                perp_in.push_back(perp_src[i]);
        double valley_half = 2.0 * cell;
        if (!perp_in.empty())
        {
            std::sort(perp_in.begin(), perp_in.end());
            const std::size_t idx = static_cast<std::size_t>(
                0.99 * (perp_in.size() - 1));
            valley_half = std::max<double>(perp_in[idx], 2.0 * cell);
        }
        const double max_width = 2.0 * valley_half;
        const double march_limit =
            std::max(request.max_half_length_metres, 4.0 * cell);
        double L = 0.5 * request.half_length_multiplier * max_width;
        L = std::max(L, request.min_half_length_metres);
        L = std::min(L, march_limit);

        result.channel_source_ids.push_back(ch.source_id);
        result.channel_half_length_metres.push_back(L);
        result.channel_max_valley_width_metres.push_back(max_width);

        // Corridor = valid pixels within L (+ a small pad) of the centreline.
        // Centreline pixels are held fixed so D stays anchored to the channel.
        const std::size_t npix = elev.size();
        std::vector<std::uint8_t> corridor(npix, 0);
        std::vector<std::uint8_t> fixed(npix, 0);
        const double corridor_lim = L + 3.0 * cell;
        for (std::size_t i = 0; i < npix; ++i)
            if (i < perp_msk.size() && perp_msk[i] && along_msk[i] &&
                perp_src[i] <= corridor_lim)
                corridor[i] = 1;
        for (std::size_t i = 0; i < nn; ++i)
        {
            const std::int64_t f = ch.centreline_flats[i];
            if (f >= 0 && static_cast<std::size_t>(f) < npix)
            {
                if (along_msk[static_cast<std::size_t>(f)])
                {
                    corridor[static_cast<std::size_t>(f)] = 1;
                    fixed[static_cast<std::size_t>(f)] = 1;
                }
            }
        }

        // Anchored Jacobi (Laplacian) smoothing of D over the corridor.
        std::vector<float> Dcur(along_src.begin(), along_src.end());
        std::vector<float> Dnext(Dcur);
        for (int pass = 0; pass < smooth_passes; ++pass)
        {
            for (int r = 0; r < n_rows; ++r)
            {
                for (int c = 0; c < n_cols; ++c)
                {
                    const std::size_t i =
                        static_cast<std::size_t>(r) * n_cols + c;
                    if (!corridor[i] || fixed[i]) continue;
                    double sum = 0.0; int cnt = 0;
                    if (r > 0)          { const std::size_t j = i - n_cols; if (corridor[j]) { sum += Dcur[j]; ++cnt; } }
                    if (r < n_rows - 1) { const std::size_t j = i + n_cols; if (corridor[j]) { sum += Dcur[j]; ++cnt; } }
                    if (c > 0)          { const std::size_t j = i - 1;      if (corridor[j]) { sum += Dcur[j]; ++cnt; } }
                    if (c < n_cols - 1) { const std::size_t j = i + 1;      if (corridor[j]) { sum += Dcur[j]; ++cnt; } }
                    if (cnt > 0)
                        Dnext[i] = static_cast<float>(sum / cnt);
                }
            }
            Dcur.swap(Dnext);
            Dnext = Dcur;
        }
        const std::vector<float>& D = Dcur;   // smoothed streamwise field

        // --- Trace both sides of every station as D-isocontours to length L --
        const std::size_t nst = station_along.size();
        std::vector<curved_side> left(nst), right(nst);
        std::vector<lsdtt_xy>    centre(nst);
        std::vector<char>        valid(nst, 0);
        for (std::size_t s = 0; s < nst; ++s)
        {
            centre[s] = { geom[s].cx, geom[s].cy };
            double level;
            if (!bilinear_raster(meta, D, corridor,
                                 geom[s].cx, geom[s].cy, level))
                continue;
            valid[s] = 1;
            left[s]  = trace_iso_side(meta, D, corridor, geom[s].cx, geom[s].cy,
                                      level, geom[s].nx, geom[s].ny, -1, step, L);
            right[s] = trace_iso_side(meta, D, corridor, geom[s].cx, geom[s].cy,
                                      level, geom[s].nx, geom[s].ny, +1, step, L);
            bool lf = false, rf = false;
            left[s].wall = find_valley_side(
                meta, vmask, geom[s].cx, geom[s].cy,
                geom[s].nx, geom[s].ny, -1, step, march_limit, gap_tol, lf);
            right[s].wall = find_valley_side(
                meta, vmask, geom[s].cx, geom[s].cy,
                geom[s].nx, geom[s].ny, +1, step, march_limit, gap_tol, rf);
            left[s].wall_found = lf; right[s].wall_found = rf;
        }

        auto floor_for = [&](double wall) -> double
        {
            return std::min(std::max((1.0 + clear) * wall, 2.0 * cell), L);
        };

        // Hierarchical (coarse-to-fine) acceptance. A single greedy sweep in
        // ALONG order lets a section bent around a bend clip out every station
        // downstream of it, leaving dense clusters next to empty gaps. Instead
        // lay down a sparse backbone first (stations every
        // curved_coarse_spacing_metres) then fill in the midpoints level by
        // level (coarse, coarse/2, coarse/4, ... down to the base spacing).
        // Each candidate is clipped at the Voronoi bisector of the ALREADY
        // accepted sections and dropped if it can no longer clear the valley or
        // would cross an accepted section. Coverage stays even: bends keep only
        // the coarse sections that fit; straight reaches fill in dense.
        struct curved_candidate
        {
            std::size_t s = 0, lk = 0, rk = 0;
            double lreach = 0.0, rreach = 0.0;
            std::vector<lsdtt_xy> test_poly;   // subsampled, for crossing tests
        };
        std::vector<curved_candidate> accepted_secs;
        const double clip_radius = 2.0 * L;

        // Dyadic insertion order over the base-spaced station indices: the
        // coarse grid first, then successive midpoints. C = coarse/base rounded
        // DOWN to a power of two so the halving is exact.
        std::size_t C = 1;
        {
            const double ratio = (spacing > 0.0)
                ? request.curved_coarse_spacing_metres / spacing : 1.0;
            const std::size_t craw =
                (ratio > 1.0) ? static_cast<std::size_t>(ratio + 0.5) : 1;
            while (C * 2 <= craw) C *= 2;
        }
        std::vector<std::size_t> order;
        order.reserve(nst);
        {
            std::vector<char> seen(nst, 0);
            auto add = [&](std::size_t i)
            { if (i < nst && !seen[i]) { seen[i] = 1; order.push_back(i); } };
            for (std::size_t i = 0; i < nst; i += C) add(i);   // coarse backbone
            for (std::size_t step = C; step > 1; step /= 2)    // fill midpoints
                for (std::size_t i = step / 2; i < nst; i += step) add(i);
            for (std::size_t i = 0; i < nst; ++i) add(i);      // any stragglers
        }

        for (std::size_t oi = 0; oi < order.size(); ++oi)
        {
            const std::size_t s = order[oi];
            if (!valid[s]) continue;

            // Clip both sides at the Voronoi bisector of nearby ACCEPTED
            // sections; drop if a side can no longer reach its wall floor.
            std::vector<lsdtt_xy> nbs;
            for (const curved_candidate& acc : accepted_secs)
            {
                const lsdtt_xy& q = centre[acc.s];
                if (std::hypot(q.x - centre[s].x, q.y - centre[s].y)
                        <= clip_radius)
                    nbs.push_back(q);
            }
            const std::size_t lc = nbs.empty()
                ? left[s].pts.size()
                : bisector_cross(left[s].pts, centre[s], nbs);
            const std::size_t rc = nbs.empty()
                ? right[s].pts.size()
                : bisector_cross(right[s].pts, centre[s], nbs);
            const std::size_t lk =
                (lc == 0) ? 0 : std::min(lc - 1, left[s].pts.size() - 1);
            const std::size_t rk =
                (rc == 0) ? 0 : std::min(rc - 1, right[s].pts.size() - 1);
            const double lreach = left[s].arc[lk];
            const double rreach = right[s].arc[rk];
            if (lreach < floor_for(left[s].wall) ||
                rreach < floor_for(right[s].wall))
                continue;

            // Subsampled polyline for the crossing test.
            curved_candidate cc;
            cc.s = s; cc.lk = lk; cc.rk = rk;
            cc.lreach = lreach; cc.rreach = rreach;
            const std::size_t stride = 3;
            for (std::size_t i = lk + 1; i-- > 0; )
                if (i == lk || i % stride == 0)
                    cc.test_poly.push_back(left[s].pts[i]);
            for (std::size_t i = 1; i <= rk; ++i)
                if (i == rk || i % stride == 0)
                    cc.test_poly.push_back(right[s].pts[i]);

            // Reject if it crosses any nearby accepted section (the isocontour
            // sides can converge at a bend's centre of curvature where no centre
            // clips them).
            bool crosses = false;
            for (std::size_t b = 0; b < accepted_secs.size() && !crosses; ++b)
            {
                const double dx = centre[cc.s].x - centre[accepted_secs[b].s].x;
                const double dy = centre[cc.s].y - centre[accepted_secs[b].s].y;
                if (std::hypot(dx, dy) > clip_radius) continue;
                if (sections_cross(cc.test_poly, accepted_secs[b].test_poly))
                    crosses = true;
            }
            if (crosses) continue;

            accepted_secs.push_back(std::move(cc));
        }

        // Emit accepted sections in ALONG order (station_seq increases
        // downstream); accepted_secs is currently in dyadic insertion order.
        std::sort(accepted_secs.begin(), accepted_secs.end(),
                  [](const curved_candidate& a, const curved_candidate& b)
                  { return a.s < b.s; });

        // Sampling pass: build the full curved cross section for each survivor.
        int station_seq_out = 0;
        for (std::size_t a = 0; a < accepted_secs.size(); ++a)
        {
            const std::size_t s  = accepted_secs[a].s;
            const std::size_t lk = accepted_secs[a].lk;
            const std::size_t rk = accepted_secs[a].rk;
            const double lreach = accepted_secs[a].lreach;
            const double rreach = accepted_secs[a].rreach;

            valley_cross_section xs;
            xs.source_id   = ch.source_id;
            xs.station_seq = static_cast<std::int32_t>(station_seq_out++);
            xs.along_distance_metres = station_along[s];
            xs.centre_x  = geom[s].cx;  xs.centre_y  = geom[s].cy;
            xs.normal_dx = geom[s].nx;  xs.normal_dy = geom[s].ny;
            xs.left_found  = left[s].wall_found;
            xs.right_found = right[s].wall_found;
            xs.left_offset_metres  = -lreach;
            xs.right_offset_metres =  rreach;
            xs.left_x  = left[s].pts[lk].x;   xs.left_y  = left[s].pts[lk].y;
            xs.right_x = right[s].pts[rk].x;  xs.right_y = right[s].pts[rk].y;
            xs.valley_width_metres = left[s].wall + right[s].wall;
            xs.left_truncated_by_neighbour  = (lk + 1 < left[s].pts.size());
            xs.right_truncated_by_neighbour = (rk + 1 < right[s].pts.size());

            // Assemble the curved polyline: left far-end -> centre -> right
            // far-end, with signed arc-length as the offset coordinate.
            xs.samples.reserve(lk + rk + 1);
            for (std::size_t i = lk + 1; i-- > 1; )
            {
                valley_cross_section_sample sp;
                sp.offset_metres = -left[s].arc[i];
                sp.x = left[s].pts[i].x; sp.y = left[s].pts[i].y;
                xs.samples.push_back(sp);
            }
            {
                valley_cross_section_sample sp;
                sp.offset_metres = 0.0;
                sp.x = geom[s].cx; sp.y = geom[s].cy;
                xs.samples.push_back(sp);
            }
            for (std::size_t i = 1; i <= rk; ++i)
            {
                valley_cross_section_sample sp;
                sp.offset_metres = right[s].arc[i];
                sp.x = right[s].pts[i].x; sp.y = right[s].pts[i].y;
                xs.samples.push_back(sp);
            }

            // Sample the DEM (and any optional rasters) along the polyline.
            for (valley_cross_section_sample& sp : xs.samples)
            {
                float e = 0.0f;
                sp.elevation_valid =
                    sample_dem(meta, elev, emask, sp.x, sp.y, e) ? 1 : 0;
                sp.elevation = e;
                sp.in_valley = in_valley_at(meta, vmask, sp.x, sp.y) ? 1 : 0;
                float cv = 0.0f;
                if (curvature && sample_raster(meta, *curvature, sp.x, sp.y, cv))
                { sp.curvature = cv; sp.curvature_valid = 1; }
                float hv = 0.0f;
                if (hillshade && sample_raster(meta, *hillshade, sp.x, sp.y, hv))
                { sp.hillshade = hv; sp.hillshade_valid = 1; }
                float chv = 0.0f;
                if (channel_mask &&
                    sample_raster(meta, *channel_mask, sp.x, sp.y, chv))
                    sp.is_channel = (chv > 0.5f) ? 1 : 0;
            }

            // --- Per-section feature detection (same as the straight path) ---
            const int ns = static_cast<int>(xs.samples.size());
            if (channel_mask)
            {
                int run_lo = -1;
                for (int j = 0; j <= ns; ++j)
                {
                    const bool cf = (j < ns) && xs.samples[j].is_channel;
                    if (cf && run_lo < 0) run_lo = j;
                    else if (!cf && run_lo >= 0)
                    {
                        const int mid = (run_lo + j - 1) / 2;
                        xs.channel_crossing_offsets.push_back(
                            xs.samples[mid].offset_metres);
                        run_lo = -1;
                    }
                }
            }
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

            xs.hydraulic_edges = hydraulic_edges_from_samples(
                xs, request,
                request.compute_hydraulic_depth ? &xs.hydraulic_curve
                                                 : nullptr,
                &xs.hydraulic_break_stages);

            result.sections.push_back(std::move(xs));
        }
    }

    return result;
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
    if (request.curved_sections)
        return extract_curved_impl(dem, centrelines, request,
                                   curvature, hillshade, channel_mask);

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
            // edges) from the sampled elevation profile. When requested, also
            // (f) keep the full BankfullMapper hydraulic-depth curve.
            xs.hydraulic_edges = hydraulic_edges_from_samples(
                xs, request,
                request.compute_hydraulic_depth ? &xs.hydraulic_curve
                                                 : nullptr,
                &xs.hydraulic_break_stages);

            result.sections.push_back(std::move(xs));
        }
    }

    return result;
}

// -----------------------------------------------------------------------------
// cross_section_informed_valley_mask
// -----------------------------------------------------------------------------
cross_section_informed_valley_result cross_section_informed_valley_mask(
    const lsdtt_raster_cpu&            dem,
    const valley_centreline_result&    centrelines,
    const valley_cross_section_result& sections)
{
    const lsdtt_raster_meta& meta = dem.get_meta();
    const std::vector<float>&        elev  = dem.get_data();
    const std::vector<std::uint8_t>& emask = dem.get_mask();
    const std::size_t n = elev.size();

    cross_section_informed_valley_result result;
    std::vector<float>        vm_data(n, 0.0f);
    std::vector<std::uint8_t> vm_mask(n, 0u);

    for (const valley_channel_centreline& ch : centrelines.channels)
    {
        const std::vector<std::int64_t>& nodes = ch.centreline_flats;
        if (nodes.size() < 2) continue;
        if (ch.centreline_baseline_elevation.size() != nodes.size()) continue;
        const std::vector<float>&        reach_d = ch.reach_id.get_data();
        const std::vector<std::uint8_t>& reach_m = ch.reach_id.get_mask();
        if (reach_d.size() != n || reach_m.empty()) continue;

        // Floodplain HAND threshold for this channel = the MEDIAN of the lowest
        // hydraulic break (break 1 = stages[1] - stages[0], the bankfull /
        // floodplain top) over all its sections. A single robust level per
        // channel avoids the lateral streaks and plateau leaks that a
        // per-section (interpolated) break produces where an individual section
        // has an anomalously high break. The second break is not used.
        std::vector<double> b1_hand;
        for (const valley_cross_section& xs : sections.sections)
        {
            if (xs.source_id != ch.source_id) continue;
            const std::vector<double>& bs = xs.hydraulic_break_stages;
            if (bs.size() >= 2) b1_hand.push_back(bs[1] - bs[0]);
        }
        if (b1_hand.empty()) continue;
        std::sort(b1_hand.begin(), b1_hand.end());
        const double break1_hand = b1_hand[b1_hand.size() / 2];  // median

        // Inflate the nearest centreline NODE INDEX across the corridor (the
        // corridor domain = the reach_id mask), so every corridor pixel gets a
        // baseline (thalweg) elevation and an along-channel position.
        std::vector<float>        sv(n, 0.0f);
        std::vector<std::uint8_t> sm(n, 0u);
        std::vector<float>        pd(n, 1.0f);
        std::vector<std::uint8_t> pm(n, 0u);
        for (std::size_t k = 0; k < nodes.size(); ++k)
        {
            const std::size_t f = static_cast<std::size_t>(nodes[k]);
            if (f < n) { sv[f] = static_cast<float>(k); sm[f] = 1u; }
        }
        for (std::size_t i = 0; i < n; ++i)
            pm[i] = (reach_m[i] || sm[i]) ? 1u : 0u;

        lsdtt_raster_cpu inflated = inflate_values_to_nearest_source(
            lsdtt_raster_cpu(std::move(sv), std::move(sm), meta),
            lsdtt_raster_cpu(std::move(pd), std::move(pm), meta));
        const std::vector<float>&        idx_d = inflated.get_data();
        const std::vector<std::uint8_t>& idx_m = inflated.get_mask();

        const int last = static_cast<int>(nodes.size()) - 1;
        std::vector<float>        cd(n, 0.0f);
        std::vector<std::uint8_t> cm(n, 0u);
        for (std::size_t i = 0; i < n; ++i)
        {
            if (!reach_m[i] || !idx_m[i]) continue;
            if (i >= emask.size() || !emask[i]) continue;
            const int k = std::clamp(
                static_cast<int>(std::llround(idx_d[i])), 0, last);
            const double baseline =
                static_cast<double>(ch.centreline_baseline_elevation[
                    static_cast<std::size_t>(k)]);
            if (!std::isfinite(baseline)) continue;
            const double rel = static_cast<double>(elev[i]) - baseline;

            // Floodplain iff within the channel's median lowest-break HAND.
            const std::uint8_t v = (rel <= break1_hand) ? 1 : 0;
            cd[i] = static_cast<float>(v);
            cm[i] = 1u;                                    // corridor pixel

            if (v != 0 && vm_data[i] == 0.0f)              // first claim wins
            {
                vm_data[i] = static_cast<float>(v);
                vm_mask[i] = 1u;
            }
        }

        result.channel_source_ids.push_back(ch.source_id);
        result.channel_masks.emplace_back(
            std::move(cd), std::move(cm), meta);
    }

    result.valley_mask =
        lsdtt_raster_cpu(std::move(vm_data), std::move(vm_mask), meta);
    return result;
}

} // namespace lsdtt
