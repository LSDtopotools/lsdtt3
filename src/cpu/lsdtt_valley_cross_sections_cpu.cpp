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
        struct station_geom { double cx, cy, nx, ny; std::size_t node; };
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

        // Pass 2: sample each transect over [-half_len, +half_len].
        for (std::size_t s = 0; s < station_along.size(); ++s)
        {
            valley_cross_section xs;
            xs.source_id   = ch.source_id;
            xs.station_seq = static_cast<std::int32_t>(s);
            xs.along_distance_metres = station_along[s];
            xs.centre_x = geom[s].cx;
            xs.centre_y = geom[s].cy;
            xs.normal_dx = geom[s].nx;
            xs.normal_dy = geom[s].ny;

            xs.left_found  = left_ok[s];
            xs.right_found = right_ok[s];
            xs.left_offset_metres  = -left_off[s];
            xs.right_offset_metres =  right_off[s];
            xs.left_x  = geom[s].cx - left_off[s]  * geom[s].nx;
            xs.left_y  = geom[s].cy - left_off[s]  * geom[s].ny;
            xs.right_x = geom[s].cx + right_off[s] * geom[s].nx;
            xs.right_y = geom[s].cy + right_off[s] * geom[s].ny;
            xs.valley_width_metres = left_off[s] + right_off[s];

            const int n_side =
                static_cast<int>(std::floor(half_len / step + 0.5));
            xs.samples.reserve(2 * n_side + 1);
            for (int j = -n_side; j <= n_side; ++j)
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

            result.sections.push_back(std::move(xs));
        }
    }

    return result;
}

} // namespace lsdtt
