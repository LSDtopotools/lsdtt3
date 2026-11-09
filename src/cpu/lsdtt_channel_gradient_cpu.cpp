// lsdtt_channel_gradient_cpu.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
//
// -----------------------------------------------------------------------------
// See lsdtt_channel_gradient_cpu.hpp.
// -----------------------------------------------------------------------------

#include "lsdtt_channel_gradient_cpu.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>

namespace lsdtt
{

// Least-squares slope of z against along-channel distance d, over the path
// cells within +/- half_length_m of path index `idx`. Returns false (slope
// untouched) if the fit is not usable; see the header for the criteria.
// Elevation decreases downstream while d increases, so the positive
// down-channel slope is -cov(d, z) / var(d).
bool fit_along_path_slope(
    const std::vector<double>&       d,
    const std::vector<double>&       z,
    const std::vector<std::uint8_t>& ok,
    std::size_t                      idx,
    double                           half_length_m,
    double&                          slope,
    double*                          fitted_z)
{
    const double d0 = d[idx];

    const std::size_t lo = static_cast<std::size_t>(
        std::lower_bound(d.begin(), d.end(), d0 - half_length_m) - d.begin());
    const std::size_t hi = static_cast<std::size_t>(
        std::upper_bound(d.begin(), d.end(), d0 + half_length_m) - d.begin());

    std::size_t n = 0;
    double sum_d = 0.0, sum_z = 0.0;
    double d_min = std::numeric_limits<double>::max();
    double d_max = std::numeric_limits<double>::lowest();

    for (std::size_t i = lo; i < hi; ++i)
    {
        if (!ok[i]) continue;
        const double di = d[i] - d0;
        sum_d += di;
        sum_z += z[i];
        d_min = std::min(d_min, di);
        d_max = std::max(d_max, di);
        ++n;
    }

    if (n < 3) return false;
    if ((d_max - d_min) < half_length_m) return false;

    const double mean_d = sum_d / static_cast<double>(n);
    const double mean_z = sum_z / static_cast<double>(n);

    double sxx = 0.0, sxz = 0.0;
    for (std::size_t i = lo; i < hi; ++i)
    {
        if (!ok[i]) continue;
        const double dd = (d[i] - d0) - mean_d;
        sxx += dd * dd;
        sxz += dd * (z[i] - mean_z);
    }

    if (!(sxx > 0.0)) return false;

    const double s = -sxz / sxx;
    if (!std::isfinite(s) || s <= 0.0) return false;

    slope = s;
    // z = mean_z - s * (d - d0 - mean_d); at d = d0:
    if (fitted_z) *fitted_z = mean_z + s * mean_d;
    return true;
}

lsdtt_channel_slope_result channel_slope_at(
    const lsdtt_flow_base&  flow,
    const lsdtt_raster_cpu& filled_dem,
    std::int64_t            flat,
    double                  half_length_m,
    std::uint64_t           min_contributing_pixels)
{
    const double cs = flow.cell_size();
    if (!(half_length_m >= 2.0 * cs))
        throw std::invalid_argument(
            "channel_slope_at: the half-length must be at least two cell sizes (" +
            std::to_string(2.0 * cs) + " m); got " + std::to_string(half_length_m) + ".");
    if (!flow.is_flat_in_range(flat) || !flow.is_valid_cell(flat))
        throw std::invalid_argument("channel_slope_at: flat " + std::to_string(flat) +
                                    " is not a valid flow cell.");
    const int nc = flow.n_cols();
    const int nr = flow.n_rows();

    // Upstream: the largest-area donor at or above the threshold.
    std::vector<std::int64_t> up;          // upstream cells, nearest first
    std::vector<double> up_dist;           // distance upstream from `flat`
    {
        std::int64_t cur = flat;
        double dist = 0.0;
        for (std::int64_t guard = 0; guard < flow.n_valid_cells(); ++guard)
        {
            const int r = static_cast<int>(cur / nc), c = static_cast<int>(cur % nc);
            std::int64_t best = -1;
            std::uint64_t best_n = 0;
            for (int dr = -1; dr <= 1; ++dr)
                for (int dc = -1; dc <= 1; ++dc)
                {
                    if (!dr && !dc) continue;
                    const int rr = r + dr, cc = c + dc;
                    if (rr < 0 || rr >= nr || cc < 0 || cc >= nc) continue;
                    const std::int64_t v = static_cast<std::int64_t>(rr) * nc + cc;
                    if (!flow.is_valid_cell(v) || flow.is_outlet(v) || flow.receiver_of(v) != cur) continue;
                    const std::uint64_t n = flow.s_vector_size(v);
                    if (n < min_contributing_pixels) continue;
                    if (best < 0 || n > best_n || (n == best_n && v < best)) { best = v; best_n = n; }
                }
            if (best < 0) break;
            dist += d8_step_length(best, cur, nc, cs);
            if (dist > half_length_m) break;
            up.push_back(best);
            up_dist.push_back(dist);
            cur = best;
        }
    }
    // Downstream: receivers.
    std::vector<std::int64_t> down;
    std::vector<double> down_dist;
    {
        std::int64_t cur = flat;
        double dist = 0.0;
        for (std::int64_t guard = 0; guard < flow.n_valid_cells(); ++guard)
        {
            if (flow.is_outlet(cur)) break;
            const std::int64_t rec = flow.receiver_of(cur);
            if (rec < 0 || rec == cur || !flow.is_valid_cell(rec)) break;
            dist += d8_step_length(cur, rec, nc, cs);
            if (dist > half_length_m) break;
            down.push_back(rec);
            down_dist.push_back(dist);
            cur = rec;
        }
    }

    // Path in downstream order, distance from its upstream end.
    std::vector<std::int64_t> path;
    std::vector<double> d;
    const double top = up_dist.empty() ? 0.0 : up_dist.back();
    for (std::size_t i = up.size(); i-- > 0;)
    {
        path.push_back(up[i]);
        d.push_back(top - up_dist[i]);
    }
    const std::size_t idx = path.size();
    path.push_back(flat);
    d.push_back(top);
    for (std::size_t i = 0; i < down.size(); ++i)
    {
        path.push_back(down[i]);
        d.push_back(top + down_dist[i]);
    }

    const std::vector<float>& zd = filled_dem.get_data();
    const std::vector<std::uint8_t>& md = filled_dem.get_mask();
    std::vector<double> z(path.size(), 0.0);
    std::vector<std::uint8_t> ok(path.size(), 0);
    lsdtt_channel_slope_result out;
    double d_min = std::numeric_limits<double>::max(), d_max = std::numeric_limits<double>::lowest();
    for (std::size_t i = 0; i < path.size(); ++i)
    {
        const std::size_t f = static_cast<std::size_t>(path[i]);
        if (f >= zd.size() || !md[f] || !std::isfinite(zd[f])) continue;
        z[i] = zd[f];
        ok[i] = 1;
        ++out.n_cells;
        d_min = std::min(d_min, d[i]);
        d_max = std::max(d_max, d[i]);
    }
    if (out.n_cells > 0) out.window_length_m = d_max - d_min;
    double slope = 0.0, zfit = 0.0;
    out.slope_valid = fit_along_path_slope(d, z, ok, idx, half_length_m, slope, &zfit);
    if (out.slope_valid)
    {
        out.slope_m_per_m = slope;
        out.smoothed_elevation_m = zfit;
    }
    return out;
}

} // namespace lsdtt
