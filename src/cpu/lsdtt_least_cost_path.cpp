// lsdtt_least_cost_path.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html

#include "lsdtt_least_cost_path.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <queue>
#include <stdexcept>
#include <unordered_set>
#include <vector>

#include "../core/lsdtt_raster_meta.hpp"

namespace lsdtt
{

// ---------------------------------------------------------------------------
// raster_least_cost_path
// ---------------------------------------------------------------------------
lcp_result raster_least_cost_path(
    const lsdtt_raster_cpu& dem,
    int                     start_row,
    int                     start_col,
    int                     goal_row,
    int                     goal_col,
    const lcp_params&       params)
{
    const lsdtt_raster_meta& meta = dem.get_meta();
    const int    nr = meta.n_rows;
    const int    nc = meta.n_cols;
    const double cs = meta.cell_size;

    const std::vector<float>&   elev = dem.get_data();
    const std::vector<uint8_t>& mask = dem.get_mask();

    // Resolve window bounds
    const int r_min = (params.r_min >= 0) ? params.r_min : 0;
    const int r_max = (params.r_max >= 0) ? std::min(params.r_max, nr-1) : nr-1;
    const int c_min = (params.c_min >= 0) ? params.c_min : 0;
    const int c_max = (params.c_max >= 0) ? std::min(params.c_max, nc-1) : nc-1;

    // Validate start
    auto in_bounds = [&](int r, int c)
    {
        return r >= 0 && r < nr && c >= 0 && c < nc;
    };

    if (!in_bounds(start_row, start_col))
        throw std::invalid_argument("raster_least_cost_path: start pixel out of range");

    const int start_flat = start_row * nc + start_col;
    if (!mask[start_flat])
        throw std::invalid_argument("raster_least_cost_path: start pixel is nodata");

    // Build goal set (multi-target) or validate single goal.
    // Multi-target goals are filtered to valid, in-window DEM pixels.
    const bool multi_target_requested = !params.goal_flats.empty();
    std::unordered_set<int> goal_set;

    auto in_window = [&](int r, int c)
    {
        return r >= r_min && r <= r_max && c >= c_min && c <= c_max;
    };

    if (!in_window(start_row, start_col))
        throw std::invalid_argument(
            "raster_least_cost_path: start pixel lies outside search window");

    const std::int64_t N64 = static_cast<std::int64_t>(nr) * nc;

    if (multi_target_requested)
    {
        goal_set.reserve(params.goal_flats.size());
        for (std::int64_t gf64 : params.goal_flats)
        {
            if (gf64 < 0 || gf64 >= N64) continue;
            const int gf = static_cast<int>(gf64);
            const int gr = gf / nc;
            const int gc = gf % nc;
            if (!in_window(gr, gc)) continue;
            if (!mask[gf])          continue;
            goal_set.insert(gf);
        }
        if (goal_set.empty())
            return lcp_result{};
    }
    else
    {
        if (!in_bounds(goal_row, goal_col))
            throw std::invalid_argument("raster_least_cost_path: goal pixel out of range");
        if (!in_window(goal_row, goal_col))
            throw std::invalid_argument(
                "raster_least_cost_path: goal pixel lies outside search window");
    }

    const bool multi_target = !goal_set.empty();

    // Channel score pointers (null when feature disabled)
    const std::vector<float>* scores = params.channel_score;
    const double              alpha  = (scores && params.channel_score_alpha > 0.0)
                                       ? params.channel_score_alpha : 0.0;

    const int N = nr * nc;
    std::vector<double> dist(N, std::numeric_limits<double>::infinity());
    std::vector<int>    prev(N, -1);

    dist[start_flat] = 0.0;

    using Node = std::pair<double, int>;
    std::priority_queue<Node, std::vector<Node>, std::greater<Node>> pq;
    pq.push({0.0, start_flat});

    static constexpr int DR[8] = {-1,-1,-1, 0, 0, 1, 1, 1};
    static constexpr int DC[8] = {-1, 0, 1,-1, 1,-1, 0, 1};

    const double gr2    = static_cast<double>(params.goal_radius)
                        * static_cast<double>(params.goal_radius);
    const double uphill = params.uphill_mult;

    int    goal_flat = -1;
    double goal_cost = std::numeric_limits<double>::infinity();

    while (!pq.empty())
    {
        auto [cd, cf] = pq.top(); pq.pop();
        if (cd > dist[cf]) continue;

        const int cr = cf / nc;
        const int cc = cf % nc;

        // Goal check
        if (multi_target)
        {
            if (goal_set.count(cf))
            {
                goal_flat = cf;
                goal_cost = cd;
                break;
            }
        }
        else
        {
            const double dgr = cr - goal_row;
            const double dgc = cc - goal_col;
            if (dgr*dgr + dgc*dgc <= gr2)
            {
                goal_flat = cf;
                goal_cost = cd;
                break;
            }
        }

        for (int k = 0; k < 8; ++k)
        {
            const int nr2 = cr + DR[k];
            const int nc2 = cc + DC[k];
            if (nr2 < r_min || nr2 > r_max || nc2 < c_min || nc2 > c_max) continue;

            const int nf = nr2 * nc + nc2;
            if (!mask[nf]) continue;

            const bool   is_diag = (DR[k] != 0 && DC[k] != 0);
            const double step    = is_diag ? cs * 1.41421356237 : cs;

            // Channel-score weighting: pixels with high channel score are cheaper.
            double channel_factor = 1.0;
            if (alpha > 0.0 && scores)
            {
                const double s = static_cast<double>((*scores)[nf]);
                channel_factor = 1.0 + alpha * (1.0 - s);
            }

            const double ed    = static_cast<double>(elev[nf])
                               - static_cast<double>(elev[cf]);
            const double ncost = cd + step * channel_factor
                               + std::max(0.0, ed) * uphill;

            if (ncost < dist[nf])
            {
                dist[nf] = ncost;
                prev[nf] = cf;
                pq.push({ncost, nf});
            }
        }
    }

    lcp_result result;
    if (goal_flat < 0)
        return result;  // not found

    // Backtrack path and compute max_climb
    const double start_elev = elev[start_flat];
    double max_climb = 0.0;

    std::vector<std::int64_t> path;
    for (int f = goal_flat; f >= 0; f = prev[f])
    {
        path.push_back(static_cast<std::int64_t>(f));
        const double climb = static_cast<double>(elev[f]) - start_elev;
        if (climb > max_climb) max_climb = climb;
    }
    std::reverse(path.begin(), path.end());

    result.path       = std::move(path);
    result.total_cost = goal_cost;
    result.max_climb  = max_climb;
    result.found      = true;
    return result;
}

// ---------------------------------------------------------------------------
// burn_monotonic_path
// ---------------------------------------------------------------------------
void burn_monotonic_path(
    lsdtt_raster_cpu&                   dem,
    const std::vector<std::int64_t>&    path,
    double                              min_slope)
{
    if (path.empty()) return;

    const lsdtt_raster_meta& meta = dem.get_meta();
    std::vector<float>& elev = dem.get_data();

    auto step_length = [&](std::int64_t a, std::int64_t b) -> double
    {
        const int ar = static_cast<int>(a / meta.n_cols);
        const int ac = static_cast<int>(a % meta.n_cols);
        const int br = static_cast<int>(b / meta.n_cols);
        const int bc = static_cast<int>(b % meta.n_cols);
        const int dr = std::abs(ar - br);
        const int dc = std::abs(ac - bc);
        return (dr == 1 && dc == 1) ? meta.cell_size * 1.4142135623730950488
                                    : meta.cell_size;
    };

    const double e0 = elev[static_cast<std::size_t>(path.front())];
    const double eN = elev[static_cast<std::size_t>(path.back())];

    auto scan = [&](auto begin_it, auto end_it)
    {
        if (begin_it == end_it) return;
        std::int64_t prev_flat = *begin_it;
        double allowed = elev[static_cast<std::size_t>(prev_flat)];

        for (auto it = begin_it; it != end_it; ++it)
        {
            const std::int64_t flat = *it;
            if (it != begin_it)
            {
                const double ds = step_length(prev_flat, flat);
                allowed -= std::max(0.0, min_slope) * ds;
            }
            const double z = elev[static_cast<std::size_t>(flat)];
            if (z < allowed)
                allowed = z;
            else if (z > allowed)
                elev[static_cast<std::size_t>(flat)] = static_cast<float>(allowed);
            prev_flat = flat;
        }
    };

    if (e0 >= eN)
        scan(path.begin(), path.end());
    else
        scan(path.rbegin(), path.rend());
}

} // namespace lsdtt
