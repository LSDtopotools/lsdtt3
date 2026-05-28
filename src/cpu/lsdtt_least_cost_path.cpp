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

    // Build goal set (multi-target) or validate single goal
    const bool multi_target = !params.goal_flats.empty();
    std::unordered_set<int> goal_set;

    if (multi_target)
    {
        goal_set.reserve(params.goal_flats.size());
        for (std::int64_t gf : params.goal_flats)
            goal_set.insert(static_cast<int>(gf));
    }
    else
    {
        if (!in_bounds(goal_row, goal_col))
            throw std::invalid_argument("raster_least_cost_path: goal pixel out of range");
    }

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
    const std::vector<std::int64_t>&    path)
{
    if (path.empty()) return;

    std::vector<float>& elev = dem.get_data();

    const double e0 = elev[path.front()];
    const double eN = elev[path.back()];

    if (e0 >= eN)
    {
        // Forward: path[0] is the upstream (high) end
        double rmin = e0;
        for (std::int64_t f : path)
        {
            const double e = elev[f];
            if (e < rmin)
                rmin = e;
            else if (e > rmin)
                elev[f] = static_cast<float>(rmin);
        }
    }
    else
    {
        // Reverse: path[N-1] is the upstream (high) end
        double rmin = eN;
        for (int i = static_cast<int>(path.size()) - 1; i >= 0; --i)
        {
            const double e = elev[path[i]];
            if (e < rmin)
                rmin = e;
            else if (e > rmin)
                elev[path[i]] = static_cast<float>(rmin);
        }
    }
}

} // namespace lsdtt
