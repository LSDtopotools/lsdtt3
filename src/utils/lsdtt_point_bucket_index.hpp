// lsdtt_point_bucket_index.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// A uniform-grid bucket index over a fixed set of 2-D points, for nearest-point
// and within-radius queries where a full kd-tree is not worth it: the points are
// polyline nodes on a raster, so they are roughly evenly spread and a bucket of a
// few cells holds a handful of them.
//
// nearest(x, y)  : index of the nearest point (ties: lower index), searching
//                  rings of buckets outward and stopping once the ring is
//                  farther than the best distance found. Exact.
// within(x, y, r, out) : appends the indices of every point within r (exact).
// -----------------------------------------------------------------------------

#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <vector>

namespace lsdtt
{

class point_bucket_index
{
public:
    // bucket_size: bucket edge length in map units (> 0).
    point_bucket_index(const std::vector<double>& x, const std::vector<double>& y,
                       double bucket_size)
        : x_(x), y_(y), size_(bucket_size > 0.0 ? bucket_size : 1.0)
    {
        if (x_.empty()) return;
        x0_ = *std::min_element(x_.begin(), x_.end());
        y0_ = *std::min_element(y_.begin(), y_.end());
        const double x1 = *std::max_element(x_.begin(), x_.end());
        const double y1 = *std::max_element(y_.begin(), y_.end());
        nx_ = static_cast<long>(std::floor((x1 - x0_) / size_)) + 1;
        ny_ = static_cast<long>(std::floor((y1 - y0_) / size_)) + 1;
        start_.assign(static_cast<std::size_t>(nx_ * ny_) + 1, 0);
        std::vector<long> b(x_.size());
        for (std::size_t i = 0; i < x_.size(); ++i)
        {
            b[i] = bucket_of(x_[i], y_[i]);
            ++start_[static_cast<std::size_t>(b[i]) + 1];
        }
        for (std::size_t k = 1; k < start_.size(); ++k) start_[k] += start_[k - 1];
        items_.resize(x_.size());
        std::vector<std::size_t> fill(start_.begin(), start_.end() - 1);
        for (std::size_t i = 0; i < x_.size(); ++i)
            items_[fill[static_cast<std::size_t>(b[i])]++] = i;
    }

    bool empty() const noexcept { return x_.empty(); }

    std::size_t nearest(double x, double y) const
    {
        std::size_t best = 0;
        double bd = std::numeric_limits<double>::infinity();
        if (x_.empty()) return best;
        const long cx = std::clamp(static_cast<long>(std::floor((x - x0_) / size_)), 0L, nx_ - 1);
        const long cy = std::clamp(static_cast<long>(std::floor((y - y0_) / size_)), 0L, ny_ - 1);
        // Distance from the query to the grid, for a query outside it.
        const double gx = std::max({x0_ - x, 0.0, x - (x0_ + nx_ * size_)});
        const double gy = std::max({y0_ - y, 0.0, y - (y0_ + ny_ * size_)});
        const long max_ring = std::max(nx_, ny_);
        for (long ring = 0; ring <= max_ring; ++ring)
        {
            // Every point in this ring or beyond is at least this far away
            // (the larger of the distance to the grid and ring - 1 buckets).
            const double ring_min = std::max(std::hypot(gx, gy),
                                             static_cast<double>(std::max(0L, ring - 1)) * size_);
            if (ring_min * ring_min > bd) break;
            for (long by = cy - ring; by <= cy + ring; ++by)
            {
                if (by < 0 || by >= ny_) continue;
                for (long bx = cx - ring; bx <= cx + ring; ++bx)
                {
                    if (bx < 0 || bx >= nx_) continue;
                    if (std::max(std::labs(bx - cx), std::labs(by - cy)) != ring) continue;
                    const std::size_t k = static_cast<std::size_t>(by * nx_ + bx);
                    for (std::size_t q = start_[k]; q < start_[k + 1]; ++q)
                    {
                        const std::size_t i = items_[q];
                        const double d = (x_[i] - x) * (x_[i] - x) + (y_[i] - y) * (y_[i] - y);
                        if (d < bd || (d == bd && i < best)) { bd = d; best = i; }
                    }
                }
            }
        }
        return best;
    }

    void within(double x, double y, double r, std::vector<std::size_t>& out) const
    {
        if (x_.empty()) return;
        const long bx0 = std::max(0L, static_cast<long>(std::floor((x - r - x0_) / size_)));
        const long bx1 = std::min(nx_ - 1, static_cast<long>(std::floor((x + r - x0_) / size_)));
        const long by0 = std::max(0L, static_cast<long>(std::floor((y - r - y0_) / size_)));
        const long by1 = std::min(ny_ - 1, static_cast<long>(std::floor((y + r - y0_) / size_)));
        const double r2 = r * r;
        for (long by = by0; by <= by1; ++by)
            for (long bx = bx0; bx <= bx1; ++bx)
            {
                const std::size_t k = static_cast<std::size_t>(by * nx_ + bx);
                for (std::size_t q = start_[k]; q < start_[k + 1]; ++q)
                {
                    const std::size_t i = items_[q];
                    if ((x_[i] - x) * (x_[i] - x) + (y_[i] - y) * (y_[i] - y) <= r2)
                        out.push_back(i);
                }
            }
    }

private:
    long bucket_of(double x, double y) const
    {
        const long bx = std::clamp(static_cast<long>(std::floor((x - x0_) / size_)), 0L, nx_ - 1);
        const long by = std::clamp(static_cast<long>(std::floor((y - y0_) / size_)), 0L, ny_ - 1);
        return by * nx_ + bx;
    }

    std::vector<double> x_, y_;
    double size_;
    double x0_ = 0.0, y0_ = 0.0;
    long nx_ = 0, ny_ = 0;
    std::vector<std::size_t> start_;
    std::vector<std::size_t> items_;
};

} // namespace lsdtt
