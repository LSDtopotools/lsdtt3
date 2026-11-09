// lsdtt_sea_mask_cpu.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
//
// -----------------------------------------------------------------------------
// See lsdtt_sea_mask_cpu.hpp.
// -----------------------------------------------------------------------------

#include "lsdtt_sea_mask_cpu.hpp"

#include <algorithm>
#include <cmath>
#include <deque>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>

namespace lsdtt
{

sea_side_choice parse_sea_side_choice(const std::string& t, const std::string& key)
{
    if (t == "auto") return sea_side_choice::automatic;
    if (t == "sea_inside") return sea_side_choice::sea_inside;
    if (t == "land_inside") return sea_side_choice::land_inside;
    if (t == "sea_left") return sea_side_choice::sea_left;
    if (t == "sea_right") return sea_side_choice::sea_right;
    throw std::invalid_argument(key + " must be one of auto, sea_inside, land_inside, sea_left, sea_right; got '" +
                                t + "'");
}

const char* sea_side_choice_name(sea_side_choice c)
{
    switch (c)
    {
    case sea_side_choice::automatic: return "auto";
    case sea_side_choice::sea_inside: return "sea_inside";
    case sea_side_choice::land_inside: return "land_inside";
    case sea_side_choice::sea_left: return "sea_left";
    case sea_side_choice::sea_right: return "sea_right";
    }
    return "auto";
}

namespace
{
struct side_stats
{
    long n = 0, n_nodata = 0;
    std::vector<float> z;
    double nodata_fraction() const { return n > 0 ? static_cast<double>(n_nodata) / n : 0.0; }
    double median()
    {
        if (z.empty()) return std::numeric_limits<double>::quiet_NaN();
        const size_t k = z.size() / 2;
        std::nth_element(z.begin(), z.begin() + static_cast<long>(k), z.end());
        return z[k];
    }
};

// +1: side A is sea; -1: side B is sea; 0: undecided. `why` gets the evidence.
int decide(side_stats& a, side_stats& b, const sea_mask_options& o, const std::string& name_a,
           const std::string& name_b, std::string& why)
{
    std::ostringstream w;
    w.precision(4);
    const double fa = a.nodata_fraction(), fb = b.nodata_fraction();
    if (a.n > 0 && b.n > 0 && std::fabs(fa - fb) >= o.nodata_evidence_fraction)
    {
        w << "NoData " << 100.0 * fa << "% " << name_a << " against " << 100.0 * fb << "% " << name_b;
        why = w.str();
        return fa > fb ? +1 : -1;
    }
    const double ma = a.median(), mb = b.median();
    if (std::isfinite(ma) && std::isfinite(mb) && std::fabs(ma - mb) >= o.min_contrast_metres)
    {
        w << "median elevation " << ma << " m " << name_a << " against " << mb << " m " << name_b;
        why = w.str();
        return ma < mb ? +1 : -1;
    }
    w << "no contrast: NoData " << 100.0 * fa << "% / " << 100.0 * fb << "%, median " << ma << " / " << mb
      << " m (" << name_a << " / " << name_b << ")";
    why = w.str();
    return 0;
}

struct grid
{
    int nr = 0, nc = 0;
    double x0 = 0, ymax = 0, cs = 1;
    lsdtt_xy centre(long f) const
    {
        return {x0 + (static_cast<double>(f % nc) + 0.5) * cs, ymax - (static_cast<double>(f / nc) + 0.5) * cs};
    }
};

// Nearest segment within `radius` of each cell: distance, side (+1 left, -1
// right, 0 on the line) and the part. Cells farther away keep part -1.
struct nearest_field
{
    std::vector<double> dist;
    std::vector<signed char> side;
    std::vector<int> part;
};

nearest_field nearest_side(const grid& g, const std::vector<std::vector<lsdtt_xy>>& lines, double radius)
{
    const size_t n = static_cast<size_t>(g.nr) * g.nc;
    nearest_field nf;
    nf.dist.assign(n, std::numeric_limits<double>::infinity());
    nf.side.assign(n, 0);
    nf.part.assign(n, -1);
    for (size_t p = 0; p < lines.size(); ++p)
        for (size_t i = 0; i + 1 < lines[p].size(); ++i)
        {
            const lsdtt_xy a = lines[p][i], b = lines[p][i + 1];
            const double dx = b.x - a.x, dy = b.y - a.y, L2 = dx * dx + dy * dy;
            if (!(L2 > 0.0)) continue;
            const int c0 = std::max(0, static_cast<int>(std::floor((std::min(a.x, b.x) - radius - g.x0) / g.cs)));
            const int c1 = std::min(g.nc - 1, static_cast<int>(std::floor((std::max(a.x, b.x) + radius - g.x0) / g.cs)));
            const int r0 = std::max(0, static_cast<int>(std::floor((g.ymax - std::max(a.y, b.y) - radius) / g.cs)));
            const int r1 = std::min(g.nr - 1, static_cast<int>(std::floor((g.ymax - std::min(a.y, b.y) + radius) / g.cs)));
            for (int r = r0; r <= r1; ++r)
                for (int c = c0; c <= c1; ++c)
                {
                    const long f = static_cast<long>(r) * g.nc + c;
                    const lsdtt_xy q = g.centre(f);
                    const double t = std::clamp(((q.x - a.x) * dx + (q.y - a.y) * dy) / L2, 0.0, 1.0);
                    const double d = std::hypot(q.x - (a.x + t * dx), q.y - (a.y + t * dy));
                    if (d > radius || d >= nf.dist[static_cast<size_t>(f)]) continue;
                    const double cross = dx * (q.y - a.y) - dy * (q.x - a.x);
                    nf.dist[static_cast<size_t>(f)] = d;
                    nf.side[static_cast<size_t>(f)] = cross > 0.0 ? 1 : (cross < 0.0 ? -1 : 0);
                    nf.part[static_cast<size_t>(f)] = static_cast<int>(p);
                }
        }
    return nf;
}

long count_valid(const std::vector<std::uint8_t>& sea, const std::vector<std::uint8_t>& mask)
{
    long n = 0;
    for (size_t i = 0; i < sea.size(); ++i) n += (sea[i] && mask[i]) ? 1 : 0;
    return n;
}

void refuse_if_everything(const sea_mask_result& r, const std::vector<std::uint8_t>& mask)
{
    long n_valid = 0;
    for (auto m : mask) n_valid += m ? 1 : 0;
    if (n_valid > 0 && r.n_sea_valid_cells == n_valid)
        throw std::invalid_argument("the sea mask would remove every valid cell of the DEM (side " + r.side_used +
                                    "); check sea_mask_side_choice");
}
} // namespace

sea_mask_result sea_mask_from_polygons(const lsdtt_raster_cpu& dem, const std::vector<std::uint8_t>& cover,
                                       const sea_mask_options& o)
{
    const auto& mask = dem.get_mask();
    const auto& z = dem.get_data();
    if (cover.size() != mask.size()) throw std::invalid_argument("sea_mask_from_polygons: cover does not match the DEM");
    sea_mask_result r;
    r.geometry = "polygons";
    long n_in = 0;
    for (auto c : cover) n_in += c ? 1 : 0;

    bool sea_in = true;
    switch (o.side)
    {
    case sea_side_choice::sea_inside: sea_in = true; r.side_used = "sea_inside"; break;
    case sea_side_choice::land_inside: sea_in = false; r.side_used = "land_inside"; break;
    case sea_side_choice::sea_left:
    case sea_side_choice::sea_right:
        throw std::invalid_argument("sea_mask_side_choice " + std::string(sea_side_choice_name(o.side)) +
                                    " is for lines; the sea mask file holds polygons (use sea_inside, land_inside or auto)");
    case sea_side_choice::automatic:
    {
        if (n_in == 0)
        {
            r.sea.assign(cover.size(), 0);
            r.side_used = "sea_inside";
            r.evidence = "the polygons cover no cell of the DEM";
            r.notes.push_back("the sea mask polygons do not cover the DEM: nothing removed");
            return r;
        }
        side_stats in, out;
        for (size_t i = 0; i < cover.size(); ++i)
        {
            side_stats& s = cover[i] ? in : out;
            ++s.n;
            if (!mask[i]) ++s.n_nodata;
            else s.z.push_back(z[i]);
        }
        std::string why;
        const int d = decide(in, out, o, "inside", "outside", why);
        if (d == 0)
            throw std::invalid_argument("sea_mask_side_choice auto cannot tell sea from land (" + why +
                                        "); give sea_inside or land_inside");
        sea_in = d > 0;
        r.side_used = sea_in ? "sea_inside" : "land_inside";
        r.evidence = why;
        break;
    }
    }
    r.sea.resize(cover.size());
    for (size_t i = 0; i < cover.size(); ++i) r.sea[i] = (cover[i] != 0) == sea_in ? 1 : 0;
    r.n_sea_valid_cells = count_valid(r.sea, mask);
    refuse_if_everything(r, mask);
    return r;
}

sea_mask_result sea_mask_from_lines(const lsdtt_raster_cpu& dem, const std::vector<std::uint8_t>& barrier,
                                    const std::vector<std::vector<lsdtt_xy>>& lines, const sea_mask_options& o)
{
    const auto& meta = dem.get_meta();
    const auto& mask = dem.get_mask();
    const auto& z = dem.get_data();
    if (barrier.size() != mask.size()) throw std::invalid_argument("sea_mask_from_lines: line cells do not match the DEM");
    if (o.side == sea_side_choice::sea_inside || o.side == sea_side_choice::land_inside)
        throw std::invalid_argument("sea_mask_side_choice " + std::string(sea_side_choice_name(o.side)) +
                                    " is for polygons; the sea mask file holds lines (use sea_left, sea_right or auto)");
    if (!(o.max_distance_metres > 0.0)) throw std::invalid_argument("sea_mask_max_distance_metres must be > 0");
    if (o.band_cells < 1) throw std::invalid_argument("sea_mask: band_cells must be >= 1");

    grid g;
    g.nr = meta.n_rows; g.nc = meta.n_cols; g.x0 = meta.x_min; g.ymax = meta.y_max(); g.cs = meta.cell_size;
    const size_t n = mask.size();
    sea_mask_result r;
    r.geometry = "lines";
    r.sea.assign(n, 0);

    const double radius = std::max(1.5, static_cast<double>(o.band_cells)) * g.cs;
    const nearest_field nf = nearest_side(g, lines, radius);
    const size_t np = lines.size();

    // sea side per part: +1 left, -1 right
    std::vector<int> sea_side(np, 0);
    if (o.side == sea_side_choice::sea_left) { std::fill(sea_side.begin(), sea_side.end(), +1); r.side_used = "sea_left"; }
    else if (o.side == sea_side_choice::sea_right) { std::fill(sea_side.begin(), sea_side.end(), -1); r.side_used = "sea_right"; }
    else
    {
        std::vector<side_stats> left(np), right(np);
        for (size_t i = 0; i < n; ++i)
        {
            const int p = nf.part[i];
            if (p < 0 || barrier[i] || nf.side[i] == 0 || nf.dist[i] > o.band_cells * g.cs) continue;
            side_stats& s = nf.side[i] > 0 ? left[static_cast<size_t>(p)] : right[static_cast<size_t>(p)];
            ++s.n;
            if (!mask[i]) ++s.n_nodata;
            else s.z.push_back(z[i]);
        }
        long w_left = 0, w_right = 0;
        int n_undecided = 0;
        std::string first_why;
        for (size_t p = 0; p < np; ++p)
        {
            if (left[p].n == 0 && right[p].n == 0) continue;
            std::string why;
            const int d = decide(left[p], right[p], o, "left", "right", why);
            if (first_why.empty() && d != 0) first_why = why;
            sea_side[p] = d;
            if (d > 0) w_left += left[p].n + right[p].n;
            else if (d < 0) w_right += left[p].n + right[p].n;
            else ++n_undecided;
        }
        if (w_left == 0 && w_right == 0)
            throw std::invalid_argument("sea_mask_side_choice auto cannot tell sea from land beside the lines (" +
                                        (first_why.empty() ? std::string("no cells beside them") : first_why) +
                                        "); give sea_left or sea_right");
        const int majority = w_left >= w_right ? +1 : -1;
        int n_flipped = 0;
        for (size_t p = 0; p < np; ++p)
        {
            if (sea_side[p] == 0) sea_side[p] = majority;
            else if (sea_side[p] != majority) ++n_flipped;
        }
        std::ostringstream w;
        w << (majority > 0 ? "sea on the left" : "sea on the right") << " of most lines (" << first_why << ")";
        if (n_flipped) w << "; " << n_flipped << " line part(s) with the sea on the other side";
        if (n_undecided) w << "; " << n_undecided << " part(s) without contrast took the majority";
        r.evidence = w.str();
        r.side_used = n_flipped ? "auto (per line)" : (majority > 0 ? "sea_left" : "sea_right");
    }

    auto is_seed = [&](size_t i, bool& sea)
    {
        if (!mask[i] || barrier[i] || nf.part[i] < 0 || nf.side[i] == 0) return false;
        const long f = static_cast<long>(i);
        const int rr = static_cast<int>(f / g.nc), cc = static_cast<int>(f % g.nc);
        bool beside = false;
        const int dr[4] = {-1, 1, 0, 0}, dc[4] = {0, 0, -1, 1};
        for (int k = 0; k < 4 && !beside; ++k)
        {
            const int r2 = rr + dr[k], c2 = cc + dc[k];
            if (r2 < 0 || c2 < 0 || r2 >= g.nr || c2 >= g.nc) continue;
            beside = barrier[static_cast<size_t>(r2) * g.nc + c2] != 0;
        }
        if (!beside) return false;
        sea = nf.side[i] == sea_side[static_cast<size_t>(nf.part[i])];
        return true;
    };

    std::vector<std::uint8_t> flooded(n, 0), land_seed(n, 0);
    std::deque<long> q;
    long n_sea_seeds = 0;
    for (size_t i = 0; i < n; ++i)
    {
        bool sea = false;
        if (!is_seed(i, sea)) continue;
        if (sea) { flooded[i] = 1; q.push_back(static_cast<long>(i)); ++n_sea_seeds; }
        else land_seed[i] = 1;
    }
    if (n_sea_seeds == 0)
    {
        r.notes.push_back("no valid cell beside the sea side of the lines: nothing removed");
        return r;
    }
    long n_leak = 0;
    while (!q.empty())
    {
        const long f = q.front();
        q.pop_front();
        if (land_seed[static_cast<size_t>(f)]) ++n_leak;
        const int rr = static_cast<int>(f / g.nc), cc = static_cast<int>(f % g.nc);
        const int dr[4] = {-1, 1, 0, 0}, dc[4] = {0, 0, -1, 1};
        for (int k = 0; k < 4; ++k)
        {
            const int r2 = rr + dr[k], c2 = cc + dc[k];
            if (r2 < 0 || c2 < 0 || r2 >= g.nr || c2 >= g.nc) continue;
            const size_t v = static_cast<size_t>(r2) * g.nc + c2;
            if (flooded[v] || !mask[v] || barrier[v]) continue;
            flooded[v] = 1;
            q.push_back(static_cast<long>(v));
        }
    }

    if (n_leak == 0)
    {
        for (size_t i = 0; i < n; ++i)
        {
            if (flooded[i]) r.sea[i] = 1;
            else if (barrier[i] && nf.part[i] >= 0 && nf.side[i] != 0 &&
                     nf.side[i] == sea_side[static_cast<size_t>(nf.part[i])])
                r.sea[i] = 1;
        }
    }
    else
    {
        // A gap, or lines that do not cross the DEM: the flood reached the
        // land side. Fall back to the sea side within max_distance_metres.
        r.leaked = true;
        r.fell_back = true;
        const nearest_field band = nearest_side(g, lines, o.max_distance_metres);
        for (size_t i = 0; i < n; ++i)
            if (band.part[i] >= 0 && band.side[i] != 0 && band.side[i] == sea_side[static_cast<size_t>(band.part[i])])
                r.sea[i] = 1;
        std::ostringstream w;
        w << "the flood from the sea side reached " << n_leak << " cell(s) beside the land side: the lines have a gap "
          << "or do not cross the DEM. Fell back to the sea side within " << o.max_distance_metres << " m of the lines";
        r.notes.push_back(w.str());
    }
    r.n_sea_valid_cells = count_valid(r.sea, mask);
    refuse_if_everything(r, mask);
    return r;
}

std::size_t apply_sea_mask(lsdtt_raster_cpu& dem, const std::vector<std::uint8_t>& sea)
{
    auto& d = dem.get_data_array();
    auto& m = dem.get_mask_array();
    const float nd = dem.get_meta().no_data_value;
    std::size_t n = 0;
    for (int r = 0; r < dem.get_meta().n_rows; ++r)
        for (int c = 0; c < dem.get_meta().n_cols; ++c)
        {
            const size_t i = static_cast<size_t>(r) * dem.get_meta().n_cols + c;
            if (i < sea.size() && sea[i] && m(r, c))
            {
                d(r, c) = nd;
                m(r, c) = 0;
                ++n;
            }
        }
    return n;
}

} // namespace lsdtt
