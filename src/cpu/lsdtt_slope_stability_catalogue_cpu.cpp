// lsdtt_slope_stability_catalogue_cpu.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// Landslide catalogues. Grouping rules and per-landslide quantities: see the header.

#include "lsdtt_slope_stability_catalogue_cpu.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

#include "../utils/lsdtt_connected_components.hpp"

namespace lsdtt
{

namespace
{

struct accumulator
{
    std::int64_t n = 0;
    double sum_h = 0.0, max_h = 0.0, sum_f = 0.0;
    double resist = 0.0, drive = 0.0;
    double sum_x = 0.0, sum_y = 0.0, sum_xx = 0.0, sum_yy = 0.0, sum_xy = 0.0;   // in cells
    std::int64_t boundary_edges = 0;
};

void finish(landslide_record& rec, const accumulator& a, const lsdtt_raster_meta& meta)
{
    const double b = meta.cell_size;
    rec.n_cells = a.n;
    rec.area_m2 = static_cast<double>(a.n) * b * b;
    rec.mean_depth_m = a.sum_h / static_cast<double>(a.n);
    rec.max_depth_m = a.max_h;
    rec.volume_m3 = rec.area_m2 * rec.mean_depth_m;
    rec.mean_safety_factor = a.sum_f / static_cast<double>(a.n);
    rec.landslide_safety_factor = (a.drive > 0.0) ? a.resist / a.drive : std::numeric_limits<double>::infinity();

    const double n = static_cast<double>(a.n);
    const double mc = a.sum_x / n, mr = a.sum_y / n;               // column, row
    rec.centroid_x = meta.x_min + (mc + 0.5) * b;
    rec.centroid_y = meta.y_min + (static_cast<double>(meta.n_rows) - mr - 0.5) * b;
    // Second central moments in cells, plus 1/12 for the cell's own extent.
    const double uxx = a.sum_xx / n - mc * mc + 1.0 / 12.0;
    const double uyy = a.sum_yy / n - mr * mr + 1.0 / 12.0;
    const double uxy = a.sum_xy / n - mc * mr;
    const double common = std::sqrt((uxx - uyy) * (uxx - uyy) + 4.0 * uxy * uxy);
    rec.length_m = 2.0 * std::sqrt(2.0) * std::sqrt(std::max(0.0, uxx + uyy + common)) * b;
    rec.width_m  = 2.0 * std::sqrt(2.0) * std::sqrt(std::max(0.0, uxx + uyy - common)) * b;
    rec.length_to_width = (rec.width_m > 0.0) ? rec.length_m / rec.width_m : std::numeric_limits<double>::infinity();
    rec.perimeter_m = static_cast<double>(a.boundary_edges) * b;
}

} // namespace

landslide_catalogues build_landslide_catalogues(
    const rupture_rasters& rupture,
    const lsdtt_raster_meta& meta,
    const std::vector<std::int32_t>& region,
    const std::map<std::int32_t, rupture_parameters>& by_region,
    int min_daylight_neighbours)
{
    const int nr = meta.n_rows, nc = meta.n_cols;
    const std::size_t n = static_cast<std::size_t>(nr) * static_cast<std::size_t>(nc);
    if (rupture.daylight.size() != n || rupture.safety_factor.get_data().size() != n)
        throw std::invalid_argument("build_landslide_catalogues: the rupture result differs from the grid");
    if (!region.empty() && region.size() != n)
        throw std::invalid_argument("build_landslide_catalogues: the region grid differs from the grid");
    if (min_daylight_neighbours < 0 || min_daylight_neighbours > 8)
        throw std::invalid_argument("build_landslide_catalogues: min_daylight_neighbours must be 0 to 8");
    if (by_region.empty())
        throw std::invalid_argument("build_landslide_catalogues: no parameter set given");

    const std::vector<float>& F = rupture.safety_factor.get_data();
    const std::vector<std::uint8_t>& Fm = rupture.safety_factor.get_mask();
    const std::vector<float>& H = rupture.depth.get_data();
    const std::vector<float>& Srup = rupture.rupture_slope.get_data();

    // Step 2: labels with "add outlets".
    std::vector<std::int64_t> lab(rupture.daylight);
    std::vector<std::uint8_t> is_daylight(n, 0);
    for (std::size_t i = 0; i < n; ++i)
        if (rupture.daylight[i] >= 0) is_daylight[static_cast<std::size_t>(rupture.daylight[i])] = 1;
    for (std::size_t d = 0; d < n; ++d)
        if (is_daylight[d]) lab[d] = static_cast<std::int64_t>(d);

    landslide_catalogues out;

    // Step 3: the 1-D rule.
    std::vector<std::uint8_t> keep_d(n, 0);
    for (std::size_t d = 0; d < n; ++d)
    {
        if (!is_daylight[d]) continue;
        ++out.n_daylight_points;
        const int r = static_cast<int>(d / nc), c = static_cast<int>(d % nc);
        int count = 0;
        for (int dr = -1; dr <= 1; ++dr)
            for (int dc = -1; dc <= 1; ++dc)
            {
                if (dr == 0 && dc == 0) continue;
                const int rr = r + dr, cc = c + dc;
                if (rr < 0 || cc < 0 || rr >= nr || cc >= nc) continue;
                if (lab[static_cast<std::size_t>(rr) * nc + cc] == static_cast<std::int64_t>(d)) ++count;
            }
        if (count >= min_daylight_neighbours) keep_d[d] = 1;
        else ++out.n_daylight_points_removed;
    }

    // Step 4: landslide cells.
    std::vector<std::uint8_t> slide(n, 0);
    for (std::size_t i = 0; i < n; ++i)
    {
        if (lab[i] < 0 || !keep_d[static_cast<std::size_t>(lab[i])]) continue;
        if (!Fm[i]) continue;
        const double f = F[i];
        if (std::isfinite(f) && f <= 1.0) { slide[i] = 1; ++out.n_landslide_cells; }
    }

    // Step 5: by-daylight ids in increasing flat index of D.
    std::vector<std::int32_t> id_of_d(n, 0);
    std::int32_t next = 0;
    {
        std::vector<std::uint8_t> used(n, 0);
        for (std::size_t i = 0; i < n; ++i)
            if (slide[i]) used[static_cast<std::size_t>(lab[i])] = 1;
        for (std::size_t d = 0; d < n; ++d)
            if (used[d]) id_of_d[d] = ++next;
    }
    out.label.assign(n, 0);
    for (std::size_t i = 0; i < n; ++i)
        if (slide[i]) out.label[i] = id_of_d[static_cast<std::size_t>(lab[i])];

    // Step 6: amalgamation.
    std::vector<std::int64_t> cc_size;
    const std::int32_t n_cc = label_connected_components(
        nr, nc, 8, [&](std::size_t i) { return slide[i] != 0; }, out.label_amalgamated, cc_size);

    // Per-landslide accumulation.
    std::vector<accumulator> acc_d(static_cast<std::size_t>(next) + 1), acc_a(static_cast<std::size_t>(n_cc) + 1);
    const rupture_parameters* single = region.empty() ? &by_region.begin()->second : nullptr;
    for (std::size_t i = 0; i < n; ++i)
    {
        if (!slide[i]) continue;
        const int r = static_cast<int>(i / nc), c = static_cast<int>(i % nc);
        const rupture_parameters* pp = single;
        if (!pp)
        {
            auto it = by_region.find(region[i]);
            if (it == by_region.end())
                throw std::invalid_argument("build_landslide_catalogues: a landslide cell's region has no parameter row");
            pp = &it->second;
        }
        const double h = H[i];
        const double t = std::atan(static_cast<double>(Srup[i]));
        const double rho_g = pp->density_kg_m3 * pp->gravity_m_s2;
        const double resist = pp->cohesion_pa + pp->friction_coefficient * rho_g * h * std::cos(t) * std::cos(t);
        const double drive = rho_g * h * std::cos(t) * std::sin(t);
        for (int which = 0; which < 2; ++which)
        {
            const std::int32_t lid = (which == 0) ? out.label[i] : out.label_amalgamated[i];
            accumulator& a = (which == 0) ? acc_d[static_cast<std::size_t>(lid)] : acc_a[static_cast<std::size_t>(lid)];
            const std::vector<std::int32_t>& L = (which == 0) ? out.label : out.label_amalgamated;
            ++a.n;
            a.sum_h += h;
            a.max_h = std::max(a.max_h, h);
            a.sum_f += F[i];
            a.resist += resist;
            a.drive += drive;
            a.sum_x += c; a.sum_y += r;
            a.sum_xx += static_cast<double>(c) * c;
            a.sum_yy += static_cast<double>(r) * r;
            a.sum_xy += static_cast<double>(c) * r;
            const int d4r[4] = {-1, 1, 0, 0}, d4c[4] = {0, 0, -1, 1};
            for (int k = 0; k < 4; ++k)
            {
                const int rr = r + d4r[k], cc2 = c + d4c[k];
                if (rr < 0 || cc2 < 0 || rr >= nr || cc2 >= nc || L[static_cast<std::size_t>(rr) * nc + cc2] != lid)
                    ++a.boundary_edges;
            }
        }
    }

    const double b = meta.cell_size;
    out.by_daylight.resize(static_cast<std::size_t>(next));
    for (std::size_t d = 0; d < n; ++d)
    {
        const std::int32_t lid = id_of_d[d];
        if (!lid) continue;
        landslide_record& rec = out.by_daylight[static_cast<std::size_t>(lid - 1)];
        rec.id = lid;
        rec.daylight_flat = static_cast<std::int64_t>(d);
        const int r = static_cast<int>(d / nc), c = static_cast<int>(d % nc);
        rec.daylight_x = meta.x_min + (c + 0.5) * b;
        rec.daylight_y = meta.y_min + (static_cast<double>(nr) - r - 0.5) * b;
        finish(rec, acc_d[static_cast<std::size_t>(lid)], meta);
    }
    // Pieces per by-daylight landslide.
    {
        // Flood fill within one id at a time: cells connect only to the same id.
        std::vector<std::uint8_t> seen(n, 0);
        for (std::size_t s = 0; s < n; ++s)
        {
            if (!slide[s] || seen[s]) continue;
            const std::int32_t lid = out.label[s];
            ++out.by_daylight[static_cast<std::size_t>(lid - 1)].n_pieces;
            std::vector<std::size_t> stack{s};
            seen[s] = 1;
            while (!stack.empty())
            {
                const std::size_t cur = stack.back();
                stack.pop_back();
                const int r = static_cast<int>(cur / nc), c = static_cast<int>(cur % nc);
                for (int dr = -1; dr <= 1; ++dr)
                    for (int dc = -1; dc <= 1; ++dc)
                    {
                        const int rr = r + dr, cc2 = c + dc;
                        if (rr < 0 || cc2 < 0 || rr >= nr || cc2 >= nc) continue;
                        const std::size_t j = static_cast<std::size_t>(rr) * nc + cc2;
                        if (!seen[j] && slide[j] && out.label[j] == lid) { seen[j] = 1; stack.push_back(j); }
                    }
            }
        }
    }
    out.amalgamated.resize(static_cast<std::size_t>(n_cc));
    {
        std::vector<std::vector<std::int32_t>> ids(static_cast<std::size_t>(n_cc) + 1);
        for (std::size_t i = 0; i < n; ++i)
            if (slide[i]) ids[static_cast<std::size_t>(out.label_amalgamated[i])].push_back(out.label[i]);
        for (std::int32_t k = 1; k <= n_cc; ++k)
        {
            landslide_record& rec = out.amalgamated[static_cast<std::size_t>(k - 1)];
            rec.id = k;
            finish(rec, acc_a[static_cast<std::size_t>(k)], meta);
            std::vector<std::int32_t>& v = ids[static_cast<std::size_t>(k)];
            std::sort(v.begin(), v.end());
            rec.n_daylight_points = static_cast<std::int32_t>(std::unique(v.begin(), v.end()) - v.begin());
            rec.n_pieces = 1;
        }
    }
    return out;
}

} // namespace lsdtt
