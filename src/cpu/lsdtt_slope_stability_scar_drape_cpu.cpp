// lsdtt_slope_stability_scar_drape_cpu.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html

#include "lsdtt_slope_stability_scar_drape_cpu.hpp"
#include "../utils/lsdtt_thread_utils.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace lsdtt {
namespace {

// Cell states in a feature's local window.
constexpr std::uint8_t k_outside = 0;  // not part of this feature's problem
constexpr std::uint8_t k_free    = 1;  // cloth free (mapped or buffer, valid DEM)
constexpr std::uint8_t k_rim     = 2;  // cloth held at the DEM
constexpr std::uint8_t k_invalid = 3;  // nodata or off the grid: no-flux

struct local_solution
{
    int r0 = 0, c0 = 0, h = 0, w = 0;     // window in grid rows / columns
    std::vector<std::uint8_t> state;      // per window cell
    std::vector<std::uint8_t> mapped;     // 1 where the polygon covers a valid cell
    std::vector<double>       z;          // cloth (free cells)
};

local_solution solve_one(const lsdtt_raster_cpu& dem, const scar_drape_feature& f,
                         const scar_drape_options& opts, scar_drape_record& rec)
{
    const lsdtt_raster_meta& meta = dem.get_meta();
    const int nr = meta.n_rows, nc = meta.n_cols;
    const double cs = meta.cell_size;
    const std::vector<float>& zd = dem.get_data();
    const std::vector<std::uint8_t>& zm = dem.get_mask();

    local_solution s;
    rec.id = f.id;
    rec.feature_index = f.feature_index;
    rec.used_all_touched = f.used_all_touched;
    if (f.cells.empty()) { rec.status = "no_cells"; return s; }

    int rmin = nr, rmax = -1, cmin = nc, cmax = -1;
    for (const std::int64_t i : f.cells)
    {
        const int r = static_cast<int>(i / nc), c = static_cast<int>(i % nc);
        rmin = std::min(rmin, r); rmax = std::max(rmax, r);
        cmin = std::min(cmin, c); cmax = std::max(cmax, c);
    }
    const int pad = opts.rim_buffer_pixels + 1;
    s.r0 = rmin - pad; s.c0 = cmin - pad;
    s.h = rmax - rmin + 1 + 2 * pad; s.w = cmax - cmin + 1 + 2 * pad;
    const int h = s.h, w = s.w;
    const auto lid = [w](int r, int c) { return static_cast<std::size_t>(r) * w + c; };
    const auto on_grid = [&](int lr, int lc) {
        const int r = s.r0 + lr, c = s.c0 + lc;
        return r >= 0 && r < nr && c >= 0 && c < nc;
    };
    const auto gid = [&](int lr, int lc) {
        return static_cast<std::size_t>(s.r0 + lr) * nc + (s.c0 + lc);
    };

    // Mapped cells, then the dilation that defines the free cloth.
    std::vector<std::uint8_t> region(static_cast<std::size_t>(h) * w, 0u);
    for (const std::int64_t i : f.cells)
        region[lid(static_cast<int>(i / nc) - s.r0, static_cast<int>(i % nc) - s.c0)] = 1u;
    s.mapped.assign(region.size(), 0u);
    for (int it = 0; it < opts.rim_buffer_pixels; ++it)
    {
        std::vector<std::uint8_t> next = region;
        for (int r = 0; r < h; ++r)
            for (int c = 0; c < w; ++c)
            {
                if (region[lid(r, c)]) continue;
                for (int dr = -1; dr <= 1 && !next[lid(r, c)]; ++dr)
                    for (int dc = -1; dc <= 1; ++dc)
                    {
                        const int rr = r + dr, cc = c + dc;
                        if (rr >= 0 && rr < h && cc >= 0 && cc < w && region[lid(rr, cc)]) { next[lid(r, c)] = 1u; break; }
                    }
            }
        region.swap(next);
    }

    s.state.assign(region.size(), k_outside);
    for (const std::int64_t i : f.cells)
    {
        const std::size_t l = lid(static_cast<int>(i / nc) - s.r0, static_cast<int>(i % nc) - s.c0);
        if (zm[static_cast<std::size_t>(i)]) s.mapped[l] = 1u;
    }
    for (int r = 0; r < h; ++r)
        for (int c = 0; c < w; ++c)
        {
            if (!region[lid(r, c)]) continue;
            if (!on_grid(r, c)) { s.state[lid(r, c)] = k_invalid; rec.touches_grid_edge = true; continue; }
            if (!zm[gid(r, c)]) { s.state[lid(r, c)] = k_invalid; ++rec.n_nodata_inside; continue; }
            s.state[lid(r, c)] = k_free;
        }
    // The held ring: 4-neighbours of free cells that are not in the region.
    static const int dr4[4] = {-1, 1, 0, 0}, dc4[4] = {0, 0, -1, 1};
    for (int r = 0; r < h; ++r)
        for (int c = 0; c < w; ++c)
        {
            if (s.state[lid(r, c)] != k_free) continue;
            for (int d = 0; d < 4; ++d)
            {
                const int rr = r + dr4[d], cc = c + dc4[d];
                if (region[lid(rr, cc)] || s.state[lid(rr, cc)] != k_outside) continue;
                if (!on_grid(rr, cc)) { s.state[lid(rr, cc)] = k_invalid; ++rec.n_rim_invalid; rec.touches_grid_edge = true; }
                else if (!zm[gid(rr, cc)]) { s.state[lid(rr, cc)] = k_invalid; ++rec.n_rim_invalid; }
                else { s.state[lid(rr, cc)] = k_rim; ++rec.n_rim_cells; }
            }
        }

    double sx = 0.0, sy = 0.0;
    double rim_sum = 0.0;
    for (int r = 0; r < h; ++r)
        for (int c = 0; c < w; ++c)
        {
            const std::uint8_t st = s.state[lid(r, c)];
            if (st == k_free)
            {
                ++rec.n_free_cells;
                if (s.mapped[lid(r, c)])
                {
                    ++rec.n_mapped_cells;
                    sx += meta.x_min + (s.c0 + c + 0.5) * cs;
                    sy += meta.y_min + (nr - (s.r0 + r) - 0.5) * cs;
                }
            }
            else if (st == k_rim) rim_sum += zd[gid(r, c)];
        }
    const double area = cs * cs;
    rec.mapped_area_m2 = rec.n_mapped_cells * area;
    rec.free_area_m2 = rec.n_free_cells * area;
    if (rec.n_mapped_cells > 0) { rec.centroid_x = sx / rec.n_mapped_cells; rec.centroid_y = sy / rec.n_mapped_cells; }
    if (rec.n_free_cells == 0) { rec.status = "no_cells"; return s; }
    if (rec.n_rim_cells == 0) { rec.status = "no_valid_rim"; return s; }

    // Start the cloth above the ground (any start converges to the same answer).
    const double z_start = rim_sum / rec.n_rim_cells;
    s.z.assign(region.size(), 0.0);
    for (int r = 0; r < h; ++r)
        for (int c = 0; c < w; ++c)
        {
            const std::uint8_t st = s.state[lid(r, c)];
            if (st == k_free) s.z[lid(r, c)] = std::max<double>(zd[gid(r, c)], z_start);
            else if (st == k_rim) s.z[lid(r, c)] = zd[gid(r, c)];
        }

    // Projected SOR. Neighbours that are off the problem (k_invalid) carry no
    // flux; the 5-point update with k valid neighbours is the mean of those
    // neighbours minus the sag term.
    const double sag = opts.sag_curvature_per_metre * cs * cs / 4.0;
    const int n_side = std::max(h, w);
    const double pi = 3.14159265358979323846;
    const double omega = 2.0 / (1.0 + std::sin(pi / std::max(2, n_side)));
    double change = 0.0;
    int it = 0;
    for (it = 0; it < opts.max_iterations; ++it)
    {
        change = 0.0;
        for (int r = 0; r < h; ++r)
            for (int c = 0; c < w; ++c)
            {
                const std::size_t l = lid(r, c);
                if (s.state[l] != k_free) continue;
                double sum = 0.0;
                int k = 0;
                for (int d = 0; d < 4; ++d)
                {
                    const std::uint8_t st = s.state[lid(r + dr4[d], c + dc4[d])];
                    if (st == k_free || st == k_rim) { sum += s.z[lid(r + dr4[d], c + dc4[d])]; ++k; }
                }
                if (k == 0) continue;
                const double target = sum / k - sag;
                const double ground = zd[gid(r, c)];
                const double znew = std::max(ground, s.z[l] + omega * (target - s.z[l]));
                change = std::max(change, std::abs(znew - s.z[l]));
                s.z[l] = znew;
            }
        if (change < opts.tolerance_metres) { ++it; break; }
    }
    rec.iterations = it;
    rec.final_change_m = change;
    rec.status = change < opts.tolerance_metres ? "ok" : "not_converged";

    for (int r = 0; r < h; ++r)
        for (int c = 0; c < w; ++c)
        {
            const std::size_t l = lid(r, c);
            if (s.state[l] != k_free) continue;
            const double dz = std::max(0.0, s.z[l] - zd[gid(r, c)]);
            rec.volume_m3 += dz * area;
            if (s.mapped[l])
            {
                rec.mapped_volume_m3 += dz * area;
                rec.max_depth_m = std::max(rec.max_depth_m, dz);
                if (dz < 1.0e-6) ++rec.n_contact_cells;
            }
        }
    if (rec.mapped_area_m2 > 0.0) rec.mean_depth_m = rec.mapped_volume_m3 / rec.mapped_area_m2;
    return s;
}

} // namespace

void validate_scar_drape_options(const scar_drape_options& o)
{
    if (o.rim_buffer_pixels < 0)
        throw std::invalid_argument("scar drape: rim buffer must be >= 0 pixels");
    if (!(o.sag_curvature_per_metre >= 0.0) || !std::isfinite(o.sag_curvature_per_metre))
        throw std::invalid_argument("scar drape: sag curvature must be >= 0 (1/m)");
    if (!(o.tolerance_metres > 0.0))
        throw std::invalid_argument("scar drape: tolerance must be > 0 m");
    if (o.max_iterations < 1)
        throw std::invalid_argument("scar drape: max iterations must be >= 1");
}

scar_drape_result drape_landslide_scars(
    const lsdtt_raster_cpu&                dem,
    const std::vector<scar_drape_feature>& features,
    const scar_drape_options&              opts)
{
    validate_scar_drape_options(opts);
    const lsdtt_raster_meta& meta = dem.get_meta();
    const std::size_t n = static_cast<std::size_t>(meta.n_rows) * meta.n_cols;
    const float nodata = static_cast<float>(meta.no_data_value);

    const int nf = static_cast<int>(features.size());
    std::vector<scar_drape_record> records(nf);
    std::vector<local_solution> sols(nf);
    const int n_threads = resolve_thread_count(opts.n_threads);
    #pragma omp parallel for schedule(dynamic, 1) num_threads(n_threads)
    for (int k = 0; k < nf; ++k)
        sols[k] = solve_one(dem, features[k], opts, records[k]);

    // Merge: where scars overlap keep the deeper cloth.
    std::vector<float> surf(n, nodata), depth(n, nodata), label(n, 0.0f);
    std::vector<std::uint8_t> smask(n, 0u), dmask(n, 0u), lmask(dem.get_mask());
    std::vector<std::uint8_t> count(n, 0u);
    const std::vector<float>& zd = dem.get_data();
    for (int k = 0; k < nf; ++k)
    {
        const local_solution& s = sols[k];
        if (s.z.empty()) continue;
        for (int r = 0; r < s.h; ++r)
            for (int c = 0; c < s.w; ++c)
            {
                const std::size_t l = static_cast<std::size_t>(r) * s.w + c;
                if (s.state[l] != 1u) continue;
                const std::size_t g = static_cast<std::size_t>(s.r0 + r) * meta.n_cols + (s.c0 + c);
                const float dz = static_cast<float>(std::max(0.0, s.z[l] - zd[g]));
                if (count[g] < 255u) ++count[g];
                if (!dmask[g] || dz > depth[g])
                {
                    depth[g] = dz; dmask[g] = 1u;
                    surf[g] = static_cast<float>(s.z[l]); smask[g] = 1u;
                    label[g] = static_cast<float>(features[k].feature_index + 1);
                }
            }
    }
    scar_drape_result out;
    for (std::size_t i = 0; i < n; ++i) if (count[i] > 1u) ++out.n_overlap_cells;
    out.records = std::move(records);
    out.surface = lsdtt_raster_cpu(std::move(surf), std::move(smask), meta);
    out.depth = lsdtt_raster_cpu(std::move(depth), std::move(dmask), meta);
    lsdtt_raster_meta lm = meta;
    out.label = lsdtt_raster_cpu(std::move(label), std::move(lmask), lm);
    return out;
}

scar_volume_area_fit fit_scar_volume_area(const std::vector<scar_drape_record>& records)
{
    std::vector<double> x, y;
    for (const auto& r : records)
        if (r.status == "ok" && r.mapped_volume_m3 > 0.0 && r.mapped_area_m2 > 0.0)
        {
            x.push_back(std::log10(r.mapped_area_m2));
            y.push_back(std::log10(r.mapped_volume_m3));
        }
    scar_volume_area_fit f;
    f.n = static_cast<int>(x.size());
    if (f.n < 3) return f;
    double mx = 0, my = 0;
    for (int i = 0; i < f.n; ++i) { mx += x[i]; my += y[i]; }
    mx /= f.n; my /= f.n;
    double sxx = 0, sxy = 0, syy = 0;
    for (int i = 0; i < f.n; ++i)
    {
        sxx += (x[i] - mx) * (x[i] - mx);
        sxy += (x[i] - mx) * (y[i] - my);
        syy += (y[i] - my) * (y[i] - my);
    }
    if (!(sxx > 0.0)) return f;
    f.gamma = sxy / sxx;
    f.alpha = std::pow(10.0, my - f.gamma * mx);
    double sse = 0;
    for (int i = 0; i < f.n; ++i)
    {
        const double e = y[i] - (my + f.gamma * (x[i] - mx));
        sse += e * e;
    }
    f.r_squared = syy > 0.0 ? 1.0 - sse / syy : 0.0;
    f.gamma_standard_error = std::sqrt(sse / (f.n - 2) / sxx);
    return f;
}

} // namespace lsdtt
