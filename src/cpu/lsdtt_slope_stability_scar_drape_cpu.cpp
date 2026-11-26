// lsdtt_slope_stability_scar_drape_cpu.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html

#include "lsdtt_slope_stability_scar_drape_cpu.hpp"
#include "lsdtt_slope_stability_dem_change_cpu.hpp"
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

// Least squares z = sum coef_j phi_j(u, v) with phi = 1, u, v (n_basis 3) or
// 1, u, v, u^2, uv, v^2 (n_basis 6), u and v centred and scaled by the caller.
// Normal equations by Gaussian elimination with partial pivoting; false if a
// pivot is below 1e-10 of the largest diagonal (points degenerate for the basis).
void rim_basis(double u, double v, int n_basis, double* phi)
{
    phi[0] = 1.0; phi[1] = u; phi[2] = v;
    if (n_basis == 6) { phi[3] = u * u; phi[4] = u * v; phi[5] = v * v; }
}

bool fit_rim_surface(const std::vector<double>& u, const std::vector<double>& v, const std::vector<double>& z,
                     int n_basis, double* coef)
{
    if (static_cast<int>(u.size()) < n_basis) return false;
    double a[6][7] = {};
    double phi[6];
    for (std::size_t k = 0; k < u.size(); ++k)
    {
        rim_basis(u[k], v[k], n_basis, phi);
        for (int i = 0; i < n_basis; ++i)
        {
            for (int j = 0; j < n_basis; ++j) a[i][j] += phi[i] * phi[j];
            a[i][n_basis] += phi[i] * z[k];
        }
    }
    double dmax = 0.0;
    for (int i = 0; i < n_basis; ++i) dmax = std::max(dmax, a[i][i]);
    for (int col = 0; col < n_basis; ++col)
    {
        int piv = col;
        for (int r = col + 1; r < n_basis; ++r) if (std::abs(a[r][col]) > std::abs(a[piv][col])) piv = r;
        if (!(std::abs(a[piv][col]) > 1e-10 * dmax)) return false;
        if (piv != col) for (int j = 0; j <= n_basis; ++j) std::swap(a[col][j], a[piv][j]);
        for (int r = col + 1; r < n_basis; ++r)
        {
            const double f = a[r][col] / a[col][col];
            for (int j = col; j <= n_basis; ++j) a[r][j] -= f * a[col][j];
        }
    }
    for (int i = n_basis - 1; i >= 0; --i)
    {
        double t = a[i][n_basis];
        for (int j = i + 1; j < n_basis; ++j) t -= a[i][j] * coef[j];
        coef[i] = t / a[i][i];
    }
    return true;
}

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
    const bool fitted_surface = opts.surface != scar_drape_surface_type::membrane;
    const int pad = opts.rim_buffer_pixels + 1 + (fitted_surface ? opts.rim_fit_width_pixels : 0);
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
    const auto dilate8 = [&](std::vector<std::uint8_t>& region, int n_steps) {
    for (int it = 0; it < n_steps; ++it)
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
    };
    dilate8(region, opts.rim_buffer_pixels);

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

    // Compactness of the mapped cells: 4 pi A / P^2, P along cell edges.
    if (rec.n_mapped_cells > 0)
    {
        std::int64_t edges = 0;
        for (int r = 0; r < h; ++r)
            for (int c = 0; c < w; ++c)
            {
                if (!s.mapped[lid(r, c)]) continue;
                for (int d = 0; d < 4; ++d)
                    if (!s.mapped[lid(r + dr4[d], c + dc4[d])]) ++edges;
            }
        const double pi = 3.14159265358979323846;
        rec.compactness = 4.0 * pi * static_cast<double>(rec.n_mapped_cells) /
                          (static_cast<double>(edges) * static_cast<double>(edges));
    }

    // Rim descriptors: relief, and the least-squares plane z = a + b x + c y
    // through the rim cells (coordinates centred for conditioning).
    if (rec.n_rim_cells >= 3)
    {
        std::vector<double> rx, ry, rz;
        double zmin = std::numeric_limits<double>::max(), zmax = -std::numeric_limits<double>::max();
        for (int r = 0; r < h; ++r)
            for (int c = 0; c < w; ++c)
                if (s.state[lid(r, c)] == k_rim)
                {
                    const double z = zd[gid(r, c)];
                    rx.push_back(c * cs); ry.push_back(-r * cs); rz.push_back(z);
                    zmin = std::min(zmin, z); zmax = std::max(zmax, z);
                }
        rec.rim_relief_m = zmax - zmin;
        const double n = static_cast<double>(rx.size());
        double mx = 0, my = 0, mz = 0;
        for (std::size_t k = 0; k < rx.size(); ++k) { mx += rx[k]; my += ry[k]; mz += rz[k]; }
        mx /= n; my /= n; mz /= n;
        double sxx = 0, syy = 0, sxy = 0, sxz = 0, syz = 0;
        for (std::size_t k = 0; k < rx.size(); ++k)
        {
            const double x = rx[k] - mx, y = ry[k] - my, z = rz[k] - mz;
            sxx += x * x; syy += y * y; sxy += x * y; sxz += x * z; syz += y * z;
        }
        const double det = sxx * syy - sxy * sxy;
        if (std::abs(det) > 1e-12 * std::max(1.0, sxx * syy))
        {
            const double b = (sxz * syy - syz * sxy) / det;
            const double cc = (syz * sxx - sxz * sxy) / det;
            double ss = 0;
            for (std::size_t k = 0; k < rx.size(); ++k)
            {
                const double e = (rz[k] - mz) - b * (rx[k] - mx) - cc * (ry[k] - my);
                ss += e * e;
            }
            rec.rim_plane_gradient = std::sqrt(b * b + cc * cc);
            rec.rim_plane_residual_m = std::sqrt(ss / n);
        }
    }

    if (rec.n_free_cells == 0) { rec.status = "no_cells"; return s; }
    if (rec.n_rim_cells == 0) { rec.status = "no_valid_rim"; return s; }

    const double area_cell = cs * cs;
    const auto integrate = [&]() {
        for (int r = 0; r < h; ++r)
            for (int c = 0; c < w; ++c)
            {
                const std::size_t l = lid(r, c);
                if (s.state[l] != k_free) continue;
                const double dz = std::max(0.0, s.z[l] - zd[gid(r, c)]);
                rec.volume_m3 += dz * area_cell;
                if (s.mapped[l])
                {
                    rec.mapped_volume_m3 += dz * area_cell;
                    rec.max_depth_m = std::max(rec.max_depth_m, dz);
                    if (dz < 1.0e-6) ++rec.n_contact_cells;
                }
            }
        if (rec.mapped_area_m2 > 0.0) rec.mean_depth_m = rec.mapped_volume_m3 / rec.mapped_area_m2;
    };

    // Fitted surfaces: least squares through the band of rim_fit_width cells
    // outside the free area, then kept no lower than the ground.
    if (fitted_surface)
    {
        std::vector<std::uint8_t> band = region;
        dilate8(band, opts.rim_fit_width_pixels);
        std::vector<double> bu, bv, bz;
        for (int r = 0; r < h; ++r)
            for (int c = 0; c < w; ++c)
                if (band[lid(r, c)] && !region[lid(r, c)] && on_grid(r, c) && zm[gid(r, c)])
                {
                    bu.push_back(c * cs); bv.push_back(-r * cs); bz.push_back(zd[gid(r, c)]);
                }
        rec.surface_fit_n_cells = static_cast<std::int64_t>(bu.size());
        const int n_basis = opts.surface == scar_drape_surface_type::rim_plane ? 3 : 6;
        double mu = 0, mv = 0, mz = 0;
        for (std::size_t k = 0; k < bu.size(); ++k) { mu += bu[k]; mv += bv[k]; mz += bz[k]; }
        const double nb = std::max<double>(1.0, static_cast<double>(bu.size()));
        mu /= nb; mv /= nb; mz /= nb;
        double scale = 0.0;
        for (std::size_t k = 0; k < bu.size(); ++k) scale = std::max({scale, std::abs(bu[k] - mu), std::abs(bv[k] - mv)});
        if (!(scale > 0.0)) scale = 1.0;
        for (std::size_t k = 0; k < bu.size(); ++k) { bu[k] = (bu[k] - mu) / scale; bv[k] = (bv[k] - mv) / scale; bz[k] -= mz; }
        double coef[6] = {};
        if (!fit_rim_surface(bu, bv, bz, n_basis, coef)) { rec.status = "rim_fit_failed"; return s; }
        double phi[6];
        const auto eval = [&](double u, double v) {
            rim_basis(u, v, n_basis, phi);
            double t = 0.0;
            for (int j = 0; j < n_basis; ++j) t += coef[j] * phi[j];
            return t;
        };
        double ss = 0.0;
        for (std::size_t k = 0; k < bu.size(); ++k) { const double e = bz[k] - eval(bu[k], bv[k]); ss += e * e; }
        rec.surface_fit_residual_m = std::sqrt(ss / nb);
        s.z.assign(region.size(), 0.0);
        for (int r = 0; r < h; ++r)
            for (int c = 0; c < w; ++c)
            {
                const std::size_t l = lid(r, c);
                if (s.state[l] == k_free)
                    s.z[l] = std::max<double>(zd[gid(r, c)], mz + eval((c * cs - mu) / scale, (-r * cs - mv) / scale));
                else if (s.state[l] == k_rim)
                    s.z[l] = zd[gid(r, c)];
            }
        rec.iterations = 0;
        rec.final_change_m = 0.0;
        rec.status = "ok";
        integrate();
        return s;
    }

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
    integrate();
    return s;
}

// splitmix64 of (seed, sample): independent streams for neighbouring samples.
std::uint64_t sample_seed(std::uint64_t seed, std::uint64_t k)
{
    std::uint64_t z = seed + 0x9E3779B97F4A7C15ull * (k + 1);
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ull;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBull;
    return z ^ (z >> 31);
}

std::vector<scar_drape_record> solve_records(const lsdtt_raster_cpu& dem, const std::vector<scar_drape_feature>& features,
                                             const scar_drape_options& opts)
{
    const int nf = static_cast<int>(features.size());
    std::vector<scar_drape_record> records(nf);
    const int n_threads = resolve_thread_count(opts.n_threads);
    #pragma omp parallel for schedule(dynamic, 1) num_threads(n_threads)
    for (int k = 0; k < nf; ++k)
        (void)solve_one(dem, features[k], opts, records[k]);
    return records;
}

} // namespace

scar_drape_surface_type parse_scar_drape_surface_type(const std::string& name)
{
    if (name == "membrane") return scar_drape_surface_type::membrane;
    if (name == "rim_plane") return scar_drape_surface_type::rim_plane;
    if (name == "rim_quadric") return scar_drape_surface_type::rim_quadric;
    throw std::invalid_argument("scar drape: surface must be membrane, rim_plane or rim_quadric, not '" + name + "'");
}

std::string scar_drape_surface_type_name(scar_drape_surface_type t)
{
    switch (t)
    {
        case scar_drape_surface_type::rim_plane: return "rim_plane";
        case scar_drape_surface_type::rim_quadric: return "rim_quadric";
        default: return "membrane";
    }
}

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
    if (o.surface == scar_drape_surface_type::rim_plane && o.rim_fit_width_pixels < 1)
        throw std::invalid_argument("scar drape: rim_plane needs a rim fit width of at least 1 pixel");
    if (o.surface == scar_drape_surface_type::rim_quadric && o.rim_fit_width_pixels < 2)
        throw std::invalid_argument("scar drape: rim_quadric needs a rim fit width of at least 2 pixels "
                                    "(a one-cell ring barely constrains the curvature)");
}

void scar_drape_uncertainty(
    const lsdtt_raster_cpu&                dem,
    const std::vector<scar_drape_feature>& features,
    const scar_drape_options&              opts,
    const scar_drape_error_model&          error,
    int                                    n_samples,
    std::uint64_t                          seed,
    std::vector<scar_drape_record>&        records)
{
    validate_scar_drape_options(opts);
    const lsdtt_raster_meta& meta = dem.get_meta();
    const std::size_t n = static_cast<std::size_t>(meta.n_rows) * meta.n_cols;
    if (n_samples < 2) throw std::invalid_argument("scar_drape_uncertainty: need at least 2 samples");
    if (records.size() != features.size())
        throw std::invalid_argument("scar_drape_uncertainty: records and features differ in length");
    if (!error.sigma_m.empty() && error.sigma_m.size() != n)
        throw std::invalid_argument("scar_drape_uncertainty: sigma_m needs one value per cell");
    if (error.sigma_m.empty() && !(error.sigma_constant_m > 0.0 && std::isfinite(error.sigma_constant_m)))
        throw std::invalid_argument("scar_drape_uncertainty: the constant sigma must be > 0 m");
    if (!(error.nugget_fraction >= 0.0 && error.nugget_fraction <= 1.0) || !(error.correlation_range_m >= 0.0))
        throw std::invalid_argument("scar_drape_uncertainty: nugget must be in [0, 1] and range >= 0");

    const std::size_t nf = features.size();
    std::vector<int> count(nf, 0);
    std::vector<double> mean(nf, 0.0), m2(nf, 0.0);
    const std::vector<float>& z0 = dem.get_data();
    const std::vector<std::uint8_t>& mask = dem.get_mask();
    for (int smp = 0; smp < n_samples; ++smp)
    {
        const std::vector<float> noise = simulate_dem_error_noise(
            meta.n_rows, meta.n_cols, meta.cell_size, error.nugget_fraction, error.correlation_range_m,
            sample_seed(seed, static_cast<std::uint64_t>(smp)));
        std::vector<float> z(z0);
        for (std::size_t i = 0; i < n; ++i)
        {
            if (!mask[i]) continue;
            const double sig = error.sigma_m.empty() ? error.sigma_constant_m : static_cast<double>(error.sigma_m[i]);
            if (std::isfinite(sig) && sig > 0.0) z[i] = static_cast<float>(z[i] + sig * noise[i]);
        }
        const lsdtt_raster_cpu noisy(std::move(z), std::vector<std::uint8_t>(mask), meta);
        const std::vector<scar_drape_record> rs = solve_records(noisy, features, opts);
        for (std::size_t k = 0; k < nf; ++k)
        {
            if (rs[k].status != "ok") continue;
            const double x = rs[k].mapped_volume_m3;
            ++count[k];
            const double d = x - mean[k];
            mean[k] += d / count[k];
            m2[k] += d * (x - mean[k]);
        }
    }
    for (std::size_t k = 0; k < nf; ++k)
    {
        records[k].n_uncertainty_samples = count[k];
        if (count[k] >= 2) records[k].volume_sigma_m3 = std::sqrt(m2[k] / (count[k] - 1));
        if (count[k] >= 1 && records[k].status == "ok") records[k].volume_noise_bias_m3 = mean[k] - records[k].mapped_volume_m3;
    }
}

void scar_drape_buffer_spread(
    const lsdtt_raster_cpu&                dem,
    const std::vector<scar_drape_feature>& features,
    const scar_drape_options&              opts,
    int                                    max_buffer_pixels,
    std::vector<scar_drape_record>&        records)
{
    if (max_buffer_pixels < 1) throw std::invalid_argument("scar_drape_buffer_spread: max buffer must be >= 1");
    if (records.size() != features.size())
        throw std::invalid_argument("scar_drape_buffer_spread: records and features differ in length");
    const std::size_t nf = features.size();
    std::vector<int> n_ok(nf, 0);
    std::vector<double> lo(nf, std::numeric_limits<double>::infinity()), hi(nf, -std::numeric_limits<double>::infinity());
    for (int b = 0; b <= max_buffer_pixels; ++b)
    {
        scar_drape_options o = opts;
        o.rim_buffer_pixels = b;
        validate_scar_drape_options(o);
        const std::vector<scar_drape_record> rs = solve_records(dem, features, o);
        for (std::size_t k = 0; k < nf; ++k)
        {
            if (rs[k].status != "ok") continue;
            ++n_ok[k];
            lo[k] = std::min(lo[k], rs[k].mapped_volume_m3);
            hi[k] = std::max(hi[k], rs[k].mapped_volume_m3);
        }
    }
    for (std::size_t k = 0; k < nf; ++k)
    {
        if (n_ok[k] != max_buffer_pixels + 1) continue;   // every buffer must solve
        records[k].volume_buffer_min_m3 = lo[k];
        records[k].volume_buffer_max_m3 = hi[k];
        records[k].volume_buffer_spread_m3 = hi[k] - lo[k];
    }
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

power_law_scaling_fit fit_scar_volume_area(
    const std::vector<scar_drape_record>& records, int n_bootstrap, std::uint64_t seed)
{
    std::vector<double> a, v;
    for (const auto& r : records)
        if (r.status == "ok") { a.push_back(r.mapped_area_m2); v.push_back(r.mapped_volume_m3); }
    return fit_power_law_scaling(a, v, n_bootstrap, seed);
}

std::vector<scar_reference_volume> scar_reference_volumes(
    const lsdtt_raster_cpu&                post_dem,
    const lsdtt_raster_cpu&                reference_dem,
    const std::vector<scar_drape_feature>& features)
{
    const lsdtt_raster_meta& a = post_dem.get_meta();
    const lsdtt_raster_meta& b = reference_dem.get_meta();
    if (a.n_rows != b.n_rows || a.n_cols != b.n_cols || std::abs(a.x_min - b.x_min) > 1e-6 * a.cell_size ||
        std::abs(a.y_min - b.y_min) > 1e-6 * a.cell_size || std::abs(a.cell_size - b.cell_size) > 1e-9 * a.cell_size)
        throw std::invalid_argument("scar_reference_volumes: the reference DEM is not on the post-event DEM's grid");
    const auto& pz = post_dem.get_data();
    const auto& pm = post_dem.get_mask();
    const auto& rz = reference_dem.get_data();
    const auto& rm = reference_dem.get_mask();
    const double area = a.cell_size * a.cell_size;
    const std::int64_t n = static_cast<std::int64_t>(a.n_rows) * a.n_cols;
    std::vector<scar_reference_volume> out(features.size());
    for (std::size_t k = 0; k < features.size(); ++k)
        for (const std::int64_t i : features[k].cells)
        {
            if (i < 0 || i >= n || !pm[i] || !rm[i]) continue;
            const double d = static_cast<double>(rz[i]) - static_cast<double>(pz[i]);
            ++out[k].n_cells;
            out[k].net_volume_m3 += d * area;
            if (d > 0) out[k].loss_volume_m3 += d * area;
            else out[k].gain_volume_m3 -= d * area;
        }
    return out;
}

scar_drape_validation validate_scar_drape(
    const std::vector<scar_drape_record>&    drape,
    const std::vector<scar_drape_record>&    control,
    const std::vector<scar_reference_volume>& reference,
    double min_reference_volume_m3, int n_bootstrap, std::uint64_t seed)
{
    if (drape.size() != control.size() || drape.size() != reference.size())
        throw std::invalid_argument("validate_scar_drape: drape, control and reference differ in length");
    scar_drape_validation v;
    v.n_scars = drape.size();
    std::vector<double> ratio, la, lv, a_d, v_d, a_r, v_r;
    for (std::size_t k = 0; k < drape.size(); ++k)
    {
        const auto& d = drape[k];
        if (d.status != "ok") continue;
        v.drape_total_m3 += d.mapped_volume_m3;
        if (control[k].status == "ok") v.control_total_m3 += control[k].mapped_volume_m3;
        v.reference_net_total_m3 += reference[k].net_volume_m3;
        v.reference_loss_total_m3 += reference[k].loss_volume_m3;
        a_d.push_back(d.mapped_area_m2); v_d.push_back(d.mapped_volume_m3);
        a_r.push_back(d.mapped_area_m2); v_r.push_back(reference[k].net_volume_m3);
        if (reference[k].net_volume_m3 > min_reference_volume_m3)
        {
            ratio.push_back(d.mapped_volume_m3 / reference[k].net_volume_m3);
            if (d.mapped_volume_m3 > 0)
            {
                la.push_back(std::log10(reference[k].net_volume_m3));
                lv.push_back(std::log10(d.mapped_volume_m3));
            }
        }
    }
    v.n_compared = ratio.size();
    if (v.reference_net_total_m3 > 0)
    {
        v.total_ratio = v.drape_total_m3 / v.reference_net_total_m3;
        v.control_fraction = v.control_total_m3 / v.reference_net_total_m3;
    }
    if (!ratio.empty())
    {
        std::sort(ratio.begin(), ratio.end());
        const auto q = [&](double p) {
            const double pos = p * static_cast<double>(ratio.size() - 1);
            const std::size_t lo = static_cast<std::size_t>(std::floor(pos));
            const std::size_t hi = std::min(lo + 1, ratio.size() - 1);
            return ratio[lo] + (pos - static_cast<double>(lo)) * (ratio[hi] - ratio[lo]);
        };
        v.ratio_p25 = q(0.25); v.median_ratio = q(0.5); v.ratio_p75 = q(0.75);
    }
    if (la.size() >= 3)
    {
        double mx = 0, my = 0;
        for (std::size_t k = 0; k < la.size(); ++k) { mx += la[k]; my += lv[k]; }
        mx /= la.size(); my /= la.size();
        double sxx = 0, syy = 0, sxy = 0;
        for (std::size_t k = 0; k < la.size(); ++k)
        {
            sxx += (la[k] - mx) * (la[k] - mx); syy += (lv[k] - my) * (lv[k] - my); sxy += (la[k] - mx) * (lv[k] - my);
        }
        if (sxx > 0 && syy > 0) v.log_volume_correlation = sxy / std::sqrt(sxx * syy);
    }
    v.drape_fit = fit_power_law_scaling(a_d, v_d, n_bootstrap, seed);
    v.reference_fit = fit_power_law_scaling(a_r, v_r, n_bootstrap, seed + 1);
    return v;
}

} // namespace lsdtt
