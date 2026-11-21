// lsdtt_slope_stability_columns_cpu.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// 3-D ellipsoids by the method of columns. Equations and conventions: see the header.

#include "lsdtt_slope_stability_columns_cpu.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>

#include "../utils/lsdtt_numeric_utils.hpp"
#include "../utils/lsdtt_random.hpp"

#ifdef _OPENMP
#include <omp.h>
#endif

namespace lsdtt
{

namespace
{
constexpr double k_inf = std::numeric_limits<double>::infinity();
constexpr double k_max_gradient = 10.0;

inline double map_x(const lsdtt_raster_meta& m, int col) { return m.x_min + (col + 0.5) * m.cell_size; }
inline double map_y(const lsdtt_raster_meta& m, int row) { return m.y_min + (m.n_rows - row - 0.5) * m.cell_size; }
}

ellipsoid_geometry make_ellipsoid_through(double xR, double yR, double zR, double xD, double yD, double zD,
                                          double c_axis, double b_axis, double z_offset)
{
    ellipsoid_geometry g;
    const double dx = xD - xR, dy = yD - yR, dz = zD - zR;
    const double L = std::sqrt(dx * dx + dy * dy + dz * dz), Lh = std::hypot(dx, dy);
    if (!(L > 0.0) || !(Lh > 0.0) || !(c_axis > 0.0) || !(b_axis > 0.0)) return g;
    g.a = 0.5 * L;
    g.b = b_axis;
    g.c = c_axis;
    g.u[0] = dx / L; g.u[1] = dy / L; g.u[2] = dz / L;
    g.v[0] = -dy / Lh; g.v[1] = dx / Lh; g.v[2] = 0.0;
    // w = u x v, flipped to point up.
    g.w[0] = g.u[1] * g.v[2] - g.u[2] * g.v[1];
    g.w[1] = g.u[2] * g.v[0] - g.u[0] * g.v[2];
    g.w[2] = g.u[0] * g.v[1] - g.u[1] * g.v[0];
    if (g.w[2] < 0.0) for (double& q : g.w) q = -q;
    g.cx = 0.5 * (xR + xD);
    g.cy = 0.5 * (yR + yD);
    g.cz = 0.5 * (zR + zD) + z_offset;
    g.valid = true;
    return g;
}

bool ellipsoid_lower_surface(const ellipsoid_geometry& g, double x, double y, double& z, double& dzdx, double& dzdy)
{
    if (!g.valid) return false;
    const double px = x - g.cx, py = y - g.cy;
    const double* ax[3] = {g.u, g.v, g.w};
    const double s2[3] = {g.a * g.a, g.b * g.b, g.c * g.c};
    double A = 0.0, B = 0.0, C = -1.0;
    for (int k = 0; k < 3; ++k)
    {
        const double hk = px * ax[k][0] + py * ax[k][1];
        A += ax[k][2] * ax[k][2] / s2[k];
        B += 2.0 * hk * ax[k][2] / s2[k];
        C += hk * hk / s2[k];
    }
    const double disc = B * B - 4.0 * A * C;
    if (!(disc > 0.0)) return false;
    const double dzl = (-B - std::sqrt(disc)) / (2.0 * A);
    z = g.cz + dzl;
    double Gx = 0.0, Gy = 0.0, Gz = 0.0;
    for (int k = 0; k < 3; ++k)
    {
        const double P = px * ax[k][0] + py * ax[k][1] + dzl * ax[k][2];
        Gx += 2.0 * P * ax[k][0] / s2[k];
        Gy += 2.0 * P * ax[k][1] / s2[k];
        Gz += 2.0 * P * ax[k][2] / s2[k];
    }
    if (!(Gz < 0.0))
    {
        dzdx = dzdy = 0.0;
        return false;
    }
    dzdx = -Gx / Gz;
    dzdy = -Gy / Gz;
    const double gm = std::hypot(dzdx, dzdy);
    if (gm > k_max_gradient) { dzdx *= k_max_gradient / gm; dzdy *= k_max_gradient / gm; }
    return true;
}

double column_area(double dzdx, double dzdy, double cell_size)
{
    return cell_size * cell_size * std::sqrt(1.0 + dzdx * dzdx + dzdy * dzdy);
}

double apparent_dip_gradient(double dzdx, double dzdy, double bearing_degrees)
{
    const double t = bearing_degrees * M_PI / 180.0;
    // Descent along the unit horizontal vector (sin t east, cos t north).
    return -(dzdx * std::sin(t) + dzdy * std::cos(t));
}

void add_column(column_sums& s, double dzdx, double dzdy, double cell_size, double weight, double u,
                double bearing, column_method method)
{
    const double A = column_area(dzdx, dzdy, cell_size);
    const double tan_c = std::hypot(dzdx, dzdy);
    const double cos_c = 1.0 / std::sqrt(1.0 + tan_c * tan_c);
    // tan beta_m = tan beta_c cos(alpha_c - alpha) = the apparent dip along alpha.
    const double tan_m = apparent_dip_gradient(dzdx, dzdy, bearing);
    const double cos_m = 1.0 / std::sqrt(1.0 + tan_m * tan_m), sin_m = tan_m * cos_m;
    const double f = (method == column_method::xie_revised) ? cos_m : 1.0;
    s.cohesion += f * A;
    s.normal += f * (weight * cos_c - u * A);
    s.driving += f * weight * sin_m;
}

double column_safety_factor(const column_sums& s, double c, double mu)
{
    if (!(s.driving > 0.0)) return k_inf;
    return (c * s.cohesion + mu * s.normal) / s.driving;
}

ellipsoid_evaluation evaluate_ellipsoid(const lsdtt_raster_cpu& dem, const ellipsoid_geometry& g, double bearing,
                                        const rupture_parameters& p, column_method method,
                                        const std::vector<double>& wtd, const std::vector<double>* layer,
                                        bool keep_cells)
{
    ellipsoid_evaluation out;
    if (!g.valid) return out;
    const lsdtt_raster_meta& m = dem.get_meta();
    const int nr = m.n_rows, nc = m.n_cols;
    const double b = m.cell_size;
    const std::vector<float>& z = dem.get_data();
    const std::vector<std::uint8_t>& zm = dem.get_mask();
    const double rho_g = p.density_kg_m3 * p.gravity_m_s2, rw_g = p.water_density_kg_m3 * p.gravity_m_s2;
    const double R = std::max({g.a, g.b, g.c}) + b;
    const int c0 = std::max(0, static_cast<int>(std::floor((g.cx - R - m.x_min) / b)));
    const int c1 = std::min(nc - 1, static_cast<int>(std::ceil((g.cx + R - m.x_min) / b)));
    const int r0 = std::max(0, static_cast<int>(std::floor((m.y_min + nr * b - (g.cy + R)) / b)));
    const int r1 = std::min(nr - 1, static_cast<int>(std::ceil((m.y_min + nr * b - (g.cy - R)) / b)));
    auto layer_base = [&](int r, int c, double& zl) {
        const std::size_t i = static_cast<std::size_t>(r) * nc + c;
        if (!zm[i] || !std::isfinite((*layer)[i])) return false;
        zl = static_cast<double>(z[i]) - (*layer)[i];
        return true;
    };
    for (int r = r0; r <= r1; ++r)
        for (int c = c0; c <= c1; ++c)
        {
            const std::size_t i = static_cast<std::size_t>(r) * nc + c;
            if (!zm[i]) continue;
            double zs, gx, gy;
            if (!ellipsoid_lower_surface(g, map_x(m, c), map_y(m, r), zs, gx, gy)) continue;
            const double zg = static_cast<double>(z[i]);
            if (layer)
            {
                double zl;
                if (layer_base(r, c, zl) && zl > zs)
                {
                    zs = zl;
                    // Layer gradient by central differences (one-sided at edges or nodata).
                    double zl_e = 0.0, zl_w = 0.0, zl_n = 0.0, zl_s = 0.0;
                    const bool e = c + 1 < nc && layer_base(r, c + 1, zl_e), wv = c > 0 && layer_base(r, c - 1, zl_w);
                    const bool nn = r > 0 && layer_base(r - 1, c, zl_n), sv = r + 1 < nr && layer_base(r + 1, c, zl_s);
                    gx = (e && wv) ? (zl_e - zl_w) / (2 * b) : e ? (zl_e - zl) / b : wv ? (zl - zl_w) / b : 0.0;
                    gy = (nn && sv) ? (zl_n - zl_s) / (2 * b) : nn ? (zl_n - zl) / b : sv ? (zl - zl_s) / b : 0.0;
                }
            }
            const double t = zg - zs;
            if (!(t > 0.0)) continue;
            const double W = rho_g * b * b * t;
            double u = 0.0;
            if (!wtd.empty() && std::isfinite(wtd[i]) && rw_g > 0.0) u = rw_g * std::max(0.0, zg - wtd[i] - zs);
            add_column(out.sums, gx, gy, b, W, u, bearing, method);
            ++out.n_columns;
            out.volume_m3 += t * b * b;
            out.max_depth_m = std::max(out.max_depth_m, t);
            if (keep_cells) { out.cells.push_back(static_cast<std::int64_t>(i)); out.depth.push_back(t); }
        }
    out.admissible = out.n_columns > 0 && out.sums.driving > 0.0;
    return out;
}

void validate_ellipsoid_options(const ellipsoid_options& o, std::size_t n_cells)
{
    if (o.search_points < 2) throw std::invalid_argument("ellipsoids: search points must be at least 2");
    if (o.random_search && o.random_draws < 1) throw std::invalid_argument("ellipsoids: random draws must be at least 1");
    if (!(o.min_half_axis_pixels > 0.0)) throw std::invalid_argument("ellipsoids: the minimum half-axis must be > 0 pixels");
    if (std::isnan(o.screen_safety_factor)) throw std::invalid_argument("ellipsoids: the screen safety factor must not be NaN");
    for (const auto& l : o.layer_depths)
        if (l.size() != n_cells) throw std::invalid_argument("ellipsoids: a layer-depth raster differs from the grid");
}

namespace
{

struct crown_geometry
{
    double xR, yR, zR, xD, yD, zD, a, bearing;
};

// F of the ellipsoid (c, b, z_offset) for the crown, lowest over truncations.
double crown_F(const lsdtt_raster_cpu& dem, const crown_geometry& cg, double c_axis, double b_axis, double z_off,
               const rupture_parameters& p, const ellipsoid_options& o, const std::vector<double>& wtd,
               int* layer_out, column_sums* sums_out, int only_layer = -2)
{
    const ellipsoid_geometry g = make_ellipsoid_through(cg.xR, cg.yR, cg.zR, cg.xD, cg.yD, cg.zD, c_axis, b_axis, z_off);
    double best = k_inf;
    for (int L = -1; L < static_cast<int>(o.layer_depths.size()); ++L)
    {
        if (only_layer != -2 && L != only_layer) continue;
        const ellipsoid_evaluation e = evaluate_ellipsoid(dem, g, cg.bearing, p, o.method, wtd,
                                                          L < 0 ? nullptr : &o.layer_depths[static_cast<std::size_t>(L)], false);
        if (!e.admissible) continue;
        const double F = column_safety_factor(e.sums, p.cohesion_pa, p.friction_coefficient);
        if (F < best)
        {
            best = F;
            if (layer_out) *layer_out = L;
            if (sums_out) *sums_out = e.sums;
        }
    }
    return best;
}

} // namespace

ellipsoid_surfaces compute_ellipsoid_surfaces(
    const lsdtt_raster_cpu& dem,
    const std::vector<std::int32_t>& region,
    const std::map<std::int32_t, rupture_parameters>& by_region,
    const rupture_rasters& rupture,
    const ellipsoid_options& o,
    int n_threads,
    const std::vector<double>& wtd)
{
    const lsdtt_raster_meta& meta = dem.get_meta();
    const int nc = meta.n_cols;
    const std::size_t n = static_cast<std::size_t>(meta.n_rows) * static_cast<std::size_t>(nc);
    validate_ellipsoid_options(o, n);
    if (rupture.daylight.size() != n)
        throw std::invalid_argument("compute_ellipsoid_surfaces: the rupture result differs from the grid");
    if (!region.empty() && region.size() != n)
        throw std::invalid_argument("compute_ellipsoid_surfaces: the region grid differs from the grid");
    if (!wtd.empty() && wtd.size() != n)
        throw std::invalid_argument("compute_ellipsoid_surfaces: the water-table grid differs from the grid");
    if (by_region.empty() || (region.empty() && by_region.size() != 1))
        throw std::invalid_argument("compute_ellipsoid_surfaces: one parameter set without regions, or one per region");
    const rupture_parameters* single = region.empty() ? &by_region.begin()->second : nullptr;
    const std::vector<float>& z = dem.get_data();
    const std::vector<float>& rf = rupture.safety_factor.get_data();
    const std::vector<std::uint8_t>& rm = rupture.safety_factor.get_mask();
    const double b = meta.cell_size, b_min = o.min_half_axis_pixels * b;

    ellipsoid_surfaces out;
    out.crown.assign(n, ellipsoid_crown{});
    slip_surface_rasters& sl = out.slip;
    sl.toe.assign(n, -1);
    sl.shape.assign(n, slip_shape{});
    const float nd = static_cast<float>(meta.no_data_value);
    std::vector<float> fs(n, nd), hc(n, nd), hmax(n, nd), bw(n, nd), len(n, nd), uns(n, nd), ev(n, nd);
    std::vector<std::uint8_t> m_eval(n, 0), m_flag(n, 0);
    std::int64_t n_eval = 0, n_screen = 0, n_unst = 0, n_inadm = 0, n_small = 0;

#ifdef _OPENMP
    const int threads = (n_threads > 0) ? n_threads : omp_get_max_threads();
#pragma omp parallel for schedule(dynamic, 8) num_threads(threads) reduction(+ : n_eval, n_screen, n_unst, n_inadm, n_small)
#else
    (void)n_threads;
#endif
    for (std::int64_t ii = 0; ii < static_cast<std::int64_t>(n); ++ii)
    {
        const std::size_t i = static_cast<std::size_t>(ii);
        if (!rm[i] || rupture.daylight[i] < 0) continue;
        if (!single && region[i] == k_rupture_region_nodata) continue;
        m_flag[i] = 1;
        if (!(static_cast<double>(rf[i]) <= o.screen_safety_factor)) { ev[i] = 0.0f; ++n_screen; continue; }
        ev[i] = 1.0f;
        ++n_eval;
        const rupture_parameters& p = single ? *single : by_region.at(region[i]);
        const std::int64_t D = rupture.daylight[i];
        crown_geometry cg;
        cg.xR = map_x(meta, static_cast<int>(ii % nc)); cg.yR = map_y(meta, static_cast<int>(ii / nc));
        cg.zR = static_cast<double>(z[i]);
        cg.xD = map_x(meta, static_cast<int>(D % nc)); cg.yD = map_y(meta, static_cast<int>(D / nc));
        cg.zD = static_cast<double>(z[static_cast<std::size_t>(D)]);
        cg.a = 0.5 * std::sqrt((cg.xD - cg.xR) * (cg.xD - cg.xR) + (cg.yD - cg.yR) * (cg.yD - cg.yR)
                               + (cg.zD - cg.zR) * (cg.zD - cg.zR));
        cg.bearing = grid_bearing_degrees(ii, D, nc);
        if (!(cg.a >= b_min)) { ++n_small; ++n_inadm; continue; }
        const double b_hi = std::max(b_min, 2.0 * cg.a);
        // Each option (untruncated, then each layer) is searched on its own, and
        // the best kept: never above the untruncated search.
        double best_f = k_inf, best_cf = 0.0, best_b = b_min, best_zo = 0.0;
        int best_layer = -1;
        const int N = o.search_points;
        for (int L = -1; L < static_cast<int>(o.layer_depths.size()); ++L)
        {
            auto F2 = [&](double cf, double bb, double zo) {
                return crown_F(dem, cg, cf * cg.a, bb, zo, p, o, wtd, nullptr, nullptr, L);
            };
            double lf = k_inf, lcf = 0.0, lb_ = b_min, lzo = 0.0;
            if (o.random_search)
            {
                lsdtt::random::splitmix64_stream rng(o.seed * 0x9E3779B97F4A7C15ULL + static_cast<std::uint64_t>(ii));
                for (int k = 0; k < o.random_draws; ++k)
                {
                    const double cf = 1e-3 + (1.0 - 1e-3) * rng.uniform01();
                    const double bb = b_min + (b_hi - b_min) * rng.uniform01();
                    const double zo = (rng.uniform01() - 0.5) * cf * cg.a;
                    const double v = F2(cf, bb, zo);
                    if (v < lf) { lf = v; lcf = cf; lb_ = bb; lzo = zo; }
                }
            }
            else
            {
                auto b_of = [&](int j) { return (N > 1) ? b_min * std::pow(b_hi / b_min, static_cast<double>(j) / (N - 1)) : b_min; };
                for (int ic = 0; ic < N; ++ic)
                    for (int jb = 0; jb < N; ++jb)
                    {
                        const double cf = (ic + 0.5) / N, bb = b_of(jb);
                        const double v = F2(cf, bb, 0.0);
                        if (v < lf) { lf = v; lcf = cf; lb_ = bb; }
                    }
                if (std::isfinite(lf))
                {
                    double dc = 1.0 / N, lb = std::log(b_hi / b_min) / std::max(1, N - 1);
                    for (int round = 0; round < 2; ++round)
                    {
                        const double bb0 = lb_;
                        const brent_minimum rc = minimise_brent_bounded([&](double cf) { return F2(cf, bb0, 0.0); },
                            std::max(1e-3, lcf - dc), std::min(1.0, lcf + dc), 1e-5, 60);
                        if (rc.f < lf) { lf = rc.f; lcf = rc.x; }
                        const double cf0 = lcf;
                        const double lo = std::max(b_min, lb_ * std::exp(-lb)), hi = std::min(b_hi, lb_ * std::exp(lb));
                        const brent_minimum rb = minimise_brent_bounded([&](double bb) { return F2(cf0, bb, 0.0); },
                            lo, std::max(lo, hi), 1e-5 * b_hi, 60);
                        if (rb.f < lf) { lf = rb.f; lb_ = rb.x; }
                        dc *= 0.5;
                        lb *= 0.5;
                    }
                }
            }
            if (lf < best_f) { best_f = lf; best_cf = lcf; best_b = lb_; best_zo = lzo; best_layer = L; }
        }
        if (!(best_f < k_inf)) { ++n_inadm; continue; }
        ellipsoid_crown& cr = out.crown[i];
        cr.valid = true;
        cr.c_axis = best_cf * cg.a;
        cr.b_axis = best_b;
        cr.z_offset = best_zo;
        cr.bearing = cg.bearing;
        cr.layer = best_layer;
        (void)crown_F(dem, cg, cr.c_axis, cr.b_axis, cr.z_offset, p, o, wtd, nullptr, &cr.sums, best_layer);
        const ellipsoid_geometry g = make_ellipsoid_through(cg.xR, cg.yR, cg.zR, cg.xD, cg.yD, cg.zD, cr.c_axis, cr.b_axis, cr.z_offset);
        const ellipsoid_evaluation e = evaluate_ellipsoid(dem, g, cg.bearing, p, o.method, wtd,
            cr.layer < 0 ? nullptr : &o.layer_depths[static_cast<std::size_t>(cr.layer)], false);
        m_eval[i] = 1;
        fs[i] = static_cast<float>(best_f);
        hc[i] = static_cast<float>(cr.c_axis);
        hmax[i] = static_cast<float>(e.max_depth_m);
        bw[i] = static_cast<float>(cr.b_axis);
        len[i] = static_cast<float>(std::hypot(cg.xD - cg.xR, cg.yD - cg.yR));
        uns[i] = (best_f <= 1.0) ? 1.0f : 0.0f;
        if (best_f <= 1.0) ++n_unst;
        sl.toe[i] = D;
    }

    // Bodies (unstable) and the cover F (all evaluated), serially in crown order.
    std::vector<double> env(n, 0.0), cover(n, k_inf);
    sl.body_offset.push_back(0);
    for (std::size_t i = 0; i < n; ++i)
    {
        if (!m_eval[i]) continue;
        const ellipsoid_crown& cr = out.crown[i];
        const rupture_parameters& p = single ? *single : by_region.at(region[i]);
        const std::int64_t D = sl.toe[i];
        const ellipsoid_geometry g = make_ellipsoid_through(
            map_x(meta, static_cast<int>(i % nc)), map_y(meta, static_cast<int>(i / nc)), z[i],
            map_x(meta, static_cast<int>(D % nc)), map_y(meta, static_cast<int>(D / nc)), z[static_cast<std::size_t>(D)],
            cr.c_axis, cr.b_axis, cr.z_offset);
        const ellipsoid_evaluation e = evaluate_ellipsoid(dem, g, cr.bearing, p, o.method, wtd,
            cr.layer < 0 ? nullptr : &o.layer_depths[static_cast<std::size_t>(cr.layer)], true);
        const double F = fs[i];
        for (std::int64_t c : e.cells) cover[static_cast<std::size_t>(c)] = std::min(cover[static_cast<std::size_t>(c)], F);
        if (!(F <= 1.0)) continue;
        sl.body_crown.push_back(static_cast<std::int64_t>(i));
        for (std::size_t q = 0; q < e.cells.size(); ++q)
        {
            sl.body_cell.push_back(e.cells[q]);
            sl.body_cell_depth.push_back(e.depth[q]);
            env[static_cast<std::size_t>(e.cells[q])] = std::max(env[static_cast<std::size_t>(e.cells[q])], e.depth[q]);
        }
        sl.body_offset.push_back(static_cast<std::int64_t>(sl.body_cell.size()));
    }
    std::vector<float> bd(n, nd), cv(n, nd);
    std::vector<std::uint8_t> bm(n, 0), cm(n, 0);
    for (std::size_t i = 0; i < n; ++i)
    {
        if (env[i] > 0.0) { bd[i] = static_cast<float>(env[i]); bm[i] = 1; }
        if (cover[i] < k_inf) { cv[i] = static_cast<float>(cover[i]); cm[i] = 1; }
    }
    auto make = [&](std::vector<float>&& v, const std::vector<std::uint8_t>& mk) {
        return lsdtt_raster_cpu(std::move(v), std::vector<std::uint8_t>(mk), meta);
    };
    const std::vector<std::uint8_t> none(n, 0);
    sl.safety_factor = make(std::move(fs), m_eval);
    sl.crown_depth = make(std::move(hc), m_eval);
    sl.max_depth = make(std::move(hmax), m_eval);
    sl.sagitta = make(std::move(bw), m_eval);
    sl.basal_gradient = make(std::vector<float>(n, nd), none);
    sl.length = make(std::move(len), m_eval);
    sl.sinuosity = make(std::vector<float>(n, nd), none);
    sl.unstable = make(std::move(uns), m_eval);
    sl.evaluated = make(std::move(ev), m_flag);
    sl.body_depth = make(std::move(bd), bm);
    sl.cover_safety_factor = make(std::move(cv), cm);
    sl.n_evaluated = n_eval;
    sl.n_screened = n_screen;
    sl.n_unstable = n_unst;
    sl.n_inadmissible = n_inadm;
    out.n_too_small = n_small;
    return out;
}

ellipsoid_probability compute_ellipsoid_failure_probability(
    const ellipsoid_surfaces& s,
    const lsdtt_raster_cpu& dem,
    const std::vector<std::int32_t>& region,
    const std::map<std::int32_t, rupture_parameters>& by_region,
    const std::map<std::int32_t, std::vector<rupture_parameter_sample>>& samples_by_region,
    const ellipsoid_options& o)
{
    const std::vector<double> wtd;   // footprints do not depend on the water table
    const lsdtt_raster_meta& meta = dem.get_meta();
    const int nc = meta.n_cols;
    const std::size_t n = static_cast<std::size_t>(meta.n_rows) * static_cast<std::size_t>(nc);
    if (s.crown.size() != n) throw std::invalid_argument("compute_ellipsoid_failure_probability: the result differs from the grid");
    if (samples_by_region.empty() || (region.empty() && samples_by_region.size() != 1))
        throw std::invalid_argument("compute_ellipsoid_failure_probability: one sample list without regions, or one per region");
    std::size_t max_samples = 0;
    for (const auto& kv : samples_by_region)
    {
        double w = 0.0;
        for (const auto& q : kv.second) w += q.weight;
        if (kv.second.empty() || std::fabs(w - 1.0) > 1e-9)
            throw std::invalid_argument("compute_ellipsoid_failure_probability: empty sample list or weights not summing to 1 (region "
                                        + std::to_string(kv.first) + ")");
        max_samples = std::max(max_samples, kv.second.size());
    }
    const std::vector<rupture_parameter_sample>* single = region.empty() ? &samples_by_region.begin()->second : nullptr;
    const rupture_parameters* single_p = region.empty() ? &by_region.begin()->second : nullptr;
    const std::vector<float>& ev = s.slip.evaluated.get_data();
    const std::vector<std::uint8_t>& em = s.slip.evaluated.get_mask();
    const std::vector<float>& z = dem.get_data();
    const float nd = static_cast<float>(meta.no_data_value);
    std::vector<float> pf(n, nd), mf(n, nd), sf(n, nd);
    std::vector<std::uint8_t> m(n, 0), m_stats(n, 0);
    std::vector<double> cover(n, -1.0);
    std::int64_t n_src = 0, n_cert = 0, n_poss = 0;
    for (std::size_t i = 0; i < n; ++i)
    {
        if (!em[i]) continue;
        if (ev[i] == 0.0f) { pf[i] = 0.0f; m[i] = 1; continue; }
        const ellipsoid_crown& cr = s.crown[i];
        if (!cr.valid) continue;
        const auto it = single ? samples_by_region.end() : samples_by_region.find(region[i]);
        if (!single && it == samples_by_region.end()) continue;
        const std::vector<rupture_parameter_sample>& smp = single ? *single : it->second;
        double p_fail = 0.0, mean = 0.0;
        std::vector<double> F(smp.size());
        for (std::size_t q = 0; q < smp.size(); ++q)
        {
            F[q] = column_safety_factor(cr.sums, smp[q].parameters.cohesion_pa, smp[q].parameters.friction_coefficient);
            if (F[q] <= 1.0) p_fail += smp[q].weight;
            mean += smp[q].weight * F[q];
        }
        double var = 0.0;
        for (std::size_t q = 0; q < smp.size(); ++q) var += smp[q].weight * (F[q] - mean) * (F[q] - mean);
        ++n_src;
        if (p_fail >= 1.0 - 1e-12) ++n_cert;
        if (p_fail > 0.0) ++n_poss;
        m[i] = 1;
        m_stats[i] = 1;
        pf[i] = static_cast<float>(std::min(1.0, p_fail));
        mf[i] = static_cast<float>(mean);
        sf[i] = static_cast<float>(std::sqrt(var));
        // Cover: the footprint of this fixed ellipsoid.
        const rupture_parameters& p = single_p ? *single_p : by_region.at(region[i]);
        const std::int64_t D = s.slip.toe[i];
        const ellipsoid_geometry g = make_ellipsoid_through(
            map_x(meta, static_cast<int>(i % nc)), map_y(meta, static_cast<int>(i / nc)), z[i],
            map_x(meta, static_cast<int>(D % nc)), map_y(meta, static_cast<int>(D / nc)), z[static_cast<std::size_t>(D)],
            cr.c_axis, cr.b_axis, cr.z_offset);
        const ellipsoid_evaluation e = evaluate_ellipsoid(dem, g, cr.bearing, p, o.method, wtd,
            cr.layer < 0 ? nullptr : &o.layer_depths[static_cast<std::size_t>(cr.layer)], true);
        for (std::int64_t c : e.cells) cover[static_cast<std::size_t>(c)] = std::max(cover[static_cast<std::size_t>(c)], p_fail);
    }
    std::vector<float> cv(n, nd);
    std::vector<std::uint8_t> cm(n, 0);
    for (std::size_t i = 0; i < n; ++i)
        if (cover[i] >= 0.0) { cv[i] = static_cast<float>(std::min(1.0, cover[i])); cm[i] = 1; }
    ellipsoid_probability out;
    out.crown.failure_probability = lsdtt_raster_cpu(std::move(pf), std::move(m), meta);
    out.crown.mean_safety_factor = lsdtt_raster_cpu(std::move(mf), std::vector<std::uint8_t>(m_stats), meta);
    out.crown.sd_safety_factor = lsdtt_raster_cpu(std::move(sf), std::move(m_stats), meta);
    out.crown.n_sources = n_src;
    out.crown.n_certain_failure = n_cert;
    out.crown.n_possible_failure = n_poss;
    out.crown.max_samples = max_samples;
    out.cover_failure_probability = lsdtt_raster_cpu(std::move(cv), std::move(cm), meta);
    return out;
}

} // namespace lsdtt
