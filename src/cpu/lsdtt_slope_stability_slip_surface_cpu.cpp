// lsdtt_slope_stability_slip_surface_cpu.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// Slip surfaces by the method of slices. Equations and conventions: see the header.

#include "lsdtt_slope_stability_slip_surface_cpu.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>

#include "../utils/lsdtt_numeric_utils.hpp"
#include "lsdtt_slope_stability_catalogue_cpu.hpp"

#ifdef _OPENMP
#include <omp.h>
#endif

namespace lsdtt
{

namespace
{
constexpr double k_inf = std::numeric_limits<double>::infinity();
}

slip_section build_slip_section(const lsdtt_flow_base& flow, const lsdtt_raster_cpu& geometry,
                                const std::vector<std::uint8_t>& channel, std::int64_t crown,
                                const std::vector<double>& water_table_depth)
{
    const lsdtt_raster_meta& meta = geometry.get_meta();
    const int nc = meta.n_cols;
    const std::vector<float>& z = geometry.get_data();
    const std::vector<std::uint8_t>& zm = geometry.get_mask();
    const double b = meta.cell_size;
    const std::int64_t max_steps = static_cast<std::int64_t>(z.size());
    slip_section s;
    s.cell_size = b;
    auto push = [&](std::int64_t flat, double x) {
        const std::size_t i = static_cast<std::size_t>(flat);
        s.x.push_back(x);
        s.z.push_back(static_cast<double>(z[i]));
        s.water_table_depth.push_back(water_table_depth.empty() ? std::numeric_limits<double>::quiet_NaN()
                                                                : water_table_depth[i]);
        s.cell.push_back(flat);
    };
    push(crown, 0.0);
    std::int64_t prev = crown, p = flow.receiver_of_unchecked(crown);
    double x = 0.0;
    for (std::int64_t step = 0; p >= 0 && step < max_steps; ++step)
    {
        const std::size_t pi = static_cast<std::size_t>(p);
        if (channel[pi] || !zm[pi]) break;
        x += b * std::hypot(static_cast<double>(p / nc - prev / nc), static_cast<double>(p % nc - prev % nc));
        push(p, x);
        prev = p;
        p = flow.receiver_of_unchecked(p);
    }
    return s;
}

double slip_section_sinuosity(const slip_section& s, int k, int nc)
{
    if (k <= 0 || k >= static_cast<int>(s.x.size())) return std::numeric_limits<double>::quiet_NaN();
    const std::int64_t a = s.cell[0], d = s.cell[static_cast<std::size_t>(k)];
    const double straight = s.cell_size * std::hypot(static_cast<double>(d / nc - a / nc),
                                                     static_cast<double>(d % nc - a % nc));
    return s.x[static_cast<std::size_t>(k)] / straight;
}

double slip_surface_geometry::base(double x) const
{
    if (straight) return z_crown_base + (z_toe - z_crown_base) * x / x_toe;
    const double dx = x - cx;
    return cz - std::sqrt(std::max(0.0, r * r - dx * dx));
}

slip_surface_geometry make_slip_surface(double x_toe, double z_crown_base, double z_toe, double sagitta)
{
    slip_surface_geometry g;
    g.x_toe = x_toe;
    g.z_crown_base = z_crown_base;
    g.z_toe = z_toe;
    if (!(x_toe > 0.0) || !std::isfinite(z_crown_base) || !std::isfinite(z_toe) || !(sagitta >= 0.0)) return g;
    const double dzv = z_toe - z_crown_base;
    const double L = std::hypot(x_toe, dzv);
    if (sagitta <= 1e-12 * L)
    {
        g.admissible = true;
        g.straight = true;
        return g;
    }
    // Down normal to the chord, centre r - s from the chord on the other side.
    const double nx = dzv / L, nz = -x_toe / L;
    const double r = (0.25 * L * L + sagitta * sagitta) / (2.0 * sagitta);
    const double mx = 0.5 * x_toe, mz = 0.5 * (z_crown_base + z_toe);
    g.straight = false;
    g.r = r;
    g.cx = mx - (r - sagitta) * nx;
    g.cz = mz - (r - sagitta) * nz;
    // Both ends on the lower half-circle, or the surface turns past vertical.
    g.admissible = z_crown_base <= g.cz && z_toe <= g.cz;
    return g;
}


double grid_bearing_degrees(std::int64_t from, std::int64_t to, int nc)
{
    const double de = static_cast<double>(to % nc - from % nc);
    const double dn = -static_cast<double>(to / nc - from / nc);   // rows grow southwards
    if (de == 0.0 && dn == 0.0) return std::numeric_limits<double>::quiet_NaN();
    double b = std::atan2(de, dn) * 180.0 / M_PI;
    if (b < 0.0) b += 360.0;
    return b;
}

double rankine_back_gradient(double mu)
{
    return std::tan(0.25 * M_PI + 0.5 * std::atan(mu));
}

namespace
{

// A slip surface made concrete on a section: ground nodes up to the toe, the
// surface elevation z_s(x), and the kinks that must be slice edges.
struct surface_def
{
    bool ok = false;
    std::vector<double> nx, nz, nd;   // ground nodes (x, z, water-table depth)
    slip_surface_geometry arc;        // plane / circle / structural plane from A to the toe
    bool bilinear = false;
    double z0 = 0.0, zA = 0.0, x_toe = 0.0, z_toe = 0.0;
    double back_g = 0.0, basal_g = 0.0;
    bool scarp = false;
    double scarp_g = 0.0;
    int toe_index = -1;
    std::vector<double> brk;          // kinks inside (0, x_toe)

    double base(double x) const
    {
        if (bilinear) return std::max(zA - back_g * x, z_toe + basal_g * (x_toe - x));
        return arc.base(x);
    }
    double zs(double x) const
    {
        const double b = base(x);
        return scarp ? std::max(b, z0 - scarp_g * x) : b;
    }
};

surface_def build_def(const slip_section& sec, const slip_shape& sh)
{
    surface_def d;
    const int N = static_cast<int>(sec.x.size());
    if (N < 2 || !(sh.depth > 0.0)) return d;
    d.z0 = sec.z[0];
    d.zA = d.z0 - sh.depth;
    d.scarp = sh.rankine_scarp && sh.family != slip_surface_family::bilinear;
    d.scarp_g = sh.back_gradient;
    if (d.scarp && !(d.scarp_g > 0.0)) return d;

    if (sh.family == slip_surface_family::structural)
    {
        // Daylight of the plane from A at gradient g (with the scarp if any).
        const double g = sh.basal_gradient;
        auto zs_plain = [&](double x) {
            const double b = d.zA - g * x;
            return d.scarp ? std::max(b, d.z0 - d.scarp_g * x) : b;
        };
        double t_prev = sec.z[0] - zs_plain(0.0);
        int j_day = -1;
        for (int j = 1; j < N; ++j)
        {
            const double t = sec.z[static_cast<std::size_t>(j)] - zs_plain(sec.x[static_cast<std::size_t>(j)]);
            if (t <= 0.0) { j_day = j; break; }
            t_prev = t;
        }
        if (j_day < 0 || !(t_prev > 0.0)) return d;
        const std::size_t jp = static_cast<std::size_t>(j_day - 1), jj = static_cast<std::size_t>(j_day);
        const double tj = sec.z[jj] - zs_plain(sec.x[jj]);
        const double f = t_prev / (t_prev - tj);
        d.x_toe = sec.x[jp] + f * (sec.x[jj] - sec.x[jp]);
        d.z_toe = sec.z[jp] + f * (sec.z[jj] - sec.z[jp]);
        for (std::size_t q = 0; q <= jp; ++q)
        {
            d.nx.push_back(sec.x[q]); d.nz.push_back(sec.z[q]); d.nd.push_back(sec.water_table_depth[q]);
        }
        const double dp = sec.water_table_depth[jp], dj = sec.water_table_depth[jj];
        d.nx.push_back(d.x_toe); d.nz.push_back(d.z_toe); d.nd.push_back(dp + f * (dj - dp));
        // The plane passes through A and the daylight point.
        d.arc = make_slip_surface(d.x_toe, d.zA, d.zA - g * d.x_toe, 0.0);
        d.toe_index = j_day;
        if (!d.arc.admissible || !(d.z0 - d.z_toe > 0.0)) return d;
    }
    else
    {
        const int k = sh.toe_index;
        if (k <= 0 || k >= N) return d;
        const std::size_t kk = static_cast<std::size_t>(k);
        d.x_toe = sec.x[kk];
        d.z_toe = sec.z[kk];
        const double H = d.z0 - d.z_toe;
        if (!(d.x_toe > 0.0) || !(H > 0.0) || !(sh.depth < H)) return d;
        for (std::size_t q = 0; q <= kk; ++q)
        {
            d.nx.push_back(sec.x[q]); d.nz.push_back(sec.z[q]); d.nd.push_back(sec.water_table_depth[q]);
        }
        d.toe_index = k;
        if (sh.family == slip_surface_family::bilinear)
        {
            const double chord_g = (H - sh.depth) / d.x_toe;
            d.bilinear = true;
            d.back_g = sh.back_gradient;
            d.basal_g = sh.basal_gradient;
            if (!(d.basal_g >= 0.0) || d.basal_g > chord_g * (1.0 + 1e-12) || !(d.back_g > d.basal_g)) return d;
            const double xk = (d.zA - d.z_toe - d.basal_g * d.x_toe) / (d.back_g - d.basal_g);
            if (!(xk < d.x_toe)) return d;
            if (xk > 0.0) d.brk.push_back(xk);
            d.arc = make_slip_surface(d.x_toe, d.zA, d.z_toe, 0.0);   // unused, admissible
        }
        else
        {
            const double s = (sh.family == slip_surface_family::circular) ? sh.sagitta : 0.0;
            d.arc = make_slip_surface(d.x_toe, d.zA, d.z_toe, s);
            if (!d.arc.admissible) return d;
        }
    }
    if (d.scarp)
    {
        // First meeting of the scarp and the base: a slice edge.
        auto gap = [&](double x) { return d.base(x) - (d.z0 - d.scarp_g * x); };
        const int n_scan = 64;
        double xa = 0.0, ga = gap(0.0);
        for (int q = 1; q <= n_scan; ++q)
        {
            const double xb = d.x_toe * q / n_scan, gb = gap(xb);
            if (ga < 0.0 && gb >= 0.0)
            {
                const brent_root r = find_root_brent(gap, xa, xb, 1e-12 * d.x_toe, 200);
                if (r.converged && r.x > 0.0 && r.x < d.x_toe) d.brk.push_back(r.x);
                break;
            }
            xa = xb; ga = gb;
        }
    }
    std::sort(d.brk.begin(), d.brk.end());
    d.ok = true;
    return d;
}

struct slice
{
    double w, W, sa, ca, u, xm, zm;
};

double ordinary_F(const std::vector<slice>& sl, double c, double mu, double drive)
{
    double num = 0.0;
    for (const slice& s : sl) num += c * s.w / s.ca + mu * (s.W * s.ca - s.u * s.w / s.ca);
    return num / drive;
}

double spencer_F(const std::vector<slice>& sl, double c, double mu, double f_ord, double min_ma)
{
    const double F_lo = (std::isfinite(f_ord) && f_ord > 0.0) ? 0.01 * f_ord : 1e-3;
    const double F_hi = (std::isfinite(f_ord) && f_ord > 0.0) ? 100.0 * f_ord : 1e3;
    // Sum of Q (which = 0) or of its moment (which = 1) at (F, theta); NaN if invalid.
    auto sums = [&](double F, double th, int which) {
        double s_f = 0.0, s_m = 0.0;
        const double st = std::sin(th), ct = std::cos(th);
        for (const slice& s : sl)
        {
            const double a = std::atan2(s.sa, s.ca);
            const double den = std::cos(a + th) + std::sin(a + th) * mu / F;
            if (den < min_ma) return std::numeric_limits<double>::quiet_NaN();
            const double l = s.w / s.ca;
            const double Q = (c * l / F + mu * (s.W * s.ca - s.u * l) / F - s.W * s.sa) / den;
            s_f += Q;
            s_m += Q * (s.xm * st - s.zm * ct);
        }
        return which == 0 ? s_f : s_m;
    };
    auto solve_F = [&](double th, int which) {
        const int n = 20;
        double Fa = F_lo, fa = sums(Fa, th, which);
        for (int q = 1; q <= n; ++q)
        {
            const double Fb = F_lo * std::pow(F_hi / F_lo, static_cast<double>(q) / n);
            const double fb = sums(Fb, th, which);
            if (std::isfinite(fa) && std::isfinite(fb) && ((fa > 0.0) != (fb > 0.0)))
            {
                const brent_root r = find_root_brent([&](double F) { return sums(F, th, which); }, Fa, Fb,
                                                     1e-12 * Fb, 200);
                if (r.converged) return r.x;
            }
            Fa = Fb; fa = fb;
        }
        return std::numeric_limits<double>::quiet_NaN();
    };
    // A planar surface satisfies moment equilibrium identically when the
    // interslice forces are parallel to the base (theta = -alpha): Spencer's
    // solution is then the force solution, the ordinary F.
    double a_min = k_inf, a_max = -k_inf;
    for (const slice& s : sl) { const double a = std::atan2(s.sa, s.ca); a_min = std::min(a_min, a); a_max = std::max(a_max, a); }
    if (a_max - a_min < 1e-9) return f_ord;
    auto gap = [&](double th) { return solve_F(th, 0) - solve_F(th, 1); };
    // Scan theta coarsely over [-1.3, 1.3] and finely around -mean(alpha), where
    // nearly planar surfaces have their root; keep the crossing nearest 0.
    double a_mean = 0.0, w_sum = 0.0;
    for (const slice& s : sl) { a_mean += s.w * std::atan2(s.sa, s.ca); w_sum += s.w; }
    a_mean /= w_sum;
    std::vector<double> ths;
    for (int q = 0; q <= 26; ++q) ths.push_back(-1.3 + 2.6 * q / 26);
    for (int q = 0; q <= 24; ++q) ths.push_back(-a_mean - 0.3 + 0.6 * q / 24);
    std::sort(ths.begin(), ths.end());
    double best_abs = k_inf, best_a = 0.0, best_b = 0.0;
    double ta = ths[0], ga = gap(ta);
    for (std::size_t q = 1; q < ths.size(); ++q)
    {
        const double tb = ths[q], gb = gap(tb);
        if (tb > ta && std::isfinite(ga) && std::isfinite(gb) && ((ga > 0.0) != (gb > 0.0)))
        {
            const double mid = 0.5 * (ta + tb);
            if (std::fabs(mid) < best_abs) { best_abs = std::fabs(mid); best_a = ta; best_b = tb; }
        }
        ta = tb; ga = gb;
    }
    if (!std::isfinite(best_abs)) return k_inf;
    const brent_root r = find_root_brent(gap, best_a, best_b, 1e-10, 200);
    if (!r.converged) return k_inf;
    const double F = solve_F(r.x, 0);
    return std::isfinite(F) ? F : k_inf;
}

double evaluate_def(const surface_def& d, const rupture_parameters& p, const slip_evaluation_options& opt)
{
    if (!d.ok) return k_inf;
    const double rho_g = p.density_kg_m3 * p.gravity_m_s2;
    const double rw_g = p.water_density_kg_m3 * p.gravity_m_s2;
    const double c = p.cohesion_pa, mu = p.friction_coefficient;
    const double H = d.z0 - d.z_toe;
    const double tol = 1e-9 * H;
    const int m = std::max(1, opt.slices_per_cell);
    const std::size_t N = d.nx.size();

    thread_local std::vector<slice> sl;
    thread_local std::vector<double> edges;
    sl.clear();
    double drive = 0.0, drive_tan = 0.0;
    const double z_ref = d.zs(0.0);
    for (std::size_t j = 0; j + 1 < N; ++j)
    {
        const double x0 = d.nx[j], x1 = d.nx[j + 1], z0 = d.nz[j], z1 = d.nz[j + 1];
        const double seg = x1 - x0;
        if (!(seg > 0.0)) continue;
        const bool wet = std::isfinite(d.nd[j]) && std::isfinite(d.nd[j + 1]);
        const double zw0 = z0 - d.nd[j], zw1 = z1 - d.nd[j + 1];
        edges.clear();
        for (int q = 0; q <= m; ++q) edges.push_back((q == m) ? x1 : x0 + seg * q / m);
        for (double b : d.brk)
            if (b > x0 && b < x1) edges.push_back(b);
        std::sort(edges.begin(), edges.end());
        for (std::size_t e = 0; e + 1 < edges.size(); ++e)
        {
            const double xa = edges[e], xb = edges[e + 1], w = xb - xa;
            if (!(w > 1e-12 * seg)) continue;
            const double fa = (xa - x0) / seg, fb = (xb - x0) / seg;
            const double ga = z0 + (z1 - z0) * fa, gb = z0 + (z1 - z0) * fb;
            const double zba = d.zs(xa);
            const double zbb = (j + 2 == N && e + 2 == edges.size()) ? d.z_toe : d.zs(xb);
            const double ta = ga - zba, tb = gb - zbb;
            if (ta < -tol || tb < -tol) return k_inf;
            const double W = rho_g * w * 0.5 * (std::max(ta, 0.0) + std::max(tb, 0.0));
            const double alpha = std::atan((zba - zbb) / w);
            const double sa = std::sin(alpha), ca = std::cos(alpha);
            double u = 0.0;
            if (wet && rw_g > 0.0)
            {
                const double fm = 0.5 * (fa + fb);
                const double zw_mid = zw0 + (zw1 - zw0) * fm;
                u = rw_g * std::max(0.0, zw_mid - 0.5 * (zba + zbb));
            }
            drive += W * sa;
            drive_tan += W * sa / ca;
            sl.push_back({w, W, sa, ca, u, 0.5 * (xa + xb), 0.5 * (zba + zbb) - z_ref});
        }
    }
    if (!(drive > 0.0)) return k_inf;
    const double f_ord = ordinary_F(sl, c, mu, drive);
    switch (opt.method)
    {
    case limit_equilibrium_method::ordinary:
        return f_ord;
    case limit_equilibrium_method::spencer:
        return spencer_F(sl, c, mu, f_ord, opt.min_m_alpha);
    case limit_equilibrium_method::bishop:
    case limit_equilibrium_method::janbu:
        break;
    }
    const bool janbu = opt.method == limit_equilibrium_method::janbu;
    if (janbu && !(drive_tan > 0.0)) return k_inf;
    double F = (std::isfinite(f_ord) && f_ord > 0.0) ? f_ord : 1.0;
    for (int it = 0; it < opt.bishop_max_iterations; ++it)
    {
        double num = 0.0;
        for (const slice& s : sl)
        {
            const double ma = s.ca + mu * s.sa / F;
            if (ma < opt.min_m_alpha) return k_inf;
            const double a = c * s.w + mu * (s.W - s.u * s.w);
            num += janbu ? a / (s.ca * ma) : a / ma;
        }
        const double F_new = num / (janbu ? drive_tan : drive);
        if (!(F_new > 0.0)) return F_new;   // no strength left: unstable
        const bool done = std::fabs(F_new - F) <= opt.bishop_tolerance * std::max(1.0, std::fabs(F));
        F = F_new;
        if (done) break;
    }
    if (janbu && opt.janbu_correction)
    {
        // d / L of the surface below the line from its upper end to the toe.
        const double xs = 0.0, zs0 = d.zs(0.0), xe = d.nx.back(), ze = d.z_toe;
        const double L = std::hypot(xe - xs, ze - zs0);
        double dmax = 0.0;
        const int n_probe = 200;
        for (int q = 1; q < n_probe; ++q)
        {
            const double x = xe * q / n_probe;
            const double z_line = zs0 + (ze - zs0) * x / xe;
            dmax = std::max(dmax, (z_line - d.zs(x)) * (xe / L));   // perpendicular distance
        }
        const double b1 = (c > 0.0 && mu > 0.0) ? 0.50 : (mu > 0.0 ? 0.31 : 0.69);
        const double r = dmax / L;
        F *= 1.0 + b1 * (r - 1.4 * r * r);
    }
    return F;
}

// Grid then Brent on one variable; returns (x, f).
template <class Fn>
std::pair<double, double> grid_brent(Fn&& f, double lo, double hi, int n)
{
    double best_x = 0.5 * (lo + hi), best_f = k_inf;
    int best_i = -1;
    for (int i = 0; i < n; ++i)
    {
        const double x = lo + (hi - lo) * (i + 0.5) / n;
        const double v = f(x);
        if (v < best_f) { best_f = v; best_x = x; best_i = i; }
    }
    if (best_i < 0) return {best_x, k_inf};
    const double step = (hi - lo) / n;
    const double a = std::max(lo, best_x - step), b = std::min(hi, best_x + step);
    const brent_minimum r = minimise_brent_bounded(f, a, b, 1e-6 * (hi - lo), 100);
    if (r.f < best_f) return {r.x, r.f};
    return {best_x, best_f};
}

} // namespace

double evaluate_slip_shape(const slip_section& sec, const slip_shape& sh, const rupture_parameters& p,
                           const slip_evaluation_options& opt)
{
    return evaluate_def(build_def(sec, sh), p, opt);
}

namespace
{
// Body thickness (ground minus surface, >= 0) at x on the def's ground nodes.
double def_thickness(const surface_def& d, double x)
{
    std::size_t j = 0;
    while (j + 2 < d.nx.size() && x > d.nx[j + 1]) ++j;
    const double f = (x - d.nx[j]) / (d.nx[j + 1] - d.nx[j]);
    const double g = d.nz[j] + (d.nz[j + 1] - d.nz[j]) * f;
    return std::max(0.0, g - d.zs(x));
}
}

std::vector<double> slip_shape_thickness(const slip_section& sec, const slip_shape& sh)
{
    const surface_def d = build_def(sec, sh);
    std::vector<double> t;
    if (!d.ok) return t;
    // Each section point before the toe gets the MEAN thickness over its span
    // along the section (half-way to its neighbours, clipped to the body), by a
    // 32-point midpoint rule: a thin wedge between cell centres still counts.
    const std::size_t N = d.nx.size();
    for (std::size_t j = 0; j + 1 < N; ++j)
    {
        const double lo = (j == 0) ? 0.0 : 0.5 * (d.nx[j - 1] + d.nx[j]);
        const double hi = std::min(d.x_toe, 0.5 * (d.nx[j] + d.nx[j + 1]));
        if (!(hi > lo)) { t.push_back(0.0); continue; }
        const int q = 32;
        double sum = 0.0;
        for (int k = 0; k < q; ++k) sum += def_thickness(d, lo + (hi - lo) * (k + 0.5) / q);
        t.push_back(sum / q);
    }
    return t;
}

double slip_shape_max_thickness(const slip_section& sec, const slip_shape& sh)
{
    const surface_def d = build_def(sec, sh);
    if (!d.ok) return 0.0;
    double m = 0.0;
    const int q = 400;
    for (int k = 0; k <= q; ++k) m = std::max(m, def_thickness(d, d.x_toe * k / q));
    for (double b : d.brk) m = std::max(m, def_thickness(d, b));
    return m;
}

double evaluate_slip_surface(const slip_section& sec, int k, double h, double sagitta,
                             const rupture_parameters& p, const slip_evaluation_options& opt)
{
    slip_shape sh;
    sh.family = (sagitta > 0.0) ? slip_surface_family::circular : slip_surface_family::plane_integrated;
    sh.toe_index = k;
    sh.depth = h;
    sh.sagitta = sagitta;
    return evaluate_slip_shape(sec, sh, p, opt);
}

slip_search_result search_slip_shape(const slip_section& sec, int k, slip_surface_family family, slip_back_scarp scarp,
                                     const rupture_parameters& p, const slip_evaluation_options& opt, int search_points,
                                     double apparent_dip_gradient)
{
    slip_search_result out;
    const int N = static_cast<int>(sec.x.size());
    if (N < 2) return out;
    const int n = std::max(2, search_points);
    slip_shape base;
    base.family = family;
    base.back_gradient = rankine_back_gradient(p.friction_coefficient);
    base.rankine_scarp = scarp == slip_back_scarp::rankine && family != slip_surface_family::bilinear;
    auto F_of = [&](const slip_shape& s) { return evaluate_slip_shape(sec, s, p, opt); };

    slip_shape best = base;
    double best_f = k_inf;

    if (family == slip_surface_family::structural)
    {
        base.basal_gradient = apparent_dip_gradient;
        double zmin = sec.z[0];
        for (double z : sec.z) zmin = std::min(zmin, z);
        const double Hmax = sec.z[0] - zmin;
        if (!(Hmax > 0.0)) return out;
        auto f = [&](double h) { slip_shape s = base; s.depth = h; return F_of(s); };
        const auto r = grid_brent(f, 1e-9 * Hmax, Hmax * (1.0 - 1e-9), n);
        best = base;
        best.depth = r.first;
        best_f = r.second;
        if (!(best_f < k_inf)) return out;
        best.toe_index = build_def(sec, best).toe_index;
    }
    else
    {
        if (k <= 0 || k >= N) return out;
        const std::size_t kk = static_cast<std::size_t>(k);
        const double X = sec.x[kk], H = sec.z[0] - sec.z[kk];
        if (!(X > 0.0) || !(H > 0.0)) return out;
        base.toe_index = k;
        const double h_lo = 1e-9 * H, h_hi = H * (1.0 - 1e-9);

        // The plane first (s = 0, tan(beta) = chord), for every family.
        slip_shape plane = base;
        plane.family = slip_surface_family::plane_integrated;
        auto f_plane = [&](double h) { slip_shape s = plane; s.depth = h; return F_of(s); };
        const auto pl = grid_brent(f_plane, h_lo, h_hi, n);
        best = base;
        best.depth = pl.first;
        best_f = pl.second;
        // Shape parameter sg in [0, sg_hi]: circular s / L; bilinear tan(beta) / chord gradient.
        const bool two_d = family == slip_surface_family::circular || family == slip_surface_family::bilinear;
        if (two_d)
        {
            const bool circ = family == slip_surface_family::circular;
            auto chord = [&](double h) { return std::hypot(X, H - h); };
            auto make = [&](double h, double sg) {
                slip_shape s = base;
                s.depth = h;
                if (circ) s.sagitta = sg * chord(h);
                else s.basal_gradient = (1.0 - sg) * (H - h) / X;   // sg = 0 is the plane
                return s;
            };
            if (!circ) best = make(best.depth, 0.0);
            double best_sg = 0.0;
            const double sg_hi = circ ? 0.5 : 1.0;
            for (int i = 0; i < n; ++i)
            {
                const double h = h_lo + (h_hi - h_lo) * (i + 0.5) / n;
                for (int j = 1; j < n; ++j)
                {
                    const double sg = sg_hi * j / (n - 1);
                    const double v = F_of(make(h, sg));
                    if (v < best_f) { best_f = v; best = make(h, sg); best_sg = sg; }
                }
            }
            if (std::isfinite(best_f))
            {
                double dsg = sg_hi / (n - 1), dh = (h_hi - h_lo) / n;
                for (int round = 0; round < 2; ++round)
                {
                    const double hb = best.depth;
                    const brent_minimum rs = minimise_brent_bounded(
                        [&](double sg) { return F_of(make(hb, sg)); },
                        std::max(0.0, best_sg - dsg), std::min(sg_hi, best_sg + dsg), 1e-6 * sg_hi, 100);
                    if (rs.f < best_f) { best_f = rs.f; best_sg = rs.x; best = make(hb, rs.x); }
                    const double sb = best_sg;
                    const brent_minimum rh = minimise_brent_bounded(
                        [&](double h) { return F_of(make(h, sb)); },
                        std::max(h_lo, best.depth - dh), std::min(h_hi, best.depth + dh), 1e-6 * H, 100);
                    if (rh.f < best_f) { best_f = rh.f; best = make(rh.x, sb); }
                    dsg *= 0.5;
                    dh *= 0.5;
                }
            }
        }
    }
    if (!(best_f < k_inf)) return out;
    out.admissible = true;
    out.safety_factor = best_f;
    out.shape = best;
    out.depth = best.depth;
    out.sagitta = best.sagitta;
    out.toe_index = best.toe_index;
    return out;
}

slip_search_result search_slip_surface(const slip_section& sec, int k, slip_surface_family family,
                                       const rupture_parameters& p, const slip_evaluation_options& opt,
                                       int search_points)
{
    return search_slip_shape(sec, k, family, slip_back_scarp::vertical, p, opt, search_points);
}

void validate_slip_surface_options(const slip_surface_options& o)
{
    if (o.family == slip_surface_family::ellipsoid)
        throw std::invalid_argument("slip surfaces: the ellipsoid family is computed by compute_ellipsoid_surfaces");
    if (o.evaluation.slices_per_cell < 1)
        throw std::invalid_argument("slip surfaces: slices per cell must be at least 1");
    if (o.search_points < 2)
        throw std::invalid_argument("slip surfaces: search points must be at least 2");
    if (o.toe_window < 0)
        throw std::invalid_argument("slip surfaces: the toe window must be >= 0");
    if (!(o.max_length_metres >= 0.0))
        throw std::invalid_argument("slip surfaces: the maximum length must be >= 0 (0 = no cap)");
    if (std::isnan(o.screen_safety_factor))
        throw std::invalid_argument("slip surfaces: the screen safety factor must not be NaN");
    if (o.family == slip_surface_family::bilinear && o.back_scarp == slip_back_scarp::rankine)
        throw std::invalid_argument("slip surfaces: bilinear already has a Rankine back segment; use the vertical back scarp");
    if (o.family == slip_surface_family::structural)
    {
        if (!std::isfinite(o.structural_dip_gradient) || !std::isfinite(o.structural_dip_direction_degrees))
            throw std::invalid_argument("slip surfaces: the structural dip gradient and direction must be finite");
    }
}

slip_surface_rasters compute_slip_surfaces(
    const lsdtt_flow_base& flow,
    const lsdtt_raster_cpu& geometry,
    const std::vector<std::uint8_t>& channel,
    const std::vector<std::int32_t>& region,
    const std::map<std::int32_t, rupture_parameters>& by_region,
    const rupture_rasters& rupture,
    const slip_surface_options& options,
    int n_threads,
    const std::vector<double>& water_table_depth)
{
    validate_slip_surface_options(options);
    const lsdtt_raster_meta& meta = geometry.get_meta();
    const int nc = meta.n_cols;
    const std::size_t n = static_cast<std::size_t>(meta.n_rows) * static_cast<std::size_t>(nc);
    if (rupture.daylight.size() != n || channel.size() != n)
        throw std::invalid_argument("compute_slip_surfaces: the rupture result or channel mask differs from the grid");
    if (!region.empty() && region.size() != n)
        throw std::invalid_argument("compute_slip_surfaces: the region grid differs from the grid");
    if (!water_table_depth.empty() && water_table_depth.size() != n)
        throw std::invalid_argument("compute_slip_surfaces: the water-table grid differs from the grid");
    if (by_region.empty() || (region.empty() && by_region.size() != 1))
        throw std::invalid_argument("compute_slip_surfaces: one parameter set without regions, or one per region");
    const bool structural = options.family == slip_surface_family::structural;
    if (structural && ((!options.structural_dip_gradient_raster.empty() && options.structural_dip_gradient_raster.size() != n)
                       || (!options.structural_dip_direction_raster.empty() && options.structural_dip_direction_raster.size() != n)))
        throw std::invalid_argument("compute_slip_surfaces: a structural raster differs from the grid");
    const rupture_parameters* single = region.empty() ? &by_region.begin()->second : nullptr;

    const float nd = static_cast<float>(meta.no_data_value);
    std::vector<float> fs(n, nd), hc(n, nd), hmax(n, nd), sg(n, nd), bg(n, nd), len(n, nd), sin_(n, nd), uns(n, nd), ev(n, nd);
    std::vector<std::uint8_t> m_eval(n, 0), m_flag(n, 0);
    slip_surface_rasters out;
    out.toe.assign(n, -1);
    out.shape.assign(n, slip_shape{});
    const std::vector<float>& rf = rupture.safety_factor.get_data();
    const std::vector<std::uint8_t>& rm = rupture.safety_factor.get_mask();
    std::int64_t n_eval = 0, n_screen = 0, n_unst = 0, n_inadm = 0;

#ifdef _OPENMP
    const int threads = (n_threads > 0) ? n_threads : omp_get_max_threads();
#pragma omp parallel for schedule(dynamic, 16) num_threads(threads) reduction(+ : n_eval, n_screen, n_unst, n_inadm)
#else
    (void)n_threads;
#endif
    for (std::int64_t ii = 0; ii < static_cast<std::int64_t>(n); ++ii)
    {
        const std::size_t i = static_cast<std::size_t>(ii);
        if (!rm[i] || rupture.daylight[i] < 0) continue;
        if (!single && region[i] == k_rupture_region_nodata) continue;
        m_flag[i] = 1;
        if (!(static_cast<double>(rf[i]) <= options.screen_safety_factor))
        {
            ev[i] = 0.0f;
            ++n_screen;
            continue;
        }
        ev[i] = 1.0f;
        ++n_eval;
        const rupture_parameters& p = single ? *single : by_region.at(region[i]);
        const slip_section sec = build_slip_section(flow, geometry, channel, ii, water_table_depth);
        const int K = static_cast<int>(sec.x.size()) - 1;
        slip_search_result best;
        if (structural)
        {
            const double dip_g = options.structural_dip_gradient_raster.empty()
                                     ? options.structural_dip_gradient : options.structural_dip_gradient_raster[i];
            const double dip_dir = options.structural_dip_direction_raster.empty()
                                       ? options.structural_dip_direction_degrees : options.structural_dip_direction_raster[i];
            const double bearing = grid_bearing_degrees(ii, rupture.daylight[i], nc);
            if (std::isfinite(dip_g) && std::isfinite(dip_dir) && std::isfinite(bearing))
            {
                const double g = dip_g * std::cos((bearing - dip_dir) * M_PI / 180.0);
                best = search_slip_shape(sec, -1, options.family, options.back_scarp, p, options.evaluation,
                                         options.search_points, g);
            }
        }
        else
        {
            int k_lo = 1, k_hi = K;
            if (!options.full_search)
            {
                int ks = -1;
                for (int k = 1; k <= K; ++k)
                    if (sec.cell[static_cast<std::size_t>(k)] == rupture.daylight[i]) { ks = k; break; }
                if (ks < 0) { ++n_inadm; continue; }
                k_lo = std::max(1, ks - options.toe_window);
                k_hi = std::min(K, ks + options.toe_window);
            }
            for (int k = k_lo; k <= k_hi; ++k)
            {
                if (options.full_search && options.max_length_metres > 0.0
                    && sec.x[static_cast<std::size_t>(k)] > options.max_length_metres) break;
                const slip_search_result r = search_slip_shape(sec, k, options.family, options.back_scarp, p,
                                                               options.evaluation, options.search_points);
                if (r.admissible && r.safety_factor < best.safety_factor) best = r;
            }
        }
        if (!best.admissible) { ++n_inadm; continue; }
        const double tmax = slip_shape_max_thickness(sec, best.shape);
        const int kt = std::min(best.shape.toe_index, K);
        const surface_def d = build_def(sec, best.shape);
        m_eval[i] = 1;
        fs[i] = static_cast<float>(best.safety_factor);
        hc[i] = static_cast<float>(best.shape.depth);
        hmax[i] = static_cast<float>(tmax);
        sg[i] = static_cast<float>(best.shape.sagitta);
        bg[i] = static_cast<float>(best.shape.basal_gradient);
        len[i] = static_cast<float>(d.x_toe);
        sin_[i] = static_cast<float>(slip_section_sinuosity(sec, kt, nc));
        uns[i] = (best.safety_factor <= 1.0) ? 1.0f : 0.0f;
        if (best.safety_factor <= 1.0) ++n_unst;
        out.toe[i] = sec.cell[static_cast<std::size_t>(kt)];
        out.shape[i] = best.shape;
    }

    // Bodies of unstable crowns, in crown order (deterministic).
    std::vector<double> env(n, 0.0);
    out.body_offset.push_back(0);
    for (std::size_t i = 0; i < n; ++i)
    {
        if (!m_eval[i] || !(fs[i] <= 1.0f)) continue;
        const slip_section sec = build_slip_section(flow, geometry, channel, static_cast<std::int64_t>(i), {});
        const std::vector<double> t = slip_shape_thickness(sec, out.shape[i]);
        out.body_crown.push_back(static_cast<std::int64_t>(i));
        for (std::size_t j = 0; j < t.size(); ++j)
        {
            if (!(t[j] > 0.0)) continue;
            const std::size_t c = static_cast<std::size_t>(sec.cell[j]);
            out.body_cell.push_back(sec.cell[j]);
            out.body_cell_depth.push_back(t[j]);
            env[c] = std::max(env[c], t[j]);
        }
        out.body_offset.push_back(static_cast<std::int64_t>(out.body_cell.size()));
    }
    std::vector<float> bd(n, nd);
    std::vector<std::uint8_t> bm(n, 0);
    for (std::size_t i = 0; i < n; ++i)
        if (env[i] > 0.0) { bd[i] = static_cast<float>(env[i]); bm[i] = 1; }

    auto make = [&](std::vector<float>&& v, const std::vector<std::uint8_t>& m) {
        return lsdtt_raster_cpu(std::move(v), std::vector<std::uint8_t>(m), meta);
    };
    out.safety_factor = make(std::move(fs), m_eval);
    out.crown_depth = make(std::move(hc), m_eval);
    out.max_depth = make(std::move(hmax), m_eval);
    out.sagitta = make(std::move(sg), m_eval);
    out.basal_gradient = make(std::move(bg), m_eval);
    out.length = make(std::move(len), m_eval);
    out.sinuosity = make(std::move(sin_), m_eval);
    out.unstable = make(std::move(uns), m_eval);
    out.evaluated = make(std::move(ev), m_flag);
    out.body_depth = make(std::move(bd), bm);
    out.n_evaluated = n_eval;
    out.n_screened = n_screen;
    out.n_unstable = n_unst;
    out.n_inadmissible = n_inadm;
    return out;
}

rupture_probability_rasters compute_slip_failure_probability(
    const slip_surface_rasters& slip,
    const lsdtt_flow_base& flow,
    const lsdtt_raster_cpu& geometry,
    const std::vector<std::uint8_t>& channel,
    const std::vector<std::int32_t>& region,
    const std::map<std::int32_t, std::vector<rupture_parameter_sample>>& samples_by_region,
    const slip_evaluation_options& evaluation,
    int n_threads,
    const std::vector<double>& water_table_depth)
{
    const lsdtt_raster_meta& meta = geometry.get_meta();
    const std::size_t n = static_cast<std::size_t>(meta.n_rows) * static_cast<std::size_t>(meta.n_cols);
    if (slip.toe.size() != n || slip.shape.size() != n)
        throw std::invalid_argument("compute_slip_failure_probability: the slip result differs from the grid");
    if (samples_by_region.empty() || (region.empty() && samples_by_region.size() != 1))
        throw std::invalid_argument("compute_slip_failure_probability: one sample list without regions, or one per region");
    std::size_t max_samples = 0;
    for (const auto& kv : samples_by_region)
    {
        double w = 0.0;
        for (const auto& s : kv.second) w += s.weight;
        if (kv.second.empty() || std::fabs(w - 1.0) > 1e-9)
            throw std::invalid_argument("compute_slip_failure_probability: empty sample list or weights not summing to 1 (region "
                                        + std::to_string(kv.first) + ")");
        max_samples = std::max(max_samples, kv.second.size());
    }
    const std::vector<rupture_parameter_sample>* single = region.empty() ? &samples_by_region.begin()->second : nullptr;
    const std::vector<float>& ev = slip.evaluated.get_data();
    const std::vector<std::uint8_t>& em = slip.evaluated.get_mask();
    const std::vector<std::uint8_t>& fm = slip.safety_factor.get_mask();
    const float nd = static_cast<float>(meta.no_data_value);
    std::vector<float> pf(n, nd), mf(n, nd), sf(n, nd);
    std::vector<std::uint8_t> m(n, 0), m_stats(n, 0);
    std::int64_t n_src = 0, n_cert = 0, n_poss = 0;

#ifdef _OPENMP
    const int threads = (n_threads > 0) ? n_threads : omp_get_max_threads();
#pragma omp parallel for schedule(dynamic, 16) num_threads(threads) reduction(+ : n_src, n_cert, n_poss)
#else
    (void)n_threads;
#endif
    for (std::int64_t ii = 0; ii < static_cast<std::int64_t>(n); ++ii)
    {
        const std::size_t i = static_cast<std::size_t>(ii);
        if (!em[i]) continue;
        if (ev[i] == 0.0f) { pf[i] = 0.0f; m[i] = 1; continue; }   // screened out: Pf 0, flagged by `evaluated`
        if (!fm[i]) continue;                                      // no admissible surface
        const auto it = single ? samples_by_region.end() : samples_by_region.find(region[i]);
        if (!single && it == samples_by_region.end()) continue;
        const std::vector<rupture_parameter_sample>& smp = single ? *single : it->second;
        const slip_section sec = build_slip_section(flow, geometry, channel, ii, water_table_depth);
        std::vector<double> F(smp.size());
        double p_fail = 0.0, mean = 0.0;
        for (std::size_t s = 0; s < smp.size(); ++s)
        {
            F[s] = evaluate_slip_shape(sec, slip.shape[i], smp[s].parameters, evaluation);
            if (F[s] <= 1.0) p_fail += smp[s].weight;
            mean += smp[s].weight * F[s];
        }
        double var = 0.0;
        for (std::size_t s = 0; s < smp.size(); ++s) var += smp[s].weight * (F[s] - mean) * (F[s] - mean);
        ++n_src;
        if (p_fail >= 1.0 - 1e-12) ++n_cert;
        if (p_fail > 0.0) ++n_poss;
        m[i] = 1;
        m_stats[i] = 1;
        pf[i] = static_cast<float>(std::min(1.0, p_fail));
        mf[i] = static_cast<float>(mean);
        sf[i] = static_cast<float>(std::sqrt(var));
    }
    rupture_probability_rasters out;
    out.failure_probability = lsdtt_raster_cpu(std::move(pf), std::move(m), meta);
    out.mean_safety_factor = lsdtt_raster_cpu(std::move(mf), std::vector<std::uint8_t>(m_stats), meta);
    out.sd_safety_factor = lsdtt_raster_cpu(std::move(sf), std::move(m_stats), meta);
    out.n_sources = n_src;
    out.n_certain_failure = n_cert;
    out.n_possible_failure = n_poss;
    out.max_samples = max_samples;
    return out;
}

lsdtt_raster_cpu compute_slip_spencer_check(
    const slip_surface_rasters& slip,
    const lsdtt_flow_base& flow,
    const lsdtt_raster_cpu& geometry,
    const std::vector<std::uint8_t>& channel,
    const std::vector<std::int32_t>& region,
    const std::map<std::int32_t, rupture_parameters>& by_region,
    const slip_evaluation_options& evaluation,
    int n_threads,
    const std::vector<double>& water_table_depth,
    std::int64_t* n_no_solution)
{
    const lsdtt_raster_meta& meta = geometry.get_meta();
    const std::size_t n = static_cast<std::size_t>(meta.n_rows) * static_cast<std::size_t>(meta.n_cols);
    if (slip.shape.size() != n) throw std::invalid_argument("compute_slip_spencer_check: the slip result differs from the grid");
    const rupture_parameters* single = region.empty() ? &by_region.begin()->second : nullptr;
    slip_evaluation_options opt = evaluation;
    opt.method = limit_equilibrium_method::spencer;
    const std::vector<std::uint8_t>& fm = slip.safety_factor.get_mask();
    const float nd = static_cast<float>(meta.no_data_value);
    std::vector<float> v(n, nd);
    std::vector<std::uint8_t> m(n, 0);
    std::int64_t none = 0;
#ifdef _OPENMP
    const int threads = (n_threads > 0) ? n_threads : omp_get_max_threads();
#pragma omp parallel for schedule(dynamic, 16) num_threads(threads) reduction(+ : none)
#else
    (void)n_threads;
#endif
    for (std::int64_t ii = 0; ii < static_cast<std::int64_t>(n); ++ii)
    {
        const std::size_t i = static_cast<std::size_t>(ii);
        if (!fm[i]) continue;
        const rupture_parameters& p = single ? *single : by_region.at(region[i]);
        const slip_section sec = build_slip_section(flow, geometry, channel, ii, water_table_depth);
        const double F = evaluate_slip_shape(sec, slip.shape[i], p, opt);
        if (!std::isfinite(F)) { ++none; continue; }
        v[i] = static_cast<float>(F);
        m[i] = 1;
    }
    if (n_no_solution) *n_no_solution = none;
    return lsdtt_raster_cpu(std::move(v), std::move(m), meta);
}

slip_body_catalogues build_slip_body_catalogues(const slip_surface_rasters& slip, const lsdtt_raster_meta& meta)
{
    const int nr = meta.n_rows, nc = meta.n_cols;
    const std::size_t n = static_cast<std::size_t>(nr) * static_cast<std::size_t>(nc);
    const double b = meta.cell_size, area = b * b;
    slip_body_catalogues out;
    const std::vector<float>& fs = slip.safety_factor.get_data();
    const std::size_t n_bodies = slip.body_crown.size();

    // By toe: groups of unstable crowns with a non-empty body.
    std::map<std::int64_t, std::vector<std::size_t>> groups;   // toe -> body indices
    for (std::size_t k = 0; k < n_bodies; ++k)
        if (slip.body_offset[k + 1] > slip.body_offset[k])
            groups[slip.toe[static_cast<std::size_t>(slip.body_crown[k])]].push_back(k);
    std::int32_t id = 0;
    std::map<std::int64_t, double> depth;   // cell -> max depth within the group
    for (const auto& kv : groups)
    {
        depth.clear();
        slip_body_record rec;
        rec.id = ++id;
        rec.toe_x = meta.x_min + (static_cast<double>(kv.first % nc) + 0.5) * b;
        rec.toe_y = meta.y_min + (static_cast<double>(nr - kv.first / nc) - 0.5) * b;
        rec.min_safety_factor = k_inf;
        double sum_f = 0.0;
        for (std::size_t k : kv.second)
        {
            const double f = fs[static_cast<std::size_t>(slip.body_crown[k])];
            rec.min_safety_factor = std::min(rec.min_safety_factor, f);
            sum_f += f;
            for (std::int64_t q = slip.body_offset[k]; q < slip.body_offset[k + 1]; ++q)
            {
                double& d = depth[slip.body_cell[static_cast<std::size_t>(q)]];
                d = std::max(d, slip.body_cell_depth[static_cast<std::size_t>(q)]);
            }
        }
        rec.n_crowns = static_cast<std::int64_t>(kv.second.size());
        rec.mean_safety_factor = sum_f / static_cast<double>(rec.n_crowns);
        cell_shape_moments mom;
        for (const auto& cd : depth)
        {
            mom.add(static_cast<int>(cd.first / nc), static_cast<int>(cd.first % nc));
            rec.volume_m3 += cd.second * area;
            rec.max_depth_m = std::max(rec.max_depth_m, cd.second);
        }
        rec.n_cells = static_cast<std::int64_t>(depth.size());
        rec.area_m2 = static_cast<double>(rec.n_cells) * area;
        rec.mean_depth_m = rec.volume_m3 / rec.area_m2;
        const cell_shape sh = shape_from_moments(mom, meta);
        rec.centroid_x = sh.centroid_x; rec.centroid_y = sh.centroid_y;
        rec.length_m = sh.length_m; rec.width_m = sh.width_m; rec.length_to_width = sh.length_to_width;
        out.by_toe.push_back(rec);
    }

    // Amalgamated: 8-connected components of the envelope; each body joins the
    // component of its first cell.
    const std::vector<float>& env = slip.body_depth.get_data();
    const std::vector<std::uint8_t>& em = slip.body_depth.get_mask();
    out.label_amalgamated.assign(n, 0);
    std::vector<std::size_t> stack;
    std::int32_t cid = 0;
    for (std::size_t s0 = 0; s0 < n; ++s0)
    {
        if (!em[s0] || out.label_amalgamated[s0]) continue;
        slip_body_record rec;
        rec.id = ++cid;
        cell_shape_moments mom;
        stack.assign(1, s0);
        out.label_amalgamated[s0] = cid;
        while (!stack.empty())
        {
            const std::size_t c = stack.back();
            stack.pop_back();
            const int r = static_cast<int>(c / static_cast<std::size_t>(nc)), col = static_cast<int>(c % static_cast<std::size_t>(nc));
            mom.add(r, col);
            ++rec.n_cells;
            rec.volume_m3 += static_cast<double>(env[c]) * area;
            rec.max_depth_m = std::max(rec.max_depth_m, static_cast<double>(env[c]));
            for (int dr = -1; dr <= 1; ++dr)
                for (int dc = -1; dc <= 1; ++dc)
                {
                    if (!dr && !dc) continue;
                    const int rr = r + dr, cc = col + dc;
                    if (rr < 0 || rr >= nr || cc < 0 || cc >= nc) continue;
                    const std::size_t q = static_cast<std::size_t>(rr) * nc + cc;
                    if (em[q] && !out.label_amalgamated[q]) { out.label_amalgamated[q] = cid; stack.push_back(q); }
                }
        }
        rec.area_m2 = static_cast<double>(rec.n_cells) * area;
        rec.mean_depth_m = rec.volume_m3 / rec.area_m2;
        rec.min_safety_factor = k_inf;
        const cell_shape sh = shape_from_moments(mom, meta);
        rec.centroid_x = sh.centroid_x; rec.centroid_y = sh.centroid_y;
        rec.length_m = sh.length_m; rec.width_m = sh.width_m; rec.length_to_width = sh.length_to_width;
        out.amalgamated.push_back(rec);
    }
    std::vector<std::vector<std::int64_t>> toes(out.amalgamated.size());
    std::vector<double> sum_f(out.amalgamated.size(), 0.0);
    for (std::size_t k = 0; k < n_bodies; ++k)
    {
        if (!(slip.body_offset[k + 1] > slip.body_offset[k])) continue;
        const std::int32_t lab = out.label_amalgamated[static_cast<std::size_t>(slip.body_cell[static_cast<std::size_t>(slip.body_offset[k])])];
        if (lab <= 0) continue;
        slip_body_record& rec = out.amalgamated[static_cast<std::size_t>(lab - 1)];
        const std::size_t crown = static_cast<std::size_t>(slip.body_crown[k]);
        ++rec.n_crowns;
        sum_f[static_cast<std::size_t>(lab - 1)] += fs[crown];
        rec.min_safety_factor = std::min(rec.min_safety_factor, static_cast<double>(fs[crown]));
        toes[static_cast<std::size_t>(lab - 1)].push_back(slip.toe[crown]);
    }
    for (std::size_t q = 0; q < out.amalgamated.size(); ++q)
    {
        slip_body_record& rec = out.amalgamated[q];
        std::sort(toes[q].begin(), toes[q].end());
        rec.n_toes = static_cast<std::int64_t>(std::unique(toes[q].begin(), toes[q].end()) - toes[q].begin());
        rec.mean_safety_factor = rec.n_crowns ? sum_f[q] / static_cast<double>(rec.n_crowns) : std::numeric_limits<double>::quiet_NaN();
        if (!rec.n_crowns) rec.min_safety_factor = std::numeric_limits<double>::quiet_NaN();
    }
    return out;
}

} // namespace lsdtt
