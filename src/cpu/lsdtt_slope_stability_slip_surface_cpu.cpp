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

double evaluate_slip_surface(const slip_section& sec, int k, double h, double sagitta,
                             const rupture_parameters& p, const slip_evaluation_options& opt)
{
    if (k <= 0 || k >= static_cast<int>(sec.x.size())) return k_inf;
    const std::size_t kk = static_cast<std::size_t>(k);
    const double X = sec.x[kk], H = sec.z[0] - sec.z[kk];
    if (!(X > 0.0) || !(H > 0.0) || !(h > 0.0) || !(h < H)) return k_inf;
    const slip_surface_geometry g = make_slip_surface(X, sec.z[0] - h, sec.z[kk], sagitta);
    if (!g.admissible) return k_inf;

    const double rho_g = p.density_kg_m3 * p.gravity_m_s2;
    const double rw_g = p.water_density_kg_m3 * p.gravity_m_s2;
    const double c = p.cohesion_pa, mu = p.friction_coefficient;
    const double tol = 1e-9 * H;
    const int m = std::max(1, opt.slices_per_cell);
    const bool bishop = opt.method == limit_equilibrium_method::bishop;

    // Slice terms. For Bishop: a_i = c w + mu (W - u w), with sin, cos kept.
    double drive = 0.0, ordinary_num = 0.0;
    thread_local std::vector<double> s_a, s_sin, s_cos;
    if (bishop) { s_a.clear(); s_sin.clear(); s_cos.clear(); }
    for (std::size_t j = 0; j < kk; ++j)
    {
        const double x0 = sec.x[j], x1 = sec.x[j + 1], z0 = sec.z[j], z1 = sec.z[j + 1];
        const double d0 = sec.water_table_depth[j], d1 = sec.water_table_depth[j + 1];
        const bool wet = std::isfinite(d0) && std::isfinite(d1);
        const double zw0 = z0 - d0, zw1 = z1 - d1;
        const double seg = x1 - x0;
        double xa = x0, zba = g.base(x0), ga = z0;
        for (int q = 1; q <= m; ++q)
        {
            const double f = static_cast<double>(q) / m;
            const double xb = (q == m) ? x1 : x0 + seg * f;
            const double gb = z0 + (z1 - z0) * f;
            const double zbb = (j + 1 == kk && q == m) ? sec.z[kk] : g.base(xb);
            const double ta = ga - zba, tb = gb - zbb;
            if (ta < -tol || tb < -tol) return k_inf;
            const double w = xb - xa;
            const double W = rho_g * w * 0.5 * (std::max(ta, 0.0) + std::max(tb, 0.0));
            const double alpha = std::atan((zba - zbb) / w);
            const double sa = std::sin(alpha), ca = std::cos(alpha);
            double u = 0.0;
            if (wet && rw_g > 0.0)
            {
                const double fm = (static_cast<double>(q) - 0.5) / m;
                const double zw_mid = zw0 + (zw1 - zw0) * fm;
                u = rw_g * std::max(0.0, zw_mid - 0.5 * (zba + zbb));
            }
            drive += W * sa;
            ordinary_num += c * w / ca + mu * (W * ca - u * w / ca);
            if (bishop)
            {
                s_a.push_back(c * w + mu * (W - u * w));
                s_sin.push_back(sa);
                s_cos.push_back(ca);
            }
            xa = xb; zba = zbb; ga = gb;
        }
    }
    if (!(drive > 0.0)) return k_inf;
    const double f_ord = ordinary_num / drive;
    if (!bishop) return f_ord;

    double F = (std::isfinite(f_ord) && f_ord > 0.0) ? f_ord : 1.0;
    for (int it = 0; it < opt.bishop_max_iterations; ++it)
    {
        double num = 0.0;
        for (std::size_t i = 0; i < s_a.size(); ++i)
        {
            const double ma = s_cos[i] + mu * s_sin[i] / F;
            if (ma < opt.min_m_alpha) return k_inf;
            num += s_a[i] / ma;
        }
        const double F_new = num / drive;
        if (!(F_new > 0.0)) return F_new;   // no strength left: unstable
        const bool done = std::fabs(F_new - F) <= opt.bishop_tolerance * std::max(1.0, std::fabs(F));
        F = F_new;
        if (done) break;
    }
    return F;
}

namespace
{

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

slip_search_result search_slip_surface(const slip_section& sec, int k, slip_surface_family family,
                                       const rupture_parameters& p, const slip_evaluation_options& opt,
                                       int search_points)
{
    slip_search_result out;
    if (k <= 0 || k >= static_cast<int>(sec.x.size())) return out;
    const std::size_t kk = static_cast<std::size_t>(k);
    const double X = sec.x[kk], H = sec.z[0] - sec.z[kk];
    if (!(X > 0.0) || !(H > 0.0)) return out;
    const int n = std::max(2, search_points);
    const double h_lo = 1e-9 * H, h_hi = H * (1.0 - 1e-9);

    // The plane (s = 0) first, for both families.
    auto F_plane = [&](double h) { return evaluate_slip_surface(sec, k, h, 0.0, p, opt); };
    const auto pl = grid_brent(F_plane, h_lo, h_hi, n);
    double best_f = pl.second, best_h = pl.first, best_sig = 0.0;

    if (family == slip_surface_family::circular)
    {
        auto chord = [&](double h) { return std::hypot(X, H - h); };
        auto F2 = [&](double h, double sig) { return evaluate_slip_surface(sec, k, h, sig * chord(h), p, opt); };
        const double sig_hi = 0.5;
        for (int i = 0; i < n; ++i)
        {
            const double h = h_lo + (h_hi - h_lo) * (i + 0.5) / n;
            for (int j = 1; j < n; ++j)
            {
                const double sig = sig_hi * j / (n - 1);
                const double v = F2(h, sig);
                if (v < best_f) { best_f = v; best_h = h; best_sig = sig; }
            }
        }
        if (std::isfinite(best_f))
        {
            double dsig = sig_hi / (n - 1), dh = (h_hi - h_lo) / n;
            for (int round = 0; round < 2; ++round)
            {
                const double s_a = std::max(0.0, best_sig - dsig), s_b = std::min(sig_hi, best_sig + dsig);
                const brent_minimum rs = minimise_brent_bounded([&](double sg) { return F2(best_h, sg); },
                                                                s_a, s_b, 1e-6 * sig_hi, 100);
                if (rs.f < best_f) { best_f = rs.f; best_sig = rs.x; }
                const double h_a = std::max(h_lo, best_h - dh), h_b = std::min(h_hi, best_h + dh);
                const brent_minimum rh = minimise_brent_bounded([&](double h) { return F2(h, best_sig); },
                                                                h_a, h_b, 1e-6 * H, 100);
                if (rh.f < best_f) { best_f = rh.f; best_h = rh.x; }
                dsig *= 0.5;
                dh *= 0.5;
            }
        }
        out.sagitta = best_sig * chord(best_h);
    }
    if (!(best_f < k_inf)) return out;
    out.admissible = true;
    out.safety_factor = best_f;
    out.depth = best_h;
    out.toe_index = k;
    return out;
}

void validate_slip_surface_options(const slip_surface_options& o)
{
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
    const rupture_parameters* single = region.empty() ? &by_region.begin()->second : nullptr;

    const float nd = static_cast<float>(meta.no_data_value);
    std::vector<float> fs(n, nd), hc(n, nd), hmax(n, nd), sg(n, nd), len(n, nd), sin_(n, nd), uns(n, nd), ev(n, nd);
    std::vector<std::uint8_t> m_eval(n, 0), m_flag(n, 0);
    slip_surface_rasters out;
    out.toe.assign(n, -1);
    out.toe_index.assign(n, -1);
    out.depth_h.assign(n, 0.0);
    out.sagitta_s.assign(n, 0.0);
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
        slip_search_result best;
        for (int k = k_lo; k <= k_hi; ++k)
        {
            if (options.full_search && options.max_length_metres > 0.0
                && sec.x[static_cast<std::size_t>(k)] > options.max_length_metres) break;
            const slip_search_result r = search_slip_surface(sec, k, options.family, p, options.evaluation,
                                                             options.search_points);
            if (r.admissible && r.safety_factor < best.safety_factor) best = r;
        }
        if (!best.admissible) { ++n_inadm; continue; }
        const std::size_t kk = static_cast<std::size_t>(best.toe_index);
        const slip_surface_geometry g = make_slip_surface(sec.x[kk], sec.z[0] - best.depth, sec.z[kk], best.sagitta);
        double tmax = 0.0;
        for (std::size_t j = 0; j < kk; ++j) tmax = std::max(tmax, sec.z[j] - g.base(sec.x[j]));
        m_eval[i] = 1;
        fs[i] = static_cast<float>(best.safety_factor);
        hc[i] = static_cast<float>(best.depth);
        hmax[i] = static_cast<float>(tmax);
        sg[i] = static_cast<float>(best.sagitta);
        len[i] = static_cast<float>(sec.x[kk]);
        sin_[i] = static_cast<float>(slip_section_sinuosity(sec, best.toe_index, nc));
        uns[i] = (best.safety_factor <= 1.0) ? 1.0f : 0.0f;
        if (best.safety_factor <= 1.0) ++n_unst;
        out.toe[i] = sec.cell[kk];
        out.toe_index[i] = best.toe_index;
        out.depth_h[i] = best.depth;
        out.sagitta_s[i] = best.sagitta;
    }

    // Bodies of unstable crowns, in crown order (deterministic).
    std::vector<double> env(n, 0.0);
    out.body_offset.push_back(0);
    for (std::size_t i = 0; i < n; ++i)
    {
        if (!m_eval[i] || !(fs[i] <= 1.0f)) continue;
        const slip_section sec = build_slip_section(flow, geometry, channel, static_cast<std::int64_t>(i), {});
        const std::size_t kk = static_cast<std::size_t>(out.toe_index[i]);
        const slip_surface_geometry g = make_slip_surface(sec.x[kk], sec.z[0] - out.depth_h[i], sec.z[kk], out.sagitta_s[i]);
        out.body_crown.push_back(static_cast<std::int64_t>(i));
        for (std::size_t j = 0; j < kk; ++j)
        {
            const double t = sec.z[j] - g.base(sec.x[j]);
            if (!(t > 0.0)) continue;
            const std::size_t c = static_cast<std::size_t>(sec.cell[j]);
            out.body_cell.push_back(sec.cell[j]);
            out.body_cell_depth.push_back(t);
            env[c] = std::max(env[c], t);
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
    if (slip.toe.size() != n) throw std::invalid_argument("compute_slip_failure_probability: the slip result differs from the grid");
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
    const float nd = static_cast<float>(meta.no_data_value);
    std::vector<float> pf(n, nd), mf(n, nd), sf(n, nd);
    std::vector<std::uint8_t> m(n, 0);
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
        if (slip.toe_index[i] < 0) continue;                       // no admissible surface
        const auto it = single ? samples_by_region.end() : samples_by_region.find(region[i]);
        if (!single && it == samples_by_region.end())
            continue;
        const std::vector<rupture_parameter_sample>& smp = single ? *single : it->second;
        const slip_section sec = build_slip_section(flow, geometry, channel, ii, water_table_depth);
        std::vector<double> F(smp.size());
        double p_fail = 0.0, mean = 0.0;
        for (std::size_t s = 0; s < smp.size(); ++s)
        {
            F[s] = evaluate_slip_surface(sec, slip.toe_index[i], slip.depth_h[i], slip.sagitta_s[i], smp[s].parameters, evaluation);
            if (F[s] <= 1.0) p_fail += smp[s].weight;
            mean += smp[s].weight * F[s];
        }
        double var = 0.0;
        for (std::size_t s = 0; s < smp.size(); ++s) var += smp[s].weight * (F[s] - mean) * (F[s] - mean);
        ++n_src;
        if (p_fail >= 1.0 - 1e-12) ++n_cert;
        if (p_fail > 0.0) ++n_poss;
        m[i] = 1;
        pf[i] = static_cast<float>(std::min(1.0, p_fail));
        mf[i] = static_cast<float>(mean);
        sf[i] = static_cast<float>(std::sqrt(var));
    }
    // Mean and SD only where evaluated.
    std::vector<std::uint8_t> m_stats(n, 0);
    for (std::size_t i = 0; i < n; ++i) m_stats[i] = m[i] && em[i] && ev[i] != 0.0f;
    rupture_probability_rasters out;
    out.failure_probability = lsdtt_raster_cpu(std::move(pf), std::vector<std::uint8_t>(m), meta);
    out.mean_safety_factor = lsdtt_raster_cpu(std::move(mf), std::vector<std::uint8_t>(m_stats), meta);
    out.sd_safety_factor = lsdtt_raster_cpu(std::move(sf), std::move(m_stats), meta);
    out.n_sources = n_src;
    out.n_certain_failure = n_cert;
    out.n_possible_failure = n_poss;
    out.max_samples = max_samples;
    return out;
}

slip_body_catalogues build_slip_body_catalogues(const slip_surface_rasters& slip, const lsdtt_raster_meta& meta)
{
    const int nr = meta.n_rows, nc = meta.n_cols;
    const std::size_t n = static_cast<std::size_t>(nr) * static_cast<std::size_t>(nc);
    const double b = meta.cell_size, area = b * b;
    slip_body_catalogues out;
    const std::vector<float>& fs = slip.safety_factor.get_data();

    // By toe: groups of unstable crowns.
    std::map<std::int64_t, std::vector<std::size_t>> groups;   // toe -> body indices
    for (std::size_t k = 0; k < slip.body_crown.size(); ++k)
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
        if (rec.n_cells == 0) { --id; continue; }
        rec.area_m2 = static_cast<double>(rec.n_cells) * area;
        rec.mean_depth_m = rec.volume_m3 / rec.area_m2;
        const cell_shape sh = shape_from_moments(mom, meta);
        rec.centroid_x = sh.centroid_x; rec.centroid_y = sh.centroid_y;
        rec.length_m = sh.length_m; rec.width_m = sh.width_m; rec.length_to_width = sh.length_to_width;
        out.by_toe.push_back(rec);
    }

    // Amalgamated: 8-connected components of the envelope.
    const std::vector<float>& env = slip.body_depth.get_data();
    const std::vector<std::uint8_t>& em = slip.body_depth.get_mask();
    // Crowns and toes per cell, for the component counts.
    std::vector<std::int64_t> crown_toe(n, -1);
    for (std::int64_t c : slip.body_crown) crown_toe[static_cast<std::size_t>(c)] = slip.toe[static_cast<std::size_t>(c)];
    out.label_amalgamated.assign(n, 0);
    std::vector<std::size_t> stack;
    std::int32_t cid = 0;
    for (std::size_t s0 = 0; s0 < n; ++s0)
    {
        if (!em[s0] || out.label_amalgamated[s0]) continue;
        slip_body_record rec;
        rec.id = ++cid;
        rec.min_safety_factor = k_inf;
        cell_shape_moments mom;
        double sum_f = 0.0;
        std::vector<std::int64_t> toes;
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
            if (crown_toe[c] >= 0)
            {
                ++rec.n_crowns;
                sum_f += fs[c];
                rec.min_safety_factor = std::min(rec.min_safety_factor, static_cast<double>(fs[c]));
                toes.push_back(crown_toe[c]);
            }
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
        std::sort(toes.begin(), toes.end());
        rec.n_toes = static_cast<std::int64_t>(std::unique(toes.begin(), toes.end()) - toes.begin());
        rec.mean_safety_factor = rec.n_crowns ? sum_f / static_cast<double>(rec.n_crowns) : std::numeric_limits<double>::quiet_NaN();
        if (!rec.n_crowns) rec.min_safety_factor = std::numeric_limits<double>::quiet_NaN();
        rec.area_m2 = static_cast<double>(rec.n_cells) * area;
        rec.mean_depth_m = rec.volume_m3 / rec.area_m2;
        const cell_shape sh = shape_from_moments(mom, meta);
        rec.centroid_x = sh.centroid_x; rec.centroid_y = sh.centroid_y;
        rec.length_m = sh.length_m; rec.width_m = sh.width_m; rec.length_to_width = sh.length_to_width;
        out.amalgamated.push_back(rec);
    }
    return out;
}

} // namespace lsdtt
