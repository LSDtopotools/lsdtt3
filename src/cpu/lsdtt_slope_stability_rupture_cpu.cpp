// lsdtt_slope_stability_rupture_cpu.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// Planar rupture surfaces after scaLr. Equations and conventions: see the header.

#include "lsdtt_slope_stability_rupture_cpu.hpp"

#include <algorithm>
#include <cmath>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>

#include "../utils/lsdtt_numeric_utils.hpp"

#ifdef _OPENMP
#include <omp.h>
#endif

namespace lsdtt
{

namespace
{
constexpr double k_inf = std::numeric_limits<double>::infinity();
}

void validate_rupture_parameters(const rupture_parameters& p)
{
    auto finite = [](double v) { return std::isfinite(v); };
    if (!finite(p.cohesion_pa) || p.cohesion_pa < 0.0)
        throw std::invalid_argument("rupture: cohesion must be finite and >= 0 (Pa)");
    if (!finite(p.friction_coefficient) || p.friction_coefficient < 0.0)
        throw std::invalid_argument("rupture: the friction coefficient must be finite and >= 0");
    if (!finite(p.density_kg_m3) || !(p.density_kg_m3 > 0.0))
        throw std::invalid_argument("rupture: density must be finite and > 0 (kg/m3)");
    if (!finite(p.gravity_m_s2) || !(p.gravity_m_s2 > 0.0))
        throw std::invalid_argument("rupture: gravity must be finite and > 0 (m/s2)");
    if (!finite(p.water_density_kg_m3) || p.water_density_kg_m3 < 0.0)
        throw std::invalid_argument("rupture: water density must be finite and >= 0 (kg/m3)");
}

double rupture_safety_factor(double h, double dx, double dz, const rupture_parameters& p)
{
    if (!(dx > 0.0) || !(h > 0.0) || !(h < dz)) return k_inf;
    const double s = dz - h;
    const double rho_g = p.density_kg_m3 * p.gravity_m_s2;
    return p.cohesion_pa / (rho_g * h) * (dx * dx + s * s) / (dx * s)
         + p.friction_coefficient * dx / s;
}

rupture_pair rupture_optimal_depth(double dx, double dz, const rupture_parameters& p)
{
    rupture_pair out;
    if (!(dx > 0.0) || !(dz > 0.0)) return out;
    out.admissible = true;
    const double c = p.cohesion_pa;
    if (c == 0.0)
    {
        // Cohesionless: F decreases towards h = 0, the chord's infinite slope.
        out.depth = 0.0;
        out.safety_factor = p.friction_coefficient * dx / dz;
        return out;
    }
    const double rho_g = p.density_kg_m3 * p.gravity_m_s2;
    const double L2 = dx * dx + dz * dz;
    const double root = std::sqrt(c * (c + p.friction_coefficient * rho_g * dz) * L2);
    out.depth = c * dz * L2 / (c * L2 + dx * root);
    out.safety_factor = rupture_safety_factor(out.depth, dx, dz, p);
    return out;
}

double rupture_safety_factor_wet(double h, double dx, double dz, double d,
                                 const rupture_parameters& p)
{
    const double dry = rupture_safety_factor(h, dx, dz, p);
    if (!std::isfinite(dry) || !std::isfinite(d) || h <= d) return dry;
    const double s = dz - h;
    return dry - p.friction_coefficient * (p.water_density_kg_m3 / p.density_kg_m3)
                 * ((h - d) / h) * (dx * dx + s * s) / (dx * s);
}

rupture_pair rupture_optimal_depth_wet(double dx, double dz, double d, const rupture_parameters& p)
{
    if (!std::isfinite(d) || d >= dz) return rupture_optimal_depth(dx, dz, p);
    rupture_pair out;
    if (!(dx > 0.0) || !(dz > 0.0)) return out;
    out.admissible = true;
    auto F = [&](double h) { return rupture_safety_factor_wet(h, dx, dz, d, p); };

    // Dry part, (0, d]: the dry optimum if it lies above the water table,
    // otherwise the boundary h = d (F is continuous there).
    if (d > 0.0)
    {
        const rupture_pair dry = rupture_optimal_depth(dx, dz, p);
        const double h = (dry.depth > 0.0 && dry.depth < d) ? dry.depth : d;
        const double f = (dry.depth > 0.0 && dry.depth < d) ? dry.safety_factor : F(d);
        if (f < out.safety_factor) { out.safety_factor = f; out.depth = h; }
    }
    // Wet part, (max(d, 0), dz): scan, then refine around the best scan point.
    const double lo = std::max(d, 0.0), hi = dz;
    const int n_scan = 32;
    const double step = (hi - lo) / n_scan;
    double best_f = std::numeric_limits<double>::infinity();
    int best_k = -1;
    for (int k = 1; k < n_scan; ++k)
    {
        const double f = F(lo + k * step);
        if (f < best_f) { best_f = f; best_k = k; }
    }
    if (best_k > 0)
    {
        const double a = lo + (best_k - 1) * step, b = lo + (best_k + 1) * step;
        const brent_minimum m = minimise_brent_bounded(F, std::max(a, lo + 1e-12 * dz), std::min(b, hi * (1.0 - 1e-12)),
                                                       1e-9 * dz, 200);
        const double f = (m.f < best_f) ? m.f : best_f;
        const double h = (m.f < best_f) ? m.x : lo + best_k * step;
        if (f < out.safety_factor) { out.safety_factor = f; out.depth = h; }
    }
    return out;
}

rupture_rasters compute_rupture_planes(
    const lsdtt_flow_base& flow,
    const lsdtt_raster_cpu& geometry,
    const std::vector<std::uint8_t>& channel,
    const std::vector<std::int32_t>& region,
    const std::map<std::int32_t, rupture_parameters>& by_region,
    int n_threads,
    const std::vector<double>& water_table_depth)
{
    const lsdtt_raster_meta& meta = geometry.get_meta();
    const int nr = meta.n_rows, nc = meta.n_cols;
    const std::size_t n = static_cast<std::size_t>(nr) * static_cast<std::size_t>(nc);
    if (flow.n_rows() != nr || flow.n_cols() != nc)
        throw std::invalid_argument("compute_rupture_planes: the flow and geometry grids differ");
    if (channel.size() != n)
        throw std::invalid_argument("compute_rupture_planes: the channel mask differs from the grid");
    if (!region.empty() && region.size() != n)
        throw std::invalid_argument("compute_rupture_planes: the region grid differs from the grid");
    if (by_region.empty())
        throw std::invalid_argument("compute_rupture_planes: no parameter set given");
    if (!water_table_depth.empty() && water_table_depth.size() != n)
        throw std::invalid_argument("compute_rupture_planes: the water-table grid differs from the grid");
    if (region.empty() && by_region.size() != 1)
        throw std::invalid_argument("compute_rupture_planes: without a region grid exactly one parameter set is used");
    for (const auto& kv : by_region)
    {
        try { validate_rupture_parameters(kv.second); }
        catch (const std::invalid_argument& e)
        {
            throw std::invalid_argument(std::string(e.what()) + " (region " + std::to_string(kv.first) + ")");
        }
    }

    const std::vector<float>& z = geometry.get_data();
    const std::vector<std::uint8_t>& zm = geometry.get_mask();
    auto is_source = [&](std::size_t i)
    {
        return zm[i] && flow.is_valid_cell(static_cast<std::int64_t>(i)) && !channel[i];
    };

    if (!region.empty())
    {
        std::set<std::int32_t> missing;
        for (std::size_t i = 0; i < n; ++i)
            if (is_source(i) && region[i] != k_rupture_region_nodata && !by_region.count(region[i]))
                missing.insert(region[i]);
        if (!missing.empty())
        {
            std::ostringstream msg;
            msg << "compute_rupture_planes: region ids with no parameter row:";
            int shown = 0;
            for (std::int32_t id : missing)
            {
                if (shown++ == 20) { msg << " ..."; break; }
                msg << ' ' << id;
            }
            throw std::invalid_argument(msg.str());
        }
    }

    const float nd = static_cast<float>(meta.no_data_value);
    std::vector<float> fs(n, nd), dep(n, nd), srup(n, nd), stopo(n, nd), len(n, nd), uns(n, nd);
    std::vector<std::uint8_t> m(n, 0);
    std::vector<std::int64_t> daylight(n, -1);
    const double b = meta.cell_size;
    const rupture_parameters* single = region.empty() ? &by_region.begin()->second : nullptr;

    std::int64_t n_sources = 0, n_unstable = 0, n_channel = 0;
    const std::int64_t max_steps = static_cast<std::int64_t>(n);   // guard against a cyclic graph

#ifdef _OPENMP
    const int threads = (n_threads > 0) ? n_threads : omp_get_max_threads();
#pragma omp parallel for schedule(dynamic, 64) num_threads(threads) reduction(+ : n_sources, n_unstable, n_channel)
#else
    (void)n_threads;
#endif
    for (std::int64_t ii = 0; ii < static_cast<std::int64_t>(n); ++ii)
    {
        const std::size_t i = static_cast<std::size_t>(ii);
        if (zm[i] && channel[i]) { ++n_channel; continue; }
        if (!is_source(i)) continue;
        if (!single && region[i] == k_rupture_region_nodata) continue;
        const rupture_parameters& p = single ? *single : by_region.at(region[i]);

        const int ri = static_cast<int>(ii / nc), ci = static_cast<int>(ii % nc);
        const double zi = static_cast<double>(z[i]);
        const double dw = water_table_depth.empty() ? std::numeric_limits<double>::quiet_NaN()
                                                    : water_table_depth[i];
        double best_f = k_inf, best_h = 0.0, best_dx = 0.0, best_dz = 0.0;
        std::int64_t best_d = -1;

        std::int64_t pflat = flow.receiver_of_unchecked(ii);
        for (std::int64_t step = 0; pflat >= 0 && step < max_steps; ++step)
        {
            const std::size_t pi = static_cast<std::size_t>(pflat);
            if (channel[pi] || !zm[pi]) break;     // walk stays on the hillslope
            const int rp = static_cast<int>(pflat / nc), cp = static_cast<int>(pflat % nc);
            const double dx = b * std::hypot(static_cast<double>(rp - ri), static_cast<double>(cp - ci));
            const double dz = zi - static_cast<double>(z[pi]);
            const rupture_pair rp_res = rupture_optimal_depth_wet(dx, dz, dw, p);
            if (rp_res.admissible && rp_res.safety_factor < best_f)
            {
                best_f = rp_res.safety_factor;
                best_h = rp_res.depth;
                best_dx = dx;
                best_dz = dz;
                best_d = pflat;
            }
            pflat = flow.receiver_of_unchecked(pflat);
        }
        if (best_d < 0) continue;

        ++n_sources;
        m[i] = 1;
        fs[i]    = static_cast<float>(best_f);
        dep[i]   = static_cast<float>(best_h);
        srup[i]  = static_cast<float>((best_dz - best_h) / best_dx);
        stopo[i] = static_cast<float>(best_dz / best_dx);
        len[i]   = static_cast<float>(best_dx);
        uns[i]   = (best_f <= 1.0) ? 1.0f : 0.0f;
        if (best_f <= 1.0) ++n_unstable;
        daylight[i] = best_d;
    }

    auto make = [&](std::vector<float>&& v)
    {
        return lsdtt_raster_cpu(std::move(v), std::vector<std::uint8_t>(m), meta);
    };
    rupture_rasters out;
    out.safety_factor     = make(std::move(fs));
    out.depth             = make(std::move(dep));
    out.rupture_slope     = make(std::move(srup));
    out.topographic_slope = make(std::move(stopo));
    out.length            = make(std::move(len));
    out.unstable          = make(std::move(uns));
    out.daylight          = std::move(daylight);
    out.n_sources = n_sources;
    out.n_unstable = n_unstable;
    out.n_channel = n_channel;
    return out;
}

} // namespace lsdtt
