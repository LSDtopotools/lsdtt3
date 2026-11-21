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

double rupture_wet_lower_bound(double dx, double dz, double d, const rupture_parameters& p)
{
    // F_wet = (A + mu) dx / s + A s / dx on the wet part, s = dz - h, with
    // A(h) = k / h - mu r, k = c / (rho g) + mu r d, r = rho_w / rho.
    // With k >= 0, A >= A_min = k / dz - mu r; s < dz - lo. Header.
    const double lo = std::max(d, 0.0), smax = dz - lo;
    if (!(dx > 0.0) || !(smax > 0.0)) return k_inf;
    const double mu = p.friction_coefficient, r = p.water_density_kg_m3 / p.density_kg_m3;
    const double k = p.cohesion_pa / (p.density_kg_m3 * p.gravity_m_s2) + mu * r * d;
    if (k < 0.0) return -k_inf;
    const double a_min = k / dz - mu * r;
    if (a_min + mu < 0.0) return -k_inf;
    return (a_min + mu) * dx / smax + std::min(a_min, 0.0) * smax / dx;
}

rupture_pair rupture_optimal_depth_wet(double dx, double dz, double d, const rupture_parameters& p,
                                       double cutoff)
{
    if (!std::isfinite(d) || d >= dz) return rupture_optimal_depth(dx, dz, p);
    rupture_pair out;
    if (!(dx > 0.0) || !(dz > 0.0)) return out;
    out.admissible = true;
    auto F = [&](double h) { return rupture_safety_factor_wet(h, dx, dz, d, p); };

    // Dry part, (0, d]: the dry optimum if it lies above the water table,
    // otherwise the boundary h = d (F is continuous there). With c = 0 the dry
    // optimum is the limit h -> 0 (depth 0, F = mu dx / dz), which counts.
    if (d > 0.0)
    {
        const rupture_pair dry = rupture_optimal_depth(dx, dz, p);
        const bool inside = dry.depth >= 0.0 && dry.depth < d;
        const double h = inside ? dry.depth : d;
        const double f = inside ? dry.safety_factor : F(d);
        if (f < out.safety_factor) { out.safety_factor = f; out.depth = h; }
    }
    // Wet part, (max(d, 0), dz): skipped when its lower bound cannot beat the
    // dry part or the caller's cutoff; otherwise Brent alone, F_wet being
    // unimodal there (header).
    if (rupture_wet_lower_bound(dx, dz, d, p) >= std::min(cutoff, out.safety_factor)) return out;
    const double lo = std::max(d, 0.0), hi = dz;
    const brent_minimum m = minimise_brent_bounded(F, lo + 1e-12 * dz, hi * (1.0 - 1e-12), 1e-7 * dz, 200);
    if (m.f < out.safety_factor) { out.safety_factor = m.f; out.depth = m.x; }
    return out;
}

namespace
{

// Checks shared by compute_rupture_planes and compute_rupture_failure_probability.
// has_id(id): a parameter row / sample list exists for region id.
template <class HasId>
void check_rupture_inputs(const lsdtt_flow_base& flow, const lsdtt_raster_cpu& geometry,
                          const std::vector<std::uint8_t>& channel, const std::vector<std::int32_t>& region,
                          std::size_t n_sets, const std::vector<double>& water_table_depth, HasId has_id,
                          const char* fn)
{
    const lsdtt_raster_meta& meta = geometry.get_meta();
    const int nr = meta.n_rows, nc = meta.n_cols;
    const std::size_t n = static_cast<std::size_t>(nr) * static_cast<std::size_t>(nc);
    const std::string f(fn);
    if (flow.n_rows() != nr || flow.n_cols() != nc)
        throw std::invalid_argument(f + ": the flow and geometry grids differ");
    if (channel.size() != n)
        throw std::invalid_argument(f + ": the channel mask differs from the grid");
    if (!region.empty() && region.size() != n)
        throw std::invalid_argument(f + ": the region grid differs from the grid");
    if (n_sets == 0)
        throw std::invalid_argument(f + ": no parameter set given");
    if (!water_table_depth.empty() && water_table_depth.size() != n)
        throw std::invalid_argument(f + ": the water-table grid differs from the grid");
    if (region.empty() && n_sets != 1)
        throw std::invalid_argument(f + ": without a region grid exactly one parameter set is used");
    if (region.empty()) return;
    const std::vector<std::uint8_t>& zm = geometry.get_mask();
    std::set<std::int32_t> missing;
    for (std::size_t i = 0; i < n; ++i)
        if (zm[i] && flow.is_valid_cell(static_cast<std::int64_t>(i)) && !channel[i]
            && region[i] != k_rupture_region_nodata && !has_id(region[i]))
            missing.insert(region[i]);
    if (!missing.empty())
    {
        std::ostringstream msg;
        msg << f << ": region ids with no parameter row:";
        int shown = 0;
        for (std::int32_t id : missing)
        {
            if (shown++ == 20) { msg << " ..."; break; }
            msg << ' ' << id;
        }
        throw std::invalid_argument(msg.str());
    }
}

// Visit every daylight candidate of source i: visit(dx, dz, flat index of D).
// The walk follows receivers and stops before a channel or nodata cell.
template <class Visit>
void walk_daylight_candidates(const lsdtt_flow_base& flow, const std::vector<float>& z,
                              const std::vector<std::uint8_t>& zm, const std::vector<std::uint8_t>& channel,
                              int nc, double b, std::int64_t ii, std::int64_t max_steps, Visit visit)
{
    const int ri = static_cast<int>(ii / nc), ci = static_cast<int>(ii % nc);
    const double zi = static_cast<double>(z[static_cast<std::size_t>(ii)]);
    std::int64_t pflat = flow.receiver_of_unchecked(ii);
    for (std::int64_t step = 0; pflat >= 0 && step < max_steps; ++step)
    {
        const std::size_t pi = static_cast<std::size_t>(pflat);
        if (channel[pi] || !zm[pi]) break;     // walk stays on the hillslope
        const int rp = static_cast<int>(pflat / nc), cp = static_cast<int>(pflat % nc);
        const double dx = b * std::hypot(static_cast<double>(rp - ri), static_cast<double>(cp - ci));
        const double dz = zi - static_cast<double>(z[pi]);
        visit(dx, dz, pflat);
        pflat = flow.receiver_of_unchecked(pflat);
    }
}

} // namespace

rupture_rasters compute_rupture_planes(
    const lsdtt_flow_base& flow,
    const lsdtt_raster_cpu& geometry,
    const std::vector<std::uint8_t>& channel,
    const std::vector<std::int32_t>& region,
    const std::map<std::int32_t, rupture_parameters>& by_region,
    int n_threads,
    const std::vector<double>& water_table_depth)
{
    check_rupture_inputs(flow, geometry, channel, region, by_region.size(), water_table_depth,
                         [&](std::int32_t id) { return by_region.count(id) > 0; }, "compute_rupture_planes");
    for (const auto& kv : by_region)
    {
        try { validate_rupture_parameters(kv.second); }
        catch (const std::invalid_argument& e)
        {
            throw std::invalid_argument(std::string(e.what()) + " (region " + std::to_string(kv.first) + ")");
        }
    }
    const lsdtt_raster_meta& meta = geometry.get_meta();
    const int nc = meta.n_cols;
    const std::size_t n = static_cast<std::size_t>(meta.n_rows) * static_cast<std::size_t>(nc);
    const std::vector<float>& z = geometry.get_data();
    const std::vector<std::uint8_t>& zm = geometry.get_mask();
    auto is_source = [&](std::size_t i)
    {
        return zm[i] && flow.is_valid_cell(static_cast<std::int64_t>(i)) && !channel[i];
    };

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
        const double dw = water_table_depth.empty() ? std::numeric_limits<double>::quiet_NaN()
                                                    : water_table_depth[i];
        double best_f = k_inf, best_h = 0.0, best_dx = 0.0, best_dz = 0.0;
        std::int64_t best_d = -1;
        walk_daylight_candidates(flow, z, zm, channel, nc, b, ii, max_steps,
            [&](double dx, double dz, std::int64_t pflat) {
                // best_f as the cutoff: a pair whose wet bound cannot beat it skips
                // its Brent search (the result is then not used either way).
                const rupture_pair rp_res = rupture_optimal_depth_wet(dx, dz, dw, p, best_f);
                if (rp_res.admissible && rp_res.safety_factor < best_f)
                {
                    best_f = rp_res.safety_factor;
                    best_h = rp_res.depth;
                    best_dx = dx;
                    best_dz = dz;
                    best_d = pflat;
                }
            });
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

std::vector<rupture_parameter_sample> rupture_parameter_samples(
    const rupture_parameters& central,
    const equal_density_distribution& cohesion, int n_cohesion,
    const equal_density_distribution& friction, int n_friction)
{
    const equal_density_sample sc = equal_density_samples(cohesion, n_cohesion);
    const equal_density_sample sm = equal_density_samples(friction, n_friction);
    std::vector<rupture_parameter_sample> out;
    out.reserve(sc.values.size() * sm.values.size());
    for (std::size_t a = 0; a < sc.values.size(); ++a)
        for (std::size_t k = 0; k < sm.values.size(); ++k)
        {
            rupture_parameter_sample s;
            s.parameters = central;
            s.parameters.cohesion_pa = sc.values[a];
            s.parameters.friction_coefficient = sm.values[k];
            s.weight = sc.weights[a] * sm.weights[k];
            validate_rupture_parameters(s.parameters);
            out.push_back(s);
        }
    return out;
}

rupture_probability_rasters compute_rupture_failure_probability(
    const lsdtt_flow_base& flow,
    const lsdtt_raster_cpu& geometry,
    const std::vector<std::uint8_t>& channel,
    const std::vector<std::int32_t>& region,
    const std::map<std::int32_t, std::vector<rupture_parameter_sample>>& samples_by_region,
    int n_threads,
    const std::vector<double>& water_table_depth)
{
    check_rupture_inputs(flow, geometry, channel, region, samples_by_region.size(), water_table_depth,
                         [&](std::int32_t id) { return samples_by_region.count(id) > 0; },
                         "compute_rupture_failure_probability");
    std::size_t max_samples = 0;
    for (const auto& kv : samples_by_region)
    {
        const std::string where = " (region " + std::to_string(kv.first) + ")";
        if (kv.second.empty())
            throw std::invalid_argument("compute_rupture_failure_probability: empty sample list" + where);
        double wsum = 0.0;
        for (const rupture_parameter_sample& smp : kv.second)
        {
            try { validate_rupture_parameters(smp.parameters); }
            catch (const std::invalid_argument& e) { throw std::invalid_argument(std::string(e.what()) + where); }
            if (!std::isfinite(smp.weight) || smp.weight < 0.0)
                throw std::invalid_argument("compute_rupture_failure_probability: negative or non-finite weight" + where);
            wsum += smp.weight;
        }
        if (std::fabs(wsum - 1.0) > 1e-9)
            throw std::invalid_argument("compute_rupture_failure_probability: weights do not sum to 1" + where);
        max_samples = std::max(max_samples, kv.second.size());
    }

    const lsdtt_raster_meta& meta = geometry.get_meta();
    const int nc = meta.n_cols;
    const std::size_t n = static_cast<std::size_t>(meta.n_rows) * static_cast<std::size_t>(nc);
    const std::vector<float>& z = geometry.get_data();
    const std::vector<std::uint8_t>& zm = geometry.get_mask();
    const float nd = static_cast<float>(meta.no_data_value);
    std::vector<float> pf(n, nd), mean_f(n, nd), sd_f(n, nd);
    std::vector<std::uint8_t> m(n, 0);
    const double b = meta.cell_size;
    const std::vector<rupture_parameter_sample>* single = region.empty() ? &samples_by_region.begin()->second : nullptr;
    const std::int64_t max_steps = static_cast<std::int64_t>(n);
    std::int64_t n_sources = 0, n_certain = 0, n_possible = 0;

#ifdef _OPENMP
    const int threads = (n_threads > 0) ? n_threads : omp_get_max_threads();
#pragma omp parallel num_threads(threads) reduction(+ : n_sources, n_certain, n_possible)
#else
    (void)n_threads;
#endif
    {
        std::vector<double> best(max_samples);
#ifdef _OPENMP
#pragma omp for schedule(dynamic, 64)
#endif
        for (std::int64_t ii = 0; ii < static_cast<std::int64_t>(n); ++ii)
        {
            const std::size_t i = static_cast<std::size_t>(ii);
            if (!zm[i] || channel[i] || !flow.is_valid_cell(ii)) continue;
            if (!single && region[i] == k_rupture_region_nodata) continue;
            const std::vector<rupture_parameter_sample>& smp = single ? *single : samples_by_region.at(region[i]);
            const std::size_t ns = smp.size();
            const double dw = water_table_depth.empty() ? std::numeric_limits<double>::quiet_NaN()
                                                        : water_table_depth[i];
            std::fill(best.begin(), best.begin() + static_cast<std::ptrdiff_t>(ns), k_inf);
            bool any = false;
            walk_daylight_candidates(flow, z, zm, channel, nc, b, ii, max_steps,
                [&](double dx, double dz, std::int64_t) {
                    if (!(dx > 0.0) || !(dz > 0.0)) return;   // inadmissible for every sample
                    any = true;
                    for (std::size_t s = 0; s < ns; ++s)
                    {
                        const rupture_pair r = rupture_optimal_depth_wet(dx, dz, dw, smp[s].parameters, best[s]);
                        if (r.safety_factor < best[s]) best[s] = r.safety_factor;
                    }
                });
            if (!any) continue;
            double p_fail = 0.0, mean = 0.0;
            for (std::size_t s = 0; s < ns; ++s)
            {
                if (best[s] <= 1.0) p_fail += smp[s].weight;
                mean += smp[s].weight * best[s];
            }
            double var = 0.0;
            for (std::size_t s = 0; s < ns; ++s)
                var += smp[s].weight * (best[s] - mean) * (best[s] - mean);
            ++n_sources;
            if (p_fail >= 1.0 - 1e-12) ++n_certain;
            if (p_fail > 0.0) ++n_possible;
            m[i] = 1;
            pf[i] = static_cast<float>(std::min(1.0, p_fail));
            mean_f[i] = static_cast<float>(mean);
            sd_f[i] = static_cast<float>(std::sqrt(var));
        }
    }

    auto make = [&](std::vector<float>&& v)
    {
        return lsdtt_raster_cpu(std::move(v), std::vector<std::uint8_t>(m), meta);
    };
    rupture_probability_rasters out;
    out.failure_probability = make(std::move(pf));
    out.mean_safety_factor  = make(std::move(mean_f));
    out.sd_safety_factor    = make(std::move(sd_f));
    out.n_sources = n_sources;
    out.n_certain_failure = n_certain;
    out.n_possible_failure = n_possible;
    out.max_samples = max_samples;
    return out;
}

} // namespace lsdtt
