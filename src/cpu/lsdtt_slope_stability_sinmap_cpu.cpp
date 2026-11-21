// lsdtt_slope_stability_sinmap_cpu.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// SINMAP stability index. Equations and conventions: see the header.

#include "lsdtt_slope_stability_sinmap_cpu.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <set>
#include <sstream>
#include <stdexcept>

#include "../utils/lsdtt_distributions.hpp"

#ifdef _OPENMP
#include <omp.h>
#endif

namespace lsdtt
{

namespace
{

constexpr double k_inf = std::numeric_limits<double>::infinity();

void check_bounds(double lo, double hi, const char* name, bool positive)
{
    if (!std::isfinite(lo) || !std::isfinite(hi))
        throw std::invalid_argument(std::string("SINMAP: ") + name + " bounds must be finite");
    if (lo > hi)
        throw std::invalid_argument(std::string("SINMAP: ") + name + " minimum exceeds its maximum");
    if (positive ? !(lo > 0.0) : (lo < 0.0))
        throw std::invalid_argument(std::string("SINMAP: ") + name
                                    + (positive ? " must be positive" : " must not be negative"));
}

} // namespace

void validate_sinmap_parameters(const sinmap_parameters& p)
{
    check_bounds(p.dimensionless_cohesion_min, p.dimensionless_cohesion_max,
                 "dimensionless cohesion", false);
    check_bounds(p.friction_coefficient_min, p.friction_coefficient_max,
                 "friction coefficient", false);
    check_bounds(p.transmissivity_over_recharge_min_metres,
                 p.transmissivity_over_recharge_max_metres,
                 "transmissivity over recharge", true);
    if (!(p.water_to_soil_density_ratio >= 0.0 && p.water_to_soil_density_ratio <= 1.0))
        throw std::invalid_argument("SINMAP: the water-to-soil density ratio must be in [0, 1]");
}

int sinmap_stability_class(double si)
{
    if (si > 1.5)  return 1;
    if (si > 1.25) return 2;
    if (si > 1.0)  return 3;
    if (si > 0.5)  return 4;
    if (si > 0.0)  return 5;
    return 6;
}

sinmap_cell_result sinmap_evaluate_cell(
    double slope_gradient,
    double a,
    const sinmap_parameters& p)
{
    sinmap_cell_result out;
    if (!(slope_gradient > 0.0))
    {
        out.flat = true;
        out.fs_min = out.fs_max = k_inf;
        out.probability_stable = 1.0;
        out.stability_index = k_inf;
        out.stability_class = 1;
        out.saturation_region = k_sinmap_always_saturated;
        out.saturation = 3.0;
        return out;
    }

    const double theta = std::atan(slope_gradient);
    const double s = std::sin(theta);
    const double c = std::cos(theta);
    const double r = p.water_to_soil_density_ratio;
    const double c1 = p.dimensionless_cohesion_min, c2 = p.dimensionless_cohesion_max;
    const double t1 = p.friction_coefficient_min,   t2 = p.friction_coefficient_max;
    const double x1 = 1.0 / p.transmissivity_over_recharge_max_metres;   // R/T, smallest
    const double x2 = 1.0 / p.transmissivity_over_recharge_min_metres;   // R/T, largest

    // Wetness at the two ends of the R/T range.
    const double q1 = x1 * a / s;
    const double q2 = x2 * a / s;

    out.fs_min = (c1 + c * (1.0 - r * std::min(q2, 1.0)) * t1) / s;
    out.fs_max = (c2 + c * (1.0 - r * std::min(q1, 1.0)) * t2) / s;

    if (q1 >= 1.0)
    {
        out.saturation_region = k_sinmap_always_saturated;
        out.saturation = 3.0;
    }
    else if (q2 > 1.0)
    {
        out.saturation_region = k_sinmap_threshold_saturated;
        out.saturation = 2.0;
    }
    else
    {
        out.saturation_region = k_sinmap_never_saturated;
        out.saturation = q2;
    }

    // P(FS <= 1), FS = X + Y B with X = C/s, B = t, Y = (c/s)(1 - r w).
    const double k = c / s;
    const double xl = c1 / s, xu = c2 / s;
    const double y_saturated = k * (1.0 - r);
    double p_fail = 0.0;
    if (out.saturation_region == k_sinmap_always_saturated)
    {
        p_fail = cdf_sum_with_product_of_uniforms(xl, xu, y_saturated, y_saturated, t1, t2, 1.0);
    }
    else if (out.saturation_region == k_sinmap_never_saturated)
    {
        // w = x a / s is uniform on [q1, q2]; Y decreases with w.
        p_fail = cdf_sum_with_product_of_uniforms(
            xl, xu, k * (1.0 - r * q2), k * (1.0 - r * q1), t1, t2, 1.0);
    }
    else
    {
        // w = 1 with probability (q2 - 1)/(q2 - q1); otherwise uniform on [q1, 1].
        const double p_sat = (q2 - 1.0) / (q2 - q1);
        const double f_sat = cdf_sum_with_product_of_uniforms(
            xl, xu, y_saturated, y_saturated, t1, t2, 1.0);
        const double f_unsat = cdf_sum_with_product_of_uniforms(
            xl, xu, y_saturated, k * (1.0 - r * q1), t1, t2, 1.0);
        p_fail = p_sat * f_sat + (1.0 - p_sat) * f_unsat;
    }
    out.probability_stable = std::min(1.0, std::max(0.0, 1.0 - p_fail));

    if (out.fs_min > 1.0)      out.stability_index = out.fs_min;
    else if (out.fs_max < 1.0) out.stability_index = 0.0;
    else                       out.stability_index = out.probability_stable;
    out.stability_class = sinmap_stability_class(out.stability_index);
    return out;
}

// -----------------------------------------------------------------------------
// Raster pass
// -----------------------------------------------------------------------------
sinmap_rasters compute_sinmap(
    const lsdtt_raster_cpu& slope,
    const lsdtt_raster_cpu& drainage_area,
    const std::vector<std::int32_t>& region,
    const std::map<std::int32_t, sinmap_parameters>& by_region,
    const std::vector<std::uint8_t>& exclude,
    int n_threads,
    const std::vector<double>& contour_width)
{
    const lsdtt_raster_meta& meta = slope.get_meta();
    const lsdtt_raster_meta& am = drainage_area.get_meta();
    if (meta.n_rows != am.n_rows || meta.n_cols != am.n_cols)
        throw std::invalid_argument("compute_sinmap: slope and drainage area grids differ");
    const std::size_t n = static_cast<std::size_t>(meta.n_rows) * static_cast<std::size_t>(meta.n_cols);
    if (!region.empty() && region.size() != n)
        throw std::invalid_argument("compute_sinmap: the region grid differs from the slope grid");
    if (!exclude.empty() && exclude.size() != n)
        throw std::invalid_argument("compute_sinmap: the exclusion grid differs from the slope grid");
    if (!contour_width.empty() && contour_width.size() != n)
        throw std::invalid_argument("compute_sinmap: the contour-width grid differs from the slope grid");
    if (by_region.empty())
        throw std::invalid_argument("compute_sinmap: no parameter set given");
    if (region.empty() && by_region.size() != 1)
        throw std::invalid_argument(
            "compute_sinmap: without a region grid exactly one parameter set is used");
    for (const auto& kv : by_region)
    {
        try { validate_sinmap_parameters(kv.second); }
        catch (const std::invalid_argument& e)
        {
            throw std::invalid_argument(std::string(e.what()) + " (region "
                                        + std::to_string(kv.first) + ")");
        }
    }

    const std::vector<float>& sv = slope.get_data();
    const std::vector<std::uint8_t>& sm = slope.get_mask();
    const std::vector<float>& av = drainage_area.get_data();
    const std::vector<std::uint8_t>& amask = drainage_area.get_mask();

    // Region ids on valid cells must all have a parameter set.
    if (!region.empty())
    {
        std::set<std::int32_t> missing;
        for (std::size_t i = 0; i < n; ++i)
            if (sm[i] && amask[i] && region[i] != k_sinmap_region_nodata
                && by_region.find(region[i]) == by_region.end())
                missing.insert(region[i]);
        if (!missing.empty())
        {
            std::ostringstream msg;
            msg << "compute_sinmap: region ids with no parameter row:";
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
    auto blank = [&]() { return std::vector<float>(n, nd); };
    std::vector<float> si = blank(), fmin = blank(), fmax = blank(), ps = blank(),
                       cls = blank(), sat = blank(), sca = blank();
    std::vector<std::uint8_t> m_cont(n, 0), m_all(n, 0);   // continuous (no flats) / every answer

    const sinmap_parameters* single = region.empty() ? &by_region.begin()->second : nullptr;
    const double b = meta.cell_size;

    std::int64_t n_eval = 0, n_flat = 0, n_excl = 0;

#ifdef _OPENMP
    const int threads = (n_threads > 0) ? n_threads : omp_get_max_threads();
#pragma omp parallel for schedule(static) num_threads(threads) reduction(+ : n_eval, n_flat, n_excl)
#else
    (void)n_threads;
#endif
    for (int row = 0; row < meta.n_rows; ++row)
    {
        for (int col = 0; col < meta.n_cols; ++col)
        {
            const std::size_t i = static_cast<std::size_t>(row) * static_cast<std::size_t>(meta.n_cols)
                                + static_cast<std::size_t>(col);
            if (!sm[i] || !amask[i]) continue;
            const double S = static_cast<double>(sv[i]);
            const double A = static_cast<double>(av[i]);
            if (!std::isfinite(S) || !std::isfinite(A) || !(A > 0.0)) continue;
            if (!single && region[i] == k_sinmap_region_nodata) continue;
            if (!exclude.empty() && exclude[i]) { ++n_excl; continue; }

            const sinmap_parameters& p = single ? *single : by_region.at(region[i]);
            const double w = contour_width.empty() ? b : contour_width[i];
            if (!(w > 0.0) || !std::isfinite(w)) continue;
            const double a = A / w;
            const sinmap_cell_result cell = sinmap_evaluate_cell(S, a, p);

            ++n_eval;
            m_all[i] = 1;
            ps[i]  = static_cast<float>(cell.probability_stable);
            cls[i] = static_cast<float>(cell.stability_class);
            sat[i] = static_cast<float>(cell.saturation);
            sca[i] = static_cast<float>(a);
            if (cell.flat) { ++n_flat; continue; }
            m_cont[i] = 1;
            si[i]   = static_cast<float>(cell.stability_index);
            fmin[i] = static_cast<float>(cell.fs_min);
            fmax[i] = static_cast<float>(cell.fs_max);
        }
    }

    auto make = [&](std::vector<float>&& v, const std::vector<std::uint8_t>& m)
    {
        return lsdtt_raster_cpu(std::move(v), std::vector<std::uint8_t>(m), meta);
    };
    sinmap_rasters out;
    out.stability_index         = make(std::move(si), m_cont);
    out.fs_min                  = make(std::move(fmin), m_cont);
    out.fs_max                  = make(std::move(fmax), m_cont);
    out.probability_stable      = make(std::move(ps), m_all);
    out.stability_class         = make(std::move(cls), m_all);
    out.saturation              = make(std::move(sat), m_all);
    out.specific_catchment_area = make(std::move(sca), m_all);
    out.n_evaluated = n_eval;
    out.n_flat = n_flat;
    out.n_excluded = n_excl;
    return out;
}

// -----------------------------------------------------------------------------
// Contour width
// -----------------------------------------------------------------------------
std::vector<double> dinf_contour_width(const lsdtt_flowgraph_cpu& graph)
{
    const lsdtt_raster_result dir = graph.dinf_flow_direction();
    const double b = dir.meta.cell_size;
    const double deg = 3.14159265358979323846 / 180.0;
    std::vector<double> w(dir.data.size(), b);
    for (std::size_t i = 0; i < w.size(); ++i)
    {
        if (!dir.mask[i]) continue;
        const double beta = static_cast<double>(dir.data[i]) * deg;
        w[i] = b * std::max(std::fabs(std::sin(beta)), std::fabs(std::cos(beta)));
    }
    return w;
}

// -----------------------------------------------------------------------------
// Edge contamination
// -----------------------------------------------------------------------------
std::vector<std::uint8_t> dinf_edge_contaminated_mask(
    const lsdtt_flowgraph_cpu& graph,
    const lsdtt_raster_meta& meta,
    const std::vector<std::uint8_t>& valid)
{
    const int nr = meta.n_rows, nc = meta.n_cols;
    const std::size_t n = static_cast<std::size_t>(nr) * static_cast<std::size_t>(nc);
    if (valid.size() != n)
        throw std::invalid_argument("dinf_edge_contaminated_mask: mask size differs from the grid");

    // Weight 1 on every valid cell that is on the raster edge or has a nodata
    // 8-neighbour; the accumulated weight is positive exactly where the
    // contributing area contains such a cell.
    std::vector<double> weight(n, 0.0);
    for (int row = 0; row < nr; ++row)
        for (int col = 0; col < nc; ++col)
        {
            const std::size_t i = static_cast<std::size_t>(row) * nc + col;
            if (!valid[i]) continue;
            bool edge = (row == 0 || col == 0 || row == nr - 1 || col == nc - 1);
            for (int dr = -1; dr <= 1 && !edge; ++dr)
                for (int dc = -1; dc <= 1 && !edge; ++dc)
                {
                    if (dr == 0 && dc == 0) continue;
                    const std::size_t j = static_cast<std::size_t>(row + dr) * nc + (col + dc);
                    if (!valid[j]) edge = true;
                }
            if (edge) weight[i] = 1.0;
        }

    const std::vector<double> reached = graph.accumulate(weight);
    std::vector<std::uint8_t> out(n, 0);
    for (std::size_t i = 0; i < n; ++i)
        out[i] = (valid[i] && reached[i] > 0.0) ? 1 : 0;
    return out;
}

// -----------------------------------------------------------------------------
// Areas by class
// -----------------------------------------------------------------------------
std::vector<sinmap_class_area_row> sinmap_class_areas(
    const sinmap_rasters& result,
    const std::vector<std::int32_t>& region)
{
    const std::vector<float>& cls = result.stability_class.get_data();
    const std::vector<std::uint8_t>& m = result.stability_class.get_mask();
    const double cell_area = result.stability_class.get_meta().cell_size
                           * result.stability_class.get_meta().cell_size;

    std::map<std::int32_t, std::map<int, std::int64_t>> counts;
    std::map<std::int32_t, std::int64_t> totals;
    for (std::size_t i = 0; i < cls.size(); ++i)
    {
        if (!m[i]) continue;
        const std::int32_t id = region.empty() ? 0 : region[i];
        ++counts[id][static_cast<int>(cls[i])];
        ++totals[id];
    }
    std::vector<sinmap_class_area_row> rows;
    for (const auto& by_id : counts)
        for (const auto& by_class : by_id.second)
        {
            sinmap_class_area_row r;
            r.region_id = by_id.first;
            r.stability_class = by_class.first;
            r.n_cells = by_class.second;
            r.area_m2 = static_cast<double>(by_class.second) * cell_area;
            r.fraction_of_region = static_cast<double>(by_class.second)
                                 / static_cast<double>(totals[by_id.first]);
            rows.push_back(r);
        }
    return rows;
}

// -----------------------------------------------------------------------------
// SA-plot curves
// -----------------------------------------------------------------------------
double sinmap_saturated_slope_for_safety_factor(
    double target, double cohesion, double friction_coefficient,
    double water_to_soil_density_ratio)
{
    // target sin - C = cos k, k = (1 - r) t. With u = cos:
    //   (target^2 + k^2) u^2 + 2 C k u + C^2 - target^2 = 0,
    //   u = [-C k + target sqrt(target^2 + k^2 - C^2)] / (target^2 + k^2).
    const double k = (1.0 - water_to_soil_density_ratio) * friction_coefficient;
    const double C = cohesion;
    if (!(target > 0.0)) return k_inf;
    if (C >= target) return k_inf;                 // FS > target at every slope below vertical
    const double denom = target * target + k * k;
    const double disc = denom - C * C;
    if (disc < 0.0) return k_inf;
    const double u = (-C * k + target * std::sqrt(disc)) / denom;
    if (!(u > 0.0) || !(u < 1.0)) return k_inf;
    const double theta = std::acos(u);
    // Squaring admits target sin - C = -k cos as well; keep only the true root.
    const double residual = target * std::sin(theta) - C - k * u;
    if (std::fabs(residual) > 1e-9 * std::max(1.0, target)) return k_inf;
    return std::tan(theta);
}

std::vector<sinmap_curve_point> sinmap_slope_area_curves(
    const sinmap_parameters& p,
    double slope_max,
    double area_min,
    double area_max,
    int n_slopes)
{
    validate_sinmap_parameters(p);
    if (!(slope_max > 0.0) || !(area_min > 0.0) || !(area_max > area_min) || n_slopes < 2)
        throw std::invalid_argument("sinmap_slope_area_curves: need slope_max > 0, 0 < area_min < area_max, n_slopes >= 2");

    const double r = p.water_to_soil_density_ratio;
    const double x1 = 1.0 / p.transmissivity_over_recharge_max_metres;
    const double x2 = 1.0 / p.transmissivity_over_recharge_min_metres;

    std::vector<double> slopes(static_cast<std::size_t>(n_slopes));
    for (int i = 0; i < n_slopes; ++i)
        slopes[static_cast<std::size_t>(i)] = slope_max * (i + 1) / n_slopes;

    std::vector<sinmap_curve_point> pts;
    auto add = [&](const std::string& curve, double S, double a)
    {
        if (a >= area_min && a <= area_max && S > 0.0 && S <= slope_max)
            pts.push_back({curve, S, a});
    };

    // Saturation lines (A-10).
    for (double S : slopes)
    {
        const double s = std::sin(std::atan(S));
        add("always_saturated_above", S, s / x1);
        add("never_saturated_below", S, s / x2);
    }

    // FS curves.
    struct fs_curve { const char* name; double target, C, t, x; };
    const fs_curve curves[] = {
        {"fs_min_1.5", 1.5, p.dimensionless_cohesion_min, p.friction_coefficient_min, x2},
        {"fs_min_1.0", 1.0, p.dimensionless_cohesion_min, p.friction_coefficient_min, x2},
        {"fs_max_1.0", 1.0, p.dimensionless_cohesion_max, p.friction_coefficient_max, x1},
    };
    for (const fs_curve& fc : curves)
    {
        // Wetness at which FS(theta, w) = target:  w = (1 - (target s - C)/(t c)) / r.
        auto wetness_needed = [&](double S) -> double
        {
            const double th = std::atan(S);
            const double s = std::sin(th), c = std::cos(th);
            if (fc.t <= 0.0 || r <= 0.0) return std::numeric_limits<double>::quiet_NaN();
            return (1.0 - (fc.target * s - fc.C) / (fc.t * c)) / r;
        };
        const double S_sat = sinmap_saturated_slope_for_safety_factor(fc.target, fc.C, fc.t, r);
        if (std::isfinite(S_sat) && S_sat > 0.0 && S_sat <= slope_max)
        {
            // Vertical part: from the top of the plot down to the saturation line.
            const double s = std::sin(std::atan(S_sat));
            add(fc.name, S_sat, area_max);
            add(fc.name, S_sat, std::min(area_max, s / fc.x));
        }
        for (double S : slopes)
        {
            if (std::isfinite(S_sat) && S <= S_sat) continue;
            const double w = wetness_needed(S);
            if (!(w > 0.0) || w >= 1.0) continue;
            const double s = std::sin(std::atan(S));
            add(fc.name, S, w * s / fc.x);
        }
    }
    return pts;
}

} // namespace lsdtt
