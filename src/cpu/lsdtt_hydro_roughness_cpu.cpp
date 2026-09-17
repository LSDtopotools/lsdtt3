// lsdtt_hydro_roughness_cpu.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html

#include "lsdtt_hydro_roughness_cpu.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>
#include <utility>

#include "../utils/lsdtt_array2d.hpp"

namespace lsdtt
{

namespace
{

constexpr double k_gravity = 9.80665;

void require_positive(double v, const char* name, const char* who)
{
    if (!(std::isfinite(v) && v > 0.0))
        throw std::invalid_argument(
            std::string(who) + ": " + name + " must be finite and > 0 (got " +
            std::to_string(v) + ")");
}

void require_same_grid(const lsdtt_raster_cpu& r, const lsdtt_raster_meta& meta,
                       const char* name, const char* who)
{
    const lsdtt_raster_meta& m = r.get_meta();
    if (m.n_rows != meta.n_rows || m.n_cols != meta.n_cols)
        throw std::invalid_argument(
            std::string(who) + ": " + name + " is " + std::to_string(m.n_rows) +
            "x" + std::to_string(m.n_cols) + " but the grid is " +
            std::to_string(meta.n_rows) + "x" + std::to_string(meta.n_cols));
}

void require_input(const lsdtt_raster_cpu* r, const char* name,
                   const char* mode, const char* who)
{
    if (r == nullptr)
        throw std::invalid_argument(
            std::string(who) + ": " + name + " is required in " + mode + " mode");
}

bool cell_ok(const lsdtt_raster_cpu* r, int row, int col)
{
    return r != nullptr && r->get_mask_array()(row, col) != 0;
}

double cell_value(const lsdtt_raster_cpu* r, int row, int col)
{
    return static_cast<double>(r->get_data_array()(row, col));
}

} // namespace

double strickler_roughness(double d50_metres, double strickler_coefficient) noexcept
{
    if (!(d50_metres > 0.0)) return 0.0;
    return strickler_coefficient * std::pow(d50_metres, 1.0 / 6.0);
}

double jarrett_roughness(double slope, double hydraulic_radius_metres) noexcept
{
    if (!(slope > 0.0) || !(hydraulic_radius_metres > 0.0)) return 0.0;
    return 0.39 * std::pow(slope, 0.38) * std::pow(hydraulic_radius_metres, -0.16);
}

double ferguson_sqrt_8_over_f(double hydraulic_radius_metres, double d84_metres,
                              double a1, double a2) noexcept
{
    if (!(hydraulic_radius_metres > 0.0) || !(d84_metres > 0.0)) return 0.0;
    const double hd  = hydraulic_radius_metres / d84_metres;
    const double num = a1 * a2 * hd;
    const double den = std::sqrt(a1 * a1 + a2 * a2 * std::pow(hd, 5.0 / 3.0));
    if (!(den > 0.0)) return 0.0;
    return num / den;
}

double ferguson_vpe_velocity(double hydraulic_radius_metres, double d84_metres,
                             double slope, double a1, double a2) noexcept
{
    if (!(slope > 0.0)) return 0.0;
    return ferguson_sqrt_8_over_f(hydraulic_radius_metres, d84_metres, a1, a2) *
           std::sqrt(k_gravity * hydraulic_radius_metres * slope);
}

double ferguson_equivalent_manning_n(double hydraulic_radius_metres,
                                     double d84_metres, double slope,
                                     double a1, double a2) noexcept
{
    const double v = ferguson_vpe_velocity(hydraulic_radius_metres, d84_metres,
                                           slope, a1, a2);
    if (!(v > 0.0)) return 0.0;
    return std::pow(hydraulic_radius_metres, 2.0 / 3.0) * std::sqrt(slope) / v;
}

namespace
{

// Median over channel cells, used when reference_ksn is left at zero so the
// scaled field is centred on n_ref rather than on an arbitrary absolute value.
double median_of(std::vector<double> v)
{
    if (v.empty()) return 0.0;
    const std::size_t mid = v.size() / 2;
    std::nth_element(v.begin(), v.begin() + static_cast<std::ptrdiff_t>(mid), v.end());
    return v[mid];
}

double ksn_at(const roughness_field_inputs& in,
              const channel_roughness_options& opts,
              int row, int col)
{
    if (opts.mode == channel_roughness_mode::ksn)
        return cell_value(in.ksn, row, col);

    const double s = std::max(cell_value(in.slope, row, col), opts.min_slope);
    const double a = cell_value(in.drainage_area, row, col);
    return s * std::pow(a, opts.chi_m_over_n);
}

} // namespace

lsdtt_raster_cpu compute_channel_roughness(
    const roughness_field_inputs&    inputs,
    const lsdtt_raster_meta&         meta,
    const channel_roughness_options& opts)
{
    const char* who = "compute_channel_roughness";

    require_positive(opts.channel_roughness, "channel_roughness", who);
    require_positive(opts.n_min, "n_min", who);
    require_positive(opts.n_max, "n_max", who);
    require_positive(opts.min_slope, "min_slope", who);
    if (opts.n_min > opts.n_max)
        throw std::invalid_argument(
            std::string(who) + ": n_min exceeds n_max");

    switch (opts.mode)
    {
        case channel_roughness_mode::constant:
            break;
        case channel_roughness_mode::drainage_area:
            require_input(inputs.drainage_area, "drainage_area", "drainage_area", who);
            require_positive(opts.reference_area_m2, "reference_area_m2", who);
            break;
        case channel_roughness_mode::ksn:
            require_input(inputs.ksn, "ksn", "ksn", who);
            break;
        case channel_roughness_mode::slope_area:
            require_input(inputs.slope, "slope", "slope_area", who);
            require_input(inputs.drainage_area, "drainage_area", "slope_area", who);
            break;
        case channel_roughness_mode::grain_size:
            require_input(inputs.d50, "d50", "grain_size", who);
            require_positive(opts.strickler_coefficient, "strickler_coefficient", who);
            break;
        case channel_roughness_mode::jarrett:
            require_input(inputs.slope, "slope", "jarrett", who);
            require_input(inputs.hydraulic_radius, "hydraulic_radius", "jarrett", who);
            break;
        case channel_roughness_mode::ferguson_vpe:
            require_input(inputs.slope, "slope", "ferguson_vpe", who);
            require_input(inputs.hydraulic_radius, "hydraulic_radius", "ferguson_vpe", who);
            require_input(inputs.d84, "d84", "ferguson_vpe", who);
            break;
    }

    const lsdtt_raster_cpu* grids[] = {
        inputs.drainage_area, inputs.slope, inputs.hydraulic_radius,
        inputs.d50, inputs.d84, inputs.ksn, inputs.landcover_roughness,
        inputs.channel_mask, inputs.surface_class
    };
    const char* names[] = {
        "drainage_area", "slope", "hydraulic_radius", "d50", "d84", "ksn",
        "landcover_roughness", "channel_mask", "surface_class"
    };
    for (std::size_t i = 0; i < sizeof(grids) / sizeof(grids[0]); ++i)
        if (grids[i]) require_same_grid(*grids[i], meta, names[i], who);

    const int   n_rows = meta.n_rows;
    const int   n_cols = meta.n_cols;
    const float ndv    = static_cast<float>(meta.no_data_value);

    auto is_channel = [&](int row, int col) {
        if (inputs.channel_mask == nullptr) return true;
        if (inputs.channel_mask->get_mask_array()(row, col) == 0) return false;
        return cell_value(inputs.channel_mask, row, col) != 0.0;
    };

    // Reference ksn: the median over channel cells unless one was given.
    double ref_ksn = opts.reference_ksn;
    const bool needs_ksn = (opts.mode == channel_roughness_mode::ksn ||
                            opts.mode == channel_roughness_mode::slope_area);
    if (needs_ksn && !(ref_ksn > 0.0))
    {
        std::vector<double> vals;
        for (int row = 0; row < n_rows; ++row)
            for (int col = 0; col < n_cols; ++col)
            {
                if (!is_channel(row, col)) continue;
                if (opts.mode == channel_roughness_mode::ksn)
                {
                    if (!cell_ok(inputs.ksn, row, col)) continue;
                }
                else
                {
                    if (!cell_ok(inputs.slope, row, col)) continue;
                    if (!cell_ok(inputs.drainage_area, row, col)) continue;
                }
                const double k = ksn_at(inputs, opts, row, col);
                if (k > 0.0 && std::isfinite(k)) vals.push_back(k);
            }
        ref_ksn = median_of(std::move(vals));
        if (!(ref_ksn > 0.0))
            throw std::invalid_argument(
                std::string(who) + ": no channel cell has a positive ksn, so no "
                "reference could be formed; supply reference_ksn explicitly");
    }

    lsdtt_array2d_f  n_data(n_rows, n_cols, ndv);
    lsdtt_array2d_u8 n_mask(n_rows, n_cols, std::uint8_t(0));

    for (int row = 0; row < n_rows; ++row)
        for (int col = 0; col < n_cols; ++col)
        {
            if (!is_channel(row, col)) continue;

            double n = 0.0;
            bool   ok = true;

            switch (opts.mode)
            {
                case channel_roughness_mode::constant:
                    n = opts.channel_roughness;
                    break;

                case channel_roughness_mode::drainage_area:
                {
                    if (!cell_ok(inputs.drainage_area, row, col)) { ok = false; break; }
                    const double a = cell_value(inputs.drainage_area, row, col);
                    if (!(a > 0.0)) { ok = false; break; }
                    n = opts.channel_roughness *
                        std::pow(a / opts.reference_area_m2, -opts.area_exponent);
                    break;
                }

                case channel_roughness_mode::ksn:
                case channel_roughness_mode::slope_area:
                {
                    if (opts.mode == channel_roughness_mode::ksn)
                    {
                        if (!cell_ok(inputs.ksn, row, col)) { ok = false; break; }
                    }
                    else
                    {
                        if (!cell_ok(inputs.slope, row, col) ||
                            !cell_ok(inputs.drainage_area, row, col)) { ok = false; break; }
                    }
                    const double k = ksn_at(inputs, opts, row, col);
                    if (!(k > 0.0)) { ok = false; break; }
                    n = opts.channel_roughness *
                        std::pow(k / ref_ksn, opts.ksn_exponent);
                    break;
                }

                case channel_roughness_mode::grain_size:
                {
                    if (!cell_ok(inputs.d50, row, col)) { ok = false; break; }
                    n = strickler_roughness(cell_value(inputs.d50, row, col),
                                            opts.strickler_coefficient);
                    if (!(n > 0.0)) ok = false;
                    break;
                }

                case channel_roughness_mode::jarrett:
                {
                    if (!cell_ok(inputs.slope, row, col) ||
                        !cell_ok(inputs.hydraulic_radius, row, col)) { ok = false; break; }
                    const double s = std::max(cell_value(inputs.slope, row, col),
                                              opts.min_slope);
                    const double r = cell_value(inputs.hydraulic_radius, row, col);
                    n = jarrett_roughness(s, r);
                    if (!(n > 0.0)) ok = false;
                    break;
                }

                case channel_roughness_mode::ferguson_vpe:
                {
                    if (!cell_ok(inputs.slope, row, col) ||
                        !cell_ok(inputs.hydraulic_radius, row, col) ||
                        !cell_ok(inputs.d84, row, col)) { ok = false; break; }
                    const double s = std::max(cell_value(inputs.slope, row, col),
                                              opts.min_slope);
                    n = ferguson_equivalent_manning_n(
                        cell_value(inputs.hydraulic_radius, row, col),
                        cell_value(inputs.d84, row, col),
                        s, opts.ferguson_a1, opts.ferguson_a2);
                    if (!(n > 0.0)) ok = false;
                    break;
                }
            }

            if (!ok || !std::isfinite(n)) continue;

            n = std::min(std::max(n, opts.n_min), opts.n_max);
            n_data(row, col) = static_cast<float>(n);
            n_mask(row, col) = std::uint8_t(1);
        }

    return lsdtt_raster_cpu(std::move(n_data), std::move(n_mask), meta);
}

roughness_field_result compute_roughness_field(
    const roughness_field_inputs&      inputs,
    const lsdtt_raster_meta&           meta,
    const roughness_partition_options& opts)
{
    const char* who = "compute_roughness_field";
    require_positive(opts.landcover_roughness_default,
                     "landcover_roughness_default", who);

    const lsdtt_raster_cpu channel_n =
        compute_channel_roughness(inputs, meta, opts.channel);

    const int   n_rows = meta.n_rows;
    const int   n_cols = meta.n_cols;
    const float ndv    = static_cast<float>(meta.no_data_value);

    lsdtt_array2d_f  n_data(n_rows, n_cols, ndv);
    lsdtt_array2d_u8 n_mask(n_rows, n_cols, std::uint8_t(0));

    roughness_partition_stats stats;

    const bool jarrett = (opts.channel.mode == channel_roughness_mode::jarrett);

    for (int row = 0; row < n_rows; ++row)
        for (int col = 0; col < n_cols; ++col)
        {
            // Channel wins: its roughness comes from the bed, not the bank.
            if (channel_n.get_mask_array()(row, col) != 0)
            {
                n_data(row, col) = channel_n.get_data_array()(row, col);
                n_mask(row, col) = std::uint8_t(1);
                ++stats.n_channel;

                if (jarrett && cell_ok(inputs.slope, row, col) &&
                    cell_ok(inputs.hydraulic_radius, row, col))
                {
                    const double s = cell_value(inputs.slope, row, col);
                    const double r = cell_value(inputs.hydraulic_radius, row, col);
                    if (s < opts.channel.jarrett_min_slope ||
                        s > opts.channel.jarrett_max_slope ||
                        r < opts.channel.jarrett_min_radius ||
                        r > opts.channel.jarrett_max_radius)
                        ++stats.n_outside_jarrett_range;
                }
                continue;
            }

            // Everything else takes land cover. Floodplain and hillslope draw
            // from the SAME source today -- the split is counted, not applied.
            // See the header on why it is carried anyway.
            const bool is_floodplain =
                cell_ok(inputs.surface_class, row, col) &&
                cell_value(inputs.surface_class, row, col) == opts.floodplain_class;

            double n = opts.landcover_roughness_default;
            bool   have_any = (inputs.landcover_roughness == nullptr);

            if (cell_ok(inputs.landcover_roughness, row, col))
            {
                const double v = cell_value(inputs.landcover_roughness, row, col);
                if (v > 0.0) { n = v; have_any = true; }
            }

            // A cell with no land cover AND no default source is unclassified:
            // left nodata rather than given a number nobody chose.
            if (!have_any && inputs.landcover_roughness != nullptr)
            {
                ++stats.n_unclassified;
                continue;
            }

            n = std::min(std::max(n, opts.channel.n_min), opts.channel.n_max);
            n_data(row, col) = static_cast<float>(n);
            n_mask(row, col) = std::uint8_t(1);

            if (is_floodplain) ++stats.n_floodplain; else ++stats.n_hillslope;
        }

    roughness_field_result out;
    out.roughness = lsdtt_raster_cpu(std::move(n_data), std::move(n_mask), meta);
    out.stats     = stats;
    return out;
}

} // namespace lsdtt
