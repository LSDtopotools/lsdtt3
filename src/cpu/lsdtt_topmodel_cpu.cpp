// lsdtt_topmodel_cpu.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html

#include "lsdtt_topmodel_cpu.hpp"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <stdexcept>
#include <string>

namespace lsdtt
{

namespace
{

constexpr double k_ln_100 = 4.605170185988091;   // cm -> m on the depth integral
constexpr double k_ln_24  = 3.1780538303479458;  // per day -> per hour

void require_positive(double v, const char* name, const char* who)
{
    if (!(std::isfinite(v) && v > 0.0))
        throw std::invalid_argument(
            std::string(who) + ": " + name + " must be finite and > 0 (got " +
            std::to_string(v) + ")");
}

} // namespace

double sti_ln_kd_to_ln_t0_m2_per_day(double sti_ln_kd) noexcept
{
    return sti_ln_kd - k_ln_100;
}

double sti_ln_kd_to_ln_t0_m2_per_hour(double sti_ln_kd) noexcept
{
    return sti_ln_kd - k_ln_100 - k_ln_24;
}

topmodel_column_fluxes topmodel_column_step(
    double rainfall, double potential_et, double local_deficit,
    const topmodel_parameters& params,
    double& s_rz, double& s_uz, double& truncated)
{
    const double p  = std::max(0.0, rainfall);
    const double ep = std::max(0.0, potential_et);
    const double d_local = std::max(0.0, local_deficit);

    topmodel_column_fluxes out;

    // Root zone: rainfall fills the deficit, any surplus passes down.
    s_rz -= p;
    if (s_rz < 0.0) { s_uz -= s_rz; s_rz = 0.0; }

    // The unsaturated store cannot exceed the local deficit; the surplus is
    // saturation excess.
    if (s_uz > d_local)
    {
        out.saturation_excess = s_uz - d_local;
        s_uz = d_local;
    }

    // Gravity drainage to the saturated zone.
    if (d_local > 0.0)
    {
        double uz = s_uz / (d_local * params.td_hours_per_metre * params.dt_hours);
        if (uz > s_uz) uz = s_uz;
        s_uz -= uz;
        // 95.02's truncation deletes water rather than routing it, so the
        // deleted volume is handed back to the caller to account for.
        if (s_uz < 1.0e-7) { truncated += s_uz; s_uz = 0.0; }
        out.drainage = uz;
    }

    // Evapotranspiration from the root zone, linear in how wet it is.
    if (ep > 0.0)
    {
        double ea = ep * (1.0 - s_rz / params.srmax_metres);
        if (ea > params.srmax_metres - s_rz) ea = params.srmax_metres - s_rz;
        if (ea < 0.0) ea = 0.0;
        s_rz += ea;
        out.actual_et = ea;
    }

    return out;
}

lsdtt_topmodel_index_cpu::lsdtt_topmodel_index_cpu(
    std::vector<topmodel_unit> units,
    const topmodel_parameters& params,
    topmodel_excess_rule       rule)
    : units_(std::move(units)), params_(params), rule_(rule)
{
    const char* who = "lsdtt_topmodel_index_cpu";

    if (units_.empty())
        throw std::invalid_argument(std::string(who) + ": unit list is empty");

    require_positive(params_.szm_metres, "szm_metres", who);
    require_positive(params_.srmax_metres, "srmax_metres", who);
    require_positive(params_.td_hours_per_metre, "td_hours_per_metre", who);
    require_positive(params_.dt_hours, "dt_hours", who);
    require_positive(params_.q0_metres_per_step, "q0_metres_per_step", who);

    if (!(std::isfinite(params_.sr0_metres) && params_.sr0_metres >= 0.0))
        throw std::invalid_argument(
            std::string(who) + ": sr0_metres must be finite and >= 0");

    double area_sum = 0.0;
    for (const topmodel_unit& u : units_)
    {
        if (!std::isfinite(u.topographic_index))
            throw std::invalid_argument(
                std::string(who) + ": a unit has a non-finite topographic index");
        if (!(u.area_fraction >= 0.0) || !std::isfinite(u.area_fraction))
            throw std::invalid_argument(
                std::string(who) + ": a unit has a negative or non-finite area fraction");
        area_sum += u.area_fraction;
    }
    if (std::fabs(area_sum - 1.0) > 1.0e-6)
        throw std::invalid_argument(
            std::string(who) + ": area fractions sum to " + std::to_string(area_sum) +
            ", not 1");

    // Areal mean of the topographic index. For the trapezoid (histogram)
    // representation this must be the trapezoidal integral between adjacent bin
    // edges, exactly as 95.02's TREAD forms TL; for per-unit (raster) cells it
    // is the ordinary area-weighted mean. Using the wrong one here shifts every
    // local deficit by a constant and silently rescales the whole model.
    if (std::isfinite(params_.lambda_override))
    {
        lambda_ = params_.lambda_override;
    }
    else if (rule_ == topmodel_excess_rule::trapezoid && units_.size() >= 2)
    {
        lambda_ = 0.0;
        for (std::size_t i = 1; i < units_.size(); ++i)
            lambda_ += units_[i].area_fraction *
                       (units_[i].topographic_index + units_[i - 1].topographic_index) * 0.5;
    }
    else
    {
        lambda_ = 0.0;
        for (const topmodel_unit& u : units_)
            lambda_ += u.area_fraction * u.topographic_index;
    }

    // szq = exp(ln(T0) + ln(dt) - lambda), i.e. the discharge at zero deficit.
    szq_ = std::exp(params_.ln_t0 + std::log(params_.dt_hours) - lambda_);
    if (!(szq_ > 0.0) || !std::isfinite(szq_))
        throw std::invalid_argument(
            std::string(who) + ": szq is not finite and positive; check ln_t0 "
            "and the topographic index values");

    s_rz_.assign(units_.size(), params_.sr0_metres);
    s_uz_.assign(units_.size(), 0.0);

    // S_bar from the initial discharge: q0 = szq * exp(-sbar/m).
    sbar_ = -params_.szm_metres * std::log(params_.q0_metres_per_step / szq_);

    balance_.initial_storage = total_storage();
}

double lsdtt_topmodel_index_cpu::total_storage() const noexcept
{
    // Deficits are negative storage; S_uz is positive storage.
    double s = -sbar_;
    for (std::size_t i = 0; i < units_.size(); ++i)
        s += units_[i].area_fraction * (s_uz_[i] - s_rz_[i]);
    return s;
}

std::vector<double> lsdtt_topmodel_index_cpu::local_deficits() const
{
    std::vector<double> d(units_.size());
    for (std::size_t i = 0; i < units_.size(); ++i)
        d[i] = std::max(0.0,
            sbar_ + params_.szm_metres * (lambda_ - units_[i].topographic_index));
    return d;
}

topmodel_step_result lsdtt_topmodel_index_cpu::step(double rainfall_metres,
                                                    double potential_et_metres)
{
    return step_impl(&rainfall_metres, &potential_et_metres, 0);
}

topmodel_step_result lsdtt_topmodel_index_cpu::step(
    const std::vector<double>& rainfall_per_unit,
    const std::vector<double>& potential_et_per_unit)
{
    const char* who = "lsdtt_topmodel_index_cpu::step";
    if (rainfall_per_unit.size() != units_.size() ||
        potential_et_per_unit.size() != units_.size())
        throw std::invalid_argument(
            std::string(who) + ": per-unit forcing must have " +
            std::to_string(units_.size()) + " entries (got " +
            std::to_string(rainfall_per_unit.size()) + " rainfall and " +
            std::to_string(potential_et_per_unit.size()) + " PET)");
    return step_impl(rainfall_per_unit.data(), potential_et_per_unit.data(), 1);
}

topmodel_step_result lsdtt_topmodel_index_cpu::step_impl(
    const double* rainfall, const double* potential_et, std::ptrdiff_t stride)
{
    const std::size_t n = units_.size();

    topmodel_step_result out;

    // With per-unit forcing the reported rainfall and PET are area-weighted
    // means, so the water balance still closes on domain totals.
    double p_mean = 0.0, ep_mean = 0.0;
    for (std::size_t i = 0; i < n; ++i)
    {
        const std::ptrdiff_t k = static_cast<std::ptrdiff_t>(i) * stride;
        p_mean  += units_[i].area_fraction * std::max(0.0, rainfall[k]);
        ep_mean += units_[i].area_fraction * std::max(0.0, potential_et[k]);
    }
    out.rainfall     = p_mean;
    out.potential_et = ep_mean;

    double quz       = 0.0;
    double qof       = 0.0;
    double sae       = 0.0;
    double acm       = 0.0;
    double truncated = 0.0;

    std::vector<double> ex(n, 0.0);

    for (std::size_t i = 0; i < n; ++i)
    {
        const double acf = units_[i].area_fraction;

        const std::ptrdiff_t k = static_cast<std::ptrdiff_t>(i) * stride;
        const double p  = std::max(0.0, rainfall[k]);
        const double ep = std::max(0.0, potential_et[k]);

        // Local deficit from the CURRENT S_bar -- start-of-step, as 95.02 does.
        const double d_local = std::max(0.0,
            sbar_ + params_.szm_metres * (lambda_ - units_[i].topographic_index));

        double trunc_here = 0.0;
        const topmodel_column_fluxes cf = topmodel_column_step(
            p, ep, d_local, params_, s_rz_[i], s_uz_[i], trunc_here);

        truncated += trunc_here * acf;
        ex[i] = cf.saturation_excess;
        if (d_local > 0.0) quz += cf.drainage * acf;
        sae += cf.actual_et * acf;

        // --- saturation excess ---
        if (rule_ == topmodel_excess_rule::per_unit)
        {
            // A cell is saturated or it is not.
            if (ex[i] > 0.0)
            {
                qof += ex[i] * acf;
                acm += acf;
            }
        }
        else if (i > 0)
        {
            // 95.02's histogram rule. A bin spans a RANGE of index values, so a
            // partly saturated bin contributes a linearly interpolated wedge
            // rather than all-or-nothing. See the header on why this does not
            // equal the per-unit rule at a finite bin count.
            double acf_bin = acf;
            double of      = 0.0;
            if (ex[i] > 0.0)
            {
                of = units_[i].area_fraction * (ex[i - 1] + ex[i]) * 0.5;
                acm += acf_bin;
            }
            else if (ex[i - 1] > 0.0)
            {
                acf_bin *= ex[i - 1] / (ex[i - 1] - ex[i]);
                of = acf_bin * ex[i - 1] * 0.5;
                acm += acf_bin;
            }
            qof += of;
        }
    }

    // Baseflow from the exponential store, then update S_bar. Both use
    // start-of-step S_bar, which is 95.02's explicit scheme.
    const double qb = szq_ * std::exp(-sbar_ / params_.szm_metres);
    sbar_ += (-quz + qb);

    out.actual_et            = sae;
    out.saturation_excess    = qof;
    out.unsaturated_drainage = quz;
    out.baseflow             = qb;
    out.total_runoff         = qb + qof;
    out.mean_deficit         = sbar_;
    out.contributing_area    = acm;

    balance_.total_rainfall   += p_mean;
    balance_.total_actual_et  += sae;
    balance_.total_runoff     += out.total_runoff;
    balance_.truncation_loss  += truncated;
    balance_.final_storage     = total_storage();

    return out;
}

std::vector<topmodel_unit> make_units_from_histogram(
    const std::vector<double>& area_fractions,
    const std::vector<double>& topographic_index_values)
{
    const char* who = "make_units_from_histogram";

    if (area_fractions.size() != topographic_index_values.size())
        throw std::invalid_argument(
            std::string(who) + ": area_fractions has " +
            std::to_string(area_fractions.size()) + " entries but there are " +
            std::to_string(topographic_index_values.size()) + " index values");

    if (area_fractions.size() < 2)
        throw std::invalid_argument(
            std::string(who) + ": need at least two bins to form a histogram");

    double total = 0.0;
    for (double a : area_fractions)
    {
        if (!(a >= 0.0) || !std::isfinite(a))
            throw std::invalid_argument(
                std::string(who) + ": a bin has a negative or non-finite area");
        total += a;
    }
    if (!(total > 0.0))
        throw std::invalid_argument(std::string(who) + ": bin areas sum to zero");

    std::vector<topmodel_unit> units(area_fractions.size());
    for (std::size_t i = 0; i < area_fractions.size(); ++i)
    {
        units[i].area_fraction     = area_fractions[i] / total;
        units[i].topographic_index = topographic_index_values[i];
    }
    return units;
}

double topmodel_lambda_from_histogram(
    const std::vector<double>& area_fractions,
    const std::vector<double>& topographic_index_values)
{
    const char* who = "topmodel_lambda_from_histogram";
    if (area_fractions.size() != topographic_index_values.size())
        throw std::invalid_argument(
            std::string(who) + ": areas and index values differ in length");
    if (area_fractions.size() < 2)
        throw std::invalid_argument(std::string(who) + ": need at least two bins");

    double total = 0.0;
    for (double a : area_fractions) total += a;
    if (!(total > 0.0))
        throw std::invalid_argument(std::string(who) + ": bin areas sum to zero");

    double lambda = 0.0;
    for (std::size_t j = 1; j < area_fractions.size(); ++j)
        lambda += (area_fractions[j] / total) *
                  (topographic_index_values[j] + topographic_index_values[j - 1]) * 0.5;
    return lambda;
}

std::vector<topmodel_unit> make_units_from_twi(
    const std::vector<float>&        twi_values,
    const std::vector<std::uint8_t>& twi_mask)
{
    const char* who = "make_units_from_twi";

    if (twi_values.size() != twi_mask.size())
        throw std::invalid_argument(
            std::string(who) + ": values and mask differ in length");

    std::size_t n_valid = 0;
    for (std::size_t i = 0; i < twi_mask.size(); ++i)
        if (twi_mask[i] != 0 && std::isfinite(twi_values[i])) ++n_valid;

    if (n_valid == 0)
        throw std::invalid_argument(
            std::string(who) + ": the TWI raster has no valid cells");

    const double w = 1.0 / static_cast<double>(n_valid);

    std::vector<topmodel_unit> units;
    units.reserve(n_valid);
    for (std::size_t i = 0; i < twi_mask.size(); ++i)
    {
        if (twi_mask[i] == 0 || !std::isfinite(twi_values[i])) continue;
        topmodel_unit u;
        u.topographic_index = static_cast<double>(twi_values[i]);
        u.area_fraction     = w;
        units.push_back(u);
    }
    return units;
}

} // namespace lsdtt
