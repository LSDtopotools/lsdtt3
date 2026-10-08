// lsdtt_infiltration_excess_cpu.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html

#include "lsdtt_infiltration_excess_cpu.hpp"
#include "lsdtt_distributed_hydrology_argument_checks_cpu.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>

namespace lsdtt
{

namespace
{

void require_non_negative(double v, const char* name, const char* who)
{
    if (!(std::isfinite(v) && v >= 0.0))
        throw std::invalid_argument(
            std::string(who) + ": " + name + " must be finite and >= 0 (got " +
            std::to_string(v) + ")");
}

} // namespace

lsdtt_infiltration_excess_cpu::lsdtt_infiltration_excess_cpu(
    const infiltration_excess_parameters& params)
    : params_(params)
{
    const char* who = "lsdtt_infiltration_excess_cpu";

    require_positive_finite(params_.xk0_metres_per_hour, "xk0_metres_per_hour", who);
    require_positive_finite(params_.szm_metres, "szm_metres", who);
    require_positive_finite(params_.tolerance, "tolerance", who);
    require_non_negative(params_.wetting_front_suction_metres,
                         "wetting_front_suction_metres", who);
    require_non_negative(params_.delta_theta, "delta_theta", who);

    if (params_.max_iterations < 1)
        throw std::invalid_argument(
            std::string(who) + ": max_iterations must be >= 1");
    if (params_.n_series_terms < 1)
        throw std::invalid_argument(
            std::string(who) + ": n_series_terms must be >= 1");

    // Everything that depends only on the parameters, computed ONCE. In the
    // original these four sit inside the per-call path, and exp(szf*cd) inside
    // the Newton loop itself.
    szf_        = 1.0 / params_.szm_metres;
    cd_         = params_.wetting_front_suction_metres * params_.delta_theta;
    exp_szf_cd_ = std::exp(szf_ * cd_);
    xkf_szf_    = params_.xk0_metres_per_hour * szf_;
}

double lsdtt_infiltration_excess_cpu::infiltration_capacity(double f) const noexcept
{
    // r(f) = -xk0 szf (cd + f) / (1 - exp(szf f))
    //
    // For f > 0 the denominator is negative and the whole expression positive.
    // At f = 0 it is 0/0; the limit is xk0 (1 + cd szf), which is what a dry
    // surface should give, so it is returned rather than a NaN.
    const double e = std::exp(szf_ * f);
    const double den = 1.0 - e;
    if (std::fabs(den) < 1.0e-300)
        return params_.xk0_metres_per_hour * (1.0 + cd_ * szf_);
    return -xkf_szf_ * (cd_ + f) / den;
}

double lsdtt_infiltration_excess_cpu::series_sum(
    double x, double* last_term_relative) const noexcept
{
    // SUM(j=1..n) x^j / (j * j!)  by the recurrence
    //     a_1 = x,  a_{j+1} = a_j * x * j / (j+1)^2
    //
    // The original evaluates pow(x, j) and a running factorial for every term.
    // This is algebraically the same series with no pow, no factorial, and no
    // overflow risk in j!.
    double term = x;
    double sum  = x;
    for (int j = 1; j < params_.n_series_terms; ++j)
    {
        const double jd = static_cast<double>(j);
        term *= x * jd / ((jd + 1.0) * (jd + 1.0));
        sum  += term;
    }
    if (last_term_relative)
        *last_term_relative = (std::fabs(sum) > 0.0)
            ? std::fabs(term) / std::fabs(sum) : 0.0;
    return sum;
}

double lsdtt_infiltration_excess_cpu::integral_term(
    double f, bool* poorly_converged) const noexcept
{
    const double fc = f + cd_;
    if (!(fc > 0.0))
    {
        if (poorly_converged) *poorly_converged = false;
        return 0.0;
    }
    double rel = 0.0;
    const double sum = series_sum(fc * szf_, &rel);
    if (poorly_converged) *poorly_converged = (rel > 1.0e-6);

    // log(fc) once, not twice as in the original.
    const double log_fc = std::log(fc);
    return log_fc - (log_fc + sum) / exp_szf_cd_;
}

infiltration_step_result lsdtt_infiltration_excess_cpu::step(
    infiltration_state& st,
    double rainfall_rate_metres_per_hour,
    double dt_hours,
    int    timestep_index) const
{
    const char* who = "lsdtt_infiltration_excess_cpu::step";
    require_positive_finite(dt_hours, "dt_hours", who);

    const double rint = std::max(0.0, rainfall_rate_metres_per_hour);
    const double dt   = dt_hours;
    const double t_end   = static_cast<double>(timestep_index) * dt;
    const double t_start = (static_cast<double>(timestep_index) - 1.0) * dt;

    infiltration_step_result out;

    // No rain, nothing to do. The original reaches the same place through its
    // caller rather than here.
    if (!(rint > 0.0))
    {
        st.ponded = false;
        out.ponded = false;
        return out;
    }

    // "No ponding this step": all the rain goes in. The original's label 20.
    auto no_ponding = [&]() {
        st.ponded = false;
        const double df = rint * dt;
        st.cumulative_infiltration += df;
        out.infiltration = df;
        out.excess       = 0.0;
        out.ponded       = false;
        return out;
    };

    bool go_to_newton = st.ponded;
    double f = st.f_estimate;

    if (!go_to_newton)
    {
        bool have_ponding_start = false;

        if (st.cumulative_infiltration > 0.0)
        {
            // Does the surface already pond at the start of the step?
            if (infiltration_capacity(st.cumulative_infiltration) < rint)
            {
                st.ponding_time_hours = t_start;
                f = st.cumulative_infiltration;
                have_ponding_start = true;
            }
        }

        if (!have_ponding_start)
        {
            // Would the step's rain exceed capacity by the end of it?
            const double f2_end = st.cumulative_infiltration + dt * rint;
            if (!(f2_end > 0.0)) return no_ponding();
            if (infiltration_capacity(f2_end) > rint) return no_ponding();

            // Bisect for the depth at which capacity equals the rainfall rate.
            // Kept as bisection, not Newton: the capacity function is monotone
            // but its derivative is awkward near f = 0, and twenty bisections
            // on an already-bracketed root are cheap next to the Newton solve
            // that follows.
            double f1 = st.cumulative_infiltration;
            double f2 = f2_end;
            f = st.cumulative_infiltration + infiltration_capacity(f2_end) * dt;

            for (int i = 0; i < params_.max_iterations; ++i)
            {
                const double r2 = infiltration_capacity(f);
                if (r2 > rint)
                {
                    f1 = f;
                    f  = (f2 + f) * 0.5;
                    if (std::fabs(f - f1) < params_.tolerance) break;
                }
                else
                {
                    f2 = f;
                    f  = (f1 + f) * 0.5;
                    if (std::fabs(f - f2) < params_.tolerance) break;
                }
            }

            st.ponding_time_hours =
                t_start + (f - st.cumulative_infiltration) / rint;
            if (st.ponding_time_hours > t_end) return no_ponding();
        }

        // Ponding begins in this step: fix the definite-integral constant at
        // the ponding depth, then advance the guess by half the remaining
        // rainfall, as the original does.
        bool poor = false;
        st.integral_constant = integral_term(f, &poor);
        out.series_poorly_converged = poor;

        st.ponded = true;
        f += 0.5 * rint * (t_end - st.ponding_time_hours);
    }

    // --- Newton-Raphson for F(t) ---
    int it = 0;
    double df_step = 0.0;
    for (; it < params_.max_iterations; ++it)
    {
        const double fc = f + cd_;
        if (!(fc > 0.0)) break;

        bool poor = false;
        const double integral = integral_term(f, &poor);
        if (poor) out.series_poorly_converged = true;

        const double func = -(integral - st.integral_constant) / xkf_szf_
                          - (t_end - st.ponding_time_hours);
        const double dfunc = (std::exp(szf_ * f) - 1.0) / (xkf_szf_ * fc);
        if (!(std::fabs(dfunc) > 0.0)) break;

        df_step = -func / dfunc;
        f += df_step;
        if (std::fabs(df_step) <= params_.tolerance) { ++it; break; }
    }
    out.iterations = it;
    out.converged  = (std::fabs(df_step) <= params_.tolerance);

    // Did the solve actually put less water in than the rain delivered? If not,
    // there is no excess and the step is an ordinary all-infiltrates step.
    //
    // The original compares F against CUMF + RINT, a depth against a depth plus
    // a RATE. `use_dimensional_ponding_test` selects CUMF + RINT*dt instead; the
    // default reproduces 95.02 exactly, including at dt != 1 where the two
    // genuinely differ.
    const double limit = params_.use_dimensional_ponding_test
        ? st.cumulative_infiltration + rint * dt
        : st.cumulative_infiltration + rint;

    if (f < limit)
    {
        const double df = f - st.cumulative_infiltration;
        st.cumulative_infiltration = f;
        st.f_estimate = f + df;          // the original's guess for the next step
        out.infiltration = std::max(0.0, df);
        out.excess       = std::max(0.0, rint * dt - out.infiltration);
        out.ponded       = true;
        return out;
    }

    return no_ponding();
}

} // namespace lsdtt
