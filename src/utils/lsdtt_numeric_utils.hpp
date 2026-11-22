// lsdtt_numeric_utils.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html

#pragma once

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

namespace lsdtt {

// Generates start, start + step, ... through the inclusive end, allowing the
// final accumulated value to lie up to half a step beyond end. This preserves
// the established sweep convention used by parameterised numeric analyses.
inline std::vector<double> inclusive_numeric_range(
    double start,
    double end,
    double step)
{
    if (!(std::isfinite(start) && std::isfinite(end) && std::isfinite(step)))
        throw std::invalid_argument("inclusive_numeric_range: values must be finite");
    if (step <= 0.0)
        throw std::invalid_argument("inclusive_numeric_range: step must be positive");
    if (end < start)
        throw std::invalid_argument("inclusive_numeric_range: end must not precede start");

    std::vector<double> values;
    for (double value = start; value <= end + 0.5 * step; value += step)
        values.push_back(value);
    return values;
}

// -----------------------------------------------------------------------------
// Bounded one-dimensional minimisation: Brent's method.
//
// Brent, R. P. (1973). Algorithms for Minimization without Derivatives,
// Prentice-Hall, ch. 5 ("localmin"); the same algorithm as Netlib's fmin.
// Golden-section steps, accelerated by successive parabolic interpolation
// when the parabola is trustworthy. No derivatives.
//
// CONTRACT
//   * Finds a LOCAL minimum of f strictly inside (lower, upper). The end points
//     themselves are never evaluated; a function that decreases all the way to
//     an end converges to within about x_tolerance of that end. A caller whose
//     function can have several minima must split the interval and compare
//     (for example the scaLr pore-pressure case, which has a kink at the water
//     table and is minimised on each side separately).
//   * The x accuracy is 2 * (sqrt(machine epsilon) * |x| + x_tolerance / 3), as
//     in Netlib fmin. Asking for better than sqrt(epsilon) relative accuracy
//     buys nothing: near a smooth minimum f is flat to that order.
//   * A NaN from f is treated as +infinity, so regions where f is undefined
//     (inadmissible geometry, say) are avoided rather than poisoning the
//     search. +infinity itself is allowed.
//   * lower == upper returns that point after one evaluation.
//
// Throws std::invalid_argument for non-finite or reversed bounds, a
// non-positive tolerance, or max_iterations < 1.
// -----------------------------------------------------------------------------
struct brent_minimum
{
    double x = 0.0;           // abscissa of the minimum found
    double f = 0.0;           // f(x)
    int iterations = 0;       // function evaluations after the first
    bool converged = false;   // false: stopped at max_iterations
};

template <class Function>
brent_minimum minimise_brent_bounded(
    Function&& f,
    double lower,
    double upper,
    double x_tolerance = 1.0e-10,
    int max_iterations = 200)
{
    if (!(std::isfinite(lower) && std::isfinite(upper)))
        throw std::invalid_argument("minimise_brent_bounded: bounds must be finite");
    if (lower > upper)
        throw std::invalid_argument("minimise_brent_bounded: lower exceeds upper");
    if (!(x_tolerance > 0.0))
        throw std::invalid_argument("minimise_brent_bounded: x_tolerance must be positive");
    if (max_iterations < 1)
        throw std::invalid_argument("minimise_brent_bounded: max_iterations must be >= 1");

    auto evaluate = [&f](double x) -> double
    {
        const double v = static_cast<double>(f(x));
        return std::isnan(v) ? std::numeric_limits<double>::infinity() : v;
    };

    brent_minimum result;
    if (lower == upper)
    {
        result.x = lower;
        result.f = evaluate(lower);
        result.converged = true;
        return result;
    }

    const double golden = 0.5 * (3.0 - std::sqrt(5.0));   // 0.381966...
    const double eps = std::sqrt(std::numeric_limits<double>::epsilon());

    double a = lower;
    double b = upper;
    double v = a + golden * (b - a);
    double w = v;
    double x = v;
    double e = 0.0;     // distance moved on the step before last
    double d = 0.0;     // distance moved on the last step
    double fx = evaluate(x);
    double fv = fx;
    double fw = fx;

    for (int iteration = 0; iteration < max_iterations; ++iteration)
    {
        const double xm   = 0.5 * (a + b);
        const double tol1 = eps * std::fabs(x) + x_tolerance / 3.0;
        const double tol2 = 2.0 * tol1;

        if (std::fabs(x - xm) <= tol2 - 0.5 * (b - a))
        {
            result.converged = true;
            break;
        }

        bool golden_step = true;
        if (std::fabs(e) > tol1)
        {
            // Trial parabola through (v, fv), (w, fw), (x, fx).
            double r = (x - w) * (fx - fv);
            double q = (x - v) * (fx - fw);
            double p = (x - v) * q - (x - w) * r;
            q = 2.0 * (q - r);
            if (q > 0.0) p = -p;
            q = std::fabs(q);
            r = e;
            e = d;

            // Accept only if it falls inside (a, b) and moves less than half the
            // step before last (so the bracket keeps shrinking).
            if (std::fabs(p) < std::fabs(0.5 * q * r)
                && p > q * (a - x) && p < q * (b - x))
            {
                d = p / q;
                const double u = x + d;
                if ((u - a) < tol2 || (b - u) < tol2)
                    d = (x < xm) ? tol1 : -tol1;
                golden_step = false;
            }
        }
        if (golden_step)
        {
            e = (x < xm) ? (b - x) : (a - x);
            d = golden * e;
        }

        const double u  = (std::fabs(d) >= tol1) ? x + d : x + ((d > 0.0) ? tol1 : -tol1);
        const double fu = evaluate(u);
        ++result.iterations;

        if (fu <= fx)
        {
            if (u >= x) a = x; else b = x;
            v = w; fv = fw;
            w = x; fw = fx;
            x = u; fx = fu;
        }
        else
        {
            if (u < x) a = u; else b = u;
            if (fu <= fw || w == x)
            {
                v = w; fv = fw;
                w = u; fw = fu;
            }
            else if (fu <= fv || v == x || v == w)
            {
                v = u; fv = fu;
            }
        }
    }

    result.x = x;
    result.f = fx;
    return result;
}

// -----------------------------------------------------------------------------
// Bracketed root of a one-dimensional function: Brent's method.
//
// Brent, R. P. (1973). Algorithms for Minimization without Derivatives,
// Prentice-Hall, ch. 4 ("zero"); the same algorithm as Netlib's zeroin.
// Bisection, accelerated by secant and inverse quadratic interpolation when
// they stay inside the bracket. No derivatives.
//
// CONTRACT
//   * f(lower) and f(upper) must be finite and of opposite sign (or one of them
//     zero); otherwise bracketed = false and x is NaN.
//   * Converges to |x - root| <= 2 eps |x| + x_tolerance / 2.
//   * Throws std::invalid_argument for non-finite or reversed bounds, a
//     non-positive tolerance, or max_iterations < 1.
// -----------------------------------------------------------------------------
struct brent_root
{
    double x = std::numeric_limits<double>::quiet_NaN();
    double f = std::numeric_limits<double>::quiet_NaN();
    int iterations = 0;
    bool bracketed = false;
    bool converged = false;
};

template <class Function>
brent_root find_root_brent(
    Function&& f,
    double lower,
    double upper,
    double x_tolerance = 1.0e-12,
    int max_iterations = 200)
{
    if (!(std::isfinite(lower) && std::isfinite(upper)))
        throw std::invalid_argument("find_root_brent: bounds must be finite");
    if (lower > upper)
        throw std::invalid_argument("find_root_brent: lower exceeds upper");
    if (!(x_tolerance > 0.0))
        throw std::invalid_argument("find_root_brent: x_tolerance must be positive");
    if (max_iterations < 1)
        throw std::invalid_argument("find_root_brent: max_iterations must be >= 1");

    brent_root out;
    double a = lower, b = upper;
    double fa = static_cast<double>(f(a)), fb = static_cast<double>(f(b));
    if (!std::isfinite(fa) || !std::isfinite(fb) || (fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0))
        return out;
    out.bracketed = true;
    if (fa == 0.0) { out.x = a; out.f = 0.0; out.converged = true; return out; }
    if (fb == 0.0) { out.x = b; out.f = 0.0; out.converged = true; return out; }
    const double eps = std::numeric_limits<double>::epsilon();
    double c = a, fc = fa, d = b - a, e = d;
    for (int it = 0; it < max_iterations; ++it)
    {
        if ((fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0)) { c = a; fc = fa; d = b - a; e = d; }
        if (std::fabs(fc) < std::fabs(fb))
        {
            a = b; b = c; c = a;
            fa = fb; fb = fc; fc = fa;
        }
        const double tol = 2.0 * eps * std::fabs(b) + 0.5 * x_tolerance;
        const double m = 0.5 * (c - b);
        if (std::fabs(m) <= tol || fb == 0.0)
        {
            out.converged = true;
            break;
        }
        if (std::fabs(e) >= tol && std::fabs(fa) > std::fabs(fb))
        {
            double p, q, r;
            const double s = fb / fa;
            if (a == c)
            {
                p = 2.0 * m * s;
                q = 1.0 - s;
            }
            else
            {
                q = fa / fc;
                r = fb / fc;
                p = s * (2.0 * m * q * (q - r) - (b - a) * (r - 1.0));
                q = (q - 1.0) * (r - 1.0) * (s - 1.0);
            }
            if (p > 0.0) q = -q; else p = -p;
            if (2.0 * p < std::min(3.0 * m * q - std::fabs(tol * q), std::fabs(e * q)))
            {
                e = d;
                d = p / q;
            }
            else
            {
                d = m;
                e = m;
            }
        }
        else
        {
            d = m;
            e = m;
        }
        a = b;
        fa = fb;
        b += (std::fabs(d) > tol) ? d : (m > 0.0 ? tol : -tol);
        fb = static_cast<double>(f(b));
        ++out.iterations;
        if (!std::isfinite(fb)) { out.converged = false; break; }
    }
    out.x = b;
    out.f = fb;
    return out;
}

} // namespace lsdtt
