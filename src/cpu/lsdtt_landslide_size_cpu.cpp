// lsdtt_landslide_size_cpu.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
//
// See lsdtt_landslide_size_cpu.hpp for the design and the sources.

#include "lsdtt_crn_units.hpp"
#include "lsdtt_landslide_size_cpu.hpp"

#include <string>

#include <algorithm>
#include <cmath>
#include <limits>
#include <sstream>
#include <stdexcept>

namespace lsdtt
{

namespace
{

constexpr double kInf = std::numeric_limits<double>::infinity();

// Regularised LOWER incomplete gamma P(a,x), by series for x < a+1 and by
// continued fraction otherwise. Written out rather than taken from a library
// because the only alternative in the toolchain is std::tgamma, which gives the
// complete function only. Numerical Recipes formulation; the two branches meet
// at x = a+1 where both converge quickly.
double gamma_p(double a, double x)
{
    if (!(a > 0.0) || x < 0.0)
        return std::numeric_limits<double>::quiet_NaN();
    if (x == 0.0)
        return 0.0;

    const double gln = std::lgamma(a);

    if (x < a + 1.0)
    {
        // Series representation.
        double ap  = a;
        double sum = 1.0 / a;
        double del = sum;
        for (int n = 0; n < 500; ++n)
        {
            ap  += 1.0;
            del *= x / ap;
            sum += del;
            if (std::fabs(del) < std::fabs(sum) * 1.0e-15)
                break;
        }
        return sum * std::exp(-x + a * std::log(x) - gln);
    }

    // Continued fraction for Q(a,x); return 1 - Q.
    const double tiny = 1.0e-300;
    double b = x + 1.0 - a;
    double c = 1.0 / tiny;
    double d = 1.0 / b;
    double h = d;
    for (int i = 1; i <= 500; ++i)
    {
        const double an = -i * (i - a);
        b += 2.0;
        d  = an * d + b;
        if (std::fabs(d) < tiny) d = tiny;
        c  = b + an / c;
        if (std::fabs(c) < tiny) c = tiny;
        d  = 1.0 / d;
        const double del = d * c;
        h *= del;
        if (std::fabs(del - 1.0) < 1.0e-15)
            break;
    }
    const double q = std::exp(-x + a * std::log(x) - gln) * h;
    return 1.0 - q;
}

// Adaptive Simpson on [lo,hi]. The densities here are smooth but span many
// orders of magnitude, so the integration is done in LOG space by the callers.
double simpson(double (*f)(double, const void*), const void* ctx,
               double lo, double hi, int n_intervals)
{
    if (n_intervals % 2 != 0) ++n_intervals;
    const double h = (hi - lo) / n_intervals;
    double sum = f(lo, ctx) + f(hi, ctx);
    for (int i = 1; i < n_intervals; ++i)
    {
        const double x = lo + i * h;
        sum += (i % 2 == 0 ? 2.0 : 4.0) * f(x, ctx);
    }
    return sum * h / 3.0;
}

struct IgCtx { double rho, a, s, gamma_av, alpha_av; };

// Integrands in log space: int f(A) dA = int f(e^u) e^u du
double ig_number_logspace(double u, const void* vctx)
{
    const auto* c = static_cast<const IgCtx*>(vctx);
    const double A = std::exp(u);
    return inverse_gamma_pdf(A, c->rho, c->a, c->s) * A;
}
double ig_volume_logspace(double u, const void* vctx)
{
    const auto* c = static_cast<const IgCtx*>(vctx);
    const double A = std::exp(u);
    return c->alpha_av * std::pow(A, c->gamma_av)
         * inverse_gamma_pdf(A, c->rho, c->a, c->s) * A;
}

}  // namespace


// -----------------------------------------------------------------------------
// Model name helpers
// -----------------------------------------------------------------------------
std::string to_string(LandslideSizeModel model)
{
    switch (model)
    {
        case LandslideSizeModel::InverseGamma:    return "inverse_gamma";
        case LandslideSizeModel::TruncatedPareto: return "truncated_pareto";
    }
    return "unknown";
}

bool parse_landslide_size_model(const std::string& name, LandslideSizeModel& out)
{
    if (name == "inverse_gamma")    { out = LandslideSizeModel::InverseGamma;    return true; }
    if (name == "truncated_pareto") { out = LandslideSizeModel::TruncatedPareto; return true; }
    return false;
}


// -----------------------------------------------------------------------------
// Area-to-volume presets (Larsen et al. 2010 Supplementary Table S1)
// -----------------------------------------------------------------------------
// Densities: soil/regolith scars ~1500 kg/m3 (dry-loose 1300 to moist 1700);
// bedrock 2650. The global and Hovius rows are fitted across materials and
// carry no defensible density, so they are flagged rather than guessed.
AreaVolumeScaling AreaVolumeScaling::larsen_global()
{ return {0.1459, 1.332, 2650.0, false}; }   // log a -0.836, n=4231, R2 0.95
AreaVolumeScaling AreaVolumeScaling::larsen_soil_scar()
{ return {0.2244, 1.262, 1500.0, true}; }    // log a -0.649, n=1617, R2 0.92
AreaVolumeScaling AreaVolumeScaling::larsen_bedrock_scar()
{ return {0.2344, 1.410, 2650.0, true}; }    // log a -0.630, n=168,  R2 0.97
AreaVolumeScaling AreaVolumeScaling::larsen_soil_all()
{ return {0.3631, 1.145, 1500.0, true}; }    // log a -0.440, n=2136 (scar+deposit)
AreaVolumeScaling AreaVolumeScaling::larsen_bedrock_all()
{ return {0.1862, 1.350, 2650.0, true}; }    // log a -0.730, n=604  (scar+deposit)
AreaVolumeScaling AreaVolumeScaling::hovius_slab()
{ return {0.05,  1.5, 2650.0, false}; }      // d = 0.05*sqrt(A), uniform depth
AreaVolumeScaling AreaVolumeScaling::hovius_parabolic()
{ return {0.025, 1.5, 2650.0, false}; }      // Niemi Eq. A12, parabolic scar

bool AreaVolumeScaling::from_name(const std::string& name, AreaVolumeScaling& out)
{
    if (name == "larsen_global")       { out = larsen_global();       return true; }
    if (name == "larsen_soil_scar")    { out = larsen_soil_scar();    return true; }
    if (name == "larsen_bedrock_scar") { out = larsen_bedrock_scar(); return true; }
    if (name == "larsen_soil_all")     { out = larsen_soil_all();     return true; }
    if (name == "larsen_bedrock_all")  { out = larsen_bedrock_all();  return true; }
    if (name == "hovius_slab")         { out = hovius_slab();         return true; }
    if (name == "hovius_parabolic")    { out = hovius_parabolic();    return true; }
    return false;
}

double AreaVolumeScaling::volume_m3(double area_m2) const
{
    if (!(area_m2 > 0.0)) return 0.0;
    return alpha * std::pow(area_m2, gamma);
}

double AreaVolumeScaling::scar_effective_depth_g_cm2(double depth_m) const
{
    return metres_to_g_cm2(depth_m, scar_density_kg_m3);
}

double AreaVolumeScaling::rock_equivalent_depth_m(double depth_m,
                                                  double rock_density_kg_m3) const
{
    if (!(rock_density_kg_m3 > 0.0)) return depth_m;
    return depth_m * scar_density_kg_m3 / rock_density_kg_m3;
}

std::string area_volume_density_note(const AreaVolumeScaling& av)
{
    if (av.density_is_material_specific) return {};
    return "the selected area-volume preset is fitted across BOTH soil and "
           "bedrock (or carries no material at all), so it has no defensible "
           "scar density; " + std::to_string(static_cast<int>(av.scar_density_kg_m3))
         + " kg/m3 is a placeholder. Set the scar density explicitly, or pick a "
           "material-specific preset (larsen_soil_scar / larsen_bedrock_scar). "
           "A soil scar carried at rock density overstates its mass depth by "
           "about a factor of two.";
}


// -----------------------------------------------------------------------------
// Continuous densities
// -----------------------------------------------------------------------------
double inverse_gamma_pdf(double area_m2, double rho, double a_param, double s_param)
{
    const double x = area_m2 - s_param;
    if (!(x > 0.0) || !(a_param > 0.0) || !(rho > 0.0))
        return 0.0;
    // 1/(a G(rho)) * (a/x)^(rho+1) * exp(-a/x), evaluated in logs: the
    // (a/x)^(rho+1) factor underflows for large A at double precision.
    const double log_term = -std::log(a_param) - std::lgamma(rho)
                          + (rho + 1.0) * (std::log(a_param) - std::log(x))
                          - a_param / x;
    return std::exp(log_term);
}

double inverse_gamma_survival(double area_m2, double rho, double a_param, double s_param)
{
    const double x = area_m2 - s_param;
    if (!(x > 0.0)) return 1.0;
    // For X ~ InvGamma(rho, a):  P(X <= x) = Q(rho, a/x) = 1 - P(rho, a/x),
    // so the survival P(X > x) is the regularised LOWER incomplete gamma.
    return gamma_p(rho, a_param / x);
}

double inverse_gamma_mean_area_m2(double rho, double a_param, double s_param)
{
    if (!(rho > 1.0)) return kInf;   // first moment diverges
    return a_param / (rho - 1.0) + s_param;
}

double truncated_pareto_pdf(double area_m2, double alpha_pdf,
                            double a_min_m2, double a_max_m2)
{
    if (area_m2 < a_min_m2 || area_m2 > a_max_m2) return 0.0;
    if (!(a_min_m2 > 0.0) || !(a_max_m2 > a_min_m2)) return 0.0;

    if (std::fabs(alpha_pdf - 1.0) < 1.0e-9)
    {
        // Normalisation becomes 1/ln(a_max/a_min).
        return 1.0 / (area_m2 * std::log(a_max_m2 / a_min_m2));
    }
    const double e = 1.0 - alpha_pdf;
    const double norm = e / (std::pow(a_max_m2, e) - std::pow(a_min_m2, e));
    return norm * std::pow(area_m2, -alpha_pdf);
}

double truncated_pareto_sample(double u01, double alpha_pdf,
                               double a_min_m2, double a_max_m2)
{
    u01 = std::clamp(u01, 0.0, 1.0);
    if (std::fabs(alpha_pdf - 1.0) < 1.0e-9)
        return a_min_m2 * std::pow(a_max_m2 / a_min_m2, u01);

    const double e   = 1.0 - alpha_pdf;
    const double lo  = std::pow(a_min_m2, e);
    const double hi  = std::pow(a_max_m2, e);
    return std::pow(lo + u01 * (hi - lo), 1.0 / e);
}

double truncated_pareto_mean_area_m2(double alpha_pdf,
                                     double a_min_m2, double a_max_m2)
{
    if (!(a_min_m2 > 0.0) || !(a_max_m2 > a_min_m2)) return 0.0;

    // alpha == 2 makes the numerator exponent zero -> logarithm.
    if (std::fabs(alpha_pdf - 2.0) < 1.0e-9)
    {
        const double e = 1.0 - alpha_pdf;   // = -1
        return std::log(a_max_m2 / a_min_m2)
             * (e / (std::pow(a_max_m2, e) - std::pow(a_min_m2, e)));
    }
    if (std::fabs(alpha_pdf - 1.0) < 1.0e-9)
        return (a_max_m2 - a_min_m2) / std::log(a_max_m2 / a_min_m2);

    const double e1 = 1.0 - alpha_pdf;
    const double e2 = 2.0 - alpha_pdf;
    return (e1 / e2) * (std::pow(a_max_m2, e2) - std::pow(a_min_m2, e2))
                     / (std::pow(a_max_m2, e1) - std::pow(a_min_m2, e1));
}


// -----------------------------------------------------------------------------
// Validation
// -----------------------------------------------------------------------------
bool LandslideSizeParams::validate(double cell_size_m, std::string& error) const
{
    std::ostringstream os;

    if (!(cell_size_m > 0.0))
    { os << "cell size must be > 0 (got " << cell_size_m << ")"; error = os.str(); return false; }

    if (n_min_pixels < 1)
    { os << "n_min_pixels must be >= 1 (got " << n_min_pixels << ")"; error = os.str(); return false; }

    if (!(area_max_m2 > 0.0))
    { os << "area_max_m2 must be > 0 (got " << area_max_m2 << ")"; error = os.str(); return false; }

    const double cell_area = cell_size_m * cell_size_m;
    const double n_max = std::floor(area_max_m2 / cell_area);
    if (n_max < static_cast<double>(n_min_pixels))
    {
        os << "area_max_m2 (" << area_max_m2 << " m2 = " << n_max
           << " pixels at " << cell_size_m << " m) is smaller than n_min_pixels ("
           << n_min_pixels << "). A landslide cannot be smaller than its minimum.";
        error = os.str();
        return false;
    }

    if (model == LandslideSizeModel::InverseGamma)
    {
        if (!(rho > 1.0))
        {
            os << "inverse gamma rho must be > 1 for a finite mean area (got "
               << rho << ")";
            error = os.str(); return false;
        }
        if (!(a_param > 0.0))
        { os << "inverse gamma a_param must be > 0 (got " << a_param << ")"; error = os.str(); return false; }
        // s is expected to be negative (Malamud: -132 m2) but a positive s is
        // legal -- it just shifts the support. Only s >= smallest bin is fatal.
        if (s_param >= cell_area * n_min_pixels)
        {
            os << "inverse gamma s_param (" << s_param
               << " m2) must be below the smallest landslide bin ("
               << cell_area * n_min_pixels << " m2); the density is zero at or below s";
            error = os.str(); return false;
        }
    }
    else
    {
        if (!(pareto_alpha_pdf > 0.0))
        {
            os << "pareto_alpha_pdf must be > 0 (got " << pareto_alpha_pdf
               << "). NOTE this is the PDF exponent, not the cumulative one: "
                  "alpha_pdf = b_cumulative + 1.";
            error = os.str(); return false;
        }
    }

    if (!(av.gamma > 0.0))
    { os << "area-volume gamma must be > 0 (got " << av.gamma << ")"; error = os.str(); return false; }
    if (!(av.alpha > 0.0))
    { os << "area-volume alpha must be > 0 (got " << av.alpha << ")"; error = os.str(); return false; }

    error.clear();
    return true;
}


// -----------------------------------------------------------------------------
// Sampler
// -----------------------------------------------------------------------------
LandslideSizeSampler::LandslideSizeSampler(const LandslideSizeParams& params,
                                           double cell_size_m)
  : params_(params)
{
    std::string err;
    if (!params_.validate(cell_size_m, err))
        throw std::invalid_argument("LandslideSizeSampler: " + err);

    cell_area_m2_ = cell_size_m * cell_size_m;
    n_min_        = params_.n_min_pixels;
    n_max_        = static_cast<int>(std::floor(params_.area_max_m2 / cell_area_m2_));

    const int n_bins = n_max_ - n_min_ + 1;
    std::vector<double> w(static_cast<std::size_t>(n_bins), 0.0);

    // Bin weights: integrate the continuous density across each pixel-count
    // bin, [ (n-1/2)*cell, (n+1/2)*cell ). This is the step that makes the
    // sampler correct at coarse resolution -- sampling continuously and
    // rounding would misweight the small bins, which is where the rollover is.
    const double a_min_edge = (static_cast<double>(n_min_) - 0.5) * cell_area_m2_;
    const double a_max_edge = (static_cast<double>(n_max_) + 0.5) * cell_area_m2_;

    for (int i = 0; i < n_bins; ++i)
    {
        const int    n  = n_min_ + i;
        const double lo = std::max((static_cast<double>(n) - 0.5) * cell_area_m2_, 1.0e-9);
        const double hi = (static_cast<double>(n) + 0.5) * cell_area_m2_;

        if (params_.model == LandslideSizeModel::InverseGamma)
        {
            IgCtx ctx{params_.rho, params_.a_param, params_.s_param,
                      params_.av.gamma, params_.av.alpha};
            // Integrate in log space: bins span orders of magnitude at the
            // large end, and the density is far smoother in log A.
            w[static_cast<std::size_t>(i)] =
                simpson(ig_number_logspace, &ctx, std::log(lo), std::log(hi), 8);
        }
        else
        {
            // Analytic bin mass from the Pareto CDF -- no quadrature needed.
            const double amin = a_min_edge;
            const double amax = a_max_edge;
            const double e    = 1.0 - params_.pareto_alpha_pdf;
            if (std::fabs(params_.pareto_alpha_pdf - 1.0) < 1.0e-9)
            {
                w[static_cast<std::size_t>(i)] =
                    std::log(hi / lo) / std::log(amax / amin);
            }
            else
            {
                const double denom = std::pow(amax, e) - std::pow(amin, e);
                w[static_cast<std::size_t>(i)] =
                    (std::pow(hi, e) - std::pow(lo, e)) / denom;
            }
        }
        if (!(w[static_cast<std::size_t>(i)] > 0.0) ||
            !std::isfinite(w[static_cast<std::size_t>(i)]))
            w[static_cast<std::size_t>(i)] = 0.0;
    }

    // Normalise and accumulate.
    double total = 0.0;
    for (double v : w) total += v;
    if (!(total > 0.0))
        throw std::invalid_argument(
            "LandslideSizeSampler: the distribution has zero mass over "
            "[n_min, n_max]; check the parameters against the cell size");

    cdf_.resize(static_cast<std::size_t>(n_bins));
    double acc = 0.0;
    for (int i = 0; i < n_bins; ++i)
    {
        acc += w[static_cast<std::size_t>(i)] / total;
        cdf_[static_cast<std::size_t>(i)] = acc;
    }
    cdf_.back() = 1.0;   // guard against accumulated round-off

    // Moments of the TRUNCATED, DISCRETISED distribution. Computed from the pmf
    // rather than in closed form: it is one pass, it is exact for what is
    // actually sampled, and it removes a whole class of algebra error.
    mean_area_m2_   = 0.0;
    mean_volume_m3_ = 0.0;
    for (int i = 0; i < n_bins; ++i)
    {
        const double p = w[static_cast<std::size_t>(i)] / total;
        const double A = static_cast<double>(n_min_ + i) * cell_area_m2_;
        mean_area_m2_   += p * A;
        mean_volume_m3_ += p * params_.av.volume_m3(A);
    }

    // How much of the UNTRUNCATED distribution lies above area_max, by number
    // and by volume. The volume fraction is expected to be large; it is
    // reported so the truncation is visible rather than silently absorbed.
    if (params_.model == LandslideSizeModel::InverseGamma)
    {
        trunc_frac_n_ = inverse_gamma_survival(a_max_edge, params_.rho,
                                               params_.a_param, params_.s_param);

        IgCtx ctx{params_.rho, params_.a_param, params_.s_param,
                  params_.av.gamma, params_.av.alpha};
        // 1e14 m^2 stands in for infinity: the volume integrand decays as
        // A^(gamma-rho-1), which for the published parameters is ~A^-1.068, so
        // the integral converges only marginally and a finite ceiling must be
        // stated rather than pretended away.
        const double lo_all = std::log(1.0e-3);
        const double hi_all = std::log(1.0e14);
        const double v_all  = simpson(ig_volume_logspace, &ctx, lo_all, hi_all, 4000);
        const double v_above = simpson(ig_volume_logspace, &ctx,
                                       std::log(a_max_edge), hi_all, 2000);
        trunc_frac_v_ = (v_all > 0.0) ? std::clamp(v_above / v_all, 0.0, 1.0) : 0.0;
    }
    else
    {
        // The truncated Pareto IS the distribution -- there is nothing above
        // a_max by construction.
        trunc_frac_n_ = 0.0;
        trunc_frac_v_ = 0.0;
    }

    // Can this cell size represent the rollover at all? The mode of the inverse
    // gamma sits near a/(rho+1) + s. If the smallest representable landslide is
    // above that, every sampled area is on the power-law tail and the choice of
    // family is inoperative.
    if (params_.model == LandslideSizeModel::InverseGamma)
    {
        const double mode = params_.a_param / (params_.rho + 1.0) + params_.s_param;
        rollover_unresolved_ =
            (static_cast<double>(n_min_) * cell_area_m2_ >= mode);
    }
}

int LandslideSizeSampler::sample_pixels(double u01) const
{
    u01 = std::clamp(u01, 0.0, std::nextafter(1.0, 0.0));
    const auto it = std::lower_bound(cdf_.begin(), cdf_.end(), u01);
    const auto idx = (it == cdf_.end())
        ? cdf_.size() - 1
        : static_cast<std::size_t>(it - cdf_.begin());
    return n_min_ + static_cast<int>(idx);
}

double LandslideSizeSampler::sample_area_m2(double u01) const
{
    return static_cast<double>(sample_pixels(u01)) * cell_area_m2_;
}

double LandslideSizeSampler::pmf(int n_pixels) const
{
    if (n_pixels < n_min_ || n_pixels > n_max_) return 0.0;
    const std::size_t i = static_cast<std::size_t>(n_pixels - n_min_);
    return (i == 0) ? cdf_[0] : cdf_[i] - cdf_[i - 1];
}

}  // namespace lsdtt
