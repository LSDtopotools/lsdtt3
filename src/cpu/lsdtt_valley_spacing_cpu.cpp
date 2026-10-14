// =============================================================================
// lsdtt_valley_spacing_cpu.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// See the header for the three measures and why all three are returned.
// =============================================================================

#include "lsdtt_valley_spacing_cpu.hpp"
#include "../utils/lsdtt_grid_neighbours.hpp"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <sstream>
#include <vector>

#ifdef _OPENMP
#  include <omp.h>
#endif

namespace lsdtt
{

// -----------------------------------------------------------------------------
double channelization_index(double K, double D, double U,
                            double domain_length_m, double m, double n)
{
    if (!(K > 0.0) || !(D > 0.0) || !(domain_length_m > 0.0)) return 0.0;

    // CI = K l^(m+n) / (D^n U^(1-n)). For n == 1 the U term vanishes, which is
    // why Bonetti et al. note the index is uplift-independent there.
    const double num = K * std::pow(domain_length_m, m + n);
    const double den = std::pow(D, n) *
                       ((std::fabs(n - 1.0) < 1e-12) ? 1.0
                                                     : std::pow(U, 1.0 - n));
    return (den > 0.0) ? num / den : 0.0;
}


namespace
{

// A direct real DFT of one line, accumulating |X_k|^2.
//
// Deliberately not FFTW: FFTW is an OPTIONAL dependency here, and a 1-D
// transform of a few hundred samples is cheap enough that carrying a conditional
// code path would cost more in confusion than it saves in time. It is also exact
// on a periodic domain, where zero-padding to a power of two would not be.
void accumulate_power(const std::vector<double>& line, std::vector<double>& acc)
{
    const std::size_t n  = line.size();
    const std::size_t nf = acc.size();
    const double two_pi_over_n = 2.0 * 3.14159265358979323846 / static_cast<double>(n);

    for (std::size_t k = 0; k < nf; ++k)
    {
        double re = 0.0, im = 0.0;
        for (std::size_t i = 0; i < n; ++i)
        {
            const double ang = two_pi_over_n * static_cast<double>(k * i % n);
            re += line[i] * std::cos(ang);
            im -= line[i] * std::sin(ang);
        }
        acc[k] += (re * re + im * im) / static_cast<double>(n * n);
    }
}

} // namespace


// -----------------------------------------------------------------------------
lsdtt_valley_spacing_result compute_valley_spacing(
    const lsdtt_raster_cpu&             dem,
    const lsdtt_valley_spacing_options& opts)
{
    lsdtt_valley_spacing_result out;

    const lsdtt_raster_meta& meta = dem.get_meta();
    const int nrows = meta.n_rows;
    const int ncols = meta.n_cols;
    const double dx = meta.cell_size;
    const std::vector<float>&   z    = dem.get_data();
    const std::vector<uint8_t>& mask = dem.get_mask();

    const bool along_rows =
        (opts.axis == lsdtt_valley_spacing_options::measurement_axis::along_rows);

    // n_along  = samples per measurement line
    // n_lines  = number of lines, i.e. the flow direction
    const int n_along = along_rows ? ncols : nrows;
    const int n_lines = along_rows ? nrows : ncols;

    if (n_along < 8 || n_lines < 8)
    {
        out.warning = "raster too small for a valley-spacing measurement";
        return out;
    }

    const int margin = std::max(0, std::min(opts.margin_cells, n_lines / 2 - 2));
    const int lo = margin, hi = n_lines - margin;

    auto at = [&](int line, int i) -> std::size_t
    {
        return along_rows
            ? static_cast<std::size_t>(line) * ncols + i
            : static_cast<std::size_t>(i) * ncols + line;
    };

    out.domain_length_m = (opts.domain_length_m > 0.0)
        ? opts.domain_length_m
        : static_cast<double>(along_rows ? nrows : ncols) * dx;

    // -------------------------------------------------------------------------
    // Spectral measure.
    //
    // Each line is detrended by its own mean -- that is the "remove the profile
    // in the flow direction" step: it strips the base-level gradient, which
    // otherwise dominates the spectrum, without touching the across-flow
    // structure the measurement is about.
    // -------------------------------------------------------------------------
    const std::size_t nf = static_cast<std::size_t>(n_along / 2) + 1;
    std::vector<double> acc(nf, 0.0);
    long n_lines_used = 0;

    for (int line = lo; line < hi; ++line)
    {
        std::vector<double> v(static_cast<std::size_t>(n_along));
        bool ok = true;
        double sum = 0.0;
        for (int i = 0; i < n_along; ++i)
        {
            const std::size_t f = at(line, i);
            if (!mask[f]) { ok = false; break; }        // a gap leaks everywhere
            v[static_cast<std::size_t>(i)] = static_cast<double>(z[f]);
            sum += v[static_cast<std::size_t>(i)];
        }
        if (!ok) continue;

        const double mean = sum / static_cast<double>(n_along);
        for (double& x : v) x -= mean;

        accumulate_power(v, acc);
        ++n_lines_used;
    }

    if (n_lines_used == 0)
    {
        out.warning = "no complete lines free of nodata; cannot form a spectrum";
        return out;
    }

    for (double& p : acc) p /= static_cast<double>(n_lines_used);
    acc[0] = 0.0;                                        // DC carries no spacing

    out.frequency_cycles_per_m.resize(nf);
    for (std::size_t k = 0; k < nf; ++k)
        out.frequency_cycles_per_m[k] =
            static_cast<double>(k) / (static_cast<double>(n_along) * dx);
    out.power = acc;

    // Local maxima only: Perron et al. average over PEAKS, not over the spectrum.
    std::vector<std::size_t> peaks;
    for (std::size_t k = 1; k + 1 < nf; ++k)
        if (acc[k] > acc[k - 1] && acc[k] > acc[k + 1]) peaks.push_back(k);

    out.n_spectral_peaks = static_cast<int>(peaks.size());

    if (!peaks.empty())
    {
        double wsum = 0.0, fsum = 0.0, best = -1.0;
        std::size_t kbest = peaks.front();
        for (std::size_t k : peaks)
        {
            wsum += acc[k];
            fsum += out.frequency_cycles_per_m[k] * acc[k];
            if (acc[k] > best) { best = acc[k]; kbest = k; }
        }
        const double f_bar = (wsum > 0.0) ? fsum / wsum : 0.0;
        out.perron_wavelength_m = (f_bar > 0.0) ? 1.0 / f_bar : 0.0;
        const double f_max = out.frequency_cycles_per_m[kbest];
        out.peak_wavelength_m = (f_max > 0.0) ? 1.0 / f_max : 0.0;
    }

    // -------------------------------------------------------------------------
    // Convex half-width: the robust measure.
    // -------------------------------------------------------------------------
    long   n_runs = 0;
    double run_total = 0.0;
    for (int line = lo; line < hi; ++line)
    {
        int cur = 0;
        for (int i = 1; i + 1 < n_along; ++i)
        {
            const std::size_t f  = at(line, i);
            const std::size_t fm = at(line, i - 1);
            const std::size_t fp = at(line, i + 1);
            if (!mask[f] || !mask[fm] || !mask[fp])
            {
                if (cur) { run_total += cur; ++n_runs; cur = 0; }
                continue;
            }
            const double curv = (static_cast<double>(z[fm]) + static_cast<double>(z[fp])
                               - 2.0 * static_cast<double>(z[f])) / (dx * dx);
            if (curv < 0.0) ++cur;                        // convex = hillslope
            else if (cur)   { run_total += cur; ++n_runs; cur = 0; }
        }
        if (cur) { run_total += cur; ++n_runs; }
    }

    out.n_convex_runs = n_runs;
    if (n_runs > 0)
    {
        out.convex_halfwidth_m  = (run_total / static_cast<double>(n_runs)) * dx;
        out.convex_wavelength_m = 2.0 * out.convex_halfwidth_m;
    }

    // -------------------------------------------------------------------------
    // Hillslope length by distance transform -- the branching-robust measure.
    //
    // Two-pass chamfer distance (3-4 weights, ~2% of true Euclidean, which is far
    // inside the scatter of the quantity itself) seeded on the concave cells.
    // Orientation-free, so oblique tributaries cost it nothing.
    // -------------------------------------------------------------------------
    {
        constexpr double k_big = 1e30;
        const double wa = 1.0 * dx;         // orthogonal step
        const double wb = grid::k_sqrt2_d * dx;  // diagonal step

        std::vector<double> dist(static_cast<std::size_t>(nrows) * ncols, k_big);
        std::vector<uint8_t> is_hill(dist.size(), 0);

        auto curv2d = [&](int r, int c) -> double
        {
            const std::size_t f  = static_cast<std::size_t>(r) * ncols + c;
            const std::size_t fn = f - ncols, fs = f + ncols, fw = f - 1, fe = f + 1;
            if (!mask[f] || !mask[fn] || !mask[fs] || !mask[fw] || !mask[fe])
                return 0.0;
            return (static_cast<double>(z[fn]) + static_cast<double>(z[fs])
                  + static_cast<double>(z[fw]) + static_cast<double>(z[fe])
                  - 4.0 * static_cast<double>(z[f])) / (dx * dx);
        };

        for (int r = 1; r < nrows - 1; ++r)
            for (int c = 1; c < ncols - 1; ++c)
            {
                const std::size_t f = static_cast<std::size_t>(r) * ncols + c;
                if (!mask[f]) continue;
                const double k = curv2d(r, c);
                if (k > 0.0) dist[f] = 0.0;           // concave: a valley cell
                else if (k < 0.0) is_hill[f] = 1;     // convex: a hillslope cell
            }

        auto relax = [&](std::size_t f, std::size_t g, double w)
        {
            if (dist[g] + w < dist[f]) dist[f] = dist[g] + w;
        };

        for (int r = 1; r < nrows - 1; ++r)               // forward pass
            for (int c = 1; c < ncols - 1; ++c)
            {
                const std::size_t f = static_cast<std::size_t>(r) * ncols + c;
                if (!mask[f]) continue;
                relax(f, f - ncols, wa); relax(f, f - 1, wa);
                relax(f, f - ncols - 1, wb); relax(f, f - ncols + 1, wb);
            }
        for (int r = nrows - 2; r >= 1; --r)              // backward pass
            for (int c = ncols - 2; c >= 1; --c)
            {
                const std::size_t f = static_cast<std::size_t>(r) * ncols + c;
                if (!mask[f]) continue;
                relax(f, f + ncols, wa); relax(f, f + 1, wa);
                relax(f, f + ncols + 1, wb); relax(f, f + ncols - 1, wb);
            }

        std::vector<double> d;
        d.reserve(dist.size() / 2);
        for (int r = lo; r < hi && along_rows; ++r)
            for (int c = 1; c < ncols - 1; ++c)
            {
                const std::size_t f = static_cast<std::size_t>(r) * ncols + c;
                if (is_hill[f] && dist[f] < k_big) d.push_back(dist[f]);
            }
        if (!along_rows)
            for (int r = 1; r < nrows - 1; ++r)
                for (int c = lo; c < hi; ++c)
                {
                    const std::size_t f = static_cast<std::size_t>(r) * ncols + c;
                    if (is_hill[f] && dist[f] < k_big) d.push_back(dist[f]);
                }

        out.n_hillslope_cells = static_cast<long>(d.size());
        if (!d.empty())
        {
            out.hillslope_length_mean_m =
                std::accumulate(d.begin(), d.end(), 0.0) / static_cast<double>(d.size());
            std::sort(d.begin(), d.end());
            out.hillslope_length_median_m = d[d.size() / 2];
            out.hillslope_length_p90_m    = d[static_cast<std::size_t>(0.9 * (d.size() - 1))];
        }
    }

    // -------------------------------------------------------------------------
    // Spectral exponent alpha, fitted BELOW the measured valley wavelength.
    // -------------------------------------------------------------------------
    const double lam = (out.convex_wavelength_m > 0.0) ? out.convex_wavelength_m
                                                       : out.perron_wavelength_m;
    if (lam > 0.0)
    {
        const double lam_hi = opts.alpha_fit_max_wavelength_frac * lam;
        const double lam_lo = opts.alpha_fit_min_wavelength_frac * lam;
        std::vector<double> lx, ly;
        for (std::size_t k = 1; k < nf; ++k)
        {
            const double f = out.frequency_cycles_per_m[k];
            if (f <= 0.0 || acc[k] <= 0.0) continue;
            const double w = 1.0 / f;
            if (w <= lam_hi && w >= lam_lo)
            {
                lx.push_back(std::log(f));
                ly.push_back(std::log(acc[k]));
            }
        }
        if (lx.size() >= 5)
        {
            const double n_pt = static_cast<double>(lx.size());
            const double sx  = std::accumulate(lx.begin(), lx.end(), 0.0);
            const double sy  = std::accumulate(ly.begin(), ly.end(), 0.0);
            double sxx = 0.0, sxy = 0.0;
            for (std::size_t i = 0; i < lx.size(); ++i)
            { sxx += lx[i] * lx[i]; sxy += lx[i] * ly[i]; }
            const double den = n_pt * sxx - sx * sx;
            if (std::fabs(den) > 0.0) out.spectral_exponent = (n_pt * sxy - sx * sy) / den;
        }
    }
    if (opts.m > 0.0) out.spectral_exponent_predicted = 2.0 * opts.m - 3.0;

    // -------------------------------------------------------------------------
    // Context and caveats.
    // -------------------------------------------------------------------------
    out.channelization_index = channelization_index(
        opts.K, opts.D, opts.U, out.domain_length_m, opts.m, opts.n);

    const double across_extent = static_cast<double>(n_along) * dx;
    if (out.convex_wavelength_m > 0.0)
        out.valleys_in_domain = across_extent / out.convex_wavelength_m;

    std::ostringstream w;
    if (out.valleys_in_domain > 0.0 && out.valleys_in_domain < 8.0)
    {
        out.spectrum_under_sampled = true;
        w << "only ~" << static_cast<int>(out.valleys_in_domain + 0.5)
          << " valleys fit across this domain, so the SPECTRAL wavelength can only "
             "land on domain harmonics (extent/n) and should not be read as a "
             "trend. The convex half-width is unaffected and is the measure to "
             "use here; a domain holding ~10 valleys would make both usable.";
    }
    if (out.channelization_index > 0.0 &&
        out.channelization_index < k_channelization_index_critical_high)
    {
        if (!w.str().empty()) w << " ";
        w << "CI = " << out.channelization_index << " is at or below the critical "
             "range " << k_channelization_index_critical_low << "-"
          << k_channelization_index_critical_high
          << " (Bonetti et al. 2020), so this landscape is at or before first "
             "channelization and regular valleys may not have formed at all.";
    }
    out.warning = w.str();

    return out;
}

} // namespace lsdtt
