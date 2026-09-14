// lsdtt_niemi_ensemble_cpu.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// See the header for what the two rates are and why the sampling is stratified.
// -----------------------------------------------------------------------------

#include "lsdtt_niemi_ensemble_cpu.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <stdexcept>

#include "lsdtt_cosmo_raster_cpu.hpp"

namespace lsdtt
{

namespace
{

double percentile(std::vector<double> v, double p)
{
    if (v.empty()) return 0.0;
    std::sort(v.begin(), v.end());
    if (v.size() == 1) return v[0];
    const double x = p * static_cast<double>(v.size() - 1);
    const std::size_t i = static_cast<std::size_t>(std::floor(x));
    const double frac = x - static_cast<double>(i);
    if (i + 1 >= v.size()) return v.back();
    return v[i] * (1.0 - frac) + v[i + 1] * frac;
}

}  // namespace


// -----------------------------------------------------------------------------
std::vector<std::int64_t> sample_basin_outlets_stratified(
    const lsdtt_flow_base& flow,
    std::int64_t           min_pixels,
    std::int64_t           max_pixels,
    int                    per_band,
    int                    bands_per_decade)
{
    if (!flow.has_s_vector())
        throw std::invalid_argument(
            "sample_basin_outlets_stratified: the flow object needs an s_vector.");
    if (per_band <= 0 || bands_per_decade <= 0)
        return {};
    if (min_pixels < 1) min_pixels = 1;

    // Collect every candidate with its upslope count, then band by log10(area).
    struct Cand { std::int64_t flat; std::int64_t n; };
    std::vector<Cand> cands;

    const std::int64_t n_total =
        static_cast<std::int64_t>(flow.n_rows()) * flow.n_cols();
    for (std::int64_t flat = 0; flat < n_total; ++flat)
    {
        if (!flow.is_valid_cell(flat)) continue;
        const std::int64_t n =
            static_cast<std::int64_t>(flow.s_vector_size_unchecked(flat));
        if (n < min_pixels || n > max_pixels) continue;
        cands.push_back({flat, n});
    }
    if (cands.empty()) return {};

    std::sort(cands.begin(), cands.end(),
              [](const Cand& a, const Cand& b) { return a.n < b.n; });

    const double lo = std::log10(static_cast<double>(cands.front().n));
    const double hi = std::log10(static_cast<double>(cands.back().n));
    const double band_w = 1.0 / static_cast<double>(bands_per_decade);
    const int n_bands = std::max(1,
        static_cast<int>(std::ceil((hi - lo) / band_w)) + 1);

    std::vector<std::vector<std::size_t>> by_band(
        static_cast<std::size_t>(n_bands));
    for (std::size_t i = 0; i < cands.size(); ++i)
    {
        int b = static_cast<int>(
            std::floor((std::log10(static_cast<double>(cands[i].n)) - lo) / band_w));
        b = std::clamp(b, 0, n_bands - 1);
        by_band[static_cast<std::size_t>(b)].push_back(i);
    }

    // Spread the picks through each band rather than taking the first few,
    // which would cluster them at the small end of the band. Nested basins are
    // unavoidable in a drainage network and are not filtered: a large basin
    // legitimately contains the small ones, and excluding them would bias the
    // large end towards whichever tributaries happened to be dropped.
    std::vector<std::int64_t> out;
    for (const auto& band : by_band)
    {
        if (band.empty()) continue;
        const int take = std::min<int>(per_band, static_cast<int>(band.size()));
        for (int k = 0; k < take; ++k)
        {
            const std::size_t idx = (take == 1)
                ? band.size() / 2
                : static_cast<std::size_t>(
                      static_cast<double>(k) * static_cast<double>(band.size() - 1)
                      / static_cast<double>(take - 1));
            out.push_back(cands[band[idx]].flat);
        }
    }
    std::sort(out.begin(), out.end());
    out.erase(std::unique(out.begin(), out.end()), out.end());
    return out;
}


// -----------------------------------------------------------------------------
std::vector<NiemiBasinSample> evaluate_niemi_basins(
    const lsdtt_flowinfo&            flow,
    const std::vector<std::int64_t>& outlet_flats,
    const lsdtt_raster_cpu&          concentration,
    const lsdtt_raster_cpu&          weight_raster,
    const lsdtt_raster_cpu&          true_rate_g_cm2_yr,
    const lsdtt_raster_cpu&          excavation_depth_m,
    const lsdtt_raster_cpu&          prod_raster,
    const lsdtt_raster_cpu&          topo_shield,
    const lsdtt_raster_cpu&          snow_shield,
    const lsdtt_raster_cpu&          self_shield,
    const std::string&               nuclide,
    const std::string&               muon_scheme,
    int                              realisation)
{
    const lsdtt_raster_meta& meta = concentration.meta();
    const double cell_area_m2 = meta.cell_size * meta.cell_size;

    const auto& cd = concentration.get_data_array();
    const auto& cm = concentration.get_mask_array();
    const auto& wd = weight_raster.get_data_array();
    const auto& wm = weight_raster.get_mask_array();
    const auto& td = true_rate_g_cm2_yr.get_data_array();
    const auto& tm = true_rate_g_cm2_yr.get_mask_array();
    const auto& ed = excavation_depth_m.get_data_array();
    const auto& em = excavation_depth_m.get_mask_array();
    const float endv = static_cast<float>(excavation_depth_m.meta().no_data_value);

    std::vector<NiemiBasinSample> out;
    out.reserve(outlet_flats.size());

    for (std::int64_t outlet : outlet_flats)
    {
        NiemiBasinSample s;
        s.realisation = realisation;
        s.outlet_flat = outlet;

        if (!flow.is_valid_cell(outlet)) { out.push_back(s); continue; }

        const lsdtt_flat_index_slice basin = flow.upslope_nodes(outlet);

        double num = 0.0, den = 0.0, true_sum = 0.0;
        std::int64_t n_true = 0;

        for (std::int64_t flat : basin)
        {
            const int r = flow.row_from_flat_unchecked(flat);
            const int c = flow.col_from_flat_unchecked(flat);
            ++s.n_pixels;

            if (em(r, c) && ed(r, c) != endv && ed(r, c) > 0.0f)
                ++s.n_failed_pixels;

            if (cm(r, c) && wm(r, c))
            {
                const double w = static_cast<double>(wd(r, c));
                num += static_cast<double>(cd(r, c)) * w;
                den += w;
            }
            if (tm(r, c))
            {
                true_sum += static_cast<double>(td(r, c));
                ++n_true;
            }
        }

        s.area_m2  = static_cast<double>(s.n_pixels) * cell_area_m2;
        s.area_km2 = s.area_m2 * 1.0e-6;
        s.fraction_failed = s.n_pixels > 0
            ? static_cast<double>(s.n_failed_pixels) / static_cast<double>(s.n_pixels)
            : 0.0;

        if (!(den > 0.0) || n_true == 0)
        {
            // No quartz, or no valid rate: a real outcome, reported as invalid
            // rather than as a number.
            out.push_back(s);
            continue;
        }

        s.accumulated_conc     = num / den;
        s.true_erate_g_cm2_yr  = true_sum / static_cast<double>(n_true);

        double predicted = 0.0;
        s.apparent_erate_g_cm2_yr = calculate_eff_erate_from_conc(
            s.accumulated_conc, nuclide, muon_scheme,
            prod_raster, topo_shield, snow_shield, self_shield,
            flow, outlet, predicted);

        if (s.true_erate_g_cm2_yr > 0.0 && s.apparent_erate_g_cm2_yr > 0.0)
        {
            s.ratio = s.apparent_erate_g_cm2_yr / s.true_erate_g_cm2_yr;
            s.valid = true;
        }
        out.push_back(s);
    }
    return out;
}


// -----------------------------------------------------------------------------
std::vector<NiemiAreaBin> bin_niemi_by_area(
    const std::vector<NiemiBasinSample>& samples,
    int                                  bands_per_decade,
    int                                  min_per_bin)
{
    std::vector<NiemiBasinSample> valid;
    for (const auto& s : samples) if (s.valid && s.area_km2 > 0.0) valid.push_back(s);
    if (valid.empty()) return {};

    double lo = 1e300, hi = -1e300;
    for (const auto& s : valid)
    {
        lo = std::min(lo, std::log10(s.area_km2));
        hi = std::max(hi, std::log10(s.area_km2));
    }
    const double band_w = 1.0 / static_cast<double>(std::max(1, bands_per_decade));
    const int n_bands = std::max(1,
        static_cast<int>(std::ceil((hi - lo) / band_w)) + 1);

    std::vector<std::vector<double>> ratios(static_cast<std::size_t>(n_bands));
    for (const auto& s : valid)
    {
        int b = static_cast<int>((std::log10(s.area_km2) - lo) / band_w);
        b = std::clamp(b, 0, n_bands - 1);
        ratios[static_cast<std::size_t>(b)].push_back(s.ratio);
    }

    std::vector<NiemiAreaBin> out;
    for (int b = 0; b < n_bands; ++b)
    {
        const auto& v = ratios[static_cast<std::size_t>(b)];
        if (static_cast<int>(v.size()) < min_per_bin) continue;

        NiemiAreaBin bin;
        bin.area_lo_km2 = std::pow(10.0, lo + b * band_w);
        bin.area_hi_km2 = std::pow(10.0, lo + (b + 1) * band_w);
        bin.n = static_cast<int>(v.size());
        bin.median = percentile(v, 0.50);
        bin.p25 = percentile(v, 0.25);
        bin.p75 = percentile(v, 0.75);
        bin.p05 = percentile(v, 0.05);
        bin.p95 = percentile(v, 0.95);
        bin.iqr = bin.p75 - bin.p25;

        int w10 = 0, w20 = 0;
        for (double x : v)
        {
            if (std::fabs(x - 1.0) <= 0.10) ++w10;
            if (std::fabs(x - 1.0) <= 0.20) ++w20;
        }
        bin.frac_within_10pc = static_cast<double>(w10) / static_cast<double>(v.size());
        bin.frac_within_20pc = static_cast<double>(w20) / static_cast<double>(v.size());
        out.push_back(bin);
    }
    return out;
}


// -----------------------------------------------------------------------------
double niemi_averaging_area_km2(const std::vector<NiemiAreaBin>& bins,
                                double tolerance, double fraction)
{
    // Walk DOWN from the largest band: the answer is the smallest area at and
    // above which the criterion holds everywhere, so a single noisy small bin
    // cannot claim convergence that later bands contradict.
    double answer = -1.0;
    for (auto it = bins.rbegin(); it != bins.rend(); ++it)
    {
        const double frac = (tolerance <= 0.10) ? it->frac_within_10pc
                                                : it->frac_within_20pc;
        if (frac >= fraction) answer = it->area_lo_km2;
        else break;
    }
    return answer;
}


// -----------------------------------------------------------------------------
void write_niemi_basins_csv(const std::string& path,
                            const std::vector<NiemiBasinSample>& samples)
{
    std::ofstream f(path);
    if (!f) return;
    f << "realisation,outlet_flat,n_pixels,area_km2,n_failed_pixels,"
         "fraction_failed,accumulated_conc_atoms_g,apparent_erate_g_cm2_yr,"
         "true_erate_g_cm2_yr,ratio,valid\n";
    f << std::setprecision(10);
    for (const auto& s : samples)
    {
        f << s.realisation << "," << s.outlet_flat << "," << s.n_pixels << ","
          << s.area_km2 << "," << s.n_failed_pixels << "," << s.fraction_failed
          << "," << s.accumulated_conc << "," << s.apparent_erate_g_cm2_yr
          << "," << s.true_erate_g_cm2_yr << "," << s.ratio << ","
          << (s.valid ? 1 : 0) << "\n";
    }
}


void write_niemi_area_bins_csv(const std::string& path,
                               const std::vector<NiemiAreaBin>& bins)
{
    std::ofstream f(path);
    if (!f) return;
    f << "area_lo_km2,area_hi_km2,n,median_ratio,p05,p25,p75,p95,iqr,"
         "frac_within_10pc,frac_within_20pc\n";
    f << std::setprecision(10);
    for (const auto& b : bins)
    {
        f << b.area_lo_km2 << "," << b.area_hi_km2 << "," << b.n << ","
          << b.median << "," << b.p05 << "," << b.p25 << "," << b.p75 << ","
          << b.p95 << "," << b.iqr << "," << b.frac_within_10pc << ","
          << b.frac_within_20pc << "\n";
    }
}

}  // namespace lsdtt
