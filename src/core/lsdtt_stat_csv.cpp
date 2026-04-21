// lsdtt_stat_csv.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html

#include "lsdtt_stat_csv.hpp"

#include <cmath>
#include <limits>
#include <string>
#include <vector>

namespace lsdtt {

// -----------------------------------------------------------------------------
void write_qq_csv(
    lsdtt_data_loader&    loader,
    const std::string&    out_path,
    const lsdtt_qq_data&  qq)
{
    if (qq.standard_normal_variate.empty()) return;

    std::vector<std::vector<std::string>> rows;
    rows.reserve(qq.standard_normal_variate.size());
    for (std::size_t i = 0; i < qq.standard_normal_variate.size(); ++i)
    {
        rows.push_back({
            std::to_string(qq.standard_normal_variate[i]),
            std::to_string(qq.data_quantile[i]),
            std::to_string(qq.normal_reference[i])
        });
    }
    loader.write_csv(out_path,
                     {"standard_normal_variate", "data_quantile", "normal_reference"},
                     rows);
}

// -----------------------------------------------------------------------------
void write_smoothed_cdf_csv(
    lsdtt_data_loader&                         loader,
    const std::string&                         out_path,
    const lsdtt_smoothed_cdf&                  scdf,
    const std::vector<lsdtt_distribution_fit>& fits)
{
    if (scdf.empty()) return;

    const lsdtt_distribution_fit* best = nullptr;
    for (const auto& f : fits)
    {
        if (!f.valid || !std::isfinite(f.ks_statistic)) continue;
        if (!best || f.ks_statistic < best->ks_statistic) best = &f;
    }

    std::vector<std::vector<std::string>> rows;
    rows.reserve(scdf.bin_centres.size());
    for (std::size_t i = 0; i < scdf.bin_centres.size(); ++i)
    {
        const double x   = scdf.bin_centres[i];
        const double fit = best ? eval_cdf(*best, x)
                                : std::numeric_limits<double>::quiet_NaN();
        rows.push_back({ std::to_string(x),
                         std::to_string(scdf.cdf[i]),
                         std::to_string(fit) });
    }

    const std::string fit_col = best
        ? std::string("fitted_cdf_") + distribution_family_name(best->family)
        : "fitted_cdf";

    loader.write_csv(out_path, {"bin_centre", "cdf", fit_col}, rows);
}

// -----------------------------------------------------------------------------
void write_dist_fits_csv(
    lsdtt_data_loader&                         loader,
    const std::string&                         out_path,
    const std::vector<lsdtt_distribution_fit>& fits)
{
    if (fits.empty()) return;

    std::vector<std::vector<std::string>> rows;
    rows.reserve(fits.size());
    for (const auto& f : fits)
    {
        const double ks_crit = (f.n_values > 0)
            ? 1.36 / std::sqrt(static_cast<double>(f.n_values))
            : std::numeric_limits<double>::quiet_NaN();

        rows.push_back({
            distribution_family_name(f.family),
            std::to_string(f.param1),
            std::to_string(f.param2),
            std::to_string(f.ks_statistic),
            std::to_string(f.max_divergence_x),
            std::to_string(f.n_values),
            std::to_string(ks_crit),
            f.valid ? "true" : "false"
        });
    }
    loader.write_csv(out_path,
        {"family","param1","param2","ks_statistic","max_divergence_x",
         "n_values","ks_critical_value_05","valid"},
        rows);
}

} // namespace lsdtt
