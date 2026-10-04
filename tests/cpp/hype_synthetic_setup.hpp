// hype_synthetic_setup.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// Synthetic HYPE set-ups and forcing, built in memory. Shared, header-only, by
// performance_hype_model (test_hype_performance.cpp) and
// unit_hype_thread_invariance (test_hype_thread_invariance.cpp). Nothing is
// read from or written to disk.
//
// SET-UP (make_hype_synthetic_setup)
//   * n_subbasins sub-basins in a binary tree: sub-basin index s (subid s+1)
//     drains to index (s-1)/2; index 0 is the single domain outlet. The model
//     derives its own upstream-first order (hype_subbasin_downstream_order).
//   * n_classes classes cycling through 3 land uses, 4 soils and 1-3 layers.
//     Every sub-basin has every class, so n_units = n_subbasins x n_classes.
//   * Per-sub-basin area, elevation, slope, river lengths, class fractions and
//     class elevation offsets come from a fixed integer hash (splitmix64), so
//     the set-up is identical on every platform and compiler.
//   * Parameters are chosen so every kept process is active: snow (lapse rate,
//     sub-zero winters), PET/AET, diversion to surface runoff and macropores,
//     percolation, soil runoff with the slope term, saturation excess (a thin
//     soil with little drainable porosity), and local + main river routing
//     with translation times of 0-2 whole steps plus a fractional part.
//
// FORCING (make_hype_synthetic_forcing)
//   The formulas of tests/scripts/generate_hype_oracle.py (forcing_of):
//   precipitation is default_all's pattern repeated, rounded to 3 decimals;
//   temperature is mean + amp x sin(2 pi (step - phase) / period), rounded to
//   2 decimals, with a 365-step period (coldest around step 15). Day of year
//   is (step mod 365) + 1.
// -----------------------------------------------------------------------------

#pragma once

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <vector>

#include "../../src/cpu/lsdtt_hype_parameters.hpp"

namespace hype_synthetic
{

// Deterministic value in [0, 1) from two integer keys (splitmix64 finaliser).
inline double hash01(std::uint64_t a, std::uint64_t b)
{
    std::uint64_t z = a * 0x9E3779B97F4A7C15ULL + b + 0x632BE59BD9B4E019ULL;
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    z = z ^ (z >> 31);
    return static_cast<double>(z >> 11) * (1.0 / 9007199254740992.0);   // 53 bits
}

// Keys for hash01's second argument, one per derived quantity.
enum : std::uint64_t
{
    k_key_area = 1, k_key_elev, k_key_slope, k_key_locriv, k_key_riv, k_key_frac = 100, k_key_dh = 200
};

inline lsdtt::hype_setup make_hype_synthetic_setup(std::size_t n_subbasins, int n_classes)
{
    lsdtt::hype_setup s;

    // General parameters: default_all's (tests/expected/hype/cases/default_all.json),
    // with a slower river so translation spans whole steps.
    s.general.ttpd     = 0.5;
    s.general.ttpi     = 1.5;
    s.general.tcalt    = 0.6;
    s.general.lp       = 0.75;
    s.general.epotdist = 4.0;
    s.general.fsceff   = 0.0;
    s.general.whcsnow  = 0.08;
    s.general.rrcs3    = 0.002;
    s.general.rivvel   = 0.1;
    s.general.damp     = 0.4;
    s.general.cevpam   = 0.3;
    s.general.cevpph   = 20.0;

    const double ttmp[3]   = {0.0, 0.5, -0.5};
    const double cmlt[3]   = {3.5, 2.5, 3.0};
    const double cevp[3]   = {0.20, 0.25, 0.15};
    const double srrcs[3]  = {0.15, 0.10, 0.20};
    const double pcluse[3] = {0.05, 0.15, 0.10};
    for (int l = 0; l < 3; ++l)
    {
        lsdtt::hype_landuse_parameters p;
        p.code = l + 1;
        p.ttmp = ttmp[l];
        p.cmlt = cmlt[l];
        p.cevp = cevp[l];
        p.srrcs = srrcs[l];
        p.pcluse = pcluse[l];
        s.landuse.push_back(p);
    }

    // Soil 3 is thin-pored and slow-draining, so saturation excess occurs.
    const double wcwp[4]     = {0.08, 0.15, 0.05, 0.10};
    const double wcfc[4]     = {0.22, 0.25, 0.15, 0.20};
    const double wcep[4]     = {0.14, 0.08, 0.03, 0.10};
    const double rrcs1[4]    = {0.25, 0.15, 0.30, 0.20};
    const double rrcs2[4]    = {0.02, 0.01, 0.03, 0.015};
    const double mperc1[4]   = {20.0, 8.0, 4.0, 12.0};
    const double mperc2[4]   = {8.0, 3.0, 2.0, 5.0};
    const double macrate[4]  = {0.15, 0.05, 0.20, 0.10};
    const double srrate[4]   = {0.10, 0.20, 0.30, 0.15};
    const double mactrinf[4] = {12.0, 6.0, 4.0, 8.0};
    const double mactrsm[4]  = {0.60, 0.70, 0.50, 0.65};
    for (int k = 0; k < 4; ++k)
    {
        lsdtt::hype_soil_parameters p;
        p.code = k + 1;
        p.wcwp = wcwp[k];
        p.wcfc = wcfc[k];
        p.wcep = wcep[k];
        p.rrcs1 = rrcs1[k];
        p.rrcs2 = rrcs2[k];
        p.mperc1 = mperc1[k];
        p.mperc2 = mperc2[k];
        p.macrate = macrate[k];
        p.srrate = srrate[k];
        p.mactrinf = mactrinf[k];
        p.mactrsm = mactrsm[k];
        s.soils.push_back(p);
    }

    for (int c = 0; c < n_classes; ++c)
    {
        lsdtt::hype_class_definition cd;
        cd.class_id = c + 1;
        cd.landuse  = c % 3 + 1;
        cd.soil     = c % 4 + 1;
        cd.n_layers = (c / 2) % 3 + 1;
        const double l1 = 0.10 + 0.05 * (c % 3);
        const double l2 = l1 + 0.40 + 0.10 * (c % 2);
        const double l3 = l2 + 0.80;
        cd.layer_bottom_m = {l1, cd.n_layers >= 2 ? l2 : 0.0, cd.n_layers >= 3 ? l3 : 0.0};
        cd.stream_depth_m = 0.25 * (1 + c % 4);
        s.classes.push_back(cd);
    }

    for (std::size_t i = 0; i < n_subbasins; ++i)
    {
        lsdtt::hype_subbasin sb;
        sb.subid    = static_cast<int>(i) + 1;
        sb.maindown = (i == 0) ? 0 : static_cast<int>((i - 1) / 2) + 1;
        sb.area_m2            = 5.0e6 + 45.0e6 * hash01(i, k_key_area);
        sb.elev_mean_m        = 100.0 + 700.0 * hash01(i, k_key_elev);
        sb.slope_mean_percent = 1.0 + 24.0 * hash01(i, k_key_slope);
        sb.loc_rivlen_m       = 1000.0 + 7000.0 * hash01(i, k_key_locriv);
        sb.rivlen_m           = 3000.0 + 27000.0 * hash01(i, k_key_riv);
        s.subbasins.push_back(sb);

        std::vector<double> w(static_cast<std::size_t>(n_classes));
        double total = 0.0;
        for (int c = 0; c < n_classes; ++c)
        {
            w[static_cast<std::size_t>(c)] = 0.2 + hash01(i, k_key_frac + static_cast<std::uint64_t>(c));
            total += w[static_cast<std::size_t>(c)];
        }
        for (int c = 0; c < n_classes; ++c)
        {
            lsdtt::hype_subbasin_class sc;
            sc.subid    = sb.subid;
            sc.class_id = c + 1;
            sc.fraction = w[static_cast<std::size_t>(c)] / total;
            sc.dh_m     = -100.0 + 200.0 * hash01(i, k_key_dh + static_cast<std::uint64_t>(c));
            s.subbasin_classes.push_back(sc);
        }
    }
    return s;
}

struct hype_synthetic_forcing
{
    std::vector<double> prec_mm;
    std::vector<double> temp_c;
    std::vector<int>    doy;
};

inline hype_synthetic_forcing make_hype_synthetic_forcing(int n_steps)
{
    // default_all's precip_pattern_mm.
    static const double pattern[] = {0, 4, 11, 0, 0, 2, 0, 18, 6, 0, 0, 0, 1, 0, 27, 3,
                                     0, 0, 0, 0, 9, 0, 0, 14, 0, 0, 0, 0, 0, 35, 2};
    const std::size_t n_pat = sizeof(pattern) / sizeof(pattern[0]);
    const double pi = 3.14159265358979323846;
    const double temp_mean = 5.0, temp_amp = 9.0, period = 365.0, phase = 106.0;

    hype_synthetic_forcing f;
    for (int i = 0; i < n_steps; ++i)
    {
        f.prec_mm.push_back(std::round(pattern[static_cast<std::size_t>(i) % n_pat] * 1000.0) / 1000.0);
        const double t = temp_mean + temp_amp * std::sin(2.0 * pi * (i - phase) / period);
        f.temp_c.push_back(std::round(t * 100.0) / 100.0);
        f.doy.push_back(i % 365 + 1);
    }
    return f;
}

} // namespace hype_synthetic
