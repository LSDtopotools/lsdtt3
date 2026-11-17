// test_chi_split_hierarchy.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
//
// -----------------------------------------------------------------------------
// The split hierarchy of the chi fit and the knickpoint measures built on it
// (docs/plan_docs/plan_chi_knickpoints.md, Part A), on synthetic one-row
// channels fitted by fit_chi_network.
//
// 1. Reach change with noise: chi = 0.02 i (i = 0..300), k_sn 100 up to
//    chi 4 (node 200) and 25 above, plus five 6 m steps over one node below
//    (nodes 30, 60, 90, 120, 150). Critical divergence 2 m. Expected: node 200
//    is the branch's first split (rank 0, depth 0) and has the largest
//    |split_delta_ksn_metres|; some step split has a larger raw
//    |split_delta_ksn| than node 200 (why the raw change is not the measure);
//    split records exist exactly at the split points.
// 2. Knickzone: k_sn 20 to chi 2 (node 100), 200 to chi 3 (node 150), 20 to
//    chi 6, plus 3 m one-node steps every 12 nodes. Threshold
//    median_plus_3mad. Expected: one knickzone, foot at chi 2, lip at chi 3,
//    k_sn below ~20, inside ~200, above ~20; both rows is_major; is_major is
//    set on exactly the within-branch rows whose |split_delta_ksn_normalised|
//    exceeds the abs_split_delta_ksn_normalised threshold of the method.
// 3. Outlet on a knickzone: k_sn 200 to chi 1 (node 50), 20 to chi 6, plus
//    6 m one-node steps every 12 nodes. Expected: one knickzone, foot_type
//    basin_outlet, foot at the outlet (chi 0, flat 0), lip within 0.1 of
//    chi 1 (the step at node 54 moves it to chi 1.08), k_sn
//    inside ~200, ksn_below and the foot split fields NaN.
// 4. Foot below the lip's threshold: k_sn 100 to chi 2, 200 to chi 3, 10 to
//    chi 6, 6 m one-node steps every 12 nodes, method otsu. The lip (chi 3)
//    is major; the foot (chi 2) is not (|split_delta_ksn_metres| 94.8, otsu
//    95.2) but passes the p95 (7.2). Expected: one knickzone, foot_type
//    split, foot chi 2, lip chi 3 (not a knickzone from the outlet).
// 5. An unknown threshold method throws.
//
// Usage: lsdtt-test-chi-split-hierarchy
// -----------------------------------------------------------------------------

#include <cmath>
#include <functional>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "../../src/cpu/lsdtt_chi_knickpoints_cpu.hpp"
#include "../../src/cpu/lsdtt_chi_network_analysis_cpu.hpp"
#include "../../src/cpu/lsdtt_raster_cpu.hpp"

namespace
{
int g_n = 0, g_fail = 0;
void check(const std::string& name, bool ok, const std::string& detail = "")
{
    ++g_n;
    if (!ok) ++g_fail;
    std::cout << (ok ? "PASS" : "FAIL") << "  " << name;
    if (!ok && !detail.empty()) std::cout << "   [" << detail << "]";
    std::cout << "\n";
}
std::string s(double v) { return std::to_string(v); }

using namespace lsdtt;

struct channel
{
    lsdtt_chi_basin_result basin;
    lsdtt_raster_cpu dem;
};

// A one-row channel: node i at column i, draining to i - 1; chi = dchi * i.
channel make_channel(int n, double dchi, const std::function<double(int)>& z_of)
{
    lsdtt_raster_meta m;
    m.n_rows = 1;
    m.n_cols = n;
    m.cell_size = 10.0;
    m.x_min = 500000.0;
    m.y_min = 5000000.0;
    m.set_crs_from_epsg(32632);
    std::vector<float> z(static_cast<std::size_t>(n));
    for (int i = 0; i < n; ++i) z[static_cast<std::size_t>(i)] = static_cast<float>(z_of(i));
    channel c{{}, lsdtt_raster_cpu(std::move(z), std::vector<std::uint8_t>(static_cast<std::size_t>(n), 1u), m)};
    c.basin.outlet_flat = 0;
    c.basin.outlet_chi = 0.0f;
    for (int i = 0; i < n; ++i)
    {
        c.basin.flats.push_back(i);
        c.basin.chi.push_back(static_cast<float>(dchi * i));
        c.basin.receiver_flats.push_back(i == 0 ? k_outlet_receiver : i - 1);
    }
    return c;
}

const lsdtt_chi_knickpoint_row* row_at(const lsdtt_chi_tributary_knickpoint_tables& t, std::int64_t flat)
{
    for (const auto& k : t.knickpoints)
        if (k.flat == flat && !k.at_junction) return &k;
    return nullptr;
}
} // namespace

int main()
{
    // ---- 1. reach change with noise ------------------------------------------
    {
        const double dchi = 0.02;
        auto z_of = [&](int i)
        {
            const double chi = dchi * i;
            double z = chi <= 4.0 ? 100.0 * chi : 400.0 + 25.0 * (chi - 4.0);
            for (const int st : {30, 60, 90, 120, 150})
                if (i >= st) z += 6.0;
            return z;
        };
        channel c = make_channel(301, dchi, z_of);
        lsdtt_chi_network_fit_options o;
        o.critical_divergence_m = 2.0f;
        const lsdtt_chi_network_fit_result r = fit_chi_network(c.basin, c.dem, o);
        bool records_match = r.split_records.size() == r.flats.size();
        for (std::size_t i = 0; records_match && i < r.flats.size(); ++i)
            records_match = (r.split_records[i].rank >= 0) == (r.is_split_point[i] != 0);
        check("split records exist exactly at split points", records_match);
        const auto& top = r.split_records[200];
        check("reach change at node 200 is the first split", top.rank == 0 && top.depth == 0,
              "rank " + s(top.rank) + " depth " + s(top.depth));
        check("first split chords: k_sn below > above", top.ksn_below > top.ksn_above,
              s(top.ksn_below) + " / " + s(top.ksn_above));

        lsdtt_chi_network_fit_collection_result fits;
        fits.basins.push_back(r);
        lsdtt_chi_collection_result nets;
        nets.basins.push_back(c.basin);
        const auto t = build_chi_tributaries_and_knickpoints(fits, nets, {}, c.dem.get_meta());
        const lsdtt_chi_knickpoint_row* k200 = row_at(t, 200);
        check("node 200 has a knickpoint row with split fields",
              k200 && k200->split_rank == 0 && std::isfinite(k200->split_delta_ksn_metres));
        if (k200)
        {
            double max_other_m = 0.0, max_other_raw = 0.0;
            for (const auto& k : t.knickpoints)
            {
                if (k.at_junction || k.flat == 200 || !std::isfinite(k.split_delta_ksn_metres)) continue;
                max_other_m = std::max(max_other_m, std::fabs(k.split_delta_ksn_metres));
                max_other_raw = std::max(max_other_raw, std::fabs(k.split_delta_ksn));
            }
            check("node 200 has the largest |split_delta_ksn_metres|",
                  std::fabs(k200->split_delta_ksn_metres) > max_other_m,
                  s(k200->split_delta_ksn_metres) + " vs " + s(max_other_m));
            check("a step split has a larger raw |split_delta_ksn| than node 200",
                  max_other_raw > std::fabs(k200->split_delta_ksn),
                  s(max_other_raw) + " vs " + s(k200->split_delta_ksn));
            const double w = std::sqrt(1.0 / (k200->split_dchi_below * k200->split_dchi_below) +
                                       1.0 / (k200->split_dchi_above * k200->split_dchi_above));
            check("split_delta_ksn_metres = delta / sqrt(1/db^2 + 1/da^2)",
                  std::fabs(k200->split_delta_ksn_metres - k200->split_delta_ksn / w) < 1e-6);
            check("normalised = metres / branch chi span (6)",
                  std::fabs(k200->split_delta_ksn_normalised - k200->split_delta_ksn_metres / 6.0) < 1e-4,
                  s(k200->split_delta_ksn_normalised));
        }
        bool has_rows = false;
        for (const auto& th : t.thresholds)
            has_rows = has_rows || (th.measure == "abs_split_delta_ksn_metres" && th.method == "p95");
        check("thresholds include abs_split_delta_ksn_metres", has_rows);
    }

    // ---- 2. knickzone --------------------------------------------------------
    {
        const double dchi = 0.02;
        auto z_of = [&](int i)
        {
            const double chi = dchi * i;
            double z = chi <= 2.0 ? 20.0 * chi
                     : chi <= 3.0 ? 40.0 + 200.0 * (chi - 2.0)
                                  : 240.0 + 20.0 * (chi - 3.0);
            for (int st = 6; st < 300; st += 12)
                if (i >= st) z += 3.0;
            return z;
        };
        channel c = make_channel(301, dchi, z_of);
        lsdtt_chi_network_fit_options o;
        o.critical_divergence_m = 2.0f;
        lsdtt_chi_network_fit_collection_result fits;
        fits.basins.push_back(fit_chi_network(c.basin, c.dem, o));
        lsdtt_chi_collection_result nets;
        nets.basins.push_back(c.basin);
        const auto t = build_chi_tributaries_and_knickpoints(fits, nets, {}, c.dem.get_meta(),
                                                             "median_plus_3mad");
        check("one knickzone", t.knickzones.size() == 1, s(static_cast<double>(t.knickzones.size())));
        if (!t.knickzones.empty())
        {
            const auto& z = t.knickzones[0];
            check("foot at chi 2, lip at chi 3",
                  std::fabs(z.foot_chi - 2.0) < 0.05 && std::fabs(z.lip_chi - 3.0) < 0.05,
                  s(z.foot_chi) + " / " + s(z.lip_chi));
            check("k_sn inside ~200", std::fabs(z.ksn_inside - 200.0) < 25.0, s(z.ksn_inside));
            check("k_sn below and above ~20 (+ the steps' 150 per chi)",
                  z.ksn_below < 100.0 && z.ksn_above < 100.0, s(z.ksn_below) + " / " + s(z.ksn_above));
            check("threshold method recorded", z.threshold_method == "median_plus_3mad");
            check("foot_type split", z.foot_type == "split", z.foot_type);
            const auto* foot = row_at(t, z.foot_flat);
            const auto* lip = row_at(t, z.lip_flat);
            check("foot and lip rows are major", foot && lip && foot->is_major && lip->is_major);
        }
        double thr = std::nan("");
        for (const auto& th : t.thresholds)
            if (th.measure == "abs_split_delta_ksn_normalised" && th.method == "median_plus_3mad")
                thr = th.threshold;
        bool major_matches = std::isfinite(thr);
        long n_major = 0;
        for (const auto& k : t.knickpoints)
        {
            const bool expect = !k.at_junction && std::isfinite(k.split_delta_ksn_normalised) &&
                                std::fabs(k.split_delta_ksn_normalised) > thr;
            major_matches = major_matches && (expect == (k.is_major != 0));
            n_major += k.is_major;
        }
        check("is_major = |split_delta_ksn_normalised| above the normalised threshold",
              major_matches && n_major >= 2, "threshold " + s(thr) + ", n_major " + s(n_major));
    }

    // ---- 3. outlet on a knickzone -------------------------------------------
    {
        const double dchi = 0.02;
        auto z_of = [&](int i)
        {
            const double chi = dchi * i;
            double z = chi <= 1.0 ? 200.0 * chi : 200.0 + 20.0 * (chi - 1.0);
            for (int st = 6; st < 300; st += 12)
                if (i >= st) z += 6.0;
            return z;
        };
        channel c = make_channel(301, dchi, z_of);
        lsdtt_chi_network_fit_options o;
        o.critical_divergence_m = 2.0f;
        lsdtt_chi_network_fit_collection_result fits;
        fits.basins.push_back(fit_chi_network(c.basin, c.dem, o));
        lsdtt_chi_collection_result nets;
        nets.basins.push_back(c.basin);
        const auto t = build_chi_tributaries_and_knickpoints(fits, nets, {}, c.dem.get_meta(),
                                                             "median_plus_3mad");
        check("outlet knickzone: one knickzone", t.knickzones.size() == 1,
              s(static_cast<double>(t.knickzones.size())));
        if (!t.knickzones.empty())
        {
            const auto& z = t.knickzones[0];
            check("outlet knickzone: foot_type basin_outlet", z.foot_type == "basin_outlet", z.foot_type);
            check("outlet knickzone: foot at the outlet, lip at chi 1",
                  z.foot_flat == 0 && std::fabs(z.foot_chi) < 1e-9 && std::fabs(z.lip_chi - 1.0) < 0.1,
                  s(static_cast<double>(z.foot_flat)) + " / " + s(z.foot_chi) + " / " + s(z.lip_chi));
            check("outlet knickzone: k_sn inside ~200", std::fabs(z.ksn_inside - 200.0) < 25.0, s(z.ksn_inside));
            check("outlet knickzone: ksn_below and foot split fields NaN",
                  std::isnan(z.ksn_below) && std::isnan(z.foot_split_delta_ksn_metres) &&
                      std::isnan(z.foot_split_delta_ksn_normalised));
        }
    }

    // ---- 4. foot below the lip's threshold --------------------------------------
    {
        const double dchi = 0.02;
        auto z_of = [&](int i)
        {
            const double chi = dchi * i;
            double z = chi <= 2.0 ? 100.0 * chi
                     : chi <= 3.0 ? 200.0 + 200.0 * (chi - 2.0)
                                  : 400.0 + 10.0 * (chi - 3.0);
            for (int st = 6; st < 300; st += 12)
                if (i >= st) z += 6.0;
            return z;
        };
        channel c = make_channel(301, dchi, z_of);
        lsdtt_chi_network_fit_options o;
        o.critical_divergence_m = 2.0f;
        lsdtt_chi_network_fit_collection_result fits;
        fits.basins.push_back(fit_chi_network(c.basin, c.dem, o));
        lsdtt_chi_collection_result nets;
        nets.basins.push_back(c.basin);
        const auto t = build_chi_tributaries_and_knickpoints(fits, nets, {}, c.dem.get_meta(), "otsu");
        check("sub-threshold foot: one knickzone", t.knickzones.size() == 1,
              s(static_cast<double>(t.knickzones.size())));
        if (!t.knickzones.empty())
        {
            const auto& z = t.knickzones[0];
            check("sub-threshold foot: foot_type split, foot chi 2, lip chi 3",
                  z.foot_type == "split" && std::fabs(z.foot_chi - 2.0) < 0.05 && std::fabs(z.lip_chi - 3.0) < 0.05,
                  z.foot_type + " " + s(z.foot_chi) + " / " + s(z.lip_chi));
            const auto* foot = row_at(t, z.foot_flat);
            const auto* lip = row_at(t, z.lip_flat);
            check("sub-threshold foot: lip major, foot not major but above the p95 of "
                  "|split_delta_ksn_metres|",
                  foot && lip && lip->is_major && !foot->is_major &&
                      std::fabs(foot->split_delta_ksn_metres) > z.foot_threshold_split_delta_ksn_metres,
                  foot ? s(foot->split_delta_ksn_metres) + " vs " + s(z.foot_threshold_split_delta_ksn_metres) : "no foot row");
        }
    }

    // ---- 5. bad method -------------------------------------------------------
    {
        bool threw = false;
        try
        {
            channel c = make_channel(10, 0.1, [](int i) { return 10.0 * i; });
            lsdtt_chi_network_fit_collection_result fits;
            fits.basins.push_back(fit_chi_network(c.basin, c.dem, {}));
            lsdtt_chi_collection_result nets;
            nets.basins.push_back(c.basin);
            (void)build_chi_tributaries_and_knickpoints(fits, nets, {}, c.dem.get_meta(), "p99");
        }
        catch (const std::invalid_argument&)
        {
            threw = true;
        }
        check("unknown threshold method throws", threw);
    }

    std::cout << "\n" << (g_n - g_fail) << " / " << g_n << " checks passed\n";
    return g_fail == 0 ? 0 : 1;
}
