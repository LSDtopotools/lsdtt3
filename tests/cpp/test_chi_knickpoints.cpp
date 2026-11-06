// test_chi_knickpoints.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
//
// -----------------------------------------------------------------------------
// build_chi_tributaries_and_knickpoints and the threshold methods, on a
// hand-built fit (no DEM).
//
// Grid 10 columns, 10 m cells. Main stem (branch 0): flats 0..6 along row 0,
// outlet 0, chi = column. Segments: s0 (0-2, k_sn 50), s1 (2-4, 100, parent
// s0), s2 (4-6, 25, parent s1). Tributary (branch 1) joins at flat 3 (inside
// s1): flats 13, 23, 33 (column 3, rows 1-3), chi 3.5, 4, 4.5; segments
// s3 (3-23, k_sn 80, no parent: a branch's first segment) and s4 (23-33,
// 160, parent s3). Drainage area: 1e6 m2 at the outlet, 2e5 at flat 13.
//
// Expected:
//   knickpoints: s1/s0 at flat 2 (delta +50, log10 2), s2/s1 at flat 4
//   (delta -75), s4/s3 at flat 23 (delta +80), and s3 at the junction (flat 3,
//   parent-branch segment s1 covers it, delta -20, at_junction 1);
//   tributary: parent 0, order 1, junction flat 3, area 2e5, fraction 0.2,
//   rank 1, source_key of its head (flat 33), n_segments 2, k_sn max 160;
//   main stem: n_descendants 1, fraction 1.
// Thresholds: percentiles, median + 3 MAD and Otsu on constructed values.
//
// Usage: lsdtt-test-chi-knickpoints
// -----------------------------------------------------------------------------

#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include "../../src/cpu/lsdtt_chi_knickpoints_cpu.hpp"

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
bool near(double a, double b, double tol = 1e-9) { return std::fabs(a - b) <= tol; }
std::string s(double v) { return std::to_string(v); }
} // namespace

int main()
{
    using namespace lsdtt;
    lsdtt_raster_meta meta;
    meta.n_rows = 10;
    meta.n_cols = 10;
    meta.x_min = 0.0;
    meta.y_min = 0.0;
    meta.cell_size = 10.0;

    // nodes: main 0..6, tributary 13, 23, 33
    const std::vector<std::int64_t> flats = {0, 1, 2, 3, 4, 5, 6, 13, 23, 33};
    const std::vector<std::int64_t> recv = {k_outlet_receiver, 0, 1, 2, 3, 4, 5, 3, 13, 23};
    const std::vector<float> chi = {0, 1, 2, 3, 4, 5, 6, 3.5f, 4, 4.5f};
    const std::vector<std::int32_t> branch = {0, 0, 0, 0, 0, 0, 0, 1, 1, 1};

    lsdtt_chi_network_fit_result r;
    r.outlet_flat = 0;
    r.flats = flats;
    r.receiver_flats = recv;
    r.chi = chi;
    r.observed_elevation = {100, 101, 102, 103, 104, 105, 106, 110, 112, 114};
    r.branch_id_of_node = branch;
    r.is_source = {0, 0, 0, 0, 0, 0, 1, 0, 0, 1};
    auto seg = [](int id, int parent, int br, std::int64_t dn, std::int64_t up, float cd, float cu, float k, int nn)
    {
        lsdtt_chi_network_fit_segment g;
        g.segment_id = id;
        g.parent_segment_id = parent;
        g.branch_id = br;
        g.downstream_flat = dn;
        g.upstream_flat = up;
        g.downstream_chi = cd;
        g.upstream_chi = cu;
        g.k_sn = k;
        g.n_nodes = nn;
        return g;
    };
    r.segments = {seg(0, -1, 0, 0, 2, 0, 2, 50, 3), seg(1, 0, 0, 2, 4, 2, 4, 100, 3),
                  seg(2, 1, 0, 4, 6, 4, 6, 25, 3), seg(3, -1, 1, 3, 23, 3, 4, 80, 3),
                  seg(4, 3, 1, 23, 33, 4, 4.5f, 160, 2)};
    lsdtt_chi_network_fit_collection_result fits;
    fits.basins.push_back(r);

    lsdtt_chi_collection_result chi_net;
    lsdtt_chi_basin_result b;
    b.outlet_flat = 0;
    b.flats = flats;
    b.receiver_flats = recv;
    b.chi = chi;
    chi_net.basins.push_back(b);

    std::vector<lsdtt_chi_basin_source_keys> keys(1);
    keys[0].basin_key = 0;
    keys[0].source_key = {0, 0, 0, 0, 0, 0, 0, 1, 1, 1};
    keys[0].is_source = {0, 0, 0, 0, 0, 0, 1, 0, 0, 1};
    keys[0].flow_distance_m = {0, 10, 20, 30, 40, 50, 60, 40, 50, 60};
    keys[0].drainage_area_m2 = {1e6, 9e5, 8e5, 7e5, 3e5, 2e5, 1e5, 2e5, 1e5, 5e4};

    const auto t = build_chi_tributaries_and_knickpoints(fits, chi_net, keys, meta);

    // ---- knickpoints ---------------------------------------------------------
    check("four knickpoint rows (3 within branches, 1 at the junction)", t.knickpoints.size() == 4,
          s(static_cast<double>(t.knickpoints.size())));
    auto find = [&](int up) -> const lsdtt_chi_knickpoint_row*
    {
        for (const auto& k : t.knickpoints) if (k.upstream_segment_id == up) return &k;
        return nullptr;
    };
    const auto* k1 = find(1);
    const auto* k2 = find(2);
    const auto* k3 = find(3);
    const auto* k4 = find(4);
    check("s1/s0 at flat 2: delta +50, log10 2, within branch",
          k1 && k1->flat == 2 && k1->downstream_segment_id == 0 && near(k1->delta_ksn, 50.0) &&
              near(k1->log_ratio_ksn, std::log10(2.0)) && k1->at_junction == 0);
    check("s2/s1 at flat 4: delta -75", k2 && k2->flat == 4 && near(k2->delta_ksn, -75.0));
    check("s4/s3 at flat 23: delta +80, chi 4, flow distance 50, area 1e5",
          k4 && k4->flat == 23 && near(k4->delta_ksn, 80.0) && near(k4->chi, 4.0) &&
              near(k4->flow_distance_m, 50.0) && near(k4->drainage_area_m2, 1e5));
    check("s3 at the junction: flat 3, downstream segment s1, delta -20, at_junction",
          k3 && k3->flat == 3 && k3->downstream_segment_id == 1 && near(k3->delta_ksn, -20.0) &&
              k3->at_junction == 1);
    check("knickpoint x / y are cell centres (flat 23: col 3, row 2)",
          k4 && near(k4->x, 35.0) && near(k4->y, 75.0));
    check("chi_from_branch_head and chi_normalised (flat 23: head chi 4.5, basin max 6)",
          k4 && near(k4->chi_from_branch_head, 0.5) && near(k4->chi_normalised, 4.0 / 6.0));
    check("knickpoint source_key is its branch's", k4 && k4->source_key == 1 && k1 && k1->source_key == 0);

    // ---- tributaries -----------------------------------------------------------
    check("two tributary rows", t.tributaries.size() == 2);
    const lsdtt_chi_tributary_row* m = nullptr;
    const lsdtt_chi_tributary_row* tr = nullptr;
    for (const auto& row : t.tributaries) (row.branch_id == 0 ? m : tr) = &row;
    check("tributary: parent 0, order 1, junction flat 3, head flat 33",
          tr && tr->parent_branch_id == 0 && tr->branch_order == 1 && tr->junction_flat == 3 && tr->head_flat == 33);
    check("tributary: area 2e5 before the junction, fraction 0.2, rank 1",
          tr && near(tr->junction_drainage_area_m2, 2e5) && near(tr->area_fraction, 0.2) && tr->area_rank_in_basin == 1,
          tr ? s(tr->junction_drainage_area_m2) + " " + s(tr->area_fraction) : "");
    check("tributary: source_key 1, length 30 m (60 - 30), chi span 1.5",
          tr && tr->source_key == 1 && near(tr->length_m, 30.0) && near(tr->chi_span, 1.5),
          tr ? s(tr->length_m) + " " + s(tr->chi_span) : "");
    check("tributary: 2 segments, k_sn max 160 (segment 4, mid chi 4.25), median 120",
          tr && tr->n_segments == 2 && near(tr->ksn_max, 160.0) && tr->ksn_max_segment_id == 4 &&
              near(tr->ksn_max_chi, 4.25) && near(tr->ksn_median, 120.0));
    check("main stem: order 0, rank 0, 1 descendant, fraction 1, junction -1",
          m && m->branch_order == 0 && m->area_rank_in_basin == 0 && m->n_descendants == 1 &&
              near(m->area_fraction, 1.0) && m->junction_flat == -1);

    // ---- thresholds ----------------------------------------------------------
    std::vector<double> v;
    for (int i = 0; i < 90; ++i) v.push_back(0.05 + 0.001 * (i % 10));
    for (int i = 0; i < 10; ++i) v.push_back(0.5 + 0.01 * i);
    const double otsu = chi_otsu_threshold(v);
    check("Otsu splits 90 small from 10 large values", otsu > 0.06 && otsu <= 0.5, s(otsu));
    check("percentile: p50 of 1..5 is 3, p90 of 1..5 is 4.6",
          near(chi_percentile_linear({1, 2, 3, 4, 5}, 50.0), 3.0) &&
              near(chi_percentile_linear({1, 2, 3, 4, 5}, 90.0), 4.6));
    // 1,2,3,4,100: median 3, |dev| 2,1,0,1,97 -> median 1 -> 3 + 3 * 1.4826
    check("median + 3 MAD (scaled)", near(chi_median_plus_3mad({1, 2, 3, 4, 100}), 3.0 + 3.0 * 1.4826));
    check("thresholds NaN with fewer than two values",
          std::isnan(chi_median_plus_3mad({1.0})) && std::isnan(chi_otsu_threshold({1.0, 1.0})));
    bool junction_left_out = false, has_ksn = false;
    for (const auto& th : t.thresholds)
    {
        if (th.measure == "abs_delta_ksn" && th.method == "p90") junction_left_out = th.n_total == 3;
        if (th.measure == "ksn" && th.method == "q3_plus_3iqr") has_ksn = th.n_total == 5;
    }
    check("thresholds use the 3 within-branch rows (junction left out)", junction_left_out);
    check("k_sn outlier threshold over all 5 segments", has_ksn);

    std::cout << "\n" << (g_n - g_fail) << " / " << g_n << " checks passed\n";
    return g_fail == 0 ? 0 : 1;
}
