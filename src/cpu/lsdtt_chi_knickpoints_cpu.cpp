// lsdtt_chi_knickpoints_cpu.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
//
// -----------------------------------------------------------------------------
// See lsdtt_chi_knickpoints_cpu.hpp.
// -----------------------------------------------------------------------------

#include "lsdtt_chi_knickpoints_cpu.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <map>
#include <sstream>
#include <stdexcept>
#include <unordered_map>

namespace lsdtt
{

namespace
{

std::vector<double> finite_only(const std::vector<double>& v)
{
    std::vector<double> out;
    out.reserve(v.size());
    for (const double x : v)
        if (std::isfinite(x)) out.push_back(x);
    return out;
}

std::string num(double v)
{
    if (!std::isfinite(v)) return "";
    std::ostringstream o;
    o << std::setprecision(10) << v;
    return o.str();
}

std::ofstream open_csv(const std::string& filename)
{
    std::ofstream out(filename);
    if (!out) throw std::runtime_error("cannot open " + filename + " for writing");
    return out;
}

} // namespace

// -----------------------------------------------------------------------------
// threshold methods
// -----------------------------------------------------------------------------
double chi_percentile_linear(std::vector<double> v, double p)
{
    v = finite_only(v);
    if (v.empty()) return k_chi_nan;
    std::sort(v.begin(), v.end());
    const double pos = std::clamp(p, 0.0, 100.0) / 100.0 * static_cast<double>(v.size() - 1);
    const std::size_t lo = static_cast<std::size_t>(std::floor(pos));
    const std::size_t hi = std::min(lo + 1, v.size() - 1);
    return v[lo] + (v[hi] - v[lo]) * (pos - static_cast<double>(lo));
}

double chi_median_plus_3mad(const std::vector<double>& v_in)
{
    const std::vector<double> v = finite_only(v_in);
    if (v.size() < 2) return k_chi_nan;
    const double med = chi_percentile_linear(v, 50.0);
    std::vector<double> dev(v.size());
    for (std::size_t i = 0; i < v.size(); ++i) dev[i] = std::fabs(v[i] - med);
    const double mad = 1.4826 * chi_percentile_linear(dev, 50.0);
    return med + 3.0 * mad;
}

double chi_otsu_threshold(const std::vector<double>& v_in, int n_bins)
{
    const std::vector<double> v = finite_only(v_in);
    if (v.size() < 2 || n_bins < 2) return k_chi_nan;
    const auto [mn_it, mx_it] = std::minmax_element(v.begin(), v.end());
    const double mn = *mn_it, mx = *mx_it;
    if (!(mx > mn)) return k_chi_nan;
    const double w = (mx - mn) / n_bins;
    std::vector<double> hist(static_cast<std::size_t>(n_bins), 0.0);
    for (const double x : v)
    {
        int b = static_cast<int>((x - mn) / w);
        b = std::clamp(b, 0, n_bins - 1);
        hist[static_cast<std::size_t>(b)] += 1.0;
    }
    const double total = static_cast<double>(v.size());
    double sum_all = 0.0;
    for (int b = 0; b < n_bins; ++b) sum_all += (mn + (b + 0.5) * w) * hist[static_cast<std::size_t>(b)];
    double w0 = 0.0, sum0 = 0.0, best = -1.0;
    int best_k = 0;
    for (int k = 0; k < n_bins - 1; ++k)
    {
        w0 += hist[static_cast<std::size_t>(k)];
        sum0 += (mn + (k + 0.5) * w) * hist[static_cast<std::size_t>(k)];
        const double w1 = total - w0;
        if (w0 <= 0.0 || w1 <= 0.0) continue;
        const double m0 = sum0 / w0, m1 = (sum_all - sum0) / w1;
        const double between = w0 * w1 * (m0 - m1) * (m0 - m1);
        if (between > best) { best = between; best_k = k; }
    }
    // The upper edge of the lower class.
    return mn + (best_k + 1) * w;
}

std::vector<lsdtt_chi_threshold_row> chi_knickpoint_thresholds(
    const std::vector<lsdtt_chi_knickpoint_row>& knickpoints,
    const std::vector<double>&                   segment_ksn)
{
    std::vector<double> abs_delta, abs_log;
    for (const auto& k : knickpoints)
    {
        if (k.at_junction) continue;
        if (std::isfinite(k.delta_ksn)) abs_delta.push_back(std::fabs(k.delta_ksn));
        if (std::isfinite(k.log_ratio_ksn)) abs_log.push_back(std::fabs(k.log_ratio_ksn));
    }
    std::vector<lsdtt_chi_threshold_row> rows;
    auto add = [&](const std::string& measure, const std::string& method, double thr,
                   const std::vector<double>& values)
    {
        lsdtt_chi_threshold_row r;
        r.measure = measure;
        r.method = method;
        r.threshold = thr;
        r.n_total = static_cast<long>(values.size());
        if (std::isfinite(thr))
            for (const double x : values) if (x > thr) ++r.n_above;
        rows.push_back(r);
    };
    for (const auto& [measure, values] : {std::make_pair(std::string("abs_delta_ksn"), abs_delta),
                                          std::make_pair(std::string("abs_log_ratio_ksn"), abs_log)})
    {
        add(measure, "p90", chi_percentile_linear(values, 90.0), values);
        add(measure, "p95", chi_percentile_linear(values, 95.0), values);
        add(measure, "median_plus_3mad", chi_median_plus_3mad(values), values);
        add(measure, "otsu", chi_otsu_threshold(values), values);
    }
    const std::vector<double> ksn = finite_only(segment_ksn);
    double q3i = k_chi_nan;
    if (ksn.size() >= 2)
    {
        const double q1 = chi_percentile_linear(ksn, 25.0), q3 = chi_percentile_linear(ksn, 75.0);
        q3i = q3 + 3.0 * (q3 - q1);
    }
    add("ksn", "q3_plus_3iqr", q3i, ksn);
    return rows;
}

// -----------------------------------------------------------------------------
// tables
// -----------------------------------------------------------------------------
lsdtt_chi_tributary_knickpoint_tables build_chi_tributaries_and_knickpoints(
    const lsdtt_chi_network_fit_collection_result& fits,
    const lsdtt_chi_collection_result&             chi_networks,
    const std::vector<lsdtt_chi_basin_source_keys>& keys,
    const lsdtt_raster_meta&                       meta)
{
    lsdtt_chi_tributary_knickpoint_tables out;
    const int nc = meta.n_cols, nr = meta.n_rows;
    const double cs = meta.cell_size;
    auto x_of = [&](std::int64_t f) { return meta.x_min + (static_cast<double>(f % nc) + 0.5) * cs; };
    auto y_of = [&](std::int64_t f) { return meta.y_min + (nr - static_cast<double>(f / nc) - 0.5) * cs; };
    std::vector<double> all_ksn;

    for (std::size_t bi = 0; bi < fits.basins.size(); ++bi)
    {
        const auto& r = fits.basins[bi];
        const std::size_t n = r.flats.size();
        if (n == 0) continue;
        std::unordered_map<std::int64_t, std::size_t> fit_idx;
        fit_idx.reserve(n);
        for (std::size_t i = 0; i < n; ++i) fit_idx[r.flats[i]] = i;
        // per-node area, flow distance and key, by flat (keys follow chi_networks order)
        std::unordered_map<std::int64_t, std::size_t> key_idx;
        const bool have_keys = bi < keys.size() && bi < chi_networks.basins.size();
        if (have_keys)
        {
            const auto& fl = chi_networks.basins[bi].flats;
            key_idx.reserve(fl.size());
            for (std::size_t i = 0; i < fl.size(); ++i) key_idx[fl[i]] = i;
        }
        auto from_keys = [&](std::int64_t f, int what) -> double
        {
            if (!have_keys) return k_chi_nan;
            const auto it = key_idx.find(f);
            if (it == key_idx.end()) return k_chi_nan;
            const auto& k = keys[bi];
            if (what == 0) return it->second < k.drainage_area_m2.size() ? k.drainage_area_m2[it->second] : k_chi_nan;
            if (what == 1) return it->second < k.flow_distance_m.size() ? k.flow_distance_m[it->second] : k_chi_nan;
            return it->second < k.source_key.size() ? static_cast<double>(k.source_key[it->second]) : -1.0;
        };
        const std::int32_t basin_key = have_keys ? keys[bi].basin_key : -1;
        double chi_max = 0.0;
        for (const float c : r.chi) chi_max = std::max(chi_max, static_cast<double>(c));

        // ---- branches -----------------------------------------------------------
        std::map<std::int32_t, std::vector<std::size_t>> branch_nodes;
        for (std::size_t i = 0; i < n; ++i) branch_nodes[r.branch_id_of_node[i]].push_back(i);
        std::map<std::int32_t, lsdtt_chi_tributary_row> trib;
        std::map<std::int32_t, std::size_t> head_of;
        for (const auto& [b, nodes] : branch_nodes)
        {
            lsdtt_chi_tributary_row t;
            t.basin_index = static_cast<std::int32_t>(bi);
            t.basin_key = basin_key;
            t.branch_id = b;
            std::size_t head = nodes.front(), low = nodes.front();
            for (const std::size_t i : nodes)
            {
                if (r.chi[i] > r.chi[head]) head = i;
                if (r.chi[i] < r.chi[low]) low = i;
            }
            head_of[b] = head;
            t.head_flat = r.flats[head];
            t.head_x = x_of(t.head_flat);
            t.head_y = y_of(t.head_flat);
            t.source_key = static_cast<std::int32_t>(from_keys(t.head_flat, 2));
            const std::int64_t rec = r.receiver_flats[low];
            const auto jt = fit_idx.find(rec);
            double chi_down = r.chi[low], fd_down = from_keys(r.flats[low], 1);
            if (b != 0 && jt != fit_idx.end())
            {
                t.junction_flat = rec;
                t.junction_x = x_of(rec);
                t.junction_y = y_of(rec);
                t.parent_branch_id = r.branch_id_of_node[jt->second];
                chi_down = r.chi[jt->second];
                fd_down = from_keys(rec, 1);
                t.junction_drainage_area_m2 = from_keys(r.flats[low], 0);
            }
            else
            {
                t.junction_drainage_area_m2 = from_keys(r.flats[low], 0);
            }
            t.chi_span = static_cast<double>(r.chi[head]) - chi_down;
            t.length_m = from_keys(t.head_flat, 1) - fd_down;
            trib[b] = t;
        }
        const auto oit = fit_idx.find(r.outlet_flat);
        const double basin_area = from_keys(oit != fit_idx.end() ? r.outlet_flat : r.flats[0], 0);
        // segments per branch
        std::map<std::int32_t, std::vector<const lsdtt_chi_network_fit_segment*>> segs_of;
        std::unordered_map<std::int32_t, const lsdtt_chi_network_fit_segment*> seg_by_id;
        for (const auto& s : r.segments)
        {
            segs_of[s.branch_id].push_back(&s);
            seg_by_id[s.segment_id] = &s;
            all_ksn.push_back(s.k_sn);
        }
        for (auto& [b, t] : trib)
        {
            t.basin_drainage_area_m2 = basin_area;
            t.area_fraction = (std::isfinite(basin_area) && basin_area > 0.0)
                                  ? t.junction_drainage_area_m2 / basin_area : k_chi_nan;
            const auto sit = segs_of.find(b);
            if (sit == segs_of.end()) continue;
            std::vector<double> ks;
            for (const auto* s : sit->second)
            {
                ks.push_back(s->k_sn);
                if (t.ksn_max_segment_id < 0 || s->k_sn > t.ksn_max)
                {
                    t.ksn_max = s->k_sn;
                    t.ksn_max_segment_id = s->segment_id;
                    t.ksn_max_chi = 0.5 * (static_cast<double>(s->downstream_chi) + s->upstream_chi);
                }
            }
            t.n_segments = static_cast<std::int32_t>(ks.size());
            t.ksn_min = *std::min_element(ks.begin(), ks.end());
            t.ksn_median = chi_percentile_linear(ks, 50.0);
        }
        // order and descendants
        for (auto& [b, t] : trib)
        {
            std::int32_t p = t.parent_branch_id, depth = 0;
            for (int guard = 0; p >= 0 && guard < static_cast<int>(trib.size()) + 1; ++guard)
            {
                ++depth;
                ++trib[p].n_descendants;
                p = trib[p].parent_branch_id;
            }
            t.branch_order = depth;
        }
        // area rank among tributaries
        std::vector<std::int32_t> tribs;
        for (const auto& [b, t] : trib) if (b != 0) tribs.push_back(b);
        std::stable_sort(tribs.begin(), tribs.end(), [&](std::int32_t a, std::int32_t c)
        {
            const double xa = trib[a].junction_drainage_area_m2, xc = trib[c].junction_drainage_area_m2;
            return (std::isfinite(xa) ? xa : -1.0) > (std::isfinite(xc) ? xc : -1.0);
        });
        for (std::size_t q = 0; q < tribs.size(); ++q) trib[tribs[q]].area_rank_in_basin = static_cast<std::int32_t>(q + 1);
        for (const auto& [b, t] : trib) out.tributaries.push_back(t);

        // ---- knickpoints ----------------------------------------------------------
        // node -> segment covering it: from the segment's upstream end down to,
        // but not including, its downstream end.
        std::unordered_map<std::int64_t, const lsdtt_chi_network_fit_segment*> seg_of_node;
        for (const auto& s : r.segments)
        {
            std::int64_t cur = s.upstream_flat;
            for (int guard = 0; guard <= s.n_nodes + 1 && cur != s.downstream_flat; ++guard)
            {
                seg_of_node.emplace(cur, &s);
                const auto it = fit_idx.find(cur);
                if (it == fit_idx.end()) break;
                cur = r.receiver_flats[it->second];
            }
        }
        for (const auto& s : r.segments)
        {
            const lsdtt_chi_network_fit_segment* p = nullptr;
            bool junction = false;
            if (s.parent_segment_id >= 0)
            {
                const auto it = seg_by_id.find(s.parent_segment_id);
                if (it != seg_by_id.end()) p = it->second;
            }
            else if (s.branch_id != 0)
            {
                const auto it = seg_of_node.find(s.downstream_flat);
                if (it != seg_of_node.end() && it->second->branch_id != s.branch_id)
                {
                    p = it->second;
                    junction = true;
                }
            }
            if (!p) continue;
            lsdtt_chi_knickpoint_row k;
            k.basin_index = static_cast<std::int32_t>(bi);
            k.basin_key = basin_key;
            k.branch_id = s.branch_id;
            k.source_key = trib.count(s.branch_id) ? trib[s.branch_id].source_key : -1;
            k.upstream_segment_id = s.segment_id;
            k.downstream_segment_id = p->segment_id;
            k.flat = s.downstream_flat;
            k.x = x_of(k.flat);
            k.y = y_of(k.flat);
            const auto nit = fit_idx.find(k.flat);
            if (nit != fit_idx.end())
            {
                k.chi = r.chi[nit->second];
                k.elevation = r.observed_elevation[nit->second];
            }
            k.flow_distance_m = from_keys(k.flat, 1);
            k.drainage_area_m2 = from_keys(k.flat, 0);
            k.ksn_upstream = s.k_sn;
            k.ksn_downstream = p->k_sn;
            k.delta_ksn = k.ksn_upstream - k.ksn_downstream;
            if (k.ksn_upstream > 0.0 && k.ksn_downstream > 0.0)
                k.log_ratio_ksn = std::log10(k.ksn_upstream / k.ksn_downstream);
            k.at_junction = junction ? 1 : 0;
            const auto hit = head_of.find(s.branch_id);
            if (hit != head_of.end() && std::isfinite(k.chi))
                k.chi_from_branch_head = static_cast<double>(r.chi[hit->second]) - k.chi;
            if (chi_max > 0.0 && std::isfinite(k.chi)) k.chi_normalised = k.chi / chi_max;
            out.knickpoints.push_back(k);
        }
    }
    out.thresholds = chi_knickpoint_thresholds(out.knickpoints, all_ksn);
    return out;
}

// -----------------------------------------------------------------------------
// writers
// -----------------------------------------------------------------------------
void write_chi_tributaries_csv(const std::string& filename,
                               const std::vector<lsdtt_chi_tributary_row>& rows)
{
    std::ofstream out = open_csv(filename);
    out << "basin_index,basin_key,source_key,branch_id,parent_branch_id,branch_order,head_flat,"
           "junction_flat,head_x,head_y,junction_x,junction_y,junction_drainage_area_m2,"
           "basin_drainage_area_m2,area_fraction,area_rank_in_basin,length_m,chi_span,n_segments,"
           "ksn_min,ksn_median,ksn_max,ksn_max_segment_id,ksn_max_chi,n_descendants\n";
    for (const auto& t : rows)
        out << t.basin_index << ',' << t.basin_key << ',' << t.source_key << ',' << t.branch_id << ','
            << t.parent_branch_id << ',' << t.branch_order << ',' << t.head_flat << ',' << t.junction_flat << ','
            << num(t.head_x) << ',' << num(t.head_y) << ',' << num(t.junction_x) << ',' << num(t.junction_y) << ','
            << num(t.junction_drainage_area_m2) << ',' << num(t.basin_drainage_area_m2) << ','
            << num(t.area_fraction) << ',' << t.area_rank_in_basin << ',' << num(t.length_m) << ','
            << num(t.chi_span) << ',' << t.n_segments << ',' << num(t.ksn_min) << ',' << num(t.ksn_median) << ','
            << num(t.ksn_max) << ',' << t.ksn_max_segment_id << ',' << num(t.ksn_max_chi) << ','
            << t.n_descendants << '\n';
}

void write_chi_knickpoints_csv(const std::string& filename,
                               const std::vector<lsdtt_chi_knickpoint_row>& rows)
{
    std::ofstream out = open_csv(filename);
    out << "basin_index,basin_key,source_key,branch_id,upstream_segment_id,downstream_segment_id,flat,"
           "x,y,chi,elevation,flow_distance_m,drainage_area_m2,ksn_upstream,ksn_downstream,delta_ksn,"
           "log_ratio_ksn,at_junction,chi_from_branch_head,chi_normalised\n";
    for (const auto& k : rows)
        out << k.basin_index << ',' << k.basin_key << ',' << k.source_key << ',' << k.branch_id << ','
            << k.upstream_segment_id << ',' << k.downstream_segment_id << ',' << k.flat << ','
            << num(k.x) << ',' << num(k.y) << ',' << num(k.chi) << ',' << num(k.elevation) << ','
            << num(k.flow_distance_m) << ',' << num(k.drainage_area_m2) << ',' << num(k.ksn_upstream) << ','
            << num(k.ksn_downstream) << ',' << num(k.delta_ksn) << ',' << num(k.log_ratio_ksn) << ','
            << static_cast<int>(k.at_junction) << ',' << num(k.chi_from_branch_head) << ','
            << num(k.chi_normalised) << '\n';
}

void write_chi_knickpoint_thresholds_csv(const std::string& filename,
                                         const std::vector<lsdtt_chi_threshold_row>& rows)
{
    std::ofstream out = open_csv(filename);
    out << "measure,method,threshold,n_above,n_total\n";
    for (const auto& r : rows)
        out << r.measure << ',' << r.method << ',' << num(r.threshold) << ',' << r.n_above << ','
            << r.n_total << '\n';
}

} // namespace lsdtt
