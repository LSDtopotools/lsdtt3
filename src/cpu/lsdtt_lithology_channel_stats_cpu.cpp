// lsdtt_lithology_channel_stats_cpu.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html

#include "lsdtt_lithology_channel_stats_cpu.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <map>
#include <numeric>
#include <set>
#include <sstream>
#include <stdexcept>
#include <unordered_map>

#include "../utils/lsdtt_statistics.hpp"
#include "lsdtt_proximity_cpu.hpp"

namespace lsdtt
{

namespace
{

constexpr double k_nan = std::numeric_limits<double>::quiet_NaN();
constexpr std::int32_t k_code_unmapped = -1;
constexpr std::int32_t k_code_all = -2;

std::string csv_text(const std::string& s)
{
    if (s.find_first_of(",\"\n") == std::string::npos) return s;
    std::string q = "\"";
    for (char c : s) { if (c == '"') q += "\"\""; else q += c; }
    return q + "\"";
}

std::string unit_name(const lsdtt_lithology_units_on_grid& units, std::int32_t code)
{
    if (code == k_code_unmapped) return "unmapped";
    if (code == k_code_all) return "all";
    for (const auto& u : units.units)
        if (u.code == code) return u.name;
    return std::to_string(code);
}

// The codes that get a row, in order: every unit, unmapped, all.
std::vector<std::int32_t> row_codes(const lsdtt_lithology_units_on_grid& units)
{
    std::vector<std::int32_t> c;
    for (const auto& u : units.units) c.push_back(u.code);
    c.push_back(k_code_unmapped);
    c.push_back(k_code_all);
    return c;
}

bool in_row(std::int32_t row_code, std::int32_t node_code)
{
    return row_code == k_code_all || row_code == node_code;
}

// The lower weighted quantile, NaN when there is no weight to take it from.
double weighted_quantile_or_nan(const std::vector<double>& v, const std::vector<double>& w, double q)
{
    double total = 0.0;
    for (double x : w) total += x;
    if (v.empty() || !(total > 0.0)) return k_nan;
    return compute_weighted_quantile(v, w, q);
}

struct reach_quantiles
{
    double median = k_nan, p25 = k_nan, p75 = k_nan, length_weighted_median = k_nan;
};

// Unweighted (D2 convention) and length-weighted quantiles of reach k_sn.
reach_quantiles quantiles_of_reaches(const std::vector<double>& k, const std::vector<double>& length)
{
    reach_quantiles q;
    if (k.empty()) return q;
    std::vector<double> a = k;
    q.median = compute_median_in_place(a);
    a = k;
    q.p25 = compute_quantile_interpolated_in_place(a, 0.25);
    a = k;
    q.p75 = compute_quantile_interpolated_in_place(a, 0.75);
    q.length_weighted_median = weighted_quantile_or_nan(k, length, 0.5);
    return q;
}

double median_or_nan(std::vector<double> v)
{
    return v.empty() ? k_nan : compute_median_in_place(v);
}

} // namespace

lsdtt_lithology_bootstrap_cluster parse_lithology_bootstrap_cluster(const std::string& s)
{
    if (s == "reach") return lsdtt_lithology_bootstrap_cluster::reach;
    if (s == "basin") return lsdtt_lithology_bootstrap_cluster::basin;
    throw std::invalid_argument("lithology_bootstrap_cluster must be reach or basin, not '" + s + "'");
}

std::string lithology_bootstrap_cluster_name(lsdtt_lithology_bootstrap_cluster c)
{
    return c == lsdtt_lithology_bootstrap_cluster::basin ? "basin" : "reach";
}

std::string chi_lithology_reach_end_name(lsdtt_chi_lithology_reach_end e)
{
    switch (e)
    {
        case lsdtt_chi_lithology_reach_end::contact: return "contact";
        case lsdtt_chi_lithology_reach_end::junction: return "junction";
        case lsdtt_chi_lithology_reach_end::source: return "source";
        case lsdtt_chi_lithology_reach_end::outlet: return "outlet";
        case lsdtt_chi_lithology_reach_end::buffer: return "buffer";
    }
    return "unknown";
}

std::vector<std::int32_t> assign_chi_fit_nodes_to_segments(const lsdtt_chi_network_fit_result& r)
{
    std::vector<std::int32_t> seg(r.flats.size(), -1);
    std::unordered_map<std::int64_t, std::size_t> idx;
    idx.reserve(r.flats.size());
    for (std::size_t i = 0; i < r.flats.size(); ++i) idx[r.flats[i]] = i;

    for (std::size_t s = 0; s < r.segments.size(); ++s)
    {
        const auto& sg = r.segments[s];
        std::int64_t cur = sg.upstream_flat;
        for (std::size_t guard = 0; guard <= r.flats.size(); ++guard)
        {
            const auto it = idx.find(cur);
            if (it == idx.end()) break;
            if (cur == sg.downstream_flat) break;   // belongs to the next segment down
            seg[it->second] = static_cast<std::int32_t>(s);
            const std::int64_t rec = r.receiver_flats[it->second];
            if (rec == k_outlet_receiver) break;
            cur = rec;
        }
    }
    // A downstream end with no segment below it (the outlet) belongs to the
    // segment ending there.
    for (std::size_t s = 0; s < r.segments.size(); ++s)
    {
        const auto it = idx.find(r.segments[s].downstream_flat);
        if (it != idx.end() && seg[it->second] < 0) seg[it->second] = static_cast<std::int32_t>(s);
    }
    return seg;
}

std::vector<float> compute_lithology_contact_distance(const lsdtt_raster_cpu& grid,
                                                      const lsdtt_lithology_units_on_grid& units,
                                                      int n_threads)
{
    const lsdtt_raster_meta& m = grid.get_meta();
    const int nr = m.n_rows, nc = m.n_cols;
    const std::size_t n = static_cast<std::size_t>(nr) * nc;
    if (units.code_of_pixel.size() != n)
        throw std::invalid_argument("compute_lithology_contact_distance: units are not on this grid");
    const std::vector<std::uint8_t>& mask = grid.get_mask();
    std::vector<float> src(n, 0.0f);
    std::size_t n_contact = 0;
    for (int r = 0; r < nr; ++r)
        for (int c = 0; c < nc; ++c)
        {
            const std::size_t i = static_cast<std::size_t>(r) * nc + c;
            if (!mask[i]) continue;
            const std::int32_t me = units.code_of_pixel[i];
            const int dr[4] = {-1, 1, 0, 0}, dc[4] = {0, 0, -1, 1};
            for (int k = 0; k < 4; ++k)
            {
                const int rr = r + dr[k], cc = c + dc[k];
                if (rr < 0 || rr >= nr || cc < 0 || cc >= nc) continue;
                const std::size_t j = static_cast<std::size_t>(rr) * nc + cc;
                if (mask[j] && units.code_of_pixel[j] != me) { src[i] = 1.0f; break; }
            }
            if (src[i] > 0.5f) ++n_contact;
        }
    std::vector<float> out(n, std::numeric_limits<float>::quiet_NaN());
    if (n_contact == 0) return out;
    const lsdtt_raster_cpu d = euclidean_distance_transform(
        lsdtt_raster_cpu(std::move(src), std::vector<std::uint8_t>(mask), m), n_threads);
    for (std::size_t i = 0; i < n; ++i)
        if (mask[i] && d.get_mask()[i]) out[i] = d.get_data()[i];
    return out;
}

lsdtt_chi_lithology_node_table build_chi_lithology_node_table(
    const lsdtt_chi_network_fit_collection_result& fits, const lsdtt_raster_meta& m,
    const lsdtt_lithology_units_on_grid& units, const std::vector<float>& contact_distance)
{
    const std::size_t n = static_cast<std::size_t>(m.n_rows) * m.n_cols;
    if (units.code_of_pixel.size() != n || contact_distance.size() != n)
        throw std::invalid_argument("build_chi_lithology_node_table: units or contact distances "
                                    "are not on the DEM grid");
    lsdtt_chi_lithology_node_table t;
    for (std::size_t b = 0; b < fits.basins.size(); ++b)
    {
        const auto& r = fits.basins[b];
        const std::vector<std::int32_t> seg = assign_chi_fit_nodes_to_segments(r);
        for (std::size_t i = 0; i < r.flats.size(); ++i)
        {
            lsdtt_chi_lithology_node_row row;
            row.basin_index = static_cast<std::int32_t>(b);
            row.node_flat = r.flats[i];
            row.row = static_cast<std::int32_t>(r.flats[i] / m.n_cols);
            row.col = static_cast<std::int32_t>(r.flats[i] % m.n_cols);
            row.x = m.x_min + (row.col + 0.5) * m.cell_size;
            row.y = m.y_min + (m.n_rows - row.row - 0.5) * m.cell_size;
            row.chi = r.chi[i];
            row.elevation = r.observed_elevation[i];
            row.segment_id = seg[i] >= 0 ? r.segments[static_cast<std::size_t>(seg[i])].segment_id : -1;
            row.k_sn = seg[i] >= 0 ? r.segments[static_cast<std::size_t>(seg[i])].k_sn
                                   : std::numeric_limits<float>::quiet_NaN();
            const std::int64_t rec = r.receiver_flats[i];
            if (rec != k_outlet_receiver && rec >= 0)
            {
                const double dr = static_cast<double>(rec / m.n_cols - row.row);
                const double dc = static_cast<double>(rec % m.n_cols - row.col);
                row.flow_length_to_receiver_m = static_cast<float>(std::hypot(dr, dc) * m.cell_size);
            }
            if (i < r.branch_id_of_node.size()) row.branch_id = r.branch_id_of_node[i];
            const std::size_t f = static_cast<std::size_t>(r.flats[i]);
            row.lith_code = units.code_of_pixel[f];
            row.distance_to_contact_m = contact_distance[f];
            t.rows.push_back(row);
        }
    }
    return t;
}

lsdtt_chi_ksn_lithology_summary summarise_chi_ksn_by_lithology(
    const lsdtt_chi_lithology_node_table& nodes, const lsdtt_lithology_units_on_grid& units,
    double m_over_n, double A_0)
{
    lsdtt_chi_ksn_lithology_summary s;
    s.m_over_n = m_over_n;
    s.A_0 = A_0;
    for (const auto& r : nodes.rows)
    {
        if (r.segment_id < 0) continue;
        ++s.n_nodes_total;
        if (r.lith_code == k_code_unmapped) ++s.n_nodes_unmapped;
    }
    for (std::int32_t code : row_codes(units))
    {
        lsdtt_chi_ksn_lithology_summary_row row;
        row.lith_code = code;
        row.lith_name = unit_name(units, code);
        std::vector<double> k;
        std::set<std::pair<std::int32_t, std::int32_t>> segs;
        std::set<std::int32_t> basins;
        double wsum = 0.0, wk = 0.0;
        for (const auto& r : nodes.rows)
        {
            if (r.segment_id < 0 || !in_row(code, r.lith_code)) continue;
            k.push_back(r.k_sn);
            segs.insert({r.basin_index, r.segment_id});
            basins.insert(r.basin_index);
            wsum += r.flow_length_to_receiver_m;
            wk += r.flow_length_to_receiver_m * r.k_sn;
        }
        row.n_nodes = k.size();
        row.n_segments = segs.size();
        row.n_basins = basins.size();
        row.channel_length_m = wsum;
        row.low_sample = row.n_segments < 5;
        if (k.empty())
        {
            row.k_sn_median = row.k_sn_q25 = row.k_sn_q75 = k_nan;
            row.k_sn_mean = row.k_sn_min = row.k_sn_max = row.k_sn_length_weighted_mean = k_nan;
        }
        else
        {
            row.k_sn_median = compute_median(k);
            row.k_sn_q25 = compute_percentile(k, 25.0);
            row.k_sn_q75 = compute_percentile(k, 75.0);
            row.k_sn_mean = compute_mean(k);
            row.k_sn_min = *std::min_element(k.begin(), k.end());
            row.k_sn_max = *std::max_element(k.begin(), k.end());
            row.k_sn_length_weighted_mean = wsum > 0.0 ? wk / wsum : k_nan;
        }
        s.rows.push_back(row);
    }
    return s;
}

lsdtt_chi_ksn_histogram make_chi_ksn_histogram_by_lithology(
    const lsdtt_chi_lithology_node_table& nodes, const lsdtt_lithology_units_on_grid& units,
    double m_over_n, double A_0, const lsdtt_chi_ksn_histogram_options& opts)
{
    if (opts.n_bins < 1) throw std::invalid_argument("k_sn histogram: n_bins must be >= 1");
    lsdtt_chi_ksn_histogram h;
    h.m_over_n = m_over_n;
    h.A_0 = A_0;

    auto tr = [&](double v) { return opts.log_bins ? std::log10(v) : v; };
    std::vector<double> pooled;
    for (const auto& r : nodes.rows)
    {
        if (r.segment_id < 0 || !std::isfinite(r.k_sn)) continue;
        if (opts.log_bins && !(r.k_sn > 0.0f)) continue;
        pooled.push_back(tr(r.k_sn));
    }
    double lo = k_nan, hi = k_nan;
    if (!pooled.empty())
    {
        lo = compute_percentile(pooled, 1.0);
        hi = compute_percentile(pooled, 99.0);
        if (!(hi > lo))   // one value: a unit-wide bin around it
        {
            lo -= 0.5;
            hi += 0.5;
        }
    }
    h.range_lower = opts.log_bins ? std::pow(10.0, lo) : lo;
    h.range_upper = opts.log_bins ? std::pow(10.0, hi) : hi;
    const double w = (hi - lo) / opts.n_bins;
    auto edge = [&](int b) {
        const double e = lo + b * w;
        return opts.log_bins ? std::pow(10.0, e) : e;
    };

    for (std::int32_t code : row_codes(units))
    {
        std::vector<std::size_t> count(static_cast<std::size_t>(opts.n_bins) + 2, 0);
        std::vector<double> length(count.size(), 0.0);
        for (const auto& r : nodes.rows)
        {
            if (r.segment_id < 0 || !std::isfinite(r.k_sn) || !in_row(code, r.lith_code)) continue;
            std::size_t slot;   // 0 below, 1..n_bins, n_bins+1 above
            if (pooled.empty() || (opts.log_bins && !(r.k_sn > 0.0f))) slot = 0;
            else
            {
                const double v = tr(r.k_sn);
                if (v < lo) slot = 0;
                else if (v > hi) slot = count.size() - 1;
                else slot = 1 + std::min<std::size_t>(static_cast<std::size_t>(opts.n_bins) - 1,
                                                      static_cast<std::size_t>((v - lo) / w));
            }
            ++count[slot];
            length[slot] += r.flow_length_to_receiver_m;
        }
        const std::string name = unit_name(units, code);
        for (std::size_t slot = 0; slot < count.size(); ++slot)
        {
            lsdtt_chi_ksn_histogram_row row;
            row.lith_code = code;
            row.lith_name = name;
            row.bin_index = static_cast<std::int32_t>(slot) - 1;
            const double inf = std::numeric_limits<double>::infinity();
            if (slot == 0) { row.bin_lower = -inf; row.bin_upper = h.range_lower; }
            else if (slot == count.size() - 1) { row.bin_lower = h.range_upper; row.bin_upper = inf; }
            else
            {
                row.bin_lower = edge(static_cast<int>(slot) - 1);
                row.bin_upper = edge(static_cast<int>(slot));
            }
            row.n_nodes = count[slot];
            row.channel_length_m = length[slot];
            h.rows.push_back(row);
        }
    }
    return h;
}

lsdtt_lithology_tag summarise_lithology_of_flats(const std::vector<std::int64_t>& flats,
                                                 const std::vector<double>& weights,
                                                 const std::vector<std::int32_t>& code_of_pixel)
{
    if (flats.size() != weights.size())
        throw std::invalid_argument("summarise_lithology_of_flats: flats and weights differ in size");
    lsdtt_lithology_tag tag;
    if (flats.empty()) return tag;
    double total = 0.0;
    for (double w : weights) total += w;
    const bool count_nodes = !(total > 0.0);
    std::map<std::int32_t, double> by_code;   // ascending, so ties go to the lower code
    for (std::size_t i = 0; i < flats.size(); ++i)
    {
        const auto f = static_cast<std::size_t>(flats[i]);
        if (f >= code_of_pixel.size())
            throw std::invalid_argument("summarise_lithology_of_flats: a node is off the grid");
        by_code[code_of_pixel[f]] += count_nodes ? 1.0 : weights[i];
    }
    if (count_nodes) total = static_cast<double>(flats.size());
    double best = -1.0;
    for (const auto& [code, w] : by_code)
        if (w > best)
        {
            best = w;
            tag.code = code;
        }
    tag.majority_fraction = best / total;
    tag.n_units = static_cast<std::int32_t>(by_code.size());
    return tag;
}

lsdtt_chi_lithology_tables build_chi_lithology_tables(
    const lsdtt_chi_network_fit_collection_result& fits, const lsdtt_raster_meta& m,
    const lsdtt_lithology_units_on_grid& units, const std::vector<float>& contact_distance,
    const std::vector<float>& drainage_area_m2, const lsdtt_chi_lithology_options& opts)
{
    const std::size_t n_grid = static_cast<std::size_t>(m.n_rows) * m.n_cols;
    if (!drainage_area_m2.empty() && drainage_area_m2.size() != n_grid)
        throw std::invalid_argument("build_chi_lithology_tables: drainage area is not on the DEM grid");
    if (!(opts.contact_buffer_m >= 0.0) || !std::isfinite(opts.contact_buffer_m))
        throw std::invalid_argument("build_chi_lithology_tables: the contact buffer must be >= 0");
    if (opts.min_reach_pixels < 1)
        throw std::invalid_argument("build_chi_lithology_tables: min_reach_pixels must be >= 1");

    lsdtt_chi_lithology_tables t;
    t.nodes = build_chi_lithology_node_table(fits, m, units, contact_distance);
    t.contact_buffer_m = opts.contact_buffer_m;
    t.min_reach_pixels = opts.min_reach_pixels;
    using end_t = lsdtt_chi_lithology_reach_end;

    std::size_t base = 0;   // first node row of this basin (rows follow r.flats)
    for (std::size_t b = 0; b < fits.basins.size(); ++b)
    {
        const auto& r = fits.basins[b];
        const std::size_t n = r.flats.size();
        auto row = [&](std::size_t i) -> lsdtt_chi_lithology_node_row& { return t.nodes.rows[base + i]; };

        std::unordered_map<std::int64_t, std::size_t> idx;
        idx.reserve(n);
        for (std::size_t i = 0; i < n; ++i) idx[r.flats[i]] = i;

        std::vector<std::uint8_t> near(n, 0);
        for (std::size_t i = 0; i < n; ++i)
        {
            auto& nr = row(i);
            if (!drainage_area_m2.empty())
                nr.drainage_area_m2 = drainage_area_m2[static_cast<std::size_t>(r.flats[i])];
            nr.near_contact = opts.contact_buffer_m > 0.0 && std::isfinite(nr.distance_to_contact_m) &&
                              nr.distance_to_contact_m < opts.contact_buffer_m;
            near[i] = nr.near_contact ? 1 : 0;
        }

        // Receiver (index in this network, or -1) and the join to it.
        const std::size_t none = static_cast<std::size_t>(-1);
        std::vector<std::size_t> rec(n, none), down(n, none), up(n, none);
        std::vector<std::uint8_t> has_donor(n, 0), donor_near(n, 0), donor_other_code(n, 0);
        for (std::size_t i = 0; i < n; ++i)
        {
            const std::int64_t rf = r.receiver_flats[i];
            if (rf == k_outlet_receiver || rf < 0) continue;
            const auto it = idx.find(rf);
            if (it == idx.end()) continue;
            const std::size_t j = it->second;
            rec[i] = j;
            has_donor[j] = 1;
            if (near[i]) donor_near[j] = 1;
            else if (row(i).lith_code != row(j).lith_code) donor_other_code[j] = 1;
            if (!near[i] && !near[j] && row(i).lith_code == row(j).lith_code &&
                row(i).branch_id == row(j).branch_id && up[j] == none)
            {
                down[i] = j;
                up[j] = i;
            }
        }

        // Reaches, numbered in the order of their downstream ends.
        for (std::size_t i = 0; i < n; ++i)
        {
            if (near[i] || down[i] != none) continue;
            std::vector<std::size_t> chain;   // downstream end first
            for (std::size_t k = i; k != none; k = up[k]) chain.push_back(k);
            const std::size_t d = chain.front(), u = chain.back();

            lsdtt_chi_lithology_reach_row rr;
            rr.reach_id = static_cast<std::int32_t>(t.reaches.size());
            rr.basin_index = static_cast<std::int32_t>(b);
            rr.branch_id = row(d).branch_id;
            rr.lith_code = row(d).lith_code;
            rr.n_nodes = chain.size();
            rr.upstream_flat = r.flats[u];
            rr.downstream_flat = r.flats[d];
            rr.chi_span = static_cast<double>(r.chi[u]) - r.chi[d];
            rr.dz_m = static_cast<double>(r.observed_elevation[u]) - r.observed_elevation[d];

            if (rec[d] == none) rr.downstream_end = end_t::outlet;
            else if (near[rec[d]]) rr.downstream_end = end_t::buffer;
            else if (row(rec[d]).lith_code != rr.lith_code) rr.downstream_end = end_t::contact;
            else rr.downstream_end = end_t::junction;

            if (!has_donor[u]) rr.upstream_end = end_t::source;
            else if (donor_near[u]) rr.upstream_end = end_t::buffer;
            else if (donor_other_code[u]) rr.upstream_end = end_t::contact;
            else rr.upstream_end = end_t::junction;

            double mx = 0.0, my = 0.0, wk = 0.0, wsum = 0.0;
            std::vector<double> areas;
            for (std::size_t k : chain)
            {
                auto& nr = row(k);
                nr.reach_id = rr.reach_id;
                rr.length_m += nr.flow_length_to_receiver_m;
                mx += r.chi[k];
                my += r.observed_elevation[k];
                if (std::isfinite(nr.drainage_area_m2)) areas.push_back(nr.drainage_area_m2);
                const double dchi = rec[k] == none ? 0.0 : static_cast<double>(r.chi[k]) - r.chi[rec[k]];
                if (std::isfinite(nr.k_sn) && dchi > 0.0)
                {
                    wk += dchi * nr.k_sn;
                    wsum += dchi;
                }
            }
            rr.median_drainage_area_m2 = median_or_nan(std::move(areas));
            if (wsum > 0.0) rr.k_sn_fit_chi_weighted = wk / wsum;

            // Ordinary least squares of elevation on chi, about the means.
            const double nn = static_cast<double>(chain.size());
            mx /= nn;
            my /= nn;
            double sxx = 0.0, sxy = 0.0, syy = 0.0;
            for (std::size_t k : chain)
            {
                const double dx = r.chi[k] - mx, dy = r.observed_elevation[k] - my;
                sxx += dx * dx;
                sxy += dx * dy;
                syy += dy * dy;
            }
            if (sxx > 0.0)
            {
                const double slope = sxy / sxx, icpt = my - slope * mx;
                double ssr = 0.0;
                for (std::size_t k : chain)
                {
                    const double e = r.observed_elevation[k] - (icpt + slope * r.chi[k]);
                    ssr += e * e;
                }
                rr.k_sn_ols = slope;
                rr.r_squared = syy > 0.0 ? 1.0 - ssr / syy : k_nan;
                if (chain.size() > 2) rr.k_sn_ols_se = std::sqrt(ssr / (nn - 2.0) / sxx);
            }
            rr.qualifies = static_cast<int>(rr.n_nodes) >= opts.min_reach_pixels &&
                           rr.chi_span > 0.0 && std::isfinite(rr.k_sn_ols);
            t.reaches.push_back(rr);
        }

        // Segments: the nodes each owns, length-weighted unit.
        const std::vector<std::int32_t> seg_of_node = assign_chi_fit_nodes_to_segments(r);
        std::vector<std::vector<std::int64_t>> seg_flats(r.segments.size());
        std::vector<std::vector<double>> seg_w(r.segments.size());
        for (std::size_t i = 0; i < n; ++i)
        {
            if (seg_of_node[i] < 0) continue;
            const auto sg = static_cast<std::size_t>(seg_of_node[i]);
            seg_flats[sg].push_back(r.flats[i]);
            seg_w[sg].push_back(row(i).flow_length_to_receiver_m);
        }
        for (std::size_t sg = 0; sg < r.segments.size(); ++sg)
        {
            lsdtt_chi_lithology_segment_row sr;
            sr.basin_index = static_cast<std::int32_t>(b);
            sr.segment_id = r.segments[sg].segment_id;
            sr.branch_id = r.segments[sg].branch_id;
            sr.k_sn = r.segments[sg].k_sn;
            sr.n_nodes = static_cast<std::int32_t>(seg_flats[sg].size());
            for (double w : seg_w[sg]) sr.length_m += w;
            sr.lith = summarise_lithology_of_flats(seg_flats[sg], seg_w[sg], units.code_of_pixel);
            t.segments.push_back(sr);
        }
        base += n;
    }
    return t;
}

std::mt19937_64 lithology_bootstrap_engine(std::int64_t seed, double m_over_n)
{
    const long long mn = std::llround(m_over_n * 10000.0);
    std::seed_seq seq{static_cast<std::uint32_t>(seed), static_cast<std::uint32_t>(mn)};
    return std::mt19937_64(seq);
}

lsdtt_lithology_median_ci bootstrap_lithology_median_ci(const std::vector<double>& values,
                                                        const std::vector<std::int32_t>& cluster_of_value,
                                                        int n_replicates, std::mt19937_64& rng)
{
    lsdtt_lithology_median_ci ci;
    if (values.empty() || n_replicates <= 0) return ci;
    if (!cluster_of_value.empty() && cluster_of_value.size() != values.size())
        throw std::invalid_argument("bootstrap_lithology_median_ci: one cluster per value");

    std::vector<std::vector<double>> clusters;
    if (cluster_of_value.empty())
        for (double v : values) clusters.push_back({v});
    else
    {
        std::map<std::int32_t, std::vector<double>> by;   // ascending cluster id
        for (std::size_t i = 0; i < values.size(); ++i) by[cluster_of_value[i]].push_back(values[i]);
        for (auto& kv : by) clusters.push_back(std::move(kv.second));
    }

    std::uniform_int_distribution<std::size_t> pick(0, clusters.size() - 1);
    std::vector<double> medians;
    medians.reserve(static_cast<std::size_t>(n_replicates));
    std::vector<double> pooled;
    for (int rep = 0; rep < n_replicates; ++rep)
    {
        pooled.clear();
        for (std::size_t c = 0; c < clusters.size(); ++c)
        {
            const auto& cl = clusters[pick(rng)];
            pooled.insert(pooled.end(), cl.begin(), cl.end());
        }
        medians.push_back(compute_median_in_place(pooled));
    }
    std::vector<double> a = medians;
    ci.low = compute_quantile_interpolated_in_place(a, 0.025);
    a = medians;
    ci.high = compute_quantile_interpolated_in_place(a, 0.975);
    ci.n_replicates = n_replicates;
    return ci;
}

lsdtt_chi_ksn_lithology_summary summarise_chi_ksn_by_lithology(
    const lsdtt_chi_lithology_tables& t, const lsdtt_lithology_units_on_grid& units,
    double m_over_n, double A_0, const lsdtt_chi_lithology_options& opts, std::mt19937_64& rng)
{
    lsdtt_chi_ksn_lithology_summary s = summarise_chi_ksn_by_lithology(t.nodes, units, m_over_n, A_0);
    s.contact_buffer_m = t.contact_buffer_m;
    s.min_reach_pixels = t.min_reach_pixels;
    s.bootstrap_cluster = opts.bootstrap_cluster;

    for (auto& row : s.rows)
    {
        const std::int32_t code = row.lith_code;
        std::vector<double> nk, nw, area;
        for (const auto& r : t.nodes.rows)
        {
            if (r.segment_id < 0 || !in_row(code, r.lith_code)) continue;
            if (r.near_contact)
            {
                ++row.n_nodes_excluded_near_contact;
                continue;
            }
            if (std::isfinite(r.k_sn))
            {
                nk.push_back(r.k_sn);
                nw.push_back(r.flow_length_to_receiver_m);
            }
            if (std::isfinite(r.drainage_area_m2)) area.push_back(r.drainage_area_m2);
        }
        row.median_drainage_area_m2 = median_or_nan(std::move(area));
        row.node_k_sn_length_weighted_median = weighted_quantile_or_nan(nk, nw, 0.5);
        row.node_k_sn_length_weighted_p25 = weighted_quantile_or_nan(nk, nw, 0.25);
        row.node_k_sn_length_weighted_p75 = weighted_quantile_or_nan(nk, nw, 0.75);

        std::vector<double> rk, rl;
        std::vector<std::int32_t> rb;
        for (const auto& r : t.reaches)
        {
            if (!in_row(code, r.lith_code)) continue;
            ++row.n_reaches;
            if (!r.qualifies) continue;
            rk.push_back(r.k_sn_ols);
            rl.push_back(r.length_m);
            rb.push_back(r.basin_index);
        }
        row.n_reaches_qualifying = rk.size();
        const reach_quantiles q = quantiles_of_reaches(rk, rl);
        row.reach_k_sn_median = q.median;
        row.reach_k_sn_p25 = q.p25;
        row.reach_k_sn_p75 = q.p75;
        row.reach_k_sn_length_weighted_median = q.length_weighted_median;
        row.low_sample = rk.size() < 5;
        if (!row.low_sample)
        {
            const std::vector<std::int32_t> none;
            const auto ci = bootstrap_lithology_median_ci(
                rk, opts.bootstrap_cluster == lsdtt_lithology_bootstrap_cluster::basin ? rb : none,
                opts.bootstrap_replicates, rng);
            row.reach_k_sn_median_ci_low = ci.low;
            row.reach_k_sn_median_ci_high = ci.high;
            row.n_bootstrap = ci.n_replicates;
        }
    }

    // Unit x basin: every (unit, basin) with a channel node; units, then unmapped.
    for (std::int32_t code : row_codes(units))
    {
        if (code == k_code_all) continue;
        std::map<std::int32_t, lsdtt_chi_ksn_lithology_basin_row> by_basin;
        std::map<std::int32_t, std::vector<double>> k_all, area, rk, rl;
        for (const auto& r : t.nodes.rows)
        {
            if (r.segment_id < 0 || r.lith_code != code) continue;
            auto& br = by_basin[r.basin_index];
            ++br.n_nodes;
            br.channel_length_m += r.flow_length_to_receiver_m;
            k_all[r.basin_index].push_back(r.k_sn);
            if (r.near_contact) ++br.n_nodes_excluded_near_contact;
            else if (std::isfinite(r.drainage_area_m2)) area[r.basin_index].push_back(r.drainage_area_m2);
        }
        for (const auto& r : t.reaches)
        {
            if (r.lith_code != code) continue;
            auto it = by_basin.find(r.basin_index);
            if (it == by_basin.end()) continue;
            ++it->second.n_reaches;
            if (!r.qualifies) continue;
            rk[r.basin_index].push_back(r.k_sn_ols);
            rl[r.basin_index].push_back(r.length_m);
        }
        for (auto& [basin, br] : by_basin)
        {
            br.lith_code = code;
            br.lith_name = unit_name(units, code);
            br.basin_index = basin;
            br.k_sn_median = compute_median(k_all[basin]);
            br.median_drainage_area_m2 = median_or_nan(area[basin]);
            br.n_reaches_qualifying = rk[basin].size();
            const reach_quantiles q = quantiles_of_reaches(rk[basin], rl[basin]);
            br.reach_k_sn_median = q.median;
            br.reach_k_sn_p25 = q.p25;
            br.reach_k_sn_p75 = q.p75;
            br.reach_k_sn_length_weighted_median = q.length_weighted_median;
            br.low_sample = br.n_reaches_qualifying < 5;
            s.basin_rows.push_back(br);
        }
    }
    return s;
}

// ---------------------------------------------------------------------------
// Channel-network records
// ---------------------------------------------------------------------------
std::vector<double> lithology_node_step_lengths(const lsdtt_flow_base& flow,
                                                const std::vector<std::int64_t>& flats)
{
    std::vector<double> w(flats.size(), 0.0);
    const double dx = flow.cell_size();
    for (std::size_t i = 0; i < flats.size(); ++i)
    {
        const std::int64_t f = flats[i];
        if (!flow.is_valid_cell(f) || flow.is_outlet(f)) continue;
        const std::int64_t rc = flow.receiver_of(f);
        if (rc < 0 || rc == f) continue;
        const int dr = flow.row_from_flat(rc) - flow.row_from_flat(f);
        const int dc = flow.col_from_flat(rc) - flow.col_from_flat(f);
        w[i] = (dr != 0 && dc != 0) ? dx * std::sqrt(2.0) : dx;
    }
    return w;
}

lsdtt_lithology_tag lithology_of_channel_nodes(const lsdtt_flow_base& flow,
                                               const std::vector<std::int64_t>& flats,
                                               const lsdtt_lithology_units_on_grid& units)
{
    return summarise_lithology_of_flats(flats, lithology_node_step_lengths(flow, flats),
                                        units.code_of_pixel);
}

std::vector<lsdtt_lithology_tag> tag_channel_link_bearings_with_lithology(
    const lsdtt_flow_base& flow, const std::vector<lsdtt_channel_link_bearing>& records,
    const lsdtt_lithology_units_on_grid& units)
{
    std::vector<lsdtt_lithology_tag> tags(records.size());
    for (std::size_t i = 0; i < records.size(); ++i)
        if (!records[i].empty()) tags[i] = lithology_of_channel_nodes(flow, records[i].sampled_flats, units);
    return tags;
}

std::vector<lsdtt_junction_lithology> tag_junction_angles_with_lithology(
    const lsdtt_flow_base& flow, const lsdtt_junction_network_base& network,
    const std::vector<lsdtt_junction_angle_record>& records,
    const lsdtt_lithology_units_on_grid& units, double min_majority_fraction)
{
    std::vector<lsdtt_junction_lithology> tags(records.size());
    for (std::size_t i = 0; i < records.size(); ++i)
    {
        const auto& rec = records[i];
        if (rec.empty()) continue;
        auto& t = tags[i];
        t.valid = true;
        const std::int64_t f = network.flat_of_junction(rec.junction);
        t.junction_code = (f >= 0 && static_cast<std::size_t>(f) < units.code_of_pixel.size())
                              ? units.code_of_pixel[static_cast<std::size_t>(f)] : -1;
        t.donor_1 = lithology_of_channel_nodes(flow, rec.donor_1_sampled_flats, units);
        t.donor_2 = lithology_of_channel_nodes(flow, rec.donor_2_sampled_flats, units);
        t.receiver = lithology_of_channel_nodes(flow, rec.receiver_sampled_flats, units);
        bool agree = true;
        for (const lsdtt_lithology_tag* a : {&t.donor_1, &t.donor_2, &t.receiver})
            agree = agree && a->n_units > 0 && a->code == t.junction_code &&
                    a->majority_fraction >= min_majority_fraction;
        t.all_arms_agree = agree;
    }
    return tags;
}

namespace
{

std::int32_t category_of_code(std::int32_t code)
{
    return code == k_code_unmapped ? k_lithology_category_unmapped : code;
}

std::map<int, std::string> category_names(const lsdtt_lithology_units_on_grid& units)
{
    std::map<int, std::string> names;
    for (const auto& u : units.units) names[u.code] = u.name;
    names[k_lithology_category_unmapped] = "unmapped";
    return names;
}

} // namespace

std::string lithology_mixed_category_name(double min_majority_fraction)
{
    std::ostringstream os;
    os << std::setprecision(6) << 100.0 * min_majority_fraction;
    return "mixed (<" + os.str() + "% one unit)";
}

channel_record_categories lithology_categories_of_links(const std::vector<lsdtt_lithology_tag>& tags,
                                                        const lsdtt_lithology_units_on_grid& units,
                                                        double min_majority_fraction,
                                                        const std::string& key)
{
    channel_record_categories c;
    c.key = key;
    c.names = category_names(units);
    c.names[k_lithology_category_mixed] = lithology_mixed_category_name(min_majority_fraction);
    c.descriptive_only.insert(k_lithology_category_mixed);
    c.category.assign(tags.size(), -1);
    for (std::size_t i = 0; i < tags.size(); ++i)
    {
        if (tags[i].n_units == 0) continue;
        c.category[i] = tags[i].majority_fraction >= min_majority_fraction
                            ? category_of_code(tags[i].code)
                            : k_lithology_category_mixed;
    }
    return c;
}

std::map<std::int32_t, std::size_t> lithology_mixed_links_by_majority_unit(
    const std::vector<lsdtt_lithology_tag>& tags, double min_majority_fraction)
{
    std::map<std::int32_t, std::size_t> out;
    for (const auto& t : tags)
        if (t.n_units > 0 && !(t.majority_fraction >= min_majority_fraction)) ++out[t.code];
    return out;
}

channel_record_categories lithology_categories_of_junctions(
    const std::vector<lsdtt_junction_lithology>& tags, const lsdtt_lithology_units_on_grid& units,
    bool interior_only, const std::string& key)
{
    channel_record_categories c;
    c.key = key;
    c.names = category_names(units);
    if (interior_only)
    {
        c.names[k_lithology_category_arms_disagree] = "arms disagree";
        c.descriptive_only.insert(k_lithology_category_arms_disagree);
    }
    c.category.assign(tags.size(), -1);
    for (std::size_t i = 0; i < tags.size(); ++i)
    {
        if (!tags[i].valid) continue;
        c.category[i] = (!interior_only || tags[i].all_arms_agree)
                            ? category_of_code(tags[i].junction_code)
                            : k_lithology_category_arms_disagree;
    }
    return c;
}

lsdtt_channel_record_extra_columns lithology_columns_of_links(const std::vector<lsdtt_lithology_tag>& tags,
                                                             const lsdtt_lithology_units_on_grid& units)
{
    lsdtt_channel_record_extra_columns x;
    x.names = {"lith_code", "lith_name", "lith_majority_fraction"};
    x.vector_names = {"lith_code", "lith_name", "lith_frac"};
    x.types = {lsdtt_vector_field_type::int64_type, lsdtt_vector_field_type::string_type,
               lsdtt_vector_field_type::double_type};
    for (const auto& t : tags)
        x.values.push_back({static_cast<std::int64_t>(t.code),
                            t.n_units > 0 ? unit_name(units, t.code) : std::string(),
                            t.majority_fraction});
    return x;
}

lsdtt_channel_record_extra_columns lithology_columns_of_junctions(
    const std::vector<lsdtt_junction_lithology>& tags, const lsdtt_lithology_units_on_grid& units)
{
    using ft = lsdtt_vector_field_type;
    lsdtt_channel_record_extra_columns x;
    x.names = {"junction_lith_code", "junction_lith_name", "donor_1_lith_code", "donor_1_lith_fraction",
               "donor_2_lith_code", "donor_2_lith_fraction", "receiver_lith_code",
               "receiver_lith_fraction", "lith_all_arms_agree"};
    x.vector_names = {"j_lith", "j_lithname", "d1_lith", "d1_frac", "d2_lith", "d2_frac",
                      "r_lith", "r_frac", "lith_agree"};
    x.types = {ft::int64_type, ft::string_type, ft::int64_type, ft::double_type, ft::int64_type,
               ft::double_type, ft::int64_type, ft::double_type, ft::int64_type};
    for (const auto& t : tags)
        x.values.push_back({static_cast<std::int64_t>(t.junction_code),
                            t.valid ? unit_name(units, t.junction_code) : std::string(),
                            static_cast<std::int64_t>(t.donor_1.code), t.donor_1.majority_fraction,
                            static_cast<std::int64_t>(t.donor_2.code), t.donor_2.majority_fraction,
                            static_cast<std::int64_t>(t.receiver.code), t.receiver.majority_fraction,
                            static_cast<std::int64_t>(t.all_arms_agree ? 1 : 0)});
    return x;
}

std::vector<std::uint8_t> pixels_draining_through(const lsdtt_flow_base& flow,
                                                  const std::vector<std::int64_t>& seed_flats)
{
    const std::size_t n = flow.full_grid_size();
    // 0 unknown, 1 drains through a seed, 2 does not.
    std::vector<std::uint8_t> state(n, 0);
    for (const std::int64_t f : seed_flats)
        if (f >= 0 && static_cast<std::size_t>(f) < n) state[static_cast<std::size_t>(f)] = 1;
    std::vector<std::int64_t> path;
    for (std::size_t p = 0; p < n; ++p)
    {
        if (state[p] != 0) continue;
        path.clear();
        std::int64_t cur = static_cast<std::int64_t>(p);
        std::uint8_t result = 2;
        while (true)
        {
            const auto c = static_cast<std::size_t>(cur);
            if (state[c] != 0) { result = state[c]; break; }
            path.push_back(cur);
            if (!flow.is_valid_cell(cur) || flow.is_outlet(cur)) break;
            const std::int64_t rc = flow.receiver_of(cur);
            if (rc < 0 || rc == cur || static_cast<std::size_t>(rc) >= n) break;
            cur = rc;
        }
        for (const std::int64_t q : path) state[static_cast<std::size_t>(q)] = result;
    }
    std::vector<std::uint8_t> out(n, 0);
    for (std::size_t p = 0; p < n; ++p) out[p] = state[p] == 1 && flow.is_valid_cell(static_cast<std::int64_t>(p)) ? 1 : 0;
    return out;
}

void write_chi_lithology_nodes_csv(const std::string& path,
                                   const lsdtt_chi_lithology_node_table& nodes,
                                   const lsdtt_lithology_units_on_grid& units)
{
    std::ofstream out(path);
    if (!out) throw std::runtime_error("cannot write " + path);
    std::map<std::int32_t, std::string> names;
    for (const auto& u : units.units) names[u.code] = csv_text(u.name);
    names[k_code_unmapped] = "unmapped";
    out << "basin_index,node_flat,row,col,x,y,chi,elevation,segment_id,k_sn,"
           "flow_length_to_receiver_m,lith_code,lith_name,distance_to_contact_m,"
           "latitude,longitude,drainage_area_m2,branch_id,near_contact,reach_id\n";
    out << std::setprecision(10);
    for (const auto& r : nodes.rows)
        out << r.basin_index << "," << r.node_flat << "," << r.row << "," << r.col << ","
            << r.x << "," << r.y << "," << r.chi << "," << r.elevation << "," << r.segment_id
            << "," << r.k_sn << "," << r.flow_length_to_receiver_m << "," << r.lith_code << ","
            << names[r.lith_code] << "," << r.distance_to_contact_m << "," << r.latitude << ","
            << r.longitude << "," << r.drainage_area_m2 << "," << r.branch_id << ","
            << (r.near_contact ? 1 : 0) << "," << r.reach_id << "\n";
}

void write_chi_ksn_lithology_summary_csv(const std::string& path,
                                         const lsdtt_chi_ksn_lithology_summary& s)
{
    std::ofstream out(path);
    if (!out) throw std::runtime_error("cannot write " + path);
    out << "lith_code,lith_name,n_nodes,n_segments,n_basins,channel_length_m,k_sn_median,"
           "k_sn_q25,k_sn_q75,k_sn_mean,k_sn_min,k_sn_max,k_sn_length_weighted_mean,low_sample,"
           "m_over_n,A_0,contact_buffer_metres,min_reach_pixels,n_nodes_excluded_near_contact,"
           "n_reaches,n_reaches_qualifying,median_drainage_area_m2,reach_k_sn_median,"
           "reach_k_sn_p25,reach_k_sn_p75,reach_k_sn_length_weighted_median,"
           "reach_k_sn_median_ci_low,reach_k_sn_median_ci_high,n_bootstrap,bootstrap_cluster,"
           "node_k_sn_length_weighted_median,node_k_sn_length_weighted_p25,"
           "node_k_sn_length_weighted_p75\n";
    out << std::setprecision(10);
    const std::string cluster = lithology_bootstrap_cluster_name(s.bootstrap_cluster);
    for (const auto& r : s.rows)
        out << r.lith_code << "," << csv_text(r.lith_name) << "," << r.n_nodes << ","
            << r.n_segments << "," << r.n_basins << "," << r.channel_length_m << ","
            << r.k_sn_median << "," << r.k_sn_q25 << "," << r.k_sn_q75 << "," << r.k_sn_mean
            << "," << r.k_sn_min << "," << r.k_sn_max << "," << r.k_sn_length_weighted_mean << ","
            << (r.low_sample ? 1 : 0) << "," << s.m_over_n << "," << s.A_0 << ","
            << s.contact_buffer_m << "," << s.min_reach_pixels << ","
            << r.n_nodes_excluded_near_contact << "," << r.n_reaches << ","
            << r.n_reaches_qualifying << "," << r.median_drainage_area_m2 << ","
            << r.reach_k_sn_median << "," << r.reach_k_sn_p25 << "," << r.reach_k_sn_p75 << ","
            << r.reach_k_sn_length_weighted_median << "," << r.reach_k_sn_median_ci_low << ","
            << r.reach_k_sn_median_ci_high << "," << r.n_bootstrap << "," << cluster << ","
            << r.node_k_sn_length_weighted_median << "," << r.node_k_sn_length_weighted_p25 << ","
            << r.node_k_sn_length_weighted_p75 << "\n";
}

void write_chi_ksn_histogram_csv(const std::string& path, const lsdtt_chi_ksn_histogram& h)
{
    std::ofstream out(path);
    if (!out) throw std::runtime_error("cannot write " + path);
    out << "lith_code,lith_name,bin_index,bin_lower,bin_upper,n_nodes,channel_length_m,"
           "m_over_n,A_0\n";
    out << std::setprecision(10);
    for (const auto& r : h.rows)
        out << r.lith_code << "," << csv_text(r.lith_name) << "," << r.bin_index << ","
            << r.bin_lower << "," << r.bin_upper << "," << r.n_nodes << "," << r.channel_length_m
            << "," << h.m_over_n << "," << h.A_0 << "\n";
}

void write_chi_lithology_reaches_csv(const std::string& path, const lsdtt_chi_lithology_tables& t,
                                     const lsdtt_lithology_units_on_grid& units)
{
    std::ofstream out(path);
    if (!out) throw std::runtime_error("cannot write " + path);
    out << "reach_id,basin_index,branch_id,lith_code,lith_name,n_nodes,length_metres,chi_span,"
           "dz_metres,median_drainage_area_m2,k_sn_ols,k_sn_ols_se,r_squared,"
           "k_sn_fit_chi_weighted,upstream_flat,downstream_flat,upstream_end,downstream_end,"
           "qualifies\n";
    out << std::setprecision(10);
    for (const auto& r : t.reaches)
        out << r.reach_id << "," << r.basin_index << "," << r.branch_id << "," << r.lith_code
            << "," << csv_text(unit_name(units, r.lith_code)) << "," << r.n_nodes << ","
            << r.length_m << "," << r.chi_span << "," << r.dz_m << "," << r.median_drainage_area_m2
            << "," << r.k_sn_ols << "," << r.k_sn_ols_se << "," << r.r_squared << ","
            << r.k_sn_fit_chi_weighted << "," << r.upstream_flat << "," << r.downstream_flat << ","
            << chi_lithology_reach_end_name(r.upstream_end) << ","
            << chi_lithology_reach_end_name(r.downstream_end) << "," << (r.qualifies ? 1 : 0)
            << "\n";
}

void write_chi_lithology_segments_csv(const std::string& path, const lsdtt_chi_lithology_tables& t,
                                      const lsdtt_lithology_units_on_grid& units)
{
    std::ofstream out(path);
    if (!out) throw std::runtime_error("cannot write " + path);
    out << "basin_index,segment_id,branch_id,k_sn,length_metres,lith_code,lith_name,"
           "lith_majority_fraction,n_lith_units,n_nodes\n";
    out << std::setprecision(10);
    for (const auto& r : t.segments)
        out << r.basin_index << "," << r.segment_id << "," << r.branch_id << "," << r.k_sn << ","
            << r.length_m << "," << r.lith.code << "," << csv_text(unit_name(units, r.lith.code))
            << "," << r.lith.majority_fraction << "," << r.lith.n_units << "," << r.n_nodes << "\n";
}

void write_chi_ksn_lithology_by_basin_csv(const std::string& path,
                                          const lsdtt_chi_ksn_lithology_summary& s)
{
    std::ofstream out(path);
    if (!out) throw std::runtime_error("cannot write " + path);
    out << "lith_code,lith_name,basin_index,n_nodes,n_nodes_excluded_near_contact,channel_length_m,"
           "k_sn_median,n_reaches,n_reaches_qualifying,median_drainage_area_m2,reach_k_sn_median,"
           "reach_k_sn_p25,reach_k_sn_p75,reach_k_sn_length_weighted_median,low_sample,m_over_n,"
           "A_0,contact_buffer_metres,min_reach_pixels\n";
    out << std::setprecision(10);
    for (const auto& r : s.basin_rows)
        out << r.lith_code << "," << csv_text(r.lith_name) << "," << r.basin_index << ","
            << r.n_nodes << "," << r.n_nodes_excluded_near_contact << "," << r.channel_length_m
            << "," << r.k_sn_median << "," << r.n_reaches << "," << r.n_reaches_qualifying << ","
            << r.median_drainage_area_m2 << "," << r.reach_k_sn_median << "," << r.reach_k_sn_p25
            << "," << r.reach_k_sn_p75 << "," << r.reach_k_sn_length_weighted_median << ","
            << (r.low_sample ? 1 : 0) << "," << s.m_over_n << "," << s.A_0 << ","
            << s.contact_buffer_m << "," << s.min_reach_pixels << "\n";
}

} // namespace lsdtt
