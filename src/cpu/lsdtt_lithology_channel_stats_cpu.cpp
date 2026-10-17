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
#include <set>
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

} // namespace

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
           "flow_length_to_receiver_m,lith_code,lith_name,distance_to_contact_m\n";
    out << std::setprecision(10);
    for (const auto& r : nodes.rows)
        out << r.basin_index << "," << r.node_flat << "," << r.row << "," << r.col << ","
            << r.x << "," << r.y << "," << r.chi << "," << r.elevation << "," << r.segment_id
            << "," << r.k_sn << "," << r.flow_length_to_receiver_m << "," << r.lith_code << ","
            << names[r.lith_code] << "," << r.distance_to_contact_m << "\n";
}

void write_chi_ksn_lithology_summary_csv(const std::string& path,
                                         const lsdtt_chi_ksn_lithology_summary& s)
{
    std::ofstream out(path);
    if (!out) throw std::runtime_error("cannot write " + path);
    out << "lith_code,lith_name,n_nodes,n_segments,n_basins,channel_length_m,k_sn_median,"
           "k_sn_q25,k_sn_q75,k_sn_mean,k_sn_min,k_sn_max,k_sn_length_weighted_mean,low_sample,"
           "m_over_n,A_0\n";
    out << std::setprecision(10);
    for (const auto& r : s.rows)
        out << r.lith_code << "," << csv_text(r.lith_name) << "," << r.n_nodes << ","
            << r.n_segments << "," << r.n_basins << "," << r.channel_length_m << ","
            << r.k_sn_median << "," << r.k_sn_q25 << "," << r.k_sn_q75 << "," << r.k_sn_mean
            << "," << r.k_sn_min << "," << r.k_sn_max << "," << r.k_sn_length_weighted_mean << ","
            << (r.low_sample ? 1 : 0) << "," << s.m_over_n << "," << s.A_0 << "\n";
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

} // namespace lsdtt
