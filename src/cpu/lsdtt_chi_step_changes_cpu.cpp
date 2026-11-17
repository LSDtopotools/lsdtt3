// lsdtt_chi_step_changes_cpu.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
//
// -----------------------------------------------------------------------------
// See lsdtt_chi_step_changes_cpu.hpp.
// -----------------------------------------------------------------------------

#include "lsdtt_chi_step_changes_cpu.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>
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

std::string num(double v)
{
    if (!std::isfinite(v)) return "";
    std::ostringstream o;
    o << std::setprecision(10) << v;
    return o.str();
}

} // namespace

// -----------------------------------------------------------------------------
// one path
// -----------------------------------------------------------------------------
std::vector<lsdtt_chi_step_reach_on_path> detect_chi_step_reaches(
    const std::vector<double>&    dz,
    const std::vector<double>&    length_metres,
    const lsdtt_chi_step_options& o)
{
    if (!(o.min_gradient > 0.0))
        throw std::invalid_argument("step_change_min_gradient must be > 0");
    if (!(o.max_window_length_metres > 0.0))
        throw std::invalid_argument("step_change_max_window_length_metres must be > 0");
    if (!(o.min_drop_metres >= 0.0))
        throw std::invalid_argument("step_change_min_drop_metres must be >= 0");
    if (dz.size() != length_metres.size())
        throw std::invalid_argument("detect_chi_step_reaches: dz and length_metres differ in size");

    const std::size_t n = dz.size();
    std::vector<lsdtt_chi_step_reach_on_path> out;
    if (n == 0) return out;

    // Steep windows: start i, end k, flow length <= the maximum (but always
    // the single pair). Marked by a difference array.
    std::vector<long> mark(n + 1, 0);
    auto for_each_window = [&](std::size_t lo, std::size_t hi, auto&& fn)
    {
        for (std::size_t i = lo; i <= hi; ++i)
        {
            double sl = 0.0, sz = 0.0;
            for (std::size_t k = i; k <= hi; ++k)
            {
                sl += length_metres[k];
                if (k > i && sl > o.max_window_length_metres) break;
                sz += dz[k];
                fn(i, k, sz / sl);
            }
        }
    };
    for_each_window(0, n - 1, [&](std::size_t i, std::size_t k, double g)
    {
        if (g >= o.min_gradient)
        {
            ++mark[i];
            --mark[k + 1];
        }
    });

    long run = 0;
    std::size_t q = 0;
    std::vector<std::uint8_t> marked(n, 0);
    for (q = 0; q < n; ++q)
    {
        run += mark[q];
        marked[q] = run > 0 ? 1 : 0;
    }

    auto pair_gradient = [&](std::size_t i) { return dz[i] / length_metres[i]; };
    for (q = 0; q < n;)
    {
        if (!marked[q]) { ++q; continue; }
        std::size_t a = q, b = q;
        while (b + 1 < n && marked[b + 1]) ++b;
        q = b + 1;
        // trim gentle end pairs (a steep window holds at least one steep
        // pair, so the reach does not empty)
        while (a < b && pair_gradient(a) < o.min_gradient) ++a;
        while (b > a && pair_gradient(b) < o.min_gradient) --b;
        lsdtt_chi_step_reach_on_path r;
        r.first_pair = a;
        r.last_pair = b;
        r.n_pairs = static_cast<std::int32_t>(b - a + 1);
        for (std::size_t i = a; i <= b; ++i)
        {
            r.drop_metres += dz[i];
            r.flow_length_metres += length_metres[i];
            if (dz[i] <= 0.0) r.contains_uphill_pairs = true;
        }
        if (r.drop_metres < o.min_drop_metres) continue;
        r.mean_gradient = r.drop_metres / r.flow_length_metres;
        double gmax = -std::numeric_limits<double>::infinity();
        for_each_window(a, b, [&](std::size_t, std::size_t, double g) { gmax = std::max(gmax, g); });
        r.max_window_gradient = gmax;
        out.push_back(r);
    }
    return out;
}

// -----------------------------------------------------------------------------
// the network
// -----------------------------------------------------------------------------
std::vector<lsdtt_chi_step_row> build_chi_step_changes(
    const lsdtt_chi_network_fit_collection_result&  fits,
    const lsdtt_chi_collection_result&              chi_networks,
    const std::vector<lsdtt_chi_basin_source_keys>& keys,
    const lsdtt_raster_meta&                        meta,
    const lsdtt_chi_step_options&                   options)
{
    // validate once, also for an empty network
    (void)detect_chi_step_reaches({}, {}, options);

    std::vector<lsdtt_chi_step_row> out;
    const std::int64_t nc = meta.n_cols;
    const int nr = meta.n_rows;
    const double cs = meta.cell_size;
    auto x_of = [&](std::int64_t f) { return meta.x_min + (static_cast<double>(f % nc) + 0.5) * cs; };
    auto y_of = [&](std::int64_t f) { return meta.y_min + (nr - static_cast<double>(f / nc) - 0.5) * cs; };
    auto d8_length = [&](std::int64_t a, std::int64_t b)
    {
        const std::int64_t dr = std::llabs(a / nc - b / nc), dc = std::llabs(a % nc - b % nc);
        return cs * std::sqrt(static_cast<double>(dr * dr + dc * dc));
    };

    for (std::size_t bi = 0; bi < fits.basins.size(); ++bi)
    {
        const auto& r = fits.basins[bi];
        const std::size_t n = r.flats.size();
        if (n == 0) continue;
        std::unordered_map<std::int64_t, std::size_t> fit_idx;
        fit_idx.reserve(n);
        for (std::size_t i = 0; i < n; ++i) fit_idx[r.flats[i]] = i;

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
            if (!have_keys) return k_chi_step_nan;
            const auto it = key_idx.find(f);
            if (it == key_idx.end()) return k_chi_step_nan;
            const auto& k = keys[bi];
            if (what == 0) return it->second < k.drainage_area_m2.size() ? k.drainage_area_m2[it->second] : k_chi_step_nan;
            if (what == 1) return it->second < k.flow_distance_m.size() ? k.flow_distance_m[it->second] : k_chi_step_nan;
            return it->second < k.source_key.size() ? static_cast<double>(k.source_key[it->second]) : -1.0;
        };
        const std::int32_t basin_key = have_keys ? keys[bi].basin_key : -1;

        // branch heads: the highest-chi node of each branch
        std::map<std::int32_t, std::size_t> head_of;
        for (std::size_t i = 0; i < n; ++i)
        {
            const std::int32_t b = r.branch_id_of_node[i];
            const auto it = head_of.find(b);
            if (it == head_of.end() || r.chi[i] > r.chi[it->second]) head_of[b] = i;
        }

        for (const auto& [b, head] : head_of)
        {
            // walk down from the head; the path ends at the outlet, or at the
            // first node of another branch (the junction), which is included
            std::vector<std::size_t> path;
            std::size_t cur = head;
            for (std::size_t guard = 0; guard <= n; ++guard)
            {
                path.push_back(cur);
                const std::int64_t rec = r.receiver_flats[cur];
                const auto it = fit_idx.find(rec);
                if (rec == r.flats[cur] || it == fit_idx.end()) break;
                if (r.branch_id_of_node[it->second] != b)
                {
                    path.push_back(it->second);
                    break;
                }
                cur = it->second;
            }
            std::reverse(path.begin(), path.end());     // base first
            if (path.size() < 2) continue;
            const bool tributary_base = r.branch_id_of_node[path.front()] != b;

            std::vector<double> dz(path.size() - 1), len(path.size() - 1), below(path.size(), 0.0);
            for (std::size_t j = 0; j + 1 < path.size(); ++j)
            {
                dz[j] = static_cast<double>(r.observed_elevation[path[j + 1]]) - r.observed_elevation[path[j]];
                len[j] = d8_length(r.flats[path[j + 1]], r.flats[path[j]]);
                below[j + 1] = below[j] + len[j];
            }
            const double sk = from_keys(r.flats[head], 2);
            const std::int32_t source_key = std::isfinite(sk) ? static_cast<std::int32_t>(sk) : -1;
            for (const auto& s : detect_chi_step_reaches(dz, len, options))
            {
                const std::size_t lo = path[s.first_pair], hi = path[s.last_pair + 1];
                lsdtt_chi_step_row row;
                row.step_id = static_cast<std::int32_t>(out.size());
                row.basin_index = static_cast<std::int32_t>(bi);
                row.basin_key = basin_key;
                row.source_key = source_key;
                row.branch_id = b;
                row.bottom_flat = r.flats[lo];
                row.top_flat = r.flats[hi];
                row.bottom_x = x_of(row.bottom_flat);
                row.bottom_y = y_of(row.bottom_flat);
                row.top_x = x_of(row.top_flat);
                row.top_y = y_of(row.top_flat);
                row.bottom_elevation = r.observed_elevation[lo];
                row.top_elevation = r.observed_elevation[hi];
                row.drop_m = s.drop_metres;
                row.flow_length_m = s.flow_length_metres;
                row.mean_gradient = s.mean_gradient;
                row.max_window_gradient = s.max_window_gradient;
                row.n_pairs = s.n_pairs;
                row.bottom_chi = r.chi[lo];
                row.top_chi = r.chi[hi];
                row.bottom_flow_distance_m = from_keys(row.bottom_flat, 1);
                row.bottom_drainage_area_m2 = from_keys(row.bottom_flat, 0);
                row.distance_to_branch_base_m = below[s.first_pair];
                const bool near_base = row.distance_to_branch_base_m <= options.max_window_length_metres;
                row.at_tributary_mouth = (near_base && tributary_base) ? 1 : 0;
                row.at_basin_outlet = (near_base && !tributary_base) ? 1 : 0;
                row.contains_uphill_pairs = s.contains_uphill_pairs ? 1 : 0;
                out.push_back(row);
            }
        }
    }
    return out;
}

void write_chi_steps_csv(const std::string& filename, const std::vector<lsdtt_chi_step_row>& rows)
{
    std::ofstream out(filename);
    if (!out) throw std::runtime_error("cannot open " + filename + " for writing");
    out << "step_id,basin_index,basin_key,source_key,branch_id,bottom_flat,top_flat,bottom_x,bottom_y,"
           "top_x,top_y,bottom_elevation,top_elevation,drop_m,flow_length_m,mean_gradient,"
           "max_window_gradient,n_pairs,bottom_chi,top_chi,bottom_flow_distance_m,"
           "bottom_drainage_area_m2,distance_to_branch_base_m,at_tributary_mouth,at_basin_outlet,"
           "contains_uphill_pairs\n";
    for (const auto& s : rows)
        out << s.step_id << ',' << s.basin_index << ',' << s.basin_key << ',' << s.source_key << ','
            << s.branch_id << ',' << s.bottom_flat << ',' << s.top_flat << ',' << num(s.bottom_x) << ','
            << num(s.bottom_y) << ',' << num(s.top_x) << ',' << num(s.top_y) << ','
            << num(s.bottom_elevation) << ',' << num(s.top_elevation) << ',' << num(s.drop_m) << ','
            << num(s.flow_length_m) << ',' << num(s.mean_gradient) << ',' << num(s.max_window_gradient)
            << ',' << s.n_pairs << ',' << num(s.bottom_chi) << ',' << num(s.top_chi) << ','
            << num(s.bottom_flow_distance_m) << ',' << num(s.bottom_drainage_area_m2) << ','
            << num(s.distance_to_branch_base_m) << ',' << int(s.at_tributary_mouth) << ','
            << int(s.at_basin_outlet) << ',' << int(s.contains_uphill_pairs) << '\n';
}

} // namespace lsdtt
