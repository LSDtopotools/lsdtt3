// lsdtt_lithology_relative_erodibility_cpu.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html

#include "lsdtt_lithology_relative_erodibility_cpu.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>

#include "../utils/lsdtt_statistics.hpp"

namespace lsdtt
{

namespace
{

constexpr double k_nan = std::numeric_limits<double>::quiet_NaN();
constexpr std::size_t k_min_reaches_for_interval = 5;   // low_sample below this, as the k_sn summary

std::string csv_text(const std::string& s)
{
    if (s.find_first_of(",\"\n") == std::string::npos) return s;
    std::string q = "\"";
    for (char c : s) { if (c == '"') q += "\"\""; else q += c; }
    return q + "\"";
}

std::string lower(std::string s)
{
    for (char& c : s) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    return s;
}

std::string num(double v, int precision = 10)
{
    std::ostringstream s;
    s << std::setprecision(precision) << v;
    return s.str();
}

// Names that are groups, not units, in the lithology tables of either CLI.
bool is_group_name(const std::string& name)
{
    const std::string l = lower(name);
    return l == "unmapped" || l == "all" || l == "arms disagree" || l.rfind("mixed", 0) == 0;
}

[[noreturn]] void throw_unknown_reference(const std::string& requested,
                                          const std::vector<std::string>& names)
{
    if (is_group_name(requested))
        throw std::invalid_argument("relative_erodibility_reference_unit '" + requested +
                                    "' is not a lithologic unit (unmapped, all, mixed and arms "
                                    "disagree are never a reference)");
    std::string list;
    for (const auto& n : names) list += (list.empty() ? "" : ", ") + n;
    throw std::invalid_argument("relative_erodibility_reference_unit '" + requested +
                                "' is not a unit of the lithology map; the units are: " + list);
}

// 2.5 and 97.5 % of the finite values; NaN with 0 when there are none.
struct interval
{
    double low = k_nan, high = k_nan;
    int n = 0;
};

interval percentile_interval(const std::vector<double>& values)
{
    std::vector<double> v;
    v.reserve(values.size());
    for (double x : values)
        if (std::isfinite(x)) v.push_back(x);
    interval out;
    if (v.empty()) return out;
    std::vector<double> a = v;
    out.low = compute_quantile_interpolated_in_place(a, 0.025);
    a = v;
    out.high = compute_quantile_interpolated_in_place(a, 0.975);
    out.n = static_cast<int>(v.size());
    return out;
}

double ratio_of(double k_ref, double k_i, double n)
{
    if (!(k_ref > 0.0) || !(k_i > 0.0) || !std::isfinite(k_ref) || !std::isfinite(k_i)) return k_nan;
    return std::pow(k_ref / k_i, n);
}

double absolute_k(double E, double A0, double m, double k_i, double n)
{
    return stream_power_erodibility_from_erosion_rate(E, k_i, m, n, A0);
}

} // namespace

double stream_power_erodibility_from_erosion_rate(double E_m_per_yr, double k_sn, double m,
                                                  double n, double A0)
{
    if (!(E_m_per_yr > 0.0) || !(k_sn > 0.0) || !std::isfinite(k_sn)) return k_nan;
    return E_m_per_yr / (std::pow(A0, m) * std::pow(k_sn, n));
}

std::mt19937_64 lithology_relative_erodibility_bootstrap_engine(std::int64_t seed, double m_over_n)
{
    const long long mn = std::llround(m_over_n * 10000.0);
    std::seed_seq seq{static_cast<std::uint32_t>(seed), static_cast<std::uint32_t>(mn), 1u};
    return std::mt19937_64(seq);
}

void check_relative_erodibility_reference_name(const lsdtt_lithology_units_on_grid& units,
                                               const std::string& requested_name)
{
    if (requested_name.empty()) return;
    std::vector<std::string> names;
    for (const auto& u : units.units)
    {
        if (u.name == requested_name) return;
        names.push_back(u.name);
    }
    throw_unknown_reference(requested_name, names);
}

lsdtt_relative_erodibility_reference choose_relative_erodibility_reference_unit(
    const lsdtt_chi_ksn_lithology_summary& summary, const std::string& requested_name)
{
    lsdtt_relative_erodibility_reference ref;
    const lsdtt_chi_ksn_lithology_summary_row* best = nullptr;
    if (!requested_name.empty())
    {
        std::vector<std::string> names;
        for (const auto& r : summary.rows)
        {
            if (r.lith_code < 0) continue;
            names.push_back(r.lith_name);
            if (r.lith_name == requested_name && !best) best = &r;
        }
        if (!best) throw_unknown_reference(requested_name, names);
        if (best->n_reaches_qualifying == 0)
            throw std::invalid_argument("relative_erodibility_reference_unit '" + requested_name +
                                        "' has no qualifying lithologic reach, so it has no k_sn");
        ref.reason = "named by relative_erodibility_reference_unit; " +
                     std::to_string(best->n_reaches_qualifying) + " qualifying reaches";
    }
    else
    {
        // Most qualifying reaches, then the larger channel length, then the
        // lower code (rows are in ascending code, so the first wins a full tie).
        std::size_t n_tied = 0;
        for (const auto& r : summary.rows)
        {
            if (r.lith_code < 0 || r.n_reaches_qualifying == 0) continue;
            if (!best || r.n_reaches_qualifying > best->n_reaches_qualifying ||
                (r.n_reaches_qualifying == best->n_reaches_qualifying &&
                 r.channel_length_m > best->channel_length_m))
                best = &r;
        }
        if (!best)
            throw std::invalid_argument("no lithologic unit has a qualifying reach, so there is no "
                                        "reference for relative erodibility");
        bool length_tie = false;
        for (const auto& r : summary.rows)
        {
            if (r.lith_code < 0 || &r == best) continue;
            if (r.n_reaches_qualifying == best->n_reaches_qualifying)
            {
                ++n_tied;
                if (r.channel_length_m == best->channel_length_m) length_tie = true;
            }
        }
        std::ostringstream why;
        why << "default: the most qualifying reaches (" << best->n_reaches_qualifying << ")";
        if (n_tied > 0)
        {
            why << ", tied with " << n_tied << " other unit(s); ";
            why << (length_tie ? "equal channel length too, so the lower code"
                               : "the larger channel length (" + num(best->channel_length_m, 8) + " m)");
        }
        ref.reason = why.str();
    }
    ref.code = best->lith_code;
    ref.name = best->lith_name;
    return ref;
}

lsdtt_relative_erodibility_result compute_lithology_relative_erodibility(
    const lsdtt_chi_lithology_tables& tables, const lsdtt_chi_ksn_lithology_summary& summary,
    const lsdtt_relative_erodibility_options& opts, std::mt19937_64& rng)
{
    if (opts.n_values.empty())
        throw std::invalid_argument("relative_erodibility_n_values is empty");
    for (double n : opts.n_values)
        if (!std::isfinite(n) || !(n > 0.0))
            throw std::invalid_argument("relative_erodibility_n_values must all be finite and > 0");
    if (!std::isfinite(opts.erosion_rate_m_per_yr) || opts.erosion_rate_m_per_yr < 0.0)
        throw std::invalid_argument("relative_erodibility_erosion_rate_m_per_yr must be >= 0 "
                                    "(0 = not given)");

    lsdtt_relative_erodibility_result res;
    res.m_over_n = summary.m_over_n;
    res.A_0 = summary.A_0;
    res.bootstrap_cluster = opts.bootstrap_cluster;
    res.erosion_rate_m_per_yr = opts.erosion_rate_m_per_yr;
    res.n_values = opts.n_values;
    res.reference = choose_relative_erodibility_reference_unit(summary, opts.reference_unit);
    const bool with_E = opts.erosion_rate_m_per_yr > 0.0;

    // The units (code >= 0) in ascending code, with their qualifying reaches.
    struct unit_data
    {
        const lsdtt_chi_ksn_lithology_summary_row* row = nullptr;
        std::vector<double> k;
        std::vector<std::int32_t> basin;
        bool low_sample = true;
        std::vector<double> replicate;   // replicate medians, NaN = none
    };
    std::vector<unit_data> units;
    std::map<std::int32_t, std::size_t> index_of_code;
    for (const auto& r : summary.rows)
    {
        if (r.lith_code < 0) continue;
        index_of_code[r.lith_code] = units.size();
        unit_data u;
        u.row = &r;
        units.push_back(std::move(u));
    }
    for (const auto& r : tables.reaches)
    {
        if (!r.qualifies) continue;
        const auto it = index_of_code.find(r.lith_code);
        if (it == index_of_code.end()) continue;
        units[it->second].k.push_back(r.k_sn_ols);
        units[it->second].basin.push_back(r.basin_index);
    }
    for (auto& u : units) u.low_sample = u.k.size() < k_min_reaches_for_interval;
    unit_data& ref = units[index_of_code.at(res.reference.code)];

    // Bootstrap: one set of replicate medians per unit, shared by every n.
    const int n_rep = std::max(0, opts.bootstrap_replicates);
    if (n_rep > 0)
    {
        for (auto& u : units)
            if (!u.low_sample) u.replicate.assign(static_cast<std::size_t>(n_rep), k_nan);
        std::vector<double> pooled;
        if (opts.bootstrap_cluster == lsdtt_lithology_bootstrap_cluster::reach)
        {
            // The reference first, then the other units in ascending code.
            std::vector<unit_data*> order;
            if (!ref.low_sample) order.push_back(&ref);
            for (auto& u : units)
                if (&u != &ref && !u.low_sample) order.push_back(&u);
            for (int rep = 0; rep < n_rep; ++rep)
                for (unit_data* u : order)
                {
                    std::uniform_int_distribution<std::size_t> pick(0, u->k.size() - 1);
                    pooled.clear();
                    for (std::size_t j = 0; j < u->k.size(); ++j) pooled.push_back(u->k[pick(rng)]);
                    u->replicate[static_cast<std::size_t>(rep)] = compute_median_in_place(pooled);
                }
        }
        else
        {
            // Basins with a qualifying reach of any unit; a draw brings every
            // unit's reaches in that basin.
            std::set<std::int32_t> basin_set;
            for (const auto& u : units) basin_set.insert(u.basin.begin(), u.basin.end());
            const std::vector<std::int32_t> basins(basin_set.begin(), basin_set.end());
            if (!basins.empty())
            {
                std::map<std::int32_t, std::size_t> slot;
                for (std::size_t b = 0; b < basins.size(); ++b) slot[basins[b]] = b;
                // reaches_in[u][b]: unit u's reach k_sn in basin slot b.
                std::vector<std::vector<std::vector<double>>> reaches_in(units.size());
                for (std::size_t ui = 0; ui < units.size(); ++ui)
                {
                    reaches_in[ui].resize(basins.size());
                    for (std::size_t j = 0; j < units[ui].k.size(); ++j)
                        reaches_in[ui][slot[units[ui].basin[j]]].push_back(units[ui].k[j]);
                }
                std::uniform_int_distribution<std::size_t> pick(0, basins.size() - 1);
                std::vector<std::size_t> drawn(basins.size());
                for (int rep = 0; rep < n_rep; ++rep)
                {
                    for (auto& d : drawn) d = pick(rng);
                    for (std::size_t ui = 0; ui < units.size(); ++ui)
                    {
                        if (units[ui].low_sample) continue;
                        pooled.clear();
                        for (std::size_t d : drawn)
                            pooled.insert(pooled.end(), reaches_in[ui][d].begin(), reaches_in[ui][d].end());
                        if (!pooled.empty())
                            units[ui].replicate[static_cast<std::size_t>(rep)] = compute_median_in_place(pooled);
                    }
                }
            }
        }
    }

    const double k_ref = ref.row->reach_k_sn_median;
    for (double n : opts.n_values)
    {
        const double m = res.m_over_n * n;
        for (const auto& u : units)
        {
            lsdtt_relative_erodibility_row row;
            row.lith_code = u.row->lith_code;
            row.lith_name = u.row->lith_name;
            row.n = n;
            row.m = m;
            row.n_reaches_qualifying = u.row->n_reaches_qualifying;
            row.reach_k_sn_median = u.row->reach_k_sn_median;
            row.low_sample = u.low_sample;
            row.is_reference = &u == &ref;
            row.k_ratio = row.is_reference && std::isfinite(k_ref) ? 1.0 : ratio_of(k_ref, row.reach_k_sn_median, n);
            if (with_E) row.K_abs = absolute_k(opts.erosion_rate_m_per_yr, res.A_0, m, row.reach_k_sn_median, n);

            if (!u.replicate.empty())
            {
                if (!ref.replicate.empty())
                {
                    std::vector<double> ratios(u.replicate.size());
                    for (std::size_t b = 0; b < ratios.size(); ++b)
                        ratios[b] = row.is_reference
                                        ? (ref.replicate[b] > 0.0 ? 1.0 : k_nan)
                                        : ratio_of(ref.replicate[b], u.replicate[b], n);
                    const interval ci = percentile_interval(ratios);
                    row.k_ratio_ci_low = ci.low;
                    row.k_ratio_ci_high = ci.high;
                    row.n_bootstrap = ci.n;
                }
                if (with_E)
                {
                    std::vector<double> ks(u.replicate.size());
                    for (std::size_t b = 0; b < ks.size(); ++b)
                        ks[b] = absolute_k(opts.erosion_rate_m_per_yr, res.A_0, m, u.replicate[b], n);
                    const interval ci = percentile_interval(ks);
                    row.K_abs_ci_low = ci.low;
                    row.K_abs_ci_high = ci.high;
                }
            }
            res.rows.push_back(row);
        }
    }
    return res;
}

std::string lithology_erodibility_units_string(double m)
{
    double e = 1.0 - 2.0 * m;
    if (std::fabs(e) < 1e-12) e = 0.0;
    std::ostringstream s;
    s << "m^(" << std::setprecision(6) << e << ") yr^-1";
    return s.str();
}

std::string relative_erodibility_n_suffix(double n)
{
    std::ostringstream s;
    s << std::fixed << std::setprecision(3) << n;
    std::string t = s.str();
    while (!t.empty() && t.back() == '0') t.pop_back();
    if (!t.empty() && t.back() == '.') t.pop_back();
    for (char& c : t) if (c == '.') c = 'p';
    return t;
}

void write_lithology_relative_erodibility_csv(const std::string& path,
                                              const lsdtt_relative_erodibility_result& r)
{
    std::ofstream out(path);
    if (!out) throw std::runtime_error("cannot write " + path);
    const bool with_E = r.erosion_rate_m_per_yr > 0.0;
    out << "lith_code,lith_name,m_over_n,A_0,n,m,reference_unit,n_reaches_qualifying,"
           "reach_k_sn_median,k_ratio,k_ratio_ci_low,k_ratio_ci_high,n_bootstrap,"
           "bootstrap_cluster,low_sample";
    if (with_E) out << ",erosion_rate_m_per_yr,K_abs,K_abs_ci_low,K_abs_ci_high,K_units";
    out << "\n" << std::setprecision(10);
    const std::string cluster = lithology_bootstrap_cluster_name(r.bootstrap_cluster);
    const std::string ref = csv_text(r.reference.name);
    for (const auto& w : r.rows)
    {
        out << w.lith_code << "," << csv_text(w.lith_name) << "," << r.m_over_n << "," << r.A_0 << ","
            << w.n << "," << w.m << "," << ref << "," << w.n_reaches_qualifying << ","
            << w.reach_k_sn_median << "," << w.k_ratio << "," << w.k_ratio_ci_low << ","
            << w.k_ratio_ci_high << "," << w.n_bootstrap << "," << cluster << ","
            << (w.low_sample ? 1 : 0);
        if (with_E)
            out << "," << r.erosion_rate_m_per_yr << "," << w.K_abs << "," << w.K_abs_ci_low << ","
                << w.K_abs_ci_high << "," << lithology_erodibility_units_string(w.m);
        out << "\n";
    }
}

std::size_t write_lithology_erodibility_property_table(
    const std::string& path, double m, double n, const std::string& provenance,
    const std::vector<std::string>& extra_comment_lines,
    const std::vector<lsdtt_lithology_erodibility_table_row>& units,
    const std::string& not_estimated_note)
{
    std::vector<const lsdtt_lithology_erodibility_table_row*> rows, missing;
    for (const auto& u : units)
    {
        if (std::isfinite(u.K) && u.K > 0.0 && u.code >= 0 && u.code <= 65534)
            rows.push_back(&u);
        else
            missing.push_back(&u);
    }
    std::map<std::string, std::string> seen;
    for (const auto* w : rows)
    {
        const auto it = seen.find(lower(w->name));
        if (it != seen.end())
            throw std::invalid_argument("units '" + it->second + "' and '" + w->name +
                                        "' have the same name ignoring case; the property table "
                                        "reader refuses that (name the units apart in a crosswalk)");
        seen[lower(w->name)] = w->name;
    }

    std::ofstream out(path);
    if (!out) throw std::runtime_error("cannot write " + path);
    out << "# lsdtt lithology property table v1\n"
        << "# erodibility_convention: map_area\n"
        << "# m_exponent: " << num(m) << "\n"
        << "# n_exponent: " << num(n) << "\n"
        << "# provenance: " << provenance << "\n";
    for (const auto& line : extra_comment_lines) out << "# " << line << "\n";
    if (!missing.empty())
    {
        out << "# not_estimated:";
        for (std::size_t i = 0; i < missing.size(); ++i)
            out << (i == 0 ? " " : "; ") << missing[i]->code << " " << missing[i]->name;
        out << " (" << not_estimated_note << ")\n";
    }
    out << "code,name,erodibility,reference\n";
    for (const auto* w : rows)
        out << w->code << "," << csv_text(w->name) << "," << num(w->K) << ","
            << csv_text(w->reference) << "\n";
    return rows.size();
}

std::size_t write_lithology_relative_erodibility_property_table(
    const std::string& path, const lsdtt_relative_erodibility_result& r, double n,
    const lsdtt_relative_erodibility_provenance& p)
{
    if (!(r.erosion_rate_m_per_yr > 0.0))
        throw std::invalid_argument("a draft property table needs an erosion rate > 0");
    const double m = r.m_over_n * n;

    std::ostringstream prov;
    prov << "draft from lsdtt-chi-analysis relative erodibility; DEM " << p.dem
         << "; lithology map " << p.lithology_map;
    if (!p.field.empty()) prov << " (field " << p.field << ")";
    if (!p.crosswalk.empty()) prov << " with crosswalk " << p.crosswalk;
    prov << "; m/n " << num(r.m_over_n) << ", A0 " << num(r.A_0) << " m^2, n " << num(n) << ", m "
         << num(m) << "; erosion rate " << num(r.erosion_rate_m_per_yr)
         << " m/yr (given, applied everywhere); reference unit " << r.reference.name
         << "; K = E / (A0^m k_sn^n) from each unit's reach-median k_sn; assumes uniform erosion and "
            "steady state";
    const std::vector<std::string> extra{
        "K units: " + lithology_erodibility_units_string(m) + " (map area in m^2)",
        "bootstrap: " + lithology_bootstrap_cluster_name(r.bootstrap_cluster) +
            " clusters; the interval in each reference column is the 95 % percentile interval of K"};

    std::vector<lsdtt_lithology_erodibility_table_row> units;
    for (const auto& w : r.rows)
    {
        if (w.n != n) continue;
        lsdtt_lithology_erodibility_table_row u;
        u.code = w.lith_code;
        u.name = w.lith_name;
        u.K = w.K_abs;
        std::ostringstream ref;
        ref << "reach median k_sn " << num(w.reach_k_sn_median, 8) << " from "
            << w.n_reaches_qualifying << " qualifying reaches; K/K_ref " << num(w.k_ratio, 8);
        if (w.low_sample)
            ref << "; low_sample (fewer than 5 qualifying reaches): no interval";
        else
            ref << "; K 95 % interval [" << num(w.K_abs_ci_low, 8) << ", " << num(w.K_abs_ci_high, 8)
                << "]";
        u.reference = ref.str();
        units.push_back(std::move(u));
    }
    return write_lithology_erodibility_property_table(
        path, m, n, prov.str(), extra, units,
        "no qualifying reach or k_sn <= 0; add rows by hand before mapping a region that "
        "contains them");
}

} // namespace lsdtt
