// lsdtt_surface_training_cpu.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html

#include "lsdtt_surface_training_cpu.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <map>
#include <stdexcept>
#include <unordered_map>

#include "lsdtt_vector_rasterise_cpu.hpp"

namespace lsdtt
{

namespace
{

// Linear-interpolation quantile of a sorted, non-empty sample.
double quantile_sorted(const std::vector<double>& v, double q)
{
    const double pos = q * static_cast<double>(v.size() - 1);
    const std::size_t lo = static_cast<std::size_t>(std::floor(pos));
    const std::size_t hi = std::min(lo + 1, v.size() - 1);
    return v[lo] + (pos - static_cast<double>(lo)) * (v[hi] - v[lo]);
}

bool all_digits(const std::string& s)
{
    return std::all_of(s.begin(), s.end(),
                       [](unsigned char c) { return std::isdigit(c) != 0; });
}

} // namespace

int parse_training_class(const std::string& label, int& number)
{
    std::string s;
    for (const unsigned char c : label)
        if (!std::isspace(c) && c != '_' && c != '-')
            s.push_back(static_cast<char>(std::tolower(c)));
    number = 0;
    if (s == "floodplain" || s == "fp") return 0;
    std::string digits;
    if (s.rfind("terrace", 0) == 0)     digits = s.substr(7);
    else if (s.size() > 1 && s[0] == 't') digits = s.substr(1);
    else return -1;
    if (!all_digits(digits)) return -1;
    if (!digits.empty()) number = std::stoi(digits);
    return 1;
}

training_set read_training_surfaces(const lsdtt_raster_cpu& reference,
                                    const std::string&      vector_path,
                                    const std::string&      layer_name,
                                    const std::string&      class_field)
{
    const polygon_field_rasterise_result r = rasterise_polygon_field_values(
        reference, vector_path, layer_name, class_field, 4,
        /*fill_coverage_gaps=*/false);

    training_set set;
    set.polygon_of_pixel = r.feature_index;
    set.polygons.resize(r.feature_values.size());
    std::vector<std::string> bad;
    for (std::size_t k = 0; k < r.feature_values.size(); ++k)
    {
        training_polygon& p = set.polygons[k];
        p.label = r.feature_values[k];
        p.cls   = parse_training_class(p.label, p.number);
        if (p.cls < 0 &&
            std::find(bad.begin(), bad.end(), p.label) == bad.end())
            bad.push_back(p.label);
    }
    if (!bad.empty())
    {
        std::string msg = "training polygons: class field '" + class_field +
                          "' has values that are neither floodplain (floodplain, "
                          "fp) nor terrace (terraceN, tN):";
        for (const std::string& b : bad) msg += " '" + b + "'";
        throw std::invalid_argument(msg);
    }
    for (const std::int32_t f : set.polygon_of_pixel)
        if (f >= 0) ++set.polygons[static_cast<std::size_t>(f)].n_cells;
    return set;
}

void compute_training_polygon_stats(training_set&           set,
                                    const lsdtt_raster_cpu& field,
                                    double                  spread_limit)
{
    const std::vector<float>&   fd = field.get_data();
    const std::vector<uint8_t>& fm = field.get_mask();
    if (fd.size() != set.polygon_of_pixel.size())
        throw std::invalid_argument(
            "compute_training_polygon_stats: field is not on the training grid");
    std::vector<std::vector<double>> vals(set.polygons.size());
    for (std::size_t i = 0; i < fd.size(); ++i)
    {
        const std::int32_t f = set.polygon_of_pixel[i];
        if (f < 0 || !fm[i] || !std::isfinite(fd[i])) continue;
        vals[static_cast<std::size_t>(f)].push_back(static_cast<double>(fd[i]));
    }
    for (std::size_t k = 0; k < set.polygons.size(); ++k)
    {
        training_polygon& p = set.polygons[k];
        std::vector<double>& v = vals[k];
        p.n_with_field = static_cast<std::int64_t>(v.size());
        if (v.empty()) continue;
        std::sort(v.begin(), v.end());
        p.median  = quantile_sorted(v, 0.5);
        p.p10     = quantile_sorted(v, 0.1);
        p.p90     = quantile_sorted(v, 0.9);
        p.flagged = spread_limit > 0.0 && (p.p90 - p.p10) > spread_limit;
    }
}

training_split best_separation(std::vector<double> lower,
                               std::vector<double> upper)
{
    if (lower.empty() || upper.empty())
        throw std::invalid_argument("best_separation: an empty sample");
    std::sort(lower.begin(), lower.end());
    std::sort(upper.begin(), upper.end());
    const double nl = static_cast<double>(lower.size());
    const double nu = static_cast<double>(upper.size());

    // Distinct candidate values, ascending. J is constant on [v_j, v_{j+1}).
    std::vector<double> cand;
    cand.reserve(lower.size() + upper.size());
    std::merge(lower.begin(), lower.end(), upper.begin(), upper.end(),
               std::back_inserter(cand));
    cand.erase(std::unique(cand.begin(), cand.end()), cand.end());

    std::vector<double> j(cand.size());
    std::size_t il = 0, iu = 0;
    for (std::size_t c = 0; c < cand.size(); ++c)
    {
        while (il < lower.size() && lower[il] <= cand[c]) ++il;
        while (iu < upper.size() && upper[iu] <= cand[c]) ++iu;
        j[c] = static_cast<double>(il) / nl +
               static_cast<double>(upper.size() - iu) / nu - 1.0;
    }
    const double jmax = *std::max_element(j.begin(), j.end());
    // The first run of maximal candidates; the threshold is the midpoint of
    // the interval they cover, [cand[a], cand[b + 1]).
    std::size_t a = 0;
    while (j[a] < jmax - 1e-12) ++a;
    std::size_t b = a;
    while (b + 1 < cand.size() && j[b + 1] >= jmax - 1e-12) ++b;
    const double hi = (b + 1 < cand.size()) ? cand[b + 1] : cand[b];

    training_split s;
    s.threshold = 0.5 * (cand[a] + hi);
    s.youden_j  = jmax;
    s.n_lower   = static_cast<std::int64_t>(lower.size());
    s.n_upper   = static_cast<std::int64_t>(upper.size());
    s.lower_below = std::upper_bound(lower.begin(), lower.end(), s.threshold) -
                    lower.begin();
    s.upper_above = upper.end() -
                    std::upper_bound(upper.begin(), upper.end(), s.threshold);
    return s;
}

training_window_fit fit_training_window(
    const std::vector<std::pair<std::int32_t, double>>& cells,
    const training_set&                                 set,
    double                                              required_gap,
    int                                                 min_cells_per_class)
{
    training_window_fit fit;
    const std::size_t min_n =
        static_cast<std::size_t>(std::max(1, min_cells_per_class));

    // Trimmed cores, per polygon, of the cells in this window.
    std::map<std::int32_t, std::vector<double>> core;
    for (const auto& [f, v] : cells)
    {
        if (f < 0) continue;
        const training_polygon& p = set.polygons[static_cast<std::size_t>(f)];
        if (p.n_with_field == 0 || v < p.p10 || v > p.p90) continue;
        core[f].push_back(v);
    }

    std::vector<double> fp_cells, terrace_cells;
    std::vector<std::int32_t> fp_ids;
    std::vector<std::int32_t> level_ids;   // single-level terrace polygons
    for (const auto& [f, v] : core)
    {
        const training_polygon& p = set.polygons[static_cast<std::size_t>(f)];
        if (p.cls == 0)
        {
            fp_cells.insert(fp_cells.end(), v.begin(), v.end());
            fp_ids.push_back(f);
        }
        else
        {
            terrace_cells.insert(terrace_cells.end(), v.begin(), v.end());
            if (!p.flagged) level_ids.push_back(f);
        }
    }

    auto add_note = [&](const std::string& s)
    {
        if (!fit.note.empty()) fit.note += "; ";
        fit.note += s;
    };

    const bool fp_ok = fp_cells.size() >= min_n;
    if (fp_ok && terrace_cells.size() >= min_n)
    {
        fit.has_bound = true;
        fit.bound = best_separation(fp_cells, terrace_cells);
    }
    else if (!fp_cells.empty() || !terrace_cells.empty())
        add_note("bound: fewer than " + std::to_string(min_n) +
                 " core cells of floodplain or terrace");

    // Levels: floodplain first, then terrace polygons grouped by median.
    std::vector<std::vector<std::int32_t>> groups;
    if (fp_ok) groups.push_back(fp_ids);
    std::sort(level_ids.begin(), level_ids.end(),
              [&](std::int32_t x, std::int32_t y)
              {
                  return set.polygons[static_cast<std::size_t>(x)].median <
                         set.polygons[static_cast<std::size_t>(y)].median;
              });
    std::vector<std::vector<std::int32_t>> tgroups;
    double last_median = 0.0;
    for (const std::int32_t f : level_ids)
    {
        const double m = set.polygons[static_cast<std::size_t>(f)].median;
        if (tgroups.empty() || m - last_median >= required_gap)
            tgroups.emplace_back();
        tgroups.back().push_back(f);
        last_median = m;
    }
    auto group_cells = [&](const std::vector<std::int32_t>& g)
    {
        std::vector<double> out;
        for (const std::int32_t f : g)
            out.insert(out.end(), core[f].begin(), core[f].end());
        return out;
    };
    for (auto& g : tgroups)
    {
        if (group_cells(g).size() >= min_n) groups.push_back(std::move(g));
        else add_note("a terrace level with fewer than " +
                      std::to_string(min_n) + " core cells dropped");
    }

    // Breaks between consecutive levels; a break not above the one below it
    // means the two upper populations do not separate in this field, so the
    // upper level is merged into the one below.
    std::vector<double> lower = groups.empty() ? std::vector<double>{}
                                               : group_cells(groups[0]);
    if (!groups.empty()) fit.levels.push_back(groups[0]);
    for (std::size_t g = 1; g < groups.size(); ++g)
    {
        const std::vector<double> upper = group_cells(groups[g]);
        const training_split s = best_separation(lower, upper);
        if (!fit.breaks.empty() && !(s.threshold > fit.breaks.back()))
        {
            add_note("level " + std::to_string(g) +
                     " does not separate from the level below; merged");
            fit.levels.back().insert(fit.levels.back().end(),
                                     groups[g].begin(), groups[g].end());
            lower.insert(lower.end(), upper.begin(), upper.end());
            continue;
        }
        fit.breaks.push_back(s.threshold);
        fit.break_fits.push_back(s);
        fit.levels.push_back(groups[g]);
        lower = upper;
    }
    return fit;
}

donor_reach_coordinate donor_reach_coordinates(
    const std::vector<std::int64_t>& donor_flat,
    const lsdtt_raster_cpu&          channel_flow_distance,
    const std::vector<std::int32_t>& channel_source_id)
{
    const std::vector<float>&   fd  = channel_flow_distance.get_data();
    const std::vector<uint8_t>& fdm = channel_flow_distance.get_mask();
    const std::size_t n = fd.size();
    if (donor_flat.size() != n || channel_source_id.size() != n)
        throw std::invalid_argument("donor_reach_coordinates: inputs differ in size");

    donor_reach_coordinate out;
    out.source_of_pixel.assign(n, -1);
    std::vector<float>   dd(n, 0.0f);
    std::vector<uint8_t> dm(n, 0);
    for (std::size_t i = 0; i < n; ++i)
    {
        const std::int64_t df = donor_flat[i];
        if (df < 0) continue;
        const std::size_t dfi = static_cast<std::size_t>(df);
        if (dfi >= n || !fdm[dfi] || !std::isfinite(fd[dfi])) continue;
        out.source_of_pixel[i] = channel_source_id[dfi];
        dd[i] = fd[dfi];
        dm[i] = 1;
    }
    out.distance = lsdtt_raster_cpu(std::move(dd), std::move(dm),
                                    channel_flow_distance.get_meta());
    return out;
}

reach_window_set assign_reach_windows(
    const std::vector<std::int32_t>& source_of_pixel,
    const lsdtt_raster_cpu&          coordinate,
    const std::vector<std::uint8_t>& defined,
    double                           reach_window_metres)
{
    const std::vector<float>&   cd = coordinate.get_data();
    const std::vector<uint8_t>& cm = coordinate.get_mask();
    const std::size_t n = cd.size();
    if (source_of_pixel.size() != n || defined.size() != n)
        throw std::invalid_argument("assign_reach_windows: inputs differ in size");
    const double cell = coordinate.get_meta().cell_size;

    auto usable = [&](std::size_t i)
    {
        return defined[i] && cm[i] && std::isfinite(cd[i]) &&
               source_of_pixel[i] >= 0;
    };

    // The ladder's windowing (build_reference_field_ladder), source by source.
    std::map<std::int32_t, std::pair<double, double>> drange;
    for (std::size_t i = 0; i < n; ++i)
    {
        if (!usable(i)) continue;
        const double d = static_cast<double>(cd[i]);
        auto it = drange.find(source_of_pixel[i]);
        if (it == drange.end()) drange[source_of_pixel[i]] = {d, d};
        else
        {
            it->second.first  = std::min(it->second.first, d);
            it->second.second = std::max(it->second.second, d);
        }
    }

    struct bins { double dmin, wstep; int n_windows; std::int32_t first; };
    std::unordered_map<std::int32_t, bins> per_source;
    reach_window_set out;
    for (const auto& [s, r] : drange)
    {
        const double span = std::max(r.second - r.first, 1.0);
        double w = reach_window_metres;
        if (w <= 0.0) w = std::max(2000.0, span / 8.0);
        bins b{r.first, std::max(w * 0.5, cell),
               0, static_cast<std::int32_t>(out.windows.size())};
        b.n_windows = static_cast<int>(span / b.wstep) + 1;
        for (int k = 0; k < b.n_windows; ++k)
        {
            reach_window_set::window win;
            win.source   = s;
            win.index    = k;
            win.centre_m = b.dmin + (k + 0.5) * b.wstep;
            for (int q = std::max(k - 1, 0); q <= std::min(k + 1, b.n_windows - 1); ++q)
                win.neighbours.push_back(b.first + q);
            out.windows.push_back(std::move(win));
        }
        per_source.emplace(s, b);
    }

    out.window_of_pixel.assign(n, -1);
    for (std::size_t i = 0; i < n; ++i)
    {
        if (!usable(i)) continue;
        const bins& b = per_source.at(source_of_pixel[i]);
        const int w = std::clamp(
            static_cast<int>((static_cast<double>(cd[i]) - b.dmin) / b.wstep),
            0, b.n_windows - 1);
        out.window_of_pixel[i] = b.first + w;
    }
    return out;
}

} // namespace lsdtt
