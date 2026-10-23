// lsdtt_erodibility_calibration_cpu.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html

#include "lsdtt_erodibility_calibration_cpu.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <map>
#include <numeric>
#include <set>
#include <sstream>
#include <stdexcept>
#include <unordered_map>

#include "../core/lsdtt_data_loader.hpp"
#include "../utils/lsdtt_statistics.hpp"
#include "lsdtt_basin_cpu.hpp"
#include "lsdtt_basin_outline_cpu.hpp"
#include "lsdtt_flow_snapping_cpu.hpp"
#include "lsdtt_lithology_relative_erodibility_cpu.hpp"

namespace lsdtt
{

namespace
{

std::string lower_trimmed(const std::string& s)
{
    std::size_t a = 0, b = s.size();
    while (a < b && std::isspace(static_cast<unsigned char>(s[a]))) ++a;
    while (b > a && std::isspace(static_cast<unsigned char>(s[b - 1]))) --b;
    std::string out = s.substr(a, b - a);
    for (char& c : out) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    return out;
}

// The attribute of a record whose key equals `column` ignoring case; nullptr
// when there is none.
const lsdtt_point_attribute_value* find_attribute(const lsdtt_input_point_with_attributes& rec,
                                                  const std::string& column)
{
    if (column.empty()) return nullptr;
    const auto exact = rec.attributes.find(column);
    if (exact != rec.attributes.end()) return &exact->second;
    const std::string want = lower_trimmed(column);
    for (const auto& [k, v] : rec.attributes)
        if (lower_trimmed(k) == want) return &v;
    return nullptr;
}

bool has_column(const std::vector<std::string>& columns, const std::string& column)
{
    if (column.empty()) return false;
    const std::string want = lower_trimmed(column);
    for (const auto& c : columns)
        if (lower_trimmed(c) == want) return true;
    return false;
}

std::string join_columns(const std::vector<std::string>& columns)
{
    std::string out;
    for (const auto& c : columns) out += (out.empty() ? "" : ", ") + c;
    return out.empty() ? "(none)" : out;
}

double numeric_or_nan(const lsdtt_point_attribute_value* v)
{
    double d = std::numeric_limits<double>::quiet_NaN();
    if (v == nullptr || !v->is_present()) return d;
    if (!v->as_double(d)) return std::numeric_limits<double>::quiet_NaN();
    return d;
}

} // namespace

std::string accepted_erosion_rate_units_list()
{
    return "mm_per_kyr, m_per_myr, cm_per_kyr, mm_per_yr, m_per_yr";
}

std::string canonical_erosion_rate_units(const std::string& units)
{
    const std::string u = lower_trimmed(units);
    if (u == "mm_per_kyr" || u == "m_per_myr") return "mm_per_kyr";
    if (u == "cm_per_kyr") return "cm_per_kyr";
    if (u == "mm_per_yr") return "mm_per_yr";
    if (u == "m_per_yr") return "m_per_yr";
    const bool mass_flux = u.find("g_per") != std::string::npos ||
                           u.find("g/cm") != std::string::npos ||
                           u.find("g cm") != std::string::npos ||
                           u.find("gcm") != std::string::npos;
    if (mass_flux)
        throw std::invalid_argument("erosion_rate_units '" + units +
                                    "' is a mass flux; convert it to a rate first by dividing by "
                                    "the rock density (e.g. 2.65 g cm^-3), then give one of " +
                                    accepted_erosion_rate_units_list());
    throw std::invalid_argument("erosion_rate_units must be one of " +
                                accepted_erosion_rate_units_list() + " (got '" + units + "')");
}

double erosion_rate_units_to_m_per_yr(const std::string& units)
{
    const std::string c = canonical_erosion_rate_units(units);
    if (c == "mm_per_kyr") return 1e-6;
    if (c == "cm_per_kyr") return 1e-5;
    if (c == "mm_per_yr") return 1e-3;
    return 1.0;
}

lsdtt_erosion_rate_read_result read_erosion_rate_samples(const lsdtt_raster_meta& meta,
                                                         const lsdtt_erosion_rate_input_options& o)
{
    lsdtt_erosion_rate_read_result res;
    res.units = canonical_erosion_rate_units(o.units);
    res.units_to_m_per_yr = erosion_rate_units_to_m_per_yr(o.units);
    if (o.path.empty() || !std::filesystem::exists(o.path))
        throw std::invalid_argument("erosion_rate_fname does not exist: " + o.path);

    lsdtt_point_input_options po;
    po.vector_layer_name = o.layer_name;
    po.csv_latitude_field_name = o.latitude_column;
    po.csv_longitude_field_name = o.longitude_column;
    po.keep_points_outside_raster = true;
    const lsdtt_point_attribute_read_result read = load_point_layer_with_attributes(o.path, meta, {}, po);
    res.columns = read.fields_found;

    if (!has_column(res.columns, o.rate_column))
        throw std::invalid_argument("erosion_rate_column '" + o.rate_column + "' is not a column of " +
                                    o.path + "; the columns are: " + join_columns(res.columns));
    if (!has_column(res.columns, o.uncertainty_column))
        throw std::invalid_argument("erosion_rate_uncertainty_column '" + o.uncertainty_column +
                                    "' is not a column of " + o.path + "; the columns are: " +
                                    join_columns(res.columns));
    res.has_sample_id = has_column(res.columns, o.sample_id_column);
    res.has_reported_area = has_column(res.columns, o.reported_area_column);
    res.has_units_column = has_column(res.columns, "erosion_rate_units");

    int row = 0;
    for (const auto& rec : read.records)
    {
        ++row;
        lsdtt_erosion_rate_sample s;
        s.input_row = row;
        s.point = rec.point;

        std::string id;
        const lsdtt_point_attribute_value* idv = find_attribute(rec, o.sample_id_column);
        if (idv && idv->is_present()) (void)idv->as_source_string(id);
        s.sample_id = lower_trimmed(id).empty() ? "row_" + std::to_string(row) : id;

        if (res.has_units_column)
        {
            std::string u;
            const lsdtt_point_attribute_value* uv = find_attribute(rec, "erosion_rate_units");
            if (uv && uv->is_present()) (void)uv->as_string(u);
            std::string canon;
            try { canon = canonical_erosion_rate_units(u); }
            catch (const std::invalid_argument&) { canon.clear(); }
            if (canon != res.units)
                throw std::invalid_argument(
                    "row " + std::to_string(row) + " (sample " + s.sample_id +
                    ") of " + o.path + ": its erosion_rate_units column says '" + u +
                    "', which disagrees with erosion_rate_units: " + o.units);
        }

        s.rate_input = numeric_or_nan(find_attribute(rec, o.rate_column));
        s.uncertainty_input = numeric_or_nan(find_attribute(rec, o.uncertainty_column));
        if (std::isfinite(s.rate_input)) s.E_m_per_yr = s.rate_input * res.units_to_m_per_yr;
        if (std::isfinite(s.uncertainty_input))
            s.sigma_E_m_per_yr = s.uncertainty_input * res.units_to_m_per_yr;
        if (res.has_reported_area)
        {
            const double a = numeric_or_nan(find_attribute(rec, o.reported_area_column));
            if (std::isfinite(a) && a > 0.0) s.reported_area_m2 = a * 1.0e6;
        }

        if (!std::isfinite(s.rate_input) || !(s.rate_input > 0.0))
        {
            s.status = "bad_rate";
            s.status_detail = "erosion rate missing, not a number or <= 0";
        }
        else if (!std::isfinite(s.uncertainty_input) || !(s.uncertainty_input > 0.0))
        {
            s.status = "missing_uncertainty";
            s.status_detail = "uncertainty missing, not a number or <= 0";
        }
        res.samples.push_back(std::move(s));
    }
    return res;
}

// ===========================================================================
// Basins, k_sn and calibration
// ===========================================================================
namespace
{

constexpr double k_cal_nan = std::numeric_limits<double>::quiet_NaN();

std::string cal_num(double v, int precision = 10)
{
    if (std::isnan(v)) return "nan";
    std::ostringstream s;
    s << std::setprecision(precision) << v;
    return s.str();
}

std::string cal_csv(const std::string& s)
{
    if (s.find_first_of(",\"\n") == std::string::npos) return s;
    std::string q = "\"";
    for (char c : s) { if (c == '"') q += "\"\""; else q += c; }
    return q + "\"";
}

// 2.5 and 97.5 % of the finite values (compute_quantile_interpolated_in_place).
struct cal_interval
{
    double low = k_cal_nan, high = k_cal_nan;
    int n = 0;
};

cal_interval cal_percentile_interval(const std::vector<double>& values)
{
    std::vector<double> v;
    v.reserve(values.size());
    for (double x : values)
        if (std::isfinite(x)) v.push_back(x);
    cal_interval out;
    if (v.empty()) return out;
    std::vector<double> a = v;
    out.low = compute_quantile_interpolated_in_place(a, 0.025);
    a = v;
    out.high = compute_quantile_interpolated_in_place(a, 0.975);
    out.n = static_cast<int>(v.size());
    return out;
}

// Widen a positive interval around `point` by a relative systematic term s:
// each half width in ln becomes sqrt(h^2 + (1.96 s)^2). s = 0 changes nothing.
void widen_for_systematic(double point, double& low, double& high, double s)
{
    if (!(s > 0.0) || !(point > 0.0) || !(low > 0.0) || !(high > 0.0)) return;
    const double add = 1.96 * s;
    const double hl = std::log(point) - std::log(low), hh = std::log(high) - std::log(point);
    low = point * std::exp(-std::sqrt(hl * hl + add * add));
    high = point * std::exp(std::sqrt(hh * hh + add * add));
}

// Union-find over sample indices.
struct cal_union_find
{
    std::vector<int> parent;
    explicit cal_union_find(std::size_t n) : parent(n) { std::iota(parent.begin(), parent.end(), 0); }
    int find(int a) { while (parent[a] != a) a = parent[a] = parent[parent[a]]; return a; }
    void unite(int a, int b) { a = find(a); b = find(b); if (a != b) parent[std::max(a, b)] = std::min(a, b); }
};

bool status_is_kept(const std::string& s) { return s.empty() || s == "missing_uncertainty"; }

// One sample of a fit or estimator: x = ln k1, y = ln E, their sigmas, the
// nesting group.
struct cal_point
{
    double x = 0.0, y = 0.0, sx = 0.0, sy = 0.0;
    int group = -1;
};

// Weighted log mean at fixed n (estimator B) over `idx` of `pts`.
struct cal_fixed_n
{
    double lnK = k_cal_nan, se = k_cal_nan, dispersion = k_cal_nan, corr = k_cal_nan;
};

cal_fixed_n fixed_n_estimate(const std::vector<cal_point>& pts, const std::vector<std::size_t>& idx,
                             double n)
{
    cal_fixed_n out;
    if (idx.empty()) return out;
    double sw = 0.0, swk = 0.0;
    for (std::size_t i : idx)
    {
        const double w = 1.0 / (pts[i].sy * pts[i].sy + n * n * pts[i].sx * pts[i].sx);
        sw += w;
        swk += w * (pts[i].y - n * pts[i].x);
    }
    if (!(sw > 0.0) || !std::isfinite(sw)) return out;
    out.lnK = swk / sw;
    out.se = 1.0 / std::sqrt(sw);
    double ss = 0.0;
    for (std::size_t i : idx)
    {
        const double w = 1.0 / (pts[i].sy * pts[i].sy + n * n * pts[i].sx * pts[i].sx);
        const double d = pts[i].y - n * pts[i].x - out.lnK;
        ss += w * d * d;
    }
    out.dispersion = std::sqrt(ss / sw);
    if (idx.size() >= 3)
    {
        double mx = 0.0, mk = 0.0;
        for (std::size_t i : idx) { mx += pts[i].x; mk += pts[i].y - n * pts[i].x; }
        mx /= static_cast<double>(idx.size());
        mk /= static_cast<double>(idx.size());
        double sxx = 0.0, skk = 0.0, sxk = 0.0;
        for (std::size_t i : idx)
        {
            const double dx = pts[i].x - mx, dk = pts[i].y - n * pts[i].x - mk;
            sxx += dx * dx; skk += dk * dk; sxk += dx * dk;
        }
        if (sxx > 0.0 && skk > 0.0) out.corr = sxk / std::sqrt(sxx * skk);
    }
    return out;
}

// The distinct groups of idx and, per group, its members.
std::map<int, std::vector<std::size_t>> members_by_group(const std::vector<cal_point>& pts,
                                                        const std::vector<std::size_t>& idx)
{
    std::map<int, std::vector<std::size_t>> m;
    for (std::size_t i : idx) m[pts[i].group].push_back(i);
    return m;
}

// R replicates of a group bootstrap: each replicate lists the drawn members.
std::vector<std::vector<std::size_t>> draw_group_replicates(
    const std::map<int, std::vector<std::size_t>>& groups, int n_rep, std::mt19937_64& rng)
{
    std::vector<std::vector<std::size_t>> reps;
    if (groups.empty() || n_rep <= 0) return reps;
    std::vector<const std::vector<std::size_t>*> g;
    for (const auto& kv : groups) g.push_back(&kv.second);
    std::uniform_int_distribution<std::size_t> pick(0, g.size() - 1);
    reps.resize(static_cast<std::size_t>(n_rep));
    for (auto& r : reps)
        for (std::size_t k = 0; k < g.size(); ++k)
        {
            const auto* members = g[pick(rng)];
            r.insert(r.end(), members->begin(), members->end());
        }
    return reps;
}

std::size_t distinct_groups(const std::vector<cal_point>& pts, const std::vector<std::size_t>& idx)
{
    std::set<int> g;
    for (std::size_t i : idx) g.insert(pts[i].group);
    return g.size();
}

// Replicate medians of a reach bootstrap.
std::vector<double> reach_median_replicates(const std::vector<double>& k, int n_rep, std::mt19937_64& rng)
{
    std::vector<double> out;
    if (k.empty() || n_rep <= 0) return out;
    out.resize(static_cast<std::size_t>(n_rep));
    std::uniform_int_distribution<std::size_t> pick(0, k.size() - 1);
    std::vector<double> pooled(k.size());
    for (auto& v : out)
    {
        for (auto& p : pooled) p = k[pick(rng)];
        v = compute_median_in_place(pooled);
    }
    return out;
}

double sigma_ln_from_interval(double low, double high)
{
    if (!(low > 0.0) || !(high > 0.0)) return k_cal_nan;
    return (std::log(high) - std::log(low)) / 3.92;
}

} // namespace

std::mt19937_64 erodibility_calibration_bootstrap_engine(std::int64_t seed, double m_over_n)
{
    const long long mn = std::llround(m_over_n * 10000.0);
    std::seed_seq seq{static_cast<std::uint32_t>(seed), static_cast<std::uint32_t>(mn), 2u};
    return std::mt19937_64(seq);
}

std::size_t count_erosion_rate_samples_with_status(const std::vector<lsdtt_erosion_rate_sample>& samples,
                                                   const std::string& status)
{
    return static_cast<std::size_t>(std::count_if(samples.begin(), samples.end(),
        [&](const lsdtt_erosion_rate_sample& s) { return s.status == status; }));
}

// ---------------------------------------------------------------------------
// Basins
// ---------------------------------------------------------------------------
lsdtt_erosion_rate_basins_result delineate_erosion_rate_basins(
    const lsdtt_flow_base& flow, std::vector<lsdtt_erosion_rate_sample>& samples,
    const lsdtt_erosion_rate_basin_options& o, const lsdtt_lithology_units_on_grid* units)
{
    const std::string& method = o.snapping_method;
    if (method != "reported_area" && method != "expanding_square" && method != "downstream" &&
        method != "largest_area" && method != "none")
        throw std::invalid_argument("erosion_rate_snapping_method must be one of reported_area, "
                                    "expanding_square, downstream, largest_area, none (got '" +
                                    method + "')");
    if (o.nested_policy != "keep" && o.nested_policy != "outermost" && o.nested_policy != "innermost")
        throw std::invalid_argument("erosion_rate_nested_basins must be one of keep, outermost, "
                                    "innermost (got '" + o.nested_policy + "')");
    if (!(o.area_ratio_tolerance > 1.0) || !std::isfinite(o.area_ratio_tolerance))
        throw std::invalid_argument("erosion_rate_area_mismatch_tolerance must be > 1");
    if (!(o.snap_max_distance_m >= 0.0) || !std::isfinite(o.snap_max_distance_m))
        throw std::invalid_argument("erosion_rate_snap_max_distance_metres must be >= 0");
    if (o.channel_threshold_pixels == 0)
        throw std::invalid_argument("the channel threshold must be > 0 pixels");

    lsdtt_erosion_rate_basins_result res;
    res.options = o;
    res.has_lithology = units != nullptr;
    const lsdtt_raster_meta& meta = flow.meta();
    const double dx = flow.cell_size();
    res.search_radius_pixels = static_cast<int>(std::ceil(o.snap_max_distance_m / dx - 1e-9));
    res.basins.resize(samples.size());

    lsdtt_square_snap_options sq;
    sq.contributing_pixels_threshold = o.channel_threshold_pixels;
    sq.min_search_radius_pixels = 0;
    sq.max_search_radius_pixels = res.search_radius_pixels;
    lsdtt_downstream_snap_options ds;
    ds.contributing_pixels_threshold = o.channel_threshold_pixels;
    lsdtt_reported_area_snap_options ao;
    ao.contributing_pixels_threshold = o.channel_threshold_pixels;
    ao.max_search_radius_pixels = res.search_radius_pixels;
    ao.area_ratio_tolerance = o.area_ratio_tolerance;

    const double tol_ln = std::log(o.area_ratio_tolerance);
    for (std::size_t i = 0; i < samples.size(); ++i)
    {
        auto& smp = samples[i];
        auto& b = res.basins[i];
        if (smp.point.status != lsdtt_point_input_status::success || smp.point.row < 0)
        {
            smp.status = "outside_dem";
            smp.status_detail = "the point is not on the DEM grid";
            ++res.n_outside;
            continue;
        }
        const bool has_area = std::isfinite(smp.reported_area_m2) && smp.reported_area_m2 > 0.0;
        b.snapping_method = method == "reported_area" && !has_area ? "expanding_square" : method;
        const lsdtt_snap_result snap = snap_pixel_by_method(
            flow, smp.point.row, smp.point.col, b.snapping_method, sq, ds,
            has_area ? smp.reported_area_m2 : 1.0, ao);
        if (snap.status != lsdtt_snap_status::success || snap.snapped_flat < 0)
        {
            smp.status = "snap_failed";
            smp.status_detail = std::string("snapping: ") + snap_status_to_string(snap.status);
            ++res.n_snap_failed;
            continue;
        }
        b.outlet_flat = snap.snapped_flat;
        b.outlet_row = snap.snapped_row;
        b.outlet_col = snap.snapped_col;
        b.snapped_x = meta.x_min + (b.outlet_col + 0.5) * dx;
        b.snapped_y = meta.y_min + (meta.n_rows - b.outlet_row - 0.5) * dx;
        b.snap_distance_m = std::hypot(b.snapped_x - smp.point.x, b.snapped_y - smp.point.y);
        const lsdtt_basin_cpu basin(b.outlet_flat, flow);
        b.dem_area_m2 = static_cast<double>(basin.get_number_of_basin_pixels()) * dx * dx;
        b.complete = basin.check_basin_completeness();
        if (has_area)
        {
            b.area_ratio = b.dem_area_m2 / smp.reported_area_m2;
            b.area_check = std::fabs(std::log(b.area_ratio)) <= tol_ln + 1e-12 ? "ok" : "mismatch";
        }
        else
        {
            b.area_check = "no_reported_area";
            ++res.n_no_reported_area;
        }
        if (!b.complete)
        {
            smp.status = "incomplete_basin";
            smp.status_detail = "the basin reaches the DEM edge or nodata";
            ++res.n_incomplete;
            continue;
        }
        if (smp.status == "bad_rate")
        {
            ++res.n_bad_rate;
            continue;
        }
        if (b.area_check == "mismatch")
        {
            if (o.area_mismatch_excludes)
            {
                smp.status = "area_mismatch";
                smp.status_detail = "DEM area / reported area = " + cal_num(b.area_ratio, 6) +
                                    ", outside 1/" + cal_num(o.area_ratio_tolerance, 6) + ".." +
                                    cal_num(o.area_ratio_tolerance, 6);
                ++res.n_area_mismatch;
                continue;
            }
            ++res.n_area_flagged;
        }
    }

    // Nesting among the samples still kept.
    std::vector<int> kept;
    for (std::size_t i = 0; i < samples.size(); ++i)
        if (status_is_kept(samples[i].status) && res.basins[i].complete) kept.push_back(static_cast<int>(i));
    std::unordered_map<std::int64_t, std::vector<int>> samples_at;
    std::vector<std::int64_t> distinct;
    for (int i : kept)
    {
        auto& v = samples_at[res.basins[i].outlet_flat];
        if (v.empty()) distinct.push_back(res.basins[i].outlet_flat);
        v.push_back(i);
    }
    for (const auto& kv : samples_at)
        for (int i : kv.second)
            for (int j : kv.second)
                if (i != j) res.basins[i].duplicates.push_back(j);
    for (std::int64_t oa : distinct)
    {
        for (const std::int64_t f : flow.upslope_nodes(oa))
        {
            if (f == oa) continue;
            const auto it = samples_at.find(f);
            if (it == samples_at.end()) continue;
            for (int b_inner : it->second)
            {
                for (int a_outer : samples_at[oa])
                {
                    res.basins[b_inner].nested_in.push_back(a_outer);
                    ++res.basins[a_outer].n_samples_inside;
                }
            }
        }
    }
    for (int i : kept)
    {
        auto& b = res.basins[i];
        std::sort(b.nested_in.begin(), b.nested_in.end());
        std::sort(b.duplicates.begin(), b.duplicates.end());
        if (!b.nested_in.empty()) ++res.n_nested;
        if (!b.duplicates.empty()) ++res.n_duplicates;
    }

    // The nested-basin policy.
    if (o.nested_policy != "keep")
        for (int i : kept)
        {
            const auto& b = res.basins[i];
            const bool drop = o.nested_policy == "outermost" ? !b.nested_in.empty() : b.n_samples_inside > 0;
            if (!drop) continue;
            samples[i].status = "nested_excluded";
            samples[i].status_detail = o.nested_policy == "outermost"
                                           ? "nested in another sample (erosion_rate_nested_basins: outermost)"
                                           : "holds another sample (erosion_rate_nested_basins: innermost)";
            ++res.n_nested_excluded;
        }

    // Nesting groups over the samples still kept.
    cal_union_find uf(samples.size());
    std::vector<int> final_kept;
    for (int i : kept)
        if (status_is_kept(samples[i].status)) final_kept.push_back(i);
    std::vector<std::uint8_t> is_final(samples.size(), 0);
    for (int i : final_kept) is_final[static_cast<std::size_t>(i)] = 1;
    for (int i : final_kept)
    {
        for (int j : res.basins[i].duplicates)
            if (is_final[static_cast<std::size_t>(j)]) uf.unite(i, j);
        for (int j : res.basins[i].nested_in)
            if (is_final[static_cast<std::size_t>(j)]) uf.unite(i, j);
    }
    std::map<int, int> group_of_root;
    for (int i : final_kept)
    {
        const int root = uf.find(i);
        auto it = group_of_root.find(root);
        if (it == group_of_root.end())
            it = group_of_root.emplace(root, static_cast<int>(group_of_root.size())).first;
        res.basins[i].nesting_group = it->second;
        if (std::find(res.outlets.begin(), res.outlets.end(), res.basins[i].outlet_flat) == res.outlets.end())
            res.outlets.push_back(res.basins[i].outlet_flat);
    }
    res.n_groups = group_of_root.size();

    // Lithologic composition of every complete basin.
    if (units != nullptr)
    {
        if (units->code_of_pixel.size() != static_cast<std::size_t>(meta.n_rows) * meta.n_cols)
            throw std::invalid_argument("the lithologic units are not on the flow grid");
        for (std::size_t i = 0; i < samples.size(); ++i)
        {
            auto& b = res.basins[i];
            if (!b.complete || b.outlet_flat < 0) continue;
            std::map<std::int32_t, lsdtt_erosion_rate_basin_unit> by;
            std::vector<std::int64_t> channel;
            double total_area = 0.0;
            for (const std::int64_t f : flow.upslope_nodes(b.outlet_flat))
            {
                const std::int32_t c = units->code_of_pixel[static_cast<std::size_t>(f)];
                auto& u = by[c];
                u.code = c;
                u.area_m2 += dx * dx;
                total_area += dx * dx;
                if (flow.s_vector_size(f) >= o.channel_threshold_pixels) channel.push_back(f);
            }
            const std::vector<double> len = lithology_node_step_lengths(flow, channel);
            for (std::size_t k = 0; k < channel.size(); ++k)
            {
                auto& u = by[units->code_of_pixel[static_cast<std::size_t>(channel[k])]];
                u.channel_length_m += len[k];
                b.channel_length_m += len[k];
            }
            double best = -1.0;
            for (const auto& kv : by)
            {
                b.units.push_back(kv.second);
                if (kv.first >= 0 && kv.second.area_m2 > best)
                {
                    best = kv.second.area_m2;
                    b.dominant_code = kv.first;
                    b.dominant_area_fraction = kv.second.area_m2 / total_area;
                    b.dominant_channel_length_fraction =
                        b.channel_length_m > 0.0 ? kv.second.channel_length_m / b.channel_length_m : k_cal_nan;
                }
            }
            const auto un = by.find(-1);
            b.unmapped_area_fraction = un == by.end() ? 0.0 : un->second.area_m2 / total_area;
        }
    }
    return res;
}

// ---------------------------------------------------------------------------
// Basin k_sn
// ---------------------------------------------------------------------------
std::vector<lsdtt_erosion_rate_basin_ksn> compute_erosion_rate_basin_ksn(
    const lsdtt_chi_network_fit_collection_result& fits, const lsdtt_raster_meta& meta,
    const std::vector<float>& drainage_area_m2, const lsdtt_erodibility_calibration_options& opts,
    std::mt19937_64& rng)
{
    const std::string& stat = opts.ksn_statistic;
    if (stat != "reach_median" && stat != "chi_z_slope" && stat != "segment_median")
        throw std::invalid_argument("erosion_rate_ksn_statistic must be one of reach_median, "
                                    "chi_z_slope, segment_median (got '" + stat + "')");
    if (opts.min_reach_pixels < 2)
        throw std::invalid_argument("erosion_rate_min_reach_pixels must be >= 2");

    std::vector<lsdtt_erosion_rate_basin_ksn> out(fits.basins.size());
    if (fits.basins.empty()) return out;

    // Reaches split at confluences only: one unit everywhere.
    const std::size_t npix = static_cast<std::size_t>(meta.n_rows) * meta.n_cols;
    lsdtt_lithology_units_on_grid single;
    single.code_of_pixel.assign(npix, 0);
    lsdtt_lithology_unit_label label;
    label.code = 0;
    label.name = "basin";
    single.units.push_back(label);
    const std::vector<float> no_contact(npix, std::numeric_limits<float>::quiet_NaN());
    lsdtt_chi_lithology_options lo;
    lo.contact_buffer_m = 0.0;
    lo.min_reach_pixels = opts.min_reach_pixels;
    lo.bootstrap_replicates = 0;
    const lsdtt_chi_lithology_tables t =
        build_chi_lithology_tables(fits, meta, single, no_contact, drainage_area_m2, lo);
    const lsdtt_chi_ksn_basin_summary_table summary = summarise_chi_ksn_by_basin(fits);

    std::vector<std::vector<double>> reach_k(fits.basins.size());
    std::vector<std::size_t> n_reaches(fits.basins.size(), 0);
    for (const auto& r : t.reaches)
    {
        if (r.basin_index < 0 || static_cast<std::size_t>(r.basin_index) >= out.size()) continue;
        ++n_reaches[static_cast<std::size_t>(r.basin_index)];
        if (r.qualifies && std::isfinite(r.k_sn_ols))
            reach_k[static_cast<std::size_t>(r.basin_index)].push_back(r.k_sn_ols);
    }

    for (std::size_t b = 0; b < fits.basins.size(); ++b)
    {
        const auto& f = fits.basins[b];
        auto& k = out[b];
        k.outlet_flat = f.outlet_flat;
        k.n_reaches = n_reaches[b];
        k.n_reaches_qualifying = reach_k[b].size();
        if (!reach_k[b].empty())
        {
            std::vector<double> v = reach_k[b];
            k.reach_median = compute_median_in_place(v);
            v = reach_k[b];
            const double q25 = compute_quantile_interpolated_in_place(v, 0.25);
            v = reach_k[b];
            const double q75 = compute_quantile_interpolated_in_place(v, 0.75);
            if (k.reach_median != 0.0) k.iqr_over_median = (q75 - q25) / k.reach_median;
        }
        // chi-z slope through the origin.
        double z0 = k_cal_nan;
        for (std::size_t i = 0; i < f.flats.size(); ++i)
            if (f.flats[i] == f.outlet_flat) z0 = f.observed_elevation[i];
        if (std::isfinite(z0) && !f.flats.empty())
        {
            double sxz = 0.0, sxx = 0.0, sz = 0.0;
            for (std::size_t i = 0; i < f.flats.size(); ++i)
            {
                const double c = f.chi[i], dz = f.observed_elevation[i] - z0;
                sxz += c * dz;
                sxx += c * c;
                sz += dz;
            }
            if (sxx > 0.0)
            {
                k.chi_z_slope = sxz / sxx;
                const double mz = sz / static_cast<double>(f.flats.size());
                double ssr = 0.0, sst = 0.0;
                for (std::size_t i = 0; i < f.flats.size(); ++i)
                {
                    const double dz = f.observed_elevation[i] - z0;
                    const double r = dz - k.chi_z_slope * f.chi[i];
                    ssr += r * r;
                    sst += (dz - mz) * (dz - mz);
                }
                if (sst > 0.0) k.chi_z_r_squared = 1.0 - ssr / sst;
            }
        }
        if (b < summary.rows.size() && summary.rows[b].n_segments > 0)
        {
            k.segment_median = summary.rows[b].median_k_sn;
            k.chi_weighted_mean = summary.rows[b].chi_weighted_mean_k_sn;
        }

        k.k_sn = stat == "reach_median" ? k.reach_median
                 : stat == "chi_z_slope" ? k.chi_z_slope
                                         : k.segment_median;
        // Reach replicates (drawn for every basin with reaches, whatever the
        // statistic, so the stream does not depend on it), scaled to the
        // chosen statistic.
        std::vector<double> reps = reach_median_replicates(reach_k[b], opts.bootstrap_replicates, rng);
        if (!reps.empty() && std::isfinite(k.k_sn) && k.reach_median > 0.0)
        {
            const double scale = stat == "reach_median" ? 1.0 : k.k_sn / k.reach_median;
            if (scale != 1.0)
                for (auto& r : reps) r *= scale;
            const cal_interval ci = cal_percentile_interval(reps);
            k.ci_low = ci.low;
            k.ci_high = ci.high;
            k.sigma_ln_k = sigma_ln_from_interval(ci.low, ci.high);
            k.replicate_k = std::move(reps);
        }
    }
    return out;
}

// ---------------------------------------------------------------------------
// Calibration
// ---------------------------------------------------------------------------
lsdtt_erodibility_calibration_result calibrate_erodibility_from_erosion_rates(
    std::vector<lsdtt_erosion_rate_sample>& samples, const lsdtt_erosion_rate_basins_result& basins,
    const std::vector<lsdtt_erosion_rate_basin_ksn>& ksn,
    const lsdtt_chi_network_fit_collection_result& fits,
    const lsdtt_chi_lithology_tables* unit_tables, const lsdtt_lithology_units_on_grid* units,
    double m_over_n, double A_0, const lsdtt_erodibility_calibration_options& opts,
    std::mt19937_64& rng)
{
    if (opts.n_values.empty()) throw std::invalid_argument("stream_power_n_values is empty");
    for (double n : opts.n_values)
        if (!std::isfinite(n) || !(n > 0.0))
            throw std::invalid_argument("stream_power_n_values must all be finite and > 0");
    if (opts.min_reaches_per_basin < 1 || opts.min_unit_reaches_per_basin < 1 ||
        opts.min_groups_for_n_fit < 3 || opts.bootstrap_replicates < 0)
        throw std::invalid_argument("erodibility calibration counts must be positive "
                                    "(min_groups_for_n_fit >= 3, bootstrap_replicates >= 0)");
    if (!(opts.mono_min_area_fraction > 0.0 && opts.mono_min_area_fraction <= 1.0))
        throw std::invalid_argument("erosion_rate_mono_lithology_min_area_fraction must be in (0, 1]");
    if (!(opts.systematic_relative_uncertainty >= 0.0) || !std::isfinite(opts.systematic_relative_uncertainty))
        throw std::invalid_argument("erosion_rate_systematic_relative_uncertainty must be >= 0");
    if (samples.size() != basins.basins.size())
        throw std::invalid_argument("samples and basins differ in size");

    lsdtt_erodibility_calibration_result res;
    res.m_over_n = m_over_n;
    res.A_0 = A_0;
    res.ksn_statistic = opts.ksn_statistic;
    res.mono_min_area_fraction = opts.mono_min_area_fraction;
    const int R = opts.bootstrap_replicates;
    const double a0_factor = std::pow(A_0, m_over_n);   // k1 = k * A0^(m/n)
    const double s_sys = opts.systematic_relative_uncertainty;

    std::unordered_map<std::int64_t, std::size_t> ksn_of_outlet;
    for (std::size_t i = 0; i < ksn.size(); ++i) ksn_of_outlet[ksn[i].outlet_flat] = i;
    auto ksn_for = [&](std::size_t si) -> const lsdtt_erosion_rate_basin_ksn* {
        const auto it = ksn_of_outlet.find(basins.basins[si].outlet_flat);
        return it == ksn_of_outlet.end() ? nullptr : &ksn[it->second];
    };

    // too_few_reaches / used.
    for (std::size_t i = 0; i < samples.size(); ++i)
    {
        auto& smp = samples[i];
        if (!status_is_kept(smp.status) || !basins.basins[i].complete) continue;
        const auto* k = ksn_for(i);
        if (k == nullptr || k->n_reaches_qualifying < static_cast<std::size_t>(opts.min_reaches_per_basin) ||
            !(k->k_sn > 0.0))
        {
            smp.status = "too_few_reaches";
            smp.status_detail = k == nullptr
                                    ? "no channel network at the outlet"
                                    : std::to_string(k->n_reaches_qualifying) + " qualifying reach(es), fewer than " +
                                          std::to_string(opts.min_reaches_per_basin);
            ++res.n_too_few_reaches;
            continue;
        }
        if (smp.status.empty()) smp.status = "used";
    }

    // E draws per used sample.
    std::vector<std::vector<double>> e_draws(samples.size());
    for (std::size_t i = 0; i < samples.size(); ++i)
    {
        const auto& smp = samples[i];
        if (smp.status != "used" || R <= 0) continue;
        std::normal_distribution<double> nd(smp.E_m_per_yr, smp.sigma_E_m_per_yr);
        auto& d = e_draws[i];
        d.resize(static_cast<std::size_t>(R));
        for (auto& e : d)
        {
            double v = nd(rng);
            int guard = 0;
            while (!(v > 0.0) && guard++ < 1000) v = nd(rng);
            e = v > 0.0 ? v : k_cal_nan;
        }
    }

    // (A) per sample and n.
    for (double n : opts.n_values)
    {
        const double m = m_over_n * n;
        for (std::size_t i = 0; i < samples.size(); ++i)
        {
            const auto& smp = samples[i];
            if (smp.status != "used" && !(smp.status == "missing_uncertainty" && ksn_for(i) &&
                                          basins.basins[i].complete && ksn_for(i)->k_sn > 0.0 &&
                                          ksn_for(i)->n_reaches_qualifying >=
                                              static_cast<std::size_t>(opts.min_reaches_per_basin)))
                continue;
            const auto* k = ksn_for(i);
            lsdtt_erodibility_basin_K_row row;
            row.sample = static_cast<int>(i);
            row.n = n;
            row.m = m;
            row.k_sn = k->k_sn;
            row.k_sn_ci_low = k->ci_low;
            row.k_sn_ci_high = k->ci_high;
            row.K = stream_power_erodibility_from_erosion_rate(smp.E_m_per_yr, k->k_sn, m, n, A_0);
            if (!e_draws[i].empty() && k->replicate_k.size() == e_draws[i].size())
            {
                std::vector<double> Ks(e_draws[i].size());
                for (std::size_t r = 0; r < Ks.size(); ++r)
                    Ks[r] = stream_power_erodibility_from_erosion_rate(e_draws[i][r], k->replicate_k[r], m, n, A_0);
                const cal_interval ci = cal_percentile_interval(Ks);
                row.K_ci_low = ci.low;
                row.K_ci_high = ci.high;
                row.n_bootstrap = ci.n;
                widen_for_systematic(row.K, row.K_ci_low, row.K_ci_high, s_sys);
            }
            res.basin_K.push_back(row);
        }
    }

    // The points of the fits.
    std::vector<cal_point> pts;
    std::vector<std::size_t> all;
    std::vector<double> k1s;
    for (std::size_t i = 0; i < samples.size(); ++i)
    {
        if (samples[i].status != "used") continue;
        const auto* k = ksn_for(i);
        cal_point p;
        p.x = std::log(k->k_sn * a0_factor);
        p.y = std::log(samples[i].E_m_per_yr);
        p.sx = std::isfinite(k->sigma_ln_k) ? k->sigma_ln_k : 0.0;
        p.sy = samples[i].sigma_E_m_per_yr / samples[i].E_m_per_yr;
        p.group = basins.basins[i].nesting_group;
        all.push_back(pts.size());
        pts.push_back(p);
        k1s.push_back(k->k_sn * a0_factor);
        res.fit_samples.push_back(static_cast<int>(i));
    }
    res.n_used = pts.size();
    res.n_groups = distinct_groups(pts, all);
    double k_range = k_cal_nan, median_k1 = k_cal_nan;
    if (!k1s.empty())
    {
        k_range = *std::max_element(k1s.begin(), k1s.end()) / *std::min_element(k1s.begin(), k1s.end());
        std::vector<double> v = k1s;
        median_k1 = compute_median_in_place(v);
    }
    const auto groups = members_by_group(pts, all);
    const auto reps = draw_group_replicates(groups, R, rng);

    // (B) fixed n.
    for (double n : opts.n_values)
    {
        lsdtt_erodibility_fit_row f;
        f.fit_type = "fixed_n";
        f.method = "weighted_log_mean";
        f.N_basins = pts.size();
        f.N_groups = res.n_groups;
        f.ksn_range_ratio = k_range;
        f.n = n;
        f.m = m_over_n * n;
        if (pts.empty())
        {
            f.status = "no_basins";
            res.fits.push_back(f);
            continue;
        }
        const cal_fixed_n e = fixed_n_estimate(pts, all, n);
        f.status = "ok";
        f.K = std::exp(e.lnK);
        f.ln_K_se = e.se;
        f.ln_K_dispersion = e.dispersion;
        f.corr_lnK_lnk = e.corr;
        f.E_at_median_ksn_m_per_yr = f.K * std::pow(median_k1, n);
        std::vector<double> Kr;
        for (const auto& r : reps)
        {
            const cal_fixed_n er = fixed_n_estimate(pts, r, n);
            Kr.push_back(std::isfinite(er.lnK) ? std::exp(er.lnK) : k_cal_nan);
        }
        const cal_interval ci = cal_percentile_interval(Kr);
        f.K_ci_low = ci.low;
        f.K_ci_high = ci.high;
        f.n_bootstrap = ci.n;
        widen_for_systematic(f.K, f.K_ci_low, f.K_ci_high, s_sys);
        res.fits.push_back(f);
    }

    // (C) fitted n.
    {
        lsdtt_erodibility_fit_row f;
        f.fit_type = "fitted_n";
        f.method = "york";
        f.N_basins = pts.size();
        f.N_groups = res.n_groups;
        f.ksn_range_ratio = k_range;
        // First-order SE(n) = s_res / (s_x sqrt(N)) from an OLS of y on x.
        if (pts.size() >= 3)
        {
            double mx = 0.0, my = 0.0;
            for (const auto& p : pts) { mx += p.x; my += p.y; }
            mx /= static_cast<double>(pts.size());
            my /= static_cast<double>(pts.size());
            double sxx = 0.0, sxy = 0.0;
            for (const auto& p : pts) { sxx += (p.x - mx) * (p.x - mx); sxy += (p.x - mx) * (p.y - my); }
            if (sxx > 0.0)
            {
                const double b = sxy / sxx;
                double ssr = 0.0;
                for (const auto& p : pts) { const double r = p.y - my - b * (p.x - mx); ssr += r * r; }
                const double nn = static_cast<double>(pts.size());
                const double s_res = std::sqrt(ssr / (nn - 2.0));
                const double s_x = std::sqrt(sxx / nn);
                f.n_needed_estimate_se = s_res / (s_x * std::sqrt(nn));
            }
        }
        if (pts.empty())
            f.status = "no_basins";
        else if (!opts.fit_n)
            f.status = "not_fitted_disabled";
        else if (res.n_groups < static_cast<std::size_t>(opts.min_groups_for_n_fit))
            f.status = "not_fitted_too_few_groups";
        else if (!(k_range >= opts.min_ksn_range_ratio))
            f.status = "not_fitted_narrow_ksn_range";
        else
        {
            std::vector<double> x, y, sx, sy;
            for (const auto& p : pts) { x.push_back(p.x); y.push_back(p.y); sx.push_back(p.sx); sy.push_back(p.sy); }
            const lsdtt_york_fit yf = fit_york_line(x, y, sx, sy);
            if (!yf.converged)
                f.status = "not_converged";
            else
            {
                f.status = "ok";
                f.n = yf.slope;
                f.m = m_over_n * yf.slope;
                f.K = std::exp(yf.intercept);
                f.n_se_york = yf.se_slope;
                f.n_se_york_mswd_scaled = yf.se_slope_scaled;
                f.ln_K_se = yf.se_intercept;
                f.ln_K_se_mswd_scaled = yf.se_intercept_scaled;
                f.corr_n_lnK = yf.corr_slope_intercept;
                f.mswd = yf.mswd;
                f.E_at_median_ksn_m_per_yr = f.K * std::pow(median_k1, f.n);
                const cal_fixed_n e = fixed_n_estimate(pts, all, f.n);
                f.ln_K_dispersion = e.dispersion;
                f.corr_lnK_lnk = e.corr;
                std::vector<double> nr, Kr;
                for (const auto& r : reps)
                {
                    if (distinct_groups(pts, r) < 3) continue;
                    std::vector<double> rx, ry, rsx, rsy;
                    for (std::size_t i : r) { rx.push_back(pts[i].x); ry.push_back(pts[i].y); rsx.push_back(pts[i].sx); rsy.push_back(pts[i].sy); }
                    const lsdtt_york_fit rf = fit_york_line(rx, ry, rsx, rsy);
                    if (!rf.converged) continue;
                    nr.push_back(rf.slope);
                    Kr.push_back(std::exp(rf.intercept));
                }
                const cal_interval cn = cal_percentile_interval(nr);
                const cal_interval ck = cal_percentile_interval(Kr);
                f.n_ci_low = cn.low;
                f.n_ci_high = cn.high;
                f.K_ci_low = ck.low;
                f.K_ci_high = ck.high;
                f.n_bootstrap = cn.n;
                widen_for_systematic(f.K, f.K_ci_low, f.K_ci_high, s_sys);
            }
        }
        res.fits.push_back(f);
    }

    // (D) per lithologic unit.
    if (units != nullptr)
    {
        const bool want_mono = std::find(opts.unit_k_methods.begin(), opts.unit_k_methods.end(), "mono") !=
                               opts.unit_k_methods.end();
        const bool want_uib = unit_tables != nullptr &&
                              std::find(opts.unit_k_methods.begin(), opts.unit_k_methods.end(),
                                        "unit_in_basin") != opts.unit_k_methods.end();
        std::vector<lsdtt_erodibility_unit_row> mono_rows, uib_rows;

        auto estimate_unit = [&](const std::vector<cal_point>& up, const std::string& method,
                                 const lsdtt_lithology_unit_label& u, std::vector<lsdtt_erodibility_unit_row>& out) {
            std::vector<std::size_t> idx(up.size());
            std::iota(idx.begin(), idx.end(), std::size_t{0});
            const auto g = members_by_group(up, idx);
            const auto ureps = draw_group_replicates(g, up.empty() ? 0 : R, rng);
            for (double n : opts.n_values)
            {
                lsdtt_erodibility_unit_row row;
                row.lith_code = u.code;
                row.lith_name = u.name;
                row.method = method;
                row.n = n;
                row.m = m_over_n * n;
                row.N_basins = up.size();
                row.N_groups = g.size();
                row.low_sample = g.size() < 3;
                if (up.empty())
                {
                    row.status = "no_basins";
                    out.push_back(row);
                    continue;
                }
                row.status = "ok";
                const cal_fixed_n e = fixed_n_estimate(up, idx, n);
                row.K = std::exp(e.lnK);
                std::vector<double> Kr;
                for (const auto& r : ureps)
                {
                    const cal_fixed_n er = fixed_n_estimate(up, r, n);
                    Kr.push_back(std::isfinite(er.lnK) ? std::exp(er.lnK) : k_cal_nan);
                }
                const cal_interval ci = cal_percentile_interval(Kr);
                row.K_ci_low = ci.low;
                row.K_ci_high = ci.high;
                row.n_bootstrap = ci.n;
                widen_for_systematic(row.K, row.K_ci_low, row.K_ci_high, s_sys);
                out.push_back(row);
            }
        };

        if (want_mono)
            for (const auto& u : units->units)
            {
                std::vector<cal_point> up;
                for (std::size_t j = 0; j < res.fit_samples.size(); ++j)
                {
                    const auto& b = basins.basins[static_cast<std::size_t>(res.fit_samples[j])];
                    if (b.dominant_code == u.code && b.dominant_area_fraction >= opts.mono_min_area_fraction)
                        up.push_back(pts[j]);
                }
                estimate_unit(up, "mono", u, mono_rows);
            }

        if (want_uib)
        {
            // Qualifying lithologic reaches per (unit, fit basin).
            std::map<std::pair<std::int32_t, std::int32_t>, std::vector<double>> reach_k;
            for (const auto& r : unit_tables->reaches)
                if (r.qualifies && r.lith_code >= 0 && std::isfinite(r.k_sn_ols))
                    reach_k[{r.lith_code, r.basin_index}].push_back(r.k_sn_ols);
            // Used samples by outlet.
            std::unordered_map<std::int64_t, std::vector<std::size_t>> used_at;
            for (std::size_t j = 0; j < res.fit_samples.size(); ++j)
                used_at[basins.basins[static_cast<std::size_t>(res.fit_samples[j])].outlet_flat].push_back(j);
            for (const auto& u : units->units)
            {
                std::vector<cal_point> up;
                for (std::size_t bi = 0; bi < fits.basins.size(); ++bi)
                {
                    const auto it = reach_k.find({u.code, static_cast<std::int32_t>(bi)});
                    if (it == reach_k.end() ||
                        it->second.size() < static_cast<std::size_t>(opts.min_unit_reaches_per_basin))
                        continue;
                    const auto ua = used_at.find(fits.basins[bi].outlet_flat);
                    if (ua == used_at.end()) continue;
                    std::vector<double> v = it->second;
                    const double kmed = compute_median_in_place(v);
                    const std::vector<double> kr = reach_median_replicates(it->second, R, rng);
                    const cal_interval ci = cal_percentile_interval(kr);
                    const double sxl = sigma_ln_from_interval(ci.low, ci.high);
                    if (!(kmed > 0.0)) continue;
                    for (std::size_t j : ua->second)
                    {
                        cal_point p = pts[j];
                        p.x = std::log(kmed * a0_factor);
                        p.sx = std::isfinite(sxl) ? sxl : 0.0;
                        up.push_back(p);
                    }
                }
                estimate_unit(up, "unit_in_basin", u, uib_rows);
            }
        }
        // Order: method, n, code.
        for (const auto* rows : {&mono_rows, &uib_rows})
            for (double n : opts.n_values)
                for (const auto& r : *rows)
                    if (r.n == n) res.units.push_back(r);
    }
    return res;
}

// ---------------------------------------------------------------------------
// Writers
// ---------------------------------------------------------------------------
namespace
{

std::string join_ids(const std::vector<int>& idx, const std::vector<lsdtt_erosion_rate_sample>& samples)
{
    std::string out;
    for (int i : idx) out += (out.empty() ? "" : ";") + samples[static_cast<std::size_t>(i)].sample_id;
    return out;
}

std::string unit_name_of(const lsdtt_lithology_units_on_grid* units, std::int32_t code)
{
    if (code < 0) return "unmapped";
    if (units)
        for (const auto& u : units->units)
            if (u.code == code) return u.name;
    return "";
}

std::ofstream open_out(const std::string& path)
{
    std::ofstream out(path);
    if (!out) throw std::runtime_error("cannot write " + path);
    return out;
}

} // namespace

void write_erodibility_calibration_samples_csv(
    const std::string& path, const std::vector<lsdtt_erosion_rate_sample>& samples,
    const lsdtt_erosion_rate_basins_result& basins, const std::vector<lsdtt_erosion_rate_basin_ksn>& ksn,
    const lsdtt_erodibility_calibration_result& cal, const std::string& input_units,
    const lsdtt_lithology_units_on_grid* units)
{
    std::ofstream out = open_out(path);
    out << "sample_id,input_row,input_latitude,input_longitude,snapped_latitude,snapped_longitude,"
           "snapping_method,snap_distance_m,outlet_flat,status,status_detail,reported_area_km2,"
           "dem_area_km2,area_ratio,area_check,nesting_group,nested_in,n_samples_inside,"
           "duplicate_outlet_of,erosion_rate_input,erosion_rate_uncertainty_input,"
           "erosion_rate_input_units,erosion_rate_m_per_yr,erosion_rate_uncertainty_m_per_yr,"
           "m_over_n,A_0,k_sn_statistic,k_sn,k_sn_ci_low,k_sn_ci_high,n_reaches,"
           "n_reaches_qualifying,k_sn_reach_iqr_over_median,k_sn_chi_z_slope,chi_z_r_squared,"
           "k_sn_segment_median,k_sn_chi_weighted_mean,dominant_lith_code,dominant_lith_name,"
           "dominant_area_fraction,dominant_channel_length_fraction,unmapped_area_fraction\n";
    std::unordered_map<std::int64_t, std::size_t> ksn_of;
    for (std::size_t i = 0; i < ksn.size(); ++i) ksn_of[ksn[i].outlet_flat] = i;
    for (std::size_t i = 0; i < samples.size(); ++i)
    {
        const auto& s = samples[i];
        const auto& b = basins.basins[i];
        const auto it = b.outlet_flat >= 0 ? ksn_of.find(b.outlet_flat) : ksn_of.end();
        const lsdtt_erosion_rate_basin_ksn* k = it == ksn_of.end() ? nullptr : &ksn[it->second];
        const lsdtt_erosion_rate_basin_ksn none;
        const auto& kk = k ? *k : none;
        out << cal_csv(s.sample_id) << "," << s.input_row << "," << cal_num(b.input_latitude, 12) << ","
            << cal_num(b.input_longitude, 12) << "," << cal_num(b.snapped_latitude, 12) << ","
            << cal_num(b.snapped_longitude, 12) << "," << b.snapping_method << ","
            << cal_num(b.snap_distance_m) << "," << b.outlet_flat << "," << s.status << ","
            << cal_csv(s.status_detail) << "," << cal_num(s.reported_area_m2 / 1e6) << ","
            << cal_num(b.dem_area_m2 / 1e6) << "," << cal_num(b.area_ratio) << "," << b.area_check << ","
            << b.nesting_group << "," << cal_csv(join_ids(b.nested_in, samples)) << ","
            << b.n_samples_inside << "," << cal_csv(join_ids(b.duplicates, samples)) << ","
            << cal_num(s.rate_input) << "," << cal_num(s.uncertainty_input) << "," << input_units << ","
            << cal_num(s.E_m_per_yr) << "," << cal_num(s.sigma_E_m_per_yr) << ","
            << cal_num(cal.m_over_n) << "," << cal_num(cal.A_0) << "," << cal.ksn_statistic << ","
            << cal_num(kk.k_sn) << "," << cal_num(kk.ci_low) << "," << cal_num(kk.ci_high) << ","
            << kk.n_reaches << "," << kk.n_reaches_qualifying << "," << cal_num(kk.iqr_over_median) << ","
            << cal_num(kk.chi_z_slope) << "," << cal_num(kk.chi_z_r_squared) << ","
            << cal_num(kk.segment_median) << "," << cal_num(kk.chi_weighted_mean) << ","
            << b.dominant_code << "," << cal_csv(basins.has_lithology && b.complete
                                                     ? unit_name_of(units, b.dominant_code) : "")
            << "," << cal_num(b.dominant_area_fraction) << ","
            << cal_num(b.dominant_channel_length_fraction) << "," << cal_num(b.unmapped_area_fraction)
            << "\n";
    }
}

void write_erodibility_calibration_basin_K_csv(
    const std::string& path, const std::vector<lsdtt_erosion_rate_sample>& samples,
    const lsdtt_erosion_rate_basins_result& basins, const lsdtt_erodibility_calibration_result& cal)
{
    std::ofstream out = open_out(path);
    out << "sample_id,status,nesting_group,m_over_n,A_0,n,m,erosion_rate_m_per_yr,"
           "erosion_rate_uncertainty_m_per_yr,k_sn,k_sn_ci_low,k_sn_ci_high,K,K_ci_low,K_ci_high,"
           "n_bootstrap,K_units\n";
    for (const auto& r : cal.basin_K)
    {
        const auto& s = samples[static_cast<std::size_t>(r.sample)];
        out << cal_csv(s.sample_id) << "," << s.status << ","
            << basins.basins[static_cast<std::size_t>(r.sample)].nesting_group << ","
            << cal_num(cal.m_over_n) << "," << cal_num(cal.A_0) << "," << cal_num(r.n) << ","
            << cal_num(r.m) << "," << cal_num(s.E_m_per_yr) << "," << cal_num(s.sigma_E_m_per_yr) << ","
            << cal_num(r.k_sn) << "," << cal_num(r.k_sn_ci_low) << "," << cal_num(r.k_sn_ci_high) << ","
            << cal_num(r.K) << "," << cal_num(r.K_ci_low) << "," << cal_num(r.K_ci_high) << ","
            << r.n_bootstrap << "," << lithology_erodibility_units_string(r.m) << "\n";
    }
}

void write_erodibility_calibration_fit_csv(const std::string& path,
                                           const lsdtt_erodibility_calibration_result& cal)
{
    std::ofstream out = open_out(path);
    out << "fit_type,method,status,N_basins,N_groups,ksn_range_ratio,m_over_n,A_0,n,n_ci_low,"
           "n_ci_high,n_se_york,n_se_york_mswd_scaled,m,K,K_ci_low,K_ci_high,ln_K_se,"
           "ln_K_se_mswd_scaled,corr_n_lnK,E_at_median_ksn_m_per_yr,mswd,ln_K_dispersion,"
           "corr_lnK_lnk,n_bootstrap,bootstrap_cluster,K_units\n";
    for (const auto& f : cal.fits)
        out << f.fit_type << "," << f.method << "," << f.status << "," << f.N_basins << ","
            << f.N_groups << "," << cal_num(f.ksn_range_ratio) << "," << cal_num(cal.m_over_n) << ","
            << cal_num(cal.A_0) << "," << cal_num(f.n) << "," << cal_num(f.n_ci_low) << ","
            << cal_num(f.n_ci_high) << "," << cal_num(f.n_se_york) << ","
            << cal_num(f.n_se_york_mswd_scaled) << "," << cal_num(f.m) << "," << cal_num(f.K) << ","
            << cal_num(f.K_ci_low) << "," << cal_num(f.K_ci_high) << "," << cal_num(f.ln_K_se) << ","
            << cal_num(f.ln_K_se_mswd_scaled) << "," << cal_num(f.corr_n_lnK) << ","
            << cal_num(f.E_at_median_ksn_m_per_yr) << "," << cal_num(f.mswd) << ","
            << cal_num(f.ln_K_dispersion) << "," << cal_num(f.corr_lnK_lnk) << "," << f.n_bootstrap
            << ",nesting_group," << (std::isfinite(f.m) ? lithology_erodibility_units_string(f.m) : "")
            << "\n";
}

void write_erodibility_calibration_by_lithology_csv(const std::string& path,
                                                    const lsdtt_erodibility_calibration_result& cal)
{
    std::ofstream out = open_out(path);
    out << "lith_code,lith_name,method,m_over_n,A_0,n,m,N_basins,N_groups,K,K_ci_low,K_ci_high,"
           "n_bootstrap,low_sample,status,purity_threshold,K_units\n";
    for (const auto& r : cal.units)
        out << r.lith_code << "," << cal_csv(r.lith_name) << "," << r.method << ","
            << cal_num(cal.m_over_n) << "," << cal_num(cal.A_0) << "," << cal_num(r.n) << ","
            << cal_num(r.m) << "," << r.N_basins << "," << r.N_groups << "," << cal_num(r.K) << ","
            << cal_num(r.K_ci_low) << "," << cal_num(r.K_ci_high) << "," << r.n_bootstrap << ","
            << (r.low_sample ? 1 : 0) << "," << r.status << ","
            << (r.method == "mono" ? cal_num(cal.mono_min_area_fraction, 6) : std::string("nan")) << ","
            << lithology_erodibility_units_string(r.m) << "\n";
}

void write_erodibility_calibration_basin_lithology_csv(
    const std::string& path, const std::vector<lsdtt_erosion_rate_sample>& samples,
    const lsdtt_erosion_rate_basins_result& basins, const lsdtt_lithology_units_on_grid& units,
    double /*cell_size*/)
{
    std::ofstream out = open_out(path);
    out << "sample_id,lith_code,lith_name,area_km2,area_fraction,channel_length_m,"
           "channel_length_fraction\n";
    for (std::size_t i = 0; i < samples.size(); ++i)
    {
        const auto& b = basins.basins[i];
        if (!b.complete || b.units.empty()) continue;
        double total = 0.0;
        for (const auto& u : b.units) total += u.area_m2;
        for (const auto& u : b.units)
            out << cal_csv(samples[i].sample_id) << "," << u.code << ","
                << cal_csv(unit_name_of(&units, u.code)) << "," << cal_num(u.area_m2 / 1e6) << ","
                << cal_num(u.area_m2 / total) << "," << cal_num(u.channel_length_m) << ","
                << cal_num(b.channel_length_m > 0.0 ? u.channel_length_m / b.channel_length_m : k_cal_nan)
                << "\n";
    }
}

std::size_t write_erodibility_calibration_basin_polygons(
    const std::string& path, const lsdtt_flow_base& flow,
    const std::vector<lsdtt_erosion_rate_sample>& samples,
    const lsdtt_erosion_rate_basins_result& basins, const std::string& crs_wkt,
    vector_output_format format)
{
    lsdtt_vector_schema schema;
    schema.field_names = {"sample_id", "status", "nesting_group", "dem_area_km2", "reported_area_km2",
                          "area_ratio"};
    schema.field_types = {lsdtt_vector_field_type::string_type, lsdtt_vector_field_type::string_type,
                          lsdtt_vector_field_type::int64_type, lsdtt_vector_field_type::double_type,
                          lsdtt_vector_field_type::double_type, lsdtt_vector_field_type::double_type};
    lsdtt_data_loader loader;
    loader.begin_vector_write(path, "erodibility_calibration_basins", schema, crs_wkt, format, false);
    std::size_t n = 0;
    auto dv = [](double v) -> lsdtt_vector_field_value {
        if (std::isfinite(v)) return v;
        return std::monostate{};
    };
    for (std::size_t i = 0; i < samples.size(); ++i)
    {
        const auto& b = basins.basins[i];
        if (!b.complete || b.outlet_flat < 0) continue;
        const lsdtt_basin_cpu basin(b.outlet_flat, flow);
        lsdtt_vector_record rec = basin.get_basin_outline().get_outline_as_polygon_record(static_cast<std::int32_t>(i));
        rec.field_values = {samples[i].sample_id, samples[i].status.empty() ? std::string("kept") : samples[i].status,
                            static_cast<std::int64_t>(b.nesting_group), dv(b.dem_area_m2 / 1e6),
                            dv(samples[i].reported_area_m2 / 1e6), dv(b.area_ratio)};
        if (write_polygon_record(loader, rec)) ++n;
    }
    loader.finish_vector_write();
    return n;
}

std::size_t write_erodibility_calibration_property_table(
    const std::string& path, const lsdtt_erodibility_calibration_result& cal, double n,
    const std::string& method, const lsdtt_lithology_units_on_grid& units,
    const std::string& provenance)
{
    std::vector<lsdtt_lithology_erodibility_table_row> rows;
    for (const auto& u : units.units)
    {
        lsdtt_lithology_erodibility_table_row r;
        r.code = u.code;
        r.name = u.name;
        for (const auto& w : cal.units)
            if (w.lith_code == u.code && w.method == method && w.n == n)
            {
                r.K = w.K;
                std::ostringstream ref;
                ref << method << " from " << w.N_basins << " basin entries in " << w.N_groups
                    << " nesting groups; K 95 % interval [" << cal_num(w.K_ci_low, 8) << ", "
                    << cal_num(w.K_ci_high, 8) << "]" << (w.low_sample ? "; low_sample" : "");
                r.reference = ref.str();
            }
        rows.push_back(r);
    }
    const double m = cal.m_over_n * n;
    return write_lithology_erodibility_property_table(
        path, m, n, provenance,
        {"K units: " + lithology_erodibility_units_string(m) + " (map area in m^2)",
         "bootstrap: nesting_group clusters; the interval in each reference column is the 95 % "
         "percentile interval of K"},
        rows, "no basin with enough reaches in this unit; add rows by hand before mapping a region "
              "that contains them");
}

} // namespace lsdtt
