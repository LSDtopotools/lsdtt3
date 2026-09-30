// lsdtt_manning_hydraulics_cpu.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html

#include "lsdtt_manning_hydraulics_cpu.hpp"
#include "../core/lsdtt_data_loader.hpp"
#include "../core/lsdtt_vector_utils.hpp"
#include "../core/lsdtt_flow_base.hpp"
#include "../core/lsdtt_types.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <iterator>
#include <limits>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <unordered_map>

namespace lsdtt
{
namespace
{

// Quote-aware CSV line splitter (handles "a,b" and "" escapes).
std::vector<std::string> split_csv_line(const std::string& line)
{
    std::vector<std::string> out;
    std::string field;
    bool in_quotes = false;

    auto finish = [&out, &field]()
    {
        while (!field.empty() && (field.back() == '\r' || field.back() == ' '))
            field.pop_back();
        std::size_t lead = 0;
        while (lead < field.size() && field[lead] == ' ') ++lead;
        out.push_back(field.substr(lead));
        field.clear();
    };

    for (std::size_t i = 0; i < line.size(); ++i)
    {
        const char c = line[i];
        if (in_quotes)
        {
            if (c == '"')
            {
                if (i + 1 < line.size() && line[i + 1] == '"')
                {
                    field.push_back('"');
                    ++i;
                }
                else
                {
                    in_quotes = false;
                }
            }
            else
            {
                field.push_back(c);
            }
        }
        else if (c == '"')
        {
            in_quotes = true;
        }
        else if (c == ',')
        {
            finish();
        }
        else
        {
            field.push_back(c);
        }
    }
    finish();
    return out;
}

int find_col(const std::vector<std::string>& header, const std::string& name)
{
    for (std::size_t i = 0; i < header.size(); ++i)
        if (header[i] == name) return static_cast<int>(i);
    throw std::runtime_error(
        "lsdtt_manning_hydraulics_cpu: missing required column: " + name);
}

double parse_double_field(
    const std::string& text, const std::string& where, std::size_t line_no)
{
    try
    {
        std::size_t pos = 0;
        const double v = std::stod(text, &pos);
        if (pos != text.size()) throw std::invalid_argument("trailing characters");
        return v;
    }
    catch (const std::exception&)
    {
        throw std::runtime_error(
            where + ": line " + std::to_string(line_no) +
            ": cannot parse '" + text + "' as a number.");
    }
}

long long parse_int_field(
    const std::string& text, const std::string& where, std::size_t line_no)
{
    try
    {
        std::size_t pos = 0;
        const long long v = std::stoll(text, &pos);
        if (pos != text.size()) throw std::invalid_argument("trailing characters");
        return v;
    }
    catch (const std::exception&)
    {
        throw std::runtime_error(
            where + ": line " + std::to_string(line_no) +
            ": cannot parse '" + text + "' as an integer.");
    }
}

// Extract a sorted, valid (s, z, x, y) profile from a cross section's samples.
bool extract_and_sort_profile(
    const lsdtt_cross_section& section,
    std::vector<double>& s,
    std::vector<double>& z,
    std::vector<double>& x,
    std::vector<double>& y)
{
    s.clear(); z.clear(); x.clear(); y.clear();

    for (const auto& sample : section.samples)
    {
        if (sample.valid == 0) continue;
        if (!std::isfinite(static_cast<double>(sample.elevation))) continue;

        s.push_back(sample.offset_metres);
        z.push_back(static_cast<double>(sample.elevation));
        x.push_back(sample.x);
        y.push_back(sample.y);
    }
    if (s.size() < 3) return false;

    std::vector<std::size_t> order(s.size());
    std::iota(order.begin(), order.end(), std::size_t{0});
    std::sort(order.begin(), order.end(),
        [&](std::size_t a, std::size_t b) { return s[a] < s[b]; });

    std::vector<double> s2(s.size()), z2(s.size()), x2(s.size()), y2(s.size());
    for (std::size_t i = 0; i < order.size(); ++i)
    {
        s2[i] = s[order[i]]; z2[i] = z[order[i]];
        x2[i] = x[order[i]]; y2[i] = y[order[i]];
    }
    s.swap(s2); z.swap(z2); x.swap(x2); y.swap(y2);
    return true;
}

// Index of the lowest sample within +/- half_width_m of the section centre
// (offset 0). If no sample lies inside the window, falls back to the sample
// nearest the centre.
std::size_t find_thalweg_index(
    const std::vector<double>& s,
    const std::vector<double>& z,
    double half_width_m)
{
    std::size_t best = s.size();
    for (std::size_t i = 0; i < s.size(); ++i)
    {
        if (std::abs(s[i]) > half_width_m) continue;
        if (best == s.size() || z[i] < z[best] ||
            (z[i] == z[best] && std::abs(s[i]) < std::abs(s[best])))
            best = i;
    }

    if (best == s.size())
    {
        best = 0;
        for (std::size_t i = 1; i < s.size(); ++i)
            if (std::abs(s[i]) < std::abs(s[best])) best = i;
    }
    return best;
}

// Finds the first crossing of the water-surface elevation moving outward
// from the thalweg, interpolating both the along-section offset (s) and the
// raster-CRS coordinates (x, y). direction: +1 = right, -1 = left.
lsdtt_bank_crossing find_bank_crossing(
    const std::vector<double>& s, const std::vector<double>& z,
    const std::vector<double>& x, const std::vector<double>& y,
    std::size_t thalweg_idx, double wse, int direction)
{
    lsdtt_bank_crossing c;
    const std::size_t n = s.size();

    if (direction > 0)
    {
        for (std::size_t i = thalweg_idx; i + 1 < n; ++i)
        {
            const bool wet0 = z[i] < wse, wet1 = z[i + 1] < wse;
            if (wet0 && !wet1)
            {
                const double t = (wse - z[i]) / (z[i + 1] - z[i]);
                c.found = true;
                c.s = s[i] + t * (s[i + 1] - s[i]);
                c.x = x[i] + t * (x[i + 1] - x[i]);
                c.y = y[i] + t * (y[i + 1] - y[i]);
                return c;
            }
        }
    }
    else
    {
        for (std::size_t i = thalweg_idx; i > 0; --i)
        {
            const std::size_t j = i - 1;
            const bool wet_i = z[i] < wse, wet_j = z[j] < wse;
            if (wet_i && !wet_j)
            {
                const double t = (wse - z[i]) / (z[j] - z[i]);
                c.found = true;
                c.s = s[i] + t * (s[j] - s[i]);
                c.x = x[i] + t * (x[j] - x[i]);
                c.y = y[i] + t * (y[j] - y[i]);
                return c;
            }
        }
    }
    return c;
}

// The D8 channel path from source to outlet, in downstream order.
// distance_m is strictly increasing, which the slope window relies on.
struct channel_path
{
    std::vector<std::int64_t> flats;
    std::vector<double>       distance_m;
    std::unordered_map<std::int64_t, std::size_t> index_by_flat;
};

// Walks downstream from source_flat via the D8 receiver chain, accumulating
// true along-channel distance (cardinal = cell_size, diagonal = cell_size *
// sqrt(2)). Direction codes (lsdtt_flowinfo_cpu):
//   7  0  1
//   6  x  2
//   5  4  3
// Odd codes are diagonal. Throws if the walk exceeds the number of valid
// cells (corrupt, cyclic graph).
channel_path build_channel_path(
    const lsdtt_flow_base& flow,
    std::int64_t           source_flat)
{
    channel_path path;

    const double cell_step = flow.cell_size();
    const double diag_step = cell_step * std::sqrt(2.0);

    auto push = [&path](std::int64_t flat, double distance)
    {
        path.index_by_flat[flat] = path.flats.size();
        path.flats.push_back(flat);
        path.distance_m.push_back(distance);
    };

    double distance = 0.0;
    std::int64_t current = source_flat;
    const std::int64_t max_steps = flow.n_valid_cells() + 1;
    std::int64_t steps = 0;

    push(current, 0.0);

    while (true)
    {
        const std::int64_t receiver = flow.receiver_of_unchecked(current);

        if (receiver == lsdtt::k_nodata_receiver ||
            receiver == lsdtt::k_outlet_receiver)
        {
            break;
        }

        const std::uint8_t dir = flow.receiver_dir_of_unchecked(current);
        const bool is_diagonal = (dir % 2 == 1);
        distance += is_diagonal ? diag_step : cell_step;

        push(receiver, distance);
        current = receiver;

        if (++steps > max_steps)
        {
            throw std::runtime_error(
                "build_channel_path: downstream walk exceeded the number "
                "of valid cells -- the flow graph contains a cycle.");
        }
    }

    return path;
}

// Least-squares slope of z against along-channel distance d, over the path
// cells within +/- half_length_m of path index `idx`. Returns false (slope
// untouched) if the fit is not usable; see the header for the criteria.
// Elevation decreases downstream while d increases, so the positive
// down-channel slope is -cov(d, z) / var(d).
bool fit_along_path_slope(
    const std::vector<double>&       d,
    const std::vector<double>&       z,
    const std::vector<std::uint8_t>& ok,
    std::size_t                      idx,
    double                           half_length_m,
    double&                          slope)
{
    const double d0 = d[idx];

    const std::size_t lo = static_cast<std::size_t>(
        std::lower_bound(d.begin(), d.end(), d0 - half_length_m) - d.begin());
    const std::size_t hi = static_cast<std::size_t>(
        std::upper_bound(d.begin(), d.end(), d0 + half_length_m) - d.begin());

    std::size_t n = 0;
    double sum_d = 0.0, sum_z = 0.0;
    double d_min = std::numeric_limits<double>::max();
    double d_max = std::numeric_limits<double>::lowest();

    for (std::size_t i = lo; i < hi; ++i)
    {
        if (!ok[i]) continue;
        const double di = d[i] - d0;
        sum_d += di;
        sum_z += z[i];
        d_min = std::min(d_min, di);
        d_max = std::max(d_max, di);
        ++n;
    }

    if (n < 3) return false;
    if ((d_max - d_min) < half_length_m) return false;

    const double mean_d = sum_d / static_cast<double>(n);
    const double mean_z = sum_z / static_cast<double>(n);

    double sxx = 0.0, sxz = 0.0;
    for (std::size_t i = lo; i < hi; ++i)
    {
        if (!ok[i]) continue;
        const double dd = (d[i] - d0) - mean_d;
        sxx += dd * dd;
        sxz += dd * (z[i] - mean_z);
    }

    if (!(sxx > 0.0)) return false;

    const double s = -sxz / sxx;
    if (!std::isfinite(s) || s <= 0.0) return false;

    slope = s;
    return true;
}

// Bisection on a discharge function that is non-decreasing in depth.
// Returns false if the target is not bracketed by [min_depth, max_depth].
template <class discharge_fn>
bool bisect_depth(
    discharge_fn&&                      q_at_depth,
    double                              target_discharge_m3s,
    const lsdtt_manning_solver_options& opts,
    double&                             solved_depth_m)
{
    double lo = opts.min_depth_m, hi = opts.max_depth_m;
    const double q_lo = q_at_depth(lo);
    const double q_hi = q_at_depth(hi);
    if (!(q_lo <= target_discharge_m3s && target_discharge_m3s <= q_hi))
        return false;

    for (int it = 0; it < opts.max_iterations; ++it)
    {
        const double mid = 0.5 * (lo + hi);
        const double q_mid = q_at_depth(mid);
        if (q_mid < target_discharge_m3s) lo = mid; else hi = mid;
        if ((hi - lo) < opts.depth_tolerance_m)
        {
            solved_depth_m = 0.5 * (lo + hi);
            return true;
        }
    }
    solved_depth_m = 0.5 * (lo + hi);
    return true;
}

} // unnamed namespace

// -----------------------------------------------------------------------------
// Reach and profile I/O
// -----------------------------------------------------------------------------
std::vector<lsdtt_manning_reach> read_manning_reaches_csv(
    const std::string& csv_path)
{
    const std::string where = "read_manning_reaches_csv(" + csv_path + ")";

    std::ifstream in(csv_path);
    if (!in.is_open())
        throw std::runtime_error("read_manning_reaches_csv: could not open " + csv_path);

    std::string header_line;
    std::getline(in, header_line);
    const std::vector<std::string> header = split_csv_line(header_line);

    const int c_start = find_col(header, "start_distance_m");
    const int c_end   = find_col(header, "end_distance_m");
    const int c_n     = find_col(header, "manning_n");
    const int c_q     = find_col(header, "target_discharge_m3s");

    std::vector<lsdtt_manning_reach> reaches;
    std::string line;
    std::size_t line_no = 1;
    while (std::getline(in, line))
    {
        ++line_no;
        if (line.empty() || line == "\r") continue;
        const std::vector<std::string> t = split_csv_line(line);
        if (t.size() < header.size())
        {
            throw std::runtime_error(
                where + ": line " + std::to_string(line_no) + " has " +
                std::to_string(t.size()) + " fields, expected " +
                std::to_string(header.size()) + ".");
        }

        lsdtt_manning_reach r;
        r.start_distance_m     = parse_double_field(t[c_start], where, line_no);
        r.end_distance_m       = parse_double_field(t[c_end], where, line_no);
        r.manning_n            = parse_double_field(t[c_n], where, line_no);
        r.target_discharge_m3s = parse_double_field(t[c_q], where, line_no);
        reaches.push_back(r);
    }
    return reaches;
}

void validate_manning_reaches(const std::vector<lsdtt_manning_reach>& reaches)
{
    if (reaches.empty())
        throw std::invalid_argument("validate_manning_reaches: no reaches supplied.");

    for (std::size_t i = 0; i < reaches.size(); ++i)
    {
        const auto& r = reaches[i];
        const std::string tag = "reach " + std::to_string(i + 1) + ": ";
        if (!std::isfinite(r.start_distance_m) || !std::isfinite(r.end_distance_m) ||
            !(r.start_distance_m < r.end_distance_m))
            throw std::invalid_argument(
                "validate_manning_reaches: " + tag + "start_distance_m must be "
                "finite and strictly less than end_distance_m.");
        if (!(r.manning_n > 0.0) || !std::isfinite(r.manning_n))
            throw std::invalid_argument(
                "validate_manning_reaches: " + tag + "manning_n must be > 0.");
        if (!(r.target_discharge_m3s > 0.0) || !std::isfinite(r.target_discharge_m3s))
            throw std::invalid_argument(
                "validate_manning_reaches: " + tag + "target_discharge_m3s must be > 0.");
    }

    std::vector<std::size_t> order(reaches.size());
    std::iota(order.begin(), order.end(), std::size_t{0});
    std::sort(order.begin(), order.end(),
        [&](std::size_t a, std::size_t b)
        { return reaches[a].start_distance_m < reaches[b].start_distance_m; });

    for (std::size_t k = 1; k < order.size(); ++k)
    {
        const auto& prev = reaches[order[k - 1]];
        const auto& next = reaches[order[k]];
        if (next.start_distance_m < prev.end_distance_m)
            throw std::invalid_argument(
                "validate_manning_reaches: reach " + std::to_string(order[k] + 1) +
                " overlaps reach " + std::to_string(order[k - 1] + 1) +
                ". Reaches are half-open [start, end) and must not overlap.");
    }
}

lsdtt_channel_profile_lookup read_channel_profile_csv(
    const std::string& csv_path)
{
    const std::string where = "read_channel_profile_csv(" + csv_path + ")";

    std::ifstream in(csv_path);
    if (!in.is_open())
        throw std::runtime_error("read_channel_profile_csv: could not open " + csv_path);

    std::string header_line;
    std::getline(in, header_line);
    const std::vector<std::string> header = split_csv_line(header_line);

    const int c_id    = find_col(header, "cross_id");
    const int c_dist  = find_col(header, "distance_from_source_m");
    const int c_area  = find_col(header, "drainage_area_m2");
    const int c_slope = find_col(header, "slope_m_per_m");
    const int c_valid = find_col(header, "slope_valid");

    lsdtt_channel_profile_lookup lookup;
    std::string line;
    std::size_t line_no = 1;
    while (std::getline(in, line))
    {
        ++line_no;
        if (line.empty() || line == "\r") continue;
        const std::vector<std::string> t = split_csv_line(line);
        if (t.size() < header.size())
        {
            throw std::runtime_error(
                where + ": line " + std::to_string(line_no) + " has " +
                std::to_string(t.size()) + " fields, expected " +
                std::to_string(header.size()) + ".");
        }

        const long long cross_id = parse_int_field(t[c_id], where, line_no);
        lookup.distance_by_cross_id[cross_id] =
            parse_double_field(t[c_dist], where, line_no);
        lookup.drainage_area_by_cross_id[cross_id] =
            parse_double_field(t[c_area], where, line_no);

        // A row with slope_valid == 0 is deliberately left OUT of
        // slope_by_cross_id so the solve flags it instead of using a bad value.
        if (parse_int_field(t[c_valid], where, line_no) != 0)
        {
            lookup.slope_by_cross_id[cross_id] =
                parse_double_field(t[c_slope], where, line_no);
        }
    }
    return lookup;
}

lsdtt_channel_profile_lookup make_channel_profile_lookup(
    const std::vector<lsdtt_channel_profile_row>& rows)
{
    lsdtt_channel_profile_lookup lookup;
    for (const auto& r : rows)
    {
        lookup.distance_by_cross_id[r.cross_section_id]      = r.distance_from_source_m;
        lookup.drainage_area_by_cross_id[r.cross_section_id] = r.drainage_area_m2;
        if (r.slope_valid)
            lookup.slope_by_cross_id[r.cross_section_id] = r.slope_m_per_m;
    }
    return lookup;
}

lsdtt_manning_profile_coverage compute_manning_profile_coverage(
    const std::vector<lsdtt_cross_section>& sections,
    const lsdtt_channel_profile_lookup&     profile)
{
    lsdtt_manning_profile_coverage c;
    c.n_sections = sections.size();
    for (const auto& sec : sections)
    {
        if (profile.distance_by_cross_id.find(sec.cross_section_id) ==
            profile.distance_by_cross_id.end())
            ++c.n_missing_distance;
        else if (profile.slope_by_cross_id.find(sec.cross_section_id) ==
                 profile.slope_by_cross_id.end())
            ++c.n_missing_slope;
    }
    return c;
}

// -----------------------------------------------------------------------------
// Geometry and Manning discharge
// -----------------------------------------------------------------------------
lsdtt_wetted_geometry compute_wetted_geometry(
    const std::vector<double>& s,
    const std::vector<double>& z,
    double water_surface_elevation)
{
    lsdtt_wetted_geometry g;
    double area = 0.0, perimeter = 0.0;

    for (std::size_t i = 0; i + 1 < s.size(); ++i)
    {
        const double ds = s[i + 1] - s[i];
        if (ds <= 0.0) continue;
        const double h0 = std::max(0.0, water_surface_elevation - z[i]);
        const double h1 = std::max(0.0, water_surface_elevation - z[i + 1]);
        area += 0.5 * (h0 + h1) * ds;
        if (h0 > 0.0 || h1 > 0.0)
            perimeter += std::sqrt(ds * ds + (z[i + 1] - z[i]) * (z[i + 1] - z[i]));
    }

    if (area <= 0.0 || perimeter <= 0.0) return g;
    g.area_m2 = area;
    g.wetted_perimeter_m = perimeter;
    g.hydraulic_radius_m = area / perimeter;
    g.valid = true;
    return g;
}

// Walks outward from the thalweg in each direction and stops at the first
// sample at or above the water surface -- the same rule find_bank_crossing
// uses -- so area, wetted perimeter and reported bank crossings all describe
// one body of water. The final segment on each side is clipped at the
// interpolated crossing. Monotonic in water-surface elevation, but it can
// jump when the surface overtops a crest into a lower area.
lsdtt_wetted_geometry compute_connected_wetted_geometry(
    const std::vector<double>& s,
    const std::vector<double>& z,
    std::size_t thalweg_idx,
    double water_surface_elevation)
{
    lsdtt_wetted_geometry g;
    const std::size_t n = s.size();
    const double wse = water_surface_elevation;

    if (thalweg_idx >= n || !(z[thalweg_idx] < wse)) return g;

    double area = 0.0, perimeter = 0.0;

    // Right of the thalweg (increasing s).
    for (std::size_t i = thalweg_idx; i + 1 < n; ++i)
    {
        const double ds = s[i + 1] - s[i];
        const double dz = z[i + 1] - z[i];
        const double h0 = wse - z[i];

        if (z[i + 1] < wse)
        {
            if (ds <= 0.0) continue;
            const double h1 = wse - z[i + 1];
            area      += 0.5 * (h0 + h1) * ds;
            perimeter += std::sqrt(ds * ds + dz * dz);
        }
        else
        {
            if (ds > 0.0)
            {
                const double t = h0 / dz;          // dz >= h0 > 0
                area      += 0.5 * h0 * (t * ds);
                perimeter += t * std::sqrt(ds * ds + dz * dz);
            }
            break;
        }
    }

    // Left of the thalweg (decreasing s).
    for (std::size_t i = thalweg_idx; i > 0; --i)
    {
        const std::size_t j = i - 1;
        const double ds = s[i] - s[j];
        const double dz = z[j] - z[i];
        const double h0 = wse - z[i];

        if (z[j] < wse)
        {
            if (ds <= 0.0) continue;
            const double h1 = wse - z[j];
            area      += 0.5 * (h0 + h1) * ds;
            perimeter += std::sqrt(ds * ds + dz * dz);
        }
        else
        {
            if (ds > 0.0)
            {
                const double t = h0 / dz;
                area      += 0.5 * h0 * (t * ds);
                perimeter += t * std::sqrt(ds * ds + dz * dz);
            }
            break;
        }
    }

    if (area <= 0.0 || perimeter <= 0.0) return g;
    g.area_m2            = area;
    g.wetted_perimeter_m = perimeter;
    g.hydraulic_radius_m = area / perimeter;
    g.valid = true;
    return g;
}

double manning_discharge(
    double area_m2, double hydraulic_radius_m,
    double slope, double manning_n)
{
    if (manning_n <= 0.0 || slope <= 0.0 || area_m2 <= 0.0) return 0.0;
    return (1.0 / manning_n) * area_m2 *
           std::pow(hydraulic_radius_m, 2.0 / 3.0) * std::sqrt(slope);
}

bool solve_connected_bankfull_depth(
    const std::vector<double>& s,
    const std::vector<double>& z,
    std::size_t thalweg_idx,
    double slope,
    double manning_n,
    double target_discharge_m3s,
    const lsdtt_manning_solver_options& opts,
    double& solved_depth_m)
{
    if (thalweg_idx >= s.size() || s.size() != z.size()) return false;
    const double thalweg_z = z[thalweg_idx];

    auto q_at_depth = [&](double depth) -> double
    {
        const lsdtt_wetted_geometry g =
            compute_connected_wetted_geometry(s, z, thalweg_idx, thalweg_z + depth);
        if (!g.valid) return 0.0;
        return manning_discharge(g.area_m2, g.hydraulic_radius_m, slope, manning_n);
    };

    return bisect_depth(q_at_depth, target_discharge_m3s, opts, solved_depth_m);
}

// -----------------------------------------------------------------------------
// compute_manning_widths_for_reach
//
// Cross sections are matched to the reach by DISTANCE FROM SOURCE (looked up
// per cross section from `profile`), not by cross_id, which is assigned by
// raster row-major scan order and is unrelated to position along a sinuous
// channel. Membership is half-open: start <= d < end.
// -----------------------------------------------------------------------------
std::vector<lsdtt_manning_width_record> compute_manning_widths_for_reach(
    const std::vector<lsdtt_cross_section>& sections,
    const lsdtt_manning_reach&              reach,
    const lsdtt_channel_profile_lookup&     profile,
    const lsdtt_manning_solver_options&     opts)
{
    std::vector<lsdtt_manning_width_record> out;

    for (const auto& sec : sections)
    {
        const auto dist_it = profile.distance_by_cross_id.find(sec.cross_section_id);
        if (dist_it == profile.distance_by_cross_id.end())
            continue;   // no distance -> cannot be assigned to any reach
        const double distance_m = dist_it->second;

        if (!(distance_m >= reach.start_distance_m &&
              distance_m <  reach.end_distance_m))
            continue;

        lsdtt_manning_width_record rec;
        rec.cross_section_id = sec.cross_section_id;
        rec.centre_x = sec.centre_x;
        rec.centre_y = sec.centre_y;
        rec.distance_from_source_m = distance_m;

        const auto slope_it = profile.slope_by_cross_id.find(sec.cross_section_id);
        if (slope_it == profile.slope_by_cross_id.end())
        {
            rec.quality_flag = "no_slope_in_channel_profile_csv";
            out.push_back(rec);
            continue;
        }
        rec.slope_m_per_m = slope_it->second;

        const auto area_it = profile.drainage_area_by_cross_id.find(sec.cross_section_id);
        if (area_it != profile.drainage_area_by_cross_id.end())
            rec.drainage_area_m2 = area_it->second;

        std::vector<double> s, z, x, y;
        if (!extract_and_sort_profile(sec, s, z, x, y))
        {
            rec.quality_flag = "invalid_profile";
            out.push_back(rec);
            continue;
        }

        const std::size_t th_idx =
            find_thalweg_index(s, z, opts.thalweg_search_half_width_m);
        rec.thalweg_s = s[th_idx];
        rec.thalweg_z = z[th_idx];

        double depth = 0.0;
        if (!solve_connected_bankfull_depth(
                s, z, th_idx, rec.slope_m_per_m, reach.manning_n,
                reach.target_discharge_m3s, opts, depth))
        {
            rec.quality_flag = "target_discharge_not_bracketed";
            out.push_back(rec);
            continue;
        }

        rec.solved_depth_m = depth;
        rec.water_surface_elevation_m = rec.thalweg_z + depth;

        const lsdtt_wetted_geometry g =
            compute_connected_wetted_geometry(s, z, th_idx, rec.water_surface_elevation_m);
        if (!g.valid)
        {
            rec.quality_flag = "dry_at_solved_depth";
            out.push_back(rec);
            continue;
        }

        const lsdtt_bank_crossing left  =
            find_bank_crossing(s, z, x, y, th_idx, rec.water_surface_elevation_m, -1);
        const lsdtt_bank_crossing right =
            find_bank_crossing(s, z, x, y, th_idx, rec.water_surface_elevation_m, +1);

        if (!left.found || !right.found)
        {
            rec.quality_flag = "bank_crossing_not_found";
            out.push_back(rec);
            continue;
        }

        rec.left_bank_s  = left.s;  rec.left_bank_x  = left.x;  rec.left_bank_y  = left.y;
        rec.right_bank_s = right.s; rec.right_bank_x = right.x; rec.right_bank_y = right.y;
        rec.width_m = rec.right_bank_s - rec.left_bank_s;

        rec.area_m2            = g.area_m2;
        rec.wetted_perimeter_m = g.wetted_perimeter_m;
        rec.hydraulic_radius_m = g.hydraulic_radius_m;
        rec.discharge_at_solution_m3s =
            manning_discharge(g.area_m2, g.hydraulic_radius_m, rec.slope_m_per_m, reach.manning_n);

        // Connected area jumps when the water surface overtops a bank crest
        // into a lower area. If the target Q sits inside that jump, bisection
        // converges on the crest and the solved depth does not reproduce Q.
        // Flag it rather than report a width for the wrong discharge.
        if (reach.target_discharge_m3s > 0.0)
        {
            const double rel_err =
                std::abs(rec.discharge_at_solution_m3s - reach.target_discharge_m3s)
                / reach.target_discharge_m3s;
            if (rel_err > opts.max_discharge_relative_error)
            {
                rec.quality_flag = "spill_discontinuity";
                out.push_back(rec);
                continue;
            }
        }

        rec.valid = true;
        rec.quality_flag = "ok";
        out.push_back(rec);
    }

    return out;
}

std::vector<lsdtt_manning_width_record> compute_manning_widths_all_reaches(
    const std::vector<lsdtt_cross_section>& sections,
    const std::vector<lsdtt_manning_reach>& reaches,
    const lsdtt_channel_profile_lookup&     profile,
    const lsdtt_manning_solver_options&     opts)
{
    validate_manning_reaches(reaches);

    std::vector<lsdtt_manning_width_record> all;
    for (const auto& reach : reaches)
    {
        auto reach_results = compute_manning_widths_for_reach(sections, reach, profile, opts);
        all.insert(all.end(),
            std::make_move_iterator(reach_results.begin()),
            std::make_move_iterator(reach_results.end()));
    }
    std::sort(all.begin(), all.end(),
        [](const lsdtt_manning_width_record& a, const lsdtt_manning_width_record& b)
        { return a.distance_from_source_m < b.distance_from_source_m; });
    return all;
}

// -----------------------------------------------------------------------------
// Writers
// -----------------------------------------------------------------------------
void write_manning_bank_points(
    const std::vector<lsdtt_manning_width_record>& results,
    const lsdtt_raster_meta&                        meta,
    const std::string&                              out_path,
    point_output_format                             format,
    bool                                             spatial_index,
    const std::string&                               layer_name)
{
    if (format == point_output_format::csv)
    {
        throw std::invalid_argument(
            "write_manning_bank_points: csv is not supported for bank points "
            "(a point CSV must carry EPSG:4326 latitude/longitude). Use "
            "flatgeobuf or gpkg, and write_manning_widths_csv for tabular output.");
    }

    lsdtt_data_loader loader;

    lsdtt_vector_schema schema;
    schema.field_names = {"cross_id", "side", "width_m", "solved_depth_m", "quality_flag"};
    schema.field_types = {
        lsdtt_vector_field_type::int64_type,
        lsdtt_vector_field_type::string_type,
        lsdtt_vector_field_type::double_type,
        lsdtt_vector_field_type::double_type,
        lsdtt_vector_field_type::string_type
    };

    const vector_output_format out_format = convert_point_to_vector_output_format(format);

    loader.begin_vector_write(out_path, layer_name, schema, meta.crs_wkt, out_format, spatial_index);

    for (const auto& r : results)
    {
        if (!r.valid) continue;

        loader.write_point(r.left_bank_x, r.left_bank_y,
            { static_cast<std::int64_t>(r.cross_section_id), std::string("left"),
              r.width_m, r.solved_depth_m, r.quality_flag });

        loader.write_point(r.right_bank_x, r.right_bank_y,
            { static_cast<std::int64_t>(r.cross_section_id), std::string("right"),
              r.width_m, r.solved_depth_m, r.quality_flag });
    }

    loader.finish_vector_write();
}

std::size_t write_manning_widths_csv(
    const std::vector<lsdtt_manning_width_record>& results,
    const std::string&                              out_path)
{
    std::ofstream out(out_path);
    if (!out.is_open())
        throw std::runtime_error("write_manning_widths_csv: could not open " + out_path);

    out << "cross_id,distance_from_source_m,centre_x,centre_y,slope_m_per_m,drainage_area_m2,"
        << "thalweg_s,thalweg_z,solved_depth_m,"
        << "water_surface_elevation_m,left_bank_s,right_bank_s,"
        << "left_bank_x,left_bank_y,right_bank_x,right_bank_y,width_m,"
        << "area_m2,wetted_perimeter_m,hydraulic_radius_m,"
        << "discharge_at_solution_m3s,valid,quality_flag\n";
    out << std::fixed << std::setprecision(4);

    std::size_t n_valid = 0;
    for (const auto& r : results)
    {
        // Slope spans orders of magnitude; fixed 4-decimal output would
        // truncate low-gradient channels, so format it separately.
        char slope_buf[32];
        std::snprintf(slope_buf, sizeof(slope_buf), "%.8g", r.slope_m_per_m);

        out << r.cross_section_id << "," << r.distance_from_source_m << ","
            << r.centre_x << "," << r.centre_y << ","
            << slope_buf << "," << r.drainage_area_m2 << ","
            << r.thalweg_s << "," << r.thalweg_z << "," << r.solved_depth_m << ","
            << r.water_surface_elevation_m << "," << r.left_bank_s << "," << r.right_bank_s << ","
            << r.left_bank_x << "," << r.left_bank_y << ","
            << r.right_bank_x << "," << r.right_bank_y << ","
            << r.width_m << "," << r.area_m2 << "," << r.wetted_perimeter_m << ","
            << r.hydraulic_radius_m << "," << r.discharge_at_solution_m3s << ","
            << (r.valid ? 1 : 0) << "," << r.quality_flag << "\n";
        if (r.valid) ++n_valid;
    }
    return n_valid;
}

// -----------------------------------------------------------------------------
// Channel profile
// -----------------------------------------------------------------------------
std::vector<lsdtt_channel_profile_row> compute_channel_profile(
    const std::vector<lsdtt_cross_section>& sections,
    const lsdtt_flow_base&                  flow,
    const lsdtt_raster_cpu&                 filled_dem,
    double                                   slope_window_half_length_metres,
    std::int64_t                             source_flat,
    double                                   thalweg_search_half_width_metres)
{
    if (!(slope_window_half_length_metres > 0.0) ||
        slope_window_half_length_metres < 2.0 * flow.cell_size())
    {
        throw std::runtime_error(
            "compute_channel_profile: slope_window_half_length_metres must be "
            "at least two cell sizes (" + std::to_string(2.0 * flow.cell_size()) +
            " m); got " + std::to_string(slope_window_half_length_metres) + ".");
    }

    const channel_path path = build_channel_path(flow, source_flat);

    // Drainage area from the D8 flow object already built.
    // NOTE: verify the unit this raster carries (pixel counts vs m2) -- the
    // output column is labelled drainage_area_m2.
    const lsdtt_raster_result drainage_area_result = flow.get_drainage_area_raster();
    const auto& da_data = drainage_area_result.data;
    const int   n_cols  = drainage_area_result.meta.n_cols;

    const int n_rows_dem = filled_dem.n_rows();
    const int n_cols_dem = filled_dem.n_cols();
    const auto& z_data = filled_dem.get_data_array();
    const auto& z_mask = filled_dem.get_mask_array();

    std::vector<double>       path_z(path.flats.size(), 0.0);
    std::vector<std::uint8_t> path_ok(path.flats.size(), 0);

    for (std::size_t i = 0; i < path.flats.size(); ++i)
    {
        const std::int64_t row = path.flats[i] / n_cols_dem;
        const std::int64_t col = path.flats[i] % n_cols_dem;
        if (row < 0 || row >= n_rows_dem) continue;

        if (z_mask(static_cast<int>(row), static_cast<int>(col)) == 0) continue;
        const double z = static_cast<double>(
            z_data(static_cast<int>(row), static_cast<int>(col)));
        if (!std::isfinite(z)) continue;

        path_z[i]  = z;
        path_ok[i] = 1;
    }

    std::vector<lsdtt_channel_profile_row> rows;
    rows.reserve(sections.size());

    for (const auto& sec : sections)
    {
        lsdtt_channel_profile_row row;
        row.cross_section_id = sec.cross_section_id;

        const auto it = path.index_by_flat.find(sec.centre_flat);
        if (it == path.index_by_flat.end())
        {
            throw std::runtime_error(
                "compute_channel_profile: cross_id " +
                std::to_string(sec.cross_section_id) +
                " is not on the path from source_flat to the outlet. "
                "Was the channel network built from more than one source "
                "point? compute_channel_profile assumes a single mainstem.");
        }
        const std::size_t path_idx = it->second;
        row.distance_from_source_m = path.distance_m[path_idx];

        // Same rule as find_thalweg_index: lowest FINITE valid sample in window.
        double thalweg_z = std::numeric_limits<double>::quiet_NaN();
        for (const auto& sample : sec.samples)
        {
            if (sample.valid == 0) continue;
            if (!std::isfinite(static_cast<double>(sample.elevation))) continue;
            if (std::abs(sample.offset_metres) > thalweg_search_half_width_metres) continue;
            if (std::isnan(thalweg_z) || sample.elevation < thalweg_z)
                thalweg_z = sample.elevation;
        }
        row.thalweg_elevation_m = thalweg_z;

        const std::size_t flat_idx =
            static_cast<std::size_t>(sec.centre_row) * static_cast<std::size_t>(n_cols)
            + static_cast<std::size_t>(sec.centre_col);
        row.drainage_area_m2 = static_cast<double>(da_data.at(flat_idx));

        double slope = 0.0;
        row.slope_valid = fit_along_path_slope(
            path.distance_m, path_z, path_ok, path_idx,
            slope_window_half_length_metres, slope);
        row.slope_m_per_m = row.slope_valid
            ? slope
            : std::numeric_limits<double>::quiet_NaN();

        rows.push_back(row);
    }

    return rows;
}

void write_channel_profile_csv(
    const std::vector<lsdtt_channel_profile_row>& rows,
    const std::string&                            out_path)
{
    std::ofstream out(out_path);
    if (!out.is_open())
        throw std::runtime_error("write_channel_profile_csv: could not open " + out_path);

    out << "cross_id,distance_from_source_m,drainage_area_m2,"
        << "thalweg_elevation_m,slope_m_per_m,slope_valid\n";
    out << std::fixed << std::setprecision(4);

    for (const auto& r : rows)
    {
        // Slope spans orders of magnitude (1e-4 .. 1e-1); format it separately.
        char slope_buf[32];
        std::snprintf(slope_buf, sizeof(slope_buf), "%.8g", r.slope_m_per_m);

        out << r.cross_section_id << ","
            << r.distance_from_source_m << ","
            << r.drainage_area_m2 << ","
            << r.thalweg_elevation_m << ","
            << slope_buf << ","
            << (r.slope_valid ? 1 : 0) << "\n";
    }
}

} // namespace lsdtt