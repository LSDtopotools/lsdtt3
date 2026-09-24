// lsdtt_manning_hydraulics_cpu.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html

#include "lsdtt_manning_hydraulics_cpu.hpp"
#include "lsdtt_raster_surface_metrics_cpu.hpp"
#include "../core/lsdtt_data_loader.hpp"
#include "../core/lsdtt_vector_utils.hpp"
#include "../core/lsdtt_flow_base.hpp"
#include "../core/lsdtt_types.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <unordered_map>

namespace lsdtt
{
namespace
{

std::vector<std::string> split_csv_line(const std::string& line)
{
    std::vector<std::string> out;
    std::stringstream ss(line);
    std::string field;
    while (std::getline(ss, field, ','))
    {
        while (!field.empty() && (field.back() == '\r' || field.back() == ' '))
            field.pop_back();
        out.push_back(field);
    }
    return out;
}

int find_col(const std::vector<std::string>& header, const std::string& name)
{
    for (std::size_t i = 0; i < header.size(); ++i)
        if (header[i] == name) return static_cast<int>(i);
    throw std::runtime_error(
        "lsdtt_manning_hydraulics_cpu: missing required column: " + name);
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
    for (std::size_t i = 0; i < order.size(); ++i) order[i] = i;
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

std::size_t find_thalweg_index(const std::vector<double>& z)
{
    return static_cast<std::size_t>(
        std::distance(z.begin(), std::min_element(z.begin(), z.end())));
}

// Finds the first crossing of the water-surface elevation moving outward
// from the thalweg, interpolating both the along-section offset (s) and
// the real raster-CRS coordinates (x, y) between the two straddling
// samples. direction: +1 = search right (increasing s), -1 = search left.
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

// Walks downstream from source_flat via the D8 receiver chain,
// accumulating true along-channel distance (cardinal steps = cell_size,
// diagonal steps = cell_size * sqrt(2)), recording the cumulative
// distance at every flat visited along the way.
//
// Direction code convention (from lsdtt_flowinfo_cpu):
//   7  0  1
//   6  x  2
//   5  4  3
// Odd codes (1,3,5,7) are diagonal; even codes (0,2,4,6) are cardinal.
//
// A step cap guards against an unexpected cycle in a corrupted graph
// rather than looping forever.
std::unordered_map<std::int64_t, double> build_distance_from_source_map(
    const lsdtt_flow_base& flow,
    std::int64_t            source_flat)
{
    std::unordered_map<std::int64_t, double> distance_by_flat;

    const double cell_step = flow.cell_size();
    const double diag_step = cell_step * std::sqrt(2.0);

    double distance = 0.0;
    std::int64_t current = source_flat;
    const std::int64_t max_steps = flow.n_valid_cells() + 1;
    std::int64_t steps = 0;

    distance_by_flat[current] = 0.0;

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

        distance_by_flat[receiver] = distance;
        current = receiver;
        ++steps;

        if (steps > max_steps)
        {
            // Should not happen on an acyclic D8 graph -- flags a corrupt
            // network rather than hanging.
            break;
        }
    }

    return distance_by_flat;
}

} // unnamed namespace

// -----------------------------------------------------------------------------
// read_manning_reaches_csv
//
// Header: start_distance_m,end_distance_m,manning_n,target_discharge_m3s
//
// Reach bounds are along-channel distance from source, in metres. This
// is deliberate: cross_section_id is assigned during extraction by a
// row-major scan of the raster and has no relationship to position
// along a (typically sinuous) channel, so a cross_id range cannot
// correctly select "every cross section between these two points along
// the channel" -- distance can.
// -----------------------------------------------------------------------------
std::vector<lsdtt_manning_reach> read_manning_reaches_csv(
    const std::string& csv_path)
{
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
    while (std::getline(in, line))
    {
        if (line.empty()) continue;
        const std::vector<std::string> t = split_csv_line(line);
        if (t.size() < header.size()) continue;

        lsdtt_manning_reach r;
        r.start_distance_m     = std::stod(t[c_start]);
        r.end_distance_m       = std::stod(t[c_end]);
        r.manning_n            = std::stod(t[c_n]);
        r.target_discharge_m3s = std::stod(t[c_q]);
        reaches.push_back(r);
    }
    return reaches;
}

// -----------------------------------------------------------------------------
// read_channel_profile_csv
//
// Loads per-cross-section distance, slope, and drainage area from the
// CSV written by write_channel_profile_csv. distance_by_cross_id is
// required for matching cross sections to reaches; a cross_id whose
// slope_valid column is 0 is deliberately left OUT of
// slope_by_cross_id, so the Manning solve can flag it
// ("no_slope_in_channel_profile_csv") rather than silently using a bad
// or zero slope value.
// -----------------------------------------------------------------------------
lsdtt_channel_profile_lookup read_channel_profile_csv(
    const std::string& csv_path)
{
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
    while (std::getline(in, line))
    {
        if (line.empty()) continue;
        const std::vector<std::string> t = split_csv_line(line);
        if (t.size() < header.size()) continue;

        const long long cross_id = std::stoll(t[c_id]);
        lookup.distance_by_cross_id[cross_id]      = std::stod(t[c_dist]);
        lookup.drainage_area_by_cross_id[cross_id] = std::stod(t[c_area]);

        const bool slope_ok = std::stoi(t[c_valid]) != 0;
        if (slope_ok)
        {
            lookup.slope_by_cross_id[cross_id] = std::stod(t[c_slope]);
        }
    }
    return lookup;
}

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

double manning_discharge(
    double area_m2, double hydraulic_radius_m,
    double slope, double manning_n)
{
    if (manning_n <= 0.0 || slope <= 0.0 || area_m2 <= 0.0) return 0.0;
    return (1.0 / manning_n) * area_m2 *
           std::pow(hydraulic_radius_m, 2.0 / 3.0) * std::sqrt(slope);
}

bool solve_bankfull_depth(
    const std::vector<double>& s,
    const std::vector<double>& z,
    double thalweg_elevation,
    double slope,
    double manning_n,
    double target_discharge_m3s,
    const manning_solver_options& opts,
    double& solved_depth_m)
{
    auto q_at_depth = [&](double depth) -> double
    {
        const lsdtt_wetted_geometry g =
            compute_wetted_geometry(s, z, thalweg_elevation + depth);
        if (!g.valid) return 0.0;
        return manning_discharge(g.area_m2, g.hydraulic_radius_m, slope, manning_n);
    };

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

// -----------------------------------------------------------------------------
// compute_manning_widths_for_reach
//
// Cross sections are matched to the reach by DISTANCE FROM SOURCE
// (looked up per cross section from `profile`), not by cross_id. This
// is the key correctness fix: cross_id is assigned during extraction by
// raster row-major scan order and is unrelated to position along a
// (generally sinuous) channel, so a cross_id numeric range cannot
// reliably select a contiguous stretch of channel -- only distance can.
//
// reach.manning_n and reach.target_discharge_m3s are held fixed across
// every matched cross section. Slope is looked up per cross section
// from `profile` -- depth is solved independently at each XS using that
// XS's own slope and geometry, not a reach-averaged slope, so Q is
// constant across the reach but solved depth varies with local slope
// and channel shape.
// -----------------------------------------------------------------------------
std::vector<lsdtt_manning_width_record> compute_manning_widths_for_reach(
    const std::vector<lsdtt_cross_section>& sections,
    const lsdtt_manning_reach&              reach,
    const lsdtt_channel_profile_lookup&     profile,
    const manning_solver_options&           opts)
{
    std::vector<lsdtt_manning_width_record> out;

    for (const auto& sec : sections)
    {
        // --- Distance lookup determines reach membership. ---
        const auto dist_it = profile.distance_by_cross_id.find(sec.cross_section_id);
        if (dist_it == profile.distance_by_cross_id.end())
        {
            // Not present in the channel-profile CSV at all -- cannot be
            // assigned to any reach without a distance value. Silently
            // skip rather than emit a spurious record for it.
            continue;
        }
        const double distance_m = dist_it->second;

        if (distance_m < reach.start_distance_m ||
            distance_m > reach.end_distance_m)
            continue;

        lsdtt_manning_width_record rec;
        rec.cross_section_id = sec.cross_section_id;
        rec.centre_x = sec.centre_x;
        rec.centre_y = sec.centre_y;
        rec.distance_from_source_m = distance_m;

        // --- Per-XS slope lookup. ---
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
        {
            rec.drainage_area_m2 = area_it->second;
        }

        std::vector<double> s, z, x, y;
        if (!extract_and_sort_profile(sec, s, z, x, y))
        {
            rec.quality_flag = "invalid_profile";
            out.push_back(rec);
            continue;
        }

        const std::size_t th_idx = find_thalweg_index(z);
        rec.thalweg_s = s[th_idx];
        rec.thalweg_z = z[th_idx];

        double depth = 0.0;
        if (!solve_bankfull_depth(
                s, z, rec.thalweg_z, rec.slope_m_per_m, reach.manning_n,
                reach.target_discharge_m3s, opts, depth))
        {
            rec.quality_flag = "target_discharge_not_bracketed";
            out.push_back(rec);
            continue;
        }

        rec.solved_depth_m = depth;
        rec.water_surface_elevation_m = rec.thalweg_z + depth;

        const lsdtt_wetted_geometry g =
            compute_wetted_geometry(s, z, rec.water_surface_elevation_m);
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
    const manning_solver_options&           opts)
{
    std::vector<lsdtt_manning_width_record> all;
    for (const auto& reach : reaches)
    {
        auto reach_results = compute_manning_widths_for_reach(sections, reach, profile, opts);
        all.insert(all.end(),
            std::make_move_iterator(reach_results.begin()),
            std::make_move_iterator(reach_results.end()));
    }
    // Sort by distance-from-source, not cross_id, so output rows read
    // in sensible along-channel order.
    std::sort(all.begin(), all.end(),
        [](const lsdtt_manning_width_record& a, const lsdtt_manning_width_record& b)
        { return a.distance_from_source_m < b.distance_from_source_m; });
    return all;
}

// -----------------------------------------------------------------------------
void write_manning_bank_points(
    const std::vector<lsdtt_manning_width_record>& results,
    const lsdtt_raster_meta&                        meta,
    const std::string&                              out_path,
    point_output_format                             format,
    bool                                             spatial_index,
    const std::string&                               layer_name)
{
    lsdtt_data_loader loader;

    if (format == point_output_format::csv)
    {
        write_manning_widths_csv(results, out_path);
        return;
    }

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
        out << r.cross_section_id << "," << r.distance_from_source_m << ","
            << r.centre_x << "," << r.centre_y << ","
            << r.slope_m_per_m << "," << r.drainage_area_m2 << ","
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
std::vector<lsdtt_channel_profile_row> compute_channel_profile(
    const std::vector<lsdtt_cross_section>& sections,
    const lsdtt_flow_base&                  flow,
    const lsdtt_raster_cpu&                 dem,
    float                                    slope_window_radius_metres,
    std::int64_t                             source_flat,
    int                                       n_threads)
{
    // --- Distance from source: single downstream walk, O(path length) ---
    const std::unordered_map<std::int64_t, double> distance_by_flat =
        build_distance_from_source_map(flow, source_flat);

    // --- Drainage area raster (D8, matches the flow object already built) ---
    const lsdtt_raster_result drainage_area_result = flow.get_drainage_area_raster();
    const auto& da_data = drainage_area_result.data;
    const int   n_cols  = drainage_area_result.meta.n_cols;

    // --- Slope raster via direct_polyfit, sampled at each XS centre pixel ---
    poly_fit_coeffs coeffs = direct_polyfit(dem, slope_window_radius_metres, n_threads);

    polyfit_metric_request slope_req;
    slope_req.metrics   = { polyfit_metric::slope };
    slope_req.n_threads = n_threads;

    polyfit_metric_result slope_res = compute_polyfit_metrics(coeffs, slope_req);
    const lsdtt_raster_cpu& slope_raster = slope_res.at(polyfit_metric::slope);
    const auto& slope_data = slope_raster.get_data_array();
    const auto& slope_mask = slope_raster.get_mask_array();

    // --- Assemble one row per cross section ---
    std::vector<lsdtt_channel_profile_row> rows;
    rows.reserve(sections.size());

    for (const auto& sec : sections)
    {
        lsdtt_channel_profile_row row;
        row.cross_section_id = sec.cross_section_id;

        const auto it = distance_by_flat.find(sec.centre_flat);
        if (it == distance_by_flat.end())
        {
            throw std::runtime_error(
                "compute_channel_profile: cross_id " +
                std::to_string(sec.cross_section_id) +
                " is not on the path from source_flat to the outlet. "
                "Was the channel network built from more than one source "
                "point? compute_channel_profile assumes a single mainstem.");
        }
        row.distance_from_source_m = it->second;

        double thalweg_z = std::numeric_limits<double>::quiet_NaN();
        for (const auto& sample : sec.samples)
        {
            if (sample.valid == 0) continue;
            if (std::isnan(thalweg_z) || sample.elevation < thalweg_z)
                thalweg_z = sample.elevation;
        }
        row.thalweg_elevation_m = thalweg_z;

        const std::size_t flat_idx =
            static_cast<std::size_t>(sec.centre_row) * static_cast<std::size_t>(n_cols)
            + static_cast<std::size_t>(sec.centre_col);
        row.drainage_area_m2 = static_cast<double>(da_data.at(flat_idx));

        const bool slope_ok = slope_mask(sec.centre_row, sec.centre_col) != 0;
        row.slope_valid = slope_ok;
        row.slope_m_per_m = slope_ok
            ? static_cast<double>(slope_data(sec.centre_row, sec.centre_col))
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
        out << r.cross_section_id << ","
            << r.distance_from_source_m << ","
            << r.drainage_area_m2 << ","
            << r.thalweg_elevation_m << ","
            << r.slope_m_per_m << ","
            << (r.slope_valid ? 1 : 0) << "\n";
    }
}

} // namespace lsdtt