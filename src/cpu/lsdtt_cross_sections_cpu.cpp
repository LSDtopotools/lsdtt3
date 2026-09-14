// lsdtt_cross_sections_cpu.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html

#include "lsdtt_cross_sections_cpu.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "../core/lsdtt_vector_utils.hpp"
#include "../utils/lsdtt_regression.hpp"
#include "../utils/lsdtt_raster_geometry.hpp"

namespace lsdtt
{
namespace
{

constexpr std::array<int, 8> d8_row{{-1, -1,  0,  1,  1,  1,  0, -1}};
constexpr std::array<int, 8> d8_col{{ 0,  1,  1,  1,  0, -1, -1, -1}};

constexpr double k_sqrt_2    = 1.4142135623730950488;
constexpr double k_rad_to_deg = 57.2957795130823208768;

struct sampled_raster_value
{
    float elevation = 0.0f;
    std::uint8_t valid = 0;
};

// -----------------------------------------------------------------------------
// distance_between_flats
// -----------------------------------------------------------------------------
double distance_between_flats(
    const lsdtt_flow_base& flow,
    std::int64_t a,
    std::int64_t b)
{
    const int row_a = flow.row_from_flat_unchecked(a);
    const int col_a = flow.col_from_flat_unchecked(a);
    const int row_b = flow.row_from_flat_unchecked(b);
    const int col_b = flow.col_from_flat_unchecked(b);

    const int drow = std::abs(row_a - row_b);
    const int dcol = std::abs(col_a - col_b);

    if (drow == 0 && dcol == 0)
    {
        return 0.0;
    }

    if (drow > 1 || dcol > 1)
    {
        throw std::runtime_error(
            "lsdtt_cross_sections_cpu: non-local step encountered while "
            "computing along-channel distance"
        );
    }

    return (drow == 1 && dcol == 1)
        ? flow.cell_size() * k_sqrt_2
        : flow.cell_size();
}

// -----------------------------------------------------------------------------
// compass_bearing_degrees_from_xy
//
// Converts a Cartesian vector (dx, dy) into a compass bearing in degrees
// clockwise from north.
// -----------------------------------------------------------------------------
double compass_bearing_degrees_from_xy(double dx, double dy)
{
    double deg = 90.0 - k_rad_to_deg * std::atan2(dy, dx);

    while (deg < 0.0)
    {
        deg += 360.0;
    }
    while (deg >= 360.0)
    {
        deg -= 360.0;
    }

    return deg;
}

// -----------------------------------------------------------------------------
// collect_direct_channel_donors
// -----------------------------------------------------------------------------
void collect_direct_channel_donors(
    const lsdtt_flow_base&             flow,
    const lsdtt_junction_network_base& network,
    std::int64_t                       flat,
    std::vector<std::int64_t>&         donors)
{
    donors.clear();

    const int row = flow.row_from_flat_unchecked(flat);
    const int col = flow.col_from_flat_unchecked(flat);

    for (int dir = 0; dir < 8; ++dir)
    {
        const int nr = row + d8_row[dir];
        const int nc = col + d8_col[dir];

        if (!flow.in_bounds(nr, nc))
        {
            continue;
        }

        const std::int64_t nflat = flow.flat_index_unchecked(nr, nc);

        if (!flow.is_valid_cell(nflat))
        {
            continue;
        }

        if (!network.is_channel(nflat))
        {
            continue;
        }

        if (flow.receiver_of_unchecked(nflat) == flat)
        {
            donors.push_back(nflat);
        }
    }

    std::sort(donors.begin(), donors.end());
}

// -----------------------------------------------------------------------------
// pick_mainstem_donor
//
// At confluences, choose the upslope donor with the largest contributing area,
// using s_vector_size() as the mainstem proxy.
// -----------------------------------------------------------------------------
std::int64_t pick_mainstem_donor(
    const lsdtt_flow_base&             flow,
    const lsdtt_junction_network_base& network,
    std::int64_t                       flat)
{
    std::vector<std::int64_t> donors;
    collect_direct_channel_donors(flow, network, flat, donors);

    if (donors.empty())
    {
        return -1;
    }

    std::int64_t best_donor = donors.front();
    std::uint64_t best_area = flow.s_vector_size(best_donor);

    for (std::size_t i = 1; i < donors.size(); ++i)
    {
        const std::int64_t donor = donors[i];
        const std::uint64_t area = flow.s_vector_size(donor);

        if (area > best_area ||
            (area == best_area && donor < best_donor))
        {
            best_donor = donor;
            best_area  = area;
        }
    }

    return best_donor;
}

// -----------------------------------------------------------------------------
// collect_bearing_window_flats
//
// Builds an along-channel window centred on centre_flat.
// Upstream follows the inferred mainstem donor path.
// Downstream follows the unique channel receiver path.
// -----------------------------------------------------------------------------
std::vector<std::int64_t> collect_bearing_window_flats(
    const lsdtt_flow_base&             flow,
    const lsdtt_junction_network_base& network,
    std::int64_t                       centre_flat,
    double                             window_length_metres)
{
    std::vector<std::int64_t> upstream;
    std::vector<std::int64_t> downstream;

    const double half_window = 0.5 * window_length_metres;

    // -------------------------------------------------------------------------
    // Upstream: follow the mainstem donor path.
    // -------------------------------------------------------------------------
    {
        std::int64_t current = centre_flat;
        double cumulative = 0.0;

        while (cumulative < half_window)
        {
            const std::int64_t donor =
                pick_mainstem_donor(flow, network, current);

            if (donor < 0)
            {
                break;
            }

            cumulative += distance_between_flats(flow, donor, current);
            upstream.push_back(donor);
            current = donor;
        }
    }

    // -------------------------------------------------------------------------
    // Downstream: follow the receiver path.
    // -------------------------------------------------------------------------
    {
        std::int64_t current = centre_flat;
        double cumulative = 0.0;

        while (cumulative < half_window)
        {
            if (flow.is_outlet_unchecked(current))
            {
                break;
            }

            const std::int64_t rec = flow.receiver_of_unchecked(current);

            if (rec < 0 || !network.is_channel(rec))
            {
                break;
            }

            cumulative += distance_between_flats(flow, current, rec);
            downstream.push_back(rec);
            current = rec;
        }
    }

    std::vector<std::int64_t> window_flats;
    window_flats.reserve(upstream.size() + 1 + downstream.size());

    for (std::int64_t i = static_cast<std::int64_t>(upstream.size()) - 1;
         i >= 0; --i)
    {
        window_flats.push_back(upstream[static_cast<std::size_t>(i)]);
    }

    window_flats.push_back(centre_flat);

    for (const std::int64_t flat : downstream)
    {
        window_flats.push_back(flat);
    }

    return window_flats;
}

// -----------------------------------------------------------------------------
// fit_tangent_from_window
//
// Fits an orthogonal linear regression to the channel-pixel centres in the
// local window. The fitted direction is then oriented downstream where possible.
// -----------------------------------------------------------------------------
bool fit_tangent_from_window(
    const lsdtt_raster_meta&         meta,
    const lsdtt_flow_base&           flow,
    const std::vector<std::int64_t>& window_flats,
    std::int64_t                     centre_flat,
    double&                          tangent_dx,
    double&                          tangent_dy)
{
    if (window_flats.size() < 2)
    {
        return false;
    }

    std::vector<double> x;
    std::vector<double> y;
    x.reserve(window_flats.size());
    y.reserve(window_flats.size());

    for (const std::int64_t flat : window_flats)
    {
        const int row = flow.row_from_flat_unchecked(flat);
        const int col = flow.col_from_flat_unchecked(flat);
        const lsdtt_xy xy = cell_center_from_row_col(meta, row, col);

        x.push_back(xy.x);
        y.push_back(xy.y);
    }

    const lsdtt_orthogonal_linear_regression_result fit =
        compute_orthogonal_linear_regression(x, y);

    if (!fit.has_defined_direction)
    {
        return false;
    }

    tangent_dx = fit.direction_x;
    tangent_dy = fit.direction_y;

    // -------------------------------------------------------------------------
    // Orient the tangent downstream if possible, using a local reference vector.
    // -------------------------------------------------------------------------
    std::size_t centre_index = window_flats.size();
    for (std::size_t i = 0; i < window_flats.size(); ++i)
    {
        if (window_flats[i] == centre_flat)
        {
            centre_index = i;
            break;
        }
    }

    if (centre_index == window_flats.size())
    {
        throw std::runtime_error(
            "lsdtt_cross_sections_cpu: centre flat not found in bearing window"
        );
    }

    const int centre_row = flow.row_from_flat_unchecked(centre_flat);
    const int centre_col = flow.col_from_flat_unchecked(centre_flat);
    const lsdtt_xy centre_xy =
        cell_center_from_row_col(meta, centre_row, centre_col);

    double ref_dx = 0.0;
    double ref_dy = 0.0;
    bool   have_reference = false;

    if (centre_index + 1 < window_flats.size())
    {
        const std::int64_t flat_ds = window_flats[centre_index + 1];
        const int row_ds = flow.row_from_flat_unchecked(flat_ds);
        const int col_ds = flow.col_from_flat_unchecked(flat_ds);
        const lsdtt_xy xy_ds =
            cell_center_from_row_col(meta, row_ds, col_ds);

        ref_dx = xy_ds.x - centre_xy.x;
        ref_dy = xy_ds.y - centre_xy.y;
        have_reference = true;
    }
    else if (centre_index > 0)
    {
        const std::int64_t flat_us = window_flats[centre_index - 1];
        const int row_us = flow.row_from_flat_unchecked(flat_us);
        const int col_us = flow.col_from_flat_unchecked(flat_us);
        const lsdtt_xy xy_us =
            cell_center_from_row_col(meta, row_us, col_us);

        ref_dx = centre_xy.x - xy_us.x;
        ref_dy = centre_xy.y - xy_us.y;
        have_reference = true;
    }

    if (have_reference)
    {
        const double dot = tangent_dx * ref_dx + tangent_dy * ref_dy;

        if (dot < 0.0)
        {
            tangent_dx = -tangent_dx;
            tangent_dy = -tangent_dy;
        }
    }

    return true;
}

// -----------------------------------------------------------------------------
// bilinear_sample_raster
//
// Generic bilinear sampling for DEM elevation or any other float raster.
// Returns invalid if sampling support contains nodata or falls outside the raster.
// -----------------------------------------------------------------------------
sampled_raster_value bilinear_sample_raster(
    const lsdtt_raster_cpu& raster,
    double                  x,
    double                  y)
{
    const lsdtt_raster_meta& meta = raster.get_meta();
    const lsdtt_array2d_f& data = raster.get_data_array();
    const lsdtt_array2d_u8& mask = raster.get_mask_array();

    sampled_raster_value out;
    out.elevation = static_cast<float>(meta.no_data_value);
    out.valid = 0;

    if (x < meta.x_min || x > meta.x_max() ||
        y < meta.y_min || y > meta.y_max())
    {
        return out;
    }

    const double col_f =
        (x - meta.x_min) / meta.cell_size - 0.5;

    const double row_f =
        (meta.y_max() - y) / meta.cell_size - 0.5;

    if (col_f < 0.0 || row_f < 0.0 ||
        col_f > static_cast<double>(meta.n_cols - 1) ||
        row_f > static_cast<double>(meta.n_rows - 1))
    {
        return out;
    }

    const int col0 = static_cast<int>(std::floor(col_f));
    const int row0 = static_cast<int>(std::floor(row_f));

    // Edge fallback: nearest-neighbour if a full 2x2 bilinear stencil
    // is not available.
    if (col0 >= meta.n_cols - 1 || row0 >= meta.n_rows - 1)
    {
        const int col_nn = static_cast<int>(std::llround(col_f));
        const int row_nn = static_cast<int>(std::llround(row_f));

        if (row_nn < 0 || row_nn >= meta.n_rows ||
            col_nn < 0 || col_nn >= meta.n_cols)
        {
            return out;
        }

        if (mask(row_nn, col_nn) == 0)
        {
            return out;
        }

        out.elevation = data(row_nn, col_nn);
        out.valid = 1;
        return out;
    }

    const int row1 = row0 + 1;
    const int col1 = col0 + 1;

    if (mask(row0, col0) == 0 ||
        mask(row0, col1) == 0 ||
        mask(row1, col0) == 0 ||
        mask(row1, col1) == 0)
    {
        return out;
    }

    const double dc = col_f - static_cast<double>(col0);
    const double dr = row_f - static_cast<double>(row0);

    const double z00 = static_cast<double>(data(row0, col0));
    const double z01 = static_cast<double>(data(row0, col1));
    const double z10 = static_cast<double>(data(row1, col0));
    const double z11 = static_cast<double>(data(row1, col1));

    const double z0 = (1.0 - dc) * z00 + dc * z01;
    const double z1 = (1.0 - dc) * z10 + dc * z11;
    const double z  = (1.0 - dr) * z0  + dr * z1;

    out.elevation = static_cast<float>(z);
    out.valid = 1;
    return out;
}

// -----------------------------------------------------------------------------
// build_one_cross_section
// -----------------------------------------------------------------------------
bool build_one_cross_section(
    const lsdtt_raster_cpu&            dem,
    const lsdtt_raster_cpu*            vertical_excess_curvature_raster,
    const lsdtt_flow_base&             flow,
    const lsdtt_junction_network_base& network,
    std::int64_t                       centre_flat,
    std::int64_t                       cross_section_id,
    const lsdtt_cross_section_options& options,
    lsdtt_cross_section&               section)
{
    const lsdtt_raster_meta& meta = dem.get_meta();

    std::vector<std::int64_t> window_flats =
        collect_bearing_window_flats(
            flow,
            network,
            centre_flat,
            options.bearing_window_length_metres
        );

    if (static_cast<int>(window_flats.size()) <
        options.minimum_bearing_window_pixels)
    {
        return false;
    }

    double tangent_dx = 0.0;
    double tangent_dy = 0.0;

    if (!fit_tangent_from_window(
            meta,
            flow,
            window_flats,
            centre_flat,
            tangent_dx,
            tangent_dy))
    {
        return false;
    }

    const double tangent_norm =
        std::sqrt(tangent_dx * tangent_dx + tangent_dy * tangent_dy);

    if (!(tangent_norm > 0.0))
    {
        return false;
    }

    tangent_dx /= tangent_norm;
    tangent_dy /= tangent_norm;

    const double normal_dx = -tangent_dy;
    const double normal_dy =  tangent_dx;

    const int centre_row = flow.row_from_flat_unchecked(centre_flat);
    const int centre_col = flow.col_from_flat_unchecked(centre_flat);
    const lsdtt_xy centre_xy =
        cell_center_from_row_col(meta, centre_row, centre_col);

    section = lsdtt_cross_section{};
    section.cross_section_id = cross_section_id;
    section.centre_flat      = centre_flat;
    section.centre_row       = centre_row;
    section.centre_col       = centre_col;
    section.centre_x         = centre_xy.x;
    section.centre_y         = centre_xy.y;

    section.tangent_dx = tangent_dx;
    section.tangent_dy = tangent_dy;
    section.normal_dx  = normal_dx;
    section.normal_dy  = normal_dy;

    section.tangent_bearing_degrees =
        compass_bearing_degrees_from_xy(tangent_dx, tangent_dy);

    section.normal_bearing_degrees =
        compass_bearing_degrees_from_xy(normal_dx, normal_dy);

    section.bearing_window_length_metres =
        options.bearing_window_length_metres;

    section.half_length_metres =
        options.cross_section_half_length_metres;

    section.sample_spacing_metres =
        options.cross_section_sample_spacing_metres;

    section.start_xy =
        lsdtt_xy{
            centre_xy.x - section.half_length_metres * normal_dx,
            centre_xy.y - section.half_length_metres * normal_dy
        };

    section.end_xy =
        lsdtt_xy{
            centre_xy.x + section.half_length_metres * normal_dx,
            centre_xy.y + section.half_length_metres * normal_dy
        };

    section.samples.clear();

    std::int32_t sample_index = 0;
    const double half_length = options.cross_section_half_length_metres;
    const double spacing     = options.cross_section_sample_spacing_metres;

    for (double offset = -half_length;
         offset <= half_length + 1.0e-9;
         offset += spacing)
    {
        lsdtt_cross_section_sample sample;
        sample.sample_index  = sample_index++;
        sample.offset_metres = offset;
        sample.x = centre_xy.x + offset * normal_dx;
        sample.y = centre_xy.y + offset * normal_dy;

        // DEM elevation
        const sampled_raster_value z =
            bilinear_sample_raster(dem, sample.x, sample.y);

        sample.elevation = z.elevation;
        sample.valid     = z.valid;

        // Optional vertical excess curvature
        if (vertical_excess_curvature_raster != nullptr)
        {
            const sampled_raster_value v =
                bilinear_sample_raster(
                    *vertical_excess_curvature_raster,
                    sample.x,
                    sample.y
                );

            sample.vertical_excess_curvature = v.elevation;
            sample.vertical_excess_curvature_valid = v.valid;
        }
        else
        {
            sample.vertical_excess_curvature =
                std::numeric_limits<float>::quiet_NaN();
            sample.vertical_excess_curvature_valid = 0;
        }

        section.samples.push_back(sample);
    }

    return true;
}

} // unnamed namespace

// =============================================================================
// extract_cross_sections
// =============================================================================
std::vector<lsdtt_cross_section> extract_cross_sections(
    const lsdtt_raster_cpu&            dem,
    const lsdtt_raster_cpu*            vertical_excess_curvature_raster,
    const lsdtt_flow_base&             flow,
    const lsdtt_junction_network_base& network,
    const lsdtt_cross_section_options& options)
{
    if (flow.empty())
    {
        throw std::runtime_error(
            "extract_cross_sections: flow object is empty"
        );
    }

    if (network.empty())
    {
        throw std::runtime_error(
            "extract_cross_sections: junction network is empty"
        );
    }

    if (options.cross_section_half_length_metres <= 0.0)
    {
        throw std::runtime_error(
            "extract_cross_sections: cross_section_half_length_metres must be > 0"
        );
    }

    if (options.cross_section_sample_spacing_metres <= 0.0)
    {
        throw std::runtime_error(
            "extract_cross_sections: cross_section_sample_spacing_metres must be > 0"
        );
    }

    if (options.bearing_window_length_metres <= 0.0)
    {
        throw std::runtime_error(
            "extract_cross_sections: bearing_window_length_metres must be > 0"
        );
    }

    if (options.cross_section_thinning <= 0)
    {
        throw std::runtime_error(
            "extract_cross_sections: cross_section_thinning must be > 0"
        );
    }

    if (options.minimum_bearing_window_pixels < 2)
    {
        throw std::runtime_error(
            "extract_cross_sections: minimum_bearing_window_pixels must be >= 2"
        );
    }

    if (vertical_excess_curvature_raster != nullptr)
    {
        const lsdtt_raster_meta& vec_meta =
            vertical_excess_curvature_raster->get_meta();

        const lsdtt_raster_meta& dem_meta = dem.get_meta();

        if (vec_meta.n_rows    != dem_meta.n_rows ||
            vec_meta.n_cols    != dem_meta.n_cols ||
            vec_meta.cell_size != dem_meta.cell_size ||
            vec_meta.x_min     != dem_meta.x_min ||
            vec_meta.y_min     != dem_meta.y_min)
        {
            throw std::runtime_error(
                "extract_cross_sections: vertical_excess_curvature_raster is "
                "not element-wise compatible with the DEM"
            );
        }
    }

    std::vector<lsdtt_cross_section> sections;
    std::int64_t next_id = 0;
    std::int64_t channel_counter = 0;

    for (int row = 0; row < network.n_rows(); ++row)
    {
        for (int col = 0; col < network.n_cols(); ++col)
        {
            const std::int64_t flat =
                flow.flat_index_unchecked(row, col);

            if (!network.is_channel(flat))
            {
                continue;
            }

            if (options.skip_junction_pixels &&
                network.is_junction_flat(flat))
            {
                continue;
            }

            if ((channel_counter % options.cross_section_thinning) != 0)
            {
                ++channel_counter;
                continue;
            }

            ++channel_counter;

            lsdtt_cross_section section;
            const bool ok =
                build_one_cross_section(
                    dem,
                    vertical_excess_curvature_raster,
                    flow,
                    network,
                    flat,
                    next_id,
                    options,
                    section
                );

            if (ok)
            {
                sections.push_back(std::move(section));
                ++next_id;
            }
        }
    }

    return sections;
}

// =============================================================================
// write_cross_section_lines
//
// Simple geometry output: line + cross_id only.
// =============================================================================
void write_cross_section_lines(
    const std::vector<lsdtt_cross_section>& sections,
    const lsdtt_raster_meta&                meta,
    const std::string&                      out_path,
    vector_output_format                    format,
    bool                                    spatial_index,
    const std::string&                      layer_name)
{
    lsdtt_vector_schema schema;
    schema.field_names = {
        "cross_id"
    };
    schema.field_types = {
        lsdtt_vector_field_type::int64_type
    };

    lsdtt_data_loader loader;
    loader.begin_vector_write(
        out_path,
        layer_name,
        schema,
        meta.crs_wkt,
        format,
        spatial_index
    );

    for (const auto& section : sections)
    {
        loader.write_line(
            {section.start_xy, section.end_xy},
            {
                section.cross_section_id
            }
        );
    }

    loader.finish_vector_write();
}

// =============================================================================
// write_cross_section_sample_points
// =============================================================================
void write_cross_section_sample_points(
    const std::vector<lsdtt_cross_section>& sections,
    const lsdtt_raster_meta&                meta,
    const std::string&                      out_path,
    point_output_format                     format,
    bool                                    spatial_index,
    const std::string&                      layer_name)
{
    lsdtt_data_loader loader;

    // -------------------------------------------------------------------------
    // CSV output in EPSG:4326
    // -------------------------------------------------------------------------
    if (format == point_output_format::csv)
    {
        if (meta.crs_wkt.empty())
        {
            throw std::runtime_error(
                "write_cross_section_sample_points: CSV point output requires "
                "a valid raster CRS so points can be transformed to EPSG:4326"
            );
        }

        const lsdtt_coordinate_transform transform =
            lsdtt_coordinate_transform::from_epsg(meta.crs_wkt, 4326);

        std::vector<std::string> headers = {
            "latitude",
            "longitude",
            "cross_id",
            "centre_flat",
            "sample_index",
            "offset_m",
            "elevation",
            "valid",
            "vertical_excess_curvature",
            "vec_valid"
        };

        auto precise_string = [](double value) -> std::string
        {
            std::ostringstream ss;
            ss.precision(17);
            ss << value;
            return ss.str();
        };

        std::vector<std::vector<std::string>> rows;

        for (const auto& section : sections)
        {
            for (const auto& sample : section.samples)
            {
                lsdtt_lonlat ll =
                    xy_to_lonlat(sample.x, sample.y, transform);

                rows.push_back({
                    precise_string(ll.latitude),
                    precise_string(ll.longitude),
                    std::to_string(section.cross_section_id),
                    std::to_string(section.centre_flat),
                    std::to_string(sample.sample_index),
                    precise_string(sample.offset_metres),
                    precise_string(static_cast<double>(sample.elevation)),
                    std::to_string(static_cast<int>(sample.valid)),
                    precise_string(
                        static_cast<double>(sample.vertical_excess_curvature)
                    ),
                    std::to_string(
                        static_cast<int>(sample.vertical_excess_curvature_valid)
                    )
                });
            }
        }

        loader.write_csv(out_path, headers, rows);
        return;
    }

    // -------------------------------------------------------------------------
    // Georeferenced point output
    // -------------------------------------------------------------------------
    lsdtt_vector_schema schema;
    schema.field_names = {
        "cross_id",
        "centre_flat",
        "sample_index",
        "offset_m",
        "elevation",
        "valid",
        "vertical_excess_curvature",
        "vec_valid"
    };
    schema.field_types = {
        lsdtt_vector_field_type::int64_type,
        lsdtt_vector_field_type::int64_type,
        lsdtt_vector_field_type::int64_type,
        lsdtt_vector_field_type::double_type,
        lsdtt_vector_field_type::double_type,
        lsdtt_vector_field_type::int64_type,
        lsdtt_vector_field_type::double_type,
        lsdtt_vector_field_type::int64_type
    };

    const vector_output_format out_format =
        convert_point_to_vector_output_format(format);

    loader.begin_vector_write(
        out_path,
        layer_name,
        schema,
        meta.crs_wkt,
        out_format,
        spatial_index
    );

    for (const auto& section : sections)
    {
        for (const auto& sample : section.samples)
        {
            loader.write_point(
                sample.x,
                sample.y,
                {
                    section.cross_section_id,
                    section.centre_flat,
                    static_cast<std::int64_t>(sample.sample_index),
                    sample.offset_metres,
                    static_cast<double>(sample.elevation),
                    static_cast<std::int64_t>(sample.valid),
                    static_cast<double>(sample.vertical_excess_curvature),
                    static_cast<std::int64_t>(
                        sample.vertical_excess_curvature_valid
                    )
                }
            );
        }
    }

    loader.finish_vector_write();
}

// =============================================================================
// write_cross_section_samples_csv
// =============================================================================
void write_cross_section_samples_csv(
    const std::vector<lsdtt_cross_section>& sections,
    const std::string&                      out_path)
{
    lsdtt_data_loader loader;

    std::vector<std::string> headers = {
        "cross_id",
        "centre_flat",
        "centre_row",
        "centre_col",
        "centre_x",
        "centre_y",
        "sample_index",
        "offset_m",
        "sample_x",
        "sample_y",
        "elevation",
        "valid",
        "vertical_excess_curvature",
        "vec_valid",
        "tangent_bearing_deg",
        "normal_bearing_deg"
    };

    auto precise_string = [](double value) -> std::string
    {
        std::ostringstream ss;
        ss.precision(17);
        ss << value;
        return ss.str();
    };

    std::vector<std::vector<std::string>> rows;

    for (const auto& section : sections)
    {
        for (const auto& sample : section.samples)
        {
            rows.push_back({
                std::to_string(section.cross_section_id),
                std::to_string(section.centre_flat),
                std::to_string(section.centre_row),
                std::to_string(section.centre_col),
                precise_string(section.centre_x),
                precise_string(section.centre_y),
                std::to_string(sample.sample_index),
                precise_string(sample.offset_metres),
                precise_string(sample.x),
                precise_string(sample.y),
                precise_string(static_cast<double>(sample.elevation)),
                std::to_string(static_cast<int>(sample.valid)),
                precise_string(
                    static_cast<double>(sample.vertical_excess_curvature)
                ),
                std::to_string(
                    static_cast<int>(sample.vertical_excess_curvature_valid)
                ),
                precise_string(section.tangent_bearing_degrees),
                precise_string(section.normal_bearing_degrees)
            });
        }
    }

    loader.write_csv(out_path, headers, rows);
}

} // namespace lsdtt
