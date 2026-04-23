// lsdtt_burn_streams_cpu.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html

#include "lsdtt_burn_streams_cpu.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "../core/lsdtt_data_loader.hpp"
#include "../core/lsdtt_raster_meta.hpp"
#include "lsdtt_proximity_cpu.hpp"

namespace lsdtt
{

// ---------------------------------------------------------------------------
// resolve_stream_burn_mode
// ---------------------------------------------------------------------------
stream_burn_mode resolve_stream_burn_mode(const std::string& s)
{
    if (s == "burn_in_place")    return stream_burn_mode::burn_in_place;
    if (s == "snap_to_network")  return stream_burn_mode::snap_to_network;
    std::cerr << "Warning: unknown stream_burn_mode '" << s
              << "'; defaulting to burn_in_place.\n";
    return stream_burn_mode::burn_in_place;
}

// ---------------------------------------------------------------------------
// burn_streams_in_place — internal helpers
// ---------------------------------------------------------------------------

namespace {

// Apply a running-minimum scan from the upstream end (last element, highest
// flow_distance) toward the outlet (first element, flow_distance = 0).
// After this call every node satisfies:
//   burned_elevation[i] <= burned_elevation[i+1]   (non-increasing downstream)
//   burned_elevation[i] <= dem_elevation[i]         (never raised above DEM)
void apply_downhill_constraint(std::vector<burned_stream_node>& nodes)
{
    if (nodes.empty()) return;

    // Nodes are sorted ascending by flow_distance (index 0 = outlet).
    // Scan from upstream (n-1) toward outlet (0), tracking the running minimum.
    double running_min = nodes.back().dem_elevation;
    for (std::size_t i = nodes.size(); i-- > 0; )
    {
        running_min = std::min(nodes[i].dem_elevation, running_min);
        nodes[i].burned_elevation = running_min;
    }
}

} // anonymous namespace

// ---------------------------------------------------------------------------
// burn_streams_in_place
// ---------------------------------------------------------------------------
std::vector<burned_stream_profile> burn_streams_in_place(
    const lsdtt_raster_cpu& dem,
    const std::string&      vector_path,
    int                     n_threads
)
{
    // -----------------------------------------------------------------------
    // Load vector features
    // -----------------------------------------------------------------------
    lsdtt_data_loader vec_loader;
    try
    {
        vec_loader.open_vector(vector_path);
    }
    catch (const std::exception& e)
    {
        throw std::runtime_error(
            std::string("burn_streams_in_place: cannot open vector file: ")
            + e.what());
    }

    const auto layer_names = vec_loader.get_layer_names();
    if (layer_names.empty())
        throw std::runtime_error(
            "burn_streams_in_place: vector file has no layers: " + vector_path);

    const lsdtt_vector_layer vlayer =
        vec_loader.read_vector_layer(layer_names[0]);
    vec_loader.close_vector();

    if (vlayer.empty())
        throw std::runtime_error(
            "burn_streams_in_place: vector layer is empty: " + vector_path);

    // -----------------------------------------------------------------------
    // DEM metadata
    // -----------------------------------------------------------------------
    const lsdtt_raster_meta& meta      = dem.get_meta();
    const int                n_rows    = meta.n_rows;
    const int                n_cols    = meta.n_cols;
    const double             cell_size = meta.cell_size;

    const std::vector<float>&   dem_data = dem.get_data();
    const std::vector<uint8_t>& dem_mask = dem.get_mask();

    // A pixel is considered "on the channel" when the perpendicular distance
    // to the polyline is within half the pixel diagonal.
    const float perp_threshold =
        static_cast<float>(std::sqrt(2.0) * 0.5 * cell_size);

    // -----------------------------------------------------------------------
    // Process each polyline feature independently
    // -----------------------------------------------------------------------
    std::vector<burned_stream_profile> profiles;
    profiles.reserve(vlayer.size());

    for (std::size_t fi = 0; fi < vlayer.size(); ++fi)
    {
        const auto& record = vlayer.records[fi];

        // Extract polyline vertices
        const lsdtt_line_geom* line_ptr =
            std::get_if<lsdtt_line_geom>(&record.geometry);
        if (!line_ptr)
        {
            std::cerr << "Warning: feature " << fi
                      << " is not a polyline; skipping.\n";
            continue;
        }
        if (line_ptr->points.size() < 2)
        {
            std::cerr << "Warning: feature " << fi
                      << " has fewer than 2 vertices; skipping.\n";
            continue;
        }

        // -------------------------------------------------------------------
        // Swath distances: along-axis and perpendicular for every DEM pixel
        // -------------------------------------------------------------------
        swath_distances swath;
        try
        {
            swath = compute_swath_distances(dem, line_ptr->points, n_threads);
        }
        catch (const std::exception& e)
        {
            std::cerr << "Warning: feature " << fi
                      << ": swath computation failed (" << e.what()
                      << "); skipping.\n";
            continue;
        }

        const std::vector<float>&   along_data = swath.along_axis.get_data();
        const std::vector<uint8_t>& along_mask = swath.along_axis.get_mask();
        const std::vector<float>&   perp_data  = swath.perp_distance.get_data();

        // -------------------------------------------------------------------
        // Collect pixels beneath the polyline
        // -------------------------------------------------------------------
        struct raw_node
        {
            double along_axis;
            double dem_elev;
            double x, y;
        };
        std::vector<raw_node> raw_nodes;

        for (int row = 0; row < n_rows; ++row)
        {
            for (int col = 0; col < n_cols; ++col)
            {
                const auto idx = static_cast<std::size_t>(row * n_cols + col);
                if (!dem_mask[idx])   continue;  // nodata
                if (!along_mask[idx]) continue;  // swath nodata (shouldn't happen)
                if (perp_data[idx] > perp_threshold) continue;  // not on channel

                const double x =
                    meta.x_min + (col + 0.5) * cell_size;
                const double y =
                    meta.y_min + (n_rows - row - 0.5) * cell_size;

                raw_nodes.push_back({
                    static_cast<double>(along_data[idx]),
                    static_cast<double>(dem_data[idx]),
                    x, y
                });
            }
        }

        if (raw_nodes.empty())
        {
            std::cerr << "Warning: feature " << fi
                      << ": no valid DEM pixels found beneath the polyline; "
                         "skipping.\n";
            continue;
        }

        // -------------------------------------------------------------------
        // Sort by along-axis distance (polyline start = 0)
        // -------------------------------------------------------------------
        std::sort(raw_nodes.begin(), raw_nodes.end(),
                  [](const raw_node& a, const raw_node& b)
                  { return a.along_axis < b.along_axis; });

        // -------------------------------------------------------------------
        // Determine which end is the outlet (lower mean elevation)
        // Compare mean elevation of the first and last 10% of nodes.
        // -------------------------------------------------------------------
        const std::size_t n_check =
            std::max(std::size_t{1}, raw_nodes.size() / 10);

        double sum_start = 0.0, sum_end = 0.0;
        for (std::size_t k = 0; k < n_check; ++k)
            sum_start += raw_nodes[k].dem_elev;
        for (std::size_t k = raw_nodes.size() - n_check;
             k < raw_nodes.size(); ++k)
            sum_end += raw_nodes[k].dem_elev;

        // outlet_at_end = true  → polyline end (high along_axis) is downstream
        // outlet_at_end = false → polyline start (along_axis ≈ 0) is downstream
        const bool outlet_at_end = (sum_end < sum_start);
        const double max_along   = raw_nodes.back().along_axis;

        // -------------------------------------------------------------------
        // Build profile nodes with flow_distance (0 = outlet, increasing ↑)
        // -------------------------------------------------------------------
        burned_stream_profile profile;
        profile.feature_id = static_cast<std::int32_t>(fi);
        profile.nodes.reserve(raw_nodes.size());

        for (const auto& rn : raw_nodes)
        {
            burned_stream_node node;
            node.flow_distance =
                outlet_at_end ? (max_along - rn.along_axis) : rn.along_axis;
            node.dem_elevation  = rn.dem_elev;
            node.burned_elevation = 0.0;  // filled below
            node.x = rn.x;
            node.y = rn.y;
            profile.nodes.push_back(node);
        }

        // Sort ascending by flow_distance so index 0 = outlet
        std::sort(profile.nodes.begin(), profile.nodes.end(),
                  [](const burned_stream_node& a, const burned_stream_node& b)
                  { return a.flow_distance < b.flow_distance; });

        // -------------------------------------------------------------------
        // Apply downhill constraint: running minimum from upstream to outlet
        // -------------------------------------------------------------------
        apply_downhill_constraint(profile.nodes);

        std::cout << "  Feature " << fi << ": " << profile.nodes.size()
                  << " channel pixels, outlet elevation = "
                  << profile.nodes.front().dem_elevation
                  << " m, headwater elevation = "
                  << profile.nodes.back().dem_elevation << " m\n";

        profiles.push_back(std::move(profile));
    }

    return profiles;
}

// ---------------------------------------------------------------------------
// write_burned_stream_profiles_csv
// ---------------------------------------------------------------------------
void write_burned_stream_profiles_csv(
    const std::vector<burned_stream_profile>& profiles,
    const std::string&                        filename
)
{
    std::ofstream ofs(filename);
    if (!ofs)
        throw std::runtime_error(
            "write_burned_stream_profiles_csv: cannot open output file: "
            + filename);

    ofs << "feature_id,flow_distance,x,y,dem_elevation,burned_elevation\n";
    ofs.precision(6);
    ofs << std::fixed;

    for (const auto& profile : profiles)
    {
        for (const auto& node : profile.nodes)
        {
            ofs << profile.feature_id << ','
                << node.flow_distance    << ','
                << node.x               << ','
                << node.y               << ','
                << node.dem_elevation   << ','
                << node.burned_elevation << '\n';
        }
    }
}

// ---------------------------------------------------------------------------
// write_burned_stream_profiles_points
// ---------------------------------------------------------------------------
void write_burned_stream_profiles_points(
    const std::vector<burned_stream_profile>& profiles,
    const lsdtt_raster_meta&                  meta,
    const std::string&                        filename,
    point_output_format                       fmt,
    bool                                      spatial_index
)
{
    if (fmt == point_output_format::csv)
        throw std::invalid_argument(
            "write_burned_stream_profiles_points: use "
            "write_burned_stream_profiles_csv for CSV output.");

    lsdtt_vector_schema schema;
    schema.field_names = {
        "feature_id", "flow_distance", "dem_elevation", "burned_elevation"
    };
    schema.field_types = {
        lsdtt_vector_field_type::int64_type,
        lsdtt_vector_field_type::double_type,
        lsdtt_vector_field_type::double_type,
        lsdtt_vector_field_type::double_type
    };

    const vector_output_format vfmt =
        convert_point_to_vector_output_format(fmt);

    lsdtt_data_loader writer;
    writer.begin_vector_write(filename, "burned_stream_profiles",
                              schema, meta.crs_wkt, vfmt, spatial_index);

    for (const auto& profile : profiles)
    {
        for (const auto& node : profile.nodes)
        {
            writer.write_point(
                node.x, node.y,
                {
                    static_cast<std::int64_t>(profile.feature_id),
                    node.flow_distance,
                    node.dem_elevation,
                    node.burned_elevation
                }
            );
        }
    }

    writer.finish_vector_write();
}

} // namespace lsdtt
