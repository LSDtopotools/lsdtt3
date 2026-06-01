// lsdtt_side_tributary_basins_cpu.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html

#include "lsdtt_side_tributary_basins_cpu.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <utility>
#include <vector>
#include <variant>
#include <utility>


#include <ogr_spatialref.h>

#include "../core/lsdtt_vector_utils.hpp"
#include "../utils/lsdtt_raster_geometry.hpp"
#include "lsdtt_raster_cpu.hpp"
#include "lsdtt_raster_binary_mask_vectorise_cpu.hpp"


namespace lsdtt
{
namespace
{

// -----------------------------------------------------------------------------
// String helpers
// -----------------------------------------------------------------------------
std::string precise_string(double value)
{
    std::ostringstream ss;
    ss.precision(17);
    ss << value;
    return ss.str();
}

// -----------------------------------------------------------------------------
// Map projected x/y to raster row/col.
//
// lsdtt_raster_meta stores x_min/y_min as the lower-left raster corner. Rows
// are indexed from the top, so row = floor((y_max - y) / cell_size).
// -----------------------------------------------------------------------------
bool xy_to_row_col(
    const lsdtt_raster_meta& meta,
    double                   x,
    double                   y,
    int&                     row,
    int&                     col
)
{
    if (meta.cell_size <= 0.0)
    {
        return false;
    }

    const double y_max = meta.y_max();

    col = static_cast<int>(std::floor((x - meta.x_min) / meta.cell_size));
    row = static_cast<int>(std::floor((y_max - y) / meta.cell_size));

    return row >= 0 && row < meta.n_rows &&
           col >= 0 && col < meta.n_cols;
}

// -----------------------------------------------------------------------------
// Build a coordinate transformation from vector CRS to raster CRS.
//
// If the vector layer has no CRS, this assumes coordinates are already in the
// raster CRS.
// -----------------------------------------------------------------------------
class local_coordinate_transform
{
public:
    local_coordinate_transform(
        const std::string& source_wkt,
        const std::string& target_wkt
    )
    {
        if (source_wkt.empty() || target_wkt.empty())
        {
            enabled_ = false;
            return;
        }

        source_srs_.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);
        target_srs_.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);

        if (source_srs_.importFromWkt(source_wkt.c_str()) != OGRERR_NONE)
        {
            throw std::runtime_error(
                "side tributary basins: failed to import input point layer CRS"
            );
        }

        if (target_srs_.importFromWkt(target_wkt.c_str()) != OGRERR_NONE)
        {
            throw std::runtime_error(
                "side tributary basins: failed to import raster CRS"
            );
        }

        if (source_srs_.IsSame(&target_srs_) != 0)
        {
            enabled_ = false;
            return;
        }

        transform_ = OGRCreateCoordinateTransformation(
            &source_srs_,
            &target_srs_
        );

        if (!transform_)
        {
            throw std::runtime_error(
                "side tributary basins: failed to create coordinate transformation"
            );
        }

        enabled_ = true;
    }

    ~local_coordinate_transform()
    {
        if (transform_)
        {
            OCTDestroyCoordinateTransformation(transform_);
        }
    }

    local_coordinate_transform(const local_coordinate_transform&) = delete;
    local_coordinate_transform& operator=(const local_coordinate_transform&) = delete;

    bool transform(double& x, double& y) const
    {
        if (!enabled_)
        {
            return true;
        }

        return transform_->Transform(1, &x, &y) != 0;
    }

private:
    bool enabled_ = false;
    OGRSpatialReference source_srs_;
    OGRSpatialReference target_srs_;
    OGRCoordinateTransformation* transform_ = nullptr;
};

// -----------------------------------------------------------------------------
// Load point vector and map points to valid flow flats.
// -----------------------------------------------------------------------------
std::vector<std::int64_t> load_point_vector_as_flow_flats(
    const lsdtt_flow_base& flow,
    const std::string&     vector_path,
    const std::string&     requested_layer_name
)
{
    if (vector_path.empty() || vector_path == "NULL")
    {
        throw std::runtime_error(
            "side tributary basins: input_points_path is empty"
        );
    }

    lsdtt_data_loader loader;
    loader.open_vector(vector_path);

    const std::vector<std::string> layer_names = loader.get_layer_names();

    if (layer_names.empty())
    {
        throw std::runtime_error(
            "side tributary basins: input vector contains no layers: " +
            vector_path
        );
    }

    const std::string layer_name =
        requested_layer_name.empty() ? layer_names.front() : requested_layer_name;

    const std::string vector_crs_wkt =
        loader.get_vector_layer_crs_wkt(layer_name);

    const lsdtt_vector_layer layer =
        loader.read_vector_layer(layer_name);

    if (layer.empty())
    {
        throw std::runtime_error(
            "side tributary basins: input point layer is empty: " + layer_name
        );
    }

    local_coordinate_transform transform(
        vector_crs_wkt,
        flow.meta().crs_wkt
    );

    std::vector<std::int64_t> flats;
    flats.reserve(layer.records.size());

    for (const lsdtt_vector_record& record : layer.records)
    {
        if (!std::holds_alternative<lsdtt_point_geom>(record.geometry))
        {
            continue;
        }

        const lsdtt_point_geom point =
            std::get<lsdtt_point_geom>(record.geometry);

        double x = point.x;
        double y = point.y;

        if (!transform.transform(x, y))
        {
            continue;
        }

        int row = -1;
        int col = -1;

        if (!xy_to_row_col(flow.meta(), x, y, row, col))
        {
            continue;
        }

        const std::int64_t flat = flow.flat_index_unchecked(row, col);

        if (!flow.is_valid_cell(flat))
        {
            continue;
        }

        flats.push_back(flat);
    }

    std::sort(flats.begin(), flats.end());
    flats.erase(std::unique(flats.begin(), flats.end()), flats.end());

    if (flats.empty())
    {
        throw std::runtime_error(
            "side tributary basins: no input points mapped to valid DEM cells"
        );
    }

    return flats;
}

// -----------------------------------------------------------------------------
// Check whether any pixel in an upslope basin is affected by nodata or an edge.
// -----------------------------------------------------------------------------
bool basin_is_complete(
    const lsdtt_flow_base&       flow,
    const lsdtt_flat_index_slice basin_slice
)
{
    for (const std::int64_t* it = basin_slice.begin();
         it != basin_slice.end();
         ++it)
    {
        if (flow.is_influenced_by_invalid(*it))
        {
            return false;
        }
    }

    return true;
}

// -----------------------------------------------------------------------------
// Check if the upslope basin intersects the already traced channel mask.
//
// The outlet itself should not be a traced-channel pixel by construction.
// This helper is defensive: it prevents accepting a side basin that somehow
// contains an input trace in its upslope set.
// -----------------------------------------------------------------------------
bool basin_intersects_traced_channel(
    const lsdtt_flat_index_slice& basin_slice,
    const std::vector<std::uint8_t>& traced_channel_mask
)
{
    for (const std::int64_t* it = basin_slice.begin();
         it != basin_slice.end();
         ++it)
    {
        const std::size_t idx = static_cast<std::size_t>(*it);

        if (idx < traced_channel_mask.size() &&
            traced_channel_mask[idx] != 0)
        {
            return true;
        }
    }

    return false;
}

// -----------------------------------------------------------------------------
// Attribute schema for side outlet and basin outline outputs.
// -----------------------------------------------------------------------------
lsdtt_vector_schema side_basin_schema()
{
    lsdtt_vector_schema schema;

    schema.field_names =
    {
        "basin_id",
        "outlet_flat",
        "outlet_row",
        "outlet_col",
        "receiver_channel_flat",
        "receiver_channel_row",
        "receiver_channel_col",
        "n_pixels",
        "drainage_area",
        "complete"
    };

    schema.field_types =
    {
        lsdtt_vector_field_type::int64_type,
        lsdtt_vector_field_type::int64_type,
        lsdtt_vector_field_type::int64_type,
        lsdtt_vector_field_type::int64_type,
        lsdtt_vector_field_type::int64_type,
        lsdtt_vector_field_type::int64_type,
        lsdtt_vector_field_type::int64_type,
        lsdtt_vector_field_type::int64_type,
        lsdtt_vector_field_type::double_type,
        lsdtt_vector_field_type::int64_type
    };

    return schema;
}

std::vector<lsdtt_vector_field_value> side_basin_field_values(
    const lsdtt_side_tributary_basin_record& basin,
    std::int64_t                             basin_id
)
{
    return
    {
        basin_id,
        basin.outlet_flat,
        static_cast<std::int64_t>(basin.outlet_row),
        static_cast<std::int64_t>(basin.outlet_col),
        basin.receiver_channel_flat,
        static_cast<std::int64_t>(basin.receiver_channel_row),
        static_cast<std::int64_t>(basin.receiver_channel_col),
        static_cast<std::int64_t>(basin.n_pixels),
        basin.drainage_area,
        static_cast<std::int64_t>(basin.complete ? 1 : 0)
    };
}

} // unnamed namespace

// =============================================================================
// extract_side_tributary_basins
// =============================================================================
lsdtt_side_tributary_basin_result extract_side_tributary_basins(
    const lsdtt_flow_base& flow,
    const lsdtt_side_tributary_basin_options& opts
)
{
    if (flow.empty())
    {
        throw std::runtime_error(
            "extract_side_tributary_basins: flow object is empty"
        );
    }

    if (!flow.has_s_vector())
    {
        throw std::runtime_error(
            "extract_side_tributary_basins: flow object must have s_vector. "
            "Set lsdtt_flowinfo_options::build_s_vector = true."
        );
    }

    if (opts.minimum_contributing_pixels <= 0)
    {
        throw std::invalid_argument(
            "extract_side_tributary_basins: minimum_contributing_pixels must be positive"
        );
    }

    const std::size_t n = flow.full_grid_size();

    lsdtt_side_tributary_basin_result result;
    result.traced_channel_mask.assign(n, std::uint8_t{0});

    const std::vector<std::int64_t> seed_flats =
        load_point_vector_as_flow_flats(
            flow,
            opts.input_points_path,
            opts.input_points_layer
        );

    // -------------------------------------------------------------------------
    // Trace D8 downstream from each seed point. If a trace reaches an already
    // traced pixel, stop; all downstream overlap has already been recorded.
    // -------------------------------------------------------------------------
    for (const std::int64_t seed : seed_flats)
    {
        if (!flow.is_valid_cell(seed))
        {
            continue;
        }

        std::int64_t current = seed;
        std::int64_t steps   = 0;

        const std::int64_t max_steps =
            static_cast<std::int64_t>(flow.full_grid_size());

        while (flow.is_valid_cell(current))
        {
            const std::size_t current_idx =
                static_cast<std::size_t>(current);

            if (result.traced_channel_mask[current_idx] != 0)
            {
                break;
            }

            result.traced_channel_mask[current_idx] = 1u;

            if (flow.is_outlet_unchecked(current))
            {
                break;
            }

            const std::int64_t receiver = flow.receiver_of_unchecked(current);

            if (receiver < 0)
            {
                break;
            }

            current = receiver;
            ++steps;

            if (steps > max_steps)
            {
                throw std::runtime_error(
                    "extract_side_tributary_basins: downstream trace exceeded raster size"
                );
            }
        }
    }

    // -------------------------------------------------------------------------
    // Find side tributary outlet pixels:
    // non-channel valid cells whose D8 receiver lies on the traced channel.
    // -------------------------------------------------------------------------
    std::vector<std::int64_t> candidate_outlets;
    candidate_outlets.reserve(1024);

    for (std::int64_t flat = 0;
         flat < static_cast<std::int64_t>(n);
         ++flat)
    {
        if (!flow.is_valid_cell(flat))
        {
            continue;
        }

        if (result.traced_channel_mask[static_cast<std::size_t>(flat)] != 0)
        {
            continue;
        }

        const std::int64_t receiver = flow.receiver_of_unchecked(flat);

        if (receiver < 0)
        {
            continue;
        }

        if (result.traced_channel_mask[static_cast<std::size_t>(receiver)] == 0)
        {
            continue;
        }

        if (flow.s_vector_size(flat) <
            static_cast<std::uint64_t>(opts.minimum_contributing_pixels))
        {
            continue;
        }

        candidate_outlets.push_back(flat);
    }

    std::sort(candidate_outlets.begin(), candidate_outlets.end());
    candidate_outlets.erase(
        std::unique(candidate_outlets.begin(), candidate_outlets.end()),
        candidate_outlets.end()
    );

    // -------------------------------------------------------------------------
    // Completion test and final record construction.
    // -------------------------------------------------------------------------
    result.accepted_basins.reserve(candidate_outlets.size());

    for (const std::int64_t outlet : candidate_outlets)
    {
        const lsdtt_flat_index_slice basin_slice =
            flow.upslope_nodes(outlet);

        if (basin_intersects_traced_channel(
                basin_slice,
                result.traced_channel_mask))
        {
            continue;
        }

        const bool complete = basin_is_complete(flow, basin_slice);

        if (opts.reject_incomplete_basins && !complete)
        {
            continue;
        }

        const std::int64_t receiver =
            flow.receiver_of_unchecked(outlet);

        if (receiver < 0)
        {
            continue;
        }

        lsdtt_side_tributary_basin_record record;

        record.outlet_flat = outlet;
        record.outlet_row  = flow.row_from_flat_unchecked(outlet);
        record.outlet_col  = flow.col_from_flat_unchecked(outlet);

        record.receiver_channel_flat = receiver;
        record.receiver_channel_row  = flow.row_from_flat_unchecked(receiver);
        record.receiver_channel_col  = flow.col_from_flat_unchecked(receiver);

        record.n_pixels = flow.s_vector_size(outlet);

        record.drainage_area =
            static_cast<double>(record.n_pixels) *
            flow.cell_size() *
            flow.cell_size();

        record.complete = complete;

        result.accepted_basins.push_back(record);
    }

    return result;
}

// =============================================================================
// write_side_tributary_outlets
// =============================================================================
void write_side_tributary_outlets(
    const lsdtt_flow_base& flow,
    const std::vector<lsdtt_side_tributary_basin_record>& basins,
    const std::string& out_path,
    point_output_format format,
    const std::string& layer_name,
    bool spatial_index
)
{
    if (basins.empty())
    {
        throw std::runtime_error(
            "write_side_tributary_outlets: no accepted side tributary basins to write"
        );
    }

    lsdtt_data_loader loader;

    if (format == point_output_format::csv)
    {
        if (flow.meta().crs_wkt.empty())
        {
            throw std::runtime_error(
                "write_side_tributary_outlets: CSV output requires raster CRS"
            );
        }

        const lsdtt_coordinate_transform transform =
            lsdtt_coordinate_transform::from_epsg(flow.meta().crs_wkt, 4326);

        std::vector<std::string> headers =
        {
            "latitude",
            "longitude",
            "basin_id",
            "outlet_flat",
            "outlet_row",
            "outlet_col",
            "receiver_channel_flat",
            "receiver_channel_row",
            "receiver_channel_col",
            "n_pixels",
            "drainage_area",
            "complete"
        };

        std::vector<std::vector<std::string>> rows;
        rows.reserve(basins.size());

        for (std::size_t i = 0; i < basins.size(); ++i)
        {
            const auto& basin = basins[i];

            const lsdtt_lonlat ll =
                cell_center_to_lonlat(
                    flow.meta(),
                    basin.outlet_row,
                    basin.outlet_col,
                    transform
                );

            rows.push_back(
                {
                    precise_string(ll.latitude),
                    precise_string(ll.longitude),
                    std::to_string(i),
                    std::to_string(basin.outlet_flat),
                    std::to_string(basin.outlet_row),
                    std::to_string(basin.outlet_col),
                    std::to_string(basin.receiver_channel_flat),
                    std::to_string(basin.receiver_channel_row),
                    std::to_string(basin.receiver_channel_col),
                    std::to_string(basin.n_pixels),
                    precise_string(basin.drainage_area),
                    std::to_string(basin.complete ? 1 : 0)
                }
            );
        }

        loader.write_csv(out_path, headers, rows);
        return;
    }

    const vector_output_format vector_format =
        convert_point_to_vector_output_format(format);

    loader.begin_vector_write(
        out_path,
        layer_name,
        side_basin_schema(),
        flow.meta().crs_wkt,
        vector_format,
        spatial_index
    );

    for (std::size_t i = 0; i < basins.size(); ++i)
    {
        const auto& basin = basins[i];

        const lsdtt_xy xy =
            cell_center_from_row_col(
                flow.meta(),
                basin.outlet_row,
                basin.outlet_col
            );

        loader.write_point(
            xy.x,
            xy.y,
            side_basin_field_values(
                basin,
                static_cast<std::int64_t>(i)
            )
        );
    }

    loader.finish_vector_write();
}

// =============================================================================
// write_side_tributary_basin_outlines
//
// This writer builds a labelled basin raster. For each basin, all pixels in its
// upslope set receive basin_id + 1. Then each basin is vectorised independently
// as a binary mask so that attributes can be attached cleanly.
// =============================================================================
void write_side_tributary_basin_outlines(
    const lsdtt_flow_base& flow,
    const std::vector<lsdtt_side_tributary_basin_record>& basins,
    const std::string& out_path,
    vector_output_format format,
    const std::string& layer_name,
    bool spatial_index
)
{
    if (basins.empty())
    {
        throw std::runtime_error(
            "write_side_tributary_basin_outlines: no accepted side tributary basins to write"
        );
    }

    if (!flow.has_s_vector())
    {
        throw std::runtime_error(
            "write_side_tributary_basin_outlines: flow object must have s_vector"
        );
    }

    const std::size_t n = flow.full_grid_size();

    // -------------------------------------------------------------------------
    // Build a labelled raster:
    //   label = basin_id + 1
    //   0     = background / nodata
    //
    // Side basins should be non-overlapping by construction. If an overlap is
    // encountered, the first basin written wins.
    // -------------------------------------------------------------------------
    lsdtt_raster_meta label_meta = flow.meta();
    label_meta.data_type = GDT_UInt32;
    label_meta.no_data_value = 0.0;

    std::vector<float> label_data(n, 0.0f);
    std::vector<std::uint8_t> label_mask(n, std::uint8_t{0});

    for (std::size_t basin_id = 0; basin_id < basins.size(); ++basin_id)
    {
        const auto& basin = basins[basin_id];

        const float label_value =
            static_cast<float>(static_cast<std::int64_t>(basin_id) + 1);

        const lsdtt_flat_index_slice slice =
            flow.upslope_nodes(basin.outlet_flat);

        for (const std::int64_t* it = slice.begin();
             it != slice.end();
             ++it)
        {
            const std::size_t idx = static_cast<std::size_t>(*it);

            if (idx >= n)
            {
                continue;
            }

            // First basin wins if there is an unexpected overlap.
            if (label_mask[idx] == 0)
            {
                label_data[idx] = label_value;
                label_mask[idx] = std::uint8_t{1};
            }
        }
    }

    const lsdtt_raster_cpu label_raster(
        std::move(label_data),
        std::move(label_mask),
        label_meta
    );

    // -------------------------------------------------------------------------
    // Vectorise labelled components to polygon records.
    //
    // The vectoriser expects the schema:
    //   label_id int64
    //   area_m2  double
    //
    // We then rewrite each record with the richer side-basin schema.
    // -------------------------------------------------------------------------
    lsdtt_vector_schema component_schema;
    component_schema.field_names = {
        "label_id",
        "area_m2"
    };
    component_schema.field_types = {
        lsdtt_vector_field_type::int64_type,
        lsdtt_vector_field_type::double_type
    };

    const std::vector<lsdtt_vector_record> component_records =
        binary_components_to_polygon_records_with_holes(
            label_raster,
            nullptr,
            component_schema
        );

    lsdtt_data_loader loader;

    loader.begin_vector_write(
        out_path,
        layer_name,
        side_basin_schema(),
        flow.meta().crs_wkt,
        format,
        spatial_index
    );

    for (const lsdtt_vector_record& component_record : component_records)
    {
        if (!std::holds_alternative<lsdtt_polygon_geom>(
                component_record.geometry))
        {
            continue;
        }

        if (component_record.field_values.empty())
        {
            continue;
        }

        std::int64_t label_id = -1;

        if (std::holds_alternative<std::int64_t>(
                component_record.field_values[0]))
        {
            label_id = std::get<std::int64_t>(
                component_record.field_values[0]
            );
        }
        else if (std::holds_alternative<double>(
                     component_record.field_values[0]))
        {
            label_id = static_cast<std::int64_t>(
                std::llround(std::get<double>(
                    component_record.field_values[0]
                ))
            );
        }
        else
        {
            continue;
        }

        if (label_id <= 0)
        {
            continue;
        }

        const std::size_t basin_index =
            static_cast<std::size_t>(label_id - 1);

        if (basin_index >= basins.size())
        {
            continue;
        }

        const lsdtt_polygon_geom polygon =
            std::get<lsdtt_polygon_geom>(component_record.geometry);

        loader.write_polygon(
            polygon.exterior,
            polygon.holes,
            side_basin_field_values(
                basins[basin_index],
                static_cast<std::int64_t>(basin_index)
            )
        );
    }

    loader.finish_vector_write();
}


} // namespace lsdtt
