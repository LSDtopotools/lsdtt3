#include <algorithm>
#include <cmath>
#include <cstdint>
#include <memory>
#include <limits>
#include <map>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <utility>
#include <vector>

#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>

#include <cpu/lsdtt_basin_cpu.hpp>
#include <cpu/lsdtt_junction_basin_manager_cpu.hpp>
#include <cpu/lsdtt_raster_surface_metrics_cpu.hpp>
#include <cpu/lsdtt_raster_spectral_cpu.hpp>
#include <cpu/lsdtt_side_tributary_basins_cpu.hpp>
#include <utils/lsdtt_statistics.hpp>
#include <utils/lsdtt_polygon_geometry.hpp>

#include <core/lsdtt_raster_warp.hpp>

#include "binding_common.hpp"

namespace lsdtt3_python {
namespace {

std::vector<std::int64_t> outlet_flats(const FlowInfo& flow, nb::object outlets,
                                       std::int64_t min_upslope_pixels)
{
    if (min_upslope_pixels < 0)
        throw std::invalid_argument("min_upslope_pixels must be non-negative");
    if (outlets.is_none())
        return lsdtt::get_all_basin_outlet_flats(flow.native(), min_upslope_pixels);

    nb::object array = nb::module_::import_("numpy").attr("asarray")(outlets, nb::arg("dtype") = "int64");
    if (nb::cast<int>(array.attr("ndim")) != 2 || nb::cast<int>(array.attr("shape")[1]) != 2)
        throw std::invalid_argument("outlets must be an N by 2 row, col array");
    const auto pairs = nb::cast<std::vector<std::vector<std::int64_t>>>(array.attr("tolist")());
    std::unordered_set<std::int64_t> seen;
    std::vector<std::int64_t> result;
    for (const auto& pair : pairs) {
        const auto row = pair[0]; const auto col = pair[1];
        if (row < 0 || row >= flow.native().n_rows() || col < 0 || col >= flow.native().n_cols())
            throw std::invalid_argument("outlets contains an out-of-grid pixel");
        const auto flat = row * static_cast<std::int64_t>(flow.native().n_cols()) + col;
        if (!flow.native().is_valid_cell(flat)) throw std::invalid_argument("outlets selects a nodata pixel");
        if (seen.insert(flat).second) result.push_back(flat);
    }
    return result;
}

class BasinCollection {
public:
    BasinCollection(const BasinCollection&) = delete;
    BasinCollection& operator=(const BasinCollection&) = delete;
    BasinCollection(BasinCollection&&) = default;
    BasinCollection& operator=(BasinCollection&&) = default;
    BasinCollection(std::shared_ptr<lsdtt::lsdtt_flowinfo> flow, std::vector<std::int64_t> outlets)
        : flow_(std::move(flow))
    {
        basins_.reserve(outlets.size());
        for (const auto outlet : outlets)
            basins_.push_back(std::make_unique<lsdtt::lsdtt_basin_cpu>(outlet, *flow_));
        parent_.assign(basins_.size(), -1);
        for (std::size_t child = 0; child < basins_.size(); ++child) {
            const auto child_outlet = basins_[child]->get_basin_outlet_flat_index();
            std::int32_t best_parent = -1;
            std::int32_t best_size = std::numeric_limits<std::int32_t>::max();
            for (std::size_t candidate = 0; candidate < basins_.size(); ++candidate) {
                if (candidate == child) continue;
                const auto parent_outlet = basins_[candidate]->get_basin_outlet_flat_index();
                if (!flow_->is_upslope_of(child_outlet, parent_outlet)) continue;
                const auto n = basins_[candidate]->pixel_count();
                if (n < best_size) { best_size = n; best_parent = static_cast<std::int32_t>(candidate); }
            }
            parent_[child] = best_parent;
        }
        depth_.assign(basins_.size(), 0);
        for (std::size_t i = 0; i < basins_.size(); ++i) {
            for (auto parent = parent_[i]; parent >= 0; parent = parent_[static_cast<std::size_t>(parent)])
                ++depth_[i];
        }
    }

    std::size_t size() const noexcept { return basins_.size(); }
    const lsdtt::lsdtt_flowinfo& flow() const noexcept { return *flow_; }
    const lsdtt::lsdtt_basin_cpu& basin(std::size_t id) const {
        if (id >= basins_.size()) throw std::out_of_range("basin_id is out of range");
        return *basins_[id];
    }
    std::int32_t parent(std::size_t id) const { if (id >= parent_.size()) throw std::out_of_range("basin_id is out of range"); return parent_[id]; }
    std::int32_t depth(std::size_t id) const { if (id >= depth_.size()) throw std::out_of_range("basin_id is out of range"); return depth_[id]; }
    std::vector<std::int32_t> memberships(std::int64_t flat) const {
        std::vector<std::int32_t> ids;
        for (std::size_t i = 0; i < basins_.size(); ++i)
            if (flow_->is_upslope_of(flat, basins_[i]->get_basin_outlet_flat_index())) ids.push_back(static_cast<std::int32_t>(i));
        return ids;
    }
private:
    std::shared_ptr<lsdtt::lsdtt_flowinfo> flow_;
    std::vector<std::unique_ptr<lsdtt::lsdtt_basin_cpu>> basins_;
    std::vector<std::int32_t> parent_, depth_;
};

Raster basin_mask(const BasinCollection& collection, std::size_t id)
{
    const auto& flow = collection.flow();
    const std::size_t n = flow.full_grid_size();
    std::vector<float> data(n, 0.0f);
    std::vector<std::uint8_t> mask(flow.mask_array().data(), flow.mask_array().data() + n);
    for (const auto flat : collection.basin(id).get_slice()) data[static_cast<std::size_t>(flat)] = 1.0f;
    return Raster(lsdtt::lsdtt_raster_cpu(std::move(data), std::move(mask), flow.meta()));
}

Raster basin_labels(const BasinCollection& collection)
{
    const auto& flow = collection.flow(); const std::size_t n = flow.full_grid_size();
    std::vector<float> data(n, 0.0f);
    std::vector<std::uint8_t> mask(flow.mask_array().data(), flow.mask_array().data() + n);
    for (std::size_t id = 0; id < collection.size(); ++id)
        for (const auto flat : collection.basin(id).get_slice()) data[static_cast<std::size_t>(flat)] = static_cast<float>(id + 1);
    return Raster(lsdtt::lsdtt_raster_cpu(std::move(data), std::move(mask), flow.meta()));
}

nb::object basin_pixels(const BasinCollection& collection, std::size_t id)
{
    const auto& flow = collection.flow(); const auto& slice = collection.basin(id).get_slice();
    std::vector<std::int64_t> values; values.reserve(slice.size() * 2);
    for (const auto flat : slice) { values.push_back(flow.row_from_flat_unchecked(flat)); values.push_back(flow.col_from_flat_unchecked(flat)); }
    return nb::module_::import_("numpy").attr("array")(values, nb::arg("dtype") = "int64")
        .attr("reshape")(nb::make_tuple(static_cast<int>(slice.size()), 2));
}

nb::object basin_table(const BasinCollection& collection)
{
    const auto& flow = collection.flow(); nb::list rows;
    for (std::size_t id = 0; id < collection.size(); ++id) {
        const auto& basin = collection.basin(id); const auto outlet = basin.get_basin_outlet_flat_index();
        nb::dict row; row["basin_id"] = id; row["outlet_flat"] = outlet;
        row["outlet_row"] = flow.row_from_flat_unchecked(outlet); row["outlet_col"] = flow.col_from_flat_unchecked(outlet);
        row["n_pixels"] = basin.pixel_count(); row["area_metres2"] = basin.pixel_count() * flow.meta().cell_size * flow.meta().cell_size;
        row["complete"] = basin.check_basin_completeness(); row["outlet_on_dem_edge"] = lsdtt::is_flat_on_dem_edge(flow, outlet);
        row["parent_id"] = collection.parent(id); row["nesting_depth"] = collection.depth(id); rows.append(std::move(row));
    }
    return nb::module_::import_("pandas").attr("DataFrame")(rows);
}

nb::object basin_membership(const BasinCollection& collection, int row, int col)
{
    const auto& flow = collection.flow(); const auto flat = flow.flat_index(row, col);
    if (!flow.is_valid_cell(flat)) throw std::invalid_argument("row, col selects nodata");
    nb::object numpy = nb::module_::import_("numpy");
    return numpy.attr("array")(collection.memberships(flat), nb::arg("dtype") = "int32");
}

BasinCollection make_basins(const FlowInfo& flow, nb::object outlets, std::int64_t min_upslope_pixels)
{
    return BasinCollection(flow.shared_native(), outlet_flats(flow, outlets, min_upslope_pixels));
}

class JunctionBasinCollection {
public:
    JunctionBasinCollection(std::shared_ptr<lsdtt::lsdtt_flowinfo> flow,
                            std::shared_ptr<lsdtt::lsdtt_junction_network_cpu> network,
                            const std::string& selection, std::vector<std::int64_t> outlets)
        : flow_(std::move(flow)), network_(std::move(network)),
          manager_(std::make_unique<lsdtt::lsdtt_junction_basin_manager_cpu>(*flow_, *network_))
    {
        if (selection == "all") manager_->select_all_junction_basins();
        else if (selection == "outlets") manager_->select_junction_basins_from_basin_outlets(outlets);
        else if (selection == "nested") manager_->select_junction_basins_from_basin_outlets_nested(outlets);
        else throw std::invalid_argument("selection must be 'all', 'outlets', or 'nested'");
    }
    std::size_t size() const noexcept { return manager_->size(); }
    const lsdtt::lsdtt_junction_basin_manager_cpu& native() const noexcept { return *manager_; }
    const lsdtt::lsdtt_flowinfo& flow() const noexcept { return *flow_; }
    JunctionBasinCollection& filter_incomplete() { manager_->filter_incomplete(); return *this; }
    JunctionBasinCollection& filter_area(std::int32_t low, std::int32_t high) { manager_->filter_by_area(low, high); return *this; }
    JunctionBasinCollection& filter_stream_order(std::int32_t low, std::int32_t high) { manager_->filter_by_stream_order(low, high); return *this; }
    JunctionBasinCollection& filter_nested() { manager_->filter_nested(); return *this; }
    JunctionBasinCollection& filter_nested_by_stream_order() { manager_->filter_nested_by_stream_order(); return *this; }
    JunctionBasinCollection& filter_largest() { manager_->filter_largest(); return *this; }
private:
    std::shared_ptr<lsdtt::lsdtt_flowinfo> flow_;
    std::shared_ptr<lsdtt::lsdtt_junction_network_cpu> network_;
    std::unique_ptr<lsdtt::lsdtt_junction_basin_manager_cpu> manager_;
};

JunctionBasinCollection make_junction_basins(const ChannelNetwork& network, const std::string& selection,
                                             nb::object outlets)
{
    std::vector<std::int64_t> flats;
    if (selection != "all") {
        nb::object array = nb::module_::import_("numpy").attr("asarray")(outlets, nb::arg("dtype") = "int64");
        if (nb::cast<int>(array.attr("ndim")) != 2 || nb::cast<int>(array.attr("shape")[1]) != 2)
            throw std::invalid_argument("outlets must be an N by 2 row, col array");
        const auto values = nb::cast<std::vector<std::int64_t>>(array.attr("ravel")());
        const auto& flow = network.flow();
        for (std::size_t i = 0; i < values.size(); i += 2) {
            const int row = static_cast<int>(values[i]), col = static_cast<int>(values[i + 1]);
            if (!flow.in_bounds(row, col) || !flow.is_valid_cell(flow.flat_index(row, col)))
                throw std::invalid_argument("outlets contains an invalid flow pixel");
            flats.push_back(flow.flat_index(row, col));
        }
    }
    return JunctionBasinCollection(network.shared_flow(), network.shared_native(), selection, std::move(flats));
}

Raster junction_basin_labels(const JunctionBasinCollection& collection)
{
    auto result = collection.native().get_junction_basins_raster();
    return Raster(lsdtt::lsdtt_raster_cpu(std::move(result.data), std::move(result.mask), result.meta));
}

Raster junction_basin_mask(const JunctionBasinCollection& collection, std::size_t id)
{
    if (id >= collection.size()) throw std::out_of_range("basin_id is out of range");
    const auto& flow = collection.flow(); std::vector<float> data(flow.full_grid_size(), 0.0f);
    std::vector<std::uint8_t> mask(flow.mask_array().data(), flow.mask_array().data() + flow.full_grid_size());
    for (const auto flat : collection.native().basin_at(id).get_slice()) data[static_cast<std::size_t>(flat)] = 1.0f;
    return Raster(lsdtt::lsdtt_raster_cpu(std::move(data), std::move(mask), flow.meta()));
}

nb::object junction_basin_table(const JunctionBasinCollection& collection)
{
    const auto& flow = collection.flow(); nb::list rows;
    for (std::size_t id = 0; id < collection.size(); ++id) {
        const auto& basin = collection.native().basin_at(id); const auto flat = basin.get_basin_outlet_flat_index();
        nb::dict row; row["basin_id"] = id; row["junction_id"] = basin.get_junction_number(); row["stream_order"] = basin.get_stream_order();
        row["outlet_flat"] = flat; row["outlet_row"] = flow.row_from_flat_unchecked(flat); row["outlet_col"] = flow.col_from_flat_unchecked(flat);
        row["n_pixels"] = basin.pixel_count(); row["area_metres2"] = basin.pixel_count() * flow.cell_size() * flow.cell_size(); row["complete"] = basin.check_basin_completeness(); rows.append(std::move(row));
    }
    return nb::module_::import_("pandas").attr("DataFrame")(rows);
}

nb::object basin_statistics(const BasinCollection& collection, const Raster& dem, bool include_slope, bool include_flow_distance,
                            bool include_shape_metrics, double perimeter_simplify_cells)
{
    const auto& flow = collection.flow();
    if (dem.native().get_meta().n_rows != flow.n_rows() || dem.native().get_meta().n_cols != flow.n_cols())
        throw std::invalid_argument("dem must have the same grid shape as the basin flow topology");
    const auto slope_raster = include_slope ? lsdtt::compute_gradient_magnitude_3x3_square(dem.native()) : lsdtt::lsdtt_raster_cpu();
    const auto distance_raster = (include_flow_distance || include_shape_metrics) ? lsdtt::lsdtt_raster_cpu(flow.get_distance_from_outlet_raster()) : lsdtt::lsdtt_raster_cpu();
    nb::list rows;
    for (std::size_t id = 0; id < collection.size(); ++id) {
        const auto& basin = collection.basin(id); const auto z = basin.get_basin_values(dem.native()); const auto zs = lsdtt::compute_descriptive_stats(z);
        nb::dict row; row["basin_id"] = id; row["n_pixels"] = basin.pixel_count(); row["area_metres2"] = basin.pixel_count() * flow.cell_size() * flow.cell_size();
        row["elevation_min"] = zs.min; row["elevation_max"] = zs.max; row["elevation_mean"] = zs.mean; row["elevation_median"] = zs.median; row["elevation_stddev"] = zs.stddev_population;
        if (include_slope) { const auto s = basin.get_basin_values(slope_raster); row["slope_25th"] = lsdtt::compute_percentile(s, 25.0); row["slope_50th"] = lsdtt::compute_percentile(s, 50.0); row["slope_75th"] = lsdtt::compute_percentile(s, 75.0); }
        if (include_flow_distance) { const auto d = basin.get_basin_values(distance_raster); row["flow_distance_25th"] = lsdtt::compute_percentile(d, 25.0); row["flow_distance_50th"] = lsdtt::compute_percentile(d, 50.0); row["flow_distance_75th"] = lsdtt::compute_percentile(d, 75.0); }
        if (include_shape_metrics) {
            const auto outline = basin.get_basin_outline(); const auto& meta = outline.divide_grid_meta(); std::vector<lsdtt::lsdtt_xy> ring; ring.reserve(outline.outline_divide_rc().size() + 1);
            for (const auto& rc : outline.outline_divide_rc()) ring.push_back({meta.x_min + rc.col * meta.cell_size, meta.y_max() - rc.row * meta.cell_size});
            if (!ring.empty() && (ring.front().x != ring.back().x || ring.front().y != ring.back().y)) ring.push_back(ring.front());
            const auto d = basin.get_basin_values(distance_raster); const auto shape = lsdtt::compute_basin_shape(ring, d.empty() ? 0.0 : *std::max_element(d.begin(), d.end()), zs.min, zs.mean, zs.max, perimeter_simplify_cells * flow.cell_size());
            double sum_x = 0.0, sum_y = 0.0; for (const auto flat : basin.get_slice()) { const auto rr = flow.row_from_flat_unchecked(flat); const auto cc = flow.col_from_flat_unchecked(flat); sum_x += flow.meta().x_min + (cc + 0.5) * flow.cell_size(); sum_y += flow.meta().y_max() - (rr + 0.5) * flow.cell_size(); }
            const double n = static_cast<double>(basin.pixel_count()); row["centroid_x"] = sum_x / n; row["centroid_y"] = sum_y / n; row["basin_length_metres"] = shape.basin_length; row["elongation_ratio"] = shape.elongation_ratio; row["form_factor"] = shape.form_factor; row["circularity_ratio"] = shape.circularity_ratio; row["gravelius_compactness"] = shape.gravelius_compactness; row["solidity"] = shape.solidity; row["eccentricity"] = shape.eccentricity; row["orientation_degrees"] = shape.orientation_deg; row["hypsometric_integral"] = shape.hypsometric_integral;
        }
        rows.append(std::move(row));
    }
    return nb::module_::import_("pandas").attr("DataFrame")(rows);
}

nb::object basin_polygons(const BasinCollection& collection)
{
    const auto& flow = collection.flow(); nb::list geometry; std::vector<std::size_t> ids; std::vector<std::int64_t> outlets; nb::object polygon = nb::module_::import_("shapely.geometry").attr("Polygon");
    for (std::size_t id = 0; id < collection.size(); ++id) { const auto outline = collection.basin(id).get_basin_outline(); const auto& meta = outline.divide_grid_meta(); nb::list coordinates; for (const auto& rc : outline.outline_divide_rc()) coordinates.append(nb::make_tuple(meta.x_min + rc.col * meta.cell_size, meta.y_max() - rc.row * meta.cell_size)); ids.push_back(id); outlets.push_back(collection.basin(id).get_basin_outlet_flat_index()); geometry.append(polygon(coordinates)); }
    nb::dict data; data["basin_id"] = ids; data["outlet_flat"] = outlets; nb::object ctor = nb::module_::import_("geopandas").attr("GeoDataFrame"); return flow.meta().crs_wkt.empty() ? ctor(data, nb::arg("geometry") = geometry) : ctor(data, nb::arg("geometry") = geometry, nb::arg("crs") = flow.meta().crs_wkt);
}

nb::object basin_outline_ridgelines(const BasinCollection& collection, const Raster& dem,
                                    const std::string& split_method)
{
    const auto& flow = collection.flow();
    const auto& meta = dem.native().get_meta();
    if (meta.n_rows != flow.n_rows() || meta.n_cols != flow.n_cols() ||
        meta.x_min != flow.meta().x_min || meta.y_min != flow.meta().y_min ||
        meta.cell_size != flow.meta().cell_size ||
        !lsdtt::crs_equal_semantically(meta.crs_wkt, flow.meta().crs_wkt))
        throw std::invalid_argument("dem must share the BasinCollection flow grid");
    if (split_method != "longest_channel" && split_method != "maximum_elevation" &&
        split_method != "equal_length")
        throw std::invalid_argument(
            "split_method must be 'longest_channel', 'maximum_elevation', or 'equal_length'");

    auto distance_result = flow.get_distance_from_outlet_raster();
    const lsdtt::lsdtt_raster_cpu distance(
        std::move(distance_result.data), std::move(distance_result.mask),
        std::move(distance_result.meta));
    nb::object line_string = nb::module_::import_("shapely.geometry").attr("LineString");
    nb::list geometry;
    std::vector<std::int32_t> basin_id, ridgeline_id, n_nodes;
    std::vector<std::int64_t> outlet_flat;
    std::vector<double> length_metres;
    std::vector<float> elevation_min, elevation_max;

    for (std::size_t basin = 0; basin < collection.size(); ++basin) {
        const auto outline = collection.basin(basin).get_basin_outline();
        std::vector<lsdtt::lsdtt_ridgeline_cpu> ridges;
        if (split_method == "longest_channel")
            ridges = outline.split_outline_into_ridgelines_by_longest_channel(distance);
        else if (split_method == "maximum_elevation")
            ridges = outline.split_outline_into_ridgelines_by_max_elevation(dem.native());
        else
            ridges = outline.split_outline_into_ridgelines_by_equal_length();

        for (std::size_t ridge = 0; ridge < ridges.size(); ++ridge) {
            auto& native = ridges[ridge];
            native.fill_elevation(dem.native());
            nb::list coordinates;
            for (const auto& xy : native.positions())
                coordinates.append(nb::make_tuple(xy.x, xy.y));
            if (nb::len(coordinates) < 2) continue;
            const auto& elevation = native.elevations();
            basin_id.push_back(static_cast<std::int32_t>(basin));
            ridgeline_id.push_back(static_cast<std::int32_t>(ridge));
            outlet_flat.push_back(collection.basin(basin).get_basin_outlet_flat_index());
            n_nodes.push_back(static_cast<std::int32_t>(native.size()));
            length_metres.push_back(native.total_length());
            const auto finite = std::minmax_element(elevation.begin(), elevation.end());
            elevation_min.push_back(elevation.empty() ? std::numeric_limits<float>::quiet_NaN() : *finite.first);
            elevation_max.push_back(elevation.empty() ? std::numeric_limits<float>::quiet_NaN() : *finite.second);
            geometry.append(line_string(coordinates));
        }
    }
    nb::dict data;
    data["basin_id"] = basin_id; data["ridgeline_id"] = ridgeline_id;
    data["outlet_flat"] = outlet_flat; data["n_nodes"] = n_nodes;
    data["length_metres"] = length_metres; data["elevation_min"] = elevation_min;
    data["elevation_max"] = elevation_max;
    nb::object ctor = nb::module_::import_("geopandas").attr("GeoDataFrame");
    return flow.meta().crs_wkt.empty() ? ctor(data, nb::arg("geometry") = geometry)
                                      : ctor(data, nb::arg("geometry") = geometry, nb::arg("crs") = flow.meta().crs_wkt);
}

nb::object basin_ridgeline_extrema(const BasinCollection& collection,
                                   const ChannelNetwork& network, const Raster& dem,
                                   const std::string& split_method,
                                   float snout_smoothing_fraction,
                                   int detrend_smoothing_pixels)
{
    if (&network.flow() != &collection.flow())
        throw std::invalid_argument("network and basins must belong to the same FlowInfo");
    if (snout_smoothing_fraction <= 0.0f || snout_smoothing_fraction > 1.0f ||
        detrend_smoothing_pixels < 1)
        throw std::invalid_argument("invalid ridgeline smoothing controls");
    const auto& flow = collection.flow();
    const auto& meta = dem.native().get_meta();
    if (meta.n_rows != flow.n_rows() || meta.n_cols != flow.n_cols() ||
        meta.x_min != flow.meta().x_min || meta.y_min != flow.meta().y_min ||
        meta.cell_size != flow.meta().cell_size ||
        !lsdtt::crs_equal_semantically(meta.crs_wkt, flow.meta().crs_wkt))
        throw std::invalid_argument("dem must share the BasinCollection flow grid");
    if (split_method != "longest_channel" && split_method != "maximum_elevation" &&
        split_method != "equal_length")
        throw std::invalid_argument("split_method must be 'longest_channel', 'maximum_elevation', or 'equal_length'");

  auto distance_result = flow.get_distance_from_outlet_raster();
  const lsdtt::lsdtt_raster_cpu distance(
      std::move(distance_result.data), std::move(distance_result.mask),
      std::move(distance_result.meta));
    nb::list rows, geometry;
    nb::object point = nb::module_::import_("shapely.geometry").attr("Point");
    for (std::size_t basin = 0; basin < collection.size(); ++basin) {
        const auto outline = collection.basin(basin).get_basin_outline();
        std::vector<lsdtt::lsdtt_ridgeline_cpu> ridges;
        if (split_method == "longest_channel") ridges = outline.split_outline_into_ridgelines_by_longest_channel(distance);
        else if (split_method == "maximum_elevation") ridges = outline.split_outline_into_ridgelines_by_max_elevation(dem.native());
        else ridges = outline.split_outline_into_ridgelines_by_equal_length();
        for (std::size_t ridge_id = 0; ridge_id < ridges.size(); ++ridge_id) {
            auto& ridge = ridges[ridge_id];
            if (ridge.size() < 3) continue;
            ridge.fill_elevation(dem.native());
            ridge.segment_ridgeline_snout_by_curvature(snout_smoothing_fraction);
            ridge.fit_ridgeline_profile(false, false);
            ridge.smoothen_and_derive_detrended_ridgeline(detrend_smoothing_pixels);
            const auto extrema = ridge.get_ridgeline_extrema(flow, network.native());
            for (const auto& extremum : extrema) {
                if (extremum.node_index < 0 || static_cast<std::size_t>(extremum.node_index) >= ridge.positions().size()) continue;
                const auto& xy = ridge.positions()[static_cast<std::size_t>(extremum.node_index)];
                nb::dict row;
                row["basin_id"] = basin; row["ridgeline_id"] = ridge_id;
                row["node_index"] = extremum.node_index; row["type"] = extremum.type;
                row["distance_to_next_metres"] = extremum.distance_to_next;
                row["relief_to_next_metres"] = extremum.relief_to_next;
                row["distance_to_tdj_metres"] = extremum.distance_to_tdj;
                row["ridgeline_slope_fit"] = ridge.slope_fit();
                row["ridgeline_r_squared_fit"] = ridge.r_squared_fit();
                rows.append(std::move(row)); geometry.append(point(xy.x, xy.y));
            }
        }
    }
    nb::object table = nb::module_::import_("pandas").attr("DataFrame")(rows);
    nb::object ctor = nb::module_::import_("geopandas").attr("GeoDataFrame");
    return flow.meta().crs_wkt.empty() ? ctor(table, nb::arg("geometry") = geometry)
                                      : ctor(table, nb::arg("geometry") = geometry, nb::arg("crs") = flow.meta().crs_wkt);
}

nb::object basin_channel_lengths(const BasinCollection& collection, const ChannelNetwork& network, int max_stream_order)
{
    if (max_stream_order < 1) throw std::invalid_argument("max_stream_order must be positive");
    if (&network.flow() != &collection.flow())
        throw std::invalid_argument("network and basins must belong to the same FlowInfo");

    const auto& flow = collection.flow(); const auto& native_network = network.native();
    const auto& channel = native_network.channel_mask_array(); nb::list rows;
    for (std::size_t id = 0; id < collection.size(); ++id) {
        const auto& basin = collection.basin(id); std::unordered_set<std::int64_t> inside;
        for (const auto flat : basin.get_slice()) inside.insert(flat);
        std::vector<double> lengths(static_cast<std::size_t>(max_stream_order) + 1, 0.0);
        for (const auto flat : basin.get_slice()) {
            const auto index = static_cast<std::size_t>(flat);
            if (!channel.data()[index]) continue;
            const auto receiver = flow.receiver_of_unchecked(flat);
            if (receiver < 0 || inside.find(receiver) == inside.end() || !channel.data()[static_cast<std::size_t>(receiver)]) continue;
            const auto junction = native_network.junction_of_flat(flat);
            if (junction < 0) continue;
            const int order = native_network.stream_order_of_junction(junction);
            if (order < 1 || order > max_stream_order) continue;
            const int row = flow.row_from_flat_unchecked(flat), col = flow.col_from_flat_unchecked(flat);
            const int receiver_row = flow.row_from_flat_unchecked(receiver), receiver_col = flow.col_from_flat_unchecked(receiver);
            lengths[static_cast<std::size_t>(order)] += flow.cell_size() * std::hypot(
                static_cast<double>(receiver_row - row), static_cast<double>(receiver_col - col));
        }
        nb::dict row; row["basin_id"] = id;
        for (int order = 1; order <= max_stream_order; ++order) {
            const auto name = "channel_length_order_" + std::to_string(order) + "_metres";
            row[name.c_str()] = lengths[static_cast<std::size_t>(order)];
        }
        rows.append(std::move(row));
    }
    return nb::module_::import_("pandas").attr("DataFrame")(rows);
}

class SideTributaryBasins {
public:
    SideTributaryBasins(std::shared_ptr<lsdtt::lsdtt_flowinfo> flow, lsdtt::lsdtt_side_tributary_basin_result result)
        : flow_(std::move(flow)), result_(std::move(result))
    {
        std::vector<std::int64_t> outlets;
        for (const auto& basin : result_.accepted_basins) outlets.push_back(basin.outlet_flat);
        basins_ = std::make_unique<BasinCollection>(flow_, std::move(outlets));
    }
    const BasinCollection& basins() const noexcept { return *basins_; }
    const lsdtt::lsdtt_side_tributary_basin_result& native() const noexcept { return result_; }
    const lsdtt::lsdtt_flowinfo& flow() const noexcept { return *flow_; }
private:
    std::shared_ptr<lsdtt::lsdtt_flowinfo> flow_;
    lsdtt::lsdtt_side_tributary_basin_result result_;
    std::unique_ptr<BasinCollection> basins_;
};

SideTributaryBasins side_tributary_basins(const FlowInfo& flow, const std::string& points_path,
                                          int minimum_contributing_pixels, bool reject_incomplete,
                                          bool single_channel, bool include_stopping_point,
                                          int stopping_point_window_pixels)
{
    if (minimum_contributing_pixels < 1) throw std::invalid_argument("minimum_contributing_pixels must be positive");
    lsdtt::lsdtt_side_tributary_basin_options options;
    options.input_points_path = points_path; options.minimum_contributing_pixels = minimum_contributing_pixels;
    options.reject_incomplete_basins = reject_incomplete; options.include_stopping_point = include_stopping_point;
    options.stopping_point_window_pixels = stopping_point_window_pixels;
    auto result = single_channel ? lsdtt::extract_single_channel_side_basins(flow.native(), options)
                                 : lsdtt::extract_side_tributary_basins(flow.native(), options);
    return SideTributaryBasins(flow.shared_native(), std::move(result));
}

Raster side_channel_mask(const SideTributaryBasins& result)
{
    const auto& flow = result.flow(); std::vector<float> data(flow.full_grid_size(), 0.0f);
    for (std::size_t i = 0; i < data.size(); ++i) data[i] = result.native().traced_channel_mask[i] ? 1.0f : 0.0f;
    std::vector<std::uint8_t> mask(flow.mask_array().data(), flow.mask_array().data() + flow.full_grid_size());
    return Raster(lsdtt::lsdtt_raster_cpu(std::move(data), std::move(mask), flow.meta()));
}

nb::object side_tributary_table(const SideTributaryBasins& result)
{
    const auto& r = result.native(); std::vector<std::int64_t> flat, receiver; std::vector<int> row, col, receiver_row, receiver_col;
    std::vector<std::uint64_t> n; std::vector<double> area; std::vector<bool> complete;
    for (const auto& basin : r.accepted_basins) { flat.push_back(basin.outlet_flat); receiver.push_back(basin.receiver_channel_flat); row.push_back(basin.outlet_row); col.push_back(basin.outlet_col); receiver_row.push_back(basin.receiver_channel_row); receiver_col.push_back(basin.receiver_channel_col); n.push_back(basin.n_pixels); area.push_back(basin.drainage_area); complete.push_back(basin.complete); }
    nb::dict data; data["outlet_flat"] = flat; data["receiver_channel_flat"] = receiver; data["outlet_row"] = row; data["outlet_col"] = col; data["receiver_channel_row"] = receiver_row; data["receiver_channel_col"] = receiver_col; data["n_pixels"] = n; data["drainage_area_metres2"] = area; data["complete"] = complete;
    return nb::module_::import_("pandas").attr("DataFrame")(data);
}

BasinCollection line_crossing_basins(const ChannelNetwork& network, nb::object line, int minimum_contributing_pixels)
{
    if (minimum_contributing_pixels < 1) throw std::invalid_argument("minimum_contributing_pixels must be positive");
    nb::object geometry = line;
    if (nb::hasattr(line, "geometry")) geometry = line.attr("geometry").attr("unary_union");
    const auto& flow = network.flow(); std::vector<std::int64_t> flats;
    nb::object point = nb::module_::import_("shapely.geometry").attr("Point");
    const double tolerance = 0.71 * flow.cell_size();
    for (std::size_t i = 0; i < flow.full_grid_size(); ++i) {
        const auto flat = static_cast<std::int64_t>(i);
        if (!network.native().channel_mask_array().data()[i] || flow.s_vector_size(flat) < minimum_contributing_pixels) continue;
        const int row = flow.row_from_flat_unchecked(flat), col = flow.col_from_flat_unchecked(flat);
        const double x = flow.meta().x_min + (col + 0.5) * flow.cell_size(), y = flow.meta().y_max() - (row + 0.5) * flow.cell_size();
        if (nb::cast<double>(geometry.attr("distance")(point(x, y))) <= tolerance) flats.push_back(flat);
    }
    return BasinCollection(network.shared_flow(), std::move(flats));
}

nb::dict basin_surface_metrics(const BasinCollection& collection, const Raster& dem,
                               const std::vector<double>& scales, const std::vector<std::string>& metric_names,
                               bool include_outline, int n_threads)
{
    if (scales.empty() || metric_names.empty()) throw std::invalid_argument("scales and metrics must not be empty");
    if (dem.native().get_meta().n_rows != collection.flow().n_rows() || dem.native().get_meta().n_cols != collection.flow().n_cols())
        throw std::invalid_argument("dem must have the same grid shape as the basin flow topology");
    const std::map<std::string, lsdtt::polyfit_metric> lookup = {
        {"slope", lsdtt::polyfit_metric::slope}, {"slope_degrees", lsdtt::polyfit_metric::slope_degrees},
        {"maximum_curvature", lsdtt::polyfit_metric::maximum_curvature}, {"minimum_curvature", lsdtt::polyfit_metric::minimum_curvature},
        {"mean_curvature", lsdtt::polyfit_metric::mean_curvature}, {"gaussian_curvature", lsdtt::polyfit_metric::gaussian_curvature},
        {"profile_curvature", lsdtt::polyfit_metric::profile_curvature}, {"plan_curvature", lsdtt::polyfit_metric::plan_curvature},
        {"tangential_curvature", lsdtt::polyfit_metric::tangential_curvature}, {"aspect_degrees", lsdtt::polyfit_metric::aspect_degrees}};
    std::vector<lsdtt::polyfit_metric> metrics; for (const auto& name : metric_names) { const auto found = lookup.find(name); if (found == lookup.end()) throw std::invalid_argument("unsupported basin surface metric: " + name); metrics.push_back(found->second); }
    nb::list interior_rows, outline_rows, flip_rows, slope_area_rows, slope_area_summary_rows;
    for (const double scale : scales) {
        if (scale <= 0.0) throw std::invalid_argument("scales must be positive");
        const auto fit = lsdtt::direct_polyfit(dem.native(), scale, n_threads); lsdtt::polyfit_metric_request request; request.metrics = metrics; request.n_threads = n_threads;
        const auto result = lsdtt::compute_polyfit_metrics(fit, request);
        for (std::size_t id = 0; id < collection.size(); ++id) {
            const auto& basin = collection.basin(id);
            for (std::size_t metric_id = 0; metric_id < metrics.size(); ++metric_id) {
                const auto& raster = result.at(metrics[metric_id]); const auto values = basin.get_basin_values(raster); const auto stats = lsdtt::compute_descriptive_stats(values);
                nb::dict row; row["basin_id"] = id; row["scale_metres"] = scale; row["metric"] = metric_names[metric_id]; row["n_values"] = stats.n_values; row["min"] = stats.min; row["mean"] = stats.mean; row["median"] = stats.median; row["max"] = stats.max; interior_rows.append(std::move(row));
                if (include_outline) {
                    const auto outline = basin.get_basin_outline(); const auto sampled = outline.interpolate_pixel_values_on_outline(raster); double distance = 0.0; const auto& rc = outline.outline_divide_rc(); const double cell = dem.native().get_meta().cell_size; int flips = 0; float previous = std::numeric_limits<float>::quiet_NaN();
                    for (std::size_t node = 0; node < sampled.size(); ++node) { if (node > 0) distance += cell * std::hypot(static_cast<double>(rc[node].row - rc[node - 1].row), static_cast<double>(rc[node].col - rc[node - 1].col)); nb::dict sample; sample["basin_id"] = id; sample["scale_metres"] = scale; sample["metric"] = metric_names[metric_id]; sample["distance_metres"] = distance; sample["value"] = sampled[node]; outline_rows.append(std::move(sample)); if (metric_names[metric_id] == "maximum_curvature" && std::isfinite(previous) && std::isfinite(sampled[node]) && previous * sampled[node] < 0.0f) ++flips; if (std::isfinite(sampled[node])) previous = sampled[node]; }
                    if (metric_names[metric_id] == "maximum_curvature") { nb::dict f; f["basin_id"] = id; f["scale_metres"] = scale; f["n_flips"] = flips; f["mean_flip_spacing_metres"] = flips > 0 ? distance / flips : std::numeric_limits<double>::quiet_NaN(); flip_rows.append(std::move(f)); }
                }
            }

            // A compact, reproducible slope--area summary.  Area is the D8
            // contributing area at each basin pixel; the 0.1 log10(A) bins
            // keep the returned table useful without prescribing a plot.
            const auto slope_position = std::find(metric_names.begin(), metric_names.end(), "slope");
            if (slope_position != metric_names.end()) {
                const auto& slope_raster = result.at(lsdtt::polyfit_metric::slope);
                const auto& slope_data = slope_raster.get_data();
                const auto& slope_mask = slope_raster.get_mask();
                std::map<std::int64_t, std::vector<float>> bins;
                for (const auto flat : basin.get_slice()) {
                    const auto index = static_cast<std::size_t>(flat);
                    if (!slope_mask[index] || !std::isfinite(slope_data[index])) continue;
                    const double area = static_cast<double>(collection.flow().s_vector_size(flat))
                                      * collection.flow().cell_size() * collection.flow().cell_size();
                    if (area <= 0.0) continue;
                    bins[static_cast<std::int64_t>(std::floor(10.0 * std::log10(area)))].push_back(slope_data[index]);
                }
                double highest_median = -std::numeric_limits<double>::infinity();
                double area_at_highest_median = std::numeric_limits<double>::quiet_NaN();
                for (const auto& [bin, values] : bins) {
                    if (values.empty()) continue;
                    const double log10_area = (static_cast<double>(bin) + 0.5) / 10.0;
                    const double median_slope = lsdtt::compute_median(values);
                    nb::dict row; row["basin_id"] = id; row["scale_metres"] = scale;
                    row["log10_area_midpoint"] = log10_area; row["area_metres2"] = std::pow(10.0, log10_area);
                    row["n_pixels"] = values.size(); row["median_slope"] = median_slope;
                    slope_area_rows.append(std::move(row));
                    if (median_slope > highest_median) { highest_median = median_slope; area_at_highest_median = std::pow(10.0, log10_area); }
                }
                nb::dict summary; summary["basin_id"] = id; summary["scale_metres"] = scale;
                summary["area_at_max_median_slope_metres2"] = area_at_highest_median;
                summary["max_median_slope"] = std::isfinite(highest_median) ? highest_median : std::numeric_limits<double>::quiet_NaN();
                slope_area_summary_rows.append(std::move(summary));
            }
        }
    }
    nb::dict output; const auto pandas = nb::module_::import_("pandas"); output["interior"] = pandas.attr("DataFrame")(interior_rows); output["outline"] = pandas.attr("DataFrame")(outline_rows); output["curvature_flips"] = pandas.attr("DataFrame")(flip_rows); output["slope_area"] = pandas.attr("DataFrame")(slope_area_rows); output["slope_area_summary"] = pandas.attr("DataFrame")(slope_area_summary_rows); return output;
}

} // namespace

void bind_basins(nb::module_& module, nb::class_<FlowInfo>& flow)
{
    nb::class_<BasinCollection>(module, "BasinCollection", R"doc(
Selected D8 drainage basins (possibly nested) on a ``FlowInfo`` topology.

Made by ``FlowInfo.basins(outlets)``. Basin ``basin_id`` (0-based) is the
whole upslope area of its outlet pixel; basins may nest (a basin whose
outlet lies inside another basin). Each basin is stored as a view of the
flow topology, so collections are cheap. ``table`` describes each basin,
``polygons`` gives outlines, ``labels`` and ``mask()`` give rasters, and
``statistics(dem)`` summarises elevation, slope and shape.
)doc")
        .def_prop_ro("size", &BasinCollection::size, "Number of basins.")
        .def_prop_ro("labels", &basin_labels, R"doc(
Raster of basin membership: ``basin_id + 1``, 0 outside all basins.

Basins are painted in ``basin_id`` order, so a pixel inside several nested
basins carries the largest ``basin_id + 1`` among them. Nodata where the DEM
is nodata.
)doc")
        .def_prop_ro("table", &basin_table, nb::sig("def table(self) -> pandas.DataFrame"), R"doc(
One row per basin (pandas.DataFrame).

Columns: ``basin_id`` (0-based), ``outlet_flat``, ``outlet_row``,
``outlet_col``, ``n_pixels``, ``area_metres2`` (m^2), ``complete`` (False if
any basin pixel, including the outlet, is next to a nodata pixel or to an
open or base-level DEM edge, so the basin may be cut off), ``outlet_on_dem_edge``,
``parent_id`` (the smallest other selected basin containing this outlet,
-1 if none), ``nesting_depth`` (0 for a basin with no parent).
)doc")
        .def_prop_ro("polygons", &basin_polygons, nb::sig("def polygons(self) -> geopandas.GeoDataFrame"), R"doc(
Basin outline polygons as a GeoDataFrame.

Columns ``basin_id``, ``outlet_flat``, ``geometry`` (polygon traced along
pixel edges), in the DEM CRS.
)doc")
        .def("outline_ridgelines", &basin_outline_ridgelines, nb::arg("dem"),
             nb::arg("split_method") = "longest_channel", R"doc(
Split every selected basin outline into two ridgelines.

Parameters
----------
dem : Raster
    Elevation (m) on the flow grid.
split_method : {"longest_channel", "maximum_elevation", "equal_length"}, default="longest_channel"
    How the outline is split in two: at the
    head of the longest flow path, at the highest outline point, or halfway
    round the perimeter.

Returns
-------
geopandas.GeoDataFrame
    Left/right ridge lines with columns ``basin_id``, ``ridgeline_id``,
    ``outlet_flat``, ``n_nodes``, ``length_metres``, ``elevation_min``,
    ``elevation_max`` (m), ``geometry``.
)doc")
        .def("mask", &basin_mask, nb::arg("basin_id"), R"doc(
Return a 0/1 Raster of one basin.

Parameters
----------
basin_id : int
    Zero-based basin index (a ``table.basin_id`` value).

Returns
-------
Raster
    1 inside the basin, 0 on other valid pixels.
)doc")
        .def("pixels", &basin_pixels, nb::arg("basin_id"), R"doc(
Return the pixels of one basin.

Parameters
----------
basin_id : int
    Zero-based basin index.

Returns
-------
numpy.ndarray
    N x 2 int64 array of ``(row, col)``; the first row is the outlet.
)doc")
        .def("membership", &basin_membership, nb::arg("row"), nb::arg("col"), R"doc(
Return all selected basin IDs containing this pixel.

Parameters
----------
row, col : int
    Zero-based pixel position of a valid pixel.

Returns
-------
numpy.ndarray
    Int32 array of zero-based basin IDs (empty if the pixel is in none).
)doc")
        .def("statistics", &basin_statistics, nb::arg("dem"), nb::arg("include_slope") = true,
             nb::arg("include_flow_distance") = true,
             nb::arg("include_shape_metrics") = false, nb::arg("perimeter_simplify_cells") = 1.0,
             R"doc(
Return per-basin elevation, slope, flow-distance, and optional shape statistics.

Parameters
----------
dem : Raster
    Elevation (m) on the flow grid.
include_slope : bool, default=True
    Add ``slope_25th``, ``slope_50th``, ``slope_75th``: percentiles of the
    3 x 3 gradient (m/m, as ``slope()``).
include_flow_distance : bool, default=True
    Add ``flow_distance_25th/50th/75th``: percentiles of D8 flow distance
    to the outlet (m).
include_shape_metrics : bool, default=False
    Add planform shape metrics: ``centroid_x``, ``centroid_y``,
    ``basin_length_metres``, ``elongation_ratio``, ``form_factor``,
    ``circularity_ratio``, ``gravelius_compactness``, ``solidity``,
    ``eccentricity``, ``orientation_degrees`` (compass, 0-180),
    ``hypsometric_integral``.
perimeter_simplify_cells : float, default=1.0
    Simplification tolerance of the outline, in cells, used for the
    perimeter in the shape metrics.

Returns
-------
pandas.DataFrame
    One row per basin: ``basin_id``, ``n_pixels``, ``area_metres2``,
    ``elevation_min``, ``elevation_max``, ``elevation_mean``,
    ``elevation_median``, ``elevation_stddev`` (m), plus the optional
    columns above.
)doc");
    nb::class_<JunctionBasinCollection>(module, "JunctionBasinCollection", R"doc(
Junction-aware D8 basins selected from a ``ChannelNetwork``.

Each basin drains to one channel junction and carries that junction's
Strahler stream order. The ``filter_*`` methods remove basins in place and
return the collection, so they can be chained:
``junction_basins(N).filter_incomplete().filter_stream_order(2, 2)``.
Created with :func:`junction_basins`.
)doc")
        .def_prop_ro("size", &JunctionBasinCollection::size, "Number of selected basins.")
        .def_prop_ro("table", &junction_basin_table, nb::sig("def table(self) -> pandas.DataFrame"), R"doc(
One row per basin (pandas.DataFrame).

Columns ``basin_id`` (0-based), ``junction_id``, ``stream_order``,
``outlet_flat``, ``outlet_row``, ``outlet_col``, ``n_pixels``,
``area_metres2`` (m^2), ``complete`` (False if the basin touches nodata or
an open edge).
)doc")
        .def_prop_ro("labels", &junction_basin_labels,
            "Raster of each basin's junction ID; nested basins overwrite their parents; nodata outside all basins.")
        .def("mask", &junction_basin_mask, nb::arg("basin_id"), R"doc(
Return a 0/1 Raster of one basin.

Parameters
----------
basin_id : int
    Zero-based basin index (a ``table.basin_id`` value).

Returns
-------
Raster
)doc")
        .def("filter_incomplete", &JunctionBasinCollection::filter_incomplete, nb::rv_policy::reference_internal,
            "Remove basins that touch nodata or an open DEM edge; returns this collection.")
        .def("filter_area", &JunctionBasinCollection::filter_area, nb::arg("min_pixels") = -1, nb::arg("max_pixels") = -1, nb::rv_policy::reference_internal, R"doc(
Keep basins whose pixel count is within [min_pixels, max_pixels].

Parameters
----------
min_pixels, max_pixels : int, default=-1
    Bounds in pixels (multiply by ``metadata.pixel_area`` for m^2); -1
    leaves that end open.

Returns
-------
JunctionBasinCollection
    This collection, filtered in place.
)doc")
        .def("filter_stream_order", &JunctionBasinCollection::filter_stream_order, nb::arg("min_order") = -1, nb::arg("max_order") = -1, nb::rv_policy::reference_internal, R"doc(
Keep basins whose Strahler order is within [min_order, max_order].

Parameters
----------
min_order, max_order : int, default=-1
    Bounds on stream order; -1 leaves that end open.

Returns
-------
JunctionBasinCollection
    This collection, filtered in place.
)doc")
        .def("filter_nested", &JunctionBasinCollection::filter_nested, nb::rv_policy::reference_internal,
            "Remove basins whose outlet lies upslope of another selected basin's outlet; returns this collection.")
        .def("filter_nested_by_stream_order", &JunctionBasinCollection::filter_nested_by_stream_order, nb::rv_policy::reference_internal,
            "Remove nested basins only where they share the enclosing basin's stream order; returns this collection.")
        .def("filter_largest", &JunctionBasinCollection::filter_largest, nb::rv_policy::reference_internal,
            "Keep only the basin with the most pixels; returns this collection.");
    module.def("junction_basins", &make_junction_basins, nb::arg("network"), nb::arg("selection") = "all", nb::arg("outlets") = nb::none(), nb::keep_alive<0, 1>(), R"doc(
Select junction-aware basins from a channel network.

Parameters
----------
network : ChannelNetwork
    Channel network whose junctions define the basins.
selection : {"all", "outlets", "nested"}, default="all"
    ``"all"``: one basin per junction, draining to the pixel just upstream
    of the next junction downstream. ``"outlets"``: one basin per supplied
    outlet. ``"nested"``: each supplied outlet's basin plus every junction
    basin nested inside it.
outlets : array-like of int or None, default=None
    N x 2 array of ``(row, col)`` outlet pixels (valid pixels); required
    for ``"outlets"`` and ``"nested"``, ignored for ``"all"``.

Returns
-------
JunctionBasinCollection
    Basins that can be filtered by area, stream order and nesting.
)doc");
    nb::class_<SideTributaryBasins>(module, "SideTributaryBasins", R"doc(
Side basins draining into a traced D8 channel.

Made by :func:`side_tributary_basins`.
)doc")
        .def_prop_ro("table", &side_tributary_table, nb::sig("def table(self) -> pandas.DataFrame"),
            "One row per accepted side basin: outlet_flat, receiver_channel_flat, outlet_row, outlet_col, receiver_channel_row, receiver_channel_col, n_pixels, drainage_area_metres2, complete.")
        .def_prop_ro("channel_mask", &side_channel_mask, "0/1 Raster of the traced trunk channel.")
        .def_prop_ro("basins", &SideTributaryBasins::basins, nb::rv_policy::reference_internal,
            "The accepted side basins as a BasinCollection.");
    module.def("side_tributary_basins", &side_tributary_basins, nb::arg("flow"), nb::arg("points_path"),
        nb::arg("minimum_contributing_pixels") = 1000, nb::arg("reject_incomplete") = true,
        nb::arg("single_channel") = false, nb::arg("include_stopping_point") = false,
        nb::arg("stopping_point_window_pixels") = 5, R"doc(
Extract the side basins that drain into channels traced from input points.

Parameters
----------
flow : FlowInfo
    Flow topology.
points_path : str
    Vector file of the points from which trunk channels are traced D8
    downstream.
minimum_contributing_pixels : int, default=1000
    Smallest side basin kept, in contributing pixels (positive).
reject_incomplete : bool, default=True
    Drop side basins that touch nodata or an open DEM edge.
single_channel : bool, default=False
    Trace a single main channel: the file holds one point (traced to the
    DEM edge) or two points (traced from the upstream one to where the
    other's flow path joins it). If False every point is a source.
include_stopping_point : bool, default=False
    Treat the last point in the file as a stopping point: each trace ends
    when it first enters a square window around it.
stopping_point_window_pixels : int, default=5
    Half-width (pixels) of the stopping-point window.

Returns
-------
SideTributaryBasins
)doc");
    module.def("line_crossing_basins", &line_crossing_basins, nb::arg("network"), nb::arg("line"),
        nb::arg("minimum_contributing_pixels") = 500, R"doc(
Create basins whose outlets are channel pixels crossed by a line.

Every network channel pixel within 0.71 cell sizes of the line, with at least
``minimum_contributing_pixels`` contributing pixels, becomes an outlet.

Parameters
----------
network : ChannelNetwork
    Channel network on the line's grid.
line : shapely geometry or geopandas.GeoDataFrame
    Line (for example a mountain front or a road) in the DEM CRS; a
    GeoDataFrame is merged into one geometry.
minimum_contributing_pixels : int, default=500
    Smallest channel crossing kept, in contributing pixels (positive).

Returns
-------
BasinCollection
)doc");
    module.def("basin_channel_lengths", &basin_channel_lengths, nb::arg("basins"), nb::arg("network"),
        nb::arg("max_stream_order") = 6, R"doc(
Return total in-basin D8 channel length by Strahler stream order.

``basins`` and ``network`` must have been created from the same ``FlowInfo``.

Parameters
----------
basins : BasinCollection
    Basins to measure.
network : ChannelNetwork
    Channel network from the same ``FlowInfo``.
max_stream_order : int, default=6
    Highest Strahler order reported (positive).

Returns
-------
pandas.DataFrame
    One row per basin: ``basin_id`` and ``channel_length_order_<k>_metres``
    for k = 1 .. ``max_stream_order`` (D8 channel length in m).
)doc");
    module.def("basin_surface_metrics", &basin_surface_metrics, nb::arg("basins"), nb::arg("dem"),
        nb::arg("scales"), nb::arg("metrics"), nb::arg("include_outline") = true, nb::arg("n_threads") = 0, R"doc(
Summarise polyfit surface metrics inside each basin at several fitting scales.

Parameters
----------
basins : BasinCollection
    Basins to summarise.
dem : Raster
    Elevation (m) on the basins' flow grid.
scales : list of float
    Polyfit length scales ``lambda_metres`` (m, positive), as in
    ``polyfit_metrics``.
metrics : list of str
    Any of ``"slope"`` (m/m), ``"slope_degrees"``, ``"aspect_degrees"``,
    ``"maximum_curvature"``, ``"minimum_curvature"``, ``"mean_curvature"``,
    ``"gaussian_curvature"``, ``"profile_curvature"``, ``"plan_curvature"``,
    ``"tangential_curvature"``.
include_outline : bool, default=True
    Also sample each metric along the basin outline.
n_threads : int, default=0
    OpenMP thread count; 0 uses the runtime default.

Returns
-------
dict of str to pandas.DataFrame
    ``"interior"`` (basin_id, scale_metres, metric, n_values, min, mean,
    median, max), ``"outline"`` (basin_id, scale_metres, metric,
    distance_metres, value), ``"curvature_flips"`` (sign changes of
    maximum_curvature along the outline), ``"slope_area"`` (median slope in
    0.1 log10 drainage-area bins) and ``"slope_area_summary"``.
)doc");
    module.def("basin_ridgeline_extrema", &basin_ridgeline_extrema,
        nb::arg("basins"), nb::arg("network"), nb::arg("dem"),
        nb::arg("split_method") = "longest_channel",
        nb::arg("snout_smoothing_fraction") = 0.1f,
        nb::arg("detrend_smoothing_pixels") = 5, R"doc(
Locate topology-aware maxima and minima along selected basin-outline ridgelines.

Parameters
----------
basins : BasinCollection
    Basin outlines to split into ridgelines.
network : ChannelNetwork
    Same-flow topology used to identify triple/quadruple divide junctions.
dem : Raster
    Elevation raster on the shared flow grid.
split_method : {"longest_channel", "maximum_elevation", "equal_length"}, optional
    Rule used to divide each basin outline into ridgelines.
snout_smoothing_fraction : float, optional
    Smoothing window, as a fraction of ridge length, used to find the ridge
    snout from its curvature.
detrend_smoothing_pixels : int, optional
    Smoothing window (pixels) used to detrend each ridge profile.

Returns
-------
geopandas.GeoDataFrame
    Extrema points with ridge fitting and topological-distance diagnostics.

Authors
-------
B. Gailleton
)doc");
    flow.def("basins", &make_basins, nb::arg("outlets") = nb::none(), nb::arg("min_upslope_pixels") = 0,
        R"doc(
Create a ``BasinCollection`` from explicit outlet pixels or all qualifying D8 outlets.

Explicit outlets may be any valid pixels, including channel pixels inside
another basin, which creates nested basins.

Parameters
----------
outlets : array-like of int or None, default=None
    N x 2 integer array of ``(row, col)`` outlet pixels; basin ``i`` is the
    upslope area of the ``i``-th distinct outlet. Duplicate outlets are
    silently dropped (so IDs follow first appearance). Out-of-grid or
    nodata pixels raise ``ValueError``. None selects every D8 outlet of the
    flow topology (pixels that drain off the grid or into nodata).
min_upslope_pixels : int, default=0
    With ``outlets=None``, keep only outlets with at least this many
    contributing pixels (including the outlet). Ignored when ``outlets`` is
    given. Non-negative.

Returns
-------
BasinCollection
)doc");
}
} // namespace lsdtt3_python
