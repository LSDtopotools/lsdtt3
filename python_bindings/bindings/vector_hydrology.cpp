#include <algorithm>
#include <cstdint>
#include <cstring>
#include <stdexcept>
#include <string>
#include <variant>
#include <vector>

#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>

#include <cpl_vsi.h>
#include <ogr_spatialref.h>

#include <cpu/lsdtt_burn_streams_cpu.hpp>
#include <cpu/lsdtt_vector_line_repair_cpu.hpp>
#include <cpu/lsdtt_water_area_centreline_cpu.hpp>
#include <utils/lsdtt_raster_geometry.hpp>

#include "binding_common.hpp"

namespace lsdtt3_python {
namespace {

std::string geodataframe_wkt(nb::handle frame, const char* name)
{
    nb::object crs = frame.attr("crs");
    if (crs.is_none()) throw std::invalid_argument(std::string(name) + " must define a CRS");
    return nb::cast<std::string>(crs.attr("to_wkt")());
}

void require_same_crs(const std::string& a, const std::string& b, const char* message)
{
    OGRSpatialReference first, second;
    if (a.empty() || b.empty() || first.SetFromUserInput(a.c_str()) != OGRERR_NONE ||
        second.SetFromUserInput(b.c_str()) != OGRERR_NONE || !first.IsSame(&second))
        throw std::invalid_argument(message);
}

std::string geodataframe_geojson(nb::handle frame, const char* name)
{
    (void)geodataframe_wkt(frame, name);
    const std::string geojson = nb::cast<std::string>(frame.attr("to_json")());
    static std::uint64_t counter = 0;
    const std::string path = "/vsimem/lsdtt3-vector-hydrology-" + std::to_string(++counter) + ".geojson";
    auto* bytes = static_cast<GByte*>(CPLMalloc(geojson.size()));
    if (!bytes) throw std::bad_alloc();
    std::memcpy(bytes, geojson.data(), geojson.size());
    if (!VSIFileFromMemBuffer(path.c_str(), bytes, static_cast<vsi_l_offset>(geojson.size()), TRUE)) {
        CPLFree(bytes);
        throw std::runtime_error("could not create in-memory GeoJSON");
    }
    return path;
}

lsdtt::lsdtt_vector_layer read_vector_layer(const std::string& path)
{
    lsdtt::lsdtt_data_loader loader;
    loader.open_vector(path);
    try {
        const auto names = loader.get_layer_names();
        if (names.empty()) throw std::runtime_error("GeoDataFrame has no vector layer");
        auto layer = loader.read_vector_layer(names.front());
        loader.close_vector();
        return layer;
    } catch (...) {
        loader.close_vector();
        throw;
    }
}

nb::object line_geometry(const lsdtt::lsdtt_geom& geometry)
{
    nb::object shapely = nb::module_::import_("shapely.geometry");
    const auto points_to_coordinates = [](const std::vector<lsdtt::lsdtt_xy>& points) {
        nb::list coordinates;
        for (const auto& point : points) coordinates.append(nb::make_tuple(point.x, point.y));
        return coordinates;
    };
    if (const auto* line = std::get_if<lsdtt::lsdtt_line_geom>(&geometry))
        return shapely.attr("LineString")(points_to_coordinates(line->points));
    if (const auto* multiline = std::get_if<lsdtt::lsdtt_multiline_geom>(&geometry)) {
        nb::list parts;
        for (const auto& part : multiline->parts) parts.append(points_to_coordinates(part));
        return shapely.attr("MultiLineString")(parts);
    }
    return nb::none();
}

nb::object layer_geodataframe(const lsdtt::lsdtt_vector_layer& layer, const std::string& crs_wkt)
{
    nb::dict data;
    std::vector<nb::list> columns(layer.schema.size());
    nb::list geometry;
    for (const auto& record : layer.records) {
        geometry.append(line_geometry(record.geometry));
        for (std::size_t i = 0; i < columns.size(); ++i) {
            if (i >= record.field_values.size() || std::holds_alternative<std::monostate>(record.field_values[i])) columns[i].append(nb::none());
            else if (const auto* value = std::get_if<std::int64_t>(&record.field_values[i])) columns[i].append(*value);
            else if (const auto* value = std::get_if<double>(&record.field_values[i])) columns[i].append(*value);
            else columns[i].append(std::get<std::string>(record.field_values[i]));
        }
    }
    for (std::size_t i = 0; i < columns.size(); ++i) data[nb::str(layer.schema.field_names[i].c_str())] = columns[i];
    return nb::module_::import_("geopandas").attr("GeoDataFrame")(data, nb::arg("geometry") = geometry, nb::arg("crs") = crs_wkt);
}

int infer_epsg(nb::handle frame, int epsg)
{
    if (epsg > 0) return epsg;
    nb::object inferred = frame.attr("crs").attr("to_epsg")();
    if (inferred.is_none()) throw std::invalid_argument("could not infer EPSG from CRS; pass epsg explicitly");
    return nb::cast<int>(inferred);
}

struct WaterCentrelineResult {
    nb::object lines;
    int n_polygons = 0, n_lochs = 0, n_centrelines = 0, n_connectors = 0, n_tributaries = 0;
};

struct WaterLineRepairResult {
    nb::object lines, tips, bridges;
    int n_dangling_tips = 0, n_candidate_pairs = 0, n_bridges_accepted = 0;
};

struct StreamBurnResult {
    Raster raster;
    nb::object profiles;
};

WaterCentrelineResult water_centrelines(nb::object areas, nb::object tributaries, int epsg,
                                        double cell_size_metres, double padding_metres,
                                        double tributary_snap_tolerance_metres,
                                        double centreline_snap_max_metres,
                                        double simplify_tolerance_metres,
                                        double min_component_length_metres)
{
    const std::string area_wkt = geodataframe_wkt(areas, "areas");
    require_same_crs(area_wkt, geodataframe_wkt(tributaries, "tributaries"), "areas and tributaries must share a CRS");
    const std::string areas_path = geodataframe_geojson(areas, "areas");
    const std::string tributaries_path = geodataframe_geojson(tributaries, "tributaries");
    try {
        lsdtt::water_centreline_options options;
        options.epsg = infer_epsg(areas, epsg);
        options.cell_size_m = cell_size_metres;
        options.padding_m = padding_metres;
        options.tributary_snap_tolerance_m = tributary_snap_tolerance_metres;
        options.centreline_snap_max_m = centreline_snap_max_metres;
        options.dp_tolerance_m = simplify_tolerance_metres;
        options.min_centreline_component_m = min_component_length_metres;
        const auto result = lsdtt::extract_water_area_centrelines(read_vector_layer(areas_path), read_vector_layer(tributaries_path), options);
        VSIUnlink(areas_path.c_str()); VSIUnlink(tributaries_path.c_str());
        return {layer_geodataframe(result.output_layer, area_wkt), result.n_polygons, result.n_lochs,
                result.n_centrelines, result.n_connectors, result.n_tributaries};
    } catch (...) { VSIUnlink(areas_path.c_str()); VSIUnlink(tributaries_path.c_str()); throw; }
}

WaterLineRepairResult repair_water_lines(nb::object water_lines, nb::object roads,
                                         double snap_tolerance_metres, double max_gap_length_metres,
                                         double max_nonroad_gap_metres, double road_intersection_tolerance_metres,
                                         bool require_mutual_best)
{
    const std::string water_wkt = geodataframe_wkt(water_lines, "water_lines");
    require_same_crs(water_wkt, geodataframe_wkt(roads, "roads"), "water_lines and roads must share a CRS");
    const std::string water_path = geodataframe_geojson(water_lines, "water_lines");
    const std::string roads_path = geodataframe_geojson(roads, "roads");
    try {
        lsdtt::water_line_repair_options options;
        options.snap_tolerance_metres = snap_tolerance_metres;
        options.max_gap_length_metres = max_gap_length_metres;
        options.max_nonroad_gap_metres = max_nonroad_gap_metres;
        options.road_intersection_tolerance_metres = road_intersection_tolerance_metres;
        options.require_mutual_best = require_mutual_best;
        const auto result = lsdtt::repair_water_lines(read_vector_layer(water_path), read_vector_layer(roads_path), options);
        VSIUnlink(water_path.c_str()); VSIUnlink(roads_path.c_str());
        nb::list tip_line, tip_start, tip_x, tip_y;
        for (const auto& tip : result.tips) { tip_line.append(tip.line_index); tip_start.append(tip.at_start); tip_x.append(tip.xy.x); tip_y.append(tip.xy.y); }
        nb::dict tip_data; tip_data["line_index"] = tip_line; tip_data["at_start"] = tip_start; tip_data["x"] = tip_x; tip_data["y"] = tip_y;
        nb::list a, b, length, score, crosses_road, road_name;
        for (const auto& bridge : result.accepted_bridges) { a.append(bridge.tip_index_a); b.append(bridge.tip_index_b); length.append(bridge.gap_length_metres); score.append(bridge.score); crosses_road.append(bridge.crosses_road); road_name.append(bridge.crossed_road_name); }
        nb::dict bridge_data; bridge_data["tip_index_a"] = a; bridge_data["tip_index_b"] = b; bridge_data["gap_length_metres"] = length; bridge_data["score"] = score; bridge_data["crosses_road"] = crosses_road; bridge_data["crossed_road_name"] = road_name;
        nb::object pandas = nb::module_::import_("pandas");
        return {layer_geodataframe(result.output_layer, water_wkt), pandas.attr("DataFrame")(tip_data), pandas.attr("DataFrame")(bridge_data),
                result.n_dangling_tips, result.n_candidate_pairs, result.n_bridges_accepted};
    } catch (...) { VSIUnlink(water_path.c_str()); VSIUnlink(roads_path.c_str()); throw; }
}

StreamBurnResult burn_streams(const Raster& dem, nb::object lines, int n_threads)
{
    const auto& native = dem.native();
    const std::string line_wkt = geodataframe_wkt(lines, "lines");
    require_same_crs(native.get_meta().crs_wkt, line_wkt, "lines CRS must match the DEM CRS; reproject explicitly first");
    const std::string path = geodataframe_geojson(lines, "lines");
    try {
        const auto profiles = lsdtt::burn_streams_in_place(native, path, n_threads);
        VSIUnlink(path.c_str());
        auto data = native.get_data(); auto mask = native.get_mask(); const auto meta = native.get_meta();
        nb::list feature_id, flow_distance, original_elevation, burned_elevation, x, y;
        for (const auto& profile : profiles) for (const auto& node : profile.nodes) {
            feature_id.append(profile.feature_id); flow_distance.append(node.flow_distance); original_elevation.append(node.dem_elevation); burned_elevation.append(node.burned_elevation); x.append(node.x); y.append(node.y);
            if (!lsdtt::point_within_raster(meta, node.x, node.y)) continue;
            const int row = lsdtt::row_from_y(meta, node.y), col = lsdtt::col_from_x(meta, node.x);
            const std::size_t flat = static_cast<std::size_t>(row) * meta.n_cols + col;
            if (mask[flat]) data[flat] = std::min(data[flat], static_cast<float>(node.burned_elevation));
        }
        nb::dict table; table["feature_id"] = feature_id; table["flow_distance"] = flow_distance; table["original_elevation"] = original_elevation; table["burned_elevation"] = burned_elevation; table["x"] = x; table["y"] = y;
        return {Raster(lsdtt::lsdtt_raster_cpu(std::move(data), std::move(mask), meta)), nb::module_::import_("pandas").attr("DataFrame")(table)};
    } catch (...) { VSIUnlink(path.c_str()); throw; }
}

} // namespace

void bind_vector_hydrology(nb::module_& module)
{
    nb::class_<WaterCentrelineResult>(module, "WaterCentrelineResult", "Vector water-area centreline output.\n\nAuthors\n-------\nB. Gailleton")
        .def_ro("lines", &WaterCentrelineResult::lines, "Centreline, connector, and tributary GeoDataFrame.")
        .def_prop_ro("summary", [](const WaterCentrelineResult& self) { nb::dict d; d["n_polygons"] = self.n_polygons; d["n_lochs"] = self.n_lochs; d["n_centrelines"] = self.n_centrelines; d["n_connectors"] = self.n_connectors; d["n_tributaries"] = self.n_tributaries; return d; });
    nb::class_<WaterLineRepairResult>(module, "WaterLineRepairResult", "Repaired water-line output and diagnostics.\n\nAuthors\n-------\nB. Gailleton")
        .def_ro("lines", &WaterLineRepairResult::lines).def_ro("tips", &WaterLineRepairResult::tips).def_ro("bridges", &WaterLineRepairResult::bridges)
        .def_prop_ro("summary", [](const WaterLineRepairResult& self) { nb::dict d; d["n_dangling_tips"] = self.n_dangling_tips; d["n_candidate_pairs"] = self.n_candidate_pairs; d["n_bridges_accepted"] = self.n_bridges_accepted; return d; });
    nb::class_<StreamBurnResult>(module, "StreamBurnResult", "DEM conditioned along mapped stream profiles.\n\nAuthors\n-------\nB. Gailleton")
        .def_ro("raster", &StreamBurnResult::raster, "Copied DEM with channel cells lowered where required.")
        .def_ro("profiles", &StreamBurnResult::profiles, "Per-profile monotone burn diagnostics.");

    module.def("water_centrelines", &water_centrelines,
        nb::arg("areas"), nb::arg("tributaries"), nb::arg("epsg") = 0, nb::arg("cell_size_metres") = 1.0,
        nb::arg("padding_metres") = 5.0, nb::arg("tributary_snap_tolerance_metres") = 2.0,
        nb::arg("centreline_snap_max_metres") = 100.0, nb::arg("simplify_tolerance_metres") = 1.0,
        nb::arg("min_component_length_metres") = 10.0,
        "Extract central flow paths from water polygons and connect tributary lines.\n\nParameters\n----------\nareas, tributaries : geopandas.GeoDataFrame\n    Projected water polygons and linework in the same CRS.\nepsg : int, default=0\n    Output CRS identifier; zero infers it from areas.\ncell_size_metres : float, default=1\n    Internal skeletonisation-grid cell size.\n\nReturns\n-------\nWaterCentrelineResult\n    Connected linework and count diagnostics.\n\nAuthors\n-------\nB. Gailleton");
    module.def("repair_water_lines", &repair_water_lines,
        nb::arg("water_lines"), nb::arg("roads"), nb::arg("snap_tolerance_metres") = 0.5,
        nb::arg("max_gap_length_metres") = 25.0, nb::arg("max_nonroad_gap_metres") = 30.0,
        nb::arg("road_intersection_tolerance_metres") = 1.0, nb::arg("require_mutual_best") = true,
        "Repair short broken water-line gaps, preferentially across roads.\n\nParameters\n----------\nwater_lines, roads : geopandas.GeoDataFrame\n    Projected line layers in the same CRS.\nmax_gap_length_metres : float, default=25\n    Maximum road-associated bridge length.\n\nReturns\n-------\nWaterLineRepairResult\n    Original linework plus accepted bridge segments and diagnostics.\n\nAuthors\n-------\nB. Gailleton");
    module.def("burn_streams", &burn_streams, nb::arg("dem"), nb::arg("lines"), nb::arg("n_threads") = 0,
        "Condition a copied DEM along mapped stream lines with monotone profiles.\n\nParameters\n----------\ndem : Raster\n    Elevation raster.\nlines : geopandas.GeoDataFrame\n    Stream LineStrings in the DEM CRS.\nn_threads : int, default=0\n    OpenMP thread count; zero uses the runtime default.\n\nReturns\n-------\nStreamBurnResult\n    Conditioned raster and profile diagnostics.\n\nNotes\n-----\nOnly cells whose monotone stream elevation is lower are changed.\n\nAuthors\n-------\nB. Gailleton");
}

} // namespace lsdtt3_python
