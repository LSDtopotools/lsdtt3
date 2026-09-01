#include <cstdint>
#include <stdexcept>
#include <string>
#include <utility>
#include <variant>
#include <vector>

#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>

#include <cpu/lsdtt_internal_drainage_cpu.hpp>
#include <cpu/lsdtt_depression_hierarchy_cpu.hpp>

#include "binding_common.hpp"

namespace lsdtt3_python {
namespace {

nb::object dataframe(nb::list rows)
{
    return nb::module_::import_("pandas").attr("DataFrame")(rows);
}

nb::object geodataframe(nb::object data, nb::list geometry,
                        const lsdtt::lsdtt_raster_meta& meta)
{
    nb::object make = nb::module_::import_("geopandas").attr("GeoDataFrame");
    if (meta.crs_wkt.empty()) return make(data, nb::arg("geometry") = geometry);
    return make(data, nb::arg("geometry") = geometry, nb::arg("crs") = meta.crs_wkt);
}

nb::object polygon_geometry(const lsdtt::lsdtt_geom& geometry)
{
    const auto* polygon = std::get_if<lsdtt::lsdtt_polygon_geom>(&geometry);
    if (!polygon) return nb::none();
    nb::list exterior, holes;
    for (const auto& point : polygon->exterior) exterior.append(nb::make_tuple(point.x, point.y));
    for (const auto& hole : polygon->holes) {
        nb::list ring;
        for (const auto& point : hole) ring.append(nb::make_tuple(point.x, point.y));
        holes.append(ring);
    }
    return nb::module_::import_("shapely.geometry").attr("Polygon")(exterior, nb::arg("holes") = holes);
}

nb::object candidates_gdf(const std::vector<lsdtt::sink_candidate>& sinks,
                          const lsdtt::lsdtt_raster_meta& meta)
{
    nb::list rows, geometry;
    const std::vector<const char*> columns = {
        "sink_id", "routing_sink_id", "pit_flat", "pit_row", "pit_col",
        "latitude", "longitude", "pit_elevation", "spill_elevation", "depth_m",
        "area_m2", "volume_m3", "n_pixels", "tagged_drainage_area_m2",
        "tagged_n_pixels", "hierarchy_level", "parent_sink_id", "accepted"
    };
    nb::object point = nb::module_::import_("shapely.geometry").attr("Point");
    for (const auto& sink : sinks) {
        nb::dict row;
        row["sink_id"] = sink.sink_id; row["routing_sink_id"] = sink.routing_sink_id;
        row["pit_flat"] = sink.pit_flat; row["pit_row"] = sink.pit_row; row["pit_col"] = sink.pit_col;
        row["latitude"] = sink.latitude; row["longitude"] = sink.longitude;
        row["pit_elevation"] = sink.pit_elevation; row["spill_elevation"] = sink.spill_elevation;
        row["depth_m"] = sink.depth_m; row["area_m2"] = sink.area_m2; row["volume_m3"] = sink.volume_m3;
        row["n_pixels"] = sink.n_pixels; row["tagged_drainage_area_m2"] = sink.tagged_drainage_area_m2;
        row["tagged_n_pixels"] = sink.tagged_n_pixels; row["hierarchy_level"] = sink.hierarchy_level;
        row["parent_sink_id"] = sink.parent_sink_id; row["accepted"] = sink.accepted;
        rows.append(std::move(row)); geometry.append(point(sink.pit_x, sink.pit_y));
    }
    nb::object table = dataframe(rows);
    if (sinks.empty()) {
        nb::list names;
        for (const char* column : columns) names.append(column);
        table = table.attr("reindex")(nb::arg("columns") = names);
    }
    return geodataframe(table, geometry, meta);
}

nb::object polygons_gdf(const lsdtt::depression_polygons_result& polygons,
                        const lsdtt::lsdtt_raster_meta& meta)
{
    nb::list rows, geometry;
    for (const auto& record : polygons.records) {
        nb::dict row;
        for (std::size_t i = 0; i < polygons.schema.size(); ++i) {
            const auto key = nb::str(polygons.schema.field_names[i].c_str());
            if (i >= record.field_values.size() || std::holds_alternative<std::monostate>(record.field_values[i])) row[key] = nb::none();
            else if (const auto* value = std::get_if<std::int64_t>(&record.field_values[i])) row[key] = *value;
            else if (const auto* value = std::get_if<double>(&record.field_values[i])) row[key] = *value;
            else row[key] = std::get<std::string>(record.field_values[i]);
        }
        rows.append(std::move(row)); geometry.append(polygon_geometry(record.geometry));
    }
    return geodataframe(dataframe(rows), geometry, meta);
}

struct InternalDrainageResult {
    nb::object candidates;
    nb::object polygons;
    int n_candidates = 0;
};

class DepressionHierarchyResult {
public:
    DepressionHierarchyResult(lsdtt::lsdtt_raster_cpu labels,
                              lsdtt::lsdtt_raster_cpu flow_direction,
                              nb::object nodes)
        : labels_(std::move(labels)), flow_direction_(std::move(flow_direction)),
          nodes_(std::move(nodes)) {}
    Raster labels() const { return Raster(labels_); }
    Raster flow_direction() const { return Raster(flow_direction_); }
    nb::object nodes() const { return nodes_; }
    int n_nodes() const { return nb::cast<int>(nodes_.attr("shape")[0]); }
private:
    lsdtt::lsdtt_raster_cpu labels_, flow_direction_;
    nb::object nodes_;
};

DepressionHierarchyResult depression_hierarchy(const Raster& dem,
                                               const std::string& boundary_conditions,
                                               bool ocean_from_boundary,
                                               bool nodata_is_ocean,
                                               nb::object ocean_level)
{
    if (boundary_conditions.size() != 4)
        throw std::invalid_argument("boundary_conditions must contain four NESW codes");
    lsdtt::depression_hierarchy_options options;
    options.boundary_conditions = boundary_conditions;
    options.ocean_from_boundary = ocean_from_boundary;
    options.nodata_is_ocean = nodata_is_ocean;
    if (!ocean_level.is_none()) options.ocean_level = nb::cast<float>(ocean_level);

    std::vector<lsdtt::dep_label_t> labels;
    std::vector<std::uint8_t> directions;
    auto hierarchy = lsdtt::build_depression_hierarchy(dem.native(), labels, directions, options);
    const auto& meta = dem.native().get_meta();
    const auto& valid = dem.native().get_mask();
    std::vector<float> label_values(labels.size()), direction_values(directions.size());
    for (std::size_t i = 0; i < labels.size(); ++i) {
        label_values[i] = static_cast<float>(labels[i]);
        direction_values[i] = directions[i] == lsdtt::k_dh_no_flow
            ? static_cast<float>(meta.no_data_value) : static_cast<float>(directions[i]);
    }
    std::vector<std::uint8_t> direction_mask(valid);
    for (std::size_t i = 0; i < directions.size(); ++i)
        if (directions[i] == lsdtt::k_dh_no_flow) direction_mask[i] = 0;

    nb::list rows;
    for (const auto& node : hierarchy) {
        nb::dict row;
        row["depression_id"] = node.dep_label; row["pit_flat"] = node.pit_cell;
        row["outlet_flat"] = node.out_cell; row["parent_id"] = node.parent == lsdtt::k_no_parent ? -1 : static_cast<std::int64_t>(node.parent);
        row["overflow_depression_id"] = node.odep == lsdtt::k_no_label ? -1 : static_cast<std::int64_t>(node.odep);
        row["geographic_overflow_id"] = node.geolink == lsdtt::k_no_label ? -1 : static_cast<std::int64_t>(node.geolink);
        row["left_child_id"] = node.lchild == lsdtt::k_no_label ? -1 : static_cast<std::int64_t>(node.lchild);
        row["right_child_id"] = node.rchild == lsdtt::k_no_label ? -1 : static_cast<std::int64_t>(node.rchild);
        row["pit_elevation"] = node.pit_elev; row["spill_elevation"] = node.out_elev;
        row["n_pixels"] = node.cell_count; row["volume_m3"] = node.dep_vol;
        row["water_volume_m3"] = node.water_vol; row["ocean_parent"] = node.ocean_parent;
        row["is_leaf"] = node.lchild == lsdtt::k_no_label && node.rchild == lsdtt::k_no_label;
        rows.append(std::move(row));
    }
    return DepressionHierarchyResult(
        lsdtt::lsdtt_raster_cpu(std::move(label_values), std::vector<std::uint8_t>(valid), meta),
        lsdtt::lsdtt_raster_cpu(std::move(direction_values), std::move(direction_mask), meta),
        dataframe(rows));
}

InternalDrainageResult find_internal_drainage(const Raster& dem, const std::string& method,
                                              const std::string& boundary_conditions,
                                              int min_pixels, float min_depth_metres,
                                              double min_area_metres2, double min_volume_metres3,
                                              bool leaves_only, bool include_polygons, int n_threads)
{
    lsdtt::sink_detection_options options;
    if (method == "hierarchy") options.method = lsdtt::sink_detection_method::depression_hierarchy;
    else if (method == "fill_difference") options.method = lsdtt::sink_detection_method::fill_difference;
    else throw std::invalid_argument("method must be 'hierarchy' or 'fill_difference'");
    options.boundary_conditions = boundary_conditions;
    options.min_pixels_in_patch = min_pixels;
    options.min_depth_m = min_depth_metres;
    options.min_area_m2 = min_area_metres2;
    options.min_volume_m3 = min_volume_metres3;
    options.leaves_only = leaves_only;
    options.n_threads = n_threads;
    const auto sinks = lsdtt::find_internal_drainage(dem.native(), options);
    const auto& meta = dem.native().get_meta();
    nb::object polygons = nb::none();
    if (include_polygons) polygons = polygons_gdf(lsdtt::sink_depressions_to_polygons(dem.native(), sinks, options), meta);
    return {candidates_gdf(sinks, meta), polygons, static_cast<int>(sinks.size())};
}

} // namespace

void bind_internal_drainage(nb::module_& module)
{
    nb::class_<InternalDrainageResult>(module, "InternalDrainageResult", R"doc(
Closed-depression candidates and, optionally, their GIS-ready footprints.

Candidates are point GeoDataFrames with hierarchy, depth, area, and storage
metadata. Footprints are polygons keyed by ``sink_id``.

Authors
-------
B. Gailleton
)doc")
        .def_prop_ro("candidates", [](const InternalDrainageResult& self) { return self.candidates; })
        .def_prop_ro("polygons", [](const InternalDrainageResult& self) { return self.polygons; })
        .def_prop_ro("n_candidates", [](const InternalDrainageResult& self) { return self.n_candidates; });

    nb::class_<DepressionHierarchyResult>(module, "DepressionHierarchyResult", R"doc(
Complete Barnes depression hierarchy with raster labels and a node table.

Authors
-------
B. Gailleton
)doc")
        .def_prop_ro("labels", &DepressionHierarchyResult::labels)
        .def_prop_ro("flow_direction", &DepressionHierarchyResult::flow_direction)
        .def_prop_ro("nodes", &DepressionHierarchyResult::nodes)
        .def_prop_ro("n_nodes", &DepressionHierarchyResult::n_nodes);

    module.def("depression_hierarchy", &depression_hierarchy,
        nb::arg("dem"), nb::arg("boundary_conditions") = "oooo",
        nb::arg("ocean_from_boundary") = true, nb::arg("nodata_is_ocean") = false,
        nb::arg("ocean_level") = nb::none(), R"doc(
Build the complete nested depression hierarchy for a DEM.

Parameters
----------
dem : Raster
    Projected elevation raster.
boundary_conditions : str, optional
    Open edges seed the ocean when ``ocean_from_boundary=True``.
ocean_from_boundary, nodata_is_ocean : bool, optional
    Define ocean seeds from open edges and/or nodata cells.
ocean_level : float, optional
    Treat valid cells at this elevation as ocean seeds.

Returns
-------
DepressionHierarchyResult
    Leaf labels, D8 spill-routing directions, and leaf/meta-depression nodes.

Authors
-------
B. Gailleton
)doc");

    module.def("find_internal_drainage", &find_internal_drainage,
        nb::arg("dem"), nb::arg("method") = "hierarchy", nb::arg("boundary_conditions") = "oooo",
        nb::arg("min_pixels") = 1, nb::arg("min_depth_metres") = 0.0f,
        nb::arg("min_area_metres2") = 0.0, nb::arg("min_volume_metres3") = 0.0,
        nb::arg("leaves_only") = true, nb::arg("include_polygons") = true, nb::arg("n_threads") = 0,
        R"doc(Detect internally drained topographic depressions in a projected DEM.

Parameters
----------
dem : Raster
    Projected elevation raster.
method : {"hierarchy", "fill_difference"}, optional
    Nested Barnes depression hierarchy, or the legacy filled-minus-raw method.
boundary_conditions : str, optional
    Four-character north/east/south/west boundary condition string.
min_pixels : int, optional
    Minimum depression footprint size in raster cells.
min_depth_metres, min_area_metres2, min_volume_metres3 : float, optional
    Post-detection filters on depression properties.
leaves_only : bool, optional
    For ``method="hierarchy"``, return only innermost depressions.
include_polygons : bool, optional
    Derive footprint polygons as well as candidate points.
n_threads : int, optional
    OpenMP thread count; zero uses the runtime default.

Returns
-------
InternalDrainageResult
    Point candidates and optional polygon footprints, both in the DEM CRS.

Notes
-----
Use the candidate ``pit_row`` and ``pit_col`` columns with
``fill_pits_with_sinks`` and ``FlowInfo(..., forced_outlets=...)`` when
modelling selected depressions as drainage exits.

Authors
-------
B. Gailleton
)doc");
}
} // namespace lsdtt3_python
