#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <limits>
#include <optional>
#include <sstream>
#include <cstring>
#include <stdexcept>
#include <string>
#include <variant>
#include <vector>

#include <nanobind/nanobind.h>
#include <nanobind/stl/optional.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>

#include <cpl_vsi.h>
#include <ogr_spatialref.h>

#include <core/lsdtt_param_parser.hpp>
#include <core/lsdtt_river_line_fetch.hpp>
#include <cpu/lsdtt_burn_corridor_cpu.hpp>
#include <cpu/lsdtt_burn_streams_cpu.hpp>
#include <cpu/lsdtt_raster_flow_preprocessing_cpu.hpp>
#include <cpu/lsdtt_river_line_prepare_cpu.hpp>
#include <cpu/lsdtt_vector_to_junction_network.hpp>
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

// ---- burn_to_corridor (lsdtt-stream-burning stream_burn_mode burn_to_corridor) ----
struct BurnCorridorResult {
    Raster dem;            // the final, conditioned DEM (<write_fname>_burned.tif)
    Raster burned;         // after the burn, before the final conditioning
    Raster burn_depth;     // original minus final (<write_fname>_burn_depth.tif)
    nb::object profiles;   // _burned_profiles.csv
    nb::object channels;   // _corridor_channels
    nb::object exit_routes;
    nb::object cross_path_receivers;
    nb::dict report;
};

nb::object nan_or(double v) { return std::isfinite(v) ? nb::object(nb::float_(v)) : nb::object(nb::float_(std::numeric_limits<double>::quiet_NaN())); }

BurnCorridorResult burn_to_corridor(
    const Raster& dem, nb::object lines, const std::string& vector_name_field,
    double vector_corridor_half_width_metres, int burn_profile_max_pit_pixels,
    double burn_profile_max_cut_metres, const std::string& burn_profile_fit_algorithm_choice,
    double m_over_n, double chi_A0, double fill_min_slope, const std::string& boundary_conditions,
    bool carve_before_fill, const std::string& breach_algorithm_choice, int breach_max_dist_pixels,
    double breach_max_cost, bool breach_minimize_dist, double breach_flat_increment,
    const std::string& breach_edge_policy, double min_elevation, int n_threads)
{
    // The same option values and checks as lsdtt-stream-burning.
    lsdtt::burn_corridor_options co;
    co.corridor_half_width_metres = vector_corridor_half_width_metres;
    co.min_gradient = std::max(1.0e-6, fill_min_slope);
    co.max_pit_pixels = burn_profile_max_pit_pixels;
    co.max_cut_metres = burn_profile_max_cut_metres;
    co.fit = lsdtt::parse_burn_profile_fit(burn_profile_fit_algorithm_choice, "burn_profile_fit_algorithm_choice");
    co.m_over_n = m_over_n;
    co.chi_A0 = chi_A0;
    co.boundary_conditions = boundary_conditions;
    co.breach.max_dist = breach_max_dist_pixels;
    co.breach.max_cost = breach_max_cost > 0.0 ? breach_max_cost : std::numeric_limits<double>::infinity();
    co.breach.minimize_dist = breach_minimize_dist;
    co.breach.flat_increment = breach_flat_increment;
    co.breach.fill_min_slope = static_cast<float>(fill_min_slope);
    co.breach.edge_policy = lsdtt::parse_breach_lc_edge_policy(breach_edge_policy, "breach_edge_policy");
    lsdtt::validate_burn_corridor_options(co);

    lsdtt::lsdtt_raster_cpu input = dem.native();
    long n_removed = 0;
    if (min_elevation > -9998.0) n_removed = static_cast<long>(lsdtt::apply_min_elevation(input, static_cast<float>(min_elevation)));

    // The lines, in the DEM's CRS, through the CLI's clip.
    const std::string& dem_wkt = input.get_meta().crs_wkt;
    nb::object frame = lines;
    if (!frame.attr("crs").is_none() && !dem_wkt.empty()) frame = frame.attr("to_crs")(dem_wkt);
    (void)geodataframe_wkt(frame, "lines");
    if (!vector_name_field.empty() &&
        !nb::cast<bool>(frame.attr("columns").attr("__contains__")(vector_name_field)))
        throw std::invalid_argument("vector_name_field '" + vector_name_field + "' is not a column of lines");
    const std::string path = geodataframe_geojson(frame, "lines");
    lsdtt::clipped_vector_lines cl;
    try { cl = lsdtt::clip_vector_lines_to_dem_with_names(path, input, vector_name_field); }
    catch (...) { VSIUnlink(path.c_str()); throw; }
    VSIUnlink(path.c_str());

    lsdtt::lsdtt_raster_cpu burned;
    const lsdtt::burn_corridor_result cr = lsdtt::burn_vector_network_to_corridor(input, burned, cl.lines, cl.names, co);
    lsdtt::corridor_final_conditioning_options fo;
    fo.carve_before_fill = carve_before_fill;
    fo.breach_algorithm_choice = breach_algorithm_choice;
    fo.breach = co.breach;
    fo.fill_min_slope = static_cast<float>(fill_min_slope);
    fo.n_threads = n_threads;
    lsdtt::lsdtt_raster_cpu final_dem = lsdtt::condition_after_corridor_burn(burned, fo);
    const lsdtt::corridor_burn_checks ck = lsdtt::check_corridor_burn(cr, burned, final_dem);

    nb::object pandas = nb::module_::import_("pandas");
    nb::object geopandas = nb::module_::import_("geopandas");
    nb::object shp = nb::module_::import_("shapely.geometry");
    const std::string& crs = final_dem.get_meta().crs_wkt;
    auto coords = [](const std::vector<lsdtt::lsdtt_xy>& xy) { nb::list c; for (const auto& q : xy) c.append(nb::make_tuple(q.x, q.y)); return c; };

    // profiles: the CSV's columns
    const auto rows = lsdtt::burned_profile_rows(cr, final_dem);
    const auto& names = lsdtt::burned_profile_columns();
    std::vector<nb::list> col(names.size());
    for (const auto& r : rows)
    {
        size_t k = 0;
        col[k++].append(r.path_id); col[k++].append(r.path_rank); col[k++].append(r.reach_id); col[k++].append(r.name);
        col[k++].append(r.point_index); col[k++].append(nan_or(r.distance_m)); col[k++].append(nan_or(r.flow_distance_m));
        col[k++].append(r.x); col[k++].append(r.y); col[k++].append(nan_or(r.z_original)); col[k++].append(nan_or(r.z_cleaned));
        col[k++].append(nan_or(r.z_profile)); col[k++].append(nan_or(r.z_final)); col[k++].append(r.pit_filled);
        col[k++].append(r.pit_id); col[k++].append(r.is_record_low); col[k++].append(r.is_rejected_low);
        col[k++].append(nan_or(r.drainage_area_m2)); col[k++].append(nan_or(r.chi)); col[k++].append(r.joins_path_id);
        col[k++].append(nan_or(r.join_distance_m)); col[k++].append(r.starts_on_path_id);
    }
    nb::dict pd; for (size_t k = 0; k < names.size(); ++k) pd[nb::str(names[k].c_str())] = col[k];
    nb::object profiles = pandas.attr("DataFrame")(pd);

    // channels: _corridor_channels
    nb::dict ch; nb::list g;
    const char* ch_names[] = {"reach_id", "from_node", "to_node", "path_id", "name", "status", "skip_reason", "length_m",
                              "upstream_length_m", "max_separation_m", "max_cut_m", "max_raise_m", "n_record_lows",
                              "n_rejected_lows", "line_index", "piece", "n_pieces"};
    std::vector<nb::list> cc(17);
    for (const auto& rc : cr.reaches)
    {
        const auto& geom = rc.channel.size() >= 2 ? rc.channel
                           : (rc.piece_line.size() >= 2 ? rc.piece_line : cl.lines[static_cast<size_t>(rc.line_index)]);
        if (geom.size() < 2) continue;
        g.append(shp.attr("LineString")(coords(geom)));
        size_t k = 0;
        cc[k++].append(rc.reach_id); cc[k++].append(rc.from_node); cc[k++].append(rc.to_node); cc[k++].append(rc.path_id);
        cc[k++].append(rc.name); cc[k++].append(rc.status); cc[k++].append(rc.skip_reason); cc[k++].append(rc.length_m);
        cc[k++].append(rc.upstream_length_m); cc[k++].append(rc.max_separation_m); cc[k++].append(rc.max_cut_m);
        cc[k++].append(rc.max_raise_m); cc[k++].append(rc.n_record_lows); cc[k++].append(rc.n_rejected_lows);
        cc[k++].append(rc.line_index); cc[k++].append(rc.piece); cc[k++].append(rc.n_pieces);
    }
    for (size_t k = 0; k < 17; ++k) ch[ch_names[k]] = cc[k];
    nb::object channels = geopandas.attr("GeoDataFrame")(ch, nb::arg("geometry") = g, nb::arg("crs") = crs);

    // exit routes
    nb::dict xr; nb::list xg;
    std::vector<nb::list> xc(7);
    for (const auto& r : cr.exit_routes)
    {
        if (r.xy.size() < 2) continue;
        xg.append(shp.attr("LineString")(coords(r.xy)));
        xc[0].append(r.path_id); xc[1].append(r.exit_type); xc[2].append(r.to_path_id); xc[3].append(r.length_m);
        xc[4].append(static_cast<long>(r.cells.size())); xc[5].append(r.n_cells_lowered); xc[6].append(r.max_cut_m);
    }
    const char* x_names[] = {"path_id", "exit_type", "to_path_id", "length_m", "n_cells", "n_cells_lowered", "max_cut_m"};
    for (size_t k = 0; k < 7; ++k) xr[x_names[k]] = xc[k];
    nb::object exits = geopandas.attr("GeoDataFrame")(xr, nb::arg("geometry") = xg, nb::arg("crs") = crs);

    // cross-path receivers
    nb::dict pr; nb::list pg;
    std::vector<nb::list> pc(9);
    for (const auto& c : cr.cross_path_receivers)
    {
        pg.append(shp.attr("Point")(c.xy.x, c.xy.y));
        pc[0].append(c.path_id); pc[1].append(c.index); pc[2].append(c.receiver_path_id); pc[3].append(c.reason);
        pc[4].append(c.z_m); pc[5].append(c.z_own_next_m); pc[6].append(c.z_receiver_m);
        pc[7].append(c.receiver_xy.x); pc[8].append(c.receiver_xy.y);
    }
    const char* p_names[] = {"path_id", "index", "receiver_path_id", "reason", "z_m", "z_own_next_m", "z_receiver_m",
                             "receiver_x", "receiver_y"};
    for (size_t k = 0; k < 9; ++k) pr[p_names[k]] = pc[k];
    nb::object cross = geopandas.attr("GeoDataFrame")(pr, nb::arg("geometry") = pg, nb::arg("crs") = crs);

    nb::dict rep;
    rep["min_elevation_cells_removed"] = n_removed;
    rep["n_lines"] = static_cast<long>(cl.lines.size());
    rep["n_nodes"] = cr.n_nodes; rep["n_reaches"] = static_cast<long>(cr.reaches.size());
    rep["n_fitted"] = cr.n_fitted; rep["n_skipped"] = cr.n_skipped; rep["n_paths"] = static_cast<long>(cr.paths.size());
    rep["n_lines_split_at_nodata"] = cr.n_lines_split_at_nodata; rep["n_pieces"] = cr.n_pieces;
    rep["n_pits_cleaned"] = cr.n_pits_cleaned; rep["n_pit_cells"] = cr.n_pit_cells;
    rep["n_record_lows"] = cr.n_record_lows; rep["n_rejected_lows"] = cr.n_rejected_lows;
    rep["n_pin_conflicts"] = cr.n_pin_conflicts; rep["n_cut_cap_exceeded"] = cr.n_cut_cap_exceeded;
    rep["n_line_cells"] = cr.n_line_cells; rep["n_line_cells_cut"] = cr.n_line_cells_cut;
    rep["n_line_cells_raised"] = cr.n_line_cells_raised;
    rep["n_ends_drained"] = cr.n_ends_drained; rep["n_cells_daylighted"] = cr.n_cells_daylighted;
    rep["n_offline_cells_carved"] = cr.n_offline_cells_carved; rep["n_offline_cells_raised"] = cr.n_offline_cells_raised;
    rep["n_bank_cells"] = cr.n_bank_cells; rep["n_cross_path_raised"] = cr.n_cross_path_raised;
    rep["n_cross_path_unresolved"] = cr.n_cross_path_unresolved;
    rep["final_line_cells_changed"] = ck.n_line_cells_changed; rep["final_other_cells_changed"] = ck.n_other_cells_changed;
    rep["n_paths_rising_after_burn"] = ck.n_rising_after_burn; rep["n_paths_rising_final"] = ck.n_rising_final;
    nb::list notes; for (const auto& n : cr.notes) notes.append(n);
    rep["notes"] = notes;

    lsdtt::lsdtt_raster_cpu depth = lsdtt::corridor_burn_depth(input, final_dem);
    return {Raster(std::move(final_dem)), Raster(std::move(burned)), Raster(std::move(depth)), profiles, channels,
            exits, cross, rep};
}

// ---- fetch_river_lines (lsdtt-fetch-river-lines) ----
struct RiverLinesResult {
    nb::object lines, nodes, topology;
    std::string attribution;
    nb::dict report;
};

std::string degrees_text(std::optional<double> v)
{
    if (!v) return "";
    std::ostringstream o; o << std::setprecision(17) << *v; return o.str();
}

RiverLinesResult fetch_river_lines(
    const std::string& river_line_source, std::optional<double> west, std::optional<double> south,
    std::optional<double> east, std::optional<double> north, const std::string& bounds_fname,
    const std::string& bounds_layer, int output_epsg, double clip_buffer_metres,
    const std::string& river_lines_fname, const std::string& river_lines_layer,
    const std::string& feature_class_field_name, const std::string& direction_field_name,
    const std::string& direction_field_mapping, bool digitised_downstream,
    const std::string& start_node_field_name, const std::string& end_node_field_name,
    const std::string& culvert_field_name, const std::string& river_feature_classes,
    double node_snap_tolerance_metres, double min_line_length_metres,
    const std::string& direction_repair_algorithm_choice, const std::string& gap_bridge_algorithm_choice,
    double max_gap_length_metres, double max_nonroad_gap_metres, double max_facing_angle_degrees,
    double max_heading_misalignment_degrees, double road_max_facing_angle_degrees,
    double tangent_fit_length_metres, double road_intersection_tolerance_metres, bool require_mutual_best,
    const std::string& road_vector_fname, const std::string& road_layer, const std::string& road_name_field,
    const std::string& river_cache_dir, int http_timeout_seconds, int http_connect_timeout_seconds,
    int http_max_retry)
{
    // The CLI's parameter names and the CLI's option builder.
    lsdtt::lsdtt_river_line_parameters p;
    p.s = {{"river_line_source", river_line_source}, {"west", degrees_text(west)}, {"south", degrees_text(south)},
           {"east", degrees_text(east)}, {"north", degrees_text(north)}, {"bounds_fname", bounds_fname},
           {"bounds_layer", bounds_layer}, {"river_lines_fname", river_lines_fname},
           {"river_lines_layer", river_lines_layer}, {"feature_class_field_name", feature_class_field_name},
           {"direction_field_name", direction_field_name}, {"direction_field_mapping", direction_field_mapping},
           {"start_node_field_name", start_node_field_name}, {"end_node_field_name", end_node_field_name},
           {"culvert_field_name", culvert_field_name}, {"river_feature_classes", river_feature_classes},
           {"direction_repair_algorithm_choice", direction_repair_algorithm_choice},
           {"gap_bridge_algorithm_choice", gap_bridge_algorithm_choice}, {"road_vector_fname", road_vector_fname},
           {"road_layer", road_layer}, {"road_name_field", road_name_field}, {"river_cache_dir", river_cache_dir}};
    p.f = {{"clip_buffer_metres", clip_buffer_metres}, {"node_snap_tolerance_metres", node_snap_tolerance_metres},
           {"min_line_length_metres", min_line_length_metres}, {"max_gap_length_metres", max_gap_length_metres},
           {"max_nonroad_gap_metres", max_nonroad_gap_metres}, {"max_facing_angle_degrees", max_facing_angle_degrees},
           {"max_heading_misalignment_degrees", max_heading_misalignment_degrees},
           {"road_max_facing_angle_degrees", road_max_facing_angle_degrees},
           {"tangent_fit_length_metres", tangent_fit_length_metres},
           {"road_intersection_tolerance_metres", road_intersection_tolerance_metres}};
    p.i = {{"output_epsg", output_epsg}, {"http_timeout_seconds", http_timeout_seconds},
           {"http_connect_timeout_seconds", http_connect_timeout_seconds}, {"http_max_retry", http_max_retry}};
    p.b = {{"digitised_downstream", digitised_downstream}, {"require_mutual_best", require_mutual_best}};
    std::vector<std::string> warnings;
    std::string area_note;
    lsdtt::lsdtt_river_line_fetch_options o =
        lsdtt::build_river_line_fetch_options(p, [](const std::string& f) { return f; }, warnings, &area_note);
    o.http.user_agent = "lsdtt3 python fetch_river_lines/" + lsdtt::lsdtt3_package_version();

    // The CLI's writer, into memory, read back: the same fields as the files.
    static std::uint64_t counter = 0;
    const std::string stem = "/vsimem/lsdtt3-river-lines-" + std::to_string(++counter);
    o.out_lines_path = stem + "_lines.fgb";
    o.out_nodes_path = stem + "_nodes.fgb";
    o.out_format = lsdtt::vector_output_format::flatgeobuf;
    o.spatial_index = false;
    lsdtt::lsdtt_river_line_fetch_result r;
    nb::object lines_gdf, nodes_gdf;
    try
    {
        r = lsdtt::fetch_river_lines(o);
        lines_gdf = layer_geodataframe(read_vector_layer(o.out_lines_path), r.crs_wkt);
        lsdtt::lsdtt_vector_layer nl = read_vector_layer(o.out_nodes_path);
        nb::object shp = nb::module_::import_("shapely.geometry");
        nb::dict d; nb::list g;
        std::vector<nb::list> cols(nl.schema.size());
        for (const auto& rec : nl.records)
        {
            if (const auto* pt = std::get_if<lsdtt::lsdtt_point_geom>(&rec.geometry)) g.append(shp.attr("Point")(pt->x, pt->y));
            else g.append(nb::none());
            for (size_t k = 0; k < cols.size(); ++k)
            {
                const auto& v = rec.field_values[k];
                if (std::holds_alternative<std::monostate>(v)) cols[k].append(nb::none());
                else if (const auto* x = std::get_if<std::int64_t>(&v)) cols[k].append(*x);
                else if (const auto* x = std::get_if<double>(&v)) cols[k].append(*x);
                else cols[k].append(std::get<std::string>(v));
            }
        }
        for (size_t k = 0; k < cols.size(); ++k) d[nb::str(nl.schema.field_names[k].c_str())] = cols[k];
        nodes_gdf = nb::module_::import_("geopandas").attr("GeoDataFrame")(d, nb::arg("geometry") = g, nb::arg("crs") = r.crs_wkt);
    }
    catch (...) { VSIUnlink(o.out_lines_path.c_str()); VSIUnlink(o.out_nodes_path.c_str()); throw; }
    VSIUnlink(o.out_lines_path.c_str()); VSIUnlink(o.out_nodes_path.c_str());

    // topology: the CSV's columns
    const char* t_names[] = {"component_id", "n_lines", "length_m", "n_tips", "n_outlets", "n_clip_edges", "n_loops",
                             "n_bifurcations", "n_conflicts", "outlet_node", "n_against_network"};
    std::vector<nb::list> tc(11);
    for (const auto& c : r.network.components)
    {
        tc[0].append(c.component_id); tc[1].append(c.n_lines); tc[2].append(c.length_m); tc[3].append(c.n_tips);
        tc[4].append(c.n_outlets); tc[5].append(c.n_clip_edges); tc[6].append(c.n_loops); tc[7].append(c.n_bifurcations);
        tc[8].append(c.n_conflicts); tc[9].append(c.outlet_node); tc[10].append(c.n_against_network);
    }
    nb::dict td; for (size_t k = 0; k < 11; ++k) td[t_names[k]] = tc[k];
    nb::object topology = nb::module_::import_("pandas").attr("DataFrame")(td);

    nb::dict rep;
    rep["crs_description"] = r.crs_description; rep["epsg"] = r.epsg; rep["read_path"] = r.read_path;
    rep["n_features_read"] = r.n_features_read; rep["n_non_line_features"] = r.n_non_line_features;
    rep["clip_buffer_metres"] = r.clip_buffer_metres; rep["node_snap_tolerance_metres"] = r.node_snap_tolerance_metres;
    rep["area"] = nb::make_tuple(r.area_x_min, r.area_y_min, r.area_x_max, r.area_y_max);
    nb::dict classes; for (const auto& [k, v] : r.classes) classes[k.c_str()] = nb::make_tuple(v.n_read, v.n_kept);
    rep["classes"] = classes;
    nb::list notes; if (!area_note.empty()) notes.append(area_note);
    for (const auto& n : r.notes) notes.append(n);
    for (const auto& n : r.network.notes) notes.append(n);
    rep["notes"] = notes;
    nb::list warns; for (const auto& w : warnings) warns.append(w); for (const auto& w : r.warnings) warns.append(w);
    rep["warnings"] = warns;
    rep["n_bridges_gap"] = r.gap_bridges.n_bridges_gap; rep["n_bridges_road"] = r.gap_bridges.n_bridges_road;
    return {lines_gdf, nodes_gdf, topology, r.attribution, rep};
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

    nb::class_<BurnCorridorResult>(module, "BurnCorridorResult",
        "Result of :func:`burn_to_corridor` (lsdtt-stream-burning burn_to_corridor): the rasters and the tables the CLI writes.")
        .def_ro("dem", &BurnCorridorResult::dem, "The final DEM: burned, then breached / filled (the CLI's _burned.tif).")
        .def_ro("burned", &BurnCorridorResult::burned, "The DEM after the burn, before the final conditioning.")
        .def_ro("burn_depth", &BurnCorridorResult::burn_depth,
            "Original minus final (m): positive where lowered, negative where raised, nodata where unchanged (_burn_depth.tif).")
        .def_ro("profiles", &BurnCorridorResult::profiles,
            "pandas.DataFrame, one row per path cell, the columns of _burned_profiles.csv (lsdviztools3 lsdviz3-burned-profiles reads them).")
        .def_ro("channels", &BurnCorridorResult::channels,
            "GeoDataFrame, one line per reach: the channel used, or the vector line if skipped (_corridor_channels).")
        .def_ro("exit_routes", &BurnCorridorResult::exit_routes,
            "GeoDataFrame, the routes lowered to drain channel ends inside the DEM (_exit_routes).")
        .def_ro("cross_path_receivers", &BurnCorridorResult::cross_path_receivers,
            "GeoDataFrame of points: line cells whose flow leaves the line for another path's cell (_cross_path_receivers).")
        .def_ro("report", &BurnCorridorResult::report,
            "dict of counts: reaches fitted / skipped, paths, pits, minima, ends drained, banks, cross-path cells, the final-pass changes and the downhill checks (all 0 when the burn drains everywhere), notes.");
    nb::class_<RiverLinesResult>(module, "RiverLinesResult",
        "Result of :func:`fetch_river_lines` (lsdtt-fetch-river-lines): the prepared lines, nodes and topology.")
        .def_ro("lines", &RiverLinesResult::lines, "GeoDataFrame of the prepared lines (the _river_lines file's fields), oriented downstream.")
        .def_ro("nodes", &RiverLinesResult::nodes, "GeoDataFrame of the nodes (_river_nodes).")
        .def_ro("topology", &RiverLinesResult::topology, "pandas.DataFrame, one row per connected network (_river_topology.csv).")
        .def_ro("attribution", &RiverLinesResult::attribution, "Licence and citation of the source (_ATTRIBUTION.txt).")
        .def_ro("report", &RiverLinesResult::report, "dict: output CRS, area, counts by feature class, notes and warnings.");

    module.def("burn_to_corridor", &burn_to_corridor, nb::arg("dem"), nb::arg("lines"), nb::kw_only(),
        nb::arg("vector_name_field") = "", nb::arg("vector_corridor_half_width_metres") = 0.0,
        nb::arg("burn_profile_max_pit_pixels") = 16, nb::arg("burn_profile_max_cut_metres") = 2000.0,
        nb::arg("burn_profile_fit_algorithm_choice") = "linear", nb::arg("m_over_n") = 0.45, nb::arg("chi_A0") = 1.0,
        nb::arg("fill_min_slope") = 0.0001, nb::arg("boundary_conditions") = "oooo",
        nb::arg("carve_before_fill") = false, nb::arg("breach_algorithm_choice") = "lindsay2016",
        nb::arg("breach_max_dist_pixels") = 1000, nb::arg("breach_max_cost") = 0.0,
        nb::arg("breach_minimize_dist") = false, nb::arg("breach_flat_increment") = 0.0,
        nb::arg("breach_edge_policy") = "outlet", nb::arg("min_elevation") = -9999.0, nb::arg("n_threads") = 0,
        R"doc(
Burn a vector river network into a DEM, taking the vector as the truth.

lsdtt-stream-burning with ``stream_burn_mode: burn_to_corridor``, through the
same C++ routines and with its parameter names. Each line's channel is found
within ``vector_corridor_half_width_metres`` of it, and burned to a profile
joining the minima of the input DEM under it (small pits filled first;
minima joined linearly, or in chi). The rest of the DEM is then carved and
filled to drain to the channels, and the whole DEM is conditioned once more
(breach if ``carve_before_fill``, then fill), so flow routing on the result
follows the network.

Parameters
----------
dem : Raster
    The input DEM (not pit-filled).
lines : geopandas.GeoDataFrame
    River lines, oriented downstream (e.g. the ``lines`` of
    :func:`fetch_river_lines`), with a CRS; reprojected to the DEM's.
vector_name_field : str, default=""
    Column holding the river name, copied to the outputs; empty for none.
vector_corridor_half_width_metres : float, default=0.0
    X: how far (m) the channel may move from the line. 0 burns the line's
    own cells.
burn_profile_max_pit_pixels : int, default=16
    Pits of fewer cells are filled before the minima are taken.
burn_profile_max_cut_metres : float, default=2000.0
    A minimum that would cut deeper than this (m) is not used.
burn_profile_fit_algorithm_choice : str, default="linear"
    How minima are joined: ``linear`` (in distance) or ``chi``.
m_over_n : float, default=0.45
    Concavity for the chi fit.
chi_A0 : float, default=1.0
    Reference area (m^2) for the chi fit.
fill_min_slope : float, default=0.0001
    Minimum gradient (m/m) of the burned profile and of filled flats.
boundary_conditions : str, default="oooo"
    NESW edge codes for the pit cleaning.
carve_before_fill : bool, default=False
    Breach before the final fill.
breach_algorithm_choice : str, default="lindsay2016"
    ``lindsay2016`` or ``least_cost``, for the final breach.
breach_max_dist_pixels : int, default=1000
    Least-cost breach: longest breach path (cells).
breach_max_cost : float, default=0.0
    Least-cost breach: largest cost; 0 is unbounded.
breach_minimize_dist : bool, default=False
    Least-cost breach: weight the cost by distance.
breach_flat_increment : float, default=0.0
    Least-cost breach: drop per cell along a breach; 0 is automatic.
breach_edge_policy : str, default="outlet"
    Least-cost breach: ``outlet`` or ``whitebox``.
min_elevation : float, default=-9999.0
    Cells at or below this (m) become nodata first; -9999 disables. Leave it
    off for DTMs with land below 0 m.
n_threads : int, default=0
    Threads for the final conditioning; 0 is the default.

Returns
-------
BurnCorridorResult
    ``dem`` (final), ``burned``, ``burn_depth``, ``profiles``,
    ``channels``, ``exit_routes``, ``cross_path_receivers`` and ``report``.
)doc");

    module.def("fetch_river_lines", &fetch_river_lines, nb::kw_only(), nb::arg("river_line_source"),
        nb::arg("west") = nb::none(), nb::arg("south") = nb::none(), nb::arg("east") = nb::none(),
        nb::arg("north") = nb::none(), nb::arg("bounds_fname") = "", nb::arg("bounds_layer") = "",
        nb::arg("output_epsg") = 0, nb::arg("clip_buffer_metres") = -1.0, nb::arg("river_lines_fname") = "",
        nb::arg("river_lines_layer") = "", nb::arg("feature_class_field_name") = "",
        nb::arg("direction_field_name") = "", nb::arg("direction_field_mapping") = "",
        nb::arg("digitised_downstream") = false, nb::arg("start_node_field_name") = "",
        nb::arg("end_node_field_name") = "", nb::arg("culvert_field_name") = "",
        nb::arg("river_feature_classes") = "river,stream,lake_connector", nb::arg("node_snap_tolerance_metres") = -1.0,
        nb::arg("min_line_length_metres") = 0.0, nb::arg("direction_repair_algorithm_choice") = "network_outlet",
        nb::arg("gap_bridge_algorithm_choice") = "none", nb::arg("max_gap_length_metres") = 25.0,
        nb::arg("max_nonroad_gap_metres") = 30.0, nb::arg("max_facing_angle_degrees") = 45.0,
        nb::arg("max_heading_misalignment_degrees") = 50.0, nb::arg("road_max_facing_angle_degrees") = 60.0,
        nb::arg("tangent_fit_length_metres") = 15.0, nb::arg("road_intersection_tolerance_metres") = 10.0,
        nb::arg("require_mutual_best") = true, nb::arg("road_vector_fname") = "", nb::arg("road_layer") = "",
        nb::arg("road_name_field") = "roadName", nb::arg("river_cache_dir") = "",
        nb::arg("http_timeout_seconds") = 60, nb::arg("http_connect_timeout_seconds") = 15,
        nb::arg("http_max_retry") = 3,
        R"doc(
Fetch or read river centrelines for an area and prepare them for burning.

lsdtt-fetch-river-lines, through the same C++ routines, with its parameter
names and the same checks: the lines are clipped to the area (plus a buffer),
noded, oriented downstream (a direction vote per network) and given a
topology. Nothing is written to disk except the download cache.

Parameters
----------
river_line_source : str, keyword-only
    ``local_file`` (any OGR line file) or ``gb_os_open_rivers`` (OS Open
    Rivers, downloaded keylessly and cached).
west, south, east, north : float or None, default=None
    The area in WGS84 degrees; or give ``bounds_fname``.
bounds_fname : str, default=""
    A raster or vector file whose extent is the area (a raster also sets the
    output CRS and the clip box).
bounds_layer : str, default=""
    Layer of a vector ``bounds_fname``; empty for the first.
output_epsg : int, default=0
    Output CRS; 0 is the bounds raster's CRS, else the UTM zone of the area.
clip_buffer_metres : float, default=-1.0
    Margin (m) around the area; negative is 5% of the shorter side, at least
    200 m.
river_lines_fname : str, default=""
    ``local_file``: the line file.
river_lines_layer : str, default=""
    ``local_file``: its layer; empty for the first.
feature_class_field_name : str, default=""
    ``local_file``: field giving the feature class.
direction_field_name : str, default=""
    ``local_file``: field giving the flow direction.
direction_field_mapping : str, default=""
    ``local_file``: values of that field and their meaning, as
    ``value=1`` (digitised downstream) / ``value=-1`` (upstream) pairs,
    comma-separated, e.g. ``"in_direction=1,in_opposite_direction=-1"``;
    spaces in a value are written as underscores.
digitised_downstream : bool, default=False
    ``local_file``: the lines are digitised downstream.
start_node_field_name : str, default=""
    ``local_file``: field of the start node id (with ``end_node_field_name``).
end_node_field_name : str, default=""
    ``local_file``: field of the end node id.
culvert_field_name : str, default=""
    ``local_file``: field flagging culverts.
river_feature_classes : str, default="river,stream,lake_connector"
    Feature classes kept, comma-separated (river, stream, lake_connector,
    canal, drain, tidal, other).
node_snap_tolerance_metres : float, default=-1.0
    Ends within this distance (m) are one node; negative is 1 m.
min_line_length_metres : float, default=0.0
    Shorter dangling lines are dropped.
direction_repair_algorithm_choice : str, default="network_outlet"
    ``network_outlet`` flips lines against their network's flow; ``none``
    only flags them.
gap_bridge_algorithm_choice : str, default="none"
    ``facing_tips`` joins dangling ends that face each other across a gap.
max_gap_length_metres : float, default=25.0
    Gap bridging: longest road-associated bridge (m).
max_nonroad_gap_metres : float, default=30.0
    Gap bridging: longest bridge not crossing a road (m).
max_facing_angle_degrees : float, default=45.0
    Gap bridging: largest angle between a tip and the bridge.
max_heading_misalignment_degrees : float, default=50.0
    Gap bridging: largest misalignment of the two tips' headings.
road_max_facing_angle_degrees : float, default=60.0
    Gap bridging: facing angle allowed across a road.
tangent_fit_length_metres : float, default=15.0
    Gap bridging: length (m) of line used for a tip's heading.
road_intersection_tolerance_metres : float, default=10.0
    Gap bridging: distance (m) within which a bridge crosses a road.
require_mutual_best : bool, default=True
    Gap bridging: only join tips that are each other's best candidate.
road_vector_fname : str, default=""
    Gap bridging: road centrelines (any OGR line file).
road_layer : str, default=""
    Gap bridging: the roads layer; empty for the first.
road_name_field : str, default="roadName"
    Gap bridging: the roads' name field.
river_cache_dir : str, default=""
    Download cache; empty is the user cache directory.
http_timeout_seconds : int, default=60
    Download timeout (s).
http_connect_timeout_seconds : int, default=15
    Connection timeout (s).
http_max_retry : int, default=3
    Download retries.

Returns
-------
RiverLinesResult
    ``lines``, ``nodes``, ``topology``, ``attribution`` and ``report``.
)doc");

    module.def("water_centrelines", &water_centrelines,
        nb::arg("areas"), nb::arg("tributaries"), nb::arg("epsg") = 0, nb::arg("cell_size_metres") = 1.0,
        nb::arg("padding_metres") = 5.0, nb::arg("tributary_snap_tolerance_metres") = 2.0,
        nb::arg("centreline_snap_max_metres") = 100.0, nb::arg("simplify_tolerance_metres") = 1.0,
        nb::arg("min_component_length_metres") = 10.0,
        "Extract central flow paths from water polygons and connect tributary lines.\n\nParameters\n----------\nareas, tributaries : geopandas.GeoDataFrame\n    Projected water polygons and linework in the same CRS.\nepsg : int, default=0\n    Output CRS identifier; zero infers it from areas.\ncell_size_metres : float, default=1\n    Internal skeletonisation-grid cell size (m).\npadding_metres : float, default=5\n    Margin (m) added around the polygons on the internal grid.\ntributary_snap_tolerance_metres : float, default=2\n    Distance (m) within which tributary ends snap to a water polygon.\ncentreline_snap_max_metres : float, default=100\n    Largest distance (m) a tributary is extended to reach the centreline.\nsimplify_tolerance_metres : float, default=1\n    Line simplification tolerance (m).\nmin_component_length_metres : float, default=10\n    Centreline pieces shorter than this (m) are dropped.\n\nReturns\n-------\nWaterCentrelineResult\n    Connected linework and count diagnostics.\n\nAuthors\n-------\nB. Gailleton");
    module.def("repair_water_lines", &repair_water_lines,
        nb::arg("water_lines"), nb::arg("roads"), nb::arg("snap_tolerance_metres") = 0.5,
        nb::arg("max_gap_length_metres") = 25.0, nb::arg("max_nonroad_gap_metres") = 30.0,
        nb::arg("road_intersection_tolerance_metres") = 1.0, nb::arg("require_mutual_best") = true,
        "Repair short broken water-line gaps, preferentially across roads.\n\nParameters\n----------\nwater_lines, roads : geopandas.GeoDataFrame\n    Projected line layers in the same CRS.\nsnap_tolerance_metres : float, default=0.5\n    Distance (m) within which line ends are treated as connected.\nmax_gap_length_metres : float, default=25\n    Maximum road-associated bridge length (m).\nmax_nonroad_gap_metres : float, default=30\n    Maximum length (m) of a bridge that does not cross a road.\nroad_intersection_tolerance_metres : float, default=1\n    Distance (m) within which a bridge counts as crossing a road.\nrequire_mutual_best : bool, default=True\n    Only bridge two tips that are each other's best candidate.\n\nReturns\n-------\nWaterLineRepairResult\n    Original linework plus accepted bridge segments and diagnostics.\n\nAuthors\n-------\nB. Gailleton");
    module.def("burn_streams", &burn_streams, nb::arg("dem"), nb::arg("lines"), nb::arg("n_threads") = 0,
        "Condition a copied DEM along mapped stream lines with monotone profiles.\n\nParameters\n----------\ndem : Raster\n    Elevation raster.\nlines : geopandas.GeoDataFrame\n    Stream LineStrings in the DEM CRS.\nn_threads : int, default=0\n    OpenMP thread count; zero uses the runtime default.\n\nReturns\n-------\nStreamBurnResult\n    Conditioned raster and profile diagnostics.\n\nNotes\n-----\nOnly cells whose monotone stream elevation is lower are changed.\n\nAuthors\n-------\nB. Gailleton");
}

} // namespace lsdtt3_python
