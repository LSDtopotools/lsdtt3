#include <cmath>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>

#include <cpu/lsdtt_tidal_channel_cpu.hpp>
#include <cpu/lsdtt_tidal_skeleton_metrics_cpu.hpp>

#include "binding_common.hpp"

namespace lsdtt3_python {
namespace {

nb::object optional_raster(const lsdtt::lsdtt_raster_cpu& raster)
{
    return raster.empty() ? nb::none() : nb::cast(Raster(raster));
}

void require_same_grid(const Raster& first, const Raster& second, const char* name)
{
    const auto& a = first.native().get_meta();
    const auto& b = second.native().get_meta();
    if (a.n_rows != b.n_rows || a.n_cols != b.n_cols || a.x_min != b.x_min ||
        a.y_min != b.y_min || a.cell_size != b.cell_size || a.crs_wkt != b.crs_wkt)
        throw std::invalid_argument(std::string(name) + " must share the DEM grid");
}

lsdtt::tidal_outlet_mode outlet_mode(const std::string& value)
{
    if (value == "edge_nodes") return lsdtt::tidal_outlet_mode::edge_nodes;
    if (value == "endpoint_nodes") return lsdtt::tidal_outlet_mode::endpoint_nodes;
    if (value == "outlet_mask") return lsdtt::tidal_outlet_mode::outlet_mask;
    if (value == "lowest_endpoint") return lsdtt::tidal_outlet_mode::lowest_endpoint;
    throw std::invalid_argument("outlet_mode must be 'edge_nodes', 'endpoint_nodes', 'outlet_mask', or 'lowest_endpoint'");
}

class TidalChannelResult {
public:
    explicit TidalChannelResult(lsdtt::tidal_channel_result result) : result_(std::move(result)) {}
    Raster channel_mask() const { return Raster(result_.channel_mask); }
    nb::object upl() const { return optional_raster(result_.upl); }
    nb::object skeleton() const { return optional_raster(result_.channel_skeleton); }
    nb::object connected_skeleton() const { return optional_raster(result_.connected_skeleton); }
    nb::object skeleton_bridge_mask() const { return optional_raster(result_.skeleton_bridge_mask); }
    nb::list median_rasters() const { nb::list out; for (const auto& r : result_.median_rasters) out.append(Raster(r)); return out; }
    nb::list residual_rasters() const { nb::list out; for (const auto& r : result_.residual_rasters) out.append(Raster(r)); return out; }
    nb::list threshold_masks() const { nb::list out; for (const auto& r : result_.threshold_masks) out.append(Raster(r)); return out; }
    nb::object scales_metres() const { return nb::cast(result_.scales_metres); }
    nb::object thresholds_metres() const { return nb::cast(result_.thresholds_metres); }
private:
    lsdtt::tidal_channel_result result_;
};

TidalChannelResult find_tidal_channels(const Raster& dem, bool use_flood_mask,
                                       double flood_elevation_metres,
                                       std::vector<double> median_window_radii_metres,
                                       std::vector<double> residual_thresholds_metres,
                                       const std::string& scale_mode,
                                       bool apply_binary_closing, double binary_closing_radius_metres,
                                       bool fill_holes, int min_channel_component_pixels,
                                       bool extract_skeleton, int skeleton_spur_max_pixels,
                                       int skeleton_terminal_loop_max_pixels, int skeleton_min_component_pixels,
                                       bool connect_skeleton_segments, double connect_max_gap_pixels,
                                       double connect_max_facing_angle_degrees,
                                       double connect_max_heading_misalignment_degrees,
                                       double connect_max_candidate_score,
                                       double connect_route_corridor_half_width_pixels,
                                       bool compute_upl, int n_threads)
{
    lsdtt::tidal_channel_request request;
    request.use_flood_mask = use_flood_mask;
    request.flood_elevation_metres = flood_elevation_metres;
    request.median_window_radii_metres = std::move(median_window_radii_metres);
    request.residual_thresholds_metres = std::move(residual_thresholds_metres);
    if (scale_mode == "any") request.scale_mode = lsdtt::tidal_channel_scale_mode::any_scale_pass;
    else if (scale_mode == "all") request.scale_mode = lsdtt::tidal_channel_scale_mode::all_scales_pass;
    else throw std::invalid_argument("scale_mode must be 'any' or 'all'");
    request.apply_binary_closing = apply_binary_closing;
    request.binary_closing_radius_metres = binary_closing_radius_metres;
    request.fill_holes = fill_holes;
    request.min_channel_component_pixels = min_channel_component_pixels;
    request.extract_skeleton = extract_skeleton;
    request.skeleton_spur_max_pixels = skeleton_spur_max_pixels;
    request.skeleton_terminal_loop_max_pixels = skeleton_terminal_loop_max_pixels;
    request.skeleton_min_component_pixels = skeleton_min_component_pixels;
    request.connect_skeleton_segments = connect_skeleton_segments;
    request.connect_max_gap_pixels = connect_max_gap_pixels;
    request.connect_max_facing_angle_degrees = connect_max_facing_angle_degrees;
    request.connect_max_heading_misalignment_degrees = connect_max_heading_misalignment_degrees;
    request.connect_max_candidate_score = connect_max_candidate_score;
    request.connect_route_corridor_half_width_pixels = connect_route_corridor_half_width_pixels;
    request.compute_upl = compute_upl;
    request.n_threads = n_threads;
    return TidalChannelResult(lsdtt::find_tidal_channels(dem.native(), request));
}

nb::object station_geodataframe(const lsdtt::lsdtt_tidal_skeleton_metrics_result& result,
                                const std::vector<std::string>& strip_names, const std::string& crs_wkt)
{
    nb::dict data;
    nb::list id, section, row, col, x, y, elevation, along, downstream, strahler, hack, width, depth, area, local_area, upstream_area;
    std::vector<nb::list> local_values(strip_names.size()), upstream_values(strip_names.size()); nb::list geometry;
    for (const auto& station : result.stations) {
        id.append(station.station_id); section.append(station.section_id); row.append(station.row); col.append(station.col);
        x.append(station.x); y.append(station.y); elevation.append(station.elevation);
        along.append(station.along_section_distance_metres); downstream.append(station.downstream_distance_metres);
        strahler.append(station.strahler_order); hack.append(station.hack_order);
        width.append(station.cross_section.valid ? nb::cast(station.cross_section.width_metres) : nb::none());
        depth.append(station.cross_section.valid ? nb::cast(station.cross_section.depth_metres) : nb::none());
        area.append(station.cross_section.valid ? nb::cast(station.cross_section.area_m2) : nb::none());
        local_area.append(station.local_strip_area_m2); upstream_area.append(station.upstream_strip_area_m2);
        geometry.append(nb::module_::import_("shapely.geometry").attr("Point")(station.x, station.y));
        for (std::size_t i = 0; i < strip_names.size(); ++i) {
            local_values[i].append(i < station.local_strip_mean_values.size() ? nb::cast(station.local_strip_mean_values[i]) : nb::none());
            upstream_values[i].append(i < station.upstream_strip_mean_values.size() ? nb::cast(station.upstream_strip_mean_values[i]) : nb::none());
        }
    }
    data["station_id"] = id; data["section_id"] = section; data["row"] = row; data["col"] = col;
    data["x"] = x; data["y"] = y; data["elevation"] = elevation; data["along_section_distance_metres"] = along;
    data["downstream_distance_metres"] = downstream; data["strahler_order"] = strahler; data["hack_order"] = hack;
    data["width_metres"] = width; data["depth_metres"] = depth; data["area_m2"] = area;
    data["local_strip_area_m2"] = local_area; data["upstream_strip_area_m2"] = upstream_area;
    for (std::size_t i = 0; i < strip_names.size(); ++i) {
        data[nb::str(("local_" + strip_names[i] + "_mean").c_str())] = local_values[i];
        data[nb::str(("upstream_" + strip_names[i] + "_mean").c_str())] = upstream_values[i];
    }
    return nb::module_::import_("geopandas").attr("GeoDataFrame")(data, nb::arg("geometry") = geometry, nb::arg("crs") = crs_wkt);
}

nb::object node_geodataframe(const lsdtt::lsdtt_tidal_skeleton_metrics_result& result, const std::string& crs_wkt)
{
    nb::dict data; nb::list id, row, col, elevation, degree, endpoint, junction, outlet, distance, geometry;
    nb::object point = nb::module_::import_("shapely.geometry").attr("Point");
    for (const auto& node : result.nodes) { id.append(node.node_id); row.append(node.row); col.append(node.col); elevation.append(node.elevation); degree.append(node.degree); endpoint.append(node.is_endpoint); junction.append(node.is_junction); outlet.append(node.is_outlet); distance.append(node.downstream_distance_metres); geometry.append(point(node.x, node.y)); }
    data["node_id"] = id; data["row"] = row; data["col"] = col; data["elevation"] = elevation; data["degree"] = degree; data["is_endpoint"] = endpoint; data["is_junction"] = junction; data["is_outlet"] = outlet; data["downstream_distance_metres"] = distance;
    return nb::module_::import_("geopandas").attr("GeoDataFrame")(data, nb::arg("geometry") = geometry, nb::arg("crs") = crs_wkt);
}

nb::object section_geodataframe(const lsdtt::lsdtt_tidal_skeleton_metrics_result& result, const std::string& crs_wkt)
{
    nb::dict data; nb::list id, upstream, downstream, length, drop, strahler, hack, loop, distance, geometry;
    nb::object line = nb::module_::import_("shapely.geometry").attr("LineString");
    for (const auto& section : result.sections) { nb::list coordinates; for (const auto& xy : section.xy) coordinates.append(nb::make_tuple(xy.x, xy.y)); id.append(section.section_id); upstream.append(section.upstream_node_id); downstream.append(section.downstream_node_id); length.append(section.total_length_metres); drop.append(section.elevation_drop_metres); strahler.append(section.strahler_order); hack.append(section.hack_order); loop.append(section.is_loop); distance.append(section.downstream_distance_at_downstream_node_metres); geometry.append(line(coordinates)); }
    data["section_id"] = id; data["upstream_node_id"] = upstream; data["downstream_node_id"] = downstream; data["length_metres"] = length; data["elevation_drop_metres"] = drop; data["strahler_order"] = strahler; data["hack_order"] = hack; data["is_loop"] = loop; data["downstream_distance_metres"] = distance;
    return nb::module_::import_("geopandas").attr("GeoDataFrame")(data, nb::arg("geometry") = geometry, nb::arg("crs") = crs_wkt);
}

class TidalSkeletonMetricsResult {
public:
    TidalSkeletonMetricsResult(lsdtt::lsdtt_tidal_skeleton_metrics_result result, std::vector<std::string> strip_names, std::string crs_wkt) : result_(std::move(result)), strip_names_(std::move(strip_names)), crs_wkt_(std::move(crs_wkt)) {}
    Raster node_ids() const { return Raster(result_.node_id_raster); }
    Raster section_ids() const { return Raster(result_.section_id_raster); }
    Raster strahler_order() const { return Raster(result_.section_strahler_order_raster); }
    Raster hack_order() const { return Raster(result_.section_hack_order_raster); }
    Raster endpoints() const { return Raster(result_.endpoint_mask); }
    Raster junctions() const { return Raster(result_.junction_mask); }
    nb::object nearest_station_ids() const { return optional_raster(result_.nearest_station_id_raster); }
    nb::object downstream_distance() const { return optional_raster(result_.station_downstream_distance_raster); }
    nb::object nodes() const { return node_geodataframe(result_, crs_wkt_); }
    nb::object sections() const { return section_geodataframe(result_, crs_wkt_); }
    nb::object stations() const { return station_geodataframe(result_, strip_names_, crs_wkt_); }
    int n_nodes() const { return static_cast<int>(result_.nodes.size()); }
    int n_sections() const { return static_cast<int>(result_.sections.size()); }
    int n_stations() const { return static_cast<int>(result_.stations.size()); }
private:
    lsdtt::lsdtt_tidal_skeleton_metrics_result result_;
    std::vector<std::string> strip_names_;
    std::string crs_wkt_;
};

TidalSkeletonMetricsResult analyse_tidal_skeleton(const Raster& dem, const Raster& channel_mask,
                                                   const Raster& skeleton, nb::object outlet_mask,
                                                   const std::string& selected_outlet_mode,
                                                   double station_spacing_metres, bool compute_cross_sections,
                                                   int tangent_support_points, double cross_section_half_length_metres,
                                                   double cross_section_sample_spacing_metres, double max_width_ratio,
                                                   bool compute_station_strip_metrics, nb::object strip_values, int n_threads)
{
    require_same_grid(dem, channel_mask, "channel_mask"); require_same_grid(dem, skeleton, "skeleton");
    const Raster* outlet = outlet_mask.is_none() ? nullptr : &nb::cast<const Raster&>(outlet_mask);
    if (outlet) require_same_grid(dem, *outlet, "outlet_mask");
    lsdtt::lsdtt_tidal_skeleton_metrics_request request;
    request.outlet_mode = outlet_mode(selected_outlet_mode);
    request.station_spacing_metres = station_spacing_metres;
    request.compute_cross_sections = compute_cross_sections;
    request.tangent_support_points = tangent_support_points;
    request.cross_section_half_length_metres = cross_section_half_length_metres;
    request.cross_section_sample_spacing_metres = cross_section_sample_spacing_metres;
    request.max_width_ratio = max_width_ratio;
    request.compute_station_strip_metrics = compute_station_strip_metrics;
    request.n_threads = n_threads;
    std::vector<std::string> strip_names;
    if (!strip_values.is_none()) {
        const nb::dict values = nb::cast<nb::dict>(strip_values);
        for (const auto& item : values) {
            const std::string name = nb::cast<std::string>(item.first);
            const Raster& value = nb::cast<const Raster&>(item.second);
            require_same_grid(dem, value, ("strip_values['" + name + "']").c_str());
            request.strip_value_rasters.push_back(&value.native()); request.strip_value_names.push_back(name); strip_names.push_back(name);
        }
    }
    return TidalSkeletonMetricsResult(lsdtt::analyse_tidal_skeleton(dem.native(), channel_mask.native(), skeleton.native(), outlet ? &outlet->native() : nullptr, request), std::move(strip_names), dem.native().get_meta().crs_wkt);
}

} // namespace

void bind_tidal_workflows(nb::module_& module)
{
    nb::class_<TidalChannelResult>(module, "TidalChannelResult", "In-memory tidal-channel extraction products.")
        .def_prop_ro("channel_mask", &TidalChannelResult::channel_mask)
        .def_prop_ro("upl", &TidalChannelResult::upl)
        .def_prop_ro("skeleton", &TidalChannelResult::skeleton)
        .def_prop_ro("connected_skeleton", &TidalChannelResult::connected_skeleton)
        .def_prop_ro("skeleton_bridge_mask", &TidalChannelResult::skeleton_bridge_mask)
        .def_prop_ro("median_rasters", &TidalChannelResult::median_rasters)
        .def_prop_ro("residual_rasters", &TidalChannelResult::residual_rasters)
        .def_prop_ro("threshold_masks", &TidalChannelResult::threshold_masks)
        .def_prop_ro("scales_metres", &TidalChannelResult::scales_metres)
        .def_prop_ro("thresholds_metres", &TidalChannelResult::thresholds_metres);
    module.def("find_tidal_channels", &find_tidal_channels, nb::arg("dem"),
        nb::arg("use_flood_mask") = false, nb::arg("flood_elevation_metres") = 0.,
        nb::arg("median_window_radii_metres") = std::vector<double>{5., 10., 20.},
        nb::arg("residual_thresholds_metres") = std::vector<double>{.1}, nb::arg("scale_mode") = "any",
        nb::arg("apply_binary_closing") = false, nb::arg("binary_closing_radius_metres") = 0.,
        nb::arg("fill_holes") = false, nb::arg("min_channel_component_pixels") = 1,
        nb::arg("extract_skeleton") = true, nb::arg("skeleton_spur_max_pixels") = 10,
        nb::arg("skeleton_terminal_loop_max_pixels") = 10, nb::arg("skeleton_min_component_pixels") = 1,
        nb::arg("connect_skeleton_segments") = false, nb::arg("connect_max_gap_pixels") = 20.,
        nb::arg("connect_max_facing_angle_degrees") = 35., nb::arg("connect_max_heading_misalignment_degrees") = 45.,
        nb::arg("connect_max_candidate_score") = 1., nb::arg("connect_route_corridor_half_width_pixels") = 12.,
        nb::arg("compute_upl") = true, nb::arg("n_threads") = 0,
        R"doc(Extract tidal channels as multiscale local topographic lows or a flood mask.

At each scale the residual ``moving-window median(dem) - dem`` is computed;
pixels lying at least the residual threshold below their local median are
channel (a TidalGeoPro-style channelness signal). Optional cleanup,
skeletonisation and gap bridging follow.

Parameters
----------
dem : Raster
    Elevation raster (m) of a marsh or tidal flat.
use_flood_mask : bool, default=False
    Skip the median workflow: pixels with elevation <= the flood elevation
    are channel.
flood_elevation_metres : float, default=0.0
    Flood elevation (m) for ``use_flood_mask``.
median_window_radii_metres : list of float, default=[5.0, 10.0, 20.0]
    Radii (m) of the median windows, one per scale.
residual_thresholds_metres : list of float, default=[0.1]
    Residual (m) a pixel must reach below the median; one value for all
    scales or one per scale.
scale_mode : {"any", "all"}, default="any"
    A pixel is channel if it passes at any scale or at every scale.
apply_binary_closing : bool, default=False
    Close the channel mask with a disk.
binary_closing_radius_metres : float, default=0.0
    Radius (m) of that disk.
fill_holes : bool, default=False
    Fill enclosed holes in the channel mask.
min_channel_component_pixels : int, default=1
    Smallest channel component kept (pixels).
extract_skeleton : bool, default=True
    Thin the channel mask to a one-pixel skeleton.
skeleton_spur_max_pixels : int, default=10
    Remove skeleton spurs up to this length (pixels).
skeleton_terminal_loop_max_pixels : int, default=10
    Remove terminal skeleton loops up to this size (pixels).
skeleton_min_component_pixels : int, default=1
    Smallest skeleton component kept (pixels).
connect_skeleton_segments : bool, default=False
    Bridge gaps between facing skeleton end points.
connect_max_gap_pixels : float, default=20.0
    Largest gap (pixels) bridged.
connect_max_facing_angle_degrees : float, default=35.0
    Largest angle (degrees) between an end point's heading and the gap.
connect_max_heading_misalignment_degrees : float, default=45.0
    Largest misalignment (degrees) between the two end headings.
connect_max_candidate_score : float, default=1.0
    Largest combined (normalised) score of an accepted bridge.
connect_route_corridor_half_width_pixels : float, default=12.0
    Half-width (pixels) of the corridor in which a bridge is routed.
compute_upl : bool, default=True
    Compute ``upl``: the Euclidean distance (m) from each pixel to the
    nearest channel pixel.
n_threads : int, default=0
    OpenMP thread count; 0 uses the runtime default.

Returns
-------
TidalChannelResult
)doc");
    nb::class_<TidalSkeletonMetricsResult>(module, "TidalSkeletonMetricsResult", "Tidal skeleton topology and station metrics.")
        .def_prop_ro("node_ids", &TidalSkeletonMetricsResult::node_ids)
        .def_prop_ro("section_ids", &TidalSkeletonMetricsResult::section_ids)
        .def_prop_ro("strahler_order", &TidalSkeletonMetricsResult::strahler_order)
        .def_prop_ro("hack_order", &TidalSkeletonMetricsResult::hack_order)
        .def_prop_ro("endpoints", &TidalSkeletonMetricsResult::endpoints)
        .def_prop_ro("junctions", &TidalSkeletonMetricsResult::junctions)
        .def_prop_ro("nearest_station_ids", &TidalSkeletonMetricsResult::nearest_station_ids)
        .def_prop_ro("downstream_distance", &TidalSkeletonMetricsResult::downstream_distance)
        .def_prop_ro("nodes", &TidalSkeletonMetricsResult::nodes)
        .def_prop_ro("sections", &TidalSkeletonMetricsResult::sections)
        .def_prop_ro("stations", &TidalSkeletonMetricsResult::stations)
        .def_prop_ro("n_nodes", &TidalSkeletonMetricsResult::n_nodes)
        .def_prop_ro("n_sections", &TidalSkeletonMetricsResult::n_sections)
        .def_prop_ro("n_stations", &TidalSkeletonMetricsResult::n_stations);
    module.def("analyse_tidal_skeleton", &analyse_tidal_skeleton, nb::arg("dem"), nb::arg("channel_mask"), nb::arg("skeleton"),
        nb::arg("outlet_mask") = nb::none(), nb::arg("outlet_mode") = "edge_nodes", nb::arg("station_spacing_metres") = 5.,
        nb::arg("compute_cross_sections") = true, nb::arg("tangent_support_points") = 2,
        nb::arg("cross_section_half_length_metres") = 25., nb::arg("cross_section_sample_spacing_metres") = 1.,
        nb::arg("max_width_ratio") = 3., nb::arg("compute_station_strip_metrics") = true,
        nb::arg("strip_values") = nb::none(), nb::arg("n_threads") = 0,
        R"doc(Analyse an existing tidal-channel skeleton into sections, orders, and station metrics.

Parameters
----------
dem : Raster
    Elevation raster (m).
channel_mask : Raster
    0/1 channel mask on the same grid (``TidalChannelResult.channel_mask``).
skeleton : Raster
    0/1 skeleton on the same grid (``TidalChannelResult.skeleton``).
outlet_mask : Raster or None, default=None
    0/1 outlet mask, needed for ``outlet_mode="outlet_mask"``.
outlet_mode : {"edge_nodes", "endpoint_nodes", "outlet_mask", "lowest_endpoint"}, default="edge_nodes"
    Outlets for downstream distances: skeleton end points touching the
    raster edge, all end points, the masked nodes, or the single lowest
    end point.
station_spacing_metres : float, default=5.0
    Spacing (m) of measurement stations along each section (positive).
compute_cross_sections : bool, default=True
    Cast cross sections at each station to measure width.
tangent_support_points : int, default=2
    Neighbouring stations used to estimate the local channel direction.
cross_section_half_length_metres : float, default=25.0
    Half-length (m) of each cross section.
cross_section_sample_spacing_metres : float, default=1.0
    Sample spacing (m) along a cross section.
max_width_ratio : float, default=3.0
    Reject cross sections implausibly long relative to local width; <= 0
    disables.
compute_station_strip_metrics : bool, default=True
    Assign every pixel to its nearest station and compute strip metrics.
strip_values : dict of str to Raster or None, default=None
    Extra co-registered rasters integrated over the station strips.
n_threads : int, default=0
    OpenMP thread count; 0 uses the runtime default.

Returns
-------
TidalSkeletonMetricsResult
)doc");
}
} // namespace lsdtt3_python
