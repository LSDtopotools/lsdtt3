#include <algorithm>
#include <stdexcept>
#include <string>
#include <cstdint>

#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>

#include <gdal_priv.h>

#include <core/lsdtt_secret.hpp>
#include <core/lsdtt_dem_grouping.hpp>
#include <cpu/lsdtt_outlet_basin_verification_cpu.hpp>

#include "binding_common.hpp"

namespace lsdtt3_python {
namespace {

lsdtt::lsdtt_dem_source source_from_name(const std::string& source)
{
    lsdtt::lsdtt_dem_source parsed;
    if (!lsdtt::parse_dem_source(source, parsed))
    {
        throw std::invalid_argument(
            "source must be cop30_aws, opentopography, or planetary_computer"
        );
    }
    return parsed;
}

Raster fetch_dem(
    double west, double south, double east, double north,
    const std::string& source, int target_epsg,
    double grid_spacing, const std::string& resampling, double nodata_value,
    int size_pixels, const std::string& cache_dir, const std::string& ot_dem_type,
    const std::string& pc_collection, const std::string& pc_asset,
    int http_timeout_seconds, int http_connect_timeout_seconds, int http_max_retry)
{
    GDALAllRegister();
    lsdtt::lsdtt_fetch_options options;
    options.west = west;
    options.south = south;
    options.east = east;
    options.north = north;
    options.source = source_from_name(source);
    options.target_epsg = target_epsg;
    options.grid_spacing = grid_spacing;
    options.resampling = resampling;
    options.no_data_value = nodata_value;
    options.size_px = size_pixels;
    options.cache_dir = cache_dir;
    options.ot_dem_type = ot_dem_type;
    options.pc_collection = pc_collection;
    options.pc_asset = pc_asset;
    options.http_timeout_seconds = http_timeout_seconds;
    options.http_connect_timeout_seconds = http_connect_timeout_seconds;
    options.http_max_retry = http_max_retry;

    if (options.source == lsdtt::lsdtt_dem_source::opentopography)
    {
        const lsdtt::lsdtt_secret_result secret = lsdtt::resolve_secret({
            "OT_API_KEY", "OT_API_KEY_FILE", "OpenTopography API key", false
        });
        if (secret.found) options.api_key = secret.value;
    }
    static std::uint64_t request_id = 0;
    const std::string temporary_path =
        "/vsimem/lsdtt3-fetch-" + std::to_string(++request_id) + ".tif";
    const lsdtt::lsdtt_fetch_result result =
        lsdtt::fetch_dem(options, temporary_path);
    if (!result.ok)
    {
        throw std::runtime_error(result.message);
    }

    try
    {
        Raster raster = Raster::read(temporary_path);
        VSIUnlink(temporary_path.c_str());
        return raster;
    }
    catch (...)
    {
        VSIUnlink(temporary_path.c_str());
        throw;
    }
}

Raster fetch_with_options(const lsdtt::lsdtt_fetch_options& options)
{
    static std::uint64_t request_id = 0;
    const std::string temporary_path =
        "/vsimem/lsdtt3-outlet-fetch-" + std::to_string(++request_id) + ".tif";
    const lsdtt::lsdtt_fetch_result result =
        lsdtt::fetch_dem(options, temporary_path);
    if (!result.ok) throw std::runtime_error(result.message);
    try
    {
        Raster raster = Raster::read(temporary_path);
        VSIUnlink(temporary_path.c_str());
        return raster;
    }
    catch (...)
    {
        VSIUnlink(temporary_path.c_str());
        throw;
    }
}

Raster basin_label_raster(
    const Raster& reference,
    const lsdtt::lsdtt_outlet_basin_verification_result& verification)
{
    const auto& meta = reference.native().get_meta();
    std::vector<float> labels(verification.basin_labels.size(), -1.0F);
    std::vector<std::uint8_t> mask(verification.basin_labels.size(), 0);
    for (std::size_t i = 0; i < labels.size(); ++i)
        if (verification.basin_labels[i] >= 0)
        {
            labels[i] = static_cast<float>(verification.basin_labels[i]);
            mask[i] = 1;
        }
    auto label_meta = meta;
    label_meta.no_data_value = -1.0;
    return Raster(lsdtt::lsdtt_raster_cpu(std::move(labels), std::move(mask), label_meta));
}

Raster mask_to_selected_basins(
    const Raster& raster,
    const lsdtt::lsdtt_outlet_basin_verification_result& verification)
{
    auto data = raster.native().get_data();
    auto mask = raster.native().get_mask();
    for (std::size_t i = 0; i < mask.size(); ++i)
        if (verification.basin_labels[i] < 0)
        {
            data[i] = static_cast<float>(raster.native().get_meta().no_data_value);
            mask[i] = 0;
        }
    return Raster(lsdtt::lsdtt_raster_cpu(
        std::move(data), std::move(mask), raster.native().get_meta()));
}

std::vector<lsdtt::lsdtt_outlet_point_wgs84> wgs84_outlet_points(nb::object outlets)
{
    nb::object geopandas = nb::module_::import_("geopandas");
    if (!nb::isinstance(outlets, geopandas.attr("GeoDataFrame")))
        throw std::invalid_argument("outlets must be a GeoDataFrame of point geometries");
    if (outlets.attr("crs").is_none())
        throw std::invalid_argument("outlets must define a CRS");

    nb::object wgs84 = outlets.attr("to_crs")("EPSG:4326");
    nb::object geometry = wgs84.attr("geometry");
    std::vector<lsdtt::lsdtt_outlet_point_wgs84> points;
    const std::size_t n = nb::cast<std::size_t>(wgs84.attr("shape")[0]);
    points.reserve(n);
    for (std::size_t i = 0; i < n; ++i)
    {
        nb::object point = geometry[i];
        if (point.is_none() || nb::cast<std::string>(point.attr("geom_type")) != "Point")
            throw std::invalid_argument("outlets must contain only non-empty Point geometries");
        points.push_back({static_cast<int>(i), nb::cast<double>(point.attr("x")),
                          nb::cast<double>(point.attr("y"))});
    }
    if (points.empty()) throw std::invalid_argument("outlets must not be empty");
    return points;
}

class OutletFetchResult
{
public:
    OutletFetchResult(std::vector<Raster> rasters,
                      lsdtt::lsdtt_outlet_bounds_discovery_result discovery,
                      std::vector<lsdtt::lsdtt_dem_group> groups,
                      std::vector<lsdtt::lsdtt_outlet_basin_verification_result> verification,
                      std::vector<Raster> basin_labels)
        : rasters_(std::move(rasters)), discovery_(std::move(discovery)), groups_(std::move(groups)),
          verification_(std::move(verification)), basin_labels_(std::move(basin_labels)) {}

    const std::vector<Raster>& rasters() const noexcept { return rasters_; }
    const lsdtt::lsdtt_outlet_bounds_discovery_result& discovery() const noexcept { return discovery_; }
    const std::vector<lsdtt::lsdtt_dem_group>& groups() const noexcept { return groups_; }
    const std::vector<lsdtt::lsdtt_outlet_basin_verification_result>& verification() const noexcept
    { return verification_; }
    const std::vector<Raster>& basin_labels() const noexcept { return basin_labels_; }

private:
    std::vector<Raster> rasters_;
    lsdtt::lsdtt_outlet_bounds_discovery_result discovery_;
    std::vector<lsdtt::lsdtt_dem_group> groups_;
    std::vector<lsdtt::lsdtt_outlet_basin_verification_result> verification_;
    std::vector<Raster> basin_labels_;
};

nb::object snapped_outlet_table(const OutletFetchResult& result)
{
    nb::dict columns;
    nb::list group_ids, outlet_ids, rows, cols, x, y;
    for (std::size_t i = 0; i < result.verification().size(); ++i)
    {
        const auto& meta = result.rasters()[i].native().get_meta();
        for (const auto& outlet : result.verification()[i].snapped_outlets)
        {
            group_ids.append(result.groups()[i].group_id); outlet_ids.append(outlet.outlet_id);
            rows.append(outlet.snapped_row); cols.append(outlet.snapped_col);
            x.append(meta.x_min + (static_cast<double>(outlet.snapped_col) + 0.5) * meta.cell_size);
            y.append(meta.y_max() - (static_cast<double>(outlet.snapped_row) + 0.5) * meta.cell_size);
        }
    }
    columns["group_id"] = std::move(group_ids); columns["outlet_id"] = std::move(outlet_ids);
    columns["row"] = std::move(rows); columns["col"] = std::move(cols);
    columns["x"] = std::move(x); columns["y"] = std::move(y);
    return nb::module_::import_("pandas").attr("DataFrame")(columns);
}

nb::object outlet_table(const OutletFetchResult& result)
{
    nb::dict columns;
    nb::list ids, longitude, latitude, status, iterations, west, south, east, north, message;
    for (const auto& record : result.discovery().records)
    {
        ids.append(record.point_id); longitude.append(record.longitude); latitude.append(record.latitude);
        status.append(record.status); iterations.append(record.iterations); west.append(record.west);
        south.append(record.south); east.append(record.east); north.append(record.north);
        message.append(record.message);
    }
    columns["outlet_id"] = std::move(ids); columns["longitude"] = std::move(longitude);
    columns["latitude"] = std::move(latitude); columns["status"] = std::move(status);
    columns["iterations"] = std::move(iterations); columns["west"] = std::move(west);
    columns["south"] = std::move(south); columns["east"] = std::move(east);
    columns["north"] = std::move(north); columns["message"] = std::move(message);
    return nb::module_::import_("pandas").attr("DataFrame")(columns);
}

nb::object group_table(const OutletFetchResult& result)
{
    nb::dict columns;
    nb::list ids, west, south, east, north, area, outlet_ids, verified, messages;
    for (std::size_t i = 0; i < result.groups().size(); ++i)
    {
        const auto& group = result.groups()[i];
        const auto& check = result.verification()[i];
        nb::list members;
        for (std::size_t index : group.record_indices)
            members.append(result.discovery().records[index].point_id);
        ids.append(group.group_id); west.append(group.west); south.append(group.south);
        east.append(group.east); north.append(group.north); area.append(group.area_km2);
        outlet_ids.append(std::move(members));
        verified.append(check.ok); messages.append(check.message);
    }
    columns["group_id"] = std::move(ids); columns["west"] = std::move(west);
    columns["south"] = std::move(south); columns["east"] = std::move(east);
    columns["north"] = std::move(north); columns["area_km2"] = std::move(area);
    columns["outlet_ids"] = std::move(outlet_ids);
    columns["basins_verified"] = std::move(verified);
    columns["verification_message"] = std::move(messages);
    return nb::module_::import_("pandas").attr("DataFrame")(columns);
}

OutletFetchResult fetch_dem_for_outlets(
    nb::object outlets, const std::string& source, int target_epsg, double grid_spacing,
    double initial_radius_km, bool auto_split_disconnected_basins,
    double split_cluster_gap_km, double max_final_dem_area_km2, double split_overlap_km,
    int max_basins_per_dem, const std::string& cache_dir, const std::string& resampling,
    double nodata_value, const std::string& snapping_method,
    int contributing_pixels_threshold, int edge_buffer_pixels,
    int nodata_buffer_pixels, int tight_bbox_buffer_pixels,
    double verify_expand_factor, int verify_max_iters, bool tight_bbox_refetch,
    int min_search_radius_pixels, int max_search_radius_pixels, bool mask_outside_basins)
{
    GDALAllRegister();
    const lsdtt::lsdtt_dem_source parsed_source = source_from_name(source);
    lsdtt::lsdtt_outlet_bounds_discovery_options discovery_options;
    discovery_options.initial_radius_km = initial_radius_km;
    const auto discovery = lsdtt::discover_outlet_bounds(
        wgs84_outlet_points(outlets), discovery_options);
    if (!discovery.ok) throw std::runtime_error(discovery.message);

    lsdtt::lsdtt_dem_grouping_options grouping_options;
    grouping_options.auto_split_disconnected_basins = auto_split_disconnected_basins;
    grouping_options.split_cluster_gap_km = split_cluster_gap_km;
    grouping_options.max_final_dem_area_km2 = max_final_dem_area_km2;
    grouping_options.split_overlap_km = split_overlap_km;
    grouping_options.max_basins_per_dem = max_basins_per_dem;
    auto groups = lsdtt::group_outlet_bounds(discovery.records, grouping_options);
    if (groups.empty()) throw std::runtime_error("could not form any outlet DEM groups");

    std::vector<Raster> rasters;
    std::vector<lsdtt::lsdtt_outlet_basin_verification_result> verification;
    std::vector<Raster> basin_labels;
    rasters.reserve(groups.size());
    verification.reserve(groups.size());
    basin_labels.reserve(groups.size());
    for (auto& group : groups)
    {
        lsdtt::lsdtt_fetch_options options;
        options.west = group.west; options.south = group.south;
        options.east = group.east; options.north = group.north;
        options.source = parsed_source;
        options.target_epsg = target_epsg;
        options.grid_spacing = grid_spacing;
        options.resampling = resampling;
        options.no_data_value = nodata_value;
        options.cache_dir = cache_dir;
        if (parsed_source == lsdtt::lsdtt_dem_source::opentopography)
        {
            const lsdtt::lsdtt_secret_result secret = lsdtt::resolve_secret({
                "OT_API_KEY", "OT_API_KEY_FILE", "OpenTopography API key", false
            });
            if (secret.found) options.api_key = secret.value;
        }
        if (contributing_pixels_threshold <= 0 || verify_max_iters < 1 ||
            min_search_radius_pixels < 0 ||
            max_search_radius_pixels < min_search_radius_pixels ||
            !(verify_expand_factor > 1.0))
            throw std::invalid_argument(
                "invalid verification or snapping options");
        lsdtt::lsdtt_outlet_basin_verification_options verify_options;
        verify_options.snapping_method = snapping_method;
        verify_options.contributing_pixels_threshold =
            static_cast<std::uint64_t>(contributing_pixels_threshold);
        verify_options.min_search_radius_pixels = min_search_radius_pixels;
        verify_options.max_search_radius_pixels = max_search_radius_pixels;
        verify_options.required_edge_buffer_pixels = edge_buffer_pixels;
        verify_options.required_nodata_buffer_pixels = nodata_buffer_pixels;
        verify_options.tight_bbox_buffer_pixels = tight_bbox_buffer_pixels;

        Raster raster = fetch_with_options(options);
        lsdtt::lsdtt_outlet_basin_verification_result check;
        for (int attempt = 1; attempt <= verify_max_iters; ++attempt)
        {
            check = lsdtt::verify_outlet_group_basins(
                raster.native(), group, discovery.records, verify_options);
            if (check.ok) break;
            if (attempt == verify_max_iters)
                throw std::runtime_error("outlet group " + std::to_string(group.group_id) +
                    " is not safely contained after " + std::to_string(attempt) +
                    " fetch attempts: " + check.message);
            const double centre_x = 0.5 * (options.west + options.east);
            const double centre_y = 0.5 * (options.south + options.north);
            const double half_width = 0.5 * (options.east - options.west) * verify_expand_factor;
            const double half_height = 0.5 * (options.north - options.south) * verify_expand_factor;
            options.west = std::max(-180.0, centre_x - half_width);
            options.east = std::min(180.0, centre_x + half_width);
            options.south = std::max(-90.0, centre_y - half_height);
            options.north = std::min(90.0, centre_y + half_height);
            raster = fetch_with_options(options);
        }
        if (tight_bbox_refetch && check.tight_bbox_valid)
        {
            const double current_area = lsdtt::approximate_bbox_area_km2(
                options.west, options.south, options.east, options.north);
            const double tight_area = lsdtt::approximate_bbox_area_km2(
                check.tight_west, check.tight_south, check.tight_east, check.tight_north);
            bool contains_input_outlets = true;
            for (const std::size_t index : group.record_indices)
            {
                const auto& outlet = discovery.records[index];
                contains_input_outlets = contains_input_outlets &&
                    outlet.longitude >= check.tight_west && outlet.longitude <= check.tight_east &&
                    outlet.latitude >= check.tight_south && outlet.latitude <= check.tight_north;
            }
            if (contains_input_outlets && tight_area < current_area * 0.98)
            {
                options.west = check.tight_west; options.south = check.tight_south;
                options.east = check.tight_east; options.north = check.tight_north;
                raster = fetch_with_options(options);
                check = lsdtt::verify_outlet_group_basins(
                    raster.native(), group, discovery.records, verify_options);
                if (!check.ok) throw std::runtime_error(
                    "tight outlet DEM refetch is not safely contained: " + check.message);
            }
        }
        group.west = options.west; group.south = options.south;
        group.east = options.east; group.north = options.north;
        group.area_km2 = lsdtt::approximate_bbox_area_km2(
            group.west, group.south, group.east, group.north);
        basin_labels.push_back(basin_label_raster(raster, check));
        if (mask_outside_basins) raster = mask_to_selected_basins(raster, check);
        rasters.push_back(std::move(raster));
        verification.push_back(std::move(check));
    }
    return OutletFetchResult(
        std::move(rasters), discovery, groups, std::move(verification), std::move(basin_labels));
}

nb::list opentopography_datasets()
{
    nb::list result;
    for (const auto& dataset : lsdtt::ot_dataset_catalogue())
    {
        nb::dict item;
        item["demtype"] = dataset.demtype;
        item["name"] = dataset.name;
        item["resolution"] = dataset.resolution;
        item["coverage"] = dataset.coverage;
        item["description"] = dataset.description;
        item["citation"] = dataset.citation;
        item["max_area_km2"] = dataset.max_area_km2;
        result.append(std::move(item));
    }
    return result;
}

} // namespace

void bind_raster_fetch(nb::module_& module)
{
    nb::class_<OutletFetchResult>(module, "OutletFetchResult", R"doc(
In-memory DEM groups fetched around outlet points.

``rasters`` and ``groups`` share order: raster ``i`` is the DEM for group
``groups.iloc[i]``. ``outlets`` records the initial WGS84 envelopes used to
build those groups. ``basin_labels`` and ``snapped_outlets`` describe the
verified selected watersheds and their snapped channel locations.
)doc")
        .def_prop_ro("rasters", [](const OutletFetchResult& result) {
            nb::list rasters;
            for (const Raster& raster : result.rasters()) rasters.append(raster);
            return rasters;
        })
        .def_prop_ro("basin_labels", [](const OutletFetchResult& result) {
            nb::list labels;
            for (const Raster& raster : result.basin_labels()) labels.append(raster);
            return labels;
        })
        .def_prop_ro("outlets", &outlet_table)
        .def_prop_ro("snapped_outlets", &snapped_outlet_table)
        .def_prop_ro("groups", &group_table);

    nb::class_<lsdtt::lsdtt_bounds_from_vector_result>(module, "VectorBounds", R"doc(
WGS84 bounds derived from a vector layer.

The input may contain polygons, or exactly two points defining opposite box
corners. Geometry is transformed to WGS84 before its bounds are calculated.
)doc")
        .def_ro("ok", &lsdtt::lsdtt_bounds_from_vector_result::ok)
        .def_ro("west", &lsdtt::lsdtt_bounds_from_vector_result::west)
        .def_ro("south", &lsdtt::lsdtt_bounds_from_vector_result::south)
        .def_ro("east", &lsdtt::lsdtt_bounds_from_vector_result::east)
        .def_ro("north", &lsdtt::lsdtt_bounds_from_vector_result::north)
        .def_ro("layer_name", &lsdtt::lsdtt_bounds_from_vector_result::layer_name)
        .def_ro("geometry_kind", &lsdtt::lsdtt_bounds_from_vector_result::geometry_kind)
        .def_ro("n_features", &lsdtt::lsdtt_bounds_from_vector_result::n_features)
        .def_ro("source_crs", &lsdtt::lsdtt_bounds_from_vector_result::source_crs)
        .def_ro("warning", &lsdtt::lsdtt_bounds_from_vector_result::warning)
        .def_ro("message", &lsdtt::lsdtt_bounds_from_vector_result::message)
        .def_prop_ro("bbox", [](const lsdtt::lsdtt_bounds_from_vector_result& result) {
            return nb::make_tuple(result.west, result.south, result.east, result.north);
        });

    module.def("fetch_dem", &fetch_dem,
        nb::arg("west"), nb::arg("south"), nb::arg("east"), nb::arg("north"),
        nb::arg("source") = "cop30_aws",
        nb::arg("target_epsg") = 0, nb::arg("grid_spacing") = 30.0,
        nb::arg("resampling") = "cubic", nb::arg("nodata_value") = -9999.0,
        nb::arg("size_pixels") = 0, nb::arg("cache_dir") = "",
        nb::arg("ot_dem_type") = "COP30", nb::arg("pc_collection") = "cop-dem-glo-30",
        nb::arg("pc_asset") = "data", nb::arg("http_timeout_seconds") = 60,
        nb::arg("http_connect_timeout_seconds") = 15, nb::arg("http_max_retry") = 3,
        R"doc(
Fetch a DEM for a WGS84 bounding box as an in-memory projected raster.

Parameters
----------
west, south, east, north : float
    WGS84 bounds in decimal degrees. Antimeridian-crossing requests are not
    supported: ``west`` must be less than ``east``.
source : {"cop30_aws", "opentopography", "planetary_computer"}, optional
    ``cop30_aws`` is a keyless global 30 m default. OpenTopography reads its
    key only from ``OT_API_KEY`` or ``OT_API_KEY_FILE``.
target_epsg : int, optional
    Output EPSG code. Zero selects UTM from the bounding-box centre.
grid_spacing : float, optional
    Output grid spacing in projected units.

cache_dir : str, optional
    Persistent cache directory. This is the only optional disk output; the
    returned raster itself is kept in memory.

Returns
-------
Raster
    A projected raster. Network, GDAL, and invalid-request failures raise
    ``RuntimeError``.
)doc");

    module.def("fetch_dem_for_outlets", &fetch_dem_for_outlets,
        nb::arg("outlets"), nb::arg("source") = "cop30_aws",
        nb::arg("target_epsg") = 0, nb::arg("grid_spacing") = 30.0,
        nb::arg("initial_radius_km") = 100.0,
        nb::arg("auto_split_disconnected_basins") = true,
        nb::arg("split_cluster_gap_km") = 25.0,
        nb::arg("max_final_dem_area_km2") = 2500.0,
        nb::arg("split_overlap_km") = 2.0, nb::arg("max_basins_per_dem") = 0,
        nb::arg("cache_dir") = "", nb::arg("resampling") = "cubic",
        nb::arg("nodata_value") = -9999.0,
        nb::arg("snapping_method") = "downstream",
        nb::arg("contributing_pixels_threshold") = 1000,
        nb::arg("edge_buffer_pixels") = 3,
        nb::arg("nodata_buffer_pixels") = 1,
        nb::arg("tight_bbox_buffer_pixels") = 3,
        nb::arg("verify_expand_factor") = 1.25,
        nb::arg("verify_max_iters") = 3,
        nb::arg("tight_bbox_refetch") = true,
        nb::arg("min_search_radius_pixels") = 0,
        nb::arg("max_search_radius_pixels") = 25,
        nb::arg("mask_outside_basins") = false,
        R"doc(
Fetch one or more in-memory DEMs around WGS84 outlet points.

The initial radius around each point is converted to WGS84 bounds, nearby
outlets are grouped, and large groups are split before each group is fetched.
The input GeoDataFrame is reprojected to WGS84 internally. Each group is
verified with D8 basin delineation, expanded if a basin reaches an edge or
nodata, and optionally refetched to its tight basin envelope.

Parameters
----------
outlets : geopandas.GeoDataFrame
    Point geometries defining outlet locations, in any declared CRS.
initial_radius_km : float, optional
    Initial half-width around each outlet before grouping.
max_final_dem_area_km2 : float, optional
    Largest allowed grouped initial extent before splitting it.
snapping_method : {"none", "expanding_square", "downstream"}, optional
    Outlet snapping used before each basin delineation.
max_search_radius_pixels : int, optional
    Largest expanding-square snap radius when ``snapping_method`` is
    ``"expanding_square"``.
edge_buffer_pixels, nodata_buffer_pixels : int, optional
    Required basin clearance from the DEM edge and invalid cells.
tight_bbox_refetch : bool, optional
    Refetch a smaller final DEM around the verified basin footprint.
mask_outside_basins : bool, optional
    Set cells outside the selected basin union to nodata in returned DEMs.

Returns
-------
OutletFetchResult
    Grouped raster collection and copied outlet/group diagnostic tables.
)doc");

    module.def("bounds_from_vector", [](const std::string& path, const std::string& layer) {
        GDALAllRegister();
        return lsdtt::bounds_from_vector_file(path, layer);
    }, nb::arg("path"), nb::arg("layer") = "", R"doc(
Derive WGS84 ``(west, south, east, north)`` bounds from a vector file.

The layer must contain polygons or exactly two points. Use the result's
``bbox`` directly as keyword arguments to :func:`fetch_dem`.
)doc");

    module.def("http_available", &lsdtt::http_available,
        "Return whether the bundled GDAL has HTTP support.");
    module.def("utm_epsg_for_point", &lsdtt::utm_epsg_for_point,
        nb::arg("longitude"), nb::arg("latitude"),
        "Return the UTM EPSG code containing a WGS84 point.");
    module.def("cop30_tiles_for_bbox", &lsdtt::cop30_tiles_for_bbox,
        nb::arg("west"), nb::arg("south"), nb::arg("east"), nb::arg("north"));
    module.def("opentopography_datasets", &opentopography_datasets,
        "Return OpenTopography DEM catalogue records as dictionaries.");
}

} // namespace lsdtt3_python
