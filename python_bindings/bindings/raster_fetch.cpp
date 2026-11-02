#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <cstdint>
#include <optional>

#include <nanobind/stl/optional.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>

#include <gdal_priv.h>
#include <Python.h>

#include <core/lsdtt_gdal_init.hpp>
#include <core/lsdtt_secret.hpp>
#include <core/lsdtt_dem_grouping.hpp>
#include <cpu/lsdtt_outlet_basin_verification_cpu.hpp>

#include "binding_common.hpp"
#include "rotation_frames.hpp"

namespace lsdtt3_python {
namespace {

lsdtt::lsdtt_dem_source source_from_name(const std::string& source)
{
    lsdtt::lsdtt_dem_source parsed;
    if (!lsdtt::parse_dem_source(source, parsed))
    {
        throw std::invalid_argument(
            "source must be one of " + lsdtt::dem_source_name_list()
        );
    }
    return parsed;
}

// The catalogue sources (national and polar DEMs), as opposed to the frozen
// Copernicus / OpenTopography / Planetary Computer paths.
bool is_catalogue_source(lsdtt::lsdtt_dem_source source)
{
    return !lsdtt::dem_source_info(source).legacy_path;
}

void warn(const std::string& text)
{
    if (PyErr_WarnEx(PyExc_UserWarning, text.c_str(), 2) < 0) throw nb::python_error();
}

// The limits lsdtt-fetch-raster applies to every source.
void set_fetch_limits(lsdtt::lsdtt_fetch_options& options, std::int64_t max_output_pixels,
                      double max_download_megabytes, int max_tiles,
                      int request_delay_milliseconds)
{
    options.max_output_pixels = max_output_pixels;
    options.max_download_megabytes = max_download_megabytes;
    options.max_tiles = max_tiles;
    options.request_delay_milliseconds = request_delay_milliseconds;
}

Raster fetch_dem(
    double west, double south, double east, double north,
    const std::string& source, int target_epsg,
    std::optional<double> grid_spacing_given, const std::string& resampling, double nodata_value,
    int size_pixels, const std::string& cache_dir, const std::string& ot_dem_type,
    const std::string& pc_collection, const std::string& pc_asset,
    int http_timeout_seconds, int http_connect_timeout_seconds, int http_max_retry,
    nb::object frame, bool fetch_from_linear_feature, const std::string& fetch_feature_fname,
    const std::string& fetch_feature_layer, double linear_feature_fetch_width,
    double linear_feature_fetch_end_buffer, bool rotate_after_fetch,
    double source_resolution_metres, const std::string& surface_model_choice,
    const std::string& source_collection_choice, std::int64_t max_output_pixels,
    double max_download_megabytes, int max_tiles, int request_delay_milliseconds)
{
    lsdtt::gdal_init();
    const lsdtt::lsdtt_dem_source parsed_source = source_from_name(source);
    // grid_spacing not given: the frame's cell size, a catalogue source's
    // native resolution (source_resolution_metres if set), as
    // lsdtt-fetch-raster does, or 30 m.
    double grid_spacing = 30.0;
    if (grid_spacing_given) grid_spacing = *grid_spacing_given;
    else if (!frame.is_none())
        grid_spacing = nb::cast<const RotationFrame&>(frame).native().frame.cell_size;
    else if (is_catalogue_source(parsed_source))
        grid_spacing = source_resolution_metres > 0.0
                     ? source_resolution_metres
                     : lsdtt::dem_source_info(parsed_source).default_resolution_metres;
    lsdtt::lsdtt_fetch_options options;
    options.west = west;
    options.south = south;
    options.east = east;
    options.north = north;

    // The area from a frame or a feature file (lsdtt-fetch-raster's
    // fetch_from_linear_feature), refusing every other way of giving one.
    const bool have_frame = !frame.is_none();
    const bool any_edge = std::isfinite(west) || std::isfinite(south) ||
                          std::isfinite(east) || std::isfinite(north);
    if (have_frame || fetch_from_linear_feature)
    {
        const char* what = have_frame ? "frame" : "fetch_from_linear_feature";
        if (have_frame && fetch_from_linear_feature)
            throw std::invalid_argument("frame and fetch_from_linear_feature cannot be combined");
        if (any_edge)
            throw std::invalid_argument(std::string(what) + " defines the area itself, so it "
                                        "cannot be combined with west/south/east/north");
        if (size_pixels != 0)
            throw std::invalid_argument(std::string(what) + " cannot be combined with size_pixels");
    }
    else if (!std::isfinite(west) || !std::isfinite(south) ||
             !std::isfinite(east) || !std::isfinite(north))
    {
        throw std::invalid_argument(
            "give west, south, east and north, or frame, or fetch_from_linear_feature");
    }
    if (!fetch_from_linear_feature && !fetch_feature_fname.empty())
    {
        if (PyErr_WarnEx(PyExc_UserWarning,
                         "fetch_feature_fname is set but fetch_from_linear_feature is False, "
                         "so it is ignored", 2) < 0)
            throw nb::python_error();
    }
    if (have_frame)
    {
        const RotationFrame& rf = nb::cast<const RotationFrame&>(frame);
        if (target_epsg != 0)
            throw std::invalid_argument(
                "target_epsg cannot be given with frame: the frame carries its rotate-back "
                "target (RotationFrame.from_feature(..., target_epsg=...))");
        if (std::abs(rf.native().frame.cell_size - grid_spacing) > 1e-9 * grid_spacing)
            throw std::invalid_argument(
                "grid_spacing (" + std::to_string(grid_spacing) + ") must equal the frame's "
                "cell_size (" + std::to_string(rf.native().frame.cell_size) + ")");
        options.use_frame = true;
        options.frame = rf.native().frame;
        options.target_epsg = rf.native().target.epsg;
    }
    else if (fetch_from_linear_feature)
    {
        if (fetch_feature_fname.empty())
            throw std::invalid_argument(
                "fetch_from_linear_feature needs fetch_feature_fname (a line, two points or "
                "polygon file)");
        const lsdtt::lsdtt_frame_from_feature_result fr = lsdtt::frame_from_feature_file(
            fetch_feature_fname, fetch_feature_layer, linear_feature_fetch_width,
            linear_feature_fetch_end_buffer, grid_spacing);
        for (const auto& w : fr.warnings)
            if (PyErr_WarnEx(PyExc_UserWarning, w.c_str(), 2) < 0) throw nb::python_error();
        if (!fr.ok) throw std::invalid_argument(fr.message);
        if (rotate_after_fetch)
        {
            options.use_frame = true;
            options.frame = fr.frame;
        }
        else
        {
            lsdtt::fetch_bbox_for_frame(fr.frame, 0.0, options.west, options.south,
                                        options.east, options.north);
        }
    }
    options.source = parsed_source;
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
    options.source_resolution_metres = source_resolution_metres;
    options.surface_model_choice = surface_model_choice;
    options.source_collection_choice = source_collection_choice;
    set_fetch_limits(options, max_output_pixels, max_download_megabytes, max_tiles,
                     request_delay_milliseconds);
    const std::string lerc_warning = lsdtt::lerc_source_warning(options.source);
    if (!lerc_warning.empty()) warn(lerc_warning);

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
    for (const std::string& w : result.warnings) warn(w);

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
    for (const std::string& w : result.warnings) warn(w);
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
    int min_search_radius_pixels, int max_search_radius_pixels, bool mask_outside_basins,
    std::int64_t max_output_pixels, double max_download_megabytes, int max_tiles,
    int request_delay_milliseconds)
{
    lsdtt::gdal_init();
    const lsdtt::lsdtt_dem_source parsed_source = source_from_name(source);
    // As lsdtt-fetch-raster's bounds_mode: outlets.
    if (is_catalogue_source(parsed_source))
        throw std::invalid_argument(
            "source " + source + " cannot be used with fetch_dem_for_outlets in this "
            "version: fetch the outlet groups with cop30_aws, then fetch each "
            "group's area with fetch_dem(..., source=\"" + source + "\")");
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
        set_fetch_limits(options, max_output_pixels, max_download_megabytes, max_tiles,
                         request_delay_milliseconds);
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

nb::list dem_sources()
{
    nb::list result;
    for (const auto& info : lsdtt::dem_source_catalogue())
    {
        nb::dict item;
        item["dem_source"] = info.name;
        item["title"] = info.title;
        item["provider"] = info.provider;
        item["description"] = info.description;
        item["coverage"] = info.coverage_text;
        item["native_resolutions_metres"] = info.native_resolutions_metres;
        item["default_resolution_metres"] = info.default_resolution_metres;
        item["surface_models"] = info.surface_models;
        item["default_surface_model"] = info.default_surface_model;
        item["vertical_datum"] = info.vertical_datum;
        nb::list collections;
        for (const auto& c : info.grid_collections) collections.append(c.label);
        item["collections"] = std::move(collections);
        item["attribution"] = info.licence_attribution;
        item["accessed"] = info.accessed;
        item["requires_lerc"] = info.requires_lerc;
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
        .def_prop_ro("outlets", &outlet_table, nb::sig("def outlets(self) -> pandas.DataFrame"),
            "Initial outlet envelopes: outlet_id, longitude, latitude, status, iterations, west, south, east, north, message.")
        .def_prop_ro("snapped_outlets", &snapped_outlet_table, nb::sig("def snapped_outlets(self) -> pandas.DataFrame"),
            "Snapped outlets: group_id, outlet_id, row, col, x, y (pixel centre in the group raster's CRS).")
        .def_prop_ro("groups", &group_table, nb::sig("def groups(self) -> pandas.DataFrame"),
            "One row per DEM: group_id, west, south, east, north (WGS84), area_km2, outlet_ids, basins_verified, verification_message.");

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

    constexpr double no_edge = std::numeric_limits<double>::quiet_NaN();
    module.def("fetch_dem", &fetch_dem,
        nb::arg("west") = no_edge, nb::arg("south") = no_edge,
        nb::arg("east") = no_edge, nb::arg("north") = no_edge,
        nb::arg("source") = "cop30_aws",
        nb::arg("target_epsg") = 0, nb::arg("grid_spacing") = nb::none(),
        nb::arg("resampling") = "cubic", nb::arg("nodata_value") = -9999.0,
        nb::arg("size_pixels") = 0, nb::arg("cache_dir") = "",
        nb::arg("ot_dem_type") = "COP30", nb::arg("pc_collection") = "cop-dem-glo-30",
        nb::arg("pc_asset") = "data", nb::arg("http_timeout_seconds") = 60,
        nb::arg("http_connect_timeout_seconds") = 15, nb::arg("http_max_retry") = 3,
        nb::arg("frame") = nb::none(), nb::arg("fetch_from_linear_feature") = false,
        nb::arg("fetch_feature_fname") = "", nb::arg("fetch_feature_layer") = "",
        nb::arg("linear_feature_fetch_width") = 0.0,
        nb::arg("linear_feature_fetch_end_buffer") = -1.0,
        nb::arg("rotate_after_fetch") = true,
        nb::arg("source_resolution_metres") = 0.0,
        nb::arg("surface_model_choice") = "source_default",
        nb::arg("source_collection_choice") = "latest",
        nb::arg("max_output_pixels") = 250000000,
        nb::arg("max_download_megabytes") = 2048.0,
        nb::arg("max_tiles") = 256, nb::arg("request_delay_milliseconds") = 200,
        R"doc(
Fetch a DEM for a WGS84 bounding box as an in-memory projected raster.

The source tiles are downloaded, mosaicked and warped to a projected grid
(UTM by default), with the same code and the same parameter names as
``lsdtt-fetch-raster``. The Copernicus DEM (the default ``cop30_aws``
source) stores the sea as exactly 0 m, not as nodata, so mask it before
flow routing: ``dem = mask_to_nodata(dem, 0.0)``. Cells outside the source
footprint are nodata.

Besides the global 30 m sources there are keyless national and polar
high-resolution DEMs (``us_3dep``, ``arctic_arcticdem``,
``antarctic_rema``, ``no_dtm``, ``nl_ahn``, ``es_mdt``, ``pl_nmt``,
``ca_hrdem``, ``ch_swissalti3d``, ``nz_lidar``, ``gb_sct_lidar``,
``jp_gsi``); ``dem_sources()`` lists each one's resolutions, surface models
(DSM or DTM), vertical datum, coverage and attribution. Warnings from the
fetch (partial coverage, overviews read, a LERC-compressed source this GDAL
cannot read, older Scottish collections covering the area) are issued as
``UserWarning``.

Or, with ``frame`` (a ``RotationFrame``) or ``fetch_from_linear_feature``,
fetch a strip along a fault, river or range warped straight from the source
tiles into an oblique frame (one resampling pass), with the feature along
+x. The return type is always ``Raster``; keep the ``RotationFrame`` to
rotate products back with ``rotate_back``.

Parameters
----------
west, south, east, north : float, optional
    WGS84 bounds in decimal degrees (longitude west/east, latitude
    south/north). Antimeridian-crossing requests are not supported:
    ``west`` must be less than ``east``. ``VectorBounds.bbox`` from
    ``bounds_from_vector`` can be unpacked here: ``fetch_dem(*b.bbox)``.
    Required unless ``frame`` or ``fetch_from_linear_feature`` gives the
    area; refused with them.
source : str, optional
    ``cop30_aws`` (default) is a keyless global 30 m Copernicus DEM.
    OpenTopography reads its key only from the ``OT_API_KEY`` or
    ``OT_API_KEY_FILE`` environment variables (never pass a key as an
    argument). ``planetary_computer`` is keyless. The national and polar
    sources are listed by ``dem_sources()``.
target_epsg : int, optional
    Output EPSG code. Zero selects the UTM zone of the bounding-box centre.
grid_spacing : float or None, optional
    Output cell size in target CRS units (m for UTM). None (default) uses
    the frame's cell size with ``frame``, the source's native resolution
    (``source_resolution_metres`` if set) for a national or polar source,
    and 30 otherwise. Output cells are always square.
resampling : str, optional
    GDAL warp resampling method (``"cubic"`` default, ``"bilinear"``,
    ``"near"``, ...).
nodata_value : float, optional
    Nodata value of the output (default -9999).
size_pixels : int, optional
    If > 0, crop an exact ``size_pixels`` x ``size_pixels`` square from the
    centre of the warped grid; 0 keeps the full extent.
cache_dir : str, optional
    Persistent cache directory. This is the only optional disk output; the
    returned raster itself is kept in memory. Empty means no caching.
ot_dem_type : str, optional
    OpenTopography dataset, for example ``"COP30"``, ``"SRTMGL1"``,
    ``"NASADEM"`` (see ``opentopography_datasets()``); other sources ignore it.
pc_collection, pc_asset : str, optional
    Planetary Computer STAC collection and asset key (defaults
    ``"cop-dem-glo-30"`` and ``"data"``); other sources ignore them.
http_timeout_seconds : int, optional
    Total time limit per HTTP request in seconds.
http_connect_timeout_seconds : int, optional
    Connection time limit per HTTP request in seconds.
http_max_retry : int, optional
    Number of retries of a failed HTTP request.
frame : RotationFrame or None, optional
    Fetch the frame's rectangle, warped into the frame. ``grid_spacing``,
    if given, must equal ``frame.cell_size``; ``target_epsg`` and ``size_pixels``
    are refused (the frame carries its rotate-back target).
fetch_from_linear_feature : bool, optional
    Build the frame from ``fetch_feature_fname`` here, as lsdtt-fetch-raster
    does; the same as ``frame=RotationFrame.from_feature(fetch_feature_fname,
    linear_feature_fetch_width, linear_feature_fetch_end_buffer,
    grid_spacing, fetch_feature_layer)``, which is the frame to keep for
    ``rotate_back``.
fetch_feature_fname : str, optional
    Line, two-point or polygon file (any OGR format and CRS).
fetch_feature_layer : str, optional
    Layer of the feature file; empty reads the first.
linear_feature_fetch_width : float, optional
    Full strip width in metres; required for a line or two points.
linear_feature_fetch_end_buffer : float, optional
    Length beyond each axis end in metres (all sides for a polygon);
    negative means 10 cells.
rotate_after_fetch : bool, optional
    With ``fetch_from_linear_feature``: True (default) returns the strip in
    the oblique frame; False returns the frame's WGS84 envelope in UTM, as a
    plain fetch.
source_resolution_metres : float, optional
    Native resolution (m) to fetch from a source that offers several (for
    example 1, 2 or 10 for ``arctic_arcticdem``); 0 (default) uses the
    source's default. Refused for a source that does not offer it.
surface_model_choice : {"source_default", "dtm", "dsm"}, optional
    Terrain (DTM) or surface (DSM) model. ``source_default`` takes the DTM
    where the source has one. Asking for a model the source lacks raises
    ``RuntimeError``; the Copernicus sources are DSMs only.
source_collection_choice : str, optional
    ``gb_sct_lidar`` only: ``latest`` (default) uses the newest collection
    with tiles in the area (a warning names older ones that also cover it);
    a collection label (``phase-1`` ... ``phase-6``,
    ``national-lidar-programme``) uses that one. Collections are never mixed.
max_output_pixels : int, optional
    Largest output grid (cells) allowed, checked before any download; for
    the WCS sources the native cells requested are counted too.
max_download_megabytes : float, optional
    Largest download allowed (MB): an estimate checked before the fetch,
    and a running count that stops the transfer where lsdtt3 does the HTTP.
max_tiles : int, optional
    Largest number of tiles, source files or service requests allowed.
request_delay_milliseconds : int, optional
    Pause between requests to tile services (ms).

Returns
-------
Raster
    A projected raster of elevation (m): UTM (or ``target_epsg``), or the
    oblique frame. Network, GDAL, and invalid-request failures (including
    an area outside the source's coverage or over a size limit) raise
    ``RuntimeError``; an inconsistent combination of arguments raises
    ``ValueError``.
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
        nb::arg("max_output_pixels") = 250000000,
        nb::arg("max_download_megabytes") = 2048.0,
        nb::arg("max_tiles") = 256, nb::arg("request_delay_milliseconds") = 200,
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
source : {"cop30_aws", "opentopography", "planetary_computer"}, optional
    DEM source, as in ``fetch_dem``. The national and polar sources are
    refused here, as in ``lsdtt-fetch-raster`` outlets mode: find the groups
    with ``cop30_aws``, then ``fetch_dem`` each group's area.
target_epsg : int, optional
    Output EPSG code; 0 selects UTM from each group's centre.
grid_spacing : float, optional
    Output cell size in target CRS units (m for UTM).
initial_radius_km : float, optional
    Initial half-width (km) of the box around each outlet before grouping.
auto_split_disconnected_basins : bool, optional
    Split groups whose outlets form separate clusters.
split_cluster_gap_km : float, optional
    Gap (km) between outlet clusters above which a group is split.
max_final_dem_area_km2 : float, optional
    Largest allowed grouped initial extent (km^2) before splitting it.
split_overlap_km : float, optional
    Overlap (km) between split extents.
max_basins_per_dem : int, optional
    Largest number of outlets per DEM; 0 means no limit.
cache_dir : str, optional
    Persistent download cache directory; empty means none.
resampling : str, optional
    GDAL warp resampling method.
nodata_value : float, optional
    Nodata value of the output DEMs.
snapping_method : {"none", "expanding_square", "downstream", "largest_area"}, optional
    Outlet snapping used before each basin delineation; the same rules as
    ``FlowInfo.snap_points`` and the CLI ``basin_outlet_snapping_method``.
contributing_pixels_threshold : int, optional
    Channel threshold (contributing pixels) that outlets are snapped to.
edge_buffer_pixels, nodata_buffer_pixels : int, optional
    Required basin clearance (pixels) from the DEM edge and invalid cells.
tight_bbox_buffer_pixels : int, optional
    Margin (pixels) added around the verified basins for the tight refetch.
verify_expand_factor : float, optional
    Factor (> 1) by which an extent grows when a basin reaches its edge.
verify_max_iters : int, optional
    Largest number of fetch-and-verify attempts per group.
tight_bbox_refetch : bool, optional
    Refetch a smaller final DEM around the verified basin footprint.
min_search_radius_pixels : int, optional
    Smallest expanding-square snap radius (pixels).
max_search_radius_pixels : int, optional
    Largest expanding-square snap radius when ``snapping_method`` is
    ``"expanding_square"``, and the half-width of the square searched for
    ``"largest_area"``.
mask_outside_basins : bool, optional
    Set cells outside the selected basin union to nodata in returned DEMs.
max_output_pixels, max_download_megabytes, max_tiles, request_delay_milliseconds
    Size limits and request pacing for each group's fetch, as in
    ``fetch_dem``.

Returns
-------
OutletFetchResult
    Grouped raster collection and copied outlet/group diagnostic tables.
)doc");

    module.def("bounds_from_vector", [](const std::string& path, const std::string& layer) {
        lsdtt::gdal_init();
        return lsdtt::bounds_from_vector_file(path, layer);
    }, nb::arg("path"), nb::arg("layer") = "", R"doc(
Derive WGS84 ``(west, south, east, north)`` bounds from a vector file.

The layer must contain polygons or exactly two points (opposite corners of
a box). Geometry is transformed to WGS84 first. ``bbox`` is a tuple
``(west, south, east, north)``, so ``fetch_dem(*b.bbox)`` fetches the area.

Parameters
----------
path : str
    Vector file readable by GDAL/OGR (GeoPackage, shapefile, GeoJSON, ...).
layer : str, default=""
    Layer name; empty uses the first layer.

Returns
-------
VectorBounds
    Bounds in decimal degrees plus diagnostics (``ok``, ``message``).
)doc");

    module.def("http_available", &lsdtt::http_available,
        R"doc(Return whether the bundled GDAL has HTTP support.

Returns
-------
bool
    True if remote downloads (``fetch_dem``) can work.
)doc");
    module.def("utm_epsg_for_point", &lsdtt::utm_epsg_for_point,
        nb::arg("longitude"), nb::arg("latitude"), R"doc(
Return the UTM EPSG code containing a WGS84 point.

Parameters
----------
longitude, latitude : float
    WGS84 position in decimal degrees (longitude first).

Returns
-------
int
    EPSG code of the WGS84 UTM zone: 326zz in the northern hemisphere,
    327zz in the southern (for example 32633 for UTM 33N).
)doc");
    module.def("cop30_tiles_for_bbox", &lsdtt::cop30_tiles_for_bbox,
        nb::arg("west"), nb::arg("south"), nb::arg("east"), nb::arg("north"), R"doc(
Return the Copernicus 30 m DEM tile names covering a WGS84 box.

Parameters
----------
west, south, east, north : float
    WGS84 bounds in decimal degrees.

Returns
-------
list of str
    Tile identifiers of the 1 x 1 degree tiles intersecting the box.
)doc");
    module.def("dem_sources", &dem_sources,
        R"doc(Return the ``fetch_dem`` DEM source catalogue as dictionaries.

The same rows ``lsdtt-fetch-raster -h dem_source`` describes.

Returns
-------
list of dict
    One dict per source with keys ``dem_source`` (the ``source`` value for
    ``fetch_dem``), ``title``, ``provider``, ``description``, ``coverage``,
    ``native_resolutions_metres``, ``default_resolution_metres``,
    ``surface_models`` (``"dtm"`` and/or ``"dsm"``),
    ``default_surface_model``, ``vertical_datum``, ``collections`` (the
    ``source_collection_choice`` labels; empty for most sources),
    ``attribution`` (to cite when you use the data), ``accessed`` (when the
    access details were last checked) and ``requires_lerc``.
)doc");
    module.def("opentopography_datasets", &opentopography_datasets,
        R"doc(Return OpenTopography DEM catalogue records as dictionaries.

Returns
-------
list of dict
    One dict per dataset with keys ``demtype`` (the ``ot_dem_type`` value
    for ``fetch_dem``), ``name``, ``resolution``, ``coverage``,
    ``description``, ``citation`` and ``max_area_km2``.
)doc");
}

} // namespace lsdtt3_python
