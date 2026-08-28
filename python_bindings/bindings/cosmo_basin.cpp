#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>

#include <cpu/lsdtt_cosmo_raster_cpu.hpp>
#include <cpu/lsdtt_flow_snapping_cpu.hpp>

#include "binding_common.hpp"

namespace lsdtt3_python {
namespace {

void require_flow_grid(const FlowInfo& flow, const Raster& raster, const char* name)
{
    const auto& a = flow.native().meta();
    const auto& b = raster.native().get_meta();
    if (a.n_rows != b.n_rows || a.n_cols != b.n_cols || a.cell_size != b.cell_size ||
        a.x_min != b.x_min || a.y_min != b.y_min)
        throw std::invalid_argument(std::string(name) + " must share the FlowInfo grid");
}

const lsdtt::lsdtt_raster_cpu& optional_or_uniform(const FlowInfo& flow, const Raster* supplied,
                                                    float value, lsdtt::lsdtt_raster_cpu& storage,
                                                    const char* name)
{
    if (supplied) {
        require_flow_grid(flow, *supplied, name);
        return supplied->native();
    }
    std::vector<float> data(flow.native().full_grid_size(), value);
    std::vector<std::uint8_t> mask(flow.native().mask_array().data(),
                                   flow.native().mask_array().data() + flow.native().full_grid_size());
    storage = lsdtt::lsdtt_raster_cpu(std::move(data), std::move(mask), flow.native().meta());
    return storage;
}

double x_from_col(const lsdtt::lsdtt_raster_meta& meta, int col)
{ return meta.x_min + (static_cast<double>(col) + 0.5) * meta.cell_size; }

double y_from_row(const lsdtt::lsdtt_raster_meta& meta, int row)
{ return meta.y_max() - (static_cast<double>(row) + 0.5) * meta.cell_size; }

nb::object crn_invert_basin_erosion(
    const FlowInfo& flow_owner, nb::object samples, const Raster& production_scaling,
    const std::string& concentration_column, const std::string& nuclide,
    const std::string& muon_scheme, int channel_threshold_pixels, int snap_radius_pixels,
    const Raster* topographic_shielding, const Raster* snow_shielding,
    const Raster* self_shielding, double rock_density_kg_m3)
{
    if (channel_threshold_pixels < 1 || snap_radius_pixels < 0 || rock_density_kg_m3 <= 0.0)
        throw std::invalid_argument("channel_threshold_pixels must be positive, snap_radius_pixels non-negative, and rock_density_kg_m3 positive");
    require_flow_grid(flow_owner, production_scaling, "production_scaling");

    nb::object geopandas = nb::module_::import_("geopandas");
    if (!nb::isinstance(samples, geopandas.attr("GeoDataFrame")))
        throw std::invalid_argument("samples must be a GeoDataFrame of Point geometries");
    nb::object result = samples.attr("copy")();
    const auto& meta = flow_owner.native().meta();
    if (!result.attr("crs").is_none() && !meta.crs_wkt.empty()) result = result.attr("to_crs")(meta.crs_wkt);
    if (result.attr("crs").is_none() && !meta.crs_wkt.empty())
        throw std::invalid_argument("samples must define a CRS when the flow grid has a CRS");
    if (!nb::cast<bool>(result.attr("columns").attr("__contains__")(concentration_column)))
        throw std::invalid_argument("samples is missing concentration column '" + concentration_column + "'");

    const std::size_t n = nb::cast<std::size_t>(result.attr("shape")[0]);
    nb::object geometry = result.attr("geometry");
    nb::object concentrations = result.attr("__getitem__")(concentration_column);
    lsdtt::lsdtt_raster_cpu topo, snow, self;
    const auto& topo_raster = optional_or_uniform(flow_owner, topographic_shielding, 1.0f, topo, "topographic_shielding");
    const auto& snow_raster = optional_or_uniform(flow_owner, snow_shielding, 0.0f, snow, "snow_shielding");
    const auto& self_raster = optional_or_uniform(flow_owner, self_shielding, 0.0f, self, "self_shielding");
    lsdtt::lsdtt_square_snap_options snap_options;
    snap_options.contributing_pixels_threshold = static_cast<std::uint64_t>(channel_threshold_pixels);
    snap_options.max_search_radius_pixels = snap_radius_pixels;
    snap_options.min_search_radius_pixels = 0;

    const double nan = std::numeric_limits<double>::quiet_NaN();
    std::vector<std::string> status(n);
    std::vector<int> input_row(n, -1), input_col(n, -1), snapped_row(n, -1), snapped_col(n, -1), search_radius(n, -1);
    std::vector<std::int64_t> outlet_flat(n, -1);
    std::vector<std::uint64_t> contributing_pixels(n, 0), basin_pixels(n, 0);
    std::vector<double> input_x(n, nan), input_y(n, nan), snapped_x(n, nan), snapped_y(n, nan);
    std::vector<double> erosion_rate(n, nan), erosion_rate_mm_kyr(n, nan), predicted_concentration(n, nan), basin_area_m2(n, nan);

    const auto& flow = flow_owner.native();
    for (std::size_t i = 0; i < n; ++i) {
        nb::object point = geometry[static_cast<int>(i)];
        if (point.is_none() || nb::cast<std::string>(point.attr("geom_type")) != "Point") { status[i] = "invalid_geometry"; continue; }
        input_x[i] = nb::cast<double>(point.attr("x")); input_y[i] = nb::cast<double>(point.attr("y"));
        const double measured = nb::cast<double>(concentrations[static_cast<int>(i)]);
        if (!std::isfinite(measured) || measured <= 0.0) { status[i] = "invalid_concentration"; continue; }
        input_row[i] = lsdtt::row_from_y(meta, input_y[i]); input_col[i] = lsdtt::col_from_x(meta, input_x[i]);
        if (!flow.in_bounds(input_row[i], input_col[i]) || !flow.is_valid_cell(input_row[i], input_col[i])) { status[i] = "outside_flow_grid"; continue; }
        const auto snapped = lsdtt::snap_pixel_to_contributing_pixels_threshold_expanding_square(flow, input_row[i], input_col[i], snap_options);
        if (snapped.status != lsdtt::lsdtt_snap_status::success) { status[i] = "channel_not_found"; continue; }
        snapped_row[i] = snapped.snapped_row; snapped_col[i] = snapped.snapped_col; outlet_flat[i] = snapped.snapped_flat;
        contributing_pixels[i] = snapped.snapped_contributing_pixels; search_radius[i] = snapped.search_radius_pixels;
        snapped_x[i] = x_from_col(meta, snapped.snapped_col); snapped_y[i] = y_from_row(meta, snapped.snapped_row);
        const auto upslope = flow.upslope_nodes(snapped.snapped_flat);
        basin_pixels[i] = static_cast<std::uint64_t>(upslope.size());
        basin_area_m2[i] = static_cast<double>(upslope.size()) * meta.cell_size * meta.cell_size;
        erosion_rate[i] = lsdtt::calculate_eff_erate_from_conc(measured, nuclide, muon_scheme,
            production_scaling.native(), topo_raster, snow_raster, self_raster, flow,
            snapped.snapped_flat, predicted_concentration[i]);
        erosion_rate_mm_kyr[i] = erosion_rate[i] * 1.0e7 / rock_density_kg_m3;
        status[i] = "ok";
    }

    const auto assign = [&result](const char* name, const auto& values) { result.attr("__setitem__")(name, nb::cast(values)); };
    assign("crn_status", status); assign("input_x", input_x); assign("input_y", input_y); assign("input_row", input_row); assign("input_col", input_col);
    assign("snapped_row", snapped_row); assign("snapped_col", snapped_col); assign("snapped_x", snapped_x); assign("snapped_y", snapped_y); assign("outlet_flat", outlet_flat);
    assign("channel_contributing_pixels", contributing_pixels); assign("snap_radius_pixels", search_radius); assign("basin_pixels", basin_pixels); assign("basin_area_m2", basin_area_m2);
    assign("erosion_rate_g_cm2_yr", erosion_rate); assign("erosion_rate_mm_kyr", erosion_rate_mm_kyr); assign("predicted_concentration_atoms_g", predicted_concentration);
    return result;
}

} // namespace

void bind_cosmo_basin(nb::module_& module)
{
    module.def("crn_invert_basin_erosion", &crn_invert_basin_erosion,
        nb::arg("flow"), nb::arg("samples"), nb::arg("production_scaling"), nb::arg("concentration_column") = "concentration",
        nb::arg("nuclide") = "Be10", nb::arg("muon_scheme") = "BraucherBorchers", nb::arg("channel_threshold_pixels") = 1000,
        nb::arg("snap_radius_pixels") = 25, nb::arg("topographic_shielding") = nb::none(), nb::arg("snow_shielding") = nb::none(),
        nb::arg("self_shielding") = nb::none(), nb::arg("rock_density_kg_m3") = 2650.0,
        R"doc(Invert basin-mean erosion rates from a point GeoDataFrame.

Points are projected to the flow-grid CRS, snapped to the nearest D8 channel
meeting ``channel_threshold_pixels``, then inverted with the native CAIRN
solver over each snapped outlet's exact upslope pixel set. Failed rows remain
in the returned GeoDataFrame and are identified by ``crn_status``.)doc");
}

} // namespace lsdtt3_python
