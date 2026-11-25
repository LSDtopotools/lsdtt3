// hillslope_channel_tagging.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
//
// -----------------------------------------------------------------------------
// lsdtt-hillslope-channel-tagging's routines: the same C++ calls as the CLI,
// with its parameter names less the hillslope_ and channel_attribute_ prefixes.
// A DataFrame of channel attributes is written to a temporary file and read by
// the CLI's reader, so that the two cannot join differently.
// -----------------------------------------------------------------------------

#include <cmath>
#include <cstdint>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>

#include <core/lsdtt_channel_attribute_input.hpp>
#include <core/lsdtt_vector_utils.hpp>
#include <cpu/lsdtt_hillslope_channel_tagging_cpu.hpp>

#include "binding_common.hpp"

namespace lsdtt3_python {
namespace {

using meta_t = lsdtt::lsdtt_raster_meta;
using i64_array = nb::ndarray<const std::int64_t, nb::ndim<1>, nb::c_contig>;

void require_grid(const meta_t& a, const meta_t& b, const std::string& what)
{
    if (a.n_rows != b.n_rows || a.n_cols != b.n_cols || a.x_min != b.x_min || a.y_min != b.y_min ||
        a.cell_size != b.cell_size)
        throw std::invalid_argument(what + " must share the grid");
}

std::vector<std::uint8_t> channel_mask_of(const Raster& channels)
{
    const auto& d = channels.native().get_data();
    const auto& m = channels.native().get_mask();
    std::vector<std::uint8_t> out(d.size(), 0u);
    for (std::size_t i = 0; i < d.size(); ++i)
        out[i] = (m[i] && std::isfinite(d[i]) && d[i] != 0.0f) ? 1u : 0u;
    return out;
}

Raster raster_of(const lsdtt::lsdtt_raster_result& r)
{
    return Raster(lsdtt::lsdtt_raster_cpu(std::vector<float>(r.data), std::vector<std::uint8_t>(r.mask), r.meta));
}

nb::object numpy_of(std::vector<std::int64_t> values)
{
    auto* owned = new std::vector<std::int64_t>(std::move(values));
    nb::capsule owner(owned, [](void* p) noexcept { delete static_cast<std::vector<std::int64_t>*>(p); });
    const std::size_t shape[1] = {owned->size()};
    return nb::cast(nb::ndarray<nb::numpy, std::int64_t, nb::ndim<1>>(owned->data(), 1, shape, owner));
}

std::vector<std::int64_t> flat_vector(nb::object channel_flat, std::size_t n)
{
    nb::object flat = nb::module_::import_("numpy").attr("asarray")(channel_flat, "int64").attr("ravel")();
    const i64_array a = nb::cast<i64_array>(flat);
    if (a.size() != n)
        throw std::invalid_argument("channel_flat must have one value per grid cell");
    return std::vector<std::int64_t>(a.data(), a.data() + a.size());
}

void add_tag(nb::dict& out, const lsdtt::hillslope_tag_result& t)
{
    out["tag"] = raster_of(t.tag);
    out["n_channel_tagged"] = t.n_channel_tagged;
    out["n_hillslope_tagged"] = t.n_hillslope_tagged;
    out["n_no_value"] = t.n_no_value;
    out["n_excluded"] = t.n_excluded;
    out["n_no_channel"] = t.n_no_channel;
}

nb::dict hillslope_relief(const FlowInfo& owner, const Raster& elevation, const Raster& channels,
                          const std::string& ridge_method, const std::string& relief_grouping,
                          nb::object aspect_surface, const std::string& aspect_method,
                          double polyfit_lambda_metres, int n_threads)
{
    const lsdtt::lsdtt_flowinfo& flow = owner.native();
    require_grid(elevation.native().get_meta(), flow.meta(), "elevation");
    require_grid(channels.native().get_meta(), flow.meta(), "channels");

    lsdtt::hillslope_relief_options options;
    options.ridge = ridge_method;
    options.grouping = relief_grouping;
    options.n_threads = n_threads;
    lsdtt::lsdtt_raster_result aspect;
    if (ridge_method == "aspect_ray")
    {
        if (aspect_method != "polyfit" && aspect_method != "dinf")
            throw std::invalid_argument("aspect_method must be polyfit or dinf");
        const lsdtt::lsdtt_raster_cpu* surface = &elevation.native();
        if (aspect_method == "polyfit")
        {
            if (aspect_surface.is_none())
                throw std::invalid_argument(
                    "ridge_method='aspect_ray' with polyfit needs aspect_surface: pass the input DEM, "
                    "as the CLI does, or the routing surface for polyfit_on_routing_surface");
            const Raster& s = nb::cast<const Raster&>(aspect_surface);
            require_grid(s.native().get_meta(), flow.meta(), "aspect_surface");
            surface = &s.native();
        }
        aspect = lsdtt::compute_hillslope_ray_aspect(*surface, aspect_method, polyfit_lambda_metres,
                                                     flow.boundary_conditions(), n_threads);
        options.ray_aspect = &aspect;
    }

    lsdtt::hillslope_relief_result r = lsdtt::compute_hillslope_relief(
        flow, elevation.native().get_data(), elevation.native().get_mask(), channel_mask_of(channels), options);

    std::vector<float> status(r.ridge_status.size(), static_cast<float>(flow.meta().no_data_value));
    const auto& valid = elevation.native().get_mask();
    for (std::size_t i = 0; i < status.size(); ++i)
        if (valid[i]) status[i] = static_cast<float>(r.ridge_status[i]);

    nb::dict out;
    out["channel_flat"] = numpy_of(std::move(r.channel_flat));
    out["entry_flat"] = numpy_of(std::move(r.entry_flat));
    out["ridge_flat"] = numpy_of(std::move(r.ridge_flat));
    out["ridge_status"] = Raster(lsdtt::lsdtt_raster_cpu(std::move(status), std::vector<std::uint8_t>(valid), flow.meta()));
    out["hand"] = raster_of(r.hand);
    out["ridge_relief"] = raster_of(r.ridge_relief);
    out["relative_height"] = raster_of(r.relative_height);
    out["channel_distance"] = raster_of(r.channel_distance);
    out["ridge_distance"] = raster_of(r.ridge_distance);
    out["n_channel"] = r.n_channel;
    out["n_hillslope"] = r.n_hillslope;
    out["n_no_channel"] = r.n_no_channel;
    out["n_zero_relief"] = r.n_zero_relief;
    out["n_above_ridge"] = r.n_above_ridge;
    out["n_crest_divide"] = r.n_crest_divide;
    out["n_ray_channel"] = r.n_ray_channel;
    out["n_ray_boundary"] = r.n_ray_boundary;
    out["n_no_aspect"] = r.n_no_aspect;
    return out;
}

nb::dict tag_from_channel(nb::object channel_flat, const Raster& channel_values, const Raster& dem)
{
    const meta_t& meta = dem.native().get_meta();
    require_grid(channel_values.native().get_meta(), meta, "channel_values");
    const std::size_t n = dem.native().get_data().size();
    const lsdtt::hillslope_tag_result t = lsdtt::tag_hillslopes_from_channel(
        meta, flat_vector(channel_flat, n), dem.native().get_mask(), channel_values.native().get_data(),
        channel_values.native().get_mask());
    nb::dict out;
    add_tag(out, t);
    return out;
}

nb::dict hillslope_chi_tags(const FlowInfo& owner, const Raster& routing_surface, nb::object channel_flat,
                            nb::object outlet_flats, nb::object chi_max_outlet_elevation_metres,
                            int threshold_contributing_pixels, double m_over_n, double chi_A0,
                            double fit_critical_divergence_metres)
{
    const lsdtt::lsdtt_flowinfo& flow = owner.native();
    require_grid(routing_surface.native().get_meta(), flow.meta(), "routing_surface");
    const std::size_t n = routing_surface.native().get_data().size();

    lsdtt::hillslope_tag_basin_options basin_options;
    basin_options.minimum_contributing_pixels = threshold_contributing_pixels;
    if (!outlet_flats.is_none())
        basin_options.outlet_flats = nb::cast<std::vector<std::int64_t>>(
            nb::module_::import_("numpy").attr("asarray")(outlet_flats, "int64").attr("ravel")().attr("tolist")());
    if (!chi_max_outlet_elevation_metres.is_none())
        basin_options.max_outlet_elevation_m = nb::cast<double>(chi_max_outlet_elevation_metres);

    const lsdtt::hillslope_tag_basins basins =
        lsdtt::select_hillslope_tag_basins(flow, routing_surface.native().get_data(), basin_options);

    lsdtt::hillslope_chi_tag_options options;
    options.chi.m_over_n = m_over_n;
    options.chi.A_0 = chi_A0;
    options.critical_divergence_m = static_cast<float>(fit_critical_divergence_metres);
    options.minimum_contributing_pixels = threshold_contributing_pixels;
    const lsdtt::hillslope_chi_tag_result tags = lsdtt::compute_hillslope_chi_tags(
        flow, routing_surface.native(), flat_vector(channel_flat, n), basins, options);

    nb::dict chi, k_sn, basin_report, out;
    add_tag(chi, tags.chi);
    add_tag(k_sn, tags.k_sn);
    basin_report["outlets"] = numpy_of(basins.outlets);
    basin_report["n_candidates"] = basins.n_candidates;
    basin_report["n_excluded_elevation"] = basins.n_excluded_elevation;
    basin_report["n_excluded_grid_edge"] = basins.n_excluded_grid_edge;
    basin_report["n_excluded_too_small"] = basins.n_excluded_too_small;
    basin_report["n_excluded_nested"] = basins.n_excluded_nested;
    out["chi"] = chi;
    out["k_sn"] = k_sn;
    out["basins"] = basin_report;
    out["n_fitted_nodes"] = tags.n_fitted_nodes;
    out["n_segments"] = tags.n_segments;
    return out;
}

nb::dict channel_attributes_to_raster(const FlowInfo& owner, nb::object source, const Raster& channels,
                                      const std::vector<std::string>& columns, double max_snap_distance_metres,
                                      double max_gap_fill_distance_metres, double nodata_value)
{
    const lsdtt::lsdtt_flowinfo& flow = owner.native();
    require_grid(channels.native().get_meta(), flow.meta(), "channels");

    // A table goes through a temporary file and the CLI's reader.
    lsdtt::channel_attribute_table table;
    const bool is_table = nb::hasattr(source, "to_csv");
    if (is_table)
    {
        nb::object directory = nb::module_::import_("tempfile").attr("TemporaryDirectory")();
        const std::string root = nb::cast<std::string>(directory.attr("name"));
        std::string path;
        if (nb::hasattr(source, "to_file") && nb::hasattr(source, "geometry"))
        {
            path = root + "/channel_attributes.fgb";
            source.attr("to_file")(path, nb::arg("driver") = "FlatGeobuf");
        }
        else
        {
            path = root + "/channel_attributes.csv";
            source.attr("to_csv")(path, nb::arg("index") = false);
        }
        table = lsdtt::read_channel_attribute_table(path, flow.meta(), columns, nodata_value);
        directory.attr("cleanup")();
    }
    else
    {
        const std::string path = nb::cast<std::string>(nb::module_::import_("os").attr("fspath")(source));
        table = lsdtt::read_channel_attribute_table(path, flow.meta(), columns, nodata_value);
    }

    const std::vector<std::uint8_t> all_valid(channels.native().get_data().size(), 1u);
    const lsdtt::channel_attribute_join_result joined = lsdtt::join_channel_attributes(
        flow, all_valid, channel_mask_of(channels), table, max_snap_distance_metres, max_gap_fill_distance_metres);

    nb::dict values, counts, out;
    for (const auto& column : joined.columns)
    {
        values[column.name.c_str()] = Raster(lsdtt::lsdtt_raster_cpu(
            std::vector<float>(column.values), std::vector<std::uint8_t>(column.mask), flow.meta()));
        nb::dict c;
        c["n_direct"] = column.n_direct;
        c["n_multiple_records"] = column.n_multiple_records;
        c["n_filled"] = column.n_filled;
        c["n_missing"] = column.n_missing;
        counts[column.name.c_str()] = c;
    }
    out["values"] = values;
    out["counts"] = counts;
    out["n_channel_pixels"] = joined.n_channel_pixels;
    out["n_records"] = joined.n_records;
    out["n_records_matched"] = joined.n_records_matched;
    out["n_records_dropped"] = joined.n_records_dropped;
    out["n_records_snapped_beyond_one_cell"] = joined.n_records_snapped_beyond_one_cell;
    out["n_exact_samples_off_channel"] = joined.n_exact_samples_off_channel;
    return out;
}

nb::object summarise_by_channel_pixel(nb::object channel_flat, const Raster& values)
{
    const meta_t& meta = values.native().get_meta();
    const std::size_t n = values.native().get_data().size();
    lsdtt::lsdtt_raster_result result;
    result.meta = meta;
    result.data = values.native().get_data();
    result.mask = values.native().get_mask();
    const lsdtt::hillslope_channel_pixel_summary s =
        lsdtt::summarise_hillslopes_by_channel_pixel(flat_vector(channel_flat, n), result);

    const std::size_t m = s.flat.size();
    std::vector<std::int64_t> row(m), col(m);
    std::vector<double> latitude(m), longitude(m);
    const lsdtt::lsdtt_coordinate_transform to_lonlat =
        lsdtt::lsdtt_coordinate_transform::from_epsg(meta.crs_wkt, 4326);
    for (std::size_t k = 0; k < m; ++k)
    {
        row[k] = s.flat[k] / meta.n_cols;
        col[k] = s.flat[k] % meta.n_cols;
        const lsdtt::lsdtt_lonlat ll = lsdtt::cell_center_to_lonlat(
            meta, static_cast<int>(row[k]), static_cast<int>(col[k]), to_lonlat);
        latitude[k] = ll.latitude;
        longitude[k] = ll.longitude;
    }
    nb::dict d;
    d["flat"] = s.flat;
    d["row"] = row;
    d["col"] = col;
    d["latitude"] = latitude;
    d["longitude"] = longitude;
    d["n_hillslope_pixels"] = s.n_hillslope;
    d["n_with_value"] = s.n_with_value;
    d["median"] = s.median;
    return nb::module_::import_("pandas").attr("DataFrame")(d);
}

} // namespace

void bind_hillslope_channel_tagging(nb::module_& module, nb::class_<FlowInfo>& flow)
{
    flow.def("hillslope_relief", &hillslope_relief, nb::arg("elevation"), nb::arg("channels"), nb::kw_only(),
        nb::arg("ridge_method") = "farthest", nb::arg("relief_grouping") = "channel_pixel",
        nb::arg("aspect_surface") = nb::none(), nb::arg("aspect_method") = "polyfit",
        nb::arg("polyfit_lambda_metres") = 70.0, nb::arg("n_threads") = 0,
        nb::sig("def hillslope_relief(self, elevation: Raster, channels: Raster, *, ridge_method: str = 'farthest', "
                "relief_grouping: str = 'channel_pixel', aspect_surface: Raster | None = None, "
                "aspect_method: str = 'polyfit', polyfit_lambda_metres: float = 70.0, n_threads: int = 0) -> dict"),
        R"doc(
HAND, ridge relief and relative height of every hillslope pixel.

lsdtt-hillslope-channel-tagging, the same routine; the keyword names are the
CLI's less the ``hillslope_`` prefix. Each hillslope pixel is linked to the
first channel pixel on its D8 path, whichever ridge method is used, and
relative height is HAND / ridge relief. With ``farthest`` and a shared group
it can exceed 1: it is not clipped, and is counted in ``n_above_ridge``.

Parameters
----------
elevation : Raster
    The routing surface this FlowInfo was built from.
channels : Raster
    Nonzero = channel pixel.
ridge_method : str, default="farthest", keyword-only
    ``farthest``, ``highest`` or ``aspect_ray``.
relief_grouping : str, default="channel_pixel", keyword-only
    Pixels sharing a ridge, for farthest and highest: ``channel_pixel``,
    ``entry_pixel`` or ``pixel``.
aspect_surface : Raster or None, default=None, keyword-only
    Surface for the polyfit aspect of ``aspect_ray``, required there. The CLI
    uses the input DEM by default.
aspect_method : str, default="polyfit", keyword-only
    ``polyfit`` or ``dinf`` (from ``elevation``), for ``aspect_ray``.
polyfit_lambda_metres : float, default=70.0, keyword-only
    Polyfit half-wavelength for the aspect.
n_threads : int, default=0, keyword-only
    OpenMP threads; the result does not depend on it.

Returns
-------
dict
    ``channel_flat``, ``entry_flat``, ``ridge_flat`` (1-D int64 arrays, -1 where
    undefined); Rasters ``ridge_status``, ``hand``, ``ridge_relief``,
    ``relative_height``, ``channel_distance``, ``ridge_distance``; and the
    counts ``n_channel``, ``n_hillslope``, ``n_no_channel``, ``n_zero_relief``,
    ``n_above_ridge``, ``n_crest_divide``, ``n_ray_channel``,
    ``n_ray_boundary``, ``n_no_aspect``.
)doc");

    flow.def("hillslope_chi_tags", &hillslope_chi_tags, nb::arg("routing_surface"), nb::arg("channel_flat"),
        nb::kw_only(), nb::arg("outlet_flats") = nb::none(), nb::arg("chi_max_outlet_elevation_metres") = nb::none(),
        nb::arg("threshold_contributing_pixels") = 1000, nb::arg("m_over_n") = 0.45, nb::arg("chi_A0") = 1.0,
        nb::arg("fit_critical_divergence_metres") = 5.0,
        nb::sig("def hillslope_chi_tags(self, routing_surface: Raster, channel_flat: numpy.ndarray, *, "
                "outlet_flats: numpy.ndarray | None = None, chi_max_outlet_elevation_metres: float | None = None, "
                "threshold_contributing_pixels: int = 1000, m_over_n: float = 0.45, chi_A0: float = 1.0, "
                "fit_critical_divergence_metres: float = 5.0) -> dict"),
        R"doc(
Tag hillslope pixels with the chi and k_sn of their channel pixel.

lsdtt-hillslope-channel-tagging ``tag_chi`` and ``tag_k_sn``, the same routine.
Give exactly one of ``outlet_flats`` and ``chi_max_outlet_elevation_metres``,
because chi needs a true baselevel. Basins with a pixel on the grid edge are
always excluded; their hillslopes are nodata and counted.

Parameters
----------
routing_surface : Raster
    The routing surface this FlowInfo was built from.
channel_flat : numpy.ndarray
    From ``hillslope_relief``.
outlet_flats : numpy.ndarray or None, default=None, keyword-only
    Flat indices of the valid basins' outlets.
chi_max_outlet_elevation_metres : float or None, default=None, keyword-only
    Every baselevel basin with its outlet at or below this elevation.
threshold_contributing_pixels : int, default=1000, keyword-only
    Channel source threshold of the fitted network.
m_over_n : float, default=0.45, keyword-only
    Concavity index.
chi_A0 : float, default=1.0, keyword-only
    Reference drainage area (m^2).
fit_critical_divergence_metres : float, default=5.0, keyword-only
    Critical elevation divergence used to split chi-elevation branches.

Returns
-------
dict
    ``chi`` and ``k_sn``, each a dict with ``tag`` (Raster) and the counts of
    ``tag_hillslopes_from_channel``; ``basins``, a dict with ``outlets`` and the
    numbers excluded; ``n_fitted_nodes`` and ``n_segments``.
)doc");

    flow.def("channel_attributes_to_raster", &channel_attributes_to_raster, nb::arg("source"), nb::arg("channels"),
        nb::arg("columns"), nb::kw_only(), nb::arg("max_snap_distance_metres") = 60.0,
        nb::arg("max_gap_fill_distance_metres") = 150.0, nb::arg("nodata_value") = -9999.0,
        nb::sig("def channel_attributes_to_raster(self, source: str | os.PathLike | pandas.DataFrame, channels: Raster, "
                "columns: list[str], *, max_snap_distance_metres: float = 60.0, "
                "max_gap_fill_distance_metres: float = 150.0, nodata_value: float = -9999.0) -> dict"),
        R"doc(
Join channel attributes from points or lines onto channel pixels.

lsdtt-hillslope-channel-tagging ``channel_attributes_fname``, the same reader
and join; the keyword names are the CLI's less the ``channel_attribute_``
prefix. Records with ``row``, ``col`` and ``flat`` join exactly; others snap to
the nearest channel pixel. Several records on a pixel give their median. A
channel pixel with no record takes the nearest record on its own
junction-to-junction link, and filling never crosses a junction.

Parameters
----------
source : str, os.PathLike, pandas.DataFrame or geopandas.GeoDataFrame
    A CSV or vector file, or a table. A DataFrame needs ``row``, ``col`` and
    ``flat``, or ``latitude`` and ``longitude``. A GeoDataFrame needs a CRS.
channels : Raster
    Nonzero = channel pixel.
columns : list of str
    Numeric columns to join.
max_snap_distance_metres : float, default=60.0, keyword-only
    Furthest a record can be from the channel pixel it joins.
max_gap_fill_distance_metres : float, default=150.0, keyword-only
    Furthest along the channel a value can fill a gap; 0 turns filling off.
nodata_value : float, default=-9999.0, keyword-only
    Attribute value meaning no data. Empty and NaN are always no data.

Returns
-------
dict
    ``values`` (column name to Raster, valid on channel pixels), ``counts``
    (column name to ``n_direct``, ``n_multiple_records``, ``n_filled``,
    ``n_missing``), ``n_channel_pixels``, ``n_records``, ``n_records_matched``,
    ``n_records_dropped``, ``n_records_snapped_beyond_one_cell`` and
    ``n_exact_samples_off_channel``.
)doc");

    module.def("tag_hillslopes_from_channel", &tag_from_channel, nb::arg("channel_flat"), nb::arg("channel_values"),
        nb::arg("dem"),
        nb::sig("def tag_hillslopes_from_channel(channel_flat: numpy.ndarray, channel_values: Raster, dem: Raster) -> dict"),
        R"doc(
Copy each channel pixel's value to the hillslope pixels draining to it.

lsdtt-hillslope-channel-tagging's tag, the same routine. Channel pixels carry
their own value. A pixel whose channel pixel has no value is nodata.

Parameters
----------
channel_flat : numpy.ndarray
    From ``FlowInfo.hillslope_relief``.
channel_values : Raster
    The value at channel pixels; nodata where there is none.
dem : Raster
    Its valid cells are the cells tagged and counted.

Returns
-------
dict
    ``tag`` (Raster), ``n_channel_tagged``, ``n_hillslope_tagged``,
    ``n_no_value``, ``n_excluded`` and ``n_no_channel``.
)doc");

    module.def("summarise_hillslopes_by_channel_pixel", &summarise_by_channel_pixel, nb::arg("channel_flat"),
        nb::arg("values"),
        nb::sig("def summarise_hillslopes_by_channel_pixel(channel_flat: numpy.ndarray, values: Raster) -> pandas.DataFrame"),
        R"doc(
One row per channel pixel: its hillslope pixels and their median value.

lsdtt-hillslope-channel-tagging ``write_channel_pixel_summary``, the same
routine, for any hillslope raster. With ``relative_height`` the CLI's
``n_with_relative_height`` and ``median_relative_height`` are ``n_with_value``
and ``median`` here.

Parameters
----------
channel_flat : numpy.ndarray
    From ``FlowInfo.hillslope_relief``.
values : Raster
    A hillslope raster, for example ``relative_height``.

Returns
-------
pandas.DataFrame
    ``flat``, ``row``, ``col``, ``latitude``, ``longitude``,
    ``n_hillslope_pixels``, ``n_with_value`` and ``median`` (NaN with no value).
)doc");
}

} // namespace lsdtt3_python
