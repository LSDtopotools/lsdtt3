#include <cmath>
#include <cstdint>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <nanobind/nanobind.h>
#include <nanobind/stl/optional.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>

#include <core/lsdtt_common_parameters.hpp>
#include <cpu/lsdtt_flow_snapping_cpu.hpp>

#include "binding_common.hpp"
#include "flow_products.hpp"

namespace lsdtt3_python {
namespace {

std::string validated_boundaries(const std::string& boundaries)
{
    if (boundaries.size() != 4)
        throw std::invalid_argument("boundary_conditions must contain four NESW codes");
    for (const char code : boundaries)
        if (code != 'o' && code != 'n' && code != 'p' && code != 'b')
            throw std::invalid_argument("boundary_conditions uses only 'o', 'n', 'p', or 'b'");
    if ((boundaries[0] == 'p') != (boundaries[2] == 'p') ||
        (boundaries[1] == 'p') != (boundaries[3] == 'p'))
        throw std::invalid_argument("periodic boundaries must occur in opposite pairs");
    return boundaries;
}

std::vector<std::int64_t> forced_outlet_flats(const Raster& dem, nb::object points)
{
    if (points.is_none()) return {};
    nb::object numpy = nb::module_::import_("numpy");
    nb::object array = numpy.attr("asarray")(points, nb::arg("dtype") = "int64");
    if (nb::cast<std::size_t>(array.attr("size")) == 0) return {};
    if (nb::cast<int>(array.attr("ndim")) != 2)
        throw std::invalid_argument("forced_outlets must be an N x 2 array of (row, col)");
    const nb::tuple shape = nb::cast<nb::tuple>(array.attr("shape"));
    if (nb::cast<int>(shape[1]) != 2)
        throw std::invalid_argument("forced_outlets must be an N x 2 array of (row, col)");
    const auto& meta = dem.native().get_meta();
    const auto& mask = dem.native().get_mask();
    const auto values = nb::cast<std::vector<std::int64_t>>(array.attr("ravel")());
    std::vector<std::int64_t> output;
    output.reserve(values.size() / 2);
    for (std::size_t i = 0; i < values.size(); i += 2)
    {
        const auto row = values[i], col = values[i + 1];
        if (row < 0 || row >= meta.n_rows || col < 0 || col >= meta.n_cols)
            throw std::invalid_argument("forced_outlets contains a point outside the raster");
        const auto flat = row * static_cast<std::int64_t>(meta.n_cols) + col;
        if (mask[static_cast<std::size_t>(flat)] == 0)
            throw std::invalid_argument("forced_outlets selects a nodata pixel");
        output.push_back(flat);
    }
    return output;
}

lsdtt::lsdtt_flowinfo make_native_flow_info(const Raster& dem, const std::string& boundary_conditions,
                                            bool interior_nodata_is_outlet, nb::object forced_outlets,
                                            bool allow_pits_as_outlets, int n_threads)
{
    lsdtt::lsdtt_flowinfo_options options;
    options.boundary_conditions = validated_boundaries(boundary_conditions);
    options.interior_nodata_is_outlet = interior_nodata_is_outlet;
    options.forced_outlet_flats = forced_outlet_flats(dem, forced_outlets);
    options.allow_pits_as_outlets = allow_pits_as_outlets;
    options.n_threads = n_threads;
    options.build_topological_order = true;
    options.build_s_vector = true;
    return lsdtt::lsdtt_flowinfo(
        dem.native().get_data_array(), dem.native().get_mask_array(),
        dem.native().get_meta(), options);
}

template <class T>
nb::object grid_array(const T* data, std::size_t size,
                      const lsdtt::lsdtt_raster_meta& meta, const char* dtype)
{
    nb::object numpy = nb::module_::import_("numpy");
    std::vector<T> values(data, data + size);
    return numpy.attr("array")(values, nb::arg("dtype") = dtype)
        .attr("reshape")(nb::make_tuple(meta.n_rows, meta.n_cols));
}

nb::object outlet_pixels(const FlowInfo& flow)
{
    const auto& native = flow.native();
    std::vector<std::int64_t> output;
    output.reserve(native.outlets().size() * 2);
    for (const auto flat : native.outlets())
    {
        output.push_back(native.row_from_flat_unchecked(flat));
        output.push_back(native.col_from_flat_unchecked(flat));
    }
    nb::object numpy = nb::module_::import_("numpy");
    return numpy.attr("array")(output, nb::arg("dtype") = "int64")
        .attr("reshape")(nb::make_tuple(static_cast<int>(native.outlets().size()), 2));
}

nb::object receiver_pixel(const FlowInfo& flow, int row, int col)
{
    const auto receiver = flow.native().receiver_of(flow.native().flat_index(row, col));
    if (receiver == lsdtt::lsdtt_flowinfo::outlet_receiver) return nb::none();
    return nb::make_tuple(flow.native().row_from_flat(receiver), flow.native().col_from_flat(receiver));
}

nb::object upslope_pixels(const FlowInfo& flow, int row, int col)
{
    const auto slice = flow.native().upslope_nodes(row, col);
    std::vector<std::int64_t> output;
    output.reserve(slice.size() * 2);
    for (const auto* it = slice.begin(); it != slice.end(); ++it)
    {
        output.push_back(flow.native().row_from_flat_unchecked(*it));
        output.push_back(flow.native().col_from_flat_unchecked(*it));
    }
    nb::object numpy = nb::module_::import_("numpy");
    return numpy.attr("array")(output, nb::arg("dtype") = "int64")
        .attr("reshape")(nb::make_tuple(static_cast<int>(slice.size()), 2));
}
std::vector<double> coordinate_values(nb::object values, const char* name)
{
    nb::object numpy = nb::module_::import_("numpy");
    nb::object array = numpy.attr("asarray")(values, nb::arg("dtype") = "float64");
    if (nb::cast<int>(array.attr("ndim")) > 1)
        throw std::invalid_argument(std::string("snap_points: ") + name + " must be a number or a 1-D sequence");
    return nb::cast<std::vector<double>>(array.attr("ravel")());
}

// Map coordinates (in the flow grid's CRS) of the points given to snap_points
// by exactly one of x/y, longitude/latitude or a point GeoDataFrame.
void snap_input_coordinates(const lsdtt::lsdtt_raster_meta& meta, nb::object x, nb::object y,
                            nb::object longitude, nb::object latitude, nb::object points,
                            const std::optional<std::string>& crs,
                            std::vector<double>& xs, std::vector<double>& ys)
{
    const bool have_xy = !x.is_none() || !y.is_none();
    const bool have_lonlat = !longitude.is_none() || !latitude.is_none();
    const bool have_points = !points.is_none();
    if (static_cast<int>(have_xy) + static_cast<int>(have_lonlat) + static_cast<int>(have_points) != 1)
        throw std::invalid_argument("snap_points: give exactly one of x and y, longitude and latitude, or points");
    if (have_points)
    {
        if (crs) throw std::invalid_argument("snap_points: crs is not used with points; the GeoDataFrame's own CRS is used");
        nb::object frame = points;
        if (!frame.attr("crs").is_none() && !meta.crs_wkt.empty())
            frame = frame.attr("to_crs")(meta.crs_wkt);
        else if (frame.attr("crs").is_none() && !meta.crs_wkt.empty())
            throw std::invalid_argument("snap_points: points must have a CRS (points.crs is None)");
        nb::object geometry = frame.attr("geometry");
        // .any() gives numpy.bool_, which nanobind's bool caster rejects.
        nb::object any_other = geometry.attr("geom_type").attr("__ne__")("Point").attr("any")();
        if (PyObject_IsTrue(any_other.ptr()) != 0)
            throw std::invalid_argument("snap_points: every geometry in points must be a Point");
        xs = nb::cast<std::vector<double>>(geometry.attr("x").attr("to_numpy")(nb::arg("dtype") = "float64"));
        ys = nb::cast<std::vector<double>>(geometry.attr("y").attr("to_numpy")(nb::arg("dtype") = "float64"));
        return;
    }
    if (have_xy && (x.is_none() || y.is_none()))
        throw std::invalid_argument("snap_points: give both x and y");
    if (have_lonlat && (longitude.is_none() || latitude.is_none()))
        throw std::invalid_argument("snap_points: give both longitude and latitude");
    xs = coordinate_values(have_xy ? x : longitude, have_xy ? "x" : "longitude");
    ys = coordinate_values(have_xy ? y : latitude, have_xy ? "y" : "latitude");
    if (xs.size() != ys.size())
        throw std::invalid_argument("snap_points: the two coordinate sequences differ in length");
    if (have_lonlat)
        transform_points_to_raster_crs(xs, ys, crs.value_or("EPSG:4326"), meta);
    else if (crs)
        transform_points_to_raster_crs(xs, ys, *crs, meta);
}

// Contributing-pixel threshold of snap_points: threshold_contributing_pixels
// as given, or the smallest whole number of pixels N with
// N * pixel_area >= threshold_drainage_area_m2. The C++ snapping routines
// accept a pixel when its contributing pixels are >= the threshold, so a pixel
// qualifies exactly when its drainage area (contributing pixels * pixel area)
// is >= threshold_drainage_area_m2. With neither, the CLIs' default
// basin_outlet_contributing_pixels_threshold.
std::uint64_t snap_threshold_pixels(const lsdtt::lsdtt_raster_meta& meta,
                                    std::optional<std::int64_t> threshold_contributing_pixels,
                                    std::optional<double> threshold_drainage_area_m2)
{
    if (threshold_contributing_pixels && threshold_drainage_area_m2)
        throw std::invalid_argument(
            "snap_points: give threshold_contributing_pixels or threshold_drainage_area_m2, not both");
    if (threshold_contributing_pixels)
    {
        if (*threshold_contributing_pixels < 1)
            throw std::invalid_argument("snap_points: threshold_contributing_pixels must be >= 1");
        return static_cast<std::uint64_t>(*threshold_contributing_pixels);
    }
    if (!threshold_drainage_area_m2)
        return static_cast<std::uint64_t>(lsdtt::default_basin_outlet_contributing_pixels_threshold);
    const double area = *threshold_drainage_area_m2;
    const double pixel_area = meta.pixel_area();
    if (!std::isfinite(area) || area <= 0.0)
        throw std::invalid_argument("snap_points: threshold_drainage_area_m2 must be positive");
    const double ratio = std::ceil(area / pixel_area);
    if (!(ratio < 9.0e15))
        throw std::invalid_argument("snap_points: threshold_drainage_area_m2 is too large for the grid");
    // ceil(area / pixel_area) can be one off when the division rounds; settle
    // it on the product that defines the threshold.
    auto n = static_cast<std::uint64_t>(ratio);
    while (n > 1 && static_cast<double>(n - 1) * pixel_area >= area) --n;
    while (static_cast<double>(n) * pixel_area < area) ++n;
    return n < 1 ? 1 : n;
}

nb::object snap_points(const FlowInfo& flow_owner, nb::object x, nb::object y,
                       nb::object longitude, nb::object latitude, nb::object points,
                       const std::optional<std::string>& crs, const std::string& method,
                       std::optional<std::int64_t> threshold_contributing_pixels,
                       std::optional<double> threshold_drainage_area_m2,
                       int min_search_radius_pixels, int max_search_radius_pixels,
                       std::optional<double> max_search_radius_metres)
{
    const auto& flow = flow_owner.native();
    const auto& meta = flow.meta();
    if (method != "downstream" && method != "expanding_square" && method != "largest_area")
        throw std::invalid_argument(
            "snap_points: method must be 'downstream', 'expanding_square' or 'largest_area'");
    if (min_search_radius_pixels < 0)
        throw std::invalid_argument("snap_points: min_search_radius_pixels must be non-negative");
    if (max_search_radius_pixels < 0)
        throw std::invalid_argument("snap_points: max_search_radius_pixels must be non-negative");
    int max_radius = max_search_radius_pixels;
    if (max_search_radius_metres)
    {
        if (!std::isfinite(*max_search_radius_metres) || *max_search_radius_metres < 0.0)
            throw std::invalid_argument("snap_points: max_search_radius_metres must be non-negative");
        // As erosion_rate_snap_max_distance_metres in the CLIs.
        max_radius = static_cast<int>(std::ceil(*max_search_radius_metres / meta.cell_size - 1e-9));
    }
    if (method == "expanding_square" && max_radius < min_search_radius_pixels)
        throw std::invalid_argument(
            "snap_points: the maximum search radius must be >= min_search_radius_pixels");

    const std::uint64_t threshold_pixels =
        snap_threshold_pixels(meta, threshold_contributing_pixels, threshold_drainage_area_m2);

    std::vector<double> xs, ys;
    snap_input_coordinates(meta, std::move(x), std::move(y), std::move(longitude),
                           std::move(latitude), std::move(points), crs, xs, ys);

    // The same options and the same dispatcher as the CLIs'
    // basin_outlet_snapping_method, basin_outlet_contributing_pixels_threshold
    // and basin_outlet_snap_min/max_search_radius_pixels.
    lsdtt::lsdtt_square_snap_options square;
    square.contributing_pixels_threshold = threshold_pixels;
    square.min_search_radius_pixels = min_search_radius_pixels;
    square.max_search_radius_pixels = max_radius;
    lsdtt::lsdtt_downstream_snap_options downstream;
    downstream.contributing_pixels_threshold = threshold_pixels;

    std::vector<std::int64_t> output;
    output.reserve(xs.size() * 2);
    std::vector<std::string> failures;
    for (std::size_t i = 0; i < xs.size(); ++i)
    {
        int row = -1, col = -1;
        std::ostringstream why;
        why.precision(12);
        if (!pixel_of_xy(meta, xs[i], ys[i], row, col))
            why << "point " << i << " (x=" << xs[i] << ", y=" << ys[i] << " in the grid CRS): "
                << lsdtt::snap_status_to_string(lsdtt::lsdtt_snap_status::input_out_of_bounds);
        else
        {
            const auto snapped = lsdtt::snap_pixel_by_method(flow, row, col, method, square, downstream);
            if (snapped.status != lsdtt::lsdtt_snap_status::success)
            {
                why << "point " << i << " at pixel (" << row << ", " << col << "): "
                    << lsdtt::snap_status_to_string(snapped.status) << " (" << method
                    << ", threshold " << threshold_pixels << " contributing pixels";
                if (method != "downstream") why << ", max_search_radius_pixels " << max_radius;
                why << ")";
            }
            else { row = snapped.snapped_row; col = snapped.snapped_col; }
        }
        if (!why.str().empty()) failures.push_back(why.str());
        output.push_back(row);
        output.push_back(col);
    }
    if (!failures.empty())
    {
        std::ostringstream message;
        message << "snap_points: " << failures.size() << " of " << xs.size() << " points could not be snapped: ";
        for (std::size_t k = 0; k < failures.size() && k < 10; ++k) message << (k ? "; " : "") << failures[k];
        if (failures.size() > 10) message << "; and " << failures.size() - 10 << " more";
        throw std::invalid_argument(message.str());
    }
    nb::object numpy = nb::module_::import_("numpy");
    return numpy.attr("array")(output, nb::arg("dtype") = "int64")
        .attr("reshape")(nb::make_tuple(static_cast<int>(xs.size()), 2));
}
} // namespace

FlowInfo::FlowInfo(const Raster& dem, const std::string& boundary_conditions,
                   bool interior_nodata_is_outlet, nb::object forced_outlets,
                   bool allow_pits_as_outlets, int n_threads)
    : flow_(std::make_shared<lsdtt::lsdtt_flowinfo>(
          make_native_flow_info(dem, boundary_conditions, interior_nodata_is_outlet,
                                std::move(forced_outlets), allow_pits_as_outlets, n_threads)))
{}

void bind_flow_info(nb::module_& module)
{
    nb::class_<FlowInfo> flow(module, "FlowInfo", R"doc(
D8 drainage topology (steepest-descent receivers) of a conditioned DEM.

Each valid pixel drains to the one of its eight neighbours with the steepest
downhill gradient, or is an outlet. ``FlowInfo`` stores this routing (and
each pixel's upslope pixels), not the elevations, so later calls that need
elevations take the DEM again. Fill or breach depressions first (for example
``fill_pits``); a pixel with no downhill neighbour and no allowed exit raises
an error unless ``allow_pits_as_outlets=True``. Pixels are addressed by
``(row, col)`` (row 0 is the northern edge) or by flat index
``row * n_cols + col``.

Typical workflow::

    F = lsdtt3.FlowInfo(lsdtt3.fill_pits(dem))
    scope = F.channel_scope(threshold_contributing_pixels=1000)
    N = F.channel_network(scope)
    B = F.basins(outlets, min_upslope_pixels=0)

Parameters
----------
dem : Raster
    Depression-conditioned elevation raster (m) in a projected metre CRS.
boundary_conditions : str, default='oooo'
    Four characters for the north, east, south and west edges: ``o`` open
    (an edge pixel may drain out of the grid when it has no lower in-grid
    neighbour), ``n`` no-flux (flow cannot leave), ``p`` periodic (wraps to
    the opposite edge; must be given for both opposite edges), ``b``
    base level (every pixel on the edge is an outlet). See
    ``make_boundary_conditions``.
interior_nodata_is_outlet : bool, default=True
    If True, a pixel next to a nodata pixel may drain into it (become an
    outlet) when it has no lower valid neighbour, as at a coastline or the
    edge of a clipped DEM. If False such pixels must drain over valid pixels.
forced_outlets : array-like of int or None, default=None
    N x 2 array of ``(row, col)`` pixels that are made outlets regardless of
    their neighbours (moulins, swallow holes, playas). Prepare the DEM with
    ``fill_pits_with_sinks`` using the same pixels. Out-of-grid or nodata
    pixels raise ``ValueError``; duplicates are ignored.
allow_pits_as_outlets : bool, default=False
    If True, a pixel that has no downhill neighbour and no allowed exit
    becomes an outlet instead of raising an error (useful in landscape
    evolution runs where transient pits form).
n_threads : int, default=0
    OpenMP thread count; 0 uses the runtime default.
)doc");
    flow.def(nb::init<const Raster&, const std::string&, bool, nb::object, bool, int>(), nb::arg("dem"),
        nb::arg("boundary_conditions") = "oooo",
        nb::arg("interior_nodata_is_outlet") = true,
        nb::arg("forced_outlets") = nb::none(),
        nb::arg("allow_pits_as_outlets") = false, nb::arg("n_threads") = 0,
        R"doc(Build D8 flow topology from a depression-conditioned DEM.

``boundary_conditions`` is a north/east/south/west string using open (``o``),
no-flux (``n``), periodic (``p``), or base-level (``b``) edges. Forced outlets
are an N by 2 array of ``(row, col)`` surface-flow terminations.)doc")
        .def_prop_ro("metadata", [](const FlowInfo& self) -> const lsdtt::lsdtt_raster_meta& { return self.native().meta(); }, nb::rv_policy::reference_internal,
            "Grid metadata (RasterMeta) of the DEM the topology was built from.")
        .def_prop_ro("shape", [](const FlowInfo& self) { return nb::make_tuple(self.native().n_rows(), self.native().n_cols()); },
            nb::sig("def shape(self) -> tuple[int, int]"), "Grid shape as (n_rows, n_cols).")
        .def_prop_ro("boundary_conditions", [](const FlowInfo& self) { return self.native().boundary_conditions(); },
            "The four-character NESW boundary-condition string used.")
        .def_prop_ro("n_valid_cells", [](const FlowInfo& self) { return self.native().n_valid_cells(); },
            "Number of valid (routed) pixels.")
        .def_prop_ro("n_outlets", [](const FlowInfo& self) { return self.native().n_outlets(); },
            "Number of outlet pixels (pixels with no receiver).")
        .def_prop_ro("outlets", &outlet_pixels, nb::sig("def outlets(self) -> numpy.ndarray"),
            "Outlet pixels as an N x 2 int64 array of (row, col), in flat-index order.")
        .def_prop_ro("valid_mask", [](const FlowInfo& self) { const auto& f = self.native(); return grid_array(f.mask_array().data(), f.full_grid_size(), f.meta(), "bool"); },
            nb::sig("def valid_mask(self) -> numpy.ndarray"),
            "Bool array (n_rows, n_cols), True for valid pixels.")
        .def_prop_ro("receiver", [](const FlowInfo& self) { const auto& f = self.native(); return grid_array(f.receiver_array().data(), f.full_grid_size(), f.meta(), "int64"); },
            nb::sig("def receiver(self) -> numpy.ndarray"),
            "Int64 array (n_rows, n_cols) of each pixel's receiver flat index; -2 at outlets, -1 at nodata.")
        .def_prop_ro("direction", [](const FlowInfo& self) { const auto& f = self.native(); return grid_array(f.receiver_dir_array().data(), f.full_grid_size(), f.meta(), "uint8"); },
            nb::sig("def direction(self) -> numpy.ndarray"),
            "Uint8 array (n_rows, n_cols) of internal D8 direction codes; see flow_direction() for a raster.")
        .def_prop_ro("donor_count", [](const FlowInfo& self) { const auto& f = self.native(); return grid_array(f.donor_count_array().data(), f.full_grid_size(), f.meta(), "uint8"); },
            nb::sig("def donor_count(self) -> numpy.ndarray"),
            "Uint8 array (n_rows, n_cols): number of neighbours draining directly into each pixel.")
        .def("flat_index", [](const FlowInfo& self, int row, int col) { return self.native().flat_index(row, col); }, nb::arg("row"), nb::arg("col"), R"doc(
Return the flat index ``row * n_cols + col`` of a pixel.

Parameters
----------
row, col : int
    Zero-based pixel position; row 0 is the northern edge.

Returns
-------
int
    Flat index, as used in the ``flat`` columns of lsdtt3 tables. Raises
    ``RuntimeError`` if the pixel is outside the grid.
)doc")
        .def("row_col", [](const FlowInfo& self, std::int64_t flat) { return nb::make_tuple(self.native().row_from_flat(flat), self.native().col_from_flat(flat)); }, nb::arg("flat_index"), R"doc(
Return the ``(row, col)`` of a flat index.

Parameters
----------
flat_index : int
    Flat index ``row * n_cols + col`` (for example a ``flat`` or
    ``outlet_flat`` table value).

Returns
-------
tuple of (int, int)
    ``(row, col)``; equivalent to ``divmod(flat_index, n_cols)``. Raises
    ``RuntimeError`` if the index is out of range.
)doc")
        .def("is_outlet", [](const FlowInfo& self, int row, int col) { return self.native().is_outlet(self.native().flat_index(row, col)); }, nb::arg("row"), nb::arg("col"), R"doc(
Return True if the pixel is an outlet (has no downstream receiver).

Parameters
----------
row, col : int
    Zero-based pixel position of a valid pixel.

Returns
-------
bool
    Raises ``RuntimeError`` for an out-of-grid or nodata pixel.
)doc")
        .def("receiver_of", &receiver_pixel, nb::arg("row"), nb::arg("col"), R"doc(
Return the pixel this pixel drains to.

Parameters
----------
row, col : int
    Zero-based pixel position of a valid pixel.

Returns
-------
tuple of (int, int) or None
    ``(row, col)`` of the D8 downstream neighbour, or None if the pixel is an
    outlet. Raises ``RuntimeError`` for an out-of-grid or nodata pixel.
)doc")
        .def("upslope_pixels", &upslope_pixels, nb::arg("row"), nb::arg("col"), R"doc(
Return every pixel that drains through this pixel, including itself.

Parameters
----------
row, col : int
    Zero-based pixel position of a valid pixel.

Returns
-------
numpy.ndarray
    N x 2 int64 array of ``(row, col)``; the first row is the pixel itself.
    N is the pixel's contributing-pixel count.
)doc")
        .def("is_upslope", [](const FlowInfo& self, int cr, int cc, int tr, int tc) { return self.native().is_upslope_of(self.native().flat_index(cr, cc), self.native().flat_index(tr, tc)); }, nb::arg("candidate_row"), nb::arg("candidate_col"), nb::arg("target_row"), nb::arg("target_col"), R"doc(
Return True if the candidate pixel drains through the target pixel.

Parameters
----------
candidate_row, candidate_col : int
    Pixel to test.
target_row, target_col : int
    Pixel whose upslope area is searched. A pixel is upslope of itself.

Returns
-------
bool
)doc")
        .def("snap_points", &snap_points, nb::arg("x") = nb::none(), nb::arg("y") = nb::none(),
            nb::kw_only(), nb::arg("longitude") = nb::none(), nb::arg("latitude") = nb::none(),
            nb::arg("points") = nb::none(), nb::arg("crs") = nb::none(),
            nb::arg("method") = "downstream",
            nb::arg("threshold_contributing_pixels") = nb::none(),
            nb::arg("threshold_drainage_area_m2") = nb::none(),
            nb::arg("min_search_radius_pixels") = 0, nb::arg("max_search_radius_pixels") = 25,
            nb::arg("max_search_radius_metres") = nb::none(), R"doc(
Snap points (for example basin outlets or sample sites) onto the flow network.

Each point is placed in the pixel that contains it and then moved by
``method``. The result is the ``(row, col)`` array that
``FlowInfo.basins(outlets)``, ``ChannelNetwork.chi(outlets=...)`` and other
``outlets`` arguments take. Give the points as exactly one of ``x`` and
``y``, ``longitude`` and ``latitude``, or ``points``.

The methods, their names and their rules are those of the command-line
programs' ``basin_outlet_snapping_method`` (the same C++ routines): a
``method``, ``threshold_contributing_pixels`` and the two radii here give
the same pixel as ``basin_outlet_snapping_method``,
``basin_outlet_contributing_pixels_threshold`` and
``basin_outlet_snap_min/max_search_radius_pixels`` with the same values. A
pixel "qualifies" when its contributing-pixel count (the pixel itself and
all pixels draining through it) is at least the threshold (``>=``).

* ``"downstream"``: follow the D8 receivers from the point's pixel to the
  first qualifying pixel. The point must be on a valid pixel and upstream
  of a qualifying pixel. The radii are not used.
* ``"expanding_square"``: search square rings of half-width
  ``min_search_radius_pixels`` to ``max_search_radius_pixels`` around the
  point's pixel, one ring at a time, and take the nearest qualifying pixel
  (straight-line distance) in the first ring that has one; ties go to the
  larger drainage area, then the smaller flat index.
* ``"largest_area"``: the pixel with the largest drainage area in the
  square of half-width ``max_search_radius_pixels`` (ties: the nearer
  pixel, then the smaller flat index), provided it qualifies. This is the
  GIS "snap pour point" rule. With a wide square it can jump onto a bigger
  neighbouring river, so keep the square to a few pixels.
  ``min_search_radius_pixels`` is not used.

Error policy: if any point is outside the grid or cannot be snapped,
``ValueError`` is raised naming each such point (by its zero-based position
in the input) and the snapping status (``input_out_of_bounds``,
``input_invalid_flow_cell``, ``threshold_not_found_in_search_window`` or
``outlet_reached_before_threshold``); no partial result is returned.

Parameters
----------
x, y : float or array-like of float, optional
    Easting and northing of the points, in the flow grid's CRS
    (``F.metadata.crs_wkt``, metres) unless ``crs`` is given.
longitude, latitude : float or array-like of float, optional, keyword-only
    Point coordinates in decimal degrees, in ``crs`` (default WGS84,
    ``"EPSG:4326"``). Longitude is the x coordinate: pass both by keyword.
points : geopandas.GeoDataFrame or geopandas.GeoSeries, optional, keyword-only
    Point geometries with a CRS; they are transformed to the grid's CRS.
crs : str or None, default=None, keyword-only
    CRS of ``x``/``y`` or ``longitude``/``latitude`` as an ``"EPSG:<code>"``
    string, WKT or PROJ string. None means the grid's CRS for ``x``/``y`` and
    ``"EPSG:4326"`` for ``longitude``/``latitude``. Not used with ``points``.
method : {"downstream", "expanding_square", "largest_area"}, default="downstream", keyword-only
    Snapping rule, as described above.
threshold_contributing_pixels : int or None, default=None, keyword-only
    Contributing pixels (>= 1) a pixel needs to qualify; the CLIs'
    ``basin_outlet_contributing_pixels_threshold``. Give this or
    ``threshold_drainage_area_m2``, not both. With neither, every method
    uses that key's CLI default, 1000 contributing pixels (as a CLI run
    without the key does); pass ``threshold_contributing_pixels=1`` to
    ``"largest_area"`` for no minimum.
threshold_drainage_area_m2 : float or None, default=None, keyword-only
    The same threshold as a drainage area (m^2, > 0). It is converted to the
    smallest whole number of contributing pixels N with
    N * ``F.metadata.pixel_area`` >= ``threshold_drainage_area_m2``, so a
    pixel qualifies exactly when its drainage area (contributing pixels *
    pixel area) is >= ``threshold_drainage_area_m2``.
min_search_radius_pixels : int, default=0, keyword-only
    Smallest ring half-width (pixels, >= 0) for ``"expanding_square"``; 0
    tests the point's own pixel first. As
    ``basin_outlet_snap_min_search_radius_pixels``.
max_search_radius_pixels : int, default=25, keyword-only
    Largest ring half-width for ``"expanding_square"`` and the square's
    half-width for ``"largest_area"`` (pixels, >= 0). As
    ``basin_outlet_snap_max_search_radius_pixels``.
max_search_radius_metres : float or None, default=None, keyword-only
    ``max_search_radius_pixels`` in metres; if given it replaces it and is
    rounded up to whole pixels.

Returns
-------
numpy.ndarray
    N x 2 int64 array of snapped ``(row, col)`` pixels, one row per input
    point and in input order (row 0 is the north edge). Points may snap to
    the same pixel; ``FlowInfo.basins`` then keeps only the first.

Examples
--------
>>> F = lsdtt3.FlowInfo(lsdtt3.fill_pits(dem))
>>> outlets = F.snap_points(longitude=[15.79008, 16.06060], latitude=[37.94234, 38.06143],
...                         method="downstream", threshold_contributing_pixels=500)
>>> B = F.basins(outlets)
)doc");
    bind_flow_products(flow);
    bind_flow_network(module, flow);
    bind_basins(module, flow);
    bind_channel_heads(module, flow);
}
} // namespace lsdtt3_python
