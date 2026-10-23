#include <cstdint>
#include <cctype>
#include <cstring>
#include <cmath>
#include <cstdlib>
#include <limits>
#include <optional>
#include <random>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <nanobind/nanobind.h>
#include <nanobind/stl/optional.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>

#include <cpl_vsi.h>
#include <gdal_utils.h>
#include <gdal_priv.h>
#include <ogr_spatialref.h>

#include <core/lsdtt_data_loader.hpp>
#include <core/lsdtt_raster_meta.hpp>
#include <cpu/lsdtt_raster_algebra.hpp>
#include <cpu/lsdtt_raster_flow_preprocessing_cpu.hpp>
#include <cpu/lsdtt_raster_breach_least_cost_cpu.hpp>
#include <cpu/lsdtt_raster_nodata_filling_cpu.hpp>
#include <cpu/lsdtt_raster_cpu.hpp>
#include <cpu/lsdtt_raster_surface_metrics_cpu.hpp>
#include <cpu/lsdtt_raster_mosaic_cpu.hpp>
#include <cpu/lsdtt_vector_rasterise_cpu.hpp>


#include "binding_common.hpp"

namespace lsdtt3_python {
namespace {
lsdtt::priority_flood_queue_backend priority_flood_backend(const std::string& name)
{
    if (name == "binary") return lsdtt::priority_flood_queue_backend::binary_heap;
    if (name == "radix") return lsdtt::priority_flood_queue_backend::radix_hybrid;
    throw std::invalid_argument("priority_flood_backend must be 'binary' or 'radix'");
}

void require_compatible(const Raster& first, const Raster& second, const char* name)
{
    const auto& a = first.native().get_meta();
    const auto& b = second.native().get_meta();
    if (a.n_rows != b.n_rows || a.n_cols != b.n_cols ||
        a.x_min != b.x_min || a.y_min != b.y_min ||
        a.cell_size != b.cell_size || a.crs_wkt != b.crs_wkt)
    {
        throw std::invalid_argument(std::string(name) + " must share the raster grid");
    }
}

std::vector<std::int64_t> selected_flat_indices(
    const Raster& raster,
    const Raster& selection,
    const char* name)
{
    require_compatible(raster, selection, name);
    const auto& source_mask = raster.native().get_mask();
    const auto& selection_mask = selection.native().get_mask();
    const auto& selection_data = selection.native().get_data();
    std::vector<std::int64_t> indices;
    for (std::size_t i = 0; i < selection_data.size(); ++i)
    {
        if (selection_mask[i] != 0 && selection_data[i] > 0.5f)
        {
            if (source_mask[i] == 0)
            {
                throw std::invalid_argument(std::string(name) + " selects a nodata pixel");
            }
            indices.push_back(static_cast<std::int64_t>(i));
        }
    }
    if (indices.empty())
    {
        throw std::invalid_argument(std::string(name) + " contains no selected pixels");
    }
    return indices;
}

std::vector<std::int64_t> row_col_points_to_flat_indices(
    const Raster& raster,
    nb::object points,
    const char* name)
{
    nb::object numpy = nb::module_::import_("numpy");
    nb::object array = numpy.attr("asarray")(points, nb::arg("dtype") = "int64");
    if (nb::cast<int>(array.attr("ndim")) != 2)
    {
        throw std::invalid_argument(std::string(name) + " must be an N x 2 array of (row, col)");
    }
    nb::tuple shape = nb::cast<nb::tuple>(array.attr("shape"));
    if (nb::cast<int>(shape[1]) != 2)
    {
        throw std::invalid_argument(std::string(name) + " must be an N x 2 array of (row, col)");
    }

    const auto& metadata = raster.native().get_meta();
    const auto& source_mask = raster.native().get_mask();
    std::vector<std::int64_t> coordinates = nb::cast<std::vector<std::int64_t>>(
        array.attr("ravel")()
    );
    std::vector<std::int64_t> indices;
    indices.reserve(coordinates.size() / 2);
    for (std::size_t i = 0; i < coordinates.size(); i += 2)
    {
        const std::int64_t row = coordinates[i];
        const std::int64_t col = coordinates[i + 1];
        if (row < 0 || row >= metadata.n_rows || col < 0 || col >= metadata.n_cols)
        {
            throw std::invalid_argument(std::string(name) + " contains a point outside the raster");
        }
        const std::int64_t flat = row * metadata.n_cols + col;
        if (source_mask[static_cast<std::size_t>(flat)] == 0)
        {
            throw std::invalid_argument(std::string(name) + " selects a nodata pixel");
        }
        indices.push_back(flat);
    }
    if (indices.empty())
    {
        throw std::invalid_argument(std::string(name) + " contains no points");
    }
    return indices;
}

Raster make_outlet_mask(
    const Raster& raster,
    bool north,
    bool east,
    bool south,
    bool west)
{
    const auto& metadata = raster.native().get_meta();
    std::vector<float> data(static_cast<std::size_t>(metadata.n_pixels()), 0.0f);
    std::vector<std::uint8_t> mask(static_cast<std::size_t>(metadata.n_pixels()), 1);
    const auto set_outlet = [&](int row, int col) {
        data[static_cast<std::size_t>(row) * metadata.n_cols + col] = 1.0f;
    };
    if (north) for (int col = 0; col < metadata.n_cols; ++col) set_outlet(0, col);
    if (east) for (int row = 0; row < metadata.n_rows; ++row) set_outlet(row, metadata.n_cols - 1);
    if (south) for (int col = 0; col < metadata.n_cols; ++col) set_outlet(metadata.n_rows - 1, col);
    if (west) for (int row = 0; row < metadata.n_rows; ++row) set_outlet(row, 0);
    return Raster(lsdtt::lsdtt_raster_cpu(std::move(data), std::move(mask), metadata));
}

Raster mask_from_points(const Raster& raster, nb::object points)
{
    const auto& metadata = raster.native().get_meta();
    std::vector<float> data(static_cast<std::size_t>(metadata.n_pixels()), 0.0f);
    std::vector<std::uint8_t> mask(static_cast<std::size_t>(metadata.n_pixels()), 1);
    for (const std::int64_t index : row_col_points_to_flat_indices(raster, points, "points"))
    {
        data[static_cast<std::size_t>(index)] = 1.0f;
    }
    return Raster(lsdtt::lsdtt_raster_cpu(std::move(data), std::move(mask), metadata));
}

char boundary_character(const std::string& value)
{
    std::string lower;
    lower.reserve(value.size());
    for (const char character : value)
    {
        lower.push_back(static_cast<char>(std::tolower(static_cast<unsigned char>(character))));
    }
    if (lower == "open" || lower == "o") return 'o';
    if (lower == "closed" || lower == "no_flux" || lower == "n") return 'n';
    if (lower == "periodic" || lower == "p") return 'p';
    if (lower == "base_level" || lower == "base-level" || lower == "b") return 'b';
    throw std::invalid_argument("boundary mode must be open, closed, periodic, or base_level");
}

std::string make_boundary_conditions(
    const std::string& north,
    const std::string& east,
    const std::string& south,
    const std::string& west)
{
    return std::string{
        boundary_character(north), boundary_character(east),
        boundary_character(south), boundary_character(west)
    };
}

std::string validate_boundary_conditions(const std::string& boundary_conditions)
{
    if (boundary_conditions.size() != 4)
    {
        throw std::invalid_argument("boundary_conditions must have four NESW characters");
    }
    std::string validated;
    validated.reserve(4);
    for (const char character : boundary_conditions)
    {
        validated.push_back(boundary_character(std::string(1, character)));
    }
    const bool periodic_ns = validated[0] == 'p' || validated[2] == 'p';
    const bool periodic_ew = validated[1] == 'p' || validated[3] == 'p';
    if ((validated[0] == 'p') != (validated[2] == 'p') ||
        (validated[1] == 'p') != (validated[3] == 'p'))
    {
        throw std::invalid_argument("periodic boundaries must occur in opposite pairs");
    }
    if (periodic_ns && periodic_ew)
    {
        throw std::invalid_argument("both axes cannot be periodic because no drainage exit exists");
    }
    return validated;
}

Raster fill_pits_with_outlets(
    const Raster& raster,
    const Raster& outlet_mask,
    float min_slope,
    int n_threads,
    const std::string& backend)
{
    const std::vector<std::int64_t> outlets = selected_flat_indices(
        raster, outlet_mask, "outlet_mask"
    );
    const auto& metadata = raster.native().get_meta();
    std::vector<int> rows;
    std::vector<int> cols;
    rows.reserve(outlets.size());
    cols.reserve(outlets.size());
    for (const std::int64_t index : outlets)
    {
        rows.push_back(static_cast<int>(index / metadata.n_cols));
        cols.push_back(static_cast<int>(index % metadata.n_cols));
    }
    return Raster(lsdtt::fill_pits_with_boundary(
        raster.native(), min_slope, {rows}, {cols}, {"outlet"}, "outlet", n_threads,
        priority_flood_backend(backend)
    ));
}

Raster fill_pits_with_sinks_mask(
    const Raster& raster,
    const Raster& sink_mask,
    float min_slope,
    int n_threads,
    const std::string& backend)
{
    return Raster(lsdtt::fill_pits_with_sinks(
        raster.native(), selected_flat_indices(raster, sink_mask, "sink_mask"), min_slope,
        n_threads, priority_flood_backend(backend)
    ));
}

Raster fill_pits_with_sinks_points(
    const Raster& raster,
    nb::object points,
    float min_slope,
    int n_threads,
    const std::string& backend)
{
    return Raster(lsdtt::fill_pits_with_sinks(
        raster.native(), row_col_points_to_flat_indices(raster, points, "sinks"), min_slope,
        n_threads, priority_flood_backend(backend)
    ));
}

Raster fill_pits_boundary_aware(
    const Raster& raster,
    const std::string& boundary_conditions,
    float min_slope,
    int n_threads,
    const std::string& backend)
{
    const std::string boundary = validate_boundary_conditions(boundary_conditions);
    const bool periodic = boundary.find('p') != std::string::npos;
    if (periodic)
    {
        for (const char character : boundary)
        {
            if (character != 'p' && character != 'o')
            {
                throw std::invalid_argument(
                    "periodic filling currently supports open boundaries on the other axis"
                );
            }
        }
        return Raster(lsdtt::fill_pits(
            raster.native(), min_slope, boundary, n_threads,
            priority_flood_backend(backend)));
    }

    // Without a closed edge this is the fill the CLIs run: the raster border
    // and every valid cell touching nodata are spill points.
    if (boundary.find('n') == std::string::npos)
    {
        return Raster(lsdtt::fill_pits(
            raster.native(), min_slope, boundary, n_threads,
            priority_flood_backend(backend)));
    }

    // Closed edges: seed the valid pixels of the other edges, plus every valid
    // cell touching nodata (a closed edge included), because nodata in a
    // downloaded DEM marks where the data ends, not a wall.
    const auto& metadata = raster.native().get_meta();
    const auto& mask = raster.native().get_mask();
    const int n_rows = metadata.n_rows;
    const int n_cols = metadata.n_cols;
    const auto valid = [&](int row, int col) {
        return mask[static_cast<std::size_t>(row) * n_cols + col] != 0;
    };
    std::vector<int> rows;
    std::vector<int> cols;
    for (int row = 0; row < n_rows; ++row)
    {
        for (int col = 0; col < n_cols; ++col)
        {
            if (!valid(row, col)) continue;
            bool seed = (row == 0 && boundary[0] != 'n') ||
                        (col == n_cols - 1 && boundary[1] != 'n') ||
                        (row == n_rows - 1 && boundary[2] != 'n') ||
                        (col == 0 && boundary[3] != 'n');
            for (int dr = -1; dr <= 1 && !seed; ++dr)
            {
                for (int dc = -1; dc <= 1 && !seed; ++dc)
                {
                    const int r = row + dr;
                    const int c = col + dc;
                    if (r >= 0 && r < n_rows && c >= 0 && c < n_cols && !valid(r, c))
                    {
                        seed = true;
                    }
                }
            }
            if (seed)
            {
                rows.push_back(row);
                cols.push_back(col);
            }
        }
    }
    if (rows.empty())
    {
        throw std::invalid_argument(
            "boundary_conditions leave no outlet: every edge is closed and no "
            "valid pixel touches nodata"
        );
    }
    return Raster(lsdtt::fill_pits_with_boundary(
        raster.native(), min_slope, {rows}, {cols}, {"outlet"}, "outlet", n_threads,
        priority_flood_backend(backend)
    ));
}

lsdtt::breach_least_cost_options least_cost_options(
    std::optional<int> max_distance,
    float max_cost,
    bool prefer_shortest_path,
    float flat_increment,
    bool fill_remainder,
    const std::string& edge_policy,
    float fill_min_slope)
{
    lsdtt::breach_least_cost_options options;
    options.max_cost       = (max_cost >= (std::numeric_limits<float>::max)())
                           ? std::numeric_limits<double>::infinity()
                           : static_cast<double>(max_cost);
    options.minimize_dist  = prefer_shortest_path;
    options.flat_increment = static_cast<double>(flat_increment);
    options.fill_remainder = fill_remainder;
    options.fill_min_slope = fill_min_slope;
    options.edge_policy    = lsdtt::parse_breach_lc_edge_policy(edge_policy, "edge_policy");
    // None: the edge policy's default (1000 outlet, to bridge large water
    // bodies; 100 whitebox), as in lsdtt-raster-preprocessing.
    options.max_dist       = max_distance.has_value()
                           ? *max_distance
                           : lsdtt::default_breach_lc_max_dist(options.edge_policy);
    return options;
}

Raster breach_least_cost(
    const Raster& raster,
    std::optional<int> max_distance,
    float max_cost,
    bool prefer_shortest_path,
    float flat_increment,
    bool fill_remainder,
    int /*n_threads*/,
    const std::string& edge_policy,
    float fill_min_slope)
{
    return Raster(lsdtt::breach_depressions_least_cost(
        raster.native(),
        least_cost_options(max_distance, max_cost, prefer_shortest_path, flat_increment,
                           fill_remainder, edge_policy, fill_min_slope)));
}

Raster least_cost_carve_depth(
    const Raster& raster,
    std::optional<int> max_distance,
    float max_cost,
    bool prefer_shortest_path,
    float flat_increment,
    bool fill_remainder,
    int /*n_threads*/,
    const std::string& edge_policy,
    float fill_min_slope)
{
    const lsdtt::lsdtt_raster_cpu breached = lsdtt::breach_depressions_least_cost(
        raster.native(),
        least_cost_options(max_distance, max_cost, prefer_shortest_path, flat_increment,
                           fill_remainder, edge_policy, fill_min_slope));
    return Raster(lsdtt::algebra::raster_positive_difference(raster.native(), breached));
}

Raster fill_pits(const Raster& raster, float min_slope, int n_threads,
                 const std::string& backend)
{
    return Raster(lsdtt::fill_pits(
        raster.native(), min_slope, n_threads, priority_flood_backend(backend)));
}

Raster fill_pits_priority_flood(
    const Raster& raster, float min_slope, int n_threads, bool use_radix_heap)
{
    return Raster(lsdtt::fill_pits_priority_flood(
        raster.native(), min_slope, n_threads, use_radix_heap));
}

Raster fill_depth(const Raster& raster, float min_slope, int n_threads)
{
    const lsdtt::lsdtt_raster_cpu filled = lsdtt::fill_pits(
        raster.native(), min_slope, n_threads
    );
    return Raster(lsdtt::algebra::raster_positive_difference(filled, raster.native()));
}

Raster breach_depressions(const Raster& raster, int n_threads)
{
    return Raster(lsdtt::breach_depressions(raster.native(), n_threads));
}

Raster carve_depth(const Raster& raster, int n_threads)
{
    const lsdtt::lsdtt_raster_cpu breached = lsdtt::breach_depressions(
        raster.native(), n_threads
    );
    return Raster(lsdtt::algebra::raster_positive_difference(raster.native(), breached));
}
} // namespace

void bind_raster_pit_processing(nb::module_& m)
{
    m.def("fill_pits", &fill_pits, nb::arg("raster"),
        nb::arg("min_slope") = 0.0001f, nb::arg("n_threads") = 0,
        nb::arg("priority_flood_backend") = "binary", R"doc(
Fill depressions while imposing an optional minimum drainage slope.

Priority-Flood filling (Barnes et al. 2014): every pit and closed depression
is raised to its spill level, so all water can drain to the raster edge or
to nodata. The raster edge and pixels next to nodata act as outlets. Use the
result to build a ``FlowInfo``.

Parameters
----------
raster : Raster
    Elevation raster (m).
min_slope : float, default=0.0001
    Gradient (m/m) imposed across filled flats: each filled pixel is raised
    at least ``min_slope * cell_size`` (``* sqrt(2)`` on a diagonal step)
    above the pixel it drains to, so every pixel has a downhill route. 0
    fills to a flat surface.
n_threads : int, default=0
    OpenMP thread count; 0 uses the runtime default.
priority_flood_backend : {"binary", "radix"}, default="binary"
    Priority-queue implementation: ``"binary"`` heap or ``"radix"``
    (hybrid radix heap). Both give the same result; only speed differs.

Returns
-------
Raster
    Filled DEM (m).
)doc");
    m.def("fill_pits_priority_flood", &fill_pits_priority_flood,
        nb::arg("raster"), nb::arg("min_slope") = 0.0001f,
        nb::arg("n_threads") = 0, nb::arg("use_radix_heap") = false, R"doc(
Explicit Improved Priority-Flood; optionally use its radix-heap queue.

As ``fill_pits`` but always uses the Priority-Flood algorithm, even when
``min_slope`` is 0.

Parameters
----------
raster : Raster
    Elevation raster (m).
min_slope : float, default=0.0001
    Gradient (m/m) imposed across filled flats: each filled pixel is raised
    at least ``min_slope * cell_size`` (``* sqrt(2)`` on a diagonal step)
    above the pixel it drains to, so every pixel has a downhill route. 0
    fills to a flat surface.
n_threads : int, default=0
    OpenMP thread count; 0 uses the runtime default.
use_radix_heap : bool, default=False
    Use the radix-heap priority queue instead of a binary heap.

Returns
-------
Raster
    Filled DEM (m).
)doc");
    m.def("fill_depth", &fill_depth, nb::arg("raster"),
        nb::arg("min_slope") = 0.0001f, nb::arg("n_threads") = 0, R"doc(
Return the non-negative elevation added by pit filling.

Parameters
----------
raster : Raster
    Elevation raster (m).
min_slope : float, default=0.0001
    Gradient (m/m) imposed across filled flats: each filled pixel is raised
    at least ``min_slope * cell_size`` (``* sqrt(2)`` on a diagonal step)
    above the pixel it drains to, so every pixel has a downhill route. 0
    fills to a flat surface.
n_threads : int, default=0
    OpenMP thread count; 0 uses the runtime default.

Returns
-------
Raster
    ``fill_pits(raster) - raster`` in metres (0 outside depressions).
)doc");
    m.def("breach_depressions", &breach_depressions, nb::arg("raster"),
        nb::arg("n_threads") = 0, R"doc(
Breach depressions by carving drainage paths.

Hybrid breaching (Lindsay 2016): each pit is drained by lowering a path out
of it rather than raising the depression, which changes less of the DEM
than filling in incised terrain. See also ``breach_least_cost``.

Parameters
----------
raster : Raster
    Elevation raster (m).
n_threads : int, default=0
    OpenMP thread count; 0 uses the runtime default.

Returns
-------
Raster
    Breached DEM (m).
)doc");
    m.def("carve_depth", &carve_depth, nb::arg("raster"), nb::arg("n_threads") = 0, R"doc(
Return the non-negative elevation removed by depression breaching.

Parameters
----------
raster : Raster
    Elevation raster (m).
n_threads : int, default=0
    OpenMP thread count; 0 uses the runtime default.

Returns
-------
Raster
    ``raster - breach_depressions(raster)`` in metres.
)doc");
    m.def("make_boundary_conditions", &make_boundary_conditions,
        nb::arg("north") = "open", nb::arg("east") = "open",
        nb::arg("south") = "open", nb::arg("west") = "open", R"doc(
Create an LSDTopoTools NESW boundary-condition string.

Parameters
----------
north, east, south, west : {"open", "closed", "periodic", "base_level"}
    Boundary treatment of each edge (``"no_flux"`` and the one-letter codes
    ``o``, ``n``, ``p``, ``b`` are also accepted). ``open``: flow may leave
    across the edge; ``closed``: no flow leaves; ``periodic``: the edge
    wraps to the opposite one (give both opposite edges); ``base_level``:
    every edge pixel is an outlet.

Returns
-------
str
    Four-character NESW code such as ``"onon"``, for ``FlowInfo``,
    ``fill_pits_boundary_aware`` or ``run_graphflood``.
)doc");
    m.def("make_outlet_mask", &make_outlet_mask, nb::arg("raster"),
        nb::arg("north") = false, nb::arg("east") = false,
        nb::arg("south") = false, nb::arg("west") = false, R"doc(
Create a valid binary mask with selected raster edges set to one.

Parameters
----------
raster : Raster
    Raster whose grid the mask uses.
north, east, south, west : bool, default=False
    Set the pixels of that edge (row 0, last column, last row, column 0) to 1.

Returns
-------
Raster
    0/1 raster with every pixel valid, for ``fill_pits_with_outlets`` or as
    an outlet mask for ``run_graphflood``.
)doc");
    m.def("mask_from_points", &mask_from_points, nb::arg("raster"), nb::arg("points"), R"doc(
Create a binary mask from zero-based ``(row, col)`` point coordinates.

Parameters
----------
raster : Raster
    Raster whose grid the mask uses.
points : array-like, shape (n, 2)
    Integer ``(row, col)`` pixels (row 0 is the northern edge); each must be
    a valid pixel of ``raster``.

Returns
-------
Raster
    0/1 raster, 1 at the points, every pixel valid.
)doc");
    m.def("fill_pits_with_outlets", &fill_pits_with_outlets, nb::arg("raster"),
        nb::arg("outlet_mask"), nb::arg("min_slope") = 0.0001f, nb::arg("n_threads") = 0,
        nb::arg("priority_flood_backend") = "binary", R"doc(
Fill pits using positive cells in outlet_mask as drainage outlets.

Only the selected pixels act as outlets (the raster edge does not, unless
selected), so a DEM can be filled towards chosen edges.

Parameters
----------
raster : Raster
    Elevation raster (m).
outlet_mask : Raster
    0/1 raster on the same grid; pixels > 0.5 are outlets (for example from
    ``make_outlet_mask``). Must select at least one valid pixel.
min_slope : float, default=0.0001
    Gradient (m/m) imposed across filled flats: each filled pixel is raised
    at least ``min_slope * cell_size`` (``* sqrt(2)`` on a diagonal step)
    above the pixel it drains to, so every pixel has a downhill route. 0
    fills to a flat surface.
n_threads : int, default=0
    OpenMP thread count; 0 uses the runtime default.
priority_flood_backend : {"binary", "radix"}, default="binary"
    Priority-queue implementation: ``"binary"`` heap or ``"radix"``
    (hybrid radix heap). Both give the same result; only speed differs.

Returns
-------
Raster
    Filled DEM (m).
)doc");
    m.def("fill_pits_with_sinks", &fill_pits_with_sinks_mask, nb::arg("raster"),
        nb::arg("sinks"), nb::arg("min_slope") = 0.0001f, nb::arg("n_threads") = 0,
        nb::arg("priority_flood_backend") = "binary", R"doc(
Fill pits while retaining positive cells in sinks as declared sinks.

Depressions draining to a declared sink (moulin, swallow hole, playa) are
not filled; use the same pixels as ``forced_outlets`` of ``FlowInfo``.

Parameters
----------
raster : Raster
    Elevation raster (m).
sinks : Raster
    0/1 raster on the same grid; pixels > 0.5 are sinks.
min_slope : float, default=0.0001
    Gradient (m/m) imposed across filled flats: each filled pixel is raised
    at least ``min_slope * cell_size`` (``* sqrt(2)`` on a diagonal step)
    above the pixel it drains to, so every pixel has a downhill route. 0
    fills to a flat surface.
n_threads : int, default=0
    OpenMP thread count; 0 uses the runtime default.
priority_flood_backend : {"binary", "radix"}, default="binary"
    Priority-queue implementation: ``"binary"`` heap or ``"radix"``
    (hybrid radix heap). Both give the same result; only speed differs.

Returns
-------
Raster
    Filled DEM (m).
)doc");
    m.def("fill_pits_with_sinks", &fill_pits_with_sinks_points, nb::arg("raster"),
        nb::arg("sinks"), nb::arg("min_slope") = 0.0001f, nb::arg("n_threads") = 0,
        nb::arg("priority_flood_backend") = "binary", R"doc(
Fill pits while retaining zero-based (row, col) sinks.

Parameters
----------
raster : Raster
    Elevation raster (m).
sinks : array-like of int
    N x 2 array of ``(row, col)`` sink pixels (valid pixels).
min_slope : float, default=0.0001
    Gradient (m/m) imposed across filled flats: each filled pixel is raised
    at least ``min_slope * cell_size`` (``* sqrt(2)`` on a diagonal step)
    above the pixel it drains to, so every pixel has a downhill route. 0
    fills to a flat surface.
n_threads : int, default=0
    OpenMP thread count; 0 uses the runtime default.
priority_flood_backend : {"binary", "radix"}, default="binary"
    Priority-queue implementation: ``"binary"`` heap or ``"radix"``
    (hybrid radix heap). Both give the same result; only speed differs.

Returns
-------
Raster
    Filled DEM (m).
)doc");
    m.def("fill_pits_boundary_aware", &fill_pits_boundary_aware, nb::arg("raster"),
        nb::arg("boundary_conditions") = "oooo", nb::arg("min_slope") = 0.0001f,
        nb::arg("n_threads") = 0, nb::arg("priority_flood_backend") = "binary", R"doc(
Fill pits with NESW open, closed, periodic, or base-level boundaries.

Closed (``n``) edges are not outlets; periodic (``p``) edges wrap, so
depressions spilling across the seam are filled to their true spill level.
Fill with the same ``boundary_conditions`` you will give ``FlowInfo``.

Every valid pixel touching nodata is also an outlet, including one on a
closed edge, so DEMs with a nodata frame or voids (such as those from
``fetch_dem``) drain out where the data ends. Nodata pixels on an open edge
are skipped. Without a closed edge this is the same fill as ``fill_pits``.

Parameters
----------
raster : Raster
    Elevation raster (m).
boundary_conditions : str, default='oooo'
    Four NESW codes ``o``, ``n``, ``p``, ``b`` (see
    ``make_boundary_conditions``). Both axes cannot be periodic; with a
    periodic axis the other edges must be open.
min_slope : float, default=0.0001
    Gradient (m/m) imposed across filled flats: each filled pixel is raised
    at least ``min_slope * cell_size`` (``* sqrt(2)`` on a diagonal step)
    above the pixel it drains to, so every pixel has a downhill route. 0
    fills to a flat surface.
n_threads : int, default=0
    OpenMP thread count; 0 uses the runtime default.
priority_flood_backend : {"binary", "radix"}, default="binary"
    Priority-queue implementation: ``"binary"`` heap or ``"radix"``
    (hybrid radix heap). Both give the same result; only speed differs.

Returns
-------
Raster
    Filled DEM (m).
)doc");
    m.def("breach_least_cost", &breach_least_cost, nb::arg("raster"),
        nb::arg("max_distance") = nb::none(),
        nb::arg("max_cost") = (std::numeric_limits<float>::max)(),
        nb::arg("prefer_shortest_path") = false, nb::arg("flat_increment") = 0.0f,
        nb::arg("fill_remainder") = true, nb::arg("n_threads") = 0,
        nb::arg("edge_policy") = "outlet", nb::arg("fill_min_slope") = 0.0001f, R"doc(
Breach depressions along least-cost paths (Lindsay 2020).

A port of WhiteboxTools BreachDepressionsLeastCost (whitebox_workflows 2.0.6).

Parameters
----------
raster : Raster
    Elevation raster (m).
max_distance : int or None, default=None
    Maximum breach path length in cells, counted in steps from the pit (a
    diagonal step counts as one). 1 to 32767. None selects the edge policy's
    default: 1000 for ``outlet`` (long enough to bridge large water bodies)
    and 100 for ``whitebox`` (the WhiteboxTools default).
max_cost : float, optional
    A pit whose cheapest path costs more is left unbreached. The cost is the
    depth each path cell is lowered, summed: metres, or metres x cells when
    ``prefer_shortest_path`` is true. Default unlimited.
prefer_shortest_path : bool, default=False
    Weight each step's cost by its length in cells, favouring shorter paths.
flat_increment : float, default=0
    Elevation drop (m) per cell along a breach channel. 0 or less selects twice
    the float32 spacing at the largest absolute elevation.
fill_remainder : bool, default=True
    Fill what the breach leaves (Priority-Flood with ``fill_min_slope``) so
    every interior cell has a lower neighbour.
n_threads : int, default=0
    Ignored; the algorithm is serial.
edge_policy : {"outlet", "whitebox"}, default="outlet"
    ``outlet``: the raster edge is an outlet, like a nodata cell.
    ``whitebox``: as in WhiteboxTools, a breach must reach an edge cell that
    is already low enough.
fill_min_slope : float, default=0.0001
    Gradient (m/m) imposed by the remainder fill.
)doc");
    m.def("least_cost_carve_depth", &least_cost_carve_depth, nb::arg("raster"),
        nb::arg("max_distance") = nb::none(),
        nb::arg("max_cost") = (std::numeric_limits<float>::max)(),
        nb::arg("prefer_shortest_path") = false, nb::arg("flat_increment") = 0.0f,
        nb::arg("fill_remainder") = true, nb::arg("n_threads") = 0,
        nb::arg("edge_policy") = "outlet", nb::arg("fill_min_slope") = 0.0001f,
        R"doc(
Return the non-negative elevation removed by least-cost breaching (parameters as breach_least_cost).

Parameters
----------
raster : Raster
    Elevation raster (m).
max_distance, max_cost, prefer_shortest_path, flat_increment, fill_remainder, n_threads, edge_policy, fill_min_slope
    As in ``breach_least_cost``.

Returns
-------
Raster
    ``raster - breach_least_cost(raster, ...)`` in metres.
)doc");
}
} // namespace lsdtt3_python
