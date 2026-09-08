#include <cstdint>
#include <cctype>
#include <cstring>
#include <cmath>
#include <cstdlib>
#include <limits>
#include <random>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <nanobind/nanobind.h>
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

    Raster outlets = make_outlet_mask(
        raster, boundary[0] != 'n', boundary[1] != 'n',
        boundary[2] != 'n', boundary[3] != 'n'
    );
    return fill_pits_with_outlets(raster, outlets, min_slope, n_threads, backend);
}

Raster breach_least_cost(
    const Raster& raster,
    int max_distance,
    float max_cost,
    bool prefer_shortest_path,
    float flat_increment,
    bool fill_remainder,
    int n_threads)
{
    return Raster(lsdtt::breach_depressions_least_cost(
        raster.native(), max_distance, max_cost, prefer_shortest_path,
        flat_increment, fill_remainder, n_threads
    ));
}

Raster least_cost_carve_depth(
    const Raster& raster,
    int max_distance,
    float max_cost,
    bool prefer_shortest_path,
    float flat_increment,
    bool fill_remainder,
    int n_threads)
{
    const lsdtt::lsdtt_raster_cpu breached = lsdtt::breach_depressions_least_cost(
        raster.native(), max_distance, max_cost, prefer_shortest_path,
        flat_increment, fill_remainder, n_threads
    );
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
        nb::arg("priority_flood_backend") = "binary",
        "Fill depressions while imposing an optional minimum drainage slope.");
    m.def("fill_pits_priority_flood", &fill_pits_priority_flood,
        nb::arg("raster"), nb::arg("min_slope") = 0.0001f,
        nb::arg("n_threads") = 0, nb::arg("use_radix_heap") = false,
        "Explicit Improved Priority-Flood; optionally use its radix-heap queue.");
    m.def("fill_depth", &fill_depth, nb::arg("raster"),
        nb::arg("min_slope") = 0.0001f, nb::arg("n_threads") = 0,
        "Return the non-negative elevation added by pit filling.");
    m.def("breach_depressions", &breach_depressions, nb::arg("raster"),
        nb::arg("n_threads") = 0,
        "Breach depressions by carving drainage paths.");
    m.def("carve_depth", &carve_depth, nb::arg("raster"), nb::arg("n_threads") = 0,
        "Return the non-negative elevation removed by depression breaching.");
    m.def("make_boundary_conditions", &make_boundary_conditions,
        nb::arg("north") = "open", nb::arg("east") = "open",
        nb::arg("south") = "open", nb::arg("west") = "open", R"doc(
Create an LSDTopoTools NESW boundary-condition string.

Parameters
----------
north, east, south, west : {"open", "closed", "periodic", "base_level"}
    Boundary treatment. The returned string has NESW order.
)doc");
    m.def("make_outlet_mask", &make_outlet_mask, nb::arg("raster"),
        nb::arg("north") = false, nb::arg("east") = false,
        nb::arg("south") = false, nb::arg("west") = false,
        "Create a valid binary mask with selected raster edges set to one.");
    m.def("mask_from_points", &mask_from_points, nb::arg("raster"), nb::arg("points"), R"doc(
Create a binary mask from zero-based ``(row, col)`` point coordinates.

Parameters
----------
points : array-like, shape (n, 2)
    Row-major NumPy coordinates. Points must be valid raster cells.
)doc");
    m.def("fill_pits_with_outlets", &fill_pits_with_outlets, nb::arg("raster"),
        nb::arg("outlet_mask"), nb::arg("min_slope") = 0.0001f, nb::arg("n_threads") = 0,
        nb::arg("priority_flood_backend") = "binary",
        "Fill pits using positive cells in outlet_mask as drainage outlets.");
    m.def("fill_pits_with_sinks", &fill_pits_with_sinks_mask, nb::arg("raster"),
        nb::arg("sinks"), nb::arg("min_slope") = 0.0001f, nb::arg("n_threads") = 0,
        nb::arg("priority_flood_backend") = "binary",
        "Fill pits while retaining positive cells in sinks as declared sinks.");
    m.def("fill_pits_with_sinks", &fill_pits_with_sinks_points, nb::arg("raster"),
        nb::arg("sinks"), nb::arg("min_slope") = 0.0001f, nb::arg("n_threads") = 0,
        nb::arg("priority_flood_backend") = "binary",
        "Fill pits while retaining zero-based (row, col) sinks.");
    m.def("fill_pits_boundary_aware", &fill_pits_boundary_aware, nb::arg("raster"),
        nb::arg("boundary_conditions") = "oooo", nb::arg("min_slope") = 0.0001f,
        nb::arg("n_threads") = 0, nb::arg("priority_flood_backend") = "binary",
        "Fill pits with NESW open, closed, periodic, or base-level boundaries.");
    m.def("breach_least_cost", &breach_least_cost, nb::arg("raster"),
        nb::arg("max_distance") = 1000,
        nb::arg("max_cost") = (std::numeric_limits<float>::max)(),
        nb::arg("prefer_shortest_path") = true, nb::arg("flat_increment") = 0.0f,
        nb::arg("fill_remainder") = true, nb::arg("n_threads") = 0, R"doc(
Breach depressions using configurable least-cost paths.

Parameters
----------
max_distance : int, default=1000
    Maximum search distance in cells.
max_cost : float, optional
    Maximum accepted path cost.
prefer_shortest_path, fill_remainder : bool, optional
    Control path selection and whether unbreached pits are filled.
flat_increment : float, optional
    Elevation decrement used to route across flats.
)doc");
    m.def("least_cost_carve_depth", &least_cost_carve_depth, nb::arg("raster"),
        nb::arg("max_distance") = 1000,
        nb::arg("max_cost") = (std::numeric_limits<float>::max)(),
        nb::arg("prefer_shortest_path") = true, nb::arg("flat_increment") = 0.0f,
        nb::arg("fill_remainder") = true, nb::arg("n_threads") = 0,
        "Return the non-negative elevation removed by least-cost breaching.");
}
} // namespace lsdtt3_python
