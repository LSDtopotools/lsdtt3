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

namespace nb = nanobind;

namespace {

class Raster
{
public:
    explicit Raster(lsdtt::lsdtt_raster_cpu raster)
        : raster_(std::move(raster))
    {
    }

    static Raster read(const std::string& path)
    {
        lsdtt::lsdtt_data_loader loader;
        loader.open_raster(path);
        return Raster(lsdtt::lsdtt_raster_cpu(loader));
    }

    void write(const std::string& path, const std::string& format) const
    {
        lsdtt::lsdtt_data_loader loader;
        loader.write_full(
            raster_.get_data(), raster_.get_meta(), path,
            lsdtt::resolve_raster_output_format(format)
        );
    }

    const lsdtt::lsdtt_raster_cpu& native() const noexcept
    {
        return raster_;
    }

private:
    lsdtt::lsdtt_raster_cpu raster_;
};

nb::object as_numpy(const std::vector<float>& data, const lsdtt::lsdtt_raster_meta& meta)
{
    nb::object numpy = nb::module_::import_("numpy");
    return numpy.attr("array")(data, nb::arg("dtype") = "float32")
        .attr("reshape")(nb::make_tuple(meta.n_rows, meta.n_cols));
}

nb::object mask_as_numpy(const Raster& raster)
{
    nb::object numpy = nb::module_::import_("numpy");
    const auto& meta = raster.native().get_meta();
    return numpy.attr("array")(raster.native().get_mask(), nb::arg("dtype") = "bool")
        .attr("reshape")(nb::make_tuple(meta.n_rows, meta.n_cols));
}

Raster add_uniform_noise(const Raster& raster, std::uint64_t seed, float magnitude)
{
    if (magnitude < 0.0f)
    {
        throw std::invalid_argument("magnitude must be non-negative");
    }

    const lsdtt::lsdtt_raster_cpu& source = raster.native();
    std::vector<float> data = source.get_data();
    std::vector<std::uint8_t> mask = source.get_mask();
    std::mt19937_64 generator(seed);
    std::uniform_real_distribution<float> distribution(-magnitude, magnitude);

    for (std::size_t i = 0; i < data.size(); ++i)
    {
        if (mask[i] != 0)
        {
            data[i] += distribution(generator);
        }
    }

    return Raster(lsdtt::lsdtt_raster_cpu(
        std::move(data), std::move(mask), source.get_meta()
    ));
}

Raster set_nodata_rectangle(
    const Raster& raster,
    int row_start,
    int row_stop,
    int col_start,
    int col_stop)
{
    const lsdtt::lsdtt_raster_cpu& source = raster.native();
    const auto& meta = source.get_meta();
    if (row_start < 0 || row_stop > meta.n_rows || row_start >= row_stop ||
        col_start < 0 || col_stop > meta.n_cols || col_start >= col_stop)
    {
        throw std::invalid_argument(
            "rectangle must be non-empty and lie within the raster"
        );
    }

    std::vector<float> data = source.get_data();
    std::vector<std::uint8_t> mask = source.get_mask();
    const float nodata = static_cast<float>(meta.no_data_value);
    for (int row = row_start; row < row_stop; ++row)
    {
        for (int col = col_start; col < col_stop; ++col)
        {
            const std::size_t index = static_cast<std::size_t>(row) * meta.n_cols + col;
            data[index] = nodata;
            mask[index] = 0;
        }
    }

    return Raster(lsdtt::lsdtt_raster_cpu(
        std::move(data), std::move(mask), meta
    ));
}

lsdtt::lsdtt_raster_meta make_raster_meta(
    int n_rows,
    int n_cols,
    double x_min,
    double y_min,
    double cell_size,
    double nodata_value,
    std::string crs_wkt,
    int epsg_code)
{
    lsdtt::lsdtt_raster_meta metadata;
    metadata.n_rows = n_rows;
    metadata.n_cols = n_cols;
    metadata.x_min = x_min;
    metadata.y_min = y_min;
    metadata.cell_size = cell_size;
    metadata.no_data_value = nodata_value;
    metadata.crs_wkt = std::move(crs_wkt);
    metadata.epsg_code = epsg_code;
    if (!metadata.is_valid())
    {
        throw std::invalid_argument(
            "n_rows, n_cols, and cell_size must be positive"
        );
    }
    return metadata;
}

Raster raster_from_numpy(
    nb::object values,
    nb::object mask,
    const lsdtt::lsdtt_raster_meta& metadata)
{
    if (!metadata.is_valid())
    {
        throw std::invalid_argument("metadata must describe a non-empty raster");
    }

    nb::object numpy = nb::module_::import_("numpy");
    nb::object data_array = numpy.attr("asarray")(
        values, nb::arg("dtype") = "float32"
    );
    nb::object mask_array = numpy.attr("asarray")(
        mask, nb::arg("dtype") = "bool"
    );

    const auto shape_matches_metadata = [&metadata](const nb::object& array) {
        if (nb::cast<int>(array.attr("ndim")) != 2)
        {
            return false;
        }
        nb::tuple shape = nb::cast<nb::tuple>(array.attr("shape"));
        return nb::cast<int>(shape[0]) == metadata.n_rows &&
               nb::cast<int>(shape[1]) == metadata.n_cols;
    };
    if (!shape_matches_metadata(data_array) || !shape_matches_metadata(mask_array))
    {
        throw std::invalid_argument(
            "data and mask must be two-dimensional and match metadata dimensions"
        );
    }

    std::vector<float> data = nb::cast<std::vector<float>>(
        data_array.attr("ravel")()
    );
    std::vector<std::uint8_t> valid = nb::cast<std::vector<std::uint8_t>>(
        mask_array.attr("astype")("uint8").attr("ravel")()
    );
    const float nodata = static_cast<float>(metadata.no_data_value);
    for (std::size_t i = 0; i < data.size(); ++i)
    {
        if (valid[i] == 0)
        {
            data[i] = nodata;
        }
    }

    return Raster(lsdtt::lsdtt_raster_cpu(
        std::move(data), std::move(valid), metadata
    ));
}

Raster fill_interior_nodata(
    const Raster& raster,
    const std::string& method,
    int max_iter,
    float tolerance,
    int n_threads)
{
    lsdtt::nodata_fill_method fill_method;
    if (method == "laplacian")
    {
        fill_method = lsdtt::nodata_fill_method::laplacian;
    }
    else if (method == "iterative")
    {
        fill_method = lsdtt::nodata_fill_method::iterative_weighted_average;
    }
    else
    {
        throw std::invalid_argument(
            "method must be 'laplacian' or 'iterative'"
        );
    }

    if (max_iter <= 0)
    {
        throw std::invalid_argument("max_iter must be positive");
    }
    if (tolerance <= 0.0f)
    {
        throw std::invalid_argument("tolerance must be positive");
    }

    return Raster(lsdtt::fill_interior_nodata(
        raster.native(), fill_method, max_iter, tolerance, n_threads
    ));
}

Raster clip_to_maximum_rectangle(const Raster& raster)
{
    return Raster(lsdtt::clip_to_maximum_rectangle(raster.native()));
}

Raster mask_to_nodata(const Raster& raster, float threshold, bool below)
{
    return Raster(lsdtt::mask_to_nodata_using_threshold(
        raster.native(), threshold, below
    ));
}

Raster threshold_above(const Raster& raster, float threshold)
{
    return Raster(lsdtt::algebra::raster_binary_true_above_threshold(
        raster.native(), threshold
    ));
}

Raster threshold_below(const Raster& raster, float threshold)
{
    return Raster(lsdtt::algebra::raster_binary_true_below_threshold(
        raster.native(), threshold
    ));
}

Raster invert_binary_mask(const Raster& raster)
{
    return Raster(lsdtt::algebra::invert_binary_mask(raster.native()));
}

Raster mask_by_binary(const Raster& values, const Raster& keep_mask)
{
    return Raster(lsdtt::algebra::raster_mask_by_binary_keep(
        values.native(), keep_mask.native()
    ));
}

Raster add(const Raster& a, const Raster& b)
{
    return Raster(lsdtt::algebra::raster_add(a.native(), b.native()));
}

Raster subtract(const Raster& a, const Raster& b)
{
    return Raster(lsdtt::algebra::raster_subtract(a.native(), b.native()));
}

Raster multiply(const Raster& a, const Raster& b)
{
    return Raster(lsdtt::algebra::raster_multiply(a.native(), b.native()));
}

Raster divide(const Raster& a, const Raster& b)
{
    return Raster(lsdtt::algebra::raster_divide(a.native(), b.native()));
}

Raster scale(const Raster& raster, float factor)
{
    return Raster(lsdtt::algebra::raster_scale(raster.native(), factor));
}

Raster normalize_unit_interval(const Raster& raster)
{
    return Raster(lsdtt::algebra::raster_normalise_unit_interval(raster.native()));
}

Raster fill_pits(const Raster& raster, float min_slope, int n_threads)
{
    return Raster(lsdtt::fill_pits(raster.native(), min_slope, n_threads));
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

Raster hillshade(const Raster& raster, double azimuth, double altitude, double z_factor)
{
    return Raster(lsdtt::compute_hillshade(
        raster.native(), azimuth, altitude, z_factor
    ));
}

Raster multidirectional_hillshade(
    const Raster& raster,
    int n_directions,
    double altitude,
    double z_factor,
    double start_azimuth)
{
    return Raster(lsdtt::compute_multidirectional_hillshade(
        raster.native(), n_directions, altitude, z_factor, start_azimuth
    ));
}

Raster slope(const Raster& raster)
{
    return Raster(lsdtt::compute_gradient_magnitude_3x3_square(raster.native()));
}

Raster plan_curvature(const Raster& raster)
{
    return Raster(lsdtt::compute_plan_curvature_3x3_square(raster.native()));
}

Raster total_curvature(const Raster& raster, double radius_metres)
{
    return Raster(lsdtt::moving_window_total_curvature(
        raster.native(), radius_metres
    ));
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
    int n_threads)
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
        raster.native(), min_slope, {rows}, {cols}, {"outlet"}, "outlet", n_threads
    ));
}

Raster fill_pits_with_sinks_mask(
    const Raster& raster,
    const Raster& sink_mask,
    float min_slope,
    int n_threads)
{
    return Raster(lsdtt::fill_pits_with_sinks(
        raster.native(), selected_flat_indices(raster, sink_mask, "sink_mask"), min_slope, n_threads
    ));
}

Raster fill_pits_with_sinks_points(
    const Raster& raster,
    nb::object points,
    float min_slope,
    int n_threads)
{
    return Raster(lsdtt::fill_pits_with_sinks(
        raster.native(), row_col_points_to_flat_indices(raster, points, "sinks"), min_slope, n_threads
    ));
}

Raster fill_pits_boundary_aware(
    const Raster& raster,
    const std::string& boundary_conditions,
    float min_slope,
    int n_threads)
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
        return Raster(lsdtt::fill_pits(raster.native(), min_slope, boundary, n_threads));
    }

    Raster outlets = make_outlet_mask(
        raster, boundary[0] != 'n', boundary[1] != 'n',
        boundary[2] != 'n', boundary[3] != 'n'
    );
    return fill_pits_with_outlets(raster, outlets, min_slope, n_threads);
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

GDALDataset* make_gdal_dataset(const Raster& raster)
{
    GDALAllRegister();
    const auto& meta = raster.native().get_meta();
    GDALDriver* driver = GetGDALDriverManager()->GetDriverByName("MEM");
    GDALDataset* dataset = driver ? driver->Create("", meta.n_cols, meta.n_rows, 1, GDT_Float32, nullptr) : nullptr;
    if (!dataset) throw std::runtime_error("could not create GDAL in-memory raster");
    double transform[6] = {meta.x_min, meta.cell_size, 0.0, meta.y_max(), 0.0, -meta.cell_size};
    dataset->SetGeoTransform(transform);
    dataset->SetProjection(meta.crs_wkt.c_str());
    GDALRasterBand* band = dataset->GetRasterBand(1);
    band->SetNoDataValue(meta.no_data_value);
    if (band->RasterIO(GF_Write, 0, 0, meta.n_cols, meta.n_rows,
                       const_cast<float*>(raster.native().get_data().data()), meta.n_cols, meta.n_rows,
                       GDT_Float32, 0, 0, nullptr) != CE_None)
    {
        GDALClose(dataset);
        throw std::runtime_error("could not populate GDAL in-memory raster");
    }
    return dataset;
}

Raster raster_from_gdal_dataset(GDALDataset* dataset)
{
    if (!dataset || dataset->GetRasterCount() < 1) throw std::runtime_error("GDAL warp produced no raster");
    const int ncols = dataset->GetRasterXSize();
    const int nrows = dataset->GetRasterYSize();
    double transform[6];
    if (dataset->GetGeoTransform(transform) != CE_None || transform[1] <= 0.0 || transform[5] >= 0.0)
        throw std::runtime_error("warped raster has an unsupported geotransform");
    lsdtt::lsdtt_raster_meta meta;
    meta.n_rows = nrows; meta.n_cols = ncols; meta.x_min = transform[0];
    meta.y_min = transform[3] + nrows * transform[5]; meta.cell_size = transform[1];
    meta.crs_wkt = dataset->GetProjectionRef();
    GDALRasterBand* band = dataset->GetRasterBand(1);
    int has_nodata = FALSE;
    meta.no_data_value = band->GetNoDataValue(&has_nodata);
    if (!has_nodata) meta.no_data_value = -9999.0;
    std::vector<float> data(static_cast<std::size_t>(meta.n_pixels()));
    std::vector<std::uint8_t> mask(data.size());
    if (band->RasterIO(GF_Read, 0, 0, ncols, nrows, data.data(), ncols, nrows, GDT_Float32, 0, 0, nullptr) != CE_None)
        throw std::runtime_error("could not read warped raster");
    std::vector<std::uint8_t> byte_mask(mask.size());
    if (band->GetMaskBand()->RasterIO(GF_Read, 0, 0, ncols, nrows, byte_mask.data(), ncols, nrows, GDT_Byte, 0, 0, nullptr) != CE_None)
        throw std::runtime_error("could not read warped raster mask");
    for (std::size_t i = 0; i < mask.size(); ++i) {
        mask[i] = byte_mask[i] ? 1 : 0;
        if (!mask[i]) data[i] = static_cast<float>(meta.no_data_value);
    }
    OGRSpatialReference srs;
    if (!meta.crs_wkt.empty() && srs.SetFromUserInput(meta.crs_wkt.c_str()) == OGRERR_NONE) {
        const char* code = srs.GetAuthorityCode(nullptr);
        if (code) meta.epsg_code = std::atoi(code);
    }
    return Raster(lsdtt::lsdtt_raster_cpu(std::move(data), std::move(mask), meta));
}

Raster gdal_warp(const Raster& raster, const std::vector<std::string>& arguments)
{
    GDALDataset* source = make_gdal_dataset(raster);
    static std::uint64_t counter = 0;
    const std::string destination_path = "/vsimem/lsdtt3-warp-" + std::to_string(++counter) + ".tif";
    std::vector<std::string> all_args = {"-of", "GTiff"};
    all_args.insert(all_args.end(), arguments.begin(), arguments.end());
    std::vector<char*> argv;
    for (auto& argument : all_args) argv.push_back(argument.data());
    argv.push_back(nullptr);
    GDALWarpAppOptions* options = GDALWarpAppOptionsNew(argv.data(), nullptr);
    int usage_error = FALSE;
    GDALDatasetH source_handle = source;
    GDALDatasetH destination = GDALWarp(destination_path.c_str(), nullptr, 1, &source_handle, options, &usage_error);
    GDALWarpAppOptionsFree(options);
    if (!destination || usage_error) {
        GDALClose(source);
        VSIUnlink(destination_path.c_str());
        throw std::runtime_error("GDAL raster warp failed");
    }
    try {
        Raster output = raster_from_gdal_dataset(static_cast<GDALDataset*>(destination));
        GDALClose(destination);
        GDALClose(source);
        VSIUnlink(destination_path.c_str());
        return output;
    } catch (...) {
        GDALClose(destination);
        GDALClose(source);
        VSIUnlink(destination_path.c_str());
        throw;
    }
}

Raster reproject_to_grid(const Raster& raster, const lsdtt::lsdtt_raster_meta& target, const std::string& resampling)
{
    if (!target.is_valid() || target.crs_wkt.empty()) throw std::invalid_argument("target metadata must define a valid CRS grid");
    return gdal_warp(raster, {"-t_srs", target.crs_wkt, "-te", std::to_string(target.x_min),
        std::to_string(target.y_min), std::to_string(target.x_max()), std::to_string(target.y_max()),
        "-ts", std::to_string(target.n_cols), std::to_string(target.n_rows), "-r", resampling,
        "-dstnodata", std::to_string(target.no_data_value)});
}

Raster reproject_to_epsg(const Raster& raster, int epsg, double cell_size, const std::string& resampling)
{
    OGRSpatialReference srs;
    if (srs.importFromEPSG(epsg) != OGRERR_NONE) throw std::invalid_argument("invalid target EPSG code");
    char* wkt = nullptr;
    srs.exportToWkt(&wkt);
    std::vector<std::string> arguments = {"-t_srs", wkt ? wkt : "", "-r", resampling};
    CPLFree(wkt);
    if (cell_size > 0.0) {
        arguments.push_back("-tr"); arguments.push_back(std::to_string(cell_size)); arguments.push_back(std::to_string(cell_size));
    }
    return gdal_warp(raster, arguments);
}

Raster mosaic(const std::vector<Raster>& rasters)
{
    if (rasters.empty()) throw std::invalid_argument("rasters must not be empty");
    const auto& first = rasters.front().native().get_meta();
    double xmin = first.x_min, ymin = first.y_min, xmax = first.x_max(), ymax = first.y_max();
    for (const auto& raster : rasters) {
        const auto& meta = raster.native().get_meta();
        if (meta.cell_size != first.cell_size || meta.crs_wkt != first.crs_wkt)
            throw std::invalid_argument("all rasters must have matching CRS and cell size");
        const double x_offset = (meta.x_min - first.x_min) / first.cell_size;
        const double y_offset = (meta.y_min - first.y_min) / first.cell_size;
        if (std::abs(x_offset - std::round(x_offset)) > 1e-9 ||
            std::abs(y_offset - std::round(y_offset)) > 1e-9)
            throw std::invalid_argument("all rasters must be aligned on the same grid");
        xmin = std::min(xmin, meta.x_min); ymin = std::min(ymin, meta.y_min);
        xmax = std::max(xmax, meta.x_max()); ymax = std::max(ymax, meta.y_max());
    }
    lsdtt::lsdtt_raster_meta meta = first;
    meta.x_min = xmin; meta.y_min = ymin;
    meta.n_cols = static_cast<int>(std::llround((xmax - xmin) / meta.cell_size));
    meta.n_rows = static_cast<int>(std::llround((ymax - ymin) / meta.cell_size));
    std::vector<float> data(static_cast<std::size_t>(meta.n_pixels()), static_cast<float>(meta.no_data_value));
    std::vector<std::uint8_t> mask(data.size(), 0);
    for (const auto& raster : rasters) {
        const auto& source = raster.native(); const auto& source_meta = source.get_meta();
        const int col_offset = static_cast<int>(std::llround((source_meta.x_min - xmin) / meta.cell_size));
        const int row_offset = static_cast<int>(std::llround((ymax - source_meta.y_max()) / meta.cell_size));
        for (int row = 0; row < source_meta.n_rows; ++row) for (int col = 0; col < source_meta.n_cols; ++col) {
            const std::size_t from = static_cast<std::size_t>(row) * source_meta.n_cols + col;
            if (!source.get_mask()[from]) continue;
            const std::size_t to = static_cast<std::size_t>(row + row_offset) * meta.n_cols + col + col_offset;
            data[to] = source.get_data()[from]; mask[to] = 1;
        }
    }
    return Raster(lsdtt::lsdtt_raster_cpu(std::move(data), std::move(mask), meta));
}

std::string geodataframe_geojson(const Raster& raster, nb::object geodataframe)
{
    nb::object crs = geodataframe.attr("crs");
    if (crs.is_none()) throw std::invalid_argument("GeoDataFrame must define a CRS");
    const std::string vector_wkt = nb::cast<std::string>(crs.attr("to_wkt")());
    const std::string& raster_wkt = raster.native().get_meta().crs_wkt;
    OGRSpatialReference vector_srs, raster_srs;
    if (raster_wkt.empty() || vector_srs.SetFromUserInput(vector_wkt.c_str()) != OGRERR_NONE ||
        raster_srs.SetFromUserInput(raster_wkt.c_str()) != OGRERR_NONE || !vector_srs.IsSame(&raster_srs))
        throw std::invalid_argument("GeoDataFrame CRS must match the raster CRS; reproject explicitly first");
    const std::string geojson = nb::cast<std::string>(geodataframe.attr("to_json")());
    static std::uint64_t counter = 0;
    const std::string path = "/vsimem/lsdtt3-geodataframe-" + std::to_string(++counter) + ".geojson";
    auto* bytes = static_cast<GByte*>(CPLMalloc(geojson.size()));
    std::memcpy(bytes, geojson.data(), geojson.size());
    if (!VSIFileFromMemBuffer(path.c_str(), bytes, static_cast<vsi_l_offset>(geojson.size()), TRUE)) {
        CPLFree(bytes); throw std::runtime_error("could not create in-memory GeoJSON");
    }
    return path;
}

lsdtt::polygon_rasterise_rule polygon_rule(const std::string& rule)
{
    if (rule == "any" || rule == "any_coverage") return lsdtt::polygon_rasterise_rule::any_coverage;
    if (rule == "majority" || rule == "majority_coverage") return lsdtt::polygon_rasterise_rule::majority_coverage;
    if (rule == "full" || rule == "full_coverage") return lsdtt::polygon_rasterise_rule::full_coverage;
    throw std::invalid_argument("coverage must be any, majority, or full");
}

Raster rasterize_polygons(const Raster& raster, nb::object polygons, const std::string& coverage)
{
    const std::string path = geodataframe_geojson(raster, polygons);
    try {
        Raster result(lsdtt::rasterise_polygons(raster.native(), path, polygon_rule(coverage)));
        VSIUnlink(path.c_str());
        return result;
    } catch (...) { VSIUnlink(path.c_str()); throw; }
}

Raster clip_to_polygons(const Raster& raster, nb::object polygons, const std::string& coverage, bool crop)
{
    nb::object bounds_object = polygons.attr("total_bounds");
    std::vector<double> bounds = nb::cast<std::vector<double>>(bounds_object);
    if (bounds.size() != 4) throw std::invalid_argument("GeoDataFrame has invalid bounds");
    Raster base = crop ? Raster(lsdtt::clip_raster_to_bbox(raster.native(), bounds[0], bounds[1], bounds[2], bounds[3])) : Raster(raster.native());
    Raster polygon_mask = rasterize_polygons(base, polygons, coverage);
    return Raster(lsdtt::algebra::raster_mask_by_binary_keep(base.native(), polygon_mask.native()));
}

nb::object sample_points(const Raster& raster, nb::object points, const std::string& column)
{
    const std::string path = geodataframe_geojson(raster, points);
    VSIUnlink(path.c_str()); // CRS validation only; point coordinates are read directly below.
    nb::object geometry = points.attr("geometry");
    std::vector<double> xs = nb::cast<std::vector<double>>(geometry.attr("x").attr("tolist")());
    std::vector<double> ys = nb::cast<std::vector<double>>(geometry.attr("y").attr("tolist")());
    const auto& source = raster.native(); const auto& meta = source.get_meta();
    nb::list values;
    for (std::size_t i = 0; i < xs.size(); ++i) {
        const int col = static_cast<int>(std::floor((xs[i] - meta.x_min) / meta.cell_size));
        const int row = static_cast<int>(std::floor((meta.y_max() - ys[i]) / meta.cell_size));
        if (row < 0 || row >= meta.n_rows || col < 0 || col >= meta.n_cols) { values.append(nb::none()); continue; }
        const std::size_t flat = static_cast<std::size_t>(row) * meta.n_cols + col;
        if (!source.get_mask()[flat]) values.append(nb::none()); else values.append(source.get_data()[flat]);
    }
    nb::object output = points.attr("copy")();
    output.attr("__setitem__")(column, values);
    return output;
}

Raster aspect(
    const Raster& raster,
    double lambda_metres,
    int n_threads,
    double flat_slope_tolerance)
{
    const double lambda = lambda_metres > 0.0
        ? lambda_metres
        : 3.0 * raster.native().get_meta().cell_size;
    const lsdtt::poly_fit_coeffs coefficients = lsdtt::direct_polyfit(
        raster.native(), lambda, n_threads
    );
    lsdtt::polyfit_metric_request request;
    request.metrics = {lsdtt::polyfit_metric::aspect_degrees};
    request.n_threads = n_threads;
    request.flat_slope_tolerance = flat_slope_tolerance;
    return Raster(lsdtt::compute_polyfit_metrics(coefficients, request).at(
        lsdtt::polyfit_metric::aspect_degrees
    ));
}

} // namespace

NB_MODULE(_core, module)
{
    module.doc() = "Python bindings for lsdtt3 raster operations.";

    nb::class_<lsdtt::lsdtt_raster_meta>(module, "RasterMeta", R"doc(
Metadata describing a north-up, square-cell raster grid.

Use ``RasterMeta.from_shape`` to construct an instance.
)doc")
        .def_static("from_shape", &make_raster_meta,
            nb::arg("n_rows"), nb::arg("n_cols"), nb::arg("x_min") = 0.0,
            nb::arg("y_min") = 0.0, nb::arg("cell_size") = 1.0,
            nb::arg("nodata_value") = -9999.0, nb::arg("crs_wkt") = "",
            nb::arg("epsg_code") = 0, R"doc(
Create raster metadata from a grid shape and lower-left origin.

Parameters
----------
n_rows, n_cols : int
    Positive grid dimensions.
x_min, y_min, cell_size : float
    Lower-left origin and square cell size.
nodata_value : float, optional
    Value stored for invalid pixels.
crs_wkt : str, optional
    Coordinate reference system in WKT form.
epsg_code : int, optional
    EPSG identifier when known.
)doc")
        .def_prop_ro("n_rows", [](const lsdtt::lsdtt_raster_meta& meta) { return meta.n_rows; })
        .def_prop_ro("n_cols", [](const lsdtt::lsdtt_raster_meta& meta) { return meta.n_cols; })
        .def_prop_ro("x_min", [](const lsdtt::lsdtt_raster_meta& meta) { return meta.x_min; })
        .def_prop_ro("x_max", &lsdtt::lsdtt_raster_meta::x_max)
        .def_prop_ro("y_min", [](const lsdtt::lsdtt_raster_meta& meta) { return meta.y_min; })
        .def_prop_ro("y_max", &lsdtt::lsdtt_raster_meta::y_max)
        .def_prop_ro("cell_size", [](const lsdtt::lsdtt_raster_meta& meta) { return meta.cell_size; })
        .def_prop_ro("nodata_value", [](const lsdtt::lsdtt_raster_meta& meta) { return meta.no_data_value; })
        .def_prop_ro("crs_wkt", [](const lsdtt::lsdtt_raster_meta& meta) { return meta.crs_wkt; })
        .def_prop_ro("epsg_code", [](const lsdtt::lsdtt_raster_meta& meta) { return meta.epsg_code; })
        .def_prop_ro("n_pixels", &lsdtt::lsdtt_raster_meta::n_pixels)
        .def_prop_ro("pixel_area", &lsdtt::lsdtt_raster_meta::pixel_area)
        .def_prop_ro("is_valid", &lsdtt::lsdtt_raster_meta::is_valid);

    nb::class_<Raster>(module, "Raster", R"doc(
An in-memory float raster with an explicit validity mask.

Instances are value objects: functions in this module return new rasters and
do not modify their inputs.
)doc")
        .def_static("read", &Raster::read, nb::arg("path"), R"doc(
Read a single-band GDAL-compatible raster.

Parameters
----------
path : path-like
    Raster file to load.

Returns
-------
Raster
    Raster values, metadata, and nodata mask.
)doc")
        .def_static(
            "from_numpy", &raster_from_numpy,
            nb::arg("data"), nb::arg("mask"), nb::arg("metadata"), R"doc(
Create a raster from copied two-dimensional NumPy-compatible arrays.

Parameters
----------
data : array-like of float
    Values with shape ``metadata.n_rows, metadata.n_cols``.
mask : array-like of bool
    True for valid pixels. False pixels are stored as ``nodata_value``.
metadata : RasterMeta
    Grid metadata matching both array dimensions.
)doc")
        .def("write", &Raster::write, nb::arg("path"), nb::arg("format") = "cog", R"doc(
Write this raster using an LSDTopoTools GDAL output format.

Parameters
----------
path : path-like
    Destination file.
format : str, default="cog"
    Output format, for example ``"cog"`` or ``"gtiff"``.
)doc")
        .def_prop_ro("metadata", [](const Raster& raster) -> const lsdtt::lsdtt_raster_meta& {
            return raster.native().get_meta();
        }, nb::rv_policy::reference_internal)
        .def_prop_ro("shape", [](const Raster& raster) {
            const auto& meta = raster.native().get_meta();
            return nb::make_tuple(meta.n_rows, meta.n_cols);
        })
        .def("to_numpy", [](const Raster& raster) {
            return as_numpy(raster.native().get_data(), raster.native().get_meta());
        }, R"doc(
Return a copied ``float32`` array of raster values.

Invalid cells contain the raster nodata value; use ``mask_numpy`` to identify
them.
)doc")
        .def("mask_numpy", &mask_as_numpy, R"doc(
Return a copied Boolean array, True where raster values are valid.
)doc")
        .def_prop_ro("valid_pixel_count", [](const Raster& raster) {
            return raster.native().count_valid_pixels();
        });

    module.def(
        "add_uniform_noise", &add_uniform_noise,
        nb::arg("raster"), nb::arg("seed"), nb::arg("magnitude"),
        R"doc(
Add reproducible uniform noise to valid pixels.

Parameters
----------
raster : Raster
seed : int
    Random-generator seed.
magnitude : float
    Maximum absolute perturbation. Must be non-negative.

Returns
-------
Raster
    New raster; nodata pixels are unchanged.
)doc"
    );
    module.def(
        "fill_interior_nodata", &fill_interior_nodata,
        nb::arg("raster"), nb::arg("method") = "laplacian",
        nb::arg("max_iter") = 10000, nb::arg("tolerance") = 1e-6f,
        nb::arg("n_threads") = 0,
        R"doc(
Fill nodata regions that are not connected to the raster edge.

Parameters
----------
method : {"laplacian", "iterative"}, default="laplacian"
    Interpolation method.
max_iter : int, default=10000
    Maximum solver iterations.
tolerance : float, default=1e-6
    Solver convergence tolerance.

Returns
-------
Raster
    New raster. Edge-connected nodata remains invalid.
)doc"
    );
    module.def(
        "clip_to_maximum_rectangle", &clip_to_maximum_rectangle,
        nb::arg("raster"),
        R"doc(
Crop to the largest axis-aligned rectangle containing only valid pixels.

Returns
-------
Raster
    New raster with adjusted metadata.
)doc"
    );
    module.def(
        "set_nodata_rectangle", &set_nodata_rectangle,
        nb::arg("raster"), nb::arg("row_start"), nb::arg("row_stop"),
        nb::arg("col_start"), nb::arg("col_stop"),
        R"doc(
Set a half-open row/column rectangle to nodata.

Parameters
----------
row_start, row_stop, col_start, col_stop : int
    Bounds following NumPy slicing semantics.
)doc"
    );
    module.def("mask_to_nodata", &mask_to_nodata, nb::arg("raster"),
        nb::arg("threshold"), nb::arg("below") = true,
        "Set valid pixels below (or above) a threshold to nodata.");
    module.def("threshold_above", &threshold_above, nb::arg("raster"), nb::arg("threshold"),
        "Return a binary mask that is True where valid values exceed threshold.");
    module.def("threshold_below", &threshold_below, nb::arg("raster"), nb::arg("threshold"),
        "Return a binary mask that is True where valid values are below threshold.");
    module.def("invert_binary_mask", &invert_binary_mask, nb::arg("raster"),
        "Invert a binary raster mask while preserving nodata.");
    module.def("mask_by_binary", &mask_by_binary, nb::arg("values"), nb::arg("keep_mask"),
        "Set values to nodata where keep_mask is false or nodata.");
    module.def("add", &add, nb::arg("a"), nb::arg("b"), "Add two grid-compatible rasters.");
    module.def("subtract", &subtract, nb::arg("a"), nb::arg("b"), "Subtract two grid-compatible rasters.");
    module.def("multiply", &multiply, nb::arg("a"), nb::arg("b"), "Multiply two grid-compatible rasters.");
    module.def("divide", &divide, nb::arg("a"), nb::arg("b"), "Divide two grid-compatible rasters.");
    module.def("scale", &scale, nb::arg("raster"), nb::arg("factor"),
        "Multiply valid raster values by a scalar factor.");
    module.def("normalize_unit_interval", &normalize_unit_interval, nb::arg("raster"),
        "Linearly rescale valid values to the closed interval [0, 1].");
    module.def("fill_pits", &fill_pits, nb::arg("raster"),
        nb::arg("min_slope") = 0.0001f, nb::arg("n_threads") = 0,
        "Fill depressions while imposing an optional minimum drainage slope.");
    module.def("fill_depth", &fill_depth, nb::arg("raster"),
        nb::arg("min_slope") = 0.0001f, nb::arg("n_threads") = 0,
        "Return the non-negative elevation added by pit filling.");
    module.def("breach_depressions", &breach_depressions, nb::arg("raster"),
        nb::arg("n_threads") = 0,
        "Breach depressions by carving drainage paths.");
    module.def("carve_depth", &carve_depth, nb::arg("raster"), nb::arg("n_threads") = 0,
        "Return the non-negative elevation removed by depression breaching.");
    module.def("hillshade", &hillshade, nb::arg("raster"), nb::arg("azimuth") = 315.0,
        nb::arg("altitude") = 45.0, nb::arg("z_factor") = 1.0,
        "Compute a single-illumination hillshade in the range 0 to 255.");
    module.def("multidirectional_hillshade", &multidirectional_hillshade,
        nb::arg("raster"), nb::arg("n_directions") = 8, nb::arg("altitude") = 45.0,
        nb::arg("z_factor") = 1.0, nb::arg("start_azimuth") = 0.0,
        "Compute an averaged multi-illumination hillshade in the range 0 to 255.");
    module.def("slope", &slope, nb::arg("raster"),
        "Compute 3-by-3 gradient-magnitude slope in rise/run units.");
    module.def("aspect", &aspect, nb::arg("raster"), nb::arg("lambda_metres") = 0.0,
        nb::arg("n_threads") = 0, nb::arg("flat_slope_tolerance") = 1e-12,
        "Compute polyfit aspect in degrees; lambda_metres defaults to three cells.");
    module.def("plan_curvature", &plan_curvature, nb::arg("raster"),
        "Compute 3-by-3 plan curvature; concave-up is positive.");
    module.def("total_curvature", &total_curvature, nb::arg("raster"), nb::arg("radius_metres"),
        "Compute moving-window total curvature at radius_metres.");
    module.def("make_boundary_conditions", &make_boundary_conditions,
        nb::arg("north") = "open", nb::arg("east") = "open",
        nb::arg("south") = "open", nb::arg("west") = "open", R"doc(
Create an LSDTopoTools NESW boundary-condition string.

Parameters
----------
north, east, south, west : {"open", "closed", "periodic", "base_level"}
    Boundary treatment. The returned string has NESW order.
)doc");
    module.def("make_outlet_mask", &make_outlet_mask, nb::arg("raster"),
        nb::arg("north") = false, nb::arg("east") = false,
        nb::arg("south") = false, nb::arg("west") = false,
        "Create a valid binary mask with selected raster edges set to one.");
    module.def("mask_from_points", &mask_from_points, nb::arg("raster"), nb::arg("points"), R"doc(
Create a binary mask from zero-based ``(row, col)`` point coordinates.

Parameters
----------
points : array-like, shape (n, 2)
    Row-major NumPy coordinates. Points must be valid raster cells.
)doc");
    module.def("fill_pits_with_outlets", &fill_pits_with_outlets, nb::arg("raster"),
        nb::arg("outlet_mask"), nb::arg("min_slope") = 0.0001f, nb::arg("n_threads") = 0,
        "Fill pits using positive cells in outlet_mask as drainage outlets.");
    module.def("fill_pits_with_sinks", &fill_pits_with_sinks_mask, nb::arg("raster"),
        nb::arg("sinks"), nb::arg("min_slope") = 0.0001f, nb::arg("n_threads") = 0,
        "Fill pits while retaining positive cells in sinks as declared sinks.");
    module.def("fill_pits_with_sinks", &fill_pits_with_sinks_points, nb::arg("raster"),
        nb::arg("sinks"), nb::arg("min_slope") = 0.0001f, nb::arg("n_threads") = 0,
        "Fill pits while retaining zero-based (row, col) sinks.");
    module.def("fill_pits_boundary_aware", &fill_pits_boundary_aware, nb::arg("raster"),
        nb::arg("boundary_conditions") = "oooo", nb::arg("min_slope") = 0.0001f,
        nb::arg("n_threads") = 0,
        "Fill pits with NESW open, closed, periodic, or base-level boundaries.");
    module.def("breach_least_cost", &breach_least_cost, nb::arg("raster"),
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
    module.def("least_cost_carve_depth", &least_cost_carve_depth, nb::arg("raster"),
        nb::arg("max_distance") = 1000,
        nb::arg("max_cost") = (std::numeric_limits<float>::max)(),
        nb::arg("prefer_shortest_path") = true, nb::arg("flat_increment") = 0.0f,
        nb::arg("fill_remainder") = true, nb::arg("n_threads") = 0,
        "Return the non-negative elevation removed by least-cost breaching.");
    module.def("reproject_to_grid", &reproject_to_grid, nb::arg("raster"),
        nb::arg("target_metadata"), nb::arg("resampling") = "bilinear", R"doc(
Reproject and resample a raster onto an explicit target grid.

Parameters
----------
raster : Raster
target_metadata : RasterMeta
    Target extent, resolution, dimensions, and CRS.
resampling : str, default="bilinear"
    GDAL resampling method.

Returns
-------
Raster
    New raster on target_metadata's grid.
)doc");
    module.def("reproject", &reproject_to_epsg, nb::arg("raster"), nb::arg("epsg"),
        nb::arg("cell_size") = 0.0, nb::arg("resampling") = "bilinear", R"doc(
Reproject a raster to an EPSG coordinate reference system.

Parameters
----------
epsg : int
    Target EPSG identifier.
cell_size : float, default=0
    Target cell size. A non-positive value lets GDAL choose it.
resampling : str, default="bilinear"
    GDAL resampling method.
)doc");
    module.def("mosaic", &mosaic, nb::arg("rasters"), R"doc(
Mosaic loaded rasters into their union extent.

Parameters
----------
rasters : sequence of Raster
    Rasters with matching CRS, cell size, and pixel alignment. Later valid
    pixels take precedence in overlaps.

Returns
-------
Raster
    Mosaic covering the union extent.
)doc");
    module.def("rasterize_polygons", &rasterize_polygons, nb::arg("raster"),
        nb::arg("polygons"), nb::arg("coverage") = "majority", R"doc(
Rasterize GeoPandas polygon geometries onto a raster grid.

Parameters
----------
polygons : geopandas.GeoDataFrame
    Polygon geometries with a CRS exactly matching raster.
coverage : {"any", "majority", "full"}, default="majority"
    Pixel-inclusion rule. ``"majority"`` tests the pixel centre.

Returns
-------
Raster
    Binary mask raster with valid pixels throughout the reference grid.

Raises
------
ValueError
    If the GeoDataFrame CRS does not match the raster CRS.
)doc");
    module.def("clip_to_polygons", &clip_to_polygons, nb::arg("raster"),
        nb::arg("polygons"), nb::arg("coverage") = "majority", nb::arg("crop") = true, R"doc(
Crop to polygon bounds and set pixels outside polygons to nodata.

Parameters
----------
polygons : geopandas.GeoDataFrame
    Polygon geometries with a CRS exactly matching raster.
coverage : {"any", "majority", "full"}, default="majority"
    Pixel-inclusion rule.
crop : bool, default=True
    Crop to polygon bounds before masking. False retains the source extent.
)doc");
    module.def("sample_points", &sample_points, nb::arg("raster"), nb::arg("points"),
        nb::arg("column") = "raster_value", R"doc(
Sample a raster at GeoPandas point geometries.

Parameters
----------
points : geopandas.GeoDataFrame
    Point geometries with a CRS exactly matching raster.
column : str, default="raster_value"
    Name of the output value column.

Returns
-------
geopandas.GeoDataFrame
    Copy of points with sampled values. Outside or nodata samples are None.
)doc");
}
