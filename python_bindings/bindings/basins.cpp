#include <algorithm>
#include <cmath>
#include <cstdint>
#include <memory>
#include <limits>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <utility>
#include <vector>

#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>

#include <cpu/lsdtt_basin_cpu.hpp>
#include <cpu/lsdtt_junction_basin_manager_cpu.hpp>
#include <cpu/lsdtt_raster_surface_metrics_cpu.hpp>
#include <cpu/lsdtt_raster_spectral_cpu.hpp>
#include <cpu/lsdtt_side_tributary_basins_cpu.hpp>
#include <utils/lsdtt_statistics.hpp>

#include "binding_common.hpp"

namespace lsdtt3_python {
namespace {

std::vector<std::int64_t> outlet_flats(const FlowInfo& flow, nb::object outlets,
                                       std::int64_t min_upslope_pixels)
{
    if (min_upslope_pixels < 0)
        throw std::invalid_argument("min_upslope_pixels must be non-negative");
    if (outlets.is_none())
        return lsdtt::get_all_basin_outlet_flats(flow.native(), min_upslope_pixels);

    nb::object array = nb::module_::import_("numpy").attr("asarray")(outlets, nb::arg("dtype") = "int64");
    if (nb::cast<int>(array.attr("ndim")) != 2 || nb::cast<int>(array.attr("shape")[1]) != 2)
        throw std::invalid_argument("outlets must be an N by 2 row, col array");
    const auto pairs = nb::cast<std::vector<std::vector<std::int64_t>>>(array.attr("tolist")());
    std::unordered_set<std::int64_t> seen;
    std::vector<std::int64_t> result;
    for (const auto& pair : pairs) {
        const auto row = pair[0]; const auto col = pair[1];
        if (row < 0 || row >= flow.native().n_rows() || col < 0 || col >= flow.native().n_cols())
            throw std::invalid_argument("outlets contains an out-of-grid pixel");
        const auto flat = row * static_cast<std::int64_t>(flow.native().n_cols()) + col;
        if (!flow.native().is_valid_cell(flat)) throw std::invalid_argument("outlets selects a nodata pixel");
        if (seen.insert(flat).second) result.push_back(flat);
    }
    return result;
}

class BasinCollection {
public:
    BasinCollection(const BasinCollection&) = delete;
    BasinCollection& operator=(const BasinCollection&) = delete;
    BasinCollection(BasinCollection&&) = default;
    BasinCollection& operator=(BasinCollection&&) = default;
    BasinCollection(std::shared_ptr<lsdtt::lsdtt_flowinfo> flow, std::vector<std::int64_t> outlets)
        : flow_(std::move(flow))
    {
        basins_.reserve(outlets.size());
        for (const auto outlet : outlets)
            basins_.push_back(std::make_unique<lsdtt::lsdtt_basin_cpu>(outlet, *flow_));
        parent_.assign(basins_.size(), -1);
        for (std::size_t child = 0; child < basins_.size(); ++child) {
            const auto child_outlet = basins_[child]->get_basin_outlet_flat_index();
            std::int32_t best_parent = -1;
            std::int32_t best_size = std::numeric_limits<std::int32_t>::max();
            for (std::size_t candidate = 0; candidate < basins_.size(); ++candidate) {
                if (candidate == child) continue;
                const auto parent_outlet = basins_[candidate]->get_basin_outlet_flat_index();
                if (!flow_->is_upslope_of(child_outlet, parent_outlet)) continue;
                const auto n = basins_[candidate]->pixel_count();
                if (n < best_size) { best_size = n; best_parent = static_cast<std::int32_t>(candidate); }
            }
            parent_[child] = best_parent;
        }
        depth_.assign(basins_.size(), 0);
        for (std::size_t i = 0; i < basins_.size(); ++i) {
            for (auto parent = parent_[i]; parent >= 0; parent = parent_[static_cast<std::size_t>(parent)])
                ++depth_[i];
        }
    }

    std::size_t size() const noexcept { return basins_.size(); }
    const lsdtt::lsdtt_flowinfo& flow() const noexcept { return *flow_; }
    const lsdtt::lsdtt_basin_cpu& basin(std::size_t id) const {
        if (id >= basins_.size()) throw std::out_of_range("basin_id is out of range");
        return *basins_[id];
    }
    std::int32_t parent(std::size_t id) const { if (id >= parent_.size()) throw std::out_of_range("basin_id is out of range"); return parent_[id]; }
    std::int32_t depth(std::size_t id) const { if (id >= depth_.size()) throw std::out_of_range("basin_id is out of range"); return depth_[id]; }
    std::vector<std::int32_t> memberships(std::int64_t flat) const {
        std::vector<std::int32_t> ids;
        for (std::size_t i = 0; i < basins_.size(); ++i)
            if (flow_->is_upslope_of(flat, basins_[i]->get_basin_outlet_flat_index())) ids.push_back(static_cast<std::int32_t>(i));
        return ids;
    }
private:
    std::shared_ptr<lsdtt::lsdtt_flowinfo> flow_;
    std::vector<std::unique_ptr<lsdtt::lsdtt_basin_cpu>> basins_;
    std::vector<std::int32_t> parent_, depth_;
};

Raster basin_mask(const BasinCollection& collection, std::size_t id)
{
    const auto& flow = collection.flow();
    const std::size_t n = flow.full_grid_size();
    std::vector<float> data(n, 0.0f);
    std::vector<std::uint8_t> mask(flow.mask_array().data(), flow.mask_array().data() + n);
    for (const auto flat : collection.basin(id).get_slice()) data[static_cast<std::size_t>(flat)] = 1.0f;
    return Raster(lsdtt::lsdtt_raster_cpu(std::move(data), std::move(mask), flow.meta()));
}

Raster basin_labels(const BasinCollection& collection)
{
    const auto& flow = collection.flow(); const std::size_t n = flow.full_grid_size();
    std::vector<float> data(n, 0.0f);
    std::vector<std::uint8_t> mask(flow.mask_array().data(), flow.mask_array().data() + n);
    for (std::size_t id = 0; id < collection.size(); ++id)
        for (const auto flat : collection.basin(id).get_slice()) data[static_cast<std::size_t>(flat)] = static_cast<float>(id + 1);
    return Raster(lsdtt::lsdtt_raster_cpu(std::move(data), std::move(mask), flow.meta()));
}

nb::object basin_pixels(const BasinCollection& collection, std::size_t id)
{
    const auto& flow = collection.flow(); const auto& slice = collection.basin(id).get_slice();
    std::vector<std::int64_t> values; values.reserve(slice.size() * 2);
    for (const auto flat : slice) { values.push_back(flow.row_from_flat_unchecked(flat)); values.push_back(flow.col_from_flat_unchecked(flat)); }
    return nb::module_::import_("numpy").attr("array")(values, nb::arg("dtype") = "int64")
        .attr("reshape")(nb::make_tuple(static_cast<int>(slice.size()), 2));
}

nb::object basin_table(const BasinCollection& collection)
{
    const auto& flow = collection.flow(); nb::list rows;
    for (std::size_t id = 0; id < collection.size(); ++id) {
        const auto& basin = collection.basin(id); const auto outlet = basin.get_basin_outlet_flat_index();
        nb::dict row; row["basin_id"] = id; row["outlet_flat"] = outlet;
        row["outlet_row"] = flow.row_from_flat_unchecked(outlet); row["outlet_col"] = flow.col_from_flat_unchecked(outlet);
        row["n_pixels"] = basin.pixel_count(); row["area_metres2"] = basin.pixel_count() * flow.meta().cell_size * flow.meta().cell_size;
        row["complete"] = basin.check_basin_completeness(); row["outlet_on_dem_edge"] = lsdtt::is_flat_on_dem_edge(flow, outlet);
        row["parent_id"] = collection.parent(id); row["nesting_depth"] = collection.depth(id); rows.append(std::move(row));
    }
    return nb::module_::import_("pandas").attr("DataFrame")(rows);
}

nb::object basin_membership(const BasinCollection& collection, int row, int col)
{
    const auto& flow = collection.flow(); const auto flat = flow.flat_index(row, col);
    if (!flow.is_valid_cell(flat)) throw std::invalid_argument("row, col selects nodata");
    nb::object numpy = nb::module_::import_("numpy");
    return numpy.attr("array")(collection.memberships(flat), nb::arg("dtype") = "int32");
}

BasinCollection make_basins(const FlowInfo& flow, nb::object outlets, std::int64_t min_upslope_pixels)
{
    return BasinCollection(flow.shared_native(), outlet_flats(flow, outlets, min_upslope_pixels));
}

class JunctionBasinCollection {
public:
    JunctionBasinCollection(std::shared_ptr<lsdtt::lsdtt_flowinfo> flow,
                            std::shared_ptr<lsdtt::lsdtt_junction_network_cpu> network,
                            const std::string& selection, std::vector<std::int64_t> outlets)
        : flow_(std::move(flow)), network_(std::move(network)),
          manager_(std::make_unique<lsdtt::lsdtt_junction_basin_manager_cpu>(*flow_, *network_))
    {
        if (selection == "all") manager_->select_all_junction_basins();
        else if (selection == "outlets") manager_->select_junction_basins_from_basin_outlets(outlets);
        else if (selection == "nested") manager_->select_junction_basins_from_basin_outlets_nested(outlets);
        else throw std::invalid_argument("selection must be 'all', 'outlets', or 'nested'");
    }
    std::size_t size() const noexcept { return manager_->size(); }
    const lsdtt::lsdtt_junction_basin_manager_cpu& native() const noexcept { return *manager_; }
    const lsdtt::lsdtt_flowinfo& flow() const noexcept { return *flow_; }
    JunctionBasinCollection& filter_incomplete() { manager_->filter_incomplete(); return *this; }
    JunctionBasinCollection& filter_area(std::int32_t low, std::int32_t high) { manager_->filter_by_area(low, high); return *this; }
    JunctionBasinCollection& filter_stream_order(std::int32_t low, std::int32_t high) { manager_->filter_by_stream_order(low, high); return *this; }
    JunctionBasinCollection& filter_nested() { manager_->filter_nested(); return *this; }
    JunctionBasinCollection& filter_nested_by_stream_order() { manager_->filter_nested_by_stream_order(); return *this; }
    JunctionBasinCollection& filter_largest() { manager_->filter_largest(); return *this; }
private:
    std::shared_ptr<lsdtt::lsdtt_flowinfo> flow_;
    std::shared_ptr<lsdtt::lsdtt_junction_network_cpu> network_;
    std::unique_ptr<lsdtt::lsdtt_junction_basin_manager_cpu> manager_;
};

JunctionBasinCollection make_junction_basins(const ChannelNetwork& network, const std::string& selection,
                                             nb::object outlets)
{
    std::vector<std::int64_t> flats;
    if (selection != "all") {
        nb::object array = nb::module_::import_("numpy").attr("asarray")(outlets, nb::arg("dtype") = "int64");
        if (nb::cast<int>(array.attr("ndim")) != 2 || nb::cast<int>(array.attr("shape")[1]) != 2)
            throw std::invalid_argument("outlets must be an N by 2 row, col array");
        const auto values = nb::cast<std::vector<std::int64_t>>(array.attr("ravel")());
        const auto& flow = network.flow();
        for (std::size_t i = 0; i < values.size(); i += 2) {
            const int row = static_cast<int>(values[i]), col = static_cast<int>(values[i + 1]);
            if (!flow.in_bounds(row, col) || !flow.is_valid_cell(flow.flat_index(row, col)))
                throw std::invalid_argument("outlets contains an invalid flow pixel");
            flats.push_back(flow.flat_index(row, col));
        }
    }
    return JunctionBasinCollection(network.shared_flow(), network.shared_native(), selection, std::move(flats));
}

Raster junction_basin_labels(const JunctionBasinCollection& collection)
{
    auto result = collection.native().get_junction_basins_raster();
    return Raster(lsdtt::lsdtt_raster_cpu(std::move(result.data), std::move(result.mask), result.meta));
}

Raster junction_basin_mask(const JunctionBasinCollection& collection, std::size_t id)
{
    if (id >= collection.size()) throw std::out_of_range("basin_id is out of range");
    const auto& flow = collection.flow(); std::vector<float> data(flow.full_grid_size(), 0.0f);
    std::vector<std::uint8_t> mask(flow.mask_array().data(), flow.mask_array().data() + flow.full_grid_size());
    for (const auto flat : collection.native().basin_at(id).get_slice()) data[static_cast<std::size_t>(flat)] = 1.0f;
    return Raster(lsdtt::lsdtt_raster_cpu(std::move(data), std::move(mask), flow.meta()));
}

nb::object junction_basin_table(const JunctionBasinCollection& collection)
{
    const auto& flow = collection.flow(); nb::list rows;
    for (std::size_t id = 0; id < collection.size(); ++id) {
        const auto& basin = collection.native().basin_at(id); const auto flat = basin.get_basin_outlet_flat_index();
        nb::dict row; row["basin_id"] = id; row["junction_id"] = basin.get_junction_number(); row["stream_order"] = basin.get_stream_order();
        row["outlet_flat"] = flat; row["outlet_row"] = flow.row_from_flat_unchecked(flat); row["outlet_col"] = flow.col_from_flat_unchecked(flat);
        row["n_pixels"] = basin.pixel_count(); row["area_metres2"] = basin.pixel_count() * flow.cell_size() * flow.cell_size(); row["complete"] = basin.check_basin_completeness(); rows.append(std::move(row));
    }
    return nb::module_::import_("pandas").attr("DataFrame")(rows);
}

nb::object basin_statistics(const BasinCollection& collection, const Raster& dem, bool include_slope, bool include_flow_distance)
{
    const auto& flow = collection.flow();
    if (dem.native().get_meta().n_rows != flow.n_rows() || dem.native().get_meta().n_cols != flow.n_cols())
        throw std::invalid_argument("dem must have the same grid shape as the basin flow topology");
    const auto slope_raster = include_slope ? lsdtt::compute_gradient_magnitude_3x3_square(dem.native()) : lsdtt::lsdtt_raster_cpu();
    const auto distance_raster = include_flow_distance ? lsdtt::lsdtt_raster_cpu(flow.get_distance_from_outlet_raster()) : lsdtt::lsdtt_raster_cpu();
    nb::list rows;
    for (std::size_t id = 0; id < collection.size(); ++id) {
        const auto& basin = collection.basin(id); const auto z = basin.get_basin_values(dem.native()); const auto zs = lsdtt::compute_descriptive_stats(z);
        nb::dict row; row["basin_id"] = id; row["n_pixels"] = basin.pixel_count(); row["area_metres2"] = basin.pixel_count() * flow.cell_size() * flow.cell_size();
        row["elevation_min"] = zs.min; row["elevation_max"] = zs.max; row["elevation_mean"] = zs.mean; row["elevation_median"] = zs.median; row["elevation_stddev"] = zs.stddev_population;
        if (include_slope) { const auto s = basin.get_basin_values(slope_raster); row["slope_25th"] = lsdtt::compute_percentile(s, 25.0); row["slope_50th"] = lsdtt::compute_percentile(s, 50.0); row["slope_75th"] = lsdtt::compute_percentile(s, 75.0); }
        if (include_flow_distance) { const auto d = basin.get_basin_values(distance_raster); row["flow_distance_25th"] = lsdtt::compute_percentile(d, 25.0); row["flow_distance_50th"] = lsdtt::compute_percentile(d, 50.0); row["flow_distance_75th"] = lsdtt::compute_percentile(d, 75.0); }
        rows.append(std::move(row));
    }
    return nb::module_::import_("pandas").attr("DataFrame")(rows);
}

class SideTributaryBasins {
public:
    SideTributaryBasins(std::shared_ptr<lsdtt::lsdtt_flowinfo> flow, lsdtt::lsdtt_side_tributary_basin_result result)
        : flow_(std::move(flow)), result_(std::move(result))
    {
        std::vector<std::int64_t> outlets;
        for (const auto& basin : result_.accepted_basins) outlets.push_back(basin.outlet_flat);
        basins_ = std::make_unique<BasinCollection>(flow_, std::move(outlets));
    }
    const BasinCollection& basins() const noexcept { return *basins_; }
    const lsdtt::lsdtt_side_tributary_basin_result& native() const noexcept { return result_; }
    const lsdtt::lsdtt_flowinfo& flow() const noexcept { return *flow_; }
private:
    std::shared_ptr<lsdtt::lsdtt_flowinfo> flow_;
    lsdtt::lsdtt_side_tributary_basin_result result_;
    std::unique_ptr<BasinCollection> basins_;
};

SideTributaryBasins side_tributary_basins(const FlowInfo& flow, const std::string& points_path,
                                          int minimum_contributing_pixels, bool reject_incomplete,
                                          bool single_channel, bool include_stopping_point,
                                          int stopping_point_window_pixels)
{
    if (minimum_contributing_pixels < 1) throw std::invalid_argument("minimum_contributing_pixels must be positive");
    lsdtt::lsdtt_side_tributary_basin_options options;
    options.input_points_path = points_path; options.minimum_contributing_pixels = minimum_contributing_pixels;
    options.reject_incomplete_basins = reject_incomplete; options.include_stopping_point = include_stopping_point;
    options.stopping_point_window_pixels = stopping_point_window_pixels;
    auto result = single_channel ? lsdtt::extract_single_channel_side_basins(flow.native(), options)
                                 : lsdtt::extract_side_tributary_basins(flow.native(), options);
    return SideTributaryBasins(flow.shared_native(), std::move(result));
}

Raster side_channel_mask(const SideTributaryBasins& result)
{
    const auto& flow = result.flow(); std::vector<float> data(flow.full_grid_size(), 0.0f);
    for (std::size_t i = 0; i < data.size(); ++i) data[i] = result.native().traced_channel_mask[i] ? 1.0f : 0.0f;
    std::vector<std::uint8_t> mask(flow.mask_array().data(), flow.mask_array().data() + flow.full_grid_size());
    return Raster(lsdtt::lsdtt_raster_cpu(std::move(data), std::move(mask), flow.meta()));
}

nb::object side_tributary_table(const SideTributaryBasins& result)
{
    const auto& r = result.native(); std::vector<std::int64_t> flat, receiver; std::vector<int> row, col, receiver_row, receiver_col;
    std::vector<std::uint64_t> n; std::vector<double> area; std::vector<bool> complete;
    for (const auto& basin : r.accepted_basins) { flat.push_back(basin.outlet_flat); receiver.push_back(basin.receiver_channel_flat); row.push_back(basin.outlet_row); col.push_back(basin.outlet_col); receiver_row.push_back(basin.receiver_channel_row); receiver_col.push_back(basin.receiver_channel_col); n.push_back(basin.n_pixels); area.push_back(basin.drainage_area); complete.push_back(basin.complete); }
    nb::dict data; data["outlet_flat"] = flat; data["receiver_channel_flat"] = receiver; data["outlet_row"] = row; data["outlet_col"] = col; data["receiver_channel_row"] = receiver_row; data["receiver_channel_col"] = receiver_col; data["n_pixels"] = n; data["drainage_area_metres2"] = area; data["complete"] = complete;
    return nb::module_::import_("pandas").attr("DataFrame")(data);
}

BasinCollection line_crossing_basins(const ChannelNetwork& network, nb::object line, int minimum_contributing_pixels)
{
    if (minimum_contributing_pixels < 1) throw std::invalid_argument("minimum_contributing_pixels must be positive");
    nb::object geometry = line;
    if (nb::hasattr(line, "geometry")) geometry = line.attr("geometry").attr("unary_union");
    const auto& flow = network.flow(); std::vector<std::int64_t> flats;
    nb::object point = nb::module_::import_("shapely.geometry").attr("Point");
    const double tolerance = 0.71 * flow.cell_size();
    for (std::size_t i = 0; i < flow.full_grid_size(); ++i) {
        const auto flat = static_cast<std::int64_t>(i);
        if (!network.native().channel_mask_array().data()[i] || flow.s_vector_size(flat) < minimum_contributing_pixels) continue;
        const int row = flow.row_from_flat_unchecked(flat), col = flow.col_from_flat_unchecked(flat);
        const double x = flow.meta().x_min + (col + 0.5) * flow.cell_size(), y = flow.meta().y_max() - (row + 0.5) * flow.cell_size();
        if (nb::cast<double>(geometry.attr("distance")(point(x, y))) <= tolerance) flats.push_back(flat);
    }
    return BasinCollection(network.shared_flow(), std::move(flats));
}

nb::dict basin_surface_metrics(const BasinCollection& collection, const Raster& dem,
                               const std::vector<double>& scales, const std::vector<std::string>& metric_names,
                               bool include_outline, int n_threads)
{
    if (scales.empty() || metric_names.empty()) throw std::invalid_argument("scales and metrics must not be empty");
    if (dem.native().get_meta().n_rows != collection.flow().n_rows() || dem.native().get_meta().n_cols != collection.flow().n_cols())
        throw std::invalid_argument("dem must have the same grid shape as the basin flow topology");
    const std::map<std::string, lsdtt::polyfit_metric> lookup = {
        {"slope", lsdtt::polyfit_metric::slope}, {"slope_degrees", lsdtt::polyfit_metric::slope_degrees},
        {"maximum_curvature", lsdtt::polyfit_metric::maximum_curvature}, {"minimum_curvature", lsdtt::polyfit_metric::minimum_curvature},
        {"mean_curvature", lsdtt::polyfit_metric::mean_curvature}, {"gaussian_curvature", lsdtt::polyfit_metric::gaussian_curvature},
        {"profile_curvature", lsdtt::polyfit_metric::profile_curvature}, {"plan_curvature", lsdtt::polyfit_metric::plan_curvature},
        {"tangential_curvature", lsdtt::polyfit_metric::tangential_curvature}, {"aspect_degrees", lsdtt::polyfit_metric::aspect_degrees}};
    std::vector<lsdtt::polyfit_metric> metrics; for (const auto& name : metric_names) { const auto found = lookup.find(name); if (found == lookup.end()) throw std::invalid_argument("unsupported basin surface metric: " + name); metrics.push_back(found->second); }
    nb::list interior_rows, outline_rows, flip_rows;
    for (const double scale : scales) {
        if (scale <= 0.0) throw std::invalid_argument("scales must be positive");
        const auto fit = lsdtt::direct_polyfit(dem.native(), scale, n_threads); lsdtt::polyfit_metric_request request; request.metrics = metrics; request.n_threads = n_threads;
        const auto result = lsdtt::compute_polyfit_metrics(fit, request);
        for (std::size_t id = 0; id < collection.size(); ++id) {
            const auto& basin = collection.basin(id);
            for (std::size_t metric_id = 0; metric_id < metrics.size(); ++metric_id) {
                const auto& raster = result.at(metrics[metric_id]); const auto values = basin.get_basin_values(raster); const auto stats = lsdtt::compute_descriptive_stats(values);
                nb::dict row; row["basin_id"] = id; row["scale_metres"] = scale; row["metric"] = metric_names[metric_id]; row["n_values"] = stats.n_values; row["min"] = stats.min; row["mean"] = stats.mean; row["median"] = stats.median; row["max"] = stats.max; interior_rows.append(std::move(row));
                if (include_outline) {
                    const auto outline = basin.get_basin_outline(); const auto sampled = outline.interpolate_pixel_values_on_outline(raster); double distance = 0.0; const auto& rc = outline.outline_divide_rc(); const double cell = dem.native().get_meta().cell_size; int flips = 0; float previous = std::numeric_limits<float>::quiet_NaN();
                    for (std::size_t node = 0; node < sampled.size(); ++node) { if (node > 0) distance += cell * std::hypot(static_cast<double>(rc[node].row - rc[node - 1].row), static_cast<double>(rc[node].col - rc[node - 1].col)); nb::dict sample; sample["basin_id"] = id; sample["scale_metres"] = scale; sample["metric"] = metric_names[metric_id]; sample["distance_metres"] = distance; sample["value"] = sampled[node]; outline_rows.append(std::move(sample)); if (metric_names[metric_id] == "maximum_curvature" && std::isfinite(previous) && std::isfinite(sampled[node]) && previous * sampled[node] < 0.0f) ++flips; if (std::isfinite(sampled[node])) previous = sampled[node]; }
                    if (metric_names[metric_id] == "maximum_curvature") { nb::dict f; f["basin_id"] = id; f["scale_metres"] = scale; f["n_flips"] = flips; f["mean_flip_spacing_metres"] = flips > 0 ? distance / flips : std::numeric_limits<double>::quiet_NaN(); flip_rows.append(std::move(f)); }
                }
            }
        }
    }
    nb::dict output; const auto pandas = nb::module_::import_("pandas"); output["interior"] = pandas.attr("DataFrame")(interior_rows); output["outline"] = pandas.attr("DataFrame")(outline_rows); output["curvature_flips"] = pandas.attr("DataFrame")(flip_rows); return output;
}

} // namespace

void bind_basins(nb::module_& module, nb::class_<FlowInfo>& flow)
{
    nb::class_<BasinCollection>(module, "BasinCollection", R"doc(
Selected nested D8 drainage basins retaining their parent ``FlowInfo`` topology.

Each basin is defined by an outlet pixel and stores an efficient native
upslope slice; no raster-sized copy is retained per basin.
)doc")
        .def_prop_ro("size", &BasinCollection::size)
        .def_prop_ro("labels", &basin_labels, "Basin IDs (one-based); nested basins overwrite parents.")
        .def_prop_ro("table", &basin_table)
        .def("mask", &basin_mask, nb::arg("basin_id"))
        .def("pixels", &basin_pixels, nb::arg("basin_id"))
        .def("membership", &basin_membership, nb::arg("row"), nb::arg("col"), "Return all selected basin IDs containing this pixel.")
        .def("statistics", &basin_statistics, nb::arg("dem"), nb::arg("include_slope") = true,
             nb::arg("include_flow_distance") = true,
             "Return per-basin elevation, slope, and flow-distance statistics.");
    nb::class_<JunctionBasinCollection>(module, "JunctionBasinCollection", R"doc(
Junction-aware D8 basins selected from a ``ChannelNetwork``.

The collection retains the channel-network stream order and supports in-place
selection filters. It is created with :func:`junction_basins`.
)doc")
        .def_prop_ro("size", &JunctionBasinCollection::size)
        .def_prop_ro("table", &junction_basin_table)
        .def_prop_ro("labels", &junction_basin_labels)
        .def("mask", &junction_basin_mask, nb::arg("basin_id"))
        .def("filter_incomplete", &JunctionBasinCollection::filter_incomplete, nb::rv_policy::reference_internal)
        .def("filter_area", &JunctionBasinCollection::filter_area, nb::arg("min_pixels") = -1, nb::arg("max_pixels") = -1, nb::rv_policy::reference_internal)
        .def("filter_stream_order", &JunctionBasinCollection::filter_stream_order, nb::arg("min_order") = -1, nb::arg("max_order") = -1, nb::rv_policy::reference_internal)
        .def("filter_nested", &JunctionBasinCollection::filter_nested, nb::rv_policy::reference_internal)
        .def("filter_nested_by_stream_order", &JunctionBasinCollection::filter_nested_by_stream_order, nb::rv_policy::reference_internal)
        .def("filter_largest", &JunctionBasinCollection::filter_largest, nb::rv_policy::reference_internal);
    module.def("junction_basins", &make_junction_basins, nb::arg("network"), nb::arg("selection") = "all", nb::arg("outlets") = nb::none(), nb::keep_alive<0, 1>(), R"doc(
Select junction-aware basins from a channel network.

``selection`` is ``"all"``, ``"outlets"``, or ``"nested"``. The latter two
require outlet pixels as an N-by-2 ``(row, col)`` array.
)doc");
    nb::class_<SideTributaryBasins>(module, "SideTributaryBasins", "Side basins draining into a traced D8 channel.")
        .def_prop_ro("table", &side_tributary_table)
        .def_prop_ro("channel_mask", &side_channel_mask)
        .def_prop_ro("basins", &SideTributaryBasins::basins, nb::rv_policy::reference_internal);
    module.def("side_tributary_basins", &side_tributary_basins, nb::arg("flow"), nb::arg("points_path"),
        nb::arg("minimum_contributing_pixels") = 1000, nb::arg("reject_incomplete") = true,
        nb::arg("single_channel") = false, nb::arg("include_stopping_point") = false,
        nb::arg("stopping_point_window_pixels") = 5);
    module.def("line_crossing_basins", &line_crossing_basins, nb::arg("network"), nb::arg("line"),
        nb::arg("minimum_contributing_pixels") = 500);
    module.def("basin_surface_metrics", &basin_surface_metrics, nb::arg("basins"), nb::arg("dem"),
        nb::arg("scales"), nb::arg("metrics"), nb::arg("include_outline") = true, nb::arg("n_threads") = 0);
    flow.def("basins", &make_basins, nb::arg("outlets") = nb::none(), nb::arg("min_upslope_pixels") = 0,
        R"doc(
Create a ``BasinCollection`` from explicit outlet pixels or all qualifying D8 outlets.

Explicit outlets may be internal channel pixels, which creates nested basins.
)doc");
}
} // namespace lsdtt3_python
