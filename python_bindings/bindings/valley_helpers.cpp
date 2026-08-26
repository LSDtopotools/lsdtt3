#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

#include <nanobind/nanobind.h>
#include <nanobind/stl/vector.h>

#include <cpu/lsdtt_channel_scope_cpu.hpp>
#include <cpu/lsdtt_valley_centreline_cpu.hpp>

#include "binding_common.hpp"

namespace lsdtt3_python {
namespace {

std::vector<float> fit_valley_baseline(nb::object along_distance, nb::object pin_elevation)
{
    const auto distance = nb::cast<std::vector<double>>(
        nb::module_::import_("numpy").attr("asarray")(
            along_distance, nb::arg("dtype") = "float64").attr("ravel")());
    const auto elevation = nb::cast<std::vector<float>>(
        nb::module_::import_("numpy").attr("asarray")(
            pin_elevation, nb::arg("dtype") = "float32").attr("ravel")());
    if (distance.size() != elevation.size())
        throw std::invalid_argument("along_distance and pin_elevation must have equal length");
    return lsdtt::fit_pinned_monotonic_baseline(distance, elevation);
}

lsdtt::channel_scope_result make_scope(const Raster& dem, nb::object channels)
{
    const auto& meta = dem.native().get_meta();
    const std::size_t n = dem.native().get_data().size();
    lsdtt::channel_scope_result scope;
    scope.mode = lsdtt::channel_scope_mode::user_sources;
    scope.channel_mask.assign(n, 0);
    scope.channel_source_id.assign(n, -1);
    const std::size_t n_channels = nb::len(channels);
    if (n_channels == 0) throw std::invalid_argument("channels must contain at least one path");
    scope.channel_node_lists.reserve(n_channels);
    scope.resolved_sources.reserve(n_channels);
    scope.resolved_outlets.reserve(n_channels);
    for (std::size_t source = 0; source < n_channels; ++source) {
        nb::object path = channels[nb::int_(source)];
        const auto pixels = nb::cast<std::vector<std::vector<int>>>(
            nb::module_::import_("numpy").attr("asarray")(path).attr("tolist")());
        if (pixels.size() < 2) throw std::invalid_argument("each channel path must contain at least two row, col pairs");
        std::vector<std::int64_t> nodes;
        nodes.reserve(pixels.size());
        for (const auto& pixel : pixels) {
            if (pixel.size() != 2 || pixel[0] < 0 || pixel[0] >= meta.n_rows ||
                pixel[1] < 0 || pixel[1] >= meta.n_cols)
                throw std::invalid_argument("channel pixels must be in-grid row, col pairs");
            const std::int64_t flat = static_cast<std::int64_t>(pixel[0]) * meta.n_cols + pixel[1];
            nodes.push_back(flat);
            scope.channel_mask[static_cast<std::size_t>(flat)] = 1;
            if (scope.channel_source_id[static_cast<std::size_t>(flat)] < 0)
                scope.channel_source_id[static_cast<std::size_t>(flat)] = static_cast<std::int32_t>(source);
        }
        scope.resolved_sources.push_back(nodes.front());
        scope.resolved_outlets.push_back(nodes.back());
        scope.channel_node_lists.push_back(std::move(nodes));
    }
    return scope;
}

class ValleyCentrelineResult {
public:
    explicit ValleyCentrelineResult(lsdtt::valley_centreline_result result)
        : result_(std::move(result)) {}
    const lsdtt::valley_centreline_result& native() const noexcept { return result_; }
private:
    lsdtt::valley_centreline_result result_;
};

ValleyCentrelineResult extract_valley_centrelines(
    const Raster& dem, nb::object channels, double corridor_width_metres,
    double fallback_valley_hand_metres, int valley_dilation_pixels,
    double trough_scaling_factor, int centreline_loops, int n_threads)
{
    lsdtt::valley_centreline_request request;
    request.corridor_width_metres = corridor_width_metres;
    request.fallback_valley_hand_metres = fallback_valley_hand_metres;
    request.valley_dilation_pixels = valley_dilation_pixels;
    request.trough_scaling_factor = trough_scaling_factor;
    request.centreline_loops = centreline_loops;
    request.n_threads = n_threads;
    return ValleyCentrelineResult(lsdtt::extract_valley_centrelines(
        dem.native(), make_scope(dem, channels), request));
}

nb::object centreline_table(const ValleyCentrelineResult& result)
{
    nb::list rows;
    for (const auto& channel : result.native().channels) {
        nb::dict row;
        row["source_id"] = channel.source_id;
        row["valley_hand_break_metres"] = channel.valley_hand_break_metres;
        row["n_valley_pixels"] = channel.n_valley_pixels;
        row["routing_succeeded"] = channel.routing_succeeded;
        row["trough_retries"] = channel.trough_retries;
        row["n_centreline_nodes"] = channel.centreline_flats.size();
        rows.append(std::move(row));
    }
    return nb::module_::import_("pandas").attr("DataFrame")(rows);
}

nb::object centreline_pixels(const ValleyCentrelineResult& result, std::size_t source_id)
{
    const auto& channels = result.native().channels;
    if (source_id >= channels.size()) throw std::out_of_range("source_id is out of range");
    const auto& meta = result.native().valley_mask.get_meta();
    nb::list pixels;
    for (const auto flat : channels[source_id].centreline_flats)
        pixels.append(nb::make_tuple(flat / meta.n_cols, flat % meta.n_cols));
    return nb::module_::import_("numpy").attr("array")(pixels, nb::arg("dtype") = "int64");
}

} // namespace

void bind_valley_helpers(nb::module_& m)
{
    m.def("fit_valley_baseline", &fit_valley_baseline,
        nb::arg("along_distance"), nb::arg("pin_elevation"), R"doc(
Interpolate a non-increasing valley baseline through finite elevation pins.

Use ``numpy.nan`` for unpinned centreline nodes.
)doc");
    nb::class_<ValleyCentrelineResult>(m, "ValleyCentrelineResult", R"doc(
Valley masks and centreline diagnostics for supplied channel paths.
)doc")
        .def_prop_ro("valley_mask", [](const ValleyCentrelineResult& r) { return Raster(r.native().valley_mask); })
        .def_prop_ro("table", &centreline_table)
        .def("centreline_pixels", &centreline_pixels, nb::arg("source_id"));
    m.def("extract_valley_centrelines", &extract_valley_centrelines,
        nb::arg("dem"), nb::arg("channels"),
        nb::arg("corridor_width_metres") = 5000.0,
        nb::arg("fallback_valley_hand_metres") = 10.0,
        nb::arg("valley_dilation_pixels") = 4,
        nb::arg("trough_scaling_factor") = 0.1,
        nb::arg("centreline_loops") = 5, nb::arg("n_threads") = 0, R"doc(
Identify valley floors and route centrelines from ordered channel paths.

``channels`` is a sequence of upstream-to-downstream ``(row, col)`` arrays.
)doc");
}

} // namespace lsdtt3_python
