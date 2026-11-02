#include <cmath>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>

#include <cpu/lsdtt_channel_scope_cpu.hpp>
#include <cpu/lsdtt_valley_centreline_cpu.hpp>
#include <cpu/lsdtt_valley_cross_sections_cpu.hpp>
#include <cpu/lsdtt_valley_relative_elevation_cpu.hpp>
#include <cpu/lsdtt_valley_terrace_profile_cpu.hpp>
#include <cpu/lsdtt_valley_width_network_cpu.hpp>

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

// The centrelines, and what the relative elevation and terrace profile need
// beyond them: each channel's pixels (head -> downstream), the flow whose
// drainage area ranks the channels, and the supplied valley floor.
class ValleyCentrelineResult {
public:
    explicit ValleyCentrelineResult(lsdtt::valley_centreline_result result)
        : result_(std::move(result)) {}
    ValleyCentrelineResult(lsdtt::valley_centreline_result result,
                           std::vector<std::vector<std::int64_t>> channel_flats,
                           std::shared_ptr<lsdtt::lsdtt_flowinfo> flow,
                           lsdtt::lsdtt_raster_cpu floor_mask)
        : result_(std::move(result)), channel_flats_(std::move(channel_flats)),
          flow_(std::move(flow)), floor_mask_(std::move(floor_mask)) {}
    const lsdtt::valley_centreline_result& native() const noexcept { return result_; }
    const std::vector<std::vector<std::int64_t>>& channel_flats() const noexcept { return channel_flats_; }
    const lsdtt::lsdtt_flowinfo* flow() const noexcept { return flow_.get(); }
    const lsdtt::lsdtt_raster_cpu& floor_mask() const noexcept { return floor_mask_; }
private:
    lsdtt::valley_centreline_result result_;
    std::vector<std::vector<std::int64_t>> channel_flats_;
    std::shared_ptr<lsdtt::lsdtt_flowinfo> flow_;
    lsdtt::lsdtt_raster_cpu floor_mask_;
};

class ValleyWidthNetwork {
public:
    ValleyWidthNetwork(lsdtt::valley_width_network_result result,
                       std::shared_ptr<lsdtt::lsdtt_flowinfo> flow,
                       lsdtt::lsdtt_raster_cpu floor_mask)
        : result_(std::move(result)), flow_(std::move(flow)), floor_mask_(std::move(floor_mask)) {}
    const lsdtt::valley_width_network_result& native() const noexcept { return result_; }
    ValleyCentrelineResult centrelines() const
    {
        std::vector<std::vector<std::int64_t>> lists;
        for (const auto& seg : result_.segments) lists.push_back(seg.channel_flats);
        return ValleyCentrelineResult(result_.centrelines, std::move(lists), flow_, floor_mask_);
    }
private:
    lsdtt::valley_width_network_result result_;
    std::shared_ptr<lsdtt::lsdtt_flowinfo> flow_;
    lsdtt::lsdtt_raster_cpu floor_mask_;
};

// Re-pin each centreline's baseline as lsdtt-valley-centreline does with
// centreline_baseline_algorithm_choice (the engine pins projected_monotonic).
void repin_baselines(lsdtt::valley_centreline_result& result,
                     const std::vector<std::vector<std::int64_t>>& channel_flats,
                     const Raster& dem, const std::string& choice_name)
{
    const lsdtt::centreline_baseline_choice choice =
        lsdtt::parse_centreline_baseline_choice(choice_name);
    if (choice == lsdtt::centreline_baseline_choice::projected_monotonic) return;
    for (std::size_t s = 0; s < result.channels.size() && s < channel_flats.size(); ++s)
        result.channels[s].centreline_baseline_elevation = lsdtt::compute_centreline_baseline(
            result.channels[s].centreline_flats, channel_flats[s], dem.native().get_data(),
            dem.native().get_meta(), choice);
}

class ValleyCrossSections {
public:
    explicit ValleyCrossSections(lsdtt::valley_cross_section_result result)
        : result_(std::move(result)) {}
    const lsdtt::valley_cross_section_result& native() const noexcept { return result_; }
private:
    lsdtt::valley_cross_section_result result_;
};

lsdtt::valley_centreline_request make_centreline_request(
    const Raster* valley_mask, double corridor_width_metres,
    double fallback_valley_hand_metres, int valley_dilation_pixels,
    double trough_scaling_factor, int centreline_loops, int n_threads)
{
    lsdtt::valley_centreline_request request;
    request.supplied_valley_mask =
        valley_mask == nullptr ? nullptr : &valley_mask->native();
    request.corridor_width_metres = corridor_width_metres;
    request.fallback_valley_hand_metres = fallback_valley_hand_metres;
    request.valley_dilation_pixels = valley_dilation_pixels;
    request.trough_scaling_factor = trough_scaling_factor;
    request.centreline_loops = centreline_loops;
    request.n_threads = n_threads;
    return request;
}

void require_valley_mask_grid(const Raster& dem, const Raster& valley_mask)
{
    const auto& a = dem.native().get_meta();
    const auto& b = valley_mask.native().get_meta();
    if (a.n_rows != b.n_rows || a.n_cols != b.n_cols || a.x_min != b.x_min ||
        a.y_min != b.y_min || a.cell_size != b.cell_size || a.crs_wkt != b.crs_wkt)
        throw std::invalid_argument("valley_mask must share the DEM grid");
}

ValleyCentrelineResult extract_valley_centrelines(
    const Raster& dem, nb::object channels, double corridor_width_metres,
    nb::object fallback_valley_hand_metres, int valley_dilation_pixels,
    double trough_scaling_factor, int centreline_loops, int n_threads,
    const Raster& valley_mask, const std::string& centreline_baseline_algorithm_choice)
{
    if (!fallback_valley_hand_metres.is_none() &&
        PyErr_WarnEx(PyExc_DeprecationWarning,
            "'fallback_valley_hand_metres' is deprecated and ignored when a valley mask is supplied", 1) < 0)
        throw nb::python_error();
    require_valley_mask_grid(dem, valley_mask);
    const auto request = make_centreline_request(
        &valley_mask, corridor_width_metres, 10.0,
        valley_dilation_pixels, trough_scaling_factor, centreline_loops,
        n_threads);
    const lsdtt::channel_scope_result scope = make_scope(dem, channels);
    auto result = lsdtt::extract_valley_centrelines(dem.native(), scope, request);
    repin_baselines(result, scope.channel_node_lists, dem, centreline_baseline_algorithm_choice);
    return ValleyCentrelineResult(std::move(result), scope.channel_node_lists, nullptr,
                                  valley_mask.native());
}

ValleyCentrelineResult extract_valley_centrelines_from_scope(
    const Raster& dem, const ChannelScope& scope, double corridor_width_metres,
    nb::object fallback_valley_hand_metres, int valley_dilation_pixels,
    double trough_scaling_factor, int centreline_loops, int n_threads,
    const Raster& valley_mask, const std::string& centreline_baseline_algorithm_choice)
{
    if (!fallback_valley_hand_metres.is_none() &&
        PyErr_WarnEx(PyExc_DeprecationWarning,
            "'fallback_valley_hand_metres' is deprecated and ignored when a valley mask is supplied", 1) < 0)
        throw nb::python_error();
    require_valley_mask_grid(dem, valley_mask);
    const auto& a = dem.native().get_meta();
    const auto& b = scope.flow().meta();
    if (a.n_rows != b.n_rows || a.n_cols != b.n_cols || a.x_min != b.x_min ||
        a.y_min != b.y_min || a.cell_size != b.cell_size || a.crs_wkt != b.crs_wkt)
        throw std::invalid_argument("dem must share the ChannelScope flow grid");
    const auto request = make_centreline_request(
        &valley_mask, corridor_width_metres, 10.0,
        valley_dilation_pixels, trough_scaling_factor, centreline_loops,
        n_threads);
    auto result = lsdtt::extract_valley_centrelines(dem.native(), scope.native(), request);
    repin_baselines(result, scope.native().channel_node_lists, dem,
                    centreline_baseline_algorithm_choice);
    return ValleyCentrelineResult(std::move(result), scope.native().channel_node_lists,
                                  scope.shared_flow(), valley_mask.native());
}

ValleyWidthNetwork build_valley_width_network(
    const Raster& dem, const ChannelScope& scope, const Raster& valley_mask,
    const std::string& network_centreline_algorithm_choice,
    double network_centreline_window_margin_multiple,
    double network_centreline_smoothing_width_multiple, int min_centreline_nodes,
    const std::string& centreline_baseline_algorithm_choice)
{
    require_valley_mask_grid(dem, valley_mask);
    const auto& a = dem.native().get_meta();
    const auto& b = scope.flow().meta();
    if (a.n_rows != b.n_rows || a.n_cols != b.n_cols || a.x_min != b.x_min ||
        a.y_min != b.y_min || a.cell_size != b.cell_size || a.crs_wkt != b.crs_wkt)
        throw std::invalid_argument("dem must share the ChannelScope flow grid");
    lsdtt::valley_width_network_request request;
    request.centreline = lsdtt::parse_network_centreline_choice(network_centreline_algorithm_choice);
    request.window_margin_multiple = network_centreline_window_margin_multiple;
    request.smoothing_width_multiple = network_centreline_smoothing_width_multiple;
    request.min_segment_pixels = min_centreline_nodes;
    request.baseline = lsdtt::parse_centreline_baseline_choice(centreline_baseline_algorithm_choice);
    std::vector<std::uint8_t> channel_mask(scope.native().channel_mask.begin(),
                                           scope.native().channel_mask.end());
    return ValleyWidthNetwork(lsdtt::build_valley_width_network(
        dem.native(), scope.flow(), channel_mask, valley_mask.native(), request),
        scope.shared_flow(), valley_mask.native());
}

nb::object network_segment_table(const ValleyWidthNetwork& network)
{
    nb::list rows;
    for (const auto& seg : network.native().segments) {
        nb::dict row;
        row["segment_id"] = seg.segment_id;
        row["parent_segment_id"] = seg.parent_segment_id;
        row["nesting_depth"] = seg.nesting_depth;
        row["n_channel_pixels"] = seg.channel_flats.size();
        row["drains_to_edge"] = seg.drains_to_edge;
        row["catchment_touches_edge"] = seg.catchment_touches_edge;
        row["valley_width_estimate_metres"] = seg.valley_width_estimate_metres;
        row["medial_path_failed"] = seg.medial_path_failed;
        row["n_samples_snapped_to_mask_edge"] = seg.n_samples_snapped_to_mask_edge;
        rows.append(std::move(row));
    }
    return nb::module_::import_("pandas").attr("DataFrame")(rows);
}

// The centrelines the relative elevation measures against, built by the
// routine lsdtt-valley-centreline uses.
std::vector<lsdtt::valley_rer_centreline> rer_centrelines(const ValleyCentrelineResult& c)
{
    if (c.flow() == nullptr)
        throw std::invalid_argument(
            "the relative elevation ranks channels by drainage area, so it needs centrelines "
            "from extract_valley_centrelines(dem, scope, ...) or build_valley_width_network");
    return lsdtt::make_rer_centrelines(c.native(), c.channel_flats(),
                                       c.flow()->get_drainage_area_raster().data,
                                       c.native().valley_mask.get_meta());
}

lsdtt::lsdtt_raster_cpu whole_dem_domain(const Raster& dem)
{
    return lsdtt::lsdtt_raster_cpu(std::vector<float>(dem.native().get_data().size(), 1.0f),
                                   std::vector<std::uint8_t>(dem.native().get_mask()),
                                   dem.native().get_meta());
}

nb::tuple relative_elevation(const ValleyCentrelineResult& centrelines, const Raster& dem,
                             const std::string& relative_elevation_domain_algorithm_choice,
                             nb::object domain_mask,
                             double relative_elevation_max_distance_metres)
{
    require_valley_mask_grid(dem, Raster(centrelines.floor_mask()));
    const auto rc = rer_centrelines(centrelines);
    const lsdtt::lsdtt_raster_cpu& floor = centrelines.floor_mask();
    lsdtt::lsdtt_raster_cpu domain = floor;
    if (relative_elevation_domain_algorithm_choice == "supplied") {
        if (domain_mask.is_none())
            throw std::invalid_argument(
                "relative_elevation_domain_algorithm_choice supplied needs domain_mask");
        const Raster& d = nb::cast<const Raster&>(domain_mask);
        require_valley_mask_grid(dem, d);
        domain = d.native();
    } else if (relative_elevation_domain_algorithm_choice == "whole_dem") {
        domain = whole_dem_domain(dem);
    } else if (relative_elevation_domain_algorithm_choice != "valley_mask") {
        throw std::invalid_argument("unknown relative_elevation_domain_algorithm_choice '" +
            relative_elevation_domain_algorithm_choice + "'. Valid: valley_mask, supplied, whole_dem.");
    }
    lsdtt::valley_rer_request request;
    request.max_distance_metres = relative_elevation_max_distance_metres;
    const lsdtt::valley_rer_result rr = lsdtt::compute_valley_relative_elevation(
        dem.native(), floor, domain, rc, request);
    // Outside the domain and the floor is nodata, as the CLI writes it.
    const auto& meta = dem.native().get_meta();
    auto clipped = [&](const lsdtt::lsdtt_raster_cpu& r) {
        std::vector<float> d = r.get_data();
        std::vector<std::uint8_t> m = r.get_mask();
        for (std::size_t i = 0; i < d.size(); ++i)
            if (!((domain.get_mask()[i] && domain.get_data()[i] > 0.5f) ||
                  (floor.get_mask()[i] && floor.get_data()[i] > 0.5f))) {
                d[i] = static_cast<float>(meta.no_data_value);
                m[i] = 0;
            }
        return Raster(lsdtt::lsdtt_raster_cpu(std::move(d), std::move(m), r.get_meta()));
    };
    return nb::make_tuple(clipped(rr.relative_elevation), clipped(rr.owner));
}

nb::tuple terrace_profile(const ValleyCentrelineResult& centrelines, const Raster& dem,
                          const Raster& surface_patch_labels,
                          double relative_elevation_max_distance_metres,
                          double terrace_profile_bin_metres)
{
    require_valley_mask_grid(dem, surface_patch_labels);
    const auto rc = rer_centrelines(centrelines);
    // Every patch pixel must reach a node, so ownership is computed over the
    // whole DEM, as lsdtt-valley-centreline does for the profile.
    lsdtt::valley_rer_request request;
    request.max_distance_metres = relative_elevation_max_distance_metres;
    const lsdtt::valley_rer_result rr = lsdtt::compute_valley_relative_elevation(
        dem.native(), centrelines.floor_mask(), whole_dem_domain(dem), rc, request);
    lsdtt::terrace_profile_request treq;
    treq.bin_metres = terrace_profile_bin_metres;
    const lsdtt::terrace_profile_result tp = lsdtt::compute_terrace_profile(
        dem.native(), surface_patch_labels.native(), rc, rr, treq);
    nb::list bins, patches;
    const char* stat_names[5] = {"min", "p25", "median", "p75", "max"};
    for (const auto& b : tp.bins) {
        nb::dict row;
        row["centreline_id"] = b.centreline_id;
        row["patch_label"] = b.patch_label;
        row["bin_centre_metres"] = b.bin_centre_metres;
        row["n_pixels"] = b.n_pixels;
        for (int k = 0; k < 5; ++k) {
            row[(std::string("elevation_") + stat_names[k]).c_str()] = b.elevation[k];
            row[(std::string("relative_elevation_") + stat_names[k]).c_str()] = b.relative_elevation[k];
        }
        bins.append(std::move(row));
    }
    for (const auto& p : tp.patches) {
        nb::dict row;
        row["patch_label"] = p.patch_label;
        row["centreline_id"] = p.centreline_id;
        row["n_pixels"] = p.n_pixels;
        row["n_pixels_reassigned"] = p.n_pixels_reassigned;
        row["along_min_metres"] = p.along_min_metres;
        row["along_max_metres"] = p.along_max_metres;
        row["elevation_median"] = p.elevation_median;
        row["relative_elevation_median"] = p.relative_elevation_median;
        patches.append(std::move(row));
    }
    nb::object df = nb::module_::import_("pandas").attr("DataFrame");
    return nb::make_tuple(df(bins), df(patches));
}

#ifdef LSDTT_PYTHON_COMPATIBILITY_PROBE
ValleyCentrelineResult extract_valley_centrelines_with_mask_probe(
    const Raster& dem, nb::object channels, const Raster& valley_mask,
    double corridor_width_metres, double fallback_valley_hand_metres,
    int valley_dilation_pixels, double trough_scaling_factor,
    int centreline_loops, int n_threads)
{
    const auto request = make_centreline_request(
        &valley_mask, corridor_width_metres, fallback_valley_hand_metres,
        valley_dilation_pixels, trough_scaling_factor, centreline_loops,
        n_threads);
    return ValleyCentrelineResult(lsdtt::extract_valley_centrelines(
        dem.native(), make_scope(dem, channels), request));
}

ValleyCentrelineResult extract_valley_centrelines_from_scope_with_mask_probe(
    const Raster& dem, const ChannelScope& scope, const Raster& valley_mask,
    double corridor_width_metres, double fallback_valley_hand_metres,
    int valley_dilation_pixels, double trough_scaling_factor,
    int centreline_loops, int n_threads)
{
    const auto& a = dem.native().get_meta();
    const auto& b = scope.flow().meta();
    if (a.n_rows != b.n_rows || a.n_cols != b.n_cols || a.x_min != b.x_min ||
        a.y_min != b.y_min || a.cell_size != b.cell_size || a.crs_wkt != b.crs_wkt)
        throw std::invalid_argument("dem must share the ChannelScope flow grid");
    const auto request = make_centreline_request(
        &valley_mask, corridor_width_metres, fallback_valley_hand_metres,
        valley_dilation_pixels, trough_scaling_factor, centreline_loops,
        n_threads);
    return ValleyCentrelineResult(lsdtt::extract_valley_centrelines(
        dem.native(), scope.native(), request));
}
#endif

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

ValleyCrossSections extract_cross_sections(const Raster& dem, const ValleyCentrelineResult& centrelines,
                                           double spacing_metres, double half_length_multiplier,
                                           double min_half_length_metres, double max_half_length_metres,
                                           int tangent_window_nodes, double sample_step_metres,
                                           bool curved_sections, nb::object curvature,
                                           nb::object hillshade, nb::object channel_mask,
                                           int side_gap_tolerance_steps,
                                           double feature_peak_min_prominence_fraction,
                                           int feature_peak_edge_margin_samples,
                                           bool clip_overlapping_sections,
                                           int neighbour_clip_window,
                                           double min_section_reach_fraction,
                                           double curved_wall_clearance_fraction,
                                           int curved_streamwise_smoothing_passes,
                                           double curved_coarse_spacing_metres,
                                           double hydraulic_stage_step_metres,
                                           double hydraulic_min_step_height_metres,
                                           int hydraulic_min_step_width_pixels,
                                           int hydraulic_midpoint_search_samples,
                                           int hydraulic_min_crossings_preferred,
                                           bool compute_hydraulic_depth)
{
    if (spacing_metres <= 0.0 || half_length_multiplier <= 0.0)
        throw std::invalid_argument("spacing_metres and half_length_multiplier must be positive");
    lsdtt::valley_cross_section_request request;
    request.spacing_metres = spacing_metres; request.half_length_multiplier = half_length_multiplier;
    request.min_half_length_metres = min_half_length_metres; request.max_half_length_metres = max_half_length_metres;
    request.tangent_window_nodes = tangent_window_nodes; request.sample_step_metres = sample_step_metres;
    request.curved_sections = curved_sections;
    request.side_gap_tolerance_steps = side_gap_tolerance_steps;
    request.feature_peak_min_prominence_fraction = feature_peak_min_prominence_fraction;
    request.feature_peak_edge_margin_samples = feature_peak_edge_margin_samples;
    request.clip_overlapping_sections = clip_overlapping_sections;
    request.neighbour_clip_window = neighbour_clip_window;
    request.min_section_reach_fraction = min_section_reach_fraction;
    request.curved_wall_clearance_fraction = curved_wall_clearance_fraction;
    request.curved_streamwise_smoothing_passes = curved_streamwise_smoothing_passes;
    request.curved_coarse_spacing_metres = curved_coarse_spacing_metres;
    request.hydraulic_stage_step_metres = hydraulic_stage_step_metres;
    request.hydraulic_min_step_height_metres = hydraulic_min_step_height_metres;
    request.hydraulic_min_step_width_pixels = hydraulic_min_step_width_pixels;
    request.hydraulic_midpoint_search_samples = hydraulic_midpoint_search_samples;
    request.hydraulic_min_crossings_preferred = hydraulic_min_crossings_preferred;
    request.compute_hydraulic_depth = compute_hydraulic_depth;
    const auto optional_raster = [](nb::object value) -> const lsdtt::lsdtt_raster_cpu* {
        return value.is_none() ? nullptr : &nb::cast<const Raster&>(value).native();
    };
    return ValleyCrossSections(lsdtt::extract_valley_cross_sections(
        dem.native(), centrelines.native(), request, optional_raster(curvature),
        optional_raster(hillshade), optional_raster(channel_mask)));
}

nb::object cross_section_stations(const ValleyCrossSections& sections)
{
    std::vector<std::int32_t> source, station; std::vector<double> along, width, cx, cy, left, right;
    std::vector<bool> left_found, right_found;
    for (const auto& s : sections.native().sections) { source.push_back(s.source_id); station.push_back(s.station_seq); along.push_back(s.along_distance_metres); width.push_back(s.valley_width_metres); cx.push_back(s.centre_x); cy.push_back(s.centre_y); left.push_back(s.left_offset_metres); right.push_back(s.right_offset_metres); left_found.push_back(s.left_found); right_found.push_back(s.right_found); }
    nb::dict data; data["source_id"] = source; data["station_id"] = station; data["along_distance_metres"] = along; data["valley_width_metres"] = width; data["centre_x"] = cx; data["centre_y"] = cy; data["left_offset_metres"] = left; data["right_offset_metres"] = right; data["left_found"] = left_found; data["right_found"] = right_found;
    return nb::module_::import_("pandas").attr("DataFrame")(data);
}

nb::object cross_section_samples(const ValleyCrossSections& sections)
{
    std::vector<std::int32_t> source, station; std::vector<double> offset, x, y; std::vector<float> elevation, curvature, hillshade;
    std::vector<bool> valid, in_valley, curvature_valid, hillshade_valid, is_channel;
    for (const auto& section : sections.native().sections) for (const auto& sample : section.samples) {
        source.push_back(section.source_id); station.push_back(section.station_seq); offset.push_back(sample.offset_metres);
        x.push_back(sample.x); y.push_back(sample.y); elevation.push_back(sample.elevation); valid.push_back(sample.elevation_valid);
        in_valley.push_back(sample.in_valley); curvature.push_back(sample.curvature); curvature_valid.push_back(sample.curvature_valid);
        hillshade.push_back(sample.hillshade); hillshade_valid.push_back(sample.hillshade_valid); is_channel.push_back(sample.is_channel);
    }
    nb::dict data; data["source_id"] = source; data["station_id"] = station; data["offset_metres"] = offset; data["x"] = x; data["y"] = y; data["elevation"] = elevation; data["elevation_valid"] = valid; data["in_valley"] = in_valley; data["curvature"] = curvature; data["curvature_valid"] = curvature_valid; data["hillshade"] = hillshade; data["hillshade_valid"] = hillshade_valid; data["is_channel"] = is_channel;
    return nb::module_::import_("pandas").attr("DataFrame")(data);
}

nb::object cross_section_features(const ValleyCrossSections& sections)
{
    std::vector<std::int32_t> source, station; std::vector<std::string> kind; std::vector<double> offset;
    for (const auto& section : sections.native().sections) {
        const auto append = [&](const char* feature, const std::vector<double>& offsets) {
            for (const double value : offsets) {
                source.push_back(section.source_id); station.push_back(section.station_seq);
                kind.emplace_back(feature); offset.push_back(value);
            }
        };
        append("channel_crossing", section.channel_crossing_offsets);
        append("curvature_peak", section.curvature_peak_offsets);
        append("hillshade_peak", section.hillshade_peak_offsets);
    }
    nb::dict data; data["source_id"] = source; data["station_id"] = station;
    data["feature"] = kind; data["offset_metres"] = offset;
    return nb::module_::import_("pandas").attr("DataFrame")(data);
}

nb::object cross_section_hydraulic_edges(const ValleyCrossSections& sections)
{
    std::vector<std::int32_t> source, station, step; std::vector<double> stage, offset;
    for (const auto& section : sections.native().sections) for (const auto& edge : section.hydraulic_edges) {
        source.push_back(section.source_id); station.push_back(section.station_seq); step.push_back(edge.step);
        stage.push_back(edge.stage_metres); offset.push_back(edge.offset_metres);
    }
    nb::dict data; data["source_id"] = source; data["station_id"] = station; data["step"] = step;
    data["stage_metres"] = stage; data["offset_metres"] = offset;
    return nb::module_::import_("pandas").attr("DataFrame")(data);
}

nb::object cross_section_hydraulic_breaks(const ValleyCrossSections& sections)
{
    std::vector<std::int32_t> source, station, break_id; std::vector<double> stage, depth;
    for (const auto& section : sections.native().sections) {
        if (section.hydraulic_break_stages.empty()) continue;
        const double thalweg = section.hydraulic_break_stages.front();
        for (std::size_t i = 0; i < section.hydraulic_break_stages.size(); ++i) {
            source.push_back(section.source_id); station.push_back(section.station_seq); break_id.push_back(static_cast<std::int32_t>(i));
            stage.push_back(section.hydraulic_break_stages[i]); depth.push_back(section.hydraulic_break_stages[i] - thalweg);
        }
    }
    nb::dict data; data["source_id"] = source; data["station_id"] = station; data["break_id"] = break_id;
    data["stage_metres"] = stage; data["depth_above_thalweg_metres"] = depth;
    return nb::module_::import_("pandas").attr("DataFrame")(data);
}

nb::object cross_section_hydraulic_curve(const ValleyCrossSections& sections)
{
    std::vector<std::int32_t> source, station; std::vector<double> stage, depth, width, area, perimeter, hydraulic_depth;
    for (const auto& section : sections.native().sections) for (const auto& sample : section.hydraulic_curve) {
        source.push_back(section.source_id); station.push_back(section.station_seq); stage.push_back(sample.stage_metres);
        depth.push_back(sample.depth_above_thalweg_metres); width.push_back(sample.width_metres); area.push_back(sample.area_metres2);
        perimeter.push_back(sample.wetted_perimeter_metres); hydraulic_depth.push_back(sample.hydraulic_depth_metres);
    }
    nb::dict data; data["source_id"] = source; data["station_id"] = station; data["stage_metres"] = stage;
    data["depth_above_thalweg_metres"] = depth; data["width_metres"] = width; data["area_metres2"] = area;
    data["wetted_perimeter_metres"] = perimeter; data["hydraulic_depth_metres"] = hydraulic_depth;
    return nb::module_::import_("pandas").attr("DataFrame")(data);
}

Raster cross_section_valley_mask(const Raster& dem, const ValleyCentrelineResult& centrelines,
                                 const ValleyCrossSections& sections)
{
    return Raster(lsdtt::cross_section_informed_valley_mask(
        dem.native(), centrelines.native(), sections.native()).valley_mask);
}

nb::object cross_section_levels(const ValleyCrossSections& sections,
                                const std::string& ladder_break_algorithm_choice,
                                int ladder_max_level_count, double vertical_accuracy_metres,
                                double ladder_resolvability_kappa,
                                double cross_section_level_min_section_fraction,
                                bool include_dropped)
{
    // The ladder request lsdtt-valley-centreline builds for
    // write_cross_section_levels_csv.
    lsdtt::ladder_request request;
    if (ladder_break_algorithm_choice == "gap") request.break_algorithm = lsdtt::ladder_break_choice::gap;
    else if (ladder_break_algorithm_choice == "antimode") request.break_algorithm = lsdtt::ladder_break_choice::antimode;
    else if (ladder_break_algorithm_choice == "jenks") request.break_algorithm = lsdtt::ladder_break_choice::jenks;
    else throw std::invalid_argument("unknown ladder_break_algorithm_choice '" +
        ladder_break_algorithm_choice + "'. Valid: gap, antimode, jenks.");
    request.max_level_count = ladder_max_level_count;
    request.vertical_accuracy_metres = vertical_accuracy_metres;
    request.resolvability_kappa = ladder_resolvability_kappa;
    const lsdtt::cross_section_levels_result result = lsdtt::cluster_cross_section_levels(
        sections.native(), request, cross_section_level_min_section_fraction);
    static const char* status_names[] = {"resolved", "unresolved", "forced", "degenerate", "supplied"};
    nb::list rows;
    auto append = [&](const std::vector<lsdtt::cross_section_level>& levels, bool kept) {
        for (const auto& lv : levels) {
            nb::dict row;
            row["source_id"] = lv.source_id;
            row["window_centre_metres"] = lv.window_centre_m;
            row["level_index"] = lv.level_index;
            row["level_hand_metres"] = lv.level_hand_m;
            row["n_breaks"] = lv.n_breaks;
            row["n_sections"] = lv.n_sections;
            row["spread_metres"] = lv.spread_m;
            row["status"] = status_names[static_cast<int>(lv.status)];
            row["required_gap_metres"] = lv.required_gap_m;
            row["min_level_gap_metres"] = lv.min_level_gap_m;
            row["kept"] = kept;
            rows.append(std::move(row));
        }
    };
    append(result.levels, true);
    if (include_dropped) append(result.dropped, false);
    return nb::module_::import_("pandas").attr("DataFrame")(rows);
}

nb::object network_section_exclusions(const ValleyWidthNetwork& network,
                                      const ValleyCrossSections& sections,
                                      bool run_tributary_stop_at_parent,
                                      double tributary_stop_width_fraction,
                                      double tributary_parent_width_window_multiple,
                                      bool run_tributary_crossing_removal)
{
    lsdtt::network_section_rules_request request;
    request.stop_tributaries_at_parent = run_tributary_stop_at_parent;
    request.tributary_stop_width_fraction = tributary_stop_width_fraction;
    request.parent_width_window_multiple = tributary_parent_width_window_multiple;
    request.remove_tributary_crossing_sections = run_tributary_crossing_removal;
    const auto excluded = lsdtt::apply_network_section_rules(
        network.native(), sections.native(), network.native().centrelines.valley_mask.get_meta(),
        request);
    std::vector<std::int32_t> source, station;
    std::vector<std::string> reason;
    for (std::size_t i = 0; i < excluded.size(); ++i) {
        source.push_back(sections.native().sections[i].source_id);
        station.push_back(sections.native().sections[i].station_seq);
        reason.emplace_back(lsdtt::to_string(excluded[i]));
    }
    nb::dict data; data["source_id"] = source; data["station_id"] = station; data["exclusion"] = reason;
    return nb::module_::import_("pandas").attr("DataFrame")(data);
}

} // namespace

void bind_valley_helpers(nb::module_& m)
{
    m.def("fit_valley_baseline", &fit_valley_baseline,
        nb::arg("along_distance"), nb::arg("pin_elevation"), R"doc(
Interpolate a non-increasing valley baseline through finite elevation pins.

Use ``numpy.nan`` for unpinned centreline nodes.

Parameters
----------
along_distance : array-like of float
    Distance (m) of each centreline node along the centreline.
pin_elevation : array-like of float
    Elevation (m) pinned at each node, NaN where unpinned.

Returns
-------
numpy.ndarray
    Baseline elevation (m) at every node, never increasing downstream.
)doc");
    nb::class_<ValleyCentrelineResult>(m, "ValleyCentrelineResult", R"doc(
Valley masks and centreline diagnostics for supplied channel paths.

Made by :func:`extract_valley_centrelines`; pass it to
:func:`extract_valley_cross_sections`.
)doc")
        .def_prop_ro("valley_mask", [](const ValleyCentrelineResult& r) { return Raster(r.native().valley_mask); },
            "0/1 valley-floor Raster used for routing.")
        .def_prop_ro("table", &centreline_table, nb::sig("def table(self) -> pandas.DataFrame"),
            "One row per channel: source_id, valley_hand_break_metres, n_valley_pixels, routing_succeeded, trough_retries, n_centreline_nodes.")
        .def("centreline_pixels", &centreline_pixels, nb::arg("source_id"), R"doc(
Return the valley centreline of one channel, upstream to downstream.

Parameters
----------
source_id : int
    Zero-based channel index (``table.source_id``).

Returns
-------
numpy.ndarray
    N x 2 int64 array of ``(row, col)``.
)doc")
        .def("relative_elevation", &relative_elevation, nb::arg("dem"),
            nb::arg("relative_elevation_domain_algorithm_choice") = "valley_mask",
            nb::arg("domain_mask") = nb::none(),
            nb::arg("relative_elevation_max_distance_metres") = 0.0, R"doc(
Relative elevation above the valley centrelines' baselines (REM).

lsdtt-valley-centreline ``write_relative_elevation_raster``. On the valley
floor each cell is measured against the centreline node of the largest
channel whose floor it lies on; elsewhere against the node reached by the
shortest floor-width-scaled ground distance. Needs centrelines from a
``ChannelScope`` or :func:`build_valley_width_network` (the channels are
ranked by drainage area).

Parameters
----------
dem : Raster
    Elevation (m), the DEM the centrelines were made on.
relative_elevation_domain_algorithm_choice : {"valley_mask", "supplied", "whole_dem"}, default="valley_mask"
    Where the field is computed: the valley floor, ``domain_mask``, or
    every valid cell. The floor is always included.
domain_mask : Raster or None, default=None
    Cells > 0.5, for ``supplied``.
relative_elevation_max_distance_metres : float, default=0.0
    Largest ground distance (m) from a cell to its node; 0 means none.

Returns
-------
tuple of Raster
    ``(relative_elevation, owner)``: elevation above the owning node's
    baseline (m) and the owning centreline's ``source_id``; nodata outside
    the domain and the floor.
)doc")
        .def("terrace_profile", &terrace_profile, nb::arg("dem"),
            nb::arg("surface_patch_labels"),
            nb::arg("relative_elevation_max_distance_metres") = 0.0,
            nb::arg("terrace_profile_bin_metres") = 100.0, R"doc(
Surface patches along the valley centrelines: the terrace long profile.

lsdtt-valley-centreline ``write_terrace_profile_csv``. Every patch pixel is
tied to a centreline node (ownership computed over the whole DEM) and
binned by distance along that centreline.

Parameters
----------
dem : Raster
    Elevation (m).
surface_patch_labels : Raster
    Patch label per cell (0 = none), for example
    ``FlatSurfaceClassification.banded_labels``.
relative_elevation_max_distance_metres : float, default=0.0
    Largest ground distance (m) from a cell to its node; 0 means none.
terrace_profile_bin_metres : float, default=100.0
    Bin width (m) along the centreline.

Returns
-------
tuple of pandas.DataFrame
    ``(bins, patches)``. ``bins``: centreline_id, patch_label,
    bin_centre_metres, n_pixels, and elevation_* and
    relative_elevation_* five-number summaries (min, p25, median, p75,
    max). ``patches``: patch_label, centreline_id, n_pixels,
    n_pixels_reassigned, along_min_metres, along_max_metres,
    elevation_median, relative_elevation_median.
)doc");
    nb::class_<ValleyCrossSections>(m, "ValleyCrossSections", R"doc(
Valley cross sections made by :func:`extract_valley_cross_sections`.

All tables are pandas DataFrames keyed by ``source_id`` (channel) and
``station_id`` (section); distances are metres.
)doc")
        .def_prop_ro("n_sections", [](const ValleyCrossSections& s) { return s.native().sections.size(); },
            "Number of cross sections.")
        .def_prop_ro("stations", &cross_section_stations, nb::sig("def stations(self) -> pandas.DataFrame"),
            "One row per section: source_id, station_id, along_distance_metres, valley_width_metres, centre_x, centre_y, left_offset_metres, right_offset_metres, left_found, right_found.")
        .def_prop_ro("samples", &cross_section_samples, nb::sig("def samples(self) -> pandas.DataFrame"),
            "One row per transect sample: source_id, station_id, offset_metres, x, y, elevation, elevation_valid, in_valley, curvature, curvature_valid, hillshade, hillshade_valid, is_channel.")
        .def_prop_ro("features", &cross_section_features,
nb::sig("def features(self) -> pandas.DataFrame"),
            "Return detected channel crossings, curvature peaks, and hillshade peaks.")
        .def_prop_ro("hydraulic_edges", &cross_section_hydraulic_edges,
nb::sig("def hydraulic_edges(self) -> pandas.DataFrame"),
            "Return hydraulic step-boundary picks for every cross section.")
        .def_prop_ro("hydraulic_breaks", &cross_section_hydraulic_breaks,
nb::sig("def hydraulic_breaks(self) -> pandas.DataFrame"),
            "Return hydraulic break stages and local heights above the thalweg.")
        .def_prop_ro("hydraulic_curve", &cross_section_hydraulic_curve,
nb::sig("def hydraulic_curve(self) -> pandas.DataFrame"),
            "Return optional BankfullMapper-style hydraulic geometry curves.")
        .def("levels", &cross_section_levels,
            nb::arg("ladder_break_algorithm_choice") = "antimode",
            nb::arg("ladder_max_level_count") = 12, nb::arg("vertical_accuracy_metres") = 0.0,
            nb::arg("ladder_resolvability_kappa") = 3.0,
            nb::arg("cross_section_level_min_section_fraction") = 0.5,
            nb::arg("include_dropped") = false, R"doc(
Surface levels per channel from every section's hydraulic breaks.

lsdtt-valley-centreline ``write_cross_section_levels_csv``: the breaks of
all of a channel's sections are pooled and split by the ladder's rule. The
result can be passed, grouped by ``source_id``, as ``supplied_levels`` to
``classify_flat_surfaces_with_scope``. The breaks are those in
``hydraulic_breaks``.

Parameters
----------
ladder_break_algorithm_choice : {"antimode", "jenks", "gap"}, default="antimode"
    How the pooled breaks are split.
ladder_max_level_count : int, default=12
    Most levels tried.
vertical_accuracy_metres : float, default=0.0
    DEM vertical accuracy E_v (m).
ladder_resolvability_kappa : float, default=3.0
    Levels must be at least kappa x E_v apart.
cross_section_level_min_section_fraction : float, default=0.5
    A level found on fewer than this fraction of the channel's sections is
    dropped.
include_dropped : bool, default=False
    Also return the dropped levels (``kept`` False).

Returns
-------
pandas.DataFrame
    One row per level: source_id, window_centre_metres, level_index,
    level_hand_metres, n_breaks, n_sections, spread_metres, status,
    required_gap_metres, min_level_gap_metres, kept.
)doc");
#ifdef LSDTT_PYTHON_COMPATIBILITY_PROBE
    m.def("_extract_valley_centrelines_with_mask_probe",
        &extract_valley_centrelines_from_scope_with_mask_probe,
        nb::arg("dem"), nb::arg("scope"), nb::arg("valley_mask"),
        nb::arg("corridor_width_metres") = 5000.0,
        nb::arg("fallback_valley_hand_metres") = 10.0,
        nb::arg("valley_dilation_pixels") = 4,
        nb::arg("trough_scaling_factor") = 0.1,
        nb::arg("centreline_loops") = 5, nb::arg("n_threads") = 0);
    m.def("_extract_valley_centrelines_with_mask_probe",
        &extract_valley_centrelines_with_mask_probe,
        nb::arg("dem"), nb::arg("channels"), nb::arg("valley_mask"),
        nb::arg("corridor_width_metres") = 5000.0,
        nb::arg("fallback_valley_hand_metres") = 10.0,
        nb::arg("valley_dilation_pixels") = 4,
        nb::arg("trough_scaling_factor") = 0.1,
        nb::arg("centreline_loops") = 5, nb::arg("n_threads") = 0);
#endif
    m.def("extract_valley_centrelines", &extract_valley_centrelines_from_scope,
        nb::arg("dem"), nb::arg("scope"),
        nb::arg("corridor_width_metres") = 5000.0,
        nb::arg("fallback_valley_hand_metres") = nb::none(),
        nb::arg("valley_dilation_pixels") = 4,
        nb::arg("trough_scaling_factor") = 0.1,
        nb::arg("centreline_loops") = 5, nb::arg("n_threads") = 0,
        nb::kw_only(), nb::arg("valley_mask"),
        nb::arg("centreline_baseline_algorithm_choice") = "projected_monotonic", R"doc(
Route centrelines through a supplied valley mask using a resolved ``ChannelScope``.
)doc");
    m.def("extract_valley_centrelines", &extract_valley_centrelines,
        nb::arg("dem"), nb::arg("channels"),
        nb::arg("corridor_width_metres") = 5000.0,
        nb::arg("fallback_valley_hand_metres") = nb::none(),
        nb::arg("valley_dilation_pixels") = 4,
        nb::arg("trough_scaling_factor") = 0.1,
        nb::arg("centreline_loops") = 5, nb::arg("n_threads") = 0,
        nb::kw_only(), nb::arg("valley_mask"),
        nb::arg("centreline_baseline_algorithm_choice") = "projected_monotonic", R"doc(
Route centrelines through a supplied valley mask from ordered channel paths.

``channels`` is a sequence of upstream-to-downstream ``(row, col)`` arrays.
A synthetic trough is carved into the valley floor, deepest along the
valley axis, and flow is routed down it, so the centreline follows the
middle of the valley rather than the channel. This function has two call
forms: pass a ``ChannelScope`` as ``scope`` (its traced channels are used) or
the ordered channel paths as ``channels``.

Parameters
----------
dem : Raster
    Elevation raster (m).
scope : ChannelScope
    First call form: channels resolved by ``FlowInfo.channel_scope``.
channels : sequence of numpy.ndarray
    Second call form: N x 2 ``(row, col)`` arrays, upstream to downstream.
corridor_width_metres : float, default=5000.0
    Half-width (m) of the corridor around each channel searched for valley
    pixels.
fallback_valley_hand_metres : float or None, default=None
    Deprecated and ignored when ``valley_mask`` is supplied (a
    ``DeprecationWarning`` is emitted if set).
valley_dilation_pixels : int, default=4
    Pixels by which valley patches are dilated (towards the channel only) to
    fill holes.
trough_scaling_factor : float, default=0.1
    Gradient (m/m) of the synthetic trough carved towards the valley axis:
    the trough is lowered by this times the distance to the valley wall
    (Clubb et al. 2022 method); doubled on retries.
centreline_loops : int, default=5
    Number of trough-carving and routing iterations.
n_threads : int, default=0
    OpenMP thread count; 0 uses the runtime default.
valley_mask : Raster, keyword-only
    0/1 valley-floor raster on the DEM grid, for example from
    ``find_valleys`` or ``valley_mask_from_cross_sections``.
centreline_baseline_algorithm_choice : {"projected_monotonic", "channel_crossings"}, default="projected_monotonic"
    How each centreline's baseline elevation is pinned to its channel, as
    in lsdtt-valley-centreline: channel elevations projected onto the
    centreline and made non-increasing, or pinned only where the centreline
    crosses the channel.

Returns
-------
ValleyCentrelineResult
)doc");

    nb::class_<ValleyWidthNetwork>(m, "ValleyWidthNetwork", R"doc(
A valley centreline for every channel segment of a network, from
:func:`build_valley_width_network`.

``centrelines`` feeds :func:`extract_valley_cross_sections` as usual;
``segments`` describes the tree (each tributary stops at its parent's
channel); ``section_exclusions`` applies the network's rules to the cross
sections.
)doc")
        .def_prop_ro("centrelines", &ValleyWidthNetwork::centrelines,
            "ValleyCentrelineResult with one channel per segment (source_id = segment_id).")
        .def_prop_ro("segments", &network_segment_table, nb::sig("def segments(self) -> pandas.DataFrame"),
            "One row per segment: segment_id, parent_segment_id (-1 for a trunk), nesting_depth, n_channel_pixels, drains_to_edge, catchment_touches_edge, valley_width_estimate_metres, medial_path_failed, n_samples_snapped_to_mask_edge.")
        .def("section_exclusions", &network_section_exclusions, nb::arg("sections"),
            nb::arg("run_tributary_stop_at_parent") = true,
            nb::arg("tributary_stop_width_fraction") = 1.0,
            nb::arg("tributary_parent_width_window_multiple") = 2.0,
            nb::arg("run_tributary_crossing_removal") = true, R"doc(
Which cross sections the network rules exclude, and why.

As lsdtt-valley-centreline with ``run_valley_width_network``.

Parameters
----------
sections : ValleyCrossSections
    Sections made from ``centrelines``.
run_tributary_stop_at_parent : bool, default=True
    Exclude a tributary's sections within its parent's valley.
tributary_stop_width_fraction : float, default=1.0
    Multiple of the parent's half-width within which a tributary section
    is excluded.
tributary_parent_width_window_multiple : float, default=2.0
    Window (multiples of the parent's width) the parent's width is taken
    over.
run_tributary_crossing_removal : bool, default=True
    Exclude sections that cross a tributary's channel.

Returns
-------
pandas.DataFrame
    One row per section: source_id, station_id, exclusion (``none``,
    ``not_in_valley``, ``side_not_found``, ``near_parent_channel`` or
    ``tributary_crossing``).
)doc");

    m.def("build_valley_width_network", &build_valley_width_network,
        nb::arg("dem"), nb::arg("scope"), nb::arg("valley_mask"),
        nb::arg("network_centreline_algorithm_choice") = "medial_path",
        nb::arg("network_centreline_window_margin_multiple") = 3.0,
        nb::arg("network_centreline_smoothing_width_multiple") = 1.0,
        nb::arg("min_centreline_nodes") = 10,
        nb::arg("centreline_baseline_algorithm_choice") = "projected_monotonic", R"doc(
Valley centrelines for every channel of a network, for valley widths.

lsdtt-valley-centreline ``run_valley_width_network``: the channel network
is split into segments at confluences (each tributary's trace stops at its
parent), and each segment gets a centreline along the valley's medial axis
or a smoothed channel.

Parameters
----------
dem : Raster
    The routing DEM; its elevations pin the baselines.
scope : ChannelScope
    The channels (every channel of the scope is used).
valley_mask : Raster
    Valley floor (cells > 0.5) on the DEM grid.
network_centreline_algorithm_choice : {"medial_path", "smoothed_channel"}, default="medial_path"
    Medial axis of the valley mask, or the channel smoothed over a valley
    width.
network_centreline_window_margin_multiple : float, default=3.0
    ``medial_path``: walk window padding, as a multiple of the largest
    distance-to-wall at the segment's channel pixels.
network_centreline_smoothing_width_multiple : float, default=1.0
    ``smoothed_channel``: Gaussian sigma as a multiple of the valley width.
min_centreline_nodes : int, default=10
    Segments with fewer channel pixels get no centreline.
centreline_baseline_algorithm_choice : {"projected_monotonic", "channel_crossings"}, default="projected_monotonic"
    How each baseline is pinned to its channel.

Returns
-------
ValleyWidthNetwork
)doc");
    m.def("extract_valley_cross_sections", &extract_cross_sections,
        nb::arg("dem"), nb::arg("centrelines"), nb::arg("spacing_metres") = 500.0,
        nb::arg("half_length_multiplier") = 1.5, nb::arg("min_half_length_metres") = 0.0,
        nb::arg("max_half_length_metres") = 5000.0, nb::arg("tangent_window_nodes") = 8,
        nb::arg("sample_step_metres") = 0.0, nb::arg("curved_sections") = false,
        nb::arg("curvature") = nb::none(), nb::arg("hillshade") = nb::none(),
        nb::arg("channel_mask") = nb::none(), nb::arg("side_gap_tolerance_steps") = 1,
        nb::arg("feature_peak_min_prominence_fraction") = 0.15,
        nb::arg("feature_peak_edge_margin_samples") = 2,
        nb::arg("clip_overlapping_sections") = true, nb::arg("neighbour_clip_window") = 2,
        nb::arg("min_section_reach_fraction") = 0.35,
        nb::arg("curved_wall_clearance_fraction") = 0.10,
        nb::arg("curved_streamwise_smoothing_passes") = 40,
        nb::arg("curved_coarse_spacing_metres") = 700.0,
        nb::arg("hydraulic_stage_step_metres") = 0.1,
        nb::arg("hydraulic_min_step_height_metres") = 2.0,
        nb::arg("hydraulic_min_step_width_pixels") = 3,
        nb::arg("hydraulic_midpoint_search_samples") = 81,
        nb::arg("hydraulic_min_crossings_preferred") = 2,
        nb::arg("compute_hydraulic_depth") = false, R"doc(
Extract straight or boundary-fitted valley cross sections.

Optional curvature, hillshade, and channel-mask rasters enable feature picks.
Hydraulic break picks are always available; set ``compute_hydraulic_depth`` to
retain the full per-stage hydraulic-geometry curve.

Stations are placed along each centreline and a transect is cast across the
valley at each; all transects of a channel have the same length (a multiple
of its widest valley width).

Parameters
----------
dem : Raster
    Elevation raster (m).
centrelines : ValleyCentrelineResult
    From :func:`extract_valley_centrelines`.
spacing_metres : float, default=500.0
    Spacing (m) of stations along the centreline (positive).
half_length_multiplier : float, default=1.5
    Transect total length as a multiple of the channel's widest bounded
    valley width (positive).
min_half_length_metres : float, default=0.0
    Floor (m) on the transect half-length; 0 means none.
max_half_length_metres : float, default=5000.0
    Cap (m) on the outward search for a valley side and on the half-length.
tangent_window_nodes : int, default=8
    Centreline nodes on each side used to smooth the local direction.
sample_step_metres : float, default=0.0
    Sample spacing (m) along a transect; 0 uses the cell size.
curved_sections : bool, default=False
    Trace boundary-fitted sections that curve through bends instead of
    straight transects.
curvature, hillshade : Raster or None, default=None
    Optional rasters sampled along transects for feature picks.
channel_mask : Raster or None, default=None
    Optional 0/1 channel raster for channel-crossing picks.
side_gap_tolerance_steps : int, default=1
    Out-of-valley steps tolerated before a valley side is declared.
feature_peak_min_prominence_fraction : float, default=0.15
    Minimum peak prominence, as a fraction of the profile range.
feature_peak_edge_margin_samples : int, default=2
    Ignore peaks this close (samples) to the transect ends.
clip_overlapping_sections : bool, default=True
    Clip transects where they would cross neighbouring sections.
neighbour_clip_window : int, default=2
    Number of neighbouring stations on each side considered for clipping.
min_section_reach_fraction : float, default=0.35
    Drop a station whose clipped length is below this fraction of the full
    transect length.
curved_wall_clearance_fraction : float, default=0.10
    Curved sections only: a clipped side still extends at least this
    fraction beyond the valley wall.
curved_streamwise_smoothing_passes : int, default=40
    Curved sections only: smoothing passes of the along-valley coordinate.
curved_coarse_spacing_metres : float, default=700.0
    Curved sections only: spacing (m) of the first, coarse set of stations
    before midpoints are added.
hydraulic_stage_step_metres : float, default=0.1
    Flood-stage increment (m) used to build width against stage.
hydraulic_min_step_height_metres : float, default=2.0
    Minimum height (m) between accepted hydraulic step boundaries.
hydraulic_min_step_width_pixels : int, default=3
    Minimum width change (pixels) between accepted step boundaries.
hydraulic_midpoint_search_samples : int, default=81
    Search resolution for placing each step midpoint.
hydraulic_min_crossings_preferred : int, default=2
    Preferred number of profile crossings at a step midpoint (2 = one bank
    pair).
compute_hydraulic_depth : bool, default=False
    Store the full stage / width / area / hydraulic-depth curve of each
    section (``hydraulic_curve``; Delchiaro et al. 2025 BankfullMapper).

Returns
-------
ValleyCrossSections
)doc");
    m.def("valley_mask_from_cross_sections", &cross_section_valley_mask,
        nb::arg("dem"), nb::arg("centrelines"), nb::arg("sections"), R"doc(
Build a valley-floor mask from the median lowest hydraulic break per channel.

This complements the corridor/HAND mask returned by ``extract_valley_centrelines``.

Parameters
----------
dem : Raster
    Elevation raster (m).
centrelines : ValleyCentrelineResult
    Centrelines the sections were cut from.
sections : ValleyCrossSections
    From :func:`extract_valley_cross_sections`.

Returns
-------
Raster
    0/1 valley-floor mask.
)doc");
}

} // namespace lsdtt3_python
