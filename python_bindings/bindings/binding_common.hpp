#pragma once
#include <cmath>
#include <string>
#include <vector>
#include <utility>
#include <memory>
#include <optional>
#include <nanobind/nanobind.h>
#include <core/lsdtt_data_loader.hpp>
#include <core/lsdtt_raster_fetch.hpp>
#include <core/lsdtt_raster_meta.hpp>
#include <cpu/lsdtt_raster_cpu.hpp>
#include <cpu/lsdtt_raster_resampling_cpu.hpp>
#include <cpu/lsdtt_flowinfo_cpu.hpp>
#include <cpu/lsdtt_channel_scope_cpu.hpp>
#include <cpu/lsdtt_junction_network_cpu.hpp>
#include <cpu/lsdtt_flowgraph_cpu.hpp>
#include <cpu/lsdtt_hydraulic_graph_cpu.hpp>
#include <cpu/lsdtt_chi_extraction_cpu.hpp>
#include <cpu/lsdtt_chi_network_analysis_cpu.hpp>
#include <cpu/lsdtt_chi_ksn_inversion_cpu.hpp>
#include <utils/lsdtt_histogram.hpp>
namespace nb = nanobind;
namespace lsdtt3_python {
inline lsdtt::raster_resampling_method raster_resampling_method(
    const std::string& method)
{
    if (method == "nearest") return lsdtt::raster_resampling_method::nearest;
    if (method == "linear" || method == "bilinear") return lsdtt::raster_resampling_method::linear;
    if (method == "mean" || method == "average") return lsdtt::raster_resampling_method::mean;
    if (method == "sum") return lsdtt::raster_resampling_method::sum;
    if (method == "min" || method == "minimum") return lsdtt::raster_resampling_method::minimum;
    if (method == "max" || method == "maximum") return lsdtt::raster_resampling_method::maximum;
    if (method == "cubic") return lsdtt::raster_resampling_method::cubic;
    if (method == "mode") return lsdtt::raster_resampling_method::mode;
    if (method == "conservative_mean") return lsdtt::raster_resampling_method::conservative_mean;
    throw std::invalid_argument(
        "method must be nearest, linear, mean, sum, min, max, cubic, mode, or conservative_mean");
}
// Pixel containing map point (x, y) in the grid's own CRS, using the
// RasterMeta convention: row 0 is the north edge, a pixel holds
// [x_min + col*cs, x_min + (col+1)*cs) x (y_max - (row+1)*cs, y_max - row*cs].
// Returns false (row = col = -1) for a point off the grid or a non-finite one.
inline bool pixel_of_xy(const lsdtt::lsdtt_raster_meta& meta, double x, double y, int& row, int& col)
{
    row = col = -1;
    if (!std::isfinite(x) || !std::isfinite(y)) return false;
    const double fcol = std::floor((x - meta.x_min) / meta.cell_size);
    const double frow = std::floor((meta.y_max() - y) / meta.cell_size);
    if (fcol < 0.0 || frow < 0.0 || fcol >= meta.n_cols || frow >= meta.n_rows) return false;
    row = static_cast<int>(frow); col = static_cast<int>(fcol);
    return true;
}
// Transform points in place from source_crs (any GDAL user CRS string, such
// as "EPSG:4326" or WKT; lon/lat axis order for geographic CRSs) into the
// CRS of meta. Throws std::invalid_argument when meta carries no CRS or
// source_crs is not understood (raster_types.cpp).
void transform_points_to_raster_crs(std::vector<double>& x, std::vector<double>& y,
                                    const std::string& source_crs,
                                    const lsdtt::lsdtt_raster_meta& meta);
class Raster { public: explicit Raster(lsdtt::lsdtt_raster_cpu raster): raster_(std::move(raster)) {} static Raster read(const std::string& path) { lsdtt::lsdtt_data_loader loader; loader.open_raster(path); return Raster(lsdtt::lsdtt_raster_cpu(loader)); } void write(const std::string& path, const std::string& format) const { lsdtt::lsdtt_data_loader loader; loader.write_full(raster_.get_data(), raster_.get_meta(), path, lsdtt::resolve_raster_output_format(format)); } const lsdtt::lsdtt_raster_cpu& native() const noexcept { return raster_; } private: lsdtt::lsdtt_raster_cpu raster_; };
class FlowInfo { public: explicit FlowInfo(lsdtt::lsdtt_flowinfo flow): flow_(std::make_shared<lsdtt::lsdtt_flowinfo>(std::move(flow))) {} FlowInfo(const Raster& dem, const std::string& boundary_conditions, bool interior_nodata_is_outlet, nb::object forced_outlets, bool allow_pits_as_outlets, int n_threads); const lsdtt::lsdtt_flowinfo& native() const noexcept { return *flow_; } std::shared_ptr<lsdtt::lsdtt_flowinfo> shared_native() const noexcept { return flow_; } private: std::shared_ptr<lsdtt::lsdtt_flowinfo> flow_; };
class ChannelHeads { public: ChannelHeads(std::vector<std::int64_t> flats, std::string method, std::shared_ptr<lsdtt::lsdtt_flowinfo> flow, std::optional<lsdtt::lsdtt_raster_cpu> filtered_dem = std::nullopt, std::optional<lsdtt::lsdtt_raster_cpu> tangential_curvature = std::nullopt, std::optional<lsdtt::lsdtt_raster_cpu> channel_mask = std::nullopt, std::optional<lsdtt::lsdtt_raster_cpu> skeleton = std::nullopt, lsdtt::lsdtt_qq_data qq_data = {}): flats_(std::move(flats)), method_(std::move(method)), flow_(std::move(flow)), filtered_dem_(std::move(filtered_dem)), tangential_curvature_(std::move(tangential_curvature)), channel_mask_(std::move(channel_mask)), skeleton_(std::move(skeleton)), qq_data_(std::move(qq_data)) {} const std::vector<std::int64_t>& flats() const noexcept { return flats_; } const std::string& method() const noexcept { return method_; } const lsdtt::lsdtt_flowinfo& flow() const noexcept { return *flow_; } std::shared_ptr<lsdtt::lsdtt_flowinfo> shared_flow() const noexcept { return flow_; } const std::optional<lsdtt::lsdtt_raster_cpu>& filtered_dem() const noexcept { return filtered_dem_; } const std::optional<lsdtt::lsdtt_raster_cpu>& tangential_curvature() const noexcept { return tangential_curvature_; } const std::optional<lsdtt::lsdtt_raster_cpu>& channel_mask() const noexcept { return channel_mask_; } const std::optional<lsdtt::lsdtt_raster_cpu>& skeleton() const noexcept { return skeleton_; } const lsdtt::lsdtt_qq_data& qq_data() const noexcept { return qq_data_; } private: std::vector<std::int64_t> flats_; std::string method_; std::shared_ptr<lsdtt::lsdtt_flowinfo> flow_; std::optional<lsdtt::lsdtt_raster_cpu> filtered_dem_, tangential_curvature_, channel_mask_, skeleton_; lsdtt::lsdtt_qq_data qq_data_; };
class ChannelScope { public: ChannelScope(lsdtt::channel_scope_result result, std::shared_ptr<lsdtt::lsdtt_flowinfo> flow): result_(std::move(result)), flow_(std::move(flow)) {} const lsdtt::channel_scope_result& native() const noexcept { return result_; } const lsdtt::lsdtt_flowinfo& flow() const noexcept { return *flow_; } std::shared_ptr<lsdtt::lsdtt_flowinfo> shared_flow() const noexcept { return flow_; } private: lsdtt::channel_scope_result result_; std::shared_ptr<lsdtt::lsdtt_flowinfo> flow_; };
class ChannelNetwork { public: ChannelNetwork(std::shared_ptr<lsdtt::lsdtt_flowinfo> flow, std::vector<std::int64_t> sources): flow_(std::move(flow)), network_(std::make_shared<lsdtt::lsdtt_junction_network_cpu>(*flow_, sources)) {} const lsdtt::lsdtt_junction_network_cpu& native() const noexcept { return *network_; } const lsdtt::lsdtt_flowinfo& flow() const noexcept { return *flow_; } std::shared_ptr<lsdtt::lsdtt_flowinfo> shared_flow() const noexcept { return flow_; } std::shared_ptr<lsdtt::lsdtt_junction_network_cpu> shared_native() const noexcept { return network_; } private: std::shared_ptr<lsdtt::lsdtt_flowinfo> flow_; std::shared_ptr<lsdtt::lsdtt_junction_network_cpu> network_; };
class ChiCollection { public: ChiCollection(lsdtt::lsdtt_chi_collection_result result, std::shared_ptr<lsdtt::lsdtt_flowinfo> flow): result_(std::move(result)), flow_(std::move(flow)) {} const lsdtt::lsdtt_chi_collection_result& native() const noexcept { return result_; } lsdtt::lsdtt_chi_collection_result& mutable_native() noexcept { return result_; } const lsdtt::lsdtt_flowinfo& flow() const noexcept { return *flow_; } std::shared_ptr<lsdtt::lsdtt_flowinfo> shared_flow() const noexcept { return flow_; } private: lsdtt::lsdtt_chi_collection_result result_; std::shared_ptr<lsdtt::lsdtt_flowinfo> flow_; };
class ChiFit { public: ChiFit(lsdtt::lsdtt_chi_network_fit_collection_result result, std::shared_ptr<lsdtt::lsdtt_flowinfo> flow): result_(std::move(result)), flow_(std::move(flow)) {} const lsdtt::lsdtt_chi_network_fit_collection_result& native() const noexcept { return result_; } const lsdtt::lsdtt_flowinfo& flow() const noexcept { return *flow_; } private: lsdtt::lsdtt_chi_network_fit_collection_result result_; std::shared_ptr<lsdtt::lsdtt_flowinfo> flow_; };
class ChiInversion { public: ChiInversion(lsdtt::lsdtt_chi_ksn_inversion_collection_result result): result_(std::move(result)) {} const lsdtt::lsdtt_chi_ksn_inversion_collection_result& native() const noexcept { return result_; } private: lsdtt::lsdtt_chi_ksn_inversion_collection_result result_; };
class FlowGraph { public: FlowGraph(const Raster& dem, const std::string& method, const std::string& boundary_conditions, bool interior_nodata_is_outlet, float mfd_exponent, int n_threads); const lsdtt::lsdtt_flowgraph_cpu& native() const noexcept { return graph_; } private: lsdtt::lsdtt_flowgraph_cpu graph_; };
void bind_raster_types(nb::module_&); void bind_raster_basic_ops(nb::module_&); void bind_raster_conditioning(nb::module_&); void bind_raster_resampling(nb::module_&); void bind_raster_hierarchy(nb::module_&); void bind_raster_pit_processing(nb::module_&); void bind_raster_surface_metrics(nb::module_&); void bind_raster_spectral(nb::module_&); void bind_raster_swath(nb::module_&); void bind_raster_mountain_front(nb::module_&); void bind_raster_feature_finder(nb::module_&); void bind_raster_flat_surfaces(nb::module_&); void bind_fluvial_surfaces(nb::module_&); void bind_raster_marsh(nb::module_&); void bind_raster_ridgeline(nb::module_&); void bind_valley_helpers(nb::module_&); void bind_flow_info(nb::module_&); void bind_channel_heads(nb::module_&, nb::class_<FlowInfo>&); void bind_flow_graph(nb::module_&); void bind_hydraulic_routing(nb::module_&); void bind_flow_indices(nb::module_&); void bind_basins(nb::module_&, nb::class_<FlowInfo>&); void bind_divide_migration(nb::module_&); void bind_raster_geospatial(nb::module_&); void bind_geospatial_runtime(nb::module_&); void bind_raster_fetch(nb::module_&); void bind_rotation_frames(nb::module_&);
void bind_cosmo_column(nb::module_&);
void bind_cosmo_raster(nb::module_&);
void bind_cosmo_basin(nb::module_&);
void bind_statistical_comparison(nb::module_&);
void bind_vector_hydrology(nb::module_&);
void bind_internal_drainage(nb::module_&);
void bind_hillslope_transport(nb::module_&);
void bind_landscape_model(nb::module_&);
void bind_sediment_routing(nb::module_&);
void bind_landslides(nb::module_&);
void bind_synthetic_workflows(nb::module_&);
void bind_tidal_workflows(nb::module_&);
} // namespace lsdtt3_python
