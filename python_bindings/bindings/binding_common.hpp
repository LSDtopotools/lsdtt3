#pragma once
#include <string>
#include <vector>
#include <utility>
#include <memory>
#include <nanobind/nanobind.h>
#include <core/lsdtt_data_loader.hpp>
#include <core/lsdtt_raster_meta.hpp>
#include <cpu/lsdtt_raster_cpu.hpp>
#include <cpu/lsdtt_flowinfo_cpu.hpp>
#include <cpu/lsdtt_channel_scope_cpu.hpp>
#include <cpu/lsdtt_junction_network_cpu.hpp>
#include <cpu/lsdtt_flowgraph_cpu.hpp>
#include <cpu/lsdtt_chi_extraction_cpu.hpp>
#include <cpu/lsdtt_chi_network_analysis_cpu.hpp>
namespace nb = nanobind;
namespace lsdtt3_python {
class Raster { public: explicit Raster(lsdtt::lsdtt_raster_cpu raster): raster_(std::move(raster)) {} static Raster read(const std::string& path) { lsdtt::lsdtt_data_loader loader; loader.open_raster(path); return Raster(lsdtt::lsdtt_raster_cpu(loader)); } void write(const std::string& path, const std::string& format) const { lsdtt::lsdtt_data_loader loader; loader.write_full(raster_.get_data(), raster_.get_meta(), path, lsdtt::resolve_raster_output_format(format)); } const lsdtt::lsdtt_raster_cpu& native() const noexcept { return raster_; } private: lsdtt::lsdtt_raster_cpu raster_; };
class FlowInfo { public: explicit FlowInfo(lsdtt::lsdtt_flowinfo flow): flow_(std::make_shared<lsdtt::lsdtt_flowinfo>(std::move(flow))) {} FlowInfo(const Raster& dem, const std::string& boundary_conditions, bool interior_nodata_is_outlet, nb::object forced_outlets, bool allow_pits_as_outlets, int n_threads); const lsdtt::lsdtt_flowinfo& native() const noexcept { return *flow_; } std::shared_ptr<lsdtt::lsdtt_flowinfo> shared_native() const noexcept { return flow_; } private: std::shared_ptr<lsdtt::lsdtt_flowinfo> flow_; };
class ChannelScope { public: ChannelScope(lsdtt::channel_scope_result result, std::shared_ptr<lsdtt::lsdtt_flowinfo> flow): result_(std::move(result)), flow_(std::move(flow)) {} const lsdtt::channel_scope_result& native() const noexcept { return result_; } const lsdtt::lsdtt_flowinfo& flow() const noexcept { return *flow_; } std::shared_ptr<lsdtt::lsdtt_flowinfo> shared_flow() const noexcept { return flow_; } private: lsdtt::channel_scope_result result_; std::shared_ptr<lsdtt::lsdtt_flowinfo> flow_; };
class ChannelNetwork { public: ChannelNetwork(std::shared_ptr<lsdtt::lsdtt_flowinfo> flow, std::vector<std::int64_t> sources): flow_(std::move(flow)), network_(std::make_unique<lsdtt::lsdtt_junction_network_cpu>(*flow_, sources)) {} const lsdtt::lsdtt_junction_network_cpu& native() const noexcept { return *network_; } const lsdtt::lsdtt_flowinfo& flow() const noexcept { return *flow_; } std::shared_ptr<lsdtt::lsdtt_flowinfo> shared_flow() const noexcept { return flow_; } private: std::shared_ptr<lsdtt::lsdtt_flowinfo> flow_; std::unique_ptr<lsdtt::lsdtt_junction_network_cpu> network_; };
class ChiCollection { public: ChiCollection(lsdtt::lsdtt_chi_collection_result result, std::shared_ptr<lsdtt::lsdtt_flowinfo> flow): result_(std::move(result)), flow_(std::move(flow)) {} const lsdtt::lsdtt_chi_collection_result& native() const noexcept { return result_; } const lsdtt::lsdtt_flowinfo& flow() const noexcept { return *flow_; } std::shared_ptr<lsdtt::lsdtt_flowinfo> shared_flow() const noexcept { return flow_; } private: lsdtt::lsdtt_chi_collection_result result_; std::shared_ptr<lsdtt::lsdtt_flowinfo> flow_; };
class ChiFit { public: ChiFit(lsdtt::lsdtt_chi_network_fit_collection_result result, std::shared_ptr<lsdtt::lsdtt_flowinfo> flow): result_(std::move(result)), flow_(std::move(flow)) {} const lsdtt::lsdtt_chi_network_fit_collection_result& native() const noexcept { return result_; } const lsdtt::lsdtt_flowinfo& flow() const noexcept { return *flow_; } private: lsdtt::lsdtt_chi_network_fit_collection_result result_; std::shared_ptr<lsdtt::lsdtt_flowinfo> flow_; };
class FlowGraph { public: FlowGraph(const Raster& dem, const std::string& method, const std::string& boundary_conditions, bool interior_nodata_is_outlet, float mfd_exponent, int n_threads); const lsdtt::lsdtt_flowgraph_cpu& native() const noexcept { return graph_; } private: lsdtt::lsdtt_flowgraph_cpu graph_; };
void bind_raster_types(nb::module_&); void bind_raster_basic_ops(nb::module_&); void bind_raster_pit_processing(nb::module_&); void bind_raster_surface_metrics(nb::module_&); void bind_raster_spectral(nb::module_&); void bind_raster_swath(nb::module_&); void bind_raster_mountain_front(nb::module_&); void bind_raster_feature_finder(nb::module_&); void bind_raster_flat_surfaces(nb::module_&); void bind_raster_marsh(nb::module_&); void bind_raster_ridgeline(nb::module_&); void bind_valley_helpers(nb::module_&); void bind_flow_info(nb::module_&); void bind_flow_graph(nb::module_&); void bind_flow_indices(nb::module_&); void bind_basins(nb::module_&, nb::class_<FlowInfo>&); void bind_raster_geospatial(nb::module_&);
} // namespace lsdtt3_python
