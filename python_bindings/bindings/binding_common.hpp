#pragma once
#include <string>
#include <vector>
#include <utility>
#include <nanobind/nanobind.h>
#include <core/lsdtt_data_loader.hpp>
#include <core/lsdtt_raster_meta.hpp>
#include <cpu/lsdtt_raster_cpu.hpp>
namespace nb = nanobind;
namespace lsdtt3_python {
class Raster { public: explicit Raster(lsdtt::lsdtt_raster_cpu raster): raster_(std::move(raster)) {} static Raster read(const std::string& path) { lsdtt::lsdtt_data_loader loader; loader.open_raster(path); return Raster(lsdtt::lsdtt_raster_cpu(loader)); } void write(const std::string& path, const std::string& format) const { lsdtt::lsdtt_data_loader loader; loader.write_full(raster_.get_data(), raster_.get_meta(), path, lsdtt::resolve_raster_output_format(format)); } const lsdtt::lsdtt_raster_cpu& native() const noexcept { return raster_; } private: lsdtt::lsdtt_raster_cpu raster_; };
void bind_raster_types(nb::module_&); void bind_raster_basic_ops(nb::module_&); void bind_raster_pit_processing(nb::module_&); void bind_raster_surface_metrics(nb::module_&); void bind_raster_spectral(nb::module_&); void bind_raster_swath(nb::module_&); void bind_raster_mountain_front(nb::module_&); void bind_raster_feature_finder(nb::module_&); void bind_raster_geospatial(nb::module_&);
} // namespace lsdtt3_python
