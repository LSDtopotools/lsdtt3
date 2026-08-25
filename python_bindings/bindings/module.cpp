#include "binding_common.hpp"
NB_MODULE(_core, module) { module.doc() = "Python bindings for lsdtt3 raster operations."; lsdtt3_python::bind_raster_types(module); lsdtt3_python::bind_raster_basic_ops(module); lsdtt3_python::bind_raster_pit_processing(module); lsdtt3_python::bind_raster_surface_metrics(module); lsdtt3_python::bind_raster_geospatial(module); }
