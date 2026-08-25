nanobind_add_module(_core
    module.cpp
    raster_types.cpp
    raster_basic_ops.cpp
    raster_pit_processing.cpp
    raster_surface_metrics.cpp
    raster_geospatial.cpp
)

target_link_libraries(_core PRIVATE lsdtt3::runtime)
target_compile_features(_core PRIVATE cxx_std_17)
set_target_properties(_core PROPERTIES
    BUILD_RPATH "${LSDTT_PYTHON_RUNTIME_PATH}"
    INSTALL_RPATH "${LSDTT_PYTHON_RUNTIME_PATH}"
)
install(TARGETS _core LIBRARY DESTINATION lsdtt3)
