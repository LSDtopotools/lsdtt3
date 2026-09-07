nanobind_add_module(_core
    module.cpp
    raster_types.cpp
    raster_basic_ops.cpp
    raster_resampling.cpp
    raster_hierarchy.cpp
    raster_pit_processing.cpp
    raster_surface_metrics.cpp
    raster_spectral.cpp
    raster_swath.cpp
    raster_mountain_front.cpp
    raster_feature_finder.cpp
    raster_flat_surfaces.cpp
    raster_marsh.cpp
    raster_ridgeline.cpp
    valley_helpers.cpp
    flow_info.cpp
    flow_products.cpp
    channel_heads.cpp
    flow_network.cpp
    flow_graph.cpp
    hydraulic_routing.cpp
    flow_indices.cpp
    basins.cpp
    divide_migration.cpp
    raster_geospatial.cpp
    raster_fetch.cpp
    cosmo_column.cpp
    cosmo_raster.cpp
    cosmo_basin.cpp
    statistical_comparison.cpp
    vector_hydrology.cpp
    internal_drainage.cpp
    hillslope_transport.cpp
    landscape_model.cpp
    sediment_routing.cpp
    landslides.cpp
    synthetic_workflows.cpp
    tidal_workflows.cpp
)

target_link_libraries(_core PRIVATE lsdtt3::runtime)
target_compile_features(_core PRIVATE cxx_std_17)
set_target_properties(_core PROPERTIES
    BUILD_RPATH "${LSDTT_PYTHON_RUNTIME_PATH}"
    INSTALL_RPATH "${LSDTT_PYTHON_RUNTIME_PATH}"
)
install(TARGETS _core LIBRARY DESTINATION lsdtt3)
