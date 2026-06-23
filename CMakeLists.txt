cmake_minimum_required(VERSION 3.15)
project(lsdtt3 VERSION 0.2.0 LANGUAGES CXX)

# =============================================================================
# C++ standard
# =============================================================================
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# =============================================================================
# Build type default
# =============================================================================
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

# =============================================================================
# Compiler warnings — all targets inherit these
# =============================================================================
add_compile_options(
    -Wall
    -Wextra
    -Wpedantic
    -Wno-unused-parameter   # OGR headers trigger this; suppress globally
)

# =============================================================================
# Options
# =============================================================================
option(LSDTT_GPU "Build GPU backends (requires CUDA)" OFF)

# =============================================================================
# Required: GDAL
# =============================================================================
find_package(GDAL REQUIRED)
if(NOT GDAL_FOUND)
    message(FATAL_ERROR "GDAL not found. Install libgdal-dev and try again.")
endif()
message(STATUS "GDAL version : ${GDAL_VERSION}")
message(STATUS "GDAL include : ${GDAL_INCLUDE_DIRS}")
message(STATUS "GDAL libs    : ${GDAL_LIBRARIES}")

# gdal_utils.h (needed by lsdtt_raster_reproject for GDALWarpAppOptions)
# is part of the standard GDAL dev package but lives in a subdirectory on
# some distributions.  Add the GDAL include path; CMake's GDAL finder
# already puts it in GDAL_INCLUDE_DIRS.

# =============================================================================
# Optional: OpenMP
# =============================================================================
find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    message(STATUS "OpenMP found — parallel CPU operations enabled")
else()
    message(WARNING "OpenMP not found — CPU operations will be single-threaded")
endif()

# =============================================================================
# Optional: FFTW3 (double, non-threaded)
# =============================================================================
find_package(PkgConfig)
set(LSDTT_HAS_FFTW OFF)
if(PkgConfig_FOUND)
    pkg_check_modules(FFTW3 fftw3)
    if(FFTW3_FOUND)
        set(LSDTT_HAS_FFTW ON)
        message(STATUS "FFTW3 found — spectral analysis enabled")
    else()
        message(WARNING "FFTW3 not found — spectral analysis will throw at runtime")
    endif()
endif()

# =============================================================================
# Optional: netCDF-C (required by lsdtt_crn_parameters ERA40 / paleomagnetic)
# =============================================================================
set(LSDTT_HAS_NETCDF OFF)
find_path(NETCDF_INCLUDE_DIR NAMES netcdf.h)
find_library(NETCDF_LIBRARY  NAMES netcdf)
if(NETCDF_INCLUDE_DIR AND NETCDF_LIBRARY)
    set(LSDTT_HAS_NETCDF ON)
    message(STATUS "netCDF-C found: ${NETCDF_LIBRARY} — ERA40 / paleomagnetic correction enabled")
else()
    message(STATUS "netCDF-C not found — ERA40 atmosphere and paleomagnetic correction disabled")
endif()

# =============================================================================
# Optional: OpenCV (core + imgproc — required by lsdtt_raster_opencv_cpu)
# =============================================================================
option(LSDTT_OPENCV "Build OpenCV image-analysis routines" ON)
set(LSDTT_HAS_OPENCV OFF)
if(LSDTT_OPENCV)
    # Search /usr/local first (locally-built OpenCV), then system paths.
    # Not QUIET so cmake prints the reason when it cannot find OpenCV.
    find_package(OpenCV
        COMPONENTS core imgproc
        HINTS /usr/local/lib/cmake/opencv4
              /usr/lib/x86_64-linux-gnu/cmake/opencv4
    )
    if(OpenCV_FOUND)
        set(LSDTT_HAS_OPENCV ON)
        message(STATUS "OpenCV found (${OpenCV_VERSION}) — image analysis enabled")
    else()
        message(STATUS "OpenCV not found — lsdtt_raster_opencv_cpu will not be built "
                       "(install libopencv-dev or pass -DLSDTT_OPENCV=OFF to silence)")
    endif()
endif()

# =============================================================================
# Optional: Eigen3 (required by lsdtt_chi_ksn_inversion_cpu)
# =============================================================================
option(LSDTT_EIGEN "Build Eigen-based inversion and sparse linear-algebra routines" ON)
set(LSDTT_HAS_EIGEN OFF)

if(LSDTT_EIGEN)
    find_package(Eigen3 QUIET NO_MODULE)
    if(Eigen3_FOUND)
        set(LSDTT_HAS_EIGEN ON)
        message(STATUS "Eigen3 found — chi ksn inversion enabled")
    else()
        message(WARNING "Eigen3 not found — lsdtt_chi_ksn_inversion_cpu will not be built "
                        "(install libeigen3-dev or pass -DLSDTT_EIGEN=OFF to silence)")
    endif()
endif()




# =============================================================================
# Optional: CUDA
# =============================================================================
if(LSDTT_GPU)
    enable_language(CUDA)
    find_package(CUDAToolkit REQUIRED)
    message(STATUS "CUDA toolkit: ${CUDAToolkit_VERSION}")
endif()

# =============================================================================
# lsdtt3_core  (static library)
#
# Abstract interfaces, factories, shared types, I/O.
# Depends only on GDAL.  No CPU arrays.  No CUDA.
# =============================================================================
set(LSDTT3_CORE_SOURCES
    src/core/lsdtt_raster_base.cpp
    src/core/lsdtt_data_loader.cpp
    src/core/lsdtt_backend_factory.cpp
    src/core/lsdtt_param_parser.cpp
    src/core/lsdtt_crs_suggest.cpp
    src/core/lsdtt_raster_reproject.cpp
    src/core/lsdtt_vector_utils.cpp
    src/core/lsdtt_raster_warp.cpp
    src/core/lsdtt_point_input.cpp
    src/core/lsdtt_stat_csv.cpp
)

add_library(lsdtt3_core STATIC ${LSDTT3_CORE_SOURCES})

target_include_directories(lsdtt3_core
    PUBLIC
        ${CMAKE_SOURCE_DIR}/src
        ${CMAKE_SOURCE_DIR}/src/core
        ${GDAL_INCLUDE_DIRS}
)

target_link_libraries(lsdtt3_core
    PUBLIC
        ${GDAL_LIBRARIES}
)

# =============================================================================
# lsdtt3_utils  (static library)
#
# Generic mathematical / utility code that is not GDAL-specific and not tied
# to a particular backend implementation.
# Header-only utilities in src/utils/ remain header-only; only compiled utility
# sources are listed here.
# =============================================================================
set(LSDTT3_UTILS_SOURCES
    src/utils/lsdtt_statistics.cpp
    src/utils/lsdtt_histogram.cpp
    src/utils/lsdtt_statistical_tests.cpp
    src/utils/lsdtt_distributions.cpp
    src/utils/lsdtt_splines.cpp
    src/utils/lsdtt_regression.cpp
    src/utils/lsdtt_binary_array2d_tracing.cpp
    src/utils/lsdtt_signal1d.cpp
)

add_library(lsdtt3_utils STATIC ${LSDTT3_UTILS_SOURCES})

target_include_directories(lsdtt3_utils
    PUBLIC
        ${CMAKE_SOURCE_DIR}/src
        ${CMAKE_SOURCE_DIR}/src/utils
)


# =============================================================================
# lsdtt3_cpu  (static library)
#
# All CPU-side implementations.
# Depends on lsdtt3_core, OpenMP (optional), FFTW (optional).
# =============================================================================
set(LSDTT3_CPU_SOURCES
    src/cpu/lsdtt_raster_cpu.cpp
    src/cpu/lsdtt_raster_algebra.cpp
    src/cpu/lsdtt_raster_surface_metrics_cpu.cpp
    src/cpu/lsdtt_surface_texture_cpu.cpp
    src/cpu/lsdtt_raster_surface_fitting_cpu.cpp
    src/cpu/lsdtt_raster_flow_preprocessing_cpu.cpp
    src/cpu/lsdtt_raster_nodata_filling_cpu.cpp
    src/cpu/lsdtt_raster_spectral_cpu.cpp
    src/cpu/lsdtt_raster_synthetic_surfaces_cpu.cpp
    src/cpu/lsdtt_synthetic_network_cpu.cpp
    src/cpu/lsdtt_synthetic_domain_cpu.cpp
    src/cpu/lsdtt_flowinfo_cpu.cpp
    src/cpu/lsdtt_flow_base.cpp
    src/cpu/lsdtt_junction_network_cpu.cpp
    src/cpu/lsdtt_junction_network_base.cpp
    src/cpu/lsdtt_raster_feature_finder_cpu.cpp
    src/cpu/lsdtt_mountain_front_cpu.cpp
    src/cpu/lsdtt_depression_hierarchy_cpu.cpp
    src/cpu/lsdtt_fill_spill_merge_cpu.cpp
    src/cpu/lsdtt_flow_snapping_cpu.cpp
    src/cpu/lsdtt_basin_cpu.cpp
    src/cpu/lsdtt_junction_basin_cpu.cpp
    src/cpu/lsdtt_junction_basin_manager_cpu.cpp
    src/cpu/lsdtt_basin_outline_cpu.cpp
    src/cpu/lsdtt_flat_surface_finder.cpp
    src/cpu/lsdtt_slope_break_finder_cpu.cpp
    src/cpu/lsdtt_flat_surface_plane_analysis.cpp
    src/cpu/lsdtt_vector_rasterise_cpu.cpp
    src/cpu/lsdtt_proximity_cpu.cpp
    src/cpu/lsdtt_swath_cpu.cpp
    src/cpu/lsdtt_streamlines_cpu.cpp
    src/cpu/lsdtt_flowgraph_cpu.cpp
    src/cpu/lsdtt_raster_vectorise_cpu.cpp
    src/cpu/lsdtt_raster_mosaic_cpu.cpp
    src/cpu/lsdtt_connectivity_cpu.cpp
    src/cpu/lsdtt_hydraulic_graph_cpu.cpp
    src/cpu/lsdtt_chi_extraction_cpu.cpp
    src/cpu/lsdtt_chi_network_analysis_cpu.cpp
    src/cpu/lsdtt_chi_disorder_cpu.cpp
    src/cpu/lsdtt_chi_data_maps_cpu.cpp
    src/cpu/lsdtt_channel_network_analysis_cpu.cpp
    src/cpu/lsdtt_strahler_tokunaga_analysis_cpu.cpp
    src/cpu/lsdtt_raster_skeleton_graph_cpu.cpp
    src/cpu/lsdtt_raster_binary_mask_vectorise_cpu.cpp
    src/cpu/lsdtt_raster_binary_mask_cpu.cpp
    src/cpu/lsdtt_raster_linear_features_cpu.cpp
    src/cpu/lsdtt_burn_streams_cpu.cpp
    src/cpu/lsdtt_least_cost_path.cpp
    src/cpu/lsdtt_terrace_edge_tracer_cpu.cpp
    src/cpu/lsdtt_raster_contours_cpu.cpp
    src/cpu/lsdtt_labelled_raster_stats_cpu.cpp
    src/cpu/lsdtt_cross_sections_cpu.cpp
    src/cpu/lsdtt_crn_parameters.cpp
    src/cpu/lsdtt_crn_column.cpp
    src/cpu/lsdtt_cosmo_raster_cpu.cpp
    src/cpu/lsdtt_vector_to_junction_network.cpp
    src/cpu/lsdtt_vector_ops_cpu.cpp
    src/cpu/lsdtt_vector_line_repair_cpu.cpp
    src/cpu/lsdtt_vector_network_d8_clean_cpu.cpp
    src/cpu/lsdtt_burn_vectormap_cpu.cpp
    src/cpu/lsdtt_water_area_centreline_cpu.cpp
    src/cpu/lsdtt_marsh_platform_cpu.cpp
    src/cpu/lsdtt_tidal_channel_cpu.cpp
    src/cpu/lsdtt_tidal_skeleton_metrics_cpu.cpp
    src/cpu/lsdtt_landscape_drape_cpu.cpp
    src/cpu/lsdtt_stream_power_cpu.cpp
    src/cpu/lsdtt_hillslope_transport_cpu.cpp
    src/cpu/lsdtt_lem_cpu.cpp
    src/cpu/lsdtt_steady_state_cpu.cpp
    src/cpu/lsdtt_hillslope_steady_state_cpu.cpp
    src/cpu/lsdtt_sediment_transport_cpu.cpp
    src/cpu/lsdtt_giuh_cpu.cpp
    src/cpu/lsdtt_side_tributary_basins_cpu.cpp
    src/cpu/lsdtt_ridgeline_cpu.cpp
)

if(LSDTT_HAS_OPENCV)
    list(APPEND LSDTT3_CPU_SOURCES src/cpu/lsdtt_raster_opencv_cpu.cpp)
endif()

if(LSDTT_HAS_EIGEN)
    list(APPEND LSDTT3_CPU_SOURCES
        src/cpu/lsdtt_chi_ksn_inversion_cpu.cpp
    )
endif()


add_library(lsdtt3_cpu STATIC ${LSDTT3_CPU_SOURCES})

target_include_directories(lsdtt3_cpu
    PUBLIC
        ${CMAKE_SOURCE_DIR}/src
        ${CMAKE_SOURCE_DIR}/src/core
        ${CMAKE_SOURCE_DIR}/src/cpu
        ${GDAL_INCLUDE_DIRS}
)

target_link_libraries(lsdtt3_cpu
    PUBLIC
        lsdtt3_core
        lsdtt3_utils
)

if(OpenMP_CXX_FOUND)
    target_link_libraries(lsdtt3_cpu PUBLIC OpenMP::OpenMP_CXX)
endif()

if(LSDTT_HAS_FFTW)
    # Give FFTW to the CPU library
    target_compile_definitions(lsdtt3_cpu PUBLIC LSDTT_HAS_FFTW=1)
    target_include_directories(lsdtt3_cpu PUBLIC ${FFTW3_INCLUDE_DIRS})
    target_link_libraries(lsdtt3_cpu PUBLIC ${FFTW3_LIBRARIES})

    # Give FFTW to the Utils library (where signal1d code lives)
    target_compile_definitions(lsdtt3_utils PUBLIC LSDTT_HAS_FFTW=1)
    target_include_directories(lsdtt3_utils PUBLIC ${FFTW3_INCLUDE_DIRS})
    target_link_libraries(lsdtt3_utils PUBLIC ${FFTW3_LIBRARIES})
endif()

if(LSDTT_HAS_NETCDF)
    target_compile_definitions(lsdtt3_cpu PUBLIC LSDTT_HAS_NETCDF=1)
    target_include_directories(lsdtt3_cpu PUBLIC ${NETCDF_INCLUDE_DIR})
    target_link_libraries(lsdtt3_cpu PUBLIC ${NETCDF_LIBRARY})
endif()

if(LSDTT_HAS_OPENCV)
    # OpenCV's CMake targets propagate their own include directories, so no
    # manual target_include_directories call is needed here.
    target_compile_definitions(lsdtt3_cpu PUBLIC LSDTT_HAS_OPENCV=1)
    target_link_libraries(lsdtt3_cpu PUBLIC ${OpenCV_LIBS})
endif()

if(LSDTT_HAS_EIGEN)
    target_compile_definitions(lsdtt3_cpu PUBLIC LSDTT_HAS_EIGEN=1)
    target_link_libraries(lsdtt3_cpu PUBLIC Eigen3::Eigen)
endif()



# =============================================================================
# lsdtt3_gpu  (static library, LSDTT_GPU=ON only)
#
# GPU compute kernels for the selective-dispatch path (currently the surface-
# fitting kernels). Self-contained: a plain-C++ POD interface (lsdtt_surface_gpu
# .hpp) backed by CUDA implementations, with no GDAL or lsdtt3_core dependency.
#
# This is the GPU *compute kernel* axis (LSDTT_GPU_KERNELS), distinct from the
# future GPU *raster-engine backend* used by make_raster_engine().
# =============================================================================
if(LSDTT_GPU)
    set(LSDTT3_GPU_SOURCES
        src/gpu/lsdtt_surface_gpu.cu
    )
    add_library(lsdtt3_gpu STATIC ${LSDTT3_GPU_SOURCES})
    target_include_directories(lsdtt3_gpu
        PUBLIC
            ${CMAKE_SOURCE_DIR}/src
    )
    target_link_libraries(lsdtt3_gpu
        PUBLIC
            CUDA::cudart
    )

    # The dispatch seam lives in lsdtt3_cpu (e.g. direct_polyfit). Enable the
    # GPU branch there and link the kernel library so its symbols resolve and
    # propagate to every CLI that links lsdtt3_cpu.
    target_compile_definitions(lsdtt3_cpu PUBLIC LSDTT_GPU_KERNELS=1)
    target_link_libraries(lsdtt3_cpu PUBLIC lsdtt3_gpu)
endif()

# =============================================================================
# Helper macro — add a CLI executable
#
# Usage: lsdtt_cli(target_name source_file)
# Links against lsdtt3_cpu and lsdtt3_core automatically.
# =============================================================================
macro(lsdtt_cli TARGET SOURCE)
    add_executable(${TARGET} ${SOURCE})
    target_include_directories(${TARGET}
        PRIVATE
            ${CMAKE_SOURCE_DIR}/src
            ${CMAKE_SOURCE_DIR}/src/core
            ${CMAKE_SOURCE_DIR}/src/cpu
            ${GDAL_INCLUDE_DIRS}
    )
    target_link_libraries(${TARGET}
        PRIVATE
            lsdtt3_cpu
            lsdtt3_core
    )
    if(LSDTT_GPU)
        target_link_libraries(${TARGET} PRIVATE lsdtt3_gpu)
    endif()
    # Install to bin/
    install(TARGETS ${TARGET} RUNTIME DESTINATION bin)
endmacro()

# =============================================================================
# CLI executables
# =============================================================================
lsdtt_cli(lsdtt-raster-preprocessing   cli/lsdtt-raster-preprocessing.cpp)
lsdtt_cli(lsdtt-raster-synthetic-tests cli/lsdtt-raster-synthetic-tests.cpp)
lsdtt_cli(lsdtt-synthetic-networks     cli/lsdtt-synthetic-networks.cpp)
lsdtt_cli(lsdtt-surface-metrics        cli/lsdtt-surface-metrics.cpp)
lsdtt_cli(lsdtt-channel-extraction     cli/lsdtt-channel-extraction.cpp)
lsdtt_cli(lsdtt-drainage-analysis      cli/lsdtt-drainage-analysis.cpp)
lsdtt_cli(lsdtt-feature-finder         cli/lsdtt-feature-finder.cpp)
lsdtt_cli(lsdtt-mountain-front         cli/lsdtt-mountain-front.cpp)
lsdtt_cli(lsdtt-basin-stats            cli/lsdtt-basin-stats.cpp)
lsdtt_cli(lsdtt-basin-extraction       cli/lsdtt-basin-extraction.cpp)
lsdtt_cli(lsdtt-flat-surface-finder    cli/lsdtt-flat-surface-finder.cpp)
lsdtt_cli(lsdtt-hydraulic-routing      cli/lsdtt-hydraulic-routing.cpp)
lsdtt_cli(lsdtt-channel-network-analysis cli/lsdtt-channel-network-analysis.cpp)
lsdtt_cli(lsdtt-cross-sections         cli/lsdtt-cross-sections.cpp)
lsdtt_cli(lsdtt-stream-burning         cli/lsdtt-stream-burning.cpp)
lsdtt_cli(lsdtt-water-line-repair      cli/lsdtt-water-line-repair.cpp)
lsdtt_cli(lsdtt-water-centreline       cli/lsdtt-water-centreline.cpp)
lsdtt_cli(lsdtt-cosmo-tool             cli/lsdtt-cosmo-tool.cpp)
lsdtt_cli(lsdtt-cosmo-prediction       cli/lsdtt-cosmo-prediction.cpp)
lsdtt_cli(lsdtt-crn-column-tool        cli/lsdtt-crn-column-tool.cpp)
lsdtt_cli(lsdtt-marsh-topography       cli/lsdtt-marsh-topography.cpp)
lsdtt_cli(lsdtt-linear-stats           cli/lsdtt-linear-stats.cpp)
lsdtt_cli(lsdtt-bearing-comparison     cli/lsdtt-bearing-comparison.cpp)
lsdtt_cli(lsdtt-swath-profile          cli/lsdtt-swath-profile.cpp)
lsdtt_cli(lsdtt-lem                    cli/lsdtt-lem.cpp)
lsdtt_cli(lsdtt-sediment-routing       cli/lsdtt-sediment-routing.cpp)
lsdtt_cli(lsdtt-ridgeline-analysis     cli/lsdtt-ridgeline-analysis.cpp)

if(LSDTT_HAS_FFTW)
    lsdtt_cli(lsdtt-spectral-analysis cli/lsdtt-spectral-analysis.cpp)
endif()

if(LSDTT_HAS_EIGEN)
    lsdtt_cli(lsdtt-chi-analysis cli/lsdtt-chi-analysis.cpp)
endif()



# =============================================================================
# Install libraries and headers (optional — useful for packaging)
# =============================================================================
install(TARGETS lsdtt3_core lsdtt3_utils lsdtt3_cpu
    ARCHIVE DESTINATION lib
)


install(DIRECTORY src/core/ DESTINATION include/lsdtt3
    FILES_MATCHING PATTERN "*.hpp"
)

install(DIRECTORY src/utils/ DESTINATION include/lsdtt3/utils
    FILES_MATCHING PATTERN "*.hpp"
)


# =============================================================================
# Summary
# =============================================================================
message(STATUS "")
message(STATUS "=== lsdtt3 build configuration ===")
message(STATUS "  Build type  : ${CMAKE_BUILD_TYPE}")
message(STATUS "  C++ std     : ${CMAKE_CXX_STANDARD}")
message(STATUS "  GDAL        : ${GDAL_VERSION}")
message(STATUS "  OpenMP      : ${OpenMP_CXX_FOUND}")
message(STATUS "  FFTW3       : ${LSDTT_HAS_FFTW}")
message(STATUS "  OpenCV      : ${LSDTT_HAS_OPENCV}")
message(STATUS "  GPU (CUDA)  : ${LSDTT_GPU}")
message(STATUS "  Eigen3      : ${LSDTT_HAS_EIGEN}")
message(STATUS "===================================")
message(STATUS "")
