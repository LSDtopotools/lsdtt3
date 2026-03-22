cmake_minimum_required(VERSION 3.15)
project(lsdtt3 VERSION 3.0.0 LANGUAGES CXX)

# -----------------------------------------------------------------------------
# C++ standard
# -----------------------------------------------------------------------------
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# -----------------------------------------------------------------------------
# Build type — default to Release if not specified
# -----------------------------------------------------------------------------
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

# -----------------------------------------------------------------------------
# Compiler warnings
# -----------------------------------------------------------------------------
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    add_compile_options(-Wall -Wextra -Wpedantic)
endif()

# -----------------------------------------------------------------------------
# Optional GPU backend
# -----------------------------------------------------------------------------
option(LSDTT_GPU "Enable GPU compute backend (requires CUDA)" OFF)

if(LSDTT_GPU)
    enable_language(CUDA)
    find_package(CUDAToolkit REQUIRED)
    add_compile_definitions(LSDTT_GPU_ENABLED)
    message(STATUS "GPU backend: enabled")
else()
    message(STATUS "GPU backend: disabled")
endif()

# -----------------------------------------------------------------------------
# Dependencies
# -----------------------------------------------------------------------------

# GDAL — required
find_package(GDAL REQUIRED)
if(GDAL_FOUND)
    message(STATUS "GDAL version: ${GDAL_VERSION}")
    message(STATUS "GDAL include: ${GDAL_INCLUDE_DIRS}")
else()
    message(FATAL_ERROR "GDAL not found. Install GDAL or set GDAL_DIR.")
endif()

# OpenMP — optional but strongly recommended for multi-core performance
find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    message(STATUS "OpenMP version: ${OpenMP_CXX_VERSION}")
else()
    message(WARNING "OpenMP not found — building without multi-core support.")
endif()

# pkg-config — optional helper for FFTW discovery
find_package(PkgConfig QUIET)

# FFTW — optional, enables spectral analysis routines
# This configuration expects the spectral code to use the double FFTW API:
#   fftw_*
# and not the threaded FFTW API.
set(LSDTT_HAS_FFTW OFF)

if(PkgConfig_FOUND)
    pkg_search_module(FFTW QUIET IMPORTED_TARGET fftw3)

    if(FFTW_FOUND)
        set(LSDTT_HAS_FFTW ON)
        message(STATUS "FFTW support: enabled")
        message(STATUS "  FFTW version: ${FFTW_VERSION}")
        message(STATUS "  FFTW libs   : ${FFTW_LIBRARIES}")
    else()
        message(WARNING
            "FFTW (pkg-config module 'fftw3') not found — "
            "lsdtt_raster_spectral_cpu.cpp will build as a stub and "
            "wavelet_quadratic_fit() should throw at runtime."
        )
    endif()
else()
    message(WARNING
        "PkgConfig not found — FFTW auto-detection disabled. "
        "lsdtt_raster_spectral_cpu.cpp will build as a stub and "
        "wavelet_quadratic_fit() should throw at runtime."
    )
endif()

# -----------------------------------------------------------------------------
# Output directories
# -----------------------------------------------------------------------------
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)

# -----------------------------------------------------------------------------
# lsdtt3_core
# Backend-agnostic infrastructure — always compiled
# -----------------------------------------------------------------------------
add_library(lsdtt3_core STATIC
    src/core/lsdtt_raster_base.cpp
    src/core/lsdtt_data_loader.cpp
    src/core/lsdtt_backend_factory.cpp
    src/core/lsdtt_param_parser.cpp
)

target_include_directories(lsdtt3_core
    PUBLIC
        ${CMAKE_SOURCE_DIR}/src
        ${GDAL_INCLUDE_DIRS}
)

target_link_libraries(lsdtt3_core
    PUBLIC
        ${GDAL_LIBRARIES}
)

# -----------------------------------------------------------------------------
# lsdtt3_cpu
# CPU raster class, CPU analysis free functions, and CPU flow routing
# -----------------------------------------------------------------------------
add_library(lsdtt3_cpu STATIC
    src/cpu/lsdtt_raster_cpu.cpp
    src/cpu/lsdtt_flowinfo_cpu.cpp
    src/cpu/lsdtt_flow_base.cpp
    src/cpu/lsdtt_raster_surface_metrics_cpu.cpp
    src/cpu/lsdtt_raster_flow_preprocessing_cpu.cpp
    src/cpu/lsdtt_raster_nodata_filling_cpu.cpp
    src/cpu/lsdtt_raster_spectral_cpu.cpp
    src/cpu/lsdtt_raster_synthetic_surfaces_cpu.cpp
)

target_include_directories(lsdtt3_cpu
    PUBLIC
        ${CMAKE_SOURCE_DIR}/src
)

target_link_libraries(lsdtt3_cpu
    PUBLIC
        lsdtt3_core
)

# OpenMP — CPU-only dependency
if(OpenMP_CXX_FOUND)
    target_link_libraries(lsdtt3_cpu PUBLIC OpenMP::OpenMP_CXX)
endif()

# FFTW — CPU-only dependency
if(LSDTT_HAS_FFTW)
    target_link_libraries(lsdtt3_cpu PUBLIC PkgConfig::FFTW)
    target_compile_definitions(lsdtt3_cpu PRIVATE LSDTT_HAS_FFTW=1)
endif()

# -----------------------------------------------------------------------------
# GPU backend — only compiled if LSDTT_GPU=ON
# -----------------------------------------------------------------------------
if(LSDTT_GPU)
    add_library(lsdtt3_gpu STATIC
        src/gpu/lsdtt_raster_gpu.cu
    )

    target_include_directories(lsdtt3_gpu
        PUBLIC
            ${CMAKE_SOURCE_DIR}/src
    )

    target_link_libraries(lsdtt3_gpu
        PUBLIC
            lsdtt3_core
            CUDA::cudart
    )
endif()

# -----------------------------------------------------------------------------
# CLI executables
# -----------------------------------------------------------------------------

# Raster preprocessing CLI
add_executable(lsdtt-raster-preprocessing
    cli/lsdtt-raster-preprocessing.cpp
)

target_link_libraries(lsdtt-raster-preprocessing
    PRIVATE
        lsdtt3_cpu
        lsdtt3_core
)

# Synthetic raster generation and wavelet test CLI
add_executable(lsdtt-raster-synthetic-tests
    cli/lsdtt-raster-synthetic-tests.cpp
)

target_link_libraries(lsdtt-raster-synthetic-tests
    PRIVATE
        lsdtt3_cpu
        lsdtt3_core
)

# D8 flow routing test and diagnostic raster output CLI
add_executable(lsdtt-test-flowinfo
    cli/lsdtt-test-flowinfo.cpp
)

target_link_libraries(lsdtt-test-flowinfo
    PRIVATE
        lsdtt3_cpu
        lsdtt3_core
)

# -----------------------------------------------------------------------------
# Installation
# -----------------------------------------------------------------------------
install(TARGETS
    lsdtt-raster-preprocessing
    lsdtt-raster-synthetic-tests
    lsdtt-test-flowinfo
    RUNTIME DESTINATION bin
)

# -----------------------------------------------------------------------------
# Summary
# -----------------------------------------------------------------------------
message(STATUS "")
message(STATUS "lsdtt3 configuration summary:")
message(STATUS "  Version     : ${PROJECT_VERSION}")
message(STATUS "  Build type  : ${CMAKE_BUILD_TYPE}")
message(STATUS "  C++ standard: ${CMAKE_CXX_STANDARD}")
message(STATUS "  GPU backend : ${LSDTT_GPU}")
message(STATUS "  OpenMP      : ${OpenMP_CXX_FOUND}")
message(STATUS "  FFTW        : ${LSDTT_HAS_FFTW}")
message(STATUS "  GDAL        : ${GDAL_VERSION}")
message(STATUS "  Output dir  : ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}")
message(STATUS "  Libraries   : lsdtt3_core lsdtt3_cpu")
message(STATUS "")
