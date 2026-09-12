cmake_minimum_required(VERSION 3.18)

# This is deliberately a standalone build. It compiles the LSDTopoTools
# implementation directly from ../src and creates no targets in the root
# project's CMake build.
project(lsdtt3_python LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

set(LSDTT3_SOURCE_DIR "${CMAKE_CURRENT_LIST_DIR}/..")
get_filename_component(LSDTT3_SOURCE_DIR "${LSDTT3_SOURCE_DIR}" REALPATH)
set(LSDTT3_SRC_DIR "${LSDTT3_SOURCE_DIR}/src")

option(LSDTT_PYTHON_OPENMP "Enable OpenMP in the bundled LSDTT runtime" ON)
option(LSDTT_PYTHON_FFTW "Enable FFTW spectral routines when available" ON)
option(LSDTT_PYTHON_NETCDF "Enable netCDF routines when available" ON)
option(LSDTT_PYTHON_OPENCV "Enable OpenCV image-analysis routines when available" ON)
option(LSDTT_PYTHON_EIGEN "Enable Eigen inversion routines when available" ON)
option(LSDTT_PYTHON_RELEASE_WHEEL
    "Stage GDAL and PROJ runtime data for a self-contained release wheel" OFF)
option(LSDTT_PYTHON_CLI
    "Build the selected vanilla LSDTopoTools command-line drivers for the wheel" ON)
option(LSDTT_PYTHON_COMPATIBILITY_PROBE
    "Build private keyword-migration probes for API contract tests" OFF)

# These optional overrides make a reproducible release build independent of
# how CMake located GDAL and PROJ on the build host.
set(LSDTT_PYTHON_GDAL_DATA_DIR "" CACHE PATH
    "Directory containing GDAL runtime data (for example gdalvrt.xsd)")
set(LSDTT_PYTHON_PROJ_DATA_DIR "" CACHE PATH
    "Directory containing PROJ runtime data (containing proj.db)")

# A direct CMake invocation does not automatically add a pip-installed
# nanobind package to CMAKE_PREFIX_PATH. Ask the active Python interpreter
# for it; scikit-build-core already supplies the same package in wheel builds.
find_package(Python COMPONENTS Interpreter Development.Module REQUIRED)
execute_process(
    COMMAND "${Python_EXECUTABLE}" -m nanobind --cmake_dir
    RESULT_VARIABLE NANOBIND_CMAKE_DIR_RESULT
    OUTPUT_VARIABLE NANOBIND_CMAKE_DIR
    OUTPUT_STRIP_TRAILING_WHITESPACE
    ERROR_QUIET
)
if(NANOBIND_CMAKE_DIR_RESULT EQUAL 0)
    list(PREPEND CMAKE_PREFIX_PATH "${NANOBIND_CMAKE_DIR}")
endif()
find_package(nanobind CONFIG REQUIRED)
if(nanobind_VERSION VERSION_LESS "3.0" OR NOT nanobind_VERSION VERSION_LESS "3.1")
    message(FATAL_ERROR "lsdtt3 Python bindings require nanobind >=3.0,<3.1")
endif()
find_package(GDAL REQUIRED)

if(LSDTT_PYTHON_OPENMP AND NOT WIN32)
    find_package(OpenMP QUIET)
endif()

if(LSDTT_PYTHON_FFTW)
    # Prefer pkg-config: it works with both distro packages and the Conda
    # package, whose CMake config can lack a Release import location.
    find_package(PkgConfig QUIET)
    if(PkgConfig_FOUND)
        pkg_check_modules(FFTW3 QUIET IMPORTED_TARGET fftw3)
        if(FFTW3_FOUND)
            set(LSDTT_PYTHON_FFTW_TARGET PkgConfig::FFTW3)
        endif()
    endif()
    if(NOT LSDTT_PYTHON_FFTW_TARGET)
        find_package(FFTW3 CONFIG QUIET)
        if(FFTW3_FOUND AND TARGET FFTW3::fftw3)
            set(LSDTT_PYTHON_FFTW_TARGET FFTW3::fftw3)
        endif()
    endif()
endif()

if(LSDTT_PYTHON_NETCDF)
    find_path(NETCDF_INCLUDE_DIR NAMES netcdf.h)
    find_library(NETCDF_LIBRARY NAMES netcdf)
endif()

if(LSDTT_PYTHON_OPENCV)
    find_package(OpenCV QUIET COMPONENTS core imgproc)
endif()

if(LSDTT_PYTHON_EIGEN)
    find_package(Eigen3 QUIET NO_MODULE)
endif()

# Keep this source selection aligned with the root CMake build. The two
# conditional implementations are only compiled after their dependency has
# been found; all other core, utility, and CPU sources belong to the runtime.
file(GLOB LSDTT3_CORE_SOURCES CONFIGURE_DEPENDS "${LSDTT3_SRC_DIR}/core/*.cpp")
file(GLOB LSDTT3_UTILS_SOURCES CONFIGURE_DEPENDS "${LSDTT3_SRC_DIR}/utils/*.cpp")
file(GLOB LSDTT3_CPU_SOURCES CONFIGURE_DEPENDS "${LSDTT3_SRC_DIR}/cpu/*.cpp")
list(REMOVE_ITEM LSDTT3_CPU_SOURCES
    "${LSDTT3_SRC_DIR}/cpu/lsdtt_raster_opencv_cpu.cpp"
    "${LSDTT3_SRC_DIR}/cpu/lsdtt_chi_ksn_inversion_cpu.cpp"
)
if(OpenCV_FOUND)
    list(APPEND LSDTT3_CPU_SOURCES "${LSDTT3_SRC_DIR}/cpu/lsdtt_raster_opencv_cpu.cpp")
endif()
if(Eigen3_FOUND)
    list(APPEND LSDTT3_CPU_SOURCES "${LSDTT3_SRC_DIR}/cpu/lsdtt_chi_ksn_inversion_cpu.cpp")
endif()

if(WIN32)
    # Keep the LSDTT implementation inside the extension on Windows. The
    # project does not annotate its C++ API for DLL export, so a separate
    # runtime DLL would not produce the import library required by _core.
    set(LSDTT_PYTHON_RUNTIME_LIBRARY_TYPE STATIC)
else()
    set(LSDTT_PYTHON_RUNTIME_LIBRARY_TYPE SHARED)
endif()
add_library(lsdtt3_runtime ${LSDTT_PYTHON_RUNTIME_LIBRARY_TYPE}
    ${LSDTT3_CORE_SOURCES}
    ${LSDTT3_UTILS_SOURCES}
    ${LSDTT3_CPU_SOURCES}
)
add_library(lsdtt3::runtime ALIAS lsdtt3_runtime)

target_include_directories(lsdtt3_runtime
    PUBLIC
        "${LSDTT3_SRC_DIR}"
        "${LSDTT3_SRC_DIR}/core"
        "${LSDTT3_SRC_DIR}/cpu"
        "${LSDTT3_SRC_DIR}/utils"
        ${GDAL_INCLUDE_DIRS}
)
target_link_libraries(lsdtt3_runtime PUBLIC ${GDAL_LIBRARIES})
target_compile_options(lsdtt3_runtime PRIVATE
    $<$<CXX_COMPILER_ID:GNU,Clang,AppleClang>:-Wall;-Wextra;-Wpedantic;-Wno-unused-parameter>
)

if(OpenMP_CXX_FOUND)
    target_link_libraries(lsdtt3_runtime PUBLIC OpenMP::OpenMP_CXX)
    target_compile_definitions(lsdtt3_runtime PUBLIC LSDTT_HAS_OPENMP=1)
endif()

if(LSDTT_PYTHON_FFTW_TARGET)
    target_link_libraries(lsdtt3_runtime PUBLIC ${LSDTT_PYTHON_FFTW_TARGET})
    target_compile_definitions(lsdtt3_runtime PUBLIC LSDTT_HAS_FFTW=1)
endif()

if(NETCDF_INCLUDE_DIR AND NETCDF_LIBRARY)
    target_include_directories(lsdtt3_runtime PUBLIC ${NETCDF_INCLUDE_DIR})
    target_link_libraries(lsdtt3_runtime PUBLIC ${NETCDF_LIBRARY})
    target_compile_definitions(lsdtt3_runtime PUBLIC LSDTT_HAS_NETCDF=1)
endif()

if(OpenCV_FOUND)
    target_link_libraries(lsdtt3_runtime PUBLIC ${OpenCV_LIBS})
    target_compile_definitions(lsdtt3_runtime PUBLIC LSDTT_HAS_OPENCV=1)
endif()

if(Eigen3_FOUND)
    target_link_libraries(lsdtt3_runtime PUBLIC Eigen3::Eigen)
    target_compile_definitions(lsdtt3_runtime PUBLIC LSDTT_HAS_EIGEN=1)
endif()

# Both shared objects are installed into the Python package. The platform
# runtime path makes _core locate liblsdtt3_runtime beside it, independently
# of the root build.
if(APPLE)
    set(LSDTT_PYTHON_RUNTIME_PATH "@loader_path")
else()
    set(LSDTT_PYTHON_RUNTIME_PATH "$ORIGIN")
endif()

add_subdirectory(bindings)

set_target_properties(lsdtt3_runtime PROPERTIES
    BUILD_RPATH "${LSDTT_PYTHON_RUNTIME_PATH}"
    INSTALL_RPATH "${LSDTT_PYTHON_RUNTIME_PATH}"
)
install(TARGETS lsdtt3_runtime LIBRARY DESTINATION lsdtt3)

# A wheel can carry native command-line programs just like it carries _core.
# The drivers are compiled unchanged from ../cli and link the same runtime as
# the bindings, so the pip and standalone interfaces remain behaviourally
# identical.
if(LSDTT_PYTHON_CLI)
    if(APPLE)
        set(LSDTT_PYTHON_CLI_RUNTIME_PATH "@loader_path/..")
    elseif(WIN32)
        set(LSDTT_PYTHON_CLI_RUNTIME_PATH "")
    else()
        set(LSDTT_PYTHON_CLI_RUNTIME_PATH "$ORIGIN/..")
    endif()

    function(lsdtt_python_cli TARGET SOURCE)
        add_executable(${TARGET} "${LSDTT3_SOURCE_DIR}/${SOURCE}")
        target_include_directories(${TARGET} PRIVATE
            "${LSDTT3_SRC_DIR}"
            "${LSDTT3_SRC_DIR}/core"
            "${LSDTT3_SRC_DIR}/cpu"
            "${LSDTT3_SRC_DIR}/utils"
            ${GDAL_INCLUDE_DIRS}
        )
        target_link_libraries(${TARGET} PRIVATE lsdtt3::runtime)
        target_compile_features(${TARGET} PRIVATE cxx_std_17)
        if(LSDTT_PYTHON_CLI_RUNTIME_PATH)
            set_target_properties(${TARGET} PROPERTIES
                BUILD_RPATH "${LSDTT_PYTHON_CLI_RUNTIME_PATH}"
                INSTALL_RPATH "${LSDTT_PYTHON_CLI_RUNTIME_PATH}"
            )
        endif()
        install(TARGETS ${TARGET} RUNTIME DESTINATION lsdtt3/bin)
    endfunction()

    lsdtt_python_cli(lsdtt-raster-preprocessing "cli/lsdtt-raster-preprocessing.cpp")
    lsdtt_python_cli(lsdtt-synthetic-surfaces "cli/lsdtt-synthetic-surfaces.cpp")
    lsdtt_python_cli(lsdtt-synthetic-networks "cli/lsdtt-synthetic-networks.cpp")
    lsdtt_python_cli(lsdtt-surface-metrics "cli/lsdtt-surface-metrics.cpp")
    lsdtt_python_cli(lsdtt-channel-extraction "cli/lsdtt-channel-extraction.cpp")
    lsdtt_python_cli(lsdtt-internal-drainage "cli/lsdtt-internal-drainage.cpp")
    lsdtt_python_cli(lsdtt-hillslope-channel-coupling "cli/lsdtt-hillslope-channel-coupling.cpp")
    lsdtt_python_cli(lsdtt-drainage-analysis "cli/lsdtt-drainage-analysis.cpp")
    lsdtt_python_cli(lsdtt-feature-finder "cli/lsdtt-feature-finder.cpp")
    lsdtt_python_cli(lsdtt-mountain-front "cli/lsdtt-mountain-front.cpp")
    lsdtt_python_cli(lsdtt-basin-stats "cli/lsdtt-basin-stats.cpp")
    lsdtt_python_cli(lsdtt-basin-extraction "cli/lsdtt-basin-extraction.cpp")
    lsdtt_python_cli(lsdtt-flat-surface-finder "cli/lsdtt-flat-surface-finder.cpp")
    lsdtt_python_cli(lsdtt-fluvial-surfaces "cli/lsdtt-fluvial-surfaces.cpp")
    lsdtt_python_cli(lsdtt-valley-centreline "cli/lsdtt-valley-centreline.cpp")
    lsdtt_python_cli(lsdtt-chi-analysis "cli/lsdtt-chi-analysis.cpp")
    lsdtt_python_cli(lsdtt-graphflood "cli/lsdtt-graphflood.cpp")
    lsdtt_python_cli(lsdtt-channel-network-analysis "cli/lsdtt-channel-network-analysis.cpp")
    lsdtt_python_cli(lsdtt-cross-sections "cli/lsdtt-cross-sections.cpp")
    lsdtt_python_cli(lsdtt-stream-burning "cli/lsdtt-stream-burning.cpp")
    lsdtt_python_cli(lsdtt-water-line-repair "cli/lsdtt-water-line-repair.cpp")
    lsdtt_python_cli(lsdtt-water-centreline "cli/lsdtt-water-centreline.cpp")
    lsdtt_python_cli(lsdtt-cosmo-tool "cli/lsdtt-cosmo-tool.cpp")
    lsdtt_python_cli(lsdtt-cosmo-prediction "cli/lsdtt-cosmo-prediction.cpp")
    lsdtt_python_cli(lsdtt-cosmo-landslides "cli/lsdtt-cosmo-landslides.cpp")
    lsdtt_python_cli(lsdtt-cosmo-column-tool "cli/lsdtt-cosmo-column-tool.cpp")
    lsdtt_python_cli(lsdtt-marsh-topography "cli/lsdtt-marsh-topography.cpp")
    lsdtt_python_cli(lsdtt-linear-stats "cli/lsdtt-linear-stats.cpp")
    lsdtt_python_cli(lsdtt-bearing-comparison "cli/lsdtt-bearing-comparison.cpp")
    lsdtt_python_cli(lsdtt-swath-profile "cli/lsdtt-swath-profile.cpp")
    lsdtt_python_cli(lsdtt-lem "cli/lsdtt-lem.cpp")
    lsdtt_python_cli(lsdtt-sediment-routing "cli/lsdtt-sediment-routing.cpp")
    lsdtt_python_cli(lsdtt-ridgeline-analysis "cli/lsdtt-ridgeline-analysis.cpp")
    lsdtt_python_cli(lsdtt-basin-surface-metrics "cli/lsdtt-basin-surface-metrics.cpp")
    lsdtt_python_cli(lsdtt-divide-migration "cli/lsdtt-divide-migration.cpp")
    lsdtt_python_cli(lsdtt-raster-calculator "cli/lsdtt-raster-calculator.cpp")
    lsdtt_python_cli(lsdtt-fetch-raster "cli/lsdtt-fetch-raster.cpp")
    if(LSDTT_PYTHON_FFTW_TARGET)
        lsdtt_python_cli(lsdtt-spectral-analysis "cli/lsdtt-spectral-analysis.cpp")
    endif()
endif()

install(
    DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/src/lsdtt3/"
    DESTINATION lsdtt3
    FILES_MATCHING PATTERN "*.py"
)

# Vanilla cosmogenic drivers locate this relative to their executable as
# ../share/lsdtt3/crn_parameters.nc.  Keeping that established layout avoids
# a Python-only path convention and leaves their CLI behaviour unchanged.
if(LSDTT_PYTHON_CLI)
    install(FILES "${LSDTT3_SOURCE_DIR}/data/crn_parameters.nc"
        DESTINATION lsdtt3/share/lsdtt3)
endif()

if(LSDTT_PYTHON_RELEASE_WHEEL)
    if(NOT LSDTT_PYTHON_GDAL_DATA_DIR)
        find_program(LSDTT_PYTHON_GDAL_CONFIG_EXECUTABLE NAMES gdal-config)
        if(LSDTT_PYTHON_GDAL_CONFIG_EXECUTABLE)
            execute_process(
                COMMAND "${LSDTT_PYTHON_GDAL_CONFIG_EXECUTABLE}" --datadir
                RESULT_VARIABLE LSDTT_PYTHON_GDAL_CONFIG_RESULT
                OUTPUT_VARIABLE LSDTT_PYTHON_GDAL_DATA_DIR
                OUTPUT_STRIP_TRAILING_WHITESPACE
            )
            if(NOT LSDTT_PYTHON_GDAL_CONFIG_RESULT EQUAL 0)
                unset(LSDTT_PYTHON_GDAL_DATA_DIR)
            endif()
        endif()
    endif()
    if(NOT LSDTT_PYTHON_GDAL_DATA_DIR)
        unset(LSDTT_PYTHON_GDAL_DATA_DIR CACHE)
        find_path(LSDTT_PYTHON_GDAL_DATA_DIR
            NAMES gdalvrt.xsd
            HINTS
                "$ENV{CONDA_PREFIX}"
                "$ENV{VCPKG_INSTALLED_DIR}/$ENV{VCPKG_TARGET_TRIPLET}"
                "$ENV{VCPKG_ROOT}/installed/$ENV{VCPKG_DEFAULT_TRIPLET}"
                "$ENV{VCPKG_ROOT}/installed/x64-windows"
            PATH_SUFFIXES share/gdal
        )
    endif()

    if(NOT LSDTT_PYTHON_PROJ_DATA_DIR)
        unset(LSDTT_PYTHON_PROJ_DATA_DIR CACHE)
        find_path(LSDTT_PYTHON_PROJ_DATA_DIR
            NAMES proj.db
            HINTS
                "$ENV{CONDA_PREFIX}"
                "$ENV{VCPKG_INSTALLED_DIR}/$ENV{VCPKG_TARGET_TRIPLET}"
                "$ENV{VCPKG_ROOT}/installed/$ENV{VCPKG_DEFAULT_TRIPLET}"
                "$ENV{VCPKG_ROOT}/installed/x64-windows"
            PATH_SUFFIXES share/proj
        )
    endif()

    foreach(LSDTT_PYTHON_RUNTIME_DATA
            LSDTT_PYTHON_GDAL_DATA_DIR LSDTT_PYTHON_PROJ_DATA_DIR)
        if(NOT IS_DIRECTORY "${${LSDTT_PYTHON_RUNTIME_DATA}}")
            message(FATAL_ERROR
                "${LSDTT_PYTHON_RUNTIME_DATA} is required when "
                "LSDTT_PYTHON_RELEASE_WHEEL=ON. Set it explicitly to the "
                "corresponding runtime-data directory.")
        endif()
    endforeach()

    install(DIRECTORY "${LSDTT_PYTHON_GDAL_DATA_DIR}/"
        DESTINATION lsdtt3/gdal_data)
    install(DIRECTORY "${LSDTT_PYTHON_PROJ_DATA_DIR}/"
        DESTINATION lsdtt3/proj_data)
endif()
