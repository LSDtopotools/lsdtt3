# -----------------------------------------------------------------------------
# Test-dependency options
#
# The native C++ tests are always registered when BUILD_TESTING is ON: they need
# nothing beyond the C++ build. The Python-driven tests are opt-in, so a machine
# with only the documented C++ packages can configure, build and run the suite:
#
#   LSDTT_ENABLE_PYTHON_TESTS      the CLI harness and the other Python checks
#                                  (standard library only; about a quarter of
#                                  them also need GDAL's command-line programs)
#   LSDTT_ENABLE_PYTHON_GEO_TESTS  the few Python checks that import third-party
#                                  geospatial packages (GDAL bindings, NumPy,
#                                  Rasterio, GeoPandas, Shapely)
#
# Python is looked for ONLY when LSDTT_ENABLE_PYTHON_TESTS is ON.
# -----------------------------------------------------------------------------
option(LSDTT_ENABLE_PYTHON_TESTS "Register Python-driven CTest tests" OFF)
option(LSDTT_ENABLE_PYTHON_GEO_TESTS "Register Python tests needing geospatial packages" OFF)

# Tests that download data: network_*, and the cosmo_realworld_* tier (each
# fetches a Copernicus DEM and runs a full basin inversion, one minute to half
# an hour apiece). Declared here, before its first use below.
option(LSDTT_ENABLE_NETWORK_TESTS "Enable tests that download data" OFF)

if(LSDTT_ENABLE_PYTHON_GEO_TESTS AND NOT LSDTT_ENABLE_PYTHON_TESTS)
    message(FATAL_ERROR
        "LSDTT_ENABLE_PYTHON_GEO_TESTS=ON needs LSDTT_ENABLE_PYTHON_TESTS=ON as "
        "well. Configure with both:\n"
        "  -DLSDTT_ENABLE_PYTHON_TESTS=ON -DLSDTT_ENABLE_PYTHON_GEO_TESTS=ON")
endif()

# 3.10 because the test scripts use `X | None` annotations without
# `from __future__ import annotations`, and `match` statements: an older
# interpreter configures fine and then fails every Python test at import.
if(LSDTT_ENABLE_PYTHON_TESTS)
    find_package(Python3 3.10 COMPONENTS Interpreter REQUIRED)
    message(STATUS "Test Python   : ${Python3_EXECUTABLE} (${Python3_VERSION})")

    # GDAL's command-line programs, which about a quarter of the Python tests
    # run (label needs_gdal_tools). Looked for on PATH, and then beside the
    # interpreter, which is where a Pixi or conda environment keeps them.
    get_filename_component(_lsdtt_py_bindir "${Python3_EXECUTABLE}" DIRECTORY)
    set(_lsdtt_missing_tools "")
    set(_lsdtt_tools_beside_python FALSE)
    # Overridable ONLY so a configure-contract test can make them "missing".
    if(NOT DEFINED LSDTT_TEST_GDAL_PROGRAMS)
        set(LSDTT_TEST_GDAL_PROGRAMS gdalinfo ogrinfo gdal_translate)
    endif()
    foreach(_tool IN LISTS LSDTT_TEST_GDAL_PROGRAMS)
        unset(_lsdtt_tool_path CACHE)
        find_program(_lsdtt_tool_path ${_tool})
        if(NOT _lsdtt_tool_path)
            unset(_lsdtt_tool_path CACHE)
            find_program(_lsdtt_tool_path ${_tool} HINTS "${_lsdtt_py_bindir}" NO_DEFAULT_PATH)
            if(_lsdtt_tool_path)
                set(_lsdtt_tools_beside_python TRUE)
            else()
                list(APPEND _lsdtt_missing_tools ${_tool})
            endif()
        endif()
        unset(_lsdtt_tool_path CACHE)
    endforeach()
    if(_lsdtt_missing_tools)
        string(REPLACE ";" ", " _lsdtt_missing_tools "${_lsdtt_missing_tools}")
        message(FATAL_ERROR
            "LSDTT_ENABLE_PYTHON_TESTS=ON, but GDAL's command-line programs are "
            "missing (${_lsdtt_missing_tools}); about a quarter of the Python tests "
            "run them. Either install them (Ubuntu: sudo apt install gdal-bin), or "
            "use the Pixi test environment, which provides them:\n"
            "  pixi install -e test\n"
            "  cmake ... -DPython3_EXECUTABLE=$(pixi run -e test python -c 'import sys; print(sys.executable)')\n"
            "  pixi run -e test ctest --test-dir <build> ...")
    endif()
    if(_lsdtt_tools_beside_python)
        message(STATUS "Test Python   : GDAL programs found beside the interpreter, "
                       "not on PATH -- run ctest inside that environment "
                       "(e.g. pixi run -e test ctest ...)")
    endif()
endif()

# The geospatial Python tests' imports, checked with the interpreter CMake
# selected, so a missing package is a configure error naming it rather than a
# ModuleNotFoundError halfway through a run.
if(LSDTT_ENABLE_PYTHON_GEO_TESTS)
    set(_lsdtt_missing_modules "")
    foreach(_mod osgeo.gdal numpy rasterio geopandas shapely pytest)
        execute_process(COMMAND "${Python3_EXECUTABLE}" -c "import ${_mod}"
                        RESULT_VARIABLE _lsdtt_rc OUTPUT_QUIET ERROR_QUIET)
        if(NOT _lsdtt_rc EQUAL 0)
            list(APPEND _lsdtt_missing_modules ${_mod})
        endif()
    endforeach()
    if(_lsdtt_missing_modules)
        string(REPLACE ";" ", " _lsdtt_missing_modules "${_lsdtt_missing_modules}")
        message(FATAL_ERROR
            "LSDTT_ENABLE_PYTHON_GEO_TESTS=ON, but ${Python3_EXECUTABLE} cannot "
            "import: ${_lsdtt_missing_modules}. Use the Pixi test environment:\n"
            "  pixi install -e test\n"
            "  cmake ... -DPython3_EXECUTABLE=$(pixi run -e test python -c 'import sys; print(sys.executable)')\n"
            "  pixi run -e test ctest --test-dir <build> -L '^tier_python_geo$'\n"
            "or configure with -DLSDTT_ENABLE_PYTHON_GEO_TESTS=OFF.")
    endif()
    message(STATUS "Test Python   : geospatial imports OK (osgeo.gdal, numpy, rasterio, geopandas, shapely, pytest)")
endif()

# Which GDAL the C++ build uses, relative to any active Pixi/conda environment.
# Built against the SYSTEM GDAL and tested inside Pixi, a CLI would read Pixi's
# PROJ and GDAL data (GDAL_DATA, PROJ_DATA, GDAL_DRIVER_PATH) for a different
# GDAL/PROJ version; the test launcher (scripts/lsdtt_test_env.py) strips them
# for lsdtt3 children using LSDTT_TEST_GDAL_LIBRARY, set on every test below.
# Building AGAINST a Pixi/conda GDAL is a separate, supported profile -- the
# CLIs then need that environment's libraries and data -- so it is reported,
# not refused.
set(_lsdtt_gdal_in_env FALSE)
if(GDAL_LIBRARY MATCHES "/\\.pixi/envs/")
    set(_lsdtt_gdal_in_env TRUE)
elseif(DEFINED ENV{CONDA_PREFIX} AND NOT "$ENV{CONDA_PREFIX}" STREQUAL "")
    # A literal prefix test: a regex built from a path breaks on '+' or '.'.
    string(FIND "${GDAL_LIBRARY}" "$ENV{CONDA_PREFIX}/" _lsdtt_pos)
    if(_lsdtt_pos EQUAL 0)
        set(_lsdtt_gdal_in_env TRUE)
    endif()
endif()
if(_lsdtt_gdal_in_env)
    message(WARNING
        "The C++ build uses a Pixi/conda GDAL (${GDAL_LIBRARY}). That is the "
        "Pixi-built profile: run the CLIs and the tests inside that environment. "
        "For a build against the system GDAL, configure outside Pixi (no "
        "`pixi run cmake`, no activated shell) and pass the test interpreter "
        "with -DPython3_EXECUTABLE.")
endif()

set(LSDTT_TEST_SCRIPT "${CMAKE_CURRENT_SOURCE_DIR}/scripts/run_cli_smoke.py")

# -----------------------------------------------------------------------------
# Runtime bookkeeping for the dependency tiers (tier_native / tier_python_stdlib
# / tier_python_geo / tier_pixi), assigned at the foot of this file.
#
# CMake cannot read a test's COMMAND back, so every registration that starts
# something other than a native executable marks itself right after add_test:
#
#   lsdtt_mark_test_runtime(python <test>)   command is ${Python3_EXECUTABLE}
#   lsdtt_mark_test_runtime(pixi   <test>)   command is `pixi run ...`
#
# An unmarked test is classified tier_native. A Python test that forgets to mark
# itself is therefore mislabelled rather than rejected; the inventory audit
# is what catches that.
#
# lsdtt_mark_test_needs(gdal_tools <test>) records that a test shells out to the
# GDAL command-line programs (gdalinfo, ogrinfo, gdal_translate). On Ubuntu they
# come from gdal-bin, which libgdal-dev does not pull in.
# -----------------------------------------------------------------------------
function(lsdtt_mark_test_runtime RUNTIME TEST_NAME)
    if(NOT RUNTIME MATCHES "^(python|pixi)$")
        message(FATAL_ERROR "lsdtt_mark_test_runtime: unknown runtime '${RUNTIME}'")
    endif()
    set_property(GLOBAL APPEND PROPERTY LSDTT_RUNTIME_${RUNTIME}_TESTS "${TEST_NAME}")
endfunction()

# A Python test that was NOT registered because its option is OFF. Recorded so
# the configure summary can say how many were left out, and so the name checks
# below can tell "switched off" from "misspelt".
function(lsdtt_mark_test_unregistered RUNTIME TEST_NAME)
    set_property(GLOBAL APPEND PROPERTY LSDTT_UNREGISTERED_${RUNTIME}_TESTS "${TEST_NAME}")
    set_property(GLOBAL APPEND PROPERTY LSDTT_UNREGISTERED_TESTS "${TEST_NAME}")
endfunction()

# lsdtt_set_tests_properties() for tests that may be switched off. A name that is
# registered gets the properties; a name recorded by
# lsdtt_mark_test_unregistered() is skipped; any other name is a typo and stops
# configure, exactly as lsdtt_set_tests_properties() itself would.
function(lsdtt_set_tests_properties)
    # Walk ARGV<n> rather than ARGN: a quoted list value such as LABELS "a;b"
    # stays ONE argument there, where ARGN would split it into two.
    math(EXPR _last "${ARGC} - 1")
    set(_names "")
    set(_props_at -1)
    foreach(_i RANGE 0 ${_last})
        if("${ARGV${_i}}" STREQUAL "PROPERTIES")
            set(_props_at ${_i})
            break()
        endif()
        list(APPEND _names "${ARGV${_i}}")
    endforeach()
    if(_props_at LESS 0)
        message(FATAL_ERROR "lsdtt_set_tests_properties: no PROPERTIES keyword")
    endif()

    get_property(_off GLOBAL PROPERTY LSDTT_UNREGISTERED_TESTS)
    set(_live "")
    foreach(_n IN LISTS _names)
        if(TEST "${_n}")
            list(APPEND _live "${_n}")
        elseif(NOT "${_n}" IN_LIST _off)
            message(FATAL_ERROR "lsdtt_set_tests_properties: no test named '${_n}'")
        endif()
    endforeach()
    if(NOT _live)
        return()
    endif()

    math(EXPR _i "${_props_at} + 1")
    while(_i LESS ARGC)
        math(EXPR _j "${_i} + 1")
        if(NOT _j LESS ARGC)
            message(FATAL_ERROR "lsdtt_set_tests_properties: property '${ARGV${_i}}' has no value")
        endif()
        set_tests_properties(${_live} PROPERTIES "${ARGV${_i}}" "${ARGV${_j}}")
        math(EXPR _i "${_i} + 2")
    endwhile()
endfunction()

function(lsdtt_mark_test_needs NEED TEST_NAME)
    set_property(GLOBAL APPEND PROPERTY LSDTT_NEEDS_${NEED}_TESTS "${TEST_NAME}")
endfunction()

# The run_cli_smoke.py assertions that shell out to GDAL's command-line tools.
set(LSDTT_HARNESS_GDAL_TOOL_ARGS
    --expect-valid-percent-100 --expect-valid-percent --expect-raster-size
    --expect-raster-min --expect-raster-stat --expect-raster-checksum
    --expect-raster-values --expect-vector-feature-min)

# Marks a run_cli_smoke.py registration: Python runtime, plus gdal_tools when
# its arguments use one of the assertions above.
#
# It also records the test's EVIDENCE and LEVEL labels from the same arguments
# (tests/cmake/evidence_labels.cmake writes them, and its hand table adds to
# them). Every harness run is a system test. It is `negative` when it expects
# the CLI to fail, and otherwise `smoke` (exit 0, the expected files, the run
# record), plus:
#   characterisation  a stored GDAL checksum or sha256 of an output
#   property          an output with no nodata or NaN: --expect-valid-percent-100,
#                     --expect-raster-min (which also requires 100 % valid), or
#                     --expect-file-not-contains <file>:nan / :-9999
function(lsdtt_mark_harness_test TEST_NAME)
    lsdtt_mark_test_runtime(python ${TEST_NAME})
    foreach(_arg IN LISTS ARGN)
        if("${_arg}" IN_LIST LSDTT_HARNESS_GDAL_TOOL_ARGS)
            lsdtt_mark_test_needs(gdal_tools ${TEST_NAME})
            break()
        endif()
    endforeach()

    set(_negative FALSE)
    set(_evidence "")
    set(_previous "")
    foreach(_arg IN LISTS ARGN)
        if(_arg STREQUAL "--expect-failure")
            set(_negative TRUE)
        elseif(_previous STREQUAL "--expect-exit-code" AND NOT _arg STREQUAL "0")
            set(_negative TRUE)
        elseif(_arg MATCHES "^--expect-(raster-checksum|sha256)$")
            list(APPEND _evidence characterisation)
        elseif(_arg MATCHES "^--expect-(valid-percent-100|raster-min)$")
            list(APPEND _evidence property)
        elseif(_previous STREQUAL "--expect-file-not-contains" AND _arg MATCHES ":(nan|-9999)$")
            list(APPEND _evidence property)
        endif()
        set(_previous "${_arg}")
    endforeach()
    if(_negative)
        list(APPEND _evidence negative)
    else()
        list(APPEND _evidence smoke)
    endif()
    list(REMOVE_DUPLICATES _evidence)
    lsdtt_record_evidence(${TEST_NAME} system ${_evidence})
endfunction()

# Records evidence and level labels for a test at registration; see
# tests/cmake/evidence_labels.cmake, which is the only writer of those labels.
function(lsdtt_record_evidence TEST_NAME LEVEL)
    set_property(GLOBAL APPEND PROPERTY LSDTT_EVIDENCE_${TEST_NAME} ${ARGN})
    set_property(GLOBAL PROPERTY LSDTT_LEVEL_${TEST_NAME} ${LEVEL})
endfunction()

# A run_cli_smoke.py test registered directly rather than through
# add_lsdtt_cli_smoke_test: ARGN is the harness's arguments, exactly as they
# would follow "${LSDTT_TEST_SCRIPT}" in add_test. Marks it like the helper.
function(lsdtt_add_harness_test TEST_NAME)
    add_test(NAME ${TEST_NAME}
             COMMAND ${Python3_EXECUTABLE} "${LSDTT_TEST_SCRIPT}" ${ARGN})
    lsdtt_mark_harness_test(${TEST_NAME} ${ARGN})
endfunction()
set(LSDTT_TEST_DATA_DIR "${CMAKE_CURRENT_SOURCE_DIR}/data")
set(LSDTT_TEST_PARAM_DIR "${CMAKE_CURRENT_SOURCE_DIR}/params")

# A real proj.db for unit_gdal_init_pinned to stage into a fake bundle. Taken
# from the environment when the build host exports it, otherwise searched in the
# usual prefixes. The test skips (exit 77) when this stays empty, so a host
# without PROJ data does not fail the suite.
if(DEFINED ENV{PROJ_DATA} AND EXISTS "$ENV{PROJ_DATA}/proj.db")
    set(LSDTT_TEST_PROJ_DB "$ENV{PROJ_DATA}/proj.db")
elseif(DEFINED ENV{PROJ_LIB} AND EXISTS "$ENV{PROJ_LIB}/proj.db")
    set(LSDTT_TEST_PROJ_DB "$ENV{PROJ_LIB}/proj.db")
else()
    find_file(LSDTT_TEST_PROJ_DB_FOUND proj.db
        PATHS /usr/share/proj /usr/local/share/proj /opt/homebrew/share/proj
              "${CMAKE_PREFIX_PATH}/share/proj"
        NO_DEFAULT_PATH)
    if(LSDTT_TEST_PROJ_DB_FOUND)
        set(LSDTT_TEST_PROJ_DB "${LSDTT_TEST_PROJ_DB_FOUND}")
    else()
        set(LSDTT_TEST_PROJ_DB "")
    endif()
endif()
message(STATUS "Test proj.db  : ${LSDTT_TEST_PROJ_DB}")

# Written by lsdtt-test-divide-migration (the unit_divide_migration CTest
# fixture) and read by the divide-migration CLI tests.
set(LSDTT_DIVIDE_MIGRATION_FIXTURE_DIR
    "${CMAKE_CURRENT_BINARY_DIR}/divide_migration_fixture")

# Written by lsdtt-test-paired-flowpaths (the unit_paired_flowpaths CTest
# fixture): a projected 120x120 diamond-square DEM.  tests/data/tiny30_dem.tif
# is geographic and lsdtt-ridgeline-analysis refuses it outright, so the paired
# flow-path CLI tests need a projected surface of their own.
set(LSDTT_PAIRED_FLOWPATHS_FIXTURE_DIR
    "${CMAKE_CURRENT_BINARY_DIR}/paired_flowpaths_fixture")

# Written by lsdtt-test-synthetic-networks (the unit_synthetic_networks CTest
# fixture): the 3000 m square polygon in EPSG:32758 that the polygon/true_area
# CLI test rasterises.
set(LSDTT_SYNTHETIC_NETWORKS_FIXTURE_DIR
    "${CMAKE_CURRENT_BINARY_DIR}/synthetic_networks_output")

set(LSDTT_INTERNAL_DRAINAGE_FIXTURE_DIR
    "${CMAKE_CURRENT_BINARY_DIR}/internal_drainage_fixture")

set(LSDTT_HILLTOP_COUPLING_FIXTURE_DIR
    "${CMAKE_CURRENT_BINARY_DIR}/hilltop_coupling_fixture")

# Written by lsdtt-test-fluvial-surfaces (the fixture_fluvial_surfaces_cli
# CTest fixture): a projected 200x200 linear-range DEM with flat piedmont and
# a rough mountain zone, plus the channel source point that
# lsdtt-valley-centreline requires.
set(LSDTT_STREAM_BURNING_FIXTURE_DIR
    "${CMAKE_CURRENT_BINARY_DIR}/stream_burning_fixture")

set(LSDTT_BASIN_STATS_FIXTURE_DIR
    "${CMAKE_CURRENT_BINARY_DIR}/basin_stats_fixture")

# Terrace-valley DEM (lsdtt-synthetic-surfaces) plus the source-point and reach
# CSVs, read by the lsdtt-cross-sections and lsdtt-bankfull-width tests.
set(LSDTT_CROSS_SECTIONS_FIXTURE_DIR
    "${CMAKE_CURRENT_BINARY_DIR}/cross_sections_fixture")

# Synthetic marsh DEM written by tests/scripts/make_marsh_topography_fixture.py.
set(LSDTT_MARSH_TOPOGRAPHY_FIXTURE_DIR
    "${CMAKE_CURRENT_BINARY_DIR}/marsh_topography_fixture")

# Reference and shifted DEM pair written by
# tests/scripts/make_dem_coregistration_fixture.py.
set(LSDTT_DEM_COREGISTRATION_FIXTURE_DIR
    "${CMAKE_CURRENT_BINARY_DIR}/dem_coregistration_fixture")

# Written by lsdtt-test-hype-cli-fixture (the fixture_hype_cli CTest fixture):
# a 60x60 projected Y-network DEM, land-use and soil rasters, the four HYPE
# parameter tables (plus three broken ones), 120 days of forcing and an outlet
# CSV, read by the column_model: hype tests of lsdtt-distributed-hydrology.
set(LSDTT_HYPE_CLI_FIXTURE_DIR
    "${CMAKE_CURRENT_BINARY_DIR}/hype_cli_fixture")

set(LSDTT_FLUVIAL_SURFACES_FIXTURE_DIR
    "${CMAKE_CURRENT_BINARY_DIR}/fluvial_surfaces_fixture")

function(add_lsdtt_cli_smoke_test TEST_NAME TARGET_NAME PARAM_TEMPLATE EXPECTED_OUTPUT)
    if(NOT LSDTT_ENABLE_PYTHON_TESTS)
        lsdtt_mark_test_unregistered(python ${TEST_NAME})
        return()
    endif()
    set(work_dir "${CMAKE_CURRENT_BINARY_DIR}/${TEST_NAME}")

    add_test(
        NAME ${TEST_NAME}
        COMMAND
            ${Python3_EXECUTABLE} "${LSDTT_TEST_SCRIPT}"
            --exe "$<TARGET_FILE:${TARGET_NAME}>"
            --cli-src-dir "${CMAKE_SOURCE_DIR}/cli"
            --param-template "${LSDTT_TEST_PARAM_DIR}/${PARAM_TEMPLATE}"
            --read-path "${LSDTT_TEST_DATA_DIR}"
            --work-dir "${work_dir}"
            --expect "${EXPECTED_OUTPUT}"
            ${ARGN}
    )
    lsdtt_mark_harness_test(${TEST_NAME} ${ARGN})

    # Evidence and level labels (smoke, negative, ...) are written by
    # tests/cmake/evidence_labels.cmake, never here.
    lsdtt_set_tests_properties(${TEST_NAME} PROPERTIES LABELS "cli")
endfunction()

# =============================================================================
# lsdtt-synthetic-networks (OCN generator)
#
# The engine properties (forest invariants, energy, D8 round-trip fidelity,
# cross-tool agreement with the LEM steady state) live in the C++ driver
# lsdtt-test-synthetic-networks; these CLI tests carry only what is genuinely
# CLI: parameter plumbing, output filenames, the failure paths, the polygon /
# CRS path, and the characterisation ORACLES on the shipped artefacts.
#
# Artefact classification (docs/developer_notes/testing_and_benchmarking.md): the GDAL checksum
# of *_ocn_drainage_cells.tif and the sha256 of *_ocn_summary.csv are ORACLES
# (refactor tripwires, re-freezable by a deliberate act). The --expect-raster-stat
# bands on *_ocn_dem.tif are BANDS, not oracles, because float chi accumulation
# is pow-implementation sensitive. There is NO numerics sentinel for this module
# and there must not be one: the OCN has no validated numeric standard.
#
# Every param file pins n_iter explicitly. Leaving it at the 40*n_rows*n_cols
# default would couple every pinned artefact to a tuning knob.
# =============================================================================
add_lsdtt_cli_smoke_test(
    smoke_synthetic_networks_rect
    lsdtt-synthetic-networks
    synthetic_networks_smoke.param.in
    ocnsmoke_ocn_dem.tif
    --expect ocnsmoke_ocn_drainage_cells.tif
    --expect ocnsmoke_ocn_summary.csv
    --expect-raster-size ocnsmoke_ocn_dem.tif:32x32
    --expect-valid-percent-100 ocnsmoke_ocn_dem.tif
    --expect-file-contains "ocnsmoke_ocn_summary.csv:n_rows,n_cols,cell_size,gamma,n_iter,cooling_rate,temperature_scale,initial_no_cooling_fraction,seed,initial_energy,final_energy,routing_fidelity"
    --stdout-contains "D8 round-trip fidelity (interior)"
)

# The parameter declarations must all render; -h needs no parameter file.
add_test(NAME smoke_synthetic_networks_help
         COMMAND lsdtt-synthetic-networks -h)
lsdtt_set_tests_properties(smoke_synthetic_networks_help PROPERTIES
    LABELS "cli;synthetic_networks"
    PASS_REGULAR_EXPRESSION "energy_exponent_gamma")

add_lsdtt_cli_smoke_test(
    smoke_synthetic_networks_phase5_defaults
    lsdtt-synthetic-networks
    synthetic_networks_phase5_defaults.param.in
    ocndefaults5_ocn_dem.tif
    --expect ocndefaults5_ocn_summary.csv
    --expect-raster-size ocndefaults5_ocn_dem.tif:128x128
    --expect-file-contains "ocndefaults5_ocn_summary.csv:128,128,1,0.5,1,1,1,0,33,"
)

add_lsdtt_cli_smoke_test(
    smoke_synthetic_networks_phase5_legacy
    lsdtt-synthetic-networks
    synthetic_networks_phase5_legacy.param.in
    ocnlegacy5_ocn_dem.tif
    --expect ocnlegacy5_ocn_summary.csv
    --expect-raster-size ocnlegacy5_ocn_dem.tif:8x8
    --expect-file-contains "ocnlegacy5_ocn_summary.csv:8,8,99,"
    --expect-file-contains "ocnlegacy5_ingested_params.param:cell_size_metres: 99"
    --expect-file-contains "ocnlegacy5_ingested_params.param:coastline_blur_passes: 2"
    --expect-file-contains "ocnlegacy5_ingested_params.param:incision_threshold_area_m2: 20"
    --expect-file-contains "ocnlegacy5_ingested_params.param:polygon_fname: unused_legacy_polygon.geojson"
    --expect-file-contains "ocnlegacy5_ingested_params.param:steady_state_tolerance_dz_metres: 0.02"
    --stderr-contains "'cell_size' is deprecated"
    --stderr-contains "'noise_blur_passes' is deprecated"
    --stderr-contains "'incision_threshold_area' is deprecated"
)

add_lsdtt_cli_smoke_test(
    smoke_synthetic_networks_phase5_alias_conflict
    lsdtt-synthetic-networks
    synthetic_networks_phase5_alias_conflict.param.in
    ocnconflict5_ocn_dem.tif
    --expect-file-contains "ocnconflict5_ocn_summary.csv:8,8,10,"
    --stderr-contains "'cell_size_metres' was also supplied"
)

# noise_island + outlet_mode auto. The valid percent is the assertion that the
# MASK actually reached the raster: 890 land pixels of 1600 = 55.62 %. A run
# that silently ignored domain_type would come back 100 % valid.
add_lsdtt_cli_smoke_test(
    smoke_synthetic_networks_island
    lsdtt-synthetic-networks
    synthetic_networks_island.param.in
    ocnisl_ocn_dem.tif
    --expect ocnisl_ocn_drainage_cells.tif
    --expect ocnisl_ocn_summary.csv
    --expect-raster-size ocnisl_ocn_dem.tif:40x40
    --expect-valid-percent ocnisl_ocn_dem.tif:55.62
    --stdout-contains "Domain: noise_island - 890 land pixels of 1600"
)

# polygon + true_area: the grid is DERIVED from the polygon extent and
# cell_size, not taken from n_rows/n_cols, and the polygon's CRS is inherited.
# The fixture is the 3000 m square in EPSG:32758 written by
# lsdtt-test-synthetic-networks, so the expected 14x14 is exact arithmetic
# (ceil(3000/250) = 12 plus a one-cell margin each side) rather than a band.
if(LSDTT_ENABLE_PYTHON_TESTS)
    lsdtt_add_harness_test(smoke_synthetic_networks_polygon_area
            --exe "$<TARGET_FILE:lsdtt-synthetic-networks>"
            --param-template "${LSDTT_TEST_PARAM_DIR}/synthetic_networks_polygon.param.in"
            --read-path "${LSDTT_SYNTHETIC_NETWORKS_FIXTURE_DIR}"
            --work-dir "${CMAKE_CURRENT_BINARY_DIR}/smoke_synthetic_networks_polygon_area"
            --expect ocnpoly_ocn_dem.tif
            --expect ocnpoly_ocn_drainage_cells.tif
        --expect-raster-size ocnpoly_ocn_dem.tif:14x14
        --stdout-contains "Polygon true-area grid: 14 x 14 @ 250 m"
        --stdout-contains "Domain: polygon - 144 land pixels of 196"
        --stderr-contains "'polygon_file' is deprecated"
    )
else()
    lsdtt_mark_test_unregistered(python smoke_synthetic_networks_polygon_area)
endif()
lsdtt_set_tests_properties(smoke_synthetic_networks_polygon_area PROPERTIES
    LABELS "cli;synthetic_networks"
    FIXTURES_REQUIRED synthetic_networks_fixture)

# B5: OCN -> equilibrate in one invocation. Asserts the stage actually ran and
# produced its two extra artefacts, not merely that the CLI exited 0 -- an
# equilibration that silently no-ops would still write the OCN DEM and pass a
# file-existence check.
#
# The cycle number is pinned deliberately. MEASURED 2026-08-14: this fixture
# reaches a period-1 fixed point at cycle 10 (first seen at cycle 9), with
# max_cycles set to 20 so there is headroom to notice it taking LONGER rather
# than the run simply hitting the cap. An earlier draft used max_cycles 8 and
# did not converge at all (max|dz| still 47.7 m), which is why the cap matters:
# a 32x32 OCN at 1000 m relief needs about ten cycles, not a handful.
add_lsdtt_cli_smoke_test(
    smoke_synthetic_networks_equilibrate
    lsdtt-synthetic-networks
    synthetic_networks_equilibrate.param.in
    ocneq_ocn_equilibrated_dem.tif
    --expect ocneq_ocn_dem.tif
    --expect ocneq_ocn_snap_convergence.csv
    --expect-raster-size ocneq_ocn_equilibrated_dem.tif:32x32
    --expect-valid-percent-100 ocneq_ocn_equilibrated_dem.tif
    --stdout-contains "Equilibrating (max 20 cycles, convergence checking on)"
    --stdout-contains "Converged at cycle 5: flow graph is a fixed point"
    --stdout-contains "Equilibrated DEM vs the ORIGINAL OCN graph"
)

# Closing the south edge on a single-south OCN leaves no path to base level. The
# router would catch this on its own (since 2026-08-14 it reports unresolvable
# cells as a normal exception rather than terminating), but only as "could not
# assign a receiver at row R, col C" -- which does not tell the user that the
# cause is a parameter they set. Pinned so the up-front, parameter-naming check
# is not quietly dropped in favour of the generic message.
add_lsdtt_cli_smoke_test(
    failure_synthetic_networks_closed_outlet
    lsdtt-synthetic-networks
    synthetic_networks_closed_outlet.param.in
    ocnclosed_ocn_dem.tif
    --expect-failure
    --stderr-contains "FATAL: boundary_conditions 'nnnn' closes the south edge"
)

# An unknown domain_type is DOCUMENTED to warn and fall back to the rectangle.
# Pinning it makes a future decision to reject it a deliberate one.
add_lsdtt_cli_smoke_test(
    smoke_synthetic_networks_unknown_domain
    lsdtt-synthetic-networks
    synthetic_networks_unknown_domain.param.in
    ocnbanana_ocn_dem.tif
    --expect-raster-size ocnbanana_ocn_dem.tif:16x16
    --expect-valid-percent-100 ocnbanana_ocn_dem.tif
    --stderr-contains "Warning: unknown domain_type 'banana' - using rectangle."
)

# ORACLE for the DEFAULT (length-weighted) energy, measured 2026-08-14.
#
# E = sum(A^gamma * L) is the default from 2026-08-14 (Strong & Mudd 2022 PNAS
# eq. 4; Goren & Shelef 2024 ESurf 12:1347 eq. 5). It changes the networks
# substantially -- drainage checksum 3275 -> 3905, final energy 2384 -> 295472
# (the units change, since energy now carries a length) -- so it needs its own
# oracle. The OCNet-parity oracle below is unchanged and still pins the old
# numbers exactly, which is what guarantees the parity switch really is parity.
add_lsdtt_cli_smoke_test(
    rigorous_synthetic_networks_oracle_lengthweighted
    lsdtt-synthetic-networks
    synthetic_networks_lengthweighted.param.in
    ocnlw_ocn_drainage_cells.tif
    --expect ocnlw_ocn_summary.csv
    --expect-raster-checksum ocnlw_ocn_drainage_cells.tif:3905
    --expect-sha256 ocnlw_ocn_summary.csv:9adc345d2046991532144521c61ec19adf7bac09cb7d2a9dfe5bcbb0df43e168
)

# ORACLE (refactor tripwire), measured 2026-08-13 with GDAL from the repo
# toolchain on the same param file as smoke_synthetic_networks_rect. The
# drainage-area field is integer-valued and exactly representable in float32, so
# its GDAL checksum is a complete fingerprint of the whole tree; the summary CSV
# prints energies at 6 significant figures, robust to last-bit drift. Both are
# re-freezable by a deliberate act, NOT numerics sentinels.
#
# THIS ORACLE NOW PINS OCNet PARITY, not the default. Its param file sets
# energy_length_weighted: false. Its values are UNCHANGED across the 2026-08-14
# switch to a length-weighted default -- same checksum 3275, same sha256 -- which
# is the evidence that the parity mode is genuinely bit-for-bit OCNet and not
# merely "close". If this test ever moves, parity has been broken.
#
# RE-FROZEN 2026-08-14, a deliberate act with a stated reason: B5 appended six
# equilibration columns to the summary CSV (equilibrated, eq_cycles_run,
# eq_converged, eq_limit_cycle_period, eq_divide_step_p90_m, eq_ocn_agreement).
# ONLY the CSV hash moved -- the drainage-area raster checksum is unchanged at
# 3275, which is the evidence that the generation path itself was not touched.
# With run_equilibration defaulting to false the new columns read
# "0,0,0,0,-1,-1", so the default output is unchanged in substance.
add_lsdtt_cli_smoke_test(
    rigorous_synthetic_networks_oracle
    lsdtt-synthetic-networks
    synthetic_networks_smoke.param.in
    ocnsmoke_ocn_drainage_cells.tif
    --expect ocnsmoke_ocn_summary.csv
    --expect-raster-checksum ocnsmoke_ocn_drainage_cells.tif:3275
    --expect-sha256 ocnsmoke_ocn_summary.csv:0f9d15a74aef81fc89b2087ddd01bc2be66633bb0ce175b4513f92cd9e03c5fc
    --expect-file-contains "ocnsmoke_ocn_summary.csv:32,32,100,0.5,20000,1,0,0,7,3118.46,2384.26,0.87"
)

# BANDS, not an oracle: float chi accumulation is pow-implementation sensitive,
# so a hash on the DEM would be flaky across libm versions. Two-sided on
# purpose - a flat surface would satisfy a one-sided bound. Measured 2026-08-13:
# min 0, max 1000.0000610352, mean 435.14471304324.
add_lsdtt_cli_smoke_test(
    rigorous_synthetic_networks_dem_band
    lsdtt-synthetic-networks
    synthetic_networks_smoke.param.in
    ocnsmoke_ocn_dem.tif
    --expect-raster-stat ocnsmoke_ocn_dem.tif:min:0.0:0.0
    --expect-raster-stat ocnsmoke_ocn_dem.tif:max:999.9:1000.1
    --expect-raster-stat ocnsmoke_ocn_dem.tif:mean:430.0:441.0
)

# Failure tier - exact messages from the engine guards and the CLI guards.
add_lsdtt_cli_smoke_test(
    failure_synthetic_networks_bad_gamma
    lsdtt-synthetic-networks
    synthetic_networks_bad_gamma.param.in
    ocnbad_gamma_ocn_dem.tif
    --expect-failure
    --stderr-contains "FATAL: network generation failed: lsdtt_synthetic_network_cpu: gamma must be in (0, 1)"
)

add_lsdtt_cli_smoke_test(
    failure_synthetic_networks_bad_gamma_zero
    lsdtt-synthetic-networks
    synthetic_networks_bad_gamma_zero.param.in
    ocnbad_gamma_zero_ocn_dem.tif
    --expect-failure
    --stderr-contains "FATAL: network generation failed: lsdtt_synthetic_network_cpu: gamma must be in (0, 1)"
)

add_lsdtt_cli_smoke_test(
    failure_synthetic_networks_bad_grid
    lsdtt-synthetic-networks
    synthetic_networks_bad_grid.param.in
    ocnbadgrid_ocn_dem.tif
    --expect-failure
    --stderr-contains "FATAL: network generation failed: lsdtt_synthetic_network_cpu: n_rows and n_cols must be >= 2"
)

add_lsdtt_cli_smoke_test(
    failure_synthetic_networks_bad_cell_size
    lsdtt-synthetic-networks
    synthetic_networks_bad_cell_size.param.in
    ocnbadcell_ocn_dem.tif
    --expect-failure
    --stderr-contains "FATAL: network generation failed: lsdtt_synthetic_network_cpu: cell_size must be finite and > 0"
)

add_lsdtt_cli_smoke_test(
    failure_synthetic_networks_polygon_no_file
    lsdtt-synthetic-networks
    synthetic_networks_polygon_no_file.param.in
    ocnpolynofile_ocn_dem.tif
    --expect-failure
    --stderr-contains "FATAL: domain_type=polygon requires polygon_fname."
)

add_lsdtt_cli_smoke_test(
    failure_synthetic_networks_polygon_missing_file
    lsdtt-synthetic-networks
    synthetic_networks_polygon_missing_file.param.in
    ocnpolymissing_ocn_dem.tif
    --expect-failure
    --stderr-contains "FATAL: domain construction failed:"
)

add_lsdtt_cli_smoke_test(
    failure_synthetic_networks_bad_baselevel
    lsdtt-synthetic-networks
    synthetic_networks_bad_baselevel.param.in
    ocnbadbase_ocn_dem.tif
    --expect-failure
    --stderr-contains "FATAL: baselevel_elevation must be > -9000."
)

add_lsdtt_cli_smoke_test(
    failure_synthetic_networks_lat_without_lon
    lsdtt-synthetic-networks
    synthetic_networks_lat_without_lon.param.in
    ocnlatonly_ocn_dem.tif
    --expect-failure
    --stderr-contains "FATAL: centre_latitude and centre_longitude must be set together."
)

add_lsdtt_cli_smoke_test(
    smoke_raster_preprocessing
    lsdtt-raster-preprocessing
    raster_preprocessing_smoke.param.in
    tiny_fill.tif
)

# The legacy replace_pixels_from_csv fixture sets listed cells to new values, and a value equal to
# the raster nodata must also clear the validity mask rather than merely writing
# the number. The fixture covers all three outcomes -- a real replacement, a
# blank-to-nodata, and a point off the grid, which must warn and be skipped
# rather than fail the run.
add_lsdtt_cli_smoke_test(
    smoke_raster_preprocessing_replace_pixels
    lsdtt-raster-preprocessing
    raster_preprocessing_replace_pixels.param.in
    tiny_replaced_pixels_replaced.tif
    --stdout-contains "Replaced 1 pixel(s)\; blanked 1 to nodata."
    --stderr-contains "1 listed pixel(s) fall outside the raster"
    --stderr-contains "'replace_pixels_from_csv' is deprecated"
)

add_lsdtt_cli_smoke_test(
    smoke_raster_preprocessing_tif
    lsdtt-raster-preprocessing
    raster_preprocessing_tif_smoke.param.in
    tiny30_fill.tif
)

# Least-cost breaching (Lindsay 2020, the Whitebox 2.0.6 port) through the CLI.
# The numerics are asserted bitwise against Whitebox in unit_breach_least_cost;
# these pin the parameter plumbing. tiny_breach_lc_dem.asc is a 7x7 bowl, pit
# 1.0 at (3,3) inside a ring of 2s and a rim of 6s, with a 4 saddle at (3,5) and
# a 3 on the east border at (3,6); breach_lc_flat_increment 0.125 keeps every
# value dyadic.
#   outlet:   the pit is raised to 1.875 and breached NE then E off the grid:
#             (2,4) 1.75, (3,5) 1.625, border (3,6) 1.5.
#   whitebox: the edge is not an outlet and no border cell is low enough, so
#             the pit is unsolved and the remainder fill raises the bowl.
# Neither file sets breach_lc_max_dist, so each pins the edge policy's default
# (1000 outlet, 100 whitebox); the max_dist failure test sets it explicitly.
add_lsdtt_cli_smoke_test(
    smoke_raster_preprocessing_breach_least_cost
    lsdtt-raster-preprocessing
    raster_preprocessing_breach_lc_outlet.param.in
    breachlc_outlet_breached.tif
    --stdout-contains "1 pits found, 0 skipped (no longer pits), 1 breached, 0 unsolved"
)

add_lsdtt_cli_smoke_test(
    rigorous_raster_preprocessing_breach_least_cost
    lsdtt-raster-preprocessing
    raster_preprocessing_breach_lc_outlet.param.in
    breachlc_outlet_breached.tif
    --stdout-contains "breach_lc_max_dist = 1000 (default for this edge policy), breach_lc_edge_policy = outlet"
    --stdout-contains "1 pits found, 0 skipped (no longer pits), 1 breached, 0 unsolved. flat_increment = 0.1250, edge_policy = outlet"
    --expect-raster-values "breachlc_outlet_breached.tif:1.5,1.625,1.75,1.875,2,6,9"
    --expect-raster-checksum breachlc_outlet_breached.tif:305
)

add_lsdtt_cli_smoke_test(
    rigorous_raster_preprocessing_breach_lc_whitebox_edge
    lsdtt-raster-preprocessing
    raster_preprocessing_breach_lc_whitebox.param.in
    breachlc_whitebox_breached.tif
    --stdout-contains "breach_lc_max_dist = 100 (default for this edge policy), breach_lc_edge_policy = whitebox"
    --stdout-contains "1 pits found, 0 skipped (no longer pits), 0 breached, 1 unsolved. flat_increment = 0.1250, edge_policy = whitebox"
    --stdout-contains "1 interior cells without descent - filling remainder with Priority-Flood, min_slope = 0.0001, flat resolution = inline_gradient."
    --expect-raster-checksum breachlc_whitebox_breached.tif:326
)

add_lsdtt_cli_smoke_test(
    failure_raster_preprocessing_breach_lc_max_dist
    lsdtt-raster-preprocessing
    raster_preprocessing_breach_lc_bad_max_dist.param.in
    breachlc_bad_breached.tif
    --expect-failure
    --stderr-contains "FATAL: lsdtt::breach_depressions_least_cost: max_dist must be between 1 and 32767 (maximum breach path length in cells)\; got 0"
)

add_lsdtt_cli_smoke_test(
    failure_raster_preprocessing_breach_lc_edge_policy
    lsdtt-raster-preprocessing
    raster_preprocessing_breach_lc_bad_edge_policy.param.in
    breachlc_badedge_breached.tif
    --expect-failure
    --stderr-contains "FATAL: breach_lc_edge_policy must be one of outlet, whitebox\; got 'Outlet'"
)

add_lsdtt_cli_smoke_test(
    smoke_channel_extraction
    lsdtt-channel-extraction
    channel_extraction_smoke.param.in
    tiny_channel_network_lines.fgb
)

add_lsdtt_cli_smoke_test(
    rigorous_channel_lines_feature_count
    lsdtt-channel-extraction
    channel_extraction_smoke.param.in
    tiny_channel_network_lines.fgb
    --expect-vector-feature-min tiny_channel_network_lines.fgb:1
)

add_lsdtt_cli_smoke_test(
    smoke_basin_extraction
    lsdtt-basin-extraction
    basin_extraction_smoke.param.in
    tiny_junction_basin_outlines.fgb
)

# Outlets snapped below threshold_contributing_pixels land off the channel
# network. Each must be skipped with a warning, not abort the run (it aborted
# with an uncaught exception until 2026-10-14).
add_lsdtt_cli_smoke_test(
    smoke_basin_extraction_outlets_off_network
    lsdtt-basin-extraction
    basin_extraction_outlets_off_network.param.in
    bx_junction_basins.tif
    --read-path "${LSDTT_BASIN_STATS_FIXTURE_DIR}/out"
    --stderr-contains "is not on the channel network"
    --stderr-contains "Set basin_outlet_contributing_pixels_threshold to at least threshold_contributing_pixels."
    --stdout-contains "Selected basins: 4"
)
lsdtt_set_tests_properties(smoke_basin_extraction_outlets_off_network PROPERTIES
    FIXTURES_REQUIRED basin_stats_fixture)

add_lsdtt_cli_smoke_test(
    rigorous_basin_extraction_legacy_threshold
    lsdtt-basin-extraction
    basin_extraction_legacy_threshold.param.in
    basin_legacy_junction_basin_outlines.fgb
    --stdout-contains "Contributing-pixel threshold: 1"
    --stderr-contains "'source_pixels_threshold' is deprecated"
    --stderr-contains "use 'threshold_contributing_pixels'"
)

add_lsdtt_cli_smoke_test(
    rigorous_basin_outlines_feature_count
    lsdtt-basin-extraction
    basin_extraction_smoke.param.in
    tiny_junction_basin_outlines.fgb
    --expect-vector-feature-min tiny_junction_basin_outlines.fgb:1
)

add_lsdtt_cli_smoke_test(
    smoke_nodata_fill
    lsdtt-raster-preprocessing
    raster_preprocessing_nodata_smoke.param.in
    tiny_nodata_nodata_filled.tif
    --expect-valid-percent-100 tiny_nodata_nodata_filled.tif
)

# interior_nodata_method (shared by lsdtt-raster-preprocessing and the analysis
# CLIs). tiny_dem_nodata has one interior nodata pixel; tiny_dem_edge_nodata has
# one on the edge, which is never filled. sink is the default and unchanged
# behaviour: the pixel stays nodata and flow ends there.
add_lsdtt_cli_smoke_test(
    smoke_interior_nodata_fill_drainage_analysis
    lsdtt-drainage-analysis
    drainage_analysis_interior_nodata_fill.param.in
    da_interior_fill_d8_drainage_area.tif
    --stdout-contains "Interior nodata: 1 pixel(s) filled by laplacian"
    --expect-valid-percent-100 da_interior_fill_d8_drainage_area.tif
)
add_lsdtt_cli_smoke_test(
    smoke_interior_nodata_sink_drainage_analysis
    lsdtt-drainage-analysis
    drainage_analysis_interior_nodata_sink.param.in
    da_interior_sink_d8_drainage_area.tif
    --stdout-contains "Interior nodata: 1 pixel(s), kept as sinks"
    --expect-valid-percent da_interior_sink_d8_drainage_area.tif:96
)
add_lsdtt_cli_smoke_test(
    smoke_interior_nodata_edge_not_filled_drainage_analysis
    lsdtt-drainage-analysis
    drainage_analysis_edge_nodata_fill.param.in
    da_edge_fill_d8_drainage_area.tif
    --stdout-contains "Interior nodata: none (interior_nodata_method: fill)"
)
add_lsdtt_cli_smoke_test(
    failure_interior_nodata_invalid_drainage_analysis
    lsdtt-drainage-analysis
    drainage_analysis_interior_nodata_invalid.param.in
    da_interior_bad_d8_drainage_area.tif
    --expect-failure
    --stderr-contains "interior_nodata_method must be sink or fill, not 'filled'"
)
# In preprocessing, fill also feeds the later steps: the pit-filled DEM has no
# interior void (write_filled_nodata_raster alone only writes a file).
add_lsdtt_cli_smoke_test(
    smoke_interior_nodata_fill_preprocessing
    lsdtt-raster-preprocessing
    raster_preprocessing_interior_nodata_fill.param.in
    tiny_interior_fill_fill.tif
    --stdout-contains "Interior nodata: 1 pixel(s) filled by laplacian"
    --expect-valid-percent-100 tiny_interior_fill_fill.tif
)
# A negative min_elevation is a threshold, not "disabled" (only -9999 or lower
# disables): on tiny_dem lowered by 160 m (-88 to -60), -80 masks 7 cells.
add_lsdtt_cli_smoke_test(
    smoke_nodata_min_elevation_negative_preprocessing
    lsdtt-raster-preprocessing
    raster_preprocessing_min_elevation_negative.param.in
    tiny_below_sea_fill.tif
    --stdout-contains "Threshold masking changed 7 pixel(s)"
    --expect-raster-stat tiny_below_sea_fill.tif:min:-79:-79
)

add_lsdtt_cli_smoke_test(
    smoke_nodata_edge_unfilled
    lsdtt-raster-preprocessing
    raster_preprocessing_nodata_edge_smoke.param.in
    tiny_edge_nodata_filled.tif
    --expect-valid-percent tiny_edge_nodata_filled.tif:96
)

# Oblique frames in lsdtt-raster-preprocessing: rotation by a line, a polygon
# and the valid-data footprint, the rotate-back mode, and the refusals. The
# frame geometry and the round trip onto the source pixels are asserted in
# unit_oblique_frame_preprocess; these pin the plumbing, file names and frame
# file keys. Fixtures in tests/data/rotation: WGS84 GeoJSON line and polygon
# over flat_resolution_lakeland_crop_dem.tif (UTM 15N, 128 x 128 at 30 m), and
# that DEM rotated by the line (lakeland_line_oblique.tif with its frame file),
# written by this CLI.
add_lsdtt_cli_smoke_test(
    smoke_raster_preprocessing_rotate_line
    lsdtt-raster-preprocessing
    raster_preprocessing_rotate_line.param.in
    rotline_oblique.tif
    --expect rotline_oblique_rotation_frame.txt
    --expect-raster-size rotline_oblique.tif:144x40
    --expect-raster-size rotline_fill.tif:144x40
    --expect-file-contains "rotline_oblique_rotation_frame.txt:frame_source: line"
    --expect-file-contains "rotline_oblique_rotation_frame.txt:feature_fname: 'lakeland_line.geojson'"
    --expect-file-contains "rotline_oblique_rotation_frame.txt:target_crs_epsg: 32615"
    --expect-file-contains "rotline_oblique_rotation_frame.txt:target_origin_x: 564930.0"
    --expect-file-contains "rotline_oblique_rotation_frame.txt:target_origin_y: 5600760.0"
    --expect-file-contains "rotline_oblique_rotation_frame.txt:source_dem_fname: 'flat_resolution_lakeland_crop_dem.tif'"
    --stdout-contains "Size    : 4320 m x 1200 m (144 x 40 cells of 30 m)"
)

add_lsdtt_cli_smoke_test(
    smoke_raster_preprocessing_rotate_polygon
    lsdtt-raster-preprocessing
    raster_preprocessing_rotate_polygon.param.in
    rotpoly_oblique.tif
    --expect-raster-size rotpoly_oblique.tif:122x52
    --expect-file-contains "rotpoly_oblique_rotation_frame.txt:frame_source: polygon"
    --expect-file-contains "rotpoly_oblique_rotation_frame.txt:buffer_all_sides: true"
)

add_lsdtt_cli_smoke_test(
    smoke_raster_preprocessing_rotate_footprint
    lsdtt-raster-preprocessing
    raster_preprocessing_rotate_footprint.param.in
    rotfoot_oblique.tif
    --expect-raster-size rotfoot_oblique.tif:150x150
    --expect-file-contains "rotfoot_oblique_rotation_frame.txt:frame_source: data_footprint"
    --stderr-contains "axis is ill-defined"
)

add_lsdtt_cli_smoke_test(
    smoke_raster_preprocessing_rotate_back
    lsdtt-raster-preprocessing
    raster_preprocessing_rotate_back.param.in
    lakeland_line_oblique_rotated_back.tif
    --expect-raster-size lakeland_line_oblique_rotated_back.tif:136x125
    --stdout-contains "Rotated back (near)"
)

foreach(_rp_case IN ITEMS
        "rotate_both_switches|rotate_to_linear_feature and rotate_to_data_footprint cannot be combined"
        "rotate_force_reproject|cannot be combined with force_reproject"
        "rotate_missing_fname|needs rotation_feature_fname"
        "rotate_back_with_fill|cannot be combined with: write_fill_raster"
        "rotate_back_no_list|needs rotate_back_file_list"
        "rotate_back_not_in_frame|not the frame CRS")
    string(REPLACE "|" ";" _rp_parts "${_rp_case}")
    list(GET _rp_parts 0 _rp_name)
    list(GET _rp_parts 1 _rp_message)
    add_lsdtt_cli_smoke_test(
        failure_raster_preprocessing_${_rp_name}
        lsdtt-raster-preprocessing
        raster_preprocessing_${_rp_name}.param.in
        unused_expected_output
        --expect-failure
        --stderr-contains "${_rp_message}"
    )
endforeach()

add_test(NAME smoke_raster_preprocessing_help_rotation
         COMMAND lsdtt-raster-preprocessing -h rotate_to_linear_feature)
lsdtt_set_tests_properties(smoke_raster_preprocessing_help_rotation PROPERTIES
    LABELS "cli"
    PASS_REGULAR_EXPRESSION "minimum-area")

# bearings_relative_to_true_north. Fixtures in tests/data/true_north: two
# planar DEMs dipping to TRUE azimuth 200 deg at their centre, one in the
# oblique frame of rotation/lakeland_line_oblique.tif (true north 39.87 deg
# clockwise of grid up) and one in UTM 15N (0.73 deg anticlockwise), with the
# grid bearing set from pyproj's meridian convergence; six outlets and the
# lakeland line in the frame CRS for the basin tools. Raster outputs are
# asserted against 200 (aspect) or 110 (contour orientation); table outputs
# pin the plumbing (start-up line, columns), and the per-row conversion is
# asserted in unit_true_north.
set(_tn_on "bearings are converted to true north at each output's location")
set(_tn_frame "Grid north is 39.9 deg anticlockwise of true north at the grid centre")
add_lsdtt_cli_smoke_test(
    smoke_surface_metrics_true_north_frame
    lsdtt-surface-metrics
    surface_metrics_true_north_frame.param.in
    tnframe_direct_polyfit_aspect_degrees.tif
    --expect-raster-stat tnframe_direct_polyfit_aspect_degrees.tif:min:199.98:200.0
    --expect-raster-stat tnframe_direct_polyfit_aspect_degrees.tif:max:200.0:200.02
    --stdout-contains "${_tn_frame}"
    --stdout-contains "${_tn_on}"
)
add_lsdtt_cli_smoke_test(
    smoke_surface_metrics_grid_north_frame
    lsdtt-surface-metrics
    surface_metrics_grid_north_frame.param.in
    gnframe_direct_polyfit_aspect_degrees.tif
    --expect-raster-stat gnframe_direct_polyfit_aspect_degrees.tif:min:239.86:239.88
    --expect-raster-stat gnframe_direct_polyfit_aspect_degrees.tif:max:239.86:239.88
    --stdout-contains "bearings are relative to grid north (bearings_relative_to_true_north is false)"
)
add_lsdtt_cli_smoke_test(
    smoke_surface_metrics_true_north_utm
    lsdtt-surface-metrics
    surface_metrics_true_north_utm.param.in
    tnutm_direct_polyfit_aspect_degrees.tif
    --expect-raster-stat tnutm_direct_polyfit_aspect_degrees.tif:min:199.98:200.0
    --expect-raster-stat tnutm_direct_polyfit_aspect_degrees.tif:max:200.0:200.02
    --stdout-contains "Grid north is 0.7 deg clockwise of true north at the grid centre"
)
add_lsdtt_cli_smoke_test(
    smoke_feature_finder_true_north_structure_tensor
    lsdtt-feature-finder
    feature_finder_true_north_structure_tensor.param.in
    tnst_orientation_deg.tif
    --expect-raster-stat tnst_orientation_deg.tif:min:105.0:110.0
    --expect-raster-stat tnst_orientation_deg.tif:max:110.0:115.0
    --stdout-contains "${_tn_on}"
)
add_lsdtt_cli_smoke_test(
    smoke_ridgeline_true_north
    lsdtt-ridgeline-analysis
    ridgeline_true_north.param.in
    tnrl_lambda_90m_direct_polyfit_aspect_degrees.tif
    --expect-raster-stat tnrl_lambda_90m_direct_polyfit_aspect_degrees.tif:min:199.98:200.0
    --expect-raster-stat tnrl_lambda_90m_direct_polyfit_aspect_degrees.tif:max:200.0:200.02
    --stdout-contains "${_tn_on}"
)
add_lsdtt_cli_smoke_test(
    smoke_basin_surface_metrics_true_north
    lsdtt-basin-surface-metrics
    basin_surface_metrics_true_north.param.in
    tnbsm_aspect_degrees_lambda90.tif
    --stdout-contains "${_tn_frame}"
    --stdout-contains "${_tn_on}"
)
add_lsdtt_cli_smoke_test(
    smoke_basin_stats_true_north
    lsdtt-basin-stats
    basin_stats_true_north.param.in
    tnbs_basin_stats.csv
    --expect-file-contains "tnbs_basin_stats.csv:com_bearing_deg"
    --expect-file-contains "tnbs_basin_stats.csv:orientation_deg"
    --stdout-contains "${_tn_on}"
)
add_lsdtt_cli_smoke_test(
    smoke_basin_extraction_true_north
    lsdtt-basin-extraction
    basin_extraction_true_north.param.in
    tnbe_line_crossing_basins.csv
    --expect-file-contains "tnbe_line_crossing_basins.csv:outlet_to_com_bearing_deg"
    --stdout-contains "(6 basins)"
    --stdout-contains "${_tn_on}"
)
add_lsdtt_cli_smoke_test(
    smoke_channel_network_true_north
    lsdtt-channel-network-analysis
    channel_network_true_north.param.in
    tncn_channel_link_bearings.csv
    --expect tncn_junction_angles.csv
    --expect tncn_junction_angle_bearing_vectors.fgb
    --expect tncn_junction_com_bearings.csv
    --expect tncn_channel_bearing_group_summary.csv
    --stdout-contains "${_tn_on}"
)
add_lsdtt_cli_smoke_test(
    smoke_flat_surface_finder_true_north
    lsdtt-flat-surface-finder
    flat_surface_finder_true_north.param.in
    tnfsf_plane_fits.csv
    --stdout-contains "${_tn_on}"
)
add_lsdtt_cli_smoke_test(
    smoke_fluvial_surfaces_true_north
    lsdtt-fluvial-surfaces
    fluvial_surfaces_true_north.param.in
    tnfs_surface_patches.csv
    --expect-file-contains "tnfs_surface_patches.csv:terrain_aspect_degrees"
    --stdout-contains "${_tn_on}"
)
add_lsdtt_cli_smoke_test(
    smoke_spectral_true_north
    lsdtt-spectral-analysis
    spectral_true_north.param.in
    tnsp_significant_frequencies_2d.csv
    --expect-file-contains "tnsp_significant_frequencies_2d.csv:# grid_convergence_deg_at_grid_centre,39.8727"
    --stdout-contains "${_tn_on}"
)

add_lsdtt_cli_smoke_test(
    smoke_drainage_analysis
    lsdtt-drainage-analysis
    drainage_analysis_smoke.param.in
    tiny_d8_drainage_area.tif
)

add_lsdtt_cli_smoke_test(
    rigorous_drainage_analysis_canonical_threshold
    lsdtt-drainage-analysis
    drainage_analysis_canonical_threshold.param.in
    drainage_canonical_height_above_channel.tif
    --stdout-contains "Channel threshold: 3 contributing pixels"
    --stderr-not-contains "deprecated"
    --expect-raster-stat drainage_canonical_height_above_channel.tif:min:0:0
    --expect-raster-stat drainage_canonical_height_above_channel.tif:max:1:1000
    --expect-file-contains "drainage_canonical_ingested_params.param:threshold_contributing_pixels: 3"
    --expect-file-not-contains "drainage_canonical_ingested_params.param:channel_threshold_contributing_pixels"
)

add_lsdtt_cli_smoke_test(
    rigorous_drainage_analysis_legacy_threshold
    lsdtt-drainage-analysis
    drainage_analysis_legacy_threshold.param.in
    drainage_legacy_height_above_channel.tif
    --stdout-contains "Channel threshold: 3 contributing pixels"
    --stderr-contains "'channel_threshold_contributing_pixels' is deprecated"
    --stderr-contains "use 'threshold_contributing_pixels'"
    --expect-file-contains "drainage_legacy_ingested_params.param:threshold_contributing_pixels: 3"
    --expect-raster-stat drainage_legacy_height_above_channel.tif:min:0:0
    --expect-raster-stat drainage_legacy_height_above_channel.tif:max:1:1000
)

lsdtt_set_tests_properties(rigorous_drainage_analysis_canonical_threshold PROPERTIES
    FIXTURES_SETUP drainage_analysis_canonical_threshold)
lsdtt_set_tests_properties(rigorous_drainage_analysis_legacy_threshold PROPERTIES
    FIXTURES_SETUP drainage_analysis_legacy_threshold)
# `cmake -E compare_files` itself, but on two Python-harness outputs, so it is
# registered with the Python tests: without them its fixtures have no setup.
if(LSDTT_ENABLE_PYTHON_TESTS)
    add_test(
        NAME rigorous_drainage_analysis_threshold_equivalence
        COMMAND ${CMAKE_COMMAND} -E compare_files
            "${CMAKE_CURRENT_BINARY_DIR}/rigorous_drainage_analysis_canonical_threshold/out/drainage_canonical_height_above_channel.tif"
            "${CMAKE_CURRENT_BINARY_DIR}/rigorous_drainage_analysis_legacy_threshold/out/drainage_legacy_height_above_channel.tif"
    )
else()
    lsdtt_mark_test_unregistered(python rigorous_drainage_analysis_threshold_equivalence)
endif()
lsdtt_set_tests_properties(rigorous_drainage_analysis_threshold_equivalence PROPERTIES
    LABELS "cli"
    FIXTURES_REQUIRED "drainage_analysis_canonical_threshold;drainage_analysis_legacy_threshold")

# C9 extraction guard: unlike drainage_analysis_smoke.param.in, this reaches
# both STI P95 Jenks call sites. The raster and diagnostic oracles were frozen
# against the original CLI-local classify_jenks_categories implementation.
# Pins the TWI raster byte-for-byte across the lift of the TWI computation out
# of cli/lsdtt-drainage-analysis.cpp and into lsdtt_wetness_index_cpu.
# The checksum was captured from the
# pre-refactor CLI; the move is a pure refactor and must not change it.
#
# rigorous_drainage_sti_jenks below guards the same computation indirectly,
# since STI consumes TWI, but it pins derived Jenks breaks rather than the TWI
# raster itself. This test is the direct guard.
add_lsdtt_cli_smoke_test(
    rigorous_drainage_twi_refactor_identical
    lsdtt-drainage-analysis
    drainage_analysis_twi_pin.param.in
    tiny_twi.tif
    --expect-raster-checksum tiny_twi.tif:65360
    --expect-raster-stat tiny_twi.tif:min:-0.9906:-0.9904
    --expect-raster-stat tiny_twi.tif:max:-0.2974:-0.2973
)

# The STI parameter names were migrated to the documented <what>_fname form and
# to spelled-out units. This asserts the old names still work, still warn, and
# still produce the SAME classified output as the canonical sibling below --
# which is the whole promise of set_deprecated_aliases().
add_lsdtt_cli_smoke_test(
    rigorous_drainage_sti_deprecated_names
    lsdtt-drainage-analysis
    drainage_analysis_sti_deprecated_names.param.in
    drainage_sti_sti_p95_jenks5.tif
    --expect-raster-values drainage_sti_sti_p95_jenks5.tif:1
    --stderr-contains "'sti_soil_vector_path' is deprecated"
    --stderr-contains "'sti_ksat_csv_path' is deprecated"
    --stderr-contains "use 'sti_soil_fname'"
    --stdout-contains "STI ksd Jenks upper breaks               : -4.90252"
)

# The same STI run with the polygons in EPSG:4326: the only end-to-end test of
# the layer-to-raster CRS transform (consolidation C03). Same pins as the UTM
# sibling below, because reprojecting the inputs must not change the result.
add_lsdtt_cli_smoke_test(
    rigorous_drainage_sti_reprojected_layers
    lsdtt-drainage-analysis
    drainage_analysis_sti_reprojected.param.in
    drainage_sti_rp_sti_p95_jenks5.tif
    --expect-raster-values drainage_sti_rp_sti_p95_jenks5.tif:1
    --stdout-contains "STI ksd Jenks classified pixels          : 3 / 3"
    --stdout-contains "STI ksd Jenks upper breaks               : -4.90252"
)

add_lsdtt_cli_smoke_test(
    rigorous_drainage_sti_jenks
    lsdtt-drainage-analysis
    drainage_analysis_sti_jenks.param.in
    drainage_sti_sti_p95_jenks5.tif
    --expect drainage_sti_sti_proxy_p95_jenks5.tif
    --expect-raster-values drainage_sti_sti_p95_jenks5.tif:1
    --expect-raster-values drainage_sti_sti_proxy_p95_jenks5.tif:1
    --stdout-contains "STI ksd Jenks classes requested/effective: 5 / 1"
    --stdout-contains "STI ksd Jenks classified pixels          : 3 / 3"
    --stdout-contains "STI ksd Jenks GVF                        : 1"
    --stdout-contains "STI ksd Jenks upper breaks               : -4.90252"
    --stdout-contains "STI proxy Jenks classes requested/effective: 5 / 1"
    --stdout-contains "STI proxy Jenks classified pixels          : 1 / 1"
    --stdout-contains "STI proxy Jenks GVF                        : 1"
    --stdout-contains "STI proxy Jenks upper breaks               : 0.906619"
)

# =============================================================================
# lsdtt-distributed-hydrology
#
# The physics is all in the modules and is tested by the C++ drivers
# (unit_topmodel_column, unit_hydro_velocity_field, unit_hydro_routing,
# unit_canopy_interception, unit_infiltration_excess, unit_channel_extent,
# unit_hydro_roughness, unit_landcover_scenario and unit_topmodel_dynamic).
# These CLI tests carry only what is genuinely CLI: that the
# pipeline assembles, that the outputs appear, and that the failure paths fail
# with a message that says what to do.
#
# tiny_dem.asc is 5x5 with 1 m cells, so twi_window_radius_metres is 2.0
# everywhere here: that leaves the interior 3x3 valid, the same nine pixels the
# TWI refactor guard uses.
# =============================================================================

# The parameter declarations must all render; -h needs no parameter file.
add_test(NAME smoke_distributed_hydrology_help
         COMMAND lsdtt-distributed-hydrology -h)
lsdtt_set_tests_properties(smoke_distributed_hydrology_help PROPERTIES
    LABELS "cli;distributed_hydrology"
    PASS_REGULAR_EXPRESSION "saturated_zone_mode")

# The end-to-end run.
#
# The TWI checksum is the point of the test, not decoration: it is the same
# 65360 that rigorous_drainage_twi_refactor_identical pins for
# lsdtt-drainage-analysis on the same DEM at the same window radius. The two
# CLIs share one TWI implementation (lsdtt_wetness_index_cpu), and a second,
# divergent TWI would be a defect rather than a variant -- so this asserts they
# really do produce the same raster, which is the only mechanical check that
# the sharing has not quietly come apart.
add_lsdtt_cli_smoke_test(
    smoke_distributed_hydrology_run
    lsdtt-distributed-hydrology
    distributed_hydrology_smoke.param.in
    tinyhyd_discharge.csv
    --expect tinyhyd_water_balance.csv
    --expect tinyhyd_iuh.csv
    --expect tinyhyd_twi.tif
    --expect tinyhyd_velocity.tif
    --expect tinyhyd_traveltime.tif
    --expect tinyhyd_saturation.tif
    --expect tinyhyd_water_table.tif
    --expect-raster-checksum tinyhyd_twi.tif:65360
    --expect-raster-size tinyhyd_velocity.tif:5x5
    --expect-file-contains "tinyhyd_discharge.csv:timestep,time_seconds,rainfall_metres,pet_metres,canopy_evaporation_metres,actual_et_metres,hortonian_runoff_metres,saturation_excess_metres,baseflow_metres,total_runoff_metres,mean_deficit_metres,contributing_area_fraction,discharge_m3_per_second"
    --expect-file-contains "tinyhyd_water_balance.csv:closure_error"
    # No not-contains guard for "inf": the substring lives inside "rainfall",
    # so the check can never pass. "nan" has no such clash and is the one that
    # matters anyway -- a nodata sentinel reaching a data column prints as nan.
    --expect-file-not-contains "tinyhyd_discharge.csv:nan"
    --stdout-contains "TOPMODEL (index)"
    --stdout-contains "Pure translation through the travel-time distribution"
    --stdout-contains "DERIVED from the geometry, not chosen"
)

add_lsdtt_cli_smoke_test(
    failure_distributed_hydrology_bad_routing_method
    lsdtt-distributed-hydrology
    distributed_hydrology_bad_routing.param.in
    badrouting_discharge.csv
    --expect-failure
    --stderr-contains "FATAL: routing_method must be 'width_function'"
)

add_lsdtt_cli_smoke_test(
    failure_distributed_hydrology_missing_rainfall
    lsdtt-distributed-hydrology
    distributed_hydrology_missing_rainfall.param.in
    norain_discharge.csv
    --expect-failure
    --stderr-contains "FATAL: rainfall_fname is required"
)

# The message is the substance of this one, not the throw. A table with an
# unknown field throws too, so a bare "it failed" check would still pass with
# the dedicated branch deleted -- and the branch exists precisely to tell
# someone arriving from SD-TOPMODEL's parameterisation where conductivity
# actually comes from.
add_lsdtt_cli_smoke_test(
    failure_distributed_hydrology_landcover_table_ks
    lsdtt-distributed-hydrology
    distributed_hydrology_landcover_ks.param.in
    lcks_discharge.csv
    --expect-failure
    --stderr-contains "FATAL: landcover_table is malformed"
    --stderr-contains "Conductivity is NOT a land-cover property here"
    --stderr-contains "set 'soil=' instead"
)

# --- column_model: hype ------------------------------------------------------
#
# The HYPE numerics are tested by the C++ drivers (unit_hype_parameters,
# unit_hype_column, unit_hype_routing, unit_hype_model against the HYPE 5.36.0
# oracle; unit_subbasin_partition and unit_hype_raster_setup for the set-up from
# rasters). These carry only what is CLI: the parameter plumbing, the outlet
# modes, the output tables and the failure paths.
#
# They read the generated input set of fixture_hype_cli (a 60x60 Y-network DEM,
# land use, soil, the four parameter tables and 120 days of forcing with snow),
# so --read-path overrides the tests/data default the function supplies:
# run_cli_smoke.py keeps the last --read-path it is given.

# User outlets: a gauge on the west arm and one on the trunk, snapped downstream,
# make a chain of two sub-basins (1 drains to 2). The trunk gauge is ten rows
# above the southern edge, so the 9 rows below it (540 pixels) are outside the
# model and must be reported. Observed PET. Also writes the sub-basin raster
# (ids 1 and 2, nodata elsewhere) and the HYPE-format GeoData/GeoClass.
add_lsdtt_cli_smoke_test(
    smoke_distributed_hydrology_hype_outlets
    lsdtt-distributed-hydrology
    distributed_hydrology_hype_outlets.param.in
    hypeout_hype_discharge.csv
    --read-path "${LSDTT_HYPE_CLI_FIXTURE_DIR}"
    --expect hypeout_hype_subbasins.csv
    --expect hypeout_hype_water_balance.csv
    --expect hypeout_hype_subbasins.tif
    --expect hypeout_hype_geodata/GeoData.txt
    --expect hypeout_hype_geodata/GeoClass.txt
    --expect-file-contains "hypeout_hype_discharge.csv:timestep,date,cout_1,cout_2,crun_1,crun_2"
    # 120 daily rows from 2001-01-01: the date column is the day count, not a copy.
    --expect-file-contains "hypeout_hype_discharge.csv:120,2001-04-30,"
    --expect-file-contains "hypeout_hype_subbasins.csv:subid,maindown,outlet_x,outlet_y,area_m2,elev_mean_m,elev_std_m,slope_mean_percent,loc_rivlen_m,rivlen_m,n_classes"
    --expect-file-contains "hypeout_hype_water_balance.csv:quantity,value,units"
    --expect-file-contains "hypeout_hype_water_balance.csv:closure_error,"
    --expect-file-not-contains "hypeout_hype_discharge.csv:nan"
    --expect-file-not-contains "hypeout_hype_subbasins.csv:nan"
    --expect-raster-size hypeout_hype_subbasins.tif:60x60
    --expect-raster-values hypeout_hype_subbasins.tif:1,2
    --stdout-contains "Column model: HYPE"
    --stdout-contains "Sub-basins         : 2"
    --stdout-contains "Unassigned pixels  : 540 of 3600 valid"
    --stdout-contains "PET                : observed, 120 days"
    --stdout-contains "HYPE water balance closure error:"
)
lsdtt_set_tests_properties(smoke_distributed_hydrology_hype_outlets PROPERTIES
    LABELS "cli;distributed_hydrology"
    FIXTURES_REQUIRED hype_cli_fixture)

# Strahler order: every first-order junction basin, i.e. the two arms of the Y,
# each its own domain outlet. HYPE's calculated PET (no pet_fname).
add_lsdtt_cli_smoke_test(
    smoke_distributed_hydrology_hype_order
    lsdtt-distributed-hydrology
    distributed_hydrology_hype_order.param.in
    hypeorder_hype_discharge.csv
    --read-path "${LSDTT_HYPE_CLI_FIXTURE_DIR}"
    --expect hypeorder_hype_subbasins.csv
    --expect hypeorder_hype_water_balance.csv
    --expect-file-contains "hypeorder_hype_discharge.csv:timestep,date,cout_1,cout_2,crun_1,crun_2"
    --expect-file-contains "hypeorder_hype_subbasins.csv:subid,maindown,outlet_x,outlet_y,area_m2"
    --expect-file-contains "hypeorder_hype_water_balance.csv:closure_error,"
    --expect-file-not-contains "hypeorder_hype_discharge.csv:nan"
    --expect-file-not-contains "hypeorder_hype_subbasins.csv:nan"
    --stdout-contains "Column model: HYPE"
    --stdout-contains "Outlet mode        : select_all_basins_of_stream_order (1)"
    --stdout-contains "Sub-basins         : 2"
    --stdout-contains "Domain outlets     : 2"
    --stdout-contains "PET                : calculated by HYPE"
    --stdout-contains "HYPE water balance closure error:"
)
lsdtt_set_tests_properties(smoke_distributed_hydrology_hype_order PROPERTIES
    LABELS "cli;distributed_hydrology"
    FIXTURES_REQUIRED hype_cli_fixture)

# Strahler order range 1..: the trunk reaches the DEM edge, so the selection
# holds a base-level outlet whose own sub-basin would be only the edge strip
# (all-nodata gradient, NaN slope). The outlet helper merges it with the
# sub-basin above (decided 2026-10-04), so this must RUN, with finite values.
# The param file leaves dt_seconds unset: the HYPE path defaults to daily.
add_lsdtt_cli_smoke_test(
    smoke_distributed_hydrology_hype_order_range
    lsdtt-distributed-hydrology
    distributed_hydrology_hype_order_range.param.in
    hypeorderrange_hype_discharge.csv
    --read-path "${LSDTT_HYPE_CLI_FIXTURE_DIR}"
    --expect hypeorderrange_hype_subbasins.csv
    --expect hypeorderrange_hype_water_balance.csv
    --expect-file-not-contains "hypeorderrange_hype_discharge.csv:nan"
    --expect-file-not-contains "hypeorderrange_hype_subbasins.csv:nan"
    --stdout-contains "dt_seconds not set: using 86400 (daily)"
    --stdout-contains "HYPE water balance closure error:"
)
lsdtt_set_tests_properties(smoke_distributed_hydrology_hype_order_range PROPERTIES
    LABELS "cli;distributed_hydrology"
    FIXTURES_REQUIRED hype_cli_fixture)

# Snapping. hype_outlets_offset.csv puts both gauges 3 columns east of the
# channel. expanding_square moves each to the nearest channel pixel; downstream
# follows flow, which lands lower on the west arm (row 18 rather than 16), so
# its sub-basin 1 is larger (54200 m2 against 45300 m2). Snapping 'none' is a
# failure case below.
add_lsdtt_cli_smoke_test(
    smoke_distributed_hydrology_hype_snap_square
    lsdtt-distributed-hydrology
    distributed_hydrology_hype_snap_square.param.in
    hypesnapsq_hype_discharge.csv
    --read-path "${LSDTT_HYPE_CLI_FIXTURE_DIR}"
    --stdout-contains "Outlet 1: row=16 col=21"
    --stdout-contains "Outlet 2: row=50 col=30"
    --stdout-contains "Sub-basins         : 2"
    --expect-file-contains "hypesnapsq_hype_subbasins.csv:,45300,"
    --expect-file-not-contains "hypesnapsq_hype_discharge.csv:nan"
)
lsdtt_set_tests_properties(smoke_distributed_hydrology_hype_snap_square PROPERTIES
    LABELS "cli;distributed_hydrology"
    FIXTURES_REQUIRED hype_cli_fixture)

add_lsdtt_cli_smoke_test(
    smoke_distributed_hydrology_hype_snap_downstream
    lsdtt-distributed-hydrology
    distributed_hydrology_hype_snap_downstream.param.in
    hypesnapds_hype_discharge.csv
    --read-path "${LSDTT_HYPE_CLI_FIXTURE_DIR}"
    --stdout-contains "Outlet 1: row=18 col=22"
    --stdout-contains "Outlet 2: row=50 col=30"
    --stdout-contains "Sub-basins         : 2"
    --expect-file-contains "hypesnapds_hype_subbasins.csv:,54200,"
    --expect-file-not-contains "hypesnapds_hype_discharge.csv:nan"
)
lsdtt_set_tests_properties(smoke_distributed_hydrology_hype_snap_downstream PROPERTIES
    LABELS "cli;distributed_hydrology"
    FIXTURES_REQUIRED hype_cli_fixture)

# largest_area with a 3-pixel half-width takes the pixel with the most
# contributing pixels in each 7 x 7 window. On the west arm that is the
# downstream pixel's (18,22), so sub-basin 1 matches downstream's 54200 m2;
# outlet 2 moves 3 rows down the trunk to (53,30). Both pixels were checked
# against a brute-force window maximum of the contributing-pixel counts
# (2026-10-23).
add_lsdtt_cli_smoke_test(
    smoke_distributed_hydrology_hype_snap_largest_area
    lsdtt-distributed-hydrology
    distributed_hydrology_hype_snap_largest_area.param.in
    hypesnapla_hype_discharge.csv
    --read-path "${LSDTT_HYPE_CLI_FIXTURE_DIR}"
    --stdout-contains "(largest_area snapping)"
    --stdout-contains "Outlet 1: row=18 col=22"
    --stdout-contains "Outlet 2: row=53 col=30"
    --stdout-contains "Sub-basins         : 2"
    --expect-file-contains "hypesnapla_hype_subbasins.csv:,54200,"
    --expect-file-not-contains "hypesnapla_hype_discharge.csv:nan"
)
lsdtt_set_tests_properties(smoke_distributed_hydrology_hype_snap_largest_area PROPERTIES
    LABELS "cli;distributed_hydrology;snapping"
    FIXTURES_REQUIRED hype_cli_fixture)

# Nodata in a class raster. A 10 x 5 block of nodata soil straddles the soil
# boundary inside the trunk sub-basin (25 pixels each of classes 2 and 3). The
# 50 pixels stay in the sub-basin area (262500 m2) but leave the class
# fractions, which renormalise over 2575 pixels. The expected fractions were
# checked against an independent pixel count of the rasters (2026-10-04):
# 465/2575, 875/2575, 605/2575, 630/2575.
add_lsdtt_cli_smoke_test(
    smoke_distributed_hydrology_hype_class_holes
    lsdtt-distributed-hydrology
    distributed_hydrology_hype_class_holes.param.in
    hypeholes_hype_discharge.csv
    --read-path "${LSDTT_HYPE_CLI_FIXTURE_DIR}"
    --stdout-contains "Pixels with no class: 50"
    --expect-file-contains "hypeholes_hype_subbasins.csv:,262500,"
    --expect-file-contains "hypeholes_hype_geodata/GeoData.txt:0.1805825243\t0.3398058252\t0.2349514563\t0.2446601942"
    --expect-file-not-contains "hypeholes_hype_discharge.csv:nan"
)
lsdtt_set_tests_properties(smoke_distributed_hydrology_hype_class_holes PROPERTIES
    LABELS "cli;distributed_hydrology"
    FIXTURES_REQUIRED hype_cli_fixture)

# n_timesteps shorter than the forcing runs only those days.
add_lsdtt_cli_smoke_test(
    smoke_distributed_hydrology_hype_short_run
    lsdtt-distributed-hydrology
    distributed_hydrology_hype_short_run.param.in
    hypeshort_hype_discharge.csv
    --read-path "${LSDTT_HYPE_CLI_FIXTURE_DIR}"
    --stdout-contains "Steps to run       : 30"
    --expect-file-contains "hypeshort_hype_discharge.csv:30,2001-01-30,"
    --expect-file-not-contains "hypeshort_hype_discharge.csv:31,2001-01-31,"
)
lsdtt_set_tests_properties(smoke_distributed_hydrology_hype_short_run PROPERTIES
    LABELS "cli;distributed_hydrology"
    FIXTURES_REQUIRED hype_cli_fixture)

# Failure paths. Each asserts the FATAL text, not just the exit code: the
# messages are what tell the user which restriction of the first version
# they hit, and the library errors must arrive with their own
# message (file, column, pair) intact.
set(_hype_failures
    "missing_temperature|FATAL: temperature_fname is required for column_model hype"
    "hourly|FATAL: column_model hype requires dt_seconds: 86400"
    "fsceff|but the first version supports only fsceff = 0"
    "unknown_parameter|unknown column \"cmelt\""
    "missing_class_pair|has no class for these (land use, soil) pairs found in the rasters: (1, 2)"
    "bad_column_model|FATAL: column_model must be 'topmodel' or 'hype'. Got: hbv"
    "no_outlet_mode|FATAL: column_model hype needs exactly one outlet mode"
    "snap_none_off_channel|basin outlet at row 15, col 23 is not on the channel network"
    "class_nodata_subbasin|sub-basin 1 (outlet flat 920) has no pixel with a valid class"
    "forcing_length|FATAL: the forcing series differ in length (rainfall 120, temperature 119, PET 120 rows)"
    "n_timesteps_too_long|FATAL: n_timesteps (121) is longer than the forcing (120 rows)"
    "bad_date_invalid|FATAL: hype_start_date must be the date of the first forcing row as YYYY-MM-DD (e.g. 2001-10-01). Got: '2001-02-30'"
    "bad_date_format|FATAL: hype_start_date must be the date of the first forcing row as YYYY-MM-DD (e.g. 2001-10-01). Got: '01/01/2001'"
)
foreach(_case IN LISTS _hype_failures)
    string(FIND "${_case}" "|" _bar)
    string(SUBSTRING "${_case}" 0 ${_bar} _name)
    math(EXPR _msg_at "${_bar} + 1")
    string(SUBSTRING "${_case}" ${_msg_at} -1 _msg)
    add_lsdtt_cli_smoke_test(
        failure_distributed_hydrology_hype_${_name}
        lsdtt-distributed-hydrology
        distributed_hydrology_hype_${_name}.param.in
        hypefail_${_name}_hype_discharge.csv
        --read-path "${LSDTT_HYPE_CLI_FIXTURE_DIR}"
        --expect-failure
        --stderr-contains "${_msg}"
    )
    lsdtt_set_tests_properties(failure_distributed_hydrology_hype_${_name} PROPERTIES
        LABELS "cli;distributed_hydrology"
        FIXTURES_REQUIRED hype_cli_fixture)
endforeach()

# lsdtt-mountain-front and lsdtt-feature-finder had NO tests at all until
# 2026-09-10. That was found when the B7 extraction (write_masked_pixels_as_points)
# refactored the outlet writer of one and both point writers of the other, and
# the change had to be verified by hand against a stashed binary because nothing
# in the suite would have caught a regression. feature-finder is the largest CLI
# in the repo at 174 parameters. These two cover the parameter plumbing and one
# output each; they are a floor, not adequate coverage.
# =============================================================================
# lsdtt-basin-surface-metrics
#
# No tests until 2026-09-10. C12 -- the outline curvature sign-flip counting and
# the slope-area analysis, ~160 lines the helper audit calls "the tool's actual
# science" -- is still inline in main(), so these tests are the guard for that
# extraction.
#
# Two cases rather than one, because the per-basin outline work and the
# whole-DEM polyfit raster path are different branches; one case with every
# switch on could not tell a failure in one from a failure in the other.
#
# Both reuse the basin-stats fixture DEM. The outlet set is committed separately
# and is INTERIOR ONLY: see basin_surface_metrics_smoke.param.in for why.
# =============================================================================
if(LSDTT_ENABLE_PYTHON_TESTS)
    lsdtt_add_harness_test(smoke_basin_surface_metrics
            --exe "$<TARGET_FILE:lsdtt-basin-surface-metrics>"
            --cli-src-dir "${CMAKE_SOURCE_DIR}/cli"
            --param-template "${LSDTT_TEST_PARAM_DIR}/basin_surface_metrics_smoke.param.in"
            --read-path "${LSDTT_BASIN_STATS_FIXTURE_DIR}/out"
            --work-dir "${CMAKE_CURRENT_BINARY_DIR}/smoke_basin_surface_metrics"
            --expect bsm_outline_curvature_flips.csv
            --expect bsm_slope_area_stats.csv
            --expect bsm_basin_polygons.fgb
            --expect-vector-feature-min bsm_basin_polygons.fgb:20
            # Both lambda scales must appear as their own columns: the per-lambda
            # expansion is part of what C12 moves, and a single-scale run would not
            # exercise it.
            --expect-file-contains "bsm_outline_curvature_flips.csv:N_flips_lambda60"
            --expect-file-contains "bsm_outline_curvature_flips.csv:N_flips_lambda120"
            --expect-file-contains "bsm_slope_area_stats.csv:A_max_slope_lambda120"
            # The warn-and-skip path: one outlet in this set cannot be snapped, and
            # the run must carry on rather than abort.
            --expect-file-contains "../stderr.log:could not be snapped downstream"
            # ORACLES for C12, not numeric standards: they say the sign-flip counting
            # and the slope-area analysis did not change when those ~160 lines moved
            # out of main(). Re-freeze deliberately if a change is intended.
            --expect-sha256 bsm_outline_curvature_flips.csv:994b77c5706e32d20ba4a9f7f1eba1c0516d9ad9d489113d3155151ec74f791f
            --expect-sha256 bsm_slope_area_stats.csv:2abfaff6b6482336c50a141f5d3aad70e974f52b25441028b474cc498bdaf949)
else()
    lsdtt_mark_test_unregistered(python smoke_basin_surface_metrics)
endif()
lsdtt_set_tests_properties(smoke_basin_surface_metrics PROPERTIES
    LABELS "cli;basin_surface_metrics"
    FIXTURES_REQUIRED basin_stats_fixture)

if(LSDTT_ENABLE_PYTHON_TESTS)
    lsdtt_add_harness_test(smoke_basin_surface_metrics_polyfit
            --exe "$<TARGET_FILE:lsdtt-basin-surface-metrics>"
            --cli-src-dir "${CMAKE_SOURCE_DIR}/cli"
            --param-template "${LSDTT_TEST_PARAM_DIR}/basin_surface_metrics_polyfit.param.in"
            --read-path "${LSDTT_BASIN_STATS_FIXTURE_DIR}/out"
            --work-dir "${CMAKE_CURRENT_BINARY_DIR}/smoke_basin_surface_metrics_polyfit"
            --expect bsmp_slope_lambda60.tif
            --expect bsmp_maximum_curvature_lambda60.tif
            --expect bsmp_dinf_area.tif
            # One raster per requested metric, named for the lambda: that naming is
            # what a caller depends on to find the output of a multi-scale run.
            --expect-raster-size bsmp_slope_lambda60.tif:200x200)
else()
    lsdtt_mark_test_unregistered(python smoke_basin_surface_metrics_polyfit)
endif()
lsdtt_set_tests_properties(smoke_basin_surface_metrics_polyfit PROPERTIES
    LABELS "cli;basin_surface_metrics"
    FIXTURES_REQUIRED basin_stats_fixture)

if(LSDTT_ENABLE_PYTHON_TESTS)
    lsdtt_add_harness_test(smoke_basin_surface_metrics_expanding_square
            --exe "$<TARGET_FILE:lsdtt-basin-surface-metrics>"
            --cli-src-dir "${CMAKE_SOURCE_DIR}/cli"
            --param-template "${LSDTT_TEST_PARAM_DIR}/basin_surface_metrics_expanding_square.param.in"
            --read-path "${LSDTT_BASIN_STATS_FIXTURE_DIR}/out"
            --work-dir "${CMAKE_CURRENT_BINARY_DIR}/smoke_basin_surface_metrics_expanding_square"
            --expect bsm_square_outline_curvature_flips.csv)
else()
    lsdtt_mark_test_unregistered(python smoke_basin_surface_metrics_expanding_square)
endif()
lsdtt_set_tests_properties(smoke_basin_surface_metrics_expanding_square PROPERTIES
    LABELS "cli;basin_surface_metrics;snapping"
    FIXTURES_REQUIRED basin_stats_fixture)

if(LSDTT_ENABLE_PYTHON_TESTS)
    lsdtt_add_harness_test(rigorous_basin_surface_metrics_legacy_threshold
            --exe "$<TARGET_FILE:lsdtt-basin-surface-metrics>"
            --cli-src-dir "${CMAKE_SOURCE_DIR}/cli"
            --param-template "${LSDTT_TEST_PARAM_DIR}/basin_surface_metrics_legacy_threshold.param.in"
            --read-path "${LSDTT_BASIN_STATS_FIXTURE_DIR}/out"
            --work-dir "${CMAKE_CURRENT_BINARY_DIR}/rigorous_basin_surface_metrics_legacy_threshold"
            --expect bsm_legacy_outline_curvature_flips.csv
            --expect-sha256 bsm_legacy_outline_curvature_flips.csv:f2f2537d2cb61a60da2fc2d4a1c862645424298d475ffcd1d4abf74010894181
            --stderr-contains "'source_pixels_threshold' is deprecated"
            --stderr-contains "use 'threshold_contributing_pixels'"
            --stderr-contains "'basin_outlet_csv' is deprecated"
            --stderr-contains "use 'basin_outlet_fname'")
else()
    lsdtt_mark_test_unregistered(python rigorous_basin_surface_metrics_legacy_threshold)
endif()
lsdtt_set_tests_properties(rigorous_basin_surface_metrics_legacy_threshold PROPERTIES
    LABELS "cli;basin_surface_metrics;compatibility"
    FIXTURES_REQUIRED basin_stats_fixture)

add_lsdtt_cli_smoke_test(
    failure_basin_surface_metrics_unknown_snap
    lsdtt-basin-surface-metrics
    basin_surface_metrics_unknown_snap.param.in
    unused
    --expect-failure
    --stderr-contains "FATAL: basin_outlet_snapping_method must be one of: none, downstream, expanding_square, largest_area."
)

# =============================================================================
# lsdtt-ridgeline-analysis: the range_defined_ridgelines path (C2 and C3)
#
# The three pre-existing ridgeline tests are all paired_flowpaths and reach
# neither resolve_dem_corners/extract_boundary_segments (C2) nor
# run_raw_dem_3x3_polyfit_with_bcs (C3). These two do.
#
# oooo and opop are both run because C3's substance is the PERIODIC wrap-around
# padding: with fully open boundaries there is nothing to wrap and the branch is
# never taken. The two produce different slope rasters.
# =============================================================================
if(LSDTT_ENABLE_PYTHON_TESTS)
    lsdtt_add_harness_test(smoke_ridgeline_range_open
            --exe "$<TARGET_FILE:lsdtt-ridgeline-analysis>"
            --cli-src-dir "${CMAKE_SOURCE_DIR}/cli"
            --param-template "${LSDTT_TEST_PARAM_DIR}/ridgeline_range_oooo.param.in"
            --read-path "${LSDTT_BASIN_STATS_FIXTURE_DIR}/out"
            --work-dir "${CMAKE_CURRENT_BINARY_DIR}/smoke_ridgeline_range_open"
            --expect rg_oooo_lambda_60m_raw_dem_3x3_polyfit_slope.tif
            --expect rg_oooo_lambda_60m_raw_dem_3x3_polyfit_maximum_curvature.tif
            --expect-raster-size rg_oooo_lambda_60m_raw_dem_3x3_polyfit_slope.tif:200x200
            # C2 must actually resolve corners and segment the boundary.
            --expect-file-contains "../stdout.log:Resolving DEM corners"
            --expect-file-contains "../stdout.log:Boundary segment corner 3 assigned"
            # ORACLE for C2+C3 under open boundaries.
            --expect-sha256 rg_oooo_lambda_60m_raw_dem_3x3_polyfit_slope.tif:c6a0d7322aa1d1f510db0fac8a2544c052238e1782d74bd68805e04a37b8b9c4)
else()
    lsdtt_mark_test_unregistered(python smoke_ridgeline_range_open)
endif()
lsdtt_set_tests_properties(smoke_ridgeline_range_open PROPERTIES
    LABELS "cli;ridgeline"
    FIXTURES_REQUIRED basin_stats_fixture)

if(LSDTT_ENABLE_PYTHON_TESTS)
    lsdtt_add_harness_test(smoke_ridgeline_range_periodic
            --exe "$<TARGET_FILE:lsdtt-ridgeline-analysis>"
            --cli-src-dir "${CMAKE_SOURCE_DIR}/cli"
            --param-template "${LSDTT_TEST_PARAM_DIR}/ridgeline_range_opop.param.in"
            --read-path "${LSDTT_BASIN_STATS_FIXTURE_DIR}/out"
            --work-dir "${CMAKE_CURRENT_BINARY_DIR}/smoke_ridgeline_range_periodic"
            --expect rg_opop_lambda_60m_raw_dem_3x3_polyfit_slope.tif
            --expect-raster-size rg_opop_lambda_60m_raw_dem_3x3_polyfit_slope.tif:200x200
            --expect-file-contains "../stdout.log:Resolving DEM corners"
            # ORACLE for the PERIODIC padding branch. This hash differs from the
            # open-boundary one above, which is the evidence that the wrap-around
            # padding in run_raw_dem_3x3_polyfit_with_bcs actually does something --
            # if a refactor made the two equal, both oracles could not hold.
            --expect-sha256 rg_opop_lambda_60m_raw_dem_3x3_polyfit_slope.tif:f7beffb5d95fb65e7ca6108624baac5f8c7ae60399622c2f235acb94a2bd9eae)
else()
    lsdtt_mark_test_unregistered(python smoke_ridgeline_range_periodic)
endif()
lsdtt_set_tests_properties(smoke_ridgeline_range_periodic PROPERTIES
    LABELS "cli;ridgeline"
    FIXTURES_REQUIRED basin_stats_fixture)

# lsdtt-chi-analysis: the outlet-FILE path (C5). The other chi tests only reach
# resolve_outlet_flats through its auto-detect branch.
if(LSDTT_ENABLE_PYTHON_TESTS)
    lsdtt_add_harness_test(smoke_chi_analysis_outlet_file
            --exe "$<TARGET_FILE:lsdtt-chi-analysis>"
            --cli-src-dir "${CMAKE_SOURCE_DIR}/cli"
            --param-template "${LSDTT_TEST_PARAM_DIR}/chi_analysis_outlet_file.param.in"
            --read-path "${LSDTT_BASIN_STATS_FIXTURE_DIR}/out"
            --work-dir "${CMAKE_CURRENT_BINARY_DIR}/smoke_chi_analysis_outlet_file"
            --expect chi_chi_basin_outlets.csv
            --expect chi_chi_data_maps.csv
            # The outlet FILE branch, not auto-detect.
            --expect-file-contains "../stdout.log:Reading outlet file"
            --expect-file-contains "../stdout.log:Using explicit m_over_n list with 1 value(s): 0.45"
            # ... and the per-outlet skip path within it.
            --expect-file-contains "../stderr.log:could not be snapped downstream"
            --stderr-contains "'basin_outlet_file' is deprecated"
            --stderr-contains "use 'basin_outlet_fname'"
            # ORACLE for C5: the resolved outlet set, which is what
            # resolve_outlet_flats exists to produce. Checked deterministic before
            # freezing.
            --expect-sha256 chi_chi_basin_outlets.csv:f2a980bf03e6870a6dc754386822b4ce6731017bacc97ff23ffa9be202f29621)
else()
    lsdtt_mark_test_unregistered(python smoke_chi_analysis_outlet_file)
endif()
lsdtt_set_tests_properties(smoke_chi_analysis_outlet_file PROPERTIES
    LABELS "cli;chi"
    FIXTURES_REQUIRED basin_stats_fixture)

# lsdtt-linear-stats: a floor, not coverage. See the param template for why the
# line fixture must stay clear of the raster edge.
add_test(NAME fixture_linear_stats_lines
    COMMAND ${CMAKE_COMMAND} -E copy
        "${LSDTT_TEST_DATA_DIR}/linear_stats_lines.fgb"
        "${LSDTT_BASIN_STATS_FIXTURE_DIR}/out/linear_stats_lines.fgb")
lsdtt_set_tests_properties(fixture_linear_stats_lines PROPERTIES
    LABELS "fixture;linear_stats"
    FIXTURES_SETUP basin_stats_fixture
    DEPENDS fixture_basin_stats_dem)

if(LSDTT_ENABLE_PYTHON_TESTS)
    lsdtt_add_harness_test(smoke_linear_stats
            --exe "$<TARGET_FILE:lsdtt-linear-stats>"
            --cli-src-dir "${CMAKE_SOURCE_DIR}/cli"
            --param-template "${LSDTT_TEST_PARAM_DIR}/linear_stats_smoke.param.in"
            --read-path "${LSDTT_BASIN_STATS_FIXTURE_DIR}/out"
            --work-dir "${CMAKE_CURRENT_BINARY_DIR}/smoke_linear_stats"
            --expect ls_linear_stats.fgb
            --expect-vector-feature-min ls_linear_stats.fgb:2
            --expect-sha256 ls_linear_stats.fgb:5089ab341d0a212acc0fee3cba76c3e414aaa89bcbf7e017aabd93fa3988224a)
else()
    lsdtt_mark_test_unregistered(python smoke_linear_stats)
endif()
lsdtt_set_tests_properties(smoke_linear_stats PROPERTIES
    LABELS "cli;linear_stats"
    FIXTURES_REQUIRED basin_stats_fixture)

# =============================================================================
# lsdtt-stream-burning
#
# No tests at all until 2026-09-10. The largest CLI in the repo (~2981 lines)
# and the highest-risk extraction of the CLI clean-up (C4, the
# inward-draining edge carve), which was to be done "with the burn tests
# green" -- so C4 was blocked on this fixture existing.
#
# The vector network is COMMITTED rather than derived from the DEM. A network
# extracted from the DEM has interior channel heads, so the carve stage runs and
# carves nothing; tests/data/stream_burning_network.fgb puts four of its five
# lines on the raster boundary so the carve actually executes. Its coordinates
# are tied to the fixed synthetic georeference of the generated DEM.
# =============================================================================
if(LSDTT_ENABLE_PYTHON_TESTS)
    lsdtt_add_harness_test(fixture_stream_burning_dem
            --exe "$<TARGET_FILE:lsdtt-synthetic-surfaces>"
            --cli-src-dir "${CMAKE_SOURCE_DIR}/cli"
            --param-template "${LSDTT_TEST_PARAM_DIR}/stream_burning_fixture_dem.param.in"
            --read-path "${LSDTT_TEST_DATA_DIR}"
            --work-dir "${LSDTT_STREAM_BURNING_FIXTURE_DIR}"
            --expect sb_diamond_square.tif)
else()
    lsdtt_mark_test_unregistered(python fixture_stream_burning_dem)
endif()
lsdtt_set_tests_properties(fixture_stream_burning_dem PROPERTIES
    LABELS "fixture;stream_burning"
    FIXTURES_SETUP stream_burning_fixture)

add_test(NAME fixture_stream_burning_network
    COMMAND ${CMAKE_COMMAND} -E copy
        "${LSDTT_TEST_DATA_DIR}/stream_burning_network.fgb"
        "${LSDTT_STREAM_BURNING_FIXTURE_DIR}/out/stream_burning_network.fgb")
lsdtt_set_tests_properties(fixture_stream_burning_network PROPERTIES
    LABELS "fixture;stream_burning"
    FIXTURES_SETUP stream_burning_fixture
    DEPENDS fixture_stream_burning_dem)

if(LSDTT_ENABLE_PYTHON_TESTS)
    lsdtt_add_harness_test(smoke_stream_burning
            --exe "$<TARGET_FILE:lsdtt-stream-burning>"
            --cli-src-dir "${CMAKE_SOURCE_DIR}/cli"
            --param-template "${LSDTT_TEST_PARAM_DIR}/stream_burning_smoke.param.in"
            --read-path "${LSDTT_STREAM_BURNING_FIXTURE_DIR}/out"
            --work-dir "${CMAKE_CURRENT_BINARY_DIR}/smoke_stream_burning"
            --expect burn_burned.tif
            --expect burn_burned_network_lines.fgb
            --expect burn_merged_lines.fgb
            --expect burn_dem_network_lines.fgb
            --expect burn_cleaned_network.fgb
            # The carve must actually CARVE. Without this the stage still runs and
            # prints "Edge sources carved inward: 0", and the C4 code path would be
            # untested while the test passed. --expect-file-contains resolves against
            # the out/ directory, so ../stdout.log is the run log.
            --expect-file-contains "../stdout.log:Edge sources carved inward: 2"
            # ORACLES frozen before the C4 extraction. The stdout assertion above
            # proves the carve RAN; these prove its RESULT did not change when the
            # ~160 lines moved out of main(). Verified deterministic across repeat
            # runs before freezing.
            --expect-sha256 burn_burned.tif:f63a442d8b42bbb6a6a4c857b0c046bda2bb889bc8c319f4cc105a8ad4b8648b
            --expect-sha256 burn_dem_network_lines.fgb:b675ef0fe65bd863eff7ee81fac93d45595366fe4b76926f051bc85788ed05a2)
else()
    lsdtt_mark_test_unregistered(python smoke_stream_burning)
endif()
lsdtt_set_tests_properties(smoke_stream_burning PROPERTIES
    LABELS "cli;stream_burning"
    FIXTURES_REQUIRED stream_burning_fixture)

# Stream-burning source seeding (river-lines plan SB1): a source over the NoData
# margin of a warped DEM is moved downstream along its edge to the first valid
# cell instead of being skipped. Exact answers from a hand-built DEM.
add_executable(lsdtt-test-vector-source-seed cpp/test_vector_source_seed.cpp)
target_link_libraries(lsdtt-test-vector-source-seed PRIVATE lsdtt3_cpu)
add_test(NAME unit_stream_burning_source_seed
         COMMAND lsdtt-test-vector-source-seed "${CMAKE_CURRENT_BINARY_DIR}/stream_burning_source_seed")
lsdtt_set_tests_properties(unit_stream_burning_source_seed PROPERTIES LABELS "stream_burning")

# Corridor burning (stream_burn_mode burn_to_corridor, river-lines plan VT):
# offset valley at three corridor widths, and a ridge across the valley.
add_executable(lsdtt-test-burn-corridor cpp/test_burn_corridor.cpp)
target_link_libraries(lsdtt-test-burn-corridor PRIVATE lsdtt3_cpu)
add_test(NAME unit_stream_burning_corridor COMMAND lsdtt-test-burn-corridor)
lsdtt_set_tests_properties(unit_stream_burning_corridor PROPERTIES LABELS "stream_burning")

# Loop labels on the vector junction network (river-lines plan SB2):
# loop_reason and closes_cycle, exact on a hand-built DEM (a stolen receiver
# slot that cuts the network, and a real braid).
add_executable(lsdtt-test-vector-loop-labels cpp/test_vector_loop_labels.cpp)
target_link_libraries(lsdtt-test-vector-loop-labels PRIVATE lsdtt3_cpu)
add_test(NAME unit_stream_burning_loop_labels
         COMMAND lsdtt-test-vector-loop-labels "${CMAKE_CURRENT_BINARY_DIR}/stream_burning_loop_labels")
lsdtt_set_tests_properties(unit_stream_burning_loop_labels PROPERTIES LABELS "stream_burning")

# lsdtt-fetch-river-lines -> lsdtt-stream-burning (river-lines plan, phase 1b):
# the fixture network, and a reversed / near-miss / multipart copy of it, each
# prepared by lsdtt-fetch-river-lines, must give stream burning's reference
# vector summary and the burned DEM frozen for smoke_stream_burning above.
if(LSDTT_ENABLE_PYTHON_TESTS)
    add_test(NAME smoke_fetch_river_lines_to_stream_burning
        COMMAND ${Python3_EXECUTABLE}
            "${CMAKE_CURRENT_SOURCE_DIR}/scripts/run_fetch_river_lines_to_stream_burning.py"
            --fetch-exe "$<TARGET_FILE:lsdtt-fetch-river-lines>"
            --burn-exe "$<TARGET_FILE:lsdtt-stream-burning>"
            --burn-param-template "${LSDTT_TEST_PARAM_DIR}/stream_burning_smoke.param.in"
            --fixture-dir "${LSDTT_STREAM_BURNING_FIXTURE_DIR}/out"
            --work-dir "${CMAKE_CURRENT_BINARY_DIR}/smoke_fetch_river_lines_to_stream_burning"
            --expect-burned-sha256 f63a442d8b42bbb6a6a4c857b0c046bda2bb889bc8c319f4cc105a8ad4b8648b)
    lsdtt_mark_test_runtime(python smoke_fetch_river_lines_to_stream_burning)
else()
    lsdtt_mark_test_unregistered(python smoke_fetch_river_lines_to_stream_burning)
endif()
lsdtt_set_tests_properties(smoke_fetch_river_lines_to_stream_burning PROPERTIES
    LABELS "cli;fetch;stream_burning"
    FIXTURES_REQUIRED stream_burning_fixture)

# Least-cost-path burn: regression test for the use-after-free segfault fixed on
# 2026-10-07, and the only test that runs raster_least_cost_path. See the
# parameter file.
add_lsdtt_cli_smoke_test(
    smoke_stream_burning_lcp_burn
    lsdtt-stream-burning
    stream_burning_lcp_burn.param.in
    lcpb_lcp_burn_paths.fgb
    --read-path "${LSDTT_STREAM_BURNING_FIXTURE_DIR}/out"
    --expect lcpb_burn_segments.fgb
    --expect lcpb_lcp_burned.tif
    --expect lcpb_connected_network_lines.fgb
    --stdout-contains "LCP burn: 12 gaps (case1/2/3), 10 D8-connected on carve"
    --stdout-contains "Candidate detection complete."
)
lsdtt_set_tests_properties(smoke_stream_burning_lcp_burn PROPERTIES
    FIXTURES_REQUIRED stream_burning_fixture)

add_lsdtt_cli_smoke_test(
    rigorous_stream_burning_legacy_threshold
    lsdtt-stream-burning
    stream_burning_legacy_threshold.param.in
    burn_legacy_cleaned_network.fgb
    --read-path "${LSDTT_STREAM_BURNING_FIXTURE_DIR}/out"
    --stdout-contains "D8 channel cells (>= 3 px):"
    --stderr-contains "'d8_channel_threshold_pixels' is deprecated"
    --stderr-contains "use 'threshold_contributing_pixels'"
    --stderr-contains "'do_stream_burn' is deprecated"
    --stderr-contains "use 'run_stream_burn'"
    --expect-file-contains "burn_legacy_ingested_params.param:threshold_contributing_pixels: 3"
    --expect-file-contains "burn_legacy_ingested_params.param:run_stream_burn: true"
)
lsdtt_set_tests_properties(rigorous_stream_burning_legacy_threshold PROPERTIES
    FIXTURES_REQUIRED stream_burning_fixture)

# =============================================================================
# lsdtt-cross-sections, lsdtt-bankfull-width, lsdtt-marsh-topography,
# lsdtt-water-centreline
#
# None of these four CLIs had a test until 2026-10-07. Each is the only CLI
# route to modules the src/utils + src/cpu consolidation changes
# (lsdtt_cross_sections_cpu, lsdtt_manning_hydraulics_cpu,
# lsdtt_marsh_platform_cpu, lsdtt_tidal_channel_cpu,
# lsdtt_tidal_skeleton_metrics_cpu, lsdtt_water_area_centreline_cpu), and that
# work gates on hashing every output the CLI tests write
# (tests/scripts/capture_consolidation_e2e_hashes.py), so a CLI without a test
# had no gate at all. These smoke tests pin what each run must report; the
# numbers themselves are guarded by the e2e hashes, not frozen here.
# =============================================================================
if(LSDTT_ENABLE_PYTHON_TESTS)
    lsdtt_add_harness_test(fixture_cross_sections_dem
            --exe "$<TARGET_FILE:lsdtt-synthetic-surfaces>"
            --cli-src-dir "${CMAKE_SOURCE_DIR}/cli"
            --param-template "${LSDTT_TEST_PARAM_DIR}/cross_sections_fixture_dem.param.in"
            --read-path "${LSDTT_TEST_DATA_DIR}"
            --work-dir "${LSDTT_CROSS_SECTIONS_FIXTURE_DIR}"
            --expect xs_terrace_staircase.tif)
else()
    lsdtt_mark_test_unregistered(python fixture_cross_sections_dem)
endif()
lsdtt_set_tests_properties(fixture_cross_sections_dem PROPERTIES
    LABELS "fixture;cross_sections"
    FIXTURES_SETUP cross_sections_fixture)

add_test(NAME fixture_cross_sections_inputs
    COMMAND ${CMAKE_COMMAND} -E copy
        "${LSDTT_TEST_DATA_DIR}/cross_sections_source.csv"
        "${LSDTT_TEST_DATA_DIR}/bankfull_width_reaches.csv"
        "${LSDTT_CROSS_SECTIONS_FIXTURE_DIR}/out/")
lsdtt_set_tests_properties(fixture_cross_sections_inputs PROPERTIES
    LABELS "fixture;cross_sections"
    FIXTURES_SETUP cross_sections_fixture
    DEPENDS fixture_cross_sections_dem)

add_lsdtt_cli_smoke_test(
    smoke_cross_sections
    lsdtt-cross-sections
    cross_sections_smoke.param.in
    xsec_cross_sections.fgb
    --read-path "${LSDTT_CROSS_SECTIONS_FIXTURE_DIR}/out"
    --expect xsec_channel_network_lines.fgb
    --expect xsec_cross_section_sample_points.fgb
    --expect xsec_vertical_excess_curvature.tif
    --stdout-contains "Points mapped       : 1"
    --stdout-contains "Channel cells : 140"
    --stdout-contains "Cross sections written from 14 channel centres."
)
lsdtt_set_tests_properties(smoke_cross_sections PROPERTIES
    FIXTURES_REQUIRED cross_sections_fixture)

add_lsdtt_cli_smoke_test(
    smoke_bankfull_width
    lsdtt-bankfull-width
    bankfull_width_smoke.param.in
    bf_manning_widths.csv
    --read-path "${LSDTT_CROSS_SECTIONS_FIXTURE_DIR}/out"
    --expect bf_channel_profile.csv
    --expect bf_cross_section_sample_points.fgb
    --expect bf_manning_bank_points.fgb
    --stdout-contains "Cross sections with valid slope: 14 of 14"
    --stdout-contains "Reaches loaded                          : 2"
    --stdout-contains "Manning widths: 14 valid out of 14 sections matched to a reach by distance."
    --expect-file-not-contains "bf_manning_widths.csv:nan"
)
lsdtt_set_tests_properties(smoke_bankfull_width PROPERTIES
    FIXTURES_REQUIRED cross_sections_fixture)

# The point-file key renamed 2026-10-26 (source_points_path ->
# source_points_fname, resolved against read_path as in lsdtt-fluvial-surfaces):
# an old file must still run, to the same result, with a warning.
add_lsdtt_cli_smoke_test(
    rigorous_cross_sections_legacy_points_key
    lsdtt-cross-sections
    cross_sections_legacy_points_key.param.in
    xsec_cross_sections.fgb
    --read-path "${LSDTT_CROSS_SECTIONS_FIXTURE_DIR}/out"
    --stdout-contains "Cross sections written from 14 channel centres."
    --stderr-contains "'source_points_path' is deprecated"
    --stderr-contains "use 'source_points_fname'"
)
lsdtt_set_tests_properties(rigorous_cross_sections_legacy_points_key PROPERTIES
    FIXTURES_REQUIRED cross_sections_fixture)

add_lsdtt_cli_smoke_test(
    rigorous_bankfull_width_legacy_points_key
    lsdtt-bankfull-width
    bankfull_width_legacy_points_key.param.in
    bf_manning_widths.csv
    --read-path "${LSDTT_CROSS_SECTIONS_FIXTURE_DIR}/out"
    --stdout-contains "Manning widths: 14 valid out of 14 sections matched to a reach by distance."
    --stderr-contains "'source_points_path' is deprecated"
    --stderr-contains "use 'source_points_fname'"
)
lsdtt_set_tests_properties(rigorous_bankfull_width_legacy_points_key PROPERTIES
    FIXTURES_REQUIRED cross_sections_fixture)

if(LSDTT_ENABLE_PYTHON_TESTS)
    add_test(NAME fixture_marsh_topography_dem
        COMMAND ${Python3_EXECUTABLE}
            "${CMAKE_CURRENT_SOURCE_DIR}/scripts/make_marsh_topography_fixture.py"
            "${LSDTT_MARSH_TOPOGRAPHY_FIXTURE_DIR}")
    lsdtt_mark_test_runtime(python fixture_marsh_topography_dem)
else()
    lsdtt_mark_test_unregistered(python fixture_marsh_topography_dem)
endif()
lsdtt_set_tests_properties(fixture_marsh_topography_dem PROPERTIES
    LABELS "fixture;marsh"
    FIXTURES_SETUP marsh_topography_fixture)

add_lsdtt_cli_smoke_test(
    smoke_marsh_topography
    lsdtt-marsh-topography
    marsh_topography_smoke.param.in
    marsh_marsh_platform_mask.tif
    --read-path "${LSDTT_MARSH_TOPOGRAPHY_FIXTURE_DIR}"
    --expect marsh_marsh_platform_labels.tif
    --expect marsh_marsh_scarp_mask.tif
    --expect marsh_marsh_platforms.fgb
    --expect-vector-feature-min marsh_marsh_platforms.fgb:1
    --expect-vector-feature-min marsh_tidal_channel_polygons.fgb:1
    --expect marsh_tidal_channel_mask.tif
    --expect marsh_tidal_channel_skeleton.tif
    --expect marsh_tidal_upl.tif
    --expect marsh_tidal_sections.fgb
    --expect marsh_tidal_stations.fgb
    --expect marsh_tidal_stations.csv
    --stdout-contains "Building nodes and sections (8 segments)..."
    --stdout-contains "Computing cross-sections (35 stations)..."
)
lsdtt_set_tests_properties(smoke_marsh_topography PROPERTIES
    FIXTURES_REQUIRED marsh_topography_fixture)

add_lsdtt_cli_smoke_test(
    smoke_water_centreline
    lsdtt-water-centreline
    water_centreline_smoke.param.in
    wc_water_centreline.fgb
    --stdout-contains "Water areas: 2   Water lines: 5"
    --stdout-contains "treated as loch    : 1"
    --stdout-contains "Tributary connectors : 4"
    --stdout-contains "Pass-through lines    : 4"
)

# =============================================================================
# lsdtt-advection-diffusion
#
# No test until 2026-10-14. The CLI generates its own red-noise surface, so it
# needs no input data. With random_seed and n_threads fixed the outputs are
# byte-identical between runs and across thread counts (checked for 1, 4 and
# 0), so the surfaces are pinned by hash.
#
# smoke_advection_diffusion is also the fixture for the resume and initial-grid
# tests, which read its out/ directory. The resume test caught a real bug:
# a resumed run re-applied fill_initial_surface to the mid-run grid and
# diverged from an uninterrupted run by up to 4 m. Its hash is the
# uninterrupted 40-step run's, so the two must now agree exactly.
# =============================================================================
set(LSDTT_ADVECTION_DIFFUSION_RUN_DIR
    "${CMAKE_CURRENT_BINARY_DIR}/smoke_advection_diffusion/out")

add_lsdtt_cli_smoke_test(
    smoke_advection_diffusion
    lsdtt-advection-diffusion
    advection_diffusion_smoke.param.in
    ad_modelled_topography.tif
    --expect ad_initial_topography.tif
    --expect ad_model_iteration_10.tif
    --expect ad_model_iteration_20.tif
    --expect ad_initial_hillshade.tif
    --expect ad_initial_hillshade.png
    --expect ad_model_iteration_10_hillshade.png
    --expect ad_model_iteration_10_log10_drainage_area.tif
    --expect ad_model_iteration_20_log10_drainage_area.png
    --expect ad_modelled_hillshade.tif
    --expect ad_modelled_log10_drainage_area.png
    --expect-raster-size ad_modelled_topography.tif:50x40
    --expect-sha256 ad_initial_topography.tif:abc51334ceca974a70a49068686ac00a1c3234dedcd8766e5604af3295b8f846
    --expect-sha256 ad_modelled_topography.tif:bde15e44e80df6d1d6e5279fed20d53f6325c71061506ca2ecc2e728b8b8680b
)
lsdtt_set_tests_properties(smoke_advection_diffusion PROPERTIES
    FIXTURES_SETUP advection_diffusion_run)

add_lsdtt_cli_smoke_test(
    smoke_advection_diffusion_resume
    lsdtt-advection-diffusion
    advection_diffusion_resume.param.in
    ../../smoke_advection_diffusion/out/ad_model_iteration_40.tif
    --read-path "${LSDTT_ADVECTION_DIFFUSION_RUN_DIR}"
    --stdout-contains "Resuming from intermediate output at iteration 20"
    --expect ../../smoke_advection_diffusion/out/ad_model_iteration_30.tif
    --expect-sha256 ../../smoke_advection_diffusion/out/ad_model_iteration_40.tif:ce85d0bfb321c5ea2d1ed27b3e42b59825279131595e222feb2846972b3ed1f7
)
lsdtt_set_tests_properties(smoke_advection_diffusion_resume PROPERTIES
    FIXTURES_REQUIRED advection_diffusion_run)

add_lsdtt_cli_smoke_test(
    smoke_advection_diffusion_initial_grid
    lsdtt-advection-diffusion
    advection_diffusion_initial_grid.param.in
    ad_modelled_topography.tif
    --read-path "${LSDTT_ADVECTION_DIFFUSION_RUN_DIR}"
    --stdout-contains "Loading existing grid for a new run"
    --expect ad_initial_topography.tif
    --expect-raster-size ad_modelled_topography.tif:50x40
)
lsdtt_set_tests_properties(smoke_advection_diffusion_initial_grid PROPERTIES
    FIXTURES_REQUIRED advection_diffusion_run)

add_lsdtt_cli_smoke_test(
    smoke_advection_diffusion_d8
    lsdtt-advection-diffusion
    advection_diffusion_d8.param.in
    ad_d8_modelled_topography.tif
    --expect ad_d8_model_iteration_20.tif
    --expect-raster-size ad_d8_modelled_topography.tif:50x40
)

foreach(_ad_bad
        "kernel|drainage_area_moving_average_kernel_size_pixels must be a positive odd number"
        "boundary_conditions|boundary_conditions must contain four characters"
        "diagnostics_format|diagnostics_write_output_format must be cog, png, or both"
        "incision_model|fluvial_incision_model must be d8 or dinf")
    string(REPLACE "|" ";" _ad_bad "${_ad_bad}")
    list(GET _ad_bad 0 _ad_case)
    list(GET _ad_bad 1 _ad_message)
    add_lsdtt_cli_smoke_test(
        failure_advection_diffusion_bad_${_ad_case}
        lsdtt-advection-diffusion
        advection_diffusion_bad_${_ad_case}.param.in
        ad_modelled_topography.tif
        --expect-failure
        --stderr-contains "FATAL: ${_ad_message}"
    )
endforeach()

# =============================================================================
# lsdtt-slope-area
#
# No test until 2026-10-14. Runs on the basin_stats fixture DEM (200 x 200, 30 m,
# UTM 58S). Outputs are identical for n_threads 1 and 4; the summary rows are
# pinned as text rather than by hash. The clip test also guards a fix: a bare
# basin_mask_fname used to be opened relative to the working directory rather
# than read_path.
# =============================================================================
if(LSDTT_ENABLE_PYTHON_TESTS)
    lsdtt_add_harness_test(fixture_slope_area_basin_mask
            --exe "$<TARGET_FILE:lsdtt-basin-extraction>"
            --cli-src-dir "${CMAKE_SOURCE_DIR}/cli"
            --param-template "${LSDTT_TEST_PARAM_DIR}/slope_area_basin_mask.param.in"
            --read-path "${LSDTT_BASIN_STATS_FIXTURE_DIR}/out"
            --work-dir "${CMAKE_CURRENT_BINARY_DIR}/fixture_slope_area_basin_mask"
            --expect ../../basin_stats_fixture/out/sa_basins_junction_basins.tif)
else()
    lsdtt_mark_test_unregistered(python fixture_slope_area_basin_mask)
endif()
lsdtt_set_tests_properties(fixture_slope_area_basin_mask PROPERTIES
    LABELS "fixture;slope_area"
    FIXTURES_SETUP slope_area_mask
    FIXTURES_REQUIRED basin_stats_fixture)

# tiny_dem.asc (5 x 5, UTM 30N, nowhere near the fixture DEM) next to the
# basin_stats DEM: a mismatched grid for lsdtt-slope-area and a non-overlapping
# target for lsdtt-outlet-matching.
add_test(NAME fixture_basin_stats_tiny_dem
    COMMAND ${CMAKE_COMMAND} -E copy
        "${LSDTT_TEST_DATA_DIR}/tiny_dem.asc"
        "${LSDTT_TEST_DATA_DIR}/tiny_dem.prj"
        "${LSDTT_BASIN_STATS_FIXTURE_DIR}/out/")
lsdtt_set_tests_properties(fixture_basin_stats_tiny_dem PROPERTIES
    LABELS "fixture;basin_stats"
    FIXTURES_SETUP basin_stats_tiny_dem
    FIXTURES_REQUIRED basin_stats_fixture)

add_lsdtt_cli_smoke_test(
    smoke_slope_area
    lsdtt-slope-area
    slope_area_smoke.param.in
    sa_maximum_d50_slope_area.csv
    --read-path "${LSDTT_BASIN_STATS_FIXTURE_DIR}/out"
    --expect sa_lambda_90m.csv
    --expect sa_lambda_150m.csv
    --expect sa_lambda_90m_slope.tif
    --expect sa_lambda_150m_slope.tif
    --expect sa_dinf_drainage_area.tif
    --expect-absent sa_d8_flow_direction.tif
    --expect-file-contains "sa_maximum_d50_slope_area.csv:90,3962.762065"
    --expect-file-contains "sa_maximum_d50_slope_area.csv:150,3952.876248"
    --expect-valid-percent sa_lambda_90m_slope.tif:96.04
)
lsdtt_set_tests_properties(smoke_slope_area PROPERTIES
    FIXTURES_REQUIRED basin_stats_fixture)

add_lsdtt_cli_smoke_test(
    smoke_slope_area_d8
    lsdtt-slope-area
    slope_area_d8.param.in
    sa_maximum_d50_slope_area.csv
    --read-path "${LSDTT_BASIN_STATS_FIXTURE_DIR}/out"
    --expect sa_d8_drainage_area.tif
    --expect sa_d8_flow_direction.tif
    --expect-absent sa_dinf_drainage_area.tif
    --expect-file-contains "sa_maximum_d50_slope_area.csv:90,3965.578391"
)
lsdtt_set_tests_properties(smoke_slope_area_d8 PROPERTIES
    FIXTURES_REQUIRED basin_stats_fixture)

add_lsdtt_cli_smoke_test(
    smoke_slope_area_no_flowinfo
    lsdtt-slope-area
    slope_area_no_flowinfo.param.in
    sa_maximum_d50_slope_area.csv
    --read-path "${LSDTT_BASIN_STATS_FIXTURE_DIR}/out"
    --expect-absent sa_dinf_drainage_area.tif
    --expect-file-contains "sa_maximum_d50_slope_area.csv:90,3962.762065"
)
lsdtt_set_tests_properties(smoke_slope_area_no_flowinfo PROPERTIES
    FIXTURES_REQUIRED basin_stats_fixture)

add_lsdtt_cli_smoke_test(
    smoke_slope_area_clip
    lsdtt-slope-area
    slope_area_clip.param.in
    sa_maximum_d50_slope_area.csv
    --read-path "${LSDTT_BASIN_STATS_FIXTURE_DIR}/out"
    --expect-valid-percent sa_lambda_90m_slope.tif:53.79
    --expect-file-contains "sa_maximum_d50_slope_area.csv:90,3906.624295"
)
lsdtt_set_tests_properties(smoke_slope_area_clip PROPERTIES
    FIXTURES_REQUIRED "basin_stats_fixture;slope_area_mask")

foreach(_sa_bad
        "bad_flow_algorithm|flow_algorithm must be dinf or d8"
        "bad_clip_mode|basin_clip_mode must be none, mask, or polygon"
        "empty_lambdas|polyfit_lambda_metres_list cannot be empty"
        "mask_grid_mismatch|Basin clip does not have exactly the same grid settings as the DEM")
    string(REPLACE "|" ";" _sa_bad "${_sa_bad}")
    list(GET _sa_bad 0 _sa_case)
    list(GET _sa_bad 1 _sa_message)
    add_lsdtt_cli_smoke_test(
        failure_slope_area_${_sa_case}
        lsdtt-slope-area
        slope_area_${_sa_case}.param.in
        sa_maximum_d50_slope_area.csv
        --read-path "${LSDTT_BASIN_STATS_FIXTURE_DIR}/out"
        --expect-failure
        --stderr-contains "Error: ${_sa_message}"
    )
    lsdtt_set_tests_properties(failure_slope_area_${_sa_case} PROPERTIES
        FIXTURES_REQUIRED "basin_stats_fixture;basin_stats_tiny_dem")
endforeach()

add_lsdtt_cli_smoke_test(
    failure_slope_area_geographic_dem
    lsdtt-slope-area
    slope_area_geographic_dem.param.in
    sa_maximum_d50_slope_area.csv
    --expect-failure
    --stderr-contains "only works on projected data"
)

# =============================================================================
# lsdtt-outlet-matching
#
# No test until 2026-10-14. Reference = the basin_stats fixture DEM (30 m) and
# its 29 outlets; target = the same DEM (identity) or resampled to 60 m. Rows
# are pinned by their integer columns (target row, col, flat, areas), not the
# 12-digit coordinates. The CLI runs with its working directory outside
# write_path, which guards the fix that made default output names honour
# write_path.
# =============================================================================
if(LSDTT_ENABLE_PYTHON_TESTS)
    lsdtt_add_harness_test(fixture_outlet_matching_target
            --exe "$<TARGET_FILE:lsdtt-raster-preprocessing>"
            --cli-src-dir "${CMAKE_SOURCE_DIR}/cli"
            --param-template "${LSDTT_TEST_PARAM_DIR}/outlet_matching_target_dem.param.in"
            --read-path "${LSDTT_BASIN_STATS_FIXTURE_DIR}/out"
            --work-dir "${CMAKE_CURRENT_BINARY_DIR}/fixture_outlet_matching_target"
            --expect ../../basin_stats_fixture/out/om_target_utm.tif)
else()
    lsdtt_mark_test_unregistered(python fixture_outlet_matching_target)
endif()
lsdtt_set_tests_properties(fixture_outlet_matching_target PROPERTIES
    LABELS "fixture;outlet_matching"
    FIXTURES_SETUP outlet_matching_target
    FIXTURES_REQUIRED basin_stats_fixture)

add_lsdtt_cli_smoke_test(
    smoke_outlet_matching_identity
    lsdtt-outlet-matching
    outlet_matching_identity.param.in
    om_matched_outlets.csv
    --read-path "${LSDTT_BASIN_STATS_FIXTURE_DIR}/out"
    --expect om_outlet_matching_diagnostics.csv
    --expect-file-not-contains "om_matched_outlets.csv:no_match"
    --expect-file-not-contains "om_matched_outlets.csv:invalid_reference_point"
    --expect-file-contains "om_matched_outlets.csv:,30,30,6030,10800,10800,0,0,"
    --expect-file-contains "om_matched_outlets.csv:,65,30,13030,5400,5400,0,0,"
    --expect-file-contains "om_matched_outlets.csv:,89,18,17818,325800,299700,26100,"
)
lsdtt_set_tests_properties(smoke_outlet_matching_identity PROPERTIES
    FIXTURES_REQUIRED basin_stats_fixture)

add_lsdtt_cli_smoke_test(
    smoke_outlet_matching_resampled
    lsdtt-outlet-matching
    outlet_matching_resampled.param.in
    om_matched_outlets.csv
    --read-path "${LSDTT_BASIN_STATS_FIXTURE_DIR}/out"
    --expect-file-not-contains "om_matched_outlets.csv:no_match"
    --expect-file-contains "om_matched_outlets.csv:,14,16,1416,10800,10800,0,0,"
    --expect-file-contains "om_matched_outlets.csv:,32,15,3215,3600,5400,-1800,"
    --expect-file-contains "om_matched_outlets.csv:,42,7,4207,291600,299700,-8100,"
)
lsdtt_set_tests_properties(smoke_outlet_matching_resampled PROPERTIES
    FIXTURES_REQUIRED "basin_stats_fixture;outlet_matching_target")

add_lsdtt_cli_smoke_test(
    failure_outlet_matching_no_overlap
    lsdtt-outlet-matching
    outlet_matching_no_overlap.param.in
    om_matched_outlets.csv
    --read-path "${LSDTT_BASIN_STATS_FIXTURE_DIR}/out"
    --expect-failure
    --stderr-contains "FATAL: Reference and target DEM extents do not overlap after CRS transformation."
)
lsdtt_set_tests_properties(failure_outlet_matching_no_overlap PROPERTIES
    FIXTURES_REQUIRED "basin_stats_fixture;basin_stats_tiny_dem")

add_lsdtt_cli_smoke_test(
    failure_outlet_matching_bad_radius
    lsdtt-outlet-matching
    outlet_matching_bad_radius.param.in
    om_matched_outlets.csv
    --read-path "${LSDTT_BASIN_STATS_FIXTURE_DIR}/out"
    --expect-failure
    --stderr-contains "local_search_radius_metres must be greater than zero."
)
lsdtt_set_tests_properties(failure_outlet_matching_bad_radius PROPERTIES
    FIXTURES_REQUIRED basin_stats_fixture)

# raster_is_filled: true on an unfilled DEM: routing fails, and the error must
# point at raster_is_filled rather than leave the user with a routing message.
add_lsdtt_cli_smoke_test(
    failure_outlet_matching_unfilled_dem
    lsdtt-outlet-matching
    outlet_matching_unfilled_dem.param.in
    om_matched_outlets.csv
    --read-path "${LSDTT_BASIN_STATS_FIXTURE_DIR}/out"
    --expect-failure
    --stderr-contains "flow routing failed on the reference DEM"
    --stderr-contains "Set raster_is_filled: false to fill both DEMs first."
)
lsdtt_set_tests_properties(failure_outlet_matching_unfilled_dem PROPERTIES
    FIXTURES_REQUIRED basin_stats_fixture)

# An outlet whose only candidate is edge-influenced is a no_match, reported with
# its reason; the run still succeeds for the others.
add_lsdtt_cli_smoke_test(
    smoke_outlet_matching_edge_only
    lsdtt-outlet-matching
    outlet_matching_edge_only.param.in
    om_matched_outlets.csv
    --read-path "${LSDTT_BASIN_STATS_FIXTURE_DIR}/out"
    --stdout-contains "Outlet 27: no_match. The only candidate drains area from the target DEM's edge"
    --stdout-contains "Outlets matched: 28, no_match: 1, invalid_reference_point: 0"
)
lsdtt_set_tests_properties(smoke_outlet_matching_edge_only PROPERTIES
    FIXTURES_REQUIRED basin_stats_fixture)

# =============================================================================
# lsdtt-water-line-repair
#
# Line-end joining by distance (lsdtt_point_snapper, river-lines plan D9):
# exact cases, including the two the old grid-rounding key got wrong.
add_executable(lsdtt-test-point-snapper cpp/test_point_snapper.cpp)
target_link_libraries(lsdtt-test-point-snapper PRIVATE lsdtt3_cpu)
add_test(NAME unit_water_line_repair_point_snapper COMMAND lsdtt-test-point-snapper)
lsdtt_set_tests_properties(unit_water_line_repair_point_snapper PROPERTIES LABELS "water")

# No test until 2026-10-14. Pure vector: tests/data/water_line_repair_fixture
# lays out one case per bridging rule (see water_line_repair_smoke.param.in).
# The bridge count in stdout is the assertion; each variant changes one rule
# and the count moves accordingly.
# =============================================================================
add_lsdtt_cli_smoke_test(
    smoke_water_line_repair
    lsdtt-water-line-repair
    water_line_repair_smoke.param.in
    wlr_water_line_repaired.fgb
    --expect wlr_repair_tips.fgb
    --stdout-contains "Dangling tips        : 24"
    --stdout-contains "Bridges accepted     : 3"
    --stdout-contains "Output features      : 15"
    --expect-vector-feature-min wlr_water_line_repaired.fgb:15
    --expect-vector-feature-min wlr_repair_tips.fgb:24
)

foreach(_wlr_case "mutual_best_on|3" "mutual_best_off|4" "wide_road_gap|4"
                  "separate_roads|3")
    string(REPLACE "|" ";" _wlr_case "${_wlr_case}")
    list(GET _wlr_case 0 _wlr_name)
    list(GET _wlr_case 1 _wlr_bridges)
    add_lsdtt_cli_smoke_test(
        smoke_water_line_repair_${_wlr_name}
        lsdtt-water-line-repair
        water_line_repair_${_wlr_name}.param.in
        wlr_water_line_repaired.fgb
        --stdout-contains "Bridges accepted     : ${_wlr_bridges}"
    )
endforeach()

add_lsdtt_cli_smoke_test(
    failure_water_line_repair_missing_layer
    lsdtt-water-line-repair
    water_line_repair_missing_layer.param.in
    wlr_water_line_repaired.fgb
    --expect-failure
    --stderr-contains "FATAL: water_layer and road_layer are required."
)

# =============================================================================
# lsdtt-swath-profile
#
# Shipped in the core profile with no test at all (it predates the rule in
# tools/cli_profiles.txt). Tested 2026-10-14 against analytic answers.
# =============================================================================
add_lsdtt_cli_smoke_test(
    smoke_swath_profile_plane
    lsdtt-swath-profile
    swath_profile_plane.param.in
    sw_swath_profile.csv
    --expect sw_swath_along_axis.tif
    --expect sw_swath_perp_distance.tif
    --expect sw_swath_signed_perp.tif
    --expect-file-contains "sw_swath_profile.csv:50,200,5,0.5,"
    --expect-file-contains "sw_swath_profile.csv:250,200,25,20.5,"
    --expect-file-contains "sw_swath_profile.csv:450,200,45,40.5,"
)

# =============================================================================
# lsdtt-basin-stats
#
# This CLI had NO tests and NO fixtures until 2026-09-10. It needed both before
# the C1 extraction moved ~840 lines (the record struct, the per-basin
# accumulation loop, and the two writers) out of its main() into
# src/cpu/lsdtt_basin_stats_cpu.
#
# The DEM is GENERATED rather than shipped, because tests/data/dem_cache is
# gitignored and nothing in the suite may depend on a fetched artefact. A fixed
# seed and the fixed synthetic georeference make it byte-reproducible, which is
# what lets the smoke test below carry a sha256 ORACLE on the stats CSV --
# a refactor tripwire, re-freezable by a deliberate act, in the sense used in
# docs/developer_notes/testing_and_benchmarking.md.
# =============================================================================
if(LSDTT_ENABLE_PYTHON_TESTS)
    lsdtt_add_harness_test(fixture_basin_stats_dem
            --exe "$<TARGET_FILE:lsdtt-synthetic-surfaces>"
            --cli-src-dir "${CMAKE_SOURCE_DIR}/cli"
            --param-template "${LSDTT_TEST_PARAM_DIR}/basin_stats_fixture_dem.param.in"
            --read-path "${LSDTT_TEST_DATA_DIR}"
            --work-dir "${LSDTT_BASIN_STATS_FIXTURE_DIR}"
            --expect basin_stats_diamond_square.tif)
else()
    lsdtt_mark_test_unregistered(python fixture_basin_stats_dem)
endif()
lsdtt_set_tests_properties(fixture_basin_stats_dem PROPERTIES
    LABELS "fixture;basin_stats"
    FIXTURES_SETUP basin_stats_fixture)

# The harness wipes its work dir, so the outlet CSV is copied in AFTER the DEM
# is generated rather than at configure time. Both are FIXTURES_SETUP for the
# same fixture; DEPENDS fixes the order between them.
add_test(NAME fixture_basin_stats_outlets
    COMMAND ${CMAKE_COMMAND} -E copy
        "${LSDTT_TEST_DATA_DIR}/basin_stats_outlets.csv"
        "${LSDTT_BASIN_STATS_FIXTURE_DIR}/out/basin_stats_outlets.csv")
lsdtt_set_tests_properties(fixture_basin_stats_outlets PROPERTIES
    LABELS "fixture;basin_stats"
    FIXTURES_SETUP basin_stats_fixture
    DEPENDS fixture_basin_stats_dem)

if(LSDTT_ENABLE_PYTHON_TESTS)
    lsdtt_add_harness_test(smoke_basin_stats
            --exe "$<TARGET_FILE:lsdtt-basin-stats>"
            --cli-src-dir "${CMAKE_SOURCE_DIR}/cli"
            --param-template "${LSDTT_TEST_PARAM_DIR}/basin_stats_smoke.param.in"
            --read-path "${LSDTT_BASIN_STATS_FIXTURE_DIR}/out"
            --work-dir "${CMAKE_CURRENT_BINARY_DIR}/smoke_basin_stats"
            --expect bs_basin_stats.csv
            --expect bs_basin_polygons.fgb
            --expect bs_junction_area_ratios.csv
            --expect-vector-feature-min bs_basin_polygons.fgb:10
            # ORACLE, not a numeric standard: a tripwire for the C1 extraction of
            # the record, the per-basin accumulation and the two writers out of
            # main(). It says "the output did not change", nothing about whether
            # the output is right. Re-freeze deliberately if a change is intended.
            # Re-frozen 2026-10-25: the table is written at 17 significant
            # digits (was the stream default of 6, ~5 m in latitude/longitude),
            # and as CSV only because the template asks for it. Check made then:
            # the new file reprinted at 6 significant digits hashes to the old
            # oracle (6b812973...), so only the precision changed.
            --expect-sha256 bs_basin_stats.csv:3917286425c4931152f04b9dbc97f52740a6d929e7d88e5eb5701f9170fbb3ea
            # Canonical-default regression: the parameter is deliberately absent
            # from the fixture, so this pins the shared default and its behaviour.
            --stdout-contains "threshold = 1000 px"
            --stderr-contains "'basin_outlet_csv' is deprecated and 'basin_outlet_fname' was also supplied")
else()
    lsdtt_mark_test_unregistered(python smoke_basin_stats)
endif()
lsdtt_set_tests_properties(smoke_basin_stats PROPERTIES
    LABELS "cli;basin_stats"
    FIXTURES_REQUIRED basin_stats_fixture)

add_lsdtt_cli_smoke_test(
    smoke_mountain_front
    lsdtt-mountain-front
    mountain_front_smoke.param.in
    mf_relief.tif
    --expect mf_mountain_front_outlets.csv
    --stdout-contains "2 outlet(s) written."
    --expect-sha256 mf_mountain_front_outlets.csv:5367371b46baac70a366d9bf1648293f5aaf665aa27ebfcc41f5341661780c56
)

# mountain_front_smoke under the names used before 2026-10-26 (a window side in
# cells, adaptive_relief_threshold, relief_threshold_metres): the same outlets.
add_lsdtt_cli_smoke_test(
    rigorous_mountain_front_legacy_keys
    lsdtt-mountain-front
    mountain_front_legacy_keys.param.in
    mf_relief.tif
    --expect-sha256 mf_mountain_front_outlets.csv:5367371b46baac70a366d9bf1648293f5aaf665aa27ebfcc41f5341661780c56
    --stderr-contains "'relief_window_side_pixels' is deprecated"
    --stderr-contains "use 'relief_window_radius_metres'"
    --stderr-contains "'adaptive_relief_threshold' is deprecated"
    --stderr-contains "use 'relief_threshold_algorithm_choice'"
    --stderr-contains "'relief_threshold_metres' is deprecated"
    --stderr-contains "use 'max_relief_metres'"
    --expect-file-contains "mf_ingested_params.param:relief_window_radius_metres: 1"
    --expect-file-contains "mf_ingested_params.param:relief_threshold_algorithm_choice: fixed"
)

# The CDF-divergence relief bound, pinned 2026-10-26 before its code moved into
# the shared threshold-selection utility (lsdtt_threshold_selection). The smoke
# test above uses a fixed bound, so nothing else exercises this path.
add_lsdtt_cli_smoke_test(
    rigorous_mountain_front_cdf_divergence
    lsdtt-mountain-front
    mountain_front_cdf_divergence.param.in
    mf_relief.tif
    --stdout-contains "Best-fit distribution: normal"
    --expect-raster-checksum mf_mountain_mask.tif:21
    --expect-sha256 mf_relief_cdf.csv:ded637956d805c709491cfb3ac9e22eae6729a193e4b6b05ec6aa7ca4b843c43
    --expect-sha256 mf_mountain_front_outlets.csv:5367371b46baac70a366d9bf1648293f5aaf665aa27ebfcc41f5341661780c56
)

# The Wiener source extraction's QQ table: the only pin on compute_qq, which
# the threshold-selection utility extends with percentile arguments. Pinned
# 2026-10-26 before that change.
add_lsdtt_cli_smoke_test(
    rigorous_channel_extraction_wiener_qq
    lsdtt-channel-extraction
    channel_extraction_wiener_qq.param.in
    wq_wiener_qq.csv
    --read-path "${LSDTT_FLUVIAL_SURFACES_FIXTURE_DIR}"
    --stdout-contains "Number of source pixels: 1"
    --expect-sha256 wq_wiener_qq.csv:6ceabf0a02d502387331ab90dcfce254e68db44b4c3b446dcacc63030356ab35
)
lsdtt_set_tests_properties(rigorous_channel_extraction_wiener_qq PROPERTIES
    FIXTURES_REQUIRED fluvial_surfaces_fixture)

add_lsdtt_cli_smoke_test(
    smoke_feature_finder
    lsdtt-feature-finder
    feature_finder_smoke.param.in
    ff_channel_candidate_score.tif
)

# feature_type find_risers on the dissected staircase (fluvial-surface methods
# plan, phase 6). With astar the six risers come back bridged across the graded
# gaps, plus twelve short segments on the gaps' up- and down-valley edges,
# which are real across-valley scarps in this fixture (step 8.7-8.9 m).
add_lsdtt_cli_smoke_test(
    rigorous_feature_finder_risers_staircase
    lsdtt-feature-finder
    feature_finder_risers_staircase.param.in
    ffr_riser_centreline.tif
    --read-path "${LSDTT_FLUVIAL_SURFACES_FIXTURE_DIR}"
    --expect ffr_riser_segments.fgb
    --expect ffr_riser_endpoints.fgb
    --expect ffr_riser_gap_closing.tif
    --expect ffr_riser_edge_field.tif
    --expect-raster-values ffr_riser_centreline.tif:0,1
    --expect-vector-feature-min ffr_riser_segments.fgb:18
    --stdout-contains "Risers kept          : 18"
    --expect-file-contains "ffr_ingested_params.param:riser_gap_closing_algorithm_choice: astar"
)
lsdtt_set_tests_properties(rigorous_feature_finder_risers_staircase PROPERTIES
    FIXTURES_REQUIRED fluvial_surfaces_fixture)

# mapped_risers_fname on the same staircase: the six tread boundaries as
# mapped lines, scored in the northern half (staircase_mapped_extent_north).
# The finder's line sits 2-3 cells (20-30 m) down onto the lower tread, at the
# foot of the riser where the 100 m polyfit puts the curvature peak, and the
# boundary line burns on one side of the step. The unmatched detections are
# the graded gaps' real up- and down-valley scarps, which the map does not
# draw. Pinned 2026-10-31 for the mean_std threshold (k = 0.5), the new
# default: 788 detected cells (905 under the 8-bit gaussian threshold);
# at 20 m precision 307/788, recall 307/600 (gaussian: 297/905, 297/600); at
# 30 m precision 423/788, recall 444/600 (gaussian: 513/905, 537/600). The
# narrower binary puts the skeleton nearer the curvature peak, further from
# the mapped line on the 30 m flank.
add_lsdtt_cli_smoke_test(
    rigorous_feature_finder_risers_mapped
    lsdtt-feature-finder
    feature_finder_risers_mapped.param.in
    ffrm_mapped_riser_scores.csv
    --read-path "${LSDTT_FLUVIAL_SURFACES_FIXTURE_DIR}"
    --expect ffrm_mapped_risers.tif
    --expect-raster-values ffrm_mapped_risers.tif:0,1
    --stdout-contains "Mapped risers (extent, 20 m): precision 0.389594 (307/788 cells), recall 0.511667 (307/600 cells)"
    --expect-file-contains "ffrm_mapped_riser_scores.csv:extent,20,10,30000,788,307,600,307,"
)
lsdtt_set_tests_properties(rigorous_feature_finder_risers_mapped PROPERTIES
    FIXTURES_REQUIRED fluvial_surfaces_fixture)

add_lsdtt_cli_smoke_test(
    rigorous_feature_finder_risers_mapped_30m
    lsdtt-feature-finder
    feature_finder_risers_mapped_30m.param.in
    ffrm30_mapped_riser_scores.csv
    --read-path "${LSDTT_FLUVIAL_SURFACES_FIXTURE_DIR}"
    --stdout-contains "Mapped risers (extent, 30 m): precision 0.536802 (423/788 cells), recall 0.74 (444/600 cells)"
)
lsdtt_set_tests_properties(rigorous_feature_finder_risers_mapped_30m PROPERTIES
    FIXTURES_REQUIRED fluvial_surfaces_fixture)

add_lsdtt_cli_smoke_test(
    rigorous_feature_finder_flow_discrimination_canonical
    lsdtt-feature-finder
    feature_finder_flow_discrimination_canonical.param.in
    ff_flow_canonical_flow_channel_mask.tif
    --expect-raster-checksum ff_flow_canonical_flow_channel_mask.tif:15
    --expect-raster-values ff_flow_canonical_flow_channel_mask.tif:0,1
    --stdout-contains "Accum. threshold     : 3 px"
    --stdout-contains "DEM is filled        : yes"
    --expect-file-contains "ff_flow_canonical_ingested_params.param:threshold_contributing_pixels: 3"
)

# Phase 4 readiness: reach the flow-discrimination engine with the legacy keys
# explicitly set. The header assertions distinguish both values from their old
# defaults; the mask oracle proves the threshold reached D8 channel selection.
add_lsdtt_cli_smoke_test(
    rigorous_feature_finder_flow_discrimination_legacy_keys
    lsdtt-feature-finder
    feature_finder_flow_discrimination.param.in
    ff_flow_flow_channel_mask.tif
    --expect ff_flow_flow_contributing_pixels.tif
    --expect-raster-checksum ff_flow_flow_channel_mask.tif:15
    --expect-raster-values ff_flow_flow_channel_mask.tif:0,1
    --stdout-contains "--- Flow-pathway road/channel discrimination ---"
    --stdout-contains "Accum. threshold     : 3 px"
    --stdout-contains "DEM is filled        : yes"
    --stderr-contains "'flow_dem_is_filled' is deprecated"
    --stderr-contains "use 'raster_is_filled'"
    --stderr-contains "'flow_accumulation_threshold_pixels' is deprecated"
    --stderr-contains "use 'threshold_contributing_pixels'"
    --expect-file-contains "ff_flow_ingested_params.param:threshold_contributing_pixels: 3"
)

# Drainage-area exceedance scaling and the Tokunaga summary: no test reached
# either before 2026-10-07. See the parameter file for the reported values.
add_lsdtt_cli_smoke_test(
    smoke_channel_network_analysis_area_scaling
    lsdtt-channel-network-analysis
    channel_network_analysis_area_scaling.param.in
    cnas_area_scaling_summary.csv
    --read-path "${LSDTT_STREAM_BURNING_FIXTURE_DIR}/out"
    --expect cnas_area_scaling_curve.csv
    --expect cnas_tokunaga_summary_fit.csv
    --expect cnas_tokunaga_summary_pair_counts.csv
    --stdout-contains "--- Drainage-area exceedance scaling ---"
    --stdout-contains "beta (MLE)                  : 0.6195"
)
lsdtt_set_tests_properties(smoke_channel_network_analysis_area_scaling PROPERTIES
    FIXTURES_REQUIRED stream_burning_fixture)

# Phase 4 readiness: this old spelling remains explicit after aliases land. The
# diagnostics and Strahler/Tokunaga CSV oracle prove that threshold-defined
# sources reach channel-network construction and a real network analysis.
add_lsdtt_cli_smoke_test(
    rigorous_channel_network_analysis_legacy_threshold
    lsdtt-channel-network-analysis
    channel_network_analysis_legacy_threshold.param.in
    cna_legacy_tokunaga_channel_data.csv
    --stdout-contains "Contributing-pixel threshold: 3"
    --stdout-contains "Number of source pixels: 5"
    --stdout-contains "--- Channel network construction ---"
    --stdout-contains "Channel cells : 15"
    --stdout-contains "Junctions     : 9"
    --stdout-contains "--- Strahler / Tokunaga analysis ---"
    --stdout-contains "Strahler channels           : 6"
    --stderr-contains "'source_pixels_threshold' is deprecated"
    --stderr-contains "use 'threshold_contributing_pixels'"
    --expect-sha256 cna_legacy_tokunaga_channel_data.csv:dd48f4ace30af5b591dc788c6930bb084da2f8274d6a5403e718f8e65b04d923
)

# =============================================================================
# Channel-link bearing and junction-angle statistics
# (lsdtt_channel_bearing_statistics_cpu, lsdtt_channel_network_grouping_cpu,
#  lsdtt-channel-network-analysis write_channel_bearing_statistics /
#  write_junction_angle_statistics)
#
# The driver modes are native and run in every configure. The first one holds
# the analytic circular two-sample cases that the retired lsdtt-bearing-comparison
# CLI tests pinned (same 8 bearings, same answers), now with no Python needed.
# The trellis and grouping modes build synthetic networks in memory: a trellis
# whose trunks run N-S and whose tributaries run E-W, so the expected bearings,
# axial/directional contrast and junction angles follow from the geometry.
# =============================================================================
add_executable(lsdtt-test-channel-bearing-statistics cpp/test_channel_bearing_statistics.cpp)
target_link_libraries(lsdtt-test-channel-bearing-statistics PRIVATE lsdtt3_cpu)

foreach(_cbs_mode
        "circular|--circular-two-sample"
        "summaries|--group-summaries"
        "comparisons|--comparisons"
        "junction_angles|--junction-angle-stats"
        "trellis|--trellis-network"
        "grouping|--grouping-keys"
        "isotropic_null|--isotropic-null")
    string(REPLACE "|" ";" _cbs_mode "${_cbs_mode}")
    list(GET _cbs_mode 0 _cbs_name)
    list(GET _cbs_mode 1 _cbs_flag)
    add_test(NAME unit_channel_bearing_statistics_${_cbs_name}
             COMMAND lsdtt-test-channel-bearing-statistics ${_cbs_flag})
    lsdtt_set_tests_properties(unit_channel_bearing_statistics_${_cbs_name} PROPERTIES
        LABELS "network")
endforeach()

# CLI fixture: the trellis DEM and its transpose, a north/south region code
# raster on the trellis grid, an off-grid region raster with no CRS, and the
# committed reference CSV copied beside them.
set(LSDTT_CHANNEL_BEARING_STATISTICS_FIXTURE_DIR
    "${CMAKE_CURRENT_BINARY_DIR}/channel_bearing_statistics_fixture")
add_test(NAME fixture_channel_bearing_statistics
         COMMAND lsdtt-test-channel-bearing-statistics --write-fixtures
                 "${LSDTT_CHANNEL_BEARING_STATISTICS_FIXTURE_DIR}"
                 "${LSDTT_TEST_DATA_DIR}/channel_bearing_statistics/ns_trunk_bearings.csv")
lsdtt_set_tests_properties(fixture_channel_bearing_statistics PROPERTIES
    LABELS "fixture;network"
    FIXTURES_SETUP channel_bearing_statistics)

function(lsdtt_add_channel_bearing_cli_test TEST_NAME TEMPLATE EXPECTED)
    add_lsdtt_cli_smoke_test(${TEST_NAME} lsdtt-channel-network-analysis
        ${TEMPLATE} ${EXPECTED}
        --read-path "${LSDTT_CHANNEL_BEARING_STATISTICS_FIXTURE_DIR}" ${ARGN})
    if(TEST ${TEST_NAME})
        lsdtt_set_tests_properties(${TEST_NAME} PROPERTIES
            FIXTURES_REQUIRED channel_bearing_statistics)
    endif()
endfunction()

lsdtt_add_channel_bearing_cli_test(
    smoke_channel_network_analysis_bearing_statistics_order
    channel_network_analysis_bearing_statistics_order.param.in
    cbso_channel_bearing_group_summary.csv
    --expect cbso_channel_bearing_comparisons.csv
    --expect cbso_channel_bearing_rose.csv
    --expect cbso_channel_bearing_group_membership.csv
    --expect cbso_junction_angle_group_summary.csv
    --expect cbso_junction_angle_comparisons.csv
    --stdout-contains "Bearing group strahler_order=2 [axial]: n="
    --stdout-contains "Junction-angle group strahler_order=2 [alpha]: n="
    --expect-file-contains "cbso_channel_bearing_group_summary.csv:measurand,group_key,group_value,stratum_key,stratum_value,angular_type,period_deg,n,n_blocks,total_length_m,circular_mean_deg_count,mean_resultant_length_count,circular_variance_count,circular_sd_deg_count,circular_mean_deg_length,mean_resultant_length_length,n_eff_length,rayleigh_z,rayleigh_p_iid_approx,fraction_on_d8_multiples,meets_min_group_size"
    --expect-file-contains "cbso_channel_bearing_comparisons.csv:measurand,stratum_key,stratum_value,angular_type,comparison_mode,group_a,group_b,n_a,n_b,n_blocks_a,n_blocks_b,circular_wasserstein_deg,rotation_invariant,optimal_rotation_deg,kuiper_V,kuiper_p_iid_approx,kuiper_p_holm,n_comparisons_in_family,permutation_p,note"
)

# The rotated trellis (trunks draining west) against a committed reference of
# N-S trunk bearings (180, order 2) read through the order filter: axial W1 90.
lsdtt_add_channel_bearing_cli_test(
    rigorous_channel_network_analysis_bearing_statistics_reference
    channel_network_analysis_bearing_statistics_reference.param.in
    cbsr_channel_bearing_comparisons.csv
    --stdout-contains "Reference bearings           : 4 from"
    --stdout-contains "Bearing comparison strahler_order=2 vs reference [axial]: n="
    --stdout-contains "W1=90.0 deg"
)

lsdtt_add_channel_bearing_cli_test(
    rigorous_channel_network_analysis_bearing_statistics_region
    channel_network_analysis_bearing_statistics_region.param.in
    cbsg_channel_bearing_group_summary.csv
    --stdout-contains "Bearing group region=1 (strahler_order=2) [axial]: n="
    --stdout-contains "Bearing group region=2 (strahler_order=2) [axial]: n="
    --stdout-contains "Bearing comparison region=1 vs 2 (strahler_order=2) [axial]:"
    --stdout-contains "W1=0.0 deg"
    --expect-file-contains "cbsg_channel_bearing_group_membership.csv:upstream_junction,downstream_junction,stream_order,total_length_m,bearing_deg,axial_bearing_deg,group_value,stratum_value,basin_id,block_id,region_value,region_purity"
)

lsdtt_add_channel_bearing_cli_test(
    failure_channel_network_analysis_bearing_stratify_equals_group
    channel_network_analysis_bearing_stratify_equals_group.param.in
    unused.csv --expect-exit-code 1
    --stderr-contains "FATAL: channel_bearing_stratify_by must differ from channel_bearing_group_by."
)
lsdtt_add_channel_bearing_cli_test(
    failure_channel_network_analysis_bearing_region_missing
    channel_network_analysis_bearing_region_missing.param.in
    unused.csv --expect-exit-code 1
    --stderr-contains "FATAL: channel_region_fname is required when channel_bearing_group_by or channel_bearing_stratify_by is region."
)
lsdtt_add_channel_bearing_cli_test(
    failure_channel_network_analysis_bearing_unknown_key
    channel_network_analysis_bearing_unknown_key.param.in
    unused.csv --expect-exit-code 1
    --stderr-contains "FATAL: channel_bearing_group_by must be one of: strahler_order, basin, region, lithology."
)
lsdtt_add_channel_bearing_cli_test(
    failure_channel_network_analysis_bearing_donor_pair
    channel_network_analysis_bearing_donor_pair.param.in
    unused.csv --expect-exit-code 1
    --stderr-contains "FATAL: donor_order_pair is only valid for junction_angle_group_by."
)
lsdtt_add_channel_bearing_cli_test(
    failure_channel_network_analysis_bearing_region_grid
    channel_network_analysis_bearing_region_grid.param.in
    unused.csv --expect-exit-code 1
    --stderr-contains "is not on the DEM grid"
)
lsdtt_add_channel_bearing_cli_test(
    failure_channel_network_analysis_bearing_reference_column
    channel_network_analysis_bearing_reference_column.param.in
    unused.csv --expect-exit-code 1
    --stderr-contains "Column 'nope' not found in"
)

# Phase 4 readiness: preserve the legacy spelling as a compatibility test after
# aliases land. Diagnostics pin threshold-defined channel initiation and the
# raster oracle pins the result of an actual sediment-routing run.
add_lsdtt_cli_smoke_test(
    rigorous_sediment_routing_legacy_threshold
    lsdtt-sediment-routing
    sediment_routing_legacy_threshold.param.in
    sediment_legacy_cumulative_deposition.tif
    --expect sediment_legacy_d50_al.tif
    --stdout-contains "Extracting channel sources (threshold = 3 px)..."
    --stdout-contains "Source pixels: 5"
    --stdout-contains "Building channel network..."
    --stdout-contains "Junctions (reaches): 9"
    --stdout-contains "Running sediment routing (2 timesteps, 3600 s each)..."
    --stderr-contains "'source_pixels_threshold' is deprecated"
    --stderr-contains "use 'threshold_contributing_pixels'"
    --expect-valid-percent sediment_legacy_d50_al.tif:60
    --expect-raster-checksum sediment_legacy_d50_al.tif:65386
)

add_lsdtt_cli_smoke_test(
    rigorous_channel_extraction_legacy_threshold
    lsdtt-channel-extraction
    channel_extraction_legacy_threshold.param.in
    channel_legacy_channel_network_lines.fgb
    --stdout-contains "Contributing-pixel threshold: 1"
    --stderr-contains "'source_pixels_threshold' is deprecated"
    --stderr-contains "use 'threshold_contributing_pixels'"
)

add_lsdtt_cli_smoke_test(
    rigorous_chi_analysis_legacy_threshold
    lsdtt-chi-analysis
    chi_analysis_legacy_threshold.param.in
    chi_legacy_chi_fit_nodes_mn_0p500.csv
    --stderr-contains "'source_pixels_threshold' is deprecated"
    --stderr-contains "use 'threshold_contributing_pixels'"
)

add_lsdtt_cli_smoke_test(
    smoke_hydraulic_routing
    lsdtt-graphflood
    hydraulic_routing_smoke.param.in
    tiny_hydraulic_water_depth.tif
    --expect tiny_hydraulic_water_surface.tif
)

# Active node states on nodata DEM cells: a FATAL message with the count and
# the first cell, not an uncaught exception from the solver's constructor.
add_lsdtt_cli_smoke_test(
    failure_hydraulic_routing_active_on_nodata
    lsdtt-graphflood
    hydraulic_routing_active_on_nodata.param.in
    tiny_hydraulic_nodata_water_depth.tif
    --expect-failure
    --stderr-contains "are active in hydraulic_node_state_fname (1 or 3) but nodata in the DEM"
)

add_lsdtt_cli_smoke_test(
    smoke_chi_analysis
    lsdtt-chi-analysis
    chi_analysis_smoke.param.in
    tiny_chi_fit_nodes_mn_0p500.csv
)

add_lsdtt_cli_smoke_test(
    smoke_raster_calculator_difference
    lsdtt-raster-calculator
    raster_calculator_difference.param.in
    tiny_and_tiny_dem_plus1_difference.tif
)

# DEM coregistration (Nuth and Kaab) had no test until 2026-10-07. The DEM pair
# is generated: see tests/params/raster_calculator_coregistration.param.in.
if(LSDTT_ENABLE_PYTHON_TESTS)
    add_test(NAME fixture_raster_calculator_coregistration_dems
        COMMAND ${Python3_EXECUTABLE}
            "${CMAKE_CURRENT_SOURCE_DIR}/scripts/make_dem_coregistration_fixture.py"
            "${LSDTT_DEM_COREGISTRATION_FIXTURE_DIR}")
    lsdtt_mark_test_runtime(python fixture_raster_calculator_coregistration_dems)
else()
    lsdtt_mark_test_unregistered(python fixture_raster_calculator_coregistration_dems)
endif()
lsdtt_set_tests_properties(fixture_raster_calculator_coregistration_dems PROPERTIES
    LABELS "fixture;raster"
    FIXTURES_SETUP dem_coregistration_fixture)

add_lsdtt_cli_smoke_test(
    smoke_raster_calculator_coregistration
    lsdtt-raster-calculator
    raster_calculator_coregistration.param.in
    coreg_and_coreg_shifted_coreg_difference.tif
    --read-path "${LSDTT_DEM_COREGISTRATION_FIXTURE_DIR}"
    --expect coreg_shifted_coregistered_to_coreg.tif
    --expect coreg_and_coreg_shifted_coreg_diagnostics.csv
    --stdout-contains "DEM coregistration offsets [m]: dx=4.01"
    --stdout-contains "DEM coregistration iterations: 3, converged: yes"
)
lsdtt_set_tests_properties(smoke_raster_calculator_coregistration PROPERTIES
    FIXTURES_REQUIRED dem_coregistration_fixture)

add_lsdtt_cli_smoke_test(
    failure_channel_extraction_invalid_threshold
    lsdtt-channel-extraction
    channel_extraction_invalid_threshold.param.in
    tiny_invalid_channel_network_lines.fgb
    --expect-failure
    --stderr-contains "FATAL: threshold_contributing_pixels must be > 0."
)

add_lsdtt_cli_smoke_test(
    failure_raster_calculator_missing_second
    lsdtt-raster-calculator
    raster_calculator_missing_second.param.in
    tiny_missing_second_difference.tif
    --expect-failure
    --stderr-contains "FATAL: second_raster_read_fname must be set."
)

add_lsdtt_cli_smoke_test(
    failure_divide_migration_missing_outlets
    lsdtt-divide-migration
    divide_migration_missing_outlets.param.in
    unused_expected_output
    --expect-failure
    --stderr-contains "FATAL: basin_outlet_fname must be set."
)

add_lsdtt_cli_smoke_test(
    failure_divide_migration_bad_integer_cluster_mode
    lsdtt-divide-migration
    divide_migration_bad_integer_cluster_mode.param.in
    unused_expected_output
    --expect-failure
    --stderr-contains "FATAL: outlet cluster value 'alpha' is not an integer"
)

# --- divide-migration end to end -------------------------------
#
# These run the whole pipeline on the synthetic three-valley DEM written by
# unit_divide_migration, so they need the fixture directory rather than
# tests/data as their read path. run_cli_smoke.py stores --read-path (it does
# not append), so the value passed here through ARGN overrides the default the
# function supplies.
#
# The assertions are on CONTENT, not exit status. The divide-segment count is
# the topology: three parallel valleys give exactly two divides (1|2 and 2|3)
# and never a 1|3, so a run that matched heads across non-adjacent units, or
# that ignored the cluster field, cannot produce these numbers.
add_lsdtt_cli_smoke_test(
    smoke_divide_migration_no_clusters
    lsdtt-divide-migration
    divide_migration_smoke.param.in
    divmig_nocluster_divides.fgb
    --read-path "${LSDTT_DIVIDE_MIGRATION_FIXTURE_DIR}"
    --expect divmig_nocluster_channel_heads.fgb
    --expect divmig_nocluster_divide_migration_metrics.csv
    --expect-vector-feature-min divmig_nocluster_divides.fgb:2
    --expect-vector-feature-min divmig_nocluster_channel_heads.fgb:69
    --expect-file-contains "divmig_nocluster_divide_migration_metrics.csv:head_id,outlet_id,unit_id,pair_id_nearest,head_x,head_y,head_elevation_m,chi,nearest_divide_x,nearest_divide_y,divide_elevation_m,distance_to_divide_m,relief_m,gradient_to_divide"
    --stdout-contains "mapped units     : 3"
    --stdout-contains "selected basins  : 3"
    --stdout-contains "divide segments  : 2"
    --stdout-contains "channel heads    : 69"
    --stdout-contains "heads with divide: 69"
    --stdout-contains "heads without    : 0"
)
lsdtt_set_tests_properties(smoke_divide_migration_no_clusters PROPERTIES
    LABELS "cli;divide_migration"
    FIXTURES_REQUIRED divide_migration_fixture)

add_lsdtt_cli_smoke_test(
    rigorous_divide_migration_legacy_threshold
    lsdtt-divide-migration
    divide_migration_legacy_threshold.param.in
    divmig_legacy_divides.fgb
    --read-path "${LSDTT_DIVIDE_MIGRATION_FIXTURE_DIR}"
    --expect divmig_legacy_channel_heads.fgb
    --stdout-contains "channel heads    : 69"
    --stderr-contains "'source_pixels_threshold' is deprecated"
    --stderr-contains "use 'threshold_contributing_pixels'"
    --stderr-contains "'basin_outlet_file' is deprecated"
    --stderr-contains "use 'basin_outlet_fname'"
)
lsdtt_set_tests_properties(rigorous_divide_migration_legacy_threshold PROPERTIES
    FIXTURES_REQUIRED divide_migration_fixture)

# Grouping outlets 1 and 2 into one unit makes the 1-2 divide internal, so the
# two-divide topology above must collapse to exactly one segment.
add_lsdtt_cli_smoke_test(
    smoke_divide_migration_clusters
    lsdtt-divide-migration
    divide_migration_clusters_smoke.param.in
    divmig_cluster_divides.fgb
    --read-path "${LSDTT_DIVIDE_MIGRATION_FIXTURE_DIR}"
    --expect divmig_cluster_channel_heads.fgb
    --expect divmig_cluster_divide_migration_metrics.csv
    --expect-vector-feature-min divmig_cluster_divides.fgb:1
    --stdout-contains "cluster field    : cluster_id"
    --stdout-contains "cluster groups   : 2"
    --stdout-contains "selected basins  : 3"
    --stdout-contains "divide segments  : 1"
    --stdout-contains "heads with divide: 69"
)
lsdtt_set_tests_properties(smoke_divide_migration_clusters PROPERTIES
    LABELS "cli;divide_migration"
    FIXTURES_REQUIRED divide_migration_fixture)

# The same clustered run, but the cluster labels arrive in a VECTOR outlet file
# instead of a CSV. Regression test: read_vector_layer("") resolves through
# GetLayerByName(), which never matches an empty string, so this path used to
# fail outright for every vector outlet file. Asserting the same topology as the
# CSV run is the point -- the container must not change the answer.
add_lsdtt_cli_smoke_test(
    smoke_divide_migration_clusters_vector
    lsdtt-divide-migration
    divide_migration_clusters_vector.param.in
    divmig_cluster_vec_divides.fgb
    --read-path "${LSDTT_DIVIDE_MIGRATION_FIXTURE_DIR}"
    --expect divmig_cluster_vec_channel_heads.fgb
    --expect divmig_cluster_vec_divide_migration_metrics.csv
    --expect-vector-feature-min divmig_cluster_vec_divides.fgb:1
    --stdout-contains "cluster field    : cluster_id"
    --stdout-contains "cluster groups   : 2"
    --stdout-contains "selected basins  : 3"
    --stdout-contains "divide segments  : 1"
    --stdout-contains "heads with divide: 69"
)
lsdtt_set_tests_properties(smoke_divide_migration_clusters_vector PROPERTIES
    LABELS "cli;divide_migration"
    FIXTURES_REQUIRED divide_migration_fixture)

add_lsdtt_cli_smoke_test(
    failure_divide_migration_vector_field_case
    lsdtt-divide-migration
    divide_migration_clusters_vector_wrong_case.param.in
    unused_expected_output
    --read-path "${LSDTT_DIVIDE_MIGRATION_FIXTURE_DIR}"
    --expect-failure
    --stderr-contains "cluster field not found in vector layer: Cluster_ID"
)
lsdtt_set_tests_properties(failure_divide_migration_vector_field_case PROPERTIES
    LABELS "cli;divide_migration"
    FIXTURES_REQUIRED divide_migration_fixture)

# A clustered run whose outlets all collapse into ONE unit has no unit pair to
# compare. Counting snapped flats rather than units let this through and simply
# produced an empty divide set; the guard now rejects it and names both counts.
add_lsdtt_cli_smoke_test(
    failure_divide_migration_single_unit
    lsdtt-divide-migration
    divide_migration_single_unit_cluster.param.in
    unused_expected_output
    --read-path "${LSDTT_DIVIDE_MIGRATION_FIXTURE_DIR}"
    --expect-failure
    --stderr-contains "FATAL: fewer than two valid snapped outlets in two distinct units remain for unit mapping (outlets=3, units=1)."
    --stderr-contains "Divide extraction needs at least two units to compare."
)
lsdtt_set_tests_properties(failure_divide_migration_single_unit PROPERTIES
    LABELS "cli;divide_migration"
    FIXTURES_REQUIRED divide_migration_fixture)

# --- divide-fragment merging end to end -------------------------
#
# A second fixture DEM: one E-W ridge with a single valley draining north and
# two draining south, so the northern basin faces BOTH southern basins. The
# three-valley fixture above cannot test merging -- every unit pair there is
# adjacent along exactly one stretch of divide.
#
# The pair of tests is the assertion. Unclustered, all three basins are pairwise
# adjacent and the run emits THREE segments. Clustering the two southern outlets
# into one unit leaves that unit bordering the northern unit along TWO fragments
# that meet on the main ridge, and merging must turn those into ONE segment. A
# run that did not merge would report 2 here, and a fixture where the second
# fragment did not exist would fail the unclustered control.
add_lsdtt_cli_smoke_test(
    smoke_divide_migration_merge_unmerged
    lsdtt-divide-migration
    divide_migration_merge_unmerged.param.in
    divmig_merge_nc_divides.fgb
    --read-path "${LSDTT_DIVIDE_MIGRATION_FIXTURE_DIR}"
    --expect-vector-feature-min divmig_merge_nc_divides.fgb:3
    --stdout-contains "mapped units     : 3"
    --stdout-contains "selected basins  : 3"
    --stdout-contains "divide segments  : 3"
)
lsdtt_set_tests_properties(smoke_divide_migration_merge_unmerged PROPERTIES
    LABELS "cli;divide_migration"
    FIXTURES_REQUIRED divide_migration_fixture)

add_lsdtt_cli_smoke_test(
    smoke_divide_migration_merge_fragments
    lsdtt-divide-migration
    divide_migration_merge_fragments.param.in
    divmig_merge_divides.fgb
    --read-path "${LSDTT_DIVIDE_MIGRATION_FIXTURE_DIR}"
    --expect divmig_merge_divide_migration_metrics.csv
    --expect-vector-feature-min divmig_merge_divides.fgb:1
    --stdout-contains "cluster field    : cluster_id"
    --stdout-contains "cluster groups   : 2"
    --stdout-contains "selected basins  : 3"
    --stdout-contains "divide segments  : 1"
    --stdout-contains "heads without    : 0"
)
lsdtt_set_tests_properties(smoke_divide_migration_merge_fragments PROPERTIES
    LABELS "cli;divide_migration"
    FIXTURES_REQUIRED divide_migration_fixture)

# --- pair-summary statistics end to end ------------------------------------
#
# The same ridge fixture and clustering as the merge test above, so the topology
# is known: exactly two units bordering along one merged divide. The pair
# summary must therefore be exactly ONE row, on pair 1|2, with heads on both
# sides -- which is what makes the header assertion worth anything.
#
# The header is pinned in full because it is the published schema (§6.3 plus the
# per-metric valid-value counts this work added). The row assertion pins the
# leading pair_id/unit_a/unit_b, and the trailing uncertainty_method proves the
# parameter was read rather than defaulted: the CLI is run with bootstrap, not
# the std_err default.
add_lsdtt_cli_smoke_test(
    smoke_divide_migration_pair_summary
    lsdtt-divide-migration
    divide_migration_pair_summary.param.in
    divmig_pairsum_divide_pair_summary.csv
    --read-path "${LSDTT_DIVIDE_MIGRATION_FIXTURE_DIR}"
    --expect divmig_pairsum_divides.fgb
    --expect divmig_pairsum_divide_migration_metrics.csv
    --expect-file-contains "divmig_pairsum_divide_pair_summary.csv:pair_id,unit_a,unit_b,n_heads_a,n_heads_b,mean_head_elevation_a,mean_head_elevation_b,delta_head_elevation,n_head_elevation_a,n_head_elevation_b,mean_gradient_a,mean_gradient_b,delta_gradient,n_gradient_a,n_gradient_b,mean_relief_a,mean_relief_b,delta_relief,n_relief_a,n_relief_b,mean_chi_a,mean_chi_b,delta_chi,n_chi_a,n_chi_b,uncertainty_method,uncertainty_head_elevation,uncertainty_gradient,uncertainty_relief,uncertainty_chi"
    --expect-file-contains "divmig_pairsum_divide_pair_summary.csv:1|2,1,2,"
    --expect-file-contains "divmig_pairsum_divide_pair_summary.csv:,bootstrap,"
    --stdout-contains "divide segments  : 1"
    --stdout-contains "pair summaries   : 1"
    --stdout-contains "uncertainty      : bootstrap"
    --stdout-contains "Pair summaries report side A minus side B"
    --stdout-contains "heads without    : 0"
)
lsdtt_set_tests_properties(smoke_divide_migration_pair_summary PROPERTIES
    LABELS "cli;divide_migration"
    FIXTURES_REQUIRED divide_migration_fixture)

# --- lsdtt-lem: the in-LEM coupled hillslope path -----------------------------
#
# Why this exists (docs/developer_notes/testing_and_benchmarking.md §8.4): the semi-implicit fast
# path applied the DAMPED Picard omega to its single IMEX solve, under-diffusing
# every step. The standalone driver (unit_hillslope_implicit) could not see it --
# its steady-state tests all finish on full Picard -- while coupled LEM steady
# states came out about 40 % HIGH. A fast path that is only taken in regimes the
# unit tests do not finish in is a fast path the unit tests cannot validate.
#
# This run spends essentially all of its steps in that fast path: U*L/D = 0.2 is
# a quarter of Sc = 0.8, well under the max-face-ratio < 0.5 semi-implicit
# threshold. It therefore exercises exactly the path the unit driver structurally
# cannot reach.
#
# The assertion is the CONVERGED STEADY STATE, two-sided:
#   peak 4.9429 m -- the recorded in-LEM strip value (§8.3), against the analytic
#                    Roering peak of ~4.85 and the explicit model's 4.91
#   mean 3.2719 m -- pins the profile SHAPE, not just its crest
# Verified to fail on perturbation: reinstating the omega bug (step_omega =
# cur_omega for a lone semi-implicit solve) moves the peak to 6.862 m (+38.8 %)
# and the mean to 4.516 m, both far outside the bands. Starting from a FLAT
# initial profile is what stops "the solver did nothing" from landing inside the
# band -- inaction leaves the strip at 0 m.
#
# ~0.5 s, so it stays in the default suite.
add_lsdtt_cli_smoke_test(
    smoke_lem_strip_fixture
    lsdtt-synthetic-surfaces
    lem_strip_fixture.param.in
    lemstrip_hillslope_strip.tif
    --expect lemstrip_hillslope_strip_channel_mask.tif
    --expect-raster-size lemstrip_hillslope_strip.tif:101x11
)
lsdtt_set_tests_properties(smoke_lem_strip_fixture PROPERTIES
    LABELS "cli;hillslope"
    FIXTURES_SETUP lem_strip_fixture)

add_lsdtt_cli_smoke_test(
    smoke_synthetic_surfaces_phase5_defaults
    lsdtt-synthetic-surfaces
    synthetic_surfaces_phase5_defaults.param.in
    rasterdefaults5_uniform_random.tif
    --expect-raster-size rasterdefaults5_uniform_random.tif:128x128
    # Before migration: 256x256, seed 1, checksum 65320. The new oracle is
    # deliberately different in both shape and stochastic realisation.
    --expect-raster-checksum rasterdefaults5_uniform_random.tif:65532
)

add_lsdtt_cli_smoke_test(
    smoke_synthetic_surfaces_phase5_old_defaults_negative_control
    lsdtt-synthetic-surfaces
    synthetic_surfaces_phase5_old_defaults.param.in
    rasterolddefaults5_uniform_random.tif
    --expect-raster-size rasterolddefaults5_uniform_random.tif:256x256
    --expect-raster-checksum rasterolddefaults5_uniform_random.tif:65320
)

add_lsdtt_cli_smoke_test(
    smoke_synthetic_surfaces_phase5_legacy
    lsdtt-synthetic-surfaces
    synthetic_surfaces_phase5_legacy.param.in
    rasterlegacy5_hillslope_strip.tif
    --expect rasterlegacy5_hillslope_strip_channel_mask.tif
    --expect-raster-size rasterlegacy5_hillslope_strip.tif:9x7
    --expect-file-contains "rasterlegacy5_ingested_params.param:cell_size_metres: 2"
    --expect-file-contains "rasterlegacy5_ingested_params.param:random_seed: 19"
    --expect-file-contains "rasterlegacy5_ingested_params.param:strip_half_width_metres: 4"
    --expect-file-contains "rasterlegacy5_ingested_params.param:fan_lobe_half_angle_degrees: 21"
    --expect-file-contains "rasterlegacy5_ingested_params.param:fan_lobe_inset_metres: 3"
    --expect-file-contains "rasterlegacy5_ingested_params.param:fan_noise_amplitude_metres: 0.2"
    --expect-file-contains "rasterlegacy5_ingested_params.param:terrace_noise_amplitude_metres: 0.3"
    --expect-file-contains "rasterlegacy5_ingested_params.param:terrace_riser_gap_width_metres: 22"
    --expect-file-contains "rasterlegacy5_ingested_params.param:terrace_riser_height_metres: 7"
    --expect-file-contains "rasterlegacy5_ingested_params.param:terrace_tread_width_metres: 44"
    --expect-file-contains "rasterlegacy5_ingested_params.param:polyfit_lambda_metres: 77"
    --expect-file-contains "rasterlegacy5_ingested_params.param:write_alluvial_fan_lobe_id_raster: false"
    --expect-file-contains "rasterlegacy5_ingested_params.param:write_terrace_staircase_tread_id_raster: false"
    --stderr-contains "'cell_size' is deprecated"
    --stderr-contains "'seed' is deprecated"
    --stderr-contains "'strip_half_width' is deprecated"
)

add_lsdtt_cli_smoke_test(
    smoke_synthetic_surfaces_phase5_alias_conflict
    lsdtt-synthetic-surfaces
    synthetic_surfaces_phase5_alias_conflict.param.in
    rasterconflict5_uniform_random.tif
    --expect-raster-size rasterconflict5_uniform_random.tif:9x7
    --expect-file-contains "rasterconflict5_ingested_params.param:random_seed: 33"
    --stderr-contains "'random_seed' was also supplied"
)

add_lsdtt_cli_smoke_test(
    smoke_lem_phase5_legacy
    lsdtt-lem
    lem_phase5_legacy.param.in
    lemlegacy5_initial_topography.tif
    --expect-raster-size lemlegacy5_initial_topography.tif:6x9
    --expect-file-contains "lemlegacy5_ingested_params.param:n_rows: 9"
    --expect-file-contains "lemlegacy5_ingested_params.param:n_cols: 6"
    --expect-file-contains "lemlegacy5_ingested_params.param:cell_size_metres: 5"
    --expect-file-contains "lemlegacy5_ingested_params.param:noise_type: red_noise"
    --expect-file-contains "lemlegacy5_ingested_params.param:noise_target_stddev_metres: 0.25"
    --expect-file-contains "lemlegacy5_ingested_params.param:hillslope_implicit_picard_tolerance_metres: 0.0002"
    --expect-file-contains "lemlegacy5_ingested_params.param:hillslope_implicit_error_tolerance_metres: 0.0003"
    --expect-file-contains "lemlegacy5_ingested_params.param:steady_state_tolerance_dz_metres: 0.004"
    --expect-file-contains "lemlegacy5_ingested_params.param:incision_threshold_specific_area: 12.5"
    --expect-file-contains "lemlegacy5_ingested_params.param:channel_mask_threshold_specific_area: 8.5"
    --stderr-contains "'nrows' is deprecated"
    --stderr-contains "'ncols' is deprecated"
    --stderr-contains "'grid_resolution' is deprecated"
    --stderr-contains "'initial_surface_type' is deprecated"
    --stderr-contains "'noise_target_stddev' is deprecated"
)

add_lsdtt_cli_smoke_test(
    smoke_lem_phase5_alias_conflict
    lsdtt-lem
    lem_phase5_alias_conflict.param.in
    lemconflict5_initial_topography.tif
    --expect-raster-size lemconflict5_initial_topography.tif:6x4
    --stderr-contains "'n_rows' was also supplied"
)

add_lsdtt_cli_smoke_test(
    smoke_lem_phase5_defaults
    lsdtt-lem
    lem_phase5_defaults.param.in
    lemdefaults5_initial_topography.tif
    --expect-raster-size lemdefaults5_initial_topography.tif:400x200
    --expect-raster-checksum lemdefaults5_initial_topography.tif:10290
)

# Time-stepping snaps the initial surface to the fluvial steady state by default,
# so a run does not spend R/U building its relief from a flat start. Asserts the
# snap runs with no parameter asking for it, and actually produces relief.
add_lsdtt_cli_smoke_test(
    smoke_lem_snap_default
    lsdtt-lem
    lem_snap_default.param.in
    lemsnapdefault_initial_topography.tif
    --stdout-contains "--- Snap initial surface to fluvial steady state ---"
    --stdout-contains "Snapped relief"
)
lsdtt_set_tests_properties(smoke_lem_snap_default PROPERTIES LABELS "cli;lem")

# =============================================================================
# lsdtt-lem: spatially varying erodibility and the model grid (CLI)
#
# Fixtures (DEM, K rasters, placement sources and the expected placement
# results) are generated by the erodibility test driver into one directory,
# which every test below reads through --read-path. Value checks compare
# written rasters with --compare-rasters (dims, cell size, origin, CRS, mask,
# data bitwise) or --compare-raster-values (dims, mask, data), never against
# stored checksums. All runs: oooo, n_threads 1, snap off unless tested.
# =============================================================================
set(LEMK_FIXTURE_DIR "${CMAKE_CURRENT_BINARY_DIR}/lem_erodibility_cli_fixture")
add_test(NAME fixture_lem_erodibility_cli
         COMMAND lsdtt-test-lem-erodibility --write-fixtures "${LEMK_FIXTURE_DIR}")
lsdtt_set_tests_properties(fixture_lem_erodibility_cli PROPERTIES
    LABELS "fixture;lem"
    FIXTURES_SETUP lem_erodibility_cli)

function(lsdtt_add_lem_k_cli_test TEST_NAME TEMPLATE EXPECTED)
    add_lsdtt_cli_smoke_test(${TEST_NAME} lsdtt-lem ${TEMPLATE} ${EXPECTED}
        --read-path "${LEMK_FIXTURE_DIR}" ${ARGN})
    lsdtt_set_tests_properties(${TEST_NAME} PROPERTIES
        LABELS "cli;lem"
        FIXTURES_REQUIRED lem_erodibility_cli
        FIXTURES_SETUP ${TEST_NAME})
endfunction()

# A native comparison of two rasters written by the tests above. Registered
# only with the Python tests, whose outputs it reads.
#
# Evidence: when B is an expected raster the fixture writer computed from the
# fixture's layout, exact_solution; when B is another CLI run's output (uniform
# K field against the scalar), metamorphic. A system test either way: both
# rasters were written by lsdtt-lem or lsdtt-lithology-erodibility.
function(lsdtt_add_lem_k_compare TEST_NAME MODE A B)
    if(NOT LSDTT_ENABLE_PYTHON_TESTS)
        lsdtt_mark_test_unregistered(python ${TEST_NAME})
        return()
    endif()
    add_test(NAME ${TEST_NAME} COMMAND lsdtt-test-lem-erodibility ${MODE} "${A}" "${B}")
    lsdtt_set_tests_properties(${TEST_NAME} PROPERTIES
        LABELS "lem"
        FIXTURES_REQUIRED "${ARGN}")
    # LITHK_FIXTURE_DIR is set further down, after the first calls; an unset
    # directory must not turn into a "/" prefix that every path matches.
    set(_expected_from_layout FALSE)
    foreach(_dir IN ITEMS "${LEMK_FIXTURE_DIR}" "${LITHK_FIXTURE_DIR}")
        if(NOT _dir STREQUAL "")
            string(FIND "${B}" "${_dir}/" _at)
            if(_at EQUAL 0)
                set(_expected_from_layout TRUE)
            endif()
        endif()
    endforeach()
    if(_expected_from_layout)
        lsdtt_record_evidence(${TEST_NAME} system exact_solution)
    else()
        lsdtt_record_evidence(${TEST_NAME} system metamorphic)
    endif()
endfunction()

# C1-C4: a uniform K raster equal to erodibility is bitwise the scalar run.
lsdtt_add_lem_k_cli_test(smoke_lem_erodibility_scalar_d8
    lem_erodibility_scalar_d8.param.in lemkscalar_modelled_topography.tif)
lsdtt_add_lem_k_cli_test(smoke_lem_erodibility_uniform_d8
    lem_erodibility_uniform_d8.param.in lemkuniform_modelled_topography.tif
    --stdout-contains "Erodibility raster"
    --stdout-contains "m^(0.55) yr^-1")
lsdtt_add_lem_k_compare(regression_lem_erodibility_uniform_vs_scalar_d8 --compare-rasters
    "${CMAKE_CURRENT_BINARY_DIR}/smoke_lem_erodibility_scalar_d8/out/lemkscalar_modelled_topography.tif"
    "${CMAKE_CURRENT_BINARY_DIR}/smoke_lem_erodibility_uniform_d8/out/lemkuniform_modelled_topography.tif"
    "smoke_lem_erodibility_scalar_d8;smoke_lem_erodibility_uniform_d8")
lsdtt_add_lem_k_cli_test(smoke_lem_erodibility_scalar_dinf
    lem_erodibility_scalar_dinf.param.in lemkscalardinf_modelled_topography.tif)
lsdtt_add_lem_k_cli_test(smoke_lem_erodibility_uniform_dinf
    lem_erodibility_uniform_dinf.param.in lemkuniformdinf_modelled_topography.tif
    --stdout-contains "Erodibility raster")
lsdtt_add_lem_k_compare(regression_lem_erodibility_uniform_vs_scalar_dinf --compare-rasters
    "${CMAKE_CURRENT_BINARY_DIR}/smoke_lem_erodibility_scalar_dinf/out/lemkscalardinf_modelled_topography.tif"
    "${CMAKE_CURRENT_BINARY_DIR}/smoke_lem_erodibility_uniform_dinf/out/lemkuniformdinf_modelled_topography.tif"
    "smoke_lem_erodibility_scalar_dinf;smoke_lem_erodibility_uniform_dinf")

# C5: the K written (final and per frame) is the K that drove the step.
lsdtt_add_lem_k_cli_test(smoke_lem_erodibility_two_region
    lem_erodibility_two_region.param.in lemktwo_erodibility.tif
    --expect lemktwo_model_iteration_10_erodibility.tif
    --expect-raster-size lemktwo_erodibility.tif:20x40
    --expect-file-contains lemktwo_model_elevation_stats.csv:erodibility_field_mean)
lsdtt_add_lem_k_compare(regression_lem_erodibility_written_k --compare-rasters
    "${CMAKE_CURRENT_BINARY_DIR}/smoke_lem_erodibility_two_region/out/lemktwo_erodibility.tif"
    "${LEMK_FIXTURE_DIR}/lemk_k_two_region.tif"
    "smoke_lem_erodibility_two_region")
lsdtt_add_lem_k_compare(regression_lem_erodibility_written_k_frame --compare-rasters
    "${CMAKE_CURRENT_BINARY_DIR}/smoke_lem_erodibility_two_region/out/lemktwo_model_iteration_10_erodibility.tif"
    "${LEMK_FIXTURE_DIR}/lemk_k_two_region.tif"
    "smoke_lem_erodibility_two_region")

# C6: relief-mode steady state scales the K raster.
lsdtt_add_lem_k_cli_test(smoke_lem_erodibility_steady_state_relief
    lem_erodibility_steady_relief.param.in lemkss_steady_state_erodibility.tif
    --stdout-contains "erodibility_field_scale")

# C7: a forcing K step scales the raster exactly when it takes effect.
lsdtt_add_lem_k_cli_test(smoke_lem_erodibility_forcing_scale
    lem_erodibility_forcing_scale.param.in lemkforce_erodibility.tif
    --stdout-contains "x0.5 on the erodibility raster")
lsdtt_add_lem_k_compare(regression_lem_erodibility_forcing_before --compare-raster-values
    "${CMAKE_CURRENT_BINARY_DIR}/smoke_lem_erodibility_forcing_scale/out/lemkforce_model_iteration_10_erodibility.tif"
    "${LEMK_FIXTURE_DIR}/lemk_k_two_region.tif"
    "smoke_lem_erodibility_forcing_scale")
lsdtt_add_lem_k_compare(regression_lem_erodibility_forcing_after --compare-raster-values
    "${CMAKE_CURRENT_BINARY_DIR}/smoke_lem_erodibility_forcing_scale/out/lemkforce_model_iteration_15_erodibility.tif"
    "${LEMK_FIXTURE_DIR}/lemk_k_two_region_half.tif"
    "smoke_lem_erodibility_forcing_scale")

# C9-C13: refused K rasters.
lsdtt_add_lem_k_cli_test(failure_lem_erodibility_grid_mismatch
    lem_erodibility_grid_mismatch.param.in unused.tif --expect-exit-code 1
    --stderr-contains "erodibility raster grid" --stderr-contains "input_raster_placement")
lsdtt_add_lem_k_cli_test(failure_lem_erodibility_wrong_dims
    lem_erodibility_wrong_dims.param.in unused.tif --expect-exit-code 1
    --stderr-contains "40x20")
lsdtt_add_lem_k_cli_test(failure_lem_erodibility_nodata
    lem_erodibility_nodata.param.in unused.tif --expect-exit-code 1
    --stderr-contains "nodata" --stderr-contains "DEM-valid")
lsdtt_add_lem_k_cli_test(failure_lem_erodibility_negative
    lem_erodibility_negative.param.in unused.tif --expect-exit-code 1
    --stderr-contains "K must be > 0")
lsdtt_add_lem_k_cli_test(failure_lem_erodibility_zero
    lem_erodibility_zero.param.in unused.tif --expect-exit-code 1
    --stderr-contains "K must be > 0")

# C14: calibration and parameter checks.
lsdtt_add_lem_k_cli_test(failure_lem_erodibility_m_mismatch
    lem_erodibility_m_mismatch.param.in unused.tif --expect-exit-code 1
    --stderr-contains "m_exponent")
lsdtt_add_lem_k_cli_test(smoke_lem_erodibility_specific_area_rescaled
    lem_erodibility_specific_area_rescaled.param.in lemkspec_modelled_topography.tif
    --stdout-contains "calibration cell size 5 to 10")
lsdtt_add_lem_k_cli_test(smoke_lem_erodibility_map_area_conversion
    lem_erodibility_map_area_conversion.param.in lemkmap_modelled_topography.tif
    --stdout-contains "factor dx^m = 2.81838")
lsdtt_add_lem_k_cli_test(failure_lem_erodibility_forcing_without_reference
    lem_erodibility_forcing_without_reference.param.in unused.tif --expect-exit-code 1
    --stderr-contains "supply erodibility")
lsdtt_add_lem_k_cli_test(failure_lem_erodibility_bad_area_convention
    lem_erodibility_bad_area_convention.param.in unused.tif --expect-exit-code 1
    --stderr-contains "map_area")
lsdtt_add_lem_k_cli_test(failure_lem_bad_input_raster_placement
    lem_bad_input_raster_placement.param.in unused.tif --expect-exit-code 1
    --stderr-contains "input_raster_placement")
lsdtt_add_lem_k_cli_test(failure_lem_bad_model_grid_source
    lem_bad_model_grid_source.param.in unused.tif --expect-exit-code 1
    --stderr-contains "model_grid_source")

# C15: an error inside a time step is a clean FATAL, exit 1 (P0b).
lsdtt_add_lem_k_cli_test(failure_lem_erodibility_step_error_is_fatal
    lem_erodibility_step_error_is_fatal.param.in unused.tif --expect-exit-code 1
    --stderr-contains "FATAL: time step failed" --stderr-contains "erodibility"
    --stderr-not-contains "terminate called")

# C16: the snap with a uniform raster is bitwise the scalar snap; D7 uses the
# forcing history's K at t = 0 for the snap.
lsdtt_add_lem_k_cli_test(smoke_lem_erodibility_snap_uniform
    lem_erodibility_snap_uniform.param.in lemksnapu_initial_topography.tif)
lsdtt_add_lem_k_cli_test(smoke_lem_erodibility_snap_scalar
    lem_erodibility_snap_scalar.param.in lemksnaps_initial_topography.tif)
lsdtt_add_lem_k_compare(regression_lem_erodibility_snap_uniform_vs_scalar --compare-rasters
    "${CMAKE_CURRENT_BINARY_DIR}/smoke_lem_erodibility_snap_scalar/out/lemksnaps_initial_topography.tif"
    "${CMAKE_CURRENT_BINARY_DIR}/smoke_lem_erodibility_snap_uniform/out/lemksnapu_initial_topography.tif"
    "smoke_lem_erodibility_snap_scalar;smoke_lem_erodibility_snap_uniform")
lsdtt_add_lem_k_cli_test(smoke_lem_erodibility_snap_forcing_t0
    lem_erodibility_snap_forcing_t0.param.in lemksnapt0_initial_topography.tif
    --stdout-contains "K=0.0005")

# C17-C19: mode 2, the model grid from the K raster.
lsdtt_add_lem_k_cli_test(smoke_lem_grid_from_erodibility_raster
    lem_grid_from_erodibility_raster.param.in lemkgrid2_initial_topography.tif
    --expect-raster-size lemkgrid2_initial_topography.tif:80x100)
lsdtt_add_lem_k_compare(regression_lem_grid_from_raster --compare-rasters
    "${CMAKE_CURRENT_BINARY_DIR}/smoke_lem_grid_from_erodibility_raster/out/lemkgrid2_erodibility.tif"
    "${LEMK_FIXTURE_DIR}/lemk_k_src_ids.tif"
    "smoke_lem_grid_from_erodibility_raster")
lsdtt_add_lem_k_cli_test(failure_lem_grid_from_geographic_raster
    lem_grid_from_geographic_raster.param.in unused.tif --expect-exit-code 1
    --stderr-contains "geographic")
lsdtt_add_lem_k_cli_test(failure_lem_grid_source_conflict
    lem_grid_source_conflict.param.in unused.tif --expect-exit-code 1
    --stderr-contains "n_rows" --stderr-contains "model_grid_source")

# C20-C24: mode 1, placing a raster from anywhere on the synthetic grid.
lsdtt_add_lem_k_cli_test(smoke_lem_erodibility_relabel_whole_raster
    lem_erodibility_relabel_whole_raster.param.in lemkwhole_erodibility.tif
    --stdout-contains "relabelled" --stdout-contains "WGS 84 / UTM zone 32N")
lsdtt_add_lem_k_compare(regression_lem_erodibility_relabel_whole_raster --compare-raster-values
    "${CMAKE_CURRENT_BINARY_DIR}/smoke_lem_erodibility_relabel_whole_raster/out/lemkwhole_erodibility.tif"
    "${LEMK_FIXTURE_DIR}/lemk_k_src_ids.tif"
    "smoke_lem_erodibility_relabel_whole_raster")
lsdtt_add_lem_k_cli_test(smoke_lem_erodibility_centre_clip
    lem_erodibility_centre_clip.param.in lemkcentre_erodibility.tif)
lsdtt_add_lem_k_compare(regression_lem_erodibility_centre_clip --compare-raster-values
    "${CMAKE_CURRENT_BINARY_DIR}/smoke_lem_erodibility_centre_clip/out/lemkcentre_erodibility.tif"
    "${LEMK_FIXTURE_DIR}/lemk_expect_centre_clip.tif"
    "smoke_lem_erodibility_centre_clip")
lsdtt_add_lem_k_cli_test(smoke_lem_erodibility_lower_left_clip
    lem_erodibility_lower_left_clip.param.in lemkll_erodibility.tif)
lsdtt_add_lem_k_compare(regression_lem_erodibility_lower_left_clip --compare-raster-values
    "${CMAKE_CURRENT_BINARY_DIR}/smoke_lem_erodibility_lower_left_clip/out/lemkll_erodibility.tif"
    "${LEMK_FIXTURE_DIR}/lemk_expect_lower_left_clip.tif"
    "smoke_lem_erodibility_lower_left_clip")
lsdtt_add_lem_k_cli_test(failure_lem_erodibility_clip_off_source
    lem_erodibility_clip_off_source.param.in unused.tif --expect-exit-code 1
    --stderr-contains "runs off the source")
lsdtt_add_lem_k_cli_test(smoke_lem_erodibility_clip_resampled
    lem_erodibility_clip_resampled.param.in lemkres_erodibility.tif
    --stdout-contains "cell size 4 to 10")
lsdtt_add_lem_k_compare(regression_lem_erodibility_clip_resampled --compare-raster-values
    "${CMAKE_CURRENT_BINARY_DIR}/smoke_lem_erodibility_clip_resampled/out/lemkres_erodibility.tif"
    "${LEMK_FIXTURE_DIR}/lemk_expect_resampled_clip.tif"
    "smoke_lem_erodibility_clip_resampled")

# C25-C26: mode 2 with a DEM (warp), and relabel refused on a DEM.
lsdtt_add_lem_k_cli_test(smoke_lem_erodibility_warp_to_dem
    lem_erodibility_warp_to_dem.param.in lemkwarp_erodibility.tif
    --stdout-contains "warped")
lsdtt_add_lem_k_compare(regression_lem_erodibility_warp_to_dem --compare-raster-values
    "${CMAKE_CURRENT_BINARY_DIR}/smoke_lem_erodibility_warp_to_dem/out/lemkwarp_erodibility.tif"
    "${LEMK_FIXTURE_DIR}/lemk_expect_warped.tif"
    "smoke_lem_erodibility_warp_to_dem")
lsdtt_add_lem_k_cli_test(failure_lem_erodibility_relabel_on_dem
    lem_erodibility_relabel_on_dem.param.in unused.tif --expect-exit-code 1
    --stderr-contains "relabel" --stderr-contains "DEM")

# C27-C28 (D1c): the channel mask follows the same grid rules.
lsdtt_add_lem_k_cli_test(failure_lem_channel_mask_grid_mismatch
    lem_channel_mask_grid_mismatch.param.in unused.tif --expect-exit-code 1
    --stderr-contains "channel mask raster grid")
lsdtt_add_lem_k_cli_test(smoke_lem_channel_mask_relabel_synthetic
    lem_channel_mask_relabel_synthetic.param.in lemkmaskrel_modelled_topography.tif
    --stdout-contains "Loading channel mask:" --stdout-contains "relabelled")

# C29 (Q3): mode 2 footprint: the model domain is the K raster's footprint.
lsdtt_add_lem_k_cli_test(smoke_lem_grid_from_erodibility_raster_footprint
    lem_grid_from_erodibility_raster_footprint.param.in lemkfoot_initial_topography.tif
    --expect-valid-percent lemkfoot_initial_topography.tif:75
    --expect-valid-percent lemkfoot_modelled_topography.tif:75
    --stdout-contains "footprint" --stdout-contains "2000")
lsdtt_add_lem_k_compare(regression_lem_grid_footprint --compare-rasters
    "${CMAKE_CURRENT_BINARY_DIR}/smoke_lem_grid_from_erodibility_raster_footprint/out/lemkfoot_erodibility.tif"
    "${LEMK_FIXTURE_DIR}/lemk_k_src_footprint.tif"
    "smoke_lem_grid_from_erodibility_raster_footprint")

# C30-C31 (Q2): the scalar erodibility, and any explicit forcing K, must be > 0.
lsdtt_add_lem_k_cli_test(failure_lem_erodibility_scalar_zero
    lem_erodibility_scalar_zero.param.in unused.tif --expect-exit-code 1
    --stderr-contains "erodibility must be > 0")
lsdtt_add_lem_k_cli_test(failure_lem_erodibility_scalar_zero_no_fluvial
    lem_erodibility_scalar_zero_no_fluvial.param.in unused.tif --expect-exit-code 1
    --stderr-contains "erodibility must be > 0")
lsdtt_add_lem_k_cli_test(failure_lem_forcing_erodibility_zero
    lem_forcing_erodibility_zero.param.in unused.tif --expect-exit-code 1
    --stderr-contains "erodibility must be > 0" --stderr-contains "row")

# =============================================================================
# lsdtt-lithology-erodibility: lithology maps to K rasters for lsdtt-lem
#
# Fixtures (polygon maps, a code raster, a DEM, property tables, a crosswalk
# and the expected K rasters) come from lsdtt-test-lithology-erodibility
# --write-fixtures. The map's unit edges lie on pixel boundaries, so every
# expected K raster is computed from the layout, not by the code under test.
# =============================================================================
set(LITHK_FIXTURE_DIR "${CMAKE_CURRENT_BINARY_DIR}/lem_lithology_erodibility_cli_fixture")
add_test(NAME fixture_lem_lithology_erodibility_cli
         COMMAND lsdtt-test-lithology-erodibility --write-fixtures "${LITHK_FIXTURE_DIR}")
lsdtt_set_tests_properties(fixture_lem_lithology_erodibility_cli PROPERTIES
    LABELS "fixture;lem"
    FIXTURES_SETUP lem_lithology_erodibility_cli)

function(lsdtt_add_lithk_cli_test TEST_NAME TEMPLATE EXPECTED)
    add_lsdtt_cli_smoke_test(${TEST_NAME} lsdtt-lithology-erodibility ${TEMPLATE} ${EXPECTED}
        --read-path "${LITHK_FIXTURE_DIR}" ${ARGN})
    lsdtt_set_tests_properties(${TEST_NAME} PROPERTIES
        LABELS "cli;lem"
        FIXTURES_REQUIRED lem_lithology_erodibility_cli
        FIXTURES_SETUP ${TEST_NAME})
endfunction()

lsdtt_add_lithk_cli_test(smoke_lem_lithology_erodibility_own_grid
    lem_lithology_erodibility_own_grid.param.in lithkown_erodibility.tif
    --expect lithkown_lithology_code.tif
    --expect lithkown_lithology_legend.csv
    --expect lithkown_lem_erodibility_parameters.txt
    --expect-raster-size lithkown_erodibility.tif:80x100
    --expect-file-contains "lithkown_lem_erodibility_parameters.txt:erodibility_raster_area_convention: map_area"
    --stdout-contains "its footprint is the model domain")
lsdtt_add_lem_k_compare(regression_lem_lithology_erodibility_own_grid --compare-rasters
    "${CMAKE_CURRENT_BINARY_DIR}/smoke_lem_lithology_erodibility_own_grid/out/lithkown_erodibility.tif"
    "${LITHK_FIXTURE_DIR}/lithk_expect_own_grid.tif"
    "smoke_lem_lithology_erodibility_own_grid")
lsdtt_add_lithk_cli_test(smoke_lem_lithology_erodibility_dem_grid
    lem_lithology_erodibility_dem_grid.param.in lithkdem_erodibility.tif)
lsdtt_add_lem_k_compare(regression_lem_lithology_erodibility_dem_grid --compare-rasters
    "${CMAKE_CURRENT_BINARY_DIR}/smoke_lem_lithology_erodibility_dem_grid/out/lithkdem_erodibility.tif"
    "${LITHK_FIXTURE_DIR}/lithk_expect_dem_grid.tif"
    "smoke_lem_lithology_erodibility_dem_grid")
lsdtt_add_lithk_cli_test(smoke_lem_lithology_erodibility_centre_clip
    lem_lithology_erodibility_centre_clip.param.in lithkclip_erodibility.tif
    --stdout-contains "relabelled onto the synthetic grid")
lsdtt_add_lem_k_compare(regression_lem_lithology_erodibility_centre_clip --compare-raster-values
    "${CMAKE_CURRENT_BINARY_DIR}/smoke_lem_lithology_erodibility_centre_clip/out/lithkclip_erodibility.tif"
    "${LITHK_FIXTURE_DIR}/lithk_expect_centre_clip.tif"
    "smoke_lem_lithology_erodibility_centre_clip")
lsdtt_add_lithk_cli_test(smoke_lem_lithology_erodibility_code_raster
    lem_lithology_erodibility_code_raster.param.in lithkcode_erodibility.tif)
lsdtt_add_lem_k_compare(regression_lem_lithology_erodibility_code_raster --compare-rasters
    "${CMAKE_CURRENT_BINARY_DIR}/smoke_lem_lithology_erodibility_code_raster/out/lithkcode_erodibility.tif"
    "${LITHK_FIXTURE_DIR}/lithk_expect_own_grid.tif"
    "smoke_lem_lithology_erodibility_code_raster")
lsdtt_add_lithk_cli_test(smoke_lem_lithology_erodibility_crosswalk
    lem_lithology_erodibility_crosswalk.param.in lithkcw_erodibility.tif)
lsdtt_add_lem_k_compare(regression_lem_lithology_erodibility_crosswalk --compare-rasters
    "${CMAKE_CURRENT_BINARY_DIR}/smoke_lem_lithology_erodibility_crosswalk/out/lithkcw_erodibility.tif"
    "${LITHK_FIXTURE_DIR}/lithk_expect_own_grid.tif"
    "smoke_lem_lithology_erodibility_crosswalk")
lsdtt_add_lithk_cli_test(smoke_lem_lithology_erodibility_unmapped_code
    lem_lithology_erodibility_unmapped_code.param.in lithkfill_erodibility.tif
    --stdout-contains "given code 30")

lsdtt_add_lithk_cli_test(failure_lem_lithology_erodibility_unknown_value
    lem_lithology_erodibility_unknown_value.param.in unused.tif --expect-exit-code 1
    --stderr-contains "could not be given a code" --stderr-contains "'basalt'")
lsdtt_add_lithk_cli_test(failure_lem_lithology_erodibility_hole
    lem_lithology_erodibility_hole.param.in unused.tif --expect-exit-code 1
    --stderr-contains "inside the model domain have no lithology")
lsdtt_add_lithk_cli_test(failure_lem_lithology_erodibility_table_zero
    lem_lithology_erodibility_table_zero.param.in unused.tif --expect-exit-code 1
    --stderr-contains "K must be > 0")
lsdtt_add_lithk_cli_test(failure_lem_lithology_erodibility_table_no_m
    lem_lithology_erodibility_table_no_m.param.in unused.tif --expect-exit-code 1
    --stderr-contains "m_exponent")
lsdtt_add_lithk_cli_test(failure_lem_lithology_erodibility_vector_no_cell_size
    lem_lithology_erodibility_vector_no_cell_size.param.in unused.tif --expect-exit-code 1
    --stderr-contains "cell_size_metres")
lsdtt_add_lithk_cli_test(failure_lem_lithology_erodibility_relabel_on_dem
    lem_lithology_erodibility_relabel_on_dem.param.in unused.tif --expect-exit-code 1
    --stderr-contains "model_grid_source: parameters")
lsdtt_add_lithk_cli_test(failure_lem_lithology_erodibility_geographic
    lem_lithology_erodibility_geographic.param.in unused.tif --expect-exit-code 1
    --stderr-contains "geographic")

# The chain: lsdtt-lem runs on the K raster the own-grid test made, with the
# parameter lines it printed, on the map's own grid.
add_lsdtt_cli_smoke_test(
    smoke_lem_lithology_erodibility_chain
    lsdtt-lem
    lem_lithology_erodibility_chain.param.in
    lithkchain_modelled_topography.tif
    --read-path "${CMAKE_CURRENT_BINARY_DIR}/smoke_lem_lithology_erodibility_own_grid/out"
    --stdout-contains "converted map-area K to specific-area K at cell size 10"
    --expect-raster-size lithkchain_modelled_topography.tif:80x100)
lsdtt_set_tests_properties(smoke_lem_lithology_erodibility_chain PROPERTIES
    LABELS "cli;lem"
    FIXTURES_REQUIRED smoke_lem_lithology_erodibility_own_grid)

# Synthetic georeference. With no synthetic_crs_wkt the grid is centred on
# Watopia (UTM zone 58S); a supplied WKT is kept. The CRS and the grid centre
# are checked on the written raster by the erodibility test driver.
add_lsdtt_cli_smoke_test(
    smoke_lem_synthetic_default_georeference
    lsdtt-lem
    lem_synthetic_default_georeference.param.in
    lemsyngeoref_modelled_topography.tif
    --stdout-contains "default: Watopia, UTM zone 58S"
)
lsdtt_set_tests_properties(smoke_lem_synthetic_default_georeference PROPERTIES
    LABELS "cli;lem"
    FIXTURES_SETUP lem_synthetic_default_georeference)

add_lsdtt_cli_smoke_test(
    smoke_lem_synthetic_supplied_crs
    lsdtt-lem
    lem_synthetic_supplied_crs.param.in
    lemsyncrs_modelled_topography.tif
)
lsdtt_set_tests_properties(smoke_lem_synthetic_supplied_crs PROPERTIES
    LABELS "cli;lem"
    FIXTURES_SETUP lem_synthetic_supplied_crs)

if(LSDTT_ENABLE_PYTHON_TESTS)
    add_test(NAME regression_lem_synthetic_default_crs
             COMMAND lsdtt-test-lem-erodibility --expect-raster-epsg
                     "${CMAKE_CURRENT_BINARY_DIR}/smoke_lem_synthetic_default_georeference/out/lemsyngeoref_modelled_topography.tif"
                     32758 --expect-centre-near 30)
    lsdtt_set_tests_properties(regression_lem_synthetic_default_crs PROPERTIES
        LABELS "lem"
        FIXTURES_REQUIRED lem_synthetic_default_georeference)

    add_test(NAME regression_lem_synthetic_supplied_crs
             COMMAND lsdtt-test-lem-erodibility --expect-raster-epsg
                     "${CMAKE_CURRENT_BINARY_DIR}/smoke_lem_synthetic_supplied_crs/out/lemsyncrs_modelled_topography.tif"
                     32630)
    lsdtt_set_tests_properties(regression_lem_synthetic_supplied_crs PROPERTIES
        LABELS "lem"
        FIXTURES_REQUIRED lem_synthetic_supplied_crs)
else()
    lsdtt_mark_test_unregistered(python regression_lem_synthetic_default_crs)
    lsdtt_mark_test_unregistered(python regression_lem_synthetic_supplied_crs)
endif()

add_lsdtt_cli_smoke_test(
    smoke_lem_roering_implicit_strip
    lsdtt-lem
    lem_roering_implicit_strip.param.in
    lemstrip_run_modelled_topography.tif
    --read-path "${CMAKE_CURRENT_BINARY_DIR}/smoke_lem_strip_fixture/out"
    --expect lemstrip_run_model_elevation_stats.csv
    --expect-raster-stat lemstrip_run_modelled_topography.tif:max:4.90:4.99
    --expect-raster-stat lemstrip_run_modelled_topography.tif:mean:3.24:3.31
    --expect-raster-stat lemstrip_run_modelled_topography.tif:min:0.0:0.0
    --expect-valid-percent-100 lemstrip_run_modelled_topography.tif
    --stdout-contains "Hillslope transp. : on  (model=roering_implicit"
    --stdout-contains "Loading channel mask:"
)
lsdtt_set_tests_properties(smoke_lem_roering_implicit_strip PROPERTIES
    LABELS "cli;hillslope"
    FIXTURES_REQUIRED lem_strip_fixture)

# run_hillslope_on_loaded_dem: solve hillslopes WITHOUT re-snapping the fluvial
# field, so a transient (a step-change knickpoint) survives.
#
# The band is the point of the test, not decoration. On tiny_dem.asc (a 5x5
# planar ramp, 72-100 m) the mode must return elevations in that range with the
# channel cell untouched. The defect this guards against is real and was caught
# here: enforce_baselevel() pins every boundary pixel to baselevel_elevation,
# which is right for a synthetic LEM run but OVERWRITES a loaded DEM. Before the
# guard this run came back at 0.02-0.04 m everywhere -- the whole surface
# relaxed to a zeroed rim. A min band of 71-73 m fails instantly in that state.
add_lsdtt_cli_smoke_test(
    smoke_lem_hillslope_on_loaded_dem
    lsdtt-lem
    lem_hillslope_on_loaded_dem.param.in
    tinyhsonly_hillslope_on_loaded_dem.tif
    --expect-raster-size tinyhsonly_hillslope_on_loaded_dem.tif:5x5
    --expect-raster-stat tinyhsonly_hillslope_on_loaded_dem.tif:min:71.0:73.0
    --expect-raster-stat tinyhsonly_hillslope_on_loaded_dem.tif:max:85.0:101.0
    --stdout-contains "Hillslope solve on the loaded DEM (fluvial field untouched)"
    --stdout-contains "converged: yes"
)
lsdtt_set_tests_properties(smoke_lem_hillslope_on_loaded_dem PROPERTIES
    LABELS "cli;hillslope")

# critical_slope takes a different path: imposed by upslope integration, not
# iterated, so it reports relief and step count rather than iterations.
add_lsdtt_cli_smoke_test(
    smoke_lem_hillslope_on_loaded_dem_critical_slope
    lsdtt-lem
    lem_hillslope_on_loaded_dem_cs.param.in
    tinyhscs_hillslope_on_loaded_dem.tif
    --expect-raster-size tinyhscs_hillslope_on_loaded_dem.tif:5x5
    --expect-raster-stat tinyhscs_hillslope_on_loaded_dem.tif:min:71.0:73.0
    --stdout-contains "critical_slope (imposed by upslope integration, not iterated)"
    --stdout-contains "max relief above channel"
)
lsdtt_set_tests_properties(smoke_lem_hillslope_on_loaded_dem_critical_slope PROPERTIES
    LABELS "cli;hillslope")

# E2's coupled landscape has no threshold-defined channels. The critical-slope
# clamp is path integration rather than an elliptic solve, so open boundaries are
# valid anchors; this must not weaken the channel requirement for other solvers.
add_lsdtt_cli_smoke_test(
    smoke_lem_hillslope_on_loaded_dem_cs_coupled
    lsdtt-lem
    lem_hillslope_on_loaded_dem_cs_coupled.param.in
    tinyhscscoupled_hillslope_on_loaded_dem.tif
    --expect-raster-size tinyhscscoupled_hillslope_on_loaded_dem.tif:5x5
    --stdout-contains "Open-boundary anchor pixels"
    --stdout-contains "critical_slope (imposed by upslope integration, not iterated)"
)
lsdtt_set_tests_properties(smoke_lem_hillslope_on_loaded_dem_cs_coupled PROPERTIES
    LABELS "cli;hillslope")

add_lsdtt_cli_smoke_test(
    smoke_lem_cellular_particles
    lsdtt-lem
    lem_cellular_particles.param.in
    cellular_modelled_topography.tif
    --expect cellular_cellular_particles.csv
    --expect-file-contains "cellular_cellular_particles.csv:id,cell,weight,source_id,age_yr,distance_cells,status,mixing_depth,seed,substep_index"
    --stdout-contains "Hillslope transp. : on  (model=cellular"
    --stdout-contains "Writing cellular particles:"
)
lsdtt_set_tests_properties(smoke_lem_cellular_particles PROPERTIES
    LABELS "cli;hillslope;cellular"
    FIXTURES_SETUP cellular_particle_fixture)

if(LSDTT_ENABLE_PYTHON_TESTS)
    lsdtt_add_harness_test(smoke_lem_cellular_particles_restart
                --exe "$<TARGET_FILE:lsdtt-lem>"
                --param-template "${LSDTT_TEST_PARAM_DIR}/lem_cellular_particles_restart.param.in"
                --read-path "${CMAKE_CURRENT_BINARY_DIR}/smoke_lem_cellular_particles/out"
                --work-dir "${CMAKE_CURRENT_BINARY_DIR}/smoke_lem_cellular_particles_restart"
                --expect cellularrestart_modelled_topography.tif
                --expect cellularrestart_cellular_particles.csv
                --stdout-contains "Loading cellular particles:"
    )
else()
    lsdtt_mark_test_unregistered(python smoke_lem_cellular_particles_restart)
endif()
lsdtt_set_tests_properties(smoke_lem_cellular_particles_restart PROPERTIES
    LABELS "cli;hillslope;cellular"
    FIXTURES_REQUIRED cellular_particle_fixture)

add_lsdtt_cli_smoke_test(
    failure_lem_cellular_particles_wrong_model
    lsdtt-lem
    lem_cellular_particles_wrong_model.param.in
    unused.tif
    --expect-failure
    --stderr-contains "enable_cellular_particles requires hillslope_model: cellular"
)
add_lsdtt_cli_smoke_test(
    failure_lem_unknown_hillslope_model
    lsdtt-lem
    lem_unknown_hillslope_model.param.in
    unused.tif
    --expect-failure
    --stderr-contains "unknown hillslope_model 'typo_cellular'"
)

# The mode selects WHICH surface the solve runs on; it does not switch the solve
# on. Asking for one without the other is a configuration error, not a silent
# no-op.
add_lsdtt_cli_smoke_test(
    failure_lem_hillslope_on_loaded_dem_no_solve
    lsdtt-lem
    lem_hillslope_on_loaded_dem_no_solve.param.in
    tinyhsbad_hillslope_on_loaded_dem.tif
    --expect-failure
    --stderr-contains "run_hillslope_on_loaded_dem requires run_hillslope_steady_state: true"
)
lsdtt_set_tests_properties(failure_lem_hillslope_on_loaded_dem_no_solve PROPERTIES
    LABELS "cli;hillslope")

add_lsdtt_cli_smoke_test(
    smoke_lem_equilibration_cycle
    lsdtt-lem
    lem_equilibration_cycle_smoke.param.in
    tiny_eqcycle_steady_state_elevation.tif
    --expect tiny_eqcycle_snap_convergence.csv
    --expect-file-contains "tiny_eqcycle_snap_convergence.csv:max_abs_dz_m"
    --expect-file-contains "tiny_eqcycle_snap_convergence.csv:divide_step_median_m"
    --expect-file-contains "tiny_eqcycle_snap_convergence.csv:divide_gradient_median"
    --expect-file-contains "tiny_eqcycle_snap_convergence.csv:divide_gradient_p90"
    --expect-file-contains "tiny_eqcycle_snap_convergence.csv:divide_gradient_target_used"
    --expect-file-contains "tiny_eqcycle_snap_convergence.csv:divide_gradient_target_margin"
    --expect-file-contains "tiny_eqcycle_snap_convergence.csv:divide_gradient_target_met"
    --stdout-contains "Path           : landscape_equilibration_cycle"
)

add_lsdtt_cli_smoke_test(
    smoke_lem_equilibration_divide_step
    lsdtt-lem
    lem_equilibration_divide_step_smoke.param.in
    tiny_eqdiv_steady_state_elevation.tif
    --expect tiny_eqdiv_divide_step.tif
    --expect tiny_eqdiv_snap_convergence.csv
    --expect-raster-size tiny_eqdiv_divide_step.tif:5x5
    --expect-raster-stat tiny_eqdiv_divide_step.tif:max:0.0001:1000000.0
    --expect-file-contains "tiny_eqdiv_snap_convergence.csv:divide_step_median_m"
)

add_lsdtt_cli_smoke_test(
    smoke_lem_equilibration_divide_step_nodiag
    lsdtt-lem
    lem_equilibration_divide_step_nodiag_smoke.param.in
    tiny_eqdiv_nodiag_steady_state_elevation.tif
    --expect tiny_eqdiv_nodiag_divide_step.tif
    --expect tiny_eqdiv_nodiag_snap_convergence.csv
    --expect-file-contains "tiny_eqdiv_nodiag_snap_convergence.csv:divide_step_median_m"
)

add_lsdtt_cli_smoke_test(
    smoke_lem_equilibration_divide_target
    lsdtt-lem
    lem_equilibration_divide_target_smoke.param.in
    tiny_eqtarget_steady_state_elevation.tif
    --expect tiny_eqtarget_snap_convergence.csv
    --stdout-contains "Divide grad med:"
    --stdout-contains ", met)"
    --expect-file-contains "tiny_eqtarget_snap_convergence.csv:divide_gradient_median"
    --expect-file-contains "tiny_eqtarget_snap_convergence.csv:divide_gradient_p90"
    --expect-file-contains "tiny_eqtarget_snap_convergence.csv:divide_gradient_target_used"
    --expect-file-contains "tiny_eqtarget_snap_convergence.csv:divide_gradient_target_margin"
    --expect-file-contains "tiny_eqtarget_snap_convergence.csv:divide_gradient_target_met"
)

add_lsdtt_cli_smoke_test(
    smoke_lem_equilibration_divide_target_nodiag
    lsdtt-lem
    lem_equilibration_divide_target_nodiag_smoke.param.in
    tiny_eqtarget_nodiag_steady_state_elevation.tif
    --expect tiny_eqtarget_nodiag_snap_convergence.csv
    --stdout-contains "Divide grad med:"
    --stdout-contains ", met)"
    --expect-file-contains "tiny_eqtarget_nodiag_snap_convergence.csv:divide_gradient_target_met"
)

add_lsdtt_cli_smoke_test(
    smoke_lem_equilibration_divide_target_notmet
    lsdtt-lem
    lem_equilibration_divide_target_notmet_smoke.param.in
    tiny_eqtarget_notmet_steady_state_elevation.tif
    --expect tiny_eqtarget_notmet_snap_convergence.csv
    --stdout-contains "Divide grad med:"
    --stdout-contains "not met"
    --expect-file-contains "tiny_eqtarget_notmet_snap_convergence.csv:divide_gradient_median"
    --expect-file-contains "tiny_eqtarget_notmet_snap_convergence.csv:divide_gradient_p90"
    --expect-file-contains "tiny_eqtarget_notmet_snap_convergence.csv:divide_gradient_target_used"
    --expect-file-contains "tiny_eqtarget_notmet_snap_convergence.csv:divide_gradient_target_margin"
    --expect-file-contains "tiny_eqtarget_notmet_snap_convergence.csv:divide_gradient_target_met"
)

add_lsdtt_cli_smoke_test(
    smoke_lem_equilibration_convergence
    lsdtt-lem
    lem_equilibration_convergence_smoke.param.in
    tiny_eqconv_steady_state_elevation.tif
    --expect tiny_eqconv_snap_convergence.csv
    --stdout-contains "Converged at cycle"
)

add_lsdtt_cli_smoke_test(
    smoke_lem_equilibration_fixed_cycles
    lsdtt-lem
    lem_equilibration_fixed_cycles_smoke.param.in
    tiny_eqfixed_steady_state_elevation.tif
    --expect tiny_eqfixed_snap_convergence.csv
    --expect-file-contains "tiny_eqfixed_snap_convergence.csv:cycle"
    --stdout-contains "Cycles completed: 4"
)

add_lsdtt_cli_smoke_test(
    smoke_lem_equilibration_hillslope
    lsdtt-lem
    lem_equilibration_hillslope_smoke.param.in
    tiny_eqhs_steady_state_elevation.tif
    --expect tiny_eqhs_hillslope_steady_state_elevation.tif
    --expect tiny_eqhs_snap_convergence.csv
    --expect-file-contains "tiny_eqhs_snap_convergence.csv:frac_receiver_changed_hillslope"
    --expect-file-contains "tiny_eqhs_snap_convergence.csv:n_receiver_changed_hillslope"
    --expect-file-contains "tiny_eqhs_snap_convergence.csv:n_hillslope_compared"
    --expect-file-contains "tiny_eqhs_snap_convergence.csv:n_channel_pixels"
    --expect-file-contains "tiny_eqhs_snap_convergence.csv:hillslope_iterations"
    --expect-file-contains "tiny_eqhs_snap_convergence.csv:hillslope_residual"
    --expect-file-contains "tiny_eqhs_snap_convergence.csv:hillslope_converged"
    --stdout-contains "Hillslope      : linear (every cycle"
    --stdout-contains "skipping the one-shot post-step"
)

add_lsdtt_cli_smoke_test(
    failure_lem_equilibration_bad_divide_steps
    lsdtt-lem
    lem_equilibration_bad_divide_steps.param.in
    unused_expected_output
    --expect-failure
    --stderr-contains "divide_merge_search_steps must be > 0"
)

# -----------------------------------------------------------------------------
# hillslope_channel_coupling: naming the two coupling regimes.
#
# lsdtt3 has always contained BOTH couplings, and which one ran was an
# unreported side effect of incision_threshold_area. These pin the naming, the
# override, and the five configurations that are contradictions rather than
# settings.
# -----------------------------------------------------------------------------

# THE BACK-COMPATIBILITY GUARD. Parameter files across this project already
# depend on the implicit behaviour, so `auto` must still infer exactly what it
# always did: a positive threshold means channels are boundaries.
add_lsdtt_cli_smoke_test(
    smoke_lem_coupling_auto
    lsdtt-lem
    lem_coupling_auto.param.in
    tiny_coup_auto_modelled_topography.tif
    --stdout-contains "Coupling regime   : channels_as_boundaries (inferred"
    --stdout-contains "drainage density is IMPOSED"
)

# The parameter OVERRIDES the threshold rather than merely describing it: the
# same file with a positive threshold runs coupled when asked to.
add_lsdtt_cli_smoke_test(
    smoke_lem_coupling_coupled_overrides_threshold
    lsdtt-lem
    lem_coupling_coupled_overrides.param.in
    tiny_coup_override_modelled_topography.tif
    --stdout-contains "Coupling regime   : coupled (requested)"
    --stdout-contains "drainage density EMERGES"
    --stdout-contains "selected D-infinity routing"
    --stdout-contains "Dinf fluvial erosion   : on"
)

# The most numerically robust hillslope solver is the one that cannot run the
# emergent-drainage-density experiment. Fatal, and the message names the four
# that can.
add_lsdtt_cli_smoke_test(
    failure_lem_coupling_ri_steady_no_dirichlet
    lsdtt-lem
    lem_coupling_ri_steady_no_dirichlet.param.in
    unused_expected_output
    --expect-failure
    --stderr-contains "needs Dirichlet data somewhere"
)

# The combination this replaces (coupled + roering_implicit) used to be a FATAL
# error on the reasoning that the implicit solver needs a Dirichlet CHANNEL set.
# It does not -- it needs boundary conditions, and in the coupled regime the
# fluvial term is the sink. The restriction ruled the most robust solver out of
# the emergent-drainage-density experiment it is most useful for. What remains
# genuinely singular is a STEADY solve with no Dirichlet data anywhere, which the
# failure test above pins; the transient case must simply RUN.
add_lsdtt_cli_smoke_test(
    smoke_lem_coupling_ri_transient
    lsdtt-lem
    lem_coupling_ri_transient.param.in
    tiny_coup_ri_trans_modelled_topography.tif
    --stdout-contains "Model run complete"
)

# D8 on a diverging hillslope can manufacture channel heads, so the drainage
# density may be an artefact of the router. Since 2026-10-16 (3d9f37f) an
# explicit enable_d8_fluvial_erosion: true with coupled is honoured with a
# warning rather than refused; this pins that the run completes AND the
# warning is printed.
add_lsdtt_cli_smoke_test(
    smoke_lem_coupling_d8_warning
    lsdtt-lem
    lem_coupling_d8_conflict.param.in
    tiny_coup_d8_modelled_topography.tif
    --stdout-contains "Model run complete"
    --stderr-contains "enable_d8_fluvial_erosion: true"
    --stderr-contains "D8 can manufacture channel heads"
    --stderr-contains "Continuing because D8/coupled was explicitly"
)

# The analytical steady-state path snaps channels then solves hillslopes against
# them, so it IS the boundary-condition model by construction.
add_lsdtt_cli_smoke_test(
    failure_lem_coupling_steady_conflict
    lsdtt-lem
    lem_coupling_steady_conflict.param.in
    unused_expected_output
    --expect-failure
    --stderr-contains "run_steady_state"
    --stderr-contains "boundary-condition model"
)

# Asking to impose a drainage density while supplying no channels.
add_lsdtt_cli_smoke_test(
    failure_lem_coupling_no_channels
    lsdtt-lem
    lem_coupling_no_channels.param.in
    unused_expected_output
    --expect-failure
    --stderr-contains "channels_as_boundaries needs a channel set"
)

# A typo must not quietly select one of two different experiments.
add_lsdtt_cli_smoke_test(
    failure_lem_coupling_typo
    lsdtt-lem
    lem_coupling_typo.param.in
    unused_expected_output
    --expect-failure
    --stderr-contains "unrecognised hillslope_channel_coupling"
)

add_lsdtt_cli_smoke_test(
    failure_lem_equilibration_dinf_conflict
    lsdtt-lem
    lem_equilibration_dinf_conflict.param.in
    unused_expected_output
    --expect-failure
    --stderr-contains "enable_dinf_fluvial_erosion"
    --stderr-contains "use_landscape_equilibration_cycle"
)

add_lsdtt_cli_smoke_test(
    failure_lem_equilibration_bad_patience
    lsdtt-lem
    lem_equilibration_bad_patience.param.in
    unused_expected_output
    --expect-failure
    --stderr-contains "steady_state_convergence_patience"
)

add_lsdtt_cli_smoke_test(
    failure_hydraulic_routing_invalid_mode
    lsdtt-graphflood
    hydraulic_routing_invalid_mode.param.in
    unused_expected_output
    --expect-failure
    --stderr-contains "FATAL: flow_routing_method must be 'mfd' or 'dinf'."
)

add_lsdtt_cli_smoke_test(
    rigorous_tif_fill_raster_size
    lsdtt-raster-preprocessing
    raster_preprocessing_tif_smoke.param.in
    tiny30_fill.tif
    --expect-raster-size tiny30_fill.tif:32x39
)

add_lsdtt_cli_smoke_test(
    rigorous_raster_calculator_difference_checksum
    lsdtt-raster-calculator
    raster_calculator_difference.param.in
    tiny_and_tiny_dem_plus1_difference.tif
    --expect-raster-checksum tiny_and_tiny_dem_plus1_difference.tif:65511
)

add_lsdtt_cli_smoke_test(
    rigorous_drainage_d8_checksum
    lsdtt-drainage-analysis
    drainage_analysis_smoke.param.in
    tiny_d8_drainage_area.tif
    --expect-raster-checksum tiny_d8_drainage_area.tif:112
)

add_lsdtt_cli_smoke_test(
    rigorous_chi_nodes_sha256
    lsdtt-chi-analysis
    chi_analysis_smoke.param.in
    tiny_chi_fit_nodes_mn_0p500.csv
    # Re-pinned 2026-09-08 when the chi integral moved from a left-endpoint
    # rectangle rule to a TRAPEZOIDAL one. Chi is lower everywhere on this
    # fixture (max 4.0297 -> 3.6297) because the old rule evaluated the
    # integrand at the upstream end of every segment, where it is largest.
    # Verified against TopoToolbox before re-pinning: the per-basin regression
    # slope went from 1.00377 to 1.00000. See chi_increment() and
    # docs/developer_notes/testing_and_benchmarking.md.
    --expect-sha256 tiny_chi_fit_nodes_mn_0p500.csv:b64aa1266abc09369b423534d44ec5ed8a3ceac52d7d31e9f8637c8893c26ad3
)

# -----------------------------------------------------------------------------
# Flow-routing benchmark vs RichDEM (docs/developer_notes/testing_and_benchmarking.md)
#
# OFF by default: these tests need network access for the fixtures and an
# optional RichDEM install, neither of which may turn CI red.  When enabled but
# unavailable, the driver exits 77 and CTest records a skip.
#
# RichDEM comes from conda-forge via pixi -- NOT pip, whose richdem 0.3.4 pins
# numpy<2 and will downgrade numpy in whatever environment it lands in:
#   pixi run fetch-benchmark-dem
#   pixi run -e benchmark flow-benchmark --exe build/lsdtt-raster-preprocessing
# -----------------------------------------------------------------------------
option(LSDTT_ENABLE_BENCHMARKS "Enable RichDEM flow-routing benchmark tests" OFF)

if(LSDTT_ENABLE_BENCHMARKS)
    # Launch through `pixi run`, not the environment's python directly.  pixi
    # exports PROJ_DATA / GDAL_DATA when it activates the environment; invoking
    # .pixi/envs/benchmark/bin/python straight from CTest skips activation and
    # GDAL then fails to open the PROJ database at write time.
    find_program(LSDTT_PIXI_EXECUTABLE pixi
                 HINTS "$ENV{HOME}/.pixi/bin" "$ENV{PIXI_HOME}/bin")

    if(NOT LSDTT_PIXI_EXECUTABLE)
        message(WARNING
            "LSDTT_ENABLE_BENCHMARKS is ON but pixi was not found. "
            "Install pixi (https://pixi.sh) or set LSDTT_PIXI_EXECUTABLE. "
            "Benchmark tests will not be registered.")
    else()
        set(LSDTT_BENCHMARK_LAUNCHER
            "${LSDTT_PIXI_EXECUTABLE}" run
            --manifest-path "${CMAKE_SOURCE_DIR}/pixi.toml"
            -e benchmark python)

    # accum_d8 is included because it is asserted against TopoToolbox, whose D8
    # uses the same steepest-descent rule as lsdtt3. The RichDEM row for the
    # same comparison is recorded but not asserted -- see REPORT_ONLY in the
    # driver and docs/developer_notes/testing_and_benchmarking.md.
    #
    # flats_float_step is a VALUE comparison: fill_flat_resolution_choice:
    # barnes2014b_float_step against RichDEM FillDepressions + ResolveFlats,
    # asserted bit for bit on float32 outside lsdtt3's |z| < 1 m near-zero
    # guard and the cells lsdtt3 reports changing after its first
    # (RichDEM-equivalent) round, which must all lie above RichDEM's, with no
    # cell left without a lower neighbour -- oracle O1
    # (benchmark_flat_resolution_richdem) on real terrain.
    #
    # breach_least_cost is a VALUE comparison too: lsdtt3's least-cost breach
    # (Lindsay 2020, ported from Whitebox) against whitebox_workflows 2.0.6
    # under matched parameters (edge policy whitebox, no remainder fill,
    # max_dist 100, minimize_dist false, flat increment float32(0.001)),
    # asserted bit for bit on float32.
    foreach(bench_fixture alpine_small lakeland_small)
        foreach(bench_case fill fill_epsilon breach breach_least_cost accum_d8 flats_float_step)
            add_test(
                NAME benchmark_${bench_case}_${bench_fixture}
                COMMAND ${LSDTT_BENCHMARK_LAUNCHER}
                        "${CMAKE_CURRENT_SOURCE_DIR}/scripts/run_flow_benchmark.py"
                        --exe "$<TARGET_FILE:lsdtt-raster-preprocessing>"
                        --fixture ${bench_fixture}
                        --comparison ${bench_case}
                        --work-dir "${CMAKE_CURRENT_BINARY_DIR}/flow_benchmark"
                        --check
            )
            lsdtt_mark_test_runtime(pixi benchmark_${bench_case}_${bench_fixture})
            lsdtt_set_tests_properties(benchmark_${bench_case}_${bench_fixture}
                PROPERTIES LABELS "benchmark;flow" SKIP_RETURN_CODE 77)
        endforeach()
    endforeach()

    # Least-cost breach fuzz against Whitebox 2.0.6: 200 seeded random DEMs
    # (integer and 0.125-step values, some with nodata blobs) x minimize_dist
    # {false,true} x max_dist {3,10,100} x max_cost {unlimited,2.5}, every run
    # bitwise on float32. No network: the DEMs are generated in the driver.
    # Exits 77 without whitebox_workflows.
    add_test(
        NAME benchmark_breach_least_cost_fuzz
        COMMAND ${LSDTT_BENCHMARK_LAUNCHER}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/check_breach_least_cost_vs_whitebox.py"
                --exe "$<TARGET_FILE:lsdtt-raster-preprocessing>"
                --work-dir "${CMAKE_CURRENT_BINARY_DIR}/flow_benchmark"
    )
    lsdtt_mark_test_runtime(pixi benchmark_breach_least_cost_fuzz)
    lsdtt_set_tests_properties(benchmark_breach_least_cost_fuzz
        PROPERTIES LABELS "benchmark;flow" SKIP_RETURN_CODE 77 TIMEOUT 1800)

    # Surface-metrics SPEED benchmark vs TopoToolbox.  TopoToolbox is a PERFORMANCE reference here, not a
    # correctness one -- correctness comes from the Tier A analytic oracle in
    # lsdtt-test-surface-metrics, which needs no network.  So --check enforces
    # ONLY the rank-correlation guard: it fails when the two tools turn out not
    # to be computing the same kind of quantity, which would make the timing
    # meaningless.  No wall-clock threshold is ever asserted.
    #
    # The _medium fixtures, not the _small ones the flow benchmark uses: at
    # 256x256 the polyfit is ~5 ms and every timing is noise.
    foreach(surf_fixture alpine_medium lakeland_medium)
        foreach(surf_case slope curvature_all4 hillshade)
            add_test(
                NAME benchmark_surface_${surf_case}_${surf_fixture}
                COMMAND ${LSDTT_BENCHMARK_LAUNCHER}
                        "${CMAKE_CURRENT_SOURCE_DIR}/scripts/run_surface_benchmark.py"
                        --exe "$<TARGET_FILE:lsdtt-surface-metrics>"
                        --fixture ${surf_fixture}
                        --comparison ${surf_case}
                        --repeats 5
                        --skip-studies
                        --work-dir "${CMAKE_CURRENT_BINARY_DIR}/surface_benchmark"
                        --check
            )
            lsdtt_mark_test_runtime(pixi benchmark_surface_${surf_case}_${surf_fixture})
            lsdtt_set_tests_properties(benchmark_surface_${surf_case}_${surf_fixture}
                PROPERTIES LABELS "benchmark;surface_metrics" SKIP_RETURN_CODE 77
                           TIMEOUT 1800)
        endforeach()
    endforeach()

    # -------------------------------------------------------------------------
    # Chi VALUE comparison against TopoToolbox.
    #
    # This is the odd one out in this block, and deliberately so: everything
    # else here is a TIMING benchmark that happens to guard comparability,
    # whereas this asserts VALUES and nothing else. Both tools integrate
    # chi = int (A0/A)^(m/n) ds; with A0 = 1 m^2, m/n matched and (since
    # 2026-09-08) the same trapezoidal quadrature on both sides, they agree to
    # float32 rounding. The pinned slope window is [0.999999, 1.000001].
    #
    # It is registered here rather than with the unit tests because it needs
    # the pixi benchmark environment for TopoToolbox, and exits 77 without it.
    #
    # NO NETWORK AND NO COMMITTED RASTER. The fixture is a synthetic island the
    # driver generates itself; lsdtt-synthetic-networks is deterministic under
    # a fixed seed, and the driver verifies the DEM's sha256 before applying any
    # tolerance. A crop of real terrain cannot be used at all: lsdtt3 drops any
    # basin whose upslope area touches invalid data, which on a cropped DEM is
    # every basin, so there would be no chi to compare.
    #
    # Timings from this driver are never asserted -- they would be meaningless
    # on a shared or thermally limited machine. --check enforces value metrics
    # only.
    # -------------------------------------------------------------------------
    add_test(
        NAME benchmark_chi_value_topotoolbox
        COMMAND ${LSDTT_BENCHMARK_LAUNCHER}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/run_svector_benchmark.py"
                --build-dir "$<TARGET_FILE_DIR:lsdtt-chi-analysis>"
                --cases chi_value
                --threshold 100
                --out-dir "${CMAKE_CURRENT_BINARY_DIR}/svector_benchmark"
                --check
    )
    lsdtt_mark_test_runtime(pixi benchmark_chi_value_topotoolbox)
    lsdtt_set_tests_properties(benchmark_chi_value_topotoolbox
        PROPERTIES LABELS "benchmark;chi" SKIP_RETURN_CODE 77 TIMEOUT 1800)

    # -------------------------------------------------------------------------
    # Flat resolution VALUE comparison against RichDEM (oracle O1).
    #
    # The first round of fill_flat_resolution_choice: barnes2014b_float_step
    # is meant to equal RichDEM's ResolveFlats on float32 bit for bit. The
    # driver fills its own synthetic fixtures (no network) and the committed
    # lakeland crop with vincent1993 at fill_min_slope 0, resolves the flats
    # with both tools and asserts equal float32 bits on every valid cell,
    # except (a) on flats with |z| < 1 m, where lsdtt3's near-zero guard
    # replaces subnormal steps by design (the driver asserts the differences
    # are exactly those cells, in RichDEM's order), and (b) on the cells
    # lsdtt3 reports changing after round 1, which must all be above RichDEM.
    # It asserts lsdtt3 leaves no cell without a lower neighbour, that RichDEM
    # does on the fixtures chosen for it, and the round-1 PAIR counts.
    # Elevations are compared, not D8 directions: RichDEM's D8 is not steepest
    # descent (see run_flow_benchmark.py). Exits 77 without richdem.
    # -------------------------------------------------------------------------
    add_test(
        NAME benchmark_flat_resolution_richdem
        COMMAND ${LSDTT_BENCHMARK_LAUNCHER} -I
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/check_flat_resolution_vs_richdem.py"
                --bin-dir "$<TARGET_FILE_DIR:lsdtt-raster-preprocessing>"
                --work-dir "${CMAKE_CURRENT_BINARY_DIR}/flat_resolution_richdem"
    )
    lsdtt_mark_test_runtime(pixi benchmark_flat_resolution_richdem)
    lsdtt_set_tests_properties(benchmark_flat_resolution_richdem
        PROPERTIES LABELS "benchmark;flow" SKIP_RETURN_CODE 77 TIMEOUT 600)
    endif()

    # Hydraulic DAGGER parity benchmark.
    # In --check mode it enforces mask agreement and pinned parity tolerances.
    # Prefer launching through pixi's benchmark env when available so GDAL/PROJ
    # and Python deps match the documented benchmark path.
    set(LSDTT_DAGGER_ROOT_DEFAULT "${CMAKE_SOURCE_DIR}/../DAGGER")
    set(LSDTT_DAGGER_REF_DEFAULT
        "${CMAKE_CURRENT_SOURCE_DIR}/scripts/dagger_graphflood_reference.py")

    if(LSDTT_PIXI_EXECUTABLE)
        set(LSDTT_HYDRAULIC_BENCHMARK_LAUNCHER ${LSDTT_BENCHMARK_LAUNCHER})
    elseif(LSDTT_ENABLE_PYTHON_TESTS)
        set(LSDTT_HYDRAULIC_BENCHMARK_LAUNCHER ${Python3_EXECUTABLE})
    else()
        set(LSDTT_HYDRAULIC_BENCHMARK_LAUNCHER "")
        message(STATUS "Hydraulic DAGGER benchmarks not registered: they need pixi "
                       "or LSDTT_ENABLE_PYTHON_TESTS=ON.")
    endif()

    if(LSDTT_HYDRAULIC_BENCHMARK_LAUNCHER)
        foreach(hyd_fixture tiny_dem.asc alpine_small)
            foreach(hyd_scenario stationary_fixed transient_dynamic)
                add_test(
                    NAME benchmark_hydraulic_dagger_${hyd_scenario}_${hyd_fixture}
                    COMMAND ${LSDTT_HYDRAULIC_BENCHMARK_LAUNCHER}
                            "${CMAKE_CURRENT_SOURCE_DIR}/scripts/run_hydraulic_dagger_benchmark.py"
                            --exe "$<TARGET_FILE:lsdtt-graphflood>"
                            --fetch-exe "$<TARGET_FILE:lsdtt-fetch-raster>"
                            --dagger-root "${LSDTT_DAGGER_ROOT_DEFAULT}"
                            --dagger-ref "${LSDTT_DAGGER_REF_DEFAULT}"
                            --fixture ${hyd_fixture}
                            --scenario ${hyd_scenario}
                            --work-dir "${CMAKE_CURRENT_BINARY_DIR}/hydraulic_benchmark/${hyd_scenario}_${hyd_fixture}"
                            --write-json
                            --check
                )
                if(LSDTT_PIXI_EXECUTABLE)
                    lsdtt_mark_test_runtime(pixi benchmark_hydraulic_dagger_${hyd_scenario}_${hyd_fixture})
                else()
                    lsdtt_mark_test_runtime(python benchmark_hydraulic_dagger_${hyd_scenario}_${hyd_fixture})
                endif()
                lsdtt_set_tests_properties(
                    benchmark_hydraulic_dagger_${hyd_scenario}_${hyd_fixture}
                    PROPERTIES LABELS "benchmark;hydraulic" SKIP_RETURN_CODE 77
                               TIMEOUT 1800)
            endforeach()
        endforeach()
    endif()
endif()

# -----------------------------------------------------------------------------
# C++ unit tests
# -----------------------------------------------------------------------------

# Surface-metrics validation: Tier A analytic oracle (global quadratic and
# periodic-wave surfaces with exact Taylor coefficients, sphere, cylinder,
# plane) plus Tier C property tests (algebraic identities, SMAP classification,
# rotation invariance, thread determinism, MADk2 invariants).  Synthetic
# fixtures only: no network and no committed DEM, so it runs in the default
# suite.
#
# Image-filter primitives and splitmix64: border mapping, box mean with
# OpenCV's rounding rules, Gaussian kernels, separable and fixed-point filters,
# Otsu, Sobel, Bresenham.
# TOPMODEL column physics, index mode.
add_executable(lsdtt-test-topmodel-column cpp/test_topmodel_column.cpp)
target_link_libraries(lsdtt-test-topmodel-column PRIVATE lsdtt3_cpu)
add_test(NAME unit_topmodel_column COMMAND lsdtt-test-topmodel-column)
lsdtt_set_tests_properties(unit_topmodel_column PROPERTIES LABELS "network")

# HYPE port. The oracle under
# tests/expected/hype is generated by tests/scripts/generate_hype_oracle.py;
# these tests only read committed files.
add_executable(lsdtt-test-distributed-hydrology-forcing cpp/test_distributed_hydrology_forcing.cpp)
target_link_libraries(lsdtt-test-distributed-hydrology-forcing PRIVATE lsdtt3_cpu)
add_test(NAME unit_distributed_hydrology_forcing
         COMMAND lsdtt-test-distributed-hydrology-forcing
                 "${CMAKE_CURRENT_BINARY_DIR}/distributed_hydrology_forcing_output")
lsdtt_set_tests_properties(unit_distributed_hydrology_forcing PROPERTIES LABELS "distributed_hydrology")

add_executable(lsdtt-test-hype-parameters cpp/test_hype_parameters.cpp)
target_link_libraries(lsdtt-test-hype-parameters PRIVATE lsdtt3_cpu)
add_test(NAME unit_hype_parameters
         COMMAND lsdtt-test-hype-parameters "${CMAKE_CURRENT_SOURCE_DIR}/expected/hype"
                 "${CMAKE_CURRENT_BINARY_DIR}/hype_parameters_output")
lsdtt_set_tests_properties(unit_hype_parameters PROPERTIES LABELS "distributed_hydrology")

add_executable(lsdtt-test-hype-column cpp/test_hype_column.cpp)
target_link_libraries(lsdtt-test-hype-column PRIVATE lsdtt3_cpu)
add_test(NAME unit_hype_column
         COMMAND lsdtt-test-hype-column "${CMAKE_CURRENT_SOURCE_DIR}/expected/hype")
lsdtt_set_tests_properties(unit_hype_column PROPERTIES LABELS "distributed_hydrology")

add_executable(lsdtt-test-hype-routing cpp/test_hype_routing.cpp)
target_link_libraries(lsdtt-test-hype-routing PRIVATE lsdtt3_cpu)
add_test(NAME unit_hype_routing
         COMMAND lsdtt-test-hype-routing "${CMAKE_CURRENT_SOURCE_DIR}/expected/hype")
lsdtt_set_tests_properties(unit_hype_routing PROPERTIES LABELS "distributed_hydrology")

add_executable(lsdtt-test-hype-model cpp/test_hype_model.cpp)
target_link_libraries(lsdtt-test-hype-model PRIVATE lsdtt3_cpu)
add_test(NAME unit_hype_model
         COMMAND lsdtt-test-hype-model "${CMAKE_CURRENT_SOURCE_DIR}/expected/hype")
lsdtt_set_tests_properties(unit_hype_model PROPERTIES LABELS "distributed_hydrology")

add_executable(lsdtt-test-hype-raster-setup cpp/test_hype_raster_setup.cpp)
target_link_libraries(lsdtt-test-hype-raster-setup PRIVATE lsdtt3_cpu)
add_test(NAME unit_hype_raster_setup
         COMMAND lsdtt-test-hype-raster-setup "${CMAKE_CURRENT_BINARY_DIR}/hype_raster_setup_output")
lsdtt_set_tests_properties(unit_hype_raster_setup PROPERTIES LABELS "distributed_hydrology")

# The generated input set for the column_model: hype CLI tests (smoke and
# failure tests under lsdtt-distributed-hydrology). Generated in C++, never
# committed.
add_executable(lsdtt-test-hype-cli-fixture cpp/test_hype_cli_fixture.cpp)
target_link_libraries(lsdtt-test-hype-cli-fixture PRIVATE lsdtt3_cpu)
add_test(NAME fixture_hype_cli
         COMMAND lsdtt-test-hype-cli-fixture "${LSDTT_HYPE_CLI_FIXTURE_DIR}")
lsdtt_set_tests_properties(fixture_hype_cli PROPERTIES
    LABELS "fixture;distributed_hydrology"
    FIXTURES_SETUP hype_cli_fixture)

# HYPE thread-count invariance: 1 thread
# vs N threads must be bitwise identical. Guards the Stage B parallel steps.
add_executable(lsdtt-test-hype-thread-invariance cpp/test_hype_thread_invariance.cpp)
target_link_libraries(lsdtt-test-hype-thread-invariance PRIVATE lsdtt3_cpu)
add_test(NAME unit_hype_thread_invariance COMMAND lsdtt-test-hype-thread-invariance)
lsdtt_set_tests_properties(unit_hype_thread_invariance PROPERTIES LABELS "distributed_hydrology")

# HYPE timing harness. Built always so it
# can be run by hand; registered only under LSDTT_ENABLE_PERF_TESTS below.
add_executable(lsdtt-test-hype-performance cpp/test_hype_performance.cpp)
target_link_libraries(lsdtt-test-hype-performance PRIVATE lsdtt3_cpu)

add_executable(lsdtt-test-subbasin-partition cpp/test_subbasin_partition.cpp)
target_link_libraries(lsdtt-test-subbasin-partition PRIVATE lsdtt3_cpu)
add_test(NAME unit_subbasin_partition COMMAND lsdtt-test-subbasin-partition)
lsdtt_set_tests_properties(unit_subbasin_partition PROPERTIES LABELS "distributed_hydrology")

# The channel mask ladder. Takes an output
# directory for the small GeoJSON fixture it writes at run time, so nothing is
# committed to tests/data.
add_executable(lsdtt-test-channel-extent cpp/test_channel_extent.cpp)
target_link_libraries(lsdtt-test-channel-extent PRIVATE lsdtt3_cpu)
add_test(NAME unit_channel_extent
         COMMAND lsdtt-test-channel-extent "${CMAKE_CURRENT_BINARY_DIR}/channel_extent_output")
lsdtt_set_tests_properties(unit_channel_extent PROPERTIES LABELS "network")
file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/channel_extent_output")

# The routing ladder.
add_executable(lsdtt-test-hydro-routing cpp/test_hydro_routing.cpp)
target_link_libraries(lsdtt-test-hydro-routing PRIVATE lsdtt3_cpu)
add_test(NAME unit_hydro_routing COMMAND lsdtt-test-hydro-routing)
lsdtt_set_tests_properties(unit_hydro_routing PROPERTIES LABELS "network")

# Canopy interception.
add_executable(lsdtt-test-canopy-interception cpp/test_canopy_interception.cpp)
target_link_libraries(lsdtt-test-canopy-interception PRIVATE lsdtt3_cpu)
add_test(NAME unit_canopy_interception COMMAND lsdtt-test-canopy-interception)
lsdtt_set_tests_properties(unit_canopy_interception PROPERTIES LABELS "network")

# The roughness partition and the channel-n
# ladder.
add_executable(lsdtt-test-hydro-roughness cpp/test_hydro_roughness.cpp)
target_link_libraries(lsdtt-test-hydro-roughness PRIVATE lsdtt3_cpu)
add_test(NAME unit_hydro_roughness COMMAND lsdtt-test-hydro-roughness)
lsdtt_set_tests_properties(unit_hydro_roughness PROPERTIES LABELS "network")

# Land-cover scenarios.
add_executable(lsdtt-test-landcover-scenario cpp/test_landcover_scenario.cpp)
target_link_libraries(lsdtt-test-landcover-scenario PRIVATE lsdtt3_cpu)
add_test(NAME unit_landcover_scenario COMMAND lsdtt-test-landcover-scenario)
lsdtt_set_tests_properties(unit_landcover_scenario PROPERTIES LABELS "network")

# Infiltration excess (exp-K Green-Ampt).
# Series, closed-form capacity, water balance and dimensional-form checks; no
# external reference input.
add_executable(lsdtt-test-infiltration-excess cpp/test_infiltration_excess.cpp)
target_link_libraries(lsdtt-test-infiltration-excess PRIVATE lsdtt3_cpu)
add_test(NAME unit_infiltration_excess COMMAND lsdtt-test-infiltration-excess)
lsdtt_set_tests_properties(unit_infiltration_excess PROPERTIES LABELS "network")

# The dynamic (prognostic) saturated zone.
add_executable(lsdtt-test-topmodel-dynamic cpp/test_topmodel_dynamic.cpp)
target_link_libraries(lsdtt-test-topmodel-dynamic PRIVATE lsdtt3_cpu)
add_test(NAME unit_topmodel_dynamic COMMAND lsdtt-test-topmodel-dynamic)
lsdtt_set_tests_properties(unit_topmodel_dynamic PROPERTIES LABELS "network")

# Velocity field, travel time, width
# function. "hydro" is claimed by the wetness_index -> network domain rule
# below, so this stays in the network compartment with its Phase 0 sibling.
add_executable(lsdtt-test-hydro-velocity-field cpp/test_hydro_velocity_field.cpp)
target_link_libraries(lsdtt-test-hydro-velocity-field PRIVATE lsdtt3_cpu)
add_test(NAME unit_hydro_velocity_field COMMAND lsdtt-test-hydro-velocity-field)
lsdtt_set_tests_properties(unit_hydro_velocity_field PROPERTIES LABELS "network")

# The shared TWI. Guards the core function
# directly, including the float precision contract that the CLI checksum test
# was MEASURED not to discriminate -- see the header of test_wetness_index.cpp.
add_executable(lsdtt-test-wetness-index cpp/test_wetness_index.cpp)
target_link_libraries(lsdtt-test-wetness-index PRIVATE lsdtt3_cpu)
add_test(NAME unit_wetness_index COMMAND lsdtt-test-wetness-index)
lsdtt_set_tests_properties(unit_wetness_index PROPERTIES LABELS "network")

add_executable(lsdtt-test-image-filters cpp/test_image_filters.cpp)
target_link_libraries(lsdtt-test-image-filters PRIVATE lsdtt3_utils)
add_test(NAME unit_image_filters COMMAND lsdtt-test-image-filters)
lsdtt_set_tests_properties(unit_image_filters PROPERTIES LABELS "image_ops")

# Native raster image operations vs outputs frozen from OpenCV 4.13 before the
# OpenCV dependency was removed (tests/expected/image_ops/), plus invariants.
add_executable(lsdtt-test-image-ops cpp/test_image_ops.cpp)
target_link_libraries(lsdtt-test-image-ops PRIVATE lsdtt3_cpu)
add_test(NAME unit_image_ops
         COMMAND lsdtt-test-image-ops ${CMAKE_CURRENT_SOURCE_DIR}/expected/image_ops)
lsdtt_set_tests_properties(unit_image_ops PROPERTIES LABELS "image_ops")

# NOTE: two MADk2 checks are KNOWN FAILURES (the anisotropy/DMC estimator treats
# axial directions as vectors).  They are reported as KNOWN-FAIL and do not fail
# the run; flip k_madk2_axial_is_fixed in the driver once the estimator is
# corrected.
add_executable(lsdtt-test-surface-metrics cpp/test_surface_metrics.cpp)
target_link_libraries(lsdtt-test-surface-metrics PRIVATE lsdtt3_cpu)

file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/surface_metrics_output)
add_test(NAME unit_surface_metrics
         COMMAND lsdtt-test-surface-metrics
                 ${CMAKE_CURRENT_BINARY_DIR}/surface_metrics_output)
lsdtt_set_tests_properties(unit_surface_metrics PROPERTIES LABELS "surface_metrics")

# Vertical units: the feet/metres trap and the plausibility warning that makes
# it discoverable (docs/developer_notes/testing_and_benchmarking.md).
add_executable(lsdtt-test-vertical-units cpp/test_vertical_units.cpp)
target_link_libraries(lsdtt-test-vertical-units PRIVATE lsdtt3_cpu)

add_test(NAME unit_vertical_units COMMAND lsdtt-test-vertical-units)
lsdtt_set_tests_properties(unit_vertical_units PROPERTIES LABELS "raster")

# -----------------------------------------------------------------------------
# lsdtt::gdal_init() -- PROJ/GDAL runtime data pinning.
#
# set_bundled_data_root() only bites before the first gdal_init(), and
# gdal_init() runs once per process, so each scenario is its own invocation.
# See tests/cpp/test_gdal_init.cpp for what each mode asserts.
# -----------------------------------------------------------------------------
add_executable(lsdtt-test-gdal-init cpp/test_gdal_init.cpp)
target_link_libraries(lsdtt-test-gdal-init PRIVATE lsdtt3_core)
target_include_directories(lsdtt-test-gdal-init PRIVATE ${GDAL_INCLUDE_DIRS})

add_test(NAME unit_gdal_init_idempotent
         COMMAND lsdtt-test-gdal-init idempotent)
lsdtt_set_tests_properties(unit_gdal_init_idempotent PROPERTIES LABELS "gdal")

# The case the module exists for: a poisoned PROJ_DATA must be ignored when a
# bundled proj.db is present. Needs a real proj.db to stage; skips without one.
add_test(NAME unit_gdal_init_pinned
         COMMAND lsdtt-test-gdal-init pinned
                 "${CMAKE_CURRENT_BINARY_DIR}/gdal_init_staging")
lsdtt_set_tests_properties(unit_gdal_init_pinned PROPERTIES
    LABELS "gdal"
    SKIP_RETURN_CODE 77
    ENVIRONMENT "LSDTT_TEST_PROJ_DB=${LSDTT_TEST_PROJ_DB}")

# The conda-forge case: nothing bundled, so the user's environment must still
# be honoured. This is the half that a naive "always override" would break.
add_test(NAME unit_gdal_init_fallthrough
         COMMAND lsdtt-test-gdal-init fallthrough)
lsdtt_set_tests_properties(unit_gdal_init_fallthrough PROPERTIES LABELS "gdal")

# Phase 1 of the GDAL/PROJ friction work: GDAL, OGR and CPL headers are confined
# to .cpp files under src/. This is a source check, not a runtime one -- a
# violation compiles perfectly well, it just quietly undoes the facade, so
# nothing but an explicit check will catch it.
# Chunk 2a of the GDAL/PROJ friction work: the format contract. Asserts every
# driver in the required set is registered and that lsdtt3 can read a raster
# written by each. Written against the system GDAL; re-run unchanged against the
# vendored GDAL in Phase 2, where it is the acceptance gate for the trimmed
# driver set.
add_executable(lsdtt-test-format-support cpp/test_format_support.cpp)
target_link_libraries(lsdtt-test-format-support PRIVATE lsdtt3_core)
target_include_directories(lsdtt-test-format-support PRIVATE ${GDAL_INCLUDE_DIRS})
add_test(NAME unit_gdal_format_support
         COMMAND lsdtt-test-format-support
                 "${CMAKE_CURRENT_BINARY_DIR}/format_support_work")
lsdtt_set_tests_properties(unit_gdal_format_support PROPERTIES LABELS "gdal")

# The loader refuses rotated, non-square-pixel and south-up rasters, which its
# meta (origin + one cell size) would otherwise misread for every CLI.
add_executable(lsdtt-test-data-loader-geotransform-guard cpp/test_data_loader_geotransform_guard.cpp)
target_link_libraries(lsdtt-test-data-loader-geotransform-guard PRIVATE lsdtt3_core)
target_include_directories(lsdtt-test-data-loader-geotransform-guard PRIVATE ${GDAL_INCLUDE_DIRS})
add_test(NAME unit_raster_loader_geotransform_guard
         COMMAND lsdtt-test-data-loader-geotransform-guard
                 "${CMAKE_CURRENT_BINARY_DIR}/geotransform_guard_work")
lsdtt_set_tests_properties(unit_raster_loader_geotransform_guard PROPERTIES LABELS "gdal")

# giuh_discharge: volume conservation and timestep invariance. Added 2026-09-16
# with the fix for the unused dt_rainfall_s defect; the GIUH code had no tests
# at all before this, which is how the defect survived.
add_executable(lsdtt-test-giuh-discharge cpp/test_giuh_discharge.cpp)
target_link_libraries(lsdtt-test-giuh-discharge PRIVATE lsdtt3_cpu)
add_test(NAME unit_giuh_discharge COMMAND lsdtt-test-giuh-discharge)
lsdtt_set_tests_properties(unit_giuh_discharge PROPERTIES LABELS "giuh")

if(LSDTT_ENABLE_PYTHON_TESTS)
    add_test(NAME unit_gdal_init_boundary
        COMMAND ${Python3_EXECUTABLE}
            "${CMAKE_CURRENT_SOURCE_DIR}/scripts/check_gdal_boundary.py")
    lsdtt_mark_test_runtime(python unit_gdal_init_boundary)
else()
    lsdtt_mark_test_unregistered(python unit_gdal_init_boundary)
endif()
lsdtt_set_tests_properties(unit_gdal_init_boundary PROPERTIES
    LABELS "gdal;conventions")

# -----------------------------------------------------------------------------
# Phase 2 of the CLI consistency clean-up.
#
# Each of these pins a helper that was extracted out of one or more cli/*.cpp
# files into a library home. They are registered here, ahead of the extraction,
# so that no Phase 2 agent has to edit this file and collide with another.
# -----------------------------------------------------------------------------

# B1 to_lower_ascii + Category A planar_distance / polyline_length /
# cell_centre_x,y / bearing_from_north_degrees.
add_executable(lsdtt-test-cli-helper-string-utils cpp/test_cli_helper_string_utils.cpp)
target_link_libraries(lsdtt-test-cli-helper-string-utils PRIVATE lsdtt3_cpu)
add_test(NAME unit_cli_helper_string_utils COMMAND lsdtt-test-cli-helper-string-utils)

# B2 output_extension(raster_output_format, read_ext) and B3 csv_number.
add_executable(lsdtt-test-cli-helper-formats cpp/test_cli_helper_formats.cpp)
target_link_libraries(lsdtt-test-cli-helper-formats PRIVATE lsdtt3_core)
add_test(NAME unit_cli_helper_formats COMMAND lsdtt-test-cli-helper-formats)

# D5 quote-aware lsdtt::csv_split_row and load_csv_table.
add_executable(lsdtt-test-csv-table-quoting cpp/test_csv_table_quoting.cpp)
target_link_libraries(lsdtt-test-csv-table-quoting PRIVATE lsdtt3_core)
add_test(NAME unit_csv_table_quoting COMMAND lsdtt-test-csv-table-quoting)

# Strict field converters shared by the hand-written CSV readers (Manning,
# crn_column erosion history, particle restart): file, line and text on error.
add_executable(lsdtt-test-csv-table-field-parsers cpp/test_csv_table_field_parsers.cpp)
target_link_libraries(lsdtt-test-csv-table-field-parsers PRIVATE lsdtt3_cpu)
add_test(NAME unit_csv_table_field_parsers
         COMMAND lsdtt-test-csv-table-field-parsers
                 ${CMAKE_CURRENT_BINARY_DIR}/csv_table_field_parsers_output)

# B4 the polyfit metric name/suffix table beside the polyfit_metric enum.
add_executable(lsdtt-test-surface-metrics-polyfit-metric-table
               cpp/test_surface_metrics_polyfit_metric_table.cpp)
target_link_libraries(lsdtt-test-surface-metrics-polyfit-metric-table PRIVATE lsdtt3_cpu)
add_test(NAME unit_surface_metrics_polyfit_metric_table
         COMMAND lsdtt-test-surface-metrics-polyfit-metric-table)
lsdtt_set_tests_properties(unit_surface_metrics_polyfit_metric_table
                     PROPERTIES LABELS "surface_metrics")

add_executable(lsdtt-test-basin-outline-metrics
               cpp/test_basin_outline_metrics.cpp)
target_link_libraries(lsdtt-test-basin-outline-metrics PRIVATE lsdtt3_cpu)
add_test(NAME unit_basin_outline_metrics
         COMMAND lsdtt-test-basin-outline-metrics)
lsdtt_set_tests_properties(unit_basin_outline_metrics
                     PROPERTIES LABELS "basin_surface_metrics")

add_executable(lsdtt-test-param-parser-help cpp/test_param_parser_help.cpp)
target_link_libraries(lsdtt-test-param-parser-help PRIVATE lsdtt3_core)

add_test(NAME unit_param_parser_help COMMAND lsdtt-test-param-parser-help)

add_executable(lsdtt-test-common-parameters cpp/test_common_parameters.cpp)
target_link_libraries(lsdtt-test-common-parameters PRIVATE lsdtt3_core)

add_test(NAME unit_common_parameters COMMAND lsdtt-test-common-parameters)
lsdtt_set_tests_properties(unit_common_parameters PROPERTIES LABELS "core;parameters")

# The shared threshold picker (fixed | qq_gaussian | cdf_divergence): the
# Clubb et al. (2017) QQ rule against samples with a known answer, and the
# persistence rule on written QQ tables. cdf_divergence is pinned end to end by
# rigorous_flat_surface_finder_cdf_divergence and
# rigorous_mountain_front_cdf_divergence.
add_executable(lsdtt-test-threshold-selection cpp/test_threshold_selection.cpp)
target_link_libraries(lsdtt-test-threshold-selection PRIVATE lsdtt3_utils)
add_test(NAME unit_fluvial_surface_threshold_selection
         COMMAND lsdtt-test-threshold-selection)
lsdtt_set_tests_properties(unit_fluvial_surface_threshold_selection PROPERTIES
    LABELS "fluvial_surfaces")

# Cosmogenic physics driver (docs/developer_notes/testing_and_benchmarking.md).  Writes its
# plottable CSVs into the test working directory.
add_executable(lsdtt-test-cosmo-physics cpp/test_cosmo_physics.cpp)
target_link_libraries(lsdtt-test-cosmo-physics PRIVATE lsdtt3_cpu)

add_executable(lsdtt-test-cosmo-ratio cpp/test_cosmo_ratio.cpp)
target_link_libraries(lsdtt-test-cosmo-ratio PRIVATE lsdtt3_cpu)

add_executable(lsdtt-test-cosmo-units cpp/test_cosmo_units.cpp)
target_link_libraries(lsdtt-test-cosmo-units PRIVATE lsdtt3_cpu)

add_executable(lsdtt-test-raster-transform cpp/test_raster_transform.cpp)
target_link_libraries(lsdtt-test-raster-transform PRIVATE lsdtt3_cpu)

# Analytical step-change rasters: the Royden & Perron / Mitchell & Yanites
# transient solution in chi, plus the imposed critical-slope hillslopes.
# Synthetic ridge-and-valley fixture, no network, no external data.
add_executable(lsdtt-test-step-change cpp/test_step_change.cpp)
target_link_libraries(lsdtt-test-step-change PRIVATE lsdtt3_cpu)

# Landslide size distribution: the imposed magnitude-frequency
# distribution and the area-volume scaling, with no landscape involved.
# Synthetic only -- no DEM, no network.
add_executable(lsdtt-test-landslide-size cpp/test_landslide_size.cpp)
target_link_libraries(lsdtt-test-landslide-size PRIVATE lsdtt3_cpu)

# Reproducing Niemi et al. (2005): the ensemble of basin samples, the
# apparent-versus-true erosion rate, and the area dependence of the scatter.
add_executable(lsdtt-test-niemi-ensemble cpp/test_niemi_ensemble.cpp)
target_link_libraries(lsdtt-test-niemi-ensemble PRIVATE lsdtt3_cpu)

# Quartz fertility as a mixing weight, and the ice mask that shares its
# machinery. Synthetic rasters only.
add_executable(lsdtt-test-quartz-fertility cpp/test_quartz_fertility.cpp)
target_link_libraries(lsdtt-test-quartz-fertility PRIVATE lsdtt3_cpu)

# Landslide scars superimposed on a step change, plus the density convention
# (real metres vs mass depth, and the per-material scar density). Synthetic
# rasters only -- no DEM, no flow routing.
add_executable(lsdtt-test-landslide-step-change cpp/test_landslide_step_change.cpp)
target_link_libraries(lsdtt-test-landslide-step-change PRIVATE lsdtt3_cpu)

# Landslide placement: eligible hosts and the s_vector footprint.
# Synthetic ridge-and-valley fixture, no network.
add_executable(lsdtt-test-landslide-placement cpp/test_landslide_placement.cpp)
target_link_libraries(lsdtt-test-landslide-placement PRIVATE lsdtt3_cpu)

# Landslide event: budget -> N -> placed scars, with mass balance.
add_executable(lsdtt-test-landslide-event cpp/test_landslide_event.cpp)
target_link_libraries(lsdtt-test-landslide-event PRIVATE lsdtt3_cpu)

# Landslide -> CRN coupling (T16). Composes existing crn_column
# calls only; no new cosmogenic physics.
add_executable(lsdtt-test-landslide-crn cpp/test_landslide_crn.cpp)
target_link_libraries(lsdtt-test-landslide-crn PRIVATE lsdtt3_cpu)

# Landslide series: multi-event integration, plus the
# Niemi rate-closure and Yanites depth-regime reproductions.
add_executable(lsdtt-test-landslide-series cpp/test_landslide_series.cpp)
target_link_libraries(lsdtt-test-landslide-series PRIVATE lsdtt3_cpu)

# Flow routing / depression handling property tests
# (docs/developer_notes/testing_and_benchmarking.md).  Synthetic fixtures only: no
# network, no RichDEM, so this runs in the default CTest suite.
add_executable(lsdtt-test-flow-properties cpp/test_flow_properties.cpp)
target_link_libraries(lsdtt-test-flow-properties PRIVATE lsdtt3_cpu)

add_executable(lsdtt-test-flow-snapping cpp/test_flow_snapping.cpp)
target_link_libraries(lsdtt-test-flow-snapping PRIVATE lsdtt3_cpu)

add_executable(lsdtt-test-mountain-front-outlet-pruning
               cpp/test_mountain_front_outlet_pruning.cpp)
target_link_libraries(lsdtt-test-mountain-front-outlet-pruning PRIVATE lsdtt3_cpu)

add_executable(lsdtt-test-synthetic-networks cpp/test_synthetic_networks.cpp)
target_link_libraries(lsdtt-test-synthetic-networks PRIVATE lsdtt3_cpu)

add_test(NAME unit_synthetic_networks
         COMMAND lsdtt-test-synthetic-networks
                 ${LSDTT_SYNTHETIC_NETWORKS_FIXTURE_DIR})
# Both synthetic-networks runs write into the same directory (the same binary
# writes ocn_square_polygon.fgb in both modes), so they must not run at once:
# under ctest -j the polygon write in one could fail mid-write of the other
# (P1 failed in ~1 run in 4, measured 2026-10-07). RESOURCE_LOCK serialises
# just these two and keeps every output path unchanged.
lsdtt_set_tests_properties(unit_synthetic_networks PROPERTIES
    LABELS "synthetic_networks"
    FIXTURES_SETUP synthetic_networks_fixture
    RESOURCE_LOCK synthetic_networks_output)
file(MAKE_DIRECTORY ${LSDTT_SYNTHETIC_NETWORKS_FIXTURE_DIR})

# Reported, never asserted: the fidelity / energy table over temperature at
# FIXED n_iter (anneal length dominates temperature, so a sweep that let n_iter
# float would confound the two). Writes ocn_fidelity_vs_temperature.csv, which
# is the table the downstream science reads to choose a temperature.
# MEASURED 2026-08-13: 1.3 s for 15 anneals of 64x64 x 163840 iterations, so it
# is NOT in LSDTT_HEAVY_TESTS - it is cheap enough for the inner loop and the
# table is then always fresh.
add_test(NAME unit_synthetic_networks_sweep
         COMMAND lsdtt-test-synthetic-networks
                 ${LSDTT_SYNTHETIC_NETWORKS_FIXTURE_DIR} --sweep)
lsdtt_set_tests_properties(unit_synthetic_networks_sweep PROPERTIES
    LABELS "synthetic_networks" TIMEOUT 1800
    RESOURCE_LOCK synthetic_networks_output)

add_test(NAME unit_flow_properties
         COMMAND lsdtt-test-flow-properties
                 ${CMAKE_CURRENT_BINARY_DIR}/flow_properties_output)
lsdtt_set_tests_properties(unit_flow_properties PROPERTIES LABELS "flow")
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/flow_properties_output)

# BYTE-IDENTITY PIN for the D8 flow object. Pins were generated from the code
# BEFORE the section-9.4 speed refactor of lsdtt_flowinfo (commit 7e81a39); a
# speed change must never need them regenerated. Integer-terrain fixtures only,
# so the hashes are reproducible on any machine. See the file header.
add_executable(lsdtt-test-flow-d8-identity cpp/test_flow_d8_identity.cpp)
target_link_libraries(lsdtt-test-flow-d8-identity PRIVATE lsdtt3_cpu)
add_test(NAME unit_flow_d8_identity COMMAND lsdtt-test-flow-d8-identity)
lsdtt_set_tests_properties(unit_flow_d8_identity PROPERTIES LABELS "flow")

# Barnes et al. (2014b) flat resolution. Hermetic: in-process inputs only.
add_executable(lsdtt-test-flat-resolution cpp/test_flat_resolution.cpp)
target_link_libraries(lsdtt-test-flat-resolution PRIVATE lsdtt3_cpu)
add_test(NAME unit_flat_resolution COMMAND lsdtt-test-flat-resolution)
lsdtt_set_tests_properties(unit_flat_resolution PROPERTIES LABELS "flow")
# The same driver on the committed real-terrain crop (Step 11d): a Copernicus
# GLO-30 lakeland window on which a single Barnes pass leaves cells the D8
# router refuses. Reads tests/data only; no network, no skip path.
add_test(NAME unit_flat_resolution_real_terrain
         COMMAND lsdtt-test-flat-resolution --real-terrain
                 "${LSDTT_TEST_DATA_DIR}/flat_resolution_lakeland_crop_dem.tif")
lsdtt_set_tests_properties(unit_flat_resolution_real_terrain PROPERTIES LABELS "flow")

# Least-cost breaching (Lindsay 2020) against frozen float32 outputs of
# Whitebox next-gen BreachDepressionsLeastCost (whitebox_workflows 2.0.6),
# bitwise, on hand-built fixtures F0-F15, with edge_policy whitebox; plus the
# Whitebox counters, max_dist validation and edge_policy outlet property checks.
# Hermetic: in-process inputs only.
add_executable(lsdtt-test-breach-least-cost cpp/test_breach_least_cost.cpp)
target_link_libraries(lsdtt-test-breach-least-cost PRIVATE lsdtt3_cpu)
add_test(NAME unit_breach_least_cost COMMAND lsdtt-test-breach-least-cost)
lsdtt_set_tests_properties(unit_breach_least_cost PROPERTIES LABELS "flow")
# The same driver on the committed lakeland crop: float32 output under matched
# parameters equals Whitebox 2.0.6's (digests computed from Whitebox itself),
# and lsdtt3's default output routes (every interior cell descends; D8
# flowinfo validates). Reads tests/data only; no network, no skip path.
add_test(NAME unit_breach_least_cost_real_terrain
         COMMAND lsdtt-test-breach-least-cost --real-terrain
                 "${LSDTT_TEST_DATA_DIR}/flat_resolution_lakeland_crop_dem.tif")
lsdtt_set_tests_properties(unit_breach_least_cost_real_terrain PROPERTIES LABELS "flow")

add_test(NAME unit_flow_snapping COMMAND lsdtt-test-flow-snapping)
lsdtt_set_tests_properties(unit_flow_snapping PROPERTIES LABELS "flow;snapping")

add_test(NAME unit_mountain_front_outlet_pruning
         COMMAND lsdtt-test-mountain-front-outlet-pruning)
lsdtt_set_tests_properties(unit_mountain_front_outlet_pruning PROPERTIES
                     LABELS "flow;mountain_front")

# Depression hierarchy, phase 0: pins the
# dangling-reference defect in outlet_link -- copy semantics, outlet_db across a
# rehash, and the algorithm on a pit field dense enough to force one.
#
# NOTE: on unmodified HEAD this test FAILS BY DESIGN. It is written before the
# phase 1 fix so that the fix is what makes it pass. Do not "repair" it by
# weakening the assertions.
# Synthetic surfaces generated in process, so there is no fixture on disk.
add_executable(lsdtt-test-depression-hierarchy cpp/test_depression_hierarchy.cpp)
target_link_libraries(lsdtt-test-depression-hierarchy PRIVATE lsdtt3_cpu)

# RUN_SERIAL because the driver's complexity guard times two grid sizes against
# each other. Under `ctest -j4` the larger grid slows far more than the smaller
# one (memory contention from neighbouring tests), so the exponent drifted to
# the 1.6 bound: 1.609 on 2026-09-18 even with interleaved best-of-five timing,
# against 1.04 alone. The whole driver takes about a second.
add_test(NAME unit_depression_hierarchy
         COMMAND lsdtt-test-depression-hierarchy)
lsdtt_set_tests_properties(unit_depression_hierarchy PROPERTIES
    LABELS "depression_hierarchy"
    RUN_SERIAL TRUE)

# Oracle comparison against the upstream Barnes implementation. The expected
# files
# were generated ONCE from upstream and are committed; upstream is NOT built by
# lsdtt3 and is not needed to run this. Provenance, including the commit hashes
# pinned, is in tests/expected/depression_hierarchy/PROVENANCE.md.
add_test(NAME unit_depression_hierarchy_oracle
         COMMAND lsdtt-test-depression-hierarchy --oracle
                 "${LSDTT_TEST_DATA_DIR}/depression_hierarchy"
                 "${CMAKE_CURRENT_SOURCE_DIR}/expected/depression_hierarchy")
lsdtt_set_tests_properties(unit_depression_hierarchy_oracle PROPERTIES
    LABELS "depression_hierarchy")

# Internal drainage, phases 0, 1 and 3: the crater-field
# oracle, fill_pits_with_sinks, forced_outlet_flats, hole punching and sink
# tagging in both routing modes.
# Synthetic surfaces generated in process, so there is no fixture on disk.
add_executable(lsdtt-test-internal-drainage cpp/test_internal_drainage.cpp)
target_link_libraries(lsdtt-test-internal-drainage PRIVATE lsdtt3_cpu)

add_test(NAME unit_internal_drainage
         COMMAND lsdtt-test-internal-drainage)
lsdtt_set_tests_properties(unit_internal_drainage PROPERTIES
    LABELS "internal_drainage")

add_test(NAME fixture_internal_drainage_cli
         COMMAND lsdtt-test-internal-drainage --write-cli-fixture
                 ${LSDTT_INTERNAL_DRAINAGE_FIXTURE_DIR})
lsdtt_set_tests_properties(fixture_internal_drainage_cli PROPERTIES
    LABELS "fixture;internal_drainage"
    FIXTURES_SETUP internal_drainage_cli_fixture)

function(add_internal_drainage_cli_test TEST_NAME PARAM_TEMPLATE EXPECTED_OUTPUT)
    if(NOT LSDTT_ENABLE_PYTHON_TESTS)
        lsdtt_mark_test_unregistered(python ${TEST_NAME})
        return()
    endif()
    set(work_dir "${CMAKE_CURRENT_BINARY_DIR}/${TEST_NAME}")
    add_test(
        NAME ${TEST_NAME}
        COMMAND ${Python3_EXECUTABLE} "${LSDTT_TEST_SCRIPT}"
                --exe "$<TARGET_FILE:lsdtt-internal-drainage>"
                --param-template "${LSDTT_TEST_PARAM_DIR}/${PARAM_TEMPLATE}"
                --read-path "${LSDTT_INTERNAL_DRAINAGE_FIXTURE_DIR}"
                --work-dir "${work_dir}"
                --expect "${EXPECTED_OUTPUT}"
                ${ARGN})
    lsdtt_mark_harness_test(${TEST_NAME} ${ARGN})
    lsdtt_set_tests_properties(${TEST_NAME} PROPERTIES
        LABELS "cli;internal_drainage"
        FIXTURES_REQUIRED internal_drainage_cli_fixture)
endfunction()

add_internal_drainage_cli_test(
    smoke_internal_drainage_detect
    internal_drainage_detect.param.in
    phase4_detect_sink_candidates.fgb
    --expect phase4_detect_sink_candidates.csv
    --expect phase4_detect_sink_depressions.fgb
    --expect phase4_detect_sink_depressions.csv
    --expect phase4_detect_fill_diff.tif
    --expect-vector-feature-min phase4_detect_sink_candidates.fgb:3
    --stdout-contains "Stage: detect")

add_internal_drainage_cli_test(
    smoke_internal_drainage_both
    internal_drainage_both.param.in
    phase4_both_sink_id.tif
    --expect phase4_both_accepted_sinks.fgb
    --expect phase4_both_sink_stats.csv
    --expect phase4_both_sink_basins.fgb
    --expect phase4_both_sink_channels.fgb
    --expect phase4_both_sink_channel_nodes.csv
    --expect phase4_both_sink_junctions.fgb
    --expect phase4_both_sink_chi.tif
    --expect phase4_both_sink_drainage_area.tif
    --expect phase4_both_sink_flow_direction.tif
    --expect-file-contains "phase4_both_sink_stats.csv:routing_sink_id,sink_id"
    --expect-file-contains "phase4_both_sink_stats.csv:0,0"
    --expect-file-contains "phase4_both_sink_channel_nodes.csv:chi"
    --expect-vector-feature-min phase4_both_accepted_sinks.fgb:3
    --expect-vector-feature-min phase4_both_sink_basins.fgb:3
    --expect-vector-feature-min phase4_both_sink_channels.fgb:1
    --expect-raster-size phase4_both_sink_id.tif:300x300
    --expect-raster-stat phase4_both_sink_id.tif:min:0:0
    --expect-raster-stat phase4_both_sink_id.tif:max:2:2
    --expect-raster-stat phase4_both_sink_chi.tif:min:0:0
    --stdout-contains "Stage: route")

add_internal_drainage_cli_test(
    smoke_internal_drainage_route
    internal_drainage_route.param.in
    phase4_route_sink_id.tif
    --expect phase4_route_accepted_sinks.csv
    --expect phase4_route_sinks_punched.tif
    --expect phase4_route_sink_stats.csv
    --expect phase4_route_sink_channels.fgb
    --expect phase4_route_sink_channel_nodes.csv
    --expect phase4_route_sink_chi.tif
    --expect-file-contains "phase4_route_sink_stats.csv:0,0"
    --expect-file-contains "phase4_route_accepted_sinks.csv:1015.9000244140625"
    --expect-raster-size phase4_route_sink_id.tif:300x300
    --expect-raster-stat phase4_route_sink_id.tif:min:0:0
    --expect-raster-stat phase4_route_sink_id.tif:max:2:2
    --stdout-contains "Stage: route")

add_internal_drainage_cli_test(
    smoke_internal_drainage_moved_route
    internal_drainage_moved_route.param.in
    phase4_moved_route_accepted_sinks.csv
    --expect phase4_moved_route_sink_stats.csv
    --expect-file-contains "phase4_moved_route_accepted_sinks.csv:1033.9000244140625"
    --expect-file-contains "phase4_moved_route_accepted_sinks.csv:nan"
    --stdout-contains "Stage: route"
    --stderr-contains "'sink_candidate_format' is deprecated")

add_internal_drainage_cli_test(
    failure_internal_drainage_breach
    internal_drainage_breach_failure.param.in
    unused.tif
    --expect-failure
    --stderr-contains "breach_depressions_with_sinks")

add_internal_drainage_cli_test(
    failure_internal_drainage_png
    internal_drainage_png_failure.param.in
    unused.tif
    --expect-failure
    --stderr-contains "png is not supported")

# sink_leaves_only means nothing to fill_difference: each filled patch is one
# sink, with no nesting to take the leaves of. Asking for it must WARN and carry
# on, not fail and not be silently ignored.
add_internal_drainage_cli_test(
    smoke_internal_drainage_leaves_only_warning
    internal_drainage_leaves_only_warning.param.in
    phase4_leaves_warning_sink_candidates.fgb
    --stderr-contains "sink_leaves_only: false has no effect")

add_internal_drainage_cli_test(
    failure_internal_drainage_nested_both
    internal_drainage_nested_both_failure.param.in
    unused.tif
    --expect-failure
    --stderr-contains "both requires sink_leaves_only: true")

add_internal_drainage_cli_test(
    smoke_internal_drainage_empty_channels
    internal_drainage_empty_channels.param.in
    phase4_empty_channels_sink_channel_nodes.csv
    --expect-file-contains "phase4_empty_channels_sink_channel_nodes.csv:routing_sink_id,sink_id"
    --stderr-contains "FlatGeobuf cannot store an empty layer")

# Canonical-default regression. The parameter template deliberately omits
# threshold_contributing_pixels so this pins the shared default's behaviour.
add_internal_drainage_cli_test(
    regression_internal_drainage_default_threshold
    internal_drainage_default_threshold.param.in
    phase4_default_threshold_sink_channel_nodes.csv
    --expect-sha256 phase4_default_threshold_sink_channel_nodes.csv:02e7dd6cfb9fe505ad2f4083e166caa131f494d18e315dad2bafc866e6267a5e)

add_test(NAME smoke_internal_drainage_help
         COMMAND lsdtt-internal-drainage -h)
lsdtt_set_tests_properties(smoke_internal_drainage_help PROPERTIES
    LABELS "cli;internal_drainage"
    PASS_REGULAR_EXPRESSION "sink_leaves_only")

# Internal drainage, phase 2b: the detection module src/cpu/lsdtt_internal_drainage_cpu.  Covers both detection
# methods against the crater generator's analytic oracle, the binary-tree node
# count of a nested depression group, the shared pit cell of a meta-depression,
# the four size filters, the symmetric CSV / FlatGeobuf / GeoPackage review
# round-trip including a write with the FlatGeobuf spatial index ON (which
# reorders features), and the depression footprint polygons.  Kept
# separate from unit_internal_drainage because it writes vector fixtures and
# therefore needs an output directory; the two together are the internal
# drainage compartment.  Fixtures are written at run time, so nothing is
# committed.
add_executable(lsdtt-test-internal-drainage-detect
               cpp/test_internal_drainage_detect.cpp)
target_link_libraries(lsdtt-test-internal-drainage-detect PRIVATE lsdtt3_cpu)

add_test(NAME unit_internal_drainage_detect
         COMMAND lsdtt-test-internal-drainage-detect
                 ${CMAKE_CURRENT_BINARY_DIR}/internal_drainage_detect_output)
lsdtt_set_tests_properties(unit_internal_drainage_detect PROPERTIES
    LABELS "internal_drainage")
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/internal_drainage_detect_output)

# Internal drainage, phase 2a: the point-input attribute reader in src/core/lsdtt_point_input.  Guards the
# load_drainage_injection_points contract (message text, silent skipping,
# point_id counting skipped features) for its two live callers, cross-format
# equivalence of load_point_layer_with_attributes across .fgb / .gpkg / .csv,
# fields_missing yielding kind::absent rather than throwing, the as_bool
# spellings, and the deliberate divergences recorded in 6.4.  Fixtures are
# written into the test output directory at run time, so nothing is committed.
add_executable(lsdtt-test-point-input cpp/test_point_input.cpp)
target_link_libraries(lsdtt-test-point-input PRIVATE lsdtt3_cpu)

add_test(NAME unit_point_input
         COMMAND lsdtt-test-point-input
                 ${CMAKE_CURRENT_BINARY_DIR}/point_input_output)
lsdtt_set_tests_properties(unit_point_input PROPERTIES
    LABELS "point_input")
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/point_input_output)

# Characterisation tests for the src/utils + src/cpu consolidation (stage 1).
# Each embeds the private copies it replaced, verbatim, as the reference, and
# stays after those copies are deleted.
add_executable(lsdtt-test-raster-cell-predicates cpp/test_raster_cell_predicates.cpp)
target_link_libraries(lsdtt-test-raster-cell-predicates PRIVATE lsdtt3_cpu)
add_test(NAME unit_raster_cell_predicates COMMAND lsdtt-test-raster-cell-predicates)
lsdtt_set_tests_properties(unit_raster_cell_predicates PROPERTIES
    LABELS "raster")

add_executable(lsdtt-test-distributed-hydrology-argument-checks
    cpp/test_distributed_hydrology_argument_checks.cpp)
target_link_libraries(lsdtt-test-distributed-hydrology-argument-checks PRIVATE lsdtt3_cpu)
add_test(NAME unit_distributed_hydrology_argument_checks
         COMMAND lsdtt-test-distributed-hydrology-argument-checks)
lsdtt_set_tests_properties(unit_distributed_hydrology_argument_checks PROPERTIES
    LABELS "distributed_hydrology")

add_executable(lsdtt-test-drainage-injection-weights cpp/test_drainage_injection_weights.cpp)
target_link_libraries(lsdtt-test-drainage-injection-weights PRIVATE lsdtt3_cpu)
add_test(NAME unit_drainage_injection_weights
         COMMAND lsdtt-test-drainage-injection-weights
                 ${CMAKE_CURRENT_BINARY_DIR}/drainage_injection_weights_output)
lsdtt_set_tests_properties(unit_drainage_injection_weights PROPERTIES
    LABELS "point_input;flow")
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/drainage_injection_weights_output)

# Stage 3: bilinear raster sampling (C07) and point_in_ring (C22).
add_executable(lsdtt-test-raster-sampling cpp/test_raster_sampling.cpp)
target_link_libraries(lsdtt-test-raster-sampling PRIVATE lsdtt3_cpu)
add_test(NAME unit_raster_sampling COMMAND lsdtt-test-raster-sampling)
lsdtt_set_tests_properties(unit_raster_sampling PROPERTIES
    LABELS "raster")

# Stage 4: junction-graph topological order and s_vector preorder (C05).
add_executable(lsdtt-test-junction-graph-orders cpp/test_junction_graph_orders.cpp)
target_link_libraries(lsdtt-test-junction-graph-orders PRIVATE lsdtt3_cpu)
add_test(NAME unit_junction_graph_orders COMMAND lsdtt-test-junction-graph-orders)
lsdtt_set_tests_properties(unit_junction_graph_orders PROPERTIES
    LABELS "network")

# Stage 5: neighbour tables (C01) and sqrt(2) (C19). The table checks are
# static_asserts, so a mapping error fails the build of this test.
add_executable(lsdtt-test-grid-neighbours cpp/test_grid_neighbours.cpp)
target_link_libraries(lsdtt-test-grid-neighbours PRIVATE lsdtt3_cpu)
add_test(NAME unit_grid_neighbours COMMAND lsdtt-test-grid-neighbours)
lsdtt_set_tests_properties(unit_grid_neighbours PROPERTIES
    LABELS "raster")

# Stage 6: shared layer-to-raster CRS transform and spatial filter (C03).
add_executable(lsdtt-test-vector-layer-transform cpp/test_vector_layer_transform.cpp)
target_link_libraries(lsdtt-test-vector-layer-transform PRIVATE lsdtt3_cpu)
add_test(NAME unit_vector_layer_transform COMMAND lsdtt-test-vector-layer-transform)
lsdtt_set_tests_properties(unit_vector_layer_transform PROPERTIES
    LABELS "gdal")

# Stage 7: connected-component labelling (C08) and grid Dijkstra (C11).
# trace_labelled_polygons: written polygons read back and compared with their
# pixels (area, bounding box). Every file it wrote was empty until 2026-10-27.
add_executable(lsdtt-test-raster-vectorise cpp/test_raster_vectorise.cpp)
target_link_libraries(lsdtt-test-raster-vectorise PRIVATE lsdtt3_cpu)
add_test(NAME unit_raster_vectorise_polygons COMMAND lsdtt-test-raster-vectorise)
lsdtt_set_tests_properties(unit_raster_vectorise_polygons PROPERTIES
    LABELS "raster")

add_executable(lsdtt-test-connected-components cpp/test_connected_components.cpp)
target_link_libraries(lsdtt-test-connected-components PRIVATE lsdtt3_cpu)
add_test(NAME unit_connected_components COMMAND lsdtt-test-connected-components)
lsdtt_set_tests_properties(unit_connected_components PROPERTIES
    LABELS "raster")

add_executable(lsdtt-test-grid-geodesic cpp/test_grid_geodesic.cpp)
target_link_libraries(lsdtt-test-grid-geodesic PRIVATE lsdtt3_cpu)
add_test(NAME unit_grid_geodesic COMMAND lsdtt-test-grid-geodesic)
lsdtt_set_tests_properties(unit_grid_geodesic PROPERTIES
    LABELS "raster")

# Stage 8: boundary-condition string parsing (C09).
add_executable(lsdtt-test-boundary-conditions cpp/test_boundary_conditions.cpp)
target_link_libraries(lsdtt-test-boundary-conditions PRIVATE lsdtt3_cpu)
add_test(NAME unit_boundary_conditions COMMAND lsdtt-test-boundary-conditions)
lsdtt_set_tests_properties(unit_boundary_conditions PROPERTIES
    LABELS "flow")

# Stage 8: open-for-write helper (C14).
add_executable(lsdtt-test-output-file cpp/test_output_file.cpp)
add_test(NAME unit_output_file COMMAND lsdtt-test-output-file
         ${CMAKE_CURRENT_BINARY_DIR}/output_file_test)
lsdtt_set_tests_properties(unit_output_file PROPERTIES
    LABELS "core")
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/output_file_test)

# Stage 8: splitmix64 copies replaced by lsdtt_random.hpp (C15).
add_executable(lsdtt-test-splitmix-copies cpp/test_splitmix_copies.cpp)
add_test(NAME unit_splitmix_copies COMMAND lsdtt-test-splitmix-copies)
lsdtt_set_tests_properties(unit_splitmix_copies PROPERTIES
    LABELS "core")

# Stage 9: thread-count invariance of the kernels whose OpenMP thread setting
# moved from omp_set_num_threads to num_threads clauses (C16, decision D4).
add_executable(lsdtt-test-thread-count-invariance cpp/test_thread_count_invariance.cpp)
target_link_libraries(lsdtt-test-thread-count-invariance PRIVATE lsdtt3_cpu)
add_test(NAME unit_thread_count_invariance COMMAND lsdtt-test-thread-count-invariance)
lsdtt_set_tests_properties(unit_thread_count_invariance PROPERTIES
    LABELS "lem")

# Stage 10: landslide xorshift64* copies replaced by lsdtt_random (L3).
add_executable(lsdtt-test-landslide-rng cpp/test_landslide_rng.cpp)
add_test(NAME unit_landslide_rng COMMAND lsdtt-test-landslide-rng)
lsdtt_set_tests_properties(unit_landslide_rng PROPERTIES
    LABELS "landslide")

# Stage 10: thread setting of the hillslope solvers moved from
# omp_set_num_threads to num_threads clauses (C16-H, decision D4).
add_executable(lsdtt-test-hillslope-thread-count-invariance cpp/test_hillslope_thread_count_invariance.cpp)
target_link_libraries(lsdtt-test-hillslope-thread-count-invariance PRIVATE lsdtt3_cpu)
add_test(NAME unit_hillslope_thread_count_invariance COMMAND lsdtt-test-hillslope-thread-count-invariance)
lsdtt_set_tests_properties(unit_hillslope_thread_count_invariance PROPERTIES
    LABELS "hillslope")

# Stage 10: bilinear_regular takes its axis weights from regular_axis_weight (L1).
add_executable(lsdtt-test-cosmo-lsdn-grid-lookup cpp/test_cosmo_lsdn_grid_lookup.cpp)
target_link_libraries(lsdtt-test-cosmo-lsdn-grid-lookup PRIVATE lsdtt3_cpu)
add_test(NAME unit_cosmo_lsdn_grid_lookup COMMAND lsdtt-test-cosmo-lsdn-grid-lookup)
lsdtt_set_tests_properties(unit_cosmo_lsdn_grid_lookup PROPERTIES
    LABELS "cosmo")

# The production / pressure rasters' WGS84 conversion must not depend on the
# axis order the DEM's CRS declares (northing-first CRSs such as NZTM).
add_executable(lsdtt-test-cosmo-raster-axis-order cpp/test_cosmo_raster_axis_order.cpp)
target_link_libraries(lsdtt-test-cosmo-raster-axis-order PRIVATE lsdtt3_cpu)
add_test(NAME unit_cosmo_raster_axis_order COMMAND lsdtt-test-cosmo-raster-axis-order)
lsdtt_set_tests_properties(unit_cosmo_raster_axis_order PROPERTIES
    LABELS "cosmo")

# Stage 11: four-cardinal centred-difference gradient shared by the
# vertical-units check and the LEM timestep diagnosis (C18).
add_executable(lsdtt-test-raster-central-difference cpp/test_raster_central_difference.cpp)
target_link_libraries(lsdtt-test-raster-central-difference PRIVATE lsdtt3_cpu)
add_test(NAME unit_raster_central_difference COMMAND lsdtt-test-raster-central-difference)
lsdtt_set_tests_properties(unit_raster_central_difference PROPERTIES
    LABELS "raster")

# Stage 11: Manning's equation shared by hydro_velocity_field and
# manning_hydraulics (C20).
add_executable(lsdtt-test-hydro-velocity-manning-formula cpp/test_hydro_velocity_manning_formula.cpp)
target_link_libraries(lsdtt-test-hydro-velocity-manning-formula PRIVATE lsdtt3_cpu)
add_test(NAME unit_hydro_velocity_manning_formula COMMAND lsdtt-test-hydro-velocity-manning-formula)
lsdtt_set_tests_properties(unit_hydro_velocity_manning_formula PROPERTIES
    LABELS "network")

# Stage 2: the in-place medians and quantiles in src/utils/lsdtt_statistics.
add_executable(lsdtt-test-order-statistics cpp/test_order_statistics.cpp)
target_link_libraries(lsdtt-test-order-statistics PRIVATE lsdtt3_cpu)
add_test(NAME unit_order_statistics COMMAND lsdtt-test-order-statistics)
lsdtt_set_tests_properties(unit_order_statistics PROPERTIES
    LABELS "statistics")

# Stage 2 part B: records the D7 change from the private ols_fit to
# compute_linear_regression in the drainage-area exceedance scaling.
add_executable(lsdtt-test-drainage-area-exceedance-ols cpp/test_drainage_area_exceedance_ols.cpp)
target_link_libraries(lsdtt-test-drainage-area-exceedance-ols PRIVATE lsdtt3_cpu)
add_test(NAME unit_drainage_area_exceedance_ols COMMAND lsdtt-test-drainage-area-exceedance-ols)
lsdtt_set_tests_properties(unit_drainage_area_exceedance_ols PROPERTIES
    LABELS "statistics")

add_executable(lsdtt-test-hydraulic-routing cpp/test_hydraulic_routing.cpp)
target_link_libraries(lsdtt-test-hydraulic-routing PRIVATE lsdtt3_cpu)

add_test(NAME unit_hydraulic_routing
         COMMAND lsdtt-test-hydraulic-routing)
lsdtt_set_tests_properties(unit_hydraulic_routing PROPERTIES LABELS "hydraulic")

add_executable(lsdtt-test-raster-resampling cpp/test_raster_resampling.cpp)
target_link_libraries(lsdtt-test-raster-resampling PRIVATE lsdtt3_cpu)

add_test(NAME unit_raster_resampling
         COMMAND lsdtt-test-raster-resampling)
lsdtt_set_tests_properties(unit_raster_resampling PROPERTIES LABELS "raster")

# Oblique analysis frames: geometry, axis convention, datum, line and polygon
# builders, checked against an independent Vincenty geodesic.
add_executable(lsdtt-test-oblique-frame cpp/test_oblique_frame.cpp)
target_link_libraries(lsdtt-test-oblique-frame PRIVATE lsdtt3_cpu)

add_test(NAME unit_oblique_frame
         COMMAND lsdtt-test-oblique-frame)
lsdtt_set_tests_properties(unit_oblique_frame PROPERTIES LABELS "raster")

# Feature files that define an oblique frame (line, two points, polygon) and
# every refusal; fixtures are GeoPackages written by the test itself.
add_executable(lsdtt-test-oblique-frame-feature cpp/test_oblique_frame_feature.cpp)
target_link_libraries(lsdtt-test-oblique-frame-feature PRIVATE lsdtt3_cpu)

add_test(NAME unit_oblique_frame_feature
         COMMAND lsdtt-test-oblique-frame-feature
                 "${CMAKE_CURRENT_BINARY_DIR}/oblique_frame_feature_work")
lsdtt_set_tests_properties(unit_oblique_frame_feature PROPERTIES LABELS "raster")

# The rotation frame file: exact write/read round trip, YAML-safe writing
# rules, consistency checks and refusals, and the frame_from_crs_wkt fallback.
add_executable(lsdtt-test-oblique-frame-file cpp/test_oblique_frame_file.cpp)
target_link_libraries(lsdtt-test-oblique-frame-file PRIVATE lsdtt3_cpu)

add_test(NAME unit_oblique_frame_file
         COMMAND lsdtt-test-oblique-frame-file
                 "${CMAKE_CURRENT_BINARY_DIR}/oblique_frame_file_work")
lsdtt_set_tests_properties(unit_oblique_frame_file PROPERTIES LABELS "raster")

# Warping a DEM into an oblique frame and rotating products back onto the
# source or -tap lattice; fixtures are GeoTIFFs written by the test itself.
add_executable(lsdtt-test-oblique-frame-warp cpp/test_oblique_frame_warp.cpp)
target_link_libraries(lsdtt-test-oblique-frame-warp PRIVATE lsdtt3_cpu)

add_test(NAME unit_oblique_frame_warp
         COMMAND lsdtt-test-oblique-frame-warp
                 "${CMAKE_CURRENT_BINARY_DIR}/oblique_frame_warp_work")
lsdtt_set_tests_properties(unit_oblique_frame_warp PROPERTIES LABELS "raster")

# Rotating an existing DEM into an oblique frame as lsdtt-raster-preprocessing
# does (line, polygon, valid-data footprint; NAD83 source) and rotating its
# products back onto the source pixels; fixtures are written by the test.
add_executable(lsdtt-test-oblique-frame-preprocess cpp/test_oblique_frame_preprocess.cpp)
target_link_libraries(lsdtt-test-oblique-frame-preprocess PRIVATE lsdtt3_cpu)

add_test(NAME unit_oblique_frame_preprocess
         COMMAND lsdtt-test-oblique-frame-preprocess
                 "${CMAKE_CURRENT_BINARY_DIR}/oblique_frame_preprocess_work")
lsdtt_set_tests_properties(unit_oblique_frame_preprocess PROPERTIES LABELS "raster")

# Grid convergence and bearings_relative_to_true_north: convergence against
# pyproj, wrapping, a plane dipping to a known true azimuth in UTM and in
# frames at several azimuths, raster lattice interpolation, start-up message.
add_executable(lsdtt-test-true-north cpp/test_true_north.cpp)
target_link_libraries(lsdtt-test-true-north PRIVATE lsdtt3_cpu)

add_test(NAME unit_true_north COMMAND lsdtt-test-true-north)
lsdtt_set_tests_properties(unit_true_north PROPERTIES LABELS "raster")

# CRS matching by meaning, not WKT text: an lsdtt3-written DEM and a gdalwarp
# mask in one oblique Mercator frame must be accepted together.
add_executable(lsdtt-test-crs-semantic-matching cpp/test_crs_semantic_matching.cpp)
target_link_libraries(lsdtt-test-crs-semantic-matching PRIVATE lsdtt3_cpu)

add_test(NAME unit_crs_semantic_matching
         COMMAND lsdtt-test-crs-semantic-matching)
lsdtt_set_tests_properties(unit_crs_semantic_matching PROPERTIES LABELS "raster")

add_executable(lsdtt-test-raster-hierarchy cpp/test_raster_hierarchy.cpp)
target_link_libraries(lsdtt-test-raster-hierarchy PRIVATE lsdtt3_cpu)

add_test(NAME unit_raster_hierarchy
         COMMAND lsdtt-test-raster-hierarchy)
lsdtt_set_tests_properties(unit_raster_hierarchy PROPERTIES LABELS "raster")

add_executable(lsdtt-test-raster-jenks-classification
               cpp/test_raster_jenks_classification.cpp)
target_link_libraries(lsdtt-test-raster-jenks-classification PRIVATE lsdtt3_cpu)

add_test(NAME unit_raster_jenks_classification
         COMMAND lsdtt-test-raster-jenks-classification)
lsdtt_set_tests_properties(unit_raster_jenks_classification PROPERTIES
    LABELS "raster;jenks")

add_executable(lsdtt-test-raster-active-masks cpp/test_raster_active_masks.cpp)
target_link_libraries(lsdtt-test-raster-active-masks PRIVATE lsdtt3_cpu)

add_test(NAME unit_raster_active_masks COMMAND lsdtt-test-raster-active-masks)
lsdtt_set_tests_properties(unit_raster_active_masks PROPERTIES
    LABELS "raster;mask")

# Divide-migration compute layer.  Covers the routines in src/cpu/lsdtt_divide_migration_cpu.cpp:
# nearest-divide geometry (mid-segment projection, zero-distance guard, negative
# relief, nodata sampling, and the PER-UNIT index that stops a head being matched
# to a divide between two other units), the common-base-level chi rebasing,
# fragment merging, the pair summaries (which are driven by the DIVIDES, so a
# pair with no matched heads still gets a row), and the eight-connected basin
# adjacency pass that decides which pairs are worth a ridgeline extraction.
#
# The same binary writes the synthetic three-valley DEM and outlet layer used by
# the divide-migration CLI end-to-end tests, so the fixture and the in-memory
# geometry come from one generator.  Test domains are generated in C++, never in
# Python (docs/developer_notes/testing_and_benchmarking.md §8).
add_executable(lsdtt-test-divide-migration cpp/test_divide_migration.cpp)
target_link_libraries(lsdtt-test-divide-migration PRIVATE lsdtt3_cpu)

add_executable(lsdtt-test-divide-diagnostics cpp/test_divide_diagnostics.cpp)
target_link_libraries(lsdtt-test-divide-diagnostics PRIVATE lsdtt3_cpu)

add_executable(lsdtt-test-paired-flowpaths cpp/test_paired_flowpaths.cpp)
target_link_libraries(lsdtt-test-paired-flowpaths PRIVATE lsdtt3_cpu)

add_executable(lsdtt-test-landscape-equilibration cpp/test_landscape_equilibration.cpp)
target_link_libraries(lsdtt-test-landscape-equilibration PRIVATE lsdtt3_cpu)

# The factor-of-two ambiguity in "mean hillslope length":
# ridge_mean is L where the all-pixel mean is L/2.
add_executable(lsdtt-test-hillslope-length-statistic cpp/test_hillslope_length_statistic.cpp)
target_link_libraries(lsdtt-test-hillslope-length-statistic PRIVATE lsdtt3_cpu)

add_test(NAME unit_hillslope_length_statistic
         COMMAND lsdtt-test-hillslope-length-statistic)
lsdtt_set_tests_properties(unit_hillslope_length_statistic PROPERTIES
    LABELS "hillslope")

add_executable(lsdtt-test-critical-slope cpp/test_critical_slope.cpp)
target_link_libraries(lsdtt-test-critical-slope PRIVATE lsdtt3_cpu)

add_test(NAME unit_divide_diagnostics
         COMMAND lsdtt-test-divide-diagnostics)
lsdtt_set_tests_properties(unit_divide_diagnostics PROPERTIES
    LABELS "divide_diagnostics")

add_test(NAME unit_paired_flowpaths
         COMMAND lsdtt-test-paired-flowpaths
                 "${LSDTT_PAIRED_FLOWPATHS_FIXTURE_DIR}")
lsdtt_set_tests_properties(unit_paired_flowpaths PROPERTIES
    LABELS "paired_flowpaths"
    FIXTURES_SETUP paired_flowpaths_fixture)

# =============================================================================
# lsdtt-ridgeline-analysis -- paired flow paths through the CLI
#
# There was NO end-to-end coverage of this CLI at all.  The unit cases in
# lsdtt-test-paired-flowpaths exercise the kernel only, and the failure mode
# that silently changed what the network-geometry sweeps measured -- a
# parameter READ via bf[...] but never DECLARED in bool_default_map, which
# default-constructs to false with no parse error and no warning -- is
# invisible to them because it lives in the CLI's parameter plumbing.  These
# run the binary.
#
# The straight-link pair is the switch assertion: the summary CSV must gain a
# non-NaN median_delta_L_straight_link when the switch is on and keep NaN when
# it is off.  Each assertion fails on the other arm's output, so the pair
# discriminates rather than merely passing.
# =============================================================================
add_lsdtt_cli_smoke_test(
    smoke_ridgeline_paired_flowpaths_sl_on
    lsdtt-ridgeline-analysis
    ridgeline_paired_flowpaths_sl_on.param.in
    pfslon_paired_flowpaths_summary.csv
    --read-path "${LSDTT_PAIRED_FLOWPATHS_FIXTURE_DIR}"
    --expect pfslon_paired_flowpaths.csv
    --expect-file-contains "pfslon_paired_flowpaths_summary.csv:,870,0."
    --stdout-contains "straight-link pairs = 870"
    --stdout-contains "rejected, paths never rejoin:        0"
)
lsdtt_set_tests_properties(smoke_ridgeline_paired_flowpaths_sl_on PROPERTIES
    LABELS "cli;paired_flowpaths"
    FIXTURES_REQUIRED paired_flowpaths_fixture)

add_lsdtt_cli_smoke_test(
    rigorous_ridgeline_legacy_threshold
    lsdtt-ridgeline-analysis
    ridgeline_legacy_threshold.param.in
    ridge_legacy_paired_flowpaths_summary.csv
    --read-path "${LSDTT_PAIRED_FLOWPATHS_FIXTURE_DIR}"
    --stdout-contains "straight-link pairs = 870"
    --stderr-contains "'source_pixels_threshold' is deprecated"
    --stderr-contains "use 'threshold_contributing_pixels'"
    --stderr-contains "'print_paired_flowpath_statistics' is deprecated"
    --stderr-contains "use 'write_paired_flowpath_statistics'"
    --expect-file-contains "ridge_legacy_ingested_params.param:write_paired_flowpath_statistics: true"
)
lsdtt_set_tests_properties(rigorous_ridgeline_legacy_threshold PROPERTIES
    FIXTURES_REQUIRED paired_flowpaths_fixture)

add_lsdtt_cli_smoke_test(
    smoke_ridgeline_paired_flowpaths_sl_off
    lsdtt-ridgeline-analysis
    ridgeline_paired_flowpaths_sl_off.param.in
    pfsloff_paired_flowpaths_summary.csv
    --read-path "${LSDTT_PAIRED_FLOWPATHS_FIXTURE_DIR}"
    --expect-file-contains "pfsloff_paired_flowpaths_summary.csv:,0,NaN,NaN,0"
    --stdout-contains "rejected, paths never rejoin:        0"
)
lsdtt_set_tests_properties(smoke_ridgeline_paired_flowpaths_sl_off PROPERTIES
    LABELS "cli;paired_flowpaths"
    FIXTURES_REQUIRED paired_flowpaths_fixture)

# The contrast that makes the two zeroes above mean something: the same fixture
# selected the OLD way (all junction basins, all pairs) rejects 32291 of 38852
# divide edges as non-rejoining.  Without this row, "0 rejected" could be an
# artefact of a fixture that simply has no non-sibling pairs to reject.
add_lsdtt_cli_smoke_test(
    smoke_ridgeline_paired_flowpaths_all_basins
    lsdtt-ridgeline-analysis
    ridgeline_paired_flowpaths_all_basins.param.in
    pfsall_paired_flowpaths_summary.csv
    --read-path "${LSDTT_PAIRED_FLOWPATHS_FIXTURE_DIR}"
    --stdout-contains "rejected, paths never rejoin:        32291"
    --expect-file-contains "pfsall_paired_flowpaths_summary.csv:6561,0."
)
lsdtt_set_tests_properties(smoke_ridgeline_paired_flowpaths_all_basins PROPERTIES
    LABELS "cli;paired_flowpaths"
    FIXTURES_REQUIRED paired_flowpaths_fixture)

add_test(NAME unit_divide_migration
         COMMAND lsdtt-test-divide-migration
                 "${LSDTT_DIVIDE_MIGRATION_FIXTURE_DIR}")
lsdtt_set_tests_properties(unit_divide_migration PROPERTIES
    LABELS "divide_migration"
    FIXTURES_SETUP divide_migration_fixture)

add_test(NAME unit_landscape_equilibration
         COMMAND lsdtt-test-landscape-equilibration)
lsdtt_set_tests_properties(unit_landscape_equilibration PROPERTIES
    LABELS "landscape_equilibration")

add_test(NAME unit_hillslope_critical_slope
         COMMAND lsdtt-test-critical-slope)
lsdtt_set_tests_properties(unit_hillslope_critical_slope PROPERTIES
    LABELS "hillslope")

# Kinematic streamline tracing (Lea 1992 + Tarboton 1997 edge routing).
# Phase 0 of the hillslope-channel coupling work: the module shipped as a draft
# with no tests and no caller, and the hilltop-flow-routing layer is to be built
# on it, so this driver characterises it before anything depends on it.
# Analytic fixtures only (tilted planes with an exact central-difference aspect,
# and a cone), so there is no committed DEM and it runs in the default suite.
add_executable(lsdtt-test-streamlines cpp/test_streamlines.cpp)
target_link_libraries(lsdtt-test-streamlines PRIVATE lsdtt3_cpu)

add_test(NAME unit_streamlines
         COMMAND lsdtt-test-streamlines)
lsdtt_set_tests_properties(unit_streamlines PROPERTIES
    LABELS "streamlines")

# Ridge extraction for hilltop flow routing (Phase 1). The divide-grid corner
# convention is the one claim everything downstream rests on, so the fixtures
# are hand-drawn 5x5 grids with the expected pixels written out explicitly,
# including a periodic-seam case and the non-4-connected step.
add_executable(lsdtt-test-hilltop-flow-routing cpp/test_hilltop_flow_routing.cpp)
target_link_libraries(lsdtt-test-hilltop-flow-routing PRIVATE lsdtt3_cpu)

add_test(NAME unit_hilltop_flow_routing
         COMMAND lsdtt-test-hilltop-flow-routing)
lsdtt_set_tests_properties(unit_hilltop_flow_routing PROPERTIES
    LABELS "hilltop_flow_routing")

# The integration half of Phase 1: run ridge extraction against REAL basin
# outlines rather than hand-made ridgelines, and MEASURE what fraction of real
# outline steps are 4-connected. If real outlines were full of diagonal steps
# the edge-based derivation would be the wrong model, so this is the check that
# decides whether the approach is sound.
add_executable(lsdtt-test-hilltop-integration cpp/test_hilltop_integration.cpp)
target_link_libraries(lsdtt-test-hilltop-integration PRIVATE lsdtt3_cpu)

add_test(NAME unit_hilltop_flow_routing_integration
         COMMAND lsdtt-test-hilltop-integration)
lsdtt_set_tests_properties(unit_hilltop_flow_routing_integration PROPERTIES
    LABELS "hilltop_flow_routing")

# --------------------------------------------------------------------------
# lsdtt-hillslope-channel-coupling
#
# The CLI runs on the SAME ridge-and-valley surface the integration driver
# asserts against, so the tool is exercised on a landscape whose behaviour is
# already characterised rather than on an unrelated fixture.
# --------------------------------------------------------------------------
file(MAKE_DIRECTORY ${LSDTT_HILLTOP_COUPLING_FIXTURE_DIR})

add_test(NAME fixture_hilltop_coupling_cli
         COMMAND lsdtt-test-hilltop-integration --write-cli-fixture
                 ${LSDTT_HILLTOP_COUPLING_FIXTURE_DIR})
lsdtt_set_tests_properties(fixture_hilltop_coupling_cli PROPERTIES
    LABELS "fixture;hilltop_flow_routing"
    FIXTURES_SETUP hilltop_coupling_cli_fixture)

function(add_hilltop_coupling_cli_test TEST_NAME PARAM_TEMPLATE EXPECTED_OUTPUT)
    if(NOT LSDTT_ENABLE_PYTHON_TESTS)
        lsdtt_mark_test_unregistered(python ${TEST_NAME})
        return()
    endif()
    set(work_dir "${CMAKE_CURRENT_BINARY_DIR}/${TEST_NAME}")
    add_test(
        NAME ${TEST_NAME}
        COMMAND ${Python3_EXECUTABLE} "${LSDTT_TEST_SCRIPT}"
                --exe "$<TARGET_FILE:lsdtt-hillslope-channel-coupling>"
                --param-template "${LSDTT_TEST_PARAM_DIR}/${PARAM_TEMPLATE}"
                --read-path "${LSDTT_HILLTOP_COUPLING_FIXTURE_DIR}"
                --work-dir "${work_dir}"
                --expect "${EXPECTED_OUTPUT}"
                ${ARGN})
    lsdtt_mark_harness_test(${TEST_NAME} ${ARGN})
    lsdtt_set_tests_properties(${TEST_NAME} PROPERTIES
        LABELS "cli;hilltop_flow_routing"
        FIXTURES_REQUIRED hilltop_coupling_cli_fixture)
endfunction()

add_hilltop_coupling_cli_test(
    smoke_hillslope_channel_coupling
    hillslope_channel_coupling_smoke.param.in
    hcc_hilltop_points.csv
    --expect hcc_ridgelines.fgb
    --stdout-contains "never traced (DEM border)")

# The simple floodplain detector, pinned under its 2026-10-26 names and under
# the old ones. Phase 4 of the fluvial-surface consolidation moves it onto the
# shared floodplain engine; the mask must not move.
add_hilltop_coupling_cli_test(
    rigorous_hillslope_channel_coupling_floodplain
    hillslope_channel_coupling_floodplain.param.in
    hccfp_floodplain_mask.tif
    --expect-raster-checksum hccfp_floodplain_mask.tif:448)

add_hilltop_coupling_cli_test(
    rigorous_hillslope_channel_coupling_floodplain_legacy_keys
    hillslope_channel_coupling_floodplain_legacy_keys.param.in
    hccfp_floodplain_mask.tif
    --expect-raster-checksum hccfp_floodplain_mask.tif:448
    --stderr-contains "'floodplain_hand_threshold' is deprecated"
    --stderr-contains "use 'max_floodplain_hand_metres'"
    --stderr-contains "'floodplain_gradient_threshold' is deprecated"
    --stderr-contains "use 'max_floodplain_gradient'")

# The floodplain bounds chosen from the data through the shared threshold
# picker: CDF divergence finds both and the chosen values are recorded; the QQ
# rule finds no HAND bound on this small DEM and must stop, not return an
# empty floodplain.
add_hilltop_coupling_cli_test(
    rigorous_hillslope_channel_coupling_floodplain_cdf
    hillslope_channel_coupling_floodplain_cdf.param.in
    hcccdf_floodplain_mask.tif
    --stdout-contains "Floodplain bounds (cdf_divergence): HAND < 5.81376123 m, gradient < 0.452985913"
    --expect-raster-checksum hcccdf_floodplain_mask.tif:1208
    --expect-file-contains "hcccdf_ingested_params.param:max_floodplain_hand_metres: 5.81376123")

add_hilltop_coupling_cli_test(
    failure_hillslope_channel_coupling_floodplain_qq_not_found
    hillslope_channel_coupling_floodplain_qq.param.in
    unused.tif
    --expect-failure
    --stderr-contains "found no HAND bound")

add_hilltop_coupling_cli_test(
    rigorous_hillslope_channel_coupling_legacy_threshold
    hillslope_channel_coupling_legacy_threshold.param.in
    hcc_legacy_hilltop_points.csv
    --expect hcc_legacy_ridgelines.fgb
    --stdout-contains "Channel pixels: 423"
    --stderr-contains "'source_pixels_threshold' is deprecated"
    --stderr-contains "use 'threshold_contributing_pixels'")

# E* and R* must be undefined, not zero, where the trace never ran.
add_hilltop_coupling_cli_test(
    rigorous_hillslope_channel_coupling_null_handling
    hillslope_channel_coupling_smoke.param.in
    hcc_hilltop_points.csv
    --expect-file-not-contains "hcc_hilltop_points.csv:nan"
    --expect-file-not-contains "hcc_hilltop_points.csv:-9999")

# Phase 5: k_sn must reach the hilltop table, and must stay EMPTY (not zero)
# on rows whose trace never found a channel.
add_hilltop_coupling_cli_test(
    smoke_hillslope_channel_coupling_chi
    hillslope_channel_coupling_chi.param.in
    hccchi_hilltop_points.csv
    --expect-file-contains "hccchi_hilltop_points.csv:channel_k_sn"
    --expect-file-not-contains "hccchi_hilltop_points.csv:nan"
    --expect-file-not-contains "hccchi_hilltop_points.csv:-9999"
    --stdout-contains "Joined k_sn onto")

add_hilltop_coupling_cli_test(
    failure_hillslope_channel_coupling_bad_sc
    hillslope_channel_coupling_bad_sc.param.in
    ""
    --expect-failure
    --stderr-contains "critical_hillslope_slope must be > 0")

add_hilltop_coupling_cli_test(
    failure_hillslope_channel_coupling_bad_ridge_choice
    hillslope_channel_coupling_bad_ridge_choice.param.in
    ""
    --expect-failure
    --stderr-contains "ridge_extraction_algorithm_choice must be")

add_hilltop_coupling_cli_test(
    failure_hillslope_channel_coupling_no_channels
    hillslope_channel_coupling_no_channels.param.in
    ""
    --expect-failure
    --stderr-contains "no channel")

# SCAFFOLDING for the E*-R* validation workflow: runs the
# hilltop pipeline on a supplied DEM. Not a test - it is built but never
# registered with CTest, and it exists only until the Phase 4 CLI replaces it.
# What incision_threshold_area is compared against: specific area in METRES
# under divide_area_by_dx, map area in m^2 without it.
add_executable(lsdtt-test-channel-threshold-units cpp/test_channel_threshold_units.cpp)
target_link_libraries(lsdtt-test-channel-threshold-units PRIVATE lsdtt3_cpu)

add_test(NAME unit_channel_threshold_units
         COMMAND lsdtt-test-channel-threshold-units)
lsdtt_set_tests_properties(unit_channel_threshold_units PROPERTIES
    LABELS "lem")

# Pathway A: the channel
# threshold search, and the ways its objective can be plausibly wrong.
add_executable(lsdtt-test-hillslope-length-search cpp/test_hillslope_length_search.cpp)
target_link_libraries(lsdtt-test-hillslope-length-search PRIVATE lsdtt3_cpu)

add_test(NAME unit_hillslope_length_search
         COMMAND lsdtt-test-hillslope-length-search)
lsdtt_set_tests_properties(unit_hillslope_length_search PROPERTIES
    LABELS "hillslope")

# Phase 5 sweep: every vector-ingesting CLI against multipart input.
# A SCRIPT rather than a C++ test because it drives the shipped binaries -- the
# defect it guards was invisible at the library level and only showed up in what
# a CLI actually wrote.
# Needs GeoPandas, NumPy, Rasterio and Shapely, so it is a geospatial Python
# test. (It used to sit behind its own quiet find_package(Python3), which meant
# a native-only configure still went looking for Python.)
if(LSDTT_ENABLE_PYTHON_GEO_TESTS)
    file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/multipart_sweep_output)
    add_test(NAME smoke_multipart_cli_sweep
             COMMAND ${Python3_EXECUTABLE}
                     ${CMAKE_CURRENT_SOURCE_DIR}/scripts/multipart_sweep.py
                     $<TARGET_FILE_DIR:lsdtt-linear-stats>
                     ${CMAKE_CURRENT_BINARY_DIR}/multipart_sweep_output)
    lsdtt_mark_test_runtime(python smoke_multipart_cli_sweep)
else()
    lsdtt_mark_test_unregistered(python_geo smoke_multipart_cli_sweep)
endif()
lsdtt_set_tests_properties(smoke_multipart_cli_sweep PROPERTIES
    LABELS "cli;raster" TIMEOUT 1200)

# Longitudinal profile from the points output.
add_executable(lsdtt-test-channel-profile-attributes cpp/test_channel_profile_attributes.cpp)
target_link_libraries(lsdtt-test-channel-profile-attributes PRIVATE lsdtt3_cpu)

file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/profile_output)
add_test(NAME unit_channel_profile_attributes
         COMMAND lsdtt-test-channel-profile-attributes
                 ${CMAKE_CURRENT_BINARY_DIR}/profile_output)
lsdtt_set_tests_properties(unit_channel_profile_attributes PROPERTIES LABELS "network")

# link_length under periodic boundaries.
add_executable(lsdtt-test-link-length-periodic cpp/test_link_length_periodic.cpp)
target_link_libraries(lsdtt-test-link-length-periodic PRIVATE lsdtt3_cpu)

file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/link_length_output)
add_test(NAME unit_link_length_periodic
         COMMAND lsdtt-test-link-length-periodic
                 ${CMAKE_CURRENT_BINARY_DIR}/link_length_output)
lsdtt_set_tests_properties(unit_link_length_periodic PROPERTIES LABELS "network")

# Multipart geometry on read.
add_executable(lsdtt-test-multipart-geometry cpp/test_multipart_geometry.cpp)
target_link_libraries(lsdtt-test-multipart-geometry PRIVATE lsdtt3_core)

file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/multipart_output)
add_test(NAME unit_multipart_geometry
         COMMAND lsdtt-test-multipart-geometry
                 ${CMAKE_CURRENT_BINARY_DIR}/multipart_output)
lsdtt_set_tests_properties(unit_multipart_geometry PROPERTIES LABELS "raster")

# Pathway A from the command line, while it is not yet a lsdtt-lem parameter.
add_executable(find-threshold-for-hillslope-length cpp/find_threshold_for_hillslope_length_driver.cpp)
target_link_libraries(find-threshold-for-hillslope-length PRIVATE lsdtt3_cpu)

# Phase 6: hillslope length on a COUPLED landscape, where there is no channel
# threshold. Valleys come from the curvature inflection instead.
add_executable(emergent-hillslope-length-driver cpp/emergent_hillslope_length_driver.cpp)
target_link_libraries(emergent-hillslope-length-driver PRIVATE lsdtt3_cpu)

# Pathway B: K and D from a target relief and hillslope length.
add_executable(lsdtt-test-hillslope-peclet cpp/test_hillslope_peclet.cpp)
target_link_libraries(lsdtt-test-hillslope-peclet PRIVATE lsdtt3_cpu)

add_test(NAME unit_hillslope_peclet COMMAND lsdtt-test-hillslope-peclet)
lsdtt_set_tests_properties(unit_hillslope_peclet PROPERTIES LABELS "hillslope")

# Boundary conditions for the explicit Roering transport solver (periodic + 'o').
add_executable(lsdtt-test-hillslope-transport-bc cpp/test_hillslope_transport_bc.cpp)
target_link_libraries(lsdtt-test-hillslope-transport-bc PRIVATE lsdtt3_cpu)

add_test(NAME unit_hillslope_transport_bc COMMAND lsdtt-test-hillslope-transport-bc)
lsdtt_set_tests_properties(unit_hillslope_transport_bc PROPERTIES LABELS "hillslope")

# Coupled-LEM timestep control: an accuracy limit that no stability check catches.
add_executable(lsdtt-test-lem-timestep cpp/test_lem_timestep.cpp)
target_link_libraries(lsdtt-test-lem-timestep PRIVATE lsdtt3_cpu)

add_test(NAME unit_lem_timestep COMMAND lsdtt-test-lem-timestep)
lsdtt_set_tests_properties(unit_lem_timestep PROPERTIES LABELS "lem")

# Process-domain coupling: every hillslope model must hold the channel pixels at
# their fluvial elevations and relax the hillslopes toward them, not toward the
# outlet. adi and cidre used to miss the channel set.
add_executable(lsdtt-test-lem-channel-coupling cpp/test_lem_channel_coupling.cpp)
target_link_libraries(lsdtt-test-lem-channel-coupling PRIVATE lsdtt3_cpu)

add_test(NAME unit_lem_channel_coupling COMMAND lsdtt-test-lem-channel-coupling)

# Forcing history: optional erodibility / hillslope_diffusivity columns, whose
# absent value must mean "keep the parameter file", never 0.
add_executable(lsdtt-test-lem-forcing-history cpp/test_lem_forcing_history.cpp)
target_link_libraries(lsdtt-test-lem-forcing-history PRIVATE lsdtt3_cpu)
add_test(NAME unit_lem_forcing_history
         COMMAND lsdtt-test-lem-forcing-history "${CMAKE_CURRENT_BINARY_DIR}")
lsdtt_set_tests_properties(unit_lem_forcing_history PROPERTIES LABELS "lem")

# Lithology maps to erodibility: property table, crosswalk, the categorical
# polygon rasteriser and the K raster (lsdtt-lithology-erodibility's library).
add_executable(lsdtt-test-lithology-erodibility cpp/test_lithology_erodibility.cpp)
target_link_libraries(lsdtt-test-lithology-erodibility PRIVATE lsdtt3_cpu)
add_test(NAME unit_lem_lithology_erodibility
         COMMAND lsdtt-test-lithology-erodibility
                 "${CMAKE_CURRENT_BINARY_DIR}/lem_lithology_erodibility_work")
lsdtt_set_tests_properties(unit_lem_lithology_erodibility PROPERTIES LABELS "lem")

# k_sn by lithologic unit: contact distance, node -> segment, bookkeeping, and
# the known-answer steady states (two units, n = 1 and n = 2).
add_executable(lsdtt-test-lithology-channel-stats cpp/test_lithology_channel_stats.cpp)
target_link_libraries(lsdtt-test-lithology-channel-stats PRIVATE lsdtt3_cpu)
add_test(NAME unit_chi_lithology_channel_stats
         COMMAND lsdtt-test-lithology-channel-stats
                 "${CMAKE_CURRENT_BINARY_DIR}/chi_lithology_channel_stats_work")
lsdtt_set_tests_properties(unit_chi_lithology_channel_stats PROPERTIES LABELS "chi")

# Erodibility calibration from basin erosion rates: the York regression
# (Pearson's data with York's weights), the erosion-rate sample reader and its
# units, the K formula and property-table seams shared with relative K, the
# declared parameter groups, and basin completeness at an open grid edge.
add_executable(lsdtt-test-erodibility-calibration cpp/test_erodibility_calibration.cpp)
target_link_libraries(lsdtt-test-erodibility-calibration PRIVATE lsdtt3_cpu)
add_test(NAME unit_chi_erodibility_calibration
         COMMAND lsdtt-test-erodibility-calibration
                 "${CMAKE_CURRENT_BINARY_DIR}/chi_erodibility_calibration_work")
lsdtt_set_tests_properties(unit_chi_erodibility_calibration PROPERTIES LABELS "chi")

# lsdtt-chi-analysis with erosion-rate samples. The fixture writer produces the
# F1 (one K, a different uplift per basin) and F2 (two units) steady states, their
# sample CSVs (mm/kyr and m/yr) and the F2 unit polygons. The CLI routes the
# steady states itself (raster_is_filled: true), which moves a few divide pixels,
# so its pins are measured separately from the in-process ones.
set(LSDTT_CHI_ERODIBILITY_CALIBRATION_FIXTURE_DIR "${CMAKE_CURRENT_BINARY_DIR}/chi_erodibility_calibration_cli_fixture")
add_test(NAME fixture_chi_erodibility_calibration_cli
         COMMAND lsdtt-test-erodibility-calibration --write-fixtures
                 "${LSDTT_CHI_ERODIBILITY_CALIBRATION_FIXTURE_DIR}")
lsdtt_set_tests_properties(fixture_chi_erodibility_calibration_cli PROPERTIES
    LABELS "fixture;chi"
    FIXTURES_SETUP chi_erodibility_calibration_cli)

add_lsdtt_cli_smoke_test(
    smoke_chi_analysis_erodibility_calibration
    lsdtt-chi-analysis
    chi_erodibility_calibration.param.in
    ercal_erodibility_calibration_samples_mn_0p450.csv
    --read-path "${LSDTT_CHI_ERODIBILITY_CALIBRATION_FIXTURE_DIR}"
    --expect ercal_erodibility_calibration_basin_K_mn_0p450.csv
    --expect ercal_erodibility_calibration_fit_mn_0p450.csv
    --expect ercal_erodibility_calibration_basins.fgb
    --expect-vector-feature-min ercal_erodibility_calibration_basins.fgb:13
    --expect-file-contains "ercal_erodibility_calibration_samples_mn_0p450.csv:sample_id,input_row,input_latitude,input_longitude,snapped_latitude,snapped_longitude,snapping_method,snap_distance_m,outlet_flat,status,status_detail,reported_area_km2,dem_area_km2,area_ratio,area_check,nesting_group,nested_in,n_samples_inside,duplicate_outlet_of,erosion_rate_input,erosion_rate_uncertainty_input,erosion_rate_input_units,erosion_rate_m_per_yr,erosion_rate_uncertainty_m_per_yr,m_over_n,A_0,k_sn_statistic,k_sn,k_sn_ci_low,k_sn_ci_high,n_reaches,n_reaches_qualifying,k_sn_reach_iqr_over_median,k_sn_chi_z_slope,chi_z_r_squared,k_sn_segment_median,k_sn_chi_weighted_mean,dominant_lith_code,dominant_lith_name,dominant_area_fraction,dominant_channel_length_fraction,unmapped_area_fraction"
    --expect-file-contains "ercal_erodibility_calibration_fit_mn_0p450.csv:fit_type,method,status,N_basins,N_groups,ksn_range_ratio,m_over_n,A_0,n,n_ci_low,n_ci_high,n_se_york,n_se_york_mswd_scaled,m,K,K_ci_low,K_ci_high,ln_K_se,ln_K_se_mswd_scaled,corr_n_lnK,E_at_median_ksn_m_per_yr,mswd,ln_K_dispersion,corr_lnK_lnk,n_bootstrap,bootstrap_cluster,K_units"
    --expect-file-contains "ercal_erodibility_calibration_basin_K_mn_0p450.csv:sample_id,status,nesting_group,m_over_n,A_0,n,m,erosion_rate_m_per_yr,erosion_rate_uncertainty_m_per_yr,k_sn,k_sn_ci_low,k_sn_ci_high,K,K_ci_low,K_ci_high,n_bootstrap,K_units"
    --expect-file-contains "ercal_erodibility_calibration_samples_mn_0p450.csv:,area_mismatch,"
    --expect-file-contains "ercal_erodibility_calibration_samples_mn_0p450.csv:,outside_dem,"
    --expect-file-contains "ercal_erodibility_calibration_fit_mn_0p450.csv:fixed_n,weighted_log_mean,ok,12,10,"
    --expect-file-contains "ercal_erodibility_calibration_basin_K_mn_0p450.csv:m^(-0.8) yr^-1"
    --stdout-contains "Erodibility calibration (m/n = 0.45, A0 = 1): 12 basin(s) used in 10 independent group(s)"
    --stdout-contains "area mismatch (outside 1/1.25..1.25 of the reported area): 1 (exclude)"
    --stdout-contains "fitted: n = "
)
lsdtt_set_tests_properties(smoke_chi_analysis_erodibility_calibration PROPERTIES
    FIXTURES_REQUIRED chi_erodibility_calibration_cli
    FIXTURES_SETUP smoke_chi_analysis_erodibility_calibration)

add_lsdtt_cli_smoke_test(
    smoke_chi_analysis_erodibility_calibration_m_per_yr
    lsdtt-chi-analysis
    chi_erodibility_calibration_m_per_yr.param.in
    ercalm_erodibility_calibration_fit_mn_0p450.csv
    --read-path "${LSDTT_CHI_ERODIBILITY_CALIBRATION_FIXTURE_DIR}"
    --stdout-contains "(units m_per_yr -> m/yr x 1)"
)
lsdtt_set_tests_properties(smoke_chi_analysis_erodibility_calibration_m_per_yr PROPERTIES
    FIXTURES_REQUIRED chi_erodibility_calibration_cli
    FIXTURES_SETUP smoke_chi_analysis_erodibility_calibration_m_per_yr)

add_lsdtt_cli_smoke_test(
    smoke_chi_analysis_erodibility_calibration_off
    lsdtt-chi-analysis
    chi_erodibility_calibration_off.param.in
    ercaloff_chi_fit_nodes_mn_0p450.csv
    --read-path "${LSDTT_CHI_ERODIBILITY_CALIBRATION_FIXTURE_DIR}"
    --expect-absent ercaloff_erodibility_calibration_samples_mn_0p450.csv
)
lsdtt_set_tests_properties(smoke_chi_analysis_erodibility_calibration_off PROPERTIES
    FIXTURES_REQUIRED chi_erodibility_calibration_cli
    FIXTURES_SETUP smoke_chi_analysis_erodibility_calibration_off)

add_lsdtt_cli_smoke_test(
    smoke_chi_analysis_erodibility_calibration_lithology
    lsdtt-chi-analysis
    chi_erodibility_calibration_lithology.param.in
    ercallith_erodibility_calibration_by_lithology_mn_0p450.csv
    --read-path "${LSDTT_CHI_ERODIBILITY_CALIBRATION_FIXTURE_DIR}"
    --expect ercallith_erodibility_calibration_property_table_mn_0p450_n_1.txt
    --expect ercallith_erodibility_calibration_basin_lithology.csv
    --expect-file-contains "ercallith_erodibility_calibration_by_lithology_mn_0p450.csv:lith_code,lith_name,method,m_over_n,A_0,n,m,N_basins,N_groups,K,K_ci_low,K_ci_high,n_bootstrap,low_sample,status,purity_threshold,K_units"
    --expect-file-contains "ercallith_erodibility_calibration_by_lithology_mn_0p450.csv:,west,mono,"
    --expect-file-contains "ercallith_erodibility_calibration_by_lithology_mn_0p450.csv:,east,mono,"
    --expect-file-contains "ercallith_erodibility_calibration_by_lithology_mn_0p450.csv:,west,unit_in_basin,"
    --expect-file-contains "ercallith_erodibility_calibration_by_lithology_mn_0p450.csv:,east,unit_in_basin,"
    --expect-file-contains "ercallith_erodibility_calibration_basin_lithology.csv:sample_id,lith_code,lith_name,area_km2,area_fraction,channel_length_m,channel_length_fraction"
    --expect-file-contains "ercallith_erodibility_calibration_property_table_mn_0p450_n_1.txt:# erodibility_convention: map_area"
    --expect-file-contains "ercallith_erodibility_calibration_property_table_mn_0p450_n_1.txt:# m_exponent: 0.45"
    --expect-file-contains "ercallith_erodibility_calibration_property_table_mn_0p450_n_1.txt:code,name,erodibility,reference"
    --stdout-contains "by lithology (unit_in_basin"
)
lsdtt_set_tests_properties(smoke_chi_analysis_erodibility_calibration_lithology PROPERTIES
    FIXTURES_REQUIRED chi_erodibility_calibration_cli
    FIXTURES_SETUP smoke_chi_analysis_erodibility_calibration_lithology)

if(TEST smoke_chi_analysis_erodibility_calibration)
    set(_ercal_out "${CMAKE_CURRENT_BINARY_DIR}/smoke_chi_analysis_erodibility_calibration/out")
    # The CLI known answer (F1, n = 1, K = 1e-5). Measured 2026-10-19: fitted
    # n 1.04242 (abs err 0.0424; K at that n 8.342e-6), fixed-n K 9.53205e-6
    # (rel err 4.68 %). Twice the measured error exceeds the plan's CLI caps,
    # so the pins are the caps: 0.05 on n, 6.4 % on K.
    add_test(NAME rigorous_chi_erodibility_calibration_n
             COMMAND lsdtt-test-erodibility-calibration --check-erodibility-calibration-fit
                     "${_ercal_out}/ercal_erodibility_calibration_fit_mn_0p450.csv" fitted_n 0 1 0.05 1e-5 1e9)
    add_test(NAME rigorous_chi_erodibility_calibration_fixed_n_K
             COMMAND lsdtt-test-erodibility-calibration --check-erodibility-calibration-fit
                     "${_ercal_out}/ercal_erodibility_calibration_fit_mn_0p450.csv" fixed_n 1 1 0 1e-5 0.064)
    foreach(_t rigorous_chi_erodibility_calibration_n rigorous_chi_erodibility_calibration_fixed_n_K)
        lsdtt_set_tests_properties(${_t} PROPERTIES
            LABELS "chi"
            FIXTURES_REQUIRED smoke_chi_analysis_erodibility_calibration)
    endforeach()
    # The same samples in m/yr give the same basin K and fit files, byte for byte.
    add_test(NAME rigorous_chi_erodibility_calibration_units_invariant
             COMMAND ${CMAKE_COMMAND}
                 "-DPREFIX_A=${_ercal_out}/ercal"
                 "-DPREFIX_B=${CMAKE_CURRENT_BINARY_DIR}/smoke_chi_analysis_erodibility_calibration_m_per_yr/out/ercalm"
                 "-DSUFFIXES=_erodibility_calibration_basin_K_mn_0p450.csv;_erodibility_calibration_fit_mn_0p450.csv"
                 -P "${CMAKE_CURRENT_SOURCE_DIR}/cmake/compare_output_sets.cmake")
    lsdtt_set_tests_properties(rigorous_chi_erodibility_calibration_units_invariant PROPERTIES
        LABELS "cli;chi"
        FIXTURES_REQUIRED "smoke_chi_analysis_erodibility_calibration;smoke_chi_analysis_erodibility_calibration_m_per_yr")
    # The calibration only observes: the chi fit tables, k_sn basin summary and
    # disorder sweep are the same with and without erosion_rate_fname.
    add_test(NAME rigorous_chi_erodibility_calibration_observation_only
             COMMAND ${CMAKE_COMMAND}
                 "-DPREFIX_A=${_ercal_out}/ercal"
                 "-DPREFIX_B=${CMAKE_CURRENT_BINARY_DIR}/smoke_chi_analysis_erodibility_calibration_off/out/ercaloff"
                 "-DSUFFIXES=_chi_fit_nodes_mn_0p450.csv;_chi_fit_segments_mn_0p450.csv;_chi_ksn_basin_summary_mn_0p450.csv;_chi_disorder_sweep.csv"
                 -P "${CMAKE_CURRENT_SOURCE_DIR}/cmake/compare_output_sets.cmake")
    lsdtt_set_tests_properties(rigorous_chi_erodibility_calibration_observation_only PROPERTIES
        LABELS "cli;chi"
        FIXTURES_REQUIRED "smoke_chi_analysis_erodibility_calibration;smoke_chi_analysis_erodibility_calibration_off")
    # The CLI known answer per unit (F2, unit_in_basin). Measured 2026-10-19:
    # west 9.46848e-6 (rel err 5.32 %), east 3.86206e-5 (rel err 3.45 %); twice
    # those exceed the plan's 6.4 % cap, so the pin is the cap. (mono, measured
    # west 7.05 %, east 3.58 %, is not asserted: west exceeds the cap.)
    foreach(_u "west;1e-5" "east;4e-5")
        list(GET _u 0 _name)
        list(GET _u 1 _K)
        add_test(NAME rigorous_chi_erodibility_calibration_lithology_units_${_name}
                 COMMAND lsdtt-test-erodibility-calibration --check-erodibility-calibration-unit
                         "${CMAKE_CURRENT_BINARY_DIR}/smoke_chi_analysis_erodibility_calibration_lithology/out/ercallith_erodibility_calibration_by_lithology_mn_0p450.csv"
                         unit_in_basin 1 ${_name} ${_K} 0.064)
        lsdtt_set_tests_properties(rigorous_chi_erodibility_calibration_lithology_units_${_name} PROPERTIES
            LABELS "chi"
            FIXTURES_REQUIRED smoke_chi_analysis_erodibility_calibration_lithology)
    endforeach()
else()
    foreach(_t n fixed_n_K units_invariant observation_only lithology_units_west lithology_units_east)
        lsdtt_mark_test_unregistered(python rigorous_chi_erodibility_calibration_${_t})
    endforeach()
endif()

foreach(_f "no_units;erosion_rate_units must be one of"
           "units_conflict;disagrees with erosion_rate_units"
           "without_fit;erosion_rate_fname needs fit_chi_networks: true"
           "missing_column;erosion_rate_column 'rate' is not a column"
           "missing_file;erosion_rate_fname does not exist")
    list(GET _f 0 _name)
    list(GET _f 1 _msg)
    add_lsdtt_cli_smoke_test(
        failure_chi_analysis_erodibility_calibration_${_name}
        lsdtt-chi-analysis
        chi_erodibility_calibration_${_name}.param.in
        unused_expected_output
        --read-path "${LSDTT_CHI_ERODIBILITY_CALIBRATION_FIXTURE_DIR}"
        --expect-failure
        --stderr-contains "${_msg}")
    lsdtt_set_tests_properties(failure_chi_analysis_erodibility_calibration_${_name} PROPERTIES
        FIXTURES_REQUIRED chi_erodibility_calibration_cli)
endforeach()

# lsdtt-chi-analysis with a lithology map. The fixture writer produces the
# steady-state landscape (two units, k_sn about 20 and 5) with its unit polygons,
# and a two-unit code raster beside a copy of tiny_dem.asc.
set(LSDTT_CHI_LITHOLOGY_FIXTURE_DIR "${CMAKE_CURRENT_BINARY_DIR}/chi_lithology_cli_fixture")
add_test(NAME fixture_chi_lithology_cli
         COMMAND lsdtt-test-lithology-channel-stats --write-fixtures
                 "${LSDTT_CHI_LITHOLOGY_FIXTURE_DIR}" "${LSDTT_TEST_DATA_DIR}/tiny_dem.asc")
lsdtt_set_tests_properties(fixture_chi_lithology_cli PROPERTIES
    LABELS "fixture;chi"
    FIXTURES_SETUP chi_lithology_cli)

add_lsdtt_cli_smoke_test(
    smoke_chi_lithology_ksn
    lsdtt-chi-analysis
    chi_lithology_ksn.param.in
    lithks_chi_ksn_by_lithology_mn_0p450.csv
    --read-path "${LSDTT_CHI_LITHOLOGY_FIXTURE_DIR}"
    --expect lithks_chi_ksn_histogram_by_lithology_mn_0p450.csv
    --expect lithks_chi_ksn_lithology_nodes_mn_0p450.csv
    --expect lithks_chi_lithology_reaches_mn_0p450.csv
    --expect lithks_chi_lithology_segments_mn_0p450.csv
    --expect lithks_chi_ksn_by_lithology_and_basin_mn_0p450.csv
    --expect lithks_lithology_units.csv
    --expect lithks_lithology_code.tif
    --expect-raster-values lithks_lithology_code.tif:1,2
    --expect-file-contains "lithks_chi_ksn_by_lithology_mn_0p450.csv:2,shale,"
    --expect-file-contains "lithks_chi_ksn_by_lithology_mn_0p450.csv:-1,unmapped,0,"
    --expect-file-contains "lithks_chi_ksn_by_lithology_mn_0p450.csv:reach_k_sn_median,reach_k_sn_p25,reach_k_sn_p75,reach_k_sn_length_weighted_median,reach_k_sn_median_ci_low,reach_k_sn_median_ci_high,n_bootstrap,bootstrap_cluster"
    # Run-local codes 1..N from the map's own text values (no crosswalk).
    --expect-file-contains "lithks_lithology_units.csv:1,sandstone,sandstone,9600,8.64,"
    --expect-file-contains "lithks_lithology_units.csv:2,shale,shale,9600,8.64,"
    --expect-file-contains "lithks_lithology_units.csv:-1,unmapped,,0,0,0"
    --expect-file-contains "lithks_chi_lithology_reaches_mn_0p450.csv:upstream_end,downstream_end,qualifies"
    # The units separate: shale peaks in bin 0 (k_sn 4.6-5.7) where sandstone has
    # no node, sandstone in bin 14 (19.8-20.9) where shale has none.
    --expect-file-contains "lithks_chi_ksn_histogram_by_lithology_mn_0p450.csv:1,sandstone,0,4.632725239,5.716628242,0,"
    --expect-file-contains "lithks_chi_ksn_histogram_by_lithology_mn_0p450.csv:2,shale,14,19.80736728,20.89127028,0,"
    --stdout-contains "2 lithologic unit(s)"
    --stdout-contains "Reading k_sn by lithologic unit:"
    --stdout-contains "reach median k_sn and its reach-bootstrap interval"
)
lsdtt_set_tests_properties(smoke_chi_lithology_ksn PROPERTIES
    FIXTURES_REQUIRED chi_lithology_cli
    FIXTURES_SETUP smoke_chi_lithology_ksn)

# Lithology is observation only: the fit nodes CSV keeps the hash pinned by
# rigorous_chi_nodes_sha256 (chi_analysis_smoke.param.in without the map).
add_lsdtt_cli_smoke_test(
    rigorous_chi_lithology_observation_only
    lsdtt-chi-analysis
    chi_lithology_observation_only.param.in
    tiny_chi_fit_nodes_mn_0p500.csv
    --read-path "${LSDTT_CHI_LITHOLOGY_FIXTURE_DIR}"
    --expect tiny_chi_ksn_by_lithology_mn_0p500.csv
    --expect-sha256 tiny_chi_fit_nodes_mn_0p500.csv:b64aa1266abc09369b423534d44ec5ed8a3ceac52d7d31e9f8637c8893c26ad3
)
lsdtt_set_tests_properties(rigorous_chi_lithology_observation_only PROPERTIES
    FIXTURES_REQUIRED chi_lithology_cli)

add_lsdtt_cli_smoke_test(
    failure_chi_lithology_without_fit
    lsdtt-chi-analysis
    chi_lithology_without_fit.param.in
    unused_expected_output
    --read-path "${LSDTT_CHI_LITHOLOGY_FIXTURE_DIR}"
    --expect-failure
    --stderr-contains "FATAL: lithology_fname needs fit_chi_networks: true"
)
lsdtt_set_tests_properties(failure_chi_lithology_without_fit PROPERTIES
    FIXTURES_REQUIRED chi_lithology_cli)

# The CLI known answer. The CLI re-routes the steady-state DEM, so
# its routing differs from the generating tree near channel heads. Measured
# 2026-10-17: reach medians sandstone 20.40326316, shale 5.163420929, ratio
# 3.951501038, relative error 1.212 % against 4; pinned at twice that, 2.42 %.
# Native, but its input comes from a Python-harness test, so it is registered
# with the Python tests.
if(LSDTT_ENABLE_PYTHON_TESTS)
    add_test(NAME rigorous_chi_lithology_ksn_ratio
             COMMAND lsdtt-test-lithology-channel-stats --check-ksn-by-lithology
                     "${CMAKE_CURRENT_BINARY_DIR}/smoke_chi_lithology_ksn/out/lithks_chi_ksn_by_lithology_mn_0p450.csv"
                     sandstone shale 4 0.0242)
    lsdtt_set_tests_properties(rigorous_chi_lithology_ksn_ratio PROPERTIES
        LABELS "chi"
        FIXTURES_REQUIRED smoke_chi_lithology_ksn)
else()
    lsdtt_mark_test_unregistered(python rigorous_chi_lithology_ksn_ratio)
endif()

# A crosswalk with a name column, a contact buffer and basin
# clusters. Codes come from the crosswalk, names from its `name` column.
add_lsdtt_cli_smoke_test(
    smoke_chi_lithology_crosswalk
    lsdtt-chi-analysis
    chi_lithology_crosswalk.param.in
    lithkscw_chi_ksn_by_lithology_mn_0p450.csv
    --read-path "${LSDTT_CHI_LITHOLOGY_FIXTURE_DIR}"
    --expect lithkscw_lithology_units.csv
    --expect lithkscw_chi_lithology_reaches_mn_0p450.csv
    --expect-file-contains "lithkscw_lithology_units.csv:10,quartz sandstone,sandstone,9600,"
    --expect-file-contains "lithkscw_lithology_units.csv:30,grey shale,shale,9600,"
    --expect-file-contains "lithkscw_chi_ksn_by_lithology_mn_0p450.csv:10,quartz sandstone,"
    # 60 m buffer: 25 + 27 = 52 nodes set aside; 200 replicates by basin.
    --expect-file-contains "lithkscw_chi_ksn_by_lithology_mn_0p450.csv:-2,all,1339,"
    --expect-file-contains "lithkscw_chi_ksn_by_lithology_mn_0p450.csv:,60,10,52,84,26,"
    --expect-file-contains "lithkscw_chi_ksn_by_lithology_mn_0p450.csv:,200,basin,"
    --expect-file-contains "lithkscw_chi_lithology_reaches_mn_0p450.csv:,buffer,"
    --stdout-contains "basin-bootstrap interval"
)
lsdtt_set_tests_properties(smoke_chi_lithology_crosswalk PROPERTIES FIXTURES_REQUIRED chi_lithology_cli)

# Refusals: a missing map, an unknown field, a value missing from the crosswalk,
# a bootstrap cluster other than reach or basin.
add_lsdtt_cli_smoke_test(
    failure_chi_lithology_missing_file
    lsdtt-chi-analysis
    chi_lithology_missing_file.param.in
    unused_expected_output
    --read-path "${LSDTT_CHI_LITHOLOGY_FIXTURE_DIR}"
    --expect-failure
    --stderr-contains "lithology map does not exist"
)
lsdtt_set_tests_properties(failure_chi_lithology_missing_file PROPERTIES FIXTURES_REQUIRED chi_lithology_cli)

add_lsdtt_cli_smoke_test(
    failure_chi_lithology_bad_field
    lsdtt-chi-analysis
    chi_lithology_bad_field.param.in
    unused_expected_output
    --read-path "${LSDTT_CHI_LITHOLOGY_FIXTURE_DIR}"
    --expect-failure
    --stderr-contains "has no field 'formation'"
    --stderr-contains "its fields are: lith_primary"
)
lsdtt_set_tests_properties(failure_chi_lithology_bad_field PROPERTIES FIXTURES_REQUIRED chi_lithology_cli)

add_lsdtt_cli_smoke_test(
    failure_chi_lithology_crosswalk_gap
    lsdtt-chi-analysis
    chi_lithology_crosswalk_gap.param.in
    unused_expected_output
    --read-path "${LSDTT_CHI_LITHOLOGY_FIXTURE_DIR}"
    --expect-failure
    --stderr-contains "'shale' (9600 pixels)"
)
lsdtt_set_tests_properties(failure_chi_lithology_crosswalk_gap PROPERTIES FIXTURES_REQUIRED chi_lithology_cli)

add_lsdtt_cli_smoke_test(
    failure_chi_lithology_bad_cluster
    lsdtt-chi-analysis
    chi_lithology_bad_cluster.param.in
    unused_expected_output
    --read-path "${LSDTT_CHI_LITHOLOGY_FIXTURE_DIR}"
    --expect-failure
    --stderr-contains "lithology_bootstrap_cluster must be reach or basin, not 'segment'"
)
lsdtt_set_tests_properties(failure_chi_lithology_bad_cluster PROPERTIES FIXTURES_REQUIRED chi_lithology_cli)

# A lithology map that labels no valid DEM pixel (200 km east, same CRS) is
# refused with both extents and CRSs: a code raster and a vector.
foreach(_lithks_far raster vector)
    add_lsdtt_cli_smoke_test(
        failure_chi_lithology_no_overlap_${_lithks_far}
        lsdtt-chi-analysis
        chi_lithology_no_overlap_${_lithks_far}.param.in
        unused_expected_output
        --read-path "${LSDTT_CHI_LITHOLOGY_FIXTURE_DIR}"
        --expect-failure
        --stderr-contains "labels no valid DEM pixel: the map covers x 700000.0 to 704800.0, y 5000000.0 to 5003600.0 (WGS 84 / UTM zone 32N)"
        --stderr-contains "the DEM covers x 500000.0 to 504800.0, y 5000000.0 to 5003600.0 (WGS 84 / UTM zone 32N)"
    )
    lsdtt_set_tests_properties(failure_chi_lithology_no_overlap_${_lithks_far} PROPERTIES
        FIXTURES_REQUIRED chi_lithology_cli)
endforeach()

# The bootstrap stream is per m/n: the 0.45 tables of a run with m/n 0.4 and
# 0.45 equal those of the single-m/n run byte for byte.
add_lsdtt_cli_smoke_test(
    smoke_chi_lithology_multi_mn
    lsdtt-chi-analysis
    chi_lithology_multi_mn.param.in
    lithksmn_chi_ksn_by_lithology_mn_0p450.csv
    --read-path "${LSDTT_CHI_LITHOLOGY_FIXTURE_DIR}"
    --expect lithksmn_chi_ksn_by_lithology_mn_0p400.csv
)
if(TEST smoke_chi_lithology_multi_mn)
    lsdtt_set_tests_properties(smoke_chi_lithology_multi_mn PROPERTIES
        FIXTURES_REQUIRED chi_lithology_cli
        FIXTURES_SETUP smoke_chi_lithology_multi_mn)
    add_test(NAME rigorous_chi_lithology_bootstrap_per_mn
             COMMAND ${CMAKE_COMMAND}
                 "-DPREFIX_A=${CMAKE_CURRENT_BINARY_DIR}/smoke_chi_lithology_ksn/out/lithks"
                 "-DPREFIX_B=${CMAKE_CURRENT_BINARY_DIR}/smoke_chi_lithology_multi_mn/out/lithksmn"
                 "-DSUFFIXES=_chi_ksn_by_lithology_mn_0p450.csv;_chi_ksn_by_lithology_and_basin_mn_0p450.csv;_chi_lithology_reaches_mn_0p450.csv"
                 -P "${CMAKE_CURRENT_SOURCE_DIR}/cmake/compare_output_sets.cmake")
    lsdtt_set_tests_properties(rigorous_chi_lithology_bootstrap_per_mn PROPERTIES
        LABELS "cli;chi"
        FIXTURES_REQUIRED "smoke_chi_lithology_ksn;smoke_chi_lithology_multi_mn")
else()
    lsdtt_mark_test_unregistered(python rigorous_chi_lithology_bootstrap_per_mn)
endif()

# Relative erodibility from k_sn by lithologic unit (write_relative_erodibility):
# K_i / K_ref = (k_ref / k_i)^n and, with E, K_i = E / k_i^n, on the n = 1 steady
# state (K1 = 1e-5 sandstone, K2 = 4e-5 shale, U = 2e-4) with E = U. The n = 2
# rows of this run apply n = 2 to an n = 1 landscape; the n = 2 known answer is
# in the driver (R4).
add_lsdtt_cli_smoke_test(
    smoke_chi_lithology_relative_erodibility
    lsdtt-chi-analysis
    chi_lithology_relative_k.param.in
    lithkrk_chi_relative_erodibility_by_lithology_mn_0p450.csv
    --read-path "${LSDTT_CHI_LITHOLOGY_FIXTURE_DIR}"
    --expect lithkrk_chi_relative_erodibility_property_table_mn_0p450_n_1.txt
    --expect lithkrk_chi_relative_erodibility_property_table_mn_0p450_n_2.txt
    --expect lithkrk_chi_ksn_by_lithology_mn_0p450.csv
    --expect-file-contains "lithkrk_chi_relative_erodibility_by_lithology_mn_0p450.csv:lith_code,lith_name,m_over_n,A_0,n,m,reference_unit,n_reaches_qualifying,reach_k_sn_median,k_ratio,k_ratio_ci_low,k_ratio_ci_high,n_bootstrap,bootstrap_cluster,low_sample,erosion_rate_m_per_yr,K_abs,K_abs_ci_low,K_abs_ci_high,K_units"
    --expect-file-contains "lithkrk_chi_relative_erodibility_by_lithology_mn_0p450.csv:1,sandstone,0.45,1,1,0.45,sandstone,15,"
    --expect-file-contains "lithkrk_chi_relative_erodibility_by_lithology_mn_0p450.csv:2,shale,0.45,1,2,0.9,sandstone,11,"
    --expect-file-contains "lithkrk_chi_relative_erodibility_by_lithology_mn_0p450.csv:,1000,reach,0,0.0002,"
    --expect-file-contains "lithkrk_chi_relative_erodibility_by_lithology_mn_0p450.csv:m^(0.1) yr^-1"
    --expect-file-contains "lithkrk_chi_relative_erodibility_by_lithology_mn_0p450.csv:m^(-0.8) yr^-1"
    --expect-file-not-contains "lithkrk_chi_relative_erodibility_by_lithology_mn_0p450.csv:unmapped"
    --expect-file-contains "lithkrk_chi_relative_erodibility_property_table_mn_0p450_n_1.txt:# erodibility_convention: map_area"
    --expect-file-contains "lithkrk_chi_relative_erodibility_property_table_mn_0p450_n_1.txt:# m_exponent: 0.45"
    --expect-file-contains "lithkrk_chi_relative_erodibility_property_table_mn_0p450_n_1.txt:# n_exponent: 1"
    --expect-file-contains "lithkrk_chi_relative_erodibility_property_table_mn_0p450_n_1.txt:DEM lithks_steady_dem.tif"
    --expect-file-contains "lithkrk_chi_relative_erodibility_property_table_mn_0p450_n_1.txt:lithology map lithks_units.gpkg (field lith_primary)"
    --expect-file-contains "lithkrk_chi_relative_erodibility_property_table_mn_0p450_n_1.txt:assumes uniform erosion and steady state"
    --expect-file-contains "lithkrk_chi_relative_erodibility_property_table_mn_0p450_n_1.txt:code,name,erodibility,reference"
    --expect-file-contains "lithkrk_chi_relative_erodibility_property_table_mn_0p450_n_2.txt:# m_exponent: 0.9"
    --stdout-contains "reference 'sandstone' (default: the most qualifying reaches (15))"
    --stdout-contains "These values hold only under these conditions:"
)
lsdtt_set_tests_properties(smoke_chi_lithology_relative_erodibility PROPERTIES
    FIXTURES_REQUIRED chi_lithology_cli
    FIXTURES_SETUP smoke_chi_lithology_relative_erodibility)

if(TEST smoke_chi_lithology_relative_erodibility)
    # Switching relative erodibility on leaves every k_sn-by-lithology output
    # of the same run byte-identical (its bootstrap stream is separate).
    add_test(NAME rigorous_chi_lithology_relative_erodibility_observation_only
             COMMAND ${CMAKE_COMMAND}
                 "-DPREFIX_A=${CMAKE_CURRENT_BINARY_DIR}/smoke_chi_lithology_ksn/out/lithks"
                 "-DPREFIX_B=${CMAKE_CURRENT_BINARY_DIR}/smoke_chi_lithology_relative_erodibility/out/lithkrk"
                 "-DSUFFIXES=_chi_ksn_by_lithology_mn_0p450.csv;_chi_ksn_by_lithology_and_basin_mn_0p450.csv;_chi_ksn_histogram_by_lithology_mn_0p450.csv;_chi_ksn_lithology_nodes_mn_0p450.csv;_chi_lithology_reaches_mn_0p450.csv;_chi_lithology_segments_mn_0p450.csv;_lithology_units.csv"
                 -P "${CMAKE_CURRENT_SOURCE_DIR}/cmake/compare_output_sets.cmake")
    lsdtt_set_tests_properties(rigorous_chi_lithology_relative_erodibility_observation_only PROPERTIES
        LABELS "cli;chi"
        FIXTURES_REQUIRED "smoke_chi_lithology_ksn;smoke_chi_lithology_relative_erodibility")

    # The CLI known answer (n = 1, E = U). The CLI re-routes the steady-state
    # DEM. Measured 2026-10-19: shale K/K_ref 3.951501037 (rel err 1.212 %, the
    # same reach medians as rigorous_chi_lithology_ksn_ratio), K shale
    # 3.87340104e-05 (rel err 3.165 %), K sandstone 9.802353596e-06 (rel err
    # 1.976 %). Pins at twice the measured error: ratio 2.42 %, K shale 6.4 %,
    # K sandstone 4 %.
    add_test(NAME rigorous_chi_lithology_relative_erodibility_ratio
             COMMAND lsdtt-test-lithology-channel-stats --check-relative-erodibility
                     "${CMAKE_CURRENT_BINARY_DIR}/smoke_chi_lithology_relative_erodibility/out/lithkrk_chi_relative_erodibility_by_lithology_mn_0p450.csv"
                     shale 1 4 0.0242 4e-5 0.064)
    lsdtt_set_tests_properties(rigorous_chi_lithology_relative_erodibility_ratio PROPERTIES
        LABELS "chi"
        FIXTURES_REQUIRED smoke_chi_lithology_relative_erodibility)
    add_test(NAME rigorous_chi_lithology_relative_erodibility_reference
             COMMAND lsdtt-test-lithology-channel-stats --check-relative-erodibility
                     "${CMAKE_CURRENT_BINARY_DIR}/smoke_chi_lithology_relative_erodibility/out/lithkrk_chi_relative_erodibility_by_lithology_mn_0p450.csv"
                     sandstone 1 1 1e-12 1e-5 0.04)
    lsdtt_set_tests_properties(rigorous_chi_lithology_relative_erodibility_reference PROPERTIES
        LABELS "chi"
        FIXTURES_REQUIRED smoke_chi_lithology_relative_erodibility)
else()
    lsdtt_mark_test_unregistered(python rigorous_chi_lithology_relative_erodibility_observation_only)
    lsdtt_mark_test_unregistered(python rigorous_chi_lithology_relative_erodibility_ratio)
    lsdtt_mark_test_unregistered(python rigorous_chi_lithology_relative_erodibility_reference)
endif()

# The loop closes: lsdtt-lithology-erodibility reads the n = 1 draft table
# unchanged and, with the same DEM and map, writes the K raster (the two K
# values of the table, min sandstone and max shale).
add_lsdtt_cli_smoke_test(
    smoke_chi_lithology_relative_erodibility_k_raster
    lsdtt-lithology-erodibility
    chi_lithology_relative_k_raster.param.in
    lithkrkk_erodibility.tif
    --read-path "${LSDTT_CHI_LITHOLOGY_FIXTURE_DIR}"
    --expect lithkrkk_lithology_legend.csv
    --expect-raster-size lithkrkk_erodibility.tif:160x120
    --expect-raster-stat lithkrkk_erodibility.tif:min:9.80235e-06:9.80236e-06
    --expect-raster-stat lithkrkk_erodibility.tif:max:3.87340e-05:3.87341e-05
    --expect-file-contains "lithkrkk_lem_erodibility_parameters.txt:erodibility_raster_m_exponent: 0.45"
    --expect-file-contains "lithkrkk_lem_erodibility_parameters.txt:erodibility_raster_area_convention: map_area"
)
lsdtt_set_tests_properties(smoke_chi_lithology_relative_erodibility_k_raster PROPERTIES
    FIXTURES_REQUIRED "chi_lithology_cli;smoke_chi_lithology_relative_erodibility")

# Per m/n stream: with m/n 0.4 and 0.45 the 0.45 relative-erodibility outputs
# equal those of the single-m/n run byte for byte.
add_lsdtt_cli_smoke_test(
    smoke_chi_lithology_relative_erodibility_multi_mn
    lsdtt-chi-analysis
    chi_lithology_relative_k_multi_mn.param.in
    lithkrkmn_chi_relative_erodibility_by_lithology_mn_0p450.csv
    --read-path "${LSDTT_CHI_LITHOLOGY_FIXTURE_DIR}"
    --expect lithkrkmn_chi_relative_erodibility_by_lithology_mn_0p400.csv
)
if(TEST smoke_chi_lithology_relative_erodibility_multi_mn)
    lsdtt_set_tests_properties(smoke_chi_lithology_relative_erodibility_multi_mn PROPERTIES
        FIXTURES_REQUIRED chi_lithology_cli
        FIXTURES_SETUP smoke_chi_lithology_relative_erodibility_multi_mn)
    add_test(NAME rigorous_chi_lithology_relative_erodibility_per_mn
             COMMAND ${CMAKE_COMMAND}
                 "-DPREFIX_A=${CMAKE_CURRENT_BINARY_DIR}/smoke_chi_lithology_relative_erodibility/out/lithkrk"
                 "-DPREFIX_B=${CMAKE_CURRENT_BINARY_DIR}/smoke_chi_lithology_relative_erodibility_multi_mn/out/lithkrkmn"
                 "-DSUFFIXES=_chi_relative_erodibility_by_lithology_mn_0p450.csv;_chi_relative_erodibility_property_table_mn_0p450_n_1.txt;_chi_relative_erodibility_property_table_mn_0p450_n_2.txt"
                 -P "${CMAKE_CURRENT_SOURCE_DIR}/cmake/compare_output_sets.cmake")
    lsdtt_set_tests_properties(rigorous_chi_lithology_relative_erodibility_per_mn PROPERTIES
        LABELS "cli;chi"
        FIXTURES_REQUIRED "smoke_chi_lithology_relative_erodibility;smoke_chi_lithology_relative_erodibility_multi_mn")
else()
    lsdtt_mark_test_unregistered(python rigorous_chi_lithology_relative_erodibility_per_mn)
endif()

# Refusals: relative erodibility without a lithology map; a reference that is
# not a unit of the map (the units are listed).
add_lsdtt_cli_smoke_test(
    failure_chi_lithology_relative_erodibility_no_map
    lsdtt-chi-analysis
    chi_lithology_relative_k_no_map.param.in
    unused_expected_output
    --read-path "${LSDTT_CHI_LITHOLOGY_FIXTURE_DIR}"
    --expect-failure
    --stderr-contains "FATAL: write_relative_erodibility needs lithology_fname"
)
lsdtt_set_tests_properties(failure_chi_lithology_relative_erodibility_no_map PROPERTIES
    FIXTURES_REQUIRED chi_lithology_cli)

add_lsdtt_cli_smoke_test(
    failure_chi_lithology_relative_erodibility_bad_reference
    lsdtt-chi-analysis
    chi_lithology_relative_k_bad_reference.param.in
    unused_expected_output
    --read-path "${LSDTT_CHI_LITHOLOGY_FIXTURE_DIR}"
    --expect-failure
    --stderr-contains "relative_erodibility_reference_unit 'granite' is not a unit of the lithology map"
    --stderr-contains "the units are: sandstone, shale"
)
lsdtt_set_tests_properties(failure_chi_lithology_relative_erodibility_bad_reference PROPERTIES
    FIXTURES_REQUIRED chi_lithology_cli)

# Channel-link bearings and junction angles by lithologic unit (T12, T13): link
# and junction tags on the synthetic trellis with in-process unit grids, the
# axial rose by unit, the junction median bootstrap, the appended record
# columns and the analysed-area mask.
add_executable(lsdtt-test-lithology-network-stats cpp/test_lithology_network_stats.cpp)
target_link_libraries(lsdtt-test-lithology-network-stats PRIVATE lsdtt3_cpu)
add_test(NAME unit_lithology_channel_network_stats
         COMMAND lsdtt-test-lithology-network-stats
                 "${CMAKE_CURRENT_BINARY_DIR}/lithology_channel_network_stats_work")
lsdtt_set_tests_properties(unit_lithology_channel_network_stats PROPERTIES LABELS "network")

# lsdtt-channel-network-analysis with a lithology code raster (the trellis
# fixture's north/south halves). These use lsdtt_add_channel_bearing_cli_test,
# defined with the bearing-statistics tests above.
lsdtt_add_channel_bearing_cli_test(
    smoke_channel_network_analysis_lithology
    channel_network_analysis_lithology.param.in
    cbslith_channel_bearing_group_summary.csv
    --expect cbslith_lithology_units.csv
    --expect cbslith_lithology_code.tif
    --expect cbslith_junction_angle_group_summary.csv
    --expect-file-contains "cbslith_channel_link_bearings.csv:downstream_longitude,lith_code,lith_name,lith_majority_fraction"
    --expect-file-contains "cbslith_junction_angles.csv:longitude,junction_lith_code,junction_lith_name,donor_1_lith_code,donor_1_lith_fraction,donor_2_lith_code,donor_2_lith_fraction,receiver_lith_code,receiver_lith_fraction,lith_all_arms_agree"
    --expect-file-contains "cbslith_junction_angle_group_summary.csv:meets_min_group_size,median_ci_low_deg,median_ci_high_deg,n_bootstrap"
    --expect-file-contains "cbslith_junction_angle_group_summary.csv:alpha,lithology_interior,"
    --expect-file-contains "cbslith_junction_angle_comparisons.csv:group_key=lithology_interior"
    --expect-file-contains "cbslith_lithology_units.csv:-1,unmapped,,0,0,0"
    # Integer codes are their own unit names; the south half's 364 pixels below
    # the lowest junctions do not drain through an analysed junction.
    --expect-file-contains "cbslith_lithology_units.csv:1,1,1,24000,2.4,24000"
    --expect-file-contains "cbslith_lithology_units.csv:2,2,2,24000,2.4,23636"
    # The 4 links of no single unit are the group mixed: summarised, never compared.
    --stdout-contains "Link bearings, lithology: 4 link(s) with majority unit 1 are 'mixed (<80% one unit)'"
    --stdout-contains "Note: 4 link record(s) in 'mixed (<80% one unit)' and"
    --expect-file-contains "cbslith_channel_bearing_group_summary.csv:link_bearing,lithology,mixed (<80% one unit),"
    --expect-file-contains "cbslith_junction_angle_group_summary.csv:alpha,lithology_interior,arms disagree,"
    --expect-file-not-contains "cbslith_channel_bearing_comparisons.csv:mixed"
    --expect-file-contains "cbslith_channel_bearing_group_summary.csv:meets_min_group_size,block_order,block_order_mode,min_blocks,n_unblocked,note"
    --expect cbslith_channel_network_block_order_counts.csv
    --stdout-contains "2 lithologic unit(s)"
    --stdout-contains "Bearing group lithology=1 (strahler_order=2) [axial]: n="
    --stdout-contains "Junction-angle group lithology_interior=1 [alpha]: n="
)

# Lithology given but no key uses it: every statistics output equals the
# order-grouped run byte for byte.
lsdtt_add_channel_bearing_cli_test(
    smoke_channel_network_analysis_lithology_observation
    channel_network_analysis_lithology_observation.param.in
    cbsobs_channel_bearing_group_summary.csv
    --expect cbsobs_lithology_units.csv)
if(TEST smoke_channel_network_analysis_lithology_observation)
    lsdtt_set_tests_properties(smoke_channel_network_analysis_lithology_observation PROPERTIES
        FIXTURES_SETUP cna_lithology_observation)
    lsdtt_set_tests_properties(smoke_channel_network_analysis_bearing_statistics_order PROPERTIES
        FIXTURES_SETUP cna_bearing_statistics_order)
    add_test(NAME rigorous_channel_network_analysis_lithology_observation_only
             COMMAND ${CMAKE_COMMAND}
                 "-DPREFIX_A=${CMAKE_CURRENT_BINARY_DIR}/smoke_channel_network_analysis_bearing_statistics_order/out/cbso"
                 "-DPREFIX_B=${CMAKE_CURRENT_BINARY_DIR}/smoke_channel_network_analysis_lithology_observation/out/cbsobs"
                 "-DSUFFIXES=_channel_bearing_group_summary.csv;_channel_bearing_comparisons.csv;_channel_bearing_rose.csv;_channel_bearing_group_membership.csv;_junction_angle_group_summary.csv;_junction_angle_comparisons.csv"
                 -P "${CMAKE_CURRENT_SOURCE_DIR}/cmake/compare_output_sets.cmake")
    lsdtt_set_tests_properties(rigorous_channel_network_analysis_lithology_observation_only PROPERTIES
        LABELS "cli;network"
        FIXTURES_REQUIRED "cna_lithology_observation;cna_bearing_statistics_order")
else()
    lsdtt_mark_test_unregistered(python rigorous_channel_network_analysis_lithology_observation_only)
endif()

lsdtt_add_channel_bearing_cli_test(
    failure_channel_network_analysis_lithology_bad_source
    channel_network_analysis_lithology_bad_source.param.in
    unused.csv --expect-exit-code 1
    --stderr-contains "auto, vector_map or code_raster")

lsdtt_add_channel_bearing_cli_test(
    failure_channel_network_analysis_lithology_key_without_map
    channel_network_analysis_lithology_key_without_map.param.in
    unused.csv --expect-exit-code 1
    --stderr-contains "lithology_fname is required when channel_bearing_group_by or channel_bearing_stratify_by is lithology")

# Lithology and region maps that label no valid DEM pixel (200 km east, same
# CRS): refused, naming both extents and CRSs.
foreach(_cbs_far lithology_no_overlap_raster lithology_no_overlap_vector
                 region_no_overlap_raster region_no_overlap_vector)
    lsdtt_add_channel_bearing_cli_test(
        failure_channel_network_analysis_${_cbs_far}
        channel_network_analysis_${_cbs_far}.param.in
        unused.csv --expect-exit-code 1
        --stderr-contains "labels no valid DEM pixel: the map covers x 913650.1 to 915250.1"
        --stderr-contains "the DEM covers x 713650.1 to 715250.1, y 8711785.5 to 8714785.5 (WGS 84 / UTM zone 58S)")
endforeach()

# Spatially varying erodibility (K): kernels, steady state, values, grid modes.
add_executable(lsdtt-test-lem-erodibility cpp/test_lem_erodibility.cpp)
target_link_libraries(lsdtt-test-lem-erodibility PRIVATE lsdtt3_cpu)
add_test(NAME unit_lem_erodibility
         COMMAND lsdtt-test-lem-erodibility --suite fast
                 "${CMAKE_CURRENT_BINARY_DIR}/lem_erodibility_work")
lsdtt_set_tests_properties(unit_lem_erodibility PROPERTIES LABELS "lem")
# The analytic steady states (strip, two catchments, red noise, snap): 2-8 s.
add_test(NAME rigorous_lem_erodibility_analytic
         COMMAND lsdtt-test-lem-erodibility --suite analytic)
lsdtt_set_tests_properties(rigorous_lem_erodibility_analytic PROPERTIES LABELS "lem")
lsdtt_set_tests_properties(unit_lem_channel_coupling PROPERTIES LABELS "lem;hillslope")

# The Newton solver for the Roering nonlinear hillslope (steady + transient).
add_executable(lsdtt-test-hillslope-newton cpp/test_hillslope_newton.cpp)
target_link_libraries(lsdtt-test-hillslope-newton PRIVATE lsdtt3_cpu)

add_test(NAME unit_hillslope_newton COMMAND lsdtt-test-hillslope-newton)
lsdtt_set_tests_properties(unit_hillslope_newton PROPERTIES LABELS "hillslope")

# Independent equation-faithful Perron (2011) Q-imp reference solver.
add_executable(lsdtt-test-perron-qimp cpp/test_perron_qimp.cpp)
target_link_libraries(lsdtt-test-perron-qimp PRIVATE lsdtt3_cpu)

add_test(NAME unit_hillslope_perron_qimp COMMAND lsdtt-test-perron-qimp)
lsdtt_set_tests_properties(unit_hillslope_perron_qimp PROPERTIES LABELS "hillslope")

# Independent C++ port of Ren et al. (2026) modular implicit method.
add_executable(lsdtt-test-ren-modular-implicit cpp/test_ren_modular_implicit.cpp)
target_link_libraries(lsdtt-test-ren-modular-implicit PRIVATE lsdtt3_cpu)

add_test(NAME unit_hillslope_ren_modular_implicit
         COMMAND lsdtt-test-ren-modular-implicit)
lsdtt_set_tests_properties(unit_hillslope_ren_modular_implicit PROPERTIES
                     LABELS "hillslope")

# Heavy P1/P2 literature comparison driver; run explicitly with an output path.
add_executable(lsdtt-perron-ren-p1-p2 cpp/perron_ren_p1_p2_driver.cpp)
target_link_libraries(lsdtt-perron-ren-p1-p2 PRIVATE lsdtt3_cpu)

# Heavy P3 equilibrium timestep sweep; run explicitly with an output CSV.
add_executable(lsdtt-perron-ren-p3 cpp/perron_ren_p3_driver.cpp)
target_link_libraries(lsdtt-perron-ren-p3 PRIVATE lsdtt3_cpu)

# Heavy P4 transient timestep sweep; `full-p2` enables the multi-hour 1 yr 2-D reference.
add_executable(lsdtt-perron-ren-p4 cpp/perron_ren_p4_driver.cpp)
target_link_libraries(lsdtt-perron-ren-p4 PRIVATE lsdtt3_cpu)

# P5 one-step stability/robustness envelope across reference and production solvers.
add_executable(lsdtt-perron-ren-p5 cpp/perron_ren_p5_driver.cpp)
target_link_libraries(lsdtt-perron-ren-p5 PRIVATE lsdtt3_cpu)

# Heavy P6 spatial refinement; cellular canonical runs stop at 100x100 by design.
add_executable(lsdtt-perron-ren-p6 cpp/perron_ren_p6_driver.cpp)
target_link_libraries(lsdtt-perron-ren-p6 PRIVATE lsdtt3_cpu)

# Small Newton-only references for the cellular hillslope investigation. The
# normal suite uses narrow strips and a 32x32 OCN; larger OCN/scaling cases are
# performance fixtures and deliberately do not run here.
add_executable(lsdtt-test-cellular-hillslope-references
               cpp/test_cellular_hillslope_references.cpp)
target_link_libraries(lsdtt-test-cellular-hillslope-references PRIVATE lsdtt3_cpu)

add_test(NAME unit_cellular_hillslope_references
         COMMAND lsdtt-test-cellular-hillslope-references
                 ${CMAKE_CURRENT_BINARY_DIR}/cellular_hillslope_reference_output)
lsdtt_set_tests_properties(unit_cellular_hillslope_references PROPERTIES
                     LABELS "hillslope;cellular")
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/cellular_hillslope_reference_output)

# Soft-critical constitutive law and whole-grid pass-cost envelope. This is P1:
# test-local mathematics only, before a raster transport solver exists.
add_executable(lsdtt-test-cellular-hillslope-constitutive
               cpp/test_cellular_hillslope_constitutive.cpp)

add_test(NAME unit_cellular_hillslope_constitutive
         COMMAND lsdtt-test-cellular-hillslope-constitutive
                 ${CMAKE_CURRENT_BINARY_DIR}/cellular_hillslope_reference_output)
lsdtt_set_tests_properties(unit_cellular_hillslope_constitutive PROPERTIES
                     LABELS "hillslope;cellular")

# P2 global deterministic cellular oracle: Newton agreement, conservation,
# periodicity, ridgetop curvature, supercritical relaxation and small-OCN modes.
add_executable(lsdtt-test-cellular-hillslope-solver
               cpp/test_cellular_hillslope_solver.cpp)
target_link_libraries(lsdtt-test-cellular-hillslope-solver PRIVATE lsdtt3_cpu)

# Research-only E2 diagnostic, built but not registered with CTest because it
# consumes user-supplied rasters rather than a committed fixture.
add_executable(cellular-hillslope-face-stats
               cpp/cellular_hillslope_face_stats.cpp)
target_link_libraries(cellular-hillslope-face-stats PRIVATE lsdtt3_cpu)

add_executable(cellular-hillslope-step-benchmark
               cpp/cellular_hillslope_step_benchmark.cpp)
target_link_libraries(cellular-hillslope-step-benchmark PRIVATE lsdtt3_cpu)

add_executable(cellular-hillslope-kernel-benchmark
               cpp/cellular_hillslope_kernel_benchmark.cpp)
target_link_libraries(cellular-hillslope-kernel-benchmark PRIVATE lsdtt3_cpu)

add_executable(cellular-hillslope-raster-compare
               cpp/cellular_hillslope_raster_compare.cpp)
target_link_libraries(cellular-hillslope-raster-compare PRIVATE lsdtt3_cpu)

add_executable(cellular-hillslope-common-limit
               cpp/cellular_hillslope_common_limit.cpp)
target_link_libraries(cellular-hillslope-common-limit PRIVATE lsdtt3_cpu)

add_executable(cellular-hillslope-grid-scaling
               cpp/cellular_hillslope_grid_scaling.cpp)
target_link_libraries(cellular-hillslope-grid-scaling PRIVATE lsdtt3_cpu)

add_executable(cellular-hillslope-nonlocal-benchmark
               cpp/cellular_hillslope_nonlocal_benchmark.cpp)
target_link_libraries(cellular-hillslope-nonlocal-benchmark PRIVATE lsdtt3_cpu)

add_executable(cellular-hillslope-nonlocal-lem-benchmark
               cpp/cellular_hillslope_nonlocal_lem_benchmark.cpp)
target_link_libraries(cellular-hillslope-nonlocal-lem-benchmark PRIVATE lsdtt3_cpu)

add_executable(cellular-hillslope-particle-benchmark
               cpp/cellular_hillslope_particle_benchmark.cpp)
target_link_libraries(cellular-hillslope-particle-benchmark PRIVATE lsdtt3_cpu)

add_executable(lsdtt-test-hillslope-particle-io
               cpp/test_hillslope_particle_io.cpp)
target_link_libraries(lsdtt-test-hillslope-particle-io PRIVATE lsdtt3_cpu)
add_test(NAME unit_hillslope_particle_io
         COMMAND lsdtt-test-hillslope-particle-io
                 ${CMAKE_CURRENT_BINARY_DIR}/cellular_hillslope_reference_output)
lsdtt_set_tests_properties(unit_hillslope_particle_io PROPERTIES
                     LABELS "hillslope;cellular")

add_test(NAME unit_cellular_hillslope_common_limit
         COMMAND cellular-hillslope-common-limit)
lsdtt_set_tests_properties(unit_cellular_hillslope_common_limit PROPERTIES
                     LABELS "hillslope;cellular")

add_test(NAME unit_cellular_hillslope_solver
         COMMAND lsdtt-test-cellular-hillslope-solver)
lsdtt_set_tests_properties(unit_cellular_hillslope_solver PROPERTIES
                     LABELS "hillslope;cellular")

# The SOR steady-state hillslope solver: conservation, precision floor, 'o' vs 'n',
# periodic invariance, and recovery of U/D from known forcing.
add_executable(lsdtt-test-hillslope-steady-state cpp/test_hillslope_steady_state.cpp)
target_link_libraries(lsdtt-test-hillslope-steady-state PRIVATE lsdtt3_cpu)

add_test(NAME unit_hillslope_steady_state COMMAND lsdtt-test-hillslope-steady-state)
lsdtt_set_tests_properties(unit_hillslope_steady_state PROPERTIES LABELS "hillslope")

# chi at the channel head, as distinct from domain-wide MaxChi.
add_executable(lsdtt-test-chi-channel-head cpp/test_chi_channel_head.cpp)
target_link_libraries(lsdtt-test-chi-channel-head PRIVATE lsdtt3_cpu)

add_test(NAME unit_chi_channel_head COMMAND lsdtt-test-chi-channel-head)
lsdtt_set_tests_properties(unit_chi_channel_head PROPERTIES LABELS "lem")

# chi_c against domain-wide MaxChi: the size of the relief
# error relief-mode carries for a landscape that will have hillslopes.
add_executable(chi-c-ratio-driver cpp/chi_c_ratio_driver.cpp)
target_link_libraries(chi-c-ratio-driver PRIVATE lsdtt3_cpu)

# Writes the D8
# distance-to-channel field so it can be compared against the streamline
# tracer's Lh on the same source pixels. A driver, not a CTest entry.
add_executable(hillslope-length-comparison-driver cpp/hillslope_length_comparison_driver.cpp)
target_link_libraries(hillslope-length-comparison-driver PRIVATE lsdtt3_cpu)

add_executable(hilltop-workflow-driver cpp/hilltop_workflow_driver.cpp)
target_link_libraries(hilltop-workflow-driver PRIVATE lsdtt3_cpu)

# Phase 3 Step 0: measures the streamline reuse ratio and the per-trace cost.
# Reports rather than asserts - the reuse ratio is landscape-dependent, so it is
# not a pass/fail quantity. Built but not registered with CTest.
add_executable(bench-streamlines cpp/bench_streamlines.cpp)
target_link_libraries(bench-streamlines PRIVATE lsdtt3_cpu)

# Implicit (backward-Euler / TR-BDF2) Roering nonlinear hillslope solver
# (docs/developer_notes/testing_and_benchmarking.md).  Rectangular-domain validation against the
# analytic Roering steady profile, the SOR steady-state solver, and the
# conservation / Sc-enforcement / timestep-independence properties.  Synthetic
# surfaces only: no network, no committed DEM, so it runs in the default suite.
# It writes its plottable profile / grid-refinement CSVs into the test working
# directory (consumed by lsdviztools.lsdtransientvalidation).
#
# The Phase 4 linear-solve TIMING sweep in the same driver is NOT part of this
# test -- it is reported, never asserted, and is registered separately as
# performance_hillslope_implicit below.
add_executable(lsdtt-test-hillslope-implicit cpp/test_hillslope_implicit.cpp)
target_link_libraries(lsdtt-test-hillslope-implicit PRIVATE lsdtt3_cpu)

add_test(NAME unit_hillslope_implicit
         COMMAND lsdtt-test-hillslope-implicit
                 ${CMAKE_CURRENT_BINARY_DIR}/hillslope_implicit_output)
lsdtt_set_tests_properties(unit_hillslope_implicit PROPERTIES LABELS "hillslope")
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/hillslope_implicit_output)

# -----------------------------------------------------------------------------
# lsdtt-fetch-raster CLI (docs/developer_notes/testing_and_benchmarking.md)
#
# These run offline: they exercise the argument handling and validation that
# happen before any network access is attempted.
# -----------------------------------------------------------------------------
add_lsdtt_cli_smoke_test(
    failure_fetch_raster_missing_bbox
    lsdtt-fetch-raster
    fetch_raster_missing_bbox.param.in
    unused_expected_output
    --expect-failure
    --stderr-contains "FATAL: no area given"
)

add_lsdtt_cli_smoke_test(
    failure_fetch_raster_partial_bbox
    lsdtt-fetch-raster
    fetch_raster_partial_bbox.param.in
    unused_expected_output
    --expect-failure
    --stderr-contains "FATAL: north must be set"
)

add_lsdtt_cli_smoke_test(
    failure_fetch_raster_bad_source
    lsdtt-fetch-raster
    fetch_raster_bad_source.param.in
    unused_expected_output
    --expect-failure
    --stderr-contains "FATAL: unknown dem_source"
)

add_lsdtt_cli_smoke_test(
    failure_fetch_raster_bad_bounds_mode
    lsdtt-fetch-raster
    fetch_raster_bad_bounds_mode.param.in
    unused_expected_output
    --expect-failure
    --stderr-contains "FATAL: bounds_mode must be one of"
)

add_lsdtt_cli_smoke_test(
    failure_fetch_raster_outlets_missing_file
    lsdtt-fetch-raster
    fetch_raster_outlets_missing_file.param.in
    unused_expected_output
    --expect-failure
    --stderr-contains "FATAL: bounds_mode=outlets requires outlet_points_fname"
)

add_lsdtt_cli_smoke_test(
    failure_fetch_raster_outlets_conflict_area
    lsdtt-fetch-raster
    fetch_raster_outlets_conflict_area.param.in
    unused_expected_output
    --expect-failure
    --stderr-contains "FATAL: bounds_mode=outlets forbids west/south/east/north and bounds_fname"
)

add_lsdtt_cli_smoke_test(
    failure_fetch_raster_bad_prefetch_source
    lsdtt-fetch-raster
    fetch_raster_bad_prefetch_source.param.in
    unused_expected_output
    --expect-failure
    --stderr-contains "FATAL: unknown prefetch_source"
)

add_lsdtt_cli_smoke_test(
    failure_fetch_raster_bad_auth_fallback_source
    lsdtt-fetch-raster
    fetch_raster_bad_auth_fallback_source.param.in
    unused_expected_output
    --expect-failure
    --stderr-contains "FATAL: auth_fallback_source must be keyless"
)

add_lsdtt_cli_smoke_test(
    failure_fetch_raster_bad_nested_stop_level
    lsdtt-fetch-raster
    fetch_raster_bad_nested_stop_level.param.in
    unused_expected_output
    --expect-failure
    --stderr-contains "FATAL: basin_nested_stop_level must be one of: coarse, 90m, 30m"
)

add_lsdtt_cli_smoke_test(
    failure_fetch_raster_outlets_missing_dataset
    lsdtt-fetch-raster
    fetch_raster_outlets_not_implemented.param.in
    unused_expected_output
    --expect-failure
    --stderr-contains "FATAL: could not open outlet points CSV"
)

add_lsdtt_cli_smoke_test(
    failure_fetch_raster_ot_missing_key_strict
    lsdtt-fetch-raster
    fetch_raster_ot_missing_key_strict.param.in
    unused_expected_output
    --expect-failure
    --stderr-contains "FATAL: no credential found"
)

add_lsdtt_cli_smoke_test(
    failure_fetch_raster_ot_missing_key_fallback
    lsdtt-fetch-raster
    fetch_raster_ot_missing_key_fallback.param.in
    unused_expected_output
    --expect-failure
    --stderr-contains "FATAL: lsdtt::fetch_dem: west must be less than east"
)

# -h must work without a parameter file, and without touching the network.
# Fetching along a linear feature into an oblique frame
# (fetch_from_linear_feature): every other way of giving an area is refused,
# and so are a missing feature file and a line without a width, all before
# any network access. green_river_line.geojson is a three-vertex WGS84 line.
foreach(_ff_case IN ITEMS
        "with_bbox|cannot be combined with west/south/east/north"
        "with_bounds_file|cannot be combined with bounds_fname"
        "with_outlets|cannot be combined with bounds_mode vector or outlets"
        "with_size_pixels|cannot be combined with size_pixels"
        "missing_fname|needs fetch_feature_fname"
        "no_width|a line or two points needs a width > 0")
    string(REPLACE "|" ";" _ff_parts "${_ff_case}")
    list(GET _ff_parts 0 _ff_name)
    list(GET _ff_parts 1 _ff_message)
    add_lsdtt_cli_smoke_test(
        failure_fetch_raster_feature_${_ff_name}
        lsdtt-fetch-raster
        fetch_raster_feature_${_ff_name}.param.in
        unused_expected_output
        --expect-failure
        --stderr-contains "${_ff_message}"
    )
endforeach()

# Catalogue sources (STAC and WCS) and the size limits: refused before any
# request, so these run offline.
add_lsdtt_cli_smoke_test(
    failure_fetch_raster_3dep_outside_coverage
    lsdtt-fetch-raster
    fetch_raster_3dep_outside_coverage.param.in
    unused_expected_output
    --expect-failure
    --stderr-contains "outside us_3dep coverage"
)

add_lsdtt_cli_smoke_test(
    failure_fetch_raster_3dep_bad_resolution
    lsdtt-fetch-raster
    fetch_raster_3dep_bad_resolution.param.in
    unused_expected_output
    --expect-failure
    --stderr-contains "us_3dep offers 10, 30 m"
)

add_lsdtt_cli_smoke_test(
    failure_fetch_raster_3dep_surface_model
    lsdtt-fetch-raster
    fetch_raster_3dep_surface_model.param.in
    unused_expected_output
    --expect-failure
    --stderr-contains "us_3dep has no dsm product"
)

add_lsdtt_cli_smoke_test(
    failure_fetch_raster_arcticdem_outside_coverage
    lsdtt-fetch-raster
    fetch_raster_arcticdem_outside_coverage.param.in
    unused_expected_output
    --expect-failure
    --stderr-contains "outside arctic_arcticdem coverage"
)

add_lsdtt_cli_smoke_test(
    failure_fetch_raster_rema_terrain_model
    lsdtt-fetch-raster
    fetch_raster_rema_terrain_model.param.in
    unused_expected_output
    --expect-failure
    --stderr-contains "antarctic_rema has no dtm product"
)

add_lsdtt_cli_smoke_test(
    failure_fetch_raster_wcs_outside_coverage
    lsdtt-fetch-raster
    fetch_raster_wcs_outside_coverage.param.in
    unused_expected_output
    --expect-failure
    --stderr-contains "outside no_dtm coverage"
)

add_lsdtt_cli_smoke_test(
    failure_fetch_raster_wcs_too_many_native_cells
    lsdtt-fetch-raster
    fetch_raster_wcs_too_many_native_cells.param.in
    unused_expected_output
    --expect-failure
    --stderr-contains "native 0.5 m cells"
)

add_lsdtt_cli_smoke_test(
    failure_fetch_raster_swissalti3d_surface_model
    lsdtt-fetch-raster
    fetch_raster_swissalti3d_surface_model.param.in
    unused_expected_output
    --expect-failure
    --stderr-contains "ch_swissalti3d has no dsm product"
)

add_lsdtt_cli_smoke_test(
    failure_fetch_raster_hrdem_bad_resolution
    lsdtt-fetch-raster
    fetch_raster_hrdem_bad_resolution.param.in
    unused_expected_output
    --expect-failure
    --stderr-contains "ca_hrdem offers 1, 2 m"
)

add_lsdtt_cli_smoke_test(
    failure_fetch_raster_sct_lidar_outside_coverage
    lsdtt-fetch-raster
    fetch_raster_sct_lidar_outside_coverage.param.in
    unused_expected_output
    --expect-failure
    --stderr-contains "outside gb_sct_lidar coverage"
)

add_lsdtt_cli_smoke_test(
    failure_fetch_raster_nz_lidar_lerc_warning
    lsdtt-fetch-raster
    fetch_raster_nz_lidar_lerc_warning.param.in
    unused_expected_output
    --expect-failure
    --stderr-contains "LERC-compressed"
    --stderr-contains "OpenTopography (opentopography.org) with an API key"
    --stderr-contains "outside nz_lidar coverage"
)

add_lsdtt_cli_smoke_test(
    failure_fetch_raster_gsi_bad_resolution
    lsdtt-fetch-raster
    fetch_raster_gsi_bad_resolution.param.in
    unused_expected_output
    --expect-failure
    --stderr-contains "jp_gsi offers 1, 5, 10 m"
)

add_lsdtt_cli_smoke_test(
    failure_fetch_raster_sct_lidar_bad_collection
    lsdtt-fetch-raster
    fetch_raster_sct_lidar_bad_collection.param.in
    unused_expected_output
    --expect-failure
    --stderr-contains "gb_sct_lidar has no collection 'phase-9'"
)

add_lsdtt_cli_smoke_test(
    failure_fetch_raster_copernicus_surface_model
    lsdtt-fetch-raster
    fetch_raster_copernicus_surface_model.param.in
    unused_expected_output
    --expect-failure
    --stderr-contains "Copernicus GLO-30 is a DSM only"
)

add_lsdtt_cli_smoke_test(
    failure_fetch_raster_too_many_pixels
    lsdtt-fetch-raster
    fetch_raster_too_many_pixels.param.in
    unused_expected_output
    --expect-failure
    --stderr-contains "more than max_output_pixels"
)

add_lsdtt_cli_smoke_test(
    failure_fetch_raster_hires_outlets
    lsdtt-fetch-raster
    fetch_raster_hires_outlets.param.in
    unused_expected_output
    --expect-failure
    --stderr-contains "cannot be used with bounds_mode: outlets"
)

add_test(NAME smoke_fetch_raster_help_dem_source
         COMMAND lsdtt-fetch-raster -h dem_source)
lsdtt_set_tests_properties(smoke_fetch_raster_help_dem_source PROPERTIES
    LABELS "cli;fetch"
    PASS_REGULAR_EXPRESSION "USGS 3DEP seamless DEM -- a DIGITAL TERRAIN MODEL")

add_test(NAME smoke_fetch_raster_help_dem_source_polar
         COMMAND lsdtt-fetch-raster -h dem_source)
lsdtt_set_tests_properties(smoke_fetch_raster_help_dem_source_polar PROPERTIES
    LABELS "cli;fetch"
    PASS_REGULAR_EXPRESSION "ArcticDEM mosaic v4.1 -- a DIGITAL SURFACE MODEL.*REMA mosaic v2.0 -- a DIGITAL SURFACE MODEL")

add_test(NAME smoke_fetch_raster_help_dem_source_wcs
         COMMAND lsdtt-fetch-raster -h dem_source)
lsdtt_set_tests_properties(smoke_fetch_raster_help_dem_source_wcs PROPERTIES
    LABELS "cli;fetch"
    PASS_REGULAR_EXPRESSION "no_dtm.*nl_ahn.*DIGITAL SURFACE MODEL \\(dsm\\).*es_mdt.*pl_nmt")

add_test(NAME smoke_fetch_raster_help_dem_source_hrdem_swiss
         COMMAND lsdtt-fetch-raster -h dem_source)
lsdtt_set_tests_properties(smoke_fetch_raster_help_dem_source_hrdem_swiss PROPERTIES
    LABELS "cli;fetch"
    PASS_REGULAR_EXPRESSION "ca_hrdem.*CGVD2013.*ch_swissalti3d.*latest acquisition of each tile")

add_test(NAME smoke_fetch_raster_help_dem_source_grid_keyed
         COMMAND lsdtt-fetch-raster -h dem_source)
lsdtt_set_tests_properties(smoke_fetch_raster_help_dem_source_grid_keyed PROPERTIES
    LABELS "cli;fetch"
    PASS_REGULAR_EXPRESSION "nz_lidar.*NZVD2016.*LERC.*gb_sct_lidar.*source_collection_choice")

add_test(NAME smoke_fetch_raster_help_dem_source_gsi
         COMMAND lsdtt-fetch-raster -h dem_source)
lsdtt_set_tests_properties(smoke_fetch_raster_help_dem_source_gsi PROPERTIES
    LABELS "cli;fetch"
    PASS_REGULAR_EXPRESSION "jp_gsi.*DIGITAL TERRAIN MODEL.*second.*resampling")

add_test(NAME smoke_fetch_raster_help_linear_feature
         COMMAND lsdtt-fetch-raster -h fetch_from_linear_feature)
lsdtt_set_tests_properties(smoke_fetch_raster_help_linear_feature PROPERTIES
    LABELS "cli;fetch"
    PASS_REGULAR_EXPRESSION "minimum-area")

add_test(NAME smoke_fetch_raster_help
         COMMAND lsdtt-fetch-raster -h)
lsdtt_set_tests_properties(smoke_fetch_raster_help PROPERTIES
    LABELS "cli;fetch"
    PASS_REGULAR_EXPRESSION "Bounding box")

# The CLI must never grow a credential parameter: parameter files get committed
# and print_parameters() echoes every parameter at the end of a run. This guards
# the rule stated in docs/developer_notes/testing_and_benchmarking.md and in the CLI's own header.
add_test(NAME rigorous_fetch_raster_no_credential_parameters
         COMMAND lsdtt-fetch-raster -h)
lsdtt_set_tests_properties(rigorous_fetch_raster_no_credential_parameters PROPERTIES
    LABELS "cli;fetch"
    FAIL_REGULAR_EXPRESSION "api_key;API_Key;password;secret;token;credential")

add_lsdtt_cli_smoke_test(
    failure_fill_vincent_needs_zero_slope
    lsdtt-raster-preprocessing
    fill_vincent_needs_zero_slope.param.in
    unused_expected_output
    --expect-failure
    --stderr-contains "cannot impose a gradient across flats"
    --stderr-contains "or set fill_flat_resolution_choice: barnes2014b_gradient"
)

# fill_flat_resolution_choice (Barnes et al. 2014b flat resolution) in
# lsdtt-raster-preprocessing. The resolver numerics are covered by
# unit_flat_resolution; these carry the CLI wiring: both options run and
# write, the three refusals are exact, the diagnostics reach stdout, and the
# rim-pit warning reaches stderr. tiny_dem.asc has no flats once filled, so
# the first two are plumbing only; flat_valley_dem.asc (11 x 24 valley, one
# flat floor) is where the options differ from the default fill.
add_lsdtt_cli_smoke_test(
    smoke_fill_flat_resolution_float_step
    lsdtt-raster-preprocessing
    fill_flat_resolution_float_step.param.in
    tinyfs_fill.tif
)

add_lsdtt_cli_smoke_test(
    smoke_fill_flat_resolution_gradient
    lsdtt-raster-preprocessing
    fill_flat_resolution_gradient.param.in
    tinygr_fill.tif
)

add_lsdtt_cli_smoke_test(
    failure_fill_flat_resolution_unknown_choice
    lsdtt-raster-preprocessing
    fill_flat_resolution_unknown.param.in
    unused_expected_output
    --expect-failure
    --stderr-contains "FATAL: fill_flat_resolution_choice must be one of inline_gradient, barnes2014b_gradient, barnes2014b_float_step\; got 'barnes2014b_gradiant'"
)

add_lsdtt_cli_smoke_test(
    failure_fill_flat_resolution_gradient_zero_slope
    lsdtt-raster-preprocessing
    fill_flat_resolution_gradient_zero_slope.param.in
    unused_expected_output
    --expect-failure
    --stderr-contains "FATAL: fill_flat_resolution_choice = barnes2014b_gradient needs fill_min_slope > 0 (it sets the gradient imposed across flats)\; got 0"
)

add_lsdtt_cli_smoke_test(
    failure_fill_flat_resolution_fill_spill_merge
    lsdtt-raster-preprocessing
    fill_flat_resolution_fsm.param.in
    unused_expected_output
    --expect-failure
    --stderr-contains "FATAL: fill_flat_resolution_choice = barnes2014b_float_step applies to the barnes2014 and vincent1993 fills, not barnesetal2021"
)

# Bands on the fill difference, not oracles of the whole surface: the largest
# raise on the valley is 48 x 2^-17 m (float step) and 48 u_q = 0.0681 m
# (gradient), against 0.0255 m for the default fill, so each band also fails
# if the option silently falls back to the default.
add_lsdtt_cli_smoke_test(
    rigorous_fill_flat_resolution_valley_float_step
    lsdtt-raster-preprocessing
    fill_flat_resolution_valley_float_step.param.in
    valleyfs_fill.tif
    --expect valleyfs_fill_diff.tif
    --expect-raster-stat valleyfs_fill_diff.tif:max:0.000366:0.000367
    --stdout-contains "Flat regions: 1, flat cells: 198, NO_FLOW cells: 195"
    --stdout-contains "Max mask: 48"
    --stdout-contains "Rim cells without descent: 0"
    --stdout-contains "Rounds: 1 ("
    --stdout-contains "NOTE: fill_min_slope (0.0001) is ignored"
)

add_lsdtt_cli_smoke_test(
    rigorous_fill_flat_resolution_valley_gradient
    lsdtt-raster-preprocessing
    fill_flat_resolution_valley_gradient.param.in
    valleygr_fill.tif
    --expect valleygr_fill_diff.tif
    --expect-raster-stat valleygr_fill_diff.tif:max:0.0681:0.0682
    --stdout-contains "Capped flats: 0"
    --stdout-contains "Float-step fallbacks: 0"
    --stdout-contains "Min effective gradient: 0.000141907 m/m"
    --stdout-contains "Rounds: 1 ("
)

# The knob 1 ULP above the valley floor: round 1 (RichDEM) leaves it without a
# lower neighbour, round 2 refills and re-resolves it, and the result routes.
# Pinned on the first run (2026-10-09).
add_lsdtt_cli_smoke_test(
    rigorous_fill_flat_resolution_knob_rounds
    lsdtt-raster-preprocessing
    fill_flat_resolution_rim_warning.param.in
    valleyknob_fill.tif
    --stdout-contains "Rounds: 2 (cells without descent after each round: 1, 0)"
    --stdout-contains "Round 1 (RichDEM-comparable): 8 pair(s), 8 cell(s) raised above surroundings, 1 cell(s) without descent"
    --stdout-contains "Rim cells without descent: 0"
    --stdout-contains "Final priority-flood pass: not needed"
    --stderr-not-contains "WARNING:"
)

lsdtt_set_tests_properties(
    smoke_fill_flat_resolution_float_step
    smoke_fill_flat_resolution_gradient
    PROPERTIES LABELS "cli;flow")
lsdtt_set_tests_properties(
    failure_fill_flat_resolution_unknown_choice
    failure_fill_flat_resolution_gradient_zero_slope
    failure_fill_flat_resolution_fill_spill_merge
    PROPERTIES LABELS "cli;flow")
lsdtt_set_tests_properties(
    rigorous_fill_flat_resolution_valley_float_step
    rigorous_fill_flat_resolution_valley_gradient
    rigorous_fill_flat_resolution_knob_rounds
    PROPERTIES LABELS "cli;flow")

# fill_flat_resolution_choice in the routing CLIs. The parameter check, the
# fill-and-resolve path, the raster_is_filled: true resolve-only path and the
# report are one shared implementation
# (src/cpu/lsdtt_fill_flat_resolution_reporting_cpu), so one CLI carries the
# raster_is_filled cases and each CLI family gets one run with an option on
# the valley. The "Flat regions" line shows the resolver ran on the valley's
# floor; under the default fill that line is never printed.
add_lsdtt_cli_smoke_test(
    smoke_chi_analysis_flat_resolution_gradient
    lsdtt-chi-analysis
    chi_analysis_flat_resolution_gradient.param.in
    chivg_chi_fit_nodes_mn_0p500.csv
    --stdout-contains "Flat regions: 1, flat cells: 198, NO_FLOW cells: 195"
    --stdout-contains "Min effective gradient: 0.000141907 m/m"
)

# raster_is_filled: true (the input is not filled again). The valley is a
# filled surface with a true flat: a barnes2014b option resolves it and the
# run succeeds; inline_gradient routes it as it is and the router refuses it.
add_lsdtt_cli_smoke_test(
    smoke_chi_analysis_flat_resolution_filled_float_step
    lsdtt-chi-analysis
    chi_analysis_flat_resolution_filled_float_step.param.in
    chivf_chi_fit_nodes_mn_0p500.csv
    --stdout-contains "Raster is marked as already filled. Skipping fill step."
    --stdout-contains "Flat regions: 1, flat cells: 198, NO_FLOW cells: 195"
    --stdout-contains "Max raise: 0.000366211 m"
)

add_lsdtt_cli_smoke_test(
    failure_chi_analysis_flat_resolution_filled_inline
    lsdtt-chi-analysis
    chi_analysis_flat_resolution_filled_inline.param.in
    unused_expected_output
    --expect-failure
    --stderr-contains "could not assign a receiver"
)

add_lsdtt_cli_smoke_test(
    failure_chi_analysis_flat_resolution_unfilled_input
    lsdtt-chi-analysis
    chi_analysis_flat_resolution_unfilled_input.param.in
    unused_expected_output
    --expect-failure
    --stderr-contains "FATAL: flat resolution failed: lsdtt::compute_barnes2014b_flat_mask: 1 flat cell(s) in 1 flat region(s) have no outlet"
    --stderr-contains "\; set raster_is_filled: false"
)

add_lsdtt_cli_smoke_test(
    smoke_channel_extraction_flat_resolution_gradient
    lsdtt-channel-extraction
    channel_extraction_flat_resolution_gradient.param.in
    chanvg_channel_network_lines.fgb
    --expect-vector-feature-min chanvg_channel_network_lines.fgb:1
    --stdout-contains "Flat regions: 1, flat cells: 198, NO_FLOW cells: 195"
    --stdout-contains "Min effective gradient: 0.000141907 m/m"
)

add_lsdtt_cli_smoke_test(
    smoke_drainage_analysis_flat_resolution_float_step
    lsdtt-drainage-analysis
    drainage_analysis_flat_resolution_float_step.param.in
    drainvf_d8_drainage_area.tif
    --stdout-contains "Flat regions: 1, flat cells: 198, NO_FLOW cells: 195"
    --stdout-contains "Max raise: 0.000366211 m"
)

# Two CLIs call their routing fill slope fill_min_slope_for_routing; the
# gradient refusal names the parameter the user actually set.
add_lsdtt_cli_smoke_test(
    failure_fluvial_surfaces_flat_resolution_gradient_zero_slope
    lsdtt-fluvial-surfaces
    fluvial_surfaces_flat_resolution_gradient_zero_slope.param.in
    unused_expected_output
    --expect-failure
    --stderr-contains "FATAL: fill_flat_resolution_choice = barnes2014b_gradient needs fill_min_slope_for_routing > 0 (it sets the gradient imposed across flats)\; got 0"
)

lsdtt_set_tests_properties(
    smoke_chi_analysis_flat_resolution_gradient
    smoke_chi_analysis_flat_resolution_filled_float_step
    smoke_channel_extraction_flat_resolution_gradient
    smoke_drainage_analysis_flat_resolution_float_step
    PROPERTIES LABELS "cli;flow")
lsdtt_set_tests_properties(
    failure_chi_analysis_flat_resolution_filled_inline
    failure_chi_analysis_flat_resolution_unfilled_input
    failure_fluvial_surfaces_flat_resolution_gradient_zero_slope
    PROPERTIES LABELS "cli;flow")

# Real-terrain gate. The committed 128 x 128 Copernicus GLO-30 lakeland crop
# (tests/data/flat_resolution_lakeland_crop_dem.tif; provenance in
# tests/README.md) is a surface on which a single Barnes pass -- RichDEM
# ResolveFlats -- leaves rim cells without a lower neighbour, so before the
# round loop these exact commands exited 1 with "could not assign a receiver".
# Each barnes2014b run must now exit 0 with no cell left without descent; the
# round and final-pass lines are pinned (2026-10-10). The synthetic valley
# fixtures above strand nothing, which is why they missed that failure.
# unit_flat_resolution_real_terrain checks the same crop at API level.
add_lsdtt_cli_smoke_test(
    rigorous_drainage_analysis_flat_resolution_real_terrain_float_step
    lsdtt-drainage-analysis
    drainage_analysis_flat_resolution_real_terrain_float_step.param.in
    rtfs_d8_drainage_area.tif
    --stdout-contains "Rim cells without descent: 0"
    --stdout-contains "Rounds: 4 (cells without descent after each round: 117, 44, 26, 28)"
    --stdout-contains "Round 1 (RichDEM-comparable): 364 pair(s), 191 cell(s) raised above surroundings, 117 cell(s) without descent"
    --stdout-contains "Final priority-flood pass: 107 cell(s)"
    --stderr-not-contains "could not assign a receiver"
)

add_lsdtt_cli_smoke_test(
    rigorous_drainage_analysis_flat_resolution_real_terrain_gradient
    lsdtt-drainage-analysis
    drainage_analysis_flat_resolution_real_terrain_gradient.param.in
    rtgr_d8_drainage_area.tif
    --stdout-contains "Rim cells without descent: 0"
    --stdout-contains "Rounds: 4 (cells without descent after each round: 115, 56, 43, 30)"
    --stdout-contains "Round 1 (RichDEM-comparable): 362 pair(s), 189 cell(s) raised above surroundings, 115 cell(s) without descent"
    --stdout-contains "Final priority-flood pass: 108 cell(s)"
    --stderr-not-contains "could not assign a receiver"
)

# Control: the default flat handling routes the same crop with the same
# command shape (no "Flat resolution:" block on this path).
add_lsdtt_cli_smoke_test(
    smoke_drainage_analysis_flat_resolution_real_terrain_inline_gradient
    lsdtt-drainage-analysis
    drainage_analysis_flat_resolution_real_terrain_inline_gradient.param.in
    rtig_d8_drainage_area.tif
)

lsdtt_set_tests_properties(
    rigorous_drainage_analysis_flat_resolution_real_terrain_float_step
    rigorous_drainage_analysis_flat_resolution_real_terrain_gradient
    PROPERTIES LABELS "cli;flow")
lsdtt_set_tests_properties(
    smoke_drainage_analysis_flat_resolution_real_terrain_inline_gradient
    PROPERTIES LABELS "cli;flow")

# The same command on whole fetched Copernicus GLO-30 fixtures: the flow
# benchmark's alpine/lakeland small and medium (fetched with
# lsdtt-fetch-raster when not cached), plus the benchmark suite's
# alpine/lakeland 500 and 1000 rungs when they are cached. inline_gradient,
# barnes2014b_float_step and barnes2014b_gradient must each exit 0, and the
# barnes runs must end with no cell without descent. Stdlib-only script.
# NEEDS NETWORK; skips (77) when a fixture cannot be fetched.
if(LSDTT_ENABLE_PYTHON_TESTS AND LSDTT_ENABLE_NETWORK_TESTS)
    add_test(NAME smoke_flat_resolution_real_terrain_routing
             COMMAND ${Python3_EXECUTABLE} -I
                     "${CMAKE_CURRENT_SOURCE_DIR}/scripts/check_flat_resolution_real_terrain_routing.py"
                     --bin-dir "$<TARGET_FILE_DIR:lsdtt-drainage-analysis>"
                     --work-dir "${CMAKE_CURRENT_BINARY_DIR}/flat_resolution_real_terrain_routing")
    lsdtt_mark_test_runtime(python smoke_flat_resolution_real_terrain_routing)
elseif(NOT LSDTT_ENABLE_NETWORK_TESTS)
    lsdtt_mark_test_unregistered(network smoke_flat_resolution_real_terrain_routing)
else()
    lsdtt_mark_test_unregistered(python smoke_flat_resolution_real_terrain_routing)
endif()
lsdtt_set_tests_properties(smoke_flat_resolution_real_terrain_routing PROPERTIES
    LABELS "flow" SKIP_RETURN_CODE 77 TIMEOUT 600)

# One run with an option per remaining wired CLI, the two options alternating
# across CLIs. Each asserts the "Flat resolution:" block (printed only on the
# barnes path, so the test fails if the option silently falls back to the
# default) and, where the valley's flat reaches an output cheaply, a value the
# default fill does not give (stated in each param file). The extra inputs are
# flat_valley_outlets.csv (the floor cell next to the spill, row 5, col 1),
# flat_valley_source.csv (the head of the valley axis, row 5, col 22) and
# flat_valley_mask.asc (the floor plus the spill cell).
add_lsdtt_cli_smoke_test(
    smoke_basin_extraction_flat_resolution_gradient
    lsdtt-basin-extraction
    basin_extraction_flat_resolution_gradient.param.in
    bevg_junction_basin_outlines.fgb
    --stdout-contains "Flat resolution: Barnes et al. (2014b), https://doi.org/10.1016/j.cageo.2013.01.009, method barnes2014b_gradient"
    --stdout-contains "Flat regions: 1, flat cells: 198, NO_FLOW cells: 195"
    --stdout-contains "Junctions     : 105"
)

# raster_is_filled: true in a second CLI family (chi-analysis carries the
# first): the resolve-only path resolves the valley's flat and routing
# succeeds; inline_gradient leaves the flat for the router, which refuses it.
add_lsdtt_cli_smoke_test(
    smoke_basin_extraction_flat_resolution_filled_float_step
    lsdtt-basin-extraction
    basin_extraction_flat_resolution_filled_float_step.param.in
    bevf_junction_basin_outlines.fgb
    --stdout-contains "Raster is marked as already filled. Skipping fill step."
    --stdout-contains "Flat resolution: Barnes et al. (2014b), https://doi.org/10.1016/j.cageo.2013.01.009, method barnes2014b_float_step"
    --stdout-contains "Flat regions: 1, flat cells: 198, NO_FLOW cells: 195"
    --stdout-contains "Max raise: 0.000366211 m"
    --stdout-contains "Junctions     : 105"
)

add_lsdtt_cli_smoke_test(
    failure_basin_extraction_flat_resolution_filled_inline
    lsdtt-basin-extraction
    basin_extraction_flat_resolution_filled_inline.param.in
    unused_expected_output
    --expect-failure
    --stderr-contains "could not assign a receiver"
)

add_lsdtt_cli_smoke_test(
    smoke_basin_stats_flat_resolution_float_step
    lsdtt-basin-stats
    basin_stats_flat_resolution_float_step.param.in
    bsvf_basin_stats.fgb
    # The default table format: a FlatGeobuf point layer, one point per basin.
    --expect-vector-feature-min bsvf_basin_stats.fgb:1
    --stdout-contains "Flat resolution: Barnes et al. (2014b), https://doi.org/10.1016/j.cageo.2013.01.009, method barnes2014b_float_step"
    --stdout-contains "Flat regions: 1, flat cells: 198, NO_FLOW cells: 195"
    --stdout-contains "Basin pixels  : 207"
)

add_lsdtt_cli_smoke_test(
    failure_basin_stats_bad_output_format
    lsdtt-basin-stats
    basin_stats_bad_output_format.param.in
    unused_expected_output
    --expect-failure
    --stderr-contains "basin_stats_write_vector_format must be flatgeobuf, gpkg or csv"
)

# fill_before_routing: false (the two CLIs that have it). The NOTE is the
# assertion; the router's refusal shows the option really was not applied.
add_lsdtt_cli_smoke_test(
    failure_basin_stats_flat_resolution_no_fill_note
    lsdtt-basin-stats
    basin_stats_flat_resolution_no_fill.param.in
    unused_expected_output
    --expect-failure
    --stdout-contains "NOTE: fill_flat_resolution_choice (barnes2014b_gradient) is ignored because fill_before_routing is false."
    --stderr-contains "could not assign a receiver"
    --stderr-not-contains "needs fill_min_slope > 0"
)

add_lsdtt_cli_smoke_test(
    smoke_basin_surface_metrics_flat_resolution_gradient
    lsdtt-basin-surface-metrics
    basin_surface_metrics_flat_resolution_gradient.param.in
    bsmvg_slope_area_stats.csv
    --expect bsmvg_outline_curvature_flips.csv
    --stdout-contains "Flat resolution: Barnes et al. (2014b), https://doi.org/10.1016/j.cageo.2013.01.009, method barnes2014b_gradient"
    --stdout-contains "Min effective gradient: 0.000141907 m/m"
    --expect-file-contains "bsmvg_slope_area_stats.csv:121,45,"
)

add_lsdtt_cli_smoke_test(
    smoke_channel_network_analysis_flat_resolution_float_step
    lsdtt-channel-network-analysis
    channel_network_analysis_flat_resolution_float_step.param.in
    cnavf_tokunaga_summary_fit.csv
    --stdout-contains "Flat resolution: Barnes et al. (2014b), https://doi.org/10.1016/j.cageo.2013.01.009, method barnes2014b_float_step"
    --stdout-contains "Flat regions: 1, flat cells: 198, NO_FLOW cells: 195"
    --stdout-contains "Junctions     : 105"
    --stdout-contains "Strahler channels           : 80"
)

# The section geometry does not depend on the choice on this fixture (see the
# param file), so this is the block and the outputs only.
add_lsdtt_cli_smoke_test(
    smoke_cross_sections_flat_resolution_gradient
    lsdtt-cross-sections
    cross_sections_flat_resolution_gradient.param.in
    xsvg_cross_sections.fgb
    --expect xsvg_cross_section_sample_points.fgb
    --stdout-contains "Flat resolution: Barnes et al. (2014b), https://doi.org/10.1016/j.cageo.2013.01.009, method barnes2014b_gradient"
    --stdout-contains "Min effective gradient: 0.000141907 m/m"
    --stdout-contains "Cross sections written from 23 channel centres."
)

# lsdtt-bankfull-width's low-slope WARNING has one sentence per choice; one
# test per branch, each asserting the whole WARNING. The barnes runs also show
# the routing change: the second section's drainage area is 20700 m2 on the
# resolved floor against 2300 m2 under the default fill.
add_lsdtt_cli_smoke_test(
    smoke_bankfull_width_flat_resolution_inline_warning
    lsdtt-bankfull-width
    bankfull_width_flat_resolution_inline.param.in
    bfvi_channel_profile.csv
    --stderr-contains "WARNING: 19 section(s) have fitted slope below 10 x fill_min_slope (0.001). Pit filling imposes fill_min_slope on flats, so on these low-gradient reaches the slope (and Manning discharge, which scales as sqrt(slope)) is materially set by the fill. Reduce fill_min_slope or supply a pre-filled DEM."
)

add_lsdtt_cli_smoke_test(
    smoke_bankfull_width_flat_resolution_gradient_warning
    lsdtt-bankfull-width
    bankfull_width_flat_resolution_gradient.param.in
    bfvg_channel_profile.csv
    --stdout-contains "Flat resolution: Barnes et al. (2014b), https://doi.org/10.1016/j.cageo.2013.01.009, method barnes2014b_gradient"
    --stderr-contains "WARNING: 19 section(s) have fitted slope below 10 x fill_min_slope (0.001). The barnes2014b_gradient flat resolver gives flats at least fill_min_slope (more on quantised-up flats, less on capped flats\; see the flat-resolution report above), so on these low-gradient reaches the slope (and Manning discharge, which scales as sqrt(slope)) is materially set by the flat resolution. Reduce fill_min_slope."
    --expect-file-contains "bfvg_channel_profile.csv:1,210.0000,20700.0000,"
)

add_lsdtt_cli_smoke_test(
    smoke_bankfull_width_flat_resolution_float_step_warning
    lsdtt-bankfull-width
    bankfull_width_flat_resolution_float_step.param.in
    bfvf_channel_profile.csv
    --stdout-contains "Flat resolution: Barnes et al. (2014b), https://doi.org/10.1016/j.cageo.2013.01.009, method barnes2014b_float_step"
    --stderr-contains "WARNING: 19 section(s) have fitted slope below 10 x fill_min_slope (0.001). Under barnes2014b_float_step, flat slopes come from float32 steps of order ULP(z) / cell size and fill_min_slope does not apply, so on these low-gradient reaches the slope (and Manning discharge, which scales as sqrt(slope)) is set by elevation precision, not the terrain."
    --expect-file-contains "bfvf_channel_profile.csv:1,210.0000,20700.0000,"
)

add_lsdtt_cli_smoke_test(
    smoke_bankfull_width_flat_resolution_filled_float_step_warning
    lsdtt-bankfull-width
    bankfull_width_flat_resolution_filled_float_step.param.in
    bfvff_channel_profile.csv
    --stdout-contains "Raster marked as already filled. Skipping."
    --stdout-contains "Flat resolution: Barnes et al. (2014b), https://doi.org/10.1016/j.cageo.2013.01.009, method barnes2014b_float_step"
    --stderr-contains "WARNING: 19 section(s) have fitted slope below 10 x fill_min_slope (0.001). Under barnes2014b_float_step, flat slopes come from float32 steps of order ULP(z) / cell size and fill_min_slope does not apply"
)

add_lsdtt_cli_smoke_test(
    smoke_distributed_hydrology_flat_resolution_hype_gradient
    lsdtt-distributed-hydrology
    distributed_hydrology_flat_resolution_hype_gradient.param.in
    hypevg_hype_discharge.csv
    --read-path "${LSDTT_HYPE_CLI_FIXTURE_DIR}"
    --stdout-contains "Column model: HYPE"
    --stdout-contains "Flat resolution: Barnes et al. (2014b), https://doi.org/10.1016/j.cageo.2013.01.009, method barnes2014b_gradient"
    --stdout-contains "Flat regions: 0, flat cells: 0, NO_FLOW cells: 0"
    --expect-file-contains "hypevg_hype_discharge.csv:30,2001-01-30,"
)

add_lsdtt_cli_smoke_test(
    smoke_distributed_hydrology_flat_resolution_topmodel_float_step
    lsdtt-distributed-hydrology
    distributed_hydrology_flat_resolution_topmodel_float_step.param.in
    tmvf_twi.tif
    --expect tmvf_discharge.csv
    --stdout-contains "Flat resolution: Barnes et al. (2014b), https://doi.org/10.1016/j.cageo.2013.01.009, method barnes2014b_float_step"
    --stdout-contains "Flat regions: 1, flat cells: 198, NO_FLOW cells: 195"
    --expect-raster-stat tmvf_twi.tif:mean:11.50:11.57
)

add_lsdtt_cli_smoke_test(
    smoke_fluvial_surfaces_flat_resolution_float_step
    lsdtt-fluvial-surfaces
    fluvial_surfaces_flat_resolution_float_step.param.in
    fsvf_surface_class.csv
    --expect fsvf_hand.tif
    --stdout-contains "Flat resolution: Barnes et al. (2014b), https://doi.org/10.1016/j.cageo.2013.01.009, method barnes2014b_float_step"
    --stdout-contains "Flat regions: 1, flat cells: 198, NO_FLOW cells: 195"
    --stdout-contains "  3 channel source(s) in scope"
)

add_lsdtt_cli_smoke_test(
    failure_fluvial_surfaces_flat_resolution_no_fill_note
    lsdtt-fluvial-surfaces
    fluvial_surfaces_flat_resolution_no_fill.param.in
    unused_expected_output
    --expect-failure
    --stdout-contains "NOTE: fill_flat_resolution_choice (barnes2014b_gradient) is ignored because fill_before_routing is false."
    --stderr-contains "could not assign a receiver"
    --stderr-not-contains "needs fill_min_slope_for_routing > 0"
)

add_lsdtt_cli_smoke_test(
    smoke_ridgeline_flat_resolution_gradient
    lsdtt-ridgeline-analysis
    ridgeline_flat_resolution_gradient.param.in
    rlvg_basin_outlines.csv
    --stdout-contains "Flat resolution: Barnes et al. (2014b), https://doi.org/10.1016/j.cageo.2013.01.009, method barnes2014b_gradient"
    --stdout-contains "Min effective gradient: 0.000141907 m/m"
    --stdout-contains "Filtering initial selection (55 basins)..."
)

add_lsdtt_cli_smoke_test(
    smoke_valley_centreline_flat_resolution_float_step
    lsdtt-valley-centreline
    valley_centreline_flat_resolution_float_step.param.in
    vcvf_valley_centreline.fgb
    --expect vcvf_valley_centreline_nodes.csv
    --stdout-contains "Flat resolution: Barnes et al. (2014b), https://doi.org/10.1016/j.cageo.2013.01.009, method barnes2014b_float_step"
    --stdout-contains "Flat regions: 1, flat cells: 198, NO_FLOW cells: 195"
    --stdout-contains "centreline 21 nodes, 2 px upstream of the head excluded from the trough"
)

lsdtt_set_tests_properties(
    smoke_basin_extraction_flat_resolution_gradient
    smoke_basin_extraction_flat_resolution_filled_float_step
    smoke_basin_stats_flat_resolution_float_step
    smoke_basin_surface_metrics_flat_resolution_gradient
    smoke_channel_network_analysis_flat_resolution_float_step
    smoke_cross_sections_flat_resolution_gradient
    smoke_bankfull_width_flat_resolution_inline_warning
    smoke_bankfull_width_flat_resolution_gradient_warning
    smoke_bankfull_width_flat_resolution_float_step_warning
    smoke_bankfull_width_flat_resolution_filled_float_step_warning
    smoke_distributed_hydrology_flat_resolution_hype_gradient
    smoke_distributed_hydrology_flat_resolution_topmodel_float_step
    smoke_fluvial_surfaces_flat_resolution_float_step
    smoke_ridgeline_flat_resolution_gradient
    smoke_valley_centreline_flat_resolution_float_step
    PROPERTIES LABELS "cli;flow")
lsdtt_set_tests_properties(
    failure_basin_extraction_flat_resolution_filled_inline
    failure_basin_stats_flat_resolution_no_fill_note
    failure_fluvial_surfaces_flat_resolution_no_fill_note
    PROPERTIES LABELS "cli;flow")
lsdtt_set_tests_properties(smoke_distributed_hydrology_flat_resolution_hype_gradient
    PROPERTIES FIXTURES_REQUIRED hype_cli_fixture)

# lsdtt-surface-metrics.  Uses tiny_dem.asc (UTM 30N) rather than tiny30_dem.tif,
# which is EPSG:4326 and is rejected by the CLI's projected-CRS check - that
# rejection is itself covered by failure_surface_metrics_geographic_crs below.
add_lsdtt_cli_smoke_test(
    smoke_surface_metrics
    lsdtt-surface-metrics
    surface_metrics_smoke.param.in
    tiny_hillshade.tif
    --expect tiny_local_relief.tif
    --expect tiny_direct_polyfit_slope.tif
    --expect tiny_direct_polyfit_aspect_degrees.tif
    --expect tiny_direct_polyfit_mean_curvature.tif
    --expect tiny_direct_polyfit_maximum_curvature.tif
    --expect tiny_direct_polyfit_minimum_curvature.tif
    --expect tiny_direct_shape_class.tif
    --expect tiny_direct_shape_gaussian_curvature.tif
    --expect tiny_direct_shape_mean_curvature.tif
)

# Both settled defaults are deliberately omitted: stdout pins the Frangi value
# propagated into the request, while help-topic tests below pin declarations.
add_lsdtt_cli_smoke_test(
    smoke_surface_metrics_frangi_default
    lsdtt-surface-metrics
    surface_metrics_frangi_default.param.in
    frangi_default_frangi.tif
    --stdout-contains "c         : 1"
)

add_lsdtt_cli_smoke_test(
    smoke_surface_metrics_frangi_c15_negative_control
    lsdtt-surface-metrics
    surface_metrics_frangi_c15.param.in
    frangi_c15_frangi.tif
    --stdout-contains "c         : 15"
)

# Valley spacing proper: no test ran it before 2026-10-07. See the parameter
# file for why (consolidation D5 changed its chamfer diagonal weight).
add_lsdtt_cli_smoke_test(
    smoke_spectral_analysis_valley_spacing
    lsdtt-spectral-analysis
    spectral_analysis_valley_spacing.param.in
    vsp_valley_spacing.csv
    --read-path "${LSDTT_STREAM_BURNING_FIXTURE_DIR}/out"
    --expect vsp_valley_spacing_spectrum.csv
    --stdout-contains "--- Valley spacing ---"
    --stdout-contains "Hillslope length (p90)  : 30 m"
)
lsdtt_set_tests_properties(smoke_spectral_analysis_valley_spacing PROPERTIES
    FIXTURES_REQUIRED stream_burning_fixture)

# 2D significant frequencies with bearings (the legacy _significant_freq.csv)
# and background_seed reproducibility. Known-answer sinusoids pin the bearing
# convention, including a non-square grid where dfx != dfy; the seed mode pins
# bit-identical backgrounds across thread counts; --edges pins the Nyquist-edge
# pairing; --align pins the bin matching of the radial CSV's beta and R^2
# columns. --threshold and --align need no FFTW.
add_executable(lsdtt-test-spectral-significant-2d cpp/test_spectral_significant_2d.cpp)
target_link_libraries(lsdtt-test-spectral-significant-2d PRIVATE lsdtt3_cpu)
set(_sss2d_modes "threshold" "align")
if(LSDTT_HAS_FFTW)
    list(APPEND _sss2d_modes "bearings" "seed" "edges")
endif()
foreach(_sss2d_mode IN LISTS _sss2d_modes)
    add_test(NAME unit_spectral_analysis_significant_2d_${_sss2d_mode}
             COMMAND lsdtt-test-spectral-significant-2d --${_sss2d_mode})
    lsdtt_set_tests_properties(unit_spectral_analysis_significant_2d_${_sss2d_mode}
        PROPERTIES LABELS "surface_metrics")
endforeach()

add_lsdtt_cli_smoke_test(
    smoke_spectral_analysis_significant_2d
    lsdtt-spectral-analysis
    spectral_analysis_significant_2d.param.in
    sig2d_significant_frequencies_2d.csv
    --read-path "${LSDTT_STREAM_BURNING_FIXTURE_DIR}/out"
    --expect-file-contains "sig2d_significant_frequencies_2d.csv:bearing_rad,bearing_deg,frequency_cycles_per_m,wavelength_m,normalised_psd"
    --expect-file-contains "sig2d_significant_frequencies_2d.csv:# normalised_psd_threshold,9.21034"
    --expect-file-contains "sig2d_significant_frequencies_2d.csv:# background_seed,17"
    # The fixture DEM has nothing significant at the defaults, so both CSVs
    # must say so in words rather than leave an empty table.
    --expect-file-contains "sig2d_significant_frequencies_2d.csv:# n_significant_frequencies,0"
    --expect-file-contains "sig2d_significant_frequencies_2d.csv:# No significant frequencies were found in this raster"
    --expect-file-contains "sig2d_significant_frequencies.csv:# n_significant_frequencies,0"
    --expect-file-contains "sig2d_significant_frequencies.csv:# No significant frequencies were found in this raster"
)
lsdtt_set_tests_properties(smoke_spectral_analysis_significant_2d PROPERTIES
    FIXTURES_REQUIRED stream_burning_fixture)

# The same DEM at confidence 0.5 (threshold 0.69): cells must be written and the
# no-significant-frequencies line must not appear.
add_lsdtt_cli_smoke_test(
    smoke_spectral_analysis_significant_2d_found
    lsdtt-spectral-analysis
    spectral_analysis_significant_2d_found.param.in
    sig2dfound_significant_frequencies_2d.csv
    --read-path "${LSDTT_STREAM_BURNING_FIXTURE_DIR}/out"
    --expect-file-contains "sig2dfound_significant_frequencies_2d.csv:# normalised_psd_threshold,0.693147"
    --expect-file-not-contains "sig2dfound_significant_frequencies_2d.csv:# n_significant_frequencies,0"
    --expect-file-not-contains "sig2dfound_significant_frequencies_2d.csv:No significant frequencies"
)
lsdtt_set_tests_properties(smoke_spectral_analysis_significant_2d_found PROPERTIES
    FIXTURES_REQUIRED stream_burning_fixture)

add_lsdtt_cli_smoke_test(
    regression_spectral_analysis_valley_spacing_legacy_alias
    lsdtt-spectral-analysis
    spectral_valley_spacing_legacy.param.in
    spectral_legacy_psd_2d.tif
    --expect-file-contains "spectral_legacy_ingested_params.param:valley_spacing_area_exponent: 0.37"
    --stderr-contains "'valley_spacing_m' is deprecated"
)

foreach(default_help_spec IN ITEMS
        "lsdtt-surface-metrics|n_threads|default : 0"
        "lsdtt-surface-metrics|frangi_c|default : 1.0"
        "lsdtt-synthetic-surfaces|n_rows|default : 128"
        "lsdtt-synthetic-surfaces|n_cols|default : 128"
        "lsdtt-spectral-analysis|n_threads|default : 0"
        "lsdtt-spectral-analysis|valley_spacing_area_exponent|default : 0.45"
        "lsdtt-valley-centreline|write_output_format|default : cog"
        "lsdtt-feature-finder|big_channel_threshold_contributing_pixels|default : 5000"
        "lsdtt-mountain-front|basin_minimum_contributing_pixels|default : 1000"
        "lsdtt-lem|incision_threshold_specific_area|A/dx channel threshold in metres"
        "lsdtt-lem|channel_mask_threshold_specific_area|A/dx threshold in metres"
        "lsdtt-drainage-analysis|connectivity_area_threshold_m2|total map area in m\\^2")
    string(REPLACE "|" ";" default_help_fields "${default_help_spec}")
    list(GET default_help_fields 0 default_help_exe)
    list(GET default_help_fields 1 default_help_topic)
    list(GET default_help_fields 2 default_help_regex)
    string(REPLACE "lsdtt-" "" default_help_name "${default_help_exe}")
    string(REPLACE "-" "_" default_help_name "${default_help_name}")
    add_test(NAME regression_${default_help_name}_${default_help_topic}_help
             COMMAND ${default_help_exe} -h ${default_help_topic})
    lsdtt_set_tests_properties(
        regression_${default_help_name}_${default_help_topic}_help PROPERTIES
        LABELS "cli;parameters"
        PASS_REGULAR_EXPRESSION "${default_help_regex}")
endforeach()

if(LSDTT_ENABLE_PYTHON_TESTS)
    add_test(NAME unit_synthetic_surfaces_help_defaults_exact
        COMMAND ${Python3_EXECUTABLE}
            "${CMAKE_CURRENT_SOURCE_DIR}/scripts/check_synthetic_surfaces_help_defaults.py"
            --exe "$<TARGET_FILE:lsdtt-synthetic-surfaces>")
    lsdtt_mark_test_runtime(python unit_synthetic_surfaces_help_defaults_exact)
else()
    lsdtt_mark_test_unregistered(python unit_synthetic_surfaces_help_defaults_exact)
endif()
lsdtt_set_tests_properties(unit_synthetic_surfaces_help_defaults_exact PROPERTIES
    LABELS "cli;parameters")

if(LSDTT_ENABLE_PYTHON_TESTS)
    add_test(NAME unit_synthetic_cli_banners
        COMMAND ${Python3_EXECUTABLE}
            "${CMAKE_CURRENT_SOURCE_DIR}/scripts/check_synthetic_cli_banners.py"
            --surfaces "$<TARGET_FILE:lsdtt-synthetic-surfaces>"
            --networks "$<TARGET_FILE:lsdtt-synthetic-networks>"
            --lem "$<TARGET_FILE:lsdtt-lem>")
    lsdtt_mark_test_runtime(python unit_synthetic_cli_banners)
else()
    lsdtt_mark_test_unregistered(python unit_synthetic_cli_banners)
endif()
lsdtt_set_tests_properties(unit_synthetic_cli_banners PROPERTIES
    LABELS "cli;synthetic_networks;lem")

add_lsdtt_cli_smoke_test(
    smoke_surface_metrics_madk2
    lsdtt-surface-metrics
    surface_metrics_madk2_smoke.param.in
    tiny_madk2_or.tif
    --expect tiny_madk2_dmc.tif
    --expect tiny_madk2_ra.tif
)

add_lsdtt_cli_smoke_test(
    failure_surface_metrics_invalid_metric
    lsdtt-surface-metrics
    surface_metrics_invalid_metric.param.in
    unused_expected_output
    --expect-failure
    --stderr-contains "Unknown polyfit metric: not_a_metric"
)

add_lsdtt_cli_smoke_test(
    failure_surface_metrics_geographic_crs
    lsdtt-surface-metrics
    surface_metrics_geographic_crs.param.in
    unused_expected_output
    --expect-failure
    --stderr-contains "geographic coordinate system"
)

add_lsdtt_cli_smoke_test(
    failure_fetch_raster_bounds_conflict
    lsdtt-fetch-raster
    fetch_raster_bounds_conflict.param.in
    unused_expected_output
    --expect-failure
    --stderr-contains "FATAL: give either bounds_fname or west/south/east/north"
)

add_lsdtt_cli_smoke_test(
    failure_fetch_raster_bounds_three_points
    lsdtt-fetch-raster
    fetch_raster_bounds_three_points.param.in
    unused_expected_output
    --expect-failure
    --stderr-contains "exactly two"
)

if(LSDTT_ENABLE_NETWORK_TESTS)
    add_lsdtt_cli_smoke_test(
        network_fetch_raster_bounds_cli
        lsdtt-fetch-raster
        fetch_raster_bounds_network.param.in
        frompoly_DEM.tif
        --expect-raster-size frompoly_DEM.tif:128x128
    )
    lsdtt_set_tests_properties(network_fetch_raster_bounds_cli
        PROPERTIES LABELS "network;cli;fetch")

    add_lsdtt_cli_smoke_test(
        network_fetch_raster_cli
        lsdtt-fetch-raster
        fetch_raster_network.param.in
        ortler_test_DEM.tif
        --expect-raster-size ortler_test_DEM.tif:128x128
    )
    lsdtt_set_tests_properties(network_fetch_raster_cli PROPERTIES LABELS "network;cli;fetch")

    # A Green River reach fetched straight into its oblique frame: 17 880 m x
    # 2 040 m at 30 m (596 x 68), all valid, with the frame file recording the
    # line and the UTM 12N rotate-back target on the (0, 0) lattice.
    add_lsdtt_cli_smoke_test(
        network_fetch_raster_linear_feature_cli
        lsdtt-fetch-raster
        fetch_raster_feature_network.param.in
        ffgreen_oblique_DEM.tif
        --expect ffgreen_oblique_DEM_rotation_frame.txt
        --expect-raster-size ffgreen_oblique_DEM.tif:596x68
        --expect-valid-percent-100 ffgreen_oblique_DEM.tif
        --expect-file-contains "ffgreen_oblique_DEM_rotation_frame.txt:frame_source: line"
        --expect-file-contains "ffgreen_oblique_DEM_rotation_frame.txt:frame_x_min: -8940.0"
        --expect-file-contains "ffgreen_oblique_DEM_rotation_frame.txt:target_crs_epsg: 32612"
        --expect-file-contains "ffgreen_oblique_DEM_rotation_frame.txt:rotated_dem_fname: 'ffgreen_oblique_DEM.tif'"
        --stdout-contains "oblique frame (rotate back to EPSG:32612)"
    )
    lsdtt_set_tests_properties(network_fetch_raster_linear_feature_cli
        PROPERTIES LABELS "network;cli;fetch")

    add_lsdtt_cli_smoke_test(
        network_fetch_raster_outlets_cli
        lsdtt-fetch-raster
        fetch_raster_outlets_network.param.in
        outlets_net_group_0001_DEM.tif
        --expect outlets_net_group_0002_DEM.tif
        --expect outlets_net_dem_groups.csv
        --expect outlets_net_outlet_bounds_status.csv
        --expect-file-contains "outlets_net_dem_groups.csv:group_id"
        --expect-file-contains "outlets_net_dem_groups.csv:group_0002_DEM.tif"
    )
    lsdtt_set_tests_properties(network_fetch_raster_outlets_cli PROPERTIES LABELS "network;cli;fetch")
endif()

# Flow-routing performance (docs/developer_notes/testing_and_benchmarking.md).
#
# Reported, never asserted: wall-clock thresholds in CI are flaky on shared
# runners. OFF by default because the sweep takes minutes; the value is a
# record comparable across commits, each row carrying a git SHA.
#
# Timings are IN-PROCESS. The lsdtt_seconds column of benchmark_results.csv is
# whole-process CLI wall time and is NOT comparable with RichDEM's in-routine
# numbers -- that is precisely what this driver exists to fix.
option(LSDTT_ENABLE_PERF_TESTS "Enable flow-routing performance benchmarks" OFF)

add_executable(lsdtt-test-flow-performance cpp/test_flow_performance.cpp)
target_link_libraries(lsdtt-test-flow-performance PRIVATE lsdtt3_cpu)

if(LSDTT_ENABLE_PERF_TESTS)
    add_test(NAME performance_flow_routines
             COMMAND lsdtt-test-flow-performance
                     ${CMAKE_CURRENT_BINARY_DIR}/perf_output --quick)
    lsdtt_set_tests_properties(performance_flow_routines
        PROPERTIES LABELS "flow" TIMEOUT 1800)

    # Implicit hillslope solver: per-solve linear-algebra cost and grid scaling
    # (128^2 / 256^2 / 512^2) for each Eigen backend, from the SAME driver as
    # unit_hillslope_implicit via --performance.
    #
    # REPORTED, NEVER ASSERTED -- it can only fail by not running. Kept out of the
    # default suite because the 512^2 LDLT solve alone is seconds of wall time and
    # the timings are worthless (and flaky) on a loaded shared runner.
    add_test(NAME performance_hillslope_implicit
             COMMAND lsdtt-test-hillslope-implicit
                     ${CMAKE_CURRENT_BINARY_DIR}/perf_output --performance)
    lsdtt_set_tests_properties(performance_hillslope_implicit
        PROPERTIES LABELS "hillslope" SKIP_RETURN_CODE 77 TIMEOUT 1800)

    # Surface-metrics timing sweep: every comparison, plus the
    # marginal-cost-per-extra-metric and window-size (O(n)) studies.
    # Reported, never asserted -- no --check, so a slow shared runner cannot
    # turn it red.  Needs pixi + topotoolbox + a fetched fixture; skips (77)
    # when any of the three is missing.
    if(LSDTT_ENABLE_BENCHMARKS AND LSDTT_PIXI_EXECUTABLE)
        add_test(
            NAME performance_surface_metrics
            COMMAND ${LSDTT_BENCHMARK_LAUNCHER}
                    "${CMAKE_CURRENT_SOURCE_DIR}/scripts/run_surface_benchmark.py"
                    --exe "$<TARGET_FILE:lsdtt-surface-metrics>"
                    --fixture alpine_medium
                    --comparison all
                    --repeats 5
                    --work-dir "${CMAKE_CURRENT_BINARY_DIR}/surface_benchmark"
        )
        lsdtt_mark_test_runtime(pixi performance_surface_metrics)
        lsdtt_set_tests_properties(performance_surface_metrics
            PROPERTIES LABELS "surface_metrics" SKIP_RETURN_CODE 77
                       TIMEOUT 3600)
    endif()
    # Per-pixel LSDn production cost. Four cases: Stone vs
    # LSDn at age 0; LSDn age 0 vs 20 ka vs 1 Ma (1 Ma pulls in the long record);
    # the time-averaged surface vs the direct per-pixel loop, both driven
    # explicitly, with the measured crossover pixel count; and the PRE-Phase-8
    # stone2000sp_paleomag baseline. Emits lsdn_performance.csv.
    #
    # REPORTED, NEVER ASSERTED -- it can only fail by not running. That the
    # surface and the direct path AGREE is asserted in L10 (unit_cosmo_physics),
    # which is where a correctness claim belongs.
    add_test(NAME perf_cosmo_lsdn_production
             COMMAND lsdtt-test-cosmo-physics
                     ${CMAKE_CURRENT_BINARY_DIR}/perf_output --performance)
    lsdtt_set_tests_properties(perf_cosmo_lsdn_production
        PROPERTIES LABELS "cosmo" TIMEOUT 1800)

    # HYPE model timing: set-up, per-step
    # and per-unit-step cost, phase split, peak RSS on synthetic in-memory
    # workloads. REPORTED, NEVER ASSERTED -- fails only on a water-balance or
    # non-finite-flow error. --quick runs the small workload only; run the
    # binary directly, without --quick, for the three-size baseline.
    add_test(NAME performance_hype_model
             COMMAND lsdtt-test-hype-performance
                     ${CMAKE_CURRENT_BINARY_DIR}/perf_output --quick)
    lsdtt_set_tests_properties(performance_hype_model
        PROPERTIES LABELS "distributed_hydrology" TIMEOUT 1800)

    file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/perf_output)
endif()

# Credential handling (docs/developer_notes/testing_and_benchmarking.md).
#
# lsdtt_secret is the only place in lsdtt3 that touches secrets, so it gets its
# own driver. Offline and deterministic; always in the default suite.
add_executable(lsdtt-test-secret cpp/test_secret.cpp)
target_link_libraries(lsdtt-test-secret PRIVATE lsdtt3_core)
target_include_directories(lsdtt-test-secret PRIVATE ${CMAKE_SOURCE_DIR}/src/core)

add_test(NAME unit_secret
         COMMAND lsdtt-test-secret ${CMAKE_CURRENT_BINARY_DIR}/secret_test_output)
lsdtt_set_tests_properties(unit_secret PROPERTIES LABELS "fetch")
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/secret_test_output)

# The test launcher's isolation from an active Pixi/conda environment
# (scripts/lsdtt_test_env.py). Runs an existing CLI's --version, which reports
# the GDAL/PROJ loaded at run time and where their data came from. Outside
# Pixi it checks only that the run-time GDAL is the compiled one; inside it
# (`pixi run -e test ctest -R unit_test_env_isolation`) it also checks the CLI
# does not read the environment's GDAL/PROJ data, with a control showing it
# would without the launcher.
if(LSDTT_ENABLE_PYTHON_TESTS)
    add_test(NAME unit_test_env_isolation
             COMMAND ${Python3_EXECUTABLE}
                     "${CMAKE_CURRENT_SOURCE_DIR}/scripts/check_test_env_isolation.py"
                     --exe "$<TARGET_FILE:lsdtt-channel-extraction>"
                     --expect-gdal "${GDAL_VERSION}")
    lsdtt_mark_test_runtime(python unit_test_env_isolation)
else()
    lsdtt_mark_test_unregistered(python unit_test_env_isolation)
endif()
lsdtt_set_tests_properties(unit_test_env_isolation PROPERTIES LABELS "core")

# CONFIGURE CONTRACTS: each disabled or missing test tier must say so at
# configure time (a warning, or an error naming what to install), never leave a
# silently partial suite. Each case configures the whole project afresh in a
# scratch directory (cmake/configure_contract.cmake), about 20 s apiece, hence
# the configure_contract label. The two that need a Python interpreter are
# Python-tier tests.
set(_lsdtt_contract_args
    -DSOURCE_DIR=${CMAKE_SOURCE_DIR}
    "-DGENERATOR=${CMAKE_GENERATOR}"
    -DCXX_COMPILER=${CMAKE_CXX_COMPILER})
foreach(_case native_default geo_without_python)
    add_test(NAME configure_contract_${_case}
             COMMAND ${CMAKE_COMMAND} -DCASE=${_case} ${_lsdtt_contract_args}
                     -DWORK_DIR=${CMAKE_CURRENT_BINARY_DIR}/configure_contract/${_case}
                     -P "${CMAKE_CURRENT_SOURCE_DIR}/cmake/configure_contract.cmake")
    lsdtt_set_tests_properties(configure_contract_${_case} PROPERTIES
        LABELS "configure_contract")
endforeach()
foreach(_case geo_missing_modules gdal_programs_missing)
    if(LSDTT_ENABLE_PYTHON_TESTS)
        add_test(NAME configure_contract_${_case}
                 COMMAND ${CMAKE_COMMAND} -DCASE=${_case} ${_lsdtt_contract_args}
                         -DPYTHON=${Python3_EXECUTABLE}
                         -DWORK_DIR=${CMAKE_CURRENT_BINARY_DIR}/configure_contract/${_case}
                         -P "${CMAKE_CURRENT_SOURCE_DIR}/cmake/configure_contract.cmake")
        lsdtt_mark_test_runtime(python configure_contract_${_case})
    else()
        lsdtt_mark_test_unregistered(python configure_contract_${_case})
    endif()
    lsdtt_set_tests_properties(configure_contract_${_case} PROPERTIES
        LABELS "configure_contract")
endforeach()

# -----------------------------------------------------------------------------
# lsdtt-fetch-lithology CLI. The offline test reads the committed synthetic tile
# (tests/data/macrostrat_tiles) with the API lookup off; the network test is
# opt-in like the fetch_raster ones.
# -----------------------------------------------------------------------------
add_lsdtt_cli_smoke_test(
    smoke_fetch_lithology_offline
    lsdtt-fetch-lithology
    fetch_lithology_offline.param.in
    lithoff_macrostrat_units.fgb
    --expect lithoff_lithology_vocabulary.csv
    --expect lithoff_crosswalk_template.csv
    --expect lithoff_ATTRIBUTION.txt
    --expect-vector-feature-min lithoff_macrostrat_units.fgb:3
    --expect-file-contains "lithoff_lithology_vocabulary.csv:\"granite\",major"
    --expect-file-contains "lithoff_lithology_vocabulary.csv:\"sandstone and shale\",text"
    --expect-file-not-contains "lithoff_lithology_vocabulary.csv:limestone"
    --expect-file-contains "lithoff_ATTRIBUTION.txt:source_id 7"
    --stdout-contains "(EPSG:32632)"
    --stdout-contains "lithology_field_name: lith_primary"
)

# -----------------------------------------------------------------------------
# lsdtt-fetch-river-lines CLI, offline: the committed Open Rivers extract read
# as a local file with the Open Rivers field names (24 links, nothing clipped).
# -----------------------------------------------------------------------------
add_lsdtt_cli_smoke_test(
    smoke_fetch_river_lines_local
    lsdtt-fetch-river-lines
    fetch_river_lines_local.param.in
    rivloc_river_lines.fgb
    --expect rivloc_river_nodes.fgb
    --expect rivloc_river_topology.csv
    --expect rivloc_ATTRIBUTION.txt
    --expect-vector-feature-min rivloc_river_lines.fgb:24
    --stdout-contains "(EPSG:32630)"
    --stdout-contains "Network: 24 lines"
    --stdout-contains "stream_vector_path: rivloc_river_lines.fgb"
)

add_lsdtt_cli_smoke_test(
    failure_fetch_river_lines_no_source
    lsdtt-fetch-river-lines
    fetch_river_lines_no_source.param.in
    unused_expected_output
    --expect-failure
    --stderr-contains "FATAL: set river_line_source"
)

add_lsdtt_cli_smoke_test(
    failure_fetch_lithology_no_bbox
    lsdtt-fetch-lithology
    fetch_lithology_no_bbox.param.in
    unused_expected_output
    --expect-failure
    --stderr-contains "FATAL: give a bounding box"
)

add_lsdtt_cli_smoke_test(
    failure_fetch_lithology_too_many_tiles
    lsdtt-fetch-lithology
    fetch_lithology_too_many_tiles.param.in
    unused_expected_output
    --expect-failure
    --stderr-contains "more than max_tiles"
)

add_lsdtt_cli_smoke_test(
    failure_fetch_lithology_web_mercator
    lsdtt-fetch-lithology
    fetch_lithology_web_mercator.param.in
    unused_expected_output
    --expect-failure
    --stderr-contains "FATAL: output_epsg: Web Mercator"
)

add_test(NAME smoke_fetch_lithology_help
         COMMAND lsdtt-fetch-lithology -h)
lsdtt_set_tests_properties(smoke_fetch_lithology_help PROPERTIES
    LABELS "cli;fetch"
    PASS_REGULAR_EXPRESSION "macrostrat_tile_zoom")

if(LSDTT_ENABLE_NETWORK_TESTS)
    add_lsdtt_cli_smoke_test(
        network_fetch_lithology_cli
        lsdtt-fetch-lithology
        fetch_lithology_network.param.in
        lithlive_macrostrat_units.fgb
        --expect lithlive_lithology_vocabulary.csv
        --expect-vector-feature-min lithlive_macrostrat_units.fgb:1
        --expect-file-contains "lithlive_ATTRIBUTION.txt:doi:10.1029/2018GC007467"
        --stdout-contains "(EPSG:32630)"
    )
    lsdtt_set_tests_properties(network_fetch_lithology_cli PROPERTIES LABELS "network;cli;fetch")
endif()

# -----------------------------------------------------------------------------
# lsdtt-fetch-erosion-rates CLI (OCTOPUS). The offline tests read the committed
# WV responses (tests/data/octopus_erosion_rates); the network test is opt-in
# like the other fetch ones and asserts minima only, because OCTOPUS grows.
# -----------------------------------------------------------------------------
add_lsdtt_cli_smoke_test(
    smoke_fetch_erosion_rates_offline
    lsdtt-fetch-erosion-rates
    fetch_erosion_rates_offline.param.in
    wvoff_erosion_rates.csv
    --expect wvoff_octopus_basins.fgb
    --expect wvoff_octopus_dem_bounds.fgb
    --expect wvoff_ATTRIBUTION.txt
    --expect-vector-feature-min wvoff_octopus_basins.fgb:28
    --expect-file-contains "wvoff_erosion_rates.csv:sample_id,latitude,longitude,erosion_rate,erosion_rate_uncertainty,erosion_rate_units,erosion_rate_uncertainty_kind,erosion_rate_internal_uncertainty,reported_area_km2,rate_scheme,rate_field,nuclide,obsid2,study_id,author,pub_year,ref_doi,mean_elevation_m,outlet_in_bbox"
    --expect-file-contains "wvoff_erosion_rates.csv:S186WTS046,38.88961897,-79.40200545,12.88,2.62,mm_per_kyr,total_1sigma,"
    --expect-file-not-contains "wvoff_erosion_rates.csv:S186WTS001"
    --expect-file-not-contains "wvoff_erosion_rates.csv:-9999"
    --expect-file-contains "wvoff_ATTRIBUTION.txt:10.5194/gchron-7-113-2025"
    --expect-file-contains "wvoff_ATTRIBUTION.txt:S186, Portenga, 2019, 10.1130/b31840.1, 17"
    --stdout-contains "Read 28 OCTOPUS basins (23 with an outlet in the box)"
    --stdout-contains "samples written: 23"
    --stdout-contains "(EPSG:32617)"
    --stdout-contains "erosion_rate_units: mm_per_kyr"
)

add_lsdtt_cli_smoke_test(
    smoke_fetch_erosion_rates_offline_lsdn
    lsdtt-fetch-erosion-rates
    fetch_erosion_rates_offline_lsdn.param.in
    wvofflsdn_erosion_rates.csv
    --expect-file-contains "wvofflsdn_erosion_rates.csv:S186WTS046,38.88961897,-79.40200545,12.8,0.83,mm_per_kyr,external_1sigma,0.29,20.14,lsdn,EBEL,be10,"
    --stdout-contains "scheme lsdn (EBEL, mm/kyr)"
)

add_lsdtt_cli_smoke_test(
    smoke_fetch_erosion_rates_offline_empty
    lsdtt-fetch-erosion-rates
    fetch_erosion_rates_offline_empty.param.in
    wvempty_erosion_rates.csv
    --expect wvempty_ATTRIBUTION.txt
    --stderr-contains "no OCTOPUS basins intersect the box"
    --stdout-contains "(header only)"
)

add_lsdtt_cli_smoke_test(
    smoke_fetch_erosion_rates_raster_bounds
    lsdtt-fetch-erosion-rates
    fetch_erosion_rates_raster_bounds.param.in
    wvrast_erosion_rates.csv
    --stdout-contains "(raster, EPSG:32632)"
    --stdout-contains "samples written: 0"
)

add_lsdtt_cli_smoke_test(
    smoke_fetch_erosion_rates_vector_bounds
    lsdtt-fetch-erosion-rates
    fetch_erosion_rates_vector_bounds.param.in
    wvvec_erosion_rates.csv
    --stdout-contains "(polygon, EPSG:4326)"
)

add_lsdtt_cli_smoke_test(
    failure_fetch_erosion_rates_no_bbox
    lsdtt-fetch-erosion-rates
    fetch_erosion_rates_no_bbox.param.in
    unused_expected_output
    --expect-failure
    --stderr-contains "FATAL: give a bounding box"
)

add_lsdtt_cli_smoke_test(
    failure_fetch_erosion_rates_too_many
    lsdtt-fetch-erosion-rates
    fetch_erosion_rates_too_many.param.in
    unused_expected_output
    --expect-failure
    --stderr-contains "more than max_features"
)

add_lsdtt_cli_smoke_test(
    failure_fetch_erosion_rates_bad_scheme
    lsdtt-fetch-erosion-rates
    fetch_erosion_rates_bad_scheme.param.in
    unused_expected_output
    --expect-failure
    --stderr-contains "st, lsdn, lsdn_no_topo_shielding, published"
)

add_lsdtt_cli_smoke_test(
    failure_fetch_erosion_rates_web_mercator
    lsdtt-fetch-erosion-rates
    fetch_erosion_rates_web_mercator.param.in
    unused_expected_output
    --expect-failure
    --stderr-contains "FATAL: output_epsg: Web Mercator"
)

add_lsdtt_cli_smoke_test(
    failure_fetch_erosion_rates_exception
    lsdtt-fetch-erosion-rates
    fetch_erosion_rates_exception.param.in
    unused_expected_output
    --expect-failure
    --stderr-contains "OCTOPUS WFS error:"
)

add_test(NAME smoke_fetch_erosion_rates_help
         COMMAND lsdtt-fetch-erosion-rates -h)
lsdtt_set_tests_properties(smoke_fetch_erosion_rates_help PROPERTIES
    LABELS "cli;fetch"
    PASS_REGULAR_EXPRESSION "octopus_rate_scheme")

if(LSDTT_ENABLE_NETWORK_TESTS)
    add_lsdtt_cli_smoke_test(
        network_fetch_erosion_rates_cli
        lsdtt-fetch-erosion-rates
        fetch_erosion_rates_network.param.in
        wvlive_erosion_rates.csv
        --expect wvlive_octopus_basins.fgb
        --expect-file-contains "wvlive_erosion_rates.csv:S186WTS046"
        --expect-vector-feature-min wvlive_octopus_basins.fgb:5
        --expect-file-contains "wvlive_ATTRIBUTION.txt:10.5194/gchron-7-113-2025"
        --stdout-contains "with an outlet in the box)"
    )
    lsdtt_set_tests_properties(network_fetch_erosion_rates_cli PROPERTIES LABELS "network;cli;fetch")
endif()

# Native DEM fetching (docs/developer_notes/testing_and_benchmarking.md).
#
# The offline half needs no network and runs in the default suite: it covers UTM
# zone derivation, tile naming, option validation, and -- most importantly --
# that an unobtainable source degrades to "skipped" rather than failing.
add_executable(lsdtt-test-raster-fetch cpp/test_raster_fetch.cpp)
target_link_libraries(lsdtt-test-raster-fetch PRIVATE lsdtt3_core)
target_include_directories(lsdtt-test-raster-fetch PRIVATE ${CMAKE_SOURCE_DIR}/src/core)

add_test(NAME unit_raster_fetch
         COMMAND lsdtt-test-raster-fetch --offline --data-dir=${CMAKE_CURRENT_SOURCE_DIR}/data/bounds)
lsdtt_set_tests_properties(unit_raster_fetch PROPERTIES LABELS "fetch")

add_executable(lsdtt-test-dem-grouping cpp/test_dem_grouping.cpp)
target_link_libraries(lsdtt-test-dem-grouping PRIVATE lsdtt3_core)
target_include_directories(lsdtt-test-dem-grouping PRIVATE ${CMAKE_SOURCE_DIR}/src/core)

add_test(NAME unit_dem_grouping
         COMMAND lsdtt-test-dem-grouping "${CMAKE_CURRENT_BINARY_DIR}/dem_grouping_output")
lsdtt_set_tests_properties(unit_dem_grouping PROPERTIES LABELS "fetch")

# Macrostrat lithology fetch (lsdtt-fetch-lithology), offline: tile maths, the
# API parsers, the lith_primary rule, and a fetch from a committed synthetic
# vector tile (tests/data/macrostrat_tiles). Skips (77) without GDAL's MVT driver.
add_executable(lsdtt-test-macrostrat-fetch cpp/test_macrostrat_fetch.cpp)
target_link_libraries(lsdtt-test-macrostrat-fetch PRIVATE lsdtt3_core)
target_include_directories(lsdtt-test-macrostrat-fetch PRIVATE ${CMAKE_SOURCE_DIR}/src/core)
add_test(NAME unit_fetch_lithology
         COMMAND lsdtt-test-macrostrat-fetch "${CMAKE_CURRENT_SOURCE_DIR}/data/macrostrat_tiles"
                 "${CMAKE_CURRENT_BINARY_DIR}/fetch_lithology_work")
lsdtt_set_tests_properties(unit_fetch_lithology PROPERTIES LABELS "fetch" SKIP_RETURN_CODE 77)

# OCTOPUS erosion-rate fetch (lsdtt-fetch-erosion-rates), offline: the WFS
# request recipe (latitude-first bbox, hits then paging, never maxFeatures), the
# sentinel predicate, the rate-field map, and reads of the committed WV
# responses (tests/data/octopus_erosion_rates). No HTTP request is made.
add_executable(lsdtt-test-octopus-fetch cpp/test_octopus_fetch.cpp)
target_link_libraries(lsdtt-test-octopus-fetch PRIVATE lsdtt3_core)
target_include_directories(lsdtt-test-octopus-fetch PRIVATE ${CMAKE_SOURCE_DIR}/src/core)
add_test(NAME unit_fetch_erosion_rates
         COMMAND lsdtt-test-octopus-fetch "${CMAKE_CURRENT_SOURCE_DIR}/data/octopus_erosion_rates"
                 "${CMAKE_CURRENT_BINARY_DIR}/fetch_erosion_rates_work")
lsdtt_set_tests_properties(unit_fetch_erosion_rates PROPERTIES LABELS "fetch")

# River-line preparation (lsdtt-fetch-river-lines), in process and GDAL-free:
# hand-built polylines covering T-junctions, near misses, reversal, braids,
# multipart features, clipping, node-id merging, conflicts and fragments.
add_executable(lsdtt-test-river-line-prepare cpp/test_river_line_prepare.cpp)
target_link_libraries(lsdtt-test-river-line-prepare PRIVATE lsdtt3_cpu)
target_include_directories(lsdtt-test-river-line-prepare PRIVATE ${CMAKE_SOURCE_DIR}/src/core ${CMAKE_SOURCE_DIR}/src/cpu)
add_test(NAME unit_fetch_river_lines_prepare COMMAND lsdtt-test-river-line-prepare)
lsdtt_set_tests_properties(unit_fetch_river_lines_prepare PROPERTIES LABELS "fetch")

# River-line fetch adapters, offline: the committed Open Rivers extract
# (tests/data/river_lines) zipped into a cache at run time, checked against a
# reference built from its node ids; local_file and reversed-copy relations.
add_executable(lsdtt-test-river-line-fetch cpp/test_river_line_fetch.cpp)
target_link_libraries(lsdtt-test-river-line-fetch PRIVATE lsdtt3_cpu)
target_include_directories(lsdtt-test-river-line-fetch PRIVATE ${CMAKE_SOURCE_DIR}/src/core ${CMAKE_SOURCE_DIR}/src/cpu)
add_test(NAME unit_fetch_river_lines_adapters
         COMMAND lsdtt-test-river-line-fetch "${CMAKE_CURRENT_SOURCE_DIR}/data/river_lines"
                 "${CMAKE_CURRENT_BINARY_DIR}/fetch_river_lines_work")
lsdtt_set_tests_properties(unit_fetch_river_lines_adapters PROPERTIES LABELS "fetch")

add_executable(lsdtt-test-coastal-mask cpp/test_coastal_mask.cpp)
target_link_libraries(lsdtt-test-coastal-mask PRIVATE lsdtt3_core)
target_include_directories(lsdtt-test-coastal-mask PRIVATE ${CMAKE_SOURCE_DIR}/src/core)

add_test(NAME unit_coastal_mask
         COMMAND lsdtt-test-coastal-mask "${CMAKE_CURRENT_BINARY_DIR}/coastal_mask_output")
lsdtt_set_tests_properties(unit_coastal_mask PROPERTIES LABELS "fetch")

# The network half actually downloads. Off by default (LSDTT_ENABLE_NETWORK_TESTS,
# declared at the top of this file); skips (77) rather than failing when the
# source is unreachable.
if(LSDTT_ENABLE_NETWORK_TESTS)
    add_test(NAME network_raster_fetch
             COMMAND lsdtt-test-raster-fetch --network
                     ${CMAKE_CURRENT_BINARY_DIR}/fetch_test_output)
    lsdtt_set_tests_properties(network_raster_fetch
        PROPERTIES LABELS "network;fetch" SKIP_RETURN_CODE 77)
    file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/fetch_test_output)
endif()

# Arguments: output dir, Stone sentinel, LSDn sentinel (a SEPARATE file with its
# own gated regeneration flag), and the tests/data fixture dir (the
# pre-LSDn crn_parameters.nc used by L5).
# Step-change rasters. Checks the analytic algebra (ksn, celerity, knickpoint
# position, the A_0 dependence) and then the rasters themselves: time-since-change
# equals the step age at the outlet and falls to zero at the knickpoint, erosion
# rates split old/new across it, the chi profile is piecewise linear, and n < 1
# produces a rarefaction fan. Emits step_change_profiles.csv.
add_test(NAME unit_step_change
         COMMAND lsdtt-test-step-change
                 ${CMAKE_CURRENT_BINARY_DIR}/step_change_test_output)
lsdtt_set_tests_properties(unit_step_change PROPERTIES LABELS "step_change")
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/step_change_test_output)

# Landslide size distribution T1-T5b. Emits landslide_truncation.csv, which
# tabulates how much of the distribution the A_max cutoff removes by number and
# by volume.
add_test(NAME unit_landslide_size
         COMMAND lsdtt-test-landslide-size
                 ${CMAKE_CURRENT_BINARY_DIR}/landslide_test_output)
lsdtt_set_tests_properties(unit_landslide_size PROPERTIES LABELS "landslide")
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/landslide_test_output)

# Niemi ensemble N1-N6. The CONTROL (N1) is the one that matters most: with no
# landslides the inversion must return the imposed rate, or nothing else here
# means anything. Emits niemi_basins.csv, niemi_area_bins.csv, niemi_ensemble.csv.
# Its OWN output directory: sharing landslide_test_output with the other
# landslide tests lets them race on the CSVs under `ctest -j`.
add_test(NAME unit_niemi_ensemble
         COMMAND lsdtt-test-niemi-ensemble
                 ${CMAKE_CURRENT_BINARY_DIR}/niemi_test_output)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/niemi_test_output)
lsdtt_set_tests_properties(unit_niemi_ensemble PROPERTIES
    LABELS "cosmo;landslide" TIMEOUT 1800)

# Quartz fertility Q1-Q6 and the ice mask I1-I3. Pins the case that makes
# fertility a correctness issue rather than a refinement: a quartz-free unit
# eroding ten times faster than its neighbour must contribute nothing.
add_test(NAME unit_quartz_fertility COMMAND lsdtt-test-quartz-fertility)
lsdtt_set_tests_properties(unit_quartz_fertility PROPERTIES LABELS "cosmo")

add_test(NAME unit_raster_transform COMMAND lsdtt-test-raster-transform)
lsdtt_set_tests_properties(unit_raster_transform PROPERTIES LABELS "raster_algebra;cosmo")

# Landslide scar coupling and densities. Emits landslide_scar_recovery.csv,
# which shows how a scar's concentration recovers with age at three depths and
# what carrying it at the wrong material density does to the answer.
add_test(NAME unit_landslide_step_change
         COMMAND lsdtt-test-landslide-step-change
                 ${CMAKE_CURRENT_BINARY_DIR}/landslide_test_output)
lsdtt_set_tests_properties(unit_landslide_step_change PROPERTIES
    LABELS "landslide;step_change")

# Landslide placement T6-T10d. Emits landslide_host_fill.csv, which records how
# completely landslides of each size fill their host basin.
add_test(NAME unit_landslide_placement
         COMMAND lsdtt-test-landslide-placement
                 ${CMAKE_CURRENT_BINARY_DIR}/landslide_test_output)
lsdtt_set_tests_properties(unit_landslide_placement PROPERTIES LABELS "landslide")

# Landslide event T11-T13b. Emits landslide_budget_vs_amax.csv.
add_test(NAME unit_landslide_event
         COMMAND lsdtt-test-landslide-event
                 ${CMAKE_CURRENT_BINARY_DIR}/landslide_test_output)
lsdtt_set_tests_properties(unit_landslide_event PROPERTIES LABELS "landslide")

# Landslide CRN response T16. Emits landslide_crn_response.csv.
add_test(NAME unit_landslide_crn
         COMMAND lsdtt-test-landslide-crn
                 ${CMAKE_CURRENT_BINARY_DIR}/landslide_test_output)
lsdtt_set_tests_properties(unit_landslide_crn PROPERTIES LABELS "landslide;cosmo")

# Landslide series: Phase 3 multi-event + Phase 5 T14/T15. Emits
# landslide_niemi_closure.csv and landslide_yanites_regimes.csv.
add_test(NAME unit_landslide_series
         COMMAND lsdtt-test-landslide-series
                 ${CMAKE_CURRENT_BINARY_DIR}/landslide_test_output)
lsdtt_set_tests_properties(unit_landslide_series PROPERTIES LABELS "landslide" TIMEOUT 900)

# The isotope-ratio output policy. One division per pixel, and every failure
# mode is silent: an inf written into a raster a basin average then consumes, a
# nodata numerator divided and emitted as a number, a mask disagreeing with its
# own data. R8 pins the 2026-08-05 promotion from a lambda in
# lsdtt-cosmo-prediction as BITWISE faithful on the input shape it was written
# for, so the refactor is proved at unit level and not only through the CLI
# oracle.
# The metre <-> g/cm2 conversions. Four-line functions, and worth pinning
# because the density split of 2026-08-04 exists precisely because one of the
# eight hand-written `x * rho * 0.1` sites chose the wrong density. U7 asserts
# the promoted functions reproduce those hand-written expressions BITWISE --
# `0.1 * rho * m` and `m * rho * 0.1` are not the same double, and the refactor
# oracle compares file hashes.
add_test(NAME unit_cosmo_units COMMAND lsdtt-test-cosmo-units)
lsdtt_set_tests_properties(unit_cosmo_units PROPERTIES LABELS "cosmo")

add_test(NAME unit_cosmo_ratio COMMAND lsdtt-test-cosmo-ratio)
lsdtt_set_tests_properties(unit_cosmo_ratio PROPERTIES LABELS "cosmo")

add_test(NAME unit_cosmo_physics
         COMMAND lsdtt-test-cosmo-physics
                 ${CMAKE_CURRENT_BINARY_DIR}/cosmo_test_output
                 ${CMAKE_CURRENT_SOURCE_DIR}/expected/cosmo/numerics_sentinel.csv
                 ${CMAKE_CURRENT_SOURCE_DIR}/expected/cosmo/numerics_sentinel_lsdn.csv
                 ${CMAKE_CURRENT_SOURCE_DIR}/data)
lsdtt_set_tests_properties(unit_cosmo_physics PROPERTIES LABELS "cosmo")
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/cosmo_test_output)

# -----------------------------------------------------------------------------
# Cosmo CLI integration tests (docs/developer_notes/testing_and_benchmarking.md)
#
# Offline: they exercise parameter-file parsing, help rendering, and the
# projected-CRS guard on the three cosmo CLIs. The help tests also drive the
# bundled-crn_parameters.nc auto-discovery path (src/core/lsdtt_exe_path.hpp)
# indirectly; the direct auto-discovery check lives in unit_cosmo_physics.
# -----------------------------------------------------------------------------

# Robustness: a geographic (non-projected) DEM must be a clean FATAL, not a
# crash or a silent wrong-units run. tiny30_dem.tif is EPSG:4326.
add_lsdtt_cli_smoke_test(
    failure_cosmo_tool_geographic_crs
    lsdtt-cosmo-tool
    cosmo_tool_geographic_crs.param.in
    unused_expected_output
    --expect-failure
    --stderr-contains "DEM must be in a projected CRS"
)
lsdtt_set_tests_properties(failure_cosmo_tool_geographic_crs
    PROPERTIES LABELS "cli;cosmo")

# --- LSDn scaling scheme through the CLI (L8) -----------------
#
# tiny_dem.asc is a 5x5 tilted plane in EPSG:32630 (UTM 30N) with its origin at
# easting/northing 0, which places it near the equator at about 7.489 deg W;
# cosmo_tiny_points.csv holds one sample at the WGS84 position of the (4.5, 0.5)
# cell centre, i.e. the low corner the plane drains to. threshold_contributing_
# pixels is dropped to 2 so a 25-pixel DEM has a channel to snap to.
#
# The run must SUCCEED, write _PROD.tif with every pixel valid and positive, and
# report the scheme in its run header. The header assertion is the point: a
# silent fall back to Stone would otherwise pass every other check here.
add_lsdtt_cli_smoke_test(
    smoke_cosmo_tool_lsdn
    lsdtt-cosmo-tool
    cosmo_tool_lsdn_smoke.param.in
    cosmo_lsdn_PROD.tif
    --expect cosmo_lsdn_CN_ERATE.csv
    --expect-raster-min cosmo_lsdn_PROD.tif:0.0001
    --stdout-contains "Scaling scheme: LSDn"
    --stdout-contains "LSDn scaling (Be10)"
)
lsdtt_set_tests_properties(smoke_cosmo_tool_lsdn
    PROPERTIES LABELS "cli;cosmo")

# self_consistent_age: the sample age and the erosion rate are solved together,
# t = Gamma_sp/eps.  The assertions are chosen so the test cannot pass if the
# iteration silently degenerates into the fixed-age path:
#   - "self-consistent iter 2" proves it went round MORE THAN ONCE.  A one-shot
#     loop, or one that converged on the first comparison, would not print it.
#   - the CSV must carry the three diagnostic columns, so the converged age is
#     recorded per sample rather than only narrated on stdout.
# Verified to fail on perturbation: forcing the loop to a single pass drops the
# "iter 2" line, and dropping the CSV columns fails the file-contains checks.
add_lsdtt_cli_smoke_test(
    smoke_cosmo_tool_self_consistent_age
    lsdtt-cosmo-tool
    cosmo_tool_self_consistent_age.param.in
    cosmo_sca_CN_ERATE.csv
    --stdout-contains "self-consistent iter 1"
    --stdout-contains "self-consistent iter 2"
    --stdout-contains "self-consistent: t = Gamma_sp/eps per sample"
    --expect-file-contains "cosmo_sca_CN_ERATE.csv:self_consistent_t_yr"
    --expect-file-contains "cosmo_sca_CN_ERATE.csv:self_consistent_n_iter"
    --expect-file-contains "cosmo_sca_CN_ERATE.csv:self_consistent_converged"
    --stderr-contains "'points_filename' is deprecated"
)
lsdtt_set_tests_properties(smoke_cosmo_tool_self_consistent_age
    PROPERTIES LABELS "cli;cosmo")

# A sample that runs out of iterations must be REPORTED, not quietly accepted.
# max_iter 1 cannot converge from the modern field (the first step moves the rate
# ~7 %), so this pins the warning text and the converged=0 flag in the CSV.  The
# run still SUCCEEDS -- non-convergence is a caveat on a number, not a failure.
add_lsdtt_cli_smoke_test(
    smoke_cosmo_tool_self_consistent_no_converge
    lsdtt-cosmo-tool
    cosmo_tool_self_consistent_no_converge.param.in
    cosmo_sca_nc_CN_ERATE.csv
    --stderr-contains "did not converge in 1 self_consistent_age iterations"
    --stderr-contains "Treat it as unconverged, not as a result"
    --expect-file-contains "cosmo_sca_nc_CN_ERATE.csv:,1,0"
)
lsdtt_set_tests_properties(smoke_cosmo_tool_self_consistent_no_converge
    PROPERTIES LABELS "cli;cosmo")

# A non-positive tolerance would make the loop exit on its first comparison and
# report a "converged" age that is nothing of the kind, so it is FATAL up front.
add_lsdtt_cli_smoke_test(
    failure_cosmo_tool_bad_self_consistent_tol
    lsdtt-cosmo-tool
    cosmo_tool_bad_self_consistent_tol.param.in
    unused_expected_output
    --expect-failure
    --stderr-contains "self_consistent_age_tolerance must be > 0"
)
lsdtt_set_tests_properties(failure_cosmo_tool_bad_self_consistent_tol
    PROPERTIES LABELS "cli;cosmo")

# The forward -> inverse CLI BRIDGE, plus knob
# coverage. lsdtt-cosmo-prediction forward-models a concentration for a known
# erosion rate and writes the points CSV it advertises as ingestible;
# lsdtt-cosmo-tool inverts that CSV; the imposed rate must come back. This is
# the only test that drives the two CLIs against each other, so it is the only
# one that can catch plumbing faults -- a column renamed on one side, a knob
# honoured asymmetrically, a shielding term applied twice -- all of which
# produce plausible numbers and pass every single-CLI test.
#
# Offline and fast (a 5x5 fixture DEM), so it runs in the default suite.
# Verified to fail on perturbation: ignoring the standardisation factor, swapping
# the +/-1 sigma bounds, and dropping topographic shielding on the inverse side
# each turn it red.
if(LSDTT_ENABLE_PYTHON_TESTS)
    add_test(NAME rigorous_cosmo_cli_bridge
             COMMAND ${Python3_EXECUTABLE}
                     "${CMAKE_CURRENT_SOURCE_DIR}/scripts/run_cosmo_bridge_test.py"
                     --pred-exe "$<TARGET_FILE:lsdtt-cosmo-prediction>"
                     --tool-exe "$<TARGET_FILE:lsdtt-cosmo-tool>"
                     --data-dir "${LSDTT_TEST_DATA_DIR}")
    lsdtt_mark_test_runtime(python rigorous_cosmo_cli_bridge)
else()
    lsdtt_mark_test_unregistered(python rigorous_cosmo_cli_bridge)
endif()
lsdtt_set_tests_properties(rigorous_cosmo_cli_bridge PROPERTIES
    LABELS "cli;cosmo" SKIP_RETURN_CODE 77 TIMEOUT 900)

# SA2 rigorous check: in multi-isotope step-change + accumulation mode, sampled
# ratios must be computed as ratios of accumulated basin means (not as means of
# per-pixel ratios).
if(LSDTT_ENABLE_PYTHON_TESTS)
    add_test(NAME rigorous_cosmo_multi_accum_ratio
             COMMAND ${Python3_EXECUTABLE}
                     "${CMAKE_CURRENT_SOURCE_DIR}/scripts/run_cosmo_multi_accum_ratio_check.py"
                     --pred-exe "$<TARGET_FILE:lsdtt-cosmo-prediction>"
                     --data-dir "${LSDTT_TEST_DATA_DIR}")
    lsdtt_mark_test_runtime(python rigorous_cosmo_multi_accum_ratio)
else()
    lsdtt_mark_test_unregistered(python rigorous_cosmo_multi_accum_ratio)
endif()
lsdtt_set_tests_properties(rigorous_cosmo_multi_accum_ratio PROPERTIES
    LABELS "cli;cosmo" SKIP_RETURN_CODE 77 TIMEOUT 900)

# Multi-isotope step change: Be-10, Al-26 and C-14 through ONE scenario, plus
# the two ratio rasters. The isotope RATIOS are what carry a transient signal,
# and a ratio is only meaningful if both members saw an identical scenario --
# hence one particle evolution per pixel, not three separate runs.
#
# Bitwise equality of each nuclide against the single-nuclide kernel is asserted
# in the unit tier (M1); this pins the CLI plumbing: all five rasters appear,
# every pixel valid, and the ratio rasters carry positive, physically sane
# values. Verified to fail on perturbation.
add_lsdtt_cli_smoke_test(
    smoke_cosmo_prediction_multi_isotope
    lsdtt-cosmo-prediction
    cosmo_prediction_multi_isotope.param.in
    mi_Be10_step_crn_conc.tif
    --expect mi_Al26_step_crn_conc.tif
    --expect mi_C14_step_crn_conc.tif
    --expect mi_R2610_step.tif
    --expect mi_R1410_step.tif
    --expect-valid-percent-100 mi_R2610_step.tif
    --expect-valid-percent-100 mi_R1410_step.tif
    --expect-raster-min mi_R2610_step.tif:1.0
    --expect-raster-min mi_R1410_step.tif:0.01
    --stdout-contains "multi-isotope step-change rasters (Be10, Al26, C14)"
    --stdout-contains "Al26/Be10 ratio raster"
    --stdout-contains "C14/Be10 ratio raster"
)
lsdtt_set_tests_properties(smoke_cosmo_prediction_multi_isotope
    PROPERTIES LABELS "cli;cosmo")

# SA2 end-to-end path: one multi-isotope step-change run produces accumulated
# basin means, accumulated ratio rasters and sampled point ratios without
# rerunning for each isotope.
add_lsdtt_cli_smoke_test(
    smoke_cosmo_prediction_multi_isotope_accum_sampling
    lsdtt-cosmo-prediction
    cosmo_prediction_multi_isotope_accum_sampling.param.in
    mi_acc_sampled_crn_conc.csv
    --expect mi_acc_Be10_crn_acc.tif
    --expect mi_acc_Al26_crn_acc.tif
    --expect mi_acc_C14_crn_acc.tif
    --expect mi_acc_R2610_acc.tif
    --expect mi_acc_R1410_acc.tif
    --expect mi_acc_sampled_crn_conc_Be10.csv
    --expect mi_acc_sampled_crn_conc_Al26.csv
    --expect mi_acc_sampled_crn_conc_C14.csv
    --expect-file-contains "mi_acc_sampled_crn_conc.csv:Be10_CONC"
    --expect-file-contains "mi_acc_sampled_crn_conc.csv:R2610_ACC"
    --expect-file-contains "mi_acc_sampled_crn_conc.csv:R1410_ACC"
    --expect-file-contains "mi_acc_sampled_crn_conc_Be10.csv:Be10_CONC"
    --expect-file-contains "mi_acc_sampled_crn_conc_Al26.csv:Al26_CONC"
    --expect-file-contains "mi_acc_sampled_crn_conc_C14.csv:C14_CONC"
    --stdout-contains "Accumulating multi-isotope concentration downstream"
)
lsdtt_set_tests_properties(smoke_cosmo_prediction_multi_isotope_accum_sampling
    PROPERTIES LABELS "cli;cosmo")

# Multi-isotope step change under LSDn scaling: smoke that the prediction CLI
# reaches the LSDn path and writes the same 3 nuclide + 2 ratio rasters.
add_lsdtt_cli_smoke_test(
    smoke_cosmo_prediction_multi_isotope_lsdn
    lsdtt-cosmo-prediction
    cosmo_prediction_multi_isotope_lsdn.param.in
    mi_lsdn_Be10_step_crn_conc.tif
    --expect mi_lsdn_Al26_step_crn_conc.tif
    --expect mi_lsdn_C14_step_crn_conc.tif
    --expect mi_lsdn_R2610_step.tif
    --expect mi_lsdn_R1410_step.tif
    --expect-valid-percent-100 mi_lsdn_R2610_step.tif
    --expect-valid-percent-100 mi_lsdn_R1410_step.tif
    --expect-raster-min mi_lsdn_R2610_step.tif:1.0
    --expect-raster-min mi_lsdn_R1410_step.tif:0.01
    --stdout-contains "Scaling scheme: LSDn"
    --stdout-contains "LSDn scaling (Be10)"
)
lsdtt_set_tests_properties(smoke_cosmo_prediction_multi_isotope_lsdn
    PROPERTIES LABELS "cli;cosmo")

# Multi-isotope ARBITRARY erosion history: the general case of the step-change
# mode. Bitwise equality with the step-change kernel for a two-row history is
# asserted in the unit tier (H1); this pins the CLI plumbing -- the history file
# is read and validated, all five rasters appear, and the ratios are sane.
add_lsdtt_cli_smoke_test(
    smoke_cosmo_prediction_multi_transient
    lsdtt-cosmo-prediction
    cosmo_prediction_multi_transient.param.in
    tr_Be10_step_crn_conc.tif
    --expect tr_Al26_step_crn_conc.tif
    --expect tr_C14_step_crn_conc.tif
    --expect tr_R2610_step.tif
    --expect tr_R1410_step.tif
    --expect-valid-percent-100 tr_R2610_step.tif
    --expect-raster-min tr_R2610_step.tif:1.0
    --stdout-contains "multi-isotope TRANSIENT rasters (Be10, Al26, C14)"
    --stdout-contains "every pixel walks the whole history"
    --stderr-contains "'transient_history_file' is deprecated"
)
lsdtt_set_tests_properties(smoke_cosmo_prediction_multi_transient
    PROPERTIES LABELS "cli;cosmo")

# An unknown scheme is FATAL, and the message names both allowed values.
add_lsdtt_cli_smoke_test(
    failure_cosmo_tool_bad_scaling_scheme
    lsdtt-cosmo-tool
    cosmo_tool_bad_scaling_scheme.param.in
    unused_expected_output
    --expect-failure
    --stderr-contains "unknown scaling_scheme"
    --stderr-contains "stone"
    --stderr-contains "lsdn"
)
lsdtt_set_tests_properties(failure_cosmo_tool_bad_scaling_scheme
    PROPERTIES LABELS "cli;cosmo")

# scaling_scheme: lsdn against an nc with NO LSDn grids must fail hard and name
# the nc it searched -- never silently degrade to Stone. crn_parameters_pre_lsdn.nc
# is the tracked pre-LSDn nc, which is exactly that situation. Library-level
# counterpart: L5(e) in unit_cosmo_physics.
add_lsdtt_cli_smoke_test(
    failure_cosmo_tool_lsdn_without_grids
    lsdtt-cosmo-tool
    cosmo_tool_lsdn_no_grids.param.in
    unused_expected_output
    --expect-failure
    --stderr-contains "requires the LSDn scaling grids"
    --stderr-contains "crn_parameters_pre_lsdn.nc"
)
lsdtt_set_tests_properties(failure_cosmo_tool_lsdn_without_grids
    PROPERTIES LABELS "cli;cosmo")

# Help must render every parameter without a parameter file and without
# touching the DEM or the network. Guards that the parameter declarations for
# each cosmo CLI stay well-formed.
# lsdtt-cosmo-landslides.
add_test(NAME smoke_cosmo_landslides_help COMMAND lsdtt-cosmo-landslides -h)
lsdtt_set_tests_properties(smoke_cosmo_landslides_help PROPERTIES
    LABELS "cli;landslide"
    PASS_REGULAR_EXPRESSION "long_term_erosion_rate_mm_kyr")

# End-to-end on the tiny fixture: an event is placed and the inventory written.
add_lsdtt_cli_smoke_test(
    smoke_cosmo_landslides_single_event
    lsdtt-cosmo-landslides
    cosmo_landslides_single_event.param.in
    lssmoke_LS_INVENTORY.csv
    --stdout-contains "landslide event:"
    --expect-file-contains "lssmoke_LS_INVENTORY.csv:n_pixels"
)

# C10 characterization: these three paths were previously untested. Counts and
# output checksums are frozen from the inline implementation before extraction.
add_lsdtt_cli_smoke_test(
    regression_cosmo_landslides_active_values
    lsdtt-cosmo-landslides
    cosmo_landslides_active_values.param.in
    lsactive_LS_TFAIL.tif
    --stdout-contains "Mask values   : 4 -> 19 host pixels"
    --stdout-contains "Active area   : 18 m2 (18 px)"
    --expect-valid-percent lsactive_LS_TFAIL.tif:0
    --expect-raster-checksum lsactive_LS_TFAIL.tif:65255
)

add_lsdtt_cli_smoke_test(
    regression_cosmo_landslides_ice_mask
    lsdtt-cosmo-landslides
    cosmo_landslides_ice_mask.param.in
    lsice_LS_TFAIL.tif
    --stdout-contains "4 pixel(s) excluded from landsliding"
    --stdout-contains "Active area   : 4 m2 (4 px)"
    --expect-valid-percent lsice_LS_TFAIL.tif:0
    --expect-raster-checksum lsice_LS_TFAIL.tif:65255
)

add_lsdtt_cli_smoke_test(
    regression_cosmo_landslides_active_values_ice
    lsdtt-cosmo-landslides
    cosmo_landslides_active_values_ice.param.in
    lsactiveice_LS_TFAIL.tif
    --stdout-contains "Mask values   : 4 -> 19 host pixels"
    --stdout-contains "4 pixel(s) excluded from landsliding"
    --stdout-contains "Active area   : 2 m2 (2 px)"
    --expect-valid-percent lsactiveice_LS_TFAIL.tif:0
    --expect-raster-checksum lsactiveice_LS_TFAIL.tif:65255
)

# R2c: CRN concentrations in lsdtt-cosmo-landslides can now build a real site
# (per-pixel production + shielding) instead of only the scalar convenience
# path.
add_lsdtt_cli_smoke_test(
    smoke_cosmo_landslides_crn_site_mode
    lsdtt-cosmo-landslides
    cosmo_landslides_crn_site_mode.param.in
    lscrn_LS_CONC.tif
    --stdout-contains "CRN site mode: per-pixel"
)

add_lsdtt_cli_smoke_test(
    smoke_cosmo_landslides_sample_ensemble
    lsdtt-cosmo-landslides
    cosmo_landslides_sample_ensemble.param.in
    lsens_LS_SAMPLE_EFFECT.csv
    --expect-file-contains "lsens_LS_SAMPLE_EFFECT.csv:ratio_to_control_median"
    --expect-file-contains "lsens_LS_SAMPLE_EFFECT.csv:ORDER_FIRST,"
    --expect-file-contains "lsens_LS_SAMPLE_EFFECT.csv:ORDER_SECOND,-7.4887036,4.5097e-06,duplicate_of:ORDER_FIRST"
    --expect-file-contains "lsens_LS_SAMPLE_EFFECT.csv:OUTSIDE,10,10,failed_no_channel"
    --stderr-contains "'sample_points_csv' is deprecated"
)

# A geographic-CRS DEM must exit cleanly with a diagnosis, NOT abort on an
# uncaught throw. The deeper check throws; the CLI catches it.
add_lsdtt_cli_smoke_test(
    failure_cosmo_landslides_geographic_crs
    lsdtt-cosmo-landslides
    cosmo_landslides_geographic.param.in
    lsgeo_LS_INVENTORY.csv
    --expect-failure
    --stderr-contains "projected"
)

# ---------------------------------------------------------------------------
# Cosmogenic parameter validation in lsdtt-cosmo-landslides.
#
# Until 2026-08-05 a bad muon_scheme or nuclide was caught inside the CRN block,
# downgraded to a stderr warning, and the tool exited **0** having written no
# concentration raster. To a scripted workflow that is a success that silently
# produced nothing -- a worse failure mode than a crash, because nothing
# notices. All three cases below were measured on the shipped binary before the
# fix; the first two exited 0.
add_lsdtt_cli_smoke_test(
    failure_cosmo_landslides_bad_muon_scheme
    lsdtt-cosmo-landslides
    cosmo_landslides_bad_muon_scheme.param.in
    lsbadmuon_LS_CONC.tif
    --expect-failure
    --stderr-contains "unknown muon_scheme"
)

add_lsdtt_cli_smoke_test(
    failure_cosmo_landslides_bad_scaling_scheme
    lsdtt-cosmo-landslides
    cosmo_landslides_bad_scaling_scheme.param.in
    lsbadscale_LS_CONC.tif
    --expect-failure
    --stderr-contains "unknown scaling_scheme"
)

# VALIDITY IS NOT CAPABILITY. is_valid_muon_scheme() accepts the CRONUS v3
# emulator schemes, but they are pressure- and nuclide-dependent and this tool
# reaches the column through reset_scaling(), which throws for them. A shared
# validator that only checked the name would accept something that fails later.
# The message must also point at the tool that DOES support them.
add_lsdtt_cli_smoke_test(
    failure_cosmo_landslides_emulator_scheme
    lsdtt-cosmo-landslides
    cosmo_landslides_emulator_scheme.param.in
    lsemul_LS_CONC.tif
    --expect-failure
    --stderr-contains "lsdtt-cosmo-prediction"
)

# The other side of the guard, and the reason it is conditional: muon_scheme is
# INERT with no cosmogenic output requested, so a stale line must not fail the
# run. Without this test the natural "fix" is to validate unconditionally, which
# would break every pure-geomorphology parameter file carrying a cosmo line.
add_lsdtt_cli_smoke_test(
    smoke_cosmo_landslides_geomorph_only_ignores_muon_scheme
    lsdtt-cosmo-landslides
    cosmo_landslides_geomorph_only.param.in
    lsgeomorph_LS_INVENTORY.csv
    --expect lsgeomorph_LS_TFAIL.tif
)

add_test(NAME smoke_cosmo_tool_help COMMAND lsdtt-cosmo-tool -h)
lsdtt_set_tests_properties(smoke_cosmo_tool_help PROPERTIES
    LABELS "cli;cosmo"
    PASS_REGULAR_EXPRESSION "sample_age_yr")

add_test(NAME smoke_divide_migration_help COMMAND lsdtt-divide-migration -h)
lsdtt_set_tests_properties(smoke_divide_migration_help PROPERTIES
    LABELS "cli"
    PASS_REGULAR_EXPRESSION "basin_outlet_fname")

# Generate the step change from the analytical fluvial solution and feed the
# generated time-since-change field straight into the concentration kernel.
# Checks that the derived K is reported (it is back-calculated, so the whole
# configuration hangs off it) and that the tag raster is written.
add_lsdtt_cli_smoke_test(
    smoke_cosmo_prediction_step_change
    lsdtt-cosmo-prediction
    cosmo_prediction_step_change.param.in
    scgen_Be10_step_crn_conc.tif
    --stdout-contains "K (back-calculated)"
    --stdout-contains "wave form"
    --expect scgen_SC_STEPTIME.tif
    --expect scgen_SC_TAG.tif
)

# The generated time-since-change field is consumed in memory. This asserts
# that route agrees with writing it to disk and reading it back in a second
# run -- the workflow the legacy tool forced. See the script header for why the
# comparison is a checksum plus a 1e-6 statistics tolerance rather than bitwise.
if(LSDTT_ENABLE_PYTHON_TESTS)
    add_test(NAME regression_step_change_disk_roundtrip
             COMMAND ${Python3_EXECUTABLE}
                     ${CMAKE_CURRENT_SOURCE_DIR}/scripts/run_step_change_roundtrip.py
                     --exe $<TARGET_FILE:lsdtt-cosmo-prediction>
                     --data-dir ${LSDTT_TEST_DATA_DIR}
                     --work-dir ${CMAKE_CURRENT_BINARY_DIR}/step_change_roundtrip)
    lsdtt_mark_test_runtime(python regression_step_change_disk_roundtrip)
else()
    lsdtt_mark_test_unregistered(python regression_step_change_disk_roundtrip)
endif()
lsdtt_set_tests_properties(regression_step_change_disk_roundtrip PROPERTIES
    LABELS "cli;cosmo;step_change")

# The scar coupling end to end: lsdtt-cosmo-landslides writes the failure-time and
# excavation rasters, lsdtt-cosmo-prediction consumes them. Asserts that
# switching the coupling on leaves unfailed ground BITWISE unchanged -- the
# property that stops a landslide run silently perturbing the whole raster.
if(LSDTT_ENABLE_PYTHON_TESTS)
    add_test(NAME regression_landslide_scar_coupling
             COMMAND ${Python3_EXECUTABLE}
                     ${CMAKE_CURRENT_SOURCE_DIR}/scripts/run_landslide_scar_coupling.py
                     --landslides-exe $<TARGET_FILE:lsdtt-cosmo-landslides>
                     --prediction-exe $<TARGET_FILE:lsdtt-cosmo-prediction>
                     --data-dir ${LSDTT_TEST_DATA_DIR}
                     --work-dir ${CMAKE_CURRENT_BINARY_DIR}/landslide_scar_coupling)
    lsdtt_mark_test_runtime(python regression_landslide_scar_coupling)
else()
    lsdtt_mark_test_unregistered(python regression_landslide_scar_coupling)
endif()
lsdtt_set_tests_properties(regression_landslide_scar_coupling PROPERTIES
    LABELS "cli;cosmo;landslide;step_change")

# SC path smoke/regression: the combined multi-isotope step-change + landslide
# coupling route that was previously blocked by a fatal guard.
if(LSDTT_ENABLE_PYTHON_TESTS)
    add_test(NAME regression_cosmo_multi_landslide_coupling
             COMMAND ${Python3_EXECUTABLE}
                     ${CMAKE_CURRENT_SOURCE_DIR}/scripts/run_cosmo_multi_landslide_coupling.py
                     --landslides-exe $<TARGET_FILE:lsdtt-cosmo-landslides>
                     --prediction-exe $<TARGET_FILE:lsdtt-cosmo-prediction>
                     --data-dir ${LSDTT_TEST_DATA_DIR}
                     --work-dir ${CMAKE_CURRENT_BINARY_DIR}/cosmo_multi_landslide_coupling)
    lsdtt_mark_test_runtime(python regression_cosmo_multi_landslide_coupling)
else()
    lsdtt_mark_test_unregistered(python regression_cosmo_multi_landslide_coupling)
endif()
lsdtt_set_tests_properties(regression_cosmo_multi_landslide_coupling PROPERTIES
    LABELS "cli;cosmo;landslide;step_change")

# ---------------------------------------------------------------------------
# The REFACTOR ORACLE.
#
# tests/data/cosmo_oracle.csv holds a bitwise hash of every raster and CSV the
# four cosmogenic CLIs produce across a matrix of scaling scheme, muon scheme,
# shielding path and prediction mode. It exists because the production/pressure/
# shielding setup is about to be extracted from the CLIs into a shared
# component, and THE ORACLE FOR A REFACTOR IS THE PRE-REFACTOR BINARY -- which
# stops existing the moment the refactor lands.
#
# Every subsequent change is one of two kinds and must declare which:
#   behaviour-PRESERVING  this test must stay green. That is the proof.
#   behaviour-CHANGING    re-freeze deliberately (--freeze) and say in the
#                         commit message which artefacts moved and why.
#
# Do NOT re-freeze to make a red test pass.
#
# Verified to fail: a 1-part-in-1e6 perturbation of the production raster moves
# 16 of the 22 cheap-tier artefacts, including lsdtt-cosmo-tool's inverted
# erosion-rate CSVs -- so the oracle reaches through the whole forward/inverse
# chain, not just the raster it perturbs.
if(LSDTT_ENABLE_PYTHON_TESTS)
    add_test(NAME regression_cosmo_oracle
             COMMAND ${Python3_EXECUTABLE}
                     ${CMAKE_CURRENT_SOURCE_DIR}/scripts/freeze_cosmo_oracle.py
                     --exe lsdtt-cosmo-prediction=$<TARGET_FILE:lsdtt-cosmo-prediction>
                     --exe lsdtt-cosmo-tool=$<TARGET_FILE:lsdtt-cosmo-tool>
                     --exe lsdtt-cosmo-column-tool=$<TARGET_FILE:lsdtt-cosmo-column-tool>
                     --exe lsdtt-cosmo-landslides=$<TARGET_FILE:lsdtt-cosmo-landslides>
                     --read-path ${LSDTT_TEST_DATA_DIR}
                     --oracle ${LSDTT_TEST_DATA_DIR}/cosmo_oracle.csv
                     --work-dir ${CMAKE_CURRENT_BINARY_DIR}/cosmo_oracle
                     --tier cheap)
    lsdtt_mark_test_runtime(python regression_cosmo_oracle)
else()
    lsdtt_mark_test_unregistered(python regression_cosmo_oracle)
endif()
lsdtt_set_tests_properties(regression_cosmo_oracle PROPERTIES
    LABELS "cli;cosmo")

# The LSDn and CRONUS-emulator rows. Separated because they need the LSDn
# scaling grids in crn_parameters.nc and carry the most setup state -- which is
# exactly why they matter to the extraction, and why they are worth the label
# rather than being dropped.
if(LSDTT_ENABLE_PYTHON_TESTS)
    add_test(NAME regression_cosmo_oracle_heavy
             COMMAND ${Python3_EXECUTABLE}
                     ${CMAKE_CURRENT_SOURCE_DIR}/scripts/freeze_cosmo_oracle.py
                     --exe lsdtt-cosmo-prediction=$<TARGET_FILE:lsdtt-cosmo-prediction>
                     --exe lsdtt-cosmo-tool=$<TARGET_FILE:lsdtt-cosmo-tool>
                     --exe lsdtt-cosmo-column-tool=$<TARGET_FILE:lsdtt-cosmo-column-tool>
                     --exe lsdtt-cosmo-landslides=$<TARGET_FILE:lsdtt-cosmo-landslides>
                     --read-path ${LSDTT_TEST_DATA_DIR}
                     --oracle ${LSDTT_TEST_DATA_DIR}/cosmo_oracle.csv
                     --work-dir ${CMAKE_CURRENT_BINARY_DIR}/cosmo_oracle_heavy
                     --tier all)
    lsdtt_mark_test_runtime(python regression_cosmo_oracle_heavy)
else()
    lsdtt_mark_test_unregistered(python regression_cosmo_oracle_heavy)
endif()
lsdtt_set_tests_properties(regression_cosmo_oracle_heavy PROPERTIES
    LABELS "cli;cosmo;heavy")

# Two prediction modes that write the SAME filenames must not both run: the
# transient goes first and the step change overwrote all five of its outputs,
# so the run reported success having produced half its results. Found while
# promoting the ratio helper; guarded rather than renamed, because
# changing shipped output filenames is a separate decision.
add_lsdtt_cli_smoke_test(
    failure_cosmo_prediction_both_multi_modes
    lsdtt-cosmo-prediction
    cosmo_prediction_both_multi_modes.param.in
    bothmulti_Be10_step_crn_conc.tif
    --expect-failure
    --stderr-contains "silently overwrite"
)

add_test(NAME smoke_cosmo_prediction_help COMMAND lsdtt-cosmo-prediction -h)
lsdtt_set_tests_properties(smoke_cosmo_prediction_help PROPERTIES
    LABELS "cli;cosmo"
    PASS_REGULAR_EXPRESSION "sample_age_yr")

add_test(NAME regression_cosmo_prediction_builtin_io_help
         COMMAND lsdtt-cosmo-prediction -h read_path)
lsdtt_set_tests_properties(regression_cosmo_prediction_builtin_io_help PROPERTIES
    LABELS "cli;cosmo;parameters"
    PASS_REGULAR_EXPRESSION "Defaults to the directory holding the parameter file")

# Functional coverage for lsdtt-cosmo-column-tool, which until 2026-07-27 had
# only the -h test below. It is the natural place to compare lsdtt3 against
# POINT calculators (CRONUS v3, CRONUScalc): no DEM, no flow routing, no basin
# averaging, so nothing can mask a scaling or muon discrepancy.
#
# The checks are LIMIT IDENTITIES of the column model, needing no external
# reference: step_change at t=0 must equal steady state at the old rate and at
# t->inf the new rate; a constant-rate full_transient must relax to steady
# state; monotonicity in depth and rate; production_scaling exactly linear in P0
# while total_shielding is not (Vermeesch redistribution); C-14 relaxing faster
# than Be-10. Verified to fail on perturbation.
if(LSDTT_ENABLE_PYTHON_TESTS)
    add_test(NAME rigorous_cosmo_column_tool
             COMMAND ${Python3_EXECUTABLE}
                     "${CMAKE_CURRENT_SOURCE_DIR}/scripts/run_cosmo_column_test.py"
                     --exe "$<TARGET_FILE:lsdtt-cosmo-column-tool>")
    lsdtt_mark_test_runtime(python rigorous_cosmo_column_tool)
else()
    lsdtt_mark_test_unregistered(python rigorous_cosmo_column_tool)
endif()
lsdtt_set_tests_properties(rigorous_cosmo_column_tool PROPERTIES
    LABELS "cli;cosmo" SKIP_RETURN_CODE 77 TIMEOUT 600)

add_test(NAME smoke_cosmo_column_tool_help COMMAND lsdtt-cosmo-column-tool -h)
lsdtt_set_tests_properties(smoke_cosmo_column_tool_help PROPERTIES
    LABELS "cli;cosmo"
    PASS_REGULAR_EXPRESSION "nuclide")

# -----------------------------------------------------------------------------
# Real-landscape cosmogenic validation tier
# (tests/expected/cosmo/published/README.md)
#
# "download -> run -> delete": run_cosmo_validation.py fetches a DEM KEYLESS via
# lsdtt-fetch-raster (Copernicus GLO-30 / Planetary Computer -- no API key),
# clips it, runs lsdtt-cosmo-tool, compares to the OCTOPUS-harvested published
# rates in the tracked sample table, then deletes the temporary working directory
# INCLUDING ON FAILURE (cleanup trap). Without --check a run FREEZES the summary
# (append); the ctest test below uses --check to REGRESS a fresh run against the
# frozen summary WITHOUT modifying the tracked file.
# -----------------------------------------------------------------------------

# Harness self-test: OFFLINE. Exercises the cleanup-trap contract (the trap fires
# and removes the workdir on a simulated failed run), the manifest/sample-table
# schema guards, and the summary-append round-trip -- with a mock workdir, no
# network, no DEM. This is the cheap tracked guard that the harness logic stays
# correct; it runs in the default unit;cosmo suite.
if(LSDTT_ENABLE_PYTHON_TESTS)
    add_test(NAME unit_cosmo_validation_harness
             COMMAND ${Python3_EXECUTABLE}
                     "${CMAKE_CURRENT_SOURCE_DIR}/scripts/run_cosmo_validation.py"
                     --self-test)
    lsdtt_mark_test_runtime(python unit_cosmo_validation_harness)
else()
    lsdtt_mark_test_unregistered(python unit_cosmo_validation_harness)
endif()
lsdtt_set_tests_properties(unit_cosmo_validation_harness PROPERTIES LABELS "cosmo")

# Riversand-comparison harness self-test (Phase 5 carried-forward item §6c(B)).
# OFFLINE scaffold check only: manifest/sample preflight + cleanup-trap contract.
if(LSDTT_ENABLE_PYTHON_TESTS)
    add_test(NAME unit_cosmo_riversand_harness
             COMMAND ${Python3_EXECUTABLE}
                     "${CMAKE_CURRENT_SOURCE_DIR}/scripts/run_cosmo_riversand_comparison.py"
                     --self-test)
    lsdtt_mark_test_runtime(python unit_cosmo_riversand_harness)
else()
    lsdtt_mark_test_unregistered(python unit_cosmo_riversand_harness)
endif()
lsdtt_set_tests_properties(unit_cosmo_riversand_harness PROPERTIES LABELS "cosmo")

# Real download -> run -> delete for the pilot dataset, in --check (regression)
# mode: it reproduces the fresh inverted rates and compares them to the FROZEN
# tracked summary without modifying it. NEEDS NETWORK and time, so every
# cosmo_realworld_* test is registered only with LSDTT_ENABLE_NETWORK_TESTS=ON
# (and the Python tests); run them with `ctest -L cosmo_realworld`. Until
# 2026-09-18 they were registered unconditionally, so a plain `ctest` downloaded
# DEMs for all eight. Skips (exit 77) when the network/DEM source is
# unavailable OR the manifest/sample table is still a template, so it can never
# turn CI red.
if(LSDTT_ENABLE_PYTHON_TESTS AND LSDTT_ENABLE_NETWORK_TESTS)
    add_test(NAME cosmo_realworld_dethier2014
             COMMAND ${Python3_EXECUTABLE}
                     "${CMAKE_CURRENT_SOURCE_DIR}/scripts/run_cosmo_validation.py"
                     dethier2014 --check
                     --exe "$<TARGET_FILE:lsdtt-cosmo-tool>"
                     --fetch-exe "$<TARGET_FILE:lsdtt-fetch-raster>")
    lsdtt_mark_test_runtime(python cosmo_realworld_dethier2014)
elseif(NOT LSDTT_ENABLE_NETWORK_TESTS)
    lsdtt_mark_test_unregistered(network cosmo_realworld_dethier2014)
else()
    lsdtt_mark_test_unregistered(python cosmo_realworld_dethier2014)
endif()
lsdtt_set_tests_properties(cosmo_realworld_dethier2014 PROPERTIES
    LABELS "cosmo_realworld" SKIP_RETURN_CODE 77 TIMEOUT 1800)

# Additional real-landscape datasets (harvested with tests/scripts/harvest_octopus.py),
# chosen to CONTRAST with the moderate-relief Front Range pilot across terrain and
# elevation: clapp2002 = Yuma Wash arid LOW-relief / LOW-elevation (OCTOPUS S014);
# dibiase2010 = San Gabriel STEEP / HIGH-relief (OCTOPUS S055). Each is fully
# manifest-driven (bbox, target UTM, run knobs) and frozen from its own accepted
# networked run. Same --check regression mode, same skip-77 / out-of-CI wiring.
if(LSDTT_ENABLE_PYTHON_TESTS AND LSDTT_ENABLE_NETWORK_TESTS)
    add_test(NAME cosmo_realworld_clapp2002
             COMMAND ${Python3_EXECUTABLE}
                     "${CMAKE_CURRENT_SOURCE_DIR}/scripts/run_cosmo_validation.py"
                     clapp2002 --check
                     --exe "$<TARGET_FILE:lsdtt-cosmo-tool>"
                     --fetch-exe "$<TARGET_FILE:lsdtt-fetch-raster>")
    lsdtt_mark_test_runtime(python cosmo_realworld_clapp2002)
elseif(NOT LSDTT_ENABLE_NETWORK_TESTS)
    lsdtt_mark_test_unregistered(network cosmo_realworld_clapp2002)
else()
    lsdtt_mark_test_unregistered(python cosmo_realworld_clapp2002)
endif()
lsdtt_set_tests_properties(cosmo_realworld_clapp2002 PROPERTIES
    LABELS "cosmo_realworld" SKIP_RETURN_CODE 77 TIMEOUT 1800)

if(LSDTT_ENABLE_PYTHON_TESTS AND LSDTT_ENABLE_NETWORK_TESTS)
    add_test(NAME cosmo_realworld_dibiase2010
             COMMAND ${Python3_EXECUTABLE}
                     "${CMAKE_CURRENT_SOURCE_DIR}/scripts/run_cosmo_validation.py"
                     dibiase2010 --check
                     --exe "$<TARGET_FILE:lsdtt-cosmo-tool>"
                     --fetch-exe "$<TARGET_FILE:lsdtt-fetch-raster>")
    lsdtt_mark_test_runtime(python cosmo_realworld_dibiase2010)
elseif(NOT LSDTT_ENABLE_NETWORK_TESTS)
    lsdtt_mark_test_unregistered(network cosmo_realworld_dibiase2010)
else()
    lsdtt_mark_test_unregistered(python cosmo_realworld_dibiase2010)
endif()
lsdtt_set_tests_properties(cosmo_realworld_dibiase2010 PROPERTIES
    LABELS "cosmo_realworld" SKIP_RETURN_CODE 77 TIMEOUT 1800)

# LSDn tier of the same three datasets: the
# END-TO-END PER-PIXEL LSDn run, (a) per-pixel LSDn basin mean vs (c) OCTOPUS EBEL.
# Each <name>_lsdn manifest REUSES its St counterpart's bbox, run knobs and SAMPLE
# TABLE -- there is exactly one copy of every measured concentration in the repo --
# and changes only scaling_scheme (lsdn) and the published target (EBE_* -> EBEL,
# via published_erate_override reading the tracked octopus_lsdn_point_reference.csv).
# The frozen band is a MONITORING BAND, NOT AN ORACLE: it mixes the kernel residual
# (pinned DEM-free by L12), the per-pixel-vs-lumped difference (pinned by L13, and a
# result rather than an error) and the DEM/delineation difference, and it runs at the
# CLI default sample_age_yr = 0 while EBEL is time-averaged over Gamma/eps. Attribute
# it with L12/L13 (see each summary header); never tune anything to shrink it.
# Same --check regression mode, same skip-77 / out-of-CI wiring as the St tier: these
# compare the fresh lsdtt3 rates against the frozen rows and must NOT rewrite them.
if(LSDTT_ENABLE_PYTHON_TESTS AND LSDTT_ENABLE_NETWORK_TESTS)
    add_test(NAME cosmo_realworld_dethier2014_lsdn
             COMMAND ${Python3_EXECUTABLE}
                     "${CMAKE_CURRENT_SOURCE_DIR}/scripts/run_cosmo_validation.py"
                     dethier2014_lsdn --check
                     --exe "$<TARGET_FILE:lsdtt-cosmo-tool>"
                     --fetch-exe "$<TARGET_FILE:lsdtt-fetch-raster>")
    lsdtt_mark_test_runtime(python cosmo_realworld_dethier2014_lsdn)
elseif(NOT LSDTT_ENABLE_NETWORK_TESTS)
    lsdtt_mark_test_unregistered(network cosmo_realworld_dethier2014_lsdn)
else()
    lsdtt_mark_test_unregistered(python cosmo_realworld_dethier2014_lsdn)
endif()
lsdtt_set_tests_properties(cosmo_realworld_dethier2014_lsdn PROPERTIES
    LABELS "cosmo_realworld" SKIP_RETURN_CODE 77 TIMEOUT 1800)

if(LSDTT_ENABLE_PYTHON_TESTS AND LSDTT_ENABLE_NETWORK_TESTS)
    add_test(NAME cosmo_realworld_clapp2002_lsdn
             COMMAND ${Python3_EXECUTABLE}
                     "${CMAKE_CURRENT_SOURCE_DIR}/scripts/run_cosmo_validation.py"
                     clapp2002_lsdn --check
                     --exe "$<TARGET_FILE:lsdtt-cosmo-tool>"
                     --fetch-exe "$<TARGET_FILE:lsdtt-fetch-raster>")
    lsdtt_mark_test_runtime(python cosmo_realworld_clapp2002_lsdn)
elseif(NOT LSDTT_ENABLE_NETWORK_TESTS)
    lsdtt_mark_test_unregistered(network cosmo_realworld_clapp2002_lsdn)
else()
    lsdtt_mark_test_unregistered(python cosmo_realworld_clapp2002_lsdn)
endif()
lsdtt_set_tests_properties(cosmo_realworld_clapp2002_lsdn PROPERTIES
    LABELS "cosmo_realworld" SKIP_RETURN_CODE 77 TIMEOUT 1800)

if(LSDTT_ENABLE_PYTHON_TESTS AND LSDTT_ENABLE_NETWORK_TESTS)
    add_test(NAME cosmo_realworld_dibiase2010_lsdn
             COMMAND ${Python3_EXECUTABLE}
                     "${CMAKE_CURRENT_SOURCE_DIR}/scripts/run_cosmo_validation.py"
                     dibiase2010_lsdn --check
                     --exe "$<TARGET_FILE:lsdtt-cosmo-tool>"
                     --fetch-exe "$<TARGET_FILE:lsdtt-fetch-raster>")
    lsdtt_mark_test_runtime(python cosmo_realworld_dibiase2010_lsdn)
elseif(NOT LSDTT_ENABLE_NETWORK_TESTS)
    lsdtt_mark_test_unregistered(network cosmo_realworld_dibiase2010_lsdn)
else()
    lsdtt_mark_test_unregistered(python cosmo_realworld_dibiase2010_lsdn)
endif()
lsdtt_set_tests_properties(cosmo_realworld_dibiase2010_lsdn PROPERTIES
    LABELS "cosmo_realworld" SKIP_RETURN_CODE 77 TIMEOUT 1800)

# §6e production-vs-BEPROD isolation (decisive test for the ~+4% Dethier residual).
# At fixed concentration inferred erosion is proportional to basin-mean production;
# this fetches the keyless COP30 DEM, runs lsdtt-cosmo-tool MATCHED to OCTOPUS
# (Braucher / NCEP / paleomag-off) with the reporting-only basin-production CSV, and
# compares lsdtt's basin-mean Stone-2000 production scaling to OCTOPUS BEPROD/BETOTS
# and its basin-mean elevation to ELEV_AVE (harvested live from the OCTOPUS WFS).
# NEEDS NETWORK; registered only with LSDTT_ENABLE_NETWORK_TESTS=ON (label
# cosmo_realworld, `ctest -L cosmo_realworld`). Skips (exit 77) when the network is unavailable.
# Runs in check-only mode: it compares against the tracked frozen artefact and must
# NOT rewrite it. Re-freezing is a deliberate act -- run the script with --freeze by
# hand and review the diff.
if(LSDTT_ENABLE_PYTHON_TESTS AND LSDTT_ENABLE_NETWORK_TESTS)
    add_test(NAME cosmo_realworld_production_isolation
             COMMAND ${Python3_EXECUTABLE}
                     "${CMAKE_CURRENT_SOURCE_DIR}/scripts/run_cosmo_production_isolation.py"
                     --exe "$<TARGET_FILE:lsdtt-cosmo-tool>"
                     --fetch-exe "$<TARGET_FILE:lsdtt-fetch-raster>")
    lsdtt_mark_test_runtime(python cosmo_realworld_production_isolation)
elseif(NOT LSDTT_ENABLE_NETWORK_TESTS)
    lsdtt_mark_test_unregistered(network cosmo_realworld_production_isolation)
else()
    lsdtt_mark_test_unregistered(python cosmo_realworld_production_isolation)
endif()
lsdtt_set_tests_properties(cosmo_realworld_production_isolation PROPERTIES
    LABELS "cosmo_realworld" SKIP_RETURN_CODE 77 TIMEOUT 1800)

# L13 GENERATION STEP -- the PER-PIXEL side of per-pixel-vs-lumped LSDn
# (test L13). Fetches each of the three
# tracked datasets' keyless COP30 DEMs and records lsdtt-cosmo-tool's basin-mean
# production scaling (Stone, and LSDn at three sample ages) over the exact upslope
# node set the inversion uses.
#
# Deliberately a TWO-STEP design: this expensive networked step is frozen into
# tests/expected/cosmo/lsdn_pixel_vs_lumped_perpixel.csv, and the CHEAP comparison
# against the lumped LSDn kernel is then regression-checked on every build inside
# unit_cosmo_physics (default unit;cosmo suite), which is where the assertions and
# the lsdn_pixel_vs_lumped.csv product live.
#
# NEEDS NETWORK; registered only with LSDTT_ENABLE_NETWORK_TESTS=ON (label cosmo_realworld). Skips
# (exit 77) when the DEM source is unreachable. Runs check-only: it compares
# against the tracked frozen artefact and must NOT rewrite it. Re-freezing is a
# deliberate act -- run the script with --freeze by hand and review the diff.
if(LSDTT_ENABLE_PYTHON_TESTS AND LSDTT_ENABLE_NETWORK_TESTS)
    add_test(NAME cosmo_realworld_lsdn_pixel_vs_lumped
             COMMAND ${Python3_EXECUTABLE}
                     "${CMAKE_CURRENT_SOURCE_DIR}/scripts/run_cosmo_lsdn_pixel_vs_lumped.py"
                     all
                     --exe "$<TARGET_FILE:lsdtt-cosmo-tool>"
                     --fetch-exe "$<TARGET_FILE:lsdtt-fetch-raster>")
    lsdtt_mark_test_runtime(python cosmo_realworld_lsdn_pixel_vs_lumped)
elseif(NOT LSDTT_ENABLE_NETWORK_TESTS)
    lsdtt_mark_test_unregistered(network cosmo_realworld_lsdn_pixel_vs_lumped)
else()
    lsdtt_mark_test_unregistered(python cosmo_realworld_lsdn_pixel_vs_lumped)
endif()
lsdtt_set_tests_properties(cosmo_realworld_lsdn_pixel_vs_lumped PROPERTIES
    LABELS "cosmo_realworld" SKIP_RETURN_CODE 77 TIMEOUT 3600)

# =============================================================================
# REDUNDANCY GUARD: no hand-written metre <-> g/cm2 conversion
#
# `x * rho * 0.1` was open-coded at nine sites across the library and two CLIs,
# each independently choosing WHICH density to use. That is not a style
# complaint: the density split of 2026-08-04 exists because one of those choices
# was wrong, and a soil scar carried at bedrock density is 1.77x too deep in
# mass -- a 2.05x error in spallation production at the scar base.
#
# lsdtt_crn_units.hpp now owns the arithmetic, with a REQUIRED density argument.
# This guard stops the pattern coming back one call site at a time, which is how
# it accumulated in the first place. It is a CONFIGURE ERROR, matching the
# "unclassified test" guard below: a lint nobody runs is a lint that does not
# exist.
#
# EXEMPT, deliberately:
#   * src/cpu/lsdtt_crn_column.{hpp,cpp} and lsdtt_crn_particle.hpp -- LOCKED
#     NUMERICS. to_eff(), remove_mass() and the four depth_m recoveries live
#     inside the validated core and must not be touched. They are excluded by
#     NOT being in the glob below; that is deliberate, not an oversight.
#     See feedback_cosmo_numerics_locked.
#   * lsdtt_crn_units.hpp itself, which is where the arithmetic now lives.
# =============================================================================

set(_lsdtt_unit_guard_files "")
file(GLOB_RECURSE _lsdtt_conv_candidates
     "${CMAKE_SOURCE_DIR}/src/cpu/lsdtt_landslide_*.cpp"
     "${CMAKE_SOURCE_DIR}/src/cpu/lsdtt_cosmo_*.cpp"
     "${CMAKE_SOURCE_DIR}/src/cpu/lsdtt_quartz_*.cpp"
     "${CMAKE_SOURCE_DIR}/src/cpu/lsdtt_step_change_*.cpp"
     "${CMAKE_SOURCE_DIR}/src/cpu/lsdtt_niemi_*.cpp"
     "${CMAKE_SOURCE_DIR}/src/cpu/lsdtt_landslide_*.hpp"
     "${CMAKE_SOURCE_DIR}/src/cpu/lsdtt_cosmo_*.hpp"
     "${CMAKE_SOURCE_DIR}/src/cpu/lsdtt_quartz_*.hpp"
     "${CMAKE_SOURCE_DIR}/src/cpu/lsdtt_step_change_*.hpp"
     "${CMAKE_SOURCE_DIR}/src/cpu/lsdtt_niemi_*.hpp"
     "${CMAKE_SOURCE_DIR}/cli/lsdtt-cosmo-*.cpp")

foreach(_f IN LISTS _lsdtt_conv_candidates)
    file(STRINGS "${_f}" _hits REGEX "(density|rho)[a-z_0-9]* *\\* *0\\.1[^0-9]")
    foreach(_h IN LISTS _hits)
        # Skip comments and help/doc strings -- they describe the conversion,
        # which is exactly what we want them to keep doing.
        string(STRIP "${_h}" _stripped)
        if(NOT _stripped MATCHES "^//" AND NOT _stripped MATCHES "^\\*" AND
           NOT _stripped MATCHES "^\"")
            list(APPEND _lsdtt_unit_guard_files "${_f}: ${_stripped}")
        endif()
    endforeach()
endforeach()

if(_lsdtt_unit_guard_files)
    string(REPLACE ";" "\n    " _guard_report "${_lsdtt_unit_guard_files}")
    message(FATAL_ERROR
        "Hand-written metre <-> g/cm2 conversion found:\n    ${_guard_report}\n\n"
        "Use lsdtt_crn_units.hpp (metres_to_g_cm2 / m_yr_to_eff_erate / "
        "g_cm2_to_metres).\nThe density argument is REQUIRED there, which is the "
        "whole point: the bug this\nprevents was never the arithmetic, it was "
        "using 2650 for a soil scar.")
endif()

# =============================================================================
# lsdtt-flat-surface-finder and lsdtt-valley-centreline
#
# Phase 0a of the fluvial-surfaces split. Both CLIs
# shipped with no coverage of any kind, so there was nothing to refactor
# against. These tests exist to pin CURRENT behaviour before the split moves it.
#
# The fixture is a projected 200x200 linear-range DEM: rough mountains, flat
# piedmont, and a channel network shed onto the piedmont. tests/data has no
# projected DEM large enough (tiny_dem.asc is 5x5, tiny30_dem.tif is geographic
# and both CLIs refuse it via check_projected_metres()).
# =============================================================================

# find_risers (lsdtt_riser_finder_cpu), behind lsdtt-feature-finder
# feature_type find_risers: staircase risers recovered at their true 7.4 m
# step, astar gap closing across a road cut out of the search mask, no
# segments on a ramp or a noisy flat.
add_executable(lsdtt-test-riser-finder cpp/test_riser_finder.cpp)
target_link_libraries(lsdtt-test-riser-finder PRIVATE lsdtt3_cpu)
add_test(NAME unit_feature_finder_risers COMMAND lsdtt-test-riser-finder)
lsdtt_set_tests_properties(unit_feature_finder_risers PROPERTIES
    LABELS "surface_metrics")

# compare_surface_classes / agreement bits: the comparison behind running
# several methods in one lsdtt-fluvial-surfaces call, on class maps
# small enough to check by hand.
add_executable(lsdtt-test-surface-method-comparison cpp/test_surface_method_comparison.cpp)
target_link_libraries(lsdtt-test-surface-method-comparison PRIVATE lsdtt3_cpu)
add_test(NAME unit_fluvial_surface_method_comparison COMMAND lsdtt-test-surface-method-comparison)
lsdtt_set_tests_properties(unit_fluvial_surface_method_comparison PROPERTIES
    LABELS "fluvial_surfaces")

# Valley widths for every channel (plan_valley_width_network): the channel
# tree, both centreline choices, the tributary stop and the tributary-crossing
# exclusion, nested, on synthetic valleys of known floor width.
add_executable(lsdtt-test-valley-width-network cpp/test_valley_width_network.cpp)
target_link_libraries(lsdtt-test-valley-width-network PRIVATE lsdtt3_cpu)
add_test(NAME unit_valley_width_network COMMAND lsdtt-test-valley-width-network)
lsdtt_set_tests_properties(unit_valley_width_network PROPERTIES
    LABELS "fluvial_surfaces")

# Relative elevation raster (lsdtt_valley_relative_elevation_cpu): the two
# centreline baselines on a hand-built zigzag channel, and the larger channel
# winning the pixels its valley discs cover.
add_executable(lsdtt-test-valley-relative-elevation cpp/test_valley_relative_elevation.cpp)
target_link_libraries(lsdtt-test-valley-relative-elevation PRIVATE lsdtt3_cpu)
add_test(NAME unit_valley_relative_elevation COMMAND lsdtt-test-valley-relative-elevation)
lsdtt_set_tests_properties(unit_valley_relative_elevation PROPERTIES
    LABELS "fluvial_surfaces")

# Perimeter sources (lsdtt_perimeter_sources_cpu, plan phases P1-P4): the
# trace on rectangles, frames, a rotated square, a neck, components and holes;
# gentle runs, gap closing, the low-spot test and the seam; outlet runs and
# sources on valley, tributary, perched-notch, basin-clipped, two-network
# and coast DEMs.
add_executable(lsdtt-test-perimeter-sources cpp/test_perimeter_sources.cpp)
target_link_libraries(lsdtt-test-perimeter-sources PRIVATE lsdtt3_cpu)
add_test(NAME unit_perimeter_sources COMMAND lsdtt-test-perimeter-sources)
lsdtt_set_tests_properties(unit_perimeter_sources PROPERTIES
    LABELS "fluvial_surfaces")

# CLI pins on the trunk-and-tributary and two-network shapes the unit
# driver writes: 2 sources and 1 source, and boundary conditions other than
# oooo refused.
set(LSDTT_PERIMETER_SOURCES_FIXTURE_DIR
    "${CMAKE_CURRENT_BINARY_DIR}/perimeter_sources_fixture")
add_test(NAME fixture_perimeter_sources_cli
         COMMAND lsdtt-test-perimeter-sources --write-cli-fixture
                 ${LSDTT_PERIMETER_SOURCES_FIXTURE_DIR})
lsdtt_set_tests_properties(fixture_perimeter_sources_cli PROPERTIES
    LABELS "fixture"
    FIXTURES_SETUP perimeter_sources_fixture)
add_lsdtt_cli_smoke_test(
    smoke_channel_extraction_perimeter_sources_tributary
    lsdtt-channel-extraction
    channel_extraction_perimeter_sources_tributary.param.in
    pst_perimeter_sources.fgb
    --read-path "${LSDTT_PERIMETER_SOURCES_FIXTURE_DIR}"
    --expect pst_perimeter_outlets.fgb
    --expect pst_perimeter_crossings.fgb
    --expect pst_perimeter_profile.csv
    --stdout-contains "Sources found                 : 2"
    --stdout-contains "Gentle runs by verdict        : outlet_crossing 1"
)
lsdtt_set_tests_properties(smoke_channel_extraction_perimeter_sources_tributary PROPERTIES
    FIXTURES_REQUIRED perimeter_sources_fixture)
add_lsdtt_cli_smoke_test(
    smoke_channel_extraction_perimeter_sources_two_network
    lsdtt-channel-extraction
    channel_extraction_perimeter_sources_two_network.param.in
    pss_perimeter_sources.fgb
    --read-path "${LSDTT_PERIMETER_SOURCES_FIXTURE_DIR}"
    --stdout-contains "Sources found                 : 1"
    --stdout-contains "outlet rank 1 (selected)"
)
lsdtt_set_tests_properties(smoke_channel_extraction_perimeter_sources_two_network PROPERTIES
    FIXTURES_REQUIRED perimeter_sources_fixture)
add_lsdtt_cli_smoke_test(
    failure_channel_extraction_perimeter_sources_boundary
    lsdtt-channel-extraction
    channel_extraction_perimeter_sources_boundary.param.in
    psb_perimeter_sources.fgb
    --read-path "${LSDTT_PERIMETER_SOURCES_FIXTURE_DIR}"
    --expect-failure
    --stderr-contains "FATAL: run_perimeter_sources requires boundary_conditions: oooo (got bbbb)."
)
lsdtt_set_tests_properties(failure_channel_extraction_perimeter_sources_boundary PROPERTIES
    FIXTURES_REQUIRED perimeter_sources_fixture)

# Training polygons: class labels, best separation, a window fit
# with a mislabelled terrace number, a flagged polygon and a shared level, and
# the reach windows, on populations small enough to check by hand.
add_executable(lsdtt-test-surface-training cpp/test_surface_training.cpp)
target_link_libraries(lsdtt-test-surface-training PRIVATE lsdtt3_cpu)
add_test(NAME unit_fluvial_surface_training COMMAND lsdtt-test-surface-training)
lsdtt_set_tests_properties(unit_fluvial_surface_training PROPERTIES
    LABELS "fluvial_surfaces")

add_executable(lsdtt-test-fluvial-surfaces cpp/test_fluvial_surfaces.cpp)
target_link_libraries(lsdtt-test-fluvial-surfaces PRIVATE lsdtt3_cpu)

add_test(NAME unit_fluvial_surfaces_fixture
         COMMAND lsdtt-test-fluvial-surfaces)

# Phase 0c. Tests T2, T4, T5 and T7: the resolvability floor
# reports UNRESOLVED rather than fabricating levels, labels are a globally
# unique running counter whose band map is single-valued, k = 2 reproduces the
# legacy valley cut bit-for-bit, and k selection is a constraint with an
# extremum rather than a trade-off.
add_test(NAME unit_fluvial_surface_ladder
         COMMAND lsdtt-test-fluvial-surfaces --ladder)

# Phase 2. T3: the fitted-surface reference field. Both arms run the SAME ladder
# with the SAME auto k selection over the SAME admitted set, and only the banded
# raster differs, so what the test measures is the field rather than the tuning.
# HAND shreds the fan into radial arcs; the residual to a fitted cone does not.
add_test(NAME unit_fluvial_surface_fitted_field
         COMMAND lsdtt-test-fluvial-surfaces --fitted-surface)

# Phase 7b. Perona-Malik is adopted on the claim that it removes within-facet
# noise while leaving the riser standing, and BOTH halves need testing: a no-op
# preserves every riser and a heavy blur removes every bit of noise, so either
# assertion alone is trivially satisfiable. A 3x3 mean filter runs alongside as
# the control - if the conductance were ever turned off, P-M and the mean filter
# would agree and PM4 would catch it. Also pins that nodata is not diffused into
# and that the lambda stability limit is a refusal rather than a warning.
add_test(NAME unit_fluvial_surface_perona_malik
         COMMAND lsdtt-test-fluvial-surfaces --perona-malik)
lsdtt_set_tests_properties(unit_fluvial_surface_perona_malik PROPERTIES
    LABELS "fluvial_surfaces")

# Phase 7c. The justification for shipping a SECOND smoother is that it
# preserves the riser for a different reason - Weickert keeps it because it is
# linear, Perona-Malik because it is steep - so the test is a riser BROKEN BY
# NOTCHES and the discriminating statistic is the step height in the notches.
# Both leave the intact riser standing; only the orientation-aware one closes
# the break. If a change ever made the two behave alike, CED3/CED4 fail.
# Phase 7d. E_v is resolved PER WINDOW, so the fixture is a channel whose
# admitted ground is gentle over its first half and steep over its second with
# everything else held identical: if the resolution ever collapsed back to one
# global number both halves would report the same required_gap_m. The two
# degenerate paths are pinned too, because each is a way for the feature to look
# like it works while doing nothing - `fixed` must ignore a supplied gradient
# raster, and `slope_dependent` with none must fall back AND say so.
# Phase 7e. The bootstrap and R_bar both exist to say whether a number already
# in the output is worth quoting, so each is tested at BOTH ends: a population
# that should be stable and one that should not, a patch facing one way and one
# facing every way. A test with only the positive half would pass on an
# implementation that returned 1.0 unconditionally. RB5 additionally pins the
# circular mean against the arithmetic one - averaging 350 and 10 degrees gives
# 180, the opposite direction - and against reporting 360 rather than 0.
add_test(NAME unit_fluvial_surface_bootstrap_rbar
         COMMAND lsdtt-test-fluvial-surfaces --bootstrap-rbar)
lsdtt_set_tests_properties(unit_fluvial_surface_bootstrap_rbar PROPERTIES
    LABELS "fluvial_surfaces")

add_test(NAME unit_fluvial_surface_slope_dependent_ev
         COMMAND lsdtt-test-fluvial-surfaces --slope-dependent-ev)
lsdtt_set_tests_properties(unit_fluvial_surface_slope_dependent_ev PROPERTIES
    LABELS "fluvial_surfaces")

# Fluvial-surface methods phase 3: the Strahler channel definition (orders
# against a brute-force recursion, the scope filter in mode 1 and a user-source
# mode, HAND on that channel against a receiver walk including the undefined
# case) and the swath-baseline relief against an all-pairs nearest search.
add_test(NAME unit_fluvial_surface_channel_order_swath
         COMMAND lsdtt-test-fluvial-surfaces --channel-order-swath)
lsdtt_set_tests_properties(unit_fluvial_surface_channel_order_swath PROPERTIES
    LABELS "fluvial_surfaces")

# Phase 4: the floodplain / terrace engine on a hand-built grid where every
# count is known - channel contact, size, seeding, holes, the gradient bound,
# the terrace window, the D2 overlap refusals, and the hillslope-coupling rule.
add_test(NAME unit_fluvial_surface_floodplain_thresholds
         COMMAND lsdtt-test-fluvial-surfaces --floodplain-thresholds)
lsdtt_set_tests_properties(unit_fluvial_surface_floodplain_thresholds PROPERTIES
    LABELS "fluvial_surfaces")

# Phase 5: cross-section levels - the ladder's rule on a population
# (partition_population), clustering section breaks with the support filter,
# the levels file, and the ladder banding on supplied breaks (status supplied,
# the floor still measured, an empty band, a missing channel, source_id -1).
add_test(NAME unit_fluvial_surface_cross_section_levels
         COMMAND lsdtt-test-fluvial-surfaces --cross-section-levels
                 "${CMAKE_CURRENT_BINARY_DIR}/unit_fluvial_surface_cross_section_levels")
lsdtt_set_tests_properties(unit_fluvial_surface_cross_section_levels PROPERTIES
    LABELS "fluvial_surfaces")

add_test(NAME unit_fluvial_surface_coherence_diffusion
         COMMAND lsdtt-test-fluvial-surfaces --coherence-diffusion)
lsdtt_set_tests_properties(unit_fluvial_surface_coherence_diffusion PROPERTIES
    LABELS "fluvial_surfaces")

add_test(NAME fixture_fluvial_surfaces_cli
         COMMAND lsdtt-test-fluvial-surfaces --write-cli-fixture
                 ${LSDTT_FLUVIAL_SURFACES_FIXTURE_DIR})
lsdtt_set_tests_properties(fixture_fluvial_surfaces_cli PROPERTIES
    LABELS "fixture"
    FIXTURES_SETUP fluvial_surfaces_fixture)

# The generators are reachable from the user-facing CLI, not only from the test
# binary. Five rasters, because the ground-truth tread_id / lobe_id layers are
# what every later assertion about purity is scored against.
add_lsdtt_cli_smoke_test(
    smoke_fluvial_surfaces_generators
    lsdtt-synthetic-surfaces
    synthetic_terrace_fan.param.in
    synfs_terrace_staircase.tif
    --expect synfs_terrace_staircase_tread_id.tif
    --expect synfs_terrace_staircase_channel_mask.tif
    --expect synfs_alluvial_fan.tif
    --expect synfs_alluvial_fan_lobe_id.tif
    --expect-raster-size synfs_terrace_staircase.tif:300x200
    --expect-raster-stat synfs_terrace_staircase_tread_id.tif:max:4.0:4.0
    --expect-raster-stat synfs_alluvial_fan_lobe_id.tif:max:2.0:2.0
)

add_test(NAME smoke_flat_surface_finder_help
         COMMAND lsdtt-flat-surface-finder -h)
lsdtt_set_tests_properties(smoke_flat_surface_finder_help PROPERTIES
    LABELS "cli"
    PASS_REGULAR_EXPRESSION "max_relief_metres")

add_test(NAME smoke_valley_centreline_help
         COMMAND lsdtt-valley-centreline -h)
lsdtt_set_tests_properties(smoke_valley_centreline_help PROPERTIES
    LABELS "cli"
    PASS_REGULAR_EXPRESSION "valley_mask_fname")

add_lsdtt_cli_smoke_test(
    smoke_flat_surface_finder
    lsdtt-flat-surface-finder
    flat_surface_finder_smoke.param.in
    fsf_labelled_flat_mask.tif
    --read-path "${LSDTT_FLUVIAL_SURFACES_FIXTURE_DIR}"
    --expect fsf_combined_flat_mask.tif
    --expect fsf_relief_cdf.csv
)
lsdtt_set_tests_properties(smoke_flat_surface_finder PROPERTIES
    FIXTURES_REQUIRED fluvial_surfaces_fixture)

# Riser bridging over many endpoints at once: moved from
# smoke_flat_surface_finder_slope_break_connect with the riser finder
# (2026-10-26). See the parameter file.
add_lsdtt_cli_smoke_test(
    smoke_feature_finder_risers_connect
    lsdtt-feature-finder
    feature_finder_risers_connect.param.in
    ffrc_riser_centreline.tif
    --read-path "${LSDTT_FLUVIAL_SURFACES_FIXTURE_DIR}"
    --expect ffrc_riser_gap_closing.tif
    --stdout-contains "endpoints found : 39"
    --stdout-contains "candidates      : 27"
    --stdout-contains "bridges accepted: 12"
)
lsdtt_set_tests_properties(smoke_feature_finder_risers_connect PROPERTIES
    FIXTURES_REQUIRED fluvial_surfaces_fixture)

# Phase 1. The reordering through the CLI: the ladder must resolve the flight,
# every patch must carry a band, and the run must report which windows resolved.
# Phase 1's reordering, now driven through lsdtt-fluvial-surfaces: phase 4 moved
# the whole height-reference half of the pipeline out of CLI 1.
add_lsdtt_cli_smoke_test(
    smoke_fluvial_surfaces_band_then_label
    lsdtt-fluvial-surfaces
    fluvial_surfaces_band_then_label.param.in
    btl_banded_labels.tif
    --read-path "${LSDTT_FLUVIAL_SURFACES_FIXTURE_DIR}"
    --expect btl_band_index.tif
    --expect btl_surface_patches.csv
    --expect btl_ladder_windows.csv
    --expect btl_surface_class.tif
    --stdout-contains "Ladder and band-then-label"
    --expect-raster-values "btl_band_index.tif:0,1,2,3"
    --expect-file-contains "btl_ladder_windows.csv:source_id,window_centre_m,n_pixels,k_selected,status,used_global_breaks,breaks_m,band_centres_m,min_band_centre_gap_m,required_gap_m,k_bootstrap_mode,k_stability,k_ci_lo,k_ci_hi,n_bootstrap,vertical_accuracy_m,median_gradient,goodness_of_variance_fit,separation_ratio"
    --expect-file-contains "btl_surface_patches.csv:label,band,class,source_id,n_pixels,area_m2,median_field_m,aspect_degrees,terrain_aspect_degrees,aspect_r_bar,window_centre_m"
)
lsdtt_set_tests_properties(smoke_fluvial_surfaces_band_then_label PROPERTIES
    FIXTURES_REQUIRED fluvial_surfaces_fixture)

# The label is the join key for everything band-then-label emits, and that is exactly the property no per-file assertion can see: the
# smoke test pins each table's SCHEMA, but a schema says nothing about whether a
# label in the raster has a row in the patch table, or whether the two agree
# about which band it is in. So these two run the CLI and then check ACROSS the
# outputs.
#
# The join arm also carries the only coverage of
# write_reference_field_profile_csv. It is off by default and no other fixture
# turns it on, so its label and band columns had never been checked against the
# rasters they claim to describe.
if(LSDTT_ENABLE_PYTHON_GEO_TESTS)
    add_test(NAME rigorous_fluvial_surfaces_label_band_identity
        COMMAND ${Python3_EXECUTABLE}
            "${CMAKE_CURRENT_SOURCE_DIR}/scripts/check_fluvial_surface_joins.py"
            --exe "$<TARGET_FILE:lsdtt-fluvial-surfaces>"
            --param-template "${LSDTT_TEST_PARAM_DIR}/fluvial_surfaces_join_check.param.in"
            --read-path "${LSDTT_FLUVIAL_SURFACES_FIXTURE_DIR}"
            --work-dir "${CMAKE_CURRENT_BINARY_DIR}/rigorous_fluvial_surfaces_label_band_identity"
            --stem join)
    lsdtt_mark_test_runtime(python rigorous_fluvial_surfaces_label_band_identity)
else()
    lsdtt_mark_test_unregistered(python_geo rigorous_fluvial_surfaces_label_band_identity)
endif()
lsdtt_set_tests_properties(rigorous_fluvial_surfaces_label_band_identity PROPERTIES
    LABELS "cli"
    FIXTURES_REQUIRED fluvial_surfaces_fixture)

# The same join checks on clubb2017 output: the fixed_thresholds path builds its
# own patch and class records, and they must join exactly as the ladder's do.
if(LSDTT_ENABLE_PYTHON_GEO_TESTS)
    add_test(NAME rigorous_fluvial_surfaces_clubb2017_joins
        COMMAND ${Python3_EXECUTABLE}
            "${CMAKE_CURRENT_SOURCE_DIR}/scripts/check_fluvial_surface_joins.py"
            --exe "$<TARGET_FILE:lsdtt-fluvial-surfaces>"
            --param-template "${LSDTT_TEST_PARAM_DIR}/fluvial_surfaces_clubb2017.param.in"
            --read-path "${LSDTT_FLUVIAL_SURFACES_FIXTURE_DIR}"
            --work-dir "${CMAKE_CURRENT_BINARY_DIR}/rigorous_fluvial_surfaces_clubb2017_joins"
            --stem clubb)
    lsdtt_mark_test_runtime(python rigorous_fluvial_surfaces_clubb2017_joins)
else()
    lsdtt_mark_test_unregistered(python_geo rigorous_fluvial_surfaces_clubb2017_joins)
endif()
lsdtt_set_tests_properties(rigorous_fluvial_surfaces_clubb2017_joins PROPERTIES
    LABELS "cli"
    FIXTURES_REQUIRED fluvial_surfaces_fixture)

# A table with no rows and a table that was never written look identical to a
# downstream reader. This pins that the empty result is the former.
if(LSDTT_ENABLE_PYTHON_GEO_TESTS)
    add_test(NAME rigorous_fluvial_surfaces_empty_result
        COMMAND ${Python3_EXECUTABLE}
            "${CMAKE_CURRENT_SOURCE_DIR}/scripts/check_fluvial_surface_joins.py"
            --exe "$<TARGET_FILE:lsdtt-fluvial-surfaces>"
            --param-template "${LSDTT_TEST_PARAM_DIR}/fluvial_surfaces_empty_result.param.in"
            --read-path "${LSDTT_FLUVIAL_SURFACES_FIXTURE_DIR}"
            --work-dir "${CMAKE_CURRENT_BINARY_DIR}/rigorous_fluvial_surfaces_empty_result"
            --stem empty
            --expect-empty-tables)
    lsdtt_mark_test_runtime(python rigorous_fluvial_surfaces_empty_result)
else()
    lsdtt_mark_test_unregistered(python_geo rigorous_fluvial_surfaces_empty_result)
endif()
lsdtt_set_tests_properties(rigorous_fluvial_surfaces_empty_result PROPERTIES
    LABELS "cli"
    FIXTURES_REQUIRED fluvial_surfaces_fixture)

# Phase 7a. `admission_algorithm_choice` exists because local relief and a
# gradient bound reject different ground; if a refactor ever collapsed them into
# one test, or wired relief_or_gradient back to the relief mask alone, every
# count would still look plausible. The synthetic DEM is built to make the two
# disagree - a smooth 8 percent ramp that only the gradient bound admits, a
# rough flat that only relief could, a smooth flat that both take and a
# gradient-1 slope that neither does - and the assertions are on the SET
# RELATIONS between the four modes, which survive the morphological
# post-processing that the raw counts do not. It also runs the pre-2026-10-26
# admission keys and requires byte-identical output to the canonical ones.
if(LSDTT_ENABLE_PYTHON_GEO_TESTS)
    add_test(NAME rigorous_flat_surface_finder_admission_modes
        COMMAND ${Python3_EXECUTABLE}
            "${CMAKE_CURRENT_SOURCE_DIR}/scripts/check_admission_modes.py"
            --exe "$<TARGET_FILE:lsdtt-flat-surface-finder>"
            --work-dir "${CMAKE_CURRENT_BINARY_DIR}/rigorous_flat_surface_finder_admission_modes"
            --cap-degrees 15.0)
    lsdtt_mark_test_runtime(python rigorous_flat_surface_finder_admission_modes)
else()
    lsdtt_mark_test_unregistered(python_geo rigorous_flat_surface_finder_admission_modes)
endif()
lsdtt_set_tests_properties(rigorous_flat_surface_finder_admission_modes PROPERTIES
    LABELS "cli;fluvial_surfaces")

# The trough's upstream-of-head exclusion. Terrain beyond the channel head stays
# in the valley MASK but must not be carved into the trough, or the pocket
# behind the head becomes a distance-to-wall bowl that spills over the upstream
# tile edge.
#
# This needs its own fixture because every smoke fixture has its channel head at
# a DEM edge, where there is no ground beyond the head at all and the exclusion
# count is legitimately zero - so none of them can tell a working guard from a
# dead one. The guard went inert at phase 5d and nothing noticed.
add_test(NAME unit_fluvial_surfaces_trough_upstream
         COMMAND lsdtt-test-fluvial-surfaces --trough-upstream)

add_test(NAME unit_fluvial_surfaces_band_then_label
         COMMAND lsdtt-test-fluvial-surfaces --band-then-label)

# Phase 4. The two retirement groups must behave DIFFERENTLY and visibly so: a
# key that turned an analysis on stops the run, a key that merely tuned a
# vanished mechanism warns and continues. If these two ever became confusable
# the whole two-group design would be pointless.
add_lsdtt_cli_smoke_test(
    failure_flat_surface_finder_removed_parameter
    lsdtt-flat-surface-finder
    flat_surface_finder_removed_parameter.param.in
    unused.tif
    --read-path "${LSDTT_FLUVIAL_SURFACES_FIXTURE_DIR}"
    --expect-failure
    --stderr-contains "find_floodplain has moved to lsdtt-fluvial-surfaces"
)
lsdtt_set_tests_properties(failure_flat_surface_finder_removed_parameter PROPERTIES
    FIXTURES_REQUIRED fluvial_surfaces_fixture)

# The riser finder's keys after its move to lsdtt-feature-finder: an output
# switch stops the run naming the new home, a tuning key warns and continues.
add_lsdtt_cli_smoke_test(
    failure_flat_surface_finder_removed_riser_parameter
    lsdtt-flat-surface-finder
    flat_surface_finder_removed_riser_parameter.param.in
    unused.tif
    --read-path "${LSDTT_FLUVIAL_SURFACES_FIXTURE_DIR}"
    --expect-failure
    --stderr-contains "write_slope_break_centreline_raster has moved"
    --stderr-contains "lsdtt-feature-finder (feature_type: find_risers)"
)
lsdtt_set_tests_properties(failure_flat_surface_finder_removed_riser_parameter PROPERTIES
    FIXTURES_REQUIRED fluvial_surfaces_fixture)

add_lsdtt_cli_smoke_test(
    smoke_flat_surface_finder_ignored_riser_parameter
    lsdtt-flat-surface-finder
    flat_surface_finder_ignored_riser_parameter.param.in
    ignr_combined_flat_mask.tif
    --read-path "${LSDTT_FLUVIAL_SURFACES_FIXTURE_DIR}"
    --stderr-contains "vec_adaptive_C is ignored"
    --stderr-contains "where it is riser_vec_adaptive_C"
    --stderr-contains "where it is riser_connect_search_radius_metres"
)
lsdtt_set_tests_properties(smoke_flat_surface_finder_ignored_riser_parameter PROPERTIES
    FIXTURES_REQUIRED fluvial_surfaces_fixture)

add_lsdtt_cli_smoke_test(
    smoke_flat_surface_finder_ignored_parameter
    lsdtt-flat-surface-finder
    flat_surface_finder_ignored_parameter.param.in
    ign_combined_flat_mask.tif
    --read-path "${LSDTT_FLUVIAL_SURFACES_FIXTURE_DIR}"
    --stderr-contains "enable_hand_banding is ignored"
    --stderr-contains "the reference-field ladder in lsdtt-fluvial-surfaces is now unconditional"
    --stderr-contains "terrace_min_height_metres has moved to lsdtt-fluvial-surfaces"
    --stderr-contains "'write_combined_flat_mask' is deprecated"
)
lsdtt_set_tests_properties(smoke_flat_surface_finder_ignored_parameter PROPERTIES
    FIXTURES_REQUIRED fluvial_surfaces_fixture)

# =============================================================================
# Parameter key policy (phase 3b): deprecated, removed and unknown keys
#
# The prerequisite for phase 4, where CLI 1 first sheds parameters. Until this
# existed a removed key was silently a no-op - the typed parses iterate over the
# DEFAULTS MAPS, so a key in the file that no map declares is never looked at by
# anything, and an old parameter file would run happily having done none of what
# it asked for.
#
# The driver re-executes itself as a child for each arm, because the removed-key
# case exits the process.
# =============================================================================
add_executable(lsdtt-test-param-deprecation cpp/test_param_deprecation.cpp)
target_link_libraries(lsdtt-test-param-deprecation PRIVATE lsdtt3_core)

add_test(NAME unit_param_parser_deprecation
         COMMAND lsdtt-test-param-deprecation
                 "${CMAKE_CURRENT_BINARY_DIR}/param_deprecation")
# The domain label is explicit: the compartmentalisation block at the foot of
# this file assigns domains by NAME match, and this test is called
# param_parser rather than fluvial_surface. It is nonetheless the guard on a
# facility built for this workstream, so it belongs in that compartment.
lsdtt_set_tests_properties(unit_param_parser_deprecation PROPERTIES
    LABELS "fluvial_surfaces")

# =============================================================================
# unit_param_parser_strict
#
# Where the parser finds the parameter file (a filename with its own directory,
# absolute or relative, sets read_path to that directory) and that the
# placeholder constructor every CLI declares prints nothing.
# =============================================================================
add_executable(lsdtt-test-param-parser-strict cpp/test_param_parser_strict.cpp)
target_link_libraries(lsdtt-test-param-parser-strict PRIVATE lsdtt3_core)

add_test(NAME unit_param_parser_strict
         COMMAND lsdtt-test-param-parser-strict
                 "${CMAKE_CURRENT_BINARY_DIR}/param_parser_strict")

# =============================================================================
# unit_run_record
#
# How every CLI run ends (lsdtt_run_record): the finished-OK line and run
# summary on success; the RUN FAILED box, ERROR line and failed summary after a
# return 1, a std::exit(1) and an uncaught exception; which files the summary
# lists as outputs.
# =============================================================================
add_executable(lsdtt-test-run-record cpp/test_run_record.cpp)
target_link_libraries(lsdtt-test-run-record PRIVATE lsdtt3_core)

add_test(NAME unit_run_record
         COMMAND lsdtt-test-run-record "${CMAKE_CURRENT_BINARY_DIR}/run_record")

# =============================================================================
# lsdtt-fluvial-surfaces
#
# The tool between lsdtt-flat-surface-finder (which ground is flat?) and
# lsdtt-valley-centreline (trace this valley): it does everything that needs a
# height reference, and it writes the valley mask CLI 3 consumes.
# =============================================================================
add_test(NAME smoke_fluvial_surfaces_help
         COMMAND lsdtt-fluvial-surfaces -h)
lsdtt_set_tests_properties(smoke_fluvial_surfaces_help PROPERTIES
    LABELS "cli"
    PASS_REGULAR_EXPRESSION "fill_before_routing")

foreach(threshold_help_target IN ITEMS
        drainage-analysis
        fluvial-surfaces
        valley-centreline
        stream-burning
        feature-finder)
    string(REPLACE "-" "_" threshold_help_slug "${threshold_help_target}")
    add_test(NAME smoke_${threshold_help_slug}_threshold_help
             COMMAND lsdtt-${threshold_help_target} -h)
    lsdtt_set_tests_properties(smoke_${threshold_help_slug}_threshold_help PROPERTIES
        LABELS "cli"
        PASS_REGULAR_EXPRESSION "threshold_contributing_pixels")
endforeach()

add_test(NAME smoke_feature_finder_fill_help
         COMMAND lsdtt-feature-finder -h)
lsdtt_set_tests_properties(smoke_feature_finder_fill_help PROPERTIES
    LABELS "cli"
    PASS_REGULAR_EXPRESSION "raster_is_filled")

# End to end, with no flat_mask_fname, so the built-in fallback admission is
# exercised too: the tool has to work on its own as well as downstream of CLI 1.
add_lsdtt_cli_smoke_test(
    smoke_fluvial_surfaces
    lsdtt-fluvial-surfaces
    fluvial_surfaces_smoke.param.in
    fsurf_band_index.tif
    --read-path "${LSDTT_FLUVIAL_SURFACES_FIXTURE_DIR}"
    --expect fsurf_banded_labels.tif
    --expect fsurf_valley_mask.tif
    --expect fsurf_surface_class.tif
    --expect fsurf_surface_patches.csv
    --expect fsurf_ladder_windows.csv
    --expect fsurf_reference_field.tif
    --stdout-contains "Ladder and band-then-label"
    --stdout-contains "window(s) resolved"
    --expect-raster-values "fsurf_band_index.tif:0,1,2,3"
    --expect-file-contains "fsurf_surface_patches.csv:label,band,class,source_id,n_pixels,area_m2,median_field_m,aspect_degrees,terrain_aspect_degrees,aspect_r_bar,window_centre_m"
    --expect-file-contains "fsurf_ladder_windows.csv:source_id,window_centre_m,n_pixels,k_selected,status,used_global_breaks,breaks_m,band_centres_m,min_band_centre_gap_m,required_gap_m,k_bootstrap_mode,k_stability,k_ci_lo,k_ci_hi,n_bootstrap,vertical_accuracy_m,median_gradient,goodness_of_variance_fit,separation_ratio"
)
lsdtt_set_tests_properties(smoke_fluvial_surfaces PROPERTIES
    FIXTURES_REQUIRED fluvial_surfaces_fixture)

# The warning names both spellings, while the routing diagnostic is emitted
# immediately before the branch and proves the aliased false value skips fill.
add_lsdtt_cli_smoke_test(
    rigorous_fluvial_surfaces_legacy_fill_key
    lsdtt-fluvial-surfaces
    fluvial_surfaces_fill_legacy.param.in
    fsurf_fill_legacy_band_index.tif
    --read-path "${LSDTT_FLUVIAL_SURFACES_FIXTURE_DIR}"
    --stdout-contains "Fill before routing  : no"
    --stdout-contains "Routing on the unfilled input DEM."
    --stderr-contains "'fill_dem_for_flow_routing' is deprecated"
    --stderr-contains "use 'fill_before_routing'"
    --stderr-contains "'channel_threshold_contributing_pixels' is deprecated"
    --stderr-contains "use 'threshold_contributing_pixels'"
    --expect-file-contains "fsurf_fill_legacy_ingested_params.param:threshold_contributing_pixels: 500"
)
lsdtt_set_tests_properties(rigorous_fluvial_surfaces_legacy_fill_key PROPERTIES
    FIXTURES_REQUIRED fluvial_surfaces_fixture)

# THE composability guarantee: HAND is the expensive stage and the level count
# is what gets iterated on, so a supplied HAND must reproduce the computed path
# exactly or every iteration is a different experiment from the final run.
# Two runs compared by checksum, which no single-run assertion can express.
if(LSDTT_ENABLE_PYTHON_TESTS)
    add_test(NAME smoke_fluvial_surfaces_hand_input
        COMMAND ${Python3_EXECUTABLE}
            "${CMAKE_CURRENT_SOURCE_DIR}/scripts/check_fluvial_surfaces_hand_roundtrip.py"
            --exe "$<TARGET_FILE:lsdtt-fluvial-surfaces>"
            --param-template "${LSDTT_TEST_PARAM_DIR}/fluvial_surfaces_hand_input.param.in"
            --read-path "${LSDTT_FLUVIAL_SURFACES_FIXTURE_DIR}"
            --work-dir "${CMAKE_CURRENT_BINARY_DIR}/smoke_fluvial_surfaces_hand_input"
            --stem fsin)
    lsdtt_mark_test_runtime(python smoke_fluvial_surfaces_hand_input)
else()
    lsdtt_mark_test_unregistered(python smoke_fluvial_surfaces_hand_input)
endif()
lsdtt_set_tests_properties(smoke_fluvial_surfaces_hand_input PROPERTIES
    LABELS "cli"
    FIXTURES_REQUIRED fluvial_surfaces_fixture)

# The fitted-surface field through the CLI. The param file deliberately does not
# name a break algorithm, so this also pins that selecting fitted_surface
# switches the finder to jenks on its own and says so - antimode bins on log10
# and cannot band a signed residual.
add_lsdtt_cli_smoke_test(
    smoke_fluvial_surfaces_fitted_surface
    lsdtt-fluvial-surfaces
    fluvial_surfaces_fitted_surface.param.in
    fsfan_band_index.tif
    --read-path "${LSDTT_FLUVIAL_SURFACES_FIXTURE_DIR}"
    --expect fsfan_banded_labels.tif
    --expect fsfan_valley_mask.tif
    --expect fsfan_surface_patches.csv
    --expect fsfan_reference_field.tif
    --stdout-contains "the break finder is jenks"
    --stdout-contains "cone fit"
    --expect-raster-values "fsfan_band_index.tif:0,1"
)
lsdtt_set_tests_properties(smoke_fluvial_surfaces_fitted_surface PROPERTIES
    FIXTURES_REQUIRED fluvial_surfaces_fixture)

# Output pins for lsdtt-fluvial-surfaces, recorded 2026-10-26 before the method
# consolidation (Clubb 2017, threshold pickers, parameter renames) touches the
# CLI. The smoke tests above pin schemas and band VALUES; these pin the numbers,
# so a later change can show it moved nothing it did not mean to. Each run was
# checked identical at n_threads 1 and 8 before pinning.
#
# A moved pin is not automatically a regression, but it must be explained
# before it is re-blessed. The float outputs (HAND, the fitted residual, the
# CSV medians and gaps) are the ones a different compiler or libm could move by
# an ulp; the integer rasters (band, label, class, valley mask) should not move
# unless a decision moved. The README is not pinned: it carries a timestamp.
#
# Re-blessed once, 2026-10-26: the three _ladder_windows.csv hashes, for the
# column rename median_slope_degrees -> median_gradient (gradients in m/m). With
# the header renamed back each file is byte-identical to its first pin; no
# raster moved.

# Fallback admission, HAND ladder, the staircase. Same template as
# smoke_fluvial_surfaces.
add_lsdtt_cli_smoke_test(
    rigorous_fluvial_surfaces_checksums
    lsdtt-fluvial-surfaces
    fluvial_surfaces_smoke.param.in
    fsurf_band_index.tif
    --read-path "${LSDTT_FLUVIAL_SURFACES_FIXTURE_DIR}"
    --expect-raster-checksum fsurf_band_index.tif:57493
    --expect-raster-checksum fsurf_banded_labels.tif:55857
    --expect-raster-checksum fsurf_surface_class.tif:37966
    --expect-raster-checksum fsurf_valley_mask.tif:972
    --expect-raster-checksum fsurf_hand.tif:57881
    --expect-raster-checksum fsurf_reference_field.tif:57881
    --expect-sha256 fsurf_ladder_windows.csv:b8113840463e5647c5b3b289de66fad93fa267c33e04ba47b042289bfff92cf1
    --expect-sha256 fsurf_surface_patches.csv:8a259c130ce3da3b69756e48c300e6059fc0b80348d0b7b29b4e753bdce4008f
    --expect-sha256 fsurf_surface_class.csv:e6f63f2c5908ed0a36a1e30a61224b1ee3a23da7e1630a3913d817ef495e0cf4
)
lsdtt_set_tests_properties(rigorous_fluvial_surfaces_checksums PROPERTIES
    FIXTURES_REQUIRED fluvial_surfaces_fixture)

# Fitted-surface (cone) reference field on the fan. Same template as
# smoke_fluvial_surfaces_fitted_surface.
add_lsdtt_cli_smoke_test(
    rigorous_fluvial_surfaces_fitted_surface_checksums
    lsdtt-fluvial-surfaces
    fluvial_surfaces_fitted_surface.param.in
    fsfan_band_index.tif
    --read-path "${LSDTT_FLUVIAL_SURFACES_FIXTURE_DIR}"
    --expect-raster-checksum fsfan_band_index.tif:6394
    --expect-raster-checksum fsfan_banded_labels.tif:26071
    --expect-raster-checksum fsfan_valley_mask.tif:24826
    --expect-raster-checksum fsfan_reference_field.tif:2938
    --expect-sha256 fsfan_ladder_windows.csv:e98845c4c205158dc633b3f7f067a3295736abff2b3c18b20bc13027eba4e120
    --expect-sha256 fsfan_surface_patches.csv:ec68e8291333913c801c26b7ec251decaa0c374e1478857c58f9dced03ceb0c3
)
lsdtt_set_tests_properties(rigorous_fluvial_surfaces_fitted_surface_checksums PROPERTIES
    FIXTURES_REQUIRED fluvial_surfaces_fixture)

# The two-step workflow: lsdtt-flat-surface-finder writes the flat mask, then
# lsdtt-fluvial-surfaces bands it. This is the default admission path
# (flat_mask_fname set), which no other test reaches. Stage 1 is a fixture of
# stage 2, so `ctest -R rigorous_fluvial_surfaces_flat_mask_input` runs both.
add_lsdtt_cli_smoke_test(
    rigorous_flat_surface_finder_staircase_mask
    lsdtt-flat-surface-finder
    flat_surface_finder_staircase_mask.param.in
    fsfst_combined_flat_mask.tif
    --read-path "${LSDTT_FLUVIAL_SURFACES_FIXTURE_DIR}"
    --expect-raster-checksum fsfst_combined_flat_mask.tif:21440
)
lsdtt_set_tests_properties(rigorous_flat_surface_finder_staircase_mask PROPERTIES
    FIXTURES_REQUIRED fluvial_surfaces_fixture
    FIXTURES_SETUP fluvial_surfaces_staircase_flat_mask)

add_lsdtt_cli_smoke_test(
    rigorous_fluvial_surfaces_flat_mask_input
    lsdtt-fluvial-surfaces
    fluvial_surfaces_flat_mask_input.param.in
    fsfm_band_index.tif
    --read-path "${LSDTT_FLUVIAL_SURFACES_FIXTURE_DIR}"
    --stdout-contains "Loading flat mask"
    --expect-raster-checksum fsfm_band_index.tif:55597
    --expect-raster-checksum fsfm_banded_labels.tif:57592
    --expect-raster-checksum fsfm_surface_class.tif:37252
    --expect-raster-checksum fsfm_valley_mask.tif:1128
    --expect-raster-checksum fsfm_reference_field.tif:57881
    --expect-sha256 fsfm_ladder_windows.csv:e8986691f8cbbd2c50f3cff34009f57c603d8a6fd30d510cc9129120da304b96
    --expect-sha256 fsfm_surface_patches.csv:20dae3570717df72dfb17d6b40086a6df84d5afff05b470fc1f2323d4a600e3b
    --expect-sha256 fsfm_surface_class.csv:dbf403d3098b67472eefb257cfd05bf9937ad39688d584c57978eabc86a7b9bb
)
lsdtt_set_tests_properties(rigorous_fluvial_surfaces_flat_mask_input PROPERTIES
    FIXTURES_REQUIRED "fluvial_surfaces_fixture;fluvial_surfaces_staircase_flat_mask")

# The renames of 2026-10-26 must carry every value: fluvial_surfaces_smoke under
# the old names, held to exactly the pins of rigorous_fluvial_surfaces_checksums,
# with the converted values recorded in the ingested file so a rerun from it
# reproduces the run without the old keys.
add_lsdtt_cli_smoke_test(
    rigorous_fluvial_surfaces_legacy_keys
    lsdtt-fluvial-surfaces
    fluvial_surfaces_legacy_keys.param.in
    fsurf_band_index.tif
    --read-path "${LSDTT_FLUVIAL_SURFACES_FIXTURE_DIR}"
    --expect-raster-checksum fsurf_band_index.tif:57493
    --expect-raster-checksum fsurf_banded_labels.tif:55857
    --expect-raster-checksum fsurf_surface_class.tif:37966
    --expect-raster-checksum fsurf_valley_mask.tif:972
    --expect-raster-checksum fsurf_hand.tif:57881
    --expect-sha256 fsurf_ladder_windows.csv:b8113840463e5647c5b3b289de66fad93fa267c33e04ba47b042289bfff92cf1
    --expect-sha256 fsurf_surface_patches.csv:8a259c130ce3da3b69756e48c300e6059fc0b80348d0b7b29b4e753bdce4008f
    --expect-sha256 fsurf_surface_class.csv:e6f63f2c5908ed0a36a1e30a61224b1ee3a23da7e1630a3913d817ef495e0cf4
    --stderr-contains "'fallback_admission_algorithm_choice' is deprecated"
    --stderr-contains "use 'admission_algorithm_choice'"
    --stderr-contains "'fallback_relief_window_radius_pixels' is deprecated"
    --stderr-contains "use 'relief_window_radius_metres'"
    --stderr-contains "'vertical_accuracy_break_degrees' is deprecated"
    --stderr-contains "use 'vertical_accuracy_break_gradient'. Converted with gradient = tan(angle)."
    --stderr-contains "'min_component_pixels' is deprecated"
    --stderr-contains "use 'min_band_component_pixels'"
    --stdout-contains "vertical_accuracy_break_degrees 11.31 -> vertical_accuracy_break_gradient 0.200001225"
    --expect-file-contains "fsurf_ingested_params.param:relief_window_radius_metres: 50"
    --expect-file-contains "fsurf_ingested_params.param:admission_algorithm_choice: relief"
    --expect-file-contains "fsurf_ingested_params.param:vertical_accuracy_break_gradient: 0.200001225"
    --expect-file-not-contains "fsurf_ingested_params.param:fallback_"
    --expect-file-not-contains "fsurf_ingested_params.param:legacy_relief_window_radius_pixels"
    --expect-file-not-contains "fsurf_ingested_params.param:legacy_vertical_accuracy_break_degrees"
    --expect-file-not-contains "fsurf_ingested_params.param:legacy_fallback_admission_algorithm_choice"
)
lsdtt_set_tests_properties(rigorous_fluvial_surfaces_legacy_keys PROPERTIES
    FIXTURES_REQUIRED fluvial_surfaces_fixture)

# ---- Fluvial-surface methods phase 4: Clubb et al. (2017) as clubb2017 ------
# On the linear-range surface: a branching network (order up to 5 at 50 px)
# draining rough mountains onto a flat piedmont. Pins recorded 2026-10-26,
# identical at n_threads 1 and 8. The run log's stage table is asserted so a
# preset that silently stopped applying would fail here.
# Re-blessed 2026-10-27: fill_floodplain_holes became true under the preset;
# 10 enclosed cells join the floodplain (mask 8344 -> 8354, the valley mask's
# value). The QQ tables, HAND and the terrace field did not move.
add_lsdtt_cli_smoke_test(
    rigorous_fluvial_surfaces_clubb2017
    lsdtt-fluvial-surfaces
    fluvial_surfaces_clubb2017.param.in
    clubb_band_index.tif
    --read-path "${LSDTT_FLUVIAL_SURFACES_FIXTURE_DIR}"
    --stdout-contains "floodplain_thresholds   (preset clubb2017)"
    --stdout-contains "fixed_thresholds   (preset clubb2017)"
    --stdout-contains "true   (preset clubb2017)"
    --stdout-contains "low relief test: 0.383436 of the DEM's floodplain gradient is below 0.1 (switch above 0.8) -> hand_and_gradient"
    --expect-file-contains "clubb_ingested_params.param:fill_floodplain_holes: true"
    --stdout-contains "channel = Strahler order >= 3: 669 channel pixel(s)"
    --stdout-contains "HAND: quantile 2379 of 10000"
    --stdout-contains "gradient: quantile 2396 of 10000"
    --stdout-contains "kept 3 component(s), dropped 1 small and 0 not touching the channel"
    --expect-raster-values "clubb_band_index.tif:0,1"
    --expect-raster-checksum clubb_band_index.tif:34254
    --expect-raster-checksum clubb_banded_labels.tif:18367
    --expect-raster-checksum clubb_surface_class.tif:8718
    --expect-raster-checksum clubb_floodplain_mask.tif:8354
    --expect-raster-checksum clubb_valley_mask.tif:8354
    --expect-raster-checksum clubb_terrace_field.tif:45462
    --expect-raster-checksum clubb_reference_field.tif:14348
    --expect-sha256 clubb_surface_patches.csv:66ea47fdfaa5eecb8853e65ebd0ff7e485ff248c7c661c9e1050ce8dcf626bea
    --expect-sha256 clubb_surface_class.csv:207af4fac576515f48c55a72590f0cecf05ac2290661272ffaf3fbffe438bc91
    --expect-sha256 clubb_floodplain_hand_qq.csv:a67f79bbbc740d057f0f8021d3db835fa29154deb6b62051d2980c219084939d
    --expect-sha256 clubb_floodplain_gradient_qq.csv:58535ca9625e09c762d14a9af6aa693ea1c37dab381e295426901c5bdffc498c
    --expect-file-contains "clubb_ingested_params.param:max_floodplain_hand_metres: 0.667831123"
    --expect-file-contains "clubb_ingested_params.param:max_floodplain_gradient: 0.0158595201"
)
lsdtt_set_tests_properties(rigorous_fluvial_surfaces_clubb2017 PROPERTIES
    FIXTURES_REQUIRED fluvial_surfaces_fixture)

# run_terrace_classification false: the floodplain alone. Same run otherwise,
# so the floodplain mask and the QQ tables carry the pins above; no terrace
# field is built or written and no cell is classed terrace.
add_lsdtt_cli_smoke_test(
    rigorous_fluvial_surfaces_clubb2017_floodplain_only
    lsdtt-fluvial-surfaces
    fluvial_surfaces_clubb2017_floodplain_only.param.in
    clubbfp_floodplain_mask.tif
    --read-path "${LSDTT_FLUVIAL_SURFACES_FIXTURE_DIR}"
    --stdout-contains "run_terrace_classification false: floodplain only"
    --expect-raster-checksum clubbfp_floodplain_mask.tif:8354
    --expect-raster-values "clubbfp_band_index.tif:0"
    --expect-absent clubbfp_terrace_field.tif
    --expect-sha256 clubbfp_floodplain_hand_qq.csv:a67f79bbbc740d057f0f8021d3db835fa29154deb6b62051d2980c219084939d
    --expect-sha256 clubbfp_floodplain_gradient_qq.csv:58535ca9625e09c762d14a9af6aa693ea1c37dab381e295426901c5bdffc498c
)
lsdtt_set_tests_properties(rigorous_fluvial_surfaces_clubb2017_floodplain_only PROPERTIES
    FIXTURES_REQUIRED fluvial_surfaces_fixture)

# floodplain_second_pass_algorithm_choice qq_upper_departure: a second HAND
# bound from the upper departure of the QQ over the first-pass floodplain, the
# floodplain rebuilt at it and filled. The first-pass mask keeps the clubb2017
# floodplain pin (8354); the pass and polygon counts were recorded 2026-10-27.
add_lsdtt_cli_smoke_test(
    rigorous_fluvial_surfaces_clubb2017_second_pass
    lsdtt-fluvial-surfaces
    fluvial_surfaces_clubb2017_second_pass.param.in
    clubbsp_floodplain_mask.tif
    --read-path "${LSDTT_FLUVIAL_SURFACES_FIXTURE_DIR}"
    --stdout-contains "second pass (qq_upper_departure): HAND bound 0.310797 m (quantile 7710 of 10000 over 8354 first-pass floodplain cell(s))"
    --expect-raster-checksum clubbsp_floodplain_first_pass_mask.tif:8354
    --expect-raster-checksum clubbsp_floodplain_mask.tif:6472
    --expect-sha256 clubbsp_floodplain_hand_qq_second_pass.csv:f21f44785d1893a4a4df38845a88833ee2e39fc4b435daf2a8e657e0eb7b2475
    --expect-vector-feature-min clubbsp_floodplain_polygons.fgb:4
    --expect-vector-feature-min clubbsp_floodplain_first_pass_polygons.fgb:3
)
lsdtt_set_tests_properties(rigorous_fluvial_surfaces_clubb2017_second_pass PROPERTIES
    FIXTURES_REQUIRED fluvial_surfaces_fixture)

# reference_field_algorithm_choice relative_elevation on the staircase: the
# centrelines run through staircase_floor_mask.tif (channel + lowest tread) and
# the field reaches the treads above. HAND there is ~9 m on treads 2-4 alike
# (they drain sideways into small channels); the relative elevation puts them at
# 15.6 / 23.7 / 31.8 m (checked 2026-10-30). Pinned by the assignment line.
add_lsdtt_cli_smoke_test(
    rigorous_fluvial_surfaces_relative_elevation
    lsdtt-fluvial-surfaces
    fluvial_surfaces_relative_elevation.param.in
    fsrer_reference_field.tif
    --read-path "${LSDTT_FLUVIAL_SURFACES_FIXTURE_DIR}"
    --stdout-contains "--- Reference field: relative elevation ---"
    --stdout-contains "9 channel segment(s), 8 centreline(s)"
    --stdout-contains "8349 floor cell(s) by the larger-channel discs, 51651 by floor-width-scaled ground distance, 0 left out"
)
lsdtt_set_tests_properties(rigorous_fluvial_surfaces_relative_elevation PROPERTIES
    FIXTURES_REQUIRED fluvial_surfaces_fixture)

# lsdtt-valley-centreline relative elevation over the whole DEM (1500 m cap):
# centrelines through the staircase's valley floor, the field on the treads.
add_lsdtt_cli_smoke_test(
    rigorous_valley_centreline_relative_elevation_whole_dem
    lsdtt-valley-centreline
    valley_centreline_relative_elevation_whole_dem.param.in
    vcrer_relative_elevation.tif
    --read-path "${LSDTT_FLUVIAL_SURFACES_FIXTURE_DIR}"
    --stdout-contains "domain (whole_dem): 60000 pixel(s): 8349 on the floor by the larger-channel discs, 46780 by floor-width-scaled ground distance, 4871 reaching no centreline"
)
lsdtt_set_tests_properties(rigorous_valley_centreline_relative_elevation_whole_dem PROPERTIES
    FIXTURES_REQUIRED fluvial_surfaces_fixture)

# Terrace profile on the staircase: tread ids as patch labels, each tread
# profiled against the trunk centreline (median relative elevations 7.7 / 15.8
# / 23.8 / 31.8 m, checked 2026-10-30), binned every 200 m.
add_lsdtt_cli_smoke_test(
    rigorous_valley_centreline_terrace_profile
    lsdtt-valley-centreline
    valley_centreline_terrace_profile.param.in
    vctp_terrace_profile.csv
    --read-path "${LSDTT_FLUVIAL_SURFACES_FIXTURE_DIR}"
    --stdout-contains "4 patch(es) in 41 bin(s)"
    --expect-file-contains "vctp_terrace_patches.csv:4,,,0,9840,612,0.000000,2825.563492,129.699997,31.840004"
)
lsdtt_set_tests_properties(rigorous_valley_centreline_terrace_profile PROPERTIES
    FIXTURES_REQUIRED fluvial_surfaces_fixture)

# floodplain_bound_algorithm_choice: the staircase has 0.709 of its floodplain
# gradient below 0.1; with the switch at 0.7 auto takes hand_only, records it
# for a rerun, and gives exactly the explicit hand_only mask (4068 cells).
add_lsdtt_cli_smoke_test(
    rigorous_fluvial_surfaces_clubb2017_low_relief_auto
    lsdtt-fluvial-surfaces
    fluvial_surfaces_clubb2017_low_relief_auto.param.in
    lr_auto_floodplain_mask.tif
    --read-path "${LSDTT_FLUVIAL_SURFACES_FIXTURE_DIR}"
    --stdout-contains "low relief test: 0.709447 of the DEM's floodplain gradient is below 0.1 (switch above 0.7) -> hand_only"
    --stdout-contains "floodplain gradient bound : none (hand_only)"
    --expect-raster-checksum lr_auto_floodplain_mask.tif:4068
    --expect-file-contains "lr_auto_ingested_params.param:floodplain_bound_algorithm_choice: hand_only"
)
lsdtt_set_tests_properties(rigorous_fluvial_surfaces_clubb2017_low_relief_auto PROPERTIES
    FIXTURES_REQUIRED fluvial_surfaces_fixture)

add_lsdtt_cli_smoke_test(
    rigorous_fluvial_surfaces_clubb2017_low_relief_hand_only
    lsdtt-fluvial-surfaces
    fluvial_surfaces_clubb2017_low_relief_hand_only.param.in
    lr_hand_only_floodplain_mask.tif
    --read-path "${LSDTT_FLUVIAL_SURFACES_FIXTURE_DIR}"
    --stdout-contains "floodplain gradient bound : none (hand_only)"
    --expect-raster-checksum lr_hand_only_floodplain_mask.tif:4068
)
lsdtt_set_tests_properties(rigorous_fluvial_surfaces_clubb2017_low_relief_hand_only PROPERTIES
    FIXTURES_REQUIRED fluvial_surfaces_fixture)

# write_surface_patch_polygons: one polygon per labelled patch with its band.
# Four jenks levels on the staircase give 32 patches (8 per level), as many
# polygons as rows in _surface_patches.csv.
add_lsdtt_cli_smoke_test(
    rigorous_fluvial_surfaces_patch_polygons
    lsdtt-fluvial-surfaces
    fluvial_surfaces_patch_polygons.param.in
    fspp_surface_patch_polygons.fgb
    --read-path "${LSDTT_FLUVIAL_SURFACES_FIXTURE_DIR}"
    --expect-vector-feature-min fspp_surface_patch_polygons.fgb:32
)
lsdtt_set_tests_properties(rigorous_fluvial_surfaces_patch_polygons PROPERTIES
    FIXTURES_REQUIRED fluvial_surfaces_fixture)

# Several methods in one run: ladder,clubb2017 with the settings
# of rigorous_fluvial_surfaces_clubb2017. The clubb2017 class raster carries
# that test's pin (8718), i.e. running inside a list changes nothing; verified
# 2026-10-27 that every per-method raster and CSV is byte-identical to its
# single-method run. Comparison outputs pinned; class polygons 3 floodplain +
# 4 terrace areas.
add_lsdtt_cli_smoke_test(
    rigorous_fluvial_surfaces_methods_list
    lsdtt-fluvial-surfaces
    fluvial_surfaces_methods_list.param.in
    ml_method_comparison.csv
    --read-path "${LSDTT_FLUVIAL_SURFACES_FIXTURE_DIR}"
    --stdout-contains "method 2 of 2: clubb2017"
    --stdout-contains "bit k = method k in this order: 0=ladder 1=clubb2017"
    --expect ml_ladder_method_record.csv
    --expect ml_clubb2017_method_record.csv
    --expect-raster-checksum ml_clubb2017_surface_class.tif:8718
    --expect-raster-checksum ml_ladder_surface_class.tif:8735
    --expect-raster-checksum ml_method_agreement_floodplain.tif:25443
    --expect-raster-checksum ml_method_agreement_terrace.tif:364
    --expect-sha256 ml_method_comparison.csv:80188349e53db8c72137a0dfb895607809ef3f0e39a1d0ad4d4159bca5dbaa7a
    --expect-sha256 ml_method_summary.csv:1c5f0a5c749491d8da02c796adb3f3f620516a50133cbbc4265aed44c2458937
    --expect-vector-feature-min ml_clubb2017_surface_class_polygons.fgb:7
)
lsdtt_set_tests_properties(rigorous_fluvial_surfaces_methods_list PROPERTIES
    FIXTURES_REQUIRED fluvial_surfaces_fixture)

# fill_admitted_holes: the one enclosed gap a strict relief admission leaves
# on the linear-range surface is admitted (11966 -> 11967 cells).
add_lsdtt_cli_smoke_test(
    rigorous_fluvial_surfaces_fill_admitted_holes
    lsdtt-fluvial-surfaces
    fluvial_surfaces_fill_admitted_holes.param.in
    fah_band_index.tif
    --read-path "${LSDTT_FLUVIAL_SURFACES_FIXTURE_DIR}"
    --stdout-contains "fill_admitted_holes: 1 enclosed cell(s) admitted"
    --stdout-contains "11967 admitted cell(s)"
)
lsdtt_set_tests_properties(rigorous_fluvial_surfaces_fill_admitted_holes PROPERTIES
    FIXTURES_REQUIRED fluvial_surfaces_fixture)

# Training polygons on the staircase: fp / t1 / T5 drawn on treads
# 1-3. The trained breaks fall midway between the polygon cores, p25..p75
# (8.14 | 15.88 m and 16.40 | 23.88 m); the ladder's own search adds the break above them (k 4,
# the third break 24.88 m in the ladder windows CSV); under
# clubb2017 the trained bound replaces the QQ bound and the second pass is
# skipped. Ladder floodplain = the lowest level touching the channel: the six
# band-0 patches within 100 m of it; two band-0 patches 286 m away stay terrace.
# Patch level clustering finds the four treads (p at its 0.001 floor).
add_lsdtt_cli_smoke_test(
    rigorous_fluvial_surfaces_training
    lsdtt-fluvial-surfaces
    fluvial_surfaces_training.param.in
    trn_ladder_training_levels.tif
    --read-path "${LSDTT_FLUVIAL_SURFACES_FIXTURE_DIR}"
    --expect trn_ladder_training_fit.csv
    --expect trn_ladder_training_polygons.csv
    --expect trn_clubb2017_training_fit.csv
    --expect trn_clubb2017_training_levels.tif
    --stdout-contains "training: 2 of 2 window(s) banded with trained breaks"
    --stdout-contains "over all training: 12.01 20.14 m"
    --stdout-contains "floodplain / terrace HAND bound 12.01 m over all training (J 1)"
    --stdout-contains "the second pass is skipped"
    --stdout-contains "6 floodplain, 26 terrace surface(s)"
    --expect-file-contains "trn_ladder_ladder_windows.csv:0,500.000000,21458,4,resolved,0,12.009998"
    --expect-file-contains "trn_ladder_training_fit.csv:window,0,500,window,break_2,20.1399956,1,809,820,809,820,"
    --expect-file-contains "trn_ladder_patch_level_clusters.csv:all,-1,0,32,1.5,4,0.906580494,0.430937805,0.001,7.94000244 16.2199936 24.2399979 32.2499962,8 8 8 8"
    --stdout-contains "per window (>= 5 patches): 2 of 2 window(s) at p <= 0.05"
    --expect-file-contains "trn_ladder_training_polygons.csv:2,T5,terrace,5,1600,1600,24.1999969"
    --expect-file-contains "trn_ladder_training_polygons.csv:0.559997559,0,0,500,2,0"
    --expect-file-contains "trn_clubb2017_training_polygons.csv:0.519996643,0,0,500,1,0"
)
lsdtt_set_tests_properties(rigorous_fluvial_surfaces_training PROPERTIES
    FIXTURES_REQUIRED fluvial_surfaces_fixture)

# Mapped surfaces on the staircase, scored against its own tread map
# (tread 1 floodplain, treads 2-4 terrace, full valley length on both flanks).
# The burn equals the generator's tread map outside the three riser-gap ramps;
# FlatGeobuf, GeoPackage and EPSG:4326 copies burn byte-identically (one sha);
# a surveyed extent over the northern half scores 30000 of 60000 cells; an
# unknown class stops the run before any method. Scores are pinned as
# measured: the ladder's floodplain precision 1, recall 0.364 (the relief
# admission leaves the tread margins and gaps out); clubb2017 scores 0 here
# (its gradient bound admits almost none of the treads).
add_lsdtt_cli_smoke_test(
    rigorous_fluvial_surfaces_mapped
    lsdtt-fluvial-surfaces
    fluvial_surfaces_mapped.param.in
    map_mapped_surface_class.tif
    --read-path "${LSDTT_FLUVIAL_SURFACES_FIXTURE_DIR}"
    --expect map_mapped_surface_scores.csv
    --expect-sha256 map_mapped_surface_class.tif:a2addfd5238435bb7ecdb64295fed92b1f4dd8d70287ddcaab527b58cc29ff29
    --stdout-contains "2 floodplain, 6 terrace polygon(s)"
    --stdout-contains "the whole DEM counts as surveyed"
    --expect-file-contains "map_mapped_surface_scores.csv:ladder,60000,100,1,0.36362069,0.36362069,0.949250958,0.454111111,0.443347616,12400,0,0,6508,4218,874,19652,0,16348"
    --expect-file-contains "map_mapped_surface_scores.csv:clubb2017,60000,100,0,0,0,nan,0,0,12243,157,0,11600,0,0,36000,0,0"
)
lsdtt_set_tests_properties(rigorous_fluvial_surfaces_mapped PROPERTIES
    FIXTURES_REQUIRED fluvial_surfaces_fixture)

add_lsdtt_cli_smoke_test(
    rigorous_fluvial_surfaces_mapped_gpkg
    lsdtt-fluvial-surfaces
    fluvial_surfaces_mapped_gpkg.param.in
    mapgp_mapped_surface_class.tif
    --read-path "${LSDTT_FLUVIAL_SURFACES_FIXTURE_DIR}"
    --expect-sha256 mapgp_mapped_surface_class.tif:a2addfd5238435bb7ecdb64295fed92b1f4dd8d70287ddcaab527b58cc29ff29
)
lsdtt_set_tests_properties(rigorous_fluvial_surfaces_mapped_gpkg PROPERTIES
    FIXTURES_REQUIRED fluvial_surfaces_fixture)

add_lsdtt_cli_smoke_test(
    rigorous_fluvial_surfaces_mapped_4326
    lsdtt-fluvial-surfaces
    fluvial_surfaces_mapped_4326.param.in
    map4326_mapped_surface_class.tif
    --read-path "${LSDTT_FLUVIAL_SURFACES_FIXTURE_DIR}"
    --expect-sha256 map4326_mapped_surface_class.tif:a2addfd5238435bb7ecdb64295fed92b1f4dd8d70287ddcaab527b58cc29ff29
)
lsdtt_set_tests_properties(rigorous_fluvial_surfaces_mapped_4326 PROPERTIES
    FIXTURES_REQUIRED fluvial_surfaces_fixture)

add_lsdtt_cli_smoke_test(
    rigorous_fluvial_surfaces_mapped_extent
    lsdtt-fluvial-surfaces
    fluvial_surfaces_mapped_extent.param.in
    mapext_mapped_surface_scores.csv
    --read-path "${LSDTT_FLUVIAL_SURFACES_FIXTURE_DIR}"
    --stdout-contains "30000 cell(s) scored"
    --stdout-contains "ladder vs mapped (30000 cells)"
)
lsdtt_set_tests_properties(rigorous_fluvial_surfaces_mapped_extent PROPERTIES
    FIXTURES_REQUIRED fluvial_surfaces_fixture)

add_lsdtt_cli_smoke_test(
    rigorous_fluvial_surfaces_mapped_extent_multi
    lsdtt-fluvial-surfaces
    fluvial_surfaces_mapped_extent_multi.param.in
    mapextm_mapped_surface_scores.csv
    --read-path "${LSDTT_FLUVIAL_SURFACES_FIXTURE_DIR}"
    --stdout-contains "47600 cell(s) scored"
)
lsdtt_set_tests_properties(rigorous_fluvial_surfaces_mapped_extent_multi PROPERTIES
    FIXTURES_REQUIRED fluvial_surfaces_fixture)

add_lsdtt_cli_smoke_test(
    failure_fluvial_surfaces_mapped_bad_class
    lsdtt-fluvial-surfaces
    fluvial_surfaces_mapped_bad_class.param.in
    unused.tif
    --read-path "${LSDTT_FLUVIAL_SURFACES_FIXTURE_DIR}"
    --expect-failure
    --stderr-contains "'bench'"
)
lsdtt_set_tests_properties(failure_fluvial_surfaces_mapped_bad_class PROPERTIES
    FIXTURES_REQUIRED fluvial_surfaces_fixture)

add_lsdtt_cli_smoke_test(
    failure_fluvial_surfaces_training_bad_class
    lsdtt-fluvial-surfaces
    fluvial_surfaces_training_bad_class.param.in
    unused.tif
    --read-path "${LSDTT_FLUVIAL_SURFACES_FIXTURE_DIR}"
    --expect-failure
    --stderr-contains "training polygons: surface polygons: class field 'class' has values that are neither"
    --stderr-contains "'bench'"
)
lsdtt_set_tests_properties(failure_fluvial_surfaces_training_bad_class PROPERTIES
    FIXTURES_REQUIRED fluvial_surfaces_fixture)

add_lsdtt_cli_smoke_test(
    failure_fluvial_surfaces_training_bad_field
    lsdtt-fluvial-surfaces
    fluvial_surfaces_training_bad_field.param.in
    unused.tif
    --read-path "${LSDTT_FLUVIAL_SURFACES_FIXTURE_DIR}"
    --expect-failure
    --stderr-contains "field 'kind' not found"
)
lsdtt_set_tests_properties(failure_fluvial_surfaces_training_bad_field PROPERTIES
    FIXTURES_REQUIRED fluvial_surfaces_fixture)

add_lsdtt_cli_smoke_test(
    failure_fluvial_surfaces_methods_list_duplicate
    lsdtt-fluvial-surfaces
    fluvial_surfaces_methods_list_duplicate.param.in
    unused.tif
    --read-path "${LSDTT_FLUVIAL_SURFACES_FIXTURE_DIR}"
    --expect-failure
    --stderr-contains "lists 'ladder' twice"
)
lsdtt_set_tests_properties(failure_fluvial_surfaces_methods_list_duplicate PROPERTIES
    FIXTURES_REQUIRED fluvial_surfaces_fixture)

# The same, counting every channel pixel: the Strahler definition must reach
# the scope (9627 pixels at order >= 1 against 669 at order >= 3).
add_lsdtt_cli_smoke_test(
    rigorous_fluvial_surfaces_clubb2017_order1
    lsdtt-fluvial-surfaces
    fluvial_surfaces_clubb2017_order1.param.in
    clubb_band_index.tif
    --read-path "${LSDTT_FLUVIAL_SURFACES_FIXTURE_DIR}"
    --stdout-contains "channel = Strahler order >= 1: 9627 channel pixel(s)"
)
lsdtt_set_tests_properties(rigorous_fluvial_surfaces_clubb2017_order1 PROPERTIES
    FIXTURES_REQUIRED fluvial_surfaces_fixture)

# One stage overridden: the named stage wins and is reported as such; the rest
# still come from the preset.
add_lsdtt_cli_smoke_test(
    smoke_fluvial_surfaces_clubb2017_ladder_override
    lsdtt-fluvial-surfaces
    fluvial_surfaces_clubb2017_ladder.param.in
    clubb_band_index.tif
    --read-path "${LSDTT_FLUVIAL_SURFACES_FIXTURE_DIR}"
    --stdout-contains "ladder   (parameter file)"
    --stdout-contains "floodplain_thresholds   (preset clubb2017)"
    --stdout-contains "Ladder and band-then-label"
)
lsdtt_set_tests_properties(smoke_fluvial_surfaces_clubb2017_ladder_override PROPERTIES
    FIXTURES_REQUIRED fluvial_surfaces_fixture)

# A clubb2017-only parameter under the ladder preset is reported, not ignored.
add_lsdtt_cli_smoke_test(
    smoke_fluvial_surfaces_unused_parameter
    lsdtt-fluvial-surfaces
    fluvial_surfaces_unused_parameter.param.in
    unused_band_index.tif
    --read-path "${LSDTT_FLUVIAL_SURFACES_FIXTURE_DIR}"
    --stderr-contains "min_channel_stream_order is set but the resolved method stages do not use it"
)
lsdtt_set_tests_properties(smoke_fluvial_surfaces_unused_parameter PROPERTIES
    FIXTURES_REQUIRED fluvial_surfaces_fixture)

add_lsdtt_cli_smoke_test(
    failure_fluvial_surfaces_clubb2017_fitted_surface
    lsdtt-fluvial-surfaces
    fluvial_surfaces_clubb2017_fitted.param.in
    unused.tif
    --read-path "${LSDTT_FLUVIAL_SURFACES_FIXTURE_DIR}"
    --expect-failure
    --stderr-contains "needs reference_field_algorithm_choice hand"
)
lsdtt_set_tests_properties(failure_fluvial_surfaces_clubb2017_fitted_surface PROPERTIES
    FIXTURES_REQUIRED fluvial_surfaces_fixture)

add_lsdtt_cli_smoke_test(
    failure_fluvial_surfaces_clubb2017_negative_terrace_height
    lsdtt-fluvial-surfaces
    fluvial_surfaces_clubb2017_negative_terrace.param.in
    unused.tif
    --read-path "${LSDTT_FLUVIAL_SURFACES_FIXTURE_DIR}"
    --expect-failure
    --stderr-contains "min_terrace_height_above_floodplain_metres must be >= 0"
)
lsdtt_set_tests_properties(failure_fluvial_surfaces_clubb2017_negative_terrace_height PROPERTIES
    FIXTURES_REQUIRED fluvial_surfaces_fixture)

# ---- Fluvial-surface methods phase 5: cross-section levels ------------------
# The three-stage workflow on the staircase, checked against the known treads:
# levels within 1 m of 8/16/24/32 m, every window on the supplied levels, and
# each band at least 95 % the tread it should be.
if(LSDTT_ENABLE_PYTHON_GEO_TESTS)
    add_test(NAME rigorous_fluvial_surfaces_cross_section_levels_workflow
        COMMAND ${Python3_EXECUTABLE}
            "${CMAKE_CURRENT_SOURCE_DIR}/scripts/check_cross_section_levels_workflow.py"
            --surfaces-exe "$<TARGET_FILE:lsdtt-fluvial-surfaces>"
            --centreline-exe "$<TARGET_FILE:lsdtt-valley-centreline>"
            --param-dir "${LSDTT_TEST_PARAM_DIR}"
            --read-path "${LSDTT_FLUVIAL_SURFACES_FIXTURE_DIR}"
            --work-dir "${CMAKE_CURRENT_BINARY_DIR}/rigorous_fluvial_surfaces_cross_section_levels_workflow")
    lsdtt_mark_test_runtime(python rigorous_fluvial_surfaces_cross_section_levels_workflow)
else()
    lsdtt_mark_test_unregistered(python_geo rigorous_fluvial_surfaces_cross_section_levels_workflow)
endif()
lsdtt_set_tests_properties(rigorous_fluvial_surfaces_cross_section_levels_workflow PROPERTIES
    LABELS "cli"
    FIXTURES_REQUIRED fluvial_surfaces_fixture)

# segmentation_algorithm_choice whole_mask against within_bands on the staircase
# With closed risers whole_mask keeps
# every patch on one tread; with the find_risers centreline, open at the graded
# gaps, treads join into one patch - the documented failure, asserted.
if(LSDTT_ENABLE_PYTHON_GEO_TESTS)
    add_test(NAME rigorous_fluvial_surfaces_whole_mask_segmentation
        COMMAND ${Python3_EXECUTABLE}
            "${CMAKE_CURRENT_SOURCE_DIR}/scripts/check_whole_mask_segmentation.py"
            --surfaces-exe "$<TARGET_FILE:lsdtt-fluvial-surfaces>"
            --feature-finder-exe "$<TARGET_FILE:lsdtt-feature-finder>"
            --param-dir "${LSDTT_TEST_PARAM_DIR}"
            --read-path "${LSDTT_FLUVIAL_SURFACES_FIXTURE_DIR}"
            --work-dir "${CMAKE_CURRENT_BINARY_DIR}/rigorous_fluvial_surfaces_whole_mask_segmentation")
    lsdtt_mark_test_runtime(python rigorous_fluvial_surfaces_whole_mask_segmentation)
else()
    lsdtt_mark_test_unregistered(python_geo rigorous_fluvial_surfaces_whole_mask_segmentation)
endif()
lsdtt_set_tests_properties(rigorous_fluvial_surfaces_whole_mask_segmentation PROPERTIES
    LABELS "cli"
    FIXTURES_REQUIRED fluvial_surfaces_fixture)

# A levels file written by hand (the survey route, no valley-centreline run):
# source_id -1 applies to every channel. Pinned 2026-10-26; the same banding
# as the cross-section route gives on this fixture.
add_lsdtt_cli_smoke_test(
    rigorous_fluvial_surfaces_levels_hand
    lsdtt-fluvial-surfaces
    fluvial_surfaces_levels_hand.param.in
    xsh_band_index.tif
    --read-path "${LSDTT_FLUVIAL_SURFACES_FIXTURE_DIR}"
    --stdout-contains "supplied_levels   (preset cross_section_levels)"
    --stdout-contains "2 with supplied levels (0 below the resolvability floor)"
    --expect-raster-values "xsh_band_index.tif:0,1,2,3"
    --expect-raster-checksum xsh_band_index.tif:59958
)
lsdtt_set_tests_properties(rigorous_fluvial_surfaces_levels_hand PROPERTIES
    FIXTURES_REQUIRED fluvial_surfaces_fixture)

add_lsdtt_cli_smoke_test(
    failure_fluvial_surfaces_levels_wrong_source
    lsdtt-fluvial-surfaces
    fluvial_surfaces_levels_wrong_source.param.in
    unused.tif
    --read-path "${LSDTT_FLUVIAL_SURFACES_FIXTURE_DIR}"
    --expect-failure
    --stderr-contains "the supplied levels have no entry for channel source_id 0"
)
lsdtt_set_tests_properties(failure_fluvial_surfaces_levels_wrong_source PROPERTIES
    FIXTURES_REQUIRED fluvial_surfaces_fixture)

add_lsdtt_cli_smoke_test(
    failure_fluvial_surfaces_bad_reference_field
    lsdtt-fluvial-surfaces
    fluvial_surfaces_bad_reference_field.param.in
    unused.tif
    --read-path "${LSDTT_FLUVIAL_SURFACES_FIXTURE_DIR}"
    --expect-failure
    --stderr-contains "unknown reference_field_algorithm_choice 'elevation'. Valid: hand, relative_elevation, fitted_surface."
)
lsdtt_set_tests_properties(failure_fluvial_surfaces_bad_reference_field PROPERTIES
    FIXTURES_REQUIRED fluvial_surfaces_fixture)

# Phase 5e. The two-command pipeline, which is the only way either tool can be
# exercised end to end now: CLI 3 will not start without a mask, and CLI 2 is
# the only supported source of one. What this pins is the HANDOVER - that stage
# 2 says identification was supplied, which is the one line distinguishing "read
# the mask" from "quietly derived its own".
if(LSDTT_ENABLE_PYTHON_TESTS)
    add_test(NAME smoke_fluvial_surfaces_pipeline
        COMMAND ${Python3_EXECUTABLE}
            "${CMAKE_CURRENT_SOURCE_DIR}/scripts/run_fluvial_pipeline.py"
            --surfaces-exe "$<TARGET_FILE:lsdtt-fluvial-surfaces>"
            --centreline-exe "$<TARGET_FILE:lsdtt-valley-centreline>"
            --surfaces-param "${LSDTT_TEST_PARAM_DIR}/fluvial_pipeline_surfaces.param.in"
            --centreline-param "${LSDTT_TEST_PARAM_DIR}/valley_centreline_smoke.param.in"
            --read-path "${LSDTT_FLUVIAL_SURFACES_FIXTURE_DIR}"
            --work-dir "${CMAKE_CURRENT_BINARY_DIR}/smoke_fluvial_surfaces_pipeline"
            --expect vcin_valley_mask.tif
            --expect vc_valley_mask_channels.tif
            --expect vc_valley_centreline.fgb
            --stage2-stdout-contains "Valley identification is supplied")
    lsdtt_mark_test_runtime(python smoke_fluvial_surfaces_pipeline)
else()
    lsdtt_mark_test_unregistered(python smoke_fluvial_surfaces_pipeline)
endif()
lsdtt_set_tests_properties(smoke_fluvial_surfaces_pipeline PROPERTIES
    LABELS "cli"
    FIXTURES_REQUIRED fluvial_surfaces_fixture)

if(LSDTT_ENABLE_PYTHON_TESTS)
    add_test(NAME rigorous_valley_centreline_legacy_threshold
        COMMAND ${Python3_EXECUTABLE}
            "${CMAKE_CURRENT_SOURCE_DIR}/scripts/run_fluvial_pipeline.py"
            --surfaces-exe "$<TARGET_FILE:lsdtt-fluvial-surfaces>"
            --centreline-exe "$<TARGET_FILE:lsdtt-valley-centreline>"
            --surfaces-param "${LSDTT_TEST_PARAM_DIR}/fluvial_pipeline_surfaces.param.in"
            --centreline-param "${LSDTT_TEST_PARAM_DIR}/valley_centreline_legacy_threshold.param.in"
            --read-path "${LSDTT_FLUVIAL_SURFACES_FIXTURE_DIR}"
            --work-dir "${CMAKE_CURRENT_BINARY_DIR}/rigorous_valley_centreline_legacy_threshold"
            --expect vc_legacy_valley_mask_channels.tif
            --expect vc_legacy_valley_centreline.fgb
            --stage2-stdout-contains "Channel threshold: 200 contributing pixels"
            --stage2-stderr-contains "'channel_threshold_contributing_pixels' is deprecated"
            --stage2-stderr-contains "use 'threshold_contributing_pixels'")
    lsdtt_mark_test_runtime(python rigorous_valley_centreline_legacy_threshold)
else()
    lsdtt_mark_test_unregistered(python rigorous_valley_centreline_legacy_threshold)
endif()
lsdtt_set_tests_properties(rigorous_valley_centreline_legacy_threshold PROPERTIES
    LABELS "cli"
    FIXTURES_REQUIRED fluvial_surfaces_fixture)

# The point-file key renamed 2026-10-26 (source_points_path ->
# source_points_fname, the same name and read_path resolution as
# lsdtt-fluvial-surfaces): an old file must still run, with a warning.
if(LSDTT_ENABLE_PYTHON_TESTS)
    add_test(NAME rigorous_valley_centreline_legacy_points_key
        COMMAND ${Python3_EXECUTABLE}
            "${CMAKE_CURRENT_SOURCE_DIR}/scripts/run_fluvial_pipeline.py"
            --surfaces-exe "$<TARGET_FILE:lsdtt-fluvial-surfaces>"
            --centreline-exe "$<TARGET_FILE:lsdtt-valley-centreline>"
            --surfaces-param "${LSDTT_TEST_PARAM_DIR}/fluvial_pipeline_surfaces.param.in"
            --centreline-param "${LSDTT_TEST_PARAM_DIR}/valley_centreline_legacy_points_key.param.in"
            --read-path "${LSDTT_FLUVIAL_SURFACES_FIXTURE_DIR}"
            --work-dir "${CMAKE_CURRENT_BINARY_DIR}/rigorous_valley_centreline_legacy_points_key"
            --expect vcold_valley_centreline.fgb
            --stage2-stderr-contains "'source_points_path' is deprecated"
            --stage2-stderr-contains "use 'source_points_fname'")
    lsdtt_mark_test_runtime(python rigorous_valley_centreline_legacy_points_key)
else()
    lsdtt_mark_test_unregistered(python rigorous_valley_centreline_legacy_points_key)
endif()
lsdtt_set_tests_properties(rigorous_valley_centreline_legacy_points_key PROPERTIES
    LABELS "cli"
    FIXTURES_REQUIRED fluvial_surfaces_fixture)

# The valley is now a required input, and the message must name where it comes
# from - a user with an old parameter file needs to be told the tool, not just
# that something is missing.
add_lsdtt_cli_smoke_test(
    failure_valley_centreline_no_mask
    lsdtt-valley-centreline
    valley_centreline_no_mask.param.in
    unused.tif
    --read-path "${LSDTT_FLUVIAL_SURFACES_FIXTURE_DIR}"
    --expect-failure
    --stderr-contains "valley_mask_fname is required"
    --stderr-contains "lsdtt-fluvial-surfaces"
)
lsdtt_set_tests_properties(failure_valley_centreline_no_mask PROPERTIES
    FIXTURES_REQUIRED fluvial_surfaces_fixture)

# rigorous_valley_centreline_break_pin was RETIRED at phase 5d, with the code it
# pinned. It asserted the breaks lsdtt-valley-centreline derived for itself, and
# that tool no longer derives any - identification moved to
# lsdtt-fluvial-surfaces. Deleting a test along with the path it covers is the
# intended outcome, not a lost guard.
#
# The C5 concern it existed for - that find_population_breaks is shared and must
# not move underneath the ladder - is now carried by T5 in
# unit_fluvial_surface_ladder, which calls that function directly and asserts the
# ladder's k = 2 reproduces its first break bit-for-bit. T5 must be green in the
# same run for this retirement to be safe.

# The phase 1 baseline. Band-then-label rewrites exactly this raster: today the
# labelling runs once over the whole broken flat mask, and after the reordering
# it runs once per band. Without this checksum recorded BEFORE the change there
# is nothing to say what the reordering actually did.
add_lsdtt_cli_smoke_test(
    rigorous_flat_surface_finder_label_baseline
    lsdtt-flat-surface-finder
    flat_surface_finder_smoke.param.in
    fsf_labelled_flat_mask.tif
    --read-path "${LSDTT_FLUVIAL_SURFACES_FIXTURE_DIR}"
    --expect-raster-checksum fsf_labelled_flat_mask.tif:20987
    --expect-raster-checksum fsf_combined_flat_mask.tif:13879
)
lsdtt_set_tests_properties(rigorous_flat_surface_finder_label_baseline PROPERTIES
    FIXTURES_REQUIRED fluvial_surfaces_fixture)

# The CDF-divergence relief bound in this tool, pinned 2026-10-26 before its
# code moved into the shared threshold-selection utility.
add_lsdtt_cli_smoke_test(
    rigorous_flat_surface_finder_cdf_divergence
    lsdtt-flat-surface-finder
    flat_surface_finder_cdf_divergence.param.in
    fsf_combined_flat_mask.tif
    --read-path "${LSDTT_FLUVIAL_SURFACES_FIXTURE_DIR}"
    --stdout-contains "Final threshold   : 8.56527 m"
    --expect-raster-checksum fsf_combined_flat_mask.tif:18337
    --expect-sha256 fsf_relief_adaptive_cdf.csv:22a90923f965f5084afb8a07772ac094e910554aad93425fce2155b434b6817c
)
lsdtt_set_tests_properties(rigorous_flat_surface_finder_cdf_divergence PROPERTIES
    FIXTURES_REQUIRED fluvial_surfaces_fixture)

# The relief bound from the Clubb et al. (2017) QQ rule, in both tools that
# admit flat ground by relief. Same DEM, same 5 x 5 relief window: the two QQ
# tables carry the same hash, because the name, the aim and the code are shared.
add_lsdtt_cli_smoke_test(
    rigorous_flat_surface_finder_relief_qq
    lsdtt-flat-surface-finder
    flat_surface_finder_relief_qq.param.in
    fsf_combined_flat_mask.tif
    --read-path "${LSDTT_FLUVIAL_SURFACES_FIXTURE_DIR}"
    --stdout-contains "QQ quantile       : 2400 of 10000"
    --stdout-contains "Final threshold   : 0.420538 m"
    --expect-raster-checksum fsf_combined_flat_mask.tif:9602
    --expect-sha256 fsf_relief_qq.csv:a74a5824fd5f2c3f10c6cf442cdbbd95dee0b8a109e088fba522a822a850590f
)
lsdtt_set_tests_properties(rigorous_flat_surface_finder_relief_qq PROPERTIES
    FIXTURES_REQUIRED fluvial_surfaces_fixture)

add_lsdtt_cli_smoke_test(
    rigorous_mountain_front_relief_qq
    lsdtt-mountain-front
    mountain_front_relief_qq.param.in
    mf_relief.tif
    --read-path "${LSDTT_FLUVIAL_SURFACES_FIXTURE_DIR}"
    --stdout-contains "QQ quantile 2400 of 10000: threshold 0.420538 m"
    --expect-raster-checksum mf_mountain_mask.tif:30398
    --expect-sha256 mf_relief_qq.csv:a74a5824fd5f2c3f10c6cf442cdbbd95dee0b8a109e088fba522a822a850590f
)
lsdtt_set_tests_properties(rigorous_mountain_front_relief_qq PROPERTIES
    FIXTURES_REQUIRED fluvial_surfaces_fixture)

# The durable guard on parameter naming. Walks every built lsdtt-* CLI, reads its declared
# parameters from the help registry (`<cli> -h json`) and applies the mechanical
# naming rules -- no `_m`/`_px`/`_deg` suffixes, no `print_*` that
# writes a file, bool raster outputs ending `_raster`/`_rasters`, one raster and
# one vector format parameter per CLI, no `minimum_` (outside the documented
# `<subject>_minimum_contributing_pixels` family) and no `tol`.
#
# Report-only against a baseline allow-list, so it PASSES on today's tree and
# fails only on a NEWLY introduced violation. The allow-list must only ever
# shrink; a stale entry is reported, not failed, so no later phase is blocked by
# it. See the script's header comment.
#
# No build dependency beyond the CLIs themselves: it runs the binaries that are
# already there. The binary directory comes from a generator expression rather
# than a hard-coded path.
if(LSDTT_ENABLE_PYTHON_TESTS)
    add_test(NAME unit_cli_parameter_conventions
        COMMAND ${Python3_EXECUTABLE}
            "${CMAKE_CURRENT_SOURCE_DIR}/scripts/check_cli_parameter_conventions.py"
            --bin-dir "$<TARGET_FILE_DIR:lsdtt-channel-extraction>"
            --cli-src-dir "${CMAKE_SOURCE_DIR}/cli"
            --baseline "${LSDTT_TEST_DATA_DIR}/cli_parameter_conventions_baseline.json")
    lsdtt_mark_test_runtime(python unit_cli_parameter_conventions)
else()
    lsdtt_mark_test_unregistered(python unit_cli_parameter_conventions)
endif()
lsdtt_set_tests_properties(unit_cli_parameter_conventions PROPERTIES
    LABELS "cli;conventions"
    TIMEOUT 900)

# Phase 6: generated Markdown parameter sections must exactly match the ordered
# JSON help registries for every CLI.
if(LSDTT_ENABLE_PYTHON_TESTS)
    add_test(NAME unit_cli_reference_is_current
        COMMAND ${Python3_EXECUTABLE}
            "${CMAKE_SOURCE_DIR}/tools/generate_cli_reference.py"
            --check
            --bin-dir "$<TARGET_FILE_DIR:lsdtt-channel-extraction>"
            --repo-root "${CMAKE_SOURCE_DIR}")
    lsdtt_mark_test_runtime(python unit_cli_reference_is_current)
else()
    lsdtt_mark_test_unregistered(python unit_cli_reference_is_current)
endif()
lsdtt_set_tests_properties(unit_cli_reference_is_current PROPERTIES
    LABELS "cli;docs"
    TIMEOUT 900)

# The pip wheel declares its CLIs in four hand-kept places, which must match the
# root lsdtt_cli() targets. python_bindings/tests/test_cli_packaging.py checks
# this too but needs pytest, so it never ran here, and four CLIs went missing
# from the wheel unnoticed. Reads source files only; needs no build.
if(LSDTT_ENABLE_PYTHON_TESTS)
    add_test(NAME unit_wheel_cli_inventory
        COMMAND ${Python3_EXECUTABLE}
            "${CMAKE_CURRENT_SOURCE_DIR}/scripts/check_wheel_cli_inventory.py"
            --repo-root "${CMAKE_SOURCE_DIR}")
    lsdtt_mark_test_runtime(python unit_wheel_cli_inventory)
else()
    lsdtt_mark_test_unregistered(python unit_wheel_cli_inventory)
endif()
lsdtt_set_tests_properties(unit_wheel_cli_inventory PROPERTIES
    LABELS "cli;packaging"
    TIMEOUT 60)

# One version number for the whole repository: project(lsdtt3 VERSION ...) is
# the truth, and pyproject.toml, the QGIS metadata.txt, the conda recipe,
# pixi.toml and debian/changelog carry checked copies (tools/set_version.py
# rewrites them together). Also checks that the PyPI README explains the
# versions that have no wheel. Reads source files only; needs no build.
if(LSDTT_ENABLE_PYTHON_TESTS)
    add_test(NAME unit_version_consistency
        COMMAND ${Python3_EXECUTABLE}
            "${CMAKE_CURRENT_SOURCE_DIR}/scripts/check_version_consistency.py"
            --repo-root "${CMAKE_SOURCE_DIR}")
    lsdtt_mark_test_runtime(python unit_version_consistency)
else()
    lsdtt_mark_test_unregistered(python unit_version_consistency)
endif()
lsdtt_set_tests_properties(unit_version_consistency PROPERTIES
    LABELS "packaging"
    TIMEOUT 60)

# Builds the Python bindings once (pip, CLIs off) into a venv and runs the
# whole python_bindings/tests suite against them: the cheap Linux check between
# manual runs of the full wheel matrix. tools/make_release.sh runs the same
# script on every release. Compiles all of src/ and needs network access (or a
# pip cache) for nanobind, numpy, geopandas and friends, hence heavy and
# serial. The interpreter needs the venv module (Debian/Ubuntu: python3-venv);
# without it the test is skipped. Point LSDTT_PYTHON_SMOKE_PYTHON at another
# interpreter to use that one instead.
if(LSDTT_ENABLE_PYTHON_TESTS)
    set(LSDTT_PYTHON_SMOKE_PYTHON "${Python3_EXECUTABLE}" CACHE FILEPATH
        "Interpreter for the python_bindings_smoke test (needs the venv module)")
    add_test(NAME python_bindings_smoke
        COMMAND bash
            "${CMAKE_CURRENT_SOURCE_DIR}/scripts/check_python_bindings_smoke.sh"
            --repo-root "${CMAKE_SOURCE_DIR}"
            --work-dir "${CMAKE_BINARY_DIR}/python_smoke"
            --python "${LSDTT_PYTHON_SMOKE_PYTHON}"
            --reuse)
    lsdtt_mark_test_runtime(python python_bindings_smoke)
else()
    lsdtt_mark_test_unregistered(python python_bindings_smoke)
endif()
lsdtt_set_tests_properties(python_bindings_smoke PROPERTIES
    LABELS "python;packaging"
    SKIP_RETURN_CODE 77
    TIMEOUT 3600)

# =============================================================================
# COMPARTMENTALISATION: domain and cost labels
#
# Running all 100+ tests after every change is wasteful: most of them cannot be
# affected by a given edit, and a handful are computationally heavy enough to
# dominate the wall clock (the implicit hillslope solver saturates every core
# for ~9 minutes on its own; the real-world cosmogenic tiers take longer still).
#
# So every test carries a DOMAIN label saying which part of the codebase it
# exercises, and the expensive ones carry a COST label so they can be excluded.
# The intended workflow:
#
#   ctest -L cosmo -LE heavy      # I changed cosmogenic code
#   ctest -L network -LE heavy    # I changed flow routing / chi / basins
#   ctest -LE heavy               # everything cheap; the normal inner loop
#   ctest                         # everything; before a commit
#
# The domain is derived from the TEST NAME rather than declared at each of the
# hundred-odd call sites, so there is one table to read and one place to change.
# A test whose name matches no rule is a CONFIGURE ERROR, not a silently
# unlabelled test -- that is what stops the scheme rotting as tests are added.
#
# DOMAIN RULES, first match wins. Add a rule here when you add a subsystem.
# =============================================================================

set(LSDTT_DOMAIN_RULES
    # Cross-package benchmarks, registered only under LSDTT_ENABLE_BENCHMARKS.
    # FIRST, because first match wins and this must claim every benchmark_*
    # test regardless of what else its name contains.
    #
    # Without this rule the guard below rejects EVERY benchmark test and
    # `cmake -DLSDTT_ENABLE_BENCHMARKS=ON` fails to configure at all -- which
    # is what it did between the guard landing and 2026-09-08, because the
    # guard was added while benchmarks were OFF and nobody reconfigured with
    # them ON. The default-OFF option hid it.
    #
    # The subsystem labels these tests set for themselves ("benchmark;flow",
    # "benchmark;surface_metrics", "benchmark;chi") are MERGED with this
    # domain further down, so `ctest -L flow` still reaches the flow benchmark.
    "^benchmark_"            "benchmark"
    "synthetic_(surfaces|networks|cli)|_ocn_" "synthetic_networks"
    "niemi"                  "cosmo"
    "quartz|fertility"       "quartz"
    "cosmo"                  "cosmo"
    "landslide"              "landslide"
    "step_change"            "step_change"
    "landscape_equilibration|lem_equilibration" "landscape_equilibration"
    "streamlines"            "streamlines"
    "hilltop_flow_routing|hilltop_coupling|hillslope_channel_coupling" "hilltop_flow_routing"
    "paired_flowpaths"       "paired_flowpaths"
    # The range_defined_ridgelines path (C2/C3) is its own subsystem: boundary
    # tracing and BC-aware surface fitting, not the paired-flowpath work the
    # other ridgeline tests cover.
    "ridgeline"              "ridgeline"
    "divide_diagnostics"     "divide_diagnostics"
    "divide_migration"       "divide_migration"
    "fetch"                  "fetch"
    "hydraulic"              "hydraulic"
    "giuh"                   "giuh"
    "surface_metrics"        "surface_metrics"
    "spectral_analysis"      "surface_metrics"
    "hillslope"              "hillslope"
    "lem"                    "lem"
    "advection_diffusion"    "lem"
    # Ahead of the "channel" rule below on purpose: CLI 2 of the fluvial-surface
    # split writes a channel network, so tests named for that would otherwise be
    # claimed by "network" and disappear from a fluvial_surfaces compartment run.
    "fluvial_surface|flat_surface|valley_centreline|valley_width|valley_relative|perimeter_sources" "fluvial_surfaces"
    # The shared Topographic Wetness Index,
    # pointed at "network" ON PURPOSE rather than given its own domain: the
    # existing TWI/STI tests already live in the network compartment, and
    # splitting them off now would quietly remove them from `ctest -L network`
    # for anyone working on drainage. When lsdtt-distributed-hydrology lands and
    # there is a real cluster of tests on this module, flip this one line to
    # "wetness_index" and move the STI tests with it, deliberately.
    "wetness_index|hydro_velocity|hydro_roughness|hydro_routing|topmodel|canopy|channel_extent|infiltration|landcover|distributed_hydrology|hype"          "network"
    "chi|channel|drainage|basin|flow_properties|flow_snapping|flow_d8|link_length|flat_resolution"  "network"
    # performance_flow_routines (LSDTT_ENABLE_PERF_TESTS) matched no rule, so
    # configuring with that option ON failed until 2026-09-18 -- the same trap
    # as the benchmark rule at the top, hidden by the option's OFF default.
    "flow_routines"          "network"
    # Depression breaching is flow preprocessing. Ahead of the raster rule
    # because "fill" would otherwise claim breach tests that mention filling.
    "breach"                 "network"
    "raster|nodata|fill|tif|coastal|dem_grouping|vertical_units|multipart" "raster"
    "param|secret|cli_reference|wheel_cli_inventory|version_consistency|python_bindings_smoke|test_env|configure_contract|run_record" "core"
    "depression_hierarchy"   "network"
    "point_input"            "core"
    # Phase 2 of the CLI consistency remediation: the shared helper homes
    # extracted out of cli/ (to_lower_ascii, output_extension, csv_number,
    # planar_distance, quote-aware csv_split_row). They are library-level
    # string/format hygiene, so they belong in the core compartment.
    "cli_helper|csv_table"   "core"
    # Added 2026-09-10 with the first tests these two CLIs have ever had.
    # mountain_front is landform/relief work and feature_finder is multi-scale
    # surface response, so both sit in the surface compartment.
    "mountain_front|feature_finder|swath_profile" "surface_metrics"
    # Added 2026-10-07 with the first tests these four CLIs have had. Cross
    # sections, bankfull width and the water-body centreline are channel /
    # hydrography work; marsh platform and tidal channel extraction is landform
    # extraction, beside feature_finder.
    "cross_sections|bankfull_width|water_centreline|water_line_repair|slope_area|outlet_matching" "network"
    "marsh_topography"       "surface_metrics"
    "basin_stats"            "network"
    "junction_graph"         "network"
    "grid_neighbours"        "raster"
    "vector_layer_transform" "gdal"
    "connected_components|grid_geodesic|oblique_frame|true_north" "raster"
    "boundary_conditions|output_file|splitmix|thread_count" "core"
    "stream_burning"         "network"
    "sediment_routing"       "network"
    "linear_stats"           "network"
    "basin_surface_metrics"  "surface_metrics"
    # Added 2026-09-14: the native image-filter primitives and raster image
    # operations (Otsu, adaptive threshold, Canny, Hough) that replaced OpenCV.
    "image_filters|image_ops" "image_ops"
    # Added 2026-09-16 with the GDAL/PROJ friction work. Its own compartment
    # rather than "core" because it is the geospatial runtime boundary: PROJ/GDAL
    # data pinning (Phase 0), the header facade (Phase 1) and the format contract
    # (Phase 2 chunk 2a). Anyone touching that boundary wants `ctest -L gdal`.
    "gdal_"                  "gdal"
    # Library-level statistics (src/utils/lsdtt_statistics and relatives).
    "order_statistics"       "core"
)

# COST: tests that are expensive enough to keep out of the inner loop. The
# figure is wall-clock when the test has the machine to itself.
#
#   unit_hillslope_implicit    ~9 min (2026-08-04, 8-core machine), saturates
#                              all 8 cores; 586 s on 2026-09-18 (8-core laptop)
#   cosmo_realworld_*          67 s to 30 min, and most need the network
#
# smoke_lem_roering_implicit_strip was listed here (~5 min on 2026-08-04). It
# ran in 2 s on 2026-09-18 with its steady-state assertions passing, and
# CTestCostData.txt records 0.1 s, so it came off this list and the serial list
# below.
set(LSDTT_HEAVY_TESTS
    unit_hillslope_implicit
    cosmo_realworld_clapp2002
    cosmo_realworld_clapp2002_lsdn
    cosmo_realworld_dethier2014
    cosmo_realworld_dethier2014_lsdn
    cosmo_realworld_dibiase2010
    cosmo_realworld_dibiase2010_lsdn
    cosmo_realworld_lsdn_pixel_vs_lumped
    cosmo_realworld_production_isolation
    python_bindings_smoke
)

# Tests that use every core they can get. Under `ctest -j` these must not share
# the machine: doing so does not merely slow them, it starves them past their
# timeout, which is how unit_hillslope_implicit came to "fail" at 1500 s while
# taking 9 minutes when run alone.
set(LSDTT_SERIAL_TESTS
    unit_hillslope_implicit
    python_bindings_smoke
)

get_property(_lsdtt_tests DIRECTORY PROPERTY TESTS)
set(_lsdtt_unclassified "")

foreach(_test IN LISTS _lsdtt_tests)
    set(_domains "")
    list(LENGTH LSDTT_DOMAIN_RULES _n_rules)
    math(EXPR _last "${_n_rules} - 1")
    foreach(_i RANGE 0 ${_last} 2)
        if(_domains STREQUAL "")
            list(GET LSDTT_DOMAIN_RULES ${_i} _pattern)
            math(EXPR _j "${_i} + 1")
            list(GET LSDTT_DOMAIN_RULES ${_j} _domain)
            if(_test MATCHES "${_pattern}")
                set(_domains "${_domain}")
            endif()
        endif()
    endforeach()

    if(_domains STREQUAL "")
        list(APPEND _lsdtt_unclassified "${_test}")
        continue()
    endif()

    # The implicit solvers are hillslope physics wherever they are driven from,
    # so an LEM test that runs one belongs to both compartments.
    if(_test MATCHES "roering|implicit" AND NOT _domains STREQUAL "hillslope")
        list(APPEND _domains "hillslope")
    endif()

    # lsdtt-cosmo-landslides is BOTH a cosmogenic tool and the landslide
    # module's CLI. "cosmo" precedes "landslide" in the rules above and first
    # match wins, so without this its tests would silently leave the landslide
    # compartment at the 2026-08-05 rename and `ctest -L landslide` would stop
    # seeing the CLI for the module it is named after.
    if(_test MATCHES "landslide" AND NOT _domains STREQUAL "landslide")
        list(APPEND _domains "landslide")
    endif()

    # Lithology cuts across compartments: the map readers live with the LEM,
    # the Macrostrat fetch with fetch, and the analyses by lithologic unit with
    # chi and network. `ctest -L lithology` reaches all of them.
    if(_test MATCHES "lithology")
        list(APPEND _domains "lithology")
    endif()

    # The erodibility calibration from basin erosion rates spans chi (the
    # calibration maths and the sample reader) and fetch (the OCTOPUS rates).
    # `ctest -L erodibility_calibration` runs the whole thread.
    if(_test MATCHES "erodibility_calibration|erosion_rates")
        list(APPEND _domains "erodibility_calibration")
    endif()

    if("${_test}" IN_LIST LSDTT_HEAVY_TESTS)
        list(APPEND _domains "heavy")
    endif()

    get_test_property(${_test} LABELS _existing)
    if(_existing)
        list(APPEND _domains ${_existing})
    endif()
    list(REMOVE_DUPLICATES _domains)
    set_tests_properties(${_test} PROPERTIES LABELS "${_domains}")

    if("${_test}" IN_LIST LSDTT_SERIAL_TESTS)
        set_tests_properties(${_test} PROPERTIES RUN_SERIAL TRUE)
        get_test_property(${_test} TIMEOUT _has_timeout)
        if(NOT _has_timeout)
            set_tests_properties(${_test} PROPERTIES TIMEOUT 3600)
        endif()
    endif()
endforeach()

if(NOT _lsdtt_unclassified STREQUAL "")
    string(REPLACE ";" "\n    " _pretty "${_lsdtt_unclassified}")
    message(FATAL_ERROR
        "These tests match no domain rule and would be unreachable by a "
        "compartment run:\n    ${_pretty}\n"
        "Add a rule to LSDTT_DOMAIN_RULES above, or rename the test so it "
        "carries its subsystem in the name. Do not leave a test unlabelled -- "
        "it will silently stop being run by anyone working on that subsystem.")
endif()

# =============================================================================
# BUILD CURRENCY GUARD
#
# A build system will happily run yesterday's binary, and a stale one produces a
# confident green. That is not hypothetical here: the parameter-convention
# baseline was once generated against a tree with ten stale binaries and one
# missing, and a smoke test passed against a CLI built before a signature change
# while a full build was failing on that very file.
#
# The CLI smoke harness checks its own executable (assert_binary_is_current in
# run_cli_smoke.py). That leaves the ~60 C++ test drivers, which nothing checks:
# a stale driver runs the old object code and reports success.
#
# So: a manifest of every executable target and its sources, written at generate
# time so $<TARGET_FILE:...> resolves, and one test that compares mtimes. Built
# from CMake's own target list rather than a hand-written one, so a driver added
# later is covered without anyone remembering to add it.
# =============================================================================
set(_currency_lines "")
foreach(_dir "${CMAKE_SOURCE_DIR}" "${CMAKE_CURRENT_SOURCE_DIR}")
    get_property(_targets DIRECTORY "${_dir}" PROPERTY BUILDSYSTEM_TARGETS)
    foreach(_tgt ${_targets})
        get_target_property(_type ${_tgt} TYPE)
        if(NOT _type STREQUAL "EXECUTABLE")
            continue()
        endif()
        get_target_property(_srcs ${_tgt} SOURCES)
        get_target_property(_tgt_src_dir ${_tgt} SOURCE_DIR)
        set(_abs_srcs "")
        foreach(_s ${_srcs})
            if(IS_ABSOLUTE "${_s}")
                list(APPEND _abs_srcs "${_s}")
            else()
                list(APPEND _abs_srcs "${_tgt_src_dir}/${_s}")
            endif()
        endforeach()
        string(REPLACE ";" "|" _joined "${_abs_srcs}")
        list(APPEND _currency_lines "${_tgt}\t$<TARGET_FILE:${_tgt}>\t${_joined}")
    endforeach()
endforeach()
list(JOIN _currency_lines "\n" _currency_body)

file(GENERATE
     OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/build_currency_manifest.tsv"
     CONTENT "${_currency_body}\n")

add_test(NAME unit_build_is_current
    COMMAND ${CMAKE_COMMAND}
        -DMANIFEST=${CMAKE_CURRENT_BINARY_DIR}/build_currency_manifest.tsv
        -P "${CMAKE_CURRENT_SOURCE_DIR}/cmake/check_build_is_current.cmake")
set_tests_properties(unit_build_is_current PROPERTIES
    LABELS "core;conventions"
    TIMEOUT 300)

# =============================================================================
# TEST INVENTORY AUDIT
#
# Checks every registered test's tier / flag / component labels and its
# environment and timeout against the command it actually runs, and that every
# test script launches processes through lsdtt_test_env.py. The rules are in
# cmake/audit_test_inventory.cmake. It reads the generated CTestTestfile.cmake
# rather than calling `ctest --show-only`, because any ctest run inside a
# running suite empties that suite's LastTest.log.
#
# The self-test runs the same audit over a canned test file with one planted
# violation per check (cmake/audit_selftest/) and passes only if exactly those
# are reported, so a broken audit cannot pass silently.
# =============================================================================
add_test(NAME unit_test_inventory_audit
    COMMAND ${CMAKE_COMMAND}
        -DTESTFILE=${CMAKE_CURRENT_BINARY_DIR}/CTestTestfile.cmake
        -DSCRIPTS_DIR=${CMAKE_CURRENT_SOURCE_DIR}/scripts
        -DCTEST_CONFIGURATION_TYPE=$<CONFIG>
        -P "${CMAKE_CURRENT_SOURCE_DIR}/cmake/audit_test_inventory.cmake")
set_tests_properties(unit_test_inventory_audit PROPERTIES
    LABELS "core;conventions"
    TIMEOUT 120)

# The test commands the documentation tells people to type: every -DLSDTT_*
# option it names exists, every anchored -L/-LE label is one this file can
# assign, and the test_report target and run_tests.cmake it names exist
# (cmake/check_doc_test_commands.cmake). Reads source files only.
add_test(NAME unit_docs_test_commands
    COMMAND ${CMAKE_COMMAND} -DSOURCE_DIR=${CMAKE_SOURCE_DIR}
        -P "${CMAKE_CURRENT_SOURCE_DIR}/cmake/check_doc_test_commands.cmake")
set_tests_properties(unit_docs_test_commands PROPERTIES
    LABELS "core;docs"
    TIMEOUT 60)

add_test(NAME unit_test_inventory_audit_selftest
    COMMAND ${CMAKE_COMMAND}
        -DSELFTEST_DIR=${CMAKE_CURRENT_SOURCE_DIR}/cmake/audit_selftest
        -P "${CMAKE_CURRENT_SOURCE_DIR}/cmake/audit_test_inventory.cmake")
set_tests_properties(unit_test_inventory_audit_selftest PROPERTIES
    LABELS "core;conventions"
    TIMEOUT 60)

# The evidence and level pytest markers of python_bindings/tests, the Python
# counterpart of evidence_labels.cmake. The checker reads each test file with
# `ast`, so it needs neither pytest nor a built binding: a test added without
# markers is caught here, in every Python-tier configure, rather than first in
# a release wheel build, where python_bindings/tests/conftest.py would stop
# collection. It also checks that its marker list equals LSDTT_EVIDENCE_LABELS
# and LSDTT_LEVEL_LABELS. The self-test runs it on built-in snippets with known
# verdicts, so a checker that accepts everything cannot pass.
if(LSDTT_ENABLE_PYTHON_TESTS)
    add_test(NAME unit_python_test_markers
        COMMAND ${Python3_EXECUTABLE}
            "${CMAKE_CURRENT_SOURCE_DIR}/scripts/check_python_test_markers.py"
            --repo-root "${CMAKE_SOURCE_DIR}"
            --cmake-labels "${CMAKE_CURRENT_SOURCE_DIR}/cmake/evidence_labels.cmake")
    lsdtt_mark_test_runtime(python unit_python_test_markers)
    add_test(NAME unit_python_test_markers_selftest
        COMMAND ${Python3_EXECUTABLE}
            "${CMAKE_CURRENT_SOURCE_DIR}/scripts/check_python_test_markers.py"
            --self-test)
    lsdtt_mark_test_runtime(python unit_python_test_markers_selftest)
else()
    lsdtt_mark_test_unregistered(python unit_python_test_markers)
    lsdtt_mark_test_unregistered(python unit_python_test_markers_selftest)
endif()
lsdtt_set_tests_properties(unit_python_test_markers unit_python_test_markers_selftest
    PROPERTIES
    LABELS "core;conventions;python"
    TIMEOUT 60)

# The python_bindings/tests files that never import lsdtt3 (API examples,
# naming, reference and binding-source contracts, CLI packaging), run with
# pytest WITHOUT a built binding, so they are checked in every geo-tier
# configure rather than only when the binding is built. pytest is part of the
# geo tier (the configure import check above, pixi feature.test).
set(LSDTT_PYTHON_BINDING_FREE_TEST_FILES
    test_api_examples.py test_api_naming.py test_api_reference_contract.py
    test_binding_source_contract.py test_cli_packaging.py)
if(LSDTT_ENABLE_PYTHON_GEO_TESTS)
    set(_lsdtt_bf_files "")
    foreach(_f IN LISTS LSDTT_PYTHON_BINDING_FREE_TEST_FILES)
        list(APPEND _lsdtt_bf_files "${CMAKE_SOURCE_DIR}/python_bindings/tests/${_f}")
    endforeach()
    add_test(NAME unit_python_binding_free_tests
        COMMAND ${Python3_EXECUTABLE} -m pytest -q -p no:cacheprovider ${_lsdtt_bf_files}
        WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}")
    lsdtt_mark_test_runtime(python unit_python_binding_free_tests)
else()
    lsdtt_mark_test_unregistered(python_geo unit_python_binding_free_tests)
endif()
lsdtt_set_tests_properties(unit_python_binding_free_tests PROPERTIES
    LABELS "core;conventions;python"
    TIMEOUT 300)

# docs/developer_notes/test_inventory.md (tests by evidence category, domain,
# CLI and Python file) must match the registered tests and the pytest markers.
# The inventory describes one configuration (CANONICAL in the generator); with
# any other the test skips, since a different option set registers different
# tests. Regenerate with tools/generate_test_inventory.py --build-dir <build> --update.
if(LSDTT_ENABLE_PYTHON_TESTS)
    add_test(NAME unit_test_inventory_is_current
        COMMAND ${Python3_EXECUTABLE}
            "${CMAKE_SOURCE_DIR}/tools/generate_test_inventory.py"
            --check --build-dir "${CMAKE_BINARY_DIR}" --ctest "${CMAKE_CTEST_COMMAND}")
    lsdtt_mark_test_runtime(python unit_test_inventory_is_current)
else()
    lsdtt_mark_test_unregistered(python unit_test_inventory_is_current)
endif()
lsdtt_set_tests_properties(unit_test_inventory_is_current PROPERTIES
    LABELS "core;docs;python"
    SKIP_RETURN_CODE 77
    TIMEOUT 300)

# =============================================================================
# TEST REPORT AND LOGGED RUNS
#
#   cmake -P <build>/run_tests.cmake [ctest args]
#       ctest through `tee` into <build>/test-logs/ctest_<timestamp>.log, with
#       the arguments recorded in test-logs/last_run.txt; a killed run still
#       leaves the log of every test that finished.
#
#   cmake --build <build> --target test_report
#       one plain-text file for a user to send when tests fail: build, machine,
#       geospatial run time, test environment, run summary and the failing
#       tests' output, with the home directory and user name masked and no
#       environment dump (cmake/make_test_report.cmake).
#
# Both are CMake scripts, so they need no Python. The run-time GDAL/PROJ facts
# come from an existing CLI's --version; no program is written for the report.
# =============================================================================
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/cmake/run_tests.cmake.in"
               "${CMAKE_BINARY_DIR}/run_tests.cmake" @ONLY)

set(_lsdtt_report_options "")
foreach(_opt BUILD_TESTING LSDTT_ENABLE_PYTHON_TESTS LSDTT_ENABLE_PYTHON_GEO_TESTS
             LSDTT_ENABLE_NETWORK_TESTS LSDTT_ENABLE_BENCHMARKS
             LSDTT_ENABLE_PERF_TESTS LSDTT_GPU)
    if(NOT DEFINED ${_opt})
        set(_v OFF)
    else()
        set(_v "${${_opt}}")
    endif()
    string(APPEND _lsdtt_report_options "${_opt}=${_v} ")
endforeach()
set(_lsdtt_report_python "")
if(LSDTT_ENABLE_PYTHON_TESTS)
    set(_lsdtt_report_python "${Python3_EXECUTABLE}")
endif()
file(GENERATE OUTPUT "${CMAKE_BINARY_DIR}/lsdtt3_build_info.cmake" CONTENT
"# Generated at CMake generate time; read by tests/cmake/make_test_report.cmake.
set(LSDTT_REPORT_VERSION [==[${LSDTT3_FULL_VERSION}]==])
set(LSDTT_REPORT_BUILD_TYPE [==[${CMAKE_BUILD_TYPE}]==])
set(LSDTT_REPORT_CXX_ID [==[${CMAKE_CXX_COMPILER_ID}]==])
set(LSDTT_REPORT_CXX_VERSION [==[${CMAKE_CXX_COMPILER_VERSION}]==])
set(LSDTT_REPORT_CXX_STANDARD [==[${CMAKE_CXX_STANDARD}]==])
set(LSDTT_REPORT_GENERATOR [==[${CMAKE_GENERATOR}]==])
set(LSDTT_REPORT_CMAKE_VERSION [==[${CMAKE_VERSION}]==])
set(LSDTT_REPORT_OPTIONS [==[${_lsdtt_report_options}]==])
set(LSDTT_REPORT_HAS_OPENMP [==[${LSDTT_HAS_OPENMP}]==])
set(LSDTT_REPORT_HAS_FFTW [==[${LSDTT_HAS_FFTW}]==])
set(LSDTT_REPORT_HAS_EIGEN [==[${LSDTT_HAS_EIGEN}]==])
set(LSDTT_REPORT_HAS_NETCDF [==[${LSDTT_HAS_NETCDF}]==])
set(LSDTT_REPORT_GDAL_VERSION [==[${GDAL_VERSION}]==])
set(LSDTT_REPORT_GDAL_LIBRARY [==[${GDAL_LIBRARY}]==])
set(LSDTT_REPORT_CLI [==[$<TARGET_FILE:lsdtt-channel-extraction>]==])
set(LSDTT_REPORT_CTEST [==[${CMAKE_CTEST_COMMAND}]==])
set(LSDTT_REPORT_PYTHON [==[${_lsdtt_report_python}]==])
set(LSDTT_REPORT_GEO_TESTS [==[${LSDTT_ENABLE_PYTHON_GEO_TESTS}]==])
")

add_custom_target(test_report
    COMMAND ${CMAKE_COMMAND} -DBUILD_DIR=${CMAKE_BINARY_DIR}
            -P "${CMAKE_CURRENT_SOURCE_DIR}/cmake/make_test_report.cmake"
    USES_TERMINAL
    COMMENT "Writing the lsdtt3 test report")

# The report built from a canned build tree (cmake/report_selftest/), with an
# API key planted in the environment: every section, the failure and timeout
# entries, the harness logs, home-directory masking, and the key's absence.
add_test(NAME unit_test_report_selftest
    COMMAND ${CMAKE_COMMAND} -E env OT_API_KEY=planted-secret-8d3f
            ${CMAKE_COMMAND}
            -DSELFTEST_DIR=${CMAKE_CURRENT_SOURCE_DIR}/cmake/report_selftest
            -DSELFTEST_WORK=${CMAKE_CURRENT_BINARY_DIR}/test_report_selftest
            -P "${CMAKE_CURRENT_SOURCE_DIR}/cmake/make_test_report.cmake")
set_tests_properties(unit_test_report_selftest PROPERTIES
    LABELS "core;conventions"
    TIMEOUT 120)

# =============================================================================
# DEPENDENCY TIERS, ENVIRONMENT FLAGS AND THE LEM COMPONENT
#
# Every test gets exactly one RUNTIME TIER label, saying what it needs to start:
#
#   tier_native         a compiled lsdtt3 executable or `cmake -E`; no Python
#   tier_python_stdlib  the Python interpreter and its standard library only
#   tier_python_geo     Python plus third-party geospatial packages
#   tier_pixi           `pixi run` into the benchmark environment
#
# plus zero or more ENVIRONMENT FLAGS:
#
#   hermetic            tier_native, no skip path, no network, no host data --
#                       the reproducible floor: `ctest -L '^hermetic$'`
#   needs_network       downloads data (skips, exit 77, when offline)
#   needs_gdal_tools    shells out to gdalinfo / ogrinfo / gdal_translate
#                       (Ubuntu: gdal-bin, which libgdal-dev does not pull in)
#
# and, for the landscape evolution model, the COMPONENT label lem_component
# (every lsdtt-lem CLI test and the unit tests of the engines it drives).
#
# These are ORTHOGONAL to the domain labels above and to the evidence and
# level labels (tests/cmake/evidence_labels.cmake). Select with
# ANCHORED patterns: ctest -L/-LE take regexes, and `-L network` already
# matches synthetic_networks as well as network.
#
#   ctest -L '^hermetic$'                      native, offline, reproducible
#   ctest -LE '^lem_component$'                everything except the LEM
#   ctest -L '^lem_component$' -LE '^heavy$'   the LEM, cheap part
#
# The runtime tier comes from lsdtt_mark_test_runtime() at registration (top of
# this file), and is then PROPAGATED through fixtures: a test that requires a
# fixture produced by a Python test cannot run without Python, whatever its own
# command is (rigorous_drainage_analysis_threshold_equivalence is `cmake -E
# compare_files` on two Python-harness outputs). needs_network and
# needs_gdal_tools propagate the same way.
#
# The lists below come from an inventory of every registered command on
# 2026-09-18: an AST scan of each test script and the local modules it imports
# for third-party packages, and a scan for the GDAL programs it runs.
# =============================================================================

# Python tests whose scripts import third-party packages.
set(LSDTT_PYTHON_GEO_TESTS
    smoke_multipart_cli_sweep                       # geopandas numpy rasterio shapely
    rigorous_fluvial_surfaces_label_band_identity   # osgeo
    rigorous_fluvial_surfaces_empty_result          # osgeo
    rigorous_flat_surface_finder_admission_modes    # osgeo, in build_dem()
    rigorous_fluvial_surfaces_clubb2017_joins       # osgeo
    rigorous_fluvial_surfaces_cross_section_levels_workflow  # osgeo
    rigorous_fluvial_surfaces_whole_mask_segmentation        # osgeo
    unit_python_binding_free_tests                  # pytest
)

# Scripts other than run_cli_smoke.py that run GDAL's command-line programs.
# (run_cli_smoke.py tests are marked from their arguments at registration.)
set(LSDTT_GDAL_TOOL_SCRIPT_TESTS
    rigorous_cosmo_multi_accum_ratio                # gdalinfo
    regression_step_change_disk_roundtrip           # gdalinfo
    smoke_fluvial_surfaces_hand_input               # gdalinfo
    regression_landslide_scar_coupling              # gdal_translate
    regression_cosmo_oracle                         # gdal_translate
    regression_cosmo_oracle_heavy                   # gdal_translate
)

# Tests that download. The benchmarks and performance_surface_metrics fetch
# their DEM fixtures with lsdtt-fetch-raster when they are not cached.
set(LSDTT_NEEDS_NETWORK_PATTERN "^(network_|cosmo_realworld_|benchmark_)")
set(LSDTT_NEEDS_NETWORK_TESTS
    performance_surface_metrics
    smoke_flat_resolution_real_terrain_routing   # fetches its fixtures when not cached
)

# Native tests that are nevertheless not hermetic. Tests with a
# SKIP_RETURN_CODE are excluded automatically; list here only what that rule
# cannot see.
set(LSDTT_NONHERMETIC_TESTS
    unit_gdal_init_pinned       # stages the host's proj.db
)

# lem_component: every test named <kind>_lem_* (the lsdtt-lem CLI tests), plus
# the unit tests of the engines lsdtt-lem drives, whose names share no prefix.
# Deliberately NOT members: unit_hillslope_length_search and
# unit_hillslope_peclet (hillslope analysis), the hillslope_channel_coupling
# CLI tests, and unit_synthetic_cli_banners (a cross-CLI banner check that
# happens to run lsdtt-lem).
set(LSDTT_LEM_NAME_PATTERN "^(smoke|failure|regression|rigorous)_lem_")
set(LSDTT_LEM_ENGINE_TESTS
    # lsdtt_lem_cpu
    unit_lem_channel_coupling
    unit_lem_timestep
    unit_lem_forcing_history
    unit_lem_erodibility
    rigorous_lem_erodibility_analytic
    # landscape equilibration
    unit_landscape_equilibration
    unit_hillslope_length_statistic
    # hillslope solvers
    unit_hillslope_implicit
    unit_hillslope_newton
    unit_hillslope_steady_state
    unit_hillslope_transport_bc
    unit_hillslope_perron_qimp
    unit_hillslope_ren_modular_implicit
    unit_hillslope_critical_slope
    unit_hillslope_thread_count_invariance
    performance_hillslope_implicit
    # cellular hillslope
    unit_cellular_hillslope_references
    unit_cellular_hillslope_constitutive
    unit_cellular_hillslope_solver
    unit_cellular_hillslope_common_limit
    unit_hillslope_particle_io
    # stream power / steady-state channel pieces of the LEM
    unit_channel_threshold_units
    unit_chi_channel_head
)

get_property(_lsdtt_class_tests DIRECTORY PROPERTY TESTS)
get_property(_lsdtt_py_tests GLOBAL PROPERTY LSDTT_RUNTIME_python_TESTS)
get_property(_lsdtt_pixi_tests GLOBAL PROPERTY LSDTT_RUNTIME_pixi_TESTS)
get_property(_lsdtt_tool_tests GLOBAL PROPERTY LSDTT_NEEDS_gdal_tools_TESTS)
list(APPEND _lsdtt_tool_tests ${LSDTT_GDAL_TOOL_SCRIPT_TESTS})

get_property(_lsdtt_off_tests GLOBAL PROPERTY LSDTT_UNREGISTERED_TESTS)

# A name in a hand-written list that is not a registered test is a typo, and a
# typo here silently mislabels a test. Tests switched off by an option are
# exempt: the Python ones recorded by lsdtt_mark_test_unregistered(), and the
# benchmark / performance ones by name.
foreach(_listed IN LISTS LSDTT_PYTHON_GEO_TESTS LSDTT_GDAL_TOOL_SCRIPT_TESTS
                         LSDTT_NONHERMETIC_TESTS LSDTT_LEM_ENGINE_TESTS
                         _lsdtt_py_tests _lsdtt_pixi_tests)
    if(NOT "${_listed}" IN_LIST _lsdtt_class_tests
       AND NOT "${_listed}" IN_LIST _lsdtt_off_tests
       AND NOT _listed MATCHES "^(benchmark_|perf)")
        message(FATAL_ERROR
            "Test classification names '${_listed}', which is not a registered "
            "test. Fix the name in the DEPENDENCY TIERS lists.")
    endif()
endforeach()

# Pass 1: each test's own runtime tier and flags; which fixtures it sets up.
# Ranks let propagation take the most demanding tier: native < stdlib < geo < pixi.
set(_lsdtt_tier_names tier_native tier_python_stdlib tier_python_geo tier_pixi)
foreach(_test IN LISTS _lsdtt_class_tests)
    if("${_test}" IN_LIST _lsdtt_pixi_tests)
        set(_rank 3)
    elseif("${_test}" IN_LIST LSDTT_PYTHON_GEO_TESTS)
        set(_rank 2)
    elseif("${_test}" IN_LIST _lsdtt_py_tests)
        set(_rank 1)
    else()
        set(_rank 0)
    endif()
    set(_lsdtt_rank_${_test} ${_rank})

    set(_net FALSE)
    if(_test MATCHES "${LSDTT_NEEDS_NETWORK_PATTERN}" OR "${_test}" IN_LIST LSDTT_NEEDS_NETWORK_TESTS)
        set(_net TRUE)
    endif()
    set(_lsdtt_net_${_test} ${_net})

    set(_tools FALSE)
    if("${_test}" IN_LIST _lsdtt_tool_tests)
        set(_tools TRUE)
    endif()
    set(_lsdtt_tools_${_test} ${_tools})

    get_test_property(${_test} FIXTURES_SETUP _setups)
    if(_setups)
        foreach(_fx IN LISTS _setups)
            list(APPEND _lsdtt_fixture_setup_${_fx} ${_test})
        endforeach()
    endif()
endforeach()

# A test that requires a fixture no registered test sets up would run anyway
# (CTest does not refuse it) against outputs that were never produced. That is
# what switching the Python tests off would do to a native consumer of a
# Python fixture, so it is a configure error rather than a confusing failure.
foreach(_test IN LISTS _lsdtt_class_tests)
    get_test_property(${_test} FIXTURES_REQUIRED _required)
    if(NOT _required)
        continue()
    endif()
    foreach(_fx IN LISTS _required)
        if(NOT _lsdtt_fixture_setup_${_fx})
            message(FATAL_ERROR
                "Test '${_test}' requires fixture '${_fx}', which no registered "
                "test sets up in this configuration. If its setup is a Python "
                "test, register '${_test}' under LSDTT_ENABLE_PYTHON_TESTS too.")
        endif()
    endforeach()
endforeach()

# Pass 2: propagate through FIXTURES_REQUIRED until nothing changes (fixtures
# can chain, so one sweep is not enough).
set(_changed TRUE)
while(_changed)
    set(_changed FALSE)
    foreach(_test IN LISTS _lsdtt_class_tests)
        get_test_property(${_test} FIXTURES_REQUIRED _required)
        if(NOT _required)
            continue()
        endif()
        foreach(_fx IN LISTS _required)
            foreach(_setup IN LISTS _lsdtt_fixture_setup_${_fx})
                if(_lsdtt_rank_${_setup} GREATER _lsdtt_rank_${_test})
                    set(_lsdtt_rank_${_test} ${_lsdtt_rank_${_setup}})
                    set(_changed TRUE)
                endif()
                if(_lsdtt_net_${_setup} AND NOT _lsdtt_net_${_test})
                    set(_lsdtt_net_${_test} TRUE)
                    set(_changed TRUE)
                endif()
                if(_lsdtt_tools_${_setup} AND NOT _lsdtt_tools_${_test})
                    set(_lsdtt_tools_${_test} TRUE)
                    set(_changed TRUE)
                endif()
            endforeach()
        endforeach()
    endforeach()
endwhile()

# Pass 3: write the labels.
set(_lsdtt_tier_counts 0 0 0 0)
set(_lsdtt_n_hermetic 0)
set(_lsdtt_n_lem 0)
foreach(_test IN LISTS _lsdtt_class_tests)
    set(_rank ${_lsdtt_rank_${_test}})
    list(GET _lsdtt_tier_names ${_rank} _tier)
    set(_add ${_tier})

    if(_lsdtt_net_${_test})
        list(APPEND _add needs_network)
    endif()
    if(_lsdtt_tools_${_test})
        list(APPEND _add needs_gdal_tools)
    endif()

    get_test_property(${_test} SKIP_RETURN_CODE _skip)
    if(_rank EQUAL 0 AND NOT _skip AND NOT _lsdtt_net_${_test}
       AND NOT _test MATCHES "^(benchmark_|perf)"
       AND NOT "${_test}" IN_LIST LSDTT_NONHERMETIC_TESTS)
        list(APPEND _add hermetic)
        math(EXPR _lsdtt_n_hermetic "${_lsdtt_n_hermetic} + 1")
    endif()

    if(_test MATCHES "${LSDTT_LEM_NAME_PATTERN}" OR "${_test}" IN_LIST LSDTT_LEM_ENGINE_TESTS)
        list(APPEND _add lem_component)
        math(EXPR _lsdtt_n_lem "${_lsdtt_n_lem} + 1")
    endif()

    set_property(TEST ${_test} APPEND PROPERTY LABELS ${_add})

    list(GET _lsdtt_tier_counts ${_rank} _c)
    math(EXPR _c "${_c} + 1")
    list(REMOVE_AT _lsdtt_tier_counts ${_rank})
    list(INSERT _lsdtt_tier_counts ${_rank} ${_c})
endforeach()

list(GET _lsdtt_tier_counts 0 _n0)
list(GET _lsdtt_tier_counts 1 _n1)
list(GET _lsdtt_tier_counts 2 _n2)
list(GET _lsdtt_tier_counts 3 _n3)
message(STATUS "Test tiers    : native ${_n0} (hermetic ${_lsdtt_n_hermetic}), "
               "python_stdlib ${_n1}, python_geo ${_n2}, pixi ${_n3}; "
               "lem_component ${_lsdtt_n_lem}")

# Say what an OFF option left out, so a reconfigured build that has quietly
# shrunk (the options default OFF) is visible rather than just "all passed".
get_property(_lsdtt_off_py GLOBAL PROPERTY LSDTT_UNREGISTERED_python_TESTS)
get_property(_lsdtt_off_geo GLOBAL PROPERTY LSDTT_UNREGISTERED_python_geo_TESTS)
list(LENGTH _lsdtt_off_py _n_off_py)
list(LENGTH _lsdtt_off_geo _n_off_geo)
if(_n_off_py GREATER 0)
    message(WARNING
        "${_n_off_py} Python-driven tests are not registered "
        "(LSDTT_ENABLE_PYTHON_TESTS=OFF). They cover most of the CLIs. "
        "Enable with -DLSDTT_ENABLE_PYTHON_TESTS=ON (needs python3; some also "
        "need GDAL's command-line programs, Ubuntu package gdal-bin).")
endif()
if(_n_off_geo GREATER 0)
    message(STATUS "Test tiers    : ${_n_off_geo} geospatial Python tests not "
                   "registered (LSDTT_ENABLE_PYTHON_GEO_TESTS=OFF)")
endif()
get_property(_lsdtt_off_net GLOBAL PROPERTY LSDTT_UNREGISTERED_network_TESTS)
list(LENGTH _lsdtt_off_net _n_off_net)
if(NOT LSDTT_ENABLE_NETWORK_TESTS)
    message(STATUS "Test tiers    : tests that download data not registered "
                   "(LSDTT_ENABLE_NETWORK_TESTS=OFF), including the "
                   "${_n_off_net} cosmo_realworld_* tests")
endif()

# =============================================================================
# EXECUTION CONTRACT: no test reads the terminal, and every test has a timeout
#
# CTest hands the invoking terminal to each test as its stdin. So "stdin is a
# terminal" is TRUE inside an interactive `ctest` run, and any code path that
# prompts on a terminal blocks the whole suite with its prompt swallowed by
# the harness. Agent and CI runs have no terminal, so they pass and never see
# it. On 2026-09-18 failure_fetch_raster_ot_missing_key_strict sat at a hidden
# OpenTopography key prompt for 40 minutes that way; unit_secret's prompt
# check had the same exposure.
#
# LSDTT_NONINTERACTIVE=1 makes lsdtt3's credential lookup refuse to prompt
# (lsdtt_secret.hpp). It is APPENDED, so a test's own ENVIRONMENT survives.
# The Python harnesses also give every child /dev/null as stdin
# (scripts/lsdtt_test_env.py); this covers native tests and any script that
# does not.
#
# OMP_WAIT_POLICY=passive is set on every correctness test except benchmark_*
# and perf* (whose timings it would change) and LSDTT_SERIAL_TESTS (which run
# alone, so have nothing to contend with, and whose tight solver loops are
# slowed badly by it: unit_hillslope_implicit, ~205 s mean in
# CTestCostData.txt, was still running after ~30 min at ~120% CPU with it).
# It changes only how idle OpenMP threads wait, not thread counts or results. Without it, OpenMP tests run
# under `ctest -j` spin-wait against each other on oversubscribed cores:
# measured 2026-09-18 on an 8-core laptop, four concurrent copies of
# lsdtt-test-hillslope-steady-state took >120 s (killed) with the default
# policy and 7.6 s with passive, against 5.6 s for one copy alone; under
# `-j4`, unit_hillslope_steady_state and unit_landscape_equilibration (1 s and
# 16 s alone) exceeded 300 s.
#
# A test with no explicit TIMEOUT gets one here, so a hang fails in minutes
# rather than at CTest's 1500 s default. Explicit timeouts are left alone.
# Basis (build/Testing/Temporary/CTestCostData.txt, 2026-09-18): every
# failure_* test finishes in about 1 s, and the slowest non-heavy,
# non-performance test averages under 20 s.
#
#   failure_*                 60 s   argument and validation paths
#   heavy                   3600 s   as the serial block above
#   benchmark_* / perf*     unset    their options own their budgets
#   everything else          300 s
#
# This block is LAST so that it also reaches tests registered after the
# classification pass (unit_build_is_current).
# =============================================================================
# =============================================================================
# Re-run from the ingested parameters
#
# Every run writes <write_fname>_ingested_params.param: every value it used,
# read or default. Running the tool on that file (write_path moved to a fresh
# directory) must write the same files with the same bytes, or the record does
# not reproduce the run (it did not: float defaults were written with six
# decimals, so a tolerance of 1e-12 came back as 0). One smoke test per CLI
# does the re-run (run_cli_smoke.py, LSDTT_RERUN_FROM_INGESTED). Not covered,
# for want of a run_cli_smoke.py test to reuse: lsdtt-graphflood,
# lsdtt-lithology-erodibility, lsdtt-ridgeline-analysis,
# lsdtt-sediment-routing and lsdtt-cosmo-column-tool.
# =============================================================================
set(LSDTT_RERUN_FROM_INGESTED_TESTS
    smoke_advection_diffusion smoke_bankfull_width smoke_basin_extraction
    smoke_basin_stats smoke_basin_surface_metrics smoke_channel_extraction
    smoke_channel_network_analysis_area_scaling smoke_chi_analysis
    smoke_cosmo_landslides_crn_site_mode smoke_cosmo_prediction_multi_isotope
    smoke_cosmo_tool_lsdn smoke_cross_sections
    smoke_distributed_hydrology_flat_resolution_topmodel_float_step
    smoke_divide_migration_clusters smoke_drainage_analysis
    smoke_feature_finder smoke_fetch_erosion_rates_offline
    smoke_fetch_lithology_offline smoke_flat_surface_finder
    smoke_fluvial_surfaces smoke_hillslope_channel_coupling
    smoke_internal_drainage_both smoke_lem_cellular_particles
    smoke_linear_stats smoke_marsh_topography smoke_mountain_front
    smoke_outlet_matching_identity smoke_raster_calculator_difference
    smoke_raster_preprocessing smoke_slope_area
    smoke_spectral_analysis_valley_spacing smoke_stream_burning
    smoke_surface_metrics smoke_swath_profile_plane
    smoke_synthetic_networks_island smoke_synthetic_surfaces_phase5_defaults
    smoke_valley_centreline_flat_resolution_float_step
    smoke_water_centreline smoke_water_line_repair)
if(LSDTT_ENABLE_PYTHON_TESTS)
    foreach(_test IN LISTS LSDTT_RERUN_FROM_INGESTED_TESTS)
        if(NOT TEST ${_test})
            message(FATAL_ERROR "LSDTT_RERUN_FROM_INGESTED_TESTS names ${_test}, "
                "which is not a test. Pick another smoke test for that CLI.")
        endif()
        set_property(TEST ${_test} APPEND PROPERTY ENVIRONMENT
                     "LSDTT_RERUN_FROM_INGESTED=1")
    endforeach()
endif()

# Evidence (oracle) and level labels, and the guard that every test has them.
# After the re-run list, which makes its tests metamorphic, and after every
# add_test, so the guard sees the whole inventory.
include("${CMAKE_CURRENT_SOURCE_DIR}/cmake/evidence_labels.cmake")

get_property(_lsdtt_all_tests DIRECTORY PROPERTY TESTS)
foreach(_test IN LISTS _lsdtt_all_tests)
    set_property(TEST ${_test} APPEND PROPERTY ENVIRONMENT "LSDTT_NONINTERACTIVE=1")
    # Which GDAL the lsdtt3 binaries were built against, so the test launcher
    # can tell whether an active Pixi/conda environment's GDAL/PROJ data belong
    # to them (see the note beside find_package(Python3) at the top).
    set_property(TEST ${_test} APPEND PROPERTY ENVIRONMENT "LSDTT_TEST_GDAL_LIBRARY=${GDAL_LIBRARY}")

    if(NOT _test MATCHES "^(benchmark_|perf)" AND NOT "${_test}" IN_LIST LSDTT_SERIAL_TESTS)
        set_property(TEST ${_test} APPEND PROPERTY ENVIRONMENT "OMP_WAIT_POLICY=passive")
    endif()

    get_test_property(${_test} TIMEOUT _timeout)
    if(_timeout)
        continue()
    endif()
    if(_test MATCHES "^(benchmark_|perf)")
        continue()
    elseif(_test MATCHES "^failure_")
        set_tests_properties(${_test} PROPERTIES TIMEOUT 60)
    elseif("${_test}" IN_LIST LSDTT_HEAVY_TESTS)
        set_tests_properties(${_test} PROPERTIES TIMEOUT 3600)
    else()
        set_tests_properties(${_test} PROPERTIES TIMEOUT 300)
    endif()
endforeach()

# =============================================================================
# CONFIGURE SUMMARY
#
# What is registered, what is not and why, and the commands to run it -- with
# counts taken from the final labels, never typed in.
# =============================================================================
get_property(_lsdtt_sum_tests DIRECTORY PROPERTY TESTS)
set(_lsdtt_sum_native 0)
set(_lsdtt_sum_hermetic 0)
set(_lsdtt_sum_stdlib 0)
set(_lsdtt_sum_geo 0)
set(_lsdtt_sum_pixi 0)
set(_lsdtt_sum_lem 0)
set(_lsdtt_sum_lem_heavy 0)
set(_lsdtt_sum_net 0)
set(_lsdtt_sum_heavy 0)
foreach(_test IN LISTS _lsdtt_sum_tests)
    get_test_property(${_test} LABELS _l)
    foreach(_pair "tier_native;native" "hermetic;hermetic" "tier_python_stdlib;stdlib"
                  "tier_python_geo;geo" "tier_pixi;pixi" "lem_component;lem"
                  "needs_network;net" "heavy;heavy")
        list(GET _pair 0 _label)
        list(GET _pair 1 _key)
        if("${_label}" IN_LIST _l)
            math(EXPR _lsdtt_sum_${_key} "${_lsdtt_sum_${_key}} + 1")
        endif()
    endforeach()
    if("lem_component" IN_LIST _l AND "heavy" IN_LIST _l)
        math(EXPR _lsdtt_sum_lem_heavy "${_lsdtt_sum_lem_heavy} + 1")
    endif()
endforeach()
list(LENGTH _lsdtt_sum_tests _lsdtt_sum_total)

macro(_lsdtt_tier_line _name _count _option _n_off)
    if(${_option})
        message(STATUS "  ${_name}${_count} tests")
    else()
        message(STATUS "  ${_name}disabled (${_option}=OFF; ${_n_off} not registered)")
    endif()
endmacro()

set(_b "${CMAKE_BINARY_DIR}")
message(STATUS "")
message(STATUS "=== lsdtt3 test configuration (${_lsdtt_sum_total} tests) ===")
message(STATUS "  tier_native        : ${_lsdtt_sum_native} tests (hermetic ${_lsdtt_sum_hermetic})")
_lsdtt_tier_line("tier_python_stdlib : " "${_lsdtt_sum_stdlib}" LSDTT_ENABLE_PYTHON_TESTS "${_n_off_py}")
_lsdtt_tier_line("tier_python_geo    : " "${_lsdtt_sum_geo}" LSDTT_ENABLE_PYTHON_GEO_TESTS "${_n_off_geo}")
if(LSDTT_ENABLE_BENCHMARKS)
    message(STATUS "  tier_pixi          : ${_lsdtt_sum_pixi} tests (benchmarks)")
else()
    message(STATUS "  tier_pixi          : disabled (LSDTT_ENABLE_BENCHMARKS=OFF)")
endif()
if(LSDTT_ENABLE_NETWORK_TESTS)
    message(STATUS "  needs_network      : ${_lsdtt_sum_net} tests (downloads; skip when offline)")
else()
    message(STATUS "  needs_network      : disabled (LSDTT_ENABLE_NETWORK_TESTS=OFF; includes ${_n_off_net} cosmo_realworld_*)")
endif()
message(STATUS "  lem_component      : ${_lsdtt_sum_lem} tests (heavy ${_lsdtt_sum_lem_heavy})")
message(STATUS "  heavy              : ${_lsdtt_sum_heavy} tests")
message(STATUS "  evidence           : ${_lsdtt_ev_summary}more than one ${_lsdtt_ev_n_multi}")
message(STATUS "  level              : unit ${_lsdtt_ev_count_unit}, system ${_lsdtt_ev_count_system}")
if(LSDTT_ENABLE_PERF_TESTS)
    message(STATUS "  performance        : ON")
else()
    message(STATUS "  performance        : OFF (LSDTT_ENABLE_PERF_TESTS)")
endif()
if(LSDTT_ENABLE_PYTHON_TESTS)
    message(STATUS "  Python             : ${Python3_EXECUTABLE} (${Python3_VERSION})")
else()
    message(STATUS "  Python             : not used")
endif()
message(STATUS "Quick    : ctest --test-dir ${_b} -L '^hermetic$' -LE '^heavy$' --output-on-failure")
message(STATUS "Logged   : cmake -P ${_b}/run_tests.cmake -LE '^heavy$' -j4")
message(STATUS "No LEM   : ctest --test-dir ${_b} -LE '^lem_component$' --output-on-failure")
message(STATUS "LEM only : ctest --test-dir ${_b} -L '^lem_component$' --output-on-failure")
message(STATUS "Report   : cmake --build ${_b} --target test_report   (after a failing run)")
message(STATUS "================================================")
