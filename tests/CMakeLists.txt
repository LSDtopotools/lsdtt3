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
    foreach(_mod osgeo.gdal numpy rasterio geopandas shapely)
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
    message(STATUS "Test Python   : geospatial imports OK (osgeo.gdal, numpy, rasterio, geopandas, shapely)")
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
function(lsdtt_mark_harness_test TEST_NAME)
    lsdtt_mark_test_runtime(python ${TEST_NAME})
    foreach(_arg IN LISTS ARGN)
        if("${_arg}" IN_LIST LSDTT_HARNESS_GDAL_TOOL_ARGS)
            lsdtt_mark_test_needs(gdal_tools ${TEST_NAME})
            break()
        endif()
    endforeach()
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

    set(test_labels "cli")
    if(TEST_NAME MATCHES "^rigorous_")
        list(APPEND test_labels "rigorous")
    elseif(TEST_NAME MATCHES "^failure_")
        list(APPEND test_labels "failure")
    else()
        list(APPEND test_labels "smoke")
    endif()

    lsdtt_set_tests_properties(${TEST_NAME} PROPERTIES LABELS "${test_labels}")
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
    LABELS "smoke;cli;synthetic_networks"
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
    add_test(
        NAME smoke_synthetic_networks_polygon_area
        COMMAND
            ${Python3_EXECUTABLE} "${LSDTT_TEST_SCRIPT}"
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
    lsdtt_mark_test_runtime(python smoke_synthetic_networks_polygon_area)
    lsdtt_mark_test_needs(gdal_tools smoke_synthetic_networks_polygon_area)
else()
    lsdtt_mark_test_unregistered(python smoke_synthetic_networks_polygon_area)
endif()
lsdtt_set_tests_properties(smoke_synthetic_networks_polygon_area PROPERTIES
    LABELS "smoke;cli;synthetic_networks"
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

add_lsdtt_cli_smoke_test(
    smoke_nodata_edge_unfilled
    lsdtt-raster-preprocessing
    raster_preprocessing_nodata_edge_smoke.param.in
    tiny_edge_nodata_filled.tif
    --expect-valid-percent tiny_edge_nodata_filled.tif:96
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
    LABELS "cli;rigorous"
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
    LABELS "smoke;cli;distributed_hydrology"
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
    LABELS "smoke;cli;distributed_hydrology"
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
    LABELS "smoke;cli;distributed_hydrology"
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
    LABELS "smoke;cli;distributed_hydrology"
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
    LABELS "smoke;cli;distributed_hydrology"
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
    LABELS "smoke;cli;distributed_hydrology"
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
    LABELS "smoke;cli;distributed_hydrology"
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
    LABELS "smoke;cli;distributed_hydrology"
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
        LABELS "failure;cli;distributed_hydrology"
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
    add_test(NAME smoke_basin_surface_metrics
        COMMAND ${Python3_EXECUTABLE} "${LSDTT_TEST_SCRIPT}"
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
    lsdtt_mark_test_runtime(python smoke_basin_surface_metrics)
    lsdtt_mark_test_needs(gdal_tools smoke_basin_surface_metrics)
else()
    lsdtt_mark_test_unregistered(python smoke_basin_surface_metrics)
endif()
lsdtt_set_tests_properties(smoke_basin_surface_metrics PROPERTIES
    LABELS "smoke;cli;basin_surface_metrics"
    FIXTURES_REQUIRED basin_stats_fixture)

if(LSDTT_ENABLE_PYTHON_TESTS)
    add_test(NAME smoke_basin_surface_metrics_polyfit
        COMMAND ${Python3_EXECUTABLE} "${LSDTT_TEST_SCRIPT}"
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
    lsdtt_mark_test_runtime(python smoke_basin_surface_metrics_polyfit)
    lsdtt_mark_test_needs(gdal_tools smoke_basin_surface_metrics_polyfit)
else()
    lsdtt_mark_test_unregistered(python smoke_basin_surface_metrics_polyfit)
endif()
lsdtt_set_tests_properties(smoke_basin_surface_metrics_polyfit PROPERTIES
    LABELS "smoke;cli;basin_surface_metrics"
    FIXTURES_REQUIRED basin_stats_fixture)

if(LSDTT_ENABLE_PYTHON_TESTS)
    add_test(NAME smoke_basin_surface_metrics_expanding_square
        COMMAND ${Python3_EXECUTABLE} "${LSDTT_TEST_SCRIPT}"
            --exe "$<TARGET_FILE:lsdtt-basin-surface-metrics>"
            --cli-src-dir "${CMAKE_SOURCE_DIR}/cli"
            --param-template "${LSDTT_TEST_PARAM_DIR}/basin_surface_metrics_expanding_square.param.in"
            --read-path "${LSDTT_BASIN_STATS_FIXTURE_DIR}/out"
            --work-dir "${CMAKE_CURRENT_BINARY_DIR}/smoke_basin_surface_metrics_expanding_square"
            --expect bsm_square_outline_curvature_flips.csv)
    lsdtt_mark_test_runtime(python smoke_basin_surface_metrics_expanding_square)
else()
    lsdtt_mark_test_unregistered(python smoke_basin_surface_metrics_expanding_square)
endif()
lsdtt_set_tests_properties(smoke_basin_surface_metrics_expanding_square PROPERTIES
    LABELS "smoke;cli;basin_surface_metrics;snapping"
    FIXTURES_REQUIRED basin_stats_fixture)

if(LSDTT_ENABLE_PYTHON_TESTS)
    add_test(NAME rigorous_basin_surface_metrics_legacy_threshold
        COMMAND ${Python3_EXECUTABLE} "${LSDTT_TEST_SCRIPT}"
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
    lsdtt_mark_test_runtime(python rigorous_basin_surface_metrics_legacy_threshold)
else()
    lsdtt_mark_test_unregistered(python rigorous_basin_surface_metrics_legacy_threshold)
endif()
lsdtt_set_tests_properties(rigorous_basin_surface_metrics_legacy_threshold PROPERTIES
    LABELS "rigorous;cli;basin_surface_metrics;compatibility"
    FIXTURES_REQUIRED basin_stats_fixture)

add_lsdtt_cli_smoke_test(
    failure_basin_surface_metrics_unknown_snap
    lsdtt-basin-surface-metrics
    basin_surface_metrics_unknown_snap.param.in
    unused
    --expect-failure
    --stderr-contains "FATAL: basin_outlet_snapping_method must be one of: none, downstream, expanding_square."
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
    add_test(NAME smoke_ridgeline_range_open
        COMMAND ${Python3_EXECUTABLE} "${LSDTT_TEST_SCRIPT}"
            --exe "$<TARGET_FILE:lsdtt-ridgeline-analysis>"
            --cli-src-dir "${CMAKE_SOURCE_DIR}/cli"
            --param-template "${LSDTT_TEST_PARAM_DIR}/ridgeline_range_oooo.param.in"
            --read-path "${LSDTT_BASIN_STATS_FIXTURE_DIR}/out"
            --work-dir "${CMAKE_CURRENT_BINARY_DIR}/smoke_ridgeline_range_open"
            --expect rg_oooo_raw_dem_3x3_polyfit_slope.tif
            --expect rg_oooo_raw_dem_3x3_polyfit_maximum_curvature.tif
            --expect-raster-size rg_oooo_raw_dem_3x3_polyfit_slope.tif:200x200
            # C2 must actually resolve corners and segment the boundary.
            --expect-file-contains "../stdout.log:Resolving DEM corners"
            --expect-file-contains "../stdout.log:Boundary segment corner 3 assigned"
            # ORACLE for C2+C3 under open boundaries.
            --expect-sha256 rg_oooo_raw_dem_3x3_polyfit_slope.tif:c6a0d7322aa1d1f510db0fac8a2544c052238e1782d74bd68805e04a37b8b9c4)
    lsdtt_mark_test_runtime(python smoke_ridgeline_range_open)
    lsdtt_mark_test_needs(gdal_tools smoke_ridgeline_range_open)
else()
    lsdtt_mark_test_unregistered(python smoke_ridgeline_range_open)
endif()
lsdtt_set_tests_properties(smoke_ridgeline_range_open PROPERTIES
    LABELS "smoke;cli;ridgeline"
    FIXTURES_REQUIRED basin_stats_fixture)

if(LSDTT_ENABLE_PYTHON_TESTS)
    add_test(NAME smoke_ridgeline_range_periodic
        COMMAND ${Python3_EXECUTABLE} "${LSDTT_TEST_SCRIPT}"
            --exe "$<TARGET_FILE:lsdtt-ridgeline-analysis>"
            --cli-src-dir "${CMAKE_SOURCE_DIR}/cli"
            --param-template "${LSDTT_TEST_PARAM_DIR}/ridgeline_range_opop.param.in"
            --read-path "${LSDTT_BASIN_STATS_FIXTURE_DIR}/out"
            --work-dir "${CMAKE_CURRENT_BINARY_DIR}/smoke_ridgeline_range_periodic"
            --expect rg_opop_raw_dem_3x3_polyfit_slope.tif
            --expect-raster-size rg_opop_raw_dem_3x3_polyfit_slope.tif:200x200
            --expect-file-contains "../stdout.log:Resolving DEM corners"
            # ORACLE for the PERIODIC padding branch. This hash differs from the
            # open-boundary one above, which is the evidence that the wrap-around
            # padding in run_raw_dem_3x3_polyfit_with_bcs actually does something --
            # if a refactor made the two equal, both oracles could not hold.
            --expect-sha256 rg_opop_raw_dem_3x3_polyfit_slope.tif:f7beffb5d95fb65e7ca6108624baac5f8c7ae60399622c2f235acb94a2bd9eae)
    lsdtt_mark_test_runtime(python smoke_ridgeline_range_periodic)
    lsdtt_mark_test_needs(gdal_tools smoke_ridgeline_range_periodic)
else()
    lsdtt_mark_test_unregistered(python smoke_ridgeline_range_periodic)
endif()
lsdtt_set_tests_properties(smoke_ridgeline_range_periodic PROPERTIES
    LABELS "smoke;cli;ridgeline"
    FIXTURES_REQUIRED basin_stats_fixture)

# lsdtt-chi-analysis: the outlet-FILE path (C5). The other chi tests only reach
# resolve_outlet_flats through its auto-detect branch.
if(LSDTT_ENABLE_PYTHON_TESTS)
    add_test(NAME smoke_chi_analysis_outlet_file
        COMMAND ${Python3_EXECUTABLE} "${LSDTT_TEST_SCRIPT}"
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
    lsdtt_mark_test_runtime(python smoke_chi_analysis_outlet_file)
else()
    lsdtt_mark_test_unregistered(python smoke_chi_analysis_outlet_file)
endif()
lsdtt_set_tests_properties(smoke_chi_analysis_outlet_file PROPERTIES
    LABELS "smoke;cli;chi"
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
    add_test(NAME smoke_linear_stats
        COMMAND ${Python3_EXECUTABLE} "${LSDTT_TEST_SCRIPT}"
            --exe "$<TARGET_FILE:lsdtt-linear-stats>"
            --cli-src-dir "${CMAKE_SOURCE_DIR}/cli"
            --param-template "${LSDTT_TEST_PARAM_DIR}/linear_stats_smoke.param.in"
            --read-path "${LSDTT_BASIN_STATS_FIXTURE_DIR}/out"
            --work-dir "${CMAKE_CURRENT_BINARY_DIR}/smoke_linear_stats"
            --expect ls_linear_stats.fgb
            --expect-vector-feature-min ls_linear_stats.fgb:2
            --expect-sha256 ls_linear_stats.fgb:5089ab341d0a212acc0fee3cba76c3e414aaa89bcbf7e017aabd93fa3988224a)
    lsdtt_mark_test_runtime(python smoke_linear_stats)
    lsdtt_mark_test_needs(gdal_tools smoke_linear_stats)
else()
    lsdtt_mark_test_unregistered(python smoke_linear_stats)
endif()
lsdtt_set_tests_properties(smoke_linear_stats PROPERTIES
    LABELS "smoke;cli;linear_stats"
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
    add_test(NAME fixture_stream_burning_dem
        COMMAND ${Python3_EXECUTABLE} "${LSDTT_TEST_SCRIPT}"
            --exe "$<TARGET_FILE:lsdtt-synthetic-surfaces>"
            --cli-src-dir "${CMAKE_SOURCE_DIR}/cli"
            --param-template "${LSDTT_TEST_PARAM_DIR}/stream_burning_fixture_dem.param.in"
            --read-path "${LSDTT_TEST_DATA_DIR}"
            --work-dir "${LSDTT_STREAM_BURNING_FIXTURE_DIR}"
            --expect sb_diamond_square.tif)
    lsdtt_mark_test_runtime(python fixture_stream_burning_dem)
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
    add_test(NAME smoke_stream_burning
        COMMAND ${Python3_EXECUTABLE} "${LSDTT_TEST_SCRIPT}"
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
    lsdtt_mark_test_runtime(python smoke_stream_burning)
else()
    lsdtt_mark_test_unregistered(python smoke_stream_burning)
endif()
lsdtt_set_tests_properties(smoke_stream_burning PROPERTIES
    LABELS "smoke;cli;stream_burning"
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
    add_test(NAME fixture_cross_sections_dem
        COMMAND ${Python3_EXECUTABLE} "${LSDTT_TEST_SCRIPT}"
            --exe "$<TARGET_FILE:lsdtt-synthetic-surfaces>"
            --cli-src-dir "${CMAKE_SOURCE_DIR}/cli"
            --param-template "${LSDTT_TEST_PARAM_DIR}/cross_sections_fixture_dem.param.in"
            --read-path "${LSDTT_TEST_DATA_DIR}"
            --work-dir "${LSDTT_CROSS_SECTIONS_FIXTURE_DIR}"
            --expect xs_terrace_staircase.tif)
    lsdtt_mark_test_runtime(python fixture_cross_sections_dem)
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
    add_test(NAME fixture_basin_stats_dem
        COMMAND ${Python3_EXECUTABLE} "${LSDTT_TEST_SCRIPT}"
            --exe "$<TARGET_FILE:lsdtt-synthetic-surfaces>"
            --cli-src-dir "${CMAKE_SOURCE_DIR}/cli"
            --param-template "${LSDTT_TEST_PARAM_DIR}/basin_stats_fixture_dem.param.in"
            --read-path "${LSDTT_TEST_DATA_DIR}"
            --work-dir "${LSDTT_BASIN_STATS_FIXTURE_DIR}"
            --expect basin_stats_diamond_square.tif)
    lsdtt_mark_test_runtime(python fixture_basin_stats_dem)
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
    add_test(NAME smoke_basin_stats
        COMMAND ${Python3_EXECUTABLE} "${LSDTT_TEST_SCRIPT}"
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
            --expect-sha256 bs_basin_stats.csv:6b8129739e73739716eb565ba822a85c4952efd55eb8ff0e3f41a1e10efa75a5
            # Canonical-default regression: the parameter is deliberately absent
            # from the fixture, so this pins the shared default and its behaviour.
            --stdout-contains "threshold = 1000 px"
            --stderr-contains "'basin_outlet_csv' is deprecated and 'basin_outlet_fname' was also supplied")
    lsdtt_mark_test_runtime(python smoke_basin_stats)
    lsdtt_mark_test_needs(gdal_tools smoke_basin_stats)
else()
    lsdtt_mark_test_unregistered(python smoke_basin_stats)
endif()
lsdtt_set_tests_properties(smoke_basin_stats PROPERTIES
    LABELS "smoke;cli;basin_stats"
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

add_lsdtt_cli_smoke_test(
    smoke_feature_finder
    lsdtt-feature-finder
    feature_finder_smoke.param.in
    ff_channel_candidate_score.tif
)

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
    LABELS "smoke;cli;divide_migration"
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
    LABELS "smoke;cli;divide_migration"
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
    LABELS "smoke;cli;divide_migration"
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
    LABELS "failure;cli;divide_migration"
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
    LABELS "failure;cli;divide_migration"
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
    LABELS "smoke;cli;divide_migration"
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
    LABELS "smoke;cli;divide_migration"
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
    LABELS "smoke;cli;divide_migration"
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
    LABELS "smoke;cli;hillslope"
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
lsdtt_set_tests_properties(smoke_lem_snap_default PROPERTIES LABELS "smoke;cli;lem")

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
    LABELS "smoke;cli;hillslope"
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
    LABELS "smoke;cli;hillslope")

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
    LABELS "smoke;cli;hillslope")

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
    LABELS "smoke;cli;hillslope")

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
    LABELS "smoke;cli;hillslope;cellular"
    FIXTURES_SETUP cellular_particle_fixture)

if(LSDTT_ENABLE_PYTHON_TESTS)
    add_test(
        NAME smoke_lem_cellular_particles_restart
        COMMAND ${Python3_EXECUTABLE} "${LSDTT_TEST_SCRIPT}"
                --exe "$<TARGET_FILE:lsdtt-lem>"
                --param-template "${LSDTT_TEST_PARAM_DIR}/lem_cellular_particles_restart.param.in"
                --read-path "${CMAKE_CURRENT_BINARY_DIR}/smoke_lem_cellular_particles/out"
                --work-dir "${CMAKE_CURRENT_BINARY_DIR}/smoke_lem_cellular_particles_restart"
                --expect cellularrestart_modelled_topography.tif
                --expect cellularrestart_cellular_particles.csv
                --stdout-contains "Loading cellular particles:"
    )
    lsdtt_mark_test_runtime(python smoke_lem_cellular_particles_restart)
else()
    lsdtt_mark_test_unregistered(python smoke_lem_cellular_particles_restart)
endif()
lsdtt_set_tests_properties(smoke_lem_cellular_particles_restart PROPERTIES
    LABELS "smoke;cli;hillslope;cellular"
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
    LABELS "failure;cli;hillslope")

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

# D8 on a diverging hillslope manufactures channel heads, so the drainage
# density would be an artefact of the router.
add_lsdtt_cli_smoke_test(
    failure_lem_coupling_d8_conflict
    lsdtt-lem
    lem_coupling_d8_conflict.param.in
    unused_expected_output
    --expect-failure
    --stderr-contains "enable_d8_fluvial_erosion"
    --stderr-contains "manufactures channel"
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
    foreach(bench_fixture alpine_small lakeland_small)
        foreach(bench_case fill fill_epsilon breach accum_d8)
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
lsdtt_set_tests_properties(unit_topmodel_column PROPERTIES LABELS "unit;network")

# HYPE port. The oracle under
# tests/expected/hype is generated by tests/scripts/generate_hype_oracle.py;
# these tests only read committed files.
add_executable(lsdtt-test-distributed-hydrology-forcing cpp/test_distributed_hydrology_forcing.cpp)
target_link_libraries(lsdtt-test-distributed-hydrology-forcing PRIVATE lsdtt3_cpu)
add_test(NAME unit_distributed_hydrology_forcing
         COMMAND lsdtt-test-distributed-hydrology-forcing
                 "${CMAKE_CURRENT_BINARY_DIR}/distributed_hydrology_forcing_output")
lsdtt_set_tests_properties(unit_distributed_hydrology_forcing PROPERTIES LABELS "unit;distributed_hydrology")

add_executable(lsdtt-test-hype-parameters cpp/test_hype_parameters.cpp)
target_link_libraries(lsdtt-test-hype-parameters PRIVATE lsdtt3_cpu)
add_test(NAME unit_hype_parameters
         COMMAND lsdtt-test-hype-parameters "${CMAKE_CURRENT_SOURCE_DIR}/expected/hype"
                 "${CMAKE_CURRENT_BINARY_DIR}/hype_parameters_output")
lsdtt_set_tests_properties(unit_hype_parameters PROPERTIES LABELS "unit;distributed_hydrology")

add_executable(lsdtt-test-hype-column cpp/test_hype_column.cpp)
target_link_libraries(lsdtt-test-hype-column PRIVATE lsdtt3_cpu)
add_test(NAME unit_hype_column
         COMMAND lsdtt-test-hype-column "${CMAKE_CURRENT_SOURCE_DIR}/expected/hype")
lsdtt_set_tests_properties(unit_hype_column PROPERTIES LABELS "unit;distributed_hydrology")

add_executable(lsdtt-test-hype-routing cpp/test_hype_routing.cpp)
target_link_libraries(lsdtt-test-hype-routing PRIVATE lsdtt3_cpu)
add_test(NAME unit_hype_routing
         COMMAND lsdtt-test-hype-routing "${CMAKE_CURRENT_SOURCE_DIR}/expected/hype")
lsdtt_set_tests_properties(unit_hype_routing PROPERTIES LABELS "unit;distributed_hydrology")

add_executable(lsdtt-test-hype-model cpp/test_hype_model.cpp)
target_link_libraries(lsdtt-test-hype-model PRIVATE lsdtt3_cpu)
add_test(NAME unit_hype_model
         COMMAND lsdtt-test-hype-model "${CMAKE_CURRENT_SOURCE_DIR}/expected/hype")
lsdtt_set_tests_properties(unit_hype_model PROPERTIES LABELS "unit;distributed_hydrology")

add_executable(lsdtt-test-hype-raster-setup cpp/test_hype_raster_setup.cpp)
target_link_libraries(lsdtt-test-hype-raster-setup PRIVATE lsdtt3_cpu)
add_test(NAME unit_hype_raster_setup
         COMMAND lsdtt-test-hype-raster-setup "${CMAKE_CURRENT_BINARY_DIR}/hype_raster_setup_output")
lsdtt_set_tests_properties(unit_hype_raster_setup PROPERTIES LABELS "unit;distributed_hydrology")

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
lsdtt_set_tests_properties(unit_hype_thread_invariance PROPERTIES LABELS "unit;distributed_hydrology")

# HYPE timing harness. Built always so it
# can be run by hand; registered only under LSDTT_ENABLE_PERF_TESTS below.
add_executable(lsdtt-test-hype-performance cpp/test_hype_performance.cpp)
target_link_libraries(lsdtt-test-hype-performance PRIVATE lsdtt3_cpu)

add_executable(lsdtt-test-subbasin-partition cpp/test_subbasin_partition.cpp)
target_link_libraries(lsdtt-test-subbasin-partition PRIVATE lsdtt3_cpu)
add_test(NAME unit_subbasin_partition COMMAND lsdtt-test-subbasin-partition)
lsdtt_set_tests_properties(unit_subbasin_partition PROPERTIES LABELS "unit;distributed_hydrology")

# The channel mask ladder. Takes an output
# directory for the small GeoJSON fixture it writes at run time, so nothing is
# committed to tests/data.
add_executable(lsdtt-test-channel-extent cpp/test_channel_extent.cpp)
target_link_libraries(lsdtt-test-channel-extent PRIVATE lsdtt3_cpu)
add_test(NAME unit_channel_extent
         COMMAND lsdtt-test-channel-extent "${CMAKE_CURRENT_BINARY_DIR}/channel_extent_output")
lsdtt_set_tests_properties(unit_channel_extent PROPERTIES LABELS "unit;network")
file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/channel_extent_output")

# The routing ladder.
add_executable(lsdtt-test-hydro-routing cpp/test_hydro_routing.cpp)
target_link_libraries(lsdtt-test-hydro-routing PRIVATE lsdtt3_cpu)
add_test(NAME unit_hydro_routing COMMAND lsdtt-test-hydro-routing)
lsdtt_set_tests_properties(unit_hydro_routing PROPERTIES LABELS "unit;network")

# Canopy interception.
add_executable(lsdtt-test-canopy-interception cpp/test_canopy_interception.cpp)
target_link_libraries(lsdtt-test-canopy-interception PRIVATE lsdtt3_cpu)
add_test(NAME unit_canopy_interception COMMAND lsdtt-test-canopy-interception)
lsdtt_set_tests_properties(unit_canopy_interception PROPERTIES LABELS "unit;network")

# The roughness partition and the channel-n
# ladder.
add_executable(lsdtt-test-hydro-roughness cpp/test_hydro_roughness.cpp)
target_link_libraries(lsdtt-test-hydro-roughness PRIVATE lsdtt3_cpu)
add_test(NAME unit_hydro_roughness COMMAND lsdtt-test-hydro-roughness)
lsdtt_set_tests_properties(unit_hydro_roughness PROPERTIES LABELS "unit;network")

# Land-cover scenarios.
add_executable(lsdtt-test-landcover-scenario cpp/test_landcover_scenario.cpp)
target_link_libraries(lsdtt-test-landcover-scenario PRIVATE lsdtt3_cpu)
add_test(NAME unit_landcover_scenario COMMAND lsdtt-test-landcover-scenario)
lsdtt_set_tests_properties(unit_landcover_scenario PROPERTIES LABELS "unit;network")

# Infiltration excess (exp-K Green-Ampt).
# Series, closed-form capacity, water balance and dimensional-form checks; no
# external reference input.
add_executable(lsdtt-test-infiltration-excess cpp/test_infiltration_excess.cpp)
target_link_libraries(lsdtt-test-infiltration-excess PRIVATE lsdtt3_cpu)
add_test(NAME unit_infiltration_excess COMMAND lsdtt-test-infiltration-excess)
lsdtt_set_tests_properties(unit_infiltration_excess PROPERTIES LABELS "unit;network")

# The dynamic (prognostic) saturated zone.
add_executable(lsdtt-test-topmodel-dynamic cpp/test_topmodel_dynamic.cpp)
target_link_libraries(lsdtt-test-topmodel-dynamic PRIVATE lsdtt3_cpu)
add_test(NAME unit_topmodel_dynamic COMMAND lsdtt-test-topmodel-dynamic)
lsdtt_set_tests_properties(unit_topmodel_dynamic PROPERTIES LABELS "unit;network")

# Velocity field, travel time, width
# function. "hydro" is claimed by the wetness_index -> network domain rule
# below, so this stays in the network compartment with its Phase 0 sibling.
add_executable(lsdtt-test-hydro-velocity-field cpp/test_hydro_velocity_field.cpp)
target_link_libraries(lsdtt-test-hydro-velocity-field PRIVATE lsdtt3_cpu)
add_test(NAME unit_hydro_velocity_field COMMAND lsdtt-test-hydro-velocity-field)
lsdtt_set_tests_properties(unit_hydro_velocity_field PROPERTIES LABELS "unit;network")

# The shared TWI. Guards the core function
# directly, including the float precision contract that the CLI checksum test
# was MEASURED not to discriminate -- see the header of test_wetness_index.cpp.
add_executable(lsdtt-test-wetness-index cpp/test_wetness_index.cpp)
target_link_libraries(lsdtt-test-wetness-index PRIVATE lsdtt3_cpu)
add_test(NAME unit_wetness_index COMMAND lsdtt-test-wetness-index)
lsdtt_set_tests_properties(unit_wetness_index PROPERTIES LABELS "unit;network")

add_executable(lsdtt-test-image-filters cpp/test_image_filters.cpp)
target_link_libraries(lsdtt-test-image-filters PRIVATE lsdtt3_utils)
add_test(NAME unit_image_filters COMMAND lsdtt-test-image-filters)
lsdtt_set_tests_properties(unit_image_filters PROPERTIES LABELS "unit;image_ops")

# Native raster image operations vs outputs frozen from OpenCV 4.13 before the
# OpenCV dependency was removed (tests/expected/image_ops/), plus invariants.
add_executable(lsdtt-test-image-ops cpp/test_image_ops.cpp)
target_link_libraries(lsdtt-test-image-ops PRIVATE lsdtt3_cpu)
add_test(NAME unit_image_ops
         COMMAND lsdtt-test-image-ops ${CMAKE_CURRENT_SOURCE_DIR}/expected/image_ops)
lsdtt_set_tests_properties(unit_image_ops PROPERTIES LABELS "unit;image_ops")

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
lsdtt_set_tests_properties(unit_surface_metrics PROPERTIES LABELS "unit;surface_metrics")

# Vertical units: the feet/metres trap and the plausibility warning that makes
# it discoverable (docs/developer_notes/testing_and_benchmarking.md).
add_executable(lsdtt-test-vertical-units cpp/test_vertical_units.cpp)
target_link_libraries(lsdtt-test-vertical-units PRIVATE lsdtt3_cpu)

add_test(NAME unit_vertical_units COMMAND lsdtt-test-vertical-units)
lsdtt_set_tests_properties(unit_vertical_units PROPERTIES LABELS "unit;raster")

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
lsdtt_set_tests_properties(unit_gdal_init_idempotent PROPERTIES LABELS "unit;gdal")

# The case the module exists for: a poisoned PROJ_DATA must be ignored when a
# bundled proj.db is present. Needs a real proj.db to stage; skips without one.
add_test(NAME unit_gdal_init_pinned
         COMMAND lsdtt-test-gdal-init pinned
                 "${CMAKE_CURRENT_BINARY_DIR}/gdal_init_staging")
lsdtt_set_tests_properties(unit_gdal_init_pinned PROPERTIES
    LABELS "unit;gdal"
    SKIP_RETURN_CODE 77
    ENVIRONMENT "LSDTT_TEST_PROJ_DB=${LSDTT_TEST_PROJ_DB}")

# The conda-forge case: nothing bundled, so the user's environment must still
# be honoured. This is the half that a naive "always override" would break.
add_test(NAME unit_gdal_init_fallthrough
         COMMAND lsdtt-test-gdal-init fallthrough)
lsdtt_set_tests_properties(unit_gdal_init_fallthrough PROPERTIES LABELS "unit;gdal")

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
lsdtt_set_tests_properties(unit_gdal_format_support PROPERTIES LABELS "unit;gdal")

# giuh_discharge: volume conservation and timestep invariance. Added 2026-09-16
# with the fix for the unused dt_rainfall_s defect; the GIUH code had no tests
# at all before this, which is how the defect survived.
add_executable(lsdtt-test-giuh-discharge cpp/test_giuh_discharge.cpp)
target_link_libraries(lsdtt-test-giuh-discharge PRIVATE lsdtt3_cpu)
add_test(NAME unit_giuh_discharge COMMAND lsdtt-test-giuh-discharge)
lsdtt_set_tests_properties(unit_giuh_discharge PROPERTIES LABELS "unit;giuh")

if(LSDTT_ENABLE_PYTHON_TESTS)
    add_test(NAME unit_gdal_init_boundary
        COMMAND ${Python3_EXECUTABLE}
            "${CMAKE_CURRENT_SOURCE_DIR}/scripts/check_gdal_boundary.py")
    lsdtt_mark_test_runtime(python unit_gdal_init_boundary)
else()
    lsdtt_mark_test_unregistered(python unit_gdal_init_boundary)
endif()
lsdtt_set_tests_properties(unit_gdal_init_boundary PROPERTIES
    LABELS "unit;gdal;conventions")

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
lsdtt_set_tests_properties(unit_cli_helper_string_utils PROPERTIES LABELS "unit")

# B2 output_extension(raster_output_format, read_ext) and B3 csv_number.
add_executable(lsdtt-test-cli-helper-formats cpp/test_cli_helper_formats.cpp)
target_link_libraries(lsdtt-test-cli-helper-formats PRIVATE lsdtt3_core)
add_test(NAME unit_cli_helper_formats COMMAND lsdtt-test-cli-helper-formats)
lsdtt_set_tests_properties(unit_cli_helper_formats PROPERTIES LABELS "unit")

# D5 quote-aware lsdtt::csv_split_row and load_csv_table.
add_executable(lsdtt-test-csv-table-quoting cpp/test_csv_table_quoting.cpp)
target_link_libraries(lsdtt-test-csv-table-quoting PRIVATE lsdtt3_core)
add_test(NAME unit_csv_table_quoting COMMAND lsdtt-test-csv-table-quoting)
lsdtt_set_tests_properties(unit_csv_table_quoting PROPERTIES LABELS "unit")

# B4 the polyfit metric name/suffix table beside the polyfit_metric enum.
add_executable(lsdtt-test-surface-metrics-polyfit-metric-table
               cpp/test_surface_metrics_polyfit_metric_table.cpp)
target_link_libraries(lsdtt-test-surface-metrics-polyfit-metric-table PRIVATE lsdtt3_cpu)
add_test(NAME unit_surface_metrics_polyfit_metric_table
         COMMAND lsdtt-test-surface-metrics-polyfit-metric-table)
lsdtt_set_tests_properties(unit_surface_metrics_polyfit_metric_table
                     PROPERTIES LABELS "unit;surface_metrics")

add_executable(lsdtt-test-basin-outline-metrics
               cpp/test_basin_outline_metrics.cpp)
target_link_libraries(lsdtt-test-basin-outline-metrics PRIVATE lsdtt3_cpu)
add_test(NAME unit_basin_outline_metrics
         COMMAND lsdtt-test-basin-outline-metrics)
lsdtt_set_tests_properties(unit_basin_outline_metrics
                     PROPERTIES LABELS "unit;basin_surface_metrics")

add_executable(lsdtt-test-param-parser-help cpp/test_param_parser_help.cpp)
target_link_libraries(lsdtt-test-param-parser-help PRIVATE lsdtt3_core)

add_test(NAME unit_param_parser_help COMMAND lsdtt-test-param-parser-help)
lsdtt_set_tests_properties(unit_param_parser_help PROPERTIES LABELS "unit")

add_executable(lsdtt-test-common-parameters cpp/test_common_parameters.cpp)
target_link_libraries(lsdtt-test-common-parameters PRIVATE lsdtt3_core)

add_test(NAME unit_common_parameters COMMAND lsdtt-test-common-parameters)
lsdtt_set_tests_properties(unit_common_parameters PROPERTIES LABELS "unit;core;parameters")

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
lsdtt_set_tests_properties(unit_synthetic_networks PROPERTIES
    LABELS "unit;synthetic_networks"
    FIXTURES_SETUP synthetic_networks_fixture)
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
    LABELS "unit;synthetic_networks" TIMEOUT 1800)

add_test(NAME unit_flow_properties
         COMMAND lsdtt-test-flow-properties
                 ${CMAKE_CURRENT_BINARY_DIR}/flow_properties_output)
lsdtt_set_tests_properties(unit_flow_properties PROPERTIES LABELS "unit;flow")
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/flow_properties_output)

# BYTE-IDENTITY PIN for the D8 flow object. Pins were generated from the code
# BEFORE the section-9.4 speed refactor of lsdtt_flowinfo (commit 7e81a39); a
# speed change must never need them regenerated. Integer-terrain fixtures only,
# so the hashes are reproducible on any machine. See the file header.
add_executable(lsdtt-test-flow-d8-identity cpp/test_flow_d8_identity.cpp)
target_link_libraries(lsdtt-test-flow-d8-identity PRIVATE lsdtt3_cpu)
add_test(NAME unit_flow_d8_identity COMMAND lsdtt-test-flow-d8-identity)
lsdtt_set_tests_properties(unit_flow_d8_identity PROPERTIES LABELS "unit;flow")

add_test(NAME unit_flow_snapping COMMAND lsdtt-test-flow-snapping)
lsdtt_set_tests_properties(unit_flow_snapping PROPERTIES LABELS "unit;flow;snapping")

add_test(NAME unit_mountain_front_outlet_pruning
         COMMAND lsdtt-test-mountain-front-outlet-pruning)
lsdtt_set_tests_properties(unit_mountain_front_outlet_pruning PROPERTIES
                     LABELS "unit;flow;mountain_front")

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
    LABELS "unit;depression_hierarchy"
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
    LABELS "unit;depression_hierarchy")

# Internal drainage, phases 0, 1 and 3: the crater-field
# oracle, fill_pits_with_sinks, forced_outlet_flats, hole punching and sink
# tagging in both routing modes.
# Synthetic surfaces generated in process, so there is no fixture on disk.
add_executable(lsdtt-test-internal-drainage cpp/test_internal_drainage.cpp)
target_link_libraries(lsdtt-test-internal-drainage PRIVATE lsdtt3_cpu)

add_test(NAME unit_internal_drainage
         COMMAND lsdtt-test-internal-drainage)
lsdtt_set_tests_properties(unit_internal_drainage PROPERTIES
    LABELS "unit;internal_drainage")

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
        LABELS "smoke;cli;internal_drainage"
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
    LABELS "smoke;cli;internal_drainage"
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
    LABELS "unit;internal_drainage")
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
    LABELS "unit;point_input")
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/point_input_output)

# Characterisation tests for the src/utils + src/cpu consolidation (stage 1).
# Each embeds the private copies it replaced, verbatim, as the reference, and
# stays after those copies are deleted.
add_executable(lsdtt-test-raster-cell-predicates cpp/test_raster_cell_predicates.cpp)
target_link_libraries(lsdtt-test-raster-cell-predicates PRIVATE lsdtt3_cpu)
add_test(NAME unit_raster_cell_predicates COMMAND lsdtt-test-raster-cell-predicates)
lsdtt_set_tests_properties(unit_raster_cell_predicates PROPERTIES
    LABELS "unit;raster")

add_executable(lsdtt-test-distributed-hydrology-argument-checks
    cpp/test_distributed_hydrology_argument_checks.cpp)
target_link_libraries(lsdtt-test-distributed-hydrology-argument-checks PRIVATE lsdtt3_cpu)
add_test(NAME unit_distributed_hydrology_argument_checks
         COMMAND lsdtt-test-distributed-hydrology-argument-checks)
lsdtt_set_tests_properties(unit_distributed_hydrology_argument_checks PROPERTIES
    LABELS "unit;distributed_hydrology")

add_executable(lsdtt-test-drainage-injection-weights cpp/test_drainage_injection_weights.cpp)
target_link_libraries(lsdtt-test-drainage-injection-weights PRIVATE lsdtt3_cpu)
add_test(NAME unit_drainage_injection_weights
         COMMAND lsdtt-test-drainage-injection-weights
                 ${CMAKE_CURRENT_BINARY_DIR}/drainage_injection_weights_output)
lsdtt_set_tests_properties(unit_drainage_injection_weights PROPERTIES
    LABELS "unit;point_input;flow")
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/drainage_injection_weights_output)

# Stage 2: the in-place medians and quantiles in src/utils/lsdtt_statistics.
add_executable(lsdtt-test-order-statistics cpp/test_order_statistics.cpp)
target_link_libraries(lsdtt-test-order-statistics PRIVATE lsdtt3_cpu)
add_test(NAME unit_order_statistics COMMAND lsdtt-test-order-statistics)
lsdtt_set_tests_properties(unit_order_statistics PROPERTIES
    LABELS "unit;statistics")

add_executable(lsdtt-test-hydraulic-routing cpp/test_hydraulic_routing.cpp)
target_link_libraries(lsdtt-test-hydraulic-routing PRIVATE lsdtt3_cpu)

add_test(NAME unit_hydraulic_routing
         COMMAND lsdtt-test-hydraulic-routing)
lsdtt_set_tests_properties(unit_hydraulic_routing PROPERTIES LABELS "unit;hydraulic")

add_executable(lsdtt-test-raster-resampling cpp/test_raster_resampling.cpp)
target_link_libraries(lsdtt-test-raster-resampling PRIVATE lsdtt3_cpu)

add_test(NAME unit_raster_resampling
         COMMAND lsdtt-test-raster-resampling)
lsdtt_set_tests_properties(unit_raster_resampling PROPERTIES LABELS "unit;raster")

add_executable(lsdtt-test-raster-hierarchy cpp/test_raster_hierarchy.cpp)
target_link_libraries(lsdtt-test-raster-hierarchy PRIVATE lsdtt3_cpu)

add_test(NAME unit_raster_hierarchy
         COMMAND lsdtt-test-raster-hierarchy)
lsdtt_set_tests_properties(unit_raster_hierarchy PROPERTIES LABELS "unit;raster")

add_executable(lsdtt-test-raster-jenks-classification
               cpp/test_raster_jenks_classification.cpp)
target_link_libraries(lsdtt-test-raster-jenks-classification PRIVATE lsdtt3_cpu)

add_test(NAME unit_raster_jenks_classification
         COMMAND lsdtt-test-raster-jenks-classification)
lsdtt_set_tests_properties(unit_raster_jenks_classification PROPERTIES
    LABELS "unit;raster;jenks")

add_executable(lsdtt-test-raster-active-masks cpp/test_raster_active_masks.cpp)
target_link_libraries(lsdtt-test-raster-active-masks PRIVATE lsdtt3_cpu)

add_test(NAME unit_raster_active_masks COMMAND lsdtt-test-raster-active-masks)
lsdtt_set_tests_properties(unit_raster_active_masks PROPERTIES
    LABELS "unit;raster;mask")

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
    LABELS "unit;hillslope")

add_executable(lsdtt-test-critical-slope cpp/test_critical_slope.cpp)
target_link_libraries(lsdtt-test-critical-slope PRIVATE lsdtt3_cpu)

add_test(NAME unit_divide_diagnostics
         COMMAND lsdtt-test-divide-diagnostics)
lsdtt_set_tests_properties(unit_divide_diagnostics PROPERTIES
    LABELS "unit;divide_diagnostics")

add_test(NAME unit_paired_flowpaths
         COMMAND lsdtt-test-paired-flowpaths
                 "${LSDTT_PAIRED_FLOWPATHS_FIXTURE_DIR}")
lsdtt_set_tests_properties(unit_paired_flowpaths PROPERTIES
    LABELS "unit;paired_flowpaths"
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
    LABELS "smoke;cli;paired_flowpaths"
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
    LABELS "smoke;cli;paired_flowpaths"
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
    LABELS "smoke;cli;paired_flowpaths"
    FIXTURES_REQUIRED paired_flowpaths_fixture)

add_test(NAME unit_divide_migration
         COMMAND lsdtt-test-divide-migration
                 "${LSDTT_DIVIDE_MIGRATION_FIXTURE_DIR}")
lsdtt_set_tests_properties(unit_divide_migration PROPERTIES
    LABELS "unit;divide_migration"
    FIXTURES_SETUP divide_migration_fixture)

add_test(NAME unit_landscape_equilibration
         COMMAND lsdtt-test-landscape-equilibration)
lsdtt_set_tests_properties(unit_landscape_equilibration PROPERTIES
    LABELS "unit;landscape_equilibration")

add_test(NAME unit_hillslope_critical_slope
         COMMAND lsdtt-test-critical-slope)
lsdtt_set_tests_properties(unit_hillslope_critical_slope PROPERTIES
    LABELS "unit;hillslope")

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
    LABELS "unit;streamlines")

# Ridge extraction for hilltop flow routing (Phase 1). The divide-grid corner
# convention is the one claim everything downstream rests on, so the fixtures
# are hand-drawn 5x5 grids with the expected pixels written out explicitly,
# including a periodic-seam case and the non-4-connected step.
add_executable(lsdtt-test-hilltop-flow-routing cpp/test_hilltop_flow_routing.cpp)
target_link_libraries(lsdtt-test-hilltop-flow-routing PRIVATE lsdtt3_cpu)

add_test(NAME unit_hilltop_flow_routing
         COMMAND lsdtt-test-hilltop-flow-routing)
lsdtt_set_tests_properties(unit_hilltop_flow_routing PROPERTIES
    LABELS "unit;hilltop_flow_routing")

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
    LABELS "unit;hilltop_flow_routing")

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
    set(lbl "cli;hilltop_flow_routing")
    if(TEST_NAME MATCHES "^failure_")
        list(APPEND lbl "failure")
    else()
        list(APPEND lbl "smoke")
    endif()
    lsdtt_set_tests_properties(${TEST_NAME} PROPERTIES
        LABELS "${lbl}"
        FIXTURES_REQUIRED hilltop_coupling_cli_fixture)
endfunction()

add_hilltop_coupling_cli_test(
    smoke_hillslope_channel_coupling
    hillslope_channel_coupling_smoke.param.in
    hcc_hilltop_points.csv
    --expect hcc_ridgelines.fgb
    --stdout-contains "never traced (DEM border)")

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
    LABELS "unit;lem")

# Pathway A: the channel
# threshold search, and the ways its objective can be plausibly wrong.
add_executable(lsdtt-test-hillslope-length-search cpp/test_hillslope_length_search.cpp)
target_link_libraries(lsdtt-test-hillslope-length-search PRIVATE lsdtt3_cpu)

add_test(NAME unit_hillslope_length_search
         COMMAND lsdtt-test-hillslope-length-search)
lsdtt_set_tests_properties(unit_hillslope_length_search PROPERTIES
    LABELS "unit;hillslope")

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
    LABELS "smoke;cli;raster" TIMEOUT 1200)

# Longitudinal profile from the points output.
add_executable(lsdtt-test-channel-profile-attributes cpp/test_channel_profile_attributes.cpp)
target_link_libraries(lsdtt-test-channel-profile-attributes PRIVATE lsdtt3_cpu)

file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/profile_output)
add_test(NAME unit_channel_profile_attributes
         COMMAND lsdtt-test-channel-profile-attributes
                 ${CMAKE_CURRENT_BINARY_DIR}/profile_output)
lsdtt_set_tests_properties(unit_channel_profile_attributes PROPERTIES LABELS "unit;network")

# link_length under periodic boundaries.
add_executable(lsdtt-test-link-length-periodic cpp/test_link_length_periodic.cpp)
target_link_libraries(lsdtt-test-link-length-periodic PRIVATE lsdtt3_cpu)

file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/link_length_output)
add_test(NAME unit_link_length_periodic
         COMMAND lsdtt-test-link-length-periodic
                 ${CMAKE_CURRENT_BINARY_DIR}/link_length_output)
lsdtt_set_tests_properties(unit_link_length_periodic PROPERTIES LABELS "unit;network")

# Multipart geometry on read.
add_executable(lsdtt-test-multipart-geometry cpp/test_multipart_geometry.cpp)
target_link_libraries(lsdtt-test-multipart-geometry PRIVATE lsdtt3_core)

file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/multipart_output)
add_test(NAME unit_multipart_geometry
         COMMAND lsdtt-test-multipart-geometry
                 ${CMAKE_CURRENT_BINARY_DIR}/multipart_output)
lsdtt_set_tests_properties(unit_multipart_geometry PROPERTIES LABELS "unit;raster")

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
lsdtt_set_tests_properties(unit_hillslope_peclet PROPERTIES LABELS "unit;hillslope")

# Boundary conditions for the explicit Roering transport solver (periodic + 'o').
add_executable(lsdtt-test-hillslope-transport-bc cpp/test_hillslope_transport_bc.cpp)
target_link_libraries(lsdtt-test-hillslope-transport-bc PRIVATE lsdtt3_cpu)

add_test(NAME unit_hillslope_transport_bc COMMAND lsdtt-test-hillslope-transport-bc)
lsdtt_set_tests_properties(unit_hillslope_transport_bc PROPERTIES LABELS "unit;hillslope")

# Coupled-LEM timestep control: an accuracy limit that no stability check catches.
add_executable(lsdtt-test-lem-timestep cpp/test_lem_timestep.cpp)
target_link_libraries(lsdtt-test-lem-timestep PRIVATE lsdtt3_cpu)

add_test(NAME unit_lem_timestep COMMAND lsdtt-test-lem-timestep)
lsdtt_set_tests_properties(unit_lem_timestep PROPERTIES LABELS "unit;lem")

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
lsdtt_set_tests_properties(unit_lem_forcing_history PROPERTIES LABELS "unit;lem")
lsdtt_set_tests_properties(unit_lem_channel_coupling PROPERTIES LABELS "unit;lem;hillslope")

# The Newton solver for the Roering nonlinear hillslope (steady + transient).
add_executable(lsdtt-test-hillslope-newton cpp/test_hillslope_newton.cpp)
target_link_libraries(lsdtt-test-hillslope-newton PRIVATE lsdtt3_cpu)

add_test(NAME unit_hillslope_newton COMMAND lsdtt-test-hillslope-newton)
lsdtt_set_tests_properties(unit_hillslope_newton PROPERTIES LABELS "unit;hillslope")

# Independent equation-faithful Perron (2011) Q-imp reference solver.
add_executable(lsdtt-test-perron-qimp cpp/test_perron_qimp.cpp)
target_link_libraries(lsdtt-test-perron-qimp PRIVATE lsdtt3_cpu)

add_test(NAME unit_hillslope_perron_qimp COMMAND lsdtt-test-perron-qimp)
lsdtt_set_tests_properties(unit_hillslope_perron_qimp PROPERTIES LABELS "unit;hillslope")

# Independent C++ port of Ren et al. (2026) modular implicit method.
add_executable(lsdtt-test-ren-modular-implicit cpp/test_ren_modular_implicit.cpp)
target_link_libraries(lsdtt-test-ren-modular-implicit PRIVATE lsdtt3_cpu)

add_test(NAME unit_hillslope_ren_modular_implicit
         COMMAND lsdtt-test-ren-modular-implicit)
lsdtt_set_tests_properties(unit_hillslope_ren_modular_implicit PROPERTIES
                     LABELS "unit;hillslope")

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
                     LABELS "unit;hillslope;cellular")
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/cellular_hillslope_reference_output)

# Soft-critical constitutive law and whole-grid pass-cost envelope. This is P1:
# test-local mathematics only, before a raster transport solver exists.
add_executable(lsdtt-test-cellular-hillslope-constitutive
               cpp/test_cellular_hillslope_constitutive.cpp)

add_test(NAME unit_cellular_hillslope_constitutive
         COMMAND lsdtt-test-cellular-hillslope-constitutive
                 ${CMAKE_CURRENT_BINARY_DIR}/cellular_hillslope_reference_output)
lsdtt_set_tests_properties(unit_cellular_hillslope_constitutive PROPERTIES
                     LABELS "unit;hillslope;cellular")

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
                     LABELS "unit;hillslope;cellular")

add_test(NAME unit_cellular_hillslope_common_limit
         COMMAND cellular-hillslope-common-limit)
lsdtt_set_tests_properties(unit_cellular_hillslope_common_limit PROPERTIES
                     LABELS "unit;hillslope;cellular")

add_test(NAME unit_cellular_hillslope_solver
         COMMAND lsdtt-test-cellular-hillslope-solver)
lsdtt_set_tests_properties(unit_cellular_hillslope_solver PROPERTIES
                     LABELS "unit;hillslope;cellular")

# The SOR steady-state hillslope solver: conservation, precision floor, 'o' vs 'n',
# periodic invariance, and recovery of U/D from known forcing.
add_executable(lsdtt-test-hillslope-steady-state cpp/test_hillslope_steady_state.cpp)
target_link_libraries(lsdtt-test-hillslope-steady-state PRIVATE lsdtt3_cpu)

add_test(NAME unit_hillslope_steady_state COMMAND lsdtt-test-hillslope-steady-state)
lsdtt_set_tests_properties(unit_hillslope_steady_state PROPERTIES LABELS "unit;hillslope")

# chi at the channel head, as distinct from domain-wide MaxChi.
add_executable(lsdtt-test-chi-channel-head cpp/test_chi_channel_head.cpp)
target_link_libraries(lsdtt-test-chi-channel-head PRIVATE lsdtt3_cpu)

add_test(NAME unit_chi_channel_head COMMAND lsdtt-test-chi-channel-head)
lsdtt_set_tests_properties(unit_chi_channel_head PROPERTIES LABELS "unit;lem")

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
lsdtt_set_tests_properties(unit_hillslope_implicit PROPERTIES LABELS "unit;hillslope")
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
add_test(NAME smoke_fetch_raster_help
         COMMAND lsdtt-fetch-raster -h)
lsdtt_set_tests_properties(smoke_fetch_raster_help PROPERTIES
    LABELS "smoke;cli;fetch"
    PASS_REGULAR_EXPRESSION "Bounding box")

# The CLI must never grow a credential parameter: parameter files get committed
# and print_parameters() echoes every parameter at the end of a run. This guards
# the rule stated in docs/developer_notes/testing_and_benchmarking.md and in the CLI's own header.
add_test(NAME rigorous_fetch_raster_no_credential_parameters
         COMMAND lsdtt-fetch-raster -h)
lsdtt_set_tests_properties(rigorous_fetch_raster_no_credential_parameters PROPERTIES
    LABELS "rigorous;cli;fetch"
    FAIL_REGULAR_EXPRESSION "api_key;API_Key;password;secret;token;credential")

add_lsdtt_cli_smoke_test(
    failure_fill_vincent_needs_zero_slope
    lsdtt-raster-preprocessing
    fill_vincent_needs_zero_slope.param.in
    unused_expected_output
    --expect-failure
    --stderr-contains "cannot impose a gradient across flats"
)

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
        LABELS "regression;cli;parameters"
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
    LABELS "unit;cli;parameters")

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
    LABELS "unit;cli;synthetic_networks;lem")

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
        PROPERTIES LABELS "performance;flow" TIMEOUT 1800)

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
        PROPERTIES LABELS "performance;hillslope" SKIP_RETURN_CODE 77 TIMEOUT 1800)

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
            PROPERTIES LABELS "performance;surface_metrics" SKIP_RETURN_CODE 77
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
        PROPERTIES LABELS "performance;cosmo" TIMEOUT 1800)

    # HYPE model timing: set-up, per-step
    # and per-unit-step cost, phase split, peak RSS on synthetic in-memory
    # workloads. REPORTED, NEVER ASSERTED -- fails only on a water-balance or
    # non-finite-flow error. --quick runs the small workload only; run the
    # binary directly, without --quick, for the three-size baseline.
    add_test(NAME performance_hype_model
             COMMAND lsdtt-test-hype-performance
                     ${CMAKE_CURRENT_BINARY_DIR}/perf_output --quick)
    lsdtt_set_tests_properties(performance_hype_model
        PROPERTIES LABELS "performance;distributed_hydrology" TIMEOUT 1800)

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
lsdtt_set_tests_properties(unit_secret PROPERTIES LABELS "unit;fetch")
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
lsdtt_set_tests_properties(unit_test_env_isolation PROPERTIES LABELS "unit;core")

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
        LABELS "unit;configure_contract")
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
        LABELS "unit;configure_contract")
endforeach()

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
lsdtt_set_tests_properties(unit_raster_fetch PROPERTIES LABELS "unit;fetch")

add_executable(lsdtt-test-dem-grouping cpp/test_dem_grouping.cpp)
target_link_libraries(lsdtt-test-dem-grouping PRIVATE lsdtt3_core)
target_include_directories(lsdtt-test-dem-grouping PRIVATE ${CMAKE_SOURCE_DIR}/src/core)

add_test(NAME unit_dem_grouping
         COMMAND lsdtt-test-dem-grouping "${CMAKE_CURRENT_BINARY_DIR}/dem_grouping_output")
lsdtt_set_tests_properties(unit_dem_grouping PROPERTIES LABELS "unit;fetch")

add_executable(lsdtt-test-coastal-mask cpp/test_coastal_mask.cpp)
target_link_libraries(lsdtt-test-coastal-mask PRIVATE lsdtt3_core)
target_include_directories(lsdtt-test-coastal-mask PRIVATE ${CMAKE_SOURCE_DIR}/src/core)

add_test(NAME unit_coastal_mask
         COMMAND lsdtt-test-coastal-mask "${CMAKE_CURRENT_BINARY_DIR}/coastal_mask_output")
lsdtt_set_tests_properties(unit_coastal_mask PROPERTIES LABELS "unit;fetch")

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
lsdtt_set_tests_properties(unit_step_change PROPERTIES LABELS "unit;step_change")
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/step_change_test_output)

# Landslide size distribution T1-T5b. Emits landslide_truncation.csv, which
# tabulates how much of the distribution the A_max cutoff removes by number and
# by volume.
add_test(NAME unit_landslide_size
         COMMAND lsdtt-test-landslide-size
                 ${CMAKE_CURRENT_BINARY_DIR}/landslide_test_output)
lsdtt_set_tests_properties(unit_landslide_size PROPERTIES LABELS "unit;landslide")
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
    LABELS "unit;cosmo;landslide" TIMEOUT 1800)

# Quartz fertility Q1-Q6 and the ice mask I1-I3. Pins the case that makes
# fertility a correctness issue rather than a refinement: a quartz-free unit
# eroding ten times faster than its neighbour must contribute nothing.
add_test(NAME unit_quartz_fertility COMMAND lsdtt-test-quartz-fertility)
lsdtt_set_tests_properties(unit_quartz_fertility PROPERTIES LABELS "unit;cosmo")

add_test(NAME unit_raster_transform COMMAND lsdtt-test-raster-transform)
lsdtt_set_tests_properties(unit_raster_transform PROPERTIES LABELS "unit;raster_algebra;cosmo")

# Landslide scar coupling and densities. Emits landslide_scar_recovery.csv,
# which shows how a scar's concentration recovers with age at three depths and
# what carrying it at the wrong material density does to the answer.
add_test(NAME unit_landslide_step_change
         COMMAND lsdtt-test-landslide-step-change
                 ${CMAKE_CURRENT_BINARY_DIR}/landslide_test_output)
lsdtt_set_tests_properties(unit_landslide_step_change PROPERTIES
    LABELS "unit;landslide;step_change")

# Landslide placement T6-T10d. Emits landslide_host_fill.csv, which records how
# completely landslides of each size fill their host basin.
add_test(NAME unit_landslide_placement
         COMMAND lsdtt-test-landslide-placement
                 ${CMAKE_CURRENT_BINARY_DIR}/landslide_test_output)
lsdtt_set_tests_properties(unit_landslide_placement PROPERTIES LABELS "unit;landslide")

# Landslide event T11-T13b. Emits landslide_budget_vs_amax.csv.
add_test(NAME unit_landslide_event
         COMMAND lsdtt-test-landslide-event
                 ${CMAKE_CURRENT_BINARY_DIR}/landslide_test_output)
lsdtt_set_tests_properties(unit_landslide_event PROPERTIES LABELS "unit;landslide")

# Landslide CRN response T16. Emits landslide_crn_response.csv.
add_test(NAME unit_landslide_crn
         COMMAND lsdtt-test-landslide-crn
                 ${CMAKE_CURRENT_BINARY_DIR}/landslide_test_output)
lsdtt_set_tests_properties(unit_landslide_crn PROPERTIES LABELS "unit;landslide;cosmo")

# Landslide series: Phase 3 multi-event + Phase 5 T14/T15. Emits
# landslide_niemi_closure.csv and landslide_yanites_regimes.csv.
add_test(NAME unit_landslide_series
         COMMAND lsdtt-test-landslide-series
                 ${CMAKE_CURRENT_BINARY_DIR}/landslide_test_output)
lsdtt_set_tests_properties(unit_landslide_series PROPERTIES LABELS "unit;landslide" TIMEOUT 900)

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
lsdtt_set_tests_properties(unit_cosmo_units PROPERTIES LABELS "unit;cosmo")

add_test(NAME unit_cosmo_ratio COMMAND lsdtt-test-cosmo-ratio)
lsdtt_set_tests_properties(unit_cosmo_ratio PROPERTIES LABELS "unit;cosmo")

add_test(NAME unit_cosmo_physics
         COMMAND lsdtt-test-cosmo-physics
                 ${CMAKE_CURRENT_BINARY_DIR}/cosmo_test_output
                 ${CMAKE_CURRENT_SOURCE_DIR}/expected/cosmo/numerics_sentinel.csv
                 ${CMAKE_CURRENT_SOURCE_DIR}/expected/cosmo/numerics_sentinel_lsdn.csv
                 ${CMAKE_CURRENT_SOURCE_DIR}/data)
lsdtt_set_tests_properties(unit_cosmo_physics PROPERTIES LABELS "unit;cosmo")
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
    PROPERTIES LABELS "failure;cli;cosmo")

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
    PROPERTIES LABELS "smoke;cli;cosmo")

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
    PROPERTIES LABELS "smoke;cli;cosmo")

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
    PROPERTIES LABELS "smoke;cli;cosmo")

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
    PROPERTIES LABELS "failure;cli;cosmo")

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
    LABELS "rigorous;cli;cosmo" SKIP_RETURN_CODE 77 TIMEOUT 900)

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
    LABELS "rigorous;cli;cosmo" SKIP_RETURN_CODE 77 TIMEOUT 900)

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
    PROPERTIES LABELS "smoke;cli;cosmo")

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
    PROPERTIES LABELS "smoke;cli;cosmo")

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
    PROPERTIES LABELS "smoke;cli;cosmo")

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
    PROPERTIES LABELS "smoke;cli;cosmo")

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
    PROPERTIES LABELS "failure;cli;cosmo")

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
    PROPERTIES LABELS "failure;cli;cosmo")

# Help must render every parameter without a parameter file and without
# touching the DEM or the network. Guards that the parameter declarations for
# each cosmo CLI stay well-formed.
# lsdtt-cosmo-landslides.
add_test(NAME smoke_cosmo_landslides_help COMMAND lsdtt-cosmo-landslides -h)
lsdtt_set_tests_properties(smoke_cosmo_landslides_help PROPERTIES
    LABELS "smoke;cli;landslide"
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
    LABELS "smoke;cli;cosmo"
    PASS_REGULAR_EXPRESSION "sample_age_yr")

add_test(NAME smoke_divide_migration_help COMMAND lsdtt-divide-migration -h)
lsdtt_set_tests_properties(smoke_divide_migration_help PROPERTIES
    LABELS "smoke;cli"
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
    LABELS "cli;cosmo;regression")

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
    LABELS "cli;cosmo;regression;heavy")

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
    LABELS "smoke;cli;cosmo"
    PASS_REGULAR_EXPRESSION "sample_age_yr")

add_test(NAME regression_cosmo_prediction_builtin_io_help
         COMMAND lsdtt-cosmo-prediction -h read_path)
lsdtt_set_tests_properties(regression_cosmo_prediction_builtin_io_help PROPERTIES
    LABELS "regression;cli;cosmo;parameters"
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
    LABELS "rigorous;cli;cosmo" SKIP_RETURN_CODE 77 TIMEOUT 600)

add_test(NAME smoke_cosmo_column_tool_help COMMAND lsdtt-cosmo-column-tool -h)
lsdtt_set_tests_properties(smoke_cosmo_column_tool_help PROPERTIES
    LABELS "smoke;cli;cosmo"
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
lsdtt_set_tests_properties(unit_cosmo_validation_harness PROPERTIES LABELS "unit;cosmo")

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
lsdtt_set_tests_properties(unit_cosmo_riversand_harness PROPERTIES LABELS "unit;cosmo")

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

add_executable(lsdtt-test-fluvial-surfaces cpp/test_fluvial_surfaces.cpp)
target_link_libraries(lsdtt-test-fluvial-surfaces PRIVATE lsdtt3_cpu)

add_test(NAME unit_fluvial_surfaces_fixture
         COMMAND lsdtt-test-fluvial-surfaces)
lsdtt_set_tests_properties(unit_fluvial_surfaces_fixture PROPERTIES
    LABELS "unit")

# Phase 0c. Tests T2, T4, T5 and T7: the resolvability floor
# reports UNRESOLVED rather than fabricating levels, labels are a globally
# unique running counter whose band map is single-valued, k = 2 reproduces the
# legacy valley cut bit-for-bit, and k selection is a constraint with an
# extremum rather than a trade-off.
add_test(NAME unit_fluvial_surface_ladder
         COMMAND lsdtt-test-fluvial-surfaces --ladder)
lsdtt_set_tests_properties(unit_fluvial_surface_ladder PROPERTIES
    LABELS "unit")

# Phase 2. T3: the fitted-surface reference field. Both arms run the SAME ladder
# with the SAME auto k selection over the SAME admitted set, and only the banded
# raster differs, so what the test measures is the field rather than the tuning.
# HAND shreds the fan into radial arcs; the residual to a fitted cone does not.
add_test(NAME unit_fluvial_surface_fitted_field
         COMMAND lsdtt-test-fluvial-surfaces --fitted-surface)
lsdtt_set_tests_properties(unit_fluvial_surface_fitted_field PROPERTIES
    LABELS "unit")

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
    LABELS "unit;fluvial_surfaces")

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
    LABELS "unit;fluvial_surfaces")

add_test(NAME unit_fluvial_surface_slope_dependent_ev
         COMMAND lsdtt-test-fluvial-surfaces --slope-dependent-ev)
lsdtt_set_tests_properties(unit_fluvial_surface_slope_dependent_ev PROPERTIES
    LABELS "unit;fluvial_surfaces")

add_test(NAME unit_fluvial_surface_coherence_diffusion
         COMMAND lsdtt-test-fluvial-surfaces --coherence-diffusion)
lsdtt_set_tests_properties(unit_fluvial_surface_coherence_diffusion PROPERTIES
    LABELS "unit;fluvial_surfaces")

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
    LABELS "smoke;cli"
    PASS_REGULAR_EXPRESSION "relief_threshold_metres")

add_test(NAME smoke_valley_centreline_help
         COMMAND lsdtt-valley-centreline -h)
lsdtt_set_tests_properties(smoke_valley_centreline_help PROPERTIES
    LABELS "smoke;cli"
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
    --expect-file-contains "btl_ladder_windows.csv:source_id,window_centre_m,n_pixels,k_selected,status,used_global_breaks,breaks_m,band_centres_m,min_band_centre_gap_m,required_gap_m,k_bootstrap_mode,k_stability,k_ci_lo,k_ci_hi,n_bootstrap,vertical_accuracy_m,median_slope_degrees,goodness_of_variance_fit,separation_ratio"
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
    LABELS "rigorous;cli"
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
    LABELS "rigorous;cli"
    FIXTURES_REQUIRED fluvial_surfaces_fixture)

# Phase 7a. `admission_mode` exists because local relief and a gradient cap
# reject different ground; if a refactor ever collapsed them into one test, or
# wired relief_or_slope back to the relief mask alone, every count would still
# look plausible. The synthetic DEM is built to make the two disagree - a smooth
# 8 percent ramp that only the cap admits, a rough flat that only relief could,
# a smooth flat that both take and a 45 degree slope that neither does - and the
# assertions are on the SET RELATIONS between the four modes, which survive the
# morphological post-processing that the raw counts do not.
if(LSDTT_ENABLE_PYTHON_GEO_TESTS)
    add_test(NAME rigorous_flat_surface_finder_admission_modes
        COMMAND ${Python3_EXECUTABLE}
            "${CMAKE_CURRENT_SOURCE_DIR}/scripts/check_admission_modes.py"
            --exe "$<TARGET_FILE:lsdtt-flat-surface-finder>"
            --work-dir "${CMAKE_CURRENT_BINARY_DIR}/rigorous_flat_surface_finder_admission_modes"
            --cap 15.0)
    lsdtt_mark_test_runtime(python rigorous_flat_surface_finder_admission_modes)
else()
    lsdtt_mark_test_unregistered(python_geo rigorous_flat_surface_finder_admission_modes)
endif()
lsdtt_set_tests_properties(rigorous_flat_surface_finder_admission_modes PROPERTIES
    LABELS "rigorous;cli;fluvial_surfaces")

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
lsdtt_set_tests_properties(unit_fluvial_surfaces_trough_upstream PROPERTIES
    LABELS "unit")

add_test(NAME unit_fluvial_surfaces_band_then_label
         COMMAND lsdtt-test-fluvial-surfaces --band-then-label)
lsdtt_set_tests_properties(unit_fluvial_surfaces_band_then_label PROPERTIES
    LABELS "unit")

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
    LABELS "unit;fluvial_surfaces")

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
    LABELS "smoke;cli"
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
        LABELS "smoke;cli"
        PASS_REGULAR_EXPRESSION "threshold_contributing_pixels")
endforeach()

add_test(NAME smoke_feature_finder_fill_help
         COMMAND lsdtt-feature-finder -h)
lsdtt_set_tests_properties(smoke_feature_finder_fill_help PROPERTIES
    LABELS "smoke;cli"
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
    --expect-file-contains "fsurf_ladder_windows.csv:source_id,window_centre_m,n_pixels,k_selected,status,used_global_breaks,breaks_m,band_centres_m,min_band_centre_gap_m,required_gap_m,k_bootstrap_mode,k_stability,k_ci_lo,k_ci_hi,n_bootstrap,vertical_accuracy_m,median_slope_degrees,goodness_of_variance_fit,separation_ratio"
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
    LABELS "smoke;cli"
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

add_lsdtt_cli_smoke_test(
    failure_fluvial_surfaces_bad_reference_field
    lsdtt-fluvial-surfaces
    fluvial_surfaces_bad_reference_field.param.in
    unused.tif
    --read-path "${LSDTT_FLUVIAL_SURFACES_FIXTURE_DIR}"
    --expect-failure
    --stderr-contains "unknown reference_field_algorithm_choice 'elevation'. Valid: hand, fitted_surface."
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
    LABELS "smoke;cli"
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
    LABELS "rigorous;cli"
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
    LABELS "unit;cli;conventions"
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
    LABELS "unit;cli;docs"
    TIMEOUT 900)

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
    # Ahead of the "channel" rule below on purpose: CLI 2 of the fluvial-surface
    # split writes a channel network, so tests named for that would otherwise be
    # claimed by "network" and disappear from a fluvial_surfaces compartment run.
    "fluvial_surface|flat_surface|valley_centreline" "fluvial_surfaces"
    # The shared Topographic Wetness Index,
    # pointed at "network" ON PURPOSE rather than given its own domain: the
    # existing TWI/STI tests already live in the network compartment, and
    # splitting them off now would quietly remove them from `ctest -L network`
    # for anyone working on drainage. When lsdtt-distributed-hydrology lands and
    # there is a real cluster of tests on this module, flip this one line to
    # "wetness_index" and move the STI tests with it, deliberately.
    "wetness_index|hydro_velocity|hydro_roughness|hydro_routing|topmodel|canopy|channel_extent|infiltration|landcover|distributed_hydrology|hype"          "network"
    "chi|channel|drainage|basin|flow_properties|flow_snapping|flow_d8|link_length"  "network"
    # performance_flow_routines (LSDTT_ENABLE_PERF_TESTS) matched no rule, so
    # configuring with that option ON failed until 2026-09-18 -- the same trap
    # as the benchmark rule at the top, hidden by the option's OFF default.
    "flow_routines"          "network"
    "raster|nodata|fill|tif|coastal|dem_grouping|vertical_units|multipart" "raster"
    "param|secret|cli_reference|test_env|configure_contract" "core"
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
    "mountain_front|feature_finder" "surface_metrics"
    # Added 2026-10-07 with the first tests these four CLIs have had. Cross
    # sections, bankfull width and the water-body centreline are channel /
    # hydrography work; marsh platform and tidal channel extraction is landform
    # extraction, beside feature_finder.
    "cross_sections|bankfull_width|water_centreline" "network"
    "marsh_topography"       "surface_metrics"
    "basin_stats"            "network"
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
)

# Tests that use every core they can get. Under `ctest -j` these must not share
# the machine: doing so does not merely slow them, it starves them past their
# timeout, which is how unit_hillslope_implicit came to "fail" at 1500 s while
# taking 9 minutes when run alone.
set(LSDTT_SERIAL_TESTS
    unit_hillslope_implicit
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
    LABELS "unit;core;conventions"
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
    LABELS "unit;core;conventions"
    TIMEOUT 120)

# The test commands the documentation tells people to type: every -DLSDTT_*
# option it names exists, every anchored -L/-LE label is one this file can
# assign, and the test_report target and run_tests.cmake it names exist
# (cmake/check_doc_test_commands.cmake). Reads source files only.
add_test(NAME unit_docs_test_commands
    COMMAND ${CMAKE_COMMAND} -DSOURCE_DIR=${CMAKE_SOURCE_DIR}
        -P "${CMAKE_CURRENT_SOURCE_DIR}/cmake/check_doc_test_commands.cmake")
set_tests_properties(unit_docs_test_commands PROPERTIES
    LABELS "unit;core;docs"
    TIMEOUT 60)

add_test(NAME unit_test_inventory_audit_selftest
    COMMAND ${CMAKE_COMMAND}
        -DSELFTEST_DIR=${CMAKE_CURRENT_SOURCE_DIR}/cmake/audit_selftest
        -P "${CMAKE_CURRENT_SOURCE_DIR}/cmake/audit_test_inventory.cmake")
set_tests_properties(unit_test_inventory_audit_selftest PROPERTIES
    LABELS "unit;core;conventions"
    TIMEOUT 60)

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
    LABELS "unit;core;conventions"
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
# These are ORTHOGONAL to the kind and domain labels above. Select with
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
