# configure_contract.cmake
#
# Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
#
# Copyright (C) 2026 Simon M. Mudd
# SPDX-License-Identifier: AGPL-3.0-or-later
#
# Configures the WHOLE project in a scratch directory for one named case and
# checks what configure says: that each disabled or missing test tier produces
# the documented message (or error) rather than a silent partial suite.
#
#   cmake -DCASE=<case> -DSOURCE_DIR=<src> -DWORK_DIR=<scratch>
#         -DGENERATOR=<gen> -DCXX_COMPILER=<c++> [-DPYTHON=<python>]
#         -P configure_contract.cmake
#
# Cases:
#   native_default         defaults: configures, warns that the Python tests are
#                          not registered, and the summary shows Python and the
#                          network tier as off
#   geo_without_python     geo ON, Python OFF: a configure error naming
#                          LSDTT_ENABLE_PYTHON_TESTS
#   geo_missing_modules    geo ON with a bare venv of PYTHON: a configure error
#                          naming the imports that fail
#   gdal_programs_missing  Python ON with the GDAL program names overridden to
#                          one that does not exist: a configure error naming it
#                          and gdal-bin
#
# Each case is a full fresh configure (~20 s), so these tests carry the label
# configure_contract.

cmake_minimum_required(VERSION 3.15)

foreach(_req CASE SOURCE_DIR WORK_DIR GENERATOR CXX_COMPILER)
    if(NOT DEFINED ${_req})
        message(FATAL_ERROR "configure_contract: ${_req} not set")
    endif()
endforeach()

file(REMOVE_RECURSE "${WORK_DIR}")
file(MAKE_DIRECTORY "${WORK_DIR}")
set(_build "${WORK_DIR}/build")

set(_opts -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=ON)
set(_expect_rc 0)
set(_expect "")

if(CASE STREQUAL "native_default")
    list(APPEND _expect
         "Python-driven tests are not registered"
         "tier_python_stdlib : disabled"
         "needs_network      : disabled"
         "Python             : not used")

elseif(CASE STREQUAL "geo_without_python")
    list(APPEND _opts -DLSDTT_ENABLE_PYTHON_TESTS=OFF -DLSDTT_ENABLE_PYTHON_GEO_TESTS=ON)
    set(_expect_rc 1)
    list(APPEND _expect "needs LSDTT_ENABLE_PYTHON_TESTS=ON")

elseif(CASE STREQUAL "geo_missing_modules")
    if(NOT DEFINED PYTHON)
        message(FATAL_ERROR "configure_contract: geo_missing_modules needs -DPYTHON")
    endif()
    # A venv without pip or site-packages: the standard library and nothing else.
    execute_process(COMMAND "${PYTHON}" -m venv --without-pip "${WORK_DIR}/venv"
                    RESULT_VARIABLE _rc OUTPUT_VARIABLE _o ERROR_VARIABLE _o)
    if(NOT _rc EQUAL 0)
        message(FATAL_ERROR "configure_contract: could not create a venv: ${_o}")
    endif()
    if(WIN32)
        set(_bare "${WORK_DIR}/venv/Scripts/python.exe")
    else()
        set(_bare "${WORK_DIR}/venv/bin/python")
    endif()
    list(APPEND _opts -DLSDTT_ENABLE_PYTHON_TESTS=ON -DLSDTT_ENABLE_PYTHON_GEO_TESTS=ON
                      "-DPython3_EXECUTABLE=${_bare}")
    set(_expect_rc 1)
    list(APPEND _expect "cannot import: osgeo.gdal" "pixi install -e test")

elseif(CASE STREQUAL "gdal_programs_missing")
    if(NOT DEFINED PYTHON)
        message(FATAL_ERROR "configure_contract: gdal_programs_missing needs -DPYTHON")
    endif()
    list(APPEND _opts -DLSDTT_ENABLE_PYTHON_TESTS=ON -DLSDTT_ENABLE_PYTHON_GEO_TESTS=OFF
                      "-DPython3_EXECUTABLE=${PYTHON}"
                      "-DLSDTT_TEST_GDAL_PROGRAMS=lsdtt_no_such_gdal_program")
    set(_expect_rc 1)
    list(APPEND _expect "GDAL's command-line programs are" "lsdtt_no_such_gdal_program" "gdal-bin")

else()
    message(FATAL_ERROR "configure_contract: unknown CASE '${CASE}'")
endif()

execute_process(
    COMMAND "${CMAKE_COMMAND}" -S "${SOURCE_DIR}" -B "${_build}"
            -G "${GENERATOR}" "-DCMAKE_CXX_COMPILER=${CXX_COMPILER}" ${_opts}
    RESULT_VARIABLE _rc
    OUTPUT_VARIABLE _out ERROR_VARIABLE _out)
file(WRITE "${WORK_DIR}/configure_output.txt" "${_out}")

# CMake's own line wrapping can split a message; compare on collapsed whitespace.
string(REGEX REPLACE "[ \t\r\n]+" " " _flat "${_out}")

set(_fail "")
if(_expect_rc EQUAL 0 AND NOT _rc EQUAL 0)
    list(APPEND _fail "configure failed (exit ${_rc}) but should have succeeded")
elseif(NOT _expect_rc EQUAL 0 AND _rc EQUAL 0)
    list(APPEND _fail "configure succeeded but should have failed")
endif()
foreach(_want IN LISTS _expect)
    string(REGEX REPLACE "[ \t\r\n]+" " " _w "${_want}")
    string(FIND "${_flat}" "${_w}" _pos)
    if(_pos LESS 0)
        list(APPEND _fail "configure output lacks: ${_want}")
    endif()
endforeach()

if(_fail)
    string(REPLACE ";" "\n  " _fail "${_fail}")
    message(FATAL_ERROR "configure contract '${CASE}' broken:\n  ${_fail}\n"
                        "(full output: ${WORK_DIR}/configure_output.txt)")
endif()
list(LENGTH _expect _n_expect)
message("configure contract '${CASE}': exit ${_rc} as expected; ${_n_expect} expected message(s) present")
file(REMOVE_RECURSE "${_build}")
