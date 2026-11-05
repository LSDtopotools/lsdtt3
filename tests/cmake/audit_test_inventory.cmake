# audit_test_inventory.cmake
#
# Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
#
# Copyright (C) 2026 Simon M. Mudd
# SPDX-License-Identifier: AGPL-3.0-or-later
#
# Checks every registered CTest test's labels and properties against the
# command it actually runs.
#
#   cmake -DTESTFILE=<build>/tests/CTestTestfile.cmake
#         -DSCRIPTS_DIR=<source>/tests/scripts
#         [-DCTEST_CONFIGURATION_TYPE=<config>]
#         -P audit_test_inventory.cmake
#
#   cmake -DSELFTEST_DIR=<source>/tests/cmake/audit_selftest -P audit_test_inventory.cmake
#
# WHY THIS READS CTestTestfile.cmake AND DOES NOT CALL `ctest --show-only`
#
# The obvious implementation runs `ctest --show-only=json-v1`. But this runs AS
# a test, and any `ctest` invocation inside a running suite rewrites
# Testing/Temporary/LastTest.log(.tmp): the suite's own log comes out empty
# (checked 2026-09-18). So instead the generated CTestTestfile.cmake -- which is
# itself CMake code, one add_test() and one set_tests_properties() per test,
# with generator expressions already resolved -- is include()d with those two
# commands redefined to record what they are given. Nothing is written.
#
# WHAT IS CHECKED (each violation is reported with its check letter)
#
#   A  every test has exactly one runtime tier label (tier_native,
#      tier_python_stdlib, tier_python_geo, tier_pixi)
#   B  the tier matches the command: a python* argv[0] is a Python tier, pixi
#      is tier_pixi; a native command may carry a Python/Pixi tier only by
#      requiring a fixture whose setup test has that tier (propagation), or
#      by passing a Python interpreter as an argument (e.g. -DPYTHON=... to a
#      cmake -P driver that creates a venv from it)
#   C  no tier_native test passes a .py file anywhere in its command
#   D  hermetic implies tier_native, no SKIP_RETURN_CODE, not needs_network
#   E  every test has LSDTT_NONINTERACTIVE=1 in ENVIRONMENT, and a TIMEOUT
#      unless it is benchmark_* / perf* (their options own their budgets)
#   F  a test that runs lsdtt-lem has lem_component
#      (unit_synthetic_cli_banners excepted: a cross-CLI banner check)
#   G  needs_gdal_tools where the test runs GDAL's command-line programs:
#      run_cli_smoke.py with an assertion that shells out to them, or any
#      other script naming gdalinfo / ogrinfo / gdal_translate
#   H  a Python test whose script (or a local module it imports) imports a
#      third-party package is tier_python_geo, and a tier_python_geo test does
#      import one
#   I  no script in tests/scripts starts a process except through
#      lsdtt_test_env.py (whose wrappers keep children off the terminal)
#
# In SELFTEST mode the audit runs over a canned test file and scripts
# directory in which each check is violated on purpose, and passes only if it
# reports exactly the planted violations -- so the audit cannot pass vacuously.

cmake_minimum_required(VERSION 3.15)

if(DEFINED SELFTEST_DIR)
    # Not named CTestTestfile.cmake: .gitignore ignores that name everywhere,
    # and under it this fixture was never committed (found 2026-09-18 when the
    # self-test failed on a fresh clone).
    set(TESTFILE "${SELFTEST_DIR}/canned_CTestTestfile.cmake.in")
    set(SCRIPTS_DIR "${SELFTEST_DIR}/scripts")
    foreach(_f "${TESTFILE}" "${SELFTEST_DIR}/expected.txt" "${SCRIPTS_DIR}")
        if(NOT EXISTS "${_f}")
            message(FATAL_ERROR "audit self-test fixture missing: ${_f}\n"
                "If this is a fresh clone, the file was probably never committed: "
                "check `git check-ignore -v <file>` against .gitignore.")
        endif()
    endforeach()
endif()

foreach(_req TESTFILE SCRIPTS_DIR)
    if(NOT DEFINED ${_req} OR NOT EXISTS "${${_req}}")
        message(FATAL_ERROR "audit_test_inventory: ${_req} not found: '${${_req}}'")
    endif()
endforeach()

# -----------------------------------------------------------------------------
# Record what the generated test file registers.
# -----------------------------------------------------------------------------
set(_audit_tests "")

function(add_test _name)
    set_property(GLOBAL APPEND PROPERTY AUDIT_TESTS "${_name}")
    math(EXPR _last "${ARGC} - 1")
    set(_argv "")
    if(_last GREATER 0)
        foreach(_i RANGE 1 ${_last})
            list(APPEND _argv "${ARGV${_i}}")
        endforeach()
    endif()
    set_property(GLOBAL PROPERTY "AUDIT_ARGV_${_name}" "${_argv}")
endfunction()

function(set_tests_properties _name)
    # set_tests_properties(<name> PROPERTIES k1 "v1" k2 "v2" ...). Values are
    # read through ARGV<n> so a quoted list ("a;b") stays one value.
    math(EXPR _i "2")
    while(_i LESS ARGC)
        math(EXPR _j "${_i} + 1")
        set_property(GLOBAL PROPERTY "AUDIT_PROP_${_name}_${ARGV${_i}}" "${ARGV${_j}}")
        math(EXPR _i "${_i} + 2")
    endwhile()
endfunction()

function(subdirs)
endfunction()

include("${TESTFILE}")
get_property(_tests GLOBAL PROPERTY AUDIT_TESTS)
list(LENGTH _tests _n_tests)

macro(_prop _test _key _out)
    get_property(${_out} GLOBAL PROPERTY "AUDIT_PROP_${_test}_${_key}")
endmacro()

set(_violations "")
# A ';' inside a message would split it into two list entries, so it becomes ','.
macro(_violate _check _msg)
    string(REPLACE ";" "," _vmsg "${_msg}")
    list(APPEND _violations "${_check}: ${_vmsg}")
endmacro()

set(_tiers tier_native tier_python_stdlib tier_python_geo tier_pixi)
set(_third_party osgeo numpy rasterio geopandas shapely pytest pyproj scipy pandas
                 matplotlib fiona richdem whitebox xarray netCDF4 h5py riversand)
set(_harness_gdal_args
    --expect-valid-percent-100 --expect-valid-percent --expect-raster-size
    --expect-raster-min --expect-raster-stat --expect-raster-checksum
    --expect-raster-values --expect-vector-feature-min)

# Third-party imports of a script, following local modules in SCRIPTS_DIR.
function(_script_third_party _path _out)
    set(_found "")
    set(_queue "${_path}")
    set(_seen "")
    while(_queue)
        list(POP_FRONT _queue _file)
        if("${_file}" IN_LIST _seen OR NOT EXISTS "${_file}")
            continue()
        endif()
        list(APPEND _seen "${_file}")
        # An import directly under `try:` is an optional dependency the script
        # guards (run_cosmo_riversand_comparison.py and riversand), not a
        # requirement, so it is skipped.
        file(STRINGS "${_file}" _lines REGEX "^[ \t]*((import|from)[ \t]+[A-Za-z_]|try:)")
        set(_after_try FALSE)
        foreach(_line IN LISTS _lines)
            if(_line MATCHES "^[ \t]*try:[ \t]*$")
                set(_after_try TRUE)
                continue()
            endif()
            if(_after_try)
                set(_after_try FALSE)
                continue()
            endif()
            if(_line MATCHES "^[ \t]*(import|from)[ \t]+([A-Za-z_][A-Za-z0-9_]*)")
                set(_mod "${CMAKE_MATCH_2}")
                if(EXISTS "${SCRIPTS_DIR}/${_mod}.py")
                    list(APPEND _queue "${SCRIPTS_DIR}/${_mod}.py")
                elseif("${_mod}" IN_LIST _third_party)
                    list(APPEND _found "${_mod}")
                endif()
            endif()
        endforeach()
    endwhile()
    list(REMOVE_DUPLICATES _found)
    set(${_out} "${_found}" PARENT_SCOPE)
endfunction()

# Fixture setup tests by fixture name, and each test's tier, for check B.
foreach(_t IN LISTS _tests)
    _prop(${_t} LABELS _labels)
    set(_tier "")
    foreach(_l IN LISTS _labels)
        if("${_l}" IN_LIST _tiers)
            list(APPEND _tier "${_l}")
        endif()
    endforeach()
    set(_tier_of_${_t} "${_tier}")
    _prop(${_t} FIXTURES_SETUP _setups)
    foreach(_fx IN LISTS _setups)
        list(APPEND _setup_of_${_fx} "${_t}")
    endforeach()
endforeach()

# -----------------------------------------------------------------------------
# Per-test checks A-H.
# -----------------------------------------------------------------------------
foreach(_t IN LISTS _tests)
    get_property(_argv GLOBAL PROPERTY "AUDIT_ARGV_${_t}")
    list(GET _argv 0 _exe)
    get_filename_component(_exe_name "${_exe}" NAME)
    _prop(${_t} LABELS _labels)
    _prop(${_t} ENVIRONMENT _env)
    _prop(${_t} TIMEOUT _timeout)
    _prop(${_t} SKIP_RETURN_CODE _skip)
    _prop(${_t} FIXTURES_REQUIRED _required)
    set(_tier "${_tier_of_${_t}}")

    # A
    list(LENGTH _tier _n_tier)
    if(NOT _n_tier EQUAL 1)
        _violate(A "${_t} has ${_n_tier} tier labels (${_tier})")
        continue()
    endif()

    # B
    if(_exe_name MATCHES "^python[0-9.]*$")
        if(NOT _tier MATCHES "^tier_python_(stdlib|geo)$")
            _violate(B "${_t} runs ${_exe_name} but is ${_tier}")
        endif()
    elseif(_exe_name STREQUAL "pixi")
        if(NOT _tier STREQUAL "tier_pixi")
            _violate(B "${_t} runs pixi but is ${_tier}")
        endif()
    elseif(NOT _tier STREQUAL "tier_native")
        set(_inherited FALSE)
        if(_tier MATCHES "^tier_python_")
            foreach(_a IN LISTS _argv)
                string(REGEX REPLACE "^-D[A-Za-z0-9_]+=" "" _av "${_a}")
                get_filename_component(_avn "${_av}" NAME)
                if(_avn MATCHES "^python[0-9.]*(\\.exe)?$")
                    set(_inherited TRUE)
                endif()
            endforeach()
        endif()
        foreach(_fx IN LISTS _required)
            foreach(_s IN LISTS _setup_of_${_fx})
                if("${_tier_of_${_s}}" STREQUAL "${_tier}"
                   OR (_tier STREQUAL "tier_python_stdlib" AND "${_tier_of_${_s}}" STREQUAL "tier_python_geo"))
                    set(_inherited TRUE)
                endif()
            endforeach()
        endforeach()
        if(NOT _inherited)
            _violate(B "${_t} runs native ${_exe_name} but is ${_tier}, and no required fixture explains it")
        endif()
    endif()

    # C
    if(_tier STREQUAL "tier_native")
        foreach(_a IN LISTS _argv)
            if(_a MATCHES "\\.py$")
                _violate(C "${_t} is tier_native but its command includes ${_a}")
            endif()
        endforeach()
    endif()

    # D
    if("hermetic" IN_LIST _labels)
        if(NOT _tier STREQUAL "tier_native")
            _violate(D "${_t} is hermetic but ${_tier}")
        endif()
        if(NOT "${_skip}" STREQUAL "")
            _violate(D "${_t} is hermetic but has SKIP_RETURN_CODE ${_skip}")
        endif()
        if("needs_network" IN_LIST _labels)
            _violate(D "${_t} is hermetic but needs_network")
        endif()
    endif()

    # E
    if(NOT "LSDTT_NONINTERACTIVE=1" IN_LIST _env)
        _violate(E "${_t} lacks LSDTT_NONINTERACTIVE=1")
    endif()
    if("${_timeout}" STREQUAL "" AND NOT _t MATCHES "^(benchmark_|perf)")
        _violate(E "${_t} has no TIMEOUT")
    endif()

    # F
    if(NOT _t STREQUAL "unit_synthetic_cli_banners" AND NOT "lem_component" IN_LIST _labels)
        foreach(_a IN LISTS _argv)
            get_filename_component(_an "${_a}" NAME)
            if(_an STREQUAL "lsdtt-lem")
                _violate(F "${_t} runs lsdtt-lem but lacks lem_component")
                break()
            endif()
        endforeach()
    endif()

    # G and H need the script, when there is one.
    set(_script "")
    foreach(_a IN LISTS _argv)
        if(_a MATCHES "\\.py$")
            set(_script "${_a}")
            break()
        endif()
    endforeach()
    if(_script)
        get_filename_component(_script_name "${_script}" NAME)
        set(_uses_tools FALSE)
        if(_script_name STREQUAL "run_cli_smoke.py")
            foreach(_a IN LISTS _argv)
                if("${_a}" IN_LIST _harness_gdal_args)
                    set(_uses_tools TRUE)
                endif()
            endforeach()
        elseif(EXISTS "${_script}")
            file(STRINGS "${_script}" _tool_lines REGEX "[\"'](gdalinfo|ogrinfo|gdal_translate)[\"']")
            if(_tool_lines)
                set(_uses_tools TRUE)
            endif()
        endif()
        if(_uses_tools AND NOT "needs_gdal_tools" IN_LIST _labels)
            _violate(G "${_t} runs GDAL's command-line programs but lacks needs_gdal_tools")
        endif()

        if(_tier MATCHES "^tier_python_")
            _script_third_party("${_script}" _imports)
            if(_imports AND _tier STREQUAL "tier_python_stdlib")
                _violate(H "${_t} is tier_python_stdlib but ${_script_name} imports ${_imports}")
            elseif(NOT _imports AND _tier STREQUAL "tier_python_geo")
                _violate(H "${_t} is tier_python_geo but ${_script_name} imports no third-party package")
            endif()
        endif()
    endif()
endforeach()

# -----------------------------------------------------------------------------
# Check I: every process launch in tests/scripts goes through lsdtt_test_env.
# -----------------------------------------------------------------------------
file(GLOB _scripts "${SCRIPTS_DIR}/*.py")
set(_n_scripts 0)
foreach(_s IN LISTS _scripts)
    get_filename_component(_sn "${_s}" NAME)
    if(_sn STREQUAL "lsdtt_test_env.py")
        continue()
    endif()
    math(EXPR _n_scripts "${_n_scripts} + 1")
    file(STRINGS "${_s}" _direct
         REGEX "subprocess\\.(run|Popen|call|check_call|check_output)\\(|os\\.(system|popen)\\(")
    foreach(_line IN LISTS _direct)
        string(STRIP "${_line}" _line)
        _violate(I "${_sn} launches a process directly (${_line}); use lsdtt_test_env")
    endforeach()
endforeach()

# -----------------------------------------------------------------------------
# Report.
# -----------------------------------------------------------------------------
list(LENGTH _violations _n_viol)
message("test inventory audit: ${_n_tests} tests, ${_n_scripts} scripts, "
        "${_n_viol} violation(s)")

if(DEFINED SELFTEST_DIR)
    # Each planted violation is named in expected.txt as "<check>: <test or
    # script>"; the audit must report exactly that set, no more and no less.
    file(STRINGS "${SELFTEST_DIR}/expected.txt" _expected REGEX "^[A-I]: ")
    set(_got "")
    foreach(_v IN LISTS _violations)
        if(_v MATCHES "^([A-I]): ([^ ]+)")
            list(APPEND _got "${CMAKE_MATCH_1}: ${CMAKE_MATCH_2}")
        endif()
    endforeach()
    list(SORT _got)
    list(SORT _expected)
    foreach(_v IN LISTS _violations)
        message("  ${_v}")
    endforeach()
    if(NOT "${_got}" STREQUAL "${_expected}")
        string(REPLACE ";" "\n    " _g "${_got}")
        string(REPLACE ";" "\n    " _e "${_expected}")
        message(FATAL_ERROR "audit self-test: reported violations differ from the planted ones.\n"
                            "  planted:\n    ${_e}\n  reported:\n    ${_g}")
    endif()
    message("audit self-test: every planted violation reported, nothing else")
    return()
endif()

if(_n_viol GREATER 0)
    string(REPLACE ";" "\n  " _report "${_violations}")
    message(FATAL_ERROR
        "\n  ${_report}\n\n"
        "The DEPENDENCY TIERS block and the execution-contract block at the foot "
        "of tests/CMakeLists.txt set these labels and properties; fix the "
        "registration (or the lists there), not this audit.")
endif()
