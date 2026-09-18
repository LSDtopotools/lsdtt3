# make_test_report.cmake
#
# Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
#
# Copyright (C) 2026 Simon M. Mudd
# SPDX-License-Identifier: AGPL-3.0-or-later
#
# Builds ONE plain-text file a user can send when tests fail: what machine,
# toolchain, dependency versions and build configuration produced the failures,
# and the failing tests' output. Run it through the build target:
#
#   cmake --build <build> --target test_report
#
# which calls
#
#   cmake -DBUILD_DIR=<build> [-DOUT=<file>] -P make_test_report.cmake
#
# It works after a completed run and after a killed one: it reads
# Testing/Temporary/LastTest.log, or LastTest.log.tmp when that is newer (the
# only record CTest leaves of an interrupted run). Test status (passed, failed,
# timeout, skipped) comes from the console log written by run_tests.cmake,
# because LastTest.log cannot tell a timeout or a skip from a failure or a
# pass; without that log, timeouts are inferred from duration against TIMEOUT.
#
# PRIVACY AND SECRETS
#
#   * The home directory becomes "~" and the user name "<user>" everywhere.
#   * There is NO environment dump. Only the variables named in
#     _report_env_vars below are reported, and nothing matching KEY, TOKEN,
#     SECRET or PASSWORD is ever read.
#   * The report ends by asking the user to read it before sending it.
#
# No program is written for this: the run-time GDAL/PROJ facts come from an
# existing CLI's --version (CLIs exist for analysis only).
#
# SELF-TEST: -DSELFTEST_DIR=<tests/cmake/report_selftest> -DSELFTEST_WORK=<dir>
# copies a canned build tree (placeholders @HOME@ and @BUILD@ substituted),
# builds a report from it, and checks the sections, the failure and timeout
# entries, the masking of the home directory, and that the value of a planted
# OT_API_KEY in the environment appears nowhere.

cmake_minimum_required(VERSION 3.15)

# -----------------------------------------------------------------------------
# Self-test set-up.
# -----------------------------------------------------------------------------
if(DEFINED SELFTEST_DIR)
    file(REMOVE_RECURSE "${SELFTEST_WORK}")
    file(COPY "${SELFTEST_DIR}/build/" DESTINATION "${SELFTEST_WORK}")
    file(GLOB_RECURSE _canned LIST_DIRECTORIES FALSE "${SELFTEST_WORK}/*")
    foreach(_f IN LISTS _canned)
        file(READ "${_f}" _c)
        string(REPLACE "@HOME@" "$ENV{HOME}" _c "${_c}")
        string(REPLACE "@BUILD@" "${SELFTEST_WORK}" _c "${_c}")
        file(WRITE "${_f}" "${_c}")
    endforeach()
    set(BUILD_DIR "${SELFTEST_WORK}")
    set(OUT "${SELFTEST_WORK}/report.txt")
endif()

if(NOT DEFINED BUILD_DIR OR NOT EXISTS "${BUILD_DIR}")
    message(FATAL_ERROR "make_test_report: BUILD_DIR not found: '${BUILD_DIR}'")
endif()
set(_info "${BUILD_DIR}/lsdtt3_build_info.cmake")
if(NOT EXISTS "${_info}")
    message(FATAL_ERROR "make_test_report: ${_info} is missing; re-run cmake on the build tree.")
endif()
include("${_info}")

# -----------------------------------------------------------------------------
# Helpers. Text read from files is "protected" before it is split into a CMake
# list, because ';', '[', ']' and '\' change how lists split.
# -----------------------------------------------------------------------------
macro(_protect _v)
    string(REPLACE "\\" "<<BS>>" ${_v} "${${_v}}")
    string(REPLACE ";" "<<SC>>" ${_v} "${${_v}}")
    string(REPLACE "[" "<<LB>>" ${_v} "${${_v}}")
    string(REPLACE "]" "<<RB>>" ${_v} "${${_v}}")
endmacro()
macro(_unprotect _v)
    string(REPLACE "<<RB>>" "]" ${_v} "${${_v}}")
    string(REPLACE "<<LB>>" "[" ${_v} "${${_v}}")
    string(REPLACE "<<SC>>" ";" ${_v} "${${_v}}")
    string(REPLACE "<<BS>>" "\\" ${_v} "${${_v}}")
endmacro()

# Lines of a file as a protected list.
function(_read_lines _path _out)
    file(READ "${_path}" _t)
    _protect(_t)
    string(REPLACE "\r" "" _t "${_t}")
    string(REPLACE "\n" ";" _t "${_t}")
    set(${_out} "${_t}" PARENT_SCOPE)
endfunction()

# The last N lines of a file, as text.
function(_tail _path _n _out)
    if(NOT EXISTS "${_path}")
        set(${_out} "    (no ${_path})\n" PARENT_SCOPE)
        return()
    endif()
    # The modification time says whether the file belongs to this run: a
    # harness that fails before launching the CLI leaves an older run's logs.
    file(TIMESTAMP "${_path}" _mtime "%Y-%m-%d %H:%M:%S")
    _read_lines("${_path}" _l)
    list(LENGTH _l _len)
    # A trailing newline leaves an empty last element; drop it.
    if(_len GREATER 0)
        math(EXPR _lastidx "${_len} - 1")
        list(GET _l ${_lastidx} _lastline)
        if(_lastline STREQUAL "")
            list(REMOVE_AT _l ${_lastidx})
            math(EXPR _len "${_len} - 1")
        endif()
    endif()
    set(_start 0)
    if(_len GREATER _n)
        math(EXPR _start "${_len} - ${_n}")
    endif()
    set(_t "    [file modified ${_mtime}]\n")
    if(_start GREATER 0)
        string(APPEND _t "    ... (${_start} earlier lines omitted)\n")
    endif()
    if(_len GREATER 0)
        math(EXPR _last "${_len} - 1")
        foreach(_i RANGE ${_start} ${_last})
            list(GET _l ${_i} _line)
            string(APPEND _t "    ${_line}\n")
        endforeach()
    endif()
    _unprotect(_t)
    set(${_out} "${_t}" PARENT_SCOPE)
endfunction()

set(_R "")
macro(_line _text)
    string(APPEND _R "${_text}\n")
endmacro()
macro(_section _title)
    string(APPEND _R "\n== ${_title} ==\n")
endmacro()

string(TIMESTAMP _now "%Y-%m-%d %H:%M:%S")
_line("lsdtt3 test report")
_line("generated: ${_now}")

# -----------------------------------------------------------------------------
# 1. Build and source
# -----------------------------------------------------------------------------
_section("Build")
_line("lsdtt3 version    : ${LSDTT_REPORT_VERSION}")
_line("build type        : ${LSDTT_REPORT_BUILD_TYPE}")
_line("compiler          : ${LSDTT_REPORT_CXX_ID} ${LSDTT_REPORT_CXX_VERSION}, C++${LSDTT_REPORT_CXX_STANDARD}")
_line("generator         : ${LSDTT_REPORT_GENERATOR} (CMake ${LSDTT_REPORT_CMAKE_VERSION})")
_line("options           : ${LSDTT_REPORT_OPTIONS}")
_line("found             : OpenMP=${LSDTT_REPORT_HAS_OPENMP} FFTW3=${LSDTT_REPORT_HAS_FFTW} Eigen3=${LSDTT_REPORT_HAS_EIGEN} NetCDF=${LSDTT_REPORT_HAS_NETCDF}")
_line("GDAL (compiled)   : ${LSDTT_REPORT_GDAL_VERSION} (${LSDTT_REPORT_GDAL_LIBRARY})")

# -----------------------------------------------------------------------------
# 2. Machine
# -----------------------------------------------------------------------------
_section("Machine")
cmake_host_system_information(RESULT _os QUERY OS_NAME OS_RELEASE OS_PLATFORM)
string(REPLACE ";" " " _os "${_os}")
_line("OS                : ${_os}")
if(EXISTS "/etc/os-release")
    file(STRINGS "/etc/os-release" _pretty REGEX "^PRETTY_NAME=")
    string(REGEX REPLACE "^PRETTY_NAME=\"?([^\"]*)\"?$" "\\1" _pretty "${_pretty}")
    _line("distribution      : ${_pretty}")
endif()
set(_wsl "no")
if(EXISTS "/proc/version")
    file(READ "/proc/version" _pv)
    string(TOLOWER "${_pv}" _pv)
    if(_pv MATCHES "microsoft")
        set(_wsl "yes")
    endif()
endif()
_line("WSL               : ${_wsl}")
cmake_host_system_information(RESULT _cores QUERY NUMBER_OF_LOGICAL_CORES)
cmake_host_system_information(RESULT _pcores QUERY NUMBER_OF_PHYSICAL_CORES)
cmake_host_system_information(RESULT _mem QUERY TOTAL_PHYSICAL_MEMORY)
cmake_host_system_information(RESULT _cpu QUERY PROCESSOR_DESCRIPTION)
_line("cores             : ${_cores} logical, ${_pcores} physical")
_line("memory            : ${_mem} MiB")
_line("CPU               : ${_cpu}")

# -----------------------------------------------------------------------------
# 3. Geospatial run time
# -----------------------------------------------------------------------------
_section("Geospatial runtime")
_line("GDAL compiled     : ${LSDTT_REPORT_GDAL_VERSION}")
if(EXISTS "${LSDTT_REPORT_CLI}")
    execute_process(COMMAND "${LSDTT_REPORT_CLI}" --version
                    OUTPUT_VARIABLE _v ERROR_VARIABLE _v TIMEOUT 30)
    _protect(_v)
    string(REPLACE "\n" ";" _vl "${_v}")
    foreach(_l IN LISTS _vl)
        if(_l MATCHES "version|GDAL|PROJ|proj data|GDAL data" AND NOT _l MATCHES "^\\|\\|")
            _unprotect(_l)
            string(STRIP "${_l}" _l)
            _line("CLI --version     : ${_l}")
        endif()
    endforeach()
    _line("  (as run from the environment this report was made in)")
else()
    _line("CLI --version     : ${LSDTT_REPORT_CLI} is not built")
endif()
foreach(_prog gdalinfo ogrinfo gdal_translate)
    unset(_found)
    find_program(_found ${_prog})
    string(SUBSTRING "${_prog}                  " 0 18 _label)
    if(_found)
        _line("${_label}: ${_found}")
    else()
        _line("${_label}: not found on PATH")
    endif()
    unset(_found CACHE)
endforeach()

# -----------------------------------------------------------------------------
# 4. Test environment
# -----------------------------------------------------------------------------
_section("Test environment")
execute_process(COMMAND "${LSDTT_REPORT_CTEST}" --version OUTPUT_VARIABLE _cv ERROR_QUIET)
string(REGEX REPLACE "\n.*" "" _cv "${_cv}")
_line("CTest             : ${_cv}")
if(LSDTT_REPORT_PYTHON AND EXISTS "${LSDTT_REPORT_PYTHON}")
    execute_process(COMMAND "${LSDTT_REPORT_PYTHON}" --version
                    OUTPUT_VARIABLE _pyv ERROR_VARIABLE _pyv)
    string(STRIP "${_pyv}" _pyv)
    _line("Python            : ${LSDTT_REPORT_PYTHON} (${_pyv})")
    if(LSDTT_REPORT_GEO_TESTS)
        foreach(_mod osgeo.gdal numpy rasterio geopandas shapely)
            execute_process(COMMAND "${LSDTT_REPORT_PYTHON}" -c "import ${_mod}"
                            RESULT_VARIABLE _rc OUTPUT_QUIET ERROR_QUIET)
            if(_rc EQUAL 0)
                _line("  import ${_mod} : ok")
            else()
                _line("  import ${_mod} : FAILS")
            endif()
        endforeach()
    endif()
elseif(LSDTT_REPORT_PYTHON)
    _line("Python            : ${LSDTT_REPORT_PYTHON} (configured, but no longer exists)")
else()
    _line("Python            : not used (LSDTT_ENABLE_PYTHON_TESTS=OFF)")
endif()
# A FIXED list. Never add a name matching KEY, TOKEN, SECRET or PASSWORD.
set(_report_env_vars CONDA_PREFIX PIXI_ENVIRONMENT_NAME PIXI_PROJECT_ROOT
    GDAL_DATA PROJ_DATA PROJ_LIB GDAL_DRIVER_PATH LD_LIBRARY_PATH
    OMP_NUM_THREADS OMP_WAIT_POLICY)
foreach(_var IN LISTS _report_env_vars)
    if(DEFINED ENV{${_var}} AND NOT "$ENV{${_var}}" STREQUAL "")
        _line("${_var} : $ENV{${_var}}")
    else()
        _line("${_var} : not set")
    endif()
endforeach()

# -----------------------------------------------------------------------------
# Read the registered tests (labels, timeouts, skip codes) from the generated
# CTestTestfile, with add_test / set_tests_properties recording, as the
# inventory audit does.
# -----------------------------------------------------------------------------
function(add_test _name)
    set_property(GLOBAL APPEND PROPERTY REPORT_TESTS "${_name}")
    # The command as CTest's log would show it, for tests that never reached
    # the log (running when an interrupted run stopped).
    set(_cmd "")
    math(EXPR _last "${ARGC} - 1")
    if(_last GREATER 0)
        foreach(_i RANGE 1 ${_last})
            string(APPEND _cmd "\"${ARGV${_i}}\" ")
        endforeach()
    endif()
    string(STRIP "${_cmd}" _cmd)
    _protect(_cmd)
    set_property(GLOBAL PROPERTY "REPORT_${_name}_COMMAND" "${_cmd}")
endfunction()
function(set_tests_properties _name)
    math(EXPR _i "2")
    while(_i LESS ARGC)
        math(EXPR _j "${_i} + 1")
        set_property(GLOBAL PROPERTY "REPORT_${_name}_${ARGV${_i}}" "${ARGV${_j}}")
        math(EXPR _i "${_i} + 2")
    endwhile()
endfunction()
function(subdirs)
endfunction()
if(EXISTS "${BUILD_DIR}/tests/CTestTestfile.cmake")
    include("${BUILD_DIR}/tests/CTestTestfile.cmake")
endif()

# -----------------------------------------------------------------------------
# 5. Run summary
# -----------------------------------------------------------------------------
_section("Run summary")
set(_tmp_dir "${BUILD_DIR}/Testing/Temporary")
set(_ctest_log "")
if(EXISTS "${_tmp_dir}/LastTest.log.tmp" AND
   (NOT EXISTS "${_tmp_dir}/LastTest.log" OR
    "${_tmp_dir}/LastTest.log.tmp" IS_NEWER_THAN "${_tmp_dir}/LastTest.log"))
    set(_ctest_log "${_tmp_dir}/LastTest.log.tmp")
elseif(EXISTS "${_tmp_dir}/LastTest.log")
    set(_ctest_log "${_tmp_dir}/LastTest.log")
endif()

# run_tests.cmake's record of the last wrapped run.
set(_console_log "")
if(EXISTS "${BUILD_DIR}/test-logs/last_run.txt")
    file(STRINGS "${BUILD_DIR}/test-logs/last_run.txt" _lr)
    set(_lr_ended FALSE)
    foreach(_l IN LISTS _lr)
        _line("last wrapped run  : ${_l}")
        if(_l MATCHES "^console_log: (.+)$")
            set(_console_log "${CMAKE_MATCH_1}")
        elseif(_l MATCHES "^end: ")
            set(_lr_ended TRUE)
        endif()
    endforeach()
    if(NOT _lr_ended)
        _line("last wrapped run  : no end recorded -- the run was killed or is still going")
    endif()
else()
    _line("last wrapped run  : none (run_tests.cmake not used; ctest arguments not recorded)")
endif()

# The console log belongs to the same run only if it was written at about the
# same time as the CTest log; a later plain `ctest` makes it stale.
if(_console_log AND EXISTS "${_console_log}" AND _ctest_log)
    file(TIMESTAMP "${_console_log}" _tc "%s" UTC)
    file(TIMESTAMP "${_ctest_log}" _tl "%s" UTC)
    math(EXPR _gap "${_tl} - ${_tc}")
    if(_gap GREATER 60)
        _line("console log       : ${_console_log} is older than the CTest log; not used")
        set(_console_log "")
    endif()
elseif(_console_log AND NOT EXISTS "${_console_log}")
    set(_console_log "")
endif()

if(NOT _ctest_log)
    _line("CTest log         : none found in ${_tmp_dir} -- run the tests first")
    set(_blocks "")
else()
    _line("CTest log         : ${_ctest_log}")
endif()

# Parse the CTest log into per-test results.
set(_names "")
if(_ctest_log)
    _read_lines("${_ctest_log}" _lines)
    set(_completed FALSE)
    set(_cur "")
    set(_in_out FALSE)
    set(_out "")
    foreach(_l IN LISTS _lines)
        if(_l MATCHES "^[0-9]+/[0-9]+ Testing: (.+)$")
            set(_cur "${CMAKE_MATCH_1}")
            list(APPEND _names "${_cur}")
            set(_res_${_cur} "no result (hung or killed)")
            set(_out "")
        elseif(_cur AND _l MATCHES "^Command: (.*)$")
            set(_cmd_${_cur} "${CMAKE_MATCH_1}")
        elseif(_cur AND _l STREQUAL "Output:")
            set(_in_out TRUE)
        elseif(_in_out AND _l STREQUAL "<end of output>")
            set(_in_out FALSE)
            set(_out_${_cur} "${_out}")
        elseif(_in_out)
            if(NOT _l MATCHES "^-+$" OR NOT _out STREQUAL "")
                list(APPEND _out "${_l}")
            endif()
        elseif(_cur AND _l MATCHES "^Test time = +([0-9.]+) sec")
            set(_time_${_cur} "${CMAKE_MATCH_1}")
        elseif(_cur AND _l MATCHES "^Test (Passed|Failed)\\.$")
            set(_res_${_cur} "${CMAKE_MATCH_1}")
        elseif(_l MATCHES "^End testing:")
            set(_completed TRUE)
        endif()
    endforeach()
    if(_completed)
        _line("run completed     : yes")
    else()
        _line("run completed     : NO -- interrupted; the tests listed as 'no result' were running when it stopped")
    endif()
endif()

# Status per test: the console log's word where there is one. A test the
# console log shows starting, with no result line, was running when the run
# stopped; CTest writes a test's block to its own log only when the test ends,
# so such a test is in the console log alone.
if(_console_log)
    _read_lines("${_console_log}" _clines)
    foreach(_l IN LISTS _clines)
        if(_l MATCHES "^ *Start +[0-9]+: ([^ ]+)$")
            set(_n "${CMAKE_MATCH_1}")
            if(NOT "${_n}" IN_LIST _names)
                list(APPEND _names "${_n}")
                set(_res_${_n} "no result (hung or killed)")
            endif()
        endif()
        if(_l MATCHES "Test +#[0-9]+: ([^ ]+) [ .]*(.*[^ ]) +[0-9.]+ sec$")
            set(_n "${CMAKE_MATCH_1}")
            set(_s "${CMAKE_MATCH_2}")
            if(_s MATCHES "Timeout")
                set(_status_${_n} "Timeout")
            elseif(_s MATCHES "Skipped")
                set(_status_${_n} "Skipped")
            elseif(_s MATCHES "Not Run")
                set(_status_${_n} "Not Run")
            elseif(_s MATCHES "Exception|aborted")
                set(_status_${_n} "Exception")
            elseif(_s MATCHES "Failed")
                set(_status_${_n} "Failed")
            elseif(_s MATCHES "Passed")
                set(_status_${_n} "Passed")
            endif()
        endif()
    endforeach()
endif()

set(_counts_Passed 0)
set(_counts_Failed 0)
set(_counts_Timeout 0)
set(_counts_Skipped 0)
set(_counts_Exception 0)
set(_counts_NotRun 0)
set(_counts_NoResult 0)
set(_bad "")
foreach(_n IN LISTS _names)
    get_property(_to GLOBAL PROPERTY "REPORT_${_n}_TIMEOUT")
    if(DEFINED _status_${_n})
        set(_st "${_status_${_n}}")
    elseif(_res_${_n} STREQUAL "Failed")
        set(_st "Failed")
        if(_to AND _time_${_n})
            # Without the console log a timeout looks like a failure; infer it.
            math(EXPR _to_int "${_to}")
            string(REGEX REPLACE "\\..*" "" _t_int "${_time_${_n}}")
            if(NOT _t_int LESS _to_int)
                set(_st "Timeout (inferred from duration)")
            endif()
        endif()
    elseif(_res_${_n} STREQUAL "Passed")
        set(_st "Passed")
    else()
        set(_st "No result (was running when the run stopped)")
    endif()
    set(_final_${_n} "${_st}")
    if(_st STREQUAL "Passed")
        math(EXPR _counts_Passed "${_counts_Passed} + 1")
    elseif(_st STREQUAL "Skipped")
        math(EXPR _counts_Skipped "${_counts_Skipped} + 1")
    else()
        list(APPEND _bad "${_n}")
        if(_st MATCHES "^Timeout")
            math(EXPR _counts_Timeout "${_counts_Timeout} + 1")
        elseif(_st STREQUAL "Exception")
            math(EXPR _counts_Exception "${_counts_Exception} + 1")
        elseif(_st STREQUAL "Not Run")
            math(EXPR _counts_NotRun "${_counts_NotRun} + 1")
        elseif(_st MATCHES "^No result")
            math(EXPR _counts_NoResult "${_counts_NoResult} + 1")
        else()
            math(EXPR _counts_Failed "${_counts_Failed} + 1")
        endif()
    endif()
endforeach()
list(LENGTH _names _n_run)
_line("tests in log      : ${_n_run}")
_line("passed            : ${_counts_Passed}")
_line("failed            : ${_counts_Failed}")
_line("timed out         : ${_counts_Timeout}")
_line("skipped           : ${_counts_Skipped}")
_line("exception         : ${_counts_Exception}")
_line("not run           : ${_counts_NotRun}")
_line("no result         : ${_counts_NoResult}")

# -----------------------------------------------------------------------------
# 6. Failures
# -----------------------------------------------------------------------------
_section("Failures")
if(NOT _bad)
    _line("none")
endif()
foreach(_n IN LISTS _bad)
    get_property(_labels GLOBAL PROPERTY "REPORT_${_n}_LABELS")
    get_property(_to GLOBAL PROPERTY "REPORT_${_n}_TIMEOUT")
    string(REPLACE ";" "," _labels "${_labels}")
    string(APPEND _R "\n--- ${_n} : ${_final_${_n}}\n")
    _line("labels  : ${_labels}")
    _line("timeout : ${_to} s    duration: ${_time_${_n}} s")
    set(_cmd "${_cmd_${_n}}")
    if(_cmd STREQUAL "")
        get_property(_cmd GLOBAL PROPERTY "REPORT_${_n}_COMMAND")
    endif()
    _unprotect(_cmd)
    _line("command : ${_cmd}")
    _line("CTest output (last 200 lines):")
    set(_o "${_out_${_n}}")
    list(LENGTH _o _olen)
    set(_ostart 0)
    if(_olen GREATER 200)
        math(EXPR _ostart "${_olen} - 200")
    endif()
    set(_otext "")
    if(_olen GREATER 0)
        math(EXPR _olast "${_olen} - 1")
        foreach(_i RANGE ${_ostart} ${_olast})
            list(GET _o ${_i} _ol)
            string(APPEND _otext "    ${_ol}\n")
        endforeach()
    else()
        set(_otext "    (none)\n")
    endif()
    _unprotect(_otext)
    string(APPEND _R "${_otext}")
    # Python-harness tests keep their own logs in --work-dir.
    if(_cmd MATCHES "\"--work-dir\" \"([^\"]+)\"")
        set(_wd "${CMAKE_MATCH_1}")
        _line("harness command.txt:")
        _tail("${_wd}/command.txt" 20 _t)
        string(APPEND _R "${_t}")
        _line("harness stdout.log (last 200 lines):")
        _tail("${_wd}/stdout.log" 200 _t)
        string(APPEND _R "${_t}")
        _line("harness stderr.log (last 200 lines):")
        _tail("${_wd}/stderr.log" 200 _t)
        string(APPEND _R "${_t}")
    endif()
endforeach()

_line("")
_line("Read this report before sending it: it contains paths and command lines")
_line("from your machine (home directory shown as ~, user name as <user>).")

# -----------------------------------------------------------------------------
# Masking, then write.
# -----------------------------------------------------------------------------
if(DEFINED ENV{HOME} AND NOT "$ENV{HOME}" STREQUAL "" AND NOT "$ENV{HOME}" STREQUAL "/")
    string(REPLACE "$ENV{HOME}" "~" _R "${_R}")
endif()
foreach(_uvar USER USERNAME LOGNAME)
    if(DEFINED ENV{${_uvar}} AND NOT "$ENV{${_uvar}}" STREQUAL "")
        string(REPLACE "/$ENV{${_uvar}}/" "/<user>/" _R "${_R}")
        string(REPLACE "\\$ENV{${_uvar}}\\" "\\<user>\\" _R "${_R}")
    endif()
endforeach()

if(NOT DEFINED OUT)
    string(TIMESTAMP _stamp "%Y%m%d_%H%M%S")
    set(OUT "${BUILD_DIR}/test-logs/lsdtt3_test_report_${_stamp}.txt")
endif()
file(WRITE "${OUT}" "${_R}")
message("lsdtt3 test report written: ${OUT}")
message("Read it before sending it; it is plain text.")

# -----------------------------------------------------------------------------
# Self-test assertions.
# -----------------------------------------------------------------------------
if(DEFINED SELFTEST_DIR)
    file(READ "${OUT}" _rep)
    set(_fail "")
    foreach(_want
            "== Build ==" "== Machine ==" "== Geospatial runtime ==" "== Test environment =="
            "== Run summary ==" "== Failures =="
            "--- sel_fail : Failed" "--- sel_timeout : Timeout"
            "--- sel_running : No result (was running when the run stopped)"
            "no result         : 1"
            "planted failure output" "planted harness stderr" "~/data/dem.tif"
            "args: -j4 -LE ^heavy$" "passed            : 1" "skipped           : 1"
            "Read this report before sending it")
        string(FIND "${_rep}" "${_want}" _pos)
        if(_pos LESS 0)
            list(APPEND _fail "missing: ${_want}")
        endif()
    endforeach()
    foreach(_forbid "--- sel_pass" "--- sel_skip" "$ENV{HOME}/" "planted-secret-8d3f")
        string(FIND "${_rep}" "${_forbid}" _pos)
        if(NOT _pos LESS 0)
            list(APPEND _fail "present but must not be: ${_forbid}")
        endif()
    endforeach()
    if(NOT "$ENV{OT_API_KEY}" STREQUAL "planted-secret-8d3f")
        list(APPEND _fail "the self-test must run with OT_API_KEY=planted-secret-8d3f in its environment")
    endif()
    if(_fail)
        string(REPLACE ";" "\n  " _fail "${_fail}")
        message(FATAL_ERROR "test report self-test failed:\n  ${_fail}\n(report: ${OUT})")
    endif()
    message("test report self-test: all sections, the failure and timeout entries, "
            "home masking and secret absence checked")
endif()
