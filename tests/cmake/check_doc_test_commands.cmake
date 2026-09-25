# check_doc_test_commands.cmake
#
# Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
#
# Copyright (C) 2026 Simon M. Mudd
# SPDX-License-Identifier: AGPL-3.0-or-later
#
# Checks that the test commands the documentation tells people to type still
# mean something:
#
#   * every -DLSDTT_<NAME> option it names is declared by an option() in
#     CMakeLists.txt or tests/CMakeLists.txt;
#   * every ANCHORED label it selects (-L '^x$', -LE '^(a|b)$') is a label the
#     test CMake can assign -- a LABELS value, a domain rule, or one of the
#     tier / flag / component / cost labels;
#   * the `test_report` target and `run_tests.cmake` it names exist.
#
# Unanchored `-L x` in historical records ("ctest -L network -LE heavy
# 92/92") is left alone: those record what was run at the time.
#
#   cmake -DSOURCE_DIR=<repo> -P check_doc_test_commands.cmake

cmake_minimum_required(VERSION 3.15)

if(NOT DEFINED SOURCE_DIR OR NOT EXISTS "${SOURCE_DIR}/tests/CMakeLists.txt")
    message(FATAL_ERROR "check_doc_test_commands: SOURCE_DIR not set or not the repository")
endif()

set(_docs
    README.md
    tests/README.md
    docs/developer_notes/building_from_source.md
    docs/developer_notes/testing_and_benchmarking.md
    apptainer_def/apptainer_setup.md
    pixi.toml)

# Declared options.
set(_options "")
foreach(_cm CMakeLists.txt tests/CMakeLists.txt)
    file(STRINGS "${SOURCE_DIR}/${_cm}" _lines REGEX "option\\((LSDTT_[A-Z0-9_]+)")
    foreach(_l IN LISTS _lines)
        if(_l MATCHES "option\\((LSDTT_[A-Z0-9_]+)")
            list(APPEND _options "${CMAKE_MATCH_1}")
        endif()
    endforeach()
endforeach()

# Labels the test CMake can assign.
set(_labels tier_native tier_python_stdlib tier_python_geo tier_pixi
            hermetic needs_network needs_gdal_tools lem_component heavy)
file(READ "${SOURCE_DIR}/tests/CMakeLists.txt" _tc)
# MATCHALL returns a list, so a ';' inside a match would split it: use ','.
string(REPLACE ";" "," _tc_flat "${_tc}")
# Literal LABELS "a;b", and the helpers' label variables (test_labels, lbl).
string(REGEX MATCHALL "(LABELS[ \t\n]+|(test_labels|lbl) +)\"[^\"]*\"" _label_props "${_tc_flat}")
foreach(_lp IN LISTS _label_props)
    string(REGEX REPLACE "^[^\"]*\"([^\"]*)\"$" "\\1" _v "${_lp}")
    string(REPLACE "," ";" _v "${_v}")
    foreach(_one IN LISTS _v)
        if(_one MATCHES "^[a-z_][a-z0-9_]*$")
            list(APPEND _labels "${_one}")
        endif()
    endforeach()
endforeach()
# Domain rules: lines of the form   "<pattern>"   "<domain>"
file(STRINGS "${SOURCE_DIR}/tests/CMakeLists.txt" _rules REGEX "^[ \t]+\"[^\"]+\"[ \t]+\"[a-z_]+\"[ \t]*$")
foreach(_r IN LISTS _rules)
    if(_r MATCHES "\"([a-z_]+)\"[ \t]*$")
        list(APPEND _labels "${CMAKE_MATCH_1}")
    endif()
endforeach()
list(REMOVE_DUPLICATES _labels)

set(_problems "")
set(_n_options 0)
set(_n_labels 0)
foreach(_doc IN LISTS _docs)
    set(_path "${SOURCE_DIR}/${_doc}")
    if(NOT EXISTS "${_path}")
        list(APPEND _problems "${_doc}: file missing")
        continue()
    endif()
    file(READ "${_path}" _text)

    string(REGEX MATCHALL "-D(LSDTT_[A-Z0-9_]+)" _opts "${_text}")
    foreach(_o IN LISTS _opts)
        string(SUBSTRING "${_o}" 2 -1 _name)
        math(EXPR _n_options "${_n_options} + 1")
        if(NOT "${_name}" IN_LIST _options)
            list(APPEND _problems "${_doc}: -D${_name} is not a declared option")
        endif()
    endforeach()

    # Anchored selections: -L '^x$'  -LE '^(a|b)$'  (quotes ' or ", $ maybe escaped)
    string(REGEX MATCHALL "-LE?[ \t]+['\"]\\^[(]?[a-z_|]+[)]?\\\\?\\$['\"]" _sels "${_text}")
    foreach(_sel IN LISTS _sels)
        string(REGEX REPLACE "^-LE?[ \t]+['\"]\\^[(]?([a-z_|]+)[)]?\\\\?\\$['\"]$" "\\1" _names "${_sel}")
        string(REPLACE "|" ";" _names "${_names}")
        foreach(_lab IN LISTS _names)
            math(EXPR _n_labels "${_n_labels} + 1")
            if(NOT "${_lab}" IN_LIST _labels)
                list(APPEND _problems "${_doc}: label '${_lab}' in ${_sel} is never assigned")
            endif()
        endforeach()
    endforeach()

    if(_text MATCHES "--target test_report" AND NOT _tc MATCHES "add_custom_target\\(test_report")
        list(APPEND _problems "${_doc}: names the test_report target, which does not exist")
    endif()
    if(_text MATCHES "run_tests\\.cmake" AND NOT EXISTS "${SOURCE_DIR}/tests/cmake/run_tests.cmake.in")
        list(APPEND _problems "${_doc}: names run_tests.cmake, whose template is missing")
    endif()
endforeach()

list(LENGTH _docs _n_docs)
message("documentation test commands: ${_n_docs} files, ${_n_options} option uses, "
        "${_n_labels} anchored label uses checked")
if(_problems)
    string(REPLACE ";" "\n  " _problems "${_problems}")
    message(FATAL_ERROR "\n  ${_problems}")
endif()
