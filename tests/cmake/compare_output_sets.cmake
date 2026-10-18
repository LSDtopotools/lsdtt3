# compare_output_sets.cmake
#
# Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
#
# Copyright (C) 2026 Simon M. Mudd
# SPDX-License-Identifier: AGPL-3.0-or-later
# Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
# See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
#
# Byte-compare the same outputs of two CLI runs:
#   cmake -DPREFIX_A=<dir>/<write_fname_a> -DPREFIX_B=<dir>/<write_fname_b>
#         -DSUFFIXES="_a.csv;_b.csv" -P compare_output_sets.cmake
# Fails naming every file pair that differs or is missing.
foreach(_v PREFIX_A PREFIX_B SUFFIXES)
    if(NOT DEFINED ${_v})
        message(FATAL_ERROR "compare_output_sets: ${_v} is required")
    endif()
endforeach()
set(_bad "")
foreach(_s IN LISTS SUFFIXES)
    set(_a "${PREFIX_A}${_s}")
    set(_b "${PREFIX_B}${_s}")
    if(NOT EXISTS "${_a}" OR NOT EXISTS "${_b}")
        list(APPEND _bad "missing: ${_a} or ${_b}")
        continue()
    endif()
    execute_process(COMMAND ${CMAKE_COMMAND} -E compare_files "${_a}" "${_b}"
                    RESULT_VARIABLE _rc OUTPUT_QUIET ERROR_QUIET)
    if(_rc EQUAL 0)
        message(STATUS "identical: ${_s}")
    else()
        list(APPEND _bad "differs: ${_s}")
    endif()
endforeach()
if(_bad)
    string(REPLACE ";" "\n  " _msg "${_bad}")
    message(FATAL_ERROR "compare_output_sets:\n  ${_msg}")
endif()
