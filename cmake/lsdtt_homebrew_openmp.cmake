# SPDX-License-Identifier: AGPL-3.0-or-later
#
# lsdtt_homebrew_openmp.cmake
#
# Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
#
# Apple Clang has no OpenMP runtime or headers. Homebrew's libomp supplies
# them, but it is keg-only: it is not linked into /opt/homebrew (Apple Silicon)
# or /usr/local (Intel), so find_package(OpenMP) cannot see it and the build
# silently falls back to single-threaded. This file points FindOpenMP at the
# keg BEFORE find_package(OpenMP) runs, using the variables FindOpenMP documents
# for exactly this case (OpenMP_<lang>_FLAGS, OpenMP_<lang>_LIB_NAMES,
# OpenMP_<lib>_LIBRARY). FindOpenMP still test-compiles with them, so a wrong
# guess fails as "OpenMP not found", never as a broken build.
#
# It does nothing when:
#   - the platform is not Apple, or the compiler is not AppleClang (a conda or
#     Homebrew LLVM clang finds its own libomp);
#   - OpenMP_CXX_FLAGS is already set (the user passed their own on the
#     command line, or an earlier configure cached a working answer);
#   - no libomp keg is found.
#
# Search order for the keg: $LIBOMP_PREFIX, `brew --prefix libomp`, then the
# two standard Homebrew prefixes.
#
# Used by the top-level CMakeLists.txt and python_bindings/CMakeLists.txt.

if(APPLE AND CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang" AND NOT OpenMP_CXX_FLAGS)
    set(_lsdtt_libomp_candidates "")
    if(DEFINED ENV{LIBOMP_PREFIX})
        list(APPEND _lsdtt_libomp_candidates "$ENV{LIBOMP_PREFIX}")
    endif()
    find_program(LSDTT_BREW_EXECUTABLE brew)
    if(LSDTT_BREW_EXECUTABLE)
        # brew --prefix <formula> prints the keg path even when the formula is
        # not installed, so every candidate is checked for the files below.
        execute_process(
            COMMAND "${LSDTT_BREW_EXECUTABLE}" --prefix libomp
            RESULT_VARIABLE _lsdtt_brew_rc
            OUTPUT_VARIABLE _lsdtt_brew_libomp
            OUTPUT_STRIP_TRAILING_WHITESPACE
            ERROR_QUIET)
        if(_lsdtt_brew_rc EQUAL 0 AND _lsdtt_brew_libomp)
            list(APPEND _lsdtt_libomp_candidates "${_lsdtt_brew_libomp}")
        endif()
    endif()
    list(APPEND _lsdtt_libomp_candidates "/opt/homebrew/opt/libomp" "/usr/local/opt/libomp")

    foreach(_lsdtt_prefix IN LISTS _lsdtt_libomp_candidates)
        if(EXISTS "${_lsdtt_prefix}/include/omp.h" AND EXISTS "${_lsdtt_prefix}/lib/libomp.dylib")
            set(OpenMP_CXX_FLAGS "-Xpreprocessor -fopenmp -I${_lsdtt_prefix}/include"
                CACHE STRING "OpenMP flags for Apple Clang + Homebrew libomp")
            set(OpenMP_CXX_LIB_NAMES "omp" CACHE STRING "OpenMP libraries for Apple Clang")
            set(OpenMP_omp_LIBRARY "${_lsdtt_prefix}/lib/libomp.dylib"
                CACHE FILEPATH "Homebrew libomp")
            message(STATUS "Apple Clang: using Homebrew libomp from ${_lsdtt_prefix}")
            break()
        endif()
    endforeach()

    if(NOT OpenMP_CXX_FLAGS)
        message(STATUS "Apple Clang: no Homebrew libomp found (brew install libomp, "
                       "or set LIBOMP_PREFIX); OpenMP will not be available")
    endif()
    unset(_lsdtt_libomp_candidates)
    unset(_lsdtt_brew_rc)
    unset(_lsdtt_brew_libomp)
    unset(_lsdtt_prefix)
endif()
