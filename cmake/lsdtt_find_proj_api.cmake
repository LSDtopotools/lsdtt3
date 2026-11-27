# cmake/lsdtt_find_proj_api.cmake
#
# Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
# SPDX-License-Identifier: AGPL-3.0-or-later
#
# Find the PROJ C API (proj.h, libproj) for the grid-availability warning
# (src/core/lsdtt_proj_grid_check). It must come from the same installation as
# GDAL: libproj is searched beside GDAL's library first, then in the default
# paths without the active conda / pixi environment ($CONDA_PREFIX) unless GDAL
# itself is from it; a PROJ and a GDAL from different environments are refused. Mixing them (system GDAL with a
# conda libproj, because a conda env was active) put the conda lib folder on
# every link line: CMake then warned about libgomp for every target, and the
# system libgdal resolved libcurl against conda's, failing with undefined
# references to curl_*@CURL_OPENSSL_4 (v0.11.1 on a server).
#
# Sets LSDTT_PROJ_API_FOUND, LSDTT_PROJ_INCLUDE_DIR, LSDTT_PROJ_LIBRARY. When
# PROJ is not found beside GDAL the grid check is simply compiled out.
# Override with -DLSDTT_PROJ_LIBRARY=... -DLSDTT_PROJ_INCLUDE_DIR=... if needed.

set(LSDTT_PROJ_API_FOUND FALSE)

# The GDAL library file: GDAL_LIBRARY from FindGDAL, else the first file in
# GDAL_LIBRARIES, else the imported target's location (GDAL config mode).
set(_lsdtt_gdal_file "")
if(GDAL_LIBRARY AND EXISTS "${GDAL_LIBRARY}")
    set(_lsdtt_gdal_file "${GDAL_LIBRARY}")
else()
    foreach(_l IN LISTS GDAL_LIBRARIES)
        if(NOT _lsdtt_gdal_file AND EXISTS "${_l}")
            set(_lsdtt_gdal_file "${_l}")
        endif()
    endforeach()
    if(NOT _lsdtt_gdal_file AND TARGET GDAL::GDAL)
        foreach(_prop IMPORTED_LOCATION IMPORTED_LOCATION_RELEASE IMPORTED_IMPLIB IMPORTED_IMPLIB_RELEASE)
            get_target_property(_loc GDAL::GDAL ${_prop})
            if(NOT _lsdtt_gdal_file AND _loc AND EXISTS "${_loc}")
                set(_lsdtt_gdal_file "${_loc}")
            endif()
        endforeach()
    endif()
endif()

# The folder that really holds a library file (symlinks resolved: Ubuntu's
# /usr/lib/libgdal.so points into /usr/lib/x86_64-linux-gnu/).
function(_lsdtt_real_libdir file out_var)
    get_filename_component(_r "${file}" REALPATH)
    get_filename_component(_d "${_r}" DIRECTORY)
    set(${out_var} "${_d}" PARENT_SCOPE)
endfunction()

# True if `file` is inside the active conda / pixi environment ($CONDA_PREFIX).
function(_lsdtt_in_env file out_var)
    set(_in FALSE)
    if(DEFINED ENV{CONDA_PREFIX} AND NOT "$ENV{CONDA_PREFIX}" STREQUAL "")
        get_filename_component(_env "$ENV{CONDA_PREFIX}" REALPATH)
        get_filename_component(_f "${file}" REALPATH)
        string(FIND "${_f}/" "${_env}/" _pos)
        if(_pos EQUAL 0)
            set(_in TRUE)
        endif()
    endif()
    set(${out_var} ${_in} PARENT_SCOPE)
endfunction()

if(_lsdtt_gdal_file)
    _lsdtt_real_libdir("${_lsdtt_gdal_file}" _lsdtt_gdal_libdir)
    get_filename_component(_lsdtt_gdal_link_dir "${_lsdtt_gdal_file}" DIRECTORY)
    _lsdtt_in_env("${_lsdtt_gdal_file}" _lsdtt_gdal_in_env)
    # A value cached by an older configure that mixes environments (a conda
    # libproj with a system GDAL, or the reverse) is dropped and searched again.
    if(LSDTT_PROJ_LIBRARY AND EXISTS "${LSDTT_PROJ_LIBRARY}")
        _lsdtt_in_env("${LSDTT_PROJ_LIBRARY}" _lsdtt_cached_in_env)
        if(NOT _lsdtt_cached_in_env STREQUAL _lsdtt_gdal_in_env)
            unset(LSDTT_PROJ_LIBRARY CACHE)
            unset(LSDTT_PROJ_INCLUDE_DIR CACHE)
        endif()
    endif()
    # 1. Beside GDAL's own library.
    find_library(LSDTT_PROJ_LIBRARY NAMES proj proj_9 proj_8
        HINTS "${_lsdtt_gdal_libdir}" "${_lsdtt_gdal_link_dir}" NO_DEFAULT_PATH)
    # 2. Otherwise the default search, but never the active environment when
    #    GDAL is not from it (its PATH entry would otherwise put <env>/lib first).
    if(NOT LSDTT_PROJ_LIBRARY)
        set(_lsdtt_saved_ignore "${CMAKE_IGNORE_PREFIX_PATH}")
        set(_lsdtt_saved_ignore_path "${CMAKE_IGNORE_PATH}")
        if(NOT _lsdtt_gdal_in_env AND DEFINED ENV{CONDA_PREFIX} AND NOT "$ENV{CONDA_PREFIX}" STREQUAL "")
            list(APPEND CMAKE_IGNORE_PREFIX_PATH "$ENV{CONDA_PREFIX}")      # CMake >= 3.23
            list(APPEND CMAKE_IGNORE_PATH "$ENV{CONDA_PREFIX}/lib" "$ENV{CONDA_PREFIX}/lib64"
                 "$ENV{CONDA_PREFIX}/bin" "$ENV{CONDA_PREFIX}/Library/lib")  # any CMake
        endif()
        find_library(LSDTT_PROJ_LIBRARY NAMES proj proj_9 proj_8)
        set(CMAKE_IGNORE_PREFIX_PATH "${_lsdtt_saved_ignore}")
        set(CMAKE_IGNORE_PATH "${_lsdtt_saved_ignore_path}")
    endif()
    if(LSDTT_PROJ_LIBRARY)
        _lsdtt_real_libdir("${LSDTT_PROJ_LIBRARY}" _lsdtt_proj_libdir)
        get_filename_component(_lsdtt_proj_prefix "${_lsdtt_proj_libdir}" DIRECTORY)
        set(_lsdtt_inc_hints ${GDAL_INCLUDE_DIRS} ${GDAL_INCLUDE_DIR} "${_lsdtt_proj_prefix}/include")
        foreach(_i IN ITEMS ${GDAL_INCLUDE_DIRS} ${GDAL_INCLUDE_DIR})
            get_filename_component(_ip "${_i}" DIRECTORY)
            list(APPEND _lsdtt_inc_hints "${_ip}")   # /usr/include/gdal -> /usr/include
        endforeach()
        get_filename_component(_lsdtt_proj_prefix2 "${_lsdtt_proj_prefix}" DIRECTORY)
        list(APPEND _lsdtt_inc_hints "${_lsdtt_proj_prefix2}/include")   # lib/x86_64-linux-gnu -> /usr/include
        find_path(LSDTT_PROJ_INCLUDE_DIR proj.h HINTS ${_lsdtt_inc_hints} NO_DEFAULT_PATH)
        _lsdtt_in_env("${LSDTT_PROJ_LIBRARY}" _lsdtt_proj_in_env)
        if(NOT _lsdtt_proj_in_env STREQUAL _lsdtt_gdal_in_env)
            message(STATUS "PROJ grid check: off (libproj ${LSDTT_PROJ_LIBRARY} and GDAL ${_lsdtt_gdal_file} "
                           "come from different installations)")
        elseif(LSDTT_PROJ_INCLUDE_DIR)
            set(LSDTT_PROJ_API_FOUND TRUE)
        endif()
    endif()
endif()

if(LSDTT_PROJ_API_FOUND)
    message(STATUS "PROJ grid check: on (${LSDTT_PROJ_LIBRARY}, same installation as GDAL)")
elseif(NOT _lsdtt_gdal_file)
    message(STATUS "PROJ grid check: off (could not tell where GDAL's library is)")
else()
    message(STATUS "PROJ grid check: off (proj.h or libproj not found for GDAL ${_lsdtt_gdal_file})")
endif()
