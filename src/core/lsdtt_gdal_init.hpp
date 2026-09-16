// lsdtt_gdal_init.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// One-time GDAL/PROJ initialisation.
//
// Every entry point that touches GDAL calls gdal_init() instead of
// GDALAllRegister(). It is idempotent and thread-safe, so it is cheap to call
// defensively at the top of any function.
//
// Beyond registering the drivers, it pins PROJ's and GDAL's runtime data
// directories to the copies bundled next to the executable, when they exist:
//
//   <exe_dir>/../share/proj/proj.db        (packaged bundle: bin/ + share/)
//   <exe_dir>/share/proj/proj.db           (in-build)
//   <exe_dir>/../../share/proj/proj.db     (nested install layouts)
//
// WHY THIS EXISTS
//
// PROJ finds proj.db through the PROJ_DATA / PROJ_LIB environment variables,
// and GDAL finds its own data through GDAL_DATA. Anyone who has installed a
// second GDAL generation -- conda, apt, a Python wheel -- is likely to have one
// of these exported, pointing at a proj.db that does not match the PROJ we
// linked against. The failure is not always loud: a mismatched proj.db can
// resolve a datum differently and hand back coordinates that are plausible and
// wrong. Pinning the bundled copy removes the user's environment from the
// equation entirely.
//
// Precedence: a bundled proj.db ALWAYS wins, including over an explicitly set
// PROJ_DATA. Only when no bundled copy is found do we leave PROJ's own search
// (and therefore the user's environment) untouched -- which is the correct
// behaviour for a conda-forge install, where sharing one GDAL is the point.
//
// This header deliberately does NOT include any GDAL, OGR or CPL header, so it
// can be included from anywhere without dragging ogr_spatialref.h into the
// translation unit.
// -----------------------------------------------------------------------------

#ifndef LSDTT_GDAL_INIT_HPP
#define LSDTT_GDAL_INIT_HPP

#include <string>

namespace lsdtt
{

// Register GDAL drivers and pin the bundled PROJ/GDAL data directories.
// Idempotent and thread-safe: the work happens on the first call only.
// Call this instead of GDALAllRegister().
void gdal_init();

// Override the directory that gdal_init() searches for bundled `share/proj`
// and `share/gdal`. Must be called BEFORE the first gdal_init() to have any
// effect; returns false if initialisation has already happened.
//
// This exists for the Python extension module, where executable_path() resolves
// to the Python interpreter rather than to the module, so the default probe
// would look in the wrong place (or, worse, in a conda prefix belonging to a
// different GDAL). The module can pass its own directory here at import time.
bool set_bundled_data_root(const std::string& dir);

// ---------------------------------------------------------------------------
// Diagnostics, for the -v version banner.
//
// These run gdal_init() first, so they always describe the state that is
// actually in force.
// ---------------------------------------------------------------------------

// Human-readable description of where PROJ is getting proj.db, one of:
//   "bundled: <dir>"          - we pinned it (a bundled copy was found)
//   "PROJ_DATA=<dir>"         - not bundled; the environment names a directory
//   "PROJ_LIB=<dir>"          - ditto, legacy variable
//   "PROJ built-in search"    - not bundled and no environment override
std::string proj_data_description();

// As above for GDAL's own data directory (gdalvrt.xsd and friends).
std::string gdal_data_description();

// GDAL release name, e.g. "3.8.4".
std::string gdal_version();

// PROJ release version, e.g. "9.3.1". Empty if GDAL cannot report it.
std::string proj_version();

} // namespace lsdtt

#endif // LSDTT_GDAL_INIT_HPP
