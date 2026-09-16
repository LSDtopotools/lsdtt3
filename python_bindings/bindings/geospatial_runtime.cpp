// geospatial_runtime.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// Exposes lsdtt3's GDAL/PROJ runtime-data controls to Python.
//
// WHY THIS IS NEEDED AT ALL
//
// lsdtt::gdal_init() pins the PROJ and GDAL data directories by probing
// relative to the RUNNING EXECUTABLE. For a command-line tool that is the tool
// itself, and <exe>/../share/proj is exactly what the release bundles stage.
//
// Inside a Python process the running executable is the INTERPRETER. So the
// probe looks at <python>/../share/proj -- which in a conda, pixi or system
// Python is that environment's PROJ data, not the wheel's. Every interpreter
// tested had such a directory. Left alone, an lsdtt3 wheel that bundles its own
// PROJ would link one PROJ and read a different one's database, which is the
// "plausible but wrong coordinates" failure the pinning exists to prevent.
//
// set_bundled_data_root() redirects the probe at the installed package
// directory, where the wheel stages share/proj and share/gdal. lsdtt3/__init__
// calls it at import, before anything can touch GDAL.
// -----------------------------------------------------------------------------

#include <string>

#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>

#include <core/lsdtt_gdal_init.hpp>

namespace nb = nanobind;

namespace lsdtt3_python
{

void bind_geospatial_runtime(nb::module_& module)
{
    module.def(
        "set_bundled_data_root",
        [](const std::string& directory) {
            return lsdtt::set_bundled_data_root(directory);
        },
        nb::arg("directory"),
        R"doc(Point the PROJ/GDAL data probe at a directory.

lsdtt3 looks for ``<directory>/share/proj/proj.db`` and
``<directory>/share/gdal/gdalvrt.xsd`` and pins whichever it finds, in
preference to the ``PROJ_DATA``, ``PROJ_LIB`` and ``GDAL_DATA`` environment
variables.

This must be called BEFORE anything else touches GDAL; ``lsdtt3/__init__``
does so at import time with the installed package directory. Returns False if
initialisation has already happened, in which case the call had no effect.
)doc");

    module.def(
        "proj_data_description",
        []() { return lsdtt::proj_data_description(); },
        R"doc(Where PROJ is getting proj.db, as a human-readable string.

One of ``"bundled: <dir>"`` (lsdtt3 pinned it), ``"PROJ_DATA=<dir>"``,
``"PROJ_LIB=<dir>"``, or ``"PROJ built-in search"``.

Put this in any bug report about projections or coordinates that look subtly
wrong -- a mismatched proj.db usually returns plausible values rather than
failing.
)doc");

    module.def(
        "gdal_data_description",
        []() { return lsdtt::gdal_data_description(); },
        "As proj_data_description(), for GDAL's own data directory.");

    module.def(
        "gdal_version",
        []() { return lsdtt::gdal_version(); },
        "GDAL release name of the library actually loaded, e.g. \"3.8.4\".");

    module.def(
        "proj_version",
        []() { return lsdtt::proj_version(); },
        R"doc(PROJ version actually in use at runtime, e.g. "9.3.1".

Reported by the loaded library rather than by the headers this was compiled
against -- that distinction is the whole point.
)doc");
}

} // namespace lsdtt3_python
