// lsdtt_gdal_init.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// See lsdtt_gdal_init.hpp for why this exists.

#include "lsdtt_gdal_init.hpp"

#include "lsdtt_exe_path.hpp"

#include <atomic>
#include <cstdlib>
#include <filesystem>
#include <mutex>
#include <string>
#include <system_error>
#include <vector>

#include <cpl_conv.h>        // CPLSetConfigOption
#include <gdal.h>            // GDALAllRegister, GDALVersionInfo
#include <ogr_srs_api.h>     // OSRSetPROJSearchPaths

namespace fs = std::filesystem;

namespace lsdtt
{

namespace
{

std::once_flag      g_init_once;
std::string         g_data_root_override;
std::string         g_pinned_proj_dir;   // "" if nothing was pinned
std::string         g_pinned_gdal_dir;   // "" if nothing was pinned

// Atomic so that set_bundled_data_root()'s "have we initialised yet?" check is
// well defined even if someone ignores the documented startup-only contract.
// The override string itself is not protected: the contract is that it is set
// before any other thread can reach gdal_init(), and a mutex here would not
// make a caller who breaks that contract correct anyway.
std::atomic<bool>   g_initialised{false};

// Read an environment variable, returning "" when unset or empty.
std::string env_or_empty(const char* name)
{
    const char* v = std::getenv(name);
    return (v && *v) ? std::string(v) : std::string();
}

// Directories to probe for bundled runtime data, in priority order.
//
// The packaged bundles put the binaries in <root>/bin and the data in
// <root>/share, so <exe_dir>/../share is the layout that matters in practice
// (see tools/package_linux_backend.sh, which stages exactly that). The other
// two cover an in-build layout and a nested install prefix.
std::vector<fs::path> candidate_roots()
{
    std::vector<fs::path> roots;

    // An explicit root REPLACES the executable-relative probe, it does not
    // precede it. A caller that sets one is telling us where its data is; if it
    // is not there, there is none. Falling through to the executable would be
    // actively wrong for the caller this exists for -- a Python extension
    // module, whose "executable" is the interpreter, so <exe>/../share/proj is
    // whatever conda or pixi environment happens to be running it. Pinning that
    // is the mismatch this module exists to prevent, not a fallback.
    if (!g_data_root_override.empty())
    {
        roots.emplace_back(g_data_root_override);
        return roots;
    }

    const fs::path exe_dir = executable_dir();
    if (!exe_dir.empty())
    {
        roots.push_back(exe_dir / "..");
        roots.push_back(exe_dir);
        roots.push_back(exe_dir / ".." / "..");
    }
    return roots;
}

// Find <root>/share/<subdir> containing `sentinel`, for the first root that
// has it. Returns "" when no bundled copy exists.
std::string find_bundled_share_dir(const std::string& subdir,
                                   const std::string& sentinel)
{
    for (const fs::path& root : candidate_roots())
    {
        const fs::path dir = root / "share" / subdir;

        std::error_code ec;
        if (!fs::exists(dir / sentinel, ec) || ec)
            continue;

        const fs::path norm = fs::weakly_canonical(dir, ec);
        return (ec ? dir : norm).string();
    }
    return "";
}

void do_init()
{
    // The one place that calls GDALAllRegister() directly. Everywhere else in
    // lsdtt3 calls gdal_init(); calling it from here would re-enter
    // std::call_once and deadlock.
    GDALAllRegister();

    // ---------------------------------------------------------------------
    // PROJ. A bundled proj.db always wins, including over PROJ_DATA /
    // PROJ_LIB: paths set through OSRSetPROJSearchPaths take precedence over
    // PROJ's environment-variable search.
    // ---------------------------------------------------------------------
    g_pinned_proj_dir = find_bundled_share_dir("proj", "proj.db");
    if (!g_pinned_proj_dir.empty())
    {
        const char* paths[] = { g_pinned_proj_dir.c_str(), nullptr };
        OSRSetPROJSearchPaths(paths);
    }

    // ---------------------------------------------------------------------
    // GDAL's own data directory has the identical failure mode: a stale
    // GDAL_DATA from another installation. A config option set here overrides
    // the environment variable.
    //
    // gdalvrt.xsd is the sentinel because it has been present in every GDAL
    // data directory for the versions we support, and python_bindings already
    // uses it to identify the directory.
    // ---------------------------------------------------------------------
    g_pinned_gdal_dir = find_bundled_share_dir("gdal", "gdalvrt.xsd");
    if (!g_pinned_gdal_dir.empty())
        CPLSetConfigOption("GDAL_DATA", g_pinned_gdal_dir.c_str());

    g_initialised = true;
}

} // namespace

void gdal_init()
{
    std::call_once(g_init_once, do_init);
}

bool set_bundled_data_root(const std::string& dir)
{
    if (g_initialised)
        return false;
    g_data_root_override = dir;
    return true;
}

std::string proj_data_description()
{
    gdal_init();

    if (!g_pinned_proj_dir.empty())
        return "bundled: " + g_pinned_proj_dir;

    const std::string proj_data = env_or_empty("PROJ_DATA");
    if (!proj_data.empty())
        return "PROJ_DATA=" + proj_data;

    const std::string proj_lib = env_or_empty("PROJ_LIB");
    if (!proj_lib.empty())
        return "PROJ_LIB=" + proj_lib;

    return "PROJ built-in search";
}

std::string gdal_data_description()
{
    gdal_init();

    if (!g_pinned_gdal_dir.empty())
        return "bundled: " + g_pinned_gdal_dir;

    const std::string gdal_data = env_or_empty("GDAL_DATA");
    if (!gdal_data.empty())
        return "GDAL_DATA=" + gdal_data;

    return "GDAL built-in search";
}

std::string gdal_version()
{
    gdal_init();
    const char* v = GDALVersionInfo("RELEASE_NAME");
    return v ? std::string(v) : std::string();
}

std::string proj_version()
{
    gdal_init();

    // OSRGetPROJVersion reports the PROJ that GDAL is actually using at
    // runtime, which is the number worth printing: it is the one that has to
    // match the proj.db we pinned. Reading PROJ_VERSION_* from proj.h would
    // report the version this file was COMPILED against, which is precisely
    // the assumption this module exists to stop making.
    int major = 0, minor = 0, patch = 0;
    OSRGetPROJVersion(&major, &minor, &patch);
    if (major == 0 && minor == 0 && patch == 0)
        return "";

    return std::to_string(major) + "." + std::to_string(minor) + "." +
           std::to_string(patch);
}

} // namespace lsdtt
