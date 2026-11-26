// lsdtt_proj_grid_check.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// Warn when a coordinate transformation falls back to a less accurate method
// because a PROJ grid is not installed.
//
// Datum shifts such as OSGB36 -> ETRS89 (British National Grid) or NAD27 ->
// NAD83 are most accurate with a grid file. Without it PROJ silently uses a
// Helmert or ballpark transformation: BNG -> UTM is then off by 0.5 - 4 m. The
// release packages and wheels carry a small default set of grids
// (default_proj_grids() below); others can be added with
// tools/fetch_proj_grids.py. See docs/what_is_lsdtt3/note_on_grids.md.
//
// This is the one place lsdtt3 calls the PROJ C API directly (everything else
// goes through GDAL). GDAL offers no way to ask which grids an operation needs.
// The check uses the same PROJ library GDAL loads, and the same search path
// that gdal_init() pinned. If lsdtt3 was built without the PROJ headers, or
// the PROJ it was compiled against differs in major version from the one GDAL
// runs with, the check is off (proj_grid_check_status() says why) and nothing
// is printed.
// -----------------------------------------------------------------------------

#pragma once

#include <string>
#include <vector>

class OGRSpatialReference;

namespace lsdtt {

struct proj_grid_check_result
{
    bool        evaluated = false;           // the check ran
    std::string best_operation;              // the most accurate operation for the area
    double      best_accuracy_m = -1.0;      // -1: unknown
    std::vector<std::string> missing_grids;  // grids of the best operation not installed
    std::string used_operation;              // the best operation whose grids are all available
    double      used_accuracy_m = -1.0;      // -1: unknown (e.g. a ballpark transformation)
};

// The grids every release package and wheel carries. The same list is in
// tools/fetch_proj_grids.py (DEFAULT_GRIDS); tests/scripts/check_proj_grid_lists.py
// checks they agree.
[[nodiscard]] const std::vector<std::string>& default_proj_grids();

// For the -v banner: how many default grids are in the pinned data directory,
// and the names of any other grids found there.
[[nodiscard]] std::string proj_grids_description();

// "on" or "off: <reason>", for the -v banner.
[[nodiscard]] std::string proj_grid_check_status();

// Evaluate the operations from src to dst over lon/lat box (west, south,
// east, north, degrees), or over the CRSs' whole areas when lonlat_box is null.
// src and dst may be any definition PROJ reads (WKT, EPSG:xxxx).
[[nodiscard]] proj_grid_check_result check_proj_transformation_grids(
    const std::string& src, const std::string& dst, const double* lonlat_box = nullptr);

// Print a warning to stderr, once per (src, dst) pair in a process, if the
// most accurate operation needs a grid that is not installed. context names
// the step (e.g. "reprojecting the DEM"). Never throws.
void warn_if_proj_grids_missing(const OGRSpatialReference& src, const OGRSpatialReference& dst,
                                const std::string& context, const double* lonlat_box = nullptr);

// As above, with the data's extent in src coordinates (converted to lon/lat
// here), so that only operations covering the data are considered.
void warn_if_proj_grids_missing(const OGRSpatialReference& src, const OGRSpatialReference& dst,
                                const std::string& context, double src_min_x, double src_min_y,
                                double src_max_x, double src_max_y);

} // namespace lsdtt
