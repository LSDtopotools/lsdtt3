// lsdtt_rotation_frame_file.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// lsdtt_rotation_frame_file
//
// The rotation frame file written next to every DEM warped into an oblique
// frame (<output stem>_rotation_frame.txt). It records the frame, how it was
// defined, where it came from, and the lattice to rotate products back onto.
//
// Format: one "key: value" per line with "#" comment lines, written so the
// file is also a flat YAML mapping (yaml.safe_load reads it):
//   * text values (WKT, paths, names, version, timestamp) are single-quoted,
//     a quote inside written as '';
//   * WKT is on one line;
//   * numbers have 17 significant digits (exact double round trip) and always
//     a decimal point before any exponent (1.0e-05; YAML 1.1 reads 1e-05 as
//     a string);
//   * booleans are true / false; enum values are plain lowercase tokens.
//
// It is not read by lsdtt_param_parser, which ends a value at the first space
// and caps a line at 1023 characters. The reader here also accepts plain and
// double-quoted values (\" and \\ escapes), so a file rewritten by a YAML tool
// still reads, as long as every value stays on one line.
//
// Reading refuses a missing or newer frame_file_version, duplicate keys,
// missing required keys, bad numbers, a frame_crs_wkt that does not match the
// stored centre and azimuth, and a frame that fails the axis self-check. It
// warns on unknown keys and when the frame rebuilt from the stored end
// points, width, buffer and cell size differs from the stored one (a hand
// edit or a different PROJ); the stored values are used, since they made the
// raster.
//
// Usage
//
//   lsdtt_rotation_target target;
//   target.epsg = 32612; target.cell_size = 30.0;
//   lsdtt_rotation_provenance prov;
//   prov.created_by = "lsdtt-fetch-raster";
//   prov.rotated_dem_fname = "ridge_rot.tif";
//   write_rotation_frame("ridge_rot_rotation_frame.txt", frame, target, prov);
//
//   lsdtt_rotation_frame_file f = read_rotation_frame("ridge_rot_rotation_frame.txt");
//   if (!f.ok) { std::cerr << f.message; }
// -----------------------------------------------------------------------------

#pragma once

#include <string>
#include <vector>

#include "lsdtt_oblique_frame.hpp"

namespace lsdtt {

constexpr int rotation_frame_file_version = 1;

// The lattice that rotate-back writes onto.
struct lsdtt_rotation_target
{
    int epsg = 0;                 // 0: none, use crs_wkt
    std::string crs_wkt;          // used when epsg == 0
    double cell_size = 0.0;
    double origin_x = 0.0;        // a lattice node: (0, 0) after a fetch (-tap),
    double origin_y = 0.0;        // the source DEM's own origin after preprocessing
    std::string source_dem_fname; // preprocessing only; empty after a fetch
};

struct lsdtt_rotation_provenance
{
    std::string lsdtt3_version;   // empty on write: this build's version
    std::string created_utc;      // empty on write: now, ISO 8601 with Z
    std::string created_by;       // the tool, e.g. "lsdtt-fetch-raster"
    std::string rotated_dem_fname;
    std::string feature_fname;    // empty when the frame came from no file
    std::string feature_layer;
};

// Write the frame file. Throws std::invalid_argument for a target with
// neither an EPSG code nor WKT or with cell_size <= 0, and std::runtime_error
// if the file cannot be written.
void write_rotation_frame(const std::string& path,
                          const lsdtt_oblique_frame& frame,
                          const lsdtt_rotation_target& target,
                          const lsdtt_rotation_provenance& provenance);

struct lsdtt_rotation_frame_file
{
    bool ok = false;
    std::string message;                // why it failed, or a summary
    std::vector<std::string> warnings;
    int version = 0;
    lsdtt_oblique_frame frame;
    lsdtt_rotation_target target;
    lsdtt_rotation_provenance provenance;
};

// Read and check a frame file. Never throws.
lsdtt_rotation_frame_file read_rotation_frame(const std::string& path);

} // namespace lsdtt
