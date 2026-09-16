// lsdtt_raster_vectorise_cpu.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// Raster-to-vector conversion - CPU free functions.
//
// These functions convert labelled integer rasters to vector polygon layers.
// Each unique label value becomes one polygon feature with the label as an
// attribute.
//
// Labelled raster conventions
// ---------------------------
//   Pixel values are integer labels stored as float32 (the native lsdtt_raster_cpu
//   storage type).  Pixels with the background value (default 0) or with mask==0
//   or with the raster's no_data_value are skipped.
//
//   Label values are cast to int64 for the output attribute.  Non-integer float
//   values are rounded before casting; the caller is responsible for ensuring
//   the raster actually contains integer labels.
//
// Polygon boundary algorithm
// --------------------------
//   A double-resolution binary image is built for each label by mapping every
//   labelled pixel to a 3x3 "stamp" centred on its double-res position.
//   Neighbouring same-label pixels share edge cells in the double-res image, so
//   the Moore boundary trace of that image follows cell edges (not cell centres).
//
//   After tracing, only the odd-row x odd-col pixels of the double-res boundary
//   are kept; these correspond to corners of original pixels in the georeferenced
//   grid.  The resulting polygon ring therefore traces the actual raster cell
//   edges of the labelled region.
//
//   Interior holes (labels that form donuts) are not currently handled.  Only
//   the outermost boundary is traced.
//
// Output schema
// -------------
//   Field        Type     Description
//   label_id     int64    The integer label value
//   pixel_count  int64    Number of raster pixels carrying this label
//   area_m2      double   Pixel count x cell_size^2 in projected map units
//   source_id    int64    (optional) channel source ID from label_source_id[label]
//                         Written only when label_source_id is non-empty.
//
// Output format
// -------------
//   Prefers FlatGeobuf (.fgb) for streaming; GeoPackage (.gpkg) for
//   spatial-index-enabled writes.  See lsdtt_vectorise_options.
// -----------------------------------------------------------------------------

#pragma once

#include <cstdint>
#include <string>
#include <vector>

#include "lsdtt_raster_cpu.hpp"
#include "../core/lsdtt_data_loader.hpp"

namespace lsdtt {

// -----------------------------------------------------------------------------
// lsdtt_vectorise_options
// -----------------------------------------------------------------------------
struct lsdtt_vectorise_options
{
    std::string          layer_name    = "labelled_regions";
    vector_output_format format        = vector_output_format::flatgeobuf;
    bool                 spatial_index = false;
    float                background    = 0.0f;  // label value to treat as no-label
};

// -----------------------------------------------------------------------------
// trace_labelled_polygons
//
// Traces each unique label in a labelled integer raster to a polygon and writes
// the result to a vector file.
//
// Parameters:
//   labelled   - CPU raster whose float32 pixel values are integer labels.
//                Pixels with mask==0 or with the raster's no_data_value are
//                treated as background regardless of their value.
//   loader     - An open lsdtt_data_loader (used only for write calls and
//                to inherit the CRS WKT from the raster metadata).
//   out_path   - Full path to the output vector file (extension determines
//                nothing; format is set by options.format).
//   options    - Layer name, output format, background label, spatial index.
//
// Throws:
//   std::runtime_error if the raster metadata is invalid or not projected.
// -----------------------------------------------------------------------------
// label_source_id (optional) - if non-empty, element [L] is the channel source
// ID for flat-surface label L.  Size must be >= max_label + 1.  Pass an empty
// vector (default) to omit the source_id attribute.
void trace_labelled_polygons(
    const lsdtt_raster_cpu&              labelled,
    lsdtt_data_loader&                   loader,
    const std::string&                   out_path,
    const lsdtt_vectorise_options&       options        = {},
    const std::vector<std::int32_t>&     label_source_id = {}
);

} // namespace lsdtt
