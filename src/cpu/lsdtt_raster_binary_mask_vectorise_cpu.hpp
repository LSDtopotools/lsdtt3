// lsdtt_raster_binary_mask_vectorise_cpu.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2025 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html
//
// -----------------------------------------------------------------------------
// Vectorisation of labelled binary-mask components.
//
// Purpose
// -------
// This module converts labelled raster components into vector polygon records.
//
// It is separated from lsdtt_raster_binary_mask_cpu because polygonisation is
// not a pure binary-mask operation: it introduces vector geometry, schemas,
// fields, and output-record construction.
//
// Typical workflow
// ----------------
//   1. create or clean a binary mask
//   2. label components with label_binary_components_by_size()
//   3. convert labelled components to polygon records here
//
// What belongs here
// -----------------
//   - component-to-polygon conversion
//   - optional per-component scalar accumulation into output attributes
//
// What does NOT belong here
// -------------------------
//   - binary morphology itself
//   - DEM feature scoring
//   - skeleton graph logic
// -----------------------------------------------------------------------------

#pragma once

#include "lsdtt_raster_cpu.hpp"
#include "lsdtt_data_loader.hpp"

#include <vector>

namespace lsdtt {

// -----------------------------------------------------------------------------
// binary_components_to_polygon_records
//
// Converts a labelled component raster to one polygon record per component.
//
// Input:
//   label_raster  - UInt32 label raster; value 0 is background
//   scalar_raster - optional per-pixel scalar raster to accumulate by component
//   schema        - caller-supplied schema
//
// Required schema fields:
//   when scalar_raster == nullptr:
//     [0] label_id  (int64)
//     [1] area_m2   (double)
//
//   when scalar_raster != nullptr:
//     [0] label_id   (int64)
//     [1] area_m2    (double)
//     [2] scalar_sum (double)
//
// Geometry:
//   polygons are traced from component pixel outlines using the same
//   double-resolution boundary-tracing approach used elsewhere in lsdtt3.
//
// Returns:
//   one lsdtt_vector_record per component, ordered by ascending label.
// -----------------------------------------------------------------------------
[[nodiscard]] std::vector<lsdtt_vector_record>
binary_components_to_polygon_records(
    const lsdtt_raster_cpu&   label_raster,
    const lsdtt_raster_cpu*   scalar_raster,
    const lsdtt_vector_schema& schema
);

// -----------------------------------------------------------------------------
// binary_components_to_polygon_records_with_holes
//
// As above, but also detects and includes interior holes in the polygons.
//
// Returns:
//   one polygon record per labelled component, ordered by ascending label.
// -----------------------------------------------------------------------------
[[nodiscard]] std::vector<lsdtt_vector_record>
binary_components_to_polygon_records_with_holes(
    const lsdtt_raster_cpu&   label_raster,
    const lsdtt_raster_cpu*   scalar_raster,
    const lsdtt_vector_schema& schema
);

} // namespace lsdtt
