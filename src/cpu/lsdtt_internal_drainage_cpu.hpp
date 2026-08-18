// lsdtt_internal_drainage_cpu.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html
//
// -----------------------------------------------------------------------------
// Internal drainage: closed-depression detection and the review round-trip.
//
// An internally drained basin is a closed depression whose water never reaches
// the DEM edge: a moulin on an ice sheet, a playa, a crater, a quarry. Standard
// flow routing destroys them, because filling or breaching raises the
// depression until it spills over its lowest saddle and the water is then
// routed out through a spill point that does not exist.
//
// This module finds those depressions and turns them into a list of
// sink_candidate records that a human can review in a GIS or a spreadsheet and
// hand back. It does NOT do any routing: forced-outlet routing lives in
// lsdtt_flowinfo (lsdtt_flowinfo_options::forced_outlet_flats) and sink-seeded
// filling in lsdtt_raster_flow_preprocessing_cpu (fill_pits_with_sinks).
//
// THE TWO ROUTING OPTIONS ARE A MATCHED PAIR IN BOTH DIRECTIONS. A forced
// outlet needs a DEM prepared by fill_pits_with_sinks(), because a plain fill
// has already drowned the pit. The converse also holds: a sink-seeded fill
// deliberately leaves genuine local minima, so building an lsdtt_flowinfo on
// that surface WITHOUT forced_outlet_flats throws "could not assign a receiver".
// Neither is usable alone.
//
// Naming trap worth stating once, loudly:
//
//   sink_candidate::area_m2 is the DEPRESSION FOOTPRINT -- the planform area of
//   the closed hollow itself. It is NOT the drainage area of the sink. On any
//   tilted background the depression intercepts the upslope strip draining into
//   it, so the internal basin is strictly LARGER than the footprint (phase 1
//   measured 6492 upslope cells against a 2828-cell crater). The depression's
//   cells are a subset of the pit's upslope set, never equal to it. The tagged
//   drainage area is a different quantity, computed by the tagging step, and
//   the two must never be conflated.
// -----------------------------------------------------------------------------

#pragma once

#include <cstdint>
#include <limits>
#include <string>
#include <vector>

#include "../core/lsdtt_data_loader.hpp"
#include "../core/lsdtt_flow_base.hpp"
#include "../core/lsdtt_raster_meta.hpp"
#include "lsdtt_raster_cpu.hpp"
#include "lsdtt_raster_flow_preprocessing_cpu.hpp"   // depression_polygons_result

namespace lsdtt {

// =============================================================================
// Detection
// =============================================================================

// -----------------------------------------------------------------------------
// sink_detection_method
//
// depression_hierarchy (default)
//   Barnes et al. (2020) depression hierarchy, via build_depression_hierarchy().
//   Depressions are found as a TREE: a lake that contains three smaller hollows
//   is a parent with three children. The pit reported is the LOWEST CELL of the
//   depression.
//
// fill_difference (legacy-compatible)
//   filled minus raw, thresholded, 8-connected components with a minimum size,
//   one sink per component. The pit reported is the DEEPEST-FILL CELL.
//
// The two pits are not the same point and the methods will legitimately
// disagree by a pixel or two. The lowest cell and the cell that had to be
// raised furthest are different questions. Do not try to reconcile them; the
// difference is a measurement to report, not a defect.
// -----------------------------------------------------------------------------
enum class sink_detection_method
{
    depression_hierarchy,
    fill_difference
};

// -----------------------------------------------------------------------------
// sink_candidate
//
// One closed depression, as offered to the reviewer.
// -----------------------------------------------------------------------------
struct sink_candidate
{
    int          sink_id         = -1;     // durable identity through review
    int          routing_sink_id = -1;     // contiguous label assigned for routing
    std::int64_t pit_flat        = -1;     // the terminal cell
    int          pit_row         = -1;
    int          pit_col         = -1;
    double       pit_x           = 0.0;    // raster CRS, cell centre
    double       pit_y           = 0.0;
    double       latitude        = 0.0;    // WGS84, for the inspection step
    double       longitude       = 0.0;
    float        pit_elevation   = 0.0f;   // raw DEM at the pit
    float        spill_elevation = 0.0f;   // out_elev of the depression
    float        depth_m         = 0.0f;   // spill_elevation - pit_elevation
    double       area_m2         = 0.0;    // DEPRESSION FOOTPRINT, not the basin
    double       volume_m3       = 0.0;    // fill volume to the spill point
    int          n_pixels        = 0;      // cells in the depression footprint
    int          hierarchy_level = 0;      // 0 = leaf; height above the leaves
    int          parent_sink_id  = -1;     // nesting, -1 if top level
    bool         accepted        = true;   // set by the review round-trip
};

// -----------------------------------------------------------------------------
// sink_detection_options
// -----------------------------------------------------------------------------
struct sink_detection_options
{
    sink_detection_method method    = sink_detection_method::depression_hierarchy;
    std::string boundary_conditions = "oooo";
    int         min_pixels_in_patch = 1;      // legacy minimum_pixels_in_patch
    float       min_depth_m         = 0.0f;   // reject shallow noise
    double      min_area_m2         = 0.0;
    double      min_volume_m3       = 0.0;
    // leaves_only: a big lake often contains several smaller closed hollows
    // within it. The depression hierarchy is a tree: the lake is a parent, the
    // hollows inside it are its children, and a "leaf" is a hollow with nothing
    // nested inside it -- the innermost, smallest closed depression at that
    // spot.
    // true  = one sink per innermost hollow (right for moulins: each swallows
    //         its own water)
    // false = also emit the enclosing lake as its own sink, which overlaps its
    //         children and double-counts area unless you handle it downstream
    //
    // Ignored by fill_difference, which has no tree: its components are
    // disjoint by construction and every sink is a leaf.
    bool        leaves_only         = true;
    int         n_threads           = 0;
};

// -----------------------------------------------------------------------------
// find_internal_drainage
//
// Detects closed depressions and returns one sink_candidate per depression,
// ordered by sink_id, which is assigned sequentially from 0 in a deterministic
// order (ascending depression label, so children precede their parents).
//
// Determinism on flat floors: a depression floor that is exactly flat over many
// cells has no unique lowest cell. The pit is resolved by taking the lowest
// elevation and, among ties, the LOWEST FLAT INDEX, so the same DEM always
// gives the same pit.
//
// The size filters (min_pixels_in_patch, min_depth_m, min_area_m2,
// min_volume_m3) are applied per depression AFTER detection, so they never
// change which depressions exist, only which are reported.
// -----------------------------------------------------------------------------
[[nodiscard]] std::vector<sink_candidate>
find_internal_drainage(const lsdtt_raster_cpu&       dem,
                       const sink_detection_options& opts);

// =============================================================================
// Review round-trip: symmetric CSV / vector I/O
//
// Every sink product can be written and read in BOTH tabular (CSV) and vector
// (FlatGeobuf, GeoPackage) form. The review step is a GIS step -- the reviewer
// is looking at satellite imagery -- so FlatGeobuf is the preferred format and
// CSV is the fallback. The two carry identical records: a .csv and a .fgb
// written from one detection run parse back to identical sink_candidate
// vectors.
// =============================================================================

enum class sink_io_format { csv, flatgeobuf, gpkg };

// -----------------------------------------------------------------------------
// resolve_sink_io_format
//
// Resolves the format from an explicit request, or from the file extension when
// the request is empty or "same_as_extension": .csv -> csv, .fgb -> flatgeobuf,
// .gpkg -> gpkg. An unknown extension or an unknown request is an error, never
// a silent default.
// -----------------------------------------------------------------------------
[[nodiscard]] sink_io_format
resolve_sink_io_format(const std::string& path, const std::string& requested);

// -----------------------------------------------------------------------------
// write_sink_candidates
//
// One writer for every format, so callers never branch on format and the CSV
// and vector paths cannot drift apart.
//
// Vector output is a POINT layer, one feature per sink, geometry at the pit
// cell centre in the raster CRS. Every sink_candidate field becomes a TYPED
// attribute -- never a string. The vector attribute model has int64 and double
// only, so ids, counts and flags are int64 and everything else is double; a
// float elevation widens to double exactly and narrows back exactly.
//
// The vector layer carries the raster CRS and needs no latitude/longitude
// columns, but writes them anyway so that the two formats round-trip to
// identical records.
//
// Doubles are written to CSV at the shortest precision that reads back to the
// same value, so a CSV round-trip is exact rather than merely close.
// -----------------------------------------------------------------------------
void write_sink_candidates(lsdtt_data_loader&                 loader,
                           const std::string&                 path,
                           const std::vector<sink_candidate>& sinks,
                           const lsdtt_raster_meta&           meta,
                           sink_io_format                     format,
                           const std::string&                 layer_name = "sink_candidates",
                           bool                               spatial_index = false);

// -----------------------------------------------------------------------------
// read_sink_candidates
//
// Reads back an edited candidate file in any supported format, through
// load_point_layer_with_attributes() -- OGR is never opened here, because
// GDAL/OGR-backed readers belong in src/core/.
//
// An `accepted` attribute (1/0, true/false, yes/no) is honoured. If the column
// is ABSENT, every feature present is treated as accepted, which reproduces the
// legacy "onlyvalid.csv" semantics: a reviewer may either flip the flag or
// simply delete the features they reject.
//
// A feature with no sink_id -- a point hand-digitised in QGIS, say -- is read
// and given a fresh sink_id continuing from the largest one in the file.
//
// POSITION is taken from, in order of preference:
//   1. the point geometry of a vector layer, transformed into the raster CRS
//   2. easting / northing attributes, assumed already in the raster CRS
//   3. latitude / longitude attributes, transformed from EPSG:4326
//   4. row / col attributes
//
// so a reviewer who drags a point in QGIS moves the sink. Note what this means
// for CSV: a CSV has no geometry of its own, and its latitude / longitude
// columns ARE its position, so for CSV input tier 1 does not apply and
// easting / northing win over latitude / longitude. That is deliberate -- the
// easting / northing pair is in the raster CRS and is therefore exact, whereas
// a lat/lon pair has been through a projection -- and it is what makes the CSV
// and vector round-trips agree field for field.
//
// The resolved position is then SNAPPED to the containing cell's centre,
// because a sink is a cell and not a point, and latitude / longitude are
// recomputed from that cell. Features that cannot be placed on the raster are
// skipped.
// -----------------------------------------------------------------------------
[[nodiscard]] std::vector<sink_candidate>
read_sink_candidates(const std::string&       path,
                     const lsdtt_raster_meta& meta,
                     sink_io_format           format,
                     const std::string&       layer_name = "");

// -----------------------------------------------------------------------------
// accepted_sinks
//
// The accepted subset, sorted by durable sink_id and assigned routing_sink_id
// contiguously from zero. sink_id and parent_sink_id remain unchanged, so an
// accepted export can still join back to the reviewed candidate layer even when
// FlatGeobuf has reordered its features.
// -----------------------------------------------------------------------------
[[nodiscard]] std::vector<sink_candidate>
accepted_sinks(const std::vector<sink_candidate>& sinks);

// =============================================================================
// Depression footprints as vectors
// =============================================================================

// -----------------------------------------------------------------------------
// sink_depressions_to_polygons
//
// One POLYGON feature per depression, keyed by sink_id so that it joins 1:1 to
// the point layer and to the stats CSV. Attributes: sink_id, pit_elevation,
// spill_elevation, depth_m, area_m2, volume_m3, n_pixels, hierarchy_level,
// parent_sink_id, accepted -- taken from the sink_candidate, so the point and
// polygon layers can never disagree.
//
// Footprints are re-derived from the DEM with the same detection settings and
// matched to the supplied sinks by pit cell, so an edited or accepted-only
// subset can be passed in and gets exactly its own polygons back. `opts` must
// therefore be the options the sinks were detected with.
//
// NESTED DEPRESSIONS: with leaves_only == false a parent lake and the hollows
// inside it are both emitted. Each cell is drawn once, into the INNERMOST
// emitted depression that contains it, so a parent's polygon is the lake with
// its children punched out as holes. Its area_m2 attribute still includes the
// children, because that is the depression's area; the polygon shows only the
// part of the lake that is not some other sink. A nested depression whose cells
// are ALL claimed by its children -- which happens whenever a lake's floor is
// exactly flat, so the lake adds no cells of its own above its hollows -- has an
// empty footprint and therefore no polygon, while still appearing in the point
// layer. Measured on a flat-floored bowl containing three cones: 5 sinks, 4
// polygons.
//
// Cost note: this re-runs the detection engine (a depression hierarchy build,
// or a fill and a component labelling) because its signature takes only a DEM.
// On a large DEM prefer to call it once, next to find_internal_drainage().
// -----------------------------------------------------------------------------
[[nodiscard]] depression_polygons_result
sink_depressions_to_polygons(const lsdtt_raster_cpu&            dem,
                             const std::vector<sink_candidate>& sinks,
                             const sink_detection_options&      opts);

// =============================================================================
// Sink-terminated routing support
// =============================================================================

enum class sink_mode
{
    forced_outlet,
    nodata_hole
};

// Returns a copy of dem with each accepted sink pit changed to nodata. The
// input raster is never modified. Sinks must already have been passed through
// accepted_sinks(), so routing_sink_id is contiguous from zero.
[[nodiscard]] lsdtt_raster_cpu
punch_sink_holes(const lsdtt_raster_cpu&            dem,
                 const std::vector<sink_candidate>& sinks);

struct sink_tagging_options
{
    float max_upslope_distance_m = std::numeric_limits<float>::max();
};

struct sink_tagging_result
{
    // Valid only where a cell is tagged. Values are routing_sink_id, including
    // zero; all other cells are nodata with mask zero. The result uses -1 as
    // nodata regardless of the DEM sentinel, so label zero remains writable.
    lsdtt_raster_cpu          sink_id_raster;
    std::vector<double>       area_m2;
    std::vector<std::int64_t> n_pixels;

    // Basins ending at outlets not represented by the supplied sinks.
    double non_sink_outlet_area_m2 = 0.0;

    // Internal-basin cells omitted by max_upslope_distance_m.
    double untagged_internal_area_m2 = 0.0;
};

// Tags each sink's disjoint upslope basin in O(total valid cells). In
// forced_outlet mode the valid pit is the seed. In nodata_hole mode the pit is
// masked and the seeds are its adjacent valid outlet cells; a seed adjacent to
// more than one hole, another void, or the raster edge is rejected as ambiguous
// because the abstract flow graph does not record which exit created an outlet.
//
// The function requires an S-vector, accepted sinks, unique pit cells, and
// routing_sink_id values exactly 0..sinks.size()-1. It enforces integer-pixel
// conservation before converting counts to areas.
[[nodiscard]] sink_tagging_result
tag_upslope_of_sinks(const lsdtt_flow_base&           flow,
                     const std::vector<sink_candidate>& sinks,
                     sink_mode                         mode,
                     const sink_tagging_options&       opts = {});

} // namespace lsdtt
