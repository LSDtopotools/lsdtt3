// lsdtt_divide_migration_cpu.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html

#pragma once

#include <cstddef>
#include <cstdint>
#include <string>
#include <utility>
#include <vector>

#include "../core/lsdtt_types.hpp"
#include "../utils/lsdtt_raster_geometry.hpp"
#include "lsdtt_raster_cpu.hpp"

namespace lsdtt
{

// Defined in lsdtt_chi_extraction_cpu.hpp; only referenced here.
struct lsdtt_chi_collection_result;

struct lsdtt_outlet_record
{
    int point_id = -1;
    int row = -1;
    int col = -1;
    std::int64_t input_flat = -1;
    std::int64_t snapped_flat = -1;
    bool snapped = false;

    int unit_id = -1;
    std::string unit_key = "";
    std::string cluster_value = "";
};

struct lsdtt_unit_record
{
    int unit_id = -1;
    std::string unit_key = "";
    std::vector<int> point_ids;
    std::vector<std::int64_t> snapped_outlet_flats;
};

struct lsdtt_divide_segment_record
{
    std::string pair_id = "";
    int unit_a = -1;
    int unit_b = -1;
    // Source-basin junction numbers. These describe a single basin pair, so they
    // are only meaningful while the record holds one extracted fragment. Once
    // fragments from several basin pairs are merged (merge_divide_fragments)
    // they are set to -1 rather than left reporting one contributor's junction
    // as though it described the whole line. n_fragments says how many fragments
    // the record holds, so -1 from merging is distinguishable from -1 unset.
    std::int32_t basin_junction_a = -1;
    std::int32_t basin_junction_b = -1;
    int segment_id = -1;
    int n_fragments = 1;
    std::vector<lsdtt_xy> points;
    double length_m = 0.0;
};

struct lsdtt_channel_head_metric_record
{
    int head_id = -1;
    int outlet_id = -1;
    int unit_id = -1;
    std::string pair_id_nearest = "";

    double head_x = 0.0;
    double head_y = 0.0;
    double head_elevation_m = 0.0;
    double chi = 0.0;

    double nearest_divide_x = 0.0;
    double nearest_divide_y = 0.0;
    double divide_elevation_m = 0.0;
    double distance_to_divide_m = 0.0;
    double relief_m = 0.0;
    double gradient_to_divide = 0.0;
};

// -----------------------------------------------------------------------------
// lsdtt_pair_summary_record
//
// One row of the side-vs-side comparison across a single divide pair.
//
// SIGN CONVENTION -- read this before interpreting any delta.
//   Side A is ALWAYS the LOWER-numbered unit of the pair and side B the higher,
//   the same ordering make_unit_pair_id uses to build "<min>|<max>". Every delta
//   is therefore
//
//       delta_<metric> = mean_<metric>_a - mean_<metric>_b   (A minus B)
//
//   and the sides can never be silently swapped: unit_a < unit_b always holds.
//   Under the standard interpretation the aggressor is the side whose heads sit
//   at LOWER elevation with HIGHER gradient and relief and LOWER chi, so a
//   negative delta_head_elevation, a positive delta_gradient, a positive
//   delta_relief and a negative delta_chi all point to A as the aggressor, and
//   the divide is expected to migrate towards B.
//
// SAMPLE SIZES.
//   n_heads_a / n_heads_b are the numbers of channel heads assigned to each
//   side of the pair. They are NOT the sample size of every mean: relief,
//   gradient and divide elevation are legitimately NaN for some heads (nearest
//   divide point on nodata; zero distance gives no gradient), and those heads
//   are skipped rather than being allowed to poison the mean. Each metric
//   therefore carries its own valid-value count, n_<metric>_a / n_<metric>_b,
//   which is the n its mean AND its uncertainty were computed from. This
//   extends the original schema, which named only
//   n_heads_a / n_heads_b.
// -----------------------------------------------------------------------------
struct lsdtt_pair_summary_record
{
    std::string pair_id = "";
    int unit_a = -1;
    int unit_b = -1;

    // Heads assigned to each side, valid metric values or not.
    std::int64_t n_heads_a = 0;
    std::int64_t n_heads_b = 0;

    double mean_head_elevation_a = 0.0;
    double mean_head_elevation_b = 0.0;
    double delta_head_elevation = 0.0;
    std::int64_t n_head_elevation_a = 0;
    std::int64_t n_head_elevation_b = 0;

    double mean_gradient_a = 0.0;
    double mean_gradient_b = 0.0;
    double delta_gradient = 0.0;
    std::int64_t n_gradient_a = 0;
    std::int64_t n_gradient_b = 0;

    double mean_relief_a = 0.0;
    double mean_relief_b = 0.0;
    double delta_relief = 0.0;
    std::int64_t n_relief_a = 0;
    std::int64_t n_relief_b = 0;

    double mean_chi_a = 0.0;
    double mean_chi_b = 0.0;
    double delta_chi = 0.0;
    std::int64_t n_chi_a = 0;
    std::int64_t n_chi_b = 0;

    std::string uncertainty_method = "";
    double uncertainty_head_elevation = 0.0;
    double uncertainty_gradient = 0.0;
    double uncertainty_relief = 0.0;
    double uncertainty_chi = 0.0;
};

// -----------------------------------------------------------------------------
// Uncertainty reported alongside each pair-summary delta.
//
// DivideTools' AcrossDivide.m reports a dispersion for EACH side and then asks
// whether the two intervals overlap. The schema here has one uncertainty column
// per metric, so each option below is the equivalent statistic OF THE DELTA
// rather than of a single side. The formulae for the per-side ingredients
// (sample standard deviation with the n-1 denominator, standard error
// s/sqrt(n), a 1000-replicate bootstrap of the mean, Welch's unequal-variance
// t-test) are taken from AcrossDivide.m; combining them into one number per
// metric is ours, and is stated here rather than implied to be DivideTools'.
//
//   std_dev    Pooled standard deviation of the two sides,
//                sqrt( ((na-1) sa^2 + (nb-1) sb^2) / (na + nb - 2) ).
//              A spread of the underlying head populations, NOT of the delta:
//              it does not shrink as more heads are added.
//
//   std_err    Standard error OF THE DIFFERENCE OF MEANS (the default),
//                sqrt( sa^2/na + sb^2/nb ).
//              This is the one to compare a delta against: |delta| > 2 * this
//              is roughly a 95% separation of the two side means.
//
//   bootstrap  Half-width of the 95% confidence interval of the BOOTSTRAP
//              DISTRIBUTION OF THE DELTA. Both sides are resampled with
//              replacement bootstrap_replicates times, the delta of means is
//              recorded for each replicate, and the reported value is
//                0.5 * ( |q_0.025 - mean_boot| + |q_0.975 - mean_boot| ),
//              the same statistic AcrossDivide.m's bootCI helper returns for a
//              single side. It is a ~95% interval, so it is roughly 1.96 times
//              std_err for a large, well-behaved sample: do not compare a
//              bootstrap column with a std_err column.
//
//   ttest      The two-sided p-value of Welch's unequal-variance t-test on the
//              difference of means, matching AcrossDivide.m's
//              ttest2(...,'Vartype','unequal'). This column is a PROBABILITY,
//              not a dispersion: small means the two sides differ. The record's
//              uncertainty_method field states which of these four produced the
//              numbers, so the column is never ambiguous.
//
// SMALL SAMPLES. Every method needs a within-side spread, so all four require
// at least TWO valid values on BOTH sides of the metric in question and return
// NaN otherwise. A single-head side and an empty side therefore give a NaN
// uncertainty rather than a meaningless zero, and no method ever divides by
// zero. A pair in which one side has no heads at all still emits a row, with
// that side's counts at zero and its mean, the delta and the uncertainty NaN:
// the divide exists in the geometry output, and a row saying so is more use
// than a silent absence.
// -----------------------------------------------------------------------------
enum class lsdtt_pair_uncertainty_method
{
    std_dev,
    std_err,
    bootstrap,
    ttest
};

// Parses "std_dev" / "std_err" / "bootstrap" / "ttest". Unrecognised text
// yields std_err (the default), which the CLI never relies on because it
// validates the string first.
[[nodiscard]] lsdtt_pair_uncertainty_method parse_pair_uncertainty_method(
    const std::string& method_text);

[[nodiscard]] std::string to_string(lsdtt_pair_uncertainty_method method);

struct lsdtt_pair_summary_options
{
    lsdtt_pair_uncertainty_method uncertainty_method =
        lsdtt_pair_uncertainty_method::std_err;

    // Only read when uncertainty_method is bootstrap.
    int           bootstrap_replicates = 1000;
    std::uint32_t bootstrap_seed = 5489;
};

// -----------------------------------------------------------------------------
// compute_pair_summaries
//
// Turns per-head records into one side-vs-side row per divide pair, mirroring
// the interpretation table of DivideTools' AcrossDivide.m.
//
// THE PAIR SET COMES FROM THE DIVIDES, NOT FROM THE HEADS. Every pair present
// in divide_segments gets exactly one row, whether or not any head was matched
// to it, so this table and the divide vector layer can never disagree about
// which pairs exist. A pair that attracted no heads at all still emits a row,
// with both counts 0 and every mean, delta and uncertainty NaN -- the same
// convention an empty single side already follows, taken one step further. Two
// segments sharing a pair_id (the same unit pair meeting along two separate
// stretches of divide) still give ONE row, not two.
//
// Grouping: heads are grouped on pair_id_nearest, and only onto pairs the
// segments already established. A head with an EMPTY pair_id_nearest matched no
// divide and takes no part in any row; a head naming a pair absent from
// divide_segments never invents one, which cannot arise from
// compute_nearest_divide_metrics since it reads every pair_id_nearest straight
// off one of these same segments. Within a pair, a head belongs to side A when
// its unit_id equals the pair's lower unit and to side B when it equals the
// higher unit; the units are read from the pair_id itself ("<min>|<max>"), so
// side A is well defined even when no head falls on it. A head whose unit_id is
// neither is not counted on either side -- it cannot arise from
// compute_nearest_divide_metrics, which only ever matches a head to a divide
// bounding its own unit.
//
// See lsdtt_pair_summary_record for the A-minus-B sign convention and the
// per-metric sample counts, and lsdtt_pair_uncertainty_method for what each
// uncertainty option reports and how small samples are handled.
//
// Determinism: rows come out ordered by pair_id. The bootstrap draws from a
// std::mt19937 seeded, through a std::seed_seq, from bootstrap_seed together
// with a fixed hash of the pair_id and the metric index -- never from the
// iteration position -- so the same seed and the same input give bit-identical
// output whatever the order in which pairs or metrics are visited.
// -----------------------------------------------------------------------------
[[nodiscard]] std::vector<lsdtt_pair_summary_record> compute_pair_summaries(
    const std::vector<lsdtt_divide_segment_record>&      divide_segments,
    const std::vector<lsdtt_channel_head_metric_record>& head_metrics,
    const lsdtt_pair_summary_options&                    options);

// -----------------------------------------------------------------------------
// merge_divide_fragments
//
// A unit made of several basins produces one divide fragment per adjacent basin
// pair, and all of those fragments share one pair_id. This joins the fragments
// that are contiguous in space into single polylines, so that per-divide length
// and along-divide profiles are available without dissolving by pair_id
// downstream.
//
// Contiguity criterion: two fragments are joined when an endpoint of one lies
// within ONE CELL DIAGONAL (cell_size_m * sqrt(2)) of an endpoint of the other.
// An exact shared-vertex test was rejected: each fragment is traced from a
// different basin's outline, and two outlines meeting at a triple point need not
// stop on the same outline node, so the two fragment ends can legitimately be
// offset by a fraction of a cell. One cell diagonal is the largest gap that can
// still be a rounding/tracing artefact rather than real separation -- ends
// further apart than that have at least one intervening cell between them and
// are genuinely disjoint stretches of divide.
//
// Not over-merging is the point of this routine. Two units can meet along two or
// more separate stretches of divide (they wrap around a third unit, say) and
// those must stay separate features, so only endpoint-to-endpoint proximity
// joins fragments; sharing a pair_id never does on its own.
//
// Mechanics:
//   - fragments may arrive in any order and with either orientation, so chains
//     are grown from both ends and fragments are reversed as needed,
//   - the shared vertex at a join is written once, not twice,
//   - a chain that closes on itself (first vertex within tolerance of the last)
//     is a ring: it is closed and no longer extended, so no chain can loop,
//   - length_m and the vertex count are recomputed from the merged geometry,
//     never summed from the inputs,
//   - segment_id is reassigned over the merged set so it stays unique and dense,
//   - basin_junction_a / basin_junction_b survive only on records built from a
//     single fragment; a merged record reports -1 for both and n_fragments > 1.
//
// Fragments with different pair_ids are never joined to each other.
// -----------------------------------------------------------------------------
[[nodiscard]] std::vector<lsdtt_divide_segment_record> merge_divide_fragments(
    const std::vector<lsdtt_divide_segment_record>& fragments,
    double                                          cell_size_m);

// -----------------------------------------------------------------------------
// find_adjacent_basin_pairs
//
// Which of the supplied basins touch each other on the grid. Divide extraction
// is expensive and only ever succeeds for basins that share a boundary, so this
// replaces the naive all-vs-all pair loop over N basins with one pass over the
// basin cells.
//
// CONNECTIVITY IS EIGHT-WAY, and that is not a free choice. A divide is traced
// by walking a basin outline and asking, at each outline NODE, which of the four
// cells meeting at that node lie outside the basin and drain to the other basin
// (lsdtt_basin_outline_cpu::get_pixels_around_basin_outline_node). Two cells
// that meet only at a corner therefore both touch the same outline node, so a
// pair of basins that touch only diagonally can still yield a divide. Four-way
// contact would miss exactly those pairs.
//
// This is a SUPERSET of the pairs that can produce a divide, which is what makes
// it safe to skip everything else. If get_interfluvial_divide_as_ridgeline
// returns geometry for a pair, some cell outside basin i shares a node with a
// cell of basin i and drains into the other basin's junction, i.e. it lies in
// the other basin's upslope set; those two cells are then eight-connected, so
// the pair is found here. Pairs found here that produce no divide simply fall
// out of the extraction exactly as they did before.
//
// basin_cells[b] holds the flat indices of basin b (the zero-copy upslope slice
// of the basin), and basins may overlap: two basins covering the same cell are
// reported as adjacent, so a nested pair is still returned here and is left to
// the caller's own nesting rule to reject.
//
// Returns unordered pairs (i, j) of indices into basin_cells with i < j, sorted
// ascending by i then j -- the same order the old nested loop visited them in,
// so downstream segment ordering and segment_id assignment do not move.
// -----------------------------------------------------------------------------
[[nodiscard]] std::vector<std::pair<int, int>> find_adjacent_basin_pairs(
    const std::vector<lsdtt_flat_index_slice>& basin_cells,
    int                                        n_rows,
    int                                        n_cols,
    bool                                       periodic_ns = false,
    bool                                       periodic_ew = false);

// -----------------------------------------------------------------------------
// compute_nearest_divide_metrics
//
// For every channel head, find the nearest location on the divide geometry and
// fill the nearest-divide fields of the record:
//   pair_id_nearest, nearest_divide_x, nearest_divide_y, divide_elevation_m,
//   distance_to_divide_m, relief_m, gradient_to_divide.
//
// Rules:
//   - only divide segments whose pair involves the head's own unit_id are
//     considered, so a head is never matched to a divide between two other
//     units. A head whose unit takes part in no divide keeps its NaN fields
//     and an empty pair_id_nearest,
//   - the nearest location is the closest point on a divide *segment*
//     (point-to-line-segment projection), not merely the closest vertex,
//   - divide_elevation_m is sampled from elevation (normally the filled/routed
//     DEM) at the cell containing the nearest point. Out-of-bounds or nodata
//     cells leave the elevation-derived fields as NaN,
//   - gradient_to_divide is NaN when the head sits on (or within rounding of)
//     the divide, so the division by zero never happens. Negative relief is a
//     legitimate result and is passed through unchanged.
//
// The search uses a uniform-grid index built per unit, so cost is close to
// linear in the number of heads for well-distributed divides.
//
// Returns the number of heads that were matched to a divide.
// -----------------------------------------------------------------------------
[[nodiscard]] std::size_t compute_nearest_divide_metrics(
    const std::vector<lsdtt_divide_segment_record>& divide_segments,
    const lsdtt_raster_cpu&                         elevation,
    std::vector<lsdtt_channel_head_metric_record>&  head_metrics);

// -----------------------------------------------------------------------------
// apply_common_base_level_to_chi
//
// Impose a common base level on an outlet-seeded chi collection, mirroring the
// intent of DivideTools' SetOutlet.m: every chi network is truncated below
// base_elevation_m and chi is re-zeroed at the lowest retained node of each
// remaining channel fragment. Comparisons across basins are then made from the
// same elevation datum rather than from each basin's own outlet.
//
// Elevations are read from elevation (this must be the filled/routed DEM so
// that flow paths are monotonic and the retained node set stays upslope
// closed). basin.outlet_flat is left untouched so caller-side outlet -> unit
// maps remain valid.
//
// Returns the number of chi nodes removed.
// -----------------------------------------------------------------------------
std::size_t apply_common_base_level_to_chi(
    lsdtt_chi_collection_result& chi,
    const lsdtt_raster_cpu&      elevation,
    double                       base_elevation_m);

} // namespace lsdtt
