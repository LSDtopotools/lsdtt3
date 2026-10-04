// lsdtt_subbasin_partition_cpu.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// Incremental sub-basin partition, and the per-sub-basin set-up quantities
// derived from it.
//
// MODEL-AGNOSTIC. Written for the simplified S-HYPE column, whose spatial unit
// is the HYPE sub-basin (an incremental drainage area with one outlet and one
// downstream sub-basin), but nothing here knows about HYPE parameters. Any
// semi-distributed model that needs "the area between this gauge and the next
// one up" can use it.
//
// NESTED VERSUS INCREMENTAL
// -------------------------
// lsdtt_basin_cpu / lsdtt_junction_basin_cpu give NESTED basins: a basin holds
// everything upstream of its outlet, including any other selected basin. A
// sub-basin here is INCREMENTAL: a pixel belongs to the FIRST selected outlet
// at or downstream of it. Gauges along one river therefore form a chain of
// sub-basins, each linked to the sub-basin that contains its outlet's
// receiver ("maindown" in HYPE's GeoData).
//
// SINGLE-RECEIVER FLOW ONLY. The partition needs every pixel to drain to
// exactly one sub-basin, so it is defined on the D8 lsdtt_flow_base. D-inf and
// MFD split a pixel between receivers and cannot be used.
//
// WHAT IS HERE
// ------------
//   outlets         subbasin_outlets_from_user_flats       (outlet list)
//                   subbasin_outlets_by_stream_order       (Strahler order)
//   partition       partition_incremental_subbasins
//                   get_subbasin_id_raster
//   GeoData-like    compute_subbasin_geodata               (area, elevation
//                                                           mean/std, slope)
//   classes         compute_subbasin_class_fractions       (fraction, dh)
//   river lengths   compute_subbasin_river_lengths         (rivlen, loc_rivlen)
//
// NO ANGLES. Slope is the dimensionless gradient S (m/m) throughout, from
// compute_gradient_magnitude_3x3_square. The only conversion is HYPE's percent,
// 100 * S. No slope angle in degrees or radians is computed anywhere here.
//
// NODATA. Every per-sub-basin statistic excludes pixels that are nodata in the
// raster being summarised, and reports how many pixels it used. Area is NOT
// reduced by nodata in a summarised raster: it is the partition's pixel count,
// i.e. every valid flow cell assigned to the sub-basin. A statistic with no
// valid pixels in a sub-basin is NaN.
//
// See docs/design_and_algorithms/hydrological_modelling.md
// (S-HYPE) for how the HYPE model uses these.
// -----------------------------------------------------------------------------

#pragma once

#include <cstdint>
#include <vector>

#include "../core/lsdtt_flow_base.hpp"
#include "../core/lsdtt_junction_network_base.hpp"
#include "../core/lsdtt_raster_meta.hpp"
#include "../core/lsdtt_raster_result.hpp"
#include "lsdtt_raster_cpu.hpp"

namespace lsdtt
{

// -----------------------------------------------------------------------------
// Outlet selection
//
// Both return basin outlet flat indices, nested or not, ready for
// partition_incremental_subbasins. They are thin wrappers on
// lsdtt_junction_basin_manager_cpu so the sub-basin set-up selects outlets
// exactly as lsdtt-basin-extraction does. Snapping user outlets to the channel
// is the caller's job (basin_outlet_snapping_method in the CLI).
//
// subbasin_outlets_from_user_flats:
//   select_junction_basins_from_basin_outlets -> get_basin_outlet_flat_indices.
//   Every flat must lie on the junction network's channel; the junction basin
//   constructor throws otherwise. Order and duplicates are passed through
//   unchanged (the partition deduplicates).
//
// subbasin_outlets_by_stream_order:
//   select_all_junction_basins -> filter_by_stream_order(min, max). Pass
//   min_order == max_order for one order; -1 leaves a bound open. Note what
//   "all junction basins" means: each junction's basin outlet is the pixel
//   immediately upstream of the NEXT downstream junction, and a channel pixel
//   at base level is itself a junction whose basin outlet is that pixel. So a
//   trunk ending at base level contributes two outlets of its order: the pixel
//   above the base-level pixel, and the base-level pixel itself.
//   The base-level pixel's own sub-basin would be only the edge strip draining
//   straight to it, which is degenerate (its 3x3 gradient is all nodata). So
//   the returned list MERGES it (decided 2026-10-04): the selected outlet
//   draining directly into a base-level outlet (the largest, if several) is
//   dropped, making that area and the strip one sub-basin with its outlet at
//   the base-level pixel. A base-level outlet with nothing selected above it is
//   dropped, and its pixels are left unassigned.
// -----------------------------------------------------------------------------
[[nodiscard]] std::vector<std::int64_t> subbasin_outlets_from_user_flats(
    const lsdtt_flow_base&             flow,
    const lsdtt_junction_network_base& junction_network,
    const std::vector<std::int64_t>&   outlet_flats
);

[[nodiscard]] std::vector<std::int64_t> subbasin_outlets_by_stream_order(
    const lsdtt_flow_base&             flow,
    const lsdtt_junction_network_base& junction_network,
    std::int32_t                       min_order,
    std::int32_t                       max_order
);

// -----------------------------------------------------------------------------
// subbasin_partition_result
//
// Sub-basin index s runs 0..n_subbasins()-1, in the order of the outlets given
// to the partition after duplicates are removed (first occurrence kept).
//
// subbasin_id has full-grid flat size n_rows * n_cols:
//   >= 0  the pixel's sub-basin;
//   -1    nodata, or a valid pixel that reaches base level without passing a
//         selected outlet (counted in n_unassigned).
//
// downstream[s] is the sub-basin containing the receiver of outlet_flat[s]
// (HYPE's maindown), or -1 when that receiver is base level or unassigned.
//
// upstream_to_downstream lists every sub-basin once, each after all sub-basins
// that drain into it, so one forward sweep over it routes the network.
// -----------------------------------------------------------------------------
struct subbasin_partition_result
{
    lsdtt_raster_meta meta;

    std::vector<std::int32_t> subbasin_id;

    std::vector<std::int64_t> outlet_flat;
    std::vector<std::int64_t> n_pixels;
    std::vector<std::int32_t> downstream;
    std::vector<std::int32_t> upstream_to_downstream;

    std::int64_t n_valid_pixels     = 0;
    std::int64_t n_unassigned       = 0;
    std::int64_t n_duplicate_outlets = 0;

    [[nodiscard]] std::int32_t n_subbasins() const noexcept
    {
        return static_cast<std::int32_t>(outlet_flat.size());
    }
};

// -----------------------------------------------------------------------------
// partition_incremental_subbasins
//
// One sweep from base level upward: each base-level outlet's upslope s_vector
// is visited in preorder (receiver before donor), and every pixel takes its
// receiver's sub-basin id unless it is itself a selected outlet, which starts a
// new id. Cost O(N). Requires the flow object's s_vector.
//
// DUPLICATE OUTLETS are DEDUPLICATED, not rejected: snapping two gauges to one
// pixel is an ordinary outcome, and the two would describe the same
// sub-basin. The first occurrence keeps its position; the number dropped is
// n_duplicate_outlets. Out-of-range or nodata outlet flats throw
// std::invalid_argument.
// -----------------------------------------------------------------------------
[[nodiscard]] subbasin_partition_result partition_incremental_subbasins(
    const lsdtt_flow_base&           flow,
    const std::vector<std::int64_t>& outlet_flats
);

// Sub-basin id as a raster. Assigned pixels carry their id; unassigned and
// nodata pixels are nodata (mask 0). Ids are exact in float up to 2^24.
[[nodiscard]] lsdtt_raster_result get_subbasin_id_raster(
    const subbasin_partition_result& partition
);

// -----------------------------------------------------------------------------
// compute_subbasin_geodata
//
// HYPE GeoData columns per sub-basin, all over the INCREMENTAL pixel set:
//
//   area_m2             n_pixels * cell_size^2 (every assigned pixel)
//   elev_mean_m         mean of dem over its valid pixels
//   elev_std_m          POPULATION standard deviation (divisor n) of the same
//                       pixels. The sub-basin's pixels are the whole
//                       population being described, not a sample of it
//                       (see lsdtt_statistics.hpp). HYPE uses elev_std only
//                       for fractional snow cover.
//   slope_mean_percent  100 * mean of gradient (dimensionless, m/m) over its
//                       valid pixels. Pass the output of
//                       compute_gradient_magnitude_3x3_square, which is
//                       nodata on the grid edge and next to nodata; those
//                       pixels are excluded.
//
// n_valid_elev / n_valid_gradient record how many pixels each mean used.
// dem and gradient must be on the partition's grid (same n_rows, n_cols).
// -----------------------------------------------------------------------------
struct subbasin_geodata_result
{
    std::vector<double>       area_m2;
    std::vector<double>       elev_mean_m;
    std::vector<double>       elev_std_m;
    std::vector<double>       slope_mean_percent;
    std::vector<std::int64_t> n_valid_elev;
    std::vector<std::int64_t> n_valid_gradient;
};

[[nodiscard]] subbasin_geodata_result compute_subbasin_geodata(
    const subbasin_partition_result& partition,
    const lsdtt_raster_cpu&          dem,
    const lsdtt_raster_cpu&          gradient
);

// -----------------------------------------------------------------------------
// compute_subbasin_class_fractions
//
// Per (sub-basin, class) pair present, HYPE's slc_n and dhslc_n:
//
//   fraction  class pixels / sub-basin pixels with a valid class code.
//             Pixels that are nodata in class_raster are left out of both
//             numerator and denominator, so a sub-basin's fractions sum to 1.
//   dh_m      mean dem over the class's pixels (valid dem only) minus the
//             sub-basin elev_mean (mean dem over ALL its valid-dem pixels,
//             as compute_subbasin_geodata). NaN if none of the class's pixels
//             has a valid elevation. When class_raster and dem share a mask,
//             sum(fraction * dh_m) over a sub-basin is 0 to rounding.
//
// class_raster values must be integral (they are class codes stored as float,
// as lsdtt_landcover_scenario_cpu reads them); a non-integral valid value
// throws std::invalid_argument.
//
// Layout: compressed rows, as the HYPE model's unit array.
// Entries for sub-basin s are entries[offset[s]] .. entries[offset[s+1]-1],
// sorted by class code. Absent classes have no entry (fraction 0).
// -----------------------------------------------------------------------------
struct subbasin_class_entry
{
    std::int32_t subbasin   = -1;
    std::int32_t class_code = 0;
    std::int64_t n_pixels   = 0;
    double       fraction   = 0.0;
    double       dh_m       = 0.0;
};

struct subbasin_class_fractions_result
{
    std::vector<std::int32_t>         class_codes;   // sorted unique codes present
    std::vector<std::int64_t>         offset;        // size n_subbasins + 1
    std::vector<subbasin_class_entry> entries;
    std::vector<std::int64_t>         n_valid_class; // per sub-basin denominator
};

[[nodiscard]] subbasin_class_fractions_result compute_subbasin_class_fractions(
    const subbasin_partition_result& partition,
    const lsdtt_raster_cpu&          class_raster,
    const lsdtt_raster_cpu&          dem
);

// -----------------------------------------------------------------------------
// compute_subbasin_river_lengths
//
// FD is get_distance_from_outlet_raster(): D8 flow-path
// length to base level, cardinal step = cell_size, diagonal = cell_size*sqrt 2.
//
// TRUNK (the sub-basin's main river), from its upstream end p_in to its
// outlet p_out:
//
//   * headwater sub-basin (nothing drains into it): start at the outlet and
//     repeatedly step to the donor that lies IN THIS SUB-BASIN and has the
//     largest contributing pixel count (include-self). Stop when no such donor
//     has a count >= channel_threshold_contributing_pixels; that pixel is p_in.
//     The threshold therefore sets where the trunk's channel head is. If the
//     outlet's own in-sub-basin donors are all below threshold, the trunk is
//     the outlet alone and rivlen is 0.
//
//   * sub-basin with upstream sub-basins: p_in is the receiver of the outlet of
//     the upstream sub-basin with the largest contributing pixel count at its
//     outlet (ties: lower sub-basin index), and the trunk is the D8 path from
//     p_in down to p_out. This is the definition used, because the main
//     river's job in HYPE is to carry the upstream inflow; the threshold plays
//     no part. Where the largest upstream inflow is also the largest-area
//     branch at every confluence, this is exactly where the in-sub-basin walk
//     above stops; where an in-sub-basin tributary outgrows that inflow, the
//     two differ and the inflow definition wins.
//
// Donor ties in the walk go to the lower flat index.
//
//   rivlen_m      FD(p_in) - FD(p_out), accumulated in double along the trunk.
//   loc_rivlen_m  mean over the sub-basin's pixels of the flow distance to the
//                 nearest downstream trunk pixel (trunk pixels contribute 0).
//                 Computed as HAND with FD in place of elevation
//                 (get_height_above_channel_mask_raster with the trunk mask):
//                 FD(p) - FD(first trunk pixel on p's path). Every sub-basin
//                 pixel reaches its own outlet, which is a trunk pixel, before
//                 it can reach any other sub-basin's trunk, so one mask holding
//                 every trunk serves all sub-basins. FD is a float raster, so
//                 this carries float rounding of FD (about 1e-7 relative).
//
// trunk_mask has full-grid flat size; nonzero marks a trunk pixel of some
// sub-basin. trunk_top_flat[s] is p_in.
// -----------------------------------------------------------------------------
struct subbasin_river_lengths_result
{
    std::vector<double>       rivlen_m;
    std::vector<double>       loc_rivlen_m;
    std::vector<std::int64_t> trunk_top_flat;
    std::vector<std::int64_t> n_trunk_pixels;
    std::vector<std::uint8_t> trunk_mask;
};

[[nodiscard]] subbasin_river_lengths_result compute_subbasin_river_lengths(
    const subbasin_partition_result& partition,
    const lsdtt_flow_base&           flow,
    std::int64_t                     channel_threshold_contributing_pixels
);

} // namespace lsdtt
