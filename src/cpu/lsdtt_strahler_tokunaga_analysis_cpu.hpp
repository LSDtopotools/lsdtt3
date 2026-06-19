// -----------------------------------------------------------------------------
// lsdtt_strahler_tokunaga_analysis_cpu
//
// CPU-side analysis of Strahler channels and Tokunaga self-similarity.
//
// IMPORTANT TERMINOLOGY
// ---------------------
// The junction network stores a compact channel graph made of "compact links".
// A compact link runs from a source or upstream junction to the next downstream
// junction or outlet.
//
// This module defines a different derived object:
//
//   Strahler channel
//
// A Strahler channel is a maximal downstream chain of one or more compact links
// that all share the same Strahler order. A Strahler channel may therefore
// contain several compact links stitched together end-to-end, provided the
// downstream continuation remains in the same order.
//
// This distinction is important:
//
//   compact link      !=   Strahler channel
//
// The junction-network object remains the source of topology. This module only
// derives same-order aggregated channels and computes statistics from them.
//
// TOKUNAGA TERMINOLOGY
// --------------------
// For a Strahler channel of order i that flows into a receiving Strahler
// channel of order j:
//
//   Tokunaga pair   = (i, j)
//   Tokunaga offset = j - i
//
// Tokunaga self-similarity is usually analysed as a function of offset rather
// than by treating each order pair independently.
//
// OUTPUT PRODUCTS
// ---------------
// This module supports two different vector/tabular views of the same derived
// Strahler-channel geometry:
//
// 1. Strahler channel output
//    - one record per Strahler channel
//    - attributes emphasise Strahler metrics such as:
//         strahler_order
//         total_length_metres
//         elevation_drop_metres
//         mean_slope
//
// 2. Tokunaga channel output
//    - one record per Strahler channel
//    - geometry is the same as the Strahler output
//    - attributes emphasise Tokunaga topology such as:
//         donor_order
//         receiver_order
//         tokunaga_offset
//         receiver_channel_id
//
// In other words, the Strahler and Tokunaga vector outputs use the same channel
// line geometry but expose different thematic attribute tables.
//
// DESIGN NOTES
// ------------
// - All routines here are CPU-side free functions and plain structs.
// - These routines operate on:
//     * lsdtt_flow_base
//     * lsdtt_junction_network_base
//     * lsdtt_raster_cpu   (for elevation drop calculations)
// - CSV and vector output are supported through lsdtt_data_loader.
// - Vector output writes one line feature per Strahler channel.
// -----------------------------------------------------------------------------



#pragma once

#include <cstddef>
#include <cstdint>
#include <limits>
#include <string>
#include <vector>

#include "../core/lsdtt_data_loader.hpp"
#include "../core/lsdtt_flow_base.hpp"
#include "../core/lsdtt_junction_network_base.hpp"
#include "../core/lsdtt_types.hpp"
#include "lsdtt_raster_cpu.hpp"
#include "../utils/lsdtt_raster_geometry.hpp"

namespace lsdtt
{

// -----------------------------------------------------------------------------
// lsdtt_strahler_channel
//
// One aggregated Strahler channel.
//
// A Strahler channel is a maximal chain of one or more compact links that all
// have the same Strahler order. Geometry is stored as a full ordered pixel chain
// from the upstream endpoint to the downstream endpoint.
// -----------------------------------------------------------------------------
struct lsdtt_strahler_channel
{
    std::int32_t channel_id = -1;

    // Strahler order of this aggregated same-order channel.
    std::int32_t strahler_order = -1;

    // Endpoints of the aggregated channel.
    std::int32_t upstream_junction   = k_nodata_junction;
    std::int32_t downstream_junction = k_outlet_junction;

    std::int64_t upstream_flat   = -1;
    std::int64_t downstream_flat = -1;

    // Upstream compact-link junction ids merged into this Strahler channel.
    // Each entry identifies one compact link by its upstream junction.
    std::vector<std::int32_t> compact_link_upstream_junctions;

    // Full geometry from upstream to downstream.
    std::vector<std::int64_t> flats;
    std::vector<lsdtt_xy>     xy;
    std::vector<double>       cumulative_distance_from_upstream_metres;

    // Longitudinal metrics.
    double total_length_metres    = std::numeric_limits<double>::quiet_NaN();
    double elevation_drop_metres  = std::numeric_limits<double>::quiet_NaN();
    double mean_slope             = std::numeric_limits<double>::quiet_NaN();

    double upstream_elevation     = std::numeric_limits<double>::quiet_NaN();
    double downstream_elevation   = std::numeric_limits<double>::quiet_NaN();

    // Tokunaga receiver information.
    // receiver_channel_id < 0 means no receiving Strahler channel exists
    // (typically because this channel drains directly to the outlet).
    std::int32_t receiver_channel_id = -1;
    std::int32_t receiver_order      = -1;
    std::int32_t tokunaga_offset     = -1;

    bool flows_to_outlet = false;

    bool empty() const noexcept
    {
        return channel_id < 0 || flats.empty();
    }

    std::size_t size() const noexcept
    {
        return flats.size();
    }
};

// -----------------------------------------------------------------------------
// lsdtt_strahler_channel_network
//
// Container for all derived Strahler channels in a network.
// channel_id_of_upstream_junction maps a compact-link upstream junction id to
// the Strahler channel id it belongs to. Entries for invalid / unused junctions
// are set to -1.
// -----------------------------------------------------------------------------
struct lsdtt_strahler_channel_network
{
    std::vector<lsdtt_strahler_channel> channels;

    std::vector<std::int32_t> channel_id_of_upstream_junction;

    std::size_t size() const noexcept
    {
        return channels.size();
    }

    bool empty() const noexcept
    {
        return channels.empty();
    }
};

// -----------------------------------------------------------------------------
// lsdtt_tokunaga_pair_count
//
// Count of channels of order i flowing into channels of order j.
// -----------------------------------------------------------------------------
struct lsdtt_tokunaga_pair_count
{
    std::int32_t donor_order    = -1;
    std::int32_t receiver_order = -1;
    std::int64_t count          = 0;
};

// -----------------------------------------------------------------------------
// lsdtt_tokunaga_offset_stat
//
// Summary by Tokunaga offset k = j - i.
//
// count   : total count of all channels with this offset
// mean_t  : mean Tokunaga count for this offset after aggregating by receiver
//           order where appropriate in later analysis
// -----------------------------------------------------------------------------
struct lsdtt_tokunaga_offset_stat
{
    std::int32_t offset = -1;
    std::int64_t count  = 0;
    double mean_t       = std::numeric_limits<double>::quiet_NaN();
};

// -----------------------------------------------------------------------------
// lsdtt_tokunaga_fit
//
// Fit of the classical Tokunaga law:
//
//   T_k = a * c^(k - 1)
//
// The fitted exponent is log(c) in a log-linear regression of T_k against
// (k - 1).
// -----------------------------------------------------------------------------
struct lsdtt_tokunaga_fit
{
    bool valid = false;

    double a         = std::numeric_limits<double>::quiet_NaN();
    double c         = std::numeric_limits<double>::quiet_NaN();

    double intercept = std::numeric_limits<double>::quiet_NaN();
    double exponent  = std::numeric_limits<double>::quiet_NaN();
    double r_squared = std::numeric_limits<double>::quiet_NaN();

    std::size_t n_offsets_used = 0;
};

// -----------------------------------------------------------------------------
// lsdtt_tokunaga_self_similarity_test
//
// Full Tokunaga summary for a Strahler channel network.
// -----------------------------------------------------------------------------
struct lsdtt_tokunaga_self_similarity_test
{
    bool valid = false;

    std::vector<lsdtt_tokunaga_pair_count>  pair_counts;
    std::vector<lsdtt_tokunaga_offset_stat> offset_stats;

    lsdtt_tokunaga_fit fit;          // ordinary log-linear least squares
    lsdtt_tokunaga_fit fit_poisson;  // count-weighted Poisson GLM (IRLS)

    double rmse_log_space = std::numeric_limits<double>::quiet_NaN();
};

// -----------------------------------------------------------------------------
// build_strahler_channel_network
//
// Build maximal same-order Strahler channels by merging consecutive compact
// links of the same Strahler order. Also calculates:
//   - total channel length
//   - upstream and downstream elevation
//   - elevation drop
//   - mean slope
//
// Requires an elevation raster because Strahler-channel drop is defined from the
// upstream-most and downstream-most pixels of the aggregated channel.
// -----------------------------------------------------------------------------
lsdtt_strahler_channel_network build_strahler_channel_network(
    const lsdtt_flow_base&             flow,
    const lsdtt_junction_network_base& network,
    const lsdtt_raster_cpu&            elevation
);

// -----------------------------------------------------------------------------
// assign_tokunaga_receiver_channels
//
// Populate Tokunaga receiver information for each Strahler channel:
//
//   receiver_channel_id
//   receiver_order
//   tokunaga_offset = receiver_order - strahler_order
//
// Channels draining directly to the outlet retain receiver_channel_id < 0.
// -----------------------------------------------------------------------------
void assign_tokunaga_receiver_channels(
    lsdtt_strahler_channel_network& strahler_network
);

// -----------------------------------------------------------------------------
// calculate_tokunaga_pair_counts
//
// Count all observed Tokunaga order pairs (i, j), where i is the order of a
// Strahler channel and j is the order of the Strahler channel it flows into.
// Channels with no receiving Strahler channel are excluded.
// -----------------------------------------------------------------------------
std::vector<lsdtt_tokunaga_pair_count> calculate_tokunaga_pair_counts(
    const lsdtt_strahler_channel_network& strahler_network
);

// -----------------------------------------------------------------------------
// calculate_tokunaga_offset_statistics
//
// Group Tokunaga pair counts by offset k = j - i.
// -----------------------------------------------------------------------------
std::vector<lsdtt_tokunaga_offset_stat> calculate_tokunaga_offset_statistics(
    const std::vector<lsdtt_tokunaga_pair_count>& pair_counts
);

// -----------------------------------------------------------------------------
// fit_tokunaga_law
//
// Fit the classical Tokunaga law:
//
//   T_k = a * c^(k - 1)
//
// using the available offset statistics.
// -----------------------------------------------------------------------------
lsdtt_tokunaga_fit fit_tokunaga_law(
    const std::vector<lsdtt_tokunaga_offset_stat>& offset_stats
);

// -----------------------------------------------------------------------------
// fit_tokunaga_law_poisson
//
// Count-weighted Poisson GLM fit of T_k = a*c^(k-1) by iteratively reweighted
// least squares (log link). Each offset k contributes its observed mean T_k with
// a prior weight equal to its pair count, so sparsely-populated high offsets no
// longer dominate the slope as they do in ordinary log-linear least squares.
// This is the statistically appropriate model for T_k as a count rate.
// -----------------------------------------------------------------------------
lsdtt_tokunaga_fit fit_tokunaga_law_poisson(
    const std::vector<lsdtt_tokunaga_offset_stat>& offset_stats
);

// -----------------------------------------------------------------------------
// analyse_tokunaga_self_similarity
//
// One-stop Tokunaga analysis:
//
//   1. calculate pair counts
//   2. reduce to offset statistics
//   3. fit Tokunaga law
// -----------------------------------------------------------------------------
lsdtt_tokunaga_self_similarity_test analyse_tokunaga_self_similarity(
    const lsdtt_strahler_channel_network& strahler_network
);

// -----------------------------------------------------------------------------
// write_strahler_channels
//
// Write one record per Strahler channel.
//
// This output is the Strahler-focused thematic view of the derived channel
// network. Geometry is one line feature per Strahler channel. Attributes
// emphasise longitudinal and order metrics such as:
//
//   - strahler_order
//   - total channel length
//   - elevation drop
//   - mean slope
//
// Output format behaviour:
//   - csv         : tabular output with endpoint coordinates in EPSG:4326
//   - gpkg / fgb  : one line feature per Strahler channel in raster CRS
//
// NOTE:
// The line geometry written here is the same aggregated Strahler-channel
// geometry that is also used by write_tokunaga_channels(). The difference
// between the two writers is the attribute schema, not the geometry.
// -----------------------------------------------------------------------------

void write_strahler_channels(
    const lsdtt_strahler_channel_network& strahler_network,
    const lsdtt_flow_base&                flow,
    const std::string&                    out_path,
    point_output_format                   format,
    const std::string&                    layer_name = "strahler_channels",
    bool                                  spatial_index = true
);

// -----------------------------------------------------------------------------
// write_tokunaga_channel_data
//
// Write one CSV row per Strahler channel containing its Tokunaga receiver
// information and longitudinal metrics.
// -----------------------------------------------------------------------------
void write_tokunaga_channel_data(
    const lsdtt_strahler_channel_network& strahler_network,
    const std::string&                    out_path
);

// -----------------------------------------------------------------------------
// write_tokunaga_summary_data
//
// Write Tokunaga pair counts, offset statistics, and fit results to CSV.
//
// This function is intended for downstream plotting / statistical analysis.
// -----------------------------------------------------------------------------
void write_tokunaga_summary_data(
    const lsdtt_tokunaga_self_similarity_test& tokunaga_result,
    const std::string&                         out_path
);


// -----------------------------------------------------------------------------
// write_tokunaga_channels
//
// Write one record per Strahler channel, but with Tokunaga-focused attributes.
//
// This output uses the same line geometry as write_strahler_channels(): one
// line feature per derived Strahler channel. The difference is that the
// attribute schema is focused on Tokunaga topology rather than on Strahler
// longitudinal metrics.
//
// Typical Tokunaga-focused attributes are:
//
//   - donor_order        : order of the current Strahler channel
//   - receiver_order     : order of the Strahler channel it flows into
//   - tokunaga_offset    : receiver_order - donor_order
//   - receiver_channel_id
//   - flows_to_outlet
//
// Output format behaviour:
//   - csv         : tabular output with endpoint coordinates in EPSG:4326
//   - gpkg / fgb  : one line feature per Strahler channel in raster CRS
//
// IMPORTANT:
// This writer does NOT create a different geometry from the Strahler writer.
// It is a different thematic view of the same Strahler-channel network, with
// Tokunaga-oriented attributes.
// -----------------------------------------------------------------------------

void write_tokunaga_channels(
    const lsdtt_strahler_channel_network& strahler_network,
    const lsdtt_flow_base&                flow,
    const std::string&                    out_path,
    point_output_format                   format,
    const std::string&                    layer_name = "tokunaga_channels",
    bool                                  spatial_index = true
);


// =============================================================================
// Junction drainage-area ratios
//
// A continuous, high-sample-count complement to order-binned Tokunaga analysis.
// At every interior confluence two (or more) channels merge; the ratio of the
// two largest incoming drainage areas r = A_trib / A_main in (0,1] is recorded
// at EVERY junction, giving far more samples than the few well-populated
// Strahler-order bins. The distribution of r encodes the network's
// area-accumulation self-similarity.
//
// CAVEATS (see docs/synthetic_networks_plan.md and the project notes):
//   - This measures area accumulation, NOT side-branching topology; it does not
//     uniquely recover Tokunaga (a,c). The c_equivalent below is a
//     model-dependent estimate (c ~= R_A^lambda), not a substitute for (a,c).
//   - Junctions are nested and hence autocorrelated; the naive summary here does
//     not branch-bootstrap, so the spread is informative but the CIs are not.
//   - Sensitive to the channel-extraction threshold (truncates the small-r tail)
//     and to the drainage-area definition (here: D8 contributing pixels).
// =============================================================================

// One record per interior confluence junction.
struct lsdtt_junction_area_ratio
{
    std::int32_t junction_id     = -1;
    double       x               = std::numeric_limits<double>::quiet_NaN();
    double       y               = std::numeric_limits<double>::quiet_NaN();
    double       area_main_cells = 0.0;  // A_main (larger tributary), pixels
    double       area_trib_cells = 0.0;  // A_trib (smaller tributary), pixels
    double       ratio           = std::numeric_limits<double>::quiet_NaN(); // A_trib/A_main
    double       area_down_cells = 0.0;  // contributing pixels at the confluence
    std::int32_t order_main      = -1;   // Strahler order of larger tributary
    std::int32_t order_trib      = -1;   // Strahler order of smaller tributary
    std::int32_t order_receiver  = -1;   // Strahler order of the confluence
    std::int32_t n_tributaries   = 0;    // number of channel tributaries merging

    // Branch label for branch-bootstrap resampling: the non-mainstem tributary
    // subtree this confluence belongs to. -1 means the sample is on the main
    // stem (always retained during resampling).
    std::int32_t branch_id       = -1;
};

// One row of the summary scalars; reused for the point estimate and for each
// bootstrap replicate.
struct lsdtt_area_ratio_summary
{
    double median_ratio = std::numeric_limits<double>::quiet_NaN();
    double lambda_all   = std::numeric_limits<double>::quiet_NaN();
    double lambda_tail  = std::numeric_limits<double>::quiet_NaN();
    double r_a          = std::numeric_limits<double>::quiet_NaN();
    double c_equivalent = std::numeric_limits<double>::quiet_NaN();
};

// One row of the area-decade collapse table: summary statistics restricted to
// confluences whose downstream area falls in [area_low, area_high). If the
// network is self-similar these rows collapse (similar median_ratio / lambda).
struct lsdtt_area_ratio_decade
{
    double       area_low     = std::numeric_limits<double>::quiet_NaN();
    double       area_high    = std::numeric_limits<double>::quiet_NaN();
    std::size_t  n            = 0;
    double       median_ratio = std::numeric_limits<double>::quiet_NaN();
    double       lambda_tail  = std::numeric_limits<double>::quiet_NaN();
};

struct lsdtt_junction_area_ratio_result
{
    std::vector<lsdtt_junction_area_ratio> junctions;

    std::size_t n_junctions   = 0;
    double mean_ratio         = std::numeric_limits<double>::quiet_NaN();
    double median_ratio       = std::numeric_limits<double>::quiet_NaN();

    // Exponential rate of x = -ln(r): lambda_all over all junctions, lambda_tail
    // over the small-r tail (x above its median). A power-law tail in r is an
    // exponential tail in x.
    double lambda_all         = std::numeric_limits<double>::quiet_NaN();
    double lambda_tail        = std::numeric_limits<double>::quiet_NaN();

    // Horton area ratio R_A from a per-order regression of mean downstream area,
    // and the model-dependent Tokunaga-equivalent c = R_A^lambda_tail.
    double r_a                = std::numeric_limits<double>::quiet_NaN();
    double c_equivalent       = std::numeric_limits<double>::quiet_NaN();

    // Branch-bootstrap 95% confidence intervals (2.5 / 97.5 percentiles). NaN if
    // bootstrapping was not requested. The main stem is always retained and
    // resampling is over non-mainstem tributary branches, the correct unit of
    // independence for nested (autocorrelated) junctions.
    std::int32_t n_bootstrap        = 0;
    std::int32_t n_branches         = 0;
    double median_ratio_ci_low  = std::numeric_limits<double>::quiet_NaN();
    double median_ratio_ci_high = std::numeric_limits<double>::quiet_NaN();
    double lambda_tail_ci_low   = std::numeric_limits<double>::quiet_NaN();
    double lambda_tail_ci_high  = std::numeric_limits<double>::quiet_NaN();
    double r_a_ci_low           = std::numeric_limits<double>::quiet_NaN();
    double r_a_ci_high          = std::numeric_limits<double>::quiet_NaN();
    double c_equivalent_ci_low  = std::numeric_limits<double>::quiet_NaN();
    double c_equivalent_ci_high = std::numeric_limits<double>::quiet_NaN();

    // Area-decade collapse table (one row per populated log10 area decade).
    std::vector<lsdtt_area_ratio_decade> decades;
};

// -----------------------------------------------------------------------------
// compute_junction_area_ratios
//
// Enumerates interior confluence junctions of the network. For each, the merging
// tributaries are the D8 channel-donors of the confluence node; their D8
// contributing-pixel counts give A_main / A_trib and r. Orders come from the
// junction governing each donor cell. Returns per-junction records plus summary
// statistics, the area-decade collapse table, and (if n_bootstrap > 0)
// branch-bootstrap confidence intervals. Requires a junction network with
// Strahler orders.
// -----------------------------------------------------------------------------
lsdtt_junction_area_ratio_result compute_junction_area_ratios(
    const lsdtt_flow_base&             flow,
    const lsdtt_junction_network_base& network,
    std::int32_t                       n_bootstrap = 0,
    std::uint64_t                      random_seed = 5489u
);

// Per-junction CSV (one row per confluence; includes branch_id).
void write_junction_area_ratio_data(
    const lsdtt_junction_area_ratio_result& result,
    const std::string&                      out_path
);

// Single-row summary CSV (distribution statistics, tail exponents, and CIs).
void write_junction_area_ratio_summary(
    const lsdtt_junction_area_ratio_result& result,
    const std::string&                      out_path
);

// Area-decade collapse CSV (one row per log10 area decade).
void write_junction_area_ratio_decades(
    const lsdtt_junction_area_ratio_result& result,
    const std::string&                      out_path
);


} // namespace lsdtt
