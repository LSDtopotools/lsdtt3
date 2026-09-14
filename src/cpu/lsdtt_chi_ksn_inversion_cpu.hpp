// lsdtt_chi_ksn_inversion_cpu.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// lsdtt_chi_ksn_inversion_cpu
//
// Spatial k_sn / u* inversion with m/n model selection.
//
// This module implements the spatial chi-based inversion in which, for each
// trial m/n value:
//
//   1. chi is recomputed on a fixed channel network
//   2. a sparse design matrix A is assembled from downstream chi increments
//   3. a coarse spatial field of u* (equivalent to k_sn when A0 = 1 m^2) is
//      solved by regularized least squares
//   4. predicted elevations and residuals are computed
//   5. model misfit is recorded
//
// The best-fit m/n is the trial value with the smallest inversion misfit.
//
// IMPORTANT:
//   - The inversion operates on a FIXED channel network geometry/topology.
//   - Chi is recalculated internally for each tested m/n.
//   - The first implementation is intended to use an unconstrained sparse
//     regularized solve. Bound constraints can be added later.
//   - The network input is generic: it only needs nodes, receiver topology,
//     and observed elevations. Helper builders are provided for common
//     lsdtt3 chi-network products.
// -----------------------------------------------------------------------------

#pragma once

#include <cstdint>
#include <limits>
#include <string>
#include <vector>

#include "../core/lsdtt_flow_base.hpp"
#include "lsdtt_chi_extraction_cpu.hpp"
#include "lsdtt_chi_network_analysis_cpu.hpp"
#include "lsdtt_raster_cpu.hpp"

namespace lsdtt
{

// -----------------------------------------------------------------------------
// Generic fixed channel-network input for k_sn / u* inversion.
//
// This stores the geometric/topological support of the inversion. Chi is NOT
// stored here because chi depends on m/n and is recomputed internally for each
// tested value.
//
// Requirements:
//   - flats, receiver_flats, and observed_elevation must have equal length
//   - nodes must be ordered receiver-before-donor
//   - flats[0] should normally be the outlet node
//   - receiver_flats[0] should normally be k_outlet_receiver
// -----------------------------------------------------------------------------
struct lsdtt_chi_ksn_inversion_network
{
    std::int64_t outlet_flat = k_nodata_receiver;

    std::vector<std::int64_t> flats;
    std::vector<std::int64_t> receiver_flats;
    std::vector<float>        observed_elevation;
};

// -----------------------------------------------------------------------------
// Per-trial summary.
//
// One of these is produced for each tested m/n value.
// -----------------------------------------------------------------------------
struct lsdtt_chi_ksn_trial_result
{
    double m_over_n = std::numeric_limits<double>::quiet_NaN();

    float mean_absolute_error = std::numeric_limits<float>::quiet_NaN();
    float root_mean_square_error = std::numeric_limits<float>::quiet_NaN();

    std::int32_t n_nodes = 0;
    std::int32_t n_parameters = 0;

    // Number of coarse-grid parameters that actually receive support from at
    // least one channel-node row in the design matrix.
    std::int32_t n_constrained_parameters = 0;
};

// -----------------------------------------------------------------------------
// Inversion options.
//
// coarse_cell_size:
//   Spatial size of the inversion grid in map units (normally metres).
//   This is larger than the DEM resolution; many channel nodes can map to one
//   inversion cell.
//
// roughness_weight:
//   Strength of the Laplacian roughness penalty on the coarse-grid u* field.
//
// diagonal_jitter:
//   Small positive value optionally added to the system diagonal for numerical
//   stabilization.
//
// A_0:
//   Reference area used for chi. If A_0 = 1 m^2 then u* is equivalent to k_sn.
//
// use_mean_absolute_error_for_model_selection:
//   If true, the best m/n is selected by mean absolute elevation misfit,
//   matching the first Matlab implementation you shared. RMSE is still
//   reported for diagnostics.
// -----------------------------------------------------------------------------
struct lsdtt_chi_ksn_inversion_options
{
    double A_0 = 1.0;

    double coarse_cell_size = 2000.0;
    double roughness_weight = 1.0e4;
    double diagonal_jitter  = 1.0e-8;

    bool use_mean_absolute_error_for_model_selection = true;
};

// -----------------------------------------------------------------------------
// Main inversion result for one basin/network.
//
// The "best" fields correspond to the trial m/n that minimizes the chosen
// model-selection misfit.
//
// u_star_field stores the coarse inversion-grid solution.
// node_u_star stores the coarse-grid value assigned back to each channel node.
// chi, predicted_elevation, and residual_elevation are all node-aligned with
// the input network.
// -----------------------------------------------------------------------------
struct lsdtt_chi_ksn_inversion_result
{
    lsdtt_chi_ksn_inversion_network network;

    std::vector<double> tested_m_over_n;
    std::vector<lsdtt_chi_ksn_trial_result> trials;

    std::int32_t best_trial_index = -1;
    double best_m_over_n = std::numeric_limits<double>::quiet_NaN();

    float best_mean_absolute_error = std::numeric_limits<float>::quiet_NaN();
    float best_root_mean_square_error = std::numeric_limits<float>::quiet_NaN();

    // Best-fit nodewise fields
    std::vector<float> chi;
    std::vector<float> predicted_elevation;
    std::vector<float> residual_elevation;
    std::vector<float> node_u_star;

    // Mapping from node i to coarse-grid parameter index.
    // -1 means the node could not be assigned, though this should normally not
    // occur for valid inputs.
    std::vector<std::int32_t> coarse_cell_id_of_node;

    // Best-fit coarse-grid outputs
    lsdtt_raster_cpu u_star_field;
    lsdtt_raster_cpu support_count_field;
};

// -----------------------------------------------------------------------------
// Optional collection wrapper for later basin-by-basin workflows.
// -----------------------------------------------------------------------------
struct lsdtt_chi_ksn_inversion_collection_result
{
    std::vector<lsdtt_chi_ksn_inversion_result> basins;
};

// =============================================================================
// Helper builders for the generic inversion-network input
// =============================================================================

// Build a fixed inversion network from an extracted chi basin and a DEM.
// Elevations are sampled from the DEM at the basin-network nodes.
[[nodiscard]] lsdtt_chi_ksn_inversion_network make_chi_ksn_inversion_network(
    const lsdtt_chi_basin_result& chi_basin,
    const lsdtt_raster_cpu& dem
);

// Build a fixed inversion network from a fitted chi network.
// This reuses the observed elevations already stored in the fit result.
[[nodiscard]] lsdtt_chi_ksn_inversion_network make_chi_ksn_inversion_network(
    const lsdtt_chi_network_fit_result& fit_result
);

// =============================================================================
// Spatial k_sn / u* inversion with m/n sweep
// =============================================================================

// Area-based chi inversion using drainage area derived from the flow object.
[[nodiscard]] lsdtt_chi_ksn_inversion_result invert_chi_ksn(
    const lsdtt_flow_base& flow,
    const lsdtt_chi_ksn_inversion_network& network,
    const std::vector<double>& tested_m_over_n,
    const lsdtt_chi_ksn_inversion_options& options = {}
);

// Custom-denominator chi inversion using a supplied positive raster
// (for example discharge or effective area) on the same grid as the flow.
[[nodiscard]] lsdtt_chi_ksn_inversion_result invert_chi_ksn(
    const lsdtt_flow_base& flow,
    const lsdtt_chi_ksn_inversion_network& network,
    const lsdtt_raster_cpu& denominator_raster,
    const std::vector<double>& tested_m_over_n,
    const lsdtt_chi_ksn_inversion_options& options = {}
);

// =============================================================================
// Collection overloads
// =============================================================================

[[nodiscard]] lsdtt_chi_ksn_inversion_collection_result invert_chi_ksn(
    const lsdtt_flow_base& flow,
    const std::vector<lsdtt_chi_ksn_inversion_network>& networks,
    const std::vector<double>& tested_m_over_n,
    const lsdtt_chi_ksn_inversion_options& options = {}
);

[[nodiscard]] lsdtt_chi_ksn_inversion_collection_result invert_chi_ksn(
    const lsdtt_flow_base& flow,
    const std::vector<lsdtt_chi_ksn_inversion_network>& networks,
    const lsdtt_raster_cpu& denominator_raster,
    const std::vector<double>& tested_m_over_n,
    const lsdtt_chi_ksn_inversion_options& options = {}
);

// =============================================================================
// CSV writers for inversion results
// =============================================================================

// Write a per-trial inversion summary (one row per basin x m/n).
void write_chi_ksn_inversion_trials_csv(
    const std::string& filename,
    const lsdtt_chi_ksn_inversion_collection_result& results
);

// Write a per-basin best-trial summary (one row per basin).
void write_chi_ksn_inversion_best_summary_csv(
    const std::string& filename,
    const lsdtt_chi_ksn_inversion_collection_result& results
);

} // namespace lsdtt
