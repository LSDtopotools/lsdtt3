// lsdtt_slope_stability_hazard_area_cpu.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// Coseismic landslide hazard area: the SHALRUN-EQ index of Milledge et al.
// (2019, NHESS 19, 837-856, eqs. 6-9).
//
// The hazard area Ah of a cell (m^2 per m of cell width) is the area of cells
// steep enough to initiate a landslide whose debris can reach it without
// crossing a cell too gentle for run-out to continue:
//
//   source    steepest downslope gradient >= initiation_slope;
//   routing   multiple flow directions, each receiver weighted by its gradient
//             (Quinn et al. 1991, exponent 1);
//   stopping  a cell with steepest gradient < stopping_slope is reached, and so
//             counted, but passes nothing on.
//
// The paper's defaults are 40 and 10 degrees, here the gradients 0.839 and
// 0.176. They are effective, resolution-dependent values, not material
// properties.
//
// Eq. 7 writes the connection as binary, but the paper reports Ah of about
// 0.1 m^2/m on a 30 m grid, which a binary count (at least one cell width)
// cannot produce. This routine therefore routes MFD FRACTIONS. Both readings
// give the same set of cells with Ah = 0. Any positive Ah too small for Float32
// is written as the smallest normal float, so Ah = 0 always means unreachable.
//
// Ah is a relative index of exposure. It takes no account of shaking, strength,
// landslide size or run-out over reverse slopes. A cell whose debris sources may
// lie beyond the DEM has an Ah that is a lower bound (see edge_influenced).
// -----------------------------------------------------------------------------

#pragma once

#include <cstdint>
#include <string>
#include <vector>

#include "lsdtt_raster_cpu.hpp"

namespace lsdtt
{

inline constexpr std::uint8_t k_hazard_zone_transit  = 0;  // neither source nor stopping
inline constexpr std::uint8_t k_hazard_zone_source   = 1;  // gradient >= initiation_slope
inline constexpr std::uint8_t k_hazard_zone_stopping = 2;  // gradient < stopping_slope

struct coseismic_hazard_area_options
{
    double initiation_slope = 0.839;  // gradient (m/m), tan 40 degrees
    double stopping_slope = 0.176;    // gradient (m/m), tan 10 degrees
    std::string boundary_conditions = "oooo";
    int n_threads = 0;
};

struct coseismic_hazard_area_result
{
    lsdtt_raster_cpu hazard_area;      // Ah (m^2/m)
    lsdtt_raster_cpu steepest_slope;   // steepest downslope gradient (m/m)
    lsdtt_raster_cpu zone;             // k_hazard_zone_*
    lsdtt_raster_cpu edge_influenced;  // 1 if debris from a DEM-edge cell reaches it

    std::int64_t n_valid = 0;
    std::int64_t n_source = 0;
    std::int64_t n_stopping = 0;
    std::int64_t n_zero_hazard = 0;
    std::int64_t n_edge_influenced = 0;
    double source_area_m2 = 0.0;
};

// Throws std::invalid_argument unless 0 < stopping_slope < initiation_slope.
void validate_coseismic_hazard_area_options(const coseismic_hazard_area_options& options);

// Eq. 9 on an explicit receiver graph (compressed sparse row). `order` lists
// donors before receivers. `passes[i]` is nonzero if cell i passes material on.
// `seed` is each cell's own contribution. If `reached` is given it is set to 1
// on every cell with a positive seed or a positive-weight arrival.
[[nodiscard]] std::vector<double> route_hazard_area(
    const std::vector<std::int64_t>& order,
    const std::vector<std::int64_t>& receiver_offsets,
    const std::vector<std::int64_t>& receiver_indices,
    const std::vector<float>& receiver_weights,
    const std::vector<std::uint8_t>& passes,
    const std::vector<double>& seed,
    std::vector<std::uint8_t>* reached = nullptr);

// `routing_surface` should be free of pits (carved and/or filled).
[[nodiscard]] coseismic_hazard_area_result compute_coseismic_hazard_area(
    const lsdtt_raster_cpu& routing_surface,
    const coseismic_hazard_area_options& options = {});

}  // namespace lsdtt
