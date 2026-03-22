// lsdtt_raster_flow_preprocessing_cpu.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2025 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html
//
// -----------------------------------------------------------------------------
// Flow preprocessing — CPU free functions that condition a DEM prior to flow
// routing. All functions are stateless: they take a raster and return a new
// raster. The input is never modified.
//
// Fill algorithm:
//   Barnes, R., Lehman, C., Mulla, D. (2014) Priority-flood: An optimal
//   depression-filling and watershed-labeling algorithm for digital elevation
//   models. Computers & Geosciences, 62, 117-127.
//   DOI: 10.1016/j.cageo.2013.04.024
//
// Breach algorithm:
//   Lindsay, J.B. (2016) Efficient hybrid breaching-filling sink removal
//   methods for flow path enforcement in digital elevation models.
//   Hydrological Processes, 30(6), 846-857.
//   DOI: https://doi.org/10.1002/hyp.10648
//   Adapted from RichDEM (r-barnes/richdem, MIT licence).
//
// Least-cost breach algorithm:
//   Lindsay, J.B. (2020) Pit-centric depression removal methods.
//   Geomorphometry 2020 Conference Proceedings.
//   Implementation reference: WhiteboxTools BreachDepressionsLeastCost.
// -----------------------------------------------------------------------------

#pragma once

#include <limits>
#include <string>
#include <vector>

#include "lsdtt_raster_cpu.hpp"

namespace lsdtt {

// -----------------------------------------------------------------------------
// mask_to_nodata_using_threshold
// -----------------------------------------------------------------------------
lsdtt_raster_cpu mask_to_nodata_using_threshold(
    const lsdtt_raster_cpu& raster,
    float threshold,
    bool  is_below_threshold_nodata
);

// -----------------------------------------------------------------------------
// fill_pits
// -----------------------------------------------------------------------------
lsdtt_raster_cpu fill_pits(
    const lsdtt_raster_cpu& raster,
    float min_slope = 0.0001f,
    int   n_threads = 0
);

// -----------------------------------------------------------------------------
// fill_pits_with_boundary
// -----------------------------------------------------------------------------
lsdtt_raster_cpu fill_pits_with_boundary(
    const lsdtt_raster_cpu&              raster,
    float                                min_slope,
    const std::vector<std::vector<int>>& boundary_rows,
    const std::vector<std::vector<int>>& boundary_cols,
    const std::vector<std::string>&      boundary_conditions,
    const std::string&                   outlet_label = "outlet",
    int                                  n_threads    = 0
);

// -----------------------------------------------------------------------------
// breach_depressions
// -----------------------------------------------------------------------------
lsdtt_raster_cpu breach_depressions(
    const lsdtt_raster_cpu& raster,
    int n_threads = 0
);

// -----------------------------------------------------------------------------
// breach_depressions_least_cost
// -----------------------------------------------------------------------------
lsdtt_raster_cpu breach_depressions_least_cost(
    const lsdtt_raster_cpu& raster,
    int   max_dist       = 1000,
    float max_cost       = std::numeric_limits<float>::max(),
    bool  min_dist       = true,
    float flat_increment = 0.0f,
    bool  fill_remainder = true,
    int   n_threads      = 0
);

} // namespace lsdtt
