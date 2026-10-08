// test_grid_neighbours.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
//
// -----------------------------------------------------------------------------
// Characterisation test for src/utils/lsdtt_grid_neighbours.hpp
// (consolidation C01 and C19).
//
// Neighbour tables: every module that had its own literal offset table now binds
// its old local name to a shared table. Each removed literal is kept below,
// copied from the source with its file and line, and checked at COMPILE time to
// equal the shared table it was mapped to, element by element and in order (the
// order decides tie-breaking, so a table with the right offsets in the wrong
// order must fail). A mapping error therefore stops the build.
//
// sqrt(2): k_sqrt2_d must be the same double as std::sqrt(2.0) and every double
// literal it replaces, and k_sqrt2_f the same float as std::sqrt(2.0f) and every
// float literal it replaces. std::sqrt is not constexpr, so those are runtime
// checks, by bit pattern.
//
// Usage: lsdtt-test-grid-neighbours
// Prints one line per check and a final "RESULT: PASS|FAIL".
// -----------------------------------------------------------------------------

#include <cmath>
#include <cstdint>
#include <cstring>
#include <iostream>
#include <string>

#include "../../src/utils/lsdtt_grid_neighbours.hpp"

using namespace lsdtt;

namespace
{

template <int N>
constexpr bool same(const int (&a)[N], const int (&b)[N])
{
    for (int i = 0; i < N; ++i)
        if (a[i] != b[i]) return false;
    return true;
}

// ---- the removed literal tables, 39 sites ----------------------------

// src/cpu/lsdtt_basin_outline_cpu.cpp:39-40 (as of 2026-10-07), stage 5a
//   static constexpr std::array<int, 8> d8_row{{-1, -1,  0,  1,  1,  1,  0, -1}};
//   static constexpr std::array<int, 8> d8_col{{ 0,  1,  1,  1,  0, -1, -1, -1}};
constexpr int site0_row[8] = {-1, -1, 0, 1, 1, 1, 0, -1};
constexpr int site0_col[8] = {0, 1, 1, 1, 0, -1, -1, -1};
static_assert(same(site0_row, grid::d8_clockwise_row) && same(site0_col, grid::d8_clockwise_col),
              "lsdtt_basin_outline_cpu.cpp:39 d8_row/d8_col is not grid::d8_clockwise");

// src/cpu/lsdtt_burn_vectormap_cpu.cpp:683-684 (as of 2026-10-07), stage 5a
//   static constexpr int DR[8] = {-1, -1, 0, 1, 1, 1, 0, -1};
//   static constexpr int DC[8] = { 0,  1, 1, 1, 0, -1, -1, -1};
constexpr int site1_row[8] = {-1, -1, 0, 1, 1, 1, 0, -1};
constexpr int site1_col[8] = {0, 1, 1, 1, 0, -1, -1, -1};
static_assert(same(site1_row, grid::d8_clockwise_row) && same(site1_col, grid::d8_clockwise_col),
              "lsdtt_burn_vectormap_cpu.cpp:683 DR/DC is not grid::d8_clockwise");

// src/cpu/lsdtt_cross_sections_cpu.cpp:32-33 (as of 2026-10-07), stage 5a
//   constexpr std::array<int, 8> d8_row{{-1, -1,  0,  1,  1,  1,  0, -1}};
//   constexpr std::array<int, 8> d8_col{{ 0,  1,  1,  1,  0, -1, -1, -1}};
constexpr int site2_row[8] = {-1, -1, 0, 1, 1, 1, 0, -1};
constexpr int site2_col[8] = {0, 1, 1, 1, 0, -1, -1, -1};
static_assert(same(site2_row, grid::d8_clockwise_row) && same(site2_col, grid::d8_clockwise_col),
              "lsdtt_cross_sections_cpu.cpp:32 d8_row/d8_col is not grid::d8_clockwise");

// src/cpu/lsdtt_fill_spill_merge_cpu.cpp:64-65 (as of 2026-10-07), stage 5a
//   constexpr int k_dy[8] = { -1, -1,  0,  1,  1,  1,  0, -1 };
//   constexpr int k_dx[8] = {  0,  1,  1,  1,  0, -1, -1, -1 };
constexpr int site3_row[8] = {-1, -1, 0, 1, 1, 1, 0, -1};
constexpr int site3_col[8] = {0, 1, 1, 1, 0, -1, -1, -1};
static_assert(same(site3_row, grid::d8_clockwise_row) && same(site3_col, grid::d8_clockwise_col),
              "lsdtt_fill_spill_merge_cpu.cpp:64 k_dy/k_dx is not grid::d8_clockwise");

// src/cpu/lsdtt_internal_drainage_cpu.cpp:1646-1647 (as of 2026-10-07), stage 5a
//   static constexpr int d_row[8] = { -1, -1, 0, 1, 1, 1, 0, -1 };
//   static constexpr int d_col[8] = { 0, 1, 1, 1, 0, -1, -1, -1 };
constexpr int site4_row[8] = {-1, -1, 0, 1, 1, 1, 0, -1};
constexpr int site4_col[8] = {0, 1, 1, 1, 0, -1, -1, -1};
static_assert(same(site4_row, grid::d8_clockwise_row) && same(site4_col, grid::d8_clockwise_col),
              "lsdtt_internal_drainage_cpu.cpp:1646 d_row/d_col is not grid::d8_clockwise");

// src/cpu/lsdtt_junction_network_cpu.cpp:44-45 (as of 2026-10-07), stage 5a
//   constexpr std::array<int, 8> d8_row{{-1, -1,  0,  1,  1,  1,  0, -1}};
//   constexpr std::array<int, 8> d8_col{{ 0,  1,  1,  1,  0, -1, -1, -1}};
constexpr int site5_row[8] = {-1, -1, 0, 1, 1, 1, 0, -1};
constexpr int site5_col[8] = {0, 1, 1, 1, 0, -1, -1, -1};
static_assert(same(site5_row, grid::d8_clockwise_row) && same(site5_col, grid::d8_clockwise_col),
              "lsdtt_junction_network_cpu.cpp:44 d8_row/d8_col is not grid::d8_clockwise");

// src/cpu/lsdtt_least_cost_path.cpp:305-306 (as of 2026-10-07), stage 5a
//   static constexpr int DR[8] = {-1, -1, 0, 1, 1, 1, 0, -1};
//   static constexpr int DC[8] = { 0,  1, 1, 1, 0, -1, -1, -1};
constexpr int site6_row[8] = {-1, -1, 0, 1, 1, 1, 0, -1};
constexpr int site6_col[8] = {0, 1, 1, 1, 0, -1, -1, -1};
static_assert(same(site6_row, grid::d8_clockwise_row) && same(site6_col, grid::d8_clockwise_col),
              "lsdtt_least_cost_path.cpp:305 DR/DC is not grid::d8_clockwise");

// src/cpu/lsdtt_proximity_cpu.cpp:779-780 (as of 2026-10-07), stage 5a
//   static constexpr std::array<int,8> k_d8_row_prox{{-1,-1,0,1,1,1,0,-1}};
//   static constexpr std::array<int,8> k_d8_col_prox{{0,1,1,1,0,-1,-1,-1}};
constexpr int site7_row[8] = {-1, -1, 0, 1, 1, 1, 0, -1};
constexpr int site7_col[8] = {0, 1, 1, 1, 0, -1, -1, -1};
static_assert(same(site7_row, grid::d8_clockwise_row) && same(site7_col, grid::d8_clockwise_col),
              "lsdtt_proximity_cpu.cpp:779 k_d8_row_prox/k_d8_col_prox is not grid::d8_clockwise");

// src/cpu/lsdtt_raster_flow_preprocessing_cpu.cpp:252-253 (as of 2026-10-07), stage 5a
//   static const int k_dr[8] = { -1, -1,  0,  1,  1,  1,  0, -1 };
//   static const int k_dc[8] = {  0,  1,  1,  1,  0, -1, -1, -1 };
constexpr int site8_row[8] = {-1, -1, 0, 1, 1, 1, 0, -1};
constexpr int site8_col[8] = {0, 1, 1, 1, 0, -1, -1, -1};
static_assert(same(site8_row, grid::d8_clockwise_row) && same(site8_col, grid::d8_clockwise_col),
              "lsdtt_raster_flow_preprocessing_cpu.cpp:252 k_dr/k_dc is not grid::d8_clockwise");

// src/cpu/lsdtt_raster_nodata_filling_cpu.cpp:31-32 (as of 2026-10-07), stage 5a
//   static const int k_nd_dr[8] = { -1, -1,  0,  1,  1,  1,  0, -1 };
//   static const int k_nd_dc[8] = {  0,  1,  1,  1,  0, -1, -1, -1 };
constexpr int site9_row[8] = {-1, -1, 0, 1, 1, 1, 0, -1};
constexpr int site9_col[8] = {0, 1, 1, 1, 0, -1, -1, -1};
static_assert(same(site9_row, grid::d8_clockwise_row) && same(site9_col, grid::d8_clockwise_col),
              "lsdtt_raster_nodata_filling_cpu.cpp:31 k_nd_dr/k_nd_dc is not grid::d8_clockwise");

// src/cpu/lsdtt_raster_skeleton_graph_cpu.cpp:787-788 (as of 2026-10-07), stage 5a
//   constexpr int DR8[8] = {-1,-1, 0, 1, 1, 1, 0,-1};
//   constexpr int DC8[8] = { 0, 1, 1, 1, 0,-1,-1,-1};
constexpr int site10_row[8] = {-1, -1, 0, 1, 1, 1, 0, -1};
constexpr int site10_col[8] = {0, 1, 1, 1, 0, -1, -1, -1};
static_assert(same(site10_row, grid::d8_clockwise_row) && same(site10_col, grid::d8_clockwise_col),
              "lsdtt_raster_skeleton_graph_cpu.cpp:787 DR8/DC8 is not grid::d8_clockwise");

// src/cpu/lsdtt_subbasin_partition_cpu.cpp:94-95 (as of 2026-10-07), stage 5a
//   static constexpr int d_row[8] = {-1, -1, 0, 1, 1, 1, 0, -1};
//   static constexpr int d_col[8] = {0, 1, 1, 1, 0, -1, -1, -1};
constexpr int site11_row[8] = {-1, -1, 0, 1, 1, 1, 0, -1};
constexpr int site11_col[8] = {0, 1, 1, 1, 0, -1, -1, -1};
static_assert(same(site11_row, grid::d8_clockwise_row) && same(site11_col, grid::d8_clockwise_col),
              "lsdtt_subbasin_partition_cpu.cpp:94 d_row/d_col is not grid::d8_clockwise");

// src/cpu/lsdtt_landscape_drape_cpu.cpp:74-75 (as of 2026-10-07), stage 5b
//   static constexpr int k_drow[8] = {-1, -1, -1,  0,  0,  1,  1,  1};
//   static constexpr int k_dcol[8] = {-1,  0,  1, -1,  1, -1,  0,  1};
constexpr int site12_row[8] = {-1, -1, -1, 0, 0, 1, 1, 1};
constexpr int site12_col[8] = {-1, 0, 1, -1, 1, -1, 0, 1};
static_assert(same(site12_row, grid::d8_rowmajor_row) && same(site12_col, grid::d8_rowmajor_col),
              "lsdtt_landscape_drape_cpu.cpp:74 k_drow/k_dcol is not grid::d8_rowmajor");

// src/cpu/lsdtt_least_cost_path.cpp:121-122 (as of 2026-10-07), stage 5b
//   static constexpr int DR[8] = {-1,-1,-1, 0, 0, 1, 1, 1};
//   static constexpr int DC[8] = {-1, 0, 1,-1, 1,-1, 0, 1};
constexpr int site13_row[8] = {-1, -1, -1, 0, 0, 1, 1, 1};
constexpr int site13_col[8] = {-1, 0, 1, -1, 1, -1, 0, 1};
static_assert(same(site13_row, grid::d8_rowmajor_row) && same(site13_col, grid::d8_rowmajor_col),
              "lsdtt_least_cost_path.cpp:121 DR/DC is not grid::d8_rowmajor");

// src/cpu/lsdtt_mountain_front_cpu.cpp:421-422 (as of 2026-10-07), stage 5b
//   static const int drow[8] = {-1,-1,-1, 0, 0, 1, 1, 1};
//   static const int dcol[8] = {-1, 0, 1,-1, 1,-1, 0, 1};
constexpr int site14_row[8] = {-1, -1, -1, 0, 0, 1, 1, 1};
constexpr int site14_col[8] = {-1, 0, 1, -1, 1, -1, 0, 1};
static_assert(same(site14_row, grid::d8_rowmajor_row) && same(site14_col, grid::d8_rowmajor_col),
              "lsdtt_mountain_front_cpu.cpp:421 drow/dcol is not grid::d8_rowmajor");

// src/cpu/lsdtt_proximity_cpu.cpp:447-448 (as of 2026-10-07), stage 5b
//   static const int DR[8] = {-1,-1,-1, 0, 0, 1, 1, 1};
//   static const int DC[8] = {-1, 0, 1,-1, 1,-1, 0, 1};
constexpr int site15_row[8] = {-1, -1, -1, 0, 0, 1, 1, 1};
constexpr int site15_col[8] = {-1, 0, 1, -1, 1, -1, 0, 1};
static_assert(same(site15_row, grid::d8_rowmajor_row) && same(site15_col, grid::d8_rowmajor_col),
              "lsdtt_proximity_cpu.cpp:447 DR/DC is not grid::d8_rowmajor");

// src/cpu/lsdtt_raster_binary_mask_cpu.cpp:697-698 (as of 2026-10-07), stage 5b
//   constexpr int DR8[8] = {-1,-1,-1, 0, 0, 1, 1, 1};
//   constexpr int DC8[8] = {-1, 0, 1,-1, 1,-1, 0, 1};
constexpr int site16_row[8] = {-1, -1, -1, 0, 0, 1, 1, 1};
constexpr int site16_col[8] = {-1, 0, 1, -1, 1, -1, 0, 1};
static_assert(same(site16_row, grid::d8_rowmajor_row) && same(site16_col, grid::d8_rowmajor_col),
              "lsdtt_raster_binary_mask_cpu.cpp:697 DR8/DC8 is not grid::d8_rowmajor");

// src/cpu/lsdtt_raster_binary_mask_cpu.cpp:846-847 (as of 2026-10-07), stage 5b
//   constexpr int DR[4] = {-1, 1,  0, 0};
//   constexpr int DC[4] = { 0, 0, -1, 1};
constexpr int site17_row[4] = {-1, 1, 0, 0};
constexpr int site17_col[4] = {0, 0, -1, 1};
static_assert(same(site17_row, grid::d4_nswe_row) && same(site17_col, grid::d4_nswe_col),
              "lsdtt_raster_binary_mask_cpu.cpp:846 DR/DC is not grid::d4_nswe");

// src/cpu/lsdtt_raster_binary_mask_cpu.cpp:1244-1245 (as of 2026-10-07), stage 5b
//   constexpr int k_dr[8] = { -1, -1, -1, 0, 0, 1, 1, 1 };
//   constexpr int k_dc[8] = { -1,  0,  1, -1, 1, -1, 0, 1 };
constexpr int site18_row[8] = {-1, -1, -1, 0, 0, 1, 1, 1};
constexpr int site18_col[8] = {-1, 0, 1, -1, 1, -1, 0, 1};
static_assert(same(site18_row, grid::d8_rowmajor_row) && same(site18_col, grid::d8_rowmajor_col),
              "lsdtt_raster_binary_mask_cpu.cpp:1244 k_dr/k_dc is not grid::d8_rowmajor");

// src/cpu/lsdtt_raster_binary_mask_vectorise_cpu.cpp:142-143 (as of 2026-10-07), stage 5b
//   const int dr4[4] = {-1, 1,  0, 0};
//   const int dc4[4] = { 0, 0, -1, 1};
constexpr int site19_row[4] = {-1, 1, 0, 0};
constexpr int site19_col[4] = {0, 0, -1, 1};
static_assert(same(site19_row, grid::d4_nswe_row) && same(site19_col, grid::d4_nswe_col),
              "lsdtt_raster_binary_mask_vectorise_cpu.cpp:142 dr4/dc4 is not grid::d4_nswe");

// src/cpu/lsdtt_raster_feature_finder_cpu.cpp:2103-2104 (as of 2026-10-07), stage 5b
//   static const int DR[8] = {-1,-1,-1, 0, 0, 1, 1, 1};
//   static const int DC[8] = {-1, 0, 1,-1, 1,-1, 0, 1};
constexpr int site20_row[8] = {-1, -1, -1, 0, 0, 1, 1, 1};
constexpr int site20_col[8] = {-1, 0, 1, -1, 1, -1, 0, 1};
static_assert(same(site20_row, grid::d8_rowmajor_row) && same(site20_col, grid::d8_rowmajor_col),
              "lsdtt_raster_feature_finder_cpu.cpp:2103 DR/DC is not grid::d8_rowmajor");

// src/cpu/lsdtt_raster_flow_preprocessing_cpu.cpp:2312-2313 (as of 2026-10-07), stage 5b
//   const int dr[4] = {-1, 1, 0, 0};
//   const int dc[4] = { 0, 0,-1, 1};
constexpr int site21_row[4] = {-1, 1, 0, 0};
constexpr int site21_col[4] = {0, 0, -1, 1};
static_assert(same(site21_row, grid::d4_nswe_row) && same(site21_col, grid::d4_nswe_col),
              "lsdtt_raster_flow_preprocessing_cpu.cpp:2312 dr/dc is not grid::d4_nswe");

// src/cpu/lsdtt_raster_flow_preprocessing_cpu.cpp:2443-2444 (as of 2026-10-07), stage 5b
//   const int dr[4] = {-1, 1, 0, 0};
//   const int dc[4] = { 0, 0,-1, 1};
constexpr int site22_row[4] = {-1, 1, 0, 0};
constexpr int site22_col[4] = {0, 0, -1, 1};
static_assert(same(site22_row, grid::d4_nswe_row) && same(site22_col, grid::d4_nswe_col),
              "lsdtt_raster_flow_preprocessing_cpu.cpp:2443 dr/dc is not grid::d4_nswe");

// src/cpu/lsdtt_raster_skeleton_graph_cpu.cpp:1306-1307 (as of 2026-10-07), stage 5b
//   constexpr int DR8[8] = {-1,-1,-1, 0, 0, 1, 1, 1};
//   constexpr int DC8[8] = {-1, 0, 1,-1, 1,-1, 0, 1};
constexpr int site23_row[8] = {-1, -1, -1, 0, 0, 1, 1, 1};
constexpr int site23_col[8] = {-1, 0, 1, -1, 1, -1, 0, 1};
static_assert(same(site23_row, grid::d8_rowmajor_row) && same(site23_col, grid::d8_rowmajor_col),
              "lsdtt_raster_skeleton_graph_cpu.cpp:1306 DR8/DC8 is not grid::d8_rowmajor");

// src/cpu/lsdtt_raster_skeleton_graph_cpu.cpp:2108-2109 (as of 2026-10-07), stage 5b
//   constexpr int DR8[8] = {-1,-1,-1, 0, 0, 1, 1, 1};
//   constexpr int DC8[8] = {-1, 0, 1,-1, 1,-1, 0, 1};
constexpr int site24_row[8] = {-1, -1, -1, 0, 0, 1, 1, 1};
constexpr int site24_col[8] = {-1, 0, 1, -1, 1, -1, 0, 1};
static_assert(same(site24_row, grid::d8_rowmajor_row) && same(site24_col, grid::d8_rowmajor_col),
              "lsdtt_raster_skeleton_graph_cpu.cpp:2108 DR8/DC8 is not grid::d8_rowmajor");

// src/cpu/lsdtt_raster_skeleton_graph_cpu.cpp:3577-3578 (as of 2026-10-07), stage 5b
//   constexpr int DR8[8] = {-1,-1,-1, 0, 0, 1, 1, 1};
//   constexpr int DC8[8] = {-1, 0, 1,-1, 1,-1, 0, 1};
constexpr int site25_row[8] = {-1, -1, -1, 0, 0, 1, 1, 1};
constexpr int site25_col[8] = {-1, 0, 1, -1, 1, -1, 0, 1};
static_assert(same(site25_row, grid::d8_rowmajor_row) && same(site25_col, grid::d8_rowmajor_col),
              "lsdtt_raster_skeleton_graph_cpu.cpp:3577 DR8/DC8 is not grid::d8_rowmajor");

// src/cpu/lsdtt_reference_field_ladder_cpu.cpp:1398-1399 (as of 2026-10-07), stage 5b
//   const int  dr8[8] = {-1, -1, -1,  0, 0,  1, 1, 1};
//   const int  dc8[8] = {-1,  0,  1, -1, 1, -1, 0, 1};
constexpr int site26_row[8] = {-1, -1, -1, 0, 0, 1, 1, 1};
constexpr int site26_col[8] = {-1, 0, 1, -1, 1, -1, 0, 1};
static_assert(same(site26_row, grid::d8_rowmajor_row) && same(site26_col, grid::d8_rowmajor_col),
              "lsdtt_reference_field_ladder_cpu.cpp:1398 dr8/dc8 is not grid::d8_rowmajor");

// src/cpu/lsdtt_reference_field_ladder_cpu.cpp:1400-1401 (as of 2026-10-07), stage 5b
//   const int  dr4[4] = {-1,  0,  0,  1};
//   const int  dc4[4] = { 0, -1,  1,  0};
constexpr int site27_row[4] = {-1, 0, 0, 1};
constexpr int site27_col[4] = {0, -1, 1, 0};
static_assert(same(site27_row, grid::d4_nwes_row) && same(site27_col, grid::d4_nwes_col),
              "lsdtt_reference_field_ladder_cpu.cpp:1400 dr4/dc4 is not grid::d4_nwes");

// src/cpu/lsdtt_strahler_tokunaga_analysis_cpu.cpp:1402-1403 (as of 2026-10-07), stage 5b
//   static const int dR[8] = {-1, -1, -1,  0, 0,  1, 1, 1};
//   static const int dC[8] = {-1,  0,  1, -1, 1, -1, 0, 1};
constexpr int site28_row[8] = {-1, -1, -1, 0, 0, 1, 1, 1};
constexpr int site28_col[8] = {-1, 0, 1, -1, 1, -1, 0, 1};
static_assert(same(site28_row, grid::d8_rowmajor_row) && same(site28_col, grid::d8_rowmajor_col),
              "lsdtt_strahler_tokunaga_analysis_cpu.cpp:1402 dR/dC is not grid::d8_rowmajor");

// src/cpu/lsdtt_terrace_edge_tracer_cpu.cpp:72-73 (as of 2026-10-07), stage 5b
//   constexpr int DR8[8] = {-1,-1,-1, 0, 0, 1, 1, 1};
//   constexpr int DC8[8] = {-1, 0, 1,-1, 1,-1, 0, 1};
constexpr int site29_row[8] = {-1, -1, -1, 0, 0, 1, 1, 1};
constexpr int site29_col[8] = {-1, 0, 1, -1, 1, -1, 0, 1};
static_assert(same(site29_row, grid::d8_rowmajor_row) && same(site29_col, grid::d8_rowmajor_col),
              "lsdtt_terrace_edge_tracer_cpu.cpp:72 DR8/DC8 is not grid::d8_rowmajor");

// src/cpu/lsdtt_terrace_edge_tracer_cpu.cpp:328-329 (as of 2026-10-07), stage 5b
//   constexpr int DR8[8] = {-1,-1,-1, 0, 0, 1, 1, 1};
//   constexpr int DC8[8] = {-1, 0, 1,-1, 1,-1, 0, 1};
constexpr int site30_row[8] = {-1, -1, -1, 0, 0, 1, 1, 1};
constexpr int site30_col[8] = {-1, 0, 1, -1, 1, -1, 0, 1};
static_assert(same(site30_row, grid::d8_rowmajor_row) && same(site30_col, grid::d8_rowmajor_col),
              "lsdtt_terrace_edge_tracer_cpu.cpp:328 DR8/DC8 is not grid::d8_rowmajor");

// src/cpu/lsdtt_synthetic_domain_cpu.cpp:66-67 (as of 2026-10-07), stage 5b
//   const int dr[4] = {-1, 1, 0, 0};
//   const int dc[4] = {0, 0, -1, 1};
constexpr int site31_row[4] = {-1, 1, 0, 0};
constexpr int site31_col[4] = {0, 0, -1, 1};
static_assert(same(site31_row, grid::d4_nswe_row) && same(site31_col, grid::d4_nswe_col),
              "lsdtt_synthetic_domain_cpu.cpp:66 dr/dc is not grid::d4_nswe");

// src/cpu/lsdtt_synthetic_network_cpu.cpp:27-28 (as of 2026-10-07), stage 5b
//   static const int k_dr[8] = {-1, -1, -1,  0, 0,  1, 1, 1};
//   static const int k_dc[8] = {-1,  0,  1, -1, 1, -1, 0, 1};
constexpr int site32_row[8] = {-1, -1, -1, 0, 0, 1, 1, 1};
constexpr int site32_col[8] = {-1, 0, 1, -1, 1, -1, 0, 1};
static_assert(same(site32_row, grid::d8_rowmajor_row) && same(site32_col, grid::d8_rowmajor_col),
              "lsdtt_synthetic_network_cpu.cpp:27 k_dr/k_dc is not grid::d8_rowmajor");

// src/cpu/lsdtt_valley_centreline_cpu.cpp:44-45 (as of 2026-10-07), stage 5b
//   constexpr int k_dr[8] = { -1, -1, -1, 0, 0, 1, 1, 1 };
//   constexpr int k_dc[8] = { -1,  0,  1, -1, 1, -1, 0, 1 };
constexpr int site33_row[8] = {-1, -1, -1, 0, 0, 1, 1, 1};
constexpr int site33_col[8] = {-1, 0, 1, -1, 1, -1, 0, 1};
static_assert(same(site33_row, grid::d8_rowmajor_row) && same(site33_col, grid::d8_rowmajor_col),
              "lsdtt_valley_centreline_cpu.cpp:44 k_dr/k_dc is not grid::d8_rowmajor");

// src/cpu/lsdtt_flowinfo_cpu.cpp:48-49 (as of 2026-10-07), stage 5c
//   constexpr std::array<int, 8> d8_row{{-1, -1,  0,  1,  1,  1,  0, -1}};
//   constexpr std::array<int, 8> d8_col{{ 0,  1,  1,  1,  0, -1, -1, -1}};
constexpr int site34_row[8] = {-1, -1, 0, 1, 1, 1, 0, -1};
constexpr int site34_col[8] = {0, 1, 1, 1, 0, -1, -1, -1};
static_assert(same(site34_row, grid::d8_clockwise_row) && same(site34_col, grid::d8_clockwise_col),
              "lsdtt_flowinfo_cpu.cpp:48 d8_row/d8_col is not grid::d8_clockwise");

// src/cpu/lsdtt_depression_hierarchy_cpu.cpp:52-53 (as of 2026-10-07), stage 5c
//   constexpr std::array<int, 8> k_dy{{ -1, -1,  0,  1,  1,  1,  0, -1 }};
//   constexpr std::array<int, 8> k_dx{{  0,  1,  1,  1,  0, -1, -1, -1 }};
constexpr int site35_row[8] = {-1, -1, 0, 1, 1, 1, 0, -1};
constexpr int site35_col[8] = {0, 1, 1, 1, 0, -1, -1, -1};
static_assert(same(site35_row, grid::d8_clockwise_row) && same(site35_col, grid::d8_clockwise_col),
              "lsdtt_depression_hierarchy_cpu.cpp:52 k_dy/k_dx is not grid::d8_clockwise");

// src/cpu/lsdtt_graphflood_optimised_cpu.cpp:30-31 (as of 2026-10-07), stage 5c
//   constexpr int k_dr[8] = {-1, 0, 0, 1, -1, -1, 1, 1};
//   constexpr int k_dc[8] = {0, -1, 1, 0, -1, 1, -1, 1};
constexpr int site36_row[8] = {-1, 0, 0, 1, -1, -1, 1, 1};
constexpr int site36_col[8] = {0, -1, 1, 0, -1, 1, -1, 1};
static_assert(same(site36_row, grid::d8_cardinal_first_row) && same(site36_col, grid::d8_cardinal_first_col),
              "lsdtt_graphflood_optimised_cpu.cpp:30 k_dr/k_dc is not grid::d8_cardinal_first");

// src/cpu/lsdtt_graphflood_optimised_cpu.cpp:32-33 (as of 2026-10-07), stage 5c
//   constexpr int k_pf_dr[8] = {-1, -1, 0, 1, 1, 1, 0, -1};
//   constexpr int k_pf_dc[8] = {0, 1, 1, 1, 0, -1, -1, -1};
constexpr int site37_row[8] = {-1, -1, 0, 1, 1, 1, 0, -1};
constexpr int site37_col[8] = {0, 1, 1, 1, 0, -1, -1, -1};
static_assert(same(site37_row, grid::d8_clockwise_row) && same(site37_col, grid::d8_clockwise_col),
              "lsdtt_graphflood_optimised_cpu.cpp:32 k_pf_dr/k_pf_dc is not grid::d8_clockwise");

// src/cpu/lsdtt_flowgraph_cpu.cpp:70-71 (as of 2026-10-07), stage 5c
//   constexpr int k_nb_dr[8] = { -1,  0,  0,  1, -1, -1,  1,  1 };
//   constexpr int k_nb_dc[8] = {  0, -1,  1,  0, -1,  1, -1,  1 };
constexpr int site38_row[8] = {-1, 0, 0, 1, -1, -1, 1, 1};
constexpr int site38_col[8] = {0, -1, 1, 0, -1, 1, -1, 1};
static_assert(same(site38_row, grid::d8_cardinal_first_row) && same(site38_col, grid::d8_cardinal_first_col),
              "lsdtt_flowgraph_cpu.cpp:70 k_nb_dr/k_nb_dc is not grid::d8_cardinal_first");



// ---- stage 10 (C01-H, Simon's D6 approval): hillslope-solver face tables ----

// src/cpu/lsdtt_hillslope_newton_cpu.cpp:139-140 (as of 2026-10-08), stage 10
constexpr int h0_row[4] = { -1,  0, 1,  0 };
constexpr int h0_col[4] = {  0, 1,  0, -1 };
static_assert(same(h0_row, grid::d4_nesw_row) && same(h0_col, grid::d4_nesw_col),
              "lsdtt_hillslope_newton_cpu.cpp:139 k_dr/k_dc is not grid::d4_nesw");

// src/cpu/lsdtt_hillslope_newton_cpu.cpp:418-419 (as of 2026-10-08), stage 10
constexpr int h1_row[4] = { -1,  0, 1,  0 };
constexpr int h1_col[4] = {  0, 1,  0, -1 };
static_assert(same(h1_row, grid::d4_nesw_row) && same(h1_col, grid::d4_nesw_col),
              "lsdtt_hillslope_newton_cpu.cpp:418 k_dr/k_dc is not grid::d4_nesw");

// src/cpu/lsdtt_hillslope_cellular_cpu.cpp:453-454 (as of 2026-10-08), stage 10
constexpr int h2_row[4] = {-1, 0, 1, 0};
constexpr int h2_col[4] = {0, 1, 0, -1};
static_assert(same(h2_row, grid::d4_nesw_row) && same(h2_col, grid::d4_nesw_col),
              "lsdtt_hillslope_cellular_cpu.cpp:453 k_dr/k_dc is not grid::d4_nesw");

// src/cpu/lsdtt_hillslope_cellular_cpu.cpp:661-662 (as of 2026-10-08), stage 10
constexpr int h3_row[4] = {-1, 0, 1, 0};
constexpr int h3_col[4] = {0, 1, 0, -1};
static_assert(same(h3_row, grid::d4_nesw_row) && same(h3_col, grid::d4_nesw_col),
              "lsdtt_hillslope_cellular_cpu.cpp:661 k_dr/k_dc is not grid::d4_nesw");

// src/cpu/lsdtt_hillslope_operator_cpu.hpp:291-292 (as of 2026-10-08), stage 10
constexpr int h4_row[4] = { -1,  0, 1,  0 };
constexpr int h4_col[4] = {  0, 1,  0, -1 };
static_assert(same(h4_row, grid::d4_nesw_row) && same(h4_col, grid::d4_nesw_col),
              "lsdtt_hillslope_operator_cpu.hpp:291 k_dr/k_dc is not grid::d4_nesw");

int g_failures = 0;

void check(const std::string& what, bool ok)
{
    std::cout << (ok ? "PASS " : "FAIL ") << what << "\n";
    if (!ok) ++g_failures;
}

std::uint64_t bits(double x) { std::uint64_t b; std::memcpy(&b, &x, sizeof b); return b; }
std::uint32_t bits(float x)  { std::uint32_t b; std::memcpy(&b, &x, sizeof b); return b; }

} // namespace

int main()
{
    std::cout << "== C01 neighbour tables and C19 sqrt(2) ==\n";
    check("39 neighbour tables (plus 5 hillslope face tables) equal their literals (static_assert at compile time)", true);
    check("d8_clockwise diagonals are the odd codes",
          grid::d8_clockwise_is_diagonal(1) && grid::d8_clockwise_is_diagonal(7) &&
          !grid::d8_clockwise_is_diagonal(0) && !grid::d8_clockwise_is_diagonal(6));

    // Double: the full literal, the nearest-double literal and std::sqrt(2.0).
    volatile double two = 2.0;   // keep std::sqrt from being folded away
    check("k_sqrt2_d == std::sqrt(2.0) (bits)", bits(grid::k_sqrt2_d) == bits(std::sqrt(static_cast<double>(two))));
    check("k_sqrt2_d == 1.4142135623730950488 (bits)", bits(grid::k_sqrt2_d) == bits(1.4142135623730950488));
    check("k_sqrt2_d == 1.4142135623730951 (bits)", bits(grid::k_sqrt2_d) == bits(1.4142135623730951));
    // Float: std::sqrt(2.0f) and the truncated float literals used at float sites.
    volatile float twof = 2.0f;
    check("k_sqrt2_f == std::sqrt(2.0f) (bits)", bits(grid::k_sqrt2_f) == bits(std::sqrt(static_cast<float>(twof))));
    check("k_sqrt2_f == 1.41421356237f (bits)", bits(grid::k_sqrt2_f) == bits(1.41421356237f));
    check("k_sqrt2_f == 1.41421356f (bits)", bits(grid::k_sqrt2_f) == bits(1.41421356f));
    check("k_sqrt2_f == static_cast<float>(1.41421356237) (bits)",
          bits(grid::k_sqrt2_f) == bits(static_cast<float>(1.41421356237)));
    // The two truncated DOUBLE literals are NOT the same double (decision D5).
    check("1.41421356237 (double) differs from k_sqrt2_d, as D5 records",
          bits(1.41421356237) != bits(grid::k_sqrt2_d));
    check("1.41421356 (double) differs from k_sqrt2_d, as D5 records",
          bits(1.41421356) != bits(grid::k_sqrt2_d));

    std::cout << "RESULT: " << (g_failures == 0 ? "PASS" : "FAIL") << "\n";
    return g_failures == 0 ? 0 : 1;
}
