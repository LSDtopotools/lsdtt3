#pragma once

// lsdtt_basin_outline_cpu
//
// Basin outline traced on the pixel-corner "divide grid".
//
// Closely based on LSDBasinOutline (LSDTopoTools2):
//  - Drainage divides are defined in-between DEM pixels (pixel corners)
//  - Outline is traced CLOCKWISE
//  - Periodic boundaries are handled by first making an "unwrapped" basin mask,
//    tracing in that unwrapped space, then folding the outline back to the
//    original domain.
//  - After tracing, the outline is rotated/opened so it starts at the basin
//    outlet cut corners.
//  - Diagonal outlet flow is a special case: start==end and the outline is
//    returned closed (first==last).
//
// lsdtt3 redesign:
//  - Construction takes an lsdtt_basin_cpu (single basin object)
//  - The basin is represented by a zero-copy slice into the flow s_vector.
//
// Copyright (C) 2026
// GPL-3.0

#include <string>
#include <vector>
#include <cstdint>

#include "../core/lsdtt_raster_meta.hpp"
#include "../core/lsdtt_types.hpp"
#include "../core/lsdtt_data_loader.hpp"     
#include "../utils/lsdtt_raster_geometry.hpp" 


namespace lsdtt
{

class lsdtt_basin_cpu;
class lsdtt_flow_base;

class lsdtt_basin_outline_cpu
{
public:
  // Builds the outline immediately.
  explicit lsdtt_basin_outline_cpu(const lsdtt_basin_cpu& basin);

  // Non-owning access
  [[nodiscard]] const lsdtt_basin_cpu& basin() const noexcept { return *basin_; }
  [[nodiscard]] const lsdtt_flow_base& flow() const noexcept { return *flow_; }

  // DEM meta
  [[nodiscard]] const lsdtt_raster_meta& meta() const noexcept { return meta_; }

  // Divide grid meta (full DEM, not cropped).
  [[nodiscard]] const lsdtt_raster_meta& divide_grid_meta() const noexcept { return divide_grid_meta_; }

  // Outline nodes (divide-grid row/col).
  // - Clockwise.
  // - Non-diagonal outlet: OPEN polyline start->...->end (inclusive).
  // - Diagonal outlet: CLOSED polyline with first==last.
  [[nodiscard]] const std::vector<lsdtt_pixel_rc>& outline_divide_rc() const noexcept
  {
    return outline_divide_rc_;
  }

  // Convenience getters (derived from basin pointer)
  [[nodiscard]] std::int64_t outlet_flat() const noexcept;

  [[nodiscard]] lsdtt_pixel_rc outlet_start_corner() const noexcept { return outlet_start_corner_; }
  [[nodiscard]] lsdtt_pixel_rc outlet_end_corner() const noexcept { return outlet_end_corner_; }
  [[nodiscard]] bool outlet_is_diagonal() const noexcept { return outlet_is_diagonal_; }

  [[nodiscard]] bool is_periodic_ns() const noexcept { return periodic_ns_; }
  [[nodiscard]] bool is_periodic_ew() const noexcept { return periodic_ew_; }

  // Function to get the 4 pixels surrounding a basin outline node, 
  // distinguishing between pixels in the basin and pixels outside the basin.
  void get_pixels_around_basin_outline_node(
      std::size_t outline_node_index,
      std::vector<std::int64_t>& in_basin,
      std::vector<std::int64_t>& out_basin
  ) const;


  // Functions to write to vector 

  // Returns one point feature per outline node.
  // Will have field values for basin_id (passed on) and cw_distance (double)
  // This can be used to make an lsdtt_vector_layer by higher-level code
  [[nodiscard]] std::vector<lsdtt_vector_record> get_outline_as_point_records(std::int32_t basin_id) const;

  // Returns one single line feature for the whole outline.
  // Will have field values only for basin_id (passed on)
  // This can be used to make an lsdtt_vector_layer by higher-level code.
  [[nodiscard]] lsdtt_vector_record get_outline_as_line_record(std::int32_t basin_id) const;

  // This returns the basin outline as a polygon.
  // Will have field values only for basin_id (passed on)
  // This can be used to make an lsdtt_vector_layer by higher-level code.
  [[nodiscard]] lsdtt_vector_record get_outline_as_polygon_record(std::int32_t basin_id) const;


private:
  const lsdtt_basin_cpu* basin_ = nullptr;   // non-owning; caller guarantees lifetime
  const lsdtt_flow_base* flow_  = nullptr;   // non-owning; owned by basin

  lsdtt_raster_meta meta_;
  lsdtt_raster_meta divide_grid_meta_;

  std::string bc_;
  bool periodic_ns_ = false;
  bool periodic_ew_ = false;

  // Outlet cut corners in DIVIDE grid coordinates (wrapped to divide grid)
  lsdtt_pixel_rc outlet_start_corner_{0, 0};
  lsdtt_pixel_rc outlet_end_corner_{0, 0};
  bool           outlet_is_diagonal_ = false;

  // Final outline after rotation/opening (divide-grid indices, wrapped)
  std::vector<lsdtt_pixel_rc> outline_divide_rc_;

private:
  void edge_tracing_normal_boundaries_();
  void edge_tracing_periodic_boundaries_();

  void compute_outlet_cut_corners_();
  void rotate_open_outline_to_outlet_cut_();

  // Divide grid meta depends on periodicity (matches LSDBasinOutline periodic handling)
  [[nodiscard]] static lsdtt_raster_meta make_divide_grid_meta_(
      const lsdtt_raster_meta& meta,
      bool periodic_ns,
      bool periodic_ew);

  // Helpers
  [[nodiscard]] static bool same_rc_(const lsdtt_pixel_rc& a, const lsdtt_pixel_rc& b) noexcept;
  static void remove_consecutive_duplicates_(std::vector<lsdtt_pixel_rc>& v);

  // Like LSDBasinOutline: scan full vector so last occurrence wins
  [[nodiscard]] static int find_node_index_last_(
      const std::vector<lsdtt_pixel_rc>& v,
      const lsdtt_pixel_rc& node) noexcept;

  [[nodiscard]] static int pos_mod_(int a, int m) noexcept;

  // Wrap a divide-grid index depending on periodicity; returns -1 if invalid in non-periodic case.
  [[nodiscard]] int wrap_div_row_(int r) const noexcept;
  [[nodiscard]] int wrap_div_col_(int c) const noexcept;

  // For baselevel outlet corner finding: treat out-of-bounds as nodata for non-periodic,
  // wrap for periodic, and within-bounds invalid mask cells as nodata.
  [[nodiscard]] bool neighbour_is_nodata_cellgrid_(int row_raw, int col_raw) const noexcept;
};

} // namespace lsdtt
