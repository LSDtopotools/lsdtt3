// ignore this comment
#include "lsdtt_basin_outline_cpu.hpp"

#include <algorithm>
#include <array>
#include <cctype>
#include <cstdint>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>
#include <cmath> 

#include "lsdtt_basin_cpu.hpp"
#include "lsdtt_flowinfo_cpu.hpp"

#include "../utils/lsdtt_array2d.hpp"
#include "../utils/lsdtt_binary_array2d_tracing.hpp" // provides trace_boundary_pixels(array2d_u8, int)

namespace lsdtt
{
namespace
{
// 3x3 offsets used for stamping into the double-resolution grid.
// (Same intent as LSDBasinOutline edge tracing.)
static constexpr int stamp9[9][2] = {
  {-1, -1}, {-1, 0}, {-1, 1},
  { 0, -1}, { 0, 0}, { 0, 1},
  { 1, -1}, { 1, 0}, { 1, 1}
};

// D8 direction offsets (lsdtt internal convention):
//   7 0 1
//   6 x 2
//   5 4 3
static constexpr std::array<int, 8> d8_row{{-1, -1,  0,  1,  1,  1,  0, -1}};
static constexpr std::array<int, 8> d8_col{{ 0,  1,  1,  1,  0, -1, -1, -1}};

} // anonymous namespace

// ----------------------------------------------------------------------------
// Public convenience
// ----------------------------------------------------------------------------
std::int64_t lsdtt_basin_outline_cpu::outlet_flat() const noexcept
{
  // basin outlet is the downstream-most pixel defining this basin (not necessarily a baselevel outlet)
  return basin_->get_basin_outlet_flat_index();
}

// ----------------------------------------------------------------------------
// Helpers
// ----------------------------------------------------------------------------
bool lsdtt_basin_outline_cpu::same_rc_(const lsdtt_pixel_rc& a,
                                      const lsdtt_pixel_rc& b) noexcept
{
  return a.row == b.row && a.col == b.col;
}

void lsdtt_basin_outline_cpu::remove_consecutive_duplicates_(std::vector<lsdtt_pixel_rc>& v)
{
  if (v.empty()) return;
  std::vector<lsdtt_pixel_rc> out;
  out.reserve(v.size());
  out.push_back(v.front());
  for (std::size_t i = 1; i < v.size(); ++i)
  {
    if (!same_rc_(v[i], out.back())) out.push_back(v[i]);
  }
  v.swap(out);
}

int lsdtt_basin_outline_cpu::find_node_index_last_(
    const std::vector<lsdtt_pixel_rc>& v,
    const lsdtt_pixel_rc& node) noexcept
{
  int idx = -1;
  for (int i = 0; i < static_cast<int>(v.size()); ++i)
  {
    if (same_rc_(v[i], node)) idx = i; // last wins (LSDBasinOutline behaviour)
  }
  return idx;
}

int lsdtt_basin_outline_cpu::pos_mod_(int a, int m) noexcept
{
  const int r = a % m;
  return (r < 0) ? (r + m) : r;
}

int lsdtt_basin_outline_cpu::wrap_div_row_(int r) const noexcept
{
  const int R = divide_grid_meta_.n_rows;
  if (periodic_ns_)
  {
    return pos_mod_(r, R);
  }
  if (r < 0 || r >= R) return -1;
  return r;
}

int lsdtt_basin_outline_cpu::wrap_div_col_(int c) const noexcept
{
  const int C = divide_grid_meta_.n_cols;
  if (periodic_ew_)
  {
    return pos_mod_(c, C);
  }
  if (c < 0 || c >= C) return -1;
  return c;
}

bool lsdtt_basin_outline_cpu::neighbour_is_nodata_cellgrid_(int row_raw, int col_raw) const noexcept
{
  const int nrows = meta_.n_rows;
  const int ncols = meta_.n_cols;

  int r = row_raw;
  int c = col_raw;

  if (periodic_ns_)
  {
    r = pos_mod_(r, nrows);
  }
  else
  {
    if (r < 0 || r >= nrows) return true;
  }

  if (periodic_ew_)
  {
    c = pos_mod_(c, ncols);
  }
  else
  {
    if (c < 0 || c >= ncols) return true;
  }

  // Within bounds: nodata if invalid (mask==0)
  return !flow_->is_valid_cell(r, c);
}

// ----------------------------------------------------------------------------
// Divide-grid meta
// ----------------------------------------------------------------------------
lsdtt_raster_meta
lsdtt_basin_outline_cpu::make_divide_grid_meta_(
    const lsdtt_raster_meta& meta,
    bool periodic_ns,
    bool periodic_ew)
{
  // Basin outlines are on a grid that is one pixel bigger and shifted half a
  // cell size to accommodate for the edge tracing along pixel edges.
  // (LSDBasinOutline rationale.)

  lsdtt_raster_meta dgm = meta;

  dgm.x_min = meta.x_min - 0.5 * meta.cell_size;
  dgm.y_min = meta.y_min - 0.5 * meta.cell_size;

  // LSDBasinOutline periodic handling:
  // - Non-periodic axis: corners are n+1
  // - Periodic axis: corners wrap so we keep n (no duplicate seam corner column/row)
  dgm.n_rows = meta.n_rows + (periodic_ns ? 0 : 1);
  dgm.n_cols = meta.n_cols + (periodic_ew ? 0 : 1);

  return dgm;
}

// ----------------------------------------------------------------------------
// Constructor
// ----------------------------------------------------------------------------
lsdtt_basin_outline_cpu::lsdtt_basin_outline_cpu(const lsdtt_basin_cpu& basin)
{
  basin_ = &basin;
  flow_  = &basin.flow();

  meta_ = flow_->meta();

  // Need boundary conditions to detect periodicity (CPU flowinfo only at present).
  const auto* fi = dynamic_cast<const lsdtt_flowinfo*>(flow_);
  if (!fi)
  {
    throw std::runtime_error(
        "lsdtt_basin_outline_cpu: flow backend is not lsdtt_flowinfo (CPU)."
    );
  }

  bc_ = fi->boundary_conditions();
  if (bc_.size() != 4)
  {
    throw std::runtime_error(
        "lsdtt_basin_outline_cpu: boundary_conditions must have length 4 (N,E,S,W)."
    );
  }

  auto lc = [](char ch) -> char {
    return static_cast<char>(std::tolower(static_cast<unsigned char>(ch)));
  };

  const char bcN = lc(bc_[0]);
  const char bcE = lc(bc_[1]);
  const char bcS = lc(bc_[2]);
  const char bcW = lc(bc_[3]);

  periodic_ns_ = (bcN == 'p' && bcS == 'p');
  periodic_ew_ = (bcE == 'p' && bcW == 'p');

  // LSDBasinOutline limitation: cannot deal with both axes periodic.
  if (periodic_ns_ && periodic_ew_)
  {
    std::cerr << "lsdtt_basin_outline_cpu: both NS and EW boundaries periodic!\n";
    std::cerr << "Can't deal with that. Turning off periodicity NS boundaries.\n";
    periodic_ns_ = false;
  }

  divide_grid_meta_ = make_divide_grid_meta_(meta_, periodic_ns_, periodic_ew_);

  // Edge tracing produces an outline ring starting at an arbitrary location.
  // These functions will deal with edge tracing and result in filled outline vectors.
  if (periodic_ns_ || periodic_ew_)
  {
    edge_tracing_periodic_boundaries_();
  }
  else
  {
    edge_tracing_normal_boundaries_();
  }

  // After tracing, the vectors just start at the top left due to edge tracing
  // but we want it to start at the basin outlet. (LSDBasinOutline.)
  compute_outlet_cut_corners_();
  rotate_open_outline_to_outlet_cut_();
}

// ----------------------------------------------------------------------------
// Edge tracing on non-periodic grid (LSDBasinOutline::edge_tracing_normal_boundaries)
// ----------------------------------------------------------------------------
void lsdtt_basin_outline_cpu::edge_tracing_normal_boundaries_()
{
  const lsdtt_flat_index_slice& basin_nodes = basin_->get_slice();

  if (basin_nodes.empty())
  {
    outline_divide_rc_.clear();
    return;
  }

  // Get row/col of basin nodes and min/max extents
  int min_r = meta_.n_rows;
  int max_r = 0;
  int min_c = meta_.n_cols;
  int max_c = 0;

  std::vector<int> r, c;
  r.reserve(basin_nodes.size());
  c.reserve(basin_nodes.size());

  for (const std::int64_t flat : basin_nodes)
  {
    const int row = flow_->row_from_flat_unchecked(flat);
    const int col = flow_->col_from_flat_unchecked(flat);
    r.push_back(row);
    c.push_back(col);

    min_r = std::min(min_r, row);
    max_r = std::max(max_r, row);
    min_c = std::min(min_c, col);
    max_c = std::max(max_c, col);
  }

  // Expand by 1 pixel in every direction (LSDBasinOutline)
  min_r -= 1; max_r += 1;
  min_c -= 1; max_c += 1;

  // Build double resolution binary grid (LSDBasinOutline comment):
  // 9 pixels associated with each DEM pixel; 4 lie on corners.
  const int NRowsDouble = (max_r - min_r + 1) * 2 + 1;
  const int NColsDouble = (max_c - min_c + 1) * 2 + 1;

  lsdtt_array2d_u8 binary(NRowsDouble, NColsDouble, 0u);

  for (std::size_t k = 0; k < r.size(); ++k)
  {
    const int base_r = (r[k] - min_r) * 2;
    const int base_c = (c[k] - min_c) * 2;

    for (int i = 0; i < 9; ++i)
    {
      const int rr = base_r + stamp9[i][0];
      const int cc = base_c + stamp9[i][1];
      binary(rr, cc) = 1u;
    }
  }

  // Trace boundary pixels (4-connectivity)
  const std::vector<lsdtt_pixel_rc> boundary = trace_boundary_pixels(binary, 4);

  // Keep only corner nodes (odd/odd), map back to divide grid coords
  outline_divide_rc_.clear();
  outline_divide_rc_.reserve(boundary.size());

  for (const auto& p : boundary)
  {
    if ((p.row & 1) && (p.col & 1))
    {
      const int div_r = ((p.row + 1) >> 1) + min_r;
      const int div_c = ((p.col + 1) >> 1) + min_c;
      outline_divide_rc_.push_back(lsdtt_pixel_rc{div_r, div_c});
    }
  }

  // Remove duplicate closure if present (LSDBasinOutline drops final point)
  if (outline_divide_rc_.size() > 1 &&
      same_rc_(outline_divide_rc_.front(), outline_divide_rc_.back()))
  {
    outline_divide_rc_.pop_back();
  }

  remove_consecutive_duplicates_(outline_divide_rc_);
}

// ----------------------------------------------------------------------------
// Edge tracing on periodic grid (LSDBasinOutline::edge_tracing_periodic_boundaries)
//
// Strategy (LSDBasinOutline):
//  - Make a binary array containing the "unwrapped" basin.
//  - Trace in unwrapped space.
//  - Fold corner coordinates back onto the original periodic domain.
// ----------------------------------------------------------------------------
void lsdtt_basin_outline_cpu::edge_tracing_periodic_boundaries_()
{
  const lsdtt_flat_index_slice& basin_nodes = basin_->get_slice();
  if (basin_nodes.empty())
  {
    outline_divide_rc_.clear();
    return;
  }

  const int nrows = meta_.n_rows;
  const int ncols = meta_.n_cols;

  // Local vectors
  std::vector<std::int64_t> flats;
  flats.reserve(basin_nodes.size());

  std::vector<int> r, c;
  r.reserve(basin_nodes.size());
  c.reserve(basin_nodes.size());

  // Flat -> local index
  std::unordered_map<std::int64_t, int> idx_of;
  idx_of.reserve(basin_nodes.size() * 2);

  int i_local = 0;
  for (const std::int64_t flat : basin_nodes)
  {
    flats.push_back(flat);
    idx_of[flat] = i_local++;

    r.push_back(flow_->row_from_flat_unchecked(flat));
    c.push_back(flow_->col_from_flat_unchecked(flat));
  }

  // grid_index of outlet is 0 (LSDBasinOutline)
  std::vector<int> grid_index(flats.size(), -1);
  grid_index[0] = 0;

  // Determine grid_index for each basin pixel based on seam crossings
  for (int i = 0; i < static_cast<int>(flats.size()); ++i)
  {
    if (i == 0) continue;

    const std::int64_t flat = flats[i];
    const std::int64_t rec  = flow_->receiver_of_unchecked(flat);

    auto it = idx_of.find(rec);
    if (it == idx_of.end())
    {
      // Should not happen for a well-defined basin slice; keep it robust.
      grid_index[i] = 0;
      continue;
    }

    int thisGridIndex = grid_index[it->second];

    const int row = r[i];
    const int col = c[i];
    const std::uint8_t flow_dir = flow_->receiver_dir_of_unchecked(flat);

    if (periodic_ew_ && col == ncols - 1)
    {
      if (flow_dir == 1 || flow_dir == 2 || flow_dir == 3)
      {
        // On periodic East boundary and flow to East
        thisGridIndex -= 1;
      }
    }
    if (periodic_ew_ && col == 0)
    {
      if (flow_dir == 5 || flow_dir == 6 || flow_dir == 7)
      {
        // On periodic West boundary and flow to West
        thisGridIndex += 1;
      }
    }
    if (periodic_ns_ && row == nrows - 1)
    {
      if (flow_dir == 3 || flow_dir == 4 || flow_dir == 5)
      {
        // On periodic South boundary and flow to South
        thisGridIndex += 1;
      }
    }
    if (periodic_ns_ && row == 0)
    {
      if (flow_dir == 7 || flow_dir == 0 || flow_dir == 1)
      {
        // On periodic North boundary and flow to North
        thisGridIndex -= 1;
      }
    }

    grid_index[i] = thisGridIndex;
  }

  const int min_grid = *std::min_element(grid_index.begin(), grid_index.end());

  // Determine unwrapped extents (LSDBasinOutline)
  int unwrapped_min_row = std::numeric_limits<int>::max();
  int unwrapped_max_row = std::numeric_limits<int>::min();
  int unwrapped_min_col = std::numeric_limits<int>::max();
  int unwrapped_max_col = std::numeric_limits<int>::min();

  for (std::size_t i = 0; i < flats.size(); ++i)
  {
    const int gidx = grid_index[i];

    const int urow = periodic_ns_ ? r[i] + (gidx - min_grid) * nrows : r[i];
    const int ucol = periodic_ew_ ? c[i] + (gidx - min_grid) * ncols : c[i];

    unwrapped_min_row = std::min(unwrapped_min_row, urow);
    unwrapped_max_row = std::max(unwrapped_max_row, urow);
    unwrapped_min_col = std::min(unwrapped_min_col, ucol);
    unwrapped_max_col = std::max(unwrapped_max_col, ucol);
  }

  // Expand by 1 DEM pixel in every direction (LSDBasinOutline)
  const int r0 = unwrapped_min_row - 1;
  const int r1 = unwrapped_max_row + 1;
  const int c0 = unwrapped_min_col - 1;
  const int c1 = unwrapped_max_col + 1;

  // Build the double-resolution grid
  const int nR = (r1 - r0 + 1) * 2 + 1;
  const int nC = (c1 - c0 + 1) * 2 + 1;

  lsdtt_array2d_u8 binary2(nR, nC, 0u);

  for (std::size_t i = 0; i < flats.size(); ++i)
  {
    const int gidx = grid_index[i];

    const int uR = periodic_ns_ ? r[i] + (gidx - min_grid) * nrows : r[i];
    const int uC = periodic_ew_ ? c[i] + (gidx - min_grid) * ncols : c[i];

    const int base_r = (uR - r0) * 2;
    const int base_c = (uC - c0) * 2;

    for (int k = 0; k < 9; ++k)
    {
      const int rr = base_r + stamp9[k][0];
      const int cc = base_c + stamp9[k][1];
      if (rr >= 0 && rr < nR && cc >= 0 && cc < nC)
      {
        binary2(rr, cc) = 1u;
      }
    }
  }

  // Trace boundary pixels
  const std::vector<lsdtt_pixel_rc> boundary = trace_boundary_pixels(binary2, 4);

  // Fold back onto original periodic divide grid (LSDBasinOutline)
  outline_divide_rc_.clear();
  outline_divide_rc_.reserve(boundary.size());

  const int NRowsDiv = divide_grid_meta_.n_rows;
  const int NColsDiv = divide_grid_meta_.n_cols;

  auto wrap = [&](int idx, int M) -> int
  {
    return pos_mod_(idx, M);
  };

  for (const auto& p : boundary)
  {
    if ((p.row & 1) && (p.col & 1))
    {
      int corner_r = ((p.row + 1) >> 1) + r0;
      int corner_c = ((p.col + 1) >> 1) + c0;

      // fold back
      if (periodic_ns_) corner_r = wrap(corner_r, NRowsDiv);
      if (periodic_ew_) corner_c = wrap(corner_c, NColsDiv);

      outline_divide_rc_.push_back(lsdtt_pixel_rc{corner_r, corner_c});
    }
  }

  // Drop the last point (duplicate of first), like LSDBasinOutline
  if (outline_divide_rc_.size() > 1 &&
      same_rc_(outline_divide_rc_.front(), outline_divide_rc_.back()))
  {
    outline_divide_rc_.pop_back();
  }

  remove_consecutive_duplicates_(outline_divide_rc_);
}

// ----------------------------------------------------------------------------
// Function to get the 4 pixels surrounding a basin outline node, 
// distinguishing between pixels in the basin and pixels outside the basin.
// ----------------------------------------------------------------------------
void lsdtt_basin_outline_cpu::get_pixels_around_basin_outline_node(
    std::size_t outline_node_index,
    std::vector<std::int64_t>& in_basin,
    std::vector<std::int64_t>& out_basin
) const
{
  in_basin.clear();
  out_basin.clear();

  const int i = static_cast<int>(outline_node_index);
  const int NNodes_outline = static_cast<int>(outline_divide_rc_.size());

  // Get the nodes of the previous and next outline nodes too
  int next_i = 0, prev_i = 0;
  if (i == NNodes_outline - 1) {
      // Next node is at the start of outline
      next_i = 0; 
      // But if diagonal outlet flow, this is actually the same node,
      // so advance one more.
      if (same_rc_(outline_divide_rc_[next_i], outline_divide_rc_[i])) {
          next_i = (next_i + 1) % NNodes_outline;
      }
  } else { next_i = i + 1; }

  if (i == 0) {
      // Previous node is at the back of outline
      prev_i = NNodes_outline - 1;
      // But if diagonal outlet flow, this is actually the same node,
      // so go one further back.
      if (same_rc_(outline_divide_rc_[i], outline_divide_rc_[prev_i])) {
          prev_i = (prev_i - 1 + NNodes_outline) % NNodes_outline;
      }
  } else { prev_i = i - 1; }

  // Get the row and col of prev_i, i, next_i
  int prev_i_row = outline_divide_rc_[prev_i].row;
  int prev_i_col = outline_divide_rc_[prev_i].col;
  int i_row      = outline_divide_rc_[i].row;
  int i_col      = outline_divide_rc_[i].col;
  int next_i_row = outline_divide_rc_[next_i].row;
  int next_i_col = outline_divide_rc_[next_i].col;

  // Local lambdas for wrapped differences in case of periodic boundaries
  const int n_rows = meta_.n_rows;
  const int n_cols = meta_.n_cols;

  auto wrapped_row_diff = [&](int from, int to) -> int {
      if (!periodic_ns_) return to - from;
      int diff = to - from;
      if (diff > n_rows / 2) diff -= n_rows;
      if (diff < -n_rows / 2) diff += n_rows;
      return diff;
  };

  auto wrapped_col_diff = [&](int from, int to) -> int {
      if (!periodic_ew_) return to - from;
      int diff = to - from;
      if (diff > n_cols / 2) diff -= n_cols;
      if (diff < -n_cols / 2) diff += n_cols;
      return diff;
  };

  // Now check if the outline edges coming in and out of i have direction
  // north (0), east (1), south (2) or west (3).
  int edge_direction_in = -1, edge_direction_out = -1;

  int d_row_out = wrapped_row_diff(i_row, next_i_row);
  int d_col_out = wrapped_col_diff(i_col, next_i_col);

  if (d_row_out < 0)      edge_direction_out = 0; // North
  else if (d_row_out > 0) edge_direction_out = 2; // South
  else if (d_col_out > 0) edge_direction_out = 1; // East
  else if (d_col_out < 0) edge_direction_out = 3; // West

  int d_row_in = wrapped_row_diff(prev_i_row, i_row);
  int d_col_in = wrapped_col_diff(prev_i_col, i_col);

  if (d_row_in < 0)      edge_direction_in = 0; // North
  else if (d_row_in > 0) edge_direction_in = 2; // South
  else if (d_col_in > 0) edge_direction_in = 1; // East
  else if (d_col_in < 0) edge_direction_in = 3; // West

  // Get the x's and y's indices of the 4 pixels of the flow grid
  int raw_x1 = i_col - 1;
  int raw_y1 = i_row - 1;
  int raw_x2 = i_col;
  int raw_y2 = i_row;

  // Account for periodicity (or bounds)
  auto wrap_idx = [](int idx, int max_val, bool periodic) -> int {
      if (periodic) return (idx % max_val + max_val) % max_val;
      if (idx < 0 || idx >= max_val) return -1;
      return idx;
  };

  int x1 = wrap_idx(raw_x1, n_cols, periodic_ew_);
  int x2 = wrap_idx(raw_x2, n_cols, periodic_ew_);
  int y1 = wrap_idx(raw_y1, n_rows, periodic_ns_);
  int y2 = wrap_idx(raw_y2, n_rows, periodic_ns_);

  std::vector<std::int64_t> dem_pixels_in_basin, dem_pixels_out_basin;
  auto add_unique = [&](std::vector<std::int64_t>& vec, std::int64_t val) {
      if (val != -1 && std::find(vec.begin(), vec.end(), val) == vec.end()) {
          vec.push_back(val);
      }
  };

  auto get_node = [&](int r, int c) -> std::int64_t {
      if (r == -1 || c == -1) return -1;
      std::int64_t flat = static_cast<std::int64_t>(r) * n_cols + c;
      if (!flow_->is_valid_cell(flat)) return -1;
      return flat;
  };

  // Figure out which of the four nodes surrounding the edge are in the 
  // basin, using the fact that outlines are traced CLOCKWISE
  std::int64_t pix_inb = -1, pix_outb = -1;

  // Lets start with outgoing node
  if (edge_direction_out == 0) {
      // North: (y1,x2) in basin, (y1, x1) out of basin
      pix_inb  = get_node(y1, x2);
      pix_outb = get_node(y1, x1);
  } else if (edge_direction_out == 1) {
      // East: (y2,x2) in basin, (y1,x2) out of basin
      pix_inb  = get_node(y2, x2);
      pix_outb = get_node(y1, x2);
  } else if (edge_direction_out == 2) {
      // South: (y2,x1) in basin, (y2,x2) out of basin
      pix_inb  = get_node(y2, x1);
      pix_outb = get_node(y2, x2);
  } else if (edge_direction_out == 3) {
      // West: (y1,x1) in basin, (y2,x1) out of basin
      pix_inb  = get_node(y1, x1);
      pix_outb = get_node(y2, x1);
  }
  if (pix_inb != -1)  add_unique(dem_pixels_in_basin, pix_inb);
  if (pix_outb != -1) add_unique(dem_pixels_out_basin, pix_outb);

  // Now the ingoing node
  if (edge_direction_in == 0) {
      pix_inb  = get_node(y2, x2);
      pix_outb = get_node(y2, x1);
  } else if (edge_direction_in == 1) {
      pix_inb  = get_node(y2, x1);
      pix_outb = get_node(y1, x1);
  } else if (edge_direction_in == 2) {
      pix_inb  = get_node(y1, x1);
      pix_outb = get_node(y1, x2);
  } else if (edge_direction_in == 3) {
      pix_inb  = get_node(y1, x2);
      pix_outb = get_node(y2, x2);
  }
  if (pix_inb != -1)  add_unique(dem_pixels_in_basin, pix_inb);
  if (pix_outb != -1) add_unique(dem_pixels_out_basin, pix_outb);

  // Corner logic with guards
  if (edge_direction_in != edge_direction_out && edge_direction_in != -1 && edge_direction_out != -1) {
      auto safe_insert = [&](int y, int x, bool in_b) {
          std::int64_t corner_pix = get_node(y, x);
          if (corner_pix != -1) {
              if (in_b) add_unique(dem_pixels_in_basin, corner_pix);
              else      add_unique(dem_pixels_out_basin, corner_pix);
          }
      };
      if (edge_direction_in == 0 && edge_direction_out == 1) safe_insert(y1, x1, false);
      if (edge_direction_in == 0 && edge_direction_out == 3) safe_insert(y1, x2, true);
      if (edge_direction_in == 1 && edge_direction_out == 0) safe_insert(y2, x2, true);
      if (edge_direction_in == 1 && edge_direction_out == 2) safe_insert(y1, x2, false);
      if (edge_direction_in == 2 && edge_direction_out == 1) safe_insert(y2, x1, true);
      if (edge_direction_in == 2 && edge_direction_out == 3) safe_insert(y2, x2, false);
      if (edge_direction_in == 3 && edge_direction_out == 2) safe_insert(y1, x1, true);
      if (edge_direction_in == 3 && edge_direction_out == 0) safe_insert(y2, x1, false);
  }

  // If the outlet node was accidentally put in the out_basin set,
  // remove it and ensure it is in the in_basin set.
  std::int64_t outlet_node = basin_->get_basin_outlet_flat_index();
  auto it = std::find(dem_pixels_out_basin.begin(), dem_pixels_out_basin.end(), outlet_node);
  if (it != dem_pixels_out_basin.end()) {
      dem_pixels_out_basin.erase(it);
      add_unique(dem_pixels_in_basin, outlet_node);
  }

  out_basin = dem_pixels_out_basin;
  in_basin  = dem_pixels_in_basin;
}

// ----------------------------------------------------------------------------
// Compute outlet cut corners (LSDBasinOutline::create logic)
//
// Uses receiver direction to decide whether this is:
//  - baselevel outlet (receiver is outlet_receiver)
//  - diagonal exit (start==end)
//  - cardinal exit (two adjacent corners on the exit edge)
// ----------------------------------------------------------------------------
void lsdtt_basin_outline_cpu::compute_outlet_cut_corners_()
{
  const std::int64_t out_flat = outlet_flat();

  const int outlet_row = flow_->row_from_flat_unchecked(out_flat);
  const int outlet_col = flow_->col_from_flat_unchecked(out_flat);

  const std::int64_t rec = flow_->receiver_of_unchecked(out_flat);
  const std::uint8_t dir = flow_->receiver_dir_of_unchecked(out_flat);

  int rowchange = 0;
  int colchange = 0;

  // Baselevel: receiver is outlet_receiver (lsdtt3 convention) or no_direction.
  if (rec == lsdtt_flow_base::outlet_receiver ||
      dir == lsdtt_flow_base::no_direction)
  {
    rowchange = 0;
    colchange = 0;
  }
  else
  {
    rowchange = d8_row[dir];
    colchange = d8_col[dir];
  }

  // 1) Baselevel outlet: find entry/exit corners on DEM edge by nodata transitions
  if (rowchange == 0 && colchange == 0)
  {
    // Neighbour scan offsets (CW)
    const std::vector<std::pair<int,int>> neighbor_offsets = {
      {-1,-1},{-1, 0},{-1, 1},
      { 0, 1},{ 1, 1},{ 1, 0},
      { 1,-1},{ 0,-1}
    };

    int startrowDiv = -1, startcolDiv = -1;
    int endrowDiv   = -1, endcolDiv   = -1;

    for (std::size_t i = 0; i < neighbor_offsets.size(); ++i)
    {
      const int dr1 = neighbor_offsets[i].first;
      const int dc1 = neighbor_offsets[i].second;

      const std::size_t ni = (i + 1) % neighbor_offsets.size();
      const int dr2 = neighbor_offsets[ni].first;
      const int dc2 = neighbor_offsets[ni].second;

      const bool nodata1 = neighbour_is_nodata_cellgrid_(outlet_row + dr1, outlet_col + dc1);
      const bool nodata2 = neighbour_is_nodata_cellgrid_(outlet_row + dr2, outlet_col + dc2);

      // entry point onto DEM
      if (nodata1 && !nodata2)
      {
        const int sum_dr = dr1 + dr2;
        const int sum_dc = dc1 + dc2;

        const int raw_row = outlet_row + (sum_dr >= 0 ? 1 : 0);
        const int raw_col = outlet_col + (sum_dc >= 0 ? 1 : 0);

        startrowDiv = wrap_div_row_(raw_row);
        startcolDiv = wrap_div_col_(raw_col);
      }

      // exit point off DEM
      if (!nodata1 && nodata2)
      {
        const int sum_dr = dr1 + dr2;
        const int sum_dc = dc1 + dc2;

        const int raw_row = outlet_row + (sum_dr >= 0 ? 1 : 0);
        const int raw_col = outlet_col + (sum_dc >= 0 ? 1 : 0);

        endrowDiv = wrap_div_row_(raw_row);
        endcolDiv = wrap_div_col_(raw_col);
      }
    }

    if (startrowDiv < 0 || startcolDiv < 0 || endrowDiv < 0 || endcolDiv < 0)
    {
      throw std::runtime_error("lsdtt_basin_outline_cpu: failed to infer baselevel outlet cut corners.");
    }

    outlet_is_diagonal_ = false;
    outlet_start_corner_ = lsdtt_pixel_rc{startrowDiv, startcolDiv};
    outlet_end_corner_   = lsdtt_pixel_rc{endrowDiv,   endcolDiv};
    return;
  }

  // 2) Diagonal flow: start==end at the exit corner
  if (rowchange != 0 && colchange != 0)
  {
    const int dr = (rowchange > 0 ? 1 : 0);
    const int dc = (colchange > 0 ? 1 : 0);

    const int rr = wrap_div_row_(outlet_row + dr);
    const int cc = wrap_div_col_(outlet_col + dc);

    outlet_is_diagonal_ = true;
    outlet_start_corner_ = lsdtt_pixel_rc{rr, cc};
    outlet_end_corner_   = outlet_start_corner_;
    return;
  }

  // 3) Horizontal flow
  if (rowchange == 0 && colchange != 0)
  {
    outlet_is_diagonal_ = false;

    if (colchange > 0)
    {
      // flow â E : corners SE then NE
      outlet_start_corner_ = lsdtt_pixel_rc{
        wrap_div_row_(outlet_row + 1),
        wrap_div_col_(outlet_col + 1)
      };
      outlet_end_corner_ = lsdtt_pixel_rc{
        wrap_div_row_(outlet_row),
        wrap_div_col_(outlet_col + 1)
      };
    }
    else
    {
      // flow â W : corners NW then SW
      outlet_start_corner_ = lsdtt_pixel_rc{
        wrap_div_row_(outlet_row),
        wrap_div_col_(outlet_col)
      };
      outlet_end_corner_ = lsdtt_pixel_rc{
        wrap_div_row_(outlet_row + 1),
        wrap_div_col_(outlet_col)
      };
    }

    return;
  }

  // 4) Vertical flow
  if (rowchange != 0 && colchange == 0)
  {
    outlet_is_diagonal_ = false;

    if (rowchange > 0)
    {
      // flow â S : corners SW then SE
      outlet_start_corner_ = lsdtt_pixel_rc{
        wrap_div_row_(outlet_row + 1),
        wrap_div_col_(outlet_col)
      };
      outlet_end_corner_ = lsdtt_pixel_rc{
        wrap_div_row_(outlet_row + 1),
        wrap_div_col_(outlet_col + 1)
      };
    }
    else
    {
      // flow â N : corners NE then NW
      outlet_start_corner_ = lsdtt_pixel_rc{
        wrap_div_row_(outlet_row),
        wrap_div_col_(outlet_col + 1)
      };
      outlet_end_corner_ = lsdtt_pixel_rc{
        wrap_div_row_(outlet_row),
        wrap_div_col_(outlet_col)
      };
    }

    return;
  }

  // Fallback (should not happen)
  outlet_is_diagonal_ = false;
  outlet_start_corner_ = lsdtt_pixel_rc{wrap_div_row_(outlet_row), wrap_div_col_(outlet_col)};
  outlet_end_corner_   = lsdtt_pixel_rc{wrap_div_row_(outlet_row), wrap_div_col_(outlet_col + 1)};
}

// ----------------------------------------------------------------------------
// Rotate/open outline so it starts at outlet cut (LSDBasinOutline::create logic)
//
// Outline is traced CLOCKWISE; we keep that order and just reindex.
//
// - Diagonal case (start==end): do a full cycle then repeat the start corner
// - General case: walk from start -> end (inclusive), wrapping if needed
// ----------------------------------------------------------------------------
void lsdtt_basin_outline_cpu::rotate_open_outline_to_outlet_cut_()
{
  const int M = static_cast<int>(outline_divide_rc_.size());
  if (M == 0) return;

  const int start_idx = find_node_index_last_(outline_divide_rc_, outlet_start_corner_);
  const int end_idx   = find_node_index_last_(outline_divide_rc_, outlet_end_corner_);

  if (start_idx < 0 || end_idx < 0)
  {
    std::cerr << "ERROR lsdtt_basin_outline_cpu: could not find start or end corner in outline\n";
    std::cerr << "Outlet flat: " << outlet_flat() << "\n";
    std::cerr << "Start corner: " << outlet_start_corner_.row << "/" << outlet_start_corner_.col << "\n";
    std::cerr << "End corner:   " << outlet_end_corner_.row   << "/" << outlet_end_corner_.col << "\n";
    throw std::runtime_error("lsdtt_basin_outline_cpu: outlet cut corners not found in outline.");
  }

  std::vector<lsdtt_pixel_rc> new_outline;
  new_outline.reserve(static_cast<std::size_t>(M + 1));

  if (start_idx == end_idx)
  {
    // Diagonal case: start==end -> full cycle then close by repeating start
    for (int k = 0; k < M; ++k)
    {
      const int idx = (start_idx + k) % M;
      new_outline.push_back(outline_divide_rc_[idx]);
    }
    new_outline.push_back(outline_divide_rc_[start_idx]); // close
  }
  else
  {
    // General case: walk start -> end (inclusive), forward (clockwise), wrapping if needed
    int i = start_idx;
    while (true)
    {
      new_outline.push_back(outline_divide_rc_[i]);
      if (i == end_idx) break;
      i = (i + 1) % M;
    }
  }

  outline_divide_rc_.swap(new_outline);
}

std::vector<lsdtt_vector_record>
lsdtt_basin_outline_cpu::get_outline_as_point_records(std::int32_t basin_id) const
{
  std::vector<lsdtt_vector_record> records;
  records.reserve(outline_divide_rc_.size());

  const double dx = meta_.cell_size; // clockwise spacing along pixel edges

  for (std::size_t i = 0; i < outline_divide_rc_.size(); ++i)
  {
    const lsdtt_pixel_rc rc = outline_divide_rc_[i];

    // Divide-grid node centre is located at DEM pixel corner position
    const lsdtt_xy xy = cell_center_from_row_col(divide_grid_meta_, rc.row, rc.col);

    const double cw_distance = static_cast<double>(i) * dx;

    lsdtt_vector_record rec;
    rec.geometry = lsdtt_point_geom{xy.x, xy.y};

    // NOTE: field_values must match caller-provided schema order.
    // basin_id is int32 in API but stored as int64 in the variant.
    rec.field_values = {
      static_cast<std::int64_t>(basin_id),
      cw_distance
    };

    records.push_back(std::move(rec));
  }

  return records;
}

lsdtt_vector_record
lsdtt_basin_outline_cpu::get_outline_as_line_record(std::int32_t basin_id) const
{
  lsdtt_vector_record rec;

  // Set attributes (match schema order)
  rec.field_values = { static_cast<std::int64_t>(basin_id) };

  // Empty outline => null geometry record
  if (outline_divide_rc_.size() < 2)
  {
    rec.geometry = std::monostate{};
    return rec;
  }

  lsdtt_line_geom line;
  line.points.reserve(outline_divide_rc_.size());

  for (const lsdtt_pixel_rc& rc : outline_divide_rc_)
  {
    const lsdtt_xy xy = cell_center_from_row_col(divide_grid_meta_, rc.row, rc.col);
    line.points.push_back(xy);
  }

  rec.geometry = std::move(line);
  return rec;
}

lsdtt_vector_record
lsdtt_basin_outline_cpu::get_outline_as_polygon_record(std::int32_t basin_id) const
{
    lsdtt_vector_record rec;
    rec.field_values = { static_cast<std::int64_t>(basin_id) };

    if (outline_divide_rc_.size() < 2)
    {
        rec.geometry = std::monostate{};
        return rec;
    }

    lsdtt_polygon_geom poly;
    poly.exterior.reserve(outline_divide_rc_.size() + 1);

    for (const lsdtt_pixel_rc& rc : outline_divide_rc_)
    {
        const lsdtt_xy xy = cell_center_from_row_col(divide_grid_meta_, rc.row, rc.col);
        poly.exterior.push_back(xy);
    }

    // Close the ring â the outline starts at the outlet cut and ends there
    // too after rotate_open_outline_to_outlet_cut_, but the closing point
    // is not repeated in outline_divide_rc_ itself.
    poly.exterior.push_back(poly.exterior.front());

    rec.geometry = std::move(poly);
    return rec;
}


// ----------------------------------------------------------------------------
// Ridgeline Generation Methods
// ----------------------------------------------------------------------------

std::vector<float> lsdtt_basin_outline_cpu::interpolate_pixel_values_on_outline(
    const lsdtt_raster_cpu& data_raster) const
{
  const int NNodes_outline = static_cast<int>(outline_divide_rc_.size());
  std::vector<float> interpolated_values(NNodes_outline, data_raster.get_meta().no_data_value);

  const int n_rows_data = data_raster.get_meta().n_rows;
  const int n_cols_data = data_raster.get_meta().n_cols;
  const auto& data = data_raster.get_data_array();
  const auto& mask = data_raster.get_mask_array();

  for (int i = 0; i < NNodes_outline; ++i)
  {
    float x_interp = static_cast<float>(outline_divide_rc_[i].col) - 0.5f;
    float y_interp = static_cast<float>(outline_divide_rc_[i].row) - 0.5f;

    int raw_x1 = static_cast<int>(std::floor(x_interp));
    int raw_y1 = static_cast<int>(std::floor(y_interp));
    int raw_x2 = raw_x1 + 1;
    int raw_y2 = raw_y1 + 1;

    int x1 = raw_x1; int x2 = raw_x2;
    int y1 = raw_y1; int y2 = raw_y2;

    if (periodic_ew_) {
      x1 = pos_mod_(raw_x1, n_cols_data);
      x2 = pos_mod_(raw_x2, n_cols_data);
    } else {
      if (x1 < 0) x1 = -1;
      if (x2 >= n_cols_data) x2 = -1;
    }

    if (periodic_ns_) {
      y1 = pos_mod_(raw_y1, n_rows_data);
      y2 = pos_mod_(raw_y2, n_rows_data);
    } else {
      if (y1 < 0) y1 = -1;
      if (y2 >= n_rows_data) y2 = -1;
    }

    float Q11 = -9999.0f, Q12 = -9999.0f, Q21 = -9999.0f, Q22 = -9999.0f;
    bool valid_Q11 = (x1 != -1 && y1 != -1 && mask(y1, x1) != 0);
    bool valid_Q12 = (x1 != -1 && y2 != -1 && mask(y2, x1) != 0);
    bool valid_Q21 = (x2 != -1 && y1 != -1 && mask(y1, x2) != 0);
    bool valid_Q22 = (x2 != -1 && y2 != -1 && mask(y2, x2) != 0);

    if (valid_Q11) Q11 = data(y1, x1);
    if (valid_Q12) Q12 = data(y2, x1);
    if (valid_Q21) Q21 = data(y1, x2);
    if (valid_Q22) Q22 = data(y2, x2);

    int valid_count = valid_Q11 + valid_Q12 + valid_Q21 + valid_Q22;

    if (valid_count == 4) {
      float fxy1 = (Q11*(x2 - x_interp) + Q21*(x_interp - x1)) / (x2 - x1);
      float fxy2 = (Q12*(x2 - x_interp) + Q22*(x_interp - x1)) / (x2 - x1);
      interpolated_values[i] = (fxy1*(y2 - y_interp) + fxy2*(y_interp - y1)) / (y2 - y1);
    }
    else if (valid_count == 3) {
      float sum = 0.0f; int count = 0;
      if (valid_Q11) { sum += Q11; count++; }
      if (valid_Q12) { sum += Q12; count++; }
      if (valid_Q21) { sum += Q21; count++; }
      if (valid_Q22) { sum += Q22; count++; }
      interpolated_values[i] = sum / count;
    }
    else if (valid_count == 2) {
      if (valid_Q11 && valid_Q21) interpolated_values[i] = (Q11*(x2 - x_interp) + Q21*(x_interp - x1)) / (x2 - x1);
      else if (valid_Q12 && valid_Q22) interpolated_values[i] = (Q12*(x2 - x_interp) + Q22*(x_interp - x1)) / (x2 - x1);
      else if (valid_Q11 && valid_Q12) interpolated_values[i] = (Q11*(y2 - y_interp) + Q12*(y_interp - y1)) / (y2 - y1);
      else if (valid_Q21 && valid_Q22) interpolated_values[i] = (Q21*(y2 - y_interp) + Q22*(y_interp - y1)) / (y2 - y1);
      else continue;
    }
    else if (valid_count == 1) {
      if (valid_Q11) interpolated_values[i] = Q11;
      else if (valid_Q12) interpolated_values[i] = Q12;
      else if (valid_Q21) interpolated_values[i] = Q21;
      else if (valid_Q22) interpolated_values[i] = Q22;
    }
  }
  return interpolated_values;
}

lsdtt_ridgeline_cpu lsdtt_basin_outline_cpu::write_outline_portion_to_ridgeline(
    int start_index, int stop_index, bool write_in_reverse) const
{
  std::vector<lsdtt_pixel_rc> rl_nodes;
  if (!write_in_reverse) {
    for (int i = start_index; i < stop_index; ++i) {
      rl_nodes.push_back(outline_divide_rc_[i]);
    }
  } else {
    for (int i = start_index; i >= stop_index + 1; --i) {
      rl_nodes.push_back(outline_divide_rc_[i]);
    }
  }
  return lsdtt_ridgeline_cpu(std::move(rl_nodes), divide_grid_meta_);
}

std::vector<lsdtt_ridgeline_cpu> lsdtt_basin_outline_cpu::split_outline_into_ridgelines_by_max_elevation(
    const lsdtt_raster_cpu& elevation_raster) const
{
  std::vector<float> outline_elev = interpolate_pixel_values_on_outline(elevation_raster);
  auto max_elev_iter = std::max_element(outline_elev.begin(), outline_elev.end());
  int current_index = static_cast<int>(std::distance(outline_elev.begin(), max_elev_iter));

  std::vector<lsdtt_ridgeline_cpu> ridgelines;
  ridgelines.push_back(write_outline_portion_to_ridgeline(0, current_index + 1, false));
  ridgelines.push_back(write_outline_portion_to_ridgeline(static_cast<int>(outline_divide_rc_.size()) - 1, current_index - 1, true));
  return ridgelines;
}

std::vector<lsdtt_ridgeline_cpu> lsdtt_basin_outline_cpu::split_outline_into_ridgelines_by_longest_channel(
    const lsdtt_raster_cpu& dist_from_outlet) const
{
  std::int64_t ch_node = -1;
  float max_dist = -1.0f;
  std::int64_t outlet_flat = basin_->get_basin_outlet_flat_index();
  const auto& ddata = dist_from_outlet.get_data_array();
        
  for (std::int64_t flat : flow_->upslope_nodes(outlet_flat)) {
    float dist = ddata.data()[flat];
    if (dist > max_dist) {
      max_dist = dist;
      ch_node = flat;
    }
  }

  int ch_row = flow_->row_from_flat_unchecked(ch_node);
  int ch_col = flow_->col_from_flat_unchecked(ch_node);
        
  float current_euclidean_distance = -1.0f;
  int current_index = 0;

  for (std::size_t i = 0; i < outline_divide_rc_.size(); ++i) {
    int outline_row = outline_divide_rc_[i].row * 2;
    int outline_col = outline_divide_rc_[i].col * 2;
    int mapped_ch_row = ch_row * 2 + 1;
    int mapped_ch_col = ch_col * 2 + 1;

    float dist = std::sqrt(std::pow(outline_row - mapped_ch_row, 2) + std::pow(outline_col - mapped_ch_col, 2));
    if (i == 0 || dist < current_euclidean_distance) {
      current_euclidean_distance = dist;
      current_index = static_cast<int>(i);
    }
  }

  std::vector<lsdtt_ridgeline_cpu> ridgelines;
  ridgelines.push_back(write_outline_portion_to_ridgeline(0, current_index + 1, false));
  ridgelines.push_back(write_outline_portion_to_ridgeline(static_cast<int>(outline_divide_rc_.size()) - 1, current_index - 1, true));
  return ridgelines;
}

std::vector<lsdtt_ridgeline_cpu> lsdtt_basin_outline_cpu::split_outline_into_ridgelines_by_equal_length() const
{
  if (outline_divide_rc_.empty()) return {};

  double total_perimeter = 0.0;
  std::vector<double> dists(outline_divide_rc_.size(), 0.0);
  for (std::size_t i = 1; i < outline_divide_rc_.size(); ++i) {
    lsdtt_xy p1 = cell_center_from_row_col(divide_grid_meta_, outline_divide_rc_[i-1].row, outline_divide_rc_[i-1].col);
    lsdtt_xy p2 = cell_center_from_row_col(divide_grid_meta_, outline_divide_rc_[i].row, outline_divide_rc_[i].col);
    double dx = p2.x - p1.x;
    double dy = p2.y - p1.y;
    total_perimeter += std::sqrt(dx*dx + dy*dy);
    dists[i] = total_perimeter;
  }

  int current_index = 0;
  while (current_index < static_cast<int>(dists.size()) && dists[current_index] < total_perimeter / 2.0) {
    current_index++;
  }

  std::vector<lsdtt_ridgeline_cpu> ridgelines;
  ridgelines.push_back(write_outline_portion_to_ridgeline(0, current_index + 1, false));
  ridgelines.push_back(write_outline_portion_to_ridgeline(static_cast<int>(outline_divide_rc_.size()) - 1, current_index, true));
  return ridgelines;
}

lsdtt_ridgeline_cpu lsdtt_basin_outline_cpu::get_left_transverse_divide_as_ridgeline(
    std::int32_t basin_junction,
    const std::vector<std::int32_t>& rsb_N,
    const std::vector<std::int32_t>& rsb_S,
    const lsdtt_junction_network_base& junction_network) const
{
  std::vector<int> start_end = find_left_transverse_divide_nodes(basin_junction, rsb_N, rsb_S, junction_network);
  if (start_end.size() != 2) return lsdtt_ridgeline_cpu();
  return write_outline_portion_to_ridgeline(start_end[0], start_end[1], false);
}

lsdtt_ridgeline_cpu lsdtt_basin_outline_cpu::get_right_transverse_divide_as_ridgeline(
    std::int32_t basin_junction,
    const std::vector<std::int32_t>& rsb_N,
    const std::vector<std::int32_t>& rsb_S,
    const lsdtt_junction_network_base& junction_network) const
{
  std::vector<int> start_end = find_right_transverse_divide_nodes(basin_junction, rsb_N, rsb_S, junction_network);
  if (start_end.size() != 2) return lsdtt_ridgeline_cpu();
  return write_outline_portion_to_ridgeline(start_end[0], start_end[1], true);
}

lsdtt_ridgeline_cpu lsdtt_basin_outline_cpu::get_outline_nodes_on_main_drainage_divide_as_ridgeline(
    std::int32_t basin_junction,
    const std::vector<std::int32_t>& rsb_N,
    const std::vector<std::int32_t>& rsb_S,
    const lsdtt_junction_network_base& junction_network) const
{
  std::vector<int> start_end = find_main_drainage_divide_nodes(basin_junction, rsb_N, rsb_S, junction_network);
  if (start_end.size() != 2) return lsdtt_ridgeline_cpu();
  return write_outline_portion_to_ridgeline(start_end[0], start_end[1], false);
}

lsdtt_ridgeline_cpu lsdtt_basin_outline_cpu::get_interfluvial_divide_as_ridgeline(
    std::int32_t basin_junction,
    std::int32_t otherBasin_junction_number,
    const lsdtt_junction_network_base& junction_network) const
{
  std::vector<int> start_end = find_interfluvial_divide_nodes(basin_junction, otherBasin_junction_number, junction_network);
  if (start_end.size() != 2) return lsdtt_ridgeline_cpu();
  return write_outline_portion_to_ridgeline(start_end[0], start_end[1], false);
}

// -------------------------------------------------------------------------
// Internal topological finders for divides
// -------------------------------------------------------------------------


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// This scans the outline for triple divide junctions where the outline 
// (potentially) joins and leaves the left transverse divide.
// By iterating through the outline clockwise (forward).
//
// Returns a vector of outline indices.
// A transverse divide is always at least 1 pixel edge long, so will have a
// distinct joining and leaving node.
// If the basin is not range-spanning or if it has no transverse divide on this
// side, an emtpy vector is returned.
//
// RV 26/08/2025 
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
std::vector<int> lsdtt_basin_outline_cpu::find_left_transverse_divide_nodes(
    std::int32_t basin_junction,
    const std::vector<std::int32_t>& rsb_N,
    const std::vector<std::int32_t>& rsb_S,
    const lsdtt_junction_network_base& junction_network) const
{
  // Find out which side the current basin is in by comparing basin_junction_
  // number to range_spanning_basins_N and range_spanning_basins_S
  bool basin_is_north = std::find(rsb_N.begin(), rsb_N.end(), basin_junction) != rsb_N.end();
  bool basin_is_south = std::find(rsb_S.begin(), rsb_S.end(), basin_junction) != rsb_S.end();
  if (!basin_is_north && !basin_is_south) return {};

  const std::vector<std::int32_t>& rspb_thisside = basin_is_north ? rsb_N : rsb_S;

  auto get_baselevel = [&](std::int64_t flat) -> std::int64_t {
    std::int64_t curr = flat;
    while (true) {
      std::int64_t rec = flow_->receiver_of_unchecked(curr);
      if (rec == curr || rec < 0 || rec == lsdtt_flow_base::outlet_receiver) return curr;
      curr = rec;
    }
  };

  std::vector<int> nodes_to_return;
  // Iterate over the nodes in the basin outline, forwards
  for (std::size_t i = 0; i < outline_divide_rc_.size(); ++i) {
    std::vector<std::int64_t> px_in_basin, px_out_basin;
    get_pixels_around_basin_outline_node(i, px_in_basin, px_out_basin);

    int N_pixels_in_adjacent_rsb = 0;
    std::vector<std::int64_t> outside_nodes_not_in_adjacent_rsb;
    
    for (std::int64_t node : px_out_basin) {
      std::int64_t bln_index = get_baselevel(node);
      std::int32_t bl_junc = junction_network.junction_of_flat(bln_index);
      if (std::find(rspb_thisside.begin(), rspb_thisside.end(), bl_junc) != rspb_thisside.end()) {
        N_pixels_in_adjacent_rsb++;
      } else {
        outside_nodes_not_in_adjacent_rsb.push_back(node);
      }
    }

    if (nodes_to_return.empty() && N_pixels_in_adjacent_rsb > 0) {
      bool join = true;
      // Ok, need to deal with an edge case:
      // This outline node has a surrounding pixel in a neighbouring range 
      // spanning basin, but it is separated by a diagonal flow edge that 
      // is in neither range-spanning basin. 
      // We can identify such cases as:
      // (1) 3 pixels outside the basin
      // (2) one of which in a different range-spanning basin and
      // (3) the 2 remaining pixels are connected by a 
      //     diagonal flow edge.
      if (px_out_basin.size() == 3 && N_pixels_in_adjacent_rsb == 1 && outside_nodes_not_in_adjacent_rsb.size() == 2) {
        std::int64_t node1 = outside_nodes_not_in_adjacent_rsb[0];
        std::int64_t node2 = outside_nodes_not_in_adjacent_rsb[1];
        std::int64_t node1rec = flow_->receiver_of_unchecked(node1);
        std::int64_t node2rec = flow_->receiver_of_unchecked(node2);
        
        if (node1rec == node2 || node2rec == node1) {
          // Ok, the two nodes outside basin draining same edge
          // are each others receiver.
          // Now check if they are diagonally connected.
          int r1 = flow_->row_from_flat_unchecked(node1);
          int c1 = flow_->col_from_flat_unchecked(node1);
          int r2 = flow_->row_from_flat_unchecked(node2);
          int c2 = flow_->col_from_flat_unchecked(node2);
          if (r1 != r2 && c1 != c2) {
            // Diagonally connected. Set join to false.
            join = false;
          }
        }
      }
      if (join) nodes_to_return.push_back(static_cast<int>(i));
    }

    // Stop if the node has a surrounding pixel outside the current
    // basin that is NOT touching a range-spanning basin
    // In this case, we have just left the Transverse divide
    if (nodes_to_return.size() == 1 && N_pixels_in_adjacent_rsb == 0) {
      nodes_to_return.push_back(static_cast<int>(i));
      return nodes_to_return;
    }
  }


  // If not returned previously, a transverse divide was not found
  // on this side.
  // This can happen due to it being the first basin along a non-periodic
  // boundary, or due to interstitial unchanneled space in-between range-spanning
  // basins (that stretches all the way from mdd to boundary).
  //cout << "Basin " << basin_junction_number << " I found no left tvd!" << endl;
  return nodes_to_return;
}


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// This scans the outline for triple divide junctions where the outline 
// (potentially) joins and leaves the right transverse divide.
// By iterating through the outline counter-clockwise.
//
// Returns a vector of outline indices.
// A transverse divide is always at least 1 pixel edge long, so will have a
// distinct joining and leaving node.
// If the basin is not range-spanning or if it has no transverse divide on this
// side, an emtpy vector is returned.
//
// RV 26/08/2025 
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
std::vector<int> lsdtt_basin_outline_cpu::find_right_transverse_divide_nodes(
    std::int32_t basin_junction,
    const std::vector<std::int32_t>& rsb_N,
    const std::vector<std::int32_t>& rsb_S,
    const lsdtt_junction_network_base& junction_network) const
{
  // Find out which side the current basin is in by comparing basin_junction_
  // number to range_spanning_basins_N and range_spanning_basins_S
  bool basin_is_north = std::find(rsb_N.begin(), rsb_N.end(), basin_junction) != rsb_N.end();
  bool basin_is_south = std::find(rsb_S.begin(), rsb_S.end(), basin_junction) != rsb_S.end();
  if (!basin_is_north && !basin_is_south) return {};

  const std::vector<std::int32_t>& rspb_thisside = basin_is_north ? rsb_N : rsb_S;

  auto get_baselevel = [&](std::int64_t flat) -> std::int64_t {
    std::int64_t curr = flat;
    while (true) {
      std::int64_t rec = flow_->receiver_of_unchecked(curr);
      if (rec == curr || rec < 0 || rec == lsdtt_flow_base::outlet_receiver) return curr;
      curr = rec;
    }
  };

  std::vector<int> nodes_to_return;
  // Iterate over the nodes in the basin outline, backwards
  for (int i = static_cast<int>(outline_divide_rc_.size()) - 1; i >= 0; --i) {
    std::vector<std::int64_t> px_in_basin, px_out_basin;
    get_pixels_around_basin_outline_node(i, px_in_basin, px_out_basin);

    int N_pixels_in_adjacent_rsb = 0;
    std::vector<std::int64_t> outside_nodes_not_in_adjacent_rsb;
    
    for (std::int64_t node : px_out_basin) {
      std::int64_t bln_index = get_baselevel(node);
      std::int32_t bl_junc = junction_network.junction_of_flat(bln_index);
      if (std::find(rspb_thisside.begin(), rspb_thisside.end(), bl_junc) != rspb_thisside.end()) {
        N_pixels_in_adjacent_rsb++;
      } else {
        outside_nodes_not_in_adjacent_rsb.push_back(node);
      }
    }

    if (nodes_to_return.empty() && N_pixels_in_adjacent_rsb > 0) {
      bool join = true;
      // Ok, need to deal with an edge case:
      // This outline node has a surrounding pixel in a neighbouring range 
      // spanning basin, but it is separated by a diagonal flow edge that 
      // is in neither range-spanning basin. 
      // We can identify such cases as:
      // (1) 3 pixels outside the basin
      // (2) one of which in a different range-spanning basin and
      // (3) the 2 remaining pixels are connected by a 
      //     diagonal flow edge.
      if (px_out_basin.size() == 3 && N_pixels_in_adjacent_rsb == 1 && outside_nodes_not_in_adjacent_rsb.size() == 2) {
        std::int64_t node1 = outside_nodes_not_in_adjacent_rsb[0];
        std::int64_t node2 = outside_nodes_not_in_adjacent_rsb[1];
        std::int64_t node1rec = flow_->receiver_of_unchecked(node1);
        std::int64_t node2rec = flow_->receiver_of_unchecked(node2);
        
        if (node1rec == node2 || node2rec == node1) {
          // Ok, the two nodes outside basin draining same edge
          // are each others receiver.
          // Now check if they are diagonally connected.
          int r1 = flow_->row_from_flat_unchecked(node1);
          int c1 = flow_->col_from_flat_unchecked(node1);
          int r2 = flow_->row_from_flat_unchecked(node2);
          int c2 = flow_->col_from_flat_unchecked(node2);
          if (r1 != r2 && c1 != c2) {
            // Diagonally connected. Set join to false.
            join = false;
          }
        }
      }
      if (join) nodes_to_return.push_back(i);
    }

    if (nodes_to_return.size() == 1 && N_pixels_in_adjacent_rsb == 0) {
      nodes_to_return.push_back(i);
      return nodes_to_return;
    }
  }
  // If not returned previously, a transverse divide was not found
  // on this side.
  // This can happen due to it being the first basin along a non-periodic
  // boundary, or due to interstitial unchanneled space in-between range-spanning
  // basins (that stretches all the way from mdd to boundary).
  //cout << "Basin " << basin_junction_number << " I found no left tvd!" << endl;
  return nodes_to_return;
}


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// This scans the outline for the points where the outline 
// (potentially) joins and leaves the main drainage divide.
//
// The start/end nodes are returned in a vector. If the returned vector if empty,
// the basin outline does not touch or stretch along the MDD.
//
// RV 26/08/2025 
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
std::vector<int> lsdtt_basin_outline_cpu::find_main_drainage_divide_nodes(
    std::int32_t basin_junction,
    const std::vector<std::int32_t>& rsb_N,
    const std::vector<std::int32_t>& rsb_S,
    const lsdtt_junction_network_base& junction_network) const
{
  // Find out which side the current basin is in by taking the current basin
  // junction number and tracing it down to the baselevel junction.
  // Then compare this baselevel junction to range_spanning_basins_N and 
  // range_spanning_basins_S.
  auto get_baselevel = [&](std::int64_t flat) -> std::int64_t {
    std::int64_t curr = flat;
    while (true) {
      std::int64_t rec = flow_->receiver_of_unchecked(curr);
      if (rec == curr || rec < 0 || rec == lsdtt_flow_base::outlet_receiver) return curr;
      curr = rec;
    }
  };

  std::int64_t bl_flat = get_baselevel(junction_network.flat_of_junction(basin_junction));
  std::int32_t basin_bl_junction = junction_network.junction_of_flat(bl_flat);

  bool basin_is_north = std::find(rsb_N.begin(), rsb_N.end(), basin_bl_junction) != rsb_N.end();
  bool basin_is_south = std::find(rsb_S.begin(), rsb_S.end(), basin_bl_junction) != rsb_S.end();
  if (!basin_is_north && !basin_is_south) return {};

  const std::vector<std::int32_t>& rspb_thisside = basin_is_north ? rsb_N : rsb_S;
  const std::vector<std::int32_t>& rspb_otherside = basin_is_north ? rsb_S : rsb_N;

  // Make a vector for joining indices and leaving indices
  std::vector<int> join_indices;
  std::vector<int> leave_indices;
  bool on_mdd = false;

  // Iterate over the nodes in the basin outline
  for (std::size_t i = 0; i < outline_divide_rc_.size(); ++i) {
    // Get the DEM pixels outside of the basin surrounding outline node
    std::vector<std::int64_t> px_in_basin, px_out_basin;
    get_pixels_around_basin_outline_node(i, px_in_basin, px_out_basin);


    // Loop through the pixels outside of the basin
    // And get the number of different basins on the same side
    int N_pixels_draining_opposite_side = 0;
    std::vector<std::int64_t> nodes_outside_basin_draining_same_side;
    
    for (std::int64_t node : px_out_basin) {
      // Get the baselevel junction of this neighbouring pixel
      std::int64_t bln_index = get_baselevel(node);
      std::int32_t bl_junc = junction_network.junction_of_flat(bln_index);
       // See if any flow to opposite side
      if (std::find(rspb_otherside.begin(), rspb_otherside.end(), bl_junc) != rspb_otherside.end()) {
        N_pixels_draining_opposite_side++;
      }
      // Identify those drainaing to the same side
      if (std::find(rspb_thisside.begin(), rspb_thisside.end(), bl_junc) != rspb_thisside.end()) {
        nodes_outside_basin_draining_same_side.push_back(node);
      }
    }


    // If we are on the mdd but we find that the 
    // current node does NOT have any surounding
    // pixels in a range-spanning basin on the opposite side, we know
    // that we have left the MDD and found the leaving node (the 
    // "leaving node" is the first node not on the MDD anymore).
    if (on_mdd && N_pixels_draining_opposite_side == 0) {
      int leaving_index = static_cast<int>(i);
      
      // Here we need to treat an edge case. If the node before
      // the leaving node is on the same location as the joining node
      // then we need to mark the leaving node one earlier.
      // This becasue it has essentially "lapped back" to the node
      // where it joins, where a diagonal flow edge connects the basin
      // pixel. But when lapping back (i.e. visiting the second time)
      // you actually have that flow edge separating the divide and 
      // the basin on the other side.
      if (leaving_index - 1 >= 0 && same_rc_(outline_divide_rc_[leaving_index - 1], outline_divide_rc_[join_indices.back()])) {
        if (leaving_index - 1 != join_indices.back()) {
          leaving_index--;
        }
      }
      leave_indices.push_back(leaving_index);
      on_mdd = false;
    }


    
    // If we are not on the mdd and
    // the pixels surrounding the node outside the 
    // current basin is in a range-spanning basin on the opposite side.
    // Then, we have found the first TDJ (we are joining the MDD)
    // Push back this index to the return vector
    if (!on_mdd && N_pixels_draining_opposite_side > 0) {
      int joining_index = static_cast<int>(i);
      bool join = true;
      
      // Ok, need to deal with an edge case:
      // The case where the node has a diagonal flow edge from another
      // basin on the same side separating the node from the basin on the
      // opposite side. We can detect these if there are 
      // (1) 3 pixels outside the basin
      // (2) one of which is on the other side of the MDD and
      // (3) the 2 remaining pixels on the same side are connected by a 
      //     diagonal flow edge.
      if (px_out_basin.size() == 3 && N_pixels_draining_opposite_side == 1 && nodes_outside_basin_draining_same_side.size() == 2) {
        std::int64_t node1 = nodes_outside_basin_draining_same_side[0];
        std::int64_t node2 = nodes_outside_basin_draining_same_side[1];
        std::int64_t node1rec = flow_->receiver_of_unchecked(node1);
        std::int64_t node2rec = flow_->receiver_of_unchecked(node2);
        
        if (node1rec == node2 || node2rec == node1) {
          int r1 = flow_->row_from_flat_unchecked(node1);
          int c1 = flow_->col_from_flat_unchecked(node1);
          int r2 = flow_->row_from_flat_unchecked(node2);
          int c2 = flow_->col_from_flat_unchecked(node2);
          if (r1 != r2 && c1 != c2) join = false;
        }
      }
      if (join) {
        join_indices.push_back(joining_index);
        on_mdd = true;
      }
    }
  }

  if (join_indices.empty()) return {};
  if (leave_indices.empty()) leave_indices.push_back(static_cast<int>(outline_divide_rc_.size()));

  return {join_indices.front(), leave_indices.back()};
}


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// This scans the outline for the points where the outline 
// (potentially) joins and leaves the outline of another, specified basin.
//
// The start/end nodes are returned in a vector. If the returned vector if empty,
// the basin outline does not touch or stretch along the outline of another basin.
//
// RV 26/08/2025 
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
std::vector<int> lsdtt_basin_outline_cpu::find_interfluvial_divide_nodes(
    std::int32_t basin_junction,
    std::int32_t otherBasin_junction_number,
    const lsdtt_junction_network_base& junction_network) const
{
  // Make a vector for joining indices and leaving indices
  std::vector<int> join_indices;
  std::vector<int> leave_indices;
  bool on_interfluvial_divide = false;

  // Check if the otherBasin_junction_number is NOT nested in the current basin junction number, or
  // vice versa
  if (junction_network.is_upslope_of(basin_junction, otherBasin_junction_number) ||
      junction_network.is_upslope_of(otherBasin_junction_number, basin_junction)) {
    return {};
  }

      
  // Iterate over the nodes in the basin outline
  for (std::size_t i = 0; i < outline_divide_rc_.size(); ++i) {
    // Get the DEM pixels surrounding the outline node on the outside of the basin
    std::vector<std::int64_t> px_in_basin, px_out_basin;
    get_pixels_around_basin_outline_node(i, px_in_basin, px_out_basin);

    // Now check if any of these pixels are in the other basin
    // and record the nodes that are in neither the current and the other basin.
    int N_pixels_draining_other_basin = 0;
    std::vector<std::int64_t> nodes_in_neither_basin;

    for (std::int64_t thisNode : px_out_basin) {
      // Get the first channelised node down the d8 flow path from thisNode
      std::int64_t curr = thisNode;
      bool reached_chan = false;
      while (true) {
        if (curr < 0 || curr == lsdtt_flow_base::outlet_receiver) break;
        if (junction_network.is_channel(curr)) {
          reached_chan = true; break;
        }
        std::int64_t rec = flow_->receiver_of_unchecked(curr);
        if (rec == curr) break;
        curr = rec;
      }

      // Check for the case that there is no channelised node
      if (!reached_chan) {
        nodes_in_neither_basin.push_back(thisNode);
        continue;
      }

      // Now get the first upstream junction of that channelised node
      std::int32_t thisNodesJunction = junction_network.junction_upstream_of_flat(curr);

      // And check if this junction is upstream of or equal to the junction number of the other basin            
      if (junction_network.is_upslope_of(thisNodesJunction, otherBasin_junction_number) || thisNodesJunction == otherBasin_junction_number) {
        N_pixels_draining_other_basin++;
      } else {
        nodes_in_neither_basin.push_back(thisNode);
      }
    }

    // If we are on the interfluvial divide between this and the other basin
    // and the current node deos not have any outside pixels in the other 
    // basin, then we know that we have left the interfluvial divide and 
    // record the leaving node. The leaving node is the first node not on the
    // interfluvial divide anymore.
    if (on_interfluvial_divide && N_pixels_draining_other_basin == 0) {
      int leaving_index = static_cast<int>(i);
      
      // Here we need to treat an edge case. If the node before
      // the leaving node is on the same location as the joining node
      // then we need to mark the leaving node one earlier.
      // This becasue it has essentially "lapped back" to the node
      // where it joins, where a diagonal flow edge connects the basin
      // pixel. But when lapping back (i.e. visiting the second time)
      // you actually have that flow edge separating the divide and 
      // the basin on the other side.
      if (leaving_index - 1 >= 0 && same_rc_(outline_divide_rc_[leaving_index - 1], outline_divide_rc_[join_indices.back()])) {
        if (leaving_index - 1 != join_indices.back()) {
          leaving_index--;
        }
      }
      // Push back the leaving index
      leave_indices.push_back(leaving_index);
      on_interfluvial_divide = false;
    }


    // If we are not on the interfluvial divide yet and one or more pixels
    // surrunding the current divide node on the outside is in the other basin, 
    // then we record this as the leaving node by pushing back this index to 
    // the return vector.
    if (!on_interfluvial_divide && N_pixels_draining_other_basin > 0) {
      int joining_index = static_cast<int>(i);
      bool join = true;
      
      // Ok, need to deal with an edge case:
      // The case where the node has a diagonal flow edge from another
      // basin separating the outline node from the other basin. 
      // We can detect these if there are 
      // (1) 3 pixels outside the basin
      // (2) one of which is in the other basin
      // (3) the 2 remaining pixels on the same side are connected by a 
      //     diagonal flow edge.
      if (px_out_basin.size() == 3 && N_pixels_draining_other_basin == 1 && nodes_in_neither_basin.size() == 2) {
        std::int64_t node1 = nodes_in_neither_basin[0];
        std::int64_t node2 = nodes_in_neither_basin[1];
        std::int64_t node1rec = flow_->receiver_of_unchecked(node1);
        std::int64_t node2rec = flow_->receiver_of_unchecked(node2);
        
        if (node1rec == node2 || node2rec == node1) {
          int r1 = flow_->row_from_flat_unchecked(node1);
          int c1 = flow_->col_from_flat_unchecked(node1);
          int r2 = flow_->row_from_flat_unchecked(node2);
          int c2 = flow_->col_from_flat_unchecked(node2);
          if (r1 != r2 && c1 != c2) join = false;
        }
      }
      if (join) {
        join_indices.push_back(joining_index);
        on_interfluvial_divide = true;
      }
    }
  }

  if (join_indices.size() == 1 && leave_indices.empty()) {
    leave_indices.push_back(static_cast<int>(outline_divide_rc_.size()));
  }

  if (join_indices.empty() && leave_indices.empty()) {
    return {};
  }

  return {join_indices.front(), leave_indices.back()};
}

} // namespace lsdtt
