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
      // flow → E : corners SE then NE
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
      // flow → W : corners NW then SW
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
      // flow → S : corners SW then SE
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
      // flow → N : corners NE then NW
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

    // Close the ring — the outline starts at the outlet cut and ends there
    // too after rotate_open_outline_to_outlet_cut_, but the closing point
    // is not repeated in outline_divide_rc_ itself.
    poly.exterior.push_back(poly.exterior.front());

    rec.geometry = std::move(poly);
    return rec;
}


} // namespace lsdtt
