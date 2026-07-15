// lsdtt_raster_vectorise_cpu.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html

#include "lsdtt_raster_vectorise_cpu.hpp"

#include "../core/lsdtt_types.hpp"
#include "../utils/lsdtt_array2d.hpp"
#include "../utils/lsdtt_binary_array2d_tracing.hpp"
#include "../utils/lsdtt_raster_geometry.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <map>
#include <stdexcept>
#include <vector>

namespace lsdtt {

// =============================================================================
// trace_labelled_polygons
// =============================================================================
//
// Algorithm overview
// ------------------
//
// For each unique label L:
//
//   1. Collect all pixel (row, col) pairs that carry label L (mask valid,
//      not nodata, not background).
//
//   2. Find the tight bounding box [min_r, max_r] x [min_c, max_c].
//
//   3. Build a double-resolution binary image of size
//        (2*height + 1) x (2*width + 1)
//      where height = max_r - min_r + 1, width = max_c - min_c + 1.
//
//      Each pixel (r, c) maps to double-res centre (cr, cc):
//        cr = 2*(r - min_r) + 1
//        cc = 2*(c - min_c) + 1
//
//      A 3x3 stamp centred at (cr, cc) fills cells:
//        rows  cr-1 ... cr+1,  cols  cc-1 ... cc+1
//
//      Adjacent pixels that share an edge therefore also share the
//      intermediate double-res cells between them, so the filled binary
//      image has no gaps along shared edges.
//
//   4. Run trace_boundary_pixels() on the double-res image with
//      connectivity = 4.
//
//   5. Discard all double-res boundary pixels that are not at odd row AND
//      odd column; these are pixel-corner positions in the original grid.
//
//   6. For each kept pixel (p.row, p.col) in the double-res grid:
//
//        row_off = (p.row - 1) / 2      (integer, in [0, height])
//        col_off = (p.col - 1) / 2      (integer, in [0, width])
//
//        x = meta.x_min + (min_c + col_off) * cell_size
//        y = meta.y_min + (meta.n_rows - min_r - row_off) * cell_size
//
//      These are pixel-edge (corner) coordinates in the projected CRS -
//      not cell-centre coordinates.
//
//   7. Close the ring (append the first point) and write.
//
// Coordinate derivation
// ---------------------
//   meta.x_min / meta.y_min are the lower-left raster *edges* (not centres).
//   Cell centre at (row, col):
//     x_ctr = x_min + (col + 0.5) * cs
//     y_ctr = y_min + (n_rows - row - 0.5) * cs
//
//   Pixel corner at integer (row_off, col_off) offset from (min_r, min_c):
//     x_corner = x_min + (min_c + col_off) * cs
//     y_corner = y_min + (n_rows - min_r - row_off) * cs
//
//   row_off = 0 -> top edge of pixel min_r
//   row_off = height -> bottom edge of pixel max_r
//   col_off = 0 -> left edge of pixel min_c
//   col_off = width -> right edge of pixel max_c
//
// =============================================================================

void trace_labelled_polygons(
    const lsdtt_raster_cpu&          labelled,
    lsdtt_data_loader&               loader,
    const std::string&               out_path,
    const lsdtt_vectorise_options&   options,
    const std::vector<std::int32_t>& label_source_id
)
{
    const lsdtt_raster_meta& meta = labelled.get_meta();
    if (!meta.is_valid())
        throw std::runtime_error(
            "lsdtt::trace_labelled_polygons: invalid raster metadata");
    meta.check_projected_metres();

    const int   nrows    = meta.n_rows;
    const int   ncols    = meta.n_cols;
    const double cs      = meta.cell_size;
    const float  bg      = options.background;

    const std::vector<float>&   src = labelled.get_data();
    const std::vector<uint8_t>& msk = labelled.get_mask();

    // -------------------------------------------------------------------------
    // Pass 1 - collect pixel positions per label.
    //
    // Key: int64 label value.
    // Value: vector of (row, col) pairs encoded as int (row * ncols + col
    //        is avoided intentionally so we can slice by row/col directly).
    // -------------------------------------------------------------------------
    struct RcPair { int r; int c; };
    std::map<std::int64_t, std::vector<RcPair>> label_pixels;

    for (int r = 0; r < nrows; ++r)
    {
        for (int c = 0; c < ncols; ++c)
        {
            const std::size_t idx = static_cast<std::size_t>(r) * ncols + c;
            if (msk[idx] == 0) continue;
            const float v = src[idx];
            if (meta.is_nodata(v)) continue;
            if (v == bg) continue;

            const auto label = static_cast<std::int64_t>(std::round(
                static_cast<double>(v)));
            label_pixels[label].push_back({r, c});
        }
    }

    if (label_pixels.empty()) return;   // nothing to trace

    // -------------------------------------------------------------------------
    // Open vector output
    // -------------------------------------------------------------------------
    const bool has_source_id = !label_source_id.empty();

    lsdtt_vector_schema schema;
    schema.field_names = { "label_id", "pixel_count", "area_m2" };
    schema.field_types = {
        lsdtt_vector_field_type::int64_type,
        lsdtt_vector_field_type::int64_type,
        lsdtt_vector_field_type::double_type
    };
    if (has_source_id)
    {
        schema.field_names.push_back("source_id");
        schema.field_types.push_back(lsdtt_vector_field_type::int64_type);
    }

    loader.begin_vector_write(
        out_path,
        options.layer_name,
        schema,
        meta.crs_wkt,
        options.format,
        options.spatial_index
    );

    // -------------------------------------------------------------------------
    // Pass 2 - per-label: build double-res binary, trace, convert, write.
    // -------------------------------------------------------------------------
    for (const auto& [label, pixels] : label_pixels)
    {
        // --- Bounding box ---
        int min_r = pixels[0].r, max_r = pixels[0].r;
        int min_c = pixels[0].c, max_c = pixels[0].c;
        for (const auto& p : pixels)
        {
            if (p.r < min_r) min_r = p.r;
            if (p.r > max_r) max_r = p.r;
            if (p.c < min_c) min_c = p.c;
            if (p.c > max_c) max_c = p.c;
        }

        const int height = max_r - min_r + 1;
        const int width  = max_c - min_c + 1;

        // Double-res grid: (2*height + 1) rows x (2*width + 1) cols.
        // Row i in double-res spans original pixels at row (i-1)/2 + min_r
        // (for odd i) or the edge between rows i/2-1 and i/2 (for even i).
        const int dr_rows = 2 * height + 1;
        const int dr_cols = 2 * width  + 1;

        lsdtt_array2d_u8 binary(dr_rows, dr_cols, uint8_t{0});

        // --- Stamp each pixel ---
        // 3x3 centred at (cr, cc) = (2*(r-min_r)+1, 2*(c-min_c)+1).
        static constexpr int stamp[3][3][2] = {
            {{-1,-1},{-1, 0},{-1, 1}},
            {{ 0,-1},{ 0, 0},{ 0, 1}},
            {{ 1,-1},{ 1, 0},{ 1, 1}}
        };

        for (const auto& p : pixels)
        {
            const int cr = 2 * (p.r - min_r) + 1;
            const int cc = 2 * (p.c - min_c) + 1;
            for (int si = 0; si < 3; ++si)
                for (int sj = 0; sj < 3; ++sj)
                    binary(cr + stamp[si][sj][0],
                           cc + stamp[si][sj][1]) = uint8_t{1};
        }

        // --- Trace boundary ---
        const std::vector<lsdtt_pixel_rc> boundary =
            trace_boundary_pixels(binary, 4);

        if (boundary.empty()) continue;

        // --- Keep only odd-row x odd-col pixels (original pixel corners) ---
        std::vector<lsdtt_xy> ring;
        ring.reserve(boundary.size() / 2 + 1);

        for (const lsdtt_pixel_rc& bp : boundary)
        {
            if ((bp.row & 1) == 0 || (bp.col & 1) == 0) continue;

            // Map double-res corner to projected coordinates.
            // row_off and col_off are integer offsets from (min_r, min_c)
            // in original pixel units.
            const int row_off = (bp.row - 1) / 2;  // in [0, height]
            const int col_off = (bp.col - 1) / 2;  // in [0, width]

            const double px =
                meta.x_min + static_cast<double>(min_c + col_off) * cs;
            const double py =
                meta.y_min +
                static_cast<double>(nrows - min_r - row_off) * cs;

            ring.push_back({px, py});
        }

        if (ring.size() < 3) continue;   // degenerate - skip

        // Close the ring
        ring.push_back(ring.front());

        // --- Write feature ---
        const auto pixel_count = static_cast<std::int64_t>(pixels.size());
        const double area_m2   = static_cast<double>(pixel_count) * cs * cs;

        std::vector<lsdtt_vector_field_value> attrs = { label, pixel_count, area_m2 };
        if (has_source_id)
        {
            const std::size_t lbl_idx = static_cast<std::size_t>(label);
            const std::int64_t src_id =
                (lbl_idx < label_source_id.size())
                ? static_cast<std::int64_t>(label_source_id[lbl_idx])
                : std::int64_t{-1};
            attrs.push_back(src_id);
        }

        loader.write_polygon(ring, attrs);
    }

    loader.finish_vector_write();
}

} // namespace lsdtt
