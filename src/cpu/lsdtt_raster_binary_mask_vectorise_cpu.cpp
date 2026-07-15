// lsdtt_raster_binary_mask_vectorise_cpu.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html


#include <algorithm>
#include <limits>
#include <map>
#include <stdexcept>
#include <utility>
#include <vector>


#include "lsdtt_raster_binary_mask_vectorise_cpu.hpp"
#include "../utils/lsdtt_binary_array2d_tracing.hpp"


namespace lsdtt
{
namespace
{


// -----------------------------------------------------------------------------
// filter_odd_odd_to_xy
//
// Shared helper: given a boundary vector in double-resolution grid coordinates
// (local, offset by box_min_r / box_min_c), filter to odd/odd corners,
// map back to divide-grid row/col, convert to CRS xy, close the ring,
// and remove consecutive duplicates.
// -----------------------------------------------------------------------------
static std::vector<lsdtt_xy> filter_odd_odd_to_xy(
    const std::vector<lsdtt_pixel_rc>& boundary,
    int                                box_min_r,
    int                                box_min_c,
    const lsdtt_raster_meta&           div_meta)
{
    std::vector<lsdtt_xy> ring;
    ring.reserve(boundary.size());

    for (const auto& p : boundary)
    {
        if (!((p.row & 1) && (p.col & 1))) continue;

        const int div_r = ((p.row + 1) >> 1) + box_min_r;
        const int div_c = ((p.col + 1) >> 1) + box_min_c;

        ring.push_back(cell_center_from_row_col(div_meta, div_r, div_c));
    }

    // Remove consecutive duplicates first
    ring.erase(
        std::unique(ring.begin(), ring.end(),
                    [](const lsdtt_xy& a, const lsdtt_xy& b)
                    { return a.x == b.x && a.y == b.y; }),
        ring.end());

    // Close the ring
    if (ring.size() >= 2)
    {
        const auto& f = ring.front();
        const auto& l = ring.back();
        if (f.x != l.x || f.y != l.y)
            ring.push_back(f);
    }

    return ring;
}



// -----------------------------------------------------------------------------
// find_hole_rings
// -----------------------------------------------------------------------------
static std::vector<std::vector<lsdtt_xy>> find_hole_rings(
    const std::vector<std::size_t>& component_pixels, // pixel_indices for this component
    int                             nrows_raster,
    int                             ncols_raster,
    int                             box_min_r,
    int                             box_max_r,
    int                             box_min_c,
    int                             box_max_c,
    const lsdtt_raster_meta&        div_meta,
    double                          cell_size)
{
    // Work on the original-resolution pixel grid, cropped to the component
    // bounding box (already expanded by 1 on each side).
    const int LR = box_max_r - box_min_r + 1;
    const int LC = box_max_c - box_min_c + 1;

    // Mark component pixels as foreground in a local grid.
    std::vector<uint8_t> local(static_cast<std::size_t>(LR) *
                               static_cast<std::size_t>(LC), 0u);

    auto lidx = [&](int r, int c) -> std::size_t
    {
        return static_cast<std::size_t>(r - box_min_r) *
               static_cast<std::size_t>(LC) +
               static_cast<std::size_t>(c - box_min_c);
    };

    for (const std::size_t px : component_pixels)
    {
        const int row = static_cast<int>(px / static_cast<std::size_t>(ncols_raster));
        const int col = static_cast<int>(px % static_cast<std::size_t>(ncols_raster));
        local[lidx(row, col)] = 1u;
    }

    // Flood-fill exterior background from the border of the local grid (4-conn).
    // Because box is expanded by 1 outside the component, the border is always
    // background - no foreground pixel can sit on it.
    std::vector<uint8_t> visited(static_cast<std::size_t>(LR) *
                                 static_cast<std::size_t>(LC), 0u);

    std::vector<std::pair<int,int>> queue;
    queue.reserve(static_cast<std::size_t>(2 * (LR + LC)));

    auto try_enqueue = [&](int r, int c)
    {
        const std::size_t i = lidx(r, c);
        if (local[i] != 0 || visited[i] != 0) return;
        visited[i] = 1u;
        queue.push_back({r, c});
    };

    for (int c = box_min_c; c <= box_max_c; ++c)
    {
        try_enqueue(box_min_r, c);
        try_enqueue(box_max_r, c);
    }
    for (int r = box_min_r + 1; r < box_max_r; ++r)
    {
        try_enqueue(r, box_min_c);
        try_enqueue(r, box_max_c);
    }

    const int dr4[4] = {-1, 1,  0, 0};
    const int dc4[4] = { 0, 0, -1, 1};

    for (std::size_t qi = 0; qi < queue.size(); ++qi)
    {
        const auto [r, c] = queue[qi];
        for (int d = 0; d < 4; ++d)
        {
            const int nr = r + dr4[d];
            const int nc = c + dc4[d];
            if (nr < box_min_r || nr > box_max_r ||
                nc < box_min_c || nc > box_max_c) continue;
            try_enqueue(nr, nc);
        }
    }

    // Any unvisited background pixel is an interior hole.
    // Collect each hole island, stamp it into a double-resolution binary,
    // trace it, and filter to odd/odd corners exactly as the outer ring.
    static constexpr int stamp9[9][2] = {
        {-1,-1},{-1, 0},{-1, 1},
        { 0,-1},{ 0, 0},{ 0, 1},
        { 1,-1},{ 1, 0},{ 1, 1}
    };

    std::vector<std::vector<lsdtt_xy>> hole_rings;

    for (int r0 = box_min_r; r0 <= box_max_r; ++r0)
    {
        for (int c0 = box_min_c; c0 <= box_max_c; ++c0)
        {
            const std::size_t i0 = lidx(r0, c0);
            if (local[i0] != 0 || visited[i0] != 0) continue;

            // Collect this hole island by flood-fill (4-conn)
            std::vector<std::pair<int,int>> island;
            std::vector<std::pair<int,int>> stack;
            stack.push_back({r0, c0});
            visited[i0] = 1u;

            while (!stack.empty())
            {
                const auto [r, c] = stack.back();
                stack.pop_back();
                island.push_back({r, c});

                for (int d = 0; d < 4; ++d)
                {
                    const int nr = r + dr4[d];
                    const int nc = c + dc4[d];
                    if (nr < box_min_r || nr > box_max_r ||
                        nc < box_min_c || nc > box_max_c) continue;
                    const std::size_t ni = lidx(nr, nc);
                    if (local[ni] != 0 || visited[ni] != 0) continue;
                    visited[ni] = 1u;
                    stack.push_back({nr, nc});
                }
            }

            // Stamp island into its own double-resolution binary and trace -
            // identical pipeline to the outer ring in binary_components_to_polygon_records.
            // Bounding box of this island:
            int h_min_r = std::numeric_limits<int>::max();
            int h_max_r = std::numeric_limits<int>::min();
            int h_min_c = std::numeric_limits<int>::max();
            int h_max_c = std::numeric_limits<int>::min();
            for (const auto& [r, c] : island)
            {
                h_min_r = std::min(h_min_r, r);
                h_max_r = std::max(h_max_r, r);
                h_min_c = std::min(h_min_c, c);
                h_max_c = std::max(h_max_c, c);
            }
            // Expand by 1 so the tracer's internal border pad works correctly
            h_min_r -= 1; h_max_r += 1;
            h_min_c -= 1; h_max_c += 1;

            const int HNR = (h_max_r - h_min_r + 1) * 2 + 1;
            const int HNC = (h_max_c - h_min_c + 1) * 2 + 1;

            lsdtt_array2d_u8 hole_bin(HNR, HNC, 0u);
            for (const auto& [r, c] : island)
            {
                const int base_r = (r - h_min_r) * 2;
                const int base_c = (c - h_min_c) * 2;
                for (int k = 0; k < 9; ++k)
                    hole_bin(base_r + stamp9[k][0], base_c + stamp9[k][1]) = 1u;
            }

            const std::vector<lsdtt_pixel_rc> hb =
                trace_boundary_pixels(hole_bin, 4);

            std::vector<lsdtt_xy> ring =
                filter_odd_odd_to_xy(hb, h_min_r, h_min_c, div_meta);

            if (ring.size() < 4) continue;

            // Reverse for correct winding: see comment in previous version.
            std::reverse(ring.begin(), ring.end() - 1);
            ring.back() = ring.front();

            hole_rings.push_back(std::move(ring));
        }
    }

    return hole_rings;
}
} // anonymous namespace




// -----------------------------------------------------------------------------
// binary_components_to_polygon_records
// -----------------------------------------------------------------------------
std::vector<lsdtt_vector_record>
binary_components_to_polygon_records(
    const lsdtt_raster_cpu&    label_raster,
    const lsdtt_raster_cpu*    scalar_raster,
    const lsdtt_vector_schema& schema)
{
    const lsdtt_raster_meta& meta  = label_raster.get_meta();
    const int nrows                = meta.n_rows;
    const int ncols                = meta.n_cols;
    const double cell_size         = meta.cell_size;
    const double cell_area         = cell_size * cell_size;

    const std::vector<float>&   lab_data = label_raster.get_data();
    const std::vector<uint8_t>& lab_mask = label_raster.get_mask();

    const std::size_t npix =
        static_cast<std::size_t>(nrows) * static_cast<std::size_t>(ncols);

    // Optional scalar data
    const std::vector<float>*   sc_data = nullptr;
    const std::vector<uint8_t>* sc_mask = nullptr;
    if (scalar_raster != nullptr)
    {
        sc_data = &scalar_raster->get_data();
        sc_mask = &scalar_raster->get_mask();
    }

    // -------------------------------------------------------------------------
    // Pass 1 - collect pixel indices per label + bounding boxes + scalar sums
    // -------------------------------------------------------------------------
    struct ComponentInfo
    {
        std::vector<std::size_t> pixel_indices;
        int min_r = std::numeric_limits<int>::max();
        int max_r = std::numeric_limits<int>::min();
        int min_c = std::numeric_limits<int>::max();
        int max_c = std::numeric_limits<int>::min();
        double scalar_sum = 0.0;
    };

    std::map<int, ComponentInfo> components;

    for (std::size_t i = 0; i < npix; ++i)
    {
        if (lab_mask[i] == 0) continue;
        const int label = static_cast<int>(lab_data[i]);
        if (label <= 0) continue;  // 0 = background / nodata

        const int row = static_cast<int>(i / static_cast<std::size_t>(ncols));
        const int col = static_cast<int>(i % static_cast<std::size_t>(ncols));

        auto& ci = components[label];
        ci.pixel_indices.push_back(i);
        ci.min_r = std::min(ci.min_r, row);
        ci.max_r = std::max(ci.max_r, row);
        ci.min_c = std::min(ci.min_c, col);
        ci.max_c = std::max(ci.max_c, col);

        if (sc_data != nullptr && (*sc_mask)[i] != 0)
        {
            const float sv = (*sc_data)[i];
            if (std::isfinite(sv) && !meta.is_nodata(sv))
                ci.scalar_sum += static_cast<double>(sv);
        }
    }

    // -------------------------------------------------------------------------
    // Divide-grid meta: same construction as lsdtt_basin_outline_cpu.
    // Shifts origin by -0.5*cell_size and adds 1 row/col (non-periodic).
    // -------------------------------------------------------------------------
    lsdtt_raster_meta div_meta = meta;
    div_meta.x_min   = meta.x_min - 0.5 * cell_size;
    div_meta.y_min   = meta.y_min - 0.5 * cell_size;
    div_meta.n_rows  = nrows + 1;
    div_meta.n_cols  = ncols + 1;

    // -------------------------------------------------------------------------
    // Pass 2 - per component: stamp double-resolution grid, trace, convert
    // -------------------------------------------------------------------------
    // 3x3 stamp offsets - identical to lsdtt_basin_outline_cpu::stamp9
    static constexpr int stamp9[9][2] = {
        {-1,-1},{-1, 0},{-1, 1},
        { 0,-1},{ 0, 0},{ 0, 1},
        { 1,-1},{ 1, 0},{ 1, 1}
    };

    std::vector<lsdtt_vector_record> records;
    records.reserve(components.size());

    for (auto& [label, ci] : components)
    {
        // Expand bounding box by 1 on each side (same as basin outline)
        const int box_min_r = ci.min_r - 1;
        const int box_max_r = ci.max_r + 1;
        const int box_min_c = ci.min_c - 1;
        const int box_max_c = ci.max_c + 1;

        const int nrows_double = (box_max_r - box_min_r + 1) * 2 + 1;
        const int ncols_double = (box_max_c - box_min_c + 1) * 2 + 1;

        lsdtt_array2d_u8 binary(nrows_double, ncols_double, 0u);

        for (const std::size_t idx : ci.pixel_indices)
        {
            const int row = static_cast<int>(idx / static_cast<std::size_t>(ncols));
            const int col = static_cast<int>(idx % static_cast<std::size_t>(ncols));

            const int base_r = (row - box_min_r) * 2;
            const int base_c = (col - box_min_c) * 2;

            for (int k = 0; k < 9; ++k)
            {
                const int rr = base_r + stamp9[k][0];
                const int cc = base_c + stamp9[k][1];
                binary(rr, cc) = 1u;
            }
        }

        // Trace boundary (4-connectivity) - same call as basin outline
        const std::vector<lsdtt_pixel_rc> boundary =
            trace_boundary_pixels(binary, 4);

        // Keep only odd/odd corners and map back to divide-grid coordinates
        // (identical to the filter in edge_tracing_normal_boundaries_)
        lsdtt_polygon_geom poly;

        for (const auto& p : boundary)
        {
            if (!((p.row & 1) && (p.col & 1))) continue;

            const int div_r = ((p.row + 1) >> 1) + box_min_r;
            const int div_c = ((p.col + 1) >> 1) + box_min_c;

            const lsdtt_xy xy =
                cell_center_from_row_col(div_meta, div_r, div_c);
            poly.exterior.push_back(xy);
        }

        // Close the ring if the tracer did not return to the start
        if (poly.exterior.size() >= 2)
        {
            const auto& first = poly.exterior.front();
            const auto& last  = poly.exterior.back();
            if (first.x != last.x || first.y != last.y)
                poly.exterior.push_back(first);
        }

        // Build record
        lsdtt_vector_record rec;
        rec.geometry = std::move(poly);

        const double area_m2 =
            static_cast<double>(ci.pixel_indices.size()) * cell_area;

        if (sc_data != nullptr)
        {
            // schema: label_id, area_m2, scalar_sum
            rec.field_values = {
                static_cast<std::int64_t>(label),
                area_m2,
                ci.scalar_sum * cell_area   // e.g. m^3 when scalar is depth in metres
            };
        }
        else
        {
            // schema: label_id, area_m2
            rec.field_values = {
                static_cast<std::int64_t>(label),
                area_m2
            };
        }

        records.push_back(std::move(rec));
    }

    return records;
}



std::vector<lsdtt_vector_record>
binary_components_to_polygon_records_with_holes(
    const lsdtt_raster_cpu&    label_raster,
    const lsdtt_raster_cpu*    scalar_raster,
    const lsdtt_vector_schema& schema)
{
    label_raster.get_meta().check_projected_metres();

    const lsdtt_raster_meta& meta  = label_raster.get_meta();
    const int nrows                = meta.n_rows;
    const int ncols                = meta.n_cols;
    const double cell_size         = meta.cell_size;
    const double cell_area         = cell_size * cell_size;

    const std::vector<float>&   lab_data = label_raster.get_data();
    const std::vector<uint8_t>& lab_mask = label_raster.get_mask();

    const std::size_t npix =
        static_cast<std::size_t>(nrows) * static_cast<std::size_t>(ncols);

    const std::vector<float>*   sc_data = nullptr;
    const std::vector<uint8_t>* sc_mask = nullptr;
    if (scalar_raster != nullptr)
    {
        sc_data = &scalar_raster->get_data();
        sc_mask = &scalar_raster->get_mask();
    }

    // -------------------------------------------------------------------------
    // Pass 1 - collect pixel indices per label, bounding boxes, scalar sums
    // (unchanged from previous version)
    // -------------------------------------------------------------------------
    struct ComponentInfo
    {
        std::vector<std::size_t> pixel_indices;
        int min_r = std::numeric_limits<int>::max();
        int max_r = std::numeric_limits<int>::min();
        int min_c = std::numeric_limits<int>::max();
        int max_c = std::numeric_limits<int>::min();
        double scalar_sum = 0.0;
    };

    std::map<int, ComponentInfo> components;

    for (std::size_t i = 0; i < npix; ++i)
    {
        if (lab_mask[i] == 0) continue;
        const int label = static_cast<int>(lab_data[i]);
        if (label <= 0) continue;

        const int row = static_cast<int>(i / static_cast<std::size_t>(ncols));
        const int col = static_cast<int>(i % static_cast<std::size_t>(ncols));

        auto& ci = components[label];
        ci.pixel_indices.push_back(i);
        ci.min_r = std::min(ci.min_r, row);
        ci.max_r = std::max(ci.max_r, row);
        ci.min_c = std::min(ci.min_c, col);
        ci.max_c = std::max(ci.max_c, col);

        if (sc_data != nullptr && (*sc_mask)[i] != 0)
        {
            const float sv = (*sc_data)[i];
            if (std::isfinite(sv) && !meta.is_nodata(sv))
                ci.scalar_sum += static_cast<double>(sv);
        }
    }

    // Divide-grid meta (non-periodic: shift origin, +1 row/col)
    lsdtt_raster_meta div_meta = meta;
    div_meta.x_min  = meta.x_min - 0.5 * cell_size;
    div_meta.y_min  = meta.y_min - 0.5 * cell_size;
    div_meta.n_rows = nrows + 1;
    div_meta.n_cols = ncols + 1;

    static constexpr int stamp9[9][2] = {
        {-1,-1},{-1, 0},{-1, 1},
        { 0,-1},{ 0, 0},{ 0, 1},
        { 1,-1},{ 1, 0},{ 1, 1}
    };

    // -------------------------------------------------------------------------
    // Pass 2 - per component: stamp, trace outer ring, find holes
    // -------------------------------------------------------------------------
    std::vector<lsdtt_vector_record> records;
    records.reserve(components.size());

    for (auto& [label, ci] : components)
    {
        const int box_min_r = ci.min_r - 1;
        const int box_max_r = ci.max_r + 1;
        const int box_min_c = ci.min_c - 1;
        const int box_max_c = ci.max_c + 1;

        const int NR = (box_max_r - box_min_r + 1) * 2 + 1;
        const int NC = (box_max_c - box_min_c + 1) * 2 + 1;

        lsdtt_array2d_u8 binary(NR, NC, 0u);

        for (const std::size_t px : ci.pixel_indices)
        {
            const int row = static_cast<int>(px / static_cast<std::size_t>(ncols));
            const int col = static_cast<int>(px % static_cast<std::size_t>(ncols));
            const int base_r = (row - box_min_r) * 2;
            const int base_c = (col - box_min_c) * 2;
            for (int k = 0; k < 9; ++k)
                binary(base_r + stamp9[k][0], base_c + stamp9[k][1]) = 1u;
        }

        // --- Outer ring ---
        const std::vector<lsdtt_pixel_rc> outer_boundary =
            trace_boundary_pixels(binary, 4);

        std::vector<lsdtt_xy> exterior =
            filter_odd_odd_to_xy(outer_boundary, box_min_r, box_min_c,
                                 div_meta);

        if (exterior.size() < 4)
            continue;  // degenerate component - skip

        // --- Hole rings ---
        std::vector<std::vector<lsdtt_xy>> holes =
            find_hole_rings(
                ci.pixel_indices,
                nrows, ncols,
                box_min_r, box_max_r,
                box_min_c, box_max_c,
                div_meta,
                cell_size
            );

        // --- Assemble record ---
        lsdtt_polygon_geom poly;
        poly.exterior = std::move(exterior);
        poly.holes    = std::move(holes);

        lsdtt_vector_record rec;
        rec.geometry = std::move(poly);

        const double area_m2 =
            static_cast<double>(ci.pixel_indices.size()) * cell_area;

        if (sc_data != nullptr)
        {
            rec.field_values = {
                static_cast<std::int64_t>(label),
                area_m2,
                ci.scalar_sum * cell_area
            };
        }
        else
        {
            rec.field_values = {
                static_cast<std::int64_t>(label),
                area_m2
            };
        }

        records.push_back(std::move(rec));
    }

    return records;
}




} // namespace lsdtt