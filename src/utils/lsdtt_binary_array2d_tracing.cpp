// lsdtt_binary_array2d_tracing.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// See lsdtt_binary_array2d_tracing.hpp for documentation.

#include "lsdtt_binary_array2d_tracing.hpp"

#include <array>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <vector>

namespace lsdtt
{
namespace
{

struct tracing_tables
{
    // Only the first n_neighbours entries are used.
    std::array<int, 8> offsets{};
    std::array<int, 8> next_search_direction_table{};
    std::array<int, 8> next_direction_table{};
    int n_neighbours = 0;
};

tracing_tables make_tables(int connectivity, int n_cols_padded)
{
    tracing_tables t;

    if (connectivity == 8)
    {
        t.n_neighbours = 8;

        // Neighbour order: N, NE, E, SE, S, SW, W, NW
        // Row-major linear index: idx = r*n_cols + c
        t.offsets = {
            -n_cols_padded,         // N
            -n_cols_padded + 1,     // NE
            1,                      // E
            n_cols_padded + 1,      // SE
            n_cols_padded,          // S
            n_cols_padded - 1,      // SW
            -1,                     // W
            -n_cols_padded - 1      // NW
        };

        // If you moved dir, start next search from:
        // N/NE -> NW, E/SE -> NE, S/SW -> SE, W/NW -> SW
        t.next_search_direction_table = {7, 7, 1, 1, 3, 3, 5, 5};

        // Clockwise neighbour search progression
        t.next_direction_table = {1, 2, 3, 4, 5, 6, 7, 0};
    }
    else
    {
        t.n_neighbours = 4;

        // Neighbour order: N, E, S, W
        t.offsets = {
            -n_cols_padded, // N
            1,              // E
            n_cols_padded,  // S
            -1,             // W
            0, 0, 0, 0
        };

        // Start direction for next search based on previous movement direction.
        // (Matches the LSDTopoTools2 table you pasted.)
        t.next_search_direction_table = {3, 0, 1, 2, 0, 0, 0, 0};

        // Clockwise neighbour search progression: N -> E -> S -> W -> N
        t.next_direction_table = {1, 2, 3, 0, 0, 0, 0, 0};
    }

    return t;
}

} // anonymous namespace

std::vector<lsdtt_pixel_rc>
trace_boundary_pixels(
    const lsdtt_array2d_u8& bw,
    int connectivity
)
{
    if (connectivity != 8 && connectivity != 4){
        throw std::runtime_error(
            "lsdtt_binary_array2d_tracing: connectivity must be 4 or 8!"
        );
    }
    
    if (!bw.is_valid())
    {
        return {};
    }

    const int bw_n_rows = bw.n_rows();
    const int bw_n_cols = bw.n_cols();

    if (bw_n_rows <= 0 || bw_n_cols <= 0)
    {
        return {};
    }

    // -------------------------------------------------------------------------
    // Pad the binary image with a background border (zeros).
    // Use int so we can write sentinel values (-1, -2) exactly as in LSDT2.
    // -------------------------------------------------------------------------
    const int padded_n_rows = bw_n_rows + 2;
    const int padded_n_cols = bw_n_cols + 2;

    lsdtt_array2d<int> L(padded_n_rows, padded_n_cols, 0);

    for (int r = 0; r < bw_n_rows; ++r)
    {
        for (int c = 0; c < bw_n_cols; ++c)
        {
            L(r + 1, c + 1) = (bw(r, c) != 0u) ? 1 : 0;
        }
    }

    const tracing_tables tables = make_tables(connectivity, padded_n_cols);

    // Sentinel markers in the padded array (as per LSDT2 implementation)
    constexpr int start_pixel_value    = -1;
    constexpr int boundary_pixel_value = -2;

    // -------------------------------------------------------------------------
    // Find a starting pixel: first foreground pixel (row-major scan)
    // such that the pixel above it is background.
    // -------------------------------------------------------------------------
    int  start_row = 0;
    int  start_col = 0;
    bool found_start = false;

    for (int r = 1; r < padded_n_rows && !found_start; ++r)
    {
        for (int c = 0; c < padded_n_cols && !found_start; ++c)
        {
            if (L(r, c) == 1 && L(r - 1, c) == 0)
            {
                start_row = r;
                start_col = c;
                found_start = true;
            }
        }
    }

    if (!found_start)
    {
        return {};
    }

    L(start_row, start_col) = start_pixel_value;

    // -------------------------------------------------------------------------
    // Core Moore neighbourhood tracing loop
    // -------------------------------------------------------------------------
    std::vector<int> boundary_indices;
    boundary_indices.reserve(static_cast<std::size_t>(bw_n_rows + bw_n_cols) * 2);

    int current_pixel = start_row * padded_n_cols + start_col;
    boundary_indices.push_back(current_pixel);

    int next_search_direction = 1;

    bool done = false;

    // Finishing condition added in LSDT2 version you pasted:
    int first_boundary_pixel = -1;
    int first_direction      = -1;

    while (!done)
    {
        bool found_next_pixel = false;

        int direction = next_search_direction;

        for (int k = 0; k < tables.n_neighbours; ++k)
        {
            const int neighbour_pixel = current_pixel + tables.offsets[direction];

            const int nr = neighbour_pixel / padded_n_cols;
            const int nc = neighbour_pixel % padded_n_cols;

            // neighbour_pixel should always be in bounds due to padding,
            // but guard anyway in case of future modifications.
            if (nr < 0 || nr >= padded_n_rows || nc < 0 || nc >= padded_n_cols)
            {
                direction = tables.next_direction_table[direction];
                continue;
            }

            const int L_value = L(nr, nc);

            if (L_value != 0)
            {
                // Finishing condition: detect when we are about to retrace.
                if (first_boundary_pixel == -1 && L_value != start_pixel_value)
                {
                    first_boundary_pixel = neighbour_pixel;
                    first_direction      = direction;
                }
                else if (neighbour_pixel == first_boundary_pixel &&
                         direction == first_direction)
                {
                    done = true;
                    break;
                }

                // Take the step.
                next_search_direction =
                    tables.next_search_direction_table[direction];

                found_next_pixel = true;

                // Mark boundary pixel unless it is the start pixel.
                if (L_value != start_pixel_value)
                {
                    L(nr, nc) = boundary_pixel_value;
                }

                boundary_indices.push_back(neighbour_pixel);
                current_pixel = neighbour_pixel;
                break;
            }

            direction = tables.next_direction_table[direction];
        }

        if (!found_next_pixel && !done)
        {
            // No neighbour found: degenerate single-pixel object.
            boundary_indices.push_back(current_pixel);
            done = true;
        }
    }

    // -------------------------------------------------------------------------
    // Convert padded linear indices back to row/col in the original bw grid.
    // padded coords (r, c) map to bw coords (r-1, c-1).
    // -------------------------------------------------------------------------
    std::vector<lsdtt_pixel_rc> out;
    out.reserve(boundary_indices.size());

    for (const int idx : boundary_indices)
    {
        const int rL = idx / padded_n_cols;
        const int cL = idx % padded_n_cols;

        out.push_back(lsdtt_pixel_rc{rL - 1, cL - 1});
    }

    return out;
}

} // namespace lsdtt
