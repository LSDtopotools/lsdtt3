// lsdtt_raster_nodata_filling.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2025 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html

#include "lsdtt_raster_nodata_filling_cpu.hpp"
#include "../cpu/lsdtt_raster_cpu.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <queue>
#include <stdexcept>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace lsdtt {

// =============================================================================
// Step 1 — interior nodata classification
// =============================================================================

// 8-neighbour offsets (used throughout this file)
static const int k_nd_dr[8] = { -1, -1,  0,  1,  1,  1,  0, -1 };
static const int k_nd_dc[8] = {  0,  1,  1,  1,  0, -1, -1, -1 };

lsdtt_array2d<uint8_t> identify_interior_nodata(
    const lsdtt_raster_cpu&raster)
{
    const lsdtt_raster_meta& meta = raster.get_meta();
    if (!meta.is_valid())
    {
        throw std::runtime_error(
            "lsdtt::identify_interior_nodata: invalid raster metadata"
        );
    }

    const int nrows = meta.n_rows;
    const int ncols = meta.n_cols;

    const lsdtt_array2d<uint8_t> mask(nrows, ncols,
                                      lsdtt_array2d<uint8_t>::from_ptr,
                                      raster.get_mask().data());

    // Classification array:
    //   0 = valid pixel or edge-connected nodata (untouched)
    //   1 = interior nodata (to be filled)
    // We initialise all nodata pixels to 1 (interior) and then flood-fill
    // from the edges, setting edge-connected nodata back to 0.
    lsdtt_array2d<uint8_t> interior(nrows, ncols, uint8_t{0});

    // Mark all nodata pixels as potentially interior
    for (int r = 0; r < nrows; ++r)
    {
        for (int c = 0; c < ncols; ++c)
        {
            if (mask(r, c) == 0)
            {
                interior(r, c) = uint8_t{1};
            }
        }
    }

    // Flood fill from edge nodata pixels using a FIFO queue.
    // Any nodata pixel reachable from the edge is edge-connected -> set to 0.
    std::queue<std::pair<int,int>> q;

    auto enqueue_if_edge_nodata = [&](int r, int c)
    {
        if (interior(r, c) == 1)
        {
            interior(r, c) = uint8_t{0};  // mark as edge-connected
            q.push({r, c});
        }
    };

    // Seed from all four raster edges
    for (int c = 0; c < ncols; ++c)
    {
        enqueue_if_edge_nodata(0,        c);
        enqueue_if_edge_nodata(nrows - 1, c);
    }
    for (int r = 1; r < nrows - 1; ++r)
    {
        enqueue_if_edge_nodata(r, 0);
        enqueue_if_edge_nodata(r, ncols - 1);
    }

    // Flood fill — 8-connected through nodata pixels
    while (!q.empty())
    {
        auto [r, c] = q.front();
        q.pop();

        for (int n = 0; n < 8; ++n)
        {
            const int nr = r + k_nd_dr[n];
            const int nc = c + k_nd_dc[n];

            if (nr < 0 || nr >= nrows || nc < 0 || nc >= ncols) continue;
            if (interior(nr, nc) != 1) continue;  // already classified or valid

            interior(nr, nc) = uint8_t{0};  // edge-connected
            q.push({nr, nc});
        }
    }

    return interior;
}

// =============================================================================
// Step 2a — iterative distance-weighted average
// =============================================================================

// One pass over the interior void pixels.
// Each interior nodata pixel that has at least one valid 8-neighbour is
// assigned the distance-weighted mean of its valid neighbours.
// Returns the number of pixels filled in this pass.
static int iterative_average_pass(
    lsdtt_array2d<float>&         dem,
    lsdtt_array2d<uint8_t>&       interior,  // modified in place as pixels are filled
    const lsdtt_raster_meta&      meta,
    int nrows, int ncols,
    int n_threads)
{
    const float nodata = static_cast<float>(meta.no_data_value);

    // Collect pixels to fill this pass — those with at least one valid neighbour
    // We compute updates into a separate array to avoid order-dependence
    lsdtt_array2d<float>   new_vals(nrows, ncols, nodata);
    lsdtt_array2d<uint8_t> filled_this_pass(nrows, ncols, uint8_t{0});

    // Distances: cardinal = 1.0, diagonal = 1/sqrt(2) for inverse weighting
    static const float k_inv_dist[8] = {
        1.0f,
        1.0f / 1.41421356f,
        1.0f,
        1.0f / 1.41421356f,
        1.0f,
        1.0f / 1.41421356f,
        1.0f,
        1.0f / 1.41421356f
    };

#ifdef _OPENMP
    #pragma omp parallel for collapse(2) schedule(dynamic) num_threads(n_threads)
#endif
    for (int r = 0; r < nrows; ++r)
    {
        for (int c = 0; c < ncols; ++c)
        {
            if (interior(r, c) == 0) continue;  // not an interior void pixel

            float weight_sum = 0.0f;
            float value_sum  = 0.0f;

            for (int n = 0; n < 8; ++n)
            {
                const int nr = r + k_nd_dr[n];
                const int nc = c + k_nd_dc[n];

                if (nr < 0 || nr >= nrows || nc < 0 || nc >= ncols) continue;
                if (interior(nr, nc) != 0) continue;  // neighbour is also void
                if (meta.is_nodata(dem(nr, nc))) continue;  // neighbour is nodata

                weight_sum += k_inv_dist[n];
                value_sum  += k_inv_dist[n] * dem(nr, nc);
            }

            if (weight_sum > 0.0f)
            {
                new_vals(r, c)          = value_sum / weight_sum;
                filled_this_pass(r, c) = uint8_t{1};
            }
        }
    }

    // Apply updates and count filled pixels
    int n_filled = 0;
    for (int r = 0; r < nrows; ++r)
    {
        for (int c = 0; c < ncols; ++c)
        {
            if (filled_this_pass(r, c) != 0)
            {
                dem(r, c)      = new_vals(r, c);
                interior(r, c) = uint8_t{0};  // no longer void
                ++n_filled;
            }
        }
    }

    return n_filled;
}

static void fill_iterative_weighted_average(
    lsdtt_array2d<float>&    dem,
    lsdtt_array2d<uint8_t>&  interior,
    const lsdtt_raster_meta& meta,
    int nrows, int ncols,
    int max_iter, int n_threads)
{
    int total_filled = 0;
    int iter         = 0;

    while (iter < max_iter)
    {
        int filled = iterative_average_pass(dem, interior, meta,
                                            nrows, ncols, n_threads);
        total_filled += filled;
        ++iter;

        if (filled == 0) break;  // no more pixels can be filled
    }

    // Count any remaining unfilled pixels
    int remaining = 0;
    for (int r = 0; r < nrows; ++r)
    {
        for (int c = 0; c < ncols; ++c)
        {
            if (interior(r, c) != 0) ++remaining;
        }
    }

    std::cout << "lsdtt::fill_interior_nodata (iterative): filled "
              << total_filled << " pixels in " << iter << " passes.";
    if (remaining > 0)
    {
        std::cout << " WARNING: " << remaining
                  << " interior void pixels could not be filled "
                  << "(isolated from all valid neighbours).";
    }
    std::cout << "\n";
}

// =============================================================================
// Step 2b — Laplacian inpainting (conjugate gradient solver)
// =============================================================================
//
// Formulation:
//   For each interior void pixel p, the discrete Laplacian equation is:
//     z(p) = (1/k) * sum of z(neighbours of p)
//   where k is the number of valid or void neighbours (all 8 in the interior).
//
//   This is equivalent to solving A*x = b where:
//     - x is the vector of unknown elevations at void pixels
//     - A is the graph Laplacian matrix of the void region
//     - b encodes the known boundary values (valid pixels adjacent to the void)
//
//   We solve this using a matrix-free conjugate gradient iteration — no
//   external sparse matrix library is required.
//
// The system is symmetric positive definite so CG is guaranteed to converge.

// Index void pixels — assign each interior void pixel a sequential integer ID.
// Returns the mapping from (row, col) to void ID, and the reverse mapping.
static int index_void_pixels(
    const lsdtt_array2d<uint8_t>& interior,
    lsdtt_array2d<int>&           void_id,   // output: -1 = not void, else ID
    std::vector<int>&             void_rows, // output: row of void pixel i
    std::vector<int>&             void_cols, // output: col of void pixel i
    int nrows, int ncols)
{
    int n_void = 0;
    for (int r = 0; r < nrows; ++r)
    {
        for (int c = 0; c < ncols; ++c)
        {
            if (interior(r, c) != 0)
            {
                void_id(r, c) = n_void;
                void_rows.push_back(r);
                void_cols.push_back(c);
                ++n_void;
            }
        }
    }
    return n_void;
}

// Matrix-vector product A*x for the Laplacian system.
// A*x[i] = n_nbr[i] * x[i] - sum of x[j] for void neighbours j
//         - sum of dem[j] for known-value neighbours j (contributes to b)
// Here we compute A*x — the b contribution is handled separately.
static void laplacian_matvec(
    const std::vector<float>&  x,
    std::vector<float>&        result,
    const lsdtt_array2d<uint8_t>& interior,
    const lsdtt_array2d<int>&  void_id,
    const std::vector<int>&    void_rows,
    const std::vector<int>&    void_cols,
    int nrows, int ncols, int n_void)
{
    result.assign(n_void, 0.0f);

    for (int i = 0; i < n_void; ++i)
    {
        const int r = void_rows[i];
        const int c = void_cols[i];
        int n_nbr = 0;

        for (int n = 0; n < 8; ++n)
        {
            const int nr = r + k_nd_dr[n];
            const int nc = c + k_nd_dc[n];

            if (nr < 0 || nr >= nrows || nc < 0 || nc >= ncols) continue;
            ++n_nbr;

            if (interior(nr, nc) != 0)
            {
                // Void neighbour — contributes to A*x
                const int j = void_id(nr, nc);
                result[i] -= x[j];
            }
            // Known-value neighbours contribute to b, handled separately
        }

        result[i] += static_cast<float>(n_nbr) * x[i];
    }
}

// Build the right-hand side vector b.
// b[i] = sum of known dem values for all valid neighbours of void pixel i
static std::vector<float> build_rhs(
    const lsdtt_array2d<float>&   dem,
    const lsdtt_array2d<uint8_t>& interior,
    const std::vector<int>&       void_rows,
    const std::vector<int>&       void_cols,
    int nrows, int ncols, int n_void)
{
    std::vector<float> b(n_void, 0.0f);

    for (int i = 0; i < n_void; ++i)
    {
        const int r = void_rows[i];
        const int c = void_cols[i];

        for (int n = 0; n < 8; ++n)
        {
            const int nr = r + k_nd_dr[n];
            const int nc = c + k_nd_dc[n];

            if (nr < 0 || nr >= nrows || nc < 0 || nc >= ncols) continue;
            if (interior(nr, nc) != 0) continue;  // void neighbour — in A, not b

            // Known valid neighbour — contributes to b
            b[i] += dem(nr, nc);
        }
    }

    return b;
}

// Conjugate gradient solver for A*x = b.
// The system matrix A is never formed explicitly — it is applied via matvec.
//
// Numerical safety measures:
//   - Dot products accumulated in double to avoid float32 rounding errors
//     that cause stagnation or divergence for elevation-scale values.
//   - Warm start from b[i] / max(n_boundary_neighbours, 1) (local mean of
//     known boundary values) rather than raw b[i] sum.
//   - p_dot_Ap <= 0 triggers an early-exit restart rather than silent divergence.
//   - Residual r is recomputed from scratch every k_cg_restart_interval
//     iterations to prevent drift from accumulated rounding errors.
//   - Voids with no boundary conditions (b entirely zero) are detected and
//     returned as zero rather than attempting to solve an underdetermined system.
static constexpr int k_cg_restart_interval = 50;

static std::vector<float> conjugate_gradient(
    const std::vector<float>&     b,
    const lsdtt_array2d<uint8_t>& interior,
    const lsdtt_array2d<int>&     void_id,
    const std::vector<int>&       void_rows,
    const std::vector<int>&       void_cols,
    int nrows, int ncols,
    int n_void,
    int max_iter, float tolerance)
{
    // Check for a singular system — void with no boundary conditions at all.
    // This happens when a void pixel is entirely surrounded by other void pixels.
    // In this case b is all zero and the system is underdetermined.
    // Return zeros rather than attempting to solve.
    double b_norm_sq_d = 0.0;
    for (float v : b) b_norm_sq_d += static_cast<double>(v) * v;

    if (b_norm_sq_d < 1e-30)
    {
        std::cout << "lsdtt::fill_interior_nodata (laplacian): void region has "
                  << "no boundary conditions — returning zero fill.\n";
        return std::vector<float>(n_void, 0.0f);
    }

    // Warm start: local mean of known boundary values at each void pixel.
    // Count how many known (non-void) neighbours each pixel has.
    std::vector<float> x(n_void, 0.0f);
    for (int i = 0; i < n_void; ++i)
    {
        const int r = void_rows[i];
        const int c = void_cols[i];
        int n_known = 0;
        for (int n = 0; n < 8; ++n)
        {
            const int nr = r + k_nd_dr[n];
            const int nc = c + k_nd_dc[n];
            if (nr < 0 || nr >= nrows || nc < 0 || nc >= ncols) continue;
            if (interior(nr, nc) == 0 && b[i] != 0.0f) ++n_known;
        }
        x[i] = (n_known > 0) ? (b[i] / static_cast<float>(n_known)) : 0.0f;
    }

    // r = b - A*x  (computed in float, converted to double for accumulation)
    std::vector<float> Ax(n_void);
    laplacian_matvec(x, Ax, interior, void_id, void_rows, void_cols,
                     nrows, ncols, n_void);

    std::vector<float> r(n_void), p(n_void);
    for (int i = 0; i < n_void; ++i)
    {
        r[i] = b[i] - Ax[i];
        p[i] = r[i];
    }

    // All dot products accumulated in double to avoid float32 rounding issues
    // with elevation-scale values (typical DEMs: 0-8000m, dot products can
    // reach 10^9 for even modest void sizes, exceeding float32 precision)
    double r_dot_r = 0.0;
    for (float v : r) r_dot_r += static_cast<double>(v) * v;

    const double tol_sq = static_cast<double>(tolerance) * tolerance
                          * std::max(b_norm_sq_d, 1.0);

    std::vector<float> Ap(n_void);
    int iter = 0;

    while (r_dot_r > tol_sq && iter < max_iter)
    {
        // Periodic restart: recompute r from scratch to prevent drift from
        // accumulated floating point errors in the residual update
        if (iter > 0 && iter % k_cg_restart_interval == 0)
        {
            laplacian_matvec(x, Ax, interior, void_id, void_rows, void_cols,
                             nrows, ncols, n_void);
            r_dot_r = 0.0;
            for (int i = 0; i < n_void; ++i)
            {
                r[i]     = b[i] - Ax[i];
                p[i]     = r[i];    // restart search direction too
                r_dot_r += static_cast<double>(r[i]) * r[i];
            }
            if (r_dot_r <= tol_sq) break;
        }

        laplacian_matvec(p, Ap, interior, void_id, void_rows, void_cols,
                         nrows, ncols, n_void);

        double p_dot_Ap = 0.0;
        for (int i = 0; i < n_void; ++i)
        {
            p_dot_Ap += static_cast<double>(p[i]) * Ap[i];
        }

        // p_dot_Ap must be strictly positive for SPD systems.
        // If it is not, the solver has broken down due to numerical error.
        // Restart from the current x rather than diverging.
        if (p_dot_Ap <= 0.0)
        {
            std::cout << "lsdtt::fill_interior_nodata (laplacian): "
                      << "CG breakdown at iteration " << iter
                      << " — restarting.\n";
            laplacian_matvec(x, Ax, interior, void_id, void_rows, void_cols,
                             nrows, ncols, n_void);
            r_dot_r = 0.0;
            for (int i = 0; i < n_void; ++i)
            {
                r[i]     = b[i] - Ax[i];
                p[i]     = r[i];
                r_dot_r += static_cast<double>(r[i]) * r[i];
            }
            if (r_dot_r <= tol_sq) break;
            continue;
        }

        const double alpha = r_dot_r / p_dot_Ap;

        for (int i = 0; i < n_void; ++i)
        {
            x[i] += static_cast<float>(alpha) * p[i];
            r[i] -= static_cast<float>(alpha) * Ap[i];
        }

        const double r_dot_r_new = [&]()
        {
            double s = 0.0;
            for (float v : r) s += static_cast<double>(v) * v;
            return s;
        }();

        const double beta = r_dot_r_new / r_dot_r;
        for (int i = 0; i < n_void; ++i)
        {
            p[i] = r[i] + static_cast<float>(beta) * p[i];
        }

        r_dot_r = r_dot_r_new;
        ++iter;
    }

    std::cout << "lsdtt::fill_interior_nodata (laplacian): CG finished after "
              << iter << " iterations, relative residual = "
              << std::sqrt(r_dot_r / std::max(b_norm_sq_d, 1.0)) << "\n";

    return x;
}

static void fill_laplacian(
    lsdtt_array2d<float>&    dem,
    lsdtt_array2d<uint8_t>&  interior,
    int nrows, int ncols,
    int max_iter, float tolerance)
{
    // Index void pixels
    lsdtt_array2d<int> void_id(nrows, ncols, -1);
    std::vector<int>   void_rows, void_cols;

    int n_void = index_void_pixels(interior, void_id,
                                   void_rows, void_cols, nrows, ncols);

    if (n_void == 0) return;

    std::cout << "lsdtt::fill_interior_nodata (laplacian): solving for "
              << n_void << " void pixels.\n";

    // Build right-hand side
    std::vector<float> b = build_rhs(dem, interior, void_rows, void_cols,
                                     nrows, ncols, n_void);

    // Solve A*x = b
    std::vector<float> x = conjugate_gradient(
        b, interior, void_id, void_rows, void_cols,
        nrows, ncols, n_void, max_iter, tolerance
    );

    // Write solution back into dem
    for (int i = 0; i < n_void; ++i)
    {
        dem(void_rows[i], void_cols[i]) = x[i];
        interior(void_rows[i], void_cols[i]) = uint8_t{0};
    }
}

// =============================================================================
// fill_interior_nodata — public function
// =============================================================================

lsdtt_raster_cpu fill_interior_nodata(
    const lsdtt_raster_cpu&raster,
    nodata_fill_method method,
    int                max_iter,
    float              tolerance,
    int                n_threads)
{
    const lsdtt_raster_meta& meta = raster.get_meta();
    if (!meta.is_valid())
    {
        throw std::runtime_error(
            "lsdtt::fill_interior_nodata: invalid raster metadata"
        );
    }

    const int   nrows  = meta.n_rows;
    const int   ncols  = meta.n_cols;

    // Working copy of elevations
    lsdtt_array2d<float> dem(nrows, ncols,
                             lsdtt_array2d<float>::from_ptr,
                             raster.get_data().data());

    // Resolve thread count
    int threads = n_threads;
#ifdef _OPENMP
    if (threads <= 0) threads = omp_get_max_threads();
#else
    threads = 1;
#endif

    // Step 1 — classify nodata pixels
    lsdtt_array2d<uint8_t> interior = identify_interior_nodata(raster);

    // Count interior void pixels before filling
    int n_interior = 0;
    for (int r = 0; r < nrows; ++r)
    {
        for (int c = 0; c < ncols; ++c)
        {
            if (interior(r, c) != 0) ++n_interior;
        }
    }

    if (n_interior == 0)
    {
        std::cout << "lsdtt::fill_interior_nodata: no interior void pixels found.\n";
        // Return a copy of the input unchanged
        std::vector<float>   data_copy = raster.get_data();
        std::vector<uint8_t> mask_copy = raster.get_mask();
        return lsdtt_raster_cpu(
            std::move(data_copy),
            std::move(mask_copy),
            raster.get_meta()
        );
    }

    std::cout << "lsdtt::fill_interior_nodata: found "
              << n_interior << " interior void pixels to fill.\n";

    // Step 2 — fill using the selected method
    if (method == nodata_fill_method::iterative_weighted_average)
    {
        fill_iterative_weighted_average(dem, interior, meta,
                                        nrows, ncols, max_iter, threads);
    }
    else  // laplacian
    {
        fill_laplacian(dem, interior, nrows, ncols, max_iter, tolerance);
    }

    // Build output mask — valid wherever dem is not nodata
    lsdtt_array2d<uint8_t> out_mask(nrows, ncols, uint8_t{0});
    for (int r = 0; r < nrows; ++r)
    {
        for (int c = 0; c < ncols; ++c)
        {
            if (!meta.is_nodata(dem(r, c))) out_mask(r, c) = uint8_t{1};
        }
    }

    return lsdtt_raster_cpu(
        std::move(dem.get_vector()),
        std::move(out_mask.get_vector()),
        meta
    );
}

// =============================================================================
// clip_to_maximum_rectangle
// =============================================================================

// Largest rectangle in a histogram — classic O(n) stack algorithm.
//
// Given a histogram of bar heights, finds the largest rectangle that fits
// entirely within the bars. Returns the rectangle as {left_col, width, height}.
//
// We also need the row at which this rectangle ends (provided by the caller
// as current_row) to recover the top-left pixel in the raster.
struct hist_rect
{
    int left;    // leftmost column of the rectangle
    int width;   // width in columns
    int height;  // height in rows
    int area()   const { return width * height; }
};

static hist_rect largest_rect_in_histogram(
    const std::vector<int>& heights,
    int ncols)
{
    // Stack stores column indices of bars whose rectangles are still open
    std::vector<int> stack;
    stack.reserve(ncols + 1);

    hist_rect best{ 0, 0, 0 };

    auto update = [&](int right_boundary)
    {
        while (!stack.empty())
        {
            const int h = heights[stack.back()];
            stack.pop_back();
            const int left  = stack.empty() ? 0 : stack.back() + 1;
            const int width = right_boundary - left;
            if (h * width > best.area())
            {
                best = { left, width, h };
            }
        }
    };

    for (int c = 0; c < ncols; ++c)
    {
        // Pop bars taller than the current bar — their rectangles cannot
        // extend past this column
        while (!stack.empty() && heights[stack.back()] > heights[c])
        {
            const int h     = heights[stack.back()];
            stack.pop_back();
            const int left  = stack.empty() ? 0 : stack.back() + 1;
            const int width = c - left;
            if (h * width > best.area())
            {
                best = { left, width, h };
            }
        }
        stack.push_back(c);
    }

    // Flush remaining bars
    update(ncols);

    return best;
}

lsdtt_raster_cpu clip_to_maximum_rectangle(
    const lsdtt_raster_cpu&raster)
{
    const lsdtt_raster_meta& meta = raster.get_meta();
    if (!meta.is_valid())
    {
        throw std::runtime_error(
            "lsdtt::clip_to_maximum_rectangle: invalid raster metadata"
        );
    }

    const int   nrows  = meta.n_rows;
    const int   ncols  = meta.n_cols;

    const lsdtt_array2d<uint8_t> mask(nrows, ncols,
                                      lsdtt_array2d<uint8_t>::from_ptr,
                                      raster.get_mask().data());

    const lsdtt_array2d<float> dem(nrows, ncols,
                                   lsdtt_array2d<float>::from_ptr,
                                   raster.get_data().data());

    // Height array: height[c] = number of consecutive valid pixels from
    // the current row upward (inclusive) in column c.
    // A nodata pixel resets the height to 0.
    std::vector<int> heights(ncols, 0);

    // Track the best rectangle found so far
    struct best_rect
    {
        int row_end  = -1;  // row index where the rectangle's bottom edge is
        int col_left =  0;  // leftmost column
        int width    =  0;
        int height   =  0;
        int area()   const { return width * height; }
    } best;

    for (int r = 0; r < nrows; ++r)
    {
        // Update height array for this row
        for (int c = 0; c < ncols; ++c)
        {
            if (mask(r, c) != 0 && !meta.is_nodata(dem(r, c)))
            {
                heights[c] += 1;
            }
            else
            {
                heights[c] = 0;  // nodata resets the run
            }
        }

        // Find the largest rectangle in this row's histogram
        hist_rect candidate = largest_rect_in_histogram(heights, ncols);

        if (candidate.area() > best.area())
        {
            best.row_end  = r;
            best.col_left = candidate.left;
            best.width    = candidate.width;
            best.height   = candidate.height;
        }
    }

    if (best.area() == 0)
    {
        throw std::runtime_error(
            "lsdtt::clip_to_maximum_rectangle: no valid rectangle found — "
            "the raster may be entirely nodata"
        );
    }

    // Recover the pixel coordinates of the rectangle
    // row_end is the bottom row of the rectangle (0-indexed, row 0 = north)
    // row_start is the top row
    const int row_start = best.row_end - best.height + 1;
    const int row_end   = best.row_end;               // inclusive
    const int col_start = best.col_left;
    const int col_end   = best.col_left + best.width - 1; // inclusive

    const int out_nrows = best.height;
    const int out_ncols = best.width;

    std::cout << "lsdtt::clip_to_maximum_rectangle:\n"
              << "  Input  : " << nrows << " x " << ncols << " pixels\n"
              << "  Output : " << out_nrows << " x " << out_ncols << " pixels\n"
              << "  Rows   : " << row_start << " to " << row_end   << "\n"
              << "  Cols   : " << col_start << " to " << col_end   << "\n"
              << "  Retained : "
              << (100.0 * best.area()) / static_cast<double>(nrows * ncols)
              << "% of original pixels\n";

    // Extract the clipped data and build output mask
    lsdtt_array2d<float>   out_dem(out_nrows, out_ncols,
                                   static_cast<float>(meta.no_data_value));
    lsdtt_array2d<uint8_t> out_mask(out_nrows, out_ncols, uint8_t{0});

    for (int r = 0; r < out_nrows; ++r)
    {
        for (int c = 0; c < out_ncols; ++c)
        {
            const int src_r = row_start + r;
            const int src_c = col_start + c;
            out_dem(r, c)  = dem(src_r, src_c);
            out_mask(r, c) = mask(src_r, src_c);
        }
    }

    // Build updated metadata
    // x_min shifts right by col_start columns
    // y_min shifts up by (nrows - row_end - 1) rows from the original bottom
    //   because row 0 is the northernmost row and y_min is the southern edge
    lsdtt_raster_meta out_meta = meta;
    out_meta.n_rows   = out_nrows;
    out_meta.n_cols   = out_ncols;
    out_meta.x_min    = meta.x_min + col_start * meta.cell_size;
    out_meta.y_min    = meta.y_min + (nrows - row_end - 1) * meta.cell_size;

    return lsdtt_raster_cpu(
        std::move(out_dem.get_vector()),
        std::move(out_mask.get_vector()),
        out_meta
    );
}

} // namespace lsdtt
