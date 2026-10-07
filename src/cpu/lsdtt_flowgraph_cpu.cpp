// lsdtt_flowgraph_cpu.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// Implementation of lsdtt_flowgraph_cpu.
// See lsdtt_flowgraph_cpu.hpp for design notes.
//
// References
// ----------
//   Tarboton D.G. (1997) A new method for the determination of flow
//     directions and upslope areas in grid digital elevation models.
//     Water Resources Research 33(2): 309-319.
//
//   Quinn P., Beven K., Chevallier P., Planchon O. (1991) The prediction
//     of hillslope flow paths for distributed hydrological modelling using
//     digital terrain models. Hydrological Processes 5(1): 59-79.
//
//   Freeman T.G. (1991) Calculating catchment area with divergent flow
//     based on a regular grid. Computers & Geosciences 17(3): 413-422.
//
//   Anand S.K., Hooshyar M., Porporato A. (2020) Linear layout of multiple
//     flow-direction networks for landscape-evolution simulations.
//     Environmental Modelling & Software 133: 104804.

#include "lsdtt_flowgraph_cpu.hpp"
#include "../utils/lsdtt_grid_neighbours.hpp"

#include "../core/lsdtt_data_loader.hpp"
#include "../core/lsdtt_point_input.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <queue>
#include <stdexcept>
#include <string>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace lsdtt
{

// =============================================================================
// Local constants (file scope)
// =============================================================================
namespace
{

// Tarboton (1997) Table 1 - eight triangular facets.
//   e1 = cardinal neighbour at (row+e1_dr, col+e1_dc)
//   e2 = diagonal neighbour at (row+e2_dr, col+e2_dc)
//   ac, af: sector index and sign used to convert the within-facet angle r
//           to a compass-referenced flow direction:
//             angle = af * r + ac * (pi/2)
constexpr int k_e1_dr[8] = {  0, -1, -1,  0,  0,  1,  1,  0 };
constexpr int k_e1_dc[8] = {  1,  0,  0, -1, -1,  0,  0,  1 };
constexpr int k_e2_dr[8] = { -1, -1, -1, -1,  1,  1,  1,  1 };
constexpr int k_e2_dc[8] = {  1,  1, -1, -1, -1, -1,  1,  1 };
constexpr int k_ac[8]    = {  0,  1,  1,  2,  2,  3,  3,  4 };
constexpr int k_af[8]    = {  1, -1,  1, -1,  1, -1,  1, -1 };

// MFD 8-neighbour offsets.  Indices 0-3 are cardinal (distance = d);
// indices 4-7 are diagonal (distance = dsqrt2).
static constexpr const auto& k_nb_dr = grid::d8_cardinal_first_row;
static constexpr const auto& k_nb_dc = grid::d8_cardinal_first_col;

constexpr float k_pi_over_4  = static_cast<float>(M_PI / 4.0);
constexpr float k_sqrt2      = grid::k_sqrt2_f;
constexpr float k_rad_to_deg = static_cast<float>(180.0 / M_PI);

} // anonymous namespace


// =============================================================================
// Constructor
// =============================================================================

lsdtt_flowgraph_cpu::lsdtt_flowgraph_cpu(
    const lsdtt_array2d_f&         dem,
    const lsdtt_array2d_u8&        mask,
    const lsdtt_raster_meta&       meta,
    const lsdtt_flowgraph_options& options)
{
    meta_    = meta;
    options_ = options;
    n_rows_  = meta.n_rows;
    n_cols_  = meta.n_cols;
    n_cells_ = static_cast<std::int64_t>(n_rows_) * n_cols_;
    bc_      = parse_boundary_conditions(options.boundary_conditions);

    meta_.check_projected_metres();

    if (dem.n_rows() != n_rows_ || dem.n_cols() != n_cols_)
        throw std::invalid_argument(
            "lsdtt_flowgraph_cpu: DEM dimensions do not match meta");
    if (mask.n_rows() != n_rows_ || mask.n_cols() != n_cols_)
        throw std::invalid_argument(
            "lsdtt_flowgraph_cpu: mask dimensions do not match meta");
    if (!options_.forced_outlet_mask.empty() &&
        static_cast<std::int64_t>(options_.forced_outlet_mask.size()) != n_cells_)
    {
        throw std::invalid_argument(
            "lsdtt_flowgraph_cpu: forced_outlet_mask must be empty or have "
            "one entry per raster cell");
    }

    // Build validity vector from both mask and DEM validity
    valid_.assign(n_cells_, 0);
    n_valid_cells_ = 0;

    const float ndv = static_cast<float>(meta_.no_data_value);

    for (int row = 0; row < n_rows_; ++row)
    {
        for (int col = 0; col < n_cols_; ++col)
        {
            const std::int64_t flat =
                static_cast<std::int64_t>(row) * n_cols_ + col;

            const float z = dem(row, col);

            if (mask(row, col) != 0 &&
                z != ndv &&
                std::isfinite(z))
            {
                valid_[flat] = 1;
                ++n_valid_cells_;
            }
        }
    }

    if (!options_.forced_outlet_mask.empty())
    {
        for (std::int64_t flat = 0; flat < n_cells_; ++flat)
        {
            if (options_.forced_outlet_mask[static_cast<std::size_t>(flat)] != 0 &&
                valid_[flat] == 0)
            {
                throw std::invalid_argument(
                    "lsdtt_flowgraph_cpu: forced outlet must be a valid raster cell");
            }
        }
    }


    if (options.method == lsdtt_flowgraph_method::mfd)
    {
        if (!std::isfinite(options.mfd_exponent) ||
            options.mfd_exponent <= 0.0f)
        {
            throw std::invalid_argument(
                "lsdtt_flowgraph_cpu: mfd_exponent must be finite and > 0");
        }
    }

    if (options.method == lsdtt_flowgraph_method::d_infinity)
        build_dinf(dem);
    else
        build_mfd(dem, options.mfd_exponent);



    build_topological_order();
}


// =============================================================================
// clear
// =============================================================================

void lsdtt_flowgraph_cpu::clear() noexcept
{
    valid_.clear();
    receiver_offsets_.clear();
    receiver_indices_.clear();
    receiver_weights_.clear();
    topological_order_.clear();
    donor_count_.clear();
    dinf_slopes_.clear();
    dinf_direction_.clear();
    n_cells_       = 0;
    n_valid_cells_ = 0;
}


// =============================================================================
// options getter
// =============================================================================
const lsdtt_flowgraph_options&
lsdtt_flowgraph_cpu::options() const noexcept
{
    return options_;
}


// =============================================================================
// n_outlets
// =============================================================================

std::int64_t lsdtt_flowgraph_cpu::n_outlets() const noexcept
{
    std::int64_t count = 0;

    for (std::int64_t flat = 0; flat < n_cells_; ++flat)
    {
        if (valid_[flat] == 0) continue;

        if (receiver_offsets_[flat + 1] == receiver_offsets_[flat])
            ++count;
    }

    return count;
}


// =============================================================================
// outlet_indices
// =============================================================================

std::vector<std::int64_t> lsdtt_flowgraph_cpu::outlet_indices() const
{
    std::vector<std::int64_t> outlets;
    outlets.reserve(static_cast<std::size_t>(n_outlets()));

    for (std::int64_t flat = 0; flat < n_cells_; ++flat)
    {
        if (valid_[flat] == 0) continue;

        if (receiver_offsets_[flat + 1] == receiver_offsets_[flat])
            outlets.push_back(flat);
    }

    return outlets;
}



// =============================================================================
// parse_boundary_conditions
// =============================================================================

lsdtt_flowgraph_cpu::boundary_flags
lsdtt_flowgraph_cpu::parse_boundary_conditions(const std::string& bc) const 
{
    if (bc.size() != 4)
    {
        throw std::invalid_argument(
            "lsdtt_flowgraph_cpu: boundary_conditions must have exactly 4 characters "
            "(north, east, south, west)");
    }

    auto is_allowed = [](char c) -> bool
    {
        return c == 'o' || c == 'O' ||
               c == 'n' || c == 'N' ||
               c == 'p' || c == 'P' ||
               c == 'b' || c == 'B';
    };

    for (char c : bc)
    {
        if (!is_allowed(c))
        {
            throw std::invalid_argument(
                "lsdtt_flowgraph_cpu: boundary_conditions contains an invalid character. "
                "Allowed values are o, n, p, b.");
        }
    }

    boundary_flags f;
    auto is_periodic = [](char c) { return c == 'p' || c == 'P'; };
    auto is_forced   = [](char c) { return c == 'b' || c == 'B'; };

    f.north_periodic = is_periodic(bc[0]);
    f.north_forced   = is_forced(bc[0]);

    f.east_periodic  = is_periodic(bc[1]);
    f.east_forced    = is_forced(bc[1]);

    f.south_periodic = is_periodic(bc[2]);
    f.south_forced   = is_forced(bc[2]);

    f.west_periodic  = is_periodic(bc[3]);
    f.west_forced    = is_forced(bc[3]);

    // Periodic boundaries must be symmetric pairs:
    // north <-> south and east <-> west.
    if (f.north_periodic != f.south_periodic)
    {
        throw std::invalid_argument(
            "lsdtt_flowgraph_cpu: north and south periodic boundary flags "
            "must match");
    }

    if (f.east_periodic != f.west_periodic)
    {
        throw std::invalid_argument(
            "lsdtt_flowgraph_cpu: east and west periodic boundary flags "
            "must match");
    }

    return f;
}




// =============================================================================
// is_forced_outlet
// =============================================================================

bool lsdtt_flowgraph_cpu::is_forced_outlet(int row, int col) const noexcept
{
    if (bc_.north_forced && row == 0)             return true;
    if (bc_.south_forced && row == n_rows_ - 1)   return true;
    if (bc_.east_forced  && col == n_cols_ - 1)   return true;
    if (bc_.west_forced  && col == 0)             return true;
    return false;
}

bool lsdtt_flowgraph_cpu::is_explicit_outlet(std::int64_t flat) const noexcept
{
    return !options_.forced_outlet_mask.empty() &&
           flat >= 0 && flat < n_cells_ &&
           options_.forced_outlet_mask[static_cast<std::size_t>(flat)] != 0;
}

// =============================================================================
// touches_interior_nodata
// =============================================================================

bool lsdtt_flowgraph_cpu::touches_interior_nodata(int row, int col) const noexcept
{
    const std::int64_t flat =
        static_cast<std::int64_t>(row) * n_cols_ + col;

    if (flat < 0 || flat >= n_cells_) return false;
    if (valid_[flat] == 0) return false;

    for (int k = 0; k < 8; ++k)
    {
        const int nr = row + k_nb_dr[k];
        const int nc = col + k_nb_dc[k];

        // If the raw neighbour lies outside the domain, this is an outer edge,
        // not an interior nodata contact, so skip it.
        if (nr < 0 || nr >= n_rows_ || nc < 0 || nc >= n_cols_)
            continue;

        const std::int64_t nflat =
            static_cast<std::int64_t>(nr) * n_cols_ + nc;

        if (valid_[nflat] == 0)
            return true;
    }

    return false;
}




// =============================================================================
// resolve_neighbor
// =============================================================================

lsdtt_pixel_rc
lsdtt_flowgraph_cpu::resolve_neighbor(
    int row, int col, int drow, int dcol) const noexcept
{
    int nr = row + drow;
    int nc = col + dcol;

    // North boundary
    if (nr < 0)
    {
        if (bc_.north_periodic) nr = n_rows_ - 1;
        else return {-1, -1};
    }
    // South boundary
    if (nr >= n_rows_)
    {
        if (bc_.south_periodic) nr = 0;
        else return {-1, -1};
    }
    // West boundary
    if (nc < 0)
    {
        if (bc_.west_periodic) nc = n_cols_ - 1;
        else return {-1, -1};
    }
    // East boundary
    if (nc >= n_cols_)
    {
        if (bc_.east_periodic) nc = 0;
        else return {-1, -1};
    }

    return {nr, nc};
}


// =============================================================================
// build_dinf
//
// Implements Tarboton (1997) D-infinity routing.
//
// For each valid cell (row, col), 8 triangular facets are evaluated.
// Each facet is defined by the centre cell (e0), a cardinal neighbour (e1),
// and a diagonal neighbour (e2).
//
//   s1 = (e0 - e1) / d          - slope along cardinal edge
//   s2 = (e1 - e2) / d          - cross-slope within facet (Tarboton eq. 2)
//   r  = atan2(s2, s1)          - within-facet steepest angle [0, pi/4]
//
// Stability improvements vs. the original LSDDinfFlowInfo:
//   - atan2(s2, s1) replaces atan(s2/s1) -> no division-by-zero when s1 = 0.
//   - r is clamped to [0, pi/4] cleanly; s is then recomputed from the
//     clamped geometry rather than carried from unclamped arithmetic.
//   - No artificial s1 = 0.00001 perturbation.
//   - Proportions are derived directly from r (not from a compass-bearing
//     decomposition of the final rounded angle).
// =============================================================================

void lsdtt_flowgraph_cpu::build_dinf(
    const lsdtt_array2d_f& dem)
{
    const float cell_size = static_cast<float>(meta_.cell_size);
    const float ndv       = static_cast<float>(meta_.no_data_value);

    // Per-cell receiver data (max 2 receivers for D-infinity)
    struct cell_rec
    {
        std::int64_t rec[2] = {-1, -1};
        float        w[2]   = {0.0f, 0.0f};
        int          n      = 0;
    };

    std::vector<cell_rec> cell_data(n_cells_);

    dinf_slopes_.assign(n_cells_, ndv);
    dinf_direction_.assign(n_cells_, ndv);

    // Parallelise the expensive per-cell facet loop.
    // Each iteration writes only to its own flat index - no race conditions.
    int nthreads = options_.n_threads;
#ifdef _OPENMP
    if (nthreads == 0) nthreads = omp_get_max_threads();
    #pragma omp parallel for schedule(static) num_threads(nthreads)
#endif
    for (int row = 0; row < n_rows_; ++row)
    {
        for (int col = 0; col < n_cols_; ++col)
        {
            const std::int64_t flat = static_cast<std::int64_t>(row) * n_cols_ + col;

            if (valid_[flat] == 0) continue;
            if (is_forced_outlet(row, col) || is_explicit_outlet(flat))
                continue;  // forced outlet: no receivers

            if (options_.interior_nodata_is_outlet &&
                touches_interior_nodata(row, col))
            {
                continue;  // treat as outlet/sink
            }

            const float e0 = dem(row, col);
            if (!std::isfinite(e0) || e0 == ndv) continue;



            float        best_s      = 0.0f;  // only accept facets with s > 0
            float        best_r      = 0.0f;
            int          best_facet  = -1;
            std::int64_t best_e1_flat = -1;
            std::int64_t best_e2_flat = -1;

            for (int a = 0; a < 8; ++a)
            {
                lsdtt_pixel_rc nb1 = resolve_neighbor(row, col, k_e1_dr[a], k_e1_dc[a]);
                lsdtt_pixel_rc nb2 = resolve_neighbor(row, col, k_e2_dr[a], k_e2_dc[a]);

                if (!nb1.is_valid() || !nb2.is_valid()) continue;

                const std::int64_t nb1_flat =
                    static_cast<std::int64_t>(nb1.row) * n_cols_ + nb1.col;
                const std::int64_t nb2_flat =
                    static_cast<std::int64_t>(nb2.row) * n_cols_ + nb2.col;

                if (nb1_flat == flat) continue;   // skip self-receivers on tiny periodic grids
                if (nb2_flat == flat) continue;   // skip self-receivers on tiny periodic grids

                if (valid_[nb1_flat] == 0) continue;
                if (valid_[nb2_flat] == 0) continue;


                const float e1 = dem(nb1.row, nb1.col);
                const float e2 = dem(nb2.row, nb2.col);
                if (!std::isfinite(e1) || !std::isfinite(e2) ||
                    e1 == ndv || e2 == ndv) continue;


                // Only enter if at least one neighbour is lower than e0.
                if (e0 <= e1 && e0 <= e2) continue;

                float r, s;
                const float s1 = (e0 - e1) / cell_size;  // always > 0 when e0 > e1
                const float s2 = (e1 - e2) / cell_size;

                if (e0 > e1 && e0 > e2)
                {
                    // Both lower: facet slope as gradient magnitude.
                    // atan2 handles s1 == 0 gracefully (returns pi/2, clamped below).
                    r = std::atan2(s2, s1);

                    if (r < 0.0f)
                    {
                        // r < 0: steepest descent is outside the facet toward e1
                        r = 0.0f;
                        s = s1;  // s1 > 0 since e0 > e1
                    }
                    else if (r > k_pi_over_4)
                    {
                        // r > pi/4: steepest descent is outside the facet toward e2
                        r = k_pi_over_4;
                        s = (e0 - e2) / (cell_size * k_sqrt2);
                    }
                    else
                    {
                        s = std::sqrt(s1 * s1 + s2 * s2);
                    }
                }
                else if (e0 > e1)
                {
                    // Only e1 is lower: flow along cardinal edge
                    r = 0.0f;
                    s = s1;
                }
                else
                {
                    // Only e2 is lower: flow along diagonal edge
                    r = k_pi_over_4;
                    s = (e0 - e2) / (cell_size * k_sqrt2);
                }

                if (s > best_s)
                {
                    best_s      = s;
                    best_r      = r;
                    best_facet  = a;
                    best_e1_flat = static_cast<std::int64_t>(nb1.row) * n_cols_ + nb1.col;
                    best_e2_flat = static_cast<std::int64_t>(nb2.row) * n_cols_ + nb2.col;
                }
            }

            if (best_facet < 0)
            {
                // No downslope facet found: this cell is a sink.
                // cell_data[flat].n remains 0.
                continue;
            }

            // Flow direction angle (compass bearing in degrees).
            //   Internal angle: af * r + ac * pi/2  (result in radians, math convention)
            //   Compass bearing = 360 - (angle_deg) + 90, wrapped to [0, 360)
            const float angle_rad =
                static_cast<float>(k_af[best_facet]) * best_r
                + static_cast<float>(k_ac[best_facet]) * (static_cast<float>(M_PI) * 0.5f);
            float bearing = 360.0f - (angle_rad * k_rad_to_deg) + 90.0f;
            while (bearing >= 360.0f) bearing -= 360.0f;
            while (bearing <    0.0f) bearing += 360.0f;

            dinf_slopes_[flat]    = best_s;
            dinf_direction_[flat] = bearing;

            // Partition flow between e1 (cardinal) and e2 (diagonal).
            //   p_e1 = (pi/4 - r) / (pi/4)  =  1 - r * 4/pi
            //   p_e2 = r / (pi/4)           =  r * 4/pi
            cell_rec& cr = cell_data[flat];

            if (best_r == 0.0f)
            {
                cr.n      = 1;
                cr.rec[0] = best_e1_flat;
                cr.w[0]   = 1.0f;
            }
            else if (best_r >= k_pi_over_4)
            {
                cr.n      = 1;
                cr.rec[0] = best_e2_flat;
                cr.w[0]   = 1.0f;
            }
            else
            {
                const float p_e2 = best_r / k_pi_over_4;
                const float p_e1 = 1.0f - p_e2;

                // On tiny periodic grids the cardinal and diagonal receivers can
                // wrap to the same cell. Merge them into a single receiver.
                if (best_e1_flat == best_e2_flat)
                {
                    cr.n      = 1;
                    cr.rec[0] = best_e1_flat;
                    cr.w[0]   = 1.0f;
                }
                else
                {
                    cr.n      = 2;
                    cr.rec[0] = best_e1_flat;  cr.w[0] = p_e1;
                    cr.rec[1] = best_e2_flat;  cr.w[1] = p_e2;
                }
            }

        }
    }

    // Build CSR from per-cell data
    receiver_offsets_.resize(n_cells_ + 1, 0);
    for (std::int64_t flat = 0; flat < n_cells_; ++flat)
        receiver_offsets_[flat + 1] = receiver_offsets_[flat] + cell_data[flat].n;

    const std::int64_t n_edges = receiver_offsets_[n_cells_];
    receiver_indices_.resize(n_edges);
    receiver_weights_.resize(n_edges);

    for (std::int64_t flat = 0; flat < n_cells_; ++flat)
    {
        const auto base = receiver_offsets_[flat];
        const cell_rec& cr = cell_data[flat];
        for (int k = 0; k < cr.n; ++k)
        {
            receiver_indices_[base + k] = cr.rec[k];
            receiver_weights_[base + k] = cr.w[k];
        }
    }
}


// =============================================================================
// build_mfd
//
// Implements Quinn et al. (1991) / Freeman (1991) MFD routing.
//
// For each valid cell, slope to each of 8 neighbours is computed.
// Downslope neighbours (slope > 0) receive weighted flow proportional to
// slope^exponent.  Weights are normalised to sum to 1.
//
// exponent = 1.0 reproduces Quinn et al. (1991).
// exponent = 1.1 is Freeman (1991)'s recommendation for natural terrain.
// =============================================================================

void lsdtt_flowgraph_cpu::build_mfd(
    const lsdtt_array2d_f& dem,
    const float exponent)
{
    const float cell_size  = static_cast<float>(meta_.cell_size);
    const float diag_size  = cell_size * k_sqrt2;
    const float ndv        = static_cast<float>(meta_.no_data_value);
    const bool  unit_exp   = (std::abs(exponent - 1.0f) < 1e-6f);

    struct cell_rec
    {
        std::int64_t rec[8];
        float        w[8];
        int          n = 0;
    };

    std::vector<cell_rec> cell_data(n_cells_);

    int nthreads = options_.n_threads;
#ifdef _OPENMP
    if (nthreads == 0) nthreads = omp_get_max_threads();
    #pragma omp parallel for schedule(static) num_threads(nthreads)
#endif
    for (int row = 0; row < n_rows_; ++row)
    {
        for (int col = 0; col < n_cols_; ++col)
        {
            const std::int64_t flat = static_cast<std::int64_t>(row) * n_cols_ + col;

            if (valid_[flat] == 0) continue;
            if (is_forced_outlet(row, col) || is_explicit_outlet(flat))
                continue;

            if (options_.interior_nodata_is_outlet &&
                touches_interior_nodata(row, col))
            {
                continue;  // treat as outlet/sink
            }

            const float e0 = dem(row, col);
            if (!std::isfinite(e0) || e0 == ndv) continue;



            cell_rec& cr = cell_data[flat];
            float total_w = 0.0f;

            for (int k = 0; k < 8; ++k)
            {
                const lsdtt_pixel_rc nb =
                    resolve_neighbor(row, col, k_nb_dr[k], k_nb_dc[k]);
                if (!nb.is_valid()) continue;

                const std::int64_t nb_flat =
                    static_cast<std::int64_t>(nb.row) * n_cols_ + nb.col;

                if (nb_flat == flat) continue;   // skip self-receivers on tiny periodic grids
                if (valid_[nb_flat] == 0) continue;

                const float e_nb = dem(nb.row, nb.col);
                if (!std::isfinite(e_nb) || e_nb == ndv) continue;



                const float dist  = (k < 4) ? cell_size : diag_size;
                const float slope = (e0 - e_nb) / dist;

                if (slope > 0.0f)
                {
                    const float w = unit_exp ? slope : std::pow(slope, exponent);
                    cr.rec[cr.n] = static_cast<std::int64_t>(nb.row) * n_cols_ + nb.col;
                    cr.w[cr.n]   = w;
                    total_w += w;
                    ++cr.n;
                }
            }

            // Merge duplicate receivers that can arise on tiny periodic grids.
            // Sum their unnormalised weights before the final normalisation.
            if (cr.n > 1)
            {
                int write_n = 0;

                for (int k = 0; k < cr.n; ++k)
                {
                    const std::int64_t rec_k = cr.rec[k];
                    const float        w_k   = cr.w[k];

                    bool merged = false;
                    for (int kk = 0; kk < write_n; ++kk)
                    {
                        if (cr.rec[kk] == rec_k)
                        {
                            cr.w[kk] += w_k;
                            merged = true;
                            break;
                        }
                    }

                    if (!merged)
                    {
                        cr.rec[write_n] = rec_k;
                        cr.w[write_n]   = w_k;
                        ++write_n;
                    }
                }

                cr.n = write_n;
            }

            // Recompute total weight after any duplicate merge.
            total_w = 0.0f;
            for (int k = 0; k < cr.n; ++k)
                total_w += cr.w[k];

            // Normalise so weights sum to 1.0
            if (total_w > 0.0f)
            {
                const float inv = 1.0f / total_w;
                for (int k = 0; k < cr.n; ++k) cr.w[k] *= inv;
            }
            else
            {
                cr.n = 0;  // all neighbours higher or nodata - cell is a sink
            }

        }
    }

    // Build CSR
    receiver_offsets_.resize(n_cells_ + 1, 0);
    for (std::int64_t flat = 0; flat < n_cells_; ++flat)
        receiver_offsets_[flat + 1] = receiver_offsets_[flat] + cell_data[flat].n;

    const std::int64_t n_edges = receiver_offsets_[n_cells_];
    receiver_indices_.resize(n_edges);
    receiver_weights_.resize(n_edges);

    for (std::int64_t flat = 0; flat < n_cells_; ++flat)
    {
        const auto base = receiver_offsets_[flat];
        const cell_rec& cr = cell_data[flat];
        for (int k = 0; k < cr.n; ++k)
        {
            receiver_indices_[base + k] = cr.rec[k];
            receiver_weights_[base + k] = cr.w[k];
        }
    }
}


// =============================================================================
// build_topological_order
//
// Kahn's BFS over the receiver DAG, following the multi-receiver extension
// of Anand et al. (2020): a cell enters the queue only once all its donors
// have been processed.
//
// Result: topological_order_[0] is a source (ridge) cell.
//         topological_order_.back() is an outlet (sink) cell.
//
// Cycle detection: if the topological order is shorter than n_valid_cells_
// after BFS completes, a cycle exists in the graph (should not occur for a
// properly pit-filled DEM).
// =============================================================================

void lsdtt_flowgraph_cpu::build_topological_order()
{
    // donor_count_[flat] = number of valid donor cells that send any fraction
    // of flow to flat. Count unique donors, not raw incoming edges, so that
    // duplicate wrapped neighbours on tiny periodic grids do not overcount.
    donor_count_.assign(n_cells_, 0);

    for (std::int64_t flat = 0; flat < n_cells_; ++flat)
    {
        const auto begin = receiver_offsets_[flat];
        const auto end   = receiver_offsets_[flat + 1];

        double weight_sum = 0.0;

        for (auto e = begin; e < end; ++e)
        {
            const std::int64_t rec = receiver_indices_[e];
            const float w          = receiver_weights_[e];

            if (rec < 0 || rec >= n_cells_)
            {
                throw std::runtime_error(
                    "lsdtt_flowgraph_cpu: invalid receiver index in receiver graph");
            }

            if (rec == flat)
            {
                throw std::runtime_error(
                    "lsdtt_flowgraph_cpu: self-receiver detected in receiver graph");
            }

            if (!std::isfinite(w) || w < 0.0f)
            {
                throw std::runtime_error(
                    "lsdtt_flowgraph_cpu: invalid receiver weight in receiver graph");
            }

            weight_sum += static_cast<double>(w);

            bool already_counted = false;
            for (auto ee = begin; ee < e; ++ee)
            {
                if (receiver_indices_[ee] == rec)
                {
                    already_counted = true;
                    break;
                }
            }

            if (!already_counted && donor_count_[rec] < 255)
            {
                ++donor_count_[rec];
            }
        }

        if (end > begin)
        {
            const double tol = 1e-5;
            if (std::abs(weight_sum - 1.0) > tol)
            {
                throw std::runtime_error(
                    "lsdtt_flowgraph_cpu: receiver weights for a donor cell do not sum to 1");
            }
        }
    }




    // Working copy for BFS (Kahn's consumes in-degree counts)
    std::vector<std::uint8_t> remaining = donor_count_;

    topological_order_.reserve(n_valid_cells_);
    std::queue<std::int64_t> q;

    // Seed: valid cells with no donors (ridge / source cells)
    for (std::int64_t flat = 0; flat < n_cells_; ++flat)
        if (valid_[flat] != 0 && remaining[flat] == 0)
            q.push(flat);

    while (!q.empty())
    {
        const std::int64_t flat = q.front(); q.pop();
        topological_order_.push_back(flat);

        const auto begin = receiver_offsets_[flat];
        const auto end   = receiver_offsets_[flat + 1];
        for (auto e = begin; e < end; ++e)
        {
            const std::int64_t rec = receiver_indices_[e];
            if (--remaining[rec] == 0)
                q.push(rec);
        }
    }

    // Cycle guard
    std::int64_t expected_valid = 0;
    for (std::uint8_t v : valid_)
        if (v != 0)
            ++expected_valid;

    if (static_cast<std::int64_t>(topological_order_.size()) != expected_valid)
        throw std::runtime_error(
            "lsdtt_flowgraph_cpu: cycle detected in receiver graph or "
            "inconsistent cell validity. Ensure the DEM is pit-filled and "
            "that valid cells have finite, non-nodata elevations.");

}


// =============================================================================
// accumulate
// =============================================================================

std::vector<double> lsdtt_flowgraph_cpu::accumulate(
    const std::vector<double>& cell_weights) const
{
    const double default_w   = meta_.cell_size * meta_.cell_size;
    const bool   use_default = cell_weights.empty();

    if (!use_default &&
        static_cast<std::int64_t>(cell_weights.size()) != n_cells_)
    {
        throw std::invalid_argument(
            "lsdtt_flowgraph_cpu::accumulate: cell_weights size must equal "
            "n_rows * n_cols");
    }

    if (!use_default)
    {
        for (std::int64_t flat = 0; flat < n_cells_; ++flat)
        {
            if (valid_[flat] == 0) continue;

            if (!std::isfinite(cell_weights[flat]))
            {
                throw std::invalid_argument(
                    "lsdtt_flowgraph_cpu::accumulate: cell_weights must be finite "
                    "for all valid cells");
            }
        }
    }

    std::vector<double> accum(n_cells_, 0.0);

    // Initialise each valid cell with its own contribution
    for (std::int64_t flat : topological_order_)
        accum[flat] = use_default ? default_w : cell_weights[flat];


    // Forward push: each cell distributes its accumulated value to receivers
    for (std::int64_t flat : topological_order_)
    {
        const double a     = accum[flat];
        const auto   begin = receiver_offsets_[flat];
        const auto   end   = receiver_offsets_[flat + 1];
        for (auto e = begin; e < end; ++e)
            accum[receiver_indices_[e]] += a * static_cast<double>(receiver_weights_[e]);
    }

    return accum;
}


// =============================================================================
// accumulate_raster
// =============================================================================

lsdtt_raster_result lsdtt_flowgraph_cpu::accumulate_raster(
    const std::vector<double>& cell_weights) const
{
    const std::vector<double> accum = accumulate(cell_weights);
    const float ndv = static_cast<float>(meta_.no_data_value);

    lsdtt_raster_result result;
    result.meta = meta_;
    result.data.resize(n_cells_, ndv);
    result.mask.resize(n_cells_, 0);

    for (std::int64_t flat : topological_order_)
    {
        result.data[flat] = static_cast<float>(accum[flat]);
        result.mask[flat] = 1;
    }

    return result;
}

// =============================================================================
// get_drainage_area_raster_with_injection
// =============================================================================

lsdtt_raster_result lsdtt_flowgraph_cpu::get_drainage_area_raster_with_injection(
    const std::string& vector_path,
    const std::string& drainage_area_field) const
{
    const std::vector<lsdtt_drainage_injection_record> injections =
        load_drainage_injection_points(vector_path, meta_, drainage_area_field);

    return accumulate_raster(drainage_area_injection_weights(
        meta_, static_cast<std::size_t>(n_cells_), injections));
}

void lsdtt_flowgraph_cpu::write_drainage_area_raster_with_injection(
    const std::string& vector_path,
    const std::string& out_path,
    raster_output_format format,
    const std::string& drainage_area_field) const
{
    lsdtt_raster_result r =
        get_drainage_area_raster_with_injection(vector_path, drainage_area_field);
    lsdtt_data_loader loader;
    loader.write_full(r, out_path, format);
}

// =============================================================================
// outlet_raster
// =============================================================================

lsdtt_raster_result lsdtt_flowgraph_cpu::outlet_raster() const
{
    const float ndv = static_cast<float>(meta_.no_data_value);

    lsdtt_raster_result result;
    result.meta = meta_;
    result.data.resize(n_cells_, ndv);
    result.mask.resize(n_cells_, 0);

    for (std::int64_t flat = 0; flat < n_cells_; ++flat)
    {
        if (valid_[flat] == 0) continue;

        result.mask[flat] = 1;

        const bool is_outlet =
            (receiver_offsets_[flat + 1] == receiver_offsets_[flat]);

        result.data[flat] = is_outlet ? 1.0f : 0.0f;
    }

    return result;
}




// =============================================================================
// dinf_slopes
// =============================================================================

lsdtt_raster_result lsdtt_flowgraph_cpu::dinf_slopes() const
{
    const float ndv = static_cast<float>(meta_.no_data_value);

    lsdtt_raster_result result;
    result.meta = meta_;
    result.data.resize(n_cells_, ndv);
    result.mask.resize(n_cells_, 0);

    if (dinf_slopes_.empty()) return result;  // MFD mode: return all nodata

    for (std::int64_t flat = 0; flat < n_cells_; ++flat)
    {
        if (dinf_slopes_[flat] != ndv && valid_[flat] != 0)
        {
            result.data[flat] = dinf_slopes_[flat];
            result.mask[flat] = 1;
        }
    }

    return result;
}


// =============================================================================
// dinf_flow_direction
// =============================================================================

lsdtt_raster_result lsdtt_flowgraph_cpu::dinf_flow_direction() const
{
    const float ndv = static_cast<float>(meta_.no_data_value);

    lsdtt_raster_result result;
    result.meta = meta_;
    result.data.resize(n_cells_, ndv);
    result.mask.resize(n_cells_, 0);

    if (dinf_direction_.empty()) return result;  // MFD mode

    for (std::int64_t flat = 0; flat < n_cells_; ++flat)
    {
        if (dinf_direction_[flat] != ndv && valid_[flat] != 0)
        {
            result.data[flat] = dinf_direction_[flat];
            result.mask[flat] = 1;
        }
    }

    return result;
}

// -----------------------------------------------------------------------------
// build_reverse_levels
// -----------------------------------------------------------------------------
std::vector<std::vector<std::int64_t>>
build_reverse_levels(const lsdtt_flowgraph_cpu& flow)
{
    const auto& topo   = flow.topological_order();
    const std::int64_t N = static_cast<std::int64_t>(topo.size());

    std::vector<int> depth(static_cast<std::size_t>(
        flow.meta().n_rows * flow.meta().n_cols), -1);

    int max_depth = 0;
    for (std::int64_t i = N - 1; i >= 0; --i)
    {
        const std::int64_t flat = topo[static_cast<std::size_t>(i)];
        const std::int64_t nr   = flow.n_receivers(flat);

        if (nr == 0)
        {
            depth[static_cast<std::size_t>(flat)] = 0;
        }
        else
        {
            int d = 0;
            for (std::int64_t k = 0; k < nr; ++k)
            {
                const std::int64_t recv = flow.receiver_index_at(flat, k);
                if (recv >= 0 && depth[static_cast<std::size_t>(recv)] >= 0)
                    d = std::max(d, depth[static_cast<std::size_t>(recv)] + 1);
            }
            depth[static_cast<std::size_t>(flat)] = d;
            max_depth = std::max(max_depth, d);
        }
    }

    std::vector<std::vector<std::int64_t>> levels(
        static_cast<std::size_t>(max_depth + 1));

    for (std::int64_t i = N - 1; i >= 0; --i)
    {
        const std::int64_t flat = topo[static_cast<std::size_t>(i)];
        const int d = depth[static_cast<std::size_t>(flat)];
        if (d >= 0)
            levels[static_cast<std::size_t>(d)].push_back(flat);
    }

    return levels;
}

// -----------------------------------------------------------------------------
// build_forward_levels
// -----------------------------------------------------------------------------
std::vector<std::vector<std::int64_t>>
build_forward_levels(const lsdtt_flowgraph_cpu& flow)
{
    const auto& topo   = flow.topological_order();
    const std::int64_t N       = static_cast<std::int64_t>(topo.size());
    const std::int64_t n_cells = static_cast<std::int64_t>(
        flow.meta().n_rows) * flow.meta().n_cols;

    std::vector<int> depth(static_cast<std::size_t>(n_cells), -1);

    int max_depth = 0;
    for (std::int64_t i = 0; i < N; ++i)
    {
        const std::int64_t flat = topo[static_cast<std::size_t>(i)];
        if (depth[static_cast<std::size_t>(flat)] < 0)
            depth[static_cast<std::size_t>(flat)] = 0;

        const int d_next = depth[static_cast<std::size_t>(flat)] + 1;
        const std::int64_t nr = flow.n_receivers(flat);
        for (std::int64_t k = 0; k < nr; ++k)
        {
            const std::int64_t recv = flow.receiver_index_at(flat, k);
            if (recv >= 0 && recv < n_cells)
            {
                if (depth[static_cast<std::size_t>(recv)] < d_next)
                {
                    depth[static_cast<std::size_t>(recv)] = d_next;
                    max_depth = std::max(max_depth, d_next);
                }
            }
        }
    }

    std::vector<std::vector<std::int64_t>> levels(
        static_cast<std::size_t>(max_depth + 1));

    for (std::int64_t i = 0; i < N; ++i)
    {
        const std::int64_t flat = topo[static_cast<std::size_t>(i)];
        const int d = depth[static_cast<std::size_t>(flat)];
        if (d >= 0)
            levels[static_cast<std::size_t>(d)].push_back(flat);
    }

    return levels;
}

} // namespace lsdtt
