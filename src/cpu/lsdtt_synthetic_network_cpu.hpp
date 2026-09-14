// lsdtt_synthetic_network_cpu.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// lsdtt_synthetic_network_cpu
//
// Generates synthetic river networks with TUNABLE branching statistics
// (in particular tunable Tokunaga exponents) using the Optimal Channel Network
// (OCN) construction: a spanning tree of D8 flow directions over a rectangular
// lattice is evolved by simulated-annealing edge flips that minimise the energy
//
//     E = sum_i  A_i^gamma * L_i    (A_i = drainage area, L_i = link length)
//
// L_i is the length of the link from i to its receiver: cell_size for a cardinal
// step, cell_size*sqrt(2) for a diagonal one. This is the DEFAULT and it is the
// physically correct form -- link energy expenditure is E_i ~ L_i Q_i S_i
// (Strong & Mudd 2022 PNAS eq. 4; Goren & Shelef 2024 ESurf 12:1347 eq. 5).
// Setting energy_length_weighted = false drops L_i and reproduces OCNet exactly.
// Do not drop it casually: on an 8-connected lattice it is what stops the
// annealer buying free lateral movement, and without it the networks come out
// almost isotropic (see the option's documentation for the measurement).
//
// A candidate reconnection is accepted with the Metropolis rule
//     accept iff  E_new <= E_prev  OR  rand() < exp(-(E_new - E_prev)/T),
// so the cooling schedule (temperature) sweeps the network from random /
// Shreve-like (high T, Tokunaga c ~ 2) to strongly aggregated (T -> 0).
//
// Elevations are derived AFTER the network is fixed by the chi-based fluvial
// steady state, via the SAME shared core (compute_fluvial_steady_state_on_tree)
// that lsdtt-lem's steady-state solver uses: chi[i] = chi[receiver] + dx_eff *
// A^(-m/n) and Z_i = baselevel + f*chi_i. Relief is controlled exactly as in the
// LEM (K-mode or relief-mode). Fluvial relief only - so an OCN and a LEM steady
// state on the same graph and parameters produce the same surface.
//
// IMPORTANT round-trip caveat: the chi DEM does NOT exactly reproduce the OCN
// graph under independent D8 routing, at ANY temperature. A greedy local
// minimum of sum(A^gamma) under single-flip moves is not identically the
// steepest-descent graph of the chi surface, so some interior pixels always
// re-route: a pixel is guaranteed to be strictly higher than its OCN receiver
// (that invariant does hold everywhere), but a different neighbour may be lower
// still and D8 will then take it.
//
// Measured 2026-08-13 (48x48, gamma 0.5, seed 42, default schedule
// n_iter = 40*n_rows*n_cols): fidelity is 0.916352 at temperature_scale = 0 and
// 0.845936 at temperature_scale = 1. Relaxing further does not close the gap -
// fidelity PLATEAUS near 0.935 and never reaches 1, even at 400000 iterations.
// Anneal LENGTH dominates temperature: the same configuration truncated to
// n_iter = 20000 gives only ~0.703 at temperature_scale = 0, so any comparison
// across temperatures must hold n_iter fixed. Low temperature helps at the
// endpoints but fidelity is not monotone in the interior of the temperature
// range. The CLI reports the realised fidelity; unit_synthetic_networks pins a
// floor on it (see tests/cpp/test_synthetic_networks.cpp).
//
// Design reference: Carraro et al. (2020) OCNet; Rodriguez-Iturbe & Rinaldo
// (1997); see ../OCNet for the reference R package.
// -----------------------------------------------------------------------------

#pragma once

#include <cstdint>
#include <random>
#include <vector>

#include "../core/lsdtt_raster_meta.hpp"
#include "../core/lsdtt_raster_result.hpp"

namespace lsdtt {

// Where drainage leaves the domain.
//   single_south : one outlet (rectangle: south-edge centre; island: southmost
//                  land pixel) - the whole domain drains through one network.
//   all_coast    : every coastal land pixel (land touching sea or the grid edge)
//                  is a base level - an island draining to its whole coastline.
// single_south : one outlet on the south edge; the whole domain is one network.
// all_coast    : every coastal land pixel is a base level.
// north_south  : every land pixel on the north AND south edges is a base level,
//                and nothing else is. Paired with periodic_columns this is the
//                classic ridge-and-valley strip -- drainage leaves through two
//                opposed boundaries and the main divide is free to sit anywhere
//                between them, which is what makes divide position a meaningful
//                measurement rather than an artefact of a single outlet.
enum class ocn_outlet_mode { single_south, all_coast, north_south };

struct lsdtt_ocn_options
{
    int       n_rows    = 128;
    int       n_cols    = 128;
    double    cell_size = 1.0;

    // Optional land/sea mask (row-major, size n_rows*n_cols, 1 = land). Empty
    // means the full rectangle is land. Pixels outside the mask are excluded
    // from the network and written as nodata.
    std::vector<std::uint8_t> domain_mask;
    ocn_outlet_mode           outlet_mode = ocn_outlet_mode::single_south;

    // Column spacing of baselevel outlets in north_south mode. 1 (the default)
    // makes every boundary pixel an outlet, so basin width is set by the grid,
    // not by the optimisation. Larger values force wider basins -- the OCN
    // analogue of the diamond-square feature wavelength. Ignored in the other
    // outlet modes.
    int outlet_spacing = 1;

    // Wrap the lattice in the column direction, so column 0 and column n_cols-1
    // are D8 neighbours. Must match the 'p' flags in the boundary_conditions
    // used to route the resulting DEM: a network annealed with wrapping and then
    // routed without it (or the reverse) is not the network that was generated.
    bool periodic_columns = false;

    // Weight each node's energy contribution by the LENGTH of its downstream
    // link:  E = sum_i A_i^gamma * L_i   instead of   E = sum_i A_i^gamma.
    //
    // This is the physically correct discretisation and it is the DEFAULT.
    // The energy expenditure of a link is E_i ~ L_i Q_i S_i (Strong & Mudd
    // 2022, PNAS, eq. 4; Goren & Shelef 2024, Earth Surf. Dynam. 12, 1347-1369,
    // eq. 5), which reduces to A_i^gamma * L_i under the usual hydraulic
    // geometry scalings. On a 4-connected lattice every L_i is equal, so it
    // factors out and the literature writes it as sum(A^gamma) -- but on the
    // 8-connected lattice used here a diagonal link is sqrt(2) longer than a
    // cardinal one, and dropping L_i under-charges it by that factor.
    //
    // The consequence of dropping it is not subtle. Lateral movement becomes
    // free, so the annealer spends it aggregating sideways and the network
    // comes out almost isotropic: measured on a 400x200 strip draining north
    // and south, 45.2% of links were diagonal (50% is random) and 26.8% ran
    // pure east-west, parallel to the divide and making no progress toward any
    // outlet. Goren & Shelef make the same point from the other side: with
    // uniform l_i the energy reduces to a function of the number of nodes.
    //
    // Set false to reproduce OCNet (Carraro et al. 2020), which sums A^gamma
    // unweighted over an 8-neighbour stencil. That is a defensible choice if
    // the object of study is the flow GRAPH rather than a geometric network,
    // and it is what the regression oracle pins -- but it is not the right
    // default for a tool whose output is a DEM.
    bool energy_length_weighted = true;

    // Optional georeference override. When use_supplied_meta is true the output
    // rasters use supplied_meta verbatim (its n_rows/n_cols/cell_size must match
    // n_rows/n_cols/cell_size above) instead of the default Watopia anchor -
    // used by true-area polygon import to preserve the polygon's real extent/CRS.
    bool              use_supplied_meta = false;
    lsdtt_raster_meta supplied_meta;

    // Optional synthetic georeference anchor. Ignored when use_supplied_meta is
    // true.
    bool   use_geographic_anchor = false;
    double centre_latitude_deg   = 0.0;
    double centre_longitude_deg  = 0.0;
    int    projection_epsg       = 0;

    double    gamma     = 0.5;     // energy exponent (OCNet expEnergy)

    // Annealing schedule. n_iter <= 0 means 40 * n_rows * n_cols (OCNet default).
    long long n_iter                     = -1;
    double    cooling_rate               = 1.0;  // larger = FASTER cooling (OCNet)
    double    initial_no_cooling_fraction = 0.0;

    // Initial temperature. If initial_temperature > 0 it is used directly;
    // otherwise T0 = temperature_scale * E0 (E0 = energy of the initial tree;
    // temperature_scale = 1 reproduces OCNet's native T0 = E0).
    double    initial_temperature = 0.0;
    double    temperature_scale   = 1.0;

    // Landscape (elevation) synthesis - chi-based fluvial steady state, identical
    // to lsdtt-lem's steady-state solver (compute_fluvial_steady_state_on_tree):
    //   chi[i] = chi[receiver] + dx_eff * A^(-m/n);  z = baselevel + f*chi.
    // OCN elevations use FLUVIAL relief only (no hillslope), so the OCN and a LEM
    // steady state on the same graph and parameters produce the same surface.
    float     m_exponent          = 0.45f;   // drainage-area exponent m
    float     n_exponent          = 1.0f;    // slope exponent n
    bool      relief_mode         = true;    // true: set total relief; false: K-mode
    float     target_relief       = 1000.0f; // relief-mode target (m)
    float     uplift_rate         = 1e-4f;   // U (m/yr)
    float     erodibility         = 1e-3f;   // K (K-mode)
    float     baselevel_elevation = 0.0f;    // elevation of the outlet(s)

    std::uint64_t seed = 1;

    // Outlet on the south (bottom) edge. outlet_col < 0 => centre column.
    int       outlet_col = -1;
};

class lsdtt_synthetic_network_cpu
{
public:
    explicit lsdtt_synthetic_network_cpu(const lsdtt_ocn_options& opts);

    // Run the simulated-annealing optimisation.
    void run();

    // ---- outputs (valid after run()) ----
    [[nodiscard]] lsdtt_raster_result get_elevation_raster() const;
    [[nodiscard]] lsdtt_raster_result get_drainage_area_cells_raster() const;

    [[nodiscard]] double initial_energy() const noexcept { return energy0_; }
    [[nodiscard]] double final_energy()   const noexcept { return energy_; }

    // Receiver flat index per node (outlet maps to itself); for verification.
    [[nodiscard]] const std::vector<std::int64_t>& receivers() const noexcept
    { return receiver_; }
    [[nodiscard]] const std::vector<double>& areas() const noexcept
    { return area_; }

    [[nodiscard]] int n_rows() const noexcept { return nrows_; }
    [[nodiscard]] int n_cols() const noexcept { return ncols_; }
    [[nodiscard]] std::int64_t outlet_flat() const noexcept { return outlet_; }
    [[nodiscard]] bool is_outlet(std::int64_t f) const noexcept
    { return f >= 0 && f < n_ && is_outlet_[static_cast<std::size_t>(f)] != 0; }
    [[nodiscard]] bool is_land(std::int64_t f) const noexcept
    { return f >= 0 && f < n_ && land_[static_cast<std::size_t>(f)] != 0; }

private:
    // Column index after applying the lateral boundary. Returns -1 for a column
    // that is genuinely off-grid, so callers can keep a single rejection test:
    // with periodic_columns the index wraps and is never -1, without it the
    // out-of-range columns are rejected exactly as before.
    [[nodiscard]] double link_length_(std::int64_t node, std::int64_t recv) const;

    [[nodiscard]] int wrap_col(int c) const noexcept
    {
        if (c >= 0 && c < ncols_) return c;
        if (!opts_.periodic_columns) return -1;
        // ncols_ >= 2 is enforced in the constructor, and neighbour offsets are
        // +/-1, so a single add/subtract is sufficient here.
        if (c < 0)       return c + ncols_;
        return c - ncols_;
    }

    lsdtt_ocn_options opts_;
    int          nrows_ = 0;
    int          ncols_ = 0;
    std::int64_t n_     = 0;
    std::int64_t outlet_ = 0;        // a representative outlet (outlets_[0])

    std::vector<std::uint8_t> land_;      // 1 = land (in domain)
    std::vector<char>         is_outlet_; // 1 = base-level outlet
    std::vector<std::int64_t> outlets_;   // all outlet flats
    std::vector<std::int64_t> receiver_;  // receiver flat index (outlet -> itself; sea -> -1)
    std::vector<double>       area_;      // accumulated area in cells (>= 1 on land)
    std::vector<std::int64_t> topo_order_;// upstream-first topological order

    double energy_  = 0.0;
    double energy0_ = 0.0;

    std::mt19937_64 rng_;

    // scratch for the incremental delta-E / commit
    std::vector<std::uint32_t> mark_;
    std::uint32_t              stamp_ = 0;
    std::vector<std::int64_t>  losers_;
    std::vector<std::int64_t>  gainers_;

    lsdtt_raster_meta build_meta_() const;
    void   setup_domain_and_outlets_();   // fills land_, is_outlet_, outlets_
    void   init_forest_();                // multi-source BFS spanning forest
    void   recompute_all_areas_();        // also fills topo_order_
    double total_energy_() const;
    double receiver_length_(std::int64_t node) const;

    [[nodiscard]] bool creates_cycle_(std::int64_t node,
                                      std::int64_t cand_receiver) const;

    // Computes delta-E for reconnecting node -> cand_receiver and fills
    // losers_/gainers_. Does NOT modify state.
    double delta_energy_(std::int64_t node, std::int64_t cand_receiver);
    void   commit_(std::int64_t node, std::int64_t cand_receiver, double dE);
};

} // namespace lsdtt
