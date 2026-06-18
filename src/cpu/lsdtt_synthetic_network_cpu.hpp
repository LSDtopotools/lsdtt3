// lsdtt_synthetic_network_cpu.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html
//
// -----------------------------------------------------------------------------
// lsdtt_synthetic_network_cpu
//
// Generates synthetic river networks with TUNABLE branching statistics
// (in particular tunable Tokunaga exponents) using the Optimal Channel Network
// (OCN) construction: a spanning tree of D8 flow directions over a rectangular
// lattice is evolved by simulated-annealing edge flips that minimise the energy
//
//     E = sum_i  A_i^gamma          (A_i = drainage area at pixel i)
//
// A candidate reconnection is accepted with the Metropolis rule
//     accept iff  E_new <= E_prev  OR  rand() < exp(-(E_new - E_prev)/T),
// so the cooling schedule (temperature) sweeps the network from random /
// Shreve-like (high T, Tokunaga c ~ 2) to strongly aggregated (T -> 0).
//
// Elevations are derived AFTER the network is fixed by upstream integration of a
// slope-area law (landscape_OCN analogue): S_i = slope0 * (A_i/A_total)^(gamma-1)
// and Z_i = Z_receiver(i) + S_i * length_i.
//
// IMPORTANT round-trip caveat: the slope-area DEM reproduces the OCN graph under
// independent D8 routing only for fully-relaxed (low-T) networks, where the OCN
// local-minimum condition makes each node's receiver its steepest-descent
// neighbour. High-T / sub-optimal networks are NOT self-consistent and D8 will
// not exactly reproduce them. The CLI reports the realised fidelity.
//
// Design reference: Carraro et al. (2020) OCNet; Rodriguez-Iturbe & Rinaldo
// (1997). See docs/synthetic_networks_plan.md (and ../OCNet for the R package).
// -----------------------------------------------------------------------------

#pragma once

#include <cstdint>
#include <random>
#include <vector>

#include "../core/lsdtt_raster_meta.hpp"
#include "../core/lsdtt_raster_result.hpp"

namespace lsdtt {

struct lsdtt_ocn_options
{
    int       n_rows    = 128;
    int       n_cols    = 128;
    double    cell_size = 1.0;

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

    // Landscape (elevation) synthesis.
    double    slope0  = 1.0;
    double    z_min   = 0.0;

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

private:
    lsdtt_ocn_options opts_;
    int          nrows_ = 0;
    int          ncols_ = 0;
    std::int64_t n_     = 0;
    std::int64_t outlet_ = 0;

    std::vector<std::int64_t> receiver_;  // receiver flat index (outlet -> itself)
    std::vector<double>       area_;      // accumulated area in cells (>= 1)
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
    void   init_comb_tree_();
    void   recompute_all_areas_();   // also fills topo_order_
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
