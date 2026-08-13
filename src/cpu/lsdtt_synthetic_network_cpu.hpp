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
// Elevations are derived AFTER the network is fixed by the chi-based fluvial
// steady state, via the SAME shared core (compute_fluvial_steady_state_on_tree)
// that lsdtt-lem's steady-state solver uses: chi[i] = chi[receiver] + dx_eff *
// A^(-m/n) and Z_i = baselevel + f*chi_i. Relief is controlled exactly as in the
// LEM (K-mode or relief-mode). Fluvial relief only - so an OCN and a LEM steady
// state on the same graph and parameters produce the same surface.
//
// IMPORTANT round-trip caveat: the chi DEM reproduces the OCN graph under
// independent D8 routing only for fully-relaxed (low-T) networks, where the OCN
// local-minimum condition makes each node's receiver its steepest-descent
// neighbour. High-T / sub-optimal networks are NOT self-consistent and D8 will
// not exactly reproduce them. The CLI reports the realised fidelity.
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
enum class ocn_outlet_mode { single_south, all_coast };

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
