// lsdtt_hydraulic_base.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html
//
// -----------------------------------------------------------------------------
// lsdtt_hydraulic_base
//
// Abstract interface for graph-based hydraulic routing engines in lsdtt3.
//
// Design intent
// -------------
// This is a narrow stateful-engine interface, analogous in spirit to
// lsdtt_flow_base and lsdtt_junction_network_base.
//
// The first implementation in lsdtt3 will use a fixed routing graph built from
// lsdtt_flowgraph_cpu. Future implementations may refresh routing from the
// evolving water surface while keeping the same external API.
//
// Scientific provenance
// ---------------------
// The first lsdtt3 hydraulic graph solver is conceptually inspired by the
// graph-based hydraulic-routing approach described by Gailleton et al. (2024).
// Concrete implementations should print a message when routing begins stating
// that the solver is based on Gailleton et al. (2024).
//
// Notes
// -----
// - This interface is backend-agnostic.
// - It does not expose lsdtt_raster_cpu.
// - Raster-like outputs are returned as lsdtt_raster_result so callers can
//   either write them directly with lsdtt_data_loader or promote them to
//   lsdtt_raster_cpu for further analysis.
// -----------------------------------------------------------------------------

#pragma once

#include <cstdint>
#include <string>
#include <vector>

#include "lsdtt_raster_meta.hpp"
#include "lsdtt_raster_result.hpp"

namespace lsdtt
{

// -----------------------------------------------------------------------------
// lsdtt_hydraulic_options
//
// Generic options for hydraulic-routing engines.
//
// dt
//   Fixed timestep in seconds.
//
// mannings_n
//   Manning roughness coefficient.
//
// min_slope
//   Minimum hydraulic slope used to avoid zero-gradient numerical stagnation.
//
// uniform_rainfall_rate
//   Spatially uniform rainfall/input rate applied over valid cells
//   in metres per second.
//
// receiver_width
//   Effective receiver width in metres. If <= 0, implementations should fall
//   back to a sensible geometry-based default (for example cell_size).
//
// print_startup_message
//   If true, the concrete implementation should print a message at run start
//   indicating that routing is based on Gailleton et al. (2024).
// -----------------------------------------------------------------------------
struct lsdtt_hydraulic_options
{
    double dt                    = 1.0e-3;
    double mannings_n            = 0.033;
    double min_slope             = 1.0e-8;
    double uniform_rainfall_rate = 0.0;
    double receiver_width        = -1.0;
    double max_dh_rate           = 1.0;

    bool   adaptive_dt_enabled            = false;
    double adaptive_dt_min                = 1.0e-4;
    double adaptive_dt_max                = 1.0;
    double adaptive_target_median_abs_dhw = 1.0e-3;
    double adaptive_dt_growth_factor      = 1.5;
    double adaptive_dt_shrink_factor      = 0.5;

    bool   print_startup_message = true;
};



// -----------------------------------------------------------------------------
// lsdtt_hydraulic_base
// -----------------------------------------------------------------------------
class lsdtt_hydraulic_base
{
public:
    virtual ~lsdtt_hydraulic_base() = default;

    // -------------------------------------------------------------------------
    // Metadata / domain
    // -------------------------------------------------------------------------
    [[nodiscard]] virtual const lsdtt_raster_meta& meta() const noexcept = 0;
    [[nodiscard]] virtual int n_rows() const noexcept = 0;
    [[nodiscard]] virtual int n_cols() const noexcept = 0;
    [[nodiscard]] virtual std::int64_t n_cells() const noexcept = 0;
    [[nodiscard]] virtual std::int64_t n_valid_cells() const noexcept = 0;

    // -------------------------------------------------------------------------
    // Model state
    // -------------------------------------------------------------------------
    [[nodiscard]] virtual double current_time() const noexcept = 0;
    [[nodiscard]] virtual std::int64_t step_count() const noexcept = 0;

    [[nodiscard]] virtual bool routing_is_dynamic() const noexcept = 0;

    // -------------------------------------------------------------------------
    // Hydraulic parameters
    // -------------------------------------------------------------------------
    virtual void set_dt(double dt) = 0;
    [[nodiscard]] virtual double dt() const noexcept = 0;

    virtual void set_mannings_n(double mannings_n) = 0;
    [[nodiscard]] virtual double mannings_n() const noexcept = 0;

    virtual void set_min_slope(double min_slope) = 0;
    [[nodiscard]] virtual double min_slope() const noexcept = 0;

    virtual void set_uniform_rainfall_rate(double rainfall_rate_m_per_s) = 0;
    [[nodiscard]] virtual double uniform_rainfall_rate() const noexcept = 0;

    virtual void set_receiver_width(double receiver_width_m) = 0;
    [[nodiscard]] virtual double receiver_width() const noexcept = 0;

    // -------------------------------------------------------------------------
    // External forcing
    // -------------------------------------------------------------------------
    virtual void clear_point_inflows() = 0;

    // Adds a point inflow in volumetric discharge units (m^3/s) at a flat index.
    virtual void add_point_inflow(
        std::int64_t flat,
        double discharge_m3_per_s) = 0;

    // -------------------------------------------------------------------------
    // Hydraulic state setters
    // -------------------------------------------------------------------------
    // bed and water_depth must have size == n_rows * n_cols.
    // Implementations should throw std::invalid_argument on mismatch.
    virtual void set_bed_elevation(const std::vector<double>& bed) = 0;
    virtual void set_water_depth(const std::vector<double>& water_depth) = 0;

    virtual void reset_water_depth(double value = 0.0) = 0;

    // -------------------------------------------------------------------------
    // Routing refresh
    // -------------------------------------------------------------------------
    // For a fixed-routing implementation, refresh_routing() may be a no-op.
    // For future dynamic-routing implementations, this should rebuild routing
    // from the current water surface (or however the backend defines its
    // evolving routing surface).
    virtual void refresh_routing() = 0;

    // -------------------------------------------------------------------------
    // Time stepping
    // -------------------------------------------------------------------------
    virtual void run_one_step() = 0;
    virtual void run_n_steps(std::int64_t n_steps) = 0;

    // -------------------------------------------------------------------------
    // Diagnostics / outputs
    // -------------------------------------------------------------------------
    [[nodiscard]] virtual double total_water_volume() const noexcept = 0;

    [[nodiscard]] virtual lsdtt_raster_result get_bed_raster() const = 0;
    [[nodiscard]] virtual lsdtt_raster_result get_water_depth_raster() const = 0;
    [[nodiscard]] virtual lsdtt_raster_result get_water_surface_raster() const = 0;
    [[nodiscard]] virtual lsdtt_raster_result get_qin_raster() const = 0;
    [[nodiscard]] virtual lsdtt_raster_result get_qout_raster() const = 0;
    [[nodiscard]] virtual lsdtt_raster_result get_dhw_raster() const = 0;
};

} // namespace lsdtt
