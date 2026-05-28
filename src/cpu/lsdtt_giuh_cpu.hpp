// lsdtt_giuh_cpu.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html
//
// -----------------------------------------------------------------------------
// Geomorphic Instantaneous Unit Hydrograph (GIUH)
//
// Computes the width-function-based IUH of a D8 channel network and generates
// a discharge time series from a rainfall input.
//
// Theory
// ------
// The approach follows Rodriguez-Iturbe & Valdes (1979) and Mesa & Mifflin
// (1986). Each reach (junction link) contributes drainage area that must travel
// a certain time to the network outlet. Aggregating area by travel time gives
// the geomorphic width function W(τ). Normalised by total catchment area, W(τ)
// is the IUH: the discharge response at the outlet to a unit rainfall applied
// instantaneously over the whole catchment.
//
// Travel time per link is estimated as:
//   v_i = v_ref * (S_i / S_ref)^β   (kinematic-wave celerity)
//   τ_i = L_i / v_i
//
// where S_i is reach slope, L_i reach length (from build_sediment_reaches),
// v_ref and S_ref are reference velocity and slope, and β = 0.5 by default.
//
// Overland-flow travel time (hillslope contribution) is optionally added to
// each link's travel time as τ_hf = hillslope_length / v_hillslope.
//
// The IUH is then binned into n_bins time steps of width dt_bin_s.
//
// Discharge generation
// --------------------
// Given a rainfall time series R[t] [m/s] and runoff coefficient C [-],
// the discharge at the outlet for each output timestep t_out is:
//
//   Q[t_out] = A_catchment * sum_k ( R[t_out - k] * C * IUH[k] * dt_rain )
//
// where IUH[k] is the fraction of catchment area contributing in bin k.
//
// References
// ----------
//   Rodriguez-Iturbe I. & Valdes J.B. (1979) WRR 15(6):1409-1420
//   Mesa O.J. & Mifflin E.R. (1986) Scale problems in hydrology, Reidel.
// -----------------------------------------------------------------------------

#pragma once

#include <cstdint>
#include <string>
#include <vector>

#include "lsdtt_flowinfo_cpu.hpp"
#include "lsdtt_junction_network_cpu.hpp"
#include "lsdtt_raster_cpu.hpp"
#include "lsdtt_sediment_transport_cpu.hpp"

namespace lsdtt
{

// ---------------------------------------------------------------------------
// Options for the GIUH computation
// ---------------------------------------------------------------------------
struct lsdtt_giuh_options
{
    // Kinematic-wave velocity model: v = v_ref * (slope / slope_ref)^velocity_exponent
    double v_ref_m_per_s   = 1.0;    // reference channel velocity [m/s]
    double slope_ref       = 0.01;   // reference slope for v_ref
    double velocity_exponent = 0.5;  // exponent β; 0.5 = Manning-like

    // Hillslope travel time added to each link
    // τ_hillslope = hillslope_length_m / v_hillslope_m_per_s
    // Set hillslope_length_m <= 0 to disable.
    double hillslope_length_m       = -1.0;
    double v_hillslope_m_per_s      = 0.005; // typical overland-flow velocity

    // IUH discretisation
    int    n_bins           = 100;    // number of time bins in the IUH
    double dt_bin_s         = 3600.0; // bin width [s] (1 hour default)

    // Minimum slope used in velocity estimation (avoids zero velocity)
    double min_slope        = 1.0e-5;

    // Runoff coefficient applied to rainfall to get effective rainfall
    double runoff_coefficient = 1.0;
};

// ---------------------------------------------------------------------------
// The computed IUH — one value per time bin
// ---------------------------------------------------------------------------
struct lsdtt_giuh_result
{
    // iuh[k] = fraction of catchment area [s⁻¹] contributing in bin k
    // Integral of iuh over all bins * dt_bin_s = 1.0
    std::vector<double> iuh;

    // Travel time at the centre of each bin [s]
    std::vector<double> bin_centre_s;

    // Per-junction travel time to the outlet [s]
    std::vector<double> travel_time_s;

    // Total catchment area used in the computation [m²]
    double catchment_area_m2 = 0.0;

    double dt_bin_s = 3600.0;
    int    n_bins   = 0;

    // Peak of the IUH [s⁻¹] and time of peak [s]
    double iuh_peak         = 0.0;
    double time_of_peak_s   = 0.0;
};

// ---------------------------------------------------------------------------
// Compute the width-function GIUH from precomputed reach geometry.
// reaches must have been produced by build_sediment_reaches().
// ---------------------------------------------------------------------------
[[nodiscard]]
lsdtt_giuh_result compute_giuh(
    const lsdtt_junction_network_cpu&          jn,
    const std::vector<lsdtt_sediment_reach>&   reaches,
    const lsdtt_giuh_options&                  opts);

// ---------------------------------------------------------------------------
// Generate a discharge time series from a rainfall input and the GIUH.
//
// rainfall_m_per_s[t] — rainfall rate at output timestep t [m/s]
// dt_rainfall_s       — duration of each rainfall timestep [s]
//
// Returns Q[t] [m³/s] at the network outlet for each rainfall timestep.
// The output vector has the same length as rainfall_m_per_s plus
// (giuh.n_bins - 1) tail steps so the complete recession is captured.
// ---------------------------------------------------------------------------
[[nodiscard]]
std::vector<double> giuh_discharge(
    const lsdtt_giuh_result&      giuh,
    const std::vector<double>&    rainfall_m_per_s,
    double                        dt_rainfall_s,
    double                        runoff_coefficient = 1.0);

// ---------------------------------------------------------------------------
// Convenience: build a synthetic uniform rainfall pulse of duration
// pulse_duration_s at rate rainfall_rate_m_per_s, then compute discharge.
// dt_s is both the rainfall timestep and the output timestep.
// ---------------------------------------------------------------------------
[[nodiscard]]
std::vector<double> giuh_pulse_discharge(
    const lsdtt_giuh_result& giuh,
    double pulse_duration_s,
    double rainfall_rate_m_per_s,
    double dt_s,
    double runoff_coefficient = 1.0);

} // namespace lsdtt
