// lsdtt_giuh_cpu.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html

#include "lsdtt_giuh_cpu.hpp"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <stdexcept>

namespace lsdtt
{

// =============================================================================
// compute_giuh
// =============================================================================
lsdtt_giuh_result compute_giuh(
    const lsdtt_junction_network_cpu&        jn,
    const std::vector<lsdtt_sediment_reach>& reaches,
    const lsdtt_giuh_options&                opts)
{
    const int n_junctions = jn.n_junctions();
    if (n_junctions == 0)
        throw std::runtime_error("compute_giuh: junction network is empty");
    if (static_cast<int>(reaches.size()) != n_junctions)
        throw std::runtime_error(
            "compute_giuh: reaches.size() != jn.n_junctions()");

    const std::vector<std::int32_t>& topo   = jn.topological_order_vector();
    const std::vector<std::int32_t>& rec_jn = jn.receiver_junction_vector();

    // -------------------------------------------------------------------------
    // Step 1 — compute link-local travel time for each reach
    //   v_i = v_ref * (max(slope, min_slope) / slope_ref)^beta
    //   τ_link_i = length_i / v_i
    // -------------------------------------------------------------------------
    std::vector<double> link_tt(n_junctions, 0.0);
    double total_area = 0.0;

    for (int idx = 0; idx < n_junctions; ++idx)
    {
        const lsdtt_sediment_reach& r = reaches[idx];

        double slope = std::max(r.slope, opts.min_slope);
        double velocity = opts.v_ref_m_per_s *
            std::pow(slope / opts.slope_ref, opts.velocity_exponent);
        velocity = std::max(velocity, 1.0e-6);

        double tt = (r.length_m > 0.0) ? r.length_m / velocity : 0.0;

        if (opts.hillslope_length_m > 0.0)
            tt += opts.hillslope_length_m / opts.v_hillslope_m_per_s;

        link_tt[idx] = tt;
        total_area  += r.drainage_area_m2;
    }

    // -------------------------------------------------------------------------
    // Step 2 — accumulate travel time from each junction to the outlet
    //   Process in topological order (headwaters first).
    //   travel_time_to_outlet[j] = link_tt[j] + travel_time_to_outlet[receiver[j]]
    //   Outlet junctions (receiver == k_outlet_junction) get just link_tt[j].
    // -------------------------------------------------------------------------
    std::vector<double> tt_to_outlet(n_junctions, 0.0);

    // Process headwaters to outlets: each junction's cumulative time =
    // its own link time plus the receiver's cumulative time.
    // Because topo is headwater-first, the receiver is already set when we
    // reach it in reverse order (outlet-first pass).

    // First pass: headwaters → outlets, set each junction's cumulative tt.
    // We need outlet-to-headwater order to propagate, so reverse the topo.
    // Actually: topo[0] is a headwater. Its receiver has a larger topo index
    // (closer to outlet). So iterate in REVERSE topo order to propagate from
    // outlets toward headwaters, then forward to accumulate.

    // Simpler: iterate forward (headwater first). For each junction j,
    //   tt_to_outlet[j] = link_tt[j]      if receiver is k_outlet_junction
    //                   = link_tt[j] + tt_to_outlet[receiver[j]]  otherwise
    // But receiver[j] may not have been processed yet at forward pass time.
    // Instead, iterate in REVERSE topological order (outlet first):
    //   outlets: tt_to_outlet[j] = link_tt[j]
    //   others:  tt_to_outlet[j] = link_tt[j] + tt_to_outlet[receiver[j]]

    for (int t = static_cast<int>(topo.size()) - 1; t >= 0; --t)
    {
        int j = topo[t];
        int r = rec_jn[j];
        if (r == lsdtt::k_outlet_junction)
            tt_to_outlet[j] = link_tt[j];
        else
            tt_to_outlet[j] = link_tt[j] + tt_to_outlet[r];
    }

    // -------------------------------------------------------------------------
    // Step 3 — bin drainage areas by travel time → width function
    // -------------------------------------------------------------------------
    const int n_bins = opts.n_bins;
    const double dt  = opts.dt_bin_s;

    double max_tt = *std::max_element(tt_to_outlet.begin(), tt_to_outlet.end());
    // Extend range slightly so the last junction falls inside the last bin
    double range  = std::max(max_tt * 1.001, dt * n_bins);

    std::vector<double> width_fn(n_bins, 0.0);

    for (int j = 0; j < n_junctions; ++j)
    {
        double tt = tt_to_outlet[j];
        int bin   = static_cast<int>(tt / range * n_bins);
        if (bin >= n_bins) bin = n_bins - 1;
        if (bin < 0)       bin = 0;
        width_fn[bin] += reaches[j].drainage_area_m2;
    }

    // -------------------------------------------------------------------------
    // Step 4 — normalise to get the IUH [s⁻¹]
    //   integral of IUH * dt over all bins = 1
    //   IUH[k] = W[k] / (total_area * dt)
    // -------------------------------------------------------------------------
    lsdtt_giuh_result result;
    result.iuh.resize(n_bins, 0.0);
    result.bin_centre_s.resize(n_bins);
    result.travel_time_s   = tt_to_outlet;
    result.catchment_area_m2 = total_area;
    result.dt_bin_s        = dt;
    result.n_bins          = n_bins;

    const double bin_width_s = range / n_bins;
    double iuh_peak = 0.0;
    double time_of_peak = 0.0;

    for (int k = 0; k < n_bins; ++k)
    {
        result.bin_centre_s[k] = (k + 0.5) * bin_width_s;
        if (total_area > 0.0)
            result.iuh[k] = width_fn[k] / (total_area * bin_width_s);
        if (result.iuh[k] > iuh_peak)
        {
            iuh_peak      = result.iuh[k];
            time_of_peak  = result.bin_centre_s[k];
        }
    }

    result.iuh_peak       = iuh_peak;
    result.time_of_peak_s = time_of_peak;

    return result;
}

// =============================================================================
// giuh_discharge
// =============================================================================
std::vector<double> giuh_discharge(
    const lsdtt_giuh_result&   giuh,
    const std::vector<double>& rainfall_m_per_s,
    double                     dt_rainfall_s,
    double                     runoff_coefficient)
{
    if (giuh.iuh.empty())
        throw std::runtime_error("giuh_discharge: IUH is empty");
    if (rainfall_m_per_s.empty())
        return {};

    const int n_rain  = static_cast<int>(rainfall_m_per_s.size());
    const int n_iuh   = static_cast<int>(giuh.iuh.size());
    const int n_out   = n_rain + n_iuh - 1;

    const double A = giuh.catchment_area_m2;

    // Resample IUH to rainfall timestep using trapezoidal quadrature
    // IUH is in [s⁻¹]; multiply by A [m²] * R [m/s] * C * dt_rain [s] → m³/s
    // If dt_rainfall_s matches dt_bin_s, no resampling needed.
    // For simplicity: we directly use the IUH bins and assume dt_rain ~ dt_bin.
    // The user should set dt_bin_s ≈ dt_rainfall_s for accuracy.

    std::vector<double> Q(n_out, 0.0);
    const double bin_width_s = giuh.dt_bin_s;

    for (int t = 0; t < n_rain; ++t)
    {
        double r_eff = rainfall_m_per_s[t] * runoff_coefficient;
        if (r_eff <= 0.0) continue;

        // How many IUH bins does one rainfall step cover?
        // Simple approach: spread rainfall pulse over the IUH.
        // Q contribution at output step t + k:
        //   A [m²] * r_eff [m/s] * dt_rain [s] * IUH[k] [s⁻¹] * bin_width_s
        // = A * r_eff * dt_rain * IUH[k] * bin_width_s
        // But IUH[k] [s⁻¹] * bin_width_s is dimensionless fraction.
        // Output Q has units m³/s, so divide by dt_rainfall_s at output step.
        // Equivalently: Q[t+k] += A * r_eff * IUH[k] * bin_width_s
        //                           (volume in m³ spread over dt_rain, giving m³/s
        //                            when IUH is per-second rate)
        // The convolution:
        //   Q[t_out] = A * sum_k R[t_out - k] * C * iuh[k] * bin_width_s * dt_rain
        // Here dt_rain is absorbed into R (R already in m/s, so rainfall depth
        // per step = R * dt_rain). The IUH is [s⁻¹], and integrating over
        // bin_width_s [s] gives a fraction. The result is in m³/s.

        for (int k = 0; k < n_iuh; ++k)
        {
            int out_idx = t + k;
            if (out_idx >= n_out) break;
            // Volume of runoff from this rain step: A * r_eff * dt_rain_s [m³]
            // Distributed over output step k: multiply by IUH[k] * bin_width_s (fraction)
            // Divide by dt_rainfall_s to get m³/s
            Q[out_idx] += A * r_eff * giuh.iuh[k] * bin_width_s;
        }
    }

    return Q;
}

// =============================================================================
// giuh_pulse_discharge
// =============================================================================
std::vector<double> giuh_pulse_discharge(
    const lsdtt_giuh_result& giuh,
    double pulse_duration_s,
    double rainfall_rate_m_per_s,
    double dt_s,
    double runoff_coefficient)
{
    if (dt_s <= 0.0)
        throw std::runtime_error("giuh_pulse_discharge: dt_s must be > 0");

    int n_rain_steps = std::max(1, static_cast<int>(
        std::ceil(pulse_duration_s / dt_s)));

    std::vector<double> rain(n_rain_steps, rainfall_rate_m_per_s);

    return giuh_discharge(giuh, rain, dt_s, runoff_coefficient);
}

} // namespace lsdtt
