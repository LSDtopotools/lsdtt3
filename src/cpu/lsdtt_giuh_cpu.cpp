// lsdtt_giuh_cpu.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html

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
    // Step 1 - compute link-local travel time for each reach
    //   v_i = v_ref * (max(slope, min_slope) / slope_ref)^beta
    //   tau_link_i = length_i / v_i
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
    // Step 2 - accumulate travel time from each junction to the outlet
    //   Process in topological order (headwaters first).
    //   travel_time_to_outlet[j] = link_tt[j] + travel_time_to_outlet[receiver[j]]
    //   Outlet junctions (receiver == k_outlet_junction) get just link_tt[j].
    // -------------------------------------------------------------------------
    std::vector<double> tt_to_outlet(n_junctions, 0.0);

    // Process headwaters to outlets: each junction's cumulative time =
    // its own link time plus the receiver's cumulative time.
    // Because topo is headwater-first, the receiver is already set when we
    // reach it in reverse order (outlet-first pass).

    // First pass: headwaters -> outlets, set each junction's cumulative tt.
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
    // Step 3 - bin drainage areas by travel time -> width function
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
    // Step 4 - normalise to get the IUH [s^-^1]
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
// ---------------------------------------------------------------------------
// resample_iuh_density
//
// Re-bins a piecewise-constant IUH density from bins of width src_dt to bins of
// width dst_dt, conserving the integral exactly.
//
// The IUH is a density in [s^-1] normalised so that sum(iuh[k]) * src_dt == 1.
// Each source bin therefore carries a dimensionless mass iuh[k] * src_dt. That
// mass is split across the destination bins it overlaps, in proportion to the
// overlap, and each destination bin's accumulated mass is divided by dst_dt to
// return a density. The result satisfies sum(out[j]) * dst_dt == 1 to within
// rounding, whatever the ratio of the two bin widths.
// ---------------------------------------------------------------------------
static std::vector<double> resample_iuh_density(
    const std::vector<double>& iuh,
    double                     src_dt,
    double                     dst_dt)
{
    const int n_src = static_cast<int>(iuh.size());
    if (n_src == 0 || src_dt <= 0.0 || dst_dt <= 0.0) return {};

    const double total_s = static_cast<double>(n_src) * src_dt;
    const int    n_dst   = std::max(1, static_cast<int>(std::ceil(total_s / dst_dt)));

    std::vector<double> mass(static_cast<std::size_t>(n_dst), 0.0);

    for (int k = 0; k < n_src; ++k)
    {
        const double src_lo = static_cast<double>(k) * src_dt;
        const double src_hi = src_lo + src_dt;
        const double density = iuh[static_cast<std::size_t>(k)];
        if (density == 0.0) continue;

        int j_lo = static_cast<int>(std::floor(src_lo / dst_dt));
        int j_hi = static_cast<int>(std::ceil (src_hi / dst_dt));
        if (j_lo < 0) j_lo = 0;
        if (j_hi > n_dst) j_hi = n_dst;

        for (int j = j_lo; j < j_hi; ++j)
        {
            const double dst_lo = static_cast<double>(j) * dst_dt;
            const double dst_hi = dst_lo + dst_dt;
            const double overlap =
                std::min(src_hi, dst_hi) - std::max(src_lo, dst_lo);
            if (overlap > 0.0)
                mass[static_cast<std::size_t>(j)] += density * overlap;
        }
    }

    for (double& m : mass) m /= dst_dt;
    return mass;
}

std::vector<double> giuh_discharge(
    const lsdtt_giuh_result&   giuh,
    const std::vector<double>& rainfall_m_per_s,
    double                     dt_rainfall_s,
    double                     runoff_coefficient)
{
    if (giuh.iuh.empty())
        throw std::runtime_error("giuh_discharge: IUH is empty");
    if (dt_rainfall_s <= 0.0)
        throw std::runtime_error("giuh_discharge: dt_rainfall_s must be > 0");
    if (giuh.dt_bin_s <= 0.0)
        throw std::runtime_error("giuh_discharge: giuh.dt_bin_s must be > 0");
    if (rainfall_m_per_s.empty())
        return {};

    const double A = giuh.catchment_area_m2;

    // The convolution indexes rainfall by dt_rainfall_s and the IUH by
    // giuh.dt_bin_s. Summing Q[t + k] across the two is only coherent when they
    // are the same width, so the IUH is re-binned onto the rainfall timestep
    // first. Before 2026-09-16 this resampling was described in a comment but
    // never performed, and the convolution additionally scaled by dt_bin_s
    // where it should have scaled by dt_rainfall_s.
    //
    // The two errors compensated in the integral, so total runoff volume came
    // out right and no mass-balance check could see the defect. What was wrong
    // was the hydrograph SHAPE, whenever the two timesteps differed: peak
    // discharge was out by a factor of dt_bin_s/dt_rainfall_s and time-to-peak
    // by its reciprocal. giuh_pulse_discharge() passes the caller's own dt
    // straight through, and lsdtt-sediment-routing passes its model timestep
    // while binning the IUH with the independent giuh_dt_bin_s parameter, so
    // the mismatch was reachable from the CLI.
    // tests/cpp/test_giuh_discharge.cpp pins both halves.
    const std::vector<double>& iuh_binned = giuh.iuh;
    std::vector<double> iuh_resampled;
    if (std::abs(dt_rainfall_s - giuh.dt_bin_s) > 1.0e-9 * giuh.dt_bin_s)
    {
        iuh_resampled = resample_iuh_density(
            iuh_binned, giuh.dt_bin_s, dt_rainfall_s);
    }
    const std::vector<double>& iuh =
        iuh_resampled.empty() ? iuh_binned : iuh_resampled;

    const int n_rain = static_cast<int>(rainfall_m_per_s.size());
    const int n_iuh  = static_cast<int>(iuh.size());
    const int n_out  = n_rain + n_iuh - 1;

    std::vector<double> Q(static_cast<std::size_t>(n_out), 0.0);

    for (int t = 0; t < n_rain; ++t)
    {
        const double r_eff =
            rainfall_m_per_s[static_cast<std::size_t>(t)] * runoff_coefficient;
        if (r_eff <= 0.0) continue;

        // Effective rainfall volume produced by this step, in m^3:
        //     A [m^2] * r_eff [m/s] * dt_rainfall_s [s]
        // The IUH is a density in [s^-1], so multiplying that volume by
        // iuh[k] gives m^3/s -- a discharge -- at lag k.
        const double step_volume_m3 = A * r_eff * dt_rainfall_s;

        for (int k = 0; k < n_iuh; ++k)
        {
            const int out_idx = t + k;
            if (out_idx >= n_out) break;
            Q[static_cast<std::size_t>(out_idx)] +=
                step_volume_m3 * iuh[static_cast<std::size_t>(k)];
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
