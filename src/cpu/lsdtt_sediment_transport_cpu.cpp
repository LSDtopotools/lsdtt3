// lsdtt_sediment_transport_cpu.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html

#include "lsdtt_sediment_transport_cpu.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <numeric>
#include <stdexcept>

#include "lsdtt_types.hpp"
#include "../utils/lsdtt_raster_geometry.hpp"

namespace lsdtt
{

// =============================================================================
// Internal helpers
// =============================================================================

namespace
{

constexpr double k_g            = 9.81;   // gravitational acceleration [m/s^2]
constexpr double k_min_slope    = 1.0e-6;
constexpr double k_min_q        = 1.0e-9; // minimum discharge [m^3/s]
constexpr double k_min_volume   = 1.0e-9; // minimum sediment volume [m^3]

// ---------------------------------------------------------------------------
// Manning flow depth for a wide rectangular channel.
// h = (Q * n / (W * sqrt(S)))^(3/5)
// Returns 0 if Q, W, or S are below sensible minima.
// ---------------------------------------------------------------------------
double manning_flow_depth(double Q, double W, double n, double S)
{
    if (Q < k_min_q || W <= 0.0 || n <= 0.0) return 0.0;
    double S_use = std::max(S, k_min_slope);
    return std::pow((Q * n) / (W * std::sqrt(S_use)), 0.6);
}

// ---------------------------------------------------------------------------
// Bed shear stress [Pa] for a wide rectangular channel.
// tau = rho_w * g * h * S
// ---------------------------------------------------------------------------
double bed_shear_stress(double h, double S, double rho_w)
{
    return rho_w * k_g * h * std::max(S, k_min_slope);
}

// ---------------------------------------------------------------------------
// Compute the D50 of an active layer from its volume distribution and the
// grain-class centre diameters (psi_m, in metres).
// Returns 0 if the layer is empty.
// ---------------------------------------------------------------------------
double active_layer_d50(
    const lsdtt_active_layer&   al,
    const std::vector<double>&  psi_m)
{
    const int nc = static_cast<int>(psi_m.size());
    double total = 0.0;
    for (int i = 0; i < nc; ++i) total += al.volume_m3[i];
    if (total < k_min_volume) return 0.0;

    // Cumulative fractional volume to locate the 50th percentile
    double cumul = 0.0;
    for (int i = 0; i < nc; ++i)
    {
        double f = al.volume_m3[i] / total;
        if (cumul + f >= 0.5)
        {
            // Linear interpolation within the class
            double f_prev = 0.5 - cumul;
            return psi_m[i] * (f_prev / f) +
                   (i > 0 ? psi_m[i-1] : psi_m[0]) * (1.0 - f_prev / f);
        }
        cumul += f;
    }
    return psi_m[nc - 1];
}

// ---------------------------------------------------------------------------
// Update active layer grain-size fractions after deposition or erosion.
// volume_change[i] > 0 -> add sediment to class i
// volume_change[i] < 0 -> remove sediment from class i (clamped to zero)
// ---------------------------------------------------------------------------
void update_active_layer(
    lsdtt_active_layer&        al,
    const std::vector<double>& volume_change)
{
    const int nc = static_cast<int>(al.volume_m3.size());
    for (int i = 0; i < nc; ++i)
    {
        al.volume_m3[i] = std::max(0.0, al.volume_m3[i] + volume_change[i]);
    }
}

// ---------------------------------------------------------------------------
// Walk pixels from junction j_up to its receiver junction, accumulating the
// planimetric Euclidean length of the link.
// ---------------------------------------------------------------------------
double compute_link_length(
    std::int32_t                      j_up,
    const lsdtt_junction_network_cpu& jn,
    const lsdtt_flowinfo&         flow)
{
    const auto& jid  = jn.junction_id_of_flat_array();
    const auto& meta = jn.meta();
    const std::int64_t start_flat = jn.flat_of_junction(j_up);

    std::int64_t current   = start_flat;
    double       link_len  = 0.0;

    for (int step = 0; step < 200000; ++step)
    {
        std::int64_t rec = flow.receiver_of_unchecked(current);

        // Reached an outlet pixel
        if (rec == lsdtt::k_outlet_receiver) break;

        // Accumulate step distance
        int r0 = flow.row_from_flat_unchecked(current);
        int c0 = flow.col_from_flat_unchecked(current);
        int r1 = flow.row_from_flat_unchecked(rec);
        int c1 = flow.col_from_flat_unchecked(rec);

        lsdtt_xy p0 = cell_center_from_row_col(meta, r0, c0);
        lsdtt_xy p1 = cell_center_from_row_col(meta, r1, c1);
        double dx = p1.x - p0.x;
        double dy = p1.y - p0.y;
        link_len += std::sqrt(dx * dx + dy * dy);

        current = rec;

        // Stop when we hit the receiver junction pixel (any junction that
        // is not the start of this link)
        if (jid.data()[current] >= 0 && current != start_flat) break;
    }

    // Minimum reach length of one diagonal cell to avoid zero-slope issues
    return std::max(link_len, meta.cell_size);
}

// ---------------------------------------------------------------------------
// Grain-class fractions from the active layer.
// Returns uniform fractions if the layer is empty.
// ---------------------------------------------------------------------------
std::vector<double> active_layer_fractions(
    const lsdtt_active_layer&  al,
    const std::vector<double>& default_fractions)
{
    const int nc = static_cast<int>(al.volume_m3.size());
    double total = 0.0;
    for (int i = 0; i < nc; ++i) total += al.volume_m3[i];
    if (total < k_min_volume) return default_fractions;

    std::vector<double> fracs(nc);
    for (int i = 0; i < nc; ++i) fracs[i] = al.volume_m3[i] / total;
    return fracs;
}

// =============================================================================
// Transport capacity formulas
// =============================================================================

// ---------------------------------------------------------------------------
// Wong & Parker (2006) with Parker (1982) hiding function.
//
// For each grain class i:
//   tau*_i  = tau / ((rho_s - rho_w) * g * D_i)          Shields number
//   tau*_ri = 0.0495 * (D_i / D50)^0.33             hiding-corrected critical Shields
//   if tau*_i > tau*_ri:
//     phi*_i = 3.97 * (tau*_i - tau*_ri)^1.5
//     q_bi = f_i * phi*_i * sqrt((s-1) * g * D_i^3) [m^2/s]
//
// Q_bi [m^3/s] = q_bi * W
// ---------------------------------------------------------------------------
std::vector<double> tr_cap_wong_parker(
    const lsdtt_sediment_reach& reach,
    double                      h,
    const std::vector<double>&  psi_m,
    const std::vector<double>&  fractions,
    double                      rho_s,
    double                      rho_w)
{
    const int nc = static_cast<int>(psi_m.size());
    std::vector<double> Q_per_class(nc, 0.0);

    double tau = bed_shear_stress(h, reach.slope, rho_w);
    if (tau <= 0.0) return Q_per_class;

    double delta      = (rho_s - rho_w) / rho_w;      // submerged specific gravity - 1
    double D50        = reach.D50_m;

    for (int i = 0; i < nc; ++i)
    {
        if (fractions[i] <= 0.0) continue;
        double Di     = psi_m[i];
        double tau_si = tau / ((rho_s - rho_w) * k_g * Di);  // Shields for class i
        double tau_ri = 0.0495 * std::pow(Di / D50, 0.33);   // hiding-corrected critical

        if (tau_si > tau_ri)
        {
            double diff   = tau_si - tau_ri;
            double phi    = 3.97 * std::pow(diff, 1.5);
            double q_bi   = phi * std::sqrt(delta * k_g * Di * Di * Di);  // m^2/s
            Q_per_class[i] = fractions[i] * q_bi * reach.width_m;
        }
    }
    return Q_per_class;
}

// ---------------------------------------------------------------------------
// Wilcock & Crowe (2003).
//
// tau_r50 = (rho_s - rho_w) * g * D50 * (0.021 + 0.015 * exp(-20 * F_s))
//   F_s = total sand fraction (D < 2 mm)
//
// For each class i:
//   b_i   = 0.67 / (1 + exp(1.5 - D_i / D50))    grain-size-dependent hiding
//   tau_ri  = tau_r50 * (D_i / D50)^b_i
//   u_i   = tau / tau_ri
//   if u_i < 1.35: W*_i = 0.002 * u_i^7.5
//   else:          W*_i = 14 * (1 - 0.894 / sqrt(u_i))^4.5
//   q_bi  = f_i * W*_i * u_star^3 / ((s-1) * g)  [m^2/s]
//
// Q_bi [m^3/s] = q_bi * W
// ---------------------------------------------------------------------------
std::vector<double> tr_cap_wilcock_crowe(
    const lsdtt_sediment_reach& reach,
    double                      h,
    const std::vector<double>&  psi_m,
    const std::vector<double>&  fractions,
    double                      rho_s,
    double                      rho_w)
{
    const int nc = static_cast<int>(psi_m.size());
    std::vector<double> Q_per_class(nc, 0.0);

    double tau = bed_shear_stress(h, reach.slope, rho_w);
    if (tau <= 0.0) return Q_per_class;

    double u_star = std::sqrt(tau / rho_w);      // friction velocity [m/s]
    double delta  = (rho_s - rho_w) / rho_w;     // submerged specific gravity - 1
    double D50    = reach.D50_m;

    // Sand fraction F_s (D < 2 mm = 0.002 m)
    double F_s = 0.0;
    for (int i = 0; i < nc; ++i)
    {
        if (psi_m[i] < 0.002) F_s += fractions[i];
    }

    // Reference shear stress for D50
    double tau_r50 = (rho_s - rho_w) * k_g * D50 *
                     (0.021 + 0.015 * std::exp(-20.0 * F_s));

    for (int i = 0; i < nc; ++i)
    {
        if (fractions[i] <= 0.0) continue;
        double Di  = psi_m[i];
        double bi  = 0.67 / (1.0 + std::exp(1.5 - Di / D50));
        double tau_ri = tau_r50 * std::pow(Di / D50, bi);

        double u_i = tau / tau_ri;
        double Wsi;
        if (u_i < 1.35)
            Wsi = 0.002 * std::pow(u_i, 7.5);
        else
        {
            double term = 1.0 - 0.894 / std::sqrt(u_i);
            Wsi = 14.0 * std::pow(term, 4.5);
        }

        double q_bi = fractions[i] * Wsi * u_star * u_star * u_star /
                      (delta * k_g);   // m^2/s
        Q_per_class[i] = q_bi * reach.width_m;
    }
    return Q_per_class;
}

// ---------------------------------------------------------------------------
// Engelund & Hansen (1967).
//
// Total-load formula for D50; distributes by GSD fraction.
//
//   tau*     = tau / ((rho_s - rho_w) * g * D50)
//   q_b    = 0.05 * tau*^(5/2) * sqrt((s-1) * g * D50^3)  [m^2/s, total]
//   Q_bi   = f_i * q_b * W
// ---------------------------------------------------------------------------
std::vector<double> tr_cap_engelund_hansen(
    const lsdtt_sediment_reach& reach,
    double                      h,
    const std::vector<double>&  psi_m,
    const std::vector<double>&  fractions,
    double                      rho_s,
    double                      rho_w)
{
    const int nc = static_cast<int>(psi_m.size());
    std::vector<double> Q_per_class(nc, 0.0);

    double tau = bed_shear_stress(h, reach.slope, rho_w);
    if (tau <= 0.0) return Q_per_class;

    double D50    = reach.D50_m;
    double delta  = (rho_s - rho_w) / rho_w;
    double tau_s  = tau / ((rho_s - rho_w) * k_g * D50);  // Shields number

    double q_b_total = 0.05 * std::pow(tau_s, 2.5) *
                       std::sqrt(delta * k_g * D50 * D50 * D50);  // m^2/s

    double Q_b_total = q_b_total * reach.width_m;  // m^3/s

    for (int i = 0; i < nc; ++i)
        Q_per_class[i] = fractions[i] * Q_b_total;

    return Q_per_class;
}

} // anonymous namespace

// =============================================================================
// Public implementations
// =============================================================================

// ---------------------------------------------------------------------------
// build_sediment_reaches
// ---------------------------------------------------------------------------
std::vector<lsdtt_sediment_reach> build_sediment_reaches(
    const lsdtt_junction_network_cpu&     jn,
    const lsdtt_flowinfo&             flow,
    const lsdtt_raster_cpu&               dem,
    const lsdtt_sediment_routing_options& opts)
{
    const int nj = jn.n_junctions();
    std::vector<lsdtt_sediment_reach> reaches(nj);

    const auto& dem_data   = dem.get_data();
    const auto& dem_mask   = dem.get_mask();
    const auto& meta       = dem.get_meta();
    const double cell_size = meta.cell_size;

    for (int j = 0; j < nj; ++j)
    {
        lsdtt_sediment_reach& r = reaches[j];
        r.junction_id          = j;
        r.receiver_junction_id = jn.receiver_of_junction(j);

        // ------------------------------------------------------------------
        // Geometry from DEM
        // ------------------------------------------------------------------
        std::int64_t flat_up = jn.flat_of_junction(j);

        r.elevation_upstream_m =
            (dem_mask.data()[flat_up] != 0)
                ? static_cast<double>(dem_data.data()[flat_up])
                : 0.0;

        // Downstream elevation
        if (jn.is_outlet_junction(j))
        {
            std::int64_t flat_out = jn.basin_outlet_flat_of_junction(j);
            r.elevation_downstream_m =
                (dem_mask.data()[flat_out] != 0)
                    ? static_cast<double>(dem_data.data()[flat_out])
                    : r.elevation_upstream_m;
        }
        else
        {
            std::int64_t flat_dn =
                jn.flat_of_junction(r.receiver_junction_id);
            r.elevation_downstream_m =
                (dem_mask.data()[flat_dn] != 0)
                    ? static_cast<double>(dem_data.data()[flat_dn])
                    : r.elevation_upstream_m;
        }

        // Link length via pixel walk
        r.length_m = compute_link_length(j, jn, flow);

        // Slope: clamp to minimum to avoid divide-by-zero
        double elev_drop = r.elevation_upstream_m - r.elevation_downstream_m;
        r.slope = std::max(elev_drop / r.length_m, opts.min_slope);

        // Drainage area from s_vector pixel count
        r.drainage_area_m2 =
            static_cast<double>(flow.s_vector_size_unchecked(flat_up)) *
            cell_size * cell_size;

        // ------------------------------------------------------------------
        // Width from hydraulic geometry: W = k * A^b
        // ------------------------------------------------------------------
        r.width_m = opts.width_coefficient *
                    std::pow(r.drainage_area_m2, opts.width_exponent);
        r.width_m = std::max(r.width_m, 0.5);  // 50 cm minimum

        // ------------------------------------------------------------------
        // Grain sizes and Manning's n from global opts (user may override)
        // ------------------------------------------------------------------
        r.mannings_n   = opts.mannings_n;
        r.D50_m        = opts.D50_m;
        r.D84_m        = opts.D84_m;
        r.D90_m        = opts.D90_m;
        r.gsd_fractions = opts.gsd_fractions;
    }

    return reaches;
}

// ---------------------------------------------------------------------------
// initialise_active_layers
// ---------------------------------------------------------------------------
std::vector<lsdtt_active_layer> initialise_active_layers(
    const std::vector<lsdtt_sediment_reach>& reaches,
    const lsdtt_sediment_routing_options&    opts,
    double                                   initial_depth_m)
{
    const int nj = static_cast<int>(reaches.size());
    const int nc = static_cast<int>(opts.psi_mm.size());
    std::vector<lsdtt_active_layer> layers(nj);

    for (int j = 0; j < nj; ++j)
    {
        const lsdtt_sediment_reach& r = reaches[j];
        lsdtt_active_layer& al        = layers[j];

        // Active layer depth
        double depth;
        if (initial_depth_m > 0.0)
            depth = initial_depth_m;
        else if (opts.active_layer_depth_m > 0.0)
            depth = opts.active_layer_depth_m;
        else
            depth = 2.0 * r.D90_m;   // Parker (1991) two-grain-diameter rule

        al.depth_m = depth;

        // Initial volume per class:
        //   V_i = f_i * depth * W * length * (1 - porosity)
        // We store bulk volume (solids + pores) so omit (1-phi) here; the
        // transport capacity is also computed as bulk volume flux.
        double reach_bed_area = r.width_m * r.length_m;
        al.volume_m3.resize(nc);
        for (int i = 0; i < nc; ++i)
            al.volume_m3[i] = r.gsd_fractions[i] * depth * reach_bed_area;
    }

    return layers;
}

// ---------------------------------------------------------------------------
// run_sediment_routing
// ---------------------------------------------------------------------------
lsdtt_sediment_routing_result run_sediment_routing(
    const lsdtt_junction_network_cpu&          jn,
    const std::vector<lsdtt_sediment_reach>&   reaches,
    std::vector<lsdtt_active_layer>&           active_layers,
    const std::vector<std::vector<double>>&    discharge,
    const lsdtt_sediment_routing_options&      opts)
{
    const int nj  = jn.n_junctions();
    const int nt  = opts.n_timesteps;
    const int nc  = static_cast<int>(opts.psi_mm.size());
    const double dt = opts.timestep_s;

    // Convert psi_mm to metres for transport calculations
    std::vector<double> psi_m(nc);
    for (int i = 0; i < nc; ++i) psi_m[i] = opts.psi_mm[i] * 1.0e-3;

    const double rho_s = opts.sediment_density_kg_m3;
    const double rho_w = opts.water_density_kg_m3;

    // Build result container
    lsdtt_sediment_routing_result result;
    result.n_junctions  = nj;
    result.n_timesteps  = nt;
    result.cumulative_deposition_m3.assign(nj, 0.0);
    result.cumulative_mobilisation_m3.assign(nj, 0.0);
    result.final_d50_al_m.assign(nj, 0.0);
    result.timestep_results.resize(nt, std::vector<lsdtt_sediment_reach_result>(nj));

    // Topological order (headwaters first)
    const auto& topo = jn.topological_order_vector();

    // Per-junction outflux per grain class for the current timestep [m^3]
    // outflux_per_class[j][i] = volume leaving junction j's reach this timestep
    std::vector<std::vector<double>> outflux(nj, std::vector<double>(nc, 0.0));

    // ------------------------------------------------------------------
    // Main time loop
    // ------------------------------------------------------------------
    for (int t = 0; t < nt; ++t)
    {
        // Zero outflux for this timestep
        for (int j = 0; j < nj; ++j)
            std::fill(outflux[j].begin(), outflux[j].end(), 0.0);

        // ------------------------------------------------------------------
        // Spatial loop: process each junction in topological order
        // ------------------------------------------------------------------
        for (int ti = 0; ti < nj; ++ti)
        {
            std::int32_t j = topo[ti];
            const lsdtt_sediment_reach& reach = reaches[j];
            lsdtt_active_layer&          al    = active_layers[j];

            // --------------------------------------------------------------
            // Discharge at this junction and timestep
            // --------------------------------------------------------------
            double Q;
            if (!discharge.empty() &&
                t < static_cast<int>(discharge.size()) &&
                j < static_cast<int>(discharge[t].size()))
            {
                Q = discharge[t][j];
            }
            else
            {
                // Hydraulic geometry fallback
                Q = opts.discharge_coefficient *
                    std::pow(reach.drainage_area_m2, opts.discharge_exponent);
            }
            Q = std::max(Q, k_min_q);

            // --------------------------------------------------------------
            // Incoming sediment from upstream junctions
            // --------------------------------------------------------------
            std::vector<double> incoming(nc, 0.0);
            jn.for_each_donor_junction(j,
                [&](std::int32_t donor)
                {
                    for (int i = 0; i < nc; ++i)
                        incoming[i] += outflux[donor][i];
                });

            double incoming_total = 0.0;
            for (int i = 0; i < nc; ++i) incoming_total += incoming[i];

            // --------------------------------------------------------------
            // Transport capacity
            // --------------------------------------------------------------
            double h = manning_flow_depth(Q, reach.width_m,
                                          reach.mannings_n, reach.slope);

            std::vector<double> fracs = active_layer_fractions(
                al, reach.gsd_fractions);

            std::vector<double> Qs_per_class;  // m^3/s per class
            switch (opts.formula)
            {
            case lsdtt_tr_cap_formula::wong_parker:
                Qs_per_class = tr_cap_wong_parker(
                    reach, h, psi_m, fracs, rho_s, rho_w);
                break;
            case lsdtt_tr_cap_formula::wilcock_crowe:
                Qs_per_class = tr_cap_wilcock_crowe(
                    reach, h, psi_m, fracs, rho_s, rho_w);
                break;
            case lsdtt_tr_cap_formula::engelund_hansen:
                Qs_per_class = tr_cap_engelund_hansen(
                    reach, h, psi_m, fracs, rho_s, rho_w);
                break;
            }

            // Capacity in m^3 over this timestep
            double Qs_total = 0.0;
            for (int i = 0; i < nc; ++i) Qs_total += Qs_per_class[i];
            double capacity_m3 = Qs_total * dt;

            // --------------------------------------------------------------
            // Mobilisation / deposition
            // --------------------------------------------------------------
            std::vector<double> delta_al(nc, 0.0);  // change to active layer
            double deposited_m3  = 0.0;
            double mobilised_m3  = 0.0;

            if (incoming_total <= capacity_m3)
            {
                // Can transport everything incoming; try to mobilise more
                double deficit = capacity_m3 - incoming_total;
                double al_total = 0.0;
                for (int i = 0; i < nc; ++i) al_total += al.volume_m3[i];

                double mobilise = std::min(deficit, al_total);
                mobilised_m3   = mobilise;

                if (al_total > k_min_volume && mobilise > 0.0)
                {
                    for (int i = 0; i < nc; ++i)
                    {
                        double removed = mobilise * (al.volume_m3[i] / al_total);
                        delta_al[i]   -= removed;
                        outflux[j][i]  = incoming[i] + removed;
                    }
                }
                else
                {
                    for (int i = 0; i < nc; ++i)
                        outflux[j][i] = incoming[i];
                }
            }
            else
            {
                // Incoming exceeds capacity; deposit the excess
                double excess    = incoming_total - capacity_m3;
                deposited_m3     = excess;

                for (int i = 0; i < nc; ++i)
                {
                    double deposit_i;
                    if (incoming_total > k_min_volume)
                        deposit_i = excess * (incoming[i] / incoming_total);
                    else
                        deposit_i = 0.0;

                    delta_al[i]   += deposit_i;
                    outflux[j][i]  = incoming[i] - deposit_i;
                    outflux[j][i]  = std::max(outflux[j][i], 0.0);
                }
            }

            // Apply active layer changes
            update_active_layer(al, delta_al);

            // --------------------------------------------------------------
            // Accumulate outlet flux
            // --------------------------------------------------------------
            if (jn.is_outlet_junction(j))
            {
                for (int i = 0; i < nc; ++i)
                    result.total_flux_out_m3 += outflux[j][i];
            }

            // --------------------------------------------------------------
            // Record this timestep's reach result
            // --------------------------------------------------------------
            lsdtt_sediment_reach_result& rr = result.timestep_results[t][j];
            rr.deposited_m3               = deposited_m3;
            rr.mobilised_m3               = mobilised_m3;
            rr.d50_al_m                   = active_layer_d50(al, psi_m);
            rr.transport_capacity_m3_per_s = Qs_total;

            result.cumulative_deposition_m3[j]   += deposited_m3;
            result.cumulative_mobilisation_m3[j] += mobilised_m3;
        }
    }

    // Final D50 from the last active layer state
    for (int j = 0; j < nj; ++j)
    {
        result.final_d50_al_m[j] =
            active_layer_d50(active_layers[j], psi_m);
    }

    return result;
}

// ---------------------------------------------------------------------------
// make_deposition_raster
// ---------------------------------------------------------------------------
lsdtt_raster_cpu make_deposition_raster(
    const lsdtt_junction_network_cpu&     jn,
    const lsdtt_flowinfo&             flow,
    const lsdtt_raster_cpu&               dem,
    const lsdtt_sediment_routing_result&  result)
{
    const auto& meta  = dem.get_meta();
    const int   n_cells = meta.n_rows * meta.n_cols;

    std::vector<float>   data(n_cells, static_cast<float>(meta.no_data_value));
    std::vector<uint8_t> mask(n_cells, 0);

    const int nj = jn.n_junctions();

    // For each junction j, assign its cumulative deposition to all pixels
    // in the link from j to its receiver.
    for (int j = 0; j < nj; ++j)
    {
        float dep_val = static_cast<float>(result.cumulative_deposition_m3[j]);

        // Walk pixels from junction j to its receiver, assigning the value
        std::int64_t current = jn.flat_of_junction(j);
        const std::int64_t start_flat = current;
        const auto& jid = jn.junction_id_of_flat_array();

        for (int step = 0; step < 200000; ++step)
        {
            data[current]  = dep_val;
            mask[current]  = 1;

            std::int64_t rec = flow.receiver_of_unchecked(current);
            if (rec == lsdtt::k_outlet_receiver) break;
            current = rec;
            if (jid.data()[current] >= 0 && current != start_flat) break;
        }
    }

    lsdtt_array2d_f   arr(meta.n_rows, meta.n_cols, std::move(data));
    lsdtt_array2d_u8  msk(meta.n_rows, meta.n_cols, std::move(mask));
    return lsdtt_raster_cpu(std::move(arr), std::move(msk), meta);
}

// ---------------------------------------------------------------------------
// make_d50_raster
// ---------------------------------------------------------------------------
lsdtt_raster_cpu make_d50_raster(
    const lsdtt_junction_network_cpu&     jn,
    const lsdtt_flowinfo&             flow,
    const lsdtt_raster_cpu&               dem,
    const lsdtt_sediment_routing_result&  result)
{
    const auto& meta  = dem.get_meta();
    const int   n_cells = meta.n_rows * meta.n_cols;

    std::vector<float>   data(n_cells, static_cast<float>(meta.no_data_value));
    std::vector<uint8_t> mask(n_cells, 0);

    const int nj = jn.n_junctions();

    for (int j = 0; j < nj; ++j)
    {
        float d50_val = static_cast<float>(result.final_d50_al_m[j]);
        std::int64_t current = jn.flat_of_junction(j);
        const std::int64_t start_flat = current;
        const auto& jid = jn.junction_id_of_flat_array();

        for (int step = 0; step < 200000; ++step)
        {
            data[current]  = d50_val;
            mask[current]  = 1;

            std::int64_t rec = flow.receiver_of_unchecked(current);
            if (rec == lsdtt::k_outlet_receiver) break;
            current = rec;
            if (jid.data()[current] >= 0 && current != start_flat) break;
        }
    }

    lsdtt_array2d_f   arr(meta.n_rows, meta.n_cols, std::move(data));
    lsdtt_array2d_u8  msk(meta.n_rows, meta.n_cols, std::move(mask));
    return lsdtt_raster_cpu(std::move(arr), std::move(msk), meta);
}

} // namespace lsdtt
