// lsdtt_step_change_cpu.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html
//
// Implementation of the analytical step-change rasters. The derivation and the
// numbering of equations (1)-(4) are in the header; read that first.
// -----------------------------------------------------------------------------

#include "lsdtt_step_change_cpu.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <vector>

#include "lsdtt_chi_extraction_cpu.hpp"

namespace lsdtt
{

namespace
{

constexpr double k_sqrt_2 = 1.4142135623730950488;

// Relative tolerance for "the two uplift rates are the same". Below this there
// is no step to solve and the shock/fan algebra is 0/0.
constexpr double k_no_step_rel_tol = 1.0e-12;

// Step length to the receiver. Even direction codes are cardinal, odd are
// diagonal -- the same convention the chi extraction uses.
double step_length_to_receiver(std::uint8_t dir, double cell_size)
{
    if (dir == lsdtt_flow_base::no_direction)
    {
        throw std::runtime_error(
            "lsdtt_step_change_cpu: encountered no_direction on a non-outlet "
            "node while walking the drainage stack.");
    }
    return (dir % 2u == 0u) ? cell_size : cell_size * k_sqrt_2;
}

bool tag_is_adjusted(StepChangeTag t)
{
    return t == StepChangeTag::fluvial_new  || t == StepChangeTag::hillslope_new
        || t == StepChangeTag::fluvial_fan  || t == StepChangeTag::hillslope_fan;
}

// Hillslope counterpart of a channel tag, for inheritance down the stack.
StepChangeTag hillslope_tag_for(StepChangeTag channel_tag)
{
    switch (channel_tag)
    {
        case StepChangeTag::fluvial_new:
        case StepChangeTag::hillslope_new:   return StepChangeTag::hillslope_new;
        case StepChangeTag::fluvial_fan:
        case StepChangeTag::hillslope_fan:   return StepChangeTag::hillslope_fan;
        default:                             return StepChangeTag::hillslope_old;
    }
}

} // namespace

// -----------------------------------------------------------------------------
double step_change_ksn(double uplift_m_yr, double K,
                       double n, double m_over_n, double A_0)
{
    return std::pow(uplift_m_yr / K, 1.0 / n) * std::pow(A_0, -m_over_n);
}

// -----------------------------------------------------------------------------
// Equation (1). A_0^m with m = m_over_n * n.
double step_change_celerity(double ksn, double K,
                            double n, double m_over_n, double A_0)
{
    return n * K * std::pow(A_0, m_over_n * n) * std::pow(ksn, n - 1.0);
}

// -----------------------------------------------------------------------------
// Invert ksn_old = relief / chi_max for K.
double step_change_K_from_relief(double uplift_old_m_yr, double target_relief_m,
                                 double chi_max, double n, double m_over_n,
                                 double A_0)
{
    if (!(chi_max > 0.0))
        throw std::invalid_argument(
            "lsdtt_step_change_cpu: chi_max must be positive to back-calculate K.");

    const double ksn_old = target_relief_m / chi_max;
    // ksn = (U/K)^(1/n) A_0^(-m/n)  =>  K = U / ( ksn A_0^(m/n) )^n
    return uplift_old_m_yr / std::pow(ksn_old * std::pow(A_0, m_over_n), n);
}

// -----------------------------------------------------------------------------
std::string StepChangeResult::summary() const
{
    std::ostringstream os;
    os << "Step-change rasters\n"
       << "-------------------\n";
    os << "  basin pixels        : " << n_basin_pixels
       << "  (channels " << n_channel_pixels << ")\n";
    os << "  outlet elevation    : " << outlet_elevation_m << " m\n";
    os << "  chi_max             : " << chi_max << "\n";
    os << std::scientific << std::setprecision(6);
    os << "  K (back-calculated) : " << K << "\n";
    os << std::defaultfloat << std::setprecision(6);
    os << "  ksn_old / ksn_new   : " << ksn_old << " / " << ksn_new << "\n";

    if (no_step)
    {
        os << "  NO STEP: the two uplift rates are equal, so the landscape is\n"
              "  a single steady state and time_since_change is zero everywhere.\n";
        return os.str();
    }

    if (is_shock)
    {
        os << "  wave form           : SHOCK (a migrating knickpoint)\n"
           << "  chi_knickpoint      : " << chi_knickpoint << "\n"
           << "  z_knickpoint        : " << z_knickpoint_m << " m\n";
    }
    else
    {
        os << "  wave form           : STRETCH (a rarefaction fan)\n"
           << "  fan in chi          : " << chi_fan_lower << " to "
           << chi_fan_upper << "\n"
           << "  leading edge z      : " << z_knickpoint_m << " m\n";
    }

    os << "  adjusted pixels     : " << n_adjusted_pixels << " of "
       << n_basin_pixels;
    if (n_basin_pixels > 0)
    {
        os << "  (" << (100.0 * static_cast<double>(n_adjusted_pixels)
                        / static_cast<double>(n_basin_pixels)) << " %)";
    }
    os << "\n";
    os << "  adjusted area       : " << adjusted_area_m2 << " m^2\n";

    if (wave_left_basin)
    {
        os << "\n  WARNING: the wave has passed the top of the basin -- every\n"
              "  pixel is adjusted and nothing records the old erosion rate.\n"
              "  Reduce step_age_yr, or accept that this run carries no\n"
              "  transient contrast.\n";
    }

    os << "\n  NOTE: K is back-calculated from target_relief_old_m, so the\n"
          "  landscape's whole configuration follows from a relief you chose.\n"
          "  Changing channel_threshold_pixels or the chi parameters changes K\n"
          "  and therefore changes every number above.\n";

    return os.str();
}

// -----------------------------------------------------------------------------
StepChangeResult make_step_change_rasters(const lsdtt_flow_base&  flow,
                                          std::int64_t            outlet_flat,
                                          const StepChangeConfig& config)
{
    // ---- validate ----------------------------------------------------------
    if (!(config.n > 0.0))
        throw std::invalid_argument(
            "lsdtt_step_change_cpu: n must be positive.");
    if (!(config.uplift_old_m_yr > 0.0) || !(config.uplift_new_m_yr > 0.0))
        throw std::invalid_argument(
            "lsdtt_step_change_cpu: both uplift rates must be positive (m/yr).");
    if (!(config.A_0 > 0.0))
        throw std::invalid_argument(
            "lsdtt_step_change_cpu: A_0 must be positive.");
    if (!(config.target_relief_old_m > 0.0))
        throw std::invalid_argument(
            "lsdtt_step_change_cpu: target_relief_old_m must be positive.");
    if (!(config.step_age_yr >= 0.0))
        throw std::invalid_argument(
            "lsdtt_step_change_cpu: step_age_yr must not be negative.");
    if (!flow.has_s_vector())
        throw std::invalid_argument(
            "lsdtt_step_change_cpu: the flow object must be built with "
            "build_s_vector = true.");
    if (!flow.is_valid_cell(outlet_flat))
        throw std::invalid_argument(
            "lsdtt_step_change_cpu: outlet_flat is not a valid cell.");

    const lsdtt_raster_meta& meta = flow.meta();
    const double cell_size = flow.cell_size();

    // ---- chi over the whole upslope domain ---------------------------------
    // The existing extraction is used rather than a second implementation of
    // the chi integral. The full-domain variant is required, not the
    // channel-only one: hillslope pixels need a chi value too.
    lsdtt_chi_options chi_opts;
    chi_opts.m_over_n = config.m_over_n;
    chi_opts.A_0      = config.A_0;

    const std::vector<std::int64_t> outlets{outlet_flat};
    lsdtt_raster_cpu chi_raster =
        get_chi_raster_from_outlets_zero_base(flow, outlets, chi_opts);

    const lsdtt_flat_index_slice basin = flow.upslope_nodes(outlet_flat);
    if (basin.empty())
        throw std::runtime_error(
            "lsdtt_step_change_cpu: empty upslope set for the outlet.");

    const std::vector<float>& chi_data = chi_raster.get_data();

    double chi_max = 0.0;
    for (std::int64_t flat : basin)
    {
        const double c = static_cast<double>(chi_data[static_cast<std::size_t>(flat)]);
        if (std::isfinite(c) && c > chi_max) chi_max = c;
    }
    if (!(chi_max > 0.0))
        throw std::runtime_error(
            "lsdtt_step_change_cpu: chi_max is not positive -- the basin has no "
            "extent above its outlet.");

    // ---- the analytic solution --------------------------------------------
    StepChangeResult r;
    r.outlet_flat = outlet_flat;
    r.chi_max     = chi_max;
    r.K = step_change_K_from_relief(config.uplift_old_m_yr,
                                    config.target_relief_old_m, chi_max,
                                    config.n, config.m_over_n, config.A_0);
    r.ksn_old = step_change_ksn(config.uplift_old_m_yr, r.K,
                                config.n, config.m_over_n, config.A_0);
    r.ksn_new = step_change_ksn(config.uplift_new_m_yr, r.K,
                                config.n, config.m_over_n, config.A_0);

    const double dU = config.uplift_new_m_yr - config.uplift_old_m_yr;
    r.no_step = std::fabs(dU)
              <= k_no_step_rel_tol * std::fabs(config.uplift_old_m_yr);

    // Which way the transition resolves is decided by the CELERITIES, not by n
    // alone: if the new steepness travels faster it overtakes the old and forms
    // a shock; if it travels slower the transition stretches. Doing it this way
    // means a DECREASE in uplift -- which swaps the two cases -- comes out right
    // without a special branch.
    const double c_old = step_change_celerity(r.ksn_old, r.K, config.n,
                                              config.m_over_n, config.A_0);
    const double c_new = step_change_celerity(r.ksn_new, r.K, config.n,
                                              config.m_over_n, config.A_0);
    r.is_shock = r.no_step || (config.n == 1.0) || (c_new >= c_old);

    const double t = config.step_age_yr;

    if (r.no_step)
    {
        r.chi_knickpoint = 0.0;
        r.chi_fan_lower  = 0.0;
        r.chi_fan_upper  = 0.0;
    }
    else if (r.is_shock)
    {
        // Equation (2). The denominator cannot vanish here: ksn_new == ksn_old
        // only if the uplift rates are equal, which is the no_step branch.
        r.chi_knickpoint = dU * t / (r.ksn_new - r.ksn_old);
        r.chi_fan_lower  = r.chi_knickpoint;
        r.chi_fan_upper  = r.chi_knickpoint;
    }
    else
    {
        // Rarefaction: the fan runs from the trailing edge (new steepness,
        // faster... no -- slower, which is why it stretches) to the leading
        // edge (old steepness). c_new < c_old here by construction.
        r.chi_fan_lower  = c_new * t;
        r.chi_fan_upper  = c_old * t;
        r.chi_knickpoint = r.chi_fan_upper;
    }

    r.wave_left_basin = !r.no_step && (r.chi_knickpoint >= chi_max);

    // The DEM is NOT an input: this routine generates a synthetic topography
    // consistent with the imposed uplift history, it does not perturb a real
    // surface. The outlet is therefore the datum at z = 0 and every elevation
    // returned is relative to it. A caller wanting the real base level adds the
    // measured outlet elevation to the whole raster.
    r.outlet_elevation_m = 0.0;

    // ---- fan profile constants (only used when !is_shock) ------------------
    // Inside the fan, S(chi) = (chi / (n K A_0^m t))^(1/(n-1)), equation (4).
    // Integrating gives the elevation; both exponents are precomputed because
    // they are used per pixel.
    const double fan_scale = config.n * r.K
                           * std::pow(config.A_0, config.m_over_n * config.n) * t;
    const double fan_exp   = (config.n == 1.0)
                           ? 0.0 : 1.0 / (config.n - 1.0);   // exponent on chi
    const double fan_int_exp = fan_exp + 1.0;                // after integration

    // z at the trailing edge, from the fully adjusted reach below it.
    const double z_fan_lower = r.ksn_new * r.chi_fan_lower;

    // Integral of S over the fan, from chi_fan_lower to chi.
    auto fan_elevation = [&](double chi_here) -> double
    {
        if (fan_int_exp == 0.0)
        {
            // Would require n -> infinity; unreachable for finite n > 0 with a
            // fan, but guard rather than divide by zero.
            return z_fan_lower;
        }
        const double k = std::pow(fan_scale, -fan_exp) / fan_int_exp;
        return z_fan_lower
             + k * (std::pow(chi_here, fan_int_exp)
                    - std::pow(r.chi_fan_lower, fan_int_exp));
    };

    if (r.no_step)
        r.z_knickpoint_m = 0.0;
    else if (r.is_shock)
        r.z_knickpoint_m = r.ksn_new * r.chi_knickpoint;
    else
        r.z_knickpoint_m = fan_elevation(r.chi_fan_upper);

    // ---- allocate outputs --------------------------------------------------
    // Built as raw arrays and promoted to rasters at the end: lsdtt_raster_cpu
    // exposes no mutable mask accessor, deliberately.
    const float ndv = static_cast<float>(meta.no_data_value);
    lsdtt_array2d_f  topo_a(meta.n_rows, meta.n_cols, ndv);
    lsdtt_array2d_u8 topo_k(meta.n_rows, meta.n_cols, 0);
    lsdtt_array2d_f  er_a(meta.n_rows, meta.n_cols, ndv);
    lsdtt_array2d_u8 er_k(meta.n_rows, meta.n_cols, 0);
    lsdtt_array2d_f  ts_a(meta.n_rows, meta.n_cols, ndv);
    lsdtt_array2d_u8 ts_k(meta.n_rows, meta.n_cols, 0);
    lsdtt_array2d_f  tg_a(meta.n_rows, meta.n_cols, ndv);
    lsdtt_array2d_u8 tg_k(meta.n_rows, meta.n_cols, 0);

    std::vector<float>&   topo_d = topo_a.get_vector();
    std::vector<uint8_t>& topo_m = topo_k.get_vector();
    std::vector<float>&   er_d   = er_a.get_vector();
    std::vector<uint8_t>& er_m   = er_k.get_vector();
    std::vector<float>&   ts_d   = ts_a.get_vector();
    std::vector<uint8_t>& ts_m   = ts_k.get_vector();
    std::vector<float>&   tg_d   = tg_a.get_vector();
    std::vector<uint8_t>& tg_m   = tg_k.get_vector();

    // Tag per basin pixel, kept alongside the raster so the hillslope walk can
    // read a receiver's class without decoding a float.
    std::vector<StepChangeTag> tag_of(static_cast<std::size_t>(meta.n_rows)
                                      * static_cast<std::size_t>(meta.n_cols),
                                      StepChangeTag::outside_basin);

    // ---- pass 1: the fluvial solution on channel pixels --------------------
    // The upslope slice is a preorder walk of the drainage stack: a node always
    // appears after its receiver. Pass 1 fills channels; pass 2 (the same
    // ordering) fills hillslopes from their receivers, so one ordered sweep
    // could do both -- they are kept separate only because the hillslope rule
    // needs to know that its receiver's value is final.
    for (std::int64_t flat : basin)
    {
        const std::size_t idx = static_cast<std::size_t>(flat);
        const double chi_here = static_cast<double>(chi_data[idx]);
        if (!std::isfinite(chi_here)) continue;

        const bool is_channel =
            static_cast<int>(flow.s_vector_size_unchecked(flat))
                > config.channel_threshold_pixels;

        double z, e, ts;
        StepChangeTag tg;

        if (r.no_step)
        {
            z  = r.ksn_old * chi_here;
            e  = config.uplift_old_m_yr;
            ts = 0.0;
            tg = StepChangeTag::fluvial_old;
        }
        else if (r.is_shock)
        {
            if (chi_here < r.chi_knickpoint)
            {
                z  = r.ksn_new * chi_here;
                e  = config.uplift_new_m_yr;
                // Equation (3): chi_kp is linear in t, so the fraction of the
                // way up to the knickpoint is the fraction of the elapsed time
                // that has passed since the wave went by.
                ts = t * (1.0 - chi_here / r.chi_knickpoint);
                tg = StepChangeTag::fluvial_new;
            }
            else
            {
                z  = r.z_knickpoint_m + (chi_here - r.chi_knickpoint) * r.ksn_old;
                e  = config.uplift_old_m_yr;
                ts = 0.0;
                tg = StepChangeTag::fluvial_old;
            }
        }
        else if (chi_here < r.chi_fan_lower)
        {
            z  = r.ksn_new * chi_here;
            e  = config.uplift_new_m_yr;
            ts = t * (1.0 - chi_here / r.chi_fan_lower);
            tg = StepChangeTag::fluvial_new;
        }
        else if (chi_here <= r.chi_fan_upper)
        {
            // Inside the fan: the steepness, and hence the erosion rate, vary
            // continuously. The change arrived when the LEADING edge passed.
            const double s_here = std::pow(chi_here / fan_scale, fan_exp);
            z  = fan_elevation(chi_here);
            e  = r.K * std::pow(config.A_0, config.m_over_n * config.n)
               * std::pow(s_here, config.n);
            ts = t * (1.0 - chi_here / r.chi_fan_upper);
            tg = StepChangeTag::fluvial_fan;
        }
        else
        {
            z  = r.z_knickpoint_m + (chi_here - r.chi_fan_upper) * r.ksn_old;
            e  = config.uplift_old_m_yr;
            ts = 0.0;
            tg = StepChangeTag::fluvial_old;
        }

        ++r.n_basin_pixels;

        if (!is_channel)
        {
            // Hillslope: elevation comes from pass 2. Everything else is
            // provisional and will be overwritten by inheritance.
            tag_of[idx] = hillslope_tag_for(tg);
        }
        else
        {
            ++r.n_channel_pixels;
            tag_of[idx] = tg;
            topo_d[idx] = static_cast<float>(z);
            topo_m[idx] = 1;
        }

        er_d[idx] = static_cast<float>(e);
        er_m[idx] = 1;
        ts_d[idx] = static_cast<float>(ts);
        ts_m[idx] = 1;
    }

    // ---- pass 2: impose hillslopes -----------------------------------------
    // Walk the stack again in the same receiver-before-donor order. A hillslope
    // pixel sits one step length above its receiver at the critical gradient,
    // and inherits the receiver's erosion rate, time-since-change and class --
    // a hillslope belongs to the channel it drains to.
    bool outlet_seen = false;
    for (std::int64_t flat : basin)
    {
        const std::size_t idx = static_cast<std::size_t>(flat);
        if (tag_of[idx] == StepChangeTag::outside_basin) continue;

        if (!outlet_seen)
        {
            // The outlet is the first entry. If it is itself below the channel
            // threshold the basin is smaller than the threshold; keep its
            // fluvial elevation as the datum rather than trying to inherit.
            outlet_seen = true;
            if (topo_m[idx] == 0)
            {
                topo_d[idx] = static_cast<float>(
                    r.ksn_old * static_cast<double>(chi_data[idx]));
                topo_m[idx] = 1;
            }
            continue;
        }

        if (topo_m[idx] != 0) continue;   // already a channel

        const std::int64_t rec = flow.receiver_of_unchecked(flat);
        if (rec < 0 || topo_m[static_cast<std::size_t>(rec)] == 0)
        {
            // Receiver outside the basin or not yet solved: cannot happen given
            // the ordering, but leave the pixel as nodata rather than invent a
            // value if it ever does.
            continue;
        }

        const std::size_t ridx = static_cast<std::size_t>(rec);
        const StepChangeTag rec_tag = tag_of[ridx];
        const StepChangeTag my_tag  = hillslope_tag_for(rec_tag);

        double s_c;
        switch (my_tag)
        {
            case StepChangeTag::hillslope_new: s_c = config.S_c_new; break;
            case StepChangeTag::hillslope_fan:
                s_c = 0.5 * (config.S_c_old + config.S_c_new);  break;
            default:                           s_c = config.S_c_old; break;
        }

        const double dx = step_length_to_receiver(
            flow.receiver_dir_of_unchecked(flat), cell_size);

        topo_d[idx] = static_cast<float>(
            static_cast<double>(topo_d[ridx]) + dx * s_c);
        topo_m[idx] = 1;

        tag_of[idx] = my_tag;
        er_d[idx]   = er_d[ridx];
        ts_d[idx]   = ts_d[ridx];
        er_m[idx]   = 1;
        ts_m[idx]   = 1;
    }

    // ---- tags raster and counts --------------------------------------------
    for (std::int64_t flat : basin)
    {
        const std::size_t idx = static_cast<std::size_t>(flat);
        const StepChangeTag tg = tag_of[idx];
        if (tg == StepChangeTag::outside_basin) continue;

        tg_d[idx] = static_cast<float>(static_cast<int>(tg));
        tg_m[idx] = 1;
        if (tag_is_adjusted(tg)) ++r.n_adjusted_pixels;
    }
    r.adjusted_area_m2 = static_cast<double>(r.n_adjusted_pixels)
                       * cell_size * cell_size;

    r.topography           = lsdtt_raster_cpu(std::move(topo_a),
                                              std::move(topo_k), meta);
    r.erosion_rate_m_yr    = lsdtt_raster_cpu(std::move(er_a),
                                              std::move(er_k), meta);
    r.time_since_change_yr = lsdtt_raster_cpu(std::move(ts_a),
                                              std::move(ts_k), meta);
    r.tags                 = lsdtt_raster_cpu(std::move(tg_a),
                                              std::move(tg_k), meta);
    r.chi                  = std::move(chi_raster);

    return r;
}

} // namespace lsdtt
