// lsdtt_landscape_equilibration_cpu.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html

#include "lsdtt_landscape_equilibration_cpu.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cctype>
#include <cstring>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>

#include "lsdtt_flowinfo_cpu.hpp"
#include "lsdtt_critical_slope_cpu.hpp"
#include "lsdtt_raster_flow_preprocessing_cpu.hpp"

namespace lsdtt
{

namespace
{

double flow_step_length_to_receiver(const lsdtt_flowinfo& flow, std::int64_t flat)
{
    constexpr double k_sqrt_2 = 1.4142135623730950488;
    const std::uint8_t dir = flow.receiver_dir_of_unchecked(flat);
    return ((dir & 1u) == 0u) ? flow.cell_size() : (flow.cell_size() * k_sqrt_2);
}

double tolerance_multiplier_for_cycle(double schedule, int cycle, int max_cycles)
{
    if (!std::isfinite(schedule) || schedule <= 1.0 || max_cycles <= 1)
        return 1.0;

    const double frac_remaining =
        static_cast<double>(max_cycles - cycle) / static_cast<double>(max_cycles - 1);
    return std::pow(schedule, std::max(0.0, frac_remaining));
}

// Flow-path distance to the nearest channel, summarised. The mean is the
// obvious summary; the p90 is the one the gradient prediction needs, because
// the solver is pushed hardest by the LONGEST hillslopes, not the average one.
// Both come from the same distance sweep, so the percentile is nearly free.
struct hillslope_length_stats
{
    double mean_m = 0.0;
    double p90_m = 0.0;
};

hillslope_length_stats estimate_hillslope_length_stats(
    const lsdtt_flowinfo& flow,
    const std::vector<std::uint8_t>& is_channel)
{
    const std::size_t n = flow.full_grid_size();
    std::vector<double> dist_to_channel(n, -1.0);
    for (std::size_t s = 0; s < n; ++s)
    {
        if (s < is_channel.size() && is_channel[s] != 0)
            dist_to_channel[s] = 0.0;
    }

    double sum = 0.0;
    std::int64_t count = 0;
    std::vector<double> sampled;

    for (std::int64_t flat = 0; flat < static_cast<std::int64_t>(n); ++flat)
    {
        if (!flow.is_valid_cell(flat)) continue;
        const std::size_t s0 = static_cast<std::size_t>(flat);
        if (s0 >= is_channel.size() || is_channel[s0] != 0) continue;

        if (dist_to_channel[s0] >= 0.0)
        {
            sum += dist_to_channel[s0];
            sampled.push_back(dist_to_channel[s0]);
            ++count;
            continue;
        }

        std::vector<std::int64_t> path;
        std::vector<double> step_lens;
        std::int64_t cur = flat;
        double known_tail = 0.0;
        bool anchored = false;

        for (std::size_t guard = 0; guard <= n; ++guard)
        {
            if (!flow.is_valid_cell(cur))
            {
                anchored = true;
                known_tail = 0.0;
                break;
            }

            const std::size_t cs = static_cast<std::size_t>(cur);
            if (dist_to_channel[cs] >= 0.0)
            {
                anchored = true;
                known_tail = dist_to_channel[cs];
                break;
            }

            path.push_back(cur);

            const std::int64_t rec = flow.receiver_of_unchecked(cur);
            if (rec < 0 || !flow.is_valid_cell(rec))
            {
                anchored = true;
                known_tail = 0.0;
                break;
            }

            step_lens.push_back(flow_step_length_to_receiver(flow, cur));
            cur = rec;
        }

        if (!anchored)
            continue;

        double tail = known_tail;
        for (std::size_t p = path.size(); p-- > 0;)
        {
            if (p < step_lens.size())
                tail += step_lens[p];
            dist_to_channel[static_cast<std::size_t>(path[p])] = tail;
        }

        if (dist_to_channel[s0] >= 0.0)
        {
            sum += dist_to_channel[s0];
            sampled.push_back(dist_to_channel[s0]);
            ++count;
        }
    }

    hillslope_length_stats stats;
    if (count > 0)
    {
        stats.mean_m = sum / static_cast<double>(count);

        // Nearest-rank p90, matching the convention used by the divide-step
        // percentiles in lsdtt_divide_diagnostics_cpu.
        std::sort(sampled.begin(), sampled.end());
        const std::size_t idx = std::min(
            sampled.size() - 1,
            static_cast<std::size_t>(0.9 * static_cast<double>(sampled.size())));
        stats.p90_m = sampled[idx];
    }
    return stats;
}

// Steady-state relief of a 1-D hillslope of length L under the nonlinear
// (Roering) transport law, obtained by integrating the closed-form gradient.
// With flux q(x) = U*x measured from the divide and q = D*S/(1-(S/Sc)^2),
//
//     S(x) = Sc^2 * ( sqrt(1 + (2*a*x/Sc)^2) - 1 ) / (2*a*x),   a = U/D
//
// and the relief is the integral of S over [0, L]. As a -> 0 this tends to
// S = a*x and the relief to U*L^2/(2D), i.e. it reduces to the linear answer,
// which is the check the unit test makes. Integrated numerically (composite
// trapezium, 4096 intervals) rather than by the analytic antiderivative: this
// runs once per equilibration, and the numerical form is far easier to read
// and to verify against the linear limit.
double predict_roering_hillslope_relief(double U, double D, double Sc, double L)
{
    if (!std::isfinite(U) || !std::isfinite(D) || !std::isfinite(Sc)
        || !std::isfinite(L) || U <= 0.0 || D <= 0.0 || Sc <= 0.0 || L <= 0.0)
    {
        return 0.0;
    }

    const double a = U / D;
    constexpr int    n_steps = 4096;
    const double     h = L / static_cast<double>(n_steps);

    const auto gradient_at = [&](double x)
    {
        if (x <= 0.0)
            return 0.0;
        const double t = 2.0 * a * x / Sc;
        // Sc^2 * (sqrt(1 + t^2) - 1) / (2*a*x). Written this way to stay
        // accurate when t is small, where sqrt(1+t^2)-1 cancels badly.
        const double num = t * t / (std::sqrt(1.0 + t * t) + 1.0);
        return Sc * Sc * num / (2.0 * a * x);
    };

    double relief = 0.5 * (gradient_at(0.0) + gradient_at(L));
    for (int i = 1; i < n_steps; ++i)
        relief += gradient_at(static_cast<double>(i) * h);

    return relief * h;
}

double predict_max_hillslope_gradient(double U, double D, double Sc, double length_m)
{
    if (!std::isfinite(U) || !std::isfinite(D) || !std::isfinite(Sc) || !std::isfinite(length_m)
        || U <= 0.0 || D <= 0.0 || Sc <= 0.0 || length_m <= 0.0)
    {
        return 0.0;
    }

    const double rhs = U * length_m / D;
    double lo = 0.0;
    double hi = Sc * (1.0 - 1e-12);

    for (int it = 0; it < 80; ++it)
    {
        const double mid = 0.5 * (lo + hi);
        const double ratio = mid / Sc;
        const double lhs = mid / (1.0 - ratio * ratio);
        if (lhs < rhs) lo = mid;
        else hi = mid;
    }
    return 0.5 * (lo + hi);
}

lsdtt_snap_cycle_record make_cycle_record(
    int                           cycle,
    const lsdtt_raster_cpu&       previous_surface,
    const lsdtt_raster_cpu&       current_surface,
    const lsdtt_flowinfo*         previous_flow,
    const lsdtt_flowinfo&         current_flow,
    const lsdtt_fluvial_steady_state_result& snap_result,
    const std::vector<std::uint8_t>* previous_is_channel,
    const std::vector<std::uint8_t>& current_is_channel,
    double                        elapsed_seconds)
{
    lsdtt_snap_cycle_record rec;
    rec.cycle = cycle;
    rec.seconds = elapsed_seconds;

    const lsdtt_raster_meta& meta = current_surface.get_meta();
    const std::vector<float>& z_prev = previous_surface.get_data();
    const std::vector<float>& z_cur = current_surface.get_data();
    const std::vector<uint8_t>& mask_prev = previous_surface.get_mask();
    const std::vector<uint8_t>& mask_cur = current_surface.get_mask();

    const std::int64_t n_cells = static_cast<std::int64_t>(std::min(
        z_prev.size(), z_cur.size()));
    double sum_abs = 0.0;
    double sum_sq = 0.0;
    std::int64_t n_eval = 0;

    double max_abs = -1.0;
    std::int64_t max_flat = -1;

    for (std::int64_t flat = 0; flat < n_cells; ++flat)
    {
        const std::size_t s = static_cast<std::size_t>(flat);
        if (!mask_cur[s] || !mask_prev[s])
            continue;

        const double dz = static_cast<double>(z_cur[s])
                        - static_cast<double>(z_prev[s]);
        const double abs_dz = std::abs(dz);

        sum_abs += abs_dz;
        sum_sq += dz * dz;
        ++n_eval;

        if (abs_dz > max_abs)
        {
            max_abs = abs_dz;
            max_flat = flat;
        }
    }

    if (n_eval > 0)
    {
        rec.max_abs_dz_m = std::max(0.0, max_abs);
        rec.mean_abs_dz_m = sum_abs / static_cast<double>(n_eval);
        rec.rms_dz_m = std::sqrt(sum_sq / static_cast<double>(n_eval));
    }

    rec.max_abs_dz_flat = max_flat;
    if (max_flat >= 0)
    {
        rec.max_abs_dz_row = static_cast<int>(
            max_flat / static_cast<std::int64_t>(meta.n_cols));
        rec.max_abs_dz_col = static_cast<int>(
            max_flat % static_cast<std::int64_t>(meta.n_cols));

        rec.max_abs_dz_x = meta.x_min +
            (static_cast<double>(rec.max_abs_dz_col) + 0.5) * meta.cell_size;
        rec.max_abs_dz_y = meta.y_min +
            (static_cast<double>(meta.n_rows - rec.max_abs_dz_row) - 0.5)
            * meta.cell_size;
    }

    const double relief = std::max(
        static_cast<double>(snap_result.computed_relief), 1e-12);
    rec.max_dz_over_relief = rec.max_abs_dz_m / relief;

    rec.K_used = static_cast<double>(snap_result.K_used);
    rec.computed_relief = static_cast<double>(snap_result.computed_relief);

    if (previous_flow)
    {
        const std::vector<std::int64_t>& prev_recv =
            previous_flow->receiver_array().get_vector();
        const std::vector<std::int64_t>& cur_recv =
            current_flow.receiver_array().get_vector();
        const std::vector<uint8_t>& prev_flow_mask =
            previous_flow->mask_array().get_vector();
        const std::vector<uint8_t>& cur_flow_mask =
            current_flow.mask_array().get_vector();

        std::int64_t n_compare = 0;
        std::int64_t n_changed = 0;
        std::int64_t n_hillslope_compare = 0;
        std::int64_t n_hillslope_changed = 0;
        const std::int64_t n_flow_cells = static_cast<std::int64_t>(std::min(
            prev_recv.size(), cur_recv.size()));

        for (std::int64_t flat = 0; flat < n_flow_cells; ++flat)
        {
            const std::size_t s = static_cast<std::size_t>(flat);
            if (!prev_flow_mask[s] || !cur_flow_mask[s])
                continue;

            ++n_compare;
            const bool changed = prev_recv[s] != cur_recv[s];
            if (changed)
                ++n_changed;

            // Restrict to pixels that were hillslope on both cycles. A pixel
            // that entered or left the channel set is a channel-set movement,
            // not a hillslope routing instability.
            if (previous_is_channel != nullptr
                && s < previous_is_channel->size()
                && s < current_is_channel.size()
                && (*previous_is_channel)[s] == 0
                && current_is_channel[s] == 0)
            {
                ++n_hillslope_compare;
                if (changed)
                    ++n_hillslope_changed;
            }
        }

        rec.n_receiver_changed = n_changed;
        if (n_compare > 0)
        {
            rec.frac_receiver_changed =
                static_cast<double>(n_changed) / static_cast<double>(n_compare);
        }

        rec.n_receiver_changed_hillslope = n_hillslope_changed;
        rec.n_hillslope_compared = n_hillslope_compare;
        if (n_hillslope_compare > 0)
        {
            rec.frac_receiver_changed_hillslope =
                static_cast<double>(n_hillslope_changed)
                / static_cast<double>(n_hillslope_compare);
        }
    }
    else
    {
        // First cycle has no prior flow graph for a receiver-change comparison.
        rec.n_receiver_changed = 0;
        rec.frac_receiver_changed = 0.0;
        rec.n_receiver_changed_hillslope = 0;
        rec.frac_receiver_changed_hillslope = 0.0;
        rec.n_hillslope_compared = 0;
    }

    return rec;
}

// 64-bit FNV-1a over the raw bytes of the receiver vector. O(N) and negligible
// against a snap; used only as a cheap index into the seen-graph table. Every
// suspected repeat is confirmed with a memcmp against the stored array, so a
// hash collision cannot produce a false limit-cycle report.
std::uint64_t fnv1a_receivers(const std::vector<std::int64_t>& receivers)
{
    constexpr std::uint64_t fnv_offset_basis = 1469598103934665603ULL;
    constexpr std::uint64_t fnv_prime        = 1099511628211ULL;

    std::uint64_t hash = fnv_offset_basis;
    const unsigned char* bytes =
        reinterpret_cast<const unsigned char*>(receivers.data());
    const std::size_t n_bytes = receivers.size() * sizeof(std::int64_t);

    for (std::size_t b = 0; b < n_bytes; ++b)
    {
        hash ^= static_cast<std::uint64_t>(bytes[b]);
        hash *= fnv_prime;
    }

    return hash;
}

// One previously visited flow graph, kept so that a suspected repeat can be
// confirmed byte-for-byte and so that an orbit member can be re-delivered.
struct lsdtt_seen_graph
{
    int                       cycle = 0;
    std::vector<std::int64_t> receivers;
};

// Count the cells that take part in an orbit: those whose receiver differs from
// the orbit's first member on at least one visited state.
std::int64_t count_orbit_pixels(
    const std::vector<const std::vector<std::int64_t>*>& orbit_receivers)
{
    if (orbit_receivers.size() < 2)
        return 0;

    const std::vector<std::int64_t>& first = *orbit_receivers.front();
    std::vector<std::uint8_t> participates(first.size(), 0);

    for (std::size_t m = 1; m < orbit_receivers.size(); ++m)
    {
        const std::vector<std::int64_t>& other = *orbit_receivers[m];
        const std::size_t n = std::min(first.size(), other.size());
        for (std::size_t s = 0; s < n; ++s)
        {
            if (first[s] != other[s])
                participates[s] = 1;
        }
    }

    std::int64_t n_pixels = 0;
    for (const std::uint8_t flag : participates)
        n_pixels += flag;

    return n_pixels;
}

// The router and the hillslope solver do not share a boundary-condition
// vocabulary (see lsdtt_hillslope_steady_state_cpu.hpp:40-43). They agree on
// 'p' (periodic) and on 'n' (no-flux in the router, zero-flux Neumann in the
// solver, which is the solver's 'o' branch). They disagree on:
//   'o' - the router lets water leave the grid; the solver reflects it back.
//   'b' - baselevel in the router; the solver reflects.
// Fixing that means editing the solver's numeric core, which is protected, so
// the cycle warns instead. Returns an empty string when the two are consistent.
std::string make_hillslope_bc_warning(const std::string& boundary_conditions)
{
    static const char* edge_names[4] = {"north", "east", "south", "west"};

    std::vector<std::string> offending;
    const std::size_t n_edges =
        std::min<std::size_t>(4, boundary_conditions.size());

    for (std::size_t e = 0; e < n_edges; ++e)
    {
        const char bc = static_cast<char>(
            std::tolower(static_cast<unsigned char>(boundary_conditions[e])));
        if (bc == 'o' || bc == 'b')
        {
            offending.push_back(
                std::string(edge_names[e]) + " ('" + bc + "')");
        }
    }

    if (offending.empty())
        return std::string();

    std::ostringstream oss;
    oss << " WARNING: hillslope boundary-condition mismatch on the ";
    for (std::size_t i = 0; i < offending.size(); ++i)
    {
        if (i > 0)
            oss << (i + 1 == offending.size() ? " and " : ", ");
        oss << offending[i];
    }
    oss << " edge" << (offending.size() > 1 ? "s" : "") << ":"
        << " the flow router drains through "
        << (offending.size() > 1 ? "these edges" : "this edge")
        << " while the hillslope steady-state solver reflects at "
        << (offending.size() > 1 ? "them" : "it")
        << " (zero-flux Neumann). The two steps are not solving the same"
        << " problem there, so divides within roughly one hillslope length of"
        << " the affected edge(s) are not trustworthy. Use 'n' or 'p' edges for"
        << " hillslope work.";

    return oss.str();
}

std::string make_non_converged_reason(
    int cycle_count,
    const lsdtt_snap_cycle_record& last,
    const lsdtt_snap_convergence_options& conv)
{
    std::ostringstream oss;
    oss << "Reached max_cycles=" << cycle_count << " without convergence";

    oss << " (max|dz|=" << last.max_abs_dz_m << " m";
    if (conv.tolerance_dz_m > 0.0)
        oss << ", target<=" << conv.tolerance_dz_m << " m";
    else if (conv.tolerance_dz_over_relief > 0.0)
        oss << ", dz/relief=" << last.max_dz_over_relief
            << ", target<=" << conv.tolerance_dz_over_relief;

    oss << "; receiver_changed=" << last.frac_receiver_changed;
    if (conv.tolerance_receiver_fraction > 0.0)
        oss << ", target<=" << conv.tolerance_receiver_fraction;
    oss << ")";

    return oss.str();
}

} // namespace

std::vector<std::string> lsdtt_snap_cycle_record::csv_headers()
{
    return {
        "cycle",
        "seconds",
        "max_abs_dz_m",
        "mean_abs_dz_m",
        "rms_dz_m",
        "max_dz_over_relief",
        "n_receiver_changed",
        "frac_receiver_changed",
        "n_receiver_changed_hillslope",
        "frac_receiver_changed_hillslope",
        "n_hillslope_compared",
        "n_cross_divide",
        "frac_cross_divide",
        "divide_n_divide_pairs",
        "divide_step_median_m",
        "divide_step_p90_m",
        "divide_step_max_m",
        "K_used",
        "computed_relief",
        "divide_gradient_median",
        "divide_gradient_p90",
        "divide_gradient_target_used",
        "divide_gradient_target_margin",
        "divide_gradient_target_met",
        "receiver_hash",
        "limit_cycle_detected",
        "limit_cycle_period",
        "limit_cycle_first",
        "limit_cycle_n_pixels",
        "n_channel_pixels",
        "hillslope_iterations",
        "hillslope_residual",
        "hillslope_converged"
    };
}

std::vector<std::string> lsdtt_snap_cycle_record::csv_row() const
{
    return {
        std::to_string(cycle),
        std::to_string(seconds),
        std::to_string(max_abs_dz_m),
        std::to_string(mean_abs_dz_m),
        std::to_string(rms_dz_m),
        std::to_string(max_dz_over_relief),
        std::to_string(n_receiver_changed),
        std::to_string(frac_receiver_changed),
        std::to_string(n_receiver_changed_hillslope),
        std::to_string(frac_receiver_changed_hillslope),
        std::to_string(n_hillslope_compared),
        std::to_string(n_cross_divide),
        std::to_string(frac_cross_divide),
        std::to_string(divide.n_divide_pairs),
        std::to_string(divide.step_median_m),
        std::to_string(divide.step_p90_m),
        std::to_string(divide.step_max_m),
        std::to_string(K_used),
        std::to_string(computed_relief),
        std::to_string(divide.gradient_median),
        std::to_string(divide.gradient_p90),
        std::to_string(divide_gradient_target_used),
        std::to_string(divide_gradient_target_margin),
        std::to_string(divide_gradient_target_met),
        std::to_string(receiver_hash),
        std::to_string(limit_cycle_detected),
        std::to_string(limit_cycle_period),
        std::to_string(limit_cycle_first),
        std::to_string(limit_cycle_n_pixels),
        std::to_string(n_channel_pixels),
        std::to_string(hillslope_iterations),
        std::to_string(hillslope_residual),
        std::to_string(hillslope_converged)
    };
}

lsdtt_equilibration_result equilibrate_landscape(
    const lsdtt_raster_cpu& initial_dem,
    const lsdtt_equilibration_options& opts)
{
    if (initial_dem.empty())
    {
        throw std::invalid_argument(
            "equilibrate_landscape: initial_dem is empty");
    }

    lsdtt_equilibration_result result;

    lsdtt_fluvial_steady_state_options fluvial_opts = opts.fluvial;
    fluvial_opts.n_threads = opts.n_threads;

    lsdtt_flowinfo_options fi_opts;
    fi_opts.boundary_conditions = opts.boundary_conditions;
    fi_opts.build_s_vector = false;
    fi_opts.build_topological_order = true;
    fi_opts.n_threads = opts.n_threads;

    const int max_cycles = std::max(1, opts.max_cycles);
    const int patience = std::max(1, opts.convergence.patience_cycles);
    int n_good_cycles = 0;

    bool hillslope_requested = opts.run_hillslope;
    std::string preflight_note;

    lsdtt_raster_cpu cycle_dem = initial_dem;
    lsdtt_raster_cpu previous_delivered = initial_dem;

    lsdtt_flowinfo previous_flow;
    bool has_previous_flow = false;

    // Limit-cycle bookkeeping. surface_history[0] is the input DEM and
    // surface_history[c] is the surface delivered by cycle c, so the input to
    // any cycle b is surface_history[b - 1] and an orbit member can be replayed
    // exactly. seen_graphs indexes previously visited flow graphs by an FNV-1a
    // hash of the receiver array, keeping the array itself so that a suspected
    // repeat is confirmed by memcmp rather than trusted on the hash.
    // channel_history[c] and hillslope_ran_history[c] are the channel mask and
    // the solve flag of cycle c, index 0 being an unused placeholder that keeps
    // them aligned with surface_history. The best-of-orbit replay needs both:
    // the channel mask carries hysteresis from the preceding cycle and the
    // preceding delivered surface is the warm start, so a replay that dropped
    // them re-ran the cycle under DIFFERENT inputs and handed back a surface
    // that was not the orbit member whose divide step was reported.
    std::vector<lsdtt_raster_cpu> surface_history;
    surface_history.push_back(initial_dem);
    std::vector<std::vector<std::uint8_t>> channel_history(1);
    std::vector<std::uint8_t> hillslope_ran_history(1, 0);
    std::map<std::uint64_t, lsdtt_seen_graph> seen_graphs;
    int last_executed_cycle = 0;

    // Channel mask from a surface and its drainage areas. Shared by the cycle
    // loop and the best-of-orbit replay so the two cannot drift. A pixel that
    // was already channel stays channel until its area falls below
    // channel_hysteresis_factor * threshold, so the Dirichlet set handed to the
    // hillslope solver does not flicker from cycle to cycle.
    const auto make_is_channel = [&](const lsdtt_raster_cpu& surface,
                                     const std::vector<double>& area_d,
                                     const std::vector<std::uint8_t>* previous)
    {
        std::vector<std::uint8_t> is_channel(surface.get_mask().size(), 0);
        if (opts.channel_threshold_area <= 0.0)
            return is_channel;

        const double on_threshold = opts.channel_threshold_area;
        const double off_threshold = std::min(
            on_threshold,
            std::max(0.0, opts.channel_hysteresis_factor) * on_threshold);

        const double dx = surface.get_meta().cell_size;
        for (std::size_t s = 0; s < area_d.size(); ++s)
        {
            if (!surface.get_mask()[s])
                continue;

            double area_for_threshold = area_d[s];
            if (fluvial_opts.divide_area_by_dx)
                area_for_threshold /= dx;

            const bool was_channel = previous != nullptr
                && s < previous->size() && (*previous)[s] != 0;

            if (area_for_threshold >= on_threshold)
                is_channel[s] = 1;
            else if (was_channel && area_for_threshold >= off_threshold)
                is_channel[s] = 1;
        }
        return is_channel;
    };

    // Hillslope steady-state solve on a snapped surface with the channel set
    // pinned as Dirichlet nodes. Returns the outcome and writes the solved
    // surface into out_surface; a no-op (ran = false) leaves out_surface alone.
    struct hillslope_outcome
    {
        bool  ran = false;
        int   iterations = 0;
        float residual = 0.0f;
        bool  converged = false;
    };

    const auto run_hillslope_solve =
        [&](const lsdtt_flowinfo&             flow,
            int                               cycle,
            const lsdtt_raster_cpu&           fluvial_surface,
            const std::vector<std::uint8_t>&  is_channel,
            const lsdtt_raster_cpu*           warm_start,
            lsdtt_raster_cpu&                 out_surface) -> hillslope_outcome
    {
        hillslope_outcome outcome;

        const std::int64_t n_channel = std::count(
            is_channel.begin(), is_channel.end(), static_cast<std::uint8_t>(1));
        if (n_channel == 0)
        {
            // With no Dirichlet nodes the elliptic problem has no unique
            // solution; skip rather than hand the solver an ill-posed system.
            return outcome;
        }

        lsdtt_hillslope_steady_state_options hs_opts = opts.hillslope;
        hs_opts.boundary_conditions = opts.boundary_conditions;
        hs_opts.n_threads = opts.n_threads;
        hs_opts.tolerance = static_cast<float>(
            static_cast<double>(hs_opts.tolerance)
            * tolerance_multiplier_for_cycle(
                opts.hillslope_tolerance_schedule,
                cycle,
                max_cycles));

        std::string model = opts.hillslope_model;
        for (char& c : model)
            c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));

        if (model == "roering")
            hs_opts.nonlinear = true;
        else if (model == "linear")
            hs_opts.nonlinear = false;
        else if (model == "critical_slope" || model == "critical-slope")
        {
            lsdtt_critical_slope_options cs_opts;
            cs_opts.mode = critical_slope_mode::fixed_gradient;
            cs_opts.S_c = static_cast<double>(opts.hillslope.Sc);
            cs_opts.is_channel_override = is_channel;

            // The cycle measures the divide step itself, on the delivered
            // surface, through compute_divide_step_stats and with its own
            // use_d8_diagonals setting. impose_critical_slope_hillslopes'
            // predicted_divide_step_m would be computed over a different pair
            // population (use_d8_diagonals defaults false there, true in
            // lsdtt_divide_diagnostic_options) and is discarded here, so the
            // path-set build it costs on every cycle is pure waste.
            cs_opts.compute_predicted_divide_step = false;

            const lsdtt_critical_slope_result cs_result =
                impose_critical_slope_hillslopes(flow, fluvial_surface, cs_opts);
            out_surface = cs_result.elevation;

            outcome.ran = true;
            outcome.iterations = 1;
            outcome.residual = 0.0f;
            outcome.converged = true;
            return outcome;
        }
        else if (model != "cidre")
        {
            throw std::invalid_argument(
                "equilibrate_landscape: unknown hillslope_model '"
                + opts.hillslope_model
                + "'; expected linear, roering, cidre, or critical_slope");
        }

        // Warm start: the previous cycle's hillslope field is the best initial
        // guess available. Channel pixels keep the current snap elevations,
        // which are the Dirichlet values.
        lsdtt_raster_cpu solve_input = fluvial_surface;
        if (opts.warm_start_hillslope
            && warm_start != nullptr
            && warm_start->get_data().size() == solve_input.get_data().size())
        {
            std::vector<float>& z = solve_input.get_data();
            const std::vector<float>& z_warm = warm_start->get_data();
            const std::vector<std::uint8_t>& m = solve_input.get_mask();
            for (std::size_t s = 0; s < z.size(); ++s)
            {
                if (m[s] != 0 && is_channel[s] == 0)
                    z[s] = z_warm[s];
            }
        }

        const lsdtt_hillslope_steady_state_result hs_result =
            (model == "cidre")
                ? compute_cidre_hillslope_steady_state(
                      solve_input, is_channel, hs_opts)
                : compute_hillslope_steady_state(
                      solve_input, is_channel, hs_opts);

        out_surface = hs_result.elevation;

        outcome.ran = true;
        outcome.iterations = hs_result.n_iterations;
        outcome.residual = hs_result.final_residual;
        outcome.converged = hs_result.converged;
        return outcome;
    };

    // The router and the solver disagree about 'o' and 'b' edges. Warn once,
    // loudly, and record it; do not throw, because 'o' is the lsdtt-lem default.
    if (opts.run_hillslope)
    {
        const std::string bc_warning =
            make_hillslope_bc_warning(opts.boundary_conditions);
        if (!bc_warning.empty())
        {
            result.hillslope_bc_mismatch = true;
            std::cerr << "Warning:" << bc_warning << "\n";
        }
    }

    std::vector<std::uint8_t> previous_is_channel;
    bool has_previous_is_channel = false;
    lsdtt_raster_cpu previous_hillslope_surface;
    bool has_previous_hillslope_surface = false;
    bool preflight_done = false;

    for (int cycle = 1; cycle <= max_cycles; ++cycle)
    {
        const auto t0 = std::chrono::steady_clock::now();

        cycle_dem = fill_pits(cycle_dem, opts.fill_min_slope, opts.n_threads);

        // TODO(dinf): add D-infinity route/accumulate/snap path.
        // Milestone scaffold intentionally runs D8 only.
        (void)opts.use_dinf;
        lsdtt_flowinfo flow(
            cycle_dem.get_data_array(),
            cycle_dem.get_mask_array(),
            cycle_dem.get_meta(),
            fi_opts);

        std::vector<double> area_d = flow.accumulate();
        lsdtt_fluvial_steady_state_result snap =
            compute_d8_fluvial_steady_state(flow, fluvial_opts);

        result.final_snap = snap;
        result.fluvial_only = snap.elevation;
        result.drainage_area = area_d;

        lsdtt_raster_cpu delivered = snap.elevation;

        // The channel set is an INPUT to the hillslope solve, not just an
        // output: it is the Dirichlet set that pins the elliptic problem.
        std::vector<std::uint8_t> is_channel = make_is_channel(
            snap.elevation,
            area_d,
            has_previous_is_channel ? &previous_is_channel : nullptr);

        if (!preflight_done && opts.run_hillslope)
        {
            const hillslope_length_stats lh_stats =
                estimate_hillslope_length_stats(flow, is_channel);
            result.mean_hillslope_length_m = lh_stats.mean_m;
            result.p90_hillslope_length_m  = lh_stats.p90_m;

            const double U = static_cast<double>(opts.hillslope.U);
            const double D = static_cast<double>(opts.hillslope.D);
            const double Sc = static_cast<double>(opts.hillslope.Sc);
            const double Lh = result.mean_hillslope_length_m;
            const double Lh90 = result.p90_hillslope_length_m;

            std::string preflight_model = opts.hillslope_model;
            for (char& c : preflight_model)
                c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));

            if (U > 0.0 && D > 0.0 && Lh > 0.0)
            {
                result.predicted_hillslope_relief_linear_m =
                    U * Lh * Lh / (2.0 * D);
            }

            // Model-aware relief. The linear form badly under-predicts what
            // critical-slope path integration imposes (Sc*L, an order of
            // magnitude more on a typical fixture), which is where using it for
            // every model made the relief_ratio warning quietest in exactly the
            // case that needs it loudest. For Roering it errs the other way --
            // the nonlinear law is the more efficient transporter, so its
            // relief is strictly SMALLER than linear and saturates at Sc*L.
            // Either way the right thing to report is the relief of the model
            // that was actually asked for.
            if (preflight_model == "roering")
            {
                result.predicted_hillslope_relief_m =
                    predict_roering_hillslope_relief(U, D, Sc, Lh);
            }
            else if (preflight_model == "critical_slope"
                     || preflight_model == "critical-slope")
            {
                // Path integration imposes exactly Sc per unit flow distance.
                result.predicted_hillslope_relief_m =
                    (Sc > 0.0 && Lh > 0.0) ? (Sc * Lh) : 0.0;
            }
            else
            {
                result.predicted_hillslope_relief_m =
                    result.predicted_hillslope_relief_linear_m;
            }

            // Evaluated at the p90 length: the gradient that pushes the
            // nonlinear solver towards Sc lives on the longest hillslopes, and
            // the mean-length value systematically under-fires the warning.
            // Both are reported so the size of that bias is visible.
            result.predicted_max_hillslope_gradient =
                predict_max_hillslope_gradient(U, D, Sc, Lh90 > 0.0 ? Lh90 : Lh);
            result.predicted_max_hillslope_gradient_at_mean_length =
                predict_max_hillslope_gradient(U, D, Sc, Lh);

            std::ostringstream pf;
            const double resolved_threshold_specific = fluvial_opts.divide_area_by_dx
                ? opts.channel_threshold_area
                : (opts.channel_threshold_area / std::max(1e-12, flow.cell_size()));
            const double resolved_threshold_map_area = fluvial_opts.divide_area_by_dx
                ? (opts.channel_threshold_area * flow.cell_size())
                : opts.channel_threshold_area;

            pf << " Pre-flight: mean_hillslope_length_m="
               << result.mean_hillslope_length_m
               << ", p90_hillslope_length_m="
               << result.p90_hillslope_length_m
               << ", predicted_hillslope_relief_m="
               << result.predicted_hillslope_relief_m
               << " (model '" << preflight_model
               << "'; linear-diffusion reference "
               << result.predicted_hillslope_relief_linear_m << " m)"
               << ", predicted_max_hillslope_gradient="
               << result.predicted_max_hillslope_gradient
               << " (at the p90 hillslope length; at the mean length it would"
                  " be "
               << result.predicted_max_hillslope_gradient_at_mean_length << ")"
               << ", channel_threshold_area(specific)="
               << resolved_threshold_specific
               << ", channel_threshold_area(map)="
               << resolved_threshold_map_area << ".";

            const double relief_ratio =
                result.predicted_hillslope_relief_m
                / std::max(1e-12, static_cast<double>(snap.computed_relief));
            if (relief_ratio > 0.5)
            {
                pf << " WARNING: predicted_hillslope_relief/computed_relief="
                   << relief_ratio << " (>0.5).";
            }
            if (Sc > 0.0 && result.predicted_max_hillslope_gradient > 0.9 * Sc)
            {
                pf << " WARNING: predicted_max_hillslope_gradient="
                   << result.predicted_max_hillslope_gradient
                   << " exceeds 0.9*Sc; consider lower sor_omega (e.g. 1.3).";
            }
            if (!opts.hillslope_every_cycle && max_cycles > 1)
            {
                pf << " WARNING: hillslope_every_cycle=false with max_cycles>1"
                   << " is a speed mode, not full route∘snap∘diffuse equilibration.";
            }

            preflight_note = pf.str();
            std::cerr << "Info:" << preflight_note << "\n";
            preflight_done = true;
        }

        hillslope_outcome hs;
        if (opts.run_hillslope && opts.hillslope_every_cycle)
        {
            hs = run_hillslope_solve(
                flow,
                cycle,
                snap.elevation,
                is_channel,
                has_previous_hillslope_surface ? &previous_hillslope_surface
                                               : nullptr,
                delivered);

            if (hs.ran)
            {
                previous_hillslope_surface = delivered;
                has_previous_hillslope_surface = true;
                result.hillslope_ran = true;

                if (!hs.converged)
                {
                    // Worth shouting about: an unconverged surface is not a
                    // steady state, so every divide metric taken from it is
                    // meaningless. The nonlinear solver in particular can
                    // diverge under the auto-omega estimate; try sor_omega=1.3.
                    std::cerr << "Warning: hillslope solve did not converge on"
                              << " cycle " << cycle << " (" << hs.iterations
                              << " iterations, residual " << hs.residual
                              << " vs tolerance " << opts.hillslope.tolerance
                              << "). Divide metrics for this cycle are not"
                                 " trustworthy; try a smaller sor_omega"
                                 " (e.g. 1.3) or more iterations.\n";
                }
            }
        }

        const auto t1 = std::chrono::steady_clock::now();
        const double elapsed_seconds =
            std::chrono::duration_cast<std::chrono::duration<double>>(t1 - t0)
                .count();

        lsdtt_snap_cycle_record rec = make_cycle_record(
            cycle,
            previous_delivered,
            delivered,
            has_previous_flow ? &previous_flow : nullptr,
            flow,
            snap,
            has_previous_is_channel ? &previous_is_channel : nullptr,
            is_channel,
            elapsed_seconds);

        rec.n_channel_pixels = std::count(
            is_channel.begin(), is_channel.end(), static_cast<std::uint8_t>(1));
        if (hs.ran)
        {
            rec.hillslope_iterations = hs.iterations;
            rec.hillslope_residual = static_cast<double>(hs.residual);
            rec.hillslope_converged = hs.converged ? 1 : 0;
        }

        rec.divide = compute_divide_step_stats(
            flow,
            delivered,
            opts.convergence.divide_opts);

        if (opts.convergence.divide_gradient_target > 0.0)
        {
            rec.divide_gradient_target_used =
                opts.convergence.divide_gradient_target;
            rec.divide_gradient_target_margin =
                opts.convergence.divide_gradient_target
                - rec.divide.gradient_median;
            rec.divide_gradient_target_met =
                rec.divide.gradient_median
                <= opts.convergence.divide_gradient_target ? 1 : 0;
        }

        if (has_previous_flow)
        {
            // Cross-divide fraction: SAME cell, TWO graphs. anc_K_prev[i] is
            // where cell i drained to before the re-route, anc_K_cur[i] is where
            // it drains to after it. An unchanged flow graph therefore gives
            // elementwise equality and exactly 0.0.
            //
            // The earlier formulation compared anc_prev[s] with
            // anc_prev[receiver(s)]; that is definitionally empty because a cell
            // and its own receiver drain to the same place by construction, so
            // it measured the one-step chain offset and returned ~1 under no
            // change at all.
            const int ancestor_steps =
                opts.convergence.divide_opts.divide_merge_search_steps;
            const std::vector<std::int64_t> anc_prev =
                compute_downstream_ancestor_labels(previous_flow, ancestor_steps);
            const std::vector<std::int64_t> anc_cur =
                compute_downstream_ancestor_labels(flow, ancestor_steps);

            const std::vector<std::uint8_t>& prev_flow_mask =
                previous_flow.mask_array().get_vector();
            const std::vector<std::uint8_t>& cur_mask =
                flow.mask_array().get_vector();

            std::int64_t n_compare = 0;
            std::int64_t n_cross = 0;
            const std::int64_t n_flow_cells = static_cast<std::int64_t>(
                std::min(anc_prev.size(), anc_cur.size()));

            for (std::int64_t flat = 0; flat < n_flow_cells; ++flat)
            {
                const std::size_t s = static_cast<std::size_t>(flat);
                if (cur_mask[s] == 0 || prev_flow_mask[s] == 0)
                    continue;

                ++n_compare;
                if (anc_prev[s] != anc_cur[s])
                    ++n_cross;
            }

            rec.n_cross_divide = n_cross;
            if (n_compare > 0)
            {
                rec.frac_cross_divide =
                    static_cast<double>(n_cross) / static_cast<double>(n_compare);
            }
        }
        else
        {
            rec.n_cross_divide = 0;
            rec.frac_cross_divide = 0.0;
        }

        // Fingerprint this cycle's flow graph before it is moved into
        // previous_flow, so a repeat of an earlier graph can be recognised.
        const std::vector<std::int64_t>& cycle_recv =
            flow.receiver_array().get_vector();
        rec.receiver_hash = fnv1a_receivers(cycle_recv);

        result.elevation = delivered;
        result.delivered_cycle = cycle;
        result.is_channel = is_channel;
        surface_history.push_back(delivered);
        channel_history.push_back(is_channel);
        hillslope_ran_history.push_back(hs.ran ? 1u : 0u);
        last_executed_cycle = cycle;

        bool stop_now = false;

        if (opts.convergence.check_convergence)
        {
            bool dz_ok = true;
            if (opts.convergence.tolerance_dz_m > 0.0)
            {
                dz_ok = rec.max_abs_dz_m <= opts.convergence.tolerance_dz_m;
            }
            else if (opts.convergence.tolerance_dz_over_relief > 0.0)
            {
                dz_ok = rec.max_dz_over_relief <=
                        opts.convergence.tolerance_dz_over_relief;
            }

            bool receiver_ok = true;
            if (opts.convergence.tolerance_receiver_fraction > 0.0)
            {
                receiver_ok = rec.frac_receiver_changed <=
                    opts.convergence.tolerance_receiver_fraction;
            }

            // Cycle 1 has no prior flow graph, so frac_receiver_changed is
            // forced to 0.0 and every clause can be satisfied vacuously -- with
            // both dz tolerances disabled a patience-1 run would otherwise
            // declare convergence having compared nothing. A cycle must have a
            // predecessor to count.
            if (cycle >= 2 && dz_ok && receiver_ok)
            {
                ++n_good_cycles;
            }
            else
            {
                n_good_cycles = 0;
            }

            if (n_good_cycles >= patience)
            {
                result.converged = true;
                std::ostringstream oss;
                oss << "Converged at cycle " << cycle
                    << " after " << patience << " consecutive cycle(s)"
                    << " (max|dz|=" << rec.max_abs_dz_m
                    << " m, dz/relief=" << rec.max_dz_over_relief
                    << ", receiver_changed=" << rec.frac_receiver_changed
                    << ")";
                result.stop_reason = oss.str();
                stop_now = true;
            }
        }

        // Limit-cycle detection. A repeat of an earlier flow graph means the
        // iteration has closed an orbit: every subsequent cycle is an exact
        // replay, so there is nothing further to learn.
        if (!stop_now && !result.limit_cycle_detected)
        {
            const auto it = seen_graphs.find(rec.receiver_hash);
            const bool confirmed_repeat =
                it != seen_graphs.end()
                && it->second.receivers.size() == cycle_recv.size()
                && std::memcmp(it->second.receivers.data(),
                               cycle_recv.data(),
                               cycle_recv.size() * sizeof(std::int64_t)) == 0;

            if (confirmed_repeat)
            {
                const int first_cycle = it->second.cycle;

                // Orbit members are the states visited at cycles
                // [first_cycle, cycle - 1]; the state at `cycle` repeats the
                // one at first_cycle.
                std::vector<const std::vector<std::int64_t>*> orbit_receivers;
                orbit_receivers.push_back(&it->second.receivers);
                for (const auto& entry : seen_graphs)
                {
                    if (entry.second.cycle > first_cycle
                        && entry.second.cycle < cycle)
                    {
                        orbit_receivers.push_back(&entry.second.receivers);
                    }
                }

                result.limit_cycle_detected = true;
                result.limit_cycle_first = first_cycle;
                result.limit_cycle_period = cycle - first_cycle;
                result.limit_cycle_n_pixels = count_orbit_pixels(orbit_receivers);

                rec.limit_cycle_detected = 1;
                rec.limit_cycle_first = result.limit_cycle_first;
                rec.limit_cycle_period = result.limit_cycle_period;
                rec.limit_cycle_n_pixels = result.limit_cycle_n_pixels;

                // Stopping is a convergence decision; a fixed-cycle run still
                // reports the orbit but must run the cycles it was asked for.
                if (opts.convergence.check_convergence)
                    stop_now = true;
            }
            else if (it == seen_graphs.end())
            {
                lsdtt_seen_graph entry;
                entry.cycle = cycle;
                entry.receivers = cycle_recv;
                seen_graphs.emplace(rec.receiver_hash, std::move(entry));
            }
            // A hash hit that fails memcmp is a collision: leave the stored
            // array in place and carry on rather than mis-report an orbit.
        }

        result.cycles.push_back(rec);

        if (stop_now)
            break;

        previous_delivered = delivered;
        previous_flow = std::move(flow);
        has_previous_flow = true;
        previous_is_channel = std::move(is_channel);
        has_previous_is_channel = true;
        cycle_dem = delivered;
    }

    if (result.limit_cycle_detected && opts.convergence.check_convergence)
    {
        std::ostringstream oss;

        if (result.limit_cycle_period == 1)
        {
            // The flow graph reproduced itself exactly: a genuine fixed point.
            result.converged = true;
            oss << "Converged at cycle " << last_executed_cycle
                << ": flow graph is a fixed point (limit_cycle_period=1,"
                << " first seen at cycle " << result.limit_cycle_first << ")";
        }
        else
        {
            // An orbit: deliver the smoothest member visited, not whatever
            // phase of the orbit the loop happened to stop on.
            int    best_cycle = result.limit_cycle_first;
            double best_p90 = -1.0;
            for (const lsdtt_snap_cycle_record& r : result.cycles)
            {
                if (r.cycle < result.limit_cycle_first
                    || r.cycle > last_executed_cycle - 1)
                    continue;

                if (best_p90 < 0.0 || r.divide.step_p90_m < best_p90)
                {
                    best_p90 = r.divide.step_p90_m;
                    best_cycle = r.cycle;
                }
            }

            if (best_cycle != last_executed_cycle
                && best_cycle >= 1
                && static_cast<std::size_t>(best_cycle) < surface_history.size())
            {
                // Replay the single cycle that produced the chosen orbit member
                // so the delivered elevation, snap, areas and channel mask are
                // all the same state rather than a mix of two phases.
                lsdtt_raster_cpu replay_dem = fill_pits(
                    surface_history[static_cast<std::size_t>(best_cycle) - 1],
                    opts.fill_min_slope,
                    opts.n_threads);

                lsdtt_flowinfo replay_flow(
                    replay_dem.get_data_array(),
                    replay_dem.get_mask_array(),
                    replay_dem.get_meta(),
                    fi_opts);

                std::vector<double> replay_area = replay_flow.accumulate();
                lsdtt_fluvial_steady_state_result replay_snap =
                    compute_d8_fluvial_steady_state(replay_flow, fluvial_opts);

                // Reproduce cycle best_cycle exactly: same hysteresis input to
                // the channel mask and same warm start as the loop used. With
                // either dropped the replay is a DIFFERENT cycle and the
                // elevation handed back does not match the divide step_p90
                // reported for the chosen orbit member.
                const std::size_t prev_idx =
                    static_cast<std::size_t>(best_cycle) - 1;
                const std::vector<std::uint8_t>* replay_prev_channel =
                    (best_cycle >= 2 && prev_idx < channel_history.size()
                     && !channel_history[prev_idx].empty())
                        ? &channel_history[prev_idx]
                        : nullptr;

                result.final_snap = replay_snap;
                result.fluvial_only = replay_snap.elevation;
                result.drainage_area = std::move(replay_area);
                result.elevation = replay_snap.elevation;
                result.is_channel = make_is_channel(
                    result.elevation, result.drainage_area, replay_prev_channel);

                if (opts.run_hillslope && opts.hillslope_every_cycle)
                {
                    const lsdtt_raster_cpu* replay_warm_start =
                        (best_cycle >= 2
                         && prev_idx < hillslope_ran_history.size()
                         && hillslope_ran_history[prev_idx] != 0)
                            ? &surface_history[prev_idx]
                            : nullptr;

                    lsdtt_raster_cpu replay_delivered = result.elevation;
                    const hillslope_outcome replay_hs = run_hillslope_solve(
                        replay_flow,
                        best_cycle,
                        result.fluvial_only,
                        result.is_channel,
                        replay_warm_start,
                        replay_delivered);
                    if (replay_hs.ran)
                        result.elevation = replay_delivered;
                }
            }

            result.delivered_cycle = best_cycle;
            result.converged = false;

            oss << "Limit cycle detected at cycle " << last_executed_cycle
                << ": flow graph repeats cycle " << result.limit_cycle_first
                << " (limit_cycle_period=" << result.limit_cycle_period
                << ", " << result.limit_cycle_n_pixels
                << " pixel(s) change receiver within the orbit)."
                << " Not converged. Delivered cycle " << best_cycle
                << ", the orbit member with the smallest divide step_p90="
                << best_p90 << " m";
        }

        result.stop_reason = oss.str();
    }

    if (!result.converged && result.stop_reason.empty())
    {
        if (!opts.convergence.check_convergence)
        {
            result.stop_reason = "Completed max_cycles with convergence checking disabled.";
        }
        else if (!result.cycles.empty())
        {
            result.stop_reason = make_non_converged_reason(
                static_cast<int>(result.cycles.size()),
                result.cycles.back(),
                opts.convergence);
        }
        else
        {
            result.stop_reason = "No cycles executed.";
        }
    }

    // hillslope_every_cycle == false: one solve on the surface actually
    // delivered, after the cycle loop (and after any best-of-orbit replay) has
    // decided which surface that is. Faster, but it never proves the divide is
    // stable under re-routing, which is the whole point of solving every cycle.
    if (opts.run_hillslope && !opts.hillslope_every_cycle
        && !result.cycles.empty())
    {
        result.fluvial_only = result.elevation;
        result.is_channel = make_is_channel(
            result.elevation, result.drainage_area, nullptr);

        lsdtt_flowinfo final_flow(
            result.fluvial_only.get_data_array(),
            result.fluvial_only.get_mask_array(),
            result.fluvial_only.get_meta(),
            fi_opts);

        lsdtt_raster_cpu final_delivered = result.elevation;
        const hillslope_outcome final_hs = run_hillslope_solve(
            final_flow,
            std::max(1, result.delivered_cycle),
            result.fluvial_only,
            result.is_channel,
            nullptr,
            final_delivered);

        if (final_hs.ran)
        {
            result.elevation = final_delivered;
            result.hillslope_ran = true;

            lsdtt_snap_cycle_record& last = result.cycles.back();
            last.hillslope_iterations = final_hs.iterations;
            last.hillslope_residual = static_cast<double>(final_hs.residual);
            last.hillslope_converged = final_hs.converged ? 1 : 0;
            last.n_channel_pixels = std::count(
                result.is_channel.begin(), result.is_channel.end(),
                static_cast<std::uint8_t>(1));
        }
    }

    if (opts.convergence.divide_gradient_target > 0.0
        && !result.cycles.empty())
    {
        const lsdtt_snap_cycle_record& last = result.cycles.back();
        const bool target_met =
            last.divide.gradient_median <= opts.convergence.divide_gradient_target;
        std::ostringstream oss;
        oss << " Divide-gradient target check: final median="
            << last.divide.gradient_median
            << ", target<=" << opts.convergence.divide_gradient_target
            << " (" << (target_met ? "met" : "not met") << ").";
        result.stop_reason += oss.str();
    }

    if (hillslope_requested)
    {
        std::ostringstream oss;
        oss << " Hillslope solve: " << opts.hillslope_model
            << (opts.hillslope_every_cycle ? ", every cycle"
                                           : ", final cycle only");
        if (!result.hillslope_ran)
        {
            oss << " -- REQUESTED BUT NOT RUN: no channel pixels were"
                   " identified, so the elliptic problem had no Dirichlet set."
                   " Set channel_threshold_area.";
        }
        else if (!result.cycles.empty())
        {
            const lsdtt_snap_cycle_record& last = result.cycles.back();
            oss << " (" << last.n_channel_pixels << " channel pixels, "
                << last.hillslope_iterations << " iterations, residual "
                << last.hillslope_residual
                << (last.hillslope_converged == 1 ? ", converged"
                                                  : ", DID NOT CONVERGE")
                << ")";
        }
        oss << ".";

        result.stop_reason += oss.str();

        if (!preflight_note.empty())
            result.stop_reason += preflight_note;

        if (result.hillslope_bc_mismatch)
        {
            result.stop_reason +=
                make_hillslope_bc_warning(opts.boundary_conditions);
        }
    }

    return result;
}

} // namespace lsdtt
