// lsdtt_lem_cpu.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html

#include "lsdtt_lem_cpu.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <numeric>
#include <fstream>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <vector>

#include "lsdtt_flowgraph_cpu.hpp"
#include "lsdtt_hillslope_implicit_cpu.hpp"
#include "lsdtt_hillslope_newton_cpu.hpp"
#include "lsdtt_hillslope_transport_cpu.hpp"
#include "lsdtt_stream_power_cpu.hpp"
#include <iostream>

namespace lsdtt
{

const char* name_of(hillslope_coupling_regime r) noexcept
{
    switch (r)
    {
        case hillslope_coupling_regime::automatic:              return "auto";
        case hillslope_coupling_regime::coupled:                return "coupled";
        case hillslope_coupling_regime::channels_as_boundaries: return "channels_as_boundaries";
    }
    return "unknown";
}

hillslope_coupling_regime parse_hillslope_coupling_regime(const std::string& name)
{
    std::string k;
    for (char ch : name)
        if (!std::isspace(static_cast<unsigned char>(ch)))
            k += static_cast<char>(std::tolower(static_cast<unsigned char>(ch)));

    if (k == "auto" || k == "automatic")   return hillslope_coupling_regime::automatic;
    if (k == "coupled")                    return hillslope_coupling_regime::coupled;
    if (k == "channels_as_boundaries" || k == "channels_as_boundary")
        return hillslope_coupling_regime::channels_as_boundaries;

    throw std::runtime_error(
        "unrecognised hillslope_channel_coupling '" + name +
        "' - expected one of: auto, coupled, channels_as_boundaries");
}

hillslope_coupling_regime resolve_hillslope_coupling_regime(
    hillslope_coupling_regime requested,
    float                     incision_threshold_area,
    bool                      has_explicit_channel_mask) noexcept
{
    if (requested != hillslope_coupling_regime::automatic) return requested;

    // The historical inference, unchanged: a channel set exists if a threshold
    // was given or a mask was supplied, and if one exists the channels are
    // Dirichlet nodes.
    return (incision_threshold_area > 0.0f || has_explicit_channel_mask)
        ? hillslope_coupling_regime::channels_as_boundaries
        : hillslope_coupling_regime::coupled;
}


namespace
{

// Shared core of the channel-mask logic: the explicit mask when the caller
// supplied one, else a drainage-area threshold on an already-built flow graph.
// `flow` may be null when an explicit mask makes the flow graph unnecessary.
std::vector<uint8_t>
build_channel_mask(
    const lsdtt_flowgraph_cpu* flow,
    const lsdtt_raster_meta&   meta,
    const lsdtt_lem_options&   opts)
{
    const std::int64_t n_cells =
        static_cast<std::int64_t>(meta.n_rows) * meta.n_cols;
    const float dx = static_cast<float>(meta.cell_size);

    std::vector<uint8_t> is_channel(static_cast<std::size_t>(n_cells), 0);

    // `coupled` means drainage density EMERGES, so there is no Dirichlet set --
    // whatever incision_threshold_area or an explicit mask say. Returning early
    // here is what makes the parameter override the threshold rather than
    // merely describe it.
    if (resolve_hillslope_coupling_regime(
            opts.coupling_regime, opts.incision_threshold_area,
            !opts.explicit_channel_mask.empty()) ==
        hillslope_coupling_regime::coupled)
    {
        return is_channel;
    }

    if (opts.explicit_channel_mask.size() == static_cast<std::size_t>(n_cells))
    {
        for (std::int64_t flat = 0; flat < n_cells; ++flat)
            is_channel[static_cast<std::size_t>(flat)] =
                opts.explicit_channel_mask[static_cast<std::size_t>(flat)] ? 1 : 0;
    }
    else if (opts.incision_threshold_area > 0.0f && flow != nullptr)
    {
        // The threshold convention lives in one place now
        // (lsdtt_stream_power_cpu.hpp): `incision_threshold_area` is a SPECIFIC
        // area in metres when divide_area_by_dx is set, and a map area in m^2
        // when it is not.
        std::vector<double> area = flow->accumulate();
        for (std::int64_t flat = 0; flat < n_cells; ++flat)
        {
            if (!flow->is_valid_flat(flat)) continue;
            const double q = channel_threshold_quantity(
                area[static_cast<std::size_t>(flat)],
                static_cast<double>(dx), opts.divide_area_by_dx);
            if (q >= static_cast<double>(opts.incision_threshold_area))
                is_channel[static_cast<std::size_t>(flat)] = 1;
        }
    }
    return is_channel;
}

// The D-infinity flow graph the hillslope step builds for channel identification.
lsdtt_flowgraph_cpu
build_hillslope_flowgraph(
    const lsdtt_raster_cpu&  dem,
    const lsdtt_lem_options& opts)
{
    lsdtt_flowgraph_options fg_opts;
    fg_opts.boundary_conditions       = opts.boundary_conditions;
    fg_opts.interior_nodata_is_outlet = true;
    fg_opts.method                    = lsdtt_flowgraph_method::d_infinity;
    fg_opts.n_threads                 = opts.n_threads;

    return lsdtt_flowgraph_cpu(
        dem.get_data_array(),
        dem.get_mask_array(),
        dem.get_meta(),
        fg_opts
    );
}

} // namespace

lsdtt_lem_forcing_point
lsdtt_lem_forcing_history::at(float time_yr) const
{
    if (points.empty()) return lsdtt_lem_forcing_point{};

    // Piecewise constant: the last row whose time_yr <= t.  Before the first
    // row, the first row applies.
    lsdtt_lem_forcing_point in_force = points.front();
    for (const lsdtt_lem_forcing_point& p : points)
    {
        if (p.time_yr > time_yr) break;
        in_force = p;
    }
    return in_force;
}

lsdtt_lem_forcing_history
load_lem_forcing_history(const std::string& csv_path)
{
    std::ifstream in(csv_path);
    if (!in.is_open())
        throw std::runtime_error("could not open forcing history CSV: " + csv_path);

    auto split = [](const std::string& line)
    {
        std::vector<std::string> cells;
        std::stringstream ss(line);
        std::string cell;
        while (std::getline(ss, cell, ','))
        {
            // trim whitespace
            const auto b = cell.find_first_not_of(" \t\r\n");
            const auto e = cell.find_last_not_of(" \t\r\n");
            cells.push_back(b == std::string::npos
                                ? std::string()
                                : cell.substr(b, e - b + 1));
        }
        return cells;
    };

    std::string line;
    if (!std::getline(in, line))
        throw std::runtime_error("forcing history CSV is empty: " + csv_path);

    // Column order is taken from the header, so columns may appear in any order.
    const std::vector<std::string> header = split(line);
    int col_time = -1, col_uplift = -1, col_fall = -1;
    for (std::size_t i = 0; i < header.size(); ++i)
    {
        if (header[i] == "time_yr")             col_time   = static_cast<int>(i);
        else if (header[i] == "uplift_rate")    col_uplift = static_cast<int>(i);
        else if (header[i] == "baselevel_fall_rate") col_fall = static_cast<int>(i);
    }
    if (col_time < 0)
        throw std::runtime_error(
            "forcing history CSV needs a time_yr column: " + csv_path);

    lsdtt_lem_forcing_history hist;
    while (std::getline(in, line))
    {
        if (line.find_first_not_of(" \t\r\n") == std::string::npos) continue;

        const std::vector<std::string> cells = split(line);
        auto cell_value = [&](int col) -> float
        {
            if (col < 0 || col >= static_cast<int>(cells.size())) return 0.0f;
            if (cells[static_cast<std::size_t>(col)].empty())      return 0.0f;
            return std::stof(cells[static_cast<std::size_t>(col)]);
        };

        lsdtt_lem_forcing_point p;
        p.time_yr             = cell_value(col_time);
        p.uplift_rate         = cell_value(col_uplift);
        p.baselevel_fall_rate = cell_value(col_fall);
        hist.points.push_back(p);
    }

    std::sort(hist.points.begin(), hist.points.end(),
              [](const lsdtt_lem_forcing_point& a, const lsdtt_lem_forcing_point& b)
              { return a.time_yr < b.time_yr; });

    return hist;
}

std::vector<std::uint8_t>
compute_lem_channel_mask(
    const lsdtt_raster_cpu&  dem,
    const lsdtt_lem_options& opts)
{
    // An explicit mask needs no flow graph at all.
    const std::size_t n_cells =
        static_cast<std::size_t>(dem.get_meta().n_rows) * dem.get_meta().n_cols;
    if (opts.explicit_channel_mask.size() == n_cells)
        return opts.explicit_channel_mask;

    lsdtt_flowgraph_cpu flow = build_hillslope_flowgraph(dem, opts);
    return build_channel_mask(&flow, dem.get_meta(), opts);
}

std::pair<lsdtt_raster_cpu, lsdtt_lem_step_stats>
run_lem_step(
    const lsdtt_raster_cpu&           dem,
    const lsdtt_lem_options&          opts,
    const std::vector<std::int64_t>&  boundary_flats,
    lsdtt_hillslope_newton_workspace* hillslope_ws,
    lsdtt_hillslope_cellular_workspace* cellular_ws)
{
    const lsdtt_raster_meta& meta = dem.get_meta();
    const std::size_t n_cells = static_cast<std::size_t>(meta.n_rows) * meta.n_cols;

    // Save pre-step elevations of boundary pixels
    const std::vector<float>& z_pre = dem.get_data();
    std::vector<float> boundary_elev(boundary_flats.size());
    for (std::size_t i = 0; i < boundary_flats.size(); ++i)
    {
        const auto flat = static_cast<std::size_t>(boundary_flats[i]);
        if (flat < n_cells)
            boundary_elev[i] = z_pre[flat];
    }

    // Snapshot of pre-step elevations for statistics (valid cells only)
    const std::vector<float> z_before = z_pre;

    lsdtt_raster_cpu current = dem;

    // -------------------------------------------------------------------------
    // 1. Fluvial incision
    // -------------------------------------------------------------------------
    if (opts.enable_d8_fluvial_erosion)
    {
        // Build D8 flow graph (only when needed)
        lsdtt_flowinfo_options fi_opts;
        fi_opts.boundary_conditions       = opts.boundary_conditions;
        fi_opts.interior_nodata_is_outlet = true;
        fi_opts.allow_pits_as_outlets     = true; // LEM surface may have transient pits

        lsdtt_flowinfo flow_d8(
            dem.get_data_array(),
            dem.get_mask_array(),
            dem.get_meta(),
            fi_opts
        );

        current = compute_d8_stream_power_incision(
            current, flow_d8,
            opts.erodibility,
            opts.m_exponent,
            opts.n_exponent,
            opts.timestep,
            opts.incision_threshold_area,
            opts.divide_area_by_dx,
            opts.slope_clamp
        );
    }
    else if (opts.enable_dinf_fluvial_erosion)
    {
        // Build D-infinity flow graph for fluvial incision
        lsdtt_flowgraph_options fg_opts;
        fg_opts.boundary_conditions       = opts.boundary_conditions;
        fg_opts.interior_nodata_is_outlet = true;
        fg_opts.method                    = lsdtt_flowgraph_method::d_infinity;
        fg_opts.n_threads                 = opts.n_threads;

        lsdtt_flowgraph_cpu flow_dinf(
            dem.get_data_array(),
            dem.get_mask_array(),
            dem.get_meta(),
            fg_opts
        );

        current = compute_dinf_stream_power_incision(
            current, flow_dinf,
            opts.erodibility,
            opts.m_exponent,
            opts.n_exponent,
            opts.timestep,
            opts.incision_threshold_area,
            opts.divide_area_by_dx,
            opts.slope_clamp,
            opts.n_threads
        );
    }

    // -------------------------------------------------------------------------
    // 2. Hillslope transport
    //    Built from the post-incision topography.
    //
    //    Cidre (default): builds a D-inf flow graph and uses the topological
    //      forward sweep for both routing and channel identification.
    //    Roering: builds a D-inf flow graph only for channel identification,
    //      then applies an explicit FD stencil (no topological ordering needed).
    // -------------------------------------------------------------------------
    if (opts.enable_hillslope_transport)
    {

        if (opts.hillslope_transport_model == "adi") 
        {
            // Alternating directions implicit does not need to build a flow graph
            current = compute_adi_hillslope_transport(
                current,
                opts.hillslope_diffusivity,
                opts.timestep,
                opts.boundary_conditions
            );
        }
        else
        {
            // Cidre routes sediment on the flow graph, so it always needs one.
            // The explicit-FD models (linear / roering / roering_implicit) use the
            // flow graph only to identify channels - which an explicit channel mask
            // already does - so skip the rebuild in that case.  It is otherwise
            // rebuilt every step, and on a long fixed-network run that dominates.
            const bool have_explicit_mask =
                opts.explicit_channel_mask.size() == n_cells;
            const bool is_cidre = (opts.hillslope_transport_model == "cidre");
            const bool is_cellular =
                (opts.hillslope_transport_model == "cellular_experimental");
            const bool coupled = resolve_hillslope_coupling_regime(
                opts.coupling_regime, opts.incision_threshold_area,
                have_explicit_mask) == hillslope_coupling_regime::coupled;

            std::optional<lsdtt_flowgraph_cpu> flow_hillslope;
            if (is_cidre || (!have_explicit_mask && !(is_cellular && coupled)))
                flow_hillslope.emplace(build_hillslope_flowgraph(current, opts));

            // Build the Dirichlet channel mask once for the channel-using explicit
            // FD models: prefer an explicit channel mask when supplied (Sec8), else
            // derive it from the drainage-area threshold.
            auto build_is_channel = [&]() -> std::vector<uint8_t>
            {
                return build_channel_mask(
                    flow_hillslope ? &*flow_hillslope : nullptr, meta, opts);
            };

            if (is_cellular)
            {
                std::vector<uint8_t> is_channel = build_is_channel();
                lsdtt_hillslope_cellular_options cellular_options;
                cellular_options.D = opts.hillslope_diffusivity;
                cellular_options.U = 0.0;
                cellular_options.Sc = opts.critical_hillslope_slope;
                cellular_options.r_match = opts.cellular_r_match;
                cellular_options.cfl_safety = opts.cellular_cfl_safety;
                cellular_options.boundary_conditions = opts.boundary_conditions;
                cellular_options.n_threads = opts.n_threads;
                const lsdtt_hillslope_cellular_result cellular =
                    compute_hillslope_cellular_transport(
                        current, is_channel, opts.timestep,
                        cellular_options, cellular_ws);
                if (!cellular.completed)
                    throw std::runtime_error(
                        "experimental cellular hillslope step failed: "
                        + cellular.message);
                current = cellular.elevation;
            }
            else if (opts.hillslope_transport_model == "linear" ||
                opts.hillslope_transport_model == "roering")
            {
                // Explicit FD on the regular grid - flow graph used only for channel ID.
                std::vector<uint8_t> is_channel = build_is_channel();

                const bool nl = (opts.hillslope_transport_model == "roering")
                                && opts.enable_nonlinear_hillslopes;

                current = compute_roering_hillslope_transport(
                    current,
                    is_channel,
                    opts.hillslope_diffusivity,
                    opts.timestep,
                    nl,
                    opts.critical_hillslope_slope,
                    opts.boundary_conditions,
                    opts.n_threads
                );
            }
            else if (opts.hillslope_transport_model == "roering_implicit")
            {
                // Unconditionally-stable backward-Euler Roering solver (Sec8).
                std::vector<uint8_t> is_channel = build_is_channel();

                // NOTE (2026-08-27): this step used to be SKIPPED whenever the
                // channel set was empty, on the reasoning that "without Dirichlet
                // channel sinks the solve degenerates to sink-free diffusion
                // toward the edges, which is meaningless".
                //
                // That reasoning applies to a STEADY solve with no Dirichlet data
                // anywhere -- an elliptic problem with pure Neumann data is
                // singular, and the field drifts without bound. It does NOT apply
                // to a transient step, which is what the guard was blocking:
                //
                //   * a diffusion step over dt is well posed under any boundary
                //     condition; mass simply redistributes;
                //   * in the COUPLED regime the fluvial term is the sink, applied
                //     to every cell in the same timestep;
                //   * an 'o' edge supplies base level regardless of channels.
                //
                // The empty channel set is exactly what `coupling_regime ==
                // coupled` produces on purpose, so the guard was disabling the
                // emergent-drainage-density mode it was meant to support. The
                // hillslope operator needs BOUNDARY CONDITIONS, not channels.
                {
                    lsdtt_roering_implicit_options ri = opts.roering_implicit_opts;
                    ri.D  = opts.hillslope_diffusivity;
                    ri.U  = 0.0f;                       // uplift applied in step 3
                    ri.Sc = opts.critical_hillslope_slope;
                    ri.nonlinear = opts.enable_nonlinear_hillslopes;
                    ri.boundary_conditions = opts.boundary_conditions;
                    ri.n_threads = opts.n_threads;

                    // Transient-coupling rule 1 (Sec8): integrate exactly one LEM step
                    // (dt_total = opts.timestep), so channels are re-incised each step.
                    //
                    // NEWTON FIRST. The Picard path lags the face coefficients, so
                    // near Sc it limit-cycles and needs heavy under-relaxation;
                    // Newton differentiates alpha and converges in an iteration
                    // count that does not grow with the grid. Picard is kept as the
                    // automatic fallback -- it is the only option in an Eigen-free
                    // build, and a step Newton cannot take is worth retrying with a
                    // damped method rather than failing the run.
                    bool done = false;
                    if (!opts.hillslope_force_picard)
                    {
                        lsdtt_hillslope_newton_transient_options nt;
                        nt.newton.D  = ri.D;
                        nt.newton.U  = ri.U;
                        nt.newton.Sc = ri.Sc;
                        nt.newton.nonlinear           = ri.nonlinear;
                        nt.newton.boundary_conditions = ri.boundary_conditions;
                        nt.newton.n_threads           = ri.n_threads;
                        nt.dt_init = (ri.dt_init > 0.0f) ? ri.dt_init : opts.timestep;

                        // The workspace carries the factorisation between
                        // timesteps. For a LINEAR run with a fixed dt the matrix is
                        // identical every step, so this turns a per-step
                        // assemble-analyse-factorise-solve into a solve.
                        const lsdtt_hillslope_newton_transient_result nr =
                            compute_hillslope_newton_transient(
                                current, is_channel, opts.timestep, nt,
                                hillslope_ws);

                        if (nr.converged)
                        {
                            current = nr.elevation;
                            done    = true;
                        }
                        else
                        {
                            std::cerr << "[lem] Newton hillslope step did not "
                                         "converge (" << nr.message
                                      << "); falling back to Picard.\n";
                        }
                    }

                    if (!done)
                    {
                        lsdtt_roering_implicit_result ri_res =
                            compute_roering_hillslope_transport_implicit(
                                current, is_channel, opts.timestep, ri);
                        current = std::move(ri_res.elevation);
                    }
                }
            }
            else  // "cidre"
            {
                current = compute_cidre_hillslope_transport(
                    current, *flow_hillslope,
                    opts.hillslope_diffusivity,
                    opts.timestep,
                    opts.enable_nonlinear_hillslopes,
                    opts.critical_hillslope_slope,
                    opts.incision_threshold_area,
                    opts.divide_area_by_dx,
                    opts.deposit_on_hillslope_baselevels,
                    opts.n_threads
                );
            }
        }
    }

    // -------------------------------------------------------------------------
    // 3. Uplift (uniform)
    // -------------------------------------------------------------------------
    const std::vector<uint8_t>& mask = current.get_mask();
    std::vector<float>& z = current.get_data();

    if (opts.enable_uplift)
    {
        const float uplift_dz = opts.uplift_rate * opts.timestep;
        for (std::size_t i = 0; i < n_cells; ++i)
        {
            if (mask[i])
                z[i] += uplift_dz;
        }
    }

    // -------------------------------------------------------------------------
    // 4. Statistics (computed before boundary restoration so they reflect the
    //    full-domain change)
    // -------------------------------------------------------------------------
    lsdtt_lem_step_stats stats{};
    // -------------------------------------------------------------------------
    // Restore the fixed-elevation boundaries BEFORE computing statistics.
    //
    // The order matters and used to be wrong. Base-level cells receive uplift like
    // everything else and are then reset to their pinned value, so they never
    // actually move -- but with the restore happening after the statistics, every
    // one of them contributed a spurious |dz| of exactly U*dt.
    //
    // On a 200 x 200 grid with two base-level rows that is 400 cells, a full 1% of
    // the domain, which pinned max|dz| AND the 99.9th percentile at U*dt for the
    // entire run. A steady-state test built on either could never be satisfied: a
    // 6 Myr run sat at "max|dz| = 0.4 m" while the landscape's real p99.9 was
    // 0.028 m and no cell anywhere moved by U*dt.
    // -------------------------------------------------------------------------
    const float baselevel_dz = opts.baselevel_fall_rate * opts.timestep;

    for (std::size_t i = 0; i < boundary_flats.size(); ++i)
    {
        const auto flat = static_cast<std::size_t>(boundary_flats[i]);
        if (flat < n_cells)
            z[flat] = boundary_elev[i] - baselevel_dz;
    }

    {
        double sum_before = 0.0, sum_after = 0.0;
        double max_before = -1e30, max_after = -1e30;
        double max_local  = 0.0;
        std::vector<float> local_changes;
        local_changes.reserve(n_cells);
        double sum_erosion = 0.0;
        double max_erosion = 0.0;
        std::int64_t n_valid = 0;

        const double inv_dt_rate = 1.0 / static_cast<double>(opts.timestep);
        if (opts.compute_erosion_rate_field)
            stats.erosion_rate.assign(n_cells, 0.0f);

        for (std::size_t i = 0; i < n_cells; ++i)
        {
            if (!mask[i])
                continue;
            ++n_valid;
            const double zb = static_cast<double>(z_before[i]);
            const double za = static_cast<double>(z[i]);

            sum_before += zb;
            sum_after  += za;
            max_before  = std::max(max_before, zb);
            max_after   = std::max(max_after,  za);

            const double local_change = std::abs(za - zb);
            max_local = std::max(max_local, local_change);
            local_changes.push_back(static_cast<float>(local_change));

            // Erosion = how much lower the surface is now (before uplift, but
            // we already applied uplift, so we back it out for the rate)
            const double uplift_dz = opts.enable_uplift
                ? static_cast<double>(opts.uplift_rate * opts.timestep)
                : 0.0;
            const double erosion = zb - (za - uplift_dz);  // positive = erosion

            // The per-pixel field keeps the sign: deposition is a real, negative
            // erosion rate, and clamping it away would hide sediment accumulation.
            if (opts.compute_erosion_rate_field)
                stats.erosion_rate[i] = static_cast<float>(erosion * inv_dt_rate);

            if (erosion > 0.0)
            {
                sum_erosion += erosion;
                max_erosion  = std::max(max_erosion, erosion);
            }
        }

        const double inv_dt = 1.0 / static_cast<double>(opts.timestep);
        const double inv_n  = n_valid > 0 ? 1.0 / static_cast<double>(n_valid) : 0.0;

        stats.delta_max_elevation        = static_cast<float>(std::abs(max_after - max_before));
        stats.delta_mean_elevation       = static_cast<float>(
            std::abs((sum_after - sum_before) * inv_n));
        stats.max_local_elevation_change = static_cast<float>(max_local);

        // p99.9 by partial sort: the whole array is never needed in order.
        if (!local_changes.empty())
        {
            const std::size_t k = static_cast<std::size_t>(
                0.999 * static_cast<double>(local_changes.size() - 1));
            std::nth_element(local_changes.begin(),
                             local_changes.begin() + static_cast<std::ptrdiff_t>(k),
                             local_changes.end());
            stats.p999_local_elevation_change = local_changes[k];
        }
        else stats.p999_local_elevation_change = 0.0f;
        stats.mean_erosion_rate          = static_cast<float>(sum_erosion * inv_n * inv_dt);
        stats.max_erosion_rate           = static_cast<float>(max_erosion * inv_dt);
    }

    // -------------------------------------------------------------------------
    // 5. Restore fixed-elevation boundary pixels.
    //
    //    With a base-level fall the boundary is not restored to its pre-step
    //    elevation but lowered by fall_rate * dt, so base level drops away from
    //    the domain and a wave of incision propagates up the channel network.
    // -------------------------------------------------------------------------
    return {std::move(current), stats};
}

} // namespace lsdtt
