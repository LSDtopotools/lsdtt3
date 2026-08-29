// lsdtt_hillslope_cellular_cpu.cpp

#include "lsdtt_hillslope_cellular_cpu.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "lsdtt_hillslope_operator_cpu.hpp"

namespace lsdtt
{
namespace
{

double mobility(double r, double r_match)
{
    if (r <= r_match) return 1.0 / (1.0 - r * r);
    const double d = 1.0 - r_match * r_match;
    const double m = 1.0 / d;
    const double mp = 2.0 * r_match / (d * d);
    return m + mp * (r - r_match);
}

double flux_tangent(double r, double r_match)
{
    if (r <= r_match)
    {
        const double d = 1.0 - r * r;
        return (1.0 + r * r) / (d * d);
    }
    const double d = 1.0 - r_match * r_match;
    const double mp = 2.0 * r_match / (d * d);
    return mobility(r, r_match) + r * mp;
}

struct face_scan
{
    long active = 0;
    long r70 = 0;
    long r80 = 0;
    long r85 = 0;
    long r100 = 0;
    double max_r = 0.0;
    double max_tangent = 1.0;
};

} // namespace

lsdtt_hillslope_cellular_result compute_hillslope_cellular_transport(
    const lsdtt_raster_cpu& dem,
    const std::vector<std::uint8_t>& is_channel,
    double dt_total,
    const lsdtt_hillslope_cellular_options& options,
    lsdtt_hillslope_cellular_workspace* workspace)
{
    if (!(options.D > 0.0) || !std::isfinite(options.D))
        throw std::invalid_argument("cellular hillslope D must be finite and > 0");
    if (!std::isfinite(options.U))
        throw std::invalid_argument("cellular hillslope U must be finite");
    if (!(options.Sc > 0.0) || !std::isfinite(options.Sc))
        throw std::invalid_argument("cellular hillslope Sc must be finite and > 0");
    if (!(options.r_match > 0.0 && options.r_match < 1.0))
        throw std::invalid_argument("cellular hillslope r_match must be in (0, 1)");
    if (!(options.cfl_safety > 0.0 && options.cfl_safety <= 1.0))
        throw std::invalid_argument("cellular hillslope cfl_safety must be in (0, 1]");
    if (!(dt_total >= 0.0) || !std::isfinite(dt_total))
        throw std::invalid_argument("cellular hillslope dt_total must be finite and >= 0");
    if (options.max_substeps <= 0)
        throw std::invalid_argument("cellular hillslope max_substeps must be > 0");

    const lsdtt_raster_meta& meta = dem.get_meta();
    const int nr = meta.n_rows;
    const int nc = meta.n_cols;
    const std::size_t n = static_cast<std::size_t>(nr) * nc;
    if (!is_channel.empty() && is_channel.size() != n)
        throw std::invalid_argument("cellular hillslope channel mask size mismatch");

    const lsdtt_hillslope_domain dom = lsdtt_hillslope_domain::make(
        nr, nc, options.boundary_conditions, dem.get_mask());
    if ((dom.bc_n == 'p') != (dom.bc_s == 'p')
        || (dom.bc_e == 'p') != (dom.bc_w == 'p'))
        throw std::invalid_argument("cellular hillslope periodic boundaries must be paired");

#ifdef _OPENMP
    if (options.n_threads > 0) omp_set_num_threads(options.n_threads);
#else
    (void)options.n_threads;
#endif

    lsdtt_hillslope_cellular_workspace local;
    lsdtt_hillslope_cellular_workspace& ws = workspace ? *workspace : local;
    ws.elevation.resize(n);
#ifdef _OPENMP
    #pragma omp parallel for schedule(static)
#endif
    for (long long i = 0; i < static_cast<long long>(n); ++i)
        ws.elevation[static_cast<std::size_t>(i)] =
            static_cast<double>(dem.get_data()[static_cast<std::size_t>(i)]);
    ws.elevation_new.resize(n);
    ws.q_east.resize(n);
    ws.q_south.resize(n);

    auto fixed = [&](std::size_t i, int r, int c) -> bool
    {
        return !dem.get_mask()[i]
            || (!is_channel.empty() && is_channel[i] != 0u)
            || dom.is_dirichlet_edge(r, c);
    };

    const bool periodic_ew = dom.bc_e == 'p' && dom.bc_w == 'p';
    const bool periodic_ns = dom.bc_n == 'p' && dom.bc_s == 'p';

    auto compute_faces = [&](bool write_flux, bool collect_diagnostics) -> face_scan
    {
        long active = 0, r70 = 0, r80 = 0, r85 = 0, r100 = 0;
        double max_r = 0.0, max_tangent = 1.0;
        const lsdtt_hillslope_stencil<double> st(dom, ws.elevation.data(), meta.cell_size);

#ifdef _OPENMP
        #pragma omp parallel for schedule(static) reduction(+:active,r70,r80,r85,r100) reduction(max:max_r,max_tangent)
#endif
        for (long long flat = 0; flat < static_cast<long long>(n); ++flat)
        {
            const std::size_t i = static_cast<std::size_t>(flat);
            const int r = static_cast<int>(flat / nc);
            const int c = static_cast<int>(flat % nc);
            double qe = 0.0, qs = 0.0;

            if (dem.get_mask()[i])
            {
                const bool have_e = c + 1 < nc || (c == nc - 1 && periodic_ew);
                if (have_e)
                {
                    int er = r, ec = c;
                    dom.neighbour(r, c, 0, 1, er, ec);
                    const std::size_t ei = dom.index(er, ec);
                    if (ei != i && dem.get_mask()[ei]
                        && !(fixed(i, r, c) && fixed(ei, er, ec)))
                    {
                        const double sn = (st.at(er, ec) - st.at(r, c)) / meta.cell_size;
                        const double stg = 0.5 * (st.slope_y(r, c) + st.slope_y(er, ec));
                        const double ratio = std::sqrt(sn * sn + stg * stg) / options.Sc;
                        qe = -options.D * mobility(ratio, options.r_match) * sn;
                        if (collect_diagnostics)
                        {
                            ++active;
                            if (ratio >= 0.70) ++r70;
                            if (ratio >= 0.80) ++r80;
                            if (ratio >= 0.85) ++r85;
                            if (ratio >= 1.00) ++r100;
                            max_r = std::max(max_r, ratio);
                        }
                        max_tangent = std::max(max_tangent,
                                               flux_tangent(ratio, options.r_match));
                    }
                }

                const bool have_s = r + 1 < nr || (r == nr - 1 && periodic_ns);
                if (have_s)
                {
                    int sr = r, sc = c;
                    dom.neighbour(r, c, 1, 0, sr, sc);
                    const std::size_t si = dom.index(sr, sc);
                    if (si != i && dem.get_mask()[si]
                        && !(fixed(i, r, c) && fixed(si, sr, sc)))
                    {
                        const double sn = (st.at(sr, sc) - st.at(r, c)) / meta.cell_size;
                        const double stg = 0.5 * (st.slope_x(r, c) + st.slope_x(sr, sc));
                        const double ratio = std::sqrt(sn * sn + stg * stg) / options.Sc;
                        qs = -options.D * mobility(ratio, options.r_match) * sn;
                        if (collect_diagnostics)
                        {
                            ++active;
                            if (ratio >= 0.70) ++r70;
                            if (ratio >= 0.80) ++r80;
                            if (ratio >= 0.85) ++r85;
                            if (ratio >= 1.00) ++r100;
                            max_r = std::max(max_r, ratio);
                        }
                        max_tangent = std::max(max_tangent,
                                               flux_tangent(ratio, options.r_match));
                    }
                }
            }
            if (write_flux)
            {
                ws.q_east[i] = qe;
                ws.q_south[i] = qs;
            }
        }
        return {active, r70, r80, r85, r100, max_r, max_tangent};
    };

    lsdtt_hillslope_cellular_result result;
    result.min_substep = dt_total > 0.0
        ? std::numeric_limits<double>::infinity() : 0.0;

    const face_scan initial = compute_faces(dt_total > 0.0, true);
    result.initial_active_faces = initial.active;
    result.initial_faces_r70 = initial.r70;
    result.initial_faces_r80 = initial.r80;
    result.initial_faces_r85 = initial.r85;
    result.initial_faces_r100 = initial.r100;
    result.initial_max_slope_ratio = initial.max_r;

    double remaining = dt_total;
    face_scan scan = initial;
    while (remaining > 0.0 && result.n_substeps < options.max_substeps)
    {
        if (result.n_substeps > 0) scan = compute_faces(true, false);
        result.n_face_evaluations += initial.active;

        const double stable_dt = options.cfl_safety * meta.cell_size * meta.cell_size
                               / (4.0 * options.D * scan.max_tangent);
        const double dt = std::min(remaining, stable_dt);
        if (!(dt > 0.0) || !std::isfinite(dt)) break;
        result.min_substep = std::min(result.min_substep, dt);
        result.max_substep = std::max(result.max_substep, dt);

        double exported_height_sum = 0.0;
#ifdef _OPENMP
        #pragma omp parallel for schedule(static) reduction(+:exported_height_sum)
#endif
        for (long long flat = 0; flat < static_cast<long long>(n); ++flat)
        {
            const std::size_t i = static_cast<std::size_t>(flat);
            const int r = static_cast<int>(flat / nc);
            const int c = static_cast<int>(flat % nc);
            ws.elevation_new[i] = ws.elevation[i];
            if (fixed(i, r, c)) continue;

            const double qe = ws.q_east[i];
            const double qs = ws.q_south[i];
            const double qw = c > 0
                ? ws.q_east[i - 1]
                : (periodic_ew ? ws.q_east[static_cast<std::size_t>(r) * nc + nc - 1] : 0.0);
            const double qn = r > 0
                ? ws.q_south[i - static_cast<std::size_t>(nc)]
                : (periodic_ns ? ws.q_south[static_cast<std::size_t>(nr - 1) * nc + c] : 0.0);
            ws.elevation_new[i] += dt * (options.U + (qw - qe + qn - qs) / meta.cell_size);
            // Sum the small per-cell balance rather than subtracting two large
            // whole-domain elevation sums. The latter loses precision even
            // though every internal face is exactly equal and opposite.
            exported_height_sum += ws.elevation[i] + options.U * dt
                                 - ws.elevation_new[i];
        }
        result.net_dirichlet_export_volume +=
            exported_height_sum * meta.cell_size * meta.cell_size;
        ws.elevation.swap(ws.elevation_new);
        ++result.n_substeps;
        if (dt >= remaining) remaining = 0.0;
        else remaining -= dt;
    }

    const face_scan final_scan = compute_faces(false, true);
    result.final_max_slope_ratio = final_scan.max_r;
    result.completed = remaining == 0.0;
    result.message = result.completed ? "completed" : "maximum substeps reached";

    std::vector<float> output(n);
#ifdef _OPENMP
    #pragma omp parallel for schedule(static)
#endif
    for (long long i = 0; i < static_cast<long long>(n); ++i)
        output[static_cast<std::size_t>(i)] =
            static_cast<float>(ws.elevation[static_cast<std::size_t>(i)]);
    result.elevation = lsdtt_raster_cpu(
        std::move(output), std::vector<std::uint8_t>(dem.get_mask()), meta);
    return result;
}

lsdtt_hillslope_cellular_result compute_hillslope_cellular_transport_multirate(
    const lsdtt_raster_cpu& dem,
    const std::vector<std::uint8_t>& is_channel,
    double dt_total,
    const lsdtt_hillslope_cellular_options& options,
    lsdtt_hillslope_cellular_workspace* workspace)
{
    if (!(options.D > 0.0) || !std::isfinite(options.D))
        throw std::invalid_argument("cellular multirate D must be finite and > 0");
    if (options.U != 0.0)
        throw std::invalid_argument("cellular multirate requires U=0; uplift is a separate LEM process");
    if (!(options.Sc > 0.0) || !std::isfinite(options.Sc))
        throw std::invalid_argument("cellular multirate Sc must be finite and > 0");
    if (!(options.r_match > 0.0 && options.r_match < 1.0))
        throw std::invalid_argument("cellular multirate r_match must be in (0, 1)");
    if (!(options.cfl_safety > 0.0 && options.cfl_safety <= 1.0))
        throw std::invalid_argument("cellular multirate cfl_safety must be in (0, 1]");
    if (!(dt_total >= 0.0) || !std::isfinite(dt_total))
        throw std::invalid_argument("cellular multirate dt_total must be finite and >= 0");
    if (options.max_multirate_level < 0 || options.max_multirate_level > 20)
        throw std::invalid_argument("cellular multirate level must be in [0, 20]");

    const lsdtt_raster_meta& meta = dem.get_meta();
    const int nr = meta.n_rows;
    const int nc = meta.n_cols;
    const std::size_t n = static_cast<std::size_t>(nr) * nc;
    if (!is_channel.empty() && is_channel.size() != n)
        throw std::invalid_argument("cellular multirate channel mask size mismatch");

    const lsdtt_hillslope_domain dom = lsdtt_hillslope_domain::make(
        nr, nc, options.boundary_conditions, dem.get_mask());
    if ((dom.bc_n == 'p') != (dom.bc_s == 'p')
        || (dom.bc_e == 'p') != (dom.bc_w == 'p'))
        throw std::invalid_argument("cellular multirate periodic boundaries must be paired");

#ifdef _OPENMP
    if (options.n_threads > 0) omp_set_num_threads(options.n_threads);
#else
    (void)options.n_threads;
#endif

    lsdtt_hillslope_cellular_workspace local;
    lsdtt_hillslope_cellular_workspace& ws = workspace ? *workspace : local;
    ws.elevation.assign(dem.get_data().begin(), dem.get_data().end());
    ws.edge_a.clear();
    ws.edge_b.clear();
    ws.edge_row.clear();
    ws.edge_col.clear();
    ws.edge_direction.clear();
    ws.edge_colour.clear();
    ws.edge_level.clear();
    ws.edge_transfer.clear();
    const std::size_t edge_reserve = 2 * n;
    ws.edge_a.reserve(edge_reserve);
    ws.edge_b.reserve(edge_reserve);
    ws.edge_row.reserve(edge_reserve);
    ws.edge_col.reserve(edge_reserve);
    ws.edge_direction.reserve(edge_reserve);
    ws.edge_colour.reserve(edge_reserve);

    auto fixed = [&](std::size_t i, int r, int c) -> bool
    {
        return !dem.get_mask()[i]
            || (!is_channel.empty() && is_channel[i] != 0u)
            || dom.is_dirichlet_edge(r, c);
    };
    const bool periodic_ew = dom.bc_e == 'p' && dom.bc_w == 'p';
    const bool periodic_ns = dom.bc_n == 'p' && dom.bc_s == 'p';

    auto add_edge = [&](int r, int c, int dr, int dc,
                        std::uint8_t direction, std::uint8_t colour)
    {
        const std::size_t a = dom.index(r, c);
        int br = r, bc = c;
        dom.neighbour(r, c, dr, dc, br, bc);
        const std::size_t b = dom.index(br, bc);
        if (a == b || !dem.get_mask()[a] || !dem.get_mask()[b]) return;
        if (fixed(a, r, c) && fixed(b, br, bc)) return;
        ws.edge_a.push_back(a);
        ws.edge_b.push_back(b);
        ws.edge_row.push_back(r);
        ws.edge_col.push_back(c);
        ws.edge_direction.push_back(direction);
        ws.edge_colour.push_back(colour);
    };

    for (int r = 0; r < nr; ++r)
        for (int c = 0; c < nc; ++c)
        {
            if (c + 1 < nc)
                add_edge(r, c, 0, 1, 0u, static_cast<std::uint8_t>(c & 1));
            else if (periodic_ew)
                add_edge(r, c, 0, 1, 0u, 4u); // odd periodic seam is separate

            if (r + 1 < nr)
                add_edge(r, c, 1, 0, 1u, static_cast<std::uint8_t>(2 + (r & 1)));
            else if (periodic_ns)
                add_edge(r, c, 1, 0, 1u, 5u);
        }

    const std::size_t ne = ws.edge_a.size();
    ws.edge_level.resize(ne, 0u);
    ws.edge_transfer.resize(ne, 0.0);
    lsdtt_hillslope_cellular_result result;
    result.initial_active_faces = static_cast<long>(ne);
    if (dt_total == 0.0)
    {
        result.elevation = dem;
        result.completed = true;
        result.message = "completed";
        return result;
    }

    auto ratio_for_edge = [&](std::size_t e,
                              const lsdtt_hillslope_stencil<double>& st,
                              double* normal) -> double
    {
        const int r = ws.edge_row[e];
        const int c = ws.edge_col[e];
        int br = r, bc = c;
        if (ws.edge_direction[e] == 0u) dom.neighbour(r, c, 0, 1, br, bc);
        else                            dom.neighbour(r, c, 1, 0, br, bc);
        const double sn = (st.at(br, bc) - st.at(r, c)) / meta.cell_size;
        const double tang = ws.edge_direction[e] == 0u
            ? 0.5 * (st.slope_y(r, c) + st.slope_y(br, bc))
            : 0.5 * (st.slope_x(r, c) + st.slope_x(br, bc));
        if (normal) *normal = sn;
        return std::sqrt(sn * sn + tang * tang) / options.Sc;
    };

    // Assign with half the requested CFL safety, leaving one level of headroom
    // for local slope growth during the macro step. A runtime check below still
    // rejects the step if that headroom is insufficient.
    const lsdtt_hillslope_stencil<double> initial_st(dom, ws.elevation.data(), meta.cell_size);
    int max_level = 0;
    for (std::size_t e = 0; e < ne; ++e)
    {
        const double ratio = ratio_for_edge(e, initial_st, nullptr);
        result.initial_max_slope_ratio = std::max(result.initial_max_slope_ratio, ratio);
        if (ratio >= 0.70) ++result.initial_faces_r70;
        if (ratio >= 0.80) ++result.initial_faces_r80;
        if (ratio >= 0.85) ++result.initial_faces_r85;
        if (ratio >= 1.00) ++result.initial_faces_r100;
        const double stable = 0.5 * options.cfl_safety * meta.cell_size * meta.cell_size
                            / (4.0 * options.D * flux_tangent(ratio, options.r_match));
        int level = 0;
        while (level < options.max_multirate_level
               && dt_total / static_cast<double>(1u << level) > stable)
            ++level;
        if (dt_total / static_cast<double>(1u << level) > stable)
        {
            result.elevation = dem;
            result.message = "required multirate level exceeds configured maximum";
            return result;
        }
        ws.edge_level[e] = static_cast<std::uint8_t>(level);
        max_level = std::max(max_level, level);
    }
    result.max_multirate_level = max_level;

    const int levels = max_level + 1;
    std::vector<std::vector<std::size_t>> bins(static_cast<std::size_t>(6 * levels));
    for (std::size_t e = 0; e < ne; ++e)
        bins[static_cast<std::size_t>(ws.edge_colour[e]) * levels + ws.edge_level[e]].push_back(e);

    const unsigned fine_ticks = 1u << max_level;
    bool stable = true;
    for (unsigned tick = 0; tick < fine_ticks && stable; ++tick)
    {
        const lsdtt_hillslope_stencil<double> st(dom, ws.elevation.data(), meta.cell_size);

        // Phase 1: all due transfers read the same field. This prevents colour
        // order from changing fluxes and makes the parallel reads race-free.
        for (int level = 0; level <= max_level; ++level)
        {
            const unsigned stride = 1u << (max_level - level);
            if (tick % stride != 0u) continue;
            const double edge_dt = dt_total / static_cast<double>(1u << level);
            for (int colour = 0; colour < 6; ++colour)
            {
                const auto& bin = bins[static_cast<std::size_t>(colour) * levels + level];
                long violations = 0;
#ifdef _OPENMP
                #pragma omp parallel for schedule(static) reduction(+:violations)
#endif
                for (long long j = 0; j < static_cast<long long>(bin.size()); ++j)
                {
                    const std::size_t e = bin[static_cast<std::size_t>(j)];
                    double sn = 0.0;
                    const double ratio = ratio_for_edge(e, st, &sn);
                    const double stable_dt = options.cfl_safety * meta.cell_size * meta.cell_size
                        / (4.0 * options.D * flux_tangent(ratio, options.r_match));
                    if (edge_dt > stable_dt * (1.0 + 1e-12)) ++violations;
                    ws.edge_transfer[e] = -options.D * mobility(ratio, options.r_match)
                                        * sn * edge_dt / meta.cell_size;
                }
                result.stability_bin_violations += violations;
                result.n_face_evaluations += static_cast<long long>(bin.size());
                result.n_face_events += static_cast<long long>(bin.size());
            }
        }
        if (result.stability_bin_violations > 0)
        {
            stable = false;
            break;
        }

        // Phase 2: a colour contains no faces sharing a primary endpoint. Seam
        // colours are separate so odd-sized periodic grids remain conflict-free.
        for (int colour = 0; colour < 6; ++colour)
            for (int level = 0; level <= max_level; ++level)
            {
                const unsigned stride = 1u << (max_level - level);
                if (tick % stride != 0u) continue;
                const auto& bin = bins[static_cast<std::size_t>(colour) * levels + level];
                double export_sum = 0.0;
#ifdef _OPENMP
                #pragma omp parallel for schedule(static) reduction(+:export_sum)
#endif
                for (long long j = 0; j < static_cast<long long>(bin.size()); ++j)
                {
                    const std::size_t e = bin[static_cast<std::size_t>(j)];
                    const std::size_t a = ws.edge_a[e];
                    const std::size_t b = ws.edge_b[e];
                    const int r = ws.edge_row[e];
                    const int c = ws.edge_col[e];
                    int br = r, bc = c;
                    if (ws.edge_direction[e] == 0u) dom.neighbour(r, c, 0, 1, br, bc);
                    else                            dom.neighbour(r, c, 1, 0, br, bc);
                    const bool af = fixed(a, r, c);
                    const bool bf = fixed(b, br, bc);
                    const double dz = ws.edge_transfer[e];
                    if (!af) ws.elevation[a] -= dz;
                    else     export_sum -= dz;
                    if (!bf) ws.elevation[b] += dz;
                    else     export_sum += dz;
                }
                result.net_dirichlet_export_volume +=
                    export_sum * meta.cell_size * meta.cell_size;
            }
    }

    std::vector<float> output(n);
#ifdef _OPENMP
    #pragma omp parallel for schedule(static)
#endif
    for (long long i = 0; i < static_cast<long long>(n); ++i)
        output[static_cast<std::size_t>(i)] = static_cast<float>(ws.elevation[static_cast<std::size_t>(i)]);
    result.elevation = lsdtt_raster_cpu(
        std::move(output), std::vector<std::uint8_t>(dem.get_mask()), meta);

    const lsdtt_hillslope_stencil<double> final_st(dom, ws.elevation.data(), meta.cell_size);
    for (std::size_t e = 0; e < ne; ++e)
        result.final_max_slope_ratio = std::max(
            result.final_max_slope_ratio, ratio_for_edge(e, final_st, nullptr));
    result.n_substeps = static_cast<long>(fine_ticks);
    result.min_substep = dt_total / static_cast<double>(fine_ticks);
    result.max_substep = dt_total;
    result.completed = stable;
    result.message = stable ? "completed" : "multirate stability bin violated; reduce macro dt";
    return result;
}

} // namespace lsdtt
