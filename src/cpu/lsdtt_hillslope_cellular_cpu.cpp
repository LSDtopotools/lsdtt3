// lsdtt_hillslope_cellular_cpu.cpp

#include "lsdtt_hillslope_cellular_cpu.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

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
    ws.elevation.assign(dem.get_data().begin(), dem.get_data().end());
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

    auto compute_faces = [&](bool write_flux) -> face_scan
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
                        ++active;
                        if (ratio >= 0.70) ++r70;
                        if (ratio >= 0.80) ++r80;
                        if (ratio >= 0.85) ++r85;
                        if (ratio >= 1.00) ++r100;
                        max_r = std::max(max_r, ratio);
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
                        ++active;
                        if (ratio >= 0.70) ++r70;
                        if (ratio >= 0.80) ++r80;
                        if (ratio >= 0.85) ++r85;
                        if (ratio >= 1.00) ++r100;
                        max_r = std::max(max_r, ratio);
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

    const face_scan initial = compute_faces(dt_total > 0.0);
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
        if (result.n_substeps > 0) scan = compute_faces(true);
        result.n_face_evaluations += scan.active;

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

    const face_scan final_scan = compute_faces(false);
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

} // namespace lsdtt
