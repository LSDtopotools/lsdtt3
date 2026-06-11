// lsdtt_hillslope_transport_cpu.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html

#include "lsdtt_hillslope_transport_cpu.hpp"

#include <algorithm>
#include <cmath>
#include <vector>

#ifdef _OPENMP
#  include <omp.h>
#endif

namespace lsdtt
{


// -----------------------------------------------------------------------------
// compute_cidre_hillslope_transport
// -----------------------------------------------------------------------------
lsdtt_raster_cpu compute_cidre_hillslope_transport(
    const lsdtt_raster_cpu&    dem,
    const lsdtt_flowgraph_cpu& flow,
    float kappa,
    float dt,
    bool  nonlinear,
    float S_c,
    float incision_threshold_area,
    bool  divide_area_by_dx,
    bool  deposit_on_hillslope_baselevels,
    int   n_threads)
{
#ifdef _OPENMP
    if (n_threads > 0)
        omp_set_num_threads(n_threads);
#else
    (void)n_threads;
#endif

    const lsdtt_raster_meta& meta = dem.get_meta();
    const float dx = static_cast<float>(meta.cell_size);
    const std::int64_t n_cells = static_cast<std::int64_t>(meta.n_rows) * meta.n_cols;

    // Working elevation array
    std::vector<float> z(dem.get_data());

    // Flux buffer (elevation units). Each cell accumulates incoming flux
    // contributed by its donors via atomic adds during the forward sweep.
    std::vector<float> flux(static_cast<std::size_t>(n_cells), 0.0f);

    // D-inf slopes from the flow graph
    const std::vector<float>& slopes = flow.dinf_slopes_raw();

    // Drainage area (for channel-pixel identification)
    std::vector<double> area_d;
    if (incision_threshold_area > 0.0f)
        area_d = flow.accumulate();

    const float Sc2              = S_c * S_c;
    const float kappa_dt_over_dx = kappa * dt / dx;

    // Per-cell channel flag (channel pixels are hillslope baselevels)
    std::vector<bool> is_channel(static_cast<std::size_t>(n_cells), false);
    if (incision_threshold_area > 0.0f)
    {
        for (std::int64_t flat = 0; flat < n_cells; ++flat)
        {
            if (!flow.is_valid_flat(flat))
                continue;
            double area = area_d[static_cast<std::size_t>(flat)];
            if (divide_area_by_dx)
                area /= static_cast<double>(dx);
            if (area >= static_cast<double>(incision_threshold_area))
                is_channel[static_cast<std::size_t>(flat)] = true;
        }
    }

    auto levels = build_forward_levels(flow);

    for (auto& level : levels)
    {
        const std::int64_t level_size = static_cast<std::int64_t>(level.size());

#ifdef _OPENMP
        #pragma omp parallel for schedule(static)
#endif
        for (std::int64_t li = 0; li < level_size; ++li)
        {
            const std::int64_t flat = level[static_cast<std::size_t>(li)];

            if (!flow.is_valid_flat(flat))
                continue;

            const bool outlet  = flow.is_outlet_flat(flat);
            const bool channel = is_channel[static_cast<std::size_t>(flat)];
            const bool is_sink = outlet || channel;

            const float flux_in = flux[static_cast<std::size_t>(flat)];

            if (is_sink)
            {
                // Sink: all incoming flux deposits (unless channel + no-deposit mode)
                if (!channel || deposit_on_hillslope_baselevels)
                    z[static_cast<std::size_t>(flat)] += flux_in;
                // No outgoing flux
                continue;
            }

            // Slope from D-inf router
            float S = slopes[static_cast<std::size_t>(flat)];
            if (S < 0.0f)
                S = 0.0f;

            // Compute: deposit (stays at node), erosion (leaves substrate),
            // outgoing (sent downstream)
            //
            // Material balance: flux_in + erosion = deposit + outgoing
            // Net elevation change: dz = deposit - erosion
            float deposit  = 0.0f;
            float erosion  = 0.0f;
            float outgoing = 0.0f;

            if (!nonlinear)
            {
                // Linear Cidre (L = 1):
                //   all incoming flux deposits locally
                //   only freshly eroded material moves on
                erosion  = kappa_dt_over_dx * S;
                deposit  = flux_in;
                outgoing = erosion;
            }
            else
            {
                if (S < S_c)
                {
                    const float S2 = S * S;
                    erosion  = kappa_dt_over_dx * S;
                    deposit  = flux_in * (1.0f - S2 / Sc2);
                    outgoing = erosion + flux_in * S2 / Sc2;
                }
                else
                {
                    // Slope-restoration (Diffusion_RetourPenteCritique):
                    //   removes exactly the excess above Sc as elevation
                    erosion  = dx * (S - S_c);
                    deposit  = 0.0f;
                    outgoing = erosion + flux_in;
                }
            }

            // Update elevation: node gains deposit, loses erosion
            z[static_cast<std::size_t>(flat)] += deposit - erosion;

            // Send outgoing flux to receivers via atomic accumulation
            const std::int64_t nr = flow.n_receivers(flat);
            for (std::int64_t k = 0; k < nr; ++k)
            {
                const std::int64_t recv = flow.receiver_index_at(flat, k);
                const float        w    = flow.receiver_weight_at(flat, k);
                if (recv >= 0 && recv < n_cells)
                {
                    const float contrib = outgoing * w;
#ifdef _OPENMP
                    #pragma omp atomic
#endif
                    flux[static_cast<std::size_t>(recv)] += contrib;
                }
            }
        }
    }

    return lsdtt_raster_cpu(
        std::move(z),
        std::vector<uint8_t>(dem.get_mask()),
        meta
    );
}


// -----------------------------------------------------------------------------
// compute_roering_hillslope_transport
// -----------------------------------------------------------------------------
lsdtt_raster_cpu compute_roering_hillslope_transport(
    const lsdtt_raster_cpu&         dem,
    const std::vector<uint8_t>&     is_channel,
    float kappa,
    float dt_total,
    bool  nonlinear,
    float S_c,
    int   n_threads)
{
#ifdef _OPENMP
    if (n_threads > 0)
        omp_set_num_threads(n_threads);
#else
    (void)n_threads;
#endif

    const lsdtt_raster_meta& meta  = dem.get_meta();
    const int    nrows  = meta.n_rows;
    const int    ncols  = meta.n_cols;
    const float  h      = static_cast<float>(meta.cell_size);
    const float  h2     = h * h;
    const float  Sc2    = S_c * S_c;
    const std::int64_t n_cells = static_cast<std::int64_t>(nrows) * ncols;

    const std::vector<uint8_t>& mask = dem.get_mask();

    // Working elevation array — updated in-place across sub-steps.
    std::vector<float> z(dem.get_data());

    // Read neighbour elevation with Neumann zero-flux at domain edges and
    // at masked interior cells (reflect the donor cell back).
    auto nbr_z = [&](int r, int c, int dr, int dc) -> float
    {
        int nr = std::max(0, std::min(nrows - 1, r + dr));
        int nc = std::max(0, std::min(ncols - 1, c + dc));
        const std::size_t nidx = static_cast<std::size_t>(nr * ncols + nc);
        if (!mask[nidx])
        {
            // masked neighbour: reflect back to (r,c)
            return z[static_cast<std::size_t>(r * ncols + c)];
        }
        return z[nidx];
    };

    float dt_remaining = dt_total;

    while (dt_remaining > 0.0f)
    {
        // -----------------------------------------------------------------
        // Determine stable sub-step size.
        // For the nonlinear case scan the current field for α_max.
        // For the linear case α = 1 and dt_stable is fixed.
        // -----------------------------------------------------------------
        float dt_sub;
        if (!nonlinear)
        {
            // h²/(4D) — always larger than any realistic LEM dt
            dt_sub = std::min(dt_remaining, h2 / (4.0f * kappa));
        }
        else
        {
            float alpha_max = 1.0f;
#ifdef _OPENMP
            #pragma omp parallel for schedule(static) reduction(max : alpha_max)
#endif
            for (std::int64_t idx = 0; idx < n_cells; ++idx)
            {
                const std::size_t sidx = static_cast<std::size_t>(idx);
                if (!mask[sidx] || is_channel[sidx]) continue;

                const int row = static_cast<int>(idx / ncols);
                const int col = static_cast<int>(idx % ncols);
                const float zc = z[sidx];

                const float zN = nbr_z(row, col, -1,  0);
                const float zS = nbr_z(row, col, +1,  0);
                const float zE = nbr_z(row, col,  0, +1);
                const float zW = nbr_z(row, col,  0, -1);

                // Full gradient magnitude squared / Sc² at this cell centre
                const float gx = (zE - zW) / (2.0f * h);
                const float gy = (zN - zS) / (2.0f * h);
                const float s2 = std::min((gx*gx + gy*gy) / Sc2, 0.9999f);
                const float alpha = 1.0f / (1.0f - s2);
                if (alpha > alpha_max) alpha_max = alpha;

                (void)zc;
            }
            // dt_stable = h²/(4·D·α_max)
            const float dt_stable = h2 / (4.0f * kappa * alpha_max);
            dt_sub = std::min(dt_remaining, dt_stable);
        }

        // Safety floor to prevent infinite loops
        if (dt_sub <= 0.0f) dt_sub = dt_remaining;

        const float kdt_over_h2 = kappa * dt_sub / h2;

        // -----------------------------------------------------------------
        // Explicit FD update.
        // For the face alphas, use the same cross-slope correction as the
        // steady-state solver (see lsdtt_hillslope_steady_state_cpu.cpp).
        //   N/S faces: cross-slope = (zE-zW)/(2h)
        //   E/W faces: cross-slope = (zN-zS)/(2h)
        // -----------------------------------------------------------------
        std::vector<float> z_new(z);  // copy: update into z_new, read from z

#ifdef _OPENMP
        #pragma omp parallel for schedule(static)
#endif
        for (std::int64_t idx = 0; idx < n_cells; ++idx)
        {
            const std::size_t sidx = static_cast<std::size_t>(idx);
            if (!mask[sidx] || is_channel[sidx]) continue;

            const int row = static_cast<int>(idx / ncols);
            const int col = static_cast<int>(idx % ncols);
            const float zc = z[sidx];

            const float zN = nbr_z(row, col, -1,  0);
            const float zS = nbr_z(row, col, +1,  0);
            const float zE = nbr_z(row, col,  0, +1);
            const float zW = nbr_z(row, col,  0, -1);

            float aN = 1.0f, aS = 1.0f, aE = 1.0f, aW = 1.0f;
            if (nonlinear)
            {
                // N/S faces: face-normal slope + E-W cross-slope
                const float sx = (zE - zW) / (2.0f * h);
                const float cross_ns = sx * sx / Sc2;

                auto alpha_face = [&](float z_nbr, float cross2) -> float
                {
                    const float dz_h = (z_nbr - zc) / h;
                    const float s2 = std::min(dz_h * dz_h / Sc2 + cross2, 0.9999f);
                    return 1.0f / (1.0f - s2);
                };

                // E/W faces: face-normal slope + N-S cross-slope
                const float sy = (zN - zS) / (2.0f * h);
                const float cross_ew = sy * sy / Sc2;

                aN = alpha_face(zN, cross_ns);
                aS = alpha_face(zS, cross_ns);
                aE = alpha_face(zE, cross_ew);
                aW = alpha_face(zW, cross_ew);
            }

            // dz/dt = D/h² * [αN*(zN-zc) + αS*(zS-zc) + αE*(zE-zc) + αW*(zW-zc)]
            const float laplacian = aN * (zN - zc) + aS * (zS - zc)
                                  + aE * (zE - zc) + aW * (zW - zc);
            z_new[sidx] = zc + kdt_over_h2 * laplacian;
        }

        z.swap(z_new);
        dt_remaining -= dt_sub;
    }

    return lsdtt_raster_cpu(
        std::move(z),
        std::vector<uint8_t>(mask),
        meta
    );
}

// -----------------------------------------------------------------------------
// compute_adi_hillslope_transport
// -----------------------------------------------------------------------------
lsdtt_raster_cpu compute_adi_hillslope_transport(
    const lsdtt_raster_cpu& dem,
    float kappa,
    float dt,
    const std::string& boundary_conditions)
{
    const lsdtt_raster_meta& meta = dem.get_meta();
    const int nrows = meta.n_rows;
    const int ncols = meta.n_cols;
    const float h   = static_cast<float>(meta.cell_size);

    const std::vector<uint8_t>& mask = dem.get_mask();

    std::vector<float> z_old = dem.get_data();
    std::vector<float> z_half(z_old.size(), 0.0f);
    std::vector<float> z_new(z_old.size(), 0.0f);

    const bool per_n = (boundary_conditions.size() > 0 && boundary_conditions[0] == 'p');
    const bool per_e = (boundary_conditions.size() > 1 && boundary_conditions[1] == 'p');
    const bool per_s = (boundary_conditions.size() > 2 && boundary_conditions[2] == 'p');
    const bool per_w = (boundary_conditions.size() > 3 && boundary_conditions[3] == 'p');

    // Safe neighbour Z fetch using reflection logic for boundary edges and masked cells
    auto nbr_z = [&](int r, int c, int dr, int dc, const std::vector<float>& arr) -> float
    {
        int nr = r + dr;
        int nc = c + dc;

        if (nr < 0) { if (per_n) nr = nrows - 1; else nr = std::max(0, r - dr); }
        else if (nr >= nrows) { if (per_s) nr = 0; else nr = std::min(nrows - 1, r - dr); }

        if (nc < 0) { if (per_w) nc = ncols - 1; else nc = std::max(0, c - dc); }
        else if (nc >= ncols) { if (per_e) nc = 0; else nc = std::min(ncols - 1, c - dc); }

        std::size_t nidx = static_cast<std::size_t>(nr * ncols + nc);
        if (!mask[nidx]) 
        {
            int ir = std::max(0, std::min(nrows - 1, r - dr));
            int ic = std::max(0, std::min(ncols - 1, c - dc));
            return arr[static_cast<std::size_t>(ir * ncols + ic)];
        }
        return arr[nidx];
    };

    auto is_valid = [&](int r, int c) -> bool {
        if (r < 0 || r >= nrows || c < 0 || c >= ncols) return false;
        return mask[static_cast<std::size_t>(r * ncols + c)] > 0;
    };

    float base_gamma = kappa * dt / (2.0f * h * h);
    int max_dim = std::max(nrows, ncols);
    
    std::vector<float> A(max_dim), B(max_dim), C_arr(max_dim), D(max_dim);
    std::vector<float> c_prime(max_dim), d_prime(max_dim), x_res(max_dim);

    // -----------------------------------------------------------------
    // Sweep 1: Implicit X-direction, Explicit Y-direction
    // -----------------------------------------------------------------
    for (int r = 0; r < nrows; ++r)
    {
        for (int c = 0; c < ncols; ++c)
        {
            std::size_t idx = static_cast<std::size_t>(r * ncols + c);
            if (!mask[idx])
            {
                A[c] = 0.0f; B[c] = 1.0f; C_arr[c] = 0.0f; D[c] = z_old[idx];
            }
            else
            {
                float gN = base_gamma, gS = base_gamma;
                float gE = base_gamma, gW = base_gamma;

                // Reduce flux to zero at open boundaries and mask borders
                if (!per_n && !is_valid(r - 1, c)) gN = 0.0f;
                if (!per_s && !is_valid(r + 1, c)) gS = 0.0f;
                if (!per_e && !is_valid(r, c + 1)) gE = 0.0f;
                if (!per_w && !is_valid(r, c - 1)) gW = 0.0f;

                // For periodic bounds, push wrapped element to explicit terms to keep tridiagonal matrix intact
                if (per_w && c == 0) gW = 0.0f;
                if (per_e && c == ncols - 1) gE = 0.0f;

                A[c] = -gW;
                C_arr[c] = -gE;
                B[c] = 1.0f + gW + gE;

                float zC = z_old[idx];
                float zN = nbr_z(r, c, -1, 0, z_old);
                float zS = nbr_z(r, c, +1, 0, z_old);
                
                D[c] = zC + gN * (zN - zC) + gS * (zS - zC);
                
                // Explictly inject wrap connections
                if (per_w && c == 0) D[c] += base_gamma * (nbr_z(r, c, 0, -1, z_old) - zC);
                if (per_e && c == ncols - 1) D[c] += base_gamma * (nbr_z(r, c, 0, +1, z_old) - zC);
            }
        }

        // Thomas Algorithm (TDMA)
        c_prime[0] = C_arr[0] / B[0];
        d_prime[0] = D[0] / B[0];
        for (int i = 1; i < ncols; ++i)
        {
            float m = 1.0f / (B[i] - A[i] * c_prime[i-1]);
            c_prime[i] = C_arr[i] * m;
            d_prime[i] = (D[i] - A[i] * d_prime[i-1]) * m;
        }
        x_res[ncols - 1] = d_prime[ncols - 1];
        for (int i = ncols - 2; i >= 0; --i)
        {
            x_res[i] = d_prime[i] - c_prime[i] * x_res[i + 1];
        }

        for (int c = 0; c < ncols; ++c) {
            z_half[static_cast<std::size_t>(r * ncols + c)] = x_res[c];
        }
    }

    // -----------------------------------------------------------------
    // Sweep 2: Implicit Y-direction, Explicit X-direction
    // -----------------------------------------------------------------
    for (int c = 0; c < ncols; ++c)
    {
        for (int r = 0; r < nrows; ++r)
        {
            std::size_t idx = static_cast<std::size_t>(r * ncols + c);
            if (!mask[idx])
            {
                A[r] = 0.0f; B[r] = 1.0f; C_arr[r] = 0.0f; D[r] = z_half[idx];
            }
            else
            {
                float gN = base_gamma, gS = base_gamma;
                float gE = base_gamma, gW = base_gamma;

                // Reduce flux to zero at open boundaries and mask borders
                if (!per_n && !is_valid(r - 1, c)) gN = 0.0f;
                if (!per_s && !is_valid(r + 1, c)) gS = 0.0f;
                if (!per_e && !is_valid(r, c + 1)) gE = 0.0f;
                if (!per_w && !is_valid(r, c - 1)) gW = 0.0f;

                // For periodic bounds, push wrapped element to explicit terms
                if (per_n && r == 0) gN = 0.0f;
                if (per_s && r == nrows - 1) gS = 0.0f;

                A[r] = -gN;
                C_arr[r] = -gS;
                B[r] = 1.0f + gN + gS;

                float zC = z_half[idx];
                float zE = nbr_z(r, c, 0, +1, z_half);
                float zW = nbr_z(r, c, 0, -1, z_half);
                
                D[r] = zC + gE * (zE - zC) + gW * (zW - zC);

                // Explictly inject wrap connections
                if (per_n && r == 0) D[r] += base_gamma * (nbr_z(r, c, -1, 0, z_half) - zC);
                if (per_s && r == nrows - 1) D[r] += base_gamma * (nbr_z(r, c, +1, 0, z_half) - zC);
            }
        }

        // Thomas Algorithm (TDMA)
        c_prime[0] = C_arr[0] / B[0];
        d_prime[0] = D[0] / B[0];
        for (int i = 1; i < nrows; ++i)
        {
            float m = 1.0f / (B[i] - A[i] * c_prime[i-1]);
            c_prime[i] = C_arr[i] * m;
            d_prime[i] = (D[i] - A[i] * d_prime[i-1]) * m;
        }
        x_res[nrows - 1] = d_prime[nrows - 1];
        for (int i = nrows - 2; i >= 0; --i)
        {
            x_res[i] = d_prime[i] - c_prime[i] * x_res[i + 1];
        }

        for (int r = 0; r < nrows; ++r) {
            z_new[static_cast<std::size_t>(r * ncols + c)] = x_res[r];
        }
    }

    return lsdtt_raster_cpu(
        std::move(z_new),
        std::vector<uint8_t>(mask),
        meta
    );
}


} // namespace lsdtt
