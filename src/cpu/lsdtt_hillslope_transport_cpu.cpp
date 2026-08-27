// lsdtt_hillslope_transport_cpu.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html

#include "lsdtt_hillslope_transport_cpu.hpp"
#include "lsdtt_hillslope_operator_cpu.hpp"

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
    const std::string&              boundary_conditions,
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

    // Working elevation array - updated in-place across sub-steps.
    std::vector<float> z(dem.get_data());

    // Boundary handling comes from the SHARED operator
    // (lsdtt_hillslope_operator_cpu.hpp): 'p' wraps, 'o' is Dirichlet at base
    // level, 'n' is Neumann zero-flux, and a masked interior neighbour reflects.
    //
    // This solver previously had NO boundary-condition concept at all -- it
    // clamped indices unconditionally, so every edge was Neumann, a periodic
    // domain could not be expressed, and an 'o' edge was silently no-flux. That
    // mattered most through the implicit solver's subcycling fallback, which
    // could switch boundary condition mid-run without any record in the output.
    const lsdtt_hillslope_domain dom =
        lsdtt_hillslope_domain::make(nrows, ncols, boundary_conditions, mask);

    auto is_fixed = [&](std::size_t sidx, int row, int col) -> bool
    {
        return !mask[sidx] || is_channel[sidx] || dom.is_dirichlet_edge(row, col);
    };

    float dt_remaining = dt_total;

    while (dt_remaining > 0.0f)
    {
        // -----------------------------------------------------------------
        // Determine stable sub-step size.
        // For the nonlinear case scan the current field for the largest
        // tangent modulus.  For the linear case it is 1 and dt_stable is fixed.
        //
        // The stability limit of an explicit step is set by the TANGENT
        // modulus dq_s/dS, not by the secant (effective) diffusivity
        // D_eff = D/(1-r^2) that appears in the flux law.  Linearising
        // q_s(S) = D*S/(1-r^2), r = S/Sc, about the current slope gives
        //
        //     dq_s/dS = D * (1 + r^2) / (1 - r^2)^2
        //
        // which diverges much faster than D_eff as r -> 1 (a factor
        // (1+r^2)/(1-r^2) larger; ~13x at r = 0.93).  Using D_eff here takes
        // steps far above the true limit on near-critical slopes, which
        // overshoots, drives faces supercritical, and then runs away through
        // the singularity.  Perron (2011) makes the same point.
        // -----------------------------------------------------------------
        float dt_sub;
        if (!nonlinear)
        {
            // h^2/(4D) - always larger than any realistic LEM dt
            dt_sub = std::min(dt_remaining, h2 / (4.0f * kappa));
        }
        else
        {
            float tangent_max = 1.0f;
            const lsdtt_hillslope_stencil<float> cfl_st(dom, z.data(),
                                                        static_cast<double>(h));
#ifdef _OPENMP
            #pragma omp parallel for schedule(static) reduction(max : tangent_max)
#endif
            for (std::int64_t idx = 0; idx < n_cells; ++idx)
            {
                const std::size_t sidx = static_cast<std::size_t>(idx);
                const int row = static_cast<int>(idx / ncols);
                const int col = static_cast<int>(idx % ncols);
                if (is_fixed(sidx, row, col)) continue;

                // Gradient magnitude at the cell centre, wrapping at a periodic
                // seam rather than reflecting -- a clamped read there returns a
                // one-sided (zero) difference and UNDERSTATES the slope, which is
                // the wrong direction to err in a stability estimate.
                const float gx = static_cast<float>(cfl_st.slope_x(row, col));
                const float gy = static_cast<float>(cfl_st.slope_y(row, col));
                const float r2 = std::min((gx*gx + gy*gy) / Sc2, 0.9999f);

                const float one_minus = 1.0f - r2;
                const float tangent = (1.0f + r2) / (one_minus * one_minus);
                if (tangent > tangent_max) tangent_max = tangent;
            }
            // dt_stable = h^2 / (4 * D * max(dq_s/dS)/D)
            const float dt_stable = h2 / (4.0f * kappa * tangent_max);
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

        const lsdtt_hillslope_stencil<float> st(dom, z.data(),
                                                static_cast<double>(h));

#ifdef _OPENMP
        #pragma omp parallel for schedule(static)
#endif
        for (std::int64_t idx = 0; idx < n_cells; ++idx)
        {
            const std::size_t sidx = static_cast<std::size_t>(idx);
            const int row = static_cast<int>(idx / ncols);
            const int col = static_cast<int>(idx % ncols);
            if (is_fixed(sidx, row, col)) continue;

            // Face-symmetric coefficients and face-normal slopes, four at a time.
            // alpha is a property of the FACE, so the tangential slope is averaged
            // across it; taking it at the updating cell (which this solver used to
            // do) gives the two cells sharing a face different diffusivities and
            // makes the scheme non-conservative.
            double a[4], sn[4];
            st.face_alphas(row, col, static_cast<double>(S_c), nonlinear, a, sn);

            // sum_f alpha_f (z_nb - z_c) == h * sum_f alpha_f * S_n,f
            double laplacian = 0.0;
            for (int k = 0; k < 4; ++k) laplacian += a[k] * sn[k] * static_cast<double>(h);

            z_new[sidx] = z[sidx] + kdt_over_h2 * static_cast<float>(laplacian);
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

    // 'o' is a DIRICHLET edge (base level): the edge cell is held at its input
    // elevation, giving an identity row in the tridiagonal solve.  'n' is the
    // Neumann (zero-flux) edge, handled by zeroing the face conductance below.
    // These were previously treated alike, both as zero-flux, which silently
    // turned a declared fixed-elevation boundary into a no-flux one.
    auto bc_at = [&](std::size_t i) -> char
    {
        return (boundary_conditions.size() > i) ? boundary_conditions[i] : 'n';
    };
    const char cN = bc_at(0), cE = bc_at(1), cS = bc_at(2), cW = bc_at(3);

    auto is_o_edge = [&](int r, int c) -> bool
    {
        return (r == 0         && cN == 'o') || (r == nrows - 1 && cS == 'o')
            || (c == 0         && cW == 'o') || (c == ncols - 1 && cE == 'o');
    };

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
    std::vector<float> y_arr(max_dim), q_arr(max_dim);

    // Standard Thomas Algorithm for standard O(N) tridiagonal matrices
    auto solve_thomas = [&](int N, const float* a, const float* b, const float* c, const float* d, float* x, float* cp, float* dp) {
        cp[0] = c[0] / b[0];
        dp[0] = d[0] / b[0];
        for (int i = 1; i < N; ++i) {
            float m = 1.0f / (b[i] - a[i] * cp[i-1]);
            cp[i] = c[i] * m;
            dp[i] = (d[i] - a[i] * dp[i-1]) * m;
        }
        x[N - 1] = dp[N - 1];
        for (int i = N - 2; i >= 0; --i) {
            x[i] = dp[i] - cp[i] * x[i + 1];
        }
    };

    // Sherman-Morrison adaptation of Thomas algorithm for Cyclic Tridiagonal (Periodic) systems
    auto solve_cyclic_thomas = [&](int N, const float* a, const float* b, const float* c, const float* d, float* x, 
                                   float* cp, float* dp, float* y, float* q) {
        if (N <= 2) { solve_thomas(N, a, b, c, d, x, cp, dp); return; }
        
        float alpha = c[N-1];
        float beta  = a[0];
        float gamma = -b[0]; 
        if (std::abs(gamma) < 1e-12f) gamma = 1.0f; // Robust fallback against masked breaks
        
        float b0_mod = b[0] - gamma;
        float bn_mod = b[N-1] - alpha * beta / gamma;
        
        // Tridiagonal solve T*y = d
        cp[0] = c[0] / b0_mod;
        dp[0] = d[0] / b0_mod;
        for (int i = 1; i < N - 1; ++i) {
            float m = 1.0f / (b[i] - a[i] * cp[i-1]);
            cp[i] = c[i] * m;
            dp[i] = (d[i] - a[i] * dp[i-1]) * m;
        }
        {
            int i = N - 1;
            float m = 1.0f / (bn_mod - a[i] * cp[i-1]);
            dp[i] = (d[i] - a[i] * dp[i-1]) * m;
        }
        y[N-1] = dp[N-1];
        for (int i = N - 2; i >= 0; --i) y[i] = dp[i] - cp[i] * y[i+1];
        
        // Tridiagonal solve T*q = u
        dp[0] = gamma / b0_mod;
        for (int i = 1; i < N - 1; ++i) {
            float m = 1.0f / (b[i] - a[i] * cp[i-1]);
            dp[i] = (0.0f - a[i] * dp[i-1]) * m; 
        }
        {
            int i = N - 1;
            float m = 1.0f / (bn_mod - a[i] * cp[i-1]);
            dp[i] = (alpha - a[i] * dp[i-1]) * m;
        }
        q[N-1] = dp[N-1];
        for (int i = N - 2; i >= 0; --i) q[i] = dp[i] - cp[i] * q[i+1];
        
        // Reconstruct the unrolled ring using Sherman-Morrison logic
        float v_dot_y = y[0] + (beta / gamma) * y[N-1];
        float v_dot_q = q[0] + (beta / gamma) * q[N-1];
        float factor = v_dot_y / (1.0f + v_dot_q);
        
        for (int i = 0; i < N; ++i) x[i] = y[i] - factor * q[i];
    };

    // -----------------------------------------------------------------
    // Sweep 1: Implicit X-direction, Explicit Y-direction
    // -----------------------------------------------------------------
    for (int r = 0; r < nrows; ++r)
    {
        for (int c = 0; c < ncols; ++c)
        {
            std::size_t idx = static_cast<std::size_t>(r * ncols + c);
            if (!mask[idx] || is_o_edge(r, c))
            {
                A[c] = 0.0f; B[c] = 1.0f; C_arr[c] = 0.0f; D[c] = z_old[idx];
            }
            else
            {
                float gN = base_gamma, gS = base_gamma;
                float gE = base_gamma, gW = base_gamma;

                if (!per_n && !is_valid(r - 1, c)) gN = 0.0f;
                if (!per_s && !is_valid(r + 1, c)) gS = 0.0f;
                if (!per_e && !is_valid(r, c + 1)) gE = 0.0f;
                if (!per_w && !is_valid(r, c - 1)) gW = 0.0f;

                A[c] = -gW;
                C_arr[c] = -gE;
                B[c] = 1.0f + gW + gE;

                float zC = z_old[idx];
                float zN = nbr_z(r, c, -1, 0, z_old);
                float zS = nbr_z(r, c, +1, 0, z_old);
                
                D[c] = zC + gN * (zN - zC) + gS * (zS - zC);
            }
        }

        if (per_w) solve_cyclic_thomas(ncols, A.data(), B.data(), C_arr.data(), D.data(), x_res.data(), c_prime.data(), d_prime.data(), y_arr.data(), q_arr.data());
        else solve_thomas(ncols, A.data(), B.data(), C_arr.data(), D.data(), x_res.data(), c_prime.data(), d_prime.data());

        for (int c = 0; c < ncols; ++c) z_half[static_cast<std::size_t>(r * ncols + c)] = x_res[c];
    }

    // -----------------------------------------------------------------
    // Sweep 2: Implicit Y-direction, Explicit X-direction
    // -----------------------------------------------------------------
    for (int c = 0; c < ncols; ++c)
    {
        for (int r = 0; r < nrows; ++r)
        {
            std::size_t idx = static_cast<std::size_t>(r * ncols + c);
            if (!mask[idx] || is_o_edge(r, c))
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

                A[r] = -gN;
                C_arr[r] = -gS;
                B[r] = 1.0f + gN + gS;

                float zC = z_half[idx];
                float zE = nbr_z(r, c, 0, +1, z_half);
                float zW = nbr_z(r, c, 0, -1, z_half);
                
                D[r] = zC + gE * (zE - zC) + gW * (zW - zC);
            }
        }

        if (per_n) solve_cyclic_thomas(nrows, A.data(), B.data(), C_arr.data(), D.data(), x_res.data(), c_prime.data(), d_prime.data(), y_arr.data(), q_arr.data());
        else solve_thomas(nrows, A.data(), B.data(), C_arr.data(), D.data(), x_res.data(), c_prime.data(), d_prime.data());

        for (int r = 0; r < nrows; ++r) z_new[static_cast<std::size_t>(r * ncols + c)] = x_res[r];
    }

    return lsdtt_raster_cpu(
        std::move(z_new),
        std::vector<uint8_t>(mask),
        meta
    );
}



} // namespace lsdtt
