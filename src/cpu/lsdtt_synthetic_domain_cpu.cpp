// lsdtt_synthetic_domain_cpu.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
//
// See lsdtt_synthetic_domain_cpu.hpp for the description.
// -----------------------------------------------------------------------------

#include <cstdint>
#include <stdexcept>
#include "lsdtt_synthetic_domain_cpu.hpp"
#include "../utils/lsdtt_connected_components.hpp"
#include "../utils/lsdtt_grid_neighbours.hpp"

#include <algorithm>
#include <cmath>
#include <queue>
#include <random>
#include <vector>

#include "lsdtt_raster_synthetic_surfaces_cpu.hpp"
#include "../utils/lsdtt_polygon_geometry.hpp"

namespace lsdtt {

// -----------------------------------------------------------------------------
void keep_largest_land_component(
    std::vector<std::uint8_t>& mask, int n_rows, int n_cols)
{
    const std::size_t n = static_cast<std::size_t>(n_rows) * n_cols;

    // zero the one-pixel border so the coastline is internal
    for (int c = 0; c < n_cols; ++c)
    {
        mask[static_cast<std::size_t>(0) * n_cols + c] = 0;
        mask[static_cast<std::size_t>(n_rows - 1) * n_cols + c] = 0;
    }
    for (int r = 0; r < n_rows; ++r)
    {
        mask[static_cast<std::size_t>(r) * n_cols + 0] = 0;
        mask[static_cast<std::size_t>(r) * n_cols + (n_cols - 1)] = 0;
    }

    // The largest 4-connected land component; on a tie, the first in row-major
    // order. With no land at all best_label stays 0 and every cell is kept,
    // as before.
    std::vector<std::int32_t> label;
    std::vector<std::int64_t> size;
    const std::int32_t n_comp = label_connected_components(n_rows, n_cols, 4,
        [&](std::size_t s) { return mask[s] != 0; }, label, size);
    std::int32_t best_label = 0;
    std::int64_t best_size = 0;
    for (std::int32_t l = 1; l <= n_comp; ++l)
    {
        if (size[static_cast<std::size_t>(l)] > best_size)
        {
            best_size = size[static_cast<std::size_t>(l)];
            best_label = l;
        }
    }

    for (std::size_t s = 0; s < n; ++s)
        mask[s] = (label[s] == best_label) ? 1 : 0;
}

// -----------------------------------------------------------------------------
std::vector<std::uint8_t> make_efd_island_mask(
    int n_rows, int n_cols,
    int n_harmonics, double roughness, double fill_fraction,
    std::uint64_t seed)
{
    std::mt19937_64 rng(seed);
    std::normal_distribution<double> g(0.0, 1.0);
    std::uniform_real_distribution<double> u(0.0, 1.0);

    const int N = std::max(1, n_harmonics);
    std::vector<double> a(N + 1, 0.0), b(N + 1, 0.0),
                        c(N + 1, 0.0), d(N + 1, 0.0);

    // base ellipse (k = 1) with a random aspect ratio and orientation
    const double ecc = 0.25 * (u(rng) - 0.5);   // mild eccentricity
    a[1] = 1.0 + ecc;
    d[1] = 1.0 - ecc;
    b[1] = 0.15 * g(rng);
    c[1] = 0.15 * g(rng);

    // higher harmonics: amplitude decays as k^-roughness
    for (int k = 2; k <= N; ++k)
    {
        const double amp = 0.35 * std::pow(static_cast<double>(k), -roughness);
        a[k] = amp * g(rng);
        b[k] = amp * g(rng);
        c[k] = amp * g(rng);
        d[k] = amp * g(rng);
    }

    // sample the closed contour
    const int M = 1024;
    std::vector<lsdtt_xy> ring(M);
    double min_x = 1e30, max_x = -1e30, min_y = 1e30, max_y = -1e30;
    for (int m = 0; m < M; ++m)
    {
        const double t = 2.0 * M_PI * static_cast<double>(m) / M;
        double x = 0.0, y = 0.0;
        for (int k = 1; k <= N; ++k)
        {
            const double ck = std::cos(k * t), sk = std::sin(k * t);
            x += a[k] * ck + b[k] * sk;
            y += c[k] * ck + d[k] * sk;
        }
        ring[m] = lsdtt_xy{x, y};
        min_x = std::min(min_x, x); max_x = std::max(max_x, x);
        min_y = std::min(min_y, y); max_y = std::max(max_y, y);
    }

    // scale/translate so the contour spans fill_fraction of the smaller grid
    // dimension, centred, with a one-pixel margin guaranteed by the component step
    const double span = std::max(max_x - min_x, max_y - min_y);
    const double frac = std::min(0.95, std::max(0.1, fill_fraction));
    const double target = frac * std::min(n_rows, n_cols);
    const double scale = (span > 0.0) ? target / span : 1.0;
    const double cx = 0.5 * (min_x + max_x);
    const double cy = 0.5 * (min_y + max_y);
    const double gcx = 0.5 * n_cols;
    const double gcy = 0.5 * n_rows;
    for (auto& p : ring)
    {
        p.x = gcx + (p.x - cx) * scale;
        p.y = gcy + (p.y - cy) * scale;
    }

    std::vector<std::uint8_t> mask(
        static_cast<std::size_t>(n_rows) * n_cols, 0);
    for (int r = 0; r < n_rows; ++r)
        for (int cc = 0; cc < n_cols; ++cc)
            if (point_in_ring(lsdtt_xy{cc + 0.5, r + 0.5}, ring))
                mask[static_cast<std::size_t>(r) * n_cols + cc] = 1;

    keep_largest_land_component(mask, n_rows, n_cols);
    return mask;
}

// -----------------------------------------------------------------------------
std::vector<std::uint8_t> make_noise_island_mask(
    int n_rows, int n_cols,
    int blur_passes, double fill_fraction,
    std::uint64_t seed)
{
    const lsdtt_raster_cpu noise = generate_red_noise_surface(
        n_rows, n_cols, 1.0, static_cast<int>(seed), 0.0, 1.0,
        std::max(1, blur_passes));
    const std::vector<float>& nd = noise.get_data();

    // radial cone taper: high in the centre, falling to < 0 at the edge
    const double gcx = 0.5 * n_cols;
    const double gcy = 0.5 * n_rows;
    const double rmax = 0.5 * std::min(n_rows, n_cols);

    std::vector<double> field(static_cast<std::size_t>(n_rows) * n_cols);
    for (int r = 0; r < n_rows; ++r)
    {
        for (int c = 0; c < n_cols; ++c)
        {
            const std::size_t s = static_cast<std::size_t>(r) * n_cols + c;
            const double dx = (c + 0.5) - gcx;
            const double dy = (r + 0.5) - gcy;
            const double rr = std::sqrt(dx * dx + dy * dy) / rmax;  // 0 centre, ~1 edge
            const double taper = 1.0 - rr;          // cone
            field[s] = static_cast<double>(nd[s]) + 2.0 * taper;
        }
    }

    // threshold so that ~fill_fraction of pixels are land
    std::vector<double> sorted = field;
    std::sort(sorted.begin(), sorted.end());
    const double frac = std::min(0.9, std::max(0.05, fill_fraction));
    const std::size_t idx = static_cast<std::size_t>(
        (1.0 - frac) * static_cast<double>(sorted.size() - 1));
    const double level = sorted[idx];

    std::vector<std::uint8_t> mask(field.size(), 0);
    for (std::size_t s = 0; s < field.size(); ++s)
        mask[s] = (field[s] > level) ? 1 : 0;

    keep_largest_land_component(mask, n_rows, n_cols);
    return mask;
}

// -----------------------------------------------------------------------------
std::vector<std::uint8_t> rasterise_polygon_scaled(
    const std::vector<lsdtt_xy>& ring_in, int n_rows, int n_cols)
{
    std::vector<std::uint8_t> mask(
        static_cast<std::size_t>(n_rows) * n_cols, 0);
    if (ring_in.size() < 3) return mask;

    double min_x = 1e30, max_x = -1e30, min_y = 1e30, max_y = -1e30;
    for (const auto& p : ring_in)
    {
        min_x = std::min(min_x, p.x); max_x = std::max(max_x, p.x);
        min_y = std::min(min_y, p.y); max_y = std::max(max_y, p.y);
    }

    // scale the polygon bbox to fill 90% of the grid, centred, preserving aspect
    const double w = max_x - min_x, h = max_y - min_y;
    const double scale = 0.90 * std::min(
        (w > 0.0 ? n_cols / w : 1.0),
        (h > 0.0 ? n_rows / h : 1.0));
    const double pcx = 0.5 * (min_x + max_x);
    const double pcy = 0.5 * (min_y + max_y);
    const double gcx = 0.5 * n_cols;
    const double gcy = 0.5 * n_rows;

    std::vector<lsdtt_xy> ring(ring_in.size());
    for (std::size_t i = 0; i < ring_in.size(); ++i)
    {
        ring[i].x = gcx + (ring_in[i].x - pcx) * scale;
        // flip y so polygon north maps to small row indices (top of grid)
        ring[i].y = gcy - (ring_in[i].y - pcy) * scale;
    }

    for (int r = 0; r < n_rows; ++r)
        for (int c = 0; c < n_cols; ++c)
            if (point_in_ring(lsdtt_xy{c + 0.5, r + 0.5}, ring))
                mask[static_cast<std::size_t>(r) * n_cols + c] = 1;

    keep_largest_land_component(mask, n_rows, n_cols);
    return mask;
}

// -----------------------------------------------------------------------------
lsdtt_polygon_true_area rasterise_polygon_true_area(
    const std::vector<lsdtt_xy>& ring, double cell_size, int margin_cells)
{
    lsdtt_polygon_true_area out;
    out.cell_size = cell_size;
    if (ring.size() < 3 || !(cell_size > 0.0)) return out;

    double min_x = 1e30, max_x = -1e30, min_y = 1e30, max_y = -1e30;
    for (const auto& p : ring)
    {
        min_x = std::min(min_x, p.x); max_x = std::max(max_x, p.x);
        min_y = std::min(min_y, p.y); max_y = std::max(max_y, p.y);
    }

    const int margin = std::max(1, margin_cells);
    const int span_c = static_cast<int>(std::ceil((max_x - min_x) / cell_size));
    const int span_r = static_cast<int>(std::ceil((max_y - min_y) / cell_size));
    out.n_cols = span_c + 2 * margin;
    out.n_rows = span_r + 2 * margin;
    out.x_min  = min_x - margin * cell_size;   // lower-left corner
    out.y_min  = min_y - margin * cell_size;
    const double y_max = out.y_min + out.n_rows * cell_size;

    out.mask.assign(
        static_cast<std::size_t>(out.n_rows) * out.n_cols, 0);

    // Cell-centre real coordinates: row 0 is the top (largest y).
    for (int r = 0; r < out.n_rows; ++r)
    {
        const double cy = y_max - (r + 0.5) * cell_size;
        for (int c = 0; c < out.n_cols; ++c)
        {
            const double cx = out.x_min + (c + 0.5) * cell_size;
            if (point_in_ring(lsdtt_xy{cx, cy}, ring))
                out.mask[static_cast<std::size_t>(r) * out.n_cols + c] = 1;
        }
    }

    keep_largest_land_component(out.mask, out.n_rows, out.n_cols);
    return out;
}


lsdtt_raster_cpu apply_sinuous_domain_mask(
    const lsdtt_raster_cpu& dem,
    bool  enable_north, float amp_north_fraction, int waves_north, float phase_north,
    bool  enable_south, float amp_south_fraction, int waves_south, float phase_south,
    float baselevel_elevation
)
{
    if (amp_north_fraction < 0.0f || amp_north_fraction > 0.25f ||
        amp_south_fraction < 0.0f || amp_south_fraction > 0.25f)
    {
        throw std::runtime_error(
            "apply_sinuous_domain_mask: amplitude fractions must be in [0, 0.25]");
    }

    const lsdtt_raster_meta& meta = dem.get_meta();
    const int    nr = meta.n_rows;
    const int    nc = meta.n_cols;
    const double dx = meta.cell_size;

    const double Lx = static_cast<double>(nc) * dx;
    const double Ly = static_cast<double>(nr) * dx;

    const double A_n = Ly * static_cast<double>(amp_north_fraction);
    const double A_s = Ly * static_cast<double>(amp_south_fraction);

    std::vector<float>        data = dem.get_data();
    std::vector<std::uint8_t> mask = dem.get_mask();

    constexpr double k_pi = 3.14159265358979323846;

    // A sub-cell tolerance, so a pixel sitting exactly on the analytic boundary
    // is kept rather than being masked by floating-point noise.
    const double edge_tol = 1e-4 * dx;

    for (int r = 0; r < nr; ++r)
    {
        // Row 0 is the NORTH edge, so y counts up from the south.
        const double y = static_cast<double>(nr - 1 - r) * dx;

        for (int c = 0; c < nc; ++c)
        {
            const double x = static_cast<double>(c) * dx;
            bool mask_out = false;

            if (enable_south)
            {
                const double shift_s = A_s * std::sin(
                    2.0 * k_pi * static_cast<double>(waves_south) * (x / Lx) +
                    2.0 * k_pi * static_cast<double>(phase_south));
                const double y_min = A_s + shift_s - edge_tol;
                if (y < y_min) mask_out = true;
            }

            if (enable_north && !mask_out)
            {
                const double shift_n = A_n * std::sin(
                    2.0 * k_pi * static_cast<double>(waves_north) * (x / Lx) +
                    2.0 * k_pi * static_cast<double>(phase_north));
                const double y_max = Ly - A_n + shift_n + edge_tol;
                if (y > y_max) mask_out = true;
            }

            if (mask_out)
            {
                const std::size_t idx =
                    static_cast<std::size_t>(r) * static_cast<std::size_t>(nc) +
                    static_cast<std::size_t>(c);
                // Flattened, not just masked: see the header.
                data[idx] = baselevel_elevation;
                mask[idx] = 0;
            }
        }
    }

    return lsdtt_raster_cpu(std::move(data), std::move(mask), meta);
}

} // namespace lsdtt
