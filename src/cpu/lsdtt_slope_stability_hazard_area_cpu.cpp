// lsdtt_slope_stability_hazard_area_cpu.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html

#include "lsdtt_slope_stability_hazard_area_cpu.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <limits>
#include <stdexcept>
#include <utility>

#include "lsdtt_flowgraph_cpu.hpp"

namespace lsdtt
{

namespace
{

constexpr float k_output_nodata = -9999.0f;

// Steepest descent over a cell's graph receivers (wrapping on periodic edges).
// A cell with no receivers, such as an outlet beside nodata, uses its in-grid
// valid neighbours, so that it can still be a source or a stopping cell.
std::vector<double> steepest_gradient(const lsdtt_flowgraph_cpu& graph, const lsdtt_raster_cpu& surface)
{
    const lsdtt_raster_meta& meta = surface.get_meta();
    const int nr = meta.n_rows, nc = meta.n_cols;
    const std::size_t n = static_cast<std::size_t>(nr) * static_cast<std::size_t>(nc);
    const double dx = meta.cell_size, dd = dx * std::sqrt(2.0);
    const std::vector<float>& z = surface.get_data();
    const std::vector<std::uint8_t>& valid = graph.valid_cells();
    const std::vector<std::int64_t>& offsets = graph.receiver_offsets();
    const std::vector<std::int64_t>& receivers = graph.receiver_indices();

    std::vector<double> gradient(n, 0.0);
    for (int row = 0; row < nr; ++row)
        for (int col = 0; col < nc; ++col)
        {
            const std::size_t i = static_cast<std::size_t>(row) * nc + col;
            if (!valid[i]) continue;
            double best = 0.0;
            if (offsets[i + 1] > offsets[i])
            {
                for (std::int64_t k = offsets[i]; k < offsets[i + 1]; ++k)
                {
                    const std::int64_t r = receivers[k];
                    int dr = std::abs(static_cast<int>(r / nc) - row);
                    int dc = std::abs(static_cast<int>(r % nc) - col);
                    dr = std::min(dr, nr - dr);
                    dc = std::min(dc, nc - dc);
                    const double length = (dr != 0 && dc != 0) ? dd : dx;
                    best = std::max(best, (static_cast<double>(z[i]) - z[r]) / length);
                }
            }
            else
            {
                for (int a = -1; a <= 1; ++a)
                    for (int b = -1; b <= 1; ++b)
                    {
                        if (a == 0 && b == 0) continue;
                        const int u = row + a, v = col + b;
                        if (u < 0 || u >= nr || v < 0 || v >= nc) continue;
                        const std::size_t j = static_cast<std::size_t>(u) * nc + v;
                        if (!valid[j]) continue;
                        best = std::max(best, (static_cast<double>(z[i]) - z[j]) / ((a != 0 && b != 0) ? dd : dx));
                    }
            }
            gradient[i] = best;
        }
    return gradient;
}

lsdtt_raster_cpu output_raster(std::vector<float> values, const std::vector<std::uint8_t>& valid,
                               const lsdtt_raster_meta& meta)
{
    // The routing surface's nodata value may be 0, a legitimate Ah.
    lsdtt_raster_meta out = meta;
    out.no_data_value = k_output_nodata;
    for (std::size_t i = 0; i < values.size(); ++i)
        if (!valid[i]) values[i] = k_output_nodata;
    return lsdtt_raster_cpu(std::move(values), std::vector<std::uint8_t>(valid), out);
}

}  // namespace

void validate_coseismic_hazard_area_options(const coseismic_hazard_area_options& o)
{
    if (!std::isfinite(o.initiation_slope) || !std::isfinite(o.stopping_slope) || !(o.stopping_slope > 0.0) ||
        !(o.initiation_slope > o.stopping_slope))
        throw std::invalid_argument(
            "coseismic hazard area: need 0 < stopping slope < initiation slope (gradients, m/m)");
}

std::vector<double> route_hazard_area(const std::vector<std::int64_t>& order,
                                      const std::vector<std::int64_t>& receiver_offsets,
                                      const std::vector<std::int64_t>& receiver_indices,
                                      const std::vector<float>& receiver_weights,
                                      const std::vector<std::uint8_t>& passes,
                                      const std::vector<double>& seed,
                                      std::vector<std::uint8_t>* reached)
{
    const std::size_t n = seed.size();
    if (passes.size() != n || receiver_offsets.size() != n + 1 ||
        receiver_indices.size() != receiver_weights.size())
        throw std::invalid_argument("route_hazard_area: inconsistent graph sizes");

    std::vector<double> total = seed;
    std::vector<std::uint8_t> hit(n, 0);
    for (std::size_t i = 0; i < n; ++i) hit[i] = seed[i] > 0.0 ? 1 : 0;

    for (const std::int64_t flat : order)
    {
        const std::size_t i = static_cast<std::size_t>(flat);
        // `hit`, not total > 0: a long chain of small fractions can underflow.
        if (!passes[i] || !hit[i]) continue;
        for (std::int64_t k = receiver_offsets[i]; k < receiver_offsets[i + 1]; ++k)
        {
            if (!(receiver_weights[k] > 0.0f)) continue;
            const std::size_t r = static_cast<std::size_t>(receiver_indices[k]);
            total[r] += total[i] * static_cast<double>(receiver_weights[k]);
            hit[r] = 1;
        }
    }
    if (reached) *reached = std::move(hit);
    return total;
}

coseismic_hazard_area_result compute_coseismic_hazard_area(const lsdtt_raster_cpu& routing_surface,
                                                           const coseismic_hazard_area_options& options)
{
    validate_coseismic_hazard_area_options(options);
    const lsdtt_raster_meta& meta = routing_surface.get_meta();
    const int nr = meta.n_rows, nc = meta.n_cols;
    const std::size_t n = static_cast<std::size_t>(nr) * static_cast<std::size_t>(nc);

    lsdtt_flowgraph_options flow;
    flow.method = lsdtt_flowgraph_method::mfd;
    flow.mfd_exponent = 1.0f;
    flow.boundary_conditions = options.boundary_conditions;
    flow.n_threads = options.n_threads;
    const lsdtt_flowgraph_cpu graph(routing_surface.get_data_array(), routing_surface.get_mask_array(), meta, flow);
    const std::vector<std::uint8_t>& valid = graph.valid_cells();
    const std::vector<double> gradient = steepest_gradient(graph, routing_surface);

    // Thresholds are compared in Float32 so that a CLI value (parsed as float)
    // and the same value from Python select exactly the same cells.
    const float initiation = static_cast<float>(options.initiation_slope);
    const float stopping = static_cast<float>(options.stopping_slope);

    coseismic_hazard_area_result out;
    std::vector<double> seed(n, 0.0), edge_seed(n, 0.0);
    std::vector<std::uint8_t> passes(n, 0);
    std::vector<float> zone(n, 0.0f), slope(n, 0.0f);
    for (int row = 0; row < nr; ++row)
        for (int col = 0; col < nc; ++col)
        {
            const std::size_t i = static_cast<std::size_t>(row) * nc + col;
            if (!valid[i]) continue;
            ++out.n_valid;
            const float g = static_cast<float>(gradient[i]);
            slope[i] = g;
            passes[i] = g >= stopping ? 1 : 0;
            if (g >= initiation)
            {
                seed[i] = meta.cell_size;  // a_i / l_j
                zone[i] = k_hazard_zone_source;
                ++out.n_source;
            }
            else if (g < stopping)
            {
                zone[i] = k_hazard_zone_stopping;
                ++out.n_stopping;
            }
            bool edge = row == 0 || col == 0 || row == nr - 1 || col == nc - 1;
            for (int a = -1; a <= 1 && !edge; ++a)
                for (int b = -1; b <= 1 && !edge; ++b)
                    edge = !valid[static_cast<std::size_t>(row + a) * nc + (col + b)];
            if (edge) edge_seed[i] = 1.0;
        }
    out.source_area_m2 = static_cast<double>(out.n_source) * meta.cell_size * meta.cell_size;

    std::vector<std::uint8_t> reached, edge_reached;
    const std::vector<double> total = route_hazard_area(graph.topological_order(), graph.receiver_offsets(),
                                                        graph.receiver_indices(), graph.receiver_weights(), passes,
                                                        seed, &reached);
    (void)route_hazard_area(graph.topological_order(), graph.receiver_offsets(), graph.receiver_indices(),
                            graph.receiver_weights(), passes, edge_seed, &edge_reached);

    std::vector<float> hazard(n, 0.0f), edge(n, 0.0f);
    for (std::size_t i = 0; i < n; ++i)
    {
        if (!valid[i]) continue;
        if (reached[i])
            hazard[i] = std::max(static_cast<float>(total[i]), std::numeric_limits<float>::min());
        else
            ++out.n_zero_hazard;
        if (edge_reached[i])
        {
            edge[i] = 1.0f;
            ++out.n_edge_influenced;
        }
    }

    out.hazard_area = output_raster(std::move(hazard), valid, meta);
    out.steepest_slope = output_raster(std::move(slope), valid, meta);
    out.zone = output_raster(std::move(zone), valid, meta);
    out.edge_influenced = output_raster(std::move(edge), valid, meta);
    return out;
}

}  // namespace lsdtt
