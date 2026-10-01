// lsdtt_subbasin_partition_cpu.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// Incremental sub-basin partition and per-sub-basin set-up quantities.
// Definitions and conventions are in the header.
// -----------------------------------------------------------------------------

#include "lsdtt_subbasin_partition_cpu.hpp"

#include "lsdtt_junction_basin_manager_cpu.hpp"
#include "../utils/lsdtt_statistics.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <stdexcept>
#include <string>

namespace lsdtt
{

namespace
{

constexpr double k_nan = std::numeric_limits<double>::quiet_NaN();

void require_same_grid(
    const subbasin_partition_result& partition,
    const lsdtt_raster_cpu&          raster,
    const char*                      who,
    const char*                      name)
{
    const lsdtt_raster_meta& m = raster.get_meta();
    if (m.n_rows != partition.meta.n_rows || m.n_cols != partition.meta.n_cols ||
        raster.get_data().size() != partition.subbasin_id.size() ||
        raster.get_mask().size() != partition.subbasin_id.size())
    {
        throw std::invalid_argument(
            std::string(who) + ": " + name +
            " is not on the partition's grid (n_rows/n_cols differ)");
    }
}

// Pixels of each sub-basin in compressed rows: pixels of s are
// flats[offset[s]] .. flats[offset[s+1]-1], in ascending flat order.
struct subbasin_pixel_rows
{
    std::vector<std::int64_t> offset;
    std::vector<std::int64_t> flats;
};

subbasin_pixel_rows build_pixel_rows(const subbasin_partition_result& partition)
{
    const std::int32_t n_sb = partition.n_subbasins();

    subbasin_pixel_rows rows;
    rows.offset.assign(static_cast<std::size_t>(n_sb) + 1, 0);
    for (std::int32_t s = 0; s < n_sb; ++s)
    {
        rows.offset[static_cast<std::size_t>(s) + 1] =
            rows.offset[static_cast<std::size_t>(s)] +
            partition.n_pixels[static_cast<std::size_t>(s)];
    }

    rows.flats.resize(static_cast<std::size_t>(rows.offset.back()));
    std::vector<std::int64_t> cursor(rows.offset.begin(), rows.offset.end() - 1);

    const std::int64_t n = static_cast<std::int64_t>(partition.subbasin_id.size());
    for (std::int64_t flat = 0; flat < n; ++flat)
    {
        const std::int32_t s = partition.subbasin_id[static_cast<std::size_t>(flat)];
        if (s < 0) continue;
        rows.flats[static_cast<std::size_t>(cursor[static_cast<std::size_t>(s)]++)] = flat;
    }
    return rows;
}

// D8 donors of `flat`: neighbours whose receiver is `flat`. Neighbour indices
// wrap around the grid so that donors across a periodic boundary are found; on
// a non-periodic grid a wrapped cell can never have `flat` as its receiver.
void collect_donors(
    const lsdtt_flow_base&     flow,
    std::int64_t               flat,
    std::vector<std::int64_t>& donors)
{
    static constexpr int d_row[8] = {-1, -1, 0, 1, 1, 1, 0, -1};
    static constexpr int d_col[8] = {0, 1, 1, 1, 0, -1, -1, -1};

    donors.clear();
    const int n_rows = flow.n_rows();
    const int n_cols = flow.n_cols();
    const int row = flow.row_from_flat_unchecked(flat);
    const int col = flow.col_from_flat_unchecked(flat);

    for (int dir = 0; dir < 8; ++dir)
    {
        const int r = ((row + d_row[dir]) % n_rows + n_rows) % n_rows;
        const int c = ((col + d_col[dir]) % n_cols + n_cols) % n_cols;
        const std::int64_t nb = flow.flat_index_unchecked(r, c);

        if (nb == flat) continue;
        if (!flow.is_valid_cell(nb)) continue;
        if (flow.receiver_of_unchecked(nb) != flat) continue;
        if (std::find(donors.begin(), donors.end(), nb) != donors.end()) continue;

        donors.push_back(nb);
    }
}

} // namespace

// =============================================================================
// Outlet selection
// =============================================================================

std::vector<std::int64_t> subbasin_outlets_from_user_flats(
    const lsdtt_flow_base&             flow,
    const lsdtt_junction_network_base& junction_network,
    const std::vector<std::int64_t>&   outlet_flats)
{
    lsdtt_junction_basin_manager_cpu manager(flow, junction_network);
    manager.select_junction_basins_from_basin_outlets(outlet_flats);
    return manager.get_basin_outlet_flat_indices();
}

std::vector<std::int64_t> subbasin_outlets_by_stream_order(
    const lsdtt_flow_base&             flow,
    const lsdtt_junction_network_base& junction_network,
    std::int32_t                       min_order,
    std::int32_t                       max_order)
{
    lsdtt_junction_basin_manager_cpu manager(flow, junction_network);
    manager.select_all_junction_basins();
    manager.filter_by_stream_order(min_order, max_order);
    return manager.get_basin_outlet_flat_indices();
}

// =============================================================================
// Partition
// =============================================================================

subbasin_partition_result partition_incremental_subbasins(
    const lsdtt_flow_base&           flow,
    const std::vector<std::int64_t>& outlet_flats)
{
    static const char* who = "lsdtt::partition_incremental_subbasins";

    if (!flow.has_s_vector())
    {
        throw std::runtime_error(
            std::string(who) + ": the flow object has no s_vector "
            "(set build_s_vector = true)");
    }

    const std::size_t n = flow.full_grid_size();

    subbasin_partition_result out;
    out.meta = flow.meta();
    out.subbasin_id.assign(n, -1);

    // Which sub-basin, if any, starts at each pixel.
    std::vector<std::int32_t> starts_here(n, -1);

    for (const std::int64_t flat : outlet_flats)
    {
        if (!flow.is_flat_in_range(flat) || !flow.is_valid_cell(flat))
        {
            throw std::invalid_argument(
                std::string(who) + ": outlet flat " + std::to_string(flat) +
                " is out of range or not a valid flow cell");
        }
        if (starts_here[static_cast<std::size_t>(flat)] >= 0)
        {
            ++out.n_duplicate_outlets;
            continue;
        }
        starts_here[static_cast<std::size_t>(flat)] =
            static_cast<std::int32_t>(out.outlet_flat.size());
        out.outlet_flat.push_back(flat);
    }

    const std::int32_t n_sb = out.n_subbasins();
    out.n_pixels.assign(static_cast<std::size_t>(n_sb), 0);

    // Sweep from base level upward. The s_vector slice of a base-level outlet
    // is in preorder, so every receiver is labelled before its donors.
    for (std::int64_t base = 0; base < static_cast<std::int64_t>(n); ++base)
    {
        if (!flow.is_valid_cell(base)) continue;
        ++out.n_valid_pixels;
        if (!flow.is_outlet_unchecked(base)) continue;

        for (const std::int64_t node : flow.upslope_nodes_unchecked(base))
        {
            const std::size_t i = static_cast<std::size_t>(node);
            std::int32_t id = starts_here[i];
            if (id < 0 && node != base)
            {
                id = out.subbasin_id[static_cast<std::size_t>(
                    flow.receiver_of_unchecked(node))];
            }
            out.subbasin_id[i] = id;
        }
    }

    std::int64_t n_assigned = 0;
    for (std::size_t i = 0; i < n; ++i)
    {
        const std::int32_t s = out.subbasin_id[i];
        if (s < 0) continue;
        ++out.n_pixels[static_cast<std::size_t>(s)];
        ++n_assigned;
    }
    out.n_unassigned = out.n_valid_pixels - n_assigned;

    // Downstream links (maindown).
    out.downstream.assign(static_cast<std::size_t>(n_sb), -1);
    for (std::int32_t s = 0; s < n_sb; ++s)
    {
        const std::int64_t rec =
            flow.receiver_of_unchecked(out.outlet_flat[static_cast<std::size_t>(s)]);
        if (rec == lsdtt_flow_base::outlet_receiver) continue;
        out.downstream[static_cast<std::size_t>(s)] =
            out.subbasin_id[static_cast<std::size_t>(rec)];
    }

    // Upstream-to-downstream order (Kahn's algorithm on the sub-basin links).
    std::vector<std::int32_t> n_upstream(static_cast<std::size_t>(n_sb), 0);
    for (std::int32_t s = 0; s < n_sb; ++s)
    {
        const std::int32_t d = out.downstream[static_cast<std::size_t>(s)];
        if (d >= 0) ++n_upstream[static_cast<std::size_t>(d)];
    }

    out.upstream_to_downstream.reserve(static_cast<std::size_t>(n_sb));
    for (std::int32_t s = 0; s < n_sb; ++s)
    {
        if (n_upstream[static_cast<std::size_t>(s)] == 0)
            out.upstream_to_downstream.push_back(s);
    }
    for (std::size_t k = 0; k < out.upstream_to_downstream.size(); ++k)
    {
        const std::int32_t d =
            out.downstream[static_cast<std::size_t>(out.upstream_to_downstream[k])];
        if (d >= 0 && --n_upstream[static_cast<std::size_t>(d)] == 0)
            out.upstream_to_downstream.push_back(d);
    }

    return out;
}

lsdtt_raster_result get_subbasin_id_raster(const subbasin_partition_result& partition)
{
    const std::size_t n   = partition.subbasin_id.size();
    const float       ndv = static_cast<float>(partition.meta.no_data_value);

    lsdtt_raster_result r;
    r.meta = partition.meta;
    r.data.assign(n, ndv);
    r.mask.assign(n, 0u);

    for (std::size_t i = 0; i < n; ++i)
    {
        const std::int32_t s = partition.subbasin_id[i];
        if (s < 0) continue;
        r.data[i] = static_cast<float>(s);
        r.mask[i] = 1u;
    }
    return r;
}

// =============================================================================
// GeoData equivalents
// =============================================================================

subbasin_geodata_result compute_subbasin_geodata(
    const subbasin_partition_result& partition,
    const lsdtt_raster_cpu&          dem,
    const lsdtt_raster_cpu&          gradient)
{
    static const char* who = "lsdtt::compute_subbasin_geodata";
    require_same_grid(partition, dem, who, "dem");
    require_same_grid(partition, gradient, who, "gradient");

    const std::int32_t n_sb      = partition.n_subbasins();
    const double       cell_area = partition.meta.cell_size * partition.meta.cell_size;

    const std::vector<float>&        z      = dem.get_data();
    const std::vector<std::uint8_t>& z_mask = dem.get_mask();
    const std::vector<float>&        g      = gradient.get_data();
    const std::vector<std::uint8_t>& g_mask = gradient.get_mask();

    const subbasin_pixel_rows rows = build_pixel_rows(partition);

    subbasin_geodata_result out;
    out.area_m2.resize(static_cast<std::size_t>(n_sb));
    out.elev_mean_m.resize(static_cast<std::size_t>(n_sb));
    out.elev_std_m.resize(static_cast<std::size_t>(n_sb));
    out.slope_mean_percent.resize(static_cast<std::size_t>(n_sb));
    out.n_valid_elev.resize(static_cast<std::size_t>(n_sb));
    out.n_valid_gradient.resize(static_cast<std::size_t>(n_sb));

    std::vector<double> elev_values;
    std::vector<double> gradient_values;

    for (std::int32_t s = 0; s < n_sb; ++s)
    {
        const std::size_t si = static_cast<std::size_t>(s);

        elev_values.clear();
        gradient_values.clear();
        for (std::int64_t k = rows.offset[si]; k < rows.offset[si + 1]; ++k)
        {
            const std::size_t i = static_cast<std::size_t>(rows.flats[static_cast<std::size_t>(k)]);
            if (z_mask[i] != 0) elev_values.push_back(static_cast<double>(z[i]));
            if (g_mask[i] != 0) gradient_values.push_back(static_cast<double>(g[i]));
        }

        out.area_m2[si]          = static_cast<double>(partition.n_pixels[si]) * cell_area;
        out.n_valid_elev[si]     = static_cast<std::int64_t>(elev_values.size());
        out.n_valid_gradient[si] = static_cast<std::int64_t>(gradient_values.size());

        if (elev_values.empty())
        {
            out.elev_mean_m[si] = k_nan;
            out.elev_std_m[si]  = k_nan;
        }
        else
        {
            out.elev_mean_m[si] = compute_mean(elev_values);
            out.elev_std_m[si]  = compute_stddev_population(elev_values);
        }

        // HYPE's slope_mean is a percentage: 100 * dimensionless gradient.
        out.slope_mean_percent[si] =
            gradient_values.empty() ? k_nan : 100.0 * compute_mean(gradient_values);
    }

    return out;
}

// =============================================================================
// Class fractions and elevation offsets
// =============================================================================

subbasin_class_fractions_result compute_subbasin_class_fractions(
    const subbasin_partition_result& partition,
    const lsdtt_raster_cpu&          class_raster,
    const lsdtt_raster_cpu&          dem)
{
    static const char* who = "lsdtt::compute_subbasin_class_fractions";
    require_same_grid(partition, class_raster, who, "class_raster");
    require_same_grid(partition, dem, who, "dem");

    const std::int32_t n_sb = partition.n_subbasins();

    const std::vector<float>&        cls      = class_raster.get_data();
    const std::vector<std::uint8_t>& cls_mask = class_raster.get_mask();
    const std::vector<float>&        z        = dem.get_data();
    const std::vector<std::uint8_t>& z_mask   = dem.get_mask();

    const subbasin_pixel_rows rows = build_pixel_rows(partition);

    struct accumulator
    {
        std::int64_t n_pixels = 0;
        double       z_sum    = 0.0;
        std::int64_t n_z      = 0;
    };

    subbasin_class_fractions_result out;
    out.offset.assign(static_cast<std::size_t>(n_sb) + 1, 0);
    out.n_valid_class.assign(static_cast<std::size_t>(n_sb), 0);

    std::vector<std::int32_t> all_codes;

    for (std::int32_t s = 0; s < n_sb; ++s)
    {
        const std::size_t si = static_cast<std::size_t>(s);

        std::map<std::int32_t, accumulator> per_class;
        double       z_sum_all = 0.0;
        std::int64_t n_z_all   = 0;

        for (std::int64_t k = rows.offset[si]; k < rows.offset[si + 1]; ++k)
        {
            const std::size_t i = static_cast<std::size_t>(rows.flats[static_cast<std::size_t>(k)]);
            const bool z_ok = (z_mask[i] != 0);

            if (z_ok)
            {
                z_sum_all += static_cast<double>(z[i]);
                ++n_z_all;
            }

            if (cls_mask[i] == 0) continue;

            const double v = static_cast<double>(cls[i]);
            if (!std::isfinite(v) || v != std::floor(v) ||
                v < static_cast<double>(std::numeric_limits<std::int32_t>::min()) ||
                v > static_cast<double>(std::numeric_limits<std::int32_t>::max()))
            {
                throw std::invalid_argument(
                    std::string(who) + ": class_raster value " + std::to_string(v) +
                    " at flat " + std::to_string(i) + " is not an integer class code");
            }

            accumulator& a = per_class[static_cast<std::int32_t>(v)];
            ++a.n_pixels;
            if (z_ok)
            {
                a.z_sum += static_cast<double>(z[i]);
                ++a.n_z;
            }
            ++out.n_valid_class[si];
        }

        const double elev_mean =
            (n_z_all > 0) ? z_sum_all / static_cast<double>(n_z_all) : k_nan;
        const double n_valid = static_cast<double>(out.n_valid_class[si]);

        for (const auto& [code, a] : per_class)
        {
            subbasin_class_entry e;
            e.subbasin   = s;
            e.class_code = code;
            e.n_pixels   = a.n_pixels;
            e.fraction   = static_cast<double>(a.n_pixels) / n_valid;
            e.dh_m       = (a.n_z > 0)
                               ? a.z_sum / static_cast<double>(a.n_z) - elev_mean
                               : k_nan;
            out.entries.push_back(e);
            all_codes.push_back(code);
        }

        out.offset[si + 1] = static_cast<std::int64_t>(out.entries.size());
    }

    std::sort(all_codes.begin(), all_codes.end());
    all_codes.erase(std::unique(all_codes.begin(), all_codes.end()), all_codes.end());
    out.class_codes = std::move(all_codes);

    return out;
}

// =============================================================================
// River lengths
// =============================================================================

subbasin_river_lengths_result compute_subbasin_river_lengths(
    const subbasin_partition_result& partition,
    const lsdtt_flow_base&           flow,
    std::int64_t                     channel_threshold_contributing_pixels)
{
    static const char* who = "lsdtt::compute_subbasin_river_lengths";

    if (flow.full_grid_size() != partition.subbasin_id.size())
    {
        throw std::invalid_argument(
            std::string(who) + ": flow object and partition are on different grids");
    }

    const std::int32_t n_sb      = partition.n_subbasins();
    const std::size_t  n         = partition.subbasin_id.size();
    const double       dx        = flow.cell_size();
    const double       dxy       = dx * std::sqrt(2.0);

    subbasin_river_lengths_result out;
    out.rivlen_m.assign(static_cast<std::size_t>(n_sb), 0.0);
    out.loc_rivlen_m.assign(static_cast<std::size_t>(n_sb), k_nan);
    out.trunk_top_flat.assign(static_cast<std::size_t>(n_sb), -1);
    out.n_trunk_pixels.assign(static_cast<std::size_t>(n_sb), 0);
    out.trunk_mask.assign(n, 0u);

    // Largest upstream sub-basin of each sub-basin, by contributing pixels at
    // its outlet (ties: lower index, from the ascending scan and strict >).
    std::vector<std::int32_t> largest_upstream(static_cast<std::size_t>(n_sb), -1);
    for (std::int32_t u = 0; u < n_sb; ++u)
    {
        const std::int32_t d = partition.downstream[static_cast<std::size_t>(u)];
        if (d < 0) continue;

        std::int32_t& best = largest_upstream[static_cast<std::size_t>(d)];
        if (best < 0 ||
            flow.s_vector_size_unchecked(partition.outlet_flat[static_cast<std::size_t>(u)]) >
            flow.s_vector_size_unchecked(partition.outlet_flat[static_cast<std::size_t>(best)]))
        {
            best = u;
        }
    }

    std::vector<std::int64_t> donors;

    for (std::int32_t s = 0; s < n_sb; ++s)
    {
        const std::size_t  si     = static_cast<std::size_t>(s);
        const std::int64_t outlet = partition.outlet_flat[si];

        // ---- upstream end of the trunk, p_in ----
        std::int64_t p_in = outlet;

        if (largest_upstream[si] >= 0)
        {
            p_in = flow.receiver_of_unchecked(
                partition.outlet_flat[static_cast<std::size_t>(largest_upstream[si])]);
        }
        else
        {
            while (true)
            {
                collect_donors(flow, p_in, donors);

                std::int64_t best      = -1;
                std::uint64_t best_cnt = 0;
                for (const std::int64_t d : donors)
                {
                    if (partition.subbasin_id[static_cast<std::size_t>(d)] != s) continue;

                    const std::uint64_t cnt = flow.s_vector_size_unchecked(d);
                    if (static_cast<std::int64_t>(cnt) < channel_threshold_contributing_pixels)
                        continue;
                    if (best < 0 || cnt > best_cnt || (cnt == best_cnt && d < best))
                    {
                        best     = d;
                        best_cnt = cnt;
                    }
                }

                if (best < 0) break;
                p_in = best;
            }
        }

        out.trunk_top_flat[si] = p_in;

        // ---- walk the trunk down to the outlet ----
        double       length  = 0.0;
        std::int64_t p       = p_in;
        std::int64_t n_trunk = 1;
        out.trunk_mask[static_cast<std::size_t>(p)] = 1u;

        while (p != outlet)
        {
            const std::uint8_t dir = flow.receiver_dir_of_unchecked(p);
            length += (dir % 2u == 0u) ? dx : dxy;
            p = flow.receiver_of_unchecked(p);
            out.trunk_mask[static_cast<std::size_t>(p)] = 1u;
            ++n_trunk;
        }

        out.rivlen_m[si]       = length;
        out.n_trunk_pixels[si] = n_trunk;
    }

    // ---- local river: mean flow distance to the trunk ----
    // HAND with flow distance standing in for elevation: FD(p) - FD(trunk pixel).
    const lsdtt_raster_result fd = flow.get_distance_from_outlet_raster();
    const lsdtt_raster_result to_trunk =
        flow.get_height_above_channel_mask_raster(fd.data, fd.mask, out.trunk_mask);

    std::vector<double>       sum(static_cast<std::size_t>(n_sb), 0.0);
    std::vector<std::int64_t> cnt(static_cast<std::size_t>(n_sb), 0);

    for (std::size_t i = 0; i < n; ++i)
    {
        const std::int32_t s = partition.subbasin_id[i];
        if (s < 0 || to_trunk.mask[i] == 0) continue;
        sum[static_cast<std::size_t>(s)] += static_cast<double>(to_trunk.data[i]);
        ++cnt[static_cast<std::size_t>(s)];
    }

    for (std::int32_t s = 0; s < n_sb; ++s)
    {
        const std::size_t si = static_cast<std::size_t>(s);
        if (cnt[si] > 0) out.loc_rivlen_m[si] = sum[si] / static_cast<double>(cnt[si]);
    }

    return out;
}

} // namespace lsdtt
