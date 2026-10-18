// lsdtt_channel_network_grouping_cpu.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html

#include "lsdtt_channel_network_grouping_cpu.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>

#include "../core/lsdtt_data_loader.hpp"
#include "../core/lsdtt_raster_warp.hpp"
#include "lsdtt_lithology_map_cpu.hpp"
#include "lsdtt_vector_rasterise_cpu.hpp"

namespace lsdtt
{
namespace
{

double quiet_nan()
{
    return std::numeric_limits<double>::quiet_NaN();
}

// A CSV field, double-quoted only when it holds a comma, quote or newline
// (lithologic unit names can), so every other field is written unchanged.
std::string csv_field(const std::string& s)
{
    if (s.find_first_of(",\"\n") == std::string::npos) return s;
    std::string q = "\"";
    for (const char c : s) { if (c == '"') q += "\"\""; else q += c; }
    return q + "\"";
}

std::string format_double(double v)
{
    if (!std::isfinite(v)) return std::isnan(v) ? "nan" : (v > 0.0 ? "inf" : "-inf");
    std::ostringstream os;
    os << std::setprecision(10) << v;
    return os.str();
}

double axial_of(double bearing_deg)
{
    if (!std::isfinite(bearing_deg)) return quiet_nan();
    double a = std::fmod(bearing_deg, 180.0);
    if (a < 0.0) a += 180.0;
    if (a >= 180.0) a = 0.0;
    return a;
}

} // namespace

// =============================================================================
// Zone labels
// =============================================================================
channel_network_zone_labels make_zone_labels_from_polygon_field(
    const lsdtt_raster_cpu& reference_dem,
    const std::string&      vector_path,
    const std::string&      layer_name,
    const std::string&      field_name,
    const std::string&      key_name)
{
    if (field_name.empty())
    {
        throw std::invalid_argument(
            "a vector " + key_name + " map needs a field name (channel_region_field_name)");
    }
    const polygon_layer_info info = read_polygon_layer_info(vector_path, layer_name);
    if (std::find(info.field_names.begin(), info.field_names.end(), field_name) ==
        info.field_names.end())
    {
        std::string list;
        for (const auto& f : info.field_names) list += (list.empty() ? "" : ", ") + f;
        throw std::invalid_argument(key_name + " map " + vector_path + " has no field '" +
                                    field_name + "'; its fields are: " + list);
    }

    polygon_field_rasterise_result r =
        rasterise_polygon_field_values(reference_dem, vector_path, layer_name, field_name);

    channel_network_zone_labels out;
    out.key_name = key_name;
    out.value_index = std::move(r.value_index);
    out.values = std::move(r.values);
    out.n_unmapped_valid = r.n_unmapped_valid;
    out.notes.push_back(
        "rasterised " + std::to_string(r.n_features_read) + " features (" +
        std::to_string(r.n_features_skipped) + " skipped), " +
        std::to_string(out.values.size()) + " distinct values of '" + field_name +
        "'; " + std::to_string(r.n_gap_filled) + " sliver pixels filled by coverage; " +
        std::to_string(r.n_overlap_pixels) + " overlapping pixels; " +
        std::to_string(r.n_unmapped_valid) + " valid pixels with no " + key_name);
    std::size_t n_valid = 0;
    for (const std::uint8_t m : reference_dem.get_mask()) n_valid += m ? 1 : 0;
    check_categorical_map_coverage(key_name, vector_path, n_valid, out.n_unmapped_valid,
                                   {info.x_min, info.y_min, info.x_max, info.y_max, info.crs_wkt},
                                   reference_dem.get_meta(), out.notes, out.warnings);
    return out;
}

channel_network_zone_labels make_zone_labels_from_code_raster(
    const lsdtt_raster_cpu& reference_dem,
    const std::string&      raster_path,
    const std::string&      key_name)
{
    const lsdtt_raster_meta& gm = reference_dem.get_meta();
    const std::size_t npix = static_cast<std::size_t>(gm.n_rows) * gm.n_cols;
    const std::vector<std::uint8_t>& gmask = reference_dem.get_mask();

    lsdtt_data_loader loader;
    loader.open_raster(raster_path);
    const lsdtt_raster_meta sm = loader.get_meta();
    const double tol = 1e-6 * gm.cell_size;
    const bool on_grid = sm.n_rows == gm.n_rows && sm.n_cols == gm.n_cols &&
                         std::abs(sm.cell_size - gm.cell_size) < tol &&
                         std::abs(sm.x_min - gm.x_min) < tol &&
                         std::abs(sm.y_min - gm.y_min) < tol &&
                         (sm.crs_wkt.empty() || gm.crs_wkt.empty() || sm.crs_wkt == gm.crs_wkt);

    channel_network_zone_labels out;
    out.key_name = key_name;

    lsdtt_raster_tile tile;
    if (on_grid)
    {
        tile = loader.read_full();
        out.notes.push_back(key_name + " code raster read on the DEM grid as it is");
    }
    else
    {
        if (sm.crs_wkt.empty() || gm.crs_wkt.empty())
        {
            std::ostringstream msg;
            msg << key_name << " code raster " << raster_path
                << " is not on the DEM grid (" << sm.n_rows << "x" << sm.n_cols
                << " at " << sm.cell_size << " against " << gm.n_rows << "x" << gm.n_cols
                << " at " << gm.cell_size << ") and cannot be warped onto it: "
                << (sm.crs_wkt.empty() ? "it has no CRS" : "the DEM has no CRS");
            throw std::invalid_argument(msg.str());
        }
        lsdtt_warp_options wo;
        wo.resampling_method = "near";
        try
        {
            tile = warp_raster_to_grid_in_memory(raster_path, gm, wo);
        }
        catch (const std::exception& e)
        {
            if (std::string(e.what()).find("do not overlap") == std::string::npos) throw;
            std::size_t n_valid = 0;
            for (const std::uint8_t v : gmask) n_valid += v ? 1 : 0;
            check_categorical_map_coverage(
                key_name, raster_path, n_valid, n_valid,
                {sm.x_min, sm.y_min, sm.x_min + sm.n_cols * sm.cell_size,
                 sm.y_min + sm.n_rows * sm.cell_size, sm.crs_wkt},
                gm, out.notes, out.warnings);   // throws
            throw;
        }
        out.notes.push_back(key_name + " code raster warped onto the DEM grid (nearest neighbour)");
    }
    if (tile.data.size() != npix)
    {
        throw std::runtime_error(key_name + " code raster: read " +
                                 std::to_string(tile.data.size()) +
                                 " pixels for a DEM grid of " + std::to_string(npix));
    }

    std::map<long, std::int32_t> seen;
    std::vector<long> code(npix, -1);
    for (std::size_t i = 0; i < npix; ++i)
    {
        if (!gmask[i] || !tile.mask[i]) continue;
        const float v = tile.data[i];
        if (!(v >= 0.0f) || v != std::floor(v))
        {
            std::ostringstream msg;
            msg << key_name << " code raster: value " << v << " at pixel ("
                << i / gm.n_cols << ", " << i % gm.n_cols << ") is not an integer code >= 0";
            throw std::invalid_argument(msg.str());
        }
        code[i] = static_cast<long>(v);
        seen.emplace(code[i], 0);
    }
    for (auto& [c, idx] : seen)
    {
        idx = static_cast<std::int32_t>(out.values.size());
        out.values.push_back(std::to_string(c));
    }
    out.value_index.assign(npix, -1);
    for (std::size_t i = 0; i < npix; ++i)
    {
        if (code[i] >= 0) out.value_index[i] = seen[code[i]];
        else if (gmask[i]) ++out.n_unmapped_valid;
    }
    out.notes.push_back(std::to_string(out.values.size()) + " distinct " + key_name +
                        " codes; " + std::to_string(out.n_unmapped_valid) +
                        " valid pixels with no " + key_name);
    std::size_t n_valid = 0;
    for (const std::uint8_t m : gmask) n_valid += m ? 1 : 0;
    check_categorical_map_coverage(
        key_name, raster_path, n_valid, out.n_unmapped_valid,
        {sm.x_min, sm.y_min, sm.x_min + sm.n_cols * sm.cell_size,
         sm.y_min + sm.n_rows * sm.cell_size, sm.crs_wkt},
        gm, out.notes, out.warnings);
    return out;
}

channel_network_zone_labels read_channel_region_labels(
    const lsdtt_raster_cpu& reference_dem,
    const std::string&      path,
    const std::string&      source,
    const std::string&      layer_name,
    const std::string&      field_name)
{
    std::string resolved;
    try
    {
        resolved = lithology_source_from_path(path, source);
    }
    catch (const std::invalid_argument&)
    {
        throw std::invalid_argument(
            "channel_region_source must be auto, vector_map or code_raster (got '" +
            source + "')");
    }
    return resolved == "code_raster"
        ? make_zone_labels_from_code_raster(reference_dem, path, "region")
        : make_zone_labels_from_polygon_field(reference_dem, path, layer_name,
                                              field_name, "region");
}

channel_zone_majority majority_zone_of_flats(
    const std::vector<std::int64_t>&  flats,
    const std::vector<std::int32_t>&  value_index)
{
    channel_zone_majority out;
    out.n_flats = flats.size();
    if (flats.empty()) return out;

    std::map<std::int32_t, std::size_t> counts;
    for (const std::int64_t f : flats)
    {
        const std::int32_t v =
            (f >= 0 && static_cast<std::size_t>(f) < value_index.size())
            ? value_index[static_cast<std::size_t>(f)] : -1;
        ++counts[v < 0 ? -1 : v];
    }
    out.n_distinct = counts.size();
    std::size_t best = 0;
    for (const auto& [v, c] : counts)   // ascending, so ties keep the lower index
    {
        if (c > best)
        {
            best = c;
            out.value_index = v;
        }
    }
    out.fraction = static_cast<double>(best) / static_cast<double>(flats.size());
    return out;
}

channel_record_zone_assignment assign_zones_to_channel_links(
    const lsdtt_flow_base&                         flow,
    const lsdtt_junction_network_base&             network,
    const std::vector<lsdtt_channel_link_bearing>& records,
    const channel_network_zone_labels&             labels)
{
    channel_record_zone_assignment out;
    out.value_index.assign(records.size(), -1);
    out.purity.assign(records.size(), quiet_nan());
    for (std::size_t i = 0; i < records.size(); ++i)
    {
        const lsdtt_traced_channel_link link =
            trace_channel_link_between_junctions(flow, network, records[i].upstream_junction);
        const channel_zone_majority m = majority_zone_of_flats(link.flats, labels.value_index);
        out.value_index[i] = m.value_index;
        out.purity[i] = m.fraction;
    }
    return out;
}

channel_record_zone_assignment assign_zones_to_junctions(
    const lsdtt_junction_network_base& network,
    const std::vector<std::int32_t>&   junctions,
    const channel_network_zone_labels& labels)
{
    channel_record_zone_assignment out;
    out.value_index.assign(junctions.size(), -1);
    out.purity.assign(junctions.size(), quiet_nan());
    for (std::size_t i = 0; i < junctions.size(); ++i)
    {
        if (!network.is_valid_junction(junctions[i])) continue;
        const std::int64_t flat = network.flat_of_junction(junctions[i]);
        const std::int32_t v =
            (flat >= 0 && static_cast<std::size_t>(flat) < labels.value_index.size())
            ? labels.value_index[static_cast<std::size_t>(flat)] : -1;
        out.value_index[i] = v < 0 ? -1 : v;
        out.purity[i] = 1.0;
    }
    return out;
}

channel_record_categories categories_from_zone_assignment(
    const std::string&                    key,
    const channel_record_zone_assignment& assignment,
    const channel_network_zone_labels&    labels,
    double                                min_purity)
{
    channel_record_categories out;
    out.key = key;
    out.category.assign(assignment.value_index.size(), -1);
    for (std::size_t i = 0; i < assignment.value_index.size(); ++i)
    {
        const std::int32_t v = assignment.value_index[i];
        if (v < 0 || static_cast<std::size_t>(v) >= labels.values.size()) continue;
        const double p = assignment.purity[i];
        if (!(p >= min_purity)) continue;
        out.category[i] = v;
        out.names[v] = labels.values[static_cast<std::size_t>(v)];
    }
    return out;
}

// =============================================================================
// Orders
// =============================================================================
channel_record_categories strahler_order_categories(
    const std::vector<std::int32_t>& orders,
    int                              max_order_group,
    const std::string&               key)
{
    channel_record_categories out;
    out.key = key;
    out.category.assign(orders.size(), -1);
    for (std::size_t i = 0; i < orders.size(); ++i)
    {
        std::int32_t o = orders[i];
        if (o < 1) continue;
        if (max_order_group > 0 && o >= max_order_group)
        {
            o = max_order_group;
            out.names[o] = std::to_string(o) + "+";
        }
        else
        {
            out.names[o] = std::to_string(o);
        }
        out.category[i] = o;
    }
    return out;
}

channel_record_categories donor_order_pair_categories(
    const std::vector<std::int32_t>& donor_1_order,
    const std::vector<std::int32_t>& donor_2_order,
    const std::string&               key)
{
    if (donor_1_order.size() != donor_2_order.size())
    {
        throw std::invalid_argument("donor_order_pair_categories: sizes differ");
    }
    channel_record_categories out;
    out.key = key;
    out.category.assign(donor_1_order.size(), -1);
    for (std::size_t i = 0; i < donor_1_order.size(); ++i)
    {
        const std::int32_t a = std::min(donor_1_order[i], donor_2_order[i]);
        const std::int32_t b = std::max(donor_1_order[i], donor_2_order[i]);
        if (a < 1) continue;
        const std::int32_t c = a * 1000 + b;
        out.category[i] = c;
        out.names[c] = std::to_string(a) + "_" + std::to_string(b);
    }
    return out;
}

// =============================================================================
// Basins
// =============================================================================
std::vector<std::int32_t> outlet_basin_of_junctions(
    const lsdtt_junction_network_base& network,
    const std::vector<std::int32_t>&   outlet_junctions)
{
    const std::int32_t nj = network.n_junctions();
    std::vector<std::int32_t> basin(static_cast<std::size_t>(std::max(0, nj)), -1);

    std::vector<std::int32_t> outlets;
    for (const std::int32_t o : outlet_junctions)
    {
        if (network.is_valid_junction(o)) outlets.push_back(o);
    }
    std::sort(outlets.begin(), outlets.end());
    outlets.erase(std::unique(outlets.begin(), outlets.end()), outlets.end());

    // Largest first, so a nested (smaller) outlet overwrites its share.
    std::stable_sort(outlets.begin(), outlets.end(),
        [&](std::int32_t a, std::int32_t b)
        { return network.s_vector_size(a) > network.s_vector_size(b); });

    for (const std::int32_t o : outlets)
    {
        basin[static_cast<std::size_t>(o)] = o;
        const lsdtt_junction_index_slice up = network.upslope_junctions(o);
        for (const std::int32_t* it = up.begin(); it != up.end(); ++it)
        {
            if (*it >= 0 && *it < nj) basin[static_cast<std::size_t>(*it)] = o;
        }
    }
    return basin;
}

std::vector<std::int32_t> order_k_basin_of_junctions(
    const lsdtt_junction_network_base& network,
    int                                k)
{
    const std::int32_t nj = network.n_junctions();
    std::vector<std::int32_t> basin(static_cast<std::size_t>(std::max(0, nj)), -1);
    std::vector<std::uint8_t> done(basin.size(), 0);
    if (k < 1) return basin;

    // basin(j) = none                       if order(j) > k
    //          = basin(receiver(j))         if receiver valid and order(r) <= k
    //          = j if order(j) == k else none, otherwise.
    std::vector<std::int32_t> chain;
    for (std::int32_t start = 0; start < nj; ++start)
    {
        if (done[static_cast<std::size_t>(start)]) continue;
        chain.clear();
        std::int32_t j = start;
        std::int32_t result = -1;
        while (true)
        {
            if (done[static_cast<std::size_t>(j)])
            {
                result = basin[static_cast<std::size_t>(j)];
                break;
            }
            chain.push_back(j);
            const std::int32_t oj = network.stream_order_of_junction(j);
            if (oj > k)
            {
                result = -1;
                break;
            }
            const std::int32_t r = network.is_outlet_junction(j)
                ? k_nodata_junction : network.receiver_of_junction(j);
            if (!network.is_valid_junction(r) || r == j ||
                network.stream_order_of_junction(r) > k)
            {
                result = oj == k ? j : -1;
                break;
            }
            j = r;
        }
        // Every junction on the chain shares the result, except those of
        // order > k (only possible as the first element), which stay -1.
        for (const std::int32_t c : chain)
        {
            basin[static_cast<std::size_t>(c)] =
                network.stream_order_of_junction(c) > k ? -1 : result;
            done[static_cast<std::size_t>(c)] = 1;
        }
    }
    return basin;
}

channel_record_categories basin_categories(
    const std::vector<std::int32_t>& record_junctions,
    const std::vector<std::int32_t>& basin_of_junction,
    const std::string&               key)
{
    channel_record_categories out;
    out.key = key;
    out.category.assign(record_junctions.size(), -1);
    for (std::size_t i = 0; i < record_junctions.size(); ++i)
    {
        const std::int32_t j = record_junctions[i];
        if (j < 0 || static_cast<std::size_t>(j) >= basin_of_junction.size()) continue;
        const std::int32_t b = basin_of_junction[static_cast<std::size_t>(j)];
        if (b < 0) continue;
        out.category[i] = b;
        out.names[b] = std::to_string(b);
    }
    return out;
}

std::vector<std::int64_t> block_ids_of_records(
    const std::vector<std::int32_t>& record_junctions,
    const std::vector<std::int32_t>& basin_of_junction)
{
    std::vector<std::int64_t> out(record_junctions.size(), -1);
    for (std::size_t i = 0; i < record_junctions.size(); ++i)
    {
        const std::int32_t j = record_junctions[i];
        if (j < 0 || static_cast<std::size_t>(j) >= basin_of_junction.size()) continue;
        out[i] = basin_of_junction[static_cast<std::size_t>(j)];
    }
    return out;
}

// =============================================================================
// Membership CSV
// =============================================================================
void write_channel_bearing_group_membership_csv(
    const std::vector<lsdtt_channel_link_bearing>& records,
    const channel_record_categories&              groups,
    const channel_record_categories*              strata,
    const std::vector<std::int32_t>&              basin_of_record,
    const std::vector<std::int64_t>&              block_of_record,
    const channel_record_zone_assignment*         region,
    const channel_network_zone_labels*            region_labels,
    const std::string&                            out_path)
{
    const std::size_t n = records.size();
    if (groups.category.size() != n || basin_of_record.size() != n ||
        block_of_record.size() != n ||
        (strata != nullptr && strata->category.size() != n) ||
        (region != nullptr && region->value_index.size() != n))
    {
        throw std::invalid_argument(
            "write_channel_bearing_group_membership_csv: per-record vectors differ in size");
    }

    std::ofstream os(out_path);
    if (!os)
    {
        throw std::runtime_error(
            "write_channel_bearing_group_membership_csv: could not open " + out_path);
    }

    os << "upstream_junction,downstream_junction,stream_order,total_length_m,"
          "bearing_deg,axial_bearing_deg,group_value,stratum_value,basin_id,"
          "block_id,region_value,region_purity\n";
    for (std::size_t i = 0; i < n; ++i)
    {
        const lsdtt_channel_link_bearing& r = records[i];
        const std::int32_t g = groups.category[i];
        std::string stratum = "all";
        if (strata != nullptr)
        {
            const std::int32_t s = strata->category[i];
            stratum = s >= 0 ? strata->name_of(s) : std::string();
        }
        std::string region_value;
        double purity = quiet_nan();
        if (region != nullptr && region_labels != nullptr)
        {
            const std::int32_t v = region->value_index[i];
            if (v >= 0 && static_cast<std::size_t>(v) < region_labels->values.size())
            {
                region_value = region_labels->values[static_cast<std::size_t>(v)];
            }
            purity = region->purity[i];
        }

        os << r.upstream_junction << ',' << r.downstream_junction << ','
           << r.stream_order << ',' << format_double(r.total_length_metres) << ','
           << format_double(r.bearing_deg) << ',' << format_double(axial_of(r.bearing_deg)) << ','
           << csv_field(g >= 0 ? groups.name_of(g) : std::string()) << ',' << csv_field(stratum) << ','
           << basin_of_record[i] << ',' << block_of_record[i] << ','
           << csv_field(region_value) << ',' << format_double(purity) << '\n';
    }
}

std::vector<std::size_t> count_order_k_basins_with_records(
    const lsdtt_junction_network_base& network,
    const std::vector<std::int32_t>&   record_junctions)
{
    int max_order = 0;
    for (std::int32_t j = 0; j < network.n_junctions(); ++j)
        max_order = std::max(max_order, static_cast<int>(network.stream_order_of_junction(j)));
    std::vector<std::size_t> counts;
    for (int k = 1; k <= max_order; ++k)
    {
        const std::vector<std::int32_t> basin = order_k_basin_of_junctions(network, k);
        std::set<std::int32_t> seen;
        for (const std::int32_t j : record_junctions)
            if (j >= 0 && static_cast<std::size_t>(j) < basin.size() &&
                basin[static_cast<std::size_t>(j)] >= 0)
                seen.insert(basin[static_cast<std::size_t>(j)]);
        counts.push_back(seen.size());
    }
    return counts;
}

channel_block_order_choice choose_channel_block_order(
    const std::vector<std::size_t>& n_basins_by_order,
    int                             override_order,
    int                             min_blocks)
{
    channel_block_order_choice c;
    c.n_basins_by_order = n_basins_by_order;
    c.min_blocks = min_blocks;
    auto count_at = [&](int k) -> std::size_t {
        return k >= 1 && static_cast<std::size_t>(k) <= n_basins_by_order.size()
                   ? n_basins_by_order[static_cast<std::size_t>(k - 1)] : 0;
    };
    if (override_order > 0)
    {
        c.automatic = false;
        c.block_order = override_order;
        if (count_at(override_order) < static_cast<std::size_t>(std::max(0, min_blocks)))
            c.note = "block order " + std::to_string(override_order) + " (override) has " +
                     std::to_string(count_at(override_order)) + " basin(s) with records, fewer than " +
                     "channel_bearing_min_blocks = " + std::to_string(min_blocks);
        return c;
    }
    std::size_t largest = 0;
    for (std::size_t i = 0; i < n_basins_by_order.size(); ++i)
    {
        largest = std::max(largest, n_basins_by_order[i]);
        if (n_basins_by_order[i] >= static_cast<std::size_t>(std::max(0, min_blocks)))
            c.block_order = static_cast<int>(i + 1);
    }
    if (c.block_order == 0)
        c.note = "network too small for block permutation: largest order-k basin count " +
                 std::to_string(largest) + " < min_blocks " + std::to_string(min_blocks);
    return c;
}

void write_channel_network_block_order_counts_csv(const channel_block_order_choice& choice,
                                                  const std::string& out_path)
{
    std::ofstream os(out_path);
    if (!os) throw std::runtime_error("write_channel_network_block_order_counts_csv: could not open " + out_path);
    os << "strahler_order,n_basins_with_records,chosen\n";
    for (std::size_t i = 0; i < choice.n_basins_by_order.size(); ++i)
        os << (i + 1) << ',' << choice.n_basins_by_order[i] << ','
           << (static_cast<int>(i + 1) == choice.block_order ? 1 : 0) << '\n';
}

} // namespace lsdtt
