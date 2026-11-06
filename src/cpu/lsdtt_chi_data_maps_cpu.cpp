// lsdtt_chi_data_maps_cpu.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// Implementation of write_chi_data_maps().
// See lsdtt_chi_data_maps_cpu.hpp for API documentation.
// -----------------------------------------------------------------------------

#include "lsdtt_chi_data_maps_cpu.hpp"

#include <cstdint>
#include <iostream>
#include <limits>
#include <optional>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

#include "../core/lsdtt_data_loader.hpp"
#include "../core/lsdtt_flow_base.hpp"
#include "../core/lsdtt_types.hpp"
#include "../core/lsdtt_vector_utils.hpp"
#include "lsdtt_basin_cpu.hpp"
#include "lsdtt_basin_outline_cpu.hpp"
#include "lsdtt_chi_extraction_cpu.hpp"
#include "lsdtt_raster_cpu.hpp"

namespace lsdtt
{

namespace
{

std::string high_prec(double v)
{
    std::ostringstream ss;
    ss.precision(17);
    ss << v;
    return ss.str();
}

} // anonymous namespace

std::vector<lsdtt_chi_basin_source_keys> assign_chi_source_keys(
    const lsdtt_chi_collection_result& chi_networks,
    const lsdtt_flow_base&             flow,
    int                                n_cols,
    double                             cell_size,
    bool                               verbose)
{
    std::vector<lsdtt_chi_basin_source_keys> out(chi_networks.basins.size());
    std::int32_t next_basin_key  = 0;
    std::int32_t next_source_key = 0;
    std::size_t  n_complete      = 0;

    for (std::size_t bi = 0; bi < chi_networks.basins.size(); ++bi)
    {
        const auto& basin = chi_networks.basins[bi];
        auto&       keys  = out[bi];
        const std::size_t n = basin.flats.size();

        if (n == 0)
        {
            if (verbose) std::cout << "  basin[" << bi << "] has no nodes; skipping.\n";
            continue;
        }

        lsdtt_basin_cpu basin_obj(basin.outlet_flat, flow);
        const bool complete = basin_obj.check_basin_completeness();
        if (verbose)
        {
            if (!complete)
                std::cout << "  basin[" << bi << "] outlet_flat=" << basin.outlet_flat
                          << " is incomplete (upslope area touches nodata); skipping.\n";
            else
                std::cout << "  basin[" << bi << "]  outlet_flat=" << basin.outlet_flat
                          << "  n_nodes=" << n << "  [complete]\n";
        }

        // Local index map and donor lists. Receiver-before-donor ordering:
        // flats[0] is the outlet and every other node's receiver has a lower index.
        std::unordered_map<std::int64_t, std::size_t> local_idx;
        local_idx.reserve(n);
        for (std::size_t i = 0; i < n; ++i)
            local_idx[basin.flats[i]] = i;

        std::vector<std::vector<std::size_t>> donors(n);
        for (std::size_t i = 1; i < n; ++i)
            donors[local_idx.at(basin.receiver_flats[i])].push_back(i);

        keys.flow_distance_m.assign(n, 0.0);
        for (std::size_t i = 1; i < n; ++i)
        {
            const std::int64_t rec = basin.receiver_flats[i];
            keys.flow_distance_m[i] = keys.flow_distance_m[local_idx.at(rec)] +
                                      d8_step_length(basin.flats[i], rec, n_cols, cell_size);
        }
        keys.drainage_area_m2.resize(n);
        for (std::size_t i = 0; i < n; ++i)
            keys.drainage_area_m2[i] =
                static_cast<double>(flow.s_vector_size_unchecked(basin.flats[i])) *
                cell_size * cell_size;

        keys.is_source.assign(n, 0);
        for (std::size_t i = 0; i < n; ++i)
            if (donors[i].empty()) keys.is_source[i] = 1;

        keys.source_key.assign(n, -1);
        if (!complete) continue;
        ++n_complete;
        keys.basin_key = next_basin_key++;

        // Primary donor: longest upstream path. Upstream-first (reverse index)
        // so every donor is done before its receiver.
        constexpr std::size_t k_no_donor = std::numeric_limits<std::size_t>::max();
        std::vector<double>      max_path_up(n, 0.0);
        std::vector<std::size_t> primary_donor(n, k_no_donor);
        for (std::size_t ii = n; ii-- > 0; )
        {
            if (donors[ii].empty()) continue;
            double      best   = -1.0;
            std::size_t best_d = k_no_donor;
            for (std::size_t d : donors[ii])
            {
                const double candidate =
                    d8_step_length(basin.flats[d], basin.flats[ii], n_cols, cell_size) +
                    max_path_up[d];
                if (best_d == k_no_donor || candidate > best)
                {
                    best   = candidate;
                    best_d = d;
                }
            }
            max_path_up[ii]   = best;
            primary_donor[ii] = best_d;
        }
        for (std::size_t ii = n; ii-- > 0; )
            keys.source_key[ii] = donors[ii].empty() ? next_source_key++
                                                     : keys.source_key[primary_donor[ii]];
    }

    if (verbose)
    {
        std::cout << "  Complete basins: " << n_complete
                  << " / " << chi_networks.basins.size() << "\n";
        std::cout.flush();
    }
    return out;
}

lsdtt_chi_fit_source_key_map build_chi_fit_source_key_map(
    const lsdtt_chi_collection_result&               chi_networks,
    const std::vector<lsdtt_chi_basin_source_keys>& keys)
{
    lsdtt_chi_fit_source_key_map out(chi_networks.basins.size());
    for (std::size_t bi = 0; bi < chi_networks.basins.size() && bi < keys.size(); ++bi)
    {
        const auto& flats = chi_networks.basins[bi].flats;
        out[bi].reserve(flats.size());
        for (std::size_t i = 0; i < flats.size() && i < keys[bi].source_key.size(); ++i)
            out[bi][flats[i]] = keys[bi].source_key[i];
    }
    return out;
}

void write_chi_data_maps(
    const lsdtt_chi_collection_result& chi_networks,
    const lsdtt_flow_base&             flow,
    const lsdtt_raster_cpu&            dem,
    const std::string&                 output_prefix,
    point_output_format                fmt)
{
    const bool   is_csv  = (fmt == point_output_format::csv);
    const std::string pt_ext  = output_extension(fmt);

    // Basin polygons must be a binary vector format.
    const vector_output_format poly_fmt =
        is_csv ? vector_output_format::flatgeobuf
               : convert_point_to_vector_output_format(fmt);
    const std::string poly_ext = output_extension(poly_fmt);

    const lsdtt_raster_meta& meta      = dem.get_meta();
    const int                n_cols    = meta.n_cols;
    const double             cell_size = meta.cell_size;

    // Coordinate transform for CSV output (projected -> EPSG:4326).
    std::optional<lsdtt_coordinate_transform> ct;
    if (is_csv)
    {
        if (meta.crs_wkt.empty())
        {
            std::cerr << "FATAL (write_chi_data_maps): csv format requires a raster "
                         "CRS so coordinates can be transformed to latitude/longitude.\n";
            return;
        }
        ct.emplace(lsdtt_coordinate_transform::from_epsg(meta.crs_wkt, 4326));
    }

    // Helpers: flat -> projected (x, y) or (lat, lon).
    auto flat_to_xy = [&](std::int64_t flat) -> lsdtt_xy
    {
        const int row = static_cast<int>(flat / static_cast<std::int64_t>(n_cols));
        const int col = static_cast<int>(flat % static_cast<std::int64_t>(n_cols));
        return {
            meta.x_min + (col + 0.5) * cell_size,
            meta.y_min + (meta.n_rows - row - 0.5) * cell_size
        };
    };

    auto flat_to_latlon = [&](std::int64_t flat) -> std::pair<double, double>
    {
        const int row = static_cast<int>(flat / static_cast<std::int64_t>(n_cols));
        const int col = static_cast<int>(flat % static_cast<std::int64_t>(n_cols));
        const lsdtt_lonlat ll = cell_center_to_lonlat(meta, row, col, *ct);
        return {ll.latitude, ll.longitude};
    };

    // Intermediate row types.
    struct NodeRow
    {
        lsdtt_xy     xy;
        double       lat = 0.0, lon = 0.0;
        float        chi = 0.0f;
        double       drainage_area = 0.0;
        double       flow_distance = 0.0;
        float        elevation = 0.0f;
        std::int32_t source_key = -1;
        std::int32_t basin_key  = -1;
    };

    struct SourceRow
    {
        lsdtt_xy     xy;
        double       lat = 0.0, lon = 0.0;
        std::int32_t source_key = -1;
        std::int32_t basin_key  = -1;
    };

    struct OutletRow
    {
        lsdtt_xy     xy;
        double       lat = 0.0, lon = 0.0;
        std::int32_t basin_key = -1;
    };

    std::vector<NodeRow>   all_nodes;
    std::vector<SourceRow> all_sources;
    std::vector<OutletRow> all_outlets;

    const std::vector<lsdtt_chi_basin_source_keys> basin_keys =
        assign_chi_source_keys(chi_networks, flow, n_cols, cell_size, true);
    std::int32_t next_basin_key = 0;

    const std::vector<float>& dem_data = dem.get_data();

    for (std::size_t bi = 0; bi < chi_networks.basins.size(); ++bi)
    {
        const lsdtt_chi_basin_source_keys& keys = basin_keys[bi];
        if (keys.basin_key < 0) continue;
        const auto&        basin     = chi_networks.basins[bi];
        const std::int32_t basin_key = keys.basin_key;
        next_basin_key = basin_key + 1;
        const std::size_t  n         = basin.flats.size();
        const std::vector<double>&       flow_dist  = keys.flow_distance_m;
        const std::vector<std::int32_t>& src_key_of = keys.source_key;
        const std::vector<std::uint8_t>& is_source  = keys.is_source;

        // Outlet row.
        {
            const lsdtt_xy xy = flat_to_xy(basin.outlet_flat);
            double lat = 0.0, lon = 0.0;
            if (is_csv)
            {
                auto [la, lo] = flat_to_latlon(basin.outlet_flat);
                lat = la; lon = lo;
            }
            all_outlets.push_back({xy, lat, lon, basin_key});
        }

        // Per-node rows.
        for (std::size_t i = 0; i < n; ++i)
        {
            const std::int64_t flat = basin.flats[i];
            const lsdtt_xy     xy   = flat_to_xy(flat);

            double lat = 0.0, lon = 0.0;
            if (is_csv)
            {
                auto [la, lo] = flat_to_latlon(flat);
                lat = la; lon = lo;
            }

            const double da = keys.drainage_area_m2[i];

            const float elev = dem_data[static_cast<std::size_t>(flat)];

            all_nodes.push_back({
                xy, lat, lon,
                basin.chi[i], da, flow_dist[i], elev,
                src_key_of[i], basin_key
            });

            if (is_source[i])
                all_sources.push_back({xy, lat, lon, src_key_of[i], basin_key});
        }
    }

    std::cout << "  total nodes="  << all_nodes.size()
              << "  sources="      << all_sources.size()
              << "  valid_basins=" << next_basin_key << "\n";
    std::cout.flush();

    // ---- Write chi data maps ----
    try
    {
        const std::string fname = output_prefix + "_chi_data_maps" + pt_ext;
        std::cout << "  Writing: " << fname << "\n";
        std::cout.flush();

        if (is_csv)
        {
            const std::vector<std::string> hdrs = {
                "latitude", "longitude", "chi", "drainage_area",
                "flow_distance", "elevation", "source_key", "basin_key"
            };
            std::vector<std::vector<std::string>> rows;
            rows.reserve(all_nodes.size());
            for (const auto& nr : all_nodes)
            {
                rows.push_back({
                    high_prec(nr.lat), high_prec(nr.lon),
                    high_prec(static_cast<double>(nr.chi)),
                    high_prec(nr.drainage_area),
                    high_prec(nr.flow_distance),
                    high_prec(static_cast<double>(nr.elevation)),
                    std::to_string(nr.source_key),
                    std::to_string(nr.basin_key)
                });
            }
            lsdtt_data_loader writer;
            writer.write_csv(fname, hdrs, rows);
        }
        else
        {
            lsdtt_vector_schema sch;
            sch.field_names = {
                "chi", "drainage_area", "flow_distance",
                "elevation", "source_key", "basin_key"
            };
            sch.field_types = {
                lsdtt_vector_field_type::double_type,
                lsdtt_vector_field_type::double_type,
                lsdtt_vector_field_type::double_type,
                lsdtt_vector_field_type::double_type,
                lsdtt_vector_field_type::int64_type,
                lsdtt_vector_field_type::int64_type
            };
            const vector_output_format vfmt = convert_point_to_vector_output_format(fmt);
            lsdtt_data_loader writer;
            writer.begin_vector_write(
                fname, "chi_data_maps", sch, meta.crs_wkt, vfmt, false);
            for (const auto& nr : all_nodes)
            {
                writer.write_point(
                    nr.xy.x, nr.xy.y,
                    {
                        static_cast<double>(nr.chi),
                        nr.drainage_area,
                        nr.flow_distance,
                        static_cast<double>(nr.elevation),
                        static_cast<std::int64_t>(nr.source_key),
                        static_cast<std::int64_t>(nr.basin_key)
                    }
                );
            }
            writer.finish_vector_write();
        }
    }
    catch (const std::exception& e)
    {
        std::cerr << "Warning (write_chi_data_maps): node file failed: "
                  << e.what() << "\n";
    }

    // ---- Write basin outlet points ----
    try
    {
        const std::string fname = output_prefix + "_chi_basin_outlets" + pt_ext;
        std::cout << "  Writing: " << fname << "\n";
        std::cout.flush();

        if (is_csv)
        {
            const std::vector<std::string> hdrs = {"latitude", "longitude", "basin_key"};
            std::vector<std::vector<std::string>> rows;
            rows.reserve(all_outlets.size());
            for (const auto& or_ : all_outlets)
            {
                rows.push_back({
                    high_prec(or_.lat), high_prec(or_.lon),
                    std::to_string(or_.basin_key)
                });
            }
            lsdtt_data_loader writer;
            writer.write_csv(fname, hdrs, rows);
        }
        else
        {
            lsdtt_vector_schema sch;
            sch.field_names = {"basin_key"};
            sch.field_types = {lsdtt_vector_field_type::int64_type};
            const vector_output_format vfmt = convert_point_to_vector_output_format(fmt);
            lsdtt_data_loader writer;
            writer.begin_vector_write(
                fname, "chi_basin_outlets", sch, meta.crs_wkt, vfmt, false);
            for (const auto& or_ : all_outlets)
            {
                writer.write_point(
                    or_.xy.x, or_.xy.y,
                    {static_cast<std::int64_t>(or_.basin_key)}
                );
            }
            writer.finish_vector_write();
        }
    }
    catch (const std::exception& e)
    {
        std::cerr << "Warning (write_chi_data_maps): outlet file failed: "
                  << e.what() << "\n";
    }

    // ---- Write source points ----
    try
    {
        const std::string fname = output_prefix + "_chi_sources" + pt_ext;
        std::cout << "  Writing: " << fname << "\n";
        std::cout.flush();

        if (is_csv)
        {
            const std::vector<std::string> hdrs = {
                "latitude", "longitude", "source_key", "basin_key"
            };
            std::vector<std::vector<std::string>> rows;
            rows.reserve(all_sources.size());
            for (const auto& sr : all_sources)
            {
                rows.push_back({
                    high_prec(sr.lat), high_prec(sr.lon),
                    std::to_string(sr.source_key),
                    std::to_string(sr.basin_key)
                });
            }
            lsdtt_data_loader writer;
            writer.write_csv(fname, hdrs, rows);
        }
        else
        {
            lsdtt_vector_schema sch;
            sch.field_names = {"source_key", "basin_key"};
            sch.field_types = {
                lsdtt_vector_field_type::int64_type,
                lsdtt_vector_field_type::int64_type
            };
            const vector_output_format vfmt = convert_point_to_vector_output_format(fmt);
            lsdtt_data_loader writer;
            writer.begin_vector_write(
                fname, "chi_sources", sch, meta.crs_wkt, vfmt, false);
            for (const auto& sr : all_sources)
            {
                writer.write_point(
                    sr.xy.x, sr.xy.y,
                    {
                        static_cast<std::int64_t>(sr.source_key),
                        static_cast<std::int64_t>(sr.basin_key)
                    }
                );
            }
            writer.finish_vector_write();
        }
    }
    catch (const std::exception& e)
    {
        std::cerr << "Warning (write_chi_data_maps): sources file failed: "
                  << e.what() << "\n";
    }

    // ---- Write basin polygons ----
    try
    {
        const std::string fname = output_prefix + "_chi_basin_polygons" + poly_ext;
        std::cout << "  Writing: " << fname << "\n";
        std::cout.flush();

        lsdtt_vector_schema sch;
        sch.field_names = {"basin_key"};
        sch.field_types = {lsdtt_vector_field_type::int64_type};

        lsdtt_data_loader pw;
        pw.begin_vector_write(
            fname, "chi_basin_polygons", sch, meta.crs_wkt, poly_fmt, false);

        std::int32_t bk = 0;
        for (std::size_t bi = 0; bi < basin_keys.size(); ++bi)
        {
            if (basin_keys[bi].basin_key < 0) continue;
            try
            {
                lsdtt_basin_cpu basin_obj(chi_networks.basins[bi].outlet_flat, flow);
                lsdtt_basin_outline_cpu outline(basin_obj);
                const lsdtt_vector_record rec = outline.get_outline_as_polygon_record(bk);

                if (!std::holds_alternative<std::monostate>(rec.geometry))
                {
                    const auto& poly = std::get<lsdtt_polygon_geom>(rec.geometry);
                    if (poly.exterior.size() >= 4)
                    {
                        pw.write_polygon(
                            poly.exterior,
                            {static_cast<std::int64_t>(bk)}
                        );
                    }
                }
            }
            catch (const std::exception& e)
            {
                std::cerr << "Warning (write_chi_data_maps): polygon for basin_key="
                          << bk << " failed: " << e.what() << "\n";
            }
            ++bk;
        }

        pw.finish_vector_write();
    }
    catch (const std::exception& e)
    {
        std::cerr << "Warning (write_chi_data_maps): polygon file failed: "
                  << e.what() << "\n";
    }

    std::cout << "  Chi data maps done.\n\n";
    std::cout.flush();
}


void write_chi_data_map_lines(
    const lsdtt_chi_collection_result& chi_networks,
    const lsdtt_flow_base&             flow,
    const lsdtt_raster_cpu&            dem,
    const std::string&                 output_prefix,
    vector_output_format               fmt,
    bool                               spatial_index)
{
    const lsdtt_raster_meta& meta      = dem.get_meta();
    const int                n_cols    = meta.n_cols;
    const int                n_rows    = meta.n_rows;
    const double             x_min     = meta.x_min;
    const double             y_min     = meta.y_min;
    const double             cell_size = meta.cell_size;

    auto flat_to_xy = [&](std::int64_t flat) -> lsdtt_xy
    {
        const int row = static_cast<int>(flat / static_cast<std::int64_t>(n_cols));
        const int col = static_cast<int>(flat % static_cast<std::int64_t>(n_cols));
        return {x_min + (col + 0.5) * cell_size,
                y_min + (n_rows - row - 0.5) * cell_size};
    };

    lsdtt_vector_schema schema;
    schema.field_names = {"source_key", "basin_key"};
    schema.field_types = {
        lsdtt_vector_field_type::int64_type,
        lsdtt_vector_field_type::int64_type
    };

    const std::string ext   = output_extension(fmt);
    const std::string fname = output_prefix + "_chi_data_map_lines" + ext;

    std::cout << "  Writing chi data map lines: " << fname << "\n";
    std::cout.flush();

    lsdtt_data_loader writer;
    writer.begin_vector_write(
        fname, "chi_data_map_lines", schema, meta.crs_wkt, fmt, spatial_index);

    const std::vector<lsdtt_chi_basin_source_keys> basin_keys =
        assign_chi_source_keys(chi_networks, flow, n_cols, cell_size, false);

    for (std::size_t bi = 0; bi < chi_networks.basins.size(); ++bi)
    {
        const lsdtt_chi_basin_source_keys& keys = basin_keys[bi];
        if (keys.basin_key < 0) continue;
        const auto&        basin     = chi_networks.basins[bi];
        const std::int32_t basin_key = keys.basin_key;
        const std::size_t  n         = basin.flats.size();
        const std::vector<std::int32_t>& src_key_of = keys.source_key;

        // Collect node indices per source_key in natural (outlet-to-source) order.
        std::unordered_map<std::int32_t, std::vector<std::size_t>> nodes_by_key;
        nodes_by_key.reserve(n / 4 + 1);
        for (std::size_t i = 0; i < n; ++i)
        {
            if (src_key_of[i] >= 0)
                nodes_by_key[src_key_of[i]].push_back(i);
        }

        // Write one polyline per source_key.
        for (const auto& [sk, indices] : nodes_by_key)
        {
            if (indices.size() < 2) continue;

            std::vector<lsdtt_xy> pts;
            pts.reserve(indices.size());
            for (std::size_t idx : indices)
                pts.push_back(flat_to_xy(basin.flats[idx]));

            writer.write_line(pts, {
                static_cast<std::int64_t>(sk),
                static_cast<std::int64_t>(basin_key)
            });
        }
    }

    writer.finish_vector_write();

    std::cout << "  Chi data map lines done.\n";
    std::cout.flush();
}

} // namespace lsdtt
