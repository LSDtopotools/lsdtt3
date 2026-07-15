// lsdtt_chi_data_maps_cpu.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html
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

    std::int32_t next_basin_key  = 0;
    std::int32_t next_source_key = 0;

    // Pre-filter: keep only non-empty complete basins.
    std::vector<std::size_t> valid_basin_indices;
    valid_basin_indices.reserve(chi_networks.basins.size());

    for (std::size_t bi = 0; bi < chi_networks.basins.size(); ++bi)
    {
        const auto& basin = chi_networks.basins[bi];

        if (basin.flats.empty())
        {
            std::cout << "  basin[" << bi << "] has no nodes; skipping.\n";
            continue;
        }

        lsdtt_basin_cpu basin_obj(basin.outlet_flat, flow);
        if (!basin_obj.check_basin_completeness())
        {
            std::cout << "  basin[" << bi << "] outlet_flat=" << basin.outlet_flat
                      << " is incomplete (upslope area touches nodata); skipping.\n";
            continue;
        }

        std::cout << "  basin[" << bi << "]  outlet_flat=" << basin.outlet_flat
                  << "  n_nodes=" << basin.flats.size() << "  [complete]\n";
        valid_basin_indices.push_back(bi);
    }

    std::cout << "  Complete basins: " << valid_basin_indices.size()
              << " / " << chi_networks.basins.size() << "\n";
    std::cout.flush();

    const std::vector<float>& dem_data = dem.get_data();

    for (const std::size_t bi : valid_basin_indices)
    {
        const auto&        basin     = chi_networks.basins[bi];
        const std::int32_t basin_key = next_basin_key++;
        const std::size_t  n         = basin.flats.size();

        // Local index map (flat -> position in basin vectors).
        std::unordered_map<std::int64_t, std::size_t> local_idx;
        local_idx.reserve(n);
        for (std::size_t i = 0; i < n; ++i)
            local_idx[basin.flats[i]] = i;

        // Donor lists: needed for flow distance, max-path-up, and source keys.
        // The receiver-before-donor ordering means flats[0] is the outlet and
        // each node with index > 0 has its receiver at a lower index.
        std::vector<std::vector<std::size_t>> donors(n);
        for (std::size_t i = 1; i < n; ++i)
        {
            const std::size_t rec_idx = local_idx.at(basin.receiver_flats[i]);
            donors[rec_idx].push_back(i);
        }

        // Flow distance from outlet (receiver-before-donor: receiver's distance
        // is set before we encounter any of its donors).
        std::vector<double> flow_dist(n, 0.0);
        for (std::size_t i = 1; i < n; ++i)
        {
            const std::int64_t flat    = basin.flats[i];
            const std::int64_t rec     = basin.receiver_flats[i];
            const std::size_t  rec_idx = local_idx.at(rec);
            flow_dist[i] = flow_dist[rec_idx] +
                           d8_step_length(flat, rec, n_cols, cell_size);
        }

        // Primary-donor (longest upstream path) for source-key inheritance.
        // Processing in reverse index order (upstream-first) ensures all donors
        // of node ii are already processed before ii itself.
        constexpr std::size_t k_no_donor = std::numeric_limits<std::size_t>::max();
        std::vector<double>      max_path_up(n, 0.0);
        std::vector<std::size_t> primary_donor(n, k_no_donor);

        for (std::size_t ii = n; ii-- > 0; )
        {
            if (donors[ii].empty())
                continue;

            double      best   = -1.0;
            std::size_t best_d = k_no_donor;

            for (std::size_t d : donors[ii])
            {
                const double step =
                    d8_step_length(basin.flats[d], basin.flats[ii], n_cols, cell_size);
                const double candidate = step + max_path_up[d];
                if (best_d == k_no_donor || candidate > best)
                {
                    best   = candidate;
                    best_d = d;
                }
            }

            max_path_up[ii]   = best;
            primary_donor[ii] = best_d;
        }

        // Source-key assignment.
        // Headwaters (no donors) receive a new global key.
        // Non-headwaters inherit from their primary donor so the longest
        // tributary carries the same key continuously to the outlet.
        std::vector<std::int32_t> src_key_of(n, -1);
        std::vector<bool>         is_source(n, false);

        for (std::size_t ii = n; ii-- > 0; )
        {
            if (donors[ii].empty())
            {
                is_source[ii]  = true;
                src_key_of[ii] = next_source_key++;
            }
            else
            {
                src_key_of[ii] = src_key_of[primary_donor[ii]];
            }
        }

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

            const double da =
                static_cast<double>(flow.s_vector_size_unchecked(flat)) *
                cell_size * cell_size;

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
        for (const std::size_t bi : valid_basin_indices)
        {
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

    std::int32_t next_basin_key  = 0;
    std::int32_t next_source_key = 0;

    for (std::size_t bi = 0; bi < chi_networks.basins.size(); ++bi)
    {
        const auto& basin = chi_networks.basins[bi];

        if (basin.flats.empty()) continue;

        lsdtt_basin_cpu basin_obj(basin.outlet_flat, flow);
        if (!basin_obj.check_basin_completeness()) continue;

        const std::int32_t basin_key = next_basin_key++;
        const std::size_t  n         = basin.flats.size();

        // Build local index map and donor lists (same as write_chi_data_maps).
        std::unordered_map<std::int64_t, std::size_t> local_idx;
        local_idx.reserve(n);
        for (std::size_t i = 0; i < n; ++i)
            local_idx[basin.flats[i]] = i;

        std::vector<std::vector<std::size_t>> donors(n);
        for (std::size_t i = 1; i < n; ++i)
        {
            const std::size_t rec_idx = local_idx.at(basin.receiver_flats[i]);
            donors[rec_idx].push_back(i);
        }

        // Primary donor (longest upstream path).
        constexpr std::size_t k_no_donor = std::numeric_limits<std::size_t>::max();
        std::vector<double>      max_path_up(n, 0.0);
        std::vector<std::size_t> primary_donor(n, k_no_donor);

        for (std::size_t ii = n; ii-- > 0; )
        {
            if (donors[ii].empty()) continue;

            double best = -1.0;
            std::size_t best_d = k_no_donor;

            for (std::size_t d : donors[ii])
            {
                const double step =
                    d8_step_length(basin.flats[d], basin.flats[ii], n_cols, cell_size);
                const double candidate = step + max_path_up[d];
                if (best_d == k_no_donor || candidate > best)
                {
                    best   = candidate;
                    best_d = d;
                }
            }
            max_path_up[ii]   = best;
            primary_donor[ii] = best_d;
        }

        // Source-key assignment (same convention as write_chi_data_maps).
        std::vector<std::int32_t> src_key_of(n, -1);
        for (std::size_t ii = n; ii-- > 0; )
        {
            if (donors[ii].empty())
                src_key_of[ii] = next_source_key++;
            else
                src_key_of[ii] = src_key_of[primary_donor[ii]];
        }

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
