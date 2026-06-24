// lsdtt_vector_to_junction_network.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html

#include "lsdtt_vector_to_junction_network.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <map>
#include <memory>
#include <queue>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <gdal_priv.h>
#include <ogrsf_frmts.h>

#include "../core/lsdtt_vector_utils.hpp"
#include "../utils/lsdtt_raster_geometry.hpp"

namespace lsdtt
{
namespace
{

struct gdal_dataset_closer
{
    void operator()(GDALDataset* ds) const
    {
        if (ds) GDALClose(ds);
    }
};

struct oct_closer
{
    void operator()(OGRCoordinateTransformation* tx) const
    {
        if (tx) OCTDestroyCoordinateTransformation(tx);
    }
};

struct raw_polyline
{
    std::vector<lsdtt_xy> xy;
};

struct raw_vertex
{
    double x = 0.0;
    double y = 0.0;
    double elevation = 0.0;
    bool is_dem_edge = false;

    std::vector<std::int32_t> incident_raw_edges;
};

struct raw_edge
{
    std::int32_t raw_edge_id = -1;
    std::int32_t start_vertex = -1;
    std::int32_t end_vertex   = -1;

    double length_metres = 0.0;
    std::vector<lsdtt_xy> xy;
};

struct compact_edge_builder
{
    std::int32_t junction_a = -1;
    std::int32_t junction_b = -1;
    double length_metres = 0.0;
    std::vector<lsdtt_xy> xy;
};

struct quantised_key
{
    long long qx = 0;
    long long qy = 0;

    bool operator<(const quantised_key& other) const noexcept
    {
        if (qx != other.qx) return qx < other.qx;
        return qy < other.qy;
    }
};

OGRPolygon make_bbox_polygon(
    const lsdtt_raster_meta& meta,
    double inset_metres)
{
    OGRPolygon poly;
    OGRLinearRing ring;

    const double x0 = meta.x_min + inset_metres;
    const double x1 = meta.x_max() - inset_metres;
    const double y0 = meta.y_min + inset_metres;
    const double y1 = meta.y_max() - inset_metres;

    if (!(x1 > x0 && y1 > y0))
    {
        throw std::runtime_error(
            "lsdtt_vector_to_junction_network: DEM inset clipping polygon is invalid"
        );
    }

    ring.addPoint(x0, y0);
    ring.addPoint(x1, y0);
    ring.addPoint(x1, y1);
    ring.addPoint(x0, y1);
    ring.addPoint(x0, y0);

    poly.addRing(&ring);
    return poly;
}


double resolved_vertex_tolerance(
    const lsdtt_vector_to_junction_network_options& options,
    const lsdtt_raster_meta& meta)
{
    if (options.vertex_snap_tolerance_metres > 0.0)
    {
        return options.vertex_snap_tolerance_metres;
    }

    return 0.25 * meta.cell_size;
}

double resolved_edge_tolerance(
    const lsdtt_vector_to_junction_network_options& options,
    const lsdtt_raster_meta& meta)
{
    if (options.dem_edge_tolerance_metres > 0.0)
    {
        return options.dem_edge_tolerance_metres;
    }

    return 0.25 * meta.cell_size;
}

double resolved_minimum_clipped_edge_length(
    const lsdtt_vector_to_junction_network_options& options,
    const lsdtt_raster_meta& meta)
{
    if (options.minimum_clipped_edge_length_metres > 0.0)
    {
        return options.minimum_clipped_edge_length_metres;
    }

    return 2.0 * meta.cell_size;
}

double resolved_both_edge_min_penetration(
    const lsdtt_vector_to_junction_network_options& options,
    const lsdtt_raster_meta& meta)
{
    if (options.both_edge_min_penetration_metres > 0.0)
    {
        return options.both_edge_min_penetration_metres;
    }

    return 100.0 * meta.cell_size;
}

double resolved_clip_inset_distance(
    const lsdtt_raster_meta& meta)
{
    // Push clipped endpoints slightly inside the DEM so boundary endpoints
    // never lie just outside the raster because of exact bbox clipping or
    // floating-point roundoff.
    return 0.25 * meta.cell_size;
}



quantised_key make_key(double x, double y, double tol)
{
    quantised_key k;
    k.qx = static_cast<long long>(std::llround(x / tol));
    k.qy = static_cast<long long>(std::llround(y / tol));
    return k;
}

bool point_is_on_dem_edge(
    const lsdtt_raster_meta& meta,
    double x,
    double y,
    double tol)
{
    const bool in_x =
        x >= meta.x_min - tol && x <= meta.x_max() + tol;

    const bool in_y =
        y >= meta.y_min - tol && y <= meta.y_max() + tol;

    if (!in_x || !in_y)
    {
        return false;
    }

    return std::abs(x - meta.x_min)   <= tol ||
           std::abs(x - meta.x_max()) <= tol ||
           std::abs(y - meta.y_min)   <= tol ||
           std::abs(y - meta.y_max()) <= tol;
}

std::vector<lsdtt_xy> deduplicate_consecutive_points(
    const std::vector<lsdtt_xy>& in,
    double tol)
{
    std::vector<lsdtt_xy> out;
    out.reserve(in.size());

    for (const auto& p : in)
    {
        if (out.empty())
        {
            out.push_back(p);
            continue;
        }

        const double dx = p.x - out.back().x;
        const double dy = p.y - out.back().y;

        if (std::sqrt(dx * dx + dy * dy) > tol)
        {
            out.push_back(p);
        }
    }

    return out;
}

double polyline_length(const std::vector<lsdtt_xy>& xy)
{
    double length = 0.0;

    for (std::size_t i = 1; i < xy.size(); ++i)
    {
        const double dx = xy[i].x - xy[i - 1].x;
        const double dy = xy[i].y - xy[i - 1].y;
        length += std::sqrt(dx * dx + dy * dy);
    }

    return length;
}

// Maximum interior penetration of a polyline: the furthest any vertex gets from
// the nearest DEM boundary. A channel that hugs the DEM edge (meandering in and
// out) has a long path length but a small penetration; a genuine through- or
// diagonal channel penetrates deep. This discriminates the two, where
// polyline_length alone cannot.
double polyline_interior_penetration(
    const std::vector<lsdtt_xy>& xy,
    const lsdtt_raster_meta& meta)
{
    const double x_max = meta.x_max();
    const double y_max = meta.y_max();

    double max_penetration = 0.0;

    for (const auto& p : xy)
    {
        const double dist_to_boundary = std::min(
            std::min(p.x - meta.x_min, x_max - p.x),
            std::min(p.y - meta.y_min, y_max - p.y));

        if (dist_to_boundary > max_penetration)
        {
            max_penetration = dist_to_boundary;
        }
    }

    return max_penetration;
}

void collect_clipped_linestrings(
    OGRGeometry* geom,
    double point_tol,
    std::vector<raw_polyline>& out)
{
    if (geom == nullptr || geom->IsEmpty())
    {
        return;
    }

    const OGRwkbGeometryType gt = wkbFlatten(geom->getGeometryType());

    if (gt == wkbLineString)
    {
        OGRLineString* ls = static_cast<OGRLineString*>(geom);

        std::vector<lsdtt_xy> xy;
        xy.reserve(static_cast<std::size_t>(ls->getNumPoints()));

        for (int i = 0; i < ls->getNumPoints(); ++i)
        {
            xy.push_back({ls->getX(i), ls->getY(i)});
        }

        xy = deduplicate_consecutive_points(xy, point_tol);

        if (xy.size() >= 2 && polyline_length(xy) > 0.0)
        {
            raw_polyline rp;
            rp.xy = std::move(xy);
            out.push_back(std::move(rp));
        }

        return;
    }

    if (gt == wkbMultiLineString || gt == wkbGeometryCollection)
    {
        OGRGeometryCollection* gc =
            static_cast<OGRGeometryCollection*>(geom);

        for (int i = 0; i < gc->getNumGeometries(); ++i)
        {
            collect_clipped_linestrings(
                gc->getGeometryRef(i),
                point_tol,
                out
            );
        }
    }
}

double sample_dem_elevation_nearest_valid(
    const lsdtt_raster_cpu& dem,
    double x,
    double y)
{
    const lsdtt_raster_meta& meta = dem.get_meta();
    const std::vector<float>& data = dem.get_data();
    const std::vector<uint8_t>& mask = dem.get_mask();

    int col = static_cast<int>(
        std::floor((x - meta.x_min) / meta.cell_size)
    );

    int row = static_cast<int>(
        std::floor((meta.y_max() - y) / meta.cell_size)
    );

    if (col == meta.n_cols) col = meta.n_cols - 1;
    if (row == meta.n_rows) row = meta.n_rows - 1;

    if (row < 0 || row >= meta.n_rows || col < 0 || col >= meta.n_cols)
    {
        throw std::runtime_error(
            "lsdtt_vector_to_junction_network: attempted to sample DEM elevation outside raster extent"
        );
    }

    auto idx_of = [&](int r, int c) -> std::size_t
    {
        return static_cast<std::size_t>(r) *
                   static_cast<std::size_t>(meta.n_cols) +
               static_cast<std::size_t>(c);
    };

    for (int radius = 0; radius <= 1; ++radius)
    {
        for (int rr = row - radius; rr <= row + radius; ++rr)
        {
            for (int cc = col - radius; cc <= col + radius; ++cc)
            {
                if (rr < 0 || rr >= meta.n_rows || cc < 0 || cc >= meta.n_cols)
                {
                    continue;
                }

                const std::size_t idx = idx_of(rr, cc);
                if (mask[idx] == 0) continue;
                if (meta.is_nodata(data[idx])) continue;
                if (!std::isfinite(data[idx])) continue;

                return static_cast<double>(data[idx]);
            }
        }
    }

    throw std::runtime_error(
        "lsdtt_vector_to_junction_network: no valid DEM elevation found beneath or adjacent to vector vertex"
    );
}

std::vector<raw_polyline> load_clipped_stream_polylines(
    const std::string&      vector_path,
    const lsdtt_raster_cpu& dem,
    double                  point_tol,
    double                  clip_inset_metres)
{
    GDALAllRegister();

    std::unique_ptr<GDALDataset, gdal_dataset_closer> ds(
        static_cast<GDALDataset*>(
            GDALOpenEx(
                vector_path.c_str(),
                GDAL_OF_VECTOR | GDAL_OF_READONLY,
                nullptr, nullptr, nullptr
            )
        )
    );

    if (!ds)
    {
        throw std::runtime_error(
            "lsdtt_vector_to_junction_network: could not open vector dataset: " +
            vector_path
        );
    }

    const lsdtt_raster_meta& meta = dem.get_meta();
    OGRPolygon bbox = make_bbox_polygon(meta, clip_inset_metres);


    OGRSpatialReference raster_srs;
    OGRSpatialReference* raster_srs_ptr = nullptr;

    if (!meta.crs_wkt.empty())
    {
        raster_srs.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);
        if (raster_srs.importFromWkt(meta.crs_wkt.c_str()) == OGRERR_NONE)
        {
            raster_srs_ptr = &raster_srs;
        }
    }

    std::vector<raw_polyline> polylines;

    for (int li = 0; li < ds->GetLayerCount(); ++li)
    {
        OGRLayer* layer = ds->GetLayer(li);
        if (!layer) continue;

        std::unique_ptr<OGRCoordinateTransformation, oct_closer> tx;

        const OGRSpatialReference* layer_srs = layer->GetSpatialRef();
        if (layer_srs && raster_srs_ptr && !layer_srs->IsSame(raster_srs_ptr))
        {
            OGRSpatialReference layer_srs_copy(*layer_srs);
            layer_srs_copy.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);

            tx.reset(
                OGRCreateCoordinateTransformation(&layer_srs_copy, raster_srs_ptr)
            );

            if (!tx)
            {
                throw std::runtime_error(
                    "lsdtt_vector_to_junction_network: failed to create vector->DEM CRS transformation"
                );
            }
        }

        layer->ResetReading();

        OGRFeature* feature = nullptr;
        while ((feature = layer->GetNextFeature()) != nullptr)
        {
            OGRGeometry* geom = feature->GetGeometryRef();

            if (!geom)
            {
                OGRFeature::DestroyFeature(feature);
                continue;
            }

            std::unique_ptr<OGRGeometry> geom_work(geom->clone());

            if (tx)
            {
                if (geom_work->transform(tx.get()) != OGRERR_NONE)
                {
                    OGRFeature::DestroyFeature(feature);
                    continue;
                }
            }

            std::unique_ptr<OGRGeometry> clipped(
                geom_work->Intersection(&bbox)
            );

            if (clipped)
            {
                collect_clipped_linestrings(
                    clipped.get(),
                    point_tol,
                    polylines
                );
            }

            OGRFeature::DestroyFeature(feature);
        }
    }

    return polylines;
}

// Drop clipped fragments that are artefacts of the DEM boundary:
//   * short stubs that touch the edge and carry only a tiny in-DEM length;
//   * "pass-in/pass-out" meanders whose BOTH endpoints sit on the edge and whose
//     interior penetration is below threshold (a channel hugging the boundary,
//     dipping in and out, vs a genuine through/diagonal channel that penetrates
//     deep). Penetration discriminates these where path length cannot.
std::vector<raw_polyline> filter_boundary_polylines(
    std::vector<raw_polyline>  polylines,
    const lsdtt_raster_meta&   meta,
    double                     edge_tol,
    double                     min_clipped_edge_length,
    double                     both_edge_min_penetration)
{
    std::vector<raw_polyline> filtered;
    filtered.reserve(polylines.size());

    for (auto& poly : polylines)
    {
        const double len = polyline_length(poly.xy);

        const bool start_on_edge = point_is_on_dem_edge(
            meta, poly.xy.front().x, poly.xy.front().y, edge_tol);
        const bool end_on_edge = point_is_on_dem_edge(
            meta, poly.xy.back().x, poly.xy.back().y, edge_tol);

        const bool touches_dem_edge = start_on_edge || end_on_edge;

        if (touches_dem_edge && len < min_clipped_edge_length)
        {
            continue;
        }

        const bool both_on_edge = start_on_edge && end_on_edge;

        if (both_on_edge &&
            polyline_interior_penetration(poly.xy, meta) <
                both_edge_min_penetration)
        {
            continue;
        }

        filtered.push_back(std::move(poly));
    }

    return filtered;
}

std::string join_ids_csv(const std::vector<std::int32_t>& ids)
{
    std::ostringstream ss;

    for (std::size_t i = 0; i < ids.size(); ++i)
    {
        if (i > 0) ss << ",";
        ss << ids[i];
    }

    return ss.str();
}

std::vector<lsdtt_xy> oriented_geometry(
    const raw_edge& edge,
    std::int32_t from_vertex)
{
    if (from_vertex == edge.start_vertex)
    {
        return edge.xy;
    }

    if (from_vertex == edge.end_vertex)
    {
        std::vector<lsdtt_xy> rev = edge.xy;
        std::reverse(rev.begin(), rev.end());
        return rev;
    }

    throw std::runtime_error(
        "lsdtt_vector_to_junction_network: requested raw-edge orientation from non-incident vertex"
    );
}

} // unnamed namespace


// -----------------------------------------------------------------------------
// clip_vector_lines_to_dem
//
// Public helper: clip a vector line file to the DEM extent and return the in-DEM
// linework, reusing the exact same logic the junction-network builder uses
// (inset-bbox OGR intersection + boundary-stub and pass-in/pass-out meander
// removal). Lines that meander in and out of the raster along the boundary are
// split by the clip into separate in-DEM pieces and the shallow boundary
// meanders are dropped, so callers (e.g. the VML stream burn) operate only on
// genuine in-DEM channels.
// -----------------------------------------------------------------------------
std::vector<std::vector<lsdtt_xy>> clip_vector_lines_to_dem(
    const std::string&                              vector_path,
    const lsdtt_raster_cpu&                         dem,
    const lsdtt_vector_to_junction_network_options& options)
{
    const lsdtt_raster_meta& meta = dem.get_meta();

    const double vertex_tol = resolved_vertex_tolerance(options, meta);
    const double edge_tol   = resolved_edge_tolerance(options, meta);
    const double min_clipped_edge_length =
        resolved_minimum_clipped_edge_length(options, meta);
    const double both_edge_min_penetration =
        resolved_both_edge_min_penetration(options, meta);
    const double clip_inset = resolved_clip_inset_distance(meta);

    std::vector<raw_polyline> polylines = filter_boundary_polylines(
        load_clipped_stream_polylines(vector_path, dem, vertex_tol, clip_inset),
        meta, edge_tol, min_clipped_edge_length, both_edge_min_penetration);

    std::vector<std::vector<lsdtt_xy>> out;
    out.reserve(polylines.size());
    for (auto& p : polylines) out.push_back(std::move(p.xy));
    return out;
}


// -----------------------------------------------------------------------------
// construction
// -----------------------------------------------------------------------------
lsdtt_vector_to_junction_network::lsdtt_vector_to_junction_network(
    const std::string&                              vector_path,
    const lsdtt_raster_cpu&                         dem,
    const lsdtt_vector_to_junction_network_options& options)
{
    create(vector_path, dem, options);
}

// -----------------------------------------------------------------------------
// clear
// -----------------------------------------------------------------------------
void lsdtt_vector_to_junction_network::clear() noexcept
{
    meta_ = lsdtt_raster_meta{};

    is_built_ = false;
    has_topological_order_ = false;
    has_s_vector_ = false;

    junctions_.clear();
    edges_.clear();
    topological_order_.clear();
    outlets_.clear();
    s_vector_.clear();
}

// -----------------------------------------------------------------------------
// create
// -----------------------------------------------------------------------------
void lsdtt_vector_to_junction_network::create(
    const std::string&                              vector_path,
    const lsdtt_raster_cpu&                         dem,
    const lsdtt_vector_to_junction_network_options& options)
{
    clear();

    const lsdtt_raster_meta& meta = dem.get_meta();
    meta.check_projected_metres();

    meta_ = meta;

    const double vertex_tol = resolved_vertex_tolerance(options, meta);
    const double edge_tol   = resolved_edge_tolerance(options, meta);
    const double min_clipped_edge_length =
        resolved_minimum_clipped_edge_length(options, meta);
    const double both_edge_min_penetration =
        resolved_both_edge_min_penetration(options, meta);
    const double clip_inset_metres =
        resolved_clip_inset_distance(meta);

    std::vector<raw_polyline> polylines =
        load_clipped_stream_polylines(
            vector_path,
            dem,
            vertex_tol,
            clip_inset_metres
        );


    if (polylines.empty())
    {
        throw std::runtime_error(
            "lsdtt_vector_to_junction_network: no clipped stream polylines intersect the DEM extent"
        );
    }

    // -------------------------------------------------------------------------
    // Remove boundary stubs and pass-in/pass-out meanders (see
    // filter_boundary_polylines).
    // -------------------------------------------------------------------------
    polylines = filter_boundary_polylines(
        std::move(polylines), meta, edge_tol,
        min_clipped_edge_length, both_edge_min_penetration);

    if (polylines.empty())
    {
        throw std::runtime_error(
            "lsdtt_vector_to_junction_network: all clipped stream polylines were removed as short boundary stubs"
        );
    }



    // -------------------------------------------------------------------------
    // Build raw endpoint graph from clipped linework.
    // -------------------------------------------------------------------------
    std::map<quantised_key, std::int32_t> vertex_lookup;
    std::vector<raw_vertex> raw_vertices;
    std::vector<raw_edge>   raw_edges;

    auto get_or_make_vertex =
        [&](const lsdtt_xy& p) -> std::int32_t
        {
            const quantised_key key = make_key(p.x, p.y, vertex_tol);
            auto it = vertex_lookup.find(key);

            if (it != vertex_lookup.end())
            {
                return it->second;
            }

            raw_vertex v;
            v.x = p.x;
            v.y = p.y;
            v.elevation = sample_dem_elevation_nearest_valid(dem, p.x, p.y);
            v.is_dem_edge = point_is_on_dem_edge(meta, p.x, p.y, edge_tol);

            const std::int32_t id =
                static_cast<std::int32_t>(raw_vertices.size());

            raw_vertices.push_back(std::move(v));
            vertex_lookup[key] = id;
            return id;
        };

    for (const auto& poly : polylines)
    {
        const std::int32_t v0 = get_or_make_vertex(poly.xy.front());
        const std::int32_t v1 = get_or_make_vertex(poly.xy.back());

        raw_edge e;
        e.raw_edge_id = static_cast<std::int32_t>(raw_edges.size());
        e.start_vertex = v0;
        e.end_vertex   = v1;
        e.xy = poly.xy;
        e.length_metres = polyline_length(e.xy);

        raw_edges.push_back(e);
        raw_vertices[static_cast<std::size_t>(v0)].incident_raw_edges.push_back(e.raw_edge_id);
        raw_vertices[static_cast<std::size_t>(v1)].incident_raw_edges.push_back(e.raw_edge_id);
    }

    if (raw_edges.empty())
    {
        throw std::runtime_error(
            "lsdtt_vector_to_junction_network: no valid clipped raw edges were produced"
        );
    }

    // -------------------------------------------------------------------------
    // Mark raw junction vertices:
    //   - degree != 2
    //   - any DEM-edge vertex
    // -------------------------------------------------------------------------
    std::vector<uint8_t> is_raw_junction(
        raw_vertices.size(),
        uint8_t{0}
    );

    for (std::size_t i = 0; i < raw_vertices.size(); ++i)
    {
        const int degree = static_cast<int>(
            raw_vertices[i].incident_raw_edges.size()
        );

        if (degree != 2 || raw_vertices[i].is_dem_edge)
        {
            is_raw_junction[i] = uint8_t{1};
        }
    }

    // -------------------------------------------------------------------------
    // If a component contains no junctions at all, break it at its lowest
    // elevation vertex so that the compact graph remains traversable.
    // -------------------------------------------------------------------------
    {
        std::vector<int> component_id(raw_vertices.size(), -1);
        int next_component = 0;

        for (std::size_t seed = 0; seed < raw_vertices.size(); ++seed)
        {
            if (component_id[seed] >= 0) continue;

            std::queue<std::int32_t> q;
            q.push(static_cast<std::int32_t>(seed));
            component_id[seed] = next_component;

            while (!q.empty())
            {
                const std::int32_t v = q.front();
                q.pop();

                for (const std::int32_t eidx :
                     raw_vertices[static_cast<std::size_t>(v)].incident_raw_edges)
                {
                    const raw_edge& e = raw_edges[static_cast<std::size_t>(eidx)];
                    const std::int32_t w =
                        (e.start_vertex == v) ? e.end_vertex : e.start_vertex;

                    if (component_id[static_cast<std::size_t>(w)] >= 0) continue;
                    component_id[static_cast<std::size_t>(w)] = next_component;
                    q.push(w);
                }
            }

            ++next_component;
        }

        for (int cid = 0; cid < next_component; ++cid)
        {
            bool has_junction = false;
            std::int32_t lowest_vertex = -1;
            double lowest_elev = std::numeric_limits<double>::infinity();

            for (std::size_t i = 0; i < raw_vertices.size(); ++i)
            {
                if (component_id[i] != cid) continue;

                if (is_raw_junction[i] != 0)
                {
                    has_junction = true;
                }

                if (raw_vertices[i].elevation < lowest_elev)
                {
                    lowest_elev = raw_vertices[i].elevation;
                    lowest_vertex = static_cast<std::int32_t>(i);
                }
            }

            if (!has_junction && lowest_vertex >= 0)
            {
                is_raw_junction[static_cast<std::size_t>(lowest_vertex)] = uint8_t{1};
            }
        }
    }

    // -------------------------------------------------------------------------
    // Build compact junction list.
    // -------------------------------------------------------------------------
    std::vector<std::int32_t> raw_to_junction(
        raw_vertices.size(),
        k_nodata_junction
    );

    for (std::size_t rv = 0; rv < raw_vertices.size(); ++rv)
    {
        if (is_raw_junction[rv] == 0) continue;

        lsdtt_vector_network_junction j;
        j.junction_id = static_cast<std::int32_t>(junctions_.size());
        j.x = raw_vertices[rv].x;
        j.y = raw_vertices[rv].y;
        j.elevation = raw_vertices[rv].elevation;
        j.is_dem_edge_junction = raw_vertices[rv].is_dem_edge;

        raw_to_junction[rv] = j.junction_id;
        junctions_.push_back(std::move(j));
    }

    // -------------------------------------------------------------------------
    // Collapse degree-2 chains into compact edges.
    // -------------------------------------------------------------------------
    std::vector<uint8_t> raw_edge_visited(
        raw_edges.size(),
        uint8_t{0}
    );

    for (std::size_t rv = 0; rv < raw_vertices.size(); ++rv)
    {
        if (is_raw_junction[rv] == 0) continue;

        const std::int32_t j0 = raw_to_junction[rv];

        for (const std::int32_t seed_edge_id :
             raw_vertices[rv].incident_raw_edges)
        {
            if (raw_edge_visited[static_cast<std::size_t>(seed_edge_id)] != 0)
            {
                continue;
            }

            compact_edge_builder builder;
            builder.junction_a = j0;

            std::int32_t current_vertex = static_cast<std::int32_t>(rv);
            std::int32_t current_edge_id = seed_edge_id;

            while (true)
            {
                raw_edge_visited[static_cast<std::size_t>(current_edge_id)] = uint8_t{1};

                const raw_edge& current_edge =
                    raw_edges[static_cast<std::size_t>(current_edge_id)];

                std::vector<lsdtt_xy> geom =
                    oriented_geometry(current_edge, current_vertex);

                if (builder.xy.empty())
                {
                    builder.xy = geom;
                }
                else
                {
                    builder.xy.insert(
                        builder.xy.end(),
                        geom.begin() + 1,
                        geom.end()
                    );
                }

                builder.length_metres += current_edge.length_metres;

                const std::int32_t next_vertex =
                    (current_edge.start_vertex == current_vertex)
                        ? current_edge.end_vertex
                        : current_edge.start_vertex;

                if (is_raw_junction[static_cast<std::size_t>(next_vertex)] != 0)
                {
                    builder.junction_b =
                        raw_to_junction[static_cast<std::size_t>(next_vertex)];
                    break;
                }

                const auto& incident =
                    raw_vertices[static_cast<std::size_t>(next_vertex)].incident_raw_edges;

                if (incident.size() != 2)
                {
                    throw std::runtime_error(
                        "lsdtt_vector_to_junction_network: encountered non-junction interior vertex with degree != 2"
                    );
                }

                const std::int32_t next_edge_id =
                    (incident[0] == current_edge_id) ? incident[1] : incident[0];

                current_vertex = next_vertex;
                current_edge_id = next_edge_id;
            }

            lsdtt_vector_network_edge edge;
            edge.edge_id = static_cast<std::int32_t>(edges_.size());
            edge.junction_a = builder.junction_a;
            edge.junction_b = builder.junction_b;
            edge.length_metres = builder.length_metres;
            edge.xy = std::move(builder.xy);

            edges_.push_back(std::move(edge));
        }
    }

    if (edges_.empty())
    {
        throw std::runtime_error(
            "lsdtt_vector_to_junction_network: no compact edges were produced"
        );
    }

    // -------------------------------------------------------------------------
    // Incident-edge lists
    // -------------------------------------------------------------------------
    for (auto& j : junctions_)
    {
        j.incident_edges.clear();
    }

    for (const auto& e : edges_)
    {
        junctions_[static_cast<std::size_t>(e.junction_a)].incident_edges.push_back(e.edge_id);
        junctions_[static_cast<std::size_t>(e.junction_b)].incident_edges.push_back(e.edge_id);
    }

    // -------------------------------------------------------------------------
    // One receiver per junction, from DEM elevation.
    // -------------------------------------------------------------------------
    for (auto& j : junctions_)
    {
        double best_slope = -std::numeric_limits<double>::infinity();
        std::int32_t best_receiver = k_outlet_junction;
        std::int32_t best_edge = -1;
        int n_downhill = 0;

        for (const std::int32_t edge_id : j.incident_edges)
        {
            const auto& e = edges_[static_cast<std::size_t>(edge_id)];

            const std::int32_t other =
                (e.junction_a == j.junction_id) ? e.junction_b : e.junction_a;

            const double dz =
                j.elevation -
                junctions_[static_cast<std::size_t>(other)].elevation;

            if (dz > options.receiver_min_drop_metres)
            {
                ++n_downhill;

                const double slope =
                    dz / std::max(e.length_metres, 1.0e-12);

                if (slope > best_slope)
                {
                    best_slope = slope;
                    best_receiver = other;
                    best_edge = edge_id;
                }
            }
        }

        if (n_downhill == 0)
        {
            j.is_outlet = true;
            j.receiver_junction = k_outlet_junction;
            j.downstream_edge_id = -1;
            outlets_.push_back(j.junction_id);
        }
        else
        {
            if (options.require_unique_receiver && n_downhill > 1)
            {
                throw std::runtime_error(
                    "lsdtt_vector_to_junction_network: junction has multiple downhill receiver candidates"
                );
            }

            if (!options.choose_steepest_receiver_when_ambiguous && n_downhill > 1)
            {
                throw std::runtime_error(
                    "lsdtt_vector_to_junction_network: ambiguous downhill receiver and steepest-choice is disabled"
                );
            }

            j.is_outlet = false;
            j.receiver_junction = best_receiver;
            j.downstream_edge_id = best_edge;
        }
    }

    // -------------------------------------------------------------------------
    // Orient compact edges using receiver assignments.
    //
    // Fallback rule:
    // If neither endpoint claims the edge as its receiver edge, but one
    // endpoint has only one incident edge, then that endpoint is forced to
    // drain through this edge. This is mainly to stabilise clipped source /
    // DEM-edge junctions where DEM endpoint sampling can be slightly noisy.
    // -------------------------------------------------------------------------
    auto remove_outlet_if_present =
        [&](std::int32_t junction_id)
        {
            outlets_.erase(
                std::remove(outlets_.begin(), outlets_.end(), junction_id),
                outlets_.end()
            );
        };

    for (auto& e : edges_)
    {
        bool a_to_b =
            junctions_[static_cast<std::size_t>(e.junction_a)].receiver_junction ==
            e.junction_b;

        bool b_to_a =
            junctions_[static_cast<std::size_t>(e.junction_b)].receiver_junction ==
            e.junction_a;

        if (a_to_b && b_to_a)
        {
            throw std::runtime_error(
                "lsdtt_vector_to_junction_network: two adjacent junctions point to each other; cyclic vector receiver graph"
            );
        }

        if (!a_to_b && !b_to_a)
        {
            auto& ja = junctions_[static_cast<std::size_t>(e.junction_a)];
            auto& jb = junctions_[static_cast<std::size_t>(e.junction_b)];

            const bool a_is_leaf =
                ja.incident_edges.size() == 1;

            const bool b_is_leaf =
                jb.incident_edges.size() == 1;

            const double tol = options.receiver_min_drop_metres;

            // If one side is a degree-1 junction, force that side to drain
            // through the only available edge.
            if (a_is_leaf && !b_is_leaf)
            {
                ja.is_outlet = false;
                ja.receiver_junction = e.junction_b;
                ja.downstream_edge_id = e.edge_id;
                remove_outlet_if_present(ja.junction_id);
                a_to_b = true;
            }
            else if (b_is_leaf && !a_is_leaf)
            {
                jb.is_outlet = false;
                jb.receiver_junction = e.junction_a;
                jb.downstream_edge_id = e.edge_id;
                remove_outlet_if_present(jb.junction_id);
                b_to_a = true;
            }
            else if (a_is_leaf && b_is_leaf)
            {
                // Isolated clipped segment: force high -> low if possible.
                if (ja.elevation > jb.elevation + tol)
                {
                    ja.is_outlet = false;
                    ja.receiver_junction = e.junction_b;
                    ja.downstream_edge_id = e.edge_id;
                    jb.is_outlet = true;
                    jb.receiver_junction = k_outlet_junction;
                    jb.downstream_edge_id = -1;
                    remove_outlet_if_present(ja.junction_id);
                    if (std::find(outlets_.begin(), outlets_.end(), jb.junction_id) ==
                        outlets_.end())
                    {
                        outlets_.push_back(jb.junction_id);
                    }
                    a_to_b = true;
                }
                else if (jb.elevation > ja.elevation + tol)
                {
                    jb.is_outlet = false;
                    jb.receiver_junction = e.junction_a;
                    jb.downstream_edge_id = e.edge_id;
                    ja.is_outlet = true;
                    ja.receiver_junction = k_outlet_junction;
                    ja.downstream_edge_id = -1;
                    remove_outlet_if_present(jb.junction_id);
                    if (std::find(outlets_.begin(), outlets_.end(), ja.junction_id) ==
                        outlets_.end())
                    {
                        outlets_.push_back(ja.junction_id);
                    }
                    b_to_a = true;
                }
            }

            if (!a_to_b && !b_to_a)
            {
                if (options.break_cycles)
                {
                    // Loop-closing edge: label it and exclude it from the
                    // directed graph so the network builds despite the cycle.
                    e.is_loop_edge = true;
                    ja.is_on_loop  = true;
                    jb.is_on_loop  = true;
                    continue;
                }
                std::ostringstream ss;
                ss.precision(17);
                ss << "lsdtt_vector_to_junction_network: compact edge is not used as a "
                   << "receiver edge by either endpoint. "
                   << "edge_id=" << e.edge_id
                   << ", junction_a=" << e.junction_a
                   << ", junction_b=" << e.junction_b
                   << ", elev_a=" << ja.elevation
                   << ", elev_b=" << jb.elevation
                   << ", n_incident_a=" << ja.incident_edges.size()
                   << ", n_incident_b=" << jb.incident_edges.size()
                   << ", receiver_a=" << ja.receiver_junction
                   << ", receiver_b=" << jb.receiver_junction;
                throw std::runtime_error(ss.str());
            }
        }

        if (a_to_b)
        {
            e.upstream_junction = e.junction_a;
            e.receiver_junction = e.junction_b;
        }
        else
        {
            e.upstream_junction = e.junction_b;
            e.receiver_junction = e.junction_a;
            std::reverse(e.xy.begin(), e.xy.end());
        }
    }


    // -------------------------------------------------------------------------
    // Label redundant parallel edges as loop edges. Each junction drains through
    // exactly one edge (its downstream_edge_id); any other edge oriented OUT of
    // that junction (e.g. a second, parallel compact edge to the same receiver,
    // a bigon) is redundant and would otherwise inflate the receiver's donor
    // count above its receiver in-degree, stalling the topological sort.
    // -------------------------------------------------------------------------
    if (options.break_cycles)
    {
        for (auto& e : edges_)
        {
            if (e.is_loop_edge) continue;
            const std::int32_t up = e.upstream_junction;
            if (up < 0) continue;
            if (junctions_[static_cast<std::size_t>(up)].downstream_edge_id
                    != e.edge_id)
            {
                e.is_loop_edge = true;
                junctions_[static_cast<std::size_t>(e.junction_a)].is_on_loop = true;
                junctions_[static_cast<std::size_t>(e.junction_b)].is_on_loop = true;
            }
        }
    }

    // -------------------------------------------------------------------------
    // Break any remaining directed receiver cycles (where every edge is claimed
    // by an endpoint, so the orientation pass above could not detect them).
    // -------------------------------------------------------------------------
    if (options.break_cycles)
    {
        break_receiver_cycles();
    }

    // -------------------------------------------------------------------------
    // Donor lists
    // -------------------------------------------------------------------------
    for (auto& j : junctions_)
    {
        j.donor_junctions.clear();
    }

    for (const auto& e : edges_)
    {
        if (e.is_loop_edge) continue;   // loop edges are not in the DAG
        junctions_[static_cast<std::size_t>(e.receiver_junction)]
            .donor_junctions.push_back(e.upstream_junction);
    }

    for (auto& j : junctions_)
    {
        std::sort(j.donor_junctions.begin(), j.donor_junctions.end());
    }

    build_topological_order();
    build_outlet_keys();
    build_strahler_order();
    build_s_vector();

    is_built_ = true;

}

// -----------------------------------------------------------------------------
// require helpers
// -----------------------------------------------------------------------------
void lsdtt_vector_to_junction_network::require_built() const
{
    if (!is_built_)
    {
        throw std::runtime_error(
            "lsdtt_vector_to_junction_network: object is empty; call create() first"
        );
    }
}

void lsdtt_vector_to_junction_network::require_valid_junction(
    std::int32_t junction) const
{
    require_built();

    if (!is_valid_junction(junction))
    {
        throw std::runtime_error(
            "lsdtt_vector_to_junction_network: invalid junction index"
        );
    }
}

// -----------------------------------------------------------------------------
// accessors
// -----------------------------------------------------------------------------
std::int32_t lsdtt_vector_to_junction_network::outlet_of_junction(
    std::int32_t junction) const
{
    require_valid_junction(junction);
    return junctions_[static_cast<std::size_t>(junction)].outlet_junction;
}

std::int32_t lsdtt_vector_to_junction_network::receiver_of_junction(
    std::int32_t junction) const
{
    require_valid_junction(junction);
    return junctions_[static_cast<std::size_t>(junction)].receiver_junction;
}

const std::vector<std::int32_t>&
lsdtt_vector_to_junction_network::donor_junctions(
    std::int32_t junction) const
{
    require_valid_junction(junction);
    return junctions_[static_cast<std::size_t>(junction)].donor_junctions;
}

const std::vector<std::int32_t>&
lsdtt_vector_to_junction_network::incident_edges(
    std::int32_t junction) const
{
    require_valid_junction(junction);
    return junctions_[static_cast<std::size_t>(junction)].incident_edges;
}

std::int32_t lsdtt_vector_to_junction_network::s_vector_start(
    std::int32_t junction) const
{
    require_valid_junction(junction);

    if (!has_s_vector_)
    {
        throw std::runtime_error(
            "lsdtt_vector_to_junction_network: s_vector was not built"
        );
    }

    return junctions_[static_cast<std::size_t>(junction)].s_vector_start;
}

std::int32_t lsdtt_vector_to_junction_network::s_vector_size(
    std::int32_t junction) const
{
    require_valid_junction(junction);

    if (!has_s_vector_)
    {
        throw std::runtime_error(
            "lsdtt_vector_to_junction_network: s_vector was not built"
        );
    }

    return junctions_[static_cast<std::size_t>(junction)].s_vector_size;
}

lsdtt_junction_index_slice
lsdtt_vector_to_junction_network::upslope_junctions(
    std::int32_t junction) const
{
    require_valid_junction(junction);

    if (!has_s_vector_)
    {
        throw std::runtime_error(
            "lsdtt_vector_to_junction_network: s_vector was not built"
        );
    }

    const std::int32_t start =
        junctions_[static_cast<std::size_t>(junction)].s_vector_start;

    const std::int32_t size =
        junctions_[static_cast<std::size_t>(junction)].s_vector_size;

    return lsdtt_junction_index_slice{
        s_vector_.data() + static_cast<std::size_t>(start),
        static_cast<std::size_t>(size)
    };
}

bool lsdtt_vector_to_junction_network::is_upslope_of(
    std::int32_t candidate_junction,
    std::int32_t target_junction) const
{
    require_valid_junction(candidate_junction);
    require_valid_junction(target_junction);

    if (!has_s_vector_)
    {
        throw std::runtime_error(
            "lsdtt_vector_to_junction_network: s_vector was not built"
        );
    }

    const std::int32_t cpos =
        junctions_[static_cast<std::size_t>(candidate_junction)].s_vector_start;

    const std::int32_t tstart =
        junctions_[static_cast<std::size_t>(target_junction)].s_vector_start;

    const std::int32_t tsize =
        junctions_[static_cast<std::size_t>(target_junction)].s_vector_size;

    return cpos >= tstart && cpos < tstart + tsize;
}


void lsdtt_vector_to_junction_network::build_outlet_keys()
{
    if (!has_topological_order_)
    {
        throw std::runtime_error(
            "lsdtt_vector_to_junction_network: topological order must be built before outlet keys"
        );
    }

    for (auto& j : junctions_)
    {
        j.outlet_junction = k_nodata_junction;
    }

    // Work from downstream to upstream so each receiver already has its outlet.
    for (auto it = topological_order_.rbegin(); it != topological_order_.rend(); ++it)
    {
        const std::int32_t junction_id = *it;
        auto& j = junctions_[static_cast<std::size_t>(junction_id)];

        if (j.is_outlet)
        {
            j.outlet_junction = j.junction_id;
        }
        else
        {
            const std::int32_t rec = j.receiver_junction;

            if (rec < 0)
            {
                throw std::runtime_error(
                    "lsdtt_vector_to_junction_network: non-outlet junction has no valid receiver while building outlet keys"
                );
            }

            const std::int32_t outlet =
                junctions_[static_cast<std::size_t>(rec)].outlet_junction;

            if (outlet < 0)
            {
                throw std::runtime_error(
                    "lsdtt_vector_to_junction_network: receiver outlet key was not set before donor while building outlet keys"
                );
            }

            j.outlet_junction = outlet;
        }
    }

    for (auto& e : edges_)
    {
        if (e.is_loop_edge) { e.outlet_junction = k_nodata_junction; continue; }
        const std::int32_t up = e.upstream_junction;
        if (up < 0)
        {
            throw std::runtime_error(
                "lsdtt_vector_to_junction_network: edge has invalid upstream junction while building outlet keys"
            );
        }

        e.outlet_junction =
            junctions_[static_cast<std::size_t>(up)].outlet_junction;
    }
}


// -----------------------------------------------------------------------------
// break_receiver_cycles
//
// The receiver graph is a functional graph (each junction has one receiver, or
// is an outlet). A directed cycle a->b->...->a cannot be ordered. We detect any
// such cycle by walking receivers, and break it by demoting the highest-
// elevation junction in the cycle to an outlet: its downstream edge is labelled
// as a loop edge (so it is excluded from the DAG and flagged for the downstream
// D8 loop resolver) and its endpoints are marked is_on_loop. Repeats until no
// cycle remains. Only invoked when options.break_cycles is set.
// -----------------------------------------------------------------------------
void lsdtt_vector_to_junction_network::break_receiver_cycles()
{
    const std::size_t n = junctions_.size();

    // Iteratively: run Kahn on the receiver functional graph; the junctions that
    // never reach in-degree 0 are exactly those trapped in directed cycles.
    // Demote the highest-elevation trapped junction to an outlet (labelling its
    // downstream edge as a loop edge) and repeat until the graph is acyclic.
    for (std::size_t guard = 0; guard <= n; ++guard)
    {
        // In-degree from current receiver pointers.
        std::vector<std::int32_t> indeg(n, 0);
        for (std::size_t j = 0; j < n; ++j)
        {
            const std::int32_t rec =
                junctions_[j].receiver_junction;
            if (rec >= 0) ++indeg[static_cast<std::size_t>(rec)];
        }

        std::vector<std::int32_t> ready;
        for (std::size_t j = 0; j < n; ++j)
            if (indeg[j] == 0) ready.push_back(static_cast<std::int32_t>(j));

        std::vector<char> processed(n, 0);
        std::size_t head = 0;
        while (head < ready.size())
        {
            const std::int32_t j = ready[head++];
            processed[static_cast<std::size_t>(j)] = 1;
            const std::int32_t rec =
                junctions_[static_cast<std::size_t>(j)].receiver_junction;
            if (rec >= 0 && --indeg[static_cast<std::size_t>(rec)] == 0)
                ready.push_back(rec);
        }

        // Find the highest-elevation junction still trapped in a cycle.
        int    hi  = -1;
        double hiz = -std::numeric_limits<double>::infinity();
        for (std::size_t j = 0; j < n; ++j)
        {
            if (processed[j]) continue;
            if (junctions_[j].elevation > hiz)
            { hiz = junctions_[j].elevation; hi = static_cast<int>(j); }
        }
        if (hi < 0) return;   // acyclic: done

        auto& jh = junctions_[static_cast<std::size_t>(hi)];
        const std::int32_t de = jh.downstream_edge_id;
        if (de >= 0)
        {
            edges_[static_cast<std::size_t>(de)].is_loop_edge = true;
            jh.is_on_loop = true;
            const std::int32_t orec = jh.receiver_junction;
            if (orec >= 0)
                junctions_[static_cast<std::size_t>(orec)].is_on_loop = true;
        }
        jh.is_outlet          = true;
        jh.receiver_junction  = k_outlet_junction;
        jh.downstream_edge_id = -1;
        if (std::find(outlets_.begin(), outlets_.end(),
                      jh.junction_id) == outlets_.end())
        {
            outlets_.push_back(jh.junction_id);
        }
    }
}

// -----------------------------------------------------------------------------
// build_topological_order
// -----------------------------------------------------------------------------
void lsdtt_vector_to_junction_network::build_topological_order()
{
    topological_order_.clear();
    topological_order_.reserve(junctions_.size());

    std::vector<std::int32_t> remaining(
        junctions_.size(),
        0
    );

    for (std::size_t i = 0; i < junctions_.size(); ++i)
    {
        remaining[i] = static_cast<std::int32_t>(
            junctions_[i].donor_junctions.size()
        );
    }

    std::vector<std::int32_t> ready;
    ready.reserve(junctions_.size());

    for (std::int32_t j = 0; j < static_cast<std::int32_t>(junctions_.size()); ++j)
    {
        if (remaining[static_cast<std::size_t>(j)] == 0)
        {
            ready.push_back(j);
        }
    }

    std::size_t head = 0;
    while (head < ready.size())
    {
        const std::int32_t j = ready[head++];
        topological_order_.push_back(j);

        const std::int32_t rec = junctions_[static_cast<std::size_t>(j)].receiver_junction;
        if (rec >= 0)
        {
            --remaining[static_cast<std::size_t>(rec)];
            if (remaining[static_cast<std::size_t>(rec)] == 0)
            {
                ready.push_back(rec);
            }
        }
    }

    if (topological_order_.size() != junctions_.size())
    {
        throw std::runtime_error(
            "lsdtt_vector_to_junction_network: vector junction graph is cyclic or incomplete"
        );
    }

    has_topological_order_ = true;
}

// -----------------------------------------------------------------------------
// build_strahler_order
// -----------------------------------------------------------------------------
void lsdtt_vector_to_junction_network::build_strahler_order()
{
    if (!has_topological_order_)
    {
        throw std::runtime_error(
            "lsdtt_vector_to_junction_network: topological order must be built before Strahler order"
        );
    }

    // Reset all orders.
    for (auto& j : junctions_)
    {
        j.junction_strahler_order = 0;
        j.downstream_edge_strahler_order = 0;
    }

    for (auto& e : edges_)
    {
        e.strahler_order = 0;
    }

    // IMPORTANT:
    // Work from upstream to downstream using the topological order.
    for (const std::int32_t junction_id : topological_order_)
    {
        auto& j = junctions_[static_cast<std::size_t>(junction_id)];

        if (j.donor_junctions.empty())
        {
            j.junction_strahler_order = 1;
        }
        else
        {
            std::int32_t max_order = 0;
            int n_at_max = 0;

            for (const std::int32_t donor_id : j.donor_junctions)
            {
                const std::int32_t donor_order =
                    junctions_[static_cast<std::size_t>(donor_id)]
                        .junction_strahler_order;

                if (donor_order > max_order)
                {
                    max_order = donor_order;
                    n_at_max = 1;
                }
                else if (donor_order == max_order)
                {
                    ++n_at_max;
                }
            }

            j.junction_strahler_order =
                (n_at_max >= 2) ? (max_order + 1) : max_order;
        }
    }

    // The downstream edge leaving a junction gets that junction's Strahler order.
    for (auto& e : edges_)
    {
        if (e.is_loop_edge) { e.strahler_order = 0; continue; }
        e.strahler_order =
            junctions_[static_cast<std::size_t>(e.upstream_junction)]
                .junction_strahler_order;
    }

    for (auto& j : junctions_)
    {
        if (j.downstream_edge_id >= 0)
        {
            j.downstream_edge_strahler_order =
                edges_[static_cast<std::size_t>(j.downstream_edge_id)]
                    .strahler_order;
        }
        else
        {
            j.downstream_edge_strahler_order = 0;
        }
    }
}


// -----------------------------------------------------------------------------
// build_s_vector
// -----------------------------------------------------------------------------
void lsdtt_vector_to_junction_network::build_s_vector()
{
    s_vector_.clear();
    s_vector_.reserve(junctions_.size());

    struct dfs_frame
    {
        std::int32_t junction = k_nodata_junction;
        std::size_t next_donor = 0;
    };

    std::vector<dfs_frame> stack;
    stack.reserve(128);

    auto push_junction =
        [&](std::int32_t junction)
        {
            junctions_[static_cast<std::size_t>(junction)].s_vector_start =
                static_cast<std::int32_t>(s_vector_.size());

            junctions_[static_cast<std::size_t>(junction)].s_vector_size = 1;
            s_vector_.push_back(junction);
            stack.push_back({junction, 0});
        };

    for (const std::int32_t outlet : outlets_)
    {
        if (junctions_[static_cast<std::size_t>(outlet)].s_vector_start >= 0)
        {
            continue;
        }

        stack.clear();
        push_junction(outlet);

        while (!stack.empty())
        {
            dfs_frame& frame = stack.back();
            auto& donors =
                junctions_[static_cast<std::size_t>(frame.junction)].donor_junctions;

            if (frame.next_donor < donors.size())
            {
                const std::int32_t child = donors[frame.next_donor++];
                if (junctions_[static_cast<std::size_t>(child)].s_vector_start >= 0)
                {
                    throw std::runtime_error(
                        "lsdtt_vector_to_junction_network: duplicate visit while building s_vector"
                    );
                }

                push_junction(child);
                continue;
            }

            const std::int32_t done = frame.junction;
            const std::int32_t done_size =
                junctions_[static_cast<std::size_t>(done)].s_vector_size;

            stack.pop_back();

            if (!stack.empty())
            {
                junctions_[static_cast<std::size_t>(stack.back().junction)].s_vector_size +=
                    done_size;
            }
        }
    }

    if (s_vector_.size() != junctions_.size())
    {
        throw std::runtime_error(
            "lsdtt_vector_to_junction_network: s_vector does not cover all junctions"
        );
    }

    has_s_vector_ = true;
}

// -----------------------------------------------------------------------------
// write_junctions_fgb
// -----------------------------------------------------------------------------
void lsdtt_vector_to_junction_network::write_junctions_fgb(
    const std::string& out_path,
    bool spatial_index,
    const std::string& layer_name) const
{
    require_built();

    lsdtt_data_loader loader;

    lsdtt_vector_schema schema;
    schema.field_names = {
        "junction_id",
        "elevation",
        "is_dem_edge",
        "is_outlet",
        "receiver_junction",
        "outlet_junction",
        "downstream_edge_id",
        "downstream_edge_order",
        "junction_order",
        "n_donors",
        "donor_junctions",
        "n_incident_edges",
        "incident_edges",
        "s_vector_start",
        "s_vector_size",
        "is_on_loop"
    };
    schema.field_types = {
        lsdtt_vector_field_type::int64_type,
        lsdtt_vector_field_type::double_type,
        lsdtt_vector_field_type::int64_type,
        lsdtt_vector_field_type::int64_type,
        lsdtt_vector_field_type::int64_type,
        lsdtt_vector_field_type::int64_type,
        lsdtt_vector_field_type::int64_type,
        lsdtt_vector_field_type::int64_type,
        lsdtt_vector_field_type::int64_type,
        lsdtt_vector_field_type::int64_type,
        lsdtt_vector_field_type::string_type,
        lsdtt_vector_field_type::int64_type,
        lsdtt_vector_field_type::string_type,
        lsdtt_vector_field_type::int64_type,
        lsdtt_vector_field_type::int64_type,
        lsdtt_vector_field_type::int64_type
    };

    loader.begin_vector_write(
        out_path,
        layer_name,
        schema,
        meta_.crs_wkt,
        vector_output_format::flatgeobuf,
        spatial_index
    );

    for (const auto& j : junctions_)
    {
        loader.write_point(
            j.x,
            j.y,
            {
                static_cast<std::int64_t>(j.junction_id),
                j.elevation,
                static_cast<std::int64_t>(j.is_dem_edge_junction ? 1 : 0),
                static_cast<std::int64_t>(j.is_outlet ? 1 : 0),
                static_cast<std::int64_t>(j.receiver_junction),
                static_cast<std::int64_t>(j.outlet_junction),
                static_cast<std::int64_t>(j.downstream_edge_id),
                static_cast<std::int64_t>(j.downstream_edge_strahler_order),
                static_cast<std::int64_t>(j.junction_strahler_order),
                static_cast<std::int64_t>(j.donor_junctions.size()),
                join_ids_csv(j.donor_junctions),
                static_cast<std::int64_t>(j.incident_edges.size()),
                join_ids_csv(j.incident_edges),
                static_cast<std::int64_t>(j.s_vector_start),
                static_cast<std::int64_t>(j.s_vector_size),
                static_cast<std::int64_t>(j.is_on_loop ? 1 : 0)
            }
        );
    }

    loader.finish_vector_write();
}

// -----------------------------------------------------------------------------
// write_edges_fgb
// -----------------------------------------------------------------------------
void lsdtt_vector_to_junction_network::write_edges_fgb(
    const std::string& out_path,
    bool spatial_index,
    const std::string& layer_name) const
{
    require_built();

    lsdtt_data_loader loader;

    lsdtt_vector_schema schema;
    schema.field_names = {
        "edge_id",
        "strahler_order",
        "upstream_junction",
        "receiver_junction",
        "outlet_junction",
        "junction_a",
        "junction_b",
        "length_m",
        "is_loop_edge"
    };
    schema.field_types = {
        lsdtt_vector_field_type::int64_type,
        lsdtt_vector_field_type::int64_type,
        lsdtt_vector_field_type::int64_type,
        lsdtt_vector_field_type::int64_type,
        lsdtt_vector_field_type::int64_type,
        lsdtt_vector_field_type::int64_type,
        lsdtt_vector_field_type::int64_type,
        lsdtt_vector_field_type::double_type,
        lsdtt_vector_field_type::int64_type
    };


    loader.begin_vector_write(
        out_path,
        layer_name,
        schema,
        meta_.crs_wkt,
        vector_output_format::flatgeobuf,
        spatial_index
    );

    for (const auto& e : edges_)
    {
        loader.write_line(
            e.xy,
            {
                static_cast<std::int64_t>(e.edge_id),
                static_cast<std::int64_t>(e.strahler_order),
                static_cast<std::int64_t>(e.upstream_junction),
                static_cast<std::int64_t>(e.receiver_junction),
                static_cast<std::int64_t>(e.outlet_junction),
                static_cast<std::int64_t>(e.junction_a),
                static_cast<std::int64_t>(e.junction_b),
                e.length_metres,
                static_cast<std::int64_t>(e.is_loop_edge ? 1 : 0)
            }
        );
    }

    loader.finish_vector_write();
}

} // namespace lsdtt
