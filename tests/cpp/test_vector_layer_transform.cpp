// test_vector_layer_transform.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
//
// -----------------------------------------------------------------------------
// Characterisation test for make_layer_to_raster_transform and
// set_layer_spatial_filter_from_raster_bbox (src/core/lsdtt_vector_utils,
// consolidation C03).
//
// lsdtt_vector_rasterise_cpu and lsdtt_sti_cpu had identical copies of both
// (apart from the failure message, now a parameter, and a null-layer guard that
// only vector_rasterise had). Both removed copies are kept below verbatim.
//
//   transform: over pairs of layer SRS and raster WKT -- projected and
//   geographic, authority-compliant and traditional axis order on the layer,
//   identical CRSs, an empty raster WKT, an unparseable one, a null layer SRS --
//   the old and new functions must agree on whether a transform is made, on the
//   transformed coordinates of sample points (bit for bit), and on any exception
//   text.
//
//   spatial filter: on in-memory layers in each SRS (and none), the filter
//   rectangle set by the old and new functions must be identical (bit for bit).
//
// Usage: lsdtt-test-vector-layer-transform
// Prints one line per check and a final "RESULT: PASS|FAIL".
// -----------------------------------------------------------------------------

#include <algorithm>
#include <cstdint>
#include <cstring>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include <cpl_error.h>
#include <gdal.h>
#include <gdal_priv.h>
#include <ogr_spatialref.h>
#include <ogrsf_frmts.h>

#include "../../src/core/lsdtt_raster_meta.hpp"
#include "../../src/core/lsdtt_vector_utils.hpp"

using namespace lsdtt;

namespace
{

// ============================================================================
// Removed copies (references)
// ============================================================================

// src/cpu/lsdtt_vector_rasterise_cpu.cpp:34-42 (identical in lsdtt_sti_cpu.cpp:27-35)
struct oct_deleter
{
    void operator()(OGRCoordinateTransformation* tx) const
    {
        if (tx) OCTDestroyCoordinateTransformation(tx);
    }
};

using ct_ptr = std::unique_ptr<OGRCoordinateTransformation, oct_deleter>;

// src/cpu/lsdtt_vector_rasterise_cpu.cpp:314-343 (as of 2026-10-07)
ct_ptr make_ct_if_needed(
    const OGRSpatialReference* layer_srs,
    const std::string&   raster_wkt)
{
    if (!layer_srs || raster_wkt.empty()) return nullptr;

    OGRSpatialReference src_srs(*layer_srs);
    OGRSpatialReference raster_srs;

    src_srs.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);
    raster_srs.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);

    if (raster_srs.importFromWkt(raster_wkt.c_str()) != OGRERR_NONE)
        return nullptr;

    if (src_srs.IsSame(&raster_srs))
        return nullptr;

    OGRCoordinateTransformation* ct =
        OGRCreateCoordinateTransformation(&src_srs, &raster_srs);

    if (!ct)
    {
        throw std::runtime_error(
            "lsdtt_vector_rasterise: CRS transform failed - "
            "check that vector and raster have compatible projections");
    }

    return ct_ptr(ct);
}

// src/cpu/lsdtt_sti_cpu.cpp:218-254 (as of 2026-10-07); the vector_rasterise copy
// (:347-400) is the same body plus "if (!layer) return;".
void sti_set_spatial_filter_from_raster_bbox(OGRLayer* layer, const lsdtt_raster_meta& meta)
{
    double minx = meta.x_min;
    double miny = meta.y_min;
    double maxx = meta.x_max();
    double maxy = meta.y_max();

    const OGRSpatialReference* layer_srs = layer->GetSpatialRef();
    if (layer_srs && !meta.crs_wkt.empty())
    {
        OGRSpatialReference raster_srs;
        raster_srs.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);
        if (raster_srs.importFromWkt(meta.crs_wkt.c_str()) == OGRERR_NONE)
        {
            OGRSpatialReference dst_srs(*layer_srs);
            dst_srs.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);
            if (!raster_srs.IsSame(&dst_srs))
            {
                OGRCoordinateTransformation* raw_ct = OGRCreateCoordinateTransformation(&raster_srs, &dst_srs);
                if (raw_ct)
                {
                    ct_ptr ct(raw_ct);
                    double xs[4] = {meta.x_min, meta.x_max(), meta.x_min, meta.x_max()};
                    double ys[4] = {meta.y_min, meta.y_min, meta.y_max(), meta.y_max()};
                    if (ct->Transform(4, xs, ys))
                    {
                        minx = *std::min_element(xs, xs + 4);
                        maxx = *std::max_element(xs, xs + 4);
                        miny = *std::min_element(ys, ys + 4);
                        maxy = *std::max_element(ys, ys + 4);
                    }
                }
            }
        }
    }

    layer->SetSpatialFilterRect(minx, miny, maxx, maxy);
}

// ============================================================================
// Harness
// ============================================================================

int g_failures = 0;

void check(const std::string& what, bool ok, const std::string& detail = "")
{
    std::cout << (ok ? "PASS " : "FAIL ") << what;
    if (!ok && !detail.empty()) std::cout << " -- " << detail;
    std::cout << "\n";
    if (!ok) ++g_failures;
}

bool same_bits(double a, double b)
{
    std::uint64_t x, y;
    std::memcpy(&x, &a, sizeof x);
    std::memcpy(&y, &b, sizeof y);
    return x == y;
}

std::string wkt_of_epsg(int epsg)
{
    OGRSpatialReference s;
    s.importFromEPSG(epsg);
    char* w = nullptr;
    s.exportToWkt(&w);
    std::string out = w ? w : "";
    CPLFree(w);
    return out;
}

struct srs_case
{
    std::string name;
    int epsg;                  // 0 = no SRS
    bool traditional;          // axis mapping strategy on the layer SRS
};

// Sample points spread over each CRS's area of use (layer coordinates in
// traditional x/y order).
std::vector<std::pair<double, double>> samples_for(int epsg)
{
    switch (epsg)
    {
        case 4326:  return {{-3.2, 55.9}, {166.9, -11.6}, {0.0, 0.0}, {-120.5, 37.1}};
        case 32630: return {{400000, 6200000}, {500000, 5000000}, {700000, 100000}};
        case 32758: return {{714050, 8713685}, {500000, 9000000}};
        case 27700: return {{300000, 700000}, {530000, 180000}};
        default:    return {{1, 2}};
    }
}

} // namespace

int main()
{
    GDALAllRegister();
    // The unparseable-WKT cases make GDAL print "missing [" errors; they are
    // expected here, so keep them off the output.
    CPLPushErrorHandler(CPLQuietErrorHandler);
    std::cout << "== C03 layer-to-raster transform and spatial filter vs the removed copies ==\n";

    const std::vector<srs_case> layers = {
        {"none", 0, true},
        {"EPSG:4326 authority", 4326, false}, {"EPSG:4326 traditional", 4326, true},
        {"EPSG:32630 authority", 32630, false}, {"EPSG:32630 traditional", 32630, true},
        {"EPSG:27700 traditional", 27700, true}, {"EPSG:32758 traditional", 32758, true},
    };
    const std::vector<std::pair<std::string, std::string>> rasters = {
        {"empty", ""}, {"unparseable", "NOT A CRS"},
        {"EPSG:4326", wkt_of_epsg(4326)}, {"EPSG:32630", wkt_of_epsg(32630)},
        {"EPSG:27700", wkt_of_epsg(27700)}, {"EPSG:32758", wkt_of_epsg(32758)},
    };
    const std::string msg = "lsdtt_vector_rasterise: CRS transform failed - "
                            "check that vector and raster have compatible projections";

    // ---- transform ----
    {
        int n = 0, bad = 0, n_made = 0;
        std::string first;
        for (const srs_case& L : layers)
        {
            OGRSpatialReference srs;
            const OGRSpatialReference* srs_ptr = nullptr;
            if (L.epsg)
            {
                srs.importFromEPSG(L.epsg);
                srs.SetAxisMappingStrategy(L.traditional ? OAMS_TRADITIONAL_GIS_ORDER
                                                         : OAMS_AUTHORITY_COMPLIANT);
                srs_ptr = &srs;
            }
            for (const auto& R : rasters)
            {
                ++n;
                std::string e_old, e_new;
                ct_ptr a;
                lsdtt_ogr_transform_ptr b;
                try { a = make_ct_if_needed(srs_ptr, R.second); } catch (const std::exception& e) { e_old = e.what(); }
                try { b = make_layer_to_raster_transform(srs_ptr, R.second, msg); } catch (const std::exception& e) { e_new = e.what(); }
                bool ok = (e_old == e_new) && (static_cast<bool>(a) == static_cast<bool>(b));
                if (ok && a)
                {
                    ++n_made;
                    for (auto [x, y] : samples_for(L.epsg))
                    {
                        double xa = x, ya = y, xb = x, yb = y;
                        const int ra = a->Transform(1, &xa, &ya);
                        const int rb = b->Transform(1, &xb, &yb);
                        ok = ok && ra == rb && (!ra || (same_bits(xa, xb) && same_bits(ya, yb)));
                    }
                }
                if (!ok)
                {
                    ++bad;
                    if (first.empty()) first = L.name + " -> " + R.first;
                }
            }
        }
        check("transform made / not made, coordinates (bits) and exceptions identical over " +
                  std::to_string(n) + " SRS pairs (" + std::to_string(n_made) + " made a transform)",
              bad == 0, first);
    }

    // ---- spatial filter ----
    {
        GDALDriver* mem = GetGDALDriverManager()->GetDriverByName("Memory");
        int n = 0, bad = 0;
        std::string first;
        for (const srs_case& L : layers)
        {
            for (const auto& R : rasters)
            {
                lsdtt_raster_meta meta;
                meta.n_rows = 37;
                meta.n_cols = 53;
                meta.cell_size = (R.first == "EPSG:4326") ? 0.01 : 30.0;
                meta.x_min = (R.first == "EPSG:4326") ? -3.5 : 400000.0;
                meta.y_min = (R.first == "EPSG:4326") ? 55.5 : 6200000.0;
                meta.crs_wkt = R.second;

                double env_old[4], env_new[4];
                for (int which = 0; which < 2; ++which)
                {
                    GDALDataset* ds = mem->Create("", 0, 0, 0, GDT_Unknown, nullptr);
                    OGRSpatialReference srs;
                    OGRSpatialReference* srs_ptr = nullptr;
                    if (L.epsg)
                    {
                        srs.importFromEPSG(L.epsg);
                        srs.SetAxisMappingStrategy(L.traditional ? OAMS_TRADITIONAL_GIS_ORDER
                                                                 : OAMS_AUTHORITY_COMPLIANT);
                        srs_ptr = &srs;
                    }
                    OGRLayer* layer = ds->CreateLayer("l", srs_ptr, wkbPoint, nullptr);
                    if (which == 0) sti_set_spatial_filter_from_raster_bbox(layer, meta);
                    else            set_layer_spatial_filter_from_raster_bbox(layer, meta);
                    OGREnvelope env;
                    layer->GetSpatialFilter()->getEnvelope(&env);
                    double* out = (which == 0) ? env_old : env_new;
                    out[0] = env.MinX; out[1] = env.MinY; out[2] = env.MaxX; out[3] = env.MaxY;
                    GDALClose(ds);
                }
                ++n;
                bool ok = true;
                for (int k = 0; k < 4; ++k) ok = ok && same_bits(env_old[k], env_new[k]);
                if (!ok)
                {
                    ++bad;
                    if (first.empty()) first = L.name + " / raster " + R.first;
                }
            }
        }
        check("spatial filter rectangle identical (bits) over " + std::to_string(n) + " layer/raster pairs",
              bad == 0, first);
        set_layer_spatial_filter_from_raster_bbox(nullptr, lsdtt_raster_meta{});
        check("null layer is a no-op (the vector_rasterise guard)", true);
    }

    std::cout << "RESULT: " << (g_failures == 0 ? "PASS" : "FAIL") << "\n";
    return g_failures == 0 ? 0 : 1;
}
