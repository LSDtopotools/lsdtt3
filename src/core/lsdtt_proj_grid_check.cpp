// lsdtt_proj_grid_check.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html

#include "lsdtt_proj_grid_check.hpp"
#include "lsdtt_gdal_init.hpp"

#include <ogr_spatialref.h>

#include <cmath>
#include <filesystem>
#include <iostream>
#include <mutex>
#include <set>
#include <sstream>

#if defined(LSDTT_HAS_PROJ_API)
#include <proj.h>
#endif

namespace lsdtt {

const std::vector<std::string>& default_proj_grids()
{
    static const std::vector<std::string> grids = {
        "uk_os_OSTN15_NTv2_OSGBtoETRS.tif",           // OSGB36 (British National Grid) <-> ETRS89
        "us_nga_egm96_15.tif",                        // EGM96 geoid (vertical_datum: egm96)
        "us_noaa_conus.tif",                          // NADCON: NAD27 <-> NAD83, contiguous US
        "us_noaa_nadcon5_nad27_nad83_1986_conus.tif", // NADCON5: NAD27 <-> NAD83(1986), contiguous US
    };
    return grids;
}

std::string proj_grids_description()
{
    const std::string dir = pinned_proj_data_dir();
    if (dir.empty()) return "no bundled data directory (PROJ's own search)";
    namespace fs = std::filesystem;
    std::size_t n_default = 0;
    std::vector<std::string> missing, extra;
    for (const auto& g : default_proj_grids())
    {
        std::error_code ec;
        if (fs::is_regular_file(fs::path(dir) / g, ec)) ++n_default;
        else missing.push_back(g);
    }
    std::error_code ec;
    for (const auto& e : fs::directory_iterator(dir, ec))
    {
        const std::string name = e.path().filename().string();
        const std::string ext = e.path().extension().string();
        if (ext != ".tif" && ext != ".tiff" && ext != ".gsb" && ext != ".gtx") continue;
        bool is_default = false;
        for (const auto& g : default_proj_grids()) is_default = is_default || g == name;
        if (!is_default) extra.push_back(name);
    }
    std::ostringstream s;
    s << n_default << " of " << default_proj_grids().size() << " default";
    if (!missing.empty())
    {
        s << " (missing:";
        for (const auto& m : missing) s << ' ' << m;
        s << ')';
    }
    if (!extra.empty()) s << ", " << extra.size() << " added";
    return s.str();
}

#if defined(LSDTT_HAS_PROJ_API)

namespace {

// The check is meaningful only if the PROJ we call is the one GDAL runs with.
std::string disabled_reason()
{
    gdal_init();
    int major = 0, minor = 0, patch = 0;
    OSRGetPROJVersion(&major, &minor, &patch);
    const PJ_INFO info = proj_info();
    if (major != info.major)
        return "PROJ " + std::to_string(info.major) + " (compiled) differs from GDAL's PROJ " +
               std::to_string(major);
    return "";
}

std::string accuracy_text(double a)
{
    if (!(a >= 0.0)) return "accuracy unknown";
    std::ostringstream s;
    s << "accuracy " << a << " m";
    return s.str();
}

} // namespace

std::string proj_grid_check_status()
{
    const std::string why = disabled_reason();
    return why.empty() ? "on" : "off: " + why;
}

proj_grid_check_result check_proj_transformation_grids(
    const std::string& src, const std::string& dst, const double* lonlat_box)
{
    proj_grid_check_result out;
    if (!disabled_reason().empty()) return out;

    PJ_CONTEXT* ctx = proj_context_create();
    if (!ctx) return out;
    proj_log_level(ctx, PJ_LOG_NONE);
    const std::string pinned = pinned_proj_data_dir();
    if (!pinned.empty())
    {
        const char* paths[] = {pinned.c_str()};
        proj_context_set_search_paths(ctx, 1, paths);
    }
    PJ* s = proj_create(ctx, src.c_str());
    PJ* d = proj_create(ctx, dst.c_str());
    PJ_OPERATION_FACTORY_CONTEXT* f = proj_create_operation_factory_context(ctx, nullptr);
    PJ_OBJ_LIST* ops = nullptr;
    if (s && d && f)
    {
        proj_operation_factory_context_set_grid_availability_use(ctx, f, PROJ_GRID_AVAILABILITY_IGNORED);
        proj_operation_factory_context_set_spatial_criterion(ctx, f, PROJ_SPATIAL_CRITERION_PARTIAL_INTERSECTION);
        if (lonlat_box)
            proj_operation_factory_context_set_area_of_interest(ctx, f, lonlat_box[0], lonlat_box[1],
                                                                lonlat_box[2], lonlat_box[3]);
        ops = proj_create_operations(ctx, s, d, f);
    }
    if (ops)
    {
        const int n = proj_list_get_count(ops);
        bool have_used = false;
        for (int i = 0; i < n && !have_used; ++i)
        {
            PJ* op = proj_list_get(ctx, ops, i);
            if (!op) continue;
            std::vector<std::string> missing;
            const int ng = proj_coordoperation_get_grid_used_count(ctx, op);
            for (int j = 0; j < ng; ++j)
            {
                const char *short_name = nullptr, *full_name = nullptr, *package = nullptr, *url = nullptr;
                int direct = 0, open_licence = 0, available = 0;
                if (proj_coordoperation_get_grid_used(ctx, op, j, &short_name, &full_name, &package, &url,
                                                      &direct, &open_licence, &available) &&
                    !available && short_name)
                    missing.emplace_back(short_name);
            }
            const char* name = proj_get_name(op);
            const double acc = proj_coordoperation_get_accuracy(ctx, op);
            if (i == 0)
            {
                out.evaluated = true;
                out.best_operation = name ? name : "";
                out.best_accuracy_m = acc;
                out.missing_grids = missing;
            }
            if (missing.empty())
            {
                out.used_operation = name ? name : "";
                out.used_accuracy_m = acc;
                have_used = true;
            }
            proj_destroy(op);
        }
        if (n == 0) out.evaluated = true;   // nothing to choose between
        proj_list_destroy(ops);
    }
    if (f) proj_operation_factory_context_destroy(f);
    if (s) proj_destroy(s);
    if (d) proj_destroy(d);
    proj_context_destroy(ctx);
    return out;
}

void warn_if_proj_grids_missing(const OGRSpatialReference& src, const OGRSpatialReference& dst,
                                const std::string& context, const double* lonlat_box)
{
    try
    {
        char* a = nullptr;
        char* b = nullptr;
        const char* opts[] = {"FORMAT=WKT2_2019", nullptr};
        if (src.exportToWkt(&a, opts) != OGRERR_NONE || dst.exportToWkt(&b, opts) != OGRERR_NONE)
        {
            CPLFree(a);
            CPLFree(b);
            return;
        }
        const std::string sa = a, sb = b;
        CPLFree(a);
        CPLFree(b);
        static std::mutex m;
        static std::set<std::string> seen;
        {
            std::lock_guard<std::mutex> lock(m);
            if (!seen.insert(sa + "\n->\n" + sb).second) return;
        }
        const proj_grid_check_result r = check_proj_transformation_grids(sa, sb, lonlat_box);
        if (!r.evaluated || r.missing_grids.empty()) return;
        std::ostringstream w;
        w << "WARNING: " << context << ": the most accurate transformation (" << r.best_operation << ", "
          << accuracy_text(r.best_accuracy_m) << ") needs the PROJ grid";
        if (r.missing_grids.size() > 1) w << "s";
        for (std::size_t k = 0; k < r.missing_grids.size(); ++k) w << (k ? ", " : " ") << r.missing_grids[k];
        w << ", which is not installed. Using "
          << (r.used_operation.empty() ? std::string("a ballpark transformation") : r.used_operation) << " ("
          << accuracy_text(r.used_accuracy_m) << ") instead. To install it: python3 fetch_proj_grids.py fetch";
        for (const auto& g : r.missing_grids) w << ' ' << g;
        w << " (the script is tools/fetch_proj_grids.py in the lsdtt3 repository; see "
             "https://lsdtopotools.github.io/lsdtt3/what_is_lsdtt3/note_on_grids.html)\n";
        std::cerr << w.str();
    }
    catch (...)
    {
        // A diagnostic must never stop a run.
    }
}

void warn_if_proj_grids_missing(const OGRSpatialReference& src, const OGRSpatialReference& dst,
                                const std::string& context, double x0, double y0, double x1, double y1)
{
    try
    {
        OGRSpatialReference s(src), ll;
        ll.importFromEPSG(4326);
        s.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);
        ll.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);
        double box[4] = {0, 0, 0, 0};
        OGRCoordinateTransformation* ct = OGRCreateCoordinateTransformation(&s, &ll);
        bool ok = ct && ct->TransformBounds(x0, y0, x1, y1, &box[0], &box[1], &box[2], &box[3], 21);
        if (ct) OGRCoordinateTransformation::DestroyCT(ct);
        warn_if_proj_grids_missing(src, dst, context, ok ? box : nullptr);
    }
    catch (...)
    {
    }
}

#else  // no PROJ API

std::string proj_grid_check_status() { return "off: built without the PROJ headers"; }

proj_grid_check_result check_proj_transformation_grids(const std::string&, const std::string&, const double*)
{
    return {};
}

void warn_if_proj_grids_missing(const OGRSpatialReference&, const OGRSpatialReference&, const std::string&,
                                const double*)
{
}

void warn_if_proj_grids_missing(const OGRSpatialReference&, const OGRSpatialReference&, const std::string&,
                                double, double, double, double)
{
}

#endif

} // namespace lsdtt
