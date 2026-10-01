// lsdtt_hype_raster_setup_cpu.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// HYPE set-up from rasters. See the header.

#include "lsdtt_hype_raster_setup_cpu.hpp"

#include <cmath>
#include <fstream>
#include <iomanip>
#include <map>
#include <sstream>
#include <stdexcept>
#include <utility>

namespace lsdtt
{

namespace
{
int integer_code(float v, const char* what)
{
    const double d = static_cast<double>(v);
    if (!std::isfinite(d) || d != std::floor(d))
    {
        std::ostringstream os;
        os << "hype raster set-up: " << what << " raster holds a non-integer code (" << d << ")";
        throw std::invalid_argument(os.str());
    }
    return static_cast<int>(d);
}
} // namespace

lsdtt_raster_cpu make_hype_class_raster(const lsdtt_raster_cpu& landuse_raster,
                                        const lsdtt_raster_cpu& soil_raster,
                                        const std::vector<hype_class_definition>& classes)
{
    const auto& lm = landuse_raster.get_meta();
    const auto& sm = soil_raster.get_meta();
    if (lm.n_rows != sm.n_rows || lm.n_cols != sm.n_cols)
    {
        throw std::invalid_argument("hype raster set-up: land-use and soil rasters differ in shape");
    }
    std::map<std::pair<int, int>, int> class_of;
    for (const auto& c : classes)
    {
        if (!class_of.emplace(std::make_pair(c.landuse, c.soil), c.class_id).second)
        {
            throw std::invalid_argument("hype raster set-up: classes " +
                std::to_string(class_of.at({c.landuse, c.soil})) + " and " + std::to_string(c.class_id) +
                " both use land use " + std::to_string(c.landuse) + " with soil " + std::to_string(c.soil));
        }
    }

    const auto& ld = landuse_raster.get_data();
    const auto& lk = landuse_raster.get_mask();
    const auto& sd = soil_raster.get_data();
    const auto& sk = soil_raster.get_mask();
    std::vector<float>        data(ld.size(), static_cast<float>(lm.no_data_value));
    std::vector<std::uint8_t> mask(ld.size(), 0u);
    std::map<std::pair<int, int>, std::int64_t> missing;
    for (std::size_t i = 0; i < ld.size(); ++i)
    {
        if (!lk[i] || !sk[i]) continue;
        const std::pair<int, int> key{integer_code(ld[i], "land-use"), integer_code(sd[i], "soil")};
        const auto it = class_of.find(key);
        if (it == class_of.end())
        {
            ++missing[key];
            continue;
        }
        data[i] = static_cast<float>(it->second);
        mask[i] = 1u;
    }
    if (!missing.empty())
    {
        std::ostringstream os;
        os << "hype raster set-up: the class table has no class for these (land use, soil) pairs found in the rasters:";
        for (const auto& [key, n] : missing) os << " (" << key.first << ", " << key.second << ") on " << n << " pixels;";
        os << " add them to the class table";
        throw std::invalid_argument(os.str());
    }
    return lsdtt_raster_cpu(std::move(data), std::move(mask), lm);
}

hype_raster_setup_result build_hype_setup_from_rasters(
    const lsdtt_flow_base& flow, const lsdtt_raster_cpu& dem, const lsdtt_raster_cpu& gradient,
    const lsdtt_raster_cpu& class_raster, const std::vector<std::int64_t>& outlet_flats,
    std::int64_t channel_threshold_contributing_pixels, const hype_setup& parameters)
{
    hype_raster_setup_result r;
    r.partition     = partition_incremental_subbasins(flow, outlet_flats);
    const auto geo  = compute_subbasin_geodata(r.partition, dem, gradient);
    const auto frac = compute_subbasin_class_fractions(r.partition, class_raster, dem);
    r.river_lengths = compute_subbasin_river_lengths(r.partition, flow, channel_threshold_contributing_pixels);

    r.setup = parameters;
    r.setup.subbasins.clear();
    r.setup.subbasin_classes.clear();
    r.setup.source_dir = "(built from rasters)";

    const std::int32_t n = r.partition.n_subbasins();
    for (std::int32_t s = 0; s < n; ++s)
    {
        const std::size_t si = static_cast<std::size_t>(s);
        if (frac.n_valid_class[si] == 0)
        {
            throw std::invalid_argument("hype raster set-up: sub-basin " + std::to_string(s + 1) +
                " (outlet flat " + std::to_string(r.partition.outlet_flat[si]) +
                ") has no pixel with a valid class; its class fractions cannot sum to 1");
        }
        hype_subbasin sb;
        sb.subid              = s + 1;
        sb.maindown           = (r.partition.downstream[si] < 0) ? 0 : r.partition.downstream[si] + 1;
        sb.area_m2            = geo.area_m2[si];
        sb.elev_mean_m        = geo.elev_mean_m[si];
        sb.elev_std_m         = geo.elev_std_m[si];
        sb.slope_mean_percent = geo.slope_mean_percent[si];
        sb.loc_rivlen_m       = r.river_lengths.loc_rivlen_m[si];
        sb.rivlen_m           = r.river_lengths.rivlen_m[si];
        r.setup.subbasins.push_back(sb);
        for (std::int64_t e = frac.offset[si]; e < frac.offset[si + 1]; ++e)
        {
            const auto& ent = frac.entries[static_cast<std::size_t>(e)];
            hype_subbasin_class sc;
            sc.subid    = s + 1;
            sc.class_id = ent.class_code;
            sc.fraction = ent.fraction;
            sc.dh_m     = ent.dh_m;
            r.setup.subbasin_classes.push_back(sc);
        }
    }
    validate_hype_setup(r.setup);
    return r;
}

void write_hype_geodata_files(const hype_setup& setup, const std::string& dir)
{
    std::vector<int> class_ids;
    for (const auto& c : setup.classes) class_ids.push_back(c.class_id);
    std::map<std::pair<int, int>, const hype_subbasin_class*> sc_of;
    for (const auto& sc : setup.subbasin_classes) sc_of[{sc.subid, sc.class_id}] = &sc;

    {
        std::ofstream f(dir + "/GeoData.txt");
        if (!f) throw std::runtime_error("cannot write " + dir + "/GeoData.txt");
        f << std::setprecision(10);
        f << "subid\tmaindown\tarea\telev_mean\telev_std\tslope_mean\tloc_rivlen\trivlen";
        for (int c : class_ids) f << "\tslc_" << c;
        for (int c : class_ids) f << "\tdhslc_" << c;
        f << "\n";
        for (const auto& sb : setup.subbasins)
        {
            f << sb.subid << '\t' << sb.maindown << '\t' << sb.area_m2 << '\t' << sb.elev_mean_m << '\t'
              << sb.elev_std_m << '\t' << sb.slope_mean_percent << '\t' << sb.loc_rivlen_m << '\t' << sb.rivlen_m;
            for (int c : class_ids)
            {
                const auto it = sc_of.find({sb.subid, c});
                f << '\t' << (it == sc_of.end() ? 0.0 : it->second->fraction);
            }
            for (int c : class_ids)
            {
                const auto it = sc_of.find({sb.subid, c});
                f << '\t' << (it == sc_of.end() ? 0.0 : it->second->dh_m);
            }
            f << "\n";
        }
    }
    {
        std::ofstream f(dir + "/GeoClass.txt");
        if (!f) throw std::runtime_error("cannot write " + dir + "/GeoClass.txt");
        f << std::setprecision(10);
        f << "!class\tlanduse\tsoil\tcrop1\tcrop2\trotation\tvegtype\tspecial\ttiledepth\tstreamdepth\tnlayers\tdepth1\tdepth2\tdepth3\n";
        for (const auto& c : setup.classes)
        {
            f << c.class_id << '\t' << c.landuse << '\t' << c.soil << "\t0\t0\t0\t1\t0\t"
              << c.tile_depth_m << '\t' << c.stream_depth_m << '\t' << c.n_layers;
            for (int k = 0; k < c.n_layers; ++k) f << '\t' << c.layer_bottom_m[static_cast<std::size_t>(k)];
            f << "\n";
        }
    }
}

} // namespace lsdtt
