#include "lsdtt_sti_cpu.hpp"
#include "../core/lsdtt_gdal_init.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <limits>
#include <memory>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <vector>

#include <gdal.h>
#include <gdal_alg.h>
#include <gdal_priv.h>
#include <ogrsf_frmts.h>

namespace lsdtt {
namespace {

struct oct_deleter
{
    void operator()(OGRCoordinateTransformation* tx) const
    {
        if (tx) OCTDestroyCoordinateTransformation(tx);
    }
};

using ct_ptr = std::unique_ptr<OGRCoordinateTransformation, oct_deleter>;

std::string trim_copy(const std::string& s)
{
    const auto is_not_space = [](unsigned char c) { return !std::isspace(c); };
    auto b = std::find_if(s.begin(), s.end(), is_not_space);
    if (b == s.end()) return "";
    auto e = std::find_if(s.rbegin(), s.rend(), is_not_space).base();
    return std::string(b, e);
}

std::string to_upper_copy(std::string s)
{
    std::transform(s.begin(), s.end(), s.begin(),
                   [](unsigned char c) { return static_cast<char>(std::toupper(c)); });
    return s;
}

std::vector<std::string> parse_csv_line(const std::string& line)
{
    std::vector<std::string> out;
    std::string token;
    bool in_quotes = false;

    for (std::size_t i = 0; i < line.size(); ++i)
    {
        const char c = line[i];
        if (in_quotes)
        {
            if (c == '"')
            {
                if (i + 1 < line.size() && line[i + 1] == '"')
                {
                    token.push_back('"');
                    ++i;
                }
                else
                {
                    in_quotes = false;
                }
            }
            else
            {
                token.push_back(c);
            }
        }
        else
        {
            if (c == ',')
            {
                out.push_back(token);
                token.clear();
            }
            else if (c == '"')
            {
                in_quotes = true;
            }
            else
            {
                token.push_back(c);
            }
        }
    }

    out.push_back(token);
    return out;
}

bool parse_double_safe(const std::string& s, double& out)
{
    const std::string t = trim_copy(s);
    if (t.empty()) return false;
    try
    {
        std::size_t idx = 0;
        const double v = std::stod(t, &idx);
        if (idx != t.size()) return false;
        if (!std::isfinite(v)) return false;
        out = v;
        return true;
    }
    catch (...)
    {
        return false;
    }
}

std::string make_pair_key(const std::string& series_code, const std::string& land_use)
{
    return series_code + "|" + land_use;
}

std::pair<std::string, std::string> split_pair_key(const std::string& key)
{
    const std::size_t p = key.find('|');
    if (p == std::string::npos) return std::make_pair(key, std::string());
    return std::make_pair(key.substr(0, p), key.substr(p + 1));
}

std::vector<std::string> fallback_land_use_order(const std::string& requested)
{
    if (requested == "CULT") return {"CULT", "SEMI", "OTHER"};
    if (requested == "SEMI") return {"SEMI", "CULT", "OTHER"};
    if (requested == "OTHER") return {"OTHER", "CULT", "SEMI"};
    return {requested, "CULT", "SEMI", "OTHER"};
}

double compute_spr_term(double spr_fraction, sti_spr_transform transform, double eps, bool& clamped)
{
    clamped = false;
    double spr = spr_fraction;
    if (spr <= eps)
    {
        spr = eps;
        clamped = true;
    }
    if (spr >= 1.0 - eps)
    {
        spr = 1.0 - eps;
        clamped = true;
    }

    if (transform == sti_spr_transform::ln_spr)
    {
        return std::log(spr);
    }

    const double om = std::max(eps, 1.0 - spr);
    if (om == eps) clamped = true;
    return std::log(om);
}

void fill_window_geotransform(const lsdtt_raster_meta& meta, int row0, int col0, int oversample, double gt[6])
{
    const double raster_top = meta.y_min + static_cast<double>(meta.n_rows) * meta.cell_size;
    gt[0] = meta.x_min + static_cast<double>(col0) * meta.cell_size;
    gt[1] = meta.cell_size / static_cast<double>(oversample);
    gt[2] = 0.0;
    gt[3] = raster_top - static_cast<double>(row0) * meta.cell_size;
    gt[4] = 0.0;
    gt[5] = -meta.cell_size / static_cast<double>(oversample);
}

GDALDataset* make_mem_byte_raster(int nrows, int ncols, const double gt[6])
{
    GDALDriver* drv = GetGDALDriverManager()->GetDriverByName("MEM");
    if (!drv) throw std::runtime_error("rasterise_polygon_string_field: MEM driver unavailable");
    GDALDataset* ds = drv->Create("", ncols, nrows, 1, GDT_Byte, nullptr);
    if (!ds) throw std::runtime_error("rasterise_polygon_string_field: failed to create MEM raster");
    ds->SetGeoTransform(const_cast<double*>(gt));
    if (ds->GetRasterBand(1)->Fill(0.0) != CE_None)
    {
        GDALClose(ds);
        throw std::runtime_error("rasterise_polygon_string_field: failed to initialise MEM raster");
    }
    return ds;
}

ct_ptr make_ct_if_needed(const OGRSpatialReference* layer_srs, const std::string& raster_wkt)
{
    if (!layer_srs || raster_wkt.empty()) return nullptr;

    OGRSpatialReference src_srs(*layer_srs);
    OGRSpatialReference raster_srs;
    src_srs.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);
    raster_srs.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);
    if (raster_srs.importFromWkt(raster_wkt.c_str()) != OGRERR_NONE) return nullptr;
    if (src_srs.IsSame(&raster_srs)) return nullptr;

    OGRCoordinateTransformation* ct = OGRCreateCoordinateTransformation(&src_srs, &raster_srs);
    if (!ct)
    {
        throw std::runtime_error("rasterise_polygon_string_field: failed to create CRS transform");
    }
    return ct_ptr(ct);
}

bool is_polygon_type(OGRwkbGeometryType t)
{
    t = wkbFlatten(t);
    return t == wkbPolygon || t == wkbMultiPolygon || t == wkbGeometryCollection || t == wkbUnknown;
}

void set_spatial_filter_from_raster_bbox(OGRLayer* layer, const lsdtt_raster_meta& meta)
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

bool raster_window_from_envelope(const lsdtt_raster_meta& meta, const OGREnvelope& env,
                                 int& row0, int& row1, int& col0, int& col1)
{
    const double raster_minx = meta.x_min;
    const double raster_maxx = meta.x_max();
    const double raster_miny = meta.y_min;
    const double raster_maxy = meta.y_max();
    const double cs = meta.cell_size;

    const double minx = std::max(env.MinX, raster_minx);
    const double maxx = std::min(env.MaxX, raster_maxx);
    const double miny = std::max(env.MinY, raster_miny);
    const double maxy = std::min(env.MaxY, raster_maxy);
    if (!(minx < maxx && miny < maxy)) return false;

    col0 = std::max(0, std::min(meta.n_cols - 1, static_cast<int>(std::floor((minx - raster_minx) / cs))));
    col1 = std::max(0, std::min(meta.n_cols - 1, static_cast<int>(std::ceil((maxx - raster_minx) / cs)) - 1));
    row0 = std::max(0, std::min(meta.n_rows - 1, static_cast<int>(std::floor((raster_maxy - maxy) / cs))));
    row1 = std::max(0, std::min(meta.n_rows - 1, static_cast<int>(std::ceil((raster_maxy - miny) / cs)) - 1));

    return row0 <= row1 && col0 <= col1;
}

void csv_write_escaped(std::ofstream& ofs, const std::string& s)
{
    bool needs_quotes = false;
    for (char c : s)
    {
        if (c == ',' || c == '"' || c == '\n' || c == '\r')
        {
            needs_quotes = true;
            break;
        }
    }

    if (!needs_quotes)
    {
        ofs << s;
        return;
    }

    ofs << '"';
    for (char c : s)
    {
        if (c == '"') ofs << '"';
        ofs << c;
    }
    ofs << '"';
}

struct ksat_profile_accum
{
    std::string series_code;
    std::string land_use;
    double transmissivity_proxy = 0.0;
    double integrated_depth_cm = 0.0;
};

double median_in_place(std::vector<double>& values)
{
    if (values.empty()) return std::numeric_limits<double>::quiet_NaN();
    std::sort(values.begin(), values.end());
    const std::size_t n = values.size();
    if (n % 2 == 1) return values[n / 2];
    return 0.5 * (values[(n / 2) - 1] + values[n / 2]);
}

std::unordered_map<std::string, std::size_t> header_index_map(const std::vector<std::string>& headers)
{
    std::unordered_map<std::string, std::size_t> out;
    out.reserve(headers.size());
    for (std::size_t i = 0; i < headers.size(); ++i)
    {
        out.emplace(to_upper_copy(trim_copy(headers[i])), i);
    }
    return out;
}

std::string fetch_cell(const std::vector<std::string>& row,
                       const std::unordered_map<std::string, std::size_t>& header_map,
                       const std::string& key)
{
    auto it = header_map.find(key);
    if (it == header_map.end()) return "";
    const std::size_t idx = it->second;
    if (idx >= row.size()) return "";
    return row[idx];
}

} // namespace

// Public: see the header for why this is exported.
bool lookup_ln_kd_with_fallback(
    const sti_lnkd_table& table,
    const std::string& series_code,
    const std::string& land_use,
    double& out_ln_kd,
    std::string& out_source,
    int& out_profile_count)
{
    out_ln_kd = 0.0;
    out_profile_count = 0;
    out_source = "missing";

    if (series_code.empty() || land_use.empty()) return false;

    const auto try_land_uses = fallback_land_use_order(land_use);
    for (std::size_t i = 0; i < try_land_uses.size(); ++i)
    {
        const std::string key = make_pair_key(series_code, try_land_uses[i]);
        auto it = table.ln_kd_by_key.find(key);
        if (it == table.ln_kd_by_key.end()) continue;

        out_ln_kd = it->second;
        auto nit = table.n_profiles_by_key.find(key);
        out_profile_count = (nit == table.n_profiles_by_key.end()) ? 0 : nit->second;
        out_source = (i == 0) ? "direct" : "fallback_other_landuse";
        return true;
    }

    return false;
}


std::string normalize_sti_series_code(const std::string& raw_code)
{
    std::string code = trim_copy(raw_code);
    if (code.empty()) return "";

    bool all_digits = true;
    for (char c : code)
    {
        if (c < '0' || c > '9')
        {
            all_digits = false;
            break;
        }
    }

    if (all_digits)
    {
        std::size_t first_non_zero = code.find_first_not_of('0');
        if (first_non_zero == std::string::npos) return "0";
        return code.substr(first_non_zero);
    }

    double numeric = 0.0;
    if (parse_double_safe(code, numeric))
    {
        const double rounded = std::round(numeric);
        if (std::fabs(numeric - rounded) < 1.0e-9)
        {
            const long long iv = static_cast<long long>(rounded);
            return std::to_string(iv);
        }
    }

    return to_upper_copy(code);
}

std::string normalize_sti_land_use(const std::string& raw_land_use)
{
    const std::string u = to_upper_copy(trim_copy(raw_land_use));
    if (u.empty()) return "";
    if (u == "CULT" || u.find("CULT") != std::string::npos) return "CULT";
    if (u == "SEMI" || u.find("SEMI") != std::string::npos) return "SEMI";
    if (u == "OTHER") return "OTHER";
    return "OTHER";
}

sti_lnkd_build_result build_sti_lnkd_table_from_csv(
    const std::string& csv_path,
    const sti_lnkd_options& options)
{
    std::ifstream ifs(csv_path);
    if (!ifs.is_open())
    {
        throw std::runtime_error("build_sti_lnkd_table_from_csv: could not open: " + csv_path);
    }

    std::string header_line;
    if (!std::getline(ifs, header_line))
    {
        throw std::runtime_error("build_sti_lnkd_table_from_csv: empty CSV: " + csv_path);
    }

    const auto headers = parse_csv_line(header_line);
    const auto hmap = header_index_map(headers);
    const std::vector<std::string> required = {
        "SERIES_CODE", "SSKIB_PROFILE_ID", "LUSE", "HORZ_TOP", "HORZ_BOTTOM", "COND_0_HYPRES"
    };
    for (const std::string& key : required)
    {
        if (hmap.find(key) == hmap.end())
        {
            throw std::runtime_error("build_sti_lnkd_table_from_csv: missing required column '" + key + "'");
        }
    }

    std::unordered_map<std::string, ksat_profile_accum> profile_accum;

    sti_lnkd_build_result result;
    std::string line;
    while (std::getline(ifs, line))
    {
        if (line.empty()) continue;
        ++result.rows_read;
        const auto row = parse_csv_line(line);

        const std::string series_code = normalize_sti_series_code(fetch_cell(row, hmap, "SERIES_CODE"));
        const std::string profile_id = trim_copy(fetch_cell(row, hmap, "SSKIB_PROFILE_ID"));
        const std::string land_use = normalize_sti_land_use(fetch_cell(row, hmap, "LUSE"));
        if (series_code.empty() || profile_id.empty() || land_use.empty()) continue;

        double top_cm = 0.0;
        double bot_cm = 0.0;
        double ks_proxy = 0.0;
        if (!parse_double_safe(fetch_cell(row, hmap, "HORZ_TOP"), top_cm)) continue;
        if (!parse_double_safe(fetch_cell(row, hmap, "HORZ_BOTTOM"), bot_cm)) continue;
        if (!parse_double_safe(fetch_cell(row, hmap, "COND_0_HYPRES"), ks_proxy)) continue;

        if (bot_cm < top_cm) std::swap(top_cm, bot_cm);

        const double clipped_top = std::max(0.0, std::min(options.depth_cap_cm, top_cm));
        const double clipped_bot = std::max(0.0, std::min(options.depth_cap_cm, bot_cm));
        const double dz_cm = std::max(0.0, clipped_bot - clipped_top);
        if (dz_cm <= 0.0) continue;

        const double ks_scaled = ks_proxy * options.ks_unit_scale_to_m_per_day;
        const std::string pkey = series_code + "|" + land_use + "|" + profile_id;

        auto& acc = profile_accum[pkey];
        acc.series_code = series_code;
        acc.land_use = land_use;
        acc.integrated_depth_cm += dz_cm;
        acc.transmissivity_proxy += ks_scaled * dz_cm;
        ++result.rows_used;
    }

    result.profiles_total = profile_accum.size();

    std::unordered_map<std::string, std::vector<double>> lnkd_values_by_pair;
    for (const auto& kv : profile_accum)
    {
        const auto& acc = kv.second;
        if (acc.integrated_depth_cm < options.min_valid_depth_cm) continue;
        const double t_proxy = std::max(options.eps, acc.transmissivity_proxy);
        const double ln_kd = std::log(t_proxy);
        const std::string key = make_pair_key(acc.series_code, acc.land_use);
        lnkd_values_by_pair[key].push_back(ln_kd);
        ++result.profiles_valid;
    }

    for (auto& kv : lnkd_values_by_pair)
    {
        std::vector<double>& values = kv.second;
        if (values.empty()) continue;
        result.table.ln_kd_by_key[kv.first] = median_in_place(values);
        result.table.n_profiles_by_key[kv.first] = static_cast<int>(values.size());
    }

    result.pairs_built = result.table.ln_kd_by_key.size();
    return result;
}

sti_spr_parse_result load_sti_spr_table_from_csv(const std::string& csv_path, double eps)
{
    std::ifstream ifs(csv_path);
    if (!ifs.is_open())
    {
        throw std::runtime_error("load_sti_spr_table_from_csv: could not open: " + csv_path);
    }

    std::string header_line;
    if (!std::getline(ifs, header_line))
    {
        throw std::runtime_error("load_sti_spr_table_from_csv: empty CSV: " + csv_path);
    }

    const auto headers = parse_csv_line(header_line);
    const auto hmap = header_index_map(headers);
    if (hmap.find("SERCDE") == hmap.end() || hmap.find("SPR_PRC") == hmap.end())
    {
        throw std::runtime_error("load_sti_spr_table_from_csv: expected columns SERCDE and SPR_prc");
    }

    sti_spr_parse_result result;
    std::string line;
    while (std::getline(ifs, line))
    {
        if (line.empty()) continue;
        ++result.rows_read;
        const auto row = parse_csv_line(line);

        const std::string series_code = normalize_sti_series_code(fetch_cell(row, hmap, "SERCDE"));
        if (series_code.empty()) continue;

        double spr_prc = 0.0;
        if (!parse_double_safe(fetch_cell(row, hmap, "SPR_PRC"), spr_prc)) continue;

        double spr = spr_prc / 100.0;
        bool clamped = false;
        if (spr <= eps)
        {
            spr = eps;
            clamped = true;
        }
        if (spr >= 1.0 - eps)
        {
            spr = 1.0 - eps;
            clamped = true;
        }

        result.table.spr_fraction_by_series[series_code] = spr;
        ++result.rows_used;
        if (clamped) ++result.clamped_rows;
    }

    return result;
}

bool sti_check_xlsx_available(const std::string& xlsx_path)
{
    std::ifstream ifs(xlsx_path, std::ios::binary);
    return ifs.good();
}

sti_polygon_string_rasterise_result rasterise_polygon_string_field(
    const lsdtt_raster_cpu& reference_raster,
    const std::string& vector_path,
    const std::string& layer_name,
    const std::string& field_name,
    int oversample)
{
    if (oversample <= 0)
    {
        throw std::runtime_error("rasterise_polygon_string_field: oversample must be > 0");
    }

    const lsdtt_raster_meta& meta = reference_raster.get_meta();
    if (!meta.is_valid())
    {
        throw std::runtime_error("rasterise_polygon_string_field: reference raster metadata is invalid");
    }

    lsdtt::gdal_init();

    GDALDataset* vec_ds = static_cast<GDALDataset*>(
        GDALOpenEx(vector_path.c_str(), GDAL_OF_VECTOR | GDAL_OF_READONLY, nullptr, nullptr, nullptr));
    if (!vec_ds)
    {
        throw std::runtime_error("rasterise_polygon_string_field: cannot open vector file: " + vector_path);
    }

    OGRLayer* layer = vec_ds->GetLayerByName(layer_name.c_str());
    if (!layer)
    {
        GDALClose(vec_ds);
        throw std::runtime_error("rasterise_polygon_string_field: layer '" + layer_name + "' not found in " + vector_path);
    }

    if (!is_polygon_type(layer->GetGeomType()))
    {
        GDALClose(vec_ds);
        throw std::runtime_error("rasterise_polygon_string_field: layer is not polygon geometry");
    }

    set_spatial_filter_from_raster_bbox(layer, meta);
    OGRFeatureDefn* defn = layer->GetLayerDefn();
    if (!defn)
    {
        GDALClose(vec_ds);
        throw std::runtime_error("rasterise_polygon_string_field: invalid layer schema");
    }
    const int field_index = defn->GetFieldIndex(field_name.c_str());
    if (field_index < 0)
    {
        GDALClose(vec_ds);
        throw std::runtime_error("rasterise_polygon_string_field: field '" + field_name + "' not found in layer '" + layer_name + "'");
    }

    auto ct = make_ct_if_needed(layer->GetSpatialRef(), meta.crs_wkt);

    std::set<std::string> unique_values;
    layer->ResetReading();
    OGRFeature* feat = nullptr;
    while ((feat = layer->GetNextFeature()) != nullptr)
    {
        OGRGeometry* geom = feat->GetGeometryRef();
        if (!geom || !is_polygon_type(geom->getGeometryType()))
        {
            OGRFeature::DestroyFeature(feat);
            continue;
        }
        if (feat->IsFieldSetAndNotNull(field_index))
        {
            unique_values.insert(feat->GetFieldAsString(field_index));
        }
        OGRFeature::DestroyFeature(feat);
    }

    std::vector<std::pair<std::uint32_t, std::string>> lookup;
    std::unordered_map<std::string, std::uint32_t> lookup_map;
    lookup.reserve(unique_values.size());
    lookup_map.reserve(unique_values.size());

    std::uint32_t code = 1;
    for (const std::string& value : unique_values)
    {
        lookup.emplace_back(code, value);
        lookup_map.emplace(value, code);
        ++code;
    }

    const std::size_t npix = static_cast<std::size_t>(meta.n_rows) * meta.n_cols;
    std::vector<float> data(npix, 0.0f);
    std::vector<uint8_t> mask(npix, 0u);
    std::vector<int> best_count(npix, -1);

    const auto& ref_mask = reference_raster.get_mask();
    for (std::size_t i = 0; i < npix; ++i)
    {
        mask[i] = ref_mask[i] ? 1u : 0u;
    }

    layer->ResetReading();
    while ((feat = layer->GetNextFeature()) != nullptr)
    {
        OGRGeometry* geom = feat->GetGeometryRef();
        if (!geom || !is_polygon_type(geom->getGeometryType()))
        {
            OGRFeature::DestroyFeature(feat);
            continue;
        }
        if (!feat->IsFieldSetAndNotNull(field_index))
        {
            OGRFeature::DestroyFeature(feat);
            continue;
        }

        const std::string value = feat->GetFieldAsString(field_index);
        auto it_code = lookup_map.find(value);
        if (it_code == lookup_map.end())
        {
            OGRFeature::DestroyFeature(feat);
            continue;
        }

        OGRGeometry* clone = geom->clone();
        if (ct)
        {
            if (clone->transform(ct.get()) != OGRERR_NONE)
            {
                delete clone;
                OGRFeature::DestroyFeature(feat);
                continue;
            }
        }

        OGREnvelope env;
        clone->getEnvelope(&env);

        int row0 = 0;
        int row1 = 0;
        int col0 = 0;
        int col1 = 0;
        if (!raster_window_from_envelope(meta, env, row0, row1, col0, col1))
        {
            delete clone;
            OGRFeature::DestroyFeature(feat);
            continue;
        }

        const int coarse_rows = row1 - row0 + 1;
        const int coarse_cols = col1 - col0 + 1;
        const int burn_rows = coarse_rows * oversample;
        const int burn_cols = coarse_cols * oversample;

        double gt[6];
        fill_window_geotransform(meta, row0, col0, oversample, gt);
        GDALDataset* mem_ds = make_mem_byte_raster(burn_rows, burn_cols, gt);

        std::vector<OGRGeometryH> handles(1, reinterpret_cast<OGRGeometryH>(clone));
        int band_no = 1;
        double burn_val = 1.0;
        const CPLErr err = GDALRasterizeGeometries(
            mem_ds,
            1,
            &band_no,
            1,
            handles.data(),
            nullptr,
            nullptr,
            &burn_val,
            nullptr,
            nullptr,
            nullptr);
        if (err != CE_None)
        {
            GDALClose(mem_ds);
            delete clone;
            OGRFeature::DestroyFeature(feat);
            GDALClose(vec_ds);
            throw std::runtime_error("rasterise_polygon_string_field: GDALRasterizeGeometries failed");
        }

        std::vector<std::uint8_t> hi_res(static_cast<std::size_t>(burn_rows) * burn_cols, 0u);
        if (mem_ds->GetRasterBand(1)->RasterIO(
                GF_Read, 0, 0, burn_cols, burn_rows,
                hi_res.data(), burn_cols, burn_rows, GDT_Byte, 0, 0) != CE_None)
        {
            GDALClose(mem_ds);
            delete clone;
            OGRFeature::DestroyFeature(feat);
            GDALClose(vec_ds);
            throw std::runtime_error("rasterise_polygon_string_field: failed reading temporary burn raster");
        }

        GDALClose(mem_ds);
        delete clone;

        for (int wr = 0; wr < coarse_rows; ++wr)
        {
            const int global_r = row0 + wr;
            for (int wc = 0; wc < coarse_cols; ++wc)
            {
                const int global_c = col0 + wc;
                const std::size_t idx = static_cast<std::size_t>(global_r) * meta.n_cols + global_c;
                if (!ref_mask[idx]) continue;

                int count = 0;
                for (int sr = 0; sr < oversample; ++sr)
                {
                    const int hi_r = wr * oversample + sr;
                    const std::size_t row_base = static_cast<std::size_t>(hi_r) * burn_cols;
                    for (int sc = 0; sc < oversample; ++sc)
                    {
                        const int hi_c = wc * oversample + sc;
                        count += static_cast<int>(hi_res[row_base + hi_c]);
                    }
                }

                if (count > best_count[idx])
                {
                    best_count[idx] = count;
                    data[idx] = static_cast<float>(it_code->second);
                }
            }
        }

        OGRFeature::DestroyFeature(feat);
    }

    layer->SetSpatialFilter(nullptr);
    GDALClose(vec_ds);

    sti_polygon_string_rasterise_result result;
    result.raster_codes = lsdtt_raster_cpu(std::move(data), std::move(mask), meta);
    result.lookup = std::move(lookup);
    return result;
}

sti_compute_result compute_sti_from_rasters(
    const lsdtt_raster_cpu& twi_raster,
    const sti_polygon_string_rasterise_result& series_code_raster,
    const sti_polygon_string_rasterise_result& land_use_raster,
    const sti_lnkd_table& lnkd_table,
    const sti_spr_table& spr_table,
    const sti_compute_options& options)
{
    const lsdtt_raster_meta& meta = twi_raster.get_meta();
    const int n_rows = meta.n_rows;
    const int n_cols = meta.n_cols;
    const float nodata = static_cast<float>(meta.no_data_value);

    lsdtt_array2d_f sti_data(n_rows, n_cols, nodata);
    lsdtt_array2d_u8 sti_mask(n_rows, n_cols, uint8_t(0));
    lsdtt_array2d_f term_data(n_rows, n_cols, nodata);
    lsdtt_array2d_u8 term_mask(n_rows, n_cols, uint8_t(0));
    lsdtt_array2d_f valid_data(n_rows, n_cols, 0.0f);
    lsdtt_array2d_u8 valid_mask(n_rows, n_cols, uint8_t(1));

    std::unordered_map<std::uint32_t, std::string> series_lookup;
    std::unordered_map<std::uint32_t, std::string> land_lookup;
    series_lookup.reserve(series_code_raster.lookup.size());
    land_lookup.reserve(land_use_raster.lookup.size());
    for (const auto& kv : series_code_raster.lookup)
    {
        series_lookup.emplace(kv.first, normalize_sti_series_code(kv.second));
    }
    for (const auto& kv : land_use_raster.lookup)
    {
        land_lookup.emplace(kv.first, normalize_sti_land_use(kv.second));
    }

    const auto& twi_data = twi_raster.get_data_array();
    const auto& twi_mask = twi_raster.get_mask_array();
    const auto& series_data = series_code_raster.raster_codes.get_data_array();
    const auto& series_mask = series_code_raster.raster_codes.get_mask_array();
    const auto& land_data = land_use_raster.raster_codes.get_data_array();
    const auto& land_mask = land_use_raster.raster_codes.get_mask_array();

    sti_compute_result result;

    for (int row = 0; row < n_rows; ++row)
    {
        for (int col = 0; col < n_cols; ++col)
        {
            if (twi_mask(row, col) == 0) continue;
            ++result.stats.valid_twi_pixels;

            if (series_mask(row, col) == 0)
            {
                ++result.stats.missing_series_code_pixels;
                continue;
            }
            if (land_mask(row, col) == 0)
            {
                ++result.stats.missing_land_use_pixels;
                continue;
            }

            const std::uint32_t series_code = static_cast<std::uint32_t>(std::lround(series_data(row, col)));
            const std::uint32_t land_code = static_cast<std::uint32_t>(std::lround(land_data(row, col)));
            if (series_code == 0)
            {
                ++result.stats.missing_series_code_pixels;
                continue;
            }
            if (land_code == 0)
            {
                ++result.stats.missing_land_use_pixels;
                continue;
            }

            const auto it_series = series_lookup.find(series_code);
            const auto it_land = land_lookup.find(land_code);
            if (it_series == series_lookup.end() || it_series->second.empty())
            {
                ++result.stats.missing_series_code_pixels;
                continue;
            }
            if (it_land == land_lookup.end() || it_land->second.empty())
            {
                ++result.stats.missing_land_use_pixels;
                continue;
            }

            const std::string& series_str = it_series->second;
            const std::string& land_str = it_land->second;
            const std::string pair_key = make_pair_key(series_str, land_str);

            double term = 0.0;
            std::string source = "missing";
            bool has_term = false;

            if (options.mode == sti_term_mode::ksd)
            {
                int n_profiles = 0;
                has_term = lookup_ln_kd_with_fallback(
                    lnkd_table,
                    series_str,
                    land_str,
                    term,
                    source,
                    n_profiles);
                (void)n_profiles;
                if (source == "fallback_other_landuse") ++result.stats.fallback_pixels;
            }
            else
            {
                auto it_spr = spr_table.spr_fraction_by_series.find(series_str);
                if (it_spr != spr_table.spr_fraction_by_series.end())
                {
                    bool clamped = false;
                    term = compute_spr_term(it_spr->second, options.spr_transform, options.eps, clamped);
                    if (clamped) ++result.stats.clamped_pixels;
                    has_term = true;
                    source = "direct";
                }
            }

            result.source_by_series_landuse[pair_key] = source;

            if (!has_term)
            {
                ++result.stats.missing_lookup_pixels;
                continue;
            }

            const float twi = twi_data(row, col);
            sti_data(row, col) = twi - static_cast<float>(term);
            sti_mask(row, col) = uint8_t(1);
            term_data(row, col) = static_cast<float>(term);
            term_mask(row, col) = uint8_t(1);
            valid_data(row, col) = 1.0f;
            ++result.stats.valid_output_pixels;
        }
    }

    result.sti_raster = lsdtt_raster_cpu(std::move(sti_data), std::move(sti_mask), meta);
    result.term_raster = lsdtt_raster_cpu(std::move(term_data), std::move(term_mask), meta);
    result.valid_mask_raster = lsdtt_raster_cpu(std::move(valid_data), std::move(valid_mask), meta);
    return result;
}

std::vector<sti_lookup_row> build_sti_lookup_rows(
    const std::unordered_map<std::string, std::string>& source_by_series_landuse,
    const sti_lnkd_table& lnkd_table,
    const sti_spr_table& spr_table,
    sti_spr_transform spr_transform,
    double eps)
{
    std::vector<sti_lookup_row> rows;
    rows.reserve(source_by_series_landuse.size());

    for (const auto& kv : source_by_series_landuse)
    {
        const auto split = split_pair_key(kv.first);
        sti_lookup_row row;
        row.series_code = split.first;
        row.land_use = split.second;
        row.source = kv.second;

        auto it_ln = lnkd_table.ln_kd_by_key.find(kv.first);
        if (it_ln != lnkd_table.ln_kd_by_key.end())
        {
            row.has_ln_kd = true;
            row.ln_kd_median = it_ln->second;
            auto it_n = lnkd_table.n_profiles_by_key.find(kv.first);
            if (it_n != lnkd_table.n_profiles_by_key.end()) row.n_profiles = it_n->second;
        }
        else
        {
            double fallback_ln = 0.0;
            int n_profiles = 0;
            std::string fallback_source;
            if (lookup_ln_kd_with_fallback(lnkd_table, row.series_code, row.land_use,
                                           fallback_ln, fallback_source, n_profiles))
            {
                row.has_ln_kd = true;
                row.ln_kd_median = fallback_ln;
                row.n_profiles = n_profiles;
            }
        }

        auto it_spr = spr_table.spr_fraction_by_series.find(row.series_code);
        if (it_spr != spr_table.spr_fraction_by_series.end())
        {
            row.has_spr_prc = true;
            row.spr_prc = it_spr->second * 100.0;
            bool clamped = false;
            row.spr_term = compute_spr_term(it_spr->second, spr_transform, eps, clamped);
            row.has_spr_term = true;
        }

        rows.push_back(std::move(row));
    }

    std::sort(rows.begin(), rows.end(),
              [](const sti_lookup_row& a, const sti_lookup_row& b)
              {
                  if (a.series_code != b.series_code) return a.series_code < b.series_code;
                  return a.land_use < b.land_use;
              });

    return rows;
}

void write_sti_lookup_csv(
    const std::vector<sti_lookup_row>& rows,
    const std::string& out_csv_path)
{
    std::ofstream ofs(out_csv_path);
    if (!ofs.is_open())
    {
        throw std::runtime_error("write_sti_lookup_csv: could not open CSV for writing: " + out_csv_path);
    }

    ofs << "series_code,land_use,n_profiles,ln_kd_median,spr_prc,spr_term,source\n";
    for (const auto& row : rows)
    {
        csv_write_escaped(ofs, row.series_code);
        ofs << ',';
        csv_write_escaped(ofs, row.land_use);
        ofs << ',';
        ofs << row.n_profiles << ',';
        if (row.has_ln_kd) ofs << row.ln_kd_median;
        ofs << ',';
        if (row.has_spr_prc) ofs << row.spr_prc;
        ofs << ',';
        if (row.has_spr_term) ofs << row.spr_term;
        ofs << ',';
        csv_write_escaped(ofs, row.source);
        ofs << '\n';
    }
}

} // namespace lsdtt
