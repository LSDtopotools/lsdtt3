#include <cmath>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>

#include <cpu/lsdtt_connectivity_cpu.hpp>
#include <cpu/lsdtt_sti_cpu.hpp>

#include <core/lsdtt_raster_warp.hpp>

#include "binding_common.hpp"

namespace lsdtt3_python {
namespace {

Raster as_raster(lsdtt::lsdtt_raster_result result)
{
    return Raster(lsdtt::lsdtt_raster_cpu(
        std::move(result.data), std::move(result.mask), result.meta));
}

bool grids_match(const lsdtt::lsdtt_raster_meta& a, const lsdtt::lsdtt_raster_meta& b)
{
    return a.n_rows == b.n_rows && a.n_cols == b.n_cols &&
        a.x_min == b.x_min && a.y_min == b.y_min &&
        a.cell_size == b.cell_size && lsdtt::crs_equal_semantically(a.crs_wkt, b.crs_wkt);
}

void require_same_grid(const Raster& a, const Raster& b, const char* name)
{
    if (!grids_match(a.native().get_meta(), b.native().get_meta()))
        throw std::invalid_argument(std::string(name) + " must share the reference raster grid");
}

void require_flow_grid(const FlowInfo& flow, const Raster& raster, const char* name)
{
    if (!grids_match(flow.native().meta(), raster.native().get_meta()))
        throw std::invalid_argument(std::string(name) + " must share the FlowInfo grid");
}

lsdtt::lsdtt_raster_result as_result(const Raster& raster)
{
    return {raster.native().get_data(), raster.native().get_mask(), raster.native().get_meta()};
}

Raster topographic_wetness_index(const Raster& drainage_area, const Raster& slope)
{
    require_same_grid(drainage_area, slope, "slope");
    const auto& area = drainage_area.native();
    const auto& gradient = slope.native();
    const auto& meta = area.get_meta();
    const auto& a = area.get_data(); const auto& am = area.get_mask();
    const auto& s = gradient.get_data(); const auto& sm = gradient.get_mask();
    const double cell_area = meta.cell_size * meta.cell_size;
    if (cell_area <= 0.0) throw std::invalid_argument("drainage_area must have a positive cell size");
    std::vector<float> data(a.size(), static_cast<float>(meta.no_data_value));
    std::vector<std::uint8_t> mask(a.size(), 0u);
    for (std::size_t i = 0; i < a.size(); ++i)
    {
        if (am[i] == 0 || sm[i] == 0 || a[i] <= 0.0f || s[i] <= 0.0f) continue;
        data[i] = static_cast<float>(std::log((static_cast<double>(a[i]) / cell_area) / s[i]));
        mask[i] = 1u;
    }
    return Raster(lsdtt::lsdtt_raster_cpu(std::move(data), std::move(mask), meta));
}

class ConnectivityResult {
public:
    ConnectivityResult(lsdtt::lsdtt_raster_result index, lsdtt::lsdtt_raster_result weight,
                       lsdtt::lsdtt_raster_result upstream, lsdtt::lsdtt_raster_result downstream)
        : index_(as_raster(std::move(index))), weight_(as_raster(std::move(weight))),
          upstream_(as_raster(std::move(upstream))), downstream_(as_raster(std::move(downstream))) {}
    const Raster& index() const noexcept { return index_; }
    const Raster& weight() const noexcept { return weight_; }
    const Raster& upstream() const noexcept { return upstream_; }
    const Raster& downstream() const noexcept { return downstream_; }
private:
    Raster index_; Raster weight_; Raster upstream_; Raster downstream_;
};

ConnectivityResult sediment_connectivity_index(const FlowInfo& flow, const Raster& roughness,
                                                const Raster& slope, const Raster& drainage_area,
                                                float area_threshold_m2)
{
    if (area_threshold_m2 <= 0.0f) throw std::invalid_argument("area_threshold_m2 must be positive");
    require_flow_grid(flow, roughness, "roughness"); require_flow_grid(flow, slope, "slope");
    require_flow_grid(flow, drainage_area, "drainage_area");
    const auto weight = lsdtt::compute_W_raster(roughness.native());
    const auto area = as_result(drainage_area);
    const auto upstream = lsdtt::compute_D_up_raster(flow.native(), weight, slope.native(), area);
    const auto downstream = lsdtt::compute_D_down_raster(flow.native(), weight, slope.native(), area, area_threshold_m2);
    const auto index = lsdtt::compute_connectivity_index(flow.native(), roughness.native(), slope.native(), area, area_threshold_m2);
    return ConnectivityResult(index, weight, upstream, downstream);
}

class STILnkdTable {
public:
    explicit STILnkdTable(lsdtt::sti_lnkd_build_result result)
        : table_(std::move(result.table)), rows_read_(result.rows_read), rows_used_(result.rows_used),
          profiles_total_(result.profiles_total), profiles_valid_(result.profiles_valid), pairs_built_(result.pairs_built) {}
    const lsdtt::sti_lnkd_table& native() const noexcept { return table_; }
    nb::dict summary() const {
        nb::dict out; out["rows_read"] = rows_read_; out["rows_used"] = rows_used_;
        out["profiles_total"] = profiles_total_; out["profiles_valid"] = profiles_valid_; out["pairs_built"] = pairs_built_; return out;
    }
private:
    lsdtt::sti_lnkd_table table_; std::size_t rows_read_, rows_used_, profiles_total_, profiles_valid_, pairs_built_;
};

class STISprTable {
public:
    explicit STISprTable(lsdtt::sti_spr_parse_result result)
        : table_(std::move(result.table)), rows_read_(result.rows_read), rows_used_(result.rows_used), clamped_rows_(result.clamped_rows) {}
    const lsdtt::sti_spr_table& native() const noexcept { return table_; }
    nb::dict summary() const { nb::dict out; out["rows_read"] = rows_read_; out["rows_used"] = rows_used_; out["clamped_rows"] = clamped_rows_; return out; }
private:
    lsdtt::sti_spr_table table_; std::size_t rows_read_, rows_used_, clamped_rows_;
};

class STICategoricalRaster {
public:
    explicit STICategoricalRaster(lsdtt::sti_polygon_string_rasterise_result result) : result_(std::move(result)) {}
    Raster codes() const { return Raster(result_.raster_codes); }
    nb::object lookup() const {
        nb::object pandas = nb::module_::import_("pandas"); nb::list rows;
        for (const auto& item : result_.lookup) { nb::dict row; row["code"] = item.first; row["value"] = item.second; rows.append(std::move(row)); }
        return pandas.attr("DataFrame")(rows);
    }
    const lsdtt::sti_polygon_string_rasterise_result& native() const noexcept { return result_; }
private: lsdtt::sti_polygon_string_rasterise_result result_;
};

class SoilTopographicIndexResult {
public:
    SoilTopographicIndexResult(lsdtt::sti_compute_result result, std::vector<lsdtt::sti_lookup_row> lookup)
        : sti_(std::move(result.sti_raster)), term_(std::move(result.term_raster)), valid_mask_(std::move(result.valid_mask_raster)), stats_(result.stats), lookup_(std::move(lookup)) {}
    const Raster& sti() const noexcept { return sti_; } const Raster& term() const noexcept { return term_; }
    const Raster& valid_mask() const noexcept { return valid_mask_; }
    nb::dict stats() const {
        nb::dict out; out["valid_twi_pixels"] = stats_.valid_twi_pixels; out["valid_output_pixels"] = stats_.valid_output_pixels;
        out["missing_series_code_pixels"] = stats_.missing_series_code_pixels; out["missing_land_use_pixels"] = stats_.missing_land_use_pixels;
        out["missing_lookup_pixels"] = stats_.missing_lookup_pixels; out["fallback_pixels"] = stats_.fallback_pixels; out["clamped_pixels"] = stats_.clamped_pixels; return out;
    }
    nb::object lookup() const {
        nb::object pandas = nb::module_::import_("pandas"); nb::list rows;
        for (const auto& item : lookup_) { nb::dict row; row["series_code"] = item.series_code; row["land_use"] = item.land_use; row["n_profiles"] = item.n_profiles; row["has_ln_kd"] = item.has_ln_kd; row["ln_kd_median"] = item.ln_kd_median; row["has_spr_prc"] = item.has_spr_prc; row["spr_prc"] = item.spr_prc; row["has_spr_term"] = item.has_spr_term; row["spr_term"] = item.spr_term; row["source"] = item.source; rows.append(std::move(row)); }
        return pandas.attr("DataFrame")(rows);
    }
private: Raster sti_; Raster term_; Raster valid_mask_; lsdtt::sti_compute_stats stats_; std::vector<lsdtt::sti_lookup_row> lookup_;
};

STILnkdTable load_sti_lnkd_table(const std::string& path, double depth_cap_cm, double min_valid_depth_cm, double ks_unit_scale_to_m_per_day, double eps)
{
    if (depth_cap_cm <= 0.0 || min_valid_depth_cm <= 0.0 || ks_unit_scale_to_m_per_day <= 0.0 || eps <= 0.0) throw std::invalid_argument("STI table options must be positive");
    lsdtt::sti_lnkd_options options; options.depth_cap_cm = depth_cap_cm; options.min_valid_depth_cm = min_valid_depth_cm; options.ks_unit_scale_to_m_per_day = ks_unit_scale_to_m_per_day; options.eps = eps;
    return STILnkdTable(lsdtt::build_sti_lnkd_table_from_csv(path, options));
}

STISprTable load_sti_spr_table(const std::string& path, double eps)
{
    if (eps <= 0.0) throw std::invalid_argument("eps must be positive");
    return STISprTable(lsdtt::load_sti_spr_table_from_csv(path, eps));
}

STICategoricalRaster rasterize_sti_categories(const Raster& reference, const std::string& vector_path, const std::string& field_name, const std::string& layer_name, int oversample)
{
    if (field_name.empty()) throw std::invalid_argument("field_name must not be empty");
    if (oversample <= 0) throw std::invalid_argument("oversample must be positive");
    return STICategoricalRaster(lsdtt::rasterise_polygon_string_field(reference.native(), vector_path, layer_name, field_name, oversample));
}

SoilTopographicIndexResult soil_topographic_index(const Raster& twi, const STICategoricalRaster& series_codes, const STICategoricalRaster& land_use, const STILnkdTable& lnkd, const STISprTable& spr, const std::string& mode, const std::string& spr_transform, double eps)
{
    if (eps <= 0.0) throw std::invalid_argument("eps must be positive");
    lsdtt::sti_compute_options options; options.eps = eps;
    if (mode == "ksd") options.mode = lsdtt::sti_term_mode::ksd; else if (mode == "spr") options.mode = lsdtt::sti_term_mode::spr; else throw std::invalid_argument("mode must be 'ksd' or 'spr'");
    if (spr_transform == "ln_1_minus_spr") options.spr_transform = lsdtt::sti_spr_transform::ln_1_minus_spr; else if (spr_transform == "ln_spr") options.spr_transform = lsdtt::sti_spr_transform::ln_spr; else throw std::invalid_argument("spr_transform must be 'ln_1_minus_spr' or 'ln_spr'");
    auto result = lsdtt::compute_sti_from_rasters(twi.native(), series_codes.native(), land_use.native(), lnkd.native(), spr.native(), options);
    const auto lookup = lsdtt::build_sti_lookup_rows(result.source_by_series_landuse, lnkd.native(), spr.native(), options.spr_transform, eps);
    return SoilTopographicIndexResult(std::move(result), lookup);
}

} // namespace

void bind_flow_indices(nb::module_& module)
{
    module.def("topographic_wetness_index", &topographic_wetness_index, nb::arg("drainage_area"), nb::arg("slope"), R"doc(Compute ``log((A / cell_area) / slope)`` from drainage area and slope rasters.

Parameters
----------
drainage_area : Raster
    Drainage area (m^2), for example ``FlowInfo.drainage_area()`` or
    ``FlowGraph.drainage_area()``.
slope : Raster
    Gradient (m/m) on the same grid, for example ``slope(dem)``.

Returns
-------
Raster
    Natural-log wetness index; nodata where area or slope is not positive.
)doc");
    nb::class_<ConnectivityResult>(module, "ConnectivityResult", "Sediment-connectivity rasters following Cavalli et al. (2013).")
        .def_prop_ro("index", &ConnectivityResult::index, nb::rv_policy::reference_internal).def_prop_ro("weight", &ConnectivityResult::weight, nb::rv_policy::reference_internal)
        .def_prop_ro("upstream", &ConnectivityResult::upstream, nb::rv_policy::reference_internal).def_prop_ro("downstream", &ConnectivityResult::downstream, nb::rv_policy::reference_internal);
    module.def("sediment_connectivity_index", &sediment_connectivity_index, nb::arg("flow"), nb::arg("roughness"), nb::arg("slope"), nb::arg("drainage_area"), nb::arg("area_threshold_m2"), R"doc(Compute sediment connectivity using D8 topology; drainage area may be D8, MFD, or D-infinity.

Cavalli et al. (2013) index ``IC = log10(D_up / D_down)`` with
``D_up = mean(W) * mean(S) * sqrt(A)`` over the upslope area and
``D_down`` the sum of ``d / (W S)`` along the D8 path to the nearest sink;
``W = max(1 - R / R_max, 0.001)`` from the roughness R.

Parameters
----------
flow : FlowInfo
    D8 flow topology.
roughness : Raster
    Surface roughness on the flow grid (for example a local standard
    deviation of residual elevation, m).
slope : Raster
    Gradient (m/m); values below 0.001 are clipped in ``D_down``.
drainage_area : Raster
    Drainage area (m^2) from any routing scheme.
area_threshold_m2 : float
    Pixels with at least this drainage area (m^2) are sinks (channels).

Returns
-------
ConnectivityResult
    ``index``, ``weight`` (W), ``upstream`` (D_up) and ``downstream``
    (D_down) Rasters.
)doc");
    nb::class_<STILnkdTable>(module, "STILnkdTable", "Parsed soil hydraulic-conductivity lookup table.").def_prop_ro("summary", &STILnkdTable::summary);
    nb::class_<STISprTable>(module, "STISprTable", "Parsed soil-profile runoff lookup table.").def_prop_ro("summary", &STISprTable::summary);
    nb::class_<STICategoricalRaster>(module, "STICategoricalRaster", "Rasterized categorical STI field and code lookup.").def_prop_ro("codes", &STICategoricalRaster::codes).def_prop_ro("lookup", &STICategoricalRaster::lookup);
    nb::class_<SoilTopographicIndexResult>(module, "SoilTopographicIndexResult", "Soil-topographic-index outputs and diagnostics.")
        .def_prop_ro("sti", &SoilTopographicIndexResult::sti, nb::rv_policy::reference_internal).def_prop_ro("term", &SoilTopographicIndexResult::term, nb::rv_policy::reference_internal)
        .def_prop_ro("valid_mask", &SoilTopographicIndexResult::valid_mask, nb::rv_policy::reference_internal).def_prop_ro("stats", &SoilTopographicIndexResult::stats).def_prop_ro("lookup", &SoilTopographicIndexResult::lookup);
    module.def("load_sti_lnkd_table", &load_sti_lnkd_table, nb::arg("csv_path"), nb::arg("depth_cap_cm") = 100.0, nb::arg("min_valid_depth_cm") = 30.0, nb::arg("ks_unit_scale_to_m_per_day") = 1.0, nb::arg("eps") = 1.0e-12, R"doc(Load the KSD STI lookup table from a soil-profile CSV.

Parameters
----------
csv_path : str
    Soil-profile CSV with saturated hydraulic conductivity by horizon.
depth_cap_cm : float, default=100.0
    Depth (cm) to which profiles are integrated.
min_valid_depth_cm : float, default=30.0
    Profiles shallower than this (cm) are rejected.
ks_unit_scale_to_m_per_day : float, default=1.0
    Factor converting the CSV conductivity to m/day.
eps : float, default=1e-12
    Small positive floor that avoids log(0).

Returns
-------
STILnkdTable
)doc");
    module.def("load_sti_spr_table", &load_sti_spr_table, nb::arg("csv_path"), nb::arg("eps") = 1.0e-12, R"doc(Load the SPR STI lookup table from a CSV.

Parameters
----------
csv_path : str
    CSV of standard percentage runoff (SPR) values.
eps : float, default=1e-12
    Small positive floor that avoids log(0).

Returns
-------
STISprTable
)doc");
    module.def("rasterize_sti_categories", &rasterize_sti_categories, nb::arg("reference"), nb::arg("vector_path"), nb::arg("field_name"), nb::arg("layer_name"), nb::arg("oversample") = 8, R"doc(Rasterize a string vector attribute and retain its categorical lookup.

Parameters
----------
reference : Raster
    Raster whose grid is used.
vector_path : str
    Polygon vector file.
field_name : str
    Attribute holding the category (for example a soil series code).
layer_name : str
    Layer in the vector file.
oversample : int, default=8
    Sub-pixel oversampling used to pick each pixel's majority category.

Returns
-------
STICategoricalRaster
    Integer codes plus the code-to-category lookup.
)doc");
    module.def("soil_topographic_index", &soil_topographic_index, nb::arg("twi"), nb::arg("series_codes"), nb::arg("land_use"), nb::arg("lnkd_table"), nb::arg("spr_table"), nb::arg("mode") = "ksd", nb::arg("spr_transform") = "ln_1_minus_spr", nb::arg("eps") = 1.0e-12, R"doc(Compute STI from TWI, categorical rasters, and lookup tables.

Parameters
----------
twi : Raster
    Topographic wetness index (``topographic_wetness_index``).
series_codes : STICategoricalRaster
    Soil series categories (``rasterize_sti_categories``).
land_use : STICategoricalRaster
    Land-use categories on the same grid.
lnkd_table : STILnkdTable
    Conductivity lookup (``load_sti_lnkd_table``).
spr_table : STISprTable
    Runoff lookup (``load_sti_spr_table``).
mode : {"ksd", "spr"}, default="ksd"
    Soil term from conductivity x depth or from standard percentage runoff.
spr_transform : {"ln_1_minus_spr", "ln_spr"}, default="ln_1_minus_spr"
    Transform of the SPR term.
eps : float, default=1e-12
    Small positive floor that avoids log(0).

Returns
-------
SoilTopographicIndexResult
)doc");
}
} // namespace lsdtt3_python
