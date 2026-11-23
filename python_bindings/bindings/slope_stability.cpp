// slope_stability.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
//
// -----------------------------------------------------------------------------
// lsdtt-slope-stability's routines: SINMAP, scaLr rupture planes and their
// probability of failure, the landslide catalogues and frequency - area fit,
// slip surfaces (method of slices) and ellipsoids (method of columns), and the
// Dupuit - Boussinesq groundwater model. The same C++ calls as the CLI, with
// its parameter names less the rupture_ / slip_surface_ / groundwater_
// prefixes. Region parameters are keyed by the CLI's region-table column
// names. The CLI's routing-surface conditioning is done in Python with the
// bound conditioning routines.
// -----------------------------------------------------------------------------

#include <cmath>
#include <cstdint>
#include <limits>
#include <map>
#include <memory>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <nanobind/nanobind.h>
#include <nanobind/stl/map.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>

#include <cpu/lsdtt_groundwater_dupuit_cpu.hpp>
#include <cpu/lsdtt_slope_stability_catalogue_cpu.hpp>
#include <cpu/lsdtt_slope_stability_columns_cpu.hpp>
#include <cpu/lsdtt_slope_stability_dem_change_cpu.hpp>
#include <cpu/lsdtt_slope_stability_hazard_area_cpu.hpp>
#include <cpu/lsdtt_slope_stability_rupture_cpu.hpp>
#include <cpu/lsdtt_slope_stability_sinmap_cpu.hpp>
#include <cpu/lsdtt_slope_stability_slip_surface_cpu.hpp>

#include "binding_common.hpp"

namespace lsdtt3_python {
namespace {

using meta_t = lsdtt::lsdtt_raster_meta;

// ----------------------------------------------------------------- helpers
void require_grid(const meta_t& a, const meta_t& b, const std::string& what)
{
    if (a.n_rows != b.n_rows || a.n_cols != b.n_cols || a.x_min != b.x_min || a.y_min != b.y_min ||
        a.cell_size != b.cell_size)
        throw std::invalid_argument(what + " must share the grid");
}

const Raster& as_grid_raster(nb::handle obj, const meta_t& meta, const std::string& what)
{
    const Raster& r = nb::cast<const Raster&>(obj);
    require_grid(r.native().get_meta(), meta, what);
    return r;
}

// None -> empty; else 1 where valid and nonzero.
std::vector<std::uint8_t> mask_of(nb::object obj, const meta_t& meta, const std::string& what)
{
    std::vector<std::uint8_t> out;
    if (obj.is_none()) return out;
    const Raster& r = as_grid_raster(obj, meta, what);
    const auto& d = r.native().get_data();
    const auto& m = r.native().get_mask();
    out.assign(d.size(), 0);
    for (std::size_t i = 0; i < d.size(); ++i) out[i] = (m[i] && std::isfinite(d[i]) && d[i] != 0.0f) ? 1 : 0;
    return out;
}

// None -> empty; else the values, NaN at nodata.
std::vector<double> values_of(nb::object obj, const meta_t& meta, const std::string& what)
{
    std::vector<double> out;
    if (obj.is_none()) return out;
    const Raster& r = as_grid_raster(obj, meta, what);
    const auto& d = r.native().get_data();
    const auto& m = r.native().get_mask();
    out.assign(d.size(), std::numeric_limits<double>::quiet_NaN());
    for (std::size_t i = 0; i < d.size(); ++i)
        if (m[i] && std::isfinite(d[i])) out[i] = d[i];
    return out;
}

// float or Raster -> one value per cell (empty for a float, which the caller uses).
bool scalar_or_raster(nb::object obj, const meta_t& meta, const std::string& what, double& scalar,
                      std::vector<double>& raster)
{
    if (nb::isinstance<Raster>(obj)) { raster = values_of(obj, meta, what); return true; }
    scalar = nb::cast<double>(obj);
    return false;
}

Raster raster_of(const std::vector<double>& v, const meta_t& meta)
{
    std::vector<float> d(v.size(), static_cast<float>(meta.no_data_value));
    std::vector<std::uint8_t> m(v.size(), 0);
    for (std::size_t i = 0; i < v.size(); ++i)
        if (std::isfinite(v[i])) { d[i] = static_cast<float>(v[i]); m[i] = 1; }
    return Raster(lsdtt::lsdtt_raster_cpu(std::move(d), std::move(m), meta));
}

Raster label_raster(const std::vector<std::int32_t>& lab, const meta_t& meta)
{
    std::vector<double> v(lab.size(), std::numeric_limits<double>::quiet_NaN());
    for (std::size_t i = 0; i < lab.size(); ++i)
        if (lab[i] > 0) v[i] = lab[i];
    return raster_of(v, meta);
}

nb::object data_frame(const nb::dict& d) { return nb::module_::import_("pandas").attr("DataFrame")(d); }

nb::object numpy_int64(const std::vector<std::int64_t>& v)
{
    return nb::module_::import_("numpy").attr("asarray")(nb::cast(v), "int64");
}

// Region raster -> ids (k_*_region_nodata where nodata); integer check.
std::vector<std::int32_t> region_ids(nb::object region, const meta_t& meta)
{
    std::vector<std::int32_t> out;
    if (region.is_none()) return out;
    const Raster& r = as_grid_raster(region, meta, "region");
    const auto& d = r.native().get_data();
    const auto& m = r.native().get_mask();
    out.assign(d.size(), lsdtt::k_rupture_region_nodata);
    for (std::size_t i = 0; i < d.size(); ++i)
    {
        if (!m[i] || !std::isfinite(d[i])) continue;
        const double v = d[i];
        if (v != std::round(v)) throw std::invalid_argument("region: region ids must be integers");
        out[i] = static_cast<std::int32_t>(v);
    }
    return out;
}

using overrides_t = std::map<std::int32_t, std::map<std::string, double>>;

// region_parameters: {region_id: {CLI region-table column: value}}. Keys
// outside `allowed` (and region_name, which is ignored) are an error, as the
// CLI's region table refuses unknown columns.
overrides_t region_overrides(nb::object region, nb::object region_parameters, const std::set<std::string>& allowed)
{
    overrides_t out;
    if (region_parameters.is_none())
    {
        if (!region.is_none())
            throw std::invalid_argument("region needs region_parameters ({region_id: {parameter: value}})");
        out[0];
        return out;
    }
    if (region.is_none()) throw std::invalid_argument("region_parameters needs a region raster");
    for (auto item : nb::cast<nb::dict>(region_parameters))
    {
        const std::int32_t id = nb::cast<std::int32_t>(item.first);
        auto& row = out[id];
        for (auto kv : nb::cast<nb::dict>(item.second))
        {
            const std::string key = nb::cast<std::string>(kv.first);
            if (key == "region_name" || key == "region_id") continue;
            if (!allowed.count(key))
            {
                std::string list;
                for (const auto& a : allowed) list += (list.empty() ? "" : ", ") + a;
                throw std::invalid_argument("region_parameters[" + std::to_string(id) + "]: unknown parameter '" + key
                                            + "' (allowed: " + list + ")");
            }
            // None or NaN (an empty cell of a table read with pandas) = not set, as the
            // CLI's region table treats an empty cell.
            if (kv.second.is_none()) continue;
            const double v = nb::cast<double>(kv.second);
            if (std::isnan(v)) continue;
            row[key] = v;
        }
    }
    if (out.empty()) throw std::invalid_argument("region_parameters is empty");
    return out;
}

double value_or(const std::map<std::string, double>& row, const std::string& key, double fallback)
{
    const auto it = row.find(key);
    return it == row.end() ? fallback : it->second;
}

// ------------------------------------------------------------------- SINMAP
const std::set<std::string> k_sinmap_keys = {
    "sinmap_dimensionless_cohesion_min", "sinmap_dimensionless_cohesion_max", "sinmap_friction_coefficient_min",
    "sinmap_friction_coefficient_max", "sinmap_transmissivity_over_recharge_min_metres",
    "sinmap_transmissivity_over_recharge_max_metres", "sinmap_water_to_soil_density_ratio"};
// The CLI's region table may hold both sets; each routine ignores the other's.
const std::set<std::string> k_rupture_keys = {
    "rupture_cohesion_pa", "rupture_friction_coefficient", "rupture_density_kg_m3",
    "rupture_cohesion_pa_min", "rupture_cohesion_pa_max", "rupture_cohesion_pa_sd",
    "rupture_friction_coefficient_min", "rupture_friction_coefficient_max", "rupture_friction_coefficient_sd"};

std::set<std::string> all_region_keys()
{
    std::set<std::string> s = k_sinmap_keys;
    s.insert(k_rupture_keys.begin(), k_rupture_keys.end());
    return s;
}

lsdtt::sinmap_parameters sinmap_params(double c_min, double c_max, double t_min, double t_max, double tr_min,
                                       double tr_max, double r)
{
    lsdtt::sinmap_parameters p;
    p.dimensionless_cohesion_min = c_min;
    p.dimensionless_cohesion_max = c_max;
    p.friction_coefficient_min = t_min;
    p.friction_coefficient_max = t_max;
    p.transmissivity_over_recharge_min_metres = tr_min;
    p.transmissivity_over_recharge_max_metres = tr_max;
    p.water_to_soil_density_ratio = r;
    lsdtt::validate_sinmap_parameters(p);
    return p;
}

lsdtt::sinmap_parameters sinmap_with(const lsdtt::sinmap_parameters& b, const std::map<std::string, double>& row)
{
    return sinmap_params(value_or(row, "sinmap_dimensionless_cohesion_min", b.dimensionless_cohesion_min),
                         value_or(row, "sinmap_dimensionless_cohesion_max", b.dimensionless_cohesion_max),
                         value_or(row, "sinmap_friction_coefficient_min", b.friction_coefficient_min),
                         value_or(row, "sinmap_friction_coefficient_max", b.friction_coefficient_max),
                         value_or(row, "sinmap_transmissivity_over_recharge_min_metres",
                                  b.transmissivity_over_recharge_min_metres),
                         value_or(row, "sinmap_transmissivity_over_recharge_max_metres",
                                  b.transmissivity_over_recharge_max_metres),
                         value_or(row, "sinmap_water_to_soil_density_ratio", b.water_to_soil_density_ratio));
}

nb::dict sinmap(const Raster& slope, const Raster& drainage_area, double c_min, double c_max, double t_min,
                double t_max, double tr_min, double tr_max, double r, nb::object region, nb::object region_parameters,
                nb::object exclude, nb::object contour_width, int n_threads)
{
    const meta_t& meta = slope.native().get_meta();
    require_grid(drainage_area.native().get_meta(), meta, "drainage_area");
    const lsdtt::sinmap_parameters base = sinmap_params(c_min, c_max, t_min, t_max, tr_min, tr_max, r);
    const auto ov = region_overrides(region, region_parameters, all_region_keys());
    std::map<std::int32_t, lsdtt::sinmap_parameters> by_region;
    for (const auto& kv : ov) by_region[kv.first] = sinmap_with(base, kv.second);
    const auto ids = region_ids(region, meta);
    const auto ex = mask_of(exclude, meta, "exclude");
    std::vector<double> width = values_of(contour_width, meta, "contour_width");
    for (double& w : width)
        if (!std::isfinite(w)) w = meta.cell_size;
    const lsdtt::sinmap_rasters res = lsdtt::compute_sinmap(slope.native(), drainage_area.native(), ids, by_region, ex,
                                                            n_threads, width);
    nb::dict out;
    out["stability_index"] = Raster(res.stability_index);
    out["stability_class"] = Raster(res.stability_class);
    out["probability_stable"] = Raster(res.probability_stable);
    out["fs_min"] = Raster(res.fs_min);
    out["fs_max"] = Raster(res.fs_max);
    out["saturation"] = Raster(res.saturation);
    out["specific_catchment_area"] = Raster(res.specific_catchment_area);
    out["n_evaluated"] = res.n_evaluated;
    out["n_flat"] = res.n_flat;
    out["n_excluded"] = res.n_excluded;
    std::vector<std::int64_t> rid, cls, nc;
    std::vector<double> area, frac;
    for (const auto& row : lsdtt::sinmap_class_areas(res, ids))
    {
        rid.push_back(row.region_id);
        cls.push_back(row.stability_class);
        nc.push_back(row.n_cells);
        area.push_back(row.area_m2);
        frac.push_back(row.fraction_of_region);
    }
    nb::dict t;
    t["region_id"] = rid;
    t["stability_class"] = cls;
    t["n_cells"] = nc;
    t["area_m2"] = area;
    t["fraction_of_region"] = frac;
    out["class_areas"] = data_frame(t);
    return out;
}

nb::object sinmap_slope_area_curves(double slope_max, double area_min, double area_max, double c_min, double c_max,
                                    double t_min, double t_max, double tr_min, double tr_max, double r, int n_slopes)
{
    if (!(slope_max > 0.0) || !(area_min > 0.0) || !(area_max > area_min))
        throw std::invalid_argument("sinmap_slope_area_curves: need slope_max > 0 and 0 < area_min < area_max");
    if (n_slopes < 2) throw std::invalid_argument("sinmap_slope_area_curves: n_slopes must be at least 2");
    const auto pts = lsdtt::sinmap_slope_area_curves(sinmap_params(c_min, c_max, t_min, t_max, tr_min, tr_max, r),
                                                     slope_max, area_min, area_max, n_slopes);
    std::vector<std::string> curve;
    std::vector<double> s, a;
    for (const auto& p : pts)
    {
        curve.push_back(p.curve);
        s.push_back(p.slope);
        a.push_back(p.specific_catchment_area_metres);
    }
    nb::dict d;
    d["curve"] = curve;
    d["slope"] = s;
    d["specific_catchment_area_m"] = a;
    return data_frame(d);
}

std::vector<std::uint8_t> graph_valid(const lsdtt::lsdtt_flowgraph_cpu& g)
{
    const std::int64_t n = static_cast<std::int64_t>(g.meta().n_rows) * g.meta().n_cols;
    std::vector<std::uint8_t> v(static_cast<std::size_t>(n), 0);
    for (std::int64_t i = 0; i < n; ++i) v[static_cast<std::size_t>(i)] = g.is_valid_flat(i) ? 1 : 0;
    return v;
}

void require_dinf(const FlowGraph& g, const char* what)
{
    if (g.native().options().method != lsdtt::lsdtt_flowgraph_method::d_infinity)
        throw std::invalid_argument(std::string(what) + " needs a FlowGraph built with method='dinf'");
}

Raster dinf_contour_width(const FlowGraph& g)
{
    require_dinf(g, "dinf_contour_width");
    const auto valid = graph_valid(g.native());
    std::vector<double> w = lsdtt::dinf_contour_width(g.native());
    for (std::size_t i = 0; i < w.size(); ++i)
        if (!valid[i]) w[i] = std::numeric_limits<double>::quiet_NaN();
    return raster_of(w, g.native().meta());
}

Raster dinf_edge_contaminated_mask(const FlowGraph& g)
{
    require_dinf(g, "dinf_edge_contaminated_mask");
    const auto valid = graph_valid(g.native());
    const auto m = lsdtt::dinf_edge_contaminated_mask(g.native(), g.native().meta(), valid);
    std::vector<double> v(m.size(), std::numeric_limits<double>::quiet_NaN());
    for (std::size_t i = 0; i < m.size(); ++i)
        if (valid[i]) v[i] = m[i];
    return raster_of(v, g.native().meta());
}

// ------------------------------------------------------------ rupture planes
// Everything a rupture run was made from: shared by RupturePlanes and the
// slip surfaces built on it.
struct rupture_state
{
    std::shared_ptr<lsdtt::lsdtt_flowinfo> flow;
    lsdtt::lsdtt_raster_cpu geometry;
    std::vector<std::uint8_t> channel;
    std::vector<std::int32_t> region;
    overrides_t overrides;
    std::map<std::int32_t, lsdtt::rupture_parameters> by_region;
    std::vector<double> water_table_depth;
    lsdtt::rupture_rasters result;
    int n_threads = 0;
};

class RupturePlanes
{
public:
    explicit RupturePlanes(std::shared_ptr<rupture_state> s) : s_(std::move(s)) {}
    const rupture_state& state() const { return *s_; }
    std::shared_ptr<rupture_state> shared() const { return s_; }

private:
    std::shared_ptr<rupture_state> s_;
};

RupturePlanes rupture_planes(const FlowInfo& flow, const Raster& geometry, nb::object channel,
                             std::int64_t threshold_contributing_pixels, double cohesion_pa,
                             double friction_coefficient, double density_kg_m3, nb::object region,
                             nb::object region_parameters, nb::object water_table_depth, int n_threads)
{
    const meta_t& meta = flow.native().meta();
    require_grid(geometry.native().get_meta(), meta, "geometry");
    auto s = std::make_shared<rupture_state>();
    s->flow = flow.shared_native();
    s->geometry = geometry.native();
    s->n_threads = n_threads;
    if (channel.is_none())
    {
        if (threshold_contributing_pixels < 1)
            throw std::invalid_argument("threshold_contributing_pixels must be at least 1");
        const lsdtt::lsdtt_raster_result npix = flow.native().get_n_contributing_pixels_raster();
        s->channel.assign(npix.data.size(), 0);
        for (std::size_t i = 0; i < s->channel.size(); ++i)
            s->channel[i] = (npix.mask[i] && static_cast<double>(npix.data[i]) >=
                                                 static_cast<double>(threshold_contributing_pixels)) ? 1 : 0;
    }
    else
        s->channel = mask_of(channel, meta, "channel");
    lsdtt::rupture_parameters base;
    base.cohesion_pa = cohesion_pa;
    base.friction_coefficient = friction_coefficient;
    base.density_kg_m3 = density_kg_m3;
    lsdtt::validate_rupture_parameters(base);
    s->overrides = region_overrides(region, region_parameters, all_region_keys());
    for (const auto& kv : s->overrides)
    {
        lsdtt::rupture_parameters p = base;
        p.cohesion_pa = value_or(kv.second, "rupture_cohesion_pa", p.cohesion_pa);
        p.friction_coefficient = value_or(kv.second, "rupture_friction_coefficient", p.friction_coefficient);
        p.density_kg_m3 = value_or(kv.second, "rupture_density_kg_m3", p.density_kg_m3);
        try { lsdtt::validate_rupture_parameters(p); }
        catch (const std::exception& e)
        {
            throw std::invalid_argument("region " + std::to_string(kv.first) + ": " + e.what());
        }
        s->by_region[kv.first] = p;
    }
    s->region = region_ids(region, meta);
    s->water_table_depth = values_of(water_table_depth, meta, "water_table_depth");
    s->result = lsdtt::compute_rupture_planes(*s->flow, s->geometry, s->channel, s->region, s->by_region, n_threads,
                                              s->water_table_depth);
    return RupturePlanes(std::move(s));
}

// The equal-density (c, mu) samples of each region, as the CLI builds them.
struct uncertainty_args
{
    std::string cohesion_distribution, friction_distribution;
    double c_min, c_max, c_sd, mu_min, mu_max, mu_sd;
    int samples;
};

std::map<std::int32_t, std::vector<lsdtt::rupture_parameter_sample>> samples_by_region(const rupture_state& s,
                                                                                         const uncertainty_args& u)
{
    if (u.samples < 1) throw std::invalid_argument("uncertainty_samples_per_parameter must be at least 1");
    std::map<std::int32_t, std::vector<lsdtt::rupture_parameter_sample>> out;
    for (const auto& kv : s.by_region)
    {
        const auto& row = s.overrides.at(kv.first);
        const lsdtt::rupture_parameters& c0 = kv.second;
        try
        {
            out[kv.first] = lsdtt::rupture_parameter_samples(
                c0,
                lsdtt::rupture_parameter_distribution(
                    u.cohesion_distribution, c0.cohesion_pa, value_or(row, "rupture_cohesion_pa_min", u.c_min),
                    value_or(row, "rupture_cohesion_pa_max", u.c_max), value_or(row, "rupture_cohesion_pa_sd", u.c_sd),
                    "cohesion_pa"),
                u.samples,
                lsdtt::rupture_parameter_distribution(
                    u.friction_distribution, c0.friction_coefficient,
                    value_or(row, "rupture_friction_coefficient_min", u.mu_min),
                    value_or(row, "rupture_friction_coefficient_max", u.mu_max),
                    value_or(row, "rupture_friction_coefficient_sd", u.mu_sd), "friction_coefficient"),
                u.samples);
        }
        catch (const std::exception& e)
        {
            throw std::invalid_argument(std::string("failure_probability")
                                        + (s.region.empty() ? "" : " (region " + std::to_string(kv.first) + ")") + ": "
                                        + e.what());
        }
    }
    return out;
}

nb::dict probability_dict(const lsdtt::rupture_probability_rasters& p)
{
    nb::dict out;
    out["failure_probability"] = Raster(p.failure_probability);
    out["mean_safety_factor"] = Raster(p.mean_safety_factor);
    out["sd_safety_factor"] = Raster(p.sd_safety_factor);
    out["n_sources"] = p.n_sources;
    out["n_possible_failure"] = p.n_possible_failure;
    out["n_certain_failure"] = p.n_certain_failure;
    out["n_samples"] = static_cast<std::int64_t>(p.max_samples);
    return out;
}

nb::dict rupture_failure_probability(const RupturePlanes& rp, const std::string& cd, double c_min, double c_max,
                                     double c_sd, const std::string& fd, double mu_min, double mu_max, double mu_sd,
                                     int samples)
{
    const rupture_state& s = rp.state();
    const auto sm = samples_by_region(s, {cd, fd, c_min, c_max, c_sd, mu_min, mu_max, mu_sd, samples});
    return probability_dict(lsdtt::compute_rupture_failure_probability(*s.flow, s.geometry, s.channel, s.region, sm,
                                                                       s.n_threads, s.water_table_depth));
}

nb::dict landslide_catalogues(const RupturePlanes& rp, int min_daylight_neighbours)
{
    const rupture_state& s = rp.state();
    const meta_t& meta = s.flow->meta();
    const lsdtt::landslide_catalogues cat =
        lsdtt::build_landslide_catalogues(s.result, meta, s.region, s.by_region, min_daylight_neighbours);
    auto table = [](const std::vector<lsdtt::landslide_record>& rows, bool by_daylight) {
        nb::dict d;
        std::vector<std::int64_t> id, ndp, ncell, npiece;
        std::vector<double> dx, dy, area, md, xd, vol, msf, lsf, cx, cy, len, wid, lw, per;
        for (const auto& r : rows)
        {
            id.push_back(r.id);
            ndp.push_back(r.n_daylight_points);
            ncell.push_back(r.n_cells);
            npiece.push_back(r.n_pieces);
            dx.push_back(r.daylight_x);
            dy.push_back(r.daylight_y);
            area.push_back(r.area_m2);
            md.push_back(r.mean_depth_m);
            xd.push_back(r.max_depth_m);
            vol.push_back(r.volume_m3);
            msf.push_back(r.mean_safety_factor);
            lsf.push_back(r.landslide_safety_factor);
            cx.push_back(r.centroid_x);
            cy.push_back(r.centroid_y);
            len.push_back(r.length_m);
            wid.push_back(r.width_m);
            lw.push_back(r.length_to_width);
            per.push_back(r.perimeter_m);
        }
        d["landslide_id"] = id;
        if (by_daylight) { d["daylight_x"] = dx; d["daylight_y"] = dy; }
        else d["n_daylight_points"] = ndp;
        d["n_cells"] = ncell;
        d["area_m2"] = area;
        d["mean_depth_m"] = md;
        d["max_depth_m"] = xd;
        d["volume_m3"] = vol;
        d["mean_safety_factor"] = msf;
        d["landslide_safety_factor"] = lsf;
        d["centroid_x"] = cx;
        d["centroid_y"] = cy;
        d["length_m"] = len;
        d["width_m"] = wid;
        d["length_to_width"] = lw;
        d["perimeter_m"] = per;
        if (by_daylight) d["n_pieces"] = npiece;
        return data_frame(d);
    };
    nb::dict out;
    out["by_daylight"] = table(cat.by_daylight, true);
    out["amalgamated"] = table(cat.amalgamated, false);
    out["label"] = label_raster(cat.label, meta);
    out["label_amalgamated"] = label_raster(cat.label_amalgamated, meta);
    out["n_daylight_points"] = cat.n_daylight_points;
    out["n_daylight_points_removed"] = cat.n_daylight_points_removed;
    out["n_landslide_cells"] = cat.n_landslide_cells;
    return out;
}

// -------------------------------------------------------------- slip surfaces
class SlipSurfaces
{
public:
    std::shared_ptr<rupture_state> rupture;
    std::string family_name, method_name;
    lsdtt::slip_surface_options options;
    lsdtt::ellipsoid_options ellipsoid;
    bool is_ellipsoid = false;
    lsdtt::ellipsoid_surfaces ell;     // ellipsoid only
    lsdtt::slip_surface_rasters slip;
};

lsdtt::slip_surface_family family_of(const std::string& f)
{
    if (f == "plane_integrated") return lsdtt::slip_surface_family::plane_integrated;
    if (f == "circular") return lsdtt::slip_surface_family::circular;
    if (f == "bilinear") return lsdtt::slip_surface_family::bilinear;
    if (f == "structural") return lsdtt::slip_surface_family::structural;
    if (f == "ellipsoid") return lsdtt::slip_surface_family::ellipsoid;
    throw std::invalid_argument("family must be plane_integrated, circular, bilinear, structural or ellipsoid");
}

const char* method_name(lsdtt::limit_equilibrium_method m)
{
    using lem = lsdtt::limit_equilibrium_method;
    return m == lem::ordinary ? "ordinary" : m == lem::bishop ? "bishop" : m == lem::janbu ? "janbu" : "spencer";
}

SlipSurfaces slip_surfaces(const RupturePlanes& rp, const std::string& family, const std::string& method,
                           const std::string& search, const std::string& back_scarp, bool janbu_correction,
                           nb::object structural_dip_gradient, nb::object structural_dip_direction_degrees,
                           const std::string& column_method, const std::string& ellipsoid_search,
                           int ellipsoid_random_draws, std::int64_t ellipsoid_random_seed,
                           double ellipsoid_min_half_axis_pixels, nb::object layer_depths,
                           double screen_safety_factor, int toe_window_pixels, double max_length_metres,
                           int slices_per_cell, int search_points)
{
    SlipSurfaces out;
    out.rupture = rp.shared();
    const rupture_state& s = *out.rupture;
    const meta_t& meta = s.flow->meta();
    lsdtt::slip_surface_options& o = out.options;
    o.family = family_of(family);
    out.family_name = family;
    out.is_ellipsoid = o.family == lsdtt::slip_surface_family::ellipsoid;
    using lem = lsdtt::limit_equilibrium_method;
    if (method == "ordinary") o.evaluation.method = lem::ordinary;
    else if (method == "bishop") o.evaluation.method = lem::bishop;
    else if (method == "janbu") o.evaluation.method = lem::janbu;
    else if (method == "spencer") o.evaluation.method = lem::spencer;
    else if (method == "auto") o.evaluation.method = lsdtt::default_limit_equilibrium_method(o.family);
    else throw std::invalid_argument("method must be auto, ordinary, bishop, janbu or spencer");
    out.method_name = out.is_ellipsoid ? column_method : method_name(o.evaluation.method);
    o.evaluation.janbu_correction = janbu_correction;
    if (back_scarp == "vertical") o.back_scarp = lsdtt::slip_back_scarp::vertical;
    else if (back_scarp == "rankine") o.back_scarp = lsdtt::slip_back_scarp::rankine;
    else throw std::invalid_argument("back_scarp must be vertical or rankine");
    if (search != "seeded" && search != "full") throw std::invalid_argument("search must be seeded or full");
    o.full_search = search == "full";
    scalar_or_raster(structural_dip_gradient, meta, "structural_dip_gradient", o.structural_dip_gradient,
                     o.structural_dip_gradient_raster);
    scalar_or_raster(structural_dip_direction_degrees, meta, "structural_dip_direction_degrees",
                     o.structural_dip_direction_degrees, o.structural_dip_direction_raster);
    o.screen_safety_factor = screen_safety_factor;
    o.toe_window = toe_window_pixels;
    o.max_length_metres = max_length_metres;
    o.evaluation.slices_per_cell = slices_per_cell;
    o.search_points = search_points;
    lsdtt::ellipsoid_options& e = out.ellipsoid;
    e.screen_safety_factor = screen_safety_factor;
    e.search_points = search_points;
    if (column_method == "xie_revised") e.method = lsdtt::column_method::xie_revised;
    else if (column_method == "hovland") e.method = lsdtt::column_method::hovland;
    else throw std::invalid_argument("column_method must be xie_revised or hovland");
    if (ellipsoid_search != "seeded" && ellipsoid_search != "random")
        throw std::invalid_argument("ellipsoid_search must be seeded or random");
    e.random_search = ellipsoid_search == "random";
    e.random_draws = ellipsoid_random_draws;
    if (ellipsoid_random_seed < 0) throw std::invalid_argument("ellipsoid_random_seed must be >= 0");
    e.seed = static_cast<std::uint64_t>(ellipsoid_random_seed);
    e.min_half_axis_pixels = ellipsoid_min_half_axis_pixels;
    if (!layer_depths.is_none())
        for (nb::handle h : layer_depths)
            e.layer_depths.push_back(values_of(nb::borrow<nb::object>(h), meta, "layer_depths"));
    if (out.is_ellipsoid)
    {
        out.ell = lsdtt::compute_ellipsoid_surfaces(s.geometry, s.region, s.by_region, s.result, e, s.n_threads,
                                                    s.water_table_depth);
        out.slip = out.ell.slip;
    }
    else
    {
        if (!e.layer_depths.empty()) throw std::invalid_argument("layer_depths is for the ellipsoid family only");
        out.slip = lsdtt::compute_slip_surfaces(*s.flow, s.geometry, s.channel, s.region, s.by_region, s.result, o,
                                                s.n_threads, s.water_table_depth);
    }
    return out;
}

nb::dict slip_failure_probability(const SlipSurfaces& ss, const std::string& cd, double c_min, double c_max,
                                  double c_sd, const std::string& fd, double mu_min, double mu_max, double mu_sd,
                                  int samples)
{
    const rupture_state& s = *ss.rupture;
    const auto sm = samples_by_region(s, {cd, fd, c_min, c_max, c_sd, mu_min, mu_max, mu_sd, samples});
    if (ss.is_ellipsoid)
    {
        const lsdtt::ellipsoid_probability ep = lsdtt::compute_ellipsoid_failure_probability(
            ss.ell, s.geometry, s.region, s.by_region, sm, ss.ellipsoid);
        nb::dict out = probability_dict(ep.crown);
        out["cover_failure_probability"] = Raster(ep.cover_failure_probability);
        return out;
    }
    return probability_dict(lsdtt::compute_slip_failure_probability(ss.slip, *s.flow, s.geometry, s.channel, s.region,
                                                                    sm, ss.options.evaluation, s.n_threads,
                                                                    s.water_table_depth));
}

nb::dict slip_spencer_check(const SlipSurfaces& ss)
{
    if (ss.is_ellipsoid) throw std::invalid_argument("spencer_check is for the 2-D families, not ellipsoid");
    const rupture_state& s = *ss.rupture;
    std::int64_t no_sol = 0;
    const lsdtt::lsdtt_raster_cpu sp = lsdtt::compute_slip_spencer_check(
        ss.slip, *s.flow, s.geometry, s.channel, s.region, s.by_region, ss.options.evaluation, s.n_threads,
        s.water_table_depth, &no_sol);
    nb::dict out;
    out["safety_factor"] = Raster(sp);
    out["n_no_solution"] = no_sol;
    return out;
}

nb::dict slip_body_catalogues(const SlipSurfaces& ss)
{
    const meta_t& meta = ss.rupture->flow->meta();
    const lsdtt::slip_body_catalogues cat = lsdtt::build_slip_body_catalogues(ss.slip, meta);
    auto table = [](const std::vector<lsdtt::slip_body_record>& rows, bool by_toe) {
        std::vector<std::int64_t> id, nt, ncr, nc;
        std::vector<double> tx, ty, area, md, xd, vol, mnsf, mesf, cx, cy, len, wid, lw;
        for (const auto& r : rows)
        {
            id.push_back(r.id);
            nt.push_back(r.n_toes);
            ncr.push_back(r.n_crowns);
            nc.push_back(r.n_cells);
            tx.push_back(r.toe_x);
            ty.push_back(r.toe_y);
            area.push_back(r.area_m2);
            md.push_back(r.mean_depth_m);
            xd.push_back(r.max_depth_m);
            vol.push_back(r.volume_m3);
            mnsf.push_back(r.min_safety_factor);
            mesf.push_back(r.mean_safety_factor);
            cx.push_back(r.centroid_x);
            cy.push_back(r.centroid_y);
            len.push_back(r.length_m);
            wid.push_back(r.width_m);
            lw.push_back(r.length_to_width);
        }
        nb::dict d;
        d["body_id"] = id;
        if (by_toe) { d["toe_x"] = tx; d["toe_y"] = ty; }
        else d["n_toes"] = nt;
        d["n_crowns"] = ncr;
        d["n_cells"] = nc;
        d["area_m2"] = area;
        d["mean_depth_m"] = md;
        d["max_depth_m"] = xd;
        d["volume_m3"] = vol;
        d["min_safety_factor"] = mnsf;
        d["mean_safety_factor"] = mesf;
        d["centroid_x"] = cx;
        d["centroid_y"] = cy;
        d["length_m"] = len;
        d["width_m"] = wid;
        d["length_to_width"] = lw;
        return data_frame(d);
    };
    nb::dict out;
    out["by_toe"] = table(cat.by_toe, true);
    out["amalgamated"] = table(cat.amalgamated, false);
    out["label_amalgamated"] = label_raster(cat.label_amalgamated, meta);
    return out;
}

// ------------------------------------------------------------ frequency - area
nb::tuple landslide_frequency_area(const std::vector<double>& areas_m2, double cell_area_m2, int bins_per_decade,
                                   std::int64_t min_tail)
{
    if (!(cell_area_m2 > 0.0)) throw std::invalid_argument("cell_area_m2 must be > 0");
    if (bins_per_decade < 1) throw std::invalid_argument("bins_per_decade must be at least 1");
    const lsdtt::frequency_area_result r = lsdtt::compute_frequency_area(
        areas_m2, cell_area_m2, bins_per_decade, static_cast<std::size_t>(std::max<std::int64_t>(2, min_tail)));
    std::vector<std::int64_t> cnt(r.density.count.begin(), r.density.count.end());
    nb::dict d;
    d["area_lower_m2"] = r.density.lower;
    d["area_upper_m2"] = r.density.upper;
    d["area_centre_m2"] = r.density.centre;
    d["n_landslides"] = cnt;
    d["probability_density_per_m2"] = r.density.density;
    nb::dict f;
    f["n_landslides"] = static_cast<std::int64_t>(r.fit.n_total);
    f["method"] = "mle_discrete_ks_xmin";
    f["area_min_m2"] = r.fit.x_min;
    f["eta"] = r.eta;
    f["eta_standard_error"] = r.eta_standard_error;
    f["n_tail"] = static_cast<std::int64_t>(r.fit.n_tail);
    f["ks_distance"] = r.fit.ks_distance;
    f["binned_eta"] = r.binned_eta;
    f["binned_eta_standard_error"] = r.binned_eta_standard_error;
    f["n_bins_fitted"] = static_cast<std::int64_t>(r.binned_n_bins);
    return nb::make_tuple(data_frame(d), f);
}

// ---------------------------------------------------------------- groundwater
class GroundwaterDupuit
{
public:
    GroundwaterDupuit(const Raster& surface, const Raster& fixed_head, double regolith_thickness_metres,
                      nb::object regolith_thickness, double hydraulic_conductivity_m_per_s, double drainable_porosity,
                      const std::string& transmissivity_profile, double conductivity_decay_depth_metres,
                      int n_threads)
        : meta_(surface.native().get_meta())
    {
        require_grid(fixed_head.native().get_meta(), meta_, "fixed_head");
        const std::size_t n = surface.native().get_data().size();
        std::vector<double> z(n), zb(n);
        valid_ = surface.native().get_mask();
        std::vector<double> thick(n, regolith_thickness_metres);
        if (!regolith_thickness.is_none()) thick = values_of(regolith_thickness, meta_, "regolith_thickness");
        else if (!(regolith_thickness_metres > 0.0)) throw std::invalid_argument("regolith_thickness_metres must be > 0");
        for (std::size_t i = 0; i < n; ++i)
        {
            if (!(thick[i] > 0.0) || !std::isfinite(thick[i])) valid_[i] = 0;
            z[i] = surface.native().get_data()[i];
            zb[i] = z[i] - (std::isfinite(thick[i]) ? thick[i] : 0.0);
        }
        lsdtt::groundwater_dupuit_parameters p;
        p.hydraulic_conductivity_m_per_s = hydraulic_conductivity_m_per_s;
        p.drainable_porosity = drainable_porosity;
        if (transmissivity_profile == "uniform") p.profile = lsdtt::groundwater_transmissivity_profile::uniform;
        else if (transmissivity_profile == "exponential") p.profile = lsdtt::groundwater_transmissivity_profile::exponential;
        else throw std::invalid_argument("transmissivity_profile must be uniform or exponential");
        p.conductivity_decay_depth_metres = conductivity_decay_depth_metres;
        const auto fixed = mask_of(nb::cast(fixed_head), meta_, "fixed_head");
        g_ = std::make_unique<lsdtt::lsdtt_groundwater_dupuit>(meta_, z, zb, valid_, fixed, std::vector<double>{}, p,
                                                               n_threads);
    }

    std::vector<double> recharge(nb::object f) const
    {
        const std::size_t n = valid_.size();
        if (nb::isinstance<Raster>(f))
        {
            std::vector<double> v = values_of(f, meta_, "recharge_m_per_s");
            for (double& x : v)
                if (!std::isfinite(x)) x = 0.0;
            return v;
        }
        return std::vector<double>(n, nb::cast<double>(f));
    }

    nb::dict run_to_steady(nb::object f, const std::string& method, double tolerance, std::int64_t max_substeps,
                           int max_newton_iterations)
    {
        const std::vector<double> rf = recharge(f);
        nb::dict out;
        out["steady_method"] = method;
        if (method == "implicit")
        {
            if (max_newton_iterations < 1) throw std::invalid_argument("max_newton_iterations must be at least 1");
            lsdtt::lsdtt_groundwater_dupuit::implicit_options io;
            io.relative_tolerance = tolerance;
            io.max_steps = max_newton_iterations;
            const auto rep = g_->run_to_steady_implicit(rf, io);
            if (!rep.converged)
                throw std::runtime_error("the implicit groundwater steady state did not converge: " + rep.message
                                         + " (raise max_newton_iterations or steady_tolerance, or try "
                                           "steady_method='explicit')");
            out["steps"] = rep.substeps;
            out["linear_iterations"] = rep.linear_iterations;
            out["used_direct_solver"] = rep.used_direct_solver;
            out["n_saturated"] = rep.n_saturated;
            out["n_dry"] = rep.n_dry;
            out["residual"] = rep.residual;
            out["pseudo_time_s"] = rep.time;
        }
        else if (method == "explicit")
        {
            if (max_substeps < 1) throw std::invalid_argument("max_substeps must be at least 1");
            const auto rep = g_->run_to_steady(rf, tolerance, max_substeps);
            if (!rep.converged)
                throw std::runtime_error("the groundwater steady state did not converge within max_substeps (raise "
                                         "it, or steady_tolerance)");
            out["steps"] = rep.substeps;
            out["residual"] = rep.residual;
            out["pseudo_time_s"] = rep.time;
        }
        else throw std::invalid_argument("steady_method must be implicit or explicit");
        return out;
    }

    nb::dict run(double duration_s, nb::object f, const std::string& method, double max_dt_s, double dt_tolerance)
    {
        if (!(duration_s > 0.0)) throw std::invalid_argument("duration_s must be > 0");
        const std::vector<double> rf = recharge(f);
        nb::dict out;
        out["transient_method"] = method;
        if (method == "explicit")
            out["steps"] = g_->run(duration_s, rf);
        else if (method == "implicit")
        {
            lsdtt::lsdtt_groundwater_dupuit::implicit_options o;
            o.max_dt_s = max_dt_s;
            o.dt_tolerance_metres = dt_tolerance;
            const auto r = g_->run_implicit(duration_s, rf, o);
            if (!r.ok) throw std::runtime_error("the implicit groundwater transient failed: " + r.message);
            out["steps"] = r.steps;
            out["newton_iterations"] = r.newton_iterations;
            out["linear_iterations"] = r.linear_iterations;
            out["splits"] = r.splits;
            out["rejected_by_error"] = r.rejected_by_error;
            out["max_error_estimate_metres"] = r.max_error_estimate_metres;
        }
        else throw std::invalid_argument("transient_method must be explicit or implicit");
        return out;
    }

    void set_thickness(const Raster& h)
    {
        std::vector<double> v = values_of(nb::cast(h), meta_, "thickness");
        for (double& x : v)
            if (!std::isfinite(x)) x = 0.0;
        g_->set_thickness(v);
    }

    Raster masked(std::vector<double> v) const
    {
        const auto wt = g_->water_table_depth();
        for (std::size_t i = 0; i < v.size(); ++i)
            if (!std::isfinite(wt[i])) v[i] = std::numeric_limits<double>::quiet_NaN();
        return raster_of(v, meta_);
    }
    Raster water_table_depth() const { return raster_of(g_->water_table_depth(), meta_); }
    Raster saturated_thickness() const { return masked(g_->thickness()); }
    Raster seepage() const { return masked(g_->seepage()); }
    double storage() const { return g_->storage(); }
    nb::dict budget() const
    {
        const auto& b = g_->budget();
        nb::dict d;
        d["recharge_m3"] = b.recharge_in;
        d["storage_change_m3"] = b.storage_change;
        d["outflow_fixed_head_m3"] = b.outflow_fixed_head;
        d["seepage_m3"] = b.seepage_out;
        d["clipped_added_m3"] = b.clipped_added;
        return d;
    }
    nb::dict steady_rates() const
    {
        const auto& r = g_->steady_rates();
        nb::dict d;
        d["recharge_m3_per_s"] = r.recharge;
        d["outflow_fixed_head_m3_per_s"] = r.outflow_fixed_head;
        d["seepage_m3_per_s"] = r.seepage;
        d["imbalance_m3_per_s"] = r.imbalance;
        return d;
    }

private:
    meta_t meta_;
    std::vector<std::uint8_t> valid_;
    std::unique_ptr<lsdtt::lsdtt_groundwater_dupuit> g_;
};

// Docstring fragments shared by the probability methods.
#define LSDTT_UNCERTAINTY_ARGS                                                                                    \
    nb::kw_only(), nb::arg("cohesion_distribution") = "uniform", nb::arg("cohesion_pa_min") = 500.0,            \
        nb::arg("cohesion_pa_max") = 1500.0, nb::arg("cohesion_pa_sd") = 250.0,                                  \
        nb::arg("friction_coefficient_distribution") = "uniform", nb::arg("friction_coefficient_min") = 0.6,     \
        nb::arg("friction_coefficient_max") = 0.8, nb::arg("friction_coefficient_sd") = 0.05,                    \
        nb::arg("uncertainty_samples_per_parameter") = 5
#define LSDTT_UNCERTAINTY_SIG                                                                                     \
    "*, cohesion_distribution: str = 'uniform', cohesion_pa_min: float = 500.0, cohesion_pa_max: float = 1500.0, " \
    "cohesion_pa_sd: float = 250.0, friction_coefficient_distribution: str = 'uniform', "                         \
    "friction_coefficient_min: float = 0.6, friction_coefficient_max: float = 0.8, "                             \
    "friction_coefficient_sd: float = 0.05, uncertainty_samples_per_parameter: int = 5"
#define LSDTT_UNCERTAINTY_DOC                                                                                     \
    "cohesion_distribution : {\"uniform\", \"normal\"}, default=\"uniform\", keyword-only\n"                      \
    "    uniform on [cohesion_pa_min, cohesion_pa_max], or normal about the central cohesion with\n"             \
    "    cohesion_pa_sd, truncated to [max(0, c - 3 sd), c + 3 sd].\n"                                             \
    "cohesion_pa_min, cohesion_pa_max, cohesion_pa_sd : float, default=500.0, 1500.0, 250.0, keyword-only\n"     \
    "    Bounds (uniform) or standard deviation (normal) of cohesion (Pa).\n"                                     \
    "friction_coefficient_distribution : {\"uniform\", \"normal\"}, default=\"uniform\", keyword-only\n"         \
    "    As cohesion_distribution, for the friction coefficient.\n"                                               \
    "friction_coefficient_min, friction_coefficient_max, friction_coefficient_sd : float, default=0.6, 0.8, "    \
    "0.05, keyword-only\n"                                                                                         \
    "    Bounds or standard deviation of the friction coefficient.\n"                                             \
    "uncertainty_samples_per_parameter : int, default=5, keyword-only\n"                                          \
    "    Equal-density samples of each parameter; their tensor product is used.\n"                                \
    "\n"                                                                                                           \
    "A region's rupture_cohesion_pa_min (and the other bounds and SDs) in the\n"                                  \
    "region_parameters given to ``FlowInfo.rupture_planes`` override these.\n"

#define LSDTT_PROBABILITY_RETURNS                                                                                 \
    "dict\n"                                                                                                       \
    "    ``failure_probability`` (Raster, P(F <= 1)), ``mean_safety_factor``,\n"                                  \
    "    ``sd_safety_factor`` (Rasters, weighted over the samples), ``n_sources``,\n"                              \
    "    ``n_possible_failure`` (Pf > 0), ``n_certain_failure`` (Pf = 1), ``n_samples``.\n"

// ------------------------------------------------- DEM of difference
Raster float_raster(const std::vector<float>& values, const meta_t& meta)
{
    std::vector<double> v(values.size(), std::numeric_limits<double>::quiet_NaN());
    const float nodata = static_cast<float>(meta.no_data_value);
    for (std::size_t i = 0; i < values.size(); ++i)
        if (values[i] != nodata) v[i] = values[i];
    return raster_of(v, meta);
}

nb::dict coseismic_hazard_area(const Raster& routing_surface, double initiation_slope, double stopping_slope,
                               const std::string& boundary_conditions, int n_threads)
{
    lsdtt::coseismic_hazard_area_options o;
    o.initiation_slope = initiation_slope;
    o.stopping_slope = stopping_slope;
    o.boundary_conditions = boundary_conditions;
    o.n_threads = n_threads;
    const lsdtt::coseismic_hazard_area_result r = lsdtt::compute_coseismic_hazard_area(routing_surface.native(), o);
    nb::dict out;
    out["hazard_area"] = float_raster(r.hazard_area.get_data(), r.hazard_area.get_meta());
    out["steepest_slope"] = float_raster(r.steepest_slope.get_data(), r.steepest_slope.get_meta());
    out["zone"] = float_raster(r.zone.get_data(), r.zone.get_meta());
    out["edge_influenced"] = float_raster(r.edge_influenced.get_data(), r.edge_influenced.get_meta());
    out["n_valid_cells"] = r.n_valid;
    out["n_source_cells"] = r.n_source;
    out["n_stopping_cells"] = r.n_stopping;
    out["n_zero_hazard_cells"] = r.n_zero_hazard;
    out["n_edge_influenced_cells"] = r.n_edge_influenced;
    out["source_area_m2"] = r.source_area_m2;
    return out;
}

nb::dict dem_change_landslides(const Raster& difference, const Raster& dem, nb::object stable_mask,
                               bool is_newer_minus_older, const std::string& change_sign, double detection_sigma,
                               double seed_sigma, double patch_significance, int min_patch_pixels,
                               int aggregation_radius_pixels, double min_slope, int n_slope_bins,
                               int n_null_simulations, std::int64_t random_seed, int n_threads)
{
    const meta_t& meta = dem.native().get_meta();
    require_grid(difference.native().get_meta(), meta, "difference");
    std::vector<float> dv = difference.native().get_data();
    std::vector<std::uint8_t> dm = difference.native().get_mask();
    if (!is_newer_minus_older)
        for (std::size_t i = 0; i < dv.size(); ++i)
            if (dm[i]) dv[i] = -dv[i];
    meta_t difference_meta = meta;
    difference_meta.no_data_value = difference.native().get_meta().no_data_value;
    const lsdtt::lsdtt_raster_cpu signed_difference(std::move(dv), std::move(dm), difference_meta);
    const std::vector<std::uint8_t> stable = mask_of(stable_mask, meta, "stable_mask");

    lsdtt::dem_change_options o;
    o.change_sign = change_sign;
    o.detection_sigma = detection_sigma;
    o.seed_sigma = seed_sigma;
    o.patch_significance = patch_significance;
    o.min_patch_cells = min_patch_pixels;
    o.aggregation_radius_pixels = aggregation_radius_pixels;
    o.min_slope = min_slope;
    o.n_slope_bins = n_slope_bins;
    o.n_null_simulations = n_null_simulations;
    o.random_seed = static_cast<std::uint64_t>(random_seed);
    o.n_threads = n_threads;
    const lsdtt::dem_change_result r = lsdtt::detect_dem_change_landslides(
        signed_difference, lsdtt::dem_change_slope(dem.native()), stable_mask.is_none() ? nullptr : &stable, o);

    const lsdtt::dem_change_uncertainty& u = r.uncertainty;
    nb::dict uncertainty;
    uncertainty["used_stable_mask"] = u.used_stable_mask;
    uncertainty["n_stable_cells"] = u.n_stable;
    uncertainty["n_clipped_cells"] = u.n_clipped;
    uncertainty["bias_metres"] = u.bias_m;
    uncertainty["nmad_metres"] = u.nmad_m;
    uncertainty["standard_deviation_metres"] = u.standard_deviation_m;
    uncertainty["rmse_metres"] = u.rmse_m;
    uncertainty["p2.5_metres"] = u.p025_m;
    uncertainty["p97.5_metres"] = u.p975_m;
    uncertainty["variogram_fitted"] = u.variogram_fitted;
    uncertainty["nugget_fraction"] = u.nugget_fraction;
    uncertainty["correlation_range_metres"] = u.correlation_range_m;
    uncertainty["sill"] = u.sill;

    nb::dict by_slope;
    by_slope["slope_lower"] = u.slope_lower;
    by_slope["slope_upper"] = u.slope_upper;
    by_slope["slope_median"] = u.slope_median;
    by_slope["n_cells"] = u.slope_n;
    by_slope["bias_metres"] = u.slope_bias_m;
    by_slope["nmad_metres"] = u.slope_nmad_m;

    nb::dict variogram;
    variogram["lag_metres"] = u.lag_m;
    variogram["semivariance"] = u.semivariance;
    variogram["n_pairs"] = u.lag_n_pairs;

    std::vector<std::int64_t> id, n_cells;
    std::vector<std::string> change;
    std::vector<double> area, area_u, perimeter, cx, cy, length, width, mean, extreme, mean_s, sig, vol, vol_s, slope;
    for (const auto& p : r.patches)
    {
        id.push_back(p.id); change.push_back(p.sign < 0 ? "loss" : "gain"); n_cells.push_back(p.n_cells);
        area.push_back(p.area_m2); area_u.push_back(p.area_uncertainty_m2); perimeter.push_back(p.perimeter_m);
        cx.push_back(p.centroid_x); cy.push_back(p.centroid_y); length.push_back(p.length_m); width.push_back(p.width_m);
        mean.push_back(p.mean_change_m); extreme.push_back(p.extreme_change_m); mean_s.push_back(p.mean_change_sigma_m);
        sig.push_back(p.significance); vol.push_back(p.volume_m3); vol_s.push_back(p.volume_sigma_m3); slope.push_back(p.mean_slope);
    }
    nb::dict patches;
    patches["patch_id"] = id; patches["change"] = change; patches["n_cells"] = n_cells; patches["area_m2"] = area;
    patches["area_uncertainty_m2"] = area_u; patches["perimeter_m"] = perimeter; patches["centroid_x"] = cx;
    patches["centroid_y"] = cy; patches["length_m"] = length; patches["width_m"] = width;
    patches["mean_change_m"] = mean; patches["extreme_change_m"] = extreme; patches["mean_change_sigma_m"] = mean_s;
    patches["significance"] = sig; patches["volume_m3"] = vol; patches["volume_sigma_m3"] = vol_s;
    patches["mean_slope"] = slope;

    nb::dict coverage;
    std::int64_t n_loss = 0, n_gain = 0;
    for (const auto& p : r.patches) (p.sign < 0 ? n_loss : n_gain)++;
    const double assessed = r.valid_area_m2 > 0.0 ? r.valid_area_m2 : 1.0;
    coverage["change"] = std::vector<std::string>{"loss", "gain"};
    coverage["n_patches"] = std::vector<std::int64_t>{n_loss, n_gain};
    coverage["area_m2"] = std::vector<double>{r.loss_area_m2, r.gain_area_m2};
    coverage["area_uncertainty_m2"] = std::vector<double>{r.loss_area_uncertainty_m2, r.gain_area_uncertainty_m2};
    coverage["fraction_of_assessed_area"] = std::vector<double>{r.loss_area_m2 / assessed, r.gain_area_m2 / assessed};
    coverage["volume_m3"] = std::vector<double>{r.loss_volume_m3, r.gain_volume_m3};
    coverage["volume_sigma_m3"] = std::vector<double>{r.loss_volume_sigma_m3, r.gain_volume_sigma_m3};
    coverage["assessed_area_m2"] = std::vector<double>{r.valid_area_m2, r.valid_area_m2};

    std::vector<double> sk, ss, sa, sf, sv;
    std::vector<std::int64_t> sn;
    for (const auto& row : r.sensitivity)
    {
        sk.push_back(row.detection_sigma); ss.push_back(row.seed_sigma); sn.push_back(row.n_patches);
        sa.push_back(row.area_m2); sf.push_back(row.coverage_fraction); sv.push_back(row.volume_m3);
    }
    nb::dict sensitivity;
    sensitivity["detection_sigma"] = sk; sensitivity["seed_sigma"] = ss; sensitivity["n_patches"] = sn;
    sensitivity["area_m2"] = sa; sensitivity["fraction_of_assessed_area"] = sf; sensitivity["volume_m3"] = sv;

    const lsdtt::dem_change_null_summary& q = r.null_simulation;
    nb::dict false_detection;
    false_detection["n_simulations"] = q.n_simulations;
    false_detection["mean_false_patches"] = q.mean_n_patches;
    false_detection["p05_false_patches"] = q.p05_n_patches;
    false_detection["p95_false_patches"] = q.p95_n_patches;
    false_detection["mean_false_area_m2"] = q.mean_area_m2;
    false_detection["p05_false_area_m2"] = q.p05_area_m2;
    false_detection["p95_false_area_m2"] = q.p95_area_m2;
    false_detection["area_false_discovery_ratio"] = q.area_false_discovery_ratio;
    std::vector<double> bl, bu, bs, bf;
    std::vector<std::int64_t> bo;
    for (const auto& b : q.by_area)
    {
        bl.push_back(b.area_lower_m2); bu.push_back(b.area_upper_m2); bo.push_back(b.n_observed);
        bs.push_back(b.mean_n_simulated); bf.push_back(b.false_discovery_ratio);
    }
    nb::dict by_area;
    by_area["area_lower_m2"] = bl; by_area["area_upper_m2"] = bu; by_area["n_observed"] = bo;
    by_area["mean_n_simulated"] = bs; by_area["false_discovery_ratio"] = bf;

    std::vector<double> label(r.label.size(), std::numeric_limits<double>::quiet_NaN());
    std::vector<double> change_class(r.label.size(), std::numeric_limits<double>::quiet_NaN());
    const float nodata = static_cast<float>(difference_meta.no_data_value);
    for (std::size_t i = 0; i < label.size(); ++i)
    {
        if (r.corrected_m[i] == nodata) continue;
        change_class[i] = r.change_class[i];
        if (r.label[i] > 0) label[i] = r.label[i];
    }

    nb::dict out;
    out["uncertainty"] = uncertainty;
    out["uncertainty_by_slope"] = data_frame(by_slope);
    out["variogram"] = data_frame(variogram);
    out["patches"] = data_frame(patches);
    out["coverage"] = data_frame(coverage);
    out["sensitivity"] = data_frame(sensitivity);
    out["false_detection"] = false_detection;
    out["false_detection_by_area"] = data_frame(by_area);
    out["corrected"] = float_raster(r.corrected_m, difference_meta);
    out["sigma"] = float_raster(r.sigma_m, difference_meta);
    out["z"] = float_raster(r.z, difference_meta);
    out["change_class"] = raster_of(change_class, meta);
    out["label"] = raster_of(label, meta);
    return out;
}

} // namespace

void bind_slope_stability_graph(nb::class_<FlowGraph>& graph)
{
    graph.def("dinf_contour_width", &dinf_contour_width,
        nb::sig("def dinf_contour_width(self) -> Raster"), R"doc(
Contour width of each cell for its D-infinity flow direction.

lsdtt-slope-stability ``contour_width_method: dinf_flow_direction``:
w = b max(|sin beta|, |cos beta|) for flow bearing beta and cell size b, the
width of a D-infinity flow tube on a plane, so A / w is the flow distance
from the divide there. Cells without a direction get b. Pass the result as
``contour_width`` to ``sinmap``.

Parameters
----------
self : FlowGraph
    A FlowGraph built with ``method="dinf"``.

Returns
-------
Raster
    Contour width (m); nodata outside the DEM.
)doc");
    graph.def("dinf_edge_contaminated_mask", &dinf_edge_contaminated_mask,
        nb::sig("def dinf_edge_contaminated_mask(self) -> Raster"), R"doc(
Cells whose D-infinity contributing area reaches the raster edge or nodata.

lsdtt-slope-stability ``exclude_edge_contaminated_cells``: the area of these
cells may be truncated, so their specific catchment area and SINMAP
stability index are lower bounds. Pass the result as ``exclude`` to
``sinmap`` to leave them out.

Parameters
----------
self : FlowGraph
    A FlowGraph built with ``method="dinf"``.

Returns
-------
Raster
    1 contaminated, 0 not; nodata outside the DEM.
)doc");
}

void bind_slope_stability(nb::module_& module, nb::class_<FlowInfo>& flow)
{
    module.def("sinmap", &sinmap, nb::arg("slope"), nb::arg("drainage_area"), nb::kw_only(),
        nb::arg("dimensionless_cohesion_min") = 0.0, nb::arg("dimensionless_cohesion_max") = 0.25,
        nb::arg("friction_coefficient_min") = 0.57735026918962573, nb::arg("friction_coefficient_max") = 1.0,
        nb::arg("transmissivity_over_recharge_min_metres") = 2000.0,
        nb::arg("transmissivity_over_recharge_max_metres") = 3000.0, nb::arg("water_to_soil_density_ratio") = 0.5,
        nb::arg("region") = nb::none(), nb::arg("region_parameters") = nb::none(), nb::arg("exclude") = nb::none(),
        nb::arg("contour_width") = nb::none(), nb::arg("n_threads") = 0,
        nb::sig("def sinmap(slope: Raster, drainage_area: Raster, *, dimensionless_cohesion_min: float = 0.0, "
                "dimensionless_cohesion_max: float = 0.25, friction_coefficient_min: float = 0.5773502691896257, "
                "friction_coefficient_max: float = 1.0, transmissivity_over_recharge_min_metres: float = 2000.0, "
                "transmissivity_over_recharge_max_metres: float = 3000.0, water_to_soil_density_ratio: float = 0.5, "
                "region: Raster | None = None, region_parameters: dict | None = None, exclude: Raster | None = None, "
                "contour_width: Raster | None = None, n_threads: int = 0) -> dict"), R"doc(
SINMAP stability index, class, probability and factor-of-safety bounds.

lsdtt-slope-stability ``run_sinmap`` (Pack et al. 1998), the same routine;
the keyword names are the CLI's less the ``sinmap_`` prefix. The infinite
slope with steady wetness w = min(1, R a / (T sin theta)) is evaluated over
uniform ranges of dimensionless cohesion C, friction coefficient tan phi and
T/R. FS_min uses (C_min, t_min, (T/R)_min), FS_max the other bounds; the
stability index is FS_min where it is above 1, P(FS > 1) where FS_min < 1 <
FS_max (closed form), and 0 where FS_max < 1. Flat cells are class 1 and
nodata in the continuous outputs. The CLI uses the D-infinity slope and area
of the carved-then-filled DEM: ``FlowGraph(dem, method="dinf")``,
``.slope()`` and ``.drainage_area()``.

Parameters
----------
slope : Raster
    Gradient S (m/m), e.g. ``FlowGraph.slope()``.
drainage_area : Raster
    Contributing area (m^2) on the same grid, e.g. ``FlowGraph.drainage_area()``.
dimensionless_cohesion_min, dimensionless_cohesion_max : float, default=0.0, 0.25, keyword-only
    C = (C_r + C_s) / (h rho_s g), bounds.
friction_coefficient_min, friction_coefficient_max : float, default=tan 30, 1.0, keyword-only
    tan phi, bounds.
transmissivity_over_recharge_min_metres, transmissivity_over_recharge_max_metres : float, default=2000.0, 3000.0, keyword-only
    T / R (m), bounds.
water_to_soil_density_ratio : float, default=0.5, keyword-only
    rho_w / rho_s.
region : Raster or None, default=None, keyword-only
    Integer region ids selecting parameter sets; nodata cells are left out.
region_parameters : dict or None, default=None, keyword-only
    ``{region_id: {name: value}}`` with the CLI's region-table column names
    (``sinmap_dimensionless_cohesion_min``, ...); a missing name takes the
    keyword value. ``rupture_*`` names are accepted and ignored, so one table
    can drive both ``sinmap`` and ``FlowInfo.rupture_planes``. Needed with
    ``region``.
exclude : Raster or None, default=None, keyword-only
    Nonzero cells are left out, e.g. ``FlowGraph.dinf_edge_contaminated_mask()``.
contour_width : Raster or None, default=None, keyword-only
    Contour width (m) for a = A / w, e.g. ``FlowGraph.dinf_contour_width()``
    (the CLI's default); None uses the cell size, as SINMAP.
n_threads : int, default=0, keyword-only
    OpenMP threads (0 = default); the result does not depend on it.

Returns
-------
dict
    Rasters ``stability_index``, ``stability_class`` (1 stable to 6 defended),
    ``probability_stable``, ``fs_min``, ``fs_max``, ``saturation`` (SINMAP's
    coding: 3 always saturated, 2 partly, else the wettest w),
    ``specific_catchment_area`` (m); counts ``n_evaluated``, ``n_flat``,
    ``n_excluded``; ``class_areas``, a DataFrame of the CLI's class-area
    table (``region_id, stability_class, n_cells, area_m2,
    fraction_of_region``).
)doc");

    module.def("sinmap_slope_area_curves", &sinmap_slope_area_curves, nb::arg("slope_max"), nb::arg("area_min"),
        nb::arg("area_max"), nb::kw_only(), nb::arg("dimensionless_cohesion_min") = 0.0,
        nb::arg("dimensionless_cohesion_max") = 0.25, nb::arg("friction_coefficient_min") = 0.57735026918962573,
        nb::arg("friction_coefficient_max") = 1.0, nb::arg("transmissivity_over_recharge_min_metres") = 2000.0,
        nb::arg("transmissivity_over_recharge_max_metres") = 3000.0, nb::arg("water_to_soil_density_ratio") = 0.5,
        nb::arg("n_slopes") = 400,
        nb::sig("def sinmap_slope_area_curves(slope_max: float, area_min: float, area_max: float, *, "
                "dimensionless_cohesion_min: float = 0.0, dimensionless_cohesion_max: float = 0.25, "
                "friction_coefficient_min: float = 0.5773502691896257, friction_coefficient_max: float = 1.0, "
                "transmissivity_over_recharge_min_metres: float = 2000.0, "
                "transmissivity_over_recharge_max_metres: float = 3000.0, water_to_soil_density_ratio: float = 0.5, "
                "n_slopes: int = 400) -> pandas.DataFrame"), R"doc(
The class boundaries of SINMAP's slope - specific catchment area plot.

lsdtt-slope-stability ``write_sinmap_slope_area_curves_table`` for one
parameter set: the saturation lines (``always_saturated_above``,
``never_saturated_below``) and the FS curves ``fs_min_1.5``, ``fs_min_1.0``
and ``fs_max_1.0`` as points. The CLI uses the largest evaluated slope,
the cell size and the largest specific catchment area as the extent.

Parameters
----------
slope_max : float
    Largest gradient (m/m) of the plot.
area_min, area_max : float
    Range of specific catchment area (m); points outside it are dropped.
dimensionless_cohesion_min, dimensionless_cohesion_max, friction_coefficient_min, friction_coefficient_max, transmissivity_over_recharge_min_metres, transmissivity_over_recharge_max_metres, water_to_soil_density_ratio
    As ``sinmap``, keyword-only.
n_slopes : int, default=400, keyword-only
    Slopes sampled along each curve.

Returns
-------
pandas.DataFrame
    ``curve, slope, specific_catchment_area_m``, the CLI table's columns
    less ``region_id``.
)doc");

    nb::class_<RupturePlanes>(module, "RupturePlanes", R"doc(
Planar rupture surfaces (scaLr) on a DEM, from ``FlowInfo.rupture_planes``.

Holds the inputs of the run, so its probability of failure, landslide
catalogues and slip surfaces are computed from the same flow, geometry,
channels, regions and water table.

Members: ``safety_factor``, ``depth``, ``rupture_slope``,
``topographic_slope``, ``length``, ``unstable``, ``channel``, ``daylight``,
``n_sources``, ``n_unstable``, ``n_channel``, ``failure_probability``,
``landslide_catalogues``, ``slip_surfaces``.
)doc")
        .def_prop_ro("safety_factor", [](const RupturePlanes& r) { return Raster(r.state().result.safety_factor); },
            "Factor of safety F at the critical daylight point (the minimum over the path); nodata on channels and cells with no admissible plane.")
        .def_prop_ro("depth", [](const RupturePlanes& r) { return Raster(r.state().result.depth); },
            "Optimal rupture depth h* below the source cell (m).")
        .def_prop_ro("rupture_slope", [](const RupturePlanes& r) { return Raster(r.state().result.rupture_slope); },
            "Gradient of the rupture plane, (dz - h*) / dx (m/m).")
        .def_prop_ro("topographic_slope", [](const RupturePlanes& r) { return Raster(r.state().result.topographic_slope); },
            "Gradient from the source to its daylight point, dz / dx (m/m).")
        .def_prop_ro("length", [](const RupturePlanes& r) { return Raster(r.state().result.length); },
            "Horizontal distance dx from the source to its daylight point (m).")
        .def_prop_ro("unstable", [](const RupturePlanes& r) { return Raster(r.state().result.unstable); },
            "1 where F <= 1, 0 where F > 1.")
        .def_prop_ro("channel", [](const RupturePlanes& r) {
                std::vector<double> v(r.state().channel.begin(), r.state().channel.end());
                return raster_of(v, r.state().flow->meta());
            }, "The channel mask used (1 channel, 0 hillslope).")
        .def_prop_ro("daylight", [](const RupturePlanes& r) { return numpy_int64(r.state().result.daylight); },
            "Flat index of each cell's daylight point (numpy int64, -1 where none).")
        .def_prop_ro("n_sources", [](const RupturePlanes& r) { return r.state().result.n_sources; },
            "Hillslope cells with an admissible daylight point.")
        .def_prop_ro("n_unstable", [](const RupturePlanes& r) { return r.state().result.n_unstable; },
            "Of the sources, those with F <= 1.")
        .def_prop_ro("n_channel", [](const RupturePlanes& r) { return r.state().result.n_channel; },
            "Channel cells (not sources).")
        .def("failure_probability", &rupture_failure_probability, LSDTT_UNCERTAINTY_ARGS,
            nb::sig("def failure_probability(self, " LSDTT_UNCERTAINTY_SIG ") -> dict"),
            "Probability of failure of the rupture planes by equal-density sampling of cohesion and friction.\n\n"
            "lsdtt-slope-stability ``parameter_uncertainty: equal_density``, the same routine. For each cell and\n"
            "each (c, mu) sample the search is the full deterministic one, so the daylight point may differ\n"
            "between samples; Pf = sum of the weights of the samples with F_min <= 1. Exact for this model\n"
            "(no Monte Carlo error). Density and the water table are held fixed.\n\n"
            "Parameters\n----------\n" LSDTT_UNCERTAINTY_DOC "\nReturns\n-------\n" LSDTT_PROBABILITY_RETURNS)
        .def("landslide_catalogues", &landslide_catalogues, nb::kw_only(), nb::arg("min_daylight_neighbours") = 2,
            nb::sig("def landslide_catalogues(self, *, min_daylight_neighbours: int = 2) -> dict"), R"doc(
The potential-landslide catalogues: by daylight point and amalgamated.

lsdtt-slope-stability ``write_landslide_catalogue`` (scaLr's MATLAB): the
unstable cells sharing a daylight point form one landslide (after
re-attaching each daylight point to its source), the 1-D rule drops daylight
points with fewer than ``min_daylight_neighbours`` daylight-point neighbours,
and 8-connected landslides are merged into the amalgamated catalogue. Areas,
depths, volumes, F_L (eq. 7) and the planform ellipse are the CLI's columns.

Parameters
----------
min_daylight_neighbours : int, default=2, keyword-only
    The 1-D rule (CLI ``rupture_min_daylight_neighbours``); 0 keeps every
    daylight point (the paper's catalogue), 0 to 8.

Returns
-------
dict
    ``by_daylight`` and ``amalgamated`` (DataFrames, the CLI's
    ``_landslide_catalogue*.csv`` columns), ``label`` and
    ``label_amalgamated`` (Rasters of the ids), ``n_daylight_points``,
    ``n_daylight_points_removed``, ``n_landslide_cells``.
)doc")
        .def("slip_surfaces", &slip_surfaces, nb::kw_only(), nb::arg("family") = "circular",
            nb::arg("method") = "auto", nb::arg("search") = "seeded", nb::arg("back_scarp") = "vertical",
            nb::arg("janbu_correction") = true, nb::arg("structural_dip_gradient") = 0.3,
            nb::arg("structural_dip_direction_degrees") = 0.0, nb::arg("column_method") = "xie_revised",
            nb::arg("ellipsoid_search") = "seeded", nb::arg("ellipsoid_random_draws") = 200,
            nb::arg("ellipsoid_random_seed") = 1, nb::arg("ellipsoid_min_half_axis_pixels") = 3.0,
            nb::arg("layer_depths") = nb::none(), nb::arg("screen_safety_factor") = 1.5,
            nb::arg("toe_window_pixels") = 0, nb::arg("max_length_metres") = 0.0, nb::arg("slices_per_cell") = 4,
            nb::arg("search_points") = 12,
            nb::sig("def slip_surfaces(self, *, family: str = 'circular', method: str = 'auto', "
                    "search: str = 'seeded', back_scarp: str = 'vertical', janbu_correction: bool = True, "
                    "structural_dip_gradient: float | Raster = 0.3, "
                    "structural_dip_direction_degrees: float | Raster = 0.0, column_method: str = 'xie_revised', "
                    "ellipsoid_search: str = 'seeded', ellipsoid_random_draws: int = 200, "
                    "ellipsoid_random_seed: int = 1, ellipsoid_min_half_axis_pixels: float = 3.0, "
                    "layer_depths: list[Raster] | None = None, screen_safety_factor: float = 1.5, "
                    "toe_window_pixels: int = 0, max_length_metres: float = 0.0, slices_per_cell: int = 4, "
                    "search_points: int = 12) -> SlipSurfaces"), R"doc(
Critical slip surfaces seeded from the rupture planes.

lsdtt-slope-stability ``slip_surface_family``, the same routines; the keyword
names are the CLI's less the ``slip_surface_`` prefix. Crowns whose rupture
F is above ``screen_safety_factor`` are screened out. The 2-D families
(method of slices on the crown's D8 section) have no side resistance;
``ellipsoid`` (method of columns) does.

Parameters
----------
family : {"plane_integrated", "circular", "bilinear", "structural", "ellipsoid"}, default="circular", keyword-only
    plane_integrated: scaLr's plane with the forces summed over the body;
    circular: arcs by sagitta; bilinear: a Rankine back segment onto a
    flatter basal plane; structural: a plane at the apparent dip of a weak
    layer, ending where it daylights; ellipsoid: 3-D ellipsoids through the
    crown and the daylight point.
method : {"auto", "ordinary", "bishop", "janbu", "spencer"}, default="auto", keyword-only
    Limit-equilibrium method of the 2-D families; auto is ordinary for
    plane_integrated, bishop for circular, janbu for bilinear and structural.
search : {"seeded", "full"}, default="seeded", keyword-only
    seeded: the toe at the rupture daylight point (+/- toe_window_pixels);
    full: every toe down the path (up to max_length_metres).
back_scarp : {"vertical", "rankine"}, default="vertical", keyword-only
    The crown's back scarp: a tension crack, or a Rankine active plane
    carrying cohesion and friction (not for bilinear).
janbu_correction : bool, default=True, keyword-only
    Apply Janbu's f0.
structural_dip_gradient : float or Raster, default=0.3, keyword-only
    structural: tan(dip) of the weak layer, constant or per cell (at the crown).
structural_dip_direction_degrees : float or Raster, default=0.0, keyword-only
    structural: dip direction, degrees clockwise from north.
column_method : {"xie_revised", "hovland"}, default="xie_revised", keyword-only
    ellipsoid: the method of columns.
ellipsoid_search : {"seeded", "random"}, default="seeded", keyword-only
    ellipsoid: a grid of shapes, or random draws.
ellipsoid_random_draws : int, default=200, keyword-only
    ellipsoid random search: shapes per crown.
ellipsoid_random_seed : int, default=1, keyword-only
    ellipsoid random search: seed.
ellipsoid_min_half_axis_pixels : float, default=3.0, keyword-only
    ellipsoid: smallest length half-axis (cells); shorter crowns are "too small".
layer_depths : list of Raster or None, default=None, keyword-only
    ellipsoid: depths (m) of weak layers that truncate the ellipsoids; each
    is searched separately.
screen_safety_factor : float, default=1.5, keyword-only
    Crowns with rupture F above this are not evaluated.
toe_window_pixels : int, default=0, keyword-only
    Seeded search: toes this many cells either side of the daylight point.
max_length_metres : float, default=0.0, keyword-only
    Full search: longest crown-to-toe flow distance (0 = no cap).
slices_per_cell : int, default=4, keyword-only
    Slices per section cell.
search_points : int, default=12, keyword-only
    Points per shape parameter in the search.

Returns
-------
SlipSurfaces
    The critical surface of each evaluated crown, its bodies and counts.
)doc");

    nb::class_<SlipSurfaces>(module, "SlipSurfaces", R"doc(
Critical slip surfaces of each crown, from ``RupturePlanes.slip_surfaces``.

Members: ``family``, ``method``, ``safety_factor``, ``crown_depth``,
``max_depth``, ``sagitta``, ``basal_gradient``, ``length``, ``sinuosity``,
``unstable``, ``evaluated``, ``body_depth``, ``cover_safety_factor``,
``n_evaluated``, ``n_screened``, ``n_unstable``, ``n_inadmissible``,
``n_too_small``, ``failure_probability``, ``spencer_check``,
``body_catalogues``.
)doc")
        .def_prop_ro("family", [](const SlipSurfaces& s) { return s.family_name; }, "The slip-surface family.")
        .def_prop_ro("method", [](const SlipSurfaces& s) { return s.method_name; },
            "The limit-equilibrium method used (auto resolved), or the column method for ellipsoids.")
        .def_prop_ro("safety_factor", [](const SlipSurfaces& s) { return Raster(s.slip.safety_factor); },
            "Factor of safety of each evaluated crown's critical surface.")
        .def_prop_ro("crown_depth", [](const SlipSurfaces& s) { return Raster(s.slip.crown_depth); },
            "Depth h below the crown (m); ellipsoid: the depth half-axis.")
        .def_prop_ro("max_depth", [](const SlipSurfaces& s) { return Raster(s.slip.max_depth); },
            "Largest body thickness (m).")
        .def_prop_ro("sagitta", [](const SlipSurfaces& s) { return Raster(s.slip.sagitta); },
            "circular: the sagitta s (m); ellipsoid: the lateral half-axis (m).")
        .def_prop_ro("basal_gradient", [](const SlipSurfaces& s) { return Raster(s.slip.basal_gradient); },
            "bilinear: tan of the basal plane; structural: the apparent dip gradient.")
        .def_prop_ro("length", [](const SlipSurfaces& s) { return Raster(s.slip.length); },
            "Flow distance crown to toe (m); ellipsoid: straight distance.")
        .def_prop_ro("sinuosity", [](const SlipSurfaces& s) { return Raster(s.slip.sinuosity); },
            "2-D: flow distance over straight distance crown to toe.")
        .def_prop_ro("unstable", [](const SlipSurfaces& s) { return Raster(s.slip.unstable); }, "1 where F <= 1.")
        .def_prop_ro("evaluated", [](const SlipSurfaces& s) { return Raster(s.slip.evaluated); },
            "1 evaluated, 0 screened out (crowns with a rupture result).")
        .def_prop_ro("body_depth", [](const SlipSurfaces& s) { return Raster(s.slip.body_depth); },
            "Envelope of the unstable bodies' thickness (m).")
        .def_prop_ro("cover_safety_factor", [](const SlipSurfaces& s) { return Raster(s.slip.cover_safety_factor); },
            "ellipsoid: lowest F of the critical surfaces covering each cell.")
        .def_prop_ro("n_evaluated", [](const SlipSurfaces& s) { return s.slip.n_evaluated; }, "Crowns evaluated.")
        .def_prop_ro("n_screened", [](const SlipSurfaces& s) { return s.slip.n_screened; },
            "Crowns screened out by screen_safety_factor.")
        .def_prop_ro("n_unstable", [](const SlipSurfaces& s) { return s.slip.n_unstable; }, "Crowns with F <= 1.")
        .def_prop_ro("n_inadmissible", [](const SlipSurfaces& s) { return s.slip.n_inadmissible; },
            "Crowns with no admissible surface.")
        .def_prop_ro("n_too_small", [](const SlipSurfaces& s) { return s.ell.n_too_small; },
            "ellipsoid: crowns below ellipsoid_min_half_axis_pixels (0 for the 2-D families).")
        .def("failure_probability", &slip_failure_probability, LSDTT_UNCERTAINTY_ARGS,
            nb::sig("def failure_probability(self, " LSDTT_UNCERTAINTY_SIG ") -> dict"),
            "Probability of failure on each crown's fixed critical surface.\n\n"
            "lsdtt-slope-stability ``parameter_uncertainty`` with a slip family. The surface found with the\n"
            "central values is held fixed and evaluated for each (c, mu) sample, so Pf is a lower bound on the\n"
            "probability that some surface fails. Screened-out crowns get 0. For ellipsoids the result also\n"
            "has ``cover_failure_probability``, the highest Pf of the ellipsoids covering each cell.\n\n"
            "Parameters\n----------\n" LSDTT_UNCERTAINTY_DOC "\nReturns\n-------\n" LSDTT_PROBABILITY_RETURNS)
        .def("spencer_check", &slip_spencer_check, nb::sig("def spencer_check(self) -> dict"), R"doc(
Spencer's factor of safety on each critical 2-D surface.

lsdtt-slope-stability ``slip_surface_spencer_check``: force and moment
equilibrium with a constant interslice inclination, on the surface the
family's own method found, as a check of the simplified methods. Not for
ellipsoids.

Parameters
----------
self : SlipSurfaces
    A 2-D family's result.

Returns
-------
dict
    ``safety_factor`` (Raster, nodata where Spencer has no solution) and
    ``n_no_solution``.
)doc")
        .def("body_catalogues", &slip_body_catalogues, nb::sig("def body_catalogues(self) -> dict"), R"doc(
Catalogues of the unstable slip bodies: by toe and amalgamated.

lsdtt-slope-stability ``write_slip_body_catalogue``. By toe: the union of the
bodies of the unstable crowns sharing a toe (groups may overlap, so their
volumes need not add up). Amalgamated: 8-connected components of the body
envelope, which do not overlap.

Parameters
----------
self : SlipSurfaces
    The slip surfaces.

Returns
-------
dict
    ``by_toe`` and ``amalgamated`` (DataFrames, the CLI's
    ``_slip_body_catalogue*.csv`` columns) and ``label_amalgamated`` (Raster).
)doc");

    flow.def("rupture_planes", &rupture_planes, nb::arg("geometry"), nb::kw_only(), nb::arg("channel") = nb::none(),
        nb::arg("threshold_contributing_pixels") = 1000, nb::arg("cohesion_pa") = 1000.0,
        nb::arg("friction_coefficient") = 0.70020753820970977, nb::arg("density_kg_m3") = 2600.0,
        nb::arg("region") = nb::none(), nb::arg("region_parameters") = nb::none(),
        nb::arg("water_table_depth") = nb::none(), nb::arg("n_threads") = 0,
        nb::sig("def rupture_planes(self, geometry: Raster, *, channel: Raster | None = None, "
                "threshold_contributing_pixels: int = 1000, cohesion_pa: float = 1000.0, "
                "friction_coefficient: float = 0.7002075382097097, density_kg_m3: float = 2600.0, "
                "region: Raster | None = None, region_parameters: dict | None = None, "
                "water_table_depth: Raster | None = None, n_threads: int = 0) -> RupturePlanes"), R"doc(
Planar rupture surfaces of every hillslope cell (scaLr, Steer et al.).

lsdtt-slope-stability ``run_rupture_planes``, the same routine; the keyword
names are the CLI's less the ``rupture_`` prefix. From each source cell R the
D8 path (this FlowInfo's receivers) is walked down to the channel; for each
candidate daylight point D the plane through D at the optimal depth h*
(closed form) gives F, and the minimum over the path is kept. With a water
table the pore pressure is hydrostatic and the minimum is found numerically.
The CLI routes on the carved-then-filled DEM and measures on the input DEM.

Parameters
----------
geometry : Raster
    Elevations for dx and dz, normally the input (unconditioned) DEM.
channel : Raster or None, default=None, keyword-only
    Nonzero = channel cell (the walk stops before it). None uses
    ``threshold_contributing_pixels`` on this FlowInfo, as the CLI.
threshold_contributing_pixels : int, default=1000, keyword-only
    Channel cells have at least this many D8 contributing pixels.
cohesion_pa : float, default=1000.0, keyword-only
    Cohesion c (Pa).
friction_coefficient : float, default=tan 35, keyword-only
    mu = tan phi.
density_kg_m3 : float, default=2600.0, keyword-only
    Bulk density rho.
region : Raster or None, default=None, keyword-only
    Integer region ids selecting parameter sets; nodata cells are left out.
region_parameters : dict or None, default=None, keyword-only
    ``{region_id: {name: value}}`` with the CLI's region-table column names
    (``rupture_cohesion_pa``, ``rupture_cohesion_pa_min``, ...); a missing name
    takes the keyword value; ``sinmap_*`` names are ignored. Needed with
    ``region``.
water_table_depth : Raster or None, default=None, keyword-only
    Depth (m) of the water table below the geometry surface, e.g.
    ``GroundwaterDupuit.water_table_depth()``; None or nodata is dry.
n_threads : int, default=0, keyword-only
    OpenMP threads; the result does not depend on it.

Returns
-------
RupturePlanes
    The rasters and counts, and methods for the probability of failure,
    catalogues and slip surfaces.
)doc");

    module.def("coseismic_hazard_area", &coseismic_hazard_area, nb::arg("routing_surface"), nb::kw_only(),
        nb::arg("initiation_slope") = 0.839, nb::arg("stopping_slope") = 0.176,
        nb::arg("boundary_conditions") = "oooo", nb::arg("n_threads") = 0,
        nb::sig("def coseismic_hazard_area(routing_surface: Raster, *, initiation_slope: float = 0.839, "
                "stopping_slope: float = 0.176, boundary_conditions: str = 'oooo', n_threads: int = 0) -> dict"), R"doc(
Coseismic landslide hazard area (SHALRUN-EQ; Milledge et al. 2019).

lsdtt-slope-stability ``run_coseismic_hazard_area``, the same routine. The hazard
area of a cell is the area of cells steep enough to initiate whose debris can
reach it by slope-weighted multiple flow directions without crossing a cell
gentler than the stopping slope, per metre of cell width. A gentler cell is
reached and counted, but passes nothing on. It is a relative index of exposure.

Parameters
----------
routing_surface : Raster
    A pit-free DEM, for example carved and filled.
initiation_slope : float, default=0.839, keyword-only
    Gradient (m/m) at or above which a cell can initiate a landslide
    (CLI ``coseismic_hazard_initiation_slope``).
stopping_slope : float, default=0.176, keyword-only
    Gradient (m/m) below which run-out stops; above 0 and below
    ``initiation_slope`` (CLI ``coseismic_hazard_stopping_slope``).
boundary_conditions : str, default="oooo", keyword-only
    Flow boundary conditions (north, east, south, west).
n_threads : int, default=0, keyword-only
    OpenMP threads for building the flow graph.

Returns
-------
dict
    ``hazard_area`` (Raster, m^2/m), ``steepest_slope`` (Raster, m/m), ``zone``
    (Raster: 0 transit, 1 source, 2 stopping), ``edge_influenced`` (Raster: 1
    where debris from a DEM-edge cell arrives, so the value may be a lower
    bound), ``n_valid_cells``, ``n_source_cells``, ``n_stopping_cells``,
    ``n_zero_hazard_cells``, ``n_edge_influenced_cells``, ``source_area_m2``.
)doc");

    module.def("dem_change_landslides", &dem_change_landslides, nb::arg("difference"), nb::arg("dem"), nb::kw_only(),
        nb::arg("stable_mask") = nb::none(), nb::arg("is_newer_minus_older") = true, nb::arg("change_sign") = "loss",
        nb::arg("detection_sigma") = 1.96, nb::arg("seed_sigma") = 3.0, nb::arg("patch_significance") = 1.96,
        nb::arg("min_patch_pixels") = 5, nb::arg("aggregation_radius_pixels") = 2, nb::arg("min_slope") = 0.0,
        nb::arg("n_slope_bins") = 10, nb::arg("n_null_simulations") = 20, nb::arg("random_seed") = 33,
        nb::arg("n_threads") = 0,
        nb::sig("def dem_change_landslides(difference: Raster, dem: Raster, *, stable_mask: Raster | None = None, "
                "is_newer_minus_older: bool = True, change_sign: str = 'loss', detection_sigma: float = 1.96, "
                "seed_sigma: float = 3.0, patch_significance: float = 1.96, min_patch_pixels: int = 5, "
                "aggregation_radius_pixels: int = 2, min_slope: float = 0.0, n_slope_bins: int = 10, "
                "n_null_simulations: int = 20, random_seed: int = 33, n_threads: int = 0) -> dict"), R"doc(
Potential landslides in a DEM of difference, with their uncertainty.

lsdtt-slope-stability ``run_dem_change_landslides``, the same routine; the
keyword names are the CLI's less the ``dem_change_`` and ``dem_difference_``
prefixes. An error model from stable terrain (bias, NMAD by slope, spatial
correlation) sets a threshold for each cell. Connected cells are grouped and
aggregated into patches, each tested against the correlated error of its mean.
Noise simulations estimate false detections. The results are CANDIDATE areas of
change: vegetation, snow, water, radar penetration and DEM artefacts also change
elevation between two products.

Parameters
----------
difference : Raster
    DEM of difference (m), for example from lsdtt-raster-calculator.
dem : Raster
    The newer DEM on the same grid, used for slope.
stable_mask : Raster or None, default=None, keyword-only
    Non-zero on stable terrain. Without it every assessable cell within 4 NMAD
    of the median is treated as stable.
is_newer_minus_older : bool, default=True, keyword-only
    False if ``difference`` is older minus newer.
change_sign : str, default="loss", keyword-only
    ``loss``, ``gain`` or ``both``.
detection_sigma : float, default=1.96, keyword-only
    Multiple of the slope-dependent error at which a cell can join a patch.
seed_sigma : float, default=3.0, keyword-only
    A connected group of candidates needs one cell at this multiple.
patch_significance : float, default=1.96, keyword-only
    Minimum absolute mean change over its standard error.
min_patch_pixels : int, default=5, keyword-only
    Smallest patch.
aggregation_radius_pixels : int, default=2, keyword-only
    Join changed cells with up to this many unchanged cells between them.
min_slope : float, default=0.0, keyword-only
    Ignore cells with a gradient (m/m) below this.
n_slope_bins : int, default=10, keyword-only
    Quantile bins of slope for the error model.
n_null_simulations : int, default=20, keyword-only
    Noise simulations for false detections; 0 turns them off.
random_seed : int, default=33, keyword-only
    Seed for the variogram sample and the simulations.
n_threads : int, default=0, keyword-only
    OpenMP threads for the simulations; the result does not depend on it.

Returns
-------
dict
    ``uncertainty`` and ``false_detection`` (dicts); DataFrames
    ``uncertainty_by_slope``, ``variogram``, ``patches``, ``coverage``,
    ``sensitivity`` and ``false_detection_by_area``; and Rasters ``corrected``,
    ``sigma``, ``z``, ``change_class`` (-1 loss, 1 gain) and ``label``.
)doc");

    module.def("landslide_frequency_area", &landslide_frequency_area, nb::arg("areas_m2"), nb::arg("cell_area_m2"),
        nb::kw_only(), nb::arg("bins_per_decade") = 10, nb::arg("min_tail") = 50,
        nb::sig("def landslide_frequency_area(areas_m2: Sequence[float], cell_area_m2: float, *, "
                "bins_per_decade: int = 10, min_tail: int = 50) -> tuple[pandas.DataFrame, dict]"), R"doc(
Frequency - area of a landslide catalogue and its power-law tail.

lsdtt-slope-stability's ``_frequency_area.csv`` and ``_frequency_area_fit.csv``
for one catalogue (landslides or slip bodies). Areas are whole numbers of
cells, so the fit and the bins are discrete in cells. The tail p(A) ~ A^eta
is fitted by maximum likelihood with A_min chosen by the Kolmogorov -
Smirnov distance (Clauset, Shalizi and Newman 2009); ``binned_eta`` is the
least-squares slope of log p(A) against log A over the bins at or above
A_min holding at least 3 landslides.

Parameters
----------
areas_m2 : sequence of float
    Landslide areas (m^2), e.g. a catalogue's ``area_m2`` column.
cell_area_m2 : float
    Area of one cell (m^2).
bins_per_decade : int, default=10, keyword-only
    Logarithmic bins per decade of area (CLI ``frequency_area_bins_per_decade``).
min_tail : int, default=50, keyword-only
    Fewest landslides in the fitted tail (at least 2; CLI ``frequency_area_min_tail``).

Returns
-------
tuple of (pandas.DataFrame, dict)
    The density (``area_lower_m2, area_upper_m2, area_centre_m2,
    n_landslides, probability_density_per_m2``) and the fit (``n_landslides,
    method, area_min_m2, eta, eta_standard_error, n_tail, ks_distance,
    binned_eta, binned_eta_standard_error, n_bins_fitted``), the CLI's columns
    less ``catalogue``.
)doc");

    nb::class_<GroundwaterDupuit>(module, "GroundwaterDupuit", R"doc(
An unconfined aquifer on an impermeable base: the Dupuit - Boussinesq model.

lsdtt-slope-stability ``water_table_method: dupuit_steady`` and
``dupuit_transient`` (after Landlab's GroundwaterDupuitPercolator, Litwin et
al. 2020). The base lies ``regolith_thickness_metres`` below the surface;
water flows down the water-table gradient; ``fixed_head`` cells (the CLI's
channels) are seepage outlets with the water table at the surface; the
raster edge and nodata are no-flux. Where the regolith saturates, water
seeps out. The keyword names are the CLI's less the ``groundwater_`` prefix.

Members: ``run_to_steady``, ``run``, ``set_thickness``,
``water_table_depth``, ``saturated_thickness``, ``seepage``, ``storage``,
``budget``, ``steady_rates``.

Parameters
----------
surface : Raster
    Ground elevation (m), normally the input DEM.
fixed_head : Raster
    Nonzero = fixed-head (channel) cell, e.g. ``RupturePlanes.channel``.
regolith_thickness_metres : float, default=2.0, keyword-only
    Aquifer thickness d (m), used where ``regolith_thickness`` is None.
regolith_thickness : Raster or None, default=None, keyword-only
    Thickness per cell (m); nodata or <= 0 cells are left out.
hydraulic_conductivity_m_per_s : float, default=0.001, keyword-only
    K, or K0 at the surface for the exponential profile.
drainable_porosity : float, default=0.2, keyword-only
    n.
transmissivity_profile : {"uniform", "exponential"}, default="uniform", keyword-only
    exponential: K decays with depth over conductivity_decay_depth_metres
    (TOPMODEL's assumption).
conductivity_decay_depth_metres : float, default=1.0, keyword-only
    z_m of the exponential profile.
n_threads : int, default=0, keyword-only
    OpenMP threads; the result does not depend on it.
)doc")
        .def(nb::init<const Raster&, const Raster&, double, nb::object, double, double, const std::string&, double,
                      int>(),
            nb::arg("surface"), nb::arg("fixed_head"), nb::kw_only(), nb::arg("regolith_thickness_metres") = 2.0,
            nb::arg("regolith_thickness") = nb::none(), nb::arg("hydraulic_conductivity_m_per_s") = 1.0e-3,
            nb::arg("drainable_porosity") = 0.2, nb::arg("transmissivity_profile") = "uniform",
            nb::arg("conductivity_decay_depth_metres") = 1.0, nb::arg("n_threads") = 0,
            nb::sig("def __init__(self, surface: Raster, fixed_head: Raster, *, "
                    "regolith_thickness_metres: float = 2.0, regolith_thickness: Raster | None = None, "
                    "hydraulic_conductivity_m_per_s: float = 0.001, drainable_porosity: float = 0.2, "
                    "transmissivity_profile: str = 'uniform', conductivity_decay_depth_metres: float = 1.0, "
                    "n_threads: int = 0) -> None"))
        .def("run_to_steady", &GroundwaterDupuit::run_to_steady, nb::arg("recharge_m_per_s"), nb::kw_only(),
            nb::arg("steady_method") = "implicit", nb::arg("steady_tolerance") = 1.0e-6,
            nb::arg("max_substeps") = 10000000, nb::arg("max_newton_iterations") = 200,
            nb::sig("def run_to_steady(self, recharge_m_per_s: float | Raster, *, steady_method: str = 'implicit', "
                    "steady_tolerance: float = 1e-06, max_substeps: int = 10000000, "
                    "max_newton_iterations: int = 200) -> dict"), R"doc(
Find the steady water table under a constant recharge.

lsdtt-slope-stability ``dupuit_steady`` (and the antecedent state of
``dupuit_transient``). Raises RuntimeError if it does not converge.

Parameters
----------
recharge_m_per_s : float or Raster
    Recharge f (m/s), > 0 somewhere; nodata is 0.
steady_method : {"implicit", "explicit"}, default="implicit", keyword-only
    implicit: solved directly by Newton's method (pseudo-transient
    continuation with an active-set step for saturated and dry cells); tens
    of iterations whatever the grid. explicit: the time-dependent model
    marched until it stops changing; its steps are limited by the Courant
    and von Neumann conditions. Both solve the same discrete equations.
steady_tolerance : float, default=1e-06, keyword-only
    Converged when the residual is at most this times the largest recharge.
max_substeps : int, default=10000000, keyword-only
    explicit: most substeps.
max_newton_iterations : int, default=200, keyword-only
    implicit: most Newton (pseudo-time) steps.

Returns
-------
dict
    ``steady_method``, ``steps``, ``residual``, ``pseudo_time_s``; implicit
    also ``linear_iterations``, ``used_direct_solver``, ``n_saturated``,
    ``n_dry``. The implicit solve's budget is in ``steady_rates()``.
)doc")
        .def("run", &GroundwaterDupuit::run, nb::arg("duration_s"), nb::arg("recharge_m_per_s"), nb::kw_only(),
            nb::arg("transient_method") = "explicit", nb::arg("implicit_max_dt_s") = 3600.0,
            nb::arg("implicit_dt_tolerance_metres") = 0.0,
            nb::sig("def run(self, duration_s: float, recharge_m_per_s: float | Raster, *, "
                    "transient_method: str = 'explicit', implicit_max_dt_s: float = 3600.0, "
                    "implicit_dt_tolerance_metres: float = 0.0) -> dict"), R"doc(
Advance the water table by ``duration_s`` under a constant recharge.

lsdtt-slope-stability ``dupuit_transient`` runs its recharge series as one
call per constant segment. ``budget()`` accumulates over calls.

Parameters
----------
duration_s : float
    Time to advance (s), > 0.
recharge_m_per_s : float or Raster
    Recharge f (m/s); nodata is 0.
transient_method : {"explicit", "implicit"}, default="explicit", keyword-only
    explicit: substeps at the stability limit (set by stability, not
    accuracy: on coarse grids they can be hours long). implicit: backward
    Euler, each step solved by Newton's method; first order in the step and
    damps fast transients.
implicit_max_dt_s : float, default=3600.0, keyword-only
    implicit: longest step (s).
implicit_dt_tolerance_metres : float, default=0.0, keyword-only
    implicit: step-doubling error tolerance on the water table (m); 0 = off.

Returns
-------
dict
    ``transient_method``, ``steps``; implicit also ``newton_iterations``,
    ``linear_iterations``, ``splits``, ``rejected_by_error``,
    ``max_error_estimate_metres``.
)doc")
        .def("set_thickness", &GroundwaterDupuit::set_thickness, nb::arg("thickness"),
            nb::sig("def set_thickness(self, thickness: Raster) -> None"), R"doc(
Set the saturated thickness h of the active cells (clipped to [0, d]).

Parameters
----------
thickness : Raster
    Saturated thickness (m); nodata is 0.

Returns
-------
None
)doc")
        .def("water_table_depth", &GroundwaterDupuit::water_table_depth,
            nb::sig("def water_table_depth(self) -> Raster"), R"doc(
Depth of the water table below the surface (m): the input to ``FlowInfo.rupture_planes``.

Parameters
----------
self : GroundwaterDupuit
    The aquifer.

Returns
-------
Raster
    Water-table depth (m); 0 where saturated; nodata outside the aquifer.
)doc")
        .def("saturated_thickness", &GroundwaterDupuit::saturated_thickness,
            nb::sig("def saturated_thickness(self) -> Raster"), R"doc(
Saturated thickness h (m).

Parameters
----------
self : GroundwaterDupuit
    The aquifer.

Returns
-------
Raster
    h (m); nodata outside the aquifer.
)doc")
        .def("seepage", &GroundwaterDupuit::seepage, nb::sig("def seepage(self) -> Raster"), R"doc(
Seepage rate (m/s) from the last call.

The implicit steady state's seepage, or the mean over the last ``run``.

Parameters
----------
self : GroundwaterDupuit
    The aquifer.

Returns
-------
Raster
    Seepage (m/s); nodata outside the aquifer.
)doc")
        .def_prop_ro("storage", &GroundwaterDupuit::storage, "Water stored, n h summed over the active cells (m^3).")
        .def("budget", &GroundwaterDupuit::budget, nb::sig("def budget(self) -> dict"), R"doc(
The cumulative water budget of the explicit steady march and the runs (m^3).

recharge_m3 + clipped_added_m3 = storage_change_m3 + outflow_fixed_head_m3 +
seepage_m3 to round-off. An implicit steady solve reports ``steady_rates()``
instead and does not add to this.

Parameters
----------
self : GroundwaterDupuit
    The aquifer.

Returns
-------
dict
    ``recharge_m3, storage_change_m3, outflow_fixed_head_m3, seepage_m3,
    clipped_added_m3``.
)doc")
        .def("steady_rates", &GroundwaterDupuit::steady_rates, nb::sig("def steady_rates(self) -> dict"), R"doc(
The budget of the last implicit steady solve, as rates (m^3/s).

Parameters
----------
self : GroundwaterDupuit
    The aquifer.

Returns
-------
dict
    ``recharge_m3_per_s, outflow_fixed_head_m3_per_s, seepage_m3_per_s,
    imbalance_m3_per_s``.
)doc");
}

} // namespace lsdtt3_python
