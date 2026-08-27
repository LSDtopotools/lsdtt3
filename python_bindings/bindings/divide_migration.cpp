#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <map>
#include <optional>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>

#include <cpu/lsdtt_basin_outline_cpu.hpp>
#include <cpu/lsdtt_chi_extraction_cpu.hpp>
#include <cpu/lsdtt_divide_migration_cpu.hpp>
#include <cpu/lsdtt_flow_snapping_cpu.hpp>
#include <cpu/lsdtt_junction_basin_manager_cpu.hpp>

#include "binding_common.hpp"

namespace lsdtt3_python {
namespace {

nb::object dataframe(nb::list rows)
{
    return nb::module_::import_("pandas").attr("DataFrame")(rows);
}

nb::object geodataframe(nb::object data, nb::list geometry,
                        const lsdtt::lsdtt_raster_meta& meta)
{
    nb::object maker = nb::module_::import_("geopandas").attr("GeoDataFrame");
    if (meta.crs_wkt.empty()) return maker(data, nb::arg("geometry") = geometry);
    return maker(data, nb::arg("geometry") = geometry, nb::arg("crs") = meta.crs_wkt);
}

class DivideMigrationResult {
public:
    DivideMigrationResult(std::shared_ptr<lsdtt::lsdtt_flowinfo> flow,
                          ChannelNetwork network, ChiCollection chi, Raster labels,
                          std::vector<lsdtt::lsdtt_outlet_record> outlets,
                          std::vector<lsdtt::lsdtt_unit_record> units,
                          std::vector<lsdtt::lsdtt_divide_segment_record> divides,
                          std::vector<lsdtt::lsdtt_channel_head_metric_record> heads,
                          std::vector<lsdtt::lsdtt_pair_summary_record> summaries)
        : flow_(std::move(flow)), network_(std::move(network)), chi_(std::move(chi)),
          labels_(std::move(labels)), outlets_(std::move(outlets)), units_(std::move(units)),
          divides_(std::move(divides)), heads_(std::move(heads)), summaries_(std::move(summaries)) {}

    const ChannelNetwork& network() const noexcept { return network_; }
    const ChiCollection& chi() const noexcept { return chi_; }
    const Raster& basin_labels() const noexcept { return labels_; }
    const std::vector<lsdtt::lsdtt_outlet_record>& outlets() const noexcept { return outlets_; }
    const std::vector<lsdtt::lsdtt_unit_record>& units() const noexcept { return units_; }
    const std::vector<lsdtt::lsdtt_divide_segment_record>& divides() const noexcept { return divides_; }
    const std::vector<lsdtt::lsdtt_channel_head_metric_record>& heads() const noexcept { return heads_; }
    const std::vector<lsdtt::lsdtt_pair_summary_record>& summaries() const noexcept { return summaries_; }
private:
    std::shared_ptr<lsdtt::lsdtt_flowinfo> flow_;
    ChannelNetwork network_;
    ChiCollection chi_;
    Raster labels_;
    std::vector<lsdtt::lsdtt_outlet_record> outlets_;
    std::vector<lsdtt::lsdtt_unit_record> units_;
    std::vector<lsdtt::lsdtt_divide_segment_record> divides_;
    std::vector<lsdtt::lsdtt_channel_head_metric_record> heads_;
    std::vector<lsdtt::lsdtt_pair_summary_record> summaries_;
};

nb::object outlet_table(const DivideMigrationResult& result)
{
    nb::list rows;
    for (const auto& r : result.outlets()) {
        nb::dict row;
        row["outlet_id"] = r.point_id; row["input_row"] = r.row; row["input_col"] = r.col;
        row["input_flat"] = r.input_flat; row["snapped_flat"] = r.snapped_flat;
        row["snapped"] = r.snapped; row["unit_id"] = r.unit_id;
        row["unit_key"] = r.unit_key; row["cluster"] = r.cluster_value;
        rows.append(std::move(row));
    }
    return dataframe(rows);
}

nb::object unit_table(const DivideMigrationResult& result)
{
    nb::list rows;
    for (const auto& u : result.units()) {
        nb::dict row; row["unit_id"] = u.unit_id; row["unit_key"] = u.unit_key;
        row["outlet_ids"] = u.point_ids; row["snapped_outlet_flats"] = u.snapped_outlet_flats;
        rows.append(std::move(row));
    }
    return dataframe(rows);
}

nb::object divide_gdf(const DivideMigrationResult& result)
{
    nb::list rows, geometry;
    nb::object line = nb::module_::import_("shapely.geometry").attr("LineString");
    for (const auto& d : result.divides()) {
        nb::dict row; row["pair_id"] = d.pair_id; row["unit_a"] = d.unit_a; row["unit_b"] = d.unit_b;
        row["junction_a"] = d.basin_junction_a; row["junction_b"] = d.basin_junction_b;
        row["segment_id"] = d.segment_id; row["n_fragments"] = d.n_fragments;
        row["length_metres"] = d.length_m; rows.append(std::move(row));
        nb::list xy; for (const auto& p : d.points) xy.append(nb::make_tuple(p.x, p.y));
        geometry.append(line(xy));
    }
    return geodataframe(dataframe(rows), geometry, result.network().flow().meta());
}

nb::object head_gdf(const DivideMigrationResult& result)
{
    nb::list rows, geometry;
    nb::object point = nb::module_::import_("shapely.geometry").attr("Point");
    for (const auto& h : result.heads()) {
        nb::dict row; row["head_id"] = h.head_id; row["outlet_id"] = h.outlet_id; row["unit_id"] = h.unit_id;
        row["pair_id_nearest"] = h.pair_id_nearest; row["head_x"] = h.head_x; row["head_y"] = h.head_y;
        row["head_elevation_m"] = h.head_elevation_m;
        row["chi"] = h.chi; row["nearest_divide_x"] = h.nearest_divide_x; row["nearest_divide_y"] = h.nearest_divide_y;
        row["divide_elevation_m"] = h.divide_elevation_m; row["distance_to_divide_m"] = h.distance_to_divide_m;
        row["relief_m"] = h.relief_m; row["gradient_to_divide"] = h.gradient_to_divide;
        rows.append(std::move(row)); geometry.append(point(h.head_x, h.head_y));
    }
    return geodataframe(dataframe(rows), geometry, result.network().flow().meta());
}

nb::object summary_table(const DivideMigrationResult& result)
{
    nb::list rows;
    for (const auto& s : result.summaries()) {
        nb::dict r; r["pair_id"] = s.pair_id; r["unit_a"] = s.unit_a; r["unit_b"] = s.unit_b;
        r["n_heads_a"] = s.n_heads_a; r["n_heads_b"] = s.n_heads_b;
        r["mean_head_elevation_a"] = s.mean_head_elevation_a; r["mean_head_elevation_b"] = s.mean_head_elevation_b; r["delta_head_elevation"] = s.delta_head_elevation;
        r["mean_gradient_a"] = s.mean_gradient_a; r["mean_gradient_b"] = s.mean_gradient_b; r["delta_gradient"] = s.delta_gradient;
        r["mean_relief_a"] = s.mean_relief_a; r["mean_relief_b"] = s.mean_relief_b; r["delta_relief"] = s.delta_relief;
        r["mean_chi_a"] = s.mean_chi_a; r["mean_chi_b"] = s.mean_chi_b; r["delta_chi"] = s.delta_chi;
        r["uncertainty_method"] = s.uncertainty_method; r["uncertainty_head_elevation"] = s.uncertainty_head_elevation;
        r["uncertainty_gradient"] = s.uncertainty_gradient; r["uncertainty_relief"] = s.uncertainty_relief; r["uncertainty_chi"] = s.uncertainty_chi;
        rows.append(std::move(r));
    }
    return dataframe(rows);
}

std::vector<std::vector<std::int64_t>> outlet_pairs(nb::object outlets)
{
    nb::object a = nb::module_::import_("numpy").attr("asarray")(outlets, nb::arg("dtype") = "int64");
    if (nb::cast<int>(a.attr("ndim")) != 2 || nb::cast<int>(a.attr("shape")[1]) != 2)
        throw std::invalid_argument("outlets must be an N by 2 (row, col) array");
    return nb::cast<std::vector<std::vector<std::int64_t>>>(a.attr("tolist")());
}

std::vector<std::string> cluster_labels(nb::object clusters, std::size_t n)
{
    if (clusters.is_none()) return {};
    nb::object values = nb::module_::import_("numpy").attr("asarray")(clusters).attr("ravel")().attr("tolist")();
    nb::list list = nb::cast<nb::list>(values);
    if (list.size() != n) throw std::invalid_argument("clusters must have one value per outlet");
    std::vector<std::string> result; result.reserve(n);
    for (nb::handle value : list) result.push_back(nb::cast<std::string>(nb::str(value)));
    return result;
}

DivideMigrationResult divide_migration(const FlowInfo& flow_owner, const Raster& dem,
                                       nb::object outlets, nb::object clusters, int source_pixels,
                                       const std::string& snapping, const std::string& outlet_control,
                                       double min_outlet_elevation, double m_over_n, double reference_area,
                                       const std::string& uncertainty, int bootstrap_replicates, std::uint32_t seed)
{
    const auto& flow = flow_owner.native();
    if (source_pixels < 1 || m_over_n <= 0.0 || reference_area <= 0.0)
        throw std::invalid_argument("source_pixels, m_over_n, and reference_area must be positive");
    if (dem.native().get_meta().n_rows != flow.n_rows() || dem.native().get_meta().n_cols != flow.n_cols())
        throw std::invalid_argument("dem must share the FlowInfo grid (use the conditioned DEM)");
    const auto pairs = outlet_pairs(outlets); const auto labels = cluster_labels(clusters, pairs.size());
    if (pairs.size() < 2) throw std::invalid_argument("at least two outlets are required");
    if (snapping != "downstream" && snapping != "none") throw std::invalid_argument("snapping must be 'downstream' or 'none'");
    if (outlet_control != "none" && outlet_control != "complete_only" && outlet_control != "elevation" && outlet_control != "max_out_elevation")
        throw std::invalid_argument("outlet_control must be 'none', 'complete_only', 'elevation', or 'max_out_elevation'");

    std::map<std::string, int> unit_for_key; std::map<std::int64_t, int> unit_for_flat, outlet_for_flat;
    std::vector<lsdtt::lsdtt_outlet_record> outlet_records; std::vector<std::int64_t> seed_flats;
    int next_unit = 1;
    lsdtt::lsdtt_downstream_snap_options snap_options; snap_options.contributing_pixels_threshold = source_pixels;
    for (std::size_t i = 0; i < pairs.size(); ++i) {
        const int row = static_cast<int>(pairs[i][0]), col = static_cast<int>(pairs[i][1]);
        if (!flow.in_bounds(row, col) || !flow.is_valid_cell(row, col)) throw std::invalid_argument("outlets contains an invalid flow pixel");
        lsdtt::lsdtt_outlet_record r; r.point_id = static_cast<int>(i + 1); r.row = row; r.col = col;
        r.input_flat = flow.flat_index(row, col); r.snapped_flat = r.input_flat; r.snapped = snapping == "none";
        if (snapping == "downstream") {
            const auto s = lsdtt::snap_pixel_to_contributing_pixels_threshold_downstream(flow, row, col, snap_options);
            if (s.status != lsdtt::lsdtt_snap_status::success) { outlet_records.push_back(std::move(r)); continue; }
            r.snapped_flat = s.snapped_flat; r.snapped = true;
        }
        r.cluster_value = labels.empty() ? "" : labels[i]; r.unit_key = labels.empty() ? "outlet:" + std::to_string(r.point_id) : "cluster:" + r.cluster_value;
        auto [it, inserted] = unit_for_key.emplace(r.unit_key, next_unit);
        if (inserted) ++next_unit; r.unit_id = it->second;
        const auto existing = unit_for_flat.find(r.snapped_flat);
        if (existing != unit_for_flat.end() && existing->second != r.unit_id)
            throw std::invalid_argument("outlets from different units snap to the same pixel");
        unit_for_flat[r.snapped_flat] = r.unit_id; outlet_for_flat[r.snapped_flat] = r.point_id;
        seed_flats.push_back(r.snapped_flat); outlet_records.push_back(std::move(r));
    }
    if (unit_for_key.size() < 2) throw std::invalid_argument("fewer than two distinct snapped units remain");

    if (outlet_control == "elevation") {
        std::vector<std::int64_t> kept;
        for (const auto flat : seed_flats)
            if (dem.native().get_data()[static_cast<std::size_t>(flat)] >= min_outlet_elevation) kept.push_back(flat);
        seed_flats.swap(kept);
    }
    const auto sources = flow.get_sources_index_threshold(source_pixels);
    ChannelNetwork network(flow_owner.shared_native(), sources);
    lsdtt::lsdtt_junction_basin_manager_cpu manager(flow, network.native());
    manager.select_junction_basins_from_basin_outlets(seed_flats);
    manager.filter_duplicates();
    if (outlet_control == "complete_only") manager.filter_incomplete();
    const auto basin_outlets = manager.get_basin_outlet_flat_indices();
    if (basin_outlets.size() < 2) throw std::invalid_argument("fewer than two selected basins remain");

    std::vector<lsdtt::lsdtt_flat_index_slice> cells; cells.reserve(manager.size());
    for (std::size_t i = 0; i < manager.size(); ++i) cells.push_back(manager.basin_at(i).get_slice());
    const auto adjacent = lsdtt::find_adjacent_basin_pairs(cells, flow.n_rows(), flow.n_cols());
    std::vector<lsdtt::lsdtt_divide_segment_record> divides; int segment_id = 1;
    for (const auto [ia, ib] : adjacent) {
        const auto& a = manager.basin_at(static_cast<std::size_t>(ia)); const auto& b = manager.basin_at(static_cast<std::size_t>(ib));
        const auto flat_a = a.get_basin_outlet_flat_index(), flat_b = b.get_basin_outlet_flat_index();
        const int unit_a = unit_for_flat.at(flat_a), unit_b = unit_for_flat.at(flat_b);
        if (unit_a == unit_b || flow.is_upslope_of(flat_a, flat_b) || flow.is_upslope_of(flat_b, flat_a)) continue;
        const auto ridge = a.get_basin_outline().get_interfluvial_divide_as_ridgeline(a.get_junction_number(), b.get_junction_number(), network.native());
        if (ridge.empty()) continue;
        lsdtt::lsdtt_divide_segment_record d; d.unit_a = std::min(unit_a, unit_b); d.unit_b = std::max(unit_a, unit_b);
        d.pair_id = std::to_string(d.unit_a) + "|" + std::to_string(d.unit_b); d.basin_junction_a = a.get_junction_number(); d.basin_junction_b = b.get_junction_number();
        d.segment_id = segment_id++; d.points = ridge.positions(); d.length_m = ridge.total_length(); divides.push_back(std::move(d));
    }
    divides = lsdtt::merge_divide_fragments(divides, flow.cell_size());

    lsdtt::lsdtt_chi_options chi_options; chi_options.A_0 = reference_area; chi_options.m_over_n = m_over_n;
    ChiCollection chi(lsdtt::get_chi_from_outlets_zero_base_threshold(flow, basin_outlets, source_pixels, chi_options), flow_owner.shared_native());
    if (outlet_control == "max_out_elevation") {
        double base = -std::numeric_limits<double>::infinity();
        for (const auto& basin : chi.native().basins)
            base = std::max(base, static_cast<double>(dem.native().get_data()[static_cast<std::size_t>(basin.outlet_flat)]));
        if (std::isfinite(base))
            lsdtt::apply_common_base_level_to_chi(chi.mutable_native(), dem.native(), base);
    }
    std::vector<lsdtt::lsdtt_channel_head_metric_record> heads; int head_id = 1;
    for (const auto& basin : chi.native().basins) {
        const int unit = unit_for_flat.at(basin.outlet_flat); std::set<std::int64_t> receivers(basin.receiver_flats.begin(), basin.receiver_flats.end());
        for (std::size_t i = 0; i < basin.flats.size(); ++i) {
            const auto flat = basin.flats[i]; if (receivers.count(flat)) continue;
            const int row = flow.row_from_flat_unchecked(flat), col = flow.col_from_flat_unchecked(flat); const auto& m = flow.meta();
            lsdtt::lsdtt_channel_head_metric_record h; h.head_id = head_id++; h.outlet_id = outlet_for_flat.at(basin.outlet_flat); h.unit_id = unit;
            h.head_x = m.x_min + (col + 0.5) * m.cell_size; h.head_y = m.y_max() - (row + 0.5) * m.cell_size;
            h.head_elevation_m = dem.native().get_data()[static_cast<std::size_t>(flat)]; h.chi = basin.chi[i];
            h.nearest_divide_x = h.nearest_divide_y = h.divide_elevation_m = h.distance_to_divide_m = h.relief_m = h.gradient_to_divide = std::numeric_limits<double>::quiet_NaN();
            heads.push_back(std::move(h));
        }
    }
    const std::size_t matched_heads = lsdtt::compute_nearest_divide_metrics(divides, dem.native(), heads);
    (void) matched_heads;
    lsdtt::lsdtt_pair_summary_options summary_options; summary_options.uncertainty_method = lsdtt::parse_pair_uncertainty_method(uncertainty);
    summary_options.bootstrap_replicates = bootstrap_replicates; summary_options.bootstrap_seed = seed;
    auto summaries = lsdtt::compute_pair_summaries(divides, heads, summary_options);
    const auto labelled = manager.get_junction_basins_raster();
    std::map<int, lsdtt::lsdtt_unit_record> by_unit;
    for (const auto flat : basin_outlets) { const int unit = unit_for_flat.at(flat); auto& u = by_unit[unit]; u.unit_id = unit; u.unit_key = ""; u.snapped_outlet_flats.push_back(flat); u.point_ids.push_back(outlet_for_flat.at(flat)); }
    for (const auto& r : outlet_records) if (r.snapped && r.unit_id >= 0) by_unit[r.unit_id].unit_key = r.unit_key;
    std::vector<lsdtt::lsdtt_unit_record> units; for (auto& [_, unit] : by_unit) units.push_back(std::move(unit));
    return DivideMigrationResult(flow_owner.shared_native(), std::move(network), std::move(chi), Raster(lsdtt::lsdtt_raster_cpu(std::vector<float>(labelled.data), std::vector<std::uint8_t>(labelled.mask), labelled.meta)), std::move(outlet_records), std::move(units), std::move(divides), std::move(heads), std::move(summaries));
}

} // namespace

void bind_divide_migration(nb::module_& module)
{
    nb::class_<DivideMigrationResult>(module, "DivideMigrationResult")
        .def_prop_ro("network", &DivideMigrationResult::network, nb::rv_policy::reference_internal)
        .def_prop_ro("chi", &DivideMigrationResult::chi, nb::rv_policy::reference_internal)
        .def_prop_ro("basin_labels", &DivideMigrationResult::basin_labels, nb::rv_policy::reference_internal)
        .def_prop_ro("outlets", &outlet_table).def_prop_ro("units", &unit_table)
        .def_prop_ro("divides", &divide_gdf).def_prop_ro("channel_heads", &head_gdf)
        .def_prop_ro("pair_summary", &summary_table);
    module.def("divide_migration", &divide_migration, nb::arg("flow"), nb::arg("dem"), nb::arg("outlets"),
        nb::arg("clusters") = nb::none(), nb::arg("source_pixels") = 1000,
        nb::arg("snapping") = "downstream", nb::arg("outlet_control") = "complete_only",
        nb::arg("min_outlet_elevation") = -std::numeric_limits<double>::infinity(),
        nb::arg("m_over_n") = 0.45, nb::arg("reference_area") = 1.0,
        nb::arg("uncertainty") = "std_err", nb::arg("bootstrap_replicates") = 1000, nb::arg("seed") = 5489u,
        "Compare channel-head elevation, chi, relief, and gradient across adjacent outlet-defined basin divides.");
}
} // namespace lsdtt3_python
