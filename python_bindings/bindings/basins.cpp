#include <algorithm>
#include <cstdint>
#include <memory>
#include <limits>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <utility>
#include <vector>

#include <nanobind/nanobind.h>
#include <nanobind/stl/vector.h>

#include <cpu/lsdtt_basin_cpu.hpp>

#include "binding_common.hpp"

namespace lsdtt3_python {
namespace {

std::vector<std::int64_t> outlet_flats(const FlowInfo& flow, nb::object outlets,
                                       std::int64_t min_upslope_pixels)
{
    if (min_upslope_pixels < 0)
        throw std::invalid_argument("min_upslope_pixels must be non-negative");
    if (outlets.is_none())
        return lsdtt::get_all_basin_outlet_flats(flow.native(), min_upslope_pixels);

    nb::object array = nb::module_::import_("numpy").attr("asarray")(outlets, nb::arg("dtype") = "int64");
    if (nb::cast<int>(array.attr("ndim")) != 2 || nb::cast<int>(array.attr("shape")[1]) != 2)
        throw std::invalid_argument("outlets must be an N by 2 row, col array");
    const auto pairs = nb::cast<std::vector<std::vector<std::int64_t>>>(array.attr("tolist")());
    std::unordered_set<std::int64_t> seen;
    std::vector<std::int64_t> result;
    for (const auto& pair : pairs) {
        const auto row = pair[0]; const auto col = pair[1];
        if (row < 0 || row >= flow.native().n_rows() || col < 0 || col >= flow.native().n_cols())
            throw std::invalid_argument("outlets contains an out-of-grid pixel");
        const auto flat = row * static_cast<std::int64_t>(flow.native().n_cols()) + col;
        if (!flow.native().is_valid_cell(flat)) throw std::invalid_argument("outlets selects a nodata pixel");
        if (seen.insert(flat).second) result.push_back(flat);
    }
    return result;
}

class BasinCollection {
public:
    BasinCollection(const BasinCollection&) = delete;
    BasinCollection& operator=(const BasinCollection&) = delete;
    BasinCollection(BasinCollection&&) = default;
    BasinCollection& operator=(BasinCollection&&) = default;
    BasinCollection(std::shared_ptr<lsdtt::lsdtt_flowinfo> flow, std::vector<std::int64_t> outlets)
        : flow_(std::move(flow))
    {
        basins_.reserve(outlets.size());
        for (const auto outlet : outlets)
            basins_.push_back(std::make_unique<lsdtt::lsdtt_basin_cpu>(outlet, *flow_));
        parent_.assign(basins_.size(), -1);
        for (std::size_t child = 0; child < basins_.size(); ++child) {
            const auto child_outlet = basins_[child]->get_basin_outlet_flat_index();
            std::int32_t best_parent = -1;
            std::int32_t best_size = std::numeric_limits<std::int32_t>::max();
            for (std::size_t candidate = 0; candidate < basins_.size(); ++candidate) {
                if (candidate == child) continue;
                const auto parent_outlet = basins_[candidate]->get_basin_outlet_flat_index();
                if (!flow_->is_upslope_of(child_outlet, parent_outlet)) continue;
                const auto n = basins_[candidate]->pixel_count();
                if (n < best_size) { best_size = n; best_parent = static_cast<std::int32_t>(candidate); }
            }
            parent_[child] = best_parent;
        }
        depth_.assign(basins_.size(), 0);
        for (std::size_t i = 0; i < basins_.size(); ++i) {
            for (auto parent = parent_[i]; parent >= 0; parent = parent_[static_cast<std::size_t>(parent)])
                ++depth_[i];
        }
    }

    std::size_t size() const noexcept { return basins_.size(); }
    const lsdtt::lsdtt_flowinfo& flow() const noexcept { return *flow_; }
    const lsdtt::lsdtt_basin_cpu& basin(std::size_t id) const {
        if (id >= basins_.size()) throw std::out_of_range("basin_id is out of range");
        return *basins_[id];
    }
    std::int32_t parent(std::size_t id) const { if (id >= parent_.size()) throw std::out_of_range("basin_id is out of range"); return parent_[id]; }
    std::int32_t depth(std::size_t id) const { if (id >= depth_.size()) throw std::out_of_range("basin_id is out of range"); return depth_[id]; }
    std::vector<std::int32_t> memberships(std::int64_t flat) const {
        std::vector<std::int32_t> ids;
        for (std::size_t i = 0; i < basins_.size(); ++i)
            if (flow_->is_upslope_of(flat, basins_[i]->get_basin_outlet_flat_index())) ids.push_back(static_cast<std::int32_t>(i));
        return ids;
    }
private:
    std::shared_ptr<lsdtt::lsdtt_flowinfo> flow_;
    std::vector<std::unique_ptr<lsdtt::lsdtt_basin_cpu>> basins_;
    std::vector<std::int32_t> parent_, depth_;
};

Raster basin_mask(const BasinCollection& collection, std::size_t id)
{
    const auto& flow = collection.flow();
    const std::size_t n = flow.full_grid_size();
    std::vector<float> data(n, 0.0f);
    std::vector<std::uint8_t> mask(flow.mask_array().data(), flow.mask_array().data() + n);
    for (const auto flat : collection.basin(id).get_slice()) data[static_cast<std::size_t>(flat)] = 1.0f;
    return Raster(lsdtt::lsdtt_raster_cpu(std::move(data), std::move(mask), flow.meta()));
}

Raster basin_labels(const BasinCollection& collection)
{
    const auto& flow = collection.flow(); const std::size_t n = flow.full_grid_size();
    std::vector<float> data(n, 0.0f);
    std::vector<std::uint8_t> mask(flow.mask_array().data(), flow.mask_array().data() + n);
    for (std::size_t id = 0; id < collection.size(); ++id)
        for (const auto flat : collection.basin(id).get_slice()) data[static_cast<std::size_t>(flat)] = static_cast<float>(id + 1);
    return Raster(lsdtt::lsdtt_raster_cpu(std::move(data), std::move(mask), flow.meta()));
}

nb::object basin_pixels(const BasinCollection& collection, std::size_t id)
{
    const auto& flow = collection.flow(); const auto& slice = collection.basin(id).get_slice();
    std::vector<std::int64_t> values; values.reserve(slice.size() * 2);
    for (const auto flat : slice) { values.push_back(flow.row_from_flat_unchecked(flat)); values.push_back(flow.col_from_flat_unchecked(flat)); }
    return nb::module_::import_("numpy").attr("array")(values, nb::arg("dtype") = "int64")
        .attr("reshape")(nb::make_tuple(static_cast<int>(slice.size()), 2));
}

nb::object basin_table(const BasinCollection& collection)
{
    const auto& flow = collection.flow(); nb::list rows;
    for (std::size_t id = 0; id < collection.size(); ++id) {
        const auto& basin = collection.basin(id); const auto outlet = basin.get_basin_outlet_flat_index();
        nb::dict row; row["basin_id"] = id; row["outlet_flat"] = outlet;
        row["outlet_row"] = flow.row_from_flat_unchecked(outlet); row["outlet_col"] = flow.col_from_flat_unchecked(outlet);
        row["n_pixels"] = basin.pixel_count(); row["area_metres2"] = basin.pixel_count() * flow.meta().cell_size * flow.meta().cell_size;
        row["complete"] = basin.check_basin_completeness(); row["outlet_on_dem_edge"] = lsdtt::is_flat_on_dem_edge(flow, outlet);
        row["parent_id"] = collection.parent(id); row["nesting_depth"] = collection.depth(id); rows.append(std::move(row));
    }
    return nb::module_::import_("pandas").attr("DataFrame")(rows);
}

nb::object basin_membership(const BasinCollection& collection, int row, int col)
{
    const auto& flow = collection.flow(); const auto flat = flow.flat_index(row, col);
    if (!flow.is_valid_cell(flat)) throw std::invalid_argument("row, col selects nodata");
    nb::object numpy = nb::module_::import_("numpy");
    return numpy.attr("array")(collection.memberships(flat), nb::arg("dtype") = "int32");
}

BasinCollection make_basins(const FlowInfo& flow, nb::object outlets, std::int64_t min_upslope_pixels)
{
    return BasinCollection(flow.shared_native(), outlet_flats(flow, outlets, min_upslope_pixels));
}

} // namespace

void bind_basins(nb::module_& module, nb::class_<FlowInfo>& flow)
{
    nb::class_<BasinCollection>(module, "BasinCollection", R"doc(
Selected nested D8 drainage basins retaining their parent ``FlowInfo`` topology.

Each basin is defined by an outlet pixel and stores an efficient native
upslope slice; no raster-sized copy is retained per basin.
)doc")
        .def_prop_ro("size", &BasinCollection::size)
        .def_prop_ro("labels", &basin_labels, "Basin IDs (one-based); nested basins overwrite parents.")
        .def_prop_ro("table", &basin_table)
        .def("mask", &basin_mask, nb::arg("basin_id"))
        .def("pixels", &basin_pixels, nb::arg("basin_id"))
        .def("membership", &basin_membership, nb::arg("row"), nb::arg("col"), "Return all selected basin IDs containing this pixel.");
    flow.def("basins", &make_basins, nb::arg("outlets") = nb::none(), nb::arg("min_upslope_pixels") = 0,
        R"doc(
Create a ``BasinCollection`` from explicit outlet pixels or all qualifying D8 outlets.

Explicit outlets may be internal channel pixels, which creates nested basins.
)doc");
}
} // namespace lsdtt3_python
