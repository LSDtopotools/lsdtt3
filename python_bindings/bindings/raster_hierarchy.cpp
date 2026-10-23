#include <stdexcept>
#include <string>

#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>

#include <cpu/lsdtt_raster_hierarchy_cpu.hpp>

#include "binding_common.hpp"

namespace lsdtt3_python {
namespace {

lsdtt::raster_hierarchy_storage parse_storage(const std::string& storage)
{
    if (storage == "active_only") return lsdtt::raster_hierarchy_storage::active_only;
    if (storage == "retain_allocated") return lsdtt::raster_hierarchy_storage::retain_allocated;
    if (storage == "preallocated") return lsdtt::raster_hierarchy_storage::preallocated;
    throw std::invalid_argument(
        "storage must be active_only, retain_allocated, or preallocated");
}

class RasterHierarchy
{
public:
    RasterHierarchy(const Raster& finest, const std::string& storage,
                    bool keep_finest, std::size_t levels)
        : hierarchy_(finest.native(), parse_storage(storage), keep_finest, levels) {}

    Raster level(std::size_t index) const { return Raster(hierarchy_.level(index)); }
    Raster restrict_to(std::size_t target, const std::string& method)
    {
        return Raster(hierarchy_.restrict_to(target, raster_resampling_method(method)));
    }
    Raster prolong_to(std::size_t target, const std::string& method)
    {
        return Raster(hierarchy_.prolong_to(target, raster_resampling_method(method)));
    }
    Raster move_to(std::size_t target, const std::string& method)
    {
        return Raster(hierarchy_.move_to(target, raster_resampling_method(method)));
    }
    std::size_t add_level_for_factor(std::size_t factor)
    {
        return hierarchy_.add_level_for_factor(factor);
    }
    void set_active(const Raster& raster) { hierarchy_.set_active(raster.native()); }

    lsdtt::lsdtt_raster_hierarchy_cpu& native() noexcept { return hierarchy_; }
    const lsdtt::lsdtt_raster_hierarchy_cpu& native() const noexcept { return hierarchy_; }

private:
    lsdtt::lsdtt_raster_hierarchy_cpu hierarchy_;
};

} // namespace

void bind_raster_hierarchy(nb::module_& m)
{
    nb::class_<RasterHierarchy>(m, "RasterHierarchy", R"doc(
Resolution hierarchy for one in-memory float Raster.

Level 0 is finest. Standard levels are exact factor-two grids; explicit
integer factors can be added over the same extent.
The hierarchy stores no solver state and each call returning a Raster returns
a copy, consistent with the rest of the Python Raster API.

Parameters
----------
finest : Raster
    Finest-resolution raster (level 0).
storage : {"active_only", "retain_allocated", "preallocated"}, default="retain_allocated"
    Which level rasters are kept in memory: only the active level, every
    level once allocated, or every level from the start.
keep_finest : bool, default=True
    Always keep level 0 in memory.
levels : int, default=0
    Number of standard factor-two levels to lay out; 0 uses every exact
    factor-two level the grid allows.
)doc")
        .def(nb::init<const Raster&, const std::string&, bool, std::size_t>(),
             nb::arg("finest"), nb::arg("storage") = "retain_allocated",
             nb::arg("keep_finest") = true, nb::arg("levels") = 0)
        .def_prop_ro("n_levels", [](const RasterHierarchy& self) {
            return self.native().n_levels();
        })
        .def_prop_ro("active_level", [](const RasterHierarchy& self) {
            return self.native().active_level();
        })
        .def("level_factor", [](const RasterHierarchy& self, std::size_t level) {
            const auto& finest = self.native().level_meta(0);
            const auto& current = self.native().level_meta(level);
            return current.n_rows == 0 ? std::size_t{0} :
                static_cast<std::size_t>(finest.n_rows / current.n_rows);
        }, nb::arg("level"), "Integer coarsening factor relative to level zero.")
        .def("level_shape", [](const RasterHierarchy& self, std::size_t level) {
            const auto& meta = self.native().level_meta(level);
            return nb::make_tuple(meta.n_rows, meta.n_cols);
        }, nb::arg("level"), "Grid shape as (n_rows, n_cols).")
        .def("has_level", [](const RasterHierarchy& self, std::size_t level) {
            return self.native().has_level(level);
        }, nb::arg("level"))
        .def("allocate_to", [](RasterHierarchy& self, std::size_t level) {
            self.native().allocate_to(level);
        }, nb::arg("level"))
        .def("add_level_for_factor", &RasterHierarchy::add_level_for_factor,
             nb::arg("factor"))
        .def("level", &RasterHierarchy::level, nb::arg("level"))
        .def("restrict_to", &RasterHierarchy::restrict_to, nb::arg("level"),
             nb::arg("method") = "conservative_mean")
        .def("prolong_to", &RasterHierarchy::prolong_to, nb::arg("level"),
             nb::arg("method") = "linear")
        .def("move_to", &RasterHierarchy::move_to, nb::arg("level"),
             nb::arg("method"))
        .def("set_active", &RasterHierarchy::set_active, nb::arg("raster"))
        .def("clear_except", [](RasterHierarchy& self, std::size_t level) {
            self.native().clear_except(level);
        }, nb::arg("level"))
        .def("clear", [](RasterHierarchy& self) { self.native().clear(); });
}

} // namespace lsdtt3_python
