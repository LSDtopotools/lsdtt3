#pragma once

// lsdtt3 — lsdtt_junction_basin_selection_cpu.hpp
// Manages a collection of lsdtt_junction_basin_cpu objects: selection,
// filtering, and per-basin value extraction for use with src/utils/
// statistical functions.
//
// Selection
// ---------
// Basins can be seeded either from explicit outlet flat indices or from
// the full junction network (one basin per junction). Both paths produce
// lsdtt_junction_basin_cpu objects.
//
// Filtering
// ---------
// Filters are member functions that return *this, enabling optional
// method chaining:
//
//   manager
//       .filter_incomplete(boundary_conditions)
//       .filter_by_area(min_pixels, max_pixels)
//       .filter_by_stream_order(min_order, max_order)
//       .filter_nested()
//       .filter_duplicates();
//
// Per-basin extraction
// --------------------
// apply_per_basin() accepts any callable that maps a const
// lsdtt_junction_basin_cpu& to an arbitrary return type R. It returns
// std::vector<R>, one element per surviving basin. The return type is
// fully deduced from the extractor, so it works equally for scalars,
// structs (e.g. lsdtt_descriptive_stats), or any other type returned
// by src/utils/ routines:
//
//   auto stats = manager.apply_per_basin([&](const lsdtt_junction_basin_cpu& b)
//   {
//       auto values = b.get_basin_values(dem);
//       return lsdtt::compute_descriptive_stats(values);
//   });
//   // stats is std::vector<lsdtt_descriptive_stats>
//
// apply_per_basin is defined in this header (template body must be
// visible at the call site).
//
// Design notes
// ------------
// - lsdtt_junction_basin_selection_cpu is a class (not a struct) because
//   it owns private state (the basin vector) that must be protected from
//   direct mutation to keep the object consistent after filtering.
// - Statistics functions are not member functions of this class. Instead,
//   use apply_per_basin with a lambda, or call the named free functions in
//   lsdtt_basin_analysis_cpu.hpp, which compose this class with src/utils/.
// - Both flowinfo and junction_network must outlive this manager (and all
//   basins it holds), because lsdtt_junction_basin_cpu stores non-owning
//   pointers to them.
//
// Copyright Simon M. Mudd, GPL-3.0

#include "lsdtt_junction_basin_cpu.hpp"
#include "lsdtt_flowinfo_cpu.hpp"
#include "lsdtt_junction_network_cpu.hpp"

#include <cstdint>
#include <string>
#include <type_traits>
#include <vector>

namespace lsdtt
{

// -----------------------------------------------------------------------
// lsdtt_junction_basin_selection_cpu
// -----------------------------------------------------------------------

class lsdtt_junction_basin_selection_cpu
{
public:

    // -----------------------------------------------------------------------
    // Construction
    // -----------------------------------------------------------------------

    // Constructs an empty manager.
    // Both flowinfo and junction_network must outlive this object.
    lsdtt_junction_basin_selection_cpu(
        const lsdtt_flowinfo_cpu&         flowinfo,
        const lsdtt_junction_network_cpu& junction_network);

    // -----------------------------------------------------------------------
    // Basin selection
    // -----------------------------------------------------------------------

    // Seeds one basin per outlet flat index in the supplied vector.
    // For each outlet, the constructor of lsdtt_junction_basin_cpu walks
    // upstream to find the associated junction. Any existing basins are
    // replaced.
    void select_basins_from_outlets(const std::vector<std::int64_t>& outlet_flat_indices);

    // Seeds one basin per junction in the junction network.
    // Any existing basins are replaced.
    void select_all_junction_basins();

    // Seeds one basin per junction encountered when walking upstream from
    // each supplied outlet flat index. The downstream-most basin in each
    // chain uses the supplied flat index as its outlet, not the junction
    // pixel. Any existing basins are replaced.
    void select_junction_basins_from_outlets(
        const std::vector<std::int64_t>& outlet_flat_indices);

    // -----------------------------------------------------------------------
    // Filtering (chainable)
    // -----------------------------------------------------------------------

    // Removes basins that are hydrologically incomplete: any basin where
    // check_basin_completeness returns false is deleted.
    // boundary_conditions is a four-character string (north/east/south/west),
    // e.g. "oooo". See lsdtt_types.hpp for boundary code definitions.
    lsdtt_junction_basin_selection_cpu&
        filter_incomplete(const std::string& boundary_conditions);

    // Removes basins whose pixel count is outside [min_pixels, max_pixels].
    // Pass -1 for either bound to leave that end open.
    lsdtt_junction_basin_selection_cpu&
        filter_by_area(std::int32_t min_pixels, std::int32_t max_pixels);

    // Removes basins whose stream order is outside [min_order, max_order].
    // Pass -1 for either bound to leave that end open.
    lsdtt_junction_basin_selection_cpu&
        filter_by_stream_order(std::int32_t min_order, std::int32_t max_order);

    // Removes all basins that are nested inside another surviving basin,
    // i.e. whose outlet flat index is upstream of another basin's outlet.
    lsdtt_junction_basin_selection_cpu& filter_nested();

    // Removes duplicate basins: if two basins share the same outlet flat
    // index, only the first one encountered is kept.
    lsdtt_junction_basin_selection_cpu& filter_duplicates();

    // Retains only the single basin with the largest pixel count.
    // If the collection is empty this is a no-op.
    lsdtt_junction_basin_selection_cpu& filter_largest();

    // -----------------------------------------------------------------------
    // Per-basin extraction (template — must be defined in header)
    // -----------------------------------------------------------------------

    // Applies extractor to every surviving basin and collects the results
    // into a std::vector<R>, where R is the return type of the extractor.
    // The extractor must be callable as:
    //   R extractor(const lsdtt_junction_basin_cpu&)
    //
    // Example — descriptive statistics per basin:
    //   auto stats = manager.apply_per_basin([&](const lsdtt_junction_basin_cpu& b)
    //   {
    //       auto values = b.get_basin_values(dem);
    //       return lsdtt::compute_descriptive_stats(values);
    //   });
    //
    // Example — scalar area per basin:
    //   auto areas = manager.apply_per_basin([&](const lsdtt_junction_basin_cpu& b)
    //   {
    //       return static_cast<double>(b.get_number_of_basin_pixels()) * pixel_area;
    //   });
    template <typename Extractor>
    auto apply_per_basin(Extractor&& extractor) const
        -> std::vector<std::invoke_result_t<Extractor, const lsdtt_junction_basin_cpu&>>
    {
        using result_type = std::invoke_result_t<Extractor, const lsdtt_junction_basin_cpu&>;

        std::vector<result_type> results;
        results.reserve(basins_.size());

        for (const lsdtt_junction_basin_cpu& basin : basins_)
        {
            results.push_back(extractor(basin));
        }

        return results;
    }

    // -----------------------------------------------------------------------
    // Accessors
    // -----------------------------------------------------------------------

    // Returns the number of basins currently in the selection.
    [[nodiscard]] std::size_t size() const noexcept;

    // Returns true if no basins survive.
    [[nodiscard]] bool empty() const noexcept;

    // Read-only access to the basin at position i.
    [[nodiscard]] const lsdtt_junction_basin_cpu& basin_at(std::size_t i) const;

    // Returns the outlet flat index for every surviving basin, in order.
    // Convenient for writing results to file alongside per-basin statistics.
    [[nodiscard]] std::vector<std::int64_t> get_outlet_flat_indices() const;

    // Returns the junction number for every surviving basin, in order.
    [[nodiscard]] std::vector<std::int32_t> get_junction_numbers() const;

    // Returns the stream order for every surviving basin, in order.
    [[nodiscard]] std::vector<std::int32_t> get_stream_orders() const;

private:

    const lsdtt_flowinfo_cpu*         flowinfo_;          // non-owning
    const lsdtt_junction_network_cpu* junction_network_;  // non-owning
    std::vector<lsdtt_junction_basin_cpu> basins_;
};

} // namespace lsdtt