#pragma once

// lsdtt3 - lsdtt_basin_cpu.hpp
// CPU implementation of a single drainage basin.
// Holds a zero-copy slice into the flow-routing s_vector
// (lsdtt_flat_index_slice) - no pixel data is copied on construction.
//
// Member functions are limited to operations that require access to
// private data (the slice and outlet index). Statistical analysis on
// extracted value vectors belongs in free functions in src/cpu/,
// following the design rule: "write free functions unless you need
// private data members".
//
// Copyright Simon M. Mudd, GPL-3.0

#include "../core/lsdtt_basin_base.hpp"
#include "../core/lsdtt_flow_base.hpp"
#include "../core/lsdtt_types.hpp"
#include "lsdtt_raster_cpu.hpp"
#include "lsdtt_basin_outline_cpu.hpp"

#include <cstdint>
#include <string>
#include <vector>

namespace lsdtt
{

// -----------------------------------------------------------------------
// lsdtt_basin_cpu
// -----------------------------------------------------------------------

// Represents a single drainage basin on the CPU backend.
//
// Construction requires:
//   - the flat index of the outlet node
//   - a const reference to the lsdtt_flow_base that owns the s_vector
//     (the flow object must outlive this basin)
//
// The basin stores an lsdtt_flat_index_slice, a zero-copy view into the
// flow s_vector covering all upstream pixels. No pixel data is copied.
class lsdtt_basin_cpu : public lsdtt_basin_base
{
public:

    // Constructs a basin from an outlet flat index and flow object.
    // The flow object must outlive this basin.
    lsdtt_basin_cpu(std::int64_t basin_outlet_flat_index,
                    const lsdtt_flow_base& flow);

    // -----------------------------------------------------------------------
    // lsdtt_basin_base interface
    // -----------------------------------------------------------------------

    [[nodiscard]] std::int64_t get_basin_outlet_flat_index() const noexcept override;

    [[nodiscard]] std::int32_t get_number_of_basin_pixels() const noexcept override;

    [[nodiscard]] bool check_basin_completeness() const override;

    [[nodiscard]] std::int32_t pixel_count() const noexcept override;

    // Paint this basin into an existing raster_result and return the modified raster.
    // The caller supplies the fill value.
    [[nodiscard]] lsdtt_raster_result
        get_basin_raster(lsdtt_raster_result raster, float fill_value) const override;
    // -----------------------------------------------------------------------
    // CPU-specific accessors
    // -----------------------------------------------------------------------

    // Returns the zero-copy slice into the flow s_vector for this basin.
    [[nodiscard]] const lsdtt_flat_index_slice& get_slice() const noexcept;

    // Samples raster values for every pixel in the basin and returns them
    // in a flat vector. Nodata pixels are excluded. The order of values in
    // the returned vector follows the s_vector order (depth-first upstream).
    //
    // This is the primary bridge between the basin and src/utils/ statistics
    // functions: call get_basin_values to obtain a valid-data vector, then
    // pass that vector to any lsdtt::compute_* utility.
    [[nodiscard]] std::vector<float>
        get_basin_values(const lsdtt_raster_cpu& raster) const;

    // Returns a CPU basin outline object.
    // NOTE: lsdtt_basin_outline_cpu stores a non-owning pointer to this basin,
    // so the returned outline must not outlive the basin object.
    [[nodiscard]] lsdtt_basin_outline_cpu get_basin_outline() const;

    // CPU-only optimisation of check_basin_completeness()
    // This allows excluding basin completeness checks of part of this
    // basins svector slice, by suppling sub-slices.
    // Used by the basin manager to avoid checking nested basin pixels multiple
    // times.
    [[nodiscard]] bool check_basin_completeness_excluding_slices(
        const std::vector<lsdtt_flat_index_slice>& excluded_slices
    ) const;

    // CPU-only accessor: flow object must outlive basin (it does by construction)
    [[nodiscard]] const lsdtt_flow_base& flow() const
    {
        if (flow_ == nullptr)
        {
            throw std::runtime_error("lsdtt_basin_cpu::flow(): flow pointer is null");
        }
        return *flow_;
    }



private:

    static void validate_inputs(const lsdtt_flow_base& flow,
                                std::int64_t outlet_flat_index);

    [[nodiscard]] static lsdtt_flat_index_slice make_slice(
        const lsdtt_flow_base& flow,
        std::int64_t basin_outlet_flat_index);

    void validate_basin_raster_input(const lsdtt_raster_result& raster) const;

private:

    std::int64_t           basin_outlet_flat_index_; // note: the basin outlet
                            // flat index is NOT the same as an outlet in flow
                            // info or flow base. The basin outlet is a chosen
                            // (arbitrary) pixel denoted by a flat_index, which
                            // is the downstream most pixel in the basin. The
                            // outlet in flow info/base is just a baselevel.
    lsdtt_flat_index_slice slice_;
    const lsdtt_flow_base* flow_;   // non-owning pointer; flow outlives basin
};

} // namespace lsdtt
