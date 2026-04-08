#pragma once

// lsdtt3 — lsdtt_basin_cpu.hpp
// CPU implementation of a single drainage basin.
// Holds a zero-copy slice into the flow-routing s_vector
// (lsdtt_flat_index_slice) — no pixel data is copied on construction.
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

    [[nodiscard]] bool
        check_basin_completeness(const std::string& boundary_conditions) const override;

    [[nodiscard]] std::int32_t pixel_count() const noexcept override;
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

    // Returns the outline of the basin as a sequence of flat indices tracing
    // the basin boundary. The lsdtt_basin_outline type is defined below.
    // TODO: lsdtt_basin_outline is not yet designed — stub returns empty vector.
    [[nodiscard]] std::vector<std::int64_t> get_basin_outline() const;

private:

    static void validate_inputs(const lsdtt_flow_base& flow,
                                std::int64_t outlet_flat_index);

    [[nodiscard]] static lsdtt_flat_index_slice make_slice(
        const lsdtt_flow_base& flow,
        std::int64_t basin_outlet_flat_index);

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
