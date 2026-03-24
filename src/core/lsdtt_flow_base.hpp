// lsdtt_flow_base.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2025 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html
//
// -----------------------------------------------------------------------------
// INFRASTRUCTURE FILE — you probably do not need to include this directly.
//
// If you are writing analysis code that uses flow relationships, include
// lsdtt_flowinfo_cpu.hpp instead and work with lsdtt_flowinfo directly.
//
// If you are writing a CLI that builds a flow routing object, include this
// file and call make_flow_base(). The factory returns the correct backend
// (currently CPU; GPU when implemented) without the CLI needing to know
// which backend is in use.
//
//   #include "../core/lsdtt_flow_base.hpp"
//
//   auto flow = lsdtt::make_flow_base(dem_data, dem_mask, meta, options);
//   // flow is unique_ptr<lsdtt_flow_base>
//   flow->upslope_nodes(flat);
//   flow->n_valid_cells();
//
// The three flow infrastructure files and what each does:
//
//   lsdtt_flow_types.hpp  — shared types (std::int64_t, sentinel values,
//                           lsdtt_flowinfo_options, lsdtt_flat_index_slice).
//                           Included by both this file and lsdtt_flowinfo_cpu.hpp
//                           to avoid a circular include dependency.
//
//   lsdtt_flow_base.hpp   — THIS FILE. Abstract interface (pure virtual base
//                           class) that lsdtt_flowinfo_cpu implements.
//                           Exposes the same query API as lsdtt_flowinfo.
//
//   lsdtt_flow_base.cpp   — Lives in src/cpu/ (not src/core/) because it
//   (in src/cpu/)           instantiates lsdtt_flowinfo_cpu. Contains only
//                           make_flow_base(). If you add a GPU backend, add
//                           the gpu_available() branch there.
//
// -----------------------------------------------------------------------------
// lsdtt_flow_base
//
// Abstract base class for D8 flow routing backends. Mirrors lsdtt_raster_base
// for the flow routing layer. lsdtt_flowinfo_cpu is the CPU implementation;
// a future lsdtt_flowinfo_gpu will be the GPU implementation.
//
// std::int64_t and all sentinel values come from lsdtt_flow_types.hpp —
// the single source of truth for flow routing types.
//
// Thread safety:
//   All const methods are safe to call concurrently after construction.
// -----------------------------------------------------------------------------

#pragma once

#include <functional>
#include <memory>
#include <string>
#include <vector>

#include "lsdtt_raster_meta.hpp"
#include "lsdtt_flow_types.hpp"
#include "../utils/lsdtt_array2d.hpp"

namespace lsdtt
{

class lsdtt_flow_base
{
public:

    // -------------------------------------------------------------------------
    // Destructor
    // -------------------------------------------------------------------------
    virtual ~lsdtt_flow_base() = default;

    // -------------------------------------------------------------------------
    // Object state
    // -------------------------------------------------------------------------
    [[nodiscard]] virtual bool        empty()                const noexcept = 0;
    [[nodiscard]] virtual bool        is_acyclic()           const noexcept = 0;
    [[nodiscard]] virtual bool        has_s_vector()         const noexcept = 0;
    [[nodiscard]] virtual bool        has_topological_order() const noexcept = 0;
    [[nodiscard]] virtual std::int64_t n_valid_cells()       const noexcept = 0;
    [[nodiscard]] virtual std::size_t  n_outlets()           const noexcept = 0;
    [[nodiscard]] virtual std::size_t  full_grid_size()      const noexcept = 0;

    // -------------------------------------------------------------------------
    // Metadata
    // -------------------------------------------------------------------------
    [[nodiscard]] virtual const lsdtt_raster_meta& meta()      const noexcept = 0;
    [[nodiscard]] virtual int                       n_rows()    const noexcept = 0;
    [[nodiscard]] virtual int                       n_cols()    const noexcept = 0;
    [[nodiscard]] virtual double                    cell_size() const noexcept = 0;

    // -------------------------------------------------------------------------
    // Index conversion
    // -------------------------------------------------------------------------
    [[nodiscard]] virtual bool in_bounds(int row, int col)        const noexcept = 0;
    [[nodiscard]] virtual bool is_flat_in_range(std::int64_t)  const noexcept = 0;

    [[nodiscard]] virtual std::int64_t flat_index_unchecked(
        int row, int col) const noexcept = 0;

    [[nodiscard]] virtual int row_from_flat_unchecked(
        std::int64_t flat) const noexcept = 0;

    [[nodiscard]] virtual int col_from_flat_unchecked(
        std::int64_t flat) const noexcept = 0;

    // Checked variants — throw on out-of-bounds
    [[nodiscard]] virtual std::int64_t flat_index(int row, int col) const = 0;
    [[nodiscard]] virtual int             row_from_flat(std::int64_t)     const = 0;
    [[nodiscard]] virtual int             col_from_flat(std::int64_t)     const = 0;

    // -------------------------------------------------------------------------
    // Cell validity and routing queries
    // -------------------------------------------------------------------------
    [[nodiscard]] virtual bool is_valid_cell(std::int64_t flat) const noexcept = 0;
    [[nodiscard]] virtual bool is_valid_cell(int row, int col)     const noexcept = 0;

    [[nodiscard]] virtual bool            is_outlet_unchecked(std::int64_t)         const noexcept = 0;
    [[nodiscard]] virtual std::int64_t receiver_of_unchecked(std::int64_t)       const noexcept = 0;
    [[nodiscard]] virtual std::uint8_t    receiver_dir_of_unchecked(std::int64_t)   const noexcept = 0;
    [[nodiscard]] virtual std::uint8_t donor_count_of_unchecked(std::int64_t)   const noexcept = 0;

    // Checked variants — throw if cell is invalid or object not built
    [[nodiscard]] virtual bool             is_outlet(std::int64_t)       const = 0;
    [[nodiscard]] virtual std::int64_t  receiver_of(std::int64_t)     const = 0;
    [[nodiscard]] virtual std::uint8_t     receiver_dir_of(std::int64_t) const = 0;
    [[nodiscard]] virtual std::uint8_t donor_count_of(std::int64_t)  const = 0;

    // -------------------------------------------------------------------------
    // s_vector queries
    //
    // upslope_nodes() returns a zero-copy slice valid for the lifetime of
    // this engine object. The GPU backend may materialise the slice on demand
    // (copying from device memory) — callers should not retain slices across
    // calls or across engine lifetime.
    // -------------------------------------------------------------------------
    [[nodiscard]] virtual std::int64_t    s_vector_start_unchecked(std::int64_t) const noexcept = 0;
    [[nodiscard]] virtual std::uint64_t s_vector_size_unchecked(std::int64_t)  const noexcept = 0;

    [[nodiscard]] virtual std::int64_t    s_vector_start(std::int64_t) const = 0;
    [[nodiscard]] virtual std::uint64_t s_vector_size(std::int64_t)  const = 0;

    [[nodiscard]] virtual lsdtt_flat_index_slice upslope_nodes(std::int64_t) const = 0;
    [[nodiscard]] virtual lsdtt_flat_index_slice upslope_nodes(int row, int col) const = 0;
    [[nodiscard]] virtual lsdtt_flat_index_slice upslope_nodes_unchecked(std::int64_t) const noexcept = 0;

    [[nodiscard]] virtual bool is_upslope_of(
        std::int64_t candidate,
        std::int64_t target) const = 0;

    [[nodiscard]] virtual bool is_upslope_of_unchecked(
        std::int64_t candidate,
        std::int64_t target) const noexcept = 0;

    // -------------------------------------------------------------------------
    // Iteration helpers
    //
    // for_each_upslope_node visits every node in the upslope set of flat,
    // calling function(node_flat) for each. The order matches the s_vector
    // preorder traversal (node before its donors).
    // -------------------------------------------------------------------------
    virtual void for_each_upslope_node(
        std::int64_t flat,
        const std::function<void(std::int64_t)>& function) const = 0;

    virtual void for_each_upslope_node_unchecked(
        std::int64_t flat,
        const std::function<void(std::int64_t)>& function) const noexcept = 0;

    // channel extraction
    [[nodiscard]] virtual std::vector<std::int64_t>
        get_sources_index_threshold(int threshold) const = 0;


    // -------------------------------------------------------------------------
    // Sentinel values — defined in lsdtt_flow_types.hpp
    // -------------------------------------------------------------------------
    static constexpr std::int64_t nodata_receiver = lsdtt::k_nodata_receiver;
    static constexpr std::int64_t outlet_receiver  = lsdtt::k_outlet_receiver;
    static constexpr std::uint8_t    no_direction     = lsdtt::k_no_direction;
};

// -----------------------------------------------------------------------------
// make_flow_base
//
// Factory function. Returns a CPU lsdtt_flowinfo engine unless a GPU backend
// is available and enabled, in which case it returns a GPU engine.
//
// Parameters match lsdtt_flowinfo_options. When the GPU backend is
// implemented, additional GPU-specific options can be added to
// lsdtt_flowinfo_options or a derived options struct without changing this
// signature.
// -----------------------------------------------------------------------------
std::unique_ptr<lsdtt_flow_base> make_flow_base(
    const lsdtt_array2d_f&        dem,
    const lsdtt_array2d_u8&       mask,
    const lsdtt_raster_meta&      meta,
    const lsdtt_flowinfo_options& options = lsdtt_flowinfo_options{}
);

} // namespace lsdtt
