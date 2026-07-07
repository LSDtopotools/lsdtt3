// lsdtt_flow_base.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
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
//   lsdtt_types.hpp       — shared types (std::int64_t, sentinel values,
//                           lsdtt_flat_index_slice).
//                           Included by both this file and lsdtt_flowinfo_cpu.hpp
//                           to avoid a circular include dependency.
//
//   lsdtt_flow_base.hpp   — THIS FILE. Abstract interface (pure virtual base
//                           class) that lsdtt_flowinfo_cpu implements.
//                           Exposes the same query API as lsdtt_flowinfo.
//                           Also has a type for lsdtt_flowinfo_options.
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
// std::int64_t and all sentinel values come from lsdtt_types.hpp —
// the single source of truth for flow routing types.
//
// Thread safety:
//   All const methods are safe to call concurrently after construction.
// -----------------------------------------------------------------------------

#pragma once

#include <cmath>
#include <cstdint>
#include <functional>
#include <memory>
#include <string>
#include <vector>

#include "lsdtt_raster_meta.hpp"
#include "lsdtt_types.hpp"
#include "lsdtt_data_loader.hpp"
#include "../utils/lsdtt_array2d.hpp"



namespace lsdtt
{

struct lsdtt_flowinfo_options
{
    std::string boundary_conditions       = "oooo";
    bool        interior_nodata_is_outlet = true;
    bool        require_acyclic_graph     = true;
    bool        build_topological_order   = true;
    bool        build_s_vector            = true;
    int         n_threads                 = 0;
    bool        allow_pits_as_outlets     = false;
};

class lsdtt_flow_base
{
public:

    // -------------------------------------------------------------------------
    // Destructor
    // -------------------------------------------------------------------------
    virtual ~lsdtt_flow_base() = default;

    // -------------------------------------------------------------------------
    // Object state
    // [[nodiscard]] is there so the compiler throws and error if this function
    //  is called but the result is discarded. Used to prevent bugs
    //  that are not visible to the compilers
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

    // This checks if cell is next to an invalid cell
    [[nodiscard]] virtual bool is_influenced_by_invalid(std::int64_t flat) const = 0;

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
    // Flow-routing raster outputs
    //
    // Two conventions are used throughout:
    //
    //   "include self"  — the pixel in question is counted as part of its own
    //                     upslope basin (s_vector_size includes the cell itself).
    //                     This is the standard hydrological drainage-area
    //                     convention: every cell contributes its own area.
    //
    //   "exclude outlet" — the pixel in question is NOT counted (s_vector_size - 1).
    //                     Useful when you want only the area/count strictly
    //                     upstream of a cell, not including the cell itself.
    //                     Headwater pixels return 0 under this convention.
    //
    // get_* methods return an lsdtt_raster_result that the caller can write to
    // disk with loader.write_full() or promote to lsdtt_raster_cpu for analysis.
    //
    // write_* methods write directly to disk; they call the corresponding get_*
    // internally and exist as a convenience so the caller does not need to
    // manage a loader object.
    // -------------------------------------------------------------------------

    // D8 flow direction — lsdtt internal direction codes:
    //   7 0 1
    //   6 x 2
    //   5 4 3
    // Outlet pixels carry their direction code (the edge they drain to).
    // Nodata pixels carry meta.no_data_value.
    [[nodiscard]] virtual lsdtt_raster_result get_d8_flow_direction_raster() const = 0;
    virtual void write_d8_flow_direction_raster(
        const std::string& out_path,
        raster_output_format format = raster_output_format::cog
    ) const = 0;

    // Number of contributing pixels — include self convention.
    // Value at each valid cell = total pixels draining through it,
    // including the cell itself (s_vector_size).
    [[nodiscard]] virtual lsdtt_raster_result get_n_contributing_pixels_raster() const = 0;
    virtual void write_n_contributing_pixels_raster(
        const std::string& out_path,
        raster_output_format format = raster_output_format::cog
    ) const = 0;

    // Number of contributing pixels — exclude self convention.
    // Value = pixels strictly upstream (s_vector_size - 1).
    // Headwater pixels return 0.
    [[nodiscard]] virtual lsdtt_raster_result get_n_contributing_pixels_exclude_outlet_raster() const = 0;
    virtual void write_n_contributing_pixels_exclude_outlet_raster(
        const std::string& out_path,
        raster_output_format format = raster_output_format::cog
    ) const = 0;

    // Drainage area — include self convention.
    // Value = n_contributing_pixels * cell_size^2  (square map units).
    [[nodiscard]] virtual lsdtt_raster_result get_drainage_area_raster() const = 0;
    virtual void write_drainage_area_raster(
        const std::string& out_path,
        raster_output_format format = raster_output_format::cog
    ) const = 0;

    // Drainage area with upstream injection from a point vector file.
    // Each successfully mapped point adds its drainage_area_field attribute
    // (in raster map units squared) to that pixel's initial weight before
    // accumulation, allowing upstream area entering from outside the tile to
    // be injected at edge-crossing river pixels.
    // Requires build_topological_order = true.
    // Default implementation throws std::runtime_error; concrete subclasses
    // override this where supported.
    [[nodiscard]] virtual lsdtt_raster_result get_drainage_area_raster_with_injection(
        const std::string& vector_path,
        const std::string& drainage_area_field = "drainage_area"
    ) const;

    // Drainage area — exclude self convention.
    // Value = (n_contributing_pixels - 1) * cell_size^2  (square map units).
    // Headwater pixels return 0.0.
    [[nodiscard]] virtual lsdtt_raster_result get_drainage_area_exclude_outlet_raster() const = 0;
    virtual void write_drainage_area_exclude_outlet_raster(
        const std::string& out_path,
        raster_output_format format = raster_output_format::cog
    ) const = 0;

    // -------------------------------------------------------------------------
    // Upslope accumulation of a secondary raster
    //
    // Both methods perform a single O(N) reverse-s_vector pass (headwaters →
    // outlets).  values and valid_mask must be flat-indexed with the same
    // length as the DEM (n_rows * n_cols).
    //
    // Pixels where valid_mask[i] == 0 contribute 0 to the running sum but
    // still propagate accumulated values received from further upstream.
    //
    // include_self (default true):
    //   true  — the pixel's own value is included (standard drainage-area
    //           convention: every cell counts itself).
    //   false — strictly upstream values only; headwaters return 0.
    //
    // get_upslope_sum_raster:
    //   Output pixel = sum of values[] over its entire upslope set.
    //
    // get_upslope_mean_raster:
    //   Output pixel = sum / count of valid upslope pixels.
    //   Nodata where no valid input pixels exist in the upslope set.
    // -------------------------------------------------------------------------
    // -------------------------------------------------------------------------
    // Flow-path distance to outlet
    //
    // Each valid pixel holds the cumulative D8 flow-path length (metres) from
    // that pixel to its basin outlet.  Outlet pixels themselves return 0.
    // Cardinal steps = cell_size; diagonal steps = cell_size * sqrt(2).
    // -------------------------------------------------------------------------
    [[nodiscard]] virtual lsdtt_raster_result get_distance_from_outlet_raster() const = 0;

    // -------------------------------------------------------------------------
    // Channel-relative height and distance
    //
    // A pixel is classified as a channel when its contributing pixel count
    // (include-self convention) >= threshold_contributing_pixels.
    //
    // get_height_above_channel_raster:
    //   For each pixel: elev[pixel] - elev[nearest_downstream_channel_pixel].
    //   Channel pixels return 0.  Pixels that drain to the DEM edge without
    //   passing through a channel carry nodata.
    //   elev / elev_mask must have flat-index size == n_rows * n_cols.
    //
    // get_distance_to_channel_raster:
    //   Cumulative D8 flow-path length (metres) to the nearest downstream
    //   channel pixel.  Channel pixels return 0.  Pixels with no downstream
    //   channel carry nodata.
    // -------------------------------------------------------------------------
    [[nodiscard]] virtual lsdtt_raster_result get_height_above_channel_raster(
        const std::vector<float>&   elev,
        const std::vector<uint8_t>& elev_mask,
        int threshold_contributing_pixels) const = 0;

    // get_height_above_channel_mask_raster:
    //   As above, but the channel set is supplied explicitly as a full-grid mask
    //   (channel_mask[flat] != 0 marks a channel pixel) rather than derived from
    //   a contributing-pixel threshold. Use this to measure height above a
    //   specific (e.g. in-scope) subset of the channel network, so surfaces are
    //   referenced to the trunk channel of interest rather than any nearby rill.
    //   channel_mask must have flat-index size == n_rows * n_cols.
    [[nodiscard]] virtual lsdtt_raster_result get_height_above_channel_mask_raster(
        const std::vector<float>&    elev,
        const std::vector<uint8_t>&  elev_mask,
        const std::vector<uint8_t>&  channel_mask) const = 0;

    [[nodiscard]] virtual lsdtt_raster_result get_distance_to_channel_raster(
        int threshold_contributing_pixels) const = 0;

    // -------------------------------------------------------------------------
    // Ridge-relative height and distance
    //
    // For every pixel, find the farthest upslope headwater — the pixel in its
    // upslope set that has the longest D8 flow path back to the current pixel.
    //
    // get_height_to_farthest_ridge_raster:
    //   elev[farthest_upslope_pixel] - elev[pixel].
    //   Positive for normal hillslopes (ridge is higher than pixel).
    //   elev / elev_mask must have flat-index size == n_rows * n_cols.
    //
    // get_distance_to_farthest_ridge_raster:
    //   D8 flow-path length (metres) from the farthest upslope headwater to
    //   the pixel.  Headwater pixels return 0.
    // -------------------------------------------------------------------------
    [[nodiscard]] virtual lsdtt_raster_result get_height_to_farthest_ridge_raster(
        const std::vector<float>&   elev,
        const std::vector<uint8_t>& elev_mask) const = 0;

    [[nodiscard]] virtual lsdtt_raster_result get_distance_to_farthest_ridge_raster()
        const = 0;

    // -------------------------------------------------------------------------
    // Downstream path sum to nearest sink
    //
    // Forward O(N) pass (preorder: receiver before donor).
    // For each non-sink pixel: output = values[self] + output[receiver], stopping
    // when the receiver is in sink_mask.  Sink pixels carry nodata in the output
    // (the connectivity index is undefined at channel cells).
    //
    // If valid_mask[i] == 0 the pixel contributes 0 to the running sum but the
    // path continues through it — callers should pre-clip nodata inputs to a safe
    // value and set valid_mask only where the result should also be nodata.
    //
    // values, valid_mask and sink_mask must have flat-index size == n_rows*n_cols.
    // -------------------------------------------------------------------------
    [[nodiscard]] virtual lsdtt_raster_result get_downstream_path_sum_raster(
        const std::vector<float>&   values,
        const std::vector<uint8_t>& valid_mask,
        const std::vector<uint8_t>& sink_mask) const = 0;

    [[nodiscard]] virtual lsdtt_raster_result get_upslope_sum_raster(
        const std::vector<float>&   values,
        const std::vector<uint8_t>& valid_mask,
        bool include_self = true) const = 0;

    [[nodiscard]] virtual lsdtt_raster_result get_upslope_mean_raster(
        const std::vector<float>&   values,
        const std::vector<uint8_t>& valid_mask,
        bool include_self = true) const = 0;

    // -------------------------------------------------------------------------
    // Flow-routing vector outputs
    // -------------------------------------------------------------------------
    virtual void write_flow_lines(
        const std::string& out_path,
        vector_output_format format,
        const std::string& layer_name = "flow_lines",
        bool spatial_index = true
    ) const = 0;

    virtual void print_arcmap_d8_flowdirection(
        const std::string& out_path,
        raster_output_format format = raster_output_format::cog
    ) const = 0;



    // -------------------------------------------------------------------------
    // Sentinel values — defined in lsdtt_types.hpp
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

// -----------------------------------------------------------------------------
// d8_step_length
//
// Flow-path step length (metres) between two adjacent flat indices.
// Cardinal neighbours return cell_size; diagonal neighbours cell_size*sqrt(2).
// Callable without including lsdtt_flowinfo_cpu.hpp.
// -----------------------------------------------------------------------------
inline double d8_step_length(
    std::int64_t flat_up,
    std::int64_t flat_dn,
    int          n_cols,
    double       cell_size)
{
    const int r_up = static_cast<int>(flat_up / static_cast<std::int64_t>(n_cols));
    const int c_up = static_cast<int>(flat_up % static_cast<std::int64_t>(n_cols));
    const int r_dn = static_cast<int>(flat_dn / static_cast<std::int64_t>(n_cols));
    const int c_dn = static_cast<int>(flat_dn % static_cast<std::int64_t>(n_cols));
    const int dr   = std::abs(r_up - r_dn);
    const int dc   = std::abs(c_up - c_dn);
    return (dr == 1 && dc == 1) ? cell_size * 1.41421356237309504880 : cell_size;
}

// -----------------------------------------------------------------------------
// get_largest_basin_outlet_flat
//
// Scans all outlet cells in the flow object and returns the flat index of the
// one with the largest upslope pixel count.  Throws if no outlet is found.
// -----------------------------------------------------------------------------
[[nodiscard]] std::int64_t get_largest_basin_outlet_flat(
    const lsdtt_flow_base& flow
);

// -----------------------------------------------------------------------------
// get_all_basin_outlet_flats
//
// Scans all cells in the flow object and returns the flat indices of every
// basin outlet (a valid cell whose receiver is k_outlet_receiver) whose upslope
// pixel count is at least min_upslope_pixels. Used to seed chi extraction from
// every basin large enough to qualify, so the largest complete basins can be
// selected downstream. Pass min_upslope_pixels = source_pixels_threshold so
// only outlets that meet the channel-head threshold are returned. Returns an
// empty vector if no qualifying outlet is found.
// -----------------------------------------------------------------------------
[[nodiscard]] std::vector<std::int64_t> get_all_basin_outlet_flats(
    const lsdtt_flow_base& flow,
    std::int64_t           min_upslope_pixels = 0
);

// -----------------------------------------------------------------------------
// is_flat_on_dem_edge
//
// Returns true if the cell lies on the DEM's bounding-box border (first/last
// row or column). An outlet on the edge drains off the clipped DEM; an outlet
// that is NOT on the edge must drain to interior nodata (e.g. an island coast,
// where the sea has been masked).
// -----------------------------------------------------------------------------
[[nodiscard]] bool is_flat_on_dem_edge(
    const lsdtt_flow_base& flow,
    std::int64_t           flat
);

// -----------------------------------------------------------------------------
// move_outlet_upstream
//
// Walks up the main stem from outlet_flat by up to n_pixels steps, choosing at
// each step the donor (a neighbour that flows directly into the current cell)
// with the largest upslope pixel count. Stops early at a channel source. Used
// to lift an interior (coastal) outlet off the nodata fringe so the resulting
// basin no longer touches nodata and passes the completeness test. Returns the
// flat index of the relocated outlet (== outlet_flat if n_pixels <= 0 or no
// donor is found).
// -----------------------------------------------------------------------------
[[nodiscard]] std::int64_t move_outlet_upstream(
    const lsdtt_flow_base& flow,
    std::int64_t           outlet_flat,
    int                    n_pixels
);

} // namespace lsdtt
