// lsdtt_raster_fetch_detail.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// Internals of lsdtt_raster_fetch, exposed so tests can pin them without the
// network. Not part of the public API: callers use fetch_dem().
//
// The pure functions here are what makes the Copernicus path checkable
// offline: the cache key, the Planetary Computer search body and the warp
// argument vector fully determine what that path asks GDAL to do, so pinning
// their exact output pins its behaviour on any machine.
// -----------------------------------------------------------------------------

#pragma once

#include <cstdint>
#include <string>
#include <vector>

#include "lsdtt_raster_fetch.hpp"

namespace lsdtt {
namespace fetch_detail {

// An extent in the target CRS (projection units).
struct fetch_extent
{
    double x_min = 0.0;
    double y_min = 0.0;
    double x_max = 0.0;
    double y_max = 0.0;
};

// The cache file name for a request: everything that changes the pixels.
std::string fetch_cache_key(const lsdtt_fetch_options& options);

// The JSON body of the Planetary Computer STAC search.
std::string pc_stac_search_body(const lsdtt_fetch_options& options);

// The JSON body of a catalogue source's STAC search: its collection for the
// resolved resolution, the bbox, and the fields filter if the row has one.
// For us_3dep (no fields filter) it equals pc_stac_search_body on that
// collection.
std::string catalogue_stac_search_body(const lsdtt_fetch_options& options);

// The STAC collection a catalogue source searches at `resolution_metres`.
std::string catalogue_stac_collection(const lsdtt_dem_source_info& info,
                                      double                       resolution_metres);

// ---- WCS (catalogue rows with a wcs_url) ----

// The native-grid extent to request: the output area (the snapped target
// extent, or the bbox in frame mode) transformed to the native CRS, snapped
// outward to native pixel edges and padded for the resampling kernel.
// Throws std::runtime_error if a transform cannot be built.
fetch_extent wcs_native_extent(const lsdtt_fetch_options&   options,
                               const lsdtt_dem_source_info& info,
                               int                          target_epsg);

// Split an extent on the native grid into chunks of at most max_side pixels
// a side, row by row from the top left. Chunks abut exactly.
std::vector<fetch_extent> wcs_plan_chunks(const fetch_extent& extent, double pixel,
                                          int max_side);

// The GetCoverage URL for one chunk (WCS 1.0.0 or 2.0.1).
std::string wcs_getcoverage_url(const lsdtt_dem_source_info& info,
                                const std::string&           coverage_id,
                                const fetch_extent&          chunk);

// The message of an OGC exception report (OWS 1.1/2.0 ExceptionReport or a
// WCS 1.0 ServiceExceptionReport); empty if the text is not one.
std::string parse_ows_exception(const std::string& xml);

// The payload of the first non-XML part of a multipart MIME body; empty if
// there is none.
std::string multipart_first_payload(const std::string& body, const std::string& boundary);

// Stage one GetCoverage response as a Float32 GeoTIFF at out_path (/vsimem)
// with the native CRS (and nodata 0 where the row says so). Refuses a
// response whose grid is not exactly the requested chunk on the native
// grid: a server that resampled would add a second resampling.
bool wcs_stage_response(const lsdtt_dem_source_info& info,
                        const std::string&           body,
                        const std::string&           content_type,
                        const fetch_extent&          chunk,
                        const std::string&           out_path,
                        std::string&                 error);

// ---- Grid-keyed buckets (catalogue rows with a grid_base_url) ----

// The NZTopo50 sheet ("BQ31") holding an NZTM2000 point; empty outside the
// sheet grid. Sheets are 24 km (E) x 36 km (N): column c spans
// x = 1012000 + 24000 (c - 1), row pair r (24 letters, no I or O; AS is the
// first row) spans y = 6234000 - 36000 (r - AS) downwards.
std::string nztopo50_sheet(double x, double y);

// The Ordnance Survey grid reference of a British National Grid point:
// digits_per_axis 1 gives the 10 km square ("NH24"), 2 the 1 km square
// ("NG2004"). Empty outside the grid.
std::string os_grid_reference(double easting, double northing, int digits_per_axis);

// The 5 km quadrant reference ("NS14NE").
std::string os_grid_quadrant(double easting, double northing);

// The keys of one collection's tiles that intersect an extent (grid_epsg),
// with {model}/{MODEL} filled in, in row-major order.
std::vector<std::string> grid_keys_for_extent(
    const lsdtt_dem_source_info::grid_collection& collection, const fetch_extent& extent,
    const std::string& model);

// The area to look for tiles in, in grid_epsg: the output area padded by
// two output cells. Throws std::runtime_error if a transform fails.
fetch_extent grid_search_extent(const lsdtt_fetch_options&   options,
                                const lsdtt_dem_source_info& info,
                                int                          target_epsg);

// True if this GDAL build can decode LERC-compressed GeoTIFFs.
bool gdal_has_lerc();

// The GDALWarp argument vector. explicit_extent (target CRS) replaces -tap
// with -te on the snapped extent; nullptr gives the Copernicus arguments.
std::vector<std::string> fetch_warp_arguments(const lsdtt_fetch_options& options,
                                              int                        epsg,
                                              const fetch_extent*        explicit_extent = nullptr);

// The envelope, in `epsg`, of the WGS84 box with every edge densified
// (64 points per side), snapped outward to multiples of grid_spacing. Throws
// std::runtime_error if the transform cannot be built.
fetch_extent target_extent_for_bbox(double west, double south, double east, double north,
                                    int epsg, double grid_spacing);

// One STAC item as the catalogue path needs it.
struct stac_item
{
    std::string id;
    std::string href;        // asset href, unsigned
    double      gsd = 0.0;   // 0 if the item carries none
    std::string datetime;
    double      bbox[4] = {0.0, 0.0, 0.0, 0.0};   // WGS84 west, south, east, north
};

// Parse a STAC ItemCollection (search response). Items without the asset are
// skipped and named in `missing_asset`. Returns false (with error) only for
// unparseable JSON. next_body / next_href: the rel=next link, if any.
// asset_key empty: take the GeoTIFF asset whose gsd (or eo:gsd) equals
// asset_gsd, for hosts that key assets by file name.
bool parse_stac_items(const std::string&        json,
                      const std::string&        asset_key,
                      std::vector<stac_item>&   items,
                      std::vector<std::string>& missing_asset,
                      std::string&              next_href,
                      std::string&              next_body,
                      std::string&              error,
                      double                    asset_gsd = 0.0);

// The body for the next page: the original search body with the rel=next
// body's members laid over it. Some hosts (swisstopo) put only a cursor in
// the next body; sent alone it drops the collection and bbox filter.
std::string merge_stac_next_body(const std::string& original, const std::string& next);

// The latest item (by datetime, then id) per tile, the tile being the item
// id with its four-digit year token removed. Output in id order.
std::vector<stac_item> latest_stac_item_per_tile(const std::vector<stac_item>& items);

// Items whose bbox intersects the WGS84 box.
std::vector<stac_item> stac_items_in_bbox(const std::vector<stac_item>& items, double west,
                                          double south, double east, double north);

// Keep the items with gsd == resolution (0 = keep all), in their original
// order. The identity for resolution 0 is what keeps planetary_computer's
// default behaviour unchanged.
std::vector<stac_item> filter_stac_items_by_gsd(const std::vector<stac_item>& items,
                                                double resolution_metres);

// Sorted distinct gsd values (0 means "none given").
std::vector<double> distinct_gsds(const std::vector<stac_item>& items);

// Sort by id: STAC search order is not guaranteed, and later inputs win where
// GDALWarp inputs overlap, so the order is part of the result.
void sort_stac_items(std::vector<stac_item>& items);

// Pre-flight size estimates. Pixels of the output grid for `extent` at
// grid_spacing; bytes downloaded for `n_cells` source cells of `bytes_per_cell`.
std::int64_t estimate_output_pixels(const fetch_extent& extent, double grid_spacing);
double       estimate_download_megabytes(double n_cells, double bytes_per_cell);

// Source cells covered by a WGS84 box on a grid of `cell_degrees`.
double geographic_cell_count(double west, double south, double east, double north,
                             double cell_degrees);

// Write the provenance items into an existing GeoTIFF's default metadata
// domain. Band values, geotransform, CRS and nodata are not touched.
bool stamp_provenance(const std::string&                                       path,
                      const std::vector<std::pair<std::string, std::string>>& items,
                      std::string&                                             error);

// The catalogue path after discovery, on local or remote files: open them,
// refuse non-metre vertical units, warp onto the snapped target grid (square
// cells), check the output and stamp its provenance. options.source must be
// a catalogue source. Used by the offline tests with synthetic files.
lsdtt_fetch_result warp_catalogue_files(const lsdtt_fetch_options&      options,
                                        const std::vector<std::string>& paths,
                                        const std::string&              out_path);

} // namespace fetch_detail
} // namespace lsdtt
