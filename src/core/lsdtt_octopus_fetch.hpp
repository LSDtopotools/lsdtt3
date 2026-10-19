// lsdtt_octopus_fetch.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// Basin-averaged cosmogenic erosion rates from OCTOPUS (octopusdata.org), for a
// bounding box, in the sample format of the erodibility calibration.
//
// OCTOPUS serves its CRN Global compilation through a keyless WFS 2.0.0
// (GeoServer) at http://geoserver.octopusdata.org/geoserver/wfs. That host
// answers over HTTP only (the TLS handshake fails), so the URL is a constant
// here, in one reviewed place; no credential is ever sent, and the data are
// public. Two layers are read:
//   be10-denude:crn_int_basins   basin polygons with every published field
//   be10-denude:crn_int_outlets  outlet points (OBSID1, OBSID2, STUDYID)
// joined on OBSID1.
//
// Request rules (each learned the hard way; see the unit tests):
//   * the box goes as bbox=S,W,N,E,urn:ogc:def:crs:EPSG::4326 with
//     srsName=EPSG:4326: LATITUDE FIRST. Any other form silently returns zero
//     features;
//   * maxFeatures is ignored by the server (one probe pulled about 300 MB), so
//     the count comes from resultType=hits first and the features are paged
//     with count / startIndex, sorted by OBSID1;
//   * a bbox query returns every basin whose polygon INTERSECTS the box, so
//     large basins with outlets far outside it come back too.
//
// No-data sentinels. The field description (Zenodo 5808770) gives -9999 and
// -9999.99; fields with fewer integer digits carry the same all-nines pattern
// truncated (-99.99, -99.999, -9.999, -9.99999 in the live layer). Every
// numeric field of every basin is scanned and these values are nulled (see
// is_octopus_no_data_sentinel), and the counts per field are reported.
//
// Rate schemes (octopus_rate_fields), all in mm/kyr:
//   st                      EBE_MMKYR / EBE_ERR (CAIRN, total 1 sigma); the
//                           internal (AMS-only) uncertainty is ERRBE_AMS
//                           converted with the basin's own EBE_MMKYR / EBE_GCMYR
//   lsdn                    EBEL / EBEL_EXT, internal EBEL_INT
//   lsdn_no_topo_shielding  EBELNS / EBELNS_EXT, internal EBELNS_INT
//   published               BE10EP / BE10EP_ERR (the source papers' own rates)
// and the Al-26 equivalents (EAL_MMKYR / EAL_ERR, EALL*, EALLNS*, AL26EP*).
//
// Data are CC BY 4.0: cite OCTOPUS (Codilean et al. 2018, 2022; Codilean &
// Munack 2025) and every original study (octopus_attribution_text).
//
// GDAL-free header; nothing here prints.
// -----------------------------------------------------------------------------

#pragma once

#include <cstdint>
#include <limits>
#include <map>
#include <string>
#include <vector>

#include "lsdtt_http_fetch.hpp"

namespace lsdtt
{

// ---------------------------------------------------------------------------
// The service (constants; a test pins them)
// ---------------------------------------------------------------------------
[[nodiscard]] std::string octopus_wfs_endpoint();        // http://geoserver.octopusdata.org/geoserver/wfs
[[nodiscard]] std::string octopus_basins_layer();        // be10-denude:crn_int_basins
[[nodiscard]] std::string octopus_outlets_layer();       // be10-denude:crn_int_outlets

// ---------------------------------------------------------------------------
// Request recipe (pure)
// ---------------------------------------------------------------------------
// "S,W,N,E,urn:ogc:def:crs:EPSG::4326" (latitude first), shortest round-trip
// decimal form of each edge.
[[nodiscard]] std::string octopus_wfs_bbox_parameter(double west, double south, double east,
                                                     double north);
// resultType=hits request for a layer and a bbox parameter.
[[nodiscard]] std::string octopus_wfs_hits_url(const std::string& layer, const std::string& bbox);
// One GeoJSON page: count features from start_index, sorted by OBSID1.
[[nodiscard]] std::string octopus_wfs_page_url(const std::string& layer, const std::string& bbox,
                                               int count, std::int64_t start_index);
// The outlets with the given OBSID1 values (CQL_FILTER, no bbox), one page.
[[nodiscard]] std::string octopus_wfs_outlets_by_id_url(const std::vector<std::string>& obsid1);
// The startIndex of every page for n_matched features.
[[nodiscard]] std::vector<std::int64_t> octopus_wfs_page_starts(std::int64_t n_matched, int page_size);
// numberMatched="<n>" of a WFS 2.0 response (hits XML or the GeoJSON key);
// -1 when absent.
[[nodiscard]] std::int64_t parse_octopus_wfs_number_matched(const std::string& body);
// The ExceptionText of an OWS ExceptionReport (all of them joined by "; "),
// or "" when the body is not one.
[[nodiscard]] std::string octopus_wfs_exception_text(const std::string& body);
// True for the no-data sentinels: a negative all-nines value with at least
// three nines (-999, -9999, -9999.99, -99.99, -99.999, -9.999, -9.99999,
// ...), to 1e-6. False for every other value, e.g. -9998.9, -9.9 and 0.
[[nodiscard]] bool is_octopus_no_data_sentinel(double x);

struct lsdtt_octopus_rate_fields
{
    std::string rate;                     // mm/kyr
    std::string uncertainty;              // mm/kyr, 1 sigma
    std::string internal_uncertainty;     // mm/kyr, or the mass-flux field when converted
    std::string internal_mass_flux_rate;  // non-empty: internal = field * rate / this field
    std::string kind;                     // total_1sigma | external_1sigma | published_1sigma
};
// scheme: st | lsdn | lsdn_no_topo_shielding | published; nuclide: be10 | al26.
// Throws std::invalid_argument for anything else, listing the accepted values.
[[nodiscard]] lsdtt_octopus_rate_fields octopus_rate_fields(const std::string& scheme,
                                                            const std::string& nuclide);

// ---------------------------------------------------------------------------
// The fetch
// ---------------------------------------------------------------------------
struct lsdtt_octopus_fetch_options
{
    double west = 0.0, south = 0.0, east = 0.0, north = 0.0;   // WGS84 degrees
    std::string rate_scheme = "st";
    std::string nuclide = "be10";
    int page_size = 100;
    int max_features = 2000;
    bool include_outlets_outside_bbox = false;
    int request_delay_milliseconds = 200;
    lsdtt_http_request_options http;
    // Non-empty: read <dir>/crn_int_basins.json and <dir>/crn_int_outlets.json
    // (the GeoJSON bodies of a previous fetch) instead of the WFS, and
    // <dir>/hits_basins.xml, when present, as the basin hits response.
    std::string response_dir;
};

struct lsdtt_octopus_value
{
    enum class kind { null, integer, real, text };
    kind type = kind::null;
    std::int64_t i = 0;
    double d = 0.0;
    std::string s;
};

struct lsdtt_octopus_field
{
    std::string name;
    lsdtt_octopus_value::kind type = lsdtt_octopus_value::kind::text;   // never null
};

struct lsdtt_octopus_basin
{
    std::vector<lsdtt_octopus_value> values;   // one per result field, sentinels nulled
    std::vector<unsigned char> wkb;            // EPSG:4326, longitude first
    bool has_outlet = false;                   // the outlet was returned
    bool outlet_in_bbox = false;
    double outlet_longitude = std::numeric_limits<double>::quiet_NaN();
    double outlet_latitude = std::numeric_limits<double>::quiet_NaN();
};

struct lsdtt_octopus_csv_row
{
    static constexpr double nan = std::numeric_limits<double>::quiet_NaN();
    std::string sample_id, obsid2, study_id, author, pub_year, ref_doi;
    double latitude = nan, longitude = nan;
    double rate = nan, uncertainty = nan, internal_uncertainty = nan;   // mm/kyr; NaN = blank
    double area_km2 = nan, mean_elevation_m = nan;
    bool outlet_in_bbox = false;
};

struct lsdtt_octopus_study
{
    std::string study_id, author, pub_year, ref_doi, ref_id;   // ref_id: REFID, for studies without a DOI
    std::size_t n_basins = 0;
};

struct lsdtt_octopus_fetch_result
{
    std::vector<lsdtt_octopus_field> fields;          // the basins layer, in layer order
    std::vector<lsdtt_octopus_basin> basins;          // as returned (sorted by OBSID1 online)
    std::int64_t n_matched = 0;                       // basins numberMatched
    std::size_t n_pages = 0;                          // GeoJSON pages read (both layers)
    std::size_t n_requests = 0;                       // HTTP requests made
    std::size_t n_bytes = 0;                          // response bytes read
    std::size_t n_outlets = 0;                        // outlets returned
    std::size_t n_outlet_in_bbox = 0;                 // basins with their outlet in the box
    std::map<std::string, std::size_t> masked_by_field;   // sentinels nulled, per field
    lsdtt_octopus_rate_fields rate_fields;
    std::string rate_scheme, nuclide;
    std::vector<lsdtt_octopus_csv_row> rows;          // the calibration samples
    std::vector<std::string> dropped_no_rate;         // OBSID1 left out: no usable rate
    std::vector<std::string> no_outlet_in_bbox;       // OBSID1 left out: outlet outside / unknown
    std::vector<lsdtt_octopus_study> studies;         // of every basin returned, by study id
    std::string db_version;                           // DBVERNO values, ";"-joined
    std::string access_date;                          // UTC, YYYY-MM-DD
    std::string source;                               // the endpoint, or the response dir
    double west = 0.0, south = 0.0, east = 0.0, north = 0.0;
    std::vector<std::string> warnings;
};

// Fetches (or reads offline) the basins and outlets, masks the sentinels,
// joins and builds the calibration rows. Throws std::runtime_error with
// "OCTOPUS WFS unreachable (<url>): <error>" on a transport failure,
// "OCTOPUS WFS error: <text>" for an ExceptionReport,
// "OCTOPUS returns <n> basins for this box, more than max_features (<m>)",
// and for a body that is neither; std::invalid_argument for a bad box,
// scheme, nuclide, page size or a missing offline file. Zero basins is not
// an error (an empty result).
lsdtt_octopus_fetch_result fetch_octopus_erosion_rates(const lsdtt_octopus_fetch_options& opts);

// The value of a field of a basin ("" / NaN when null or absent).
[[nodiscard]] std::string octopus_basin_text(const lsdtt_octopus_fetch_result& r,
                                             const lsdtt_octopus_basin& b, const std::string& field);
[[nodiscard]] double octopus_basin_number(const lsdtt_octopus_fetch_result& r,
                                          const lsdtt_octopus_basin& b, const std::string& field);

// The calibration CSV (header in lsdtt_octopus_erosion_rate_csv_header()).
[[nodiscard]] std::string octopus_erosion_rate_csv_header();
void write_octopus_erosion_rate_csv(const std::string& path, const lsdtt_octopus_fetch_result& r);

// The basin polygons in EPSG:epsg (layer "octopus_basins"): every field
// returned (the WFS feature id "id" and every OCTOPUS attribute)
// with sentinels nulled, plus outlet_in_bbox, outlet_longitude and
// outlet_latitude. driver: "FlatGeobuf" or "GPKG". Throws on failure.
void write_octopus_basin_polygons(const std::string& path, const std::string& driver, int epsg,
                                  const lsdtt_octopus_fetch_result& r);

// The suggested DEM box: the envelope of the basins with an outlet in the box
// and AREA <= max_area_km2, padded by pad_degrees. ok = false when no basin
// qualifies.
struct lsdtt_octopus_dem_bounds
{
    bool ok = false;
    double west = 0.0, south = 0.0, east = 0.0, north = 0.0;
    std::size_t n_basins = 0;
};
[[nodiscard]] lsdtt_octopus_dem_bounds octopus_dem_bounds(const lsdtt_octopus_fetch_result& r,
                                                          double pad_degrees, double max_area_km2);
// One EPSG:4326 rectangle (layer "dem_bounds") in a FlatGeobuf.
void write_octopus_dem_bounds(const std::string& path, const lsdtt_octopus_dem_bounds& b);

// CC BY 4.0 attribution: OCTOPUS, the endpoint and layers, the access date
// and database version, the citations, the calculator note for the scheme,
// and one line per study.
[[nodiscard]] std::string octopus_attribution_text(const lsdtt_octopus_fetch_result& r);

} // namespace lsdtt
