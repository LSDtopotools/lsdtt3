// lsdtt_hype_parameters.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// HYPE parameter tables and the sub-basin set-up struct
//
// WP2 of the S-HYPE port (plan_potential_shype_port.md, sections 2.10, 5.2.2,
// 5.3 and 7.2). This file holds the inputs of the simplified HYPE column and
// nothing else: plain structs, CSV readers and validation. No process code.
//
// WHAT IS HERE
// ------------
//   hype_general_parameters   HYPE's general (par.txt, one value) parameters
//   hype_landuse_parameters   one row per land-use code (par.txt, per land use)
//   hype_soil_parameters      one row per soil code     (par.txt, per soil type)
//   hype_class_definition     one row per class         (GeoClass.txt)
//   hype_subbasin             one row per sub-basin     (GeoData.txt)
//   hype_subbasin_class       one row per (sub-basin, class) pair, the slc_n
//                             and dhslc_n columns of GeoData.txt
//   hype_setup                all of the above: the seam between the model and
//                             where its inputs come from. The oracle tests fill
//                             it from an lsdtt_input/ directory of CSV tables;
//                             the CLI will fill it from rasters through
//                             lsdtt_subbasin_partition_cpu.
//
// Field names are HYPE's own parameter / column names so that a value can be
// traced straight back to the Fortran. Units are HYPE's: mm, mm per timestep,
// degrees C, m. The model runs daily (plan 7.1), so "per timestep" is "per day".
//
// ABSENT MEANS ZERO
// -----------------
// In HYPE a parameter that is not in par.txt is ZERO, and the oracle cases use
// that to switch processes off. The readers keep that rule for PARAMETERS: a
// parameter missing from general_parameters.csv, or a parameter column missing
// from the land-use or soil table, reads as 0. KEY columns (codes, ids, the
// class's land use, soil and layer count, fractions, areas) are required.
//
// UNKNOWN NAMES ARE ERRORS
// ------------------------
// HYPE silently ignores a parameter name it does not recognise. Combined with
// "absent means zero", a typo such as `mperc` for `mperc1` turns a process off
// without a word. Every reader here therefore rejects a name or column it does
// not know, in every table.
//
// ERRORS
// ------
// Problems reading or parsing a file throw std::runtime_error (mostly from
// lsdtt_csv_table). Content that parses but is not a valid set-up throws
// std::invalid_argument. Every message names the file (or table), the 1-based
// data row and the column. Catch std::exception to handle both.
// -----------------------------------------------------------------------------

#pragma once

#include <array>
#include <cstddef>
#include <string>
#include <vector>

namespace lsdtt
{

// ---------------------------------------------------------------------------
// General parameters (plan 2.10). All default to 0, HYPE's value for a
// parameter that is not given.
// ---------------------------------------------------------------------------
struct hype_general_parameters
{
    double ttpd = 0.0;      // offset of the rain/snow threshold from ttmp [deg C]
    double ttpi = 0.0;      // half-width of the mixed rain/snow interval [deg C]
    double tcalt = 0.0;     // temperature lapse rate [deg C per 100 m]
    double lp = 0.0;        // fraction of field capacity below which AET is reduced [-]
    double epotdist = 0.0;  // decay coefficient of the PET depth distribution [1/m]
    double fsceff = 0.0;    // efficiency of fractional snow cover [-]
    double whcsnow = 0.0;   // liquid-water holding capacity of snow [fraction of snow]
    double rrcs3 = 0.0;     // slope dependence of the soil recession coefficient [1/day per %]
    double rivvel = 0.0;    // river celerity [m/s]
    double damp = 0.0;      // fraction of the river delay taken by the damping box [-]
    double cevpam = 0.0;    // amplitude of the seasonal PET correction [-]
    double cevpph = 0.0;    // phase of the seasonal PET correction [day]
};

// ---------------------------------------------------------------------------
// Land-use parameters, one row per land-use code.
// ---------------------------------------------------------------------------
struct hype_landuse_parameters
{
    int code = 0;           // land-use code (key)
    double ttmp = 0.0;      // threshold temperature for melt and evaporation [deg C]
    double cmlt = 0.0;      // degree-day melt factor [mm / (deg C day)]
    double cevp = 0.0;      // PET coefficient [mm / (deg C day)]
    double srrcs = 0.0;     // surface runoff recession coefficient [1/day]
    double pcluse = 0.0;    // precipitation correction, interception loss [-]
};

// ---------------------------------------------------------------------------
// Soil parameters, one row per soil code. wcwp, wcfc and wcep are STACKED
// volume fractions (plan 2.2), not cumulative contents.
// ---------------------------------------------------------------------------
struct hype_soil_parameters
{
    int code = 0;            // soil code (key)
    double wcwp = 0.0;       // water content below wilting point [-]
    double wcfc = 0.0;       // plant-available water, wilting point to field capacity [-]
    double wcep = 0.0;       // effective porosity, above field capacity [-]
    double rrcs1 = 0.0;      // soil runoff recession coefficient, top layer [1/day]
    double rrcs2 = 0.0;      // soil runoff recession coefficient, bottom layer [1/day]
    double mperc1 = 0.0;     // maximum percolation, layer 1 to 2 [mm/day]
    double mperc2 = 0.0;     // maximum percolation, layer 2 to 3 [mm/day]
    double macrate = 0.0;    // fraction of excess infiltration to macropores [-]
    double srrate = 0.0;     // fraction of excess infiltration to surface runoff [-]
    double mactrinf = 0.0;   // infiltration threshold for macropore/surface runoff [mm/day]
    double mactrsm = 0.0;    // soil-moisture threshold, fraction of field capacity [-]
};

// ---------------------------------------------------------------------------
// Class definition (a GeoClass row): a land use and a soil with a layer
// geometry. Layer bottom depths are below the ground surface; entries beyond
// n_layers are 0.
// ---------------------------------------------------------------------------
struct hype_class_definition
{
    int class_id = 0;                          // class id (key)
    int landuse = 0;                           // land-use code
    int soil = 0;                              // soil code
    int n_layers = 0;                          // number of soil layers, 1..3
    std::array<double, 3> layer_bottom_m{};    // layer bottom depths [m], strictly increasing
    double stream_depth_m = 0.0;               // stream depth d_d [m]
    double tile_depth_m = 0.0;                 // tile drain depth d_t [m]; 0 = none
};

// ---------------------------------------------------------------------------
// Sub-basin (a GeoData row, plan 5.2.2).
// ---------------------------------------------------------------------------
struct hype_subbasin
{
    int subid = 0;                     // sub-basin id (key), > 0
    int maindown = 0;                  // downstream sub-basin id; 0 = outlet of the domain
    double area_m2 = 0.0;              // area [m^2]
    double elev_mean_m = 0.0;          // mean elevation [m]
    double elev_std_m = 0.0;           // elevation standard deviation [m]; carried, unused while fsceff = 0
    double slope_mean_percent = 0.0;   // mean slope, 100 x gradient [%]
    double loc_rivlen_m = 0.0;         // local river length [m]
    double rivlen_m = 0.0;             // main river length [m]
};

// ---------------------------------------------------------------------------
// One (sub-basin, class) pair: HYPE's slc_n and dhslc_n.
// ---------------------------------------------------------------------------
struct hype_subbasin_class
{
    int subid = 0;          // sub-basin id
    int class_id = 0;       // class id
    double fraction = 0.0;  // fraction of the sub-basin area in this class [-]
    double dh_m = 0.0;      // class mean elevation minus sub-basin mean elevation [m]
};

// ---------------------------------------------------------------------------
// The whole set-up. Rows keep the order of their source (file order when read
// from CSV), so a row's index + 1 is its data-row number in validation messages.
// ---------------------------------------------------------------------------
struct hype_setup
{
    hype_general_parameters general;
    std::vector<hype_landuse_parameters> landuse;
    std::vector<hype_soil_parameters> soils;
    std::vector<hype_class_definition> classes;
    std::vector<hype_subbasin> subbasins;
    std::vector<hype_subbasin_class> subbasin_classes;

    // Directory the tables were read from. Used only to name files in
    // validation messages; empty for a set-up built in memory.
    std::string source_dir;
};

// ---------------------------------------------------------------------------
// The recognised parameter names, in the order of the structs above. These
// are the only names accepted in general_parameters.csv and the only
// parameter columns accepted in the land-use and soil tables.
// ---------------------------------------------------------------------------
[[nodiscard]] const std::vector<std::string>& hype_general_parameter_names();
[[nodiscard]] const std::vector<std::string>& hype_landuse_parameter_names();
[[nodiscard]] const std::vector<std::string>& hype_soil_parameter_names();

// ---------------------------------------------------------------------------
// Table readers. Each reads one CSV, rejects unknown names and columns,
// parses, and runs the checks that need only that table (duplicate keys,
// layer depths, fraction range, area, ...). Cross-table checks need the whole
// set-up: see validate_hype_setup.
//
//   general_parameters.csv   columns name,value; one row per parameter
//   landuse_parameters.csv   code + any of hype_landuse_parameter_names()
//   soil_parameters.csv      code + any of hype_soil_parameter_names()
//   classes.csv              class_id,landuse,soil,n_layers,layer1_bottom_m
//                            required; layer2_bottom_m, layer3_bottom_m,
//                            stream_depth_m, tile_depth_m optional. Layers
//                            beyond n_layers must be blank ("0.3,,").
//   subbasins.csv            subid,maindown,area_m2,elev_mean_m,
//                            slope_mean_percent,loc_rivlen_m,rivlen_m
//                            required; elev_std_m optional (0)
//   subbasin_classes.csv     subid,class_id,fraction required; dh_m optional (0)
// ---------------------------------------------------------------------------
[[nodiscard]] hype_general_parameters load_hype_general_parameters(const std::string& path);
[[nodiscard]] std::vector<hype_landuse_parameters> load_hype_landuse_parameters(const std::string& path);
[[nodiscard]] std::vector<hype_soil_parameters> load_hype_soil_parameters(const std::string& path);
[[nodiscard]] std::vector<hype_class_definition> load_hype_classes(const std::string& path);
[[nodiscard]] std::vector<hype_subbasin> load_hype_subbasins(const std::string& path);
[[nodiscard]] std::vector<hype_subbasin_class> load_hype_subbasin_classes(const std::string& path);

// Reads the six tables above from `dir` (an lsdtt_input/ directory) and runs
// validate_hype_setup. forcing.csv in the same directory is not read here; it
// belongs to lsdtt_distributed_hydrology_forcing
// (load_distributed_hydrology_depth_series and
// load_distributed_hydrology_temperature_series).
[[nodiscard]] hype_setup load_hype_setup_directory(const std::string& dir);

// Every check, single-table and cross-table, on a set-up however it was
// built. Throws std::invalid_argument on the first problem found:
//   * duplicate land-use code, soil code, class id, subid or (subid, class_id);
//   * a class whose land use or soil code is not in its table;
//   * n_layers outside 1..3; layer depths not positive and strictly
//     increasing; a depth set beyond n_layers;
//   * negative stream or tile depth;
//   * subid <= 0; area <= 0; maindown not 0 and not a subid; a maindown cycle;
//   * a fraction outside [0, 1]; a (sub-basin, class) row naming an unknown
//     subid or class id; a sub-basin whose fractions do not sum to 1 within
//     1e-6 (this includes a sub-basin with no class rows).
void validate_hype_setup(const hype_setup& setup);

// Sub-basins in downstream (topological) order: every sub-basin comes before
// its maindown, so one forward sweep routes the network (plan 7.3). Returns
// INDICES into `subbasins`. Ties keep input order, so the result is
// deterministic. Throws std::invalid_argument on a maindown that is neither 0
// nor a subid, or on a cycle (naming the sub-basins on it).
[[nodiscard]] std::vector<std::size_t> hype_subbasin_downstream_order(
    const std::vector<hype_subbasin>& subbasins);

} // namespace lsdtt
