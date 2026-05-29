// lsdtt-sediment-routing.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html
//
// -----------------------------------------------------------------------------
// lsdtt-sediment-routing
//
// 1D reach-scale bedload sediment routing on a D8 channel network.
// Reach geometry (length, slope, drainage area) is derived from the DEM.
// Width and grain-size distribution are specified in the parameter file.
// Discharge is estimated from drainage area via a hydraulic-geometry power law
// or from a geomorphic instantaneous unit hydrograph (GIUH).
//
// Outputs
//   <write_fname>_cumulative_deposition.tif  — cumulative deposited volume [m³]
//   <write_fname>_d50_al.tif                 — final active-layer D50 [m]
//
// Usage
//   lsdtt-sediment-routing <path> <parameter_file>
//   lsdtt-sediment-routing <parameter_file>
//   lsdtt-sediment-routing -h
// -----------------------------------------------------------------------------

#include <cstdlib>
#include <iostream>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "../src/core/lsdtt_data_loader.hpp"
#include "../src/core/lsdtt_param_parser.hpp"
#include "../src/core/lsdtt_types.hpp"
#include "../src/cpu/lsdtt_flowinfo_cpu.hpp"
#include "../src/cpu/lsdtt_giuh_cpu.hpp"
#include "../src/cpu/lsdtt_junction_network_cpu.hpp"
#include "../src/cpu/lsdtt_raster_cpu.hpp"
#include "../src/cpu/lsdtt_raster_flow_preprocessing_cpu.hpp"
#include "../src/cpu/lsdtt_sediment_transport_cpu.hpp"

static const std::string k_version  = "0.1.0";
static const std::string k_citation =
static const std::string k_tool_name = "lsdtt-sediment-routing";
    "d-CASCADE: Battista G. et al. (2022) GMD 15:4227-4267\n"
    "lsdtt3:    https://github.com/LSDtopotools/lsdtt3";

// -----------------------------------------------------------------------------
static void print_banner()
{
    std::cout
        << "=====================================================\n"
        << "|| Welcome to lsdtt-sediment-routing              ||\n"
        << "|| 1D reach-scale bedload sediment routing        ||\n"
        << "|| Part of lsdtt3 by Simon M. Mudd                ||\n"
        << "=====================================================\n"
        << "|| Version : " << k_version  << "\n"
        << "|| Citation: " << k_citation << "\n"
        << "=====================================================\n\n";
}

// =============================================================================
// main
// =============================================================================
int main(int argc, char* argv[])
{
    print_banner();

    auto [path, fname] = lsdtt::lsdtt_param_parser::parse_args(
        argc, argv, k_tool_name, k_version, k_citation);

    // -------------------------------------------------------------------------
    // Parameter file
    // -------------------------------------------------------------------------
    lsdtt::lsdtt_param_parser pp;
    try { pp = lsdtt::lsdtt_param_parser(path, fname); }
    catch (const std::exception& e)
    {
        std::cerr << "Error loading parameter file: " << e.what() << "\n";
        return 1;
    }

    // -------------------------------------------------------------------------
    // Parameter declarations
    //
    // Each parameter has its default and its help entry declared on consecutive
    // lines so the full contract for each parameter is visible in one place.
    // Parameters are grouped by workflow step.
    //
    // help_map entry format: { type, default, short_desc, long_desc }
    // -------------------------------------------------------------------------
    std::map<std::string, bool>                     bool_default_map;
    std::map<std::string, float>                    float_default_map;
    std::map<std::string, int>                      int_default_map;
    std::map<std::string, std::string>              string_default_map;
    std::map<std::string, std::vector<std::string>> help_map;

    // --- Threading ---
    int_default_map["n_threads"] = 0;
    help_map["n_threads"] = {"int","0",
        "Number of OpenMP threads.",
        "0 = use all available hardware threads. "
        "Applies to pit filling and flow-routing construction."};

    // --- DEM and flow routing ---
    string_default_map["boundary_conditions"] = "oooo";
    help_map["boundary_conditions"] = {"string","oooo",
        "Boundary conditions for D8 flow routing.",
        "Four-character string for north, east, south, west. "
        "o = open, n = no-flux, p = periodic, b = forced outlet."};

    bool_default_map["raster_is_filled"] = false;
    help_map["raster_is_filled"] = {"bool","false",
        "Whether the input raster is already pit-filled and ready for routing.",
        "If false (default), fill_pits() is called before building flow routing. "
        "If true, the DEM is used as supplied."};

    float_default_map["fill_min_slope"] = 0.0001f;
    help_map["fill_min_slope"] = {"float","0.0001",
        "Minimum slope imposed by pit filling.",
        "Used only when raster_is_filled is false. "
        "Passed to fill_pits() to gently enforce downstream descent."};

    int_default_map["source_pixels_threshold"] = 1000;
    help_map["source_pixels_threshold"] = {"int","1000",
        "Contributing-pixel threshold used to define channel source heads.",
        "Channel heads are extracted from the flow object by this threshold. "
        "Reduce to extract a denser network; increase for sparser."};

    // --- Transport formula ---
    string_default_map["transport_formula_choice"] = "wilcock_crowe";
    help_map["transport_formula_choice"] = {"string","wilcock_crowe",
        "Bedload transport capacity formula.",
        "wilcock_crowe: Wilcock & Crowe (2003) mixed grain-size formula. "
        "Accounts for sand fraction and grain-size-dependent hiding. "
        "wong_parker: Wong & Parker (2006) with Parker (1982) hiding function. "
        "Best for unimodal gravel beds. "
        "engelund_hansen: Engelund & Hansen (1967) total-load formula using D50. "
        "Simple but does not resolve grain-size classes."};

    // --- Routing time control ---
    int_default_map["n_timesteps"] = 365;
    help_map["n_timesteps"] = {"int","365",
        "Number of routing timesteps.",
        "The model marches upstream-to-downstream this many times. "
        "Use 365 with timestep_s = 86400 for one year of daily steps."};

    float_default_map["timestep_s"] = 86400.0f;
    help_map["timestep_s"] = {"float","86400",
        "Duration of each routing timestep in seconds.",
        "Default is one day (86400 s). "
        "Must be consistent with the discharge input if use_giuh_discharge is true."};

    // --- Sediment properties ---
    float_default_map["porosity"] = 0.4f;
    help_map["porosity"] = {"float","0.4",
        "Bed sediment porosity (dimensionless).",
        "Used to convert between bulk volume and solid volume when updating "
        "the active layer. Typical gravel beds: 0.35–0.45."};

    float_default_map["sediment_density_kg_m3"] = 2650.0f;
    help_map["sediment_density_kg_m3"] = {"float","2650",
        "Sediment grain density [kg/m³].",
        "Quartz density is 2650 kg/m³. Used in Shields stress and "
        "transport capacity calculations."};

    float_default_map["active_layer_depth_m"] = -1.0f;
    help_map["active_layer_depth_m"] = {"float","-1",
        "Fixed active-layer depth [m] applied to all reaches.",
        "The active layer is the mobile surface sediment that participates in "
        "exchange each timestep. "
        "Values <= 0 (default) use 2 × D90 per reach (Parker 1991), "
        "which scales with the coarsest grains."};

    // --- Grain-size distribution ---
    string_default_map["psi_mm"] = "4,8,16,32,64";
    help_map["psi_mm"] = {"string","4,8,16,32,64",
        "Grain-class centre diameters [mm], comma-separated.",
        "Defines the grain-size classes used in transport calculations. "
        "Values correspond to phi-scale class centres (4, 8, 16, 32, 64 mm "
        "span fine gravel to coarse gravel). "
        "Must have the same number of entries as gsd_fractions."};

    string_default_map["gsd_fractions"] = "0.10,0.25,0.35,0.20,0.10";
    help_map["gsd_fractions"] = {"string","0.10,0.25,0.35,0.20,0.10",
        "Mass fraction per grain class, comma-separated (must sum to 1).",
        "Applied uniformly to all reaches unless overridden. "
        "Must have the same number of entries as psi_mm."};

    float_default_map["D50_m"] = 0.016f;
    help_map["D50_m"] = {"float","0.016",
        "Median grain diameter [m] used in transport and hydraulic calculations.",
        "Used directly by the Engelund-Hansen formula and as the reference "
        "grain size for hiding-function calculations in Wong-Parker."};

    float_default_map["D84_m"] = 0.045f;
    help_map["D84_m"] = {"float","0.045",
        "D84 grain diameter [m].",
        "84th percentile grain diameter. Used in Wilcock-Crowe hiding function."};

    float_default_map["D90_m"] = 0.055f;
    help_map["D90_m"] = {"float","0.055",
        "D90 grain diameter [m].",
        "90th percentile grain diameter. Sets the default active-layer depth "
        "(2 × D90) when active_layer_depth_m <= 0."};

    float_default_map["mannings_n"] = 0.04f;
    help_map["mannings_n"] = {"float","0.04",
        "Manning's roughness coefficient for flow-depth estimation.",
        "Used to compute normal flow depth from discharge, width, and slope "
        "via Manning's equation. Typical gravel-bed rivers: 0.03–0.07."};

    // --- Channel width (hydraulic geometry) ---
    float_default_map["width_coefficient"] = 2.76f;
    help_map["width_coefficient"] = {"float","2.76",
        "Hydraulic-geometry coefficient for channel width.",
        "Width [m] = width_coefficient × drainage_area [m²]^width_exponent. "
        "Default (2.76) is from Leopold & Maddock (1953) for gravel-bed rivers."};

    float_default_map["width_exponent"] = 0.40f;
    help_map["width_exponent"] = {"float","0.40",
        "Hydraulic-geometry exponent for channel width.",
        "Typical range: 0.35–0.50 for natural rivers."};

    // --- Discharge (hydraulic geometry fallback) ---
    float_default_map["discharge_coefficient"] = 1.0e-5f;
    help_map["discharge_coefficient"] = {"float","1e-5",
        "Hydraulic-geometry coefficient for mean discharge when GIUH is not used.",
        "Q [m³/s] = discharge_coefficient × drainage_area [m²]^discharge_exponent. "
        "This is a simple power-law scaling. Use use_giuh_discharge for a "
        "physics-based estimate."};

    float_default_map["discharge_exponent"] = 0.80f;
    help_map["discharge_exponent"] = {"float","0.80",
        "Hydraulic-geometry exponent for mean discharge.",
        "Typical range: 0.7–0.9. Used only when use_giuh_discharge is false."};

    // --- GIUH discharge estimation ---
    bool_default_map["use_giuh_discharge"] = false;
    help_map["use_giuh_discharge"] = {"bool","false",
        "Use the geomorphic IUH to estimate discharge for sediment routing.",
        "When true, the width-function-based IUH (Rodriguez-Iturbe & Valdes 1979) "
        "is computed from the channel network geometry. A synthetic rainfall pulse "
        "is convolved with the IUH to produce per-reach discharge at each timestep. "
        "When false, discharge is estimated from drainage area using the "
        "discharge_coefficient and discharge_exponent power law."};

    float_default_map["giuh_v_ref_m_per_s"] = 1.0f;
    help_map["giuh_v_ref_m_per_s"] = {"float","1.0",
        "Reference channel velocity [m/s] for GIUH travel-time estimation.",
        "Travel velocity in each reach is scaled from this reference: "
        "v = v_ref × (slope / slope_ref)^velocity_exponent. "
        "Typical bankfull channel velocities are 0.5–2.0 m/s."};

    float_default_map["giuh_slope_ref"] = 0.01f;
    help_map["giuh_slope_ref"] = {"float","0.01",
        "Reference slope corresponding to giuh_v_ref_m_per_s.",
        "A reach with slope equal to giuh_slope_ref is assigned velocity "
        "giuh_v_ref_m_per_s. Steeper reaches are faster; shallower are slower."};

    float_default_map["giuh_velocity_exponent"] = 0.5f;
    help_map["giuh_velocity_exponent"] = {"float","0.5",
        "Exponent β in the kinematic-wave velocity scaling v ∝ slope^β.",
        "β = 0.5 is consistent with Manning's equation and the kinematic wave "
        "approximation. Typical range: 0.3–0.6."};

    float_default_map["giuh_hillslope_length_m"] = -1.0f;
    help_map["giuh_hillslope_length_m"] = {"float","-1",
        "Mean hillslope length [m] added as overland-flow travel time per reach.",
        "When > 0, an overland-flow travel time = hillslope_length / v_hillslope "
        "is added to each reach, delaying the IUH peak. "
        "Values <= 0 disable the hillslope contribution."};

    float_default_map["giuh_v_hillslope_m_per_s"] = 0.005f;
    help_map["giuh_v_hillslope_m_per_s"] = {"float","0.005",
        "Overland-flow velocity [m/s] for hillslope travel time.",
        "Only used when giuh_hillslope_length_m > 0. "
        "Typical sheet-flow velocities: 0.001–0.01 m/s."};

    int_default_map["giuh_n_bins"] = 200;
    help_map["giuh_n_bins"] = {"int","200",
        "Number of time bins used to discretise the IUH.",
        "More bins give a smoother IUH but use more memory. "
        "200 bins at 1-hour resolution covers 200 hours of response time."};

    float_default_map["giuh_dt_bin_s"] = 3600.0f;
    help_map["giuh_dt_bin_s"] = {"float","3600",
        "Width of each IUH time bin in seconds.",
        "Set to match or be smaller than timestep_s for accurate convolution. "
        "Default is 1 hour (3600 s)."};

    float_default_map["giuh_runoff_coefficient"] = 1.0f;
    help_map["giuh_runoff_coefficient"] = {"float","1.0",
        "Fraction of rainfall that becomes runoff (dimensionless, 0–1).",
        "Applied to the rainfall rate before convolving with the IUH. "
        "1.0 = all rain reaches the channel (no infiltration losses). "
        "Lower values reduce the discharge estimate."};

    float_default_map["giuh_rainfall_rate_m_per_s"] = 1.0e-5f;
    help_map["giuh_rainfall_rate_m_per_s"] = {"float","1e-5",
        "Rainfall rate [m/s] for the synthetic IUH pulse (≈ 36 mm/hr).",
        "A uniform pulse of this intensity lasting giuh_pulse_duration_s "
        "is convolved with the IUH to produce the discharge time series. "
        "Adjust to represent the design storm relevant to your setting."};

    float_default_map["giuh_pulse_duration_s"] = 86400.0f;
    help_map["giuh_pulse_duration_s"] = {"float","86400",
        "Duration of the synthetic rainfall pulse in seconds.",
        "The pulse is applied uniformly over this duration at giuh_rainfall_rate_m_per_s. "
        "Default is 86400 s (one day), matching the default timestep_s."};

    // --- Output ---
    string_default_map["write_output_format"] = "cog";
    help_map["write_output_format"] = {"string","cog",
        "Output raster format.",
        "cog: Cloud Optimised GeoTIFF (default). "
        "gtiff: plain compressed GeoTIFF. "
        "same_as_input: match the input raster format."};

    bool_default_map["write_deposition_raster"] = true;
    help_map["write_deposition_raster"] = {"bool","true",
        "Write cumulative sediment deposition [m³] per reach as a raster.",
        "Each channel pixel in a reach is assigned the cumulative deposition "
        "volume for that reach summed over all timesteps. "
        "Output file: <write_fname>_cumulative_deposition.tif."};

    bool_default_map["write_d50_raster"] = true;
    help_map["write_d50_raster"] = {"bool","true",
        "Write the final active-layer D50 [m] per reach as a raster.",
        "D50 is computed from the grain-size distribution of the active layer "
        "at the end of the last timestep. "
        "Output file: <write_fname>_d50_al.tif."};

    // -------------------------------------------------------------------------
    // Parse — fills maps from the parameter file, falling back to defaults
    // -------------------------------------------------------------------------
    pp.parse_all_parameters(
        float_default_map,
        int_default_map,
        bool_default_map,
        string_default_map
    );

    std::map<std::string, bool>        bf  = pp.get_bool_parameters();
    std::map<std::string, float>       ff  = pp.get_float_parameters();
    std::map<std::string, int>         ifs = pp.get_int_parameters();
    std::map<std::string, std::string> sf  = pp.get_string_parameters();

    // -------------------------------------------------------------------------
    // Build routing options
    // -------------------------------------------------------------------------
    lsdtt::lsdtt_sediment_routing_options opts;

    opts.n_timesteps            = ifs["n_timesteps"];
    opts.timestep_s             = static_cast<double>(ff["timestep_s"]);
    opts.porosity               = static_cast<double>(ff["porosity"]);
    opts.sediment_density_kg_m3 = static_cast<double>(ff["sediment_density_kg_m3"]);
    opts.active_layer_depth_m   = static_cast<double>(ff["active_layer_depth_m"]);
    opts.D50_m                  = static_cast<double>(ff["D50_m"]);
    opts.D84_m                  = static_cast<double>(ff["D84_m"]);
    opts.D90_m                  = static_cast<double>(ff["D90_m"]);
    opts.mannings_n             = static_cast<double>(ff["mannings_n"]);
    opts.width_coefficient      = static_cast<double>(ff["width_coefficient"]);
    opts.width_exponent         = static_cast<double>(ff["width_exponent"]);
    opts.discharge_coefficient  = static_cast<double>(ff["discharge_coefficient"]);
    opts.discharge_exponent     = static_cast<double>(ff["discharge_exponent"]);
    opts.n_threads              = ifs["n_threads"];

    if      (sf["transport_formula_choice"] == "wong_parker")
        opts.formula = lsdtt::lsdtt_tr_cap_formula::wong_parker;
    else if (sf["transport_formula_choice"] == "engelund_hansen")
        opts.formula = lsdtt::lsdtt_tr_cap_formula::engelund_hansen;
    else
        opts.formula = lsdtt::lsdtt_tr_cap_formula::wilcock_crowe;

    auto psi_v  = pp.parse_double_vector("psi_mm");
    auto frac_v = pp.parse_double_vector("gsd_fractions");
    if (!psi_v.empty() && psi_v.size() == frac_v.size())
    {
        opts.psi_mm        = psi_v;
        opts.gsd_fractions = frac_v;
    }
    else if (!psi_v.empty())
    {
        std::cerr << "[warning] psi_mm and gsd_fractions must have the same "
                     "number of values. Using defaults.\n";
    }

    lsdtt::raster_output_format out_fmt = lsdtt::raster_output_format::cog;
    if      (sf["write_output_format"] == "gtiff")
        out_fmt = lsdtt::raster_output_format::gtiff;
    else if (sf["write_output_format"] == "same_as_input")
        out_fmt = lsdtt::raster_output_format::same_as_input;

    // -------------------------------------------------------------------------
    // Load DEM
    // -------------------------------------------------------------------------
    const std::string input_path = pp.get_full_read_path();
    std::cout << "Opening DEM: " << input_path << "\n";

    lsdtt::lsdtt_data_loader loader;
    try { loader.open_raster(input_path); }
    catch (const std::exception& e)
    {
        std::cerr << "Error opening DEM: " << e.what() << "\n";
        return 1;
    }

    loader.print_raster_info();

    try { loader.get_meta().check_projected_metres(); }
    catch (const std::exception& e)
    {
        std::cerr << "FATAL: " << e.what() << "\n";
        return 1;
    }

    lsdtt::lsdtt_raster_cpu dem(loader);
    lsdtt::lsdtt_raster_cpu routed_dem = dem;

    // -------------------------------------------------------------------------
    // Pit filling
    // -------------------------------------------------------------------------
    if (!bf["raster_is_filled"])
    {
        std::cout << "Pit filling (min slope " << ff["fill_min_slope"] << ")...\n";
        try
        {
            routed_dem = lsdtt::fill_pits(dem,
                static_cast<double>(ff["fill_min_slope"]),
                ifs["n_threads"]);
        }
        catch (const std::exception& e)
        {
            std::cerr << "FATAL: pit filling failed: " << e.what() << "\n";
            return 1;
        }
    }

    // -------------------------------------------------------------------------
    // D8 flow routing
    // -------------------------------------------------------------------------
    std::cout << "Building D8 flow routing...\n";

    lsdtt::lsdtt_flowinfo_options fi_opts;
    fi_opts.boundary_conditions       = sf["boundary_conditions"];
    fi_opts.interior_nodata_is_outlet = true;
    fi_opts.require_acyclic_graph     = true;
    fi_opts.build_topological_order   = true;
    fi_opts.build_s_vector            = true;
    fi_opts.n_threads                 = ifs["n_threads"];

    lsdtt::lsdtt_flowinfo flow;
    try
    {
        flow = lsdtt::lsdtt_flowinfo(
            routed_dem.get_data_array(),
            routed_dem.get_mask_array(),
            routed_dem.get_meta(),
            fi_opts);
    }
    catch (const std::exception& e)
    {
        std::cerr << "FATAL: flow routing failed: " << e.what() << "\n";
        return 1;
    }

    // -------------------------------------------------------------------------
    // Source extraction
    // -------------------------------------------------------------------------
    std::cout << "Extracting channel sources (threshold = "
              << ifs["source_pixels_threshold"] << " px)...\n";

    std::vector<std::int64_t> source_flats;
    try
    {
        source_flats = flow.get_sources_index_threshold(
            ifs["source_pixels_threshold"]);
    }
    catch (const std::exception& e)
    {
        std::cerr << "FATAL: source extraction failed: " << e.what() << "\n";
        return 1;
    }

    if (source_flats.empty())
    {
        std::cerr << "FATAL: no source pixels found. "
                     "Try lowering source_pixels_threshold.\n";
        return 1;
    }
    std::cout << "  Source pixels: " << source_flats.size() << "\n";

    // -------------------------------------------------------------------------
    // Junction / channel network
    // -------------------------------------------------------------------------
    std::cout << "Building channel network...\n";

    lsdtt::lsdtt_junction_network_options jn_opts;
    jn_opts.build_topological_order = true;
    jn_opts.build_s_vector          = true;

    lsdtt::lsdtt_junction_network_cpu jn;
    try
    {
        jn = lsdtt::lsdtt_junction_network_cpu(flow, source_flats, jn_opts);
    }
    catch (const std::exception& e)
    {
        std::cerr << "FATAL: channel-network construction failed: "
                  << e.what() << "\n";
        return 1;
    }

    std::cout << "  Junctions (reaches): " << jn.n_junctions() << "\n";

    // -------------------------------------------------------------------------
    // Build sediment reaches from DEM geometry
    // -------------------------------------------------------------------------
    std::cout << "Building sediment reaches from DEM...\n";
    auto reaches = lsdtt::build_sediment_reaches(jn, flow, dem, opts);

    // -------------------------------------------------------------------------
    // Initialise active layers
    // -------------------------------------------------------------------------
    auto active_layers = lsdtt::initialise_active_layers(reaches, opts);

    // -------------------------------------------------------------------------
    // GIUH discharge estimation (optional)
    // -------------------------------------------------------------------------
    std::vector<std::vector<double>> discharge_giuh;

    if (bf["use_giuh_discharge"])
    {
        std::cout << "Computing GIUH...\n";

        lsdtt::lsdtt_giuh_options giuh_opts;
        giuh_opts.v_ref_m_per_s       = static_cast<double>(ff["giuh_v_ref_m_per_s"]);
        giuh_opts.slope_ref           = static_cast<double>(ff["giuh_slope_ref"]);
        giuh_opts.velocity_exponent   = static_cast<double>(ff["giuh_velocity_exponent"]);
        giuh_opts.hillslope_length_m  = static_cast<double>(ff["giuh_hillslope_length_m"]);
        giuh_opts.v_hillslope_m_per_s = static_cast<double>(ff["giuh_v_hillslope_m_per_s"]);
        giuh_opts.n_bins              = ifs["giuh_n_bins"];
        giuh_opts.dt_bin_s            = static_cast<double>(ff["giuh_dt_bin_s"]);
        giuh_opts.runoff_coefficient  = static_cast<double>(ff["giuh_runoff_coefficient"]);
        giuh_opts.min_slope           = opts.min_slope;

        try
        {
            lsdtt::lsdtt_giuh_result giuh = lsdtt::compute_giuh(jn, reaches, giuh_opts);

            std::cout << "  IUH peak      : " << giuh.iuh_peak << " s⁻¹\n"
                      << "  Time of peak  : " << giuh.time_of_peak_s / 3600.0
                      << " hours\n"
                      << "  Catchment area: " << giuh.catchment_area_m2 / 1.0e6
                      << " km²\n";

            // Convolve a synthetic uniform rainfall pulse with the IUH
            std::vector<double> Q_giuh = lsdtt::giuh_pulse_discharge(
                giuh,
                static_cast<double>(ff["giuh_pulse_duration_s"]),
                static_cast<double>(ff["giuh_rainfall_rate_m_per_s"]),
                opts.timestep_s,
                static_cast<double>(ff["giuh_runoff_coefficient"]));

            const int n_Q = static_cast<int>(Q_giuh.size());
            const double total_A = giuh.catchment_area_m2;

            for (int t = 0; t < opts.n_timesteps; ++t)
            {
                int q_idx = (n_Q > 0) ? (t % n_Q) : 0;
                double q  = (n_Q > 0) ? Q_giuh[q_idx] : 0.0;
                std::vector<double> row(jn.n_junctions());
                for (int j = 0; j < jn.n_junctions(); ++j)
                {
                    double frac = (total_A > 0.0)
                        ? reaches[j].drainage_area_m2 / total_A : 1.0;
                    row[j] = q * frac;
                }
                discharge_giuh.push_back(std::move(row));
            }
        }
        catch (const std::exception& e)
        {
            std::cerr << "Warning: GIUH computation failed: " << e.what()
                      << "\nFalling back to hydraulic-geometry discharge.\n";
            discharge_giuh.clear();
        }
    }

    // -------------------------------------------------------------------------
    // Run routing
    // -------------------------------------------------------------------------
    std::cout << "Running sediment routing ("
              << opts.n_timesteps << " timesteps, "
              << opts.timestep_s  << " s each)...\n";

    lsdtt::lsdtt_sediment_routing_result result;
    try
    {
        result = lsdtt::run_sediment_routing(
            jn, reaches, active_layers, discharge_giuh, opts);
    }
    catch (const std::exception& e)
    {
        std::cerr << "FATAL: sediment routing failed: " << e.what() << "\n";
        return 1;
    }

    std::cout << "  Total sediment flux out of network: "
              << result.total_flux_out_m3 << " m³\n";

    // -------------------------------------------------------------------------
    // Write outputs
    // -------------------------------------------------------------------------
    const std::string write_prefix = pp.get_full_write_path();

    if (bf["write_deposition_raster"])
    {
        std::cout << "Writing deposition raster...\n";
        try
        {
            lsdtt::lsdtt_raster_cpu dep_raster =
                lsdtt::make_deposition_raster(jn, flow, dem, result);
            const std::string dep_path = write_prefix + "_cumulative_deposition.tif";
            loader.write_full(
                lsdtt::lsdtt_raster_result{
                    dep_raster.get_data(),
                    dep_raster.get_mask(),
                    dep_raster.get_meta()},
                dep_path, out_fmt);
            std::cout << "  Written: " << dep_path << "\n";
        }
        catch (const std::exception& e)
        {
            std::cerr << "Warning: deposition raster write failed: "
                      << e.what() << "\n";
        }
    }

    if (bf["write_d50_raster"])
    {
        std::cout << "Writing D50 raster...\n";
        try
        {
            lsdtt::lsdtt_raster_cpu d50_raster =
                lsdtt::make_d50_raster(jn, flow, dem, result);
            const std::string d50_path = write_prefix + "_d50_al.tif";
            loader.write_full(
                lsdtt::lsdtt_raster_result{
                    d50_raster.get_data(),
                    d50_raster.get_mask(),
                    d50_raster.get_meta()},
                d50_path, out_fmt);
            std::cout << "  Written: " << d50_path << "\n";
        }
        catch (const std::exception& e)
        {
            std::cerr << "Warning: D50 raster write failed: "
                      << e.what() << "\n";
        }
    }

    // -------------------------------------------------------------------------
    // Help file + parameter echo
    // -------------------------------------------------------------------------
    pp.print_help(
        help_map,
        pp.get_full_write_path() + "_README",
        k_version,
        k_citation
    );

    pp.print_parameters();

    return 0;
}
