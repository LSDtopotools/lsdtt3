// lsdtt-cosmo-landslides.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
//
// -----------------------------------------------------------------------------
// Triggered landslides with an imposed size distribution.
//
// The user prescribes WHEN (a trigger time, or a recurrence) and HOW MUCH (an
// event volume, or a long-term erosion rate). This tool answers how many
// landslides that implies, what sizes they are, and where they go -- then
// excavates them and writes the fields the cosmogenic and sediment-routing
// tools consume.
//
// Sizes are IMPOSED, not emergent: they are sampled from a prescribed
// distribution (Malamud inverse gamma by default) rather than derived from
// slope stability. That is deliberate -- the purpose is to impose a known
// landslide flux and ask what it does to something else, most immediately to
// cosmogenic nuclide concentrations.
//
// Landslides come from HILLSLOPES. Initiation is restricted to pixels whose
// contributing area is below the channel initiation threshold, and a landslide
// is built from that pixel's upslope set, so no scar can reach into a channel
// and no landslide can be larger than the hillslope it comes from.
// -----------------------------------------------------------------------------

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <vector>

#include "../src/core/lsdtt_data_loader.hpp"
#include "../src/core/lsdtt_flow_base.hpp"
#include "../src/core/lsdtt_param_parser.hpp"
#include "../src/cpu/lsdtt_crn_parameters.hpp"
#include "../src/cpu/lsdtt_landslide_crn_cpu.hpp"
#include "../src/cpu/lsdtt_landslide_event_cpu.hpp"
#include "../src/cpu/lsdtt_raster_cpu.hpp"
#include "../src/cpu/lsdtt_raster_flow_preprocessing_cpu.hpp"

static const std::string k_version   = "0.1.0";
static const std::string k_citation  = "https://github.com/LSDtopotools/lsdtt3";
static const std::string k_tool_name = "lsdtt-cosmo-landslides";

// ---------------------------------------------------------------------------
static void print_banner()
{
    std::cout
        << "========================================================\n"
        << "|| Welcome to lsdtt-cosmo-landslides                        ||\n"
        << "|| Triggered landslides, imposed size distribution    ||\n"
        << "|| Part of lsdtt3 by Simon M. Mudd                    ||\n"
        << "========================================================\n"
        << "|| CLI version: " << k_version  << "\n"
        << "|| Citation   : " << k_citation << "\n"
        << "|| -v or --version : full lsdtt3 release version\n"
        << "|| -h              : help explaining the parameters\n"
        << "========================================================\n\n";
}

// ---------------------------------------------------------------------------
// Cosmogenic parameters this tool can actually service.
//
// VALIDATE UP FRONT AND FATALLY. Until 2026-08-05 a bad muon_scheme or nuclide
// was discovered only when the CRN block ran, was caught, downgraded to a
// warning on stderr, and the tool then exited **0** having written no
// concentration raster. To a scripted workflow that is a success that silently
// produced nothing -- worse than a crash, because nothing notices.
//
// CAPABILITY, NOT JUST VALIDITY. lsdtt_crn_parameters::is_valid_muon_scheme()
// accepts every scheme lsdtt3 knows, including cronusEmulatorSt and
// cronusEmulatorLSDn. Those are pressure- AND nuclide-dependent, and this tool
// reaches the column through reset_scaling(), which deliberately throws for
// them because it is handed neither. So validating against is_valid_muon_scheme
// alone would accept a name that fails later -- the same trap
// lsdtt-cosmo-column-tool documents and guards against. Checked against the
// serviceable subset here instead.
//
// TO BE REVISITED AT PLAN PHASE R2c: once this tool builds a real site with a
// pressure raster, the emulator schemes become serviceable and this rejection
// must be lifted rather than left to reject something that now works.
static bool landslide_cosmo_parameters_supported(const std::string& muon_scheme,
                                                 const std::string& nuclide)
{
    if (muon_scheme == "cronusEmulatorSt" || muon_scheme == "cronusEmulatorLSDn")
    {
        std::cerr << "FATAL: muon_scheme \"" << muon_scheme
                  << "\" is not supported by " << k_tool_name << ".\n"
                  << "       Its muon terms depend on site pressure and "
                     "nuclide, and this tool applies\n"
                  << "       a scalar production scaling with no pressure "
                     "raster. Use lsdtt-cosmo-prediction\n"
                  << "       for the CRONUS v3 emulator schemes.\n"
                  << "       Schemes available here: Braucher, "
                     "BraucherBorchers, Granger, Schaller, newCRONUS.\n";
        return false;
    }
    if (!lsdtt::lsdtt_crn_parameters::is_valid_muon_scheme(muon_scheme))
    {
        std::cerr << "FATAL: unknown muon_scheme \"" << muon_scheme
                  << "\". Allowed values: "
                  << lsdtt::lsdtt_crn_parameters::allowed_muon_schemes_csv()
                  << ".\n";
        return false;
    }
    if (nuclide != "Be10" && nuclide != "Al26" && nuclide != "C14")
    {
        std::cerr << "FATAL: unknown nuclide \"" << nuclide
                  << "\". Allowed values: Be10, Al26, C14.\n";
        return false;
    }
    return true;
}

// Build a raster matching `ref` in geometry, filled with `value` and carrying
// ref's valid mask. Both output fields need this, and neither may be built with
// an empty mask -- the event writer sets data only.
static lsdtt::lsdtt_raster_cpu like(const lsdtt::lsdtt_raster_cpu& ref, float value)
{
    const lsdtt::lsdtt_raster_meta& meta = ref.get_meta();
    lsdtt::lsdtt_array2d_f  d(meta.n_rows, meta.n_cols, value);
    lsdtt::lsdtt_array2d_u8 m(meta.n_rows, meta.n_cols, 0);
    const auto& rm = ref.get_mask_array();
    for (int r = 0; r < meta.n_rows; ++r)
        for (int c = 0; c < meta.n_cols; ++c)
            if (rm(r, c)) m(r, c) = 1;
    return lsdtt::lsdtt_raster_cpu(std::move(d), std::move(m), meta);
}


int main(int argc, char* argv[])
{
    print_banner();

    lsdtt::lsdtt_cli_args cli_args = lsdtt::lsdtt_param_parser::parse_args_deferred(
        argc, argv, k_tool_name, k_version, k_citation);

    std::map<std::string,float>       float_default_map;
    std::map<std::string,int>         int_default_map;
    std::map<std::string,bool>        bool_default_map;
    std::map<std::string,std::string> string_default_map;
    lsdtt::lsdtt_help_registry help_map;

    // --- DEM preparation ---
    bool_default_map["fill_before_routing"] = true;
    help_map["fill_before_routing"] = {"bool","true",
        "Fill DEM pits before building flow routing.",
        "Set false only if the DEM is already pit-filled.",
        "DEM preparation","param",""};

    float_default_map["min_slope_for_fill"] = 0.0001f;
    help_map["min_slope_for_fill"] = {"float","0.0001",
        "Minimum imposed slope during pit filling (m/m).","",
        "DEM preparation","param",""};

    string_default_map["boundary_conditions"] = "oooo";
    help_map["boundary_conditions"] = {"string","oooo",
        "D8 flow routing boundary conditions.",
        "Four-character string (N,E,S,W): o=open, b=forced outlet.",
        "DEM preparation","param",""};

    string_default_map["write_output_format"] = "cog";
    help_map["write_output_format"] = {"string","cog",
        "Raster output format: cog, gtiff, or same_as_input.","",
        "Output","param",""};

    // --- Where landslides may start ---
    int_default_map["threshold_contributing_pixels"] = 1000;
    help_map["threshold_contributing_pixels"] = {"int","1000",
        "Channel initiation threshold, in pixels.",
        "Landslides start only BELOW this, i.e. on hillslopes. It also caps the "
        "largest landslide, since a landslide cannot exceed the hillslope it "
        "comes from.",
        "Placement","param",""};

    string_default_map["mask_fname"] = "NULL";
    help_map["mask_fname"] = {"string","NULL",
        "Optional raster restricting where landslides may occur.",
        "Non-zero and valid = active. NULL uses the whole DEM.",
        "Placement","param",""};

    // --- Size distribution ---
    string_default_map["landslide_size_model"] = "inverse_gamma";
    help_map["landslide_size_model"] = {"string","inverse_gamma",
        "Size distribution family: inverse_gamma or truncated_pareto.",
        "inverse_gamma (Malamud et al. 2004) has a power-law tail AND a "
        "small-area rollover. truncated_pareto (Yanites et al. 2009) is a pure "
        "power law, kept to reproduce published results.",
        "Size distribution","param",""};

    float_default_map["inverse_gamma_rho"] = 1.40f;
    help_map["inverse_gamma_rho"] = {"float","1.40",
        "Inverse gamma rho. Tail exponent is -(rho+1), so 1.40 gives -2.40.",
        "Malamud et al. (2004) maximum-likelihood fit.",
        "Size distribution","param",""};

    float_default_map["inverse_gamma_a_m2"] = 1280.0f;
    help_map["inverse_gamma_a_m2"] = {"float","1280",
        "Inverse gamma a (m2). Sets the modal landslide size.",
        "Malamud et al. (2004): 1.28e-3 km2.",
        "Size distribution","param",""};

    float_default_map["inverse_gamma_s_m2"] = -132.0f;
    help_map["inverse_gamma_s_m2"] = {"float","-132",
        "Inverse gamma s (m2), NEGATIVE. Sets the small-area rollover.",
        "Malamud et al. (2004): -1.32e-4 km2.",
        "Size distribution","param",""};

    float_default_map["pareto_alpha_pdf"] = 2.3f;
    help_map["pareto_alpha_pdf"] = {"float","2.3",
        "Truncated Pareto exponent. This is the PDF exponent, NOT cumulative.",
        "alpha_pdf = b_cumulative + 1. Valagussa et al. (2019) measure 1.75-3.11 "
        "across lithologies; Niemi's b = 1.18 is alpha_pdf = 2.18.",
        "Size distribution","param",""};

    int_default_map["min_landslide_pixels"] = 1;
    help_map["min_landslide_pixels"] = {"int","1",
        "Smallest landslide, in pixels.",
        "At coarse resolution this sets a floor above the modal size and the "
        "rollover becomes unrepresentable; the tool warns when that happens.",
        "Size distribution","param",""};

    float_default_map["max_landslide_area_m2"] = 0.0f;
    help_map["max_landslide_area_m2"] = {"float","0",
        "Largest landslide (m2). 0 = use the channel threshold, which is the "
        "physical cap.",
        "May not exceed threshold_contributing_pixels x cell area: a landslide "
        "cannot be bigger than its hillslope.",
        "Size distribution","param",""};

    string_default_map["area_volume_preset"] = "larsen_global";
    help_map["area_volume_preset"] = {"string","larsen_global",
        "V = alpha*A^gamma preset (Larsen et al. 2010 Table S1).",
        "larsen_global, larsen_soil_scar, larsen_bedrock_scar, larsen_soil_all, "
        "larsen_bedrock_all, hovius_slab, hovius_parabolic. PREFER THE SCAR "
        "PRESETS: this tool excavates scars, and the commonly quoted soil/bedrock "
        "gammas are combined scar+deposit fits. hovius_* reproduce Niemi/Yanites "
        "but overestimate volume by >= 2x per Larsen.",
        "Size distribution","param",""};

    // --- Triggering ---
    string_default_map["landslide_mode"] = "single_event";
    help_map["landslide_mode"] = {"string","single_event",
        "single_event or series.",
        "single_event triggers once at event_time_yr. series runs repeatedly "
        "between series_start_yr and series_end_yr.",
        "Triggering","param",""};

    float_default_map["event_time_yr"] = 0.0f;
    help_map["event_time_yr"] = {"float","0",
        "Trigger time for single_event mode (yr, on the user's own clock).",
        "Stamped into the time-of-failure raster, which is both the annealing "
        "clock and the cosmogenic time-since-change field.",
        "Triggering","param",""};

    float_default_map["event_volume_m3"] = 0.0f;
    help_map["event_volume_m3"] = {"float","0",
        "Total volume to remove in one event (m3).","",
        "Triggering","param",""};

    float_default_map["event_volume_per_area_m"] = 0.0f;
    help_map["event_volume_per_area_m"] = {"float","0",
        "Event volume per unit ACTIVE area (m3/m2, i.e. metres of lowering).",
        "Used when event_volume_m3 is 0. Active area is the ground that can host "
        "a landslide, not the whole raster.",
        "Triggering","param",""};

    float_default_map["series_start_yr"] = 0.0f;
    help_map["series_start_yr"] = {"float","0","Series start time (yr).","",
        "Triggering","param",""};

    float_default_map["series_end_yr"] = 100000.0f;
    help_map["series_end_yr"] = {"float","100000","Series end time (yr).","",
        "Triggering","param",""};

    float_default_map["mean_recurrence_yr"] = 1000.0f;
    help_map["mean_recurrence_yr"] = {"float","1000",
        "Mean interval between events in series mode (yr).","",
        "Triggering","param",""};

    bool_default_map["poisson_timing"] = true;
    help_map["poisson_timing"] = {"bool","true",
        "Draw inter-event intervals from an exponential distribution.",
        "Yanites et al. (2009) Eq. 14. They are explicit that this is a "
        "simplification: real landslides are triggered by storms and earthquakes "
        "and are neither independent nor uncorrelated. Set false for exactly "
        "regular events.",
        "Triggering","switch",""};

    float_default_map["long_term_erosion_rate_mm_kyr"] = 0.0f;
    help_map["long_term_erosion_rate_mm_kyr"] = {"float","0",
        "Landslide erosion rate over the active area (mm/kyr), series mode.",
        "The Niemi et al. (2005) convention: prescribe the rate and let the "
        "model find the events. Overrides event_volume_m3 in series mode.",
        "Triggering","param",""};

    float_default_map["anneal_time_yr"] = 0.0f;
    help_map["anneal_time_yr"] = {"float","0",
        "A failed pixel cannot fail again for this long (yr).",
        "Stands in for regolith regeneration. UNCALIBRATED -- there is no strong "
        "literature constraint on it.",
        "Triggering","param",""};

    float_default_map["budget_tolerance"] = 0.02f;
    help_map["budget_tolerance"] = {"float","0.02",
        "Relative tolerance when solving for the number of landslides.",
        "Reachable only when the event volume is large compared with the biggest "
        "permitted landslide; see the CLI reference.",
        "Triggering","param",""};

    int_default_map["budget_max_iterations"] = 50;
    help_map["budget_max_iterations"] = {"int","50",
        "Iteration cap for the landslide-count solve.","",
        "Triggering","param",""};

    int_default_map["random_seed"] = 1;
    help_map["random_seed"] = {"int","1",
        "Seed. The same seed reproduces a run exactly.","",
        "Triggering","param",""};

    // --- Outputs ---
    bool_default_map["write_excavation_depth_raster"] = true;
    help_map["write_excavation_depth_raster"] = {"bool","true",
        "Write cumulative excavation depth (m).",
        "Filename: <write_fname>_LS_EXCAVATION.tif",
        "Output","switch",""};

    bool_default_map["write_time_of_failure_raster"] = true;
    help_map["write_time_of_failure_raster"] = {"bool","true",
        "Write the time each pixel last failed.",
        "This is the annealing clock AND the cosmogenic time-since-change field. "
        "Filename: <write_fname>_LS_TFAIL.tif",
        "Output","switch",""};

    bool_default_map["write_post_landslide_dem"] = false;
    help_map["write_post_landslide_dem"] = {"bool","false",
        "Write the DEM after excavation.",
        "Filename: <write_fname>_LS_DEM.tif",
        "Output","switch",""};

    bool_default_map["write_landslide_csv"] = true;
    help_map["write_landslide_csv"] = {"bool","true",
        "Write one row per landslide (size in pixels, area, volume).",
        "The realised size distribution, for comparison against the requested "
        "one. Filename: <write_fname>_LS_INVENTORY.csv",
        "Output","switch",""};

    // --- Optional CRN coupling ---
    bool_default_map["calculate_crn_concentration"] = false;
    help_map["calculate_crn_concentration"] = {"bool","false",
        "Compute the post-landslide surface nuclide concentration raster.",
        "Composes existing validated CRN calls: steady state for the inherited "
        "profile, transient re-accumulation since failure. Filename: "
        "<write_fname>_LS_CONC.tif",
        "Cosmogenics","switch",""};

    string_default_map["nuclide"] = "Be10";
    help_map["nuclide"] = {"string","Be10","Nuclide: Be10, Al26 or C14.","",
        "Cosmogenics","param",""};

    string_default_map["muon_scheme"] = "BraucherBorchers";
    help_map["muon_scheme"] = {"string","BraucherBorchers",
        "Muon production scheme.","",
        "Cosmogenics","param",""};

    float_default_map["background_eff_erate_g_cm2_yr"] = 0.001f;
    help_map["background_eff_erate_g_cm2_yr"] = {"float","0.001",
        "Background effective erosion rate (g/cm2/yr), before and after failure.",
        "~3.8 mm/kyr at rho = 2650.",
        "Cosmogenics","param",""};

    float_default_map["production_scaling"] = 1.0f;
    help_map["production_scaling"] = {"float","1.0",
        "Stone 2000 production scaling (1.0 at SLHL).","",
        "Cosmogenics","param",""};

    float_default_map["total_shielding"] = 1.0f;
    help_map["total_shielding"] = {"float","1.0",
        "Topographic x snow x self shielding.",
        "Applied through the Vermeesch redistribution, separately from "
        "production_scaling.",
        "Cosmogenics","param",""};

    float_default_map["rock_density_kg_m3"] = 2650.0f;
    help_map["rock_density_kg_m3"] = {"float","2650",
        "Density of the ERODING BEDROCK COLUMN (kg/m3).",
        "Converts the background erosion rate between m/yr and g/cm2/yr. This is "
        "NOT the density of the material a landslide removes -- see "
        "scar_density_kg_m3.",
        "Cosmogenics","param",""};

    string_default_map["ice_mask_fname"] = "NULL";
    help_map["ice_mask_fname"] = {"string","NULL",
        "Raster marking ice cover (non-zero = ice). No landslide is placed on "
        "ice.",
        "Intersected with mask_fname if both are given: a pixel can host a "
        "landslide only if it is inside the active mask AND not under ice. The "
        "number of pixels excluded is reported, since it changes the active "
        "area and therefore the landslide count for a prescribed erosion rate.",
        "Placement","param",""};

    float_default_map["scar_density_kg_m3"] = -1.0f;
    help_map["scar_density_kg_m3"] = {"float","-1",
        "Density of the material a landslide EXCAVATES (kg/m3). Negative takes "
        "it from area_volume_preset.",
        "This is a different quantity from rock_density_kg_m3 and conflating "
        "them is a real error. Larsen's V = alpha A^gamma is fitted to MEASURED "
        "scar volumes, so the excavation depth is a true geometric thickness of "
        "whatever failed -- soil, at roughly 1500 kg/m3, for a shallow "
        "landslide. Carrying that at 2650 overstates the mass depth by about a "
        "factor of 1.77, and since production falls exponentially with mass "
        "depth it changes the concentration of the delivered material by more "
        "than a factor of two. The material-specific presets carry a density: "
        "larsen_soil_scar 1500, larsen_bedrock_scar 2650. larsen_global and the "
        "hovius presets are fitted across materials and carry none, so leaving "
        "this negative with one of those selected raises a warning.",
        "Cosmogenics","param",""};

    float_default_map["crn_report_time_yr"] = 0.0f;
    help_map["crn_report_time_yr"] = {"float","0",
        "Time at which concentrations are reported, on the same clock as the "
        "event times.","",
        "Cosmogenics","param",""};

    // -------------------------------------------------------------------------
    if (cli_args.help_requested)
    {
        lsdtt::lsdtt_param_parser::print_terminal_help(
            k_tool_name, k_version, k_citation, help_map, cli_args);
        return 0;
    }

    lsdtt::lsdtt_param_parser pp;
    try { pp = lsdtt::lsdtt_param_parser(cli_args.path, cli_args.fname); }
    catch (const std::exception& e)
    { std::cerr << "Error loading parameter file: " << e.what() << "\n"; return 1; }

    pp.parse_all_parameters(float_default_map, int_default_map,
                            bool_default_map, string_default_map);
    auto bf  = pp.get_bool_parameters();
    auto ff  = pp.get_float_parameters();
    auto sf  = pp.get_string_parameters();
    auto ifs = pp.get_int_parameters();

    // -------------------------------------------------------------------------
    // Cross-parameter validation, BEFORE the DEM is read or a single landslide
    // is placed. Failing here costs the user a second; failing after the event
    // series costs them the run.
    //
    // Only checked when the cosmogenic output is actually requested: muon_scheme
    // and nuclide are inert for a pure-geomorphology run, and rejecting a stale
    // line in a parameter file that changes nothing would be officious.
    // -------------------------------------------------------------------------
    if (bf["calculate_crn_concentration"] &&
        !landslide_cosmo_parameters_supported(sf["muon_scheme"], sf["nuclide"]))
    {
        return 1;
    }

    // -------------------------------------------------------------------------
    // Load and prepare the DEM
    // -------------------------------------------------------------------------
    const std::string dem_path = pp.get_full_read_path();
    std::cout << "Opening DEM: " << dem_path << "\n";
    lsdtt::lsdtt_data_loader loader;
    loader.open_raster(dem_path);
    lsdtt::lsdtt_raster_cpu dem(loader);

    // Projected CRS, checked HERE with a clean message. The check exists deeper
    // in the stack too, but it throws, and an uncaught throw from a CLI is an
    // abort with a core dump rather than a diagnosis.
    try { dem.get_meta().check_projected_metres(); }
    catch (const std::exception& e)
    {
        std::cerr << "FATAL: " << e.what() << "\n"
                  << "       Landslide areas and volumes are computed from cell "
                     "size in metres, so a\n"
                  << "       geographic (degree) grid would produce meaningless "
                     "numbers. Reproject\n"
                  << "       the DEM to UTM or another projected CRS first "
                     "(lsdtt-raster-preprocessing).\n";
        return 1;
    }

    const lsdtt::lsdtt_raster_meta meta = dem.get_meta();
    const double cell = meta.cell_size;
    std::cout << "DEM: " << meta.n_rows << " x " << meta.n_cols
              << " at " << cell << " m\n";

    lsdtt::lsdtt_raster_cpu filled = bf["fill_before_routing"]
        ? lsdtt::fill_pits(dem, ff["min_slope_for_fill"])
        : dem;

    lsdtt::lsdtt_flowinfo_options fopts;
    fopts.boundary_conditions = sf["boundary_conditions"];
    std::unique_ptr<lsdtt::lsdtt_flow_base> flow =
        lsdtt::make_flow_base(filled.get_data_array(), filled.get_mask_array(),
                              meta, fopts);
    if (!flow) { std::cerr << "FATAL: flow routing failed.\n"; return 1; }

    // -------------------------------------------------------------------------
    // Size distribution
    // -------------------------------------------------------------------------
    lsdtt::LandslideSizeParams sp;
    if (!lsdtt::parse_landslide_size_model(sf["landslide_size_model"], sp.model))
    {
        std::cerr << "FATAL: unknown landslide_size_model \""
                  << sf["landslide_size_model"]
                  << "\". Options: inverse_gamma, truncated_pareto.\n";
        return 1;
    }
    sp.rho              = ff["inverse_gamma_rho"];
    sp.a_param          = ff["inverse_gamma_a_m2"];
    sp.s_param          = ff["inverse_gamma_s_m2"];
    sp.pareto_alpha_pdf = ff["pareto_alpha_pdf"];
    sp.n_min_pixels     = ifs["min_landslide_pixels"];

    if (!lsdtt::AreaVolumeScaling::from_name(sf["area_volume_preset"], sp.av))
    {
        std::cerr << "FATAL: unknown area_volume_preset \""
                  << sf["area_volume_preset"] << "\".\n";
        return 1;
    }

    const std::int64_t chan_thresh = ifs["threshold_contributing_pixels"];
    const double cap_area = static_cast<double>(chan_thresh) * cell * cell;
    sp.area_max_m2 = (ff["max_landslide_area_m2"] > 0.0f)
        ? static_cast<double>(ff["max_landslide_area_m2"]) : cap_area;

    // A landslide may not be larger than the hillslope it comes from. Fatal
    // rather than clamped: silently clamping would misreport the distribution
    // actually sampled.
    if (sp.area_max_m2 > cap_area)
    {
        std::cerr << "FATAL: max_landslide_area_m2 (" << sp.area_max_m2
                  << " m2) exceeds the channel threshold ("
                  << chan_thresh << " px = " << cap_area << " m2).\n"
                  << "       A landslide cannot be larger than the hillslope it "
                     "comes from. Lower max_landslide_area_m2 or raise "
                     "threshold_contributing_pixels.\n";
        return 1;
    }

    std::unique_ptr<lsdtt::LandslideSizeSampler> sampler;
    try { sampler = std::make_unique<lsdtt::LandslideSizeSampler>(sp, cell); }
    catch (const std::exception& e)
    { std::cerr << "FATAL: " << e.what() << "\n"; return 1; }

    std::cout << "\nSize model    : " << lsdtt::to_string(sp.model) << "\n"
              << "A_max         : " << sp.area_max_m2 << " m2 ("
              << sampler->n_max_pixels() << " px)\n"
              << "E[V]          : " << sampler->mean_volume_m3() << " m3\n"
              << "Truncated at A_max: "
              << 100.0 * sampler->truncated_fraction_number() << "% by number, "
              << 100.0 * sampler->truncated_fraction_volume() << "% by volume\n";
    if (sampler->rollover_unresolved())
        std::cout << "WARNING: cell size cannot resolve the small-area "
                     "rollover.\n"
                  << "         Every sampled landslide sits on the power-law "
                     "tail, so the choice\n"
                  << "         of distribution family has no effect at this "
                     "resolution.\n";

    // -------------------------------------------------------------------------
    // Placement
    // -------------------------------------------------------------------------
    lsdtt::lsdtt_raster_cpu mask_raster;
    const lsdtt::lsdtt_raster_cpu* mask_ptr = nullptr;
    if (sf["mask_fname"] != "NULL")
    {
        try
        {
            lsdtt::lsdtt_data_loader mask_loader;
            mask_loader.open_raster(pp.add_read_path_if_needed(sf["mask_fname"]));
            mask_raster = lsdtt::lsdtt_raster_cpu(mask_loader);
            mask_ptr = &mask_raster;
            std::cout << "Mask          : " << sf["mask_fname"] << "\n";
        }
        catch (const std::exception& e)
        { std::cerr << "FATAL: could not load mask: " << e.what() << "\n"; return 1; }
    }

    // Ice excludes ground from landsliding. The placer takes ONE mask, so the
    // ice is intersected into it here rather than threaded through as a second
    // field -- a pixel is hostable if it is inside the active mask AND not
    // under ice.
    lsdtt::lsdtt_raster_cpu combined_mask;
    if (sf["ice_mask_fname"] != "NULL")
    {
        lsdtt::lsdtt_raster_cpu ice_raster;
        try
        {
            lsdtt::lsdtt_data_loader ice_loader;
            ice_loader.open_raster(
                pp.add_read_path_if_needed(sf["ice_mask_fname"]));
            ice_raster = lsdtt::lsdtt_raster_cpu(ice_loader);
        }
        catch (const std::exception& e)
        {
            std::cerr << "FATAL: could not load ice mask: " << e.what() << "\n";
            return 1;
        }
        if (ice_raster.get_meta().n_rows != meta.n_rows ||
            ice_raster.get_meta().n_cols != meta.n_cols)
        {
            std::cerr << "FATAL: the ice mask does not match the DEM grid.\n";
            return 1;
        }

        const auto& icd = ice_raster.get_data_array();
        const auto& icm = ice_raster.get_mask_array();
        const float ice_ndv = static_cast<float>(ice_raster.get_meta().no_data_value);

        lsdtt::lsdtt_array2d_f  cd(meta.n_rows, meta.n_cols, 0.0f);
        lsdtt::lsdtt_array2d_u8 ck(meta.n_rows, meta.n_cols, 1);
        std::int64_t n_excluded = 0;
        for (int r = 0; r < meta.n_rows; ++r)
        {
            for (int c = 0; c < meta.n_cols; ++c)
            {
                const bool in_active = mask_ptr
                    ? (mask_ptr->get_mask_array()(r, c) &&
                       mask_ptr->get_data_array()(r, c) != 0.0f)
                    : true;
                const bool iced = icm(r, c) && icd(r, c) != ice_ndv
                                            && icd(r, c) != 0.0f;
                if (iced && in_active) ++n_excluded;
                cd(r, c) = (in_active && !iced) ? 1.0f : 0.0f;
            }
        }
        combined_mask = lsdtt::lsdtt_raster_cpu(std::move(cd), std::move(ck), meta);
        mask_ptr = &combined_mask;
        std::cout << "Ice mask      : " << sf["ice_mask_fname"] << " -- "
                  << n_excluded << " pixel(s) excluded from landsliding\n";
    }

    lsdtt::LandslidePlacementConfig pcfg;
    pcfg.channel_threshold_pixels = chan_thresh;
    lsdtt::LandslidePlacer placer(*flow, meta, pcfg, mask_ptr);

    const double active_area =
        static_cast<double>(placer.n_eligible_hosts(1)) * cell * cell;
    std::cout << "Active area   : " << active_area << " m2 ("
              << placer.n_eligible_hosts(1) << " px)\n"
              << "Largest hostable landslide: " << placer.max_hostable_pixels()
              << " px\n\n";
    if (placer.n_eligible_hosts(1) == 0)
    {
        std::cerr << "FATAL: no pixel can host a landslide. Every cell is at or "
                     "above the channel threshold, or the mask excludes them all.\n";
        return 1;
    }

    // -------------------------------------------------------------------------
    // Run
    // -------------------------------------------------------------------------
    lsdtt::lsdtt_raster_cpu tfail = like(filled, static_cast<float>(meta.no_data_value));
    lsdtt::lsdtt_raster_cpu exc   = like(filled, 0.0f);
    lsdtt::lsdtt_raster_cpu work  = filled;

    std::vector<int> realised;
    const std::string mode = sf["landslide_mode"];

    if (mode == "single_event")
    {
        lsdtt::LandslideEventConfig ecfg;
        ecfg.event_time_yr            = ff["event_time_yr"];
        ecfg.target_volume_m3         = ff["event_volume_m3"];
        ecfg.volume_per_area_m3_per_m2 = ff["event_volume_per_area_m"];
        ecfg.anneal_time_yr           = ff["anneal_time_yr"];
        ecfg.budget_tol               = ff["budget_tolerance"];
        ecfg.budget_max_iter          = ifs["budget_max_iterations"];
        ecfg.seed = static_cast<std::uint64_t>(ifs["random_seed"]);

        if (!(ecfg.target_volume_m3 > 0.0) &&
            !(ecfg.volume_per_area_m3_per_m2 > 0.0))
        {
            std::cerr << "FATAL: single_event mode needs event_volume_m3 or "
                         "event_volume_per_area_m.\n";
            return 1;
        }

        const lsdtt::LandslideEventResult r = lsdtt::trigger_landslide_event(
            work, tfail, exc, placer, *sampler, ecfg);
        std::cout << r.summary() << "\n";
        realised = r.realised_pixels;
    }
    else if (mode == "series")
    {
        lsdtt::LandslideSeriesConfig scfg;
        scfg.start_time_yr      = ff["series_start_yr"];
        scfg.end_time_yr        = ff["series_end_yr"];
        scfg.mean_recurrence_yr = ff["mean_recurrence_yr"];
        scfg.poisson_timing     = bf["poisson_timing"];
        scfg.target_volume_per_event_m3   = ff["event_volume_m3"];
        scfg.long_term_erosion_rate_mm_kyr = ff["long_term_erosion_rate_mm_kyr"];
        scfg.anneal_time_yr     = ff["anneal_time_yr"];
        scfg.budget_tol         = ff["budget_tolerance"];
        scfg.budget_max_iter    = ifs["budget_max_iterations"];
        scfg.seed = static_cast<std::uint64_t>(ifs["random_seed"]);

        if (!(scfg.target_volume_per_event_m3 > 0.0) &&
            !(scfg.long_term_erosion_rate_mm_kyr > 0.0))
        {
            std::cerr << "FATAL: series mode needs long_term_erosion_rate_mm_kyr "
                         "or event_volume_m3.\n";
            return 1;
        }

        const lsdtt::LandslideSeriesResult r = lsdtt::run_landslide_series(
            work, tfail, exc, placer, *sampler, scfg);
        std::cout << r.summary() << "\n";
        realised = r.realised_pixels;
    }
    else
    {
        std::cerr << "FATAL: unknown landslide_mode \"" << mode
                  << "\". Options: single_event, series.\n";
        return 1;
    }

    if (realised.empty())
        std::cerr << "\nWarning: no landslides were placed. Check the budget, "
                     "the channel threshold and the mask.\n";

    // -------------------------------------------------------------------------
    // Outputs
    // -------------------------------------------------------------------------
    lsdtt::raster_output_format out_fmt = lsdtt::raster_output_format::cog;
    if (sf["write_output_format"] == "gtiff")
        out_fmt = lsdtt::raster_output_format::gtiff;
    else if (sf["write_output_format"] == "same_as_input")
        out_fmt = lsdtt::raster_output_format::same_as_input;

    const std::string base = pp.get_full_write_path();

    if (bf["write_excavation_depth_raster"])
    {
        const std::string path = base + "_LS_EXCAVATION.tif";
        loader.write_full(exc.get_data(), exc.get_meta(), path, out_fmt);
        std::cout << "Wrote " << path << "\n";
    }
    if (bf["write_time_of_failure_raster"])
    {
        const std::string path = base + "_LS_TFAIL.tif";
        loader.write_full(tfail.get_data(), tfail.get_meta(), path, out_fmt);
        std::cout << "Wrote " << path << "\n";
    }
    if (bf["write_post_landslide_dem"])
    {
        const std::string path = base + "_LS_DEM.tif";
        loader.write_full(work.get_data(), work.get_meta(), path, out_fmt);
        std::cout << "Wrote " << path << "\n";
    }
    if (bf["write_landslide_csv"])
    {
        const std::string path = base + "_LS_INVENTORY.csv";
        std::ofstream f(path);
        if (!f.is_open())
            std::cerr << "Warning: could not write " << path << "\n";
        else
        {
            f << "landslide_id,n_pixels,area_m2,volume_m3,mean_depth_m\n";
            f.precision(8);
            for (std::size_t i = 0; i < realised.size(); ++i)
            {
                const double a = static_cast<double>(realised[i]) * cell * cell;
                const double v = sampler->volume_of_pixels(realised[i]);
                f << (i + 1) << "," << realised[i] << "," << a << "," << v
                  << "," << (a > 0.0 ? v / a : 0.0) << "\n";
            }
            std::cout << "Wrote " << path << " (" << realised.size()
                      << " landslides)\n";
        }
    }

    if (bf["calculate_crn_concentration"])
    {
        lsdtt::LandslideCrnConfig ccfg;
        ccfg.background_eff_erate_g_cm2_yr = ff["background_eff_erate_g_cm2_yr"];
        ccfg.production_scaling = ff["production_scaling"];
        ccfg.total_shielding    = ff["total_shielding"];
        ccfg.muon_scheme        = sf["muon_scheme"];
        ccfg.rock_density_kg_m3 = ff["rock_density_kg_m3"];

        // The scar density is a property of the material that failed, so it
        // comes from the area-volume preset unless the user overrides it. A
        // preset fitted across materials has no defensible value and says so.
        if (ff["scar_density_kg_m3"] > 0.0f)
        {
            ccfg.scar_density_kg_m3 = ff["scar_density_kg_m3"];
        }
        else
        {
            ccfg.scar_density_kg_m3 = sp.av.scar_density_kg_m3;
            const std::string note = lsdtt::area_volume_density_note(sp.av);
            if (!note.empty())
                std::cout << "\n  WARNING: " << note << "\n\n";
        }
        std::cout << "  Densities: bedrock column " << ccfg.rock_density_kg_m3
                  << " kg/m3, scar material " << ccfg.scar_density_kg_m3
                  << " kg/m3\n";
        ccfg.t_now_yr           = ff["crn_report_time_yr"];
        ccfg.nuclide            = sf["nuclide"];

        try
        {
            const lsdtt::lsdtt_raster_cpu conc =
                lsdtt::compute_landslide_surface_concentration(exc, tfail, ccfg);
            const std::string path = base + "_LS_CONC.tif";
            loader.write_full(conc.get_data(), conc.get_meta(), path, out_fmt);
            std::cout << "Wrote " << path << "\n"
                      << "  background (no landslide) concentration: "
                      << lsdtt::landslide_crn_background_concentration(ccfg)
                      << " atoms/g\n";
        }
        catch (const std::exception& e)
        {
            // FATAL, not a warning. The user asked for a concentration raster;
            // if it could not be produced, the run did not do what it was told.
            // This used to warn and fall through to exit 0, so a scripted
            // workflow saw a success that silently produced nothing. The
            // up-front check above catches the two known causes early and with
            // a better message; this is the backstop for the rest.
            std::cerr << "FATAL: CRN concentration failed: " << e.what() << "\n";
            return 1;
        }
    }

    pp.print_help(help_map, pp.get_write_fname() + "_README",
                  k_version, k_citation);
    return 0;
}
