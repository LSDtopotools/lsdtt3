#pragma once
// =============================================================================
//  DO NOT CHANGE THE NUMBERS THIS CODE PRODUCES — ask Simon Mudd first.
//
//  The cosmogenic calculations in this file (the CAIRN erosion-rate method and
//  the step-change concentration solutions) were checked by hand over several
//  years. To stop anyone changing those numbers by accident, the test suite
//  keeps a stored list of known-correct output values in
//  tests/expected/cosmo/numerics_sentinel.csv. The tests in
//  tests/cpp/test_cosmo_physics.cpp re-run these calculations and compare the
//  results against that stored list to 12 decimal places.
//
//  So if you edit this file and any output value changes — even from a refactor
//  or an algebra rewrite you believe is harmless — those tests will fail. When
//  that happens, do NOT edit the stored list to make the tests pass again. A
//  changed value has to be reviewed and approved by Simon Mudd before it is
//  committed.
// =============================================================================
#include "lsdtt_crn_parameters.hpp"
#include "lsdtt_crn_particle.hpp"
#include <string>
#include <vector>

namespace lsdtt
{

// 
// CrnColumnSample
// One output row from a column simulation.
// 
struct CrnColumnSample
{
    double depth_m;             // current depth below surface (m)
    double effective_depth;     // current effective depth (g/cm^2)
    double conc_10Be;           // atoms/g
    double conc_26Al;           // atoms/g
    double conc_14C;            // atoms/g
};

// 
// lsdtt_crn_column
//
// Simulates cosmogenic nuclide concentrations in a bedrock/sediment column
// undergoing a user-specified erosion history.  Mirrors the
// transient_column_calculator functionality from LSDTopoTools2.
//
// Erosion history file format (CSV, three columns with header):
//   time,rate,removal
//   - time   : time before present (yr) of the rate change; increasing order
//   - rate   : erosion rate (m/yr); 0 -> instantaneous removal event
//   - removal: depth removed (m), used only when rate == 0
//   - last row: time = -9999, rate = long-term background rate, removal = 0
//
// Coordinate convention:
//   - sampling_depths_m[i] is the depth BELOW THE CURRENT SURFACE in metres
//   - rock_density_kg_m3 converts depths to effective depths (g/cm^2):
//       d_eff = depth_m x density_kg_m3 x 0.1
// 
class lsdtt_crn_column
{
public:
    explicit lsdtt_crn_column(double rock_density_kg_m3 = 2650.0)
        : rock_density_(rock_density_kg_m3) {}

    //  Erosion history 
    void read_erosion_time_series(const std::string& csv_path);
    void set_erosion_time_series(const std::vector<double>& time_yr,
                                 const std::vector<double>& rate_m_yr,
                                 const std::vector<double>& removal_m);

    // Total depth (m) added to sampling positions to recover their starting
    // depth at the beginning of the erosion record.
    double calculate_depth_change() const;

    //  Column simulation 
    // Run the full transient simulation for a set of sampling depths.
    // production_scaling : Stone 2000 scaling factor (1.0 at SLHL)
    // muon_scheme        : parameter scheme name (see lsdtt_crn_parameters)
    // initial_10Be / 26Al / 14C : optional uniform initial concentrations
    //   (atoms/g); NaN or < 0 -> initialise at SS with background erosion rate
    std::vector<CrnColumnSample>
    run_transient(const std::vector<double>& sampling_depths_m,
                  double production_scaling,
                  const std::string& muon_scheme,
                  double initial_10Be = -1.0,
                  double initial_26Al = -1.0,
                  double initial_14C  = -1.0) const;

    //  Simple single-particle modes 
    // Steady-state concentration at a point sample (effective depth d_eff g/cm^2).
    CrnColumnSample steady_state(double effective_depth_g_cm2,
                                 double eff_erate_g_cm2_yr,
                                 double production_scaling,
                                 const std::string& muon_scheme) const;

    // Step change: particle was at eff_erate_old for a very long time, then
    // at time_since_step_yr before present changed to eff_erate_new.
    // Returns concentration at depth effective_depth_g_cm2 below the current surface.
    CrnColumnSample step_change(double effective_depth_g_cm2,
                                double eff_erate_old_g_cm2_yr,
                                double eff_erate_new_g_cm2_yr,
                                double time_since_step_yr,
                                double production_scaling,
                                const std::string& muon_scheme) const;

    //  Point modes with SHIELDING HANDLED SEPARATELY
    // ADDITIVE overloads (2026-07-27). The overloads above take a single
    // scaling factor and apply all of it through apply_production_scaling, i.e.
    // a flat P0 multiply. That is right for a SITE PRODUCTION SCALING (it may
    // exceed 1 at altitude) but WRONG for SHIELDING, which the rest of lsdtt3
    // -- calculate_crn_concentration_raster and every validated comparison --
    // puts through the Vermeesch virtual-attenuation redistribution in
    // scale_F_values, attenuating spallation and the muon pathways unequally.
    //
    // Passing production_scaling * total_shielding to the single-factor form
    // therefore gives a DIFFERENT answer from the raster pipeline, and made
    // lsdtt-cosmo-column-tool's two point modes disagree by 11.4 % at a
    // shielding of 0.5 -- in a case (t = 0) where a limit identity says they
    // must agree exactly.
    //
    // These overloads take the two separately and apply each by its own
    // mechanism. total_shielding == 1.0 reproduces the single-factor overloads
    // bit-for-bit, so this is purely additive: nothing that calls the existing
    // signatures can move, and the numerics sentinels are unaffected.
    //
    // total_shielding: topographic x snow x self, dimensionless, <= 1.
    CrnColumnSample steady_state(double effective_depth_g_cm2,
                                 double eff_erate_g_cm2_yr,
                                 double production_scaling,
                                 double total_shielding,
                                 const std::string& muon_scheme) const;

    CrnColumnSample step_change(double effective_depth_g_cm2,
                                double eff_erate_old_g_cm2_yr,
                                double eff_erate_new_g_cm2_yr,
                                double time_since_step_yr,
                                double production_scaling,
                                double total_shielding,
                                const std::string& muon_scheme) const;

    // run_transient on the raster convention too, so an arbitrary erosion
    // history can drive a pixel prediction and agree with the raster kernel.
    // The single-factor overload above is unchanged and still applies its
    // argument as a flat P0 multiply.
    std::vector<CrnColumnSample>
    run_transient(const std::vector<double>& sampling_depths_m,
                  double production_scaling,
                  double total_shielding,
                  const std::string& muon_scheme,
                  double initial_10Be = -1.0,
                  double initial_26Al = -1.0,
                  double initial_14C  = -1.0) const;

    //  Accessors 
    double rock_density() const { return rock_density_; }
    const std::vector<double>& time_vec()    const { return time_vec_; }
    const std::vector<double>& rate_vec()    const { return rate_vec_; }
    const std::vector<double>& removal_vec() const { return removal_vec_; }

private:
    double rock_density_;               // kg/m^3
    std::vector<double> time_vec_;      // yr BP, ascending; last entry = -9999
    std::vector<double> rate_vec_;      // m/yr
    std::vector<double> removal_vec_;   // m  (non-zero only when rate==0)

    void check_erosion_time_series() const;
    lsdtt_crn_parameters make_params(double production_scaling,
                                     const std::string& muon_scheme) const;
    // make_params plus the Vermeesch shielding redistribution. Guarded so that
    // total_shielding == 1.0 returns make_params' object bit-for-bit.
    lsdtt_crn_parameters make_params_shielded(double production_scaling,
                                              double total_shielding,
                                              const std::string& muon_scheme) const;
    // Shared body of the two run_transient overloads; they differ only in the
    // parameter object, so the evolution lives in exactly one place.
    std::vector<CrnColumnSample>
    run_transient_impl(const std::vector<double>& sampling_depths_m,
                       const lsdtt_crn_parameters& params,
                       double initial_10Be,
                       double initial_26Al,
                       double initial_14C) const;
    double to_eff(double depth_m) const { return depth_m * rock_density_ * 0.1; }
};

} // namespace lsdtt
