#pragma once
#include "lsdtt_crn_parameters.hpp"
#include "lsdtt_crn_particle.hpp"
#include <string>
#include <vector>

namespace lsdtt
{

// ─────────────────────────────────────────────────────────────────────────────
// CrnColumnSample
// One output row from a column simulation.
// ─────────────────────────────────────────────────────────────────────────────
struct CrnColumnSample
{
    double depth_m;             // current depth below surface (m)
    double effective_depth;     // current effective depth (g/cm²)
    double conc_10Be;           // atoms/g
    double conc_26Al;           // atoms/g
    double conc_14C;            // atoms/g
};

// ─────────────────────────────────────────────────────────────────────────────
// lsdtt_crn_column
//
// Simulates cosmogenic nuclide concentrations in a bedrock/sediment column
// undergoing a user-specified erosion history.  Mirrors the
// transient_column_calculator functionality from LSDTopoTools2.
//
// Erosion history file format (CSV, three columns with header):
//   time,rate,removal
//   - time   : time before present (yr) of the rate change; increasing order
//   - rate   : erosion rate (m/yr); 0 → instantaneous removal event
//   - removal: depth removed (m), used only when rate == 0
//   - last row: time = -9999, rate = long-term background rate, removal = 0
//
// Coordinate convention:
//   - sampling_depths_m[i] is the depth BELOW THE CURRENT SURFACE in metres
//   - rock_density_kg_m3 converts depths to effective depths (g/cm²):
//       d_eff = depth_m × density_kg_m3 × 0.1
// ─────────────────────────────────────────────────────────────────────────────
class lsdtt_crn_column
{
public:
    explicit lsdtt_crn_column(double rock_density_kg_m3 = 2650.0)
        : rock_density_(rock_density_kg_m3) {}

    // ── Erosion history ───────────────────────────────────────────────────────
    void read_erosion_time_series(const std::string& csv_path);
    void set_erosion_time_series(const std::vector<double>& time_yr,
                                 const std::vector<double>& rate_m_yr,
                                 const std::vector<double>& removal_m);

    // Total depth (m) added to sampling positions to recover their starting
    // depth at the beginning of the erosion record.
    double calculate_depth_change() const;

    // ── Column simulation ─────────────────────────────────────────────────────
    // Run the full transient simulation for a set of sampling depths.
    // production_scaling : Stone 2000 scaling factor (1.0 at SLHL)
    // muon_scheme        : parameter scheme name (see lsdtt_crn_parameters)
    // initial_10Be / 26Al / 14C : optional uniform initial concentrations
    //   (atoms/g); NaN or < 0 → initialise at SS with background erosion rate
    std::vector<CrnColumnSample>
    run_transient(const std::vector<double>& sampling_depths_m,
                  double production_scaling,
                  const std::string& muon_scheme,
                  double initial_10Be = -1.0,
                  double initial_26Al = -1.0,
                  double initial_14C  = -1.0) const;

    // ── Simple single-particle modes ──────────────────────────────────────────
    // Steady-state concentration at a point sample (effective depth d_eff g/cm²).
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

    // ── Accessors ─────────────────────────────────────────────────────────────
    double rock_density() const { return rock_density_; }
    const std::vector<double>& time_vec()    const { return time_vec_; }
    const std::vector<double>& rate_vec()    const { return rate_vec_; }
    const std::vector<double>& removal_vec() const { return removal_vec_; }

private:
    double rock_density_;               // kg/m³
    std::vector<double> time_vec_;      // yr BP, ascending; last entry = -9999
    std::vector<double> rate_vec_;      // m/yr
    std::vector<double> removal_vec_;   // m  (non-zero only when rate==0)

    void check_erosion_time_series() const;
    lsdtt_crn_parameters make_params(double production_scaling,
                                     const std::string& muon_scheme) const;
    double to_eff(double depth_m) const { return depth_m * rock_density_ * 0.1; }
};

} // namespace lsdtt
