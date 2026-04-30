#pragma once
#include <string>
#include <map>
#include <vector>

namespace lsdtt
{

// ─────────────────────────────────────────────────────────────────────────────
// lsdtt_crn_parameters
//
// Holds all physics constants for cosmogenic nuclide (CRN) production and decay.
//
// Workflow per pixel:
//   1. reset_scaling(muon_scheme)   — restore standard F values for the scheme
//   2. scale_F_values(shielding, nuclide) — fold total shielding into F arrays
//   3. Pass const& to lsdtt_crn_particle for concentration calculation
//
// Supported nuclides: Be10, Al26, C14
// Supported muon schemes: "Braucher", "BraucherBorchers" (default),
//                         "Granger", "Schaller", "newCRONUS"
// ─────────────────────────────────────────────────────────────────────────────
class lsdtt_crn_parameters
{
public:
    lsdtt_crn_parameters();

    // ── Parameter schemes ────────────────────────────────────────────────────
    void set_braucher_parameters();
    void set_braucher_borchers_parameters();   // default
    void set_granger_parameters();
    void set_schaller_parameters();
    void set_new_cronus_parameters();

    // Selects scheme by name (case-sensitive: "Braucher", "BraucherBorchers",
    // "Granger", "Schaller", "newCRONUS").  Falls back to BraucherBorchers.
    void reset_scaling(const std::string& muon_scheme);

    // ── Shielding ─────────────────────────────────────────────────────────────
    // Rescales F values for the named nuclide using a Newton-Raphson search that
    // converts total_shielding (= production_scale × topo_shield, dimensionless
    // [0,1]) into an equivalent effective depth so that each production pathway
    // is attenuated correctly.  Call reset_scaling first, then scale_F_values.
    // nuclide: "Be10", "Al26", or "C14".
    void scale_F_values(double total_shielding, const std::string& nuclide);

    // ── Atmospheric pressure ─────────────────────────────────────────────────
    // Standard-atmosphere formula (Balco's barometric equation with SLHL values).
    // Returns pressure in mbar.
    static double atmospheric_pressure_standard(double elevation_m);

    // ── Stone 2000 scaling ───────────────────────────────────────────────────
    // Returns the dimensionless Stone (2000) production-rate scaling factor.
    // lat_deg: site latitude (deg, southern hemisphere OK — sign is ignored).
    // P_mbar : atmospheric pressure in mbar.
    // Fsp    : fraction of total production from spallation (use 1.0 for default).
    double stone2000sp(double lat_deg, double P_mbar, double Fsp) const;

    // ── Standardisation ──────────────────────────────────────────────────────
    // Returns the factor to multiply a measured concentration by to convert it
    // to the reference standard (07KNSTD for Be10, KNSTD for Al26).
    double be10_standardisation_factor(const std::string& standard_name) const;
    double al26_standardisation_factor(const std::string& standard_name) const;

    // Applies standardisation to all concentrations in a vector.
    // standard_names: one per sample (empty string → reference standard assumed).
    // nuclide: "Be10" or "Al26".
    std::vector<double> apply_standardisation(
        const std::vector<double>& concentrations,
        const std::vector<std::string>& standard_names,
        const std::string& nuclide) const;

    // ── Getters (physics constants) ───────────────────────────────────────────
    double get_lambda_10Be() const { return lambda_10Be_; }
    double get_lambda_26Al() const { return lambda_26Al_; }
    double get_lambda_14C()  const { return lambda_14C_;  }
    double get_P0_10Be()     const { return P0_10Be_; }
    double get_P0_26Al()     const { return P0_26Al_; }
    double get_P0_14C()      const { return P0_14C_;  }
    double get_Gamma(int i)  const { return Gamma_[i]; }
    double get_F_10Be(int i) const { return F_10Be_[i]; }
    double get_F_26Al(int i) const { return F_26Al_[i]; }
    double get_F_14C(int i)  const { return F_14C_[i]; }

private:
    // Decay constants (yr⁻¹)
    double lambda_10Be_, lambda_26Al_, lambda_14C_, lambda_36Cl_;

    // Sea-level high-latitude production rates (atoms g⁻¹ yr⁻¹)
    double P0_10Be_, P0_26Al_, P0_14C_, P0_36Cl_;

    // Attenuation lengths (g cm⁻²):
    //   [0] spallation  [1] slow neg. muons  [2] fast muons  [3] very fast muons
    double Gamma_[4];

    // Production fractions (dimensionless); modified in-place by scale_F_values
    double F_10Be_[4], F_26Al_[4], F_14C_[4];

    // Initialise to Vermeesh (2007) defaults
    void init_defaults();
    void init_standardisation_maps();

    // Internal Newton-Raphson: finds effective depth x such that
    //   Σᵢ F[i] * exp(-x/Γ[i]) == total_shielding, then applies it to F[].
    static void scale_one_nuclide(double total_shielding,
                                  const double Gamma[4],
                                  double F[4]);

    // Linear 1-D interpolation on ordered x
    static double interp1d(const std::vector<double>& x,
                           const std::vector<double>& y,
                           double xi);

    std::map<std::string, double> be10_std_map_;
    std::map<std::string, double> al26_std_map_;
};

} // namespace lsdtt
