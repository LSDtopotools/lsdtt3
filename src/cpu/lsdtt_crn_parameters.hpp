#pragma once
#include <string>
#include <map>
#include <vector>

namespace lsdtt
{

// ─────────────────────────────────────────────────────────────────────────────
// CrnNcData
//
// ERA40 atmospheric grids and Muscheler VDM record loaded from crn_parameters.nc.
// Load once with lsdtt_crn_parameters::load_nc_data(), then pass by const-ref
// to ERA40 and paleomagnetic correction functions.
// ─────────────────────────────────────────────────────────────────────────────
struct CrnNcData
{
    // ERA40 grid axes
    std::vector<double> era40_lat;     // (73,)   deg, 90 → -90 (descending)
    std::vector<double> era40_lon;     // (145,)  deg, 0 → 360 (ascending)
    // ERA40 climatological fields — row-major [ilat][ilon]
    std::vector<float>  era40_meanP;   // (73×145) hPa
    std::vector<float>  era40_meanT;   // (73×145) K
    // Lifton COSPAR CIRA-86 lapse rate polynomial coefficients (7 terms, K m⁻¹)
    std::vector<double> lifton_lr;     // (7,)
    // Muscheler VDM record
    std::vector<double> muscheler_time; // (5970,) yr BP, ascending
    std::vector<float>  muscheler_vdm;  // (5970,) 10²² Am²
    bool loaded = false;
};

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

    // ERA40-based atmospheric pressure (Balco/Lifton, CRONUS origin).
    // Bilinear interpolation of ERA40 climatological SLP and SLT grids, then
    // applies the Lifton latitude-dependent lapse rate polynomial.
    // nc: loaded via load_nc_data(). Falls back to atmospheric_pressure_standard
    //     if nc.loaded == false.
    // Returns pressure in mbar.
    static double atmospheric_pressure_era40(double lat_deg, double lon_deg,
                                              double elevation_m,
                                              const CrnNcData& nc);

    // Load ERA40 and paleomagnetic data from a crn_parameters.nc file.
    // Returns CrnNcData with loaded=false if the file cannot be read or if
    // the code was compiled without netCDF-C support.
    static CrnNcData load_nc_data(const std::string& nc_path);

    // ── Stone 2000 scaling ───────────────────────────────────────────────────
    // Returns the dimensionless Stone (2000) production-rate scaling factor.
    // lat_deg: site latitude (deg, southern hemisphere OK — sign is ignored).
    // P_mbar : atmospheric pressure in mbar.
    // Fsp    : fraction of total production from spallation (use 1.0 for default).
    double stone2000sp(double lat_deg, double P_mbar, double Fsp) const;

    // Stone 2000 scaling with time-integrated paleomagnetic correction
    // (Stone/Lal scheme, following Martin et al. 2017 / CREp).
    // Computes the effective time-averaged spallation scaling factor over the
    // sample's exposure age using the Muscheler VDM record, then rescales the
    // total Stone factor accordingly.  Returns stone2000sp(lat,P,Fsp) if
    // nc.loaded == false or age_yr <= 0.
    // age_yr: sample exposure age in years (positive, years before present).
    double stone2000sp_paleomag(double lat_deg, double P_mbar, double Fsp,
                                double age_yr, const CrnNcData& nc) const;

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

    // Bilinear interpolation on ERA40 grid (lon ascending, lat descending).
    // Returns NaN if the grid vectors are empty.
    static double bilinear_era40(const std::vector<double>& lat_ax,
                                 const std::vector<double>& lon_ax,
                                 const std::vector<float>&  grid,
                                 double lat, double lon);

    // Stone spallation factor as a function of geomagnetic cutoff rigidity Rc [GV].
    // Builds the Rc→sf look-up table at pressure P_mbar using the Stone (2000)
    // coefficients and linearly interpolates.  Mirrors Rc_StoneFactors() in
    // PR_processes.py (Martin et al. / CREp).
    static double stone_factor_at_rc(double Rc_GV, double P_mbar);

    std::map<std::string, double> be10_std_map_;
    std::map<std::string, double> al26_std_map_;
};

} // namespace lsdtt
