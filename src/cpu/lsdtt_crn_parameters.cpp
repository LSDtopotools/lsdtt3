#include "lsdtt_crn_parameters.hpp"
#include <cmath>
#include <iostream>
#include <stdexcept>

namespace lsdtt
{

// ─────────────────────────────────────────────────────────────────────────────
lsdtt_crn_parameters::lsdtt_crn_parameters()
{
    init_defaults();
    set_braucher_borchers_parameters();
    init_standardisation_maps();
}

// ─────────────────────────────────────────────────────────────────────────────
void lsdtt_crn_parameters::init_defaults()
{
    // Vermeesh (2007) default values — overwritten by any set_*_parameters call
    lambda_10Be_ = 456e-9;
    lambda_26Al_ = 980e-9;
    lambda_14C_  = 121e-6;
    lambda_36Cl_ = 230e-8;

    P0_10Be_ = 5.11;
    P0_26Al_ = 30.31;
    P0_14C_  = 5.86;
    P0_36Cl_ = 55.45;

    Gamma_[0] = 160.0;
    Gamma_[1] = 738.6;
    Gamma_[2] = 2688.0;
    Gamma_[3] = 4360.0;

    F_10Be_[0] = 0.9724; F_10Be_[1] = 0.0186; F_10Be_[2] = 0.004;  F_10Be_[3] = 0.005;
    F_26Al_[0] = 0.9655; F_26Al_[1] = 0.0233; F_26Al_[2] = 0.005;  F_26Al_[3] = 0.0062;
    F_14C_[0]  = 0.83;   F_14C_[1]  = 0.0691; F_14C_[2]  = 0.0809; F_14C_[3]  = 0.02;
}

// ─────────────────────────────────────────────────────────────────────────────
void lsdtt_crn_parameters::init_standardisation_maps()
{
    // 10Be: reference standard is 07KNSTD (factor = 1.0)
    // Values convert measured standard → 07KNSTD equivalent
    be10_std_map_["07KNSTD"]      = 1.0;
    be10_std_map_["KNSTD"]        = 0.9042;
    be10_std_map_["NIST_Certified"] = 1.0425;
    be10_std_map_["LLNL31000"]    = 0.8761;
    be10_std_map_["LLNL10000"]    = 0.9042;
    be10_std_map_["LLNL3000"]     = 0.8644;
    be10_std_map_["LLNL1000"]     = 0.9313;
    be10_std_map_["LLNL300"]      = 0.8562;
    be10_std_map_["NIST_30000"]   = 0.9313;
    be10_std_map_["NIST_30200"]   = 0.9251;
    be10_std_map_["NIST_30300"]   = 0.9221;
    be10_std_map_["NIST_30600"]   = 0.9130;
    be10_std_map_["NIST_27900"]   = 1.0;
    be10_std_map_["S555"]         = 0.9124;
    be10_std_map_["S2007"]        = 0.9124;
    be10_std_map_["BEST433"]      = 0.9124;
    be10_std_map_["BEST433N"]     = 1.0;
    be10_std_map_["S555N"]        = 1.0;
    be10_std_map_["S2007N"]       = 1.0;

    // 26Al: reference standard is KNSTD (factor = 1.0)
    al26_std_map_["KNSTD"]  = 1.0;
    al26_std_map_["ZAL94"]  = 0.9134;
    al26_std_map_["SMAL11"] = 1.021;
    al26_std_map_["0"]      = 1.0;
    al26_std_map_["ZAL94N"] = 1.0;
    al26_std_map_["ASTER"]  = 1.021;
    al26_std_map_["Z92-0222"] = 1.0;
}

// ─────────────────────────────────────────────────────────────────────────────
void lsdtt_crn_parameters::set_braucher_parameters()
{
    // Cosmocalc v2.0 (Braucher et al 2011) — Chmeleff/Korschinek decay constant
    lambda_10Be_ = 500e-9;
    lambda_26Al_ = 980e-9;
    lambda_14C_  = 121e-6;
    lambda_36Cl_ = 230e-8;

    P0_10Be_ = 4.30;
    P0_26Al_ = 31.10;
    P0_14C_  = 15.21;
    P0_36Cl_ = 58.95;

    Gamma_[0] = 160.0; Gamma_[1] = 1500.0; Gamma_[2] = 1500.0; Gamma_[3] = 4320.0;

    F_10Be_[0] = 0.9887; F_10Be_[1] = 0.0027; F_10Be_[2] = 0.0;   F_10Be_[3] = 0.0086;
    F_26Al_[0] = 0.9699; F_26Al_[1] = 0.0275; F_26Al_[2] = 0.000; F_26Al_[3] = 0.0026;
    F_14C_[0]  = 0.83;   F_14C_[1]  = 0.15;   F_14C_[2]  = 0.0;   F_14C_[3]  = 0.02;
}

// ─────────────────────────────────────────────────────────────────────────────
void lsdtt_crn_parameters::set_braucher_borchers_parameters()
{
    // Braucher et al EPSL 2011 F values + Borchers et al 2016 P0 values,
    // compiled by Mirjam Schaller.  14C uses Lupker et al QG 2015 spallation.
    lambda_10Be_ = 500e-9;
    lambda_26Al_ = 980e-9;
    lambda_14C_  = 121.58e-6;   // updated half-life (5700 yr), IAEA 2024
    lambda_36Cl_ = 230e-8;

    P0_10Be_ = 4.061;
    P0_26Al_ = 28.851;
    P0_14C_  = 15.533;
    P0_36Cl_ = 55.705;

    Gamma_[0] = 160.0; Gamma_[1] = 1500.0; Gamma_[2] = 1500.0; Gamma_[3] = 4320.0;

    F_10Be_[0] = 0.9874; F_10Be_[1] = 0.0030; F_10Be_[2] = 0.0;   F_10Be_[3] = 0.0096;
    F_26Al_[0] = 0.9681; F_26Al_[1] = 0.0291; F_26Al_[2] = 0.000; F_26Al_[3] = 0.0028;
    F_14C_[0]  = 0.788;  F_14C_[1]  = 0.212;  F_14C_[2]  = 0.0;   F_14C_[3]  = 0.0;
}

// ─────────────────────────────────────────────────────────────────────────────
void lsdtt_crn_parameters::set_granger_parameters()
{
    lambda_10Be_ = 500e-9;
    lambda_26Al_ = 980e-9;
    lambda_14C_  = 121e-6;
    lambda_36Cl_ = 230e-8;

    P0_10Be_ = 5.11;
    P0_26Al_ = 30.31;
    P0_14C_  = 5.86;
    P0_36Cl_ = 55.45;

    Gamma_[0] = 160.0; Gamma_[1] = 738.6; Gamma_[2] = 2688.0; Gamma_[3] = 4360.0;

    F_10Be_[0] = 0.9724; F_10Be_[1] = 0.0186; F_10Be_[2] = 0.004;  F_10Be_[3] = 0.005;
    F_26Al_[0] = 0.9655; F_26Al_[1] = 0.0233; F_26Al_[2] = 0.005;  F_26Al_[3] = 0.0062;
    F_14C_[0]  = 0.83;   F_14C_[1]  = 0.0691; F_14C_[2]  = 0.0809; F_14C_[3]  = 0.02;
}

// ─────────────────────────────────────────────────────────────────────────────
void lsdtt_crn_parameters::set_schaller_parameters()
{
    // Schaller et al 2016 — updated Gamma[1] for negative muons
    lambda_10Be_ = 500e-9;
    lambda_26Al_ = 980e-9;
    lambda_14C_  = 121e-6;
    lambda_36Cl_ = 230e-8;

    P0_10Be_ = 4.01;
    P0_26Al_ = 28.84;
    P0_14C_  = 12.1;
    P0_36Cl_ = 55.45;

    Gamma_[0] = 160.0; Gamma_[1] = 1500.0; Gamma_[2] = 1500.0; Gamma_[3] = 4320.0;

    F_10Be_[0] = 0.9884; F_10Be_[1] = 0.0042; F_10Be_[2] = 0.0;    F_10Be_[3] = 0.0074;
    F_26Al_[0] = 0.9697; F_26Al_[1] = 0.0270; F_26Al_[2] = 0.0;    F_26Al_[3] = 0.0033;
    F_14C_[0]  = 0.83;   F_14C_[1]  = 0.15;   F_14C_[2]  = 0.0;    F_14C_[3]  = 0.02;
}

// ─────────────────────────────────────────────────────────────────────────────
void lsdtt_crn_parameters::set_new_cronus_parameters()
{
    // New CRONUS 10Be production curve (Shasta Marrero); all other nuclides
    // from Braucher et al 2009 / Cosmocalc v2.0.
    lambda_10Be_ = 500e-9;
    lambda_26Al_ = 980e-9;
    lambda_14C_  = 121e-6;
    lambda_36Cl_ = 230e-8;

    P0_10Be_ = 4.075213;
    P0_26Al_ = 27.9;
    P0_14C_  = 12.2;
    P0_36Cl_ = 58.95;

    Gamma_[0] = 160.0;
    Gamma_[1] = 1459.76761923;
    Gamma_[2] = 11039.2402217;
    Gamma_[3] = 4320.0;

    F_10Be_[0] = 0.98374;          F_10Be_[1] = 0.0137188126531;
    F_10Be_[2] = 0.00252519164093; F_10Be_[3] = 0.0;

    F_26Al_[0] = 0.9699; F_26Al_[1] = 0.0275; F_26Al_[2] = 0.000; F_26Al_[3] = 0.0026;
    F_14C_[0]  = 0.83;   F_14C_[1]  = 0.15;   F_14C_[2]  = 0.0;   F_14C_[3]  = 0.02;
}

// ─────────────────────────────────────────────────────────────────────────────
void lsdtt_crn_parameters::reset_scaling(const std::string& muon_scheme)
{
    if      (muon_scheme == "Braucher")          set_braucher_parameters();
    else if (muon_scheme == "BraucherBorchers")  set_braucher_borchers_parameters();
    else if (muon_scheme == "Granger")           set_granger_parameters();
    else if (muon_scheme == "Schaller")          set_schaller_parameters();
    else if (muon_scheme == "newCRONUS")         set_new_cronus_parameters();
    else
    {
        std::cout << "lsdtt_crn_parameters: unknown muon scheme \"" << muon_scheme
                  << "\". Defaulting to BraucherBorchers." << std::endl;
        set_braucher_borchers_parameters();
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// scale_one_nuclide: Newton-Raphson to find effective depth x such that
//   Σᵢ F[i]*exp(-x/Γ[i]) == total_shielding
// then rescale: F[i] *= exp(-x/Γ[i]).
// This preserves the relative attenuation shape of each production pathway
// rather than applying a flat multiplier.
void lsdtt_crn_parameters::scale_one_nuclide(double total_shielding,
                                              const double Gamma[4],
                                              double F[4])
{
    if (total_shielding >= 1.0)
        return;   // no shielding — nothing to do

    const double tolerance    = 1e-7;
    const double displace_x   = 1e-6;
    double x                  = 0.0;
    double x_change           = 1.0;

    while (std::fabs(x_change) > tolerance)
    {
        double s  = 0.0, sd = 0.0;
        for (int i = 0; i < 4; ++i)
        {
            s  += F[i] * std::exp(-x           / Gamma[i]);
            sd += F[i] * std::exp(-(x + displace_x) / Gamma[i]);
        }
        double fx        = s  - total_shielding;
        double fx_d      = sd - total_shielding;
        double deriv     = (fx_d - fx) / displace_x;
        if (deriv == 0.0) break;
        x_change = fx / deriv;
        x -= x_change;
    }

    for (int i = 0; i < 4; ++i)
        F[i] *= std::exp(-x / Gamma[i]);
}

// ─────────────────────────────────────────────────────────────────────────────
void lsdtt_crn_parameters::scale_F_values(double total_shielding,
                                           const std::string& nuclide)
{
    if (nuclide == "Be10")
        scale_one_nuclide(total_shielding, Gamma_, F_10Be_);
    else if (nuclide == "Al26")
        scale_one_nuclide(total_shielding, Gamma_, F_26Al_);
    else if (nuclide == "C14")
        scale_one_nuclide(total_shielding, Gamma_, F_14C_);
    else
    {
        std::cout << "lsdtt_crn_parameters::scale_F_values: unknown nuclide \""
                  << nuclide << "\". No scaling applied." << std::endl;
    }
}

// ─────────────────────────────────────────────────────────────────────────────
double lsdtt_crn_parameters::atmospheric_pressure_standard(double elevation_m)
{
    // Balco barometric formula with SLHL standard-atmosphere constants.
    // Returns pressure in mbar (hPa).
    const double gmr   = -0.03417;   // g·M/R (mol⁻¹ K⁻¹)
    const double dtdz  =  0.0065;    // environmental lapse rate (K m⁻¹)
    const double P_slp =  1013.25;   // standard sea-level pressure (mbar)
    const double T_slp =  288.15;    // standard sea-level temperature (K)
    double T_site = T_slp - dtdz * elevation_m;
    if (T_site <= 0.0) T_site = 1.0;   // guard against unphysical elevations
    return P_slp * std::exp((gmr / dtdz) * (std::log(T_slp) - std::log(T_site)));
}

// ─────────────────────────────────────────────────────────────────────────────
// Stone (2000) spallation scaling.  Faithfully ported from LSDTopoTools2
// (Balco / Stone 2000 Table 1 coefficients).
double lsdtt_crn_parameters::stone2000sp(double lat_deg,
                                          double P_mbar,
                                          double Fsp) const
{
    if (Fsp > 1.0) Fsp = 0.978;
    if (std::fabs(lat_deg) > 90.0)
    {
        std::cout << "stone2000sp: |latitude| > 90. Defaulting to 45 degrees." << std::endl;
        lat_deg = 45.0;
    }

    // Spallogenic production at 7 index latitudes (0–60°), from Stone 2000 Table 1.
    // S(P) = a + b·exp(-P/150) + c·P + d·P² + e·P³
    const double a[7] = { 31.8518, 34.3699, 40.3153, 42.0983, 56.7733, 69.0720, 71.8733 };
    const double b[7] = { 250.3193, 258.4759, 308.9894, 512.6857, 649.1343, 832.4566, 863.1927 };
    const double c[7] = { -0.083393, -0.089807, -0.106248, -0.120551, -0.160859, -0.199252, -0.207069 };
    const double d[7] = { 7.4260e-5, 7.9457e-5, 9.4508e-5, 1.1752e-4, 1.5463e-4, 1.9391e-4, 2.0127e-4 };
    const double e[7] = { -2.2397e-8, -2.3697e-8, -2.8234e-8, -3.8809e-8, -5.0330e-8, -6.3653e-8, -6.6043e-8 };
    const double ilat[7] = { 0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0 };

    // Muonic scaling constants at 7 index latitudes
    const double mk[7] = { 0.587, 0.600, 0.678, 0.833, 0.933, 1.000, 1.000 };

    const double P2 = P_mbar * P_mbar;
    const double P3 = P2 * P_mbar;
    const double eP = std::exp(P_mbar / (-150.0));

    std::vector<double> S_at_ilat(7), M_at_ilat(7);
    for (int k = 0; k < 7; ++k)
    {
        S_at_ilat[k] = a[k] + b[k] * eP + c[k] * P_mbar + d[k] * P2 + e[k] * P3;
        M_at_ilat[k] = mk[k] * std::exp((1013.25 - P_mbar) / 242.0);
    }

    std::vector<double> ilat_v(ilat, ilat + 7);

    lat_deg = std::fabs(lat_deg);
    if (lat_deg > 60.0) lat_deg = 60.0;

    double S = interp1d(ilat_v, S_at_ilat, lat_deg);
    double M = interp1d(ilat_v, M_at_ilat, lat_deg);

    return (S * Fsp) + (M * (1.0 - Fsp));
}

// ─────────────────────────────────────────────────────────────────────────────
double lsdtt_crn_parameters::be10_standardisation_factor(
    const std::string& standard_name) const
{
    auto it = be10_std_map_.find(standard_name);
    if (it == be10_std_map_.end())
    {
        std::cout << "be10_standardisation_factor: unknown standard \""
                  << standard_name << "\". Assuming 07KNSTD (factor=1)." << std::endl;
        return 1.0;
    }
    return it->second;
}

double lsdtt_crn_parameters::al26_standardisation_factor(
    const std::string& standard_name) const
{
    auto it = al26_std_map_.find(standard_name);
    if (it == al26_std_map_.end())
    {
        std::cout << "al26_standardisation_factor: unknown standard \""
                  << standard_name << "\". Assuming KNSTD (factor=1)." << std::endl;
        return 1.0;
    }
    return it->second;
}

// ─────────────────────────────────────────────────────────────────────────────
std::vector<double> lsdtt_crn_parameters::apply_standardisation(
    const std::vector<double>& concentrations,
    const std::vector<std::string>& standard_names,
    const std::string& nuclide) const
{
    std::vector<double> out(concentrations.size());
    bool is_10Be = (nuclide == "Be10");

    for (std::size_t i = 0; i < concentrations.size(); ++i)
    {
        std::string std_name = (i < standard_names.size()) ? standard_names[i] : "";
        double factor = 1.0;
        if (!std_name.empty())
            factor = is_10Be ? be10_standardisation_factor(std_name)
                             : al26_standardisation_factor(std_name);
        out[i] = concentrations[i] * factor;
    }
    return out;
}

// ─────────────────────────────────────────────────────────────────────────────
double lsdtt_crn_parameters::interp1d(const std::vector<double>& x,
                                       const std::vector<double>& y,
                                       double xi)
{
    // Linear interpolation on a monotone-increasing x vector.
    // Clamps to boundary values outside [x.front(), x.back()].
    int n = static_cast<int>(x.size());
    if (xi <= x[0])       return y[0];
    if (xi >= x[n - 1])   return y[n - 1];
    for (int k = 0; k < n - 1; ++k)
    {
        if (xi >= x[k] && xi <= x[k + 1])
        {
            double t = (xi - x[k]) / (x[k + 1] - x[k]);
            return y[k] + t * (y[k + 1] - y[k]);
        }
    }
    return y[n - 1];
}

} // namespace lsdtt
