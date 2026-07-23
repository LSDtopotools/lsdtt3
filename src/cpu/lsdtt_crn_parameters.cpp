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
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <algorithm>
#ifdef LSDTT_HAS_NETCDF
#include <netcdf.h>
#endif

namespace lsdtt
{

namespace {
constexpr double kPi = 3.14159265358979323846;
}

// 
lsdtt_crn_parameters::lsdtt_crn_parameters()
{
    init_defaults();
    set_braucher_borchers_parameters();
    init_standardisation_maps();
}

// 
void lsdtt_crn_parameters::init_defaults()
{
    // Vermeesh (2007) default values - overwritten by any set_*_parameters call
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

// 
void lsdtt_crn_parameters::init_standardisation_maps()
{
    // 10Be: reference standard is 07KNSTD (factor = 1.0)
    // Values convert measured standard -> 07KNSTD equivalent
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

// 
void lsdtt_crn_parameters::set_braucher_parameters()
{
    // Cosmocalc v2.0 (Braucher et al 2011) - Chmeleff/Korschinek decay constant
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

// 
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

// 
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

// 
void lsdtt_crn_parameters::set_schaller_parameters()
{
    // Schaller et al 2016 - updated Gamma[1] for negative muons
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

// 
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

// 
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

//
// scale_one_nuclide: Newton-Raphson to find effective depth x such that
//   Sigma_i F[i]*exp(-x/Gamma[i]) == total_shielding
// then rescale: F[i] *= exp(-x/Gamma[i]).
// This preserves the relative attenuation shape of each production pathway
// rather than applying a flat multiplier.
//
// total_shielding here is production scaling x topographic shielding, so it
// is routinely > 1 (Stone scaling exceeds 1 above a few hundred metres);
// the Vermeesch (2007) scheme handles that with a NEGATIVE effective depth
// (amplification, spallation amplified most).  Fix approved by Simon Mudd
// 2026-07-16: the previous ">= 1.0 -> return" guard silently dropped the
// production scaling at altitude, found by the CRONUS v3 comparison tests.
void lsdtt_crn_parameters::scale_one_nuclide(double total_shielding,
                                              const double Gamma[4],
                                              double F[4])
{
    if (total_shielding == 1.0)
        return;   // exactly neutral - nothing to do

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

// 
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

// 
double lsdtt_crn_parameters::atmospheric_pressure_standard(double elevation_m)
{
    // Balco barometric formula with SLHL standard-atmosphere constants.
    // Returns pressure in mbar (hPa).
    const double gmr   = -0.03417;   // g*M/R (mol^-^1 K^-^1)
    const double dtdz  =  0.0065;    // environmental lapse rate (K m^-^1)
    const double P_slp =  1013.25;   // standard sea-level pressure (mbar)
    const double T_slp =  288.15;    // standard sea-level temperature (K)
    double T_site = T_slp - dtdz * elevation_m;
    if (T_site <= 0.0) T_site = 1.0;   // guard against unphysical elevations
    return P_slp * std::exp((gmr / dtdz) * (std::log(T_slp) - std::log(T_site)));
}

// 
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

    // Spallogenic production at 7 index latitudes (0-60 deg), from Stone 2000 Table 1.
    // S(P) = a + b*exp(-P/150) + c*P + d*P^2 + e*P^3
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

// 
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

// 
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

// 
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

// 
// load_nc_data
// 
CrnNcData lsdtt_crn_parameters::load_nc_data(const std::string& nc_path)
{
    CrnNcData out;
#ifdef LSDTT_HAS_NETCDF
    int ncid;
    if (nc_open(nc_path.c_str(), NC_NOWRITE, &ncid) != NC_NOERR)
    {
        std::cerr << "load_nc_data: cannot open \"" << nc_path << "\"\n";
        return out;
    }

    auto read_var1d_d = [&](const char* name, std::vector<double>& v) {
        int varid;
        if (nc_inq_varid(ncid, name, &varid) != NC_NOERR) return;
        int ndims; int dimids[8];
        nc_inq_var(ncid, varid, nullptr, nullptr, &ndims, dimids, nullptr);
        std::size_t n = 1;
        for (int d = 0; d < ndims; ++d) {
            std::size_t s; nc_inq_dimlen(ncid, dimids[d], &s); n *= s;
        }
        v.resize(n);
        nc_get_var_double(ncid, varid, v.data());
    };

    auto read_var2d_f = [&](const char* name, std::vector<float>& v) {
        int varid;
        if (nc_inq_varid(ncid, name, &varid) != NC_NOERR) return;
        int ndims; int dimids[8];
        nc_inq_var(ncid, varid, nullptr, nullptr, &ndims, dimids, nullptr);
        std::size_t n = 1;
        for (int d = 0; d < ndims; ++d) {
            std::size_t s; nc_inq_dimlen(ncid, dimids[d], &s); n *= s;
        }
        v.resize(n);
        nc_get_var_float(ncid, varid, v.data());
    };

    auto read_var1d_f = [&](const char* name, std::vector<float>& v) {
        int varid;
        if (nc_inq_varid(ncid, name, &varid) != NC_NOERR) return;
        int ndims; int dimids[8];
        nc_inq_var(ncid, varid, nullptr, nullptr, &ndims, dimids, nullptr);
        std::size_t n = 1;
        for (int d = 0; d < ndims; ++d) {
            std::size_t s; nc_inq_dimlen(ncid, dimids[d], &s); n *= s;
        }
        v.resize(n);
        nc_get_var_float(ncid, varid, v.data());
    };

    read_var1d_d("era40_lat",       out.era40_lat);
    read_var1d_d("era40_lon",       out.era40_lon);
    read_var2d_f("era40_meanP",     out.era40_meanP);
    read_var2d_f("era40_meanT",     out.era40_meanT);
    read_var1d_d("lifton_lr",       out.lifton_lr);
    read_var1d_d("muscheler_time",  out.muscheler_time);
    read_var1d_f("muscheler_vdm",   out.muscheler_vdm);

    nc_close(ncid);

    out.loaded = (!out.era40_lat.empty() && !out.era40_meanP.empty() &&
                  !out.lifton_lr.empty() && !out.muscheler_vdm.empty());
    if (out.loaded)
        std::cout << "load_nc_data: loaded " << nc_path
                  << "  ERA40=" << out.era40_lat.size() << "x" << out.era40_lon.size()
                  << "  Muscheler=" << out.muscheler_time.size() << " steps\n";
    else
        std::cerr << "load_nc_data: file opened but variables missing in \"" << nc_path << "\"\n";
#else
    std::cerr << "load_nc_data: lsdtt3 was compiled without netCDF-C support.\n";
#endif
    return out;
}

// 
// bilinear_era40
// Grid axes: lat_ax descending (90->-90), lon_ax ascending (0->360).
// 
double lsdtt_crn_parameters::bilinear_era40(const std::vector<double>& lat_ax,
                                             const std::vector<double>& lon_ax,
                                             const std::vector<float>&  grid,
                                             double lat, double lon)
{
    int nlat = static_cast<int>(lat_ax.size());
    int nlon = static_cast<int>(lon_ax.size());
    if (nlat < 2 || nlon < 2) return std::numeric_limits<double>::quiet_NaN();

    // Clamp lat to grid bounds
    lat = std::max(lat_ax.back(), std::min(lat_ax.front(), lat));
    // Wrap lon to [0, 360)
    while (lon <   0.0) lon += 360.0;
    while (lon > 360.0) lon -= 360.0;

    // Find lon bracket (ascending)
    int jlo = 0;
    for (int j = 0; j < nlon - 1; ++j) {
        if (lon >= lon_ax[j] && lon <= lon_ax[j+1]) { jlo = j; break; }
        if (j == nlon - 2) jlo = j;
    }
    double tlon = (lon_ax[jlo+1] > lon_ax[jlo])
                  ? (lon - lon_ax[jlo]) / (lon_ax[jlo+1] - lon_ax[jlo])
                  : 0.0;

    // Find lat bracket (descending - lat_ax[0] is largest)
    int ilo = 0;
    for (int i = 0; i < nlat - 1; ++i) {
        if (lat <= lat_ax[i] && lat >= lat_ax[i+1]) { ilo = i; break; }
        if (i == nlat - 2) ilo = i;
    }
    double tlat = (lat_ax[ilo] > lat_ax[ilo+1])
                  ? (lat_ax[ilo] - lat) / (lat_ax[ilo] - lat_ax[ilo+1])
                  : 0.0;

    // Bilinear: corners (ilo,jlo) (ilo,jlo+1) (ilo+1,jlo) (ilo+1,jlo+1)
    auto idx = [&](int i, int j) { return static_cast<std::size_t>(i) * nlon + j; };
    double v00 = grid[idx(ilo,   jlo)];
    double v01 = grid[idx(ilo,   jlo+1)];
    double v10 = grid[idx(ilo+1, jlo)];
    double v11 = grid[idx(ilo+1, jlo+1)];

    return (1-tlat)*(1-tlon)*v00 + (1-tlat)*tlon*v01
           + tlat*(1-tlon)*v10  + tlat*tlon*v11;
}

// 
// atmospheric_pressure_era40
// Mirrors ERA40atm() from lsd.py (Balco/Lifton, CRONUS origin).
// 
double lsdtt_crn_parameters::atmospheric_pressure_era40(double lat_deg, double lon_deg,
                                                         double elevation_m,
                                                         const CrnNcData& nc)
{
    if (!nc.loaded || nc.era40_lat.empty() || nc.era40_meanP.empty())
        return atmospheric_pressure_standard(elevation_m);

    // Normalise longitude to [0, 360)
    if (lon_deg < 0.0) lon_deg += 360.0;

    double P_sl = bilinear_era40(nc.era40_lat, nc.era40_lon, nc.era40_meanP, lat_deg, lon_deg);
    double T_sl = bilinear_era40(nc.era40_lat, nc.era40_lon, nc.era40_meanT, lat_deg, lon_deg);

    // Lifton lapse rate polynomial: dtdz(lat) = -Sigma_k lr[k]*lat^k  (K m^-^1, positive upward)
    // Signed latitude to match COSPAR CIRA-86 fit (lsd.py uses site_lat signed).
    double dtdz = 0.0;
    double lat_k = 1.0;
    for (std::size_t k = 0; k < nc.lifton_lr.size(); ++k) {
        dtdz += nc.lifton_lr[k] * lat_k;
        lat_k *= lat_deg;
    }
    dtdz = -dtdz;   // sign convention: dtdz > 0 means temperature falls with altitude

    const double gmr = -0.03417;   // -g*M/R (Balco 2008)
    double T_site = T_sl - dtdz * elevation_m;
    if (T_site <= 0.0) T_site = 1.0;   // guard unphysical elevations
    return P_sl * std::exp((gmr / dtdz) * (std::log(T_sl) - std::log(T_site)));
}

// 
// stone_factor_at_rc
// Converts geomagnetic cutoff rigidity Rc [GV] to Stone 2000 spallation scaling
// factor at pressure P_mbar.  Mirrors Rc_StoneFactors() from PR_processes.py.
// 
double lsdtt_crn_parameters::stone_factor_at_rc(double Rc_GV, double P_mbar)
{
    // Stone table: latitudes, coefficients (same values as in stone2000sp)
    static const double ilat[7]  = { 0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0 };
    static const double sa[7] = { 31.8518, 34.3699, 40.3153, 42.0983, 56.7733, 69.0720, 71.8733 };
    static const double sb[7] = { 250.3193, 258.4759, 308.9894, 512.6857, 649.1343, 832.4566, 863.1927 };
    static const double sc[7] = { -0.083393,-0.089807,-0.106248,-0.120551,-0.160859,-0.199252,-0.207069 };
    static const double sd[7] = { 7.4260e-5,7.9457e-5,9.4508e-5,1.1752e-4,1.5463e-4,1.9391e-4,2.0127e-4 };
    static const double se[7] = { -2.2397e-8,-2.3697e-8,-2.8234e-8,-3.8809e-8,-5.0330e-8,-6.3653e-8,-6.6043e-8 };

    const double P  = P_mbar;
    const double P2 = P * P, P3 = P2 * P;
    const double eP = std::exp(-P / 150.0);

    // Reference Rc values from Stone table latitudes: i_Rc = 14.3 * cos^4(lat_rad)
    // and corresponding spallation factors at this pressure
    std::vector<double> i_Rc(7), i_sf(7);
    for (int k = 0; k < 7; ++k) {
        double lat_r = ilat[k] * kPi / 180.0;
        double c4 = std::pow(std::cos(lat_r), 4.0);
        i_Rc[k] = 14.3 * c4;
        i_sf[k] = sa[k] + sb[k] * eP + sc[k] * P + sd[k] * P2 + se[k] * P3;
    }

    // Extend to Rc = 0: Stone factor at 60 deg applies below this rigidity
    i_Rc.push_back(0.0);
    i_sf.push_back(i_sf[6]);

    // Extend to Rc = 25 GV using log-log fit to the first 3 points (lat 0,10,20)
    // fitting: log(sf) = slope * log(Rc) + intercept
    double slope = (std::log(i_sf[0]) - std::log(i_sf[2]))
                 / (std::log(i_Rc[0]) - std::log(i_Rc[2]));
    for (int rc_ext = 25; rc_ext > 15; --rc_ext) {
        double sf_ext = std::exp(std::log(i_sf[0])
                                 + slope * (std::log(rc_ext) - std::log(i_Rc[0])));
        i_Rc.push_back(static_cast<double>(rc_ext));
        i_sf.push_back(sf_ext);
    }

    // Sort both vectors by ascending Rc for interpolation
    std::size_t n = i_Rc.size();
    std::vector<std::size_t> idx(n);
    for (std::size_t i = 0; i < n; ++i) idx[i] = i;
    std::sort(idx.begin(), idx.end(), [&](std::size_t a, std::size_t b) {
        return i_Rc[a] < i_Rc[b];
    });
    std::vector<double> Rc_s(n), sf_s(n);
    for (std::size_t i = 0; i < n; ++i) { Rc_s[i] = i_Rc[idx[i]]; sf_s[i] = i_sf[idx[i]]; }

    return interp1d(Rc_s, sf_s, Rc_GV);
}

// 
// stone2000sp_paleomag
// Stone 2000 scaling with time-integrated paleomagnetic correction via the
// Muscheler VDM record.  Follows Rc_correctedStoneFactors() / PR_processes.py.
// 
double lsdtt_crn_parameters::stone2000sp_paleomag(double lat_deg, double P_mbar,
                                                   double Fsp, double age_yr,
                                                   const CrnNcData& nc) const
{
    if (!nc.loaded || nc.muscheler_time.empty() || age_yr <= 0.0)
        return stone2000sp(lat_deg, P_mbar, Fsp);

    const double lat_rad = std::fabs(lat_deg) * kPi / 180.0;
    const double cos4lat = std::pow(std::cos(lat_rad), 4.0);

    // Rc [GV] from VDM [10^2^2 Am^2]: dipole formula (Stormer unit)
    // Rc = (mu_0 x M x 1e22 / 4pi) / R^2 x c x cos^4(lat) / 1e9
    // = (1e-7 x M x 1e22 x 3e8) / (6.3712e6)^2 x cos^4(lat)  - but must match Python exactly:
    // Python: (p * 1e22 * 4e-7 * 3e8) / (16e9 * (6.3712e6)^2) * cos^4(lat)
    static const double Rc_c = (4e-7 * 3e8 * 1e22) / (16e9 * 6.3712e6 * 6.3712e6);

    // Locate the portion of Muscheler record within [0, age_yr]
    int n = static_cast<int>(nc.muscheler_time.size());
    int n_use = 0;
    for (int i = 0; i < n; ++i) {
        if (nc.muscheler_time[i] <= age_yr) n_use = i + 1; else break;
    }
    if (n_use < 2)
        return stone2000sp(lat_deg, P_mbar, Fsp);

    // Build Stone spallation factor at each time step
    std::vector<double> sf_t(n_use);
    for (int i = 0; i < n_use; ++i) {
        double Rc = Rc_c * static_cast<double>(nc.muscheler_vdm[i]) * cos4lat;
        sf_t[i] = stone_factor_at_rc(Rc, P_mbar);
    }

    // Trapezoidal time average of sf over [time[0], time[n_use-1]]
    double integral = 0.0;
    for (int i = 1; i < n_use; ++i) {
        double dt = nc.muscheler_time[i] - nc.muscheler_time[i-1];
        integral += 0.5 * (sf_t[i] + sf_t[i-1]) * dt;
    }
    double t_span = nc.muscheler_time[n_use-1] - nc.muscheler_time[0];
    if (t_span <= 0.0)
        return stone2000sp(lat_deg, P_mbar, Fsp);
    double sf_avg = integral / t_span;

    // Modern spallation factor (VDM at t ~= 0, i.e. index 0)
    double Rc_modern = Rc_c * static_cast<double>(nc.muscheler_vdm[0]) * cos4lat;
    double sf_modern = stone_factor_at_rc(Rc_modern, P_mbar);
    if (sf_modern <= 0.0)
        return stone2000sp(lat_deg, P_mbar, Fsp);

    // Scale the total Stone factor by the paleomagnetic correction ratio
    // (correction applies to the spallation component only; muon term unchanged)
    double correction = sf_avg / sf_modern;
    return stone2000sp(lat_deg, P_mbar, Fsp) * correction;
}

} // namespace lsdtt
