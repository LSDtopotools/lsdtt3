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
#include <string>
#include <map>
#include <vector>

namespace lsdtt
{

// 
// CrnNcData
//
// ERA40 atmospheric grids and Muscheler VDM record loaded from crn_parameters.nc.
// Load once with lsdtt_crn_parameters::load_nc_data(), then pass by const-ref
// to ERA40 and paleomagnetic correction functions.
// 
struct CrnNcData
{
    // ERA40 grid axes
    std::vector<double> era40_lat;     // (73,)   deg, 90 -> -90 (descending)
    std::vector<double> era40_lon;     // (145,)  deg, 0 -> 360 (ascending)
    // ERA40 climatological fields - row-major [ilat][ilon]
    std::vector<float>  era40_meanP;   // (73x145) hPa
    std::vector<float>  era40_meanT;   // (73x145) K
    // Lifton COSPAR CIRA-86 lapse rate polynomial coefficients (7 terms, K m^-^1)
    std::vector<double> lifton_lr;     // (7,)
    // Muscheler VDM record
    std::vector<double> muscheler_time; // (5970,) yr BP, ascending
    std::vector<float>  muscheler_vdm;  // (5970,) 10^2^2 Am^2
    bool loaded = false;
    // NCEP2 reanalysis climatology (the atmosphere CAIRN/OCTOPUS use; Balco
    // NCEPatm_2). Optional: absent in older crn_parameters.nc files. Only the
    // sea-level pressure + 1000 hPa temperature grids are needed - the legacy
    // geopotential-height arrays are never used in the pressure calculation.
    std::vector<double> ncep_lat;      // (73,)   deg, 90 -> -90 (descending)
    std::vector<double> ncep_lon;      // (145,)  deg, 0 -> 360 (ascending)
    std::vector<float>  ncep_meanP;    // (73x145) hPa (mean sea-level pressure)
    std::vector<float>  ncep_meanT;    // (73x145) degC (mean 1000 hPa temperature)
    bool ncep_loaded = false;

    // LSDn (Lifton, Sato & Dunai 2014) nuclide-specific scaling grids and the
    // paleomagnetic / solar-modulation records they are driven by. Optional:
    // absent in older crn_parameters.nc files.
    std::vector<double> lsdn_rc;       // (201,) GV, ascending 0 -> 20 step 0.1
    std::vector<double> lsdn_p;        // (380,) hPa, ascending 300 -> 1058 step 2
    std::vector<float>  lsdn_b_Be10;   // (380x201) row-major [ip][irc]
    std::vector<float>  lsdn_m_Be10;   // (380x201)
    std::vector<float>  lsdn_b_Al26, lsdn_m_Al26;
    std::vector<float>  lsdn_b_C14,  lsdn_m_C14;
    std::vector<double> lsdn_sphi;     // (91,)  MV, per Holocene time step
    double              lsdn_sphi_inf = 0.0;            // MV, long-term mean
    std::vector<double> pavon_t1min, pavon_t1max;       // (91,)   yr
    std::vector<double> pavon_mm01;                     // (91,)   M_t / M_0
    std::vector<double> pavon_pole_lat, pavon_pole_lon; // (91,)   rad
    std::vector<double> pavon_t2min, pavon_t2max;       // (1710,) yr
    std::vector<double> pavon_mm02;                     // (1710,) M_t / M_0
    // Long-record step CENTRES. Used ONLY as the branch threshold for
    // "exposure longer than the long record" (get_DipRc.m uses m.t2(end) =
    // 2.00074e7 yr, NOT t2max.back() = 2.02044e7 yr). Never a time weight -
    // the duration averaging uses the bounds above.
    std::vector<double> pavon_t2;                       // (1710,) yr
    double              pavon_tzero = 0.0;              // yr CE (2010)
    bool                lsdn_loaded = false;
};

//
// LsdnTimeAvgSurface
//
// Time-averaged LSDn scaling precomputed over a horizontal (lat, lon) grid and
// the native pressure nodes of the packed LSDn grids, for ONE nuclide and ONE
// sample age.
//
// Why it exists: the Pavon-Carrasco record has 91 Holocene steps plus up to
// 1710 long-record steps, so the direct per-pixel time loop of lsdn_scaling()
// costs up to ~1800 grid evaluations per pixel.  But the time average
//
//   Abar(lat, lon, P) = < b(Rc(t, lat, lon), P) + S(t)*m(Rc(t, lat, lon), P) >_t
//
// depends on the site only through (lat, lon, P) - the time axis, the field
// record and the solar record are common to every pixel.  So the average is
// built once over a modest horizontal grid and evaluated per pixel with a
// single trilinear interpolation.
//
// LONGITUDE IS PART OF THE SURFACE, not an optimisation to be dropped: the
// Holocene geomagnetic latitude is the great-circle distance to a MOVING
// paleopole, so Rc(t) depends on site longitude too.  A latitude-only surface
// would carry an error growing with distance from the paleopole meridian, and
// a latitude-only sensitivity sweep would never see it.
//
// THE PRESSURE AXIS IS EXACT, NOT APPROXIMATED.  b and m are linear in P
// between adjacent grid rows and the time average is a linear operator, so
// building on the 380 native pressure nodes and then interpolating linearly in
// P is identical (to floating point) to time-averaging at an arbitrary P.  The
// ONLY approximation in this object is the horizontal one.
//
// LIFETIME AND THREADING: built ONCE per raster call, before the OpenMP
// region; thereafter READ-ONLY and safe to share by const& across threads.
// Evaluation allocates nothing and touches no global or static mutable state,
// so it can never go stale against a different nuclide, age or extent.
//
struct LsdnTimeAvgSurface
{
    std::string         nuclide;            // "Be10" | "Al26" | "C14" | "Stone"
    double              age_yr   = 0.0;
    std::vector<double> lat_ax;             // ascending, deg, regular
    // Ascending, deg, regular.  Kept in the CALLER'S longitude convention (the
    // convention of the bounds it was built from), NOT forced into [0,360):
    // normalising would tear a raster that straddles the prime meridian.
    // lsdn_scaling_from_surface() shifts a query longitude by +/-360 as needed,
    // so either convention works at call time.  May legitimately hold a single
    // node for a kernel with no longitude dependence (Stone paleomag).
    std::vector<double> lon_ax;
    // LSDn: copy of nc.lsdn_p (380 native nodes).  EMPTY for the Stone-paleomag
    // flavour, which has no pressure axis at all - see avg_modern below.
    std::vector<double> p_ax;
    // LSDn: n_lat*n_lon*n_p time-averaged RAW values (NOT SLHL-normalised),
    // indexed ((ilat*n_lon) + ilon)*n_p + ip.
    // Stone paleomag: n_lat * stone_paleomag_basis_n time-averaged Stone-table
    // interpolation weights, indexed ilat*stone_paleomag_basis_n + j.
    std::vector<double> avg;
    // STONE-PALEOMAG FLAVOUR ONLY (empty for LSDn).  The modern (record step 0)
    // counterpart of avg, same layout.  It exists so that the correction's
    // numerator AND denominator are both drawn from this one object at the same
    // latitude - which is what makes the constant-VDM ratio cancel identically
    // (test L14b / G1) rather than only approximately.
    std::vector<double> avg_modern;
    double              ref      = 0.0;     // SLHL denominator for this nuclide
    bool                valid    = false;   // false => caller uses the direct path
    // Why the surface was or was not built; for the run header.
    std::string         reason;
    // Diagnostics for the run header and the performance CSV.
    std::size_t         n_nodes  = 0;
    std::size_t         n_bytes  = 0;
    double              build_seconds = 0.0;
};

// 
// lsdtt_crn_parameters
//
// Holds all physics constants for cosmogenic nuclide (CRN) production and decay.
//
// Workflow per pixel:
//   1. reset_scaling(muon_scheme)   - restore standard F values for the scheme
//   2. scale_F_values(shielding, nuclide) - fold total shielding into F arrays
//   3. Pass const& to lsdtt_crn_particle for concentration calculation
//
// Supported nuclides: Be10, Al26, C14
// Supported muon schemes: "Braucher", "BraucherBorchers" (default),
//                         "Granger", "Schaller", "newCRONUS",
//                         "balcoCRONUSv3"
// 
class lsdtt_crn_parameters
{
public:
    lsdtt_crn_parameters();

    //  Parameter schemes 
    void set_braucher_parameters();
    void set_braucher_borchers_parameters();   // default
    void set_granger_parameters();
    void set_schaller_parameters();
    void set_new_cronus_parameters();
    void set_balco_cronus_v3_parameters(double P_mbar = 1013.25);

    // Selects scheme by name (case-sensitive: "Braucher", "BraucherBorchers",
    // "Granger", "Schaller", "newCRONUS", "balcoCRONUSv3").
    // Falls back to BraucherBorchers.
    void reset_scaling(const std::string& muon_scheme);

    //  Shielding / production scaling (Vermeesch 2007 virtual attenuation)
    // Rescales F values for the named nuclide using a Newton-Raphson search
    // that converts total_shielding (= production_scale x topo_shield,
    // dimensionless) into an equivalent effective depth so that each
    // production pathway is attenuated (or amplified) correctly.
    //
    // NOTE: total_shielding is routinely > 1 - the Stone (2000) production
    // scaling exceeds 1 above a few hundred metres elevation, and the
    // Vermeesch approximation handles this with a negative effective depth
    // (spallation amplified most, muon pathways least, mirroring how the
    // pathways scale with altitude).  Handling S_tot > 1 here is required to
    // reproduce CRONUS calculator results at altitude (verified against the
    // CRONUS v3 online calculator, 2026-07-16).
    // Call reset_scaling first, then scale_F_values.
    // nuclide: "Be10", "Al26", or "C14".
    void scale_F_values(double total_shielding, const std::string& nuclide);

    //  Atmospheric pressure 
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

    // NCEP2-based atmospheric pressure - a faithful port of Balco's NCEPatm_2
    // (the atmosphere CAIRN/OCTOPUS use). Bilinear interpolation of the NCEP2
    // sea-level pressure and 1000 hPa temperature grids, then the standard-
    // atmosphere barometric formula with a single fixed lapse rate (0.0065 K/m)
    // - NOT the Lifton latitude polynomial the ERA40 path uses. This is the ONLY
    // difference in atmosphere model between the two; the climatology grids differ
    // too. Falls back to atmospheric_pressure_standard if nc.ncep_loaded == false.
    // Returns pressure in mbar.
    static double atmospheric_pressure_ncep(double lat_deg, double lon_deg,
                                             double elevation_m,
                                             const CrnNcData& nc);

    // Load ERA40 and paleomagnetic data from a crn_parameters.nc file.
    // Returns CrnNcData with loaded=false if the file cannot be read or if
    // the code was compiled without netCDF-C support.
    static CrnNcData load_nc_data(const std::string& nc_path);

    //  Stone 2000 scaling 
    // Returns the dimensionless Stone (2000) production-rate scaling factor.
    // lat_deg: site latitude (deg, southern hemisphere OK - sign is ignored).
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
    //
    // EVALUATED THROUGH THE PRECOMPUTED SURFACE (see
    // build_stone_paleomag_surface below), on the same global latitude lattice
    // the raster path uses, so this function and calculate_production_raster
    // return the same number at the same site.  The direct per-timestep
    // definition remains permanently available as
    // stone_paleomag_correction_direct() and is the oracle the surface is
    // tested against (L14).
    double stone2000sp_paleomag(double lat_deg, double P_mbar, double Fsp,
                                double age_yr, const CrnNcData& nc) const;

    // ---- Stone paleomagnetic correction: definition, surface, evaluator -----
    //
    // The correction is sf_avg(lat,P) / sf_modern(lat,P), where sf_avg is the
    // trapezoidal time average of stone_factor_at_rc(Rc_c*vdm(t)*cos^4(lat), P)
    // over the Muscheler record truncated at age_yr, and sf_modern is the same
    // quantity at record step 0.  Multiplying the time-independent Stone factor
    // by it gives stone2000sp_paleomag.
    //
    // WHY A SURFACE.  Called per pixel the direct form costs 108 evaluations of
    // stone_factor_at_rc, each of which rebuilds and sorts an 18-point table on
    // five freshly heap-allocated vectors.  Everything in that loop except
    // cos^4(lat) and the pressure is common to every pixel, so it is hoisted to
    // a coarse LATITUDE lattice built once per run.
    //
    // WHY THE PRESSURE STAYS EXACT (the analogue of the LSDn surface's native
    // 380 pressure nodes, LsdnTimeAvgSurface above).  stone_factor_at_rc is
    // interp1d over a table whose BREAKPOINTS Rc_j are independent of P (they
    // are 14.3*cos^4(ilat_k), 0, and 16...25 GV) and whose VALUES S_j(P) depend
    // on P alone.  So
    //
    //     stone_factor_at_rc(Rc, P) = sum_j beta_j(Rc) * S_j(P)
    //
    // exactly, with P-independent weights beta.  Time-averaging is linear, so
    // averaging beta and applying S(P) afterwards is algebraically the same
    // number as averaging the factors themselves.  The surface therefore stores
    // the AVERAGED WEIGHTS, and pressure is applied per pixel at full accuracy.
    // The ONLY approximation in the object is the latitude interpolation, which
    // is what test L11 measures.
    static constexpr int stone_paleomag_basis_n = 18;

    // Shipped latitude resolution of the Stone-paleomag surface, CHOSEN BY
    // MEASUREMENT (L11 sweep, 2026-07-26), not assumed.  Measured max/RMS
    // relative deviation from the direct definition over off-node latitudes
    // 0-70 deg, pressures 500-1013 hPa and ages 2e4/5e4 yr:
    //
    //     h (deg)   max rel     rms rel
    //     1.0       1.21e-03    4.37e-04
    //     0.1       2.87e-06    7.41e-07
    //     0.025     4.25e-07    8.20e-08
    //     0.01      7.89e-08    1.52e-08     <-- shipped
    //
    // 0.01 deg is chosen because 7.6e-8 is at the precision of the float32
    // production raster this feeds (float32 eps 1.19e-7, half-ulp 5.96e-8), so
    // the surface is indistinguishable from the direct path in the actual
    // output; and because the Stone surface has NO pressure axis, so its build
    // is n_lat * n_record_steps and stays negligible even at this step.  That
    // is why it can afford to be 10x finer than the LSDn surface's 0.1 deg,
    // whose build carries a factor of 380 pressure nodes.
    //
    // Convergence is deliberately NOT assumed to be second order: the
    // correction is C0 but not C1 in latitude, because each record step drags
    // Rc_i(lat) across the Stone table's fixed breakpoints, so the error curve
    // has first-order plateaus (measured ratios per halving as low as 1.71).
    static constexpr double stone_paleomag_h_step_deg = 0.01;

    // THE DEFINITION.  Direct per-timestep trapezoidal average, no surface, no
    // interpolation in latitude - exactly the arithmetic this function has
    // always used.  Kept public and permanent: it is the oracle the surface is
    // measured against (L11/L14), and "the surface is faster" is never a reason
    // to delete it.  Returns 1.0 (not an error) for every degenerate case in
    // which the correction is not defined: nc not loaded, no record, age <= 0,
    // fewer than two record points within the age, or zero time span.
    static double stone_paleomag_correction_direct(double lat_deg, double P_mbar,
                                                   double age_yr,
                                                   const CrnNcData& nc);

    // Stone-paleomag flavour of the time-averaged surface.  Same struct, same
    // ownership and lifetime rules as the LSDn surface: built ONCE per raster
    // call before the OpenMP region, thereafter read-only and safe to share by
    // const& across threads; evaluation allocates nothing.
    //
    // The axis is |latitude| (Rc depends on cos^4(lat) only, so the surface is
    // symmetric about the equator and its lon_ax degenerates to a single node);
    // bounds that straddle the equator are widened to include 0.  Nodes are
    // snapped outward to multiples of h_step_deg, so node positions depend on
    // the step alone and not on the raster's exact extent - which is what lets
    // the scalar stone2000sp_paleomag() and the raster path agree.
    //
    // n_pixels is accepted for symmetry with build_lsdn_time_avg_surface but
    // gates nothing: with no pressure axis the build costs n_lat * n_steps -
    // about the work of two pixels - so there is no raster too small to repay
    // it, and no cost model to get wrong.
    //
    // Returns valid == false (having allocated nothing) with the decision in
    // .reason when the correction is not defined at all (the same degenerate
    // cases stone_paleomag_correction_direct returns 1.0 for).
    static LsdnTimeAvgSurface build_stone_paleomag_surface(
        const CrnNcData& nc, double age_yr,
        double lat_min_deg, double lat_max_deg,
        std::size_t n_pixels = 0,
        double h_step_deg = stone_paleomag_h_step_deg,
        bool force_build = false);

    // The paleomagnetic correction at one site from a prebuilt surface: linear
    // in |latitude| on the averaged and modern weight vectors, then the exact
    // pressure basis applied to both and the ratio taken.  Equal to
    // stone_paleomag_correction_direct() at a latitude node (to rounding) and
    // to within the L11 interpolation error between nodes.
    //
    // Pure arithmetic on const data - no allocation, no locking, no static
    // state - so it is safe inside an OpenMP parallel-for.  Queries outside the
    // surface's latitude extent are clamped to it.
    // Throws std::runtime_error if the surface is not valid; that is a caller
    // bug (validity belongs outside the parallel region).
    static double stone_paleomag_correction_from_surface(
        const LsdnTimeAvgSurface& surface, double lat_deg, double P_mbar);

    // The full surface-backed Stone factor: the time-independent Stone (2000)
    // factor times the correction above.  This is what the raster path calls.
    static double stone_paleomag_from_surface(const LsdnTimeAvgSurface& surface,
                                              double lat_deg, double P_mbar,
                                              double Fsp,
                                              const lsdtt_crn_parameters& params);

    //  LSDn scaling (Lifton, Sato & Dunai 2014)
    // Dimensionless production-rate scaling factor, normalised so that it is
    // exactly 1 at SLHL (Rc = 0 GV, P = 1013.25 hPa, S = S_ref = nc.lsdn_sphi[0]
    // = 624.5717522098033 MV):
    //
    //   lsdn_scaling = <LSDn_raw(Rc(t), P, S(t))>_t / LSDn_raw(0, 1013.25, S_ref)
    //
    // Drops into the same slot as stone2000sp(); P0, F[], Gamma[] and the muon
    // pathway are unaffected. Spallation scaling only - the muon pathways
    // receive it through the Vermeesch redistribution in scale_F_values(),
    // exactly as under Stone.
    //
    // nuclide : "Be10", "Al26" or "C14".
    // lat_deg, lon_deg : site position (deg; lon may be -180..180 or 0..360).
    //           LONGITUDE IS REQUIRED, not decorative: the Holocene geomagnetic
    //           latitude is the great-circle distance to a MOVING paleopole, so
    //           Rc(t) depends on site longitude as well as latitude. (The long
    //           record is GAD and uses |latitude| only.)
    // P_mbar  : atmospheric pressure (from any of the atmospheric_pressure_* fns).
    // age_yr  : <= 0 -> present-day field and solar modulation (record step 0),
    //                   two bilinear lookups, no time loop;
    //           >  0 -> duration-weighted mean over the Pavon-Carrasco field
    //                   record and the Lifton solar record, truncated at age_yr.
    //                   The records are defined on time-step BOUNDS that abut
    //                   exactly, so the average is duration-weighted over those
    //                   bounds - NOT the trapezoidal average over points that
    //                   stone2000sp_paleomag applies to the Muscheler record.
    //                   Allocation-free: accumulates scalars, never a per-call
    //                   std::vector.
    //
    // ASSUMPTION (documented, not silent): lsdtt3 has no sample-collection-year
    // input, so "now" is taken to be the records' reference year
    // nc.pavon_tzero = 2010 CE. The v3 calculator instead shifts the record by
    // (collection year - tzero); that shift is deliberately not applied here. A
    // sample_year parameter may be added later and will default to 2010.
    //
    // This is the REFERENCE path: correctness first. For raster work the caller
    // will later build a precomputed time-averaged surface once per run and
    // evaluate that per pixel; the two must agree, and this function is that
    // comparison's oracle as well as the clearest statement of the definition.
    //
    // Throws std::runtime_error if nc.lsdn_loaded is false or the nuclide has no
    // packed grid - LSDn is never silently substituted by Stone (deliberately
    // unlike the atmosphere selector, which may fall back quietly because it
    // refines one input rather than replacing the physics).
    double lsdn_scaling(const std::string& nuclide,
                        double lat_deg, double lon_deg, double P_mbar,
                        double age_yr, const CrnNcData& nc) const;

    // Raw grid value b(Rc,P) + S*m(Rc,P) for one nuclide, bilinearly
    // interpolated on the packed LSDn scaling grids.  This is NOT SLHL-
    // normalised - the SLHL normalisation lives in the scaling function that
    // consumes it, so that lsdn_scaling == 1 at (Rc = 0, P = 1013.25 hPa,
    // S = nc.lsdn_sphi[0]).
    //
    // nuclide : "Be10", "Al26" or "C14".
    // Rc_GV   : geomagnetic cutoff rigidity, clamped to [0, 20] GV (the grid's
    //           range; the same clamp the Balco v3 calculator applies).  No
    //           low-Rc floor is applied: the grid columns below 1 GV are
    //           already flat because the generator floored Rc there.
    // P_hPa   : atmospheric pressure, clamped to the grid range [300, 1058]
    //           hPa.  v3 would return NaN outside it; clamping is the right
    //           behaviour for a raster tool, where one NaN poisons a basin mean.
    // S_MV    : solar modulation potential.
    // Throws std::runtime_error if the nuclide has no packed grid.
    //
    // Grids and their axes are from the Balco v3 calculator data distribution
    // (LSDn2015_grid.mat, LSDn2015_N{10,26,14}quartz.mat), which precomputes the
    // Sato et al. (2008) PARMA spectral integration used by Lifton, Sato &
    // Dunai (2014) EPSL 386:149-160.
    static double lsdn_grid_value(const CrnNcData& nc, const std::string& nuclide,
                                  double Rc_GV, double P_hPa, double S_MV);

    // Lifton, Sato & Dunai (2014) EPSL 386:149-160, Eq. 2: geocentric dipole
    // effective vertical cutoff rigidity [GV].
    //   MM0        = M_t / M_0, dimensionless (M_0 = 7.746e22 A m^-^2, 2010 DGRF)
    //   gm_lat_rad = geomagnetic latitude [rad]
    // Returns |Rc|: the 6th-order polynomial in cos(lambda) leaks slightly
    // negative above ~80 deg, and Lifton et al. take the absolute value.
    // Deliberately NOT the Rc = Rc_c*vdm*cos^4(lat) form used by
    // stone2000sp_paleomag - the two differ by ~+12 to +24 % at mid to high
    // latitudes and must not be interchanged.
    static double lsd_dipole_rc(double MM0, double gm_lat_rad);

    // Geomagnetic latitude [rad] of a site given a (paleo)pole position, from
    // the great-circle angular distance theta to that pole:
    //   gm_lat = |pi/2 - theta|.
    // Used with the moving Pavon-Carrasco paleopole for the Holocene record, so
    // the geomagnetic latitude depends on site LONGITUDE as well as latitude.
    // The long (GAD) record instead uses |site latitude| directly.
    static double lsd_geomagnetic_lat(double lat_rad, double lon_rad,
                                      double pole_lat_rad, double pole_lon_rad);

    // Builds the LSDn time-averaged surface for a raster's actual footprint.
    // The lat/lon bounds are the raster's WGS84 extent; the caller should pad
    // them by about one cell so that every pixel is interior to the surface
    // (queries outside are clamped, not extrapolated).  The axes are snapped
    // outward to multiples of h_step_deg, so node positions are reproducible
    // and independent of the exact bounds.
    //
    // nuclide     : "Be10", "Al26" or "C14".
    // age_yr      : the sample age driving the time average, as in lsdn_scaling.
    // n_pixels    : the number of pixels the surface would serve; used only by
    //               the small-raster guard below.
    // h_step_deg  : horizontal resolution of the surface, in degrees.  ***THE
    //               DEFAULT IS A PLACEHOLDER*** carried over from the plan; it
    //               has NOT yet been chosen by the L11 sensitivity sweep.  Do
    //               not cite it as a validated value.
    // force_build : build even when the guard below would decline.  Exists so
    //               the equivalence test (L10) can exercise the surface path on
    //               small rasters where the direct path is cheaper - the whole
    //               point of that test is that it must not depend on which path
    //               the guard happens to pick.
    //
    // Returns valid == false (having allocated nothing) when a surface is not
    // worthwhile, with the decision recorded in .reason:
    //   - the loaded nc has no LSDn data;
    //   - age_yr <= 0                : the static path is already two lookups;
    //   - n_pixels < 512             : small-raster guard.  DELIBERATELY CRUDE,
    //                                  not a tuned cost model.  Both paths give
    //                                  the same answer (L10), so this only ever
    //                                  affects run time; past the crossover the
    //                                  surface wins by 490-2180x per pixel, and
    //                                  below it the whole stake is a few ms.
    //                                  The measured crossover is not stable
    //                                  enough to fit (190-522 px and 74-157 px
    //                                  in two honest runs), and every value in
    //                                  that range is far below any real basin -
    //                                  so this is a round number, not a fitted
    //                                  one, and needs no re-measurement when
    //                                  h_step_deg or the pressure axis changes.
    // The caller then uses lsdn_scaling() per pixel.  Both paths compute the
    // same quantity (test L10), so the choice is purely one of cost.
    static LsdnTimeAvgSurface build_lsdn_time_avg_surface(
        const CrnNcData& nc, const std::string& nuclide, double age_yr,
        double lat_min_deg, double lat_max_deg,
        double lon_min_deg, double lon_max_deg,
        std::size_t n_pixels, double h_step_deg = 0.1,
        bool force_build = false);

    // SLHL-normalised LSDn scaling at one site from a prebuilt surface:
    // trilinear in (lat, lon, P) - P innermost, then lon, then lat, so that a
    // query sitting exactly on a horizontal node reproduces the pure pressure
    // interpolation exactly - then divided by surface.ref.
    //
    // Equivalent to lsdn_scaling(surface.nuclide, lat, lon, P, surface.age_yr,
    // nc) to within the surface's horizontal interpolation error; exact (to
    // rounding) in pressure and at horizontal nodes.  This is test L10.
    //
    // Pure arithmetic on const data: no allocation, no locking, no static
    // state, so it is safe to call from inside an OpenMP parallel-for.
    // Queries outside the surface's horizontal extent are clamped to it - the
    // caller is expected to have built the surface with a pad.
    // Throws std::runtime_error if the surface is not valid; that is a caller
    // bug (the validity check belongs outside the parallel region).
    static double lsdn_scaling_from_surface(const LsdnTimeAvgSurface& surface,
                                            double lat_deg, double lon_deg,
                                            double P_mbar);

    //  Production scaling 
    // Multiplies P0 for all nuclides by the Stone 2000 factor (or any other
    // site-specific production scaling).  Call after reset_scaling(), before
    // scale_F_values().  Primarily used by the column / particle model.
    void apply_production_scaling(double factor)
    {
        P0_10Be_ *= factor;
        P0_26Al_ *= factor;
        P0_14C_  *= factor;
        P0_36Cl_ *= factor;
    }

    //  Standardisation 
    // Returns the factor to multiply a measured concentration by to convert it
    // to the reference standard (07KNSTD for Be10, KNSTD for Al26).
    double be10_standardisation_factor(const std::string& standard_name) const;
    double al26_standardisation_factor(const std::string& standard_name) const;

    // Applies standardisation to all concentrations in a vector.
    // standard_names: one per sample (empty string -> reference standard assumed).
    // nuclide: "Be10" or "Al26".
    std::vector<double> apply_standardisation(
        const std::vector<double>& concentrations,
        const std::vector<std::string>& standard_names,
        const std::string& nuclide) const;

    //  Getters (physics constants) 
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
    // Decay constants (yr^-^1)
    double lambda_10Be_, lambda_26Al_, lambda_14C_, lambda_36Cl_;

    // Sea-level high-latitude production rates (atoms g^-^1 yr^-^1)
    double P0_10Be_, P0_26Al_, P0_14C_, P0_36Cl_;

    // Attenuation lengths (g cm^-^2):
    //   [0] spallation  [1] slow neg. muons  [2] fast muons  [3] very fast muons
    double Gamma_[4];

    // Production fractions (dimensionless); modified in-place by scale_F_values
    double F_10Be_[4], F_26Al_[4], F_14C_[4];

    // Initialise to Vermeesh (2007) defaults
    void init_defaults();
    void init_standardisation_maps();

    // Internal Newton-Raphson: finds effective depth x such that
    //   Sigma_i F[i] * exp(-x/Gamma[i]) == total_shielding, then applies it to F[].
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

    // Bilinear interpolation on a REGULAR ascending grid, x fast and y slow, so
    // the flat index is iy*x_ax.size() + ix.  THIS AXIS ORDER IS BINDING: the
    // packed LSDn grids are stored (p, Rc) with p slow, so LSDn callers pass
    // (x_ax = lsdn_rc, y_ax = lsdn_p, x = Rc, y = P), and any later precomputed
    // surface must store and index its grids the same way.  A transposed lookup
    // is the silent failure mode here - it is caught by evaluating
    // b + 624.5718*m at the deliberately off-node (Rc = 0, p = 1013.25), which
    // must give 4.880884 for Be-10.
    // No search loop: the bracketing
    // index comes straight from (v - v0)/dv.  Distinct from bilinear_era40,
    // whose lat axis descends, whose lon axis wraps, and which searches.
    // Both axes are clamped to the grid bounds.
    // Returns NaN if either axis has fewer than two nodes.
    static double bilinear_regular(const std::vector<double>& x_ax,
                                   const std::vector<double>& y_ax,
                                   const std::vector<float>&  grid,
                                   double x, double y);

    // Stone spallation factor as a function of geomagnetic cutoff rigidity Rc [GV].
    // Builds the Rc->sf look-up table at pressure P_mbar using the Stone (2000)
    // coefficients and linearly interpolates.  Mirrors Rc_StoneFactors() in
    // PR_processes.py (Martin et al. / CREp).
    static double stone_factor_at_rc(double Rc_GV, double P_mbar);

    std::map<std::string, double> be10_std_map_;
    std::map<std::string, double> al26_std_map_;
};

} // namespace lsdtt
