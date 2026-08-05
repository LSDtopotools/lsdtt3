// lsdtt_crn_units.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
//
// -----------------------------------------------------------------------------
// The unit conversions that sit between real geometry and the cosmogenic core.
//
// WHY THIS FILE EXISTS. The CRN physics works in MASS DEPTH, g/cm2. Metres
// appear only at the boundaries: DEM lowering, the Larsen area-volume scaling,
// burial depth of stored sediment. Every metre <-> g/cm2 conversion therefore
// carries a density -- and until 2026-08-05 the expression `x * rho * 0.1` was
// written out by hand at eight sites across the library and two CLIs, each
// independently deciding WHICH density to use.
//
// That is not a hypothetical risk. The density split of 2026-08-04 exists
// because a soil scar carried at the bedrock density (2650 rather than ~1500)
// is 1.77x too deep in mass, which changes spallation production at the scar
// base by a factor of 2.05. The arithmetic was never wrong; the density was.
//
// SO THE DENSITY ARGUMENT HAS NO DEFAULT, DELIBERATELY. A defaulted density
// would restore exactly the failure this file prevents: a call site that
// compiles, looks complete, and silently assumes bedrock. If you find yourself
// wanting a default, the question you actually need to answer is "what material
// is at this pixel", and the answer belongs at the call site.
//
// EXPRESSION ORDER IS PRESERVED, NOT TIDIED. Each function reproduces the
// multiplication order of the hand-written sites it replaced, so migrating to
// it is bit-for-bit. `0.1 * rho * m` and `m * rho * 0.1` are not the same
// double, and the refactor oracle compares hashes.
// -----------------------------------------------------------------------------

#ifndef LSDTT_CRN_UNITS_HPP
#define LSDTT_CRN_UNITS_HPP

namespace lsdtt
{

// A real thickness in metres -> mass depth in g/cm2.
//
//   1 kg/m3 * 1 m = 1 kg/m2 = 1000 g / 10^4 cm2 = 0.1 g/cm2
//
// `rho_kg_m3` is the density of the material ACTUALLY AT THAT PIXEL: bedrock
// (~2650) for an eroding rock column, soil (~1500) for a shallow landslide
// scar, alluvium (~2000-2200) for burial under channel fill.
constexpr double metres_to_g_cm2(double thickness_m, double rho_kg_m3)
{
    return thickness_m * rho_kg_m3 * 0.1;
}

// Mass depth in g/cm2 -> real thickness in metres. The exact inverse.
constexpr double g_cm2_to_metres(double mass_depth_g_cm2, double rho_kg_m3)
{
    return mass_depth_g_cm2 / (rho_kg_m3 * 0.1);
}

// A physical erosion rate in m/yr -> an effective erosion rate in g/cm2/yr.
// Same conversion as metres_to_g_cm2, named separately because a rate and a
// depth are different quantities and reading `metres_to_g_cm2(rate, rho)` at a
// call site invites the wrong review.
constexpr double m_yr_to_eff_erate(double rate_m_yr, double rho_kg_m3)
{
    return rate_m_yr * rho_kg_m3 * 0.1;
}

// Effective erosion rate in g/cm2/yr -> physical rate in m/yr.
constexpr double eff_erate_to_m_yr(double eff_g_cm2_yr, double rho_kg_m3)
{
    return eff_g_cm2_yr / (rho_kg_m3 * 0.1);
}

// mm/kyr <-> m/yr. No density: this is a pure change of units on a physical
// rate, and the two forms are both in common use in the literature.
constexpr double mm_kyr_to_m_yr(double rate_mm_kyr) { return rate_mm_kyr * 1.0e-6; }
constexpr double m_yr_to_mm_kyr(double rate_m_yr)   { return rate_m_yr   * 1.0e6;  }

// Effective erosion rate in g/cm2/yr -> mass flux in kg/m2/yr.
//
//   1 g/cm2 = 0.001 kg / 10^-4 m2 = 10 kg/m2
//
// No density: both sides are already masses per unit area, so this is a pure
// unit change. It is the conversion the delivery-mass weight uses, where the
// mixing weight must be commensurable between background erosion and a
// flushing landslide slab.
constexpr double eff_erate_to_kg_m2_yr(double eff_g_cm2_yr)
{
    return eff_g_cm2_yr * 10.0;
}

}  // namespace lsdtt

#endif  // LSDTT_CRN_UNITS_HPP
