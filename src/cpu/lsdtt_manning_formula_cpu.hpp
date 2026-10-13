// lsdtt_manning_formula_cpu.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// Manning's equation for mean flow velocity (SI units):
//
//   v = (1/n) R^(2/3) S^(1/2)
//
// n is Manning's roughness (s m^-1/3), R the hydraulic radius (m) and S the
// energy slope as a dimensionless gradient (m/m). There are no guards: a caller
// whose inputs can be zero or negative checks them first, because what such a
// cell should get (zero, nodata, an error) differs between callers.
//
// Shared by lsdtt_hydro_velocity_field_cpu (which had it as a local function)
// and lsdtt_manning_hydraulics_cpu (consolidation C20). The inverse and
// wide-channel forms in lsdtt_hydro_roughness_cpu, lsdtt_hydro_routing_cpu and
// the generalised friction law in lsdtt_hydraulic_graph_cpu are different
// formulas and stay where they are.
// -----------------------------------------------------------------------------

#pragma once

#include <cmath>

namespace lsdtt
{

inline double manning_velocity(double roughness, double hydraulic_radius, double slope)
{
    return (1.0 / roughness) * std::pow(hydraulic_radius, 2.0 / 3.0) *
           std::sqrt(slope);
}

} // namespace lsdtt
