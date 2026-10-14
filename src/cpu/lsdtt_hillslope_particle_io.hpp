// lsdtt_hillslope_particle_io.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html

#pragma once

#include <string>

#include "lsdtt_hillslope_cellular_cpu.hpp"

namespace lsdtt
{
void write_hillslope_particles_csv(
    const std::string& path,
    const lsdtt_hillslope_particle_tracker& particles);

lsdtt_hillslope_particle_tracker load_hillslope_particles_csv(
    const std::string& path);
}
