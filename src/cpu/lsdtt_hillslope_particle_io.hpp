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
