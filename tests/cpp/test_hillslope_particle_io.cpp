#include <filesystem>
#include <iostream>

#include "../../src/cpu/lsdtt_hillslope_particle_io.hpp"

int main(int argc, char** argv)
{
    const std::filesystem::path dir = argc > 1 ? argv[1] : ".";
    std::filesystem::create_directories(dir);
    const std::string path = (dir / "particle_restart_roundtrip.csv").string();
    lsdtt::lsdtt_hillslope_particle_tracker a;
    a.id = {4, 9}; a.cell = {12, 18}; a.weight = {1.25, 2.5};
    a.source_id = {3, 7}; a.age = {10.5, 20.25};
    a.distance_cells = {2.0, 8.5}; a.status = {0, 1};
    a.mixing_depth = 0.75; a.seed = 123456789; a.substep_index = 42;
    lsdtt::write_hillslope_particles_csv(path, a);
    const auto b = lsdtt::load_hillslope_particles_csv(path);
    const bool same = a.id == b.id && a.cell == b.cell && a.weight == b.weight
        && a.source_id == b.source_id && a.age == b.age
        && a.distance_cells == b.distance_cells && a.status == b.status
        && a.mixing_depth == b.mixing_depth && a.seed == b.seed
        && a.substep_index == b.substep_index;
    std::cout << (same ? "particle restart round-trip exact\n"
                       : "particle restart round-trip FAILED\n");
    return same ? 0 : 1;
}
