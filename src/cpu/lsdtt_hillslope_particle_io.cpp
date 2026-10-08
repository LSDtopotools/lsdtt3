#include "lsdtt_hillslope_particle_io.hpp"
#include "../utils/lsdtt_output_file.hpp"

#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <vector>

namespace lsdtt
{

void write_hillslope_particles_csv(
    const std::string& path,
    const lsdtt_hillslope_particle_tracker& p)
{
    std::ofstream out = open_output_file(
        path, "cannot write particle CSV: ");
    out << "id,cell,weight,source_id,age_yr,distance_cells,status,"
           "mixing_depth,seed,substep_index\n";
    out << std::setprecision(17);
    for (std::size_t i = 0; i < p.size(); ++i)
        out << p.id[i] << ',' << p.cell[i] << ',' << p.weight[i] << ','
            << p.source_id[i] << ',' << p.age[i] << ',' << p.distance_cells[i]
            << ',' << static_cast<int>(p.status[i]) << ',' << p.mixing_depth
            << ',' << p.seed << ',' << p.substep_index << '\n';
}

lsdtt_hillslope_particle_tracker load_hillslope_particles_csv(
    const std::string& path)
{
    std::ifstream in(path);
    if (!in) throw std::runtime_error("cannot open particle restart CSV: " + path);
    std::string line;
    if (!std::getline(in, line)
        || line != "id,cell,weight,source_id,age_yr,distance_cells,status,mixing_depth,seed,substep_index")
        throw std::runtime_error("invalid particle restart header: " + path);

    lsdtt_hillslope_particle_tracker p;
    bool first = true;
    while (std::getline(in, line))
    {
        if (line.empty()) continue;
        std::stringstream ss(line);
        std::vector<std::string> v;
        std::string value;
        while (std::getline(ss, value, ',')) v.push_back(value);
        if (v.size() != 10)
            throw std::runtime_error("invalid particle restart row: " + path);
        p.id.push_back(std::stoull(v[0]));
        p.cell.push_back(std::stoll(v[1]));
        p.weight.push_back(std::stod(v[2]));
        p.source_id.push_back(static_cast<std::uint32_t>(std::stoul(v[3])));
        p.age.push_back(std::stod(v[4]));
        p.distance_cells.push_back(std::stod(v[5]));
        p.status.push_back(static_cast<std::uint8_t>(std::stoul(v[6])));
        const double mixing = std::stod(v[7]);
        const std::uint64_t seed = std::stoull(v[8]);
        const std::uint64_t step = std::stoull(v[9]);
        if (first)
        {
            p.mixing_depth = mixing;
            p.seed = seed;
            p.substep_index = step;
            first = false;
        }
        else if (mixing != p.mixing_depth || seed != p.seed
                 || step != p.substep_index)
            throw std::runtime_error("inconsistent particle restart metadata: " + path);
    }
    if (p.size() == 0)
        throw std::runtime_error("particle restart contains no rows: " + path);
    return p;
}

} // namespace lsdtt
