// lsdtt_hillslope_particle_io.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html

#include "lsdtt_hillslope_particle_io.hpp"
#include "../utils/lsdtt_output_file.hpp"
#include "../core/lsdtt_csv_table.hpp"

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
    static const std::string k_header =
        "id,cell,weight,source_id,age_yr,distance_cells,status,mixing_depth,seed,substep_index";
    std::string line;
    if (!std::getline(in, line) || line != k_header)
    {
        // A restart file edited or saved on Windows gains CRLF line endings;
        // say so, because the header otherwise looks right.
        const bool crlf = !line.empty() && line.back() == '\r' &&
                          line.compare(0, line.size() - 1, k_header) == 0;
        throw std::runtime_error("invalid particle restart header: " + path +
                                 (crlf ? " (the file has Windows CRLF line endings; "
                                         "convert it to LF)"
                                       : " (expected \"" + k_header + "\")"));
    }

    // Every field must be a number in full; a bad one is reported with the file,
    // line and text rather than as a bare "stod" / "stoull".
    const std::string where = "particle restart CSV " + path;
    lsdtt_hillslope_particle_tracker p;
    bool first = true;
    std::size_t line_no = 1;
    while (std::getline(in, line))
    {
        ++line_no;
        if (line.empty()) continue;
        std::stringstream ss(line);
        std::vector<std::string> v;
        std::string value;
        while (std::getline(ss, value, ',')) v.push_back(value);
        if (v.size() != 10)
            throw std::runtime_error("invalid particle restart row: " + path + ": line " +
                                     std::to_string(line_no) + " has " +
                                     std::to_string(v.size()) + " fields, expected 10");
        p.id.push_back(csv_parse_uint_field(v[0], where, line_no));
        p.cell.push_back(csv_parse_int_field(v[1], where, line_no));
        p.weight.push_back(csv_parse_double_field(v[2], where, line_no));
        p.source_id.push_back(static_cast<std::uint32_t>(csv_parse_uint_field(v[3], where, line_no)));
        p.age.push_back(csv_parse_double_field(v[4], where, line_no));
        p.distance_cells.push_back(csv_parse_double_field(v[5], where, line_no));
        p.status.push_back(static_cast<std::uint8_t>(csv_parse_uint_field(v[6], where, line_no)));
        const double mixing = csv_parse_double_field(v[7], where, line_no);
        const std::uint64_t seed = csv_parse_uint_field(v[8], where, line_no);
        const std::uint64_t step = csv_parse_uint_field(v[9], where, line_no);
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
