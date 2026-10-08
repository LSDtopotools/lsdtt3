// lsdtt_output_file.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// open_output_file(path, failure_prefix)
//
// Opens path for writing (truncating) and returns the stream. If it cannot be
// opened, throws std::runtime_error(failure_prefix + path).
//
// The caller supplies the whole message prefix, so each writer keeps the exact
// text it always had ("write_chi_network_fit_nodes_csv: cannot open ",
// "Cannot open CRONUS v3 input file: ", ...). It replaced the same
// open / test / throw lines in 22 writers.
// -----------------------------------------------------------------------------

#pragma once

#include <fstream>
#include <stdexcept>
#include <string>

namespace lsdtt
{

inline std::ofstream open_output_file(const std::string& path,
                                      const std::string& failure_prefix)
{
    std::ofstream out(path);
    if (!out)
    {
        throw std::runtime_error(failure_prefix + path);
    }
    return out;
}

} // namespace lsdtt
