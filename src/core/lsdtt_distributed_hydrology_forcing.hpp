#pragma once
#ifndef LSDTT_DISTRIBUTED_HYDROLOGY_FORCING_HPP
#define LSDTT_DISTRIBUTED_HYDROLOGY_FORCING_HPP
// lsdtt_distributed_hydrology_forcing.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// Readers for the per-timestep CSV forcing of lsdtt-distributed-hydrology, for
// both the TOPMODEL path and the HYPE column. One CSV, one row per timestep, in
// order, one named column holding the domain-wide value for that step.
//
// This began as a static function in the CLI, which reads rainfall and PET
// depths. It moved into the library because the HYPE column also needs
// TEMPERATURE forcing, and the depth reader could not be reused for it: it
// rejects negative values, which is right for a depth and wrong for a
// temperature. Rather than grow a second copy, both readers go through one
// internal routine and differ only in a sign policy (non-negative, or any
// finite value).
//
// It is deliberately named for the one component it serves. Cosmo, the LEM and
// graphflood each have forcings of their own, with their own units and rules;
// this is not a general forcing reader.
//
// Every failure throws std::runtime_error naming the file, and the row where
// there is one. The CLI is responsible for turning that into a FATAL line.
// -----------------------------------------------------------------------------

#include <string>
#include <vector>

namespace lsdtt
{

// A non-negative depth per timestep, in metres (rainfall, PET). Reads
// column_metres if it is present, otherwise column_mm multiplied by 0.001.
//
// If both columns are present column_metres wins and column_mm is ignored
// entirely -- it is not read, not checked and not compared. Negative or
// non-finite values, a file with neither column and a file with a header but
// no data rows are rejected.
//
// Two spellings are accepted because rain gauge and reanalysis products are
// published in millimetres while everything inside the model is in metres.
// Converting at the door is better than letting a user multiply by 1000 in a
// spreadsheet and get it wrong; there is no way to detect that mistake
// downstream, because a millimetre of rain read as a metre is a perfectly
// plausible storm.
std::vector<double> load_distributed_hydrology_depth_series(
    const std::string& path,
    const std::string& column_metres,
    const std::string& column_mm);

// Air temperature per timestep in degrees C (e.g. column "temperature_c"),
// read from one column with no unit conversion. Negative values are accepted.
// Non-finite values, a missing column and a file with a header but no data
// rows are rejected: a gap in a temperature record has to be filled
// deliberately by the user, not passed through into the snowpack.
std::vector<double> load_distributed_hydrology_temperature_series(
    const std::string& path,
    const std::string& column);

} // namespace lsdtt

#endif // LSDTT_DISTRIBUTED_HYDROLOGY_FORCING_HPP
