#pragma once
#ifndef LSDTT_CSV_TABLE_HPP
#define LSDTT_CSV_TABLE_HPP
// lsdtt_csv_table.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html
//
// -----------------------------------------------------------------------------
// A minimal header-keyed CSV reader for CLI input tables.
//
// This exists because several CLIs had each grown their own copy: read a header
// row, split the rest on commas, look fields up by name, and convert to double
// with a useful error message. The copies drifted -- two different whitespace
// trims, two different "is this field missing" rules -- which is exactly the
// kind of divergence that makes one tool reject a file another accepts.
//
// SCOPE, deliberately small. This is a SIMPLE splitter: it does not handle
// quoted fields containing commas, embedded newlines, or escape sequences. It is
// for the tabular parameter/sample inputs lsdtt3 CLIs read, which are machine
// generated or hand-written flat tables. If you need real RFC 4180 parsing, this
// is the wrong tool and should not be extended into one -- reach for a parser.
//
// For raster/point inputs that must be transformed into a raster CRS, use
// lsdtt_point_input instead; this is the plain-table path.
// -----------------------------------------------------------------------------

#include <string>
#include <unordered_map>
#include <vector>

namespace lsdtt
{

// One parsed row: field name -> raw (trimmed) text. Every header column is
// present as a key, so a short data row yields empty strings rather than
// missing keys, and callers never have to distinguish "absent" from "blank".
using lsdtt_csv_row = std::unordered_map<std::string, std::string>;

struct lsdtt_csv_table
{
    std::vector<std::string> headers;
    std::vector<lsdtt_csv_row> rows;
};

// Trim leading/trailing whitespace (space, tab, CR, LF). Returns "" for an
// all-whitespace input.
std::string csv_trim(const std::string& s);

// Split one line on commas. No quote handling -- see the scope note above.
std::vector<std::string> csv_split_row(const std::string& row);

// Throws std::runtime_error if the file cannot be opened, is empty, or has no
// header columns. Blank lines in the body are skipped.
lsdtt_csv_table load_csv_table(const std::string& path);

bool csv_has_column(const lsdtt_csv_table& t, const std::string& key);

// Missing key OR empty value both yield `fallback`: for these inputs a blank
// cell and an absent column mean the same thing to the caller.
std::string csv_row_string(const lsdtt_csv_row& row,
                           const std::string& key,
                           const std::string& fallback = "");

// Throws std::runtime_error naming the file, the 1-based row and the offending
// value if the field is missing, blank, or not numeric.
double csv_row_double_required(const lsdtt_csv_row& row,
                               const std::string& key,
                               const std::string& csv_path,
                               int row_index_1based);

// As above but a missing or blank field returns `fallback`. A field that is
// PRESENT but not numeric still throws: that is a typo, not an omission.
double csv_row_double_optional(const lsdtt_csv_row& row,
                               const std::string& key,
                               double fallback,
                               const std::string& csv_path,
                               int row_index_1based);

} // namespace lsdtt

#endif // LSDTT_CSV_TABLE_HPP
