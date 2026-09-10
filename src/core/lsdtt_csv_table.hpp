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
// SCOPE, deliberately small. The splitter is quote-aware in the RFC 4180 sense
// -- a field may be wrapped in double quotes, a comma inside quotes is literal,
// and "" inside quotes is one literal quote -- but it is line-oriented and so
// does NOT support a newline embedded inside a quoted field: the reader has
// already cut the file into lines before the splitter sees it. There are no
// other escape sequences, no configurable delimiter and no type inference. It
// is for the tabular parameter/sample inputs lsdtt3 CLIs read, which are
// machine generated or hand-written flat tables. If you need multi-line fields
// or a full dialect-aware reader, this is the wrong tool -- reach for a parser
// rather than extending this one.
//
// It is also deliberately LENIENT. These files come from users, and three CLIs
// read them; refusing a file outright is a worse failure than a permissive
// parse, so no malformed input throws from the splitter. csv_split_row below
// states exactly what each malformed shape yields.
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

// Split one line into fields on commas, honouring double quotes. This does NOT
// trim: the caller does that (load_csv_table calls csv_trim on every token), so
// that a quoted field can keep leading/trailing spaces if a caller ever wants
// them. The contract, exactly:
//
//   a,b,c        -> {a, b, c}          unquoted fields are taken verbatim
//   a,"b",c      -> {a, b, c}          the wrapping quotes are removed
//   a,"b,c",d    -> {a, b,c, d}        a comma inside quotes is literal: 3 fields
//   a,"b""c",d   -> {a, b"c, d}        "" inside quotes is one literal quote
//   ,,           -> {"", "", ""}       empty fields are empty strings
//   a,b,         -> {a, b, ""}         a trailing comma yields a trailing field
//   ""           -> {""}               a quoted empty field is an empty string,
//                                      indistinguishable from an unquoted one
//   ""(the empty line) -> {}           an empty line has NO fields, not one
//    "a, b"      -> {a, b }            whitespace before the opening quote is
//                                      dropped with the quote; text after the
//                                      closing quote is kept, then trimmed by
//                                      the caller
//
// Malformed input, the lenient readings (none of these throw):
//
//   * UNTERMINATED QUOTE, e.g.  a,"b,c   -> {a, b,c}. End of line is treated as
//     an implicit closing quote, so the rest of the line is one field and the
//     commas in it stay literal. It is not re-parsed as unquoted text.
//   * QUOTE MID-FIELD in an otherwise unquoted field, e.g.  ab"cd  -> {ab"cd}.
//     Only a quote at the START of a field (after any leading whitespace) opens
//     a quoted field, so an interior quote is just a character. This is what
//     keeps every previously-working unquoted file byte-identical.
//   * TEXT AFTER A CLOSING QUOTE, e.g.  "ab"cd  -> {abcd}. The trailing run is
//     appended verbatim to the unquoted content rather than being dropped or
//     rejected. A trailing \r on a CRLF file lands here and in the plain
//     unquoted case, and is removed by csv_trim in the caller, not here.
//
// NOTE, behaviour change (decision D5): the previous getline-based splitter
// dropped a trailing empty field, so "a," gave {a}. It now gives {a, ""}, which
// means a header line ending in a comma now declares one extra, empty-named
// column.
std::vector<std::string> csv_split_row(const std::string& row);

// Split one line into fields AND trim each one -- csv_split_row followed by
// csv_trim on every field.
//
// This is what a caller reading a CSV by hand almost always wants, and it is
// what load_csv_table does internally. It exists as a named function because
// four CLIs had each grown their own splitter, and unifying them on
// csv_split_row alone simply moved the duplication: each then needed the same
// three-line trim loop. Reach for csv_split_row directly only when you have a
// reason to keep the surrounding whitespace on a field.
std::vector<std::string> csv_split_row_trimmed(const std::string& row);

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
