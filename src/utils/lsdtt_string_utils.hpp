// lsdtt_string_utils.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// Header-only string helpers shared by the CLIs and the library.
//
// Pure std::string manipulation only:
//   - no GDAL / OGR
//   - no file I/O
//   - no locale dependence
//
// This header exists because the same tiny case-folding helper had grown a
// private copy in 14 separate cli/*.cpp files, under six different names
// (to_lower_copy, to_lower, lower, lowercase_copy, plus bare inline
// std::transform loops). Parameter-string comparison is done everywhere in
// lsdtt3, so it needs exactly one implementation.
// -----------------------------------------------------------------------------

#pragma once

#include <cmath>
#include <sstream>
#include <cctype>
#include <string>

namespace lsdtt {

// -----------------------------------------------------------------------------
// to_lower_ascii
//
// Returns a lowercased copy of s. Takes the string BY VALUE and returns the
// modified copy, so a call on an rvalue moves rather than copies.
//
// ASCII-only and locale-independent by design:
//   - only the 26 bytes 'A'-'Z' are folded; every other byte is returned
//     unchanged. Digits, punctuation and whitespace pass through.
//   - bytes >= 0x80 (UTF-8 continuation bytes, Latin-1 accented characters)
//     pass through unchanged. Multi-byte UTF-8 sequences are therefore NOT
//     case-folded, which is the correct behaviour here: this is used to
//     normalise ASCII parameter keywords ("D8", "Same_As_Input", ...), not
//     human-language text.
//
// The std::tolower argument is cast to unsigned char first. Passing a plain
// char with the high bit set to std::tolower is undefined behaviour on any
// platform with a signed char, which includes x86-64 Linux; the cast is the
// whole reason this helper exists rather than a bare std::transform.
//
// Note that std::tolower consults the global C locale. For bytes < 0x80 in
// any locale lsdtt3 could plausibly run under, the mapping is the ASCII one,
// so the result is stable; we do not call std::setlocale anywhere.
// -----------------------------------------------------------------------------
inline std::string to_lower_ascii(std::string s)
{
    for (char& c : s)
    {
        c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    }
    return s;
}

// -----------------------------------------------------------------------------
// csv_number
//
// Formats a double for a CSV field. This is the single spelling of numeric
// output for every CSV lsdtt3 writes by hand; it replaces four private copies
// that had grown in cli/ (named number, format_double, fmt and precise_string)
// which between them used three precisions and three different non-finite
// spellings.
//
// Format: std::ostringstream in its default (defaultfloat) float format with
// precision(precision). That is exactly what all four copies did, so for
// finite values the output is byte-identical to what they produced at the same
// precision.
//
// NON-FINITE CONTRACT -- relied on by downstream readers:
//   NaN        -> "nan"
//   +infinity  -> "inf"
//   -infinity  -> "-inf"
// All lowercase. pandas.read_csv, numpy.loadtxt and R's read.csv all parse
// those three spellings natively, so a CSV written through this function reads
// back with its non-finite values intact and needs no na_values= argument.
// This also removes the platform-dependent "nan" / "-nan" that the two copies
// with no non-finite handling emitted by accident, since the sign of a NaN is
// not meaningful.
//
// PRECISION IS A PER-CALL-SITE PROPERTY, NOT A STYLE CHOICE. The default of 17
// is the round-trip precision for IEEE-754 double and is what new code should
// use. Existing call sites that were written at another precision must keep
// passing it explicitly (the swath-profile writers pass 9, the basin-extraction
// writers pass 10) so that their numeric output stays byte-identical to what
// they emitted before the unification. Do NOT "tidy" those call sites onto the
// default: doing so silently changes the contents of published CSVs.
// -----------------------------------------------------------------------------
inline std::string csv_number(double v, int precision = 17)
{
    if (std::isnan(v))     return "nan";
    if (!std::isfinite(v)) return (v > 0.0) ? "inf" : "-inf";

    std::ostringstream ss;
    ss.precision(precision);
    ss << v;
    return ss.str();
}

} // namespace lsdtt
