// lsdtt_string_utils.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html
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

} // namespace lsdtt
