// lsdtt_boundary_conditions.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// Parsing of the four-character boundary-condition string (north, east, south,
// west; codes such as o = open, n = no flux, p = periodic, b = base level).
//
// parse_boundary_codes(text, allowed, fold_case) checks, in this order:
//   1. the length is 4                          (status wrong_length)
//   2. each code, north to west, is in allowed  (status invalid_code; bad_index
//      is the first offender, and codes[] is filled up to and including it)
//   3. then reports whether north/south and east/west disagree on periodicity
//      (periodic_ns_mismatch, periodic_ew_mismatch).
// With fold_case each code is lower-cased (std::tolower) before it is checked
// and stored; allowed should then be lower case.
//
// It never throws: what a mismatch or a bad code MEANS differs between callers
// (lsdtt_flowinfo and lsdtt_flowgraph_cpu reject a periodic mismatch,
// lsdtt_depression_hierarchy makes the axis symmetric, each has its own
// exception type and text), so every caller keeps its own policy and messages.
// tests/cpp/test_boundary_conditions.cpp keeps the five replaced parsers as the
// reference.
//
// Header-only: string logic.
// -----------------------------------------------------------------------------

#pragma once

#include <array>
#include <cctype>
#include <cstddef>
#include <string>
#include <string_view>

namespace lsdtt
{

struct lsdtt_boundary_parse
{
    enum class status
    {
        ok,
        wrong_length,
        invalid_code
    };

    status              result               = status::ok;
    std::array<char, 4> codes                = {{0, 0, 0, 0}};
    std::size_t         bad_index            = 0;
    bool                periodic_ns_mismatch = false;
    bool                periodic_ew_mismatch = false;
};

inline lsdtt_boundary_parse parse_boundary_codes(const std::string& text,
                                                 std::string_view   allowed,
                                                 bool               fold_case)
{
    lsdtt_boundary_parse out;
    if (text.size() != 4)
    {
        out.result = lsdtt_boundary_parse::status::wrong_length;
        return out;
    }

    for (std::size_t i = 0; i < 4; ++i)
    {
        const char c = fold_case
            ? static_cast<char>(std::tolower(static_cast<unsigned char>(text[i])))
            : text[i];
        out.codes[i] = c;
        if (allowed.find(c) == std::string_view::npos)
        {
            out.result    = lsdtt_boundary_parse::status::invalid_code;
            out.bad_index = i;
            return out;
        }
    }

    out.periodic_ns_mismatch = (out.codes[0] == 'p') != (out.codes[2] == 'p');
    out.periodic_ew_mismatch = (out.codes[1] == 'p') != (out.codes[3] == 'p');
    return out;
}

} // namespace lsdtt
