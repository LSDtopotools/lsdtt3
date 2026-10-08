// test_boundary_conditions.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
//
// -----------------------------------------------------------------------------
// Characterisation test for parse_boundary_codes
// (src/utils/lsdtt_boundary_conditions.hpp, consolidation C09).
//
// Five boundary-condition parsers now call the shared parser and keep their own
// policy and messages. Each removed parser is kept below verbatim (member state
// passed as parameters) and compared with the migrated version on every
// four-character string over the alphabet "onpbONPBx" (6561 strings) plus the
// empty string and strings of length 3 and 5: the outcome (no throw, or the
// exception type and exact what()) and the parsed result must be identical.
//
// The sites that only test for 'p' (fill_pits, resampling, surface metrics,
// basin outline, synthetic surfaces) do not parse the string and are not
// migrated; the hillslope ADI parser is a stage 10 item.
//
// Usage: lsdtt-test-boundary-conditions
// Prints one line per check and a final "RESULT: PASS|FAIL".
// -----------------------------------------------------------------------------

#include <array>
#include <cctype>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "../../src/utils/lsdtt_boundary_conditions.hpp"

using namespace lsdtt;

namespace
{

// ============================================================================
// Removed parsers (references)
// ============================================================================

bool is_valid_boundary_code(char c)
{
    return c == 'o' || c == 'n' || c == 'p' || c == 'b';
}

// src/cpu/lsdtt_flowinfo_cpu.cpp:844-883 (as of 2026-10-08); boundary_conditions_
// is the std::string member (pre-sized to 4).
void old_flowinfo(const std::string& boundary_conditions, std::string& boundary_conditions_)
{
    if (boundary_conditions.size() != 4)
    {
        throw std::runtime_error(
            "lsdtt_flowinfo: boundary_conditions must have length 4 "
            "(north, east, south, west)"
        );
    }

    for (std::size_t i = 0; i < 4; ++i)
    {
        const char c = static_cast<char>(
            std::tolower(static_cast<unsigned char>(boundary_conditions[i]))
        );

        if (!is_valid_boundary_code(c))
        {
            throw std::runtime_error(
                "lsdtt_flowinfo: invalid boundary code '" +
                std::string(1, boundary_conditions[i]) + "'"
            );
        }

        boundary_conditions_[i] = c;
    }

    // Periodic boundaries must be symmetric in each axis.
    if ((boundary_conditions_[0] == 'p') != (boundary_conditions_[2] == 'p'))
    {
        throw std::runtime_error(
            "lsdtt_flowinfo: north and south periodic boundary flags must match"
        );
    }

    if ((boundary_conditions_[1] == 'p') != (boundary_conditions_[3] == 'p'))
    {
        throw std::runtime_error(
            "lsdtt_flowinfo: east and west periodic boundary flags must match"
        );
    }
}

// src/cpu/lsdtt_depression_hierarchy_cpu.cpp:63-94 (as of 2026-10-08)
std::array<char, 4> old_depression(const std::string& bc_str)
{
    if (bc_str.size() != 4)
    {
        throw std::runtime_error(
            "depression_hierarchy: boundary_conditions must be 4 characters "
            "(north, east, south, west)");
    }

    std::array<char, 4> bc{};
    for (int i = 0; i < 4; ++i)
    {
        const char c = static_cast<char>(
            std::tolower(static_cast<unsigned char>(
                bc_str[static_cast<std::size_t>(i)])));
        if (c != 'o' && c != 'b' && c != 'n' && c != 'p')
        {
            throw std::runtime_error(
                std::string("depression_hierarchy: invalid boundary code '") +
                c + "' (must be o, b, n or p)");
        }
        bc[static_cast<std::size_t>(i)] = c;
    }

    // Periodic boundaries must be symmetric on each axis.
    if (bc[0] == 'p') bc[2] = 'p';
    if (bc[2] == 'p') bc[0] = 'p';
    if (bc[1] == 'p') bc[3] = 'p';
    if (bc[3] == 'p') bc[1] = 'p';

    return bc;
}

struct boundary_flags
{
    bool north_periodic = false, north_forced = false;
    bool east_periodic = false,  east_forced = false;
    bool south_periodic = false, south_forced = false;
    bool west_periodic = false,  west_forced = false;
    bool operator==(const boundary_flags& o) const
    {
        return north_periodic == o.north_periodic && north_forced == o.north_forced &&
               east_periodic == o.east_periodic && east_forced == o.east_forced &&
               south_periodic == o.south_periodic && south_forced == o.south_forced &&
               west_periodic == o.west_periodic && west_forced == o.west_forced;
    }
};

// src/cpu/lsdtt_flowgraph_cpu.cpp:250-310 (as of 2026-10-08)
boundary_flags old_flowgraph(const std::string& bc)
{
    if (bc.size() != 4)
    {
        throw std::invalid_argument(
            "lsdtt_flowgraph_cpu: boundary_conditions must have exactly 4 characters "
            "(north, east, south, west)");
    }

    auto is_allowed = [](char c) -> bool
    {
        return c == 'o' || c == 'O' ||
               c == 'n' || c == 'N' ||
               c == 'p' || c == 'P' ||
               c == 'b' || c == 'B';
    };

    for (char c : bc)
    {
        if (!is_allowed(c))
        {
            throw std::invalid_argument(
                "lsdtt_flowgraph_cpu: boundary_conditions contains an invalid character. "
                "Allowed values are o, n, p, b.");
        }
    }

    boundary_flags f;
    auto is_periodic = [](char c) { return c == 'p' || c == 'P'; };
    auto is_forced   = [](char c) { return c == 'b' || c == 'B'; };

    f.north_periodic = is_periodic(bc[0]);
    f.north_forced   = is_forced(bc[0]);

    f.east_periodic  = is_periodic(bc[1]);
    f.east_forced    = is_forced(bc[1]);

    f.south_periodic = is_periodic(bc[2]);
    f.south_forced   = is_forced(bc[2]);

    f.west_periodic  = is_periodic(bc[3]);
    f.west_forced    = is_forced(bc[3]);

    if (f.north_periodic != f.south_periodic)
    {
        throw std::invalid_argument(
            "lsdtt_flowgraph_cpu: north and south periodic boundary flags "
            "must match");
    }

    if (f.east_periodic != f.west_periodic)
    {
        throw std::invalid_argument(
            "lsdtt_flowgraph_cpu: east and west periodic boundary flags "
            "must match");
    }

    return f;
}

// src/cpu/lsdtt_raster_flow_preprocessing_cpu.cpp:930-950 (as of 2026-10-08),
// fill_hydraulic_surface; returns {bn, be, bs, bw}.
std::array<char, 4> old_hydraulic(const std::string& boundary_conditions)
{
    if (boundary_conditions.size() != 4)
        throw std::invalid_argument(
            "fill_hydraulic_surface: boundary_conditions must contain NESW codes");

    const auto normalise_bc = [](char c) -> char
    {
        return static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    };
    const char bn = normalise_bc(boundary_conditions[0]);
    const char be = normalise_bc(boundary_conditions[1]);
    const char bs = normalise_bc(boundary_conditions[2]);
    const char bw = normalise_bc(boundary_conditions[3]);
    for (const char b : {bn, be, bs, bw})
    {
        if (b != 'o' && b != 'n' && b != 'p')
            throw std::invalid_argument(
                "fill_hydraulic_surface: hydraulic boundaries allow only o, n, or p");
    }
    if ((bn == 'p') != (bs == 'p') || (be == 'p') != (bw == 'p'))
        throw std::invalid_argument(
            "fill_hydraulic_surface: periodic boundaries must occur in opposite pairs");
    return {bn, be, bs, bw};
}

// src/cpu/lsdtt_graphflood_optimised_cpu.cpp:178-188 (as of 2026-10-08)
void old_graphflood(const std::string& bc)
{
    if (bc.size() != 4)
        throw std::invalid_argument("run_graphflood_optimised: boundary conditions must be NESW");

    for (const char code : bc)
        if (code != 'o' && code != 'n' && code != 'p')
            throw std::invalid_argument("run_graphflood_optimised: boundaries use o, n, or p");
    if ((bc[0] == 'p') != (bc[2] == 'p') || (bc[1] == 'p') != (bc[3] == 'p'))
        throw std::invalid_argument("run_graphflood_optimised: periodic boundaries must be paired");
}

// ============================================================================
// The migrated parsers, written exactly as they now appear in the sources
// ============================================================================

void new_flowinfo(const std::string& boundary_conditions, std::string& boundary_conditions_)
{
    const lsdtt_boundary_parse p = parse_boundary_codes(boundary_conditions, "onpb", true);
    if (p.result == lsdtt_boundary_parse::status::wrong_length)
    {
        throw std::runtime_error(
            "lsdtt_flowinfo: boundary_conditions must have length 4 "
            "(north, east, south, west)"
        );
    }
    if (p.result == lsdtt_boundary_parse::status::invalid_code)
    {
        // The codes before the bad one were stored before the throw, as always.
        for (std::size_t i = 0; i < p.bad_index; ++i) boundary_conditions_[i] = p.codes[i];
        throw std::runtime_error(
            "lsdtt_flowinfo: invalid boundary code '" +
            std::string(1, boundary_conditions[p.bad_index]) + "'"
        );
    }
    for (std::size_t i = 0; i < 4; ++i) boundary_conditions_[i] = p.codes[i];

    if (p.periodic_ns_mismatch)
    {
        throw std::runtime_error(
            "lsdtt_flowinfo: north and south periodic boundary flags must match"
        );
    }
    if (p.periodic_ew_mismatch)
    {
        throw std::runtime_error(
            "lsdtt_flowinfo: east and west periodic boundary flags must match"
        );
    }
}

std::array<char, 4> new_depression(const std::string& bc_str)
{
    const lsdtt_boundary_parse p = parse_boundary_codes(bc_str, "obnp", true);
    if (p.result == lsdtt_boundary_parse::status::wrong_length)
    {
        throw std::runtime_error(
            "depression_hierarchy: boundary_conditions must be 4 characters "
            "(north, east, south, west)");
    }
    if (p.result == lsdtt_boundary_parse::status::invalid_code)
    {
        throw std::runtime_error(
            std::string("depression_hierarchy: invalid boundary code '") +
            p.codes[p.bad_index] + "' (must be o, b, n or p)");
    }

    std::array<char, 4> bc = p.codes;

    // Periodic boundaries must be symmetric on each axis.
    if (bc[0] == 'p') bc[2] = 'p';
    if (bc[2] == 'p') bc[0] = 'p';
    if (bc[1] == 'p') bc[3] = 'p';
    if (bc[3] == 'p') bc[1] = 'p';

    return bc;
}

boundary_flags new_flowgraph(const std::string& bc)
{
    const lsdtt_boundary_parse p = parse_boundary_codes(bc, "onpb", true);
    if (p.result == lsdtt_boundary_parse::status::wrong_length)
    {
        throw std::invalid_argument(
            "lsdtt_flowgraph_cpu: boundary_conditions must have exactly 4 characters "
            "(north, east, south, west)");
    }
    if (p.result == lsdtt_boundary_parse::status::invalid_code)
    {
        throw std::invalid_argument(
            "lsdtt_flowgraph_cpu: boundary_conditions contains an invalid character. "
            "Allowed values are o, n, p, b.");
    }

    boundary_flags f;
    f.north_periodic = p.codes[0] == 'p';
    f.north_forced   = p.codes[0] == 'b';
    f.east_periodic  = p.codes[1] == 'p';
    f.east_forced    = p.codes[1] == 'b';
    f.south_periodic = p.codes[2] == 'p';
    f.south_forced   = p.codes[2] == 'b';
    f.west_periodic  = p.codes[3] == 'p';
    f.west_forced    = p.codes[3] == 'b';

    if (p.periodic_ns_mismatch)
    {
        throw std::invalid_argument(
            "lsdtt_flowgraph_cpu: north and south periodic boundary flags "
            "must match");
    }
    if (p.periodic_ew_mismatch)
    {
        throw std::invalid_argument(
            "lsdtt_flowgraph_cpu: east and west periodic boundary flags "
            "must match");
    }

    return f;
}

std::array<char, 4> new_hydraulic(const std::string& boundary_conditions)
{
    const lsdtt_boundary_parse p = parse_boundary_codes(boundary_conditions, "onp", true);
    if (p.result == lsdtt_boundary_parse::status::wrong_length)
        throw std::invalid_argument(
            "fill_hydraulic_surface: boundary_conditions must contain NESW codes");
    if (p.result == lsdtt_boundary_parse::status::invalid_code)
        throw std::invalid_argument(
            "fill_hydraulic_surface: hydraulic boundaries allow only o, n, or p");
    if (p.periodic_ns_mismatch || p.periodic_ew_mismatch)
        throw std::invalid_argument(
            "fill_hydraulic_surface: periodic boundaries must occur in opposite pairs");
    return p.codes;
}

void new_graphflood(const std::string& bc)
{
    const lsdtt_boundary_parse p = parse_boundary_codes(bc, "onp", false);
    if (p.result == lsdtt_boundary_parse::status::wrong_length)
        throw std::invalid_argument("run_graphflood_optimised: boundary conditions must be NESW");
    if (p.result == lsdtt_boundary_parse::status::invalid_code)
        throw std::invalid_argument("run_graphflood_optimised: boundaries use o, n, or p");
    if (p.periodic_ns_mismatch || p.periodic_ew_mismatch)
        throw std::invalid_argument("run_graphflood_optimised: periodic boundaries must be paired");
}

// ============================================================================
// Harness
// ============================================================================

int g_failures = 0;

void check(const std::string& what, bool ok, long long n, const std::string& detail = "")
{
    std::cout << (ok ? "PASS " : "FAIL ") << what << " (" << n << " strings)";
    if (!ok && !detail.empty()) std::cout << " -- " << detail;
    std::cout << "\n";
    if (!ok) ++g_failures;
}

// "none" or "<type>|<what>"
template <class F>
std::string outcome(F f)
{
    try { f(); return "none"; }
    catch (const std::invalid_argument& e) { return std::string("invalid_argument|") + e.what(); }
    catch (const std::runtime_error& e) { return std::string("runtime_error|") + e.what(); }
    catch (const std::exception& e) { return std::string("exception|") + e.what(); }
}

std::vector<std::string> all_inputs()
{
    const std::string alphabet = "onpbONPBx";
    std::vector<std::string> out = {"", "onp", "onpbo", "ONPB", "pppp", "pnpn"};
    for (char a : alphabet)
        for (char b : alphabet)
            for (char c : alphabet)
                for (char d : alphabet)
                    out.push_back(std::string{a, b, c, d});
    return out;
}

} // namespace

int main()
{
    std::cout << "== C09 boundary-condition parsers vs the removed copies ==\n";
    const std::vector<std::string> inputs = all_inputs();
    const long long n = static_cast<long long>(inputs.size());

    long long bad_fi = 0, bad_dh = 0, bad_fg = 0, bad_hy = 0, bad_gf = 0;
    std::string first;
    for (const std::string& s : inputs)
    {
        {
            std::string m_old = "????", m_new = "????";
            const std::string a = outcome([&] { old_flowinfo(s, m_old); });
            const std::string b = outcome([&] { new_flowinfo(s, m_new); });
            if (a != b || m_old != m_new) { ++bad_fi; if (first.empty()) first = "flowinfo '" + s + "'"; }
        }
        {
            std::array<char, 4> r_old{}, r_new{};
            const std::string a = outcome([&] { r_old = old_depression(s); });
            const std::string b = outcome([&] { r_new = new_depression(s); });
            if (a != b || r_old != r_new) { ++bad_dh; if (first.empty()) first = "depression '" + s + "'"; }
        }
        {
            boundary_flags r_old, r_new;
            const std::string a = outcome([&] { r_old = old_flowgraph(s); });
            const std::string b = outcome([&] { r_new = new_flowgraph(s); });
            if (a != b || !(r_old == r_new)) { ++bad_fg; if (first.empty()) first = "flowgraph '" + s + "'"; }
        }
        {
            std::array<char, 4> r_old{}, r_new{};
            const std::string a = outcome([&] { r_old = old_hydraulic(s); });
            const std::string b = outcome([&] { r_new = new_hydraulic(s); });
            if (a != b || r_old != r_new) { ++bad_hy; if (first.empty()) first = "hydraulic '" + s + "'"; }
        }
        {
            const std::string a = outcome([&] { old_graphflood(s); });
            const std::string b = outcome([&] { new_graphflood(s); });
            if (a != b) { ++bad_gf; if (first.empty()) first = "graphflood '" + s + "'"; }
        }
    }
    check("lsdtt_flowinfo: outcome, message and stored codes (incl. the prefix stored before a throw)", bad_fi == 0, n, first);
    check("depression_hierarchy: outcome, message and symmetrised codes", bad_dh == 0, n, first);
    check("lsdtt_flowgraph_cpu: outcome, message and flags", bad_fg == 0, n, first);
    check("fill_hydraulic_surface: outcome, message and codes", bad_hy == 0, n, first);
    check("graphflood_optimised: outcome and message (case-sensitive)", bad_gf == 0, n, first);

    std::cout << "RESULT: " << (g_failures == 0 ? "PASS" : "FAIL") << "\n";
    return g_failures == 0 ? 0 : 1;
}
