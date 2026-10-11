// test_breach_least_cost.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
//
// -----------------------------------------------------------------------------
// Hermetic value tests for least-cost depression breaching
// (lsdtt::breach_depressions_least_cost; Lindsay 2020), against frozen
// outputs of WhiteboxTools next-gen BreachDepressionsLeastCost as shipped in
// whitebox_workflows 2.0.6 (run in memory, max_procs 1, fill_deps = False).
//
// Every fixture is built in-process from integer / dyadic values, and its
// float32 bits are checked against a digest of the input Whitebox was given,
// so a drifting builder is reported as such and not as a breach mismatch.
// The expected outputs are Whitebox's f64 values rounded to float32 (round to
// nearest even), stored as the list of cells Whitebox changed; every other
// valid cell must be unchanged. Comparison is BITWISE on every valid cell, and
// the nodata mask must be identical. No tolerance anywhere.
//
// Fixtures (in-process):
//   F0   single-cell spike on a ramp           pit raise and `small`
//   F1   bowl with one saddle                  basic path, carve values
//   F2   two sub-bowls sharing a flat rim      order, still-pit skip, nesting
//   F3   symmetric bowl, two equal saddles     tie-break (heap + neighbour order)
//   F4   diagonal vs cardinal saddle           step distance, minimize_dist
//   F5   troughs needing 10/11/12 steps        max_dist = path steps, +1 discovery
//   F6   pit with total cut exactly 2.0        max_cost on pop, strict >
//   F7   flat 5 x 5 floor                      in-place alternating pit raise
//   F8   flat 3 x 20 lake                      relief exhaustion, comb
//   F9   nodata neighbour / hole / ring        nodata is an outlet
//   F10  only low rim is a higher border cell  edge is NOT an outlet (Whitebox)
//   F11  two equal-z pits                      equal-z order: last in scan first
//   F12  F1 at floors 0, -50, 4000             auto increment, float32 rounding
//   F13  dyadic copies of F1/F3/F7             pure float32-exact arithmetic
//   F14  3x3, 1x9, all nodata, plane           degenerate inputs
//   F15  pit beside the east edge, high        edge policy: whitebox cuts the
//        divide to a lower interior basin      divide, outlet drains to the edge
//
// Runs per fixture: cell size {1, 30} m x minimize_dist {F, T} x increment
// {0.125, Whitebox auto}, max_dist 100; plus F5 at max_dist {10, 11, 12} and
// F6 at max_cost {1.0, 1.125, 1.25, 1.9, 2.0, 2.1}. "auto" runs pass
// Whitebox's own auto increment EXPLICITLY (its exact double value), so they
// test the core, not lsdtt3's auto-increment rule.
//
// Checks:
//   D0  each in-process fixture equals the input Whitebox was run on (FNV-1a
//       64 over the float32 bits, row-major).
//   W   per fixture: every run is bitwise equal to float32(Whitebox), with
//       edge_policy = whitebox.
//   C   per fixture: Whitebox's internal counters (pits detected / searched /
//       skipped / solved / unsolved) -- once the library exposes them.
//   E   max_dist < 1 and max_dist > 32767 are rejected with a message naming
//       max_dist; max_dist 32767 is accepted. fill_min_slope < 0 / NaN / inf
//       are rejected; barnes2014b_gradient with fill_min_slope 0 is rejected
//       unless fill_remainder is false. The edge-policy spelling (outlet,
//       whitebox) and its exact error message.
//   A   lsdtt3's automatic flat increment (flat_increment <= 0):
//       2 x ULP_f32(max(|z|max, 1)), pinned on F12's three elevation ranges,
//       flat rasters, a power-of-two maximum and an all-nodata raster.
//   R   with library defaults (auto increment, remainder fill at min_slope
//       1e-4) in BOTH edge policies, every interior cell (all 8 neighbours in
//       the raster and valid) has a strictly lower neighbour: F1, F7, F8,
//       F10, F15, cell size 1 and 30 m; also F8 at increment 0.125 (16
//       unsolved pits), and a control showing fill_min_slope 0 fails it.
//       The same with fill_flat_resolution barnes2014b_gradient and
//       barnes2014b_float_step (whitebox mode), and on F8 at 0.125 each
//       differs from inline_gradient (the option reaches the fill).
//   O   edge_policy = outlet (an off-grid neighbour is treated exactly as
//       Whitebox treats a nodata neighbour: the search ends there): default
//       is outlet; equal to Whitebox on every run whose Whitebox search never
//       expands a border cell; F10 solved through the border cell; F15 drains
//       to the east edge with a smaller maximum cut than Whitebox.
//
// API. The checks use the options overload
//
//   struct lsdtt::breach_least_cost_options {
//       int    max_dist;          // path length in steps (Whitebox semantics)
//       double max_cost;          // Whitebox units, pop-and-break
//       bool   minimize_dist;
//       double flat_increment;    // MUST be double: auto runs pass 2e-14 etc.
//       bool   fill_remainder;
//       lsdtt::breach_lc_edge_policy edge_policy;   // { outlet, whitebox }
//       float  fill_min_slope;    // remainder fill gradient (default 1e-4)
//       lsdtt::flat_resolution_choice fill_flat_resolution;  // inline_gradient
//   };
//   lsdtt_raster_cpu breach_depressions_least_cost(
//       const lsdtt_raster_cpu&, const breach_least_cost_options&,
//       breach_least_cost_diagnostics* = nullptr);
//
// (C reads the counters from breach_least_cost_diagnostics; E goes through
// the positional overload, which forwards to the options overload.) The two
// guards below (LSDTT_BREACH_LC_HAS_OPTIONS / _DIAGNOSTICS) select these; set
// to 0 they fall back to the positional API only, explicit-increment W runs,
// and NOT RUN for C and O.
//
// Real terrain (--real-terrain <crop.tif>, unit_breach_least_cost_real_terrain):
//   the committed 128 x 128 lakeland crop (tests/data/
//   flat_resolution_lakeland_crop_dem.tif).
//   RT0  the crop loads (128 x 128, 30 m, every cell valid) and its float32
//        values have the digest of the input Whitebox was run on.
//   RT1  matched parameters (edge_policy whitebox, fill_remainder false,
//        flat_increment float32(0.001) widened to double -- what the CLI
//        passes for 0.001 -- max_dist 100 and 1000, minimize_dist false and
//        true): the float32 output's FNV-1a-64 digest equals that of
//        float32(Whitebox 2.0.6's f64 output) on the same input, computed
//        from Whitebox itself (2026-10-11). The pit counters are pinned from
//        the lsdtt3 run that matched.
//   RT2  lsdtt3 defaults (edge_policy outlet, max_dist 1000, auto increment,
//        remainder fill on) and whitebox mode at its default max_dist 100;
//        and max_dist 3 (pits left
//        unsolved, so the remainder fill must run) with inline_gradient and
//        each barnes2014b option: every interior cell has a strictly lower
//        neighbour.
//   RT3  the same outputs: lsdtt_flowinfo (D8, boundary oooo) builds and
//        validate() passes.
//
// Usage: lsdtt-test-breach-least-cost [--real-terrain <crop.tif>]
// Prints one PASS/FAIL/NOT RUN line per check, a per-fixture table, and a
// final "RESULT: PASS|FAIL". Exit code 1 on any failure.
// -----------------------------------------------------------------------------

#ifndef LSDTT_BREACH_LC_HAS_OPTIONS
#define LSDTT_BREACH_LC_HAS_OPTIONS 1
#endif
#ifndef LSDTT_BREACH_LC_HAS_DIAGNOSTICS
#define LSDTT_BREACH_LC_HAS_DIAGNOSTICS 1
#endif

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <string>
#include <type_traits>
#include <vector>

#include "lsdtt_data_loader.hpp"
#include "lsdtt_flowinfo_cpu.hpp"
#include "lsdtt_raster_cpu.hpp"
#include "lsdtt_raster_flow_preprocessing_cpu.hpp"
#include "lsdtt_raster_breach_least_cost_cpu.hpp"

namespace
{

// =============================================================================
// Tiny check harness (same shape as tests/cpp/test_flow_properties.cpp)
// =============================================================================
int g_n_checks  = 0;
int g_n_failed  = 0;
int g_n_not_run = 0;

void check(const std::string& name, bool ok, const std::string& detail = "")
{
    ++g_n_checks;
    if (!ok) ++g_n_failed;
    std::cout << (ok ? "PASS" : "FAIL") << "  " << name;
    if (!ok && !detail.empty()) std::cout << "   [" << detail << "]";
    std::cout << "\n";
}

void not_run(const std::string& name, const std::string& why)
{
    ++g_n_not_run;
    std::cout << "NOT RUN  " << name << "   [" << why << "]\n";
}

// The current library prints progress lines on std::cout for every call;
// silence them so the PASS/FAIL lines stay readable.
struct cout_silencer
{
    std::streambuf* old_;
    std::ostringstream sink_;
    cout_silencer() : old_(std::cout.rdbuf(sink_.rdbuf())) {}
    ~cout_silencer() { std::cout.rdbuf(old_); }
};

constexpr double k_nodata = -9999.0;
constexpr double k_inf    = std::numeric_limits<double>::infinity();


// =============================================================================
// Literal tables (generated from the frozen Whitebox 2.0.6 characterisation;
// see the header). Do not edit by hand.
// =============================================================================
struct wb_change
{
    int   r;
    int   c;
    float z;          // float32(Whitebox f64 output)
};

struct wb_counters
{
    int pits_detected;
    int pits_searched;
    int pits_skipped_not_pit;
    int pits_solved;
    int pits_unsolved_exhausted;
    int pits_unsolved_max_cost;
};

struct wb_run
{
    const char*      fixture;
    const char*      key;
    double           cell_size;
    bool             minimize_dist;
    bool             explicit_increment;  // false: Whitebox auto (passed explicitly as `small`)
    double           small;               // Whitebox's flat increment, exact
    int              max_dist;
    double           max_cost;            // k_inf = unlimited
    wb_counters      counters;
    // Whether any Whitebox search in this run expanded a raster-border cell.
    // If not, the outlet edge policy cannot change the result.
    bool             border_popped;
    // Result of a scratch outlet-mode transliteration (informational; the O1
    // check selects on border_popped).
    bool             outlet_equal;
    const wb_change* changes;
    int              n_changes;
};

struct fixture_digest
{
    const char* name;
    int         rows;
    int         cols;
    uint64_t    fnv1a64;
};

// ---- BEGIN GENERATED LITERALS ----
const fixture_digest k_fixture_digests[] = {
    {"F0", 9, 9, 0xff0cd3b4f89a6ddaULL},
    {"F1", 11, 11, 0x54bdd5864a1eb7b8ULL},
    {"F2", 11, 19, 0xc72a2fc134533f78ULL},
    {"F3", 11, 11, 0x7fcf6c65499ff2d8ULL},
    {"F4", 9, 9, 0xcc7660415e5e6f08ULL},
    {"F5_s10", 7, 15, 0x2af839836d920745ULL},
    {"F5_d10", 15, 15, 0xbed8af7661ffaf35ULL},
    {"F5_s11", 7, 16, 0x63fc3034b347de95ULL},
    {"F5_d11", 16, 16, 0x6548db0b324794a5ULL},
    {"F5_s12", 7, 17, 0x008a4b6bcfa7bbfdULL},
    {"F5_d12", 17, 17, 0xd3b1c8d17ff36875ULL},
    {"F6", 7, 7, 0x6beafe08a34a5a23ULL},
    {"F7", 11, 11, 0x5ef8af65b87667b5ULL},
    {"F8", 7, 26, 0x5c0d264099a0b1abULL},
    {"F9a", 7, 7, 0x0ba7d3b6bc9e5952ULL},
    {"F9b", 9, 9, 0x22247506b4dd78b2ULL},
    {"F9c", 11, 11, 0x21a31810436d0698ULL},
    {"F10", 7, 7, 0xaae5bedebfe2431bULL},
    {"F11", 7, 13, 0x87a910e0101ad975ULL},
    {"F12_0", 11, 11, 0x4e534c72325e0875ULL},
    {"F12_m50", 11, 11, 0xc251f2449111f56fULL},
    {"F12_4000", 11, 11, 0xb74e7c15e2e2fa90ULL},
    {"F13_F1", 11, 11, 0x31bbc8a00f75ab28ULL},
    {"F13_F3", 11, 11, 0xdd1ab7038ee046f8ULL},
    {"F13_F7", 11, 11, 0x9adac44c39a5e5f7ULL},
    {"F14_3x3", 3, 3, 0x25bda83a8ea4e1f5ULL},
    {"F14_1x9", 1, 9, 0x29ac3e08bc21c828ULL},
    {"F14_allnd", 4, 4, 0xbfebaf0a75a08925ULL},
    {"F14_nopit", 6, 6, 0x6bc9cdf9664313f5ULL},
    {"F15", 9, 12, 0xc6398a39d573db0bULL},
};

const wb_change k_chg_000[] = {
    {4,4,0x1.9c00000000000p+3f}, {5,3,0x1.9800000000000p+3f}};
const wb_change k_chg_001[] = {
    {4,4,0x1.a000000000000p+3f}, {5,3,0x1.a000000000000p+3f}};
const wb_change k_chg_002[] = {
    {4,3,0x1.9800000000000p+3f}, {4,4,0x1.9c00000000000p+3f}};
const wb_change k_chg_003[] = {
    {4,3,0x1.a000000000000p+3f}, {4,4,0x1.a000000000000p+3f}};
const wb_change k_chg_004[] = {
    {5,5,0x1.e000000000000p+0f}, {5,8,0x1.8000000000000p+0f}, {6,6,0x1.c000000000000p+0f},
    {6,7,0x1.a000000000000p+0f}};
const wb_change k_chg_005[] = {
    {5,5,0x1.0000000000000p+1f}, {5,8,0x1.0000000000000p+1f}, {6,6,0x1.0000000000000p+1f},
    {6,7,0x1.0000000000000p+1f}};
const wb_change k_chg_006[] = {
    {5,5,0x1.e000000000000p+0f}, {5,6,0x1.c000000000000p+0f}, {5,7,0x1.a000000000000p+0f},
    {5,8,0x1.8000000000000p+0f}};
const wb_change k_chg_007[] = {
    {5,5,0x1.0000000000000p+1f}, {5,6,0x1.0000000000000p+1f}, {5,7,0x1.0000000000000p+1f},
    {5,8,0x1.0000000000000p+1f}};
const wb_change k_chg_008[] = {
    {1,6,0x1.e000000000000p+0f}, {2,6,0x1.0000000000000p+1f}, {3,6,0x1.1000000000000p+1f},
    {3,9,0x1.7000000000000p+1f}, {3,10,0x1.6000000000000p+1f}, {4,6,0x1.2000000000000p+1f},
    {4,11,0x1.5000000000000p+1f}, {5,5,0x1.3000000000000p+1f}, {5,9,0x1.7000000000000p+1f},
    {5,13,0x1.3000000000000p+1f}, {5,17,0x1.e000000000000p+0f}, {6,8,0x1.6000000000000p+1f},
    {6,14,0x1.2000000000000p+1f}, {6,16,0x1.0000000000000p+1f}, {7,7,0x1.5000000000000p+1f},
    {7,8,0x1.6000000000000p+1f}, {7,9,0x1.7000000000000p+1f}, {7,15,0x1.1000000000000p+1f}};
const wb_change k_chg_009[] = {
    {3,9,0x1.4000000000000p+1f}, {3,10,0x1.4000000000000p+1f}, {3,11,0x1.4000000000000p+1f},
    {4,6,0x1.4000000000000p+1f}, {4,8,0x1.4000000000000p+1f}, {4,12,0x1.4000000000000p+1f},
    {5,5,0x1.4000000000000p+1f}, {5,7,0x1.4000000000000p+1f}, {5,9,0x1.8000000000000p+1f},
    {5,13,0x1.4000000000000p+1f}, {5,17,0x1.4000000000000p+1f}, {6,10,0x1.8000000000000p+1f},
    {6,14,0x1.4000000000000p+1f}, {6,16,0x1.4000000000000p+1f}, {7,9,0x1.8000000000000p+1f},
    {7,15,0x1.4000000000000p+1f}};
const wb_change k_chg_010[] = {
    {3,7,0x1.5000000000000p+1f}, {3,8,0x1.6000000000000p+1f}, {3,9,0x1.7000000000000p+1f},
    {5,5,0x1.3000000000000p+1f}, {5,9,0x1.7000000000000p+1f}, {5,13,0x1.3000000000000p+1f},
    {5,14,0x1.2000000000000p+1f}, {5,15,0x1.1000000000000p+1f}, {5,16,0x1.0000000000000p+1f},
    {5,17,0x1.e000000000000p+0f}, {6,5,0x1.2000000000000p+1f}, {6,8,0x1.6000000000000p+1f},
    {7,5,0x1.1000000000000p+1f}, {7,7,0x1.5000000000000p+1f}, {7,8,0x1.6000000000000p+1f},
    {7,9,0x1.7000000000000p+1f}, {8,5,0x1.0000000000000p+1f}, {9,5,0x1.e000000000000p+0f}};
const wb_change k_chg_011[] = {
    {3,9,0x1.8000000000000p+1f}, {3,10,0x1.8000000000000p+1f}, {5,5,0x1.4000000000000p+1f},
    {5,6,0x1.4000000000000p+1f}, {5,7,0x1.4000000000000p+1f}, {5,8,0x1.4000000000000p+1f},
    {5,9,0x1.4000000000000p+1f}, {5,10,0x1.4000000000000p+1f}, {5,11,0x1.4000000000000p+1f},
    {5,12,0x1.4000000000000p+1f}, {5,13,0x1.4000000000000p+1f}, {6,5,0x1.4000000000000p+1f},
    {7,5,0x1.4000000000000p+1f}, {7,9,0x1.8000000000000p+1f}, {7,10,0x1.8000000000000p+1f},
    {8,5,0x1.4000000000000p+1f}, {9,5,0x1.4000000000000p+1f}};
const wb_change k_chg_012[] = {
    {3,5,0x1.3000000000000p+2f}, {4,4,0x1.3800000000000p+2f}};
const wb_change k_chg_013[] = {
    {3,5,0x1.4000000000000p+2f}, {4,4,0x1.4000000000000p+2f}};
const wb_change k_chg_014[] = {
    {4,4,0x1.3800000000000p+2f}, {4,5,0x1.3000000000000p+2f}};
const wb_change k_chg_015[] = {
    {4,4,0x1.4000000000000p+2f}, {4,5,0x1.4000000000000p+2f}};
const wb_change k_chg_016[] = {
    {3,2,0x1.4800000000000p+2f}, {3,3,0x1.4000000000000p+2f}, {3,4,0x1.3800000000000p+2f},
    {3,5,0x1.3000000000000p+2f}, {3,6,0x1.2800000000000p+2f}, {3,7,0x1.2000000000000p+2f},
    {3,8,0x1.1800000000000p+2f}, {3,9,0x1.1000000000000p+2f}, {3,10,0x1.0800000000000p+2f},
    {3,11,0x1.0000000000000p+2f}};
const wb_change k_chg_017[] = {
    {3,2,0x1.5000000000000p+2f}, {3,3,0x1.5000000000000p+2f}, {3,4,0x1.5000000000000p+2f},
    {3,5,0x1.5000000000000p+2f}, {3,6,0x1.5000000000000p+2f}, {3,7,0x1.5000000000000p+2f},
    {3,8,0x1.5000000000000p+2f}, {3,9,0x1.5000000000000p+2f}, {3,10,0x1.5000000000000p+2f},
    {3,11,0x1.5000000000000p+2f}};
const wb_change k_chg_018[] = {
    {2,2,0x1.4800000000000p+2f}, {3,3,0x1.4000000000000p+2f}, {4,4,0x1.3800000000000p+2f},
    {5,5,0x1.3000000000000p+2f}, {6,6,0x1.2800000000000p+2f}, {7,7,0x1.2000000000000p+2f},
    {8,8,0x1.1800000000000p+2f}, {9,9,0x1.1000000000000p+2f}, {10,10,0x1.0800000000000p+2f},
    {11,11,0x1.0000000000000p+2f}};
const wb_change k_chg_019[] = {
    {2,2,0x1.5000000000000p+2f}, {3,3,0x1.5000000000000p+2f}, {4,4,0x1.5000000000000p+2f},
    {5,5,0x1.5000000000000p+2f}, {6,6,0x1.5000000000000p+2f}, {7,7,0x1.5000000000000p+2f},
    {8,8,0x1.5000000000000p+2f}, {9,9,0x1.5000000000000p+2f}, {10,10,0x1.5000000000000p+2f},
    {11,11,0x1.5000000000000p+2f}};
const wb_change k_chg_020[] = {
    {3,2,0x1.4800000000000p+2f}, {3,3,0x1.4000000000000p+2f}, {3,4,0x1.3800000000000p+2f},
    {3,5,0x1.3000000000000p+2f}, {3,6,0x1.2800000000000p+2f}, {3,7,0x1.2000000000000p+2f},
    {3,8,0x1.1800000000000p+2f}, {3,9,0x1.1000000000000p+2f}, {3,10,0x1.0800000000000p+2f},
    {3,11,0x1.0000000000000p+2f}, {3,12,0x1.f000000000000p+1f}};
const wb_change k_chg_021[] = {
    {3,2,0x1.5000000000000p+2f}, {3,3,0x1.5000000000000p+2f}, {3,4,0x1.5000000000000p+2f},
    {3,5,0x1.5000000000000p+2f}, {3,6,0x1.5000000000000p+2f}, {3,7,0x1.5000000000000p+2f},
    {3,8,0x1.5000000000000p+2f}, {3,9,0x1.5000000000000p+2f}, {3,10,0x1.5000000000000p+2f},
    {3,11,0x1.5000000000000p+2f}, {3,12,0x1.5000000000000p+2f}};
const wb_change k_chg_022[] = {
    {2,2,0x1.4800000000000p+2f}, {3,3,0x1.4000000000000p+2f}, {4,4,0x1.3800000000000p+2f},
    {5,5,0x1.3000000000000p+2f}, {6,6,0x1.2800000000000p+2f}, {7,7,0x1.2000000000000p+2f},
    {8,8,0x1.1800000000000p+2f}, {9,9,0x1.1000000000000p+2f}, {10,10,0x1.0800000000000p+2f},
    {11,11,0x1.0000000000000p+2f}, {12,12,0x1.f000000000000p+1f}};
const wb_change k_chg_023[] = {
    {2,2,0x1.5000000000000p+2f}, {3,3,0x1.5000000000000p+2f}, {4,4,0x1.5000000000000p+2f},
    {5,5,0x1.5000000000000p+2f}, {6,6,0x1.5000000000000p+2f}, {7,7,0x1.5000000000000p+2f},
    {8,8,0x1.5000000000000p+2f}, {9,9,0x1.5000000000000p+2f}, {10,10,0x1.5000000000000p+2f},
    {11,11,0x1.5000000000000p+2f}, {12,12,0x1.5000000000000p+2f}};
const wb_change k_chg_024[] = {
    {3,2,0x1.4800000000000p+2f}};
const wb_change k_chg_025[] = {
    {3,2,0x1.4800000000000p+2f}, {3,3,0x1.4000000000000p+2f}, {3,4,0x1.3800000000000p+2f},
    {3,5,0x1.3000000000000p+2f}, {3,6,0x1.2800000000000p+2f}, {3,7,0x1.2000000000000p+2f},
    {3,8,0x1.1800000000000p+2f}, {3,9,0x1.1000000000000p+2f}, {3,10,0x1.0800000000000p+2f},
    {3,11,0x1.0000000000000p+2f}, {3,12,0x1.f000000000000p+1f}, {3,13,0x1.e000000000000p+1f}};
const wb_change k_chg_026[] = {
    {3,2,0x1.5000000000000p+2f}};
const wb_change k_chg_027[] = {
    {3,2,0x1.5000000000000p+2f}, {3,3,0x1.5000000000000p+2f}, {3,4,0x1.5000000000000p+2f},
    {3,5,0x1.5000000000000p+2f}, {3,6,0x1.5000000000000p+2f}, {3,7,0x1.5000000000000p+2f},
    {3,8,0x1.5000000000000p+2f}, {3,9,0x1.5000000000000p+2f}, {3,10,0x1.5000000000000p+2f},
    {3,11,0x1.5000000000000p+2f}, {3,12,0x1.5000000000000p+2f}, {3,13,0x1.5000000000000p+2f}};
const wb_change k_chg_028[] = {
    {2,2,0x1.4800000000000p+2f}};
const wb_change k_chg_029[] = {
    {2,2,0x1.4800000000000p+2f}, {3,3,0x1.4000000000000p+2f}, {4,4,0x1.3800000000000p+2f},
    {5,5,0x1.3000000000000p+2f}, {6,6,0x1.2800000000000p+2f}, {7,7,0x1.2000000000000p+2f},
    {8,8,0x1.1800000000000p+2f}, {9,9,0x1.1000000000000p+2f}, {10,10,0x1.0800000000000p+2f},
    {11,11,0x1.0000000000000p+2f}, {12,12,0x1.f000000000000p+1f}, {13,13,0x1.e000000000000p+1f}};
const wb_change k_chg_030[] = {
    {2,2,0x1.5000000000000p+2f}};
const wb_change k_chg_031[] = {
    {2,2,0x1.5000000000000p+2f}, {3,3,0x1.5000000000000p+2f}, {4,4,0x1.5000000000000p+2f},
    {5,5,0x1.5000000000000p+2f}, {6,6,0x1.5000000000000p+2f}, {7,7,0x1.5000000000000p+2f},
    {8,8,0x1.5000000000000p+2f}, {9,9,0x1.5000000000000p+2f}, {10,10,0x1.5000000000000p+2f},
    {11,11,0x1.5000000000000p+2f}, {12,12,0x1.5000000000000p+2f}, {13,13,0x1.5000000000000p+2f}};
const wb_change k_chg_032[] = {
    {3,3,0x1.c000000000000p-1f}, {3,4,0x1.8000000000000p-1f}, {3,5,0x1.4000000000000p-1f}};
const wb_change k_chg_033[] = {
    {3,3,0x1.0000000000000p+0f}, {3,4,0x1.0000000000000p+0f}, {3,5,0x1.0000000000000p+0f}};
const wb_change k_chg_034[] = {
    {3,3,0x1.c000000000000p-1f}};
const wb_change k_chg_035[] = {
    {3,3,0x1.c000000000000p-1f}, {3,5,0x1.c000000000000p-1f}, {3,7,0x1.c000000000000p-1f},
    {4,4,0x1.8000000000000p-1f}, {4,6,0x1.8000000000000p-1f}, {4,7,0x1.8000000000000p-1f},
    {5,3,0x1.c000000000000p-1f}, {5,5,0x1.4000000000000p-1f}, {5,7,0x1.4000000000000p-1f},
    {5,9,0x1.0000000000000p-3f}, {6,4,0x1.8000000000000p-1f}, {6,6,0x1.0000000000000p-1f},
    {6,8,0x1.0000000000000p-2f}, {7,3,0x1.c000000000000p-1f}, {7,5,0x1.4000000000000p-1f},
    {7,7,0x1.8000000000000p-2f}};
const wb_change k_chg_036[] = {
    {3,3,0x1.0000000000000p+0f}, {3,5,0x1.0000000000000p+0f}, {3,7,0x1.0000000000000p+0f},
    {4,4,0x1.0000000000000p+0f}, {4,6,0x1.0000000000000p+0f}, {4,7,0x1.0000000000000p+0f},
    {5,3,0x1.0000000000000p+0f}, {5,5,0x1.0000000000000p+0f}, {5,7,0x1.0000000000000p+0f},
    {6,4,0x1.0000000000000p+0f}, {6,6,0x1.0000000000000p+0f}, {6,8,0x1.0000000000000p+0f},
    {7,3,0x1.0000000000000p+0f}, {7,5,0x1.0000000000000p+0f}, {7,7,0x1.0000000000000p+0f}};
const wb_change k_chg_037[] = {
    {3,3,0x1.c000000000000p-1f}, {3,4,0x1.8000000000000p-1f}, {3,5,0x1.4000000000000p-1f},
    {3,7,0x1.c000000000000p-1f}, {4,6,0x1.0000000000000p-1f}, {4,7,0x1.8000000000000p-1f},
    {5,3,0x1.c000000000000p-1f}, {5,5,0x1.c000000000000p-1f}, {5,7,0x1.8000000000000p-2f},
    {5,9,0x1.0000000000000p-3f}, {6,4,0x1.8000000000000p-1f}, {6,5,0x1.8000000000000p-1f},
    {6,8,0x1.0000000000000p-2f}, {7,3,0x1.c000000000000p-1f}, {7,4,0x1.8000000000000p-1f},
    {7,5,0x1.4000000000000p-1f}, {7,6,0x1.0000000000000p-1f}, {7,7,0x1.8000000000000p-2f}};
const wb_change k_chg_038[] = {
    {3,3,0x1.0000000000000p+0f}, {3,4,0x1.0000000000000p+0f}, {3,5,0x1.0000000000000p+0f},
    {3,7,0x1.0000000000000p+0f}, {4,6,0x1.0000000000000p+0f}, {4,7,0x1.0000000000000p+0f},
    {5,3,0x1.0000000000000p+0f}, {5,5,0x1.0000000000000p+0f}, {5,7,0x1.0000000000000p+0f},
    {6,4,0x1.0000000000000p+0f}, {6,5,0x1.0000000000000p+0f}, {6,8,0x1.0000000000000p+0f},
    {7,3,0x1.0000000000000p+0f}, {7,4,0x1.0000000000000p+0f}, {7,5,0x1.0000000000000p+0f},
    {7,6,0x1.0000000000000p+0f}, {7,7,0x1.0000000000000p+0f}};
const wb_change k_chg_039[] = {
    {2,2,0x1.c000000000000p-1f}, {2,4,0x1.c000000000000p-1f}, {2,6,0x1.c000000000000p-1f},
    {2,8,0x1.c000000000000p-1f}, {2,10,0x1.c000000000000p-1f}, {2,12,0x1.c000000000000p-1f},
    {2,14,0x1.c000000000000p-1f}, {2,16,0x1.c000000000000p-1f}, {2,18,0x1.c000000000000p-1f},
    {2,20,0x1.c000000000000p-1f}, {3,19,0x1.8000000000000p-1f}, {3,21,0x1.0000000000000p-1f},
    {3,22,0x1.8000000000000p-2f}, {3,23,0x1.0000000000000p-2f}, {3,24,0x1.0000000000000p-3f},
    {4,2,0x1.c000000000000p-1f}, {4,4,0x1.c000000000000p-1f}, {4,6,0x1.c000000000000p-1f},
    {4,8,0x1.c000000000000p-1f}, {4,10,0x1.c000000000000p-1f}, {4,12,0x1.c000000000000p-1f},
    {4,14,0x1.c000000000000p-1f}, {4,16,0x1.c000000000000p-1f}, {4,18,0x1.c000000000000p-1f},
    {4,20,0x1.4000000000000p-1f}};
const wb_change k_chg_040[] = {
    {2,2,0x1.0000000000000p+0f}, {2,4,0x1.0000000000000p+0f}, {2,6,0x1.0000000000000p+0f},
    {2,8,0x1.0000000000000p+0f}, {2,10,0x1.0000000000000p+0f}, {2,12,0x1.0000000000000p+0f},
    {2,14,0x1.0000000000000p+0f}, {2,16,0x1.0000000000000p+0f}, {2,18,0x1.0000000000000p+0f},
    {2,20,0x1.0000000000000p+0f}, {3,3,0x1.0000000000000p+0f}, {3,5,0x1.0000000000000p+0f},
    {3,7,0x1.0000000000000p+0f}, {3,9,0x1.0000000000000p+0f}, {3,11,0x1.0000000000000p+0f},
    {3,13,0x1.0000000000000p+0f}, {3,15,0x1.0000000000000p+0f}, {3,17,0x1.0000000000000p+0f},
    {3,19,0x1.0000000000000p+0f}, {3,21,0x1.0000000000000p+0f}, {3,22,0x1.0000000000000p+0f},
    {4,2,0x1.0000000000000p+0f}, {4,4,0x1.0000000000000p+0f}, {4,6,0x1.0000000000000p+0f},
    {4,8,0x1.0000000000000p+0f}, {4,10,0x1.0000000000000p+0f}, {4,12,0x1.0000000000000p+0f},
    {4,14,0x1.0000000000000p+0f}, {4,16,0x1.0000000000000p+0f}, {4,18,0x1.0000000000000p+0f},
    {4,20,0x1.0000000000000p+0f}};
const wb_change k_chg_041[] = {
    {2,2,0x1.c000000000000p-1f}, {2,4,0x1.c000000000000p-1f}, {2,6,0x1.c000000000000p-1f},
    {2,8,0x1.c000000000000p-1f}, {2,10,0x1.c000000000000p-1f}, {2,12,0x1.c000000000000p-1f},
    {2,14,0x1.c000000000000p-1f}, {2,16,0x1.c000000000000p-1f}, {2,18,0x1.c000000000000p-1f},
    {2,20,0x1.c000000000000p-1f}, {2,21,0x1.8000000000000p-1f}, {3,19,0x1.8000000000000p-1f},
    {3,22,0x1.8000000000000p-2f}, {3,23,0x1.0000000000000p-2f}, {3,24,0x1.0000000000000p-3f},
    {4,2,0x1.c000000000000p-1f}, {4,4,0x1.c000000000000p-1f}, {4,6,0x1.c000000000000p-1f},
    {4,8,0x1.c000000000000p-1f}, {4,10,0x1.c000000000000p-1f}, {4,12,0x1.c000000000000p-1f},
    {4,14,0x1.c000000000000p-1f}, {4,16,0x1.c000000000000p-1f}, {4,18,0x1.c000000000000p-1f},
    {4,19,0x1.8000000000000p-1f}, {4,20,0x1.4000000000000p-1f}, {4,21,0x1.0000000000000p-1f}};
const wb_change k_chg_042[] = {
    {2,2,0x1.0000000000000p+0f}, {2,4,0x1.0000000000000p+0f}, {2,6,0x1.0000000000000p+0f},
    {2,8,0x1.0000000000000p+0f}, {2,10,0x1.0000000000000p+0f}, {2,12,0x1.0000000000000p+0f},
    {2,14,0x1.0000000000000p+0f}, {2,16,0x1.0000000000000p+0f}, {2,18,0x1.0000000000000p+0f},
    {2,20,0x1.0000000000000p+0f}, {2,21,0x1.0000000000000p+0f}, {3,3,0x1.0000000000000p+0f},
    {3,4,0x1.0000000000000p+0f}, {3,6,0x1.0000000000000p+0f}, {3,8,0x1.0000000000000p+0f},
    {3,10,0x1.0000000000000p+0f}, {3,12,0x1.0000000000000p+0f}, {3,14,0x1.0000000000000p+0f},
    {3,16,0x1.0000000000000p+0f}, {3,18,0x1.0000000000000p+0f}, {3,22,0x1.0000000000000p+0f},
    {4,2,0x1.0000000000000p+0f}, {4,3,0x1.0000000000000p+0f}, {4,4,0x1.0000000000000p+0f},
    {4,5,0x1.0000000000000p+0f}, {4,6,0x1.0000000000000p+0f}, {4,7,0x1.0000000000000p+0f},
    {4,8,0x1.0000000000000p+0f}, {4,9,0x1.0000000000000p+0f}, {4,10,0x1.0000000000000p+0f},
    {4,11,0x1.0000000000000p+0f}, {4,12,0x1.0000000000000p+0f}, {4,13,0x1.0000000000000p+0f},
    {4,14,0x1.0000000000000p+0f}, {4,15,0x1.0000000000000p+0f}, {4,16,0x1.0000000000000p+0f},
    {4,17,0x1.0000000000000p+0f}, {4,18,0x1.0000000000000p+0f}, {4,19,0x1.0000000000000p+0f},
    {4,20,0x1.0000000000000p+0f}, {4,21,0x1.0000000000000p+0f}};
const wb_change k_chg_043[] = {
    {-1,-1,0.0f}};
const wb_change k_chg_044[] = {
    {4,4,0x1.3800000000000p+2f}, {5,5,0x1.3000000000000p+2f}, {5,6,0x1.2800000000000p+2f}};
const wb_change k_chg_045[] = {
    {4,4,0x1.4000000000000p+2f}, {5,5,0x1.4000000000000p+2f}, {5,6,0x1.4000000000000p+2f}};
const wb_change k_chg_046[] = {
    {4,4,0x1.3800000000000p+2f}, {4,5,0x1.3000000000000p+2f}, {4,6,0x1.2800000000000p+2f}};
const wb_change k_chg_047[] = {
    {4,4,0x1.4000000000000p+2f}, {4,5,0x1.4000000000000p+2f}, {4,6,0x1.4000000000000p+2f}};
const wb_change k_chg_048[] = {
    {2,6,0x1.4000000000000p+1f}, {3,6,0x1.5000000000000p+1f}, {4,6,0x1.6000000000000p+1f},
    {5,5,0x1.7000000000000p+1f}};
const wb_change k_chg_049[] = {
    {2,6,0x1.8000000000000p+1f}, {3,6,0x1.8000000000000p+1f}, {4,6,0x1.8000000000000p+1f},
    {5,5,0x1.8000000000000p+1f}};
const wb_change k_chg_050[] = {
    {2,5,0x1.4000000000000p+1f}, {3,5,0x1.5000000000000p+1f}, {4,5,0x1.6000000000000p+1f},
    {5,5,0x1.7000000000000p+1f}};
const wb_change k_chg_051[] = {
    {2,5,0x1.8000000000000p+1f}, {3,5,0x1.8000000000000p+1f}, {4,5,0x1.8000000000000p+1f},
    {5,5,0x1.8000000000000p+1f}};
const wb_change k_chg_052[] = {
    {1,3,0x1.e000000000000p+0f}};
const wb_change k_chg_053[] = {
    {1,3,0x1.0000000000000p+1f}};
const wb_change k_chg_054[] = {
    {3,3,0x1.7000000000000p+1f}, {3,9,0x1.1000000000000p+1f}, {4,4,0x1.6000000000000p+1f},
    {4,8,0x1.2000000000000p+1f}, {4,10,0x1.0000000000000p+1f}, {4,11,0x1.e000000000000p+0f},
    {5,5,0x1.5000000000000p+1f}, {5,7,0x1.3000000000000p+1f}, {6,6,0x1.4000000000000p+1f}};
const wb_change k_chg_055[] = {
    {3,3,0x1.8000000000000p+1f}, {3,9,0x1.8000000000000p+1f}, {4,4,0x1.8000000000000p+1f},
    {4,8,0x1.8000000000000p+1f}, {4,10,0x1.8000000000000p+1f}, {4,11,0x1.8000000000000p+1f},
    {5,5,0x1.8000000000000p+1f}, {5,7,0x1.8000000000000p+1f}, {6,6,0x1.8000000000000p+1f}};
const wb_change k_chg_056[] = {
    {3,3,0x1.7000000000000p+1f}, {3,4,0x1.6000000000000p+1f}, {3,5,0x1.5000000000000p+1f},
    {3,6,0x1.4000000000000p+1f}, {3,7,0x1.3000000000000p+1f}, {3,8,0x1.2000000000000p+1f},
    {3,9,0x1.1000000000000p+1f}, {3,10,0x1.0000000000000p+1f}, {3,11,0x1.e000000000000p+0f}};
const wb_change k_chg_057[] = {
    {3,3,0x1.8000000000000p+1f}, {3,4,0x1.8000000000000p+1f}, {3,5,0x1.8000000000000p+1f},
    {3,6,0x1.8000000000000p+1f}, {3,7,0x1.8000000000000p+1f}, {3,8,0x1.8000000000000p+1f},
    {3,9,0x1.8000000000000p+1f}, {3,10,0x1.8000000000000p+1f}, {3,11,0x1.8000000000000p+1f}};
const wb_change k_chg_058[] = {
    {5,5,0x1.c000000000000p-1f}, {5,8,0x1.0000000000000p-1f}, {6,6,0x1.8000000000000p-1f},
    {6,7,0x1.4000000000000p-1f}};
const wb_change k_chg_059[] = {
    {5,5,0x1.0000000000000p+0f}, {5,8,0x1.0000000000000p+0f}, {6,6,0x1.0000000000000p+0f},
    {6,7,0x1.0000000000000p+0f}};
const wb_change k_chg_060[] = {
    {5,5,0x1.c000000000000p-1f}, {5,6,0x1.8000000000000p-1f}, {5,7,0x1.4000000000000p-1f},
    {5,8,0x1.0000000000000p-1f}};
const wb_change k_chg_061[] = {
    {5,5,0x1.0000000000000p+0f}, {5,6,0x1.0000000000000p+0f}, {5,7,0x1.0000000000000p+0f},
    {5,8,0x1.0000000000000p+0f}};
const wb_change k_chg_062[] = {
    {5,5,-0x1.8900000000000p+5f}, {5,8,-0x1.8c00000000000p+5f}, {6,6,-0x1.8a00000000000p+5f},
    {6,7,-0x1.8b00000000000p+5f}};
const wb_change k_chg_063[] = {
    {5,5,-0x1.8800000000000p+5f}, {5,8,-0x1.8800000000000p+5f}, {6,6,-0x1.8800000000000p+5f},
    {6,7,-0x1.8800000000000p+5f}};
const wb_change k_chg_064[] = {
    {5,5,-0x1.8900000000000p+5f}, {5,6,-0x1.8a00000000000p+5f}, {5,7,-0x1.8b00000000000p+5f},
    {5,8,-0x1.8c00000000000p+5f}};
const wb_change k_chg_065[] = {
    {5,5,-0x1.8800000000000p+5f}, {5,6,-0x1.8800000000000p+5f}, {5,7,-0x1.8800000000000p+5f},
    {5,8,-0x1.8800000000000p+5f}};
const wb_change k_chg_066[] = {
    {5,5,0x1.f41c000000000p+11f}, {5,8,0x1.f410000000000p+11f}, {6,6,0x1.f418000000000p+11f},
    {6,7,0x1.f414000000000p+11f}};
const wb_change k_chg_067[] = {
    {5,5,0x1.f420000000000p+11f}, {5,8,0x1.f420000000000p+11f}, {6,6,0x1.f420000000000p+11f},
    {6,7,0x1.f420000000000p+11f}};
const wb_change k_chg_068[] = {
    {5,5,0x1.f41c000000000p+11f}, {5,6,0x1.f418000000000p+11f}, {5,7,0x1.f414000000000p+11f},
    {5,8,0x1.f410000000000p+11f}};
const wb_change k_chg_069[] = {
    {5,5,0x1.f420000000000p+11f}, {5,6,0x1.f420000000000p+11f}, {5,7,0x1.f420000000000p+11f},
    {5,8,0x1.f420000000000p+11f}};
const wb_change k_chg_070[] = {
    {5,5,0x1.4000000000000p+1f}, {5,8,0x1.1000000000000p+1f}, {6,6,0x1.3000000000000p+1f},
    {6,7,0x1.2000000000000p+1f}};
const wb_change k_chg_071[] = {
    {5,5,0x1.5000000000000p+1f}, {5,8,0x1.5000000000000p+1f}, {6,6,0x1.5000000000000p+1f},
    {6,7,0x1.5000000000000p+1f}};
const wb_change k_chg_072[] = {
    {5,5,0x1.4000000000000p+1f}, {5,6,0x1.3000000000000p+1f}, {5,7,0x1.2000000000000p+1f},
    {5,8,0x1.1000000000000p+1f}};
const wb_change k_chg_073[] = {
    {5,5,0x1.5000000000000p+1f}, {5,6,0x1.5000000000000p+1f}, {5,7,0x1.5000000000000p+1f},
    {5,8,0x1.5000000000000p+1f}};
const wb_change k_chg_074[] = {
    {3,3,0x1.6000000000000p+0f}, {3,5,0x1.6000000000000p+0f}, {3,7,0x1.6000000000000p+0f},
    {4,4,0x1.4000000000000p+0f}, {4,6,0x1.4000000000000p+0f}, {4,7,0x1.4000000000000p+0f},
    {5,3,0x1.6000000000000p+0f}, {5,5,0x1.2000000000000p+0f}, {5,7,0x1.2000000000000p+0f},
    {5,9,0x1.4000000000000p-1f}, {6,4,0x1.4000000000000p+0f}, {6,6,0x1.0000000000000p+0f},
    {6,8,0x1.8000000000000p-1f}, {7,3,0x1.6000000000000p+0f}, {7,5,0x1.2000000000000p+0f},
    {7,7,0x1.c000000000000p-1f}};
const wb_change k_chg_075[] = {
    {3,3,0x1.8000000000000p+0f}, {3,5,0x1.8000000000000p+0f}, {3,7,0x1.8000000000000p+0f},
    {4,4,0x1.8000000000000p+0f}, {4,6,0x1.8000000000000p+0f}, {4,7,0x1.8000000000000p+0f},
    {5,3,0x1.8000000000000p+0f}, {5,5,0x1.8000000000000p+0f}, {5,7,0x1.8000000000000p+0f},
    {6,4,0x1.8000000000000p+0f}, {6,6,0x1.8000000000000p+0f}, {6,8,0x1.8000000000000p+0f},
    {7,3,0x1.8000000000000p+0f}, {7,5,0x1.8000000000000p+0f}, {7,7,0x1.8000000000000p+0f}};
const wb_change k_chg_076[] = {
    {3,3,0x1.6000000000000p+0f}, {3,4,0x1.4000000000000p+0f}, {3,5,0x1.2000000000000p+0f},
    {3,7,0x1.6000000000000p+0f}, {4,6,0x1.0000000000000p+0f}, {4,7,0x1.4000000000000p+0f},
    {5,3,0x1.6000000000000p+0f}, {5,5,0x1.6000000000000p+0f}, {5,7,0x1.c000000000000p-1f},
    {5,9,0x1.4000000000000p-1f}, {6,4,0x1.4000000000000p+0f}, {6,5,0x1.4000000000000p+0f},
    {6,8,0x1.8000000000000p-1f}, {7,3,0x1.6000000000000p+0f}, {7,4,0x1.4000000000000p+0f},
    {7,5,0x1.2000000000000p+0f}, {7,6,0x1.0000000000000p+0f}, {7,7,0x1.c000000000000p-1f}};
const wb_change k_chg_077[] = {
    {3,3,0x1.8000000000000p+0f}, {3,4,0x1.8000000000000p+0f}, {3,5,0x1.8000000000000p+0f},
    {3,7,0x1.8000000000000p+0f}, {4,6,0x1.8000000000000p+0f}, {4,7,0x1.8000000000000p+0f},
    {5,3,0x1.8000000000000p+0f}, {5,5,0x1.8000000000000p+0f}, {5,7,0x1.8000000000000p+0f},
    {6,4,0x1.8000000000000p+0f}, {6,5,0x1.8000000000000p+0f}, {6,8,0x1.8000000000000p+0f},
    {7,3,0x1.8000000000000p+0f}, {7,4,0x1.8000000000000p+0f}, {7,5,0x1.8000000000000p+0f},
    {7,6,0x1.8000000000000p+0f}, {7,7,0x1.8000000000000p+0f}};
const wb_change k_chg_078[] = {
    {1,1,0x1.3800000000000p+2f}};
const wb_change k_chg_079[] = {
    {1,1,0x1.4000000000000p+2f}};
const wb_change k_chg_080[] = {
    {1,7,0x1.7000000000000p+3f}, {2,8,0x1.7400000000000p+3f}, {3,9,0x1.7800000000000p+3f},
    {4,10,0x1.7c00000000000p+3f}};
const wb_change k_chg_081[] = {
    {4,9,0x1.8000000000000p+3f}, {4,10,0x1.8000000000000p+3f}, {5,8,0x1.8000000000000p+3f},
    {6,7,0x1.8000000000000p+3f}};
const wb_change k_chg_082[] = {
    {4,7,0x1.7000000000000p+3f}, {4,8,0x1.7400000000000p+3f}, {4,9,0x1.7800000000000p+3f},
    {4,10,0x1.7c00000000000p+3f}};
const wb_change k_chg_083[] = {
    {4,7,0x1.8000000000000p+3f}, {4,8,0x1.8000000000000p+3f}, {4,9,0x1.8000000000000p+3f},
    {4,10,0x1.8000000000000p+3f}};

const wb_run k_wb_runs[] = {
    {"F0", "cs1_md0_inc0.125_d100", 1.0, false, true, 0x1.0000000000000p-3, 100, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_000, 2},
    {"F0", "cs1_md0_incauto_d100", 1.0, false, false, 0x1.c25c268497682p-43, 100, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_001, 2},
    {"F0", "cs1_md1_inc0.125_d100", 1.0, true, true, 0x1.0000000000000p-3, 100, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_002, 2},
    {"F0", "cs1_md1_incauto_d100", 1.0, true, false, 0x1.c25c268497682p-43, 100, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_003, 2},
    {"F0", "cs30_md0_inc0.125_d100", 30.0, false, true, 0x1.0000000000000p-3, 100, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_000, 2},
    {"F0", "cs30_md0_incauto_d100", 30.0, false, false, 0x1.2e95e9e115b9fp-38, 100, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_001, 2},
    {"F0", "cs30_md1_inc0.125_d100", 30.0, true, true, 0x1.0000000000000p-3, 100, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_002, 2},
    {"F0", "cs30_md1_incauto_d100", 30.0, true, false, 0x1.2e95e9e115b9fp-38, 100, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_003, 2},
    {"F1", "cs1_md0_inc0.125_d100", 1.0, false, true, 0x1.0000000000000p-3, 100, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_004, 4},
    {"F1", "cs1_md0_incauto_d100", 1.0, false, false, 0x1.6849b86a12b9bp-46, 100, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_005, 4},
    {"F1", "cs1_md1_inc0.125_d100", 1.0, true, true, 0x1.0000000000000p-3, 100, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_006, 4},
    {"F1", "cs1_md1_incauto_d100", 1.0, true, false, 0x1.6849b86a12b9bp-46, 100, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_007, 4},
    {"F1", "cs30_md0_inc0.125_d100", 30.0, false, true, 0x1.0000000000000p-3, 100, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_004, 4},
    {"F1", "cs30_md0_incauto_d100", 30.0, false, false, 0x1.e4230fce89298p-42, 100, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_005, 4},
    {"F1", "cs30_md1_inc0.125_d100", 30.0, true, true, 0x1.0000000000000p-3, 100, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_006, 4},
    {"F1", "cs30_md1_incauto_d100", 30.0, true, false, 0x1.e4230fce89298p-42, 100, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_007, 4},
    {"F2", "cs1_md0_inc0.125_d100", 1.0, false, true, 0x1.0000000000000p-3, 100, k_inf, {5, 5, 0, 5, 0, 0}, false, true, k_chg_008, 18},
    {"F2", "cs1_md0_incauto_d100", 1.0, false, false, 0x1.6849b86a12b9bp-46, 100, k_inf, {5, 2, 3, 2, 0, 0}, false, true, k_chg_009, 16},
    {"F2", "cs1_md1_inc0.125_d100", 1.0, true, true, 0x1.0000000000000p-3, 100, k_inf, {5, 5, 0, 5, 0, 0}, false, true, k_chg_010, 18},
    {"F2", "cs1_md1_incauto_d100", 1.0, true, false, 0x1.6849b86a12b9bp-46, 100, k_inf, {5, 4, 1, 4, 0, 0}, false, true, k_chg_011, 17},
    {"F2", "cs30_md0_inc0.125_d100", 30.0, false, true, 0x1.0000000000000p-3, 100, k_inf, {5, 5, 0, 5, 0, 0}, false, true, k_chg_008, 18},
    {"F2", "cs30_md0_incauto_d100", 30.0, false, false, 0x1.e4230fce89298p-42, 100, k_inf, {5, 2, 3, 2, 0, 0}, false, true, k_chg_009, 16},
    {"F2", "cs30_md1_inc0.125_d100", 30.0, true, true, 0x1.0000000000000p-3, 100, k_inf, {5, 5, 0, 5, 0, 0}, false, true, k_chg_010, 18},
    {"F2", "cs30_md1_incauto_d100", 30.0, true, false, 0x1.e4230fce89298p-42, 100, k_inf, {5, 4, 1, 4, 0, 0}, false, true, k_chg_011, 17},
    {"F3", "cs1_md0_inc0.125_d100", 1.0, false, true, 0x1.0000000000000p-3, 100, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_004, 4},
    {"F3", "cs1_md0_incauto_d100", 1.0, false, false, 0x1.6849b86a12b9bp-46, 100, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_005, 4},
    {"F3", "cs1_md1_inc0.125_d100", 1.0, true, true, 0x1.0000000000000p-3, 100, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_006, 4},
    {"F3", "cs1_md1_incauto_d100", 1.0, true, false, 0x1.6849b86a12b9bp-46, 100, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_007, 4},
    {"F3", "cs30_md0_inc0.125_d100", 30.0, false, true, 0x1.0000000000000p-3, 100, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_004, 4},
    {"F3", "cs30_md0_incauto_d100", 30.0, false, false, 0x1.e4230fce89298p-42, 100, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_005, 4},
    {"F3", "cs30_md1_inc0.125_d100", 30.0, true, true, 0x1.0000000000000p-3, 100, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_006, 4},
    {"F3", "cs30_md1_incauto_d100", 30.0, true, false, 0x1.e4230fce89298p-42, 100, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_007, 4},
    {"F4", "cs1_md0_inc0.125_d100", 1.0, false, true, 0x1.0000000000000p-3, 100, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_012, 2},
    {"F4", "cs1_md0_incauto_d100", 1.0, false, false, 0x1.c25c268497682p-43, 100, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_013, 2},
    {"F4", "cs1_md1_inc0.125_d100", 1.0, true, true, 0x1.0000000000000p-3, 100, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_014, 2},
    {"F4", "cs1_md1_incauto_d100", 1.0, true, false, 0x1.c25c268497682p-43, 100, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_015, 2},
    {"F4", "cs30_md0_inc0.125_d100", 30.0, false, true, 0x1.0000000000000p-3, 100, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_012, 2},
    {"F4", "cs30_md0_incauto_d100", 30.0, false, false, 0x1.2e95e9e115b9fp-38, 100, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_013, 2},
    {"F4", "cs30_md1_inc0.125_d100", 30.0, true, true, 0x1.0000000000000p-3, 100, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_014, 2},
    {"F4", "cs30_md1_incauto_d100", 30.0, true, false, 0x1.2e95e9e115b9fp-38, 100, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_015, 2},
    {"F5_s10", "cs1_md0_inc0.125_d10", 1.0, false, true, 0x1.0000000000000p-3, 10, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_016, 10},
    {"F5_s10", "cs1_md0_inc0.125_d11", 1.0, false, true, 0x1.0000000000000p-3, 11, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_016, 10},
    {"F5_s10", "cs1_md0_inc0.125_d12", 1.0, false, true, 0x1.0000000000000p-3, 12, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_016, 10},
    {"F5_s10", "cs1_md0_inc0.125_d100", 1.0, false, true, 0x1.0000000000000p-3, 100, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_016, 10},
    {"F5_s10", "cs1_md0_incauto_d10", 1.0, false, false, 0x1.c25c268497682p-43, 10, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_017, 10},
    {"F5_s10", "cs1_md0_incauto_d11", 1.0, false, false, 0x1.c25c268497682p-43, 11, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_017, 10},
    {"F5_s10", "cs1_md0_incauto_d12", 1.0, false, false, 0x1.c25c268497682p-43, 12, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_017, 10},
    {"F5_s10", "cs1_md0_incauto_d100", 1.0, false, false, 0x1.c25c268497682p-43, 100, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_017, 10},
    {"F5_s10", "cs1_md1_inc0.125_d10", 1.0, true, true, 0x1.0000000000000p-3, 10, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_016, 10},
    {"F5_s10", "cs1_md1_inc0.125_d11", 1.0, true, true, 0x1.0000000000000p-3, 11, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_016, 10},
    {"F5_s10", "cs1_md1_inc0.125_d12", 1.0, true, true, 0x1.0000000000000p-3, 12, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_016, 10},
    {"F5_s10", "cs1_md1_inc0.125_d100", 1.0, true, true, 0x1.0000000000000p-3, 100, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_016, 10},
    {"F5_s10", "cs1_md1_incauto_d10", 1.0, true, false, 0x1.c25c268497682p-43, 10, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_017, 10},
    {"F5_s10", "cs1_md1_incauto_d11", 1.0, true, false, 0x1.c25c268497682p-43, 11, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_017, 10},
    {"F5_s10", "cs1_md1_incauto_d12", 1.0, true, false, 0x1.c25c268497682p-43, 12, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_017, 10},
    {"F5_s10", "cs1_md1_incauto_d100", 1.0, true, false, 0x1.c25c268497682p-43, 100, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_017, 10},
    {"F5_s10", "cs30_md0_inc0.125_d10", 30.0, false, true, 0x1.0000000000000p-3, 10, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_016, 10},
    {"F5_s10", "cs30_md0_inc0.125_d11", 30.0, false, true, 0x1.0000000000000p-3, 11, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_016, 10},
    {"F5_s10", "cs30_md0_inc0.125_d12", 30.0, false, true, 0x1.0000000000000p-3, 12, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_016, 10},
    {"F5_s10", "cs30_md0_inc0.125_d100", 30.0, false, true, 0x1.0000000000000p-3, 100, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_016, 10},
    {"F5_s10", "cs30_md0_incauto_d10", 30.0, false, false, 0x1.2e95e9e115b9fp-38, 10, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_017, 10},
    {"F5_s10", "cs30_md0_incauto_d11", 30.0, false, false, 0x1.2e95e9e115b9fp-38, 11, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_017, 10},
    {"F5_s10", "cs30_md0_incauto_d12", 30.0, false, false, 0x1.2e95e9e115b9fp-38, 12, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_017, 10},
    {"F5_s10", "cs30_md0_incauto_d100", 30.0, false, false, 0x1.2e95e9e115b9fp-38, 100, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_017, 10},
    {"F5_s10", "cs30_md1_inc0.125_d10", 30.0, true, true, 0x1.0000000000000p-3, 10, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_016, 10},
    {"F5_s10", "cs30_md1_inc0.125_d11", 30.0, true, true, 0x1.0000000000000p-3, 11, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_016, 10},
    {"F5_s10", "cs30_md1_inc0.125_d12", 30.0, true, true, 0x1.0000000000000p-3, 12, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_016, 10},
    {"F5_s10", "cs30_md1_inc0.125_d100", 30.0, true, true, 0x1.0000000000000p-3, 100, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_016, 10},
    {"F5_s10", "cs30_md1_incauto_d10", 30.0, true, false, 0x1.2e95e9e115b9fp-38, 10, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_017, 10},
    {"F5_s10", "cs30_md1_incauto_d11", 30.0, true, false, 0x1.2e95e9e115b9fp-38, 11, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_017, 10},
    {"F5_s10", "cs30_md1_incauto_d12", 30.0, true, false, 0x1.2e95e9e115b9fp-38, 12, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_017, 10},
    {"F5_s10", "cs30_md1_incauto_d100", 30.0, true, false, 0x1.2e95e9e115b9fp-38, 100, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_017, 10},
    {"F5_d10", "cs1_md0_inc0.125_d10", 1.0, false, true, 0x1.0000000000000p-3, 10, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_018, 10},
    {"F5_d10", "cs1_md0_inc0.125_d11", 1.0, false, true, 0x1.0000000000000p-3, 11, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_018, 10},
    {"F5_d10", "cs1_md0_inc0.125_d12", 1.0, false, true, 0x1.0000000000000p-3, 12, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_018, 10},
    {"F5_d10", "cs1_md0_inc0.125_d100", 1.0, false, true, 0x1.0000000000000p-3, 100, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_018, 10},
    {"F5_d10", "cs1_md0_incauto_d10", 1.0, false, false, 0x1.c25c268497682p-43, 10, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_019, 10},
    {"F5_d10", "cs1_md0_incauto_d11", 1.0, false, false, 0x1.c25c268497682p-43, 11, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_019, 10},
    {"F5_d10", "cs1_md0_incauto_d12", 1.0, false, false, 0x1.c25c268497682p-43, 12, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_019, 10},
    {"F5_d10", "cs1_md0_incauto_d100", 1.0, false, false, 0x1.c25c268497682p-43, 100, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_019, 10},
    {"F5_d10", "cs1_md1_inc0.125_d10", 1.0, true, true, 0x1.0000000000000p-3, 10, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_018, 10},
    {"F5_d10", "cs1_md1_inc0.125_d11", 1.0, true, true, 0x1.0000000000000p-3, 11, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_018, 10},
    {"F5_d10", "cs1_md1_inc0.125_d12", 1.0, true, true, 0x1.0000000000000p-3, 12, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_018, 10},
    {"F5_d10", "cs1_md1_inc0.125_d100", 1.0, true, true, 0x1.0000000000000p-3, 100, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_018, 10},
    {"F5_d10", "cs1_md1_incauto_d10", 1.0, true, false, 0x1.c25c268497682p-43, 10, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_019, 10},
    {"F5_d10", "cs1_md1_incauto_d11", 1.0, true, false, 0x1.c25c268497682p-43, 11, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_019, 10},
    {"F5_d10", "cs1_md1_incauto_d12", 1.0, true, false, 0x1.c25c268497682p-43, 12, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_019, 10},
    {"F5_d10", "cs1_md1_incauto_d100", 1.0, true, false, 0x1.c25c268497682p-43, 100, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_019, 10},
    {"F5_d10", "cs30_md0_inc0.125_d10", 30.0, false, true, 0x1.0000000000000p-3, 10, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_018, 10},
    {"F5_d10", "cs30_md0_inc0.125_d11", 30.0, false, true, 0x1.0000000000000p-3, 11, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_018, 10},
    {"F5_d10", "cs30_md0_inc0.125_d12", 30.0, false, true, 0x1.0000000000000p-3, 12, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_018, 10},
    {"F5_d10", "cs30_md0_inc0.125_d100", 30.0, false, true, 0x1.0000000000000p-3, 100, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_018, 10},
    {"F5_d10", "cs30_md0_incauto_d10", 30.0, false, false, 0x1.2e95e9e115b9fp-38, 10, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_019, 10},
    {"F5_d10", "cs30_md0_incauto_d11", 30.0, false, false, 0x1.2e95e9e115b9fp-38, 11, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_019, 10},
    {"F5_d10", "cs30_md0_incauto_d12", 30.0, false, false, 0x1.2e95e9e115b9fp-38, 12, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_019, 10},
    {"F5_d10", "cs30_md0_incauto_d100", 30.0, false, false, 0x1.2e95e9e115b9fp-38, 100, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_019, 10},
    {"F5_d10", "cs30_md1_inc0.125_d10", 30.0, true, true, 0x1.0000000000000p-3, 10, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_018, 10},
    {"F5_d10", "cs30_md1_inc0.125_d11", 30.0, true, true, 0x1.0000000000000p-3, 11, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_018, 10},
    {"F5_d10", "cs30_md1_inc0.125_d12", 30.0, true, true, 0x1.0000000000000p-3, 12, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_018, 10},
    {"F5_d10", "cs30_md1_inc0.125_d100", 30.0, true, true, 0x1.0000000000000p-3, 100, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_018, 10},
    {"F5_d10", "cs30_md1_incauto_d10", 30.0, true, false, 0x1.2e95e9e115b9fp-38, 10, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_019, 10},
    {"F5_d10", "cs30_md1_incauto_d11", 30.0, true, false, 0x1.2e95e9e115b9fp-38, 11, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_019, 10},
    {"F5_d10", "cs30_md1_incauto_d12", 30.0, true, false, 0x1.2e95e9e115b9fp-38, 12, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_019, 10},
    {"F5_d10", "cs30_md1_incauto_d100", 30.0, true, false, 0x1.2e95e9e115b9fp-38, 100, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_019, 10},
    {"F5_s11", "cs1_md0_inc0.125_d10", 1.0, false, true, 0x1.0000000000000p-3, 10, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_020, 11},
    {"F5_s11", "cs1_md0_inc0.125_d11", 1.0, false, true, 0x1.0000000000000p-3, 11, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_020, 11},
    {"F5_s11", "cs1_md0_inc0.125_d12", 1.0, false, true, 0x1.0000000000000p-3, 12, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_020, 11},
    {"F5_s11", "cs1_md0_inc0.125_d100", 1.0, false, true, 0x1.0000000000000p-3, 100, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_020, 11},
    {"F5_s11", "cs1_md0_incauto_d10", 1.0, false, false, 0x1.c25c268497682p-43, 10, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_021, 11},
    {"F5_s11", "cs1_md0_incauto_d11", 1.0, false, false, 0x1.c25c268497682p-43, 11, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_021, 11},
    {"F5_s11", "cs1_md0_incauto_d12", 1.0, false, false, 0x1.c25c268497682p-43, 12, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_021, 11},
    {"F5_s11", "cs1_md0_incauto_d100", 1.0, false, false, 0x1.c25c268497682p-43, 100, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_021, 11},
    {"F5_s11", "cs1_md1_inc0.125_d10", 1.0, true, true, 0x1.0000000000000p-3, 10, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_020, 11},
    {"F5_s11", "cs1_md1_inc0.125_d11", 1.0, true, true, 0x1.0000000000000p-3, 11, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_020, 11},
    {"F5_s11", "cs1_md1_inc0.125_d12", 1.0, true, true, 0x1.0000000000000p-3, 12, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_020, 11},
    {"F5_s11", "cs1_md1_inc0.125_d100", 1.0, true, true, 0x1.0000000000000p-3, 100, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_020, 11},
    {"F5_s11", "cs1_md1_incauto_d10", 1.0, true, false, 0x1.c25c268497682p-43, 10, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_021, 11},
    {"F5_s11", "cs1_md1_incauto_d11", 1.0, true, false, 0x1.c25c268497682p-43, 11, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_021, 11},
    {"F5_s11", "cs1_md1_incauto_d12", 1.0, true, false, 0x1.c25c268497682p-43, 12, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_021, 11},
    {"F5_s11", "cs1_md1_incauto_d100", 1.0, true, false, 0x1.c25c268497682p-43, 100, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_021, 11},
    {"F5_s11", "cs30_md0_inc0.125_d10", 30.0, false, true, 0x1.0000000000000p-3, 10, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_020, 11},
    {"F5_s11", "cs30_md0_inc0.125_d11", 30.0, false, true, 0x1.0000000000000p-3, 11, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_020, 11},
    {"F5_s11", "cs30_md0_inc0.125_d12", 30.0, false, true, 0x1.0000000000000p-3, 12, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_020, 11},
    {"F5_s11", "cs30_md0_inc0.125_d100", 30.0, false, true, 0x1.0000000000000p-3, 100, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_020, 11},
    {"F5_s11", "cs30_md0_incauto_d10", 30.0, false, false, 0x1.2e95e9e115b9fp-38, 10, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_021, 11},
    {"F5_s11", "cs30_md0_incauto_d11", 30.0, false, false, 0x1.2e95e9e115b9fp-38, 11, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_021, 11},
    {"F5_s11", "cs30_md0_incauto_d12", 30.0, false, false, 0x1.2e95e9e115b9fp-38, 12, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_021, 11},
    {"F5_s11", "cs30_md0_incauto_d100", 30.0, false, false, 0x1.2e95e9e115b9fp-38, 100, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_021, 11},
    {"F5_s11", "cs30_md1_inc0.125_d10", 30.0, true, true, 0x1.0000000000000p-3, 10, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_020, 11},
    {"F5_s11", "cs30_md1_inc0.125_d11", 30.0, true, true, 0x1.0000000000000p-3, 11, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_020, 11},
    {"F5_s11", "cs30_md1_inc0.125_d12", 30.0, true, true, 0x1.0000000000000p-3, 12, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_020, 11},
    {"F5_s11", "cs30_md1_inc0.125_d100", 30.0, true, true, 0x1.0000000000000p-3, 100, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_020, 11},
    {"F5_s11", "cs30_md1_incauto_d10", 30.0, true, false, 0x1.2e95e9e115b9fp-38, 10, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_021, 11},
    {"F5_s11", "cs30_md1_incauto_d11", 30.0, true, false, 0x1.2e95e9e115b9fp-38, 11, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_021, 11},
    {"F5_s11", "cs30_md1_incauto_d12", 30.0, true, false, 0x1.2e95e9e115b9fp-38, 12, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_021, 11},
    {"F5_s11", "cs30_md1_incauto_d100", 30.0, true, false, 0x1.2e95e9e115b9fp-38, 100, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_021, 11},
    {"F5_d11", "cs1_md0_inc0.125_d10", 1.0, false, true, 0x1.0000000000000p-3, 10, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_022, 11},
    {"F5_d11", "cs1_md0_inc0.125_d11", 1.0, false, true, 0x1.0000000000000p-3, 11, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_022, 11},
    {"F5_d11", "cs1_md0_inc0.125_d12", 1.0, false, true, 0x1.0000000000000p-3, 12, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_022, 11},
    {"F5_d11", "cs1_md0_inc0.125_d100", 1.0, false, true, 0x1.0000000000000p-3, 100, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_022, 11},
    {"F5_d11", "cs1_md0_incauto_d10", 1.0, false, false, 0x1.c25c268497682p-43, 10, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_023, 11},
    {"F5_d11", "cs1_md0_incauto_d11", 1.0, false, false, 0x1.c25c268497682p-43, 11, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_023, 11},
    {"F5_d11", "cs1_md0_incauto_d12", 1.0, false, false, 0x1.c25c268497682p-43, 12, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_023, 11},
    {"F5_d11", "cs1_md0_incauto_d100", 1.0, false, false, 0x1.c25c268497682p-43, 100, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_023, 11},
    {"F5_d11", "cs1_md1_inc0.125_d10", 1.0, true, true, 0x1.0000000000000p-3, 10, k_inf, {1, 1, 0, 1, 0, 0}, true, false, k_chg_022, 11},
    {"F5_d11", "cs1_md1_inc0.125_d11", 1.0, true, true, 0x1.0000000000000p-3, 11, k_inf, {1, 1, 0, 1, 0, 0}, true, false, k_chg_022, 11},
    {"F5_d11", "cs1_md1_inc0.125_d12", 1.0, true, true, 0x1.0000000000000p-3, 12, k_inf, {1, 1, 0, 1, 0, 0}, true, false, k_chg_022, 11},
    {"F5_d11", "cs1_md1_inc0.125_d100", 1.0, true, true, 0x1.0000000000000p-3, 100, k_inf, {1, 1, 0, 1, 0, 0}, true, false, k_chg_022, 11},
    {"F5_d11", "cs1_md1_incauto_d10", 1.0, true, false, 0x1.c25c268497682p-43, 10, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_023, 11},
    {"F5_d11", "cs1_md1_incauto_d11", 1.0, true, false, 0x1.c25c268497682p-43, 11, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_023, 11},
    {"F5_d11", "cs1_md1_incauto_d12", 1.0, true, false, 0x1.c25c268497682p-43, 12, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_023, 11},
    {"F5_d11", "cs1_md1_incauto_d100", 1.0, true, false, 0x1.c25c268497682p-43, 100, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_023, 11},
    {"F5_d11", "cs30_md0_inc0.125_d10", 30.0, false, true, 0x1.0000000000000p-3, 10, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_022, 11},
    {"F5_d11", "cs30_md0_inc0.125_d11", 30.0, false, true, 0x1.0000000000000p-3, 11, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_022, 11},
    {"F5_d11", "cs30_md0_inc0.125_d12", 30.0, false, true, 0x1.0000000000000p-3, 12, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_022, 11},
    {"F5_d11", "cs30_md0_inc0.125_d100", 30.0, false, true, 0x1.0000000000000p-3, 100, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_022, 11},
    {"F5_d11", "cs30_md0_incauto_d10", 30.0, false, false, 0x1.2e95e9e115b9fp-38, 10, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_023, 11},
    {"F5_d11", "cs30_md0_incauto_d11", 30.0, false, false, 0x1.2e95e9e115b9fp-38, 11, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_023, 11},
    {"F5_d11", "cs30_md0_incauto_d12", 30.0, false, false, 0x1.2e95e9e115b9fp-38, 12, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_023, 11},
    {"F5_d11", "cs30_md0_incauto_d100", 30.0, false, false, 0x1.2e95e9e115b9fp-38, 100, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_023, 11},
    {"F5_d11", "cs30_md1_inc0.125_d10", 30.0, true, true, 0x1.0000000000000p-3, 10, k_inf, {1, 1, 0, 1, 0, 0}, true, false, k_chg_022, 11},
    {"F5_d11", "cs30_md1_inc0.125_d11", 30.0, true, true, 0x1.0000000000000p-3, 11, k_inf, {1, 1, 0, 1, 0, 0}, true, false, k_chg_022, 11},
    {"F5_d11", "cs30_md1_inc0.125_d12", 30.0, true, true, 0x1.0000000000000p-3, 12, k_inf, {1, 1, 0, 1, 0, 0}, true, false, k_chg_022, 11},
    {"F5_d11", "cs30_md1_inc0.125_d100", 30.0, true, true, 0x1.0000000000000p-3, 100, k_inf, {1, 1, 0, 1, 0, 0}, true, false, k_chg_022, 11},
    {"F5_d11", "cs30_md1_incauto_d10", 30.0, true, false, 0x1.2e95e9e115b9fp-38, 10, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_023, 11},
    {"F5_d11", "cs30_md1_incauto_d11", 30.0, true, false, 0x1.2e95e9e115b9fp-38, 11, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_023, 11},
    {"F5_d11", "cs30_md1_incauto_d12", 30.0, true, false, 0x1.2e95e9e115b9fp-38, 12, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_023, 11},
    {"F5_d11", "cs30_md1_incauto_d100", 30.0, true, false, 0x1.2e95e9e115b9fp-38, 100, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_023, 11},
    {"F5_s12", "cs1_md0_inc0.125_d10", 1.0, false, true, 0x1.0000000000000p-3, 10, k_inf, {1, 1, 0, 0, 1, 0}, true, false, k_chg_024, 1},
    {"F5_s12", "cs1_md0_inc0.125_d11", 1.0, false, true, 0x1.0000000000000p-3, 11, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_025, 12},
    {"F5_s12", "cs1_md0_inc0.125_d12", 1.0, false, true, 0x1.0000000000000p-3, 12, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_025, 12},
    {"F5_s12", "cs1_md0_inc0.125_d100", 1.0, false, true, 0x1.0000000000000p-3, 100, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_025, 12},
    {"F5_s12", "cs1_md0_incauto_d10", 1.0, false, false, 0x1.c25c268497682p-43, 10, k_inf, {1, 1, 0, 0, 1, 0}, true, false, k_chg_026, 1},
    {"F5_s12", "cs1_md0_incauto_d11", 1.0, false, false, 0x1.c25c268497682p-43, 11, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_027, 12},
    {"F5_s12", "cs1_md0_incauto_d12", 1.0, false, false, 0x1.c25c268497682p-43, 12, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_027, 12},
    {"F5_s12", "cs1_md0_incauto_d100", 1.0, false, false, 0x1.c25c268497682p-43, 100, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_027, 12},
    {"F5_s12", "cs1_md1_inc0.125_d10", 1.0, true, true, 0x1.0000000000000p-3, 10, k_inf, {1, 1, 0, 0, 1, 0}, true, false, k_chg_024, 1},
    {"F5_s12", "cs1_md1_inc0.125_d11", 1.0, true, true, 0x1.0000000000000p-3, 11, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_025, 12},
    {"F5_s12", "cs1_md1_inc0.125_d12", 1.0, true, true, 0x1.0000000000000p-3, 12, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_025, 12},
    {"F5_s12", "cs1_md1_inc0.125_d100", 1.0, true, true, 0x1.0000000000000p-3, 100, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_025, 12},
    {"F5_s12", "cs1_md1_incauto_d10", 1.0, true, false, 0x1.c25c268497682p-43, 10, k_inf, {1, 1, 0, 0, 1, 0}, true, false, k_chg_026, 1},
    {"F5_s12", "cs1_md1_incauto_d11", 1.0, true, false, 0x1.c25c268497682p-43, 11, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_027, 12},
    {"F5_s12", "cs1_md1_incauto_d12", 1.0, true, false, 0x1.c25c268497682p-43, 12, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_027, 12},
    {"F5_s12", "cs1_md1_incauto_d100", 1.0, true, false, 0x1.c25c268497682p-43, 100, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_027, 12},
    {"F5_s12", "cs30_md0_inc0.125_d10", 30.0, false, true, 0x1.0000000000000p-3, 10, k_inf, {1, 1, 0, 0, 1, 0}, true, false, k_chg_024, 1},
    {"F5_s12", "cs30_md0_inc0.125_d11", 30.0, false, true, 0x1.0000000000000p-3, 11, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_025, 12},
    {"F5_s12", "cs30_md0_inc0.125_d12", 30.0, false, true, 0x1.0000000000000p-3, 12, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_025, 12},
    {"F5_s12", "cs30_md0_inc0.125_d100", 30.0, false, true, 0x1.0000000000000p-3, 100, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_025, 12},
    {"F5_s12", "cs30_md0_incauto_d10", 30.0, false, false, 0x1.2e95e9e115b9fp-38, 10, k_inf, {1, 1, 0, 0, 1, 0}, true, false, k_chg_026, 1},
    {"F5_s12", "cs30_md0_incauto_d11", 30.0, false, false, 0x1.2e95e9e115b9fp-38, 11, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_027, 12},
    {"F5_s12", "cs30_md0_incauto_d12", 30.0, false, false, 0x1.2e95e9e115b9fp-38, 12, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_027, 12},
    {"F5_s12", "cs30_md0_incauto_d100", 30.0, false, false, 0x1.2e95e9e115b9fp-38, 100, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_027, 12},
    {"F5_s12", "cs30_md1_inc0.125_d10", 30.0, true, true, 0x1.0000000000000p-3, 10, k_inf, {1, 1, 0, 0, 1, 0}, true, false, k_chg_024, 1},
    {"F5_s12", "cs30_md1_inc0.125_d11", 30.0, true, true, 0x1.0000000000000p-3, 11, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_025, 12},
    {"F5_s12", "cs30_md1_inc0.125_d12", 30.0, true, true, 0x1.0000000000000p-3, 12, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_025, 12},
    {"F5_s12", "cs30_md1_inc0.125_d100", 30.0, true, true, 0x1.0000000000000p-3, 100, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_025, 12},
    {"F5_s12", "cs30_md1_incauto_d10", 30.0, true, false, 0x1.2e95e9e115b9fp-38, 10, k_inf, {1, 1, 0, 0, 1, 0}, true, false, k_chg_026, 1},
    {"F5_s12", "cs30_md1_incauto_d11", 30.0, true, false, 0x1.2e95e9e115b9fp-38, 11, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_027, 12},
    {"F5_s12", "cs30_md1_incauto_d12", 30.0, true, false, 0x1.2e95e9e115b9fp-38, 12, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_027, 12},
    {"F5_s12", "cs30_md1_incauto_d100", 30.0, true, false, 0x1.2e95e9e115b9fp-38, 100, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_027, 12},
    {"F5_d12", "cs1_md0_inc0.125_d10", 1.0, false, true, 0x1.0000000000000p-3, 10, k_inf, {1, 1, 0, 0, 1, 0}, true, false, k_chg_028, 1},
    {"F5_d12", "cs1_md0_inc0.125_d11", 1.0, false, true, 0x1.0000000000000p-3, 11, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_029, 12},
    {"F5_d12", "cs1_md0_inc0.125_d12", 1.0, false, true, 0x1.0000000000000p-3, 12, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_029, 12},
    {"F5_d12", "cs1_md0_inc0.125_d100", 1.0, false, true, 0x1.0000000000000p-3, 100, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_029, 12},
    {"F5_d12", "cs1_md0_incauto_d10", 1.0, false, false, 0x1.c25c268497682p-43, 10, k_inf, {1, 1, 0, 0, 1, 0}, true, false, k_chg_030, 1},
    {"F5_d12", "cs1_md0_incauto_d11", 1.0, false, false, 0x1.c25c268497682p-43, 11, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_031, 12},
    {"F5_d12", "cs1_md0_incauto_d12", 1.0, false, false, 0x1.c25c268497682p-43, 12, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_031, 12},
    {"F5_d12", "cs1_md0_incauto_d100", 1.0, false, false, 0x1.c25c268497682p-43, 100, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_031, 12},
    {"F5_d12", "cs1_md1_inc0.125_d10", 1.0, true, true, 0x1.0000000000000p-3, 10, k_inf, {1, 1, 0, 0, 1, 0}, true, false, k_chg_028, 1},
    {"F5_d12", "cs1_md1_inc0.125_d11", 1.0, true, true, 0x1.0000000000000p-3, 11, k_inf, {1, 1, 0, 1, 0, 0}, true, false, k_chg_029, 12},
    {"F5_d12", "cs1_md1_inc0.125_d12", 1.0, true, true, 0x1.0000000000000p-3, 12, k_inf, {1, 1, 0, 1, 0, 0}, true, false, k_chg_029, 12},
    {"F5_d12", "cs1_md1_inc0.125_d100", 1.0, true, true, 0x1.0000000000000p-3, 100, k_inf, {1, 1, 0, 1, 0, 0}, true, false, k_chg_029, 12},
    {"F5_d12", "cs1_md1_incauto_d10", 1.0, true, false, 0x1.c25c268497682p-43, 10, k_inf, {1, 1, 0, 0, 1, 0}, true, false, k_chg_030, 1},
    {"F5_d12", "cs1_md1_incauto_d11", 1.0, true, false, 0x1.c25c268497682p-43, 11, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_031, 12},
    {"F5_d12", "cs1_md1_incauto_d12", 1.0, true, false, 0x1.c25c268497682p-43, 12, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_031, 12},
    {"F5_d12", "cs1_md1_incauto_d100", 1.0, true, false, 0x1.c25c268497682p-43, 100, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_031, 12},
    {"F5_d12", "cs30_md0_inc0.125_d10", 30.0, false, true, 0x1.0000000000000p-3, 10, k_inf, {1, 1, 0, 0, 1, 0}, true, false, k_chg_028, 1},
    {"F5_d12", "cs30_md0_inc0.125_d11", 30.0, false, true, 0x1.0000000000000p-3, 11, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_029, 12},
    {"F5_d12", "cs30_md0_inc0.125_d12", 30.0, false, true, 0x1.0000000000000p-3, 12, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_029, 12},
    {"F5_d12", "cs30_md0_inc0.125_d100", 30.0, false, true, 0x1.0000000000000p-3, 100, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_029, 12},
    {"F5_d12", "cs30_md0_incauto_d10", 30.0, false, false, 0x1.2e95e9e115b9fp-38, 10, k_inf, {1, 1, 0, 0, 1, 0}, true, false, k_chg_030, 1},
    {"F5_d12", "cs30_md0_incauto_d11", 30.0, false, false, 0x1.2e95e9e115b9fp-38, 11, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_031, 12},
    {"F5_d12", "cs30_md0_incauto_d12", 30.0, false, false, 0x1.2e95e9e115b9fp-38, 12, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_031, 12},
    {"F5_d12", "cs30_md0_incauto_d100", 30.0, false, false, 0x1.2e95e9e115b9fp-38, 100, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_031, 12},
    {"F5_d12", "cs30_md1_inc0.125_d10", 30.0, true, true, 0x1.0000000000000p-3, 10, k_inf, {1, 1, 0, 0, 1, 0}, true, false, k_chg_028, 1},
    {"F5_d12", "cs30_md1_inc0.125_d11", 30.0, true, true, 0x1.0000000000000p-3, 11, k_inf, {1, 1, 0, 1, 0, 0}, true, false, k_chg_029, 12},
    {"F5_d12", "cs30_md1_inc0.125_d12", 30.0, true, true, 0x1.0000000000000p-3, 12, k_inf, {1, 1, 0, 1, 0, 0}, true, false, k_chg_029, 12},
    {"F5_d12", "cs30_md1_inc0.125_d100", 30.0, true, true, 0x1.0000000000000p-3, 100, k_inf, {1, 1, 0, 1, 0, 0}, true, false, k_chg_029, 12},
    {"F5_d12", "cs30_md1_incauto_d10", 30.0, true, false, 0x1.2e95e9e115b9fp-38, 10, k_inf, {1, 1, 0, 0, 1, 0}, true, false, k_chg_030, 1},
    {"F5_d12", "cs30_md1_incauto_d11", 30.0, true, false, 0x1.2e95e9e115b9fp-38, 11, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_031, 12},
    {"F5_d12", "cs30_md1_incauto_d12", 30.0, true, false, 0x1.2e95e9e115b9fp-38, 12, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_031, 12},
    {"F5_d12", "cs30_md1_incauto_d100", 30.0, true, false, 0x1.2e95e9e115b9fp-38, 100, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_031, 12},
    {"F6", "cs1_md0_inc0.125_d100", 1.0, false, true, 0x1.0000000000000p-3, 100, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_032, 3},
    {"F6", "cs1_md0_incauto_d100", 1.0, false, false, 0x1.c25c268497682p-43, 100, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_033, 3},
    {"F6", "cs1_md1_inc0.125_d100", 1.0, true, true, 0x1.0000000000000p-3, 100, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_032, 3},
    {"F6", "cs1_md1_incauto_d100", 1.0, true, false, 0x1.c25c268497682p-43, 100, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_033, 3},
    {"F6", "cs30_md0_inc0.125_d100", 30.0, false, true, 0x1.0000000000000p-3, 100, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_032, 3},
    {"F6", "cs30_md0_incauto_d100", 30.0, false, false, 0x1.2e95e9e115b9fp-38, 100, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_033, 3},
    {"F6", "cs30_md1_inc0.125_d100", 30.0, true, true, 0x1.0000000000000p-3, 100, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_032, 3},
    {"F6", "cs30_md1_incauto_d100", 30.0, true, false, 0x1.2e95e9e115b9fp-38, 100, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_033, 3},
    {"F6", "cs1_md0_inc0.125_d100_c1.0", 1.0, false, true, 0x1.0000000000000p-3, 100, 1.0, {1, 1, 0, 0, 0, 1}, false, true, k_chg_034, 1},
    {"F6", "cs1_md0_inc0.125_d100_c1.125", 1.0, false, true, 0x1.0000000000000p-3, 100, 1.125, {1, 1, 0, 0, 0, 1}, false, true, k_chg_034, 1},
    {"F6", "cs1_md0_inc0.125_d100_c1.25", 1.0, false, true, 0x1.0000000000000p-3, 100, 1.25, {1, 1, 0, 0, 0, 1}, false, true, k_chg_034, 1},
    {"F6", "cs1_md0_inc0.125_d100_c1.9", 1.0, false, true, 0x1.0000000000000p-3, 100, 1.9, {1, 1, 0, 0, 0, 1}, false, true, k_chg_034, 1},
    {"F6", "cs1_md0_inc0.125_d100_c2.0", 1.0, false, true, 0x1.0000000000000p-3, 100, 2.0, {1, 1, 0, 1, 0, 0}, false, true, k_chg_032, 3},
    {"F6", "cs1_md0_inc0.125_d100_c2.1", 1.0, false, true, 0x1.0000000000000p-3, 100, 2.1, {1, 1, 0, 1, 0, 0}, false, true, k_chg_032, 3},
    {"F6", "cs1_md1_inc0.125_d100_c1.0", 1.0, true, true, 0x1.0000000000000p-3, 100, 1.0, {1, 1, 0, 0, 0, 1}, false, true, k_chg_034, 1},
    {"F6", "cs1_md1_inc0.125_d100_c1.125", 1.0, true, true, 0x1.0000000000000p-3, 100, 1.125, {1, 1, 0, 1, 0, 0}, false, true, k_chg_032, 3},
    {"F6", "cs1_md1_inc0.125_d100_c1.25", 1.0, true, true, 0x1.0000000000000p-3, 100, 1.25, {1, 1, 0, 1, 0, 0}, false, true, k_chg_032, 3},
    {"F6", "cs1_md1_inc0.125_d100_c1.9", 1.0, true, true, 0x1.0000000000000p-3, 100, 1.9, {1, 1, 0, 1, 0, 0}, false, true, k_chg_032, 3},
    {"F6", "cs1_md1_inc0.125_d100_c2.0", 1.0, true, true, 0x1.0000000000000p-3, 100, 2.0, {1, 1, 0, 1, 0, 0}, false, true, k_chg_032, 3},
    {"F6", "cs1_md1_inc0.125_d100_c2.1", 1.0, true, true, 0x1.0000000000000p-3, 100, 2.1, {1, 1, 0, 1, 0, 0}, false, true, k_chg_032, 3},
    {"F7", "cs1_md0_inc0.125_d100", 1.0, false, true, 0x1.0000000000000p-3, 100, k_inf, {9, 6, 3, 6, 0, 0}, false, true, k_chg_035, 16},
    {"F7", "cs1_md0_incauto_d100", 1.0, false, false, 0x1.6849b86a12b9bp-46, 100, k_inf, {9, 6, 3, 6, 0, 0}, false, true, k_chg_036, 15},
    {"F7", "cs1_md1_inc0.125_d100", 1.0, true, true, 0x1.0000000000000p-3, 100, k_inf, {9, 8, 1, 8, 0, 0}, false, true, k_chg_037, 18},
    {"F7", "cs1_md1_incauto_d100", 1.0, true, false, 0x1.6849b86a12b9bp-46, 100, k_inf, {9, 8, 1, 8, 0, 0}, false, true, k_chg_038, 17},
    {"F7", "cs30_md0_inc0.125_d100", 30.0, false, true, 0x1.0000000000000p-3, 100, k_inf, {9, 6, 3, 6, 0, 0}, false, true, k_chg_035, 16},
    {"F7", "cs30_md0_incauto_d100", 30.0, false, false, 0x1.e4230fce89298p-42, 100, k_inf, {9, 6, 3, 6, 0, 0}, false, true, k_chg_036, 15},
    {"F7", "cs30_md1_inc0.125_d100", 30.0, true, true, 0x1.0000000000000p-3, 100, k_inf, {9, 8, 1, 8, 0, 0}, false, true, k_chg_037, 18},
    {"F7", "cs30_md1_incauto_d100", 30.0, true, false, 0x1.e4230fce89298p-42, 100, k_inf, {9, 8, 1, 8, 0, 0}, false, true, k_chg_038, 17},
    {"F8", "cs1_md0_inc0.125_d100", 1.0, false, true, 0x1.0000000000000p-3, 100, k_inf, {20, 18, 2, 2, 16, 0}, true, false, k_chg_039, 25},
    {"F8", "cs1_md0_incauto_d100", 1.0, false, false, 0x1.6849b86a12b9bp-46, 100, k_inf, {20, 10, 10, 10, 0, 0}, false, true, k_chg_040, 31},
    {"F8", "cs1_md1_inc0.125_d100", 1.0, true, true, 0x1.0000000000000p-3, 100, k_inf, {20, 20, 0, 4, 16, 0}, true, false, k_chg_041, 27},
    {"F8", "cs1_md1_incauto_d100", 1.0, true, false, 0x1.6849b86a12b9bp-46, 100, k_inf, {20, 20, 0, 20, 0, 0}, false, true, k_chg_042, 41},
    {"F8", "cs30_md0_inc0.125_d100", 30.0, false, true, 0x1.0000000000000p-3, 100, k_inf, {20, 18, 2, 2, 16, 0}, true, false, k_chg_039, 25},
    {"F8", "cs30_md0_incauto_d100", 30.0, false, false, 0x1.e4230fce89298p-42, 100, k_inf, {20, 10, 10, 10, 0, 0}, false, true, k_chg_040, 31},
    {"F8", "cs30_md1_inc0.125_d100", 30.0, true, true, 0x1.0000000000000p-3, 100, k_inf, {20, 20, 0, 4, 16, 0}, true, false, k_chg_041, 27},
    {"F8", "cs30_md1_incauto_d100", 30.0, true, false, 0x1.e4230fce89298p-42, 100, k_inf, {20, 20, 0, 20, 0, 0}, false, true, k_chg_042, 41},
    {"F9a", "cs1_md0_inc0.125_d100", 1.0, false, true, 0x1.0000000000000p-3, 100, k_inf, {0, 0, 0, 0, 0, 0}, false, true, k_chg_043, 0},
    {"F9a", "cs1_md0_incauto_d100", 1.0, false, false, 0x1.6849b86a12b9bp-46, 100, k_inf, {0, 0, 0, 0, 0, 0}, false, true, k_chg_043, 0},
    {"F9a", "cs1_md1_inc0.125_d100", 1.0, true, true, 0x1.0000000000000p-3, 100, k_inf, {0, 0, 0, 0, 0, 0}, false, true, k_chg_043, 0},
    {"F9a", "cs1_md1_incauto_d100", 1.0, true, false, 0x1.6849b86a12b9bp-46, 100, k_inf, {0, 0, 0, 0, 0, 0}, false, true, k_chg_043, 0},
    {"F9a", "cs30_md0_inc0.125_d100", 30.0, false, true, 0x1.0000000000000p-3, 100, k_inf, {0, 0, 0, 0, 0, 0}, false, true, k_chg_043, 0},
    {"F9a", "cs30_md0_incauto_d100", 30.0, false, false, 0x1.e4230fce89298p-42, 100, k_inf, {0, 0, 0, 0, 0, 0}, false, true, k_chg_043, 0},
    {"F9a", "cs30_md1_inc0.125_d100", 30.0, true, true, 0x1.0000000000000p-3, 100, k_inf, {0, 0, 0, 0, 0, 0}, false, true, k_chg_043, 0},
    {"F9a", "cs30_md1_incauto_d100", 30.0, true, false, 0x1.e4230fce89298p-42, 100, k_inf, {0, 0, 0, 0, 0, 0}, false, true, k_chg_043, 0},
    {"F9b", "cs1_md0_inc0.125_d100", 1.0, false, true, 0x1.0000000000000p-3, 100, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_044, 3},
    {"F9b", "cs1_md0_incauto_d100", 1.0, false, false, 0x1.6849b86a12b9bp-46, 100, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_045, 3},
    {"F9b", "cs1_md1_inc0.125_d100", 1.0, true, true, 0x1.0000000000000p-3, 100, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_046, 3},
    {"F9b", "cs1_md1_incauto_d100", 1.0, true, false, 0x1.6849b86a12b9bp-46, 100, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_047, 3},
    {"F9b", "cs30_md0_inc0.125_d100", 30.0, false, true, 0x1.0000000000000p-3, 100, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_044, 3},
    {"F9b", "cs30_md0_incauto_d100", 30.0, false, false, 0x1.e4230fce89298p-42, 100, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_045, 3},
    {"F9b", "cs30_md1_inc0.125_d100", 30.0, true, true, 0x1.0000000000000p-3, 100, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_046, 3},
    {"F9b", "cs30_md1_incauto_d100", 30.0, true, false, 0x1.e4230fce89298p-42, 100, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_047, 3},
    {"F9c", "cs1_md0_inc0.125_d100", 1.0, false, true, 0x1.0000000000000p-3, 100, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_048, 4},
    {"F9c", "cs1_md0_incauto_d100", 1.0, false, false, 0x1.6849b86a12b9bp-46, 100, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_049, 4},
    {"F9c", "cs1_md1_inc0.125_d100", 1.0, true, true, 0x1.0000000000000p-3, 100, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_050, 4},
    {"F9c", "cs1_md1_incauto_d100", 1.0, true, false, 0x1.6849b86a12b9bp-46, 100, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_051, 4},
    {"F9c", "cs30_md0_inc0.125_d100", 30.0, false, true, 0x1.0000000000000p-3, 100, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_048, 4},
    {"F9c", "cs30_md0_incauto_d100", 30.0, false, false, 0x1.e4230fce89298p-42, 100, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_049, 4},
    {"F9c", "cs30_md1_inc0.125_d100", 30.0, true, true, 0x1.0000000000000p-3, 100, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_050, 4},
    {"F9c", "cs30_md1_incauto_d100", 30.0, true, false, 0x1.e4230fce89298p-42, 100, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_051, 4},
    {"F10", "cs1_md0_inc0.125_d100", 1.0, false, true, 0x1.0000000000000p-3, 100, k_inf, {1, 1, 0, 0, 1, 0}, true, false, k_chg_052, 1},
    {"F10", "cs1_md0_incauto_d100", 1.0, false, false, 0x1.c25c268497682p-43, 100, k_inf, {1, 1, 0, 0, 1, 0}, true, true, k_chg_053, 1},
    {"F10", "cs1_md1_inc0.125_d100", 1.0, true, true, 0x1.0000000000000p-3, 100, k_inf, {1, 1, 0, 0, 1, 0}, true, false, k_chg_052, 1},
    {"F10", "cs1_md1_incauto_d100", 1.0, true, false, 0x1.c25c268497682p-43, 100, k_inf, {1, 1, 0, 0, 1, 0}, true, true, k_chg_053, 1},
    {"F10", "cs30_md0_inc0.125_d100", 30.0, false, true, 0x1.0000000000000p-3, 100, k_inf, {1, 1, 0, 0, 1, 0}, true, false, k_chg_052, 1},
    {"F10", "cs30_md0_incauto_d100", 30.0, false, false, 0x1.2e95e9e115b9fp-38, 100, k_inf, {1, 1, 0, 0, 1, 0}, true, true, k_chg_053, 1},
    {"F10", "cs30_md1_inc0.125_d100", 30.0, true, true, 0x1.0000000000000p-3, 100, k_inf, {1, 1, 0, 0, 1, 0}, true, false, k_chg_052, 1},
    {"F10", "cs30_md1_incauto_d100", 30.0, true, false, 0x1.2e95e9e115b9fp-38, 100, k_inf, {1, 1, 0, 0, 1, 0}, true, true, k_chg_053, 1},
    {"F11", "cs1_md0_inc0.125_d100", 1.0, false, true, 0x1.0000000000000p-3, 100, k_inf, {2, 2, 0, 2, 0, 0}, true, false, k_chg_054, 9},
    {"F11", "cs1_md0_incauto_d100", 1.0, false, false, 0x1.6849b86a12b9bp-46, 100, k_inf, {2, 2, 0, 2, 0, 0}, true, false, k_chg_055, 9},
    {"F11", "cs1_md1_inc0.125_d100", 1.0, true, true, 0x1.0000000000000p-3, 100, k_inf, {2, 2, 0, 2, 0, 0}, true, false, k_chg_056, 9},
    {"F11", "cs1_md1_incauto_d100", 1.0, true, false, 0x1.6849b86a12b9bp-46, 100, k_inf, {2, 2, 0, 2, 0, 0}, true, false, k_chg_057, 9},
    {"F11", "cs30_md0_inc0.125_d100", 30.0, false, true, 0x1.0000000000000p-3, 100, k_inf, {2, 2, 0, 2, 0, 0}, true, false, k_chg_054, 9},
    {"F11", "cs30_md0_incauto_d100", 30.0, false, false, 0x1.e4230fce89298p-42, 100, k_inf, {2, 2, 0, 2, 0, 0}, true, false, k_chg_055, 9},
    {"F11", "cs30_md1_inc0.125_d100", 30.0, true, true, 0x1.0000000000000p-3, 100, k_inf, {2, 2, 0, 2, 0, 0}, true, false, k_chg_056, 9},
    {"F11", "cs30_md1_incauto_d100", 30.0, true, false, 0x1.e4230fce89298p-42, 100, k_inf, {2, 2, 0, 2, 0, 0}, true, false, k_chg_057, 9},
    {"F12_0", "cs1_md0_inc0.125_d100", 1.0, false, true, 0x1.0000000000000p-3, 100, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_058, 4},
    {"F12_0", "cs1_md0_incauto_d100", 1.0, false, false, 0x1.6849b86a12b9bp-46, 100, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_059, 4},
    {"F12_0", "cs1_md1_inc0.125_d100", 1.0, true, true, 0x1.0000000000000p-3, 100, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_060, 4},
    {"F12_0", "cs1_md1_incauto_d100", 1.0, true, false, 0x1.6849b86a12b9bp-46, 100, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_061, 4},
    {"F12_0", "cs30_md0_inc0.125_d100", 30.0, false, true, 0x1.0000000000000p-3, 100, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_058, 4},
    {"F12_0", "cs30_md0_incauto_d100", 30.0, false, false, 0x1.e4230fce89298p-42, 100, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_059, 4},
    {"F12_0", "cs30_md1_inc0.125_d100", 30.0, true, true, 0x1.0000000000000p-3, 100, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_060, 4},
    {"F12_0", "cs30_md1_incauto_d100", 30.0, true, false, 0x1.e4230fce89298p-42, 100, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_061, 4},
    {"F12_m50", "cs1_md0_inc0.125_d100", 1.0, false, true, 0x1.0000000000000p-3, 100, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_062, 4},
    {"F12_m50", "cs1_md0_incauto_d100", 1.0, false, false, 0x1.c25c268497682p-43, 100, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_063, 4},
    {"F12_m50", "cs1_md1_inc0.125_d100", 1.0, true, true, 0x1.0000000000000p-3, 100, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_064, 4},
    {"F12_m50", "cs1_md1_incauto_d100", 1.0, true, false, 0x1.c25c268497682p-43, 100, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_065, 4},
    {"F12_m50", "cs30_md0_inc0.125_d100", 30.0, false, true, 0x1.0000000000000p-3, 100, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_062, 4},
    {"F12_m50", "cs30_md0_incauto_d100", 30.0, false, false, 0x1.2e95e9e115b9fp-38, 100, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_063, 4},
    {"F12_m50", "cs30_md1_inc0.125_d100", 30.0, true, true, 0x1.0000000000000p-3, 100, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_064, 4},
    {"F12_m50", "cs30_md1_incauto_d100", 30.0, true, false, 0x1.2e95e9e115b9fp-38, 100, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_065, 4},
    {"F12_4000", "cs1_md0_inc0.125_d100", 1.0, false, true, 0x1.0000000000000p-3, 100, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_066, 4},
    {"F12_4000", "cs1_md0_incauto_d100", 1.0, false, false, 0x1.5fd7fe1796495p-36, 100, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_067, 4},
    {"F12_4000", "cs1_md1_inc0.125_d100", 1.0, true, true, 0x1.0000000000000p-3, 100, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_068, 4},
    {"F12_4000", "cs1_md1_incauto_d100", 1.0, true, false, 0x1.5fd7fe1796495p-36, 100, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_069, 4},
    {"F12_4000", "cs30_md0_inc0.125_d100", 30.0, false, true, 0x1.0000000000000p-3, 100, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_066, 4},
    {"F12_4000", "cs30_md0_incauto_d100", 30.0, false, false, 0x1.d8ca3d6fb1f28p-32, 100, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_067, 4},
    {"F12_4000", "cs30_md1_inc0.125_d100", 30.0, true, true, 0x1.0000000000000p-3, 100, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_068, 4},
    {"F12_4000", "cs30_md1_incauto_d100", 30.0, true, false, 0x1.d8ca3d6fb1f28p-32, 100, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_069, 4},
    {"F13_F1", "cs1_md0_inc0.125_d100", 1.0, false, true, 0x1.0000000000000p-3, 100, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_070, 4},
    {"F13_F1", "cs1_md0_incauto_d100", 1.0, false, false, 0x1.6849b86a12b9bp-46, 100, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_071, 4},
    {"F13_F1", "cs1_md1_inc0.125_d100", 1.0, true, true, 0x1.0000000000000p-3, 100, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_072, 4},
    {"F13_F1", "cs1_md1_incauto_d100", 1.0, true, false, 0x1.6849b86a12b9bp-46, 100, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_073, 4},
    {"F13_F1", "cs30_md0_inc0.125_d100", 30.0, false, true, 0x1.0000000000000p-3, 100, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_070, 4},
    {"F13_F1", "cs30_md0_incauto_d100", 30.0, false, false, 0x1.e4230fce89298p-42, 100, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_071, 4},
    {"F13_F1", "cs30_md1_inc0.125_d100", 30.0, true, true, 0x1.0000000000000p-3, 100, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_072, 4},
    {"F13_F1", "cs30_md1_incauto_d100", 30.0, true, false, 0x1.e4230fce89298p-42, 100, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_073, 4},
    {"F13_F3", "cs1_md0_inc0.125_d100", 1.0, false, true, 0x1.0000000000000p-3, 100, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_070, 4},
    {"F13_F3", "cs1_md0_incauto_d100", 1.0, false, false, 0x1.6849b86a12b9bp-46, 100, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_071, 4},
    {"F13_F3", "cs1_md1_inc0.125_d100", 1.0, true, true, 0x1.0000000000000p-3, 100, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_072, 4},
    {"F13_F3", "cs1_md1_incauto_d100", 1.0, true, false, 0x1.6849b86a12b9bp-46, 100, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_073, 4},
    {"F13_F3", "cs30_md0_inc0.125_d100", 30.0, false, true, 0x1.0000000000000p-3, 100, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_070, 4},
    {"F13_F3", "cs30_md0_incauto_d100", 30.0, false, false, 0x1.e4230fce89298p-42, 100, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_071, 4},
    {"F13_F3", "cs30_md1_inc0.125_d100", 30.0, true, true, 0x1.0000000000000p-3, 100, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_072, 4},
    {"F13_F3", "cs30_md1_incauto_d100", 30.0, true, false, 0x1.e4230fce89298p-42, 100, k_inf, {1, 1, 0, 1, 0, 0}, false, true, k_chg_073, 4},
    {"F13_F7", "cs1_md0_inc0.125_d100", 1.0, false, true, 0x1.0000000000000p-3, 100, k_inf, {9, 6, 3, 6, 0, 0}, false, true, k_chg_074, 16},
    {"F13_F7", "cs1_md0_incauto_d100", 1.0, false, false, 0x1.6849b86a12b9bp-46, 100, k_inf, {9, 6, 3, 6, 0, 0}, false, true, k_chg_075, 15},
    {"F13_F7", "cs1_md1_inc0.125_d100", 1.0, true, true, 0x1.0000000000000p-3, 100, k_inf, {9, 8, 1, 8, 0, 0}, false, true, k_chg_076, 18},
    {"F13_F7", "cs1_md1_incauto_d100", 1.0, true, false, 0x1.6849b86a12b9bp-46, 100, k_inf, {9, 8, 1, 8, 0, 0}, false, true, k_chg_077, 17},
    {"F13_F7", "cs30_md0_inc0.125_d100", 30.0, false, true, 0x1.0000000000000p-3, 100, k_inf, {9, 6, 3, 6, 0, 0}, false, true, k_chg_074, 16},
    {"F13_F7", "cs30_md0_incauto_d100", 30.0, false, false, 0x1.e4230fce89298p-42, 100, k_inf, {9, 6, 3, 6, 0, 0}, false, true, k_chg_075, 15},
    {"F13_F7", "cs30_md1_inc0.125_d100", 30.0, true, true, 0x1.0000000000000p-3, 100, k_inf, {9, 8, 1, 8, 0, 0}, false, true, k_chg_076, 18},
    {"F13_F7", "cs30_md1_incauto_d100", 30.0, true, false, 0x1.e4230fce89298p-42, 100, k_inf, {9, 8, 1, 8, 0, 0}, false, true, k_chg_077, 17},
    {"F14_3x3", "cs1_md0_inc0.125_d100", 1.0, false, true, 0x1.0000000000000p-3, 100, k_inf, {1, 1, 0, 0, 1, 0}, true, false, k_chg_078, 1},
    {"F14_3x3", "cs1_md0_incauto_d100", 1.0, false, false, 0x1.6849b86a12b9bp-46, 100, k_inf, {1, 1, 0, 0, 1, 0}, true, true, k_chg_079, 1},
    {"F14_3x3", "cs1_md1_inc0.125_d100", 1.0, true, true, 0x1.0000000000000p-3, 100, k_inf, {1, 1, 0, 0, 1, 0}, true, false, k_chg_078, 1},
    {"F14_3x3", "cs1_md1_incauto_d100", 1.0, true, false, 0x1.6849b86a12b9bp-46, 100, k_inf, {1, 1, 0, 0, 1, 0}, true, true, k_chg_079, 1},
    {"F14_3x3", "cs30_md0_inc0.125_d100", 30.0, false, true, 0x1.0000000000000p-3, 100, k_inf, {1, 1, 0, 0, 1, 0}, true, false, k_chg_078, 1},
    {"F14_3x3", "cs30_md0_incauto_d100", 30.0, false, false, 0x1.e4230fce89298p-42, 100, k_inf, {1, 1, 0, 0, 1, 0}, true, true, k_chg_079, 1},
    {"F14_3x3", "cs30_md1_inc0.125_d100", 30.0, true, true, 0x1.0000000000000p-3, 100, k_inf, {1, 1, 0, 0, 1, 0}, true, false, k_chg_078, 1},
    {"F14_3x3", "cs30_md1_incauto_d100", 30.0, true, false, 0x1.e4230fce89298p-42, 100, k_inf, {1, 1, 0, 0, 1, 0}, true, true, k_chg_079, 1},
    {"F14_1x9", "cs1_md0_inc0.125_d100", 1.0, false, true, 0x1.0000000000000p-3, 100, k_inf, {0, 0, 0, 0, 0, 0}, false, true, k_chg_043, 0},
    {"F14_1x9", "cs1_md0_incauto_d100", 1.0, false, false, 0x1.6849b86a12b9bp-46, 100, k_inf, {0, 0, 0, 0, 0, 0}, false, true, k_chg_043, 0},
    {"F14_1x9", "cs1_md1_inc0.125_d100", 1.0, true, true, 0x1.0000000000000p-3, 100, k_inf, {0, 0, 0, 0, 0, 0}, false, true, k_chg_043, 0},
    {"F14_1x9", "cs1_md1_incauto_d100", 1.0, true, false, 0x1.6849b86a12b9bp-46, 100, k_inf, {0, 0, 0, 0, 0, 0}, false, true, k_chg_043, 0},
    {"F14_1x9", "cs30_md0_inc0.125_d100", 30.0, false, true, 0x1.0000000000000p-3, 100, k_inf, {0, 0, 0, 0, 0, 0}, false, true, k_chg_043, 0},
    {"F14_1x9", "cs30_md0_incauto_d100", 30.0, false, false, 0x1.e4230fce89298p-42, 100, k_inf, {0, 0, 0, 0, 0, 0}, false, true, k_chg_043, 0},
    {"F14_1x9", "cs30_md1_inc0.125_d100", 30.0, true, true, 0x1.0000000000000p-3, 100, k_inf, {0, 0, 0, 0, 0, 0}, false, true, k_chg_043, 0},
    {"F14_1x9", "cs30_md1_incauto_d100", 30.0, true, false, 0x1.e4230fce89298p-42, 100, k_inf, {0, 0, 0, 0, 0, 0}, false, true, k_chg_043, 0},
    {"F14_allnd", "cs1_md0_inc0.125_d100", 1.0, false, true, 0x1.0000000000000p-3, 100, k_inf, {0, 0, 0, 0, 0, 0}, false, true, k_chg_043, 0},
    {"F14_allnd", "cs1_md0_incauto_d100", 1.0, false, false, 0x1.5798ee2308c3ap-27, 100, k_inf, {0, 0, 0, 0, 0, 0}, false, true, k_chg_043, 0},
    {"F14_allnd", "cs1_md1_inc0.125_d100", 1.0, true, true, 0x1.0000000000000p-3, 100, k_inf, {0, 0, 0, 0, 0, 0}, false, true, k_chg_043, 0},
    {"F14_allnd", "cs1_md1_incauto_d100", 1.0, true, false, 0x1.5798ee2308c3ap-27, 100, k_inf, {0, 0, 0, 0, 0, 0}, false, true, k_chg_043, 0},
    {"F14_allnd", "cs30_md0_inc0.125_d100", 30.0, false, true, 0x1.0000000000000p-3, 100, k_inf, {0, 0, 0, 0, 0, 0}, false, true, k_chg_043, 0},
    {"F14_allnd", "cs30_md0_incauto_d100", 30.0, false, false, 0x1.01b2b29a4692cp-24, 100, k_inf, {0, 0, 0, 0, 0, 0}, false, true, k_chg_043, 0},
    {"F14_allnd", "cs30_md1_inc0.125_d100", 30.0, true, true, 0x1.0000000000000p-3, 100, k_inf, {0, 0, 0, 0, 0, 0}, false, true, k_chg_043, 0},
    {"F14_allnd", "cs30_md1_incauto_d100", 30.0, true, false, 0x1.01b2b29a4692cp-24, 100, k_inf, {0, 0, 0, 0, 0, 0}, false, true, k_chg_043, 0},
    {"F14_nopit", "cs1_md0_inc0.125_d100", 1.0, false, true, 0x1.0000000000000p-3, 100, k_inf, {0, 0, 0, 0, 0, 0}, false, true, k_chg_043, 0},
    {"F14_nopit", "cs1_md0_incauto_d100", 1.0, false, false, 0x1.c25c268497682p-43, 100, k_inf, {0, 0, 0, 0, 0, 0}, false, true, k_chg_043, 0},
    {"F14_nopit", "cs1_md1_inc0.125_d100", 1.0, true, true, 0x1.0000000000000p-3, 100, k_inf, {0, 0, 0, 0, 0, 0}, false, true, k_chg_043, 0},
    {"F14_nopit", "cs1_md1_incauto_d100", 1.0, true, false, 0x1.c25c268497682p-43, 100, k_inf, {0, 0, 0, 0, 0, 0}, false, true, k_chg_043, 0},
    {"F14_nopit", "cs30_md0_inc0.125_d100", 30.0, false, true, 0x1.0000000000000p-3, 100, k_inf, {0, 0, 0, 0, 0, 0}, false, true, k_chg_043, 0},
    {"F14_nopit", "cs30_md0_incauto_d100", 30.0, false, false, 0x1.2e95e9e115b9fp-38, 100, k_inf, {0, 0, 0, 0, 0, 0}, false, true, k_chg_043, 0},
    {"F14_nopit", "cs30_md1_inc0.125_d100", 30.0, true, true, 0x1.0000000000000p-3, 100, k_inf, {0, 0, 0, 0, 0, 0}, false, true, k_chg_043, 0},
    {"F14_nopit", "cs30_md1_incauto_d100", 30.0, true, false, 0x1.2e95e9e115b9fp-38, 100, k_inf, {0, 0, 0, 0, 0, 0}, false, true, k_chg_043, 0},
    {"F15", "cs1_md0_inc0.125_d100", 1.0, false, true, 0x1.0000000000000p-3, 100, k_inf, {1, 1, 0, 1, 0, 0}, true, false, k_chg_080, 4},
    {"F15", "cs1_md0_incauto_d100", 1.0, false, false, 0x1.c25c268497682p-43, 100, k_inf, {1, 1, 0, 1, 0, 0}, true, false, k_chg_081, 4},
    {"F15", "cs1_md1_inc0.125_d100", 1.0, true, true, 0x1.0000000000000p-3, 100, k_inf, {1, 1, 0, 1, 0, 0}, true, false, k_chg_082, 4},
    {"F15", "cs1_md1_incauto_d100", 1.0, true, false, 0x1.c25c268497682p-43, 100, k_inf, {1, 1, 0, 1, 0, 0}, true, false, k_chg_083, 4},
    {"F15", "cs30_md0_inc0.125_d100", 30.0, false, true, 0x1.0000000000000p-3, 100, k_inf, {1, 1, 0, 1, 0, 0}, true, false, k_chg_080, 4},
    {"F15", "cs30_md0_incauto_d100", 30.0, false, false, 0x1.2e95e9e115b9fp-38, 100, k_inf, {1, 1, 0, 1, 0, 0}, true, false, k_chg_081, 4},
    {"F15", "cs30_md1_inc0.125_d100", 30.0, true, true, 0x1.0000000000000p-3, 100, k_inf, {1, 1, 0, 1, 0, 0}, true, false, k_chg_082, 4},
    {"F15", "cs30_md1_incauto_d100", 30.0, true, false, 0x1.2e95e9e115b9fp-38, 100, k_inf, {1, 1, 0, 1, 0, 0}, true, false, k_chg_083, 4},
};
// ---- END GENERATED LITERALS ----


// =============================================================================
// Fixture builders (in-process). Values are computed in double and rounded
// once to float32, exactly as the characterisation built them.
// =============================================================================
struct grid
{
    int rows = 0;
    int cols = 0;
    std::vector<double> z;

    grid(int r, int c, double v = 0.0) : rows(r), cols(c), z(static_cast<std::size_t>(r) * c, v) {}
    double& at(int r, int c) { return z[static_cast<std::size_t>(r) * cols + c]; }
    double  at(int r, int c) const { return z[static_cast<std::size_t>(r) * cols + c]; }
};

int cheb(int r, int c, int r0, int c0)
{
    return std::max(std::abs(r - r0), std::abs(c - c0));
}

grid f0()
{
    grid g(9, 9);
    for (int r = 0; r < 9; ++r)
        for (int c = 0; c < 9; ++c)
            g.at(r, c) = 10.0 + c;
    g.at(4, 4) -= 5.0;
    return g;
}

grid bowl11(const std::vector<std::pair<int,int>>& saddles)
{
    const double lut[6] = {1.0, 2.0, 4.0, 5.0, 0.5, 0.0};
    grid g(11, 11);
    for (int r = 0; r < 11; ++r)
        for (int c = 0; c < 11; ++c)
            g.at(r, c) = lut[cheb(r, c, 5, 5)];
    for (const auto& s : saddles) g.at(s.first, s.second) = 3.0;
    return g;
}

grid f1() { return bowl11({{5, 8}}); }
grid f3() { return bowl11({{5, 8}, {5, 2}}); }

grid f2()
{
    const int rows = 11, cols = 19;
    grid g(rows, cols, 4.0);
    for (int c = 0; c < cols; ++c) { g.at(0, c) = 0.0; g.at(rows - 1, c) = 0.0; }
    for (int r = 0; r < rows; ++r) { g.at(r, 0) = 0.0; g.at(r, cols - 1) = 0.0; }
    for (int c = 1; c < cols - 1; ++c) { g.at(1, c) = 6.0; g.at(rows - 2, c) = 6.0; }
    for (int r = 1; r < rows - 1; ++r) { g.at(r, 1) = 6.0; g.at(r, cols - 2) = 6.0; }
    g.at(5, 17) = 3.5;
    for (int r = 3; r < 8; ++r)
        for (int c = 8; c < 11; ++c)
            g.at(r, c) = 3.0;
    const int    centres[2][2] = {{5, 5}, {5, 13}};
    const double floors[2]     = {1.0, 2.0};
    for (int b = 0; b < 2; ++b)
    {
        for (int r = 0; r < rows; ++r)
            for (int c = 0; c < cols; ++c)
                if (cheb(r, c, centres[b][0], centres[b][1]) == 2) g.at(r, c) = 2.75;
        for (int r = 0; r < rows; ++r)
            for (int c = 0; c < cols; ++c)
                if (cheb(r, c, centres[b][0], centres[b][1]) == 1) g.at(r, c) = 2.5;
        g.at(centres[b][0], centres[b][1]) = floors[b];
    }
    return g;
}

grid f4()
{
    grid g(9, 9);
    for (int r = 0; r < 9; ++r)
        for (int c = 0; c < 9; ++c)
            g.at(r, c) = 10.0 + cheb(r, c, 4, 4);
    g.at(4, 4) = 1.0;
    g.at(4, 5) = 5.0;
    g.at(3, 5) = 5.0;
    g.at(4, 6) = 0.0; g.at(4, 7) = -1.0; g.at(4, 8) = -2.0;
    g.at(2, 6) = 0.0; g.at(1, 7) = -1.0; g.at(0, 8) = -2.0;
    return g;
}

// Pit + ascending trough of L-1 cells (z = 5 + 0.25 k) + outlet at step L
// (z = 0) + a channel falling to the raster edge; cone 20 + d elsewhere.
grid f5(int L, bool diagonal)
{
    int rows, cols, r0, c0;
    std::vector<std::pair<int,int>> cells;
    if (diagonal)
    {
        const int n = L + 5;
        rows = cols = n;
        r0 = c0 = 2;
        for (int k = 0; k < n; ++k) cells.emplace_back(r0 + k, c0 + k);
    }
    else
    {
        rows = 7; cols = L + 5; r0 = 3; c0 = 2;
        for (int k = 0; k < cols - c0; ++k) cells.emplace_back(r0, c0 + k);
    }
    grid g(rows, cols);
    for (int r = 0; r < rows; ++r)
        for (int c = 0; c < cols; ++c)
            g.at(r, c) = 20.0 + cheb(r, c, r0, c0);
    g.at(r0, c0) = 1.0;
    for (int k = 1; k < L; ++k) g.at(cells[k].first, cells[k].second) = 5.0 + 0.25 * k;
    double zz = 0.0;
    for (std::size_t k = static_cast<std::size_t>(L); k < cells.size(); ++k)
    {
        const int r = cells[k].first, c = cells[k].second;
        if (r < rows && c < cols)
        {
            g.at(r, c) = zz;
            zz -= 1.0;
        }
    }
    return g;
}

grid f6()
{
    grid g(7, 7);
    for (int r = 0; r < 7; ++r)
        for (int c = 0; c < 7; ++c)
            g.at(r, c) = 10.0 + cheb(r, c, 3, 3);
    g.at(3, 3) = 0.0;
    g.at(3, 4) = 1.0;
    g.at(3, 5) = 2.375;
    g.at(3, 6) = 0.0;
    return g;
}

// Distance (Chebyshev) from a cell to the rectangle [r_lo, r_hi] x [c_lo, c_hi].
int rect_dist(int r, int c, int r_lo, int r_hi, int c_lo, int c_hi)
{
    const int dr = std::max(std::max(r_lo - r, r - r_hi), 0);
    const int dc = std::max(std::max(c_lo - c, c - c_hi), 0);
    return std::max(dr, dc);
}

grid f7()
{
    grid g(11, 11);
    for (int r = 0; r < 11; ++r)
        for (int c = 0; c < 11; ++c)
        {
            const int d = rect_dist(r, c, 3, 7, 3, 7);
            g.at(r, c) = (d == 0) ? 1.0 : 2.0 + d;
        }
    g.at(5, 8) = 3.0; g.at(5, 9) = 0.5; g.at(5, 10) = 0.0;
    return g;
}

grid f8()
{
    grid g(7, 26);
    for (int r = 0; r < 7; ++r)
        for (int c = 0; c < 26; ++c)
        {
            const int d = rect_dist(r, c, 2, 4, 2, 21);
            g.at(r, c) = (d == 0) ? 1.0 : 2.0 + d;
        }
    g.at(3, 22) = 1.0; g.at(3, 23) = 0.5; g.at(3, 24) = 0.25; g.at(3, 25) = 0.0;
    return g;
}

grid f9a()
{
    grid g(7, 7);
    for (int r = 0; r < 7; ++r)
        for (int c = 0; c < 7; ++c)
            g.at(r, c) = 4.0 + cheb(r, c, 3, 3);
    g.at(3, 3) = 1.0;
    g.at(2, 3) = k_nodata;
    return g;
}

grid f9b()
{
    grid g(9, 9);
    for (int r = 0; r < 9; ++r)
        for (int c = 0; c < 9; ++c)
            g.at(r, c) = 4.0 + cheb(r, c, 4, 4);
    g.at(4, 4) = 1.0;
    g.at(4, 7) = k_nodata;
    return g;
}

grid f9c()
{
    const double lut[6] = {1.0, 3.0, 4.0, 5.0, k_nodata, 2.0};
    grid g(11, 11);
    for (int r = 0; r < 11; ++r)
        for (int c = 0; c < 11; ++c)
            g.at(r, c) = lut[cheb(r, c, 5, 5)];
    return g;
}

grid f10()
{
    grid g(7, 7);
    for (int r = 0; r < 7; ++r)
        for (int c = 0; c < 7; ++c)
            g.at(r, c) = 5.0 + cheb(r, c, 1, 3);
    g.at(1, 3) = 1.0;
    g.at(0, 3) = 2.0;
    return g;
}

grid f11()
{
    grid g(7, 13);
    for (int r = 0; r < 7; ++r)
        for (int c = 0; c < 13; ++c)
            g.at(r, c) = 2.0 + std::min(cheb(r, c, 3, 3), cheb(r, c, 3, 9));
    g.at(3, 3) = 1.0;
    g.at(3, 9) = 1.0;
    g.at(3, 12) = 0.0;
    return g;
}

grid f12(double offset)
{
    grid g = f1();
    for (double& v : g.z) v += offset;
    return g;
}

grid f13(const grid& base)
{
    grid g = base;
    for (double& v : g.z) v = 1.125 * v + 0.375;
    return g;
}

grid f14_3x3()
{
    grid g(3, 3, 5.0);
    g.at(1, 1) = 0.0;
    return g;
}

grid f14_1x9()
{
    const double v[9] = {5.0, 4.0, 3.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
    grid g(1, 9);
    for (int c = 0; c < 9; ++c) g.at(0, c) = v[c];
    return g;
}

grid f14_allnd() { return grid(4, 4, k_nodata); }

grid f14_nopit()
{
    grid g(6, 6);
    for (int r = 0; r < 6; ++r)
        for (int c = 0; c < 6; ++c)
            g.at(r, c) = static_cast<double>(r + c);
    return g;
}

// F15 (edge policy): pit (4,10) = 10 beside the east edge. The east border
// column is 12, so the pit drains east over a 2 m rim once the edge is an
// outlet. Westward, cols 9/10 (22/21) and a 40 m divide (cols 7-8) separate it
// from a lowland (z = 1 + col) that drains to the west edge. With the edge NOT
// an outlet (Whitebox) the only escape is through the divide.
grid f15()
{
    grid g(9, 12);
    for (int r = 0; r < 9; ++r)
        for (int c = 0; c < 12; ++c)
        {
            double z;
            if (c <= 6)       z = 1.0 + c;
            else if (c <= 8)  z = 40.0;
            else if (c == 9)  z = 22.0;
            else if (c == 10) z = 21.0;
            else              z = 12.0;
            g.at(r, c) = z;
        }
    g.at(4, 10) = 10.0;
    return g;
}

std::map<std::string, grid> build_fixtures()
{
    std::map<std::string, grid> f;
    f.emplace("F0", f0());
    f.emplace("F1", f1());
    f.emplace("F2", f2());
    f.emplace("F3", f3());
    f.emplace("F4", f4());
    for (int L : {10, 11, 12})
    {
        f.emplace("F5_s" + std::to_string(L), f5(L, false));
        f.emplace("F5_d" + std::to_string(L), f5(L, true));
    }
    f.emplace("F6", f6());
    f.emplace("F7", f7());
    f.emplace("F8", f8());
    f.emplace("F9a", f9a());
    f.emplace("F9b", f9b());
    f.emplace("F9c", f9c());
    f.emplace("F10", f10());
    f.emplace("F11", f11());
    f.emplace("F12_0", f12(-1.0));
    f.emplace("F12_m50", f12(-51.0));
    f.emplace("F12_4000", f12(3999.0));
    f.emplace("F13_F1", f13(f1()));
    f.emplace("F13_F3", f13(f3()));
    f.emplace("F13_F7", f13(f7()));
    f.emplace("F14_3x3", f14_3x3());
    f.emplace("F14_1x9", f14_1x9());
    f.emplace("F14_allnd", f14_allnd());
    f.emplace("F14_nopit", f14_nopit());
    f.emplace("F15", f15());
    return f;
}

std::vector<float> to_f32(const grid& g)
{
    std::vector<float> out(g.z.size());
    for (std::size_t i = 0; i < g.z.size(); ++i) out[i] = static_cast<float>(g.z[i]);
    return out;
}

uint64_t fnv1a64(const std::vector<float>& v)
{
    uint64_t h = 0xcbf29ce484222325ULL;
    for (float x : v)
    {
        uint32_t bits;
        std::memcpy(&bits, &x, sizeof bits);
        for (int b = 0; b < 4; ++b)
        {
            h ^= static_cast<uint8_t>((bits >> (8 * b)) & 0xFFu);
            h *= 0x100000001b3ULL;
        }
    }
    return h;
}

lsdtt::lsdtt_raster_cpu make_raster(const grid& g, double cell_size)
{
    lsdtt::lsdtt_raster_meta meta;
    meta.n_rows        = g.rows;
    meta.n_cols        = g.cols;
    meta.cell_size     = cell_size;
    meta.no_data_value = k_nodata;
    meta.set_default_synthetic_georeference();

    std::vector<float>   data = to_f32(g);
    std::vector<uint8_t> mask(data.size(), 1);
    for (std::size_t i = 0; i < data.size(); ++i)
        if (g.z[i] == k_nodata) mask[i] = 0;

    lsdtt::lsdtt_raster_cpu dem;
    dem.set_from_vectors(std::move(data), std::move(mask), meta);
    return dem;
}


// =============================================================================
// Running the library
// =============================================================================
enum class edge_mode { whitebox, outlet };

// Can this run be executed with the API the library currently provides?
bool runnable(const wb_run& run, edge_mode mode)
{
#if LSDTT_BREACH_LC_HAS_OPTIONS
    (void)run; (void)mode;
    return true;
#else
    // Positional overload only: it has no edge policy, and its float
    // increment cannot carry Whitebox's double auto value.
    return mode == edge_mode::whitebox && run.explicit_increment;
#endif
}

#if LSDTT_BREACH_LC_HAS_OPTIONS
lsdtt::breach_least_cost_options make_options(const wb_run& run, edge_mode mode)
{
    lsdtt::breach_least_cost_options o;
    static_assert(std::is_same_v<decltype(o.flat_increment), double>,
                  "flat_increment must be double: the auto runs pass Whitebox's exact increment");
    o.max_dist       = run.max_dist;
    o.max_cost       = run.max_cost;
    o.minimize_dist  = run.minimize_dist;
    o.flat_increment = run.small;
    o.fill_remainder = false;
    o.edge_policy    = (mode == edge_mode::whitebox) ? lsdtt::breach_lc_edge_policy::whitebox
                                                     : lsdtt::breach_lc_edge_policy::outlet;
    return o;
}
#endif

lsdtt::lsdtt_raster_cpu run_breach(const lsdtt::lsdtt_raster_cpu& dem, const wb_run& run,
                                   edge_mode mode)
{
    cout_silencer quiet;
#if LSDTT_BREACH_LC_HAS_OPTIONS
    return lsdtt::breach_depressions_least_cost(dem, make_options(run, mode));
#else
    (void)mode;
    const float max_cost = std::isinf(run.max_cost)
                         ? (std::numeric_limits<float>::max)()
                         : static_cast<float>(run.max_cost);
    return lsdtt::breach_depressions_least_cost(
        dem, run.max_dist, max_cost, /*min_dist=*/run.minimize_dist,
        /*flat_increment=*/static_cast<float>(run.small), /*fill_remainder=*/false);
#endif
}

// Expected float32 output: the input with Whitebox's changed cells applied.
std::vector<float> expected_output(const grid& g, const wb_run& run)
{
    std::vector<float> e = to_f32(g);
    for (int i = 0; i < run.n_changes; ++i)
        e[static_cast<std::size_t>(run.changes[i].r) * g.cols + run.changes[i].c] = run.changes[i].z;
    return e;
}

struct compare_result
{
    bool        equal = true;
    int         n_diff = 0;
    std::string detail;
};

std::string fmt_float(float v)
{
    std::ostringstream s;
    s << std::setprecision(9) << v;
    return s.str();
}

// Bitwise on every valid cell; the mask must be identical.
compare_result compare_bitwise(const grid& g, const lsdtt::lsdtt_raster_cpu& out,
                               const std::vector<float>& expected)
{
    compare_result res;
    const auto& od = out.get_data();
    const auto& om = out.get_mask_array().get_vector();
    if (od.size() != expected.size() || om.size() != expected.size())
    {
        res.equal = false;
        res.detail = "output size " + std::to_string(od.size()) + " != " + std::to_string(expected.size());
        return res;
    }
    std::ostringstream first;
    for (std::size_t i = 0; i < expected.size(); ++i)
    {
        const bool valid_in = g.z[i] != k_nodata;
        bool bad = false;
        if ((om[i] != 0) != valid_in) bad = true;
        else if (valid_in)
        {
            uint32_t a, b;
            std::memcpy(&a, &od[i], sizeof a);
            std::memcpy(&b, &expected[i], sizeof b);
            bad = (a != b);
        }
        if (bad)
        {
            if (res.n_diff < 6)
            {
                const int r = static_cast<int>(i) / g.cols, c = static_cast<int>(i) % g.cols;
                first << " (" << r << "," << c << ") got " << fmt_float(od[i])
                      << (om[i] ? "" : "[nd]") << " want " << fmt_float(expected[i]) << ";";
            }
            ++res.n_diff;
            res.equal = false;
        }
    }
    if (!res.equal) res.detail = std::to_string(res.n_diff) + " cells differ:" + first.str();
    return res;
}

#if LSDTT_BREACH_LC_HAS_OPTIONS
// Helpers for the outlet-mode property checks.

// Changed valid cells as (index -> output value).
std::map<std::size_t, float> changed_cells(const grid& g, const lsdtt::lsdtt_raster_cpu& out)
{
    std::map<std::size_t, float> ch;
    const auto in = to_f32(g);
    const auto& od = out.get_data();
    for (std::size_t i = 0; i < in.size(); ++i)
        if (g.z[i] != k_nodata && od[i] != in[i]) ch[i] = od[i];
    return ch;
}

float max_cut(const grid& g, const lsdtt::lsdtt_raster_cpu& out)
{
    float m = 0.0f;
    const auto in = to_f32(g);
    const auto& od = out.get_data();
    for (std::size_t i = 0; i < in.size(); ++i)
        if (g.z[i] != k_nodata) m = std::max(m, in[i] - od[i]);
    return m;
}

float max_cut_expected(const grid& g, const wb_run& run)
{
    float m = 0.0f;
    const auto in = to_f32(g);
    for (int i = 0; i < run.n_changes; ++i)
    {
        const std::size_t k = static_cast<std::size_t>(run.changes[i].r) * g.cols + run.changes[i].c;
        m = std::max(m, in[k] - run.changes[i].z);
    }
    return m;
}

#endif

std::vector<const wb_run*> runs_of(const std::string& fixture)
{
    std::vector<const wb_run*> v;
    for (const wb_run& r : k_wb_runs)
        if (fixture == r.fixture) v.push_back(&r);
    return v;
}


// =============================================================================
// D0 -- the in-process fixtures are the inputs Whitebox was run on
// =============================================================================
void test_fixture_digests(const std::map<std::string, grid>& fx)
{
    for (const fixture_digest& d : k_fixture_digests)
    {
        const auto it = fx.find(d.name);
        if (it == fx.end())
        {
            check(std::string("D0 fixture builder exists: ") + d.name, false, "no builder");
            continue;
        }
        const grid& g = it->second;
        const uint64_t h = fnv1a64(to_f32(g));
        std::ostringstream det;
        det << "shape " << g.rows << "x" << g.cols << " (want " << d.rows << "x" << d.cols
            << "), fnv1a64 0x" << std::hex << h << " (want 0x" << d.fnv1a64 << ")";
        check(std::string("D0 fixture equals the Whitebox input: ") + d.name,
              g.rows == d.rows && g.cols == d.cols && h == d.fnv1a64, det.str());
    }
}


// =============================================================================
// W -- bitwise against float32(Whitebox), edge_policy = whitebox
// =============================================================================
struct fixture_summary
{
    int n_runs = 0;
    int n_run = 0;
    int n_pass = 0;
    std::vector<std::string> failed_keys;
};

std::map<std::string, fixture_summary> g_summary;

void test_whitebox_literals(const std::map<std::string, grid>& fx)
{
    for (const fixture_digest& d : k_fixture_digests)
    {
        const grid& g = fx.at(d.name);
        fixture_summary& s = g_summary[d.name];
        std::ostringstream det;
        for (const wb_run* run : runs_of(d.name))
        {
            ++s.n_runs;
            if (!runnable(*run, edge_mode::whitebox)) continue;
            ++s.n_run;
            compare_result cr;
            try
            {
                const auto out = run_breach(make_raster(g, run->cell_size), *run, edge_mode::whitebox);
                cr = compare_bitwise(g, out, expected_output(g, *run));
            }
            catch (const std::exception& e)
            {
                cr.equal = false;
                cr.detail = std::string("threw: ") + e.what();
            }
            if (cr.equal) ++s.n_pass;
            else
            {
                s.failed_keys.push_back(run->key);
                if (s.failed_keys.size() <= 2) det << run->key << ": " << cr.detail << " | ";
            }
        }
        if (s.n_run == 0)
        {
            not_run(std::string("W bitwise vs Whitebox 2.0.6 (edge whitebox): ") + d.name,
                    "no run executable with this API (auto increment needs the options overload)");
            continue;
        }
        std::ostringstream name;
        name << "W bitwise vs Whitebox 2.0.6 (edge whitebox): " << d.name << " "
             << s.n_pass << "/" << s.n_run << " runs";
        if (s.n_run < s.n_runs) name << " (" << (s.n_runs - s.n_run) << " not runnable yet)";
        check(name.str(), s.n_pass == s.n_run, det.str());
    }
}


// =============================================================================
// C -- Whitebox's internal counters
// =============================================================================
void test_counters(const std::map<std::string, grid>& fx)
{
#if LSDTT_BREACH_LC_HAS_OPTIONS && LSDTT_BREACH_LC_HAS_DIAGNOSTICS
    for (const fixture_digest& d : k_fixture_digests)
    {
        const grid& g = fx.at(d.name);
        int n_ok = 0, n = 0;
        std::ostringstream det;
        for (const wb_run* run : runs_of(d.name))
        {
            ++n;
            lsdtt::breach_least_cost_diagnostics diag;
            {
                cout_silencer quiet;
                (void)lsdtt::breach_depressions_least_cost(
                    make_raster(g, run->cell_size), make_options(*run, edge_mode::whitebox), &diag);
            }
            const wb_counters& w = run->counters;
            const bool ok = diag.pits_detected == w.pits_detected &&
                            diag.pits_searched == w.pits_searched &&
                            diag.pits_skipped_not_pit == w.pits_skipped_not_pit &&
                            diag.pits_solved == w.pits_solved &&
                            diag.pits_unsolved_exhausted == w.pits_unsolved_exhausted &&
                            diag.pits_unsolved_max_cost == w.pits_unsolved_max_cost;
            if (ok) ++n_ok;
            else det << run->key << " ";
        }
        check("C counters equal Whitebox's: " + std::string(d.name) + " " +
              std::to_string(n_ok) + "/" + std::to_string(n), n_ok == n, det.str());
    }
#else
    (void)fx;
    not_run("C counters equal Whitebox's (all fixtures)",
            "the library exposes no breach diagnostics yet");
#endif
}


// =============================================================================
// E -- parameter checking
// =============================================================================
void test_errors(const std::map<std::string, grid>& fx)
{
    const auto dem = make_raster(fx.at("F1"), 1.0);
    auto attempt = [&](int max_dist, std::string& msg) -> bool
    {
        cout_silencer quiet;
        try
        {
            (void)lsdtt::breach_depressions_least_cost(
                dem, max_dist, (std::numeric_limits<float>::max)(), false, 0.125f, false);
        }
        catch (const std::exception& e)
        {
            msg = e.what();
            return true;
        }
        return false;
    };

    for (int bad : {0, -1, 32768})
    {
        std::string msg;
        const bool threw = attempt(bad, msg);
        check("E max_dist " + std::to_string(bad) + " is rejected with a message naming max_dist",
              threw && msg.find("max_dist") != std::string::npos,
              threw ? "message: " + msg : "did not throw");
    }
    std::string msg;
    const bool threw = attempt(32767, msg);
    check("E max_dist 32767 is accepted", !threw, "threw: " + msg);

#if LSDTT_BREACH_LC_HAS_OPTIONS
    for (float bad : {-1.0e-4f, std::numeric_limits<float>::quiet_NaN(),
                      std::numeric_limits<float>::infinity()})
    {
        lsdtt::breach_least_cost_options o;
        o.fill_min_slope = bad;
        std::string m;
        bool t = false;
        {
            cout_silencer quiet;
            try { (void)lsdtt::breach_depressions_least_cost(dem, o); }
            catch (const std::exception& e) { m = e.what(); t = true; }
        }
        check("E fill_min_slope " + fmt_float(bad) + " is rejected with a message naming fill_min_slope",
              t && m.find("fill_min_slope") != std::string::npos,
              t ? "message: " + m : "did not throw");
    }

    // barnes2014b_gradient takes its gradient from fill_min_slope, so 0 is
    // refused -- but only when the remainder fill can run.
    {
        lsdtt::breach_least_cost_options o;
        o.fill_flat_resolution = lsdtt::flat_resolution_choice::barnes2014b_gradient;
        o.fill_min_slope       = 0.0f;
        std::string m;
        bool t = false;
        {
            cout_silencer quiet;
            try { (void)lsdtt::breach_depressions_least_cost(dem, o); }
            catch (const std::exception& e) { m = e.what(); t = true; }
        }
        check("E barnes2014b_gradient with fill_min_slope 0 is rejected, naming both",
              t && m.find("fill_min_slope") != std::string::npos &&
                  m.find("barnes2014b_gradient") != std::string::npos,
              t ? "message: " + m : "did not throw");

        o.fill_remainder = false;
        t = false;
        {
            cout_silencer quiet;
            try { (void)lsdtt::breach_depressions_least_cost(dem, o); }
            catch (const std::exception& e) { m = e.what(); t = true; }
        }
        check("E barnes2014b_gradient with fill_min_slope 0 is accepted when fill_remainder is false",
              !t, "threw: " + m);
    }

    // The edge-policy spelling shared by the CLIs and the binding.
    {
        bool ok = lsdtt::parse_breach_lc_edge_policy("outlet") == lsdtt::breach_lc_edge_policy::outlet &&
                  lsdtt::parse_breach_lc_edge_policy("whitebox") == lsdtt::breach_lc_edge_policy::whitebox &&
                  std::string(lsdtt::breach_lc_edge_policy_name(lsdtt::breach_lc_edge_policy::outlet)) == "outlet" &&
                  std::string(lsdtt::breach_lc_edge_policy_name(lsdtt::breach_lc_edge_policy::whitebox)) == "whitebox";
        check("E edge policy parses outlet / whitebox and names them back", ok);
        std::string m;
        bool t = false;
        try { (void)lsdtt::parse_breach_lc_edge_policy("Outlet", "breach_lc_edge_policy"); }
        catch (const std::invalid_argument& e) { m = e.what(); t = true; }
        check("E edge policy 'Outlet' is rejected with the exact message",
              t && m == "breach_lc_edge_policy must be one of outlet, whitebox; got 'Outlet'",
              t ? "message: " + m : "did not throw");
    }
#endif
}


// =============================================================================
// A -- lsdtt3's automatic flat increment
// =============================================================================
#if LSDTT_BREACH_LC_HAS_OPTIONS && LSDTT_BREACH_LC_HAS_DIAGNOSTICS
double auto_increment_of(const grid& g)
{
    lsdtt::breach_least_cost_options o;   // flat_increment 0 -> auto
    lsdtt::breach_least_cost_diagnostics diag;
    cout_silencer quiet;
    (void)lsdtt::breach_depressions_least_cost(make_raster(g, 1.0), o, &diag);
    return diag.flat_increment;
}
#endif

void test_auto_increment(const std::map<std::string, grid>& fx)
{
#if LSDTT_BREACH_LC_HAS_OPTIONS && LSDTT_BREACH_LC_HAS_DIAGNOSTICS
    // 2 x ULP_f32(a), a = max(|z|max, 1), ULP_f32(a) = nextafterf(a, +inf) - a.
    // The expected values are written out, not computed with the rule.
    struct auto_case { const char* name; grid g; double want; };
    const std::vector<auto_case> cases = {
        // F1 - 1: z in [-1, 4]; |z|max 4 -> ULP 2^-21.
        {"F12_0 (|z|max 4)",          fx.at("F12_0"),    0x1p-20},
        // F1 - 51: z in [-51, -46]; |z|max 51 (the most negative) -> ULP 2^-18.
        {"F12_m50 (|z|max 51)",       fx.at("F12_m50"),  0x1p-17},
        // F1 + 3999: z in [3999, 4004] -> ULP 2^-12.
        {"F12_4000 (|z|max 4004)",    fx.at("F12_4000"), 0x1p-11},
        // Flat rasters: no special case.
        {"flat 100",                  grid(6, 6, 100.0), 0x1p-16},
        {"flat 0 (floor of 1)",       grid(6, 6, 0.0),   0x1p-22},
        // A power of two takes the gap ABOVE it: ULP(1024) = 2^-13.
        {"flat 1024 (binade above)",  grid(6, 6, 1024.0), 0x1p-12},
        // All nodata: |z|max 0 -> the floor of 1.
        {"all nodata (floor of 1)",   fx.at("F14_allnd"), 0x1p-22},
    };
    for (const auto_case& c : cases)
    {
        const double got = auto_increment_of(c.g);
        std::ostringstream det;
        det << std::hexfloat << "got " << got << " want " << c.want;
        check(std::string("A auto flat_increment = 2 x ULP_f32: ") + c.name, got == c.want, det.str());
    }
    // An explicit increment is used as given.
    {
        lsdtt::breach_least_cost_options o;
        o.flat_increment = 0.125;
        lsdtt::breach_least_cost_diagnostics diag;
        {
            cout_silencer quiet;
            (void)lsdtt::breach_depressions_least_cost(make_raster(fx.at("F12_4000"), 1.0), o, &diag);
        }
        check("A explicit flat_increment 0.125 is used as given", diag.flat_increment == 0.125);
    }
#else
    (void)fx;
    not_run("A auto flat_increment", "breach_least_cost_options / diagnostics not available");
#endif
}


// =============================================================================
// R -- routability after the remainder fill, library defaults
// =============================================================================
#if LSDTT_BREACH_LC_HAS_OPTIONS
// Interior cells (all 8 neighbours in the raster and valid) with no strictly
// lower neighbour.
int count_without_descent(const lsdtt::lsdtt_raster_cpu& out)
{
    const auto& od = out.get_data();
    const auto& om = out.get_mask_array().get_vector();
    const int rows = static_cast<int>(out.get_meta().n_rows);
    const int cols = static_cast<int>(out.get_meta().n_cols);
    int n = 0;
    for (int r = 1; r < rows - 1; ++r)
        for (int c = 1; c < cols - 1; ++c)
        {
            const std::size_t i = static_cast<std::size_t>(r) * cols + c;
            if (!om[i]) continue;
            bool interior = true, lower = false;
            for (int dr = -1; dr <= 1; ++dr)
                for (int dc = -1; dc <= 1; ++dc)
                {
                    if (dr == 0 && dc == 0) continue;
                    const std::size_t j = static_cast<std::size_t>(r + dr) * cols + (c + dc);
                    if (!om[j]) interior = false;
                    else if (od[j] < od[i]) lower = true;
                }
            if (interior && !lower) ++n;
        }
    return n;
}
#endif

void test_routability(const std::map<std::string, grid>& fx)
{
#if LSDTT_BREACH_LC_HAS_OPTIONS && LSDTT_BREACH_LC_HAS_DIAGNOSTICS
    for (const char* name : {"F1", "F7", "F8", "F10", "F15"})
    {
        const grid& g = fx.at(name);
        for (edge_mode mode : {edge_mode::outlet, edge_mode::whitebox})
        {
            for (double cs : {1.0, 30.0})
            {
                lsdtt::breach_least_cost_options o;   // defaults
                o.edge_policy = (mode == edge_mode::whitebox) ? lsdtt::breach_lc_edge_policy::whitebox
                                                              : lsdtt::breach_lc_edge_policy::outlet;
                lsdtt::breach_least_cost_diagnostics diag;
                lsdtt::lsdtt_raster_cpu out;
                {
                    cout_silencer quiet;
                    out = lsdtt::breach_depressions_least_cost(make_raster(g, cs), o, &diag);
                }
                const int n_bad = count_without_descent(out);
                std::ostringstream lab, det;
                lab << "R defaults, " << (mode == edge_mode::whitebox ? "whitebox" : "outlet")
                    << ", cs " << cs << ": every interior cell descends: " << name
                    << " (breach left " << diag.cells_without_descent << ", remainder fill "
                    << (diag.remainder_filled ? "ran" : "not needed") << ")";
                det << n_bad << " interior cells without a strictly lower neighbour";
                check(lab.str(), n_bad == 0, det.str());
            }
        }
    }

    // A larger remainder: F8 at increment 0.125 leaves 16 pits unsolved
    // (relief exhaustion, whitebox mode); the default fill must still leave
    // every interior cell descending.
    {
        lsdtt::breach_least_cost_options o;
        o.flat_increment = 0.125;
        o.edge_policy    = lsdtt::breach_lc_edge_policy::whitebox;
        lsdtt::breach_least_cost_diagnostics diag;
        lsdtt::lsdtt_raster_cpu out;
        {
            cout_silencer quiet;
            out = lsdtt::breach_depressions_least_cost(make_raster(fx.at("F8"), 1.0), o, &diag);
        }
        const int n_bad = count_without_descent(out);
        check("R F8 whitebox, increment 0.125 (" + std::to_string(diag.pits_unsolved_exhausted) +
              " unsolved): every interior cell descends after the remainder fill",
              diag.remainder_filled && n_bad == 0,
              std::to_string(n_bad) + " interior cells without a strictly lower neighbour");
    }

    // Control: with fill_min_slope 0 (the flat fill used before) the F10
    // whitebox remainder is NOT routable, so the R checks can fail.
    {
        lsdtt::breach_least_cost_options o;
        o.edge_policy    = lsdtt::breach_lc_edge_policy::whitebox;
        o.fill_min_slope = 0.0f;
        lsdtt::lsdtt_raster_cpu out;
        {
            cout_silencer quiet;
            out = lsdtt::breach_depressions_least_cost(make_raster(fx.at("F10"), 1.0), o);
        }
        const int n_bad = count_without_descent(out);
        check("R control: F10 whitebox with fill_min_slope 0 leaves cells without descent",
              n_bad > 0, "0 cells without descent: the R checks cannot detect a flat fill");
    }

    // The remainder fill honours fill_flat_resolution: with either barnes2014b
    // option every interior cell still descends (whitebox mode, so F10 and
    // F15 leave a remainder for the fill), and on F8 at increment 0.125 (16
    // unsolved pits) the result differs from inline_gradient, so the option
    // reaches the fill.
    for (lsdtt::flat_resolution_choice fr : {lsdtt::flat_resolution_choice::barnes2014b_gradient,
                                             lsdtt::flat_resolution_choice::barnes2014b_float_step})
    {
        const std::string fr_name = lsdtt::flat_resolution_choice_name(fr);
        for (const char* name : {"F1", "F7", "F8", "F10", "F15"})
        {
            for (double cs : {1.0, 30.0})
            {
                lsdtt::breach_least_cost_options o;
                o.edge_policy          = lsdtt::breach_lc_edge_policy::whitebox;
                o.fill_flat_resolution = fr;
                lsdtt::breach_least_cost_diagnostics diag;
                lsdtt::lsdtt_raster_cpu out;
                {
                    cout_silencer quiet;
                    out = lsdtt::breach_depressions_least_cost(make_raster(fx.at(name), cs), o, &diag);
                }
                const int n_bad = count_without_descent(out);
                std::ostringstream lab;
                lab << "R " << fr_name << ", whitebox, cs " << cs
                    << ": every interior cell descends: " << name << " (remainder fill "
                    << (diag.remainder_filled ? "ran" : "not needed") << ")";
                check(lab.str(), n_bad == 0,
                      std::to_string(n_bad) + " interior cells without a strictly lower neighbour");
            }
        }

        lsdtt::breach_least_cost_options o;
        o.flat_increment = 0.125;
        o.edge_policy    = lsdtt::breach_lc_edge_policy::whitebox;
        lsdtt::lsdtt_raster_cpu inline_out, fr_out;
        lsdtt::breach_least_cost_diagnostics diag;
        {
            cout_silencer quiet;
            inline_out = lsdtt::breach_depressions_least_cost(make_raster(fx.at("F8"), 1.0), o);
            o.fill_flat_resolution = fr;
            fr_out = lsdtt::breach_depressions_least_cost(make_raster(fx.at("F8"), 1.0), o, &diag);
        }
        const int n_bad = count_without_descent(fr_out);
        check("R F8 whitebox, increment 0.125, " + fr_name +
              ": remainder fill ran, every interior cell descends, result differs from inline_gradient",
              diag.remainder_filled && n_bad == 0 && fr_out.get_data() != inline_out.get_data(),
              std::to_string(n_bad) + " interior cells without a strictly lower neighbour");
    }
#else
    (void)fx;
    not_run("R routability after the remainder fill", "breach_least_cost_options not available");
#endif
}


// =============================================================================
// O -- edge_policy = outlet (property tests; no external oracle)
// =============================================================================
void test_outlet_policy(const std::map<std::string, grid>& fx)
{
#if LSDTT_BREACH_LC_HAS_OPTIONS
    // O0: outlet is the default.
    {
        const lsdtt::breach_least_cost_options o{};
        check("O0 breach_least_cost_options default edge_policy is outlet",
              o.edge_policy == lsdtt::breach_lc_edge_policy::outlet);
        // max_dist defaults by edge policy: 1000 outlet (bridges large water
        // bodies), 100 whitebox (Whitebox's default); the struct carries the
        // outlet one.
        check("O0 default max_dist: 1000 outlet, 100 whitebox; the struct default is the outlet one",
              lsdtt::default_breach_lc_max_dist(lsdtt::breach_lc_edge_policy::outlet) == 1000 &&
                  lsdtt::default_breach_lc_max_dist(lsdtt::breach_lc_edge_policy::whitebox) == 100 &&
                  o.max_dist == 1000,
              "struct max_dist " + std::to_string(o.max_dist));
    }

    // O1: where Whitebox's searches never expand a border cell, outlet mode
    // cannot differ from Whitebox.
    for (const fixture_digest& d : k_fixture_digests)
    {
        const grid& g = fx.at(d.name);
        int n = 0, n_ok = 0;
        std::ostringstream det;
        for (const wb_run* run : runs_of(d.name))
        {
            if (run->border_popped) continue;
            ++n;
            const auto out = run_breach(make_raster(g, run->cell_size), *run, edge_mode::outlet);
            const compare_result cr = compare_bitwise(g, out, expected_output(g, *run));
            if (cr.equal) ++n_ok;
            else if (n - n_ok <= 2) det << run->key << ": " << cr.detail << " | ";
        }
        if (n == 0) continue;
        check("O1 outlet == Whitebox where no search reaches the border: " + std::string(d.name) +
              " " + std::to_string(n_ok) + "/" + std::to_string(n), n_ok == n, det.str());
    }

    // O2: F10 -- the pit's only low rim is the border cell (0,3) = 2 > pit.
    // Whitebox leaves it unsolved (W covers that); outlet mode cuts the border
    // cell one step below the raised pit and drains off the grid.
    {
        const grid& g = fx.at("F10");
        for (const wb_run* run : runs_of("F10"))
        {
            if (!run->explicit_increment) continue;   // exact values need the dyadic increment
            const auto out = run_breach(make_raster(g, run->cell_size), *run, edge_mode::outlet);
            const auto ch = changed_cells(g, out);
            const std::size_t pit = 1 * static_cast<std::size_t>(g.cols) + 3;
            const std::size_t rim = 0 * static_cast<std::size_t>(g.cols) + 3;
            const bool ok = ch.size() == 2 && ch.count(pit) && ch.count(rim) &&
                            ch.at(pit) == 1.875f && ch.at(rim) == 1.75f;
            std::ostringstream det;
            for (const auto& kv : ch)
                det << "(" << kv.first / g.cols << "," << kv.first % g.cols << ")=" << fmt_float(kv.second) << " ";
            check(std::string("O2 F10 outlet mode: pit 1.875, border cell (0,3) cut to 1.75, nothing else: ") +
                  run->key, ok, det.str());
        }
    }

    // O3: F15 -- outlet mode drains east to the edge: only the pit (raised to
    // 11.875) and ONE east-border cell beside it (cut to 11.75) change, and the
    // maximum cut is below Whitebox's divide crossing.
    {
        const grid& g = fx.at("F15");
        for (const wb_run* run : runs_of("F15"))
        {
            if (!run->explicit_increment) continue;
            const auto out = run_breach(make_raster(g, run->cell_size), *run, edge_mode::outlet);
            const auto ch = changed_cells(g, out);
            const std::size_t pit = 4 * static_cast<std::size_t>(g.cols) + 10;
            int n_border = 0;
            bool other = false;
            for (const auto& kv : ch)
            {
                const int r = static_cast<int>(kv.first) / g.cols, c = static_cast<int>(kv.first) % g.cols;
                if (kv.first == pit) continue;
                if (c == g.cols - 1 && r >= 3 && r <= 5 && kv.second == 11.75f) ++n_border;
                else other = true;
            }
            const float cut_outlet = max_cut(g, out);
            const float cut_wb     = max_cut_expected(g, *run);
            const bool ok = ch.count(pit) && ch.at(pit) == 11.875f && n_border == 1 && !other &&
                            cut_outlet < cut_wb;
            std::ostringstream det;
            for (const auto& kv : ch)
                det << "(" << kv.first / g.cols << "," << kv.first % g.cols << ")=" << fmt_float(kv.second) << " ";
            det << "max cut outlet " << cut_outlet << " vs whitebox " << cut_wb;
            check(std::string("O3 F15 outlet mode drains to the east edge, max cut < Whitebox's: ") +
                  run->key, ok, det.str());
        }
    }
#else
    (void)fx;
    for (const char* name : {"O0 default edge_policy is outlet",
                             "O1 outlet == Whitebox where no search reaches the border",
                             "O2 F10 outlet mode solved through the border cell",
                             "O3 F15 outlet mode drains to the east edge"})
        not_run(name, "breach_least_cost_options / edge_policy not available yet");
#endif
}


// =============================================================================
// RT -- real terrain: the committed lakeland crop
// =============================================================================
#if LSDTT_BREACH_LC_HAS_OPTIONS && LSDTT_BREACH_LC_HAS_DIAGNOSTICS
std::string hex64(uint64_t h)
{
    std::ostringstream o;
    o << "0x" << std::hex << std::setw(16) << std::setfill('0') << h;
    return o.str();
}

void run_real_terrain_checks(const std::string& path)
{
    std::cout << "=== Real terrain: " << path << " ===\n";

    lsdtt::lsdtt_raster_cpu dem;
    try
    {
        lsdtt::lsdtt_data_loader loader;
        loader.open_raster(path);
        lsdtt::lsdtt_raster_tile tile = loader.read_full();
        dem.set_from_vectors(std::move(tile.data), std::move(tile.mask), tile.meta);
    }
    catch (const std::exception& e)
    {
        check("RT0 the crop loads", false, e.what());
        return;
    }
    const lsdtt::lsdtt_raster_meta& meta = dem.get_meta();
    std::int64_t n_valid = 0;
    for (const std::uint8_t v : dem.get_mask()) n_valid += v;
    check("RT0 the crop is 128 x 128 at 30 m with every cell valid",
          meta.n_rows == 128 && meta.n_cols == 128 && meta.cell_size == 30.0 && n_valid == 128 * 128,
          std::to_string(meta.n_rows) + " x " + std::to_string(meta.n_cols) + ", " +
              std::to_string(n_valid) + " valid");
    // The float32 input Whitebox was run on (predictor-free copy of this file).
    const uint64_t k_input_digest = 0x5ca1d7cd28291e35ULL;
    check("RT0 the crop's float32 values are the Whitebox input: digest " + hex64(k_input_digest),
          fnv1a64(dem.get_data()) == k_input_digest, "got " + hex64(fnv1a64(dem.get_data())));

    // RT1 -- matched parameters vs Whitebox 2.0.6, frozen digests.
    struct rt1_case
    {
        int          max_dist;
        bool         minimize_dist;
        uint64_t     whitebox_digest;   // FNV-1a-64 of float32(Whitebox f64 output)
        std::int64_t detected, searched, skipped, solved, unsolved;
    };
    const std::vector<rt1_case> cases = {
        {100,  false, 0xce3742d1cccf9d7fULL, 822, 564, 258, 564, 0},
        {1000, false, 0xce3742d1cccf9d7fULL, 822, 564, 258, 564, 0},
        {100,  true,  0xd3377851c2c5f3fdULL, 822, 614, 208, 614, 0},
        {1000, true,  0xd3377851c2c5f3fdULL, 822, 614, 208, 614, 0},
    };
    for (const rt1_case& c : cases)
    {
        lsdtt::breach_least_cost_options o;
        o.max_dist       = c.max_dist;
        o.minimize_dist  = c.minimize_dist;
        o.flat_increment = static_cast<double>(0.001f);
        o.fill_remainder = false;
        o.edge_policy    = lsdtt::breach_lc_edge_policy::whitebox;
        lsdtt::breach_least_cost_diagnostics diag;
        lsdtt::lsdtt_raster_cpu out;
        {
            cout_silencer quiet;
            out = lsdtt::breach_depressions_least_cost(dem, o, &diag);
        }
        const uint64_t h = fnv1a64(out.get_data());
        std::ostringstream lab, det;
        lab << "RT1 crop, whitebox mode, max_dist " << c.max_dist << ", minimize_dist "
            << (c.minimize_dist ? "true" : "false")
            << ": float32 output == Whitebox 2.0.6 (digest " << hex64(c.whitebox_digest) << ")";
        det << "got " << hex64(h);
        check(lab.str(), h == c.whitebox_digest, det.str());

        const std::int64_t unsolved = diag.pits_unsolved_exhausted + diag.pits_unsolved_max_cost;
        std::ostringstream cl, cd;
        cl << "RT1 crop counters, max_dist " << c.max_dist << ", minimize_dist "
           << (c.minimize_dist ? "true" : "false") << ": detected/searched/skipped/solved/unsolved = "
           << c.detected << "/" << c.searched << "/" << c.skipped << "/" << c.solved << "/" << c.unsolved;
        cd << "got " << diag.pits_detected << "/" << diag.pits_searched << "/"
           << diag.pits_skipped_not_pit << "/" << diag.pits_solved << "/" << unsolved;
        check(cl.str(),
              diag.pits_detected == c.detected && diag.pits_searched == c.searched &&
                  diag.pits_skipped_not_pit == c.skipped && diag.pits_solved == c.solved &&
                  unsolved == c.unsolved,
              cd.str());
    }

    // RT2 / RT3 -- routable outputs.
    // On the crop every pit breaches within either default max_dist, so the
    // remainder fill does not run; the max_dist 3 rows leave pits unsolved
    // (required below) so each fill is exercised on real terrain.
    struct rt2_case
    {
        const char*                   name;
        int                           max_dist;
        lsdtt::breach_lc_edge_policy  edge;
        lsdtt::flat_resolution_choice flat;
    };
    const std::vector<rt2_case> rcases = {
        {"defaults (outlet, max_dist 1000, inline_gradient)", 1000,
         lsdtt::breach_lc_edge_policy::outlet, lsdtt::flat_resolution_choice::inline_gradient},
        {"whitebox, max_dist 100 (its default), inline_gradient", 100,
         lsdtt::breach_lc_edge_policy::whitebox, lsdtt::flat_resolution_choice::inline_gradient},
        {"outlet, max_dist 3, inline_gradient", 3,
         lsdtt::breach_lc_edge_policy::outlet, lsdtt::flat_resolution_choice::inline_gradient},
        {"outlet, max_dist 3, barnes2014b_gradient", 3,
         lsdtt::breach_lc_edge_policy::outlet, lsdtt::flat_resolution_choice::barnes2014b_gradient},
        {"outlet, max_dist 3, barnes2014b_float_step", 3,
         lsdtt::breach_lc_edge_policy::outlet, lsdtt::flat_resolution_choice::barnes2014b_float_step},
        {"whitebox, max_dist 3, barnes2014b_float_step", 3,
         lsdtt::breach_lc_edge_policy::whitebox, lsdtt::flat_resolution_choice::barnes2014b_float_step},
    };
    for (const rt2_case& c : rcases)
    {
        lsdtt::breach_least_cost_options o;   // defaults otherwise
        o.max_dist             = c.max_dist;
        o.edge_policy          = c.edge;
        o.fill_flat_resolution = c.flat;
        lsdtt::breach_least_cost_diagnostics diag;
        lsdtt::lsdtt_raster_cpu out;
        {
            cout_silencer quiet;
            out = lsdtt::breach_depressions_least_cost(dem, o, &diag);
        }
        const int n_bad = count_without_descent(out);
        check(std::string("RT2 crop, ") + c.name + ": every interior cell descends (breach left " +
                  std::to_string(diag.cells_without_descent) + ", remainder fill " +
                  (diag.remainder_filled ? "ran" : "not needed") + ")",
              n_bad == 0 && (c.max_dist != 3 || diag.remainder_filled),
              std::to_string(n_bad) + " interior cells without a strictly lower neighbour" +
                  (c.max_dist == 3 && !diag.remainder_filled
                       ? "; max_dist 3 left nothing for the fill to do" : ""));

        std::string err;
        try
        {
            lsdtt::lsdtt_flowinfo_options fo;
            fo.boundary_conditions = "oooo";
            const lsdtt::lsdtt_flowinfo fi(out.get_data_array(), out.get_mask_array(), out.meta(), fo);
            fi.validate();
        }
        catch (const std::exception& e)
        {
            err = e.what();
        }
        check(std::string("RT3 crop, ") + c.name + ": D8 flowinfo builds and validate() passes",
              err.empty(), err);
    }
}
#endif

} // namespace


int main(int argc, char** argv)
{
    // Real-terrain mode (unit_breach_least_cost_real_terrain): the committed
    // crop only, so the hermetic default run stays free of file input.
    if (argc == 3 && std::string(argv[1]) == "--real-terrain")
    {
#if LSDTT_BREACH_LC_HAS_OPTIONS && LSDTT_BREACH_LC_HAS_DIAGNOSTICS
        run_real_terrain_checks(argv[2]);
#else
        not_run("RT real terrain", "breach_least_cost_options / diagnostics not available");
#endif
        std::cout << "\n" << g_n_checks << " checks, " << g_n_failed << " failed\n";
        std::cout << "RESULT: " << (g_n_failed == 0 && g_n_checks > 0 ? "PASS" : "FAIL") << "\n";
        return g_n_failed == 0 && g_n_checks > 0 ? 0 : 1;
    }
    if (argc != 1)
    {
        std::cerr << "usage: lsdtt-test-breach-least-cost [--real-terrain <crop.tif>]\n";
        return 2;
    }

    std::cout << "unit_breach_least_cost: least-cost breaching vs frozen Whitebox 2.0.6\n"
              << "options overload: " << (LSDTT_BREACH_LC_HAS_OPTIONS ? "yes" : "no (positional API)")
              << ", diagnostics: " << (LSDTT_BREACH_LC_HAS_DIAGNOSTICS ? "yes" : "no") << "\n\n";

    const std::map<std::string, grid> fx = build_fixtures();

    test_fixture_digests(fx);
    test_whitebox_literals(fx);
    test_counters(fx);
    test_errors(fx);
    test_auto_increment(fx);
    test_routability(fx);
    test_outlet_policy(fx);

    std::cout << "\nPer-fixture summary (W, edge_policy whitebox):\n"
              << std::left << std::setw(12) << "fixture" << std::setw(8) << "runs"
              << std::setw(8) << "run" << std::setw(8) << "pass" << "failing run keys\n";
    for (const fixture_digest& d : k_fixture_digests)
    {
        const fixture_summary& s = g_summary[d.name];
        std::cout << std::setw(12) << d.name << std::setw(8) << s.n_runs << std::setw(8) << s.n_run
                  << std::setw(8) << s.n_pass;
        for (std::size_t i = 0; i < s.failed_keys.size(); ++i)
            std::cout << (i ? "," : "") << s.failed_keys[i];
        std::cout << "\n";
    }

    std::cout << "\n" << g_n_checks << " checks, " << g_n_failed << " failed, "
              << g_n_not_run << " not run\n"
              << "RESULT: " << (g_n_failed == 0 ? "PASS" : "FAIL") << "\n";
    return g_n_failed == 0 ? 0 : 1;
}
