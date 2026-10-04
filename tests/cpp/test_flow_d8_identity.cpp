// test_flow_d8_identity.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
//
// -----------------------------------------------------------------------------
// BYTE-IDENTITY PIN for the D8 flow object (lsdtt_flowinfo).
//
// WHY THIS EXISTS
//
// lsdtt_flowinfo::create() is being refactored for speed (interior fast
// paths in the receiver build and the s-vector donor
// search, merged passes, row-major accumulation writes, no double fills). The
// requirement is that NOT ONE BYTE of any output changes. This test was written,
// and its pins generated, from the code BEFORE that refactor (commit 7e81a39),
// so every later change is measured against the original behaviour.
//
// The junction network built on top of each fixture is audited the same way:
// lsdtt_junction_network_cpu::validate() no longer runs inside its create()
// either, so it is called here for every fixture, at a dense, a moderate and a
// sparse channel threshold, and for each junction-network option combination.
//
// It is also the replacement for lsdtt_flowinfo::validate() running on every
// build: code correctness is established here, once, in CTest, rather than paid
// for in every flow routing call. Every object built below is passed through
// validate() explicitly -- every fixture, every configuration, both thread
// counts.
//
// WHAT IS PINNED
//
// For every fixture, every option configuration, at 1 and 3 threads, a 64-bit
// FNV-1a hash of each of:
//
//   receiver   receiver_of_unchecked() for every cell (valid or not)
//   direction  receiver_dir_of_unchecked() for every cell
//   donors     donor_count_of_unchecked() for every cell
//   outlets    outlets()
//   topo       topological_order()          (configs that keep it)
//   svector    the full s-vector, rebuilt from upslope_nodes() of each outlet,
//              plus s_vector_start / s_vector_size for every cell
//   area       get_drainage_area_raster() data AND mask
//   weighted   accumulate() with non-uniform weights (the cosmogenic path)
//
// Components are pinned separately so a failure names what changed. Absent
// components (e.g. no s-vector in a config) hash to a fixed sentinel.
//
// WHY THE FIXTURES ARE INTEGER TERRAIN
//
// tests/data/dem_cache is gitignored, so real DEMs are not available on the
// server or on a fresh clone; and surfaces built with libm functions or filled
// with a float increment can differ in the last bit between libm versions or
// where a compiler contracts to FMA. Every elevation here is an integer:
//
//     z = 1000 * d + noise,   noise in {0, 250, 500, 750}
//
// where d is the BFS step distance (through the same boundary-condition rules
// the router uses) to the nearest cell that may drain out. A cell at distance d
// always has a neighbour at d - 1 at least 250 lower, so every surface is
// routable by construction, needs no fill, and has MANY exactly equal
// neighbour slopes -- which is what exercises the strict `>` first-wins
// tie-break. The only floating point in the router is division by the cell
// size and addition, both exactly reproducible under IEEE-754.
//
// REGENERATING PINS
//
// Only for a DELIBERATE behaviour change, recorded with its reason:
//     lsdtt-test-flow-d8-identity --print-pins
// and paste the output over the k_pins table. A speed refactor must never
// need this.
// -----------------------------------------------------------------------------

#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <deque>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "lsdtt_flowinfo_cpu.hpp"
#include "lsdtt_junction_network_cpu.hpp"

namespace
{

// =============================================================================
// Check harness (same shape as the other offline test drivers)
// =============================================================================
int g_n_checks = 0;
int g_n_failed = 0;

void check(const std::string& name, bool ok, const std::string& detail = "")
{
    ++g_n_checks;
    if (!ok) ++g_n_failed;
    std::cout << (ok ? "PASS" : "FAIL") << "  " << name;
    if (!detail.empty()) std::cout << "   [" << detail << "]";
    std::cout << "\n";
}

// =============================================================================
// Hashing
// =============================================================================
struct fnv64
{
    std::uint64_t h = 1469598103934665603ULL;
    void bytes(const void* p, std::size_t n)
    {
        const auto* c = static_cast<const unsigned char*>(p);
        for (std::size_t i = 0; i < n; ++i) { h ^= c[i]; h *= 1099511628211ULL; }
    }
    template <class T> void value(const T& v) { bytes(&v, sizeof(T)); }
    template <class T> void vec(const std::vector<T>& v)
    {
        value(static_cast<std::uint64_t>(v.size()));
        if (!v.empty()) bytes(v.data(), v.size() * sizeof(T));
    }
};

constexpr std::uint64_t k_absent = 0xabad1deaabad1deaULL;

enum component : int
{
    c_receiver = 0, c_direction, c_donors, c_outlets, c_topo, c_svector,
    c_area, c_weighted, c_count
};
const char* const k_component_names[c_count] = {
    "receiver", "direction", "donors", "outlets", "topo", "svector",
    "area", "weighted"
};

// =============================================================================
// Fixtures
// =============================================================================
constexpr std::array<int, 8> d8_row{{-1, -1, 0, 1, 1, 1, 0, -1}};
constexpr std::array<int, 8> d8_col{{0, 1, 1, 1, 0, -1, -1, -1}};

enum class mask_kind { none, holes, border };

struct fixture_spec
{
    std::string name;
    int         n_rows;
    int         n_cols;
    double      cell_size;
    std::string bc;                        // N E S W
    mask_kind   mask;
    bool        interior_nodata_is_outlet;
    int         n_forced_outlets;          // placed deterministically
    bool        pits_as_outlets;           // pure noise surface, no d term
    std::uint64_t seed;
};

struct fixture
{
    fixture_spec              spec;
    lsdtt::lsdtt_raster_meta  meta;
    lsdtt::lsdtt_array2d_f    dem;
    lsdtt::lsdtt_array2d_u8   mask;
    std::vector<std::int64_t> forced;
};

std::uint64_t lcg(std::uint64_t& s)
{
    s = s * 6364136223846793005ULL + 1442695040888963407ULL;
    return s >> 33;
}

// The router's own neighbour rule, restated: returns false for a no-flux
// step, sets `off` for a step off an 'o'/'b' edge, else writes (nr, nc).
bool step(const fixture_spec& f, int row, int col, int dir,
          bool& off, int& nr, int& nc)
{
    int r = row + d8_row[dir], c = col + d8_col[dir];
    off = false;
    auto edge = [&](char b, int& v, int wrap) {
        if (b == 'p') { v = wrap; return true; }
        if (b == 'o' || b == 'b') { off = true; return true; }
        return false;
    };
    if (r < 0 && !edge(f.bc[0], r, f.n_rows - 1)) return false;
    if (r >= f.n_rows && !edge(f.bc[2], r, 0)) return false;
    if (c < 0 && !edge(f.bc[3], c, f.n_cols - 1)) return false;
    if (c >= f.n_cols && !edge(f.bc[1], c, 0)) return false;
    nr = r; nc = c;
    return true;
}

fixture make_fixture(const fixture_spec& spec)
{
    fixture fx;
    fx.spec = spec;
    const int R = spec.n_rows, C = spec.n_cols;
    const std::int64_t N = static_cast<std::int64_t>(R) * C;

    fx.meta.n_rows        = R;
    fx.meta.n_cols        = C;
    fx.meta.cell_size     = spec.cell_size;
    fx.meta.no_data_value = -9999.0;
    fx.meta.set_default_synthetic_georeference();

    std::uint64_t s = spec.seed;
    std::vector<std::uint8_t> m(static_cast<std::size_t>(N), 1);
    for (int r = 0; r < R; ++r)
        for (int c = 0; c < C; ++c)
        {
            const std::size_t i = static_cast<std::size_t>(r) * C + c;
            if (spec.mask == mask_kind::holes && lcg(s) % 100 < 4) m[i] = 0;
            if (spec.mask == mask_kind::border &&
                (r < 2 || c < 3 || r >= R - 1)) m[i] = 0;
        }

    // Forced outlets on valid cells, spread deterministically.
    for (int k = 0; k < spec.n_forced_outlets; ++k)
    {
        std::int64_t f = static_cast<std::int64_t>(lcg(s) % static_cast<std::uint64_t>(N));
        while (m[static_cast<std::size_t>(f)] == 0) f = (f + 1) % N;
        if (std::find(fx.forced.begin(), fx.forced.end(), f) == fx.forced.end())
            fx.forced.push_back(f);
    }

    // BFS step distance to the nearest cell that may drain out.
    std::vector<std::int64_t> dist(static_cast<std::size_t>(N), -1);
    std::deque<std::int64_t> q;
    for (int r = 0; r < R; ++r)
        for (int c = 0; c < C; ++c)
        {
            const std::int64_t i = static_cast<std::int64_t>(r) * C + c;
            if (!m[static_cast<std::size_t>(i)]) continue;
            bool sink = std::find(fx.forced.begin(), fx.forced.end(), i) != fx.forced.end();
            if ((r == 0 && spec.bc[0] == 'b') || (c == C - 1 && spec.bc[1] == 'b') ||
                (r == R - 1 && spec.bc[2] == 'b') || (c == 0 && spec.bc[3] == 'b'))
                sink = true;
            for (int d = 0; d < 8 && !sink; ++d)
            {
                bool off; int nr = 0, nc = 0;
                if (!step(spec, r, c, d, off, nr, nc)) continue;
                if (off) { sink = true; break; }
                if (!m[static_cast<std::size_t>(nr) * C + nc] &&
                    spec.interior_nodata_is_outlet) sink = true;
            }
            if (sink) { dist[static_cast<std::size_t>(i)] = 0; q.push_back(i); }
        }
    while (!q.empty())
    {
        const std::int64_t i = q.front(); q.pop_front();
        const int r = static_cast<int>(i / C), c = static_cast<int>(i % C);
        for (int d = 0; d < 8; ++d)
        {
            bool off; int nr = 0, nc = 0;
            if (!step(spec, r, c, d, off, nr, nc) || off) continue;
            const std::size_t j = static_cast<std::size_t>(nr) * C + nc;
            if (!m[j] || dist[j] >= 0) continue;
            dist[j] = dist[static_cast<std::size_t>(i)] + 1;
            q.push_back(static_cast<std::int64_t>(j));
        }
    }

    fx.dem  = lsdtt::lsdtt_array2d_f(R, C, static_cast<float>(fx.meta.no_data_value));
    fx.mask = lsdtt::lsdtt_array2d_u8(R, C, std::uint8_t{0});
    for (int r = 0; r < R; ++r)
        for (int c = 0; c < C; ++c)
        {
            const std::size_t i = static_cast<std::size_t>(r) * C + c;
            if (!m[i]) continue;
            const std::uint64_t noise = (lcg(s) % 4) * 250;
            std::int64_t z;
            if (spec.pits_as_outlets)
                z = 1000 + static_cast<std::int64_t>(noise);
            else
            {
                if (dist[i] < 0)
                    throw std::logic_error("fixture " + spec.name +
                        " has a cell with no route out; fix the spec");
                z = 1000 * dist[i] + static_cast<std::int64_t>(noise);
            }
            fx.dem(r, c)  = static_cast<float>(z);
            fx.mask(r, c) = 1;
        }
    return fx;
}

const std::vector<fixture_spec>& fixture_specs()
{
    static const std::vector<fixture_spec> specs = {
        //  name                 rows cols cell     bc      mask            nd_out forced pits  seed
        {"open_square",           61,  67, 10.0,    "oooo", mask_kind::none,   true,  0, false, 11},
        {"open_nonint_cell",      47,  53, 30.137,  "oooo", mask_kind::none,   true,  0, false, 12},
        {"holes_nodata_outlet",   59,  71, 10.0,    "oooo", mask_kind::holes,  true,  0, false, 13},
        {"holes_nodata_blocked",  59,  71, 10.0,    "oooo", mask_kind::holes,  false, 0, false, 14},
        {"nodata_border",         43,  57, 10.0,    "oooo", mask_kind::border, true,  0, false, 15},
        {"noflux_forced",         41,  45, 10.0,    "nnnn", mask_kind::none,   true,  6, false, 16},
        {"periodic_forced",       37,  39, 10.0,    "pppp", mask_kind::none,   true,  4, false, 17},
        {"periodic_ns_open_ew",   35,  41, 10.0,    "popo", mask_kind::none,   true,  0, false, 18},
        {"baselevel_edges",       33,  36, 10.0,    "bbbb", mask_kind::none,   true,  0, false, 19},
        {"mixed_onbn",            38,  44, 10.0,    "onbn", mask_kind::holes,  true,  0, false, 20},
        {"forced_on_open",        40,  40, 10.0,    "oooo", mask_kind::none,   true,  5, false, 21},
        {"pits_as_outlets",       31,  37, 10.0,    "nnnn", mask_kind::none,   true,  0, true,  22},
        // Periodic wrap on 2 rows maps N and S to the same cell: the case the
        // s-vector donor de-duplication exists for.
        {"periodic_2_rows",        2,   9, 10.0,    "pppp", mask_kind::none,   true,  1, false, 23},
        // 1 column periodic: E and W wrap onto the cell itself (self-receiver skip).
        {"periodic_1_col",         7,   1, 10.0,    "pppp", mask_kind::none,   true,  1, false, 24},
        {"tiny_3x3_periodic",      3,   3, 10.0,    "pppp", mask_kind::none,   true,  1, false, 25},
        {"single_cell",            1,   1, 10.0,    "oooo", mask_kind::none,   true,  0, false, 26},
        {"two_by_n",               2,  11, 10.0,    "oooo", mask_kind::none,   true,  0, false, 27},
    };
    return specs;
}

// =============================================================================
// Option configurations
// =============================================================================
struct config
{
    const char* name;
    bool topo;
    bool svector;
    bool acyclic;
};
constexpr config k_configs[] = {
    {"sv_on",          true,  true,  true},   // what the CLIs build
    {"sv_off",         true,  false, true},   // build_flow_s_vector: false
    {"sv_only",        false, true,  true},   // order built internally, released
    {"receivers_only", false, false, false},  // no order, no s-vector
};

using hash_row = std::array<std::uint64_t, c_count>;

hash_row hash_object(const fixture& fx, const config& cfg, int threads,
                     std::vector<float>* area_data_out,
                     std::vector<std::uint8_t>* area_mask_out,
                     std::string* validate_error)
{
    lsdtt::lsdtt_flowinfo_options opts;
    opts.boundary_conditions       = fx.spec.bc;
    opts.interior_nodata_is_outlet = fx.spec.interior_nodata_is_outlet;
    opts.require_acyclic_graph     = cfg.acyclic;
    opts.build_topological_order   = cfg.topo;
    opts.build_s_vector            = cfg.svector;
    opts.allow_pits_as_outlets     = fx.spec.pits_as_outlets;
    opts.forced_outlet_flats       = fx.forced;
    opts.n_threads                 = threads;

    const lsdtt::lsdtt_flowinfo fi(fx.dem, fx.mask, fx.meta, opts);

    // The internal-consistency audit that create() used to run on every build.
    if (validate_error)
    {
        try { fi.validate(); }
        catch (const std::exception& e) { *validate_error = e.what(); }
    }

    const std::int64_t n = static_cast<std::int64_t>(fx.meta.n_rows) * fx.meta.n_cols;

    hash_row out{};
    out.fill(k_absent);

    fnv64 hr, hd, hc;
    for (std::int64_t f = 0; f < n; ++f)
    {
        hr.value(fi.receiver_of_unchecked(f));
        hd.value(fi.receiver_dir_of_unchecked(f));
        hc.value(fi.donor_count_of_unchecked(f));
    }
    out[c_receiver]  = hr.h;
    out[c_direction] = hd.h;
    out[c_donors]    = hc.h;

    fnv64 ho; ho.vec(fi.outlets()); out[c_outlets] = ho.h;

    if (fi.has_topological_order())
    {
        fnv64 ht; ht.vec(fi.topological_order()); out[c_topo] = ht.h;
    }

    if (fi.has_s_vector())
    {
        fnv64 hs;
        std::uint64_t total = 0;
        for (const std::int64_t o : fi.outlets())
        {
            const auto slice = fi.upslope_nodes(o);
            total += slice.size();
            hs.bytes(slice.data(), slice.size() * sizeof(std::int64_t));
        }
        hs.value(total);
        for (std::int64_t f = 0; f < n; ++f)
        {
            hs.value(fi.s_vector_start_unchecked(f));
            hs.value(fi.s_vector_size_unchecked(f));
        }
        out[c_svector] = hs.h;
    }

    if (fi.has_s_vector() || fi.has_topological_order())
    {
        const lsdtt::lsdtt_raster_result area = fi.get_drainage_area_raster();
        fnv64 ha; ha.vec(area.data); ha.vec(area.mask); out[c_area] = ha.h;
        if (area_data_out) *area_data_out = area.data;
        if (area_mask_out) *area_mask_out = area.mask;
    }

    if (fi.has_topological_order())
    {
        // Non-uniform, exactly representable weights: the path the cosmogenic
        // accumulations take. Order-of-addition changes would show here.
        std::vector<double> w(static_cast<std::size_t>(n));
        for (std::int64_t f = 0; f < n; ++f)
            w[static_cast<std::size_t>(f)] = 0.125 + 0.375 * static_cast<double>(f % 7);
        fnv64 hw; hw.vec(fi.accumulate(w)); out[c_weighted] = hw.h;
    }
    return out;
}

// Build junction networks on the fixture and audit each with validate().
// Sources need the flow object's s-vector and topological order.
void check_junction_networks_validate(const fixture& fx)
{
    lsdtt::lsdtt_flowinfo_options opts;
    opts.boundary_conditions       = fx.spec.bc;
    opts.interior_nodata_is_outlet = fx.spec.interior_nodata_is_outlet;
    opts.allow_pits_as_outlets     = fx.spec.pits_as_outlets;
    opts.forced_outlet_flats       = fx.forced;
    const lsdtt::lsdtt_flowinfo flow(fx.dem, fx.mask, fx.meta, opts);

    struct jcfg { const char* name; bool topo; bool sv; };
    constexpr jcfg jcfgs[] = {
        {"topo+sv", true, true}, {"topo_only", true, false}, {"neither", false, false}};

    const int sparse = static_cast<int>(std::max<std::int64_t>(
        1, (static_cast<std::int64_t>(fx.spec.n_rows) * fx.spec.n_cols) / 8));

    int n_built = 0, n_failed = 0;
    std::string first_error;
    for (const int threshold : {3, 12, sparse})
    {
        const std::vector<std::int64_t> sources =
            flow.get_sources_index_threshold(threshold);
        // An empty source list is an input error by design (create() rejects
        // it); a tiny grid at a sparse threshold simply has no channel.
        if (sources.empty()) continue;
        for (const jcfg& jc : jcfgs)
        {
            lsdtt::lsdtt_junction_network_options jo;
            jo.build_topological_order = jc.topo;
            jo.build_s_vector          = jc.sv;
            std::string err;
            try
            {
                const lsdtt::lsdtt_junction_network_cpu jn(flow, sources, jo);
                jn.validate();
                ++n_built;
            }
            catch (const std::exception& e) { err = e.what(); }
            if (!err.empty())
            {
                ++n_failed;
                if (first_error.empty())
                    first_error = "threshold " + std::to_string(threshold) + " "
                                + jc.name + ": " + err;
            }
        }
    }
    std::ostringstream dt;
    dt << n_built << " networks built and validated";
    if (!first_error.empty()) dt << "; first failure: " << first_error;
    check(fx.spec.name + " junction networks validate()", n_failed == 0, dt.str());
}

// How many valid non-outlet cells have two or more neighbours on the exact
// steepest slope -- i.e. how much the tie-break is actually exercised.
std::int64_t count_tied_cells(const fixture& fx)
{
    std::int64_t tied = 0;
    const int R = fx.spec.n_rows, C = fx.spec.n_cols;
    const double diag = fx.meta.cell_size * 1.4142135623730950488;
    for (int r = 0; r < R; ++r)
        for (int c = 0; c < C; ++c)
        {
            if (!fx.mask(r, c)) continue;
            double best = 0.0; int n_best = 0;
            for (int d = 0; d < 8; ++d)
            {
                bool off; int nr = 0, nc = 0;
                if (!step(fx.spec, r, c, d, off, nr, nc) || off) continue;
                if (!fx.mask(nr, nc)) continue;
                const double dz = static_cast<double>(fx.dem(r, c)) - fx.dem(nr, nc);
                if (dz <= 0.0) continue;
                const double sl = dz / ((d % 2 == 0) ? fx.meta.cell_size : diag);
                if (sl > best) { best = sl; n_best = 1; }
                else if (sl == best) ++n_best;
            }
            if (n_best >= 2) ++tied;
        }
    return tied;
}

// =============================================================================
// PINS -- generated with --print-pins from commit 7e81a39 (pre-refactor).
// Order per row: receiver, direction, donors, outlets, topo, svector, area,
// weighted. 0xabad1deaabad1dea = component absent in that configuration.
// =============================================================================
const std::map<std::string, hash_row>& pins()
{
    static const std::map<std::string, hash_row> k_pins = {
        {"baselevel_edges/receivers_only", {0x29dfb643928f825fULL, 0x6406e59b9f99b25bULL, 0x90b8a6a17476aa85ULL, 0x6cc87171295ec2daULL, 0xabad1deaabad1deaULL, 0xabad1deaabad1deaULL, 0xabad1deaabad1deaULL, 0xabad1deaabad1deaULL}},
        {"baselevel_edges/sv_off", {0x29dfb643928f825fULL, 0x6406e59b9f99b25bULL, 0x90b8a6a17476aa85ULL, 0x6cc87171295ec2daULL, 0x7cd72ec96d2a144bULL, 0xabad1deaabad1deaULL, 0x415a7c7a90353160ULL, 0xd4d69249449849ffULL}},
        {"baselevel_edges/sv_on", {0x29dfb643928f825fULL, 0x6406e59b9f99b25bULL, 0x90b8a6a17476aa85ULL, 0x6cc87171295ec2daULL, 0x7cd72ec96d2a144bULL, 0x7ddd82c83276005bULL, 0x415a7c7a90353160ULL, 0xd4d69249449849ffULL}},
        {"baselevel_edges/sv_only", {0x29dfb643928f825fULL, 0x6406e59b9f99b25bULL, 0x90b8a6a17476aa85ULL, 0x6cc87171295ec2daULL, 0xabad1deaabad1deaULL, 0x7ddd82c83276005bULL, 0x415a7c7a90353160ULL, 0xabad1deaabad1deaULL}},
        {"forced_on_open/receivers_only", {0xcd76680751b19e20ULL, 0x91ec9ac38392acb6ULL, 0x5057ad887e1ffd14ULL, 0x13c5219de7486442ULL, 0xabad1deaabad1deaULL, 0xabad1deaabad1deaULL, 0xabad1deaabad1deaULL, 0xabad1deaabad1deaULL}},
        {"forced_on_open/sv_off", {0xcd76680751b19e20ULL, 0x91ec9ac38392acb6ULL, 0x5057ad887e1ffd14ULL, 0x13c5219de7486442ULL, 0x27973022aa813c2dULL, 0xabad1deaabad1deaULL, 0xee0261880c19781bULL, 0x11e19c7ac1376925ULL}},
        {"forced_on_open/sv_on", {0xcd76680751b19e20ULL, 0x91ec9ac38392acb6ULL, 0x5057ad887e1ffd14ULL, 0x13c5219de7486442ULL, 0x27973022aa813c2dULL, 0x20cdd293519e5845ULL, 0xee0261880c19781bULL, 0x11e19c7ac1376925ULL}},
        {"forced_on_open/sv_only", {0xcd76680751b19e20ULL, 0x91ec9ac38392acb6ULL, 0x5057ad887e1ffd14ULL, 0x13c5219de7486442ULL, 0xabad1deaabad1deaULL, 0x20cdd293519e5845ULL, 0xee0261880c19781bULL, 0xabad1deaabad1deaULL}},
        {"holes_nodata_blocked/receivers_only", {0xb5e0fe65c870f1aeULL, 0x5f30dda972d2ebe7ULL, 0x012a20bee48c7a31ULL, 0x72c26e718f42f342ULL, 0xabad1deaabad1deaULL, 0xabad1deaabad1deaULL, 0xabad1deaabad1deaULL, 0xabad1deaabad1deaULL}},
        {"holes_nodata_blocked/sv_off", {0xb5e0fe65c870f1aeULL, 0x5f30dda972d2ebe7ULL, 0x012a20bee48c7a31ULL, 0x72c26e718f42f342ULL, 0xc029c1ca7a796678ULL, 0xabad1deaabad1deaULL, 0x1bb49aff268996a7ULL, 0x3a4d7bbbc415da9cULL}},
        {"holes_nodata_blocked/sv_on", {0xb5e0fe65c870f1aeULL, 0x5f30dda972d2ebe7ULL, 0x012a20bee48c7a31ULL, 0x72c26e718f42f342ULL, 0xc029c1ca7a796678ULL, 0x211c9fa4f1786fd4ULL, 0x1bb49aff268996a7ULL, 0x3a4d7bbbc415da9cULL}},
        {"holes_nodata_blocked/sv_only", {0xb5e0fe65c870f1aeULL, 0x5f30dda972d2ebe7ULL, 0x012a20bee48c7a31ULL, 0x72c26e718f42f342ULL, 0xabad1deaabad1deaULL, 0x211c9fa4f1786fd4ULL, 0x1bb49aff268996a7ULL, 0xabad1deaabad1deaULL}},
        {"holes_nodata_outlet/receivers_only", {0xfa9488d73bd0f028ULL, 0x6d5ce94e6b1753abULL, 0xce03978ea9f85507ULL, 0x5557e199a1d3e6baULL, 0xabad1deaabad1deaULL, 0xabad1deaabad1deaULL, 0xabad1deaabad1deaULL, 0xabad1deaabad1deaULL}},
        {"holes_nodata_outlet/sv_off", {0xfa9488d73bd0f028ULL, 0x6d5ce94e6b1753abULL, 0xce03978ea9f85507ULL, 0x5557e199a1d3e6baULL, 0x48d15f355a3926bdULL, 0xabad1deaabad1deaULL, 0x32f65f41b7a6566dULL, 0x2da7c306f084ba22ULL}},
        {"holes_nodata_outlet/sv_on", {0xfa9488d73bd0f028ULL, 0x6d5ce94e6b1753abULL, 0xce03978ea9f85507ULL, 0x5557e199a1d3e6baULL, 0x48d15f355a3926bdULL, 0x7c284e51fbe7cb7bULL, 0x32f65f41b7a6566dULL, 0x2da7c306f084ba22ULL}},
        {"holes_nodata_outlet/sv_only", {0xfa9488d73bd0f028ULL, 0x6d5ce94e6b1753abULL, 0xce03978ea9f85507ULL, 0x5557e199a1d3e6baULL, 0xabad1deaabad1deaULL, 0x7c284e51fbe7cb7bULL, 0x32f65f41b7a6566dULL, 0xabad1deaabad1deaULL}},
        {"mixed_onbn/receivers_only", {0x300244e60868de03ULL, 0x367fe6b3ef6f72a3ULL, 0xbef324851d890e1bULL, 0x02b1377f4f2ec1a8ULL, 0xabad1deaabad1deaULL, 0xabad1deaabad1deaULL, 0xabad1deaabad1deaULL, 0xabad1deaabad1deaULL}},
        {"mixed_onbn/sv_off", {0x300244e60868de03ULL, 0x367fe6b3ef6f72a3ULL, 0xbef324851d890e1bULL, 0x02b1377f4f2ec1a8ULL, 0x057c92c8bd7e1bdaULL, 0xabad1deaabad1deaULL, 0x82462f0c3f1ac474ULL, 0xfae5dd5af9b09da7ULL}},
        {"mixed_onbn/sv_on", {0x300244e60868de03ULL, 0x367fe6b3ef6f72a3ULL, 0xbef324851d890e1bULL, 0x02b1377f4f2ec1a8ULL, 0x057c92c8bd7e1bdaULL, 0x40f245f9703a2a9cULL, 0x82462f0c3f1ac474ULL, 0xfae5dd5af9b09da7ULL}},
        {"mixed_onbn/sv_only", {0x300244e60868de03ULL, 0x367fe6b3ef6f72a3ULL, 0xbef324851d890e1bULL, 0x02b1377f4f2ec1a8ULL, 0xabad1deaabad1deaULL, 0x40f245f9703a2a9cULL, 0x82462f0c3f1ac474ULL, 0xabad1deaabad1deaULL}},
        {"nodata_border/receivers_only", {0xc987492a69ff2dc5ULL, 0x1c036ff51837a360ULL, 0xba0fbede802130c7ULL, 0xd115368ac8fee742ULL, 0xabad1deaabad1deaULL, 0xabad1deaabad1deaULL, 0xabad1deaabad1deaULL, 0xabad1deaabad1deaULL}},
        {"nodata_border/sv_off", {0xc987492a69ff2dc5ULL, 0x1c036ff51837a360ULL, 0xba0fbede802130c7ULL, 0xd115368ac8fee742ULL, 0x61d998bfe3c14928ULL, 0xabad1deaabad1deaULL, 0xa13d36cf982d562cULL, 0xf1ea85dcad85eedfULL}},
        {"nodata_border/sv_on", {0xc987492a69ff2dc5ULL, 0x1c036ff51837a360ULL, 0xba0fbede802130c7ULL, 0xd115368ac8fee742ULL, 0x61d998bfe3c14928ULL, 0x8b31df1ee90634b8ULL, 0xa13d36cf982d562cULL, 0xf1ea85dcad85eedfULL}},
        {"nodata_border/sv_only", {0xc987492a69ff2dc5ULL, 0x1c036ff51837a360ULL, 0xba0fbede802130c7ULL, 0xd115368ac8fee742ULL, 0xabad1deaabad1deaULL, 0x8b31df1ee90634b8ULL, 0xa13d36cf982d562cULL, 0xabad1deaabad1deaULL}},
        {"noflux_forced/receivers_only", {0x1529ca3a3ccf0273ULL, 0xe18353f130318425ULL, 0xc11691161d696602ULL, 0x7d0ffe15e4555e5cULL, 0xabad1deaabad1deaULL, 0xabad1deaabad1deaULL, 0xabad1deaabad1deaULL, 0xabad1deaabad1deaULL}},
        {"noflux_forced/sv_off", {0x1529ca3a3ccf0273ULL, 0xe18353f130318425ULL, 0xc11691161d696602ULL, 0x7d0ffe15e4555e5cULL, 0xc678cf881520b9e6ULL, 0xabad1deaabad1deaULL, 0x4112d2660091305cULL, 0x1bb0082c978ac26dULL}},
        {"noflux_forced/sv_on", {0x1529ca3a3ccf0273ULL, 0xe18353f130318425ULL, 0xc11691161d696602ULL, 0x7d0ffe15e4555e5cULL, 0xc678cf881520b9e6ULL, 0x3aeae107db9bc576ULL, 0x4112d2660091305cULL, 0x1bb0082c978ac26dULL}},
        {"noflux_forced/sv_only", {0x1529ca3a3ccf0273ULL, 0xe18353f130318425ULL, 0xc11691161d696602ULL, 0x7d0ffe15e4555e5cULL, 0xabad1deaabad1deaULL, 0x3aeae107db9bc576ULL, 0x4112d2660091305cULL, 0xabad1deaabad1deaULL}},
        {"open_nonint_cell/receivers_only", {0xee3deb3821d82e49ULL, 0xb23a9c63c64bdd34ULL, 0xaf89a96d82a1a405ULL, 0x08f1ced35a42ba38ULL, 0xabad1deaabad1deaULL, 0xabad1deaabad1deaULL, 0xabad1deaabad1deaULL, 0xabad1deaabad1deaULL}},
        {"open_nonint_cell/sv_off", {0xee3deb3821d82e49ULL, 0xb23a9c63c64bdd34ULL, 0xaf89a96d82a1a405ULL, 0x08f1ced35a42ba38ULL, 0x46923fd73a4d62d3ULL, 0xabad1deaabad1deaULL, 0x3521a8df6e65a8ccULL, 0x42053156dd9e5177ULL}},
        {"open_nonint_cell/sv_on", {0xee3deb3821d82e49ULL, 0xb23a9c63c64bdd34ULL, 0xaf89a96d82a1a405ULL, 0x08f1ced35a42ba38ULL, 0x46923fd73a4d62d3ULL, 0x2bc984ca49d69791ULL, 0x3521a8df6e65a8ccULL, 0x42053156dd9e5177ULL}},
        {"open_nonint_cell/sv_only", {0xee3deb3821d82e49ULL, 0xb23a9c63c64bdd34ULL, 0xaf89a96d82a1a405ULL, 0x08f1ced35a42ba38ULL, 0xabad1deaabad1deaULL, 0x2bc984ca49d69791ULL, 0x3521a8df6e65a8ccULL, 0xabad1deaabad1deaULL}},
        {"open_square/receivers_only", {0x7aa013143f0d29dcULL, 0x55d1ac98a5bdb1f3ULL, 0xc94e4a4fba5e8d1cULL, 0x8eb6e0da03bf69f9ULL, 0xabad1deaabad1deaULL, 0xabad1deaabad1deaULL, 0xabad1deaabad1deaULL, 0xabad1deaabad1deaULL}},
        {"open_square/sv_off", {0x7aa013143f0d29dcULL, 0x55d1ac98a5bdb1f3ULL, 0xc94e4a4fba5e8d1cULL, 0x8eb6e0da03bf69f9ULL, 0xe94fee55dd57532fULL, 0xabad1deaabad1deaULL, 0x93dcea7be82928b1ULL, 0xf44a951c8977ba30ULL}},
        {"open_square/sv_on", {0x7aa013143f0d29dcULL, 0x55d1ac98a5bdb1f3ULL, 0xc94e4a4fba5e8d1cULL, 0x8eb6e0da03bf69f9ULL, 0xe94fee55dd57532fULL, 0x313f1ef91204efa7ULL, 0x93dcea7be82928b1ULL, 0xf44a951c8977ba30ULL}},
        {"open_square/sv_only", {0x7aa013143f0d29dcULL, 0x55d1ac98a5bdb1f3ULL, 0xc94e4a4fba5e8d1cULL, 0x8eb6e0da03bf69f9ULL, 0xabad1deaabad1deaULL, 0x313f1ef91204efa7ULL, 0x93dcea7be82928b1ULL, 0xabad1deaabad1deaULL}},
        {"periodic_1_col/receivers_only", {0xb6734340959c3698ULL, 0xc28801ba30c9ef4aULL, 0xb61e0ad4ea398697ULL, 0x35254a4ef7ecc9e1ULL, 0xabad1deaabad1deaULL, 0xabad1deaabad1deaULL, 0xabad1deaabad1deaULL, 0xabad1deaabad1deaULL}},
        {"periodic_1_col/sv_off", {0xb6734340959c3698ULL, 0xc28801ba30c9ef4aULL, 0xb61e0ad4ea398697ULL, 0x35254a4ef7ecc9e1ULL, 0xedb5a868aaff9ac3ULL, 0xabad1deaabad1deaULL, 0x247d2e53eea1a011ULL, 0xdc5b3f17dd499bb6ULL}},
        {"periodic_1_col/sv_on", {0xb6734340959c3698ULL, 0xc28801ba30c9ef4aULL, 0xb61e0ad4ea398697ULL, 0x35254a4ef7ecc9e1ULL, 0xedb5a868aaff9ac3ULL, 0x5d6c668d8faf04a3ULL, 0x247d2e53eea1a011ULL, 0xdc5b3f17dd499bb6ULL}},
        {"periodic_1_col/sv_only", {0xb6734340959c3698ULL, 0xc28801ba30c9ef4aULL, 0xb61e0ad4ea398697ULL, 0x35254a4ef7ecc9e1ULL, 0xabad1deaabad1deaULL, 0x5d6c668d8faf04a3ULL, 0x247d2e53eea1a011ULL, 0xabad1deaabad1deaULL}},
        {"periodic_2_rows/receivers_only", {0x029803136202be42ULL, 0xf35a33bdf93c2072ULL, 0x4e48e167dbc5cd82ULL, 0x89ebd7b27034f94cULL, 0xabad1deaabad1deaULL, 0xabad1deaabad1deaULL, 0xabad1deaabad1deaULL, 0xabad1deaabad1deaULL}},
        {"periodic_2_rows/sv_off", {0x029803136202be42ULL, 0xf35a33bdf93c2072ULL, 0x4e48e167dbc5cd82ULL, 0x89ebd7b27034f94cULL, 0x6537cfcec931f2d0ULL, 0xabad1deaabad1deaULL, 0xf0ef6cc975e9a6c1ULL, 0x155986f5f5740980ULL}},
        {"periodic_2_rows/sv_on", {0x029803136202be42ULL, 0xf35a33bdf93c2072ULL, 0x4e48e167dbc5cd82ULL, 0x89ebd7b27034f94cULL, 0x6537cfcec931f2d0ULL, 0xe36bd49e08cd6460ULL, 0xf0ef6cc975e9a6c1ULL, 0x155986f5f5740980ULL}},
        {"periodic_2_rows/sv_only", {0x029803136202be42ULL, 0xf35a33bdf93c2072ULL, 0x4e48e167dbc5cd82ULL, 0x89ebd7b27034f94cULL, 0xabad1deaabad1deaULL, 0xe36bd49e08cd6460ULL, 0xf0ef6cc975e9a6c1ULL, 0xabad1deaabad1deaULL}},
        {"periodic_forced/receivers_only", {0x118ef2ac377faceaULL, 0x32bf29504342b4b8ULL, 0x1b382a2d448454dcULL, 0x62c3b423ed30921fULL, 0xabad1deaabad1deaULL, 0xabad1deaabad1deaULL, 0xabad1deaabad1deaULL, 0xabad1deaabad1deaULL}},
        {"periodic_forced/sv_off", {0x118ef2ac377faceaULL, 0x32bf29504342b4b8ULL, 0x1b382a2d448454dcULL, 0x62c3b423ed30921fULL, 0xfaf5a910b027c2c3ULL, 0xabad1deaabad1deaULL, 0x923df63a9a3d6adfULL, 0xcfba2ffd42ced756ULL}},
        {"periodic_forced/sv_on", {0x118ef2ac377faceaULL, 0x32bf29504342b4b8ULL, 0x1b382a2d448454dcULL, 0x62c3b423ed30921fULL, 0xfaf5a910b027c2c3ULL, 0x63f96d93d07257ddULL, 0x923df63a9a3d6adfULL, 0xcfba2ffd42ced756ULL}},
        {"periodic_forced/sv_only", {0x118ef2ac377faceaULL, 0x32bf29504342b4b8ULL, 0x1b382a2d448454dcULL, 0x62c3b423ed30921fULL, 0xabad1deaabad1deaULL, 0x63f96d93d07257ddULL, 0x923df63a9a3d6adfULL, 0xabad1deaabad1deaULL}},
        {"periodic_ns_open_ew/receivers_only", {0x15389e882ffdf6aeULL, 0xd0edfcf033d0a2b0ULL, 0x3f3b8931f99e1bd5ULL, 0x21e674b29e578d74ULL, 0xabad1deaabad1deaULL, 0xabad1deaabad1deaULL, 0xabad1deaabad1deaULL, 0xabad1deaabad1deaULL}},
        {"periodic_ns_open_ew/sv_off", {0x15389e882ffdf6aeULL, 0xd0edfcf033d0a2b0ULL, 0x3f3b8931f99e1bd5ULL, 0x21e674b29e578d74ULL, 0x2e2115950cfa3f9bULL, 0xabad1deaabad1deaULL, 0xb9ec23431dfe9a11ULL, 0xdff28222f4e4bae8ULL}},
        {"periodic_ns_open_ew/sv_on", {0x15389e882ffdf6aeULL, 0xd0edfcf033d0a2b0ULL, 0x3f3b8931f99e1bd5ULL, 0x21e674b29e578d74ULL, 0x2e2115950cfa3f9bULL, 0x0732439d7d68301cULL, 0xb9ec23431dfe9a11ULL, 0xdff28222f4e4bae8ULL}},
        {"periodic_ns_open_ew/sv_only", {0x15389e882ffdf6aeULL, 0xd0edfcf033d0a2b0ULL, 0x3f3b8931f99e1bd5ULL, 0x21e674b29e578d74ULL, 0xabad1deaabad1deaULL, 0x0732439d7d68301cULL, 0xb9ec23431dfe9a11ULL, 0xabad1deaabad1deaULL}},
        {"pits_as_outlets/receivers_only", {0x8d8a525827653025ULL, 0xd1579aaab2853a35ULL, 0xe254d251bed81e97ULL, 0x26dfdf29199019f7ULL, 0xabad1deaabad1deaULL, 0xabad1deaabad1deaULL, 0xabad1deaabad1deaULL, 0xabad1deaabad1deaULL}},
        {"pits_as_outlets/sv_off", {0x8d8a525827653025ULL, 0xd1579aaab2853a35ULL, 0xe254d251bed81e97ULL, 0x26dfdf29199019f7ULL, 0xdd0271fbe652e08fULL, 0xabad1deaabad1deaULL, 0xd08747f48af86ac1ULL, 0xd626b05d09725316ULL}},
        {"pits_as_outlets/sv_on", {0x8d8a525827653025ULL, 0xd1579aaab2853a35ULL, 0xe254d251bed81e97ULL, 0x26dfdf29199019f7ULL, 0xdd0271fbe652e08fULL, 0x7dce996a15e3d8fcULL, 0xd08747f48af86ac1ULL, 0xd626b05d09725316ULL}},
        {"pits_as_outlets/sv_only", {0x8d8a525827653025ULL, 0xd1579aaab2853a35ULL, 0xe254d251bed81e97ULL, 0x26dfdf29199019f7ULL, 0xabad1deaabad1deaULL, 0x7dce996a15e3d8fcULL, 0xd08747f48af86ac1ULL, 0xabad1deaabad1deaULL}},
        {"single_cell/receivers_only", {0xbf141897110b0c3aULL, 0x44bd24d473ccebb4ULL, 0x44bd2bd473ccf799ULL, 0x5420115802dc1402ULL, 0xabad1deaabad1deaULL, 0xabad1deaabad1deaULL, 0xabad1deaabad1deaULL, 0xabad1deaabad1deaULL}},
        {"single_cell/sv_off", {0xbf141897110b0c3aULL, 0x44bd24d473ccebb4ULL, 0x44bd2bd473ccf799ULL, 0x5420115802dc1402ULL, 0x5420115802dc1402ULL, 0xabad1deaabad1deaULL, 0xbc04e78b04f846ecULL, 0x56ac64580506264bULL}},
        {"single_cell/sv_on", {0xbf141897110b0c3aULL, 0x44bd24d473ccebb4ULL, 0x44bd2bd473ccf799ULL, 0x5420115802dc1402ULL, 0x5420115802dc1402ULL, 0x0f6ed56b1666e4c3ULL, 0xbc04e78b04f846ecULL, 0x56ac64580506264bULL}},
        {"single_cell/sv_only", {0xbf141897110b0c3aULL, 0x44bd24d473ccebb4ULL, 0x44bd2bd473ccf799ULL, 0x5420115802dc1402ULL, 0xabad1deaabad1deaULL, 0x0f6ed56b1666e4c3ULL, 0xbc04e78b04f846ecULL, 0xabad1deaabad1deaULL}},
        {"tiny_3x3_periodic/receivers_only", {0x4a8f50c6a08998faULL, 0xe64a276c966ce18aULL, 0xe31e01444eec3e01ULL, 0xb110667323a9f265ULL, 0xabad1deaabad1deaULL, 0xabad1deaabad1deaULL, 0xabad1deaabad1deaULL, 0xabad1deaabad1deaULL}},
        {"tiny_3x3_periodic/sv_off", {0x4a8f50c6a08998faULL, 0xe64a276c966ce18aULL, 0xe31e01444eec3e01ULL, 0xb110667323a9f265ULL, 0x15433362b658c162ULL, 0xabad1deaabad1deaULL, 0xb64eb506c4633d71ULL, 0x684ba53da5aed8d9ULL}},
        {"tiny_3x3_periodic/sv_on", {0x4a8f50c6a08998faULL, 0xe64a276c966ce18aULL, 0xe31e01444eec3e01ULL, 0xb110667323a9f265ULL, 0x15433362b658c162ULL, 0x2603f0299ae8b8c3ULL, 0xb64eb506c4633d71ULL, 0x684ba53da5aed8d9ULL}},
        {"tiny_3x3_periodic/sv_only", {0x4a8f50c6a08998faULL, 0xe64a276c966ce18aULL, 0xe31e01444eec3e01ULL, 0xb110667323a9f265ULL, 0xabad1deaabad1deaULL, 0x2603f0299ae8b8c3ULL, 0xb64eb506c4633d71ULL, 0xabad1deaabad1deaULL}},
        {"two_by_n/receivers_only", {0x5fa60e40cf5e4eaaULL, 0x1520b125eb30429eULL, 0xd77e8119c3408d49ULL, 0x864d0b905b20586cULL, 0xabad1deaabad1deaULL, 0xabad1deaabad1deaULL, 0xabad1deaabad1deaULL, 0xabad1deaabad1deaULL}},
        {"two_by_n/sv_off", {0x5fa60e40cf5e4eaaULL, 0x1520b125eb30429eULL, 0xd77e8119c3408d49ULL, 0x864d0b905b20586cULL, 0x4f80482045319554ULL, 0xabad1deaabad1deaULL, 0xa2ffed654d4ff410ULL, 0xcb4e940532b4a017ULL}},
        {"two_by_n/sv_on", {0x5fa60e40cf5e4eaaULL, 0x1520b125eb30429eULL, 0xd77e8119c3408d49ULL, 0x864d0b905b20586cULL, 0x4f80482045319554ULL, 0xf5404cd6045265f1ULL, 0xa2ffed654d4ff410ULL, 0xcb4e940532b4a017ULL}},
        {"two_by_n/sv_only", {0x5fa60e40cf5e4eaaULL, 0x1520b125eb30429eULL, 0xd77e8119c3408d49ULL, 0x864d0b905b20586cULL, 0xabad1deaabad1deaULL, 0xf5404cd6045265f1ULL, 0xa2ffed654d4ff410ULL, 0xabad1deaabad1deaULL}},
    };
    return k_pins;
}

} // namespace

int main(int argc, char** argv)
{
    const bool print_pins = argc > 1 && std::string(argv[1]) == "--print-pins";

    std::cout << "lsdtt3 D8 flow object byte-identity pin\n\n";

    std::map<std::string, hash_row> computed;

    for (const fixture_spec& spec : fixture_specs())
    {
        const fixture fx = make_fixture(spec);

        if (!print_pins)
        {
            std::ostringstream dt;
            const std::int64_t tied = count_tied_cells(fx);
            dt << tied << " cells with a tied steepest slope";
            // Only the larger fixtures are required to exercise ties; tiny
            // grids cannot be expected to.
            if (spec.n_rows * spec.n_cols >= 400)
                check("fixture " + spec.name + " exercises the tie-break", tied > 0, dt.str());
        }

        std::vector<float> area_sv_on, area_sv_off;
        std::vector<std::uint8_t> mask_sv_on, mask_sv_off;

        for (const config& cfg : k_configs)
        {
            const std::string key = spec.name + "/" + cfg.name;
            hash_row h1, h3;
            std::string err, verr1, verr3;
            try
            {
                h1 = hash_object(fx, cfg, 1,
                                 std::string(cfg.name) == "sv_on"  ? &area_sv_on
                               : std::string(cfg.name) == "sv_off" ? &area_sv_off : nullptr,
                                 std::string(cfg.name) == "sv_on"  ? &mask_sv_on
                               : std::string(cfg.name) == "sv_off" ? &mask_sv_off : nullptr,
                                 &verr1);
                h3 = hash_object(fx, cfg, 3, nullptr, nullptr, &verr3);
            }
            catch (const std::exception& e)
            {
                err = e.what();
            }
            if (!err.empty())
            {
                check(key + " builds", false, err);
                continue;
            }
            computed[key] = h1;
            if (print_pins) continue;

            check(key + " identical at 1 and 3 threads", h1 == h3);
            check(key + " validate() passes at 1 and 3 threads",
                  verr1.empty() && verr3.empty(), verr1.empty() ? verr3 : verr1);

            const auto it = pins().find(key);
            if (it == pins().end())
            {
                check(key + " has a pin", false, "run --print-pins");
                continue;
            }
            for (int k = 0; k < c_count; ++k)
            {
                std::ostringstream dt;
                dt << std::hex << "pinned 0x" << it->second[k] << ", got 0x" << h1[k];
                check(key + " " + k_component_names[k], h1[k] == it->second[k],
                      h1[k] == it->second[k] ? "" : dt.str());
            }
        }

        if (!print_pins)
        {
            check_junction_networks_validate(fx);

            // lsdtt3 must not behave differently with the s-vector on or off.
            check(spec.name + " drainage area identical with s-vector on and off",
                  !area_sv_on.empty() && area_sv_on.size() == area_sv_off.size() &&
                  std::memcmp(area_sv_on.data(), area_sv_off.data(),
                              area_sv_on.size() * sizeof(float)) == 0 &&
                  mask_sv_on == mask_sv_off);
        }
    }

    if (print_pins)
    {
        for (const auto& [key, h] : computed)
        {
            std::cout << "        {\"" << key << "\", {";
            for (int k = 0; k < c_count; ++k)
                std::cout << (k ? ", " : "") << "0x" << std::hex << std::setw(16)
                          << std::setfill('0') << h[k] << "ULL" << std::dec;
            std::cout << "}},\n";
        }
        return 0;
    }

    std::cout << "\n" << g_n_checks - g_n_failed << " / " << g_n_checks
              << " checks passed\n";
    std::cout << "RESULT: " << (g_n_failed == 0 ? "PASS" : "FAIL") << "\n";
    return g_n_failed == 0 ? 0 : 1;
}
