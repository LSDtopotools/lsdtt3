// test_depression_hierarchy.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
//
// -----------------------------------------------------------------------------
// Unit tests for src/cpu/lsdtt_depression_hierarchy_cpu, our port of Richard
// Barnes's depression hierarchy (Barnes et al. 2020, doi:10.5194/esurf-8-431-2020).
//
// Phase 0 of docs/plan_docs/depression_hierarchy_fix_plan.md.
//
// WHAT THESE TESTS ARE FOR
// ------------------------
// They pin the defect described in section 1 of that plan: outlet_link
// (lsdtt_depression_hierarchy_cpu.hpp) derives from std::pair and adds two
// REFERENCE members aliasing first/second.  The implicitly generated copy
// constructor binds a copy's references to the SOURCE object's members, so when
// a stack-local outlet_link is copied into outlet_db the stored key's references
// dangle as soon as that loop iteration ends.  outlet_link_hash reads exactly
// those references and libstdc++ does not cache hash codes for this type, so
// every rehash hashes reclaimed stack memory.
//
// That has two consequences, and the SECOND is the serious one:
//
//   1. every element collapses into one bucket, so each find() in the BFS inner
//      loop becomes a linear scan -- genuine O(n^2);
//   2. lookups MISS keys that are present, so the "keep the lowest saddle
//      between depressions A and B" de-duplication silently fails and the
//      stored outlet is overwritten by a higher one instead of kept.
//
// WHY THIS WAS NOT CAUGHT BEFORE
// ------------------------------
// outlet_db is reserved at 3 * (pit_cell_count + 1).  The bug only bites once
// the map rehashes past that reserve, and every pre-existing fixture is far too
// small.  A correct-looking result on a small grid proves nothing about this
// code path, which is why test_outlet_db_survives_rehash below deliberately
// forces a rehash rather than trusting a realistic-looking DEM to do it.
//
// EXPECTED STATE
// --------------
// On unmodified HEAD these tests FAIL.  That is the point of phase 0: they are
// written before the fix so that the fix can be shown to be what makes them
// pass.  Once phase 1 lands they must pass and stay passing.
//
// Usage: lsdtt-test-depression-hierarchy
// Prints one PASS/FAIL line per check and a final "RESULT: PASS|FAIL".
// -----------------------------------------------------------------------------

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <iostream>
#include <limits>
#include <string>
#include <array>
#include <chrono>
#include <type_traits>
#include <unordered_map>
#include <vector>

#include "lsdtt_data_loader.hpp"
#include "lsdtt_depression_hierarchy_cpu.hpp"
#include "lsdtt_raster_cpu.hpp"
#include "lsdtt_raster_synthetic_surfaces_cpu.hpp"

namespace
{

// =============================================================================
// Tiny check harness (same shape as the other lsdtt3 test drivers)
// =============================================================================
int g_n_checks = 0;
int g_n_failed = 0;

void check(const std::string& name, bool ok, const std::string& detail = "")
{
    ++g_n_checks;
    if (!ok) ++g_n_failed;
    std::cout << (ok ? "PASS" : "FAIL") << "  " << name;
    if (!detail.empty()) std::cout << "   [" << detail << "]";
    std::cout << std::endl;
}

void report(const std::string& name, const std::string& detail)
{
    std::cout << "NOTE  " << name << "   [" << detail << "]" << std::endl;
}

lsdtt::lsdtt_raster_meta make_meta(int n_rows, int n_cols, double cell_size)
{
    lsdtt::lsdtt_raster_meta meta;
    meta.n_rows      = n_rows;
    meta.n_cols      = n_cols;
    meta.cell_size   = cell_size;
    meta.x_min       = 0.0;
    meta.y_min       = 0.0;
    meta.no_data_value = -9999.0;
    // A projected metric CRS: build_depression_hierarchy requires one.
    meta.crs_wkt = "PROJCS[\"WGS 84 / UTM zone 58S\","
                   "GEOGCS[\"WGS 84\",DATUM[\"WGS_1984\","
                   "SPHEROID[\"WGS 84\",6378137,298.257223563]],"
                   "PRIMEM[\"Greenwich\",0],UNIT[\"degree\",0.0174532925199433]],"
                   "PROJECTION[\"Transverse_Mercator\"],"
                   "PARAMETER[\"latitude_of_origin\",0],"
                   "PARAMETER[\"central_meridian\",165],"
                   "PARAMETER[\"scale_factor\",0.9996],"
                   "PARAMETER[\"false_easting\",500000],"
                   "PARAMETER[\"false_northing\",10000000],"
                   "UNIT[\"metre\",1]]";
    return meta;
}

// =============================================================================
// TEST 1 -- the mechanism, directly
//
// Builds the exact container the algorithm builds, populates it the exact way
// the BFS inner loop does (lsdtt_depression_hierarchy_cpu.cpp:399-412: construct
// outlet_link as a loop-local, insert, let it die), and forces several rehashes
// by reserving far less than the number inserted.
//
// This pins the MECHANISM rather than a symptom, so it stays meaningful after
// the fix: it is a standing guard against anyone reintroducing reference or
// pointer members into a hashed key.
//
// Two independent assertions:
//   (a) every inserted key must be findable          -- correctness
//   (b) the keys must be spread across buckets       -- performance
//
// (a) is the one that matters scientifically.  (b) is what makes the bug show
// up as a hang rather than a wrong number.
// =============================================================================
void test_outlet_db_survives_rehash()
{
    using db_t = std::unordered_map<lsdtt::outlet_link,
                                    lsdtt::depression_outlet,
                                    lsdtt::outlet_link_hash>;

    // Deliberately tiny reserve so the map rehashes many times over the loop
    // below.  The algorithm reserves 3 * (pit_cell_count + 1); on a real DEM
    // with many depressions the link count exceeds that and it rehashes for the
    // same reason.  Forcing it here keeps the test fast and deterministic
    // instead of depending on a large fixture.
    const std::size_t n_links = 20000;

    db_t db;
    db.reserve(16);

    for (std::size_t i = 0; i < n_links; ++i)
    {
        // Scoped exactly as in the BFS inner loop: constructed, used, destroyed.
        const lsdtt::outlet_link lnk(static_cast<lsdtt::dep_label_t>(i),
                                     static_cast<lsdtt::dep_label_t>(i + 1));
        db[lnk] = lsdtt::depression_outlet(
            static_cast<lsdtt::dep_label_t>(i),
            static_cast<lsdtt::dep_label_t>(i + 1),
            static_cast<std::int64_t>(i),
            static_cast<float>(i));
    }

    check("outlet_db retains every inserted link",
          db.size() == n_links,
          "size " + std::to_string(db.size()) + " of " + std::to_string(n_links));

    // (a) Correctness: every key must be findable.  Equality on outlet_link is
    // std::pair's, which compares first/second by VALUE and is therefore
    // unaffected by the bug -- only the hash is wrong.  So a miss here means the
    // key was filed in the wrong bucket, which is precisely how the "keep the
    // lowest saddle" de-duplication fails.
    std::size_t misses = 0;
    for (std::size_t i = 0; i < n_links; ++i)
    {
        const lsdtt::outlet_link lnk(static_cast<lsdtt::dep_label_t>(i),
                                     static_cast<lsdtt::dep_label_t>(i + 1));
        if (db.find(lnk) == db.end()) ++misses;
    }

    check("every outlet_link is findable after rehashing",
          misses == 0,
          std::to_string(misses) + " of " + std::to_string(n_links) +
              " lookups missed a key that is present");

    // (b) Performance: the keys must actually be distributed.  A correct hash
    // over 20000 distinct links gives a max bucket in the low single digits; the
    // dangling-reference hash puts everything in one bucket.  The bound is
    // deliberately loose -- this is a collapse detector, not a hash-quality
    // test.
    std::size_t max_bucket = 0;
    for (std::size_t b = 0; b < db.bucket_count(); ++b)
        max_bucket = std::max(max_bucket, db.bucket_size(b));

    const std::size_t bucket_limit = 64;

    check("outlet_link_hash distributes keys across buckets",
          max_bucket <= bucket_limit,
          "max bucket " + std::to_string(max_bucket) + " (limit " +
              std::to_string(bucket_limit) + ", buckets " +
              std::to_string(db.bucket_count()) + ")");

    report("bucket occupancy",
           "size " + std::to_string(db.size()) + ", buckets " +
               std::to_string(db.bucket_count()) + ", max bucket " +
               std::to_string(max_bucket));
}

// =============================================================================
// TEST 1b -- copy semantics, in isolation
//
// The narrowest possible statement of the defect, independent of any container:
// a COPY of an outlet_link must describe the same pair as the original once the
// original is gone.  Under the current reference members it describes whatever
// happens to be in the original's reclaimed storage.
//
// Kept separate from test 1 because it needs no rehash and no map, so it says
// unambiguously WHY test 1 fails rather than merely that it does.
// =============================================================================
void test_outlet_link_copy_is_self_contained()
{
    // A hashed key must be SELF-CONTAINED: its hash must depend only on its own
    // storage. The way to state that as a test is through type properties, not
    // by reading a copy after its source dies -- that read is undefined
    // behaviour, and an undefined-behaviour test can pass by luck (the dead
    // stack slot often still holds the old value), which would make this suite
    // nondeterministic. These checks are deterministic and compile-time-decidable.

    // Reference members suppress the implicit copy-assignment operator. This is
    // both a symptom and a smell: a key type that cannot be assigned has
    // something aliasing inside it.
    check("outlet_link is copy-assignable",
          std::is_copy_assignable<lsdtt::outlet_link>::value,
          "reference members delete the implicit operator=");

    // The strongest available statement of "no aliasing members". A key of two
    // integer labels has no business being anything else.
    check("outlet_link is trivially copyable",
          std::is_trivially_copyable<lsdtt::outlet_link>::value,
          "a hashed key must own its data outright");

    // Equal links must hash equally regardless of how they were built. This
    // holds trivially for a value type; it is the invariant the hash relies on.
    const lsdtt::outlet_link a(7u, 42u);
    const lsdtt::outlet_link b(42u, 7u);   // constructor normalises the order
    const lsdtt::outlet_link_hash hash;

    check("outlet_link normalises pair order",
          a.dep_a == b.dep_a && a.dep_b == b.dep_b,
          "a(" + std::to_string(a.dep_a) + "," + std::to_string(a.dep_b) +
              ") vs b(" + std::to_string(b.dep_a) + "," + std::to_string(b.dep_b) + ")");

    check("equal outlet_links hash equally",
          hash(a) == hash(b),
          std::to_string(hash(a)) + " vs " + std::to_string(hash(b)));
}

// =============================================================================
// TEST 2 -- the algorithm, on a surface large enough to rehash
//
// An egg-carton of many small closed hollows.  Two things matter about the
// fixture:
//
//   - it must contain enough ADJACENT depression pairs to push outlet_db past
//     its 3 * (pit_cell_count + 1) reserve, because below that the bug is
//     dormant (see the header note);
//   - every hollow's true spill elevation must be knowable independently, so
//     the assertion is against ground truth rather than against the code's own
//     output.
//
// The surface is a raised plateau with a grid of square pits punched into it.
// Every pit is a simple flat-floored box, so its spill elevation is exactly the
// plateau height: to leave the pit, water must rise to the rim.  That is an
// analytic oracle, not a golden value.
//
// A depression whose out_elev is BELOW the plateau means the de-duplication
// kept, or overwrote with, the wrong saddle.
// =============================================================================
void test_hierarchy_outlets_on_dense_pit_field()
{
    // Grid of pits: pit_stride controls how many.  These numbers are chosen to
    // produce several hundred depressions -- enough adjacency to force a rehash
    // -- while still running in well under a second.
    const int    n          = 400;
    const double cell_size  = 10.0;
    const int    pit_stride = 8;     // one pit every 8 cells in each direction
    const int    pit_half   = 2;     // pits are (2*2+1) = 5 cells across
    const float  plateau    = 100.0f;
    const float  pit_floor  = 90.0f;

    lsdtt::lsdtt_raster_meta meta = make_meta(n, n, cell_size);

    std::vector<float>   data(static_cast<std::size_t>(n) * n, plateau);
    std::vector<uint8_t> mask(static_cast<std::size_t>(n) * n, 1u);

    int n_pits = 0;
    for (int pr = pit_stride; pr < n - pit_stride; pr += pit_stride)
    {
        for (int pc = pit_stride; pc < n - pit_stride; pc += pit_stride)
        {
            ++n_pits;
            for (int dr = -pit_half; dr <= pit_half; ++dr)
                for (int dc = -pit_half; dc <= pit_half; ++dc)
                {
                    const int r = pr + dr;
                    const int c = pc + dc;
                    data[static_cast<std::size_t>(r) * n + c] = pit_floor;
                }
        }
    }

    lsdtt::lsdtt_raster_cpu dem;
    dem.set_from_vectors(std::move(data), std::move(mask), meta);

    std::vector<lsdtt::dep_label_t> labels;
    std::vector<uint8_t>            flowdirs;

    const lsdtt::lsdtt_depression_hierarchy hier =
        lsdtt::build_depression_hierarchy(dem, labels, flowdirs);

    report("dense pit field",
           std::to_string(n) + "x" + std::to_string(n) + " grid, " +
               std::to_string(n_pits) + " pits punched, " +
               std::to_string(hier.size()) + " hierarchy nodes");

    // IMPORTANT, and the reason this test does not currently reproduce the
    // defect: outlet_db is reserved at 3 * (pit_cell_count + 1), and
    // std::unordered_map does not rehash while the element count stays within a
    // reserve. On this fixture the outlet count is orders of magnitude below the
    // reserve, so no rehash occurs and the dangling references are never
    // re-read. See the stderr line "BFS done. N leaf depressions, M outlets"
    // above against the pit-cell count. Forcing outlets > 3 * pit_cells through
    // the public API is an OPEN PROBLEM -- see the plan, phase 0 results.
    report("rehash headroom",
           "outlet_db reserve is 3 * (pit cells + 1); a rehash needs the outlet "
           "count to exceed it. Compare the two figures on the stderr lines above.");

    check("dense pit field produced a hierarchy",
          hier.size() > 1u,
          std::to_string(hier.size()) + " nodes");

    // Every LEAF depression is a flat-floored box whose only escape is over the
    // plateau, so its outlet elevation must be exactly the plateau height.
    // Meta-depressions (nodes with children) merge boxes and may sit higher, so
    // only leaves are checked.
    std::size_t n_leaves      = 0;
    std::size_t n_bad_outlet  = 0;
    float       worst_out     = plateau;

    for (std::size_t d = 1; d < hier.size(); ++d)
    {
        const lsdtt::depression_node& node = hier[d];

        const bool is_leaf = (node.lchild == lsdtt::k_no_label &&
                              node.rchild == lsdtt::k_no_label);
        if (!is_leaf) continue;
        ++n_leaves;

        // An unset outlet (still +inf) is a separate condition and is not what
        // this check is about; skip those rather than counting them as low.
        if (!std::isfinite(node.out_elev)) continue;

        if (node.out_elev < plateau - 1e-3f)
        {
            ++n_bad_outlet;
            worst_out = std::min(worst_out, node.out_elev);
        }
    }

    check("leaf depressions have their true spill elevation",
          n_bad_outlet == 0,
          std::to_string(n_bad_outlet) + " of " + std::to_string(n_leaves) +
              " leaves spill below the plateau (" + std::to_string(plateau) +
              "); lowest seen " + std::to_string(worst_out));

    // Every pit floor cell must belong to some depression.  A cell left
    // unlabelled means the BFS lost it.
    std::size_t unlabelled_floor = 0;
    const std::vector<float>& elev = dem.get_data();
    for (std::size_t i = 0; i < elev.size(); ++i)
        if (elev[i] <= pit_floor + 1e-3f &&
            labels[i] == lsdtt::k_no_label) ++unlabelled_floor;

    check("every pit-floor cell carries a depression label",
          unlabelled_floor == 0,
          std::to_string(unlabelled_floor) + " unlabelled floor cells");
}

// =============================================================================
// COMPLEXITY GUARD  (phase 4 of the plan, section 4.3)
//
// The outlet_link defect made this algorithm QUADRATIC in cell count: 56.9 s on
// a 400x400 noise surface, and 2 h 33 min on the 17.4 Mpx Greenland DEM. Barnes
// reports it as near-linear, and after the fix it is (0.23 s at 400x400).
//
// A wall-clock ratio test would have caught that regression on day one, which is
// why one exists now. It is deliberately loose -- this is a "has it gone
// quadratic again" alarm, not a benchmark, and it must not go off because CI is
// busy. Quadratic is exponent 2.0; the bound is 1.6.
//
// Uniform noise is used because outlet DENSITY, not cell count, is what drives
// the cost: a flat-floored plateau of the same size yields ~2300 outlets against
// ~32000 depressions for noise, and stays fast even with the defect present.
// =============================================================================
void test_scaling_is_not_quadratic()
{
    struct point { double cells; double secs; };
    const int sizes[2] = {150, 300};

    std::vector<lsdtt::lsdtt_raster_cpu> dems;
    for (const int n : sizes)
    {
        lsdtt::lsdtt_raster_meta meta = make_meta(n, n, 10.0);
        dems.push_back(lsdtt::generate_uniform_random_noise_surface(meta, 1, 0.0, 10.0));
    }

    // Best of five, INTERLEAVED (small, large, small, large, ...). Each build
    // takes ~20-200 ms, so under `ctest -j` one burst of load from a
    // neighbouring test moved the exponent past the bound: 1.66 in 2 of 3 `-j4`
    // runs on 2026-09-18, passing alone. A best-of-3 per size, timed back to
    // back, still failed once in 7 runs, because a burst can cover all three
    // runs of one size. Interleaving exposes both sizes to the same load, and
    // the minimum is the estimate closest to the true cost, since contention
    // only ever adds time.
    double best[2] = {0.0, 0.0};
    for (int rep = 0; rep < 5; ++rep)
    {
        for (int k = 0; k < 2; ++k)
        {
            std::vector<lsdtt::dep_label_t> labels;
            std::vector<uint8_t>            flowdirs;

            const auto t0 = std::chrono::steady_clock::now();
            const lsdtt::lsdtt_depression_hierarchy hier =
                lsdtt::build_depression_hierarchy(dems[k], labels, flowdirs);
            const auto t1 = std::chrono::steady_clock::now();

            (void)hier;
            const double secs = std::chrono::duration<double>(t1 - t0).count();
            if (rep == 0 || secs < best[k]) best[k] = secs;
        }
    }

    std::vector<point> pts;
    for (int k = 0; k < 2; ++k)
        pts.push_back({static_cast<double>(sizes[k]) * sizes[k], best[k]});

    // Guard against a machine so fast the smaller case is unmeasurable, which
    // would make the ratio meaningless rather than failing honestly.
    const double floor_s = 1e-4;
    if (pts[0].secs < floor_s)
    {
        report("scaling exponent",
               "smaller case too fast to time reliably (" +
                   std::to_string(pts[0].secs) + " s); check skipped");
        return;
    }

    const double exponent = std::log(pts[1].secs / pts[0].secs) /
                            std::log(pts[1].cells / pts[0].cells);

    check("depression hierarchy is not quadratic in cell count",
          exponent < 1.6,
          "exponent " + std::to_string(exponent) + " (quadratic is 2.0); " +
              std::to_string(pts[0].secs) + " s at " +
              std::to_string(static_cast<long>(pts[0].cells)) + " cells vs " +
              std::to_string(pts[1].secs) + " s at " +
              std::to_string(static_cast<long>(pts[1].cells)) + " cells");
}

// =============================================================================
// SCALING HARNESS  (--scaling)
//
// A fast stand-in for the Greenland case, so a candidate fix can be judged in
// seconds instead of the ~3 hours the 17.4 Mpx DEM takes.
//
// WHAT MAKES GREENLAND SLOW IS OUTLET DENSITY, NOT SIZE.  The plateau fixture
// above is 160k cells with 2304 outlets (1.4% of cells).  Greenland is 17.4M
// cells with >1.26M outlets (7.2%).  Uniform noise reproduces that regime: every
// local minimum becomes a depression and every depression boundary a set of
// outlet links, which is what a noisy supraglacial ice surface looks like to
// this algorithm.
//
// Reports wall time and, crucially, TIME PER CELL.  If the algorithm is
// near-linear as Barnes reports, time per cell is flat as n grows.  If it grows,
// the run is superlinear and the growth rate localises the problem.
//
// Not a CTest case: it prints numbers to be read, and its runtime is
// deliberately size-dependent.  Run it by hand when working on performance.
// =============================================================================
void run_scaling_harness()
{
    std::cout << "\n=== scaling harness: uniform noise ===\n";
    std::cout << "Doubling the cell count should roughly DOUBLE the time if the\n"
                 "algorithm is near-linear, as Barnes et al. report.\n\n";
    std::cout << "    n      cells      nodes    wall_s    us/cell   ratio_vs_prev\n";

    struct row { double cells; double secs; };
    std::vector<row> rows;

    for (const int n : {100, 150, 200, 250, 300, 400})
    {
        lsdtt::lsdtt_raster_meta meta = make_meta(n, n, 10.0);

        lsdtt::lsdtt_raster_cpu dem =
            lsdtt::generate_uniform_random_noise_surface(meta, 1, 0.0, 10.0);

        std::vector<lsdtt::dep_label_t> labels;
        std::vector<uint8_t>            flowdirs;

        const auto t0 = std::chrono::steady_clock::now();
        const lsdtt::lsdtt_depression_hierarchy hier =
            lsdtt::build_depression_hierarchy(dem, labels, flowdirs);
        const auto t1 = std::chrono::steady_clock::now();

        const double secs  = std::chrono::duration<double>(t1 - t0).count();
        const double cells = static_cast<double>(n) * n;

        double ratio = 0.0;
        if (!rows.empty())
        {
            const double cell_ratio = cells / rows.back().cells;
            const double time_ratio = secs  / rows.back().secs;
            // Exponent e such that time ~ cells^e.
            ratio = std::log(time_ratio) / std::log(cell_ratio);
        }
        rows.push_back({cells, secs});

        std::printf("%5d %10.0f %10zu %9.3f %10.2f %14s\n",
                    n, cells, hier.size(), secs, 1e6 * secs / cells,
                    rows.size() == 1 ? "-"
                                     : (std::to_string(ratio).substr(0, 5) +
                                        " (exponent)").c_str());
        std::fflush(stdout);
    }

    std::cout << "\nThe exponent column is e in  time ~ cells^e.\n"
                 "  e ~ 1.0  near-linear, which is what this algorithm should be\n"
                 "  e ~ 2.0  quadratic\n";
}

// =============================================================================
// ORACLE DUMP  (--dump <dem> [--crop])
//
// Phase 2 of the plan: prints our hierarchy in a form comparable with upstream
// Barnes, so the two can be checked against each other on the 25 test DEMs.
//
// --crop strips a fully-nodata outer border.  This exists because of the
// ocean-labelling divergence recorded in plan 3.2: upstream treats NODATA cells
// as OCEAN, while lsdtt3 seeds ocean from valid cells touching an open GRID
// EDGE and ignores nodata entirely.  Every upstream fixture has a nodata ring,
// so without cropping lsdtt3 finds no ocean seed at all and throws.
// =============================================================================
void run_oracle_dump(const std::string& path, double ocean_level)
{
    lsdtt::lsdtt_data_loader loader;
    loader.open_raster(path);
    lsdtt::lsdtt_raster_tile tile = loader.read_full();
    lsdtt::lsdtt_raster_cpu  dem;
    dem.set_from_vectors(std::move(tile.data), std::move(tile.mask), tile.meta);

    const lsdtt::lsdtt_raster_meta m0 = dem.get_meta();
    std::cout << "# input " << path << "  " << m0.n_rows << "x" << m0.n_cols
              << "  nodata " << m0.no_data_value
              << "  ocean_level " << ocean_level << "\n";

    // Upstream semantics exactly (Barnes .../src/main.cpp):
    //     label(i) = OCEAN  where  topo.isNoData(i) || topo(i) == ocean_level
    // and nothing else is ocean -- so the boundary-condition rule is switched
    // off. This is what option (A) of plan 4.2a bought us: no fixture surgery,
    // no cropping, both implementations asked the same question.
    lsdtt::depression_hierarchy_options opts;
    opts.ocean_from_boundary = false;
    opts.nodata_is_ocean     = true;
    opts.ocean_level         = static_cast<float>(ocean_level);

    std::vector<lsdtt::dep_label_t> labels;
    std::vector<uint8_t>            flowdirs;

    try
    {
        const lsdtt::lsdtt_depression_hierarchy hier =
            lsdtt::build_depression_hierarchy(dem, labels, flowdirs, opts);

        std::cout << "# nodes " << hier.size() << "\n";
        std::cout << "dep,parent,lchild,rchild,pit_elev,out_elev,cell_count,"
                     "ocean_parent\n";
        for (std::size_t d = 1; d < hier.size(); ++d)
        {
            const lsdtt::depression_node& n = hier[d];
            std::printf("%zu,%u,%u,%u,%.6g,%.6g,%u,%d\n",
                        d, n.parent, n.lchild, n.rchild,
                        static_cast<double>(n.pit_elev),
                        static_cast<double>(n.out_elev),
                        n.cell_count, n.ocean_parent ? 1 : 0);
        }
    }
    catch (const std::exception& e)
    {
        std::cout << "# THREW: " << e.what() << "\n";
    }
}

// =============================================================================
// ORACLE SUITE  (--oracle <data_dir> <expected_dir>)
//
// Phase 3 of docs/plan_docs/depression_hierarchy_fix_plan.md.
//
// Compares build_depression_hierarchy() field by field against expected output
// generated ONCE from the upstream Barnes implementation.  The expected files
// are committed; upstream is NOT built by lsdtt3 and is not needed to run this.
// See the plan's 3.3 for how they were produced and which commits they pin.
//
// Both sides are asked the same question via the explicit ocean options added in
// phase 2 (plan 4.2b): ocean = NODATA cells or cells at ocean_level, exactly as
// upstream's driver defines it, with the boundary-condition rule switched off.
//
// KNOWN DIVERGENCES.  Six fixtures disagree structurally: upstream records
// depressions joined through the ocean in a Depression::ocean_linked LIST, while
// our port instead creates a binary meta-depression parenting them.  Those are
// named below and reported as KNOWN, not silently skipped -- and a fixture NOT on
// the list that disagrees is a hard failure.  A listed fixture that starts
// agreeing is also reported, so the list gets pruned when phase 5 fixes the
// underlying divergence rather than quietly rotting.
// =============================================================================
//
// ALL of these are DEGENERATE-TIE cases, and none is a detection error.
//
// For every one of the 8, the LEAF SET is identical to upstream -- same
// depressions, same pit and outlet elevations, same cell counts -- so the BFS
// and the partition of cells agree exactly. Only the tree ABOVE the leaves
// differs, and in five of the 8 every outlet in the DEM shares a single
// elevation (`kerry_test1`: nine leaves, every out_elev 6), so there is no
// well-defined nesting to find: upstream's binary tree is an artefact of its
// hash-map iteration order, and ours of our tie-break. See plan 7.2.
//
// The leaf-set equality is ASSERTED below rather than skipped, so a genuine
// detection bug cannot hide behind this list.
//
// Historical note: phase 5b removed five other entries (kerry_test, kerry_test2,
// kerry_test3, kerry_test7, testdem2) by giving the outlet sort a total order. Upstream records depressions joined through the ocean
// in a Depression::ocean_linked LIST, leaving their `parent` as OCEAN; our port
// instead links them into the binary tree, either by chaining parents or by
// inventing a meta-depression. Two symptoms, same defect:
//
//   parent chain      kerry_test3: upstream parents are 0,0,0
//                     ours are 2,3,0
//   extra meta node   testdem2:    upstream has 2 ocean-linked leaves
//                     ours adds a meta-depression parenting both
//
// The other 12 fixtures agree field for field, so this is a contained defect in
// the ocean-link handling and not a general disagreement about the algorithm.
// Diagnosing it is phase 5; these fixtures are its test set.
const char* const k_known_divergent[] = {
    "kerry_test1",  "kerry_test4",  "kerry_test5",  "kerry_test6",
    "kerry_test10", "kerry_test11", "kerry_test12", "testdem10",
};

constexpr unsigned k_no_child = 4294967295u;   // upstream NO_VALUE

struct dep_row
{
    unsigned dep = 0, parent = 0, lchild = 0, rchild = 0, cell_count = 0;
    double   pit_elev = 0.0, out_elev = 0.0;
    int      ocean_parent = 0;
};

std::vector<dep_row> parse_expected(const std::string& path, bool& ok)
{
    std::vector<dep_row> rows;
    std::ifstream in(path);
    ok = in.good();
    if (!ok) return rows;
    std::string line;
    std::getline(in, line);                     // header
    while (std::getline(in, line))
    {
        if (line.empty()) continue;

        // Split on commas and convert each field with std::stod / std::stoul.
        // NOT with operator>> on a double: libstdc++'s num_get FAILS on the
        // token "inf", which upstream writes for every meta-depression's
        // pit_elev. That silently zeroed every field after it and produced a
        // suite full of bogus cell_count mismatches.
        std::vector<std::string> tok;
        std::string              cur;
        std::istringstream       ss(line);
        while (std::getline(ss, cur, ',')) tok.push_back(cur);
        if (tok.size() < 8) continue;

        dep_row r;
        r.dep          = static_cast<unsigned>(std::stoul(tok[0]));
        r.parent       = static_cast<unsigned>(std::stoul(tok[1]));
        r.lchild       = static_cast<unsigned>(std::stoul(tok[2]));
        r.rchild       = static_cast<unsigned>(std::stoul(tok[3]));
        r.pit_elev     = std::stod(tok[4]);
        r.out_elev     = std::stod(tok[5]);
        r.cell_count   = static_cast<unsigned>(std::stoul(tok[6]));
        r.ocean_parent = std::stoi(tok[7]);
        rows.push_back(r);
    }
    return rows;
}

void run_oracle_suite(const std::string& data_dir, const std::string& expected_dir)
{
    std::cout << "=== depression hierarchy: oracle suite vs upstream Barnes ===\n";

    // Every fixture in the upstream test_cases set.
    const char* const fixtures[] = {
        "kerry_test",  "kerry_test1", "kerry_test2", "kerry_test3", "kerry_test4",
        "kerry_test5", "kerry_test6", "kerry_test7", "kerry_test8", "kerry_test9",
        "kerry_test10","kerry_test11","kerry_test12",
        "testdem1", "testdem2", "testdem3",  "testdem4",  "testdem5",  "testdem6",
        "testdem7", "testdem8", "testdem9",  "testdem10", "testdem11", "testdem12",
    };

    for (const char* name : fixtures)
    {
        const std::string dem = data_dir + "/" + name + ".dem";
        const std::string exp = expected_dir + "/" + name + ".csv";

        bool known = false;
        for (const char* k : k_known_divergent)
            if (std::string(k) == name) known = true;

        bool have_expected = false;
        const std::vector<dep_row> want = parse_expected(exp, have_expected);
        if (!have_expected)
        {
            check(std::string("oracle ") + name + ": expected file present",
                  false, "missing " + exp);
            continue;
        }

        lsdtt::lsdtt_data_loader loader;
        loader.open_raster(dem);
        lsdtt::lsdtt_raster_tile tile = loader.read_full();
        lsdtt::lsdtt_raster_cpu  raster;
        raster.set_from_vectors(std::move(tile.data), std::move(tile.mask),
                                tile.meta);

        lsdtt::depression_hierarchy_options opts;
        opts.ocean_from_boundary = false;
        opts.nodata_is_ocean     = true;
        opts.ocean_level         = -9999.0f;

        std::vector<lsdtt::dep_label_t> labels;
        std::vector<uint8_t>            flowdirs;
        lsdtt::lsdtt_depression_hierarchy hier;
        try
        {
            hier = lsdtt::build_depression_hierarchy(raster, labels, flowdirs,
                                                     opts);
        }
        catch (const std::exception& e)
        {
            check(std::string("oracle ") + name + ": builds", false, e.what());
            continue;
        }

        const std::size_t got_n = hier.size() - 1;   // exclude the ocean node

        // Field-by-field comparison. Elevations are copied, never computed, so
        // they are compared exactly; a tolerance here would hide real errors.
        std::string first_diff;
        bool        agrees = (got_n == want.size());
        if (!agrees)
            first_diff = "count " + std::to_string(got_n) + " vs " +
                         std::to_string(want.size());

        if (agrees)
        {
            for (std::size_t i = 0; i < want.size(); ++i)
            {
                const dep_row&                w = want[i];
                const lsdtt::depression_node& g = hier[i + 1];
                const char* field = nullptr;
                if      (g.parent     != w.parent)     field = "parent";
                else if (g.lchild     != w.lchild)     field = "lchild";
                else if (g.rchild     != w.rchild)     field = "rchild";
                else if (g.cell_count != w.cell_count) field = "cell_count";
                else if ((g.ocean_parent ? 1 : 0) != w.ocean_parent)
                                                       field = "ocean_parent";
                // Compare in FLOAT, not double. The fixture stores %.9g, which
                // round-trips a float exactly -- but static_cast<double>(float)
                // is not equal to std::stod of that text, so comparing as
                // doubles fails on values that are in fact bit-identical.
                // Elevations are copied, never computed, so exact float
                // equality is the right test and no tolerance is warranted.
                else if (std::isfinite(w.pit_elev) &&
                         g.pit_elev != static_cast<float>(w.pit_elev))
                                                       field = "pit_elev";
                else if (std::isfinite(w.out_elev) &&
                         g.out_elev != static_cast<float>(w.out_elev))
                                                       field = "out_elev";
                if (field)
                {
                    agrees    = false;
                    first_diff = std::string("dep ") + std::to_string(i + 1) +
                                 " field " + field;
                    break;
                }
            }
        }

        if (known)
        {
            if (agrees)
            {
                report(std::string("oracle ") + name,
                       "KNOWN DIVERGENCE NOW AGREES -- remove it from "
                       "k_known_divergent");
                continue;
            }

            report(std::string("oracle ") + name,
                   "known divergence (" + first_diff + ")");

            // The tree above the leaves is arbitrary here, but the LEAF SET is
            // not: the depressions themselves, and the partition of cells
            // between them, are fully determined. Assert that much rather than
            // skipping the fixture outright -- otherwise a real detection bug
            // could hide behind a "known divergence" label.
            std::vector<std::array<double, 3>> want_leaves, got_leaves;

            for (const dep_row& w : want)
                if (w.lchild == k_no_child && w.rchild == k_no_child)
                    want_leaves.push_back({w.pit_elev, w.out_elev,
                                           static_cast<double>(w.cell_count)});

            for (std::size_t d = 1; d < hier.size(); ++d)
            {
                const lsdtt::depression_node& g = hier[d];
                if (g.lchild == lsdtt::k_no_label &&
                    g.rchild == lsdtt::k_no_label)
                    got_leaves.push_back({static_cast<double>(g.pit_elev),
                                          static_cast<double>(g.out_elev),
                                          static_cast<double>(g.cell_count)});
            }

            std::sort(want_leaves.begin(), want_leaves.end());
            std::sort(got_leaves.begin(),  got_leaves.end());

            check(std::string("oracle ") + name +
                      ": leaf depressions match upstream exactly",
                  want_leaves == got_leaves,
                  std::to_string(got_leaves.size()) + " leaves vs " +
                      std::to_string(want_leaves.size()));
        }
        else
        {
            check(std::string("oracle ") + name + " matches upstream",
                  agrees, agrees ? std::to_string(got_n) + " depressions"
                                 : first_diff);
        }
    }
}

}  // namespace

int main(int argc, char** argv)
{
    if (argc > 1 && std::string(argv[1]) == "--scaling")
    {
        run_scaling_harness();
        return 0;
    }
    if (argc > 3 && std::string(argv[1]) == "--oracle")
    {
        run_oracle_suite(argv[2], argv[3]);
        std::cout << "\n" << g_n_checks - g_n_failed << " / " << g_n_checks
                  << " checks passed\n";
        std::cout << "RESULT: " << (g_n_failed == 0 ? "PASS" : "FAIL") << "\n";
        return (g_n_failed == 0) ? 0 : 1;
    }
    if (argc > 2 && std::string(argv[1]) == "--dump")
    {
        const double ocean_level = (argc > 3) ? std::stod(argv[3]) : -9999.0;
        run_oracle_dump(argv[2], ocean_level);
        return 0;
    }

    std::cout << "=== depression hierarchy: phase 0 ===\n";
    std::cout << "These tests are EXPECTED TO FAIL on unmodified HEAD.\n"
                 "See docs/plan_docs/depression_hierarchy_fix_plan.md phase 0.\n";

    try
    {
        std::cout << "\n--- outlet_link copy semantics ---\n";
        test_outlet_link_copy_is_self_contained();

        std::cout << "\n--- outlet_db across a rehash ---\n";
        test_outlet_db_survives_rehash();

        std::cout << "\n--- the algorithm on a dense pit field ---\n";
        test_hierarchy_outlets_on_dense_pit_field();

        std::cout << "\n--- complexity guard ---\n";
        test_scaling_is_not_quadratic();
    }
    catch (const std::exception& e)
    {
        std::cout << "FAIL  unexpected exception: " << e.what() << std::endl;
        ++g_n_checks;
        ++g_n_failed;
    }

    std::cout << "\n" << g_n_checks - g_n_failed << " / " << g_n_checks
              << " checks passed\n";
    std::cout << "RESULT: " << (g_n_failed == 0 ? "PASS" : "FAIL") << "\n";

    return (g_n_failed == 0) ? 0 : 1;
}
