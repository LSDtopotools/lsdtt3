// test_depression_hierarchy.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
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
#include <iostream>
#include <limits>
#include <string>
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

}  // namespace

int main(int argc, char** argv)
{
    if (argc > 1 && std::string(argv[1]) == "--scaling")
    {
        run_scaling_harness();
        return 0;
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
