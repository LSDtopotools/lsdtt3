// Minimal OpenMP scaling ceiling test.
//
// Standalone by design: no lsdtt3 headers, no CMake, no dependency beyond a
// compiler with OpenMP. It is meant to be copied to a machine that has no
// lsdtt3 build on it and answer, in about ten seconds, "can this machine
// parallelise anything at all?"
//
//     g++ -O3 -DNDEBUG -fopenmp omp_ceiling.cpp -o omp_ceiling && ./omp_ceiling
//
// Two kernels at opposite ends of the arithmetic-to-memory-traffic ratio:
//
//   compute-bound : many sqrt() per element, data stays in registers, so the
//                   working set is irrelevant and stays small and fixed
//   memory-bound  : streaming triad, almost no arithmetic, sized at runtime to
//                   4x the machine's last-level cache so it really does stream
//
// The point is not the absolute times. It is that a machine which parallelises
// the first kernel and not the second has no threading fault: a bandwidth-
// limited kernel has nothing for extra threads to win. Interpreting an
// application's scaling curve without this ceiling beside it is how a
// memory-bound kernel gets misread as a threading bug.
//
// It is NOT a CTest test and nothing registers it: tests/CMakeLists.txt names
// every cpp/*.cpp target explicitly and does not glob, so this file is inert
// until someone compiles it by hand. It asserts nothing and cannot fail a
// build. Reference measurements at 1 / 2 / 4 threads:
//
//   i5-1145G7, 4 cores / 8 threads, 8 MiB L3, GCC 11.4, WSL2:
//     compute   1.00x   1.76x   3.07x
//     memory    1.00x   0.97x   0.96x
//
//   EPYC 9654P, 96 cores / 192 threads, 384 MiB L3, GCC 13.3, 8-CPU Slurm
//   allocation:
//     compute   1.00x   2.00x   4.00x
//
// The laptop figures came off a thermally limited machine under WSL2 and drift
// between runs; they are a shape to recognise, not a target to hit. What they
// say is that compute scales with cores and streaming does not.
//
// READING IT
//
// The COMPUTE row is the diagnostic. It should approach cores-many speedup on
// any healthy machine. A machine that does not reach it is throttled,
// cpuset-limited or oversubscribed, and no application scaling number measured
// on it means anything until that is resolved. The same EPYC above reported
// 1.62x at BOTH 2 and 4 threads while a login-node affinity mask held it to
// `Cpus_allowed_list: 1,97` -- two SMT siblings of one physical core -- and
// that single line of output is what identified it.
//
// The MEMORY row is context, not a health check. A flat curve there is the
// correct answer and not a fault.
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <algorithm>
#include <omp.h>

static const int kThreadList[] = {1, 2, 4};
static const int kReps = 7;          // median of 7: cheap, and robust to one
                                     // scheduler hiccup on a laptop

static double median(std::vector<double> v)
{
    std::sort(v.begin(), v.end());
    return v[v.size() / 2];
}

// -----------------------------------------------------------------------------
// Compute-bound: 60 dependent sqrt per element. One load, one store, so the
// arithmetic dominates by orders of magnitude and the kernel is a near-pure
// test of how many cores actually execute.
// -----------------------------------------------------------------------------
static double compute_bound(std::vector<double>& a, int threads)
{
    const std::size_t n = a.size();
    const double t0 = omp_get_wtime();
#pragma omp parallel for num_threads(threads) schedule(static)
    for (std::size_t i = 0; i < n; ++i)
    {
        double x = a[i];
        for (int k = 0; k < 60; ++k) x = std::sqrt(x + 1.0);
        a[i] = x;
    }
    return omp_get_wtime() - t0;
}

// -----------------------------------------------------------------------------
// Memory-bound: streaming triad. Three arrays touched per element, two reads
// and one write, one multiply-add. Limited by DRAM bandwidth, which is shared
// by every core and does not grow when threads are added.
// -----------------------------------------------------------------------------
// The working set MUST exceed last-level cache or this kernel is not a memory
// test at all. A fixed size cannot do that across machines: 1M doubles is a
// 24 MB triad, which is DRAM-bound against an 8 MiB laptop L3 and entirely
// cache-resident against the 384 MiB L3 of an EPYC 9654P. The latter duly
// reported 1.97x/4.06x -- a perfect "bandwidth" curve that was really an L3
// curve, and exactly the false reassurance this file exists to prevent.
//
// So read the real L3 and size the arrays to a multiple of it.
static std::size_t last_level_cache_bytes()
{
    // sysfs, in preference to anything derived from the CPU model: it is what
    // the kernel actually sees, and it reports a VM's or container's view
    // rather than the bare metal's.
    std::size_t best = 0;
    for (int idx = 0; idx < 10; ++idx)
    {
        char path[128];
        std::snprintf(path, sizeof path,
                      "/sys/devices/system/cpu/cpu0/cache/index%d/size", idx);
        std::FILE* f = std::fopen(path, "r");
        if (!f) continue;
        long value = 0; char unit = 0;
        if (std::fscanf(f, "%ld%c", &value, &unit) >= 1)
        {
            std::size_t bytes = static_cast<std::size_t>(value);
            if (unit == 'K') bytes *= 1024;
            else if (unit == 'M') bytes *= 1024 * 1024;
            else if (unit == 'G') bytes *= 1024ull * 1024 * 1024;
            best = std::max(best, bytes);
        }
        std::fclose(f);
    }
    return best;
}

static double memory_bound(std::vector<double>& a,
                           const std::vector<double>& b,
                           const std::vector<double>& c,
                           int threads, int passes)
{
    const std::size_t n = a.size();
    const double t0 = omp_get_wtime();
    // One pass is over in a millisecond or two, short enough that timer
    // granularity and a single scheduler hiccup dominate -- measured swinging
    // between 1.05x and 0.59x on one machine. Repeating inside the timed
    // region lengthens the measurement without changing the working set.
    for (int pass = 0; pass < passes; ++pass)
    {
#pragma omp parallel for num_threads(threads) schedule(static)
        for (std::size_t i = 0; i < n; ++i)
            a[i] = b[i] + 3.0 * c[i];
    }
    return omp_get_wtime() - t0;
}

int main(int argc, char** argv)
{
    // The compute kernel is size-independent -- 60 sqrt per element, data in
    // registers -- so it keeps a small fixed array and stays quick.
    const std::size_t n_compute = 1000000;

    // The memory kernel is sized from last-level cache. 4x LLC is enough that
    // the triad cannot be served from cache even with a generous prefetcher.
    // Override with a first argument in MiB (per array) when sysfs is absent
    // or the automatic size is inconvenient.
    const std::size_t llc = last_level_cache_bytes();
    std::size_t triad_bytes = (llc ? llc * 4 : 128ull * 1024 * 1024);
    if (argc > 1)
    {
        const long mib = std::strtol(argv[1], nullptr, 10);
        if (mib > 0) triad_bytes = static_cast<std::size_t>(mib) * 1024 * 1024 * 3;
    }
    // Cap at 3 GB of triad so an auto-sized run on a huge-cache part does not
    // become an allocation experiment.
    triad_bytes = std::min<std::size_t>(triad_bytes, 3ull * 1024 * 1024 * 1024);
    const std::size_t n_memory = std::max<std::size_t>(triad_bytes / (3 * sizeof(double)),
                                                       1000000);

    // Enough passes for roughly 4 GB of traffic, so the timed region is tens
    // of milliseconds on any machine regardless of the size chosen above.
    const std::size_t bytes_per_pass = n_memory * 3 * sizeof(double);
    const int passes = static_cast<int>(
        std::clamp<std::size_t>(4ull * 1024 * 1024 * 1024 / bytes_per_pass, 3, 200));

    std::printf("logical cores reported by OpenMP: %d\n", omp_get_max_threads());
    if (llc)
        std::printf("last-level cache               : %.1f MiB\n",
                    llc / (1024.0 * 1024.0));
    else
        std::printf("last-level cache               : UNKNOWN (sysfs unreadable)\n");
    std::printf("memory-kernel working set      : %.1f MiB over %d passes\n",
                bytes_per_pass / (1024.0 * 1024.0), passes);
    if (llc && bytes_per_pass < llc * 2)
        std::printf("  !! WARNING: working set is not comfortably larger than "
                    "cache.\n  !! The memory row below measures CACHE bandwidth "
                    "and WILL appear to scale.\n");
    std::printf("\n");

    std::vector<double> b(n_memory), c(n_memory);
    for (std::size_t i = 0; i < n_memory; ++i) { b[i] = 1.0 + i % 97; c[i] = 2.0 + i % 89; }

    struct Result { double t[3]; };
    Result comp{}, mem{};

    for (int ti = 0; ti < 3; ++ti)
    {
        const int threads = kThreadList[ti];

        std::vector<double> times;
        for (int r = 0; r < kReps; ++r)
        {
            std::vector<double> a(n_compute, 1.0);
            times.push_back(compute_bound(a, threads));
        }
        comp.t[ti] = median(times);

        times.clear();
        for (int r = 0; r < kReps; ++r)
        {
            std::vector<double> a(n_memory, 0.0);
            times.push_back(memory_bound(a, b, c, threads, passes));
        }
        mem.t[ti] = median(times);
    }

    std::printf("%-16s %10s %10s %10s\n", "kernel", "1t", "2t", "4t");
    std::printf("%-16s %9.4fs %9.4fs %9.4fs\n", "compute (s)",
                comp.t[0], comp.t[1], comp.t[2]);
    std::printf("%-16s %9.2fx %9.2fx %9.2fx\n", "compute speedup",
                1.0, comp.t[0] / comp.t[1], comp.t[0] / comp.t[2]);
    std::printf("%-16s %9.4fs %9.4fs %9.4fs\n", "memory (s)",
                mem.t[0], mem.t[1], mem.t[2]);
    std::printf("%-16s %9.2fx %9.2fx %9.2fx\n", "memory speedup",
                1.0, mem.t[0] / mem.t[1], mem.t[0] / mem.t[2]);
    return 0;
}
