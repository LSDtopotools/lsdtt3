// lsdtt_reference_field_ladder_cpu.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html

#include "lsdtt_reference_field_ladder_cpu.hpp"
#include "lsdtt_raster_surface_fitting_cpu.hpp"
#include "../utils/lsdtt_random.hpp"
#include "../utils/lsdtt_grid_neighbours.hpp"
#include "../utils/lsdtt_statistics.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <queue>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

#include "lsdtt_flat_surface_channel_classify_cpu.hpp"   // find_population_breaks
#include "../utils/lsdtt_natural_breaks.hpp"

namespace lsdtt
{

namespace
{

// Float32 exact-integer ceiling. Above this a label cannot round-trip through
// the raster's backing store, so a join keyed on it silently corrupts.
constexpr std::int64_t k_max_exact_float_label = 16777216;

// -----------------------------------------------------------------------------
// band_of
//
// Band index for a value against ascending breaks: the number of breaks the
// value exceeds.
// -----------------------------------------------------------------------------
int band_of(double value, const std::vector<double>& breaks)
{
    int b = 0;
    for (const double lvl : breaks)
    {
        if (value > lvl) ++b;
        else break;
    }
    return b;
}

// -----------------------------------------------------------------------------
// evaluate_partition
//
// Band centres (medians) and the tightest gap between consecutive centres, for
// one candidate break set over one population.
//
// The gap is measured between CENTRES, not between breaks: two bands can be
// separated by a break and still hold populations whose typical elevations are
// indistinguishable, and it is the typical elevations that a claim about
// "levels" is about.
// -----------------------------------------------------------------------------
struct partition_metrics
{
    std::vector<double> band_centres;
    double              min_centre_gap = 0.0;
    bool                usable         = false;
};

partition_metrics evaluate_partition(
    const std::vector<double>& values,
    const std::vector<double>& breaks)
{
    partition_metrics pm;
    const std::size_t k = breaks.size() + 1;
    if (values.empty() || k < 1) return pm;

    std::vector<std::vector<double>> bins(k);
    for (const double v : values)
    {
        bins[static_cast<std::size_t>(band_of(v, breaks))].push_back(v);
    }

    // An empty band means the partition is not real: some level would be
    // claimed with nothing in it.
    for (const std::vector<double>& b : bins)
    {
        if (b.empty()) return pm;
    }

    pm.band_centres.reserve(k);
    for (std::vector<double>& b : bins)
    {
        pm.band_centres.push_back(compute_median_of_copy(std::move(b), 0.0));
    }
    std::sort(pm.band_centres.begin(), pm.band_centres.end());

    pm.min_centre_gap = std::numeric_limits<double>::infinity();
    for (std::size_t j = 1; j < pm.band_centres.size(); ++j)
    {
        pm.min_centre_gap = std::min(
            pm.min_centre_gap, pm.band_centres[j] - pm.band_centres[j - 1]);
    }
    if (pm.band_centres.size() == 1) pm.min_centre_gap = 0.0;

    pm.usable = true;
    return pm;
}

// -----------------------------------------------------------------------------
// jenks_breaks_for_k
// -----------------------------------------------------------------------------
std::vector<double> jenks_breaks_for_k(
    const std::vector<double>& values, int k, double& gvf_out)
{
    gvf_out = 0.0;
    std::vector<double> breaks;
    if (k < 2 || values.size() < static_cast<std::size_t>(k)) return breaks;

    lsdtt_jenks_options opts;
    opts.n_classes = k;
    const lsdtt_jenks_result jr = compute_jenks_breaks(values, opts);
    gvf_out = jr.goodness_of_variance_fit;

    // compute_jenks_breaks returns UPPER breaks including the top of the range;
    // the ladder wants only the k-1 interior cuts.
    breaks = jr.upper_breaks;
    if (!breaks.empty()) breaks.pop_back();
    if (static_cast<int>(breaks.size()) > k - 1)
    {
        breaks.resize(static_cast<std::size_t>(k - 1));
    }
    return breaks;
}

// -----------------------------------------------------------------------------
// window_histogram / separation_ratio_at
//
// The separation diagnostic. See ladder_window::separation_ratio for what the
// number means and why the resolvability floor does not replace it.
//
// The histogram is LINEAR in field units, not log10, because it is meant to be
// read against the breaks - and a break at 2 m and a break at 20 m have to sit
// on the same axis as the data. It is smoothed with a small moving average
// before the ratio is taken, for the same reason the antimode finder smooths:
// a raw histogram is bumpy enough that a single empty bin beside a break would
// otherwise read as a perfect separation.
// -----------------------------------------------------------------------------
void window_histogram(const std::vector<double>& values, int nbins,
                      double& lo, double& hi, std::vector<std::int64_t>& counts)
{
    counts.assign(static_cast<std::size_t>(std::max(nbins, 4)), 0);
    if (values.empty()) { lo = hi = 0.0; return; }

    lo = *std::min_element(values.begin(), values.end());
    hi = *std::max_element(values.begin(), values.end());
    if (!(hi > lo)) { hi = lo + 1e-6; }

    const double span = hi - lo;
    const int nb = static_cast<int>(counts.size());
    for (const double v : values)
    {
        int b = static_cast<int>((v - lo) / span * nb);
        b = std::clamp(b, 0, nb - 1);
        ++counts[static_cast<std::size_t>(b)];
    }
}

double separation_ratio_at(const std::vector<std::int64_t>& counts,
                           double lo, double hi, double break_value)
{
    const int nb = static_cast<int>(counts.size());
    if (nb < 4 || !(hi > lo)) return 1.0;

    // Smooth with a 3-bin moving average; the raw histogram is too bumpy for a
    // ratio taken at a single bin to mean anything.
    std::vector<double> s(static_cast<std::size_t>(nb), 0.0);
    for (int i = 0; i < nb; ++i)
    {
        double acc = 0.0; int m = 0;
        for (int d = -1; d <= 1; ++d)
        {
            const int j = i + d;
            if (j < 0 || j >= nb) continue;
            acc += static_cast<double>(counts[static_cast<std::size_t>(j)]);
            ++m;
        }
        s[static_cast<std::size_t>(i)] = (m > 0) ? acc / m : 0.0;
    }

    int bb = static_cast<int>((break_value - lo) / (hi - lo) * nb);
    bb = std::clamp(bb, 0, nb - 1);

    // Tallest peak on each side of the break. If either side is empty there is
    // no dip to speak of, so the break is not a separation.
    double peak_lo = 0.0, peak_hi = 0.0;
    for (int i = 0; i < bb; ++i)  peak_lo = std::max(peak_lo, s[static_cast<std::size_t>(i)]);
    for (int i = bb + 1; i < nb; ++i) peak_hi = std::max(peak_hi, s[static_cast<std::size_t>(i)]);
    const double flank = std::min(peak_lo, peak_hi);
    if (!(flank > 0.0)) return 1.0;

    return s[static_cast<std::size_t>(bb)] / flank;
}

// -----------------------------------------------------------------------------
// select_partition
//
// The k rule, and the whole reason this engine exists.
//
// A band count is ADMISSIBLE when the tightest gap between consecutive band
// CENTRES clears kappa * E_v. The answer is the LARGEST admissible k.
//
// EVERY k IN RANGE IS SCANNED, never searched downward from the top. Under
// median centres admissibility is not monotone in k: merging two well-separated
// bands can put the merged median right next to a neighbour's, so a coarser
// partition can fail where a finer one passed. A downward search would stop at
// the first failure and report a k that is not the largest admissible one.
//
// SELECTION IS OVER PARTITIONS, NOT OVER BREAKS. Dropping one break from a
// candidate set merges two bands and therefore moves EVERY centre and the
// minimum gap with them; the admissibility of the parts says nothing about the
// admissibility of the whole. So each candidate subset is rebuilt and
// re-evaluated in full.
//
// There is deliberately NO additive score of the form
//   w1 * GVF + w2 * (gap / gap_min).
// Goodness-of-variance-fit is bounded on [0, 1] while the band gap scales as
// range / k and is unbounded, so over a large vertical range the gap term
// dominates by orders of magnitude and the criterion collapses to the smallest
// k in the range regardless of what the data say. The rule here is a
// CONSTRAINT with an extremum over the feasible set. GVF is reported for
// diagnosis and is never an input to the choice.
// -----------------------------------------------------------------------------
struct selection
{
    std::vector<double> breaks;
    std::vector<double> band_centres;
    double              min_centre_gap = 0.0;
    double              gvf            = 0.0;
    int                 k              = 1;
    ladder_status       status         = ladder_status::degenerate;
};

// Budget on how many k-1 sized subsets of the candidate breaks are tried when
// the top-ranked subset is inadmissible. Pragmatic, not derived: with the
// candidate counts the antimode finder actually returns (rarely above 8) this
// is usually exhaustive.
constexpr int k_subset_budget = 16;

selection select_partition(
    const std::vector<double>& values,
    const std::vector<double>& candidate_breaks,
    const ladder_request&      request,
    double                     required_gap)
{
    selection best;

    if (values.size() < 2)
    {
        best.status = ladder_status::degenerate;
        return best;
    }

    const int max_k = std::max(1, request.max_level_count);

    // The window's own density, so a candidate partition can be tested for
    // whether its breaks sit in real gaps and not only for whether the DEM
    // could measure them. Built once per selection rather than per candidate.
    double sep_lo = 0.0, sep_hi = 0.0;
    std::vector<std::int64_t> sep_counts;
    if (request.max_separation_ratio > 0.0)
        window_histogram(values, request.window_hist_bins,
                         sep_lo, sep_hi, sep_counts);

    auto try_partition = [&](const std::vector<double>& brk) -> bool
    {
        const partition_metrics pm = evaluate_partition(values, brk);
        if (!pm.usable) return false;
        // Constraint 1: the instrument floor - can the DEM measure it.
        if (pm.min_centre_gap < required_gap) return false;
        // Constraint 2: separation - is there anything there to measure.
        // Two constraints, never blended into a score; see max_separation_ratio.
        if (request.max_separation_ratio > 0.0 && !sep_counts.empty())
        {
            for (const double b : brk)
                if (separation_ratio_at(sep_counts, sep_lo, sep_hi, b) >
                    request.max_separation_ratio)
                    return false;
        }
        return true;
    };

    // ---------------------------------------------------------------------
    // Explicit level count. Honoured, and reported honestly: an inadmissible
    // explicit k comes back with the requested k and status unresolved, never
    // silently replaced by a different number.
    // ---------------------------------------------------------------------
    if (request.level_count >= 2)
    {
        const int k = std::min(request.level_count, max_k);
        std::vector<double> brk;
        double gvf = 0.0;

        if (request.break_algorithm == ladder_break_choice::jenks)
        {
            brk = jenks_breaks_for_k(values, k, gvf);
        }
        else if (static_cast<int>(candidate_breaks.size()) >= k - 1)
        {
            // Take the k-1 strongest candidates, then test THAT partition; if
            // it fails, try other subsets under the budget before giving up.
            brk.assign(candidate_breaks.begin(),
                       candidate_breaks.begin() + (k - 1));
            if (!try_partition(brk))
            {
                std::vector<double> trial = candidate_breaks;
                int tried = 0;
                bool found = false;
                // Slide the window of k-1 consecutive candidates.
                const int n_cand = static_cast<int>(candidate_breaks.size());
                for (int start = 1;
                     start + (k - 1) <= n_cand && tried < k_subset_budget;
                     ++start, ++tried)
                {
                    std::vector<double> sub(
                        candidate_breaks.begin() + start,
                        candidate_breaks.begin() + start + (k - 1));
                    if (try_partition(sub))
                    {
                        brk = std::move(sub);
                        found = true;
                        break;
                    }
                }
                (void)found;
            }
        }
        else
        {
            brk = candidate_breaks;
        }

        const partition_metrics pm = evaluate_partition(values, brk);
        best.breaks         = brk;
        best.band_centres   = pm.band_centres;
        best.min_centre_gap = pm.usable ? pm.min_centre_gap : 0.0;
        best.gvf            = gvf;
        best.k              = static_cast<int>(brk.size()) + 1;
        best.status         = (pm.usable && try_partition(brk))
            ? ladder_status::forced
            : ladder_status::unresolved;
        return best;
    }

    // ---------------------------------------------------------------------
    // Auto: the largest admissible k. Scan every k in range.
    // ---------------------------------------------------------------------
    for (int k = 2; k <= max_k; ++k)
    {
        std::vector<double> brk;
        double gvf = 0.0;

        if (request.break_algorithm == ladder_break_choice::jenks)
        {
            brk = jenks_breaks_for_k(values, k, gvf);
        }
        else
        {
            if (static_cast<int>(candidate_breaks.size()) < k - 1) break;
            brk.assign(candidate_breaks.begin(),
                       candidate_breaks.begin() + (k - 1));
        }
        if (static_cast<int>(brk.size()) != k - 1) continue;

        const partition_metrics pm = evaluate_partition(values, brk);
        if (!pm.usable) continue;
        // BOTH constraints, via the one predicate. Checking the floor inline
        // here is how the separation constraint came to be silently bypassed on
        // the auto path: it applied to the explicit-k branch and to nothing
        // else. Any new admissibility rule belongs in try_partition, not here.
        if (!try_partition(brk)) continue;

        // Admissible. Keep it and keep going: a larger k may also be
        // admissible, and the rule wants the largest.
        best.breaks         = brk;
        best.band_centres   = pm.band_centres;
        best.min_centre_gap = pm.min_centre_gap;
        best.gvf            = gvf;
        best.k              = k;
        best.status         = ladder_status::resolved;
    }

    if (best.status != ladder_status::resolved)
    {
        // No k qualifies. One band, and say so - this is not "k = 1 is the
        // answer", it is "the DEM cannot resolve a ladder here".
        const partition_metrics pm = evaluate_partition(values, {});
        best.breaks.clear();
        best.band_centres   = pm.band_centres;
        best.min_centre_gap = 0.0;
        best.k              = 1;
        best.status         = ladder_status::unresolved;
    }
    return best;
}




// -----------------------------------------------------------------------------
// gap_break_candidates
//
// Breaks from the EMPTY SPACE between clusters, which is what a reader's eye
// uses on a HAND histogram and what neither other finder looks for. See
// ladder_break_choice in the header for why jenks and antimode both miss it.
//
// Bins linearly and emits two kinds of candidate (see the header for why both
// are needed):
//
//   EMPTY RUNS - contiguous bins below empty_fraction of the tallest bin,
//                spanning at least min_width. Break placed per gap_break_at.
//   TROUGHS    - interior local minima of the 3-bin-smoothed histogram whose
//                floor sits below gap_trough_max_ratio of the smaller flanking
//                peak. "Empty" is relative to the TALLEST bin, so where one
//                level dominates the population - a mountain floodplain, say -
//                a genuine dip above it never reads as empty and only the
//                sparse tail does. Without troughs the finder then has no
//                candidate below 100 m and bands the terraces as valley floor.
//
// Leading and trailing runs are dropped: the empty space below the lowest datum
// and above the highest is not a boundary between two populations. Ranking is
// gap_rank: lowest break first by default, so the ladder is walked bottom-up.
// -----------------------------------------------------------------------------
std::vector<double> gap_break_candidates(
    const std::vector<double>& values, const ladder_request& request,
    std::vector<double>* gap_lo, std::vector<double>* gap_hi)
{
    std::vector<double> out;
    if (gap_lo) gap_lo->clear();
    if (gap_hi) gap_hi->clear();
    if (values.size() < 2) return out;

    const int nb = std::max(request.gap_hist_bins, 8);
    const double lo = *std::min_element(values.begin(), values.end());
    const double hi = *std::max_element(values.begin(), values.end());
    if (!(hi > lo)) return out;
    const double span = hi - lo;
    const double bw   = span / static_cast<double>(nb);

    std::vector<std::int64_t> counts(static_cast<std::size_t>(nb), 0);
    for (const double v : values)
    {
        int b = static_cast<int>((v - lo) / span * nb);
        b = std::clamp(b, 0, nb - 1);
        ++counts[static_cast<std::size_t>(b)];
    }

    const std::int64_t peak = *std::max_element(counts.begin(), counts.end());
    if (peak <= 0) return out;
    const double empty_at =
        std::max(0.0, request.gap_empty_fraction) * static_cast<double>(peak);

    struct cand { double at; double width; double lo; double hi; };
    std::vector<cand> cands;

    // Contiguous runs of empty bins.
    struct run { int lo_bin; int hi_bin; };
    std::vector<run> runs;
    int start = -1;
    for (int i = 0; i < nb; ++i)
    {
        const bool empty = static_cast<double>(counts[static_cast<std::size_t>(i)]) <= empty_at;
        if (empty && start < 0) start = i;
        if (!empty && start >= 0) { runs.push_back({start, i - 1}); start = -1; }
    }
    if (start >= 0) runs.push_back({start, nb - 1});

    std::vector<char> in_run(static_cast<std::size_t>(nb), 0);
    for (const run& r : runs)
    {
        // Drop runs touching either end: empty space below the lowest datum or
        // above the highest is not a boundary between two populations.
        if (r.lo_bin == 0 || r.hi_bin == nb - 1) continue;
        const double glo = lo + bw * static_cast<double>(r.lo_bin);
        const double ghi = lo + bw * static_cast<double>(r.hi_bin + 1);
        const double w   = ghi - glo;
        if (w < request.gap_min_width_metres) continue;
        for (int i = r.lo_bin; i <= r.hi_bin; ++i)
            in_run[static_cast<std::size_t>(i)] = 1;
        double at = 0.5 * (glo + ghi);
        if (request.gap_break_at == gap_break_position::lower)      at = glo;
        else if (request.gap_break_at == gap_break_position::upper) at = ghi;
        cands.push_back({at, w, glo, ghi});
    }

    // Troughs: interior local minima of the smoothed histogram. Smoothed with
    // the same 3-bin average separation_ratio_at uses, so a candidate offered
    // here is scored on the same curve it will later be judged against.
    if (request.gap_trough_max_ratio > 0.0)
    {
        std::vector<double> s(static_cast<std::size_t>(nb), 0.0);
        for (int i = 0; i < nb; ++i)
        {
            double acc = 0.0; int m = 0;
            for (int d = -1; d <= 1; ++d)
            {
                const int j = i + d;
                if (j < 0 || j >= nb) continue;
                acc += static_cast<double>(counts[static_cast<std::size_t>(j)]);
                ++m;
            }
            s[static_cast<std::size_t>(i)] = (m > 0) ? acc / m : 0.0;
        }

        for (int i = 1; i < nb - 1; ++i)
        {
            if (in_run[static_cast<std::size_t>(i)]) continue;   // already a gap
            const double c = s[static_cast<std::size_t>(i)];
            if (c > s[static_cast<std::size_t>(i - 1)]) continue;
            if (c > s[static_cast<std::size_t>(i + 1)]) continue;
            // A flat shoulder is not a trough: require a strict fall on one side.
            if (!(c < s[static_cast<std::size_t>(i - 1)] ||
                  c < s[static_cast<std::size_t>(i + 1)])) continue;

            double peak_lo = 0.0, peak_hi = 0.0;
            for (int j = 0; j < i; ++j)
                peak_lo = std::max(peak_lo, s[static_cast<std::size_t>(j)]);
            for (int j = i + 1; j < nb; ++j)
                peak_hi = std::max(peak_hi, s[static_cast<std::size_t>(j)]);
            const double flank = std::min(peak_lo, peak_hi);
            if (!(flank > 0.0)) continue;
            if (c / flank > request.gap_trough_max_ratio) continue;

            // The trough has no span, so lo and hi are the bin it sits in and
            // gap_break_at has nothing to choose between.
            const double blo = lo + bw * static_cast<double>(i);
            const double bhi = blo + bw;
            cands.push_back({0.5 * (blo + bhi), 0.0, blo, bhi});
        }
    }

    if (request.gap_rank == gap_rank_choice::width)
    {
        std::sort(cands.begin(), cands.end(),
                  [](const cand& a, const cand& b) { return a.width > b.width; });
    }
    else
    {
        // Bottom-up. Ties (a trough sitting at the same height as a run) go to
        // the wider one, which is the better-separated statement of the two.
        std::sort(cands.begin(), cands.end(),
                  [](const cand& a, const cand& b) {
                      if (a.at != b.at) return a.at < b.at;
                      return a.width > b.width;
                  });
    }
    for (const cand& c : cands)
    {
        out.push_back(c.at);
        if (gap_lo) gap_lo->push_back(c.lo);
        if (gap_hi) gap_hi->push_back(c.hi);
    }
    return out;
}

// -----------------------------------------------------------------------------
// clamp_breaks
//
// Apply the hard min / max bounds on accepted breaks.
// -----------------------------------------------------------------------------
std::vector<double> clamp_breaks(
    std::vector<double> breaks, const ladder_request& request)
{
    std::vector<double> out;
    out.reserve(breaks.size());
    for (const double b : breaks)
    {
        if (request.min_break_metres > 0.0 && b < request.min_break_metres) continue;
        if (request.max_break_metres > 0.0 && b > request.max_break_metres) continue;
        out.push_back(b);
    }
    std::sort(out.begin(), out.end());
    out.erase(std::unique(out.begin(), out.end()), out.end());
    return out;
}

// -----------------------------------------------------------------------------
// bootstrap_k
//
// Resample a window's population and re-run the whole break search on each
// replicate, so k_selected can be reported with a stability rather than on its
// own. See ladder_request::k_bootstrap for what this costs and why it is off by
// default.
//
// Resampling is WITH REPLACEMENT at the population's own size, capped at
// max_sample. The cap is not an approximation to apologise for: every finder
// bins before it searches, so past a few thousand values the extra draws change
// no histogram this looks at, and the cost would be spent for nothing.
//
// The generator is a fixed-seed LCG rather than <random>, so a replicate count
// reproduces exactly across platforms and standard libraries. A stability
// figure that moved with the toolchain would be worse than none.
// -----------------------------------------------------------------------------
struct bootstrap_summary
{
    int    mode      = -1;
    double stability = 0.0;
    int    ci_lo     = -1;
    int    ci_hi     = -1;
    int    n_run     = 0;
};

bootstrap_summary bootstrap_k(
    const std::vector<double>& values,
    const ladder_request&      request,
    double                     required_gap,
    std::uint64_t              seed)
{
    bootstrap_summary out;
    if (request.k_bootstrap < 1 || values.size() < 2) return out;

    const std::size_t cap = static_cast<std::size_t>(
        std::max(2, request.k_bootstrap_max_sample));
    const std::size_t draw = std::min(values.size(), cap);

    // splitmix64: cheap, well-distributed, and fully specified (lsdtt_random)
    // so the sequence does not depend on the standard library.
    random::splitmix64_stream rng(seed ? seed : 1u);
    auto next = [&rng]() -> std::uint64_t { return rng.next(); };

    std::vector<int> ks;
    ks.reserve(static_cast<std::size_t>(request.k_bootstrap));
    std::vector<double> sample(draw, 0.0);

    for (int b = 0; b < request.k_bootstrap; ++b)
    {
        for (std::size_t i = 0; i < draw; ++i)
            sample[i] = values[next() % values.size()];

        std::vector<double> cand;
        if (request.break_algorithm == ladder_break_choice::gap)
        {
            cand = clamp_breaks(
                gap_break_candidates(sample, request, nullptr, nullptr),
                request);
        }
        else if (request.break_algorithm == ladder_break_choice::antimode)
        {
            const int nb = std::max(request.hist_bins, 16);
            const double floor_m = std::max(request.min_field_metres, 1e-3);
            const double lo_l = std::log10(floor_m);
            double hi_l = lo_l;
            for (const double v : sample)
                hi_l = std::max(hi_l, std::log10(std::max(v, floor_m)));
            hi_l += 1e-6;
            std::vector<std::int64_t> h(static_cast<std::size_t>(nb), 0);
            for (const double v : sample)
            {
                const double l = std::log10(std::max(v, floor_m));
                int bb = static_cast<int>((l - lo_l) / (hi_l - lo_l) * nb);
                bb = std::clamp(bb, 0, nb - 1);
                ++h[static_cast<std::size_t>(bb)];
            }
            cand = clamp_breaks(
                find_population_breaks(h, lo_l, hi_l,
                                       request.max_valley_fraction),
                request);
        }
        // jenks takes k directly, so select_partition builds its own breaks.

        const selection sel =
            select_partition(sample, cand, request, required_gap);
        ks.push_back(sel.k);
    }

    if (ks.empty()) return out;
    out.n_run = static_cast<int>(ks.size());

    std::unordered_map<int, int> tally;
    for (const int k : ks) ++tally[k];
    int best_k = ks.front(), best_n = 0;
    for (const auto& [k, cnt] : tally)
    {
        if (cnt > best_n || (cnt == best_n && k < best_k))
        {
            best_k = k; best_n = cnt;
        }
    }
    out.mode      = best_k;
    out.stability = static_cast<double>(best_n) / static_cast<double>(ks.size());

    std::sort(ks.begin(), ks.end());
    auto pct = [&](double f) -> int
    {
        std::size_t i = static_cast<std::size_t>(
            f * static_cast<double>(ks.size() - 1) + 0.5);
        if (i >= ks.size()) i = ks.size() - 1;
        return ks[i];
    };
    out.ci_lo = pct(0.025);
    out.ci_hi = pct(0.975);
    return out;
}

}  // namespace

// =============================================================================
// build_fitted_surface_reference_field
// =============================================================================
fitted_surface_field build_fitted_surface_reference_field(
    const lsdtt_raster_cpu&       dem,
    const lsdtt_raster_cpu&       fit_mask,
    const fitted_surface_request& request)
{
    static const char* fn = "lsdtt::build_fitted_surface_reference_field";

    const lsdtt_raster_meta& meta = dem.get_meta();
    const std::size_t n =
        static_cast<std::size_t>(meta.n_rows) * static_cast<std::size_t>(meta.n_cols);

    if (fit_mask.get_data().size() != n)
    {
        throw std::invalid_argument(
            std::string(fn) + ": fit_mask does not match dem");
    }

    fitted_surface_field out;
    out.algorithm_used = request.algorithm;

    if (request.algorithm == fitted_surface_choice::cone)
    {
        double apex_x = request.apex_x;
        double apex_y = request.apex_y;
        if (!request.apex_is_set &&
            !find_cone_apex(dem, fit_mask, request.apex_fraction, apex_x, apex_y))
        {
            return out;   // the mask selects nothing; valid stays false
        }

        const lsdtt_cone_fit_result fit =
            fit_cone_to_mask(dem, fit_mask, apex_x, apex_y);
        if (!fit.valid) return out;

        out.cone_coeffs   = fit.coeffs;
        out.rmse          = fit.rmse;
        out.n_pixels_used = fit.n_pixels_used;
        out.field         = residual_from_fitted_cone(
            dem, fit.coeffs, request.n_threads);
        out.valid         = true;
        return out;
    }

    if (request.algorithm == fitted_surface_choice::poly2)
    {
        const lsdtt_poly2_fit_result fit = fit_poly2_to_mask(dem, fit_mask);
        if (!fit.valid) return out;   // valid stays false; caller decides

        out.poly2_coeffs  = fit.coeffs;
        out.rmse          = fit.rmse;
        out.n_pixels_used = fit.n_pixels_used;
        out.field         = residual_from_fitted_poly2(
            dem, fit.coeffs, request.n_threads);
        out.valid         = true;
        return out;
    }

    const lsdtt_plane_fit_result fit = fit_plane_to_mask_ransac(
        dem, fit_mask,
        request.ransac_tolerance_metres,
        request.ransac_iterations,
        request.seed);
    if (!fit.valid) return out;

    out.plane_coeffs  = fit.coeffs;
    out.rmse          = fit.rmse;
    out.n_pixels_used = fit.n_pixels_used;
    out.field         = residual_from_fitted_plane(
        dem, fit.coeffs, request.n_threads);
    out.valid         = true;
    return out;
}

// =============================================================================
// build_reference_field_ladder
// =============================================================================
ladder_result build_reference_field_ladder(
    const lsdtt_raster_cpu&          reference_field,
    const lsdtt_raster_cpu&          admitted_mask,
    const ladder_request&            request,
    const std::vector<std::int32_t>* window_source_id,
    const lsdtt_raster_cpu*          window_coordinate,
    const lsdtt_raster_cpu*          slope_degrees)
{
    static const char* fn = "lsdtt::build_reference_field_ladder";

    const lsdtt_raster_meta& meta = reference_field.get_meta();
    const std::size_t n =
        static_cast<std::size_t>(meta.n_rows) * static_cast<std::size_t>(meta.n_cols);

    if (admitted_mask.get_data().size() != n)
    {
        throw std::invalid_argument(
            std::string(fn) + ": admitted_mask does not match reference_field");
    }

    const std::vector<float>&   fv = reference_field.get_data();
    const std::vector<uint8_t>& fm = reference_field.get_mask();
    const std::vector<float>&   am = admitted_mask.get_data();

    ladder_result result;
    lsdtt_raster_meta band_meta = meta;
    band_meta.no_data_value = -1.0;

    // Band index raster: -1 everywhere a pixel is not admitted or the field is
    // undefined, so a consumer can never mistake "not banded" for band 0.
    std::vector<float> band_data(n, -1.0f);
    std::vector<uint8_t> band_mask(n, 1);

    auto is_live = [&](std::size_t i) -> bool
    {
        return am[i] > 0.5f && fm[i] != 0 && std::isfinite(fv[i]);
    };

    // -------------------------------------------------------------------------
    // Population, and the global break candidates.
    // -------------------------------------------------------------------------
    std::vector<double> all_values;
    all_values.reserve(n / 4);
    double lo = std::numeric_limits<double>::infinity();
    double hi = -std::numeric_limits<double>::infinity();
    for (std::size_t i = 0; i < n; ++i)
    {
        if (!is_live(i)) continue;
        const double v = static_cast<double>(fv[i]);
        all_values.push_back(v);
        lo = std::min(lo, v);
        hi = std::max(hi, v);
    }

    // -------------------------------------------------------------------------
    // The resolvability floor. Under `fixed` this is one number for the whole
    // raster; under `slope_dependent` each window reads its own E_v off the
    // median gradient of its OWN admitted cells, because required_gap
    // constrains the separation between band CENTRES and a band centre is a
    // property of a partition, not of a place - it has no slope of its own.
    // -------------------------------------------------------------------------
    const bool want_slope_ev =
        request.vertical_accuracy_mode ==
            vertical_accuracy_choice::slope_dependent;

    const std::vector<float>*   sld = nullptr;
    const std::vector<uint8_t>* slm = nullptr;
    if (want_slope_ev && slope_degrees != nullptr &&
        slope_degrees->get_data().size() == n)
    {
        sld = &slope_degrees->get_data();
        slm = &slope_degrees->get_mask();
    }
    else if (want_slope_ev)
    {
        result.notes.push_back(
            "vertical_accuracy_mode is slope_dependent but no usable gradient "
            "raster was supplied, so the fixed vertical_accuracy_metres was "
            "used for every window.");
    }

    // E_v for a population, given the gradients of the same cells.
    auto ev_for = [&](const std::vector<double>& slopes) -> double
    {
        if (sld == nullptr || slopes.empty())
            return request.vertical_accuracy_metres;
        const double med = compute_median_of_copy(slopes);
        return (med <= request.vertical_accuracy_break_degrees)
            ? request.vertical_accuracy_gentle_metres
            : request.vertical_accuracy_steep_metres;
    };
    auto median_slope_of = [&](const std::vector<double>& slopes) -> double
    {
        if (sld == nullptr || slopes.empty())
            return std::numeric_limits<double>::quiet_NaN();
        return compute_median_of_copy(slopes);
    };
    auto gap_from_ev = [&](double ev) -> double
    {
        return std::max(0.0, request.resolvability_kappa * ev);
    };

    // Gradients of the whole admitted population, for the global partition and
    // for every window that has to fall back on it.
    std::vector<double> all_slopes;
    if (sld != nullptr)
    {
        all_slopes.reserve(all_values.size());
        for (std::size_t i = 0; i < n; ++i)
        {
            if (!is_live(i)) continue;
            if (!(*slm)[i] || !std::isfinite((*sld)[i])) continue;
            all_slopes.push_back(static_cast<double>((*sld)[i]));
        }
    }

    const double global_ev  = ev_for(all_slopes);
    const double required_gap = gap_from_ev(global_ev);

    // The antimode finder bins on log10 with a small positive floor, so every
    // negative value lands in the lowest bin. On a signed field - a
    // fitted-surface residual, where below-surface IS the signal - that folds
    // the entire below-surface population into one spike and the breaks that
    // come back mean nothing. Refuse rather than return a plausible-looking
    // ladder built on a destroyed histogram. See ladder_break_choice.
    // gap and jenks both work on the linear field, so only antimode is caught.
    if (request.break_algorithm == ladder_break_choice::antimode &&
        std::isfinite(lo) && lo < 0.0)
    {
        throw std::invalid_argument(
            std::string(fn) + ": the antimode break finder bins on log10 and "
            "cannot band a field containing negative values (minimum here is " +
            std::to_string(lo) + "). Use ladder_break_choice::jenks, which "
            "works on the linear field.");
    }

    if (all_values.size() < 2 || !(hi > lo))
    {
        // Nothing to stratify. Everything live becomes band 0.
        for (std::size_t i = 0; i < n; ++i)
        {
            if (is_live(i)) band_data[i] = 0.0f;
        }
        ladder_window w;
        w.n_pixels      = static_cast<std::int64_t>(all_values.size());
        w.status        = ladder_status::degenerate;
        w.required_gap_m = required_gap;
        w.vertical_accuracy_m = global_ev;
        result.windows.push_back(std::move(w));
        result.n_windows_degenerate = 1;
        result.band_index = lsdtt_raster_cpu(
            std::move(band_data), std::move(band_mask), band_meta);
        return result;
    }

    // Candidate breaks from the antimode finder, on a log10 histogram, exactly
    // as the valley engine does. Wrapped, never modified.
    const int    nbins   = std::max(request.hist_bins, 16);
    const double floor_m = std::max(request.min_field_metres, 1e-3);
    const double lo_log  = std::log10(floor_m);
    double       hi_log  = lo_log;
    for (const double v : all_values)
    {
        hi_log = std::max(hi_log, std::log10(std::max(v, floor_m)));
    }
    hi_log += 1e-6;

    auto log_bin = [&](double v) -> int
    {
        const double l = std::log10(std::max(v, floor_m));
        int b = static_cast<int>((l - lo_log) / (hi_log - lo_log) * nbins);
        return std::clamp(b, 0, nbins - 1);
    };

    {
        std::vector<std::int64_t> ghist(static_cast<std::size_t>(nbins), 0);
        for (const double v : all_values) ++ghist[static_cast<std::size_t>(log_bin(v))];
        result.global_breaks_m = clamp_breaks(
            request.break_algorithm == ladder_break_choice::gap
                ? gap_break_candidates(all_values, request, nullptr, nullptr)
                : find_population_breaks(ghist, lo_log, hi_log,
                                         request.max_valley_fraction),
            request);
    }
    const selection global_selection = select_partition(
        all_values, result.global_breaks_m, request, required_gap);

    // -------------------------------------------------------------------------
    // Windowing. A fan has no along-channel coordinate, so the fitted-surface
    // field is always banded globally regardless of reach_window_metres.
    // -------------------------------------------------------------------------
    const bool want_windows =
        request.reference_field == reference_field_choice::hand &&
        request.reach_window_metres >= 0.0 &&
        window_source_id != nullptr &&
        window_coordinate != nullptr &&
        window_source_id->size() == n &&
        window_coordinate->get_data().size() == n;

    if (!want_windows)
    {
        // One global ladder.
        const selection& sel = global_selection;

        for (std::size_t i = 0; i < n; ++i)
        {
            if (!is_live(i)) continue;
            band_data[i] = static_cast<float>(
                band_of(static_cast<double>(fv[i]), sel.breaks));
        }

        ladder_window w;
        w.source_id             = -1;
        w.window_centre_m       = 0.0;
        w.n_pixels              = static_cast<std::int64_t>(all_values.size());
        w.k_selected            = sel.k;
        w.status                = sel.status;
        w.used_global_breaks    = true;
        w.breaks_m              = sel.breaks;
        {
            double hlo = 0.0, hhi = 0.0;
            std::vector<std::int64_t> hc;
            window_histogram(all_values, request.window_hist_bins,
                             hlo, hhi, hc);
            for (const double b : sel.breaks)
                w.separation_ratio.push_back(
                    separation_ratio_at(hc, hlo, hhi, b));
            if (request.record_window_histograms)
            {
                w.hist_lo_m   = hlo;
                w.hist_hi_m   = hhi;
                w.hist_counts = std::move(hc);
            }
        }
        w.band_centres_m        = sel.band_centres;
        w.min_band_centre_gap_m = sel.min_centre_gap;
        w.required_gap_m        = required_gap;
        w.vertical_accuracy_m   = global_ev;
        w.median_slope_degrees  = median_slope_of(all_slopes);
        w.goodness_of_variance_fit = sel.gvf;
        {
            const bootstrap_summary bs = bootstrap_k(
                all_values, request, required_gap,
                static_cast<std::uint64_t>(request.k_bootstrap_seed));
            w.k_bootstrap_mode = bs.mode;
            w.k_stability      = bs.stability;
            w.k_ci_lo          = bs.ci_lo;
            w.k_ci_hi          = bs.ci_hi;
            w.n_bootstrap      = bs.n_run;
        }
        result.windows.push_back(std::move(w));

        result.max_band_index = std::max(0, sel.k - 1);
        if (sel.status == ladder_status::resolved ||
            sel.status == ladder_status::forced)      ++result.n_windows_resolved;
        else if (sel.status == ladder_status::unresolved) ++result.n_windows_unresolved;
        else                                              ++result.n_windows_degenerate;

        result.band_index = lsdtt_raster_cpu(
            std::move(band_data), std::move(band_mask), band_meta);
        return result;
    }

    // -------------------------------------------------------------------------
    // Reach-relative ladders. Windows of fixed length along each source's
    // channel; a window pools its own histogram with its two neighbours so a
    // sparse reach still has a population to work with.
    // -------------------------------------------------------------------------
    const std::vector<float>&   cd = window_coordinate->get_data();
    const std::vector<uint8_t>& cm = window_coordinate->get_mask();
    const double cell = meta.cell_size;

    struct source_bins
    {
        double dmin      = 0.0;
        double wstep     = 1.0;
        int    n_windows = 1;
        std::vector<std::vector<double>> window_values;
        // Gradients of the same cells, in the same order, so a window's E_v is
        // read off the ground the window actually covers.
        std::vector<std::vector<double>> window_slopes;
    };

    // Per-source coordinate range.
    std::unordered_map<std::int32_t, std::pair<double, double>> drange;
    for (std::size_t i = 0; i < n; ++i)
    {
        if (!is_live(i)) continue;
        if (!cm[i] || !std::isfinite(cd[i])) continue;
        const std::int32_t s = (*window_source_id)[i];
        if (s < 0) continue;
        const double d = static_cast<double>(cd[i]);
        auto it = drange.find(s);
        if (it == drange.end()) drange[s] = {d, d};
        else
        {
            it->second.first  = std::min(it->second.first, d);
            it->second.second = std::max(it->second.second, d);
        }
    }

    std::unordered_map<std::int32_t, source_bins> bins;
    for (const auto& [s, r] : drange)
    {
        source_bins sb;
        sb.dmin = r.first;
        const double span = std::max(r.second - r.first, 1.0);
        double w = request.reach_window_metres;
        if (w <= 0.0) w = std::max(2000.0, span / 8.0);
        sb.wstep     = std::max(w * 0.5, cell);
        sb.n_windows = static_cast<int>(span / sb.wstep) + 1;
        sb.window_values.resize(static_cast<std::size_t>(sb.n_windows));
        sb.window_slopes.resize(static_cast<std::size_t>(sb.n_windows));
        bins.emplace(s, std::move(sb));
    }

    auto window_of = [](const source_bins& sb, double d) -> int
    {
        const int w = static_cast<int>((d - sb.dmin) / sb.wstep);
        return std::clamp(w, 0, sb.n_windows - 1);
    };

    for (std::size_t i = 0; i < n; ++i)
    {
        if (!is_live(i)) continue;
        if (!cm[i] || !std::isfinite(cd[i])) continue;
        const std::int32_t s = (*window_source_id)[i];
        if (s < 0) continue;
        source_bins& sb = bins.at(s);
        const std::size_t w = static_cast<std::size_t>(
            window_of(sb, static_cast<double>(cd[i])));
        sb.window_values[w].push_back(static_cast<double>(fv[i]));
        if (sld != nullptr && (*slm)[i] && std::isfinite((*sld)[i]))
        {
            sb.window_slopes[w].push_back(static_cast<double>((*sld)[i]));
        }
    }

    // Per (source, window) selection.
    std::unordered_map<std::int32_t, std::vector<std::vector<double>>> wbreaks;
    for (auto& [s, sb] : bins)
    {
        std::vector<std::vector<double>>& wb = wbreaks[s];
        wb.resize(static_cast<std::size_t>(sb.n_windows));
        std::vector<char> has_own(static_cast<std::size_t>(sb.n_windows), 0);
        std::vector<selection> selected(static_cast<std::size_t>(sb.n_windows));
        const std::size_t first_rec = result.windows.size();

        for (int w = 0; w < sb.n_windows; ++w)
        {
            // Pool with the immediate neighbours. The gradients are pooled
            // over exactly the same windows, so E_v is read off the ground the
            // breaks were actually chosen from rather than off a wider reach.
            std::vector<double> pooled;
            std::vector<double> pooled_slopes;
            for (int k = std::max(w - 1, 0);
                 k <= std::min(w + 1, sb.n_windows - 1); ++k)
            {
                const std::vector<double>& src =
                    sb.window_values[static_cast<std::size_t>(k)];
                pooled.insert(pooled.end(), src.begin(), src.end());
                const std::vector<double>& ssrc =
                    sb.window_slopes[static_cast<std::size_t>(k)];
                pooled_slopes.insert(pooled_slopes.end(),
                                     ssrc.begin(), ssrc.end());
            }

            const double win_ev  = ev_for(pooled_slopes);
            const double win_gap = gap_from_ev(win_ev);

            ladder_window rec;
            rec.source_id       = s;
            rec.window_centre_m = sb.dmin + (w + 0.5) * sb.wstep;
            rec.n_pixels        = static_cast<std::int64_t>(pooled.size());
            rec.required_gap_m  = win_gap;
            rec.vertical_accuracy_m  = win_ev;
            rec.median_slope_degrees = median_slope_of(pooled_slopes);

            if (static_cast<int>(pooled.size()) >= request.min_window_pixels)
            {
                std::vector<std::int64_t> h(static_cast<std::size_t>(nbins), 0);
                for (const double v : pooled)
                    ++h[static_cast<std::size_t>(log_bin(v))];

                const std::vector<double> cand = clamp_breaks(
                    request.break_algorithm == ladder_break_choice::gap
                        ? gap_break_candidates(pooled, request,
                                               &rec.gap_lo_m, &rec.gap_hi_m)
                        : find_population_breaks(h, lo_log, hi_log,
                                                 request.max_valley_fraction),
                    request);

                const selection sel =
                    select_partition(pooled, cand, request, win_gap);

                rec.candidate_breaks_m = cand;
                selected[static_cast<std::size_t>(w)] = sel;
                has_own[static_cast<std::size_t>(w)] = 1;

                rec.k_selected            = sel.k;
                rec.status                = sel.status;
                rec.breaks_m              = sel.breaks;
                rec.band_centres_m        = sel.band_centres;
                rec.min_band_centre_gap_m = sel.min_centre_gap;
                rec.goodness_of_variance_fit = sel.gvf;

                // Seeded per (source, window) so a window's stability is
                // reproducible on its own, and does not shift because a
                // neighbouring reach changed.
                {
                    const std::uint64_t wseed =
                        static_cast<std::uint64_t>(request.k_bootstrap_seed) * 1000003ull
                        + static_cast<std::uint64_t>(static_cast<std::uint32_t>(s)) * 10007ull
                        + static_cast<std::uint64_t>(w) + 1ull;
                    const bootstrap_summary bs =
                        bootstrap_k(pooled, request, win_gap, wseed);
                    rec.k_bootstrap_mode = bs.mode;
                    rec.k_stability      = bs.stability;
                    rec.k_ci_lo          = bs.ci_lo;
                    rec.k_ci_hi          = bs.ci_hi;
                    rec.n_bootstrap      = bs.n_run;
                }

                // How cleanly the density actually separates at each break.
                // Computed regardless of the break algorithm: jenks has no
                // notion of a dip and will split a single mode without
                // complaint, which is exactly when this number matters.
                {
                    double hlo = 0.0, hhi = 0.0;
                    std::vector<std::int64_t> hc;
                    window_histogram(pooled, request.window_hist_bins,
                                     hlo, hhi, hc);
                    rec.separation_ratio.clear();
                    for (const double b : sel.breaks)
                        rec.separation_ratio.push_back(
                            separation_ratio_at(hc, hlo, hhi, b));
                    if (request.record_window_histograms)
                    {
                        rec.hist_lo_m    = hlo;
                        rec.hist_hi_m    = hhi;
                        // The 3-bin smoothing the separation ratio is taken on,
                        // so a reader sees the curve the decision used.
                        const int nbh = static_cast<int>(hc.size());
                        rec.hist_smoothed.assign(hc.size(), 0.0);
                        for (int q = 0; q < nbh; ++q)
                        {
                            double acc = 0.0; int m = 0;
                            for (int dq = -1; dq <= 1; ++dq)
                            {
                                const int r2 = q + dq;
                                if (r2 < 0 || r2 >= nbh) continue;
                                acc += static_cast<double>(hc[static_cast<std::size_t>(r2)]);
                                ++m;
                            }
                            rec.hist_smoothed[static_cast<std::size_t>(q)] =
                                (m > 0) ? acc / m : 0.0;
                        }
                        rec.hist_counts  = std::move(hc);
                    }
                }

                wb[static_cast<std::size_t>(w)] = sel.breaks;
            }
            else
            {
                rec.status = ladder_status::degenerate;
            }

            result.windows.push_back(std::move(rec));
        }

        // A sparse window inherits the complete selection from the nearest
        // adequately populated window on the same source, or the complete
        // global selection when that source has none. An adequately populated
        // unresolved window owns its one-band selection and is never replaced
        // by unrelated multi-band breaks.
        for (int w = 0; w < sb.n_windows; ++w)
        {
            if (has_own[static_cast<std::size_t>(w)]) continue;
            int best = -1;
            for (int k = 0; k < sb.n_windows; ++k)
            {
                if (has_own[static_cast<std::size_t>(k)] &&
                    (best < 0 || std::abs(k - w) < std::abs(best - w)))
                {
                    best = k;
                }
            }
            const selection& fallback = (best >= 0)
                ? selected[static_cast<std::size_t>(best)]
                : global_selection;
            ladder_window& rec =
                result.windows[first_rec + static_cast<std::size_t>(w)];
            rec.k_selected = fallback.k;
            rec.breaks_m = fallback.breaks;
            rec.band_centres_m = fallback.band_centres;
            rec.min_band_centre_gap_m = fallback.min_centre_gap;
            rec.goodness_of_variance_fit = fallback.gvf;
            rec.used_global_breaks = best < 0;
            wb[static_cast<std::size_t>(w)] = fallback.breaks;
        }
    }

    // Per-pixel band index from the pixel's own window. Pixels without a usable
    // donor coordinate remain admitted and use the global partition.
    // Which window each pixel drew its breaks from. Keyed on (source, window)
    // so it survives the deterministic sort of `windows` below.
    std::vector<float>   widx(n, -1.0f);
    std::vector<uint8_t> widx_mask(n, 1u);
    std::vector<std::pair<std::int32_t,int>> pixel_win(n, {-1, -1});

    for (std::size_t i = 0; i < n; ++i)
    {
        if (!is_live(i)) continue;
        const std::vector<double>* brk = &global_selection.breaks;
        if (cm[i] && std::isfinite(cd[i]))
        {
            const std::int32_t s = (*window_source_id)[i];
            const auto bit = bins.find(s);
            if (s >= 0 && bit != bins.end())
            {
                const source_bins& sb = bit->second;
                const int wi = window_of(sb, static_cast<double>(cd[i]));
                brk = &wbreaks.at(s)[static_cast<std::size_t>(wi)];
                pixel_win[i] = {s, wi};
            }
        }
        const int b = band_of(static_cast<double>(fv[i]), *brk);
        band_data[i] = static_cast<float>(b);
        result.max_band_index = std::max(result.max_band_index, b);
    }

    // The value is the window's ORDINAL WITHIN ITS SOURCE, counting downstream
    // - a reach number, not an index into the sorted window table. That keeps
    // it meaningful on its own: window 3 of channel 1 is the fourth reach down
    // channel 1, whatever order the table happens to be in.
    for (std::size_t i = 0; i < n; ++i)
        if (pixel_win[i].second >= 0)
            widx[i] = static_cast<float>(pixel_win[i].second);
    result.window_index = lsdtt_raster_cpu(
        std::move(widx), std::move(widx_mask), band_meta);

    std::sort(result.windows.begin(), result.windows.end(),
        [](const ladder_window& a, const ladder_window& b)
        {
            if (a.source_id != b.source_id) return a.source_id < b.source_id;
            return a.window_centre_m < b.window_centre_m;
        });

    for (const ladder_window& w : result.windows)
    {
        if (w.status == ladder_status::resolved ||
            w.status == ladder_status::forced)          ++result.n_windows_resolved;
        else if (w.status == ladder_status::unresolved) ++result.n_windows_unresolved;
        else                                            ++result.n_windows_degenerate;
    }

    result.band_index = lsdtt_raster_cpu(
        std::move(band_data), std::move(band_mask), band_meta);
    return result;
}

// =============================================================================
// label_components_within_bands
// =============================================================================
banded_label_result label_components_within_bands(
    const lsdtt_raster_cpu& band_index,
    const lsdtt_raster_cpu& admitted_mask,
    const lsdtt_raster_cpu& reference_field,
    int                     min_component_pixels,
    pixel_connectivity      connectivity,
    const lsdtt_raster_cpu* terrain_aspect_degrees)
{
    static const char* fn = "lsdtt::label_components_within_bands";

    const lsdtt_raster_meta& meta = band_index.get_meta();
    const int n_rows = meta.n_rows;
    const int n_cols = meta.n_cols;
    const std::size_t n =
        static_cast<std::size_t>(n_rows) * static_cast<std::size_t>(n_cols);

    if (admitted_mask.get_data().size() != n ||
        reference_field.get_data().size() != n)
    {
        throw std::invalid_argument(
            std::string(fn) + ": input rasters are not co-registered");
    }

    const std::vector<float>& bi = band_index.get_data();
    const std::vector<float>& am = admitted_mask.get_data();
    const std::vector<float>& fv = reference_field.get_data();

    int max_band = -1;
    for (std::size_t i = 0; i < n; ++i)
    {
        if (am[i] > 0.5f && bi[i] >= 0.0f)
        {
            max_band = std::max(max_band, static_cast<int>(bi[i]));
        }
    }

    std::vector<float> label_data(n, 0.0f);
    banded_label_result out;

    const int  n_off = (connectivity == pixel_connectivity::d8) ? 8 : 4;
    static constexpr const auto& dr8 = grid::d8_rowmajor_row;
    static constexpr const auto& dc8 = grid::d8_rowmajor_col;
    static constexpr const auto& dr4 = grid::d4_nwes_row;
    static constexpr const auto& dc4 = grid::d4_nwes_col;
    const int* dr = (connectivity == pixel_connectivity::d8) ? dr8 : dr4;
    const int* dc = (connectivity == pixel_connectivity::d8) ? dc8 : dc4;

    const double cell_area = meta.cell_size * meta.cell_size;

    // A GLOBAL RUNNING COUNTER, advanced band by band. The band is not packed
    // into the label; it travels in band_index and in the patch record.
    std::int64_t next_label = 1;

    for (int band = 0; band <= max_band; ++band)
    {
        for (std::size_t seed = 0; seed < n; ++seed)
        {
            if (am[seed] <= 0.5f) continue;
            if (static_cast<int>(bi[seed]) != band) continue;
            if (bi[seed] < 0.0f) continue;
            if (label_data[seed] != 0.0f) continue;

            if (next_label >= k_max_exact_float_label)
            {
                throw std::invalid_argument(
                    std::string(fn) + ": patch count reached the float32 exact-"
                    "integer ceiling (" +
                    std::to_string(k_max_exact_float_label) +
                    "); labels above it cannot round-trip through the raster "
                    "backing store and every join keyed on them would silently "
                    "corrupt");
            }

            // Flood fill this component, within this band only.
            std::vector<std::size_t> members;
            std::queue<std::size_t> q;
            q.push(seed);
            label_data[seed] = static_cast<float>(next_label);

            while (!q.empty())
            {
                const std::size_t cur = q.front();
                q.pop();
                members.push_back(cur);

                const int r = static_cast<int>(cur / static_cast<std::size_t>(n_cols));
                const int c = static_cast<int>(cur % static_cast<std::size_t>(n_cols));

                for (int k = 0; k < n_off; ++k)
                {
                    const int nr = r + dr[k];
                    const int nc = c + dc[k];
                    if (nr < 0 || nr >= n_rows || nc < 0 || nc >= n_cols) continue;
                    const std::size_t ni =
                        static_cast<std::size_t>(nr) * n_cols + nc;
                    if (label_data[ni] != 0.0f) continue;
                    if (am[ni] <= 0.5f) continue;
                    if (bi[ni] < 0.0f) continue;
                    if (static_cast<int>(bi[ni]) != band) continue;
                    label_data[ni] = static_cast<float>(next_label);
                    q.push(ni);
                }
            }

            if (static_cast<int>(members.size()) < min_component_pixels)
            {
                // Too small to keep: erase it and reuse the number.
                for (const std::size_t m : members) label_data[m] = 0.0f;
                continue;
            }

            patch_record rec;
            rec.label    = next_label;
            rec.band     = band;
            rec.n_pixels = static_cast<std::int64_t>(members.size());
            rec.area_m2  = static_cast<double>(members.size()) * cell_area;

            std::vector<double> vals;
            std::vector<double> xs, ys;
            vals.reserve(members.size());
            xs.reserve(members.size());
            ys.reserve(members.size());
            for (const std::size_t m : members)
            {
                if (!std::isfinite(fv[m])) continue;
                const int r = static_cast<int>(m / static_cast<std::size_t>(n_cols));
                const int c = static_cast<int>(m % static_cast<std::size_t>(n_cols));
                xs.push_back(static_cast<double>(c) * meta.cell_size);
                ys.push_back(-static_cast<double>(r) * meta.cell_size);
                vals.push_back(static_cast<double>(fv[m]));
            }
            rec.median_field = compute_median_of_copy(vals, 0.0);

            // Patch tilt from the shared least-squares plane (decision D7,
            // 2026-10-08): it uses a relative degeneracy test,
            // |det| < 1e-10 * max(Sxx * Syy, 1), where this code had its own
            // fit with an absolute |det| > 1e-12.
            if (vals.size() >= 3)
            {
                const lsdtt_plane_fit_result plane =
                    fit_plane_least_squares(xs, ys, vals);
                if (plane.valid)
                {
                    const double a = plane.coeffs.a;
                    const double b = plane.coeffs.b;
                    if (std::hypot(a, b) > 1e-12)
                    {
                        constexpr double rad_to_deg =
                            180.0 / 3.14159265358979323846;
                        rec.aspect_degrees = std::atan2(-a, -b) * rad_to_deg;
                        if (rec.aspect_degrees < 0.0)
                            rec.aspect_degrees += 360.0;
                    }
                }
            }

            // Terrain aspect: CIRCULAR mean plus the resultant length that
            // says whether that mean means anything. Summing unit vectors is
            // the only correct way to average a direction - the arithmetic mean
            // of 350 and 10 degrees is 180, the opposite way.
            if (terrain_aspect_degrees != nullptr &&
                terrain_aspect_degrees->get_data().size() == n)
            {
                const std::vector<float>&   ad = terrain_aspect_degrees->get_data();
                const std::vector<uint8_t>& aspect_mask =
                    terrain_aspect_degrees->get_mask();
                constexpr double deg_to_rad = 3.14159265358979323846 / 180.0;
                double sn = 0.0, cs = 0.0;
                std::int64_t na = 0;
                for (const std::size_t m : members)
                {
                    if (!aspect_mask[m] || !std::isfinite(ad[m])) continue;
                    const double a_rad = static_cast<double>(ad[m]) * deg_to_rad;
                    sn += std::sin(a_rad);
                    cs += std::cos(a_rad);
                    ++na;
                }
                if (na > 0)
                {
                    const double inv = 1.0 / static_cast<double>(na);
                    const double mn = sn * inv, mc = cs * inv;
                    rec.aspect_r_bar = std::sqrt(mn * mn + mc * mc);
                    double deg = std::atan2(mn, mc) *
                                 (180.0 / 3.14159265358979323846);
                    // fmod, not a single += 360: a mean that lands a hair below
                    // zero - which is exactly what a set straddling north
                    // gives - would otherwise be reported as 360 degrees, and a
                    // reader binning aspect into octants would put it in the
                    // wrong bin or off the end of the array.
                    deg = std::fmod(deg + 360.0, 360.0);
                    rec.terrain_aspect_degrees = deg;
                }
            }

            out.patches.push_back(std::move(rec));
            ++next_label;
        }
    }

    std::vector<uint8_t> label_mask(n, 1);
    out.labels = lsdtt_raster_cpu(
        std::move(label_data), std::move(label_mask), meta);
    return out;
}

} // namespace lsdtt
