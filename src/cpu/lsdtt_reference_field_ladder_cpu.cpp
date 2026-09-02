// lsdtt_reference_field_ladder_cpu.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html

#include "lsdtt_reference_field_ladder_cpu.hpp"

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

double median_of_copy(std::vector<double> v)
{
    if (v.empty()) return 0.0;
    const std::size_t mid = v.size() / 2;
    std::nth_element(v.begin(), v.begin() + mid, v.end());
    return v[mid];
}

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
        pm.band_centres.push_back(median_of_copy(std::move(b)));
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

    auto try_partition = [&](const std::vector<double>& brk) -> bool
    {
        const partition_metrics pm = evaluate_partition(values, brk);
        if (!pm.usable) return false;
        return pm.min_centre_gap >= required_gap;
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
        best.status         = (pm.usable && pm.min_centre_gap >= required_gap)
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
        if (pm.min_centre_gap < required_gap) continue;

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

}  // namespace

// =============================================================================
// build_reference_field_ladder
// =============================================================================
ladder_result build_reference_field_ladder(
    const lsdtt_raster_cpu&          reference_field,
    const lsdtt_raster_cpu&          admitted_mask,
    const ladder_request&            request,
    const std::vector<std::int32_t>* window_source_id,
    const lsdtt_raster_cpu*          window_coordinate)
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

    const double required_gap =
        std::max(0.0, request.resolvability_kappa * request.vertical_accuracy_metres);

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
        result.windows.push_back(std::move(w));
        result.n_windows_degenerate = 1;
        result.band_index = lsdtt_raster_cpu(
            std::move(band_data), std::move(band_mask), meta);
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
            find_population_breaks(ghist, lo_log, hi_log,
                                   request.max_valley_fraction),
            request);
    }

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
        const selection sel = select_partition(
            all_values, result.global_breaks_m, request, required_gap);

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
        w.band_centres_m        = sel.band_centres;
        w.min_band_centre_gap_m = sel.min_centre_gap;
        w.required_gap_m        = required_gap;
        w.goodness_of_variance_fit = sel.gvf;
        result.windows.push_back(std::move(w));

        result.max_band_index = std::max(0, sel.k - 1);
        if (sel.status == ladder_status::resolved ||
            sel.status == ladder_status::forced)      ++result.n_windows_resolved;
        else if (sel.status == ladder_status::unresolved) ++result.n_windows_unresolved;
        else                                              ++result.n_windows_degenerate;

        result.band_index = lsdtt_raster_cpu(
            std::move(band_data), std::move(band_mask), meta);
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
        sb.window_values[static_cast<std::size_t>(
            window_of(sb, static_cast<double>(cd[i])))]
            .push_back(static_cast<double>(fv[i]));
    }

    // Per (source, window) selection.
    std::unordered_map<std::int32_t, std::vector<std::vector<double>>> wbreaks;
    for (auto& [s, sb] : bins)
    {
        std::vector<std::vector<double>>& wb = wbreaks[s];
        wb.resize(static_cast<std::size_t>(sb.n_windows));
        std::vector<char> has_own(static_cast<std::size_t>(sb.n_windows), 0);
        const std::size_t first_rec = result.windows.size();

        for (int w = 0; w < sb.n_windows; ++w)
        {
            // Pool with the immediate neighbours.
            std::vector<double> pooled;
            for (int k = std::max(w - 1, 0);
                 k <= std::min(w + 1, sb.n_windows - 1); ++k)
            {
                const std::vector<double>& src =
                    sb.window_values[static_cast<std::size_t>(k)];
                pooled.insert(pooled.end(), src.begin(), src.end());
            }

            ladder_window rec;
            rec.source_id       = s;
            rec.window_centre_m = sb.dmin + (w + 0.5) * sb.wstep;
            rec.n_pixels        = static_cast<std::int64_t>(pooled.size());
            rec.required_gap_m  = required_gap;

            if (static_cast<int>(pooled.size()) >= request.min_window_pixels)
            {
                std::vector<std::int64_t> h(static_cast<std::size_t>(nbins), 0);
                for (const double v : pooled)
                    ++h[static_cast<std::size_t>(log_bin(v))];

                const std::vector<double> cand = clamp_breaks(
                    find_population_breaks(h, lo_log, hi_log,
                                           request.max_valley_fraction),
                    request);

                const selection sel =
                    select_partition(pooled, cand, request, required_gap);

                rec.k_selected            = sel.k;
                rec.status                = sel.status;
                rec.breaks_m              = sel.breaks;
                rec.band_centres_m        = sel.band_centres;
                rec.min_band_centre_gap_m = sel.min_centre_gap;
                rec.goodness_of_variance_fit = sel.gvf;

                if (!sel.breaks.empty())
                {
                    wb[static_cast<std::size_t>(w)] = sel.breaks;
                    has_own[static_cast<std::size_t>(w)] = 1;
                }
            }
            else
            {
                rec.status = ladder_status::degenerate;
            }

            rec.used_global_breaks = !has_own[static_cast<std::size_t>(w)];
            result.windows.push_back(std::move(rec));
        }

        // A window without its own breaks inherits the nearest window on the
        // same source that has some, and falls back to the global set only
        // when that source has none anywhere.
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
            wb[static_cast<std::size_t>(w)] = (best >= 0)
                ? wb[static_cast<std::size_t>(best)]
                : result.global_breaks_m;
            result.windows[first_rec + static_cast<std::size_t>(w)].breaks_m =
                wb[static_cast<std::size_t>(w)];
        }
    }

    // Per-pixel band index from the pixel's own window.
    for (std::size_t i = 0; i < n; ++i)
    {
        if (!is_live(i)) continue;
        if (!cm[i] || !std::isfinite(cd[i])) continue;
        const std::int32_t s = (*window_source_id)[i];
        if (s < 0) continue;
        const source_bins& sb = bins.at(s);
        const std::vector<double>& brk =
            wbreaks.at(s)[static_cast<std::size_t>(
                window_of(sb, static_cast<double>(cd[i])))];
        const int b = band_of(static_cast<double>(fv[i]), brk);
        band_data[i] = static_cast<float>(b);
        result.max_band_index = std::max(result.max_band_index, b);
    }

    for (const ladder_window& w : result.windows)
    {
        if (w.status == ladder_status::resolved ||
            w.status == ladder_status::forced)          ++result.n_windows_resolved;
        else if (w.status == ladder_status::unresolved) ++result.n_windows_unresolved;
        else                                            ++result.n_windows_degenerate;
    }

    result.band_index = lsdtt_raster_cpu(
        std::move(band_data), std::move(band_mask), meta);
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
    pixel_connectivity      connectivity)
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
    const int  dr8[8] = {-1, -1, -1,  0, 0,  1, 1, 1};
    const int  dc8[8] = {-1,  0,  1, -1, 1, -1, 0, 1};
    const int  dr4[4] = {-1,  0,  0,  1};
    const int  dc4[4] = { 0, -1,  1,  0};
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
            vals.reserve(members.size());
            for (const std::size_t m : members)
            {
                if (std::isfinite(fv[m])) vals.push_back(static_cast<double>(fv[m]));
            }
            rec.median_field = median_of_copy(std::move(vals));

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
