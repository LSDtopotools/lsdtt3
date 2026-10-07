// test_order_statistics.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
//
// -----------------------------------------------------------------------------
// Characterisation test for the in-place order statistics in
// src/utils/lsdtt_statistics (compute_median_in_place,
// compute_median_upper_in_place, compute_quantile_interpolated_in_place,
// compute_quantile_floor_in_place).
//
// Each replaced a family of private helpers. Every removed helper is kept below
// VERBATIM (the body, wrapped in a function where it was inline code), with the
// file and line it came from, and is compared with the shared function exactly
// as the migrated site now calls it:
//
//   - results compared by bit pattern (so -0.0 vs 0.0 is a difference);
//   - for sites that select in place and then use the vector again, the order
//     the vector is left in is compared too;
//   - inputs: n = 0..41 and 1000, heavy ties (values rounded to a coarse grid),
//     float and double; for the two filtering percentiles, +-inf and NaN too.
//
// Signed zeros are run as a separate regime. A sort-based helper and a
// selection-based one pick the same order statistic but not necessarily the
// same zero when +0 and -0 are tied, and this does reach the result at the five
// sites whose algorithm changed (terrace, paired_flowpaths, sti, chi_network:
// sort or two-pass selection -> one selection; raster_algebra: selection ->
// sort). Simon accepted that as D7 on 2026-10-07, so in that regime those five
// must agree in value (+0 == -0) and every other site still bit for bit.
//
// Usage: lsdtt-test-order-statistics
// Prints one line per check and a final "RESULT: PASS|FAIL".
// -----------------------------------------------------------------------------

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <functional>
#include <iostream>
#include <limits>
#include <random>
#include <string>
#include <vector>

#include "../../src/utils/lsdtt_statistics.hpp"

using namespace lsdtt;

namespace
{

const double k_nan = std::numeric_limits<double>::quiet_NaN();

// ============================================================================
// Removed copies (references)
// ============================================================================

// ---- mean of the two middle values -----------------------------------------

// src/cpu/lsdtt_terrace_edge_tracer_cpu.cpp:40-58 (as of 2026-10-07)
static double terrace_median_copy_(std::vector<double> vals)
{
    if (vals.empty())
        return std::numeric_limits<double>::quiet_NaN();

    const std::size_t n = vals.size();
    const std::size_t mid = n / 2;

    std::nth_element(vals.begin(), vals.begin() + mid, vals.end());
    double med = vals[mid];

    if (n % 2 == 0)
    {
        std::nth_element(vals.begin(), vals.begin() + mid - 1, vals.end());
        med = 0.5 * (med + vals[mid - 1]);
    }

    return med;
}

// src/cpu/lsdtt_paired_flowpaths_cpu.cpp:35-41 (as of 2026-10-07)
double paired_median_of(std::vector<double> v)
{
    if (v.empty()) return std::numeric_limits<double>::quiet_NaN();
    std::sort(v.begin(), v.end());
    const std::size_t n = v.size();
    return (n % 2) ? v[n / 2] : 0.5 * (v[n / 2 - 1] + v[n / 2]);
}

// src/cpu/lsdtt_sti_cpu.cpp:314-321 (as of 2026-10-07)
double sti_median_in_place(std::vector<double>& values)
{
    if (values.empty()) return std::numeric_limits<double>::quiet_NaN();
    std::sort(values.begin(), values.end());
    const std::size_t n = values.size();
    if (n % 2 == 1) return values[n / 2];
    return 0.5 * (values[(n / 2) - 1] + values[n / 2]);
}

// src/cpu/lsdtt_chi_network_analysis_cpu.cpp:1352-1369 (as of 2026-10-07)
double chi_network_median_of(std::vector<double> values)
{
    if (values.empty())
    {
        return 0.0;
    }

    std::sort(values.begin(), values.end());

    const std::size_t mid = values.size() / 2;

    if (values.size() % 2 == 0)
    {
        return 0.5 * (values[mid - 1] + values[mid]);
    }

    return values[mid];
}

// src/cpu/lsdtt_reference_field_ladder_cpu.cpp:34-43 (as of 2026-10-07)
double ladder_median_of_copy(std::vector<double> v)
{
    if (v.empty()) return 0.0;
    const std::size_t mid = v.size() / 2;
    std::nth_element(v.begin(), v.begin() + mid, v.end());
    const double upper = v[mid];
    if (v.size() % 2 != 0) return upper;
    const double lower = *std::max_element(v.begin(), v.begin() + mid);
    return 0.5 * (lower + upper);
}

// src/utils/lsdtt_signal1d.cpp:87-96 (as of 2026-10-07), the get_median lambda
double signal1d_get_median(std::vector<double> v)
{
    if (v.empty()) return 0.0;
    std::size_t mid_n = v.size() / 2;
    std::nth_element(v.begin(), v.begin() + mid_n, v.end());
    if (v.size() % 2 == 0) {
        auto rit = std::max_element(v.begin(), v.begin() + mid_n);
        return (*rit + v[mid_n]) / 2.0;
    }
    return v[mid_n];
}

// src/cpu/lsdtt_surface_texture_cpu.cpp:303-312 (as of 2026-10-07), inline in
// the window loop; the caller marks the window invalid when n < 3, so n >= 3.
double surface_texture_median(std::vector<float>& k2_diffs)
{
    const std::size_t mid = k2_diffs.size() / 2;
    std::nth_element(k2_diffs.begin(), k2_diffs.begin() + mid, k2_diffs.end());
    double med = static_cast<double>(k2_diffs[mid]);
    if (k2_diffs.size() % 2 == 0)
    {
        // Average with the largest value in the lower half.
        const float lo =
            *std::max_element(k2_diffs.begin(), k2_diffs.begin() + mid);
        med = (static_cast<double>(lo) + med) * 0.5;
    }
    return med;
}

// src/cpu/lsdtt_flat_surface_channel_classify_cpu.cpp:772-779 (as of
// 2026-10-07), inline; the caller skips empty h and then sums h IN ITS
// POST-SELECTION ORDER, so the order is part of the contract.
double flat_classify_median_height(std::vector<float>& h)
{
    double median_height;
    const std::size_t mid = h.size() / 2;
    std::nth_element(h.begin(), h.begin() + mid, h.end());
    median_height = h[mid];
    if (h.size() % 2 == 0)
    {
        const float upper = *std::max_element(h.begin(), h.begin() + mid);
        median_height = 0.5 * (median_height + upper);
    }
    return median_height;
}

// ---- upper middle value ----------------------------------------------------

// src/cpu/lsdtt_usgs_elevation.cpp:81-87 (as of 2026-10-07)
double usgs_median_of(std::vector<double> v)
{
    if (v.empty()) return 0.0;
    const std::size_t k = v.size() / 2;
    std::nth_element(v.begin(), v.begin() + k, v.end());
    return v[k];
}

// src/cpu/lsdtt_hydro_roughness_cpu.cpp:95-101 (as of 2026-10-07)
double hydro_roughness_median_of(std::vector<double> v)
{
    if (v.empty()) return 0.0;
    const std::size_t mid = v.size() / 2;
    std::nth_element(v.begin(), v.begin() + static_cast<std::ptrdiff_t>(mid), v.end());
    return v[mid];
}

// src/cpu/lsdtt_flat_surface_channel_classify_cpu.cpp:213-219 (as of 2026-10-07)
double flat_classify_median_of(std::vector<float>& v)
{
    if (v.empty()) return -1.0;
    const std::size_t mid = v.size() / 2;
    std::nth_element(v.begin(), v.begin() + mid, v.end());
    return static_cast<double>(v[mid]);
}

// src/cpu/lsdtt_landslide_placement_cpu.cpp:302-307 (as of 2026-10-07),
// LandslidePlacer::finalise_diagnostics
double landslide_median_host_fill_ratio(const std::vector<double>& fill_ratios_)
{
    if (fill_ratios_.empty()) { return 0.0; }
    std::vector<double> v = fill_ratios_;
    const std::size_t mid = v.size() / 2;
    std::nth_element(v.begin(), v.begin() + mid, v.end());
    return v[mid];
}

// src/cpu/lsdtt_lem_timestep_cpu.cpp:50-52 (as of 2026-10-07), the tail of
// median_interior_gradient
double lem_timestep_median(std::vector<double>& g)
{
    if (g.empty()) return 0.0;
    std::nth_element(g.begin(), g.begin() + g.size() / 2, g.end());
    return g[g.size() / 2];
}

// src/cpu/lsdtt_hydraulic_graph_cpu.cpp:1379-1381 (as of 2026-10-07); the
// caller returns early when there are no valid cells, so n >= 1.
double hydraulic_graph_median(std::vector<double>& abs_vals)
{
    auto mid_abs = abs_vals.begin() + abs_vals.size() / 2;
    std::nth_element(abs_vals.begin(), mid_abs, abs_vals.end());
    return *mid_abs;
}

// src/cpu/lsdtt_graphflood_optimised_cpu.cpp:866-871 (as of 2026-10-07);
// caller returns {0, 0} when empty.
double graphflood_median(std::vector<double>& convergence_abs_)
{
    const std::size_t middle = convergence_abs_.size() / 2;
    std::nth_element(convergence_abs_.begin(),
                     convergence_abs_.begin() + middle, convergence_abs_.end());
    return convergence_abs_[middle];
}

// src/cpu/lsdtt_reference_field_ladder_cpu.cpp:897-903 (as of 2026-10-07),
// median_slope_of lambda (ev_for at :889-895 takes the same median and
// compares it with a threshold).
double ladder_median_slope_of(const std::vector<double>& slopes)
{
    if (slopes.empty())
        return std::numeric_limits<double>::quiet_NaN();
    std::vector<double> s = slopes;
    const std::size_t mid = s.size() / 2;
    std::nth_element(s.begin(), s.begin() + static_cast<long>(mid), s.end());
    return s[mid];
}

// ---- interpolated quantiles ------------------------------------------------

// src/cpu/lsdtt_strahler_tokunaga_analysis_cpu.cpp:1374-1386 (as of 2026-10-07)
static double strahler_finite_percentile(std::vector<double> v, double p)
{
    v.erase(std::remove_if(v.begin(), v.end(),
            [](double x){ return !std::isfinite(x); }), v.end());
    if (v.empty()) return std::numeric_limits<double>::quiet_NaN();
    std::sort(v.begin(), v.end());
    const double rank = (p / 100.0) * static_cast<double>(v.size() - 1);
    const std::size_t lo = static_cast<std::size_t>(std::floor(rank));
    const std::size_t hi = static_cast<std::size_t>(std::ceil(rank));
    if (lo == hi) return v[lo];
    const double w = rank - static_cast<double>(lo);
    return v[lo] * (1.0 - w) + v[hi] * w;
}

// src/cpu/lsdtt_chi_disorder_cpu.cpp:596-637 (as of 2026-10-07); is_finite and
// quiet_nan were that file's float helpers.
float chi_disorder_percentile_of_finite(std::vector<float> values, double p)
{
    values.erase(
        std::remove_if(
            values.begin(),
            values.end(),
            [](float x) { return !std::isfinite(x); }
        ),
        values.end()
    );

    if (values.empty())
    {
        return std::numeric_limits<float>::quiet_NaN();
    }

    std::sort(values.begin(), values.end());

    if (p <= 0.0)
    {
        return values.front();
    }
    if (p >= 1.0)
    {
        return values.back();
    }

    const double position = p * static_cast<double>(values.size() - 1);
    const std::size_t lower = static_cast<std::size_t>(std::floor(position));
    const std::size_t upper = static_cast<std::size_t>(std::ceil(position));

    if (lower == upper)
    {
        return values[lower];
    }

    const double weight_upper = position - static_cast<double>(lower);
    const double weight_lower = 1.0 - weight_upper;

    return static_cast<float>(
        weight_lower * static_cast<double>(values[lower]) +
        weight_upper * static_cast<double>(values[upper])
    );
}

// src/cpu/lsdtt_raster_algebra.cpp:1571-1608 (as of 2026-10-07), after
// collect_valid (which throws on an empty valid set, so n >= 1).
float raster_algebra_percentile(std::vector<float> values, double p)
{
    const std::size_t n = values.size();

    if (n == 1) return values[0];

    const double rank = (p / 100.0) * static_cast<double>(n - 1);

    const std::size_t lower_index = static_cast<std::size_t>(std::floor(rank));
    const std::size_t upper_index = static_cast<std::size_t>(std::ceil(rank));

    std::nth_element(
        values.begin(),
        values.begin() + static_cast<std::ptrdiff_t>(upper_index),
        values.end()
    );
    const float upper_value = values[upper_index];

    if (lower_index == upper_index) return upper_value;

    std::nth_element(
        values.begin(),
        values.begin() + static_cast<std::ptrdiff_t>(lower_index),
        values.begin() + static_cast<std::ptrdiff_t>(upper_index) + 1
    );
    const float lower_value = values[lower_index];

    const double weight_upper = rank - static_cast<double>(lower_index);
    return static_cast<float>(
        (1.0 - weight_upper) * static_cast<double>(lower_value) +
        weight_upper         * static_cast<double>(upper_value)
    );
}

// ---- floor-index quantiles -------------------------------------------------

// src/cpu/lsdtt_raster_cpu.cpp:340-346 (as of 2026-10-07), the pct lambda; the
// caller returns early when slopes is empty and calls pct several times on the
// same vector.
double raster_cpu_pct(std::vector<float>& slopes, double p)
{
    const std::size_t k = static_cast<std::size_t>(
        p / 100.0 * static_cast<double>(slopes.size() - 1));
    std::nth_element(slopes.begin(), slopes.begin() + k, slopes.end());
    return static_cast<double>(slopes[k]);
}

// src/cpu/lsdtt_raster_flow_preprocessing_cpu.cpp:2399-2405 (as of
// 2026-10-07); the caller throws when vals is empty.
double flow_preprocessing_k(std::vector<float>& vals, double k_gradient_quantile)
{
    std::size_t idx = static_cast<std::size_t>(
        k_gradient_quantile * static_cast<double>(vals.size() - 1));
    idx = std::min(idx, vals.size() - 1);
    std::nth_element(vals.begin(), vals.begin() + static_cast<long>(idx),
                     vals.end());
    return static_cast<double>(vals[idx]);
}

// src/cpu/lsdtt_lem_cpu.cpp:676-683 (as of 2026-10-07); empty gives 0.
float lem_p999(std::vector<float>& local_changes)
{
    if (!local_changes.empty())
    {
        const std::size_t k = static_cast<std::size_t>(
            0.999 * static_cast<double>(local_changes.size() - 1));
        std::nth_element(local_changes.begin(),
                         local_changes.begin() + static_cast<std::ptrdiff_t>(k),
                         local_changes.end());
        return local_changes[k];
    }
    else return 0.0f;
}

// ============================================================================
// The migrated call sites, written exactly as they now appear in the sources
// ============================================================================

double new_strahler_finite_percentile(std::vector<double> v, double p)
{
    v.erase(std::remove_if(v.begin(), v.end(),
            [](double x){ return !std::isfinite(x); }), v.end());
    return compute_quantile_interpolated_in_place(v, p / 100.0);
}

float new_chi_disorder_percentile_of_finite(std::vector<float> values, double p)
{
    values.erase(
        std::remove_if(values.begin(), values.end(),
                       [](float x) { return !std::isfinite(x); }),
        values.end());
    if (values.empty()) return std::numeric_limits<float>::quiet_NaN();
    return static_cast<float>(
        compute_quantile_interpolated_in_place(values, std::clamp(p, 0.0, 1.0)));
}

float new_raster_algebra_percentile(std::vector<float> values, double p)
{
    return static_cast<float>(compute_quantile_interpolated_in_place(values, p / 100.0));
}

float new_lem_p999(std::vector<float>& local_changes)
{
    if (!local_changes.empty())
        return static_cast<float>(compute_quantile_floor_in_place(local_changes, 0.999));
    else return 0.0f;
}

// ============================================================================
// Harness
// ============================================================================

int g_failures = 0;

void check(const std::string& what, bool ok, long long n_cases, const std::string& detail = "")
{
    std::cout << (ok ? "PASS " : "FAIL ") << what << " (" << n_cases << " cases)";
    if (!ok && !detail.empty()) std::cout << " -- " << detail;
    std::cout << "\n";
    if (!ok) ++g_failures;
}

std::uint64_t bits(double x) { std::uint64_t b; std::memcpy(&b, &x, sizeof b); return b; }
std::uint32_t bits(float x)  { std::uint32_t b; std::memcpy(&b, &x, sizeof b); return b; }

template <class T>
bool same(T a, T b) { return (std::isnan(a) && std::isnan(b)) || bits(a) == bits(b); }

template <class T>
bool same_order(const std::vector<T>& a, const std::vector<T>& b)
{
    if (a.size() != b.size()) return false;
    for (std::size_t i = 0; i < a.size(); ++i) if (!same(a[i], b[i])) return false;
    return true;
}

// Generated inputs. signed_zeros=false keeps every zero +0.0; true mixes in -0.0.
template <class T>
std::vector<std::vector<T>> make_inputs(unsigned seed, bool signed_zeros)
{
    std::mt19937 rng(seed);
    std::vector<std::vector<T>> out;
    std::vector<std::size_t> sizes;
    for (std::size_t n = 0; n <= 41; ++n) sizes.push_back(n);
    sizes.push_back(1000);
    for (int rep = 0; rep < 12; ++rep)
        for (std::size_t n : sizes)
        {
            std::vector<T> v(n);
            const int mode = rep % 4;
            std::uniform_real_distribution<double> u(-50.0, 50.0);
            std::uniform_int_distribution<int> small(-3, 3);
            for (auto& x : v)
            {
                double d;
                if (mode == 0)      d = u(rng);                        // distinct
                else if (mode == 1) d = std::round(u(rng) / 10.0);     // heavy ties
                else if (mode == 2) d = small(rng) * 0.25;             // ties around zero
                else                d = std::exp(u(rng) / 10.0);       // positive, skewed
                if (d == 0.0) d = (signed_zeros && (rng() & 1u)) ? -0.0 : 0.0;
                x = static_cast<T>(d);
            }
            out.push_back(std::move(v));
        }
    return out;
}

// Compares a double-valued reference against the migrated site over inputs.
template <class T>
void compare_family(const std::string& name,
                    const std::vector<std::vector<T>>& inputs,
                    const std::function<double(std::vector<T>&)>& old_site,
                    const std::function<double(std::vector<T>&)>& new_site,
                    bool compare_order,
                    std::size_t min_n = 0,
                    bool signed_zero_may_differ = false)
{
    long long n = 0, bad = 0;
    std::string first;
    for (const auto& in : inputs)
    {
        if (in.size() < min_n) continue;
        std::vector<T> a = in, b = in;
        const double ra = old_site(a);
        const double rb = new_site(b);
        ++n;
        // D7 (Simon, 2026-10-07): where the selection algorithm changed, a zero
        // result may change sign when +0 and -0 are tied in the input.
        const bool value_ok = same(ra, rb) ||
            (signed_zero_may_differ && ra == 0.0 && rb == 0.0);
        const bool ok = value_ok && (!compare_order || same_order(a, b));
        if (!ok)
        {
            ++bad;
            if (first.empty())
                first = "n=" + std::to_string(in.size()) + " old=" + std::to_string(ra) +
                        " new=" + std::to_string(rb) +
                        (same(ra, rb) ? " (order differs)" : "");
        }
    }
    check(name, bad == 0, n, first);
}

} // namespace

int main()
{
    std::cout << "== in-place order statistics vs the removed private copies ==\n";

    using VD = std::vector<double>;
    using VF = std::vector<float>;

    for (bool signed_zeros : {false, true})
    {
        const std::string regime = signed_zeros ? " [mixed +-0]" : "";
        const auto ind = make_inputs<double>(signed_zeros ? 902u : 901u, signed_zeros);
        const auto inf = make_inputs<float>(signed_zeros ? 904u : 903u, signed_zeros);

        // ---- mean of the two middle values ----
        compare_family<double>("terrace median_copy_ == compute_median_in_place" + regime, ind,
            [](VD& v) { return terrace_median_copy_(v); },
            [](VD& v) { return compute_median_of_copy(v); }, false, 0, signed_zeros);
        compare_family<double>("paired_flowpaths median_of == compute_median_in_place" + regime, ind,
            [](VD& v) { return paired_median_of(v); },
            [](VD& v) { return compute_median_of_copy(v); }, false, 0, signed_zeros);
        compare_family<double>("sti median_in_place == compute_median_in_place" + regime, ind,
            [](VD& v) { return sti_median_in_place(v); },
            [](VD& v) { return compute_median_in_place(v); }, false, 0, signed_zeros);
        compare_family<double>("chi_network median_of == compute_median_in_place(.., 0)" + regime, ind,
            [](VD& v) { return chi_network_median_of(v); },
            [](VD& v) { return compute_median_of_copy(v, 0.0); }, false, 0, signed_zeros);
        compare_family<double>("ladder median_of_copy == compute_median_in_place(.., 0)" + regime, ind,
            [](VD& v) { return ladder_median_of_copy(v); },
            [](VD& v) { return compute_median_of_copy(v, 0.0); }, false);
        compare_family<double>("signal1d get_median == compute_median_in_place(.., 0)" + regime, ind,
            [](VD& v) { return signal1d_get_median(v); },
            [](VD& v) { return compute_median_of_copy(v, 0.0); }, false);
        compare_family<float>("surface_texture median == compute_median_in_place, incl. order" + regime, inf,
            [](VF& v) { return surface_texture_median(v); },
            [](VF& v) { return compute_median_in_place(v); }, true, 3);
        compare_family<float>("flat_classify median_height == compute_median_in_place, incl. order" + regime, inf,
            [](VF& v) { return flat_classify_median_height(v); },
            [](VF& v) { return compute_median_in_place(v); }, true, 1);

        // ---- upper middle value ----
        compare_family<double>("usgs median_of == compute_median_upper_in_place(.., 0)" + regime, ind,
            [](VD& v) { return usgs_median_of(v); },
            [](VD& v) { return compute_median_upper_of_copy(v, 0.0); }, false);
        compare_family<double>("hydro_roughness median_of == compute_median_upper_in_place(.., 0)" + regime, ind,
            [](VD& v) { return hydro_roughness_median_of(v); },
            [](VD& v) { return compute_median_upper_of_copy(v, 0.0); }, false);
        compare_family<float>("flat_classify median_of == compute_median_upper_in_place(.., -1), incl. order" + regime, inf,
            [](VF& v) { return flat_classify_median_of(v); },
            [](VF& v) { return compute_median_upper_in_place(v, -1.0); }, true);
        compare_family<double>("landslide median_host_fill_ratio == compute_median_upper_in_place(.., 0)" + regime, ind,
            [](VD& v) { return landslide_median_host_fill_ratio(v); },
            [](VD& v) { return compute_median_upper_of_copy(v, 0.0); }, false);
        compare_family<double>("lem_timestep median == compute_median_upper_in_place(.., 0), incl. order" + regime, ind,
            [](VD& v) { return lem_timestep_median(v); },
            [](VD& v) { return compute_median_upper_in_place(v, 0.0); }, true);
        compare_family<double>("hydraulic_graph median == compute_median_upper_in_place, incl. order" + regime, ind,
            [](VD& v) { return hydraulic_graph_median(v); },
            [](VD& v) { return compute_median_upper_in_place(v); }, true, 1);
        compare_family<double>("graphflood median == compute_median_upper_in_place, incl. order" + regime, ind,
            [](VD& v) { return graphflood_median(v); },
            [](VD& v) { return compute_median_upper_in_place(v); }, true, 1);
        compare_family<double>("ladder median_slope_of == compute_median_upper_in_place" + regime, ind,
            [](VD& v) { return ladder_median_slope_of(v); },
            [](VD& v) { return compute_median_upper_of_copy(v); }, false);

        // ---- interpolated quantiles ----
        for (double p : {0.0, 2.0, 10.0, 25.0, 33.3, 50.0, 75.0, 90.0, 97.5, 100.0})
        {
            const std::string tag = " p=" + std::to_string(p) + regime;
            // Non-finite values mixed into the filtering helpers' inputs.
            auto ind_nf = ind;
            auto inf_nf = inf;
            for (std::size_t i = 0; i < ind_nf.size(); i += 5)
                if (!ind_nf[i].empty())
                {
                    ind_nf[i][0] = std::numeric_limits<double>::infinity();
                    ind_nf[i].push_back(k_nan);
                }
            for (std::size_t i = 0; i < inf_nf.size(); i += 5)
                if (!inf_nf[i].empty())
                {
                    inf_nf[i][0] = -std::numeric_limits<float>::infinity();
                    inf_nf[i].push_back(std::numeric_limits<float>::quiet_NaN());
                }
            compare_family<double>("strahler finite_percentile" + tag, ind_nf,
                [p](VD& v) { return strahler_finite_percentile(v, p); },
                [p](VD& v) { return new_strahler_finite_percentile(v, p); }, false);
            compare_family<float>("chi_disorder percentile_of_finite" + tag, inf_nf,
                [p](VF& v) { return static_cast<double>(chi_disorder_percentile_of_finite(v, p / 100.0)); },
                [p](VF& v) { return static_cast<double>(new_chi_disorder_percentile_of_finite(v, p / 100.0)); }, false);
            compare_family<float>("raster_algebra raster_percentile" + tag, inf,
                [p](VF& v) { return static_cast<double>(raster_algebra_percentile(v, p)); },
                [p](VF& v) { return static_cast<double>(new_raster_algebra_percentile(v, p)); }, false, 1, signed_zeros);
        }

        // ---- floor-index quantiles ----
        compare_family<float>("raster_cpu pct sequence 50,90,99 (shared vector), incl. order" + regime, inf,
            [](VF& v) {
                const double a = raster_cpu_pct(v, 50.0), b = raster_cpu_pct(v, 90.0), c = raster_cpu_pct(v, 99.0);
                return a * 1e6 + b * 1e3 + c; },   // any change in any of the three shows
            [](VF& v) {
                const double a = compute_quantile_floor_in_place(v, 50.0 / 100.0);
                const double b = compute_quantile_floor_in_place(v, 90.0 / 100.0);
                const double c = compute_quantile_floor_in_place(v, 99.0 / 100.0);
                return a * 1e6 + b * 1e3 + c; }, true, 1);
        for (double q : {0.0, 0.5, 0.9, 0.95, 1.0})
            compare_family<float>("flow_preprocessing K quantile q=" + std::to_string(q) + ", incl. order" + regime, inf,
                [q](VF& v) { return flow_preprocessing_k(v, q); },
                [q](VF& v) { return compute_quantile_floor_in_place(v, q); }, true, 1);
        compare_family<float>("lem p99.9, incl. order" + regime, inf,
            [](VF& v) { return static_cast<double>(lem_p999(v)); },
            [](VF& v) { return static_cast<double>(new_lem_p999(v)); }, true);
    }

    std::cout << "RESULT: " << (g_failures == 0 ? "PASS" : "FAIL") << "\n";
    return g_failures == 0 ? 0 : 1;
}
