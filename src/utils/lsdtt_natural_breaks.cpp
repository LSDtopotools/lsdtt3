// lsdtt_natural_breaks.cpp

#include "lsdtt_natural_breaks.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace lsdtt {

namespace {

struct weighted_series
{
    std::vector<double> x;
    std::vector<double> w;
    std::vector<double> upper_edge;
};

weighted_series build_weighted_histogram(
    const std::vector<double>& values,
    const int max_bins,
    std::int64_t& n_values_used)
{
    if (max_bins <= 0)
    {
        throw std::runtime_error("compute_jenks_breaks: max_histogram_bins must be > 0");
    }

    double vmin = std::numeric_limits<double>::infinity();
    double vmax = -std::numeric_limits<double>::infinity();
    n_values_used = 0;
    for (const double v : values)
    {
        if (!std::isfinite(v)) continue;
        ++n_values_used;
        vmin = std::min(vmin, v);
        vmax = std::max(vmax, v);
    }

    if (n_values_used == 0)
    {
        throw std::runtime_error("compute_jenks_breaks: no finite values available");
    }

    weighted_series out;
    if (vmax <= vmin)
    {
        out.x.push_back(vmin);
        out.w.push_back(static_cast<double>(n_values_used));
        out.upper_edge.push_back(vmax);
        return out;
    }

    std::vector<double> counts(static_cast<std::size_t>(max_bins), 0.0);
    const double width = (vmax - vmin) / static_cast<double>(max_bins);
    for (const double v : values)
    {
        if (!std::isfinite(v)) continue;
        int idx = static_cast<int>(std::floor((v - vmin) / width));
        if (idx < 0) idx = 0;
        if (idx >= max_bins) idx = max_bins - 1;
        counts[static_cast<std::size_t>(idx)] += 1.0;
    }

    out.x.reserve(static_cast<std::size_t>(max_bins));
    out.w.reserve(static_cast<std::size_t>(max_bins));
    out.upper_edge.reserve(static_cast<std::size_t>(max_bins));

    for (int i = 0; i < max_bins; ++i)
    {
        const double c = counts[static_cast<std::size_t>(i)];
        if (c <= 0.0) continue;
        const double lo = vmin + static_cast<double>(i) * width;
        const double hi = (i == max_bins - 1) ? vmax : (lo + width);
        out.x.push_back(0.5 * (lo + hi));
        out.w.push_back(c);
        out.upper_edge.push_back(hi);
    }

    return out;
}

double sse_interval(
    const std::vector<double>& pref_w,
    const std::vector<double>& pref_wx,
    const std::vector<double>& pref_wx2,
    const int l,
    const int r)
{
    if (l > r) return 0.0;
    const double w = pref_w[r] - pref_w[l - 1];
    if (w <= 0.0) return 0.0;
    const double wx = pref_wx[r] - pref_wx[l - 1];
    const double wx2 = pref_wx2[r] - pref_wx2[l - 1];
    const double sse = wx2 - (wx * wx) / w;
    return sse > 0.0 ? sse : 0.0;
}

} // namespace

lsdtt_jenks_result compute_jenks_breaks(
    const std::vector<double>& values,
    const lsdtt_jenks_options& opts)
{
    if (opts.n_classes <= 0)
    {
        throw std::runtime_error("compute_jenks_breaks: n_classes must be > 0");
    }

    lsdtt_jenks_result result;
    result.requested_classes = opts.n_classes;

    std::int64_t n_values_used = 0;
    const weighted_series ws = build_weighted_histogram(
        values,
        opts.max_histogram_bins,
        n_values_used);

    result.n_values = n_values_used;
    result.n_histogram_bins_used = static_cast<std::int64_t>(ws.x.size());

    const int n = static_cast<int>(ws.x.size());
    if (n == 1)
    {
        result.effective_classes = 1;
        result.upper_breaks = {ws.upper_edge[0]};
        result.goodness_of_variance_fit = 1.0;
        return result;
    }

    const int k = std::min(opts.n_classes, n);
    result.effective_classes = k;

    std::vector<double> pref_w(static_cast<std::size_t>(n + 1), 0.0);
    std::vector<double> pref_wx(static_cast<std::size_t>(n + 1), 0.0);
    std::vector<double> pref_wx2(static_cast<std::size_t>(n + 1), 0.0);

    for (int i = 1; i <= n; ++i)
    {
        const double xi = ws.x[static_cast<std::size_t>(i - 1)];
        const double wi = ws.w[static_cast<std::size_t>(i - 1)];
        pref_w[static_cast<std::size_t>(i)] = pref_w[static_cast<std::size_t>(i - 1)] + wi;
        pref_wx[static_cast<std::size_t>(i)] = pref_wx[static_cast<std::size_t>(i - 1)] + wi * xi;
        pref_wx2[static_cast<std::size_t>(i)] = pref_wx2[static_cast<std::size_t>(i - 1)] + wi * xi * xi;
    }

    const double inf = std::numeric_limits<double>::infinity();
    std::vector<std::vector<double>> dp(
        static_cast<std::size_t>(k + 1),
        std::vector<double>(static_cast<std::size_t>(n + 1), inf));
    std::vector<std::vector<int>> back(
        static_cast<std::size_t>(k + 1),
        std::vector<int>(static_cast<std::size_t>(n + 1), -1));

    for (int i = 1; i <= n; ++i)
    {
        dp[1][i] = sse_interval(pref_w, pref_wx, pref_wx2, 1, i);
        back[1][i] = 0;
    }

    for (int c = 2; c <= k; ++c)
    {
        for (int i = c; i <= n; ++i)
        {
            double best = inf;
            int best_m = -1;
            for (int m = c - 1; m <= i - 1; ++m)
            {
                const double cand = dp[c - 1][m] + sse_interval(pref_w, pref_wx, pref_wx2, m + 1, i);
                if (cand < best)
                {
                    best = cand;
                    best_m = m;
                }
            }
            dp[c][i] = best;
            back[c][i] = best_m;
        }
    }

    std::vector<int> class_end(static_cast<std::size_t>(k), 0);
    int i = n;
    for (int c = k; c >= 1; --c)
    {
        class_end[static_cast<std::size_t>(c - 1)] = i;
        i = back[c][i];
        if (i < 0 && c > 1)
        {
            throw std::runtime_error("compute_jenks_breaks: dynamic-programming backtrace failed");
        }
    }

    result.upper_breaks.resize(static_cast<std::size_t>(k));
    for (int c = 0; c < k; ++c)
    {
        const int end_idx = class_end[static_cast<std::size_t>(c)] - 1;
        result.upper_breaks[static_cast<std::size_t>(c)] =
            ws.upper_edge[static_cast<std::size_t>(std::max(0, end_idx))];
    }
    result.upper_breaks.back() = ws.upper_edge.back();

    const double sdam = sse_interval(pref_w, pref_wx, pref_wx2, 1, n);
    const double sdcm = dp[k][n];
    if (sdam <= 0.0)
    {
        result.goodness_of_variance_fit = 1.0;
    }
    else
    {
        const double gvf = 1.0 - (sdcm / sdam);
        result.goodness_of_variance_fit = std::clamp(gvf, 0.0, 1.0);
    }

    return result;
}

int classify_jenks_value(
    const double value,
    const std::vector<double>& upper_breaks)
{
    if (!std::isfinite(value) || upper_breaks.empty()) return 0;
    for (std::size_t i = 0; i < upper_breaks.size(); ++i)
    {
        if (value <= upper_breaks[i])
        {
            return static_cast<int>(i + 1);
        }
    }
    return static_cast<int>(upper_breaks.size());
}

} // namespace lsdtt
