// lsdtt_slope_stability_dem_change_cpu.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html

#include "lsdtt_slope_stability_dem_change_cpu.hpp"

#include <algorithm>
#include <cmath>
#include <map>
#include <random>
#include <stdexcept>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "../utils/lsdtt_connected_components.hpp"
#include "../utils/lsdtt_statistics.hpp"
#include "../utils/lsdtt_thread_utils.hpp"
#include "lsdtt_flowgraph_cpu.hpp"
#include "lsdtt_raster_central_difference_cpu.hpp"
#include "lsdtt_slope_stability_catalogue_cpu.hpp"

namespace lsdtt
{

namespace
{

const double k_nmad_factor = 1.482602218505602;
const double k_pi = 3.14159265358979323846;

double nmad_of(const std::vector<double>& v)
{
    return k_nmad_factor * compute_mad(v);
}

double spherical(double h, double range)
{
    if (range <= 0.0) return 1.0;
    const double u = h / range;
    return u >= 1.0 ? 1.0 : 1.5 * u - 0.5 * u * u * u;
}

// sigma(S): linear between the slope bins' medians, constant beyond them.
double sigma_at(const dem_change_uncertainty& u, double slope)
{
    const std::vector<double>& x = u.slope_median;
    const std::vector<double>& y = u.slope_nmad_m;
    if (slope <= x.front()) return y.front();
    if (slope >= x.back()) return y.back();
    const std::size_t k = static_cast<std::size_t>(std::upper_bound(x.begin(), x.end(), slope) - x.begin());
    const double w = (slope - x[k - 1]) / (x[k] - x[k - 1]);
    return (1.0 - w) * y[k - 1] + w * y[k];
}

void estimate_variogram(dem_change_uncertainty& u, const lsdtt_raster_meta& meta,
                        const std::vector<std::size_t>& stable_cells, const std::vector<float>& z,
                        const std::vector<std::uint8_t>& is_stable, int n_pairs, std::uint64_t seed)
{
    const int nr = meta.n_rows, nc = meta.n_cols;
    const double max_lag = 0.5 * std::min(nr, nc);
    const int n_bins = 20;
    if (max_lag < 2.0 || stable_cells.size() < 100 || n_pairs < n_bins) return;

    std::mt19937_64 rng(seed);
    std::uniform_int_distribution<std::size_t> pick(0, stable_cells.size() - 1);
    std::uniform_real_distribution<double> unit(0.0, 1.0);
    const double log_max = std::log(max_lag);
    for (int b = 0; b < n_bins; ++b)
    {
        const double lo = std::exp(log_max * b / n_bins), hi = std::exp(log_max * (b + 1) / n_bins);
        const std::int64_t wanted = n_pairs / n_bins;
        std::int64_t count = 0;
        double sum = 0.0, sum_h = 0.0;
        for (std::int64_t attempt = 0; attempt < 4 * wanted && count < wanted; ++attempt)
        {
            const std::size_t i = stable_cells[pick(rng)];
            const double h = lo * std::exp(unit(rng) * std::log(hi / lo));
            const double angle = 2.0 * k_pi * unit(rng);
            const int r = static_cast<int>(i / nc), c = static_cast<int>(i % nc);
            const int r2 = r + static_cast<int>(std::lround(h * std::sin(angle)));
            const int c2 = c + static_cast<int>(std::lround(h * std::cos(angle)));
            if (r2 < 0 || r2 >= nr || c2 < 0 || c2 >= nc || (r2 == r && c2 == c)) continue;
            const std::size_t j = static_cast<std::size_t>(r2) * nc + c2;
            if (!is_stable[j]) continue;
            const double dz = static_cast<double>(z[i]) - static_cast<double>(z[j]);
            sum += 0.5 * dz * dz;
            sum_h += std::hypot(static_cast<double>(r2 - r), static_cast<double>(c2 - c));
            ++count;
        }
        if (count < 30) continue;
        u.lag_m.push_back(meta.cell_size * sum_h / static_cast<double>(count));
        u.semivariance.push_back(sum / static_cast<double>(count));
        u.lag_n_pairs.push_back(count);
    }
    if (u.lag_m.size() < 5) return;

    // Grid search on nugget and range; the sill has a closed form for each.
    double best = std::numeric_limits<double>::infinity();
    const double r_min = meta.cell_size, r_max = meta.cell_size * max_lag;
    for (int ir = 0; ir < 60; ++ir)
    {
        const double range = r_min * std::pow(r_max / r_min, ir / 59.0);
        for (int in = 0; in <= 19; ++in)
        {
            const double nugget = 0.05 * in;
            double sxy = 0.0, sxx = 0.0;
            for (std::size_t k = 0; k < u.lag_m.size(); ++k)
            {
                const double m = nugget + (1.0 - nugget) * spherical(u.lag_m[k], range);
                const double w = static_cast<double>(u.lag_n_pairs[k]);
                sxy += w * m * u.semivariance[k];
                sxx += w * m * m;
            }
            if (sxx <= 0.0) continue;
            const double sill = sxy / sxx;
            double cost = 0.0;
            for (std::size_t k = 0; k < u.lag_m.size(); ++k)
            {
                const double m = sill * (nugget + (1.0 - nugget) * spherical(u.lag_m[k], range));
                const double e = u.semivariance[k] - m;
                cost += static_cast<double>(u.lag_n_pairs[k]) * e * e;
            }
            if (cost < best)
            {
                best = cost;
                u.nugget_fraction = nugget;
                u.correlation_range_m = range;
                u.sill = sill;
            }
        }
    }
    u.variogram_fitted = true;
}

dem_change_uncertainty estimate_uncertainty(const lsdtt_raster_meta& meta, const std::vector<float>& dh,
                                            const std::vector<std::uint8_t>& assessable,
                                            const std::vector<double>& slope,
                                            const std::vector<std::uint8_t>* stable_mask,
                                            const dem_change_options& options)
{
    dem_change_uncertainty u;
    u.used_stable_mask = stable_mask != nullptr;
    const std::size_t n = dh.size();

    std::vector<double> first;
    first.reserve(n);
    for (std::size_t i = 0; i < n; ++i)
        if (assessable[i] && (!stable_mask || (*stable_mask)[i])) first.push_back(dh[i]);
    if (first.size() < 100)
        throw std::invalid_argument("detect_dem_change_landslides: fewer than 100 stable cells");

    const double raw_median = compute_median(first);
    const double raw_nmad = nmad_of(first);
    const double clip = 4.0 * std::max(raw_nmad, 1e-6);

    std::vector<std::size_t> cells;
    cells.reserve(first.size());
    std::vector<double> values;
    values.reserve(first.size());
    std::vector<std::uint8_t> is_stable(n, 0u);
    for (std::size_t i = 0; i < n; ++i)
    {
        if (!assessable[i] || (stable_mask && !(*stable_mask)[i])) continue;
        if (std::abs(static_cast<double>(dh[i]) - raw_median) > clip) { ++u.n_clipped; continue; }
        cells.push_back(i);
        values.push_back(dh[i]);
        is_stable[i] = 1u;
    }
    u.n_stable = static_cast<std::int64_t>(cells.size());
    if (cells.size() < 100)
        throw std::invalid_argument("detect_dem_change_landslides: fewer than 100 stable cells after clipping");

    u.bias_m = compute_median(values);
    u.nmad_m = nmad_of(values);
    double sum2 = 0.0, sum = 0.0;
    for (const double v : values)
    {
        sum += v;
        sum2 += v * v;
    }
    const double count = static_cast<double>(values.size());
    u.rmse_m = std::sqrt(sum2 / count);
    u.standard_deviation_m = std::sqrt(std::max(0.0, sum2 / count - (sum / count) * (sum / count)));
    u.p025_m = compute_percentile(values, 2.5);
    u.p975_m = compute_percentile(values, 97.5);

    // Quantile bins of slope.
    std::vector<std::pair<double, double>> by_slope;
    by_slope.reserve(cells.size());
    for (std::size_t k = 0; k < cells.size(); ++k) by_slope.emplace_back(slope[cells[k]], values[k]);
    // By slope only: sorting ties by value would put the lowest errors in one bin.
    std::stable_sort(by_slope.begin(), by_slope.end(),
                     [](const std::pair<double, double>& p, const std::pair<double, double>& q) {
                         return p.first < q.first;
                     });
    const int n_bins = std::max(1, options.n_slope_bins);
    std::size_t next = 0;
    for (int b = 0; b < n_bins && next < by_slope.size(); ++b)
    {
        const std::size_t lo = next;
        std::size_t hi = std::max(lo + 1, by_slope.size() * static_cast<std::size_t>(b + 1) / n_bins);
        // Equal slopes stay in one bin.
        while (hi < by_slope.size() && by_slope[hi].first == by_slope[hi - 1].first) ++hi;
        next = hi;
        std::vector<double> s, d;
        s.reserve(hi - lo);
        d.reserve(hi - lo);
        for (std::size_t k = lo; k < hi; ++k)
        {
            s.push_back(by_slope[k].first);
            d.push_back(by_slope[k].second);
        }
        const double median_slope = compute_median(s);
        u.slope_lower.push_back(s.front());
        u.slope_upper.push_back(s.back());
        u.slope_median.push_back(median_slope);
        u.slope_bias_m.push_back(compute_median(d));
        u.slope_nmad_m.push_back(std::max(nmad_of(d), 1e-6));
        u.slope_n.push_back(static_cast<std::int64_t>(d.size()));
    }

    std::vector<float> z(n, 0.0f);
    for (std::size_t k = 0; k < cells.size(); ++k)
        z[cells[k]] = static_cast<float>((values[k] - u.bias_m) / sigma_at(u, slope[cells[k]]));
    estimate_variogram(u, meta, cells, z, is_stable, options.n_variogram_pairs, options.random_seed);
    return u;
}

struct detection
{
    std::vector<std::int32_t> label;
    std::vector<dem_change_patch> patches;
};

// Everything the detection needs that does not change between runs.
struct detection_context
{
    const lsdtt_raster_meta* meta = nullptr;
    const std::vector<std::uint8_t>* assessable = nullptr;
    const std::vector<float>* sigma = nullptr;
    const std::vector<double>* slope = nullptr;
    const dem_change_uncertainty* uncertainty = nullptr;
    std::vector<std::pair<int, int>> disc;
    bool want_loss = true, want_gain = false;
    double patch_significance = 1.96;
    double min_mean_change = 0.0;
    int min_patch_cells = 5;
    bool full_records = true;
};

double mean_correlation_cached(std::map<std::int64_t, double>& cache, std::int64_t n_cells, double cell_area,
                               double range)
{
    const auto it = cache.find(n_cells);
    if (it != cache.end()) return it->second;
    const double r = std::sqrt(static_cast<double>(n_cells) * cell_area / k_pi);
    const double c = spherical_disc_mean_correlation(r, range);
    cache[n_cells] = c;
    return c;
}

void detect_one_sign(const detection_context& ctx, const std::vector<float>& corrected, int sign, double k_low,
                     double k_high, detection& out, std::map<std::int64_t, double>& correlation_cache)
{
    const lsdtt_raster_meta& meta = *ctx.meta;
    const int nr = meta.n_rows, nc = meta.n_cols;
    const std::size_t n = corrected.size();
    const std::vector<std::uint8_t>& ok = *ctx.assessable;
    const std::vector<float>& sigma = *ctx.sigma;

    std::vector<std::uint8_t> candidate(n, 0u);
    for (std::size_t i = 0; i < n; ++i)
        if (ok[i] && sign * corrected[i] >= k_low * sigma[i]) candidate[i] = 1u;

    std::vector<std::int32_t> component;
    std::vector<std::int64_t> size;
    const std::int32_t n_components = label_connected_components(
        nr, nc, 8, [&](std::size_t i) { return candidate[i] != 0u; }, component, size);
    if (n_components == 0) return;

    // Hysteresis.
    std::vector<std::uint8_t> has_seed(static_cast<std::size_t>(n_components) + 1, 0u);
    for (std::size_t i = 0; i < n; ++i)
        if (component[i] > 0 && sign * corrected[i] >= k_high * sigma[i])
            has_seed[static_cast<std::size_t>(component[i])] = 1u;

    // Aggregation: union of kept cells whose centres are within the search disc.
    std::vector<std::uint8_t> kept(n, 0u);
    std::vector<std::int32_t> index(n, -1);
    std::vector<std::size_t> kept_cells;
    for (std::size_t i = 0; i < n; ++i)
    {
        if (component[i] == 0 || !has_seed[static_cast<std::size_t>(component[i])]) continue;
        kept[i] = 1u;
        index[i] = static_cast<std::int32_t>(kept_cells.size());
        kept_cells.push_back(i);
    }
    if (kept_cells.empty()) return;
    std::vector<std::int32_t> parent(kept_cells.size());
    for (std::size_t k = 0; k < parent.size(); ++k) parent[k] = static_cast<std::int32_t>(k);
    const auto find = [&](std::int32_t x) {
        while (parent[static_cast<std::size_t>(x)] != x)
        {
            parent[static_cast<std::size_t>(x)] = parent[static_cast<std::size_t>(parent[static_cast<std::size_t>(x)])];
            x = parent[static_cast<std::size_t>(x)];
        }
        return x;
    };
    for (std::size_t k = 0; k < kept_cells.size(); ++k)
    {
        const int r = static_cast<int>(kept_cells[k] / nc), c = static_cast<int>(kept_cells[k] % nc);
        for (const auto& d : ctx.disc)
        {
            const int r2 = r + d.first, c2 = c + d.second;
            if (r2 < 0 || r2 >= nr || c2 < 0 || c2 >= nc) continue;
            const std::int32_t other = index[static_cast<std::size_t>(r2) * nc + c2];
            if (other < 0) continue;
            const std::int32_t ra = find(static_cast<std::int32_t>(k)), rb = find(other);
            if (ra != rb) parent[static_cast<std::size_t>(std::max(ra, rb))] = std::min(ra, rb);
        }
    }
    // Groups numbered 1.. in row-major order of their first cell.
    std::vector<std::int32_t> group(n, 0), group_of_root(kept_cells.size(), 0);
    std::int32_t n_groups = 0;
    for (std::size_t k = 0; k < kept_cells.size(); ++k)
    {
        const std::size_t root = static_cast<std::size_t>(find(static_cast<std::int32_t>(k)));
        if (group_of_root[root] == 0) group_of_root[root] = ++n_groups;
        group[kept_cells[k]] = group_of_root[root];
    }

    const double cell_area = meta.cell_size * meta.cell_size;
    const dem_change_uncertainty& u = *ctx.uncertainty;
    const std::size_t m = static_cast<std::size_t>(n_groups) + 1;
    std::vector<std::int64_t> count(m, 0), edges(m, 0);
    std::vector<double> sum_change(m, 0.0), sum_var(m, 0.0), extreme(m, 0.0), sum_slope(m, 0.0);
    std::vector<cell_shape_moments> moments;
    if (ctx.full_records) moments.resize(m);

    for (std::size_t i = 0; i < n; ++i)
    {
        if (!kept[i]) continue;
        const std::size_t g = static_cast<std::size_t>(group[i]);
        const double v = corrected[i];
        ++count[g];
        sum_change[g] += v;
        sum_var[g] += static_cast<double>(sigma[i]) * sigma[i];
        if (std::abs(v) > std::abs(extreme[g])) extreme[g] = v;
        if (ctx.full_records)
        {
            const int r = static_cast<int>(i / nc), c = static_cast<int>(i % nc);
            sum_slope[g] += (*ctx.slope)[i];
            moments[g].add(r, c);
            // A neighbour counts as inside only if it changed and is in this patch.
            const auto outside = [&](int r2, int c2) {
                if (r2 < 0 || r2 >= nr || c2 < 0 || c2 >= nc) return true;
                const std::size_t j = static_cast<std::size_t>(r2) * nc + c2;
                return !kept[j] || static_cast<std::size_t>(group[j]) != g;
            };
            edges[g] += (outside(r - 1, c) ? 1 : 0) + (outside(r + 1, c) ? 1 : 0) +
                        (outside(r, c - 1) ? 1 : 0) + (outside(r, c + 1) ? 1 : 0);
        }
    }

    std::vector<std::int32_t> new_id(m, 0);
    for (std::size_t g = 1; g < m; ++g)
    {
        if (count[g] < ctx.min_patch_cells) continue;
        const double c = static_cast<double>(count[g]);
        const double area = c * cell_area;
        const double correlated = u.variogram_fitted
                                      ? mean_correlation_cached(correlation_cache, count[g], cell_area, u.correlation_range_m)
                                      : 0.0;
        const double nugget = u.variogram_fitted ? u.nugget_fraction : 1.0;
        const double sigma_mean = std::sqrt((sum_var[g] / c) * (nugget * cell_area / area + (1.0 - nugget) * correlated));
        const double mean = sum_change[g] / c;
        const double significance = std::abs(mean) / std::max(sigma_mean, 1e-12);
        if (significance < ctx.patch_significance) continue;
        if (std::abs(mean) < ctx.min_mean_change) continue;

        dem_change_patch p;
        p.id = static_cast<std::int32_t>(out.patches.size()) + 1;
        p.sign = sign;
        p.n_cells = count[g];
        p.area_m2 = area;
        p.mean_change_m = mean;
        p.extreme_change_m = extreme[g];
        p.mean_change_sigma_m = sigma_mean;
        p.significance = significance;
        p.volume_m3 = mean * area;
        p.volume_sigma_m3 = sigma_mean * area;
        if (ctx.full_records)
        {
            p.perimeter_m = static_cast<double>(edges[g]) * meta.cell_size;
            p.area_uncertainty_m2 = 0.5 * p.perimeter_m * meta.cell_size;
            const cell_shape shape = shape_from_moments(moments[g], meta);
            p.centroid_x = shape.centroid_x;
            p.centroid_y = shape.centroid_y;
            p.length_m = shape.length_m;
            p.width_m = shape.width_m;
            p.mean_slope = sum_slope[g] / c;
        }
        out.patches.push_back(p);
        new_id[g] = p.id;
    }
    for (std::size_t i = 0; i < n; ++i)
        if (kept[i] && new_id[static_cast<std::size_t>(group[i])] > 0)
            out.label[i] = new_id[static_cast<std::size_t>(group[i])];
}

detection run_detection(const detection_context& ctx, const std::vector<float>& corrected, double k_low, double k_high,
                        std::map<std::int64_t, double>& correlation_cache)
{
    detection out;
    out.label.assign(corrected.size(), 0);
    if (ctx.want_loss) detect_one_sign(ctx, corrected, -1, k_low, k_high, out, correlation_cache);
    if (ctx.want_gain) detect_one_sign(ctx, corrected, +1, k_low, k_high, out, correlation_cache);
    return out;
}

// One pass of a running mean along a line, clamped at the ends.
void box_line(std::vector<double>& line, std::vector<double>& scratch, int half_width)
{
    const int n = static_cast<int>(line.size());
    const auto at = [&](int k) { return line[static_cast<std::size_t>(std::max(0, std::min(n - 1, k)))]; };
    double sum = 0.0;
    for (int k = -half_width; k <= half_width; ++k) sum += at(k);
    const double width = 2.0 * half_width + 1.0;
    for (int k = 0; k < n; ++k)
    {
        scratch[static_cast<std::size_t>(k)] = sum / width;
        sum += at(k + half_width + 1) - at(k - half_width);
    }
    line.swap(scratch);
}

// Unit-variance noise: a white part with variance `nugget`, and a part smoothed
// by three box passes, close to a Gaussian kernel. A spherical correlation falls
// to 0.5 at 0.347 L and a Gaussian-smoothed field at 1.665 sigma_k, so
// sigma_k = 0.208 L.
std::vector<float> simulate_noise(int nr, int nc, double cell, double nugget, double range, std::uint64_t seed)
{
    const std::size_t n = static_cast<std::size_t>(nr) * nc;
    std::mt19937_64 rng(seed);
    std::normal_distribution<double> gauss(0.0, 1.0);
    std::vector<double> smooth(n), white(n);
    for (std::size_t i = 0; i < n; ++i)
    {
        smooth[i] = gauss(rng);
        white[i] = gauss(rng);
    }
    const double sigma_cells = 0.208 * range / cell;
    const int half = static_cast<int>(std::lround((std::sqrt(4.0 * sigma_cells * sigma_cells + 1.0) - 1.0) / 2.0));
    if (half > 0)
    {
        std::vector<double> line(static_cast<std::size_t>(nc)), scratch(static_cast<std::size_t>(nc));
        for (int r = 0; r < nr; ++r)
        {
            std::copy_n(smooth.begin() + static_cast<std::ptrdiff_t>(r) * nc, nc, line.begin());
            for (int pass = 0; pass < 3; ++pass) box_line(line, scratch, half);
            std::copy_n(line.begin(), nc, smooth.begin() + static_cast<std::ptrdiff_t>(r) * nc);
        }
        line.resize(static_cast<std::size_t>(nr));
        scratch.resize(static_cast<std::size_t>(nr));
        for (int c = 0; c < nc; ++c)
        {
            for (int r = 0; r < nr; ++r) line[static_cast<std::size_t>(r)] = smooth[static_cast<std::size_t>(r) * nc + c];
            for (int pass = 0; pass < 3; ++pass) box_line(line, scratch, half);
            for (int r = 0; r < nr; ++r) smooth[static_cast<std::size_t>(r) * nc + c] = line[static_cast<std::size_t>(r)];
        }
    }
    double sum = 0.0, sum2 = 0.0;
    for (const double v : smooth)
    {
        sum += v;
        sum2 += v * v;
    }
    const double mean = sum / static_cast<double>(n);
    const double sd = std::sqrt(std::max(1e-12, sum2 / static_cast<double>(n) - mean * mean));
    std::vector<float> out(n);
    const double a = std::sqrt(1.0 - nugget), b = std::sqrt(nugget);
    for (std::size_t i = 0; i < n; ++i)
        out[i] = static_cast<float>(a * (smooth[i] - mean) / sd + b * white[i]);
    return out;
}

} // namespace

double spherical_disc_mean_correlation(double r, double L)
{
    if (!(r > 0.0)) return 1.0;
    if (!(L > 0.0)) return 0.0;
    // Distance h between two random points of a disc:
    //   p(h) = (4 h / (pi r^2)) [acos(q) - q sqrt(1 - q^2)], q = h / (2 r).
    const int n = 400;
    const double h_max = std::min(2.0 * r, L);
    const double dh = h_max / n;
    double total = 0.0;
    for (int k = 0; k < n; ++k)
    {
        const double h = (k + 0.5) * dh;
        const double q = h / (2.0 * r);
        const double p = 4.0 * h / (k_pi * r * r) * (std::acos(q) - q * std::sqrt(1.0 - q * q));
        total += (1.0 - spherical(h, L)) * p * dh;
    }
    return std::min(1.0, std::max(0.0, total));
}

std::vector<double> dem_change_slope(const lsdtt_raster_cpu& dem)
{
    const lsdtt_raster_meta& meta = dem.get_meta();
    const int nr = meta.n_rows, nc = meta.n_cols;
    const std::vector<float>& z = dem.get_data();
    const std::vector<std::uint8_t>& mask = dem.get_mask();
    std::vector<double> out(z.size(), std::numeric_limits<double>::quiet_NaN());
    for (int r = 1; r < nr - 1; ++r)
        for (int c = 1; c < nc - 1; ++c)
        {
            const std::size_t i = static_cast<std::size_t>(r) * nc + c;
            double gx = 0.0, gy = 0.0;
            if (mask[i] && cardinal_central_gradient(z, mask, nc, i, meta.cell_size, gx, gy))
                out[i] = gradient_magnitude(gx, gy);
        }
    return out;
}

inventory_dod_match match_inventory_to_dem_change(
    const dem_change_result&                      result,
    const lsdtt_raster_meta&                      meta,
    const std::vector<std::vector<std::int64_t>>& polygon_cells)
{
    const std::size_t n = static_cast<std::size_t>(meta.n_rows) * meta.n_cols;
    if (result.label.size() != n || result.corrected_m.size() != n)
        throw std::invalid_argument("match_inventory_to_dem_change: the result does not match the grid");
    const double a = meta.cell_size * meta.cell_size;
    const float nodata = static_cast<float>(meta.no_data_value);
    // Patch id -> position, and sign.
    std::vector<std::size_t> pos_of(result.patches.size() + 1, 0);
    inventory_dod_match out;
    out.patches.resize(result.patches.size());
    for (std::size_t k = 0; k < result.patches.size(); ++k)
    {
        const auto& p = result.patches[k];
        if (p.id < 1 || static_cast<std::size_t>(p.id) >= pos_of.size())
            throw std::invalid_argument("match_inventory_to_dem_change: unexpected patch id");
        pos_of[static_cast<std::size_t>(p.id)] = k;
        out.patches[k].patch_id = p.id;
        out.patches[k].sign = p.sign;
        if (p.sign < 0) ++out.n_loss_patches;
    }
    const auto patch_at = [&](std::size_t i) -> const dem_change_patch* {
        const std::int32_t id = result.label[i];
        if (id < 1 || static_cast<std::size_t>(id) >= pos_of.size()) return nullptr;
        return &result.patches[pos_of[static_cast<std::size_t>(id)]];
    };
    for (std::size_t i = 0; i < n; ++i)
        if (const dem_change_patch* p = patch_at(i))
        {
            ++out.patches[pos_of[static_cast<std::size_t>(p->id)]].n_cells;
            if (p->sign < 0) out.detected_loss_area_m2 += a;
        }

    std::vector<std::uint8_t> in_union(n, 0u);
    out.polygons.resize(polygon_cells.size());
    for (std::size_t q = 0; q < polygon_cells.size(); ++q)
    {
        inventory_dod_polygon_match& m = out.polygons[q];
        for (const std::int64_t ci : polygon_cells[q])
        {
            if (ci < 0 || static_cast<std::size_t>(ci) >= n)
                throw std::invalid_argument("match_inventory_to_dem_change: a polygon cell is off the grid");
            const std::size_t i = static_cast<std::size_t>(ci);
            ++m.n_cells;
            in_union[i] = 1u;
            const float c = result.corrected_m[i];
            const bool valid = c != nodata && std::isfinite(c);
            if (valid)
            {
                ++m.n_valid_cells;
                m.dod_net_change_m3 += c * a;
                if (c < 0) m.dod_loss_volume_m3 -= c * a; else m.dod_gain_volume_m3 += c * a;
            }
            if (const dem_change_patch* p = patch_at(i))
            {
                const double mag = valid ? std::abs(static_cast<double>(c)) * a : 0.0;
                if (p->sign < 0) { m.detected_loss_area_m2 += a; m.detected_loss_volume_m3 += mag; m.loss_patch_ids.push_back(p->id); }
                else { m.detected_gain_area_m2 += a; m.detected_gain_volume_m3 += mag; m.gain_patch_ids.push_back(p->id); }
                auto& pl = out.patches[pos_of[static_cast<std::size_t>(p->id)]].polygons;
                if (pl.empty() || pl.back() != static_cast<int>(q)) pl.push_back(static_cast<int>(q));
            }
        }
        m.area_m2 = static_cast<double>(m.n_cells) * a;
        for (auto* ids : {&m.loss_patch_ids, &m.gain_patch_ids})
        {
            std::sort(ids->begin(), ids->end());
            ids->erase(std::unique(ids->begin(), ids->end()), ids->end());
        }
        m.matched = !m.loss_patch_ids.empty();
        if (m.matched) ++out.n_polygons_matched;
    }
    for (std::size_t i = 0; i < n; ++i)
    {
        if (!in_union[i]) continue;
        out.polygon_union_area_m2 += a;
        if (const dem_change_patch* p = patch_at(i))
        {
            ++out.patches[pos_of[static_cast<std::size_t>(p->id)]].n_cells_in_polygons;
            if (p->sign < 0) out.detected_loss_in_polygons_m2 += a;
        }
    }
    for (const auto& p : out.patches)
        if (p.sign < 0 && p.n_cells_in_polygons > 0) ++out.n_loss_patches_matched;
    if (out.polygon_union_area_m2 > 0)
        out.mapped_area_detected_fraction = out.detected_loss_in_polygons_m2 / out.polygon_union_area_m2;
    if (out.detected_loss_area_m2 > 0)
        out.detected_area_mapped_fraction = out.detected_loss_in_polygons_m2 / out.detected_loss_area_m2;
    return out;
}

std::vector<float> simulate_dem_error_noise(
    int nr, int nc, double cell_size, double nugget_fraction, double correlation_range_m, std::uint64_t seed)
{
    if (nr < 1 || nc < 1 || !(cell_size > 0.0) || !(nugget_fraction >= 0.0 && nugget_fraction <= 1.0) ||
        !(correlation_range_m >= 0.0))
        throw std::invalid_argument("simulate_dem_error_noise: an argument is out of range");
    return simulate_noise(nr, nc, cell_size, nugget_fraction, correlation_range_m, seed);
}

dem_change_error_model estimate_dem_change_error_model(
    const lsdtt_raster_cpu&          difference,
    const std::vector<double>&       slope,
    const std::vector<std::uint8_t>* stable_mask,
    const dem_change_options&        options,
    const std::vector<std::uint8_t>* exclusion_mask)
{
    if (options.min_slope < 0.0 || options.n_slope_bins < 1 || options.n_variogram_pairs < 0)
        throw std::invalid_argument("estimate_dem_change_error_model: an option is out of range");
    const lsdtt_raster_meta& meta = difference.get_meta();
    meta.check_projected_metres();
    const std::size_t n = static_cast<std::size_t>(meta.n_rows) * meta.n_cols;
    if (slope.size() != n || (stable_mask && stable_mask->size() != n) || (exclusion_mask && exclusion_mask->size() != n))
        throw std::invalid_argument(
            "estimate_dem_change_error_model: slope, stable_mask and exclusion_mask need one value per cell");
    const std::vector<float>& dh = difference.get_data();
    const std::vector<std::uint8_t>& mask = difference.get_mask();
    std::vector<std::uint8_t> assessable(n, 0u);
    for (std::size_t i = 0; i < n; ++i)
    {
        if (!mask[i] || !std::isfinite(dh[i])) continue;
        if (exclusion_mask && (*exclusion_mask)[i]) continue;
        if (std::isfinite(slope[i]) && slope[i] >= options.min_slope) assessable[i] = 1u;
    }
    dem_change_error_model out;
    out.uncertainty = estimate_uncertainty(meta, dh, assessable, slope, stable_mask, options);
    out.sigma_m.assign(n, std::numeric_limits<float>::quiet_NaN());
    for (std::size_t i = 0; i < n; ++i)
        if (std::isfinite(slope[i])) out.sigma_m[i] = static_cast<float>(sigma_at(out.uncertainty, slope[i]));
    return out;
}

dem_change_result detect_dem_change_landslides(
    const lsdtt_raster_cpu&          difference,
    const std::vector<double>&       slope,
    const std::vector<std::uint8_t>* stable_mask,
    const dem_change_options&        options,
    const std::vector<std::uint8_t>* exclusion_mask)
{
    if (options.change_sign != "loss" && options.change_sign != "gain" && options.change_sign != "both")
        throw std::invalid_argument("detect_dem_change_landslides: change_sign must be loss, gain or both");
    if (!(options.detection_sigma > 0.0) || options.seed_sigma < options.detection_sigma)
        throw std::invalid_argument("detect_dem_change_landslides: need 0 < detection_sigma <= seed_sigma");
    if (options.patch_significance < 0.0 || options.min_patch_cells < 1 || options.aggregation_radius_pixels < 0 ||
        options.min_slope < 0.0 || options.n_slope_bins < 1 || options.n_null_simulations < 0 ||
        options.n_variogram_pairs < 0 || !(options.min_patch_mean_change_m >= 0.0))
        throw std::invalid_argument("detect_dem_change_landslides: an option is out of range");

    const lsdtt_raster_meta& meta = difference.get_meta();
    meta.check_projected_metres();
    const int nr = meta.n_rows, nc = meta.n_cols;
    const std::size_t n = static_cast<std::size_t>(nr) * nc;
    if (slope.size() != n || (stable_mask && stable_mask->size() != n) || (exclusion_mask && exclusion_mask->size() != n))
        throw std::invalid_argument(
            "detect_dem_change_landslides: slope, stable_mask and exclusion_mask need one value per cell");

    const std::vector<float>& dh = difference.get_data();
    const std::vector<std::uint8_t>& mask = difference.get_mask();
    std::vector<std::uint8_t> assessable(n, 0u);
    dem_change_result out;
    for (std::size_t i = 0; i < n; ++i)
    {
        if (!mask[i] || !std::isfinite(dh[i])) continue;
        if (exclusion_mask && (*exclusion_mask)[i]) { ++out.n_excluded_cells; continue; }
        if (std::isfinite(slope[i]) && slope[i] >= options.min_slope)
        {
            assessable[i] = 1u;
            ++out.n_valid_cells;
        }
    }
    const double cell_area = meta.cell_size * meta.cell_size;
    out.valid_area_m2 = static_cast<double>(out.n_valid_cells) * cell_area;

    out.uncertainty = estimate_uncertainty(meta, dh, assessable, slope, stable_mask, options);
    const dem_change_uncertainty& u = out.uncertainty;

    const float nodata = static_cast<float>(meta.no_data_value);
    out.corrected_m.assign(n, nodata);
    out.sigma_m.assign(n, nodata);
    out.z.assign(n, nodata);
    std::vector<float> corrected(n, 0.0f), sigma(n, 0.0f);
    for (std::size_t i = 0; i < n; ++i)
    {
        if (!mask[i] || !std::isfinite(dh[i])) continue;
        const float value = static_cast<float>(static_cast<double>(dh[i]) - u.bias_m);
        out.corrected_m[i] = value;
        corrected[i] = value;
        if (!std::isfinite(slope[i])) continue;
        const float s = static_cast<float>(sigma_at(u, slope[i]));
        sigma[i] = s;
        out.sigma_m[i] = s;
        out.z[i] = value / s;
    }

    detection_context ctx;
    ctx.meta = &meta;
    ctx.assessable = &assessable;
    ctx.sigma = &sigma;
    ctx.slope = &slope;
    ctx.uncertainty = &u;
    ctx.want_loss = options.change_sign != "gain";
    ctx.want_gain = options.change_sign != "loss";
    ctx.patch_significance = options.patch_significance;
    ctx.min_mean_change = options.min_patch_mean_change_m;
    ctx.min_patch_cells = options.min_patch_cells;
    const int radius = options.aggregation_radius_pixels;
    // Half of the neighbourhood, since a union is symmetric. Radius 0 is the
    // eight neighbours; radius g joins cells with up to g cells between them.
    const int reach = radius + 1;
    for (int r = 0; r <= reach; ++r)
        for (int c = -reach; c <= reach; ++c)
            if ((r > 0 || c > 0) && r * r + c * c <= reach * reach + 1) ctx.disc.emplace_back(r, c);

    std::map<std::int64_t, double> correlation_cache;
    detection found = run_detection(ctx, corrected, options.detection_sigma, options.seed_sigma, correlation_cache);
    out.patches = std::move(found.patches);
    out.label = std::move(found.label);
    out.change_class.assign(n, 0);
    for (std::size_t i = 0; i < n; ++i)
        if (out.label[i] > 0)
            out.change_class[i] = static_cast<std::int8_t>(out.patches[static_cast<std::size_t>(out.label[i] - 1)].sign);

    double loss_var = 0.0, gain_var = 0.0, loss_a2 = 0.0, gain_a2 = 0.0;
    for (const dem_change_patch& p : out.patches)
    {
        const double v2 = p.volume_sigma_m3 * p.volume_sigma_m3, a2 = p.area_uncertainty_m2 * p.area_uncertainty_m2;
        if (p.sign < 0) { out.loss_area_m2 += p.area_m2; out.loss_volume_m3 += p.volume_m3; loss_var += v2; loss_a2 += a2; }
        else            { out.gain_area_m2 += p.area_m2; out.gain_volume_m3 += p.volume_m3; gain_var += v2; gain_a2 += a2; }
    }
    out.loss_volume_sigma_m3 = std::sqrt(loss_var);
    out.gain_volume_sigma_m3 = std::sqrt(gain_var);
    out.loss_area_uncertainty_m2 = std::sqrt(loss_a2);
    out.gain_area_uncertainty_m2 = std::sqrt(gain_a2);

    // Sensitivity to the threshold, with the seed scaled in proportion.
    detection_context light = ctx;
    light.full_records = false;
    const double seed_ratio = options.seed_sigma / options.detection_sigma;
    for (const double k : {1.5, options.detection_sigma, 2.5, 3.0, 4.0})
    {
        bool repeated = false;
        for (const auto& row : out.sensitivity) repeated = repeated || row.detection_sigma == k;
        if (repeated) continue;
        const detection d = run_detection(light, corrected, k, k * seed_ratio, correlation_cache);
        dem_change_sensitivity_row row;
        row.detection_sigma = k;
        row.seed_sigma = k * seed_ratio;
        row.n_patches = static_cast<std::int64_t>(d.patches.size());
        for (const auto& p : d.patches)
        {
            row.area_m2 += p.area_m2;
            row.volume_m3 += p.volume_m3;
        }
        row.coverage_fraction = out.valid_area_m2 > 0.0 ? row.area_m2 / out.valid_area_m2 : 0.0;
        out.sensitivity.push_back(row);
    }
    std::sort(out.sensitivity.begin(), out.sensitivity.end(),
              [](const dem_change_sensitivity_row& a, const dem_change_sensitivity_row& b) {
                  return a.detection_sigma < b.detection_sigma;
              });

    // False detections on simulated stable terrain.
    if (options.n_null_simulations > 0)
    {
        const int n_sim = options.n_null_simulations;
        const double nugget = u.variogram_fitted ? u.nugget_fraction : 1.0;
        const double range = u.variogram_fitted ? u.correlation_range_m : 0.0;
        std::vector<std::vector<double>> sim_areas(static_cast<std::size_t>(n_sim));
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) num_threads(resolve_thread_count(options.n_threads))
#endif
        for (int s = 0; s < n_sim; ++s)
        {
            const std::vector<float> noise =
                simulate_noise(nr, nc, meta.cell_size, nugget, range, options.random_seed + 1000u + static_cast<std::uint64_t>(s));
            std::vector<float> simulated(n, 0.0f);
            for (std::size_t i = 0; i < n; ++i)
                if (assessable[i]) simulated[i] = noise[i] * sigma[i];
            std::map<std::int64_t, double> local_cache;
            const detection d = run_detection(light, simulated, options.detection_sigma, options.seed_sigma, local_cache);
            std::vector<double>& areas = sim_areas[static_cast<std::size_t>(s)];
            areas.reserve(d.patches.size());
            for (const auto& p : d.patches) areas.push_back(p.area_m2);
        }

        dem_change_null_summary& null = out.null_simulation;
        null.n_simulations = n_sim;
        std::vector<double> counts, totals;
        double largest = cell_area;
        for (const auto& areas : sim_areas)
        {
            counts.push_back(static_cast<double>(areas.size()));
            double total = 0.0;
            for (const double a : areas)
            {
                total += a;
                largest = std::max(largest, a);
            }
            totals.push_back(total);
        }
        for (const auto& p : out.patches) largest = std::max(largest, p.area_m2);
        null.mean_n_patches = compute_mean(counts);
        null.p05_n_patches = compute_percentile(counts, 5.0);
        null.p95_n_patches = compute_percentile(counts, 95.0);
        null.mean_area_m2 = compute_mean(totals);
        null.p05_area_m2 = compute_percentile(totals, 5.0);
        null.p95_area_m2 = compute_percentile(totals, 95.0);
        const double observed = out.loss_area_m2 + out.gain_area_m2;
        if (observed > 0.0) null.area_false_discovery_ratio = std::min(1.0, null.mean_area_m2 / observed);

        // Quarter-decade bins of area, starting at one cell.
        const int n_bins = static_cast<int>(std::ceil(4.0 * std::log10(largest / cell_area) + 1e-9)) + 1;
        for (int b = 0; b < n_bins; ++b)
        {
            dem_change_null_bin bin;
            bin.area_lower_m2 = cell_area * std::pow(10.0, b / 4.0);
            bin.area_upper_m2 = cell_area * std::pow(10.0, (b + 1) / 4.0);
            for (const auto& p : out.patches)
                if (p.area_m2 >= bin.area_lower_m2 && p.area_m2 < bin.area_upper_m2) ++bin.n_observed;
            std::int64_t n_simulated = 0;
            for (const auto& areas : sim_areas)
                for (const double a : areas)
                    if (a >= bin.area_lower_m2 && a < bin.area_upper_m2) ++n_simulated;
            bin.mean_n_simulated = static_cast<double>(n_simulated) / n_sim;
            if (bin.n_observed > 0)
                bin.false_discovery_ratio = std::min(1.0, bin.mean_n_simulated / static_cast<double>(bin.n_observed));
            null.by_area.push_back(bin);
        }
    }
    return out;
}

std::vector<std::uint8_t> dem_change_channel_mask(const lsdtt_raster_cpu& routing_dem,
                                                  double specific_catchment_threshold_m, int buffer_pixels,
                                                  const std::string& boundary_conditions, int n_threads)
{
    if (!(specific_catchment_threshold_m > 0.0))
        throw std::invalid_argument("dem_change_channel_mask: the specific catchment threshold must be > 0");
    if (buffer_pixels < 0) throw std::invalid_argument("dem_change_channel_mask: the buffer must be >= 0");
    const lsdtt_raster_meta& meta = routing_dem.get_meta();
    meta.check_projected_metres();
    lsdtt_flowgraph_options flow_options;
    flow_options.method = lsdtt_flowgraph_method::d_infinity;
    flow_options.boundary_conditions = boundary_conditions;
    flow_options.n_threads = n_threads;
    const lsdtt_flowgraph_cpu graph(routing_dem.get_data_array(), routing_dem.get_mask_array(), meta, flow_options);
    const std::vector<double> area = graph.accumulate();
    const std::vector<std::uint8_t>& valid = routing_dem.get_mask();
    const int nr = meta.n_rows, nc = meta.n_cols;
    const std::size_t n = static_cast<std::size_t>(nr) * nc;
    const double area_threshold = specific_catchment_threshold_m * meta.cell_size;
    std::vector<std::uint8_t> channel(n, 0u);
    for (std::size_t i = 0; i < n; ++i)
        if (valid[i] && area[i] >= area_threshold) channel[i] = 1u;
    if (buffer_pixels == 0) return channel;
    std::vector<std::pair<int, int>> disc;
    for (int dr = -buffer_pixels; dr <= buffer_pixels; ++dr)
        for (int dc = -buffer_pixels; dc <= buffer_pixels; ++dc)
            if ((dr != 0 || dc != 0) && dr * dr + dc * dc <= buffer_pixels * buffer_pixels) disc.emplace_back(dr, dc);
    std::vector<std::uint8_t> grown = channel;
    for (int r = 0; r < nr; ++r)
        for (int c = 0; c < nc; ++c)
        {
            if (!channel[static_cast<std::size_t>(r) * nc + c]) continue;
            for (const auto& d : disc)
            {
                const int r2 = r + d.first, c2 = c + d.second;
                if (r2 >= 0 && r2 < nr && c2 >= 0 && c2 < nc) grown[static_cast<std::size_t>(r2) * nc + c2] = 1u;
            }
        }
    return grown;
}

dem_change_scaling fit_dem_change_scaling(const std::vector<dem_change_patch>& patches, int sign, int n_bootstrap,
                                          std::uint64_t seed, std::size_t min_tail)
{
    if (sign != -1 && sign != 1) throw std::invalid_argument("fit_dem_change_scaling: sign must be -1 or +1");
    dem_change_scaling out;
    out.sign = sign;
    std::vector<double> area, volume, depth;
    for (const dem_change_patch& p : patches)
    {
        if (p.sign != sign) continue;
        area.push_back(p.area_m2);
        volume.push_back(std::abs(p.volume_m3));
        depth.push_back(std::abs(p.mean_change_m));
    }
    out.n_patches = area.size();
    out.volume_area = fit_power_law_scaling(area, volume, n_bootstrap, seed);
    out.depth_area = fit_power_law_scaling(area, depth, n_bootstrap, seed + 1);
    out.volume_tail = fit_power_law_tail(volume, false, std::max<std::size_t>(2, min_tail));
    return out;
}

} // namespace lsdtt
