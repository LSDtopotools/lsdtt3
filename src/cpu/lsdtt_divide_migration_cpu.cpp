// lsdtt_divide_migration_cpu.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html
//
// -----------------------------------------------------------------------------
// CPU compute routines for divide-migration analysis.
// -----------------------------------------------------------------------------

#include "lsdtt_divide_migration_cpu.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <map>
#include <random>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

#include "lsdtt_chi_extraction_cpu.hpp"

namespace lsdtt
{

namespace
{

constexpr double k_nan = std::numeric_limits<double>::quiet_NaN();

// Distance below which a head is treated as sitting on the divide, so that the
// gradient is reported as NaN instead of exploding.
constexpr double k_min_divide_distance_m = 1e-6;

// -----------------------------------------------------------------------------
// One straight piece of divide geometry, tagged with the segment record it came
// from so that pair_id can be reported for the winning piece.
// -----------------------------------------------------------------------------
struct divide_piece
{
    lsdtt_xy a{};
    lsdtt_xy b{};
    const lsdtt_divide_segment_record* owner = nullptr;
};

struct nearest_piece_result
{
    lsdtt_xy point{};
    double   distance_sq = std::numeric_limits<double>::infinity();
    const lsdtt_divide_segment_record* owner = nullptr;
};

// Closest point on the line segment a->b to p, plus the squared distance.
nearest_piece_result closest_point_on_piece(
    const lsdtt_xy& p,
    const divide_piece& piece)
{
    const double vx = piece.b.x - piece.a.x;
    const double vy = piece.b.y - piece.a.y;
    const double len_sq = vx * vx + vy * vy;

    double t = 0.0;
    if (len_sq > 0.0)
    {
        t = ((p.x - piece.a.x) * vx + (p.y - piece.a.y) * vy) / len_sq;
        t = std::min(1.0, std::max(0.0, t));
    }

    nearest_piece_result out;
    out.point.x = piece.a.x + t * vx;
    out.point.y = piece.a.y + t * vy;
    const double dx = p.x - out.point.x;
    const double dy = p.y - out.point.y;
    out.distance_sq = dx * dx + dy * dy;
    out.owner = piece.owner;
    return out;
}

// -----------------------------------------------------------------------------
// Uniform-grid spatial index over divide pieces.
//
// Each piece is stored in every bucket overlapped by its bounding box, so a
// bucket-by-bucket ring search around a query point never misses a candidate.
// -----------------------------------------------------------------------------
class divide_piece_grid
{
public:
    void build(std::vector<divide_piece> pieces)
    {
        pieces_ = std::move(pieces);
        buckets_.clear();
        if (pieces_.empty()) return;

        double x_max = -std::numeric_limits<double>::infinity();
        double y_max = -std::numeric_limits<double>::infinity();
        x_min_ = std::numeric_limits<double>::infinity();
        y_min_ = std::numeric_limits<double>::infinity();
        for (const divide_piece& piece : pieces_)
        {
            x_min_ = std::min({x_min_, piece.a.x, piece.b.x});
            y_min_ = std::min({y_min_, piece.a.y, piece.b.y});
            x_max  = std::max({x_max,  piece.a.x, piece.b.x});
            y_max  = std::max({y_max,  piece.a.y, piece.b.y});
        }

        // Aim for a handful of pieces per bucket. Degenerate (zero-area)
        // extents fall back to a single bucket.
        const double width  = std::max(0.0, x_max - x_min_);
        const double height = std::max(0.0, y_max - y_min_);
        const double area   = width * height;
        const double target_buckets =
            std::max(1.0, static_cast<double>(pieces_.size()) / 4.0);
        bucket_size_ = (area > 0.0) ? std::sqrt(area / target_buckets)
                                    : std::max(width, height);
        if (!(bucket_size_ > 0.0)) bucket_size_ = 1.0;

        n_x_ = std::max(1, static_cast<int>(std::floor(width  / bucket_size_)) + 1);
        n_y_ = std::max(1, static_cast<int>(std::floor(height / bucket_size_)) + 1);
        buckets_.resize(static_cast<std::size_t>(n_x_) * static_cast<std::size_t>(n_y_));

        for (std::size_t i = 0; i < pieces_.size(); ++i)
        {
            const divide_piece& piece = pieces_[i];
            const int i0 = bucket_x(std::min(piece.a.x, piece.b.x));
            const int i1 = bucket_x(std::max(piece.a.x, piece.b.x));
            const int j0 = bucket_y(std::min(piece.a.y, piece.b.y));
            const int j1 = bucket_y(std::max(piece.a.y, piece.b.y));
            for (int j = j0; j <= j1; ++j)
            {
                for (int i_b = i0; i_b <= i1; ++i_b)
                {
                    buckets_[bucket_index(i_b, j)].push_back(i);
                }
            }
        }
    }

    [[nodiscard]] bool empty() const { return pieces_.empty(); }

    // Nearest point on any indexed piece. Returns a result with a null owner if
    // the index holds no geometry.
    [[nodiscard]] nearest_piece_result query_nearest(const lsdtt_xy& p) const
    {
        nearest_piece_result best;
        if (pieces_.empty()) return best;

        const int cx = bucket_x(p.x);
        const int cy = bucket_y(p.y);
        const int r_max = std::max({cx, n_x_ - 1 - cx, cy, n_y_ - 1 - cy});

        for (int r = 0; r <= r_max; ++r)
        {
            // Everything within Chebyshev ring r-1 has been scanned already, so
            // anything still unseen is at least (r-1) buckets away.
            if (best.owner != nullptr && r > 0)
            {
                const double scanned = static_cast<double>(r - 1) * bucket_size_;
                if (scanned > 0.0 && best.distance_sq <= scanned * scanned) break;
            }

            for (int j = cy - r; j <= cy + r; ++j)
            {
                if (j < 0 || j >= n_y_) continue;
                const bool j_on_edge = (j == cy - r || j == cy + r);
                for (int i_b = cx - r; i_b <= cx + r; ++i_b)
                {
                    if (i_b < 0 || i_b >= n_x_) continue;
                    if (!j_on_edge && i_b != cx - r && i_b != cx + r) continue;

                    for (const std::size_t piece_idx : buckets_[bucket_index(i_b, j)])
                    {
                        const nearest_piece_result cand =
                            closest_point_on_piece(p, pieces_[piece_idx]);
                        if (cand.distance_sq < best.distance_sq)
                        {
                            best = cand;
                        }
                    }
                }
            }
        }
        return best;
    }

private:
    [[nodiscard]] int bucket_x(double x) const
    {
        const int i = static_cast<int>(std::floor((x - x_min_) / bucket_size_));
        return std::min(n_x_ - 1, std::max(0, i));
    }

    [[nodiscard]] int bucket_y(double y) const
    {
        const int j = static_cast<int>(std::floor((y - y_min_) / bucket_size_));
        return std::min(n_y_ - 1, std::max(0, j));
    }

    [[nodiscard]] std::size_t bucket_index(int i, int j) const
    {
        return static_cast<std::size_t>(j) * static_cast<std::size_t>(n_x_) +
               static_cast<std::size_t>(i);
    }

    std::vector<divide_piece> pieces_;
    std::vector<std::vector<std::size_t>> buckets_;
    double x_min_ = 0.0;
    double y_min_ = 0.0;
    double bucket_size_ = 1.0;
    int n_x_ = 1;
    int n_y_ = 1;
};

// Elevation at the cell containing (x, y), or NaN if the point falls outside
// the raster or on a nodata cell.
double sample_elevation_at_point(
    const lsdtt_raster_cpu& elevation,
    double x,
    double y)
{
    const lsdtt_raster_meta& meta = elevation.get_meta();
    if (!(meta.cell_size > 0.0)) return k_nan;

    const int col = static_cast<int>(std::floor((x - meta.x_min) / meta.cell_size));
    const int row = static_cast<int>(std::floor((meta.y_max() - y) / meta.cell_size));
    if (col < 0 || col >= meta.n_cols || row < 0 || row >= meta.n_rows) return k_nan;

    const std::size_t idx =
        static_cast<std::size_t>(row) * static_cast<std::size_t>(meta.n_cols) +
        static_cast<std::size_t>(col);

    const std::vector<float>& data = elevation.get_data();
    const std::vector<std::uint8_t>& mask = elevation.get_mask();
    if (idx >= data.size() || idx >= mask.size() || mask[idx] == 0) return k_nan;

    const double z = static_cast<double>(data[idx]);
    if (!std::isfinite(z)) return k_nan;
    return z;
}

// Elevation at a flat index, or NaN for nodata / out-of-range indices.
double sample_elevation_at_flat(
    const lsdtt_raster_cpu& elevation,
    std::int64_t flat)
{
    if (flat < 0) return k_nan;
    const std::size_t idx = static_cast<std::size_t>(flat);
    const std::vector<float>& data = elevation.get_data();
    const std::vector<std::uint8_t>& mask = elevation.get_mask();
    if (idx >= data.size() || idx >= mask.size() || mask[idx] == 0) return k_nan;

    const double z = static_cast<double>(data[idx]);
    if (!std::isfinite(z)) return k_nan;
    return z;
}

// Squared distance between two map points.
double distance_sq_between(const lsdtt_xy& a, const lsdtt_xy& b)
{
    const double dx = a.x - b.x;
    const double dy = a.y - b.y;
    return dx * dx + dy * dy;
}

// Summed euclidean length of a polyline.
double polyline_length(const std::vector<lsdtt_xy>& points)
{
    double total = 0.0;
    for (std::size_t i = 1; i < points.size(); ++i)
    {
        total += std::sqrt(distance_sq_between(points[i - 1], points[i]));
    }
    return total;
}

// -----------------------------------------------------------------------------
// Pair-summary statistics helpers
// -----------------------------------------------------------------------------

// Mean of a sample, NaN when the sample is empty. Callers only ever pass
// samples from which the NaNs have already been dropped.
double sample_mean(const std::vector<double>& v)
{
    if (v.empty()) return k_nan;
    double sum = 0.0;
    for (const double x : v) sum += x;
    return sum / static_cast<double>(v.size());
}

// Sample standard deviation with the n-1 denominator, matching MATLAB std() as
// used by DivideTools. NaN for fewer than two values.
double sample_std_dev(const std::vector<double>& v)
{
    const std::size_t n = v.size();
    if (n < 2) return k_nan;
    const double m = sample_mean(v);
    double ss = 0.0;
    for (const double x : v)
    {
        const double d = x - m;
        ss += d * d;
    }
    return std::sqrt(ss / static_cast<double>(n - 1));
}

// Continued-fraction evaluation of the incomplete beta function (Lentz's
// method), the standard Numerical Recipes betacf. Used only for the Welch
// t-test p-value.
double beta_continued_fraction(double a, double b, double x)
{
    constexpr int    k_max_iter = 300;
    constexpr double k_eps      = 3.0e-16;
    constexpr double k_tiny     = 1.0e-300;

    const double qab = a + b;
    const double qap = a + 1.0;
    const double qam = a - 1.0;

    double c = 1.0;
    double d = 1.0 - qab * x / qap;
    if (std::fabs(d) < k_tiny) d = k_tiny;
    d = 1.0 / d;
    double h = d;

    for (int m = 1; m <= k_max_iter; ++m)
    {
        const double dm = static_cast<double>(m);
        const double m2 = 2.0 * dm;

        double aa = dm * (b - dm) * x / ((qam + m2) * (a + m2));
        d = 1.0 + aa * d;
        if (std::fabs(d) < k_tiny) d = k_tiny;
        c = 1.0 + aa / c;
        if (std::fabs(c) < k_tiny) c = k_tiny;
        d = 1.0 / d;
        h *= d * c;

        aa = -(a + dm) * (qab + dm) * x / ((a + m2) * (qap + m2));
        d = 1.0 + aa * d;
        if (std::fabs(d) < k_tiny) d = k_tiny;
        c = 1.0 + aa / c;
        if (std::fabs(c) < k_tiny) c = k_tiny;
        d = 1.0 / d;
        const double del = d * c;
        h *= del;

        if (std::fabs(del - 1.0) < k_eps) break;
    }
    return h;
}

// Regularised incomplete beta function I_x(a, b).
double regularised_incomplete_beta(double a, double b, double x)
{
    if (!(x >= 0.0) || !(x <= 1.0)) return k_nan;
    if (x == 0.0) return 0.0;
    if (x == 1.0) return 1.0;

    const double log_front = std::lgamma(a + b) - std::lgamma(a) - std::lgamma(b) +
                             a * std::log(x) + b * std::log1p(-x);
    const double front = std::exp(log_front);

    if (x < (a + 1.0) / (a + b + 2.0))
    {
        return front * beta_continued_fraction(a, b, x) / a;
    }
    return 1.0 - front * beta_continued_fraction(b, a, 1.0 - x) / b;
}

// Two-sided p-value of Welch's unequal-variance t-test on the difference of
// means, as DivideTools calls ttest2(...,'Vartype','unequal'). NaN when either
// side has fewer than two values or when both sides are constant, in which case
// there is no within-side spread to test against.
double welch_t_test_p_value(
    const std::vector<double>& a,
    const std::vector<double>& b)
{
    const std::size_t na = a.size();
    const std::size_t nb = b.size();
    if (na < 2 || nb < 2) return k_nan;

    const double sa = sample_std_dev(a);
    const double sb = sample_std_dev(b);
    if (!std::isfinite(sa) || !std::isfinite(sb)) return k_nan;

    const double va = sa * sa / static_cast<double>(na);
    const double vb = sb * sb / static_cast<double>(nb);
    const double se = std::sqrt(va + vb);
    if (!(se > 0.0)) return k_nan;

    const double t = (sample_mean(a) - sample_mean(b)) / se;

    // Welch-Satterthwaite degrees of freedom.
    const double df_num = (va + vb) * (va + vb);
    const double df_den = va * va / static_cast<double>(na - 1) +
                          vb * vb / static_cast<double>(nb - 1);
    if (!(df_den > 0.0)) return k_nan;
    const double df = df_num / df_den;

    const double p = regularised_incomplete_beta(
        0.5 * df, 0.5, df / (df + t * t));
    if (!std::isfinite(p)) return k_nan;
    return std::min(1.0, std::max(0.0, p));
}

// 32-bit FNV-1a. A fixed, self-contained hash so that the bootstrap seeding
// does not depend on std::hash, which is implementation defined.
std::uint32_t fnv1a_32(const std::string& text)
{
    std::uint32_t h = 2166136261u;
    for (const char c : text)
    {
        h ^= static_cast<std::uint32_t>(static_cast<unsigned char>(c));
        h *= 16777619u;
    }
    return h;
}

// Bootstrap 95% confidence-interval half-width of the DELTA of side means.
//
// Both sides are resampled with replacement n_replicates times and the delta of
// the resampled means recorded; the reported value is the mean of the distances
// from the bootstrap mean delta to the 2.5th and 97.5th percentiles, which is
// the statistic DivideTools' bootCI helper returns for a single side.
//
// Reproducibility: the generator is a std::mt19937 seeded through a
// std::seed_seq from the caller's seed, a fixed hash of the pair_id and the
// metric index, and indices are drawn with a plain modulo of the generator
// output rather than through std::uniform_int_distribution, whose mapping is
// implementation defined. The result therefore depends on nothing but the
// inputs and the seed.
double bootstrap_delta_ci_half_width(
    const std::vector<double>& a,
    const std::vector<double>& b,
    int                        n_replicates,
    std::uint32_t              seed,
    const std::string&         pair_id,
    int                        metric_index)
{
    const std::size_t na = a.size();
    const std::size_t nb = b.size();
    if (na < 2 || nb < 2) return k_nan;
    if (n_replicates < 2) return k_nan;

    std::seed_seq seq{seed, fnv1a_32(pair_id),
                      static_cast<std::uint32_t>(metric_index)};
    std::mt19937 gen(seq);

    std::vector<double> deltas;
    deltas.reserve(static_cast<std::size_t>(n_replicates));
    for (int r = 0; r < n_replicates; ++r)
    {
        double sum_a = 0.0;
        for (std::size_t i = 0; i < na; ++i)
        {
            sum_a += a[static_cast<std::size_t>(gen()) % na];
        }
        double sum_b = 0.0;
        for (std::size_t i = 0; i < nb; ++i)
        {
            sum_b += b[static_cast<std::size_t>(gen()) % nb];
        }
        deltas.push_back(sum_a / static_cast<double>(na) -
                         sum_b / static_cast<double>(nb));
    }

    const double mean_delta = sample_mean(deltas);
    std::sort(deltas.begin(), deltas.end());

    // Nearest-rank percentiles: q(p) = sorted[ceil(p * B) - 1].
    const double b_count = static_cast<double>(deltas.size());
    const auto rank = [&](double p) {
        const long long idx =
            static_cast<long long>(std::ceil(p * b_count)) - 1;
        const long long clamped = std::min<long long>(
            static_cast<long long>(deltas.size()) - 1,
            std::max<long long>(0, idx));
        return deltas[static_cast<std::size_t>(clamped)];
    };

    const double lo = rank(0.025);
    const double hi = rank(0.975);
    return 0.5 * (std::fabs(lo - mean_delta) + std::fabs(hi - mean_delta));
}

// The uncertainty reported for one metric, given the two sides' valid values.
// Every method needs a within-side spread, so all of them require at least two
// valid values on both sides and return NaN otherwise.
double pair_uncertainty(
    const std::vector<double>&          a,
    const std::vector<double>&          b,
    const lsdtt_pair_summary_options&   options,
    const std::string&                  pair_id,
    int                                 metric_index)
{
    const std::size_t na = a.size();
    const std::size_t nb = b.size();
    if (na < 2 || nb < 2) return k_nan;

    switch (options.uncertainty_method)
    {
    case lsdtt_pair_uncertainty_method::std_dev:
    {
        // Pooled standard deviation of the two sides.
        const double sa = sample_std_dev(a);
        const double sb = sample_std_dev(b);
        if (!std::isfinite(sa) || !std::isfinite(sb)) return k_nan;
        const double num = static_cast<double>(na - 1) * sa * sa +
                           static_cast<double>(nb - 1) * sb * sb;
        const double den = static_cast<double>(na + nb - 2);
        if (!(den > 0.0)) return k_nan;
        return std::sqrt(num / den);
    }
    case lsdtt_pair_uncertainty_method::std_err:
    {
        // Standard error of the difference of means.
        const double sa = sample_std_dev(a);
        const double sb = sample_std_dev(b);
        if (!std::isfinite(sa) || !std::isfinite(sb)) return k_nan;
        return std::sqrt(sa * sa / static_cast<double>(na) +
                         sb * sb / static_cast<double>(nb));
    }
    case lsdtt_pair_uncertainty_method::bootstrap:
        return bootstrap_delta_ci_half_width(
            a, b, options.bootstrap_replicates, options.bootstrap_seed,
            pair_id, metric_index);
    case lsdtt_pair_uncertainty_method::ttest:
        return welch_t_test_p_value(a, b);
    }
    return k_nan;
}

// Side A / side B valid-value samples for one metric of one pair.
struct pair_metric_sample
{
    std::vector<double> a;
    std::vector<double> b;
};

// Everything gathered for one pair before the summary row is built.
struct pair_accumulator
{
    int unit_a = -1;
    int unit_b = -1;
    std::int64_t n_heads_a = 0;
    std::int64_t n_heads_b = 0;
    // Indexed by metric: 0 head elevation, 1 gradient, 2 relief, 3 chi. The
    // order is also the bootstrap metric index, so it must not be reshuffled.
    pair_metric_sample metrics[4];
};

// Splits "<unit_a>|<unit_b>" into its two unit ids. Returns false for anything
// that is not that shape, so a malformed pair id is skipped rather than
// silently summarised against the wrong units.
bool parse_pair_id(const std::string& pair_id, int& unit_a, int& unit_b)
{
    const std::size_t bar = pair_id.find('|');
    if (bar == std::string::npos || bar == 0 || bar + 1 >= pair_id.size())
    {
        return false;
    }
    try
    {
        std::size_t used_a = 0;
        std::size_t used_b = 0;
        const std::string text_a = pair_id.substr(0, bar);
        const std::string text_b = pair_id.substr(bar + 1);
        const int a = std::stoi(text_a, &used_a);
        const int b = std::stoi(text_b, &used_b);
        if (used_a != text_a.size() || used_b != text_b.size()) return false;
        unit_a = a;
        unit_b = b;
        return true;
    }
    catch (const std::exception&)
    {
        return false;
    }
}

} // namespace

// -----------------------------------------------------------------------------
// merge_divide_fragments
// -----------------------------------------------------------------------------
std::vector<lsdtt_divide_segment_record> merge_divide_fragments(
    const std::vector<lsdtt_divide_segment_record>& fragments,
    double                                          cell_size_m)
{
    std::vector<lsdtt_divide_segment_record> merged;
    if (fragments.empty()) return merged;

    // One cell diagonal: see the criterion note in the header.
    const double tolerance =
        (cell_size_m > 0.0) ? cell_size_m * std::sqrt(2.0) : 0.0;
    const double tolerance_sq = tolerance * tolerance;

    const auto joins = [&](const lsdtt_xy& a, const lsdtt_xy& b) {
        return distance_sq_between(a, b) <= tolerance_sq;
    };

    std::vector<char> used(fragments.size(), 0);

    for (std::size_t seed = 0; seed < fragments.size(); ++seed)
    {
        if (used[seed] != 0) continue;
        used[seed] = 1;

        const lsdtt_divide_segment_record& first = fragments[seed];

        // An empty fragment carries no geometry to chain with, so it is passed
        // straight through rather than silently dropped.
        std::vector<lsdtt_xy> chain = first.points;
        int n_used = 1;

        // A chain that has been built from two or more fragments and whose two
        // ends now meet is a closed ring: stop, so nothing is attached across
        // the closure point and the loop cannot run round again. The test is
        // deliberately not applied to a lone fragment, whose two ends can sit
        // within a cell diagonal of each other without the fragment being a
        // ring at all.
        const auto chain_is_closed = [&](int n_in_chain) {
            return n_in_chain > 1 && chain.size() > 2 &&
                   joins(chain.front(), chain.back());
        };

        // Each pass consumes at most one fragment and every fragment is marked
        // used before it is appended, so this terminates for any input.
        bool grew = !chain.empty();
        while (grew && !chain_is_closed(n_used))
        {
            grew = false;
            for (std::size_t k = 0; k < fragments.size(); ++k)
            {
                if (used[k] != 0) continue;
                const lsdtt_divide_segment_record& cand = fragments[k];
                if (cand.pair_id != first.pair_id) continue;
                if (cand.points.empty()) continue;

                const lsdtt_xy& c_front = cand.points.front();
                const lsdtt_xy& c_back  = cand.points.back();

                if (joins(chain.back(), c_front) || joins(chain.back(), c_back))
                {
                    // Append at the tail, reversing the candidate when it is
                    // the candidate's END that meets the chain. The shared
                    // vertex is dropped so it is not written twice.
                    const bool reverse = !joins(chain.back(), c_front);
                    if (reverse)
                    {
                        for (std::size_t i = cand.points.size() - 1; i-- > 0;)
                        {
                            chain.push_back(cand.points[i]);
                        }
                    }
                    else
                    {
                        for (std::size_t i = 1; i < cand.points.size(); ++i)
                        {
                            chain.push_back(cand.points[i]);
                        }
                    }
                }
                else if (joins(chain.front(), c_back) ||
                         joins(chain.front(), c_front))
                {
                    // Prepend at the head, again dropping the shared vertex.
                    const bool reverse = !joins(chain.front(), c_back);
                    std::vector<lsdtt_xy> prefix;
                    prefix.reserve(cand.points.size() - 1 + chain.size());
                    if (reverse)
                    {
                        for (std::size_t i = cand.points.size(); i-- > 1;)
                        {
                            prefix.push_back(cand.points[i]);
                        }
                    }
                    else
                    {
                        for (std::size_t i = 0; i + 1 < cand.points.size(); ++i)
                        {
                            prefix.push_back(cand.points[i]);
                        }
                    }
                    prefix.insert(prefix.end(), chain.begin(), chain.end());
                    chain = std::move(prefix);
                }
                else
                {
                    continue;
                }

                used[k] = 1;
                ++n_used;
                grew = true;
                break;
            }
        }

        lsdtt_divide_segment_record out;
        out.pair_id = first.pair_id;
        out.unit_a = first.unit_a;
        out.unit_b = first.unit_b;
        out.n_fragments = n_used;
        // Junction numbers describe one basin pair, so they only survive on a
        // record that still holds exactly one fragment.
        if (n_used == 1)
        {
            out.basin_junction_a = first.basin_junction_a;
            out.basin_junction_b = first.basin_junction_b;
        }
        out.segment_id = static_cast<int>(merged.size()) + 1;
        out.points = std::move(chain);
        out.length_m = polyline_length(out.points);
        merged.push_back(std::move(out));
    }

    return merged;
}

// -----------------------------------------------------------------------------
// find_adjacent_basin_pairs
// -----------------------------------------------------------------------------
std::vector<std::pair<int, int>> find_adjacent_basin_pairs(
    const std::vector<lsdtt_flat_index_slice>& basin_cells,
    int                                        n_rows,
    int                                        n_cols,
    bool                                       periodic_ns,
    bool                                       periodic_ew)
{
    std::vector<std::pair<int, int>> pairs;
    const std::size_t n_basins = basin_cells.size();
    if (n_basins < 2 || n_rows <= 0 || n_cols <= 0) return pairs;

    const std::size_t n_cells =
        static_cast<std::size_t>(n_rows) * static_cast<std::size_t>(n_cols);

    // Cell -> basin labels, as a CSR layout. Basins may overlap (a nested basin
    // shares every one of its cells with its parent), so a cell can carry more
    // than one label and a single label grid would hide contacts by overwriting.
    // The offsets are 64-bit because they index the concatenated basin cells,
    // not the grid.
    std::vector<std::int64_t> offsets(n_cells + 1, 0);
    for (const lsdtt_flat_index_slice& slice : basin_cells)
    {
        for (const std::int64_t flat : slice)
        {
            if (flat < 0 || static_cast<std::size_t>(flat) >= n_cells) continue;
            ++offsets[static_cast<std::size_t>(flat) + 1];
        }
    }
    for (std::size_t c = 0; c < n_cells; ++c)
    {
        offsets[c + 1] += offsets[c];
    }

    std::vector<std::int32_t> labels(static_cast<std::size_t>(offsets[n_cells]));
    std::vector<std::int64_t> cursor(offsets.begin(), offsets.end() - 1);
    for (std::size_t b = 0; b < n_basins; ++b)
    {
        for (const std::int64_t flat : basin_cells[b])
        {
            if (flat < 0 || static_cast<std::size_t>(flat) >= n_cells) continue;
            const std::size_t c = static_cast<std::size_t>(flat);
            labels[static_cast<std::size_t>(cursor[c]++)] =
                static_cast<std::int32_t>(b);
        }
    }

    // Pairs are keyed as i * n_basins + j so that a plain sort/unique
    // deduplicates them and leaves them in (i, j) ascending order -- the order
    // the old nested loop visited them in.
    std::vector<std::uint64_t> keys;
    const std::uint64_t stride = static_cast<std::uint64_t>(n_basins);

    auto record = [&](std::int32_t p, std::int32_t q) {
        if (p == q) return;
        const std::int32_t lo = (p < q) ? p : q;
        const std::int32_t hi = (p < q) ? q : p;
        keys.push_back(static_cast<std::uint64_t>(lo) * stride +
                       static_cast<std::uint64_t>(hi));
    };

    // Only the four "forward" neighbours are visited, which covers every
    // eight-connected cell pair exactly once.
    constexpr int k_d_row[4] = {0, 1, 1, 1};
    constexpr int k_d_col[4] = {1, -1, 0, 1};

    for (int row = 0; row < n_rows; ++row)
    {
        for (int col = 0; col < n_cols; ++col)
        {
            const std::size_t cell =
                static_cast<std::size_t>(row) * static_cast<std::size_t>(n_cols) +
                static_cast<std::size_t>(col);
            const std::int64_t begin = offsets[cell];
            const std::int64_t end   = offsets[cell + 1];
            if (begin == end) continue;

            // Two basins covering the SAME cell overlap, which is contact too.
            for (std::int64_t u = begin; u < end; ++u)
            {
                for (std::int64_t v = u + 1; v < end; ++v)
                {
                    record(labels[static_cast<std::size_t>(u)],
                           labels[static_cast<std::size_t>(v)]);
                }
            }

            for (int k = 0; k < 4; ++k)
            {
                int n_row = row + k_d_row[k];
                int n_col = col + k_d_col[k];

                if (periodic_ns && n_rows > 1)
                {
                    n_row = (n_row + n_rows) % n_rows;
                }
                if (periodic_ew && n_cols > 1)
                {
                    n_col = (n_col + n_cols) % n_cols;
                }
                if (n_row < 0 || n_row >= n_rows) continue;
                if (n_col < 0 || n_col >= n_cols) continue;

                const std::size_t neighbour =
                    static_cast<std::size_t>(n_row) *
                        static_cast<std::size_t>(n_cols) +
                    static_cast<std::size_t>(n_col);
                if (neighbour == cell) continue;

                for (std::int64_t u = begin; u < end; ++u)
                {
                    for (std::int64_t v = offsets[neighbour];
                         v < offsets[neighbour + 1]; ++v)
                    {
                        record(labels[static_cast<std::size_t>(u)],
                               labels[static_cast<std::size_t>(v)]);
                    }
                }
            }
        }
    }

    std::sort(keys.begin(), keys.end());
    keys.erase(std::unique(keys.begin(), keys.end()), keys.end());

    pairs.reserve(keys.size());
    for (const std::uint64_t key : keys)
    {
        pairs.emplace_back(static_cast<int>(key / stride),
                           static_cast<int>(key % stride));
    }
    return pairs;
}

// -----------------------------------------------------------------------------
// compute_nearest_divide_metrics
// -----------------------------------------------------------------------------
std::size_t compute_nearest_divide_metrics(
    const std::vector<lsdtt_divide_segment_record>& divide_segments,
    const lsdtt_raster_cpu&                         elevation,
    std::vector<lsdtt_channel_head_metric_record>&  head_metrics)
{
    // One index per unit: a divide belongs to both of the units it separates,
    // so a head only ever searches divides that bound its own unit.
    std::map<int, std::vector<divide_piece>> pieces_by_unit;
    for (const lsdtt_divide_segment_record& seg : divide_segments)
    {
        if (seg.points.empty()) continue;

        std::vector<divide_piece> pieces;
        if (seg.points.size() == 1)
        {
            // A single-vertex divide is a degenerate zero-length piece.
            pieces.push_back(divide_piece{seg.points[0], seg.points[0], &seg});
        }
        else
        {
            pieces.reserve(seg.points.size() - 1);
            for (std::size_t k = 1; k < seg.points.size(); ++k)
            {
                pieces.push_back(divide_piece{seg.points[k - 1], seg.points[k], &seg});
            }
        }

        std::vector<divide_piece>& a_list = pieces_by_unit[seg.unit_a];
        a_list.insert(a_list.end(), pieces.begin(), pieces.end());
        if (seg.unit_b != seg.unit_a)
        {
            std::vector<divide_piece>& b_list = pieces_by_unit[seg.unit_b];
            b_list.insert(b_list.end(), pieces.begin(), pieces.end());
        }
    }

    std::map<int, divide_piece_grid> grid_by_unit;
    for (auto& entry : pieces_by_unit)
    {
        grid_by_unit[entry.first].build(std::move(entry.second));
    }

    std::size_t n_matched = 0;
    for (lsdtt_channel_head_metric_record& rec : head_metrics)
    {
        rec.pair_id_nearest = "";
        rec.nearest_divide_x = k_nan;
        rec.nearest_divide_y = k_nan;
        rec.divide_elevation_m = k_nan;
        rec.distance_to_divide_m = k_nan;
        rec.relief_m = k_nan;
        rec.gradient_to_divide = k_nan;

        auto git = grid_by_unit.find(rec.unit_id);
        if (git == grid_by_unit.end() || git->second.empty()) continue;

        const lsdtt_xy head_xy{rec.head_x, rec.head_y};
        const nearest_piece_result hit = git->second.query_nearest(head_xy);
        if (hit.owner == nullptr) continue;

        rec.pair_id_nearest = hit.owner->pair_id;
        rec.nearest_divide_x = hit.point.x;
        rec.nearest_divide_y = hit.point.y;
        rec.distance_to_divide_m = std::sqrt(hit.distance_sq);
        ++n_matched;

        const double divide_z =
            sample_elevation_at_point(elevation, hit.point.x, hit.point.y);
        if (!std::isfinite(divide_z)) continue;

        rec.divide_elevation_m = divide_z;
        rec.relief_m = divide_z - rec.head_elevation_m;

        // A head sitting on its own divide has no meaningful gradient.
        if (rec.distance_to_divide_m > k_min_divide_distance_m)
        {
            rec.gradient_to_divide = rec.relief_m / rec.distance_to_divide_m;
        }
    }

    return n_matched;
}

// -----------------------------------------------------------------------------
// apply_common_base_level_to_chi
// -----------------------------------------------------------------------------
std::size_t apply_common_base_level_to_chi(
    lsdtt_chi_collection_result& chi,
    const lsdtt_raster_cpu&      elevation,
    double                       base_elevation_m)
{
    std::size_t n_removed = 0;

    for (lsdtt_chi_basin_result& basin : chi.basins)
    {
        const std::size_t n_nodes = basin.flats.size();
        if (n_nodes == 0) continue;
        if (basin.chi.size() != n_nodes || basin.receiver_flats.size() != n_nodes)
        {
            continue;
        }

        std::unordered_map<std::int64_t, std::size_t> flat_to_index;
        flat_to_index.reserve(n_nodes * 2);
        for (std::size_t i = 0; i < n_nodes; ++i)
        {
            flat_to_index[basin.flats[i]] = i;
        }

        // Retained nodes sit at or above the common base level. In a filled DEM
        // this set is upslope closed, so each retained fragment is a subtree.
        std::vector<char> keep(n_nodes, 0);
        for (std::size_t i = 0; i < n_nodes; ++i)
        {
            const double z = sample_elevation_at_flat(elevation, basin.flats[i]);
            keep[i] = (std::isfinite(z) && z >= base_elevation_m) ? 1 : 0;
        }

        // chi offset of each retained node is the chi of the root of its
        // fragment, resolved by walking downstream while the receiver is kept.
        constexpr std::size_t k_unset = static_cast<std::size_t>(-1);
        std::vector<std::size_t> root_of(n_nodes, k_unset);
        std::vector<std::size_t> walk;
        for (std::size_t i = 0; i < n_nodes; ++i)
        {
            if (keep[i] == 0 || root_of[i] != k_unset) continue;

            walk.clear();
            std::size_t cursor = i;
            std::size_t root = k_unset;
            while (true)
            {
                walk.push_back(cursor);
                const std::int64_t receiver = basin.receiver_flats[cursor];
                auto rit = (receiver < 0) ? flat_to_index.end()
                                          : flat_to_index.find(receiver);
                if (rit == flat_to_index.end() || keep[rit->second] == 0)
                {
                    root = cursor;
                    break;
                }
                if (root_of[rit->second] != k_unset)
                {
                    root = root_of[rit->second];
                    break;
                }
                cursor = rit->second;
            }
            for (const std::size_t node : walk)
            {
                root_of[node] = root;
            }
        }

        std::vector<std::int64_t> new_flats;
        std::vector<float>        new_chi;
        std::vector<std::int64_t> new_receivers;
        new_flats.reserve(n_nodes);
        new_chi.reserve(n_nodes);
        new_receivers.reserve(n_nodes);

        for (std::size_t i = 0; i < n_nodes; ++i)
        {
            if (keep[i] == 0)
            {
                ++n_removed;
                continue;
            }
            new_flats.push_back(basin.flats[i]);
            new_chi.push_back(basin.chi[i] - basin.chi[root_of[i]]);
            new_receivers.push_back(basin.receiver_flats[i]);
        }

        basin.flats = std::move(new_flats);
        basin.chi = std::move(new_chi);
        basin.receiver_flats = std::move(new_receivers);
        basin.outlet_chi = 0.0f;
    }

    return n_removed;
}

// -----------------------------------------------------------------------------
// parse_pair_uncertainty_method / to_string
// -----------------------------------------------------------------------------
lsdtt_pair_uncertainty_method parse_pair_uncertainty_method(
    const std::string& method_text)
{
    if (method_text == "std_dev")   return lsdtt_pair_uncertainty_method::std_dev;
    if (method_text == "bootstrap") return lsdtt_pair_uncertainty_method::bootstrap;
    if (method_text == "ttest")     return lsdtt_pair_uncertainty_method::ttest;
    return lsdtt_pair_uncertainty_method::std_err;
}

std::string to_string(lsdtt_pair_uncertainty_method method)
{
    switch (method)
    {
    case lsdtt_pair_uncertainty_method::std_dev:   return "std_dev";
    case lsdtt_pair_uncertainty_method::std_err:   return "std_err";
    case lsdtt_pair_uncertainty_method::bootstrap: return "bootstrap";
    case lsdtt_pair_uncertainty_method::ttest:     return "ttest";
    }
    return "std_err";
}

// -----------------------------------------------------------------------------
// compute_pair_summaries
// -----------------------------------------------------------------------------
std::vector<lsdtt_pair_summary_record> compute_pair_summaries(
    const std::vector<lsdtt_divide_segment_record>&      divide_segments,
    const std::vector<lsdtt_channel_head_metric_record>& head_metrics,
    const lsdtt_pair_summary_options&                    options)
{
    // Ordered by pair_id, so the row order of the output never depends on the
    // order the segments or the heads happened to arrive in.
    std::map<std::string, pair_accumulator> by_pair;

    // The pair set comes from the DIVIDES, not from the heads: every extracted
    // divide pair gets a row whether or not any head was matched to it, so the
    // summary table and the divide layer always agree about which pairs exist.
    // Two segments sharing a pair_id (two separate stretches of the same divide)
    // collapse onto one accumulator here and so yield exactly one row.
    for (const lsdtt_divide_segment_record& seg : divide_segments)
    {
        if (seg.pair_id.empty()) continue;

        int unit_a = -1;
        int unit_b = -1;
        if (!parse_pair_id(seg.pair_id, unit_a, unit_b))
        {
            // A pair_id that is not "<unit_a>|<unit_b>" cannot be built by
            // make_unit_pair_id; fall back to the record's own unit ids rather
            // than dropping a divide that exists in the geometry output.
            unit_a = seg.unit_a;
            unit_b = seg.unit_b;
        }
        // Side A is the lower-numbered unit, matching make_unit_pair_id, so the
        // A-minus-B sign convention cannot depend on how the pair was built.
        if (unit_b < unit_a) std::swap(unit_a, unit_b);

        pair_accumulator& acc = by_pair[seg.pair_id];
        acc.unit_a = unit_a;
        acc.unit_b = unit_b;
    }

    for (const lsdtt_channel_head_metric_record& rec : head_metrics)
    {
        // An empty pair_id_nearest means the head matched no divide, so it
        // belongs to no side of any pair.
        if (rec.pair_id_nearest.empty()) continue;

        // Heads only ever attach to a pair the divides already established. A
        // head naming a pair that is absent from divide_segments never invents a
        // row; it cannot arise from compute_nearest_divide_metrics, which reads
        // every pair_id_nearest straight off one of these same segments.
        auto it = by_pair.find(rec.pair_id_nearest);
        if (it == by_pair.end()) continue;

        pair_accumulator& acc = it->second;
        const int unit_a = acc.unit_a;
        const int unit_b = acc.unit_b;

        // A head is only ever counted on the side its own unit_id names. It is
        // never assigned to a side by elimination.
        const bool is_a = (rec.unit_id == unit_a);
        const bool is_b = (rec.unit_id == unit_b);
        if (!is_a && !is_b) continue;

        if (is_a) ++acc.n_heads_a; else ++acc.n_heads_b;

        const double values[4] = {rec.head_elevation_m, rec.gradient_to_divide,
                                  rec.relief_m, rec.chi};
        for (int m = 0; m < 4; ++m)
        {
            // NaN metrics are skipped, not propagated: one head whose nearest
            // divide point fell on nodata must not destroy the side mean.
            if (!std::isfinite(values[m])) continue;
            if (is_a) acc.metrics[m].a.push_back(values[m]);
            else      acc.metrics[m].b.push_back(values[m]);
        }
    }

    std::vector<lsdtt_pair_summary_record> out;
    out.reserve(by_pair.size());

    for (const auto& entry : by_pair)
    {
        const pair_accumulator& acc = entry.second;

        lsdtt_pair_summary_record rec;
        rec.pair_id = entry.first;
        rec.unit_a = acc.unit_a;
        rec.unit_b = acc.unit_b;
        rec.n_heads_a = acc.n_heads_a;
        rec.n_heads_b = acc.n_heads_b;
        rec.uncertainty_method = to_string(options.uncertainty_method);

        double* mean_a[4] = {&rec.mean_head_elevation_a, &rec.mean_gradient_a,
                             &rec.mean_relief_a, &rec.mean_chi_a};
        double* mean_b[4] = {&rec.mean_head_elevation_b, &rec.mean_gradient_b,
                             &rec.mean_relief_b, &rec.mean_chi_b};
        double* delta[4]  = {&rec.delta_head_elevation, &rec.delta_gradient,
                             &rec.delta_relief, &rec.delta_chi};
        std::int64_t* n_a[4] = {&rec.n_head_elevation_a, &rec.n_gradient_a,
                                &rec.n_relief_a, &rec.n_chi_a};
        std::int64_t* n_b[4] = {&rec.n_head_elevation_b, &rec.n_gradient_b,
                                &rec.n_relief_b, &rec.n_chi_b};
        double* uncertainty[4] = {&rec.uncertainty_head_elevation,
                                  &rec.uncertainty_gradient,
                                  &rec.uncertainty_relief,
                                  &rec.uncertainty_chi};

        for (int m = 0; m < 4; ++m)
        {
            const pair_metric_sample& sample = acc.metrics[m];
            *n_a[m] = static_cast<std::int64_t>(sample.a.size());
            *n_b[m] = static_cast<std::int64_t>(sample.b.size());
            *mean_a[m] = sample_mean(sample.a);
            *mean_b[m] = sample_mean(sample.b);
            // A minus B, A being the lower-numbered unit.
            *delta[m] = *mean_a[m] - *mean_b[m];
            *uncertainty[m] = pair_uncertainty(
                sample.a, sample.b, options, rec.pair_id, m);
        }

        out.push_back(std::move(rec));
    }

    return out;
}

} // namespace lsdtt
