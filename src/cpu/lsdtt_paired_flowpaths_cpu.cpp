// lsdtt_paired_flowpaths_cpu.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
//
// See lsdtt_paired_flowpaths_cpu.hpp for the method and its provenance.
// -----------------------------------------------------------------------------

#include "lsdtt_paired_flowpaths_cpu.hpp"

#include <algorithm>
#include <cmath>
#include <set>
#include <stdexcept>
#include <utility>

#include "../core/lsdtt_data_loader.hpp"

namespace lsdtt
{

namespace
{

std::string precise(double v)
{
    if (!std::isfinite(v)) return "NaN";
    std::string s = std::to_string(v);
    return s;
}

double median_of(std::vector<double> v)
{
    if (v.empty()) return std::numeric_limits<double>::quiet_NaN();
    std::sort(v.begin(), v.end());
    const std::size_t n = v.size();
    return (n % 2) ? v[n / 2] : 0.5 * (v[n / 2 - 1] + v[n / 2]);
}

// Normalised asymmetry, 2|x-y|/(x+y), the form Goren & Shelef use for both dL
// and dchi. Undefined when the sum is not positive.
double normalised_asymmetry(double x, double y)
{
    const double s = x + y;
    if (!(s > 0.0) || !std::isfinite(s)) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    return 2.0 * std::abs(x - y) / s;
}

} // namespace


std::int64_t find_common_downstream_node(
    const lsdtt_flow_base& flow,
    std::int64_t           a,
    std::int64_t           b)
{
    if (a < 0 || b < 0) return -1;
    if (a == b) return a;

    const lsdtt_raster_meta& meta = flow.meta();
    const std::int64_t n =
        static_cast<std::int64_t>(meta.n_rows) * static_cast<std::int64_t>(meta.n_cols);

    // Mark every node on a's path to its outlet. A step cap equal to the cell
    // count makes a corrupt receiver graph terminate rather than spin.
    std::vector<std::uint8_t> on_a(static_cast<std::size_t>(n), 0u);

    std::int64_t x = a;
    for (std::int64_t guard = 0; guard <= n; ++guard)
    {
        if (x < 0 || x >= n) break;
        on_a[static_cast<std::size_t>(x)] = 1u;
        const std::int64_t r = flow.receiver_of(x);
        if (r < 0 || r == x) break;          // outlet, or self-receiver
        x = r;
    }

    // Walk b until it lands on a's path.
    std::int64_t y = b;
    for (std::int64_t guard = 0; guard <= n; ++guard)
    {
        if (y < 0 || y >= n) break;
        if (on_a[static_cast<std::size_t>(y)]) return y;
        const std::int64_t r = flow.receiver_of(y);
        if (r < 0 || r == y) break;
        y = r;
    }

    return -1;
}


lsdtt_paired_flowpath_result compute_paired_flowpath_statistics(
    const lsdtt_flow_base&                  flow,
    const std::vector<lsdtt_ridgeline_cpu>& divide_ridgelines,
    const std::vector<float>&               fd_at_channel,
    const std::vector<float>&               chi_at_channel,
    const std::vector<float>&               fd_from_outlet,
    const std::vector<float>&               chi_from_outlet,
    const std::vector<std::uint8_t>&        channel_mask,
    const lsdtt_paired_flowpath_options&    opts)
{
    lsdtt_paired_flowpath_result out;

    const lsdtt_raster_meta& meta = flow.meta();
    const int nrows = meta.n_rows;
    const int ncols = meta.n_cols;
    const std::size_t n =
        static_cast<std::size_t>(nrows) * static_cast<std::size_t>(ncols);

    if (fd_at_channel.size() != n || chi_at_channel.size() != n ||
        fd_from_outlet.size() != n || chi_from_outlet.size() != n ||
        channel_mask.size() != n)
    {
        throw std::invalid_argument(
            "compute_paired_flowpath_statistics: all input rasters must have "
            "flat-index size n_rows * n_cols");
    }

    auto flat_of = [ncols](int r, int c) -> std::int64_t
    {
        return static_cast<std::int64_t>(r) * ncols + c;
    };
    auto in_grid = [nrows, ncols](int r, int c) -> bool
    {
        return r >= 0 && r < nrows && c >= 0 && c < ncols;
    };

    std::set<std::pair<std::int64_t, std::int64_t>> seen_pairs;

    for (const lsdtt_ridgeline_cpu& rl : divide_ridgelines)
    {
        const std::vector<lsdtt_pixel_rc>& nodes = rl.nodes();
        const std::vector<double>&         dist  = rl.distances();
        if (nodes.size() < 2) continue;

        const double total = rl.total_length();

        for (std::size_t k = 0; k + 1 < nodes.size(); ++k)
        {
            ++out.n_edges_examined;

            // Distance along the divide to the NEARER end. Both ends are where
            // the divide meets a junction, so proximity to either is what the
            // threshold is about.
            double d_along = std::numeric_limits<double>::quiet_NaN();
            if (k < dist.size())
            {
                const double d0 = dist[k];
                d_along = std::min(d0, total - d0);
            }

            if (opts.threshold_divide_distance_m > 0.0 &&
                std::isfinite(d_along) &&
                d_along < opts.threshold_divide_distance_m)
            {
                ++out.n_rejected_distance;
                continue;
            }

            // The two DEM pixels flanking this divide-grid edge.
            //
            // Divide-grid node (r, c) is the top-left CORNER of DEM pixel
            // (r, c). So a horizontal edge (r,c)->(r,c+1) runs along the top of
            // DEM row r and separates rows r-1 and r; a vertical edge
            // (r,c)->(r+1,c) runs down the left of DEM column c and separates
            // columns c-1 and c.
            const lsdtt_pixel_rc& p0 = nodes[k];
            const lsdtt_pixel_rc& p1 = nodes[k + 1];

            const int dr = p1.row - p0.row;
            const int dc = p1.col - p0.col;

            int ra, ca, rb, cb;
            if (dr == 0 && (dc == 1 || dc == -1))
            {
                const int c = std::min(p0.col, p1.col);
                ra = p0.row - 1; ca = c;
                rb = p0.row;     cb = c;
            }
            else if (dc == 0 && (dr == 1 || dr == -1))
            {
                const int r = std::min(p0.row, p1.row);
                ra = r; ca = p0.col - 1;
                rb = r; cb = p0.col;
            }
            else
            {
                // Diagonal or wrapped step: not a single pixel edge, so there is
                // no unambiguous flanking pair. Skipped rather than guessed.
                continue;
            }

            if (!in_grid(ra, ca) || !in_grid(rb, cb)) continue;

            const std::int64_t fa = flat_of(ra, ca);
            const std::int64_t fb = flat_of(rb, cb);
            const std::size_t  ia = static_cast<std::size_t>(fa);
            const std::size_t  ib = static_cast<std::size_t>(fb);

            if (!flow.is_valid_cell(fa) || !flow.is_valid_cell(fb)) continue;

            // Both sides must actually reach a channel.
            if (!std::isfinite(fd_at_channel[ia]) ||
                !std::isfinite(fd_at_channel[ib]))
            {
                ++out.n_rejected_no_channel;
                continue;
            }

            // The channel each side drains into.
            std::int64_t ch_a = -1, ch_b = -1;
            {
                std::int64_t x = fa;
                for (std::int64_t g = 0; g <= static_cast<std::int64_t>(n); ++g)
                {
                    if (x < 0) break;
                    if (channel_mask[static_cast<std::size_t>(x)]) { ch_a = x; break; }
                    const std::int64_t r = flow.receiver_of(x);
                    if (r < 0 || r == x) break;
                    x = r;
                }
                std::int64_t y = fb;
                for (std::int64_t g = 0; g <= static_cast<std::int64_t>(n); ++g)
                {
                    if (y < 0) break;
                    if (channel_mask[static_cast<std::size_t>(y)]) { ch_b = y; break; }
                    const std::int64_t r = flow.receiver_of(y);
                    if (r < 0 || r == y) break;
                    y = r;
                }
            }

            if (ch_a < 0 || ch_b < 0 || ch_a == ch_b)
            {
                ++out.n_rejected_no_channel;
                continue;
            }

            if (opts.unique_pairs)
            {
                const std::pair<std::int64_t, std::int64_t> key =
                    (ch_a < ch_b) ? std::make_pair(ch_a, ch_b)
                                  : std::make_pair(ch_b, ch_a);
                if (!seen_pairs.insert(key).second)
                {
                    ++out.n_rejected_duplicate;
                    continue;
                }
            }

            if (opts.only_take_channel_heads)
            {
                // A channel head has no channel pixel anywhere upslope of it.
                //
                // The s_vector answers this directly: upslope_nodes() is a
                // zero-copy slice of the node's whole upslope set in preorder,
                // with the node itself first. Testing the entire upslope set
                // rather than just the donors is not extra work -- channel-ness
                // is monotone downstream (contributing area only grows), so an
                // upslope channel pixel implies a channel DONOR and vice versa
                // -- and it costs nothing, because the scan exits on the first
                // channel found. For a trunk node that is immediately; only a
                // genuine head walks its (by definition sub-threshold) set.
                auto is_head = [&](std::int64_t node) -> bool
                {
                    if (flow.donor_count_of(node) == 0) return true;

                    const lsdtt_flat_index_slice up = flow.upslope_nodes(node);
                    for (std::size_t k = 1; k < up.size(); ++k)
                    {
                        const std::int64_t u = up[k];
                        if (u >= 0 && channel_mask[static_cast<std::size_t>(u)])
                            return false;
                    }
                    return true;
                };
                if (!is_head(ch_a) || !is_head(ch_b))
                {
                    ++out.n_rejected_not_head;
                    continue;
                }
            }

            const std::int64_t cj = find_common_downstream_node(flow, ch_a, ch_b);
            if (cj < 0)
            {
                // The two paths never rejoin: different basins. This is the
                // condition Goren & Shelef enforce with a basin raster; testing
                // the flow paths directly is equivalent and needs no extra input.
                ++out.n_rejected_no_junction;
                continue;
            }

            const std::size_t icj = static_cast<std::size_t>(cj);
            if (!std::isfinite(fd_from_outlet[icj])) continue;

            // Both measures are referred to the common junction, so each path
            // is counted only over the stretch where the two are separate.
            const double fd_cj = static_cast<double>(fd_from_outlet[icj]);
            const double L_a = static_cast<double>(fd_at_channel[ia]) - fd_cj;
            const double L_b = static_cast<double>(fd_at_channel[ib]) - fd_cj;

            const double chi_cj = static_cast<double>(chi_from_outlet[icj]);

            lsdtt_paired_flowpath rec;
            rec.pixel_a = fa;
            rec.pixel_b = fb;
            rec.channel_a = ch_a;
            rec.channel_b = ch_b;
            rec.common_junction = cj;
            rec.divide_distance_m = d_along;
            rec.L_a = L_a;
            rec.L_b = L_b;
            rec.chi_a = static_cast<double>(chi_at_channel[ia]) - chi_cj;
            rec.chi_b = static_cast<double>(chi_at_channel[ib]) - chi_cj;
            rec.delta_L = normalised_asymmetry(L_a, L_b);
            rec.delta_chi = normalised_asymmetry(rec.chi_a, rec.chi_b);

            if (!std::isfinite(rec.delta_L)) continue;

            out.pairs.push_back(rec);
        }
    }

    out.n_pairs = out.pairs.size();

    std::vector<double> dl, dchi;
    dl.reserve(out.pairs.size());
    for (const lsdtt_paired_flowpath& p : out.pairs)
    {
        dl.push_back(p.delta_L);
        if (std::isfinite(p.delta_chi)) dchi.push_back(p.delta_chi);
    }

    if (!dl.empty())
    {
        out.median_delta_L = median_of(dl);
        double s = 0.0;
        for (double v : dl) s += v;
        out.mean_delta_L = s / static_cast<double>(dl.size());
    }
    out.median_delta_chi = median_of(dchi);

    return out;
}


void write_paired_flowpath_data(
    const lsdtt_paired_flowpath_result& result,
    const lsdtt_raster_meta&            meta,
    const std::string&                  out_path)
{
    lsdtt_data_loader loader;

    std::vector<std::string> headers = {
        "pair_number",
        "row_a", "col_a", "row_b", "col_b",
        "divide_distance_m",
        "path_length_a", "path_length_b",
        "chi_a", "chi_b",
        "delta_L", "delta_chi"
    };

    const int ncols = meta.n_cols;

    std::vector<std::vector<std::string>> rows;
    rows.reserve(result.pairs.size());
    for (std::size_t i = 0; i < result.pairs.size(); ++i)
    {
        const lsdtt_paired_flowpath& p = result.pairs[i];
        rows.push_back({
            std::to_string(i),
            std::to_string(static_cast<int>(p.pixel_a / ncols)),
            std::to_string(static_cast<int>(p.pixel_a % ncols)),
            std::to_string(static_cast<int>(p.pixel_b / ncols)),
            std::to_string(static_cast<int>(p.pixel_b % ncols)),
            precise(p.divide_distance_m),
            precise(p.L_a), precise(p.L_b),
            precise(p.chi_a), precise(p.chi_b),
            precise(p.delta_L), precise(p.delta_chi)
        });
    }

    loader.write_csv(out_path, headers, rows);
}


void write_paired_flowpath_summary(
    const lsdtt_paired_flowpath_result& result,
    const std::string&                  out_path)
{
    lsdtt_data_loader loader;

    std::vector<std::string> headers = {
        "n_pairs", "median_delta_L", "mean_delta_L", "median_delta_chi",
        "n_edges_examined",
        "n_rejected_distance", "n_rejected_no_channel",
        "n_rejected_no_junction", "n_rejected_duplicate", "n_rejected_not_head"
    };

    std::vector<std::vector<std::string>> rows = {{
        std::to_string(result.n_pairs),
        precise(result.median_delta_L),
        precise(result.mean_delta_L),
        precise(result.median_delta_chi),
        std::to_string(result.n_edges_examined),
        std::to_string(result.n_rejected_distance),
        std::to_string(result.n_rejected_no_channel),
        std::to_string(result.n_rejected_no_junction),
        std::to_string(result.n_rejected_duplicate),
        std::to_string(result.n_rejected_not_head)
    }};

    loader.write_csv(out_path, headers, rows);
}

} // namespace lsdtt
