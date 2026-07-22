// lsdtt_raster_skeleton_graph_cpu.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html

#include "lsdtt_raster_skeleton_graph_cpu.hpp"

#include "../utils/lsdtt_linear_algebra.hpp"
#include "../utils/lsdtt_raster_geometry.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <queue>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>
#include <set>
#include <map>

namespace lsdtt
{
namespace
{

// Crossing number for pixel (row, col) in a binary skeleton.
//
// Counts the number of 0->1 transitions in the clockwise-ordered 8-neighbourhood.
// This is the standard digital-topology endpoint test:
//   CN = 1  ->  endpoint  (one branch leaving the pixel)
//   CN = 2  ->  chain pixel
//   CN >= 3 ->  junction
//
// Unlike the raw 8-neighbour degree count, CN correctly handles corner pixels
// where two adjacent diagonal neighbours are present (e.g. south + southeast):
// those count as ONE transition, not two branches.
template<typename IsForeground>
inline int crossing_number(int row, int col, IsForeground is_fg)
{
    // Clockwise from east: E, NE, N, NW, W, SW, S, SE
    static constexpr int dr[8] = { 0, -1, -1, -1,  0,  1,  1,  1};
    static constexpr int dc[8] = { 1,  1,  0, -1, -1, -1,  0,  1};

    int cn = 0;
    for (int i = 0; i < 8; ++i)
    {
        int v0 = is_fg(row + dr[i],       col + dc[i])       ? 1 : 0;
        int v1 = is_fg(row + dr[(i+1)%8], col + dc[(i+1)%8]) ? 1 : 0;
        if (v0 == 0 && v1 == 1) ++cn;
    }
    return cn;
}



// -----------------------------------------------------------------------------
// estimate_endpoint_tangent_one_side_
// -----------------------------------------------------------------------------
// Internal helper for one segment end.
// Returns true on success and fills `out`.
// -----------------------------------------------------------------------------
static bool estimate_endpoint_tangent_one_side_(
    const skeleton_segment&       seg,
    std::size_t                   segment_index,
    bool                          at_start,
    double                        fit_length_pixels,
    int                           min_support_pixels,
    skeleton_endpoint_tangent&    out)
{
    if (seg.is_loop) return false;
    if (seg.pixels.empty()) return false;

    const skeleton_node& node =
        at_start ? seg.start_node : seg.end_node;

    if (node.crossing_number != 1) return false;


    const std::size_t n = seg.pixels.size();
    if (n < 2) return false;

    // -------------------------------------------------------------------------
    // Collect support pixels from endpoint inward along the ordered segment.
    // -------------------------------------------------------------------------
    std::vector<skeleton_pixel> support;
    support.reserve(n);

    double acc_len = 0.0;

    auto pixel_at = [&](std::size_t i) -> const skeleton_pixel&
    {
        return at_start ? seg.pixels[i]
                        : seg.pixels[n - 1 - i];
    };

    support.push_back(pixel_at(0));

    for (std::size_t i = 1; i < n; ++i)
    {
        const skeleton_pixel& p0 = pixel_at(i - 1);
        const skeleton_pixel& p1 = pixel_at(i);

        const double dr = static_cast<double>(p1.row - p0.row);
        const double dc = static_cast<double>(p1.col - p0.col);
        acc_len += std::sqrt(dr * dr + dc * dc);

        support.push_back(p1);

        if (acc_len >= fit_length_pixels)
            break;
    }

    if (static_cast<int>(support.size()) < min_support_pixels)
        return false;

    // -------------------------------------------------------------------------
    // Endpoint and inward chord.
    // support[0] is the endpoint; support.back() is the deepest sampled pixel.
    // -------------------------------------------------------------------------
    const double er = static_cast<double>(support.front().row);
    const double ec = static_cast<double>(support.front().col);

    const double chord_r =
        static_cast<double>(support.back().row) - er;
    const double chord_c =
        static_cast<double>(support.back().col) - ec;

    const double chord_norm =
        std::sqrt(chord_r * chord_r + chord_c * chord_c);

    if (!(std::isfinite(chord_norm)) || chord_norm <= 0.0)
        return false;

    // -------------------------------------------------------------------------
    // PCA / covariance orientation of support points.
    // -------------------------------------------------------------------------
    double mean_r = 0.0;
    double mean_c = 0.0;

    for (const auto& p : support)
    {
        mean_r += static_cast<double>(p.row);
        mean_c += static_cast<double>(p.col);
    }

    mean_r /= static_cast<double>(support.size());
    mean_c /= static_cast<double>(support.size());

    double cov_rr = 0.0;
    double cov_cc = 0.0;
    double cov_rc = 0.0;

    for (const auto& p : support)
    {
        const double dr = static_cast<double>(p.row) - mean_r;
        const double dc = static_cast<double>(p.col) - mean_c;
        cov_rr += dr * dr;
        cov_cc += dc * dc;
        cov_rc += dr * dc;
    }

    cov_rr /= static_cast<double>(support.size());
    cov_cc /= static_cast<double>(support.size());
    cov_rc /= static_cast<double>(support.size());

    const double trace = cov_rr + cov_cc;
    const double disc  = std::sqrt(
        std::max(0.0, (cov_rr - cov_cc) * (cov_rr - cov_cc) +
                      4.0 * cov_rc * cov_rc)
    );

    const double lambda_max = 0.5 * (trace + disc);

    double vr = 0.0;
    double vc = 0.0;

    // Eigenvector for the dominant eigenvalue.
    if (std::abs(cov_rc) > 1e-12 ||
        std::abs(lambda_max - cov_rr) > 1e-12)
    {
        vr = cov_rc;
        vc = lambda_max - cov_rr;
    }
    else
    {
        // Diagonal covariance fallback.
        if (cov_rr >= cov_cc)
        {
            vr = 1.0;
            vc = 0.0;
        }
        else
        {
            vr = 0.0;
            vc = 1.0;
        }
    }

    double vnorm = std::sqrt(vr * vr + vc * vc);

    // Degenerate covariance fallback: use the inward chord direction.
    if (!(std::isfinite(vnorm)) || vnorm <= 0.0)
    {
        vr = chord_r / chord_norm;
        vc = chord_c / chord_norm;
        vnorm = 1.0;
    }
    else
    {
        vr /= vnorm;
        vc /= vnorm;
    }

    // -------------------------------------------------------------------------
    // Orient the tangent outward from the endpoint.
    // The chord points inward from endpoint into the segment. We want the
    // opposite direction.
    // -------------------------------------------------------------------------
    const double dot_inward = vr * chord_r + vc * chord_c;
    if (dot_inward > 0.0)
    {
        vr = -vr;
        vc = -vc;
    }

    out.segment_index         = segment_index;
    out.at_start              = at_start;
    out.row                   = node.row;
    out.col                   = node.col;
    out.d_row                 = vr;
    out.d_col                 = vc;
    out.n_support_pixels      = static_cast<int>(support.size());
    out.support_length_pixels = acc_len;

    return true;
}






// -----------------------------------------------------------------------------
// validate_endpoint_pair_request
// -----------------------------------------------------------------------------
static void validate_endpoint_pair_request(
    const endpoint_pair_request& request)
{
    if (!(std::isfinite(request.max_gap_pixels) &&
          request.max_gap_pixels > 0.0))
    {
        throw std::runtime_error(
            "lsdtt::find_skeleton_endpoint_pair_candidates: max_gap_pixels must be finite and > 0"
        );
    }

    if (!(std::isfinite(request.max_facing_angle_degrees) &&
          request.max_facing_angle_degrees > 0.0 &&
          request.max_facing_angle_degrees <= 180.0))
    {
        throw std::runtime_error(
            "lsdtt::find_skeleton_endpoint_pair_candidates: max_facing_angle_degrees must be in (0,180]"
        );
    }

    if (!(std::isfinite(request.max_heading_misalignment_degrees) &&
          request.max_heading_misalignment_degrees >= 0.0 &&
          request.max_heading_misalignment_degrees <= 180.0))
    {
        throw std::runtime_error(
            "lsdtt::find_skeleton_endpoint_pair_candidates: max_heading_misalignment_degrees must be in [0,180]"
        );
    }

    if (!(std::isfinite(request.weight_distance) &&
          std::isfinite(request.weight_facing) &&
          std::isfinite(request.weight_heading)))
    {
        throw std::runtime_error(
            "lsdtt::find_skeleton_endpoint_pair_candidates: all weights must be finite"
        );
    }

    if (request.weight_distance < 0.0 ||
        request.weight_facing   < 0.0 ||
        request.weight_heading  < 0.0)
    {
        throw std::runtime_error(
            "lsdtt::find_skeleton_endpoint_pair_candidates: all weights must be >= 0"
        );
    }

    if (!(std::isfinite(request.min_segment_mean_strength) &&
          request.min_segment_mean_strength >= 0.0))
    {
        throw std::runtime_error(
            "lsdtt::find_skeleton_endpoint_pair_candidates: min_segment_mean_strength must be finite and >= 0"
        );
    }

    const double wsum =
        request.weight_distance +
        request.weight_facing +
        request.weight_heading;

    if (wsum <= 0.0)
    {
        throw std::runtime_error(
            "lsdtt::find_skeleton_endpoint_pair_candidates: at least one weight must be > 0"
        );
    }
}


// -----------------------------------------------------------------------------
// compute_skeleton_segment_component_ids_
//
// Infers connected-component IDs for extracted skeleton segments.
//
// Two segments are adjacent if they share a node coordinate
// (start_node or end_node row/col).
//
// Loop segments that share no node with any other segment each form their own
// singleton component.
// -----------------------------------------------------------------------------
static std::vector<int> compute_skeleton_segment_component_ids_(
    const std::vector<skeleton_segment>& segments)
{
    const std::size_t nseg = segments.size();
    std::vector<int> component_id(nseg, -1);

    using node_key = std::pair<int,int>;
    std::map<node_key, std::vector<std::size_t>> node_to_segments;

    for (std::size_t i = 0; i < nseg; ++i)
    {
        const skeleton_segment& seg = segments[i];
        if (seg.pixels.empty()) continue;

        const node_key k0(seg.start_node.row, seg.start_node.col);
        const node_key k1(seg.end_node.row,   seg.end_node.col);

        node_to_segments[k0].push_back(i);
        if (k1 != k0)
            node_to_segments[k1].push_back(i);
    }

    int next_component_id = 0;

    for (std::size_t seed = 0; seed < nseg; ++seed)
    {
        if (component_id[seed] >= 0) continue;

        std::vector<std::size_t> stack;
        stack.push_back(seed);
        component_id[seed] = next_component_id;

        while (!stack.empty())
        {
            const std::size_t sidx = stack.back();
            stack.pop_back();

            const skeleton_segment& seg = segments[sidx];

            const node_key keys[2] =
            {
                node_key(seg.start_node.row, seg.start_node.col),
                node_key(seg.end_node.row,   seg.end_node.col)
            };

            for (const node_key& key : keys)
            {
                const auto it = node_to_segments.find(key);
                if (it == node_to_segments.end()) continue;

                for (const std::size_t nbr : it->second)
                {
                    if (component_id[nbr] >= 0) continue;
                    component_id[nbr] = next_component_id;
                    stack.push_back(nbr);
                }
            }
        }

        ++next_component_id;
    }

    return component_id;
}


// -----------------------------------------------------------------------------
// validate_endpoint_pair_selection_request
// -----------------------------------------------------------------------------
static void validate_endpoint_pair_selection_request(
    const endpoint_pair_selection_request& request)
{
    if (!std::isfinite(request.max_score) &&
        request.max_score != std::numeric_limits<double>::infinity())
    {
        throw std::runtime_error(
            "lsdtt::select_skeleton_endpoint_pairs: max_score must be finite or +infinity"
        );
    }
}

// -----------------------------------------------------------------------------
// endpoint_pair_candidate_less_
//
// Deterministic total ordering for candidates.
// Lower score is better.
// -----------------------------------------------------------------------------
static bool endpoint_pair_candidate_less_(
    const skeleton_endpoint_pair_candidate& lhs,
    const skeleton_endpoint_pair_candidate& rhs)
{
    if (lhs.score != rhs.score)
        return lhs.score < rhs.score;

    if (lhs.gap_distance_pixels != rhs.gap_distance_pixels)
        return lhs.gap_distance_pixels < rhs.gap_distance_pixels;

    if (lhs.tangent_index_a != rhs.tangent_index_a)
        return lhs.tangent_index_a < rhs.tangent_index_a;

    return lhs.tangent_index_b < rhs.tangent_index_b;
}





// -----------------------------------------------------------------------------
// validate_skeleton_bridge_rasterise_request
// -----------------------------------------------------------------------------
static void validate_skeleton_bridge_rasterise_request(
    const skeleton_bridge_rasterise_request& request)
{
    (void)request;
}

// -----------------------------------------------------------------------------
// rasterise_line_bresenham_
//
// Returns all raster pixels on the integer line from (r0,c0) to (r1,c1),
// including both endpoints.
// -----------------------------------------------------------------------------
static std::vector<std::pair<int,int>> rasterise_line_bresenham_(
    int r0, int c0, int r1, int c1)
{
    std::vector<std::pair<int,int>> pixels;

    int x0 = c0;
    int y0 = r0;
    int x1 = c1;
    int y1 = r1;

    const int dx = std::abs(x1 - x0);
    const int sx = (x0 < x1) ? 1 : -1;
    const int dy = -std::abs(y1 - y0);
    const int sy = (y0 < y1) ? 1 : -1;
    int err = dx + dy;

    while (true)
    {
        pixels.emplace_back(y0, x0);

        if (x0 == x1 && y0 == y1)
            break;

        const int e2 = 2 * err;

        if (e2 >= dy)
        {
            err += dy;
            x0  += sx;
        }

        if (e2 <= dx)
        {
            err += dx;
            y0  += sy;
        }
    }

    return pixels;
}



// -----------------------------------------------------------------------------
// validate_skeleton_bridge_route_request
// -----------------------------------------------------------------------------
static void validate_skeleton_bridge_route_request(
    const skeleton_bridge_route_request& request)
{
    if (!(std::isfinite(request.corridor_half_width_pixels) &&
          request.corridor_half_width_pixels > 0.0))
    {
        throw std::runtime_error(
            "lsdtt::route_skeleton_endpoint_bridges_least_cost: corridor_half_width_pixels must be finite and > 0"
        );
    }

    if (!(std::isfinite(request.base_step_cost) &&
          request.base_step_cost > 0.0))
    {
        throw std::runtime_error(
            "lsdtt::route_skeleton_endpoint_bridges_least_cost: base_step_cost must be finite and > 0"
        );
    }

    if (!(std::isfinite(request.support_penalty_weight) &&
          request.support_penalty_weight >= 0.0))
    {
        throw std::runtime_error(
            "lsdtt::route_skeleton_endpoint_bridges_least_cost: support_penalty_weight must be finite and >= 0"
        );
    }

    if (!(std::isfinite(request.foreground_cost_scale) &&
          request.foreground_cost_scale > 0.0))
    {
        throw std::runtime_error(
            "lsdtt::route_skeleton_endpoint_bridges_least_cost: foreground_cost_scale must be finite and > 0"
        );
    }
}

} // end unnamed namespace


// -----------------------------------------------------------------------------
// thin_zhang_suen
// -----------------------------------------------------------------------------
lsdtt_raster_cpu thin_zhang_suen(
    const lsdtt_raster_cpu& binary)
{
    const lsdtt_raster_meta& meta = binary.get_meta();

    if (!meta.is_valid())
    {
        throw std::runtime_error(
            "lsdtt::thin_zhang_suen: input raster metadata is invalid"
        );
    }

    const int nrows = meta.n_rows;
    const int ncols = meta.n_cols;

    const std::vector<float>&   in_data = binary.get_data();
    const std::vector<uint8_t>& in_mask = binary.get_mask();

    const std::size_t npix =
        static_cast<std::size_t>(nrows) * static_cast<std::size_t>(ncols);

    std::vector<uint8_t> valid_domain(npix, uint8_t{0});
    std::vector<uint8_t> work(npix, uint8_t{0});

    for (std::size_t i = 0; i < npix; ++i)
    {
        if (in_mask[i] == 0) continue;
        if (meta.is_nodata(in_data[i])) continue;
        if (!std::isfinite(in_data[i])) continue;

        valid_domain[i] = uint8_t{1};
        work[i] = (in_data[i] > 0.5f) ? uint8_t{1} : uint8_t{0};
    }

    // Nothing to do on tiny rasters
    if (nrows < 3 || ncols < 3)
    {
        std::vector<float>   out_data(npix, 0.0f);
        std::vector<uint8_t> out_mask(npix, uint8_t{0});

        for (std::size_t i = 0; i < npix; ++i)
        {
            if (valid_domain[i] != 0)
            {
                out_mask[i] = uint8_t{1};
                out_data[i] = (work[i] != 0) ? 1.0f : 0.0f;
            }
        }

        return lsdtt_raster_cpu(std::move(out_data), std::move(out_mask), meta);
    }

    auto idx_of = [ncols](int row, int col) -> std::size_t
    {
        return static_cast<std::size_t>(row) * static_cast<std::size_t>(ncols) +
               static_cast<std::size_t>(col);
    };

    auto pixel = [&](int row, int col) -> uint8_t
    {
        if (row < 0 || row >= nrows || col < 0 || col >= ncols)
            return uint8_t{0};

        const std::size_t idx = idx_of(row, col);
        if (valid_domain[idx] == 0)
            return uint8_t{0};

        return work[idx];
    };

    auto neighbour_count = [&](int row, int col) -> int
    {
        return static_cast<int>(pixel(row - 1, col    )) +  // p2
               static_cast<int>(pixel(row - 1, col + 1)) +  // p3
               static_cast<int>(pixel(row,     col + 1)) +  // p4
               static_cast<int>(pixel(row + 1, col + 1)) +  // p5
               static_cast<int>(pixel(row + 1, col    )) +  // p6
               static_cast<int>(pixel(row + 1, col - 1)) +  // p7
               static_cast<int>(pixel(row,     col - 1)) +  // p8
               static_cast<int>(pixel(row - 1, col - 1));   // p9
    };

    auto transition_count = [&](int row, int col) -> int
    {
        const uint8_t p2 = pixel(row - 1, col    );
        const uint8_t p3 = pixel(row - 1, col + 1);
        const uint8_t p4 = pixel(row,     col + 1);
        const uint8_t p5 = pixel(row + 1, col + 1);
        const uint8_t p6 = pixel(row + 1, col    );
        const uint8_t p7 = pixel(row + 1, col - 1);
        const uint8_t p8 = pixel(row,     col - 1);
        const uint8_t p9 = pixel(row - 1, col - 1);

        int a = 0;
        if (p2 == 0 && p3 != 0) ++a;
        if (p3 == 0 && p4 != 0) ++a;
        if (p4 == 0 && p5 != 0) ++a;
        if (p5 == 0 && p6 != 0) ++a;
        if (p6 == 0 && p7 != 0) ++a;
        if (p7 == 0 && p8 != 0) ++a;
        if (p8 == 0 && p9 != 0) ++a;
        if (p9 == 0 && p2 != 0) ++a;
        return a;
    };

    std::vector<std::size_t> to_delete;
    bool changed = true;

    while (changed)
    {
        changed = false;
        to_delete.clear();

        // Sub-iteration 1
        for (int row = 1; row < nrows - 1; ++row)
        {
            for (int col = 1; col < ncols - 1; ++col)
            {
                const std::size_t idx = idx_of(row, col);
                if (valid_domain[idx] == 0 || work[idx] == 0) continue;

                const uint8_t p2 = pixel(row - 1, col    );
                const uint8_t p4 = pixel(row,     col + 1);
                const uint8_t p6 = pixel(row + 1, col    );
                const uint8_t p8 = pixel(row,     col - 1);

                const int b = neighbour_count(row, col);
                const int a = transition_count(row, col);

                if (b < 2 || b > 6) continue;
                if (a != 1)         continue;
                if (p2 != 0 && p4 != 0 && p6 != 0) continue;
                if (p4 != 0 && p6 != 0 && p8 != 0) continue;

                to_delete.push_back(idx);
            }
        }

        if (!to_delete.empty())
        {
            changed = true;
            for (std::size_t idx : to_delete) work[idx] = uint8_t{0};
        }

        to_delete.clear();

        // Sub-iteration 2
        for (int row = 1; row < nrows - 1; ++row)
        {
            for (int col = 1; col < ncols - 1; ++col)
            {
                const std::size_t idx = idx_of(row, col);
                if (valid_domain[idx] == 0 || work[idx] == 0) continue;

                const uint8_t p2 = pixel(row - 1, col    );
                const uint8_t p4 = pixel(row,     col + 1);
                const uint8_t p6 = pixel(row + 1, col    );
                const uint8_t p8 = pixel(row,     col - 1);

                const int b = neighbour_count(row, col);
                const int a = transition_count(row, col);

                if (b < 2 || b > 6) continue;
                if (a != 1)         continue;
                if (p2 != 0 && p4 != 0 && p8 != 0) continue;
                if (p2 != 0 && p6 != 0 && p8 != 0) continue;

                to_delete.push_back(idx);
            }
        }

        if (!to_delete.empty())
        {
            changed = true;
            for (std::size_t idx : to_delete) work[idx] = uint8_t{0};
        }
    }

    std::vector<float>   out_data(npix, 0.0f);
    std::vector<uint8_t> out_mask(npix, uint8_t{0});

    for (std::size_t i = 0; i < npix; ++i)
    {
        if (valid_domain[i] != 0)
        {
            out_mask[i] = uint8_t{1};
            out_data[i] = (work[i] != 0) ? 1.0f : 0.0f;
        }
    }

    return lsdtt_raster_cpu(std::move(out_data), std::move(out_mask), meta);
}






// -----------------------------------------------------------------------------
// thin_distance_guided
// -----------------------------------------------------------------------------
lsdtt_raster_cpu thin_distance_guided(
    const lsdtt_raster_cpu& binary,
    double                  preserve_distance_pixels)
{
    const lsdtt_raster_meta& meta = binary.get_meta();

    if (!meta.is_valid())
    {
        throw std::runtime_error(
            "lsdtt::thin_distance_guided: input raster metadata is invalid"
        );
    }

    if (!std::isfinite(preserve_distance_pixels) ||
        preserve_distance_pixels < 0.0)
    {
        throw std::runtime_error(
            "lsdtt::thin_distance_guided: preserve_distance_pixels must be finite and >= 0"
        );
    }

    const int nrows = meta.n_rows;
    const int ncols = meta.n_cols;

    const std::vector<float>&   in_data = binary.get_data();
    const std::vector<uint8_t>& in_mask = binary.get_mask();

    const std::size_t npix =
        static_cast<std::size_t>(nrows) * static_cast<std::size_t>(ncols);

    std::vector<uint8_t> valid_domain(npix, uint8_t{0});
    std::vector<uint8_t> work(npix, uint8_t{0});

    for (std::size_t i = 0; i < npix; ++i)
    {
        if (in_mask[i] == 0) continue;
        if (meta.is_nodata(in_data[i])) continue;
        if (!std::isfinite(in_data[i])) continue;

        valid_domain[i] = uint8_t{1};
        work[i] = (in_data[i] > 0.5f) ? uint8_t{1} : uint8_t{0};
    }

    // Small rasters: just return the valid-domain binary mask unchanged.
    if (nrows < 3 || ncols < 3)
    {
        std::vector<float>   out_data(npix, 0.0f);
        std::vector<uint8_t> out_mask(npix, uint8_t{0});

        for (std::size_t i = 0; i < npix; ++i)
        {
            if (valid_domain[i] != 0)
            {
                out_mask[i] = uint8_t{1};
                out_data[i] = (work[i] != 0) ? 1.0f : 0.0f;
            }
        }

        return lsdtt_raster_cpu(std::move(out_data), std::move(out_mask), meta);
    }

    auto idx_of = [ncols](int row, int col) -> std::size_t
    {
        return static_cast<std::size_t>(row) * static_cast<std::size_t>(ncols) +
               static_cast<std::size_t>(col);
    };

    auto in_bounds = [nrows, ncols](int row, int col) -> bool
    {
        return row >= 0 && row < nrows && col >= 0 && col < ncols;
    };

    auto pixel = [&](int row, int col) -> uint8_t
    {
        if (row < 0 || row >= nrows || col < 0 || col >= ncols)
            return uint8_t{0};

        const std::size_t idx = idx_of(row, col);
        if (valid_domain[idx] == 0)
            return uint8_t{0};

        return work[idx];
    };

    // -------------------------------------------------------------------------
    // Step 1. Build an approximate Euclidean support field inside the binary
    // mask: distance from each valid pixel to the nearest valid background
    // pixel. We use an 8-neighbour Dijkstra metric with axial cost 1 and
    // diagonal cost sqrt(2).
    //
    // This is not an exact continuous EDT, but it is a robust centre-support
    // field that is adequate for guiding thinning.
    // -------------------------------------------------------------------------
    constexpr int DR8[8] = {-1,-1, 0, 1, 1, 1, 0,-1};
    constexpr int DC8[8] = { 0, 1, 1, 1, 0,-1,-1,-1};

    const double sqrt2 = std::sqrt(2.0);
    const double step_cost[8] =
        {1.0, sqrt2, 1.0, sqrt2, 1.0, sqrt2, 1.0, sqrt2};

    struct qnode
    {
        double      d   = std::numeric_limits<double>::infinity();
        std::size_t idx = 0;
    };

    struct qnode_greater
    {
        bool operator()(const qnode& a, const qnode& b) const
        {
            if (a.d != b.d) return a.d > b.d;
            return a.idx > b.idx;
        }
    };

    std::vector<double> dist(npix, std::numeric_limits<double>::infinity());

    std::priority_queue<
        qnode,
        std::vector<qnode>,
        qnode_greater> pq;

    int n_background_seed = 0;

    for (int row = 0; row < nrows; ++row)
    {
        for (int col = 0; col < ncols; ++col)
        {
            const std::size_t idx = idx_of(row, col);
            if (valid_domain[idx] == 0) continue;

            if (work[idx] == 0)
            {
                dist[idx] = 0.0;
                pq.push({0.0, idx});
                ++n_background_seed;
            }
        }
    }

    // Degenerate case: no valid background exists inside the valid domain.
    // In that case, fall back to the standard thinner.
    if (n_background_seed == 0)
    {
        return thin_zhang_suen(binary);
    }

    while (!pq.empty())
    {
        const qnode cur = pq.top();
        pq.pop();

        if (cur.d != dist[cur.idx])
            continue;

        const int row =
            static_cast<int>(cur.idx / static_cast<std::size_t>(ncols));
        const int col =
            static_cast<int>(cur.idx % static_cast<std::size_t>(ncols));

        for (int k = 0; k < 8; ++k)
        {
            const int rr = row + DR8[k];
            const int cc = col + DC8[k];

            if (!in_bounds(rr, cc)) continue;

            const std::size_t nidx = idx_of(rr, cc);
            if (valid_domain[nidx] == 0) continue;

            const double cand = cur.d + step_cost[k];
            if (cand < dist[nidx])
            {
                dist[nidx] = cand;
                pq.push({cand, nidx});
            }
        }
    }

    auto neighbour_count = [&](int row, int col) -> int
    {
        return static_cast<int>(pixel(row - 1, col    )) +  // p2
               static_cast<int>(pixel(row - 1, col + 1)) +  // p3
               static_cast<int>(pixel(row,     col + 1)) +  // p4
               static_cast<int>(pixel(row + 1, col + 1)) +  // p5
               static_cast<int>(pixel(row + 1, col    )) +  // p6
               static_cast<int>(pixel(row + 1, col - 1)) +  // p7
               static_cast<int>(pixel(row,     col - 1)) +  // p8
               static_cast<int>(pixel(row - 1, col - 1));   // p9
    };

    auto transition_count = [&](int row, int col) -> int
    {
        const uint8_t p2 = pixel(row - 1, col    );
        const uint8_t p3 = pixel(row - 1, col + 1);
        const uint8_t p4 = pixel(row,     col + 1);
        const uint8_t p5 = pixel(row + 1, col + 1);
        const uint8_t p6 = pixel(row + 1, col    );
        const uint8_t p7 = pixel(row + 1, col - 1);
        const uint8_t p8 = pixel(row,     col - 1);
        const uint8_t p9 = pixel(row - 1, col - 1);

        int a = 0;
        if (p2 == 0 && p3 != 0) ++a;
        if (p3 == 0 && p4 != 0) ++a;
        if (p4 == 0 && p5 != 0) ++a;
        if (p5 == 0 && p6 != 0) ++a;
        if (p6 == 0 && p7 != 0) ++a;
        if (p7 == 0 && p8 != 0) ++a;
        if (p8 == 0 && p9 != 0) ++a;
        if (p9 == 0 && p2 != 0) ++a;
        return a;
    };

    // -------------------------------------------------------------------------
    // Protect locally central pixels.
    //
    // A foreground pixel is protected from deletion if:
    //   - it has enough support distance, AND
    //   - among current foreground neighbours, it is a local distance maximum.
    //
    // This keeps the thinning biased toward the middle of the mask, especially
    // across corners and width transitions.
    // -------------------------------------------------------------------------
    auto is_protected_central_pixel =
        [&](int row, int col) -> bool
    {
        if (pixel(row, col) == 0)
            return false;

        const std::size_t idx = idx_of(row, col);
        const double d0 = dist[idx];

        if (!std::isfinite(d0) || d0 < preserve_distance_pixels)
            return false;

        bool has_fg_neighbour = false;
        bool ge_all = true;
        bool gt_any = false;

        for (int k = 0; k < 8; ++k)
        {
            const int rr = row + DR8[k];
            const int cc = col + DC8[k];
            if (pixel(rr, cc) == 0)
                continue;

            has_fg_neighbour = true;

            const double dn = dist[idx_of(rr, cc)];
            if (dn > d0 + 1.0e-12)
                ge_all = false;

            if (d0 > dn + 1.0e-12)
                gt_any = true;
        }

        return has_fg_neighbour && ge_all && gt_any;
    };

    std::vector<std::size_t> to_delete;
    bool changed = true;

    while (changed)
    {
        changed = false;
        to_delete.clear();

        // ---------------------------------------------------------------------
        // Sub-iteration 1
        // Zhang-Suen deletion conditions, but skip pixels that are locally
        // central according to the distance-support field.
        // ---------------------------------------------------------------------
        for (int row = 1; row < nrows - 1; ++row)
        {
            for (int col = 1; col < ncols - 1; ++col)
            {
                const std::size_t idx = idx_of(row, col);
                if (valid_domain[idx] == 0 || work[idx] == 0) continue;
                if (is_protected_central_pixel(row, col)) continue;

                const uint8_t p2 = pixel(row - 1, col    );
                const uint8_t p4 = pixel(row,     col + 1);
                const uint8_t p6 = pixel(row + 1, col    );
                const uint8_t p8 = pixel(row,     col - 1);

                const int b = neighbour_count(row, col);
                const int a = transition_count(row, col);

                if (b < 2 || b > 6) continue;
                if (a != 1)         continue;
                if (p2 != 0 && p4 != 0 && p6 != 0) continue;
                if (p4 != 0 && p6 != 0 && p8 != 0) continue;

                to_delete.push_back(idx);
            }
        }

        if (!to_delete.empty())
        {
            changed = true;
            for (std::size_t idx : to_delete)
                work[idx] = uint8_t{0};
        }

        to_delete.clear();

        // ---------------------------------------------------------------------
        // Sub-iteration 2
        // ---------------------------------------------------------------------
        for (int row = 1; row < nrows - 1; ++row)
        {
            for (int col = 1; col < ncols - 1; ++col)
            {
                const std::size_t idx = idx_of(row, col);
                if (valid_domain[idx] == 0 || work[idx] == 0) continue;
                if (is_protected_central_pixel(row, col)) continue;

                const uint8_t p2 = pixel(row - 1, col    );
                const uint8_t p4 = pixel(row,     col + 1);
                const uint8_t p6 = pixel(row + 1, col    );
                const uint8_t p8 = pixel(row,     col - 1);

                const int b = neighbour_count(row, col);
                const int a = transition_count(row, col);

                if (b < 2 || b > 6) continue;
                if (a != 1)         continue;
                if (p2 != 0 && p4 != 0 && p8 != 0) continue;
                if (p2 != 0 && p6 != 0 && p8 != 0) continue;

                to_delete.push_back(idx);
            }
        }

        if (!to_delete.empty())
        {
            changed = true;
            for (std::size_t idx : to_delete)
                work[idx] = uint8_t{0};
        }
    }

    std::vector<float>   out_data(npix, 0.0f);
    std::vector<uint8_t> out_mask(npix, uint8_t{0});

    for (std::size_t i = 0; i < npix; ++i)
    {
        if (valid_domain[i] != 0)
        {
            out_mask[i] = uint8_t{1};
            out_data[i] = (work[i] != 0) ? 1.0f : 0.0f;
        }
    }

    return lsdtt_raster_cpu(
        std::move(out_data),
        std::move(out_mask),
        meta
    );
}







// -----------------------------------------------------------------------------
// classify_skeleton_nodes
// -----------------------------------------------------------------------------
skeleton_node_result classify_skeleton_nodes(
    const lsdtt_raster_cpu& skeleton)
{
    const lsdtt_raster_meta& meta = skeleton.get_meta();

    if (!meta.is_valid())
    {
        throw std::runtime_error(
            "lsdtt::classify_skeleton_nodes: input raster metadata is invalid"
        );
    }

    const int nrows = meta.n_rows;
    const int ncols = meta.n_cols;
    const std::size_t npix =
        static_cast<std::size_t>(nrows) * static_cast<std::size_t>(ncols);

    const std::vector<float>&   in_data = skeleton.get_data();
    const std::vector<uint8_t>& in_mask = skeleton.get_mask();

    const float nodata = static_cast<float>(meta.no_data_value);

    std::vector<float>   degree_data(npix, nodata);
    std::vector<uint8_t> degree_mask(npix, uint8_t{0});

    std::vector<float>   endpoint_data(npix, nodata);
    std::vector<uint8_t> endpoint_mask(npix, uint8_t{0});

    std::vector<float>   junction_data(npix, nodata);
    std::vector<uint8_t> junction_mask(npix, uint8_t{0});

    std::vector<skeleton_node> endpoints;
    std::vector<skeleton_node> junctions;

    auto is_valid = [&](int row, int col) -> bool
    {
        if (row < 0 || row >= nrows || col < 0 || col >= ncols)
            return false;

        const std::size_t idx =
            static_cast<std::size_t>(row) * static_cast<std::size_t>(ncols) +
            static_cast<std::size_t>(col);

        return in_mask[idx] != 0 &&
               !meta.is_nodata(in_data[idx]) &&
               std::isfinite(in_data[idx]);
    };

    auto is_foreground = [&](int row, int col) -> bool
    {
        if (!is_valid(row, col)) return false;

        const std::size_t idx =
            static_cast<std::size_t>(row) * static_cast<std::size_t>(ncols) +
            static_cast<std::size_t>(col);

        return in_data[idx] > 0.5f;
    };

    for (int row = 0; row < nrows; ++row)
    {
        for (int col = 0; col < ncols; ++col)
        {
            const std::size_t idx =
                static_cast<std::size_t>(row) * static_cast<std::size_t>(ncols) +
                static_cast<std::size_t>(col);

            if (!is_valid(row, col))
                continue;

            degree_mask[idx]   = uint8_t{1};
            endpoint_mask[idx] = uint8_t{1};
            junction_mask[idx] = uint8_t{1};

            degree_data[idx]   = 0.0f;
            endpoint_data[idx] = 0.0f;
            junction_data[idx] = 0.0f;

            if (!is_foreground(row, col))
                continue;

            int degree = 0;

            for (int dr = -1; dr <= 1; ++dr)
            {
                for (int dc = -1; dc <= 1; ++dc)
                {
                    if (dr == 0 && dc == 0) continue;
                    if (is_foreground(row + dr, col + dc))
                        ++degree;
                }
            }

            degree_data[idx] = static_cast<float>(degree);

            const int cn = crossing_number(row, col, is_foreground);

            if (cn == 1)
            {
                endpoint_data[idx] = 1.0f;
                endpoints.push_back({row, col, degree, cn});
            }
            else if (cn >= 3)
            {
                junction_data[idx] = 1.0f;
                junctions.push_back({row, col, degree, cn});
            }
        }
    }

    skeleton_node_result result;
    result.degree_raster = lsdtt_raster_cpu(
        std::move(degree_data), std::move(degree_mask), meta);
    result.endpoint_mask = lsdtt_raster_cpu(
        std::move(endpoint_data), std::move(endpoint_mask), meta);
    result.junction_mask = lsdtt_raster_cpu(
        std::move(junction_data), std::move(junction_mask), meta);
    result.endpoints = std::move(endpoints);
    result.junctions = std::move(junctions);

    return result;
}






















// -----------------------------------------------------------------------------
// prune_skeleton_spurs
// -----------------------------------------------------------------------------
lsdtt_raster_cpu prune_skeleton_spurs(
    const lsdtt_raster_cpu& skeleton,
    int                     max_spur_pixels,
    int                     max_iterations)
{
    const lsdtt_raster_meta& meta = skeleton.get_meta();

    if (!meta.is_valid())
    {
        throw std::runtime_error(
            "lsdtt::prune_skeleton_spurs: input raster metadata is invalid"
        );
    }

    if (max_spur_pixels < 1)
    {
        throw std::runtime_error(
            "lsdtt::prune_skeleton_spurs: max_spur_pixels must be >= 1"
        );
    }

    if (max_iterations < 1)
    {
        throw std::runtime_error(
            "lsdtt::prune_skeleton_spurs: max_iterations must be >= 1"
        );
    }

    const int nrows = meta.n_rows;
    const int ncols = meta.n_cols;
    const std::size_t npix =
        static_cast<std::size_t>(nrows) * static_cast<std::size_t>(ncols);

    const std::vector<float>&   in_data = skeleton.get_data();
    const std::vector<uint8_t>& in_mask = skeleton.get_mask();

    const float nodata = static_cast<float>(meta.no_data_value);

    // -------------------------------------------------------------------------
    // Build working foreground plane over the DEM-valid domain.
    // -------------------------------------------------------------------------
    std::vector<uint8_t> valid(npix, uint8_t{0});
    std::vector<uint8_t> work (npix, uint8_t{0});

    for (std::size_t i = 0; i < npix; ++i)
    {
        if (in_mask[i] == 0) continue;
        if (meta.is_nodata(in_data[i])) continue;
        if (!std::isfinite(in_data[i])) continue;

        valid[i] = uint8_t{1};
        work[i]  = (in_data[i] > 0.5f) ? uint8_t{1} : uint8_t{0};
    }

    auto idx_of = [ncols](int row, int col) -> std::size_t
    {
        return static_cast<std::size_t>(row) * static_cast<std::size_t>(ncols) +
               static_cast<std::size_t>(col);
    };

    auto in_bounds = [nrows, ncols](int row, int col) -> bool
    {
        return row >= 0 && row < nrows && col >= 0 && col < ncols;
    };

    auto is_foreground = [&](int row, int col) -> bool
    {
        if (!in_bounds(row, col))
            return false;

        const std::size_t idx = idx_of(row, col);
        return valid[idx] != 0 && work[idx] != 0;
    };

    auto cn_of = [&](std::size_t idx) -> int
    {
        const int row = static_cast<int>(idx / static_cast<std::size_t>(ncols));
        const int col = static_cast<int>(idx % static_cast<std::size_t>(ncols));
        return crossing_number(row, col, is_foreground);
    };

    auto are_8_neighbours = [&](std::size_t a, std::size_t b) -> bool
    {
        const int ar = static_cast<int>(a / static_cast<std::size_t>(ncols));
        const int ac = static_cast<int>(a % static_cast<std::size_t>(ncols));
        const int br = static_cast<int>(b / static_cast<std::size_t>(ncols));
        const int bc = static_cast<int>(b % static_cast<std::size_t>(ncols));

        return std::max(std::abs(ar - br), std::abs(ac - bc)) == 1;
    };




    constexpr int DR8[8] = {-1,-1,-1, 0, 0, 1, 1, 1};
    constexpr int DC8[8] = {-1, 0, 1,-1, 1,-1, 0, 1};

    // -------------------------------------------------------------------------
    // Iterative pruning.
    // -------------------------------------------------------------------------
    for (int iter = 0; iter < max_iterations; ++iter)
    {
        // ---------------------------------------------------------------------
        // Recompute 8-neighbour degree on the current working skeleton.
        // ---------------------------------------------------------------------
        std::vector<int> degree(npix, 0);

        for (int row = 0; row < nrows; ++row)
        {
            for (int col = 0; col < ncols; ++col)
            {
                const std::size_t idx = idx_of(row, col);
                if (!valid[idx] || !work[idx]) continue;

                int deg = 0;
                for (int d = 0; d < 8; ++d)
                {
                    const int rr = row + DR8[d];
                    const int cc = col + DC8[d];
                    if (!in_bounds(rr, cc)) continue;

                    const std::size_t nidx = idx_of(rr, cc);
                    if (valid[nidx] && work[nidx]) ++deg;
                }
                degree[idx] = deg;
            }
        }

        // traced marks skeleton pixels already assigned to a traced endpoint-rooted
        // path in this iteration, preventing double-processing of endpoint-endpoint
        // segments from both ends.
        std::vector<uint8_t> traced(npix, uint8_t{0});

        // delete_mark accumulates all pixels to be removed in this iteration.
        std::vector<uint8_t> delete_mark(npix, uint8_t{0});
        bool any_deleted_this_iter = false;

        enum class trace_terminal
        {
            hit_junction,
            hit_endpoint,
            ambiguous,
            dead_end
        };

        for (int row0 = 0; row0 < nrows; ++row0)
        {
            for (int col0 = 0; col0 < ncols; ++col0)
            {
                const std::size_t seed = idx_of(row0, col0);

                if (!valid[seed] || !work[seed]) continue;
                if (cn_of(seed) != 1) continue;
                if (traced[seed] != 0) continue;

                std::vector<std::size_t> removable_pixels;

                std::size_t prev = static_cast<std::size_t>(-1);
                std::size_t curr = seed;

                trace_terminal terminal = trace_terminal::dead_end;

                while (true)
                {
                    traced[curr] = uint8_t{1};
                    removable_pixels.push_back(curr);

                    const int crow = static_cast<int>(curr / static_cast<std::size_t>(ncols));
                    const int ccol = static_cast<int>(curr % static_cast<std::size_t>(ncols));
                    const int ccn  = cn_of(curr);

                    if (curr != seed && ccn == 1)
                    {
                        terminal = trace_terminal::hit_endpoint;
                        break;
                    }

                    if (ccn >= 3)
                    {
                        // True topological junction.
                        removable_pixels.pop_back();
                        terminal = trace_terminal::hit_junction;
                        break;
                    }

                    std::vector<std::size_t> next_pixels;
                    next_pixels.reserve(8);

                    for (int d = 0; d < 8; ++d)
                    {
                        const int rr = crow + DR8[d];
                        const int cc = ccol + DC8[d];
                        if (!in_bounds(rr, cc)) continue;

                        const std::size_t nidx = idx_of(rr, cc);
                        if (!valid[nidx] || !work[nidx]) continue;
                        if (nidx == prev) continue;

                        next_pixels.push_back(nidx);
                    }

                    if (next_pixels.empty())
                    {
                        terminal = trace_terminal::dead_end;
                        break;
                    }

                    // Group forward candidates into 8-connected clusters.
                    // Multiple clusters means we have reached a true branch
                    // attachment / junction from the perspective of this trace.
                    std::vector<int> group(next_pixels.size(), -1);
                    int n_groups = 0;

                    for (std::size_t i = 0; i < next_pixels.size(); ++i)
                    {
                        if (group[i] >= 0) continue;

                        group[i] = n_groups;
                        std::vector<std::size_t> stack;
                        stack.push_back(i);

                        while (!stack.empty())
                        {
                            const std::size_t k = stack.back();
                            stack.pop_back();

                            for (std::size_t j = 0; j < next_pixels.size(); ++j)
                            {
                                if (group[j] >= 0) continue;
                                if (!are_8_neighbours(next_pixels[k], next_pixels[j])) continue;

                                group[j] = n_groups;
                                stack.push_back(j);
                            }
                        }

                        ++n_groups;
                    }

                    if (n_groups > 1)
                    {
                        // Preserve the current attachment pixel and prune only
                        // the endpoint-rooted prefix leading into it.
                        removable_pixels.pop_back();
                        terminal = trace_terminal::hit_junction;
                        break;
                    }

                    // Single forward cluster: continue through the nearest
                    // candidate pixel. This suppresses diagonal look-ahead on
                    // staircase lines.
                    std::size_t next = next_pixels.front();
                    int best_step2 = std::numeric_limits<int>::max();

                    for (const std::size_t cand : next_pixels)
                    {
                        const int nr = static_cast<int>(cand / static_cast<std::size_t>(ncols));
                        const int nc = static_cast<int>(cand % static_cast<std::size_t>(ncols));

                        const int dr = nr - crow;
                        const int dc = nc - ccol;
                        const int step2 = dr * dr + dc * dc;

                        if (step2 < best_step2 ||
                            (step2 == best_step2 && cand < next))
                        {
                            best_step2 = step2;
                            next = cand;
                        }
                    }

                    prev = curr;
                    curr = next;
                }

                const bool eligible =
                    terminal == trace_terminal::hit_junction ||
                    terminal == trace_terminal::hit_endpoint;

                if (!eligible)
                {
                    continue;
                }

                if (static_cast<int>(removable_pixels.size()) <= max_spur_pixels)
                {
                    any_deleted_this_iter = true;
                    for (const std::size_t idx : removable_pixels)
                    {
                        delete_mark[idx] = uint8_t{1};
                    }
                }
            }
        }

        if (!any_deleted_this_iter)
        {
            break;
        }

        for (std::size_t i = 0; i < npix; ++i)
        {
            if (delete_mark[i] != 0)
            {
                work[i] = uint8_t{0};
            }
        }
    }

    // -------------------------------------------------------------------------
    // Convert the working plane back to a raster with the original valid domain.
    // -------------------------------------------------------------------------
    std::vector<float>   out_data(npix, nodata);
    std::vector<uint8_t> out_mask(npix, uint8_t{0});

    for (std::size_t i = 0; i < npix; ++i)
    {
        if (!valid[i]) continue;

        out_mask[i] = uint8_t{1};
        out_data[i] = (work[i] != 0) ? 1.0f : 0.0f;
    }

    return lsdtt_raster_cpu(
        std::move(out_data),
        std::move(out_mask),
        meta
    );
}



// -----------------------------------------------------------------------------
// prune_terminal_loops
//
// Segment-graph version.
//
// Removes small cycle-containing clusters of extracted skeleton segments when
// those clusters are attached to the rest of the skeleton through at most one
// OUTSIDE SEGMENT COMPONENT.
//
// This is intended to catch:
//   - isolated small loops
//   - spur-with-loop artefacts
//   - small loop clusters attached to the end of a long line
//
// Strategy:
//   1. extract ordered skeleton segments
//   2. build a segment adjacency graph from shared endpoint node coordinates
//   3. find connected components of "short" segments
//   4. keep only those short-segment components that:
//        - contain a cycle
//        - have <= max_loop_pixels unique foreground pixels total
//        - touch <= 1 external segment component
//   5. delete all pixels in those components except attachment-node pixels
//      that touch the outside skeleton
// -----------------------------------------------------------------------------
lsdtt_raster_cpu prune_terminal_loops(
    const lsdtt_raster_cpu& skeleton,
    int                     max_loop_pixels)
{
    const lsdtt_raster_meta& meta = skeleton.get_meta();

    if (!meta.is_valid())
    {
        throw std::runtime_error(
            "lsdtt::prune_terminal_loops: input raster metadata is invalid"
        );
    }

    if (max_loop_pixels < 1)
    {
        throw std::runtime_error(
            "lsdtt::prune_terminal_loops: max_loop_pixels must be >= 1"
        );
    }

    const int nrows = meta.n_rows;
    const int ncols = meta.n_cols;
    const std::size_t npix =
        static_cast<std::size_t>(nrows) * static_cast<std::size_t>(ncols);

    const std::vector<float>&   in_data = skeleton.get_data();
    const std::vector<uint8_t>& in_mask = skeleton.get_mask();

    const float nodata = static_cast<float>(meta.no_data_value);

    // -------------------------------------------------------------------------
    // Build valid/original foreground planes.
    // -------------------------------------------------------------------------
    std::vector<uint8_t> valid(npix, uint8_t{0});
    std::vector<uint8_t> orig (npix, uint8_t{0});

    for (std::size_t i = 0; i < npix; ++i)
    {
        if (in_mask[i] == 0) continue;
        if (meta.is_nodata(in_data[i])) continue;
        if (!std::isfinite(in_data[i])) continue;

        valid[i] = uint8_t{1};
        orig[i]  = (in_data[i] > 0.5f) ? uint8_t{1} : uint8_t{0};
    }

    auto idx_of = [ncols](int row, int col) -> std::size_t
    {
        return static_cast<std::size_t>(row) * static_cast<std::size_t>(ncols) +
               static_cast<std::size_t>(col);
    };

    // -------------------------------------------------------------------------
    // Extract ordered segments from the current skeleton.
    // -------------------------------------------------------------------------
    const std::vector<skeleton_segment> segments =
        extract_skeleton_segments(skeleton);

    if (segments.empty())
    {
        return skeleton;
    }

    using node_key = std::pair<int,int>;

    // -------------------------------------------------------------------------
    // Build node -> incident segments map.
    // -------------------------------------------------------------------------
    std::map<node_key, std::vector<std::size_t>> node_to_segments;

    for (std::size_t i = 0; i < segments.size(); ++i)
    {
        const skeleton_segment& seg = segments[i];
        if (seg.pixels.empty()) continue;

        const node_key k0(seg.start_node.row, seg.start_node.col);
        const node_key k1(seg.end_node.row,   seg.end_node.col);

        node_to_segments[k0].push_back(i);
        if (k1 != k0)
            node_to_segments[k1].push_back(i);
    }

    // -------------------------------------------------------------------------
    // Build segment adjacency from shared endpoint nodes.
    // -------------------------------------------------------------------------
    std::vector<std::vector<std::size_t>> adj(segments.size());

    for (const auto& kv : node_to_segments)
    {
        const std::vector<std::size_t>& inc = kv.second;
        for (std::size_t a = 0; a < inc.size(); ++a)
        {
            for (std::size_t b = a + 1; b < inc.size(); ++b)
            {
                const std::size_t i = inc[a];
                const std::size_t j = inc[b];
                if (i == j) continue;

                adj[i].push_back(j);
                adj[j].push_back(i);
            }
        }
    }

    // -------------------------------------------------------------------------
    // Define "short" segments as those no longer than max_loop_pixels.
    //
    // We later ALSO require the total unique pixel count of the whole candidate
    // cluster to be <= max_loop_pixels.
    // -------------------------------------------------------------------------
    std::vector<uint8_t> short_seg(segments.size(), uint8_t{0});

    for (std::size_t i = 0; i < segments.size(); ++i)
    {
        if (!segments[i].pixels.empty() &&
            static_cast<int>(segments[i].pixels.size()) <= max_loop_pixels)
        {
            short_seg[i] = uint8_t{1};
        }
    }

    // -------------------------------------------------------------------------
    // Mark pixels to delete.
    // -------------------------------------------------------------------------
    std::vector<uint8_t> delete_mark(npix, uint8_t{0});
    std::vector<uint8_t> visited_short(segments.size(), uint8_t{0});

    // -------------------------------------------------------------------------
    // Connected components of SHORT segments only.
    // -------------------------------------------------------------------------
    for (std::size_t seed = 0; seed < segments.size(); ++seed)
    {
        if (!short_seg[seed]) continue;
        if (visited_short[seed]) continue;

        std::vector<std::size_t> comp_segments;
        std::queue<std::size_t> q;

        q.push(seed);
        visited_short[seed] = uint8_t{1};

        while (!q.empty())
        {
            const std::size_t s = q.front();
            q.pop();
            comp_segments.push_back(s);

            for (const std::size_t nbr : adj[s])
            {
                if (!short_seg[nbr]) continue;
                if (visited_short[nbr]) continue;

                visited_short[nbr] = uint8_t{1};
                q.push(nbr);
            }
        }

        if (comp_segments.empty())
            continue;

        // ---------------------------------------------------------------------
        // Collect unique endpoint nodes and unique pixels in the component.
        // ---------------------------------------------------------------------
        std::set<node_key> comp_nodes;
        std::vector<std::size_t> comp_pixels;
        comp_pixels.reserve(comp_segments.size() * 8);

        bool contains_cycle = false;

        for (const std::size_t sidx : comp_segments)
        {
            const skeleton_segment& seg = segments[sidx];
            if (seg.pixels.empty()) continue;

            comp_nodes.insert(node_key(seg.start_node.row, seg.start_node.col));
            comp_nodes.insert(node_key(seg.end_node.row,   seg.end_node.col));

            const bool self_returning =
                seg.start_node.row == seg.end_node.row &&
                seg.start_node.col == seg.end_node.col;

            if (seg.is_loop || self_returning)
                contains_cycle = true;

            for (const auto& px : seg.pixels)
                comp_pixels.push_back(idx_of(px.row, px.col));
        }

        std::sort(comp_pixels.begin(), comp_pixels.end());
        comp_pixels.erase(std::unique(comp_pixels.begin(), comp_pixels.end()),
                          comp_pixels.end());

        const int E = static_cast<int>(comp_segments.size());
        const int V = static_cast<int>(comp_nodes.size());

        // Connected undirected segment-node graph has a cycle if E >= V.
        if (!contains_cycle && E >= V)
            contains_cycle = true;

        if (!contains_cycle)
            continue;

        if (static_cast<int>(comp_pixels.size()) > max_loop_pixels)
            continue;

        // ---------------------------------------------------------------------
        // Build fast membership mask for this candidate component.
        // ---------------------------------------------------------------------
        std::vector<uint8_t> in_comp(segments.size(), uint8_t{0});
        for (const std::size_t sidx : comp_segments)
            in_comp[sidx] = uint8_t{1};

        // ---------------------------------------------------------------------
        // Find outside segments adjacent to this candidate component.
        // Also record which candidate nodes touch the outside.
        // ---------------------------------------------------------------------
        std::set<std::size_t> outside_segments_set;
        std::set<node_key>    keep_nodes;   // attachment nodes to preserve

        for (const node_key& nk : comp_nodes)
        {
            const auto it = node_to_segments.find(nk);
            if (it == node_to_segments.end()) continue;

            bool touches_outside = false;

            for (const std::size_t sidx : it->second)
            {
                if (!in_comp[sidx])
                {
                    outside_segments_set.insert(sidx);
                    touches_outside = true;
                }
            }

            if (touches_outside)
                keep_nodes.insert(nk);
        }

        // ---------------------------------------------------------------------
        // Count connected components in the OUTSIDE segment graph, with the
        // candidate component removed.
        //
        // We prune when there are 0 or 1 outside segment components:
        //   0 -> isolated small loop
        //   1 -> terminal loop cluster at line end
        //   2+ -> internal / bridging / non-terminal cycle, keep it
        // ---------------------------------------------------------------------
        int n_outside_components = 0;

        if (!outside_segments_set.empty())
        {
            std::vector<uint8_t> outside_seen(segments.size(), uint8_t{0});

            for (const std::size_t s0 : outside_segments_set)
            {
                if (outside_seen[s0]) continue;

                ++n_outside_components;
                if (n_outside_components > 1)
                    break;

                std::queue<std::size_t> oq;
                oq.push(s0);
                outside_seen[s0] = uint8_t{1};

                while (!oq.empty())
                {
                    const std::size_t s = oq.front();
                    oq.pop();

                    for (const std::size_t nbr : adj[s])
                    {
                        if (in_comp[nbr]) continue;
                        if (outside_seen[nbr]) continue;

                        outside_seen[nbr] = uint8_t{1};
                        oq.push(nbr);
                    }
                }
            }
        }

        if (n_outside_components > 1)
            continue;

        // ---------------------------------------------------------------------
        // Prune this terminal loop cluster.
        //
        // Keep only attachment-node pixels where the cluster touches the
        // outside skeleton. Delete all other pixels in the candidate cluster.
        // ---------------------------------------------------------------------
        for (const std::size_t sidx : comp_segments)
        {
            const skeleton_segment& seg = segments[sidx];

            for (const auto& px : seg.pixels)
            {
                const node_key pk(px.row, px.col);
                if (keep_nodes.find(pk) != keep_nodes.end())
                    continue;

                const std::size_t idx = idx_of(px.row, px.col);
                if (valid[idx] != 0)
                    delete_mark[idx] = uint8_t{1};
            }
        }
    }

    // -------------------------------------------------------------------------
    // Convert back to raster.
    // -------------------------------------------------------------------------
    std::vector<float>   out_data(npix, nodata);
    std::vector<uint8_t> out_mask(npix, uint8_t{0});

    for (std::size_t i = 0; i < npix; ++i)
    {
        if (!valid[i]) continue;

        out_mask[i] = uint8_t{1};
        out_data[i] = (orig[i] != 0 && delete_mark[i] == 0) ? 1.0f : 0.0f;
    }

    return lsdtt_raster_cpu(
        std::move(out_data),
        std::move(out_mask),
        meta
    );
}


// -----------------------------------------------------------------------------
// fill_dense_false_pixels_3x3
//
// Local 3x3 repair on a binary skeleton/binary mask.
//
// For each valid background pixel, count the number of foreground pixels in its
// 8-neighbourhood. If that count is >= min_true_neighbours, the centre pixel is
// turned on. This is useful for filling tiny holes / ring artefacts before
// re-skeletonisation.
//
// The operation is applied iteratively n_iterations times.
// -----------------------------------------------------------------------------
lsdtt_raster_cpu fill_dense_false_pixels_3x3(
    const lsdtt_raster_cpu& binary,
    int                     min_true_neighbours,
    int                     n_iterations)
{
    const lsdtt_raster_meta& meta = binary.get_meta();

    if (!meta.is_valid())
    {
        throw std::runtime_error(
            "lsdtt::fill_dense_false_pixels_3x3: input raster metadata is invalid"
        );
    }

    if (min_true_neighbours < 1 || min_true_neighbours > 8)
    {
        throw std::runtime_error(
            "lsdtt::fill_dense_false_pixels_3x3: min_true_neighbours must be in [1,8]"
        );
    }

    if (n_iterations < 1)
    {
        throw std::runtime_error(
            "lsdtt::fill_dense_false_pixels_3x3: n_iterations must be >= 1"
        );
    }

    const int nrows = meta.n_rows;
    const int ncols = meta.n_cols;
    const std::size_t npix =
        static_cast<std::size_t>(nrows) * static_cast<std::size_t>(ncols);

    const std::vector<float>&   in_data = binary.get_data();
    const std::vector<uint8_t>& in_mask = binary.get_mask();

    const float nodata = static_cast<float>(meta.no_data_value);

    std::vector<uint8_t> valid(npix, uint8_t{0});
    std::vector<uint8_t> work (npix, uint8_t{0});

    for (std::size_t i = 0; i < npix; ++i)
    {
        if (in_mask[i] == 0) continue;
        if (meta.is_nodata(in_data[i])) continue;
        if (!std::isfinite(in_data[i])) continue;

        valid[i] = uint8_t{1};
        work[i]  = (in_data[i] > 0.5f) ? uint8_t{1} : uint8_t{0};
    }

    auto idx_of = [ncols](int row, int col) -> std::size_t
    {
        return static_cast<std::size_t>(row) * static_cast<std::size_t>(ncols) +
               static_cast<std::size_t>(col);
    };

    auto in_bounds = [nrows, ncols](int row, int col) -> bool
    {
        return row >= 0 && row < nrows && col >= 0 && col < ncols;
    };

    for (int iter = 0; iter < n_iterations; ++iter)
    {
        std::vector<uint8_t> next = work;

        for (int row = 0; row < nrows; ++row)
        {
            for (int col = 0; col < ncols; ++col)
            {
                const std::size_t idx = idx_of(row, col);
                if (!valid[idx]) continue;

                // Only fill background pixels.
                if (work[idx] != 0) continue;

                int n_true = 0;

                for (int dr = -1; dr <= 1; ++dr)
                {
                    for (int dc = -1; dc <= 1; ++dc)
                    {
                        if (dr == 0 && dc == 0) continue;

                        const int rr = row + dr;
                        const int cc = col + dc;
                        if (!in_bounds(rr, cc)) continue;

                        const std::size_t nidx = idx_of(rr, cc);
                        if (!valid[nidx]) continue;
                        if (work[nidx] != 0) ++n_true;
                    }
                }

                if (n_true >= min_true_neighbours)
                    next[idx] = uint8_t{1};
            }
        }

        work.swap(next);
    }

    std::vector<float>   out_data(npix, nodata);
    std::vector<uint8_t> out_mask(npix, uint8_t{0});

    for (std::size_t i = 0; i < npix; ++i)
    {
        if (!valid[i]) continue;

        out_mask[i] = uint8_t{1};
        out_data[i] = (work[i] != 0) ? 1.0f : 0.0f;
    }

    return lsdtt_raster_cpu(
        std::move(out_data),
        std::move(out_mask),
        meta
    );
}





// -----------------------------------------------------------------------------
// extract_skeleton_segments
// -----------------------------------------------------------------------------
std::vector<skeleton_segment>
extract_skeleton_segments(
    const lsdtt_raster_cpu& skeleton)
{
    const lsdtt_raster_meta& meta = skeleton.get_meta();

    if (!meta.is_valid())
    {
        throw std::runtime_error(
            "lsdtt::extract_skeleton_segments: input raster metadata is invalid"
        );
    }

    const int nrows = meta.n_rows;
    const int ncols = meta.n_cols;
    const std::size_t npix =
        static_cast<std::size_t>(nrows) * static_cast<std::size_t>(ncols);

    const std::vector<float>&   in_data = skeleton.get_data();
    const std::vector<uint8_t>& in_mask = skeleton.get_mask();

    // -------------------------------------------------------------------------
    // Build working foreground plane over the DEM-valid domain.
    // -------------------------------------------------------------------------
    std::vector<uint8_t> valid(npix, uint8_t{0});
    std::vector<uint8_t> work (npix, uint8_t{0});

    for (std::size_t i = 0; i < npix; ++i)
    {
        if (in_mask[i] == 0) continue;
        if (meta.is_nodata(in_data[i])) continue;
        if (!std::isfinite(in_data[i])) continue;

        valid[i] = uint8_t{1};
        work[i]  = (in_data[i] > 0.5f) ? uint8_t{1} : uint8_t{0};
    }

    auto idx_of = [ncols](int row, int col) -> std::size_t
    {
        return static_cast<std::size_t>(row) * static_cast<std::size_t>(ncols) +
               static_cast<std::size_t>(col);
    };

    auto in_bounds = [nrows, ncols](int row, int col) -> bool
    {
        return row >= 0 && row < nrows && col >= 0 && col < ncols;
    };

    auto row_of = [ncols](std::size_t idx) -> int
    {
        return static_cast<int>(idx / static_cast<std::size_t>(ncols));
    };

    auto col_of = [ncols](std::size_t idx) -> int
    {
        return static_cast<int>(idx % static_cast<std::size_t>(ncols));
    };

    auto is_foreground = [&](int row, int col) -> bool
    {
        if (!in_bounds(row, col))
            return false;

        const std::size_t idx = idx_of(row, col);
        return valid[idx] != 0 && work[idx] != 0;
    };


    constexpr int DR8[8] = {-1,-1,-1, 0, 0, 1, 1, 1};
    constexpr int DC8[8] = {-1, 0, 1,-1, 1,-1, 0, 1};

    // -------------------------------------------------------------------------
    // Compute 8-neighbour foreground degree.
    // -------------------------------------------------------------------------
    std::vector<int> degree(npix, 0);

    for (int row = 0; row < nrows; ++row)
    {
        for (int col = 0; col < ncols; ++col)
        {
            const std::size_t idx = idx_of(row, col);
            if (!valid[idx] || !work[idx]) continue;

            int deg = 0;
            for (int d = 0; d < 8; ++d)
            {
                const int rr = row + DR8[d];
                const int cc = col + DC8[d];
                if (!in_bounds(rr, cc)) continue;

                const std::size_t nidx = idx_of(rr, cc);
                if (valid[nidx] && work[nidx]) ++deg;
            }

            degree[idx] = deg;
        }
    }

    // -------------------------------------------------------------------------
    // Compute crossing number for topological node classification.
    // Use this for endpoint / chain / junction logic; raw degree is retained
    // only as a neighbour-count diagnostic.
    // -------------------------------------------------------------------------
    std::vector<int> cn(npix, 0);

    for (int row = 0; row < nrows; ++row)
    {
        for (int col = 0; col < ncols; ++col)
        {
            const std::size_t idx = idx_of(row, col);
            if (!valid[idx] || !work[idx]) continue;
            cn[idx] = crossing_number(row, col, is_foreground);
        }
    }

    auto foreground_neighbours = [&](std::size_t idx) -> std::vector<std::size_t>
    {
        std::vector<std::size_t> nbrs;
        nbrs.reserve(8);

        const int row = row_of(idx);
        const int col = col_of(idx);

        for (int d = 0; d < 8; ++d)
        {
            const int rr = row + DR8[d];
            const int cc = col + DC8[d];
            if (!in_bounds(rr, cc)) continue;

            const std::size_t nidx = idx_of(rr, cc);
            if (valid[nidx] && work[nidx])
                nbrs.push_back(nidx);
        }

        return nbrs;
    };





    using edge_key = std::pair<std::size_t, std::size_t>;

    auto make_edge_key = [](std::size_t a, std::size_t b) -> edge_key
    {
        return (a < b) ? edge_key(a, b) : edge_key(b, a);
    };

    std::set<edge_key> visited_edges;
    std::vector<skeleton_segment> segments;

    // -------------------------------------------------------------------------
    // Step 1: degree-0 singleton pixels
    // -------------------------------------------------------------------------
    for (int row = 0; row < nrows; ++row)
    {
        for (int col = 0; col < ncols; ++col)
        {
            const std::size_t idx = idx_of(row, col);
            if (!valid[idx] || !work[idx]) continue;
            if (degree[idx] != 0) continue;

            skeleton_segment seg;
            seg.is_loop    = false;
            seg.start_node = {row, col, 0, 0};
            seg.end_node   = {row, col, 0, 0};
            seg.pixels.push_back({row, col});

            segments.push_back(std::move(seg));
        }
    }

    // -------------------------------------------------------------------------
    // Step 2: open segments starting from node pixels (degree != 2, degree > 0)
    // -------------------------------------------------------------------------
    for (int row = 0; row < nrows; ++row)
    {
        for (int col = 0; col < ncols; ++col)
        {
            const std::size_t start_idx = idx_of(row, col);
            if (!valid[start_idx] || !work[start_idx]) continue;
            if (degree[start_idx] == 0) continue;
            if (cn[start_idx] == 2) continue;


            const std::vector<std::size_t> nbrs =
                foreground_neighbours(start_idx);

            for (const std::size_t first_next : nbrs)
            {
                const edge_key first_edge =
                    make_edge_key(start_idx, first_next);

                if (visited_edges.find(first_edge) != visited_edges.end())
                    continue;

                visited_edges.insert(first_edge);

                skeleton_segment seg;
                seg.is_loop    = false;
                seg.start_node = {row, col, degree[start_idx], cn[start_idx]};
                seg.pixels.push_back({row, col});

                std::size_t prev = start_idx;
                std::size_t curr = first_next;

                while (true)
                {
                    const int crow = row_of(curr);
                    const int ccol = col_of(curr);

                    seg.pixels.push_back({crow, ccol});

                    if (cn[curr] != 2)
                    {
                        seg.end_node = {crow, ccol, degree[curr], cn[curr]};
                        break;
                    }


                    const std::vector<std::size_t> curr_nbrs =
                        foreground_neighbours(curr);

                    if (curr_nbrs.size() != 2)
                    {
                        // Conservative fallback on malformed local topology.
                        seg.end_node = {crow, ccol, degree[curr], cn[curr]};
                        break;
                    }

                    std::size_t next = static_cast<std::size_t>(-1);

                    if (curr_nbrs[0] == prev && curr_nbrs[1] != prev)
                    {
                        next = curr_nbrs[1];
                    }
                    else if (curr_nbrs[1] == prev && curr_nbrs[0] != prev)
                    {
                        next = curr_nbrs[0];
                    }
                    else
                    {
                        seg.end_node = {crow, ccol, degree[curr], cn[curr]};
                        break;
                    }

                    const edge_key e = make_edge_key(curr, next);

                    if (visited_edges.find(e) != visited_edges.end())
                    {
                        seg.end_node = {crow, ccol, degree[curr], cn[curr]};
                        break;
                    }

                    visited_edges.insert(e);
                    prev = curr;
                    curr = next;
                }

                segments.push_back(std::move(seg));
            }
        }
    }

    // -------------------------------------------------------------------------
    // Step 3: closed loops (all remaining unvisited degree-2 foreground pixels)
    // -------------------------------------------------------------------------
    for (int row = 0; row < nrows; ++row)
    {
        for (int col = 0; col < ncols; ++col)
        {
            const std::size_t seed = idx_of(row, col);
            if (!valid[seed] || !work[seed]) continue;
            if (cn[seed] != 2) continue;


            const std::vector<std::size_t> seed_nbrs =
                foreground_neighbours(seed);

            std::size_t first_next = static_cast<std::size_t>(-1);

            for (const std::size_t nidx : seed_nbrs)
            {
                if (visited_edges.find(make_edge_key(seed, nidx)) ==
                    visited_edges.end())
                {
                    first_next = nidx;
                    break;
                }
            }

            if (first_next == static_cast<std::size_t>(-1))
                continue;   // this loop already traced

            skeleton_segment seg;
            seg.is_loop    = true;
            seg.start_node = {row, col, degree[seed], cn[seed]};
            seg.end_node   = {row, col, degree[seed], cn[seed]};

            seg.pixels.push_back({row, col});

            std::size_t prev = seed;
            std::size_t curr = first_next;

            visited_edges.insert(make_edge_key(seed, curr));

            while (true)
            {
                const int crow = row_of(curr);
                const int ccol = col_of(curr);

                seg.pixels.push_back({crow, ccol});

                const std::vector<std::size_t> curr_nbrs =
                    foreground_neighbours(curr);

                if (curr_nbrs.size() != 2)
                {
                    // Unexpected malformed topology; degrade gracefully.
                    seg.is_loop  = false;
                    seg.end_node = {crow, ccol, degree[curr], cn[curr]};
                    break;
                }

                std::size_t next = static_cast<std::size_t>(-1);

                if (curr_nbrs[0] != prev)
                {
                    next = curr_nbrs[0];
                }
                else if (curr_nbrs[1] != prev)
                {
                    next = curr_nbrs[1];
                }

                if (next == static_cast<std::size_t>(-1))
                {
                    seg.is_loop  = false;
                    seg.end_node = {crow, ccol, degree[curr], cn[curr]};
                    break;
                }

                if (next == seed)
                {
                    visited_edges.insert(make_edge_key(curr, seed));
                    break;   // close loop without repeating the seed pixel
                }

                const edge_key e = make_edge_key(curr, next);

                if (visited_edges.find(e) != visited_edges.end())
                {
                    seg.is_loop  = false;
                    seg.end_node = {crow, ccol, degree[curr], cn[curr]};
                    break;
                }

                visited_edges.insert(e);
                prev = curr;
                curr = next;
            }

            segments.push_back(std::move(seg));
        }
    }

    return segments;
}


// -----------------------------------------------------------------------------
// compute_skeleton_segment_strength
// -----------------------------------------------------------------------------
std::vector<skeleton_segment_strength>
compute_skeleton_segment_strength(
    const std::vector<skeleton_segment>& segments,
    const lsdtt_raster_cpu&             strength_raster)
{
    const lsdtt_raster_meta& meta = strength_raster.get_meta();

    if (!meta.is_valid())
    {
        throw std::runtime_error(
            "lsdtt::compute_skeleton_segment_strength: strength raster metadata is invalid"
        );
    }

    const int nrows = meta.n_rows;
    const int ncols = meta.n_cols;

    const std::vector<float>&   data = strength_raster.get_data();
    const std::vector<uint8_t>& mask = strength_raster.get_mask();

    auto in_bounds = [nrows, ncols](int row, int col) -> bool
    {
        return row >= 0 && row < nrows && col >= 0 && col < ncols;
    };

    auto idx_of = [ncols](int row, int col) -> std::size_t
    {
        return static_cast<std::size_t>(row) * static_cast<std::size_t>(ncols) +
               static_cast<std::size_t>(col);
    };

    std::vector<skeleton_segment_strength> out;
    out.reserve(segments.size());

    for (std::size_t s = 0; s < segments.size(); ++s)
    {
        skeleton_segment_strength ss;
        ss.segment_index = s;
        ss.n_segment_pixels = static_cast<int>(segments[s].pixels.size());

        double sum = 0.0;
        double vmax = 0.0;
        int nvalid = 0;

        for (const auto& px : segments[s].pixels)
        {
            if (!in_bounds(px.row, px.col))
            {
                throw std::runtime_error(
                    "lsdtt::compute_skeleton_segment_strength: segment pixel lies outside raster bounds"
                );
            }

            const std::size_t idx = idx_of(px.row, px.col);

            if (mask[idx] == 0) continue;
            if (meta.is_nodata(data[idx])) continue;
            if (!std::isfinite(data[idx])) continue;

            const double v = static_cast<double>(data[idx]);
            sum += v;
            if (v > vmax) vmax = v;
            ++nvalid;
        }

        ss.n_valid_strength_pixels = nvalid;
        if (nvalid > 0)
        {
            ss.mean_strength = sum / static_cast<double>(nvalid);
            ss.max_strength  = vmax;
        }

        out.push_back(ss);
    }

    return out;
}



// -----------------------------------------------------------------------------
// compute_skeleton_topology
// -----------------------------------------------------------------------------
skeleton_topology_summary
compute_skeleton_topology(
    const std::vector<skeleton_segment>& segments)
{
    skeleton_topology_summary summary;

    if (segments.empty())
        return summary;

    // -------------------------------------------------------------------------
    // Walk every segment to:
    //   (a) classify it (singleton / loop / open),
    //   (b) register its endpoint nodes in a map keyed by (row, col),
    //   (c) accumulate segment lengths and orientation vectors.
    //
    // Node registration
    // -----------------
    // We map each unique pixel position to its crossing number (CN).  The CN is
    // a property of the raster pixel and is the same no matter which segment
    // reports it, so the first encounter is sufficient.
    //
    // Node classes:
    //   CN = 0  -> singleton isolated pixel (vertex only, not a graph edge)
    //   CN = 1  -> endpoint / dangling tip
    //   CN = 2  -> chain pixel (interior of a segment, not stored as a node)
    //   CN >= 3 -> junction / branch point
    //
    // Segment classes (after extraction by extract_skeleton_segments):
    //   Singleton : is_loop=false, both endpoints have CN=0, 1 pixel
    //   Loop      : is_loop=true,  start_node == end_node, CN=2 at both
    //   Open      : all remaining (CN=1..n at endpoints, start != end or self-
    //               returning with CN>=3 at both ends)
    // -------------------------------------------------------------------------
    using node_key = std::pair<int, int>;

    // node_cn maps unique node (row, col) -> crossing number.
    std::map<node_key, int> node_cn;

    int n_singleton_segs = 0;
    int n_loop_segs      = 0;
    int n_open_segs      = 0;

    // Segment pixel counts for length statistics.
    // Singletons are excluded because they are not graph edges.
    std::vector<double> seg_pixel_counts;
    seg_pixel_counts.reserve(segments.size());

    // -------------------------------------------------------------------------
    // Orientation accumulation via doubled-angle circular statistics.
    //
    // For undirected line orientations theta in [0, pi), the standard trick is
    // to double the angles, accumulate unit vectors in the complex plane, then
    // halve the resultant angle to recover the dominant orientation.
    //
    //   sum_k w_k exp(2i theta_k)  =  (sum_cos) + i (sum_sin)
    //
    // The magnitude of the normalised resultant:
    //   R = sqrt(sum_cos^2 + sum_sin^2) / sum_w
    // is the mean resultant length, also called the anisotropy index.
    //   R = 0 -> orientations uniformly distributed (isotropic)
    //   R = 1 -> all orientations identical (perfectly aligned)
    //
    // Each segment is weighted by its pixel count so that long segments
    // contribute more than short ones.
    //
    // Segment orientation is taken from the chord between start_node and
    // end_node.  Loop segments and zero-chord self-returning open segments are
    // excluded (they carry no net directional information).
    //
    // Reference: Mardia & Jupp (2000), Directional Statistics, Wiley, sec. 2.3.
    // -------------------------------------------------------------------------
    double sum_cos = 0.0;
    double sum_sin = 0.0;
    double sum_w   = 0.0;

    for (const skeleton_segment& seg : segments)
    {
        if (seg.pixels.empty())
            continue;

        const node_key k_start(seg.start_node.row, seg.start_node.col);
        const node_key k_end  (seg.end_node.row,   seg.end_node.col);

        // Register both endpoint node positions.
        // For singletons and pure loops start == end, so only one entry is made.
        if (node_cn.find(k_start) == node_cn.end())
            node_cn[k_start] = seg.start_node.crossing_number;

        if (node_cn.find(k_end) == node_cn.end())
            node_cn[k_end] = seg.end_node.crossing_number;

        // -----------------------------------------------------------------
        // Singleton detection:
        //   A singleton is an isolated degree-0 pixel with CN=0 on both ends
        //   and is_loop=false.  It contributes one vertex but zero graph edges.
        //   Exclude it from all per-segment statistics.
        // -----------------------------------------------------------------
        if (!seg.is_loop &&
            seg.start_node.crossing_number == 0 &&
            seg.end_node.crossing_number   == 0)
        {
            ++n_singleton_segs;
            continue;
        }

        // -----------------------------------------------------------------
        // Classify as loop or open, then accumulate length and orientation.
        // -----------------------------------------------------------------
        if (seg.is_loop)
            ++n_loop_segs;
        else
            ++n_open_segs;

        const double len =
            static_cast<double>(seg.pixels.size());

        seg_pixel_counts.push_back(len);

        // Orientation: open segments only, using chord from start to end node.
        // Loop segments have start == end so the chord is zero; self-returning
        // open segments (attached small loops at junctions) likewise have zero
        // chord.  The chord_len > 0 guard handles both cases.
        if (!seg.is_loop)
        {
            const double dr = static_cast<double>(
                seg.end_node.row - seg.start_node.row);
            const double dc = static_cast<double>(
                seg.end_node.col - seg.start_node.col);

            const double chord_len = std::sqrt(dr * dr + dc * dc);

            if (chord_len > 0.0)
            {
                // COMPASS axial orientation: 0 = north, clockwise, [0, 180).
                // dc is a column difference (east-positive) and dr a ROW
                // difference (south-positive), so north is -dr and the compass
                // form is atan2(east, north) = atan2(dc, -dr).  Using
                // atan2(dc, dr) - as this did until 2026-07-22 - mirrors the
                // result: N-S and E-W came out right by accident, but a NE-SW
                // axis reported 135 where compass says 45.
                //
                // atan2 in (-pi, pi). Fold to [0, pi) for undirected orientation.
                double theta = std::atan2(dc, -dr);
                if (theta < 0.0) theta += 3.14159265358979323846;

                const double two_theta = 2.0 * theta;
                sum_cos += len * std::cos(two_theta);
                sum_sin += len * std::sin(two_theta);
                sum_w   += len;
            }
        }
    }

    // -------------------------------------------------------------------------
    // Count endpoint and junction nodes from the registered node map.
    //
    // CN = 1  -> endpoint (dangling tip, counts as a topological free end)
    // CN >= 3 -> junction (branch point)
    // CN = 2  -> pure-loop seeds; not endpoints or junctions in the graph sense
    // CN = 0  -> singleton pixels; isolated vertices, not topological nodes
    // -------------------------------------------------------------------------
    int n_endpoint_nodes = 0;
    int n_junction_nodes = 0;

    for (const auto& kv : node_cn)
    {
        const int cn = kv.second;
        if      (cn == 1) ++n_endpoint_nodes;
        else if (cn >= 3) ++n_junction_nodes;
    }

    // -------------------------------------------------------------------------
    // Graph counts for the Euler formula.
    //
    //   V = total unique node positions (all CN classes, including singletons
    //       and loop seeds with CN=2).
    //   E = true graph edges = non-singleton segments (open + loop).
    //   C = number of connected components, computed below.
    // -------------------------------------------------------------------------
    const int V = static_cast<int>(node_cn.size());
    const int E = n_open_segs + n_loop_segs;

    // -------------------------------------------------------------------------
    // Connected components.
    //
    // compute_skeleton_segment_component_ids_ performs a BFS over the segment
    // adjacency graph: two segments are adjacent when they share a node
    // (row, col) coordinate.  Each singleton becomes its own component because
    // it shares its unique position with no other segment.
    //
    // Component IDs are assigned as 0, 1, 2, ... so the total count is
    // max_id + 1 when the list is non-empty.
    // -------------------------------------------------------------------------
    const std::vector<int> comp_ids =
        compute_skeleton_segment_component_ids_(segments);

    int max_comp_id = -1;
    for (const int c : comp_ids)
    {
        if (c > max_comp_id) max_comp_id = c;
    }

    const int C = (max_comp_id >= 0) ? max_comp_id + 1 : 0;

    // -------------------------------------------------------------------------
    // First Betti number (cycle rank, beta_1).
    //
    // For a finite graph with V vertices, E edges, and C connected components:
    //   beta_1 = E - V + C
    //
    // This equals the dimension of the cycle space - the minimum number of edges
    // that must be removed to make the graph a spanning forest.  Each independent
    // loop in the skeleton contributes exactly 1 to beta_1.
    //
    // Consistency check for special cases:
    //   Singleton alone:   E=0, V=1, C=1  -> beta_1 = 0  (correct, no loops)
    //   Isolated loop:     E=1, V=1, C=1  -> beta_1 = 1  (correct, one cycle)
    //   Open path A-B:     E=1, V=2, C=1  -> beta_1 = 0  (correct, acyclic)
    //   Triangle A-B-C-A:  E=3, V=3, C=1  -> beta_1 = 1  (correct)
    //
    // Reference: Hatcher (2002), Algebraic Topology, Cambridge UP, section 2.2.
    // -------------------------------------------------------------------------
    const int beta_1 = E - V + C;

    // Graph Euler characteristic chi = beta_0 - beta_1, where beta_0 = C.
    const int chi = C - beta_1;

    // -------------------------------------------------------------------------
    // Segment length statistics.
    //
    // Lengths are in pixels (pixel count along the ordered segment chain).
    // This is a proxy for arc length; callers can multiply by cell size
    // in metres for physical units.
    // -------------------------------------------------------------------------
    double mean_len  = 0.0;
    double std_len   = 0.0;
    double total_len = 0.0;

    if (!seg_pixel_counts.empty())
    {
        for (const double l : seg_pixel_counts)
            total_len += l;

        mean_len = total_len /
            static_cast<double>(seg_pixel_counts.size());

        double var = 0.0;
        for (const double l : seg_pixel_counts)
        {
            const double d = l - mean_len;
            var += d * d;
        }
        var /= static_cast<double>(seg_pixel_counts.size());
        std_len = std::sqrt(var);
    }

    // -------------------------------------------------------------------------
    // Orientation anisotropy and dominant direction.
    //
    // Normalise the accumulated doubled-angle sum to get the mean resultant:
    //   (rx, ry) = (sum_cos / sum_w, sum_sin / sum_w)
    //   R = sqrt(rx^2 + ry^2) in [0, 1]
    //
    // Recover the dominant undirected orientation from the resultant angle phi:
    //   phi = atan2(ry, rx) in (-pi, pi)
    //   theta* = phi / 2               (halve the doubled angle)
    //   if theta* < 0: theta* += pi    (map to [0, pi))
    //
    // The dominant orientation in degrees is theta* * 180 / pi in [0, 180).
    // -------------------------------------------------------------------------
    double anisotropy               = 0.0;
    double dominant_orientation_deg = 0.0;

    if (sum_w > 0.0)
    {
        const double rx = sum_cos / sum_w;
        const double ry = sum_sin / sum_w;

        anisotropy = std::sqrt(rx * rx + ry * ry);

        // atan2 of the doubled-angle resultant, then halve to get [0, pi).
        constexpr double pi = 3.14159265358979323846;
        double phi      = std::atan2(ry, rx);
        double theta_star = phi / 2.0;
        if (theta_star < 0.0) theta_star += pi;

        dominant_orientation_deg = theta_star * (180.0 / pi);
    }

    // -------------------------------------------------------------------------
    // Derived ratios.
    //
    // branching_ratio: proportion of all topologically significant nodes that
    //   are junctions.  Low for open/dendritic patterns, high for labyrinthine.
    //
    // loop_fraction: proportion of proper segments (non-singleton) that are
    //   closed loops.  Distinct from beta_1: a single ring has loop_fraction=1,
    //   but a figure-eight has loop_fraction=0 yet beta_1=2.
    // -------------------------------------------------------------------------
    const int topo_nodes = n_endpoint_nodes + n_junction_nodes;
    const double branching_ratio =
        topo_nodes > 0
            ? static_cast<double>(n_junction_nodes) /
              static_cast<double>(topo_nodes)
            : 0.0;

    const int n_proper_segs = n_open_segs + n_loop_segs;
    const double loop_fraction =
        n_proper_segs > 0
            ? static_cast<double>(n_loop_segs) /
              static_cast<double>(n_proper_segs)
            : 0.0;

    // -------------------------------------------------------------------------
    // Pack results.
    // -------------------------------------------------------------------------
    summary.n_segments             = n_singleton_segs + n_open_segs + n_loop_segs;
    summary.n_open_segments        = n_open_segs;
    summary.n_loop_segments        = n_loop_segs;
    summary.n_singleton_segments   = n_singleton_segs;
    summary.n_endpoint_nodes       = n_endpoint_nodes;
    summary.n_junction_nodes       = n_junction_nodes;
    summary.n_connected_components = C;
    summary.cycle_rank             = beta_1;
    summary.euler_characteristic   = chi;
    summary.branching_ratio        = branching_ratio;
    summary.loop_fraction          = loop_fraction;
    summary.mean_segment_length_pixels  = mean_len;
    summary.std_segment_length_pixels   = std_len;
    summary.total_length_pixels         = total_len;
    summary.orientation_anisotropy      = anisotropy;
    summary.dominant_orientation_degrees = dominant_orientation_deg;

    return summary;
}


// -----------------------------------------------------------------------------
// estimate_skeleton_endpoint_tangents
// -----------------------------------------------------------------------------
std::vector<skeleton_endpoint_tangent>
estimate_skeleton_endpoint_tangents(
    const std::vector<skeleton_segment>& segments,
    double                              fit_length_pixels,
    int                                 min_support_pixels)
{
    if (!std::isfinite(fit_length_pixels) || fit_length_pixels <= 0.0)
    {
        throw std::runtime_error(
            "lsdtt::estimate_skeleton_endpoint_tangents: fit_length_pixels must be finite and > 0"
        );
    }

    if (min_support_pixels < 2)
    {
        throw std::runtime_error(
            "lsdtt::estimate_skeleton_endpoint_tangents: min_support_pixels must be >= 2"
        );
    }

    std::vector<skeleton_endpoint_tangent> tangents;
    tangents.reserve(segments.size() * 2);

    for (std::size_t i = 0; i < segments.size(); ++i)
    {
        const skeleton_segment& seg = segments[i];

        skeleton_endpoint_tangent t0;
        if (estimate_endpoint_tangent_one_side_(
                seg, i, true, fit_length_pixels, min_support_pixels, t0))
        {
            tangents.push_back(std::move(t0));
        }

        skeleton_endpoint_tangent t1;
        if (estimate_endpoint_tangent_one_side_(
                seg, i, false, fit_length_pixels, min_support_pixels, t1))
        {
            tangents.push_back(std::move(t1));
        }
    }

    return tangents;
}


// -----------------------------------------------------------------------------
// find_skeleton_endpoint_pair_candidates
// -----------------------------------------------------------------------------
std::vector<skeleton_endpoint_pair_candidate>
find_skeleton_endpoint_pair_candidates(
    const std::vector<skeleton_segment>&                   segments,
    const std::vector<skeleton_endpoint_tangent>&          tangents,
    const endpoint_pair_request&                           request,
    const std::vector<skeleton_segment_strength>&          segment_strengths)
{
    validate_endpoint_pair_request(request);

    if (!segment_strengths.empty() &&
        segment_strengths.size() != segments.size())
    {
        throw std::runtime_error(
            "lsdtt::find_skeleton_endpoint_pair_candidates: segment_strengths size must match segments size"
        );
    }

    if (request.require_strong_segments && segment_strengths.empty())
    {
        throw std::runtime_error(
            "lsdtt::find_skeleton_endpoint_pair_candidates: strong-segment filtering requested but no segment strengths were supplied"
        );
    }

    auto get_segment_mean_strength = [&](std::size_t segment_index) -> double
    {
        if (segment_strengths.empty()) return 0.0;
        return segment_strengths[segment_index].mean_strength;
    };


    const std::vector<int> segment_component_id =
        compute_skeleton_segment_component_ids_(segments);

    for (const auto& t : tangents)
    {
        if (t.segment_index >= segments.size())
        {
            throw std::runtime_error(
                "lsdtt::find_skeleton_endpoint_pair_candidates: tangent segment_index out of range"
            );
        }
    }

    std::vector<skeleton_endpoint_pair_candidate> candidates;
    candidates.reserve(tangents.size());

    const double weight_sum =
        request.weight_distance +
        request.weight_facing +
        request.weight_heading;

    for (std::size_t i = 0; i < tangents.size(); ++i)
    {
        const auto& a = tangents[i];

        for (std::size_t j = i + 1; j < tangents.size(); ++j)
        {
            const auto& b = tangents[j];

            // -----------------------------------------------------------------
            // Structural exclusions
            // -----------------------------------------------------------------
            if (a.segment_index == b.segment_index)
                continue;

            const int comp_a = segment_component_id[a.segment_index];
            const int comp_b = segment_component_id[b.segment_index];

            if (request.reject_same_component && comp_a == comp_b)
                continue;

            const double seg_strength_a =
                get_segment_mean_strength(a.segment_index);
            const double seg_strength_b =
                get_segment_mean_strength(b.segment_index);

            if (request.require_strong_segments)
            {
                if (seg_strength_a < request.min_segment_mean_strength ||
                    seg_strength_b < request.min_segment_mean_strength)
                {
                    continue;
                }
            }

            // -----------------------------------------------------------------
            // Gap geometry
            // -----------------------------------------------------------------
            const double g_row = static_cast<double>(b.row - a.row);
            const double g_col = static_cast<double>(b.col - a.col);


            const double gap_distance =
                std::sqrt(g_row * g_row + g_col * g_col);

            if (!(std::isfinite(gap_distance)) || gap_distance <= 0.0)
                continue;

            if (gap_distance > request.max_gap_pixels)
                continue;

            // -----------------------------------------------------------------
            // Facing tests
            // A should point toward B; B should point toward A.
            // -----------------------------------------------------------------
            const double facing_a =
                angle_between_vectors_degrees(
                    a.d_row, a.d_col,
                    g_row,   g_col
                );

            const double facing_b =
                angle_between_vectors_degrees(
                    b.d_row, b.d_col,
                   -g_row,  -g_col
                );

            if (!(std::isfinite(facing_a) && std::isfinite(facing_b)))
                continue;

            if (facing_a > request.max_facing_angle_degrees ||
                facing_b > request.max_facing_angle_degrees)
            {
                continue;
            }

            // -----------------------------------------------------------------
            // Heading compatibility:
            // outward tangent at A should align with the opposite of outward
            // tangent at B.
            // -----------------------------------------------------------------
            const double heading_misalignment =
                angle_between_vectors_degrees(
                    a.d_row,  a.d_col,
                   -b.d_row, -b.d_col
                );

            if (!std::isfinite(heading_misalignment))
                continue;

            if (heading_misalignment > request.max_heading_misalignment_degrees)
                continue;

            // -----------------------------------------------------------------
            // Normalised penalties and score (lower is better)
            // -----------------------------------------------------------------
            const double distance_penalty =
                gap_distance / request.max_gap_pixels;

            const double facing_penalty =
                0.5 * (
                    facing_a / request.max_facing_angle_degrees +
                    facing_b / request.max_facing_angle_degrees
                );

            double heading_penalty = 0.0;
            if (request.max_heading_misalignment_degrees > 0.0)
            {
                heading_penalty =
                    heading_misalignment /
                    request.max_heading_misalignment_degrees;
            }
            else
            {
                heading_penalty = (heading_misalignment == 0.0) ? 0.0
                                                                : std::numeric_limits<double>::infinity();
            }

            if (!std::isfinite(heading_penalty))
                continue;

            const double score =
                (request.weight_distance * distance_penalty +
                 request.weight_facing   * facing_penalty   +
                 request.weight_heading  * heading_penalty) / weight_sum;

            skeleton_endpoint_pair_candidate cand;
            cand.tangent_index_a = i;
            cand.tangent_index_b = j;
            cand.segment_index_a = a.segment_index;
            cand.segment_index_b = b.segment_index;
            cand.component_id_a  = comp_a;
            cand.component_id_b  = comp_b;
            cand.row_a           = a.row;
            cand.col_a           = a.col;
            cand.row_b           = b.row;
            cand.col_b           = b.col;
            cand.gap_distance_pixels          = gap_distance;
            cand.facing_angle_a_degrees       = facing_a;
            cand.facing_angle_b_degrees       = facing_b;
            cand.heading_misalignment_degrees = heading_misalignment;
            cand.segment_mean_strength_a      = seg_strength_a;
            cand.segment_mean_strength_b      = seg_strength_b;
            cand.pair_min_segment_strength    =
                std::min(seg_strength_a, seg_strength_b);
            cand.score                        = score;


            candidates.push_back(std::move(cand));
        }
    }

    std::sort(
        candidates.begin(),
        candidates.end(),
        [](const skeleton_endpoint_pair_candidate& lhs,
           const skeleton_endpoint_pair_candidate& rhs)
        {
            if (lhs.score != rhs.score)
                return lhs.score < rhs.score;

            if (lhs.gap_distance_pixels != rhs.gap_distance_pixels)
                return lhs.gap_distance_pixels < rhs.gap_distance_pixels;

            if (lhs.tangent_index_a != rhs.tangent_index_a)
                return lhs.tangent_index_a < rhs.tangent_index_a;

            return lhs.tangent_index_b < rhs.tangent_index_b;
        }
    );

    return candidates;
}



// -----------------------------------------------------------------------------
// select_skeleton_endpoint_pairs
// -----------------------------------------------------------------------------
std::vector<skeleton_endpoint_pair_candidate>
select_skeleton_endpoint_pairs(
    const std::vector<skeleton_endpoint_pair_candidate>& candidates,
    std::size_t                                          n_tangents,
    const endpoint_pair_selection_request&               request)
{
    validate_endpoint_pair_selection_request(request);

    if (n_tangents == 0 || candidates.empty())
    {
        return {};
    }

    // -------------------------------------------------------------------------
    // Copy and sort candidates into deterministic best-first order.
    // -------------------------------------------------------------------------
    std::vector<skeleton_endpoint_pair_candidate> sorted = candidates;

    std::sort(
        sorted.begin(),
        sorted.end(),
        endpoint_pair_candidate_less_);

    // -------------------------------------------------------------------------
    // Validate indices and mark score-surviving candidates.
    // -------------------------------------------------------------------------
    std::vector<uint8_t> survives_score_filter(
        sorted.size(), uint8_t{0});

    for (std::size_t i = 0; i < sorted.size(); ++i)
    {
        const auto& c = sorted[i];

        if (c.tangent_index_a >= n_tangents ||
            c.tangent_index_b >= n_tangents)
        {
            throw std::runtime_error(
                "lsdtt::select_skeleton_endpoint_pairs: candidate tangent index out of range"
            );
        }

        if (c.tangent_index_a == c.tangent_index_b)
        {
            continue;
        }

        if (!std::isfinite(c.score))
        {
            continue;
        }

        if (c.score <= request.max_score)
        {
            survives_score_filter[i] = uint8_t{1};
        }
    }

    // -------------------------------------------------------------------------
    // Best surviving candidate incident on each tangent.
    // Because sorted[] is already in best-first order, the first surviving
    // incident candidate is the best one.
    // -------------------------------------------------------------------------
    std::vector<int> best_for_tangent(
        n_tangents, -1);

    for (std::size_t i = 0; i < sorted.size(); ++i)
    {
        if (survives_score_filter[i] == 0) continue;

        const auto& c = sorted[i];

        if (best_for_tangent[c.tangent_index_a] < 0)
            best_for_tangent[c.tangent_index_a] = static_cast<int>(i);

        if (best_for_tangent[c.tangent_index_b] < 0)
            best_for_tangent[c.tangent_index_b] = static_cast<int>(i);
    }

    // -------------------------------------------------------------------------
    // Mutual-best filter.
    // -------------------------------------------------------------------------
    std::vector<std::size_t> eligible_indices;
    eligible_indices.reserve(sorted.size());

    for (std::size_t i = 0; i < sorted.size(); ++i)
    {
        if (survives_score_filter[i] == 0) continue;

        const auto& c = sorted[i];

        bool keep = true;

        if (request.require_mutual_best)
        {
            keep =
                best_for_tangent[c.tangent_index_a] == static_cast<int>(i) &&
                best_for_tangent[c.tangent_index_b] == static_cast<int>(i);
        }

        if (keep)
        {
            eligible_indices.push_back(i);
        }
    }

    // -------------------------------------------------------------------------
    // One-to-one greedy acceptance.
    // Since eligible_indices follows sorted[] order, greedily taking unused
    // endpoints yields a consistent best-first matching.
    // -------------------------------------------------------------------------
    std::vector<skeleton_endpoint_pair_candidate> accepted;
    accepted.reserve(eligible_indices.size());

    if (request.enforce_one_to_one)
    {
        std::vector<uint8_t> tangent_used(
            n_tangents, uint8_t{0});

        for (const std::size_t i : eligible_indices)
        {
            const auto& c = sorted[i];

            if (tangent_used[c.tangent_index_a] != 0 ||
                tangent_used[c.tangent_index_b] != 0)
            {
                continue;
            }

            tangent_used[c.tangent_index_a] = uint8_t{1};
            tangent_used[c.tangent_index_b] = uint8_t{1};
            accepted.push_back(c);
        }
    }
    else
    {
        for (const std::size_t i : eligible_indices)
        {
            accepted.push_back(sorted[i]);
        }
    }

    return accepted;
}






// -----------------------------------------------------------------------------
// rasterise_skeleton_endpoint_bridges
// -----------------------------------------------------------------------------
skeleton_bridge_rasterise_result
rasterise_skeleton_endpoint_bridges(
    const lsdtt_raster_cpu&                              skeleton,
    const std::vector<skeleton_endpoint_pair_candidate>& accepted_pairs,
    const skeleton_bridge_rasterise_request&             request)
{
    validate_skeleton_bridge_rasterise_request(request);

    const lsdtt_raster_meta& meta = skeleton.get_meta();

    if (!meta.is_valid())
    {
        throw std::runtime_error(
            "lsdtt::rasterise_skeleton_endpoint_bridges: input raster metadata is invalid"
        );
    }

    const int nrows = meta.n_rows;
    const int ncols = meta.n_cols;
    const std::size_t npix =
        static_cast<std::size_t>(nrows) * static_cast<std::size_t>(ncols);

    const std::vector<float>&   in_data = skeleton.get_data();
    const std::vector<uint8_t>& in_mask = skeleton.get_mask();

    const float nodata = static_cast<float>(meta.no_data_value);

    // -------------------------------------------------------------------------
    // Build working valid/original-foreground planes over the DEM-valid domain.
    // -------------------------------------------------------------------------
    std::vector<uint8_t> valid(npix, uint8_t{0});
    std::vector<uint8_t> orig (npix, uint8_t{0});

    for (std::size_t i = 0; i < npix; ++i)
    {
        if (in_mask[i] == 0) continue;
        if (meta.is_nodata(in_data[i])) continue;
        if (!std::isfinite(in_data[i])) continue;

        valid[i] = uint8_t{1};
        orig[i]  = (in_data[i] > 0.5f) ? uint8_t{1} : uint8_t{0};
    }

    auto idx_of = [ncols](int row, int col) -> std::size_t
    {
        return static_cast<std::size_t>(row) * static_cast<std::size_t>(ncols) +
               static_cast<std::size_t>(col);
    };

    auto in_bounds = [nrows, ncols](int row, int col) -> bool
    {
        return row >= 0 && row < nrows && col >= 0 && col < ncols;
    };

    // Newly burned bridge pixels only.
    std::vector<uint8_t> bridge(npix, uint8_t{0});

    int n_requested                     = 0;
    int n_burned                        = 0;
    int n_rejected_invalid_domain_cross = 0;

    // -------------------------------------------------------------------------
    // Burn accepted straight-line bridges.
    // -------------------------------------------------------------------------
    for (const auto& cand : accepted_pairs)
    {
        ++n_requested;

        const std::vector<std::pair<int,int>> line_pixels =
            rasterise_line_bresenham_(
                cand.row_a, cand.col_a,
                cand.row_b, cand.col_b
            );

        bool line_is_fully_valid = true;

        for (const auto& rc : line_pixels)
        {
            const int row = rc.first;
            const int col = rc.second;

            if (!in_bounds(row, col))
            {
                line_is_fully_valid = false;
                break;
            }

            const std::size_t idx = idx_of(row, col);
            if (valid[idx] == 0)
            {
                line_is_fully_valid = false;
                break;
            }
        }

        if (request.require_all_bridge_pixels_valid && !line_is_fully_valid)
        {
            ++n_rejected_invalid_domain_cross;
            continue;
        }

        bool burned_this_bridge = false;

        for (const auto& rc : line_pixels)
        {
            const int row = rc.first;
            const int col = rc.second;

            if (!in_bounds(row, col))
                continue;

            const std::size_t idx = idx_of(row, col);
            if (valid[idx] == 0)
                continue;

            if (orig[idx] == 0 && bridge[idx] == 0)
            {
                bridge[idx] = uint8_t{1};
                burned_this_bridge = true;
            }
        }

        if (burned_this_bridge)
        {
            ++n_burned;
        }
    }

    // -------------------------------------------------------------------------
    // Convert bridge/original planes back to rasters.
    // -------------------------------------------------------------------------
    std::vector<float>   bridge_data(npix, nodata);
    std::vector<uint8_t> bridge_mask(npix, uint8_t{0});

    std::vector<float>   connected_data(npix, nodata);
    std::vector<uint8_t> connected_mask(npix, uint8_t{0});

    for (std::size_t i = 0; i < npix; ++i)
    {
        if (valid[i] == 0)
            continue;

        bridge_mask[i]    = uint8_t{1};
        connected_mask[i] = uint8_t{1};

        bridge_data[i] =
            (bridge[i] != 0) ? 1.0f : 0.0f;

        connected_data[i] =
            (orig[i] != 0 || bridge[i] != 0) ? 1.0f : 0.0f;
    }

    skeleton_bridge_rasterise_result result;
    result.bridge_mask = lsdtt_raster_cpu(
        std::move(bridge_data), std::move(bridge_mask), meta);
    result.connected_skeleton = lsdtt_raster_cpu(
        std::move(connected_data), std::move(connected_mask), meta);

    result.n_requested                     = n_requested;
    result.n_burned                        = n_burned;
    result.n_rejected_invalid_domain_cross = n_rejected_invalid_domain_cross;

    return result;
}



// -----------------------------------------------------------------------------
// route_skeleton_endpoint_bridges_least_cost
// -----------------------------------------------------------------------------
skeleton_bridge_route_result
route_skeleton_endpoint_bridges_least_cost(
    const lsdtt_raster_cpu&                              skeleton,
    const lsdtt_raster_cpu&                              support_raster,
    const std::vector<skeleton_endpoint_pair_candidate>& accepted_pairs,
    const skeleton_bridge_route_request&                 request)
{
    validate_skeleton_bridge_route_request(request);

    const lsdtt_raster_meta& meta = skeleton.get_meta();
    const lsdtt_raster_meta& smeta = support_raster.get_meta();

    if (!meta.is_valid())
    {
        throw std::runtime_error(
            "lsdtt::route_skeleton_endpoint_bridges_least_cost: skeleton metadata is invalid"
        );
    }

    if (!smeta.is_valid())
    {
        throw std::runtime_error(
            "lsdtt::route_skeleton_endpoint_bridges_least_cost: support raster metadata is invalid"
        );
    }

    if (!same_grid(meta, smeta))
    {
        throw std::runtime_error(
            "lsdtt::route_skeleton_endpoint_bridges_least_cost: skeleton and support raster are not co-registered"
        );
    }

    const int nrows = meta.n_rows;
    const int ncols = meta.n_cols;
    const std::size_t npix =
        static_cast<std::size_t>(nrows) * static_cast<std::size_t>(ncols);

    const std::vector<float>&   in_data = skeleton.get_data();
    const std::vector<uint8_t>& in_mask = skeleton.get_mask();

    const std::vector<float>&   sup_data = support_raster.get_data();
    const std::vector<uint8_t>& sup_mask = support_raster.get_mask();

    const float nodata = static_cast<float>(meta.no_data_value);

    // -------------------------------------------------------------------------
    // Build valid/original-skeleton planes. The skeleton valid domain is the
    // authoritative routing domain. Support nodata inside that domain simply
    // receives zero support.
    // -------------------------------------------------------------------------
    std::vector<uint8_t> valid(npix, uint8_t{0});
    std::vector<uint8_t> orig (npix, uint8_t{0});

    for (std::size_t i = 0; i < npix; ++i)
    {
        if (in_mask[i] == 0) continue;
        if (meta.is_nodata(in_data[i])) continue;
        if (!std::isfinite(in_data[i])) continue;

        valid[i] = uint8_t{1};
        orig[i]  = (in_data[i] > 0.5f) ? uint8_t{1} : uint8_t{0};
    }

    // -------------------------------------------------------------------------
    // Normalise support to [0,1] over valid pixels where support is finite.
    // Invalid support on valid skeleton pixels becomes 0 support.
    // -------------------------------------------------------------------------
    std::vector<double> support_norm(npix, 0.0);

    bool   have_support = false;
    double smin =  std::numeric_limits<double>::infinity();
    double smax = -std::numeric_limits<double>::infinity();

    for (std::size_t i = 0; i < npix; ++i)
    {
        if (valid[i] == 0) continue;
        if (sup_mask[i] == 0) continue;
        if (smeta.is_nodata(sup_data[i])) continue;
        if (!std::isfinite(sup_data[i])) continue;

        const double v = static_cast<double>(sup_data[i]);
        if (v < smin) smin = v;
        if (v > smax) smax = v;
        have_support = true;
    }

    if (have_support)
    {
        if (smax > smin)
        {
            const double inv = 1.0 / (smax - smin);
            for (std::size_t i = 0; i < npix; ++i)
            {
                if (valid[i] == 0) continue;
                if (sup_mask[i] == 0) continue;
                if (smeta.is_nodata(sup_data[i])) continue;
                if (!std::isfinite(sup_data[i])) continue;

                support_norm[i] =
                    (static_cast<double>(sup_data[i]) - smin) * inv;
            }
        }
        else
        {
            const double fill_value = (smax > 0.0) ? 1.0 : 0.0;

            for (std::size_t i = 0; i < npix; ++i)
            {
                if (valid[i] == 0) continue;
                if (sup_mask[i] == 0) continue;
                if (smeta.is_nodata(sup_data[i])) continue;
                if (!std::isfinite(sup_data[i])) continue;

                support_norm[i] = fill_value;
            }
        }
    }

    auto idx_of = [ncols](int row, int col) -> std::size_t
    {
        return static_cast<std::size_t>(row) * static_cast<std::size_t>(ncols) +
               static_cast<std::size_t>(col);
    };

    auto in_bounds = [nrows, ncols](int row, int col) -> bool
    {
        return row >= 0 && row < nrows && col >= 0 && col < ncols;
    };

    auto row_of = [ncols](std::size_t idx) -> int
    {
        return static_cast<int>(idx / static_cast<std::size_t>(ncols));
    };

    auto col_of = [ncols](std::size_t idx) -> int
    {
        return static_cast<int>(idx % static_cast<std::size_t>(ncols));
    };

    constexpr int DR8[8] = {-1,-1,-1, 0, 0, 1, 1, 1};
    constexpr int DC8[8] = {-1, 0, 1,-1, 1,-1, 0, 1};

    struct open_node
    {
        double      f   = std::numeric_limits<double>::infinity();
        double      g   = std::numeric_limits<double>::infinity();
        std::size_t idx = 0;
    };

    struct open_node_greater
    {
        bool operator()(const open_node& a, const open_node& b) const
        {
            if (a.f != b.f) return a.f > b.f;
            if (a.g != b.g) return a.g > b.g;
            return a.idx > b.idx;
        }
    };

    std::vector<uint8_t> bridge(npix, uint8_t{0});

    int n_requested               = 0;
    int n_routed                  = 0;
    int n_failed_invalid_endpoint = 0;
    int n_failed_no_path          = 0;

    const double min_cost_scale =
        std::min(1.0, request.foreground_cost_scale);

    const double heuristic_cost_per_pixel =
        request.base_step_cost * min_cost_scale;

    // -------------------------------------------------------------------------
    // Route each accepted pair independently.
    // -------------------------------------------------------------------------
    for (const auto& cand : accepted_pairs)
    {
        ++n_requested;

        if (!in_bounds(cand.row_a, cand.col_a) ||
            !in_bounds(cand.row_b, cand.col_b))
        {
            ++n_failed_invalid_endpoint;
            continue;
        }

        const std::size_t start = idx_of(cand.row_a, cand.col_a);
        const std::size_t goal  = idx_of(cand.row_b, cand.col_b);

        if (valid[start] == 0 || valid[goal] == 0)
        {
            ++n_failed_invalid_endpoint;
            continue;
        }

        // ---------------------------------------------------------------------
        // Build a local corridor around the straight endpoint-to-endpoint line.
        // ---------------------------------------------------------------------
        const int pad =
            std::max(1, static_cast<int>(
                std::ceil(request.corridor_half_width_pixels))) + 1;

        const int rmin = std::max(0, std::min(cand.row_a, cand.row_b) - pad);
        const int rmax = std::min(nrows - 1, std::max(cand.row_a, cand.row_b) + pad);
        const int cmin = std::max(0, std::min(cand.col_a, cand.col_b) - pad);
        const int cmax = std::min(ncols - 1, std::max(cand.col_a, cand.col_b) + pad);

        std::vector<uint8_t> in_corridor(npix, uint8_t{0});

        for (int row = rmin; row <= rmax; ++row)
        {
            for (int col = cmin; col <= cmax; ++col)
            {
                const std::size_t idx = idx_of(row, col);
                if (valid[idx] == 0) continue;

                const bool is_endpoint =
                    (row == cand.row_a && col == cand.col_a) ||
                    (row == cand.row_b && col == cand.col_b);

                if (is_endpoint)
                {
                    in_corridor[idx] = uint8_t{1};
                    continue;
                }

                const double dist =
                    point_to_segment_distance(
                        static_cast<double>(row),
                        static_cast<double>(col),
                        static_cast<double>(cand.row_a),
                        static_cast<double>(cand.col_a),
                        static_cast<double>(cand.row_b),
                        static_cast<double>(cand.col_b)
                    );

                if (dist <= request.corridor_half_width_pixels)
                {
                    in_corridor[idx] = uint8_t{1};
                }
            }
        }

        if (in_corridor[start] == 0 || in_corridor[goal] == 0)
        {
            ++n_failed_no_path;
            continue;
        }

        // ---------------------------------------------------------------------
        // A* over 8-neighbour pixels inside the corridor.
        // ---------------------------------------------------------------------
        std::vector<double> gscore(npix,
            std::numeric_limits<double>::infinity());
        std::vector<int>    parent(npix, -1);
        std::vector<uint8_t> closed(npix, uint8_t{0});

        auto heuristic = [&](std::size_t idx) -> double
        {
            const double dr =
                static_cast<double>(row_of(goal) - row_of(idx));
            const double dc =
                static_cast<double>(col_of(goal) - col_of(idx));

            return std::sqrt(dr * dr + dc * dc) * heuristic_cost_per_pixel;
        };

        auto pixel_traversal_cost = [&](std::size_t idx) -> double
        {
            double c =
                request.base_step_cost +
                request.support_penalty_weight * (1.0 - support_norm[idx]);

            if (orig[idx] != 0)
            {
                c *= request.foreground_cost_scale;
            }

            return c;
        };

        std::priority_queue<
            open_node,
            std::vector<open_node>,
            open_node_greater> open;

        gscore[start] = 0.0;
        open.push({heuristic(start), 0.0, start});

        bool found = false;

        while (!open.empty())
        {
            const open_node cur = open.top();
            open.pop();

            if (closed[cur.idx] != 0)
                continue;

            closed[cur.idx] = uint8_t{1};

            if (cur.idx == goal)
            {
                found = true;
                break;
            }

            const int row = row_of(cur.idx);
            const int col = col_of(cur.idx);

            for (int d = 0; d < 8; ++d)
            {
                const int rr = row + DR8[d];
                const int cc = col + DC8[d];

                if (!in_bounds(rr, cc)) continue;

                const std::size_t nidx = idx_of(rr, cc);

                if (valid[nidx] == 0)       continue;
                if (in_corridor[nidx] == 0) continue;
                if (closed[nidx] != 0)      continue;

                const double step_len =
                    (DR8[d] == 0 || DC8[d] == 0)
                        ? 1.0
                        : std::sqrt(2.0);

                const double tentative_g =
                    gscore[cur.idx] +
                    step_len * pixel_traversal_cost(nidx);

                if (tentative_g < gscore[nidx])
                {
                    gscore[nidx] = tentative_g;
                    parent[nidx] = static_cast<int>(cur.idx);

                    open.push({
                        tentative_g + heuristic(nidx),
                        tentative_g,
                        nidx
                    });
                }
            }
        }

        if (!found || !std::isfinite(gscore[goal]))
        {
            ++n_failed_no_path;
            continue;
        }

        // ---------------------------------------------------------------------
        // Reconstruct path and burn only newly inserted pixels.
        // ---------------------------------------------------------------------
        std::size_t cur = goal;
        while (true)
        {
            if (orig[cur] == 0)
            {
                bridge[cur] = uint8_t{1};
            }

            if (cur == start)
                break;

            const int p = parent[cur];
            if (p < 0)
                break;

            cur = static_cast<std::size_t>(p);
        }

        ++n_routed;
    }

    // -------------------------------------------------------------------------
    // Convert back to rasters.
    // -------------------------------------------------------------------------
    std::vector<float>   bridge_data(npix, nodata);
    std::vector<uint8_t> bridge_mask(npix, uint8_t{0});

    std::vector<float>   connected_data(npix, nodata);
    std::vector<uint8_t> connected_mask(npix, uint8_t{0});

    for (std::size_t i = 0; i < npix; ++i)
    {
        if (valid[i] == 0)
            continue;

        bridge_mask[i]    = uint8_t{1};
        connected_mask[i] = uint8_t{1};

        bridge_data[i] =
            (bridge[i] != 0) ? 1.0f : 0.0f;

        connected_data[i] =
            (orig[i] != 0 || bridge[i] != 0) ? 1.0f : 0.0f;
    }

    skeleton_bridge_route_result result;
    result.bridge_mask = lsdtt_raster_cpu(
        std::move(bridge_data), std::move(bridge_mask), meta);
    result.connected_skeleton = lsdtt_raster_cpu(
        std::move(connected_data), std::move(connected_mask), meta);

    result.n_requested               = n_requested;
    result.n_routed                  = n_routed;
    result.n_failed_invalid_endpoint = n_failed_invalid_endpoint;
    result.n_failed_no_path          = n_failed_no_path;

    return result;
}


// -----------------------------------------------------------------------------
// skeleton_segment_to_polyline
// -----------------------------------------------------------------------------
std::vector<lsdtt_xy> skeleton_segment_to_polyline(
    const skeleton_segment&  segment,
    const lsdtt_raster_meta& meta)
{
    std::vector<lsdtt_xy> out;
    out.reserve(segment.pixels.size());
    for (const skeleton_pixel& px : segment.pixels)
    {
        out.push_back(cell_center_from_row_col(meta, px.row, px.col));
    }
    return out;
}


} // end namespace lsdtt