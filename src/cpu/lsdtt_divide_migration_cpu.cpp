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

} // namespace

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

} // namespace lsdtt
