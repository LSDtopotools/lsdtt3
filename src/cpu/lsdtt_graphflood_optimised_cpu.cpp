// lsdtt_graphflood_optimised_cpu.cpp

#include "lsdtt_graphflood_optimised_cpu.hpp"
#include "../utils/lsdtt_grid_neighbours.hpp"
#include "../utils/lsdtt_statistics.hpp"

#include "lsdtt_raster_resampling_cpu.hpp"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstring>
#include <limits>
#include <optional>
#include <stdexcept>
#include <utility>

namespace lsdtt
{
namespace
{

using clock_type = std::chrono::steady_clock;

double seconds_since(clock_type::time_point start)
{
    return std::chrono::duration<double>(clock_type::now() - start).count();
}

static constexpr const auto& k_dr = grid::d8_cardinal_first_row;
static constexpr const auto& k_dc = grid::d8_cardinal_first_col;
static constexpr const auto& k_pf_dr = grid::d8_clockwise_row;
static constexpr const auto& k_pf_dc = grid::d8_clockwise_col;
constexpr float k_sqrt2 = grid::k_sqrt2_f;
constexpr std::uint32_t k_no_receiver =
    std::numeric_limits<std::uint32_t>::max();

std::uint32_t ordered_float_key(float value) noexcept
{
    // Preserve float ordering with an unsigned radix key.
    std::uint32_t bits;
    std::memcpy(&bits, &value, sizeof(bits));
    return bits ^ ((bits & 0x80000000u) ? 0xffffffffu : 0x80000000u);
}

struct flood_node
{
    float z;
    std::uint32_t flat;
};

struct higher_flood_node
{
    bool operator()(const flood_node& lhs, const flood_node& rhs) const noexcept
    {
        return lhs.z > rhs.z;
    }
};

class optimised_solver
{
public:
    optimised_solver(
        const lsdtt_raster_cpu& bed,
        const lsdtt_array2d_u8& node_state,
        const std::vector<double>& rainfall_rate,
        const std::vector<double>& friction_coefficient,
        const std::vector<double>& point_inflow,
        const std::vector<double>& initial_water_depth,
        const std::vector<double>& prescribed_qin,
        const graphflood_optimised_options& options,
        bool use_v2)
        : meta_(bed.get_meta()), n_rows_(meta_.n_rows), n_cols_(meta_.n_cols),
          n_cells_(static_cast<std::size_t>(n_rows_) * n_cols_), options_(options),
          use_v2_(use_v2)
    {
        validate_inputs(bed, node_state, rainfall_rate,
                        friction_coefficient, point_inflow,
                        initial_water_depth, prescribed_qin);

        const auto& bed_data = bed.get_data();
        state_ = node_state.get_vector();
        bed_.assign(bed_data.begin(), bed_data.end());
        surface_ = bed_;
        depth_ = initial_water_depth;
        if (use_v2_ && !prescribed_qin.empty()) prescribed_qin_ = prescribed_qin;
        qin_.assign(n_cells_, 0.0);
        qout_.assign(n_cells_, 0.0);
        dhw_.assign(n_cells_, 0.0);
        slope_.assign(n_cells_, 0.0);
        width_.assign(n_cells_, static_cast<float>(meta_.cell_size));
        receiver_distance_.assign(n_cells_, static_cast<float>(meta_.cell_size));
        best_receiver_.assign(n_cells_, k_no_receiver);
        rank_.assign(n_cells_, k_no_receiver);
        flood_state_.assign(n_cells_, -1);
        interior_.assign(n_cells_, 0);
        source_q_.assign(n_cells_, 0.0);
        specific_discharge_.assign(n_cells_, 0.0);
        active_.reserve(n_cells_);
        order_.reserve(n_cells_);
        sort_scratch_.resize(active_.size());
        flood_queue_.reserve(n_cells_);
        heap_.reserve(std::min<std::size_t>(n_cells_, 65536));
        convergence_abs_.reserve(n_cells_);
        convergence_rel_.reserve(n_cells_);

        const double area = meta_.cell_size * meta_.cell_size;
        for (std::size_t flat = 0; flat < n_cells_; ++flat)
        {
            if (state_[flat] == 0) continue;
            active_.push_back(static_cast<std::uint32_t>(flat));
            source_q_[flat] = rainfall_rate[flat] * area + point_inflow[flat];
            friction_[flat] = friction_coefficient[flat];
        }

        initialise_v2_geometry();

        if (options_.check_valid_outlet) validate_outlet_components();
        rebuild_topology();
        update_working_bytes();
    }

    graphflood_optimised_result run(
        std::int64_t max_steps,
        bool run_until_converged,
        int convergence_check_every,
        double tol_abs,
        double tol_rel)
    {
        if (max_steps < 0) throw std::invalid_argument("max_steps must be >= 0");
        if (convergence_check_every <= 0)
            throw std::invalid_argument("convergence_check_every must be > 0");

        bool converged = false;
        for (std::int64_t step = 0; step < max_steps; ++step)
        {
            run_one_step();
            ++steps_;

            if (run_until_converged && steps_ % convergence_check_every == 0)
            {
                const auto start = clock_type::now();
                const auto [median_abs, median_rel] = convergence_medians();
                performance_.convergence_seconds += seconds_since(start);
                if (median_abs <= tol_abs && median_rel <= tol_rel)
                {
                    converged = true;
                    break;
                }
            }
        }

        return make_result(converged);
    }

private:
    void validate_inputs(
        const lsdtt_raster_cpu& bed,
        const lsdtt_array2d_u8& node_state,
        const std::vector<double>& rainfall_rate,
        const std::vector<double>& friction_coefficient,
        const std::vector<double>& point_inflow,
        const std::vector<double>& initial_water_depth,
        const std::vector<double>& prescribed_qin)
    {
        if (!meta_.is_valid() || bed.empty())
            throw std::invalid_argument("run_graphflood_optimised: invalid bed raster");
        if (n_cells_ > std::numeric_limits<std::uint32_t>::max())
            throw std::invalid_argument("run_graphflood_optimised: raster exceeds uint32 indexing");
        if (node_state.n_rows() != n_rows_ || node_state.n_cols() != n_cols_)
            throw std::invalid_argument("run_graphflood_optimised: node_state shape mismatch");
        if (rainfall_rate.size() != n_cells_ ||
            friction_coefficient.size() != n_cells_ || point_inflow.size() != n_cells_ ||
            initial_water_depth.size() != n_cells_ ||
            (use_v2_ && !prescribed_qin.empty() && prescribed_qin.size() != n_cells_))
            throw std::invalid_argument("run_graphflood_optimised: forcing shape mismatch");
        if (options_.boundary_conditions.size() != 4)
            throw std::invalid_argument("run_graphflood_optimised: boundary conditions must be NESW");

        for (const char code : options_.boundary_conditions)
            if (code != 'o' && code != 'n' && code != 'p')
                throw std::invalid_argument("run_graphflood_optimised: boundaries use o, n, or p");
        if ((options_.boundary_conditions[0] == 'p') !=
                (options_.boundary_conditions[2] == 'p') ||
            (options_.boundary_conditions[1] == 'p') !=
                (options_.boundary_conditions[3] == 'p'))
            throw std::invalid_argument("run_graphflood_optimised: periodic boundaries must be paired");
        if (!(options_.mfd_exponent > 0.0f) ||
            !(options_.fill_slope_epsilon >= 0.0f) ||
            !(options_.min_slope >= 0.0) || !(options_.dt > 0.0) ||
            !(options_.friction_depth_exponent > 0.0) ||
            !(options_.friction_slope_exponent > 0.0) ||
            options_.dirty_tile_size <= 0 || options_.full_validation_every < 0)
            throw std::invalid_argument("run_graphflood_optimised: invalid numerical option");

        const auto& mask = bed.get_mask();
        const auto& values = bed.get_data();
        const auto states = node_state.get_vector();
        friction_.assign(n_cells_, 0.0);
        std::size_t outlets = 0;
        for (std::size_t flat = 0; flat < n_cells_; ++flat)
        {
            const std::uint8_t state = states[flat];
            if (state != 0 && state != 1 && state != 3 &&
                !(use_v2_ && state == 2))
                throw std::invalid_argument(use_v2_
                    ? "run_graphflood_optimised: node states must be 0, 1, 2, or 3"
                    : "run_graphflood_optimised: node states must be 0, 1, or 3");
            if (state == 0) continue;
            if (mask[flat] == 0 || !std::isfinite(values[flat]) ||
                !std::isfinite(rainfall_rate[flat]) || rainfall_rate[flat] < 0.0 ||
                !std::isfinite(friction_coefficient[flat]) ||
                friction_coefficient[flat] <= 0.0 ||
                !std::isfinite(point_inflow[flat]) || point_inflow[flat] < 0.0 ||
                !std::isfinite(initial_water_depth[flat]) ||
                initial_water_depth[flat] < 0.0)
                throw std::invalid_argument("run_graphflood_optimised: invalid active-cell input");
            if (state == 2 && (prescribed_qin.size() != n_cells_ ||
                               !std::isfinite(prescribed_qin[flat]) ||
                               prescribed_qin[flat] < 0.0))
                throw std::invalid_argument("run_graphflood_optimised: invalid ghost prescribed_qin");
            if (state == 3)
            {
                ++outlets;
            }
        }
        if (outlets == 0)
            throw std::invalid_argument("run_graphflood_optimised: at least one outlet is required");
    }

    bool neighbour(int row, int col, int direction,
                   std::uint32_t& receiver) const noexcept
    {
        return neighbour_offset(row, col, k_dr[direction], k_dc[direction], receiver);
    }

    bool neighbour_offset(int row, int col, int dr, int dc,
                          std::uint32_t& receiver) const noexcept
    {
        int nr = row + dr;
        int nc = col + dc;
        if (options_.boundary_conditions[0] == 'p')
        {
            if (nr < 0) nr += n_rows_;
            else if (nr >= n_rows_) nr -= n_rows_;
        }
        if (options_.boundary_conditions[1] == 'p')
        {
            if (nc < 0) nc += n_cols_;
            else if (nc >= n_cols_) nc -= n_cols_;
        }
        if (nr < 0 || nr >= n_rows_ || nc < 0 || nc >= n_cols_) return false;
        receiver = static_cast<std::uint32_t>(nr * n_cols_ + nc);
        return state_[receiver] != 0;
    }

    void initialise_v2_geometry()
    {
        flat_offsets_ = {
            -n_cols_, -1, 1, n_cols_,
            -n_cols_ - 1, -n_cols_ + 1, n_cols_ - 1, n_cols_ + 1};
        if (!use_v2_) return;

        for (int row = 1; row + 1 < n_rows_; ++row)
        {
            for (int col = 1; col + 1 < n_cols_; ++col)
            {
                const auto flat = static_cast<std::size_t>(row) * n_cols_ + col;
                if (state_[flat] == 0) continue;
                bool complete = true;
                for (int direction = 0; direction < 8; ++direction)
                    complete = complete &&
                        state_[flat + flat_offsets_[direction]] != 0;
                interior_[flat] = complete ? 1 : 0;
            }
        }

        tile_rows_ = (n_rows_ + options_.dirty_tile_size - 1) /
            options_.dirty_tile_size;
        tile_cols_ = (n_cols_ + options_.dirty_tile_size - 1) /
            options_.dirty_tile_size;
        dirty_tiles_.assign(
            static_cast<std::size_t>(tile_rows_) * tile_cols_, 0);
        check_tiles_.assign(dirty_tiles_.size(), 0);
    }

    std::size_t tile_of(std::uint32_t flat) const noexcept
    {
        const int row = static_cast<int>(flat / n_cols_);
        const int col = static_cast<int>(flat % n_cols_);
        return static_cast<std::size_t>(row / options_.dirty_tile_size) *
            tile_cols_ + col / options_.dirty_tile_size;
    }

    void validate_outlet_components()
    {
        std::fill(flood_state_.begin(), flood_state_.end(), std::int8_t{-1});
        flood_queue_.clear();
        for (std::uint32_t start : active_)
        {
            if (flood_state_[start] != -1) continue;
            bool has_outlet = false;
            std::size_t head = 0;
            flood_queue_.clear();
            flood_queue_.push_back({0.0f, start});
            flood_state_[start] = 0;
            while (head < flood_queue_.size())
            {
                const std::uint32_t flat = flood_queue_[head++].flat;
                has_outlet = has_outlet || state_[flat] == 3;
                const int row = static_cast<int>(flat / n_cols_);
                const int col = static_cast<int>(flat % n_cols_);
                for (int direction = 0; direction < 8; ++direction)
                {
                    std::uint32_t next;
                    if (!neighbour(row, col, direction, next) ||
                        flood_state_[next] != -1) continue;
                    flood_state_[next] = 0;
                    flood_queue_.push_back({0.0f, next});
                }
            }
            if (!has_outlet)
                throw std::invalid_argument(
                    "run_graphflood_optimised: active component has no outlet");
        }
    }

    void heap_push(flood_node node)
    {
        const std::size_t old_capacity = heap_.capacity();
        heap_.push_back(node);
        if (heap_.capacity() != old_capacity)
            ++performance_.scratch_reallocations;
        std::push_heap(heap_.begin(), heap_.end(), higher_flood_node{});
    }

    flood_node heap_pop()
    {
        std::pop_heap(heap_.begin(), heap_.end(), higher_flood_node{});
        const flood_node node = heap_.back();
        heap_.pop_back();
        return node;
    }

    static int radix_bucket(std::uint32_t key, std::uint32_t last) noexcept
    {
        const std::uint32_t different = key ^ last;
        if (different == 0) return 0;
#if defined(__GNUC__) || defined(__clang__)
        return 32 - __builtin_clz(different);
#else
        int bucket = 0;
        for (std::uint32_t value = different; value != 0; value >>= 1) ++bucket;
        return bucket;
#endif
    }

    void radix_clear()
    {
        for (auto& bucket : radix_buckets_) bucket.clear();
        radix_last_ = 0;
        radix_size_ = 0;
    }

    void radix_push(flood_node node)
    {
        const std::uint32_t key = ordered_float_key(node.z);
        if (key < radix_last_)
            throw std::runtime_error("run_graphflood_optimised_v2: non-monotone flood key");
        auto& bucket = radix_buckets_[radix_bucket(key, radix_last_)];
        const std::size_t old_capacity = bucket.capacity();
        bucket.push_back(node);
        if (bucket.capacity() != old_capacity)
            ++performance_.scratch_reallocations;
        ++radix_size_;
    }

    flood_node radix_pop()
    {
        if (radix_buckets_[0].empty())
        {
            int source = 1;
            while (source < 33 && radix_buckets_[source].empty()) ++source;
            if (source == 33)
                throw std::runtime_error("run_graphflood_optimised_v2: empty radix heap");

            std::uint32_t next = std::numeric_limits<std::uint32_t>::max();
            for (const flood_node& node : radix_buckets_[source])
                next = std::min(next, ordered_float_key(node.z));
            radix_last_ = next;
            auto& nodes = radix_buckets_[source];
            const std::size_t moved = nodes.size();
            for (const flood_node& node : nodes)
            {
                auto& destination = radix_buckets_[
                    radix_bucket(ordered_float_key(node.z), radix_last_)];
                const std::size_t old_capacity = destination.capacity();
                destination.push_back(node);
                if (destination.capacity() != old_capacity)
                    ++performance_.scratch_reallocations;
            }
            nodes.clear();
            performance_.radix_redistributions +=
                static_cast<std::int64_t>(moved);
        }
        flood_node node = radix_buckets_[0].back();
        radix_buckets_[0].pop_back();
        --radix_size_;
        return node;
    }

    void flood_push(flood_node node)
    {
        if (use_v2_) radix_push(node);
        else heap_push(node);
    }

    bool flood_empty() const noexcept
    {
        return use_v2_ ? radix_size_ == 0 : heap_.empty();
    }

    flood_node flood_pop()
    {
        return use_v2_ ? radix_pop() : heap_pop();
    }

    void rebuild_topology()
    {
        const auto start = clock_type::now();
        heap_.clear();
        radix_clear();
        flood_queue_.clear();
        order_.clear();
        std::fill(flood_state_.begin(), flood_state_.end(), std::int8_t{-1});

        for (std::uint32_t flat : active_)
        {
            surface_[flat] = static_cast<float>(bed_[flat] + depth_[flat]);
            flood_state_[flat] = 0;
        }
        for (std::uint32_t flat : active_)
        {
            if (state_[flat] != 3) continue;
            flood_push({surface_[flat], flat});
            flood_state_[flat] = 1;
        }

        std::size_t queue_head = 0;
        std::size_t closed_count = 0;
        const float cell_size = static_cast<float>(meta_.cell_size);

        while (!flood_empty() || queue_head < flood_queue_.size())
        {
            flood_node centre;
            if (queue_head < flood_queue_.size() &&
                (use_v2_ || flood_empty() ||
                 flood_queue_[queue_head].z <= heap_.front().z))
                centre = flood_queue_[queue_head++];
            else
                centre = flood_pop();

            if (flood_state_[centre.flat] == 2) continue;
            flood_state_[centre.flat] = 2;
            ++closed_count;

            const int row = static_cast<int>(centre.flat / n_cols_);
            const int col = static_cast<int>(centre.flat % n_cols_);
            for (int direction = 0; direction < 8; ++direction)
            {
                std::uint32_t next;
                if (!neighbour_offset(row, col, k_pf_dr[direction],
                                      k_pf_dc[direction], next) ||
                    flood_state_[next] != 0)
                    continue;

                float z = surface_[next];
                if (z <= centre.z)
                {
                    const float distance = direction % 2 == 0
                        ? cell_size : cell_size * k_sqrt2;
                    const float requested = centre.z +
                        options_.fill_slope_epsilon * distance;
                    const float representable = std::nextafter(
                        centre.z, std::numeric_limits<float>::infinity());
                    z = std::max(requested, representable);
                    surface_[next] = z;
                    if (use_v2_) flood_push({z, next});
                    else flood_queue_.push_back({z, next});
                }
                else
                {
                    flood_push({z, next});
                }
                flood_state_[next] = 1;
            }
        }

        if (closed_count != active_.size())
            throw std::runtime_error("run_graphflood_optimised: priority flood missed active cells");

        build_surface_order();
        std::fill(rank_.begin(), rank_.end(), k_no_receiver);
        for (std::size_t i = 0; i < order_.size(); ++i)
            rank_[order_[i]] = static_cast<std::uint32_t>(i);

        const double area = meta_.cell_size * meta_.cell_size;
        for (std::uint32_t flat : active_)
        {
            if (state_[flat] == 2) continue;
            const double old_depth = depth_[flat];
            depth_[flat] = std::max(0.0,
                static_cast<double>(surface_[flat]) - bed_[flat]);
            projection_volume_ += std::max(0.0, depth_[flat] - old_depth) * area;
        }
        topology_valid_ = true;
        std::fill(dirty_tiles_.begin(), dirty_tiles_.end(), 0);
        ++performance_.topology_rebuilds;
        performance_.projected_cells += static_cast<std::int64_t>(active_.size());
        performance_.projection_seconds += seconds_since(start);
        update_working_bytes();
    }

    std::uint32_t surface_key(std::uint32_t flat) const noexcept
    {
        return ordered_float_key(surface_[flat]);
    }

    void build_surface_order()
    {
        order_ = active_;
        sort_scratch_.resize(order_.size());
        for (unsigned shift = 0; shift < 32; shift += 8)
        {
            std::size_t counts[256] = {};
            for (std::uint32_t flat : order_)
                ++counts[(surface_key(flat) >> shift) & 0xffu];
            std::size_t offsets[256];
            std::size_t offset = 0;
            for (int bucket = 0; bucket < 256; ++bucket)
            {
                offsets[bucket] = offset;
                offset += counts[bucket];
            }
            for (std::uint32_t flat : order_)
            {
                const auto bucket = (surface_key(flat) >> shift) & 0xffu;
                sort_scratch_[offsets[bucket]++] = flat;
            }
            order_.swap(sort_scratch_);
        }
        std::reverse(order_.begin(), order_.end());
    }

    bool cell_topology_is_valid(std::uint32_t flat) const
    {
        const int row = static_cast<int>(flat / n_cols_);
        const int col = static_cast<int>(flat % n_cols_);
        bool has_lower = state_[flat] == 3;
        for (int direction = 0; direction < 8; ++direction)
        {
            std::uint32_t next;
            if (!neighbour(row, col, direction, next)) continue;
            if (state_[flat] == 3 || surface_[next] >= surface_[flat]) continue;
            has_lower = true;
            if (rank_[flat] == k_no_receiver || rank_[next] <= rank_[flat])
                return false;
        }
        return has_lower;
    }

    bool topology_still_valid()
    {
        const auto start = clock_type::now();
        if (!topology_valid_)
        {
            performance_.topology_check_seconds += seconds_since(start);
            return false;
        }

        for (std::uint32_t flat : active_)
            surface_[flat] = static_cast<float>(bed_[flat] + depth_[flat]);

        ++performance_.full_topology_checks;
        for (std::uint32_t flat : active_)
            if (!cell_topology_is_valid(flat))
            {
                performance_.topology_check_seconds += seconds_since(start);
                return false;
            }
        performance_.topology_check_seconds += seconds_since(start);
        return true;
    }

    bool dirty_topology_still_valid()
    {
        const auto start = clock_type::now();
        if (!topology_valid_)
        {
            performance_.topology_check_seconds += seconds_since(start);
            return false;
        }

        const bool force_full = options_.full_validation_every > 0 &&
            steps_ > 0 && steps_ % options_.full_validation_every == 0;
        if (force_full) return topology_still_valid();

        // Unchanged tiles cannot introduce a new receiver edge.
        std::fill(check_tiles_.begin(), check_tiles_.end(), 0);
        for (int tile_row = 0; tile_row < tile_rows_; ++tile_row)
        {
            for (int tile_col = 0; tile_col < tile_cols_; ++tile_col)
            {
                const std::size_t tile =
                    static_cast<std::size_t>(tile_row) * tile_cols_ + tile_col;
                if (dirty_tiles_[tile] == 0) continue;
                for (int tr = std::max(0, tile_row - 1);
                     tr <= std::min(tile_rows_ - 1, tile_row + 1); ++tr)
                    for (int tc = std::max(0, tile_col - 1);
                         tc <= std::min(tile_cols_ - 1, tile_col + 1); ++tc)
                        check_tiles_[static_cast<std::size_t>(tr) * tile_cols_ + tc] = 1;
            }
        }

        for (int tile_row = 0; tile_row < tile_rows_; ++tile_row)
        {
            for (int tile_col = 0; tile_col < tile_cols_; ++tile_col)
            {
                const std::size_t tile =
                    static_cast<std::size_t>(tile_row) * tile_cols_ + tile_col;
                if (check_tiles_[tile] == 0) continue;
                ++performance_.dirty_tiles_checked;
                const int row_begin = tile_row * options_.dirty_tile_size;
                const int row_end = std::min(n_rows_, row_begin + options_.dirty_tile_size);
                const int col_begin = tile_col * options_.dirty_tile_size;
                const int col_end = std::min(n_cols_, col_begin + options_.dirty_tile_size);
                for (int row = row_begin; row < row_end; ++row)
                    for (int col = col_begin; col < col_end; ++col)
                    {
                        const auto flat = static_cast<std::uint32_t>(row * n_cols_ + col);
                        if (state_[flat] != 0)
                            surface_[flat] = static_cast<float>(bed_[flat] + depth_[flat]);
                    }
            }
        }

        for (int tile_row = 0; tile_row < tile_rows_; ++tile_row)
        {
            for (int tile_col = 0; tile_col < tile_cols_; ++tile_col)
            {
                const std::size_t tile =
                    static_cast<std::size_t>(tile_row) * tile_cols_ + tile_col;
                if (check_tiles_[tile] == 0) continue;
                const int row_begin = tile_row * options_.dirty_tile_size;
                const int row_end = std::min(n_rows_, row_begin + options_.dirty_tile_size);
                const int col_begin = tile_col * options_.dirty_tile_size;
                const int col_end = std::min(n_cols_, col_begin + options_.dirty_tile_size);
                for (int row = row_begin; row < row_end; ++row)
                    for (int col = col_begin; col < col_end; ++col)
                    {
                        const auto flat = static_cast<std::uint32_t>(row * n_cols_ + col);
                        if (state_[flat] != 0 && !cell_topology_is_valid(flat))
                        {
                            performance_.topology_check_seconds += seconds_since(start);
                            return false;
                        }
                    }
            }
        }
        std::fill(dirty_tiles_.begin(), dirty_tiles_.end(), 0);
        performance_.topology_check_seconds += seconds_since(start);
        return true;
    }

    void prepare_topology()
    {
        if (options_.reuse_topology &&
            (use_v2_ ? dirty_topology_still_valid() : topology_still_valid()))
        {
            ++performance_.topology_reuses;
            return;
        }
        rebuild_topology();
    }

    void mark_dirty(std::uint32_t flat)
    {
        if (!use_v2_) return;
        const float next_surface = static_cast<float>(bed_[flat] + depth_[flat]);
        if (next_surface != surface_[flat]) dirty_tiles_[tile_of(flat)] = 1;
    }

    void accumulate_mfd()
    {
        const auto start = clock_type::now();
        for (std::uint32_t flat : active_)
        {
            qin_[flat] = source_q_[flat];
            qout_[flat] = 0.0;
            dhw_[flat] = 0.0;
            slope_[flat] = 0.0;
            specific_discharge_[flat] = 0.0;
            best_receiver_[flat] = k_no_receiver;
        }

        const double cell_size = meta_.cell_size;
        const bool unit_exponent = std::abs(options_.mfd_exponent - 1.0f) < 1.0e-6f;
        for (std::uint32_t flat : order_)
        {
            if (state_[flat] == 2) qin_[flat] = prescribed_qin_[flat];
            if (state_[flat] == 3) continue;
            std::uint32_t receivers[8];
            double weights[8];
            int count = 0;
            double weight_sum = 0.0;
            double best_slope = -1.0;

            if (use_v2_ && interior_[flat] != 0)
            {
                ++performance_.interior_routed_cells;
                for (int direction = 0; direction < 8; ++direction)
                {
                    const auto next = static_cast<std::uint32_t>(
                        static_cast<std::int64_t>(flat) + flat_offsets_[direction]);
                    if (surface_[next] >= surface_[flat]) continue;
                    const double distance = direction < 4
                        ? cell_size : cell_size * k_sqrt2;
                    const double local_slope =
                        (surface_[flat] - surface_[next]) / distance;
                    double weight = unit_exponent ? local_slope
                        : std::pow(local_slope, options_.mfd_exponent);
                    weight *= distance;
                    receivers[count] = next;
                    weights[count] = weight;
                    ++count;
                    weight_sum += weight;
                    if (local_slope > best_slope)
                    {
                        best_slope = local_slope;
                        best_receiver_[flat] = next;
                        receiver_distance_[flat] = static_cast<float>(distance);
                        width_[flat] = static_cast<float>(distance);
                    }
                }
            }
            else
            {
                if (use_v2_) ++performance_.boundary_routed_cells;
                const int row = static_cast<int>(flat / n_cols_);
                const int col = static_cast<int>(flat % n_cols_);
                for (int direction = 0; direction < 8; ++direction)
                {
                    std::uint32_t next;
                    if (!neighbour(row, col, direction, next) || next == flat ||
                        surface_[next] >= surface_[flat]) continue;
                    const double distance = direction < 4
                        ? cell_size : cell_size * k_sqrt2;
                    const double local_slope =
                        (surface_[flat] - surface_[next]) / distance;
                    double weight = unit_exponent ? local_slope
                        : std::pow(local_slope, options_.mfd_exponent);
                    weight *= distance;

                    int existing = -1;
                    for (int k = 0; k < count; ++k)
                        if (receivers[k] == next) { existing = k; break; }
                    if (existing >= 0)
                        weights[existing] += weight;
                    else
                    {
                        receivers[count] = next;
                        weights[count] = weight;
                        ++count;
                    }
                    weight_sum += weight;
                    if (local_slope > best_slope)
                    {
                        best_slope = local_slope;
                        best_receiver_[flat] = next;
                        receiver_distance_[flat] = static_cast<float>(distance);
                        width_[flat] = static_cast<float>(distance);
                    }
                }
            }

            if (weight_sum > 0.0)
                for (int k = 0; k < count; ++k)
                    qin_[receivers[k]] += qin_[flat] * weights[k] / weight_sum;
            slope_[flat] = best_slope > 0.0
                ? std::max(best_slope, options_.min_slope) : 0.0;
        }

        performance_.routed_cells += static_cast<std::int64_t>(active_.size());
        performance_.routing_seconds += seconds_since(start);
    }

    double capacity(std::uint32_t flat, double depth, double slope) const noexcept
    {
        if (depth <= 0.0 || slope <= 0.0) return 0.0;
        const double depth_term = use_v2_ &&
            std::abs(options_.friction_depth_exponent - 5.0 / 3.0) < 1.0e-12
            ? depth * std::cbrt(depth * depth)
            : std::pow(depth, options_.friction_depth_exponent);
        const double slope_term = use_v2_ &&
            std::abs(options_.friction_slope_exponent - 0.5) < 1.0e-12
            ? std::sqrt(slope)
            : std::pow(slope, options_.friction_slope_exponent);
        return width_[flat] * friction_[flat] * depth_term * slope_term;
    }

    void compute_qout()
    {
        const auto start = clock_type::now();
        const double area = meta_.cell_size * meta_.cell_size;
        for (std::uint32_t flat : active_)
        {
            if (state_[flat] == 3)
            {
                qout_[flat] = qin_[flat];
                specific_discharge_[flat] = qout_[flat] / meta_.cell_size;
                continue;
            }
            const double q_capacity = capacity(flat, depth_[flat], slope_[flat]);
            const double q_available = qin_[flat] + depth_[flat] * area / options_.dt;
            qout_[flat] = std::min(q_capacity, q_available);
            specific_discharge_[flat] = width_[flat] > 0.0f
                ? qout_[flat] / width_[flat] : 0.0;
        }
        performance_.friction_seconds += seconds_since(start);
    }

    void explicit_depth_update()
    {
        const auto start = clock_type::now();
        const double scale = options_.dt /
            (meta_.cell_size * meta_.cell_size);
        for (std::uint32_t flat : active_)
        {
            const double old_depth = depth_[flat];
            depth_[flat] = std::max(0.0,
                old_depth + scale * (qin_[flat] - qout_[flat]));
            dhw_[flat] = depth_[flat] - old_depth;
            mark_dirty(flat);
        }
        performance_.depth_update_seconds += seconds_since(start);
    }

    void run_one_step()
    {
        prepare_topology();
        accumulate_mfd();
        compute_qout();
        explicit_depth_update();
    }

    std::pair<double, double> convergence_medians()
    {
        convergence_abs_.clear();
        convergence_rel_.clear();
        for (std::uint32_t flat : active_)
        {
            const double value = std::abs(dhw_[flat]);
            convergence_abs_.push_back(value);
            convergence_rel_.push_back(value / std::max(depth_[flat], 1.0e-12));
        }
        if (convergence_abs_.empty()) return {0.0, 0.0};
        return {compute_median_in_place(convergence_abs_),
                compute_median_in_place(convergence_rel_)};
    }

    lsdtt_raster_result raster_from_float(const std::vector<float>& values) const
    {
        lsdtt_raster_result result;
        result.meta = meta_;
        result.data.assign(n_cells_, static_cast<float>(meta_.no_data_value));
        result.mask.assign(n_cells_, 0);
        for (std::uint32_t flat : active_)
        {
            result.data[flat] = values[flat];
            result.mask[flat] = 1;
        }
        return result;
    }

    lsdtt_raster_result raster_from_double(const std::vector<double>& values) const
    {
        lsdtt_raster_result result;
        result.meta = meta_;
        result.data.assign(n_cells_, static_cast<float>(meta_.no_data_value));
        result.mask.assign(n_cells_, 0);
        for (std::uint32_t flat : active_)
        {
            result.data[flat] = static_cast<float>(values[flat]);
            result.mask[flat] = 1;
        }
        return result;
    }

    void update_working_bytes()
    {
        performance_.persistent_working_bytes =
            bed_.capacity() * sizeof(float) + surface_.capacity() * sizeof(float) +
            state_.capacity() * sizeof(std::uint8_t) +
            interior_.capacity() * sizeof(std::uint8_t) +
            dirty_tiles_.capacity() * sizeof(std::uint8_t) +
            check_tiles_.capacity() * sizeof(std::uint8_t) +
            flood_state_.capacity() * sizeof(std::int8_t) +
            active_.capacity() * sizeof(std::uint32_t) +
            order_.capacity() * sizeof(std::uint32_t) +
            sort_scratch_.capacity() * sizeof(std::uint32_t) +
            rank_.capacity() * sizeof(std::uint32_t) +
            best_receiver_.capacity() * sizeof(std::uint32_t) +
            receiver_distance_.capacity() * sizeof(float) +
            width_.capacity() * sizeof(float) +
            depth_.capacity() * sizeof(double) + qin_.capacity() * sizeof(double) +
            qout_.capacity() * sizeof(double) + dhw_.capacity() * sizeof(double) +
            source_q_.capacity() * sizeof(double) + friction_.capacity() * sizeof(double) +
            prescribed_qin_.capacity() * sizeof(double) +
            slope_.capacity() * sizeof(double) +
            specific_discharge_.capacity() * sizeof(double) +
            heap_.capacity() * sizeof(flood_node) +
            flood_queue_.capacity() * sizeof(flood_node);
        for (const auto& bucket : radix_buckets_)
            performance_.persistent_working_bytes +=
                bucket.capacity() * sizeof(flood_node);
    }

    graphflood_optimised_result make_result(bool converged)
    {
        std::vector<float> final_surface(n_cells_, static_cast<float>(meta_.no_data_value));
        for (std::uint32_t flat : active_)
            final_surface[flat] = static_cast<float>(bed_[flat] + depth_[flat]);

        graphflood_optimised_result result;
        result.bed = raster_from_float(bed_);
        result.water_depth = raster_from_double(depth_);
        result.water_surface = raster_from_float(final_surface);
        result.qin = raster_from_double(qin_);
        result.qout = raster_from_double(qout_);
        result.dhw = raster_from_double(dhw_);
        result.hydraulic_slope = raster_from_double(slope_);
        result.specific_discharge = raster_from_double(specific_discharge_);
        result.total_discharge = raster_from_double(qin_);
        result.performance = performance_;
        result.converged = converged;
        result.step_count = steps_;
        result.projection_volume_m3 = projection_volume_;
        const double area = meta_.cell_size * meta_.cell_size;
        for (std::uint32_t flat : active_)
            result.total_water_volume_m3 += depth_[flat] * area;
        return result;
    }

    lsdtt_raster_meta meta_;
    int n_rows_ = 0;
    int n_cols_ = 0;
    std::size_t n_cells_ = 0;
    graphflood_optimised_options options_;
    std::vector<float> bed_, surface_, width_, receiver_distance_;
    std::vector<std::uint8_t> state_, interior_, dirty_tiles_, check_tiles_;
    std::vector<std::int8_t> flood_state_;
    std::vector<std::uint32_t> active_, order_, sort_scratch_, best_receiver_;
    std::vector<std::uint32_t> rank_;
    std::vector<double> depth_, qin_, qout_, dhw_, source_q_, friction_, prescribed_qin_;
    std::vector<double> slope_, specific_discharge_;
    std::vector<double> convergence_abs_, convergence_rel_;
    std::vector<flood_node> heap_, flood_queue_;
    std::array<std::vector<flood_node>, 33> radix_buckets_;
    std::array<std::int32_t, 8> flat_offsets_{};
    bool topology_valid_ = false;
    bool use_v2_ = false;
    int tile_rows_ = 0;
    int tile_cols_ = 0;
    std::uint32_t radix_last_ = 0;
    std::size_t radix_size_ = 0;
    std::int64_t steps_ = 0;
    double projection_volume_ = 0.0;
    graphflood_optimised_performance performance_;
};

} // namespace

graphflood_optimised_result run_graphflood_optimised_cpu(
    const lsdtt_raster_cpu& bed,
    const lsdtt_array2d_u8& node_state,
    const std::vector<double>& rainfall_rate,
    const std::vector<double>& friction_coefficient,
    const std::vector<double>& point_inflow,
    const graphflood_optimised_options& options,
    std::int64_t max_steps,
    bool run_until_converged,
    int convergence_check_every,
    double tol_median_abs_dhw,
    double tol_median_rel_dhw)
{
    std::vector<double> initial_water_depth(point_inflow.size(), 0.0);
    return run_graphflood_optimised_cpu(
        bed, node_state, rainfall_rate, friction_coefficient, point_inflow,
        initial_water_depth, options, max_steps, run_until_converged,
        convergence_check_every, tol_median_abs_dhw, tol_median_rel_dhw);
}

graphflood_optimised_result run_graphflood_optimised_cpu(
    const lsdtt_raster_cpu& bed,
    const lsdtt_array2d_u8& node_state,
    const std::vector<double>& rainfall_rate,
    const std::vector<double>& friction_coefficient,
    const std::vector<double>& point_inflow,
    const std::vector<double>& initial_water_depth,
    const graphflood_optimised_options& options,
    std::int64_t max_steps,
    bool run_until_converged,
    int convergence_check_every,
    double tol_median_abs_dhw,
    double tol_median_rel_dhw)
{
    optimised_solver solver(bed, node_state, rainfall_rate,
                            friction_coefficient, point_inflow,
                            initial_water_depth, {}, options, false);
    return solver.run(max_steps, run_until_converged, convergence_check_every,
                      tol_median_abs_dhw, tol_median_rel_dhw);
}

graphflood_optimised_result run_graphflood_optimised_v2_cpu(
    const lsdtt_raster_cpu& bed,
    const lsdtt_array2d_u8& node_state,
    const std::vector<double>& rainfall_rate,
    const std::vector<double>& friction_coefficient,
    const std::vector<double>& point_inflow,
    const std::vector<double>& initial_water_depth,
    const std::vector<double>& prescribed_qin,
    const graphflood_optimised_options& options,
    std::int64_t max_steps,
    bool run_until_converged,
    int convergence_check_every,
    double tol_median_abs_dhw,
    double tol_median_rel_dhw)
{
    optimised_solver solver(bed, node_state, rainfall_rate,
                            friction_coefficient, point_inflow,
                            initial_water_depth, prescribed_qin, options, true);
    return solver.run(max_steps, run_until_converged, convergence_check_every,
                      tol_median_abs_dhw, tol_median_rel_dhw);
}

graphflood_multigrid_result run_graphflood_optimised_v2_multigrid_cpu(
    const lsdtt_raster_cpu& bed,
    const lsdtt_array2d_u8& node_state,
    const std::vector<double>& rainfall_rate,
    const std::vector<double>& friction_coefficient,
    const std::vector<double>& point_inflow,
    const std::vector<double>& initial_water_depth,
    const graphflood_optimised_options& options,
    const graphflood_multigrid_transfer_policy& transfer_policy,
    const std::vector<int>& factors,
    const std::vector<std::int64_t>& steps_per_level,
    bool run_until_converged,
    int convergence_check_every,
    double tol_median_abs_dhw,
    double tol_median_rel_dhw)
{
    if (factors.empty() || factors.back() != 1 || factors.size() != steps_per_level.size())
        throw std::invalid_argument("multigrid factors must end at one and match steps_per_level");
    const auto& fine_meta = bed.get_meta();
    const auto n_fine = static_cast<std::size_t>(fine_meta.n_pixels());
    if (rainfall_rate.size() != n_fine || friction_coefficient.size() != n_fine ||
        point_inflow.size() != n_fine || initial_water_depth.size() != n_fine)
        throw std::invalid_argument("multigrid forcing shape mismatch");

    const auto make_raster = [&fine_meta](const std::vector<double>& values,
                                          const std::vector<std::uint8_t>& mask) {
        std::vector<float> data(values.size());
        for (std::size_t i = 0; i < values.size(); ++i) data[i] = static_cast<float>(values[i]);
        return lsdtt_raster_cpu(std::move(data), std::vector<std::uint8_t>(mask), fine_meta);
    };
    const auto fine_states = node_state.get_vector();
    std::vector<std::uint8_t> active(n_fine, 0);
    for (std::size_t i = 0; i < n_fine; ++i) active[i] = fine_states[i] != 0;
    std::vector<double> rainfall_transfer(n_fine, 0.0);
    for (std::size_t i = 0; i < n_fine; ++i)
        if (active[i] != 0) rainfall_transfer[i] = rainfall_rate[i];
    const auto rain_source = make_raster(rainfall_transfer, std::vector<std::uint8_t>(n_fine, 1));
    const auto friction_source = make_raster(friction_coefficient, active);
    const auto depth_source = make_raster(initial_water_depth, active);

    graphflood_multigrid_result output;
    output.factors = factors;
    std::optional<lsdtt_raster_cpu> previous_depth;
    for (std::size_t level = 0; level < factors.size(); ++level)
    {
        const int factor = factors[level];
        if (factor <= 0 || fine_meta.n_rows % factor != 0 || fine_meta.n_cols % factor != 0 ||
            (level > 0 && (factors[level - 1] <= factor || factors[level - 1] % factor != 0)))
            throw std::invalid_argument("multigrid factors must be descending divisors of the fine grid");

        auto meta = fine_meta;
        meta.n_rows /= factor;
        meta.n_cols /= factor;
        meta.cell_size *= factor;
        const auto n_cells = static_cast<std::size_t>(meta.n_pixels());
        std::vector<std::uint8_t> coarse_states(n_cells, 0), coarse_mask(n_cells, 0);
        for (int row = 0; row < fine_meta.n_rows; ++row)
            for (int col = 0; col < fine_meta.n_cols; ++col)
            {
                const auto fine = static_cast<std::size_t>(row) * fine_meta.n_cols + col;
                const auto coarse = static_cast<std::size_t>(row / factor) * meta.n_cols + col / factor;
                if (fine_states[fine] == 0) continue;
                coarse_mask[coarse] = 1;
                coarse_states[coarse] = fine_states[fine] == 3 ? 3 :
                    (coarse_states[coarse] == 3 ? 3 : 1);
            }
        lsdtt_array2d_u8 states(meta.n_rows, meta.n_cols, 0);
        for (int row = 0; row < meta.n_rows; ++row)
            for (int col = 0; col < meta.n_cols; ++col)
                states(row, col) = coarse_states[static_cast<std::size_t>(row) * meta.n_cols + col];

        const auto resample = [&meta](const lsdtt_raster_cpu& source,
                                      raster_resampling_method method) {
            const auto n = static_cast<std::size_t>(meta.n_pixels());
            lsdtt_raster_cpu destination(
                std::vector<float>(n, static_cast<float>(meta.no_data_value)),
                std::vector<std::uint8_t>(n, 0), meta);
            resample_raster(source, destination, method);
            return destination;
        };
        auto level_bed = factor == 1 ? bed : resample(bed, transfer_policy.static_field);
        auto level_rain = factor == 1 ? rain_source : resample(rain_source, transfer_policy.static_field);
        auto level_friction = factor == 1 ? friction_source : resample(friction_source, transfer_policy.static_field);
        const auto with_mask = [&coarse_mask](const lsdtt_raster_cpu& raster) {
            return lsdtt_raster_cpu(std::vector<float>(raster.get_data()),
                                    std::vector<std::uint8_t>(coarse_mask),
                                    raster.get_meta());
        };
        level_bed = with_mask(level_bed);
        level_rain = with_mask(level_rain);
        level_friction = with_mask(level_friction);

        lsdtt_raster_cpu level_depth = previous_depth
            ? resample(*previous_depth, transfer_policy.depth_prolong)
            : (factor == 1 ? depth_source : resample(depth_source, transfer_policy.depth_restrict));
        level_depth = with_mask(level_depth);
        std::vector<double> level_points(n_cells, 0.0);
        for (int row = 0; row < fine_meta.n_rows; ++row)
            for (int col = 0; col < fine_meta.n_cols; ++col)
                level_points[static_cast<std::size_t>(row / factor) * meta.n_cols + col / factor] +=
                    point_inflow[static_cast<std::size_t>(row) * fine_meta.n_cols + col];

        const auto to_double = [](const lsdtt_raster_cpu& raster) {
            const auto& values = raster.get_data();
            return std::vector<double>(values.begin(), values.end());
        };
        auto result = run_graphflood_optimised_v2_cpu(
            level_bed, states, to_double(level_rain), to_double(level_friction), level_points,
            to_double(level_depth), std::vector<double>(n_cells, 0.0), options,
            steps_per_level[level], run_until_converged,
            convergence_check_every, tol_median_abs_dhw, tol_median_rel_dhw);
        previous_depth.emplace(std::vector<float>(result.water_depth.data),
                               std::vector<std::uint8_t>(result.water_depth.mask),
                               result.water_depth.meta);
        output.level_results.push_back(std::move(result));
    }
    return output;
}

} // namespace lsdtt
