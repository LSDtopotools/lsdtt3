// lsdtt_synthetic_network_cpu.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
//
// See lsdtt_synthetic_network_cpu.hpp for the algorithm description.
// -----------------------------------------------------------------------------

#include "lsdtt_synthetic_network_cpu.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace lsdtt {

// D8 neighbour offsets (row, col)
static const int k_dr[8] = {-1, -1, -1,  0, 0,  1, 1, 1};
static const int k_dc[8] = {-1,  0,  1, -1, 1, -1, 0, 1};

// -----------------------------------------------------------------------------
lsdtt_synthetic_network_cpu::lsdtt_synthetic_network_cpu(
    const lsdtt_ocn_options& opts)
    : opts_(opts)
{
    if (opts_.n_rows < 2 || opts_.n_cols < 2)
    {
        throw std::runtime_error(
            "lsdtt_synthetic_network_cpu: n_rows and n_cols must be >= 2");
    }
    if (!(opts_.cell_size > 0.0) || !std::isfinite(opts_.cell_size))
    {
        throw std::runtime_error(
            "lsdtt_synthetic_network_cpu: cell_size must be finite and > 0");
    }
    if (!(opts_.gamma > 0.0) || opts_.gamma >= 1.0)
    {
        throw std::runtime_error(
            "lsdtt_synthetic_network_cpu: gamma must be in (0, 1)");
    }

    nrows_ = opts_.n_rows;
    ncols_ = opts_.n_cols;
    n_     = static_cast<std::int64_t>(nrows_) * static_cast<std::int64_t>(ncols_);

    int out_col = opts_.outlet_col;
    if (out_col < 0 || out_col >= ncols_) out_col = ncols_ / 2;
    outlet_ = static_cast<std::int64_t>(nrows_ - 1) * ncols_ + out_col;

    rng_.seed(opts_.seed);

    receiver_.assign(n_, -1);
    area_.assign(n_, 1.0);
    mark_.assign(n_, 0);

    init_comb_tree_();
    recompute_all_areas_();
    energy_  = total_energy_();
    energy0_ = energy_;
}

// -----------------------------------------------------------------------------
lsdtt_raster_meta lsdtt_synthetic_network_cpu::build_meta_() const
{
    lsdtt_raster_meta meta;
    meta.n_rows        = nrows_;
    meta.n_cols        = ncols_;
    meta.x_min         = 0.0;
    meta.y_min         = 0.0;
    meta.cell_size     = opts_.cell_size;
    meta.no_data_value = -9999.0;
    meta.data_type     = GDT_Float32;
    meta.crs_wkt       = "";
    return meta;
}

// -----------------------------------------------------------------------------
// Initial valid spanning tree: a "comb". Every pixel above the bottom row drains
// straight south; the bottom row drains horizontally toward the outlet column.
// -----------------------------------------------------------------------------
void lsdtt_synthetic_network_cpu::init_comb_tree_()
{
    const int out_col = static_cast<int>(outlet_ % ncols_);

    for (int r = 0; r < nrows_; ++r)
    {
        for (int c = 0; c < ncols_; ++c)
        {
            const std::int64_t f = static_cast<std::int64_t>(r) * ncols_ + c;

            if (f == outlet_)
            {
                receiver_[f] = outlet_;   // outlet drains to itself (sentinel)
            }
            else if (r < nrows_ - 1)
            {
                receiver_[f] = static_cast<std::int64_t>(r + 1) * ncols_ + c;
            }
            else // bottom row, not the outlet: step toward the outlet column
            {
                const int nc = (c < out_col) ? c + 1 : c - 1;
                receiver_[f] = static_cast<std::int64_t>(r) * ncols_ + nc;
            }
        }
    }
}

// -----------------------------------------------------------------------------
// O(N) flow accumulation over the tree (Kahn). Records topo_order_ with sources
// first and the outlet last.
// -----------------------------------------------------------------------------
void lsdtt_synthetic_network_cpu::recompute_all_areas_()
{
    std::vector<std::int32_t> indeg(n_, 0);
    for (std::int64_t f = 0; f < n_; ++f)
    {
        if (f != outlet_) indeg[receiver_[f]] += 1;
    }

    std::fill(area_.begin(), area_.end(), 1.0);
    topo_order_.clear();
    topo_order_.reserve(static_cast<std::size_t>(n_));

    std::vector<std::int64_t> stack;
    stack.reserve(static_cast<std::size_t>(n_));
    for (std::int64_t f = 0; f < n_; ++f)
        if (indeg[f] == 0) stack.push_back(f);

    while (!stack.empty())
    {
        const std::int64_t f = stack.back();
        stack.pop_back();
        topo_order_.push_back(f);

        if (f == outlet_) continue;

        const std::int64_t r = receiver_[f];
        area_[r] += area_[f];
        if (--indeg[r] == 0) stack.push_back(r);
    }
}

// -----------------------------------------------------------------------------
double lsdtt_synthetic_network_cpu::total_energy_() const
{
    const double g = opts_.gamma;
    double e = 0.0;
    for (std::int64_t f = 0; f < n_; ++f) e += std::pow(area_[f], g);
    return e;
}

// -----------------------------------------------------------------------------
double lsdtt_synthetic_network_cpu::receiver_length_(std::int64_t node) const
{
    const std::int64_t r = receiver_[node];
    const int dr = std::abs(static_cast<int>(node / ncols_) -
                            static_cast<int>(r / ncols_));
    const int dc = std::abs(static_cast<int>(node % ncols_) -
                            static_cast<int>(r % ncols_));
    const double diag = (dr + dc == 2);
    return diag ? opts_.cell_size * std::sqrt(2.0) : opts_.cell_size;
}

// -----------------------------------------------------------------------------
// True if reconnecting node -> cand_receiver would create a cycle, i.e. if
// cand_receiver lies in the subtree of node (following receivers from
// cand_receiver reaches node before the outlet).
// -----------------------------------------------------------------------------
bool lsdtt_synthetic_network_cpu::creates_cycle_(
    std::int64_t node, std::int64_t cand_receiver) const
{
    std::int64_t x = cand_receiver;
    while (true)
    {
        if (x == node)   return true;
        if (x == outlet_) return false;
        x = receiver_[x];
    }
}

// -----------------------------------------------------------------------------
double lsdtt_synthetic_network_cpu::delta_energy_(
    std::int64_t node, std::int64_t cand_receiver)
{
    const double g   = opts_.gamma;
    const double sub = area_[node];
    const std::int64_t old_r = receiver_[node];

    // mark the old downstream path (old_r ... outlet, inclusive)
    ++stamp_;
    for (std::int64_t x = old_r; ; x = receiver_[x])
    {
        mark_[x] = stamp_;
        if (x == outlet_) break;
    }

    // new path: walk until the confluence (first node already on the old path)
    gainers_.clear();
    std::int64_t confluence = cand_receiver;
    while (mark_[confluence] != stamp_)
    {
        gainers_.push_back(confluence);
        confluence = receiver_[confluence];
    }

    // old-path losers: old_r ... up to (not including) the confluence
    losers_.clear();
    for (std::int64_t x = old_r; x != confluence; x = receiver_[x])
        losers_.push_back(x);

    double dE = 0.0;
    for (const std::int64_t nidx : losers_)
        dE += std::pow(area_[nidx] - sub, g) - std::pow(area_[nidx], g);
    for (const std::int64_t nidx : gainers_)
        dE += std::pow(area_[nidx] + sub, g) - std::pow(area_[nidx], g);

    return dE;
}

// -----------------------------------------------------------------------------
void lsdtt_synthetic_network_cpu::commit_(
    std::int64_t node, std::int64_t cand_receiver, double dE)
{
    const double sub = area_[node];
    for (const std::int64_t nidx : losers_)  area_[nidx] -= sub;
    for (const std::int64_t nidx : gainers_) area_[nidx] += sub;
    receiver_[node] = cand_receiver;
    energy_ += dE;
}

// -----------------------------------------------------------------------------
void lsdtt_synthetic_network_cpu::run()
{
    const long long n_iter =
        (opts_.n_iter > 0) ? opts_.n_iter : 40LL * n_;

    double T0 = opts_.initial_temperature;
    if (!(T0 > 0.0)) T0 = opts_.temperature_scale * energy0_;

    const long long n0 = static_cast<long long>(
        std::llround(opts_.initial_no_cooling_fraction * n_iter));

    std::uniform_int_distribution<std::int64_t> pick_node(0, n_ - 1);
    std::uniform_int_distribution<int>           pick_dir(0, 7);
    std::uniform_real_distribution<double>       unit(0.0, 1.0);

    for (long long it = 0; it < n_iter; ++it)
    {
        // temperature for this iteration
        double T;
        if (it < n0)
        {
            T = T0;
        }
        else
        {
            const double k = static_cast<double>(it - n0);
            T = T0 * std::exp(-opts_.cooling_rate * k /
                              static_cast<double>(n_));
        }

        // propose: a random non-outlet node and a random in-bounds neighbour
        const std::int64_t node = pick_node(rng_);
        if (node == outlet_) continue;

        const int dir = pick_dir(rng_);
        const int nr = static_cast<int>(node / ncols_) + k_dr[dir];
        const int nc = static_cast<int>(node % ncols_) + k_dc[dir];
        if (nr < 0 || nr >= nrows_ || nc < 0 || nc >= ncols_) continue;

        const std::int64_t cand = static_cast<std::int64_t>(nr) * ncols_ + nc;
        if (cand == receiver_[node]) continue;
        if (creates_cycle_(node, cand))  continue;

        const double dE = delta_energy_(node, cand);

        bool accept = (dE <= 0.0);
        if (!accept && T > 0.0)
        {
            accept = (unit(rng_) < std::exp(-dE / T));
        }

        if (accept) commit_(node, cand, dE);
    }

    // refresh exact areas + topological order for the elevation pass
    recompute_all_areas_();
    energy_ = total_energy_();
}

// -----------------------------------------------------------------------------
lsdtt_raster_result lsdtt_synthetic_network_cpu::get_drainage_area_cells_raster() const
{
    lsdtt_raster_result out;
    out.meta = build_meta_();
    out.data.resize(static_cast<std::size_t>(n_));
    out.mask.assign(static_cast<std::size_t>(n_), 1);
    for (std::int64_t f = 0; f < n_; ++f)
        out.data[static_cast<std::size_t>(f)] = static_cast<float>(area_[f]);
    return out;
}

// -----------------------------------------------------------------------------
// Slope-area elevation: S_i = slope0 * (A_i / N)^(gamma-1),
// Z_i = Z_receiver(i) + S_i * length_i, integrated upstream from the outlet.
// -----------------------------------------------------------------------------
lsdtt_raster_result lsdtt_synthetic_network_cpu::get_elevation_raster() const
{
    const double g    = opts_.gamma;
    const double invN = 1.0 / static_cast<double>(n_);

    std::vector<float> z(static_cast<std::size_t>(n_), 0.0f);

    // topo_order_ is sources-first, outlet-last; iterate in reverse so each
    // node's receiver elevation is already set.
    for (auto it = topo_order_.rbegin(); it != topo_order_.rend(); ++it)
    {
        const std::int64_t node = *it;
        if (node == outlet_)
        {
            z[static_cast<std::size_t>(node)] = static_cast<float>(opts_.z_min);
            continue;
        }
        const double slope = opts_.slope0 *
                             std::pow(area_[node] * invN, g - 1.0);
        const double dz = slope * receiver_length_(node);
        z[static_cast<std::size_t>(node)] =
            z[static_cast<std::size_t>(receiver_[node])] + static_cast<float>(dz);
    }

    lsdtt_raster_result out;
    out.meta = build_meta_();
    out.data = std::move(z);
    out.mask.assign(static_cast<std::size_t>(n_), 1);
    return out;
}

} // namespace lsdtt
