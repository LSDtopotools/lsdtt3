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
#include <queue>
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
        throw std::runtime_error(
            "lsdtt_synthetic_network_cpu: n_rows and n_cols must be >= 2");
    if (!(opts_.cell_size > 0.0) || !std::isfinite(opts_.cell_size))
        throw std::runtime_error(
            "lsdtt_synthetic_network_cpu: cell_size must be finite and > 0");
    if (!(opts_.gamma > 0.0) || opts_.gamma >= 1.0)
        throw std::runtime_error(
            "lsdtt_synthetic_network_cpu: gamma must be in (0, 1)");

    nrows_ = opts_.n_rows;
    ncols_ = opts_.n_cols;
    n_     = static_cast<std::int64_t>(nrows_) * static_cast<std::int64_t>(ncols_);

    rng_.seed(opts_.seed);

    receiver_.assign(n_, -1);
    area_.assign(n_, 0.0);
    mark_.assign(n_, 0);

    setup_domain_and_outlets_();
    init_forest_();
    recompute_all_areas_();
    energy_  = total_energy_();
    energy0_ = energy_;
}

// -----------------------------------------------------------------------------
lsdtt_raster_meta lsdtt_synthetic_network_cpu::build_meta_() const
{
    if (opts_.use_supplied_meta) return opts_.supplied_meta;

    lsdtt_raster_meta meta;
    meta.n_rows        = nrows_;
    meta.n_cols        = ncols_;
    meta.cell_size     = opts_.cell_size;
    meta.no_data_value = -9999.0;
    meta.data_type     = GDT_Float32;
    meta.set_default_synthetic_georeference();
    return meta;
}

// -----------------------------------------------------------------------------
void lsdtt_synthetic_network_cpu::setup_domain_and_outlets_()
{
    land_.assign(static_cast<std::size_t>(n_), 1);
    if (!opts_.domain_mask.empty())
    {
        if (static_cast<std::int64_t>(opts_.domain_mask.size()) != n_)
            throw std::runtime_error(
                "lsdtt_synthetic_network_cpu: domain_mask size mismatch");
        for (std::int64_t f = 0; f < n_; ++f)
            land_[static_cast<std::size_t>(f)] =
                opts_.domain_mask[static_cast<std::size_t>(f)] ? 1 : 0;
    }

    is_outlet_.assign(static_cast<std::size_t>(n_), 0);
    outlets_.clear();

    std::int64_t n_land = 0;
    for (std::int64_t f = 0; f < n_; ++f)
        if (land_[static_cast<std::size_t>(f)]) ++n_land;
    if (n_land == 0)
        throw std::runtime_error(
            "lsdtt_synthetic_network_cpu: domain has no land pixels");

    auto land_at = [&](int r, int c) -> bool
    {
        if (r < 0 || r >= nrows_ || c < 0 || c >= ncols_) return false;
        return land_[static_cast<std::size_t>(r) * ncols_ + c] != 0;
    };

    if (opts_.outlet_mode == ocn_outlet_mode::all_coast)
    {
        // every land pixel adjacent (D8) to sea or the grid edge is an outlet
        for (int r = 0; r < nrows_; ++r)
        {
            for (int c = 0; c < ncols_; ++c)
            {
                if (!land_at(r, c)) continue;
                bool coastal = false;
                for (int k = 0; k < 8 && !coastal; ++k)
                    if (!land_at(r + k_dr[k], c + k_dc[k])) coastal = true;
                if (coastal)
                {
                    const std::int64_t f =
                        static_cast<std::int64_t>(r) * ncols_ + c;
                    is_outlet_[static_cast<std::size_t>(f)] = 1;
                    outlets_.push_back(f);
                }
            }
        }
    }
    else // single_south: southmost land pixel near the target column
    {
        int target_col = opts_.outlet_col;
        if (target_col < 0 || target_col >= ncols_) target_col = ncols_ / 2;

        std::int64_t best = -1;
        int best_row = -1, best_dc = ncols_ + 1;
        for (int r = 0; r < nrows_; ++r)
        {
            for (int c = 0; c < ncols_; ++c)
            {
                if (!land_at(r, c)) continue;
                const int dc = std::abs(c - target_col);
                if (r > best_row || (r == best_row && dc < best_dc))
                {
                    best_row = r; best_dc = dc;
                    best = static_cast<std::int64_t>(r) * ncols_ + c;
                }
            }
        }
        if (best >= 0)
        {
            is_outlet_[static_cast<std::size_t>(best)] = 1;
            outlets_.push_back(best);
        }
    }

    outlet_ = outlets_.empty() ? -1 : outlets_[0];
}

// -----------------------------------------------------------------------------
// Multi-source BFS from all outlets over land: each land pixel's receiver is the
// neighbour from which it was first reached, giving a valid spanning forest.
// -----------------------------------------------------------------------------
void lsdtt_synthetic_network_cpu::init_forest_()
{
    std::fill(receiver_.begin(), receiver_.end(), static_cast<std::int64_t>(-1));
    std::vector<char> visited(static_cast<std::size_t>(n_), 0);
    std::queue<std::int64_t> q;

    for (const std::int64_t o : outlets_)
    {
        receiver_[static_cast<std::size_t>(o)] = o;   // outlet drains to itself
        visited[static_cast<std::size_t>(o)] = 1;
        q.push(o);
    }

    while (!q.empty())
    {
        const std::int64_t f = q.front();
        q.pop();
        const int r = static_cast<int>(f / ncols_);
        const int c = static_cast<int>(f % ncols_);
        for (int k = 0; k < 8; ++k)
        {
            const int nr = r + k_dr[k];
            const int nc = c + k_dc[k];
            if (nr < 0 || nr >= nrows_ || nc < 0 || nc >= ncols_) continue;
            const std::int64_t nf = static_cast<std::int64_t>(nr) * ncols_ + nc;
            if (!land_[static_cast<std::size_t>(nf)] ||
                visited[static_cast<std::size_t>(nf)]) continue;
            visited[static_cast<std::size_t>(nf)] = 1;
            receiver_[static_cast<std::size_t>(nf)] = f;
            q.push(nf);
        }
    }

    // Land not reachable from any outlet (disconnected) becomes its own outlet.
    for (std::int64_t f = 0; f < n_; ++f)
    {
        if (land_[static_cast<std::size_t>(f)] &&
            !visited[static_cast<std::size_t>(f)])
        {
            receiver_[static_cast<std::size_t>(f)] = f;
            is_outlet_[static_cast<std::size_t>(f)] = 1;
            outlets_.push_back(f);
        }
    }
    if (outlet_ < 0 && !outlets_.empty()) outlet_ = outlets_[0];
}

// -----------------------------------------------------------------------------
void lsdtt_synthetic_network_cpu::recompute_all_areas_()
{
    std::vector<std::int32_t> indeg(static_cast<std::size_t>(n_), 0);
    for (std::int64_t f = 0; f < n_; ++f)
    {
        if (land_[static_cast<std::size_t>(f)] &&
            !is_outlet_[static_cast<std::size_t>(f)] &&
            receiver_[static_cast<std::size_t>(f)] >= 0)
        {
            indeg[static_cast<std::size_t>(receiver_[static_cast<std::size_t>(f)])] += 1;
        }
    }

    std::fill(area_.begin(), area_.end(), 0.0);
    for (std::int64_t f = 0; f < n_; ++f)
        if (land_[static_cast<std::size_t>(f)]) area_[static_cast<std::size_t>(f)] = 1.0;

    topo_order_.clear();
    topo_order_.reserve(static_cast<std::size_t>(n_));

    std::vector<std::int64_t> stack;
    for (std::int64_t f = 0; f < n_; ++f)
        if (land_[static_cast<std::size_t>(f)] && indeg[static_cast<std::size_t>(f)] == 0)
            stack.push_back(f);

    while (!stack.empty())
    {
        const std::int64_t f = stack.back();
        stack.pop_back();
        topo_order_.push_back(f);

        if (is_outlet_[static_cast<std::size_t>(f)]) continue;

        const std::int64_t r = receiver_[static_cast<std::size_t>(f)];
        if (r < 0) continue;
        area_[static_cast<std::size_t>(r)] += area_[static_cast<std::size_t>(f)];
        if (--indeg[static_cast<std::size_t>(r)] == 0) stack.push_back(r);
    }
}

// -----------------------------------------------------------------------------
double lsdtt_synthetic_network_cpu::total_energy_() const
{
    const double g = opts_.gamma;
    double e = 0.0;
    for (std::int64_t f = 0; f < n_; ++f)
        if (land_[static_cast<std::size_t>(f)])
            e += std::pow(area_[static_cast<std::size_t>(f)], g);
    return e;
}

// -----------------------------------------------------------------------------
double lsdtt_synthetic_network_cpu::receiver_length_(std::int64_t node) const
{
    const std::int64_t r = receiver_[static_cast<std::size_t>(node)];
    const int dr = std::abs(static_cast<int>(node / ncols_) -
                            static_cast<int>(r / ncols_));
    const int dc = std::abs(static_cast<int>(node % ncols_) -
                            static_cast<int>(r % ncols_));
    const bool diag = (dr + dc == 2);
    return diag ? opts_.cell_size * std::sqrt(2.0) : opts_.cell_size;
}

// -----------------------------------------------------------------------------
bool lsdtt_synthetic_network_cpu::creates_cycle_(
    std::int64_t node, std::int64_t cand_receiver) const
{
    std::int64_t x = cand_receiver;
    while (true)
    {
        if (x == node)                                  return true;
        if (is_outlet_[static_cast<std::size_t>(x)])    return false;
        x = receiver_[static_cast<std::size_t>(x)];
        if (x < 0)                                      return false;
    }
}

// -----------------------------------------------------------------------------
double lsdtt_synthetic_network_cpu::delta_energy_(
    std::int64_t node, std::int64_t cand_receiver)
{
    const double g   = opts_.gamma;
    const double sub = area_[static_cast<std::size_t>(node)];
    const std::int64_t old_r = receiver_[static_cast<std::size_t>(node)];

    // mark the old downstream path (old_r ... its outlet, inclusive)
    ++stamp_;
    for (std::int64_t x = old_r; ; x = receiver_[static_cast<std::size_t>(x)])
    {
        mark_[static_cast<std::size_t>(x)] = stamp_;
        if (is_outlet_[static_cast<std::size_t>(x)]) break;
    }

    // new path: walk until a marked node (confluence) or its own outlet (disjoint)
    gainers_.clear();
    std::int64_t confluence = -1;
    for (std::int64_t x = cand_receiver; ; )
    {
        if (mark_[static_cast<std::size_t>(x)] == stamp_) { confluence = x; break; }
        gainers_.push_back(x);
        if (is_outlet_[static_cast<std::size_t>(x)]) { confluence = -1; break; }
        x = receiver_[static_cast<std::size_t>(x)];
    }

    // losers: old path from old_r to the confluence (exclusive), or the whole
    // old path to its outlet if the two paths are disjoint
    losers_.clear();
    for (std::int64_t x = old_r; ; x = receiver_[static_cast<std::size_t>(x)])
    {
        if (confluence >= 0 && x == confluence) break;
        losers_.push_back(x);
        if (is_outlet_[static_cast<std::size_t>(x)]) break;
    }

    double dE = 0.0;
    for (const std::int64_t nidx : losers_)
        dE += std::pow(area_[static_cast<std::size_t>(nidx)] - sub, g)
            - std::pow(area_[static_cast<std::size_t>(nidx)], g);
    for (const std::int64_t nidx : gainers_)
        dE += std::pow(area_[static_cast<std::size_t>(nidx)] + sub, g)
            - std::pow(area_[static_cast<std::size_t>(nidx)], g);

    return dE;
}

// -----------------------------------------------------------------------------
void lsdtt_synthetic_network_cpu::commit_(
    std::int64_t node, std::int64_t cand_receiver, double dE)
{
    const double sub = area_[static_cast<std::size_t>(node)];
    for (const std::int64_t nidx : losers_)
        area_[static_cast<std::size_t>(nidx)] -= sub;
    for (const std::int64_t nidx : gainers_)
        area_[static_cast<std::size_t>(nidx)] += sub;
    receiver_[static_cast<std::size_t>(node)] = cand_receiver;
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
        double T;
        if (it < n0) T = T0;
        else
        {
            const double k = static_cast<double>(it - n0);
            T = T0 * std::exp(-opts_.cooling_rate * k / static_cast<double>(n_));
        }

        const std::int64_t node = pick_node(rng_);
        if (!land_[static_cast<std::size_t>(node)]) continue;
        if (is_outlet_[static_cast<std::size_t>(node)]) continue;

        const int dir = pick_dir(rng_);
        const int nr = static_cast<int>(node / ncols_) + k_dr[dir];
        const int nc = static_cast<int>(node % ncols_) + k_dc[dir];
        if (nr < 0 || nr >= nrows_ || nc < 0 || nc >= ncols_) continue;

        const std::int64_t cand = static_cast<std::int64_t>(nr) * ncols_ + nc;
        if (!land_[static_cast<std::size_t>(cand)]) continue;
        if (cand == receiver_[static_cast<std::size_t>(node)]) continue;
        if (creates_cycle_(node, cand)) continue;

        const double dE = delta_energy_(node, cand);

        bool accept = (dE <= 0.0);
        if (!accept && T > 0.0) accept = (unit(rng_) < std::exp(-dE / T));
        if (accept) commit_(node, cand, dE);
    }

    recompute_all_areas_();
    energy_ = total_energy_();
}

// -----------------------------------------------------------------------------
lsdtt_raster_result lsdtt_synthetic_network_cpu::get_drainage_area_cells_raster() const
{
    lsdtt_raster_result out;
    out.meta = build_meta_();
    out.data.resize(static_cast<std::size_t>(n_));
    out.mask.assign(static_cast<std::size_t>(n_), 0);
    for (std::int64_t f = 0; f < n_; ++f)
    {
        if (land_[static_cast<std::size_t>(f)])
        {
            out.data[static_cast<std::size_t>(f)] =
                static_cast<float>(area_[static_cast<std::size_t>(f)]);
            out.mask[static_cast<std::size_t>(f)] = 1;
        }
        else
        {
            out.data[static_cast<std::size_t>(f)] =
                static_cast<float>(out.meta.no_data_value);
        }
    }
    return out;
}

// -----------------------------------------------------------------------------
lsdtt_raster_result lsdtt_synthetic_network_cpu::get_elevation_raster() const
{
    const double g    = opts_.gamma;
    const double invN = 1.0 / static_cast<double>(n_);
    const lsdtt_raster_meta meta = build_meta_();

    std::vector<float> z(static_cast<std::size_t>(n_),
                         static_cast<float>(meta.no_data_value));

    for (auto it = topo_order_.rbegin(); it != topo_order_.rend(); ++it)
    {
        const std::int64_t node = *it;
        if (is_outlet_[static_cast<std::size_t>(node)])
        {
            z[static_cast<std::size_t>(node)] = static_cast<float>(opts_.z_min);
            continue;
        }
        const double slope = opts_.slope0 *
            std::pow(area_[static_cast<std::size_t>(node)] * invN, g - 1.0);
        const double dz = slope * receiver_length_(node);
        z[static_cast<std::size_t>(node)] =
            z[static_cast<std::size_t>(receiver_[static_cast<std::size_t>(node)])]
            + static_cast<float>(dz);
    }

    lsdtt_raster_result out;
    out.meta = meta;
    out.data = std::move(z);
    out.mask.assign(static_cast<std::size_t>(n_), 0);
    for (std::int64_t f = 0; f < n_; ++f)
        if (land_[static_cast<std::size_t>(f)]) out.mask[static_cast<std::size_t>(f)] = 1;
    return out;
}

} // namespace lsdtt
