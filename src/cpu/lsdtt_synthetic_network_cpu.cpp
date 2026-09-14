// lsdtt_synthetic_network_cpu.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
//
// See lsdtt_synthetic_network_cpu.hpp for the algorithm description.
// -----------------------------------------------------------------------------

#include "lsdtt_synthetic_network_cpu.hpp"

#include <algorithm>
#include <cmath>
#include <queue>
#include <limits>
#include <utility>
#include <functional>
#include <stdexcept>

#include "lsdtt_steady_state_cpu.hpp"

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
    if (opts_.use_geographic_anchor)
    {
        meta.set_synthetic_georeference(opts_.centre_latitude_deg,
                                        opts_.centre_longitude_deg,
                                        opts_.projection_epsg);
    }
    else
    {
        meta.set_default_synthetic_georeference();
    }
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
        if (r < 0 || r >= nrows_) return false;
        c = wrap_col(c);
        if (c < 0) return false;
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
    else if (opts_.outlet_mode == ocn_outlet_mode::north_south)
    {
        // Land pixels on the north and south edges, every outlet_spacing-th
        // column, and nothing else. With periodic_columns this gives a strip
        // that drains through two opposed boundaries, so the main divide is
        // free to sit anywhere between them.
        //
        // outlet_spacing is what sets the BASIN WIDTH, and it is worth being
        // explicit about why it exists. At the default of 1 every boundary
        // pixel is its own baselevel, so the tile is partitioned into as many
        // basins as there are boundary cells and the basin scale is imposed by
        // the grid rather than emerging from the optimisation. Annealing
        // temperature cannot change that: it controls how each basin BRANCHES,
        // not how many there are. Raising outlet_spacing forces the drainage to
        // organise into wider basins, which is the OCN analogue of the
        // diamond-square feature wavelength, and the continuous version of the
        // single-outlet configuration used in much of the OCN literature.
        const int spacing = (opts_.outlet_spacing > 1) ? opts_.outlet_spacing : 1;

        for (int c = 0; c < ncols_; c += spacing)
        {
            for (int r : {0, nrows_ - 1})
            {
                if (!land_at(r, c)) continue;
                const std::int64_t f =
                    static_cast<std::int64_t>(r) * ncols_ + c;
                if (is_outlet_[static_cast<std::size_t>(f)]) continue;  // nrows_ == 1
                is_outlet_[static_cast<std::size_t>(f)] = 1;
                outlets_.push_back(f);
            }
        }

        // A spaced set can miss every land pixel on a ragged coast, which would
        // leave the domain with no exit at all. Fall back to every column
        // rather than throw: the caller asked for a spacing, not for failure.
        if (outlets_.empty() && spacing > 1)
        {
            for (int c = 0; c < ncols_; ++c)
            {
                for (int r : {0, nrows_ - 1})
                {
                    if (!land_at(r, c)) continue;
                    const std::int64_t f =
                        static_cast<std::int64_t>(r) * ncols_ + c;
                    if (is_outlet_[static_cast<std::size_t>(f)]) continue;
                    is_outlet_[static_cast<std::size_t>(f)] = 1;
                    outlets_.push_back(f);
                }
            }
        }

        if (outlets_.empty())
            throw std::runtime_error(
                "lsdtt_synthetic_network_cpu: outlet_mode 'north_south' found no "
                "land pixels on the north or south edges");
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

    if (opts_.energy_length_weighted)
    {
        // Dijkstra on EUCLIDEAN link length.
        //
        // An unweighted BFS over eight neighbours is the CHEBYSHEV metric, in
        // which a diagonal step and a cardinal step are the same distance. That
        // is the wrong metric for a length-weighted energy: it seeds a network
        // whose shortest paths are already diagonal, and single-flip annealing
        // at low temperature is greedy enough not to undo it. Seeding on the
        // metric the energy actually uses removes that head start.
        //
        // Ties are broken by flat index so the seed is deterministic and
        // independent of the queue's internal ordering.
        const double kInf = std::numeric_limits<double>::infinity();
        std::vector<double> dist(static_cast<std::size_t>(n_), kInf);

        using Entry = std::pair<double, std::int64_t>;   // (distance, node)
        std::priority_queue<Entry, std::vector<Entry>, std::greater<Entry>> pq;

        for (const std::int64_t o : outlets_)
        {
            dist[static_cast<std::size_t>(o)] = 0.0;
            pq.emplace(0.0, o);
        }

        while (!pq.empty())
        {
            const auto [d, f] = pq.top();
            pq.pop();
            if (d > dist[static_cast<std::size_t>(f)]) continue;   // stale
            visited[static_cast<std::size_t>(f)] = 1;

            const int r = static_cast<int>(f / ncols_);
            const int c = static_cast<int>(f % ncols_);
            for (int k = 0; k < 8; ++k)
            {
                const int nr = r + k_dr[k];
                const int nc = wrap_col(c + k_dc[k]);
                if (nr < 0 || nr >= nrows_ || nc < 0) continue;
                const std::int64_t nf = static_cast<std::int64_t>(nr) * ncols_ + nc;
                if (!land_[static_cast<std::size_t>(nf)]) continue;

                const double nd = d + link_length_(nf, f);
                const std::size_t ni = static_cast<std::size_t>(nf);
                if (nd < dist[ni] - 1e-12)
                {
                    dist[ni] = nd;
                    receiver_[ni] = f;
                    pq.emplace(nd, nf);
                }
            }
        }

        for (std::int64_t f = 0; f < n_; ++f)
            if (land_[static_cast<std::size_t>(f)] &&
                dist[static_cast<std::size_t>(f)] < kInf)
                visited[static_cast<std::size_t>(f)] = 1;
    }
    else
    {
        while (!q.empty())
        {
            const std::int64_t f = q.front();
            q.pop();
            const int r = static_cast<int>(f / ncols_);
            const int c = static_cast<int>(f % ncols_);
            for (int k = 0; k < 8; ++k)
            {
                const int nr = r + k_dr[k];
                const int nc = wrap_col(c + k_dc[k]);
                if (nr < 0 || nr >= nrows_ || nc < 0) continue;
                const std::int64_t nf = static_cast<std::int64_t>(nr) * ncols_ + nc;
                if (!land_[static_cast<std::size_t>(nf)] ||
                    visited[static_cast<std::size_t>(nf)]) continue;
                visited[static_cast<std::size_t>(nf)] = 1;
                receiver_[static_cast<std::size_t>(nf)] = f;
                q.push(nf);
            }
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

    // E = sum A^gamma * L   (default, physically correct on an 8-connected
    // lattice) or E = sum A^gamma (OCNet parity). See energy_length_weighted.
    if (opts_.energy_length_weighted)
    {
        for (std::int64_t f = 0; f < n_; ++f)
            if (land_[static_cast<std::size_t>(f)])
                e += std::pow(area_[static_cast<std::size_t>(f)], g)
                   * receiver_length_(f);
    }
    else
    {
        for (std::int64_t f = 0; f < n_; ++f)
            if (land_[static_cast<std::size_t>(f)])
                e += std::pow(area_[static_cast<std::size_t>(f)], g);
    }
    return e;
}

// -----------------------------------------------------------------------------
double lsdtt_synthetic_network_cpu::receiver_length_(std::int64_t node) const
{
    return link_length_(node, receiver_[static_cast<std::size_t>(node)]);
}

// -----------------------------------------------------------------------------
// Length of the link from `node` to `recv`, taken from the indices alone so it
// can be evaluated for a CANDIDATE receiver that has not been committed.
//
// An outlet drains to itself (dr = dc = 0), which is not diagonal, so it is
// given the cardinal cell size. An outlet has no downstream link, so the value
// is a convention rather than a measurement -- it only matters because the
// energy sum runs over every land node, and it is the same convention on both
// sides of any comparison.
double lsdtt_synthetic_network_cpu::link_length_(std::int64_t node,
                                                 std::int64_t recv) const
{
    const int dr = std::abs(static_cast<int>(node / ncols_) -
                            static_cast<int>(recv / ncols_));
    int dc = std::abs(static_cast<int>(node % ncols_) -
                      static_cast<int>(recv % ncols_));

    // A step that wraps the periodic seam is one column apart, not ncols-1.
    // Without this a wrapped DIAGONAL step reads as dr+dc = 1+(ncols-1), fails
    // the diagonal test, and is given the cardinal length -- an under-estimate
    // of sqrt(2) on exactly the steps that cross the seam.
    if (opts_.periodic_columns && dc > ncols_ - dc) dc = ncols_ - dc;

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

    if (!opts_.energy_length_weighted)
    {
        for (const std::int64_t nidx : losers_)
            dE += std::pow(area_[static_cast<std::size_t>(nidx)] - sub, g)
                - std::pow(area_[static_cast<std::size_t>(nidx)], g);
        for (const std::int64_t nidx : gainers_)
            dE += std::pow(area_[static_cast<std::size_t>(nidx)] + sub, g)
                - std::pow(area_[static_cast<std::size_t>(nidx)], g);
        return dE;
    }

    // Length-weighted: E = sum A^gamma * L.
    //
    // Losers and gainers keep their own receivers, so their L is unchanged and
    // only their area moves -- the area delta is simply scaled by their L.
    for (const std::int64_t nidx : losers_)
        dE += (std::pow(area_[static_cast<std::size_t>(nidx)] - sub, g)
             - std::pow(area_[static_cast<std::size_t>(nidx)], g))
            * receiver_length_(nidx);
    for (const std::int64_t nidx : gainers_)
        dE += (std::pow(area_[static_cast<std::size_t>(nidx)] + sub, g)
             - std::pow(area_[static_cast<std::size_t>(nidx)], g))
            * receiver_length_(nidx);

    // The moved node itself: its AREA is unchanged but its LINK LENGTH changes,
    // because that is exactly what the flip does. This single term is the whole
    // diagonal penalty -- without it a diagonal link costs the same as a
    // cardinal one and lateral movement is free.
    const double L_old = link_length_(node, old_r);
    const double L_new = link_length_(node, cand_receiver);
    dE += std::pow(sub, g) * (L_new - L_old);

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
        const int nc = wrap_col(static_cast<int>(node % ncols_) + k_dc[dir]);
        if (nr < 0 || nr >= nrows_ || nc < 0) continue;

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
// -----------------------------------------------------------------------------
// Elevation: chi-based fluvial steady state on the OCN tree, computed by the
// SAME shared core that lsdtt-lem's steady-state solver uses, so an OCN and a
// LEM steady state on the same graph and parameters give the same surface.
// Drainage area is converted to map area (cells * cell_size^2) to match the D8
// flow router's accumulate() convention. Fluvial relief only (no hillslope).
// -----------------------------------------------------------------------------
lsdtt_raster_result lsdtt_synthetic_network_cpu::get_elevation_raster() const
{
    const lsdtt_raster_meta meta = build_meta_();
    const double cs  = opts_.cell_size;
    const double cs2 = cs * cs;

    std::vector<double>       area_m2(static_cast<std::size_t>(n_), 0.0);
    std::vector<float>        dx_eff(static_cast<std::size_t>(n_),
                                     static_cast<float>(cs));
    std::vector<std::uint8_t> valid(static_cast<std::size_t>(n_), 0);
    for (std::int64_t f = 0; f < n_; ++f)
    {
        const std::size_t s = static_cast<std::size_t>(f);
        if (!land_[s]) continue;
        valid[s]   = 1;
        area_m2[s] = area_[s] * cs2;
        if (receiver_[s] >= 0 && receiver_[s] != f)
            dx_eff[s] = static_cast<float>(receiver_length_(f));
    }

    lsdtt_fluvial_steady_state_options sopts;
    sopts.m                   = opts_.m_exponent;
    sopts.n                   = opts_.n_exponent;
    sopts.U                   = opts_.uplift_rate;
    sopts.K                   = opts_.erodibility;
    sopts.relief_mode         = opts_.relief_mode;
    sopts.target_relief       = opts_.target_relief;
    sopts.baselevel_elevation = opts_.baselevel_elevation;
    sopts.divide_area_by_dx   = true;

    const lsdtt_fluvial_steady_state_field field =
        compute_fluvial_steady_state_on_tree(
            receiver_, area_m2, dx_eff, valid, topo_order_,
            cs, meta.no_data_value, sopts);

    lsdtt_raster_result out;
    out.meta = meta;
    out.data = field.elevation;
    out.mask = std::move(valid);
    return out;
}

} // namespace lsdtt
