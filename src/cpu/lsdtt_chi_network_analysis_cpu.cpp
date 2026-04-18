// lsdtt_chi_network_analysis_cpu.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html
//
// -----------------------------------------------------------------------------
// lsdtt_chi_network_analysis_cpu
//
// Analysis tools for extracted chi channel networks.
//
// First tool implemented here:
//   Recursive piecewise-linear fitting of an extracted chi network in
//   chi-elevation space.
//
// The fitted network preserves the same nodes as the input chi network.
// Output includes fitted elevation at every node and the locations of all
// interior split points.
// -----------------------------------------------------------------------------

#include "lsdtt_chi_network_analysis_cpu.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace lsdtt
{
namespace
{

constexpr double k_sqrt_2 = 1.4142135623730950488;

// -----------------------------------------------------------------------------
// small helpers
// -----------------------------------------------------------------------------
bool nearly_equal(double a, double b)
{
    const double scale =
        std::max(1.0, std::max(std::fabs(a), std::fabs(b)));

    return std::fabs(a - b) <= 1e-9 * scale;
}

bool fitted_values_consistent(float a, float b)
{
    const double scale =
        std::max(1.0, std::max(std::fabs(static_cast<double>(a)),
                               std::fabs(static_cast<double>(b))));

    return std::fabs(static_cast<double>(a) - static_cast<double>(b))
           <= 1e-6 * scale;
}

void validate_fit_options(const lsdtt_chi_network_fit_options& options)
{
    if (!std::isfinite(options.critical_divergence_m) ||
        options.critical_divergence_m < 0.0f)
    {
        throw std::invalid_argument(
            "lsdtt_chi_network_analysis_cpu: critical_divergence_m must be "
            "finite and >= 0."
        );
    }
}

double edge_length_from_flats(
    std::int64_t upstream_flat,
    std::int64_t downstream_flat,
    int n_cols,
    double cell_size)
{
    const int up_row = static_cast<int>(upstream_flat / static_cast<std::int64_t>(n_cols));
    const int up_col = static_cast<int>(upstream_flat % static_cast<std::int64_t>(n_cols));

    const int dn_row = static_cast<int>(downstream_flat / static_cast<std::int64_t>(n_cols));
    const int dn_col = static_cast<int>(downstream_flat % static_cast<std::int64_t>(n_cols));

    const int dr = std::abs(up_row - dn_row);
    const int dc = std::abs(up_col - dn_col);

    if (dr == 0 && dc == 0)
    {
        return 0.0;
    }

    if (dr > 1 || dc > 1)
    {
        throw std::runtime_error(
            "lsdtt_chi_network_analysis_cpu: chi network receiver relation is "
            "not D8-adjacent."
        );
    }

    return (dr == 1 && dc == 1) ? cell_size * k_sqrt_2 : cell_size;
}

// =============================================================================
// chi_network_fitter
// =============================================================================
class chi_network_fitter
{
public:
    chi_network_fitter(
        const lsdtt_chi_basin_result&         chi_network,
        const lsdtt_raster_cpu&               dem,
        const lsdtt_chi_network_fit_options&  options
    )
        : chi_network_(chi_network),
          dem_(dem),
          options_(options)
    {
    }

    lsdtt_chi_network_fit_result run()
    {
        validate_fit_options(options_);
        validate_inputs();
        initialise_result();
        sample_observed_elevation();
        build_local_graph();
        compute_downstream_distances();
        compute_upstream_metrics();
        fit_entire_network();
        finalise_split_points();
        finalise_segment_parents();
        return result_;
    }

private:
    // -------------------------------------------------------------------------
    // input
    // -------------------------------------------------------------------------
    const lsdtt_chi_basin_result&        chi_network_;
    const lsdtt_raster_cpu&              dem_;
    const lsdtt_chi_network_fit_options  options_;

    // -------------------------------------------------------------------------
    // local graph / bookkeeping
    // -------------------------------------------------------------------------
    std::unordered_map<std::int64_t, int> local_index_of_flat_;
    std::vector<int>                      receiver_local_;
    std::vector<std::vector<int>>         donors_local_;
    std::vector<double>                   step_to_receiver_;
    std::vector<double>                   distance_to_outlet_;

    // for max flow-path mainstem choice
    std::vector<double>                   max_upstream_length_;
    std::vector<int>                      best_donor_by_length_;

    // for max chi mainstem choice
    std::vector<float>                    max_upstream_source_chi_;
    std::vector<int>                      best_donor_by_chi_;

    // fitted bookkeeping
    std::vector<std::uint8_t>             fitted_assigned_;
    int                                   root_local_index_ = -1;
    int                                   next_branch_id_   = 1;
    int                                   next_segment_id_  = 0;

    // -------------------------------------------------------------------------
    // output
    // -------------------------------------------------------------------------
    lsdtt_chi_network_fit_result result_;

private:
    void validate_inputs() const
    {
        dem_.get_meta().check_projected_metres();

        const std::size_t n = chi_network_.flats.size();

        if (n == 0)
        {
            throw std::invalid_argument(
                "lsdtt_chi_network_analysis_cpu: chi network is empty."
            );
        }

        if (chi_network_.chi.size() != n ||
            chi_network_.receiver_flats.size() != n)
        {
            throw std::invalid_argument(
                "lsdtt_chi_network_analysis_cpu: chi network vectors have "
                "inconsistent lengths."
            );
        }

        if (chi_network_.outlet_flat == k_nodata_receiver)
        {
            throw std::invalid_argument(
                "lsdtt_chi_network_analysis_cpu: chi network has invalid outlet_flat."
            );
        }

        if (chi_network_.flats[0] != chi_network_.outlet_flat)
        {
            throw std::invalid_argument(
                "lsdtt_chi_network_analysis_cpu: chi network must begin with "
                "the outlet node."
            );
        }

        if (chi_network_.receiver_flats[0] != k_outlet_receiver)
        {
            throw std::invalid_argument(
                "lsdtt_chi_network_analysis_cpu: outlet node must have "
                "receiver_flats[0] == k_outlet_receiver."
            );
        }

        if (!std::isfinite(chi_network_.outlet_chi))
        {
            throw std::invalid_argument(
                "lsdtt_chi_network_analysis_cpu: outlet_chi must be finite."
            );
        }

        if (!std::isfinite(chi_network_.chi[0]))
        {
            throw std::invalid_argument(
                "lsdtt_chi_network_analysis_cpu: outlet chi coordinate must be finite."
            );
        }
    }

    void initialise_result()
    {
        result_.outlet_flat        = chi_network_.outlet_flat;
        result_.flats             = chi_network_.flats;
        result_.receiver_flats    = chi_network_.receiver_flats;
        result_.chi               = chi_network_.chi;
        result_.observed_elevation.assign(chi_network_.flats.size(), 0.0f);
        result_.fitted_elevation.assign(
            chi_network_.flats.size(),
            std::numeric_limits<float>::quiet_NaN()
        );
        result_.is_split_point.assign(chi_network_.flats.size(), std::uint8_t{0});
        result_.split_flats.clear();
        result_.segments.clear();

        fitted_assigned_.assign(chi_network_.flats.size(), std::uint8_t{0});
    }

    void sample_observed_elevation()
    {
        const std::vector<float>&   dem_data = dem_.get_data();
        const std::vector<uint8_t>& dem_mask = dem_.get_mask();

        for (std::size_t i = 0; i < chi_network_.flats.size(); ++i)
        {
            const std::int64_t flat = chi_network_.flats[i];

            if (flat < 0 ||
                flat >= static_cast<std::int64_t>(dem_data.size()))
            {
                throw std::runtime_error(
                    "lsdtt_chi_network_analysis_cpu: chi network flat index "
                    "lies outside the DEM."
                );
            }

            const std::size_t idx = static_cast<std::size_t>(flat);

            if (dem_mask[idx] == 0)
            {
                throw std::runtime_error(
                    "lsdtt_chi_network_analysis_cpu: DEM is nodata at a chi "
                    "network node."
                );
            }

            const float z = dem_data[idx];
            if (!std::isfinite(z))
            {
                throw std::runtime_error(
                    "lsdtt_chi_network_analysis_cpu: DEM elevation is not "
                    "finite at a chi network node."
                );
            }

            result_.observed_elevation[i] = z;
        }
    }

    void build_local_graph()
    {
        const std::size_t n = chi_network_.flats.size();

        local_index_of_flat_.clear();
        local_index_of_flat_.reserve(n);

        for (std::size_t i = 0; i < n; ++i)
        {
            const auto inserted =
                local_index_of_flat_.emplace(chi_network_.flats[i], static_cast<int>(i));

            if (!inserted.second)
            {
                throw std::invalid_argument(
                    "lsdtt_chi_network_analysis_cpu: chi network contains "
                    "duplicate flat indices."
                );
            }
        }

        root_local_index_ = local_index_of_flat_.at(chi_network_.outlet_flat);

        receiver_local_.assign(n, -1);
        donors_local_.assign(n, {});
        step_to_receiver_.assign(n, 0.0);

        const int    n_cols    = dem_.get_meta().n_cols;
        const double cell_size = dem_.get_meta().cell_size;

        for (std::size_t i = 0; i < n; ++i)
        {
            const std::int64_t flat = chi_network_.flats[i];
            const std::int64_t rec  = chi_network_.receiver_flats[i];

            const float chi_here = chi_network_.chi[i];
            if (!std::isfinite(chi_here))
            {
                throw std::invalid_argument(
                    "lsdtt_chi_network_analysis_cpu: chi network contains "
                    "non-finite chi values."
                );
            }

            if (i == static_cast<std::size_t>(root_local_index_))
            {
                if (rec != k_outlet_receiver)
                {
                    throw std::invalid_argument(
                        "lsdtt_chi_network_analysis_cpu: outlet node does not "
                        "route to k_outlet_receiver."
                    );
                }

                receiver_local_[i]   = -1;
                step_to_receiver_[i] = 0.0;
                continue;
            }

            if (rec == k_outlet_receiver)
            {
                throw std::invalid_argument(
                    "lsdtt_chi_network_analysis_cpu: non-root node routes to "
                    "k_outlet_receiver."
                );
            }

            const auto found = local_index_of_flat_.find(rec);
            if (found == local_index_of_flat_.end())
            {
                throw std::invalid_argument(
                    "lsdtt_chi_network_analysis_cpu: receiver flat is not "
                    "present in the chi network."
                );
            }

            const int rec_local = found->second;

            if (rec_local >= static_cast<int>(i))
            {
                throw std::invalid_argument(
                    "lsdtt_chi_network_analysis_cpu: chi network is not in "
                    "receiver-before-donor order."
                );
            }

            if (!(chi_here >= chi_network_.chi[static_cast<std::size_t>(rec_local)]))
            {
                throw std::invalid_argument(
                    "lsdtt_chi_network_analysis_cpu: chi values decrease from "
                    "receiver to donor."
                );
            }

            receiver_local_[i] = rec_local;
            donors_local_[static_cast<std::size_t>(rec_local)].push_back(static_cast<int>(i));
            step_to_receiver_[i] = edge_length_from_flats(
                flat,
                rec,
                n_cols,
                cell_size
            );
        }
    }

    void compute_downstream_distances()
    {
        distance_to_outlet_.assign(chi_network_.flats.size(), 0.0);

        for (std::size_t i = 0; i < chi_network_.flats.size(); ++i)
        {
            const int rec = receiver_local_[i];

            if (rec < 0)
            {
                distance_to_outlet_[i] = 0.0;
            }
            else
            {
                distance_to_outlet_[i] =
                    distance_to_outlet_[static_cast<std::size_t>(rec)] +
                    step_to_receiver_[i];
            }
        }
    }

    void compute_upstream_metrics()
    {
        const std::size_t n = chi_network_.flats.size();

        max_upstream_length_.assign(n, 0.0);
        best_donor_by_length_.assign(n, -1);

        max_upstream_source_chi_.assign(
            n,
            -std::numeric_limits<float>::infinity()
        );
        best_donor_by_chi_.assign(n, -1);

        for (std::int64_t i = static_cast<std::int64_t>(n) - 1; i >= 0; --i)
        {
            const std::size_t idx = static_cast<std::size_t>(i);

            if (donors_local_[idx].empty())
            {
                max_upstream_length_[idx]     = 0.0;
                best_donor_by_length_[idx]    = -1;
                max_upstream_source_chi_[idx] = chi_network_.chi[idx];
                best_donor_by_chi_[idx]       = -1;
                continue;
            }

            // ---- max flow-path length ----
            double best_length = -1.0;
            int    best_donor_length = -1;

            for (const int donor : donors_local_[idx])
            {
                const double candidate =
                    step_to_receiver_[static_cast<std::size_t>(donor)] +
                    max_upstream_length_[static_cast<std::size_t>(donor)];

                if (candidate > best_length)
                {
                    best_length       = candidate;
                    best_donor_length = donor;
                }
                else if (nearly_equal(candidate, best_length))
                {
                    // Tie-break on larger source chi.
                    if (max_upstream_source_chi_[static_cast<std::size_t>(donor)] >
                        max_upstream_source_chi_[static_cast<std::size_t>(best_donor_length)])
                    {
                        best_donor_length = donor;
                    }
                }
            }

            max_upstream_length_[idx]  = best_length;
            best_donor_by_length_[idx] = best_donor_length;

            // ---- max chi ----
            float best_source_chi = -std::numeric_limits<float>::infinity();
            int   best_donor_chi  = -1;

            for (const int donor : donors_local_[idx])
            {
                const float candidate =
                    max_upstream_source_chi_[static_cast<std::size_t>(donor)];

                if (candidate > best_source_chi)
                {
                    best_source_chi = candidate;
                    best_donor_chi  = donor;
                }
                else if (candidate == best_source_chi)
                {
                    // Tie-break on longer flow path.
                    if (max_upstream_length_[static_cast<std::size_t>(donor)] >
                        max_upstream_length_[static_cast<std::size_t>(best_donor_chi)])
                    {
                        best_donor_chi = donor;
                    }
                }
            }

            max_upstream_source_chi_[idx] = best_source_chi;
            best_donor_by_chi_[idx]       = best_donor_chi;
        }
    }

    int choose_primary_donor(int local_idx) const
    {
        if (options_.mainstem_choice == lsdtt_chi_mainstem_choice::max_chi)
        {
            return best_donor_by_chi_[static_cast<std::size_t>(local_idx)];
        }

        return best_donor_by_length_[static_cast<std::size_t>(local_idx)];
    }

    double tributary_priority(int donor_local_idx) const
    {
        if (options_.mainstem_choice == lsdtt_chi_mainstem_choice::max_chi)
        {
            return static_cast<double>(
                max_upstream_source_chi_[static_cast<std::size_t>(donor_local_idx)]
            );
        }

        return step_to_receiver_[static_cast<std::size_t>(donor_local_idx)] +
               max_upstream_length_[static_cast<std::size_t>(donor_local_idx)];
    }

    std::vector<int> build_branch_path(
        int anchor_local_idx,
        int first_upstream_local_idx
    ) const
    {
        std::vector<int> path;
        path.reserve(64);

        path.push_back(anchor_local_idx);

        int current = first_upstream_local_idx;
        if (current < 0)
        {
            current = choose_primary_donor(anchor_local_idx);
        }

        while (current >= 0)
        {
            path.push_back(current);
            current = choose_primary_donor(current);
        }

        return path;
    }

    void fit_entire_network()
    {
        // The outlet elevation anchors the whole fitted network.
        result_.fitted_elevation[static_cast<std::size_t>(root_local_index_)] =
            result_.observed_elevation[static_cast<std::size_t>(root_local_index_)];
        fitted_assigned_[static_cast<std::size_t>(root_local_index_)] = 1u;

        // Main stem always starts at the outlet and follows the default primary path.
        process_branch(root_local_index_, -1, 0);
    }

    void process_branch(
        int anchor_local_idx,
        int first_upstream_local_idx,
        int branch_id)
    {
        if (anchor_local_idx < 0 ||
            anchor_local_idx >= static_cast<int>(chi_network_.flats.size()))
        {
            throw std::runtime_error(
                "lsdtt_chi_network_analysis_cpu: invalid branch anchor index."
            );
        }

        if (!fitted_assigned_[static_cast<std::size_t>(anchor_local_idx)])
        {
            // Defensive fallback. In the normal workflow all tributary anchors
            // should already have been fitted from the parent branch.
            result_.fitted_elevation[static_cast<std::size_t>(anchor_local_idx)] =
                result_.observed_elevation[static_cast<std::size_t>(anchor_local_idx)];
            fitted_assigned_[static_cast<std::size_t>(anchor_local_idx)] = 1u;
        }

        const std::vector<int> path =
            build_branch_path(anchor_local_idx, first_upstream_local_idx);

        if (path.size() >= 2)
        {
            fit_path_recursive(
                path,
                0,
                static_cast<int>(path.size()) - 1,
                branch_id
            );
        }

        // Now process all side tributaries off this fitted path, walking
        // downstream to upstream so the parent path is already available.
        for (std::size_t pos = 0; pos < path.size(); ++pos)
        {
            const int node = path[pos];

            int main_continuation = -1;
            if (pos + 1 < path.size())
            {
                main_continuation = path[pos + 1];
            }

            std::vector<int> side_donors;
            side_donors.reserve(donors_local_[static_cast<std::size_t>(node)].size());

            for (const int donor : donors_local_[static_cast<std::size_t>(node)])
            {
                if (donor != main_continuation)
                {
                    side_donors.push_back(donor);
                }
            }

            std::sort(
                side_donors.begin(),
                side_donors.end(),
                [this](int a, int b)
                {
                    return tributary_priority(a) > tributary_priority(b);
                }
            );

            for (const int donor : side_donors)
            {
                process_branch(node, donor, next_branch_id_++);
            }
        }
    }

    void fit_path_recursive(
        const std::vector<int>& path,
        int down_pos,
        int up_pos,
        int branch_id)
    {
        if (down_pos < 0 || up_pos < 0 || down_pos >= up_pos)
        {
            return;
        }

        const int down_local = path[static_cast<std::size_t>(down_pos)];
        const int up_local   = path[static_cast<std::size_t>(up_pos)];

        if (!fitted_assigned_[static_cast<std::size_t>(down_local)])
        {
            throw std::runtime_error(
                "lsdtt_chi_network_analysis_cpu: downstream fit anchor was not assigned."
            );
        }

        const float z_down =
            result_.fitted_elevation[static_cast<std::size_t>(down_local)];

        const float z_up =
            fitted_assigned_[static_cast<std::size_t>(up_local)]
                ? result_.fitted_elevation[static_cast<std::size_t>(up_local)]
                : result_.observed_elevation[static_cast<std::size_t>(up_local)];

        const double chi_down =
            static_cast<double>(result_.chi[static_cast<std::size_t>(down_local)]);

        const double chi_up =
            static_cast<double>(result_.chi[static_cast<std::size_t>(up_local)]);

        if (!(chi_up >= chi_down))
        {
            throw std::runtime_error(
                "lsdtt_chi_network_analysis_cpu: branch chi is not monotonic."
            );
        }

        if (up_pos - down_pos <= 1)
        {
            accept_segment(path, down_pos, up_pos, branch_id, z_down, z_up);
            return;
        }

        int   split_pos   = -1;
        float max_abs_res = -1.0f;

        for (int pos = down_pos + 1; pos < up_pos; ++pos)
        {
            const int local_idx = path[static_cast<std::size_t>(pos)];

            const float z_line = line_elevation(
                chi_down,
                chi_up,
                z_down,
                z_up,
                static_cast<double>(result_.chi[static_cast<std::size_t>(local_idx)])
            );

            const float z_obs =
                result_.observed_elevation[static_cast<std::size_t>(local_idx)];

            const float abs_res = std::fabs(z_obs - z_line);

            if (abs_res > max_abs_res)
            {
                max_abs_res = abs_res;
                split_pos   = pos;
            }
        }

        if (split_pos > down_pos &&
            split_pos < up_pos &&
            max_abs_res > options_.critical_divergence_m)
        {
            const int split_local = path[static_cast<std::size_t>(split_pos)];

            // Split points are constrained to the observed profile.
            if (fitted_assigned_[static_cast<std::size_t>(split_local)])
            {
                const float z_existing =
                    result_.fitted_elevation[static_cast<std::size_t>(split_local)];
                const float z_obs =
                    result_.observed_elevation[static_cast<std::size_t>(split_local)];

                if (!fitted_values_consistent(z_existing, z_obs))
                {
                    throw std::runtime_error(
                        "lsdtt_chi_network_analysis_cpu: split point already "
                        "has an inconsistent fitted elevation."
                    );
                }
            }
            else
            {
                result_.fitted_elevation[static_cast<std::size_t>(split_local)] =
                    result_.observed_elevation[static_cast<std::size_t>(split_local)];
                fitted_assigned_[static_cast<std::size_t>(split_local)] = 1u;
            }

            result_.is_split_point[static_cast<std::size_t>(split_local)] = 1u;

            fit_path_recursive(path, down_pos, split_pos, branch_id);
            fit_path_recursive(path, split_pos, up_pos, branch_id);
            return;
        }

        accept_segment(path, down_pos, up_pos, branch_id, z_down, z_up);
    }

    float line_elevation(
        double chi_down,
        double chi_up,
        float z_down,
        float z_up,
        double chi_here) const
    {
        if (nearly_equal(chi_up, chi_down))
        {
            return z_down;
        }

        const double t = (chi_here - chi_down) / (chi_up - chi_down);

        return static_cast<float>(
            static_cast<double>(z_down) +
            t * (static_cast<double>(z_up) - static_cast<double>(z_down))
        );
    }

    void accept_segment(
        const std::vector<int>& path,
        int down_pos,
        int up_pos,
        int branch_id,
        float z_down,
        float z_up)
    {
        const int down_local = path[static_cast<std::size_t>(down_pos)];
        const int up_local   = path[static_cast<std::size_t>(up_pos)];

        const double chi_down =
            static_cast<double>(result_.chi[static_cast<std::size_t>(down_local)]);
        const double chi_up =
            static_cast<double>(result_.chi[static_cast<std::size_t>(up_local)]);

        for (int pos = down_pos; pos <= up_pos; ++pos)
        {
            const int local_idx = path[static_cast<std::size_t>(pos)];

            const float z_fit = line_elevation(
                chi_down,
                chi_up,
                z_down,
                z_up,
                static_cast<double>(result_.chi[static_cast<std::size_t>(local_idx)])
            );

            if (fitted_assigned_[static_cast<std::size_t>(local_idx)])
            {
                const float z_existing =
                    result_.fitted_elevation[static_cast<std::size_t>(local_idx)];

                if (!fitted_values_consistent(z_existing, z_fit))
                {
                    throw std::runtime_error(
                        "lsdtt_chi_network_analysis_cpu: inconsistent fitted "
                        "elevation assigned to a chi network node."
                    );
                }
            }
            else
            {
                result_.fitted_elevation[static_cast<std::size_t>(local_idx)] = z_fit;
                fitted_assigned_[static_cast<std::size_t>(local_idx)] = 1u;
            }
        }

        lsdtt_chi_network_fit_segment segment;
        segment.segment_id        = next_segment_id_++;
        segment.parent_segment_id = -1;
        segment.branch_id         = branch_id;

        segment.downstream_flat = result_.flats[static_cast<std::size_t>(down_local)];
        segment.upstream_flat   = result_.flats[static_cast<std::size_t>(up_local)];

        segment.downstream_chi =
            result_.chi[static_cast<std::size_t>(down_local)];
        segment.upstream_chi =
            result_.chi[static_cast<std::size_t>(up_local)];

        segment.downstream_elevation = z_down;
        segment.upstream_elevation   = z_up;

        result_.segments.push_back(segment);
    }

    void finalise_split_points()
    {
        result_.split_flats.clear();

        for (std::size_t i = 0; i < result_.is_split_point.size(); ++i)
        {
            if (result_.is_split_point[i] != 0)
            {
                result_.split_flats.push_back(result_.flats[i]);
            }
        }
    }

    void finalise_segment_parents()
    {
        for (auto& seg : result_.segments)
        {
            seg.parent_segment_id = -1;
        }

        // For each branch, the parent segment of an upstream segment is the
        // downstream segment on the same branch whose upstream endpoint equals
        // the current segment's downstream endpoint.
        for (auto& seg : result_.segments)
        {
            for (const auto& candidate_parent : result_.segments)
            {
                if (candidate_parent.segment_id == seg.segment_id)
                {
                    continue;
                }

                if (candidate_parent.branch_id != seg.branch_id)
                {
                    continue;
                }

                if (candidate_parent.upstream_flat == seg.downstream_flat)
                {
                    seg.parent_segment_id = candidate_parent.segment_id;
                    break;
                }
            }
        }
    }
};

} // anonymous namespace

// =============================================================================
// public API
// =============================================================================
lsdtt_chi_network_fit_result fit_chi_network(
    const lsdtt_chi_basin_result& chi_network,
    const lsdtt_raster_cpu& dem,
    const lsdtt_chi_network_fit_options& options)
{
    chi_network_fitter fitter(chi_network, dem, options);
    return fitter.run();
}

} // namespace lsdtt
