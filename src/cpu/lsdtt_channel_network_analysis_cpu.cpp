// lsdtt_channel_network_analysis_cpu.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html

#include "lsdtt_channel_network_analysis_cpu.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <vector>
#include <limits>
#include <sstream>

#include "../utils/lsdtt_raster_geometry.hpp"
#include "../utils/lsdtt_polygon_geometry.hpp"
#include "../utils/lsdtt_regression.hpp"
#include "../core/lsdtt_data_loader.hpp"
#include "../core/lsdtt_vector_utils.hpp"


namespace lsdtt
{
namespace
{


double quiet_nan()
{
    return std::numeric_limits<double>::quiet_NaN();
}


double mean_total_link_length_from_bearing_records(
    const std::vector<lsdtt_channel_link_bearing>& records
)
{
    double sum = 0.0;
    std::size_t n = 0;

    for (const auto& record : records)
    {
        if (record.empty())
        {
            continue;
        }

        if (std::isfinite(record.total_length_metres) &&
            record.total_length_metres > 0.0)
        {
            sum += record.total_length_metres;
            ++n;
        }
    }

    if (n == 0)
    {
        return quiet_nan();
    }

    return sum / static_cast<double>(n);
}

double mean_total_link_length_from_junction_angle_records(
    const std::vector<lsdtt_junction_angle_record>& records
)
{
    double sum = 0.0;
    std::size_t n = 0;

    for (const auto& record : records)
    {
        if (record.empty())
        {
            continue;
        }

        const double lengths[3] = {
            record.donor_1_total_length_metres,
            record.donor_2_total_length_metres,
            record.receiver_total_length_metres
        };

        for (double length : lengths)
        {
            if (std::isfinite(length) && length > 0.0)
            {
                sum += length;
                ++n;
            }
        }
    }

    if (n == 0)
    {
        return quiet_nan();
    }

    return sum / static_cast<double>(n);
}

std::vector<lsdtt_xy> downstream_vector_geometry(
    double origin_x,
    double origin_y,
    double direction_x,
    double direction_y,
    double vector_length_metres
)
{
    std::vector<lsdtt_xy> points(2);

    points[0] = lsdtt_xy{origin_x, origin_y};
    points[1] = lsdtt_xy{
        origin_x + vector_length_metres * direction_x,
        origin_y + vector_length_metres * direction_y
    };

    return points;
}

std::vector<lsdtt_xy> donor_arm_geometry_from_inward_bearing(
    double junction_x,
    double junction_y,
    double inward_bearing_rad,
    double arm_length_metres
)
{
    const double dx = std::sin(inward_bearing_rad);
    const double dy = std::cos(inward_bearing_rad);

    std::vector<lsdtt_xy> points(2);

    // Donor bearing points toward the junction, so the arm should terminate
    // at the junction and extend only upstream of it.
    points[0] = lsdtt_xy{
        junction_x - arm_length_metres * dx,
        junction_y - arm_length_metres * dy
    };

    points[1] = lsdtt_xy{
        junction_x,
        junction_y
    };

    return points;
}

std::vector<lsdtt_xy> receiver_arm_geometry_from_outward_bearing(
    double junction_x,
    double junction_y,
    double outward_bearing_rad,
    double arm_length_metres
)
{
    const double dx = std::sin(outward_bearing_rad);
    const double dy = std::cos(outward_bearing_rad);

    std::vector<lsdtt_xy> points(2);

    // Receiver bearing points downstream away from the junction, so the arm
    // should start at the junction and extend only downstream of it.
    points[0] = lsdtt_xy{
        junction_x,
        junction_y
    };

    points[1] = lsdtt_xy{
        junction_x + arm_length_metres * dx,
        junction_y + arm_length_metres * dy
    };

    return points;
}




std::string precise_string(double value)
{
    std::ostringstream ss;
    ss.precision(17);
    ss << value;
    return ss.str();
}

lsdtt_lonlat flat_to_lonlat(
    const lsdtt_flow_base&              flow,
    std::int64_t                        flat,
    const lsdtt_coordinate_transform&   transform
)
{
    const int row = flow.row_from_flat_unchecked(flat);
    const int col = flow.col_from_flat_unchecked(flat);

    return cell_center_to_lonlat(flow.meta(), row, col, transform);
}

lsdtt_xy flat_to_xy(
    const lsdtt_flow_base& flow,
    std::int64_t           flat
)
{
    const int row = flow.row_from_flat_unchecked(flat);
    const int col = flow.col_from_flat_unchecked(flat);

    return cell_center_from_row_col(flow.meta(), row, col);
}



double wrap_angle_0_2pi(double angle)
{
    const double two_pi = 2.0 * std::acos(-1.0);

    angle = std::fmod(angle, two_pi);

    if (angle < 0.0)
    {
        angle += two_pi;
    }

    return angle;
}

double wrap_angle_minus_pi_pi(double angle)
{
    const double pi = std::acos(-1.0);
    const double two_pi = 2.0 * pi;

    angle = std::fmod(angle + pi, two_pi);

    if (angle < 0.0)
    {
        angle += two_pi;
    }

    return angle - pi;
}

double radians_to_degrees(double angle_rad)
{
    return angle_rad * 180.0 / std::acos(-1.0);
}

double drainage_area_from_flat(
    const lsdtt_flow_base& flow,
    std::int64_t           flat
)
{
    return static_cast<double>(flow.s_vector_size(flat)) *
           flow.cell_size() *
           flow.cell_size();
}

std::vector<std::int32_t> resolve_junctions_to_analyse(
    const lsdtt_junction_network_base& network,
    const std::vector<std::int32_t>&   requested_junctions
)
{
    std::vector<std::int32_t> junctions;

    if (requested_junctions.empty())
    {
        junctions.reserve(static_cast<std::size_t>(network.n_junctions()));

        for (std::int32_t j = 0; j < network.n_junctions(); ++j)
        {
            junctions.push_back(j);
        }

        return junctions;
    }

    junctions = requested_junctions;

    for (const std::int32_t junction : junctions)
    {
        if (!network.is_valid_junction(junction))
        {
            throw std::runtime_error(
                "calculate_junction_angles: requested junction is invalid"
            );
        }
    }

    std::sort(junctions.begin(), junctions.end());
    junctions.erase(
        std::unique(junctions.begin(), junctions.end()),
        junctions.end()
    );

    return junctions;
}



double bearing_clockwise_from_north(
    double dx,
    double dy
)
{
    if (dx == 0.0 && dy == 0.0)
    {
        return std::numeric_limits<double>::quiet_NaN();
    }

    double bearing = std::atan2(dx, dy);

    if (bearing < 0.0)
    {
        bearing += 2.0 * std::acos(-1.0);
    }

    return bearing;
}

std::vector<std::int32_t> resolve_upstream_junctions_to_analyse(
    const lsdtt_junction_network_base& network,
    const std::vector<std::int32_t>&   requested_junctions
)
{
    std::vector<std::int32_t> junctions;

    if (requested_junctions.empty())
    {
        junctions.reserve(static_cast<std::size_t>(network.n_junctions()));

        for (std::int32_t j = 0; j < network.n_junctions(); ++j)
        {
            if (!network.is_outlet_junction(j))
            {
                junctions.push_back(j);
            }
        }

        return junctions;
    }

    junctions = requested_junctions;

    for (const std::int32_t junction : junctions)
    {
        if (!network.is_valid_junction(junction))
        {
            throw std::runtime_error(
                "calculate_channel_link_bearings: requested junction is invalid"
            );
        }
    }

    std::sort(junctions.begin(), junctions.end());
    junctions.erase(
        std::unique(junctions.begin(), junctions.end()),
        junctions.end()
    );

    return junctions;
}


std::vector<double> build_cumulative_distance(
    const std::vector<lsdtt_xy>& xy
)
{
    std::vector<double> cumulative(xy.size(), 0.0);

    for (std::size_t i = 1; i < xy.size(); ++i)
    {
        cumulative[i] =
            cumulative[i - 1] + planar_distance(xy[i - 1], xy[i]);
    }

    return cumulative;
}

void validate_sampling_options(
    const lsdtt_link_sampling_options& options
)
{
    if (options.min_nodes < 1)
    {
        throw std::invalid_argument(
            "sample_traced_channel_link: min_nodes must be >= 1"
        );
    }

    if (options.mode == lsdtt_link_sampling_mode::fixed_n_pixels &&
        options.fixed_n_pixels <= 0)
    {
        throw std::invalid_argument(
            "sample_traced_channel_link: fixed_n_pixels mode requires "
            "fixed_n_pixels > 0"
        );
    }

    if (options.mode == lsdtt_link_sampling_mode::fixed_length &&
        options.fixed_length_metres <= 0.0)
    {
        throw std::invalid_argument(
            "sample_traced_channel_link: fixed_length mode requires "
            "fixed_length_metres > 0"
        );
    }
}

void validate_traced_link(
    const lsdtt_traced_channel_link& traced_link
)
{
    if (traced_link.flats.size() != traced_link.xy.size())
    {
        throw std::runtime_error(
            "sample_traced_channel_link: traced link flats and xy sizes differ"
        );
    }

    if (traced_link.flats.size() !=
        traced_link.cumulative_distance_from_upstream_metres.size())
    {
        throw std::runtime_error(
            "sample_traced_channel_link: traced link flats and cumulative "
            "distance sizes differ"
        );
    }

    if (!traced_link.flats.empty())
    {
        if (traced_link.cumulative_distance_from_upstream_metres.front() != 0.0)
        {
            throw std::runtime_error(
                "sample_traced_channel_link: traced link cumulative distance "
                "must start at 0"
            );
        }

        const double total =
            traced_link.cumulative_distance_from_upstream_metres.back();

        if (std::abs(total - traced_link.total_length_metres) > 1e-9)
        {
            throw std::runtime_error(
                "sample_traced_channel_link: traced link total length does not "
                "match final cumulative distance"
            );
        }
    }
}

lsdtt_sampled_channel_link make_empty_sample(
    const lsdtt_traced_channel_link& traced_link
)
{
    lsdtt_sampled_channel_link sampled;
    sampled.upstream_junction   = traced_link.upstream_junction;
    sampled.downstream_junction = traced_link.downstream_junction;
    return sampled;
}

lsdtt_sampled_channel_link make_sample_from_half_open_interval(
    const lsdtt_traced_channel_link& traced_link,
    std::size_t begin_index,
    std::size_t end_index
)
{
    lsdtt_sampled_channel_link sampled;
    sampled.upstream_junction   = traced_link.upstream_junction;
    sampled.downstream_junction = traced_link.downstream_junction;

    if (begin_index >= end_index || end_index > traced_link.size())
    {
        return sampled;
    }

    sampled.flats.assign(
        traced_link.flats.begin() + static_cast<std::ptrdiff_t>(begin_index),
        traced_link.flats.begin() + static_cast<std::ptrdiff_t>(end_index)
    );

    sampled.xy.assign(
        traced_link.xy.begin() + static_cast<std::ptrdiff_t>(begin_index),
        traced_link.xy.begin() + static_cast<std::ptrdiff_t>(end_index)
    );

    sampled.cumulative_distance_metres =
        build_cumulative_distance(sampled.xy);

    sampled.sampled_length_metres =
        sampled.cumulative_distance_metres.empty()
            ? 0.0
            : sampled.cumulative_distance_metres.back();

    return sampled;
}

} // unnamed namespace

// -----------------------------------------------------------------------------
// trace_channel_link_between_junctions
// -----------------------------------------------------------------------------
lsdtt_traced_channel_link trace_channel_link_between_junctions(
    const lsdtt_flow_base&             flow,
    const lsdtt_junction_network_base& network,
    std::int32_t                       upstream_junction
)
{
    if (flow.empty())
    {
        throw std::runtime_error(
            "trace_channel_link_between_junctions: flow object is empty"
        );
    }

    if (network.empty())
    {
        throw std::runtime_error(
            "trace_channel_link_between_junctions: junction network is empty"
        );
    }

    if (!network.is_valid_junction(upstream_junction))
    {
        throw std::runtime_error(
            "trace_channel_link_between_junctions: upstream_junction is invalid"
        );
    }

    if (network.is_outlet_junction(upstream_junction))
    {
        throw std::runtime_error(
            "trace_channel_link_between_junctions: cannot trace a downstream "
            "link from an outlet junction"
        );
    }

    const std::int32_t downstream_junction =
        network.receiver_of_junction(upstream_junction);

    if (!network.is_valid_junction(downstream_junction))
    {
        throw std::runtime_error(
            "trace_channel_link_between_junctions: receiver_of_junction did not "
            "return a valid downstream junction"
        );
    }

    const std::int64_t upstream_flat =
        network.flat_of_junction(upstream_junction);

    const std::int64_t downstream_flat =
        network.flat_of_junction(downstream_junction);

    if (!flow.is_valid_cell(upstream_flat) || !flow.is_valid_cell(downstream_flat))
    {
        throw std::runtime_error(
            "trace_channel_link_between_junctions: endpoint junction flat index "
            "is not a valid flow cell"
        );
    }

    if (!network.is_channel(upstream_flat) || !network.is_channel(downstream_flat))
    {
        throw std::runtime_error(
            "trace_channel_link_between_junctions: endpoint junction flat index "
            "is not on the channel network"
        );
    }

    lsdtt_traced_channel_link traced;
    traced.upstream_junction   = upstream_junction;
    traced.downstream_junction = downstream_junction;
    traced.upstream_flat       = upstream_flat;
    traced.downstream_flat     = downstream_flat;

    const std::int64_t max_steps =
        static_cast<std::int64_t>(flow.n_rows()) *
        static_cast<std::int64_t>(flow.n_cols());

    std::int64_t current = upstream_flat;
    traced.flats.push_back(current);

    for (std::int64_t steps = 0; steps <= max_steps; ++steps)
    {
        if (current == downstream_flat)
        {
            break;
        }

        if (flow.is_outlet_unchecked(current))
        {
            throw std::runtime_error(
                "trace_channel_link_between_junctions: reached an outlet before "
                "arriving at the declared downstream junction"
            );
        }

        const std::int64_t rec = flow.receiver_of_unchecked(current);

        if (rec < 0)
        {
            throw std::runtime_error(
                "trace_channel_link_between_junctions: encountered a negative "
                "receiver before arriving at the downstream junction"
            );
        }

        if (!network.is_channel(rec))
        {
            throw std::runtime_error(
                "trace_channel_link_between_junctions: traced path left the "
                "channel network"
            );
        }

        current = rec;
        traced.flats.push_back(current);

        if (network.is_junction_flat(current) && current != downstream_flat)
        {
            throw std::runtime_error(
                "trace_channel_link_between_junctions: encountered an unexpected "
                "intermediate junction before reaching the declared receiver"
            );
        }

        if (current == downstream_flat)
        {
            break;
        }

        if (steps == max_steps)
        {
            throw std::runtime_error(
                "trace_channel_link_between_junctions: exceeded maximum number "
                "of tracing steps"
            );
        }
    }

    if (traced.flats.empty() || traced.flats.back() != downstream_flat)
    {
        throw std::runtime_error(
            "trace_channel_link_between_junctions: failed to reach the declared "
            "downstream junction"
        );
    }

    traced.xy.reserve(traced.flats.size());
    for (const std::int64_t flat : traced.flats)
    {
        const int row = flow.row_from_flat_unchecked(flat);
        const int col = flow.col_from_flat_unchecked(flat);

        traced.xy.push_back(
            cell_center_from_row_col(flow.meta(), row, col)
        );
    }

    traced.cumulative_distance_from_upstream_metres =
        build_cumulative_distance(traced.xy);

    traced.total_length_metres =
        traced.cumulative_distance_from_upstream_metres.empty()
            ? 0.0
            : traced.cumulative_distance_from_upstream_metres.back();

    return traced;
}

// -----------------------------------------------------------------------------
// sample_traced_channel_link
// -----------------------------------------------------------------------------
lsdtt_sampled_channel_link sample_traced_channel_link(
    const lsdtt_traced_channel_link&    traced_link,
    const lsdtt_link_sampling_options&  options,
    lsdtt_link_sample_anchor            anchor
)
{
    validate_sampling_options(options);
    validate_traced_link(traced_link);

    if (traced_link.empty())
    {
        return make_empty_sample(traced_link);
    }

    if (options.mode == lsdtt_link_sampling_mode::full_link)
    {
        lsdtt_sampled_channel_link sampled =
            make_sample_from_half_open_interval(
                traced_link,
                0,
                traced_link.size()
            );

        if (sampled.size() < static_cast<std::size_t>(options.min_nodes))
        {
            return make_empty_sample(traced_link);
        }

        return sampled;
    }

    if (options.mode == lsdtt_link_sampling_mode::fixed_n_pixels)
    {
        const std::size_t n_available = traced_link.size();
        const std::size_t n_requested =
            static_cast<std::size_t>(options.fixed_n_pixels);

        const std::size_t n_used = std::min(n_available, n_requested);

        std::size_t begin_index = 0;
        std::size_t end_index   = n_used;

        if (anchor == lsdtt_link_sample_anchor::downstream_end)
        {
            begin_index = n_available - n_used;
            end_index   = n_available;
        }

        lsdtt_sampled_channel_link sampled =
            make_sample_from_half_open_interval(
                traced_link,
                begin_index,
                end_index
            );

        if (sampled.size() < static_cast<std::size_t>(options.min_nodes))
        {
            return make_empty_sample(traced_link);
        }

        return sampled;
    }

    if (options.mode == lsdtt_link_sampling_mode::fixed_length)
    {
        const double limit = options.fixed_length_metres;
        const std::vector<double>& cumulative =
            traced_link.cumulative_distance_from_upstream_metres;

        std::size_t begin_index = 0;
        std::size_t end_index   = 1;

        if (anchor == lsdtt_link_sample_anchor::upstream_end)
        {
            end_index = traced_link.size();

            for (std::size_t i = 0; i < cumulative.size(); ++i)
            {
                if (cumulative[i] <= limit)
                {
                    end_index = i + 1;
                }
                else
                {
                    break;
                }
            }
        }
        else
        {
            begin_index = traced_link.size() - 1;
            end_index   = traced_link.size();

            const double total = traced_link.total_length_metres;

            for (std::size_t i = traced_link.size(); i > 0; --i)
            {
                const std::size_t idx = i - 1;
                const double distance_from_downstream = total - cumulative[idx];

                if (distance_from_downstream <= limit)
                {
                    begin_index = idx;
                }
                else
                {
                    break;
                }
            }
        }

        lsdtt_sampled_channel_link sampled =
            make_sample_from_half_open_interval(
                traced_link,
                begin_index,
                end_index
            );

        if (sampled.size() < static_cast<std::size_t>(options.min_nodes))
        {
            return make_empty_sample(traced_link);
        }

        return sampled;
    }

    throw std::runtime_error(
        "sample_traced_channel_link: unrecognised sampling mode"
    );
}

// -----------------------------------------------------------------------------
// reverse_sampled_channel_link
// -----------------------------------------------------------------------------
void reverse_sampled_channel_link(
    lsdtt_sampled_channel_link& sampled_link
)
{
    std::swap(
        sampled_link.upstream_junction,
        sampled_link.downstream_junction
    );

    std::reverse(
        sampled_link.flats.begin(),
        sampled_link.flats.end()
    );

    std::reverse(
        sampled_link.xy.begin(),
        sampled_link.xy.end()
    );

    sampled_link.cumulative_distance_metres =
        build_cumulative_distance(sampled_link.xy);

    sampled_link.sampled_length_metres =
        sampled_link.cumulative_distance_metres.empty()
            ? 0.0
            : sampled_link.cumulative_distance_metres.back();
}


// -----------------------------------------------------------------------------
// fit_channel_bearing_to_sampled_link
// -----------------------------------------------------------------------------
lsdtt_channel_bearing_fit fit_channel_bearing_to_sampled_link(
    const lsdtt_sampled_channel_link& sampled_link
)
{
    lsdtt_channel_bearing_fit result;

    if (sampled_link.flats.size() != sampled_link.xy.size())
    {
        throw std::runtime_error(
            "fit_channel_bearing_to_sampled_link: sampled link flats and xy "
            "sizes differ"
        );
    }

    if (sampled_link.flats.size() != sampled_link.cumulative_distance_metres.size())
    {
        throw std::runtime_error(
            "fit_channel_bearing_to_sampled_link: sampled link flats and "
            "cumulative distance sizes differ"
        );
    }

    result.n_nodes_used = sampled_link.size();
    result.sampled_length_metres = sampled_link.sampled_length_metres;

    if (sampled_link.size() < 2)
    {
        return result;
    }

    std::vector<double> x;
    std::vector<double> y;
    x.reserve(sampled_link.xy.size());
    y.reserve(sampled_link.xy.size());

    for (const lsdtt_xy& point : sampled_link.xy)
    {
        x.push_back(point.x);
        y.push_back(point.y);
    }

    const lsdtt_orthogonal_linear_regression_result fit =
        compute_orthogonal_linear_regression(x, y);

    if (!fit.has_defined_direction)
    {
        return result;
    }

    const double chord_dx =
        sampled_link.xy.back().x - sampled_link.xy.front().x;

    const double chord_dy =
        sampled_link.xy.back().y - sampled_link.xy.front().y;

    if (chord_dx == 0.0 && chord_dy == 0.0)
    {
        return result;
    }

    double direction_x = fit.direction_x;
    double direction_y = fit.direction_y;

    const double dot = direction_x * chord_dx + direction_y * chord_dy;

    if (dot < 0.0)
    {
        direction_x = -direction_x;
        direction_y = -direction_y;
    }

    result.valid         = true;
    result.direction_x   = direction_x;
    result.direction_y   = direction_y;
    result.fit_linearity = fit.linearity_ratio;
    result.bearing_rad   = bearing_clockwise_from_north(direction_x, direction_y);
    result.bearing_deg   = result.bearing_rad * 180.0 / std::acos(-1.0);

    return result;
}

// -----------------------------------------------------------------------------
// calculate_channel_link_bearings
// -----------------------------------------------------------------------------
std::vector<lsdtt_channel_link_bearing> calculate_channel_link_bearings(
    const lsdtt_flow_base&                    flow,
    const lsdtt_junction_network_base&        network,
    const std::vector<std::int32_t>&          upstream_junctions_to_analyse,
    const lsdtt_link_sampling_options&        options,
    lsdtt_link_sample_anchor                  anchor
)
{
    if (flow.empty())
    {
        throw std::runtime_error(
            "calculate_channel_link_bearings: flow object is empty"
        );
    }

    if (network.empty())
    {
        throw std::runtime_error(
            "calculate_channel_link_bearings: junction network is empty"
        );
    }

    const std::vector<std::int32_t> junctions =
        resolve_upstream_junctions_to_analyse(
            network,
            upstream_junctions_to_analyse
        );

    std::vector<lsdtt_channel_link_bearing> records;
    records.reserve(junctions.size());

    for (const std::int32_t upstream_junction : junctions)
    {
        if (network.is_outlet_junction(upstream_junction))
        {
            continue;
        }

        const lsdtt_traced_channel_link traced_link =
            trace_channel_link_between_junctions(
                flow,
                network,
                upstream_junction
            );

        const lsdtt_sampled_channel_link sampled_link =
            sample_traced_channel_link(
                traced_link,
                options,
                anchor
            );

        const lsdtt_channel_bearing_fit fit =
            fit_channel_bearing_to_sampled_link(sampled_link);

        if (!fit.valid)
        {
            continue;
        }

        lsdtt_channel_link_bearing record;
        record.upstream_junction   = traced_link.upstream_junction;
        record.downstream_junction = traced_link.downstream_junction;
        record.upstream_flat       = traced_link.upstream_flat;
        record.downstream_flat     = traced_link.downstream_flat;

        record.stream_order =
            network.stream_order_of_junction(upstream_junction);

        record.drainage_area_m2 =
            static_cast<double>(flow.s_vector_size(traced_link.upstream_flat)) *
            flow.cell_size() *
            flow.cell_size();

        record.n_nodes_total       = traced_link.size();
        record.n_nodes_used        = fit.n_nodes_used;
        record.total_length_metres = traced_link.total_length_metres;
        record.sampled_length_metres = fit.sampled_length_metres;

        record.direction_x   = fit.direction_x;
        record.direction_y   = fit.direction_y;
        record.bearing_rad   = fit.bearing_rad;
        record.bearing_deg   = fit.bearing_deg;
        record.fit_linearity = fit.fit_linearity;

        if (!sampled_link.empty())
        {
            const std::size_t mid =
                sampled_link.size() / 2;

            record.representative_flat =
                sampled_link.flats[mid];

            record.representative_x =
                sampled_link.xy[mid].x;

            record.representative_y =
                sampled_link.xy[mid].y;
        }

        records.push_back(std::move(record));
    }

    return records;
}

// -----------------------------------------------------------------------------
// IMPORTANT TERMINOLOGY NOTE - JUNCTION-ANGLE NAMING CONVENTION
//
// This routine uses internally descriptive variable names for the three angular
// sectors around a two-donor junction, but these names do NOT exactly match the
// historical naming convention used in earlier LSDTopoTools work.
//
// In the naming convention used by Strong and Mudd (PNAS, 2023):
//
//   - the angle BETWEEN the two donor channels was called the
//     "junction angle" and denoted alpha
//
//   - the two angles between each donor and the receiver channel were called
//     the two "bending angles" and denoted beta_1 and beta_2
//
// In THIS implementation, the stored names are:
//
//   donor_opening_angle
//       = the angle between the two donor channels
//       = legacy "junction angle" alpha
//
//   donor_1_receiver_angle
//       = the angle between donor_1 and the receiver channel
//       = legacy bending angle beta_1
//
//   donor_2_receiver_angle
//       = the angle between donor_2 and the receiver channel
//       = legacy bending angle beta_2
//
// So the mapping is:
//
//   donor_opening_angle      <->  alpha
//   donor_1_receiver_angle   <->  beta_1
//   donor_2_receiver_angle   <->  beta_2
//
// The term "opening angle" is used here because it is geometrically explicit:
// it is simply the angular opening between the two incoming donor links.
//
// The term "receiver angle" is also used explicitly because each of the other
// two sectors is defined by one donor and the receiver link.
//
// IMPORTANT:
// This code may also report an additional metric called:
//
//   bending_angle_signed
//   bending_angle_abs
//
// This is NOT part of the original Strong and Mudd (PNAS, 2023) reported
// alpha / beta_1 / beta_2 set.
//
// Instead, it is an additional diagnostic quantity defined here as the angular
// offset between:
//   (1) the inward-pointing receiver direction, and
//   (2) the bisector of the donor opening angle.
//
// This extra bending metric should therefore be treated as a new derived
// diagnostic, not as a replacement for either beta_1 or beta_2.
//
// DONOR LABELS:
// donor_1 and donor_2 are implementation labels only. They are stable, but
// they are not guaranteed to mean "left" and "right" tributary in map view.
// In the current junction-network implementation, donor ordering is determined
// by the deterministic ordering of direct donor junction ids, not by geometry.
// Therefore:
//
//   donor_1_receiver_angle and donor_2_receiver_angle
//
// are stable and reproducible, but beta_1 and beta_2 should not be interpreted
// as left-bank versus right-bank bending angles unless an additional geometric
// re-labelling step is performed.
//
// BEARING DIRECTION CONVENTION:
// For angle calculations, donor bearings are treated as inward-pointing toward
// the analysed junction. The receiver link is traced downstream away from the
// junction, so its outward bearing is reversed internally when constructing the
// three-sector junction geometry.
//
// In summary:
//
//   internal name              historical Strong & Mudd name
//   --------------------------------------------------------
//   donor_opening_angle      = alpha
//   donor_1_receiver_angle   = beta_1
//   donor_2_receiver_angle   = beta_2
//
// This mapping is intentional and must be kept clearly documented anywhere
// junction-angle results are exposed to users.
// -----------------------------------------------------------------------------
// calculate_junction_angles
// -----------------------------------------------------------------------------
std::vector<lsdtt_junction_angle_record> calculate_junction_angles(
    const lsdtt_flow_base&                    flow,
    const lsdtt_junction_network_base&        network,
    const std::vector<std::int32_t>&          junctions_to_analyse,
    const lsdtt_link_sampling_options&        options
)
{
    if (flow.empty())
    {
        throw std::runtime_error(
            "calculate_junction_angles: flow object is empty"
        );
    }

    if (network.empty())
    {
        throw std::runtime_error(
            "calculate_junction_angles: junction network is empty"
        );
    }

    const std::vector<std::int32_t> junctions =
        resolve_junctions_to_analyse(
            network,
            junctions_to_analyse
        );

    std::vector<lsdtt_junction_angle_record> records;
    records.reserve(junctions.size());

    const double pi = std::acos(-1.0);
    const double two_pi = 2.0 * pi;

    for (const std::int32_t junction : junctions)
    {
        if (!network.is_valid_junction(junction))
        {
            continue;
        }

        if (network.is_outlet_junction(junction))
        {
            continue;
        }

        const lsdtt_junction_index_slice donors =
            network.donor_junctions(junction);

        if (donors.size() != 2)
        {
            continue;
        }

        const std::int32_t donor_1 = donors[0];
        const std::int32_t donor_2 = donors[1];
        const std::int32_t receiver = network.receiver_of_junction(junction);

        if (!network.is_valid_junction(donor_1) ||
            !network.is_valid_junction(donor_2) ||
            !network.is_valid_junction(receiver))
        {
            continue;
        }

        lsdtt_traced_channel_link donor_1_traced;
        lsdtt_traced_channel_link donor_2_traced;
        lsdtt_traced_channel_link receiver_traced;

        lsdtt_sampled_channel_link donor_1_sampled;
        lsdtt_sampled_channel_link donor_2_sampled;
        lsdtt_sampled_channel_link receiver_sampled;

        lsdtt_channel_bearing_fit donor_1_fit;
        lsdtt_channel_bearing_fit donor_2_fit;
        lsdtt_channel_bearing_fit receiver_fit;

        try
        {
            donor_1_traced =
                trace_channel_link_between_junctions(
                    flow,
                    network,
                    donor_1
                );

            donor_2_traced =
                trace_channel_link_between_junctions(
                    flow,
                    network,
                    donor_2
                );

            receiver_traced =
                trace_channel_link_between_junctions(
                    flow,
                    network,
                    junction
                );

            donor_1_sampled =
                sample_traced_channel_link(
                    donor_1_traced,
                    options,
                    lsdtt_link_sample_anchor::downstream_end
                );

            donor_2_sampled =
                sample_traced_channel_link(
                    donor_2_traced,
                    options,
                    lsdtt_link_sample_anchor::downstream_end
                );

            receiver_sampled =
                sample_traced_channel_link(
                    receiver_traced,
                    options,
                    lsdtt_link_sample_anchor::upstream_end
                );

            donor_1_fit =
                fit_channel_bearing_to_sampled_link(donor_1_sampled);

            donor_2_fit =
                fit_channel_bearing_to_sampled_link(donor_2_sampled);

            receiver_fit =
                fit_channel_bearing_to_sampled_link(receiver_sampled);
        }
        catch (const std::exception&)
        {
            continue;
        }

        if (!donor_1_fit.valid || !donor_2_fit.valid || !receiver_fit.valid)
        {
            continue;
        }

        // Donor bearings already point inward toward the analysed junction.
        const double donor_1_bearing =
            wrap_angle_0_2pi(donor_1_fit.bearing_rad);

        const double donor_2_bearing =
            wrap_angle_0_2pi(donor_2_fit.bearing_rad);

        // Receiver fit points outward/downstream away from the analysed junction.
        const double receiver_bearing_out =
            wrap_angle_0_2pi(receiver_fit.bearing_rad);

        // Reverse by pi so that the receiver also points inward toward the junction.
        const double receiver_bearing_in =
            wrap_angle_0_2pi(receiver_bearing_out + pi);

        // Rotate so that inward receiver direction is zero.
        const double donor_1_rot =
            wrap_angle_0_2pi(donor_1_bearing - receiver_bearing_in);

        const double donor_2_rot =
            wrap_angle_0_2pi(donor_2_bearing - receiver_bearing_in);

        double donor_opening_angle = std::numeric_limits<double>::quiet_NaN();
        double donor_1_receiver_angle = std::numeric_limits<double>::quiet_NaN();
        double donor_2_receiver_angle = std::numeric_limits<double>::quiet_NaN();
        double opening_midpoint_rot = std::numeric_limits<double>::quiet_NaN();

        if (donor_2_rot > donor_1_rot)
        {
            donor_1_receiver_angle = donor_1_rot;
            donor_2_receiver_angle = two_pi - donor_2_rot;
            donor_opening_angle    = donor_2_rot - donor_1_rot;
            opening_midpoint_rot   = donor_1_rot + 0.5 * donor_opening_angle;
        }
        else
        {
            donor_1_receiver_angle = two_pi - donor_1_rot;
            donor_2_receiver_angle = donor_2_rot;
            donor_opening_angle    = donor_1_rot - donor_2_rot;
            opening_midpoint_rot   = donor_2_rot + 0.5 * donor_opening_angle;
        }

        const double angle_sum =
            donor_opening_angle +
            donor_1_receiver_angle +
            donor_2_receiver_angle;

        if (!std::isfinite(angle_sum) ||
            std::abs(angle_sum - two_pi) > 1e-6)
        {
            continue;
        }

        const double bending_angle_signed =
            wrap_angle_minus_pi_pi(opening_midpoint_rot);

        const double bending_angle_abs =
            std::abs(bending_angle_signed);

        lsdtt_junction_angle_record record;
        record.valid = true;

        record.junction = junction;
        record.donor_1  = donor_1;
        record.donor_2  = donor_2;
        record.receiver = receiver;

        record.junction_order = network.stream_order_of_junction(junction);
        record.donor_1_order  = network.stream_order_of_junction(donor_1);
        record.donor_2_order  = network.stream_order_of_junction(donor_2);
        record.receiver_order = network.stream_order_of_junction(receiver);

        record.junction_area_m2 =
            drainage_area_from_flat(flow, network.flat_of_junction(junction));

        // Tributary areas just upstream of the analysed junction.
        if (donor_1_traced.size() >= 2)
        {
            record.donor_1_area_m2 =
                drainage_area_from_flat(
                    flow,
                    donor_1_traced.flats[donor_1_traced.size() - 2]
                );
        }

        if (donor_2_traced.size() >= 2)
        {
            record.donor_2_area_m2 =
                drainage_area_from_flat(
                    flow,
                    donor_2_traced.flats[donor_2_traced.size() - 2]
                );
        }

        record.donor_1_bearing_rad = donor_1_bearing;
        record.donor_2_bearing_rad = donor_2_bearing;
        record.receiver_bearing_out_rad = receiver_bearing_out;
        record.receiver_bearing_in_rad  = receiver_bearing_in;

        record.donor_1_bearing_deg = radians_to_degrees(donor_1_bearing);
        record.donor_2_bearing_deg = radians_to_degrees(donor_2_bearing);
        record.receiver_bearing_out_deg = radians_to_degrees(receiver_bearing_out);
        record.receiver_bearing_in_deg  = radians_to_degrees(receiver_bearing_in);

        record.donor_opening_angle_rad = donor_opening_angle;
        record.donor_1_receiver_angle_rad = donor_1_receiver_angle;
        record.donor_2_receiver_angle_rad = donor_2_receiver_angle;

        record.donor_opening_angle_deg = radians_to_degrees(donor_opening_angle);
        record.donor_1_receiver_angle_deg = radians_to_degrees(donor_1_receiver_angle);
        record.donor_2_receiver_angle_deg = radians_to_degrees(donor_2_receiver_angle);

        record.bending_angle_signed_rad = bending_angle_signed;
        record.bending_angle_abs_rad    = bending_angle_abs;

        record.bending_angle_signed_deg = radians_to_degrees(bending_angle_signed);
        record.bending_angle_abs_deg    = radians_to_degrees(bending_angle_abs);

        record.donor_1_linearity = donor_1_fit.fit_linearity;
        record.donor_2_linearity = donor_2_fit.fit_linearity;
        record.receiver_linearity = receiver_fit.fit_linearity;

        record.donor_1_nodes_total = donor_1_traced.size();
        record.donor_2_nodes_total = donor_2_traced.size();
        record.receiver_nodes_total = receiver_traced.size();

        record.donor_1_nodes_used = donor_1_fit.n_nodes_used;
        record.donor_2_nodes_used = donor_2_fit.n_nodes_used;
        record.receiver_nodes_used = receiver_fit.n_nodes_used;

        record.donor_1_total_length_metres = donor_1_traced.total_length_metres;
        record.donor_2_total_length_metres = donor_2_traced.total_length_metres;
        record.receiver_total_length_metres = receiver_traced.total_length_metres;

        record.donor_1_sampled_length_metres = donor_1_fit.sampled_length_metres;
        record.donor_2_sampled_length_metres = donor_2_fit.sampled_length_metres;
        record.receiver_sampled_length_metres = receiver_fit.sampled_length_metres;

        records.push_back(std::move(record));
    }

    return records;
}

// -----------------------------------------------------------------------------
// write_channel_link_bearings
// -----------------------------------------------------------------------------
void write_channel_link_bearings(
    const std::vector<lsdtt_channel_link_bearing>& records,
    const lsdtt_flow_base&                         flow,
    const lsdtt_junction_network_base&             network,
    const std::string&                             out_path,
    point_output_format                            format,
    const std::string&                             layer_name,
    bool                                           spatial_index
)
{
    if (flow.empty())
    {
        throw std::runtime_error(
            "write_channel_link_bearings: flow object is empty"
        );
    }

    if (network.empty())
    {
        throw std::runtime_error(
            "write_channel_link_bearings: junction network is empty"
        );
    }

    lsdtt_data_loader loader;

    // -------------------------------------------------------------------------
    // CSV output
    //
    // This is tabular output, not line geometry. We therefore export upstream
    // and downstream coordinates as latitude/longitude columns.
    // -------------------------------------------------------------------------
    if (format == point_output_format::csv)
    {
        if (flow.meta().crs_wkt.empty())
        {
            throw std::runtime_error(
                "write_channel_link_bearings: CSV output requires a valid raster "
                "CRS so coordinates can be transformed to EPSG:4326"
            );
        }

        const lsdtt_coordinate_transform transform =
            lsdtt_coordinate_transform::from_epsg(flow.meta().crs_wkt, 4326);

        std::vector<std::string> headers =
        {
            "upstream_junction",
            "downstream_junction",
            "stream_order",
            "drainage_area_m2",
            "n_nodes_total",
            "n_nodes_used",
            "total_length_metres",
            "sampled_length_metres",
            "direction_x",
            "direction_y",
            "bearing_rad",
            "bearing_deg",
            "fit_linearity",
            "upstream_latitude",
            "upstream_longitude",
            "downstream_latitude",
            "downstream_longitude"
        };

        std::vector<std::vector<std::string>> rows;
        rows.reserve(records.size());

        for (const auto& record : records)
        {
            if (record.empty())
            {
                continue;
            }

            const lsdtt_lonlat ll_up =
                flat_to_lonlat(flow, record.upstream_flat, transform);

            const lsdtt_lonlat ll_down =
                flat_to_lonlat(flow, record.downstream_flat, transform);

            rows.push_back({
                std::to_string(record.upstream_junction),
                std::to_string(record.downstream_junction),
                std::to_string(record.stream_order),
                precise_string(record.drainage_area_m2),
                std::to_string(record.n_nodes_total),
                std::to_string(record.n_nodes_used),
                precise_string(record.total_length_metres),
                precise_string(record.sampled_length_metres),
                precise_string(record.direction_x),
                precise_string(record.direction_y),
                precise_string(record.bearing_rad),
                precise_string(record.bearing_deg),
                precise_string(record.fit_linearity),
                precise_string(ll_up.latitude),
                precise_string(ll_up.longitude),
                precise_string(ll_down.latitude),
                precise_string(ll_down.longitude)
            });
        }

        loader.write_csv(out_path, headers, rows);
        return;
    }

    // -------------------------------------------------------------------------
    // Vector output: one line per analysed link, written in raster CRS.
    // Geometry uses the full compact channel link, not just the sampled subset.
    // -------------------------------------------------------------------------
    const vector_output_format out_format =
        convert_point_to_vector_output_format(format);

    lsdtt_vector_schema schema;
    schema.field_names = {
        "up_junction",
        "down_junction",
        "stream_order",
        "drainage_area_m2",
        "n_nodes_total",
        "n_nodes_used",
        "len_total_m",
        "len_sample_m",
        "dir_x",
        "dir_y",
        "bearing_rad",
        "bearing_deg",
        "linearity"
    };
    schema.field_types = {
        lsdtt_vector_field_type::int64_type,
        lsdtt_vector_field_type::int64_type,
        lsdtt_vector_field_type::int64_type,
        lsdtt_vector_field_type::double_type,
        lsdtt_vector_field_type::int64_type,
        lsdtt_vector_field_type::int64_type,
        lsdtt_vector_field_type::double_type,
        lsdtt_vector_field_type::double_type,
        lsdtt_vector_field_type::double_type,
        lsdtt_vector_field_type::double_type,
        lsdtt_vector_field_type::double_type,
        lsdtt_vector_field_type::double_type,
        lsdtt_vector_field_type::double_type
    };

    loader.begin_vector_write(
        out_path,
        layer_name,
        schema,
        flow.meta().crs_wkt,
        out_format,
        spatial_index
    );

    for (const auto& record : records)
    {
        if (record.empty())
        {
            continue;
        }

        const lsdtt_traced_channel_link traced =
            trace_channel_link_between_junctions(
                flow,
                network,
                record.upstream_junction
            );

        loader.write_line(
            traced.xy,
            {
                static_cast<std::int64_t>(record.upstream_junction),
                static_cast<std::int64_t>(record.downstream_junction),
                static_cast<std::int64_t>(record.stream_order),
                record.drainage_area_m2,
                static_cast<std::int64_t>(record.n_nodes_total),
                static_cast<std::int64_t>(record.n_nodes_used),
                record.total_length_metres,
                record.sampled_length_metres,
                record.direction_x,
                record.direction_y,
                record.bearing_rad,
                record.bearing_deg,
                record.fit_linearity
            }
        );
    }

    loader.finish_vector_write();
}


// -----------------------------------------------------------------------------
// write_junction_angle_data
// -----------------------------------------------------------------------------
void write_junction_angle_data(
    const std::vector<lsdtt_junction_angle_record>& records,
    const lsdtt_flow_base&                          flow,
    const lsdtt_junction_network_base&              network,
    const std::string&                              out_path,
    point_output_format                             format,
    const std::string&                              layer_name,
    bool                                            spatial_index
)
{
    if (flow.empty())
    {
        throw std::runtime_error(
            "write_junction_angle_data: flow object is empty"
        );
    }

    if (network.empty())
    {
        throw std::runtime_error(
            "write_junction_angle_data: junction network is empty"
        );
    }

    lsdtt_data_loader loader;

    // -------------------------------------------------------------------------
    // CSV output
    // -------------------------------------------------------------------------
    if (format == point_output_format::csv)
    {
        if (flow.meta().crs_wkt.empty())
        {
            throw std::runtime_error(
                "write_junction_angle_data: CSV output requires a valid raster "
                "CRS so coordinates can be transformed to EPSG:4326"
            );
        }

        const lsdtt_coordinate_transform transform =
            lsdtt_coordinate_transform::from_epsg(flow.meta().crs_wkt, 4326);

        std::vector<std::string> headers =
        {
            "junction",
            "donor_1",
            "donor_2",
            "receiver",
            "junction_order",
            "donor_1_order",
            "donor_2_order",
            "receiver_order",
            "junction_area_m2",
            "donor_1_area_m2",
            "donor_2_area_m2",

            "donor_1_bearing_rad",
            "donor_2_bearing_rad",
            "receiver_bearing_out_rad",
            "receiver_bearing_in_rad",

            "donor_1_bearing_deg",
            "donor_2_bearing_deg",
            "receiver_bearing_out_deg",
            "receiver_bearing_in_deg",

            "donor_opening_angle_rad",
            "donor_1_receiver_angle_rad",
            "donor_2_receiver_angle_rad",

            "donor_opening_angle_deg",
            "donor_1_receiver_angle_deg",
            "donor_2_receiver_angle_deg",

            "alpha_rad",
            "beta_1_rad",
            "beta_2_rad",
            "alpha_deg",
            "beta_1_deg",
            "beta_2_deg",

            "bending_angle_signed_rad",
            "bending_angle_abs_rad",
            "bending_angle_signed_deg",
            "bending_angle_abs_deg",

            "donor_1_linearity",
            "donor_2_linearity",
            "receiver_linearity",

            "donor_1_nodes_total",
            "donor_2_nodes_total",
            "receiver_nodes_total",
            "donor_1_nodes_used",
            "donor_2_nodes_used",
            "receiver_nodes_used",

            "donor_1_total_length_metres",
            "donor_2_total_length_metres",
            "receiver_total_length_metres",
            "donor_1_sampled_length_metres",
            "donor_2_sampled_length_metres",
            "receiver_sampled_length_metres",

            "latitude",
            "longitude"
        };

        std::vector<std::vector<std::string>> rows;
        rows.reserve(records.size());

        for (const auto& record : records)
        {
            if (record.empty())
            {
                continue;
            }

            const lsdtt_lonlat ll =
                flat_to_lonlat(
                    flow,
                    network.flat_of_junction(record.junction),
                    transform
                );

            rows.push_back({
                std::to_string(record.junction),
                std::to_string(record.donor_1),
                std::to_string(record.donor_2),
                std::to_string(record.receiver),
                std::to_string(record.junction_order),
                std::to_string(record.donor_1_order),
                std::to_string(record.donor_2_order),
                std::to_string(record.receiver_order),
                precise_string(record.junction_area_m2),
                precise_string(record.donor_1_area_m2),
                precise_string(record.donor_2_area_m2),

                precise_string(record.donor_1_bearing_rad),
                precise_string(record.donor_2_bearing_rad),
                precise_string(record.receiver_bearing_out_rad),
                precise_string(record.receiver_bearing_in_rad),

                precise_string(record.donor_1_bearing_deg),
                precise_string(record.donor_2_bearing_deg),
                precise_string(record.receiver_bearing_out_deg),
                precise_string(record.receiver_bearing_in_deg),

                precise_string(record.donor_opening_angle_rad),
                precise_string(record.donor_1_receiver_angle_rad),
                precise_string(record.donor_2_receiver_angle_rad),

                precise_string(record.donor_opening_angle_deg),
                precise_string(record.donor_1_receiver_angle_deg),
                precise_string(record.donor_2_receiver_angle_deg),

                precise_string(record.donor_opening_angle_rad),
                precise_string(record.donor_1_receiver_angle_rad),
                precise_string(record.donor_2_receiver_angle_rad),
                precise_string(record.donor_opening_angle_deg),
                precise_string(record.donor_1_receiver_angle_deg),
                precise_string(record.donor_2_receiver_angle_deg),

                precise_string(record.bending_angle_signed_rad),
                precise_string(record.bending_angle_abs_rad),
                precise_string(record.bending_angle_signed_deg),
                precise_string(record.bending_angle_abs_deg),

                precise_string(record.donor_1_linearity),
                precise_string(record.donor_2_linearity),
                precise_string(record.receiver_linearity),

                std::to_string(record.donor_1_nodes_total),
                std::to_string(record.donor_2_nodes_total),
                std::to_string(record.receiver_nodes_total),
                std::to_string(record.donor_1_nodes_used),
                std::to_string(record.donor_2_nodes_used),
                std::to_string(record.receiver_nodes_used),

                precise_string(record.donor_1_total_length_metres),
                precise_string(record.donor_2_total_length_metres),
                precise_string(record.receiver_total_length_metres),
                precise_string(record.donor_1_sampled_length_metres),
                precise_string(record.donor_2_sampled_length_metres),
                precise_string(record.receiver_sampled_length_metres),

                precise_string(ll.latitude),
                precise_string(ll.longitude)
            });
        }

        loader.write_csv(out_path, headers, rows);
        return;
    }

    // -------------------------------------------------------------------------
    // Vector output: one point per analysed junction
    // -------------------------------------------------------------------------
    const vector_output_format out_format =
        convert_point_to_vector_output_format(format);

    lsdtt_vector_schema schema;
    schema.field_names = {
        "junction",
        "donor_1",
        "donor_2",
        "receiver",
        "j_order",
        "d1_order",
        "d2_order",
        "r_order",
        "j_area_m2",
        "d1_area_m2",
        "d2_area_m2",

        "d1_bear_deg",
        "d2_bear_deg",
        "rout_b_deg",
        "rin_b_deg",

        "open_deg",
        "d1r_deg",
        "d2r_deg",

        "alpha_deg",
        "beta1_deg",
        "beta2_deg",

        "bend_s_deg",
        "bend_a_deg",

        "d1_lin",
        "d2_lin",
        "r_lin"
    };
    schema.field_types = {
        lsdtt_vector_field_type::int64_type,
        lsdtt_vector_field_type::int64_type,
        lsdtt_vector_field_type::int64_type,
        lsdtt_vector_field_type::int64_type,
        lsdtt_vector_field_type::int64_type,
        lsdtt_vector_field_type::int64_type,
        lsdtt_vector_field_type::int64_type,
        lsdtt_vector_field_type::int64_type,
        lsdtt_vector_field_type::double_type,
        lsdtt_vector_field_type::double_type,
        lsdtt_vector_field_type::double_type,

        lsdtt_vector_field_type::double_type,
        lsdtt_vector_field_type::double_type,
        lsdtt_vector_field_type::double_type,
        lsdtt_vector_field_type::double_type,

        lsdtt_vector_field_type::double_type,
        lsdtt_vector_field_type::double_type,
        lsdtt_vector_field_type::double_type,

        lsdtt_vector_field_type::double_type,
        lsdtt_vector_field_type::double_type,
        lsdtt_vector_field_type::double_type,

        lsdtt_vector_field_type::double_type,
        lsdtt_vector_field_type::double_type,

        lsdtt_vector_field_type::double_type,
        lsdtt_vector_field_type::double_type,
        lsdtt_vector_field_type::double_type
    };

    loader.begin_vector_write(
        out_path,
        layer_name,
        schema,
        flow.meta().crs_wkt,
        out_format,
        spatial_index
    );

    for (const auto& record : records)
    {
        if (record.empty())
        {
            continue;
        }

        const lsdtt_xy xy =
            flat_to_xy(
                flow,
                network.flat_of_junction(record.junction)
            );

        loader.write_point(
            xy.x,
            xy.y,
            {
                static_cast<std::int64_t>(record.junction),
                static_cast<std::int64_t>(record.donor_1),
                static_cast<std::int64_t>(record.donor_2),
                static_cast<std::int64_t>(record.receiver),
                static_cast<std::int64_t>(record.junction_order),
                static_cast<std::int64_t>(record.donor_1_order),
                static_cast<std::int64_t>(record.donor_2_order),
                static_cast<std::int64_t>(record.receiver_order),
                record.junction_area_m2,
                record.donor_1_area_m2,
                record.donor_2_area_m2,

                record.donor_1_bearing_deg,
                record.donor_2_bearing_deg,
                record.receiver_bearing_out_deg,
                record.receiver_bearing_in_deg,

                record.donor_opening_angle_deg,
                record.donor_1_receiver_angle_deg,
                record.donor_2_receiver_angle_deg,

                record.donor_opening_angle_deg,       // alpha
                record.donor_1_receiver_angle_deg,    // beta_1
                record.donor_2_receiver_angle_deg,    // beta_2

                record.bending_angle_signed_deg,
                record.bending_angle_abs_deg,

                record.donor_1_linearity,
                record.donor_2_linearity,
                record.receiver_linearity
            }
        );
    }

    loader.finish_vector_write();
}


// -----------------------------------------------------------------------------
// write_channel_link_bearing_vectors
// -----------------------------------------------------------------------------
void write_channel_link_bearing_vectors(
    const std::vector<lsdtt_channel_link_bearing>& records,
    const lsdtt_flow_base&                         flow,
    const std::string&                             out_path,
    vector_output_format                           format,
    const std::string&                             layer_name,
    bool                                           spatial_index
)
{
    const double mean_link_length =
        mean_total_link_length_from_bearing_records(records);

    if (!std::isfinite(mean_link_length) || mean_link_length <= 0.0)
    {
        throw std::runtime_error(
            "write_channel_link_bearing_vectors: could not determine mean link "
            "length from bearing records"
        );
    }

    lsdtt_data_loader loader;

    lsdtt_vector_schema schema;
    schema.field_names = {
        "up_junction",
        "down_junction",
        "stream_order",
        "bearing_rad",
        "bearing_deg",
        "linearity",
        "link_length_m",
        "vector_len_m"
    };
    schema.field_types = {
        lsdtt_vector_field_type::int64_type,
        lsdtt_vector_field_type::int64_type,
        lsdtt_vector_field_type::int64_type,
        lsdtt_vector_field_type::double_type,
        lsdtt_vector_field_type::double_type,
        lsdtt_vector_field_type::double_type,
        lsdtt_vector_field_type::double_type,
        lsdtt_vector_field_type::double_type
    };

    // CRS is not carried by the records, so we infer it from the representative
    // point coordinates already stored in the record. These were written in the
    // raster CRS during analysis. The caller should therefore use this writer
    // only with records generated from a valid projected DEM and channel network.
    //
    // Because the records themselves do not store WKT, this writer opens an
    // empty string CRS. If you want the file CRS embedded, the cleaner long-term
    // extension is to add a CRS string parameter here or store CRS WKT alongside
    // the records.
    loader.begin_vector_write(
        out_path,
        layer_name,
        schema,
        flow.meta().crs_wkt,
        format,
        spatial_index
    );


    const double global_vector_scale = mean_link_length / 3.0;

    for (const auto& record : records)
    {
        if (record.empty())
        {
            continue;
        }

        if (!std::isfinite(record.total_length_metres) ||
            record.total_length_metres <= 0.0)
        {
            continue;
        }

        if (!std::isfinite(record.direction_x) ||
            !std::isfinite(record.direction_y) ||
            !std::isfinite(record.representative_x) ||
            !std::isfinite(record.representative_y))
        {
            continue;
        }

        const double vector_length =
            std::min(global_vector_scale, 0.5 * record.total_length_metres);

        if (!(vector_length > 0.0))
        {
            continue;
        }

        const std::vector<lsdtt_xy> points =
            downstream_vector_geometry(
                record.representative_x,
                record.representative_y,
                record.direction_x,
                record.direction_y,
                vector_length
            );

        loader.write_line(
            points,
            {
                static_cast<std::int64_t>(record.upstream_junction),
                static_cast<std::int64_t>(record.downstream_junction),
                static_cast<std::int64_t>(record.stream_order),
                record.bearing_rad,
                record.bearing_deg,
                record.fit_linearity,
                record.total_length_metres,
                vector_length
            }
        );
    }

    loader.finish_vector_write();
}


// -----------------------------------------------------------------------------
// write_junction_angle_bearing_vectors
// -----------------------------------------------------------------------------
void write_junction_angle_bearing_vectors(
    const std::vector<lsdtt_junction_angle_record>& records,
    const lsdtt_flow_base&                          flow,
    const lsdtt_junction_network_base&              network,
    const std::string&                              out_path,
    vector_output_format                            format,
    const std::string&                              layer_name,
    bool                                            spatial_index
)
{
    if (flow.empty())
    {
        throw std::runtime_error(
            "write_junction_angle_bearing_vectors: flow object is empty"
        );
    }

    if (network.empty())
    {
        throw std::runtime_error(
            "write_junction_angle_bearing_vectors: junction network is empty"
        );
    }

    const double mean_link_length =
        mean_total_link_length_from_junction_angle_records(records);

    if (!std::isfinite(mean_link_length) || mean_link_length <= 0.0)
    {
        throw std::runtime_error(
            "write_junction_angle_bearing_vectors: could not determine mean "
            "link length from junction-angle records"
        );
    }

    lsdtt_data_loader loader;

    lsdtt_vector_schema schema;
    schema.field_names = {
        "junction",
        "arm_type",
        "arm_code",
        "bearing_rad",
        "bearing_deg",
        "display_len_m",
        "junction_order",
        "donor_1_order",
        "donor_2_order",
        "receiver_order"
    };
    schema.field_types = {
        lsdtt_vector_field_type::int64_type,
        lsdtt_vector_field_type::string_type,
        lsdtt_vector_field_type::int64_type,
        lsdtt_vector_field_type::double_type,
        lsdtt_vector_field_type::double_type,
        lsdtt_vector_field_type::double_type,
        lsdtt_vector_field_type::int64_type,
        lsdtt_vector_field_type::int64_type,
        lsdtt_vector_field_type::int64_type,
        lsdtt_vector_field_type::int64_type
    };

    loader.begin_vector_write(
        out_path,
        layer_name,
        schema,
        flow.meta().crs_wkt,
        format,
        spatial_index
    );

    const double global_display_scale = mean_link_length / 3.0;

    for (const auto& record : records)
    {
        if (record.empty())
        {
            continue;
        }

        if (!network.is_valid_junction(record.junction))
        {
            continue;
        }

        const lsdtt_xy centre =
            flat_to_xy(
                flow,
                network.flat_of_junction(record.junction)
            );

        struct arm_info
        {
            const char* arm_type;
            std::int64_t arm_code;
            bool is_receiver;
            double bearing_rad;
            double bearing_deg;
            double total_length_metres;
        };


        // Visual junction-angle QA convention:
        //
        //   - donor bearings are stored as inward-pointing toward the junction,
        //     so donor arms terminate at the junction and extend only upstream
        //
        //   - receiver_bearing_out points downstream away from the junction,
        //     so the receiver arm starts at the junction and extends only
        //     downstream
        //
        // This avoids the misleading centred rendering in which arms extend to
        // both sides of the junction.
        const arm_info arms[3] =
        {
            {
                "donor_1",
                1,
                false,
                record.donor_1_bearing_rad,
                record.donor_1_bearing_deg,
                record.donor_1_total_length_metres
            },
            {
                "donor_2",
                2,
                false,
                record.donor_2_bearing_rad,
                record.donor_2_bearing_deg,
                record.donor_2_total_length_metres
            },
            {
                "receiver",
                3,
                true,
                record.receiver_bearing_out_rad,
                record.receiver_bearing_out_deg,
                record.receiver_total_length_metres
            }
        };


        for (const arm_info& arm : arms)
        {
            if (!std::isfinite(arm.bearing_rad) ||
                !std::isfinite(arm.total_length_metres) ||
                arm.total_length_metres <= 0.0)
            {
                continue;
            }

            const double display_length =
                std::min(global_display_scale, 0.5 * arm.total_length_metres);

            if (!(display_length > 0.0))
            {
                continue;
            }

            std::vector<lsdtt_xy> points;

            if (arm.is_receiver)
            {
                points =
                    receiver_arm_geometry_from_outward_bearing(
                        centre.x,
                        centre.y,
                        arm.bearing_rad,
                        display_length
                    );
            }
            else
            {
                points =
                    donor_arm_geometry_from_inward_bearing(
                        centre.x,
                        centre.y,
                        arm.bearing_rad,
                        display_length
                    );
            }


            loader.write_line(
                points,
                {
                    static_cast<std::int64_t>(record.junction),
                    std::string(arm.arm_type),
                    arm.arm_code,
                    arm.bearing_rad,
                    arm.bearing_deg,
                    display_length,
                    static_cast<std::int64_t>(record.junction_order),
                    static_cast<std::int64_t>(record.donor_1_order),
                    static_cast<std::int64_t>(record.donor_2_order),
                    static_cast<std::int64_t>(record.receiver_order)
                }
            );
        }
    }

    loader.finish_vector_write();
}


// =============================================================================
// bearing_from_north_degrees
// =============================================================================
double bearing_from_north_degrees(double ox, double oy, double cx, double cy)
{
    const double de = cx - ox;   // easting difference
    const double dn = cy - oy;   // northing difference
    if (de == 0.0 && dn == 0.0)
        return std::numeric_limits<double>::quiet_NaN();
    double b = std::atan2(de, dn) * (180.0 / M_PI);
    if (b < 0.0) b += 360.0;
    return b;
}

// =============================================================================
// calculate_junction_com_bearings
// =============================================================================
std::vector<lsdtt_junction_com_bearing> calculate_junction_com_bearings(
    const lsdtt_flow_base&             flow,
    const lsdtt_junction_network_base& network)
{
    std::vector<lsdtt_junction_com_bearing> out;

    const lsdtt_raster_meta& meta = flow.meta();
    const int n_rows = meta.n_rows;
    const int n_cols = meta.n_cols;
    const std::size_t n =
        static_cast<std::size_t>(n_rows) * static_cast<std::size_t>(n_cols);

    // Origin-relative cell-centre coordinates as accumulation inputs (kept small
    // to preserve float precision through the O(N) upslope sum).
    std::vector<float>   vx(n, 0.0f);
    std::vector<float>   vy(n, 0.0f);
    std::vector<uint8_t> mask(n, 0);
    for (int r = 0; r < n_rows; ++r)
    {
        for (int c = 0; c < n_cols; ++c)
        {
            const std::size_t f = static_cast<std::size_t>(r) * n_cols + c;
            if (!flow.is_valid_cell(static_cast<std::int64_t>(f))) continue;
            mask[f] = 1;
            vx[f] = static_cast<float>(cell_centre_x(meta, c) - meta.x_min);
            vy[f] = static_cast<float>(cell_centre_y(meta, r) - meta.y_min);
        }
    }

    const lsdtt_raster_result sx = flow.get_upslope_sum_raster(vx, mask, true);
    const lsdtt_raster_result sy = flow.get_upslope_sum_raster(vy, mask, true);
    const lsdtt_raster_result nc = flow.get_n_contributing_pixels_raster();
    const std::vector<float>& AX = sx.data;
    const std::vector<float>& AY = sy.data;
    const std::vector<float>& NP = nc.data;

    const std::int32_t n_junc = network.n_junctions();
    out.reserve(static_cast<std::size_t>(n_junc));

    for (std::int32_t j = 0; j < n_junc; ++j)
    {
        if (!network.is_valid_junction(j)) continue;

        const std::int64_t node = network.flat_of_junction(j);
        const std::int32_t recv_j = network.receiver_of_junction(j);
        const std::int64_t recv_node =
            (recv_j >= 0) ? network.flat_of_junction(recv_j) : -1;

        // Walk downstream to the pixel just before the downstream junction
        // (or to base level for the outlet junction).
        std::int64_t outlet = node;
        while (true)
        {
            const std::int64_t nxt = flow.receiver_of(outlet);
            if (nxt < 0) break;                              // base level
            if (recv_node >= 0 && nxt == recv_node) break;   // just before downstream junction
            outlet = nxt;
        }

        const std::size_t oi = static_cast<std::size_t>(outlet);
        if (oi >= NP.size() || NP[oi] <= 0.0f) continue;

        const double area  = static_cast<double>(NP[oi]);
        const double com_x = meta.x_min + static_cast<double>(AX[oi]) / area;
        const double com_y = meta.y_min + static_cast<double>(AY[oi]) / area;

        const int orow = flow.row_from_flat_unchecked(outlet);
        const int ocol = flow.col_from_flat_unchecked(outlet);
        const lsdtt_xy oxy = cell_center_from_row_col(meta, orow, ocol);

        lsdtt_junction_com_bearing rec;
        rec.junction_id = j;
        rec.order       = network.stream_order_of_junction(j);
        rec.outlet_x    = oxy.x;
        rec.outlet_y    = oxy.y;
        rec.com_x       = com_x;
        rec.com_y       = com_y;
        rec.bearing_deg = bearing_from_north_degrees(oxy.x, oxy.y, com_x, com_y);
        rec.distance_m  = std::hypot(com_x - oxy.x, com_y - oxy.y);
        rec.n_pixels    = static_cast<std::int64_t>(area);
        out.push_back(rec);
    }

    return out;
}

// =============================================================================
// write_junction_com_bearings
// =============================================================================
void write_junction_com_bearings(
    const std::vector<lsdtt_junction_com_bearing>& records,
    const std::string&                             out_path)
{
    lsdtt_data_loader loader;

    std::vector<std::string> headers =
    {
        "junction_id", "order",
        "outlet_x", "outlet_y", "com_x", "com_y",
        "bearing_deg", "distance_m", "n_pixels"
    };

    std::vector<std::vector<std::string>> rows;
    rows.reserve(records.size());
    for (const auto& r : records)
    {
        rows.push_back({
            std::to_string(r.junction_id),
            std::to_string(r.order),
            precise_string(r.outlet_x),
            precise_string(r.outlet_y),
            precise_string(r.com_x),
            precise_string(r.com_y),
            precise_string(r.bearing_deg),
            precise_string(r.distance_m),
            std::to_string(r.n_pixels)
        });
    }

    loader.write_csv(out_path, headers, rows);
}


} // namespace lsdtt
