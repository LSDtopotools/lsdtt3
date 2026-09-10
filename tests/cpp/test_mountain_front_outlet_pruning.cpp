// Direct tests for downstream-most mountain-front outlet pruning.

#include <cstdint>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "lsdtt_flowinfo_cpu.hpp"
#include "lsdtt_mountain_front_cpu.hpp"
#include "lsdtt_raster_cpu.hpp"

namespace
{

int failures = 0;

void check(const std::string& name, const bool ok)
{
    std::cout << (ok ? "PASS" : "FAIL") << "  " << name << '\n';
    if (!ok) ++failures;
}

lsdtt::lsdtt_raster_meta make_meta(const int n_cols)
{
    lsdtt::lsdtt_raster_meta meta;
    meta.n_rows = 1;
    meta.n_cols = n_cols;
    meta.cell_size = 10.0;
    meta.no_data_value = -9999.0;
    meta.set_default_synthetic_georeference();
    return meta;
}

lsdtt::lsdtt_flowinfo make_flow(std::vector<std::int64_t> receivers)
{
    const lsdtt::lsdtt_raster_meta meta =
        make_meta(static_cast<int>(receivers.size()));
    std::vector<float> elevations(receivers.size());
    std::vector<std::uint8_t> mask(receivers.size(), std::uint8_t{1});
    for (std::size_t i = 0; i < elevations.size(); ++i)
    {
        elevations[i] = static_cast<float>(elevations.size() - i);
    }

    const lsdtt::lsdtt_raster_cpu dem(
        std::move(elevations), std::move(mask), meta);
    lsdtt::lsdtt_flowinfo flow(
        dem.get_data_array(), dem.get_mask_array(), dem.get_meta());

    // The pruning routine reads only this array and full_grid_size(). Replacing
    // it allows cycle and branch topologies that valid D8 construction rejects.
    auto& flow_receivers = const_cast<std::vector<std::int64_t>&>(
        flow.receiver_array().get_vector());
    flow_receivers = std::move(receivers);
    return flow;
}

lsdtt::lsdtt_raster_cpu make_outlets(
    std::vector<float> values,
    std::vector<std::uint8_t> mask = {})
{
    const lsdtt::lsdtt_raster_meta meta =
        make_meta(static_cast<int>(values.size()));
    if (mask.empty()) mask.assign(values.size(), std::uint8_t{1});
    return lsdtt::lsdtt_raster_cpu(std::move(values), std::move(mask), meta);
}

bool values_are(const lsdtt::lsdtt_raster_cpu& raster,
                const std::vector<float>& expected)
{
    return raster.get_data() == expected;
}

} // namespace

int main()
{
    using flow_base = lsdtt::lsdtt_flow_base;

    {
        const auto flow = make_flow({1, 2, 3, flow_base::outlet_receiver});
        const auto result = lsdtt::retain_downstream_most_mountain_front_outlets(
            make_outlets({1.0f, 0.0f, 1.0f, 0.0f}), flow);
        check("one chain retains only its downstream-most candidate",
              values_are(result, {0.0f, 0.0f, 1.0f, 0.0f}));
    }

    {
        const auto flow = make_flow({2, 2, 3, flow_base::outlet_receiver});
        const auto result = lsdtt::retain_downstream_most_mountain_front_outlets(
            make_outlets({1.0f, 1.0f, 0.0f, 0.0f}), flow);
        check("candidates on independent branches are retained",
              values_are(result, {1.0f, 1.0f, 0.0f, 0.0f}));
    }

    {
        const auto flow = make_flow({0, flow_base::outlet_receiver});
        const auto result = lsdtt::retain_downstream_most_mountain_front_outlets(
            make_outlets({1.0f, 1.0f}), flow);
        check("self-receiving candidate is pruned and outlet candidate retained",
              values_are(result, {0.0f, 1.0f}));
    }

    {
        const auto flow = make_flow({1, 2, 1});
        const auto result = lsdtt::retain_downstream_most_mountain_front_outlets(
            make_outlets({1.0f, 0.0f, 0.0f}), flow);
        check("receiver cycle is bounded and retains an unshadowed candidate",
              values_are(result, {1.0f, 0.0f, 0.0f}));
    }

    {
        const float nan = std::numeric_limits<float>::quiet_NaN();
        const auto flow = make_flow({1, 2, 3, 4, 5, flow_base::outlet_receiver});
        const auto result = lsdtt::retain_downstream_most_mountain_front_outlets(
            make_outlets({1.0f, 1.0f, nan, -9999.0f, 0.5f, 0.5001f},
                         {1, 0, 1, 1, 1, 1}),
            flow);
        check("masked, nodata, non-finite, and threshold inputs are handled",
              result.get_mask() ==
                  std::vector<std::uint8_t>({1, 0, 1, 1, 1, 1}) &&
              values_are(result,
                         {1.0f, -9999.0f, 0.0f, 0.0f, 0.0f, 1.0f}));
    }

    {
        bool threw = false;
        try
        {
            const auto flow = make_flow(
                {1, 2, flow_base::outlet_receiver});
            (void)lsdtt::retain_downstream_most_mountain_front_outlets(
                make_outlets({1.0f, 0.0f, 0.0f, 0.0f}), flow);
        }
        catch (const std::runtime_error& error)
        {
            threw = std::string(error.what()).find("grid size") !=
                    std::string::npos;
        }
        check("incompatible flow and outlet dimensions throw", threw);
    }

    std::cout << "RESULT: " << (failures == 0 ? "PASS" : "FAIL") << '\n';
    return failures == 0 ? 0 : 1;
}
