#include <cmath>
#include <cstdint>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

#include "lsdtt_raster_resampling_cpu.hpp"

namespace {

int checks = 0;
int failures = 0;

void check(const std::string& name, bool ok)
{
    ++checks;
    if (!ok) ++failures;
    std::cout << (ok ? "PASS" : "FAIL") << "  " << name << '\n';
}

lsdtt::lsdtt_raster_cpu raster(int rows, int cols, double cell_size,
                               std::vector<float> values)
{
    lsdtt::lsdtt_raster_meta meta;
    meta.n_rows = rows;
    meta.n_cols = cols;
    meta.cell_size = cell_size;
    meta.no_data_value = -9999.0;
    return lsdtt::lsdtt_raster_cpu(
        std::move(values), std::vector<std::uint8_t>(rows * cols, 1), meta);
}

lsdtt::lsdtt_raster_cpu empty(int rows, int cols, double cell_size,
                              double x_min = 0.0, double y_min = 0.0)
{
    lsdtt::lsdtt_raster_meta meta;
    meta.n_rows = rows;
    meta.n_cols = cols;
    meta.cell_size = cell_size;
    meta.x_min = x_min;
    meta.y_min = y_min;
    meta.no_data_value = -9999.0;
    return lsdtt::lsdtt_raster_cpu(
        std::vector<float>(rows * cols, -9999.0f),
        std::vector<std::uint8_t>(rows * cols, 0), meta);
}

bool equal(const std::vector<float>& actual, const std::vector<float>& expected,
           double tolerance = 1e-6)
{
    if (actual.size() != expected.size()) return false;
    for (std::size_t i = 0; i < actual.size(); ++i)
        if (std::abs(actual[i] - expected[i]) > tolerance) return false;
    return true;
}

} // namespace

int main()
{
    auto source = raster(4, 4, 1.0,
        {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16});

    auto mean = empty(2, 2, 2.0);
    lsdtt::halve_raster(source, mean, lsdtt::raster_resampling_method::mean);
    check("optimized half mean", equal(mean.get_data(), {3.5f, 5.5f, 11.5f, 13.5f}));

    auto sum = empty(2, 2, 2.0);
    lsdtt::resample_raster(source, sum, lsdtt::raster_resampling_method::sum);
    check("optimized half sum", equal(sum.get_data(), {14, 22, 46, 54}));

    auto minimum = empty(2, 2, 2.0);
    lsdtt::resample_raster(source, minimum, lsdtt::raster_resampling_method::minimum);
    check("minimum", equal(minimum.get_data(), {1, 3, 9, 11}));

    auto maximum = empty(2, 2, 2.0);
    lsdtt::resample_raster(source, maximum, lsdtt::raster_resampling_method::maximum);
    check("maximum", equal(maximum.get_data(), {6, 8, 14, 16}));

    auto doubled = empty(8, 8, 0.5);
    lsdtt::double_raster(source, doubled, lsdtt::raster_resampling_method::nearest);
    check("optimized double nearest",
          doubled.get_data_array()(0, 0) == 1.0f &&
          doubled.get_data_array()(1, 1) == 1.0f &&
          doubled.get_data_array()(7, 7) == 16.0f);

    auto split_sum = empty(8, 8, 0.5);
    lsdtt::double_raster(source, split_sum, lsdtt::raster_resampling_method::sum);
    double split_total = 0.0;
    for (float value : split_sum.get_data()) split_total += value;
    check("double sum is conservative", std::abs(split_total - 136.0) < 1e-6);

    auto thirds = empty(1, 1, 1.25, 0.5, 0.0);
    auto strip = raster(1, 2, 1.0, {0.0f, 10.0f});
    lsdtt::resample_raster(strip, thirds, lsdtt::raster_resampling_method::mean);
    check("unweighted mean", std::abs(thirds.get_data()[0] - 5.0f) < 1e-6);
    lsdtt::resample_raster(strip, thirds, lsdtt::raster_resampling_method::conservative_mean);
    check("overlap-weighted conservative mean",
          std::abs(thirds.get_data()[0] - 6.0f) < 1e-6);

    auto categories = raster(2, 2, 1.0, {2, 1, 2, 1});
    auto category = empty(1, 1, 2.0);
    lsdtt::resample_raster(categories, category, lsdtt::raster_resampling_method::mode);
    check("mode has deterministic low-value tie break", category.get_data()[0] == 1.0f);

    auto square = raster(3, 3, 1.0, {0, 1, 2, 3, 4, 5, 6, 7, 8});
    auto padded = lsdtt::pad_to_nearest_power_of_two(square);
    check("power-of-two padding", padded.n_rows() == 4 && padded.n_cols() == 4 &&
          padded.get_mask_array()(0, 0) == 0 && padded.get_data_array()(1, 0) == 0.0f);

    auto power = lsdtt::interpolate_to_nearest_power_of_two(
        square, lsdtt::power_of_two_direction::higher,
        lsdtt::raster_resampling_method::cubic);
    check("power-of-two interpolation", power.n_rows() == 4 && power.n_cols() == 4 &&
          std::abs(power.get_meta().cell_size - 0.75) < 1e-12);

    std::cout << checks - failures << "/" << checks << " checks passed\n";
    return failures == 0 ? 0 : 1;
}
