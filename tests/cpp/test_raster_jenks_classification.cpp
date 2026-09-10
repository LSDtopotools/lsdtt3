#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "lsdtt_raster_algebra.hpp"

namespace {

int checks = 0;
int failures = 0;

void check(const std::string& name, const bool ok)
{
    ++checks;
    if (!ok) ++failures;
    std::cout << (ok ? "PASS" : "FAIL") << "  " << name << '\n';
}

lsdtt::lsdtt_raster_cpu raster(std::vector<float> data,
                               std::vector<std::uint8_t> mask)
{
    lsdtt::lsdtt_raster_meta meta;
    meta.n_rows = 1;
    meta.n_cols = static_cast<int>(data.size());
    meta.x_min = 123.0;
    meta.y_min = 456.0;
    meta.cell_size = 7.5;
    meta.no_data_value = -9999.0;
    meta.crs_wkt = "test projection";
    return lsdtt::lsdtt_raster_cpu(std::move(data), std::move(mask), meta);
}

bool equal(const std::vector<float>& actual, const std::vector<float>& expected)
{
    return actual == expected;
}

} // namespace

int main()
{
    {
        const auto result = lsdtt::algebra::classify_jenks_categories(
            raster({1, 2, 100, 101}, {1, 1, 1, 1}), 2, 128);
        check("valid classification", equal(result.classified_raster.get_data(), {1, 1, 2, 2}));
        check("valid classified count", result.classified_pixels == 4 && result.jenks.n_values == 4);
        check("valid output class range", *std::min_element(result.classified_raster.get_data().begin(),
              result.classified_raster.get_data().end()) == 1.0f &&
              *std::max_element(result.classified_raster.get_data().begin(),
              result.classified_raster.get_data().end()) == 2.0f);
    }

    {
        const float nan = std::numeric_limits<float>::quiet_NaN();
        const float inf = std::numeric_limits<float>::infinity();
        const auto result = lsdtt::algebra::classify_jenks_categories(
            raster({1, 2, -9999, nan, inf, 10}, {1, 0, 1, 1, 1, 1}), 2, 16);
        check("masked/non-finite/nodata data", equal(result.classified_raster.get_data(),
              {1, -9999, -9999, -9999, -9999, 2}));
        check("masked/non-finite/nodata mask", result.classified_raster.get_mask() ==
              std::vector<std::uint8_t>({1, 0, 0, 0, 0, 1}));
        check("invalid pixels excluded from count", result.classified_pixels == 2 &&
              result.jenks.n_values == 2);
        const auto& meta = result.classified_raster.get_meta();
        check("metadata preserved", meta.n_rows == 1 && meta.n_cols == 6 &&
              meta.x_min == 123.0 && meta.y_min == 456.0 && meta.cell_size == 7.5 &&
              meta.no_data_value == -9999.0 && meta.crs_wkt == "test projection");
    }

    {
        bool threw = false;
        try
        {
            (void)lsdtt::algebra::classify_jenks_categories(
                raster({-9999, std::numeric_limits<float>::quiet_NaN()}, {0, 1}), 3, 8);
        }
        catch (const std::runtime_error& e)
        {
            threw = std::string(e.what()).find("no finite values") != std::string::npos;
        }
        check("empty valid input preserves error", threw);
    }

    {
        const auto result = lsdtt::algebra::classify_jenks_categories(
            raster({7, 7, 7}, {1, 1, 1}), 5, 8);
        check("one unique value", result.jenks.effective_classes == 1 &&
              result.jenks.upper_breaks == std::vector<double>({7.0}) &&
              equal(result.classified_raster.get_data(), {1, 1, 1}));
    }

    {
        const auto result = lsdtt::algebra::classify_jenks_categories(
            raster({1, 2, 2}, {1, 1, 1}), 5, 16);
        check("classes exceeding unique values", result.jenks.requested_classes == 5 &&
              result.jenks.effective_classes == 2 &&
              equal(result.classified_raster.get_data(), {1, 2, 2}));
    }

    {
        const auto result = lsdtt::algebra::classify_jenks_categories(
            raster({0, 1, 2, 3, 4, 5, 6, 7, 8, 9},
                   {1, 1, 1, 1, 1, 1, 1, 1, 1, 1}), 3, 3);
        const auto& data = result.classified_raster.get_data();
        check("histogram path uses capped bins", result.jenks.n_histogram_bins_used == 3 &&
              result.jenks.effective_classes == 3);
        check("histogram output range and count", result.classified_pixels == 10 &&
              *std::min_element(data.begin(), data.end()) == 1.0f &&
              *std::max_element(data.begin(), data.end()) == 3.0f);
    }

    std::cout << checks - failures << "/" << checks << " checks passed\n";
    return failures == 0 ? 0 : 1;
}
