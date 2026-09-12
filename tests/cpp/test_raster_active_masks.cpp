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
                               std::vector<std::uint8_t> mask,
                               int n_rows = 1)
{
    lsdtt::lsdtt_raster_meta meta;
    meta.n_rows = n_rows;
    meta.n_cols = static_cast<int>(data.size()) / n_rows;
    meta.x_min = 123.0;
    meta.y_min = 456.0;
    meta.cell_size = 7.5;
    meta.no_data_value = -9999.0;
    meta.crs_wkt = "test projection";
    return lsdtt::lsdtt_raster_cpu(std::move(data), std::move(mask), meta);
}

bool metadata_preserved(const lsdtt::lsdtt_raster_cpu& value)
{
    const auto& meta = value.get_meta();
    return meta.n_rows == 1 && meta.x_min == 123.0 && meta.y_min == 456.0 &&
           meta.cell_size == 7.5 && meta.no_data_value == -9999.0 &&
           meta.crs_wkt == "test projection";
}

} // namespace

int main()
{
    const float nan = std::numeric_limits<float>::quiet_NaN();
    const float inf = std::numeric_limits<float>::infinity();

    {
        const auto result = lsdtt::algebra::raster_mask_from_rounded_value_set(
            raster({-2.5f, -1.5f, -0.5f, 0.49f, 0.5f, 1.49f, 1.5f,
                    2.6f, nan, inf, -9999.0f, 4.0f},
                   {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0}),
            {-3, -2, -1, 1, 2, 3, 3});
        check("value-set lround truth table",
              result.raster.get_data() == std::vector<float>(
                  {1, 1, 1, 0, 1, 1, 1, 1, 0, 0, 0, 0}));
        check("value-set preserves source mask",
              result.raster.get_mask() == std::vector<std::uint8_t>(
                  {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0}));
        check("value-set duplicates do not alter count", result.selected_pixels == 7);
        check("value-set metadata preserved", metadata_preserved(result.raster));
    }

    {
        const auto empty = lsdtt::algebra::raster_mask_from_rounded_value_set(
            raster({-9999.0f, 2.0f}, {1, 0}), {});
        check("empty set writes zero under valid and masked pixels",
              empty.raster.get_data() == std::vector<float>({0, 0}) &&
              empty.raster.get_mask() == std::vector<std::uint8_t>({1, 0}) &&
              empty.selected_pixels == 0);

        const auto nodata_as_value =
            lsdtt::algebra::raster_mask_from_rounded_value_set(
                raster({-9999.0f}, {1}), {-9999});
        check("source mask, not nodata sentinel, controls membership",
              nodata_as_value.raster.get_data() == std::vector<float>({1}) &&
              nodata_as_value.selected_pixels == 1);
    }

    {
        const auto domain = raster(std::vector<float>(8, -9999.0f),
                                   std::vector<std::uint8_t>(8, 0));
        const auto active = raster({0, 2, -3, 1, 5, nan, 0, inf},
                                   {1, 1, 1, 0, 1, 1, 0, 1});
        const auto ice = raster({0, 0, 2, 4, -7, 8, nan, inf},
                                {1, 1, 1, 0, 1, 0, 1, 1});
        const auto result = lsdtt::algebra::raster_active_mask_and_not_excluded(
            domain, &active, ice);
        check("active AND-NOT-excluded truth table",
              result.raster.get_data() ==
                  std::vector<float>({0, 1, 0, 0, 0, 1, 0, 0}));
        check("AND-NOT output mask is all valid",
              result.raster.get_mask() == std::vector<std::uint8_t>(8, 1));
        check("AND-NOT counts active exclusions only",
              result.active_pixels == 2 && result.excluded_pixels == 3);
        check("AND-NOT uses domain metadata despite invalid domain",
              metadata_preserved(result.raster));
    }

    {
        const auto domain = raster(std::vector<float>(6, 0),
                                   std::vector<std::uint8_t>(6, 1));
        const auto ice = raster({0, -9999, 2, -4, nan, inf},
                                {1, 1, 1, 1, 1, 0});
        const auto result = lsdtt::algebra::raster_active_mask_and_not_excluded(
            domain, nullptr, ice);
        check("null active and ice nodata/non-finite truth table",
              result.raster.get_data() == std::vector<float>({1, 1, 0, 0, 0, 1}));
        check("null active counts", result.active_pixels == 3 &&
              result.excluded_pixels == 3);
    }

    {
        const auto domain = raster({0, 0}, {1, 1});
        const auto short_raster = raster({0}, {1});
        bool ice_threw = false;
        bool active_threw = false;
        try
        {
            (void)lsdtt::algebra::raster_active_mask_and_not_excluded(
                domain, nullptr, short_raster);
        }
        catch (const std::invalid_argument&) { ice_threw = true; }
        try
        {
            (void)lsdtt::algebra::raster_active_mask_and_not_excluded(
                domain, &short_raster, domain);
        }
        catch (const std::invalid_argument&) { active_threw = true; }
        check("exclusion dimension mismatch throws", ice_threw);
        check("active dimension mismatch throws", active_threw);
    }

    std::cout << checks - failures << "/" << checks << " checks passed\n";
    return failures == 0 ? 0 : 1;
}
