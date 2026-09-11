// Focused contract tests for algebra::raster_transform (Phase 3 C11).

#include <cmath>
#include <cstring>
#include <iostream>
#include <limits>
#include <string>

#include "../../src/cpu/lsdtt_crn_units.hpp"
#include "../../src/cpu/lsdtt_raster_algebra.hpp"

namespace
{
int failures = 0;

void check(bool condition, const std::string& message)
{
    if (!condition)
    {
        ++failures;
        std::cerr << "FAIL: " << message << '\n';
    }
}

bool same_float(float a, float b)
{
    return std::memcmp(&a, &b, sizeof(float)) == 0;
}

bool same_meta(const lsdtt::lsdtt_raster_meta& a,
               const lsdtt::lsdtt_raster_meta& b)
{
    return a.n_rows == b.n_rows && a.n_cols == b.n_cols
        && a.x_min == b.x_min && a.y_min == b.y_min
        && a.cell_size == b.cell_size
        && a.no_data_value == b.no_data_value
        && a.crs_wkt == b.crs_wkt && a.epsg_code == b.epsg_code
        && a.data_type == b.data_type;
}

lsdtt::lsdtt_raster_cpu old_mass_flux_loop(
    const lsdtt::lsdtt_raster_cpu& source)
{
    const auto& meta = source.get_meta();
    lsdtt::lsdtt_array2d_f data(meta.n_rows, meta.n_cols,
                                static_cast<float>(meta.no_data_value));
    lsdtt::lsdtt_array2d_u8 mask(meta.n_rows, meta.n_cols, 0u);
    for (int row = 0; row < meta.n_rows; ++row)
        for (int col = 0; col < meta.n_cols; ++col)
            if (source.get_mask_array()(row, col))
            {
                data(row, col) = static_cast<float>(
                    lsdtt::eff_erate_to_kg_m2_yr(static_cast<double>(
                        source.get_data_array()(row, col))));
                mask(row, col) = 1u;
            }
    return lsdtt::lsdtt_raster_cpu(std::move(data), std::move(mask), meta);
}
}  // namespace

int main()
{
    lsdtt::lsdtt_raster_meta meta;
    meta.n_rows = 2;
    meta.n_cols = 4;
    meta.x_min = 123.5;
    meta.y_min = -456.25;
    meta.cell_size = 30.0;
    meta.no_data_value = -9999.0;
    meta.crs_wkt = "fixture-wkt";
    meta.epsg_code = 32630;
    meta.data_type = GDT_Float64;

    lsdtt::lsdtt_array2d_f data(2, 4, 0.0f);
    lsdtt::lsdtt_array2d_u8 mask(2, 4, 1u);
    data(0, 0) = 0.001f;
    data(0, 1) = static_cast<float>(meta.no_data_value); // valid nodata payload
    data(0, 2) = std::numeric_limits<float>::quiet_NaN();
    data(0, 3) = std::numeric_limits<float>::infinity();
    data(1, 0) = -std::numeric_limits<float>::infinity();
    data(1, 1) = 123456.0f; // masked garbage
    mask(1, 1) = 0u;
    data(1, 2) = std::numeric_limits<float>::quiet_NaN(); // masked garbage
    mask(1, 2) = 0u;
    data(1, 3) = 1.2345678f;
    lsdtt::lsdtt_raster_cpu source(std::move(data), std::move(mask), meta);

    int callback_calls = 0;
    const auto custom = lsdtt::algebra::raster_transform(
        source, [&callback_calls](float value) {
            ++callback_calls;
            return value + 0.25f;
        });
    check(callback_calls == 6, "callback runs once per valid cell only");
    check(same_float(custom.get_data_array()(0, 0), 0.251f),
          "callback output is retained");
    check(same_meta(custom.get_meta(), meta), "metadata is preserved");
    check(custom.get_mask() == source.get_mask(), "valid-cell mask is preserved");
    check(same_float(custom.get_data_array()(1, 1), -9999.0f)
          && same_float(custom.get_data_array()(1, 2), -9999.0f),
          "masked garbage becomes canonical nodata");
    check(same_float(custom.get_data_array()(0, 1), -9998.75f),
          "valid nodata-valued payload reaches callback");
    check(std::isnan(custom.get_data_array()(0, 2))
          && std::isinf(custom.get_data_array()(0, 3))
          && std::isinf(custom.get_data_array()(1, 0)),
          "valid non-finite payloads reach callback");

    const auto converted = lsdtt::algebra::raster_transform(
        source, [](float value) {
            return static_cast<float>(lsdtt::eff_erate_to_kg_m2_yr(
                static_cast<double>(value)));
        });
    const auto old = old_mass_flux_loop(source);
    check(converted.get_mask() == old.get_mask(),
          "conversion mask exactly matches old scalar loop");
    for (std::size_t i = 0; i < converted.get_data().size(); ++i)
        check(same_float(converted.get_data()[i], old.get_data()[i]),
              "conversion data bitwise matches old scalar loop at cell "
              + std::to_string(i));

    int empty_calls = 0;
    const lsdtt::lsdtt_raster_cpu empty;
    const auto empty_out = lsdtt::algebra::raster_transform(
        empty, [&empty_calls](float value) { ++empty_calls; return value; });
    check(empty_out.empty() && empty_out.get_mask().empty(),
          "empty raster produces empty raster");
    check(empty_calls == 0, "empty raster does not invoke callback");
    check(same_meta(empty_out.get_meta(), empty.get_meta()),
          "empty raster metadata is preserved");

    if (failures != 0)
    {
        std::cerr << failures << " raster_transform check(s) failed\n";
        return 1;
    }
    std::cout << "PASS: raster_transform contract and old-loop agreement\n";
    return 0;
}
