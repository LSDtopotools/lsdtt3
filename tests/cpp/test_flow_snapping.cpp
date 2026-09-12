// Direct tests for the string snapping dispatcher.

#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "lsdtt_flow_snapping_cpu.hpp"
#include "lsdtt_flowinfo_cpu.hpp"
#include "lsdtt_raster_cpu.hpp"

namespace
{

int failures = 0;

void check(const std::string& name, const bool ok)
{
    std::cout << (ok ? "PASS" : "FAIL") << "  " << name << '\n';
    if (!ok) ++failures;
}

lsdtt::lsdtt_flowinfo make_flow()
{
    lsdtt::lsdtt_raster_meta meta;
    meta.n_rows = 7;
    meta.n_cols = 7;
    meta.cell_size = 10.0;
    meta.no_data_value = -9999.0;
    meta.set_default_synthetic_georeference();

    std::vector<float> data(49);
    std::vector<std::uint8_t> mask(49, 1);
    for (int row = 0; row < meta.n_rows; ++row)
    {
        for (int col = 0; col < meta.n_cols; ++col)
        {
            data[static_cast<std::size_t>(row * meta.n_cols + col)] =
                static_cast<float>(100 - row);
        }
    }

    const lsdtt::lsdtt_raster_cpu dem(std::move(data), std::move(mask), meta);
    lsdtt::lsdtt_flowinfo_options options;
    options.boundary_conditions = "oooo";
    options.build_topological_order = true;
    options.build_s_vector = true;
    return lsdtt::lsdtt_flowinfo(
        dem.get_data_array(), dem.get_mask_array(), dem.get_meta(), options);
}

} // namespace

int main()
{
    const lsdtt::lsdtt_flowinfo flow = make_flow();
    lsdtt::lsdtt_square_snap_options square_options;
    square_options.contributing_pixels_threshold = 1;
    square_options.min_search_radius_pixels = 0;
    square_options.max_search_radius_pixels = 2;

    lsdtt::lsdtt_downstream_snap_options downstream_options;
    downstream_options.contributing_pixels_threshold = 1;

    const auto square_direct =
        lsdtt::snap_pixel_to_contributing_pixels_threshold_expanding_square(
            flow, 3, 3, square_options);
    const auto square_dispatched = lsdtt::snap_pixel_by_method(
        flow, 3, 3, "expanding_square", square_options, downstream_options);
    check("expanding_square dispatch matches direct call",
          square_dispatched.status == square_direct.status &&
          square_dispatched.snapped_flat == square_direct.snapped_flat &&
          square_dispatched.search_radius_pixels == square_direct.search_radius_pixels);

    const auto downstream_direct =
        lsdtt::snap_pixel_to_contributing_pixels_threshold_downstream(
            flow, 3, 3, downstream_options);
    const auto downstream_dispatched = lsdtt::snap_pixel_by_method(
        flow, 3, 3, "downstream", square_options, downstream_options);
    check("downstream dispatch matches direct call",
          downstream_dispatched.status == downstream_direct.status &&
          downstream_dispatched.snapped_flat == downstream_direct.snapped_flat &&
          downstream_dispatched.downstream_steps == downstream_direct.downstream_steps);

    const auto none = lsdtt::snap_pixel_by_method(
        flow, 3, 3, "none", square_options, downstream_options);
    check("none returns the valid input pixel",
          none.status == lsdtt::lsdtt_snap_status::success &&
          none.snapped_row == 3 && none.snapped_col == 3 &&
          none.snapped_flat == flow.flat_index(3, 3));

    const auto invalid_none = lsdtt::snap_pixel_by_method(
        flow, -1, 3, "none", square_options, downstream_options);
    check("none validates the input pixel",
          invalid_none.status == lsdtt::lsdtt_snap_status::input_invalid_flow_cell);

    bool unknown_threw = false;
    try
    {
        (void)lsdtt::snap_pixel_by_method(
            flow, 3, 3, "sideways", square_options, downstream_options);
    }
    catch (const std::invalid_argument& e)
    {
        unknown_threw = std::string(e.what()).find("sideways") != std::string::npos;
    }
    check("unknown method throws invalid_argument naming the method", unknown_threw);

    std::cout << "RESULT: " << (failures == 0 ? "PASS" : "FAIL") << '\n';
    return failures == 0 ? 0 : 1;
}
