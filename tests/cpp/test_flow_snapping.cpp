// Direct tests for the string snapping dispatcher, for area-matched snapping
// (snap_pixel_to_reported_drainage_area_expanding_square) on a synthetic
// tributary junction, and for largest-area snapping
// (snap_pixel_to_contributing_pixels_threshold_largest_area).

#include <algorithm>
#include <cmath>
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

// A trunk running south down column 10 and a tributary running east along
// row 12 into it, on a 30 x 21 grid at 10 m. Hillslopes fall 3 m per pixel
// towards the nearer channel. The tributary meets the trunk at (13, 10), so a
// point at (13, 9) is one pixel from the tributary cell (12, 9) and one pixel
// from the trunk cell (13, 10) below the confluence.
lsdtt::lsdtt_flowinfo make_junction_flow()
{
    lsdtt::lsdtt_raster_meta meta;
    meta.n_rows = 30;
    meta.n_cols = 21;
    meta.cell_size = 10.0;
    meta.no_data_value = -9999.0;
    meta.set_default_synthetic_georeference();

    std::vector<float> data(static_cast<std::size_t>(meta.n_rows * meta.n_cols));
    std::vector<std::uint8_t> mask(data.size(), 1);
    for (int row = 0; row < meta.n_rows; ++row)
    {
        for (int col = 0; col < meta.n_cols; ++col)
        {
            const double trunk = 100.0 - row + 3.0 * std::abs(col - 10);
            double z = trunk;
            if (col < 10)
            {
                const double trib = 88.0 + (10 - col) + 3.0 * std::abs(row - 12);
                z = std::min(z, trib);
            }
            data[static_cast<std::size_t>(row * meta.n_cols + col)] = static_cast<float>(z);
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

void check_reported_area_snapping()
{
    const lsdtt::lsdtt_flowinfo flow = make_junction_flow();
    const double cell_area = 100.0;
    auto area_of = [&](int r, int c) {
        return static_cast<double>(flow.s_vector_size(flow.flat_index(r, c))) * cell_area;
    };
    const double a_trib = area_of(12, 9);
    const double a_trunk = area_of(13, 10);
    std::cout << "junction fixture: tributary (12,9) " << a_trib / cell_area
              << " px, trunk (13,10) " << a_trunk / cell_area << " px, trunk above (11,10) "
              << area_of(11, 10) / cell_area << " px, (12,10) " << area_of(12, 10) / cell_area
              << " px, (12,8) " << area_of(12, 8) / cell_area << " px\n";
    check("T-SN fixture: the trunk below the confluence is more than 1.25x the tributary",
          a_trunk > 1.25 * a_trib && a_trib >= 20.0 * cell_area);

    lsdtt::lsdtt_square_snap_options sq;
    sq.contributing_pixels_threshold = 20;
    sq.min_search_radius_pixels = 0;
    sq.max_search_radius_pixels = 3;
    const auto nearest =
        lsdtt::snap_pixel_to_contributing_pixels_threshold_expanding_square(flow, 13, 9, sq);
    check("T-SN1: expanding_square from (13,9) takes the trunk (13,10) (the failure "
          "area matching fixes)",
          nearest.status == lsdtt::lsdtt_snap_status::success &&
          nearest.snapped_row == 13 && nearest.snapped_col == 10);

    lsdtt::lsdtt_reported_area_snap_options ao;
    ao.contributing_pixels_threshold = 20;
    ao.max_search_radius_pixels = 3;
    ao.area_ratio_tolerance = 1.25;
    const auto matched =
        lsdtt::snap_pixel_to_reported_drainage_area_expanding_square(flow, 13, 9, a_trib, ao);
    check("T-SN1: reported_area with the tributary's area takes the tributary cell (12,9)",
          matched.status == lsdtt::lsdtt_snap_status::success &&
          matched.snapped_row == 12 && matched.snapped_col == 9 &&
          matched.snapped_flat == flow.flat_index(12, 9) && matched.search_radius_pixels == 1 &&
          matched.downstream_steps == -1);

    // T-SN2: nothing within the tolerance: the cell closest in area (the
    // largest in the window, the trunk at the window's bottom row) and still
    // success; the caller sees the ratio is outside the tolerance.
    const double huge = 50.0 * a_trunk;
    const auto closest =
        lsdtt::snap_pixel_to_reported_drainage_area_expanding_square(flow, 13, 9, huge, ao);
    const double ratio = closest.status == lsdtt::lsdtt_snap_status::success
                             ? static_cast<double>(closest.snapped_contributing_pixels) * cell_area / huge
                             : 0.0;
    check("T-SN2: no cell within the tolerance: success on the min |ln ratio| cell (16,10)",
          closest.status == lsdtt::lsdtt_snap_status::success && closest.snapped_row == 16 &&
          closest.snapped_col == 10 && closest.search_radius_pixels == 3);
    check("T-SN2: the snapped area ratio is outside the tolerance for the caller to flag",
          ratio < 1.0 / 1.25);

    // T-SN3: no channel cell in the window.
    lsdtt::lsdtt_reported_area_snap_options far = ao;
    far.max_search_radius_pixels = 1;
    const auto none_found =
        lsdtt::snap_pixel_to_reported_drainage_area_expanding_square(flow, 2, 2, a_trib, far);
    check("T-SN3: no channel cell within R -> threshold_not_found_in_search_window",
          none_found.status == lsdtt::lsdtt_snap_status::threshold_not_found_in_search_window &&
          none_found.snapped_flat == -1);
    const auto off_grid =
        lsdtt::snap_pixel_to_reported_drainage_area_expanding_square(flow, -5, 2, a_trib, ao);
    check("T-SN3: a point off the grid -> input_out_of_bounds",
          off_grid.status == lsdtt::lsdtt_snap_status::input_out_of_bounds);

    // Ties: from (11,9) the tributary (12,9) and the trunk (11,10) are both one
    // pixel away and both within a tolerance of 2 of 150 and of 110 px; the
    // smaller |ln ratio| decides, so the choice flips with the reported area.
    lsdtt::lsdtt_reported_area_snap_options wide = ao;
    wide.max_search_radius_pixels = 1;
    wide.area_ratio_tolerance = 2.0;
    const auto tie_trunk = lsdtt::snap_pixel_to_reported_drainage_area_expanding_square(
        flow, 11, 9, 150.0 * cell_area, wide);
    const auto tie_trib = lsdtt::snap_pixel_to_reported_drainage_area_expanding_square(
        flow, 11, 9, 110.0 * cell_area, wide);
    check("T-SN: equally near in-tolerance cells go to the smaller |ln ratio|",
          flow.s_vector_size(flow.flat_index(11, 9)) < 20 &&
          tie_trunk.status == lsdtt::lsdtt_snap_status::success && tie_trunk.snapped_row == 11 &&
          tie_trunk.snapped_col == 10 && tie_trib.status == lsdtt::lsdtt_snap_status::success &&
          tie_trib.snapped_row == 12 && tie_trib.snapped_col == 9);

    auto throws = [&](double area, double tol, std::uint64_t thr) {
        lsdtt::lsdtt_reported_area_snap_options o = ao;
        o.area_ratio_tolerance = tol;
        o.contributing_pixels_threshold = thr;
        try
        {
            (void)lsdtt::snap_pixel_to_reported_drainage_area_expanding_square(flow, 13, 9, area, o);
        }
        catch (const std::invalid_argument&)
        {
            return true;
        }
        return false;
    };
    check("T-SN: tolerance <= 1, a non-positive area and a zero threshold throw",
          throws(a_trib, 1.0, 20) && throws(0.0, 1.25, 20) && throws(std::nan(""), 1.25, 20) &&
          throws(a_trib, 1.25, 0));

    // T-SN4: the dispatcher.
    lsdtt::lsdtt_downstream_snap_options ds;
    ds.contributing_pixels_threshold = 20;
    const auto dispatched =
        lsdtt::snap_pixel_by_method(flow, 13, 9, "reported_area", sq, ds, a_trib, ao);
    check("T-SN4: reported_area dispatch matches the direct call",
          dispatched.status == matched.status && dispatched.snapped_flat == matched.snapped_flat &&
          dispatched.search_radius_pixels == matched.search_radius_pixels);
    const auto sq_dispatched =
        lsdtt::snap_pixel_by_method(flow, 13, 9, "expanding_square", sq, ds, a_trib, ao);
    check("T-SN4: the area overload passes the other methods through unchanged",
          sq_dispatched.snapped_flat == nearest.snapped_flat &&
          sq_dispatched.search_radius_pixels == nearest.search_radius_pixels);
    bool unknown = false, needs_area = false;
    try { (void)lsdtt::snap_pixel_by_method(flow, 13, 9, "sideways", sq, ds, a_trib, ao); }
    catch (const std::invalid_argument& e)
    {
        unknown = std::string(e.what()).find("reported_area") != std::string::npos;
    }
    try { (void)lsdtt::snap_pixel_by_method(flow, 13, 9, "reported_area", sq, ds); }
    catch (const std::invalid_argument& e)
    {
        needs_area = std::string(e.what()).find("reported drainage area") != std::string::npos;
    }
    check("T-SN4: an unknown method still throws (listing reported_area)", unknown);
    check("T-SN4: reported_area without an area throws", needs_area);
}

// T-LA: largest-area snapping on the junction fixture and on the south-tilted
// plane of make_flow, where every pixel drains due south so a pixel in row r
// has r + 1 contributing pixels.
void check_largest_area_snapping()
{
    const lsdtt::lsdtt_flowinfo flow = make_junction_flow();
    auto count_of = [&](int r, int c) { return flow.s_vector_size(flow.flat_index(r, c)); };

    // Brute-force maximum over the window, for comparison.
    std::uint64_t window_max = 0;
    for (int r = 10; r <= 16; ++r)
        for (int c = 6; c <= 12; ++c) window_max = std::max(window_max, count_of(r, c));

    lsdtt::lsdtt_largest_area_snap_options lo;
    lo.contributing_pixels_threshold = 20;
    lo.max_search_radius_pixels = 3;
    const auto largest =
        lsdtt::snap_pixel_to_contributing_pixels_threshold_largest_area(flow, 13, 9, lo);
    check("T-LA1: from (13,9) with R=3 the trunk at the window's bottom row (16,10) is taken",
          largest.status == lsdtt::lsdtt_snap_status::success && largest.snapped_row == 16 &&
          largest.snapped_col == 10 && largest.snapped_flat == flow.flat_index(16, 10) &&
          largest.snapped_contributing_pixels == window_max &&
          largest.search_radius_pixels == 3 && largest.downstream_steps == -1);

    // The threshold only accepts or rejects the largest cell.
    lsdtt::lsdtt_largest_area_snap_options high = lo;
    high.contributing_pixels_threshold = window_max;
    const auto at_threshold =
        lsdtt::snap_pixel_to_contributing_pixels_threshold_largest_area(flow, 13, 9, high);
    high.contributing_pixels_threshold = window_max + 1;
    const auto above_threshold =
        lsdtt::snap_pixel_to_contributing_pixels_threshold_largest_area(flow, 13, 9, high);
    check("T-LA2: a threshold equal to the largest count succeeds (>=), one above fails",
          at_threshold.status == lsdtt::lsdtt_snap_status::success &&
          at_threshold.snapped_flat == largest.snapped_flat &&
          above_threshold.status ==
              lsdtt::lsdtt_snap_status::threshold_not_found_in_search_window &&
          above_threshold.snapped_flat == -1);

    lsdtt::lsdtt_largest_area_snap_options own = lo;
    own.max_search_radius_pixels = 0;
    own.contributing_pixels_threshold = 1;
    const auto own_pixel =
        lsdtt::snap_pixel_to_contributing_pixels_threshold_largest_area(flow, 13, 9, own);
    check("T-LA3: R=0 with threshold 1 keeps the input pixel",
          own_pixel.status == lsdtt::lsdtt_snap_status::success && own_pixel.snapped_row == 13 &&
          own_pixel.snapped_col == 9 && own_pixel.search_radius_pixels == 0);

    const auto off_grid =
        lsdtt::snap_pixel_to_contributing_pixels_threshold_largest_area(flow, -5, 2, lo);
    check("T-LA4: a point off the grid -> input_out_of_bounds",
          off_grid.status == lsdtt::lsdtt_snap_status::input_out_of_bounds);

    auto throws = [&](std::uint64_t thr, int radius) {
        lsdtt::lsdtt_largest_area_snap_options o;
        o.contributing_pixels_threshold = thr;
        o.max_search_radius_pixels = radius;
        try
        {
            (void)lsdtt::snap_pixel_to_contributing_pixels_threshold_largest_area(flow, 13, 9, o);
        }
        catch (const std::invalid_argument&)
        {
            return true;
        }
        return false;
    };
    check("T-LA5: a zero threshold and a negative radius throw", throws(0, 3) && throws(20, -1));

    // Ties on count go to the nearer pixel: on the plane every bottom-row
    // pixel has 7 contributing pixels, and (6,3) is the nearest to (4,3).
    const lsdtt::lsdtt_flowinfo plane = make_flow();
    lsdtt::lsdtt_largest_area_snap_options po;
    po.contributing_pixels_threshold = 1;
    po.max_search_radius_pixels = 3;
    const auto tie = lsdtt::snap_pixel_to_contributing_pixels_threshold_largest_area(plane, 4, 3, po);
    check("T-LA6: equal counts go to the nearer pixel ((6,3) from (4,3))",
          plane.s_vector_size(plane.flat_index(6, 1)) == 7 &&
          plane.s_vector_size(plane.flat_index(6, 3)) == 7 &&
          tie.status == lsdtt::lsdtt_snap_status::success && tie.snapped_row == 6 &&
          tie.snapped_col == 3 && tie.search_radius_pixels == 2);

    // Dispatch: largest_area takes the threshold and max radius from the
    // square options and ignores min_search_radius_pixels.
    lsdtt::lsdtt_square_snap_options sq;
    sq.contributing_pixels_threshold = lo.contributing_pixels_threshold;
    sq.min_search_radius_pixels = 2;
    sq.max_search_radius_pixels = lo.max_search_radius_pixels;
    lsdtt::lsdtt_downstream_snap_options ds;
    ds.contributing_pixels_threshold = 20;
    const auto dispatched = lsdtt::snap_pixel_by_method(flow, 13, 9, "largest_area", sq, ds);
    check("T-LA7: largest_area dispatch matches the direct call",
          dispatched.status == largest.status && dispatched.snapped_flat == largest.snapped_flat &&
          dispatched.search_radius_pixels == largest.search_radius_pixels);
    lsdtt::lsdtt_reported_area_snap_options ao;
    const auto area_overload =
        lsdtt::snap_pixel_by_method(flow, 13, 9, "largest_area", sq, ds, 1.0, ao);
    check("T-LA7: the area overload passes largest_area through unchanged",
          area_overload.status == largest.status &&
          area_overload.snapped_flat == largest.snapped_flat);
    bool listed = false;
    try { (void)lsdtt::snap_pixel_by_method(flow, 13, 9, "sideways", sq, ds); }
    catch (const std::invalid_argument& e)
    {
        listed = std::string(e.what()).find("largest_area") != std::string::npos;
    }
    check("T-LA7: the unknown-method message lists largest_area", listed);
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

    check_reported_area_snapping();
    check_largest_area_snapping();

    std::cout << "RESULT: " << (failures == 0 ? "PASS" : "FAIL") << '\n';
    return failures == 0 ? 0 : 1;
}
