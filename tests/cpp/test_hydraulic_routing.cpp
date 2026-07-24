// test_hydraulic_routing.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)

#include <cmath>
#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

#include "lsdtt_array2d.hpp"
#include "lsdtt_flowgraph_cpu.hpp"
#include "lsdtt_hydraulic_graph_cpu.hpp"

namespace
{

int g_n_checks = 0;
int g_n_failed = 0;

void check(const std::string& name, bool ok, const std::string& detail = "")
{
    ++g_n_checks;
    if (!ok) ++g_n_failed;
    std::cout << (ok ? "PASS" : "FAIL") << "  " << name;
    if (!ok && !detail.empty()) std::cout << "   [" << detail << "]";
    std::cout << "\n";
}

lsdtt::lsdtt_raster_meta make_meta_3x3()
{
    lsdtt::lsdtt_raster_meta meta;
    meta.n_rows = 3;
    meta.n_cols = 3;
    meta.cell_size = 1.0;
    meta.no_data_value = -9999.0;
    meta.set_default_synthetic_georeference();
    return meta;
}

lsdtt::lsdtt_flowgraph_cpu make_flowgraph_centre_pit()
{
    const auto meta = make_meta_3x3();

    lsdtt::lsdtt_array2d_f data(3, 3, 10.0f);
    data(1, 1) = 0.0f;

    lsdtt::lsdtt_array2d_u8 mask(3, 3, 1);

    lsdtt::lsdtt_flowgraph_options opts;
    opts.boundary_conditions = "oooo";
    opts.method = lsdtt::lsdtt_flowgraph_method::mfd;
    opts.mfd_exponent = 1.0f;
    opts.n_threads = 1;

    return lsdtt::lsdtt_flowgraph_cpu(data, mask, meta, opts);
}

double run_one_step_centre_qin(lsdtt_hydraulic_run_mode mode)
{
    const auto fg = make_flowgraph_centre_pit();

    lsdtt::lsdtt_array2d_f bed(3, 3, 10.0f);
    bed(1, 1) = 0.0f;

    lsdtt::lsdtt_hydraulic_options options;
    options.dt = 0.1;
    options.mannings_n = 0.033;
    options.min_slope = 1.0e-8;
    options.uniform_rainfall_rate = 1.0e-3;
    options.receiver_width = -1.0;
    options.max_dh_rate = 1.0;
    options.print_startup_message = false;

    lsdtt::lsdtt_hydraulic_graph_cpu solver(fg, bed, options);
    solver.set_run_mode(mode);
    solver.set_reroute_every(1);

    std::vector<double> water_depth(9, 0.0);
    water_depth[4] = 50.0;
    solver.set_water_depth(water_depth);

    solver.run_one_step();
    const auto qin = solver.get_qin_raster();
    return static_cast<double>(qin.data[4]);
}

void test_transient_dynamic_dispatch_changes_qin()
{
    const double q_fixed_centre = run_one_step_centre_qin(
        lsdtt_hydraulic_run_mode::transient_fixed_graph);
    const double q_dynamic_centre = run_one_step_centre_qin(
        lsdtt_hydraulic_run_mode::transient_dynamic_graph);

    const bool changed =
        std::abs(q_fixed_centre - q_dynamic_centre) > 1.0e-9 &&
        q_dynamic_centre < q_fixed_centre;

    check("H1 transient dynamic reroutes qin", changed,
          "expected centre qin to differ and decrease under dynamic rerouting");
}

} // namespace

int main()
{
    test_transient_dynamic_dispatch_changes_qin();

    if (g_n_failed == 0)
    {
        std::cout << "RESULT: PASS (" << g_n_checks << " checks)\n";
        return 0;
    }

    std::cout << "RESULT: FAIL (" << g_n_failed << " / " << g_n_checks
              << " checks failed)\n";
    return 1;
}
