// test_cellular_hillslope_solver.cpp
//
// Correctness tests for the global deterministic cellular hillslope oracle.
// Newton-Roering is the only numerical reference.

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

#include "../../src/cpu/lsdtt_hillslope_cellular_cpu.hpp"
#include "../../src/cpu/lsdtt_hillslope_newton_cpu.hpp"
#include "../../src/cpu/lsdtt_lem_cpu.hpp"
#include "../../src/cpu/lsdtt_raster_synthetic_surfaces_cpu.hpp"
#include "../../src/cpu/lsdtt_synthetic_network_cpu.hpp"

namespace
{

constexpr double PI = 3.14159265358979323846;
constexpr double DX = 5.0;
constexpr double D = 0.1;
constexpr double U = 4.0e-4;
constexpr double SC = 0.8;

int g_checks = 0;
int g_failures = 0;

void check(bool ok, const std::string& what)
{
    ++g_checks;
    if (ok) std::cout << "  ok  : " << what << "\n";
    else
    {
        ++g_failures;
        std::cout << "  FAIL: " << what << "\n";
    }
}

lsdtt::lsdtt_raster_meta meta(int nr, int nc, double dx = DX)
{
    lsdtt::lsdtt_raster_meta m;
    m.n_rows = nr;
    m.n_cols = nc;
    m.cell_size = dx;
    return m;
}

double linf(const lsdtt::lsdtt_raster_cpu& a, const lsdtt::lsdtt_raster_cpu& b)
{
    double out = 0.0;
    for (std::size_t i = 0; i < a.get_data().size(); ++i)
        if (a.get_mask()[i] && b.get_mask()[i])
            out = std::max(out, std::fabs(
                static_cast<double>(a.get_data()[i]) - b.get_data()[i]));
    return out;
}

double elevation_sum(const lsdtt::lsdtt_raster_cpu& dem)
{
    double out = 0.0;
    for (std::size_t i = 0; i < dem.get_data().size(); ++i)
        if (dem.get_mask()[i]) out += dem.get_data()[i];
    return out;
}

std::vector<std::uint8_t> strip_channels(
    const lsdtt::lsdtt_raster_cpu& channel_raster)
{
    std::vector<std::uint8_t> out(channel_raster.get_data().size(), 0u);
    for (std::size_t i = 0; i < out.size(); ++i)
        out[i] = channel_raster.get_data()[i] != 0.0f ? 1u : 0u;
    return out;
}

std::vector<std::uint8_t> edge_channels(int nr, int nc)
{
    std::vector<std::uint8_t> out(static_cast<std::size_t>(nr) * nc, 0u);
    for (int r = 0; r < nr; ++r)
        for (int c = 0; c < nc; ++c)
            if (r == 0 || c == 0 || r == nr - 1 || c == nc - 1)
                out[static_cast<std::size_t>(r) * nc + c] = 1u;
    return out;
}

lsdtt::lsdtt_hillslope_cellular_options cellular_options(
    const std::string& bc = "nnnn")
{
    lsdtt::lsdtt_hillslope_cellular_options o;
    o.D = D;
    o.U = 0.0;
    o.Sc = SC;
    o.r_match = 0.85;
    o.cfl_safety = 0.8;
    o.boundary_conditions = bc;
    o.n_threads = 1;
    return o;
}

lsdtt::lsdtt_hillslope_newton_options newton_options(
    const std::string& bc = "nnnn")
{
    lsdtt::lsdtt_hillslope_newton_options o;
    o.D = static_cast<float>(D);
    o.U = 0.0f;
    o.Sc = static_cast<float>(SC);
    o.nonlinear = true;
    o.boundary_conditions = bc;
    o.residual_tol = 1e-11;
    o.max_iterations = 60;
    o.solver_choice = 0;
    o.n_threads = 1;
    return o;
}

lsdtt::lsdtt_raster_cpu periodic_dem(int nr, int nc)
{
    std::vector<float> z(static_cast<std::size_t>(nr) * nc, 0.0f);
    std::vector<std::uint8_t> mask(z.size(), 1u);
    for (int r = 0; r < nr; ++r)
        for (int c = 0; c < nc; ++c)
            z[static_cast<std::size_t>(r) * nc + c] = static_cast<float>(
                0.8 * std::sin(2.0 * PI * r / nr)
                + 0.5 * std::cos(2.0 * PI * c / nc));
    return lsdtt::lsdtt_raster_cpu(std::move(z), std::move(mask), meta(nr, nc));
}

} // namespace

int main()
{
    using namespace lsdtt;
    std::cout << "\n=== Cellular hillslope global oracle ===\n";

    // Probe Newton once so the supported Eigen-free build skips comparisons.
    const lsdtt_raster_cpu probe_dem = periodic_dem(5, 5);
    const auto probe = compute_hillslope_newton(
        probe_dem, edge_channels(5, 5), newton_options());
    if (!probe.had_eigen)
    {
        std::cout << "SKIP: Newton reference requires Eigen\n";
        return 0;
    }

    // C1: every internal face is equal and opposite on a closed domain.
    {
        const lsdtt_raster_cpu dem = periodic_dem(13, 17);
        const auto result = compute_hillslope_cellular_transport(
            dem, {}, 20.0, cellular_options("pppp"));
        check(result.completed, "closed periodic transport completed");
        check(std::fabs(result.net_dirichlet_export_volume) < 1e-9,
              "closed periodic transport conserves volume in double");
        check(result.n_face_evaluations > 0,
              "face-work counter records the global oracle work");
    }

    // C2: rolling a periodic field only rolls the answer.
    {
        constexpr int nr = 13, nc = 17, shift = 5;
        const lsdtt_raster_cpu dem = periodic_dem(nr, nc);
        std::vector<float> zr(dem.get_data().size(), 0.0f);
        for (int r = 0; r < nr; ++r)
            for (int c = 0; c < nc; ++c)
                zr[static_cast<std::size_t>(r) * nc + ((c + shift) % nc)] =
                    dem.get_data()[static_cast<std::size_t>(r) * nc + c];
        lsdtt_raster_cpu rolled(
            std::move(zr), std::vector<std::uint8_t>(dem.get_mask()), meta(nr, nc));

        const auto a = compute_hillslope_cellular_transport(
            dem, {}, 20.0, cellular_options("pppp"));
        const auto b = compute_hillslope_cellular_transport(
            rolled, {}, 20.0, cellular_options("pppp"));
        auto parallel_options = cellular_options("pppp");
        parallel_options.n_threads = 4;
        const auto parallel = compute_hillslope_cellular_transport(
            dem, {}, 20.0, parallel_options);
        auto multirate_one_options = cellular_options("pppp");
        multirate_one_options.n_threads = 1;
        auto multirate_four_options = multirate_one_options;
        multirate_four_options.n_threads = 4;
        const auto multirate_one = compute_hillslope_cellular_transport_multirate(
            dem, {}, 20.0, multirate_one_options);
        const auto multirate_four = compute_hillslope_cellular_transport_multirate(
            dem, {}, 20.0, multirate_four_options);
        double worst = 0.0;
        for (int r = 0; r < nr; ++r)
            for (int c = 0; c < nc; ++c)
                worst = std::max(worst, std::fabs(
                    static_cast<double>(a.elevation.get_data()[static_cast<std::size_t>(r) * nc + c])
                    - b.elevation.get_data()[static_cast<std::size_t>(r) * nc + ((c + shift) % nc)]));
        check(worst < 1e-7, "periodic transport is translation invariant");
        check(linf(a.elevation, parallel.elevation) == 0.0,
              "topography is deterministic across one and four threads");
        check(multirate_one.completed && multirate_four.completed,
              "periodic multirate transport completed");
        check(linf(multirate_one.elevation, multirate_four.elevation) == 0.0,
              "multirate topography is deterministic across one and four threads");
        check(std::fabs(multirate_one.net_dirichlet_export_volume) < 1e-9,
              "closed periodic multirate transport conserves volume");
    }

    // C3: at small dt, forward cellular and backward-Euler Newton converge.
    {
        constexpr int nr = 9, nc = 17;
        const lsdtt_raster_cpu dem = periodic_dem(nr, nc);
        const auto channels = edge_channels(nr, nc);
        constexpr double dt = 0.01;

        const auto cellular = compute_hillslope_cellular_transport(
            dem, channels, dt, cellular_options());
        const auto multirate = compute_hillslope_cellular_transport_multirate(
            dem, channels, dt, cellular_options());
        auto no = newton_options();
        no.inv_dt = 1.0 / dt;
        const auto newton = compute_hillslope_newton(
            dem, channels, no, dem.get_data());
        check(newton.converged, "small-dt Newton reference converged");
        check(cellular.initial_max_slope_ratio < 0.85,
              "small-dt comparison stays on the exact Roering branch");
        check(linf(cellular.elevation, newton.elevation) < 1e-6,
              "cellular transient converges to Newton as dt approaches zero");
        check(multirate.completed
              && linf(multirate.elevation, cellular.elevation) < 1e-6,
              "one-level multirate update agrees with the global oracle");
    }

    // C4: a narrow pseudo-1D strip relaxes to Newton steady state and U/D
    // curvature without needing a square fixture.
    {
        constexpr int nr = 5, nc = 17;
        lsdtt_hillslope_strip_surface strip = generate_hillslope_strip_surface(
            nr, nc, DX, lsdtt_strip_channel_config::both_ends,
            lsdtt_strip_init_profile::flat, true, 0.0, 0.0,
            U, D, SC, true, 0.0, 1, 1);
        const auto channels = strip_channels(strip.channel_mask);

        auto co = cellular_options("popo");
        co.U = U;
        const double half_width = 0.5 * (nc - 1) * DX;
        const double tau = half_width * half_width / D;
        const auto cellular = compute_hillslope_cellular_transport(
            strip.dem, channels, 8.0 * tau, co);

        auto no = newton_options("popo");
        no.U = static_cast<float>(U);
        const auto newton = compute_hillslope_newton(strip.dem, channels, no);
        check(cellular.completed && newton.converged,
              "pseudo-1D cellular and Newton steady solves completed");
        const double difference = linf(cellular.elevation, newton.elevation);
        check(difference < 2e-3,
              "pseudo-1D cellular field reaches Newton steady state");

        const std::vector<float>& z = cellular.elevation.get_data();
        const std::size_t centre = static_cast<std::size_t>(nr / 2) * nc + nc / 2;
        const double curvature =
            (static_cast<double>(z[centre - 1]) + z[centre + 1] - 2.0 * z[centre])
            / (DX * DX);
        check(std::fabs(curvature + U / D) / (U / D) < 0.02,
              "pseudo-1D ridgetop curvature recovers -U/D within 2 percent");
        check(cellular.net_dirichlet_export_volume > 0.0,
              "steady strip records positive export into fixed channels");
    }

    // C5: supercritical faces are finite and can persist for a short interval.
    {
        constexpr int nr = 17, nc = 33;
        std::vector<float> z(static_cast<std::size_t>(nr) * nc, 0.0f);
        std::vector<std::uint8_t> mask(z.size(), 1u);
        z[static_cast<std::size_t>(nr / 2) * nc + nc / 2] = 6.0f;
        const lsdtt_raster_cpu dem(std::move(z), std::move(mask), meta(nr, nc));
        const auto channels = edge_channels(nr, nc);

        const auto short_step = compute_hillslope_cellular_transport(
            dem, channels, 0.01, cellular_options());
        check(short_step.initial_faces_r100 > 0,
              "supercritical fixture records faces above Sc");
        check(short_step.final_max_slope_ratio > 1.0,
              "a supercritical face can persist for a finite short interval");
        check(short_step.final_max_slope_ratio < short_step.initial_max_slope_ratio,
              "unforced supercritical patch relaxes monotonically at short time");

        auto coarse_options = cellular_options();
        coarse_options.cfl_safety = 0.8;
        auto medium_options = coarse_options;
        medium_options.cfl_safety = 0.4;
        auto fine_options = coarse_options;
        fine_options.cfl_safety = 0.2;
        const auto coarse = compute_hillslope_cellular_transport(
            dem, channels, 100.0, coarse_options);
        const auto medium = compute_hillslope_cellular_transport(
            dem, channels, 100.0, medium_options);
        const auto fine = compute_hillslope_cellular_transport(
            dem, channels, 100.0, fine_options);
        const auto multirate = compute_hillslope_cellular_transport_multirate(
            dem, channels, 100.0, coarse_options);
        check(coarse.completed && medium.completed && fine.completed,
              "supercritical timestep-refinement runs completed");
        check(linf(medium.elevation, fine.elevation)
              < linf(coarse.elevation, medium.elevation),
              "supercritical relaxation converges under timestep refinement");
        check(multirate.completed,
              "localized supercritical multirate update completes without bin violation");
        check(multirate.n_face_events < coarse.n_face_evaluations,
              "multirate scheduler evaluates fewer faces than global subcycling");
    }

    // C6: the small deterministic OCN supports imposed and empty channel sets.
    {
        lsdtt_ocn_options oo;
        oo.n_rows = 32;
        oo.n_cols = 32;
        oo.cell_size = 100.0;
        oo.outlet_mode = ocn_outlet_mode::single_south;
        oo.gamma = 0.5;
        oo.n_iter = 20000;
        oo.temperature_scale = 0.0;
        oo.seed = 7;
        oo.energy_length_weighted = false;
        oo.relief_mode = true;
        oo.target_relief = 20.0f;
        lsdtt_synthetic_network_cpu engine(oo);
        engine.run();
        lsdtt_raster_cpu dem(engine.get_elevation_raster());

        std::vector<std::uint8_t> imposed(engine.areas().size(), 0u);
        for (std::size_t i = 0; i < imposed.size(); ++i)
            imposed[i] = engine.areas()[i] >= 16.0 ? 1u : 0u;

        auto co = cellular_options("nnnn");
        co.D = 0.1;
        constexpr double dt = 250.0;
        const auto cellular_imposed = compute_hillslope_cellular_transport(
            dem, imposed, dt, co);
        const auto cellular_coupled = compute_hillslope_cellular_transport(
            dem, {}, dt, co);

        auto no = newton_options("nnnn");
        no.D = 0.1f;
        no.inv_dt = 1.0 / dt;
        const auto newton_imposed = compute_hillslope_newton(
            dem, imposed, no, dem.get_data());
        const auto newton_coupled = compute_hillslope_newton(
            dem, {}, no, dem.get_data());

        check(cellular_imposed.completed && cellular_coupled.completed
              && newton_imposed.converged && newton_coupled.converged,
              "both small-OCN coupling regimes completed");
        check(cellular_imposed.initial_max_slope_ratio < 0.85
              && cellular_coupled.initial_max_slope_ratio < 0.85,
              "small OCN comparisons stay on the exact Roering branch");
        const double imposed_difference =
            linf(cellular_imposed.elevation, newton_imposed.elevation);
        const double coupled_difference =
            linf(cellular_coupled.elevation, newton_coupled.elevation);
        std::cout << "     OCN imposed Linf(cellular-Newton) = "
                  << imposed_difference << " m, substeps="
                  << cellular_imposed.n_substeps << "\n";
        std::cout << "     OCN coupled Linf(cellular-Newton) = "
                  << coupled_difference << " m, substeps="
                  << cellular_coupled.n_substeps << ", closed export="
                  << cellular_coupled.net_dirichlet_export_volume << " m3\n";
        check(imposed_difference < 2e-3,
              "imposed-channel OCN agrees with Newton for a resolved step");
        check(coupled_difference < 2e-3,
              "empty-channel OCN agrees with Newton for a resolved step");
        check(cellular_imposed.net_dirichlet_export_volume > 0.0,
              "imposed OCN channels receive cellular hillslope export");
        check(std::fabs(cellular_coupled.net_dirichlet_export_volume) < 1e-8,
              "empty-channel closed OCN conserves cellular hillslope volume");
    }

    // C7: LEM orchestration preserves both coupling regimes and process order.
    {
        constexpr int nr = 13, nc = 17;
        const lsdtt_raster_cpu dem = periodic_dem(nr, nc);
        const auto channels = edge_channels(nr, nc);

        lsdtt_lem_options imposed;
        imposed.timestep = 20.0f;
        imposed.enable_d8_fluvial_erosion = false;
        imposed.enable_dinf_fluvial_erosion = false;
        imposed.enable_hillslope_transport = true;
        imposed.enable_uplift = false;
        imposed.hillslope_transport_model = "cellular";
        imposed.hillslope_diffusivity = static_cast<float>(D);
        imposed.critical_hillslope_slope = static_cast<float>(SC);
        imposed.cellular_r_match = 0.85;
        imposed.cellular_cfl_safety = 0.8;
        imposed.coupling_regime = hillslope_coupling_regime::channels_as_boundaries;
        imposed.explicit_channel_mask = channels;
        imposed.boundary_conditions = "nnnn";
        imposed.n_threads = 2;

        lsdtt_hillslope_cellular_workspace lem_ws;
        const auto lem_result = run_lem_step(dem, imposed, {}, nullptr, &lem_ws);
        auto direct_options = cellular_options("nnnn");
        direct_options.D = imposed.hillslope_diffusivity;
        direct_options.Sc = imposed.critical_hillslope_slope;
        direct_options.r_match = imposed.cellular_r_match;
        direct_options.cfl_safety = imposed.cellular_cfl_safety;
        direct_options.n_threads = 2;
        const auto direct = compute_hillslope_cellular_transport(
            dem, channels, imposed.timestep, direct_options);
        check(linf(lem_result.first, direct.elevation) == 0.0,
              "imposed-channel LEM path equals the direct cellular step");
        bool held = true;
        for (std::size_t i = 0; i < channels.size(); ++i)
            if (channels[i] && lem_result.first.get_data()[i] != dem.get_data()[i]) held = false;
        check(held, "imposed-channel LEM path holds channel elevations fixed");

        std::vector<float> slope_z(static_cast<std::size_t>(nr) * nc, 0.0f);
        std::vector<std::uint8_t> slope_mask(slope_z.size(), 1u);
        for (int r = 0; r < nr; ++r)
            for (int c = 0; c < nc; ++c)
                slope_z[static_cast<std::size_t>(r) * nc + c] =
                    static_cast<float>(0.15 * (nr - 1 - r)
                        + 0.05 * std::sin(2.0 * PI * c / nc));
        const lsdtt_raster_cpu slope_dem(
            std::move(slope_z), std::move(slope_mask), meta(nr, nc));

        lsdtt_lem_options coupled = imposed;
        coupled.timestep = 10.0f;
        coupled.enable_d8_fluvial_erosion = true;
        coupled.enable_uplift = true;
        coupled.uplift_rate = static_cast<float>(U);
        coupled.erodibility = 1e-5f;
        coupled.incision_threshold_area = 0.0f;
        coupled.divide_area_by_dx = false;
        coupled.coupling_regime = hillslope_coupling_regime::coupled;
        coupled.explicit_channel_mask.clear();

        std::vector<std::int64_t> fixed_bottom;
        for (int c = 0; c < nc; ++c)
            fixed_bottom.push_back(static_cast<std::int64_t>(nr - 1) * nc + c);
        const auto empty_channels = compute_lem_channel_mask(slope_dem, coupled);
        check(std::count(empty_channels.begin(), empty_channels.end(), 1u) == 0,
              "everywhere-acting cellular LEM has an empty channel set");

        const auto coupled_result = run_lem_step(
            slope_dem, coupled, fixed_bottom, nullptr, &lem_ws);

        // Rebuild the orchestrated result process by process. This both pins
        // ordering and exposes a separate volume budget for each process.
        auto fluvial_options = coupled;
        fluvial_options.enable_hillslope_transport = false;
        fluvial_options.enable_uplift = false;
        const auto fluvial = run_lem_step(
            slope_dem, fluvial_options, {});

        auto coupled_cellular_options = cellular_options("nnnn");
        coupled_cellular_options.D = coupled.hillslope_diffusivity;
        coupled_cellular_options.Sc = coupled.critical_hillslope_slope;
        coupled_cellular_options.r_match = coupled.cellular_r_match;
        coupled_cellular_options.cfl_safety = coupled.cellular_cfl_safety;
        coupled_cellular_options.n_threads = coupled.n_threads;
        const auto hillslope = compute_hillslope_cellular_transport(
            fluvial.first, {}, coupled.timestep, coupled_cellular_options);

        lsdtt_raster_cpu manual = hillslope.elevation;
        for (std::size_t i = 0; i < manual.get_data().size(); ++i)
            if (manual.get_mask()[i])
                manual.get_data()[i] += coupled.uplift_rate * coupled.timestep;
        for (std::int64_t i : fixed_bottom)
            manual.get_data()[static_cast<std::size_t>(i)] =
                slope_dem.get_data()[static_cast<std::size_t>(i)];

        check(linf(coupled_result.first, fluvial.first) > 0.0,
              "cellular transport acts after fluvial incision in coupled mode");
        check(linf(coupled_result.first, manual) == 0.0,
              "manual fluvial-cellular-uplift-restore sequence equals run_lem_step");

        const double cell_area = DX * DX;
        const double fluvial_removed =
            (elevation_sum(slope_dem) - elevation_sum(fluvial.first)) * cell_area;
        const double hillslope_export =
            (elevation_sum(fluvial.first) - elevation_sum(hillslope.elevation)) * cell_area;
        check(fluvial_removed > 0.0,
              "coupled process budget records positive fluvial removal");
        check(std::fabs(hillslope_export) < 1e-3,
              "coupled closed hillslope redistribution conserves volume");
        check(std::fabs(hillslope.net_dirichlet_export_volume) < 1e-9,
              "cellular internal export diagnostic closes the hillslope budget");
        bool bottom_restored = true;
        for (std::int64_t i : fixed_bottom)
            if (coupled_result.first.get_data()[static_cast<std::size_t>(i)]
                != slope_dem.get_data()[static_cast<std::size_t>(i)])
                bottom_restored = false;
        check(bottom_restored,
              "fixed base-level row is restored after fluvial, cellular and uplift processes");
        check(std::isfinite(coupled_result.second.mean_erosion_rate),
              "coupled cellular LEM reports finite erosion statistics");
    }

    // C8: deterministic bounded nonlocal transport preserves the local law
    // below activation, conserves mass, and terminates at fixed channels.
    {
        const lsdtt_raster_cpu low = periodic_dem(13, 17);
        auto local_options = cellular_options("pppp");
        auto nonlocal_options = local_options;
        nonlocal_options.nonlocal_fraction = 0.75;
        nonlocal_options.nonlocal_activation_ratio = 0.8;
        nonlocal_options.nonlocal_deposition_probability = 0.5;
        nonlocal_options.nonlocal_max_hops = 4;
        const auto local_low = compute_hillslope_cellular_transport(
            low, {}, 20.0, local_options);
        const auto nonlocal_low = compute_hillslope_cellular_transport(
            low, {}, 20.0, nonlocal_options);
        check(nonlocal_low.nonlocal_face_events == 0,
              "nonlocal transport is inactive below its slope threshold");
        check(linf(local_low.elevation, nonlocal_low.elevation) == 0.0,
              "nonlocal mode leaves low-slope topography bit-identical");

        constexpr int nr = 17, nc = 33;
        std::vector<float> peak_z(static_cast<std::size_t>(nr) * nc, 0.0f);
        std::vector<std::uint8_t> peak_mask(peak_z.size(), 1u);
        peak_z[static_cast<std::size_t>(nr / 2) * nc + nc / 2] = 6.0f;
        const lsdtt_raster_cpu peak(
            std::move(peak_z), std::move(peak_mask), meta(nr, nc));
        auto nl_one_options = cellular_options("nnnn");
        nl_one_options.nonlocal_fraction = 0.75;
        nl_one_options.nonlocal_activation_ratio = 0.8;
        nl_one_options.nonlocal_deposition_probability = 0.5;
        nl_one_options.nonlocal_max_hops = 4;
        nl_one_options.n_threads = 1;
        auto nl_four_options = nl_one_options;
        nl_four_options.n_threads = 4;
        const auto peak_local = compute_hillslope_cellular_transport(
            peak, {}, 0.1, cellular_options("nnnn"));
        const auto peak_one = compute_hillslope_cellular_transport(
            peak, {}, 0.1, nl_one_options);
        const auto peak_four = compute_hillslope_cellular_transport(
            peak, {}, 0.1, nl_four_options);
        check(peak_one.completed && peak_one.nonlocal_face_events > 0
              && peak_one.nonlocal_routed_volume > 0.0,
              "supercritical patch activates nonlocal routing");
        check(std::fabs(peak_one.net_dirichlet_export_volume) < 1e-9
              && peak_one.nonlocal_export_volume == 0.0,
              "closed nonlocal routing conserves volume");
        check(peak_one.nonlocal_mean_hops >= 1.0
              && peak_one.nonlocal_mean_hops <= nl_one_options.nonlocal_max_hops,
              "nonlocal mean hop count remains within configured bounds");
        check(linf(peak_one.elevation, peak_four.elevation) == 0.0,
              "nonlocal topography is deterministic across thread counts");

        constexpr int rr = 5, cc = 9;
        std::vector<float> ramp_z(static_cast<std::size_t>(rr) * cc, 0.0f);
        std::vector<std::uint8_t> ramp_mask(ramp_z.size(), 1u);
        std::vector<std::uint8_t> ramp_channel(ramp_z.size(), 0u);
        for (int r = 0; r < rr; ++r)
            for (int c = 0; c < cc; ++c)
            {
                ramp_z[static_cast<std::size_t>(r) * cc + c] =
                    static_cast<float>(4.5 * (cc - 1 - c));
                if (c == cc - 1)
                    ramp_channel[static_cast<std::size_t>(r) * cc + c] = 1u;
            }
        const lsdtt_raster_cpu ramp(
            std::move(ramp_z), std::move(ramp_mask), meta(rr, cc));
        auto ramp_options = cellular_options("nnnn");
        ramp_options.nonlocal_fraction = 1.0;
        ramp_options.nonlocal_activation_ratio = 0.7;
        ramp_options.nonlocal_deposition_probability = 0.1;
        ramp_options.nonlocal_max_hops = 8;
        const auto routed = compute_hillslope_cellular_transport(
            ramp, ramp_channel, 0.01, ramp_options);
        auto hop2_options = ramp_options;
        hop2_options.nonlocal_max_hops = 2;
        auto hop4_options = ramp_options;
        hop4_options.nonlocal_max_hops = 4;
        const auto hop2 = compute_hillslope_cellular_transport(
            ramp, ramp_channel, 0.01, hop2_options);
        const auto hop4 = compute_hillslope_cellular_transport(
            ramp, ramp_channel, 0.01, hop4_options);
        auto ramp_local_options = cellular_options("nnnn");
        const auto ramp_local = compute_hillslope_cellular_transport(
            ramp, ramp_channel, 0.01, ramp_local_options);
        check(linf(routed.elevation, ramp_local.elevation) > 0.0,
              "nonlocal routing changes deposition along a continuing slope");
        check(routed.nonlocal_export_volume > 0.0,
              "nonlocal paths stop and export at the first fixed channel");
        check(hop2.nonlocal_mean_hops <= hop4.nonlocal_mean_hops
              && hop4.nonlocal_mean_hops <= routed.nonlocal_mean_hops,
              "mean travel distance increases with the bounded hop scale");
        bool channel_held = true;
        for (std::size_t i = 0; i < ramp_channel.size(); ++i)
            if (ramp_channel[i] && routed.elevation.get_data()[i] != ramp.get_data()[i])
                channel_held = false;
        check(channel_held, "nonlocal interception does not alter fixed channel cells");

        constexpr int dn = 9;
        std::vector<float> diagonal_z(static_cast<std::size_t>(dn) * dn, 0.0f);
        std::vector<std::uint8_t> diagonal_mask(diagonal_z.size(), 1u);
        std::vector<std::uint8_t> diagonal_channels(diagonal_z.size(), 0u);
        for (int r = 0; r < dn; ++r)
            for (int c = 0; c < dn; ++c)
            {
                diagonal_z[static_cast<std::size_t>(r) * dn + c] =
                    static_cast<float>(3.2 * (2 * dn - 2 - r - c));
                if (r == dn - 1 || c == dn - 1)
                    diagonal_channels[static_cast<std::size_t>(r) * dn + c] = 1u;
            }
        const lsdtt_raster_cpu diagonal(
            std::move(diagonal_z), std::move(diagonal_mask), meta(dn, dn));
        const auto diagonal_result = compute_hillslope_cellular_transport(
            diagonal, diagonal_channels, 0.01, ramp_options);
        double transpose_error = 0.0;
        for (int r = 0; r < dn; ++r)
            for (int c = 0; c < dn; ++c)
                transpose_error = std::max(transpose_error, std::fabs(
                    static_cast<double>(diagonal_result.elevation.get_data()[
                        static_cast<std::size_t>(r) * dn + c])
                    - diagonal_result.elevation.get_data()[
                        static_cast<std::size_t>(c) * dn + r]));
        check(transpose_error < 1e-6,
              "nonlocal routing preserves transpose symmetry on a diagonal slope");
    }

    // C9: sparse Lagrangian particles sample accepted flux without changing
    // topography, provenance, or thread determinism.
    {
        const lsdtt_raster_cpu density_dem = periodic_dem(100, 100);
        const auto d1000 = seed_sparse_hillslope_particles(
            density_dem, 1000, 1.0, 17);
        const auto d100 = seed_sparse_hillslope_particles(
            density_dem, 100, 1.0, 17);
        const auto d10 = seed_sparse_hillslope_particles(
            density_dem, 10, 1.0, 17);
        check(d1000.size() == 10 && d100.size() == 100 && d10.size() == 1000,
              "systematic seeding produces 1:1000, 1:100 and 1:10 densities");
        double seeded_weight = 0.0;
        for (double w : d100.weight) seeded_weight += w;
        check(std::fabs(seeded_weight - 10000.0 * DX * DX) < 1e-8,
              "particle weights exactly represent valid active-layer volume");

        constexpr int nr = 9, nc = 33;
        std::vector<float> z(static_cast<std::size_t>(nr) * nc, 0.0f);
        std::vector<std::uint8_t> mask(z.size(), 1u);
        std::vector<std::uint8_t> channels(z.size(), 0u);
        for (int r = 0; r < nr; ++r)
            for (int c = 0; c < nc; ++c)
            {
                z[static_cast<std::size_t>(r) * nc + c] =
                    static_cast<float>(4.5 * (nc - 1 - c));
                if (c == nc - 1) channels[static_cast<std::size_t>(r) * nc + c] = 1u;
            }
        const lsdtt_raster_cpu ramp(
            std::move(z), std::move(mask), meta(nr, nc));
        auto particles = seed_sparse_hillslope_particles(
            ramp, 1, 0.001, 29, channels);
        for (std::size_t p = 0; p < particles.size(); ++p)
        {
            const int c = static_cast<int>(particles.cell[p] % nc);
            particles.source_id[p] = c < nc / 2 ? 1u : 2u;
        }
        const auto initial_particles = particles;
        auto particles_four = particles;

        auto options = cellular_options("nnnn");
        options.nonlocal_fraction = 1.0;
        options.nonlocal_activation_ratio = 0.7;
        options.nonlocal_deposition_probability = 0.5;
        options.nonlocal_max_hops = 8;
        options.n_threads = 1;
        const auto untraced = compute_hillslope_cellular_transport(
            ramp, channels, 0.01, options);
        const auto traced = compute_hillslope_cellular_transport(
            ramp, channels, 0.01, options, nullptr, &particles);
        auto options_four = options;
        options_four.n_threads = 4;
        const auto traced_four = compute_hillslope_cellular_transport(
            ramp, channels, 0.01, options_four, nullptr, &particles_four);

        check(linf(untraced.elevation, traced.elevation) == 0.0,
              "enabling sparse particles does not change topography");
        check(linf(traced.elevation, traced_four.elevation) == 0.0,
              "traced topography is thread deterministic");
        check(particles.cell == particles_four.cell
              && particles.status == particles_four.status
              && particles.distance_cells == particles_four.distance_cells,
              "particle trajectories are thread deterministic");
        check(traced.tracker_moves > 0 && traced.tracker_exports > 0,
              "particles move and export through the fixed channel");
        check(*std::max_element(particles.distance_cells.begin(),
                                particles.distance_cells.end()) > 1.0,
              "particles sample nonlocal paths longer than one cell");
        check(particles.source_id == initial_particles.source_id,
              "one- and two-source provenance IDs survive transport unchanged");
        double source1_before = 0.0, source1_after = 0.0;
        double source2_before = 0.0, source2_after = 0.0;
        for (std::size_t p = 0; p < particles.size(); ++p)
        {
            if (initial_particles.source_id[p] == 1u)
                source1_before += initial_particles.weight[p];
            else source2_before += initial_particles.weight[p];
            if (particles.source_id[p] == 1u) source1_after += particles.weight[p];
            else source2_after += particles.weight[p];
        }
        check(source1_before == source1_after && source2_before == source2_after,
              "weighted two-source provenance budget closes exactly");

        lsdtt_lem_options lem;
        lem.timestep = 0.01f;
        lem.enable_d8_fluvial_erosion = false;
        lem.enable_dinf_fluvial_erosion = false;
        lem.enable_hillslope_transport = true;
        lem.enable_uplift = false;
        lem.hillslope_transport_model = "cellular";
        lem.hillslope_diffusivity = static_cast<float>(D);
        lem.critical_hillslope_slope = static_cast<float>(SC);
        lem.cellular_nonlocal_fraction = 1.0;
        lem.cellular_nonlocal_activation_ratio = 0.7;
        lem.cellular_nonlocal_deposition_probability = 0.5;
        lem.cellular_nonlocal_max_hops = 8;
        lem.coupling_regime = hillslope_coupling_regime::channels_as_boundaries;
        lem.explicit_channel_mask = channels;
        lem.boundary_conditions = "nnnn";
        lem.n_threads = 1;
        auto lem_particles = initial_particles;
        lsdtt_hillslope_cellular_workspace lem_ws;
        const auto lem_result = run_lem_step(
            ramp, lem, {}, nullptr, &lem_ws, &lem_particles);
        check(linf(lem_result.first, traced.elevation) == 0.0
              && lem_particles.cell == particles.cell
              && lem_particles.status == particles.status,
              "LEM forwards the identical sparse-particle transport step");

        // Exact Eulerian/Lagrangian export check: one active cell surrounded
        // by fixed sinks. The topographic export thickness determines the
        // well-mixed active-layer survival probability analytically.
        constexpr int en = 3;
        std::vector<float> export_z(en * en, 0.0f);
        std::vector<std::uint8_t> export_mask(en * en, 1u);
        std::vector<std::uint8_t> export_channels(en * en, 1u);
        const std::size_t centre = 4;
        export_z[centre] = 1.0f;
        export_channels[centre] = 0u;
        const lsdtt_raster_cpu export_dem(
            std::move(export_z), std::move(export_mask), meta(en, en));

        constexpr std::size_t np = 20000;
        lsdtt_hillslope_particle_tracker export_particles;
        export_particles.mixing_depth = 0.001;
        export_particles.seed = 73;
        export_particles.id.resize(np);
        export_particles.cell.assign(np, static_cast<std::int64_t>(centre));
        export_particles.weight.assign(np, DX * DX * export_particles.mixing_depth / np);
        export_particles.source_id.resize(np);
        export_particles.age.assign(np, 0.0);
        export_particles.distance_cells.assign(np, 0.0);
        export_particles.status.assign(np, static_cast<std::uint8_t>(
            lsdtt_hillslope_particle_status::active));
        for (std::size_t p = 0; p < np; ++p)
        {
            export_particles.id[p] = p;
            export_particles.source_id[p] = p < np / 2 ? 1u : 2u;
        }
        const auto export_result = compute_hillslope_cellular_transport(
            export_dem, export_channels, 0.05, cellular_options("nnnn"),
            nullptr, &export_particles);
        const double exported_thickness =
            export_result.net_dirichlet_export_volume / (DX * DX);
        const double expected_fraction = 1.0 - std::exp(
            -exported_thickness / export_particles.mixing_depth);
        const double observed_fraction =
            static_cast<double>(export_result.tracker_exports) / np;
        const double sigma = std::sqrt(expected_fraction
            * (1.0 - expected_fraction) / np);
        check(std::fabs(observed_fraction - expected_fraction) < 6.0 * sigma,
              "particle export agrees with the Eulerian active-layer budget");
        double exported_weight = 0.0;
        for (std::size_t p = 0; p < np; ++p)
            if (export_particles.status[p] == static_cast<std::uint8_t>(
                    lsdtt_hillslope_particle_status::exported))
                exported_weight += export_particles.weight[p];
        const double expected_weight = DX * DX * export_particles.mixing_depth
                                     * expected_fraction;
        check(std::fabs(exported_weight - expected_weight)
              < 6.0 * sigma * DX * DX * export_particles.mixing_depth,
              "weighted particle export closes against expected exported inventory");
    }

    std::cout << "\n" << (g_checks - g_failures) << " / " << g_checks
              << " checks passed\n";
    return g_failures == 0 ? 0 : 1;
}
