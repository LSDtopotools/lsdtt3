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
        double worst = 0.0;
        for (int r = 0; r < nr; ++r)
            for (int c = 0; c < nc; ++c)
                worst = std::max(worst, std::fabs(
                    static_cast<double>(a.elevation.get_data()[static_cast<std::size_t>(r) * nc + c])
                    - b.elevation.get_data()[static_cast<std::size_t>(r) * nc + ((c + shift) % nc)]));
        check(worst < 1e-7, "periodic transport is translation invariant");
        check(linf(a.elevation, parallel.elevation) == 0.0,
              "topography is deterministic across one and four threads");
    }

    // C3: at small dt, forward cellular and backward-Euler Newton converge.
    {
        constexpr int nr = 9, nc = 17;
        const lsdtt_raster_cpu dem = periodic_dem(nr, nc);
        const auto channels = edge_channels(nr, nc);
        constexpr double dt = 0.01;

        const auto cellular = compute_hillslope_cellular_transport(
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
        check(coarse.completed && medium.completed && fine.completed,
              "supercritical timestep-refinement runs completed");
        check(linf(medium.elevation, fine.elevation)
              < linf(coarse.elevation, medium.elevation),
              "supercritical relaxation converges under timestep refinement");
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

    std::cout << "\n" << (g_checks - g_failures) << " / " << g_checks
              << " checks passed\n";
    return g_failures == 0 ? 0 : 1;
}
