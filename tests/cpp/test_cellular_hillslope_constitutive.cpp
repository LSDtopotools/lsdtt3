// test_cellular_hillslope_constitutive.cpp
//
// Phase P1 of plan_cellular_hillslope_transport.md: test the proposed
// soft-critical mobility without implementing a raster transport solver, and
// turn its tangent modulus into an explicit whole-grid pass-cost envelope.

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

namespace
{

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

double roering_mobility(double r)
{
    return 1.0 / (1.0 - r * r);
}

double roering_mobility_derivative(double r)
{
    const double d = 1.0 - r * r;
    return 2.0 * r / (d * d);
}

double soft_mobility(double r, double r_match)
{
    if (r <= r_match) return roering_mobility(r);
    return roering_mobility(r_match)
         + roering_mobility_derivative(r_match) * (r - r_match);
}

double soft_flux(double r, double r_match)
{
    return r * soft_mobility(r, r_match);
}

// Derivative of the dimensionless flux r*M(r), which controls the explicit
// diffusion stability limit. This, not M alone, is the tangent modulus.
double soft_flux_tangent(double r, double r_match)
{
    if (r <= r_match)
    {
        const double d = 1.0 - r * r;
        return (1.0 + r * r) / (d * d);
    }
    const double mp = roering_mobility_derivative(r_match);
    return soft_mobility(r, r_match) + r * mp;
}

double roering_flux_tangent(double r)
{
    const double d = 1.0 - r * r;
    return (1.0 + r * r) / (d * d);
}

struct cost_case
{
    std::string name;
    long cells;
    double dx;
    double D;
    double dt;
    double max_r;
    double r_match;
    double newton_ms; // < 0 where the E2 Newton timing is not yet available
};

long predicted_substeps(const cost_case& c)
{
    const double tangent = soft_flux_tangent(c.max_r, c.r_match);
    return std::max(1L, static_cast<long>(std::ceil(
        4.0 * c.D * c.dt * tangent / (c.dx * c.dx))));
}

void write_law_table(const std::filesystem::path& path)
{
    std::ofstream out(path);
    if (!out) throw std::runtime_error("cannot write constitutive table: " + path.string());
    out << "r_match,slope_ratio,mobility,dimensionless_flux,flux_tangent\n";
    out << std::setprecision(17);
    for (double rm : {0.80, 0.85, 0.90})
        for (int i = 0; i <= 150; ++i)
        {
            const double r = static_cast<double>(i) / 100.0;
            out << rm << ',' << r << ',' << soft_mobility(r, rm) << ','
                << soft_flux(r, rm) << ',' << soft_flux_tangent(r, rm) << '\n';
        }
}

void write_cost_table(const std::filesystem::path& path,
                      const std::vector<cost_case>& cases)
{
    std::ofstream out(path);
    if (!out) throw std::runtime_error("cannot write cost table: " + path.string());
    out << "case,cells,dx_m,D_m2_per_yr,dt_yr,max_slope_ratio,r_match,"
           "max_flux_tangent,predicted_substeps,predicted_streaming_passes,"
           "newton_ms,allowed_ns_per_cell_pass_for_2x\n";
    out << std::setprecision(17);
    for (const cost_case& c : cases)
    {
        const long substeps = predicted_substeps(c);
        const long passes = 2L * substeps; // face compute plus cell gather
        double allowance = -1.0;
        if (c.newton_ms > 0.0)
            allowance = (0.5 * c.newton_ms * 1.0e6)
                      / (static_cast<double>(passes) * static_cast<double>(c.cells));
        out << c.name << ',' << c.cells << ',' << c.dx << ',' << c.D << ','
            << c.dt << ',' << c.max_r << ',' << c.r_match << ','
            << soft_flux_tangent(c.max_r, c.r_match) << ',' << substeps << ','
            << passes << ',' << c.newton_ms << ',' << allowance << '\n';
    }
}

} // namespace

int main(int argc, char** argv)
{
    std::cout << "\n=== Soft-critical cellular constitutive law ===\n";

    for (double rm : {0.80, 0.85, 0.90})
    {
        check(soft_mobility(0.0, rm) == 1.0,
              "r_match=" + std::to_string(rm) + " has linear mobility at zero slope");
        check(soft_flux(0.0, rm) == 0.0,
              "r_match=" + std::to_string(rm) + " has zero flux at zero slope");
        check(soft_flux_tangent(0.0, rm) == 1.0,
              "r_match=" + std::to_string(rm) + " has unit low-slope tangent");

        const double m_at = soft_mobility(rm, rm);
        const double f_at = soft_flux(rm, rm);
        const double t_at = soft_flux_tangent(rm, rm);
        const double eps = 1e-9;
        check(std::fabs(soft_mobility(rm + eps, rm) - m_at) < 1e-7,
              "mobility is continuous at r_match=" + std::to_string(rm));
        check(std::fabs(soft_flux(rm + eps, rm) - f_at) < 1e-7,
              "flux is continuous at r_match=" + std::to_string(rm));
        check(std::fabs(soft_flux_tangent(rm + eps, rm) - t_at) < 1e-6,
              "flux tangent is continuous at r_match=" + std::to_string(rm));

        bool finite = true;
        bool monotonic = true;
        double previous_flux = -1.0;
        for (int i = 0; i <= 1500; ++i)
        {
            const double r = static_cast<double>(i) / 1000.0;
            const double m = soft_mobility(r, rm);
            const double f = soft_flux(r, rm);
            const double t = soft_flux_tangent(r, rm);
            finite = finite && std::isfinite(m) && std::isfinite(f)
                    && std::isfinite(t) && m > 0.0 && t > 0.0;
            monotonic = monotonic && f >= previous_flux;
            previous_flux = f;
            if (r < rm)
                finite = finite && std::fabs(m - roering_mobility(r)) < 1e-12;
        }
        check(finite, "law is finite and positive over r=0..1.5 for r_match="
                      + std::to_string(rm));
        check(monotonic, "flux is monotonic over r=0..1.5 for r_match="
                         + std::to_string(rm));
        check(std::isfinite(soft_flux(1.0, rm)),
              "flux remains finite at Sc for r_match=" + std::to_string(rm));
    }

    const double soft_tangent_at_sc = soft_flux_tangent(1.0, 0.85);
    const double singular_tangent_near_sc = roering_flux_tangent(0.99);
    std::cout << "\nr_match=0.85: tangent at Sc = " << soft_tangent_at_sc
              << ", Roering tangent at 0.99 Sc = " << singular_tangent_near_sc
              << "\n";
    check(soft_tangent_at_sc < 30.0,
          "r_match=0.85 keeps the tangent at Sc below 30");
    check(singular_tangent_near_sc / soft_tangent_at_sc > 100.0,
          "continuation removes more than 100x of the near-singular tangent");

    // P0 cases use measured Newton medians. The E2 rows are a documented cost
    // envelope (600x600, dx=5 m, dt=5000 yr); their Newton timings remain -1
    // until the staged E2 runs produce landscapes.
    const std::vector<cost_case> cases = {
        {"p0_ocn_coupled",       1024, 100.0, 0.1000, 1000.0, 0.0393, 0.85, 26.16},
        {"p0_incision_pulse",     903,   5.0, 0.1000,  500.0, 0.7342, 0.85, 15.81},
        {"p0_supercritical_patch",561,   5.0, 0.1000,  100.0, 1.5000, 0.85, 27.78},
        {"e2_D00077_r085",     360000,   5.0, 0.0077, 5000.0, 0.8500, 0.85, -1.0},
        {"e2_D0064_r085",      360000,   5.0, 0.0640, 5000.0, 0.8500, 0.85, -1.0},
        {"e2_D0385_r085",      360000,   5.0, 0.3850, 5000.0, 0.8500, 0.85, -1.0},
        {"e2_D00077_r100",     360000,   5.0, 0.0077, 5000.0, 1.0000, 0.85, -1.0},
        {"e2_D0064_r100",      360000,   5.0, 0.0640, 5000.0, 1.0000, 0.85, -1.0},
        {"e2_D0385_r100",      360000,   5.0, 0.3850, 5000.0, 1.0000, 0.85, -1.0}
    };

    std::cout << "\n-- Global-CFL cost envelope --\n";
    for (const cost_case& c : cases)
    {
        const long substeps = predicted_substeps(c);
        std::cout << "  " << std::setw(27) << std::left << c.name
                  << " substeps=" << std::setw(6) << substeps
                  << " passes=" << 2L * substeps;
        if (c.newton_ms > 0.0)
        {
            const double allowance = (0.5 * c.newton_ms * 1.0e6)
                / (static_cast<double>(2L * substeps) * static_cast<double>(c.cells));
            std::cout << "  2x allowance=" << allowance << " ns/cell/pass";
        }
        std::cout << '\n';
    }

    check(predicted_substeps(cases[0]) == 1,
          "gentle 100 m OCN needs one global substep");
    check(predicted_substeps(cases[1]) > 40,
          "incision pulse exposes material global subcycling cost");
    check(predicted_substeps(cases[5]) > 6000,
          "high-D E2 envelope exposes prohibitive global subcycling");

    const std::filesystem::path output_dir = argc > 1
        ? std::filesystem::path(argv[1]) : std::filesystem::path(".");
    std::filesystem::create_directories(output_dir);
    write_law_table(output_dir / "cellular_hillslope_constitutive.csv");
    write_cost_table(output_dir / "cellular_hillslope_cost_envelope.csv", cases);
    std::cout << "\nTables written to " << output_dir.string() << "\n";
    std::cout << (g_checks - g_failures) << " / " << g_checks
              << " checks passed\n";
    return g_failures == 0 ? 0 : 1;
}
