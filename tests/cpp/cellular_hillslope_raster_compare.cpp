// cellular_hillslope_raster_compare.cpp
// Native field comparison for experimental E2 outputs.

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iostream>

#include "../../src/core/lsdtt_data_loader.hpp"
#include "../../src/cpu/lsdtt_raster_cpu.hpp"

namespace
{
double curvature_flip_fraction(const lsdtt::lsdtt_raster_cpu& dem)
{
    const int nr = dem.get_meta().n_rows;
    const int nc = dem.get_meta().n_cols;
    const auto& z = dem.get_data();
    const auto& mask = dem.get_mask();
    auto at = [&](int r, int c) -> double
    { return z[static_cast<std::size_t>(r) * nc + c]; };
    auto valid = [&](int r, int c) -> bool
    { return mask[static_cast<std::size_t>(r) * nc + c] != 0u; };

    std::size_t flips = 0, pairs = 0;
    for (int r = 1; r < nr - 1; ++r)
        for (int c = 1; c < nc - 2; ++c)
            if (valid(r, c - 1) && valid(r, c) && valid(r, c + 1)
                && valid(r, c + 2))
            {
                const double a = at(r, c - 1) - 2.0 * at(r, c) + at(r, c + 1);
                const double b = at(r, c) - 2.0 * at(r, c + 1) + at(r, c + 2);
                if (a != 0.0 && b != 0.0)
                {
                    ++pairs;
                    if ((a < 0.0) != (b < 0.0)) ++flips;
                }
            }
    for (int r = 1; r < nr - 2; ++r)
        for (int c = 1; c < nc - 1; ++c)
            if (valid(r - 1, c) && valid(r, c) && valid(r + 1, c)
                && valid(r + 2, c))
            {
                const double a = at(r - 1, c) - 2.0 * at(r, c) + at(r + 1, c);
                const double b = at(r, c) - 2.0 * at(r + 1, c) + at(r + 2, c);
                if (a != 0.0 && b != 0.0)
                {
                    ++pairs;
                    if ((a < 0.0) != (b < 0.0)) ++flips;
                }
            }
    return pairs ? static_cast<double>(flips) / pairs : 0.0;
}
}

int main(int argc, char** argv)
{
    if (argc != 3)
    {
        std::cerr << "usage: cellular-hillslope-raster-compare A B\n";
        return 2;
    }
    lsdtt::lsdtt_data_loader la, lb;
    la.open_raster(argv[1]);
    lb.open_raster(argv[2]);
    lsdtt::lsdtt_raster_cpu a(la), b(lb);
    if (a.get_data().size() != b.get_data().size())
    {
        std::cerr << "raster size mismatch\n";
        return 1;
    }

    double linf = 0.0, sum = 0.0, sum2 = 0.0;
    double amin = 0.0, amax = 0.0, bmin = 0.0, bmax = 0.0;
    std::size_t n = 0;
    bool first = true;
    for (std::size_t i = 0; i < a.get_data().size(); ++i)
    {
        if (!a.get_mask()[i] || !b.get_mask()[i]) continue;
        const double av = a.get_data()[i];
        const double bv = b.get_data()[i];
        const double d = av - bv;
        linf = std::max(linf, std::fabs(d));
        sum += d;
        sum2 += d * d;
        if (first)
        {
            amin = amax = av;
            bmin = bmax = bv;
            first = false;
        }
        else
        {
            amin = std::min(amin, av); amax = std::max(amax, av);
            bmin = std::min(bmin, bv); bmax = std::max(bmax, bv);
        }
        ++n;
    }
    std::cout << "valid_cells," << n << '\n'
              << "linf_m," << linf << '\n'
              << "rms_m," << (n ? std::sqrt(sum2 / n) : 0.0) << '\n'
              << "mean_a_minus_b_m," << (n ? sum / n : 0.0) << '\n'
              << "a_relief_m," << amax - amin << '\n'
              << "b_relief_m," << bmax - bmin << '\n'
              << "a_curvature_flip_fraction," << curvature_flip_fraction(a) << '\n'
              << "b_curvature_flip_fraction," << curvature_flip_fraction(b) << '\n';
    return 0;
}
