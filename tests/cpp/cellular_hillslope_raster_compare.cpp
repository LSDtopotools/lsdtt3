// cellular_hillslope_raster_compare.cpp
// Native field comparison for experimental E2 outputs.

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iostream>

#include "../../src/core/lsdtt_data_loader.hpp"
#include "../../src/cpu/lsdtt_raster_cpu.hpp"

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
              << "b_relief_m," << bmax - bmin << '\n';
    return 0;
}
