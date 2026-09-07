#include <cmath>
#include <cstdint>
#include <iostream>
#include <utility>
#include <vector>

#include "lsdtt_raster_hierarchy_cpu.hpp"

namespace {

int checks = 0;
int failures = 0;

void check(const char* name, bool ok)
{
    ++checks;
    if (!ok) ++failures;
    std::cout << (ok ? "PASS" : "FAIL") << "  " << name << '\n';
}

lsdtt::lsdtt_raster_cpu source()
{
    lsdtt::lsdtt_raster_meta meta;
    meta.n_rows = 8;
    meta.n_cols = 8;
    meta.cell_size = 1.0;
    meta.no_data_value = -9999.0;
    std::vector<float> values(64);
    for (std::size_t i = 0; i < values.size(); ++i) values[i] = static_cast<float>(i);
    return {std::move(values), std::vector<std::uint8_t>(64, 1), meta};
}

} // namespace

int main()
{
    auto dem = source();
    lsdtt::lsdtt_raster_hierarchy_cpu hierarchy(dem);

    check("dyadic layout", hierarchy.n_levels() == 4 &&
          hierarchy.level_meta(2).n_rows == 2 &&
          std::abs(hierarchy.level_meta(2).cell_size - 4.0) < 1e-12);
    check("finest retained", hierarchy.has_level(0) && !hierarchy.has_level(1));

    const auto& coarse = hierarchy.restrict_to(2);
    check("restriction activates target", hierarchy.active_level() == 2 &&
          coarse.n_rows() == 2 && coarse.n_cols() == 2 &&
          std::abs(coarse.get_data()[0] - 13.5f) < 1e-6f);

    const auto& fine = hierarchy.prolong_to(0);
    check("prolongation activates finest", hierarchy.active_level() == 0 &&
          fine.n_rows() == 8 && fine.n_cols() == 8);

    lsdtt::lsdtt_raster_hierarchy_cpu active_only(
        dem, lsdtt::raster_hierarchy_storage::active_only, false);
    active_only.restrict_to(1);
    check("active-only releases source", active_only.has_level(1) && !active_only.has_level(0));

    lsdtt::lsdtt_raster_hierarchy_cpu preallocated(
        dem, lsdtt::raster_hierarchy_storage::preallocated);
    check("preallocation", preallocated.has_level(0) && preallocated.has_level(3));

    lsdtt::lsdtt_raster_meta five_meta;
    five_meta.n_rows = 5;
    five_meta.n_cols = 5;
    five_meta.cell_size = 1.0;
    five_meta.no_data_value = -9999.0;
    lsdtt::lsdtt_raster_cpu five(
        std::vector<float>(25, 1.0f), std::vector<std::uint8_t>(25, 1), five_meta);
    lsdtt::lsdtt_raster_hierarchy_cpu arbitrary(five);
    const auto factor_five = arbitrary.add_level_for_factor(5);
    arbitrary.restrict_to(factor_five, lsdtt::raster_resampling_method::mean);
    check("dynamic integer factor", arbitrary.level(factor_five).n_rows() == 1 &&
          std::abs(arbitrary.level(factor_five).cell_size() - 5.0) < 1e-12);

    std::cout << checks - failures << '/' << checks << " checks passed\n";
    return failures == 0 ? 0 : 1;
}
