// test_canopy_interception.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
//
// -----------------------------------------------------------------------------
// Unit tests for src/cpu/lsdtt_canopy_interception_cpu -- Phase 4 of the
// distributed-hydrology work.
//
// What is under test
// ------------------
//   the LAI closed forms (capacity and gap fraction) and their limits
//   the Rutter bucket's mass balance
//   THE EVAPORATIVE RE-PARTITION, including total ET legitimately exceeding PET
//   seasonality, at the magnitude that matters for a winter UK storm
//   the land-cover -> LAI ingest, including unmatched classes
//   the composed canopy + column balance across BOTH stores
//
// THE CHECK THAT MATTERS MOST
// ---------------------------
// test_evaporative_repartition. Interception is not a skim off the top: wet
// canopy evaporation suppresses transpiration, so only the residual demand
// reaches the root zone, and the residual is PET - E_canopy/evap_factor rather
// than PET - E_canopy. Dropping the division caps total ET at PET and silently
// understates forest water use. It is the sort of error that leaves every
// balance closing and every store non-negative.
//
// Usage: lsdtt-test-canopy-interception
// -----------------------------------------------------------------------------

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

#include "lsdtt_canopy_interception_cpu.hpp"
#include "lsdtt_topmodel_cpu.hpp"

namespace
{

int g_n_checks = 0;
int g_n_failed = 0;

void check(const std::string& name, bool ok, const std::string& detail = "")
{
    ++g_n_checks;
    if (!ok) ++g_n_failed;
    std::cout << (ok ? "PASS" : "FAIL") << "  " << name;
    if (!detail.empty()) std::cout << "   [" << detail << "]";
    std::cout << std::endl;
}

bool close_rel(double a, double b, double tol)
{
    const double d = std::max(std::fabs(a), std::fabs(b));
    if (d < 1e-300) return std::fabs(a - b) < tol;
    return std::fabs(a - b) / d < tol;
}

std::vector<double> equal_areas(std::size_t n)
{
    return std::vector<double>(n, 1.0 / static_cast<double>(n));
}

// =============================================================================
// 1. The LAI closed forms
// =============================================================================
void test_lai_forms()
{
    // C = 0.2 mm per unit LAI, returned in metres.
    check("lai: capacity is k_S * LAI, converted to metres",
          close_rel(lsdtt::canopy_capacity_from_lai(5.0, 0.2), 5.0 * 0.2e-3, 1e-15),
          "C(LAI=5) = " + std::to_string(lsdtt::canopy_capacity_from_lai(5.0, 0.2)) + " m");
    check("lai: capacity is linear in LAI",
          close_rel(lsdtt::canopy_capacity_from_lai(10.0, 0.2),
                    2.0 * lsdtt::canopy_capacity_from_lai(5.0, 0.2), 1e-15));
    check("lai: bare ground has zero capacity",
          lsdtt::canopy_capacity_from_lai(0.0, 0.2) == 0.0);

    // p = exp(-kappa LAI).
    check("lai: gap fraction is exp(-kappa LAI)",
          close_rel(lsdtt::free_throughfall_from_lai(3.0, 0.6),
                    std::exp(-0.6 * 3.0), 1e-15));
    check("lai: bare ground passes ALL rain",
          lsdtt::free_throughfall_from_lai(0.0, 0.6) == 1.0);
    check("lai: a dense canopy passes almost none",
          lsdtt::free_throughfall_from_lai(10.0, 0.6) < 0.01,
          "p(LAI=10) = " + std::to_string(lsdtt::free_throughfall_from_lai(10.0, 0.6)));
    check("lai: gap fraction falls monotonically with LAI",
          lsdtt::free_throughfall_from_lai(1.0, 0.6) >
          lsdtt::free_throughfall_from_lai(2.0, 0.6));
}

// =============================================================================
// 2. Mass balance of the bucket
// =============================================================================
void test_balance()
{
    lsdtt::canopy_parameters p;
    p.canopy_storage_mm_per_lai = 0.2;
    p.lai_gap_coefficient       = 0.6;

    const std::vector<double> lai = {0.0, 1.0, 3.0, 5.0, 8.0};
    lsdtt::lsdtt_canopy_interception_cpu c(lai, equal_areas(lai.size()), p);

    for (int t = 0; t < 2000; ++t)
    {
        const bool storm = (t % 37) < 8;
        (void)c.step(storm ? 0.003 : 0.0, 0.0002);
    }

    const auto& b = c.water_balance();
    const double scale = std::max(b.total_rainfall, 1e-12);
    check("balance: the canopy store closes over 2000 steps",
          std::fabs(b.closure_error()) / scale < 1e-12,
          "relative error " + std::to_string(std::fabs(b.closure_error()) / scale) +
          ", P=" + std::to_string(b.total_rainfall) +
          " TF=" + std::to_string(b.total_throughfall) +
          " E=" + std::to_string(b.total_canopy_evaporation));

    check("balance: throughfall never exceeds rainfall in total",
          b.total_throughfall <= b.total_rainfall + 1e-12);

    bool store_bounded = true;
    for (std::size_t i = 0; i < c.n_units(); ++i)
        if (c.canopy_storage()[i] < -1e-15 ||
            c.canopy_storage()[i] > c.capacity()[i] + 1e-15) store_bounded = false;
    check("balance: every store stays within [0, capacity]", store_bounded);
}

// =============================================================================
// 3. THE evaporative re-partition
// =============================================================================
void test_evaporative_repartition()
{
    const double pet = 0.001;

    // With evap_factor = 1 the residual is simply PET - E_canopy.
    {
        lsdtt::canopy_parameters p;
        p.interception_evap_factor = 1.0;
        lsdtt::lsdtt_canopy_interception_cpu c({4.0}, {1.0}, p);
        const auto w = c.step(0.002, pet);          // wet the canopy
        check("repartition: at factor 1, residual = PET - E_canopy",
              close_rel(w.residual_pet, pet - w.canopy_evaporation, 1e-12),
              "E=" + std::to_string(w.canopy_evaporation) +
              " residual=" + std::to_string(w.residual_pet));
    }

    // With evap_factor = 2 a rough canopy evaporates at twice PET, but consumes
    // only half as much transpiration demand per unit evaporated -- so TOTAL ET
    // exceeds PET. This is the Rutter/Calder result and the reason the division
    // by evap_factor is in the code.
    {
        lsdtt::canopy_parameters p;
        p.interception_evap_factor  = 2.0;
        p.canopy_storage_mm_per_lai = 2.0;          // a big store, so E is not
        lsdtt::lsdtt_canopy_interception_cpu c({5.0}, {1.0}, p);   // storage-limited
        const auto w = c.step(0.02, pet);

        check("repartition: at factor 2, canopy evaporation reaches 2 x PET",
              close_rel(w.canopy_evaporation, 2.0 * pet, 1e-12),
              "E=" + std::to_string(w.canopy_evaporation) +
              " vs PET=" + std::to_string(pet));
        check("repartition: and transpiration demand is fully consumed",
              close_rel(w.residual_pet, 0.0, 1e-12) || w.residual_pet < 1e-15,
              "residual=" + std::to_string(w.residual_pet));

        // The headline consequence, stated as a number.
        const double total_et = w.canopy_evaporation + w.residual_pet;
        check("repartition: TOTAL ET legitimately EXCEEDS PET",
              total_et > pet * 1.5,
              "total ET " + std::to_string(total_et) +
              " against PET " + std::to_string(pet) +
              " (ratio " + std::to_string(total_et / pet) + ")");
    }

    // The case above is DEMAND-saturated: residual clamps to zero either way,
    // so it cannot tell PET - E/factor from PET - E. Measured: a mutation
    // dropping the division survived it. The discriminating case is a
    // STORAGE-limited canopy, where E < PET*factor and both forms stay
    // positive.
    {
        lsdtt::canopy_parameters p;
        p.interception_evap_factor  = 2.0;
        p.canopy_storage_mm_per_lai = 0.1;      // a small store
        lsdtt::lsdtt_canopy_interception_cpu c({5.0}, {1.0}, p);
        // Capacity is 0.5 mm; wet it fully, then ask for more than it holds.
        const auto w = c.step(0.02, pet);

        const double cap = 0.5e-3;
        check("repartition: the storage-limited case really is storage-limited",
              close_rel(w.canopy_evaporation, cap, 1e-12) &&
              w.canopy_evaporation < pet * p.interception_evap_factor,
              "E=" + std::to_string(w.canopy_evaporation) +
              " capacity=" + std::to_string(cap) +
              " demand=" + std::to_string(pet * p.interception_evap_factor));

        const double want = pet - w.canopy_evaporation / p.interception_evap_factor;
        check("repartition: residual is PET - E/factor, NOT PET - E",
              close_rel(w.residual_pet, want, 1e-12),
              "residual " + std::to_string(w.residual_pet) +
              ", PET-E/f " + std::to_string(want) +
              ", PET-E " + std::to_string(pet - w.canopy_evaporation));
        check("repartition: and the two forms are genuinely distinguishable here",
              std::fabs(want - (pet - w.canopy_evaporation)) > 1e-9,
              "they differ by " +
              std::to_string(want - (pet - w.canopy_evaporation)));
    }

    // Free throughfall must actually bypass the canopy. With a store big enough
    // to swallow the whole event, throughfall is exactly p_free * rainfall --
    // and a mutation setting direct throughfall to zero gives zero instead.
    {
        lsdtt::canopy_parameters p;
        p.lai_gap_coefficient       = 0.6;
        p.canopy_storage_mm_per_lai = 5.0;      // 25 mm at LAI 5: nothing drains
        lsdtt::lsdtt_canopy_interception_cpu c({1.0}, {1.0}, p);
        const double rain = 0.001;
        const auto r = c.step(rain, 0.0);
        const double p_free = std::exp(-0.6 * 1.0);
        check("repartition: free throughfall bypasses the canopy exactly",
              close_rel(r.throughfall, p_free * rain, 1e-12),
              "throughfall " + std::to_string(r.throughfall) +
              " want " + std::to_string(p_free * rain));
        check("repartition: and the rest is held, not passed",
              close_rel(r.canopy_storage, (1.0 - p_free) * rain, 1e-12));
    }

    // A dry canopy must pass the entire demand through untouched.
    {
        lsdtt::canopy_parameters p;
        lsdtt::lsdtt_canopy_interception_cpu c({4.0}, {1.0}, p);
        const auto d = c.step(0.0, pet);
        check("repartition: a dry canopy passes PET through unchanged",
              close_rel(d.residual_pet, pet, 1e-12) &&
              d.canopy_evaporation == 0.0);
    }

    // Bare ground is a pass-through in both directions.
    {
        lsdtt::canopy_parameters p;
        lsdtt::lsdtt_canopy_interception_cpu c({0.0}, {1.0}, p);
        const auto r = c.step(0.005, pet);
        check("repartition: bare ground intercepts nothing at all",
              close_rel(r.throughfall, 0.005, 1e-15) &&
              r.canopy_evaporation == 0.0 &&
              close_rel(r.residual_pet, pet, 1e-15));
    }
}

// =============================================================================
// 4. Seasonality -- at the magnitude that matters
// =============================================================================
void test_seasonality()
{
    const std::vector<double> lai = {5.0};   // summer deciduous woodland

    auto intercepted = [&](double season) {
        lsdtt::canopy_parameters p;
        p.lai_seasonal_factor = season;
        lsdtt::lsdtt_canopy_interception_cpu c(lai, {1.0}, p);
        double e = 0.0;
        for (int t = 0; t < 200; ++t)
            e += c.step((t % 10) < 4 ? 0.002 : 0.0, 0.0003).canopy_evaporation;
        return e;
    };

    const double summer = intercepted(1.0);
    const double winter = intercepted(0.15);   // deciduous LAI ~0.75 in winter

    check("season: winter intercepts far less than summer",
          winter < 0.5 * summer,
          "summer " + std::to_string(summer) + " m, winter " +
          std::to_string(winter) + " m (ratio " +
          std::to_string(winter / summer) + ")");

    // The capacity itself scales exactly with the factor, which is what makes
    // the effect an order of magnitude rather than a nuance.
    lsdtt::canopy_parameters ps, pw;
    pw.lai_seasonal_factor = 0.15;
    lsdtt::lsdtt_canopy_interception_cpu cs(lai, {1.0}, ps);
    lsdtt::lsdtt_canopy_interception_cpu cw(lai, {1.0}, pw);
    check("season: capacity scales linearly with the seasonal factor",
          close_rel(cw.capacity()[0], 0.15 * cs.capacity()[0], 1e-12),
          "summer C " + std::to_string(cs.capacity()[0]) +
          " m, winter C " + std::to_string(cw.capacity()[0]) + " m");
    check("season: and the gap fraction opens up in winter",
          cw.free_throughfall()[0] > cs.free_throughfall()[0],
          "summer p " + std::to_string(cs.free_throughfall()[0]) +
          ", winter p " + std::to_string(cw.free_throughfall()[0]));
}

// =============================================================================
// 5. Land cover -> LAI
// =============================================================================
void test_landcover_ingest()
{
    const std::vector<int> classes = {1, 2, 1, 3, 99};
    const std::unordered_map<int, double> table = {{1, 1.5}, {2, 5.0}, {3, 0.2}};

    const auto r = lsdtt::lai_from_landcover_classes(classes, table, 0.0);
    check("ingest: known classes map to their table LAI",
          r.lai[0] == 1.5 && r.lai[1] == 5.0 && r.lai[3] == 0.2);
    check("ingest: an unknown class takes the default and is COUNTED",
          r.lai[4] == 0.0 && r.n_unmatched == 1,
          std::to_string(r.n_unmatched) + " unmatched");

    const auto r2 = lsdtt::lai_from_landcover_classes(classes, table, 2.0);
    check("ingest: the default is used verbatim", r2.lai[4] == 2.0);

    bool threw = false;
    try { (void)lsdtt::lai_from_landcover_classes(classes, {{1, -1.0}}, 0.0); }
    catch (const std::invalid_argument&) { threw = true; }
    check("ingest: a negative LAI in the table throws", threw);
}

// =============================================================================
// 6. Composed: canopy + column, one balance across BOTH stores
// =============================================================================
void test_composed_balance()
{
    const int nb = 30;
    std::vector<double> idx(nb), ar(nb), lai(nb);
    double s = 0.0;
    for (int i = 0; i < nb; ++i)
    {
        idx[static_cast<std::size_t>(i)] = 12.0 - 10.0 * i / (nb - 1);
        const double z = (idx[static_cast<std::size_t>(i)] - 7.0) / 1.6;
        ar[static_cast<std::size_t>(i)] = std::exp(-0.5 * z * z);
        s += ar[static_cast<std::size_t>(i)];
        // Wetter units carry more vegetation, so LAI varies across the domain
        // and the per-unit forcing path is genuinely exercised.
        lai[static_cast<std::size_t>(i)] = 0.5 + 5.0 * i / (nb - 1);
    }
    for (double& x : ar) x /= s;

    lsdtt::topmodel_parameters tp;
    tp.dt_hours = 1.0;
    lsdtt::lsdtt_topmodel_index_cpu col(
        lsdtt::make_units_from_histogram(ar, idx), tp,
        lsdtt::topmodel_excess_rule::per_unit);

    lsdtt::canopy_parameters cp;
    cp.interception_evap_factor = 1.3;
    lsdtt::lsdtt_canopy_interception_cpu can(lai, ar, cp);

    double total_canopy_evap = 0.0;
    std::vector<double> tf, rp;
    for (int t = 0; t < 3000; ++t)
    {
        const bool storm = (t % 91) < 15;
        const auto cr = can.step(storm ? 0.004 : 0.0, 0.00025, tf, rp);
        total_canopy_evap += cr.canopy_evaporation;
        (void)col.step(tf, rp);
    }

    const auto& cb = can.water_balance();
    const auto& tb = col.water_balance();

    // The column's "rainfall" is the canopy's throughfall, by construction.
    check("composed: the column received exactly the canopy's throughfall",
          close_rel(tb.total_rainfall, cb.total_throughfall, 1e-9),
          "column P " + std::to_string(tb.total_rainfall) +
          " canopy TF " + std::to_string(cb.total_throughfall));

    // One balance over both stores.
    const double storage_change =
        (cb.final_storage - cb.initial_storage) +
        (tb.final_storage - tb.initial_storage);
    const double net_input = cb.total_rainfall
                           - cb.total_canopy_evaporation
                           - tb.total_actual_et
                           - tb.total_runoff
                           - tb.truncation_loss;
    const double err = std::fabs(storage_change - net_input);
    const double scale = std::max(cb.total_rainfall, 1e-12);

    check("composed: the two-store balance closes",
          err / scale < 1e-10,
          "relative error " + std::to_string(err / scale) +
          ", P=" + std::to_string(cb.total_rainfall) +
          " Ecan=" + std::to_string(cb.total_canopy_evaporation) +
          " Ecol=" + std::to_string(tb.total_actual_et) +
          " Q=" + std::to_string(tb.total_runoff));

    check("composed: interception actually removed water",
          total_canopy_evap > 0.0,
          "canopy evaporation " + std::to_string(total_canopy_evap) + " m");

    // Per-unit forcing must reach the column: a length mismatch has to throw
    // rather than silently applying one unit's forcing to another.
    bool threw = false;
    try { std::vector<double> bad(5, 0.0); (void)col.step(bad, bad); }
    catch (const std::invalid_argument&) { threw = true; }
    check("composed: mismatched per-unit forcing lengths throw", threw);
}

// =============================================================================
// 7. Validation
// =============================================================================
void test_validation()
{
    auto throws = [](auto&& fn) {
        try { fn(); } catch (const std::invalid_argument&) { return true; }
        catch (...) { return false; }
        return false;
    };
    lsdtt::canopy_parameters p;

    check("validation: an empty LAI vector throws",
          throws([&]{ lsdtt::lsdtt_canopy_interception_cpu c({}, {}, p); }));
    check("validation: mismatched LAI and area lengths throw",
          throws([&]{ lsdtt::lsdtt_canopy_interception_cpu c({1.0, 2.0}, {1.0}, p); }));
    check("validation: areas not summing to 1 throw",
          throws([&]{ lsdtt::lsdtt_canopy_interception_cpu c({1.0, 2.0}, {0.3, 0.3}, p); }));
    check("validation: a negative LAI throws",
          throws([&]{ lsdtt::lsdtt_canopy_interception_cpu c({-1.0}, {1.0}, p); }));
    check("validation: a non-positive seasonal factor throws",
          throws([&]{ auto q = p; q.lai_seasonal_factor = 0.0;
                      lsdtt::lsdtt_canopy_interception_cpu c({1.0}, {1.0}, q); }));
    check("validation: a non-positive evaporation factor throws",
          throws([&]{ auto q = p; q.interception_evap_factor = 0.0;
                      lsdtt::lsdtt_canopy_interception_cpu c({1.0}, {1.0}, q); }));
}

} // namespace

int main()
{
    std::cout << "lsdtt-test-canopy-interception\n"
              << "------------------------------\n";
    try
    {
        test_lai_forms();
        test_balance();
        test_evaporative_repartition();
        test_seasonality();
        test_landcover_ingest();
        test_composed_balance();
        test_validation();
    }
    catch (const std::exception& e)
    {
        std::cout << "FAIL  unexpected exception: " << e.what() << std::endl;
        ++g_n_failed; ++g_n_checks;
    }

    std::cout << "------------------------------\n"
              << g_n_checks - g_n_failed << " / " << g_n_checks << " checks passed\n"
              << "RESULT: " << (g_n_failed == 0 ? "PASS" : "FAIL") << std::endl;
    return g_n_failed == 0 ? 0 : 1;
}
