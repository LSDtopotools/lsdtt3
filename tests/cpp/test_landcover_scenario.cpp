// test_landcover_scenario.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
//
// -----------------------------------------------------------------------------
// Unit tests for src/cpu/lsdtt_landcover_scenario_cpu -- Phase 9.
//
// THE CENTREPIECE: transmissivity is RE-KEYED, not multiplied
// -----------------------------------------------------------
// test_transmissivity_is_rekeyed is the check that matters. SD-TOPMODEL
// represents woodland planting as `Ks x 1.5` -- one multiplier applied
// everywhere. Here, reclassifying a pixel re-keys the MEASURED (soil series,
// land use) lookup and returns the transmissivity that soil actually shows
// under the new cover.
//
// The difference is testable, and the test is built so a multiplier CANNOT
// masquerade as re-keying: the fixture gives two soil series whose measured
// response to the same land-cover change differs in RATIO. Any implementation
// that scales the old value produces one ratio for both, so the check fails.
// One of the two series is also given a response in the opposite DIRECTION,
// which no positive multiplier can reproduce at all.
//
// Usage: lsdtt-test-landcover-scenario
// -----------------------------------------------------------------------------

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "lsdtt_landcover_scenario_cpu.hpp"
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

std::string fx(double v)
{
    std::ostringstream o; o << std::fixed << std::setprecision(4) << v;
    return o.str();
}

constexpr double k_nodata = -9999.0;
constexpr double k_float_tol = 1e-6;   // rasters store float

lsdtt::lsdtt_raster_meta make_meta(int nr, int nc)
{
    lsdtt::lsdtt_raster_meta m;
    m.n_rows = nr; m.n_cols = nc; m.cell_size = 10.0;
    m.no_data_value = k_nodata;
    m.set_default_synthetic_georeference();
    return m;
}

lsdtt::lsdtt_raster_cpu make_raster(const lsdtt::lsdtt_raster_meta& meta,
                                    const std::vector<float>& v)
{
    lsdtt::lsdtt_array2d_f  d(meta.n_rows, meta.n_cols, static_cast<float>(k_nodata));
    lsdtt::lsdtt_array2d_u8 m(meta.n_rows, meta.n_cols, std::uint8_t(0));
    for (int r = 0; r < meta.n_rows; ++r)
        for (int c = 0; c < meta.n_cols; ++c)
        {
            const float x = v[static_cast<std::size_t>(r) * meta.n_cols + c];
            d(r, c) = x;
            m(r, c) = (x == static_cast<float>(k_nodata)) ? std::uint8_t(0)
                                                          : std::uint8_t(1);
        }
    return lsdtt::lsdtt_raster_cpu(std::move(d), std::move(m), meta);
}

const char* k_table_spec =
    "grazing:code=1,lai=1.5,n=0.30,soil=CULT; "
    "woodland:code=2,lai=5.0,n=0.60,soil=SEMI; "
    "hedgerow:code=3,lai=3.0,n=0.45,soil=SEMI";

// =============================================================================
// 1. Parsing
// =============================================================================
void test_parsing()
{
    const auto t = lsdtt::landcover_table::parse(k_table_spec);
    check("parse: all three classes are read", t.size() == 3,
          std::to_string(t.size()) + " classes");
    check("parse: fields are read correctly",
          t.by_name("woodland") != nullptr &&
          t.by_name("woodland")->code == 2 &&
          close_rel(t.by_name("woodland")->lai, 5.0, 1e-12) &&
          close_rel(t.by_name("woodland")->manning_n, 0.60, 1e-12) &&
          t.by_name("woodland")->soil_land_use == "SEMI");
    check("parse: lookup by code and by name agree",
          t.by_code(1) == t.by_name("grazing"));

    const auto sc = lsdtt::landcover_scenario::parse("grazing->woodland; hedgerow->woodland");
    check("parse: scenario rules are read",
          sc.rules.size() == 2 && sc.rules[0].from == "grazing" &&
          sc.rules[0].to == "woodland");

    auto throws = [](auto&& fn) {
        try { fn(); } catch (const std::invalid_argument&) { return true; }
        catch (...) { return false; }
        return false;
    };

    // The one that matters for this phase: conductivity is NOT a land-cover
    // property, and saying so must be an error rather than silently ignored.
    //
    // The MESSAGE is checked, not just that something threw. Deleting the
    // dedicated branch leaves `ks=` falling through to the generic "unknown
    // field" error, which still throws and so still passed a bare throws()
    // check -- but tells a user migrating from SD-TOPMODEL nothing about why
    // their parameter is not accepted or what to write instead.
    auto throws_saying = [](const std::string& spec,
                            const std::string& must_contain) {
        try { (void)lsdtt::landcover_table::parse(spec); }
        catch (const std::invalid_argument& e)
        {
            return std::string(e.what()).find(must_contain) != std::string::npos;
        }
        catch (...) { return false; }
        return false;
    };
    check("parse: ks= is rejected AND the error explains where Ks comes from",
          throws_saying("grazing:code=1,lai=1.5,n=0.3,ks=0.004",
                        "measured soil lookup"));
    check("parse: kv= -- SD-TOPMODEL's name for it -- is rejected the same way",
          throws_saying("grazing:code=1,lai=1.5,n=0.3,kv=25",
                        "measured soil lookup"));
    check("parse: and the message says what to write instead",
          throws_saying("grazing:code=1,lai=1.5,n=0.3,ks=0.004", "soil="));

    check("parse: a duplicate code throws",
          throws([]{ (void)lsdtt::landcover_table::parse(
              "a:code=1,lai=1,n=0.3; b:code=1,lai=2,n=0.4"); }));
    check("parse: a missing required field throws",
          throws([]{ (void)lsdtt::landcover_table::parse("a:code=1,lai=1"); }));
    check("parse: a negative LAI throws",
          throws([]{ (void)lsdtt::landcover_table::parse(
              "a:code=1,lai=-1,n=0.3"); }));
    check("parse: an unknown field throws",
          throws([]{ (void)lsdtt::landcover_table::parse(
              "a:code=1,lai=1,n=0.3,colour=green"); }));
    check("parse: a rule without '->' throws",
          throws([]{ (void)lsdtt::landcover_scenario::parse("grazing woodland"); }));
}

// =============================================================================
// 2. Applying a scenario
// =============================================================================
void test_apply()
{
    const auto table = lsdtt::landcover_table::parse(k_table_spec);
    const auto meta  = make_meta(2, 4);
    //                                 grazing grazing woodland grazing
    const auto lc = make_raster(meta, {1.f, 1.f, 2.f, 1.f,
                                       1.f, 1.f, 2.f, 1.f});
    // Intervention only on the left half.
    const auto mask = make_raster(meta, {1.f, 1.f, 0.f, 0.f,
                                         1.f, 1.f, 0.f, 0.f});

    const auto sc = lsdtt::landcover_scenario::parse("grazing->woodland");
    const auto r = lsdtt::apply_landcover_scenario(lc, &mask, table, sc);

    check("apply: only cells inside the intervention are reclassified",
          r.n_reclassified == 4 && r.n_considered == 4,
          std::to_string(r.n_reclassified) + " of " +
          std::to_string(r.n_considered) + " considered");
    check("apply: the reclassified cells carry the new code",
          r.landcover.get_data_array()(0, 0) == 2.f &&
          r.landcover.get_data_array()(1, 1) == 2.f);
    check("apply: cells outside the intervention are untouched",
          r.landcover.get_data_array()(0, 3) == 1.f &&
          r.landcover.get_data_array()(0, 2) == 2.f);
    check("apply: the transition is counted and labelled",
          r.transitions.count("grazing->woodland") == 1 &&
          r.transitions.at("grazing->woodland") == 4);

    // A class matching no rule is counted, not silently passed over.
    {
        const auto sc2 = lsdtt::landcover_scenario::parse("hedgerow->woodland");
        const auto r2 = lsdtt::apply_landcover_scenario(lc, &mask, table, sc2);
        check("apply: cells matching no rule are COUNTED as unmatched",
              r2.n_unmatched == 4 && r2.n_reclassified == 0,
              std::to_string(r2.n_unmatched) + " unmatched");
    }

    // A wildcard rule, and first-match-wins ordering.
    {
        const auto sc3 = lsdtt::landcover_scenario::parse("grazing->woodland; *->hedgerow");
        const auto r3 = lsdtt::apply_landcover_scenario(lc, nullptr, table, sc3);
        check("apply: an earlier specific rule beats a later catch-all",
              r3.landcover.get_data_array()(0, 0) == 2.f &&   // grazing -> woodland
              r3.landcover.get_data_array()(0, 2) == 3.f,     // woodland -> hedgerow
              "grazing became " +
              std::to_string(static_cast<int>(r3.landcover.get_data_array()(0, 0))) +
              ", woodland became " +
              std::to_string(static_cast<int>(r3.landcover.get_data_array()(0, 2))));
    }

    // No mask means everywhere.
    {
        const auto r4 = lsdtt::apply_landcover_scenario(lc, nullptr, table, sc);
        check("apply: with no intervention mask the scenario applies everywhere",
              r4.n_considered == 8 && r4.n_reclassified == 6,
              std::to_string(r4.n_reclassified) + " reclassified of 8");
    }

    auto throws = [](auto&& fn) {
        try { fn(); } catch (const std::invalid_argument&) { return true; }
        catch (...) { return false; }
        return false;
    };
    check("apply: an unknown class code in the raster throws",
          throws([&]{ const auto bad = make_raster(meta, {9.f,9.f,9.f,9.f,
                                                          9.f,9.f,9.f,9.f});
                      (void)lsdtt::apply_landcover_scenario(bad, nullptr, table, sc); }));
    check("apply: a rule naming an unknown class throws",
          throws([&]{ const auto s = lsdtt::landcover_scenario::parse("grazing->forest");
                      (void)lsdtt::apply_landcover_scenario(lc, nullptr, table, s); }));
}

// =============================================================================
// 3. THE CENTREPIECE: re-keying, not multiplying
// =============================================================================
void test_transmissivity_is_rekeyed()
{
    const auto table = lsdtt::landcover_table::parse(k_table_spec);
    const auto meta  = make_meta(1, 2);

    // Two cells, two DIFFERENT soil series, both starting as grazing (CULT).
    const auto lc   = make_raster(meta, {1.f, 1.f});
    const auto soil = make_raster(meta, {10.f, 20.f});
    const std::vector<std::pair<std::uint32_t, std::string>> soil_lookup = {
        {10u, "PEATYGLEY"}, {20u, "BROWNEARTH"}
    };

    // A measured lookup in which the two series respond to the SAME land-cover
    // change by different ratios -- and in opposite directions. This is what
    // makes the test able to tell re-keying from multiplication: no single
    // positive multiplier can produce both.
    lsdtt::sti_lnkd_table lnkd;
    lnkd.ln_kd_by_key["PEATYGLEY|CULT"]  = std::log(2.0);
    lnkd.ln_kd_by_key["PEATYGLEY|SEMI"]  = std::log(8.0);    // x4 -- planting helps a lot
    lnkd.ln_kd_by_key["BROWNEARTH|CULT"] = std::log(50.0);
    lnkd.ln_kd_by_key["BROWNEARTH|SEMI"] = std::log(25.0);   // x0.5 -- the OTHER way

    lsdtt::landcover_field_inputs in;
    in.soil_series_codes  = &soil;
    in.soil_series_lookup = &soil_lookup;
    in.lnkd_table         = &lnkd;

    in.landcover = &lc;
    const auto before = lsdtt::derive_landcover_fields(in, table);

    const auto sc  = lsdtt::landcover_scenario::parse("grazing->woodland");
    const auto app = lsdtt::apply_landcover_scenario(lc, nullptr, table, sc);
    in.landcover = &app.landcover;
    const auto after = lsdtt::derive_landcover_fields(in, table);

    check("rekey: every cell found a direct measured match",
          before.n_t0_direct == 2 && after.n_t0_direct == 2 &&
          before.n_t0_missing == 0 && after.n_t0_missing == 0);

    // T0 itself, recovered from ln(T0).
    auto t0 = [](const lsdtt::lsdtt_raster_cpu& r, int c) {
        return std::exp(static_cast<double>(r.get_data_array()(0, c)));
    };
    const double peat_before = t0(before.ln_t0, 0), peat_after = t0(after.ln_t0, 0);
    const double earth_before = t0(before.ln_t0, 1), earth_after = t0(after.ln_t0, 1);

    const double peat_ratio  = peat_after / peat_before;
    const double earth_ratio = earth_after / earth_before;

    check("rekey: the peaty gley's transmissivity QUADRUPLES under woodland",
          close_rel(peat_ratio, 4.0, 1e-5),
          "ratio " + fx(peat_ratio));
    check("rekey: the brown earth's HALVES -- the opposite direction",
          close_rel(earth_ratio, 0.5, 1e-5),
          "ratio " + fx(earth_ratio));

    // The discriminating statement. A multiplier gives one ratio; measured
    // re-keying gives the ratio each soil actually shows.
    check("REKEY: the two soils respond DIFFERENTLY to the same land-cover change",
          !close_rel(peat_ratio, earth_ratio, 0.1),
          "peaty gley x" + fx(peat_ratio) + ", brown earth x" + fx(earth_ratio) +
          " -- no single multiplier can produce both");
    check("REKEY: and one response is in the opposite direction to the other",
          (peat_ratio > 1.0) != (earth_ratio > 1.0),
          "which no positive multiplier can reproduce at all");

    // And the values are the MEASURED ones, not anything derived from the old
    // value: each equals the table entry for (this series, the new land use).
    check("rekey: the new T0 is exactly the measured (series, land use) value",
          close_rel(peat_after,
                    std::exp(lsdtt::sti_ln_kd_to_ln_t0_m2_per_hour(std::log(8.0))),
                    1e-5) &&
          close_rel(earth_after,
                    std::exp(lsdtt::sti_ln_kd_to_ln_t0_m2_per_hour(std::log(25.0))),
                    1e-5));

    // The unit conversion is applied: ln_kd is ln(m.cm/day), ln_t0 is ln(m2/h).
    check("rekey: the STI units conversion is applied on the way out",
          close_rel(static_cast<double>(after.ln_t0.get_data_array()(0, 0)),
                    lsdtt::sti_ln_kd_to_ln_t0_m2_per_hour(std::log(8.0)), 1e-6),
          "raw ln_kd " + fx(std::log(8.0)) + " -> ln_t0 " +
          fx(after.ln_t0.get_data_array()(0, 0)));

    // A missing pair falls back, and the substitution is counted rather than
    // passed off as measured.
    {
        lsdtt::sti_lnkd_table sparse;
        sparse.ln_kd_by_key["PEATYGLEY|CULT"]  = std::log(2.0);
        sparse.ln_kd_by_key["BROWNEARTH|CULT"] = std::log(50.0);
        // no SEMI entries at all
        lsdtt::landcover_field_inputs in2 = in;
        in2.lnkd_table = &sparse;
        const auto f = lsdtt::derive_landcover_fields(in2, table);
        check("rekey: a missing (series, land use) pair is served by fallback AND counted",
              f.n_t0_fallback == 2 && f.n_t0_direct == 0,
              std::to_string(f.n_t0_fallback) + " fallbacks, " +
              std::to_string(f.n_t0_direct) + " direct");
    }

    // An unknown soil series leaves T0 nodata and is counted, rather than
    // inventing a value.
    {
        const auto unknown = make_raster(meta, {99.f, 99.f});
        lsdtt::landcover_field_inputs in3 = in;
        in3.soil_series_codes = &unknown;
        const auto f = lsdtt::derive_landcover_fields(in3, table);
        check("rekey: an unknown soil series leaves T0 nodata and is counted",
              f.n_t0_missing == 2 &&
              f.ln_t0.get_mask_array()(0, 0) == 0);
    }
}

// =============================================================================
// 4. The three fields move together
// =============================================================================
void test_fields_move_together()
{
    const auto table = lsdtt::landcover_table::parse(k_table_spec);
    const auto meta  = make_meta(1, 2);
    const auto lc    = make_raster(meta, {1.f, 1.f});
    const auto soil  = make_raster(meta, {10.f, 10.f});
    const std::vector<std::pair<std::uint32_t, std::string>> lookup = {{10u, "S"}};

    lsdtt::sti_lnkd_table lnkd;
    lnkd.ln_kd_by_key["S|CULT"] = std::log(2.0);
    lnkd.ln_kd_by_key["S|SEMI"] = std::log(6.0);

    lsdtt::landcover_field_inputs in;
    in.landcover = &lc; in.soil_series_codes = &soil;
    in.soil_series_lookup = &lookup; in.lnkd_table = &lnkd;
    const auto before = lsdtt::derive_landcover_fields(in, table);

    const auto app = lsdtt::apply_landcover_scenario(
        lc, nullptr, table, lsdtt::landcover_scenario::parse("grazing->woodland"));
    in.landcover = &app.landcover;
    const auto after = lsdtt::derive_landcover_fields(in, table);

    check("together: LAI rises with the planting",
          close_rel(before.lai.get_data_array()(0, 0), 1.5, k_float_tol) &&
          close_rel(after.lai.get_data_array()(0, 0), 5.0, k_float_tol),
          "LAI 1.5 -> 5.0");
    check("together: roughness rises with it",
          close_rel(before.roughness.get_data_array()(0, 0), 0.30, k_float_tol) &&
          close_rel(after.roughness.get_data_array()(0, 0), 0.60, k_float_tol),
          "n 0.30 -> 0.60");
    check("together: and transmissivity rises with it",
          after.ln_t0.get_data_array()(0, 0) >
          before.ln_t0.get_data_array()(0, 0),
          "ln(T0) " + fx(before.ln_t0.get_data_array()(0, 0)) + " -> " +
          fx(after.ln_t0.get_data_array()(0, 0)));

    // All three changed: a scenario that moved only one would be a bug.
    const bool all_three =
        !close_rel(before.lai.get_data_array()(0, 0),
                   after.lai.get_data_array()(0, 0), 1e-3) &&
        !close_rel(before.roughness.get_data_array()(0, 0),
                   after.roughness.get_data_array()(0, 0), 1e-3) &&
        !close_rel(before.ln_t0.get_data_array()(0, 0),
                   after.ln_t0.get_data_array()(0, 0), 1e-3);
    check("together: ALL THREE fields moved, not just one", all_three);

    // Without soil inputs the other two still work; ln_t0 is simply absent.
    {
        lsdtt::landcover_field_inputs bare;
        bare.landcover = &lc;
        const auto f = lsdtt::derive_landcover_fields(bare, table);
        check("together: LAI and roughness work with no soil data at all",
              f.n_cells == 2 &&
              close_rel(f.lai.get_data_array()(0, 0), 1.5, k_float_tol));
    }

    auto throws = [](auto&& fn) {
        try { fn(); } catch (const std::invalid_argument&) { return true; }
        catch (...) { return false; }
        return false;
    };
    check("together: soil codes without a lookup table throws",
          throws([&]{ lsdtt::landcover_field_inputs b;
                      b.landcover = &lc; b.soil_series_codes = &soil;
                      (void)lsdtt::derive_landcover_fields(b, table); }));
    check("together: a missing land-cover raster throws",
          throws([&]{ lsdtt::landcover_field_inputs b;
                      (void)lsdtt::derive_landcover_fields(b, table); }));
}

} // namespace

int main()
{
    std::cout << "lsdtt-test-landcover-scenario\n"
              << "-----------------------------\n";
    try
    {
        test_parsing();
        test_apply();
        test_transmissivity_is_rekeyed();
        test_fields_move_together();
    }
    catch (const std::exception& e)
    {
        std::cout << "FAIL  unexpected exception: " << e.what() << std::endl;
        ++g_n_failed; ++g_n_checks;
    }

    std::cout << "-----------------------------\n"
              << g_n_checks - g_n_failed << " / " << g_n_checks << " checks passed\n"
              << "RESULT: " << (g_n_failed == 0 ? "PASS" : "FAIL") << std::endl;
    return g_n_failed == 0 ? 0 : 1;
}
