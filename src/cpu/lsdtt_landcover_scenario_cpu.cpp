// lsdtt_landcover_scenario_cpu.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html

#include "lsdtt_landcover_scenario_cpu.hpp"

#include <algorithm>
#include <cmath>
#include <sstream>
#include <stdexcept>
#include <utility>

#include "lsdtt_topmodel_cpu.hpp"
#include "../utils/lsdtt_array2d.hpp"

namespace lsdtt
{

namespace
{

std::string trim(const std::string& s)
{
    const auto a = s.find_first_not_of(" \t\r\n");
    if (a == std::string::npos) return "";
    const auto b = s.find_last_not_of(" \t\r\n");
    return s.substr(a, b - a + 1);
}

std::vector<std::string> split(const std::string& s, char sep)
{
    std::vector<std::string> out;
    std::string cell;
    std::istringstream is(s);
    while (std::getline(is, cell, sep))
    {
        const std::string t = trim(cell);
        if (!t.empty()) out.push_back(t);
    }
    return out;
}

void require_same_grid(const lsdtt_raster_cpu& r, const lsdtt_raster_meta& meta,
                       const char* name, const char* who)
{
    const lsdtt_raster_meta& m = r.get_meta();
    if (m.n_rows != meta.n_rows || m.n_cols != meta.n_cols)
        throw std::invalid_argument(
            std::string(who) + ": " + name + " is " + std::to_string(m.n_rows) +
            "x" + std::to_string(m.n_cols) + " but the land-cover grid is " +
            std::to_string(meta.n_rows) + "x" + std::to_string(meta.n_cols));
}

std::string series_name_for(
    std::uint32_t code,
    const std::vector<std::pair<std::uint32_t, std::string>>& lookup)
{
    for (const auto& kv : lookup) if (kv.first == code) return kv.second;
    return "";
}

} // namespace

landcover_table landcover_table::parse(const std::string& spec)
{
    const char* who = "landcover_table::parse";
    landcover_table out;

    for (const std::string& entry : split(spec, ';'))
    {
        const auto colon = entry.find(':');
        if (colon == std::string::npos)
            throw std::invalid_argument(
                std::string(who) + ": entry '" + entry +
                "' has no ':' separating the class name from its fields");

        landcover_class c;
        c.name = trim(entry.substr(0, colon));
        if (c.name.empty())
            throw std::invalid_argument(
                std::string(who) + ": entry '" + entry + "' has an empty name");

        bool have_code = false, have_lai = false, have_n = false;
        for (const std::string& kv : split(entry.substr(colon + 1), ','))
        {
            const auto eq = kv.find('=');
            if (eq == std::string::npos)
                throw std::invalid_argument(
                    std::string(who) + ": field '" + kv + "' in class '" +
                    c.name + "' is not key=value");
            const std::string key = trim(kv.substr(0, eq));
            const std::string val = trim(kv.substr(eq + 1));

            try
            {
                if      (key == "code") { c.code = std::stoi(val); have_code = true; }
                else if (key == "lai")  { c.lai  = std::stod(val); have_lai  = true; }
                else if (key == "n")    { c.manning_n = std::stod(val); have_n = true; }
                else if (key == "soil") { c.soil_land_use = val; }
                else if (key == "ks" || key == "kv")
                    throw std::invalid_argument(
                        std::string(who) + ": class '" + c.name + "' sets '" + key +
                        "'. Conductivity is NOT a land-cover property here: it "
                        "comes from the measured soil lookup, re-keyed by land "
                        "use. Remove it and set 'soil=' instead.");
                else
                    throw std::invalid_argument(
                        std::string(who) + ": class '" + c.name +
                        "' has unknown field '" + key + "'");
            }
            catch (const std::out_of_range&)
            {
                throw std::invalid_argument(
                    std::string(who) + ": value '" + val + "' for '" + key +
                    "' in class '" + c.name + "' is out of range");
            }
            catch (const std::invalid_argument&)
            {
                // std::stoi/stod also throw invalid_argument on garbage; keep
                // our own messages but make a parse failure legible.
                if (key == "code" || key == "lai" || key == "n")
                    throw std::invalid_argument(
                        std::string(who) + ": value '" + val + "' for '" + key +
                        "' in class '" + c.name + "' is not a number");
                throw;
            }
        }

        if (!have_code || !have_lai || !have_n)
            throw std::invalid_argument(
                std::string(who) + ": class '" + c.name +
                "' needs code=, lai= and n=");
        if (!(c.lai >= 0.0) || !std::isfinite(c.lai))
            throw std::invalid_argument(
                std::string(who) + ": class '" + c.name + "' has a negative or "
                "non-finite LAI");
        if (!(c.manning_n > 0.0) || !std::isfinite(c.manning_n))
            throw std::invalid_argument(
                std::string(who) + ": class '" + c.name + "' has a non-positive "
                "or non-finite roughness");

        for (const landcover_class& e : out.classes)
        {
            if (e.code == c.code)
                throw std::invalid_argument(
                    std::string(who) + ": code " + std::to_string(c.code) +
                    " is used by both '" + e.name + "' and '" + c.name + "'");
            if (e.name == c.name)
                throw std::invalid_argument(
                    std::string(who) + ": class '" + c.name + "' is defined twice");
        }
        out.classes.push_back(std::move(c));
    }

    if (out.classes.empty())
        throw std::invalid_argument(std::string(who) + ": no classes defined");
    return out;
}

const landcover_class* landcover_table::by_code(int code) const noexcept
{
    for (const landcover_class& c : classes) if (c.code == code) return &c;
    return nullptr;
}

const landcover_class* landcover_table::by_name(const std::string& name) const noexcept
{
    for (const landcover_class& c : classes) if (c.name == name) return &c;
    return nullptr;
}

landcover_scenario landcover_scenario::parse(const std::string& spec)
{
    const char* who = "landcover_scenario::parse";
    landcover_scenario out;

    for (const std::string& rule : split(spec, ';'))
    {
        const auto arrow = rule.find("->");
        if (arrow == std::string::npos)
            throw std::invalid_argument(
                std::string(who) + ": rule '" + rule + "' has no '->'");
        landcover_scenario_rule r;
        r.from = trim(rule.substr(0, arrow));
        r.to   = trim(rule.substr(arrow + 2));
        if (r.from.empty() || r.to.empty())
            throw std::invalid_argument(
                std::string(who) + ": rule '" + rule + "' has an empty side");
        out.rules.push_back(std::move(r));
    }
    return out;
}

scenario_application_result apply_landcover_scenario(
    const lsdtt_raster_cpu&   landcover,
    const lsdtt_raster_cpu*   intervention_mask,
    const landcover_table&    table,
    const landcover_scenario& scenario)
{
    const char* who = "apply_landcover_scenario";

    const lsdtt_raster_meta& meta = landcover.get_meta();
    if (intervention_mask)
        require_same_grid(*intervention_mask, meta, "intervention_mask", who);

    // Resolve the rules against the table up front, so an unknown class name is
    // an error at the start rather than a surprise part-way through a raster.
    struct resolved { int from; int to; std::string label; };
    std::vector<resolved> rules;
    for (const landcover_scenario_rule& r : scenario.rules)
    {
        const landcover_class* to = table.by_name(r.to);
        if (to == nullptr)
            throw std::invalid_argument(
                std::string(who) + ": rule target '" + r.to +
                "' is not a class in the table");
        int from_code = -1;
        if (r.from != "*")
        {
            const landcover_class* f = table.by_name(r.from);
            if (f == nullptr)
                throw std::invalid_argument(
                    std::string(who) + ": rule source '" + r.from +
                    "' is not a class in the table");
            from_code = f->code;
        }
        rules.push_back({from_code, to->code, r.from + "->" + r.to});
    }

    const int   n_rows = meta.n_rows;
    const int   n_cols = meta.n_cols;
    const float ndv    = static_cast<float>(meta.no_data_value);

    const auto& lc_d = landcover.get_data_array();
    const auto& lc_m = landcover.get_mask_array();

    lsdtt_array2d_f  out_d(n_rows, n_cols, ndv);
    lsdtt_array2d_u8 out_m(n_rows, n_cols, std::uint8_t(0));

    scenario_application_result out;

    for (int row = 0; row < n_rows; ++row)
        for (int col = 0; col < n_cols; ++col)
        {
            if (lc_m(row, col) == 0) continue;

            const int code = static_cast<int>(std::lround(lc_d(row, col)));
            if (table.by_code(code) == nullptr)
                throw std::invalid_argument(
                    std::string(who) + ": the raster contains class code " +
                    std::to_string(code) + ", which the table does not define");

            out_d(row, col) = static_cast<float>(code);
            out_m(row, col) = std::uint8_t(1);

            bool in_mask = true;
            if (intervention_mask)
            {
                in_mask = intervention_mask->get_mask_array()(row, col) != 0 &&
                          intervention_mask->get_data_array()(row, col) != 0.0f;
            }
            if (!in_mask) continue;

            ++out.n_considered;

            bool matched = false;
            for (const resolved& r : rules)
            {
                if (r.from != -1 && r.from != code) continue;
                // First matching rule wins, so an earlier specific rule beats a
                // later catch-all.
                out_d(row, col) = static_cast<float>(r.to);
                ++out.transitions[r.label];
                if (r.to != code) ++out.n_reclassified;
                matched = true;
                break;
            }
            if (!matched) ++out.n_unmatched;
        }

    out.landcover = lsdtt_raster_cpu(std::move(out_d), std::move(out_m), meta);
    return out;
}

landcover_field_result derive_landcover_fields(
    const landcover_field_inputs& inputs,
    const landcover_table&        table)
{
    const char* who = "derive_landcover_fields";

    if (inputs.landcover == nullptr)
        throw std::invalid_argument(std::string(who) + ": landcover is required");

    const bool want_t0 = (inputs.soil_series_codes != nullptr);
    if (want_t0)
    {
        if (inputs.soil_series_lookup == nullptr || inputs.lnkd_table == nullptr)
            throw std::invalid_argument(
                std::string(who) + ": soil_series_codes was supplied, so "
                "soil_series_lookup and lnkd_table are required too");
    }

    const lsdtt_raster_meta& meta = inputs.landcover->get_meta();
    if (want_t0) require_same_grid(*inputs.soil_series_codes, meta,
                                   "soil_series_codes", who);

    const int   n_rows = meta.n_rows;
    const int   n_cols = meta.n_cols;
    const float ndv    = static_cast<float>(meta.no_data_value);

    const auto& lc_d = inputs.landcover->get_data_array();
    const auto& lc_m = inputs.landcover->get_mask_array();

    lsdtt_array2d_f  lai_d(n_rows, n_cols, ndv);
    lsdtt_array2d_u8 lai_m(n_rows, n_cols, std::uint8_t(0));
    lsdtt_array2d_f  n_d(n_rows, n_cols, ndv);
    lsdtt_array2d_u8 n_m(n_rows, n_cols, std::uint8_t(0));
    lsdtt_array2d_f  t0_d(n_rows, n_cols, ndv);
    lsdtt_array2d_u8 t0_m(n_rows, n_cols, std::uint8_t(0));

    landcover_field_result out;

    for (int row = 0; row < n_rows; ++row)
        for (int col = 0; col < n_cols; ++col)
        {
            if (lc_m(row, col) == 0) continue;

            const int code = static_cast<int>(std::lround(lc_d(row, col)));
            const landcover_class* c = table.by_code(code);
            if (c == nullptr)
                throw std::invalid_argument(
                    std::string(who) + ": the raster contains class code " +
                    std::to_string(code) + ", which the table does not define");

            ++out.n_cells;

            lai_d(row, col) = static_cast<float>(c->lai);
            lai_m(row, col) = std::uint8_t(1);
            n_d(row, col)   = static_cast<float>(c->manning_n);
            n_m(row, col)   = std::uint8_t(1);

            if (!want_t0) continue;

            if (inputs.soil_series_codes->get_mask_array()(row, col) == 0)
            { ++out.n_t0_missing; continue; }

            const std::uint32_t scode = static_cast<std::uint32_t>(std::lround(
                inputs.soil_series_codes->get_data_array()(row, col)));
            const std::string series =
                series_name_for(scode, *inputs.soil_series_lookup);
            if (series.empty()) { ++out.n_t0_missing; continue; }

            // THE RE-KEYING. The land use comes from the (possibly
            // reclassified) class, so a scenario changes transmissivity by
            // selecting a different MEASURED profile for the same soil -- not by
            // scaling the old value.
            double ln_kd = 0.0;
            std::string source;
            int n_profiles = 0;
            const bool found = lookup_ln_kd_with_fallback(
                *inputs.lnkd_table, series, c->soil_land_use,
                ln_kd, source, n_profiles);

            if (!found) { ++out.n_t0_missing; continue; }
            if (source == "direct") ++out.n_t0_direct; else ++out.n_t0_fallback;

            // ln_kd is ln of m.cm/day; convert to ln(T0) in m^2/h so the value
            // can be used as topmodel_parameters::ln_t0 directly.
            t0_d(row, col) = static_cast<float>(
                sti_ln_kd_to_ln_t0_m2_per_hour(ln_kd));
            t0_m(row, col) = std::uint8_t(1);
        }

    out.lai       = lsdtt_raster_cpu(std::move(lai_d), std::move(lai_m), meta);
    out.roughness = lsdtt_raster_cpu(std::move(n_d),   std::move(n_m),   meta);
    if (want_t0)
        out.ln_t0 = lsdtt_raster_cpu(std::move(t0_d), std::move(t0_m), meta);

    return out;
}

} // namespace lsdtt
