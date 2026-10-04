// lsdtt_hype_parameters.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// HYPE parameter tables and the sub-basin set-up struct: readers and
// validation. See lsdtt_hype_parameters.hpp for the contract.
//
// The parameter names live in one place each, the k_*_fields tables below,
// which pair a HYPE name with the struct member it fills. The readers, the
// unknown-name check and the public *_parameter_names() lists are all driven
// from those tables, so a parameter added to a struct needs one new line here
// and cannot be accepted by one path and rejected by another.
// -----------------------------------------------------------------------------

#include "lsdtt_hype_parameters.hpp"

#include <climits>
#include <cmath>
#include <filesystem>
#include <sstream>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>
#include <utility>

#include "lsdtt_csv_table.hpp"

namespace lsdtt
{

namespace
{

// ---------------------------------------------------------------------------
// Name -> member tables
// ---------------------------------------------------------------------------
template <typename T>
struct hype_field
{
    const char* name;
    double T::*member;
};

const hype_field<hype_general_parameters> k_general_fields[] = {
    {"ttpd", &hype_general_parameters::ttpd},
    {"ttpi", &hype_general_parameters::ttpi},
    {"tcalt", &hype_general_parameters::tcalt},
    {"lp", &hype_general_parameters::lp},
    {"epotdist", &hype_general_parameters::epotdist},
    {"fsceff", &hype_general_parameters::fsceff},
    {"whcsnow", &hype_general_parameters::whcsnow},
    {"rrcs3", &hype_general_parameters::rrcs3},
    {"rivvel", &hype_general_parameters::rivvel},
    {"damp", &hype_general_parameters::damp},
    {"cevpam", &hype_general_parameters::cevpam},
    {"cevpph", &hype_general_parameters::cevpph},
};

const hype_field<hype_landuse_parameters> k_landuse_fields[] = {
    {"ttmp", &hype_landuse_parameters::ttmp},
    {"cmlt", &hype_landuse_parameters::cmlt},
    {"cevp", &hype_landuse_parameters::cevp},
    {"srrcs", &hype_landuse_parameters::srrcs},
    {"pcluse", &hype_landuse_parameters::pcluse},
};

const hype_field<hype_soil_parameters> k_soil_fields[] = {
    {"wcwp", &hype_soil_parameters::wcwp},
    {"wcfc", &hype_soil_parameters::wcfc},
    {"wcep", &hype_soil_parameters::wcep},
    {"rrcs1", &hype_soil_parameters::rrcs1},
    {"rrcs2", &hype_soil_parameters::rrcs2},
    {"mperc1", &hype_soil_parameters::mperc1},
    {"mperc2", &hype_soil_parameters::mperc2},
    {"macrate", &hype_soil_parameters::macrate},
    {"srrate", &hype_soil_parameters::srrate},
    {"mactrinf", &hype_soil_parameters::mactrinf},
    {"mactrsm", &hype_soil_parameters::mactrsm},
};

template <typename T, std::size_t N>
std::vector<std::string> field_names(const hype_field<T> (&fields)[N])
{
    std::vector<std::string> names;
    names.reserve(N);
    for (const auto& f : fields) names.emplace_back(f.name);
    return names;
}

const char* const k_layer_columns[3] = {"layer1_bottom_m", "layer2_bottom_m", "layer3_bottom_m"};

// ---------------------------------------------------------------------------
// Message helpers
// ---------------------------------------------------------------------------
std::string num(double v)
{
    std::ostringstream os;
    os.precision(10);
    os << v;
    return os.str();
}

// "<label> row <n>, column "<col>": "
std::string at(const std::string& label, std::size_t row_1based, const std::string& column)
{
    return label + " row " + std::to_string(row_1based) + ", column \"" + column + "\": ";
}

std::string join(const std::vector<std::string>& v)
{
    std::string s;
    for (std::size_t i = 0; i < v.size(); ++i)
    {
        if (i) s += ", ";
        s += v[i];
    }
    return s;
}

std::string table_path(const std::string& dir, const char* file)
{
    return (std::filesystem::path(dir) / file).string();
}

// How validate_hype_setup names a table: the file if the set-up came from a
// directory, otherwise just the table.
std::string table_label(const hype_setup& s, const char* file)
{
    if (!s.source_dir.empty()) return table_path(s.source_dir, file);
    return std::string(file) + " table";
}

// ---------------------------------------------------------------------------
// Column and cell helpers
// ---------------------------------------------------------------------------

// Rejects duplicate, unknown and missing columns. An unknown column is the
// check that matters: with "absent means zero", a misspelt parameter column
// would otherwise silently read as zero.
void check_columns(const lsdtt_csv_table& t,
                   const std::string& path,
                   const std::vector<std::string>& required,
                   const std::vector<std::string>& optional)
{
    std::unordered_set<std::string> allowed(required.begin(), required.end());
    allowed.insert(optional.begin(), optional.end());

    std::unordered_set<std::string> seen;
    for (const std::string& h : t.headers)
    {
        if (!seen.insert(h).second)
            throw std::invalid_argument(path + " header: column \"" + h + "\" appears more than once");
        if (allowed.count(h) == 0)
        {
            std::vector<std::string> all = required;
            all.insert(all.end(), optional.begin(), optional.end());
            throw std::invalid_argument(path + " header: unknown column \"" + h +
                                        "\" (accepted: " + join(all) + ")");
        }
    }
    for (const std::string& r : required)
        if (seen.count(r) == 0)
            throw std::invalid_argument(path + " header: missing required column \"" + r + "\"");
}

bool cell_blank(const lsdtt_csv_row& row, const std::string& key)
{
    return csv_row_string(row, key).empty();
}

int row_int(const lsdtt_csv_row& row, const std::string& key, const std::string& path, std::size_t r)
{
    const std::string s = csv_row_string(row, key);
    if (s.empty())
        throw std::invalid_argument(at(path, r, key) + "missing required value");
    std::size_t pos = 0;
    long long v = 0;
    try
    {
        v = std::stoll(s, &pos);
    }
    catch (...)
    {
        pos = 0;
    }
    if (pos == 0 || pos != s.size() || v < INT_MIN || v > INT_MAX)
        throw std::invalid_argument(at(path, r, key) + "\"" + s + "\" is not an integer");
    return static_cast<int>(v);
}

// A numeric cell. A required column must have a value in every row. An
// optional column that is ABSENT reads as `fallback`; an optional column that
// is PRESENT must still have a value in every row, because a blank cell in a
// column the user wrote is more likely an omission than a deliberate zero.
double row_double(const lsdtt_csv_table& t, const lsdtt_csv_row& row, const std::string& key,
                  double fallback, const std::string& path, std::size_t r)
{
    if (!csv_has_column(t, key)) return fallback;
    const double v = csv_row_double_required(row, key, path, static_cast<int>(r));
    if (!std::isfinite(v))
        throw std::invalid_argument(at(path, r, key) + "value is not finite");
    return v;
}

template <typename T, std::size_t N>
void read_fields(const hype_field<T> (&fields)[N], const lsdtt_csv_table& t,
                 const lsdtt_csv_row& row, const std::string& path, std::size_t r, T& out)
{
    for (const auto& f : fields) out.*(f.member) = row_double(t, row, f.name, 0.0, path, r);
}

// ---------------------------------------------------------------------------
// Single-table validators. Row numbers are index + 1: rows keep file order.
// ---------------------------------------------------------------------------
template <typename T>
void check_unique_codes(const std::vector<T>& rows, int T::*key, const char* column,
                        const std::string& label)
{
    std::unordered_map<int, std::size_t> first;
    for (std::size_t i = 0; i < rows.size(); ++i)
    {
        const int k = rows[i].*key;
        auto ins = first.emplace(k, i);
        if (!ins.second)
            throw std::invalid_argument(at(label, i + 1, column) + "duplicate " + column + " " +
                                        std::to_string(k) + " (first at row " +
                                        std::to_string(ins.first->second + 1) + ")");
    }
}

void validate_classes_table(const std::vector<hype_class_definition>& classes,
                            const std::string& label)
{
    check_unique_codes(classes, &hype_class_definition::class_id, "class_id", label);
    for (std::size_t i = 0; i < classes.size(); ++i)
    {
        const hype_class_definition& c = classes[i];
        const std::size_t r = i + 1;
        const std::string cid = " (class_id " + std::to_string(c.class_id) + ")";
        if (c.n_layers < 1 || c.n_layers > 3)
            throw std::invalid_argument(at(label, r, "n_layers") + "n_layers is " +
                                        std::to_string(c.n_layers) + cid + "; must be 1, 2 or 3");
        double previous = 0.0;
        for (int k = 0; k < 3; ++k)
        {
            const double d = c.layer_bottom_m[static_cast<std::size_t>(k)];
            if (k < c.n_layers)
            {
                if (!(d > previous))
                    throw std::invalid_argument(
                        at(label, r, k_layer_columns[k]) + "layer bottom " + num(d) + " m" + cid +
                        (k == 0 ? " must be positive"
                                : " must be deeper than the layer above (" + num(previous) + " m)"));
                previous = d;
            }
            else if (d != 0.0)
            {
                throw std::invalid_argument(at(label, r, k_layer_columns[k]) + "layer bottom set" +
                                            cid + " but n_layers is " + std::to_string(c.n_layers));
            }
        }
        if (!(c.stream_depth_m >= 0.0))
            throw std::invalid_argument(at(label, r, "stream_depth_m") + "stream depth " +
                                        num(c.stream_depth_m) + " m" + cid + " is negative");
        if (!(c.tile_depth_m >= 0.0))
            throw std::invalid_argument(at(label, r, "tile_depth_m") + "tile depth " +
                                        num(c.tile_depth_m) + " m" + cid + " is negative");
    }
}

// Kahn's algorithm over the maindown links. Each sub-basin has at most one
// downstream neighbour, so whatever Kahn's leaves unvisited is exactly the set
// of sub-basins on cycles (a cycle has no exit), and following maindown from
// any of them walks the cycle.
std::vector<std::size_t> downstream_order(const std::vector<hype_subbasin>& subbasins,
                                          const std::string& label)
{
    const std::size_t n = subbasins.size();
    std::unordered_map<int, std::size_t> index;
    for (std::size_t i = 0; i < n; ++i)
    {
        auto ins = index.emplace(subbasins[i].subid, i);
        if (!ins.second)
            throw std::invalid_argument(at(label, i + 1, "subid") + "duplicate subid " +
                                        std::to_string(subbasins[i].subid) + " (first at row " +
                                        std::to_string(ins.first->second + 1) + ")");
    }

    // down[i] = index of the downstream sub-basin, or n for the domain outlet.
    std::vector<std::size_t> down(n, n);
    std::vector<std::size_t> indegree(n, 0);
    for (std::size_t i = 0; i < n; ++i)
    {
        const int md = subbasins[i].maindown;
        if (md == 0) continue;
        auto it = index.find(md);
        if (it == index.end())
            throw std::invalid_argument(at(label, i + 1, "maindown") + "maindown " +
                                        std::to_string(md) + " of subid " +
                                        std::to_string(subbasins[i].subid) +
                                        " is not a subid in the table (0 marks the domain outlet)");
        down[i] = it->second;
        ++indegree[it->second];
    }

    std::vector<std::size_t> order;
    order.reserve(n);
    for (std::size_t i = 0; i < n; ++i)
        if (indegree[i] == 0) order.push_back(i);
    // `order` doubles as the FIFO queue: ties keep input order.
    for (std::size_t head = 0; head < order.size(); ++head)
    {
        const std::size_t d = down[order[head]];
        if (d != n && --indegree[d] == 0) order.push_back(d);
    }

    if (order.size() != n)
    {
        std::size_t start = 0;
        while (indegree[start] == 0) ++start;
        std::string path = std::to_string(subbasins[start].subid);
        for (std::size_t j = down[start]; ; j = down[j])
        {
            path += " -> " + std::to_string(subbasins[j].subid);
            if (j == start) break;
        }
        throw std::invalid_argument(at(label, start + 1, "maindown") +
                                    "maindown links form a cycle: " + path);
    }
    return order;
}

void validate_subbasins_table(const std::vector<hype_subbasin>& subbasins, const std::string& label)
{
    for (std::size_t i = 0; i < subbasins.size(); ++i)
    {
        const hype_subbasin& s = subbasins[i];
        if (s.subid <= 0)
            throw std::invalid_argument(at(label, i + 1, "subid") + "subid " +
                                        std::to_string(s.subid) + " must be positive");
        if (!(s.area_m2 > 0.0))
            throw std::invalid_argument(at(label, i + 1, "area_m2") + "area " + num(s.area_m2) +
                                        " m^2 of subid " + std::to_string(s.subid) +
                                        " must be positive");
    }
    // Duplicates, unknown maindown and cycles.
    (void)downstream_order(subbasins, label);
}

void validate_subbasin_classes_table(const std::vector<hype_subbasin_class>& rows,
                                     const std::string& label)
{
    std::unordered_map<long long, std::size_t> first;
    for (std::size_t i = 0; i < rows.size(); ++i)
    {
        const hype_subbasin_class& sc = rows[i];
        const std::string pair = " (subid " + std::to_string(sc.subid) + ", class_id " +
                                 std::to_string(sc.class_id) + ")";
        if (!(sc.fraction >= 0.0 && sc.fraction <= 1.0))
            throw std::invalid_argument(at(label, i + 1, "fraction") + "fraction " +
                                        num(sc.fraction) + pair + " is outside [0, 1]");
        const long long key = (static_cast<long long>(sc.subid) << 32) ^
                              static_cast<long long>(static_cast<unsigned int>(sc.class_id));
        auto ins = first.emplace(key, i);
        if (!ins.second)
            throw std::invalid_argument(at(label, i + 1, "class_id") + "duplicate row" + pair +
                                        " (first at row " +
                                        std::to_string(ins.first->second + 1) + ")");
    }
}

} // namespace

// ---------------------------------------------------------------------------
// Names
// ---------------------------------------------------------------------------
const std::vector<std::string>& hype_general_parameter_names()
{
    static const std::vector<std::string> names = field_names(k_general_fields);
    return names;
}

const std::vector<std::string>& hype_landuse_parameter_names()
{
    static const std::vector<std::string> names = field_names(k_landuse_fields);
    return names;
}

const std::vector<std::string>& hype_soil_parameter_names()
{
    static const std::vector<std::string> names = field_names(k_soil_fields);
    return names;
}

// ---------------------------------------------------------------------------
// Readers
// ---------------------------------------------------------------------------
hype_general_parameters load_hype_general_parameters(const std::string& path)
{
    const lsdtt_csv_table t = load_csv_table(path);
    check_columns(t, path, {"name", "value"}, {});

    hype_general_parameters p;
    std::unordered_map<std::string, std::size_t> seen;
    for (std::size_t i = 0; i < t.rows.size(); ++i)
    {
        const lsdtt_csv_row& row = t.rows[i];
        const std::size_t r = i + 1;
        const std::string name = csv_row_string(row, "name");
        if (name.empty())
            throw std::invalid_argument(at(path, r, "name") + "missing parameter name");

        const hype_field<hype_general_parameters>* field = nullptr;
        for (const auto& f : k_general_fields)
            if (name == f.name) field = &f;
        if (field == nullptr)
            throw std::invalid_argument(at(path, r, "name") + "unknown parameter \"" + name +
                                        "\" (accepted: " + join(hype_general_parameter_names()) +
                                        ")");

        auto ins = seen.emplace(name, r);
        if (!ins.second)
            throw std::invalid_argument(at(path, r, "name") + "parameter \"" + name +
                                        "\" given twice (first at row " +
                                        std::to_string(ins.first->second) + ")");

        p.*(field->member) = row_double(t, row, "value", 0.0, path, r);
    }
    return p;
}

std::vector<hype_landuse_parameters> load_hype_landuse_parameters(const std::string& path)
{
    const lsdtt_csv_table t = load_csv_table(path);
    check_columns(t, path, {"code"}, hype_landuse_parameter_names());

    std::vector<hype_landuse_parameters> out;
    out.reserve(t.rows.size());
    for (std::size_t i = 0; i < t.rows.size(); ++i)
    {
        hype_landuse_parameters lu;
        lu.code = row_int(t.rows[i], "code", path, i + 1);
        read_fields(k_landuse_fields, t, t.rows[i], path, i + 1, lu);
        out.push_back(lu);
    }
    check_unique_codes(out, &hype_landuse_parameters::code, "code", path);
    return out;
}

std::vector<hype_soil_parameters> load_hype_soil_parameters(const std::string& path)
{
    const lsdtt_csv_table t = load_csv_table(path);
    check_columns(t, path, {"code"}, hype_soil_parameter_names());

    std::vector<hype_soil_parameters> out;
    out.reserve(t.rows.size());
    for (std::size_t i = 0; i < t.rows.size(); ++i)
    {
        hype_soil_parameters so;
        so.code = row_int(t.rows[i], "code", path, i + 1);
        read_fields(k_soil_fields, t, t.rows[i], path, i + 1, so);
        out.push_back(so);
    }
    check_unique_codes(out, &hype_soil_parameters::code, "code", path);
    return out;
}

std::vector<hype_class_definition> load_hype_classes(const std::string& path)
{
    const lsdtt_csv_table t = load_csv_table(path);
    check_columns(t, path, {"class_id", "landuse", "soil", "n_layers", k_layer_columns[0]},
                  {k_layer_columns[1], k_layer_columns[2], "stream_depth_m", "tile_depth_m"});

    std::vector<hype_class_definition> out;
    out.reserve(t.rows.size());
    for (std::size_t i = 0; i < t.rows.size(); ++i)
    {
        const lsdtt_csv_row& row = t.rows[i];
        const std::size_t r = i + 1;
        hype_class_definition c;
        c.class_id = row_int(row, "class_id", path, r);
        c.landuse = row_int(row, "landuse", path, r);
        c.soil = row_int(row, "soil", path, r);
        c.n_layers = row_int(row, "n_layers", path, r);
        if (c.n_layers < 1 || c.n_layers > 3)
            throw std::invalid_argument(at(path, r, "n_layers") + "n_layers is " +
                                        std::to_string(c.n_layers) + " (class_id " +
                                        std::to_string(c.class_id) + "); must be 1, 2 or 3");

        // Layers beyond n_layers are blank cells ("0.3,,") or absent columns;
        // the number of filled depth cells must equal n_layers.
        for (int k = 0; k < 3; ++k)
        {
            const std::string col = k_layer_columns[k];
            const bool blank = cell_blank(row, col);
            if (k < c.n_layers && blank)
                throw std::invalid_argument(at(path, r, col) + "blank, but n_layers is " +
                                            std::to_string(c.n_layers) + " (class_id " +
                                            std::to_string(c.class_id) + ")");
            if (k >= c.n_layers && !blank)
                throw std::invalid_argument(at(path, r, col) + "set, but n_layers is " +
                                            std::to_string(c.n_layers) + " (class_id " +
                                            std::to_string(c.class_id) +
                                            "); leave unused layers blank");
            if (!blank)
                c.layer_bottom_m[static_cast<std::size_t>(k)] =
                    row_double(t, row, col, 0.0, path, r);
        }
        c.stream_depth_m = row_double(t, row, "stream_depth_m", 0.0, path, r);
        c.tile_depth_m = row_double(t, row, "tile_depth_m", 0.0, path, r);
        out.push_back(c);
    }
    validate_classes_table(out, path);
    return out;
}

std::vector<hype_subbasin> load_hype_subbasins(const std::string& path)
{
    const lsdtt_csv_table t = load_csv_table(path);
    check_columns(t, path,
                  {"subid", "maindown", "area_m2", "elev_mean_m", "slope_mean_percent",
                   "loc_rivlen_m", "rivlen_m"},
                  {"elev_std_m"});

    std::vector<hype_subbasin> out;
    out.reserve(t.rows.size());
    for (std::size_t i = 0; i < t.rows.size(); ++i)
    {
        const lsdtt_csv_row& row = t.rows[i];
        const std::size_t r = i + 1;
        hype_subbasin s;
        s.subid = row_int(row, "subid", path, r);
        s.maindown = row_int(row, "maindown", path, r);
        s.area_m2 = row_double(t, row, "area_m2", 0.0, path, r);
        s.elev_mean_m = row_double(t, row, "elev_mean_m", 0.0, path, r);
        s.elev_std_m = row_double(t, row, "elev_std_m", 0.0, path, r);
        s.slope_mean_percent = row_double(t, row, "slope_mean_percent", 0.0, path, r);
        s.loc_rivlen_m = row_double(t, row, "loc_rivlen_m", 0.0, path, r);
        s.rivlen_m = row_double(t, row, "rivlen_m", 0.0, path, r);
        out.push_back(s);
    }
    validate_subbasins_table(out, path);
    return out;
}

std::vector<hype_subbasin_class> load_hype_subbasin_classes(const std::string& path)
{
    const lsdtt_csv_table t = load_csv_table(path);
    check_columns(t, path, {"subid", "class_id", "fraction"}, {"dh_m"});

    std::vector<hype_subbasin_class> out;
    out.reserve(t.rows.size());
    for (std::size_t i = 0; i < t.rows.size(); ++i)
    {
        const lsdtt_csv_row& row = t.rows[i];
        const std::size_t r = i + 1;
        hype_subbasin_class sc;
        sc.subid = row_int(row, "subid", path, r);
        sc.class_id = row_int(row, "class_id", path, r);
        sc.fraction = row_double(t, row, "fraction", 0.0, path, r);
        sc.dh_m = row_double(t, row, "dh_m", 0.0, path, r);
        out.push_back(sc);
    }
    validate_subbasin_classes_table(out, path);
    return out;
}

hype_setup load_hype_setup_directory(const std::string& dir)
{
    if (!std::filesystem::is_directory(dir))
        throw std::runtime_error("HYPE set-up directory does not exist: " + dir);

    hype_setup s;
    s.source_dir = dir;
    s.general = load_hype_general_parameters(table_path(dir, "general_parameters.csv"));
    s.landuse = load_hype_landuse_parameters(table_path(dir, "landuse_parameters.csv"));
    s.soils = load_hype_soil_parameters(table_path(dir, "soil_parameters.csv"));
    s.classes = load_hype_classes(table_path(dir, "classes.csv"));
    s.subbasins = load_hype_subbasins(table_path(dir, "subbasins.csv"));
    s.subbasin_classes = load_hype_subbasin_classes(table_path(dir, "subbasin_classes.csv"));
    validate_hype_setup(s);
    return s;
}

// ---------------------------------------------------------------------------
// Validation
// ---------------------------------------------------------------------------
void validate_hype_setup(const hype_setup& s)
{
    const std::string landuse_label = table_label(s, "landuse_parameters.csv");
    const std::string soil_label = table_label(s, "soil_parameters.csv");
    const std::string class_label = table_label(s, "classes.csv");
    const std::string subbasin_label = table_label(s, "subbasins.csv");
    const std::string sc_label = table_label(s, "subbasin_classes.csv");

    // Single-table checks again: a set-up built in memory has not been
    // through the readers.
    check_unique_codes(s.landuse, &hype_landuse_parameters::code, "code", landuse_label);
    check_unique_codes(s.soils, &hype_soil_parameters::code, "code", soil_label);
    validate_classes_table(s.classes, class_label);
    validate_subbasins_table(s.subbasins, subbasin_label);
    validate_subbasin_classes_table(s.subbasin_classes, sc_label);

    // Parameter ranges. HYPE itself clamps or silently misbehaves on some of
    // these (e.g. a negative srrcs makes the saturation-excess term remove
    // more water than layer 1 holds, and the column stops); reject them here.
    const std::string general_label = table_label(s, "general_parameters.csv");
    auto require = [](bool ok, const std::string& where, const std::string& what) {
        if (!ok) throw std::invalid_argument(where + what);
    };
    auto general_at = [&](const char* name) { return general_label + ", parameter \"" + std::string(name) + "\": "; };
    const hype_general_parameters& g = s.general;
    require(g.ttpi >= 0.0, general_at("ttpi"), "must be >= 0, got " + num(g.ttpi));
    require(g.lp >= 0.0, general_at("lp"), "must be >= 0, got " + num(g.lp));
    require(g.epotdist >= 0.0, general_at("epotdist"), "must be >= 0, got " + num(g.epotdist));
    require(g.fsceff >= 0.0 && g.fsceff <= 1.0, general_at("fsceff"), "must be in [0, 1], got " + num(g.fsceff));
    require(g.whcsnow >= 0.0 && g.whcsnow < 1.0, general_at("whcsnow"), "must be in [0, 1), got " + num(g.whcsnow));
    require(g.rrcs3 >= 0.0, general_at("rrcs3"), "must be >= 0, got " + num(g.rrcs3));
    require(g.rivvel >= 0.0, general_at("rivvel"), "must be >= 0, got " + num(g.rivvel));
    require(g.damp >= 0.0 && g.damp <= 1.0, general_at("damp"), "must be in [0, 1], got " + num(g.damp));
    for (std::size_t i = 0; i < s.landuse.size(); ++i)
    {
        const hype_landuse_parameters& l = s.landuse[i];
        require(l.cmlt >= 0.0, at(landuse_label, i + 1, "cmlt"), "must be >= 0, got " + num(l.cmlt));
        require(l.cevp >= 0.0, at(landuse_label, i + 1, "cevp"), "must be >= 0, got " + num(l.cevp));
        require(l.srrcs >= 0.0, at(landuse_label, i + 1, "srrcs"), "must be >= 0, got " + num(l.srrcs));
        require(l.pcluse >= 0.0 && l.pcluse <= 1.0, at(landuse_label, i + 1, "pcluse"),
                "must be in [0, 1], got " + num(l.pcluse));
    }
    for (std::size_t i = 0; i < s.soils.size(); ++i)
    {
        const hype_soil_parameters& o = s.soils[i];
        const std::pair<const char*, double> nonneg[] = {
            {"wcwp", o.wcwp}, {"wcfc", o.wcfc}, {"wcep", o.wcep}, {"rrcs1", o.rrcs1}, {"rrcs2", o.rrcs2},
            {"mperc1", o.mperc1}, {"mperc2", o.mperc2}, {"macrate", o.macrate}, {"srrate", o.srrate},
            {"mactrinf", o.mactrinf}, {"mactrsm", o.mactrsm}};
        for (const auto& [name, v] : nonneg)
            require(v >= 0.0, at(soil_label, i + 1, name), "must be >= 0, got " + num(v));
        const double porosity = o.wcwp + o.wcfc + o.wcep;
        require(porosity <= 1.0, at(soil_label, i + 1, "wcwp+wcfc+wcep"),
                "the stacked volume fractions sum to " + num(porosity) + "; they cannot exceed 1");
    }

    // Sub-basin and class values must be finite. A sub-basin derived from
    // rasters can carry NaN (e.g. a slope over pixels whose gradient is all
    // nodata), and the column fails on it much later with an unrelated message.
    for (std::size_t i = 0; i < s.subbasins.size(); ++i)
    {
        const hype_subbasin& b = s.subbasins[i];
        const std::pair<const char*, double> vals[] = {
            {"area_m2", b.area_m2}, {"elev_mean_m", b.elev_mean_m}, {"elev_std_m", b.elev_std_m},
            {"slope_mean_percent", b.slope_mean_percent}, {"loc_rivlen_m", b.loc_rivlen_m},
            {"rivlen_m", b.rivlen_m}};
        for (const auto& [name, v] : vals)
            require(std::isfinite(v), at(subbasin_label, i + 1, name),
                    "subid " + std::to_string(b.subid) + " has a non-finite value (" + num(v) + ")");
    }
    for (std::size_t i = 0; i < s.subbasin_classes.size(); ++i)
    {
        const hype_subbasin_class& sc = s.subbasin_classes[i];
        require(std::isfinite(sc.fraction), at(sc_label, i + 1, "fraction"), "is not finite (" + num(sc.fraction) + ")");
        require(std::isfinite(sc.dh_m), at(sc_label, i + 1, "dh_m"), "is not finite (" + num(sc.dh_m) + ")");
    }

    // Classes -> land use and soil.
    std::unordered_set<int> landuse_codes, soil_codes, class_ids;
    for (const auto& lu : s.landuse) landuse_codes.insert(lu.code);
    for (const auto& so : s.soils) soil_codes.insert(so.code);
    for (std::size_t i = 0; i < s.classes.size(); ++i)
    {
        const hype_class_definition& c = s.classes[i];
        class_ids.insert(c.class_id);
        if (landuse_codes.count(c.landuse) == 0)
            throw std::invalid_argument(at(class_label, i + 1, "landuse") + "land-use code " +
                                        std::to_string(c.landuse) + " of class_id " +
                                        std::to_string(c.class_id) + " is not in " +
                                        landuse_label);
        if (soil_codes.count(c.soil) == 0)
            throw std::invalid_argument(at(class_label, i + 1, "soil") + "soil code " +
                                        std::to_string(c.soil) + " of class_id " +
                                        std::to_string(c.class_id) + " is not in " + soil_label);
    }

    // (Sub-basin, class) rows -> sub-basins and classes; fractions sum to 1.
    std::unordered_map<int, double> fraction_sum;
    for (const auto& sb : s.subbasins) fraction_sum.emplace(sb.subid, 0.0);
    for (std::size_t i = 0; i < s.subbasin_classes.size(); ++i)
    {
        const hype_subbasin_class& sc = s.subbasin_classes[i];
        auto it = fraction_sum.find(sc.subid);
        if (it == fraction_sum.end())
            throw std::invalid_argument(at(sc_label, i + 1, "subid") + "subid " +
                                        std::to_string(sc.subid) + " is not in " + subbasin_label);
        if (class_ids.count(sc.class_id) == 0)
            throw std::invalid_argument(at(sc_label, i + 1, "class_id") + "class_id " +
                                        std::to_string(sc.class_id) + " is not in " + class_label);
        it->second += sc.fraction;
    }
    constexpr double k_fraction_tolerance = 1e-6;
    for (const auto& sb : s.subbasins)
    {
        const double sum = fraction_sum[sb.subid];
        if (std::fabs(sum - 1.0) > k_fraction_tolerance)
            throw std::invalid_argument(sc_label + ", column \"fraction\": fractions of subid " +
                                        std::to_string(sb.subid) + " sum to " + num(sum) +
                                        "; they must sum to 1 within 1e-6");
    }
}

std::vector<std::size_t> hype_subbasin_downstream_order(const std::vector<hype_subbasin>& subbasins)
{
    return downstream_order(subbasins, "subbasins table");
}

} // namespace lsdtt
