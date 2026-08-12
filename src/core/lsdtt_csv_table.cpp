// lsdtt_csv_table.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html

#include "lsdtt_csv_table.hpp"

#include <fstream>
#include <sstream>
#include <stdexcept>

namespace lsdtt
{

// The full isspace set, not just " \t\r\n". The two CLI copies this replaces
// disagreed here: one trimmed " \t\r\n", the other used std::isspace and so
// also took \v and \f. Trimming the superset is the safe reconciliation -- a
// field containing \v or \f would previously have failed std::stod under the
// narrower rule, so nothing that used to parse can start failing.
static const char* k_csv_ws = " \t\n\v\f\r";

std::string csv_trim(const std::string& s)
{
    const std::size_t b = s.find_first_not_of(k_csv_ws);
    if (b == std::string::npos) return "";
    const std::size_t e = s.find_last_not_of(k_csv_ws);
    return s.substr(b, e - b + 1);
}

std::vector<std::string> csv_split_row(const std::string& row)
{
    std::vector<std::string> out;
    std::istringstream ss(row);
    std::string item;
    while (std::getline(ss, item, ',')) out.push_back(item);
    return out;
}

lsdtt_csv_table load_csv_table(const std::string& path)
{
    std::ifstream f(path);
    if (!f.is_open())
        throw std::runtime_error("Cannot open input CSV: " + path);

    lsdtt_csv_table t;
    std::string line;
    if (!std::getline(f, line))
        throw std::runtime_error("Input CSV is empty: " + path);

    std::istringstream hs(line);
    std::string tok;
    while (std::getline(hs, tok, ','))
        t.headers.push_back(csv_trim(tok));
    if (t.headers.empty())
        throw std::runtime_error("Input CSV has no header columns: " + path);

    while (std::getline(f, line))
    {
        if (csv_trim(line).empty()) continue;
        lsdtt_csv_row row;
        std::istringstream rs(line);
        std::size_t i = 0;
        while (std::getline(rs, tok, ',') && i < t.headers.size())
        {
            row[t.headers[i]] = csv_trim(tok);
            ++i;
        }
        // Short rows are padded so every header key exists; see the header note.
        while (i < t.headers.size())
        {
            row[t.headers[i]] = "";
            ++i;
        }
        t.rows.push_back(std::move(row));
    }
    return t;
}

bool csv_has_column(const lsdtt_csv_table& t, const std::string& key)
{
    for (const std::string& h : t.headers)
        if (h == key) return true;
    return false;
}

std::string csv_row_string(const lsdtt_csv_row& row,
                           const std::string& key,
                           const std::string& fallback)
{
    auto it = row.find(key);
    if (it == row.end()) return fallback;
    return it->second.empty() ? fallback : it->second;
}

double csv_row_double_required(const lsdtt_csv_row& row,
                              const std::string& key,
                              const std::string& csv_path,
                              const int row_index_1based)
{
    auto it = row.find(key);
    if (it == row.end() || csv_trim(it->second).empty())
        throw std::runtime_error("Missing required field \"" + key + "\" in " +
                                 csv_path + " row " + std::to_string(row_index_1based));
    try
    {
        return std::stod(it->second);
    }
    catch (...)
    {
        throw std::runtime_error("Field \"" + key + "\" is not numeric in " +
                                 csv_path + " row " + std::to_string(row_index_1based) +
                                 " (value: \"" + it->second + "\")");
    }
}

double csv_row_double_optional(const lsdtt_csv_row& row,
                              const std::string& key,
                              const double fallback,
                              const std::string& csv_path,
                              const int row_index_1based)
{
    auto it = row.find(key);
    if (it == row.end() || csv_trim(it->second).empty()) return fallback;
    try
    {
        return std::stod(it->second);
    }
    catch (...)
    {
        throw std::runtime_error("Field \"" + key + "\" is not numeric in " +
                                 csv_path + " row " + std::to_string(row_index_1based) +
                                 " (value: \"" + it->second + "\")");
    }
}

} // namespace lsdtt
