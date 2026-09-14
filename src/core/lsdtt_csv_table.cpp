// lsdtt_csv_table.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html

#include "lsdtt_csv_table.hpp"

#include <fstream>
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

// A quote-aware splitter, RFC 4180 shaped but deliberately lenient. The exact
// contract, including the malformed cases, is documented in the header; keep
// the two in step.
static bool is_csv_ws(const char c)
{
    return c == ' ' || c == '\t' || c == '\n' || c == '\v' || c == '\f' || c == '\r';
}

std::vector<std::string> csv_split_row(const std::string& row)
{
    std::vector<std::string> out;
    // An empty line has no fields at all, not one empty field. load_csv_table
    // uses this to reject a header-less file.
    if (row.empty()) return out;

    const std::size_t n = row.size();
    std::size_t i = 0;
    while (true)
    {
        std::string field;

        // Whitespace before an opening quote is dropped with the quote; if no
        // quote follows, the field is verbatim from `i` and that whitespace is
        // part of it (csv_trim, in the caller, deals with it).
        std::size_t j = i;
        while (j < n && is_csv_ws(row[j])) ++j;

        if (j < n && row[j] == '"')
        {
            std::size_t k = j + 1;
            while (k < n)
            {
                if (row[k] == '"')
                {
                    // "" inside quotes is one literal quote; a lone " closes.
                    if (k + 1 < n && row[k + 1] == '"') { field.push_back('"'); k += 2; continue; }
                    ++k;
                    break;
                }
                field.push_back(row[k]);
                ++k;
            }
            // Lenient: an unterminated quote was just closed by end of line,
            // and any text between a closing quote and the next separating
            // comma is appended verbatim rather than being an error.
            while (k < n && row[k] != ',') { field.push_back(row[k]); ++k; }
            i = k;
        }
        else
        {
            std::size_t k = i;
            while (k < n && row[k] != ',') { field.push_back(row[k]); ++k; }
            i = k;
        }

        out.push_back(field);
        if (i >= n) break;
        ++i; // step over the separating comma; a trailing one yields an empty field
    }
    return out;
}

std::vector<std::string> csv_split_row_trimmed(const std::string& row)
{
    std::vector<std::string> fields = csv_split_row(row);
    for (std::string& f : fields) f = csv_trim(f);
    return fields;
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

    for (const std::string& tok : csv_split_row_trimmed(line))
        t.headers.push_back(tok);
    if (t.headers.empty())
        throw std::runtime_error("Input CSV has no header columns: " + path);

    while (std::getline(f, line))
    {
        if (csv_trim(line).empty()) continue;
        lsdtt_csv_row row;
        // Both the header and the body go through csv_split_row so a quoted
        // field means the same thing in either -- these used to be two inline
        // getline loops that knew nothing about quotes.
        const std::vector<std::string> toks = csv_split_row(line);
        std::size_t i = 0;
        // Over-long rows are truncated to the header count.
        for (; i < toks.size() && i < t.headers.size(); ++i)
            row[t.headers[i]] = csv_trim(toks[i]);
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
