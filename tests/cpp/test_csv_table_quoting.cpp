// test_csv_table_quoting.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// The quote-aware CSV splitter and table reader (src/core/lsdtt_csv_table,
// decision D5), against the contract written in lsdtt_csv_table.hpp.
//
// Several CLIs read user tables through load_csv_table (the HYPE set-up, the
// distributed-hydrology forcing, sample and parameter tables), and until D5
// each had its own getline splitter that knew nothing about quotes. Nothing
// else tests the splitter directly, so every case the header documents is
// pinned here, with the expected fields written out by hand from the header:
//
//   S  csv_split_row, the documented cases one by one, including the lenient
//      readings of malformed input (unterminated quote, a quote mid-field,
//      text after a closing quote) and the D5 change: a trailing comma gives
//      a trailing empty field ("a," -> {a, ""}, not {a}).
//   T  csv_trim and csv_split_row_trimmed.
//   L  load_csv_table on files written here: a quoted header name, a comma
//      inside a quoted value, CRLF line ends, a blank line skipped, a short
//      row padded with empty strings, an over-long row truncated; and an
//      empty file and a missing file refused.
//
// Usage: lsdtt-test-csv-table-quoting [scratch_dir]
// Prints one line per check and a final "RESULT: PASS|FAIL".
// -----------------------------------------------------------------------------

#include "lsdtt_csv_table.hpp"

#include <cstdio>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace
{
int g_n_checks = 0;
int g_n_failed = 0;

void check(const std::string& name, bool ok)
{
    ++g_n_checks;
    if (!ok) ++g_n_failed;
    std::cout << (ok ? "PASS" : "FAIL") << "  " << name << "\n";
}

std::string show(const std::vector<std::string>& v)
{
    std::string s = "{";
    for (std::size_t i = 0; i < v.size(); ++i)
        s += (i ? ", " : "") + std::string("[") + v[i] + "]";
    return s + "}";
}

void split_is(const std::string& name, const std::string& row,
              const std::vector<std::string>& expected)
{
    const std::vector<std::string> got = lsdtt::csv_split_row(row);
    check(name + "  " + show(got), got == expected);
}

void write_file(const std::string& path, const std::string& text)
{
    std::ofstream f(path, std::ios::binary);   // binary: keep the \r\n as written
    f << text;
}

void test_split()
{
    // The documented contract, case by case.
    split_is("S1 unquoted fields verbatim", "a,b,c", {"a", "b", "c"});
    split_is("S2 wrapping quotes removed", "a,\"b\",c", {"a", "b", "c"});
    split_is("S3 comma inside quotes is literal", "a,\"b,c\",d", {"a", "b,c", "d"});
    split_is("S4 doubled quote inside quotes is one quote", "a,\"b\"\"c\",d", {"a", "b\"c", "d"});
    split_is("S5 empty fields", ",,", {"", "", ""});
    split_is("S6 trailing comma gives a trailing empty field (D5)", "a,b,", {"a", "b", ""});
    split_is("S6b single trailing comma (D5: was {a})", "a,", {"a", ""});
    split_is("S7 quoted empty field", "\"\"", {""});
    split_is("S8 empty line has no fields", "", {});
    split_is("S9 whitespace before the opening quote is dropped, text after kept",
             " \"a, b\" ", {"a, b "});
    split_is("S9b unquoted field keeps its whitespace (no trimming here)",
             " a , b ", {" a ", " b "});

    // Malformed input: the lenient readings, none of which throw.
    split_is("S10 unterminated quote runs to the end of the line",
             "a,\"b,c", {"a", "b,c"});
    split_is("S11 a quote mid-field is an ordinary character", "ab\"cd", {"ab\"cd"});
    split_is("S12 text after a closing quote is appended", "\"ab\"cd", {"abcd"});
    split_is("S13 a CRLF line keeps its \\r for the caller to trim",
             "a,b\r", {"a", "b\r"});
    split_is("S13b a CRLF line ending in a quoted field", "a,\"b\"\r", {"a", "b\r"});
}

void test_trim()
{
    check("T1 csv_trim removes spaces, tabs, CR and LF", lsdtt::csv_trim(" \t a b \r\n") == "a b");
    check("T2 csv_trim of whitespace only is empty", lsdtt::csv_trim(" \t\r\n").empty());
    const std::vector<std::string> t = lsdtt::csv_split_row_trimmed(" a , \"b, c\" ,d\r");
    check("T3 csv_split_row_trimmed = split then trim each  " + show(t),
          t == std::vector<std::string>({"a", "b, c", "d"}));
}

void test_load(const std::filesystem::path& dir)
{
    const std::string path = (dir / "quoted_table.csv").string();
    write_file(path,
               "id,\"site, name\",value\r\n"
               "1,\"Glen Doll, upper\",2.5\r\n"
               "\r\n"
               "2,plain,\r\n"
               "3,short\r\n"
               "4,\"long\",7,extra\r\n");
    lsdtt::lsdtt_csv_table t;
    bool loaded = true;
    try { t = lsdtt::load_csv_table(path); } catch (const std::exception&) { loaded = false; }
    check("L1 table loads", loaded);
    if (!loaded) return;

    check("L2 quoted header name with a comma is one column",
          t.headers == std::vector<std::string>({"id", "site, name", "value"}));
    check("L3 blank (CRLF-only) line skipped: 4 rows", t.rows.size() == 4);
    if (t.rows.size() != 4) return;
    check("L4 comma inside a quoted value kept, CR trimmed",
          t.rows[0].at("site, name") == "Glen Doll, upper" && t.rows[0].at("value") == "2.5");
    check("L5 trailing empty field is an empty string",
          t.rows[1].at("site, name") == "plain" && t.rows[1].at("value").empty());
    check("L6 short row padded: every header key present and empty",
          t.rows[2].count("value") == 1 && t.rows[2].at("value").empty());
    check("L7 over-long row truncated to the header count",
          t.rows[3].size() == 3 && t.rows[3].at("value") == "7");
    check("L8 csv_has_column finds the quoted name",
          lsdtt::csv_has_column(t, "site, name") && !lsdtt::csv_has_column(t, "site"));

    const std::string empty = (dir / "empty.csv").string();
    write_file(empty, "");
    bool threw = false;
    try { lsdtt::load_csv_table(empty); } catch (const std::runtime_error&) { threw = true; }
    check("L9 an empty file is refused", threw);

    threw = false;
    try { lsdtt::load_csv_table((dir / "no_such_file.csv").string()); }
    catch (const std::runtime_error&) { threw = true; }
    check("L10 a missing file is refused", threw);
}
}  // namespace

int main(int argc, char** argv)
{
    std::filesystem::path dir = argc > 1
        ? std::filesystem::path(argv[1])
        : std::filesystem::temp_directory_path() / "lsdtt3_test_csv_table_quoting";
    std::filesystem::create_directories(dir);

    test_split();
    test_trim();
    test_load(dir);

    std::cout << "RESULT: " << (g_n_failed == 0 ? "PASS" : "FAIL")
              << "  (" << (g_n_checks - g_n_failed) << "/" << g_n_checks
              << " checks passed)\n";
    return g_n_failed == 0 ? 0 : 1;
}
