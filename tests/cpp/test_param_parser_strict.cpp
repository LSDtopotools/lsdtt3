// test_param_parser_strict.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// unit_param_parser_strict
//
// Where the parser finds the parameter file, and what it prints doing so.
//
// LOCATION. A tool is run as `tool <parameter_file>` or `tool <path>
// <parameter_file>`, and read_path defaults to the directory holding the
// file. A coding agent passes absolute paths by default, and the one-argument
// form used to glue an absolute filename onto "./" (".//abs/dir/run.param"),
// so the file was not found. A filename carrying its own directory must be
// split, absolute or relative, and read_path must be that directory.
//
// SILENCE. Every CLI declares a placeholder parser before assigning the real
// one; the default constructor used to print "created empty parser", which
// appeared on every run and read like an error.
//
// STRICT VALUES. A value the parser cannot honour stops the run, naming the
// key, the value and the fix, instead of being half-read:
//   booleans   only true/True/TRUE/t/T/1 and false/False/FALSE/f/F/0; "yes"
//              used to mean false, silently
//   integers   the whole value; "1e3" and "1.5" used to become 1, "100m" 100
//   reals      the whole value, finite; "100m" used to become 100
//   empty      a numeric or boolean key with no value (used to crash in stoi)
//   quoted     quotes are kept by lsdtt3, so a quoted value is always wrong
//   misspelt   an unknown key that is a near miss of a declared key
// Every problem in the file is reported in the one run. The rules end the
// process (no CLI guards the parse call), so each case runs in a CHILD: the
// driver re-executes itself with --child. The clean file must still parse,
// with the accepted spellings read correctly, so the rules cannot pass by
// rejecting everything.
// -----------------------------------------------------------------------------

#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>

#ifndef _WIN32
#include <sys/wait.h>
#endif

#include "../../src/core/lsdtt_param_parser.hpp"

namespace
{

int g_n_checks = 0;
int g_n_failed = 0;

void check(bool condition, const std::string& what)
{
    ++g_n_checks;
    if (condition)
    {
        std::cout << "PASS  " << what << "\n";
    }
    else
    {
        ++g_n_failed;
        std::cout << "FAIL  " << what << "\n";
    }
}

void write_file(const std::filesystem::path& p, const std::string& text)
{
    std::filesystem::create_directories(p.parent_path());
    std::ofstream out(p);
    out << text;
}

// read_path as the parser resolved it, compared as a canonical directory so
// "a/b/" and "./a/b/" agree.
bool same_dir(const std::string& got, const std::filesystem::path& want)
{
    std::error_code ec;
    const auto a = std::filesystem::weakly_canonical(got, ec);
    const auto b = std::filesystem::weakly_canonical(want, ec);
    return a == b || a / "" == b / "";
}

// Construct a parser and report read_path, or "THREW: <message>".
std::string resolved_read_path(const std::string& path, const std::string& fname)
{
    try
    {
        lsdtt::lsdtt_param_parser pp(path, fname);
        return pp.get_read_path();
    }
    catch (const std::exception& e)
    {
        return std::string("THREW: ") + e.what();
    }
}

// -----------------------------------------------------------------------------
// child: parse with a fixed set of declared keys and print what was read
// -----------------------------------------------------------------------------
int run_child(const std::string& path, const std::string& fname)
{
    lsdtt::lsdtt_param_parser pp(path, fname);
    std::map<std::string, float>       ff = {{"erodibility", 1.0f}};
    std::map<std::string, int>         ii = {{"n_iterations", 10}};
    std::map<std::string, bool>        bb = {{"write_fill_raster", false},
                                             {"write_slope_raster", true}};
    std::map<std::string, std::string> ss = {{"method", "dinf"},
                                             {"K", ""}};
    std::map<std::string, double>      dd = {{"tolerance", 1e-6}};
    pp.parse_all_parameters(ff, ii, bb, ss, dd);

    std::cout << "CHILD erodibility=" << pp.get_float_parameters().at("erodibility")
              << "\nCHILD n_iterations=" << pp.get_int_parameters().at("n_iterations")
              << "\nCHILD write_fill_raster="
              << pp.get_bool_parameters().at("write_fill_raster")
              << "\nCHILD write_slope_raster="
              << pp.get_bool_parameters().at("write_slope_raster")
              << "\nCHILD tolerance=" << pp.get_double_parameters().at("tolerance")
              << "\nCHILD_COMPLETED\n";
    return 0;
}

struct child_run
{
    int         exit_code = 0;
    std::string output;      // stdout and stderr, merged
};

child_run run_case(const std::string& self, const std::filesystem::path& dir,
                   const std::string& fname, const std::string& contents)
{
    write_file(dir / fname, contents);
    const std::filesystem::path log = dir / (fname + ".log");
    std::ostringstream cmd;
    cmd << "\"" << self << "\" --child \"" << dir.string() << "/\" \"" << fname
        << "\" > \"" << log.string() << "\" 2>&1";
    child_run out;
    const int rc = std::system(cmd.str().c_str());
#ifdef _WIN32
    out.exit_code = rc;
#else
    out.exit_code = (rc == -1) ? -1
                  : (WIFEXITED(rc) ? WEXITSTATUS(rc) : 128 + WTERMSIG(rc));
#endif
    std::ifstream ifs(log);
    std::ostringstream buf;
    buf << ifs.rdbuf();
    out.output = buf.str();
    return out;
}

bool contains(const std::string& haystack, const std::string& needle)
{
    return haystack.find(needle) != std::string::npos;
}

const std::string k_clean = "write_fname: strict\n";

// One bad line: the run must stop, before any value is used, with every
// expected fragment in the output.
void expect_stop(const std::string& self, const std::filesystem::path& dir,
                 const std::string& name, const std::string& line,
                 const std::vector<std::string>& fragments)
{
    const child_run r = run_case(self, dir, name + ".param", k_clean + line + "\n");
    bool all = true;
    for (const auto& f : fragments) all = all && contains(r.output, f);
    check(r.exit_code == 1 && !contains(r.output, "CHILD_COMPLETED") && all,
          "stops: " + line);
    if (!all) std::cout << r.output << "\n";
}

} // namespace

int main(int argc, char* argv[])
{
    if (argc == 4 && std::string(argv[1]) == "--child")
        return run_child(argv[2], argv[3]);
    if (argc < 2)
    {
        std::cerr << "usage: lsdtt-test-param-parser-strict <work_dir>\n";
        return 2;
    }
    // Resolved before any case changes the working directory.
    const std::string self = std::filesystem::absolute(argv[0]).string();
    const std::filesystem::path work =
        std::filesystem::absolute(argv[1]);
    std::filesystem::remove_all(work);
    const std::filesystem::path sub = work / "runs" / "a";
    write_file(sub / "run.param", "write_fname: strict\nn_threads: 1\n");

    // --- location ------------------------------------------------------------
    {
        const std::string got =
            resolved_read_path("", (sub / "run.param").string());
        check(same_dir(got, sub),
              "one argument, absolute file: read_path is the file's directory "
              "(got '" + got + "')");
    }
    {
        const std::string got =
            resolved_read_path("/somewhere/else", (sub / "run.param").string());
        check(same_dir(got, sub),
              "absolute file wins over a <path> argument (got '" + got + "')");
    }
    {
        const std::string got = resolved_read_path(work.string(), "runs/a/run.param");
        check(same_dir(got, sub),
              "<path> plus a relative file with a directory: read_path is "
              "<path>/<dir> (got '" + got + "')");
    }
    {
        const std::string got = resolved_read_path(sub.string(), "run.param");
        check(same_dir(got, sub),
              "<path> plus a bare filename is unchanged (got '" + got + "')");
    }
    {
        std::filesystem::current_path(work);
        const std::string got = resolved_read_path("", "runs/a/run.param");
        check(same_dir(got, sub),
              "one argument, relative file with a directory: read_path is "
              "that directory (got '" + got + "')");
    }
    {
        const std::string got =
            resolved_read_path("", (work / "missing" / "run.param").string());
        check(got.rfind("THREW: ", 0) == 0 &&
                  got.find((work / "missing").string()) != std::string::npos,
              "a missing absolute file is reported by its real path");
    }

    // --- silence -------------------------------------------------------------
    {
        std::ostringstream captured;
        std::streambuf* old = std::cout.rdbuf(captured.rdbuf());
        { lsdtt::lsdtt_param_parser placeholder; }
        std::cout.rdbuf(old);
        check(captured.str().empty(),
              "the default (placeholder) constructor prints nothing");
    }

    // --- strict values ------------------------------------------------------
    const std::filesystem::path cases = work / "cases";
    std::filesystem::create_directories(cases);
    {
        const child_run r = run_case(self, cases, "clean.param", k_clean +
            "erodibility: 1e-4\nn_iterations: +250\nwrite_fill_raster: T\n"
            "write_slope_raster: false\ntolerance: 2.5\nmethod: d8\nK:\n");
        check(r.exit_code == 0 && contains(r.output, "CHILD_COMPLETED"),
              "a clean file parses");
        check(contains(r.output, "CHILD erodibility=0.0001") &&
                  contains(r.output, "CHILD n_iterations=250") &&
                  contains(r.output, "CHILD write_fill_raster=1") &&
                  contains(r.output, "CHILD write_slope_raster=0") &&
                  contains(r.output, "CHILD tolerance=2.5"),
              "  accepted spellings are read correctly (1e-4, +250, T, false), "
              "and an empty string value is allowed");
        check(!contains(r.output, "FATAL"), "  with no FATAL line");
    }
    expect_stop(self, cases, "bool_yes", "write_fill_raster: yes",
                {"FATAL: parameter 'write_fill_raster' value 'yes' is not true "
                 "or false"});
    expect_stop(self, cases, "int_exponent", "n_iterations: 1e3",
                {"value '1e3' is not a whole number",
                 "Write it as a whole number: 1000."});
    expect_stop(self, cases, "int_fraction", "n_iterations: 1.5",
                {"value '1.5' is not a whole number", "e.g. 1000"});
    expect_stop(self, cases, "int_units", "n_iterations: 100m",
                {"value '100m' is not a whole number"});
    expect_stop(self, cases, "int_range", "n_iterations: 99999999999",
                {"is out of range"});
    expect_stop(self, cases, "real_units", "erodibility: 100m",
                {"value '100m' is not a number", "with no units"});
    expect_stop(self, cases, "real_nan", "tolerance: nan",
                {"value 'nan' is not a number"});
    expect_stop(self, cases, "real_float_range", "erodibility: 1e300",
                {"value '1e300' is out of range"});
    expect_stop(self, cases, "empty_int", "n_iterations:",
                {"parameter 'n_iterations' has no value"});
    expect_stop(self, cases, "empty_bool", "write_fill_raster:",
                {"parameter 'write_fill_raster' has no value", "true or false"});
    expect_stop(self, cases, "quoted", "method: \"d8\"",
                {"value \"d8\" is quoted", "Write d8."});
    expect_stop(self, cases, "misspelt", "write_fil_raster: true",
                {"parameter 'write_fil_raster' is not a parameter of this tool",
                 "Did you mean 'write_fill_raster'?"});
    expect_stop(self, cases, "misspelt_case", "Erodibility: 2",
                {"Did you mean 'erodibility'?"});
    expect_stop(self, cases, "misspelt_reserved", "write_fnme: x",
                {"Did you mean 'write_fname'?"});
    {
        // A short unknown key is not "corrected" to the one-letter key K.
        const child_run r = run_case(self, cases, "short_unknown.param",
                                     k_clean + "m: 0.5\n");
        check(r.exit_code == 0 && contains(r.output, "CHILD_COMPLETED") &&
                  contains(r.output, "are not used by this tool"),
              "a short unrelated key is a warning, not a near miss of K");
    }
    {
        const child_run r = run_case(self, cases, "several.param", k_clean +
            "write_fill_raster: yes\nn_iterations: 1.5\nmethd: d8\n");
        check(r.exit_code == 1 && contains(r.output, "3 problem(s) in several.param"),
              "every problem in the file is reported in one run");
    }

    std::cout << "\n" << g_n_checks - g_n_failed << " / " << g_n_checks
              << " checks passed\n";
    return g_n_failed == 0 ? 0 : 1;
}
