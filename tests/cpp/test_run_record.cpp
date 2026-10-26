// test_run_record.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// unit_run_record
//
// lsdtt_run_record decides how every CLI run ends: the finished-OK line and a
// run summary on success; a boxed RUN FAILED message on stdout, an ERROR line
// on stderr and a failed summary otherwise. A run can end four ways, and each
// is a different code path in the record:
//
//   ok         return run_record.succeed()      (destructor)
//   returned   return 1 after a FATAL message   (destructor)
//   exited     std::exit(1), as the parser does (atexit handler)
//   threw      an exception escaping main       (terminate handler)
//
// Each runs in a CHILD process (the driver re-executes itself), and the parent
// checks the exit code, stdout, stderr and the summary file. The ok case also
// checks output discovery: a file named <write_fname>_..., a file inside a
// <write_fname>_... subdirectory, a registered file with another name, and
// that a pre-existing file and an unrelated new file are NOT listed.
// -----------------------------------------------------------------------------

#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

#ifndef _WIN32
#include <sys/wait.h>
#endif

#include "../../src/core/lsdtt_run_record.hpp"

namespace fs = std::filesystem;

namespace
{

int g_n_checks = 0;
int g_n_failed = 0;

void check(bool condition, const std::string& what)
{
    ++g_n_checks;
    std::cout << (condition ? "PASS  " : "FAIL  ") << what << "\n";
    if (!condition) ++g_n_failed;
}

bool contains(const std::string& haystack, const std::string& needle)
{
    return haystack.find(needle) != std::string::npos;
}

std::string slurp(const fs::path& p)
{
    std::ifstream in(p);
    std::ostringstream s;
    s << in.rdbuf();
    return s.str();
}

void touch(const fs::path& p, const std::string& text = "x")
{
    fs::create_directories(p.parent_path());
    std::ofstream(p) << text;
}

// -----------------------------------------------------------------------------
// child
// -----------------------------------------------------------------------------
int run_child(const std::string& mode, const fs::path& dir, int argc, char* argv[])
{
    lsdtt::lsdtt_run_record run_record("lsdtt-test-tool", argc, argv);
    lsdtt::lsdtt_run_record::note_paths((dir / "run.param").string(), dir.string(),
                                        "dem", dir.string(), "rr");
    std::cerr << "WARNING: 1 parameter(s) in run.param are not used\n"
              << "  frobnicate\n";

    if (mode == "ok")
    {
        touch(dir / "rr_hillshade.tif");
        touch(dir / "rr_geodata" / "GeoData.txt");
        touch(dir / "dem_rotated_back.tif");
        lsdtt::lsdtt_run_record::add_output((dir / "dem_rotated_back.tif").string());
        touch(dir / "unrelated.csv");
        return run_record.succeed();
    }
    if (mode == "returned")
    {
        std::cerr << "FATAL: the DEM has no valid cells\n"
                  << "  Check read_fname.\n";
        return 1;
    }
    if (mode == "exited")
    {
        std::cerr << "FATAL: parameter 'n' value '1.5' is not a whole number.\n";
        std::exit(1);
    }
    if (mode == "threw")
    {
        throw std::runtime_error("disk full while writing rr_slope.tif");
    }
    return 3;
}

// -----------------------------------------------------------------------------
// parent
// -----------------------------------------------------------------------------
struct child_run
{
    int exit_code = 0;
    std::string out, err;
};

child_run run_case(const std::string& self, const std::string& mode, const fs::path& dir)
{
    fs::remove_all(dir);
    fs::create_directories(dir);
    touch(dir / "rr_old_output.csv");   // written before the run started
    fs::last_write_time(dir / "rr_old_output.csv",
                        fs::file_time_type::clock::now() - std::chrono::hours(1));
    const fs::path out = dir.parent_path() / (mode + ".out");
    const fs::path err = dir.parent_path() / (mode + ".err");
    std::ostringstream cmd;
    cmd << "\"" << self << "\" --child " << mode << " \"" << dir.string() << "\" > \""
        << out.string() << "\" 2> \"" << err.string() << "\"";
    child_run r;
    const int rc = std::system(cmd.str().c_str());
#ifdef _WIN32
    r.exit_code = rc;
#else
    r.exit_code = (rc == -1) ? -1 : (WIFEXITED(rc) ? WEXITSTATUS(rc) : 128 + WTERMSIG(rc));
#endif
    r.out = slurp(out);
    r.err = slurp(err);
    return r;
}

} // namespace

int main(int argc, char* argv[])
{
    if (argc >= 4 && std::string(argv[1]) == "--child")
        return run_child(argv[2], argv[3], argc, argv);
    if (argc < 2)
    {
        std::cerr << "usage: lsdtt-test-run-record <work_dir>\n";
        return 2;
    }
    const std::string self = fs::absolute(argv[0]).string();
    const fs::path work = fs::absolute(argv[1]);

    // --- ok -------------------------------------------------------------------
    {
        const fs::path dir = work / "ok";
        const child_run r = run_case(self, "ok", dir);
        const std::string summary = slurp(dir / "rr_run_summary.json");
        check(r.exit_code == 0, "ok: exit code 0");
        check(contains(r.out, "lsdtt3: lsdtt-test-tool finished OK; 3 output files; summary: "),
              "ok: the last line reports success and three outputs");
        check(!contains(r.out, "RUN FAILED") && !contains(r.err, "ERROR:"),
              "ok: no failure box, no ERROR line");
        check(contains(summary, "\"status\": \"ok\""), "ok: summary status ok");
        check(contains(summary, "\"path\": \"rr_hillshade.tif\"") &&
                  contains(summary, "\"path\": \"rr_geodata/GeoData.txt\"") &&
                  contains(summary, "\"path\": \"dem_rotated_back.tif\""),
              "ok: lists a prefixed file, a file in a prefixed directory and a "
              "registered file");
        check(!contains(summary, "unrelated.csv") && !contains(summary, "rr_old_output.csv") &&
                  !contains(summary, "rr_run_summary.json"),
              "ok: does not list an unrelated file, an older file or itself");
        check(contains(summary, "\"warnings\": [\"WARNING: 1 parameter(s) in run.param "
                                "are not used\\n  frobnicate\"]"),
              "ok: a warning keeps its indented continuation line");
    }

    // --- the three ways to fail -------------------------------------------------
    struct failure { const char* mode; int exit_code; const char* error; };
    for (const failure& f : {failure{"returned", 1, "FATAL: the DEM has no valid cells"},
                             failure{"exited", 1, "is not a whole number"},
                             failure{"threw", -1, "uncaught error: disk full"}})
    {
        const std::string m = f.mode;
        const fs::path dir = work / m;
        const child_run r = run_case(self, m, dir);
        const std::string summary = slurp(dir / "rr_run_summary.json");
        check(f.exit_code < 0 ? r.exit_code != 0 : r.exit_code == f.exit_code,
              m + ": non-zero exit code");
        check(contains(r.out, "|| RUN FAILED: lsdtt-test-tool") &&
                  contains(r.out, std::string("|| ") + "") && contains(r.out, f.error),
              m + ": stdout has the RUN FAILED box with the error");
        check(contains(r.out, "|| Summary: "), m + ": the box names the summary");
        check(contains(r.err, "ERROR: lsdtt-test-tool did not finish; summary: "),
              m + ": stderr has the ERROR line");
        check(contains(summary, "\"status\": \"failed\"") && contains(summary, f.error),
              m + ": summary status failed, with the error");
        check(!contains(summary, "frobnicate\", \"FATAL") &&
                  !contains(summary, "\"errors\": [\"WARNING"),
              m + ": the warning before the error is not taken as the error");
    }
    {
        const child_run r = run_case(self, "returned", work / "returned");
        check(contains(r.out, "|| FATAL: the DEM has no valid cells\n||   Check read_fname."),
              "returned: the error's indented continuation is in the box");
    }

    std::cout << "\n" << g_n_checks - g_n_failed << " / " << g_n_checks
              << " checks passed\n";
    return g_n_failed == 0 ? 0 : 1;
}
