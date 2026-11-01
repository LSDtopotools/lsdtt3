// lsdtt_run_record.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// lsdtt_run_record: what a CLI run did, for a person and for a program.
//
// Every run ends in one of two ways:
//
//   success   the last stdout line is
//               lsdtt3: <tool> finished OK; <n> output files; summary: <path>
//   failure   a boxed RUN FAILED message on STDOUT, with the error lines the
//             run printed, plus a short ERROR line on stderr. A typical user
//             watches stdout and never looks at stderr, so the failure has to
//             be loud where they are looking.
//
// and, where write_path is usable, writes <write_fname>_run_summary.json:
// tool, versions, command line, resolved paths, start/end time, status, the
// error and warning lines, and every output file the run wrote. A coding agent
// driving the binaries reads that file instead of scraping a log.
//
// HOW IT WORKS WITHOUT TOUCHING ANY FAILURE SITE. The CLIs have about 1,100
// failure returns and 850 FATAL messages, all on stderr. The record does not
// need any of them changed:
//
//   * it puts a pass-through buffer on std::cerr that keeps the lines the run
//     prints, and takes the error lines for the box and the summary from them;
//   * the CLI ends each successful path with `return run_record.succeed();`;
//     a run that ends any other way (return 1, std::exit from the parser, an
//     uncaught exception) is a failure;
//   * the parameter parser reports its resolved paths to the live record
//     (note_paths), so no CLI has to;
//   * outputs are found after the run: every file under write_path, at any
//     depth, inside a top-level entry named <write_fname>_..., modified since
//     the run started. A tool whose output does not follow that naming
//     registers it with add_output().
//
// Use, in a CLI's main, after the help branch:
//
//     lsdtt::lsdtt_run_record run_record(k_tool_name, argc, argv);
//     ...
//     return run_record.succeed();
//
// One record per process. The record writes nothing for help, -v or -c, which
// return before it exists.
// -----------------------------------------------------------------------------

#pragma once

#include <chrono>
#include <cstdint>
#include <filesystem>
#include <memory>
#include <string>
#include <vector>
#include <utility>

namespace lsdtt {

class lsdtt_run_record
{
public:
    lsdtt_run_record(const std::string& tool_name, int argc, char* argv[]);
    ~lsdtt_run_record();

    lsdtt_run_record(const lsdtt_run_record&)            = delete;
    lsdtt_run_record& operator=(const lsdtt_run_record&) = delete;

    // The run finished. Returns 0, so a CLI ends with
    //     return run_record.succeed();
    // which stays one statement wherever the return sits.
    int succeed() { succeeded_ = true; return 0; }

    // Called by lsdtt_param_parser once it has resolved its paths. A no-op
    // when no record is live (unit tests, bindings).
    static void note_paths(const std::string& parameter_file,
                           const std::string& read_path,
                           const std::string& read_fname,
                           const std::string& write_path,
                           const std::string& write_fname);

    // An output whose name does not start with <write_fname>_ (the naming the
    // record finds outputs by). A no-op when no record is live.
    static void add_output(const std::string& path);

    // A provenance item (key -> value) for the summary's "provenance" object,
    // e.g. which data source a DEM came from. The object is written only when
    // at least one item was added, so other CLIs' summaries are unchanged.
    // A later value for the same key replaces the earlier one. A no-op when no
    // record is live.
    static void add_provenance(const std::string& key, const std::string& value);

    // Name of the summary file for a given write_fname.
    static std::string summary_filename(const std::string& write_fname);

    // Write the summary and the closing message now (once). Called by the
    // destructor, at std::exit and from the terminate handler.
    void finish();

private:
    struct paths
    {
        std::string parameter_file, read_path, read_fname, write_path, write_fname;
        bool known = false;
    };
    struct output_file
    {
        std::string path;        // relative to write_path when inside it
        std::uintmax_t bytes = 0;
    };

    std::vector<output_file> find_outputs() const;
    std::vector<std::string> error_lines() const;
    std::vector<std::string> warning_lines() const;
    std::string write_summary(const std::vector<output_file>& outputs,
                              const std::vector<std::string>& errors,
                              const std::vector<std::string>& warnings) const;

    std::string tool_name_;
    std::string command_line_;
    std::chrono::system_clock::time_point     start_wall_;
    std::chrono::steady_clock::time_point     start_steady_;
    std::filesystem::file_time_type           start_file_time_;
    paths       paths_;
    std::vector<std::string> extra_outputs_;
    std::vector<std::pair<std::string, std::string>> provenance_;
    std::vector<std::string> extra_errors_;    // e.g. an uncaught exception
    bool succeeded_ = false;
    bool finished_  = false;

    class line_capture;
    std::unique_ptr<line_capture> capture_;

    static void on_exit();
    static void on_terminate();
};

} // namespace lsdtt
