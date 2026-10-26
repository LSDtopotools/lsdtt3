// lsdtt_run_record.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// See lsdtt_run_record.hpp.

#include "lsdtt_run_record.hpp"

#include "lsdtt_gdal_init.hpp"
#include "lsdtt_param_parser.hpp"

#include <algorithm>
#include <cstdlib>
#include <ctime>
#include <deque>
#include <exception>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <streambuf>

namespace lsdtt {

namespace {

lsdtt_run_record* g_live = nullptr;
bool g_handlers_installed = false;
std::terminate_handler g_previous_terminate = nullptr;

constexpr std::size_t k_kept_lines     = 400;
constexpr std::size_t k_max_box_lines  = 15;
const char* const     k_frame =
    "=========================================================";

bool starts_with(const std::string& s, const char* prefix)
{
    return s.rfind(prefix, 0) == 0;
}

bool is_error_line(const std::string& s)
{
    return starts_with(s, "FATAL") || starts_with(s, "ERROR") ||
           starts_with(s, "Error") || starts_with(s, "error:");
}

bool is_warning_line(const std::string& s)
{
    return starts_with(s, "WARNING") || starts_with(s, "Warning");
}

std::string json_string(const std::string& s)
{
    std::ostringstream out;
    out << '"';
    for (const unsigned char c : s)
    {
        switch (c)
        {
            case '"':  out << "\\\""; break;
            case '\\': out << "\\\\"; break;
            case '\n': out << "\\n";  break;
            case '\r': out << "\\r";  break;
            case '\t': out << "\\t";  break;
            default:
                if (c < 0x20)
                    out << "\\u" << std::hex << std::setw(4) << std::setfill('0')
                        << static_cast<int>(c) << std::dec << std::setfill(' ');
                else
                    out << c;
        }
    }
    out << '"';
    return out.str();
}

std::string iso_utc(std::chrono::system_clock::time_point t)
{
    const std::time_t tt = std::chrono::system_clock::to_time_t(t);
    std::tm utc{};
#ifdef _WIN32
    gmtime_s(&utc, &tt);
#else
    gmtime_r(&tt, &utc);
#endif
    std::ostringstream out;
    out << std::put_time(&utc, "%Y-%m-%dT%H:%M:%SZ");
    return out.str();
}

std::string absolute_or_given(const std::string& p)
{
    if (p.empty()) return p;
    std::error_code ec;
    const auto a = std::filesystem::absolute(p, ec);
    return ec ? p : a.lexically_normal().string();
}

} // namespace

// -----------------------------------------------------------------------------
// line_capture: passes everything to the real std::cerr buffer and keeps the
// last k_kept_lines complete lines.
// -----------------------------------------------------------------------------
class lsdtt_run_record::line_capture : public std::streambuf
{
public:
    explicit line_capture(std::streambuf* target) : target_(target) {}

    std::streambuf* target() const { return target_; }

    std::vector<std::string> lines() const
    {
        std::vector<std::string> out(lines_.begin(), lines_.end());
        if (!current_.empty()) out.push_back(current_);
        return out;
    }

protected:
    int_type overflow(int_type ch) override
    {
        if (traits_type::eq_int_type(ch, traits_type::eof()))
            return traits_type::not_eof(ch);
        keep(traits_type::to_char_type(ch));
        return target_->sputc(traits_type::to_char_type(ch));
    }

    std::streamsize xsputn(const char* s, std::streamsize n) override
    {
        for (std::streamsize i = 0; i < n; ++i) keep(s[i]);
        return target_->sputn(s, n);
    }

    int sync() override { return target_->pubsync(); }

private:
    void keep(char c)
    {
        if (c == '\n')
        {
            lines_.push_back(current_);
            current_.clear();
            if (lines_.size() > k_kept_lines) lines_.pop_front();
        }
        else if (c != '\r')
        {
            current_ += c;
        }
    }

    std::streambuf*         target_;
    std::deque<std::string> lines_;
    std::string             current_;
};

// -----------------------------------------------------------------------------
// Lifetime
// -----------------------------------------------------------------------------
lsdtt_run_record::lsdtt_run_record(const std::string& tool_name,
                                   int argc, char* argv[])
    : tool_name_(tool_name),
      start_wall_(std::chrono::system_clock::now()),
      start_steady_(std::chrono::steady_clock::now()),
      start_file_time_(std::filesystem::file_time_type::clock::now())
{
    for (int i = 0; i < argc; ++i)
    {
        if (i) command_line_ += ' ';
        command_line_ += argv[i];
    }

    capture_ = std::make_unique<line_capture>(std::cerr.rdbuf());
    std::cerr.rdbuf(capture_.get());

    g_live = this;
    if (!g_handlers_installed)
    {
        // std::exit (the parser's FATAL path) skips this object's destructor
        // but runs atexit handlers; an uncaught exception runs neither.
        std::atexit(&lsdtt_run_record::on_exit);
        g_previous_terminate = std::set_terminate(&lsdtt_run_record::on_terminate);
        g_handlers_installed = true;
    }
}

lsdtt_run_record::~lsdtt_run_record()
{
    finish();
    if (g_live == this) g_live = nullptr;
}

void lsdtt_run_record::on_exit()
{
    if (g_live) g_live->finish();
}

void lsdtt_run_record::on_terminate()
{
    if (g_live)
    {
        std::string what = "the run stopped on an unexpected error";
        if (const std::exception_ptr e = std::current_exception())
        {
            try { std::rethrow_exception(e); }
            catch (const std::exception& ex) { what = std::string("uncaught error: ") + ex.what(); }
            catch (...) { what = "uncaught error of unknown type"; }
        }
        g_live->extra_errors_.push_back(what);
        g_live->finish();
    }
    if (g_previous_terminate) g_previous_terminate();
    std::abort();
}

// -----------------------------------------------------------------------------
// Static hooks
// -----------------------------------------------------------------------------
void lsdtt_run_record::note_paths(const std::string& parameter_file,
                                  const std::string& read_path,
                                  const std::string& read_fname,
                                  const std::string& write_path,
                                  const std::string& write_fname)
{
    if (!g_live) return;
    g_live->paths_ = {absolute_or_given(parameter_file),
                      absolute_or_given(read_path), read_fname,
                      absolute_or_given(write_path), write_fname, true};
}

void lsdtt_run_record::add_output(const std::string& path)
{
    if (g_live) g_live->extra_outputs_.push_back(path);
}

std::string lsdtt_run_record::summary_filename(const std::string& write_fname)
{
    return write_fname + "_run_summary.json";
}

// -----------------------------------------------------------------------------
// What the run did
// -----------------------------------------------------------------------------
std::vector<lsdtt_run_record::output_file> lsdtt_run_record::find_outputs() const
{
    namespace fs = std::filesystem;
    std::vector<output_file> out;
    if (!paths_.known || paths_.write_path.empty()) return out;

    const fs::path root(paths_.write_path);
    std::error_code ec;
    if (!fs::is_directory(root, ec)) return out;

    // A one-second allowance for file systems with coarse time stamps.
    const auto since = start_file_time_ - std::chrono::seconds(1);
    const std::string prefix = paths_.write_fname + "_";
    const std::string summary = summary_filename(paths_.write_fname);

    auto consider = [&](const fs::path& file)
    {
        std::error_code e;
        if (!fs::is_regular_file(file, e)) return;
        if (fs::last_write_time(file, e) < since || e) return;
        output_file f;
        f.path  = fs::relative(file, root, e).generic_string();
        if (e || f.path.empty()) f.path = file.generic_string();
        f.bytes = fs::file_size(file, e);
        if (f.path == summary) return;
        out.push_back(f);
    };

    for (const auto& entry : fs::directory_iterator(root, ec))
    {
        const std::string name = entry.path().filename().string();
        if (!paths_.write_fname.empty() && name.rfind(prefix, 0) != 0) continue;
        std::error_code e;
        if (entry.is_directory(e))
        {
            for (const auto& inner : fs::recursive_directory_iterator(entry.path(), e))
                consider(inner.path());
        }
        else
        {
            consider(entry.path());
        }
    }
    for (const std::string& extra : extra_outputs_)
    {
        fs::path p(extra);
        if (p.is_relative()) p = root / p;
        consider(p);
    }

    std::sort(out.begin(), out.end(),
              [](const output_file& a, const output_file& b) { return a.path < b.path; });
    out.erase(std::unique(out.begin(), out.end(),
                          [](const output_file& a, const output_file& b)
                          { return a.path == b.path; }),
              out.end());
    return out;
}

// The block of error lines nearest the end of what the run printed: the last
// FATAL/ERROR/Error line, with the error and indented lines around it.
std::vector<std::string> lsdtt_run_record::error_lines() const
{
    const std::vector<std::string> lines = capture_ ? capture_->lines()
                                                    : std::vector<std::string>{};
    std::vector<std::string> out;
    std::size_t last = lines.size();
    for (std::size_t i = lines.size(); i-- > 0;)
        if (is_error_line(lines[i])) { last = i; break; }

    if (last < lines.size())
    {
        auto indented = [&](std::size_t i)
        {
            return !lines[i].empty() && (lines[i][0] == ' ' || lines[i][0] == '\t');
        };
        // Back over earlier error lines and their indented continuations, but
        // not over the indented lines of a warning that came just before.
        std::size_t first = last;
        for (std::size_t i = last; i-- > 0;)
        {
            if (is_error_line(lines[i])) first = i;
            else if (!indented(i)) break;
        }
        std::size_t end = last + 1;
        while (end < lines.size() && (is_error_line(lines[end]) || indented(end))) ++end;
        out.assign(lines.begin() + static_cast<std::ptrdiff_t>(first),
                   lines.begin() + static_cast<std::ptrdiff_t>(end));
    }
    else
    {
        // No recognisable error line: the last few lines printed to stderr.
        for (std::size_t i = lines.size(); i-- > 0 && out.size() < 5;)
            if (!lines[i].empty() && !is_warning_line(lines[i]))
                out.insert(out.begin(), lines[i]);
    }
    out.insert(out.end(), extra_errors_.begin(), extra_errors_.end());
    if (out.empty()) out.push_back("(no error message was printed)");
    return out;
}

std::vector<std::string> lsdtt_run_record::warning_lines() const
{
    // Each warning with its indented continuation lines (the ignored key
    // names, for instance), joined into one entry.
    std::vector<std::string> out;
    if (!capture_) return out;
    bool in_warning = false;
    for (const std::string& line : capture_->lines())
    {
        if (is_warning_line(line))
        {
            out.push_back(line);
            in_warning = true;
        }
        else if (in_warning && !line.empty() && (line[0] == ' ' || line[0] == '\t'))
        {
            out.back() += "\n" + line;
        }
        else
        {
            in_warning = false;
        }
    }
    return out;
}

// Returns the summary's path, or "" if it could not be written.
std::string lsdtt_run_record::write_summary(const std::vector<output_file>& outputs,
                                            const std::vector<std::string>& errors,
                                            const std::vector<std::string>& warnings) const
{
    namespace fs = std::filesystem;
    if (!paths_.known || paths_.write_path.empty()) return "";
    std::error_code ec;
    if (!fs::is_directory(paths_.write_path, ec)) return "";

    const fs::path path = fs::path(paths_.write_path) /
                          summary_filename(paths_.write_fname);
    std::ofstream out(path);
    if (!out) return "";

    const auto end_wall = std::chrono::system_clock::now();
    const double seconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - start_steady_).count();
    const std::string ingested = paths_.write_fname + "_ingested_params.param";
    const bool has_ingested =
        fs::exists(fs::path(paths_.write_path) / ingested, ec);

    auto string_list = [&](const std::vector<std::string>& items)
    {
        std::string s = "[";
        for (std::size_t i = 0; i < items.size(); ++i)
            s += (i ? ", " : "") + json_string(items[i]);
        return s + "]";
    };

    out << "{\n"
        << "  \"tool\": " << json_string(tool_name_) << ",\n"
        << "  \"lsdtt3_version\": " << json_string(lsdtt3_package_version()) << ",\n"
        << "  \"gdal_version\": " << json_string(gdal_version()) << ",\n"
        << "  \"proj_version\": " << json_string(proj_version()) << ",\n"
        << "  \"command_line\": " << json_string(command_line_) << ",\n"
        << "  \"parameter_file\": " << json_string(paths_.parameter_file) << ",\n"
        << "  \"read_path\": " << json_string(paths_.read_path) << ",\n"
        << "  \"read_fname\": " << json_string(paths_.read_fname) << ",\n"
        << "  \"write_path\": " << json_string(paths_.write_path) << ",\n"
        << "  \"write_fname\": " << json_string(paths_.write_fname) << ",\n"
        << "  \"start_utc\": " << json_string(iso_utc(start_wall_)) << ",\n"
        << "  \"end_utc\": " << json_string(iso_utc(end_wall)) << ",\n"
        << "  \"wall_time_seconds\": " << std::fixed << std::setprecision(3)
        << seconds << std::defaultfloat << ",\n"
        << "  \"status\": " << json_string(succeeded_ ? "ok" : "failed") << ",\n"
        << "  \"errors\": " << string_list(errors) << ",\n"
        << "  \"warnings\": " << string_list(warnings) << ",\n"
        << "  \"ingested_parameters_file\": "
        << (has_ingested ? json_string(ingested) : std::string("null")) << ",\n"
        << "  \"outputs\": [";
    for (std::size_t i = 0; i < outputs.size(); ++i)
    {
        out << (i ? ",\n" : "\n") << "    {\"path\": " << json_string(outputs[i].path)
            << ", \"bytes\": " << outputs[i].bytes << "}";
    }
    out << (outputs.empty() ? "]\n" : "\n  ]\n") << "}\n";
    out.close();
    return out ? path.string() : "";
}

// -----------------------------------------------------------------------------
// The end of the run
// -----------------------------------------------------------------------------
void lsdtt_run_record::finish()
{
    if (finished_) return;
    finished_ = true;

    // Stop capturing before printing anything of our own.
    if (capture_ && std::cerr.rdbuf() == capture_.get())
        std::cerr.rdbuf(capture_->target());

    const std::vector<output_file> outputs = find_outputs();
    const std::vector<std::string> errors =
        succeeded_ ? std::vector<std::string>{} : error_lines();
    const std::vector<std::string> warnings = warning_lines();
    const std::string summary = write_summary(outputs, errors, warnings);

    std::cerr.flush();
    if (succeeded_)
    {
        std::cout << "lsdtt3: " << tool_name_ << " finished OK; " << outputs.size()
                  << " output file" << (outputs.size() == 1 ? "" : "s") << "; summary: "
                  << (summary.empty() ? std::string("not written") : summary) << std::endl;
        return;
    }

    std::cout << "\n" << k_frame << "\n"
              << "|| RUN FAILED: " << tool_name_ << "\n";
    const std::size_t shown = std::min(errors.size(), k_max_box_lines);
    for (std::size_t i = 0; i < shown; ++i)
        std::cout << "|| " << errors[i] << "\n";
    if (errors.size() > shown)
        std::cout << "|| ... " << errors.size() - shown << " more line(s) above\n";
    if (summary.empty())
        std::cout << "|| No summary file was written (no usable write_path).\n";
    else
        std::cout << "|| Summary: " << summary << "\n";
    std::cout << k_frame << std::endl;

    std::cerr << "ERROR: " << tool_name_ << " did not finish; "
              << (summary.empty() ? std::string("no summary file was written")
                                  : "summary: " + summary)
              << std::endl;
}

} // namespace lsdtt
