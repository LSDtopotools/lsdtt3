// lsdtt_param_parser.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html
//
// Original LSDTopoTools implementation by Simon M. Mudd and contributors.
// Redesigned for lsdtt3: removed extension logic, removed hardcoded file
// parameters, renamed keys to snake_case, removed external dependencies.
// This reads parameter files (which we sometimes call driver files)
// format is
// parameter_name: value
// lines with hash at beginning are treated as comments
#include "lsdtt_param_parser.hpp"

#include <algorithm>
#include <cctype>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>

// Package-level version injected by CMake (release: "X.Y.Z"; development
// build: "X.Y.Z-dev, git <describe>"). Distinct from each CLI's k_version.
#ifndef LSDTT3_VERSION
#define LSDTT3_VERSION "unknown"
#endif

namespace lsdtt {

// -----------------------------------------------------------------------------
// Construction
// -----------------------------------------------------------------------------

lsdtt_param_parser::lsdtt_param_parser()
{
    std::cout << "lsdtt_param_parser: created empty parser. "
              << "Provide a path and filename to parse a parameter file.\n";
}

lsdtt_param_parser::lsdtt_param_parser(
    const std::string& path,
    const std::string& filename)
{
    param_file_path_ = fix_path(path);
    param_fname_     = filename;

    std::string full_path = param_file_path_ + param_fname_;

    // Verify the file exists before parsing
    std::ifstream test(full_path);
    if (!test.good())
    {
        throw std::runtime_error(
            "lsdtt_param_parser: parameter file not found: " + full_path
        );
    }
    test.close();

    std::cout << "lsdtt_param_parser: parsing " << full_path << "\n";
    parse_file(full_path);
    parse_file_io();
    check_paths();
}

// -----------------------------------------------------------------------------
// Parsing
// -----------------------------------------------------------------------------

void lsdtt_param_parser::parse_file(const std::string& full_path)
{
    std::ifstream infile(full_path);
    if (!infile.good())
    {
        throw std::runtime_error(
            "lsdtt_param_parser::parse_file: cannot open " + full_path
        );
    }

    std::map<std::string, std::string> temp_map;
    std::string parameter, value;

    while (infile.good())
    {
        parse_line(infile, parameter, value);
        if (parameter == "NULL") continue;

        value = strip_control_chars(value);
        temp_map[parameter] = value;
    }

    parameter_map_ = temp_map;
}

void lsdtt_param_parser::parse_file_io()
{
    // read_path
    if (parameter_map_.count("read_path"))
    {
        read_path_ = strip_control_chars(parameter_map_["read_path"]);
    }
    else
    {
        read_path_ = param_file_path_;
    }

    // write_path
    if (parameter_map_.count("write_path"))
    {
        write_path_ = strip_control_chars(parameter_map_["write_path"]);
    }

    // read_fname
    if (parameter_map_.count("read_fname"))
    {
        read_fname_ = strip_control_chars(parameter_map_["read_fname"]);
    }

    // write_fname
    if (parameter_map_.count("write_fname"))
    {
        write_fname_ = strip_control_chars(parameter_map_["write_fname"]);
    }

    // boundary_conditions
    if (parameter_map_.count("boundary_conditions"))
    {
        const std::string& bc_str = parameter_map_["boundary_conditions"];
        boundary_conditions_.clear();
        for (char c : bc_str)
        {
            std::string this_bc(1, static_cast<char>(
                std::tolower(static_cast<unsigned char>(c))
            ));
            if (this_bc == "o" || this_bc == "n" || this_bc == "p")
            {
                boundary_conditions_.push_back(this_bc);
            }
            else
            {
                std::cout << "lsdtt_param_parser: '" << this_bc
                          << "' is not a valid boundary condition "
                          << "(valid: o, n, p). Defaulting to open (o).\n";
                boundary_conditions_.push_back("o");
            }
        }
    }
}

void lsdtt_param_parser::parse_all_parameters(
    const std::map<std::string, float>&       defaults_f,
    const std::map<std::string, int>&         defaults_i,
    const std::map<std::string, bool>&        defaults_b,
    const std::map<std::string, std::string>& defaults_s
)
{
    parse_float_parameters(defaults_f);
    parse_int_parameters(defaults_i);
    parse_bool_parameters(defaults_b);
    parse_string_parameters(defaults_s);
}

void lsdtt_param_parser::parse_all_parameters(
    const std::map<std::string, float>&       defaults_f,
    const std::map<std::string, int>&         defaults_i,
    const std::map<std::string, bool>&        defaults_b,
    const std::map<std::string, std::string>& defaults_s,
    const std::map<std::string, double>&      defaults_d
)
{
    parse_float_parameters(defaults_f);
    parse_int_parameters(defaults_i);
    parse_bool_parameters(defaults_b);
    parse_string_parameters(defaults_s);
    parse_double_parameters(defaults_d);
}

void lsdtt_param_parser::parse_float_parameters(
    const std::map<std::string, float>& defaults)
{
    for (const auto& [key, default_val] : defaults)
    {
        if (parameter_map_.count(key))
        {
            float_parameters_[key] = std::stof(parameter_map_[key]);
            parameters_read_map_[key] = parameter_map_[key];
        }
        else
        {
            float_parameters_[key] = default_val;
            defaults_used_map_[key] = std::to_string(default_val);
        }
    }
}

void lsdtt_param_parser::parse_int_parameters(
    const std::map<std::string, int>& defaults)
{
    for (const auto& [key, default_val] : defaults)
    {
        if (parameter_map_.count(key))
        {
            int_parameters_[key] = std::stoi(parameter_map_[key]);
            parameters_read_map_[key] = parameter_map_[key];
        }
        else
        {
            int_parameters_[key] = default_val;
            defaults_used_map_[key] = std::to_string(default_val);
        }
    }
}

void lsdtt_param_parser::parse_bool_parameters(
    const std::map<std::string, bool>& defaults)
{
    for (const auto& [key, default_val] : defaults)
    {
        if (parameter_map_.count(key))
        {
            const std::string& v = parameter_map_[key];
            bool result = (v == "true" || v == "True" || v == "TRUE" ||
                           v == "t"    || v == "T"    || v == "1");
            bool_parameters_[key] = result;
            parameters_read_map_[key] = v;
        }
        else
        {
            bool_parameters_[key] = default_val;
            defaults_used_map_[key] = default_val ? "true" : "false";
        }
    }
}

void lsdtt_param_parser::parse_string_parameters(
    const std::map<std::string, std::string>& defaults)
{
    for (const auto& [key, default_val] : defaults)
    {
        if (parameter_map_.count(key))
        {
            string_parameters_[key] = parameter_map_[key];
            parameters_read_map_[key] = parameter_map_[key];
        }
        else
        {
            string_parameters_[key] = default_val;
            defaults_used_map_[key] = default_val;
        }
    }
}

void lsdtt_param_parser::parse_double_parameters(
    const std::map<std::string, double>& defaults)
{
    for (const auto& [key, default_val] : defaults)
    {
        if (parameter_map_.count(key))
        {
            double_parameters_[key] = std::stod(parameter_map_[key]);
            parameters_read_map_[key] = parameter_map_[key];
        }
        else
        {
            double_parameters_[key] = default_val;
            defaults_used_map_[key] = std::to_string(default_val);
        }
    }
}

void lsdtt_param_parser::force_parse()
{
    // Reserved IO keys - skip these when force-parsing
    static const std::vector<std::string> io_keys = {
        "read_fname", "write_fname", "read_path", "write_path",
        "boundary_conditions", "NULL"
    };

    for (const auto& [key, value] : parameter_map_)
    {
        bool is_io = false;
        for (const auto& io : io_keys)
        {
            if (key == io) { is_io = true; break; }
        }
        if (!is_io)
        {
            parameters_read_map_[key] = value;
        }
    }
}

// -----------------------------------------------------------------------------
// Vector parameter parsing
// -----------------------------------------------------------------------------

std::vector<std::string> lsdtt_param_parser::parse_string_vector(
    const std::string& key) const
{
    std::vector<std::string> result;
    auto it = string_parameters_.find(key);
    if (it == string_parameters_.end()) return result;

    std::stringstream ss(it->second);
    std::string token;
    while (std::getline(ss, token, ','))
    {
        token.erase(std::remove_if(token.begin(), token.end(), ::isspace),
                    token.end());
        token.erase(std::remove_if(token.begin(), token.end(), ::iscntrl),
                    token.end());
        if (!token.empty()) result.push_back(token);
    }
    return result;
}

std::vector<int> lsdtt_param_parser::parse_int_vector(
    const std::string& key) const
{
    std::vector<int> result;
    for (const auto& s : parse_string_vector(key))
    {
        result.push_back(std::stoi(s));
    }
    return result;
}

std::vector<float> lsdtt_param_parser::parse_float_vector(
    const std::string& key) const
{
    std::vector<float> result;
    for (const auto& s : parse_string_vector(key))
    {
        result.push_back(std::stof(s));
    }
    return result;
}

std::vector<double> lsdtt_param_parser::parse_double_vector(
    const std::string& key) const
{
    std::vector<double> result;
    for (const auto& s : parse_string_vector(key))
    {
        result.push_back(std::stod(s));
    }
    return result;
}

// -----------------------------------------------------------------------------
// Boundary conditions
// -----------------------------------------------------------------------------

void lsdtt_param_parser::check_boundary_conditions()
{
    if (boundary_conditions_.size() != 4)
    {
        std::cout << "lsdtt_param_parser: boundary conditions not set or "
                  << "not exactly 4 - defaulting to no-flux (n) on all sides.\n";
        boundary_conditions_.assign(4, "n");
    }
}

// -----------------------------------------------------------------------------
// File I/O path helpers
// -----------------------------------------------------------------------------

void lsdtt_param_parser::check_paths()
{
    // Ensure read_path ends with /
    if (!read_path_.empty())
    {
        read_path_ = fix_path(read_path_);
    }

    // write_path defaults to read_path if not set
    if (write_path_.empty())
    {
        write_path_ = read_path_;
    }
    else
    {
        write_path_ = fix_path(write_path_);
    }

    // Split extension from read_fname if present.
    // read_fname may arrive as "my_dem.tif" - split into stem "my_dem"
    // and extension ".tif". If no extension is present read_ext_ stays empty
    // and GDAL will try to open the filename as-is (may still work for
    // extensionless files or if the full path was specified differently).
    if (!read_fname_.empty())
    {
        std::string ext = extension(read_fname_);
        if (!ext.empty())
        {
            read_ext_   = ext;
            read_fname_ = stem(read_fname_);
        }
        // If no extension found, read_ext_ remains empty -
        // get_full_read_path() will return path + fname without extension
    }
    else
    {
        // read_fname not set - default to parameter file stem
        read_fname_ = stem(param_fname_);
        read_ext_   = "";
    }

    // write_fname defaults to read_fname stem (never includes extension)
    if (write_fname_.empty())
    {
        write_fname_ = read_fname_;
    }

    std::cout << "lsdtt_param_parser: read  : "
              << read_path_ << read_fname_ << read_ext_ << "\n";
    std::cout << "lsdtt_param_parser: write prefix: "
              << write_path_ << write_fname_ << "\n";
}

std::string lsdtt_param_parser::add_read_path_if_needed(
    const std::string& this_string) const
{
    if (this_string == "NULL" || this_string == "null" || this_string == "Null")
    {
        return "NULL";
    }
    if (this_string.find('/') != std::string::npos ||
        this_string.find('\\') != std::string::npos)
    {
        return this_string;  // already has a path
    }
    return read_path_ + this_string;
}

std::string lsdtt_param_parser::stem(const std::string& filename)
{
    std::size_t pos = filename.find_last_of('.');
    if (pos == std::string::npos) return filename;
    return filename.substr(0, pos);
}

std::string lsdtt_param_parser::extension(const std::string& filename)
{
    std::size_t pos = filename.find_last_of('.');
    if (pos == std::string::npos) return "";
    return filename.substr(pos);  // includes the dot
}

// -----------------------------------------------------------------------------
// lsdtt_help_registry
// -----------------------------------------------------------------------------

std::vector<std::string>& lsdtt_help_registry::operator[](const std::string& key)
{
    auto it = entries_.find(key);
    if (it == entries_.end())
    {
        key_order_.push_back(key);
        return entries_[key];
    }
    return it->second;
}

std::string lsdtt_help_registry::field(
    const std::vector<std::string>& entry, std::size_t index)
{
    return index < entry.size() ? entry[index] : std::string();
}

// -----------------------------------------------------------------------------
// Deferred argument parsing and terminal help
// -----------------------------------------------------------------------------

namespace {

// Entry field indices (see lsdtt_help_registry documentation).
constexpr std::size_t k_help_type    = 0;
constexpr std::size_t k_help_default = 1;
constexpr std::size_t k_help_short   = 2;
constexpr std::size_t k_help_long    = 3;
constexpr std::size_t k_help_group   = 4;
constexpr std::size_t k_help_role    = 5;
constexpr std::size_t k_help_parent  = 6;

constexpr std::size_t k_help_wrap_width = 100;

// Wrap text at word boundaries; every output line starts with indent.
void write_wrapped(std::ostream& os, const std::string& text,
                   const std::string& indent)
{
    std::istringstream words(text);
    std::string word;
    std::string line = indent;
    bool line_has_word = false;
    while (words >> word)
    {
        if (line_has_word &&
            line.size() + 1 + word.size() > k_help_wrap_width)
        {
            os << line << "\n";
            line = indent;
            line_has_word = false;
        }
        if (line_has_word) line += " ";
        line += word;
        line_has_word = true;
    }
    if (line_has_word) os << line << "\n";
}

// One overview line: "<indent><key> (<type>=<default>)  <short_desc>",
// with the description wrapped onto its own indented line when too long.
void write_overview_line(std::ostream& os,
                         const std::string& indent,
                         const std::string& key,
                         const std::vector<std::string>& entry)
{
    const std::string type    = lsdtt_help_registry::field(entry, k_help_type);
    const std::string dflt    = lsdtt_help_registry::field(entry, k_help_default);
    const std::string short_d = lsdtt_help_registry::field(entry, k_help_short);

    std::string head = indent + key;
    if (!type.empty())
    {
        head += " (" + type;
        if (!dflt.empty()) head += "=" + dflt;
        head += ")";
    }

    if (!short_d.empty() &&
        head.size() + 2 + short_d.size() <= k_help_wrap_width)
    {
        os << head << "  " << short_d << "\n";
    }
    else
    {
        os << head << "\n";
        if (!short_d.empty()) write_wrapped(os, short_d, indent + "    ");
    }
}

std::string group_of(const std::vector<std::string>& entry)
{
    std::string g = lsdtt_help_registry::field(entry, k_help_group);
    return g.empty() ? std::string("General") : g;
}

std::string role_of(const std::vector<std::string>& entry)
{
    return lsdtt_help_registry::field(entry, k_help_role);
}

std::string parent_of(const std::vector<std::string>& entry)
{
    return lsdtt_help_registry::field(entry, k_help_parent);
}

// Return a copy of the registry with the parser-native file I/O keys
// (read_path, read_fname, write_path, write_fname) inserted first. These
// are handled by lsdtt_param_parser itself, so CLIs do not declare them;
// a CLI's own entry wins if it does.
lsdtt_help_registry with_builtin_io_entries(const lsdtt_help_registry& help)
{
    static const std::vector<
        std::pair<std::string, std::vector<std::string>>> builtin = {
        {"read_path", {"string", "./",
            "Directory from which input files are read.",
            "Defaults to the directory holding the parameter file.",
            "File I/O", "io", ""}},
        {"read_fname", {"string", "",
            "Input filename, with or without extension (e.g. my_dem.tif).",
            "Read from read_path. lsdtt3 uses GDAL, so any GDAL-readable "
            "raster format works. If empty, defaults to the parameter "
            "filename stem.",
            "File I/O", "io", ""}},
        {"write_path", {"string", "",
            "Directory to which output files are written.",
            "Defaults to read_path when not set.",
            "File I/O", "io", ""}},
        {"write_fname", {"string", "",
            "Output filename prefix (no path, no extension).",
            "All outputs are named <write_fname>_<descriptor>.<ext>. "
            "Defaults to the stem of read_fname when not set.",
            "File I/O", "io", ""}},
    };

    lsdtt_help_registry out;
    const auto& entries = help.entries();
    for (const auto& [key, entry] : builtin)
    {
        auto it = entries.find(key);
        out[key] = (it != entries.end()) ? it->second : entry;
    }
    for (const auto& key : help.key_order())
    {
        auto out_it = out.entries().find(key);
        if (out_it == out.entries().end())
        {
            out[key] = entries.at(key);
        }
    }
    return out;
}

// Detailed view for "-h <parameter_name>".
void write_topic_help(std::ostream& os,
                      const std::string& tool_name,
                      const lsdtt_help_registry& help,
                      const std::string& topic)
{
    const auto& entries = help.entries();
    auto it = entries.find(topic);
    if (it == entries.end())
    {
        os << "Unknown parameter: " << topic << "\n";
        std::vector<std::string> close;
        for (const auto& key : help.key_order())
        {
            if (key.find(topic) != std::string::npos) close.push_back(key);
        }
        if (!close.empty())
        {
            os << "Close matches:\n";
            for (const auto& key : close) os << "  " << key << "\n";
        }
        os << "Run \"" << tool_name << " -h\" for the full parameter list.\n";
        return;
    }

    const auto& e = it->second;
    os << topic << "\n";
    os << "  type    : " << lsdtt_help_registry::field(e, k_help_type) << "\n";
    os << "  default : " << lsdtt_help_registry::field(e, k_help_default) << "\n";
    os << "  group   : " << group_of(e) << "\n";
    const std::string role = role_of(e);
    if (!role.empty()) os << "  role    : " << role << "\n";
    const std::string parent = parent_of(e);
    if (!parent.empty()) os << "  used by : " << parent << "\n";
    os << "\n";
    write_wrapped(os, lsdtt_help_registry::field(e, k_help_short), "  ");
    const std::string long_d = lsdtt_help_registry::field(e, k_help_long);
    if (!long_d.empty())
    {
        os << "\n";
        write_wrapped(os, long_d, "  ");
    }

    if (role == "switch")
    {
        std::vector<std::string> children;
        for (const auto& key : help.key_order())
        {
            auto child = entries.find(key);
            if (child != entries.end() && parent_of(child->second) == topic)
            {
                children.push_back(key);
            }
        }
        if (!children.empty())
        {
            os << "\nParameters controlled by this switch:\n";
            for (const auto& key : children) os << "  " << key << "\n";
        }
    }
}

// Grouped overview for plain "-h".
void write_overview_help(std::ostream& os,
                         const std::string& tool_name,
                         const std::string& version,
                         const std::string& citation,
                         const lsdtt_help_registry& help)
{
    os << tool_name << " " << version << "\n";
    os << "Citation: " << citation << "\n\n";
    os << "Usage: " << tool_name << " <path> <parameter_file>\n";
    os << "       " << tool_name << " <parameter_file>\n";
    os << "Parameter file lines are \"key: value\"; # starts a comment.\n";

    const auto& entries = help.entries();
    const auto& order   = help.key_order();

    // File I/O and run control first.
    bool io_header_written = false;
    for (const auto& key : order)
    {
        const auto& e = entries.at(key);
        if (role_of(e) != "io") continue;
        if (!io_header_written)
        {
            os << "\nFile I/O and run control:\n";
            io_header_written = true;
        }
        write_overview_line(os, "  ", key, e);
    }

    // Group order = first appearance among non-io entries.
    std::vector<std::string> group_order;
    for (const auto& key : order)
    {
        const auto& e = entries.at(key);
        if (role_of(e) == "io") continue;
        const std::string g = group_of(e);
        if (std::find(group_order.begin(), group_order.end(), g) ==
            group_order.end())
        {
            group_order.push_back(g);
        }
    }

    for (const auto& group : group_order)
    {
        os << "\n" << group << ":\n";

        // Switches first, each followed by its dependent parameters.
        for (const auto& key : order)
        {
            const auto& e = entries.at(key);
            if (role_of(e) == "io" || group_of(e) != group) continue;
            if (role_of(e) != "switch") continue;
            write_overview_line(os, "  ", key, e);
            for (const auto& child_key : order)
            {
                const auto& c = entries.at(child_key);
                if (parent_of(c) == key && group_of(c) == group)
                {
                    write_overview_line(os, "    ", child_key, c);
                }
            }
        }

        // Remaining parameters in the group (no parent switch, not io).
        for (const auto& key : order)
        {
            const auto& e = entries.at(key);
            if (role_of(e) == "io" || group_of(e) != group) continue;
            if (role_of(e) == "switch" || !parent_of(e).empty()) continue;
            write_overview_line(os, "  ", key, e);
        }
    }

    os << "\nRun \"" << tool_name
       << " -h <parameter_name>\" for full guidance on one parameter.\n";
    os << "Run \"" << tool_name
       << " -h html\" to write an HTML/CSV reference to the current directory.\n";
    os << "Run \"" << tool_name
       << " --help-json\" for machine-readable output.\n";
}

// Escape minimal HTML special characters.
std::string html_escape(const std::string& s)
{
    std::string out;
    out.reserve(s.size());
    for (char c : s)
    {
        switch (c)
        {
            case '&':  out += "&amp;";  break;
            case '<':  out += "&lt;";   break;
            case '>':  out += "&gt;";   break;
            case '"':  out += "&quot;"; break;
            default:   out += c;        break;
        }
    }
    return out;
}

// Escape a string for a JSON string literal.
std::string json_escape(const std::string& s)
{
    std::string out;
    out.reserve(s.size());
    for (char c : s)
    {
        switch (c)
        {
            case '"':  out += "\\\""; break;
            case '\\': out += "\\\\"; break;
            case '\n': out += "\\n";  break;
            case '\r': out += "\\r";  break;
            case '\t': out += "\\t";  break;
            default:
                if (static_cast<unsigned char>(c) < 0x20)
                {
                    char buf[8];
                    std::snprintf(buf, sizeof(buf), "\\u%04x", c);
                    out += buf;
                }
                else
                {
                    out += c;
                }
                break;
        }
    }
    return out;
}

// Machine-readable help for "--help-json". Consumed by the QGIS plugin for
// capability gating.
void write_json_help(std::ostream& os,
                     const std::string& tool_name,
                     const std::string& version,
                     const std::string& citation,
                     const lsdtt_help_registry& help,
                     const std::vector<std::string>& features)
{
    os << "{\n";
    os << "  \"tool\": \"" << json_escape(tool_name) << "\",\n";
    os << "  \"version\": \"" << json_escape(version) << "\",\n";
    os << "  \"citation\": \"" << json_escape(citation) << "\",\n";

    os << "  \"features\": [";
    for (std::size_t i = 0; i < features.size(); ++i)
    {
        if (i > 0) os << ", ";
        os << "\"" << json_escape(features[i]) << "\"";
    }
    os << "],\n";

    os << "  \"parameters\": [\n";
    const auto& entries = help.entries();
    const auto& order   = help.key_order();
    for (std::size_t i = 0; i < order.size(); ++i)
    {
        const auto& key = order[i];
        const auto& e   = entries.at(key);
        os << "    {\"name\": \"" << json_escape(key) << "\", "
           << "\"type\": \""
           << json_escape(lsdtt_help_registry::field(e, k_help_type)) << "\", "
           << "\"default\": \""
           << json_escape(lsdtt_help_registry::field(e, k_help_default)) << "\", "
           << "\"description\": \""
           << json_escape(lsdtt_help_registry::field(e, k_help_short)) << "\", "
           << "\"guidance\": \""
           << json_escape(lsdtt_help_registry::field(e, k_help_long)) << "\", "
           << "\"group\": \""
           << json_escape(lsdtt_help_registry::field(e, k_help_group)) << "\", "
           << "\"role\": \""
           << json_escape(lsdtt_help_registry::field(e, k_help_role)) << "\", "
           << "\"parent\": \""
           << json_escape(lsdtt_help_registry::field(e, k_help_parent)) << "\"}";
        os << (i + 1 < order.size() ? ",\n" : "\n");
    }
    os << "  ]\n";
    os << "}\n";
}

// CSV-quote a field (always quoted so commas in prose are safe).
std::string csv_quote(const std::string& s)
{
    std::string out = "\"";
    for (char c : s)
    {
        if (c == '"') out += "\"\"";
        else          out += c;
    }
    out += "\"";
    return out;
}

} // anonymous namespace

lsdtt_cli_args lsdtt_param_parser::parse_args_deferred(
    int argc, char* argv[],
    const std::string& tool_name,
    const std::string& version,
    const std::string& citation)
{
    auto is_flag = [](const std::string& s,
                      std::initializer_list<const char*> flags)
    {
        for (const auto* f : flags) if (s == f) return true;
        return false;
    };

    lsdtt_cli_args args;

    for (int i = 1; i < argc; ++i)
    {
        std::string a = argv[i];

        if (is_flag(a, {"--help-json", "-help-json"}))
        {
            args.help_requested = true;
            args.help_json = true;
            return args;
        }
        if (is_flag(a, {"-h", "--help", "help", "-H", "--H"}))
        {
            args.help_requested = true;
            if (i + 1 < argc)
            {
                const std::string next = argv[i + 1];
                if (next == "html")
                {
                    args.help_html = true;
                }
                else if (next == "json")
                {
                    args.help_json = true;
                }
                else if (!next.empty() && next[0] != '-')
                {
                    args.help_topic = next;
                }
            }
            return args;
        }
        if (is_flag(a, {"-v", "--version", "version", "-V"}))
        {
            std::cout << tool_name << " version " << version
                      << " (lsdtt3 " << LSDTT3_VERSION << ")\n";
            std::exit(0);
        }
        if (is_flag(a, {"-c", "--cite", "citation", "-citation", "--citation"}))
        {
            std::cout << "Citation: " << citation << "\n";
            std::exit(0);
        }
    }

    if (argc == 1)
    {
        std::cerr << "No parameter file specified. Use -h for help.\n";
        std::exit(1);
    }
    else if (argc == 2)
    {
        args.fname = argv[1];
    }
    else
    {
        args.path  = argv[1];
        args.fname = argv[2];
    }

    return args;
}

void lsdtt_param_parser::print_terminal_help(
    const std::string& tool_name,
    const std::string& version,
    const std::string& citation,
    const lsdtt_help_registry& help,
    const lsdtt_cli_args& args,
    std::ostream& os,
    const std::vector<std::string>& features)
{
    const lsdtt_help_registry full = with_builtin_io_entries(help);

    if (args.help_json)
    {
        write_json_help(os, tool_name, version, citation, full, features);
        return;
    }

    if (args.help_html)
    {
        write_reference_files(tool_name, version, citation, full,
                              "./" + tool_name + "_help");
        os << "Reference files written to:\n"
           << "  ./" << tool_name << "_help.html\n"
           << "  ./" << tool_name << "_help.csv\n";
        return;
    }

    if (!args.help_topic.empty())
    {
        write_topic_help(os, tool_name, full, args.help_topic);
        return;
    }

    write_overview_help(os, tool_name, version, citation, full);
}

void lsdtt_param_parser::print_terminal_help(
    const std::string& tool_name,
    const std::string& version,
    const std::string& citation,
    const lsdtt_help_registry& help,
    const lsdtt_cli_args& args,
    const std::vector<std::string>& features)
{
    print_terminal_help(tool_name, version, citation, help, args,
                        std::cout, features);
}

void lsdtt_param_parser::write_reference_files(
    const std::string& tool_name,
    const std::string& version,
    const std::string& citation,
    const lsdtt_help_registry& help,
    const std::string& out_prefix)
{
    const lsdtt_help_registry full = with_builtin_io_entries(help);
    const auto& entries = full.entries();
    const auto& order   = full.key_order();

    // ---- CSV (flat, grouped columns appended) ----
    {
        std::ofstream csv(out_prefix + ".csv");
        if (csv.is_open())
        {
            csv << "name,type,default,description,guidance,group,role,parent\n";
            csv << "version," << csv_quote(version) << ",citation,"
                << csv_quote(citation) << ",,,,\n";
            for (const auto& key : order)
            {
                const auto& e = entries.at(key);
                csv << key << ","
                    << csv_quote(lsdtt_help_registry::field(e, k_help_type))    << ","
                    << csv_quote(lsdtt_help_registry::field(e, k_help_default)) << ","
                    << csv_quote(lsdtt_help_registry::field(e, k_help_short))   << ","
                    << csv_quote(lsdtt_help_registry::field(e, k_help_long))    << ","
                    << csv_quote(lsdtt_help_registry::field(e, k_help_group))   << ","
                    << csv_quote(lsdtt_help_registry::field(e, k_help_role))    << ","
                    << csv_quote(lsdtt_help_registry::field(e, k_help_parent))  << "\n";
            }
        }
    }

    // ---- HTML (grouped layout) ----
    std::ofstream html(out_prefix + ".html");
    if (!html.is_open()) return;

    html << "<!DOCTYPE html>\n<html>\n<head>\n"
         << "<title>Help: " << html_escape(tool_name) << "</title>\n"
         << "<style>table,th,td{border:1px solid black;"
         << "border-collapse:collapse;padding:4px;}"
         << "h2{margin-top:1.5em;}</style>\n"
         << "</head>\n<body>\n";

    html << "<h1>Help: " << html_escape(tool_name) << "</h1>\n";
    html << "<p>Version: " << html_escape(version) << "</p>\n";
    html << "<p>Citation: " << html_escape(citation) << "</p>\n";
    html << "<h2>Usage</h2>\n"
         << "<p><code>" << html_escape(tool_name)
         << " my_params.driver</code></p>\n"
         << "<p>Parameter file format: <code>key: value</code> "
         << "(one per line, # for comments)</p>\n";

    auto write_row = [&html](const std::string& key,
                             const std::vector<std::string>& e,
                             bool indent = false)
    {
        html << "<tr><td>" << (indent ? "&nbsp;&nbsp;" : "")
             << html_escape(key) << "</td><td>"
             << html_escape(lsdtt_help_registry::field(e, k_help_type))
             << "</td><td>"
             << html_escape(lsdtt_help_registry::field(e, k_help_default))
             << "</td><td>"
             << html_escape(lsdtt_help_registry::field(e, k_help_short))
             << "</td><td>"
             << html_escape(lsdtt_help_registry::field(e, k_help_long))
             << "</td></tr>\n";
    };

    const char* table_header =
        "<table>\n<tr><th>name</th><th>type</th><th>default</th>"
        "<th>description</th><th>guidance</th></tr>\n";

    // File I/O section.
    bool io_open = false;
    for (const auto& key : order)
    {
        const auto& e = entries.at(key);
        if (role_of(e) != "io") continue;
        if (!io_open)
        {
            html << "<h2>File I/O and run control</h2>\n" << table_header;
            io_open = true;
        }
        write_row(key, e);
    }
    if (io_open) html << "</table>\n";

    // Groups in declaration order.
    std::vector<std::string> group_order;
    for (const auto& key : order)
    {
        const auto& e = entries.at(key);
        if (role_of(e) == "io") continue;
        const std::string g = group_of(e);
        if (std::find(group_order.begin(), group_order.end(), g) ==
            group_order.end())
        {
            group_order.push_back(g);
        }
    }

    for (const auto& group : group_order)
    {
        html << "<h2>" << html_escape(group) << "</h2>\n" << table_header;
        // Switches with their children directly beneath, then loose params.
        for (const auto& key : order)
        {
            const auto& e = entries.at(key);
            if (role_of(e) == "io" || group_of(e) != group) continue;
            if (role_of(e) != "switch") continue;
            write_row(key, e);
            for (const auto& child_key : order)
            {
                const auto& c = entries.at(child_key);
                if (parent_of(c) == key && group_of(c) == group)
                {
                    write_row(child_key, c, true);
                }
            }
        }
        for (const auto& key : order)
        {
            const auto& e = entries.at(key);
            if (role_of(e) == "io" || group_of(e) != group) continue;
            if (role_of(e) == "switch" || !parent_of(e).empty()) continue;
            write_row(key, e);
        }
        html << "</table>\n";
    }

    html << "</body>\n</html>\n";
}

// -----------------------------------------------------------------------------
// Output / diagnostics
// -----------------------------------------------------------------------------

void lsdtt_param_parser::print_parameters() const
{
    print_parameters(write_fname_ + "_ingested_params.param");
}

void lsdtt_param_parser::print_parameters(const std::string& fname_prefix) const
{
    std::string fname = write_path_ + fname_prefix;
    std::ofstream out(fname);
    if (!out.is_open())
    {
        std::cerr << "lsdtt_param_parser::print_parameters: "
                  << "could not open " << fname << "\n";
        return;
    }

    out << "# lsdtt3 parameter file - ingested parameters\n";
    out << "# File I/O\n";
    out << "read_path: "   << read_path_            << "\n";
    out << "read_fname: "  << read_fname_ << read_ext_ << "\n";
    out << "write_path: "  << write_path_  << "\n";
    out << "write_fname: " << write_fname_ << "\n\n";

    out << "# Parameters read from file\n";
    for (const auto& [key, val] : parameters_read_map_)
    {
        out << key << ": " << val << "\n";
    }

    out << "\n# Default parameters used\n";
    for (const auto& [key, val] : defaults_used_map_)
    {
        out << key << ": " << val << "\n";
    }
}

void lsdtt_param_parser::print_help(
    const std::map<std::string, std::vector<std::string>>& help_map,
    const std::string& file_prefix,
    const std::string& version,
    const std::string& citation
) const
{
    // Collect all keys from read and defaults maps
    std::vector<std::string> all_keys;
    for (const auto& [k, v] : parameters_read_map_) all_keys.push_back(k);
    for (const auto& [k, v] : defaults_used_map_)   all_keys.push_back(k);
    std::sort(all_keys.begin(), all_keys.end());
    all_keys.erase(std::unique(all_keys.begin(), all_keys.end()), all_keys.end());

    // CSV help file
    std::string csv_path = write_path_ + file_prefix + ".csv";
    std::ofstream csv(csv_path);
    if (csv.is_open())
    {
        // group/role/parent columns are appended so readers of the previous
        // five-column format keep working; entries without the grouped
        // metadata leave those columns empty.
        csv << "name,type,default,description,guidance,group,role,parent\n";
        csv << "version," << version << ",citation," << citation
            << ",see code output for additional algorithm citations,,,\n";

        for (const auto& key : all_keys)
        {
            auto it = help_map.find(key);
            if (it == help_map.end() || it->second.size() < 4)
            {
                csv << key << ",not available,not available,"
                    << "not available,not available,,,\n";
            }
            else
            {
                const auto& h = it->second;
                csv << key << "," << h[0] << "," << h[1] << ","
                    << h[2] << "," << h[3] << ","
                    << lsdtt_help_registry::field(h, 4) << ","
                    << lsdtt_help_registry::field(h, 5) << ","
                    << lsdtt_help_registry::field(h, 6) << "\n";
            }
        }
        csv.close();
    }

    // HTML help file
    std::string html_path = write_path_ + file_prefix + ".html";
    std::ofstream html(html_path);
    if (html.is_open())
    {
        std::string prog = stem(file_prefix);

        html << "<!DOCTYPE html>\n<html>\n<head>\n"
             << "<title>Help: " << prog << "</title>\n"
             << "<style>table,th,td{border:1px solid black;"
             << "border-collapse:collapse;padding:4px;}</style>\n"
             << "</head>\n<body>\n";

        html << "<h1>Help: " << prog << "</h1>\n";
        html << "<p>Version: " << version << "</p>\n";
        html << "<p>Citation: <a href=\"" << citation << "\">"
             << citation << "</a></p>\n";

        html << "<h2>Usage</h2>\n"
             << "<p>Call with a parameter file: "
             << "<code>" << prog << " my_params.driver</code></p>\n"
             << "<p>Parameter file format: "
             << "<code>key: value</code> (one per line, # for comments)</p>\n"
             << "<p>All keys use snake_case. No spaces in key names.</p>\n";

        html << "<h2>Parameters</h2>\n"
             << "<table>\n"
             << "<tr><th>name</th><th>type</th><th>default</th>"
             << "<th>description</th><th>guidance</th>"
             << "<th>group</th><th>role</th><th>parent</th></tr>\n";

        for (const auto& key : all_keys)
        {
            auto it = help_map.find(key);
            if (it == help_map.end() || it->second.size() < 4)
            {
                html << "<tr><td>" << key
                     << "</td><td colspan=\"7\">not available</td></tr>\n";
            }
            else
            {
                const auto& h = it->second;
                html << "<tr><td>" << key << "</td><td>" << h[0]
                     << "</td><td>" << h[1] << "</td><td>" << h[2]
                     << "</td><td>" << h[3] << "</td><td>"
                     << lsdtt_help_registry::field(h, 4) << "</td><td>"
                     << lsdtt_help_registry::field(h, 5) << "</td><td>"
                     << lsdtt_help_registry::field(h, 6) << "</td></tr>\n";
            }
        }
        html << "</table>\n</body>\n</html>\n";
        html.close();
    }

    std::cout << "lsdtt_param_parser: help files written to:\n"
              << "  " << csv_path  << "\n"
              << "  " << html_path << "\n";
}

void lsdtt_param_parser::print_help(
    const lsdtt_help_registry& help,
    const std::string& file_prefix,
    const std::string& version,
    const std::string& citation
) const
{
    print_help(help.entries(), file_prefix, version, citation);
}

void lsdtt_param_parser::replace_and_print_parameter_file(
    const std::string& output_fname,
    const std::string& new_read_path,
    const std::string& new_read_fname,
    const std::string& new_write_path,
    const std::string& new_write_fname,
    const std::map<std::string, std::string>& replacements
) const
{
    std::ofstream out(output_fname);
    if (!out.is_open())
    {
        throw std::runtime_error(
            "lsdtt_param_parser::replace_and_print_parameter_file: "
            "cannot open " + output_fname
        );
    }

    out << "# lsdtt3 parameter file (generated)\n";
    out << "read_path: "   << fix_path(new_read_path)  << "\n";
    out << "read_fname: "  << new_read_fname            << "\n";
    out << "write_path: "  << fix_path(new_write_path) << "\n";
    out << "write_fname: " << new_write_fname           << "\n\n";

    // Write all parameters, applying replacements
    for (const auto& [key, val] : parameter_map_)
    {
        // Skip IO keys - already written above
        if (key == "read_path"  || key == "write_path" ||
            key == "read_fname" || key == "write_fname") continue;

        auto rep_it = replacements.find(key);
        if (rep_it != replacements.end())
        {
            out << key << ": " << rep_it->second << "\n";
        }
        else
        {
            out << key << ": " << val << "\n";
        }
    }
}

// -----------------------------------------------------------------------------
// Private helpers
// -----------------------------------------------------------------------------

void lsdtt_param_parser::parse_line(
    std::ifstream& infile,
    std::string&   parameter,
    std::string&   value) const
{
    char    c;
    char    buff[1024];
    int     pos  = 0;
    int     word = 0;

    parameter = "NULL";
    value     = "NULL";

    while (infile.get(c))
    {
        if (pos >= 1023)
        {
            std::cerr << "lsdtt_param_parser: buffer overrun on line starting: "
                      << buff << "\n";
            std::string line;
            std::getline(infile, line);
            return;
        }

        // Comment character - stop reading this line
        if (c == '#')
        {
            if (word == 1) value = "NULL";
            word = 2;
        }

        if ((c == ' ' || c == '\t') && pos == 0)
        {
            continue;  // leading whitespace
        }
        else if ((c == ':' && word == 0) ||
                 ((c == ' ' || c == '\n' || c == '\t') && word == 1))
        {
            // Strip trailing whitespace
            while (pos > 0 && (buff[pos-1] == ' ' || buff[pos-1] == '\t'))
            {
                --pos;
            }
            buff[pos] = '\0';
            if (word == 0)      parameter = buff;
            else if (word == 1) value     = buff;
            ++word;
            pos = 0;
        }
        else if (c == '\n' && word == 0)
        {
            parameter = "NULL";
            buff[pos] = '\0';
            value = buff;
            ++word;
        }
        else if (word < 2)
        {
            buff[pos++] = c;
        }

        if (c == '\n') break;
    }

    if (word == 0)
    {
        parameter = "NULL";
        value     = "NULL";
    }
}

std::string lsdtt_param_parser::strip_control_chars(const std::string& s)
{
    std::string result = s;
    result.erase(
        std::remove_if(result.begin(), result.end(),
                       [](unsigned char c)
                       {
                           return c == '\r' || c == '\n' || c == '\t' ||
                                  (c < 32 && c != ' ');
                       }),
        result.end()
    );
    return result;
}

std::string lsdtt_param_parser::fix_path(const std::string& path)
{
    if (path.empty()) return "./";
    if (path.back() != '/' && path.back() != '\\')
    {
        return path + "/";
    }
    return path;
}

template<typename V>
std::vector<std::string> lsdtt_param_parser::keys_of(
    const std::map<std::string, V>& m)
{
    std::vector<std::string> keys;
    keys.reserve(m.size());
    for (const auto& [k, v] : m)
    {
        keys.push_back(k);
    }
    return keys;
}

} // namespace lsdtt
