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
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>

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
    // Reserved IO keys — skip these when force-parsing
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
                  << "not exactly 4 — defaulting to no-flux (n) on all sides.\n";
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
    // read_fname may arrive as "my_dem.tif" — split into stem "my_dem"
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
        // If no extension found, read_ext_ remains empty —
        // get_full_read_path() will return path + fname without extension
    }
    else
    {
        // read_fname not set — default to parameter file stem
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

std::pair<std::string, std::string> lsdtt_param_parser::parse_args(
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

    for (int i = 1; i < argc; ++i)
    {
        std::string a = argv[i];

        if (is_flag(a, {"-h", "--help", "help", "-H", "--H"}))
        {
            std::cout
                << "Usage: " << tool_name << " <path> <parameter_file>\n"
                << "       " << tool_name << " <parameter_file>\n"
                << "       " << tool_name << " -h   (this help)\n"
                << "       " << tool_name << " -v   (version)\n"
                << "       " << tool_name << " -c   (citation)\n";
            std::exit(0);
        }
        if (is_flag(a, {"-v", "--version", "version", "-V"}))
        {
            std::cout << tool_name << " version " << version << "\n";
            std::exit(0);
        }
        if (is_flag(a, {"-c", "--cite", "citation", "-citation", "--citation"}))
        {
            std::cout << "Citation: " << citation << "\n";
            std::exit(0);
        }
    }

    std::string path = "./";
    std::string fname;

    if (argc == 1)
    {
        std::cerr << "No parameter file specified. Use -h for help.\n";
        std::exit(1);
    }
    else if (argc == 2)
    {
        fname = argv[1];
    }
    else
    {
        path  = argv[1];
        fname = argv[2];
    }

    return {path, fname};
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

    out << "# lsdtt3 parameter file — ingested parameters\n";
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
        csv << "name,type,default,description,guidance\n";
        csv << "version," << version << ",citation," << citation
            << ",see code output for additional algorithm citations\n";

        for (const auto& key : all_keys)
        {
            auto it = help_map.find(key);
            if (it == help_map.end() || it->second.size() < 4)
            {
                csv << key << ",not available,not available,"
                    << "not available,not available\n";
            }
            else
            {
                const auto& h = it->second;
                csv << key << "," << h[0] << "," << h[1] << ","
                    << h[2] << "," << h[3] << "\n";
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
             << "<th>description</th><th>guidance</th></tr>\n";

        for (const auto& key : all_keys)
        {
            auto it = help_map.find(key);
            if (it == help_map.end() || it->second.size() < 4)
            {
                html << "<tr><td>" << key
                     << "</td><td colspan=\"4\">not available</td></tr>\n";
            }
            else
            {
                const auto& h = it->second;
                html << "<tr><td>" << key << "</td><td>" << h[0]
                     << "</td><td>" << h[1] << "</td><td>" << h[2]
                     << "</td><td>" << h[3] << "</td></tr>\n";
            }
        }
        html << "</table>\n</body>\n</html>\n";
        html.close();
    }

    std::cout << "lsdtt_param_parser: help files written to:\n"
              << "  " << csv_path  << "\n"
              << "  " << html_path << "\n";
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
        // Skip IO keys — already written above
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

        // Comment character — stop reading this line
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
