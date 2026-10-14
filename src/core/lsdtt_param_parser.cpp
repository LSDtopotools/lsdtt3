// lsdtt_param_parser.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// Original LSDTopoTools implementation by Simon M. Mudd and contributors.
// Redesigned for lsdtt3: removed extension logic, removed hardcoded file
// parameters, renamed keys to snake_case, removed external dependencies.
// This reads parameter files (which we sometimes call driver files)
// format is
// parameter_name: value
// lines with hash at beginning are treated as comments
#include "lsdtt_param_parser.hpp"
#include "../utils/lsdtt_output_file.hpp"

#include "lsdtt_gdal_init.hpp"

#include <algorithm>
#include <cmath>
#include <cctype>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <locale>
#include <set>
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
            "lsdtt_param_parser: parameter file not found: " + full_path +
            "\n  The parameter file is a plain-text file of \"key: value\" "
            "lines that you create.\n  Run the tool with no arguments for "
            "its format, or with -h to list its parameters."
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

// -----------------------------------------------------------------------------
// apply_key_policy
//
// Runs against parameter_map_ - the RAW KEYS READ FROM THE FILE - not against
// the defaults maps. See the header for why that distinction is the design.
//
// Order matters. Removed keys are checked first and stop the run, because there
// is no point warning about aliases in a file that is going to be rejected
// anyway. Transformed aliases and then ordinary aliases are rewritten into
// parameter_map_ so the typed parses that follow pick values up under their new
// names. Consumed old keys are erased so they do not also surface as unknown.
// -----------------------------------------------------------------------------
void lsdtt_param_parser::set_deprecated_compatibility_aliases(
    std::map<std::string, lsdtt_deprecated_compatibility_alias> aliases)
{
    deprecated_compatibility_alias_map_ = std::move(aliases);
    compatibility_storage_keys_.clear();
    for (const auto& [old_key, policy] : deprecated_compatibility_alias_map_)
    {
        (void)old_key;
        if (!policy.storage_key.empty())
            compatibility_storage_keys_.insert(policy.storage_key);
    }
}

void lsdtt_param_parser::apply_key_policy(
    const std::vector<std::string>& declared_keys)
{
    // Keys the parser consumes itself in parse_file_io(). They never appear in
    // a CLI's defaults maps, so without this they would every one of them be
    // reported as unknown.
    static const std::set<std::string> reserved = {
        "read_path", "write_path", "read_fname", "write_fname",
        "boundary_conditions"
    };

    // --- 1. Removed keys: print every one, then stop -------------------------
    std::vector<std::string> removed_present;
    for (const auto& [old_key, message] : removed_key_message_map_)
    {
        if (parameter_map_.count(old_key)) removed_present.push_back(old_key);
    }
    if (!removed_present.empty())
    {
        // All of them, not just the first: a user migrating an old file should
        // get the whole list in one run rather than one key per attempt.
        std::cerr << "\n";
        for (const std::string& key : removed_present)
        {
            std::cerr << "FATAL: " << removed_key_message_map_[key] << "\n";
        }
        std::cerr << "  " << removed_present.size()
                  << " removed parameter(s) in " << param_fname_
                  << ". Remove or migrate them and run again.\n";
        std::exit(1);
    }

    // A key cannot have two incompatible alias policies. Check the maps, not
    // just keys present in this file, so configuration errors are deterministic.
    for (const auto& [old_key, policy] : deprecated_numeric_alias_map_)
    {
        (void)policy;
        if (!deprecated_alias_map_.count(old_key)) continue;
        std::cerr << "FATAL: deprecated key '" << old_key
                  << "' is registered as both an ordinary and transformed "
                     "numeric alias.\n";
        std::exit(1);
    }

    for (const auto& [old_key, policy] : deprecated_compatibility_alias_map_)
    {
        if (deprecated_alias_map_.count(old_key) ||
            deprecated_numeric_alias_map_.count(old_key) ||
            policy.storage_key.empty() || policy.replacement_key.empty())
        {
            std::cerr << "FATAL: invalid or duplicate compatibility alias policy for '"
                      << old_key << "'.\n";
            std::exit(1);
        }
    }

    // --- 2. Transformed numeric aliases --------------------------------------
    for (const auto& [old_key, policy] : deprecated_numeric_alias_map_)
    {
        if (!parameter_map_.count(old_key)) continue;

        auto format_number = [](double value)
        {
            std::ostringstream out;
            out.imbue(std::locale::classic());
            out << std::setprecision(std::numeric_limits<double>::max_digits10)
                << value;
            return out.str();
        };
        const std::string transform =
            "new = old * " + format_number(policy.scale) + " + " +
            format_number(policy.offset) +
            (policy.preserve_zero ? "; zero remains zero" : "");

        if (policy.new_key.empty() || !std::isfinite(policy.scale) ||
            !std::isfinite(policy.offset))
        {
            std::cerr << "FATAL: invalid transformed numeric alias policy for '"
                      << old_key << "'.\n";
            std::exit(1);
        }

        const std::string& text = parameter_map_[old_key];
        double old_value = 0.0;
        std::size_t consumed = 0;
        try
        {
            old_value = std::stod(text, &consumed);
        }
        catch (const std::exception&)
        {
            consumed = 0;
        }
        if (consumed != text.size() || !std::isfinite(old_value))
        {
            std::cerr << "FATAL: deprecated numeric parameter '" << old_key
                      << "' must be a finite number for transformation to '"
                      << policy.new_key << "' (" << transform << "); got '"
                      << text << "'.\n";
            std::exit(1);
        }

        const bool canonical_supplied = parameter_map_.count(policy.new_key) != 0;
        if (canonical_supplied)
        {
            std::cerr << "WARNING: '" << old_key << "' is deprecated and '"
                      << policy.new_key << "' was also supplied. Using '"
                      << policy.new_key << "'; the deprecated key is ignored "
                      << "and its numeric transformation (" << transform
                      << ") was not applied.\n";
        }
        else
        {
            const double new_value =
                (policy.preserve_zero && old_value == 0.0)
                    ? 0.0
                    : old_value * policy.scale + policy.offset;
            if (!std::isfinite(new_value))
            {
                std::cerr << "FATAL: transformed value for deprecated numeric "
                             "parameter '" << old_key << "' is not finite ("
                          << transform << ").\n";
                std::exit(1);
            }
            parameter_map_[policy.new_key] = format_number(new_value);
            std::cerr << "WARNING: '" << old_key << "' is deprecated; use '"
                      << policy.new_key << "'. Its numeric value has been "
                      << "transformed (" << transform << ") and applied to '"
                      << policy.new_key << "' for this run.\n";
        }
        parameter_map_.erase(old_key);
    }

    // Compatibility aliases are injective storage rewrites. The hidden storage
    // key preserves an old per-output value without advertising a non-genuine
    // public override.
    for (const auto& [old_key, policy] : deprecated_compatibility_alias_map_)
    {
        if (!parameter_map_.count(old_key)) continue;
        if (parameter_map_.count(policy.storage_key))
        {
            std::cerr << "WARNING: '" << old_key << "' is deprecated and its "
                      << "compatibility value was already supplied. Using the "
                      << "existing value; migrate to '" << policy.replacement_key
                      << "'.\n";
        }
        else
        {
            parameter_map_[policy.storage_key] = parameter_map_[old_key];
            std::cerr << "WARNING: '" << old_key << "' is deprecated; use '"
                      << policy.replacement_key << "'. Its value is preserved "
                      << "for this output in this run; review consolidated "
                      << "format settings before removing the old key.\n";
        }
        parameter_map_.erase(old_key);
    }

    // --- 3. Deprecated aliases: carry the value across, name the replacement --
    for (const auto& [old_key, new_key] : deprecated_alias_map_)
    {
        if (!parameter_map_.count(old_key)) continue;

        if (parameter_map_.count(new_key))
        {
            std::cerr << "WARNING: '" << old_key << "' is deprecated and '"
                      << new_key << "' was also supplied. Using '" << new_key
                      << "'; the deprecated key is ignored.\n";
        }
        else
        {
            parameter_map_[new_key] = parameter_map_[old_key];
            std::cerr << "WARNING: '" << old_key << "' is deprecated; use '"
                      << new_key << "'. Its value has been applied to '"
                      << new_key << "' for this run.\n";
        }
        parameter_map_.erase(old_key);
    }

    // --- 4. Ignored keys: named, explained, and survived ---------------------
    // Erased after reporting so they do not also surface as unknown: the
    // specific message is strictly better than "not recognised", and printing
    // both for one key would read as two separate problems.
    for (const auto& [old_key, message] : ignored_key_message_map_)
    {
        if (!parameter_map_.count(old_key)) continue;
        std::cerr << "WARNING: " << message << "\n";
        parameter_map_.erase(old_key);
    }

    // --- 5. Unknown keys: a warning, never fatal -----------------------------
    if (!report_unknown_keys_) return;

    const std::set<std::string> declared(declared_keys.begin(),
                                         declared_keys.end());
    std::vector<std::string> unknown;
    for (const auto& [key, value] : parameter_map_)
    {
        (void)value;
        if (reserved.count(key) || declared.count(key)) continue;
        unknown.push_back(key);
    }
    if (unknown.empty()) return;

    std::cerr << "WARNING: " << unknown.size()
              << " parameter(s) in " << param_fname_
              << " are not used by this tool and were ignored:\n";
    for (const std::string& key : unknown)
        std::cerr << "  " << key << "\n";
    std::cerr << "  Check for a typo, or run with -h to list the parameters "
                 "this tool accepts.\n";
}

void lsdtt_param_parser::parse_all_parameters(
    const std::map<std::string, float>&       defaults_f,
    const std::map<std::string, int>&         defaults_i,
    const std::map<std::string, bool>&        defaults_b,
    const std::map<std::string, std::string>& defaults_s
)
{
    std::vector<std::string> declared;
    for (const auto& kv : defaults_f) declared.push_back(kv.first);
    for (const auto& kv : defaults_i) declared.push_back(kv.first);
    for (const auto& kv : defaults_b) declared.push_back(kv.first);
    for (const auto& kv : defaults_s) declared.push_back(kv.first);
    apply_key_policy(declared);

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
    std::vector<std::string> declared;
    for (const auto& kv : defaults_f) declared.push_back(kv.first);
    for (const auto& kv : defaults_i) declared.push_back(kv.first);
    for (const auto& kv : defaults_b) declared.push_back(kv.first);
    for (const auto& kv : defaults_s) declared.push_back(kv.first);
    for (const auto& kv : defaults_d) declared.push_back(kv.first);
    apply_key_policy(declared);

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
        if (!is_io && !compatibility_storage_keys_.count(key))
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

// How to run a tool and write its parameter file. Shared by plain "-h" and by
// a run with no arguments, so that someone holding only the binary (often an
// LLM with no access to the documentation) learns that options live in a
// parameter file rather than in command-line flags, and how that file is
// written. `help` is null on the no-argument path, which runs before the CLI
// has declared its parameters; when present it supplies a real example key.
void write_parameter_file_guide(std::ostream& os,
                                const std::string& tool_name,
                                const lsdtt_help_registry* help)
{
    os << "How to run " << tool_name << ":\n"
       << "  " << tool_name << " -h                        "
          "list every parameter this tool accepts\n"
       << "  " << tool_name << " -h <parameter_name>       "
          "full guidance on one parameter\n"
       << "  " << tool_name << " <parameter_file>          "
          "run the analysis described in the file\n"
       << "  " << tool_name << " <path> <parameter_file>   "
          "same, with the parameter file in directory <path>\n"
       << "\n"
       << "All analysis options are set in a plain-text parameter file, not with\n"
       << "command-line flags. The parameter names, types and defaults are\n"
       << "listed by \"" << tool_name << " -h\".\n"
       << "\n"
       << "Parameter file format (any filename; my_params.driver is the convention):\n"
       << "  - one \"key: value\" per line; the key is a parameter name from -h\n"
       << "  - keys are case-sensitive; any parameter left out takes its default\n"
       << "  - a value must not contain spaces: lists are comma-separated with no\n"
       << "    spaces (a,b,c), and paths must not contain spaces\n"
       << "  - booleans are true or false\n"
       << "  - # starts a comment; blank lines are ignored\n"
       << "  - read_path defaults to the parameter file's directory, and\n"
       << "    write_path defaults to read_path; write_path must already exist\n"
       << "  - outputs are named <write_fname>_<descriptor>.<ext> in write_path\n"
       << "  - unrecognised keys are reported as warnings and ignored\n"
       << "\n"
       << "Example parameter file:\n"
       << "  # my_params.driver\n"
       << "  read_path: /data/my_project\n"
       << "  write_path: /data/my_project/output\n"
       << "  read_fname: my_dem.tif\n"
       << "  write_fname: my_run\n";

    // One real parameter of this tool, so the example is runnable in shape.
    // Prefer a switch (they turn analyses on); fall back to any parameter with
    // a usable default.
    std::string example_key, example_value;
    if (help)
    {
        const auto& entries = help->entries();
        for (const std::string pass : {"switch", ""})
        {
            for (const auto& key : help->key_order())
            {
                const auto& e = entries.at(key);
                if (role_of(e) == "io") continue;
                if (!pass.empty() && role_of(e) != pass) continue;
                std::string v = lsdtt_help_registry::field(e, k_help_default);
                if (lsdtt_help_registry::field(e, k_help_type) == "bool")
                    v = "true";
                if (v.empty() || v.find(' ') != std::string::npos) continue;
                example_key = key;
                example_value = v;
                break;
            }
            if (!example_key.empty()) break;
        }
    }
    if (!example_key.empty())
    {
        os << "  " << example_key << ": " << example_value << "\n";
    }
    else
    {
        os << "  <parameter_name>: <value>   # parameters listed by -h\n";
    }
    os << "\n"
       << "Run it with:  " << tool_name << " /data/my_project my_params.driver\n";
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
    os << "  set in the parameter file as \"" << topic << ": <value>\"\n";
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
    write_parameter_file_guide(os, tool_name, &help);
    os << "\nParameters accepted in the parameter file "
          "(name (type=default)  description):\n";

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

    // EVERY declared parameter must appear exactly once.
    //
    // A `parent` is a hint about presentation: it lets a dependent parameter be
    // shown indented under the switch that enables it. It must never be able to
    // delete a parameter from the documentation, and it must never be able to
    // show one twice. Both used to happen. A parented entry was rendered only
    // when its parent was a switch in the same group; every other parented
    // entry was skipped by the nesting pass and then excluded from the flat
    // pass because it had a parent, so it appeared nowhere at all. That hid 68
    // declared parameters across 8 CLIs - including every ladder_gap_* control
    // in lsdtt-fluvial-surfaces - from -h, -h html and the per-run README
    // files, while the CLI reference pages told users those outputs were
    // complete.
    //
    // The flat pass therefore skips exactly what the nesting pass WROTE,
    // recorded as it goes, rather than re-deriving that from the parent's role
    // and group. Re-deriving it is what drifted: the nesting pass descends
    // through a chain of any roles, so a predicate that only looked one level
    // up and demanded role == "switch" disagreed with it and printed the
    // middle of every chain twice. `written` cannot disagree with itself.
    std::set<std::string> written;

    // Seed with the io entries: they were listed above under "File I/O and run
    // control", and a tree would otherwise print an io child a second time
    // inside its group.
    for (const auto& key : order)
    {
        if (role_of(entries.at(key)) == "io") written.insert(key);
    }

    for (const auto& group : group_order)
    {
        os << "\n" << group << ":\n";

        // Print a switch and, beneath it, everything that depends on it,
        // recursively: a switch may itself depend on a switch (clean_valley_mask
        // under write_valley_mask_raster, with six parameters under that), and a
        // dependent parameter may have dependants of its own. `written` doubles
        // as the cycle guard for a parent chain that loops back on itself.
        std::function<void(const std::string&, const std::string&)> write_tree =
            [&](const std::string& key, const std::string& indent)
        {
            if (!written.insert(key).second) return;
            write_overview_line(os, indent, key, entries.at(key));
            for (const auto& child_key : order)
            {
                const auto& c = entries.at(child_key);
                if (parent_of(c) == key && group_of(c) == group)
                {
                    write_tree(child_key, indent + "  ");
                }
            }
        };

        // Switches that nothing else claims as a dependant, with their trees.
        for (const auto& key : order)
        {
            const auto& e = entries.at(key);
            if (role_of(e) == "io" || group_of(e) != group) continue;
            if (role_of(e) != "switch") continue;
            const std::string p = parent_of(e);
            const auto pit = entries.find(p);
            if (!p.empty() && pit != entries.end() &&
                group_of(pit->second) == group) continue;   // renders as a child
            write_tree(key, "  ");
        }

        // Everything else in the group that the trees did not already write.
        // Also via write_tree, so a dependency chain nests whatever the roles
        // along it are: a parameter may depend on an algorithm choice rather
        // than on a switch (cellular_nonlocal_max_hops on
        // cellular_nonlocal_fraction), and that dependency is worth showing.
        for (const auto& key : order)
        {
            const auto& e = entries.at(key);
            if (role_of(e) == "io" || group_of(e) != group) continue;
            if (written.count(key)) continue;
            write_tree(key, "  ");
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

            // Report the geospatial runtime as well. A mismatched proj.db can
            // return coordinates that are plausible and wrong rather than
            // failing, so naming the one actually in use is what makes such a
            // report diagnosable. See src/core/lsdtt_gdal_init.hpp.
            const std::string proj_ver = proj_version();
            std::cout << "GDAL " << gdal_version();
            if (!proj_ver.empty())
                std::cout << ", PROJ " << proj_ver;
            std::cout << "\n"
                      << "  proj data: " << proj_data_description() << "\n"
                      << "  GDAL data: " << gdal_data_description() << "\n";
            std::exit(0);
        }
        if (is_flag(a, {"-c", "--cite", "citation", "-citation", "--citation"}))
        {
            std::cout << "Citation: " << citation << "\n";
            std::exit(0);
        }
    }

    // Anything that looks like a flag but is not one of the above, or more
    // than two positional arguments, means the caller is guessing at a
    // flag-style interface (e.g. "--input dem.tif"). Left alone, the first
    // would be taken as <path> and the second as the parameter file.
    for (int i = 1; i < argc; ++i)
    {
        const std::string a = argv[i];
        if (a.size() > 1 && a[0] == '-')
        {
            std::cerr << "Unrecognised option: " << a << "\n" << tool_name
                      << " takes no option flags except -h, -v and -c; all "
                         "other options go in the parameter file.\n\n";
            write_parameter_file_guide(std::cerr, tool_name, nullptr);
            std::exit(1);
        }
    }
    if (argc > 3)
    {
        std::cerr << "Too many arguments: " << tool_name
                  << " takes at most <path> <parameter_file>.\n\n";
        write_parameter_file_guide(std::cerr, tool_name, nullptr);
        std::exit(1);
    }

    if (argc == 1)
    {
        std::cerr << "No parameter file specified: " << tool_name
                  << " reads all of its options from a parameter file.\n\n";
        write_parameter_file_guide(std::cerr, tool_name, nullptr);
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
         << " my_params.driver</code> or <code>" << html_escape(tool_name)
         << " &lt;path&gt; my_params.driver</code></p>\n"
         << "<p>All options are set in a plain-text parameter file, not with "
         << "command-line flags. Format: one <code>key: value</code> per line; "
         << "keys are case-sensitive and omitted keys take their defaults; a "
         << "value contains no spaces (lists are comma-separated, e.g. "
         << "<code>a,b,c</code>); booleans are <code>true</code>/<code>false</code>; "
         << "<code>#</code> starts a comment.</p>\n";

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

    // Exactly as print_terminal_help: skip what the trees WROTE, never what a
    // predicate over parent role and group says they should have written. See
    // the long comment there for why the two must not be derived separately.
    std::set<std::string> written;
    for (const auto& key : order)
    {
        if (role_of(entries.at(key)) == "io") written.insert(key);
    }

    for (const auto& group : group_order)
    {
        html << "<h2>" << html_escape(group) << "</h2>\n" << table_header;

        std::function<void(const std::string&, bool)> write_tree =
            [&](const std::string& key, bool indent)
        {
            if (!written.insert(key).second) return;
            write_row(key, entries.at(key), indent);
            for (const auto& child_key : order)
            {
                const auto& c = entries.at(child_key);
                if (parent_of(c) == key && group_of(c) == group)
                {
                    write_tree(child_key, true);
                }
            }
        };

        for (const auto& key : order)
        {
            const auto& e = entries.at(key);
            if (role_of(e) == "io" || group_of(e) != group) continue;
            if (role_of(e) != "switch") continue;
            const std::string p = parent_of(e);
            const auto pit = entries.find(p);
            if (!p.empty() && pit != entries.end() &&
                group_of(pit->second) == group) continue;
            write_tree(key, false);
        }
        for (const auto& key : order)
        {
            const auto& e = entries.at(key);
            if (role_of(e) == "io" || group_of(e) != group) continue;
            if (written.count(key)) continue;
            write_tree(key, false);
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
        if (compatibility_storage_keys_.count(key)) continue;
        out << key << ": " << val << "\n";
    }

    out << "\n# Default parameters used\n";
    for (const auto& [key, val] : defaults_used_map_)
    {
        if (compatibility_storage_keys_.count(key)) continue;
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
    // Parser-owned I/O never enters the typed parameter maps. Include it here
    // centrally, with this run's resolved values, so per-run README files do
    // not depend on every CLI redeclaring the same four keys.
    std::map<std::string, std::vector<std::string>> full_help = help_map;
    full_help["read_path"] = {"string", read_path_,
        "Resolved directory from which input files are read.",
        "Defaults to the directory holding the parameter file.",
        "File I/O", "io", ""};
    full_help["read_fname"] = {"string", read_fname_ + read_ext_,
        "Resolved input filename, including its extension.",
        "Read from read_path; defaults to the parameter filename stem when omitted.",
        "File I/O", "io", ""};
    full_help["write_path"] = {"string", write_path_,
        "Resolved directory to which output files are written.",
        "Defaults to read_path when omitted.",
        "File I/O", "io", ""};
    full_help["write_fname"] = {"string", write_fname_,
        "Resolved output filename prefix (no path or extension).",
        "Defaults to the stem of read_fname when omitted.",
        "File I/O", "io", ""};

    // Collect all keys from read and defaults maps
    std::vector<std::string> all_keys = {
        "read_path", "read_fname", "write_path", "write_fname"};
    for (const auto& [k, v] : parameters_read_map_)
        if (!compatibility_storage_keys_.count(k)) all_keys.push_back(k);
    for (const auto& [k, v] : defaults_used_map_)
        if (!compatibility_storage_keys_.count(k)) all_keys.push_back(k);
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
        // Every field goes through csv_quote(). Descriptions and guidance
        // routinely contain commas, and writing them raw produced a file that
        // no CSV reader could parse: the group column picked up sentence
        // fragments and every column after it shifted. The sibling writer above
        // already quoted; this one did not.
        csv << "name,type,default,description,guidance,group,role,parent\n";
        csv << "version," << csv_quote(version) << ",citation,"
            << csv_quote(citation)
            << ",see code output for additional algorithm citations,,,\n";

        for (const auto& key : all_keys)
        {
            auto it = full_help.find(key);
            if (it == full_help.end() || it->second.size() < 4)
            {
                csv << csv_quote(key) << ",not available,not available,"
                    << "not available,not available,,,\n";
            }
            else
            {
                const auto& h = it->second;
                csv << csv_quote(key) << ","
                    << csv_quote(h[0]) << ","
                    << csv_quote(h[1]) << ","
                    << csv_quote(h[2]) << ","
                    << csv_quote(h[3]) << ","
                    << csv_quote(lsdtt_help_registry::field(h, 4)) << ","
                    << csv_quote(lsdtt_help_registry::field(h, 5)) << ","
                    << csv_quote(lsdtt_help_registry::field(h, 6)) << "\n";
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
            auto it = full_help.find(key);
            if (it == full_help.end() || it->second.size() < 4)
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
    std::ofstream out = open_output_file(
        output_fname, "lsdtt_param_parser::replace_and_print_parameter_file: cannot open ");
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

            // A value ends at the first whitespace.  Anything further on the
            // line (other than a # comment) therefore used to be DISCARDED
            // SILENTLY -- so "polyfit_metrics: a, b, c" ingested as "a," and
            // the run computed one metric instead of three, exit 0, no warning.
            // Reject it instead: a value the user wrote and we cannot honour is
            // an error, not a default.  Trailing whitespace and trailing
            // comments remain valid.
            if (c == ' ' || c == '\t')
            {
                std::string tail;
                char        t;
                while (infile.get(t) && t != '\n') tail += t;

                std::size_t first = tail.find_first_not_of(" \t\r");
                if (first != std::string::npos && tail[first] != '#')
                {
                    // Strip a trailing comment before echoing the value back.
                    std::string shown = tail.substr(first);
                    const std::size_t hash = shown.find('#');
                    if (hash != std::string::npos) shown = shown.substr(0, hash);
                    while (!shown.empty() &&
                           (shown.back() == ' ' || shown.back() == '\t' ||
                            shown.back() == '\r'))
                    {
                        shown.pop_back();
                    }

                    // Suggestion: the same value with ALL whitespace removed,
                    // which is the fix in the overwhelmingly common case of a
                    // comma-separated list typed with spaces after the commas.
                    std::string suggestion = value + shown;
                    suggestion.erase(
                        std::remove_if(suggestion.begin(), suggestion.end(),
                                       [](unsigned char ch) {
                                           return ch == ' ' || ch == '\t' ||
                                                  ch == '\r';
                                       }),
                        suggestion.end());

                    throw std::runtime_error(
                        "lsdtt_param_parser: parameter '" + parameter +
                        "' has a value containing spaces: '" + value + " " +
                        shown + "'.\n"
                        "  A value must be a single unbroken token: everything "
                        "after the first space is ignored,\n"
                        "  so this would silently run with only '" + value +
                        "'.\n"
                        "  If this is a comma-separated list, remove the "
                        "spaces:\n    " + parameter + ": " + suggestion
                    );
                }
                break;   // rest of the line was blank or a comment
            }
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
