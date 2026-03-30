// lsdtt_param_parser.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html
//
// -----------------------------------------------------------------------------
// lsdtt_param_parser
//
// Reads and parses lsdtt3 parameter files. Parameter files are plain text
// with the format:
//
//   # Lines beginning with # are comments and are ignored
//   key: value
//
// Keys use snake_case with no spaces. Values are strings which are converted
// to the requested type (float, int, bool, double, string) by comparing
// against caller-supplied default maps.
//
// File I/O parameters recognised directly by the parser:
//
//   read_path   — directory from which input files are read (default: ./)
//   write_path  — directory to which output files are written (default: read_path)
//   read_fname  — input filename prefix (no path, no extension)
//   write_fname — output filename prefix (no path, no extension)
//
// File extensions are not parsed here. lsdtt3 uses GDAL for all raster I/O
// and GDAL auto-detects file formats — no extension configuration is needed.
//
// Boundary conditions:
//   boundary_conditions — a string of characters, one per boundary segment.
//                         Valid characters: o (open), n (no flux), p (periodic)
//                         Example: boundary_conditions: ono
//
// Usage:
//
//   lsdtt::lsdtt_param_parser parser("./", "my_analysis.driver");
//
//   map<string,float>  defaults_f = { {"min_slope", 0.0001f} };
//   map<string,int>    defaults_i = { {"n_threads", 1} };
//   map<string,bool>   defaults_b = { {"fill_pits", true} };
//   map<string,string> defaults_s = { {"method", "wang_liu"} };
//
//   parser.parse_all_parameters(defaults_f, defaults_i, defaults_b, defaults_s);
//
//   float slope = parser.get_float_parameters()["min_slope"];
//   string path = parser.get_read_path();
//
// Original LSDTopoTools implementation by Simon M. Mudd and contributors.
// Redesigned for lsdtt3: removed extension logic, removed hardcoded file
// parameters, renamed keys to snake_case, removed external dependencies.
// -----------------------------------------------------------------------------

#pragma once

#include <fstream>
#include <map>
#include <string>
#include <vector>

namespace lsdtt {

class lsdtt_param_parser
{
public:

    // -------------------------------------------------------------------------
    // Construction
    // -------------------------------------------------------------------------

    // Default — produces an empty parser. Always prefer the path+file constructor.
    lsdtt_param_parser();

    // Construct from a path and filename. Immediately parses the file.
    // Throws std::runtime_error if the file cannot be opened.
    lsdtt_param_parser(const std::string& path, const std::string& filename);

    // -------------------------------------------------------------------------
    // Parsing
    // -------------------------------------------------------------------------

    // Parse the file into the internal string map. Called automatically
    // by the path+file constructor.
    void parse_file(const std::string& full_path);

    // Parse the file I/O parameters (read_path, write_path, read_fname,
    // write_fname, boundary_conditions) from the string map.
    // Called automatically by the path+file constructor after parse_file().
    void parse_file_io();

    // Parse all parameter types from their respective default maps.
    // Convenience wrapper around the four typed parse functions.
    void parse_all_parameters(
        const std::map<std::string, float>&       defaults_f,
        const std::map<std::string, int>&         defaults_i,
        const std::map<std::string, bool>&        defaults_b,
        const std::map<std::string, std::string>& defaults_s
    );

    // Overload including double parameters.
    void parse_all_parameters(
        const std::map<std::string, float>&       defaults_f,
        const std::map<std::string, int>&         defaults_i,
        const std::map<std::string, bool>&        defaults_b,
        const std::map<std::string, std::string>& defaults_s,
        const std::map<std::string, double>&      defaults_d
    );

    // Parse individual parameter types from their default maps.
    // Each function looks up each key in the parsed string map; if found it
    // converts and stores the value, otherwise it stores the default.
    void parse_float_parameters (const std::map<std::string, float>&       defaults);
    void parse_int_parameters   (const std::map<std::string, int>&         defaults);
    void parse_bool_parameters  (const std::map<std::string, bool>&        defaults);
    void parse_string_parameters(const std::map<std::string, std::string>& defaults);
    void parse_double_parameters(const std::map<std::string, double>&      defaults);

    // Force all entries in the string map into the parameters_read_map.
    // Used when copying parameter files.
    void force_parse();

    // -------------------------------------------------------------------------
    // Vector parameter parsing
    //
    // Comma-separated values stored in the string parameter map can be
    // unpacked into vectors. The key must have been parsed into
    // string_parameters first (i.e. listed in the defaults_s map).
    // -------------------------------------------------------------------------

    std::vector<std::string> parse_string_vector(const std::string& key) const;
    std::vector<int>         parse_int_vector   (const std::string& key) const;
    std::vector<float>       parse_float_vector (const std::string& key) const;
    std::vector<double>      parse_double_vector(const std::string& key) const;

    // -------------------------------------------------------------------------
    // Boundary conditions
    // -------------------------------------------------------------------------

    // Check that boundary conditions have been set. If not, defaults to
    // four no-flux boundaries ("n").
    void check_boundary_conditions();

    // -------------------------------------------------------------------------
    // File I/O path helpers
    // -------------------------------------------------------------------------

    // Check paths and filenames. If write_path is empty it defaults to
    // read_path; if read_fname is empty it defaults to the parameter filename
    // stem.
    void check_paths();

    // If this_string does not contain a path separator, prepend read_path.
    // Returns "NULL" unchanged.
    std::string add_read_path_if_needed(const std::string& this_string) const;

    // Returns the part of a string before the last dot.
    // e.g. "my_dem.tif" → "my_dem"
    static std::string stem(const std::string& filename);

    // Returns the extension of a filename including the dot.
    // e.g. "my_dem.tif" → ".tif"
    // Returns empty string if no dot is found.
    static std::string extension(const std::string& filename);

    // -------------------------------------------------------------------------
    // Output / diagnostics
    // -------------------------------------------------------------------------

    // Write all parsed parameters to a file so the user can verify ingestion.
    // Output is written to write_path + write_fname + "_ingested_params.param"
    void print_parameters() const;

    // Write to a named file.
    void print_parameters(const std::string& fname_prefix) const;

    // Write a help file (CSV + HTML) from a help map.
    //   help_map    — key → {type, default, description, guidance}
    //   file_prefix — output filename prefix
    //   version     — version string to embed in the help file
    //   citation    — citation URL or string to embed
    void print_help(
        const std::map<std::string, std::vector<std::string>>& help_map,
        const std::string& file_prefix,
        const std::string& version,
        const std::string& citation
    ) const;

    // Write a parameter file with some values replaced. Used for spawning
    // multiple analyses (e.g. one per catchment).
    void replace_and_print_parameter_file(
        const std::string& output_fname,
        const std::string& new_read_path,
        const std::string& new_read_fname,
        const std::string& new_write_path,
        const std::string& new_write_fname,
        const std::map<std::string, std::string>& replacements
    ) const;

    // -------------------------------------------------------------------------
    // Getters — file I/O
    // -------------------------------------------------------------------------

    std::string get_param_file_path() const { return param_file_path_; }
    std::string get_param_fname()     const { return param_fname_; }
    std::string get_read_path()       const { return read_path_; }
    std::string get_write_path()      const { return write_path_; }
    std::string get_read_fname()      const { return read_fname_; }
    std::string get_read_ext()        const { return read_ext_; }
    std::string get_write_fname()     const { return write_fname_; }

    // Returns the full input path including extension: read_path + read_fname + read_ext
    // e.g. "/data/dems/my_dem.tif"
    std::string get_full_read_path() const
    {
        return read_path_ + read_fname_ + read_ext_;
    }

    // Returns the write prefix: write_path + write_fname (no extension)
    // Append your own suffix before the extension, e.g.:
    //   pp.get_full_write_path() + "_hillshade.tif"
    std::string get_full_write_path() const
    {
        return write_path_ + write_fname_;
    }

    // -------------------------------------------------------------------------
    // Getters — boundary conditions
    // -------------------------------------------------------------------------

    std::vector<std::string> get_boundary_conditions() const
    {
        return boundary_conditions_;
    }

    // -------------------------------------------------------------------------
    // Getters — typed parameter maps
    // -------------------------------------------------------------------------

    std::map<std::string, float>       get_float_parameters()  const { return float_parameters_; }
    std::map<std::string, int>         get_int_parameters()    const { return int_parameters_; }
    std::map<std::string, bool>        get_bool_parameters()   const { return bool_parameters_; }
    std::map<std::string, std::string> get_string_parameters() const { return string_parameters_; }
    std::map<std::string, double>      get_double_parameters() const { return double_parameters_; }

    // -------------------------------------------------------------------------
    // Setters
    // -------------------------------------------------------------------------

    void set_read_fname (const std::string& s) { read_fname_  = s; }
    void set_write_fname(const std::string& s) { write_fname_ = s; }
    void set_boundary_conditions(const std::vector<std::string>& bc)
    {
        boundary_conditions_ = bc;
    }

private:

    // -------------------------------------------------------------------------
    // Internal helpers
    // -------------------------------------------------------------------------

    // Read one line from the parameter file, splitting on ':'.
    // parameter is the key (left of colon), value is the right.
    // Sets both to "NULL" for blank or comment lines.
    void parse_line(std::ifstream& infile,
                    std::string& parameter,
                    std::string& value) const;

    // Remove control characters (carriage returns etc.) from string end.
    // Needed for parameter files created on Windows.
    static std::string strip_control_chars(const std::string& s);

    // Ensure a path string ends with '/'.
    static std::string fix_path(const std::string& path);

    // Extract keys from a map into a sorted vector.
    template<typename V>
    static std::vector<std::string> keys_of(const std::map<std::string, V>& m);

    // -------------------------------------------------------------------------
    // Data members
    // -------------------------------------------------------------------------

    // Parameter file location
    std::string param_file_path_;
    std::string param_fname_;

    // File I/O — populated by parse_file_io()
    std::string read_path_;
    std::string write_path_;
    std::string read_fname_;   // stem only, no extension
    std::string read_ext_;     // extension including dot, e.g. ".tif"
    std::string write_fname_;  // stem only, no extension

    // Boundary conditions — populated by parse_file_io()
    std::vector<std::string> boundary_conditions_;

    // Raw string map of all key:value pairs read from file
    std::map<std::string, std::string> parameter_map_;

    // Diagnostic maps — used for print_parameters() and print_help()
    std::map<std::string, std::string> defaults_used_map_;
    std::map<std::string, std::string> parameters_read_map_;

    // Typed parameter maps — populated by parse_*_parameters()
    std::map<std::string, float>       float_parameters_;
    std::map<std::string, int>         int_parameters_;
    std::map<std::string, bool>        bool_parameters_;
    std::map<std::string, std::string> string_parameters_;
    std::map<std::string, double>      double_parameters_;
};

} // namespace lsdtt
