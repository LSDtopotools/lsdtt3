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
//   read_path   - directory from which input files are read (default: ./)
//   write_path  - directory to which output files are written (default: read_path)
//   read_fname  - input filename (no path; may include extension, e.g. my_dem.tif)
//   write_fname - output filename prefix (no path, no extension)
//
// If read_fname includes an extension, the parser splits it internally into
// stem + extension. lsdtt3 uses GDAL for all raster I/O, so no separate file
// extension parameter is needed.
//
// Boundary conditions:
//   boundary_conditions - a string of characters, one per boundary segment.
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
#include <iosfwd>
#include <map>
#include <string>
#include <vector>

namespace lsdtt {

// -----------------------------------------------------------------------------
// lsdtt_cli_args
//
// Result of parse_args_deferred(). Unlike the legacy parse_args(), help
// requests do not exit immediately: the CLI declares its parameter maps
// first and then renders help from them, so -h can show real parameter
// documentation. See docs/developer_notes/lsdtt3_architecture.md (CLI design).
// -----------------------------------------------------------------------------
struct lsdtt_cli_args
{
    std::string path = "./";
    std::string fname;
    bool        help_requested = false;  // -h / --help
    bool        help_html      = false;  // -h html (write reference files)
    bool        help_json      = false;  // --help-json / -h json
    std::string help_topic;              // -h <parameter_name>
};

// -----------------------------------------------------------------------------
// lsdtt_compiled_features
//
// Returns the optional-dependency features compiled into the CALLING
// translation unit. Deliberately inline: the LSDTT_HAS_* macros are PUBLIC
// compile definitions on lsdtt3_cpu, so a CLI that links lsdtt3_cpu sees
// them here, while lsdtt3_core itself is compiled without them. Used as the
// default `features` argument of print_terminal_help(), which is evaluated
// at the call site.
// -----------------------------------------------------------------------------
inline std::vector<std::string> lsdtt_compiled_features()
{
    std::vector<std::string> features;
#ifdef LSDTT_HAS_FFTW
    features.push_back("fftw");
#endif
#ifdef LSDTT_HAS_EIGEN
    features.push_back("eigen");
#endif
#ifdef LSDTT_HAS_OPENCV
    features.push_back("opencv");
#endif
#ifdef LSDTT_HAS_NETCDF
    features.push_back("netcdf");
#endif
#ifdef LSDTT_GPU_KERNELS
    features.push_back("gpu_kernels");
#endif
    return features;
}

// -----------------------------------------------------------------------------
// lsdtt_help_registry
//
// Drop-in replacement for the plain help map
// (std::map<std::string, std::vector<std::string>>) that additionally
// records the order in which keys were first inserted. Declaration order
// follows the pipeline in every CLI, so help output can be grouped in
// pipeline order rather than alphabetically.
//
// Help entry format (positional, entries with only the first 4 elements
// remain valid):
//   0 type        "bool" / "int" / "float" / "string"
//   1 default     default value as string
//   2 short_desc  one-line description
//   3 long_desc   guidance paragraph
//   4 group       pipeline-step label (matches the // --- Section --- comment)
//   5 role        "switch" (enables an analysis/output), "param" (tunes one),
//                 "io" (paths, formats, threading)
//   6 parent      key of the switch this parameter belongs to ("" = none)
// -----------------------------------------------------------------------------
class lsdtt_help_registry
{
public:

    // Map-compatible accessor; records insertion order on first use of a key.
    std::vector<std::string>& operator[](const std::string& key);

    const std::map<std::string, std::vector<std::string>>& entries() const
    {
        return entries_;
    }

    // Keys in first-insertion order.
    const std::vector<std::string>& key_order() const { return key_order_; }

    bool empty() const { return entries_.empty(); }

    // Field access that tolerates short (4-element) entries: returns ""
    // when the index is out of range for the entry.
    static std::string field(const std::vector<std::string>& entry,
                             std::size_t index);

private:

    std::map<std::string, std::vector<std::string>> entries_;
    std::vector<std::string> key_order_;
};

class lsdtt_param_parser
{
public:

    // -------------------------------------------------------------------------
    // Construction
    // -------------------------------------------------------------------------

    // Default - produces an empty parser. Always prefer the path+file constructor.
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
    // e.g. "my_dem.tif" -> "my_dem"
    static std::string stem(const std::string& filename);

    // Returns the extension of a filename including the dot.
    // e.g. "my_dem.tif" -> ".tif"
    // Returns empty string if no dot is found.
    static std::string extension(const std::string& filename);

    // Parse command-line arguments for a standard lsdtt3 CLI
    // (see docs/developer_notes/lsdtt3_architecture.md, CLI design).
    // -v and -c print and exit. Help flags do NOT exit: they are
    // recorded in the returned struct so the CLI can declare its parameter
    // maps and then call print_terminal_help(). Accepted help forms:
    //   -h / --help            grouped parameter overview to stdout
    //   -h <parameter_name>    full guidance for one parameter
    //   -h html                write <tool>_help.html/.csv to the current dir
    //   --help-json / -h json  machine-readable JSON to stdout (tool,
    //                          version, citation, compiled features, and
    //                          all parameters with their grouped metadata)
    // When no help is requested, path/fname resolution requires a parameter
    // file (exit(1) when none is given).
    static lsdtt_cli_args parse_args_deferred(
        int argc, char* argv[],
        const std::string& tool_name,
        const std::string& version,
        const std::string& citation
    );

    // Render help for a CLI from its help registry. Static: needs no
    // parameter file and no parser instance. Behaviour depends on args:
    //   help_topic set  - detailed entry for that parameter (with close-match
    //                     suggestions when the name is unknown)
    //   help_html set   - write grouped <tool>_help.html/.csv to the current
    //                     directory
    //   otherwise       - grouped overview: io keys first, then one block per
    //                     group in declaration order, switches first with
    //                     their dependent parameters indented beneath
    //   help_json set   - machine-readable JSON to os
    // Output goes to os (default std::cout) except the html/csv files.
    // `features` defaults to lsdtt_compiled_features(), which is evaluated
    // at the call site so it reflects the CLI's compiled optional
    // dependencies; it is reported in the JSON output.
    static void print_terminal_help(
        const std::string& tool_name,
        const std::string& version,
        const std::string& citation,
        const lsdtt_help_registry& help,
        const lsdtt_cli_args& args,
        std::ostream& os,
        const std::vector<std::string>& features = lsdtt_compiled_features()
    );

    // Convenience overload writing to std::cout.
    static void print_terminal_help(
        const std::string& tool_name,
        const std::string& version,
        const std::string& citation,
        const lsdtt_help_registry& help,
        const lsdtt_cli_args& args,
        const std::vector<std::string>& features = lsdtt_compiled_features()
    );

    // Write the grouped HTML + CSV reference files for a tool to
    // out_prefix + ".html" / ".csv" (used by "-h html"; independent of any
    // run outputs).
    static void write_reference_files(
        const std::string& tool_name,
        const std::string& version,
        const std::string& citation,
        const lsdtt_help_registry& help,
        const std::string& out_prefix
    );

    // -------------------------------------------------------------------------
    // Output / diagnostics
    // -------------------------------------------------------------------------

    // Write all parsed parameters to a file so the user can verify ingestion.
    // Output is written to write_path + write_fname + "_ingested_params.param"
    void print_parameters() const;

    // Write to a named file.
    void print_parameters(const std::string& fname_prefix) const;

    // Write a help file (CSV + HTML) from a help map.
    //   help_map    - key -> {type, default, description, guidance}
    //   file_prefix - output filename prefix
    //   version     - version string to embed in the help file
    //   citation    - citation URL or string to embed
    void print_help(
        const std::map<std::string, std::vector<std::string>>& help_map,
        const std::string& file_prefix,
        const std::string& version,
        const std::string& citation
    ) const;

    // Overload taking a help registry; identical output (the end-of-run
    // files list the parameters of THIS run, flat table). Entries with
    // group/role/parent metadata get those values in appended columns.
    void print_help(
        const lsdtt_help_registry& help,
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
    // Getters - file I/O
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

    // Build a full output path: write_path + write_fname + suffix [+ ext].
    // Use when suffix already includes the extension:
    //   pp.make_output_path("_hillshade.tif")
    // or when suffix and extension are separate:
    //   pp.make_output_path("_hillshade", ".tif")
    std::string make_output_path(const std::string& suffix,
                                 const std::string& ext = "") const
    {
        return get_full_write_path() + suffix + ext;
    }

    // -------------------------------------------------------------------------
    // Getters - boundary conditions
    // -------------------------------------------------------------------------

    std::vector<std::string> get_boundary_conditions() const
    {
        return boundary_conditions_;
    }

    // -------------------------------------------------------------------------
    // Getters - typed parameter maps
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

    // File I/O - populated by parse_file_io()
    std::string read_path_;
    std::string write_path_;
    std::string read_fname_;   // stem only, no extension
    std::string read_ext_;     // extension including dot, e.g. ".tif"
    std::string write_fname_;  // stem only, no extension

    // Boundary conditions - populated by parse_file_io()
    std::vector<std::string> boundary_conditions_;

    // Raw string map of all key:value pairs read from file
    std::map<std::string, std::string> parameter_map_;

public:
    // True if the parameter file actually contained this key.
    //
    // Distinguishes "the user asked for the default value" from "the user said
    // nothing", which the typed maps cannot: both come back as the default.
    // Needed wherever a tool wants to CHOOSE a setting on the user's behalf but
    // must not override an explicit request -- selecting a router to match a
    // coupling regime, for instance, where silently overriding what was asked
    // for would produce a run that is not the experiment requested.
    [[nodiscard]] bool parameter_supplied(const std::string& key) const
    { return parameter_map_.find(key) != parameter_map_.end(); }

private:

    // Diagnostic maps - used for print_parameters() and print_help()
    std::map<std::string, std::string> defaults_used_map_;
    std::map<std::string, std::string> parameters_read_map_;

    // Typed parameter maps - populated by parse_*_parameters()
    std::map<std::string, float>       float_parameters_;
    std::map<std::string, int>         int_parameters_;
    std::map<std::string, bool>        bool_parameters_;
    std::map<std::string, std::string> string_parameters_;
    std::map<std::string, double>      double_parameters_;
};

} // namespace lsdtt
