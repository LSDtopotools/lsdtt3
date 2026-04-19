// lsdtt-chi-analysis.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html
//
// -----------------------------------------------------------------------------
// lsdtt-chi-analysis
//
// End-to-end χ workflow for testing:
//
//   1. χ extraction on threshold-defined channel networks
//   2. χ-network fitting
//   3. χ disorder metrics and uncertainty products
//   4. spatial k_sn / u* inversion with m/n sweep
//
// Outlets can be supplied as a CSV of latitude/longitude points. If no outlet
// CSV is provided, the largest basin in the DEM is used automatically.
//
// Usage:
//   lsdtt-chi-analysis <path> <parameter_file>
//   lsdtt-chi-analysis <parameter_file>   (path assumed to be ./)
//   lsdtt-chi-analysis -h                 (print help and exit)
//   lsdtt-chi-analysis -v                 (print version and exit)
//   lsdtt-chi-analysis -c                 (print citation and exit)
// -----------------------------------------------------------------------------

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <initializer_list>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "../src/core/lsdtt_data_loader.hpp"
#include "../src/core/lsdtt_flow_base.hpp"
#include "../src/core/lsdtt_param_parser.hpp"
#include "../src/core/lsdtt_point_input.hpp"
#include "../src/core/lsdtt_types.hpp"
#include "../src/cpu/lsdtt_chi_disorder_cpu.hpp"
#include "../src/cpu/lsdtt_chi_extraction_cpu.hpp"
#include "../src/cpu/lsdtt_chi_ksn_inversion_cpu.hpp"
#include "../src/cpu/lsdtt_chi_network_analysis_cpu.hpp"
#include "../src/cpu/lsdtt_flow_snapping_cpu.hpp"
#include "../src/cpu/lsdtt_raster_cpu.hpp"
#include "../src/cpu/lsdtt_raster_flow_preprocessing_cpu.hpp"

static const std::string k_version  = "0.2.0";
static const std::string k_citation = "https://github.com/LSDtopotools/lsdtt3";

// -----------------------------------------------------------------------------
// print_banner
// -----------------------------------------------------------------------------
static void print_banner()
{
    std::cout
        << "=========================================================\n"
        << "|| Welcome to lsdtt-chi-analysis                       ||\n"
        << "|| I am going to test chi extraction and analysis      ||\n"
        << "|| extraction, fitting, disorder, and inversion        ||\n"
        << "|| Part of lsdtt3 by Simon M. Mudd                     ||\n"
        << "=========================================================\n"
        << "|| Version : " << k_version  << "\n"
        << "|| Citation: " << k_citation << "\n"
        << "=========================================================\n\n";
}

// -----------------------------------------------------------------------------
// parse_args
// -----------------------------------------------------------------------------
static std::pair<std::string, std::string> parse_args(int argc, char* argv[])
{
    auto is_flag = [](const std::string& s,
                      std::initializer_list<std::string> flags)
    {
        for (const auto& f : flags) if (s == f) return true;
        return false;
    };

    for (int i = 1; i < argc; ++i)
    {
        std::string a = argv[i];

        if (is_flag(a, {"-h", "--help", "help", "-H", "--H"}))
        {
            std::cout
                << "Usage: lsdtt-chi-analysis <path> <parameter_file>\n"
                << "       lsdtt-chi-analysis <parameter_file>\n"
                << "       lsdtt-chi-analysis -h   (this help)\n"
                << "       lsdtt-chi-analysis -v   (version)\n"
                << "       lsdtt-chi-analysis -c   (citation)\n";
            std::exit(0);
        }
        if (is_flag(a, {"-v", "--version", "version", "-V"}))
        {
            std::cout << "lsdtt-chi-analysis version " << k_version << "\n";
            std::exit(0);
        }
        if (is_flag(a, {"-c", "--cite", "citation", "-citation", "--citation"}))
        {
            std::cout << "Citation: " << k_citation << "\n";
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
// out_path
// -----------------------------------------------------------------------------
static std::string out_path(
    const lsdtt::lsdtt_param_parser& pp,
    const std::string& suffix,
    const std::string& extension)
{
    return pp.get_full_write_path() + suffix + extension;
}

// -----------------------------------------------------------------------------
// helpers
// -----------------------------------------------------------------------------
static std::string to_lower_copy(std::string s)
{
    for (char& ch : s)
    {
        ch = static_cast<char>(std::tolower(static_cast<unsigned char>(ch)));
    }
    return s;
}

static std::string mn_suffix(double mn)
{
    std::ostringstream ss;
    ss << std::fixed << std::setprecision(3) << mn;
    std::string s = ss.str();
    for (char& c : s)
    {
        if (c == '.') c = 'p';
        if (c == '-') c = 'm';
    }
    return s;
}

static std::vector<double> build_mn_vector(double mn_min, double mn_max, double mn_step)
{
    if (!(std::isfinite(mn_min) && std::isfinite(mn_max) && std::isfinite(mn_step)))
    {
        throw std::invalid_argument("m/n sweep parameters must be finite.");
    }
    if (mn_step <= 0.0)
    {
        throw std::invalid_argument("mn_step must be > 0.");
    }
    if (mn_max < mn_min)
    {
        throw std::invalid_argument("mn_max must be >= mn_min.");
    }

    std::vector<double> tested;
    for (double mn = mn_min; mn <= mn_max + 0.5 * mn_step; mn += mn_step)
    {
        tested.push_back(mn);
    }
    return tested;
}

static std::int64_t largest_basin_outlet_flat(const lsdtt::lsdtt_flow_base& flow)
{
    std::int64_t best_flat = lsdtt::k_nodata_receiver;
    std::int64_t best_size = -1;

    for (int row = 0; row < flow.n_rows(); ++row)
    {
        for (int col = 0; col < flow.n_cols(); ++col)
        {
            if (!flow.is_valid_cell(row, col))
            {
                continue;
            }

            const std::int64_t flat = flow.flat_index_unchecked(row, col);
            if (flow.receiver_of_unchecked(flat) != lsdtt::k_outlet_receiver)
            {
                continue;
            }

            const std::int64_t n_up =
                static_cast<std::int64_t>(flow.s_vector_size_unchecked(flat));

            if (n_up > best_size)
            {
                best_size = n_up;
                best_flat = flat;
            }
        }
    }

    if (best_flat == lsdtt::k_nodata_receiver)
    {
        throw std::runtime_error(
            "Could not identify a largest basin outlet in the flow object."
        );
    }

    return best_flat;
}

static std::vector<std::int64_t> resolve_outlet_flats(
    const lsdtt::lsdtt_param_parser& pp,
    const lsdtt::lsdtt_raster_cpu& dem,
    const lsdtt::lsdtt_flow_base& flow,
    const std::string& outlet_csv,
    const std::string& snapping_method,
    int snap_threshold_pixels,
    int min_search_radius_pixels,
    int max_search_radius_pixels)
{
    if (outlet_csv == "NULL" || outlet_csv.empty())
    {
        const std::int64_t flat = largest_basin_outlet_flat(flow);
        std::cout << "No outlet CSV supplied; using largest basin outlet flat: "
                  << flat << "\n";
        return {flat};
    }

    const std::string csv_path = pp.add_read_path_if_needed(outlet_csv);

    std::cout << "Reading outlet CSV: " << csv_path << "\n";

    std::vector<lsdtt::lsdtt_input_point> requested_outlets =
        lsdtt::load_csv_latitude_longitude_points_in_raster_crs(
            csv_path,
            dem.get_meta()
        );

    if (requested_outlets.empty())
    {
        throw std::runtime_error("Outlet CSV contained no usable points.");
    }

    lsdtt::lsdtt_square_snap_options square_snap_options;
    square_snap_options.contributing_pixels_threshold =
        static_cast<std::uint64_t>(snap_threshold_pixels);
    square_snap_options.min_search_radius_pixels = min_search_radius_pixels;
    square_snap_options.max_search_radius_pixels = max_search_radius_pixels;

    lsdtt::lsdtt_downstream_snap_options downstream_snap_options;
    downstream_snap_options.contributing_pixels_threshold =
        static_cast<std::uint64_t>(snap_threshold_pixels);

    std::vector<std::int64_t> outlet_flats;
    outlet_flats.reserve(requested_outlets.size());

    for (const auto& pt : requested_outlets)
    {
        if (pt.status != lsdtt::lsdtt_point_input_status::success)
        {
            std::cerr << "Warning: outlet point " << pt.point_id
                      << " did not map to a valid raster cell; skipping.\n";
            continue;
        }

        int row = pt.row;
        int col = pt.col;
        std::int64_t flat = pt.flat;

        if (snapping_method == "expanding_square")
        {
            const lsdtt::lsdtt_snap_result snapped =
                lsdtt::snap_pixel_to_contributing_pixels_threshold_expanding_square(
                    flow,
                    row,
                    col,
                    square_snap_options
                );

            if (snapped.status != lsdtt::lsdtt_snap_status::success)
            {
                std::cerr << "Warning: outlet point " << pt.point_id
                          << " could not be snapped with expanding_square; skipping.\n";
                continue;
            }

            flat = snapped.snapped_flat;
        }
        else if (snapping_method == "downstream")
        {
            const lsdtt::lsdtt_snap_result snapped =
                lsdtt::snap_pixel_to_contributing_pixels_threshold_downstream(
                    flow,
                    row,
                    col,
                    downstream_snap_options
                );

            if (snapped.status != lsdtt::lsdtt_snap_status::success)
            {
                std::cerr << "Warning: outlet point " << pt.point_id
                          << " could not be snapped downstream; skipping.\n";
                continue;
            }

            flat = snapped.snapped_flat;
        }
        else
        {
            if (!flow.is_valid_cell(flat))
            {
                std::cerr << "Warning: outlet point " << pt.point_id
                          << " maps to invalid flow cell; skipping.\n";
                continue;
            }
        }

        outlet_flats.push_back(flat);
    }

    if (outlet_flats.empty())
    {
        throw std::runtime_error("No valid outlet flats were produced from the CSV.");
    }

    return outlet_flats;
}

static std::int32_t count_unique_branch_ids(
    const lsdtt::lsdtt_chi_network_fit_result& fit_result)
{
    std::vector<std::int32_t> ids = fit_result.branch_id_of_node;
    std::sort(ids.begin(), ids.end());
    ids.erase(std::unique(ids.begin(), ids.end()), ids.end());
    return static_cast<std::int32_t>(ids.size());
}

static std::int32_t count_main_stem_nodes(
    const lsdtt::lsdtt_chi_network_fit_result& fit_result)
{
    std::int32_t count = 0;
    for (std::uint8_t f : fit_result.is_main_stem)
    {
        if (f != 0u) ++count;
    }
    return count;
}

static std::string join_branch_ids(const std::vector<std::int32_t>& ids)
{
    std::ostringstream ss;
    for (std::size_t i = 0; i < ids.size(); ++i)
    {
        if (i > 0) ss << ';';
        ss << ids[i];
    }
    return ss.str();
}

// -----------------------------------------------------------------------------
// CSV writers
// -----------------------------------------------------------------------------
static void write_fit_node_table_csv(
    const std::string& filename,
    const lsdtt::lsdtt_chi_network_fit_node_table& table)
{
    std::ofstream out(filename);
    if (!out)
    {
        throw std::runtime_error("Could not open file for writing: " + filename);
    }

    out << "basin_index,outlet_flat,flat,receiver_flat,chi,observed_elevation,"
           "fitted_elevation,residual_elevation,abs_residual_elevation,branch_id,"
           "is_main_stem,is_source,is_junction,is_split_point\n";

    const std::size_t n = table.flat.size();
    for (std::size_t i = 0; i < n; ++i)
    {
        out << table.basin_index[i] << ','
            << table.outlet_flat[i] << ','
            << table.flat[i] << ','
            << table.receiver_flat[i] << ','
            << table.chi[i] << ','
            << table.observed_elevation[i] << ','
            << table.fitted_elevation[i] << ','
            << table.residual_elevation[i] << ','
            << table.abs_residual_elevation[i] << ','
            << table.branch_id[i] << ','
            << static_cast<int>(table.is_main_stem[i]) << ','
            << static_cast<int>(table.is_source[i]) << ','
            << static_cast<int>(table.is_junction[i]) << ','
            << static_cast<int>(table.is_split_point[i]) << '\n';
    }
}

static void write_fit_segment_table_csv(
    const std::string& filename,
    const lsdtt::lsdtt_chi_network_fit_segment_table& table)
{
    std::ofstream out(filename);
    if (!out)
    {
        throw std::runtime_error("Could not open file for writing: " + filename);
    }

    out << "basin_index,outlet_flat,segment_id,parent_segment_id,branch_id,"
           "downstream_flat,upstream_flat,downstream_chi,upstream_chi,"
           "downstream_elevation,upstream_elevation,max_abs_residual,n_nodes\n";

    const std::size_t n = table.segment_id.size();
    for (std::size_t i = 0; i < n; ++i)
    {
        out << table.basin_index[i] << ','
            << table.outlet_flat[i] << ','
            << table.segment_id[i] << ','
            << table.parent_segment_id[i] << ','
            << table.branch_id[i] << ','
            << table.downstream_flat[i] << ','
            << table.upstream_flat[i] << ','
            << table.downstream_chi[i] << ','
            << table.upstream_chi[i] << ','
            << table.downstream_elevation[i] << ','
            << table.upstream_elevation[i] << ','
            << table.max_abs_residual[i] << ','
            << table.n_nodes[i] << '\n';
    }
}

static void write_subset_table_csv(
    const std::string& filename,
    const lsdtt::lsdtt_chi_branch_subset_sensitivity_sample_table& table)
{
    std::ofstream out(filename);
    if (!out)
    {
        throw std::runtime_error("Could not open file for writing: " + filename);
    }

    out << "basin_index,full_network_disorder,sample_index,disorder,n_nodes,"
           "included_branch_ids\n";

    const std::size_t n = table.sample_index.size();
    for (std::size_t i = 0; i < n; ++i)
    {
        out << table.basin_index[i] << ','
            << table.full_network_disorder[i] << ','
            << table.sample_index[i] << ','
            << table.disorder[i] << ','
            << table.n_nodes[i] << ",\""
            << join_branch_ids(table.included_branch_ids[i]) << "\"\n";
    }
}

static void write_jackknife_table_csv(
    const std::string& filename,
    const lsdtt::lsdtt_chi_disorder_jackknife_table& table)
{
    std::ofstream out(filename);
    if (!out)
    {
        throw std::runtime_error("Could not open file for writing: " + filename);
    }

    out << "basin_index,full_network_disorder,replicate_index,removed_branch_id,"
           "leave_one_out_disorder,leave_one_out_n_nodes,jackknife_mean,"
           "jackknife_standard_error\n";

    const std::size_t n = table.replicate_index.size();
    for (std::size_t i = 0; i < n; ++i)
    {
        out << table.basin_index[i] << ','
            << table.full_network_disorder[i] << ','
            << table.replicate_index[i] << ','
            << table.removed_branch_id[i] << ','
            << table.leave_one_out_disorder[i] << ','
            << table.leave_one_out_n_nodes[i] << ','
            << table.jackknife_mean[i] << ','
            << table.jackknife_standard_error[i] << '\n';
    }
}

static void write_bootstrap_table_csv(
    const std::string& filename,
    const lsdtt::lsdtt_chi_disorder_bootstrap_table& table)
{
    std::ofstream out(filename);
    if (!out)
    {
        throw std::runtime_error("Could not open file for writing: " + filename);
    }

    out << "basin_index,full_network_disorder,replicate_index,bootstrap_disorder,"
           "bootstrap_n_nodes,bootstrap_mean,bootstrap_standard_deviation,"
           "percentile_2p5,percentile_97p5\n";

    const std::size_t n = table.replicate_index.size();
    for (std::size_t i = 0; i < n; ++i)
    {
        out << table.basin_index[i] << ','
            << table.full_network_disorder[i] << ','
            << table.replicate_index[i] << ','
            << table.bootstrap_disorder[i] << ','
            << table.bootstrap_n_nodes[i] << ','
            << table.bootstrap_mean[i] << ','
            << table.bootstrap_standard_deviation[i] << ','
            << table.percentile_2p5[i] << ','
            << table.percentile_97p5[i] << '\n';
    }
}

struct disorder_sweep_row
{
    double mn = std::numeric_limits<double>::quiet_NaN();
    std::int32_t basin_index = -1;
    std::int64_t outlet_flat = lsdtt::k_nodata_receiver;

    std::int32_t n_nodes = 0;
    std::int32_t n_main_stem_nodes = 0;
    std::int32_t n_branches = 0;
    std::int32_t n_side_branches = 0;

    float full_network_disorder = std::numeric_limits<float>::quiet_NaN();

    std::int32_t n_subset_sensitivity_samples = 0;

    float jackknife_mean = std::numeric_limits<float>::quiet_NaN();
    float jackknife_standard_error = std::numeric_limits<float>::quiet_NaN();

    float bootstrap_mean = std::numeric_limits<float>::quiet_NaN();
    float bootstrap_standard_deviation = std::numeric_limits<float>::quiet_NaN();
    float bootstrap_percentile_2p5 = std::numeric_limits<float>::quiet_NaN();
    float bootstrap_percentile_97p5 = std::numeric_limits<float>::quiet_NaN();
};

static void write_disorder_sweep_csv(
    const std::string& filename,
    const std::vector<disorder_sweep_row>& rows)
{
    std::ofstream out(filename);
    if (!out)
    {
        throw std::runtime_error("Could not open file for writing: " + filename);
    }

    out << "m_over_n,basin_index,outlet_flat,n_nodes,n_main_stem_nodes,n_branches,"
           "n_side_branches,full_network_disorder,n_subset_sensitivity_samples,"
           "jackknife_mean,jackknife_standard_error,bootstrap_mean,"
           "bootstrap_standard_deviation,bootstrap_percentile_2p5,"
           "bootstrap_percentile_97p5\n";

    for (const auto& r : rows)
    {
        out << r.mn << ','
            << r.basin_index << ','
            << r.outlet_flat << ','
            << r.n_nodes << ','
            << r.n_main_stem_nodes << ','
            << r.n_branches << ','
            << r.n_side_branches << ','
            << r.full_network_disorder << ','
            << r.n_subset_sensitivity_samples << ','
            << r.jackknife_mean << ','
            << r.jackknife_standard_error << ','
            << r.bootstrap_mean << ','
            << r.bootstrap_standard_deviation << ','
            << r.bootstrap_percentile_2p5 << ','
            << r.bootstrap_percentile_97p5 << '\n';
    }
}

static void write_inversion_trial_csv(
    const std::string& filename,
    const lsdtt::lsdtt_chi_ksn_inversion_collection_result& results)
{
    std::ofstream out(filename);
    if (!out)
    {
        throw std::runtime_error("Could not open file for writing: " + filename);
    }

    out << "basin_index,outlet_flat,trial_index,m_over_n,mean_absolute_error,"
           "root_mean_square_error,n_nodes,n_parameters,n_constrained_parameters\n";

    for (std::size_t b = 0; b < results.basins.size(); ++b)
    {
        const auto& basin = results.basins[b];
        for (std::size_t t = 0; t < basin.trials.size(); ++t)
        {
            const auto& tr = basin.trials[t];
            out << static_cast<std::int32_t>(b) << ','
                << basin.network.outlet_flat << ','
                << static_cast<std::int32_t>(t) << ','
                << tr.m_over_n << ','
                << tr.mean_absolute_error << ','
                << tr.root_mean_square_error << ','
                << tr.n_nodes << ','
                << tr.n_parameters << ','
                << tr.n_constrained_parameters << '\n';
        }
    }
}

static void write_inversion_best_summary_csv(
    const std::string& filename,
    const lsdtt::lsdtt_chi_ksn_inversion_collection_result& results)
{
    std::ofstream out(filename);
    if (!out)
    {
        throw std::runtime_error("Could not open file for writing: " + filename);
    }

    out << "basin_index,outlet_flat,best_trial_index,best_m_over_n,"
           "best_mean_absolute_error,best_root_mean_square_error,"
           "n_nodes_best,n_parameters_best,n_constrained_parameters_best\n";

    for (std::size_t b = 0; b < results.basins.size(); ++b)
    {
        const auto& basin = results.basins[b];

        std::int32_t n_nodes_best = 0;
        std::int32_t n_params_best = 0;
        std::int32_t n_const_best = 0;

        if (basin.best_trial_index >= 0 &&
            static_cast<std::size_t>(basin.best_trial_index) < basin.trials.size())
        {
            const auto& tr = basin.trials[static_cast<std::size_t>(basin.best_trial_index)];
            n_nodes_best = tr.n_nodes;
            n_params_best = tr.n_parameters;
            n_const_best = tr.n_constrained_parameters;
        }

        out << static_cast<std::int32_t>(b) << ','
            << basin.network.outlet_flat << ','
            << basin.best_trial_index << ','
            << basin.best_m_over_n << ','
            << basin.best_mean_absolute_error << ','
            << basin.best_root_mean_square_error << ','
            << n_nodes_best << ','
            << n_params_best << ','
            << n_const_best << '\n';
    }
}

// =============================================================================
// main
// =============================================================================
int main(int argc, char* argv[])
{
    print_banner();

    auto [path, fname] = parse_args(argc, argv);

    // -------------------------------------------------------------------------
    // Parameter parsing
    // -------------------------------------------------------------------------
    lsdtt::lsdtt_param_parser pp;
    try
    {
        pp = lsdtt::lsdtt_param_parser(path, fname);
    }
    catch (const std::exception& e)
    {
        std::cerr << "Error loading parameter file: " << e.what() << "\n";
        return 1;
    }

    std::map<std::string, bool>                        bool_default_map;
    std::map<std::string, float>                       float_default_map;
    std::map<std::string, int>                         int_default_map;
    std::map<std::string, std::string>                 string_default_map;
    std::map<std::string, std::vector<std::string>>    help_map;

    // Flow / DEM prep
    int_default_map["n_threads"] = 0;
    string_default_map["boundary_conditions"] = "oooo";
    bool_default_map["raster_is_filled"] = false;
    float_default_map["fill_min_slope"] = 0.0001f;

    // Outlet selection
    string_default_map["outlet_csv"] = "NULL";
    string_default_map["outlet_snapping_method"] = "downstream";
    int_default_map["outlet_threshold_contributing_pixels"] = 1000;
    int_default_map["outlet_min_search_radius_pixels"] = 0;
    int_default_map["outlet_max_search_radius_pixels"] = 25;

    // Chi extraction / fitting
    int_default_map["channel_threshold_contributing_pixels"] = 1000;
    float_default_map["chi_A0"] = 1.0f;
    float_default_map["mn_min"] = 0.3f;
    float_default_map["mn_max"] = 0.6f;
    float_default_map["mn_step"] = 0.02f;
    float_default_map["fit_critical_divergence_m"] = 5.0f;

    // Disorder
    bool_default_map["run_subset_sensitivity"] = true;
    int_default_map["subset_n_side_branches_per_sample"] = 3;

    bool_default_map["run_jackknife"] = true;

    bool_default_map["run_bootstrap"] = true;
    int_default_map["bootstrap_replicates"] = 1000;
    bool_default_map["bootstrap_sample_full_tributary_count"] = true;
    int_default_map["bootstrap_n_side_branches_per_replicate"] = 3;
    bool_default_map["bootstrap_with_replacement"] = true;
    int_default_map["bootstrap_random_seed"] = 5489;

    // Inversion
    bool_default_map["run_inversion"] = true;
    float_default_map["inversion_coarse_cell_size"] = 2000.0f;
    float_default_map["inversion_roughness_weight"] = 1.0e4f;
    float_default_map["inversion_diagonal_jitter"] = 1.0e-8f;
    bool_default_map["inversion_use_mae_for_model_selection"] = true;
    bool_default_map["write_inversion_rasters"] = true;

    // Outputs
    bool_default_map["write_fit_tables"] = false;
    bool_default_map["write_disorder_sample_tables"] = false;
    string_default_map["write_output_format"] = "cog";

    pp.parse_all_parameters(
        float_default_map,
        int_default_map,
        bool_default_map,
        string_default_map
    );

    auto bf = pp.get_bool_parameters();
    auto ff = pp.get_float_parameters();
    auto ifs = pp.get_int_parameters();
    auto sf = pp.get_string_parameters();

    const std::string snapping_method = to_lower_copy(sf["outlet_snapping_method"]);
    const lsdtt::raster_output_format out_fmt =
        lsdtt::resolve_raster_output_format(sf["write_output_format"]);

    // -------------------------------------------------------------------------
    // Validate params
    // -------------------------------------------------------------------------
    if (ifs["channel_threshold_contributing_pixels"] <= 0)
    {
        std::cerr << "FATAL: channel_threshold_contributing_pixels must be > 0.\n";
        return 1;
    }

    if (snapping_method != "none" &&
        snapping_method != "expanding_square" &&
        snapping_method != "downstream")
    {
        std::cerr << "FATAL: outlet_snapping_method must be one of: none, expanding_square, downstream.\n";
        return 1;
    }

    std::vector<double> tested_mn;
    try
    {
        tested_mn = build_mn_vector(
            static_cast<double>(ff["mn_min"]),
            static_cast<double>(ff["mn_max"]),
            static_cast<double>(ff["mn_step"])
        );
    }
    catch (const std::exception& e)
    {
        std::cerr << "FATAL: invalid m/n sweep: " << e.what() << "\n";
        return 1;
    }

    // -------------------------------------------------------------------------
    // Open DEM
    // -------------------------------------------------------------------------
    const std::string input_path = pp.get_full_read_path();
    std::cout << "\nOpening DEM: " << input_path << "\n";

    lsdtt::lsdtt_data_loader loader;
    try
    {
        loader.open_raster(input_path);
    }
    catch (const std::exception& e)
    {
        std::cerr << "Error opening DEM: " << e.what() << "\n";
        return 1;
    }

    loader.print_raster_info();

    try
    {
        loader.get_meta().check_projected_metres();
    }
    catch (const std::exception& e)
    {
        std::cerr << "\nFATAL: " << e.what() << "\n";
        return 1;
    }

    lsdtt::lsdtt_raster_cpu dem(loader);
    lsdtt::lsdtt_raster_cpu routed_dem = dem;

    // -------------------------------------------------------------------------
    // Optional fill
    // -------------------------------------------------------------------------
    if (!bf["raster_is_filled"])
    {
        std::cout << "\n--- Pit filling ---\n";
        try
        {
            routed_dem = lsdtt::fill_pits(
                dem,
                ff["fill_min_slope"],
                ifs["n_threads"]
            );
        }
        catch (const std::exception& e)
        {
            std::cerr << "FATAL: pit filling failed: " << e.what() << "\n";
            return 1;
        }
    }
    else
    {
        std::cout << "\nRaster is marked as already filled. Skipping fill step.\n";
    }

    // -------------------------------------------------------------------------
    // Flow routing
    // -------------------------------------------------------------------------
    std::cout << "\n--- Flow routing ---\n";

    lsdtt::lsdtt_flowinfo_options fi_opts;
    fi_opts.boundary_conditions       = sf["boundary_conditions"];
    fi_opts.interior_nodata_is_outlet = true;
    fi_opts.require_acyclic_graph     = true;
    fi_opts.build_topological_order   = true;
    fi_opts.build_s_vector            = true;
    fi_opts.n_threads                 = ifs["n_threads"];

    std::unique_ptr<lsdtt::lsdtt_flow_base> flow;
    try
    {
        flow = lsdtt::make_flow_base(
            routed_dem.get_data_array(),
            routed_dem.get_mask_array(),
            routed_dem.get_meta(),
            fi_opts
        );
    }
    catch (const std::exception& e)
    {
        std::cerr << "FATAL: flow routing failed: " << e.what() << "\n";
        return 1;
    }

    // -------------------------------------------------------------------------
    // Resolve outlets
    // -------------------------------------------------------------------------
    std::vector<std::int64_t> outlet_flats;
    try
    {
        outlet_flats = resolve_outlet_flats(
            pp,
            dem,
            *flow,
            sf["outlet_csv"],
            snapping_method,
            ifs["outlet_threshold_contributing_pixels"],
            ifs["outlet_min_search_radius_pixels"],
            ifs["outlet_max_search_radius_pixels"]
        );
    }
    catch (const std::exception& e)
    {
        std::cerr << "FATAL: outlet resolution failed: " << e.what() << "\n";
        return 1;
    }

    std::cout << "Number of outlet seeds: " << outlet_flats.size() << "\n";

    // -------------------------------------------------------------------------
    // Build inversion network support once.
    //
    // The inversion uses a fixed network geometry/topology, so we extract the
    // threshold channel network once using the first tested m/n. Chi itself is
    // recomputed internally during inversion, so only flats/receiver topology
    // and elevations matter here.
    // -------------------------------------------------------------------------
    std::vector<lsdtt::lsdtt_chi_ksn_inversion_network> inversion_networks;

    if (bf["run_inversion"])
    {
        std::cout << "\n--- Preparing inversion networks ---\n";
        try
        {
            lsdtt::lsdtt_chi_options chi_opts;
            chi_opts.m_over_n = tested_mn.front();
            chi_opts.A_0 = static_cast<double>(ff["chi_A0"]);

            lsdtt::lsdtt_chi_collection_result chi_networks =
                lsdtt::get_chi_from_outlets_zero_base_threshold(
                    *flow,
                    outlet_flats,
                    ifs["channel_threshold_contributing_pixels"],
                    chi_opts
                );

            inversion_networks.reserve(chi_networks.basins.size());
            for (const auto& basin : chi_networks.basins)
            {
                inversion_networks.push_back(
                    lsdtt::make_chi_ksn_inversion_network(basin, dem)
                );
            }
        }
        catch (const std::exception& e)
        {
            std::cerr << "FATAL: preparing inversion networks failed: "
                      << e.what() << "\n";
            return 1;
        }
    }

    // -------------------------------------------------------------------------
    // χ extraction + fitting + disorder sweep
    // -------------------------------------------------------------------------
    std::vector<disorder_sweep_row> disorder_rows;

    for (double mn : tested_mn)
    {
        std::cout << "\n=====================================================\n";
        std::cout << "Testing m/n = " << mn << "\n";
        std::cout << "=====================================================\n";

        lsdtt::lsdtt_chi_options chi_opts;
        chi_opts.m_over_n = mn;
        chi_opts.A_0 = static_cast<double>(ff["chi_A0"]);

        lsdtt::lsdtt_chi_collection_result chi_networks;
        lsdtt::lsdtt_chi_network_fit_collection_result fit_results;

        try
        {
            chi_networks = lsdtt::get_chi_from_outlets_zero_base_threshold(
                *flow,
                outlet_flats,
                ifs["channel_threshold_contributing_pixels"],
                chi_opts
            );

            lsdtt::lsdtt_chi_network_fit_options fit_opts;
            fit_opts.critical_divergence_m = ff["fit_critical_divergence_m"];
            fit_opts.mainstem_choice = lsdtt::lsdtt_chi_mainstem_choice::max_flow_length;

            fit_results = lsdtt::fit_chi_networks(
                chi_networks,
                dem,
                fit_opts
            );
        }
        catch (const std::exception& e)
        {
            std::cerr << "FATAL: chi extraction/fitting failed for m/n = "
                      << mn << ": " << e.what() << "\n";
            return 1;
        }

        if (bf["write_fit_tables"])
        {
            try
            {
                const auto node_table = lsdtt::export_chi_network_fit_nodes(fit_results);
                const auto seg_table  = lsdtt::export_chi_network_fit_segments(fit_results);

                write_fit_node_table_csv(
                    out_path(pp, "_chi_fit_nodes_mn_" + mn_suffix(mn), ".csv"),
                    node_table
                );
                write_fit_segment_table_csv(
                    out_path(pp, "_chi_fit_segments_mn_" + mn_suffix(mn), ".csv"),
                    seg_table
                );
            }
            catch (const std::exception& e)
            {
                std::cerr << "Warning: writing fit tables failed for m/n = "
                          << mn << ": " << e.what() << "\n";
            }
        }

        lsdtt::lsdtt_chi_disorder_collection_result disorder_results;
        lsdtt::lsdtt_chi_branch_subset_sensitivity_collection_result subset_results;
        lsdtt::lsdtt_chi_disorder_jackknife_collection_result jackknife_results;
        lsdtt::lsdtt_chi_disorder_bootstrap_collection_result bootstrap_results;

        try
        {
            disorder_results = lsdtt::get_chi_disorder(fit_results);

            if (bf["run_subset_sensitivity"])
            {
                lsdtt::lsdtt_chi_branch_subset_sensitivity_options opts;
                opts.n_side_branches_per_sample =
                    ifs["subset_n_side_branches_per_sample"];
                subset_results = lsdtt::get_chi_branch_subset_sensitivity(
                    fit_results, opts
                );
            }

            if (bf["run_jackknife"])
            {
                jackknife_results = lsdtt::get_chi_disorder_jackknife(fit_results);
            }

            if (bf["run_bootstrap"])
            {
                lsdtt::lsdtt_chi_disorder_bootstrap_options opts;
                opts.n_replicates = ifs["bootstrap_replicates"];
                opts.sample_full_tributary_count =
                    bf["bootstrap_sample_full_tributary_count"];
                opts.n_side_branches_per_replicate =
                    ifs["bootstrap_n_side_branches_per_replicate"];
                opts.sample_with_replacement =
                    bf["bootstrap_with_replacement"];
                opts.random_seed =
                    static_cast<std::uint64_t>(ifs["bootstrap_random_seed"]);

                bootstrap_results = lsdtt::get_chi_disorder_bootstrap(
                    fit_results, opts
                );
            }
        }
        catch (const std::exception& e)
        {
            std::cerr << "FATAL: disorder analysis failed for m/n = "
                      << mn << ": " << e.what() << "\n";
            return 1;
        }

        // Optionally write sample tables for this m/n
        if (bf["write_disorder_sample_tables"])
        {
            try
            {
                if (bf["run_subset_sensitivity"])
                {
                    write_subset_table_csv(
                        out_path(pp, "_chi_disorder_subset_samples_mn_" + mn_suffix(mn), ".csv"),
                        lsdtt::export_chi_branch_subset_sensitivity_samples(subset_results)
                    );
                }

                if (bf["run_jackknife"])
                {
                    write_jackknife_table_csv(
                        out_path(pp, "_chi_disorder_jackknife_mn_" + mn_suffix(mn), ".csv"),
                        lsdtt::export_chi_disorder_jackknife(jackknife_results)
                    );
                }

                if (bf["run_bootstrap"])
                {
                    write_bootstrap_table_csv(
                        out_path(pp, "_chi_disorder_bootstrap_mn_" + mn_suffix(mn), ".csv"),
                        lsdtt::export_chi_disorder_bootstrap(bootstrap_results)
                    );
                }
            }
            catch (const std::exception& e)
            {
                std::cerr << "Warning: writing disorder sample tables failed for m/n = "
                          << mn << ": " << e.what() << "\n";
            }
        }

        // Build basin-level sweep rows
        for (std::size_t i = 0; i < fit_results.basins.size(); ++i)
        {
            disorder_sweep_row row;
            row.mn = mn;
            row.basin_index = static_cast<std::int32_t>(i);
            row.outlet_flat = fit_results.basins[i].outlet_flat;

            row.n_nodes = static_cast<std::int32_t>(fit_results.basins[i].flats.size());
            row.n_main_stem_nodes = count_main_stem_nodes(fit_results.basins[i]);
            row.n_branches = count_unique_branch_ids(fit_results.basins[i]);
            row.n_side_branches = std::max<std::int32_t>(0, row.n_branches - 1);

            row.full_network_disorder = disorder_results.disorder[i];

            if (bf["run_subset_sensitivity"])
            {
                row.n_subset_sensitivity_samples =
                    static_cast<std::int32_t>(subset_results.basins[i].samples.size());
            }

            if (bf["run_jackknife"])
            {
                row.jackknife_mean = jackknife_results.basins[i].jackknife_mean;
                row.jackknife_standard_error =
                    jackknife_results.basins[i].jackknife_standard_error;
            }

            if (bf["run_bootstrap"])
            {
                row.bootstrap_mean = bootstrap_results.basins[i].bootstrap_mean;
                row.bootstrap_standard_deviation =
                    bootstrap_results.basins[i].bootstrap_standard_deviation;
                row.bootstrap_percentile_2p5 =
                    bootstrap_results.basins[i].percentile_2p5;
                row.bootstrap_percentile_97p5 =
                    bootstrap_results.basins[i].percentile_97p5;
            }

            disorder_rows.push_back(row);
        }
    }

    try
    {
        write_disorder_sweep_csv(
            out_path(pp, "_chi_disorder_sweep", ".csv"),
            disorder_rows
        );
    }
    catch (const std::exception& e)
    {
        std::cerr << "Warning: writing disorder sweep summary failed: "
                  << e.what() << "\n";
    }

    // -------------------------------------------------------------------------
    // Spatial k_sn / u* inversion
    // -------------------------------------------------------------------------
    if (bf["run_inversion"])
    {
        std::cout << "\n--- Running spatial chi k_sn / u* inversion ---\n";

        lsdtt::lsdtt_chi_ksn_inversion_options inv_opts;
        inv_opts.A_0 = static_cast<double>(ff["chi_A0"]);
        inv_opts.coarse_cell_size =
            static_cast<double>(ff["inversion_coarse_cell_size"]);
        inv_opts.roughness_weight =
            static_cast<double>(ff["inversion_roughness_weight"]);
        inv_opts.diagonal_jitter =
            static_cast<double>(ff["inversion_diagonal_jitter"]);
        inv_opts.use_mean_absolute_error_for_model_selection =
            bf["inversion_use_mae_for_model_selection"];

        lsdtt::lsdtt_chi_ksn_inversion_collection_result inversion_results;
        try
        {
            inversion_results = lsdtt::invert_chi_ksn(
                *flow,
                inversion_networks,
                tested_mn,
                inv_opts
            );
        }
        catch (const std::exception& e)
        {
            std::cerr << "FATAL: chi ksn inversion failed: "
                      << e.what() << "\n";
            return 1;
        }

        try
        {
            write_inversion_trial_csv(
                out_path(pp, "_chi_ksn_inversion_trials", ".csv"),
                inversion_results
            );

            write_inversion_best_summary_csv(
                out_path(pp, "_chi_ksn_inversion_best", ".csv"),
                inversion_results
            );
        }
        catch (const std::exception& e)
        {
            std::cerr << "Warning: writing inversion CSV summaries failed: "
                      << e.what() << "\n";
        }

        if (bf["write_inversion_rasters"])
        {
            for (std::size_t i = 0; i < inversion_results.basins.size(); ++i)
            {
                const auto& basin = inversion_results.basins[i];

                try
                {
                    const std::string u_file =
                        out_path(pp,
                                 "_chi_ksn_u_star_basin_" + std::to_string(i),
                                 ".tif");
                    loader.write_full(
                        basin.u_star_field.get_data(),
                        basin.u_star_field.get_meta(),
                        u_file,
                        out_fmt
                    );

                    const std::string s_file =
                        out_path(pp,
                                 "_chi_ksn_support_count_basin_" + std::to_string(i),
                                 ".tif");
                    loader.write_full(
                        basin.support_count_field.get_data(),
                        basin.support_count_field.get_meta(),
                        s_file,
                        out_fmt
                    );
                }
                catch (const std::exception& e)
                {
                    std::cerr << "Warning: writing inversion rasters failed for basin "
                              << i << ": " << e.what() << "\n";
                }
            }
        }
    }

    // -------------------------------------------------------------------------
    // Echo parameters
    // -------------------------------------------------------------------------
    pp.print_parameters();

    std::cout << "\nlsdtt-chi-analysis: done.\n";
    return 0;
}
