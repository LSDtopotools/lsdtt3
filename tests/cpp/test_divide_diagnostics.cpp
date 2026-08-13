// test_divide_diagnostics.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "lsdtt_divide_diagnostics_cpu.hpp"
#include "lsdtt_flowinfo_cpu.hpp"
#include "lsdtt_raster_cpu.hpp"
#include "lsdtt_raster_flow_preprocessing_cpu.hpp"

namespace
{

int g_n_checks = 0;
int g_n_failed = 0;

void check(const std::string& name, bool ok, const std::string& detail = "")
{
    ++g_n_checks;
    if (!ok) ++g_n_failed;

    std::cout << (ok ? "PASS" : "FAIL") << "  " << name;
    if (!ok && !detail.empty()) std::cout << "   [" << detail << "]";
    std::cout << "\n";
}

bool is_finite_non_negative(double v)
{
    return std::isfinite(v) && v >= 0.0;
}

lsdtt::lsdtt_raster_cpu make_south_ramp_dem(int n_rows, int n_cols, double cell_size)
{
    lsdtt::lsdtt_raster_meta meta;
    meta.n_rows = n_rows;
    meta.n_cols = n_cols;
    meta.cell_size = cell_size;
    meta.no_data_value = -9999.0;
    meta.set_default_synthetic_georeference();

    std::vector<float> data(static_cast<std::size_t>(n_rows) * n_cols, 0.0f);
    std::vector<std::uint8_t> mask(data.size(), 1);

    for (int row = 0; row < n_rows; ++row)
    {
        const float z = static_cast<float>(1000.0 - static_cast<double>(row));
        for (int col = 0; col < n_cols; ++col)
        {
            const std::size_t idx =
                static_cast<std::size_t>(row) * static_cast<std::size_t>(n_cols) +
                static_cast<std::size_t>(col);
            data[idx] = z;
        }
    }

    return lsdtt::lsdtt_raster_cpu(std::move(data), std::move(mask), meta);
}

lsdtt::lsdtt_flowinfo make_flow(const lsdtt::lsdtt_raster_cpu& dem)
{
    lsdtt::lsdtt_flowinfo_options fopts;
    fopts.boundary_conditions = "oooo";
    fopts.build_topological_order = true;
    fopts.build_s_vector = true;

    return lsdtt::lsdtt_flowinfo(
        dem.get_data_array(), dem.get_mask_array(), dem.get_meta(), fopts);
}

lsdtt::lsdtt_raster_cpu make_half_step_south_ramp_dem(
    int n_rows,
    int n_cols,
    double cell_size,
    double right_half_step_m)
{
    lsdtt::lsdtt_raster_meta meta;
    meta.n_rows = n_rows;
    meta.n_cols = n_cols;
    meta.cell_size = cell_size;
    meta.no_data_value = -9999.0;
    meta.set_default_synthetic_georeference();

    std::vector<float> data(static_cast<std::size_t>(n_rows) * n_cols, 0.0f);
    std::vector<std::uint8_t> mask(data.size(), 1);

    const int split_col = n_cols / 2;
    const float half_step = static_cast<float>(right_half_step_m);

    for (int row = 0; row < n_rows; ++row)
    {
        for (int col = 0; col < n_cols; ++col)
        {
            const std::size_t idx =
                static_cast<std::size_t>(row) * static_cast<std::size_t>(n_cols) +
                static_cast<std::size_t>(col);

            const float base = static_cast<float>(1000.0 - static_cast<double>(row));
            data[idx] = (col >= split_col) ? (base + half_step) : base;
        }
    }

    return lsdtt::lsdtt_raster_cpu(std::move(data), std::move(mask), meta);
}

void test_invalid_k_throws()
{
    const lsdtt::lsdtt_raster_cpu dem = make_south_ramp_dem(6, 6, 30.0);
    const lsdtt::lsdtt_flowinfo flow = make_flow(dem);

    bool threw = false;
    try
    {
        (void)lsdtt::compute_downstream_ancestor_labels(flow, 0);
    }
    catch (const std::runtime_error&)
    {
        threw = true;
    }

    check("invalid merge_search_steps throws", threw);
}

void test_idempotent_zero_change_ramp()
{
    const int n_rows = 9;
    const int n_cols = 7;
    const lsdtt::lsdtt_raster_cpu dem = make_south_ramp_dem(n_rows, n_cols, 30.0);
    const lsdtt::lsdtt_flowinfo flow = make_flow(dem);

    lsdtt::lsdtt_divide_diagnostic_options opts;
    opts.divide_merge_search_steps = 1;
    opts.use_d8_diagonals = false;

    const lsdtt::lsdtt_divide_step_stats s1 =
        lsdtt::compute_divide_step_stats(flow, dem, opts);
    const lsdtt::lsdtt_divide_step_stats s2 =
        lsdtt::compute_divide_step_stats(flow, dem, opts);

    const std::int64_t expected_pairs =
        static_cast<std::int64_t>(n_rows) * static_cast<std::int64_t>(n_cols - 1);

    check("ramp divide pairs are found",
          s1.n_divide_pairs == expected_pairs,
          "got=" + std::to_string(s1.n_divide_pairs) +
          " expected=" + std::to_string(expected_pairs));

    check("ramp has zero divide step",
          s1.step_median_m == 0.0 && s1.step_mean_m == 0.0 &&
          s1.step_p90_m == 0.0 && s1.step_max_m == 0.0);

    check("ramp stats are finite and non-negative",
          is_finite_non_negative(s1.step_median_m) &&
          is_finite_non_negative(s1.step_mean_m) &&
          is_finite_non_negative(s1.step_p90_m) &&
          is_finite_non_negative(s1.step_max_m) &&
          is_finite_non_negative(s1.gradient_median) &&
          is_finite_non_negative(s1.gradient_p90));

    check("ramp valid-cell count is non-negative",
          s1.n_valid_cells >= 0,
          "got=" + std::to_string(s1.n_valid_cells));

    check("repeated evaluation is idempotent",
          s1.n_divide_pairs == s2.n_divide_pairs &&
          s1.step_median_m == s2.step_median_m &&
          s1.step_mean_m == s2.step_mean_m &&
          s1.step_p90_m == s2.step_p90_m &&
          s1.step_max_m == s2.step_max_m);
}

void test_half_step_and_class_sanity()
{
    const int n_rows = 8;
    const int n_cols = 8;
    const lsdtt::lsdtt_raster_cpu dem =
        make_half_step_south_ramp_dem(n_rows, n_cols, 30.0, 15.0);
    const lsdtt::lsdtt_flowinfo flow = make_flow(dem);

    lsdtt::lsdtt_divide_diagnostic_options opts;
    opts.divide_merge_search_steps = 1;
    opts.use_d8_diagonals = false;
    lsdtt::lsdtt_divide_diagnostic_options opts_diag = opts;
    opts_diag.use_d8_diagonals = true;

    std::vector<float> classes(static_cast<std::size_t>(n_rows) * n_cols, 0.0f);
    for (int row = 0; row < n_rows; ++row)
    {
        for (int col = 0; col < n_cols; ++col)
        {
            const std::size_t idx =
                static_cast<std::size_t>(row) * static_cast<std::size_t>(n_cols) +
                static_cast<std::size_t>(col);
            classes[idx] = ((row + col) % 2 == 0) ? 1.0f : 0.0f;
        }
    }

    const lsdtt::lsdtt_divide_step_stats s =
        lsdtt::compute_divide_step_stats(flow, dem, opts, &classes, 1.0f, 0.0f);
    const lsdtt::lsdtt_divide_step_stats s_diag =
        lsdtt::compute_divide_step_stats(flow, dem, opts_diag, &classes, 1.0f, 0.0f);
    const lsdtt::lsdtt_raster_cpu step_raster =
        lsdtt::make_divide_step_raster(flow, dem, opts);
    const lsdtt::lsdtt_raster_cpu step_raster_diag =
        lsdtt::make_divide_step_raster(flow, dem, opts_diag);

    check("half-step case finds divide pairs",
          s.n_divide_pairs > 0,
          "got=" + std::to_string(s.n_divide_pairs));

    check("half-step step stats are ordered",
          s.step_max_m >= s.step_p90_m &&
          s.step_p90_m >= s.step_median_m &&
          s.step_median_m >= 0.0);

    check("half-step gradient stats are ordered",
          s.gradient_p90 >= s.gradient_median &&
          s.gradient_median >= 0.0);

    check("half-step stats are finite",
          std::isfinite(s.step_median_m) &&
          std::isfinite(s.step_mean_m) &&
          std::isfinite(s.step_p90_m) &&
          std::isfinite(s.step_max_m) &&
          std::isfinite(s.gradient_median) &&
          std::isfinite(s.gradient_p90));

    check("class pairs are present with class raster",
          s.n_class_pairs > 0,
          "got=" + std::to_string(s.n_class_pairs));

    check("class higher fraction is bounded",
          std::isfinite(s.class_higher_fraction) &&
          s.class_higher_fraction >= 0.0 &&
          s.class_higher_fraction <= 1.0,
          "got=" + std::to_string(s.class_higher_fraction));

    check("divide-step raster size/mask match input",
          step_raster.get_data().size() == dem.get_data().size() &&
          step_raster.get_mask() == dem.get_mask());

    check("diagonal divide pairs are not fewer than non-diagonal",
          s_diag.n_divide_pairs >= s.n_divide_pairs,
          "diag=" + std::to_string(s_diag.n_divide_pairs) +
          " non_diag=" + std::to_string(s.n_divide_pairs));

    check("diagonal divide-step raster size/mask match input",
          step_raster_diag.get_data().size() == dem.get_data().size() &&
          step_raster_diag.get_mask() == dem.get_mask());

    std::int64_t n_positive = 0;
    std::int64_t n_positive_diag = 0;
    bool any_negative_valid = false;
    bool any_negative_valid_diag = false;
    double step_sum = 0.0;
    double step_sum_diag = 0.0;
    for (std::size_t i = 0; i < step_raster.get_data().size(); ++i)
    {
        if (step_raster.get_mask()[i] == 0) continue;

        const float v = step_raster.get_data()[i];
        if (v > 0.0f) ++n_positive;
        if (v < 0.0f) any_negative_valid = true;
        step_sum += static_cast<double>(v);

        const float v_diag = step_raster_diag.get_data()[i];
        if (v_diag > 0.0f) ++n_positive_diag;
        if (v_diag < 0.0f) any_negative_valid_diag = true;
        step_sum_diag += static_cast<double>(v_diag);
    }

    check("divide-step raster has positive values on half-step fixture",
          n_positive > 0,
          "got=" + std::to_string(n_positive));

    check("divide-step raster has no negative values on valid pixels",
          !any_negative_valid);

    check("diagonal divide-step raster has no negative values on valid pixels",
          !any_negative_valid_diag);

    check("diagonal-on has at least as much divide-step signal",
          n_positive_diag >= n_positive || step_sum_diag >= step_sum,
          "diag_pos=" + std::to_string(n_positive_diag) +
          " non_diag_pos=" + std::to_string(n_positive) +
          " diag_sum=" + std::to_string(step_sum_diag) +
          " non_diag_sum=" + std::to_string(step_sum));
}

// -----------------------------------------------------------------------------
// Dendritic fixture.
//
// FIXTURE RULE: never test a divide metric on a planar ramp or a tent alone.
// Neither has any flow convergence, so a broken merge test looks perfect on
// them. Only a dendritic network -- one with real junctions where two flow
// paths meet at different distances -- discriminates between the terminal-node
// test and the path-set test.
//
// Noise on a cone, pit-filled and D8-routed, gives exactly that.
// -----------------------------------------------------------------------------
lsdtt::lsdtt_raster_cpu make_dendritic_cone_dem(int n_rows, int n_cols, double cell_size)
{
    lsdtt::lsdtt_raster_meta meta;
    meta.n_rows = n_rows;
    meta.n_cols = n_cols;
    meta.cell_size = cell_size;
    meta.no_data_value = -9999.0;
    meta.set_default_synthetic_georeference();

    std::vector<float> data(static_cast<std::size_t>(n_rows) * n_cols, 0.0f);
    std::vector<std::uint8_t> mask(data.size(), 1);

    // Deterministic hash noise: reproducible on every platform and compiler,
    // unlike <random> distributions.
    auto hash_noise = [](int row, int col) {
        std::uint32_t h = static_cast<std::uint32_t>(row) * 73856093u ^
                          static_cast<std::uint32_t>(col) * 19349663u ^ 0x9e3779b9u;
        h ^= h >> 15; h *= 0x2c1b3c6du;
        h ^= h >> 12; h *= 0x297a2d39u;
        h ^= h >> 15;
        return static_cast<double>(h % 100000u) / 100000.0 - 0.5;  // [-0.5, 0.5)
    };

    for (int row = 0; row < n_rows; ++row)
    {
        for (int col = 0; col < n_cols; ++col)
        {
            const double dr = row - 0.5 * n_rows;
            const double dc = col - 0.5 * n_cols;
            const double radius = std::sqrt(dr * dr + dc * dc);
            const std::size_t idx =
                static_cast<std::size_t>(row) * static_cast<std::size_t>(n_cols) +
                static_cast<std::size_t>(col);
            data[idx] = static_cast<float>(600.0 - 6.0 * radius + 16.0 * hash_noise(row, col));
        }
    }

    const lsdtt::lsdtt_raster_cpu raw(std::move(data), std::move(mask), meta);
    return lsdtt::fill_pits(raw, 0.0001f);
}

// Number of cells that receive flow from two or more donors. A fixture with
// zero junctions cannot discriminate the two divide tests.
std::int64_t count_junctions(const lsdtt::lsdtt_flowinfo& flow)
{
    const std::size_t n = flow.full_grid_size();
    std::vector<int> donors(n, 0);

    for (std::size_t idx = 0; idx < n; ++idx)
    {
        const std::int64_t flat = static_cast<std::int64_t>(idx);
        if (!flow.is_valid_cell(flat)) continue;

        const std::int64_t rec = flow.receiver_of_unchecked(flat);
        if (rec >= 0 && rec != flat && flow.is_valid_cell(rec))
        {
            ++donors[static_cast<std::size_t>(rec)];
        }
    }

    std::int64_t n_junctions = 0;
    for (const int d : donors) if (d >= 2) ++n_junctions;
    return n_junctions;
}

// -----------------------------------------------------------------------------
// Brute-force references. Both walk the pair enumeration independently of the
// library so they can disagree with it.
//
//   reference_path_set_pairs      -- the intended definition: the two forward
//                                    K-step paths share no node.
//   reference_terminal_node_pairs -- the old (defective for K >= 2) test:
//                                    anc_K[a] != anc_K[b]. Retained so the
//                                    K == 1 equivalence can be pinned.
// -----------------------------------------------------------------------------
std::vector<std::int64_t> reference_path(
    const lsdtt::lsdtt_flowinfo& flow, std::int64_t start, int k)
{
    std::vector<std::int64_t> path;
    path.push_back(start);

    std::int64_t cur = start;
    for (int step = 0; step < k; ++step)
    {
        const std::int64_t rec = flow.receiver_of_unchecked(cur);
        if (rec < 0 || rec == cur || !flow.is_valid_cell(rec)) break;
        cur = rec;
        path.push_back(cur);
    }

    return path;
}

bool reference_is_valid(
    std::int64_t flat,
    const lsdtt::lsdtt_flowinfo& flow,
    const lsdtt::lsdtt_raster_cpu& dem)
{
    if (!flow.is_valid_cell(flat)) return false;
    const std::size_t idx = static_cast<std::size_t>(flat);
    if (dem.get_mask()[idx] == 0) return false;
    return std::isfinite(static_cast<double>(dem.get_data()[idx]));
}

// use_path_set = true  -> intended definition
// use_path_set = false -> old terminal-node test
std::int64_t reference_pair_count(
    const lsdtt::lsdtt_flowinfo& flow,
    const lsdtt::lsdtt_raster_cpu& dem,
    int k,
    bool use_d8_diagonals,
    bool use_path_set,
    std::vector<double>* steps_out = nullptr)
{
    const int n_rows = dem.get_meta().n_rows;
    const int n_cols = dem.get_meta().n_cols;

    const std::vector<std::int64_t> anc =
        lsdtt::compute_downstream_ancestor_labels(flow, k);

    std::int64_t n_pairs = 0;

    for (int row = 0; row < n_rows; ++row)
    {
        for (int col = 0; col < n_cols; ++col)
        {
            const std::int64_t a = flow.flat_index_unchecked(row, col);
            if (!reference_is_valid(a, flow, dem)) continue;

            for (int nb = 0; nb < (use_d8_diagonals ? 4 : 2); ++nb)
            {
                int nr = row;
                int nc = col;
                if (nb == 0) nc = col + 1;
                if (nb == 1) nr = row + 1;
                if (nb == 2) { nr = row + 1; nc = col + 1; }
                if (nb == 3) { nr = row + 1; nc = col - 1; }

                if (nr < 0 || nr >= n_rows || nc < 0 || nc >= n_cols) continue;

                const std::int64_t b = flow.flat_index_unchecked(nr, nc);
                if (!reference_is_valid(b, flow, dem)) continue;

                if (flow.receiver_of_unchecked(a) == b ||
                    flow.receiver_of_unchecked(b) == a)
                {
                    continue;
                }

                bool is_divide = false;
                if (use_path_set)
                {
                    const std::vector<std::int64_t> pa = reference_path(flow, a, k);
                    const std::vector<std::int64_t> pb = reference_path(flow, b, k);

                    bool shares = false;
                    for (const std::int64_t va : pa)
                    {
                        if (std::find(pb.begin(), pb.end(), va) != pb.end())
                        {
                            shares = true;
                            break;
                        }
                    }
                    is_divide = !shares;
                }
                else
                {
                    is_divide = anc[static_cast<std::size_t>(a)] !=
                                anc[static_cast<std::size_t>(b)];
                }

                if (!is_divide) continue;

                ++n_pairs;
                if (steps_out != nullptr)
                {
                    const double za =
                        static_cast<double>(dem.get_data()[static_cast<std::size_t>(a)]);
                    const double zb =
                        static_cast<double>(dem.get_data()[static_cast<std::size_t>(b)]);
                    steps_out->push_back(std::fabs(za - zb));
                }
            }
        }
    }

    return n_pairs;
}

void test_dendritic_path_set_agreement()
{
    const int n_rows = 64;
    const int n_cols = 64;
    const lsdtt::lsdtt_raster_cpu dem = make_dendritic_cone_dem(n_rows, n_cols, 30.0);
    const lsdtt::lsdtt_flowinfo flow = make_flow(dem);

    const std::int64_t n_junctions = count_junctions(flow);
    check("dendritic fixture has real junctions",
          n_junctions > 100,
          "junctions=" + std::to_string(n_junctions));

    const std::vector<int> k_values = {1, 2, 3, 5, 10};

    for (const bool diagonals : {false, true})
    {
        std::int64_t previous_pairs = -1;

        for (const int k : k_values)
        {
            lsdtt::lsdtt_divide_diagnostic_options opts;
            opts.divide_merge_search_steps = k;
            opts.use_d8_diagonals = diagonals;

            const lsdtt::lsdtt_divide_step_stats s =
                lsdtt::compute_divide_step_stats(flow, dem, opts);

            std::vector<double> ref_steps;
            const std::int64_t ref_pairs =
                reference_pair_count(flow, dem, k, diagonals, true, &ref_steps);

            const std::string tag =
                std::string(diagonals ? "d8" : "d4") + " K=" + std::to_string(k);

            check("path-set agreement with brute force, " + tag,
                  s.n_divide_pairs == ref_pairs,
                  "impl=" + std::to_string(s.n_divide_pairs) +
                  " ref=" + std::to_string(ref_pairs));

            double ref_sum = 0.0;
            for (const double v : ref_steps) ref_sum += v;
            const double ref_mean =
                ref_steps.empty() ? 0.0 : ref_sum / static_cast<double>(ref_steps.size());

            check("path-set mean step matches brute force, " + tag,
                  std::fabs(s.step_mean_m - ref_mean) < 1.0e-9,
                  "impl=" + std::to_string(s.step_mean_m) +
                  " ref=" + std::to_string(ref_mean));

            if (previous_pairs >= 0)
            {
                check("n_divide_pairs monotone non-increasing in K, " + tag,
                      s.n_divide_pairs <= previous_pairs,
                      "this=" + std::to_string(s.n_divide_pairs) +
                      " previous=" + std::to_string(previous_pairs));
            }
            previous_pairs = s.n_divide_pairs;
        }
    }
}

void test_k1_terminal_node_regression_pin()
{
    // K == 1 is the one case where the terminal-node test and the path-set test
    // agree exactly, so the pre-change behaviour is pinned by comparing against
    // a brute-force implementation of the OLD test.
    const int n_rows = 64;
    const int n_cols = 64;
    const lsdtt::lsdtt_raster_cpu dem = make_dendritic_cone_dem(n_rows, n_cols, 30.0);
    const lsdtt::lsdtt_flowinfo flow = make_flow(dem);

    for (const bool diagonals : {false, true})
    {
        lsdtt::lsdtt_divide_diagnostic_options opts;
        opts.divide_merge_search_steps = 1;
        opts.use_d8_diagonals = diagonals;

        const lsdtt::lsdtt_divide_step_stats s =
            lsdtt::compute_divide_step_stats(flow, dem, opts);

        const std::int64_t old_pairs =
            reference_pair_count(flow, dem, 1, diagonals, false);

        check(std::string("K=1 matches the old terminal-node test, ") +
                  (diagonals ? "d8" : "d4"),
              s.n_divide_pairs == old_pairs,
              "impl=" + std::to_string(s.n_divide_pairs) +
              " old=" + std::to_string(old_pairs));
    }

    // And the same fixture must show the two tests DISAGREEING at K >= 2 --
    // otherwise the fixture is not discriminating and the agreement checks
    // above are vacuous.
    const std::int64_t new_k5 = reference_pair_count(flow, dem, 5, false, true);
    const std::int64_t old_k5 = reference_pair_count(flow, dem, 5, false, false);

    check("fixture discriminates the two tests at K=5",
          old_k5 > new_k5,
          "terminal_node=" + std::to_string(old_k5) +
          " path_set=" + std::to_string(new_k5));
}

void test_thread_count_invariance()
{
    const lsdtt::lsdtt_raster_cpu dem = make_dendritic_cone_dem(48, 48, 30.0);
    const lsdtt::lsdtt_flowinfo flow = make_flow(dem);

    lsdtt::lsdtt_divide_diagnostic_options opts_1;
    opts_1.divide_merge_search_steps = 5;
    opts_1.n_threads = 1;

    lsdtt::lsdtt_divide_diagnostic_options opts_n = opts_1;
    opts_n.n_threads = 4;

    const lsdtt::lsdtt_divide_step_stats s1 =
        lsdtt::compute_divide_step_stats(flow, dem, opts_1);
    const lsdtt::lsdtt_divide_step_stats sn =
        lsdtt::compute_divide_step_stats(flow, dem, opts_n);

    check("stats are independent of n_threads",
          s1.n_divide_pairs == sn.n_divide_pairs &&
          s1.step_median_m == sn.step_median_m &&
          s1.step_mean_m == sn.step_mean_m &&
          s1.step_p90_m == sn.step_p90_m &&
          s1.step_max_m == sn.step_max_m &&
          s1.step_max_flat_a == sn.step_max_flat_a &&
          s1.step_max_flat_b == sn.step_max_flat_b &&
          s1.gradient_median == sn.gradient_median,
          "pairs " + std::to_string(s1.n_divide_pairs) + " vs " +
          std::to_string(sn.n_divide_pairs));

    const lsdtt::lsdtt_raster_cpu r1 = lsdtt::make_divide_step_raster(flow, dem, opts_1);
    const lsdtt::lsdtt_raster_cpu rn = lsdtt::make_divide_step_raster(flow, dem, opts_n);

    check("divide-step raster is independent of n_threads",
          r1.get_data() == rn.get_data() && r1.get_mask() == rn.get_mask());
}

void test_raster_mask_never_holds_nodata()
{
    // A cell the DEM mask accepts but the flow graph rejects used to keep
    // mask == 1 while carrying the no-data value.
    const int n_rows = 32;
    const int n_cols = 32;
    lsdtt::lsdtt_raster_cpu dem = make_dendritic_cone_dem(n_rows, n_cols, 30.0);

    lsdtt::lsdtt_divide_diagnostic_options opts;
    opts.divide_merge_search_steps = 5;

    const lsdtt::lsdtt_flowinfo flow = make_flow(dem);
    const lsdtt::lsdtt_raster_cpu step_raster =
        lsdtt::make_divide_step_raster(flow, dem, opts);

    const float nodata = static_cast<float>(dem.get_meta().no_data_value);
    std::int64_t n_bad = 0;
    for (std::size_t i = 0; i < step_raster.get_data().size(); ++i)
    {
        if (step_raster.get_mask()[i] == 1 && step_raster.get_data()[i] == nodata)
        {
            ++n_bad;
        }
    }

    check("no masked-valid pixel carries the no-data value",
          n_bad == 0,
          "n_bad=" + std::to_string(n_bad));
}

void report_dendritic_pair_counts()
{
    const lsdtt::lsdtt_raster_cpu dem = make_dendritic_cone_dem(64, 64, 30.0);
    const lsdtt::lsdtt_flowinfo flow = make_flow(dem);

    std::cout << "\n  dendritic fixture (64x64, d4 adjacency): "
              << "terminal-node vs path-set pair counts\n";

    for (const int k : {1, 2, 3, 5, 10})
    {
        lsdtt::lsdtt_divide_diagnostic_options opts;
        opts.divide_merge_search_steps = k;
        opts.use_d8_diagonals = false;

        const lsdtt::lsdtt_divide_step_stats s =
            lsdtt::compute_divide_step_stats(flow, dem, opts);
        const std::int64_t old_pairs = reference_pair_count(flow, dem, k, false, false);

        std::cout << "    K=" << k
                  << "  terminal_node=" << old_pairs
                  << "  path_set=" << s.n_divide_pairs
                  << "  gradient_median=" << s.gradient_median << "\n";
    }
    std::cout << "\n";
}

} // namespace

int main()
{
    std::cout << "==============================================\n"
              << " divide diagnostics scaffold\n"
              << "==============================================\n";

    test_invalid_k_throws();
    test_idempotent_zero_change_ramp();
    test_half_step_and_class_sanity();
    test_dendritic_path_set_agreement();
    test_k1_terminal_node_regression_pin();
    test_thread_count_invariance();
    test_raster_mask_never_holds_nodata();
    report_dendritic_pair_counts();

    std::cout << "\nChecks: " << g_n_checks
              << "  Failed: " << g_n_failed << "\n";
    if (g_n_failed == 0)
    {
        std::cout << "RESULT: PASS\n";
        return 0;
    }

    std::cout << "RESULT: FAIL\n";
    return 1;
}
