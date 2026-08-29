// =============================================================================
// lsdtt_hillslope_cellular_cpu.hpp
//
// Experimental conservative face-transfer hillslope solver. This is the global
// deterministic correctness oracle for plan_cellular_hillslope_transport.md;
// it is not exposed through the LEM or CLI.
// =============================================================================

#ifndef LSDTT_HILLSLOPE_CELLULAR_CPU_HPP
#define LSDTT_HILLSLOPE_CELLULAR_CPU_HPP

#include <cstdint>
#include <string>
#include <vector>

#include "lsdtt_raster_cpu.hpp"

namespace lsdtt
{

struct lsdtt_hillslope_cellular_options
{
    double D = 1e-3;
    double U = 0.0;
    double Sc = 0.8;
    double r_match = 0.85;
    double cfl_safety = 0.9;

    // Optional deterministic nonlocal split of near-critical excess flux.
    // Zero fraction preserves the local solver bit-for-bit.
    double nonlocal_fraction = 0.0;
    double nonlocal_activation_ratio = 0.8;
    double nonlocal_deposition_probability = 0.5;
    int nonlocal_max_hops = 4;
    std::string boundary_conditions = "oooo";
    long max_substeps = 1000000;
    int max_multirate_level = 20;
    int n_threads = 0;
};

struct lsdtt_hillslope_cellular_workspace
{
    std::vector<double> elevation;
    std::vector<double> elevation_new;
    std::vector<double> q_east;
    std::vector<double> q_south;
    std::vector<double> q_nonlocal_east;
    std::vector<double> q_nonlocal_south;

    // Multirate edge schedule, retained to avoid repeated allocations.
    std::vector<std::size_t> edge_a;
    std::vector<std::size_t> edge_b;
    std::vector<int> edge_row;
    std::vector<int> edge_col;
    std::vector<std::uint8_t> edge_direction; // 0 east, 1 south
    std::vector<std::uint8_t> edge_colour;
    std::vector<std::uint8_t> edge_level;
    std::vector<double> edge_transfer;
};

struct lsdtt_hillslope_cellular_result
{
    lsdtt_raster_cpu elevation;

    long n_substeps = 0;
    long long n_face_evaluations = 0;
    long long n_face_events = 0;
    int max_multirate_level = 0;
    long stability_bin_violations = 0;
    double min_substep = 0.0;
    double max_substep = 0.0;

    long initial_active_faces = 0;
    long initial_faces_r70 = 0;
    long initial_faces_r80 = 0;
    long initial_faces_r85 = 0;
    long initial_faces_r100 = 0;
    double initial_max_slope_ratio = 0.0;
    double final_max_slope_ratio = 0.0;

    // Positive means net volume left solved cells through a channel or 'o' cell.
    // On a closed domain this should be zero to round-off.
    double net_dirichlet_export_volume = 0.0;
    long long nonlocal_face_events = 0;
    double nonlocal_routed_volume = 0.0;
    double nonlocal_export_volume = 0.0;
    double nonlocal_mean_hops = 0.0;

    bool completed = false;
    std::string message;
};

// Advance the soft-critical flux law over dt_total. Channel cells and 'o' edge
// cells are fixed. is_channel may be empty (no channel Dirichlet set) or have
// one byte per raster cell.
lsdtt_hillslope_cellular_result compute_hillslope_cellular_transport(
    const lsdtt_raster_cpu& dem,
    const std::vector<std::uint8_t>& is_channel,
    double dt_total,
    const lsdtt_hillslope_cellular_options& options,
    lsdtt_hillslope_cellular_workspace* workspace = nullptr);

// Experimental local-time-step form. Each face is assigned a power-of-two
// update level from its local tangent modulus. Transfers due at one fine tick
// are computed from the same frozen field, then applied in deterministic,
// conflict-free edge colours. U must be zero: uplift remains the separate LEM
// process and is not locally subcycled.
lsdtt_hillslope_cellular_result compute_hillslope_cellular_transport_multirate(
    const lsdtt_raster_cpu& dem,
    const std::vector<std::uint8_t>& is_channel,
    double dt_total,
    const lsdtt_hillslope_cellular_options& options,
    lsdtt_hillslope_cellular_workspace* workspace = nullptr);

} // namespace lsdtt

#endif // LSDTT_HILLSLOPE_CELLULAR_CPU_HPP
