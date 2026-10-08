// test_streamlines.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
//
// -----------------------------------------------------------------------------
// Phase 0 of the hillslope-channel coupling work: pin down lsdtt_streamlines_cpu.
//
// The module implements Lea (1992) sub-pixel kinematic routing with Tarboton
// (1997) edge routing.  It shipped as "draft" with no tests and no caller, and
// hilltop flow routing rests on it, so this driver exists to
// establish what it actually does before anything is built on top.
//
// THE ANALYTIC CONSTRUCTION
//
// On an exact tilted plane the central-difference aspect is exact, so the
// steepest-descent path is a straight line of known bearing and the trace
// length has a closed form.  Build the plane so that its downslope azimuth is a
// chosen bearing beta:
//
//   z(r,c) = z0 - S*cell*( sin(beta)*c + cos(beta)*(n_rows-1-r) )
//
// gives p = dz/dx = -S*sin(beta) and q = dz/dy = -S*cos(beta) EXACTLY under
// central differences, hence atan2(-p,-q) = beta exactly.
//
// Put the sink on a whole ROW.  A straight path at bearing beta covering a
// north-south distance D has length D/|cos(beta)|, and the trace stops on
// ENTERING the sink cell, so the north-south distance covered is
// (rs - r0 - 0.5) cells: half a cell out of the source centre, then whole cells.
// Therefore
//
//   length_m = (rs - r0 - 0.5) * cell_size / |cos(beta)|
//
// which holds for every bearing with a southward component, diagonals included,
// and does not depend on any internal branch of the implementation.  The
// column analogue with |sin(beta)| covers easterly and westerly planes.
// -----------------------------------------------------------------------------

#include <cmath>
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "lsdtt_raster_cpu.hpp"
#include "lsdtt_streamlines_cpu.hpp"

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

bool rel_close(double got, double want, double rel_tol)
{
    const double denom = std::abs(want) > 1e-12 ? std::abs(want) : 1.0;
    return std::abs(got - want) / denom <= rel_tol;
}

std::string fmt(double v)
{
    return std::to_string(v);
}

constexpr double k_deg2rad = 3.14159265358979323846 / 180.0;

lsdtt::lsdtt_raster_meta make_meta(int n_rows, int n_cols, double cell_size)
{
    lsdtt::lsdtt_raster_meta meta;
    meta.n_rows        = n_rows;
    meta.n_cols        = n_cols;
    meta.cell_size     = cell_size;
    meta.no_data_value = -9999.0;
    meta.set_default_synthetic_georeference();
    return meta;
}

// z(r,c) for a plane whose downslope azimuth is exactly `bearing_deg`.
double plane_z(int r, int c, int n_rows, double cell, double bearing_deg,
               double gradient, double z0)
{
    const double sb = std::sin(bearing_deg * k_deg2rad);
    const double cb = std::cos(bearing_deg * k_deg2rad);
    return z0 - gradient * cell *
                (sb * static_cast<double>(c) +
                 cb * static_cast<double>(n_rows - 1 - r));
}

lsdtt::lsdtt_raster_cpu make_plane(int n_rows, int n_cols, double cell,
                                   double bearing_deg, double gradient = 0.1,
                                   double z0 = 1000.0)
{
    auto meta = make_meta(n_rows, n_cols, cell);
    std::vector<float> data(static_cast<std::size_t>(n_rows) * n_cols, 0.0f);
    std::vector<std::uint8_t> mask(data.size(), 1);

    for (int r = 0; r < n_rows; ++r)
        for (int c = 0; c < n_cols; ++c)
            data[static_cast<std::size_t>(r) * n_cols + c] =
                static_cast<float>(plane_z(r, c, n_rows, cell, bearing_deg,
                                           gradient, z0));

    return lsdtt::lsdtt_raster_cpu(std::move(data), std::move(mask), meta);
}

// A mask raster that is 0 everywhere, with the given cells set to 1.
lsdtt::lsdtt_raster_cpu make_mask_raster(const lsdtt::lsdtt_raster_meta& meta)
{
    std::vector<float> data(
        static_cast<std::size_t>(meta.n_rows) * meta.n_cols, 0.0f);
    std::vector<std::uint8_t> mask(data.size(), 1);
    return lsdtt::lsdtt_raster_cpu(std::move(data), std::move(mask), meta);
}

void set_cell(lsdtt::lsdtt_raster_cpu& r, int row, int col, float v)
{
    r.get_data_array().at(row, col) = v;
}

float get_cell(const lsdtt::lsdtt_raster_cpu& r, int row, int col)
{
    return r.get_data_array().at(row, col);
}

const lsdtt::lsdtt_streamline* find_trace(
    const lsdtt::lsdtt_streamline_result& res, int row, int col)
{
    for (const auto& s : res.streamlines)
        if (s.source_row == row && s.source_col == col) return &s;
    return nullptr;
}

// =============================================================================
// A. Analytic plane: length, relief, linkage
// =============================================================================

// Sink on a whole row; expected length = (rs - r0 - 0.5)*cell / |cos(beta)|.
void test_plane_row_sink(double bearing_deg, const std::string& label,
                         double rel_tol = 2e-4)
{
    const int    nr = 61, nc = 121;
    const double cell = 10.0;
    const double grad = 0.1;
    const double z0   = 1000.0;

    auto dem = make_plane(nr, nc, cell, bearing_deg, grad, z0);

    const double cb = std::cos(bearing_deg * k_deg2rad);
    const double sb = std::sin(bearing_deg * k_deg2rad);

    // A plane descending northward must be traced from a source SOUTH of the
    // sink, and vice versa.
    const bool southward = (cb < 0.0);
    const int  r0 = southward ? 10 : 40;      // source row
    const int  rs = southward ? 40 : 10;      // sink row

    // The straight path drifts east-west as it descends; centre the source so
    // the whole path stays inside the grid rather than exiting the side.
    const double s_cells = std::abs(static_cast<double>(rs - r0)) / std::abs(cb);
    const double drift   = s_cells * sb;
    const int    c0      = static_cast<int>(60.0 - 0.5 * drift);

    auto src = make_mask_raster(dem.get_meta());
    set_cell(src, r0, c0, 1.0f);

    auto snk = make_mask_raster(dem.get_meta());
    for (int c = 0; c < nc; ++c) set_cell(snk, rs, c, 1.0f);

    lsdtt::lsdtt_streamline_options opts;
    opts.n_threads = 1;

    auto res = lsdtt::trace_streamlines(dem, src, snk, opts);

    const auto* t = find_trace(res, r0, c0);
    if (!t)
    {
        check("plane " + label + ": trace produced", false);
        return;
    }

    check("plane " + label + ": reached_sink",
          t->outcome == lsdtt::streamline_outcome::reached_sink,
          "outcome=" + std::to_string(static_cast<int>(t->outcome)));

    if (t->outcome != lsdtt::streamline_outcome::reached_sink) return;

    check("plane " + label + ": sink row is the sink row",
          t->sink_row == rs, "got " + std::to_string(t->sink_row));

    const double want_len =
        (std::abs(static_cast<double>(rs - r0)) - 0.5) * cell /
        std::abs(std::cos(bearing_deg * k_deg2rad));

    check("plane " + label + ": length_m analytic",
          rel_close(t->length_m, want_len, rel_tol),
          "got " + fmt(t->length_m) + " want " + fmt(want_len));

    // relief is measured against the sink cell the trace actually reported
    const double want_relief =
        static_cast<double>(get_cell(dem, r0, c0)) -
        static_cast<double>(get_cell(dem, t->sink_row, t->sink_col));

    check("plane " + label + ": relief_m matches DEM at reported sink",
          rel_close(t->relief_m, want_relief, 1e-5),
          "got " + fmt(t->relief_m) + " want " + fmt(want_relief));

    check("plane " + label + ": relief positive (descending)",
          t->relief_m > 0.0, "got " + fmt(t->relief_m));
}

// Sink on a whole column; expected length = |dc - 0.5|*cell / |sin(beta)|.
void test_plane_col_sink(double bearing_deg, const std::string& label,
                         double rel_tol = 2e-4)
{
    const int    nr = 121, nc = 61;
    const double cell = 10.0;

    auto dem = make_plane(nr, nc, cell, bearing_deg, 0.1, 1000.0);

    const int r0 = 60;
    const int c0 = (bearing_deg < 180.0) ? 10 : 50;   // east-flowing vs west
    const int cs = (bearing_deg < 180.0) ? 40 : 20;

    auto src = make_mask_raster(dem.get_meta());
    set_cell(src, r0, c0, 1.0f);

    auto snk = make_mask_raster(dem.get_meta());
    for (int r = 0; r < nr; ++r) set_cell(snk, r, cs, 1.0f);

    lsdtt::lsdtt_streamline_options opts;
    opts.n_threads = 1;

    auto res = lsdtt::trace_streamlines(dem, src, snk, opts);
    const auto* t = find_trace(res, r0, c0);
    if (!t) { check("plane " + label + ": trace produced", false); return; }

    check("plane " + label + ": reached_sink",
          t->outcome == lsdtt::streamline_outcome::reached_sink,
          "outcome=" + std::to_string(static_cast<int>(t->outcome)));
    if (t->outcome != lsdtt::streamline_outcome::reached_sink) return;

    check("plane " + label + ": sink col is the sink col",
          t->sink_col == cs, "got " + std::to_string(t->sink_col));

    const double want_len =
        (std::abs(static_cast<double>(cs - c0)) - 0.5) * cell /
        std::abs(std::sin(bearing_deg * k_deg2rad));

    check("plane " + label + ": length_m analytic",
          rel_close(t->length_m, want_len, rel_tol),
          "got " + fmt(t->length_m) + " want " + fmt(want_len));
}

// The header claims length is "from source pixel centre to sink pixel centre".
// Measure which of the two candidate definitions the code actually implements.
void test_length_endpoint_definition()
{
    const int    nr = 41, nc = 41;
    const double cell = 10.0;

    auto dem = make_plane(nr, nc, cell, 180.0, 0.1, 1000.0);   // due south

    const int r0 = 10, c0 = 20, rs = 30;

    auto src = make_mask_raster(dem.get_meta());
    set_cell(src, r0, c0, 1.0f);
    auto snk = make_mask_raster(dem.get_meta());
    for (int c = 0; c < nc; ++c) set_cell(snk, rs, c, 1.0f);

    lsdtt::lsdtt_streamline_options opts;
    opts.n_threads = 1;
    auto res = lsdtt::trace_streamlines(dem, src, snk, opts);
    const auto* t = find_trace(res, r0, c0);
    if (!t) { check("length endpoint: trace produced", false); return; }

    const double to_entry  = (static_cast<double>(rs - r0) - 0.5) * cell;
    const double to_centre = static_cast<double>(rs - r0) * cell;

    const bool is_entry  = rel_close(t->length_m, to_entry, 1e-6);
    const bool is_centre = rel_close(t->length_m, to_centre, 1e-6);

    check("length endpoint: is source-centre -> sink-cell-ENTRY (actual)",
          is_entry,
          "got " + fmt(t->length_m) + " entry " + fmt(to_entry) +
          " centre " + fmt(to_centre));

    check("length endpoint: header's 'to sink pixel centre' claim is WRONG "
          "(documentation defect)",
          !is_centre,
          "length matches sink-centre, so the header would be correct and this "
          "test should be deleted");
}

// =============================================================================
// B. Failure modes
// =============================================================================

void test_flat_plateau_no_aspect()
{
    const int nr = 21, nc = 21;
    auto meta = make_meta(nr, nc, 10.0);

    std::vector<float> data(static_cast<std::size_t>(nr) * nc, 500.0f);
    std::vector<std::uint8_t> mask(data.size(), 1);
    lsdtt::lsdtt_raster_cpu dem(std::move(data), std::move(mask), meta);

    auto src = make_mask_raster(meta);
    set_cell(src, 10, 10, 1.0f);
    auto snk = make_mask_raster(meta);
    set_cell(snk, 15, 10, 1.0f);

    lsdtt::lsdtt_streamline_options opts;
    opts.n_threads = 1;
    auto res = lsdtt::trace_streamlines(dem, src, snk, opts);
    const auto* t = find_trace(res, 10, 10);
    if (!t) { check("flat plateau: trace produced", false); return; }

    check("flat plateau: outcome is no_aspect",
          t->outcome == lsdtt::streamline_outcome::no_aspect,
          "outcome=" + std::to_string(static_cast<int>(t->outcome)));
    check("flat plateau: counted in n_no_aspect",
          res.n_no_aspect == 1, "got " + std::to_string(res.n_no_aspect));
}

// FINDING: the public enumeration in trace_streamlines loops r = 1..nr-2 and
// c = 1..nc-2, so a source on the outermost ring is SILENTLY DROPPED - no
// record of any kind is emitted.  The border guard inside trace_one_streamline
// (which would return hit_edge) is therefore unreachable from the public API.
// Pin the real behaviour, because a caller that hands in N hilltop pixels and
// gets back fewer than N traces needs to know why.
void test_source_on_border_is_silently_dropped()
{
    const int nr = 21, nc = 21;
    auto dem = make_plane(nr, nc, 10.0, 180.0);

    auto src = make_mask_raster(dem.get_meta());
    set_cell(src, 0, 10, 1.0f);      // top row: dropped
    set_cell(src, 5, 10, 1.0f);      // interior: traced
    auto snk = make_mask_raster(dem.get_meta());
    for (int c = 0; c < nc; ++c) set_cell(snk, 15, c, 1.0f);

    lsdtt::lsdtt_streamline_options opts;
    opts.n_threads = 1;
    auto res = lsdtt::trace_streamlines(dem, src, snk, opts);

    check("border source: only the interior source is returned "
          "(border sources are silently dropped, NOT reported as hit_edge)",
          res.streamlines.size() == 1,
          "got " + std::to_string(res.streamlines.size()) + " traces");

    check("border source: no record at all for the border pixel",
          find_trace(res, 0, 10) == nullptr);

    check("border source: the interior source is traced normally",
          find_trace(res, 5, 10) != nullptr);
}

void test_no_sink_hits_edge()
{
    const int nr = 31, nc = 31;
    auto dem = make_plane(nr, nc, 10.0, 180.0);

    auto src = make_mask_raster(dem.get_meta());
    set_cell(src, 5, 15, 1.0f);
    auto snk = make_mask_raster(dem.get_meta());   // no sinks at all

    lsdtt::lsdtt_streamline_options opts;
    opts.n_threads = 1;
    auto res = lsdtt::trace_streamlines(dem, src, snk, opts);
    const auto* t = find_trace(res, 5, 15);
    if (!t) { check("no sink: trace produced", false); return; }

    check("no sink: outcome is hit_edge",
          t->outcome == lsdtt::streamline_outcome::hit_edge,
          "outcome=" + std::to_string(static_cast<int>(t->outcome)));
    check("no sink: length is the full run to the edge",
          t->length_m > 0.0, "got " + fmt(t->length_m));
}

void test_nodata_neighbour_no_aspect()
{
    const int nr = 21, nc = 21;
    auto dem = make_plane(nr, nc, 10.0, 180.0);

    // punch a nodata hole immediately north of the source
    dem.get_mask_array().at(9, 10) = 0;

    auto src = make_mask_raster(dem.get_meta());
    set_cell(src, 10, 10, 1.0f);
    auto snk = make_mask_raster(dem.get_meta());
    for (int c = 0; c < nc; ++c) set_cell(snk, 15, c, 1.0f);

    lsdtt::lsdtt_streamline_options opts;
    opts.n_threads = 1;
    auto res = lsdtt::trace_streamlines(dem, src, snk, opts);
    const auto* t = find_trace(res, 10, 10);
    if (!t) { check("nodata neighbour: trace produced", false); return; }

    check("nodata neighbour: outcome is no_aspect",
          t->outcome == lsdtt::streamline_outcome::no_aspect,
          "outcome=" + std::to_string(static_cast<int>(t->outcome)));
}

// The main loop runs while the cell is strictly inside the outer ring, so a
// sink on the outermost ring can never be tested.  Pin that, because a user
// whose channel network touches the DEM edge needs to know.
void test_sink_on_outer_ring_is_unreachable()
{
    const int nr = 21, nc = 21;
    auto dem = make_plane(nr, nc, 10.0, 180.0);

    auto src = make_mask_raster(dem.get_meta());
    set_cell(src, 5, 10, 1.0f);

    auto snk = make_mask_raster(dem.get_meta());
    for (int c = 0; c < nc; ++c) set_cell(snk, nr - 1, c, 1.0f);  // bottom ring

    lsdtt::lsdtt_streamline_options opts;
    opts.n_threads = 1;
    auto res = lsdtt::trace_streamlines(dem, src, snk, opts);
    const auto* t = find_trace(res, 5, 10);
    if (!t) { check("outer-ring sink: trace produced", false); return; }

    check("outer-ring sink: NOT detected, reported hit_edge "
          "(sinks must be >= 1 pixel inside)",
          t->outcome == lsdtt::streamline_outcome::hit_edge,
          "outcome=" + std::to_string(static_cast<int>(t->outcome)));
}

// A closed pit: what does the tracer do?  Assert only that it terminates and
// reports a defined failure outcome; record which one.
void test_closed_pit_terminates()
{
    const int    nr = 31, nc = 31;
    const double cell = 10.0;
    auto meta = make_meta(nr, nc, cell);

    std::vector<float> data(static_cast<std::size_t>(nr) * nc, 0.0f);
    std::vector<std::uint8_t> mask(data.size(), 1);

    // cone-shaped pit centred on (15,15): z rises with radius
    for (int r = 0; r < nr; ++r)
        for (int c = 0; c < nc; ++c)
        {
            const double dr = r - 15.0, dc = c - 15.0;
            data[static_cast<std::size_t>(r) * nc + c] =
                static_cast<float>(100.0 + 2.0 * std::sqrt(dr*dr + dc*dc));
        }

    lsdtt::lsdtt_raster_cpu dem(std::move(data), std::move(mask), meta);

    auto src = make_mask_raster(meta);
    set_cell(src, 5, 15, 1.0f);
    auto snk = make_mask_raster(meta);   // no sink: the pit is the only terminus

    lsdtt::lsdtt_streamline_options opts;
    opts.n_threads = 1;
    opts.max_steps = 10000;

    auto res = lsdtt::trace_streamlines(dem, src, snk, opts);
    const auto* t = find_trace(res, 5, 15);
    if (!t) { check("closed pit: trace produced", false); return; }

    const bool defined =
        t->outcome == lsdtt::streamline_outcome::loop_detected ||
        t->outcome == lsdtt::streamline_outcome::stalled_on_flat ||
        t->outcome == lsdtt::streamline_outcome::hit_edge ||
        t->outcome == lsdtt::streamline_outcome::no_aspect;

    check("closed pit: terminates with a defined failure outcome", defined,
          "outcome=" + std::to_string(static_cast<int>(t->outcome)));

    std::cout << "      [info] closed-pit outcome = "
              << static_cast<int>(t->outcome)
              << " (0=reached_sink 1=hit_edge 2=loop_detected 3=no_aspect)"
              << ", length_m = " << t->length_m << "\n";
}

// =============================================================================
// C. Determinism
// =============================================================================

void test_thread_independence()
{
    const int    nr = 61, nc = 121;
    const double cell = 10.0;
    auto dem = make_plane(nr, nc, cell, 160.0, 0.1, 1000.0);

    auto src = make_mask_raster(dem.get_meta());
    for (int r = 5; r < 15; ++r)
        for (int c = 10; c < 50; ++c)
            set_cell(src, r, c, 1.0f);

    auto snk = make_mask_raster(dem.get_meta());
    for (int c = 0; c < nc; ++c) set_cell(snk, 45, c, 1.0f);

    lsdtt::lsdtt_streamline_result ref;
    bool all_match = true;
    std::string detail;

    for (int nt : {1, 2, 4, 8})
    {
        lsdtt::lsdtt_streamline_options opts;
        opts.n_threads = nt;
        auto res = lsdtt::trace_streamlines(dem, src, snk, opts);

        if (nt == 1) { ref = res; continue; }

        if (res.streamlines.size() != ref.streamlines.size())
        {
            all_match = false;
            detail = "size differs at nt=" + std::to_string(nt);
            break;
        }
        for (std::size_t i = 0; i < res.streamlines.size(); ++i)
        {
            const auto& a = ref.streamlines[i];
            const auto& b = res.streamlines[i];
            if (a.source_row != b.source_row || a.source_col != b.source_col ||
                a.sink_row != b.sink_row || a.sink_col != b.sink_col ||
                a.outcome != b.outcome ||
                !rel_close(a.length_m, b.length_m, 1e-12))
            {
                all_match = false;
                detail = "trace " + std::to_string(i) + " differs at nt=" +
                         std::to_string(nt);
                break;
            }
        }
        if (!all_match) break;
    }

    check("thread independence: identical results at 1/2/4/8 threads",
          all_match, detail);
    check("thread independence: outcome counts sum to the trace count",
          ref.n_reached_sink + ref.n_hit_edge + ref.n_loop_detected +
          ref.n_no_aspect == static_cast<int>(ref.streamlines.size()),
          std::to_string(ref.streamlines.size()));

    check("thread independence: the bulk of traces reached the sink",
          ref.n_reached_sink > 0.9 * static_cast<double>(ref.streamlines.size()),
          std::to_string(ref.n_reached_sink) + " of " +
          std::to_string(ref.streamlines.size()));
}

// =============================================================================
// D. Options and contracts
// =============================================================================

void test_save_nodes_geometry()
{
    const int    nr = 41, nc = 41;
    const double cell = 10.0;
    auto dem = make_plane(nr, nc, cell, 180.0, 0.1, 1000.0);

    const int r0 = 10, c0 = 20, rs = 30;

    auto src = make_mask_raster(dem.get_meta());
    set_cell(src, r0, c0, 1.0f);
    auto snk = make_mask_raster(dem.get_meta());
    for (int c = 0; c < nc; ++c) set_cell(snk, rs, c, 1.0f);

    lsdtt::lsdtt_streamline_options opts;
    opts.n_threads  = 1;
    opts.save_nodes = true;

    auto res = lsdtt::trace_streamlines(dem, src, snk, opts);
    const auto* t = find_trace(res, r0, c0);
    if (!t) { check("save_nodes: trace produced", false); return; }

    check("save_nodes: nodes populated", !t->nodes.empty(),
          "n=" + std::to_string(t->nodes.size()));
    if (t->nodes.empty()) return;

    bool monotone = true;
    for (std::size_t i = 1; i < t->nodes.size(); ++i)
        if (t->nodes[i].distance_m < t->nodes[i-1].distance_m - 1e-9)
            { monotone = false; break; }

    check("save_nodes: distance_m non-decreasing along the path", monotone);

    check("save_nodes: first node distance is zero",
          std::abs(t->nodes.front().distance_m) < 1e-9,
          fmt(t->nodes.front().distance_m));

    check("save_nodes: last node distance equals length_m",
          rel_close(t->nodes.back().distance_m, t->length_m, 1e-9),
          fmt(t->nodes.back().distance_m) + " vs " + fmt(t->length_m));

    // first node must sit at the source cell centre
    const auto& meta = dem.get_meta();
    const double want_x = meta.x_min + (c0 + 0.5) * cell;
    const double want_y = meta.y_min + (nr - r0 - 1 + 0.5) * cell;

    check("save_nodes: first node at source cell centre",
          rel_close(t->nodes.front().x, want_x, 1e-9) &&
          rel_close(t->nodes.front().y, want_y, 1e-9),
          "got (" + fmt(t->nodes.front().x) + "," + fmt(t->nodes.front().y) +
          ") want (" + fmt(want_x) + "," + fmt(want_y) + ")");

    check("save_nodes: first node elevation is the source elevation",
          rel_close(t->nodes.front().z,
                    static_cast<double>(get_cell(dem, r0, c0)), 1e-6));

    // save_nodes must not change the measured quantities
    lsdtt::lsdtt_streamline_options bare;
    bare.n_threads = 1;
    auto res_bare = lsdtt::trace_streamlines(dem, src, snk, bare);
    const auto* tb = find_trace(res_bare, r0, c0);
    check("save_nodes: does not change length_m",
          tb && rel_close(t->length_m, tb->length_m, 1e-12));
}

void test_thinning()
{
    const int nr = 41, nc = 41;
    auto dem = make_plane(nr, nc, 10.0, 180.0);

    auto src = make_mask_raster(dem.get_meta());
    int n_src = 0;
    for (int c = 10; c < 30; ++c) { set_cell(src, 10, c, 1.0f); ++n_src; }

    auto snk = make_mask_raster(dem.get_meta());
    for (int c = 0; c < nc; ++c) set_cell(snk, 30, c, 1.0f);

    lsdtt::lsdtt_streamline_options o1;
    o1.n_threads = 1;
    auto full = lsdtt::trace_streamlines(dem, src, snk, o1);

    lsdtt::lsdtt_streamline_options o2;
    o2.n_threads = 1;
    o2.thinning  = 4;
    auto thin = lsdtt::trace_streamlines(dem, src, snk, o2);

    check("thinning: full run traces every source",
          static_cast<int>(full.streamlines.size()) == n_src,
          std::to_string(full.streamlines.size()) + " of " +
          std::to_string(n_src));

    check("thinning: thinning=4 traces about a quarter",
          thin.streamlines.size() < full.streamlines.size() &&
          thin.streamlines.size() >= static_cast<std::size_t>(n_src / 4 - 1) &&
          thin.streamlines.size() <= static_cast<std::size_t>(n_src / 4 + 1),
          std::to_string(thin.streamlines.size()));
}

void test_grid_mismatch_throws()
{
    auto dem = make_plane(21, 21, 10.0, 180.0);
    auto src = make_mask_raster(make_meta(21, 21, 10.0));
    auto bad = make_mask_raster(make_meta(20, 21, 10.0));   // wrong n_rows

    bool threw = false;
    try
    {
        lsdtt::lsdtt_streamline_options opts;
        opts.n_threads = 1;
        (void)lsdtt::trace_streamlines(dem, src, bad, opts);
    }
    catch (const std::invalid_argument&) { threw = true; }
    catch (const std::exception&)        { threw = true; }

    check("grid mismatch: trace_streamlines throws", threw);
}

// =============================================================================
// E. Cone: radial descent, approximate (central differences are not exact here)
// =============================================================================

void test_cone_radial()
{
    const int    nr = 61, nc = 61;
    const double cell = 10.0;
    const double apex_r = 30.0, apex_c = 30.0;
    auto meta = make_meta(nr, nc, cell);

    std::vector<float> data(static_cast<std::size_t>(nr) * nc, 0.0f);
    std::vector<std::uint8_t> mask(data.size(), 1);

    for (int r = 0; r < nr; ++r)
        for (int c = 0; c < nc; ++c)
        {
            const double dr = r - apex_r, dc = c - apex_c;
            data[static_cast<std::size_t>(r) * nc + c] =
                static_cast<float>(1000.0 - 0.2 * cell * std::sqrt(dr*dr + dc*dc));
        }

    lsdtt::lsdtt_raster_cpu dem(std::move(data), std::move(mask), meta);

    // sink ring: every cell at radius >= 20 pixels
    auto snk = make_mask_raster(meta);
    for (int r = 1; r < nr - 1; ++r)
        for (int c = 1; c < nc - 1; ++c)
        {
            const double dr = r - apex_r, dc = c - apex_c;
            if (std::sqrt(dr*dr + dc*dc) >= 20.0) set_cell(snk, r, c, 1.0f);
        }

    // sources on a ring at radius 5
    auto src = make_mask_raster(meta);
    int n_src = 0;
    for (int r = 1; r < nr - 1; ++r)
        for (int c = 1; c < nc - 1; ++c)
        {
            const double dr = r - apex_r, dc = c - apex_c;
            const double rad = std::sqrt(dr*dr + dc*dc);
            if (rad >= 4.75 && rad <= 5.25) { set_cell(src, r, c, 1.0f); ++n_src; }
        }

    lsdtt::lsdtt_streamline_options opts;
    opts.n_threads = 1;
    auto res = lsdtt::trace_streamlines(dem, src, snk, opts);

    check("cone: sources found", n_src > 0, std::to_string(n_src));
    check("cone: all traces reached the sink ring",
          res.n_reached_sink == static_cast<int>(res.streamlines.size()) &&
          !res.streamlines.empty(),
          std::to_string(res.n_reached_sink) + " of " +
          std::to_string(res.streamlines.size()));

    // A radial path from radius 5 to the ring entry has length ~ (20 - 5 - 0.5)
    // cells.  Central differences are not exact on a cone, so allow 5 %.
    double worst = 0.0;
    for (const auto& t : res.streamlines)
    {
        if (t.outcome != lsdtt::streamline_outcome::reached_sink) continue;
        const double dr = t.source_row - apex_r, dc = t.source_col - apex_c;
        const double r_src = std::sqrt(dr*dr + dc*dc);
        const double dr2 = t.sink_row - apex_r, dc2 = t.sink_col - apex_c;
        const double r_snk = std::sqrt(dr2*dr2 + dc2*dc2);
        const double want  = (r_snk - r_src) * cell;
        const double err   = std::abs(t.length_m - want) / want;
        if (err > worst) worst = err;
    }

    check("cone: path length within 5 % of the radial distance",
          worst < 0.05, "worst relative error " + fmt(worst));
    std::cout << "      [info] cone worst relative length error = "
              << worst << "\n";
}


// -----------------------------------------------------------------------------
// The loop / stall split.
//
// A revisited cell was one outcome until 2026-08-25, and it hid a data problem
// inside what read as a tracer bug. On Indian Creek 35% of hilltop traces came
// back `loop_detected`; the elevation profiles showed every one descending
// monotonically ridge to valley floor, and they ended a median of 116 m from
// the nearest channel while successful traces ended ON one every time. The path
// was right -- there was simply no sink where it arrived. Splitting the outcome
// moved 74% of those into `stalled_on_flat` and cut the genuine loop rate from
// 35% to 9%.
//
// The two are separated by the ground the trace stops on:
//   - stops where it could still descend -> loop_detected   (routing failed)
//   - stops where it could not           -> stalled_on_flat (no sink there)
//
// FIXTURE. A hillslope draining into a closed bowl whose centre is off the
// column axis. Two details are load-bearing:
//   - a bowl, not a flat plain: a dead-flat surface has no aspect at all, so a
//     trace breaks on nodata-aspect before it can ever revisit a cell, which is
//     a different failure;
//   - off-centre: a bowl centred on a symmetry axis is descended straight down
//     that axis and never revisits anything.
// Only the bowl's GRADIENT changes between the two cases, so nothing but the
// classifier can be responsible for the different verdicts.
// -----------------------------------------------------------------------------
lsdtt::lsdtt_raster_cpu bowl_dem(const lsdtt::lsdtt_raster_meta& meta,
                                 double bowl_gradient)
{
    const int    nr = meta.n_rows, nc = meta.n_cols;
    const double cell = meta.cell_size;

    std::vector<float>        data(static_cast<std::size_t>(nr) * nc);
    std::vector<std::uint8_t> mask(data.size(), 1);

    for (int r = 0; r < nr; ++r)
        for (int c = 0; c < nc; ++c)
        {
            const double dr = r - 28.0, dc = c - 20.5;   // 20.5: off-axis
            data[static_cast<std::size_t>(r) * nc + c] =
                (r < 20) ? static_cast<float>(120.0 - 2.0 * r)
                         : static_cast<float>(80.0 + bowl_gradient * cell *
                                              std::sqrt(dr * dr + dc * dc));
        }

    return lsdtt::lsdtt_raster_cpu(std::move(data), std::move(mask), meta);
}

lsdtt::lsdtt_streamline_result trace_into_bowl(double bowl_gradient)
{
    auto meta = make_meta(41, 41, 10.0);
    auto dem  = bowl_dem(meta, bowl_gradient);

    auto src = make_mask_raster(meta);
    set_cell(src, 3, 20, 1.0f);
    auto snk = make_mask_raster(meta);   // deliberately empty: no sink anywhere

    lsdtt::lsdtt_streamline_options opts;
    opts.n_threads = 1;
    opts.max_steps = 50000;
    return lsdtt::trace_streamlines(dem, src, snk, opts);
}

void test_revisit_on_flat_ground_is_a_stall()
{
    auto res = trace_into_bowl(0.004);
    const auto* t = find_trace(res, 3, 20);
    if (!t) { check("stall: trace produced", false); return; }

    check("a revisit on ground too flat to leave is stalled_on_flat",
          t->outcome == lsdtt::streamline_outcome::stalled_on_flat,
          "outcome=" + std::to_string(static_cast<int>(t->outcome)));
    check("and is counted apart from loops",
          res.n_stalled_flat == 1 && res.n_loop_detected == 0,
          "stalled=" + std::to_string(res.n_stalled_flat) +
          " looped=" + std::to_string(res.n_loop_detected));

    // The whole point of the split: the trace DESCENDED. A stall must never be
    // read as a path that went nowhere.
    check("a stalled trace still descended the hillslope", t->relief_m > 30.0,
          "relief=" + std::to_string(t->relief_m));
}

void test_revisit_on_descending_ground_is_still_a_loop()
{
    // Same landscape, same source, ONLY the bowl gradient raised above the
    // threshold. If this also came back as a stall, the classifier would be
    // answering something other than the ground.
    auto res = trace_into_bowl(0.20);
    const auto* t = find_trace(res, 3, 20);
    if (!t) { check("loop: trace produced", false); return; }

    check("a revisit where the trace could still descend is loop_detected",
          t->outcome == lsdtt::streamline_outcome::loop_detected,
          "outcome=" + std::to_string(static_cast<int>(t->outcome)));
    check("and no stall is recorded", res.n_stalled_flat == 0,
          "stalled=" + std::to_string(res.n_stalled_flat));
}

void test_the_split_turns_on_the_threshold()
{
    // The gentle bowl, reclassified purely by moving the threshold below its
    // gradient. Nothing about the DEM changes.
    auto meta = make_meta(41, 41, 10.0);
    auto dem  = bowl_dem(meta, 0.004);
    auto src  = make_mask_raster(meta);
    set_cell(src, 3, 20, 1.0f);
    auto snk  = make_mask_raster(meta);

    lsdtt::lsdtt_streamline_options opts;
    opts.n_threads = 1;
    opts.max_steps = 50000;
    opts.flat_terminus_gradient = 0.0;   // "no ground is ever flat"

    auto res = lsdtt::trace_streamlines(dem, src, snk, opts);
    const auto* t = find_trace(res, 3, 20);
    if (!t) { check("threshold: trace produced", false); return; }

    check("with flat_terminus_gradient = 0 the same trace is a loop",
          t->outcome == lsdtt::streamline_outcome::loop_detected,
          "outcome=" + std::to_string(static_cast<int>(t->outcome)));
}

void test_a_reachable_sink_is_unaffected_by_the_split()
{
    // The split must not cost a successful trace.
    auto meta = make_meta(41, 41, 10.0);

    std::vector<float>        data(static_cast<std::size_t>(41) * 41);
    std::vector<std::uint8_t> mask(data.size(), 1);
    for (int r = 0; r < 41; ++r)
        for (int c = 0; c < 41; ++c)
            data[static_cast<std::size_t>(r) * 41 + c] =
                static_cast<float>(100.0 - 2.0 * r);

    lsdtt::lsdtt_raster_cpu dem(std::move(data), std::move(mask), meta);
    auto src = make_mask_raster(meta);
    set_cell(src, 5, 20, 1.0f);
    auto snk = make_mask_raster(meta);
    for (int c = 0; c < 41; ++c) set_cell(snk, 30, c, 1.0f);

    lsdtt::lsdtt_streamline_options opts;
    opts.n_threads = 1;

    auto res = lsdtt::trace_streamlines(dem, src, snk, opts);
    const auto* t = find_trace(res, 5, 20);
    if (!t) { check("reachable sink: trace produced", false); return; }

    check("a reachable sink still reports reached_sink",
          t->outcome == lsdtt::streamline_outcome::reached_sink,
          "outcome=" + std::to_string(static_cast<int>(t->outcome)));
    check("and no stall is recorded", res.n_stalled_flat == 0,
          "stalled=" + std::to_string(res.n_stalled_flat));
}

} // namespace

int main()
{
    std::cout << "=== lsdtt_streamlines_cpu: Phase 0 characterisation ===\n\n";

    std::cout << "-- A. analytic plane -----------------------------------\n";
    test_plane_row_sink(180.0, "due south");
    test_plane_row_sink(0.0,   "due north (source south of sink)");
    test_plane_row_sink(135.0, "SE diagonal");
    test_plane_row_sink(225.0, "SW diagonal");
    test_plane_row_sink(160.0, "S-SE oblique");
    test_plane_row_sink(200.0, "S-SW oblique");
    test_plane_col_sink(90.0,  "due east");
    test_plane_col_sink(270.0, "due west");

    std::cout << "\n-- A2. length endpoint definition ----------------------\n";
    test_length_endpoint_definition();

    std::cout << "\n-- B. failure modes ------------------------------------\n";
    test_flat_plateau_no_aspect();
    test_source_on_border_is_silently_dropped();
    test_no_sink_hits_edge();
    test_nodata_neighbour_no_aspect();
    test_sink_on_outer_ring_is_unreachable();
    test_closed_pit_terminates();

    std::cout << "\n-- B2. loop vs stall -----------------------------------\n";
    test_revisit_on_flat_ground_is_a_stall();
    test_revisit_on_descending_ground_is_still_a_loop();
    test_the_split_turns_on_the_threshold();
    test_a_reachable_sink_is_unaffected_by_the_split();

    std::cout << "\n-- C. determinism --------------------------------------\n";
    test_thread_independence();

    std::cout << "\n-- D. options and contracts ----------------------------\n";
    test_save_nodes_geometry();
    test_thinning();
    test_grid_mismatch_throws();

    std::cout << "\n-- E. cone ---------------------------------------------\n";
    test_cone_radial();

    std::cout << "\n=======================================================\n";
    std::cout << "checks: " << g_n_checks
              << "   failed: " << g_n_failed << "\n";

    return g_n_failed == 0 ? 0 : 1;
}
