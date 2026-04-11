// lsdtt_streamlines_cpu.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html
//
// Kinematic streamline tracing after Lea (1992) and Tarboton (1997).
// Adapted for hilltop-to-channel tracing and erosion rate mapping in:
//   Hurst et al. (2012) JGR-ES doi:10.1029/2011JF002057
//   Grieve et al. (2016) ESurf doi:10.5194/esurf-4-309-2016
// See lsdtt_streamlines_cpu.hpp for full reference details.

#include "lsdtt_streamlines_cpu.hpp"
#include "../utils/lsdtt_thread_utils.hpp"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace lsdtt {

// =============================================================================
// Internal constants
// =============================================================================

static constexpr float  k_nan_f     = std::numeric_limits<float>::quiet_NaN();
static constexpr double k_pi        = M_PI;
static constexpr double k_half_pi   = M_PI * 0.5;
// Epsilon used to nudge sub-pixel positions away from exact cell-wall corners.
// Prevents degenerate tangent values on the next step.
static constexpr double k_wall_eps  = 1.0e-5;

// =============================================================================
// Bearing / angle helpers
// =============================================================================

// Convert compass bearing (degrees; north=0, east=90, clockwise) to the
// east=0 radian convention used in the kinematic routing formulas.
//   north (0°)  → π/2
//   east  (90°) → 0
//   south (180°) → -π/2
//   west  (270°) → ±π
static inline double bearing_to_theta(double degs)
{
    return (90.0 - degs) * (k_pi / 180.0);
}

// =============================================================================
// precompute_aspect
//
// For every interior pixel with a fully-valid 3x3 neighbourhood, compute:
//   bearing_out  — steepest-descent bearing in [0, 360) degrees
//   tan_out      — tan(theta) where theta is the east=0 radian form
//   inv_tan_out  — 1 / tan(theta) (clamped for near-zero tan)
//
// Edge pixels and pixels with any nodata neighbour receive NaN in all three
// arrays.  Flat pixels (zero gradient in both directions) also receive NaN.
// =============================================================================
static void precompute_aspect(
    const float*   dem_data,
    const uint8_t* dem_mask,
    float*         bearing_out,
    float*         tan_out,
    float*         inv_tan_out,
    int            n_rows,
    int            n_cols,
    double         cell_size,
    int            n_threads
)
{
    const int ntot = n_rows * n_cols;

    for (int k = 0; k < ntot; ++k)
    {
        bearing_out[k] = k_nan_f;
        tan_out[k]     = k_nan_f;
        inv_tan_out[k] = k_nan_f;
    }

    #pragma omp parallel for collapse(2) schedule(static) num_threads(n_threads)
    for (int r = 1; r < n_rows - 1; ++r)
    {
        for (int c = 1; c < n_cols - 1; ++c)
        {
            int idx = r * n_cols + c;
            if (!dem_mask[idx]) continue;

            int idx_n = (r-1)*n_cols + c;
            int idx_s = (r+1)*n_cols + c;
            int idx_e = r*n_cols + (c+1);
            int idx_w = r*n_cols + (c-1);

            if (!dem_mask[idx_n] || !dem_mask[idx_s] ||
                !dem_mask[idx_e] || !dem_mask[idx_w]) continue;

            // Central-difference first derivatives
            // p = dz/dx (east positive); q = dz/dy (north positive, row-1 is north)
            double p = (dem_data[idx_e] - dem_data[idx_w]) / (2.0 * cell_size);
            double q = (dem_data[idx_n] - dem_data[idx_s]) / (2.0 * cell_size);

            if (p == 0.0 && q == 0.0) continue;  // flat pixel

            // Steepest-descent bearing: direction of -gradient = (-p, -q) in (east, north)
            // atan2(east_component, north_component) gives bearing from north, clockwise
            double raw = std::atan2(-p, -q) * (180.0 / k_pi);
            if (raw < 0.0) raw += 360.0;

            double theta = bearing_to_theta(raw);
            double t     = std::tan(theta);
            // Clamp 1/tan to avoid division by zero (near-cardinal directions)
            double it = (std::abs(t) > 1.0e-12)
                        ? 1.0 / t
                        : std::copysign(1.0e12, t);

            bearing_out[idx] = static_cast<float>(raw);
            tan_out[idx]     = static_cast<float>(t);
            inv_tan_out[idx] = static_cast<float>(it);
        }
    }
}

// =============================================================================
// trace_one_streamline
//
// Traces a single kinematic streamline from (src_row, src_col) to the nearest
// reachable sink pixel using the Lea (1992) / Tarboton (1997) algorithm.
//
// Cell coordinate convention:
//   xi increases east (toward col+1): 0 = west wall, 1 = east wall
//   yi increases north (toward row-1): 0 = south wall, 1 = north wall
//   (0.5, 0.5) is the cell centre
//
// Direction codes (consistent with old lsdtt convention):
//   dir = 1 — eastward travel  (entered cell from west  face)
//   dir = 2 — southward travel (entered cell from north face)
//   dir = 3 — westward travel  (entered cell from east  face)
//   dir = 4 — northward travel (entered cell from south face)
//
// Projected coordinates of sub-pixel position (r, c, xi, yi):
//   x = meta.x_min + (c + xi) * cell_size
//   y = meta.y_min + (n_rows - r - 1 + yi) * cell_size
// =============================================================================
static lsdtt_streamline trace_one_streamline(
    int                             src_row,
    int                             src_col,
    const float*                    dem_data,
    const uint8_t*                  dem_mask,
    const float*                    sink_data,
    const uint8_t*                  sink_mask_p,
    const float*                    bearing_arr,
    const float*                    tan_arr,
    const float*                    inv_tan_arr,
    const lsdtt_raster_meta&        meta,
    const lsdtt_streamline_options& opts
)
{
    const int    nr      = meta.n_rows;
    const int    nc      = meta.n_cols;
    const double cs      = meta.cell_size;
    const double x_min   = meta.x_min;
    const double y_min   = meta.y_min;

    // Flat-index and projected-coordinate helpers
    auto flat  = [&](int r, int c) -> int { return r * nc + c; };
    auto px_x  = [&](int c, double xi) { return x_min + (c + xi) * cs; };
    auto px_y  = [&](int r, double yi) {
        return y_min + (nr - r - 1 + yi) * cs;
    };

    lsdtt_streamline result;
    result.source_row = src_row;
    result.source_col = src_col;

    // Source pixel must be at least one pixel from every edge so that both the
    // initial step and the 3x3 neighbourhood for the aspect computation are valid.
    if (src_row <= 0 || src_row >= nr-1 || src_col <= 0 || src_col >= nc-1)
    {
        result.outcome = streamline_outcome::hit_edge;
        return result;
    }

    float src_bearing_f = bearing_arr[flat(src_row, src_col)];
    if (std::isnan(src_bearing_f))
    {
        result.outcome = streamline_outcome::no_aspect;
        return result;
    }

    const float src_elev = dem_data[flat(src_row, src_col)];

    // -------------------------------------------------------------------------
    // State
    // -------------------------------------------------------------------------
    int    a     = src_row;   // current cell row
    int    b     = src_col;   // current cell col
    double length = 0.0;      // cumulative path length in cell-widths

    double degs;              // bearing at current cell (degrees)
    double xi, yi;            // entry point into current cell
    int    dir;               // travel direction (1=E 2=S 3=W 4=N)

    std::vector<lsdtt_streamline_node> nodes;

    if (opts.save_nodes)
    {
        lsdtt_streamline_node n0;
        n0.x = px_x(b, 0.5);
        n0.y = px_y(a, 0.5);
        n0.z = src_elev;
        n0.distance_m = 0.0;
        nodes.push_back(n0);
    }

    // -------------------------------------------------------------------------
    // Initial step: exit the source cell from its centre (xi=0.5, yi=0.5)
    //
    // The four quadrant formulas for entry at the cell centre follow directly
    // from the general in-cell routing formulas evaluated at xi=0.5, yi=0.5.
    // -------------------------------------------------------------------------
    {
        degs           = src_bearing_f;
        double theta   = bearing_to_theta(degs);
        double t       = static_cast<double>(tan_arr[flat(a, b)]);
        double it      = static_cast<double>(inv_tan_arr[flat(a, b)]);
        double d;
        double xo, yo;

        if (degs >= 45.0 && degs < 135.0)
        {
            // Easterly
            xo = 1.0;  yo = (1.0 + t) * 0.5;
            d  = std::abs(1.0 / (2.0 * std::cos(theta)));
            xi = 0.0;  yi = yo;
            dir = 1;   ++b;
            if      (yi <= 0.0) yi = k_wall_eps;
            else if (yi >= 1.0) yi = 1.0 - k_wall_eps;
        }
        else if (degs >= 135.0 && degs < 225.0)
        {
            // Southerly
            xo = (1.0 - it) * 0.5;  yo = 0.0;
            d  = std::abs(1.0 / (2.0 * std::cos(k_half_pi - theta)));
            xi = xo;  yi = 1.0;
            dir = 2;   ++a;
            if      (xi <= 0.0) xi = k_wall_eps;
            else if (xi >= 1.0) xi = 1.0 - k_wall_eps;
        }
        else if (degs >= 225.0 && degs < 315.0)
        {
            // Westerly
            xo = 0.0;  yo = (1.0 - t) * 0.5;
            d  = std::abs(1.0 / (2.0 * std::cos(theta)));
            xi = 1.0;  yi = yo;
            dir = 3;   --b;
            if      (yi <= 0.0) yi = k_wall_eps;
            else if (yi >= 1.0) yi = 1.0 - k_wall_eps;
        }
        else  // northerly: >= 315 || < 45
        {
            xo = (1.0 + it) * 0.5;  yo = 1.0;
            d  = std::abs(1.0 / (2.0 * std::cos(k_half_pi - theta)));
            xi = xo;  yi = 0.0;
            dir = 4;   --a;
            if      (xi <= 0.0) xi = k_wall_eps;
            else if (xi >= 1.0) xi = 1.0 - k_wall_eps;
        }

        length += d;

        if (opts.save_nodes && a >= 0 && a < nr && b >= 0 && b < nc)
        {
            lsdtt_streamline_node nd;
            nd.x = px_x(b, xi);
            nd.y = px_y(a, yi);
            nd.z = dem_mask[flat(a,b)] ? dem_data[flat(a,b)] : src_elev;
            nd.distance_m = length * cs;
            nodes.push_back(nd);
        }
    }

    // -------------------------------------------------------------------------
    // Loop-detection state: flat-index → visit count (thread-private)
    // -------------------------------------------------------------------------
    std::unordered_map<int, int> visit_count;

    // -------------------------------------------------------------------------
    // Main trace loop
    // -------------------------------------------------------------------------
    bool flag       = true;
    bool skip_trace = false;
    int  steps      = 1;

    while (flag && a > 0 && a < nr-1 && b > 0 && b < nc-1)
    {
        int fk = flat(a, b);

        // Loop detection: increment visit count; terminate on third visit
        int& vc = visit_count[fk];
        ++vc;
        if (vc >= 3)
        {
            skip_trace = true;
            flag = false;
            break;
        }

        // Sink check
        if (sink_mask_p[fk] && sink_data[fk] > 0.5f)
        {
            result.sink_row = a;
            result.sink_col = b;
            result.outcome  = streamline_outcome::reached_sink;
            result.length_m = length * cs;
            result.relief_m = static_cast<double>(src_elev) - dem_data[fk];
            if (opts.save_nodes)
            {
                lsdtt_streamline_node nd;
                nd.x = px_x(b, 0.5);
                nd.y = px_y(a, 0.5);
                nd.z = dem_data[fk];
                nd.distance_m = length * cs;
                nodes.push_back(nd);
                result.nodes = std::move(nodes);
            }
            return result;
        }

        if (++steps > opts.max_steps) { flag = false; break; }

        float degs_new_f = bearing_arr[fk];
        if (std::isnan(degs_new_f)) { flag = false; break; }

        double degs_new = static_cast<double>(degs_new_f);
        double t_new    = static_cast<double>(tan_arr[fk]);
        double it_new   = static_cast<double>(inv_tan_arr[fk]);
        double d;

        // -----------------------------------------------------------------------
        // Perimeter flow test
        //
        // The flow is "perimeter" (non-reversed) if the new bearing has a
        // component in the same general direction as travel.  Reversed flow
        // (heading back toward the inlet face) triggers edge routing.
        // -----------------------------------------------------------------------
        bool perimeter = false;
        if      (dir == 1 && degs_new >   0.0 && degs_new < 180.0) perimeter = true;
        else if (dir == 2 && degs_new >  90.0 && degs_new < 270.0) perimeter = true;
        else if (dir == 3 && degs_new > 180.0 && degs_new < 360.0) perimeter = true;
        else if (dir == 4 && (degs_new > 270.0 || degs_new < 90.0)) perimeter = true;

        if (perimeter)
        {
            // ------------------------------------------------------------------
            // Normal kinematic routing
            //
            // Test all four possible outlet faces.  The inlet face is excluded
            // by setting its candidate value to -1 (outside [0,1]).
            //
            // temp_yo1: exit east  (xo=1), find yo = yi + (1-xi)*tan
            // temp_xo2: exit south (yo=0), find xo = xi - yi/tan
            // temp_yo2: exit west  (xo=0), find yo = yi - xi*tan
            // temp_xo1: exit north (yo=1), find xo = xi + (1-yi)/tan
            // ------------------------------------------------------------------
            double temp_yo1 = yi + (1.0 - xi) * t_new;
            double temp_xo1 = xi + (1.0 - yi) * it_new;
            double temp_yo2 = yi - xi * t_new;
            double temp_xo2 = xi - yi * it_new;

            // Exclude inlet face
            if      (dir == 1) temp_yo2 = -1.0;  // entered from west → cannot exit west
            else if (dir == 2) temp_xo1 = -1.0;  // entered from north → cannot exit north
            else if (dir == 3) temp_yo1 = -1.0;  // entered from east → cannot exit east
            else if (dir == 4) temp_xo2 = -1.0;  // entered from south → cannot exit south

            if (temp_yo1 > 0.0 && temp_yo1 <= 1.0)
            {
                double xo = 1.0, yo = temp_yo1;
                d  = std::sqrt((xo-xi)*(xo-xi) + (yo-yi)*(yo-yi));
                xi = 0.0;  yi = yo;
                dir = 1;   ++b;
                if      (yi <= 0.0) yi = k_wall_eps;
                else if (yi >= 1.0) yi = 1.0 - k_wall_eps;
            }
            else if (temp_xo2 > 0.0 && temp_xo2 <= 1.0)
            {
                double xo = temp_xo2, yo = 0.0;
                d  = std::sqrt((xo-xi)*(xo-xi) + (yo-yi)*(yo-yi));
                xi = xo;  yi = 1.0;
                dir = 2;   ++a;
                if      (xi <= 0.0) xi = k_wall_eps;
                else if (xi >= 1.0) xi = 1.0 - k_wall_eps;
            }
            else if (temp_yo2 > 0.0 && temp_yo2 <= 1.0)
            {
                double xo = 0.0, yo = temp_yo2;
                d  = std::sqrt((xo-xi)*(xo-xi) + (yo-yi)*(yo-yi));
                xi = 1.0;  yi = yo;
                dir = 3;   --b;
                if      (yi <= 0.0) yi = k_wall_eps;
                else if (yi >= 1.0) yi = 1.0 - k_wall_eps;
            }
            else if (temp_xo1 > 0.0 && temp_xo1 <= 1.0)
            {
                double xo = temp_xo1, yo = 1.0;
                d  = std::sqrt((xo-xi)*(xo-xi) + (yo-yi)*(yo-yi));
                xi = xo;  yi = 0.0;
                dir = 4;   --a;
                if      (xi <= 0.0) xi = k_wall_eps;
                else if (xi >= 1.0) xi = 1.0 - k_wall_eps;
            }
            else
            {
                // No valid outlet face found — numerical edge case with valid
                // aspect data.  Skip rather than loop indefinitely.
                skip_trace = true;
                flag = false;
                break;
            }
        }
        else
        {
            // ------------------------------------------------------------------
            // Edge routing
            //
            // The flow direction has reversed relative to the direction of
            // travel (locally uphill step on rough topography).  Redirect
            // along the cell edge toward the secondary downslope component.
            // ------------------------------------------------------------------
            if (dir == 1)
            {
                if (degs_new <= 90.0 || degs_new >= 270.0)
                {
                    // Secondary component: north
                    double xo = k_wall_eps, yo = 1.0;
                    d  = std::sqrt((xo-xi)*(xo-xi) + (yo-yi)*(yo-yi));
                    xi = xo;  yi = 1.0 - yo;
                    dir = 4;   --a;
                }
                else
                {
                    // Secondary component: south
                    double xo = k_wall_eps, yo = 0.0;
                    d  = std::sqrt((xo-xi)*(xo-xi) + (yo-yi)*(yo-yi));
                    xi = xo;  yi = 1.0 - yo;
                    dir = 2;   ++a;
                }
            }
            else if (dir == 2)
            {
                if (degs_new >= 0.0 && degs_new <= 180.0)
                {
                    // Secondary component: east
                    double xo = 1.0, yo = 1.0 - k_wall_eps;
                    d  = std::sqrt((xo-xi)*(xo-xi) + (yo-yi)*(yo-yi));
                    xi = 1.0 - xo;  yi = yo;
                    dir = 1;   ++b;
                }
                else
                {
                    // Secondary component: west
                    double xo = 0.0, yo = 1.0 - k_wall_eps;
                    d  = std::sqrt((xo-xi)*(xo-xi) + (yo-yi)*(yo-yi));
                    xi = 1.0 - xo;  yi = yo;
                    dir = 3;   --b;
                }
            }
            else if (dir == 3)
            {
                if (degs_new >= 90.0 && degs_new <= 270.0)
                {
                    // Secondary component: south
                    double xo = 1.0 - k_wall_eps, yo = 0.0;
                    d  = std::sqrt((xo-xi)*(xo-xi) + (yo-yi)*(yo-yi));
                    xi = xo;  yi = 1.0 - yo;
                    dir = 2;   ++a;
                }
                else
                {
                    // Secondary component: north
                    double xo = 1.0 - k_wall_eps, yo = 1.0;
                    d  = std::sqrt((xo-xi)*(xo-xi) + (yo-yi)*(yo-yi));
                    xi = xo;  yi = 1.0 - yo;
                    dir = 4;   --a;
                }
            }
            else  // dir == 4
            {
                if (degs_new >= 180.0 && degs_new <= 360.0)
                {
                    // Secondary component: west
                    double xo = 0.0, yo = k_wall_eps;
                    d  = std::sqrt((xo-xi)*(xo-xi) + (yo-yi)*(yo-yi));
                    xi = 1.0 - xo;  yi = yo;
                    dir = 3;   --b;
                }
                else
                {
                    // Secondary component: east
                    double xo = 1.0, yo = k_wall_eps;
                    d  = std::sqrt((xo-xi)*(xo-xi) + (yo-yi)*(yo-yi));
                    xi = 1.0 - xo;  yi = yo;
                    dir = 1;   ++b;
                }
            }
        }

        // Accumulate length only on the first visit to this step
        // (vc is 1 on first visit, 2 on second; we skip double-counting)
        if (vc < 2) length += d;

        degs = degs_new;

        if (opts.save_nodes && a >= 0 && a < nr && b >= 0 && b < nc)
        {
            lsdtt_streamline_node nd;
            nd.x = px_x(b, xi);
            nd.y = px_y(a, yi);
            nd.z = (dem_mask[flat(a,b)]) ? dem_data[flat(a,b)] : src_elev;
            nd.distance_m = length * cs;
            nodes.push_back(nd);
        }
    }  // end while

    // -------------------------------------------------------------------------
    // Classify failed trace
    // -------------------------------------------------------------------------
    if (skip_trace)
    {
        result.outcome = streamline_outcome::loop_detected;
    }
    else if (a <= 0 || a >= nr-1 || b <= 0 || b >= nc-1)
    {
        result.outcome = streamline_outcome::hit_edge;
    }
    else
    {
        // Exceeded max_steps or encountered nodata aspect mid-trace
        result.outcome = streamline_outcome::hit_edge;
    }

    result.length_m = length * cs;

    // Relief against the last-visited cell (may be near-zero on failed traces)
    {
        int ra = std::max(0, std::min(nr-1, a));
        int rb = std::max(0, std::min(nc-1, b));
        int fk = flat(ra, rb);
        result.relief_m = dem_mask[fk]
            ? static_cast<double>(src_elev) - dem_data[fk]
            : 0.0;
    }

    if (opts.save_nodes) result.nodes = std::move(nodes);
    return result;
}

// =============================================================================
// trace_streamlines — public API
// =============================================================================
lsdtt_streamline_result trace_streamlines(
    const lsdtt_raster_cpu&         dem,
    const lsdtt_raster_cpu&         source_mask,
    const lsdtt_raster_cpu&         sink_mask,
    const lsdtt_streamline_options& opts
)
{
    const lsdtt_raster_meta& meta = dem.get_meta();

    // -------------------------------------------------------------------------
    // Validate grid compatibility
    // -------------------------------------------------------------------------
    auto check_compat = [&](const lsdtt_raster_cpu& r, const char* name)
    {
        const auto& m = r.get_meta();
        bool ok = (m.n_rows == meta.n_rows)
               && (m.n_cols == meta.n_cols)
               && (std::abs(m.x_min     - meta.x_min)     < 1.0e-3)
               && (std::abs(m.y_min     - meta.y_min)     < 1.0e-3)
               && (std::abs(m.cell_size - meta.cell_size) < 1.0e-6);
        if (!ok)
        {
            throw std::invalid_argument(
                std::string("trace_streamlines: ") + name +
                " raster is not grid-compatible with the DEM"
            );
        }
    };
    check_compat(source_mask, "source_mask");
    check_compat(sink_mask,   "sink_mask");

    if (opts.thinning < 1)
        throw std::invalid_argument("trace_streamlines: thinning must be >= 1");

    const int nr   = meta.n_rows;
    const int nc   = meta.n_cols;
    const int ntot = nr * nc;

    const int nthreads = resolve_thread_count(opts.n_threads);

    const float*   dem_data   = dem.get_data().data();
    const uint8_t* dem_mask_p = dem.get_mask().data();
    const float*   src_data   = source_mask.get_data().data();
    const uint8_t* src_mask_p = source_mask.get_mask().data();
    const float*   snk_data   = sink_mask.get_data().data();
    const uint8_t* snk_mask_p = sink_mask.get_mask().data();

    // -------------------------------------------------------------------------
    // Precompute aspect arrays (parallelised)
    // -------------------------------------------------------------------------
    std::vector<float> bearing_arr(ntot);
    std::vector<float> tan_arr(ntot);
    std::vector<float> inv_tan_arr(ntot);

    precompute_aspect(
        dem_data, dem_mask_p,
        bearing_arr.data(), tan_arr.data(), inv_tan_arr.data(),
        nr, nc, meta.cell_size, nthreads
    );

    // -------------------------------------------------------------------------
    // Collect source pixels (row-major scan order)
    // -------------------------------------------------------------------------
    std::vector<std::pair<int,int>> sources;
    sources.reserve(static_cast<std::size_t>(ntot) / 10);

    for (int r = 1; r < nr-1; ++r)
    {
        for (int c = 1; c < nc-1; ++c)
        {
            int k = r * nc + c;
            if (!src_mask_p[k])     continue;
            if (src_data[k] <= 0.5f) continue;
            sources.emplace_back(r, c);
        }
    }

    // Apply thinning: retain every Nth source pixel
    if (opts.thinning > 1)
    {
        std::vector<std::pair<int,int>> thinned;
        thinned.reserve(sources.size() / static_cast<std::size_t>(opts.thinning) + 1);
        for (int k = 0; k < static_cast<int>(sources.size()); k += opts.thinning)
            thinned.push_back(sources[k]);
        sources = std::move(thinned);
    }

    // -------------------------------------------------------------------------
    // Trace all streamlines (OMP parallel over source pixels)
    //
    // Pre-allocate output slots so each thread writes to a unique index with
    // no synchronisation needed.  Progress is reported every 5% to stdout
    // using a shared atomic counter; the critical section for printing is
    // entered at most 20 times over the whole run.
    // -------------------------------------------------------------------------
    const int nsrc = static_cast<int>(sources.size());
    std::vector<lsdtt_streamline> raw(nsrc);

    std::cout << "Tracing " << nsrc << " streamlines using "
              << nthreads << " thread" << (nthreads > 1 ? "s" : "") << "...\n"
              << std::flush;

    // Report at 5% intervals (20 messages total); at least every 1 trace.
    const int report_interval = std::max(1, nsrc / 20);

    std::atomic<int> done_count{0};

    #pragma omp parallel for schedule(dynamic, 16) num_threads(nthreads)
    for (int k = 0; k < nsrc; ++k)
    {
        raw[k] = trace_one_streamline(
            sources[k].first, sources[k].second,
            dem_data, dem_mask_p,
            snk_data, snk_mask_p,
            bearing_arr.data(), tan_arr.data(), inv_tan_arr.data(),
            meta, opts
        );

        int completed = ++done_count;

        // Each value of completed is unique (atomic), so the modulo check
        // triggers at most once per interval — no double-printing.
        if (completed % report_interval == 0)
        {
            #pragma omp critical (streamline_progress)
            {
                int pct = (completed * 100) / nsrc;
                std::cout << "\r  " << pct << "% complete  ("
                          << completed << " / " << nsrc << " traces)"
                          << std::flush;
            }
        }
    }

    std::cout << "\r  100% complete  (" << nsrc << " / " << nsrc << " traces)\n"
              << std::flush;

    // -------------------------------------------------------------------------
    // Collect results and tally outcomes
    // -------------------------------------------------------------------------
    lsdtt_streamline_result result;
    result.streamlines.reserve(nsrc);

    for (auto& sl : raw)
    {
        switch (sl.outcome)
        {
            case streamline_outcome::reached_sink:  ++result.n_reached_sink;  break;
            case streamline_outcome::hit_edge:      ++result.n_hit_edge;      break;
            case streamline_outcome::loop_detected: ++result.n_loop_detected; break;
            case streamline_outcome::no_aspect:     ++result.n_no_aspect;     break;
        }
        result.streamlines.push_back(std::move(sl));
    }

    std::cout << "Streamline tracing complete.\n"
              << "  Reached sink:  " << result.n_reached_sink  << "\n"
              << "  Hit edge:      " << result.n_hit_edge       << "\n"
              << "  Loop detected: " << result.n_loop_detected  << "\n"
              << "  No aspect:     " << result.n_no_aspect      << "\n"
              << std::flush;

    return result;
}

// =============================================================================
// write_streamlines
// =============================================================================
void write_streamlines(
    const lsdtt_streamline_result& result,
    const lsdtt_raster_meta&       meta,
    lsdtt_data_loader&             loader,
    const std::string&             out_path,
    const std::string&             layer_name,
    vector_output_format           format,
    bool                           all_outcomes
)
{
    // -------------------------------------------------------------------------
    // GeoPackage warning
    // -------------------------------------------------------------------------
    if (format == vector_output_format::gpkg)
    {
        std::cerr << "WARNING: GeoPackage output selected for streamlines.\n"
                  << "  GeoPackage is backed by SQLite and can fail silently when writing\n"
                  << "  to mounted or network filesystems (e.g. WSL2's /mnt/c).\n"
                  << "  FlatGeobuf (.fgb) is the recommended format and has no such\n"
                  << "  restriction.  If the output file is empty or corrupt, retry\n"
                  << "  with format = flatgeobuf.\n"
                  << std::flush;
    }

    // -------------------------------------------------------------------------
    // Schema
    //
    // outcome integer encoding:
    //   0 = reached_sink
    //   1 = hit_edge
    //   2 = loop_detected
    //   3 = no_aspect
    // -------------------------------------------------------------------------
    lsdtt_vector_schema schema;
    schema.field_names = {
        "source_row", "source_col",
        "sink_row",   "sink_col",
        "length_m",   "relief_m",
        "outcome"
    };
    schema.field_types = {
        lsdtt_vector_field_type::int64_type,   // source_row
        lsdtt_vector_field_type::int64_type,   // source_col
        lsdtt_vector_field_type::int64_type,   // sink_row
        lsdtt_vector_field_type::int64_type,   // sink_col
        lsdtt_vector_field_type::double_type,  // length_m
        lsdtt_vector_field_type::double_type,  // relief_m
        lsdtt_vector_field_type::int64_type    // outcome
    };

    // Helper: projected coordinates of a cell centre
    const double cs    = meta.cell_size;
    const double x_min = meta.x_min;
    const double y_min = meta.y_min;
    const int    nr    = meta.n_rows;
    auto cell_x = [&](int col) { return x_min + (col + 0.5) * cs; };
    auto cell_y = [&](int row) { return y_min + (nr - row - 0.5) * cs; };

    // Helper: outcome to int
    auto outcome_int = [](streamline_outcome o) -> std::int64_t
    {
        switch (o)
        {
            case streamline_outcome::reached_sink:  return 0;
            case streamline_outcome::hit_edge:      return 1;
            case streamline_outcome::loop_detected: return 2;
            case streamline_outcome::no_aspect:     return 3;
        }
        return 3;
    };

    loader.begin_vector_write(
        out_path, layer_name, schema, meta.crs_wkt, format,
        /*spatial_index=*/false
    );

    int n_written = 0;

    for (const auto& sl : result.streamlines)
    {
        // no_aspect traces always skipped — no usable geometry
        if (sl.outcome == streamline_outcome::no_aspect) continue;

        // Filter to reached_sink only unless all_outcomes requested
        if (!all_outcomes && sl.outcome != streamline_outcome::reached_sink) continue;

        // -----------------------------------------------------------------------
        // Build line geometry
        //
        // Prefer the saved sub-pixel node list when available (save_nodes=true).
        // Fall back to a two-vertex source→sink centre line for reached_sink
        // traces when save_nodes=false.  Other outcomes with no nodes are
        // skipped because a single-vertex line is invalid.
        // -----------------------------------------------------------------------
        std::vector<lsdtt_xy> pts;

        if (!sl.nodes.empty())
        {
            pts.reserve(sl.nodes.size());
            for (const auto& nd : sl.nodes)
                pts.push_back({nd.x, nd.y});
        }
        else if (sl.outcome == streamline_outcome::reached_sink)
        {
            pts.push_back({ cell_x(sl.source_col), cell_y(sl.source_row) });
            pts.push_back({ cell_x(sl.sink_col),   cell_y(sl.sink_row)   });
        }
        else
        {
            // Non-sink outcome with no node list — cannot form a valid line
            continue;
        }

        if (pts.size() < 2) continue;

        std::vector<lsdtt_vector_field_value> fv = {
            static_cast<std::int64_t>(sl.source_row),
            static_cast<std::int64_t>(sl.source_col),
            static_cast<std::int64_t>(sl.sink_row),
            static_cast<std::int64_t>(sl.sink_col),
            sl.length_m,
            sl.relief_m,
            outcome_int(sl.outcome)
        };

        loader.write_line(pts, fv);
        ++n_written;
    }

    loader.finish_vector_write();

    std::cout << "Streamlines written: " << n_written
              << " features → " << out_path << "\n"
              << std::flush;
}

} // namespace lsdtt
