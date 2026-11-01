// lsdtt_perimeter_sources_cpu.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// See lsdtt_perimeter_sources_cpu.hpp.
// -----------------------------------------------------------------------------

#include "lsdtt_perimeter_sources_cpu.hpp"

#include <algorithm>
#include <cmath>
#include <deque>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <unordered_map>

#include "../utils/lsdtt_array2d.hpp"
#include "../utils/lsdtt_binary_array2d_tracing.hpp"

namespace lsdtt
{

std::string perimeter_run_verdict_name(perimeter_run_verdict v)
{
    switch (v)
    {
        case perimeter_run_verdict::too_narrow:          return "too_narrow";
        case perimeter_run_verdict::not_low_spot:        return "not_low_spot";
        case perimeter_run_verdict::outlet_crossing:     return "outlet_crossing";
        case perimeter_run_verdict::failed_inward:       return "failed_inward";
        case perimeter_run_verdict::drains_elsewhere:    return "drains_elsewhere";
        case perimeter_run_verdict::too_close_to_outlet: return "too_close_to_outlet";
        case perimeter_run_verdict::source:              return "source";
    }
    return "unknown";
}

namespace
{

constexpr int k_dr[8] = {-1, -1, 0, 1, 1, 1, 0, -1};
constexpr int k_dc[8] = {0, 1, 1, 1, 0, -1, -1, -1};

// Erode a mask k times: a cell goes if any 8-neighbour is off or off-grid.
std::vector<std::uint8_t> erode(std::vector<std::uint8_t> m, int rows, int cols, int k)
{
    for (int it = 0; it < k; ++it)
    {
        std::vector<std::uint8_t> out = m;
        for (int r = 0; r < rows; ++r)
            for (int c = 0; c < cols; ++c)
            {
                const std::size_t i = static_cast<std::size_t>(r) * cols + c;
                if (!m[i]) continue;
                for (int d = 0; d < 8; ++d)
                {
                    const int rr = r + k_dr[d], cc = c + k_dc[d];
                    if (rr < 0 || rr >= rows || cc < 0 || cc >= cols ||
                        !m[static_cast<std::size_t>(rr) * cols + cc])
                    { out[i] = 0; break; }
                }
            }
        m.swap(out);
    }
    return m;
}

// Connected components of the cells where `on` is set; returns labels (-1
// off) and the size of each component.
std::vector<std::int32_t> label_components(const std::vector<std::uint8_t>& on, int rows,
                                           int cols, int connectivity,
                                           std::vector<std::int64_t>& sizes)
{
    const std::size_t n = on.size();
    std::vector<std::int32_t> lab(n, -1);
    std::vector<std::size_t> stack;
    const int nd = connectivity == 8 ? 8 : 4;
    // 4-connected steps are the even entries of the D8 tables.
    for (std::size_t s = 0; s < n; ++s)
    {
        if (!on[s] || lab[s] >= 0) continue;
        const std::int32_t id = static_cast<std::int32_t>(sizes.size());
        sizes.push_back(0);
        lab[s] = id;
        stack.push_back(s);
        while (!stack.empty())
        {
            const std::size_t i = stack.back();
            stack.pop_back();
            ++sizes[static_cast<std::size_t>(id)];
            const int r = static_cast<int>(i / cols), c = static_cast<int>(i % cols);
            for (int k = 0; k < nd; ++k)
            {
                const int d = nd == 8 ? k : 2 * k;
                const int rr = r + k_dr[d], cc = c + k_dc[d];
                if (rr < 0 || rr >= rows || cc < 0 || cc >= cols) continue;
                const std::size_t j = static_cast<std::size_t>(rr) * cols + cc;
                if (on[j] && lab[j] < 0) { lab[j] = id; stack.push_back(j); }
            }
        }
    }
    return lab;
}

} // namespace

perimeter_sources_result find_perimeter_sources(
    const lsdtt_raster_cpu&          dem,
    const lsdtt_raster_cpu&          routed_dem,
    const lsdtt_flow_base&           flow,
    const perimeter_sources_request& rq)
{
    const lsdtt_raster_meta& meta = dem.get_meta();
    const int rows = meta.n_rows, cols = meta.n_cols;
    const std::size_t n = static_cast<std::size_t>(rows) * cols;
    const double cell = meta.cell_size;
    if (routed_dem.get_data().size() != n || flow.full_grid_size() != n)
        throw std::invalid_argument("find_perimeter_sources: DEM / routed DEM / flow size mismatch.");
    if (rq.n_outlets < 1)
        throw std::invalid_argument("find_perimeter_sources: n_outlets must be >= 1.");
    if (!(rq.min_valley_width_metres > 0.0) || !(rq.max_gradient > 0.0) ||
        !(rq.inward_max_gradient > 0.0) || rq.min_relief_metres < 0.0 ||
        rq.inward_check_distance_metres < 0.0 || rq.trim_pixels < 0)
        throw std::invalid_argument(
            "find_perimeter_sources: widths and gradients must be > 0; relief, check "
            "distance and trim must be >= 0.");

    perimeter_sources_result res;
    perimeter_sources_report& rep = res.report;
    rep.gradient_radius_metres = rq.gradient_radius_metres >= 0.0
        ? rq.gradient_radius_metres : std::max(5.0 * cell, 10.0);
    rep.gap_closing_metres = rq.gap_closing_metres >= 0.0
        ? rq.gap_closing_metres : std::max(2.0 * cell, 20.0);
    rep.inward_gradient_window_metres = rq.inward_gradient_window_metres > 0.0
        ? rq.inward_gradient_window_metres : std::max(5.0 * cell, 25.0);
    rep.min_source_path_metres = rq.min_source_path_metres >= 0.0
        ? rq.min_source_path_metres : 2.0 * rq.min_valley_width_metres;

    const std::vector<float>&        z  = dem.get_data();
    const std::vector<std::uint8_t>& em = dem.get_mask();
    const std::vector<float>&        zr = routed_dem.get_data();

    // -------------------------------------------------------------------------
    // Valid region: trim, keep the largest 8-connected component, find holes.
    // -------------------------------------------------------------------------
    std::vector<std::uint8_t> valid(n);
    for (std::size_t i = 0; i < n; ++i) valid[i] = em[i] ? 1 : 0;
    rep.n_valid_cells = std::count(valid.begin(), valid.end(), std::uint8_t{1});
    if (rep.n_valid_cells == 0)
        throw std::runtime_error("find_perimeter_sources: the DEM has no valid cells.");
    const std::vector<std::uint8_t> trimmed = erode(valid, rows, cols, rq.trim_pixels);

    std::vector<std::int64_t> comp_sizes;
    const std::vector<std::int32_t> comp = label_components(trimmed, rows, cols, 8, comp_sizes);
    if (comp_sizes.empty())
        throw std::runtime_error(
            "find_perimeter_sources: no valid cells left after perimeter_trim_pixels.");
    const std::int32_t main_comp = static_cast<std::int32_t>(
        std::max_element(comp_sizes.begin(), comp_sizes.end()) - comp_sizes.begin());
    rep.n_traced_component_cells = comp_sizes[static_cast<std::size_t>(main_comp)];
    rep.n_components_dropped = static_cast<int>(comp_sizes.size()) - 1;
    rep.n_cells_dropped = std::accumulate(comp_sizes.begin(), comp_sizes.end(),
                                          std::int64_t{0}) - rep.n_traced_component_cells;

    // Background, 4-connected: the parts not touching the grid edge are holes.
    std::vector<std::uint8_t> background(n);
    for (std::size_t i = 0; i < n; ++i) background[i] = comp[i] == main_comp ? 0 : 1;
    std::vector<std::int64_t> bg_sizes;
    const std::vector<std::int32_t> bg = label_components(background, rows, cols, 4, bg_sizes);
    std::vector<std::uint8_t> bg_exterior(bg_sizes.size(), 0);
    for (int r = 0; r < rows; ++r)
        for (int c = 0; c < cols; ++c)
        {
            if (r != 0 && r != rows - 1 && c != 0 && c != cols - 1) continue;
            const std::int32_t b = bg[static_cast<std::size_t>(r) * cols + c];
            if (b >= 0) bg_exterior[static_cast<std::size_t>(b)] = 1;
        }
    for (std::size_t b = 0; b < bg_sizes.size(); ++b)
        if (!bg_exterior[b]) { ++rep.n_holes; rep.n_hole_cells += bg_sizes[b]; }

    // -------------------------------------------------------------------------
    // Trace the outer boundary (the tracer pads by one background cell, so
    // data reaching the grid edge is traced like data reaching nodata).
    // -------------------------------------------------------------------------
    lsdtt_array2d_u8 bw(rows, cols, 0);
    for (std::size_t i = 0; i < n; ++i)
        if (comp[i] == main_comp) bw.data()[i] = 1;
    std::vector<lsdtt_pixel_rc> trace = trace_boundary_pixels(bw, 8);
    // The trace ends by stepping back onto its start.
    if (trace.size() > 1 && trace.back().row == trace.front().row &&
        trace.back().col == trace.front().col)
        trace.pop_back();
    const std::size_t N = trace.size();

    std::vector<perimeter_profile_visit>& prof = res.profile;
    prof.resize(N);
    for (std::size_t k = 0; k < N; ++k)
    {
        prof[k].row = trace[k].row;
        prof[k].col = trace[k].col;
        prof[k].flat = static_cast<std::int64_t>(trace[k].row) * cols + trace[k].col;
        prof[k].elevation = z[static_cast<std::size_t>(prof[k].flat)];
    }
    // Start at the highest cell, which cannot be inside a valley crossing.
    {
        std::size_t top = 0;
        for (std::size_t k = 1; k < N; ++k)
            if (prof[k].elevation > prof[top].elevation) top = k;
        std::rotate(prof.begin(), prof.begin() + static_cast<std::ptrdiff_t>(top), prof.end());
    }
    std::unordered_map<std::int64_t, std::size_t> first_visit;
    for (std::size_t k = 0; k < N; ++k)
    {
        const auto ins = first_visit.emplace(prof[k].flat, k);
        if (!ins.second) { prof[k].repeated_visit = true; ++rep.n_repeated_visits; }
    }
    auto step_metres = [&](std::size_t a, std::size_t b) {
        const int dr = std::abs(prof[a].row - prof[b].row), dc = std::abs(prof[a].col - prof[b].col);
        if (dr == 0 && dc == 0) return 0.0;
        return (dr && dc) ? cell * std::sqrt(2.0) : cell; };
    for (std::size_t k = 1; k < N; ++k)
        prof[k].s_metres = prof[k - 1].s_metres + step_metres(k - 1, k);
    const double L = N > 1 ? prof[N - 1].s_metres + step_metres(N - 1, 0) : cell;
    rep.n_profile_visits = N;
    rep.profile_length_metres = L;
    // Forward distance from visit a to visit b around the circle.
    auto fwd = [&](std::size_t a, std::size_t b) {
        return b >= a ? prof[b].s_metres - prof[a].s_metres
                      : L - prof[a].s_metres + prof[b].s_metres; };

    // -------------------------------------------------------------------------
    // Gradient: least-squares plane through the valid cells within the radius.
    // -------------------------------------------------------------------------
    {
        const double R = rep.gradient_radius_metres;
        const int rc = std::max(1, static_cast<int>(std::ceil(R / cell)));
        const float nanf = std::numeric_limits<float>::quiet_NaN();
        for (std::size_t k = 0; k < N; ++k)
        {
            if (prof[k].repeated_visit)
            { prof[k].gradient = prof[first_visit[prof[k].flat]].gradient; continue; }
            const int r0 = prof[k].row, c0 = prof[k].col;
            const double z0 = prof[k].elevation;
            double sxx = 0, sxy = 0, syy = 0, sx = 0, sy = 0, sz = 0, sxz = 0, syz = 0, m = 0;
            for (int dr = -rc; dr <= rc; ++dr)
                for (int dc = -rc; dc <= rc; ++dc)
                {
                    const int r = r0 + dr, c = c0 + dc;
                    if (r < 0 || r >= rows || c < 0 || c >= cols) continue;
                    const double x = dc * cell, y = -dr * cell;
                    if (x * x + y * y > R * R + 1e-9) continue;
                    const std::size_t i = static_cast<std::size_t>(r) * cols + c;
                    if (!em[i]) continue;
                    const double zz = z[i] - z0;
                    sxx += x * x; sxy += x * y; syy += y * y; sx += x; sy += y;
                    sz += zz; sxz += x * zz; syz += y * zz; m += 1.0;
                }
            // Normal equations for z = a x + b y + c, by Cramer's rule.
            const double det = sxx * (syy * m - sy * sy) - sxy * (sxy * m - sy * sx) +
                               sx * (sxy * sy - syy * sx);
            if (m < 3.0 || std::abs(det) < 1e-9 * std::pow(cell, 4) * m * m * m)
            { prof[k].gradient = nanf; continue; }
            const double da = sxz * (syy * m - sy * sy) - sxy * (syz * m - sy * sz) +
                              sx * (syz * sy - syy * sz);
            const double db = sxx * (syz * m - sz * sy) - sxz * (sxy * m - sy * sx) +
                              sx * (sxy * sz - syz * sx);
            prof[k].gradient = static_cast<float>(std::hypot(da / det, db / det));
        }
        // Cells with no fit take the median of their nearest fitted neighbours.
        std::vector<float> filled(N);
        for (std::size_t k = 0; k < N; ++k)
        {
            filled[k] = prof[k].gradient;
            if (std::isfinite(prof[k].gradient)) continue;
            std::vector<float> nb;
            for (int side = -1; side <= 1; side += 2)
            {
                int got = 0;
                for (std::size_t step = 1; step < N && got < 2; ++step)
                {
                    const std::size_t j = side < 0 ? (k + N - step % N) % N : (k + step) % N;
                    if (std::isfinite(prof[j].gradient)) { nb.push_back(prof[j].gradient); ++got; }
                }
            }
            ++rep.n_gradient_filled;
            if (nb.empty()) { filled[k] = std::numeric_limits<float>::infinity(); continue; }
            std::sort(nb.begin(), nb.end());
            const std::size_t h = nb.size() / 2;
            filled[k] = nb.size() % 2 ? nb[h] : 0.5f * (nb[h - 1] + nb[h]);
        }
        for (std::size_t k = 0; k < N; ++k) prof[k].gradient = filled[k];
    }

    // -------------------------------------------------------------------------
    // Gentle cells and gap closing (circular).
    // -------------------------------------------------------------------------
    std::vector<std::uint8_t> gentle(N);
    for (std::size_t k = 0; k < N; ++k)
        gentle[k] = prof[k].gradient < rq.max_gradient ? 1 : 0;
    const std::size_t n_raw_gentle = std::count(gentle.begin(), gentle.end(), std::uint8_t{1});
    if (n_raw_gentle > 0 && n_raw_gentle < N)
    {
        std::vector<std::uint8_t> closed = gentle;
        // Start just after a gentle cell so no steep run straddles the scan start.
        std::size_t g0 = 0;
        while (!gentle[g0]) ++g0;
        std::size_t k = 1;
        while (k <= N)
        {
            const std::size_t i = (g0 + k) % N;
            if (gentle[i]) { ++k; continue; }
            std::size_t len = 0;
            while (!gentle[(g0 + k + len) % N]) ++len;
            const std::size_t before = (g0 + k - 1) % N, after = (g0 + k + len) % N;
            const double gap = fwd(before, after) - cell;
            if (gap <= rep.gap_closing_metres)
                for (std::size_t q = 0; q < len; ++q) closed[(g0 + k + q) % N] = 1;
            k += len;
        }
        gentle.swap(closed);
    }
    for (std::size_t k = 0; k < N; ++k) prof[k].gentle = gentle[k] != 0;
    rep.n_gentle_visits = std::count(gentle.begin(), gentle.end(), std::uint8_t{1});

    // -------------------------------------------------------------------------
    // Gentle runs, width, lowest cell, low-spot test.
    // -------------------------------------------------------------------------
    std::vector<perimeter_run>& runs = res.runs;
    auto add_run = [&](std::size_t first, std::size_t count) {
        perimeter_run pr;
        pr.id = static_cast<std::int32_t>(runs.size());
        pr.first_visit = first;
        pr.n_visits = count;
        pr.last_visit = (first + count - 1) % N;
        pr.width_metres = count == N ? L : std::min(L, fwd(first, pr.last_visit) + cell);
        double best_mid = std::numeric_limits<double>::infinity();
        pr.min_elevation = std::numeric_limits<float>::infinity();
        for (std::size_t q = 0; q < count; ++q)
        {
            const std::size_t v = (first + q) % N;
            prof[v].run_id = pr.id;
            const double mid = std::abs(static_cast<double>(q) - 0.5 * (count - 1));
            if (prof[v].elevation < pr.min_elevation ||
                (prof[v].elevation == pr.min_elevation && mid < best_mid))
            { pr.min_elevation = prof[v].elevation; pr.lowest_visit = v; best_mid = mid; }
        }
        pr.lowest_flat = prof[pr.lowest_visit].flat;
        runs.push_back(pr);
    };
    if (rep.n_gentle_visits == N && N > 0) add_run(0, N);
    else if (rep.n_gentle_visits > 0)
    {
        std::size_t s0 = 0;
        while (!(gentle[s0] && !gentle[(s0 + N - 1) % N])) ++s0;
        std::size_t k = 0;
        while (k < N)
        {
            const std::size_t i = (s0 + k) % N;
            if (!gentle[i]) { ++k; continue; }
            std::size_t len = 0;
            while (len < N && gentle[(s0 + k + len) % N]) ++len;
            add_run(i, len);
            k += len;
        }
    }
    for (perimeter_run& pr : runs)
    {
        const double W = std::max(pr.width_metres, rq.min_valley_width_metres);
        const std::size_t outside = N - pr.n_visits;
        float zb = -std::numeric_limits<float>::infinity(), za = zb;
        for (std::size_t q = 1; q <= outside; ++q)
        {
            const std::size_t v = (pr.first_visit + N - q % N) % N;
            if (fwd(v, pr.first_visit) > W) break;
            zb = std::max(zb, prof[v].elevation);
        }
        for (std::size_t q = 1; q <= outside; ++q)
        {
            const std::size_t v = (pr.last_visit + q) % N;
            if (fwd(pr.last_visit, v) > W) break;
            za = std::max(za, prof[v].elevation);
        }
        pr.relief_before_metres = std::isfinite(zb) ? zb - pr.min_elevation : 0.0;
        pr.relief_after_metres  = std::isfinite(za) ? za - pr.min_elevation : 0.0;
        pr.relief_metres = std::min(pr.relief_before_metres, pr.relief_after_metres);
        pr.low_spot = pr.relief_metres >= rq.min_relief_metres;
        pr.verdict = pr.width_metres < rq.min_valley_width_metres
            ? perimeter_run_verdict::too_narrow
            : (pr.low_spot ? perimeter_run_verdict::source   // provisional
                           : perimeter_run_verdict::not_low_spot);
    }

    // -------------------------------------------------------------------------
    // Outlet cells: attach each outer D8 outlet (one touching exterior nodata
    // or the grid edge) to its nearest perimeter cell. With trim_pixels = 0
    // most are perimeter cells; a concave-corner cell touching nodata only
    // diagonally is one step away.
    // -------------------------------------------------------------------------
    std::unordered_map<std::int64_t, std::size_t> near_visit;
    {
        std::deque<std::pair<std::int64_t, int>> q;
        for (const auto& [flat, v] : first_visit) { near_visit.emplace(flat, v); q.emplace_back(flat, 0); }
        const int max_depth = rq.trim_pixels + 1;
        while (!q.empty())
        {
            const auto [f, d] = q.front();
            q.pop_front();
            if (d >= max_depth) continue;
            const int r = static_cast<int>(f / cols), c = static_cast<int>(f % cols);
            for (int k = 0; k < 8; ++k)
            {
                const int rr = r + k_dr[k], cc = c + k_dc[k];
                if (rr < 0 || rr >= rows || cc < 0 || cc >= cols) continue;
                const std::int64_t g = static_cast<std::int64_t>(rr) * cols + cc;
                if (!valid[static_cast<std::size_t>(g)] || near_visit.count(g)) continue;
                near_visit.emplace(g, near_visit[f]);
                q.emplace_back(g, d + 1);
            }
        }
    }
    // The run each perimeter visit groups its outlets into: the low-spot runs,
    // widened by the gradient radius at both ends. The plane fit makes the
    // edge of a valley floor read steep for about that distance, and the
    // floor-edge cells are where each side slope's flow arrives, so leaving
    // them out would split the valley's exit.
    std::vector<std::int32_t> grouping_run(N, -1);
    for (const perimeter_run& pr : runs)
    {
        if (!pr.low_spot) continue;
        for (std::size_t q = 0; q < pr.n_visits; ++q)
            grouping_run[(pr.first_visit + q) % N] = pr.id;
        const std::size_t outside = N - pr.n_visits;
        for (int side = -1; side <= 1; side += 2)
            for (std::size_t q = 1; q <= outside; ++q)
            {
                const std::size_t v = side < 0 ? (pr.first_visit + N - q % N) % N
                                               : (pr.last_visit + q) % N;
                const double d = side < 0 ? fwd(v, pr.first_visit) : fwd(pr.last_visit, v);
                if (d > rep.gradient_radius_metres || prof[v].run_id >= 0) break;
                if (grouping_run[v] < 0) grouping_run[v] = pr.id;
            }
    }
    const std::vector<float> area = flow.get_drainage_area_raster().data;
    std::vector<perimeter_outlet_run> oruns;
    std::vector<std::int32_t> run_to_orun(runs.size(), -1);
    for (std::size_t i = 0; i < n; ++i)
    {
        const std::int64_t f = static_cast<std::int64_t>(i);
        if (!em[i] || !flow.is_valid_cell(f) || !flow.is_outlet_unchecked(f)) continue;
        const auto it = near_visit.find(f);
        // An outlet one step from the perimeter but beside a hole, not the
        // outer edge, drains into the hole.
        bool beside_hole = false, beside_exterior = false;
        const int r = static_cast<int>(i / cols), c = static_cast<int>(i % cols);
        for (int k = 0; k < 8; ++k)
        {
            const int rr = r + k_dr[k], cc = c + k_dc[k];
            if (rr < 0 || rr >= rows || cc < 0 || cc >= cols) { beside_exterior = true; continue; }
            const std::int32_t b = bg[static_cast<std::size_t>(rr) * cols + cc];
            if (b < 0 || em[static_cast<std::size_t>(rr) * cols + cc]) continue;
            if (bg_exterior[static_cast<std::size_t>(b)]) beside_exterior = true;
            else beside_hole = true;
        }
        if (it == near_visit.end() || (beside_hole && !beside_exterior))
        {
            if (comp[i] == main_comp) rep.area_draining_to_holes_m2 += area[i];
            else rep.area_draining_elsewhere_m2 += area[i];
            continue;
        }
        ++rep.n_outer_outlet_cells;
        const std::int32_t run = grouping_run[it->second];
        std::int32_t o;
        if (run >= 0)
        {
            o = run_to_orun[static_cast<std::size_t>(run)];
            if (o < 0)
            {
                o = static_cast<std::int32_t>(oruns.size());
                run_to_orun[static_cast<std::size_t>(run)] = o;
                perimeter_outlet_run orr;
                orr.on_crossing = true;
                orr.run_id = run;
                orr.width_metres = runs[static_cast<std::size_t>(run)].width_metres;
                oruns.push_back(orr);
            }
        }
        else
        {
            o = static_cast<std::int32_t>(oruns.size());
            perimeter_outlet_run orr;
            orr.width_metres = cell;
            oruns.push_back(orr);
        }
        perimeter_outlet_run& orr = oruns[static_cast<std::size_t>(o)];
        orr.summed_drainage_area_m2 += area[i];
        ++orr.n_cells;
        orr.cells.push_back(f);
        // The point written for the run is its largest-area cell: where the
        // river is, not wherever a flat floor happens to be lowest.
        if (orr.outlet_flat < 0 || area[i] > area[static_cast<std::size_t>(orr.outlet_flat)])
        { orr.outlet_flat = f; orr.elevation = z[i]; }
    }
    std::stable_sort(oruns.begin(), oruns.end(),
                     [](const perimeter_outlet_run& a, const perimeter_outlet_run& b) {
                         return a.summed_drainage_area_m2 > b.summed_drainage_area_m2; });
    std::unordered_map<std::int64_t, std::int32_t> outlet_of_cell;
    for (std::size_t k = 0; k < oruns.size(); ++k)
    {
        oruns[k].id = static_cast<std::int32_t>(k);
        oruns[k].rank = static_cast<std::int32_t>(k + 1);
        oruns[k].selected = k < static_cast<std::size_t>(rq.n_outlets);
        for (const std::int64_t f : oruns[k].cells)
        {
            outlet_of_cell[f] = oruns[k].id;
            const std::size_t v = near_visit[f];
            prof[v].outlet_cell = true;
            prof[v].outlet_run_id = oruns[k].id;
        }
    }
    rep.n_outlet_runs = oruns.size();
    rep.fewer_outlets_than_requested = oruns.size() < static_cast<std::size_t>(rq.n_outlets);
    {
        const std::size_t nn = static_cast<std::size_t>(rq.n_outlets);
        if (oruns.size() > nn && oruns[nn].summed_drainage_area_m2 > 0.0)
            rep.outlet_area_ratio = oruns[nn - 1].summed_drainage_area_m2 /
                                    oruns[nn].summed_drainage_area_m2;
    }

    // -------------------------------------------------------------------------
    // Inward test and sources.
    // -------------------------------------------------------------------------
    const double W = rep.inward_gradient_window_metres;
    // D8 path from a cell to its outlet, with distance and routed elevation.
    auto trace_path = [&](std::int64_t from, std::vector<std::int64_t>& path,
                          std::vector<double>& dist, std::vector<double>& pz) {
        path.assign(1, from);
        dist.assign(1, 0.0);
        pz.assign(1, zr[static_cast<std::size_t>(from)]);
        std::int64_t cur = from;
        for (std::size_t guard = 0; guard < n && !flow.is_outlet_unchecked(cur); ++guard)
        {
            const std::int64_t nx = flow.receiver_of_unchecked(cur);
            if (nx == cur || nx < 0) break;
            const bool diag = (nx / cols != cur / cols) && (nx % cols != cur % cols);
            dist.push_back(dist.back() + (diag ? cell * std::sqrt(2.0) : cell));
            pz.push_back(zr[static_cast<std::size_t>(nx)]);
            path.push_back(nx);
            cur = nx;
        }
    };
    auto selected_outlet_of = [&](std::int64_t terminal) {
        const auto it = outlet_of_cell.find(terminal);
        return it != outlet_of_cell.end() && oruns[static_cast<std::size_t>(it->second)].selected
            ? it->second : -1; };
    std::vector<std::int64_t> path;
    std::vector<double> dist, pz;
    for (perimeter_run& pr : runs)
    {
        if (pr.verdict != perimeter_run_verdict::source) continue;
        bool holds_selected = false;
        for (const perimeter_outlet_run& orr : oruns)
            if (orr.selected && orr.run_id == pr.id)
            { holds_selected = true; pr.path_outlet_run_id = orr.id; }
        if (holds_selected) { pr.verdict = perimeter_run_verdict::outlet_crossing; continue; }

        // Start cell: the lowest floor cell (within min_relief of the run's
        // lowest) whose flow reaches a selected outlet. Usually the lowest
        // cell itself; but a few centimetres of noise can make an entering
        // channel's edge cell a local low that drains straight out, and a
        // source that drains out at once is no use to a channel tracer.
        std::vector<std::size_t> floor;
        for (std::size_t q = 0; q < pr.n_visits; ++q)
        {
            const std::size_t v = (pr.first_visit + q) % N;
            if (!prof[v].repeated_visit &&
                prof[v].elevation <= pr.min_elevation + rq.min_relief_metres)
                floor.push_back(v);
        }
        std::stable_sort(floor.begin(), floor.end(), [&](std::size_t a, std::size_t b) {
            return prof[a].elevation < prof[b].elevation; });
        std::int64_t start = -1;
        for (const std::size_t v : floor)
        {
            if (flow.is_outlet_unchecked(prof[v].flat)) continue;
            trace_path(prof[v].flat, path, dist, pz);
            if (selected_outlet_of(path.back()) >= 0) { start = prof[v].flat; break; }
        }
        if (start < 0)
        {
            // Nothing on the floor reaches a selected outlet: report the
            // lowest cell's own fate.
            if (flow.is_outlet_unchecked(pr.lowest_flat))
            {
                pr.verdict = perimeter_run_verdict::outlet_crossing;
                const auto it = outlet_of_cell.find(pr.lowest_flat);
                pr.path_outlet_run_id = it != outlet_of_cell.end() ? it->second : -1;
                continue;
            }
            trace_path(pr.lowest_flat, path, dist, pz);
            start = pr.lowest_flat;
        }
        pr.start_flat = start;
        const double total = dist.back();
        pr.inward_checked_metres = std::min(total, rq.inward_check_distance_metres);
        bool failed = false;
        std::size_t j = 0;
        for (std::size_t k = 1; k < path.size() && dist[k] <= rq.inward_check_distance_metres; ++k)
        {
            while (j + 1 < k && dist[k] - dist[j + 1] >= W) ++j;
            double g;
            if (dist[k] - dist[j] >= W) g = (pz[j] - pz[k]) / (dist[k] - dist[j]);
            else if (k + 1 == path.size()) g = (pz[0] - pz[k]) / dist[k];  // short path
            else continue;
            pr.inward_max_gradient = std::max(pr.inward_max_gradient, g);
            if (g > rq.inward_max_gradient)
            {
                failed = true;
                pr.inward_failure_flat = path[k];
                pr.inward_checked_metres = dist[k];
                break;
            }
        }
        const auto it = outlet_of_cell.find(path.back());
        pr.path_outlet_run_id = it != outlet_of_cell.end() ? it->second : -1;
        pr.path_length_to_outlet_metres = total;
        pr.drop_to_outlet_metres = pz.front() - pz.back();
        if (failed) pr.verdict = perimeter_run_verdict::failed_inward;
        else if (pr.path_outlet_run_id < 0 ||
                 !oruns[static_cast<std::size_t>(pr.path_outlet_run_id)].selected)
            pr.verdict = perimeter_run_verdict::drains_elsewhere;
        else if (total < rep.min_source_path_metres)
            pr.verdict = perimeter_run_verdict::too_close_to_outlet;
        else
        {
            perimeter_source ps;
            ps.id = static_cast<std::int32_t>(res.sources.size());
            ps.run_id = pr.id;
            ps.outlet_run_id = pr.path_outlet_run_id;
            ps.flat = start;
            ps.row = static_cast<int>(start / cols);
            ps.col = static_cast<int>(start % cols);
            ps.elevation = z[static_cast<std::size_t>(start)];
            ps.crossing_width_metres = pr.width_metres;
            ps.crossing_relief_metres = pr.relief_metres;
            ps.path_length_to_outlet_metres = total;
            ps.drop_to_outlet_metres = pr.drop_to_outlet_metres;
            res.sources.push_back(ps);
        }
    }
    res.outlet_runs = std::move(oruns);
    return res;
}

} // namespace lsdtt
