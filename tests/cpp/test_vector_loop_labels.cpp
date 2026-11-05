// test_vector_loop_labels.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
//
// -----------------------------------------------------------------------------
// Loop labels of lsdtt_vector_to_junction_network (river-lines plan SB2):
// loop_reason and closes_cycle on the labelled edges.
//
// A 50 x 30 DEM, 10 m cells, falling 1 m per cell to the east, and two
// separate networks built on it. Expected labels follow from the rule that a
// junction drains through its single steepest downhill edge:
//
//   * stolen slot. A main stem W -> A -> B -> C runs east, with tributaries
//     D -> A from the north and E -> B from the south (so B is a real
//     junction, not a degree-2 point that would be merged away). D sits in a
//     pit (elevation 0), so by DEM elevation the tributary runs uphill from D
//     to A: from A, the edge to D (84.5 m over 70 m) is steeper than the edge
//     to B (10 m over 100 m), so A drains to D. B drains to C and cannot drain
//     uphill to A. Nobody claims A-B: it is labelled "unclaimed", and since A-B
//     is the only link between the two halves, closes_cycle = 0 (the network
//     is cut). W-A, B-C and E-B are unlabelled;
//   * braid. P -> Q twice, once straight and once via a bend, with a line in
//     to P and a line out of Q. One of the two P-Q edges is labelled, and it
//     closes a real cycle: closes_cycle = 1.
//
// Usage: lsdtt-test-vector-loop-labels <work_dir>
// -----------------------------------------------------------------------------

#include <cmath>
#include <filesystem>
#include <iostream>
#include <string>
#include <vector>

#include <gdal_priv.h>
#include <ogrsf_frmts.h>

#include "../../src/cpu/lsdtt_raster_cpu.hpp"
#include "../../src/cpu/lsdtt_vector_to_junction_network.hpp"

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

constexpr double k_x0 = 500000.0;
constexpr double k_y0 = 6000000.0;
constexpr double k_cell = 10.0;
constexpr int k_cols = 50;
constexpr int k_rows = 30;

using lsdtt::lsdtt_xy;

// The edge whose end points are (near) p and q, or nullptr.
const lsdtt::lsdtt_vector_network_edge* edge_between(const lsdtt::lsdtt_vector_to_junction_network& net,
                                                     lsdtt_xy p, lsdtt_xy q, bool want_bent = false)
{
    for (const auto& e : net.edges())
    {
        if (e.xy.size() < 2) continue;
        const lsdtt_xy a = e.xy.front(), b = e.xy.back();
        const bool fwd = std::hypot(a.x - p.x, a.y - p.y) < 5 && std::hypot(b.x - q.x, b.y - q.y) < 5;
        const bool rev = std::hypot(a.x - q.x, a.y - q.y) < 5 && std::hypot(b.x - p.x, b.y - p.y) < 5;
        if (!(fwd || rev)) continue;
        // Largest distance of a vertex from the straight line p-q.
        const double len = std::hypot(q.x - p.x, q.y - p.y);
        double maxdev = 0.0;
        for (const auto& v : e.xy)
            maxdev = std::max(maxdev, std::fabs((q.x - p.x) * (p.y - v.y) - (p.x - v.x) * (q.y - p.y)) / len);
        if (want_bent == (maxdev > 10.0)) return &e;
    }
    return nullptr;
}

} // namespace

int main(int argc, char* argv[])
{
    if (argc < 2)
    {
        std::cerr << "usage: lsdtt-test-vector-loop-labels <work_dir>\n";
        return 2;
    }
    GDALAllRegister();
    const std::filesystem::path work = argv[1];
    std::filesystem::create_directories(work);

    const lsdtt_xy W{k_x0 + 50, k_y0 + 200}, A{k_x0 + 150, k_y0 + 200}, B{k_x0 + 250, k_y0 + 200},
                   C{k_x0 + 450, k_y0 + 200}, D{k_x0 + 150, k_y0 + 270},
                   E{k_x0 + 200, k_y0 + 140};
    const lsdtt_xy IN{k_x0 + 30, k_y0 + 60}, P{k_x0 + 150, k_y0 + 60}, Q{k_x0 + 350, k_y0 + 60},
                   OUT{k_x0 + 470, k_y0 + 60}, BEND{k_x0 + 250, k_y0 + 110};

    lsdtt::lsdtt_raster_meta meta;
    meta.n_rows = k_rows;
    meta.n_cols = k_cols;
    meta.x_min = k_x0;
    meta.y_min = k_y0;
    meta.cell_size = k_cell;
    meta.no_data_value = -9999.0f;
    meta.set_crs_from_epsg(32630);
    std::vector<float> data(static_cast<size_t>(k_rows * k_cols));
    std::vector<uint8_t> mask(data.size(), 1);
    for (int r = 0; r < k_rows; ++r)
        for (int c = 0; c < k_cols; ++c)
        {
            const double x = k_x0 + (c + 0.5) * k_cell;
            const double y = k_y0 + (k_rows - r - 0.5) * k_cell;
            double z = 100.0 - 0.1 * (x - k_x0);
            if (std::hypot(x - D.x, y - D.y) < 25.0) z = 0.0;   // the pit under D
            data[static_cast<size_t>(r * k_cols + c)] = static_cast<float>(z);
        }
    const lsdtt::lsdtt_raster_cpu dem(std::move(data), std::move(mask), meta);

    const std::string vec = (work / "loop_lines.fgb").string();
    {
        GDALDriver* drv = GetGDALDriverManager()->GetDriverByName("FlatGeobuf");
        std::filesystem::remove(vec);
        GDALDataset* ds = drv->Create(vec.c_str(), 0, 0, 0, GDT_Unknown, nullptr);
        OGRSpatialReference srs;
        srs.importFromEPSG(32630);
        OGRLayer* l = ds->CreateLayer("lines", &srs, wkbLineString, nullptr);
        auto put = [&](std::vector<lsdtt_xy> pts)
        {
            OGRFeature* f = OGRFeature::CreateFeature(l->GetLayerDefn());
            OGRLineString ls;
            for (const auto& p : pts) ls.addPoint(p.x, p.y);
            f->SetGeometry(&ls);
            const OGRErr err = l->CreateFeature(f);
            OGRFeature::DestroyFeature(f);
            if (err != OGRERR_NONE) throw std::runtime_error("could not write a test line");
        };
        put({W, A}); put({A, B}); put({B, C}); put({D, A}); put({E, B});
        put({IN, P}); put({P, Q}); put({P, BEND, Q}); put({Q, OUT});
        GDALClose(ds);
    }

    lsdtt::lsdtt_vector_to_junction_network_options opts;
    opts.break_cycles = true;
    lsdtt::lsdtt_vector_to_junction_network net;
    net.create(vec, dem, opts);

    const auto* ab = edge_between(net, A, B);
    check("stolen slot: edge A-B found", ab != nullptr);
    if (ab)
    {
        check("stolen slot: A-B labelled", ab->is_loop_edge);
        check("stolen slot: reason unclaimed", ab->loop_reason == "unclaimed", ab->loop_reason);
        check("stolen slot: A-B cuts the network (closes_cycle 0)", ab->closes_cycle == 0,
              std::to_string(ab->closes_cycle));
    }
    const auto* ad = edge_between(net, D, A);
    check("stolen slot: A drains to D (the steeper, uphill-along-the-data edge)",
          ad && !ad->is_loop_edge && std::hypot(ad->xy.front().x - A.x, ad->xy.front().y - A.y) < 5);
    for (auto pq : {std::make_pair(W, A), std::make_pair(B, C), std::make_pair(E, B)})
    {
        const auto* e = edge_between(net, pq.first, pq.second);
        check("stolen slot: other main-stem edges unlabelled",
              e && !e->is_loop_edge && e->loop_reason.empty() && e->closes_cycle == -1);
    }

    const auto* straight = edge_between(net, P, Q, false);
    const auto* bent = edge_between(net, P, Q, true);
    check("braid: both P-Q edges found", straight && bent);
    if (straight && bent)
    {
        const int n_lab = (straight->is_loop_edge ? 1 : 0) + (bent->is_loop_edge ? 1 : 0);
        check("braid: exactly one P-Q edge labelled", n_lab == 1, std::to_string(n_lab));
        const auto* lab = straight->is_loop_edge ? straight : bent;
        check("braid: the labelled edge closes a real cycle", lab->closes_cycle == 1,
              std::to_string(lab->closes_cycle));
        check("braid: it has a reason", !lab->loop_reason.empty());
    }

    std::cout << "\n" << (g_n_checks - g_n_failed) << " / " << g_n_checks << " checks passed\n";
    return g_n_failed == 0 ? 0 : 1;
}
