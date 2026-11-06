// test_point_snapper.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
//
// -----------------------------------------------------------------------------
// lsdtt_point_snapper and build_line_endpoint_graph (lsdtt_vector_ops_cpu):
// line ends are joined by DISTANCE, not by rounding to a grid (river-lines plan,
// decision D9). Every expected id below follows from the rule "the nearest
// point already added within the tolerance, else a new id":
//
//   * two points 0.02 m apart either side of a 1 m grid line (x = 0.49 and
//     0.51, where the old llround(x / 1) key gave 0 and 1) are joined;
//   * two points 1.2 m apart in one old grid cell (llround gave the same key)
//     are NOT joined at a 1 m tolerance;
//   * exactly at the tolerance is joined (inclusive);
//   * the nearest of two candidates is chosen;
//   * tolerance 0 joins coincident points only;
//   * the endpoint graph of two lines meeting 0.02 m apart across a grid line
//     has three nodes, one of degree 2.
//
// Usage: lsdtt-test-point-snapper
// -----------------------------------------------------------------------------

#include <iostream>
#include <string>

#include "../../src/cpu/lsdtt_vector_ops_cpu.hpp"

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

lsdtt::lsdtt_vector_record line(std::vector<lsdtt::lsdtt_xy> pts)
{
    lsdtt::lsdtt_vector_record r;
    lsdtt::lsdtt_line_geom g;
    g.points = std::move(pts);
    r.geometry = g;
    return r;
}

} // namespace

int main()
{
    {
        lsdtt::lsdtt_point_snapper s(1.0);
        const int a = s.id_for({0.49, 0.0});
        const int b = s.id_for({0.51, 0.0});
        check("straddling a grid line, 0.02 m apart: joined", a == b,
              std::to_string(a) + " vs " + std::to_string(b));
        // The old key: llround(0.49) = 0, llround(0.51) = 1.
        check("(the old key split them)",
              !(lsdtt::make_snap_key(0.49, 0.0, 1.0) == lsdtt::make_snap_key(0.51, 0.0, 1.0)));
    }
    {
        lsdtt::lsdtt_point_snapper s(1.0);
        const int a = s.id_for({0.0, 0.0});
        const int b = s.id_for({0.6, 0.6});   // 0.849 m: joined
        const int c = s.id_for({-0.45, 0.4}); // 0.602 m from (0,0): joined
        check("within the tolerance: joined", a == b && a == c);
    }
    {
        lsdtt::lsdtt_point_snapper s(1.0);
        const int a = s.id_for({-0.45, -0.45});
        const int b = s.id_for({0.40, 0.40});  // 1.20 m: not joined
        check("1.2 m apart in one old grid cell: not joined", a != b);
        check("(the old key joined them)",
              lsdtt::make_snap_key(-0.45, -0.45, 1.0) == lsdtt::make_snap_key(0.40, 0.40, 1.0));
    }
    {
        lsdtt::lsdtt_point_snapper s(1.0);
        const int a = s.id_for({10.0, 10.0});
        check("exactly at the tolerance: joined (inclusive)", s.id_for({11.0, 10.0}) == a);
        check("just beyond: new id", s.id_for({11.0 + 1e-6, 10.0}) != a);
    }
    {
        lsdtt::lsdtt_point_snapper s(1.0);
        const int a = s.id_for({0.0, 0.0});
        const int b = s.id_for({1.5, 0.0});
        check("two separate points", a != b);
        check("nearest of two candidates", s.id_for({1.0, 0.0}) == b && s.id_for({0.7, 0.0}) == a);
        check("find never adds", s.find({50.0, 50.0}) == -1 && s.size() == 2);
        check("first point added is the representative", s.point(a).x == 0.0 && s.point(b).x == 1.5);
    }
    {
        lsdtt::lsdtt_point_snapper s(0.0);
        const int a = s.id_for({3.0, 4.0});
        check("tolerance 0: coincident joined", s.id_for({3.0, 4.0}) == a);
        check("tolerance 0: 1e-9 apart not joined", s.id_for({3.0 + 1e-9, 4.0}) != a);
    }
    {
        lsdtt::lsdtt_vector_layer layer;
        layer.records.push_back(line({{-10.0, 0.0}, {0.49, 0.0}}));
        layer.records.push_back(line({{0.51, 0.0}, {10.0, 0.0}}));
        const lsdtt::lsdtt_endpoint_graph g = lsdtt::build_line_endpoint_graph(layer, 1.0);
        check("endpoint graph across a grid line: 3 nodes", g.n_nodes == 3, std::to_string(g.n_nodes));
        int deg2 = 0;
        for (int d : g.node_degree) if (d == 2) ++deg2;
        check("endpoint graph: one node of degree 2", deg2 == 1);
    }

    std::cout << "\n" << (g_n_checks - g_n_failed) << " / " << g_n_checks << " checks passed\n";
    return g_n_failed == 0 ? 0 : 1;
}
