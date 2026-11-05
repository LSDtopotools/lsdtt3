// test_river_line_prepare.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
//
// -----------------------------------------------------------------------------
// Tests for lsdtt_river_line_prepare_cpu (lsdtt-fetch-river-lines), in process.
// Every expected count, degree, type and coordinate below is worked out by hand
// from the synthetic polylines built in each case:
//
//   * a Y network with a T-junction on a segment interior, a near-miss
//     T-junction at an interior vertex, a reversed tributary and a near-miss
//     endpoint pair;
//   * a braid (two edges on a loop), a closed ring (a self-loop);
//   * a multipart feature and a degenerate part;
//   * clipping: a crossing line, a line that leaves and re-enters, a line
//     entirely outside;
//   * merging by source node id beyond the snap tolerance;
//   * direction conflicts;
//   * dropping isolated short fragments only;
//   * zero tolerance: only coincident endpoints merge;
//   * the graph-only direction vote: one reversed reach in a main stem is the
//     one line against the network, and is flipped when repair is on;
//   * a sink whose outflow lies outside the data is the outlet, not a conflict.
//
// Usage: lsdtt-test-river-line-prepare
// -----------------------------------------------------------------------------

#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include "lsdtt_river_line_prepare_cpu.hpp"

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

void check_int(const std::string& name, long got, long want)
{
    check(name, got == want, "got " + std::to_string(got) + ", want " + std::to_string(want));
}

bool near(double a, double b, double eps = 1e-9) { return std::fabs(a - b) <= eps; }

lsdtt::river_line_feature line(std::vector<lsdtt::lsdtt_xy> pts, int direction = 1,
                               const std::string& id = "")
{
    lsdtt::river_line_feature f;
    f.parts.push_back(std::move(pts));
    f.direction = direction;
    f.direction_source = direction == 0 ? "unknown" : "attribute";
    f.source_feature_id = id;
    f.feature_class = "river";
    return f;
}

// The node at (x, y), or nullptr.
const lsdtt::river_node* node_at(const lsdtt::river_line_network& net, double x, double y,
                                 double eps = 1e-6)
{
    for (const auto& n : net.nodes)
    {
        if (near(n.xy.x, x, eps) && near(n.xy.y, y, eps)) return &n;
    }
    return nullptr;
}

std::string type_at(const lsdtt::river_line_network& net, double x, double y)
{
    const lsdtt::river_node* n = node_at(net, x, y);
    return n ? n->node_type : std::string("<no node>");
}

int count_type(const lsdtt::river_line_network& net, const std::string& t)
{
    int c = 0;
    for (const auto& n : net.nodes) if (n.node_type == t) ++c;
    return c;
}

void test_y_network()
{
    // Main stem M flows west to east, (0,0) -> (100,0) -> (200,0).
    // T1 ends 0.4 m above M's interior vertex (100,0): a near-miss T-junction.
    // T2 is digitised upstream from (150,0), the middle of M's second segment.
    // T3 ends at (0.5,0.5), 0.71 m from M's start: a near-miss endpoint pair.
    std::vector<lsdtt::river_line_feature> f;
    f.push_back(line({{0, 0}, {100, 0}, {200, 0}}, 1, "M"));
    f.push_back(line({{100, 100}, {100, 0.4}}, 1, "T1"));
    f.push_back(line({{150, 0}, {150, -80}}, -1, "T2"));
    f.push_back(line({{0, 50}, {0.5, 0.5}}, 1, "T3"));

    lsdtt::river_line_prepare_options o;
    o.node_snap_tolerance_metres = 1.0;
    const auto net = lsdtt::prepare_river_lines(f, o);

    check_int("Y: lines (M in three pieces + three tributaries)", static_cast<long>(net.lines.size()), 6);
    check_int("Y: nodes", static_cast<long>(net.nodes.size()), 7);
    check_int("Y: endpoints merged", net.n_endpoints_merged, 1);
    check_int("Y: T-junction splits", net.n_t_junction_splits, 2);
    check_int("Y: parts reversed", net.n_parts_reversed, 1);
    check_int("Y: components", static_cast<long>(net.components.size()), 1);
    check_int("Y: loops", net.components.empty() ? -1 : net.components[0].n_loops, 0);
    // Joined by lsdtt_point_snapper: the node sits at the first endpoint in
    // input order, M's start (0,0); T3's end (0.5,0.5) moves onto it.
    check("Y: joined endpoint at M's start (0,0) is a link",
          type_at(net, 0, 0) == "link", type_at(net, 0, 0));
    check("Y: near-miss T-junction node keeps T1's end (100,0.4), confluence",
          type_at(net, 100, 0.4) == "confluence", type_at(net, 100, 0.4));
    check("Y: T2's junction (150,0) is a confluence", type_at(net, 150, 0) == "confluence",
          type_at(net, 150, 0));
    check("Y: reversed T2 starts at a tip (150,-80)", type_at(net, 150, -80) == "tip",
          type_at(net, 150, -80));
    check("Y: outlet at (200,0)", type_at(net, 200, 0) == "outlet", type_at(net, 200, 0));
    check_int("Y: tips", count_type(net, "tip"), 3);

    double m_len = 0.0;
    int m_pieces = 0;
    bool m_splits_ok = true;
    bool all_downstream = true;
    for (const auto& l : net.lines)
    {
        if (l.source_feature_id == "M")
        {
            ++m_pieces;
            m_len += l.length_m;
            if (l.n_splits != 2) m_splits_ok = false;
        }
        if (l.on_loop) all_downstream = false;
        // every line's last vertex is its to-node
        const auto& to = net.nodes[static_cast<size_t>(l.to_node)].xy;
        if (!near(l.vertices.back().x, to.x) || !near(l.vertices.back().y, to.y)) all_downstream = false;
    }
    check_int("Y: M pieces", m_pieces, 3);
    check("Y: each M piece records 2 splits", m_splits_ok);
    // M: (0,0)->(100,0.4)->(150,0)->(200,0); the detour to y=0.4 adds a little.
    const double want = std::hypot(100.0, 0.4) + std::hypot(50.0, 0.4) + 50.0;
    check("Y: M length after snapping", near(m_len, want, 1e-9),
          std::to_string(m_len) + " vs " + std::to_string(want));
    check("Y: no loops flagged; every line ends at its to-node", all_downstream);
}

void test_braid_and_ring()
{
    std::vector<lsdtt::river_line_feature> f;
    f.push_back(line({{-50, 0}, {0, 0}}, 1, "up"));
    f.push_back(line({{0, 0}, {100, 0}}, 1, "arm_a"));
    f.push_back(line({{0, 0}, {50, 30}, {100, 0}}, 1, "arm_b"));
    f.push_back(line({{100, 0}, {150, 0}}, 1, "down"));
    // A separate closed ring (a lake shore digitised as a line).
    f.push_back(line({{500, 500}, {510, 500}, {510, 510}, {500, 500}}, 1, "ring"));

    lsdtt::river_line_prepare_options o;
    const auto net = lsdtt::prepare_river_lines(f, o);

    bool arms_on_loop = true, others_off = true, ring_on_loop = false;
    for (const auto& l : net.lines)
    {
        if (l.source_feature_id == "arm_a" || l.source_feature_id == "arm_b") arms_on_loop &= l.on_loop;
        else if (l.source_feature_id == "ring") ring_on_loop = l.on_loop;
        else others_off &= !l.on_loop;
    }
    check("braid: both arms on_loop", arms_on_loop);
    check("braid: stem lines not on_loop", others_off);
    check("ring: self-loop on_loop", ring_on_loop);
    check("braid: (0,0) bifurcation", type_at(net, 0, 0) == "bifurcation", type_at(net, 0, 0));
    check("braid: (100,0) confluence", type_at(net, 100, 0) == "confluence", type_at(net, 100, 0));
    check_int("braid+ring: components", static_cast<long>(net.components.size()), 2);
    int loops = 0;
    for (const auto& c : net.components) loops += c.n_loops;
    check_int("braid+ring: independent cycles", loops, 2);
    check_int("braid+ring: bifurcations", count_type(net, "bifurcation"), 1);
}

void test_multipart_and_degenerate()
{
    lsdtt::river_line_feature f;
    f.parts.push_back({{0, 0}, {10, 0}});
    f.parts.push_back({{100, 0}, {110, 0}});
    f.parts.push_back({{200, 0}, {200, 0}, {200, 0}});   // degenerate
    f.start_node_id = "X";                              // ignored for multipart
    f.end_node_id = "X";
    f.direction = 1;
    f.source_feature_id = "multi";
    const auto net = lsdtt::prepare_river_lines({f}, lsdtt::river_line_prepare_options());
    check_int("multipart: parts read", net.n_input_parts, 3);
    check_int("multipart: degenerate dropped", net.n_parts_dropped_degenerate, 1);
    check_int("multipart: lines", static_cast<long>(net.lines.size()), 2);
    check_int("multipart: ids ignored, two components", static_cast<long>(net.components.size()), 2);
    bool attr = net.lines.size() == 2 && net.lines[0].source_feature_id == "multi"
                && net.lines[1].source_feature_id == "multi";
    check("multipart: attributes copied to each part", attr);
}

void test_clipping()
{
    lsdtt::river_line_prepare_options o;
    o.clip_to_box = true;
    o.clip_x_min = 0; o.clip_y_min = 0; o.clip_x_max = 100; o.clip_y_max = 100;

    std::vector<lsdtt::river_line_feature> f;
    f.push_back(line({{-50, 50}, {150, 50}}, 1, "cross"));
    f.push_back(line({{30, 80}, {30, 150}, {60, 150}, {60, 80}}, 1, "out_and_back"));
    f.push_back(line({{300, 300}, {400, 300}}, 1, "outside"));
    const auto net = lsdtt::prepare_river_lines(f, o);

    check_int("clip: parts outside the box", net.n_parts_clipped_away, 1);
    int cross = 0, oab = 0;
    double cross_len = 0, oab_len = 0;
    for (const auto& l : net.lines)
    {
        if (l.source_feature_id == "cross") { ++cross; cross_len += l.length_m; }
        if (l.source_feature_id == "out_and_back") { ++oab; oab_len += l.length_m; }
    }
    check_int("clip: crossing line, one piece", cross, 1);
    check("clip: crossing piece length 100", near(cross_len, 100.0), std::to_string(cross_len));
    check_int("clip: out-and-back line, two pieces", oab, 2);
    check("clip: out-and-back pieces 20 + 20", near(oab_len, 40.0), std::to_string(oab_len));
    check("clip: (0,50) clip_edge", type_at(net, 0, 50) == "clip_edge", type_at(net, 0, 50));
    check("clip: (100,50) clip_edge", type_at(net, 100, 50) == "clip_edge", type_at(net, 100, 50));
    check("clip: (30,100) clip_edge", type_at(net, 30, 100) == "clip_edge", type_at(net, 30, 100));
    check("clip: (60,100) clip_edge", type_at(net, 60, 100) == "clip_edge", type_at(net, 60, 100));
    check("clip: uncut start (30,80) is a tip", type_at(net, 30, 80) == "tip", type_at(net, 30, 80));
    check("clip: uncut end (60,80) is an outlet", type_at(net, 60, 80) == "outlet", type_at(net, 60, 80));
    check_int("clip: clip_edge nodes", count_type(net, "clip_edge"), 4);

    std::vector<bool> cs, ce;
    const auto pieces = lsdtt::clip_polyline_to_box({{10, 10}, {90, 90}}, 0, 0, 100, 100, cs, ce);
    check("clip: inside line unchanged, no cut flags",
          pieces.size() == 1 && pieces[0].size() == 2 && !cs[0] && !ce[0]);
}

void test_node_ids()
{
    // Two tributaries and a stem whose ends are 5 m apart but share an id.
    lsdtt::river_line_feature a = line({{0, 100}, {0, 5}}, 1, "a");
    a.end_node_id = "J";
    lsdtt::river_line_feature b = line({{100, 0}, {5, 0}}, 1, "b");
    b.end_node_id = "J";
    lsdtt::river_line_feature c = line({{0, 0}, {-100, -100}}, 1, "c");
    c.start_node_id = "J";
    lsdtt::river_line_prepare_options o;
    o.node_snap_tolerance_metres = 1.0;
    const auto net = lsdtt::prepare_river_lines({a, b, c}, o);
    check_int("ids: nodes (J merged beyond tolerance)", static_cast<long>(net.nodes.size()), 4);
    // The node sits at the first endpoint with id J in input order: a's end (0,5).
    check("ids: J at a's end (0,5), the first endpoint with that id, is a confluence",
          type_at(net, 0, 5) == "confluence", type_at(net, 0, 5));
}

void test_conflicts()
{
    // h1 and h2 run head to head into (50,0); h2 starts at (100,0), where a
    // long line also leaves for (500,0). The vote puts the outlet at (500,0)
    // (only h2's 50 m disagrees, against 400 m for an outlet at (50,0)). With
    // repair off, (50,0) has two inflows and no outflow and is not the outlet,
    // and (100,0) has two outflows and no inflow: both are conflicts.
    std::vector<lsdtt::river_line_feature> f;
    f.push_back(line({{0, 0}, {50, 0}}, 1, "h1"));
    f.push_back(line({{100, 0}, {50, 0}}, 1, "h2"));
    f.push_back(line({{100, 0}, {500, 0}}, 1, "long"));
    lsdtt::river_line_prepare_options no_repair;
    no_repair.repair_direction = false;
    const auto net = lsdtt::prepare_river_lines(f, no_repair);
    check("conflict: head-to-head (50,0)", type_at(net, 50, 0) == "conflict", type_at(net, 50, 0));
    check("conflict: two outflows, no inflow (100,0)", type_at(net, 100, 0) == "conflict",
          type_at(net, 100, 0));
    check_int("conflict: count", count_type(net, "conflict"), 2);
    check("conflict: outlet voted at (500,0)", type_at(net, 500, 0) == "outlet", type_at(net, 500, 0));
}

void test_fragments()
{
    std::vector<lsdtt::river_line_feature> f;
    f.push_back(line({{0, 0}, {100, 0}}, 1, "stem"));
    f.push_back(line({{100, 0}, {103, 0}}, 1, "short_connected"));
    f.push_back(line({{500, 0}, {503, 0}}, 1, "short_isolated"));
    lsdtt::river_line_prepare_options o;
    o.min_line_length_metres = 10.0;
    const auto net = lsdtt::prepare_river_lines(f, o);
    bool has_connected = false, has_isolated = false;
    for (const auto& l : net.lines)
    {
        if (l.source_feature_id == "short_connected") has_connected = true;
        if (l.source_feature_id == "short_isolated") has_isolated = true;
    }
    check_int("fragments: dropped", net.n_fragments_dropped, 1);
    check("fragments: connected short line kept", has_connected);
    check("fragments: isolated short line dropped", !has_isolated);
    check_int("fragments: nodes renumbered", static_cast<long>(net.nodes.size()), 3);
    bool ids_ok = true;
    for (size_t i = 0; i < net.nodes.size(); ++i) ids_ok &= net.nodes[i].node_id == static_cast<int>(i);
    check("fragments: node ids are 0..n-1", ids_ok);
}

void test_zero_tolerance()
{
    std::vector<lsdtt::river_line_feature> f;
    f.push_back(line({{0, 0}, {10, 0}}, 1, "a"));
    f.push_back(line({{10, 0}, {20, 0}}, 1, "b"));        // coincident
    f.push_back(line({{20.5, 0}, {30, 0}}, 1, "c"));      // 0.5 m gap
    lsdtt::river_line_prepare_options o;
    o.node_snap_tolerance_metres = 0.0;
    const auto net = lsdtt::prepare_river_lines(f, o);
    check_int("zero tol: nodes", static_cast<long>(net.nodes.size()), 5);
    check_int("zero tol: components", static_cast<long>(net.components.size()), 2);
    o.node_snap_tolerance_metres = 1.0;
    const auto net1 = lsdtt::prepare_river_lines(f, o);
    check_int("1 m tol: the 0.5 m gap closes", static_cast<long>(net1.components.size()), 1);
}


void test_direction_vote()
{
    // Main stem (0,0) -> (400,0) in four 100 m reaches, with tributaries
    // joining at (100,0) and (200,0). The reach (100,0)-(200,0) is digitised
    // the wrong way round (attribute says downstream, but it runs 200 -> 100).
    // Votes, by length of line disagreeing:
    //   outlet (400,0): only the reversed reach disagrees: 100 m;
    //   outlet (100,0) (a data sink: three lines end there): the two reaches
    //     below (200,0) disagree: 200 m;
    //   any tip: more. So the outlet is (400,0) and the reversed reach is the
    //   one line against the network.
    std::vector<lsdtt::river_line_feature> f;
    f.push_back(line({{0, 0}, {100, 0}}, 1, "m1"));
    f.push_back(line({{200, 0}, {100, 0}}, 1, "m2_reversed"));
    f.push_back(line({{200, 0}, {300, 0}}, 1, "m3"));
    f.push_back(line({{300, 0}, {400, 0}}, 1, "m4"));
    f.push_back(line({{100, 100}, {100, 0}}, 1, "t1"));
    f.push_back(line({{200, 100}, {200, 0}}, 1, "t2"));

    lsdtt::river_line_prepare_options off;
    off.repair_direction = false;
    const auto raw = lsdtt::prepare_river_lines(f, off);
    check_int("vote, repair off: one line against the network", raw.n_lines_against_network, 1);
    check_int("vote, repair off: nothing flipped", raw.n_lines_flipped, 0);
    bool flagged_right = true;
    for (const auto& l : raw.lines)
        flagged_right &= (l.direction_consistent == (l.source_feature_id == "m2_reversed" ? 0 : 1));
    check("vote, repair off: only m2_reversed has direction_consistent 0", flagged_right);
    check("vote, repair off: the reversed reach leaves a conflict at (100,0)",
          type_at(raw, 100, 0) == "conflict", type_at(raw, 100, 0));
    check_int("vote: the chosen outlet is (400,0)",
              raw.components.size() == 1 ? raw.components[0].outlet_node : -1,
              node_at(raw, 400, 0) ? node_at(raw, 400, 0)->node_id : -2);

    lsdtt::river_line_prepare_options on;   // repair is the default
    const auto fixed = lsdtt::prepare_river_lines(f, on);
    check_int("vote, repair on: one line flipped", fixed.n_lines_flipped, 1);
    bool m2_ok = false;
    for (const auto& l : fixed.lines)
    {
        if (l.source_feature_id != "m2_reversed") continue;
        m2_ok = near(l.vertices.front().x, 100) && near(l.vertices.back().x, 200) &&
                l.direction_source == "network" && l.direction_consistent == 0;
    }
    check("vote, repair on: m2 now runs 100 -> 200, direction_source network", m2_ok);
    check("vote, repair on: (100,0) is a confluence", type_at(fixed, 100, 0) == "confluence",
          type_at(fixed, 100, 0));
    check("vote, repair on: (200,0) is a confluence", type_at(fixed, 200, 0) == "confluence",
          type_at(fixed, 200, 0));
    check_int("vote, repair on: no conflicts or bifurcations",
              count_type(fixed, "conflict") + count_type(fixed, "bifurcation"), 0);
    check("vote, repair on: outlet at (400,0)", type_at(fixed, 400, 0) == "outlet", type_at(fixed, 400, 0));
}

void test_truncated_outlet()
{
    // Two tributaries end at (100,0) and nothing leaves it: the outflow lies
    // outside the data. The sink is the outlet (no disagreement), nothing is
    // flipped, and it is typed outlet, not conflict.
    std::vector<lsdtt::river_line_feature> f;
    f.push_back(line({{0, 50}, {100, 0}}, 1, "a"));
    f.push_back(line({{0, -50}, {100, 0}}, 1, "b"));
    f.push_back(line({{-100, 50}, {0, 50}}, 1, "c"));
    const auto net = lsdtt::prepare_river_lines(f, lsdtt::river_line_prepare_options());
    check_int("truncated: nothing flipped", net.n_lines_flipped, 0);
    check("truncated: the sink (100,0) is the outlet", type_at(net, 100, 0) == "outlet",
          type_at(net, 100, 0));
}
} // namespace

int main()
{
    test_y_network();
    test_braid_and_ring();
    test_multipart_and_degenerate();
    test_clipping();
    test_node_ids();
    test_conflicts();
    test_fragments();
    test_zero_tolerance();
    test_direction_vote();
    test_truncated_outlet();
    std::cout << "\n" << (g_n_checks - g_n_failed) << " / " << g_n_checks << " checks passed\n";
    return g_n_failed == 0 ? 0 : 1;
}
