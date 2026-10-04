// test_hype_parameters.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// Unit tests for src/cpu/lsdtt_hype_parameters (S-HYPE column).
//
// What is under test
// ------------------
//   load_hype_setup_directory()
//     every case under tests/expected/hype with an lsdtt_input/ directory
//     loads. Facts from the committed tables and cases/*.json are pinned by
//     hand below, so a reader that swaps columns or misreads a blank cell is
//     caught without a JSON parser.
//
//   "absent means zero"
//     a parameter missing from general_parameters.csv (diversion has no tcalt)
//     and parameter COLUMNS missing from the land-use and soil tables (no
//     committed case does that, so a copy is written to the output directory)
//     read as 0.
//
//   hype_subbasin_downstream_order()
//     [1, 2] on default_all and routing; the upstream-first property on a
//     branched network given in reverse order.
//
//   validation
//     a deliberately broken copy of default_all per rule, written into the
//     output directory. Each must be rejected, and the message must name the
//     file and the offending name or value. validate_hype_setup() is also
//     exercised on a set-up built in memory.
//
// Usage: lsdtt-test-hype-parameters <oracle_dir> <output_dir>
// Prints one PASS/FAIL line per check and a final "RESULT: PASS|FAIL".
// -----------------------------------------------------------------------------

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include "lsdtt_hype_parameters.hpp"

namespace fs = std::filesystem;

namespace
{

int g_n_checks = 0;
int g_n_failed = 0;

void check(const std::string& name, bool ok, const std::string& detail = "")
{
    ++g_n_checks;
    if (!ok) ++g_n_failed;
    std::cout << (ok ? "PASS" : "FAIL") << "  " << name;
    if (!detail.empty()) std::cout << "   [" << detail << "]";
    std::cout << std::endl;
}

const lsdtt::hype_subbasin* find_subbasin(const lsdtt::hype_setup& s, int subid)
{
    for (const auto& b : s.subbasins)
        if (b.subid == subid) return &b;
    return nullptr;
}

const lsdtt::hype_class_definition* find_class(const lsdtt::hype_setup& s, int id)
{
    for (const auto& c : s.classes)
        if (c.class_id == id) return &c;
    return nullptr;
}

const lsdtt::hype_landuse_parameters* find_landuse(const lsdtt::hype_setup& s, int code)
{
    for (const auto& l : s.landuse)
        if (l.code == code) return &l;
    return nullptr;
}

const lsdtt::hype_soil_parameters* find_soil(const lsdtt::hype_setup& s, int code)
{
    for (const auto& so : s.soils)
        if (so.code == code) return &so;
    return nullptr;
}

std::vector<int> order_subids(const lsdtt::hype_setup& s)
{
    std::vector<int> ids;
    for (std::size_t i : lsdtt::hype_subbasin_downstream_order(s.subbasins))
        ids.push_back(s.subbasins[i].subid);
    return ids;
}

std::string ids_str(const std::vector<int>& v)
{
    std::string out = "[";
    for (std::size_t i = 0; i < v.size(); ++i) out += (i ? ", " : "") + std::to_string(v[i]);
    return out + "]";
}

void write_file(const fs::path& p, const std::string& content)
{
    std::ofstream f(p, std::ios::trunc);
    if (!f) throw std::runtime_error("cannot write " + p.string());
    f << content;
}

// Copy default_all's tables (not forcing.csv) into <output_dir>/<name>/ and
// overwrite one of them. Returns the new directory.
fs::path make_variant(const fs::path& oracle, const fs::path& out, const std::string& name,
                      const std::string& file, const std::string& content)
{
    const fs::path src = oracle / "default_all" / "lsdtt_input";
    const fs::path dst = out / name;
    fs::create_directories(dst);
    for (const char* f : {"general_parameters.csv", "landuse_parameters.csv", "soil_parameters.csv",
                          "classes.csv", "subbasins.csv", "subbasin_classes.csv"})
        fs::copy_file(src / f, dst / f, fs::copy_options::overwrite_existing);
    write_file(dst / file, content);
    return dst;
}

// ---------------------------------------------------------------------------
// The committed cases load, with hand-pinned facts
// ---------------------------------------------------------------------------
void test_oracle_cases(const fs::path& oracle)
{
    std::vector<std::string> cases;
    for (const auto& e : fs::directory_iterator(oracle))
        if (e.is_directory() && fs::is_directory(e.path() / "lsdtt_input"))
            cases.push_back(e.path().filename().string());
    std::sort(cases.begin(), cases.end());

    const std::vector<std::string> expected_cases = {"default_all", "diversion", "et", "percolation",
                                                     "routing", "snow", "soil_runoff"};
    for (const std::string& c : expected_cases)
        check("case directory present: " + c,
              std::find(cases.begin(), cases.end(), c) != cases.end());

    for (const std::string& c : cases)
    {
        try
        {
            const lsdtt::hype_setup s =
                lsdtt::load_hype_setup_directory((oracle / c / "lsdtt_input").string());
            check("loads: " + c, !s.subbasins.empty() && !s.classes.empty());
        }
        catch (const std::exception& e)
        {
            check("loads: " + c, false, e.what());
        }
    }
}

void test_pinned_facts(const fs::path& oracle)
{
    auto load = [&](const char* c)
    { return lsdtt::load_hype_setup_directory((oracle / c / "lsdtt_input").string()); };

    // default_all
    {
        const lsdtt::hype_setup s = load("default_all");
        check("default_all: 2 sub-basins", s.subbasins.size() == 2);
        const auto* b1 = find_subbasin(s, 1);
        const auto* b2 = find_subbasin(s, 2);
        check("default_all: subid 1 maindown is 2", b1 && b1->maindown == 2);
        check("default_all: subid 2 is the outlet", b2 && b2->maindown == 0);
        check("default_all: subid 2 area 35e6 m2, elev 150, slope 4 %",
              b2 && b2->area_m2 == 35000000.0 && b2->elev_mean_m == 150.0 &&
                  b2->slope_mean_percent == 4.0);
        check("default_all: subid 1 loc_rivlen 5000, rivlen 8000",
              b1 && b1->loc_rivlen_m == 5000.0 && b1->rivlen_m == 8000.0);
        check("default_all: elev_std_m absent reads 0", b1 && b1->elev_std_m == 0.0);
        check("default_all: 3 classes", s.classes.size() == 3);
        const auto* c3 = find_class(s, 3);
        check("default_all: class 3 is 1 layer, bottom 0.3 m, stream 0.25 m",
              c3 && c3->n_layers == 1 && c3->layer_bottom_m[0] == 0.3 &&
                  c3->layer_bottom_m[1] == 0.0 && c3->layer_bottom_m[2] == 0.0 &&
                  c3->stream_depth_m == 0.25 && c3->tile_depth_m == 0.0);
        check("default_all: class 3 is land use 2, soil 2", c3 && c3->landuse == 2 && c3->soil == 2);
        const auto* c2 = find_class(s, 2);
        check("default_all: class 2 is 2 layers, 0.15 / 0.6 m",
              c2 && c2->n_layers == 2 && c2->layer_bottom_m[0] == 0.15 &&
                  c2->layer_bottom_m[1] == 0.6 && c2->layer_bottom_m[2] == 0.0);
        const auto* c1 = find_class(s, 1);
        check("default_all: class 1 is 3 layers, 0.1 / 0.5 / 1.5 m, stream 1 m",
              c1 && c1->n_layers == 3 && c1->layer_bottom_m[2] == 1.5 && c1->stream_depth_m == 1.0);
        const auto* so2 = find_soil(s, 2);
        check("default_all: soil 2 wcep 0.08, wcfc 0.25, mperc1 8",
              so2 && so2->wcep == 0.08 && so2->wcfc == 0.25 && so2->mperc1 == 8.0);
        const auto* lu2 = find_landuse(s, 2);
        check("default_all: land use 2 pcluse 0.15, ttmp 0.5",
              lu2 && lu2->pcluse == 0.15 && lu2->ttmp == 0.5);
        check("default_all: general rivvel 0.7, damp 0.3, cevpph 20, rrcs3 0.002",
              s.general.rivvel == 0.7 && s.general.damp == 0.3 && s.general.cevpph == 20.0 &&
                  s.general.rrcs3 == 0.002);
        check("default_all: 6 sub-basin/class rows", s.subbasin_classes.size() == 6);
        check("default_all: (1,1) fraction 0.5, dh 80; (2,3) dh -40",
              s.subbasin_classes.size() == 6 && s.subbasin_classes[0].subid == 1 &&
                  s.subbasin_classes[0].class_id == 1 && s.subbasin_classes[0].fraction == 0.5 &&
                  s.subbasin_classes[0].dh_m == 80.0 && s.subbasin_classes[5].dh_m == -40.0);
        const std::vector<int> ord = order_subids(s);
        check("default_all: downstream order [1, 2]", ord == std::vector<int>{1, 2}, ids_str(ord));
        check("default_all: source_dir recorded", !s.source_dir.empty());
    }

    // soil_runoff
    {
        const lsdtt::hype_setup s = load("soil_runoff");
        const auto* b1 = find_subbasin(s, 1);
        const auto* b2 = find_subbasin(s, 2);
        check("soil_runoff: slopes 2 % and 25 %",
              b1 && b2 && b1->slope_mean_percent == 2.0 && b2->slope_mean_percent == 25.0);
        check("soil_runoff: sub-basins not connected", b1 && b2 && b1->maindown == 0 && b2->maindown == 0);
        const auto* c1 = find_class(s, 1);
        const auto* c2 = find_class(s, 2);
        const auto* c3 = find_class(s, 3);
        check("soil_runoff: stream depths 0.05 / 0.3 / 1 m",
              c1 && c2 && c3 && c1->stream_depth_m == 0.05 && c2->stream_depth_m == 0.3 &&
                  c3->stream_depth_m == 1.0);
        check("soil_runoff: rrcs3 0.004, macrate 0.2",
              s.general.rrcs3 == 0.004 && !s.soils.empty() && s.soils[0].macrate == 0.2);
    }

    // routing
    {
        const lsdtt::hype_setup s = load("routing");
        const auto* b2 = find_subbasin(s, 2);
        check("routing: subid 2 rivlen 100000 m, loc_rivlen 20000 m",
              b2 && b2->rivlen_m == 100000.0 && b2->loc_rivlen_m == 20000.0);
        check("routing: rivvel 0.5, damp 0.4", s.general.rivvel == 0.5 && s.general.damp == 0.4);
        const std::vector<int> ord = order_subids(s);
        check("routing: downstream order [1, 2]", ord == std::vector<int>{1, 2}, ids_str(ord));
    }

    // diversion: tcalt, fsceff, whcsnow, rrcs3, cevpam, cevpph are not in its
    // general_parameters.csv, so they are zero.
    {
        const lsdtt::hype_setup s = load("diversion");
        check("diversion: absent general parameters read 0 (tcalt, whcsnow, rrcs3, cevpam, cevpph)",
              s.general.tcalt == 0.0 && s.general.whcsnow == 0.0 && s.general.rrcs3 == 0.0 &&
                  s.general.cevpam == 0.0 && s.general.cevpph == 0.0 && s.general.fsceff == 0.0);
        check("diversion: lp 0.9, srrate 0.3, mactrinf 8, ttmp -10",
              s.general.lp == 0.9 && s.soils[0].srrate == 0.3 && s.soils[0].mactrinf == 8.0 &&
                  s.landuse[0].ttmp == -10.0);
    }

    // et
    {
        const lsdtt::hype_setup s = load("et");
        check("et: cevpam 0.4, cevpph 30, epotdist 4, cevp 0.25",
              s.general.cevpam == 0.4 && s.general.cevpph == 30.0 && s.general.epotdist == 4.0 &&
                  s.landuse[0].cevp == 0.25);
    }

    // percolation
    {
        const lsdtt::hype_setup s = load("percolation");
        check("percolation: mperc1 6, mperc2 3, layer 3 bottom 0.9 m",
              s.soils[0].mperc1 == 6.0 && s.soils[0].mperc2 == 3.0 &&
                  s.classes[0].layer_bottom_m[2] == 0.9);
    }

    // snow
    {
        const lsdtt::hype_setup s = load("snow");
        check("snow: dh 150 m, tcalt 0.6, whcsnow 0.08, cmlt 3",
              s.subbasin_classes[0].dh_m == 150.0 && s.general.tcalt == 0.6 &&
                  s.general.whcsnow == 0.08 && s.landuse[0].cmlt == 3.0);
    }
}

// ---------------------------------------------------------------------------
// Absent parameter columns read as zero (no committed case drops one)
// ---------------------------------------------------------------------------
void test_absent_columns(const fs::path& oracle, const fs::path& out)
{
    try
    {
        // pcluse and srrcs dropped from land use.
        fs::path d = make_variant(oracle, out, "absent_columns", "landuse_parameters.csv",
                                  "code,cevp,cmlt,ttmp\n1,0.2,3.5,0\n2,0.25,2.5,0.5\n");
        // macrate, mactrinf and rrcs2 dropped from soil.
        write_file(d / "soil_parameters.csv",
                   "code,mactrsm,mperc1,mperc2,rrcs1,srrate,wcep,wcfc,wcwp\n"
                   "1,0.6,20,8,0.25,0.1,0.14,0.22,0.08\n"
                   "2,0.7,8,3,0.15,0.2,0.08,0.25,0.15\n");
        const lsdtt::hype_setup s = lsdtt::load_hype_setup_directory(d.string());
        const auto* lu2 = find_landuse(s, 2);
        check("absent land-use columns pcluse, srrcs read 0; present ones kept",
              lu2 && lu2->pcluse == 0.0 && lu2->srrcs == 0.0 && lu2->cmlt == 2.5);
        const auto* so1 = find_soil(s, 1);
        check("absent soil columns macrate, mactrinf, rrcs2 read 0; present ones kept",
              so1 && so1->macrate == 0.0 && so1->mactrinf == 0.0 && so1->rrcs2 == 0.0 &&
                  so1->mperc1 == 20.0);
    }
    catch (const std::exception& e)
    {
        check("absent parameter columns load", false, e.what());
    }
}

// ---------------------------------------------------------------------------
// Downstream order on a branched network given downstream-first
// ---------------------------------------------------------------------------
void test_downstream_order()
{
    // 1 -> 2 -> 3 (outlet); 4 -> 3; 5 isolated.
    std::vector<lsdtt::hype_subbasin> b(5);
    const int ids[5] = {3, 2, 5, 1, 4};
    const int downs[5] = {0, 3, 0, 2, 3};
    for (int i = 0; i < 5; ++i)
    {
        b[static_cast<std::size_t>(i)].subid = ids[i];
        b[static_cast<std::size_t>(i)].maindown = downs[i];
        b[static_cast<std::size_t>(i)].area_m2 = 1.0;
    }
    const std::vector<std::size_t> ord = lsdtt::hype_subbasin_downstream_order(b);
    std::vector<int> pos(6, -1);
    for (std::size_t k = 0; k < ord.size(); ++k) pos[static_cast<std::size_t>(b[ord[k]].subid)] = static_cast<int>(k);
    bool ok = ord.size() == 5;
    for (const auto& s : b)
        if (s.maindown != 0 && !(pos[static_cast<std::size_t>(s.subid)] < pos[static_cast<std::size_t>(s.maindown)])) ok = false;
    std::vector<int> ids_out;
    for (std::size_t i : ord) ids_out.push_back(b[i].subid);
    check("branched network: every sub-basin before its maindown", ok, ids_str(ids_out));
    // Kahn's with input-order ties: headwaters 5, 1, 4 first, then 2, then 3.
    check("branched network: deterministic order [5, 1, 4, 2, 3]",
          ids_out == std::vector<int>{5, 1, 4, 2, 3}, ids_str(ids_out));

    b[0].maindown = 1;   // 3 -> 1 closes 1 -> 2 -> 3 -> 1
    bool threw = false;
    std::string msg;
    try { (void)lsdtt::hype_subbasin_downstream_order(b); }
    catch (const std::invalid_argument& e) { threw = true; msg = e.what(); }
    check("in-memory cycle rejected by hype_subbasin_downstream_order",
          threw && msg.find("cycle") != std::string::npos, msg);
}

// ---------------------------------------------------------------------------
// Broken tables
// ---------------------------------------------------------------------------
struct broken_case
{
    const char* name;
    const char* file;
    std::string content;
    std::vector<std::string> must_contain;
};

void test_rejections(const fs::path& oracle, const fs::path& out)
{
    const char* classes_hdr =
        "class_id,landuse,soil,n_layers,layer1_bottom_m,layer2_bottom_m,layer3_bottom_m,"
        "stream_depth_m,tile_depth_m\n";
    const std::string cl_ok2 = "2,2,1,2,0.15,0.6,,0.5,0\n3,2,2,1,0.3,,,0.25,0\n";
    const char* sb_hdr = "subid,maindown,area_m2,elev_mean_m,slope_mean_percent,loc_rivlen_m,rivlen_m\n";
    const char* sc_hdr = "subid,class_id,fraction,dh_m\n";
    const std::string sc_sub2 = "2,1,0.2,30\n2,2,0.4,0\n2,3,0.4,-40\n";
    const char* lu_hdr = "code,cevp,cmlt,pcluse,srrcs,ttmp\n";
    const char* so_hdr = "code,macrate,mactrinf,mactrsm,mperc1,mperc2,rrcs1,rrcs2,srrate,wcep,wcfc,wcwp\n";
    const std::string so_row2 = "2,0.05,6,0.7,8,3,0.15,0.01,0.2,0.08,0.25,0.15\n";

    std::vector<broken_case> cases;

    cases.push_back({"unknown_general_name", "general_parameters.csv",
                     "name,value\nttpd,0.5\nttpx,1.5\n",
                     {"general_parameters.csv", "row 2", "unknown parameter", "ttpx"}});
    cases.push_back({"duplicate_general_name", "general_parameters.csv",
                     "name,value\nlp,0.75\nlp,0.8\n",
                     {"general_parameters.csv", "row 2", "\"lp\" given twice"}});
    cases.push_back({"general_value_not_numeric", "general_parameters.csv",
                     "name,value\nlp,abc\n", {"general_parameters.csv", "row 1", "value"}});
    cases.push_back({"unknown_landuse_column", "landuse_parameters.csv",
                     "code,cevp,cmlt,pcluse,srrcs,ttmpp\n1,0.2,3.5,0.05,0.15,0\n2,0.25,2.5,0.15,0.1,0.5\n",
                     {"landuse_parameters.csv", "unknown column", "ttmpp"}});
    cases.push_back({"unknown_soil_column", "soil_parameters.csv",
                     (std::string("code,macrate,mactrinf,mactrsm,mperc1,mperc3,rrcs1,rrcs2,srrate,wcep,wcfc,wcwp\n") +
                          "1,0.15,12,0.6,20,8,0.25,0.02,0.1,0.14,0.22,0.08\n" + so_row2),
                     {"soil_parameters.csv", "unknown column", "mperc3"}});
    cases.push_back({"duplicate_landuse_code", "landuse_parameters.csv",
                     (std::string(lu_hdr) + "1,0.2,3.5,0.05,0.15,0\n2,0.25,2.5,0.15,0.1,0.5\n1,0.2,3.5,0.05,0.15,0\n"),
                     {"landuse_parameters.csv", "row 3", "code", "duplicate code 1"}});
    cases.push_back({"duplicate_soil_code", "soil_parameters.csv",
                     (std::string(so_hdr) + "2,0.15,12,0.6,20,8,0.25,0.02,0.1,0.14,0.22,0.08\n" + so_row2),
                     {"soil_parameters.csv", "row 2", "duplicate code 2"}});
    cases.push_back({"landuse_code_not_integer", "landuse_parameters.csv",
                     (std::string(lu_hdr) + "1.5,0.2,3.5,0.05,0.15,0\n2,0.25,2.5,0.15,0.1,0.5\n"),
                     {"landuse_parameters.csv", "row 1", "\"code\"", "1.5", "not an integer"}});
    cases.push_back({"duplicate_class_id", "classes.csv",
                     (std::string(classes_hdr) + "2,1,1,3,0.1,0.5,1.5,1,0\n" + cl_ok2),
                     {"classes.csv", "row 2", "duplicate class_id 2"}});
    cases.push_back({"class_missing_landuse", "classes.csv",
                     (std::string(classes_hdr) + "1,9,1,3,0.1,0.5,1.5,1,0\n" + cl_ok2),
                     {"classes.csv", "row 1", "\"landuse\"", "land-use code 9"}});
    cases.push_back({"class_missing_soil", "classes.csv",
                     (std::string(classes_hdr) + "1,1,7,3,0.1,0.5,1.5,1,0\n" + cl_ok2),
                     {"classes.csv", "row 1", "\"soil\"", "soil code 7"}});
    cases.push_back({"n_layers_out_of_range", "classes.csv",
                     (std::string(classes_hdr) + "1,1,1,4,0.1,0.5,1.5,1,0\n" + cl_ok2),
                     {"classes.csv", "row 1", "n_layers", "4"}});
    cases.push_back({"n_layers_zero", "classes.csv",
                     (std::string(classes_hdr) + "1,1,1,0,0.1,0.5,1.5,1,0\n" + cl_ok2),
                     {"classes.csv", "row 1", "n_layers", "0"}});
    cases.push_back({"n_layers_fewer_than_depths", "classes.csv",
                     (std::string(classes_hdr) + "1,1,1,2,0.1,0.5,1.5,1,0\n" + cl_ok2),
                     {"classes.csv", "row 1", "layer3_bottom_m", "n_layers is 2"}});
    cases.push_back({"n_layers_more_than_depths", "classes.csv",
                     (std::string(classes_hdr) + "1,1,1,3,0.1,0.5,,1,0\n" + cl_ok2),
                     {"classes.csv", "row 1", "layer3_bottom_m", "blank"}});
    cases.push_back({"depths_not_increasing", "classes.csv",
                     (std::string(classes_hdr) + "1,1,1,3,0.1,0.5,0.5,1,0\n" + cl_ok2),
                     {"classes.csv", "row 1", "layer3_bottom_m", "deeper"}});
    cases.push_back({"depth_not_positive", "classes.csv",
                     (std::string(classes_hdr) + "1,1,1,3,0,0.5,1.5,1,0\n" + cl_ok2),
                     {"classes.csv", "row 1", "layer1_bottom_m", "positive"}});
    cases.push_back({"stream_depth_negative", "classes.csv",
                     (std::string(classes_hdr) + "1,1,1,3,0.1,0.5,1.5,-0.1,0\n" + cl_ok2),
                     {"classes.csv", "row 1", "stream_depth_m", "negative"}});
    cases.push_back({"duplicate_subid", "subbasins.csv",
                     (std::string(sb_hdr) + "2,0,20000000,350,12,5000,8000\n2,0,35000000,150,4,7000,25000\n"),
                     {"subbasins.csv", "row 2", "duplicate subid 2"}});
    cases.push_back({"maindown_missing", "subbasins.csv",
                     (std::string(sb_hdr) + "1,7,20000000,350,12,5000,8000\n2,0,35000000,150,4,7000,25000\n"),
                     {"subbasins.csv", "row 1", "maindown", "maindown 7"}});
    cases.push_back({"maindown_cycle", "subbasins.csv",
                     (std::string(sb_hdr) + "1,2,20000000,350,12,5000,8000\n2,1,35000000,150,4,7000,25000\n"),
                     {"subbasins.csv", "maindown", "cycle", "1 -> 2 -> 1"}});
    cases.push_back({"maindown_self", "subbasins.csv",
                     (std::string(sb_hdr) + "1,1,20000000,350,12,5000,8000\n2,0,35000000,150,4,7000,25000\n"),
                     {"subbasins.csv", "row 1", "cycle", "1 -> 1"}});
    cases.push_back({"area_zero", "subbasins.csv",
                     (std::string(sb_hdr) + "1,2,0,350,12,5000,8000\n2,0,35000000,150,4,7000,25000\n"),
                     {"subbasins.csv", "row 1", "area_m2", "subid 1"}});
    cases.push_back({"subid_zero", "subbasins.csv",
                     (std::string(sb_hdr) + "1,0,20000000,350,12,5000,8000\n0,0,35000000,150,4,7000,25000\n"),
                     {"subbasins.csv", "row 2", "subid", "positive"}});
    cases.push_back({"subbasins_missing_column", "subbasins.csv",
                     "subid,maindown,area_m2,elev_mean_m,slope_mean_percent,loc_rivlen_m\n"
                     "1,2,20000000,350,12,5000\n2,0,35000000,150,4,7000\n",
                     {"subbasins.csv", "missing required column", "rivlen_m"}});
    cases.push_back({"fraction_out_of_range", "subbasin_classes.csv",
                     (std::string(sc_hdr) + "1,1,1.5,80\n1,2,-0.3,-20\n1,3,-0.2,-60\n" + sc_sub2),
                     {"subbasin_classes.csv", "row 1", "fraction", "1.5", "[0, 1]"}});
    cases.push_back({"fractions_do_not_sum", "subbasin_classes.csv",
                     (std::string(sc_hdr) + "1,1,0.5,80\n1,2,0.3,-20\n1,3,0.1,-60\n" + sc_sub2),
                     {"subbasin_classes.csv", "fraction", "subid 1", "sum to 0.9"}});
    cases.push_back({"subbasin_without_classes", "subbasin_classes.csv",
                     (std::string(sc_hdr) + sc_sub2),
                     {"subbasin_classes.csv", "subid 1", "sum to 0"}});
    cases.push_back({"unknown_subid_in_classes", "subbasin_classes.csv",
                     (std::string(sc_hdr) + "1,1,0.5,80\n1,2,0.3,-20\n1,3,0.2,-60\n9,1,1,0\n" + sc_sub2),
                     {"subbasin_classes.csv", "row 4", "subid 9"}});
    cases.push_back({"unknown_class_id_in_classes", "subbasin_classes.csv",
                     (std::string(sc_hdr) + "1,1,0.5,80\n1,2,0.3,-20\n1,8,0.2,-60\n" + sc_sub2),
                     {"subbasin_classes.csv", "row 3", "class_id 8"}});
    cases.push_back({"duplicate_subbasin_class_pair", "subbasin_classes.csv",
                     (std::string(sc_hdr) + "1,1,0.5,80\n1,2,0.3,-20\n1,2,0.2,-60\n" + sc_sub2),
                     {"subbasin_classes.csv", "row 3", "duplicate row", "class_id 2"}});

    for (const broken_case& bc : cases)
    {
        std::string msg;
        bool threw = false;
        try
        {
            const fs::path d = make_variant(oracle, out, bc.name, bc.file, bc.content);
            (void)lsdtt::load_hype_setup_directory(d.string());
        }
        catch (const std::exception& e)
        {
            threw = true;
            msg = e.what();
        }
        bool named = threw;
        std::string missing;
        for (const std::string& m : bc.must_contain)
            if (msg.find(m) == std::string::npos)
            {
                named = false;
                missing += " '" + m + "'";
            }
        check(std::string("rejects ") + bc.name, threw && named,
              threw ? (named ? msg : "message lacks" + missing + ": " + msg) : "did not throw");
    }
}

// ---------------------------------------------------------------------------
// validate_hype_setup on a set-up built in memory
// ---------------------------------------------------------------------------
void test_in_memory_validation(const fs::path& oracle)
{
    const lsdtt::hype_setup base =
        lsdtt::load_hype_setup_directory((oracle / "default_all" / "lsdtt_input").string());

    lsdtt::hype_setup s = base;
    s.source_dir.clear();
    bool ok = true;
    try { lsdtt::validate_hype_setup(s); } catch (...) { ok = false; }
    check("in-memory copy of default_all validates", ok);

    s.classes[2].layer_bottom_m[1] = 0.6;   // class 3 has n_layers 1
    std::string msg;
    try { lsdtt::validate_hype_setup(s); } catch (const std::invalid_argument& e) { msg = e.what(); }
    check("in-memory: depth beyond n_layers rejected, table named",
          msg.find("classes.csv table") != std::string::npos &&
              msg.find("layer2_bottom_m") != std::string::npos,
          msg);

    s = base;
    s.source_dir.clear();
    s.subbasin_classes[0].fraction = 0.5 + 2e-6;   // just outside the 1e-6 tolerance
    msg.clear();
    try { lsdtt::validate_hype_setup(s); } catch (const std::invalid_argument& e) { msg = e.what(); }
    check("in-memory: fraction sum 1 + 2e-6 rejected", msg.find("subid 1") != std::string::npos, msg);

    s.subbasin_classes[0].fraction = 0.5 + 5e-7;   // inside it
    ok = true;
    try { lsdtt::validate_hype_setup(s); } catch (...) { ok = false; }
    check("in-memory: fraction sum 1 + 5e-7 accepted", ok);
}

} // namespace


// Parameter ranges (validate_hype_setup): values HYPE would clamp or that make
// the column fail are rejected, naming the table and column.
void parameter_range_checks(const fs::path& oracle_dir)
{
    std::cout << "\n--- parameter ranges ---\n";
    const auto base = lsdtt::load_hype_setup_directory((oracle_dir / "default_all" / "lsdtt_input").string());
    auto expect_reject = [&](const std::string& what, const std::function<void(lsdtt::hype_setup&)>& mutate,
                             const std::string& needle)
    {
        lsdtt::hype_setup s = base;
        mutate(s);
        std::string msg;
        bool threw = false;
        try { lsdtt::validate_hype_setup(s); }
        catch (const std::invalid_argument& e) { threw = true; msg = e.what(); }
        check("range: rejects " + what, threw && msg.find(needle) != std::string::npos, msg);
    };
    expect_reject("negative srrcs", [](auto& s) { s.landuse[0].srrcs = -10.0; }, "srrcs");
    expect_reject("pcluse > 1", [](auto& s) { s.landuse[0].pcluse = 1.5; }, "pcluse");
    expect_reject("negative wcwp", [](auto& s) { s.soils[0].wcwp = -0.1; }, "wcwp");
    expect_reject("negative mperc1", [](auto& s) { s.soils[0].mperc1 = -1.0; }, "mperc1");
    expect_reject("wcwp+wcfc+wcep > 1", [](auto& s) { s.soils[0].wcep = 0.9; }, "wcwp+wcfc+wcep");
    expect_reject("damp > 1", [](auto& s) { s.general.damp = 1.2; }, "damp");
    expect_reject("whcsnow = 1", [](auto& s) { s.general.whcsnow = 1.0; }, "whcsnow");
    expect_reject("negative rivvel", [](auto& s) { s.general.rivvel = -0.5; }, "rivvel");
    expect_reject("NaN sub-basin slope", [](auto& s) { s.subbasins[0].slope_mean_percent = std::numeric_limits<double>::quiet_NaN(); }, "slope_mean_percent");
    expect_reject("infinite class dh", [](auto& s) { s.subbasin_classes[0].dh_m = std::numeric_limits<double>::infinity(); }, "dh_m");
    bool ok = true;
    try { lsdtt::validate_hype_setup(base); } catch (...) { ok = false; }
    check("range: the unmodified default_all set-up still validates", ok);
}

int main(int argc, char** argv)
{
    std::cout << "lsdtt-test-hype-parameters\n"
              << "--------------------------\n";
    if (argc != 3)
    {
        std::cerr << "usage: lsdtt-test-hype-parameters <oracle_dir> <output_dir>\n";
        return 2;
    }
    const fs::path oracle = argv[1];
    const fs::path out = argv[2];
    try
    {
        fs::create_directories(out);
        test_oracle_cases(oracle);
        test_pinned_facts(oracle);
        test_absent_columns(oracle, out);
        test_downstream_order();
        test_rejections(oracle, out);
        test_in_memory_validation(oracle);
        parameter_range_checks(oracle);
    }
    catch (const std::exception& e)
    {
        std::cout << "FAIL  unexpected exception: " << e.what() << std::endl;
        ++g_n_failed; ++g_n_checks;
    }

    std::cout << "--------------------------\n"
              << g_n_checks - g_n_failed << " / " << g_n_checks << " checks passed\n"
              << "RESULT: " << (g_n_failed == 0 ? "PASS" : "FAIL") << std::endl;
    return g_n_failed == 0 ? 0 : 1;
}
