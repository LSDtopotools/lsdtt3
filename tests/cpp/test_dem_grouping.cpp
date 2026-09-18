// test_dem_grouping.cpp

#include <filesystem>
#include <iostream>
#include <string>
#include <vector>
#include <fstream>

#include "lsdtt_dem_grouping.hpp"
#include "lsdtt_outlet_bounds_discovery.hpp"

namespace {

int g_checks = 0;
int g_failed = 0;

void check(const std::string& name, bool ok)
{
    ++g_checks;
    if (!ok) ++g_failed;
    std::cout << (ok ? "PASS" : "FAIL") << "  " << name << "\n";
}

lsdtt::lsdtt_outlet_bounds_record make_record(
    int id, double lon, double lat,
    double west, double south, double east, double north)
{
    lsdtt::lsdtt_outlet_bounds_record r;
    r.point_id = id;
    r.longitude = lon;
    r.latitude = lat;
    r.status = "ok";
    r.west = west;
    r.south = south;
    r.east = east;
    r.north = north;
    return r;
}

void test_disconnect_split()
{
    std::vector<lsdtt::lsdtt_outlet_bounds_record> records;
    records.push_back(make_record(1, 120.0, -4.0, 119.9, -4.1, 120.1, -3.9));
    records.push_back(make_record(2, 120.1, -4.0, 120.0, -4.1, 120.2, -3.9));
    records.push_back(make_record(3, 122.0, -4.0, 121.9, -4.1, 122.1, -3.9));

    lsdtt::lsdtt_dem_grouping_options o;
    o.auto_split_disconnected_basins = true;
    o.split_cluster_gap_km = 25.0;
    o.max_final_dem_area_km2 = 1.0e9;
    o.max_basins_per_dem = 0;

    const auto groups = lsdtt::group_outlet_bounds(records, o);
    check("disconnect -> two groups", groups.size() == 2);
}

void test_max_basins_split_and_overlap()
{
    std::vector<lsdtt::lsdtt_outlet_bounds_record> records;
    records.push_back(make_record(1, 120.0, -4.0, 119.95, -4.05, 120.05, -3.95));
    records.push_back(make_record(2, 120.02, -4.0, 119.97, -4.05, 120.07, -3.95));

    lsdtt::lsdtt_dem_grouping_options o;
    o.auto_split_disconnected_basins = false;
    o.max_basins_per_dem = 1;
    o.max_final_dem_area_km2 = 1.0e9;
    o.split_overlap_km = 2.0;

    const auto groups = lsdtt::group_outlet_bounds(records, o);
    check("max_basins_per_dem=1 splits", groups.size() == 2);

    if (groups.size() == 2)
    {
        check("overlap expanded bbox west", groups[0].west < records[0].west || groups[1].west < records[0].west);
    }
}

// Where the test writes its scratch files. Passed by CTest as argv[1] (a
// directory under the build tree); falls back to <system temp>/lsdtt3_test_*.
// This used to be a hard-coded /tmp/opencode, which only existed because some
// other tool had created it on one workstation: on a clean machine every write
// failed (found 2026-09-18 by the minimal Apptainer image gate).
std::string g_work_dir;

void test_outlet_discovery_csv_single_group()
{
    const std::string csv_path = g_work_dir + "/test_outlet_points_e5.csv";
    {
        std::ofstream out(csv_path);
        out << "latitude,longitude\n";
        out << "-4.0000,120.0000\n";
        out << "-4.0020,120.0030\n";
    }

    const auto points = lsdtt::load_outlet_points_wgs84(csv_path);
    check("outlet discovery csv loads two points", points.size() == 2);

    lsdtt::lsdtt_outlet_bounds_discovery_options o;
    o.initial_radius_km = 10.0;
    o.max_iters = 2;
    const auto discovered = lsdtt::discover_outlet_bounds(points, o);
    check("outlet discovery result ok", discovered.ok);
    check("outlet discovery gives merged bbox", discovered.merged_west < discovered.merged_east && discovered.merged_south < discovered.merged_north);

    lsdtt::lsdtt_dem_grouping_options gopt;
    gopt.auto_split_disconnected_basins = true;
    gopt.split_cluster_gap_km = 25.0;
    gopt.max_final_dem_area_km2 = 1.0e9;
    const auto groups = lsdtt::group_outlet_bounds(discovered.records, gopt);
    check("nearby outlets produce one group", groups.size() == 1);

    const std::string status_csv = g_work_dir + "/test_outlet_status_e5.csv";
    std::string error;
    const bool wrote = lsdtt::write_outlet_bounds_status_csv(status_csv, discovered, error);
    check("status csv write ok", wrote && error.empty());
}

} // namespace

int main(int argc, char** argv)
{
    g_work_dir = (argc > 1)
        ? std::string(argv[1])
        : (std::filesystem::temp_directory_path() / "lsdtt3_test_dem_grouping").string();
    std::filesystem::create_directories(g_work_dir);

    test_disconnect_split();
    test_max_basins_split_and_overlap();
    test_outlet_discovery_csv_single_group();

    std::cout << "RESULT: " << (g_failed == 0 ? "PASS" : "FAIL")
              << " (" << (g_checks - g_failed) << "/" << g_checks << ")\n";
    return g_failed == 0 ? 0 : 1;
}
