// test_common_parameters.cpp
//
// Direct tests for the selective shared CLI parameter declarations.

#include "lsdtt_common_parameters.hpp"

#include <iostream>
#include <map>
#include <string>

namespace
{

int n_failed = 0;

void check(bool condition, const std::string& label)
{
    if (condition)
    {
        std::cout << "PASS: " << label << "\n";
    }
    else
    {
        std::cerr << "FAIL: " << label << "\n";
        ++n_failed;
    }
}

std::string help_field(const lsdtt::lsdtt_help_registry& help,
                       const std::string& key,
                       std::size_t field)
{
    return lsdtt::lsdtt_help_registry::field(help.entries().at(key), field);
}

void test_independent_selection()
{
    std::map<std::string, bool> bools;
    lsdtt::lsdtt_help_registry help;

    lsdtt::declare_raster_is_filled_parameter(bools, help);
    check(bools.size() == 1 && bools.count("raster_is_filled") == 1,
          "mandatory-fill declaration adds only raster_is_filled");
    check(!bools.at("raster_is_filled"),
          "raster_is_filled default is false");
    check(help.entries().count("fill_before_routing") == 0,
          "mandatory-fill declaration does not advertise optional filling");
    check(help_field(help, "raster_is_filled", 3).find("REQUIRED") !=
              std::string::npos,
          "mandatory-fill help states that filling is required");

    bools.clear();
    help = lsdtt::lsdtt_help_registry{};
    lsdtt::declare_fill_before_routing_parameter(bools, help);
    check(bools.size() == 1 && bools.count("fill_before_routing") == 1,
          "optional-fill declaration adds only fill_before_routing");
    check(bools.at("fill_before_routing"),
          "fill_before_routing default is true");
    check(help.entries().count("raster_is_filled") == 0,
          "optional-fill declaration does not advertise a filled-input assertion");
    check(help_field(help, "fill_before_routing", 3).find("OPTIONAL") !=
              std::string::npos,
          "optional-fill help states that filling is optional");
}

void test_defaults_and_composition()
{
    std::map<std::string, bool> bools;
    std::map<std::string, int> ints;
    std::map<std::string, float> floats;
    std::map<std::string, std::string> strings;
    lsdtt::lsdtt_help_registry help;

    lsdtt::declare_routing_core_parameters(ints, strings, help);
    lsdtt::declare_raster_is_filled_parameter(bools, help);
    lsdtt::declare_fill_min_slope_parameter(floats, help);
    lsdtt::declare_channel_initiation_parameters(ints, help);

    check(strings.size() == 1 && strings.at("boundary_conditions") == "oooo",
          "routing core declares canonical boundary conditions");
    check(ints.at("n_threads") == 0,
          "routing core declares canonical thread default");
    check(ints.at("threshold_contributing_pixels") == 1000,
          "channel initiation declares canonical threshold");
    check(floats.size() == 1 && floats.at("fill_min_slope") == 0.0001f,
          "fill slope declaration has canonical default");
    check(help.entries().size() == 5,
          "composed declarations contain exactly the requested keys");
    check(help_field(help, "threshold_contributing_pixels", 1) == "1000",
          "channel threshold help default matches typed default");
}

void test_prefixes()
{
    std::map<std::string, bool> bools;
    std::map<std::string, int> ints;
    std::map<std::string, float> floats;
    std::map<std::string, std::string> strings;
    lsdtt::lsdtt_help_registry help;

    lsdtt::declare_routing_core_parameters(ints, strings, help, "coarse_");
    lsdtt::declare_fill_before_routing_parameter(bools, help, "coarse_");
    lsdtt::declare_fill_min_slope_parameter(floats, help, "coarse_");
    lsdtt::declare_channel_initiation_parameters(ints, help, "coarse_");

    check(strings.count("coarse_boundary_conditions") == 1,
          "trailing-underscore prefix does not produce a double underscore");
    check(ints.count("coarse_n_threads") == 1 &&
              ints.count("coarse_threshold_contributing_pixels") == 1,
          "prefix applies to routing and channel keys");
    check(bools.count("coarse_fill_before_routing") == 1 &&
              floats.count("coarse_fill_min_slope") == 1,
          "prefix applies to optional-fill keys");
    check(help_field(help, "coarse_fill_before_routing", 4) ==
              "Pit filling (coarse)",
          "prefix appears in the help group label");
}

void test_overwrite_policy()
{
    std::map<std::string, int> ints{{"threshold_contributing_pixels", 42}};
    lsdtt::lsdtt_help_registry help;
    help["threshold_contributing_pixels"] = {
        "int", "42", "local", "local", "Local", "param", ""};

    lsdtt::declare_channel_initiation_parameters(ints, help);
    check(ints.at("threshold_contributing_pixels") == 1000,
          "shared declaration overwrites a pre-existing typed default");
    check(help_field(help, "threshold_contributing_pixels", 1) == "1000",
          "shared declaration overwrites pre-existing help");
    check(help.key_order().size() == 1,
          "overwriting help does not duplicate insertion order");

    ints["threshold_contributing_pixels"] = 250;
    help["threshold_contributing_pixels"][1] = "250";
    check(ints.at("threshold_contributing_pixels") == 250 &&
              help_field(help, "threshold_contributing_pixels", 1) == "250",
          "caller can deliberately override typed and documented defaults after declaration");
}

} // namespace

int main()
{
    test_independent_selection();
    test_defaults_and_composition();
    test_prefixes();
    test_overwrite_policy();

    if (n_failed != 0)
    {
        std::cerr << n_failed << " common-parameter check(s) failed\n";
        return 1;
    }
    std::cout << "All common-parameter checks passed\n";
    return 0;
}
