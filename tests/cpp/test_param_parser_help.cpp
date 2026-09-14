// test_param_parser_help.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
//
// Unit tests for the CLI help infrastructure
// (docs/developer_notes/lsdtt3_architecture.md, CLI design):
// lsdtt_help_registry insertion order, parse_args_deferred flag handling,
// and print_terminal_help rendering (grouping, switch/child nesting,
// 4-element entry fallback).

#include "lsdtt_param_parser.hpp"

#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace {

int g_failures = 0;

void check(bool condition, const std::string& label)
{
    if (condition)
    {
        std::cout << "PASS: " << label << "\n";
    }
    else
    {
        std::cerr << "FAIL: " << label << "\n";
        ++g_failures;
    }
}

// True if needle_a appears before needle_b in haystack (both must exist).
bool appears_before(const std::string& haystack,
                    const std::string& needle_a,
                    const std::string& needle_b)
{
    const std::size_t a = haystack.find(needle_a);
    const std::size_t b = haystack.find(needle_b);
    return a != std::string::npos && b != std::string::npos && a < b;
}

lsdtt::lsdtt_help_registry make_test_registry()
{
    lsdtt::lsdtt_help_registry help;

    // Deliberately non-alphabetical declaration order.
    help["write_fname"] = {"string","",
        "Output filename stem.", "Guidance for write_fname.",
        "", "io", ""};
    help["n_threads"] = {"int","0",
        "Number of threads.", "0 = all hardware threads.",
        "", "io", ""};

    help["write_zebra_raster"] = {"bool","false",
        "Write the zebra raster.", "Long zebra guidance.",
        "Zebra analysis", "switch", ""};
    help["zebra_stripe_width"] = {"float","2.5",
        "Stripe width in metres.", "Only used when zebra is on.",
        "Zebra analysis", "param", "write_zebra_raster"};

    help["write_aardvark_raster"] = {"bool","true",
        "Write the aardvark raster.", "Long aardvark guidance.",
        "Aardvark analysis", "switch", ""};
    help["aardvark_free_param"] = {"int","7",
        "A group-level parameter with no parent switch.", "More guidance.",
        "Aardvark analysis", "param", ""};

    // Legacy 4-element entry: no group/role/parent.
    help["legacy_param"] = {"float","1.0",
        "A legacy entry.", "Legacy guidance."};

    return help;
}

std::string render_overview(const lsdtt::lsdtt_help_registry& help)
{
    lsdtt::lsdtt_cli_args args;
    args.help_requested = true;
    std::ostringstream os;
    lsdtt::lsdtt_param_parser::print_terminal_help(
        "lsdtt-test-tool", "0.1", "https://example.org/citation",
        help, args, os);
    return os.str();
}

std::string render_topic(const lsdtt::lsdtt_help_registry& help,
                         const std::string& topic)
{
    lsdtt::lsdtt_cli_args args;
    args.help_requested = true;
    args.help_topic = topic;
    std::ostringstream os;
    lsdtt::lsdtt_param_parser::print_terminal_help(
        "lsdtt-test-tool", "0.1", "https://example.org/citation",
        help, args, os);
    return os.str();
}

void test_registry_insertion_order()
{
    lsdtt::lsdtt_help_registry help = make_test_registry();
    const auto& order = help.key_order();
    check(order.size() == 7, "registry records all keys");
    check(order.front() == "write_fname",
          "registry preserves first insertion");
    check(order[2] == "write_zebra_raster",
          "registry preserves declaration order, not alphabetical order");

    // Re-accessing an existing key must not duplicate it.
    help["write_fname"];
    check(help.key_order().size() == 7,
          "re-access does not duplicate key in order record");
}

void test_parse_args_deferred()
{
    {
        char prog[] = "tool"; char h[] = "-h";
        char* argv[] = {prog, h};
        auto args = lsdtt::lsdtt_param_parser::parse_args_deferred(
            2, argv, "tool", "0.1", "cite");
        check(args.help_requested && !args.help_html &&
              args.help_topic.empty(),
              "-h sets help_requested only");
    }
    {
        char prog[] = "tool"; char h[] = "-h"; char html[] = "html";
        char* argv[] = {prog, h, html};
        auto args = lsdtt::lsdtt_param_parser::parse_args_deferred(
            3, argv, "tool", "0.1", "cite");
        check(args.help_requested && args.help_html,
              "-h html sets help_html");
    }
    {
        char prog[] = "tool"; char h[] = "--help"; char t[] = "n_threads";
        char* argv[] = {prog, h, t};
        auto args = lsdtt::lsdtt_param_parser::parse_args_deferred(
            3, argv, "tool", "0.1", "cite");
        check(args.help_requested && args.help_topic == "n_threads",
              "--help <name> sets help_topic");
    }
    {
        char prog[] = "tool"; char f[] = "params.driver";
        char* argv[] = {prog, f};
        auto args = lsdtt::lsdtt_param_parser::parse_args_deferred(
            2, argv, "tool", "0.1", "cite");
        check(!args.help_requested && args.path == "./" &&
              args.fname == "params.driver",
              "two-arg form resolves fname with default path");
    }
    {
        char prog[] = "tool"; char p[] = "/data/"; char f[] = "params.driver";
        char* argv[] = {prog, p, f};
        auto args = lsdtt::lsdtt_param_parser::parse_args_deferred(
            3, argv, "tool", "0.1", "cite");
        check(!args.help_requested && args.path == "/data/" &&
              args.fname == "params.driver",
              "three-arg form resolves path and fname");
    }
}

void test_overview_rendering()
{
    const std::string out = render_overview(make_test_registry());

    check(out.find("lsdtt-test-tool 0.1") != std::string::npos,
          "overview shows tool name and version");
    check(out.find("File I/O and run control:") != std::string::npos,
          "overview has io section");
    check(appears_before(out, "File I/O and run control:", "Zebra analysis:"),
          "io section comes before analysis groups");
    check(appears_before(out, "Zebra analysis:", "Aardvark analysis:"),
          "groups appear in declaration order, not alphabetical");
    check(out.find("    zebra_stripe_width") != std::string::npos,
          "child parameter is indented beneath its switch");
    check(appears_before(out, "write_zebra_raster", "zebra_stripe_width"),
          "switch appears before its child parameter");
    check(out.find("aardvark_free_param") != std::string::npos,
          "group-level parameter without parent is listed");
    check(out.find("General:") != std::string::npos &&
          out.find("legacy_param") != std::string::npos,
          "legacy 4-element entry falls back to General group");
    check(out.find("Long zebra guidance.") == std::string::npos,
          "overview shows short descriptions only, not guidance text");

    // Parser-native file I/O keys are synthesised even though no CLI
    // declares them.
    for (const char* key : {"read_path", "read_fname",
                            "write_path", "write_fname"})
    {
        check(out.find(key) != std::string::npos,
              std::string("overview lists builtin io key ") + key);
    }
    check(appears_before(out, "read_path", "write_fname"),
          "builtin io keys are listed before CLI io keys");
    check(appears_before(out, "write_fname", "n_threads"),
          "builtin io keys come before declared io keys");
}

void test_builtin_io_topic_and_override()
{
    const auto help = make_test_registry();

    // -h read_fname works without the CLI declaring it.
    const std::string topic = render_topic(help, "read_fname");
    check(topic.find("Input filename") != std::string::npos,
          "topic view resolves builtin io key read_fname");

    // A CLI-declared entry for a builtin key wins over the synthesised one.
    lsdtt::lsdtt_help_registry custom = make_test_registry();
    custom["read_fname"] = {"string","",
        "Custom read_fname description.", "Custom guidance.",
        "", "io", ""};
    const std::string out = render_overview(custom);
    check(out.find("Custom read_fname description.") != std::string::npos &&
          out.find("Input filename, with or without") == std::string::npos,
          "CLI-declared builtin io entry overrides the synthesised one");
}

void test_topic_rendering()
{
    const auto help = make_test_registry();

    const std::string zebra = render_topic(help, "write_zebra_raster");
    check(zebra.find("Long zebra guidance.") != std::string::npos,
          "topic view shows the long guidance text");
    check(zebra.find("zebra_stripe_width") != std::string::npos,
          "topic view of a switch lists its dependent parameters");

    const std::string child = render_topic(help, "zebra_stripe_width");
    check(child.find("used by : write_zebra_raster") != std::string::npos,
          "topic view of a param names its parent switch");

    const std::string unknown = render_topic(help, "zebra");
    check(unknown.find("Unknown parameter") != std::string::npos,
          "unknown topic reports as unknown");
    check(unknown.find("write_zebra_raster") != std::string::npos,
          "unknown topic suggests close matches");
}

void test_json_help()
{
    {
        char prog[] = "tool"; char j[] = "--help-json";
        char* argv[] = {prog, j};
        auto args = lsdtt::lsdtt_param_parser::parse_args_deferred(
            2, argv, "tool", "0.1", "cite");
        check(args.help_requested && args.help_json,
              "--help-json sets help_json");
    }
    {
        char prog[] = "tool"; char h[] = "-h"; char j[] = "json";
        char* argv[] = {prog, h, j};
        auto args = lsdtt::lsdtt_param_parser::parse_args_deferred(
            3, argv, "tool", "0.1", "cite");
        check(args.help_requested && args.help_json,
              "-h json sets help_json");
    }

    lsdtt::lsdtt_cli_args args;
    args.help_requested = true;
    args.help_json = true;
    std::ostringstream os;
    lsdtt::lsdtt_param_parser::print_terminal_help(
        "lsdtt-test-tool", "0.1", "https://example.org/citation",
        make_test_registry(), args, os, {"fftw", "eigen"});
    const std::string out = os.str();

    check(out.find("\"tool\": \"lsdtt-test-tool\"") != std::string::npos,
          "json contains tool name");
    check(out.find("\"features\": [\"fftw\", \"eigen\"]") != std::string::npos,
          "json contains the supplied feature list");
    check(out.find("\"name\": \"write_zebra_raster\"") != std::string::npos &&
          out.find("\"role\": \"switch\"") != std::string::npos,
          "json contains parameters with role metadata");
    check(out.find("\"parent\": \"write_zebra_raster\"") != std::string::npos,
          "json contains parent metadata");
    check(out.find("\"name\": \"read_fname\"") != std::string::npos,
          "json includes the builtin io keys");
    check(out.rfind("{", 0) == 0,
          "json output starts with an opening brace");
}


// -----------------------------------------------------------------------------
// EVERY DECLARED PARAMETER RENDERS EXACTLY ONCE.
//
// A `parent` is a presentation hint. It must never be able to delete a
// parameter from the help output, and it must never be able to show one twice.
// Both used to happen: an entry was nested only when its parent was a switch in
// the SAME group, and every other parented entry fell through both passes and
// appeared nowhere. That silently hid 68 declared parameters across 8 CLIs -
// including every ladder_gap_* control in lsdtt-fluvial-surfaces - from -h,
// -h html and the per-run README files, while the CLI reference pages told
// users those outputs listed every parameter.
//
// The registry below reproduces every shape that broke it. Each case is a real
// declaration pattern found in the tree, named after where it was found.
// -----------------------------------------------------------------------------
lsdtt::lsdtt_help_registry make_awkward_parent_registry()
{
    lsdtt::lsdtt_help_registry help;

    help["write_thing_raster"] = {"bool","false",
        "A switch.", "", "Group A", "switch", ""};

    // 1. Parent is a switch in a DIFFERENT group (lsdtt-internal-drainage's
    //    sink_candidate_fname, lsdtt-hillslope-channel-coupling's trace_thinning).
    help["cross_group_child"] = {"int","1",
        "Parent is a switch, but in another group.", "",
        "Group B", "param", "write_thing_raster"};

    // 2. Parent is a "param", not a switch (lsdtt-lem's steady_state_* chain,
    //    lsdtt-fluvial-surfaces' ladder_gap_* under an algorithm choice).
    help["plain_parent"] = {"string","auto",
        "An algorithm choice, declared as a param.", "",
        "Group A", "param", ""};
    help["child_of_a_param"] = {"float","0.5",
        "Parent has role param.", "",
        "Group A", "param", "plain_parent"};

    // 3. A chain: child of a child of a param. The nesting pass descends
    //    through roles, so anything that re-derives "was this nested?" from the
    //    immediate parent's role disagrees with it and prints the middle twice.
    help["grandchild_of_a_param"] = {"int","3",
        "Two levels below a non-switch.", "",
        "Group A", "param", "child_of_a_param"};

    // 4. Parent is not a declared key at all: several CLIs put the ALLOWED
    //    VALUES or a CONDITION in the parent slot (lsdtt-raster-preprocessing's
    //    vertical_units, lsdtt-lem's cellular_* under "hillslope_model=cellular").
    help["parent_is_an_enum_list"] = {"string","metres",
        "Parent slot holds allowed values, not a key.", "",
        "Group A", "param", "metres|feet|us_survey_feet"};

    // 5. A switch nested under another switch (lsdtt-fluvial-surfaces'
    //    clean_valley_mask), itself carrying a child.
    help["nested_switch"] = {"bool","true",
        "A switch that depends on another switch.", "",
        "Group A", "switch", "write_thing_raster"};
    help["child_of_nested_switch"] = {"float","1.0",
        "Two levels down, under switches.", "",
        "Group A", "param", "nested_switch"};

    // 6. An io entry claimed as a parent. io is listed in its own block up
    //    front, so a tree must not print it a second time inside a group
    //    (lsdtt-fluvial-surfaces' fitted_surface_mask_fname).
    help["some_input_fname"] = {"string","",
        "An input path.", "", "Group A", "io", ""};
    help["child_of_an_io_key"] = {"int","2",
        "Parent has role io.", "",
        "Group A", "param", "some_input_fname"};

    return help;
}

// Count how many overview lines name `key`. The overview writes the key
// followed by " (" and its type, at some indent, so anchoring on that pair
// avoids matching the key where it is merely mentioned in prose.
int count_overview_entries(const std::string& rendered, const std::string& key)
{
    const std::string needle = key + " (";
    int n = 0;
    for (std::size_t pos = rendered.find(needle);
         pos != std::string::npos;
         pos = rendered.find(needle, pos + 1))
    {
        // Must be at the start of a line apart from indentation.
        std::size_t line_start = rendered.rfind('\n', pos);
        line_start = (line_start == std::string::npos) ? 0 : line_start + 1;
        bool only_spaces = true;
        for (std::size_t i = line_start; i < pos; ++i)
        {
            if (rendered[i] != ' ') { only_spaces = false; break; }
        }
        if (only_spaces) ++n;
    }
    return n;
}

void test_every_parameter_renders_exactly_once()
{
    const lsdtt::lsdtt_help_registry help = make_awkward_parent_registry();
    const std::string out = render_overview(help);

    for (const auto& key : help.key_order())
    {
        const int n = count_overview_entries(out, key);
        check(n == 1,
              "'" + key + "' renders exactly once (rendered " +
              std::to_string(n) + " times)");
    }

    // The nesting that DOES work must keep working: a switch's dependants stay
    // indented under it rather than being flattened out to fix the above.
    check(appears_before(out, "write_thing_raster", "nested_switch"),
          "a nested switch still renders beneath its parent switch");
    check(out.find("    nested_switch (") != std::string::npos,
          "a nested switch is indented, not promoted to the top level");
    check(out.find("      child_of_nested_switch (") != std::string::npos,
          "a grandchild is indented one level deeper again");

    // A dependency on a non-switch parameter must nest too. This is the shape
    // lsdtt-lem's cellular_nonlocal_* parameters have: they depend on
    // cellular_nonlocal_fraction, which is a float, not a switch. Before the
    // parent slot was corrected they declared "cellular_nonlocal_fraction>0" -
    // a CONDITION, not a key - and rendered nowhere at all.
    check(out.find("    child_of_a_param (") != std::string::npos,
          "a dependant of a non-switch parameter nests under it");
    check(out.find("      grandchild_of_a_param (") != std::string::npos,
          "the chain keeps nesting below a non-switch parameter");
    check(appears_before(out, "plain_parent", "child_of_a_param"),
          "a non-switch parent still precedes its dependants");
}

// -----------------------------------------------------------------------------
// THE PARENT SLOT NAMES A PARAMETER, NOT A VALUE SET OR A CONDITION.
//
// Several CLIs used to put allowed values ("metres|feet|us_survey_feet",
// "true|false") or a condition ("hillslope_model=cellular",
// "cellular_nonlocal_fraction>0") in the parent slot. Both read naturally and
// both were wrong: the renderer looks the parent up as a key, so those entries
// matched nothing and vanished from the help. They are fixed at the
// declarations - allowed values moved to the description, conditions to the
// guidance with the parent pointing at the actual key - and this test states
// the contract so the shape does not come back.
//
// A parent that names no declared key is not an error the parser can reject
// (a CLI may legitimately reference a key it declares conditionally), so the
// guarantee is the weaker but sufficient one: it never costs the parameter its
// place in the documentation.
// -----------------------------------------------------------------------------
void test_parent_slot_naming_a_non_key_still_renders()
{
    lsdtt::lsdtt_help_registry help;
    help["real_switch"] = {"bool","false",
        "A real switch.", "", "Group A", "switch", ""};
    help["values_in_parent_slot"] = {"string","metres",
        "Allowed values wrongly placed in the parent slot.", "",
        "Group A", "param", "metres|feet|us_survey_feet"};
    help["condition_in_parent_slot"] = {"float","0.85",
        "A condition wrongly placed in the parent slot.", "",
        "Group A", "param", "some_model=cellular"};

    const std::string out = render_overview(help);
    for (const auto& key : help.key_order())
    {
        check(count_overview_entries(out, key) == 1,
              "'" + key + "' survives a malformed parent slot");
    }
}

} // anonymous namespace

int main()
{
    test_registry_insertion_order();
    test_parse_args_deferred();
    test_overview_rendering();
    test_topic_rendering();
    test_builtin_io_topic_and_override();
    test_json_help();
    test_every_parameter_renders_exactly_once();
    test_parent_slot_naming_a_non_key_still_renders();

    if (g_failures > 0)
    {
        std::cerr << g_failures << " test(s) failed.\n";
        return 1;
    }
    std::cout << "All help infrastructure tests passed.\n";
    return 0;
}
