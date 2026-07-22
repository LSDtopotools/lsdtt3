// test_param_parser_help.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
//
// Unit tests for the CLI help infrastructure (docs/plan_docs/cli_help_plan.md phase 1):
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

} // anonymous namespace

int main()
{
    test_registry_insertion_order();
    test_parse_args_deferred();
    test_overview_rendering();
    test_topic_rendering();
    test_builtin_io_topic_and_override();
    test_json_help();

    if (g_failures > 0)
    {
        std::cerr << g_failures << " test(s) failed.\n";
        return 1;
    }
    std::cout << "All help infrastructure tests passed.\n";
    return 0;
}
