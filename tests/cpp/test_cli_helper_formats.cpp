// test_cli_helper_formats.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html
//
// -----------------------------------------------------------------------------
// Phase 2 CLI-consistency helper pins.
//
// Covers:
//   B2  lsdtt::output_extension(raster_output_format, read_ext)
//                                           (src/core/lsdtt_data_loader.hpp)
//       plus an overload-resolution regression on the two pre-existing
//       output_extension siblings -- adding a third overload with a defaulted
//       second argument is the one thing that could have broken them.
//   B3  lsdtt::csv_number                   (src/core/lsdtt_stat_csv.hpp)
//
// The precision-9 and precision-10 expectations below are the byte-exact
// strings that cli/lsdtt-swath-profile.cpp:64 (fmt) and
// cli/lsdtt-basin-extraction.cpp:85 (precise_string) produce today. Stage 2
// migrates those call sites onto csv_number(v, 9) and csv_number(v, 10); these
// assertions are what guarantees their published CSVs do not change.
// -----------------------------------------------------------------------------

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <string>

#include "../../src/core/lsdtt_data_loader.hpp"
#include "../../src/core/lsdtt_stat_csv.hpp"

namespace
{

int g_failures = 0;

void check_string(
    const std::string& what,
    const std::string& got,
    const std::string& expected)
{
    if (got != expected)
    {
        std::printf("FAIL %s: expected \"%s\", got \"%s\"\n",
                    what.c_str(), expected.c_str(), got.c_str());
        ++g_failures;
    }
}

void check_true(const std::string& what, bool ok)
{
    if (!ok)
    {
        std::printf("FAIL %s: expected true\n", what.c_str());
        ++g_failures;
    }
}

// -----------------------------------------------------------------------------
// B2: output_extension(raster_output_format, read_ext)
// -----------------------------------------------------------------------------
void test_raster_output_extension()
{
    using lsdtt::output_extension;
    using lsdtt::raster_output_format;

    check_string("output_extension(png)",
                 output_extension(raster_output_format::png), ".png");
    check_string("output_extension(cog)",
                 output_extension(raster_output_format::cog), ".tif");
    check_string("output_extension(gtiff)",
                 output_extension(raster_output_format::gtiff), ".tif");

    // png ignores read_ext.
    check_string("output_extension(png, \".bil\") ignores read_ext",
                 output_extension(raster_output_format::png, ".bil"), ".png");

    // same_as_input with no input extension falls back to .tif.
    check_string("output_extension(same_as_input, \"\")",
                 output_extension(raster_output_format::same_as_input), ".tif");
    check_string("output_extension(same_as_input, \"\") explicit",
                 output_extension(raster_output_format::same_as_input, ""),
                 ".tif");

    // same_as_input echoes a non-empty input extension verbatim.
    check_string("output_extension(same_as_input, \".tif\")",
                 output_extension(raster_output_format::same_as_input, ".tif"),
                 ".tif");
    check_string("output_extension(same_as_input, \".vrt\")",
                 output_extension(raster_output_format::same_as_input, ".vrt"),
                 ".vrt");
    check_string("output_extension(same_as_input, \".img\")",
                 output_extension(raster_output_format::same_as_input, ".img"),
                 ".img");
    check_string("output_extension(same_as_input, \".bil\")",
                 output_extension(raster_output_format::same_as_input, ".bil"),
                 ".bil");
}

// -----------------------------------------------------------------------------
// Overload-resolution regression: the two pre-existing siblings must still
// resolve unambiguously now that a third overload with a defaulted second
// parameter sits alongside them.
// -----------------------------------------------------------------------------
void test_output_extension_overloads()
{
    using lsdtt::output_extension;

    check_string("output_extension(vector: flatgeobuf)",
                 output_extension(lsdtt::vector_output_format::flatgeobuf),
                 ".fgb");
    check_string("output_extension(vector: gpkg)",
                 output_extension(lsdtt::vector_output_format::gpkg), ".gpkg");

    check_string("output_extension(point: csv)",
                 output_extension(lsdtt::point_output_format::csv), ".csv");
    check_string("output_extension(point: flatgeobuf)",
                 output_extension(lsdtt::point_output_format::flatgeobuf),
                 ".fgb");
    check_string("output_extension(point: gpkg)",
                 output_extension(lsdtt::point_output_format::gpkg), ".gpkg");

    // The same call through a const lvalue of each enum type, which is how the
    // CLIs actually call it.
    const lsdtt::vector_output_format vfmt = lsdtt::vector_output_format::gpkg;
    const lsdtt::point_output_format  pfmt = lsdtt::point_output_format::csv;
    const lsdtt::raster_output_format rfmt = lsdtt::raster_output_format::png;
    check_string("output_extension(const vector lvalue)",
                 output_extension(vfmt), ".gpkg");
    check_string("output_extension(const point lvalue)",
                 output_extension(pfmt), ".csv");
    check_string("output_extension(const raster lvalue)",
                 output_extension(rfmt), ".png");
}

// -----------------------------------------------------------------------------
// B3: csv_number
// -----------------------------------------------------------------------------
void test_csv_number_finite()
{
    using lsdtt::csv_number;

    check_string("csv_number(0.0)",  csv_number(0.0),  "0");
    check_string("csv_number(-0.5)", csv_number(-0.5), "-0.5");
    check_string("csv_number(-123456.789012345)",
                 csv_number(-123456.789012345), "-123456.78901234501");

    // Default precision 17 must round-trip a value that needs all 17 digits.
    {
        const double v   = 1.0 / 3.0;
        const std::string s = csv_number(v);
        check_string("csv_number(1/3) at default precision",
                     s, "0.33333333333333331");
        const double back = std::strtod(s.c_str(), nullptr);
        check_true("csv_number(1/3) round-trips bit-exactly", back == v);
    }
    {
        const double v   = 0.1;
        const double back = std::strtod(csv_number(v).c_str(), nullptr);
        check_true("csv_number(0.1) round-trips bit-exactly", back == v);
    }

    // Precision 9: byte-identical to cli/lsdtt-swath-profile.cpp:64 fmt().
    check_string("csv_number(pi, 9) == swath-profile fmt()",
                 csv_number(3.14159265358979, 9), "3.14159265");
    check_string("csv_number(2/7, 9) == swath-profile fmt()",
                 csv_number(2.0 / 7.0, 9), "0.285714286");
    check_string("csv_number(1234567890.12345, 9) == swath-profile fmt()",
                 csv_number(1234567890.12345, 9), "1.23456789e+09");

    // Precision 10: byte-identical to cli/lsdtt-basin-extraction.cpp:85
    // precise_string().
    check_string("csv_number(pi, 10) == basin-extraction precise_string()",
                 csv_number(3.14159265358979, 10), "3.141592654");
    check_string("csv_number(2/7, 10) == basin-extraction precise_string()",
                 csv_number(2.0 / 7.0, 10), "0.2857142857");
    check_string("csv_number(1234567890.12345, 10) == "
                 "basin-extraction precise_string()",
                 csv_number(1234567890.12345, 10), "1234567890");
    check_string("csv_number(-123456.789012345, 10) == "
                 "basin-extraction precise_string()",
                 csv_number(-123456.789012345, 10), "-123456.789");
}

void test_csv_number_non_finite()
{
    using lsdtt::csv_number;

    const double nan_v  = std::numeric_limits<double>::quiet_NaN();
    const double inf_v  = std::numeric_limits<double>::infinity();

    check_string("csv_number(NaN)",  csv_number(nan_v),  "nan");
    check_string("csv_number(+Inf)", csv_number(inf_v),  "inf");
    check_string("csv_number(-Inf)", csv_number(-inf_v), "-inf");

    // A negative NaN must still print "nan": the sign bit of a NaN carries no
    // meaning and the platform-dependent "-nan" was one of the things this
    // helper exists to remove.
    check_string("csv_number(-NaN)", csv_number(-nan_v), "nan");

    // The spelling must not depend on the requested precision.
    check_string("csv_number(NaN, 9)",   csv_number(nan_v, 9),   "nan");
    check_string("csv_number(-Inf, 10)", csv_number(-inf_v, 10), "-inf");
}

} // namespace

int main()
{
    std::printf("test_cli_helper_formats: B2 output_extension + B3 csv_number\n");

    test_raster_output_extension();
    test_output_extension_overloads();
    test_csv_number_finite();
    test_csv_number_non_finite();

    if (g_failures != 0)
    {
        std::printf("test_cli_helper_formats: %d FAILURE(S)\n", g_failures);
        return 1;
    }

    std::printf("test_cli_helper_formats: all checks passed\n");
    return 0;
}
