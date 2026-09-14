// test_image_ops.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// Tests for src/cpu/lsdtt_raster_image_ops_cpu (OpenCV removal plan, Phase 2).
// No OpenCV needed: OpenCV's behaviour is read from the frozen oracle in
// tests/expected/image_ops/ (written by lsdtt-test-image-ops-oracle-capture).
//
//   lsdtt-test-image-ops <path to tests/expected/image_ops>
//
//   O  oracle: native output vs OpenCV, pixel for pixel. Exact for every input
//      without nodata. For inputs with nodata the native code deliberately
//      differs (valid-pixel normalisation, plan D3); the difference is reported
//      as a measurement, and I4/I5 check the corrected behaviour directly.
//   P  primitives: the native primitives' own 8-bit / 16-bit images vs OpenCV
//      on a large random image (hash comparison; every rounding decision counts)
//   I  invariants of the native implementation (thread count, nodata, offset,
//      the valid_pixels normalisation)
// -----------------------------------------------------------------------------

#include "../../src/cpu/lsdtt_raster_image_ops_cpu.hpp"
#include "../../src/utils/lsdtt_image_filters.hpp"
#include "image_ops_oracle_inputs.hpp"

#include <cmath>
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

using lsdtt::lsdtt_raster_cpu;
namespace ops = lsdtt::image_ops_cpu;
using namespace lsdtt_oracle;

namespace {

int g_n_checks = 0;
int g_n_failed = 0;

void check(const std::string& name, bool ok, const std::string& detail = "")
{
    ++g_n_checks;
    if (!ok) ++g_n_failed;
    std::cout << (ok ? "PASS" : "FAIL") << "  " << name;
    if (!detail.empty()) std::cout << "   [" << detail << "]";
    std::cout << "\n";
}

void info(const std::string& text) { std::cout << "INFO  " << text << "\n"; }

std::string hex64(std::uint64_t v)
{
    char buf[17];
    std::snprintf(buf, sizeof buf, "%016llx", static_cast<unsigned long long>(v));
    return buf;
}

// --- manifest ---------------------------------------------------------------

struct record
{
    std::string input, input_hash, op, params, encoding;
    std::uint64_t offset = 0, nbytes = 0;
    long n_pixels = 0;
};

std::vector<std::string> split_csv_line(const std::string& line)
{
    std::vector<std::string> out;
    std::string cur;
    bool quoted = false;
    for (char ch : line)
    {
        if (ch == '"') { quoted = !quoted; continue; }
        if (ch == ',' && !quoted) { out.push_back(cur); cur.clear(); continue; }
        cur.push_back(ch);
    }
    out.push_back(cur);
    return out;
}

std::vector<record> read_manifest(const std::string& dir)
{
    std::ifstream in(dir + "/manifest.csv");
    std::vector<record> out;
    std::string line;
    std::getline(in, line);   // header
    while (std::getline(in, line))
    {
        if (line.empty()) continue;
        const auto f = split_csv_line(line);
        // input,rows,cols,input_hash,op,params,encoding,offset,nbytes,
        // reference_hash,default_hash,default_diff_pixels,n_pixels
        if (f.size() < 13) continue;
        record r;
        r.input = f[0]; r.input_hash = f[3]; r.op = f[4]; r.params = f[5]; r.encoding = f[6];
        r.offset = std::stoull(f[7]); r.nbytes = std::stoull(f[8]); r.n_pixels = std::stol(f[12]);
        out.push_back(r);
    }
    return out;
}

std::map<std::string, double> parse_params(const std::string& p)
{
    std::map<std::string, double> out;
    std::stringstream ss(p);
    std::string item;
    while (std::getline(ss, item, ';'))
    {
        const auto eq = item.find('=');
        if (eq != std::string::npos) out[item.substr(0, eq)] = std::stod(item.substr(eq + 1));
    }
    return out;
}

// Decode a record into one value per pixel (u8 value, or 0/1 for bits).
std::vector<std::uint8_t> load_record(const std::string& dir, const record& r)
{
    std::ifstream in(dir + "/" + r.input + ".bin", std::ios::binary);
    std::vector<char> bytes(r.nbytes);
    in.seekg(static_cast<std::streamoff>(r.offset));
    in.read(bytes.data(), static_cast<std::streamsize>(r.nbytes));
    std::vector<std::uint8_t> px(static_cast<std::size_t>(r.n_pixels), 0);
    if (r.encoding == "u8")
    {
        for (long i = 0; i < r.n_pixels; ++i) px[static_cast<std::size_t>(i)] = static_cast<std::uint8_t>(bytes[static_cast<std::size_t>(i)]);
    }
    else
    {
        for (long i = 0; i < r.n_pixels; ++i)
            px[static_cast<std::size_t>(i)] = (static_cast<std::uint8_t>(bytes[static_cast<std::size_t>(i / 8)]) >> (i % 8)) & 1;
    }
    return px;
}

std::vector<std::uint8_t> binary_pixels(const lsdtt_raster_cpu& r)
{
    const auto& d = r.get_data();
    std::vector<std::uint8_t> px(d.size());
    for (std::size_t i = 0; i < d.size(); ++i) px[i] = d[i] != 0.0f ? 1 : 0;
    return px;
}

long count_ones(const lsdtt_raster_cpu& r)
{
    long n = 0;
    for (float v : r.get_data()) if (v != 0.0f) ++n;
    return n;
}

bool same_raster(const lsdtt_raster_cpu& a, const lsdtt_raster_cpu& b)
{
    return a.get_data() == b.get_data() && a.get_mask() == b.get_mask();
}

// Run the native op for a record with the given normalisation.
bool run_native(const record& r, const lsdtt_raster_cpu& in,
                std::vector<std::uint8_t>& px, int n_threads = 1)
{
    const auto p = parse_params(r.params);
    if (r.op == "u8")
    {
        ops::u8_params params;
        px = ops::raster_to_u8(in, params, n_threads);
        return true;
    }
    if (r.op == "otsu")
    {
        px = binary_pixels(ops::threshold_otsu(in, n_threads));
        return true;
    }
    if (r.op == "adaptive_mean" || r.op == "adaptive_gaussian")
    {
        const int block = static_cast<int>(p.at("block"));
        const double C = p.at("C"), sigma = p.at("sigma");
        px = binary_pixels(r.op == "adaptive_mean"
            ? ops::threshold_adaptive_mean(in, block, C, sigma, n_threads)
            : ops::threshold_adaptive_gaussian(in, block, C, sigma, n_threads));
        return true;
    }
    return false;   // not ported in this phase (canny, hough)
}


// =============================================================================
// O - oracle
// =============================================================================
void test_oracle(const std::string& dir)
{
    std::cout << "--- O: native vs frozen OpenCV output ---\n";

    const std::vector<record> recs = read_manifest(dir);
    check("O0 manifest read (832 records expected)", recs.size() == 832,
          std::to_string(recs.size()) + " records");

    std::map<std::string, lsdtt_raster_cpu> inputs;
    bool hashes_ok = true;
    for (auto& in : make_oracle_inputs()) inputs.emplace(in.name, std::move(in.raster));
    for (const auto& r : recs)
    {
        const auto it = inputs.find(r.input);
        if (it == inputs.end() || hex64(hash_raster(it->second)) != r.input_hash) { hashes_ok = false; break; }
    }
    check("O1 regenerated inputs match the hashes recorded at capture", hashes_ok);
    if (!hashes_ok) return;

    struct tally { long records = 0, bad_records = 0, bad_pixels = 0; std::string first_bad; };
    std::map<std::string, tally> exact_by_op, nodata_by_op;
    std::map<std::string, long> nodata_valid_pixels;

    for (const auto& r : recs)
    {
        std::vector<std::uint8_t> got;
        const lsdtt_raster_cpu& in = inputs.at(r.input);
        if (!run_native(r, in, got)) continue;
        const std::vector<std::uint8_t> want = load_record(dir, r);
        long diff = 0;
        long valid = 0;
        for (std::size_t i = 0; i < want.size(); ++i)
        {
            if (in.get_mask()[i] != 0) ++valid;
            if (got[i] != want[i]) ++diff;
        }
        const bool has_nodata = valid < static_cast<long>(want.size());
        tally& t = has_nodata ? nodata_by_op[r.op] : exact_by_op[r.op];
        if (has_nodata) nodata_valid_pixels[r.op] += valid;
        ++t.records;
        if (diff)
        {
            if (t.bad_records++ == 0) t.first_bad = r.input + " " + r.params + " (" + std::to_string(diff) + " px)";
            t.bad_pixels += diff;
        }
    }

    for (const char* op : {"u8", "otsu", "adaptive_mean", "adaptive_gaussian"})
    {
        const tally& t = exact_by_op[op];
        std::ostringstream d;
        d << t.records << " records, " << t.bad_records << " differ, " << t.bad_pixels << " pixels";
        if (t.bad_records) d << "; first: " << t.first_bad;
        check(std::string("O2 ") + op + " identical to OpenCV on every input without nodata",
              t.records > 0 && t.bad_records == 0, d.str());
    }
    for (const char* op : {"u8", "otsu", "adaptive_mean", "adaptive_gaussian"})
    {
        const tally& t = nodata_by_op[op];
        std::ostringstream d;
        d << op << " on inputs with nodata (valid-pixel normalisation vs OpenCV): "
          << t.bad_records << " of " << t.records << " records differ, "
          << t.bad_pixels << " pixels (all pixels, incl. nodata, for u8; valid pixels total "
          << nodata_valid_pixels[op] << ")";
        info(d.str());
    }
}


// =============================================================================
// P - primitives
// =============================================================================
void test_primitives(const std::string& dir)
{
    std::cout << "\n--- P: native primitives vs OpenCV on a 512x487 random image ---\n";
    namespace img = lsdtt::image;

    const int R = k_primitive_rows, C = k_primitive_cols;
    const std::vector<std::uint8_t> u8 = make_primitive_u8();
    const std::vector<float>        f32 = make_primitive_f32();
    const std::string in_u8_hash  = hex64(fnv1a64(u8.data(), u8.size()));
    const std::string in_f32_hash = hex64(fnv1a64(reinterpret_cast<const unsigned char*>(f32.data()), f32.size() * sizeof(float)));

    std::ifstream in(dir + "/primitives.csv");
    std::string line;
    std::getline(in, line);
    int n_rows = 0;
    while (std::getline(in, line))
    {
        if (line.empty()) continue;
        const auto f = split_csv_line(line);
        if (f.size() < 8) continue;
        ++n_rows;
        const std::string name = f[0], params = f[1], input_hash = f[4], want = f[6];
        const auto p = parse_params(params);

        std::string got;
        if (name == "normalize_u8")
        {
            if (input_hash != in_f32_hash) { check("P " + name + " input hash", false); continue; }
            lsdtt::lsdtt_raster_meta meta = oracle_meta(R, C);
            lsdtt_raster_cpu r(std::vector<float>(f32), std::vector<std::uint8_t>(f32.size(), 1), meta);
            ops::u8_params prm;
            const auto out = ops::raster_to_u8(r, prm, 1);
            got = hex64(fnv1a64(out.data(), out.size()));
        }
        else
        {
            if (input_hash != in_u8_hash) { check("P " + name + " input hash", false); continue; }
            if (name == "gaussian_mean_u8")
            {
                const int block = static_cast<int>(p.at("block"));
                std::vector<float> fi(u8.begin(), u8.end()), fo;
                const auto k = img::gaussian_kernel_opencv(block, 0.0);
                img::separable_filter_f32(fi, R, C, k, k, img::border_mode::replicate, fo);
                std::vector<std::uint8_t> out(fo.size());
                for (std::size_t i = 0; i < fo.size(); ++i)
                {
                    // saturate_cast<uchar>(float): nearest, ties to even
                    float v = fo[i];
                    v = v < 0.0f ? 0.0f : (v > 255.0f ? 255.0f : v);
                    out[i] = static_cast<std::uint8_t>(std::nearbyint(v));
                }
                got = hex64(fnv1a64(out.data(), out.size()));
            }
            else if (name == "gaussian_blur_u8")
            {
                const double sigma = p.at("sigma");
                const int ksize = 2 * static_cast<int>(std::ceil(3.0 * sigma)) + 1;
                std::vector<std::uint8_t> out;
                img::gaussian_blur_u8(u8, R, C, ksize, sigma, img::border_mode::reflect_101, out);
                got = hex64(fnv1a64(out.data(), out.size()));
            }
            else if (name == "box_mean_u8")
            {
                std::vector<std::uint8_t> out;
                img::box_mean_u8(u8, R, C, static_cast<int>(p.at("block")), img::border_mode::replicate, out);
                got = hex64(fnv1a64(out.data(), out.size()));
            }
            else if (name == "sobel_s16")
            {
                const int ap = static_cast<int>(p.at("aperture"));
                const int dx = static_cast<int>(p.at("dx"));
                std::vector<std::int16_t> out;
                img::sobel_s16(u8, R, C, dx, 1 - dx, ap, ap == 7 ? 16 : 1, img::border_mode::replicate, out);
                got = hex64(fnv1a64(reinterpret_cast<const unsigned char*>(out.data()), out.size() * sizeof(std::int16_t)));
            }
            else
            {
                continue;
            }
        }
        check("P " + name + " " + params + " identical to OpenCV (reference mode)", got == want,
              got == want ? "" : "native " + got + " vs OpenCV " + want);
    }
    check("P0 primitives.csv read (22 rows expected)", n_rows == 22, std::to_string(n_rows) + " rows");
}


// =============================================================================
// I - invariants
// =============================================================================
void test_invariants()
{
    std::cout << "\n--- I: invariants ---\n";

    std::map<std::string, lsdtt_raster_cpu> in;
    for (auto& x : make_oracle_inputs()) in.emplace(x.name, std::move(x.raster));

    // I1: bit-identical across thread counts.
    {
        bool same = true;
        for (const char* name : {"peaks_nodata", "signed_nodata", "lines_full", "entropy_nodata"})
        {
            const auto& r = in.at(name);
            const auto o1 = ops::threshold_otsu(r, 1);
            const auto m1 = ops::threshold_adaptive_mean(r, 11, 0.0, 1.5, 1);
            const auto g1 = ops::threshold_adaptive_gaussian(r, 51, -2.0, 0.0, 1);
            for (int t : {2, 8})
            {
                if (!same_raster(o1, ops::threshold_otsu(r, t))) same = false;
                if (!same_raster(m1, ops::threshold_adaptive_mean(r, 11, 0.0, 1.5, t))) same = false;
                if (!same_raster(g1, ops::threshold_adaptive_gaussian(r, 51, -2.0, 0.0, t))) same = false;
            }
        }
        check("I1 outputs bit-identical for n_threads 1, 2, 8", same);
    }

    // I2: nodata pixels output 0 with mask 0; output mask equals input mask.
    {
        bool ok = true;
        const auto& r = in.at("ramp_nodata");
        for (const auto& out : {ops::threshold_otsu(r, 1),
                                ops::threshold_adaptive_mean(r, 11, -5.0, 0.0, 1),
                                ops::threshold_adaptive_gaussian(r, 11, 2.5, 1.5, 1)})
        {
            if (out.get_mask() != r.get_mask()) ok = false;
            for (std::size_t i = 0; i < out.get_data().size(); ++i)
                if (r.get_mask()[i] == 0 && out.get_data()[i] != 0.0f) ok = false;
        }
        check("I2 output mask equals input mask; nodata pixels are 0", ok);
    }

    // I3: a constant offset does not change any output in valid_pixels mode.
    // Constructed so the scaling is exact in float: integer data spanning
    // exactly 0..255 (scale = 1) and an offset of 256 (shift = -256).
    {
        const int n = 64;
        lsdtt::lsdtt_raster_meta meta = oracle_meta(n);
        std::vector<float> base(static_cast<std::size_t>(n) * n), shifted(base.size());
        std::vector<std::uint8_t> mask(base.size(), 1);
        for (int r = 0; r < n; ++r)
            for (int c = 0; c < n; ++c)
            {
                const std::size_t i = static_cast<std::size_t>(r) * n + c;
                double v = std::floor(255.0 * hash_uniform(77, r, c));
                if (r == 0 && c == 0) v = 0.0;
                if (r == 0 && c == 1) v = 255.0;
                base[i] = static_cast<float>(v);
                shifted[i] = static_cast<float>(v + 256.0);
                if (oracle_is_nodata(n, r, c)) mask[i] = (r == 0) ? 1 : 0;   // keep the 0 and 255 anchors valid
            }
        std::vector<std::uint8_t> mask2 = mask;
        lsdtt_raster_cpu a(std::vector<float>(base), std::move(mask), meta);
        lsdtt_raster_cpu b(std::vector<float>(shifted), std::move(mask2), meta);
        const bool ok = same_raster(ops::threshold_otsu(a, 1), ops::threshold_otsu(b, 1)) &&
                        same_raster(ops::threshold_adaptive_mean(a, 11, 0.0, 1.5, 1),
                                    ops::threshold_adaptive_mean(b, 11, 0.0, 1.5, 1)) &&
                        same_raster(ops::threshold_adaptive_gaussian(a, 51, -2.0, 0.0, 1),
                                    ops::threshold_adaptive_gaussian(b, 51, -2.0, 0.0, 1));
        check("I3 adding 256 to every value (exact in float) changes no output, with nodata present", ok);
    }

    // I4: valid_pixels 8-bit image on valid pixels equals that of the same data
    // with nodata removed (cropped to the valid interior of a border-only mask).
    {
        const int n = 128, border = 6, m = n - 2 * border;
        const auto& full = in.at("ramp_full");
        lsdtt::lsdtt_raster_meta meta = oracle_meta(n);
        std::vector<float> d = full.get_data();
        std::vector<std::uint8_t> mk(d.size(), 1);
        for (int r = 0; r < n; ++r)
            for (int c = 0; c < n; ++c)
                if (r < border || c < border || r >= n - border || c >= n - border)
                { mk[static_cast<std::size_t>(r) * n + c] = 0; d[static_cast<std::size_t>(r) * n + c] = -9999.0f; }
        lsdtt_raster_cpu bordered(std::move(d), std::move(mk), meta);

        lsdtt::lsdtt_raster_meta cmeta = oracle_meta(m);
        std::vector<float> cd(static_cast<std::size_t>(m) * m);
        for (int r = 0; r < m; ++r)
            for (int c = 0; c < m; ++c)
                cd[static_cast<std::size_t>(r) * m + c] = full.get_data()[static_cast<std::size_t>(r + border) * n + c + border];
        lsdtt_raster_cpu cropped(std::move(cd), std::vector<std::uint8_t>(static_cast<std::size_t>(m) * m, 1), cmeta);

        ops::u8_params p1, p2;
        const auto u_b = ops::raster_to_u8(bordered, p1, 1);
        const auto u_c = ops::raster_to_u8(cropped,  p2, 1);
        bool same = true;
        for (int r = 0; r < m; ++r)
            for (int c = 0; c < m; ++c)
                if (u_b[static_cast<std::size_t>(r + border) * n + c + border] != u_c[static_cast<std::size_t>(r) * m + c]) same = false;
        check("I4 8-bit image on valid pixels equals that of the nodata-free crop", same);

        // Otsu: same foreground on valid pixels as the crop.
        const auto o_b = ops::threshold_otsu(bordered, 1);
        const auto o_c = ops::threshold_otsu(cropped,  1);
        bool same_o = true;
        for (int r = 0; r < m; ++r)
            for (int c = 0; c < m; ++c)
                if (o_b.get_data()[static_cast<std::size_t>(r + border) * n + c + border] != o_c.get_data()[static_cast<std::size_t>(r) * m + c]) same_o = false;
        check("I5 Otsu on a raster with a nodata border equals Otsu on the crop", same_o);

    }

    // Measured: Otsu foreground on the oracle _nodata inputs.
    for (const char* name : {"ramp_nodata", "peaks_nodata", "noise_nodata", "steps_nodata"})
    {
        const auto& r = in.at(name);
        long valid = 0;
        for (auto v : r.get_mask()) if (v) ++valid;
        std::ostringstream s;
        s << "Otsu foreground on " << name << ": " << count_ones(ops::threshold_otsu(r, 1))
          << " of " << valid << " valid pixels";
        info(s.str());
    }
}

} // namespace

int main(int argc, char* argv[])
{
    if (argc != 2)
    {
        std::cerr << "usage: lsdtt-test-image-ops <tests/expected/image_ops>\n";
        return 2;
    }
    test_oracle(argv[1]);
    test_primitives(argv[1]);
    test_invariants();
    std::cout << "\n" << (g_n_checks - g_n_failed) << "/" << g_n_checks << " checks passed\n";
    return g_n_failed == 0 ? 0 : 1;
}
