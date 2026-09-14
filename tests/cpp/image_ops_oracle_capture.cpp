// image_ops_oracle_capture.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// lsdtt-test-image-ops-oracle-capture
//
// Freezes the behaviour of the OpenCV image operations (lsdtt::opencv_cpu)
// before OpenCV is removed, so the native replacements can be checked against
// it afterwards without OpenCV installed.
//
// NOT part of the CTest suite: run by hand, once, while OpenCV is linked.
//
//   lsdtt-test-image-ops-oracle-capture <output_dir>
//
// Two OpenCV execution modes are captured:
//
//   reference - cv::ipp::setUseIPP(false) and cv::setUseOptimized(false):
//               OpenCV's own portable code paths (baseline SIMD only). The
//               native port targets THIS mode; its outputs are stored in full.
//   default   - IPP and dispatched SIMD as the build enables them. Stored only
//               as a hash plus the number of pixels that differ from
//               reference, to measure how much OpenCV itself varies by build.
//
// Outputs in <output_dir>:
//   <input>.bin        concatenated records for one input (see manifest)
//   manifest.csv       one row per record: parameters, byte offset/length,
//                      encoding, hashes, default-mode difference count
//   capture_info.txt   OpenCV version, IPP status, CPU features, inputs
//   planted_segments.csv  ground-truth segments of the "lines" input
//
// Encodings:
//   u8       one byte per pixel, row-major (the 8-bit normalised image)
//   bits     binary raster data packed 1 bit per pixel, row-major, LSB first
//            (output mask always equals the input mask, so it is not stored)
// -----------------------------------------------------------------------------

// OpenCV must precede lsdtt_array2d.hpp (see lsdtt_raster_opencv_cpu.hpp).
#include "../../src/cpu/lsdtt_raster_opencv_cpu.hpp"

#include <opencv2/core/utility.hpp>

#include <cmath>
#include <cstdint>
#include <cstdio>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "image_ops_oracle_inputs.hpp"

namespace fs = std::filesystem;
using lsdtt::lsdtt_raster_cpu;
namespace ocv = lsdtt::opencv_cpu;
using namespace lsdtt_oracle;

namespace {

enum class mode { reference, deflt };

void set_mode(mode m)
{
    const bool on = (m == mode::deflt);
    cv::ipp::setUseIPP(on);
    cv::setUseOptimized(on);
}

std::vector<uint8_t> pack_bits(const std::vector<float>& d)
{
    std::vector<uint8_t> out((d.size() + 7) / 8, 0);
    for (std::size_t i = 0; i < d.size(); ++i)
        if (d[i] != 0.0f) out[i / 8] |= static_cast<uint8_t>(1u << (i % 8));
    return out;
}

std::string hex64(std::uint64_t v)
{
    char buf[17];
    std::snprintf(buf, sizeof buf, "%016llx", static_cast<unsigned long long>(v));
    return buf;
}

std::string fmt(double v)
{
    std::ostringstream s;
    s << v;
    return s.str();
}

// One operation to capture: a label, a parameter string, and a function that
// returns the record bytes (and the pixel-level vector used for diff counts).
struct op_result
{
    std::vector<uint8_t> bytes;     // stored encoding
    std::vector<uint8_t> pixels;    // one value per pixel, for diff counting
    std::string          encoding;
};

op_result from_binary_raster(const lsdtt_raster_cpu& r)
{
    const auto& d = r.get_data();
    op_result o;
    o.bytes    = pack_bits(d);
    o.pixels.resize(d.size());
    for (std::size_t i = 0; i < d.size(); ++i) o.pixels[i] = d[i] != 0.0f ? 1 : 0;
    o.encoding = "bits";
    return o;
}

struct op_spec
{
    std::string                                    op;
    std::string                                    params;
    std::function<op_result(const lsdtt_raster_cpu&, int /*opencv_threads*/)> run;
};

std::vector<op_spec> make_ops()
{
    std::vector<op_spec> ops;

    ops.push_back({"u8", "", [](const lsdtt_raster_cpu& r, int)
    {
        ocv::mat8u_params p;
        cv::Mat m = ocv::raster_to_mat8u(r, p, 1);
        op_result o;
        o.bytes.assign(m.data, m.data + m.total());
        o.pixels   = o.bytes;
        o.encoding = "u8";
        return o;
    }});

    ops.push_back({"otsu", "", [](const lsdtt_raster_cpu& r, int t)
    {
        return from_binary_raster(ocv::threshold_otsu(r, 1, t));
    }});

    for (const char* method : {"mean", "gaussian"})
    for (int block : {3, 11, 51})
    for (double C : {-5.0, 0.0, 2.5})
    for (double sigma : {0.0, 1.5})
    {
        const std::string m = method;
        std::ostringstream ps;
        ps << "block=" << block << ";C=" << fmt(C) << ";sigma=" << fmt(sigma);
        ops.push_back({"adaptive_" + m, ps.str(),
            [m, block, C, sigma](const lsdtt_raster_cpu& r, int t)
            {
                return from_binary_raster(m == "mean"
                    ? ocv::threshold_adaptive_mean(r, block, C, sigma, 1, t)
                    : ocv::threshold_adaptive_gaussian(r, block, C, sigma, 1, t));
            }});
    }

    for (double thr : {20.0, 50.0})
    for (int aperture : {3, 5, 7})
    for (double sigma : {0.0, 1.0})
    {
        std::ostringstream ps;
        ps << "threshold=" << fmt(thr) << ";ratio=2;aperture=" << aperture
           << ";sigma=" << fmt(sigma);
        ops.push_back({"canny", ps.str(),
            [thr, aperture, sigma](const lsdtt_raster_cpu& r, int t)
            {
                return from_binary_raster(
                    ocv::detect_edges_canny(r, thr, 2.0, aperture, sigma, 1, t));
            }});
    }

    // Hough: CLI defaults, and a stricter setting.
    struct hough_case { double rho, theta; int thr; double minlen, gap, canny, ratio, sigma; };
    const hough_case hcs[] = {
        {1.0, 1.0, 50, 30.0, 10.0, 50.0, 2.0, 1.0},
        {1.0, 1.0, 80, 60.0,  5.0, 50.0, 2.0, 0.0},
    };
    for (const auto& h : hcs)
    {
        std::ostringstream ps;
        ps << "rho=" << fmt(h.rho) << ";theta_deg=" << fmt(h.theta)
           << ";threshold=" << h.thr << ";min_line_length=" << fmt(h.minlen)
           << ";max_line_gap=" << fmt(h.gap) << ";canny_threshold=" << fmt(h.canny)
           << ";canny_ratio=" << fmt(h.ratio) << ";sigma=" << fmt(h.sigma);
        ops.push_back({"hough", ps.str(), [h](const lsdtt_raster_cpu& r, int t)
        {
            return from_binary_raster(ocv::extract_lines_hough(
                r, h.rho, h.theta, h.thr, h.minlen, h.gap, h.canny, h.ratio,
                h.sigma, 1, t));
        }});
    }

    return ops;
}

} // namespace

int main(int argc, char* argv[])
{
    if (argc != 2)
    {
        std::cerr << "usage: lsdtt-test-image-ops-oracle-capture <output_dir>\n";
        return 2;
    }
    const fs::path out_dir = argv[1];
    fs::create_directories(out_dir);

    const std::vector<oracle_input> inputs = make_oracle_inputs();
    const std::vector<op_spec>      ops    = make_ops();

    std::ofstream manifest(out_dir / "manifest.csv");
    manifest << "input,rows,cols,input_hash,op,params,encoding,offset,nbytes,"
                "reference_hash,default_hash,default_diff_pixels,n_pixels\n";

    long n_records = 0, n_records_differing = 0, total_diff_pixels = 0;

    for (const auto& in : inputs)
    {
        const auto& meta = in.raster.get_meta();
        const std::string input_hash = hex64(hash_raster(in.raster));
        std::ofstream bin(out_dir / (in.name + ".bin"), std::ios::binary);
        std::uint64_t offset = 0;

        for (const auto& op : ops)
        {
            set_mode(mode::reference);
            const op_result ref = op.run(in.raster, 1);

            set_mode(mode::deflt);
            const op_result def = op.run(in.raster, -1);

            long diff = 0;
            for (std::size_t i = 0; i < ref.pixels.size(); ++i)
                if (ref.pixels[i] != def.pixels[i]) ++diff;

            bin.write(reinterpret_cast<const char*>(ref.bytes.data()),
                      static_cast<std::streamsize>(ref.bytes.size()));

            manifest << in.name << ',' << meta.n_rows << ',' << meta.n_cols << ','
                     << input_hash << ',' << op.op << ",\"" << op.params << "\","
                     << ref.encoding << ',' << offset << ',' << ref.bytes.size() << ','
                     << hex64(fnv1a64(ref.bytes.data(), ref.bytes.size())) << ','
                     << hex64(fnv1a64(def.bytes.data(), def.bytes.size())) << ','
                     << diff << ',' << ref.pixels.size() << '\n';

            offset += ref.bytes.size();
            ++n_records;
            if (diff) { ++n_records_differing; total_diff_pixels += diff; }
        }
        std::cout << "captured " << in.name << " (" << ops.size() << " records)\n";
    }

    // -------------------------------------------------------------------------
    // Primitive-level records: OpenCV's own intermediate 8-bit / 16-bit images
    // on a large random image, stored as hashes (primitives.csv). These pin the
    // rounding and accumulation order of each primitive, which thresholded
    // outputs are too coarse to reveal.
    // -------------------------------------------------------------------------
    {
        const int R = k_primitive_rows, Cc = k_primitive_cols;
        std::vector<std::uint8_t> u8 = make_primitive_u8();
        std::vector<float>        f32 = make_primitive_f32();
        const cv::Mat mu8(R, Cc, CV_8UC1, u8.data());
        const cv::Mat mf32(R, Cc, CV_32FC1, f32.data());

        std::ofstream prim(out_dir / "primitives.csv");
        prim << "primitive,params,rows,cols,input_hash,dtype,reference_hash,default_hash\n";
        const std::string in_u8_hash  = hex64(fnv1a64(u8.data(), u8.size()));
        const std::string in_f32_hash = hex64(fnv1a64(reinterpret_cast<const unsigned char*>(f32.data()), f32.size() * sizeof(float)));

        auto hash_mat = [](const cv::Mat& m)
        {
            return hex64(fnv1a64(m.data, m.total() * m.elemSize()));
        };
        auto emit = [&](const std::string& name, const std::string& params, const std::string& in_hash,
                        const char* dtype, const std::function<cv::Mat()>& run)
        {
            set_mode(mode::reference);
            const std::string ref = hash_mat(run());
            set_mode(mode::deflt);
            const std::string def = hash_mat(run());
            prim << name << ",\"" << params << "\"," << R << ',' << Cc << ',' << in_hash << ','
                 << dtype << ',' << ref << ',' << def << '\n';
        };

        for (int block : {11, 21, 31, 51})
            emit("gaussian_mean_u8", "block=" + std::to_string(block), in_u8_hash, "u8", [&]()
            {
                cv::Mat f, g, out;
                mu8.convertTo(f, CV_32F);
                cv::GaussianBlur(f, g, cv::Size(block, block), 0, 0, cv::BORDER_REPLICATE | cv::BORDER_ISOLATED);
                g.convertTo(out, CV_8U);
                return out;
            });
        for (double sigma : {0.7, 1.0, 1.5, 2.3, 3.1})
            emit("gaussian_blur_u8", "sigma=" + fmt(sigma), in_u8_hash, "u8", [&]()
            {
                const int ksize = 2 * static_cast<int>(std::ceil(3.0 * sigma)) + 1;
                cv::Mat out;
                cv::GaussianBlur(mu8, out, cv::Size(ksize, ksize), sigma, sigma);
                return out;
            });
        for (int block : {3, 9, 15, 17, 31, 51})
            emit("box_mean_u8", "block=" + std::to_string(block), in_u8_hash, "u8", [&]()
            {
                cv::Mat out;
                cv::boxFilter(mu8, out, CV_8U, cv::Size(block, block), cv::Point(-1, -1), true,
                              cv::BORDER_REPLICATE | cv::BORDER_ISOLATED);
                return out;
            });
        emit("normalize_u8", "", in_f32_hash, "u8", [&]()
        {
            cv::Mat out;
            cv::normalize(mf32, out, 0, 255, cv::NORM_MINMAX, CV_8UC1);
            return out;
        });
        for (int ap : {3, 5, 7})
            for (const auto& th : {std::pair<double, double>{20.0, 40.0}, {50.0, 100.0}, {80.0, 300.0}, {120.0, 60.0}})
                emit("canny_u8", "aperture=" + std::to_string(ap) + ";low=" + fmt(th.first) + ";high=" + fmt(th.second),
                     in_u8_hash, "u8", [&]()
                     {
                         cv::Mat out;
                         cv::Canny(mu8, out, th.first, th.second, ap);
                         return out;
                     });
        // cv::line (LINE_8, thickness 1) for 2000 segments whose endpoints are
        // drawn from splitmix64 inside the image; covers every octant.
        emit("line8", "segments=2000", "none", "u8", [&]()
        {
            cv::Mat out(R, Cc, CV_8UC1, cv::Scalar(0));
            for (std::uint64_t i = 0; i < 2000; ++i)
            {
                const int x0 = static_cast<int>(splitmix64(0x11 + 4 * i) % static_cast<std::uint64_t>(Cc));
                const int y0 = static_cast<int>(splitmix64(0x12 + 4 * i) % static_cast<std::uint64_t>(R));
                const int x1 = static_cast<int>(splitmix64(0x13 + 4 * i) % static_cast<std::uint64_t>(Cc));
                const int y1 = static_cast<int>(splitmix64(0x14 + 4 * i) % static_cast<std::uint64_t>(R));
                // Accumulate counts so overlapping segments are also checked.
                cv::Mat one(R, Cc, CV_8UC1, cv::Scalar(0));
                cv::line(one, cv::Point(x0, y0), cv::Point(x1, y1), cv::Scalar(1), 1, cv::LINE_8);
                out += one;
            }
            return out;
        });
        for (int ap : {3, 5, 7})
            for (int dir = 0; dir < 2; ++dir)
                emit("sobel_s16", "aperture=" + std::to_string(ap) + ";dx=" + std::to_string(1 - dir),
                     in_u8_hash, "s16", [&]()
                     {
                         cv::Mat out;
                         cv::Sobel(mu8, out, CV_16S, 1 - dir, dir, ap, ap == 7 ? 1.0 / 16.0 : 1.0, 0,
                                   cv::BORDER_REPLICATE);
                         return out;
                     });
    }

    // -------------------------------------------------------------------------
    // Hough segment lists: cv::HoughLinesP on the Canny edges of every input
    // without nodata (OpenCV's own 8-bit conversion, pre-blur and Canny). The
    // native Canny is identical on these inputs, so a test can replay cv::RNG
    // through the native Hough and require the same segments, in order.
    // -------------------------------------------------------------------------
    {
        set_mode(mode::reference);
        std::ofstream hs(out_dir / "hough_opencv_segments.csv");
        hs << "input,params,n_segments,index,x0,y0,x1,y1\n";
        for (const auto& in : inputs)
        {
            if (in.name.find("_full") == std::string::npos) continue;
            ocv::mat8u_params prm;
            const cv::Mat base = ocv::raster_to_mat8u(in.raster, prm, 1);
            for (double sigma : {0.0, 1.0})
                for (double ct : {30.0, 50.0})
                    for (int thr : {20, 50, 80})
                    {
                        cv::Mat img = base.clone();
                        if (sigma > 0.0)
                        {
                            const int ksize = 2 * static_cast<int>(std::ceil(3.0 * sigma)) + 1;
                            cv::GaussianBlur(img, img, cv::Size(ksize, ksize), sigma, sigma);
                        }
                        cv::Mat edges;
                        cv::Canny(img, edges, ct, ct * 2.0, 3);
                        std::vector<cv::Vec4i> lines;
                        cv::HoughLinesP(edges, lines, 1.0, CV_PI / 180.0, thr, 30.0, 10.0);
                        std::ostringstream ps;
                        ps << "sigma=" << fmt(sigma) << ";canny_threshold=" << fmt(ct) << ";threshold=" << thr
                           << ";min_line_length=30;max_line_gap=10";
                        if (lines.empty())
                            hs << in.name << ",\"" << ps.str() << "\",0,-1,0,0,0,0\n";
                        for (std::size_t i = 0; i < lines.size(); ++i)
                            hs << in.name << ",\"" << ps.str() << "\"," << lines.size() << ',' << i << ','
                               << lines[i][0] << ',' << lines[i][1] << ',' << lines[i][2] << ',' << lines[i][3] << '\n';
                    }
        }
    }

    {
        std::ofstream seg(out_dir / "planted_segments.csv");
        seg << "c0,r0,c1,r1\n";
        for (const auto& s : lines_planted_segments())
            seg << s.c0 << ',' << s.r0 << ',' << s.c1 << ',' << s.r1 << '\n';
    }

    {
        set_mode(mode::reference);
        const bool ref_ipp = cv::ipp::useIPP(), ref_opt = cv::useOptimized();
        set_mode(mode::deflt);
        const bool ipp_available = cv::ipp::useIPP();
        const bool def_opt = cv::useOptimized();
        std::ofstream info(out_dir / "capture_info.txt");
        const std::time_t now = std::time(nullptr);
        char date[32];
        std::strftime(date, sizeof date, "%Y-%m-%dT%H:%M:%SZ", std::gmtime(&now));
        info << "captured_utc: " << date << "\n"
             << "opencv_version: " << cv::getVersionString() << "\n"
             << "ipp_available_in_default_mode: " << (ipp_available ? "yes" : "no") << "\n"
             << "ipp_version: " << cv::ipp::getIppVersion() << "\n"
             << "reference_mode_state: useIPP=" << ref_ipp << " useOptimized=" << ref_opt << "\n"
             << "default_mode_state: useIPP=" << ipp_available << " useOptimized=" << def_opt << "\n"
             << "cpu_features: " << cv::getCPUFeaturesLine() << "\n"
             << "records: " << n_records << "\n"
             << "records_where_default_differs_from_reference: " << n_records_differing << "\n"
             << "total_default_diff_pixels: " << total_diff_pixels << "\n"
             << "inputs: generated by tests/cpp/image_ops_oracle_inputs.hpp\n";
    }

    std::cout << "records: " << n_records
              << "; default mode differs from reference in " << n_records_differing
              << " records (" << total_diff_pixels << " pixels)\n";
    return 0;
}
