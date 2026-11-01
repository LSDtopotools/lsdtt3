// lsdtt_surface_method_comparison_cpu.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html

#include "lsdtt_surface_method_comparison_cpu.hpp"

#include <limits>
#include <stdexcept>

namespace lsdtt
{

namespace
{

void check_inputs(const std::vector<method_class_map>& methods,
                  const std::vector<std::uint8_t>&     valid,
                  const char*                          fn)
{
    for (const method_class_map& m : methods)
    {
        if (m.cls.size() != valid.size())
            throw std::invalid_argument(std::string(fn) + ": method '" + m.name +
                                        "' is not on the comparison grid");
        for (const std::uint8_t c : m.cls)
            if (c > 2)
                throw std::invalid_argument(std::string(fn) + ": method '" +
                                            m.name + "' has a class above 2");
    }
}

double ratio(std::int64_t num, std::int64_t den)
{
    return den > 0 ? static_cast<double>(num) / static_cast<double>(den)
                   : std::numeric_limits<double>::quiet_NaN();
}

} // namespace

std::vector<method_pair_comparison> compare_surface_classes(
    const std::vector<method_class_map>& methods,
    const std::vector<std::uint8_t>&     valid)
{
    check_inputs(methods, valid, "compare_surface_classes");
    std::vector<method_pair_comparison> out;
    for (std::size_t a = 0; a < methods.size(); ++a)
        for (std::size_t b = a + 1; b < methods.size(); ++b)
        {
            method_pair_comparison pc;
            pc.a = a;
            pc.b = b;
            const std::vector<std::uint8_t>& ca = methods[a].cls;
            const std::vector<std::uint8_t>& cb = methods[b].cls;
            for (std::size_t i = 0; i < valid.size(); ++i)
                if (valid[i]) ++pc.confusion[ca[i]][cb[i]];
            auto iou = [&](int k)
            {
                std::int64_t inter = pc.confusion[k][k], uni = 0;
                for (int j = 0; j < 3; ++j)
                    uni += pc.confusion[k][j] + pc.confusion[j][k];
                uni -= pc.confusion[k][k];
                return ratio(inter, uni);
            };
            pc.floodplain_iou = iou(1);
            pc.terrace_iou    = iou(2);
            std::int64_t either = 0;
            for (int i = 0; i < 3; ++i)
                for (int j = 0; j < 3; ++j)
                    if (i != 0 || j != 0) either += pc.confusion[i][j];
            pc.classified_agreement =
                ratio(pc.confusion[1][1] + pc.confusion[2][2], either);
            out.push_back(pc);
        }
    return out;
}

mapped_surface_score score_against_mapped(const method_class_map&          method,
                                          const std::vector<std::uint8_t>& mapped,
                                          const std::vector<std::uint8_t>& scored)
{
    check_inputs({method}, scored, "score_against_mapped");
    if (mapped.size() != scored.size())
        throw std::invalid_argument(
            "score_against_mapped: the mapped classes are not on the grid");
    mapped_surface_score sc;
    for (std::size_t i = 0; i < scored.size(); ++i)
    {
        if (!scored[i]) continue;
        if (mapped[i] > 2)
            throw std::invalid_argument(
                "score_against_mapped: a mapped class above 2");
        ++sc.confusion[mapped[i]][method.cls[i]];
        ++sc.n_scored;
    }
    for (int k = 1; k <= 2; ++k)
    {
        std::int64_t by_method = 0, by_map = 0;
        for (int j = 0; j < 3; ++j)
        {
            by_method += sc.confusion[j][k];
            by_map    += sc.confusion[k][j];
        }
        const std::int64_t tp = sc.confusion[k][k];
        sc.precision[k] = ratio(tp, by_method);
        sc.recall[k]    = ratio(tp, by_map);
        sc.iou[k]       = ratio(tp, by_method + by_map - tp);
    }
    return sc;
}

std::vector<std::uint8_t> surface_class_agreement_bits(
    const std::vector<method_class_map>& methods,
    const std::vector<std::uint8_t>&     valid,
    std::uint8_t                         cls)
{
    check_inputs(methods, valid, "surface_class_agreement_bits");
    if (methods.size() > 8)
        throw std::invalid_argument(
            "surface_class_agreement_bits: at most 8 methods fit one byte");
    std::vector<std::uint8_t> bits(valid.size(), 0);
    for (std::size_t k = 0; k < methods.size(); ++k)
        for (std::size_t i = 0; i < valid.size(); ++i)
            if (valid[i] && methods[k].cls[i] == cls)
                bits[i] = static_cast<std::uint8_t>(bits[i] | (1u << k));
    return bits;
}

} // namespace lsdtt
