// lsdtt_channel_network_ingest.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html

#include "lsdtt_channel_network_ingest.hpp"

#include <stdexcept>
#include <utility>
#include <variant>

#include "../utils/lsdtt_string_utils.hpp"
#include "lsdtt_csv_table.hpp"
#include "lsdtt_data_loader.hpp"

namespace lsdtt
{

lsdtt_ingested_channel_network ingest_channel_network_points(
    const std::string&     path,
    const lsdtt_flow_base& flow)
{
    const int nr = flow.meta().n_rows, nc = flow.meta().n_cols;
    lsdtt_ingested_channel_network out;
    out.mask.assign(static_cast<std::size_t>(nr) * static_cast<std::size_t>(nc), 0u);

    std::vector<std::pair<int,int>> pts;
    const std::string ext = path.size() >= 4
        ? to_lower_ascii(path.substr(path.size() - 4)) : std::string();

    if (ext == ".csv")
    {
        const lsdtt_csv_table t = load_csv_table(path);
        if (!csv_has_column(t, "row") || !csv_has_column(t, "col"))
            throw std::invalid_argument(
                "channel_network_fname " + path + " has no 'row'/'col' "
                "columns. Write it with lsdtt-channel-extraction's "
                "write_channel_network_points.");
        int line = 1;
        for (const auto& r : t.rows)
        {
            ++line;
            pts.emplace_back(
                static_cast<int>(csv_row_double_required(r, "row", path, line)),
                static_cast<int>(csv_row_double_required(r, "col", path, line)));
        }
    }
    else
    {
        lsdtt_data_loader vl;
        vl.open_vector(path);
        const auto names = vl.get_layer_names();
        if (names.empty())
            throw std::invalid_argument(
                "channel_network_fname " + path + " has no layers");
        const lsdtt_vector_layer layer = vl.read_vector_layer(names.front());
        const int i_row = layer.schema.index_of("row");
        const int i_col = layer.schema.index_of("col");
        if (i_row < 0 || i_col < 0)
            throw std::invalid_argument(
                "channel_network_fname " + path + " has no 'row'/'col' "
                "fields. Write it with lsdtt-channel-extraction's "
                "write_channel_network_points.");
        for (const auto& rec : layer.records)
            pts.emplace_back(
                static_cast<int>(std::get<std::int64_t>(rec.field_values[i_row])),
                static_cast<int>(std::get<std::int64_t>(rec.field_values[i_col])));
    }

    if (pts.empty())
        throw std::runtime_error(
            "channel_network_fname " + path + " contained no points");

    // The row/col indices are row-major on the grid the network was
    // extracted from. Reject an out-of-range index rather than wrapping
    // it into a wrong pixel -- a network from a different DEM must not
    // be silently accepted.
    std::int64_t n_bad = 0, n_invalid = 0;
    for (const auto& rc : pts)
    {
        if (rc.first < 0 || rc.first >= nr || rc.second < 0 || rc.second >= nc)
        { ++n_bad; continue; }
        const std::size_t k =
            static_cast<std::size_t>(rc.first) * nc + rc.second;
        if (!flow.is_valid_cell(static_cast<std::int64_t>(k)))
        { ++n_invalid; continue; }
        if (out.mask[k] == 0u) { out.mask[k] = 1u; ++out.n_channel; }
    }
    if (n_bad > 0)
        throw std::runtime_error(
            "channel_network_fname has " + std::to_string(n_bad) +
            " points outside this DEM (" + std::to_string(nr) + "x" +
            std::to_string(nc) + "). It was extracted on a different "
            "grid; row/col are row-major on that grid and cannot be "
            "reused here.");

    if (out.n_channel == 0)
    {
        // Two different causes, and blaming the wrong one sends the user
        // to re-extract a network that was never the problem. If NO
        // point had usable geometry the file is the issue; if points
        // were read but all landed on nodata, the DEM is.
        if (pts.empty())
            throw std::runtime_error(
                "channel_network_fname: no usable POINT geometry was read "
                "from '" + path + "'. Multipart points are supported "
                "and split automatically, so this is either an empty "
                "layer or a geometry type that is not points - check that "
                "the file is the *_channel_network_points output and not "
                "the *_channel_network_lines one.");

        throw std::runtime_error(
            "channel_network_fname: all " + std::to_string(pts.size()) +
            " ingested channel points fell on nodata in this DEM. The "
            "network was extracted on a different (or unmasked) DEM - "
            "re-extract it on the same one.");
    }

    out.n_points    = static_cast<std::int64_t>(pts.size());
    out.n_on_nodata = n_invalid;
    return out;
}

} // namespace lsdtt
