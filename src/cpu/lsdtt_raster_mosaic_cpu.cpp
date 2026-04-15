// lsdtt_raster_mosaic_cpu.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html

#include "lsdtt_raster_mosaic_cpu.hpp"

#include "../core/lsdtt_data_loader.hpp"
#include "../core/lsdtt_raster_meta.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include <ogr_spatialref.h>

namespace lsdtt {

// =============================================================================
// raster_mosaic
// =============================================================================

lsdtt_raster_cpu raster_mosaic(
    const std::vector<std::string>& paths)
{
    if (paths.empty())
        throw std::invalid_argument("raster_mosaic: path list is empty.");

    // -------------------------------------------------------------------------
    // Pass 1 — load metadata, validate, collect accepted entries.
    //
    // The first raster that can be opened sets the reference CRS, cell size,
    // and grid origin.  All subsequent rasters are checked against it.
    // -------------------------------------------------------------------------
    struct Entry
    {
        std::string       path;
        lsdtt_raster_meta meta;
    };

    std::vector<Entry> accepted;
    OGRSpatialReference ref_srs;
    double ref_cs   = 0.0;
    double ref_xmin = 0.0;
    double ref_ymin = 0.0;

    for (const auto& p : paths)
    {
        // --- Open and read metadata ---
        lsdtt_raster_meta meta;
        try
        {
            lsdtt_data_loader ld;
            ld.open_raster(p);
            meta = ld.get_meta();
        }
        catch (const std::exception& e)
        {
            std::cerr << "raster_mosaic: skipping '" << p
                      << "' — cannot open: " << e.what() << "\n";
            continue;
        }

        if (!meta.is_valid())
        {
            std::cerr << "raster_mosaic: skipping '" << p
                      << "' — invalid metadata.\n";
            continue;
        }

        // --- First accepted raster sets the reference ---
        if (accepted.empty())
        {
            if (ref_srs.importFromWkt(meta.crs_wkt.c_str()) != OGRERR_NONE)
            {
                std::cerr << "raster_mosaic: skipping '" << p
                          << "' — cannot parse CRS WKT.\n";
                continue;
            }
            ref_cs   = meta.cell_size;
            ref_xmin = meta.x_min;
            ref_ymin = meta.y_min;
            accepted.push_back({p, meta});
            std::cout << "raster_mosaic: reference raster: " << p
                      << " (cell_size=" << ref_cs << ")\n";
            continue;
        }

        // --- CRS check ---
        OGRSpatialReference this_srs;
        if (this_srs.importFromWkt(meta.crs_wkt.c_str()) != OGRERR_NONE)
        {
            std::cerr << "raster_mosaic: skipping '" << p
                      << "' — cannot parse CRS WKT.\n";
            continue;
        }
        if (!ref_srs.IsSame(&this_srs))
        {
            std::cerr << "raster_mosaic: skipping '" << p
                      << "' — CRS differs from reference.\n";
            continue;
        }

        // --- Cell size check (0.001 % relative tolerance) ---
        const double cs_rel_diff =
            std::abs(meta.cell_size - ref_cs) / ref_cs;
        if (cs_rel_diff > 1e-5)
        {
            std::cerr << "raster_mosaic: skipping '" << p
                      << "' — cell size " << meta.cell_size
                      << " differs from reference " << ref_cs
                      << " (relative diff " << cs_rel_diff << ").\n";
            continue;
        }

        // --- Pixel alignment check (0.01 % of cell_size tolerance) ---
        //
        // The x_min and y_min of this raster must sit on the reference grid,
        // i.e. the differences must be integer multiples of cell_size.
        auto misaligned = [&](double a, double b) -> bool
        {
            double rem = std::fmod(std::abs(a - b), ref_cs);
            if (rem > ref_cs * 0.5) rem = ref_cs - rem;   // nearest integer
            return rem > ref_cs * 1e-4;
        };

        if (misaligned(meta.x_min, ref_xmin) || misaligned(meta.y_min, ref_ymin))
        {
            std::cerr << "raster_mosaic: skipping '" << p
                      << "' — grid origin is not pixel-aligned with reference "
                      << "(x_min=" << meta.x_min
                      << " y_min=" << meta.y_min << ").\n";
            continue;
        }

        accepted.push_back({p, meta});
    }

    if (accepted.empty())
        throw std::runtime_error(
            "raster_mosaic: no compatible rasters found in the input list.");

    std::cout << "raster_mosaic: accepted " << accepted.size()
              << " of " << paths.size() << " rasters.\n";

    // -------------------------------------------------------------------------
    // Compute output extent (union of all accepted bounding boxes)
    // -------------------------------------------------------------------------
    const double cs = ref_cs;

    double out_xmin = accepted[0].meta.x_min;
    double out_ymin = accepted[0].meta.y_min;
    double out_xmax = accepted[0].meta.x_max();
    double out_ymax = accepted[0].meta.y_max();

    for (const auto& e : accepted)
    {
        out_xmin = std::min(out_xmin, e.meta.x_min);
        out_ymin = std::min(out_ymin, e.meta.y_min);
        out_xmax = std::max(out_xmax, e.meta.x_max());
        out_ymax = std::max(out_ymax, e.meta.y_max());
    }

    const int out_cols = static_cast<int>(std::round((out_xmax - out_xmin) / cs));
    const int out_rows = static_cast<int>(std::round((out_ymax - out_ymin) / cs));

    std::cout << "raster_mosaic: output extent"
              << "  x=[" << out_xmin << ", " << out_xmax << "]"
              << "  y=[" << out_ymin << ", " << out_ymax << "]"
              << "  size=" << out_cols << " x " << out_rows << "\n";

    // -------------------------------------------------------------------------
    // Build output metadata
    // -------------------------------------------------------------------------
    lsdtt_raster_meta out_meta = accepted[0].meta;
    out_meta.n_rows = out_rows;
    out_meta.n_cols = out_cols;
    out_meta.x_min  = out_xmin;
    out_meta.y_min  = out_ymin;

    // -------------------------------------------------------------------------
    // Allocate output arrays (initialised to nodata / mask 0)
    // -------------------------------------------------------------------------
    const std::size_t out_npix =
        static_cast<std::size_t>(out_rows) * static_cast<std::size_t>(out_cols);
    const float ndv = static_cast<float>(out_meta.no_data_value);

    std::vector<float>   out_data(out_npix, ndv);
    std::vector<uint8_t> out_mask(out_npix, uint8_t{0});

    // -------------------------------------------------------------------------
    // Pass 2 — paste each raster in list order (last valid wins)
    // -------------------------------------------------------------------------
    for (const auto& e : accepted)
    {
        const lsdtt_raster_meta& sm = e.meta;

        // Pixel offsets: where does the top-left corner of this raster land
        // in the output grid?
        //   col_off: distance in pixels from the output left edge
        //   row_off: distance in pixels from the output top edge
        //            (output row 0 corresponds to y = out_ymax)
        const int col_off = static_cast<int>(
            std::round((sm.x_min  - out_xmin) / cs));
        const int row_off = static_cast<int>(
            std::round((out_ymax  - sm.y_max()) / cs));

        // Read the entire source raster
        lsdtt_data_loader ld;
        ld.open_raster(e.path);
        const lsdtt_raster_tile tile = ld.read_full();

        const int src_rows = sm.n_rows;
        const int src_cols = sm.n_cols;

        std::size_t pasted = 0;

        for (int r = 0; r < src_rows; ++r)
        {
            const int out_r = row_off + r;
            if (out_r < 0 || out_r >= out_rows) continue;

            const std::size_t src_row_base =
                static_cast<std::size_t>(r) * static_cast<std::size_t>(src_cols);
            const std::size_t out_row_base =
                static_cast<std::size_t>(out_r) * static_cast<std::size_t>(out_cols);

            for (int c = 0; c < src_cols; ++c)
            {
                const int out_c = col_off + c;
                if (out_c < 0 || out_c >= out_cols) continue;

                const std::size_t src_idx = src_row_base + static_cast<std::size_t>(c);
                if (tile.mask[src_idx] == 0) continue;   // nodata — skip

                const std::size_t out_idx =
                    out_row_base + static_cast<std::size_t>(out_c);
                out_data[out_idx] = tile.data[src_idx];
                out_mask[out_idx] = uint8_t{1};
                ++pasted;
            }
        }

        std::cout << "  pasted: " << e.path
                  << "  col_off=" << col_off << "  row_off=" << row_off
                  << "  valid_pixels=" << pasted << "\n";
    }

    return lsdtt_raster_cpu(
        std::move(out_data), std::move(out_mask), out_meta);
}

} // namespace lsdtt
