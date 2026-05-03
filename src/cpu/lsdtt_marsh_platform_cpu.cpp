// lsdtt_marsh_platform_cpu.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html
//
// A port of:
// Goodwin, G.C.H., Mudd, S.M., Clubb, F.J., 2018. 
// Unsupervised detection of salt marsh platforms: a topographic method. 
// Earth Surface Dynamics 6, 239–255. https://doi.org/10.5194/esurf-6-239-2018
#include "lsdtt_marsh_platform_cpu.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>
#include <numeric>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "lsdtt_raster_surface_metrics_cpu.hpp"
#include "lsdtt_raster_binary_mask_cpu.hpp"
#include "../utils/lsdtt_histogram.hpp"
#include "../utils/lsdtt_raster_geometry.hpp"

namespace lsdtt
{
namespace
{

// -----------------------------------------------------------------------------
// Basic helpers
// -----------------------------------------------------------------------------
inline bool is_valid_scalar_(
    const lsdtt_raster_meta&       meta,
    const std::vector<float>&      data,
    const std::vector<uint8_t>&    mask,
    std::size_t                    idx)
{
    return mask[idx] != 0 &&
           !meta.is_nodata(data[idx]) &&
           std::isfinite(data[idx]);
}

inline bool is_binary_fg_(
    const lsdtt_raster_meta&       meta,
    const std::vector<float>&      data,
    const std::vector<uint8_t>&    mask,
    std::size_t                    idx)
{
    return is_valid_scalar_(meta, data, mask, idx) && data[idx] > 0.5f;
}

inline std::size_t idx_of_(int row, int col, int ncols)
{
    return static_cast<std::size_t>(row) * static_cast<std::size_t>(ncols) +
           static_cast<std::size_t>(col);
}

inline bool in_bounds_(int row, int col, int nrows, int ncols)
{
    return row >= 0 && row < nrows && col >= 0 && col < ncols;
}

inline double pixel_distance_(int r0, int c0, int r1, int c1)
{
    const double dr = static_cast<double>(r0 - r1);
    const double dc = static_cast<double>(c0 - c1);
    return std::sqrt(dr * dr + dc * dc);
}

static void validate_request_(const marsh_platform_request& request)
{
    if (request.histogram_n_bins < 3)
    {
        throw std::invalid_argument(
            "lsdtt::find_marsh_platform: histogram_n_bins must be >= 3"
        );
    }

    if (request.max_scarp_growth_iterations < 1)
    {
        throw std::invalid_argument(
            "lsdtt::find_marsh_platform: max_scarp_growth_iterations must be >= 1"
        );
    }

    if (request.max_platform_growth_iterations < 1)
    {
        throw std::invalid_argument(
            "lsdtt::find_marsh_platform: max_platform_growth_iterations must be >= 1"
        );
    }

    if (request.cutoff_histogram_run_length_bins < 1)
    {
        throw std::invalid_argument(
            "lsdtt::find_marsh_platform: cutoff_histogram_run_length_bins must be >= 1"
        );
    }

    if (request.scarp_clean_min_neighbour_pixels < 1)
    {
        throw std::invalid_argument(
            "lsdtt::find_marsh_platform: scarp_clean_min_neighbour_pixels must be >= 1"
        );
    }

    if (request.min_platform_component_pixels < 1)
    {
        throw std::invalid_argument(
            "lsdtt::find_marsh_platform: min_platform_component_pixels must be >= 1"
        );
    }
}

static std::vector<uint8_t> make_valid_domain_mask_(const lsdtt_raster_cpu& raster)
{
    const lsdtt_raster_meta& meta = raster.get_meta();
    const std::vector<float>& data = raster.get_data();
    const std::vector<uint8_t>& mask = raster.get_mask();

    std::vector<uint8_t> out(data.size(), uint8_t{0});
    for (std::size_t i = 0; i < data.size(); ++i)
    {
        if (is_valid_scalar_(meta, data, mask, i))
        {
            out[i] = uint8_t{1};
        }
    }
    return out;
}

static lsdtt_raster_cpu raster_from_values_on_valid_domain_(
    std::vector<float>&&           values,
    const std::vector<uint8_t>&    valid_domain,
    const lsdtt_raster_meta&       meta)
{
    if (values.size() != valid_domain.size())
    {
        throw std::runtime_error(
            "lsdtt_marsh_platform_cpu: values and valid_domain size mismatch"
        );
    }

    std::vector<float> out(values.size(), static_cast<float>(meta.no_data_value));
    std::vector<uint8_t> out_mask(values.size(), uint8_t{0});

    for (std::size_t i = 0; i < values.size(); ++i)
    {
        if (valid_domain[i] == 0) continue;
        out[i] = values[i];
        out_mask[i] = uint8_t{1};
    }

    return lsdtt_raster_cpu(std::move(out), std::move(out_mask), meta);
}

static lsdtt_raster_cpu binary_raster_from_mask_(
    const std::vector<uint8_t>&    fg,
    const std::vector<uint8_t>&    valid_domain,
    const lsdtt_raster_meta&       meta)
{
    if (fg.size() != valid_domain.size())
    {
        throw std::runtime_error(
            "lsdtt_marsh_platform_cpu: binary mask and valid_domain size mismatch"
        );
    }

    std::vector<float> out(fg.size(), static_cast<float>(meta.no_data_value));
    std::vector<uint8_t> out_mask(fg.size(), uint8_t{0});

    for (std::size_t i = 0; i < fg.size(); ++i)
    {
        if (valid_domain[i] == 0) continue;
        out[i] = (fg[i] != 0) ? 1.0f : 0.0f;
        out_mask[i] = uint8_t{1};
    }

    return lsdtt_raster_cpu(std::move(out), std::move(out_mask), meta);
}

// -----------------------------------------------------------------------------
// Relative relief
// -----------------------------------------------------------------------------
static lsdtt_raster_cpu compute_relative_relief_(
    const lsdtt_raster_cpu& raster,
    int                     n_threads)
{
    const lsdtt_raster_meta& meta = raster.get_meta();
    const std::vector<float>& data = raster.get_data();
    const std::vector<uint8_t>& mask = raster.get_mask();

    const std::size_t npix = data.size();
    const float ndv = static_cast<float>(meta.no_data_value);

    double vmin = std::numeric_limits<double>::infinity();
    double vmax = -std::numeric_limits<double>::infinity();

    for (std::size_t i = 0; i < npix; ++i)
    {
        if (!is_valid_scalar_(meta, data, mask, i)) continue;
        const double v = static_cast<double>(data[i]);
        vmin = std::min(vmin, v);
        vmax = std::max(vmax, v);
    }

    if (!std::isfinite(vmin) || !std::isfinite(vmax))
    {
        throw std::runtime_error(
            "lsdtt::find_marsh_platform: relative_relief found no valid pixels"
        );
    }

    const double denom = std::max(vmax - vmin, 1.0e-12);

    std::vector<float> out(npix, ndv);
    std::vector<uint8_t> out_mask(npix, uint8_t{0});

#ifdef _OPENMP
#pragma omp parallel for if(npix > 100000) num_threads(n_threads > 0 ? n_threads : omp_get_max_threads())
#endif
    for (std::int64_t i = 0; i < static_cast<std::int64_t>(npix); ++i)
    {
        if (!is_valid_scalar_(meta, data, mask, static_cast<std::size_t>(i))) continue;

        out[static_cast<std::size_t>(i)] =
            static_cast<float>((static_cast<double>(data[static_cast<std::size_t>(i)]) - vmin) / denom);
        out_mask[static_cast<std::size_t>(i)] = uint8_t{1};
    }

    return lsdtt_raster_cpu(std::move(out), std::move(out_mask), meta);
}

// -----------------------------------------------------------------------------
// Internal slope
//
// For marsh-platform extraction the DEM-only path should be fast and simple.
// We therefore use the optimised 3x3 gradient-magnitude slope raster when no
// external slope raster is supplied.
// -----------------------------------------------------------------------------
static lsdtt_raster_cpu compute_internal_slope_(
    const lsdtt_raster_cpu& dem)
{
    return compute_gradient_magnitude_3x3_square(dem);
}

// -----------------------------------------------------------------------------
// Histogram diagnostics
// -----------------------------------------------------------------------------
static marsh_histogram_diagnostics histogram_with_derivative_(
    const std::vector<double>& values,
    int                        n_bins,
    double                     derivative_threshold)
{
    marsh_histogram_diagnostics diag;

    if (values.empty())
    {
        return diag;
    }

    const lsdtt_histogram hist = compute_histogram(values, static_cast<std::size_t>(n_bins));
    if (hist.empty())
    {
        return diag;
    }

    diag.bin_centres   = hist.bin_centres;
    diag.probabilities = hist.probabilities;
    diag.derivative.assign(hist.probabilities.size(), 0.0);

    if (hist.probabilities.size() > 1)
    {
        const double x_step =
            (hist.bin_centres.back() - hist.bin_centres.front()) /
            static_cast<double>(hist.bin_centres.size());

        if (std::abs(x_step) > 0.0)
        {
            for (std::size_t j = 1; j < hist.probabilities.size(); ++j)
            {
                diag.derivative[j] =
                    (hist.probabilities[j] - hist.probabilities[j - 1]) / x_step;
            }
        }
    }

    double threshold = 0.0;
    for (std::size_t j = 1; j + 1 < diag.derivative.size(); ++j)
    {
        if (diag.derivative[j] < derivative_threshold &&
            diag.derivative[j + 1] >= derivative_threshold)
        {
            threshold = diag.bin_centres[j];
        }
    }

    diag.threshold_value = threshold;
    return diag;
}

struct cutoff_result_
{
    marsh_histogram_diagnostics diag;
    double dominant_peak_value = 0.0;
};

static cutoff_result_ compute_cutoff_histogram_(
    const std::vector<float>&      platform_order,
    const lsdtt_raster_cpu&        dem,
    const lsdtt_raster_cpu&        search_space,
    int                            n_bins,
    int                            run_length_bins)
{
    cutoff_result_ out;

    const lsdtt_raster_meta& meta = dem.get_meta();
    const std::vector<float>& dd = dem.get_data();
    const std::vector<uint8_t>& dm = dem.get_mask();

    const std::vector<float>& sd = search_space.get_data();
    const std::vector<uint8_t>& sm = search_space.get_mask();

    std::vector<double> platform_values;
    std::vector<double> search_values;
    platform_values.reserve(dd.size());
    search_values.reserve(dd.size());

    for (std::size_t i = 0; i < dd.size(); ++i)
    {
        if (!is_valid_scalar_(meta, dd, dm, i)) continue;

        platform_values.push_back((platform_order[i] > 0.5f) ? static_cast<double>(dd[i]) : 0.0);
        search_values.push_back((sm[i] != 0 && sd[i] > 0.5f) ? static_cast<double>(dd[i]) : 0.0);
    }

    if (platform_values.empty())
    {
        return out;
    }

    const lsdtt_histogram ph = compute_histogram(platform_values, static_cast<std::size_t>(n_bins));
    if (ph.empty())
    {
        return out;
    }

    out.diag.bin_centres   = ph.bin_centres;
    out.diag.probabilities = ph.probabilities;

    std::size_t dominant_index = ph.probabilities.size() - 1;
    double max_nonzero_peak = 0.0;

    if (ph.probabilities.size() > 1)
    {
        max_nonzero_peak = *std::max_element(ph.probabilities.begin() + 1, ph.probabilities.end());

        for (std::size_t j = 1; j + 1 < ph.probabilities.size(); ++j)
        {
            if (ph.probabilities[j] > 0.9 * max_nonzero_peak &&
                ph.probabilities[j] > ph.probabilities[j - 1] &&
                ph.probabilities[j] > ph.probabilities[j + 1])
            {
                dominant_index = j;
            }
        }
    }

    out.dominant_peak_value = ph.bin_centres[dominant_index];

    const double mean_prob =
        std::accumulate(ph.probabilities.begin(), ph.probabilities.end(), 0.0) /
        static_cast<double>(ph.probabilities.size());

    int low_run = 0;
    double cutoff_z = 0.0;

    for (std::int64_t j = static_cast<std::int64_t>(dominant_index); j >= 0; --j)
    {
        if (ph.probabilities[static_cast<std::size_t>(j)] < mean_prob)
        {
            ++low_run;
        }
        else
        {
            low_run = 0;
        }

        if (low_run > run_length_bins)
        {
            cutoff_z = ph.bin_centres[static_cast<std::size_t>(j)];
            break;
        }

        if (low_run > run_length_bins / 2)
        {
            cutoff_z = ph.bin_centres[static_cast<std::size_t>(j)];
        }
    }

    if (!search_values.empty())
    {
        const lsdtt_histogram sh = compute_histogram(search_values, static_cast<std::size_t>(n_bins));

        if (!sh.empty() && sh.probabilities.size() > 1)
        {
            const double max_prob =
                *std::max_element(sh.probabilities.begin() + 1, sh.probabilities.end());
            const double low_prop = 0.01 * max_prob;

            for (std::size_t j = 1; j < sh.probabilities.size(); ++j)
            {
                if (sh.probabilities[j] > low_prop)
                {
                    cutoff_z = std::max(cutoff_z, sh.bin_centres[j]);
                    break;
                }
            }
        }
    }

    out.diag.threshold_value = cutoff_z;
    return out;
}

// -----------------------------------------------------------------------------
// Search space
// -----------------------------------------------------------------------------
struct search_space_result_
{
    lsdtt_raster_cpu crossover;
    lsdtt_raster_cpu search_space;
    marsh_histogram_diagnostics histogram;
};

static search_space_result_ compute_search_space_(
    const lsdtt_raster_cpu& dem_rr,
    const lsdtt_raster_cpu& slope_rr,
    const marsh_platform_request& request)
{
    const lsdtt_raster_meta& meta = dem_rr.get_meta();

    const std::vector<float>& dd = dem_rr.get_data();
    const std::vector<uint8_t>& dm = dem_rr.get_mask();
    const std::vector<float>& sd = slope_rr.get_data();
    const std::vector<uint8_t>& sm = slope_rr.get_mask();

    if (!same_grid(meta, slope_rr.get_meta()))
    {
        throw std::invalid_argument(
            "lsdtt::find_marsh_platform: DEM and slope relative-relief rasters are not co-registered"
        );
    }

    const int nrows = meta.n_rows;
    const int ncols = meta.n_cols;
    const std::size_t npix = dd.size();
    const float ndv = static_cast<float>(meta.no_data_value);

    std::vector<float> crossover(npix, ndv);
    std::vector<uint8_t> crossover_mask(npix, uint8_t{0});

    std::vector<double> crossover_values;
    crossover_values.reserve(npix / 2);

#ifdef _OPENMP
#pragma omp parallel if(npix > 100000) num_threads(request.n_threads > 0 ? request.n_threads : omp_get_max_threads())
    {
        std::vector<double> local_values;
        local_values.reserve(1024);

#pragma omp for nowait
        for (std::int64_t i = 0; i < static_cast<std::int64_t>(npix); ++i)
        {
            const std::size_t idx = static_cast<std::size_t>(i);

            if (!is_valid_scalar_(meta, dd, dm, idx)) continue;
            if (!is_valid_scalar_(meta, sd, sm, idx)) continue;

            const float v = dd[idx] * sd[idx];
            crossover[idx] = v;
            crossover_mask[idx] = uint8_t{1};
            local_values.push_back(static_cast<double>(v));
        }

#pragma omp critical
        crossover_values.insert(
            crossover_values.end(),
            local_values.begin(),
            local_values.end());
    }
#else
    for (std::size_t i = 0; i < npix; ++i)
    {
        if (!is_valid_scalar_(meta, dd, dm, i)) continue;
        if (!is_valid_scalar_(meta, sd, sm, i)) continue;

        const float v = dd[i] * sd[i];
        crossover[i] = v;
        crossover_mask[i] = uint8_t{1};
        crossover_values.push_back(static_cast<double>(v));
    }
#endif

    search_space_result_ out;
    out.crossover = lsdtt_raster_cpu(
        std::move(crossover), std::move(crossover_mask), meta);

    out.histogram = histogram_with_derivative_(
        crossover_values,
        request.histogram_n_bins,
        request.search_space_histogram_derivative_threshold);

    const double threshold = out.histogram.threshold_value;

    std::vector<float> search(npix, ndv);
    std::vector<uint8_t> search_mask(npix, uint8_t{0});

#ifdef _OPENMP
#pragma omp parallel for if(npix > 100000) num_threads(request.n_threads > 0 ? request.n_threads : omp_get_max_threads())
#endif
    for (std::int64_t i = 0; i < static_cast<std::int64_t>(npix); ++i)
    {
        const std::size_t idx = static_cast<std::size_t>(i);
        if (!is_valid_scalar_(meta, out.crossover.get_data(), out.crossover.get_mask(), idx)) continue;

        search_mask[idx] = uint8_t{1};
        search[idx] = (out.crossover.get_data()[idx] > static_cast<float>(threshold)) ? 1.0f : 0.0f;
    }

    for (int col = 0; col < ncols; ++col)
    {
        const std::size_t i0 = idx_of_(0, col, ncols);
        const std::size_t i1 = idx_of_(nrows - 1, col, ncols);
        if (search_mask[i0] != 0) search[i0] = 0.0f;
        if (search_mask[i1] != 0) search[i1] = 0.0f;
    }

    for (int row = 0; row < nrows; ++row)
    {
        const std::size_t i0 = idx_of_(row, 0, ncols);
        const std::size_t i1 = idx_of_(row, ncols - 1, ncols);
        if (search_mask[i0] != 0) search[i0] = 0.0f;
        if (search_mask[i1] != 0) search[i1] = 0.0f;
    }

    out.search_space = lsdtt_raster_cpu(
        std::move(search), std::move(search_mask), meta);

    return out;
}

// -----------------------------------------------------------------------------
// Scarp seed detection
// -----------------------------------------------------------------------------
static std::vector<uint8_t> detect_scarp_seeds_(
    const lsdtt_raster_cpu& slope,
    const lsdtt_raster_cpu& search_space,
    int                     n_threads)
{
    const lsdtt_raster_meta& meta = slope.get_meta();
    const std::vector<float>& sd = slope.get_data();
    const std::vector<uint8_t>& sm = slope.get_mask();
    const std::vector<float>& qd = search_space.get_data();
    const std::vector<uint8_t>& qm = search_space.get_mask();

    const int nrows = meta.n_rows;
    const int ncols = meta.n_cols;

    std::vector<uint8_t> seeds(sd.size(), uint8_t{0});

#ifdef _OPENMP
#pragma omp parallel for collapse(2) if(sd.size() > 100000) num_threads(n_threads > 0 ? n_threads : omp_get_max_threads())
#endif
    for (int row = 1; row < nrows - 1; ++row)
    {
        for (int col = 1; col < ncols - 1; ++col)
        {
            const std::size_t idx = idx_of_(row, col, ncols);

            if (!is_valid_scalar_(meta, sd, sm, idx)) continue;
            if (!(qm[idx] != 0 && qd[idx] > 0.5f)) continue;

            const float centre = sd[idx];
            bool is_max = true;

            for (int dr = -1; dr <= 1 && is_max; ++dr)
            {
                for (int dc = -1; dc <= 1; ++dc)
                {
                    if (dr == 0 && dc == 0) continue;
                    const std::size_t nidx = idx_of_(row + dr, col + dc, ncols);
                    if (!is_valid_scalar_(meta, sd, sm, nidx)) continue;
                    if (sd[nidx] > centre)
                    {
                        is_max = false;
                        break;
                    }
                }
            }

            if (is_max)
            {
                seeds[idx] = uint8_t{1};
            }
        }
    }

    return seeds;
}

static bool choose_best_neighbour_(
    int                         row,
    int                         col,
    const lsdtt_raster_cpu&     slope,
    const lsdtt_raster_cpu&     search_space,
    const std::vector<uint8_t>& available,
    int&                        best_row,
    int&                        best_col,
    int                         avoid_row = std::numeric_limits<int>::min(),
    int                         avoid_col = std::numeric_limits<int>::min(),
    int                         grandad_row = std::numeric_limits<int>::min(),
    int                         grandad_col = std::numeric_limits<int>::min())
{
    const lsdtt_raster_meta& meta = slope.get_meta();
    const std::vector<float>& sd = slope.get_data();
    const std::vector<uint8_t>& sm = slope.get_mask();
    const std::vector<float>& qd = search_space.get_data();
    const std::vector<uint8_t>& qm = search_space.get_mask();

    const int nrows = meta.n_rows;
    const int ncols = meta.n_cols;

    float best_value = -std::numeric_limits<float>::infinity();
    bool found = false;

    for (int dr = -1; dr <= 1; ++dr)
    {
        for (int dc = -1; dc <= 1; ++dc)
        {
            if (dr == 0 && dc == 0) continue;

            const int nr = row + dr;
            const int nc = col + dc;
            if (!in_bounds_(nr, nc, nrows, ncols)) continue;

            const std::size_t nidx = idx_of_(nr, nc, ncols);

            if (available[nidx] == 0) continue;
            if (!is_valid_scalar_(meta, sd, sm, nidx)) continue;
            if (!(qm[nidx] != 0 && qd[nidx] > 0.5f)) continue;

            if (avoid_row != std::numeric_limits<int>::min())
            {
                if (pixel_distance_(nr, nc, avoid_row, avoid_col) <= std::sqrt(2.0))
                {
                    continue;
                }
            }

            if (grandad_row != std::numeric_limits<int>::min())
            {
                if (pixel_distance_(nr, nc, grandad_row, grandad_col) <= std::sqrt(2.0))
                {
                    continue;
                }
            }

            if (sd[nidx] > best_value)
            {
                best_value = sd[nidx];
                best_row = nr;
                best_col = nc;
                found = true;
            }
        }
    }

    return found;
}

// -----------------------------------------------------------------------------
// Scarp growth
// -----------------------------------------------------------------------------
static void grow_scarps_(
    const lsdtt_raster_cpu&     slope,
    const lsdtt_raster_cpu&     search_space,
    const std::vector<uint8_t>& seed_mask,
    int                         max_iterations,
    std::vector<float>&         scarp_order)
{
    const lsdtt_raster_meta& meta = slope.get_meta();
    const int nrows = meta.n_rows;
    const int ncols = meta.n_cols;
    const std::size_t npix = scarp_order.size();

    std::vector<uint8_t> available(npix, uint8_t{0});
    for (std::size_t i = 0; i < npix; ++i)
    {
        if (is_binary_fg_(search_space.get_meta(), search_space.get_data(), search_space.get_mask(), i))
        {
            available[i] = uint8_t{1};
        }
    }

    // Seeds are order 1 and become unavailable.
    for (std::size_t i = 0; i < npix; ++i)
    {
        if (seed_mask[i] != 0)
        {
            scarp_order[i] = 1.0f;
            available[i] = uint8_t{0};
        }
    }

    // Initiate ridge (order 2).
    for (int row = 0; row < nrows; ++row)
    {
        for (int col = 0; col < ncols; ++col)
        {
            const std::size_t idx = idx_of_(row, col, ncols);
            if (seed_mask[idx] == 0) continue;

            int ridge_count = 0;
            for (int dr = -1; dr <= 1; ++dr)
            {
                for (int dc = -1; dc <= 1; ++dc)
                {
                    const int nr = row + dr;
                    const int nc = col + dc;
                    if (!in_bounds_(nr, nc, nrows, ncols)) continue;
                    if (scarp_order[idx_of_(nr, nc, ncols)] > 0.5f) ++ridge_count;
                }
            }

            int r1 = -1, c1 = -1;
            if (!choose_best_neighbour_(row, col, slope, search_space, available, r1, c1))
            {
                continue;
            }

            const std::size_t i1 = idx_of_(r1, c1, ncols);
            scarp_order[i1] = 2.0f;
            available[i1] = uint8_t{0};

            if (ridge_count == 1)
            {
                int r2 = -1, c2 = -1;
                if (choose_best_neighbour_(row, col, slope, search_space, available, r2, c2, r1, c1))
                {
                    const std::size_t i2 = idx_of_(r2, c2, ncols);
                    scarp_order[i2] = 2.0f;
                    available[i2] = uint8_t{0};
                }
            }
        }
    }

    // Continue ridge.
    for (int order = 3; order <= max_iterations; ++order)
    {
        bool any_added = false;

        for (int row = 0; row < nrows; ++row)
        {
            for (int col = 0; col < ncols; ++col)
            {
                const std::size_t idx = idx_of_(row, col, ncols);
                if (std::round(scarp_order[idx]) != order - 1) continue;

                int ridge_count = 0;
                int grandad_row = std::numeric_limits<int>::min();
                int grandad_col = std::numeric_limits<int>::min();

                for (int dr = -1; dr <= 1; ++dr)
                {
                    for (int dc = -1; dc <= 1; ++dc)
                    {
                        const int nr = row + dr;
                        const int nc = col + dc;
                        if (!in_bounds_(nr, nc, nrows, ncols)) continue;

                        const std::size_t nidx = idx_of_(nr, nc, ncols);
                        if (scarp_order[nidx] > 0.5f)
                        {
                            ++ridge_count;
                            if (std::round(scarp_order[nidx]) == order - 2)
                            {
                                grandad_row = nr;
                                grandad_col = nc;
                            }
                        }
                    }
                }

                if (ridge_count != 2) continue;
                if (grandad_row == std::numeric_limits<int>::min()) continue;

                int nr = -1, nc = -1;
                if (!choose_best_neighbour_(
                        row, col, slope, search_space, available,
                        nr, nc,
                        std::numeric_limits<int>::min(),
                        std::numeric_limits<int>::min(),
                        grandad_row, grandad_col))
                {
                    continue;
                }

                const std::size_t nidx = idx_of_(nr, nc, ncols);
                scarp_order[nidx] = static_cast<float>(order);
                available[nidx] = uint8_t{0};
                any_added = true;
            }
        }

        if (!any_added) break;
    }
}

// -----------------------------------------------------------------------------
// Scarp cleaning
// -----------------------------------------------------------------------------
static std::vector<uint8_t> clean_scarps_(
    const std::vector<float>&     scarp_order,
    const lsdtt_raster_cpu&       dem,
    const marsh_platform_request& request)
{
    const lsdtt_raster_meta& meta = dem.get_meta();
    const std::vector<float>& dd = dem.get_data();
    const std::vector<uint8_t>& dm = dem.get_mask();

    const int nrows = meta.n_rows;
    const int ncols = meta.n_cols;
    const std::size_t npix = dd.size();

    std::vector<float> clipped_dem(npix, 0.0f);
    std::vector<double> valid_dem;
    valid_dem.reserve(npix);

    for (std::size_t i = 0; i < npix; ++i)
    {
        if (!is_valid_scalar_(meta, dd, dm, i)) continue;
        clipped_dem[i] = dd[i];
        valid_dem.push_back(static_cast<double>(dd[i]));
    }

    if (valid_dem.empty())
    {
        throw std::runtime_error(
            "lsdtt::find_marsh_platform: no valid DEM pixels during scarp cleaning"
        );
    }

    const double q75 = compute_percentile(valid_dem, 75.0);
    double threshold = 0.0;
    for (double v : valid_dem)
    {
        if (v < q75) threshold = std::max(threshold, v);
    }
    if (!(threshold > 0.0))
    {
        threshold = q75;
    }

    for (std::size_t i = 0; i < npix; ++i)
    {
        if (clipped_dem[i] > static_cast<float>(threshold))
        {
            clipped_dem[i] = static_cast<float>(threshold);
        }
    }

    std::vector<uint8_t> keep1(npix, uint8_t{0});

#ifdef _OPENMP
#pragma omp parallel for collapse(2) if(npix > 100000) num_threads(request.n_threads > 0 ? request.n_threads : omp_get_max_threads())
#endif
    for (int row = 0; row < nrows; ++row)
    {
        for (int col = 0; col < ncols; ++col)
        {
            const std::size_t idx = idx_of_(row, col, ncols);
            if (!(scarp_order[idx] > 0.5f)) continue;

            float local_max = 0.0f;
            for (int dr = -4; dr <= 4; ++dr)
            {
                for (int dc = -4; dc <= 4; ++dc)
                {
                    const int nr = row + dr;
                    const int nc = col + dc;
                    if (!in_bounds_(nr, nc, nrows, ncols)) continue;

                    const std::size_t nidx = idx_of_(nr, nc, ncols);
                    if (!is_valid_scalar_(meta, dd, dm, nidx)) continue;
                    local_max = std::max(local_max, clipped_dem[nidx]);
                }
            }

            if (threshold > 0.0 &&
                (static_cast<double>(local_max) / threshold) >= request.scarp_clean_relief_fraction)
            {
                keep1[idx] = uint8_t{1};
            }
        }
    }

    std::vector<uint8_t> keep2(npix, uint8_t{0});

#ifdef _OPENMP
#pragma omp parallel for collapse(2) if(npix > 100000) num_threads(request.n_threads > 0 ? request.n_threads : omp_get_max_threads())
#endif
    for (int row = 0; row < nrows; ++row)
    {
        for (int col = 0; col < ncols; ++col)
        {
            const std::size_t idx = idx_of_(row, col, ncols);
            if (keep1[idx] == 0) continue;

            int count = 0;
            for (int dr = -4; dr <= 4; ++dr)
            {
                for (int dc = -4; dc <= 4; ++dc)
                {
                    const int nr = row + dr;
                    const int nc = col + dc;
                    if (!in_bounds_(nr, nc, nrows, ncols)) continue;

                    if (keep1[idx_of_(nr, nc, ncols)] != 0) ++count;
                }
            }

            if (count >= request.scarp_clean_min_neighbour_pixels)
            {
                keep2[idx] = uint8_t{1};
            }
        }
    }

    return keep2;
}

// -----------------------------------------------------------------------------
// Precompute local 11x11 DEM maxima for platform growth
// -----------------------------------------------------------------------------
static std::vector<float> precompute_local_dem_max_11x11_(
    const lsdtt_raster_cpu& dem,
    int                     n_threads)
{
    const lsdtt_raster_meta& meta = dem.get_meta();
    const std::vector<float>& dd = dem.get_data();
    const std::vector<uint8_t>& dm = dem.get_mask();

    const int nrows = meta.n_rows;
    const int ncols = meta.n_cols;
    const std::size_t npix = dd.size();

    std::vector<float> out(npix, static_cast<float>(meta.no_data_value));

#ifdef _OPENMP
#pragma omp parallel for collapse(2) if(npix > 100000) num_threads(n_threads > 0 ? n_threads : omp_get_max_threads())
#endif
    for (int row = 0; row < nrows; ++row)
    {
        for (int col = 0; col < ncols; ++col)
        {
            const std::size_t idx = idx_of_(row, col, ncols);
            if (!is_valid_scalar_(meta, dd, dm, idx)) continue;

            float local_max = -std::numeric_limits<float>::infinity();

            for (int dr = -5; dr <= 5; ++dr)
            {
                for (int dc = -5; dc <= 5; ++dc)
                {
                    const int nr = row + dr;
                    const int nc = col + dc;
                    if (!in_bounds_(nr, nc, nrows, ncols)) continue;

                    const std::size_t nidx = idx_of_(nr, nc, ncols);
                    if (!is_valid_scalar_(meta, dd, dm, nidx)) continue;

                    local_max = std::max(local_max, dd[nidx]);
                }
            }

            out[idx] = local_max;
        }
    }

    return out;
}

// -----------------------------------------------------------------------------
// Platform seeding and growth
// -----------------------------------------------------------------------------
static std::vector<float> seed_platform_from_scarps_(
    const std::vector<uint8_t>& scarp_mask,
    const lsdtt_raster_cpu&     dem,
    int                         n_threads)
{
    const lsdtt_raster_meta& meta = dem.get_meta();
    const std::vector<float>& dd = dem.get_data();
    const std::vector<uint8_t>& dm = dem.get_mask();

    const int nrows = meta.n_rows;
    const int ncols = meta.n_cols;
    const std::size_t npix = dd.size();

    std::vector<uint8_t> seed_bin(npix, uint8_t{0});

#ifdef _OPENMP
#pragma omp parallel for collapse(2) if(npix > 100000) num_threads(n_threads > 0 ? n_threads : omp_get_max_threads())
#endif
    for (int row = 1; row < nrows - 1; ++row)
    {
        for (int col = 1; col < ncols - 1; ++col)
        {
            const std::size_t idx = idx_of_(row, col, ncols);
            if (!is_valid_scalar_(meta, dd, dm, idx)) continue;
            if (scarp_mask[idx] != 0) continue;

            bool seed = false;
            for (int dr = -1; dr <= 1 && !seed; ++dr)
            {
                for (int dc = -1; dc <= 1; ++dc)
                {
                    if (dr == 0 && dc == 0) continue;

                    const int nr = row + dr;
                    const int nc = col + dc;
                    const std::size_t nidx = idx_of_(nr, nc, ncols);

                    if (scarp_mask[nidx] == 0) continue;
                    if (!is_valid_scalar_(meta, dd, dm, nidx)) continue;

                    if (dd[idx] >= dd[nidx])
                    {
                        seed = true;
                        break;
                    }
                }
            }

            if (seed)
            {
                seed_bin[idx] = uint8_t{1};
            }
        }
    }

    // Cull sparse seed pixels.
    std::vector<float> order(npix, 0.0f);

#ifdef _OPENMP
#pragma omp parallel for collapse(2) if(npix > 100000) num_threads(n_threads > 0 ? n_threads : omp_get_max_threads())
#endif
    for (int row = 1; row < nrows - 1; ++row)
    {
        for (int col = 1; col < ncols - 1; ++col)
        {
            const std::size_t idx = idx_of_(row, col, ncols);
            if (seed_bin[idx] == 0) continue;

            int count = 0;
            for (int dr = -1; dr <= 1; ++dr)
            {
                for (int dc = -1; dc <= 1; ++dc)
                {
                    const std::size_t nidx = idx_of_(row + dr, col + dc, ncols);
                    if (seed_bin[nidx] != 0) ++count;
                }
            }

            if (count > 2)
            {
                order[idx] = 1.0f;
            }
        }
    }

    return order;
}

static void grow_platform_(
    std::vector<float>&         platform_order,
    const std::vector<uint8_t>& scarp_mask,
    const lsdtt_raster_cpu&     dem,
    const std::vector<float>&   local_dem_max_11x11,
    int                         max_iterations)
{
    const lsdtt_raster_meta& meta = dem.get_meta();
    const std::vector<float>& dd = dem.get_data();
    const std::vector<uint8_t>& dm = dem.get_mask();

    const int nrows = meta.n_rows;
    const int ncols = meta.n_cols;

    std::vector<std::size_t> frontier;
    frontier.reserve(1024);
    for (std::size_t i = 0; i < platform_order.size(); ++i)
    {
        if (std::round(platform_order[i]) == 1)
        {
            frontier.push_back(i);
        }
    }

    for (int order = 2; order <= max_iterations; ++order)
    {
        std::vector<std::size_t> next_frontier;
        next_frontier.reserve(frontier.size() * 2 + 16);

        for (std::size_t parent_idx : frontier)
        {
            const int row = static_cast<int>(parent_idx / static_cast<std::size_t>(ncols));
            const int col = static_cast<int>(parent_idx % static_cast<std::size_t>(ncols));

            std::vector<std::pair<int,int>> ridge_pos;
            std::vector<std::pair<int,int>> parent_pos;

            for (int dr = -1; dr <= 1; ++dr)
            {
                for (int dc = -1; dc <= 1; ++dc)
                {
                    const int nr = row + dr;
                    const int nc = col + dc;
                    if (!in_bounds_(nr, nc, nrows, ncols)) continue;

                    const std::size_t nidx = idx_of_(nr, nc, ncols);
                    if (scarp_mask[nidx] != 0)
                    {
                        ridge_pos.push_back({nr, nc});
                    }
                    if (std::round(platform_order[nidx]) == order - 1)
                    {
                        parent_pos.push_back({nr, nc});
                    }
                }
            }

            const float local_max = local_dem_max_11x11[parent_idx];

            for (int dr = -1; dr <= 1; ++dr)
            {
                for (int dc = -1; dc <= 1; ++dc)
                {
                    if (dr == 0 && dc == 0) continue;

                    const int nr = row + dr;
                    const int nc = col + dc;
                    if (!in_bounds_(nr, nc, nrows, ncols)) continue;

                    const std::size_t nidx = idx_of_(nr, nc, ncols);

                    if (!is_valid_scalar_(meta, dd, dm, nidx)) continue;
                    if (scarp_mask[nidx] != 0) continue;
                    if (platform_order[nidx] > 0.5f) continue;

                    if (!(dd[nidx] > local_max - 0.20f)) continue;

                    bool accept = false;

                    if (ridge_pos.empty())
                    {
                        accept = true;
                    }
                    else
                    {
                        double min_ridge_dist  = std::numeric_limits<double>::infinity();
                        double min_parent_dist = std::numeric_limits<double>::infinity();

                        for (const auto& p : ridge_pos)
                        {
                            min_ridge_dist = std::min(
                                min_ridge_dist,
                                pixel_distance_(nr, nc, p.first, p.second));
                        }

                        for (const auto& p : parent_pos)
                        {
                            min_parent_dist = std::min(
                                min_parent_dist,
                                pixel_distance_(nr, nc, p.first, p.second));
                        }

                        if (min_ridge_dist > min_parent_dist)
                        {
                            accept = true;
                        }
                    }

                    if (accept)
                    {
                        platform_order[nidx] = static_cast<float>(order);
                        next_frontier.push_back(nidx);
                    }
                }
            }
        }

        if (next_frontier.empty()) break;
        frontier = std::move(next_frontier);
    }
}

// -----------------------------------------------------------------------------
// Legacy-style fills
// -----------------------------------------------------------------------------
static void apply_cutoff_(
    std::vector<float>&         platform_order,
    const lsdtt_raster_cpu&     dem,
    double                      cutoff_elevation,
    int                         n_threads)
{
    const lsdtt_raster_meta& meta = dem.get_meta();
    const std::vector<float>& dd = dem.get_data();
    const std::vector<uint8_t>& dm = dem.get_mask();

#ifdef _OPENMP
#pragma omp parallel for if(platform_order.size() > 100000) num_threads(n_threads > 0 ? n_threads : omp_get_max_threads())
#endif
    for (std::int64_t i = 0; i < static_cast<std::int64_t>(platform_order.size()); ++i)
    {
        const std::size_t idx = static_cast<std::size_t>(i);
        if (platform_order[idx] <= 0.5f) continue;
        if (!is_valid_scalar_(meta, dd, dm, idx)) continue;

        if (static_cast<double>(dd[idx]) < cutoff_elevation)
        {
            platform_order[idx] = 0.0f;
        }
    }
}

static void fill_high_gaps_(
    std::vector<float>&         platform_order,
    const std::vector<uint8_t>& scarp_mask,
    const lsdtt_raster_cpu&     dem,
    double                      dominant_peak_elevation,
    int                         n_threads)
{
    const lsdtt_raster_meta& meta = dem.get_meta();
    const std::vector<float>& dd = dem.get_data();
    const std::vector<uint8_t>& dm = dem.get_mask();

#ifdef _OPENMP
#pragma omp parallel for if(platform_order.size() > 100000) num_threads(n_threads > 0 ? n_threads : omp_get_max_threads())
#endif
    for (std::int64_t i = 0; i < static_cast<std::int64_t>(platform_order.size()); ++i)
    {
        const std::size_t idx = static_cast<std::size_t>(i);
        if (!is_valid_scalar_(meta, dd, dm, idx)) continue;
        if (scarp_mask[idx] != 0) continue;
        if (platform_order[idx] > 0.5f) continue;

        if (static_cast<double>(dd[idx]) >= dominant_peak_elevation)
        {
            platform_order[idx] = 3.0f;
        }
    }
}

static void fill_pools_(
    std::vector<float>&         platform_order,
    const std::vector<uint8_t>& scarp_mask,
    const lsdtt_raster_cpu&     dem)
{
    const lsdtt_raster_meta& meta = dem.get_meta();
    const int nrows = meta.n_rows;
    const int ncols = meta.n_cols;

    int max_order = 0;
    for (float v : platform_order)
    {
        max_order = std::max(max_order, static_cast<int>(std::round(v)));
    }
    max_order = std::max(max_order, 100);

    for (int iter = 0; iter < 10; ++iter)
    {
        for (int counter = max_order - 1; counter > 2; --counter)
        {
            for (int row = 1; row < nrows - 1; ++row)
            {
                for (int col = 1; col < ncols - 1; ++col)
                {
                    const std::size_t idx = idx_of_(row, col, ncols);
                    if (std::round(platform_order[idx]) != counter + 1) continue;

                    int marsh_count = 0;
                    for (int dr = -1; dr <= 1; ++dr)
                    {
                        for (int dc = -1; dc <= 1; ++dc)
                        {
                            if (platform_order[idx_of_(row + dr, col + dc, ncols)] > 0.5f)
                            {
                                ++marsh_count;
                            }
                        }
                    }

                    if (marsh_count > 6)
                    {
                        for (int dr = -1; dr <= 1; ++dr)
                        {
                            for (int dc = -1; dc <= 1; ++dc)
                            {
                                const std::size_t nidx = idx_of_(row + dr, col + dc, ncols);
                                if (platform_order[nidx] <= 0.5f && scarp_mask[nidx] == 0)
                                {
                                    platform_order[nidx] = static_cast<float>(counter);
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

static void fill_small_holes_(
    std::vector<float>&         platform_order,
    const std::vector<uint8_t>& scarp_mask,
    const lsdtt_raster_cpu&     dem,
    int                         n_threads)
{
    const lsdtt_raster_meta& meta = dem.get_meta();
    const std::vector<float>& dd = dem.get_data();
    const std::vector<uint8_t>& dm = dem.get_mask();

    const int nrows = meta.n_rows;
    const int ncols = meta.n_cols;
    std::vector<uint8_t> fill(platform_order.size(), uint8_t{0});

#ifdef _OPENMP
#pragma omp parallel for collapse(2) if(platform_order.size() > 100000) num_threads(n_threads > 0 ? n_threads : omp_get_max_threads())
#endif
    for (int row = 1; row < nrows - 1; ++row)
    {
        for (int col = 1; col < ncols - 1; ++col)
        {
            const std::size_t idx = idx_of_(row, col, ncols);
            if (!is_valid_scalar_(meta, dd, dm, idx)) continue;
            if (scarp_mask[idx] != 0) continue;
            if (platform_order[idx] > 0.5f) continue;

            int count = 0;
            for (int dr = -1; dr <= 1; ++dr)
            {
                for (int dc = -1; dc <= 1; ++dc)
                {
                    if (dr == 0 && dc == 0) continue;
                    if (platform_order[idx_of_(row + dr, col + dc, ncols)] > 0.5f)
                    {
                        ++count;
                    }
                }
            }

            if (count == 8)
            {
                fill[idx] = uint8_t{1};
            }
        }
    }

    for (std::size_t i = 0; i < platform_order.size(); ++i)
    {
        if (fill[i] != 0)
        {
            platform_order[i] = 105.0f;
        }
    }
}

// -----------------------------------------------------------------------------
// Clean scarps against platform adjacency before reinserting
// -----------------------------------------------------------------------------
static std::vector<uint8_t> clean_scarps_from_platform_(
    const std::vector<uint8_t>& scarp_mask,
    const std::vector<float>&   platform_order,
    const lsdtt_raster_cpu&     dem,
    int                         min_neighbour_pixels,
    int                         n_threads)
{
    const lsdtt_raster_meta& meta = dem.get_meta();
    const int nrows = meta.n_rows;
    const int ncols = meta.n_cols;

    std::vector<uint8_t> keep1(scarp_mask.size(), uint8_t{0});

#ifdef _OPENMP
#pragma omp parallel for collapse(2) if(scarp_mask.size() > 100000) num_threads(n_threads > 0 ? n_threads : omp_get_max_threads())
#endif
    for (int row = 0; row < nrows; ++row)
    {
        for (int col = 0; col < ncols; ++col)
        {
            const std::size_t idx = idx_of_(row, col, ncols);
            if (scarp_mask[idx] == 0) continue;

            bool any_platform = false;
            for (int dr = -1; dr <= 1 && !any_platform; ++dr)
            {
                for (int dc = -1; dc <= 1; ++dc)
                {
                    const int nr = row + dr;
                    const int nc = col + dc;
                    if (!in_bounds_(nr, nc, nrows, ncols)) continue;

                    if (platform_order[idx_of_(nr, nc, ncols)] > 0.5f)
                    {
                        any_platform = true;
                        break;
                    }
                }
            }

            if (any_platform)
            {
                keep1[idx] = uint8_t{1};
            }
        }
    }

    std::vector<uint8_t> keep2(scarp_mask.size(), uint8_t{0});

#ifdef _OPENMP
#pragma omp parallel for collapse(2) if(scarp_mask.size() > 100000) num_threads(n_threads > 0 ? n_threads : omp_get_max_threads())
#endif
    for (int row = 0; row < nrows; ++row)
    {
        for (int col = 0; col < ncols; ++col)
        {
            const std::size_t idx = idx_of_(row, col, ncols);
            if (keep1[idx] == 0) continue;

            int count = 0;
            for (int dr = -4; dr <= 4; ++dr)
            {
                for (int dc = -4; dc <= 4; ++dc)
                {
                    const int nr = row + dr;
                    const int nc = col + dc;
                    if (!in_bounds_(nr, nc, nrows, ncols)) continue;

                    if (keep1[idx_of_(nr, nc, ncols)] != 0) ++count;
                }
            }

            if (count >= min_neighbour_pixels)
            {
                keep2[idx] = uint8_t{1};
            }
        }
    }

    return keep2;
}

} // unnamed namespace

// -----------------------------------------------------------------------------
// find_marsh_platform
// -----------------------------------------------------------------------------
marsh_platform_result find_marsh_platform(
    const lsdtt_raster_cpu&       dem,
    const lsdtt_raster_cpu*       supplied_slope,
    const marsh_platform_request& request)
{
    validate_request_(request);

    const lsdtt_raster_meta& meta = dem.get_meta();
    if (!meta.is_valid())
    {
        throw std::runtime_error(
            "lsdtt::find_marsh_platform: DEM metadata are invalid"
        );
    }
    meta.check_projected_metres();

    if (dem.get_data().empty())
    {
        throw std::runtime_error(
            "lsdtt::find_marsh_platform: DEM raster is empty"
        );
    }

    if (supplied_slope != nullptr && !same_grid(meta, supplied_slope->get_meta()))
    {
        throw std::invalid_argument(
            "lsdtt::find_marsh_platform: supplied slope raster is not co-registered with the DEM"
        );
    }

    marsh_platform_result result;
    const std::vector<uint8_t> dem_valid = make_valid_domain_mask_(dem);

    // -------------------------------------------------------------------------
    // Slope
    // -------------------------------------------------------------------------
    if (supplied_slope != nullptr)
    {
        result.slope = *supplied_slope;
    }
    else
    {
        result.slope = compute_internal_slope_(dem);
    }

    // -------------------------------------------------------------------------
    // Relative reliefs
    // -------------------------------------------------------------------------
    result.dem_relative_relief   = compute_relative_relief_(dem, request.n_threads);
    result.slope_relative_relief = compute_relative_relief_(result.slope, request.n_threads);

    // -------------------------------------------------------------------------
    // Search space
    // -------------------------------------------------------------------------
    {
        const search_space_result_ ss =
            compute_search_space_(
                result.dem_relative_relief,
                result.slope_relative_relief,
                request);

        result.search_crossover = ss.crossover;
        result.search_space     = ss.search_space;
        result.search_histogram = ss.histogram;
        result.search_space_threshold_value = ss.histogram.threshold_value;
    }

    // -------------------------------------------------------------------------
    // Scarp seeds
    // -------------------------------------------------------------------------
    const std::vector<uint8_t> scarp_seed_mask =
        detect_scarp_seeds_(result.slope, result.search_space, request.n_threads);

    result.scarp_seed_mask =
        binary_raster_from_mask_(scarp_seed_mask, dem_valid, meta);

    // -------------------------------------------------------------------------
    // Scarp growth
    // -------------------------------------------------------------------------
    std::vector<float> scarp_order(meta.n_pixels(), 0.0f);
    grow_scarps_(
        result.slope,
        result.search_space,
        scarp_seed_mask,
        request.max_scarp_growth_iterations,
        scarp_order);

    result.scarp_order_raster =
        raster_from_values_on_valid_domain_(
            std::move(scarp_order),
            dem_valid,
            meta);

    // -------------------------------------------------------------------------
    // Scarp cleaning
    // -------------------------------------------------------------------------
    const std::vector<uint8_t> scarp_mask =
        clean_scarps_(
            result.scarp_order_raster.get_data(),
            dem,
            request);

    result.scarp_mask =
        binary_raster_from_mask_(scarp_mask, dem_valid, meta);

    // -------------------------------------------------------------------------
    // Platform seeding
    // -------------------------------------------------------------------------
    std::vector<float> platform_order =
        seed_platform_from_scarps_(scarp_mask, dem, request.n_threads);

    {
        std::vector<uint8_t> seed_bin(platform_order.size(), uint8_t{0});
        for (std::size_t i = 0; i < platform_order.size(); ++i)
        {
            if (platform_order[i] > 0.5f) seed_bin[i] = uint8_t{1};
        }
        result.platform_seed_mask =
            binary_raster_from_mask_(seed_bin, dem_valid, meta);
    }

    // -------------------------------------------------------------------------
    // Platform growth
    // -------------------------------------------------------------------------
    const std::vector<float> local_dem_max_11x11 =
        precompute_local_dem_max_11x11_(dem, request.n_threads);

    grow_platform_(
        platform_order,
        scarp_mask,
        dem,
        local_dem_max_11x11,
        request.max_platform_growth_iterations);

    {
        std::vector<uint8_t> provisional_bin(platform_order.size(), uint8_t{0});
        for (std::size_t i = 0; i < platform_order.size(); ++i)
        {
            if (platform_order[i] > 0.5f) provisional_bin[i] = uint8_t{1};
        }
        result.platform_provisional =
            binary_raster_from_mask_(provisional_bin, dem_valid, meta);
    }

    // -------------------------------------------------------------------------
    // First cutoff
    // -------------------------------------------------------------------------
    cutoff_result_ cut1 =
        compute_cutoff_histogram_(
            platform_order,
            dem,
            result.search_space,
            request.histogram_n_bins,
            request.cutoff_histogram_run_length_bins);

    result.cutoff_histogram       = cut1.diag;
    result.cutoff_elevation_metres = cut1.diag.threshold_value;

    apply_cutoff_(
        platform_order,
        dem,
        cut1.diag.threshold_value,
        request.n_threads);

    if (request.fill_high_gaps)
    {
        fill_high_gaps_(
            platform_order,
            scarp_mask,
            dem,
            cut1.dominant_peak_value,
            request.n_threads);
    }

    if (request.fill_pools)
    {
        fill_pools_(
            platform_order,
            scarp_mask,
            dem);
    }

    if (request.fill_small_holes)
    {
        fill_small_holes_(
            platform_order,
            scarp_mask,
            dem,
            request.n_threads);
    }

    // -------------------------------------------------------------------------
    // Reinsert scarps
    // -------------------------------------------------------------------------
    std::vector<uint8_t> final_scarps = scarp_mask;

    if (request.preserve_scarps_in_final_platform)
    {
        final_scarps =
            clean_scarps_from_platform_(
                scarp_mask,
                platform_order,
                dem,
                request.scarp_clean_min_neighbour_pixels,
                request.n_threads);

        for (std::size_t i = 0; i < platform_order.size(); ++i)
        {
            if (final_scarps[i] != 0)
            {
                platform_order[i] = std::max(platform_order[i], 110.0f);
            }
        }

        if (request.fill_high_gaps)
        {
            cutoff_result_ tmp_cut =
                compute_cutoff_histogram_(
                    platform_order,
                    dem,
                    result.search_space,
                    request.histogram_n_bins,
                    request.cutoff_histogram_run_length_bins);

            fill_high_gaps_(
                platform_order,
                final_scarps,
                dem,
                tmp_cut.dominant_peak_value,
                request.n_threads);
        }

        if (request.fill_pools)
        {
            fill_pools_(
                platform_order,
                final_scarps,
                dem);
        }
    }

    // -------------------------------------------------------------------------
    // Final cutoff
    // -------------------------------------------------------------------------
    cutoff_result_ cut2 =
        compute_cutoff_histogram_(
            platform_order,
            dem,
            result.search_space,
            request.histogram_n_bins,
            request.cutoff_histogram_run_length_bins);

    result.cutoff_histogram        = cut2.diag;
    result.cutoff_elevation_metres = cut2.diag.threshold_value;

    apply_cutoff_(
        platform_order,
        dem,
        cut2.diag.threshold_value,
        request.n_threads);

    // -------------------------------------------------------------------------
    // Final binary platform mask
    // -------------------------------------------------------------------------
    {
        std::vector<uint8_t> platform_bin(platform_order.size(), uint8_t{0});
        for (std::size_t i = 0; i < platform_order.size(); ++i)
        {
            if (platform_order[i] > 0.5f) platform_bin[i] = uint8_t{1};
        }

        result.platform_mask =
            binary_raster_from_mask_(platform_bin, dem_valid, meta);
    }

    if (request.min_platform_component_pixels > 1)
    {
        result.platform_mask =
            filter_binary_components_by_size(
                result.platform_mask,
                request.min_platform_component_pixels,
                pixel_connectivity::d8);
    }

    result.platform_label_raster =
        label_binary_components_by_size(
            result.platform_mask,
            std::max(1, request.min_platform_component_pixels),
            pixel_connectivity::d8);

    return result;
}

} // namespace lsdtt
