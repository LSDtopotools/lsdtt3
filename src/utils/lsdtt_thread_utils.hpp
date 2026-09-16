// lsdtt_thread_utils.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// Thread count resolution helper.
//
// Used by CPU kernels that accept an n_threads parameter.
// If n_threads <= 0, the function returns the hardware-detected maximum.
// -----------------------------------------------------------------------------

#pragma once

#include <thread>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace lsdtt {

// Returns n_threads if > 0, otherwise returns the OpenMP or hardware maximum.
// Result is always >= 1.
inline int resolve_thread_count(int n_threads)
{
    int threads = n_threads;

    if (threads <= 0)
    {
#ifdef _OPENMP
        threads = omp_get_max_threads();
#else
        threads = static_cast<int>(std::thread::hardware_concurrency());
#endif
    }

    if (threads < 1) threads = 1;
    return threads;
}

} // namespace lsdtt
