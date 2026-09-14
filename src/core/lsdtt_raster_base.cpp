// lsdtt_raster_base.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html

#include "lsdtt_raster_base.hpp"
#include "lsdtt_backend_factory.hpp"

#include <memory>

namespace lsdtt
{

std::unique_ptr<lsdtt_raster_base> make_raster_base(
    const lsdtt_data_loader& loader
)
{
    return make_raster_engine(loader, backend_type::auto_detect);
}

} // namespace lsdtt
