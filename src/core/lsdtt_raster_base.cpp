// lsdtt_raster_base.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html

#include "lsdtt_raster_base.hpp"
#include "lsdtt_backend_factory.hpp"

namespace lsdtt
{

std::unique_ptr<lsdtt_raster_base> make_raster_base(
    const lsdtt_data_loader& loader
)
{
    return make_raster_engine(loader, backend_type::auto_detect);
}

} // namespace lsdtt
