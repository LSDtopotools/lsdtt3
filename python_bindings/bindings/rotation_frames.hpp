// SPDX-License-Identifier: AGPL-3.0-or-later
//
// RotationFrame: an oblique analysis frame and its rotate-back target, held
// as the core lsdtt_rotation_frame_file so fetch_dem can use it too.

#pragma once

#include <filesystem>
#include <string>
#include <utility>

#include <core/lsdtt_rotation_frame_file.hpp>

#include "binding_common.hpp"

namespace lsdtt3_python {

class RotationFrame
{
public:
    explicit RotationFrame(lsdtt::lsdtt_rotation_frame_file file) : file_(std::move(file)) {}

    static RotationFrame from_feature(const std::string& path, double width, double end_buffer,
                                      double grid_spacing, const std::string& layer,
                                      nb::object source, int target_epsg);
    static RotationFrame from_valid_data(const Raster& source, double end_buffer,
                                         double grid_spacing);
    static RotationFrame from_raster(const Raster& raster);
    static RotationFrame read(const std::string& path);

    void write(const std::string& path) const;
    double grid_convergence_degrees(double x, double y) const;

    const lsdtt::lsdtt_rotation_frame_file& native() const noexcept { return file_; }

private:
    lsdtt::lsdtt_rotation_frame_file file_;
};

Raster rotate_to_frame(const Raster& raster, const RotationFrame& frame,
                       const std::string& resampling);
Raster rotate_back(const Raster& raster, const RotationFrame& frame,
                   const std::string& resampling, int target_epsg);

} // namespace lsdtt3_python
