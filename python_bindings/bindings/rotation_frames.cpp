// SPDX-License-Identifier: AGPL-3.0-or-later
//
// Oblique analysis frames: RotationFrame, rotate_to_frame and rotate_back.
// Every operation calls the same core routines as lsdtt-fetch-raster and
// lsdtt-raster-preprocessing (lsdtt_oblique_frame, lsdtt_rotation_frame_file,
// lsdtt_raster_warp); in-memory rasters pass through /vsimem/ GeoTIFFs.

#include <cstdint>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>

#include <cpl_vsi.h>
#include <Python.h>

#include <core/lsdtt_gdal_init.hpp>
#include <core/lsdtt_oblique_frame.hpp>
#include <core/lsdtt_raster_warp.hpp>
#include <core/lsdtt_rotation_frame_file.hpp>
#include <core/lsdtt_true_north.hpp>

#include "binding_common.hpp"
#include "rotation_frames.hpp"

namespace lsdtt3_python {

namespace {

// Default cell size of a frame built for a fetch (fetch_dem's default).
constexpr double default_fetch_grid_spacing = 30.0;

// A raster written to a /vsimem/ GeoTIFF for the life of this object.
class VsimemRaster
{
public:
    explicit VsimemRaster(const Raster& raster)
        : path_("/vsimem/lsdtt3-rotation-" + std::to_string(++counter_) + ".tif")
    {
        lsdtt::gdal_init();
        raster.write(path_, "gtiff");
    }
    VsimemRaster()
        : path_("/vsimem/lsdtt3-rotation-" + std::to_string(++counter_) + ".tif")
    {
        lsdtt::gdal_init();
    }
    ~VsimemRaster() { VSIUnlink(path_.c_str()); }
    VsimemRaster(const VsimemRaster&) = delete;
    VsimemRaster& operator=(const VsimemRaster&) = delete;
    const std::string& path() const noexcept { return path_; }

private:
    static inline std::uint64_t counter_ = 0;
    std::string path_;
};

void warn_all(const std::vector<std::string>& warnings)
{
    for (const auto& w : warnings)
        if (PyErr_WarnEx(PyExc_UserWarning, w.c_str(), 2) < 0) throw nb::python_error();
}

const char* definition_name(lsdtt::oblique_frame_definition d)
{
    switch (d)
    {
    case lsdtt::oblique_frame_definition::two_points:     return "two_points";
    case lsdtt::oblique_frame_definition::line:           return "line";
    case lsdtt::oblique_frame_definition::polygon:        return "polygon";
    case lsdtt::oblique_frame_definition::crs:            return "crs";
    case lsdtt::oblique_frame_definition::data_footprint: return "data_footprint";
    }
    return "unknown";
}

// The rotate-back target of a frame built for a fetch: the analysis CRS of
// the frame centre (or target_epsg) at the frame cell size on the lattice
// anchored at (0, 0), as lsdtt::fetch_dem writes it.
lsdtt::lsdtt_rotation_target fetch_target(const lsdtt::lsdtt_oblique_frame& frame, int target_epsg)
{
    lsdtt::lsdtt_rotation_target t;
    t.epsg = target_epsg > 0 ? target_epsg : lsdtt::analysis_epsg_for_frame_centre(frame);
    t.cell_size = frame.cell_size;
    return t;
}

RotationFrame frame_for_source(const lsdtt::lsdtt_frame_from_feature_result& fr,
                               const Raster& source, const std::string& feature_fname)
{
    lsdtt::lsdtt_rotation_frame_file f;
    f.ok = true;
    f.version = lsdtt::rotation_frame_file_version;
    f.frame = fr.frame;
    f.target = lsdtt::rotation_target_for_raster(source.native().get_meta(), "", fr.frame);
    f.provenance.created_by = "lsdtt3 (Python)";
    f.provenance.feature_fname = feature_fname;
    f.provenance.feature_layer = fr.layer_name;
    f.warnings = fr.warnings;
    return RotationFrame(std::move(f));
}

} // namespace

RotationFrame RotationFrame::from_feature(const std::string& path, double width,
                                          double end_buffer, double grid_spacing,
                                          const std::string& layer, nb::object source,
                                          int target_epsg)
{
    lsdtt::gdal_init();
    const bool have_source = !source.is_none();
    if (have_source && target_epsg != 0)
        throw std::invalid_argument(
            "target_epsg cannot be given with source: the rotate-back target is the source "
            "raster's own grid");
    const std::string feature_fname = std::filesystem::path(path).filename().string();
    if (have_source)
    {
        const Raster& dem = nb::cast<const Raster&>(source);
        const auto& meta = dem.native().get_meta();
        if (meta.crs_wkt.empty())
            throw std::invalid_argument("source has no CRS, so it cannot be rotated into a frame");
        const double cell = lsdtt::frame_cell_size_for_raster(meta, grid_spacing);
        const auto fr = lsdtt::frame_from_feature_file(
            path, layer, width, end_buffer, cell, lsdtt::base_geographic_wkt_of(meta.crs_wkt));
        warn_all(fr.warnings);
        if (!fr.ok) throw std::invalid_argument(fr.message);
        RotationFrame out = frame_for_source(fr, dem, feature_fname);
        return out;
    }
    const double cell = grid_spacing > 0.0 ? grid_spacing : default_fetch_grid_spacing;
    const auto fr = lsdtt::frame_from_feature_file(path, layer, width, end_buffer, cell);
    warn_all(fr.warnings);
    if (!fr.ok) throw std::invalid_argument(fr.message);
    lsdtt::lsdtt_rotation_frame_file f;
    f.ok = true;
    f.version = lsdtt::rotation_frame_file_version;
    f.frame = fr.frame;
    f.target = fetch_target(fr.frame, target_epsg);
    f.provenance.created_by = "lsdtt3 (Python)";
    f.provenance.feature_fname = feature_fname;
    f.provenance.feature_layer = fr.layer_name;
    f.warnings = fr.warnings;
    return RotationFrame(std::move(f));
}

RotationFrame RotationFrame::from_valid_data(const Raster& source, double end_buffer,
                                             double grid_spacing)
{
    const auto& meta = source.native().get_meta();
    if (meta.crs_wkt.empty())
        throw std::invalid_argument("source has no CRS, so it cannot be rotated into a frame");
    const VsimemRaster file(source);
    const double cell = lsdtt::frame_cell_size_for_raster(meta, grid_spacing);
    const auto fr = lsdtt::frame_from_valid_data(file.path(), end_buffer, cell);
    warn_all(fr.warnings);
    if (!fr.ok) throw std::invalid_argument(fr.message);
    return frame_for_source(fr, source, "");
}

RotationFrame RotationFrame::from_raster(const Raster& raster)
{
    const VsimemRaster file(raster);
    lsdtt::lsdtt_rotation_frame_file f = lsdtt::rotation_frame_from_raster(file.path());
    if (!f.ok) throw std::invalid_argument(f.message);
    warn_all(f.warnings);
    return RotationFrame(std::move(f));
}

RotationFrame RotationFrame::read(const std::string& path)
{
    lsdtt::gdal_init();
    lsdtt::lsdtt_rotation_frame_file f = lsdtt::read_rotation_frame(path);
    if (!f.ok) throw std::invalid_argument(f.message);
    warn_all(f.warnings);
    return RotationFrame(std::move(f));
}

void RotationFrame::write(const std::string& path) const
{
    lsdtt::lsdtt_rotation_provenance p = file_.provenance;
    p.lsdtt3_version.clear();   // this build's version and the current time
    p.created_utc.clear();
    lsdtt::write_rotation_frame(path, file_.frame, file_.target, p);
}

double RotationFrame::grid_convergence_degrees(double x, double y) const
{
    const lsdtt::lsdtt_grid_convergence c(file_.frame.to_wkt());
    return c.degrees_at(x, y);
}

Raster rotate_to_frame(const Raster& raster, const RotationFrame& frame,
                       const std::string& resampling)
{
    if (raster.native().get_meta().crs_wkt.empty())
        throw std::invalid_argument("raster has no CRS, so it cannot be rotated into a frame");
    const VsimemRaster in(raster);
    const VsimemRaster out;
    lsdtt::lsdtt_warp_options opts;
    opts.resampling_method = resampling;
    opts.output_format = lsdtt::raster_output_format::gtiff;
    lsdtt::warp_raster_to_frame(in.path(), frame.native().frame, out.path(), opts);
    return Raster::read(out.path());
}

Raster rotate_back(const Raster& raster, const RotationFrame& frame,
                   const std::string& resampling, int target_epsg)
{
    const VsimemRaster in(raster);
    const VsimemRaster out;
    lsdtt::lsdtt_warp_options opts;
    opts.output_format = lsdtt::raster_output_format::gtiff;
    const lsdtt::lsdtt_rotate_back_result r = lsdtt::rotate_back_raster(
        in.path(), frame.native(), resampling, target_epsg, out.path(), opts);
    warn_all(r.warnings);
    return Raster::read(out.path());
}

void bind_rotation_frames(nb::module_& module)
{
    nb::class_<RotationFrame>(module, "RotationFrame", R"doc(
An oblique analysis frame aligned with a feature, and its rotate-back target.

A frame is a Hotine Oblique Mercator (variant B) CRS whose central line runs
along a fault, river, coastline or elongated range, with the feature along
+x and the frame centre at (0, 0). A raster in the frame is an ordinary
north-up raster in that CRS, so every lsdtt3 analysis runs on it unchanged,
on a grid that holds the strip rather than its bounding box. Bearings from
analyses are then clockwise from grid up (90 = along the feature) unless
``bearings_relative_to_true_north=True`` is passed.

The frame also carries the lattice that ``rotate_back`` writes products onto:
the source DEM's own grid for a frame built with ``source=``, else the UTM
(or polar stereographic) CRS of the frame centre on the lattice a fetch uses.
Bring back products (basins, channels, IDs), never the DEM.

Get one with ``RotationFrame.from_feature`` (a line, two points or polygon
file), ``RotationFrame.from_valid_data`` (the footprint of a raster's valid
cells), ``RotationFrame.read`` (a ``*_rotation_frame.txt`` written by a CLI
or by ``write``), or ``RotationFrame.from_raster`` (rebuilt from a frame
raster's own CRS when the frame file is lost). Typical use::

    frame = lsdtt3.RotationFrame.from_feature("fault.gpkg", width=10000.0,
                                              grid_spacing=30.0)
    dem = lsdtt3.fetch_dem(frame=frame, grid_spacing=30.0)
    frame.write("fault_rotation_frame.txt")
    # ... analysis on dem ...
    basins_utm = lsdtt3.rotate_back(basins, frame, resampling="near")

Property names spell out the units of the frame-file keys: ``centre_lon`` is
``centre_longitude``, ``azimuth_deg`` is ``azimuth_degrees``, ``width_m`` is
``width_metres``, ``target_crs_epsg`` is ``target_epsg``, and so on.
)doc")
        .def_static("from_feature", &RotationFrame::from_feature, nb::arg("path"),
            nb::arg("width") = 0.0, nb::arg("end_buffer") = -1.0,
            nb::arg("grid_spacing") = 0.0, nb::arg("layer") = "",
            nb::arg("source") = nb::none(), nb::arg("target_epsg") = 0, R"doc(
Build the frame of a line, two-point or polygon feature file.

A line (or a multi-part line) uses its two furthest end points as the axis;
two points are the axis ends; a polygon uses its minimum-area enclosing
rectangle, longest side as the axis, and takes its width from the polygon.
The same routine as ``fetch_from_linear_feature`` in lsdtt-fetch-raster
and ``rotate_to_linear_feature`` in lsdtt-raster-preprocessing.

Parameters
----------
path : str
    Feature file (any OGR vector: GeoPackage, GeoJSON, shapefile, ...), in
    any CRS.
width : float, default=0.0
    Full width of the strip across the axis, in metres. Required (> 0) for
    a line or two points; ignored, with a warning, for a polygon.
end_buffer : float, default=-1.0
    Length added beyond each axis end, in metres (on all four sides for a
    polygon). Negative means 10 cells.
grid_spacing : float, default=0.0
    Frame cell size in metres. Without ``source``, 0 means 30 m (the
    ``fetch_dem`` default). With ``source``, 0 means the source's own cell
    size when it is projected in metres, else the size reprojection would
    choose, as lsdtt-raster-preprocessing does.
layer : str, default=""
    Layer to read; empty reads the first layer (warned if there are several).
source : Raster or None, default=None
    The DEM this frame will rotate (``rotate_to_frame``). The frame is built
    on its geographic datum, and the rotate-back target is its own grid
    (CRS, cell size and lower-left corner) when it is projected in metres.
    None builds a frame for ``fetch_dem`` on WGS84.
target_epsg : int, default=0
    Without ``source`` only: EPSG code of the rotate-back target. 0 picks the
    UTM (or polar stereographic beyond 75 degrees) CRS of the frame centre.

Returns
-------
RotationFrame
    The frame. Reader and builder notes (several layers, an ill-defined
    polygon axis) are issued as ``UserWarning``; a file that cannot be read
    or used raises ``ValueError``.
)doc")
        .def_static("from_valid_data", &RotationFrame::from_valid_data, nb::arg("source"),
            nb::arg("end_buffer") = -1.0, nb::arg("grid_spacing") = 0.0, R"doc(
Build the frame of a raster's valid cells (its minimum-area rectangle).

The polygon of the outer corners of the first and last valid cell of every
row, framed as a polygon feature would be; the frame of
``rotate_to_data_footprint`` in lsdtt-raster-preprocessing. A footprint with
no clear long axis (a square DEM) is warned.

Parameters
----------
source : Raster
    The DEM; it must have a CRS. The rotate-back target is its own grid
    when it is projected in metres.
end_buffer : float, default=-1.0
    Length added on all four sides, in metres; negative means 10 cells.
grid_spacing : float, default=0.0
    Frame cell size in metres; 0 chooses it from ``source`` as
    ``from_feature`` does.

Returns
-------
RotationFrame
    The frame; ``ValueError`` for a raster with no CRS or no valid cells.
)doc")
        .def_static("from_raster", &RotationFrame::from_raster, nb::arg("raster"), R"doc(
Rebuild the frame of a raster that is already in an lsdtt3 oblique frame.

The fallback when the frame file is lost: the frame comes from the raster's
own CRS and extent, and the rotate-back target is the UTM (or polar
stereographic) CRS of the frame centre on the lattice anchored at (0, 0). A
``UserWarning`` says that a source DEM origin off that lattice cannot be
recovered without the frame file.

Parameters
----------
raster : Raster
    A raster in an lsdtt3 oblique Mercator frame.

Returns
-------
RotationFrame
    The rebuilt frame; ``ValueError`` for a raster in any other CRS.
)doc")
        .def_static("read", &RotationFrame::read, nb::arg("path"), R"doc(
Read a rotation frame file (``*_rotation_frame.txt``).

The files written by lsdtt-fetch-raster, lsdtt-raster-preprocessing and
``RotationFrame.write`` are flat YAML mappings; this reader checks the
stored CRS against the stored centre and azimuth and runs the frame's axis
self-check.

Parameters
----------
path : str
    The frame file.

Returns
-------
RotationFrame
    The frame and its rotate-back target; ``ValueError`` if the file is
    missing, newer than this lsdtt3, or inconsistent.
)doc")
        .def("write", &RotationFrame::write, nb::arg("path"), R"doc(
Write the frame and its rotate-back target as a rotation frame file.

The format is the one the CLIs write and read (``rotation_frame_fname`` in
lsdtt-raster-preprocessing), so a frame made in Python can rotate products
back from the command line and the other way round.

Parameters
----------
path : str
    Output file, conventionally ``<stem>_rotation_frame.txt``.

Returns
-------
None
    Nothing; ``RuntimeError`` if the file cannot be written.
)doc")
        .def("grid_convergence_degrees", &RotationFrame::grid_convergence_degrees,
            nb::arg("x"), nb::arg("y"), R"doc(
Angle of true north, clockwise from grid up, at a point of the frame.

A grid bearing ``g`` at ``(x, y)`` is the true bearing ``g - c``. At the
frame centre ``c`` is ``90 - azimuth_degrees``.

Parameters
----------
x, y : float
    Frame coordinates in metres.

Returns
-------
float
    The convergence in degrees.
)doc")
        .def_prop_ro("frame_source", [](const RotationFrame& f) {
            return std::string(definition_name(f.native().frame.definition)); },
            "How the frame was defined: line, two_points, polygon, data_footprint or crs.")
        .def_prop_ro("centre_longitude", [](const RotationFrame& f) { return f.native().frame.centre.longitude; },
            "Longitude of the frame centre, degrees in the base geographic CRS.")
        .def_prop_ro("centre_latitude", [](const RotationFrame& f) { return f.native().frame.centre.latitude; },
            "Latitude of the frame centre, degrees in the base geographic CRS.")
        .def_prop_ro("azimuth_degrees", [](const RotationFrame& f) { return f.native().frame.azimuth_degrees; },
            "Azimuth of the axis at the centre, clockwise from true north, in (-90, 90].")
        .def_prop_ro("axis_length_metres", [](const RotationFrame& f) { return f.native().frame.axis_length_metres; },
            "Geodesic length of the axis between its end points (without the end buffer), m.")
        .def_prop_ro("width_metres", [](const RotationFrame& f) { return f.native().frame.width_metres; },
            "Full width of the strip across the axis, m.")
        .def_prop_ro("end_buffer_metres", [](const RotationFrame& f) { return f.native().frame.end_buffer_metres; },
            "Length added beyond each axis end (and each long side for a polygon), m.")
        .def_prop_ro("buffer_all_sides", [](const RotationFrame& f) { return f.native().frame.buffer_all_sides; },
            "True when the end buffer is also on both long sides (polygon frames).")
        .def_prop_ro("end_points", [](const RotationFrame& f) {
            const auto& fr = f.native().frame;
            return nb::make_tuple(nb::make_tuple(fr.end1.longitude, fr.end1.latitude),
                                  nb::make_tuple(fr.end2.longitude, fr.end2.latitude)); },
            "The axis end points ((lon, lat), (lon, lat)), west one first.")
        .def_prop_ro("cell_size", [](const RotationFrame& f) { return f.native().frame.cell_size; },
            "Frame cell size, m.")
        .def_prop_ro("n_cols", [](const RotationFrame& f) { return f.native().frame.n_cols(); },
            "Number of columns of a raster filling the frame.")
        .def_prop_ro("n_rows", [](const RotationFrame& f) { return f.native().frame.n_rows(); },
            "Number of rows of a raster filling the frame.")
        .def_prop_ro("extent", [](const RotationFrame& f) {
            const auto& fr = f.native().frame;
            return nb::make_tuple(fr.x_min, fr.x_max, fr.y_min, fr.y_max); },
            "Frame extent (x_min, x_max, y_min, y_max) in frame metres.")
        .def_prop_ro("crs_wkt", [](const RotationFrame& f) { return f.native().frame.to_wkt(); },
            "The frame CRS as WKT.")
        .def_prop_ro("proj_string", [](const RotationFrame& f) { return f.native().frame.to_proj(); },
            "The frame CRS as a PROJ string.")
        .def_prop_ro("base_crs_wkt", [](const RotationFrame& f) { return f.native().frame.base_geographic_wkt; },
            "WKT of the geographic CRS (datum) the frame is built on.")
        .def_prop_ro("target_epsg", [](const RotationFrame& f) { return f.native().target.epsg; },
            "EPSG code of the rotate-back target CRS (0 when it is given as WKT).")
        .def_prop_ro("target_crs_wkt", [](const RotationFrame& f) { return f.native().target.crs_wkt; },
            "WKT of the rotate-back target CRS when it has no EPSG code, else empty.")
        .def_prop_ro("target_cell_size", [](const RotationFrame& f) { return f.native().target.cell_size; },
            "Cell size of the rotate-back target grid, in its CRS units.")
        .def_prop_ro("target_origin", [](const RotationFrame& f) {
            return nb::make_tuple(f.native().target.origin_x, f.native().target.origin_y); },
            "A node (x, y) of the rotate-back lattice: the source DEM's lower-left corner, or (0, 0).")
        .def_prop_ro("source_dem_fname", [](const RotationFrame& f) { return f.native().target.source_dem_fname; },
            "File name of the source DEM, when the frame file records one.")
        .def_prop_ro("feature_fname", [](const RotationFrame& f) { return f.native().provenance.feature_fname; },
            "File name of the feature the frame was built from, or empty.")
        .def_prop_ro("feature_layer", [](const RotationFrame& f) { return f.native().provenance.feature_layer; },
            "Layer of the feature file the frame was built from, or empty.")
        .def_prop_ro("warnings", [](const RotationFrame& f) { return f.native().warnings; },
            "Notes from building or reading the frame (also issued as warnings).")
        .def("__repr__", [](const RotationFrame& f) {
            const auto& fr = f.native().frame;
            return "RotationFrame(" + std::string(definition_name(fr.definition)) +
                   ", azimuth " + std::to_string(fr.azimuth_degrees) + " deg, " +
                   std::to_string(fr.n_cols()) + " x " + std::to_string(fr.n_rows()) +
                   " cells of " + std::to_string(fr.cell_size) + " m)"; });

    module.def("rotate_to_frame", &rotate_to_frame, nb::arg("raster"), nb::arg("frame"),
        nb::arg("resampling") = "cubic", R"doc(
Warp a raster into an oblique frame, in one resampling pass from its own CRS.

The output is a north-up raster in ``frame.crs_wkt`` covering the frame
(``frame.n_cols`` x ``frame.n_rows`` cells of ``frame.cell_size``), with the
feature along +x; cells outside the source are nodata. The routine of
``rotate_to_linear_feature`` / ``rotate_to_data_footprint`` in
lsdtt-raster-preprocessing, so the same frame and resampling give the same
raster.

Parameters
----------
raster : Raster
    The DEM (or any raster) to rotate; it must have a CRS.
frame : RotationFrame
    The frame, normally built with ``source=raster`` so that rotate-back
    lands on the raster's own pixels.
resampling : str, default="cubic"
    GDAL warp method: near, bilinear, cubic, cubicspline, lanczos, average,
    mode, max, min, med, q1, q3, sum (``resampling_method`` in
    lsdtt-raster-preprocessing).

Returns
-------
Raster
    The raster in the frame.
)doc");

    module.def("rotate_back", &rotate_back, nb::arg("raster"), nb::arg("frame"),
        nb::arg("resampling") = "bilinear", nb::arg("target_epsg") = 0, R"doc(
Warp a product in an oblique frame back onto the frame's rotate-back grid.

The routine of ``run_rotate_back`` in lsdtt-raster-preprocessing. A raster at
the frame's cell size goes to the target cell size, any other keeps its own;
the output lattice is anchored at the target origin, so after
``rotate_to_frame`` with a frame built from a projected DEM the product lands
on that DEM's own pixels. Rotate back products (basins, channel masks, IDs,
metrics), never the DEM itself: a round trip resamples twice.

Parameters
----------
raster : Raster
    A product in ``frame.crs_wkt``; any other CRS raises ``ValueError``.
frame : RotationFrame
    The frame the product was made in.
resampling : str, default="bilinear"
    GDAL warp method (as ``rotate_to_frame``); use ``"near"`` for
    categorical rasters (basin IDs, masks, junctions, lithology).
target_epsg : int, default=0
    EPSG code overriding the frame's target CRS (``rotate_back_epsg``); a
    CRS other than the frame's target is gridded on the lattice anchored at
    (0, 0), with a warning. 0 uses the frame's target.

Returns
-------
Raster
    The product on the target grid.
)doc");
}

} // namespace lsdtt3_python
