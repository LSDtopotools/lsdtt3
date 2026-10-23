#include <cstdint>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>

#include <cpu/lsdtt_raster_synthetic_surfaces_cpu.hpp>
#include <cpu/lsdtt_synthetic_network_cpu.hpp>

#include "binding_common.hpp"

namespace lsdtt3_python {
namespace {

Raster as_raster(lsdtt::lsdtt_raster_result result)
{
    return Raster(lsdtt::lsdtt_raster_cpu(
        std::move(result.data), std::move(result.mask), std::move(result.meta)));
}

nb::object array_2d(const std::vector<std::int64_t>& values, int n_rows, int n_cols)
{
    nb::object numpy = nb::module_::import_("numpy");
    return numpy.attr("array")(values, nb::arg("dtype") = "int64")
        .attr("reshape")(nb::make_tuple(n_rows, n_cols));
}

lsdtt::lsdtt_range_edge range_edge(const std::string& edge)
{
    if (edge == "north") return lsdtt::lsdtt_range_edge::north;
    if (edge == "south") return lsdtt::lsdtt_range_edge::south;
    if (edge == "east") return lsdtt::lsdtt_range_edge::east;
    if (edge == "west") return lsdtt::lsdtt_range_edge::west;
    throw std::invalid_argument("edge must be 'north', 'south', 'east', or 'west'");
}

lsdtt::lsdtt_strip_channel_config strip_channels(const std::string& value)
{
    if (value == "both_ends") return lsdtt::lsdtt_strip_channel_config::both_ends;
    if (value == "single_end_low") return lsdtt::lsdtt_strip_channel_config::single_end_low;
    if (value == "single_end_high") return lsdtt::lsdtt_strip_channel_config::single_end_high;
    throw std::invalid_argument("channel_config must be 'both_ends', 'single_end_low', or 'single_end_high'");
}

lsdtt::lsdtt_strip_init_profile strip_profile(const std::string& value)
{
    if (value == "flat") return lsdtt::lsdtt_strip_init_profile::flat;
    if (value == "analytic_parabola") return lsdtt::lsdtt_strip_init_profile::analytic_parabola;
    if (value == "raised_plateau") return lsdtt::lsdtt_strip_init_profile::raised_plateau;
    if (value == "perturbed") return lsdtt::lsdtt_strip_init_profile::perturbed;
    if (value == "linear_ramp") return lsdtt::lsdtt_strip_init_profile::linear_ramp;
    throw std::invalid_argument("init_profile must be 'flat', 'analytic_parabola', 'raised_plateau', 'perturbed', or 'linear_ramp'");
}

Raster synthetic_noise(const lsdtt::lsdtt_raster_meta& meta, const std::string& kind,
                       int seed, double mean, double stddev, int blur_passes,
                       double roughness, double beta, int n_threads)
{
    if (!meta.is_valid()) throw std::invalid_argument("metadata must describe a non-empty raster");
    if (kind == "uniform")
        return Raster(lsdtt::generate_uniform_random_noise_surface(meta, seed, mean - stddev, mean + stddev, n_threads));
    if (kind == "red")
        return Raster(lsdtt::generate_red_noise_surface(meta, seed, mean, stddev, blur_passes, n_threads));
    if (kind == "blue")
        return Raster(lsdtt::generate_blue_noise_surface(meta, seed, mean, stddev, blur_passes, n_threads));
    if (kind == "diamond_square")
        return Raster(lsdtt::generate_diamond_square_surface(meta, seed, roughness, 1.0, mean, stddev, n_threads));
    if (kind == "power_law")
        return Raster(lsdtt::generate_power_law_fractal_surface_fft(meta, seed, beta, mean, stddev, n_threads));
    throw std::invalid_argument("kind must be 'uniform', 'red', 'blue', 'diamond_square', or 'power_law'");
}

nb::dict synthetic_mountain_range(const lsdtt::lsdtt_raster_meta& meta, const std::string& geometry,
                                  const std::string& edge, double ridge_elevation,
                                  double piedmont_elevation, double front_fraction,
                                  double noise_amplitude, int seed, int n_threads)
{
    lsdtt::lsdtt_mountain_range_surface result;
    if (geometry == "linear")
        result = lsdtt::generate_linear_range_surface(meta, ridge_elevation, piedmont_elevation,
            0.5, front_fraction, noise_amplitude, noise_amplitude * 0.1, true, 0.55, seed, n_threads);
    else if (geometry == "halfspace")
        result = lsdtt::generate_halfspace_range_surface(meta, range_edge(edge), ridge_elevation,
            piedmont_elevation, front_fraction, noise_amplitude, noise_amplitude * 0.1, 0.55, seed, n_threads);
    else if (geometry == "circular")
        result = lsdtt::generate_circular_range_surface(meta, ridge_elevation, piedmont_elevation,
            0.5, 0.5, front_fraction, front_fraction, noise_amplitude, noise_amplitude * 0.1, 0.55, seed, n_threads);
    else if (geometry == "arc")
        result = lsdtt::generate_arc_range_surface(meta, range_edge(edge), ridge_elevation,
            piedmont_elevation, front_fraction, 0.5, 1.0, noise_amplitude,
            noise_amplitude * 0.1, 0.55, seed, n_threads);
    else
        throw std::invalid_argument("geometry must be 'linear', 'halfspace', 'circular', or 'arc'");
    nb::dict out; out["elevation"] = Raster(std::move(result.dem)); out["mountain_mask"] = Raster(std::move(result.front_mask)); return out;
}

nb::dict synthetic_hillslope_strip(const lsdtt::lsdtt_raster_meta& meta,
                                   const std::string& channel_config,
                                   const std::string& init_profile,
                                   bool cross_axis_is_x, double baselevel_elevation,
                                   double half_width, double uplift_rate,
                                   double diffusivity, double critical_slope,
                                   bool nonlinear, double initial_amplitude,
                                   int seed, int n_threads)
{
    const auto result = lsdtt::generate_hillslope_strip_surface(meta, strip_channels(channel_config),
        strip_profile(init_profile), cross_axis_is_x, baselevel_elevation, half_width,
        uplift_rate, diffusivity, critical_slope, nonlinear, initial_amplitude, seed, n_threads);
    nb::dict out; out["elevation"] = Raster(result.dem); out["channel_mask"] = Raster(result.channel_mask); return out;
}

class SyntheticNetworkResult {
public:
    SyntheticNetworkResult(lsdtt::lsdtt_synthetic_network_cpu network)
        : network_(std::move(network)) {}
    Raster elevation() const { return as_raster(network_.get_elevation_raster()); }
    Raster drainage_area_cells() const { return as_raster(network_.get_drainage_area_cells_raster()); }
    double initial_energy() const { return network_.initial_energy(); }
    double final_energy() const { return network_.final_energy(); }
    nb::object receivers() const { return array_2d(network_.receivers(), network_.n_rows(), network_.n_cols()); }
private:
    lsdtt::lsdtt_synthetic_network_cpu network_;
};

SyntheticNetworkResult synthetic_network(const lsdtt::lsdtt_raster_meta& meta,
                                         const std::string& outlet_mode, bool periodic_columns,
                                         int outlet_spacing, double gamma, long long n_iterations,
                                         double temperature_scale, double m_exponent, double n_exponent,
                                          double relief_metres, std::uint64_t seed)
{
    if (!meta.is_valid()) throw std::invalid_argument("metadata must describe a non-empty raster");
    lsdtt::lsdtt_ocn_options options;
    options.n_rows = meta.n_rows; options.n_cols = meta.n_cols; options.cell_size = meta.cell_size;
    options.use_supplied_meta = true; options.supplied_meta = meta;
    options.periodic_columns = periodic_columns; options.outlet_spacing = outlet_spacing;
    options.gamma = gamma; options.n_iter = n_iterations; options.temperature_scale = temperature_scale;
    options.m_exponent = static_cast<float>(m_exponent); options.n_exponent = static_cast<float>(n_exponent);
    options.relief_mode = true; options.target_relief = static_cast<float>(relief_metres); options.seed = seed;
    if (outlet_mode == "single_south") options.outlet_mode = lsdtt::ocn_outlet_mode::single_south;
    else if (outlet_mode == "all_coast") options.outlet_mode = lsdtt::ocn_outlet_mode::all_coast;
    else if (outlet_mode == "north_south") options.outlet_mode = lsdtt::ocn_outlet_mode::north_south;
    else throw std::invalid_argument("outlet_mode must be 'single_south', 'all_coast', or 'north_south'");
    lsdtt::lsdtt_synthetic_network_cpu network(options); network.run(); return SyntheticNetworkResult(std::move(network));
}

} // namespace

void bind_synthetic_workflows(nb::module_& module)
{
    module.def("synthetic_noise", &synthetic_noise, nb::arg("metadata"), nb::arg("kind") = "red",
        nb::arg("seed") = 1, nb::arg("mean") = 0.0, nb::arg("stddev") = 1.0,
        nb::arg("blur_passes") = 6, nb::arg("roughness") = 0.55, nb::arg("beta") = 2.0,
        nb::arg("n_threads") = 0, R"doc(Generate a deterministic synthetic noise raster.

Parameters
----------
metadata : RasterMeta
    Output grid definition.
kind : {"uniform", "red", "blue", "diamond_square", "power_law"}, default="red"
    Surface generator.
seed : int, default=1
    Reproducible random seed.
mean, stddev : float, default=0, 1
    Target mean and spread; for uniform noise, values span mean ± stddev.
blur_passes, roughness, beta : int, float, float
    Generator-specific correlation or spectral controls.
n_threads : int, default=0
    OpenMP thread count.

Returns
-------
Raster
    Newly allocated synthetic raster.

Authors
-------
B. Gailleton
)doc");
    module.def("synthetic_mountain_range", &synthetic_mountain_range, nb::arg("metadata"),
        nb::arg("geometry") = "linear", nb::arg("edge") = "north", nb::arg("ridge_elevation") = 2000.,
        nb::arg("piedmont_elevation") = 0., nb::arg("front_fraction") = .2,
        nb::arg("noise_amplitude") = 100., nb::arg("seed") = 1, nb::arg("n_threads") = 0,
        R"doc(Generate a synthetic mountain range and its known mountain mask.

Parameters
----------
metadata : RasterMeta
    Output grid (``RasterMeta.from_shape``).
geometry : {"linear", "halfspace", "circular", "arc"}, default="linear"
    ``linear``: a north-south Gaussian ridge through the centre with
    piedmonts on both sides; ``halfspace``: mountains against one edge;
    ``circular``: an elliptical dome in the middle; ``arc``: mountains
    against one edge with a front that curves into the piedmont.
edge : {"north", "south", "east", "west"}, default="north"
    Edge the mountains are against (``halfspace`` and ``arc`` only).
ridge_elevation : float, default=2000.0
    Crest elevation (m).
piedmont_elevation : float, default=0.0
    Base (piedmont) elevation (m).
front_fraction : float, default=0.2
    Size of the mountain zone as a fraction of the grid extent: the ridge
    half-width (``linear``), the ellipse semi-axes (``circular``) or the
    front's distance from the edge (``halfspace``, ``arc``).
noise_amplitude : float, default=100.0
    Amplitude (m) of the fractal noise added to the surface.
seed : int, default=1
    Random seed.
n_threads : int, default=0
    OpenMP thread count; 0 uses the runtime default.

Returns
-------
dict
    ``"elevation"`` (Raster, m) and ``"mountain_mask"`` (Raster, 1 =
    mountain, 0 = piedmont).
)doc");
    module.def("synthetic_hillslope_strip", &synthetic_hillslope_strip, nb::arg("metadata"),
        nb::arg("channel_config") = "both_ends", nb::arg("init_profile") = "flat",
        nb::arg("cross_axis_is_x") = true, nb::arg("baselevel_elevation") = 0., nb::arg("half_width") = 0.,
        nb::arg("uplift_rate") = 1e-3, nb::arg("diffusivity") = 1e-2, nb::arg("critical_slope") = .8,
        nb::arg("nonlinear") = false, nb::arg("initial_amplitude") = 50., nb::arg("seed") = 1,
        nb::arg("n_threads") = 0, R"doc(Generate a one-dimensional hillslope test DEM and fixed channel mask.

Elevation varies across one axis only, between channel lines at one or both
ends, for testing hillslope solvers against analytic profiles.

Parameters
----------
metadata : RasterMeta
    Output grid.
channel_config : {"both_ends", "single_end_low", "single_end_high"}, default="both_ends"
    Channels at both ends (divide in the middle) or at one end (divide at
    the other).
init_profile : {"flat", "analytic_parabola", "raised_plateau", "perturbed", "linear_ramp"}, default="flat"
    Initial hillslope surface: flat at base level, the analytic steady
    profile, a raised block, flat plus noise, or a triangular ridge.
cross_axis_is_x : bool, default=True
    True: elevation varies east-west and channels are north-south columns;
    False: it varies north-south and channels are rows.
baselevel_elevation : float, default=0.0
    Elevation (m) of the channel cells.
half_width : float, default=0.0
    Hillslope half-width L (m, divide to channel) of the analytic profiles;
    <= 0 derives it from the grid.
uplift_rate : float, default=0.001
    Uplift rate U (m/yr) of the analytic profile.
diffusivity : float, default=0.01
    Diffusivity D (m^2/yr) of the analytic profile.
critical_slope : float, default=0.8
    Critical gradient S_c (m/m) of the nonlinear analytic profile.
nonlinear : bool, default=False
    Use the nonlinear (Roering 2007) analytic profile instead of the linear
    parabola ``z = U / (2 D) (L^2 - x^2)``.
initial_amplitude : float, default=50.0
    Height (m) of the raised plateau, noise or ramp crest.
seed : int, default=1
    Random seed for ``"perturbed"``.
n_threads : int, default=0
    OpenMP thread count; 0 uses the runtime default.

Returns
-------
dict
    ``"elevation"`` (Raster, m) and ``"channel_mask"`` (Raster, 0/1).
)doc");
    nb::class_<SyntheticNetworkResult>(module, "SyntheticNetworkResult",
        "In-memory optimal-channel-network outputs.")
        .def_prop_ro("elevation", &SyntheticNetworkResult::elevation, "Steady-state elevation Raster (m) of the network.")
        .def_prop_ro("drainage_area_cells", &SyntheticNetworkResult::drainage_area_cells, "Drainage area in cells (Raster).")
        .def_prop_ro("initial_energy", &SyntheticNetworkResult::initial_energy, "Network energy sum(A^gamma L) before optimisation.")
        .def_prop_ro("final_energy", &SyntheticNetworkResult::final_energy, "Network energy after optimisation.")
        .def_prop_ro("receivers", &SyntheticNetworkResult::receivers, nb::sig("def receivers(self) -> numpy.ndarray"),
            "Int array (n_rows, n_cols) of each cell's receiver flat index.");
    module.def("synthetic_network", &synthetic_network, nb::arg("metadata"),
        nb::arg("outlet_mode") = "single_south", nb::arg("periodic_columns") = false,
        nb::arg("outlet_spacing") = 1, nb::arg("gamma") = .5, nb::arg("n_iterations") = -1,
        nb::arg("temperature_scale") = 1., nb::arg("m_exponent") = .45, nb::arg("n_exponent") = 1.,
        nb::arg("relief_metres") = 1000., nb::arg("seed") = 1,
        R"doc(Optimise a synthetic D8 optimal channel network and derive its chi steady-state DEM.

An optimal channel network (OCN) is found by simulated annealing of the
energy ``sum(A**gamma * L)`` over D8 networks, and a steady-state
stream-power DEM is built on it from chi.

Parameters
----------
metadata : RasterMeta
    Output grid.
outlet_mode : {"single_south", "all_coast", "north_south"}, default="single_south"
    One outlet at the centre of the south edge; every edge pixel an outlet;
    or outlets along the north and south edges.
periodic_columns : bool, default=False
    Wrap the grid east-west (route the DEM with ``p`` on east and west).
outlet_spacing : int, default=1
    ``"north_south"`` only: column spacing of the outlets (larger values
    force wider basins).
gamma : float, default=0.5
    Energy exponent on drainage area.
n_iterations : int, default=-1
    Annealing iterations; <= 0 uses 40 x n_rows x n_cols.
temperature_scale : float, default=1.0
    Multiplier on the annealing temperature schedule.
m_exponent, n_exponent : float, default=0.45, 1.0
    Stream-power exponents used to build the DEM from chi.
relief_metres : float, default=1000.0
    Relief (m) of the output DEM.
seed : int, default=1
    Random seed.

Returns
-------
SyntheticNetworkResult
)doc");
}
} // namespace lsdtt3_python
