#include <algorithm>
#include <cctype>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>

#include <cpu/lsdtt_raster_spectral_cpu.hpp>

#include "binding_common.hpp"

namespace lsdtt3_python {
namespace {

lsdtt::spectral_window_type parse_window(std::string name)
{
    std::transform(name.begin(), name.end(), name.begin(), [](unsigned char c) {
        return static_cast<char>(std::tolower(c));
    });
    if (name == "none") return lsdtt::spectral_window_type::none;
    if (name == "hann") return lsdtt::spectral_window_type::hann;
    if (name == "hamming") return lsdtt::spectral_window_type::hamming;
    throw std::invalid_argument(
        "window must be one of 'none', 'hann', or 'hamming'"
    );
}

nb::object numpy_vector(const std::vector<double>& values)
{
    nb::object numpy = nb::module_::import_("numpy");
    return numpy.attr("array")(values, nb::arg("dtype") = "float64");
}

nb::object numpy_array2d(const lsdtt::lsdtt_array2d_f& values)
{
    nb::object numpy = nb::module_::import_("numpy");
    return numpy.attr("array")(values.get_vector(), nb::arg("dtype") = "float32")
        .attr("reshape")(nb::make_tuple(values.n_rows(), values.n_cols()));
}

class SpectralResult {
public:
    explicit SpectralResult(lsdtt::lsdtt_spectral_analysis_result result)
        : result_(std::move(result)) {}

    const lsdtt::lsdtt_spectral_analysis_result& native() const noexcept
    {
        return result_;
    }

private:
    lsdtt::lsdtt_spectral_analysis_result result_;
};

SpectralResult spectral_analysis(
    const Raster& raster,
    double log_bin_width,
    int n_background_iterations,
    int seed,
    const std::string& window,
    bool clip_to_valid_rectangle,
    int n_threads)
{
    lsdtt::lsdtt_spectral_analysis_options options;
    options.log_bin_width = log_bin_width;
    options.n_background_iterations = n_background_iterations;
    options.background_seed = seed;
    options.window = parse_window(window);
    options.clip_to_valid_rectangle = clip_to_valid_rectangle;
    options.n_threads = n_threads;
    return SpectralResult(lsdtt::full_spectral_analysis(raster.native(), options));
}

Raster lowpass_filter(const Raster& raster, double cutoff_metres,
                      const std::string& window, bool clip_to_valid_rectangle,
                      int n_threads)
{
    return Raster(lsdtt::fftw2D_lowpass_filter(
        raster.native(), cutoff_metres, parse_window(window),
        clip_to_valid_rectangle, n_threads
    ));
}

Raster highpass_filter(const Raster& raster, double cutoff_metres,
                       const std::string& window, bool clip_to_valid_rectangle,
                       int n_threads)
{
    return Raster(lsdtt::fftw2D_highpass_filter(
        raster.native(), cutoff_metres, parse_window(window),
        clip_to_valid_rectangle, n_threads
    ));
}

Raster bandpass_filter(const Raster& raster, double short_metres,
                       double long_metres, const std::string& window,
                       bool clip_to_valid_rectangle, int n_threads)
{
    return Raster(lsdtt::fftw2D_bandpass_filter(
        raster.native(), short_metres, long_metres, parse_window(window),
        clip_to_valid_rectangle, n_threads
    ));
}

Raster wiener_filter(const Raster& raster, double signal_wavelength_long_metres,
                      double signal_wavelength_short_metres, const std::string& window,
                     bool clip_to_valid_rectangle, int n_threads)
{
    return Raster(lsdtt::fftw2D_wiener_filter(
        raster.native(), signal_wavelength_long_metres, signal_wavelength_short_metres,
        parse_window(window), clip_to_valid_rectangle, n_threads
    ));
}

Raster masked_surface_recovery(const Raster& raster, const Raster& mask,
                               double cutoff_metres, const std::string& window,
                               int n_threads)
{
    return Raster(lsdtt::fftw2D_masked_surface_recovery(
        raster.native(), mask.native(), cutoff_metres, parse_window(window), n_threads
    ));
}

nb::dict detect_directional_periodicity(const Raster& raster, int window_size_pixels,
                                        int step_pixels, double min_spacing_metres,
                                        double max_spacing_metres, double min_valid_fraction,
                                        int n_threads)
{
    const auto result = lsdtt::detect_directional_periodicity(
        raster.native(), window_size_pixels, step_pixels, min_spacing_metres,
        max_spacing_metres, min_valid_fraction, n_threads
    );
    nb::dict output;
    output["strength"] = Raster(result.strength);
    output["spacing_metres"] = Raster(result.spacing_metres);
    output["orientation_degrees"] = Raster(result.orientation_deg);
    return output;
}

std::vector<double> frequency_axis(int size, double spacing)
{
    std::vector<double> axis(static_cast<std::size_t>(size));
    const int centre = size / 2;
    for (int i = 0; i < size; ++i)
        axis[static_cast<std::size_t>(i)] = (i - centre) * spacing;
    return axis;
}

} // namespace

void bind_raster_spectral(nb::module_& m)
{
    nb::class_<SpectralResult>(m, "SpectralResult", R"doc(
Results from a full DEM spectral analysis.

Spectral arrays are returned as copied NumPy arrays. Frequencies are in
cycles per metre and wavelengths are in metres.
)doc")
        .def_prop_ro("padded_shape", [](const SpectralResult& result) {
            return nb::make_tuple(result.native().padded_rows, result.native().padded_cols);
        })
        .def_prop_ro("frequency_x", [](const SpectralResult& result) {
            return numpy_vector(frequency_axis(result.native().padded_cols, result.native().dfx));
        }, nb::sig("def frequency_x(self) -> numpy.ndarray"))
        .def_prop_ro("frequency_y", [](const SpectralResult& result) {
            return numpy_vector(frequency_axis(result.native().padded_rows, result.native().dfy));
        }, nb::sig("def frequency_y(self) -> numpy.ndarray"))
        .def_prop_ro("psd_2d", [](const SpectralResult& result) {
            return numpy_array2d(result.native().psd_2d);
        }, nb::sig("def psd_2d(self) -> numpy.ndarray"))
        .def_prop_ro("background_psd_2d", [](const SpectralResult& result) {
            return numpy_array2d(result.native().background_psd_2d);
        }, nb::sig("def background_psd_2d(self) -> numpy.ndarray"))
        .def_prop_ro("frequency", [](const SpectralResult& result) { return numpy_vector(result.native().radial_frequency); }, nb::sig("def frequency(self) -> numpy.ndarray"))
        .def_prop_ro("wavelength", [](const SpectralResult& result) { return numpy_vector(result.native().radial_wavelength); }, nb::sig("def wavelength(self) -> numpy.ndarray"))
        .def_prop_ro("psd", [](const SpectralResult& result) { return numpy_vector(result.native().radial_psd); }, nb::sig("def psd(self) -> numpy.ndarray"))
        .def_prop_ro("background_psd", [](const SpectralResult& result) { return numpy_vector(result.native().background_psd); }, nb::sig("def background_psd(self) -> numpy.ndarray"))
        .def_prop_ro("normalised_psd", [](const SpectralResult& result) { return numpy_vector(result.native().normalised_psd); }, nb::sig("def normalised_psd(self) -> numpy.ndarray"))
        .def_prop_ro("rollover_frequency", [](const SpectralResult& result) { return result.native().rollover_frequency; })
        .def_prop_ro("rollover_wavelength", [](const SpectralResult& result) {
            const double frequency = result.native().rollover_frequency;
            return frequency > 0.0 ? 1.0 / frequency : 0.0;
        })
        .def_prop_ro("rollover_beta", [](const SpectralResult& result) { return result.native().rollover_beta; })
        .def_prop_ro("sub_rollover_beta", [](const SpectralResult& result) { return result.native().sub_rollover_beta; })
        .def_prop_ro("binned", [](const SpectralResult& result) {
            const auto& binned = result.native().binned;
            nb::dict output;
            output["frequency"] = numpy_vector(binned.frequency);
            output["wavelength"] = numpy_vector(binned.wavelength);
            output["psd"] = numpy_vector(binned.psd);
            output["background_psd"] = numpy_vector(binned.background_psd);
            output["normalised_psd"] = numpy_vector(binned.normalised_psd);
            output["beta"] = numpy_vector(binned.beta_by_bin);
            output["r_squared"] = numpy_vector(binned.r_sq_by_bin);
            return output;
        });

    m.def("spectral_available", []() {
#ifdef LSDTT_HAS_FFTW
        return true;
#else
        return false;
#endif
    }, R"doc(Return whether this build includes FFTW spectral support.

Returns
-------
bool
    True if the FFT-based functions (spectral analysis, filters, Wiener
    channel heads) are available.
)doc");
    m.def("spectral_analysis", &spectral_analysis, nb::arg("raster"),
        nb::arg("log_bin_width") = 0.05, nb::arg("n_background_iterations") = 100,
        nb::arg("seed") = -1, nb::arg("window") = "hann",
        nb::arg("clip_to_valid_rectangle") = false, nb::arg("n_threads") = 0,
        R"doc(Run full FFT-based spectral analysis and return a SpectralResult.

The detrended, tapered DEM is Fourier transformed and its power spectral
density is compared with a Monte Carlo background of red-noise surfaces
(Perron et al. 2008). Needs FFTW (``spectral_available()``).

Parameters
----------
raster : Raster
    Elevation raster (m) in a projected metre CRS.
log_bin_width : float, default=0.05
    Width of the log10-frequency bins of the binned spectrum.
n_background_iterations : int, default=100
    Number of background noise surfaces averaged.
seed : int, default=-1
    Base random seed of the background; negative uses a fixed built-in
    sequence (still reproducible).
window : {"hann", "hamming", "none"}, default="hann"
    Taper applied before the FFT.
clip_to_valid_rectangle : bool, default=False
    First crop to the largest all-valid rectangle (the result is then
    smaller than the input).
n_threads : int, default=0
    OpenMP thread count; 0 uses the runtime default.

Returns
-------
SpectralResult
)doc");
    m.def("lowpass_filter", &lowpass_filter, nb::arg("raster"), nb::arg("cutoff_metres"),
        nb::arg("window") = "hann", nb::arg("clip_to_valid_rectangle") = false,
        nb::arg("n_threads") = 0, R"doc(Return the long-wavelength DEM component.

Parameters
----------
raster : Raster
    Elevation raster (m) in a projected metre CRS.
cutoff_metres : float
    Wavelengths (m) longer than this are kept.
window : {"hann", "hamming", "none"}, default="hann"
    Taper applied before the FFT.
clip_to_valid_rectangle : bool, default=False
    First crop to the largest all-valid rectangle (the result is then
    smaller than the input).
n_threads : int, default=0
    OpenMP thread count; 0 uses the runtime default.

Returns
-------
Raster
)doc");
    m.def("highpass_filter", &highpass_filter, nb::arg("raster"), nb::arg("cutoff_metres"),
        nb::arg("window") = "hann", nb::arg("clip_to_valid_rectangle") = false,
        nb::arg("n_threads") = 0, R"doc(Return the short-wavelength DEM component.

Parameters
----------
raster : Raster
    Elevation raster (m) in a projected metre CRS.
cutoff_metres : float
    Wavelengths (m) shorter than this are kept.
window : {"hann", "hamming", "none"}, default="hann"
    Taper applied before the FFT.
clip_to_valid_rectangle : bool, default=False
    First crop to the largest all-valid rectangle (the result is then
    smaller than the input).
n_threads : int, default=0
    OpenMP thread count; 0 uses the runtime default.

Returns
-------
Raster
)doc");
    m.def("bandpass_filter", &bandpass_filter, nb::arg("raster"), nb::arg("short_metres"),
        nb::arg("long_metres"), nb::arg("window") = "hann",
        nb::arg("clip_to_valid_rectangle") = false, nb::arg("n_threads") = 0,
        R"doc(Return the DEM component between two wavelengths.

Parameters
----------
raster : Raster
    Elevation raster (m) in a projected metre CRS.
short_metres, long_metres : float
    Shortest and longest wavelength (m) kept; ``short_metres`` must be less
    than ``long_metres``.
window : {"hann", "hamming", "none"}, default="hann"
    Taper applied before the FFT.
clip_to_valid_rectangle : bool, default=False
    First crop to the largest all-valid rectangle (the result is then
    smaller than the input).
n_threads : int, default=0
    OpenMP thread count; 0 uses the runtime default.

Returns
-------
Raster
)doc");
    m.def("wiener_filter", &wiener_filter, nb::arg("raster"),
        nb::arg("signal_wavelength_long_metres") = 1000.0,
        nb::arg("signal_wavelength_short_metres") = 100.0, nb::arg("window") = "hann",
        nb::arg("clip_to_valid_rectangle") = false, nb::arg("n_threads") = 0,
        R"doc(Attenuate high-frequency white noise with an adaptive Wiener filter.

A power law fitted to the spectrum between the two signal wavelengths gives
the signal; the noise level comes from the near-Nyquist spectrum; each
frequency is weighted by signal / (signal + noise).

Parameters
----------
raster : Raster
    Elevation raster (m) in a projected metre CRS.
signal_wavelength_long_metres : float, default=1000.0
    Longest wavelength (m) of the band used to fit the signal.
signal_wavelength_short_metres : float, default=100.0
    Shortest wavelength (m) of that band.
window : {"hann", "hamming", "none"}, default="hann"
    Taper applied before the FFT.
clip_to_valid_rectangle : bool, default=False
    First crop to the largest all-valid rectangle (the result is then
    smaller than the input).
n_threads : int, default=0
    OpenMP thread count; 0 uses the runtime default.

Returns
-------
Raster
    Filtered DEM (m).
)doc");
    m.def("masked_surface_recovery", &masked_surface_recovery, nb::arg("raster"),
        nb::arg("mask"), nb::arg("cutoff_metres"), nb::arg("window") = "hann",
        nb::arg("n_threads") = 0, R"doc(Recover a smooth surface from selected pixels.

The DEM is detrended, zeroed outside ``mask`` and low-pass filtered, so a
smooth long-wavelength surface is interpolated under the masked-out area
(for example the piedmont surface beneath a mountain range).

Parameters
----------
raster : Raster
    Elevation raster (m).
mask : Raster
    0/1 raster; pixels > 0.5 are used.
cutoff_metres : float
    Low-pass cut-off wavelength (m).
window : {"hann", "hamming", "none"}, default="hann"
    Taper applied before the FFT.
n_threads : int, default=0
    OpenMP thread count; 0 uses the runtime default.

Returns
-------
Raster
    Smooth surface (m), valid wherever the DEM is valid.
)doc");
    m.def("detect_directional_periodicity", &detect_directional_periodicity,
        nb::arg("raster"), nb::arg("window_size_pixels"), nb::arg("step_pixels"),
        nb::arg("min_spacing_metres"), nb::arg("max_spacing_metres"),
        nb::arg("min_valid_fraction") = 0.6, nb::arg("n_threads") = 0,
        R"doc(Detect local regularly spaced parallel features with sliding-window FFTs.

Families of evenly spaced parallel lines (ditches, tracks, cultivation
ridges) give a sharp peak in the 2D power spectrum at 1 / spacing. A
window slides over the raster and the strongest peak in the spacing band
is recorded; each pixel keeps its best window.

Parameters
----------
raster : Raster
    Input raster, for example a DEM or a curvature response.
window_size_pixels : int
    Window side (pixels, >= 8); should span several of the largest spacings.
step_pixels : int
    Window stride (pixels, >= 1), typically half the window.
min_spacing_metres, max_spacing_metres : float
    Range of spacings (m) searched; the minimum must exceed 2 cell sizes and
    the maximum be below the window extent.
min_valid_fraction : float, default=0.6
    Windows with a smaller fraction of valid pixels are skipped.
n_threads : int, default=0
    OpenMP thread count; 0 uses the runtime default.

Returns
-------
dict
    ``"strength"`` (peak prominence), ``"spacing_metres"`` and
    ``"orientation_degrees"`` (line orientation) Rasters.
)doc");
}

} // namespace lsdtt3_python
