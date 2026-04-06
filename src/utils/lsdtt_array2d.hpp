// lsdtt_array2d.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html
//
// -----------------------------------------------------------------------------
// lsdtt_array2d<T>
//
// Data from rasters is stored using GDAL format in a continuous
// flattened vector. However code is  easier to read if you can use
// row,column referencing (at least Simon thinks it is easier).
// So this little helper object is made to facilitate that. 
// Any time heavy matrix compuation is required you should pass to 
// the Eigen linear algebra library, and if you want to do things 
// with image processing you pass to OpenCV.
//
// Or, more succinctly (in AI speak):
// Lightweight 2D array wrapper over a contiguous std::vector<T>.
// Provides row/column indexing via operator()(row, col) with the same
// memory layout as GDAL (row-major, origin top-left).
//
// All conversion helpers (OpenCV, Eigen, FFTW) produce zero-copy views
// over the same underlying buffer. No data is duplicated.
//
// Storage convention:
//   Row-major (C order). Element (row, col) is at index row * n_cols + col.
//   Row 0 is the northernmost row, consistent with GDAL's convention.
//
// Type parameter T:
//   float   — standard raster storage, matches GDAL GDT_Float32
//   double  — high-precision intermediate calculations
//   uint8_t — nodata masks, binary rasters
//   int32_t — integer classifications, flow direction grids
//
// Conversion helpers are conditionally compiled:
//   OpenCV  — include <opencv2/core.hpp> before this header
//   Eigen   — include <Eigen/Core> before this header
//   FFTW    — no special include needed; use data() directly
//
// Example usage:
//
//   lsdtt::lsdtt_array2d<float> dem(nrows, ncols, -9999.0f);
//
//   // Fill from GDAL
//   band->RasterIO(GF_Read, 0, 0, ncols, nrows, dem.data(), ...);
//
//   // Access by row and column
//   float z = dem(row, col);
//   dem(row, col) = 0.0f;
//
//   // If you want to check bounds use
//   float z = dem.at(row,col)
//
//   // Pass to OpenCV (zero copy)
//   cv::Mat cv_dem = dem.to_cv_mat();
//
//   // Pass to Eigen (zero copy, row-major)
//   auto eig = dem.to_eigen_map();
//
//   // Pass raw pointer to GDAL or FFTW
//   fftwf_plan p = fftwf_plan_dft_r2c_2d(nrows, ncols, dem.data(), ...);
// -----------------------------------------------------------------------------

#pragma once

#include <vector>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <algorithm>

// Conditionally include conversion targets if already available
#ifdef CV_VERSION
#include <opencv2/core.hpp>
#endif

#ifdef EIGEN_WORLD_VERSION
#include <Eigen/Core>
#endif

namespace lsdtt {

template<typename T>
class lsdtt_array2d
{
public:

    // -------------------------------------------------------------------------
    // Types
    // -------------------------------------------------------------------------

    using value_type     = T;
    using reference      = T&;
    using const_reference = const T&;
    using iterator       = typename std::vector<T>::iterator;
    using const_iterator = typename std::vector<T>::const_iterator;

    // -------------------------------------------------------------------------
    // Construction
    // -------------------------------------------------------------------------

    // Default — empty, invalid array
    lsdtt_array2d()
        : n_rows_(0), n_cols_(0)
    {}

    // Construct with given dimensions, all elements initialised to fill_value.
    // The fill_value parameter is typed as T so integer literals passed for
    // narrow types (uint8_t, int16_t etc.) must be explicitly cast:
    //   lsdtt_array2d<uint8_t> mask(nrows, ncols, uint8_t{0});   // correct
    //   lsdtt_array2d<uint8_t> mask(nrows, ncols, 0);            // ambiguous
    lsdtt_array2d(int n_rows, int n_cols, T fill_value = T{})
        : n_rows_(n_rows)
        , n_cols_(n_cols)
        , data_(static_cast<std::size_t>(n_rows) * n_cols, fill_value)
    {
        check_dimensions(n_rows, n_cols);
    }

    // Tag type used to disambiguate the pointer constructor from the fill
    // constructor. Without this tag, passing an integer literal 0 as the
    // third argument is ambiguous — the compiler cannot decide whether it
    // is a fill value or a null pointer.
    //
    // Usage:
    //   lsdtt_array2d<float> dem(nrows, ncols, lsdtt_array2d<float>::from_ptr, ptr);
    struct from_ptr_tag {};
    static constexpr from_ptr_tag from_ptr{};

    // Construct from an existing contiguous buffer (copies data).
    // You give it the pointer and the the array2d then has the data, but need
    // the from_ptr_tag to make sure it doens't think 0 is a null pointer
    // Use the from_ptr tag to call this constructor unambiguously:
    //   lsdtt_array2d<float> dem(nrows, ncols, lsdtt_array2d<float>::from_ptr, src);
    lsdtt_array2d(int n_rows, int n_cols, from_ptr_tag, const T* src)
        : n_rows_(n_rows)
        , n_cols_(n_cols)
        , data_(src, src + static_cast<std::size_t>(n_rows) * n_cols)
    {
        check_dimensions(n_rows, n_cols);
    }

    // Construct by moving an existing vector (zero copy)
    // Vector size must equal n_rows * n_cols
    lsdtt_array2d(int n_rows, int n_cols, std::vector<T>&& vec)
        : n_rows_(n_rows)
        , n_cols_(n_cols)
        , data_(std::move(vec))
    {
        check_dimensions(n_rows, n_cols);
        if (static_cast<std::size_t>(n_rows) * n_cols != data_.size())
        {
            throw std::runtime_error(
                "lsdtt_array2d: vector size " + std::to_string(data_.size()) +
                " does not match dimensions " + std::to_string(n_rows) +
                " x " + std::to_string(n_cols)
            );
        }
    }

    // Copy and move — default implementations are correct
    lsdtt_array2d(const lsdtt_array2d&)            = default;
    lsdtt_array2d& operator=(const lsdtt_array2d&) = default;
    lsdtt_array2d(lsdtt_array2d&&)                 = default;
    lsdtt_array2d& operator=(lsdtt_array2d&&)      = default;

    // -------------------------------------------------------------------------
    // Element access
    // -------------------------------------------------------------------------

    // Bounds-checked access
    // Slower but safer
    reference at(int row, int col)
    {
        check_bounds(row, col);
        return data_[index(row, col)];
    }

    const_reference at(int row, int col) const
    {
        check_bounds(row, col);
        return data_[index(row, col)];
    }

    // Unchecked access — use in performance-critical inner loops
    // Because stuff is accessed via a pointer anything beyond bounds
    // will still result in a number, but it won't be in the raster
    // so make sure you use this where you are confident in the size of the
    // raster
    reference operator()(int row, int col)
    {
        return data_[index(row, col)];
    }

    const_reference operator()(int row, int col) const
    {
        return data_[index(row, col)];
    }

    // -------------------------------------------------------------------------
    // Dimensions
    // -------------------------------------------------------------------------

    int  n_rows()   const { return n_rows_; }
    int  n_cols()   const { return n_cols_; }
    bool is_valid() const { return n_rows_ > 0 && n_cols_ > 0; }

    std::size_t size() const { return data_.size(); }

    // Returns true if (row, col) is within the array bounds
    bool in_bounds(int row, int col) const
    {
        return row >= 0 && row < n_rows_ && col >= 0 && col < n_cols_;
    }

    // Returns true if (row, col) is an interior pixel
    // (i.e. a full 3x3 neighbourhood exists around it)
    bool is_interior(int row, int col) const
    {
        return row > 0 && row < n_rows_ - 1 &&
               col > 0 && col < n_cols_ - 1;
    }

    // -------------------------------------------------------------------------
    // Raw memory access
    // Used by GDAL RasterIO, FFTW, and the zero-copy conversion helpers below
    // -------------------------------------------------------------------------

    T*       data()       { return data_.data(); }
    const T* data() const { return data_.data(); }

    // Access the underlying vector directly (e.g. for std::move)
    std::vector<T>&       get_vector()       { return data_; }
    const std::vector<T>& get_vector() const { return data_; }

    // -------------------------------------------------------------------------
    // Iterators — for range-based for loops and STL algorithms
    // -------------------------------------------------------------------------

    iterator       begin()        { return data_.begin(); }
    iterator       end()          { return data_.end();   }
    const_iterator begin()  const { return data_.begin(); }
    const_iterator end()    const { return data_.end();   }
    const_iterator cbegin() const { return data_.cbegin(); }
    const_iterator cend()   const { return data_.cend();   }

    // -------------------------------------------------------------------------
    // Fill and resize
    // -------------------------------------------------------------------------

    void fill(T value)
    {
        std::fill(data_.begin(), data_.end(), value);
    }

    // Resize and reinitialise — existing data is discarded
    void resize(int n_rows, int n_cols, T fill_value = T{})
    {
        check_dimensions(n_rows, n_cols);
        n_rows_ = n_rows;
        n_cols_ = n_cols;
        data_.assign(static_cast<std::size_t>(n_rows) * n_cols, fill_value);
    }

    // -------------------------------------------------------------------------
    // Type conversion — create a new array with a different element type
    // Always produces a copy (explicit, intentional type promotion/demotion)
    // -------------------------------------------------------------------------

    template<typename U>
    lsdtt_array2d<U> as_type() const
    {
        lsdtt_array2d<U> result(n_rows_, n_cols_);
        std::transform(
            data_.begin(), data_.end(),
            result.get_vector().begin(),
            [](const T& v) { return static_cast<U>(v); }
        );
        return result;
    }

    // -------------------------------------------------------------------------
    // OpenCV zero-copy view
    // Requires <opencv2/core.hpp> to be included before this header.
    // The returned cv::Mat shares memory with this array — do not let the
    // lsdtt_array2d go out of scope while the cv::Mat is in use.
    // -------------------------------------------------------------------------

#ifdef CV_VERSION
    cv::Mat to_cv_mat()
    {
        return cv::Mat(n_rows_, n_cols_, cv_type(), data_.data());
    }

    cv::Mat to_cv_mat() const
    {
        // const_cast is safe here — cv::Mat is used read-only
        return cv::Mat(n_rows_, n_cols_, cv_type(),
                       const_cast<T*>(data_.data()));
    }

private:
    // Map T to the corresponding OpenCV type code
    static int cv_type()
    {
        if constexpr (std::is_same_v<T, float>)   return CV_32F;
        if constexpr (std::is_same_v<T, double>)  return CV_64F;
        if constexpr (std::is_same_v<T, uint8_t>) return CV_8U;
        if constexpr (std::is_same_v<T, int32_t>) return CV_32S;
        if constexpr (std::is_same_v<T, int16_t>) return CV_16S;
        static_assert(
            std::is_same_v<T, float> || std::is_same_v<T, double> ||
            std::is_same_v<T, uint8_t> || std::is_same_v<T, int32_t> ||
            std::is_same_v<T, int16_t>,
            "lsdtt_array2d::to_cv_mat() — unsupported type T"
        );
        return CV_32F;  // unreachable but suppresses warnings
    }
public:
#endif

    // -------------------------------------------------------------------------
    // Eigen zero-copy view (row-major)
    // Requires <Eigen/Core> to be included before this header.
    //
    // IMPORTANT: Eigen defaults to column-major storage. This map explicitly
    // uses row-major to match GDAL and OpenCV conventions. Always use the
    // returned map type rather than a plain Eigen::MatrixXf to avoid
    // accidentally re-interpreting the memory as column-major.
    //
    // The returned map shares memory with this array — do not let the
    // lsdtt_array2d go out of scope while the map is in use.
    // -------------------------------------------------------------------------

#ifdef EIGEN_WORLD_VERSION
    using eigen_map_t = Eigen::Map<
        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>
    >;
    using const_eigen_map_t = Eigen::Map<
        const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>
    >;

    eigen_map_t to_eigen_map()
    {
        return eigen_map_t(data_.data(), n_rows_, n_cols_);
    }

    const_eigen_map_t to_eigen_map() const
    {
        return const_eigen_map_t(data_.data(), n_rows_, n_cols_);
    }
#endif

private:

    // -------------------------------------------------------------------------
    // Private helpers
    // -------------------------------------------------------------------------

    int n_rows_;
    int n_cols_;
    std::vector<T> data_;

    std::size_t index(int row, int col) const
    {
        return static_cast<std::size_t>(row) * n_cols_ + col;
    }

    static void check_dimensions(int n_rows, int n_cols)
    {
        if (n_rows < 0 || n_cols < 0)
        {
            throw std::runtime_error(
                "lsdtt_array2d: negative dimensions (" +
                std::to_string(n_rows) + ", " + std::to_string(n_cols) + ")"
            );
        }
    }

    void check_bounds(int row, int col) const
    {
        if (row < 0 || row >= n_rows_ || col < 0 || col >= n_cols_)
        {
            throw std::runtime_error(
                "lsdtt_array2d: index (" + std::to_string(row) + ", " +
                std::to_string(col) + ") out of bounds for array of size (" +
                std::to_string(n_rows_) + ", " + std::to_string(n_cols_) + ")"
            );
        }
    }
};

// -----------------------------------------------------------------------------
// Convenience type aliases
// -----------------------------------------------------------------------------

using lsdtt_array2d_f  = lsdtt_array2d<float>;    // standard raster
using lsdtt_array2d_d  = lsdtt_array2d<double>;   // high precision
using lsdtt_array2d_i  = lsdtt_array2d<int32_t>;  // integer grids
using lsdtt_array2d_u8 = lsdtt_array2d<uint8_t>;  // masks, binary rasters

} // namespace lsdtt
