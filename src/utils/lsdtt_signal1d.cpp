// lsdtt_signal1d.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Rutger Vervoordeldonk and Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html

#include "lsdtt_signal1d.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace lsdtt
{

// -----------------------------------------------------------------------------
// Smoothing Algorithms
// -----------------------------------------------------------------------------

std::vector<double> compute_moving_average_smoothing(
    const std::vector<double>& y,
    int window,
    bool periodic)
{
    const int n = static_cast<int>(y.size());
    if (n == 0 || window < 1) return y;

    // Enforce odd window length for natural symmetry symmetry
    if (window % 2 == 0) ++window;
    const int half = window / 2;

    std::vector<double> out(n, 0.0);

    for (int i = 0; i < n; ++i) 
    {
        double sum = 0.0;
        int count = 0;

        for (int k = -half; k <= half; ++k) 
        {
            int j = i + k;
            if (periodic) 
            {
                j = (j % n + n) % n;
                sum += y[j];
                ++count;
            } 
            else 
            {
                if (j >= 0 && j < n) 
                {
                    sum += y[j];
                    ++count;
                }
            }
        }
        if (count > 0) out[i] = sum / static_cast<double>(count);
    }
    return out;
}

std::vector<double> compute_lowess_smoothing(
    const std::vector<double>& x,
    const std::vector<double>& y,
    double bandwidth,
    int iterations)
{
    const std::size_t n = x.size();
    if (n != y.size()) throw std::invalid_argument("LOWESS: x and y must be same length.");
    if (n == 0) return {};

    std::vector<double> smoothed_y = y; 
    std::vector<double> weights(n);
    std::vector<double> residuals(n, 0.0);

    auto get_median = [](std::vector<double> v) -> double {
        if (v.empty()) return 0.0;
        std::size_t mid_n = v.size() / 2;
        std::nth_element(v.begin(), v.begin() + mid_n, v.end());
        if (v.size() % 2 == 0) {
            auto rit = std::max_element(v.begin(), v.begin() + mid_n);
            return (*rit + v[mid_n]) / 2.0;
        }
        return v[mid_n];
    };

    auto weighted_least_squares = [](const std::vector<double>& dx, const std::vector<double>& dy, const std::vector<double>& dw, double& intercept, double& gradient) {
        double sw = 0, swx = 0, swx2 = 0, swy = 0, swxy = 0;
        for (std::size_t j = 0; j < dx.size(); ++j) {
            double wx = dw[j] * dx[j];
            sw   += dw[j];
            swx  += wx;
            swx2 += wx * dx[j];
            swy  += dw[j] * dy[j];
            swxy += wx * dy[j];
        }
        double D = sw * swx2 - swx * swx;
        if (std::abs(D) > 1e-15) {
            intercept = (swx2 * swy - swx * swxy) / D;
            gradient  = (sw * swxy - swx * swy) / D;
        } else {
            intercept = sw != 0.0 ? (swy / sw) : 0.0;
            gradient  = 0.0;
        }
    };

    for (int iter = 0; iter < iterations; ++iter) 
    {
        double scale = 0.0;
        if (iter > 0) {
            // Note: LSDStatsTools updated residuals per point natively during smoothing,
            // resulting in a slight dependency order anomaly. We calculate it correctly over the whole vector block.
            scale = 6.0 * get_median(residuals); 
        }

        std::vector<double> current_smoothed_y(n, 0.0);

        for (std::size_t i = 0; i < n; ++i) 
        {
            double x_i = x[i];

            for (std::size_t j = 0; j < n; ++j) 
            {
                double w = 0.0;
                double dist = std::abs(x_i - x[j]) / bandwidth;
                
                // Tricube distance kernel
                if (dist < 1.0) {
                    double tmp = 1.0 - (dist * dist * dist);
                    w = tmp * tmp * tmp;
                }
                
                // Bisquare robustness weighting
                if (iter > 0 && scale != 0.0) {
                    double r = residuals[j] / scale;
                    if (std::abs(r) < 1.0) {
                        double tmp = 1.0 - r * r;
                        w *= (tmp * tmp);
                    } else {
                        w = 0.0;
                    }
                }
                weights[j] = w;
            }

            double intercept, gradient;
            weighted_least_squares(x, y, weights, intercept, gradient);
            current_smoothed_y[i] = intercept + gradient * x_i;
        }

        smoothed_y = current_smoothed_y;

        // Recalculate robust residuals blockwise
        for (std::size_t i = 0; i < n; ++i) {
            residuals[i] = std::abs(y[i] - smoothed_y[i]);
        }
    }

    return smoothed_y;
}

// -----------------------------------------------------------------------------
// Profile Signal Operators
// -----------------------------------------------------------------------------

std::vector<double> compute_central_difference(
    const std::vector<double>& x,
    const std::vector<double>& y,
    bool periodic)
{
    const std::size_t n = x.size();
    if (n != y.size()) throw std::invalid_argument("Central diff: x and y must be same length.");
    std::vector<double> dydx(n, 0.0);
    if (n < 2) return dydx;

    for (std::size_t i = 0; i < n; ++i) 
    {
        std::size_t i_prev, i_next;

        if (i == 0) {
            i_prev = periodic ? n - 1 : i;
            i_next = i + 1;
        } else if (i == n - 1) {
            i_prev = i - 1;
            i_next = periodic ? 0 : i;
        } else {
            i_prev = i - 1;
            i_next = i + 1;
        }

        // We explicitly track forwards/backwards relative spacing logic
        // This is robust against coordinate wrapping in cumulative distanced profiles (lsdtt2)
        double dx_forward  = (i_next > i) ? (x[i_next] - x[i]) : (x[i] - x[i - 1]);
        double dx_backward = (i > i_prev) ? (x[i] - x[i_prev]) : (x[i + 1] - x[i]);
        double dx = dx_forward + dx_backward;

        if (std::abs(dx) > 1e-15) {
            dydx[i] = (y[i_next] - y[i_prev]) / dx;
        }
    }
    return dydx;
}

std::vector<std::size_t> find_knees(
    const std::vector<double>& x,
    const std::vector<double>& y,
    double sensitivity)
{
    const std::size_t n = x.size();
    if (n != y.size()) throw std::invalid_argument("Kneedle: x and y must be same length.");
    if (n < 3) return {};

    double min_x = *std::min_element(x.begin(), x.end());
    double max_x = *std::max_element(x.begin(), x.end());
    double min_y = *std::min_element(y.begin(), y.end());
    double max_y = *std::max_element(y.begin(), y.end());

    double range_x = (max_x > min_x) ? max_x - min_x : 1.0;
    double range_y = (max_y > min_y) ? max_y - min_y : 1.0;

    std::vector<double> x_norm(n), y_norm(n), y_diff(n);
    for (std::size_t i = 0; i < n; ++i) 
    {
        x_norm[i] = (x[i] - min_x) / range_x;
        y_norm[i] = (y[i] - min_y) / range_y;
        y_diff[i] = y_norm[i] - x_norm[i];
    }

    std::vector<std::size_t> local_maxima;
    for (std::size_t i = 1; i < n - 1; ++i) 
    {
        if (y_diff[i - 1] < y_diff[i] && y_diff[i + 1] < y_diff[i]) {
            local_maxima.push_back(i);
        }
    }

    if (local_maxima.empty()) return {};

    double avg_dx_term = sensitivity * (x_norm.back() - x_norm.front()) / static_cast<double>(n > 1 ? n - 1 : 1);

    std::vector<std::size_t> knees;
    for (std::size_t lm_idx = 0; lm_idx < local_maxima.size(); ++lm_idx) 
    {
        std::size_t lmx_i = local_maxima[lm_idx];
        double Tlmxi = y_diff[lmx_i] - avg_dx_term;

        for (std::size_t j = lmx_i + 1; j < n; ++j) 
        {
            if (lm_idx + 1 < local_maxima.size() && j == local_maxima[lm_idx + 1]) break;
            if (j + 1 < n && y_diff[j] < y_diff[j + 1]) break; // Trailing local minimum guard

            if (y_diff[j] < Tlmxi) {
                knees.push_back(lmx_i);
                break;
            }
        }
    }
    return knees;
}

} // namespace lsdtt
