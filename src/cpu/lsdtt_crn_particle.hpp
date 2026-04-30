#pragma once
#include "lsdtt_crn_parameters.hpp"
#include <cmath>
#include <iostream>

namespace lsdtt
{

// ─────────────────────────────────────────────────────────────────────────────
// lsdtt_crn_particle
//
// Lightweight particle that tracks cosmogenic nuclide concentrations.
// Header-only; all methods are inline because the formulas are compact.
//
// Before calling any update_*_ss method, the caller must have already:
//   1. Called params.reset_scaling(muon_scheme)
//   2. Called params.scale_F_values(total_shielding, nuclide)
//
// Steady-state depth-integrated formula (all pathways i = 0..3):
//
//   Surface (d_top == d_bottom == d):
//     C = P0 · Σᵢ F[i]·Γ[i] · exp(-d/Γ[i]) / (ε + Γ[i]·λ)
//
//   Depth-integrated [d_top, d_bottom]:
//     C = P0/(d_bot-d_top) · Σᵢ F[i]·Γ[i]² / (ε + Γ[i]·λ) ·
//              [exp(-d_top/Γ[i]) - exp(-d_bot/Γ[i])]
//
// F[i] values carry the shielding (via scale_F_values), so no S_t multiplier
// is needed here (S_t = 1 after reset_scaling).
//
// Units: concentration in atoms g⁻¹, erosion_rate in g cm⁻² yr⁻¹.
// ─────────────────────────────────────────────────────────────────────────────
class lsdtt_crn_particle
{
public:
    lsdtt_crn_particle() : conc_10Be_(0.0), conc_26Al_(0.0), conc_14C_(0.0) {}

    // ── Steady-state concentration updates ───────────────────────────────────
    // erosion_rate : effective erosion rate (g cm⁻² yr⁻¹, > 0)
    // params       : CRN parameters with F values already scaled by shielding
    // d_eff_top    : effective depth of sample top (g cm⁻²), 0 for surface
    // d_eff_bottom : effective depth of sample bottom; equal to d_eff_top for
    //               a point sample, > d_eff_top for a depth-integrated sample

    void update_10Be_ss(double erosion_rate,
                        const lsdtt_crn_parameters& params,
                        double d_eff_top    = 0.0,
                        double d_eff_bottom = 0.0)
    {
        conc_10Be_ = compute_ss(erosion_rate,
                                params.get_P0_10Be(), params.get_lambda_10Be(),
                                params, &lsdtt_crn_parameters::get_F_10Be,
                                d_eff_top, d_eff_bottom);
    }

    void update_26Al_ss(double erosion_rate,
                        const lsdtt_crn_parameters& params,
                        double d_eff_top    = 0.0,
                        double d_eff_bottom = 0.0)
    {
        conc_26Al_ = compute_ss(erosion_rate,
                                params.get_P0_26Al(), params.get_lambda_26Al(),
                                params, &lsdtt_crn_parameters::get_F_26Al,
                                d_eff_top, d_eff_bottom);
    }

    void update_14C_ss(double erosion_rate,
                       const lsdtt_crn_parameters& params,
                       double d_eff_top    = 0.0,
                       double d_eff_bottom = 0.0)
    {
        conc_14C_ = compute_ss(erosion_rate,
                               params.get_P0_14C(), params.get_lambda_14C(),
                               params, &lsdtt_crn_parameters::get_F_14C,
                               d_eff_top, d_eff_bottom);
    }

    // ── Getters ──────────────────────────────────────────────────────────────
    double get_conc_10Be() const { return conc_10Be_; }
    double get_conc_26Al() const { return conc_26Al_; }
    double get_conc_14C()  const { return conc_14C_;  }

private:
    double conc_10Be_, conc_26Al_, conc_14C_;

    // Generic steady-state helper — works for any nuclide via a pointer-to-member
    // getter for the F[] array.
    using FGetter = double (lsdtt_crn_parameters::*)(int) const;

    static double compute_ss(double erosion_rate,
                             double P0, double lambda,
                             const lsdtt_crn_parameters& params,
                             FGetter get_F,
                             double d_top, double d_bot)
    {
        if (d_top > d_bot)
        {
            std::cout << "lsdtt_crn_particle: depth range reversed; swapping." << std::endl;
            std::swap(d_top, d_bot);
        }

        double sum = 0.0;

        if (d_top == d_bot)
        {
            // Point sample at effective depth d_top
            for (int i = 0; i < 4; ++i)
            {
                double F = (params.*get_F)(i);
                double G = params.get_Gamma(i);
                sum += F * G * std::exp(-d_top / G) / (erosion_rate + G * lambda);
            }
        }
        else
        {
            // Depth-integrated over [d_top, d_bot]
            double denom = d_bot - d_top;
            for (int i = 0; i < 4; ++i)
            {
                double F = (params.*get_F)(i);
                double G = params.get_Gamma(i);
                double top_term = std::exp(-d_top / G);
                double bot_term = std::exp(-d_bot / G);
                sum += F * G * G * (top_term - bot_term) / ((erosion_rate + G * lambda) * denom);
            }
        }

        return P0 * sum;
    }
};

} // namespace lsdtt
