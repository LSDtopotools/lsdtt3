#pragma once
// =============================================================================
//  DO NOT CHANGE THE NUMBERS THIS CODE PRODUCES — ask Simon Mudd first.
//
//  The cosmogenic calculations in this file (the CAIRN erosion-rate method and
//  the step-change concentration solutions) were checked by hand over several
//  years. To stop anyone changing those numbers by accident, the test suite
//  keeps a stored list of known-correct output values in
//  tests/expected/cosmo/numerics_sentinel.csv. The tests in
//  tests/cpp/test_cosmo_physics.cpp re-run these calculations and compare the
//  results against that stored list to 12 decimal places.
//
//  So if you edit this file and any output value changes — even from a refactor
//  or an algebra rewrite you believe is harmless — those tests will fail. When
//  that happens, do NOT edit the stored list to make the tests pass again. A
//  changed value has to be reviewed and approved by Simon Mudd before it is
//  committed.
// =============================================================================
#include "lsdtt_crn_parameters.hpp"
#include <cmath>
#include <iostream>

namespace lsdtt
{

// 
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
//     C = P0 * Sigma_i F[i]*Gamma[i] * exp(-d/Gamma[i]) / (epsilon + Gamma[i]*lambda)
//
//   Depth-integrated [d_top, d_bottom]:
//     C = P0/(d_bot-d_top) * Sigma_i F[i]*Gamma[i]^2 / (epsilon + Gamma[i]*lambda) *
//              [exp(-d_top/Gamma[i]) - exp(-d_bot/Gamma[i])]
//
// F[i] values carry the shielding (via scale_F_values), so no S_t multiplier
// is needed here (S_t = 1 after reset_scaling).
//
// Units: concentration in atoms g^-^1, erosion_rate in g cm^-^2 yr^-^1.
// 
class lsdtt_crn_particle
{
public:
    lsdtt_crn_particle()
        : conc_10Be_(0.0), conc_26Al_(0.0), conc_14C_(0.0), d_eff_(0.0) {}

    explicit lsdtt_crn_particle(double d_eff_g_cm2)
        : conc_10Be_(0.0), conc_26Al_(0.0), conc_14C_(0.0), d_eff_(d_eff_g_cm2) {}

    //  Steady-state concentration updates 
    // erosion_rate : effective erosion rate (g cm^-^2 yr^-^1, > 0)
    // params       : CRN parameters with F values already scaled by shielding
    // d_eff_top    : effective depth of sample top (g cm^-^2), 0 for surface
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

    //  SS initialiser using tracked depth 
    // Computes steady-state for all three nuclides at the particle's current d_eff_.
    void init_ss_all(double erosion_rate, const lsdtt_crn_parameters& params)
    {
        conc_10Be_ = compute_ss(erosion_rate, params.get_P0_10Be(), params.get_lambda_10Be(),
                                params, &lsdtt_crn_parameters::get_F_10Be, d_eff_, d_eff_);
        conc_26Al_ = compute_ss(erosion_rate, params.get_P0_26Al(), params.get_lambda_26Al(),
                                params, &lsdtt_crn_parameters::get_F_26Al, d_eff_, d_eff_);
        conc_14C_  = compute_ss(erosion_rate, params.get_P0_14C(),  params.get_lambda_14C(),
                                params, &lsdtt_crn_parameters::get_F_14C,  d_eff_, d_eff_);
    }

    //  Transient update 
    // Advances concentrations by dt (yr) at constant erosion rate eff_erate (g cm^-^2 yr^-^1).
    // Uses the analytic solution:
    //   C_new = C_old * exp(-lambda*dt)
    //         + P0 * Sigma_i F_i * exp(-d/Gamma_i) * [exp(epsilon*dt/Gamma_i) - exp(-lambda*dt)] / (epsilon/Gamma_i + lambda)
    // where d = d_eff_ (the depth at the START of the time step).
    // Depth is NOT updated here - call remove_mass() after.
    void update_conc(double dt, double eff_erate, const lsdtt_crn_parameters& params)
    {
        double d0 = d_eff_;   // capture starting depth before any depth change
        conc_10Be_ = advance_conc(conc_10Be_, dt, eff_erate, d0,
                                  params.get_P0_10Be(), params.get_lambda_10Be(),
                                  params, &lsdtt_crn_parameters::get_F_10Be);
        conc_26Al_ = advance_conc(conc_26Al_, dt, eff_erate, d0,
                                  params.get_P0_26Al(), params.get_lambda_26Al(),
                                  params, &lsdtt_crn_parameters::get_F_26Al);
        conc_14C_  = advance_conc(conc_14C_,  dt, eff_erate, d0,
                                  params.get_P0_14C(),  params.get_lambda_14C(),
                                  params, &lsdtt_crn_parameters::get_F_14C);
    }

    //  Depth management 
    // Remove depth_m metres of material from the surface (particle moves up).
    // density_kg_m3: rock density; converts metres to effective depth (g/cm^2).
    void remove_mass(double depth_m, double density_kg_m3)
    {
        d_eff_ -= depth_m * density_kg_m3 * 0.1;
        if (d_eff_ < 0.0) d_eff_ = 0.0;
    }

    void set_effective_depth(double d) { d_eff_ = d; }
    void set_conc_10Be(double c) { conc_10Be_ = c; }
    void set_conc_26Al(double c) { conc_26Al_ = c; }
    void set_conc_14C(double c)  { conc_14C_  = c; }

    //  Getters 
    double get_conc_10Be()       const { return conc_10Be_; }
    double get_conc_26Al()       const { return conc_26Al_; }
    double get_conc_14C()        const { return conc_14C_;  }
    double get_effective_depth() const { return d_eff_; }

private:
    double conc_10Be_, conc_26Al_, conc_14C_;
    double d_eff_;   // current effective depth (g cm^-^2)

    using FGetter = double (lsdtt_crn_parameters::*)(int) const;

    // Advance one nuclide concentration by dt yr at erosion rate eff_erate (g/cm^2/yr).
    // d0: effective depth at the START of the time step (g/cm^2).
    // Analytic solution of dC/dt = P0*Sigma_i F_i*exp(-d(t)/Gamma_i) - lambda*C  with d(t)=d0-epsilon*t:
    //   C(dt) = C_old*exp(-lambda*dt)
    //         + P0*Sigma_i F_i*exp(-d0/Gamma_i)*[exp(epsilon*dt/Gamma_i) - exp(-lambda*dt)] / (epsilon/Gamma_i+lambda)
    // (Mudd et al. 2016 ESurf Eq. 4; knickpoint-paper Eq. 7.)
    // The growth term is grouped as exp(-(d0-epsilon*dt)/Gamma_i), the
    // end-of-step depth (>= 0 for all callers), so it cannot overflow for
    // large epsilon*dt/Gamma_i the way exp(-d0/G)*exp(epsilon*dt/G) does.
    // Fix approved by Simon Mudd 2026-07-16: the previous implementation
    // multiplied the growth term by an extra exp(-lambda*dt) (missing
    // +lambda*dt in the exponent), biasing all transients low by ~lambda*t.
    static double advance_conc(double C_old, double dt, double eff_erate, double d0,
                               double P0, double lambda,
                               const lsdtt_crn_parameters& params, FGetter get_F)
    {
        double exp_decay = std::exp(-lambda * dt);
        double sum = 0.0;
        for (int i = 0; i < 4; ++i)
        {
            double G     = params.get_Gamma(i);
            double alpha = eff_erate / G;    // epsilon/Gamma_i  (yr^-^1, positive)
            double denom = alpha + lambda;
            if (denom <= 0.0) continue;
            // exp(-d0/G) * [exp(alpha*dt) - exp(-lambda*dt)], stable grouping
            double term = std::exp(-(d0 - eff_erate * dt) / G)
                          - std::exp(-d0 / G) * exp_decay;
            sum += (params.*get_F)(i) * term / denom;
        }
        return C_old * exp_decay + P0 * sum;
    }

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
