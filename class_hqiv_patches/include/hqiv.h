/** @file hqiv.h Documented includes for HQIV module
 *
 * Horizon-Quantized Informational Vacuum (HQIV) Cosmological Framework
 * 
 * This module implements the HQIV modifications to standard ΛCDM cosmology
 * as derived in the paper "Horizon-Quantized Informational Vacuum: 
 * A Covariant Baryon-Only Cosmological Framework from Quantised Inertia"
 * 
 * Key equations:
 * - Modified Friedmann (Paper Section 5, Eq. 13): 3H² - γH = 8πG_eff(ρ_m + ρ_r)
 * - Inertia modification (Paper Section 4.2, Eq. 11): f(α,ϕ) = max(α/(α + χϕ/6), f_min)
 * - Varying G (Paper Section 2.3): G(a) = G0 × (H(a)/H0)^α
 * - Vorticity source (Paper Section 6): ∂ω/∂t + (v·∇)ω = (∂f/∂ϕ)(k × ∇ϕ)·ê_ω
 *
 * Author: HQIV Team
 * Reference: new-final.tex (February 2026)
 */

#ifndef __HQIV__
#define __HQIV__

#include "common.h"

/**
 * HQIV parameter structure
 * 
 * All parameters have theoretical origins as described in the paper:
 * - γ: fixed by Brodie's Rindler-cosmic horizon overlap integral
 * - α: derived from galaxy rotation curve requirements
 * - χ: from full light-cone average
 * - f_min: saturation of informational-energy budget
 */

struct hqiv_parameters {
  
  /** @name - Main HQIV switch */
  
  short hqiv_on;           /**< If _TRUE_, use HQIV modified expansion. If _FALSE_, standard ΛCDM */
  
  /** @name - Background parameters (Paper Section 5, Table 1) */
  
  double gamma_hqiv;       /**< Thermodynamic coefficient from Brodie's overlap integral 
                                Default: 0.40, Range: 0.35-0.45 */
  
  double alpha_hqiv;       /**< Exponent for varying G(a), derived from rotation curves
                                Default: 0.60 */
  
  double chi_hqiv;         /**< Scaling factor from light-cone average for minimum acceleration
                                Default: 0.172 */
  
  double fmin_hqiv;        /**< Thermodynamic floor from informational-energy saturation
                                Default: 0.01 */
  
  /** @name - Derived quantities computed during evolution */
  
  double phi_horizon;      /**< Horizon field φ = cH (in 1/s units, computed at each step) */
  
  double G_eff_ratio;      /**< G_eff/G0 ratio at current scale factor */
  
  double inertia_factor;   /**< Inertia reduction factor f(α,φ) at current scale factor */
  
  /** @name - Output quantities */
  
  double age_hqiv;         /**< Universe age in Gyr (CLASS-consistent; ~32 Gyr for fiducial) */
  
  double vorticity_growth_exponent;  /**< Vorticity growth exponent (target: +1.9) */
  
  double raw_acoustic_scale;         /**< Raw acoustic scale ℓ_A (target: ~340) */
  
};

/**
 * Indices for HQIV-specific background quantities
 */

enum hqiv_background_indices {
  index_hqiv_phi,          /**< Horizon field φ = cH */
  index_hqiv_G_eff,        /**< Effective gravitational coupling G_eff/G0 */
  index_hqiv_f_inertia,    /**< Inertia reduction factor f(α,φ) */
  index_hqiv_df_dphi,      /**< Derivative df/dφ for vorticity source */
  hqiv_bg_size             /**< Total size of HQIV background vector */
};

/**************************************************************/
/* @cond INCLUDE_WITH_DOXYGEN */
/*
 * Boilerplate for C++
 */
#ifdef __cplusplus
extern "C" {
#endif

  /**
   * Initialize HQIV parameters with default values
   * 
   * @param phqiv  Output: pointer to HQIV parameters structure
   * @return the error status
   */
  
  int hqiv_init(struct hqiv_parameters * phqiv);
  
  /**
   * Free HQIV parameters
   * 
   * @param phqiv  Input: pointer to HQIV parameters structure
   * @return the error status
   */
  
  int hqiv_free(struct hqiv_parameters * phqiv);
  
  /**
   * Compute HQIV-modified Hubble parameter
   * 
   * Solves the quadratic Friedmann equation (Paper Eq. 13, Appendix A.2):
   *   3H² - γH = 8πG_eff(ρ_m + ρ_r)
   * 
   * Positive root solution:
   *   H = [γ + √(γ² + 96πG_eff ρ_tot)] / 6
   * 
   * @param gamma        Input: thermodynamic coefficient
   * @param rho_tot      Input: total energy density (matter + radiation)
   * @param G_eff_ratio  Input: G_eff/G0 ratio
   * @param H            Output: Hubble parameter
   * @return the error status
   */
  
  int hqiv_Hubble_from_density(
                               double gamma,
                               double rho_tot,
                               double G_eff_ratio,
                               double * H
                               );
  
  /**
   * Compute effective gravitational coupling G_eff(a)
   * 
   * Paper Section 2.3, Eq. 6:
   *   G(a) = G0 × (H(a)/H0)^α
   * 
   * Also includes Planck-suppressed renormalisation (Paper Section 4.3):
   *   G_eff(φ) = G0 / (1 + γ ℓ_P² φ²)
   * 
   * @param H            Input: Hubble parameter at scale factor a
   * @param H0           Input: Hubble parameter today
   * @param alpha        Input: varying G exponent
   * @param phi          Input: horizon field φ = cH
   * @param gamma        Input: thermodynamic coefficient
   * @param G_eff_ratio  Output: G_eff/G0
   * @return the error status
   */
  
  int hqiv_G_eff(
                 double H,
                 double H0,
                 double alpha,
                 double phi,
                 double gamma,
                 double * G_eff_ratio
                 );
  
  /**
   * Compute horizon field φ
   * 
   * Paper Section 3, Appendix A:
   *   φ ≡ 2c²/Θ_local
   * 
   * In homogeneous FLRW background:
   *   Θ_local = 2c/H  →  φ = cH
   * 
   * @param H     Input: Hubble parameter
   * @param phi   Output: horizon field
   * @return the error status
   */
  
  int hqiv_phi_from_H(
                      double H,
                      double * phi
                      );
  
  /**
   * Compute inertia reduction factor f(α, φ)
   * 
   * Paper Section 4.2, Eq. 11:
   *   f(α, φ) = max( α / (α + χφ/6), f_min )
   * 
   * This factor modifies the inertial mass in the Euler equation:
   *   m_i = m_g × f(α, φ)
   * 
   * @param alpha      Input: varying G exponent
   * @param phi        Input: horizon field
   * @param chi        Input: acceleration scaling factor
   * @param f_min      Input: thermodynamic floor
   * @param f_inertia  Output: inertia reduction factor
   * @return the error status
   */
  
  int hqiv_inertia_factor(
                          double alpha,
                          double phi,
                          double chi,
                          double f_min,
                          double * f_inertia
                          );
  
  /**
   * Compute derivative of inertia factor with respect to φ
   * 
   * Used in vorticity source term (Paper Section 6):
   *   ∂f/∂φ = -αχ / [6(α + χφ/6)²]
   * 
   * @param alpha    Input: varying G exponent
   * @param phi      Input: horizon field
   * @param chi      Input: acceleration scaling factor
   * @param df_dphi  Output: derivative df/dφ
   * @return the error status
   */
  
  int hqiv_df_dphi(
                   double alpha,
                   double phi,
                   double chi,
                   double * df_dphi
                   );
  
  /**
   * Compute vorticity source term
   * 
   * Paper Section 6:
   *   S_ω = (∂f/∂φ) (k × ∇φ) · ê_ω
   * 
   * In linear perturbation theory, ∇φ is computed from metric perturbations.
   * 
   * @param df_dphi     Input: derivative of inertia factor
   * @param grad_phi    Input: gradient of horizon field
   * @param k_vector    Input: wave vector
   * @param source      Output: vorticity source term
   * @return the error status
   */
  
  int hqiv_vorticity_source(
                            double df_dphi,
                            double grad_phi,
                            double k_vector,
                            double * source
                            );
  
  /**
   * Output HQIV-specific background quantities
   * 
   * @param phqiv       Input: pointer to HQIV parameters
   * @param titles      Output: column titles
   * @return the error status
   */
  
  int hqiv_output_titles(
                         struct hqiv_parameters * phqiv,
                         char titles[_MAXTITLESTRINGLENGTH_]
                         );
  
  /**
   * Output HQIV-specific background data
   * 
   * @param phqiv           Input: pointer to HQIV parameters
   * @param number_of_titles Input: number of columns
   * @param data            Output: data array
   * @return the error status
   */
  
  int hqiv_output_data(
                       struct hqiv_parameters * phqiv,
                       int number_of_titles,
                       double * data
                       );

#ifdef __cplusplus
}
#endif

/**************************************************************/

/* Physical constants used in HQIV calculations */

#define _HQIV_c_ 2.99792458e8           /**< Speed of light [m/s] */
#define _HQIV_G0_ 6.67430e-11           /**< Newton constant [m³ kg⁻¹ s⁻²] */
#define _HQIV_hbar_ 1.054571817e-34     /**< Reduced Planck constant [J·s] */
#define _HQIV_lP_ 1.616255e-35          /**< Planck length [m] */

/* Default parameter values (Paper Table 1) */

#define _HQIV_GAMMA_DEFAULT_ 0.40       /**< Default thermodynamic coefficient */
#define _HQIV_ALPHA_DEFAULT_ 0.60       /**< Default varying G exponent */
#define _HQIV_CHI_DEFAULT_ 0.172        /**< Default acceleration scaling factor */
#define _HQIV_FMIN_DEFAULT_ 0.01        /**< Default thermodynamic floor */

/* Parameter ranges for validation */

#define _HQIV_GAMMA_MIN_ 0.0            /**< Minimum gamma (0 = standard ΛCDM) */
#define _HQIV_GAMMA_MAX_ 1.0            /**< Maximum gamma */
#define _HQIV_ALPHA_MIN_ 0.0            /**< Minimum alpha */
#define _HQIV_ALPHA_MAX_ 2.0            /**< Maximum alpha */
#define _HQIV_CHI_MIN_ 0.0              /**< Minimum chi */
#define _HQIV_CHI_MAX_ 1.0              /**< Maximum chi */
#define _HQIV_FMIN_MIN_ 0.001           /**< Minimum f_min */
#define _HQIV_FMIN_MAX_ 0.1             /**< Maximum f_min */

#endif /* __HQIV__ */
/* @endcond */
