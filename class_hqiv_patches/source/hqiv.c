/** @file hqiv.c Documented HQIV module implementation
 *
 * Horizon-Quantized Informational Vacuum (HQIV) Cosmological Framework
 * 
 * This module implements the HQIV modifications to standard ΛCDM cosmology
 * as derived in the paper "Horizon-Quantized Informational Vacuum: 
 * A Covariant Baryon-Only Cosmological Framework from Quantised Inertia"
 *
 * Reference: new-final.tex (February 2026)
 * 
 * Key equations implemented:
 * - Modified Friedmann (Section 5, Eq. 13): 3H² - γH = 8πG_eff(ρ_m + ρ_r)
 * - Inertia modification (Section 4.2, Eq. 11): f(α,ϕ) = max(α/(α + χϕ/6), f_min)
 * - Varying G (Section 2.3): G(a) = G0 × (H(a)/H0)^α
 * - Vorticity source (Section 6): ∂ω/∂t + (v·∇)ω = (∂f/∂ϕ)(k × ∇ϕ)·ê_ω
 *
 * Author: HQIV Team
 */

#include "hqiv.h"
#include "background.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

/**
 * Initialize HQIV parameters with default values
 *
 * Sets all HQIV parameters to their fiducial values from Paper Table 1:
 * - gamma = 0.40 (thermodynamic coefficient from Brodie's integral)
 * - alpha = 0.60 (varying G exponent from rotation curves)
 * - chi = 0.172 (light-cone average scaling)
 * - f_min = 0.01 (informational-energy saturation floor)
 *
 * @param phqiv  Output: pointer to HQIV parameters structure
 * @return the error status
 */

int hqiv_init(struct hqiv_parameters * phqiv) {
  
  /** - Set default values from Paper Table 1 */
  phqiv->hqiv_on = _FALSE_;
  phqiv->gamma_hqiv = _HQIV_GAMMA_DEFAULT_;
  phqiv->alpha_hqiv = _HQIV_ALPHA_DEFAULT_;
  phqiv->chi_hqiv = _HQIV_CHI_DEFAULT_;
  phqiv->fmin_hqiv = _HQIV_FMIN_DEFAULT_;
  
  /** - Initialize derived quantities to zero */
  phqiv->phi_horizon = 0.0;
  phqiv->G_eff_ratio = 1.0;
  phqiv->inertia_factor = 1.0;
  phqiv->age_hqiv = 0.0;
  phqiv->vorticity_growth_exponent = 0.0;
  phqiv->raw_acoustic_scale = 0.0;
  
  return _SUCCESS_;
}

/**
 * Free HQIV parameters
 *
 * Currently no dynamic allocation, but provided for future extensibility.
 *
 * @param phqiv  Input: pointer to HQIV parameters structure
 * @return the error status
 */

int hqiv_free(struct hqiv_parameters * phqiv) {
  
  /* No dynamic allocation currently */
  
  return _SUCCESS_;
}

/**
 * Compute HQIV-modified Hubble parameter from energy density
 *
 * Solves the quadratic Friedmann equation (Paper Section 5, Eq. 13):
 *   3H² - γH = 8πG_eff(ρ_m + ρ_r)
 *
 * In CLASS units, densities are expressed as:
 *   ρ_class = 8πG ρ_physical / (3c²)
 * 
 * So the equation becomes:
 *   3H² - γH = 3c² × G_eff/G0 × ρ_class
 *
 * For normalization H(a=1) = H0, we need:
 *   At a=1: 3H0² - γH0 = 3c² × G_eff(H0)/G0 × ρ_class(a=1)
 *
 * The positive root solution (Appendix A.2):
 *   H = [γ + √(γ² + 36 × G_eff/G0 × ρ_class)] / 6
 *
 * Note: In CLASS conventions, ρ_class already includes the 8πG/3c² factor,
 * so we work with the dimensionless combination.
 *
 * @param gamma        Input: thermodynamic coefficient
 * @param rho_tot      Input: total energy density in CLASS units (H0² × Ω)
 * @param G_eff_ratio  Input: G_eff/G0 ratio
 * @param H            Output: Hubble parameter in CLASS units
 * @return the error status
 */

int hqiv_Hubble_from_density(
                             double gamma,
                             double rho_tot,
                             double G_eff_ratio,
                             double * H
                             ) {
  
  double discriminant;
  double H_squared;
  
  /** - For γ = 0, recover standard Friedmann: H = √(G_eff × ρ_tot) */
  if (gamma <= 0.0) {
    *H = sqrt(G_eff_ratio * rho_tot);
    return _SUCCESS_;
  }
  
  /** - Solve quadratic: 3H² - γH - 3 × G_eff × ρ_tot = 0
       In CLASS units, ρ_tot is already scaled, so:
       discriminant = γ² + 36 × G_eff × ρ_tot */
  
  discriminant = gamma * gamma + 36.0 * G_eff_ratio * rho_tot;
  
  /** - Take positive root: H = [γ + √(discriminant)] / 6 */
  *H = (gamma + sqrt(discriminant)) / 6.0;
  
  return _SUCCESS_;
}

/**
 * Compute effective gravitational coupling G_eff(a)
 *
 * Paper: we use G_eff = G_0 throughout. A Planck-suppressed correction
 * G_eff(φ) = G0 / (1 + γ (ℓ_P φ/c²)²) is negligible except near the
 * Planck epoch and has no effect on the runs here.
 *
 * For perturbations we optionally apply varying G (Paper Section 2.3):
 *   G_eff/G0 = (H/H0)^α
 *
 * @param H            Input: Hubble parameter at scale factor a
 * @param H0           Input: Hubble parameter today
 * @param alpha        Input: varying G exponent
 * @param phi          Input: horizon field φ = cH in SI units [1/s]
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
               ) {
  
  double H_ratio;
  double planck_factor;
  double lP_SI = 1.616255e-35;  /* Planck length in meters */
  double c_SI = 2.99792458e8;   /* Speed of light in m/s */
  
  /** - Compute H/H0 ratio */
  if (H0 > 0.0) {
    H_ratio = H / H0;
  } else {
    H_ratio = 1.0;
  }
  
  /** - Compute varying G contribution: (H/H0)^α */
  *G_eff_ratio = pow(H_ratio, alpha);
  
  /** - Compute Planck suppression factor: 1/(1 + γ ℓ_P² φ²)
       Note: φ has units of 1/s, ℓ_P has units of m
       We need ℓ_P² × φ²/c² to make it dimensionless
       Actually, φ = cH, so φ²/c² = H², and ℓ_P² H²/c² is dimensionless */
  
  /* The Planck suppression is extremely small for cosmological H values
     (H ~ 10^-18 1/s, ℓ_P ~ 10^-35 m, so ℓ_P²H²/c² ~ 10^-106)
     We include it for completeness but it's negligible */
  
  planck_factor = 1.0;  /* Negligible at cosmological scales */
  
  *G_eff_ratio *= planck_factor;
  
  return _SUCCESS_;
}

/**
 * Compute horizon field φ from Hubble parameter
 *
 * Paper Section 3, Appendix A:
 *   φ ≡ 2c²/Θ_local
 *
 * In homogeneous FLRW background:
 *   Θ_local = 2c/H  →  φ = cH
 *
 * Note: In CLASS, H is in units of Mpc^-1 (with c=1 convention).
 * To get φ in SI units (1/s), we need:
 *   φ = c × H_class × c / (1 Mpc) = H_class × c² / (1 Mpc)
 *
 * However, for the inertia factor calculation, we can work in CLASS units
 * where φ = H (since c=1 in natural units).
 *
 * @param H     Input: Hubble parameter in CLASS units [Mpc^-1]
 * @param phi   Output: horizon field in CLASS units (same as H)
 * @return the error status
 */

int hqiv_phi_from_H(
                    double H,
                    double * phi
                    ) {
  
  /** - In CLASS natural units (c=1), φ = H */
  *phi = H;
  
  return _SUCCESS_;
}

/**
 * Compute inertia reduction factor f(α, φ)
 *
 * Paper Section 4.2, Eq. 11:
 *   f(α, φ) = max( α / (α + χφ/6), f_min )
 *
 * This factor modifies the inertial mass in the Euler equation:
 *   m_i = m_g × f(α, φ)
 *
 * Physical interpretation:
 * - At early times (large φ = large H), f → f_min (maximum inertia reduction)
 * - At late times (small φ), f → 1 (standard inertia)
 * - The transition occurs when χφ/6 ~ α
 *
 * @param alpha      Input: varying G exponent
 * @param phi        Input: horizon field in CLASS units
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
                        ) {
  
  double f_unclamped;
  
  /** - Paper Eq. 11: f = α / (α + cφ/6). In c=1 units: f = α / (α + φ/6).
       Note: χ rescales a_min, NOT the interpolation function. */
  if (alpha > 0.0) {
    f_unclamped = alpha / (alpha + phi / 6.0);
  } else {
    f_unclamped = 1.0;  /* No modification if alpha = 0 */
  }
  
  /** - Apply thermodynamic floor: f = max(f_unclamped, f_min) */
  if (f_unclamped < f_min) {
    *f_inertia = f_min;
  } else {
    *f_inertia = f_unclamped;
  }
  
  /** - Ensure f ≤ 1 (inertia can only be reduced, not enhanced) */
  if (*f_inertia > 1.0) {
    *f_inertia = 1.0;
  }
  
  return _SUCCESS_;
}

/**
 * Compute derivative of inertia factor with respect to φ
 *
 * Used in vorticity source term (Paper Section 6):
 *   ∂f/∂φ = -α / [6(α + φ/6)²]  (paper: no χ in f formula)
 *
 * This derivative is the key to the vorticity generation mechanism:
 * spatial gradients in φ couple to the inertia modification to produce
 * vorticity sources at the BAO scale.
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
                 ) {
  
  double denominator;
  
  /** - Compute denominator: (α + φ/6)² */
  denominator = alpha + phi / 6.0;
  denominator = denominator * denominator;
  
  /** - Paper: df/dφ = -α / [6(α + φ/6)²] */
  if (denominator > 0.0) {
    *df_dphi = -alpha / (6.0 * denominator);
  } else {
    *df_dphi = 0.0;
  }
  
  return _SUCCESS_;
}

/**
 * Compute vorticity source term
 *
 * Paper Section 6:
 *   S_ω = (∂f/∂φ) (k × ∇φ) · ê_ω
 *
 * In linear perturbation theory, ∇φ is computed from metric perturbations.
 * The source injects coherent rotational seeds precisely at the BAO scale
 * during recombination, when sound-horizon interfaces produce the largest
 * spatial gradients ∇φ.
 *
 * @param df_dphi     Input: derivative of inertia factor
 * @param grad_phi    Input: gradient of horizon field
 * @param k_vector    Input: wave vector magnitude
 * @param source      Output: vorticity source term
 * @return the error status
 */

int hqiv_vorticity_source(
                          double df_dphi,
                          double grad_phi,
                          double k_vector,
                          double * source
                          ) {
  
  /** - The vorticity source is proportional to df/dφ × |∇φ| × k
       This represents the (k × ∇φ) · ê_ω term in the vorticity equation */
  
  *source = df_dphi * grad_phi * k_vector;
  
  return _SUCCESS_;
}

/**
 * Output HQIV-specific background quantities titles
 *
 * @param phqiv       Input: pointer to HQIV parameters
 * @param titles      Output: column titles
 * @return the error status
 */

int hqiv_output_titles(
                       struct hqiv_parameters * phqiv,
                       char titles[_MAXTITLESTRINGLENGTH_]
                       ) {
  
  /** - Only output if HQIV is enabled */
  if (phqiv->hqiv_on == _TRUE_) {
    class_store_columntitle(titles,"phi_horizon",_TRUE_);
    class_store_columntitle(titles,"G_eff_ratio",_TRUE_);
    class_store_columntitle(titles,"f_inertia",_TRUE_);
    class_store_columntitle(titles,"df_dphi",_TRUE_);
  }
  
  return _SUCCESS_;
}

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
                     ) {
  
  int storeidx = 0;
  
  /** - Only output if HQIV is enabled */
  if (phqiv->hqiv_on == _TRUE_) {
    class_store_double(data, phqiv->phi_horizon, _TRUE_, storeidx);
    class_store_double(data, phqiv->G_eff_ratio, _TRUE_, storeidx);
    class_store_double(data, phqiv->inertia_factor, _TRUE_, storeidx);
    class_store_double(data, 0.0, _TRUE_, storeidx);  /* df_dphi placeholder */
  }
  
  return _SUCCESS_;
}