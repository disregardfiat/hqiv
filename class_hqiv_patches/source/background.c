/** @file background.c Documented background module
 *
 * * Julien Lesgourgues, 17.04.2011
 * * routines related to ncdm written by T. Tram in 2011
 * * new integration scheme written by N. Schoeneberg in 2020
 * * HQIV modifications added 2026
 *
 * Deals with the cosmological background evolution.
 * 
 * HQIV Modification: The Friedmann equation is modified to:
 *   3H² - γH = 8πG_eff ρ
 * 
 * Solution: H = [γ + √(γ² + 24πG_eff ρ)] / 6
 */

#include "background.h"

/**
 * Background quantities at given redshift z.
 */

int background_at_z(
                    struct background *pba,
                    double z,
                    enum vecback_format return_format,
                    enum interpolation_method inter_mode,
                    int * last_index,
                    double * pvecback
                    ) {

  int pvecback_size;
  double loga;

  loga = -log(1+z);

  class_test(loga < pba->loga_table[0],
             pba->error_message,
             "out of range: a/a_0 = %e < a_min/a_0 = %e, you should decrease the precision parameter a_ini_over_a_today_default\n",1./(1.+z),exp(pba->loga_table[0]));

  class_test(loga > pba->loga_table[pba->bt_size-1],
             pba->error_message,
             "out of range: a/a_0 = %e > a_max/a_0 = %e\n",1./(1.+z),exp(pba->loga_table[pba->bt_size-1]));

  if (return_format == normal_info) {
    pvecback_size=pba->bg_size_normal;
  }
  else {
    if (return_format == short_info) {
      pvecback_size=pba->bg_size_short;
    }
    else {
      pvecback_size=pba->bg_size;
    }
  }

  if (inter_mode == inter_normal) {
    class_call(array_interpolate_spline(
                                        pba->loga_table,
                                        pba->bt_size,
                                        pba->background_table,
                                        pba->d2background_dloga2_table,
                                        pba->bg_size,
                                        loga,
                                        last_index,
                                        pvecback,
                                        pvecback_size,
                                        pba->error_message),
               pba->error_message,
               pba->error_message);
  }
  if (inter_mode == inter_closeby) {
    class_call(array_interpolate_spline_growing_closeby(
                                                        pba->loga_table,
                                                        pba->bt_size,
                                                        pba->background_table,
                                                        pba->d2background_dloga2_table,
                                                        pba->bg_size,
                                                        loga,
                                                        last_index,
                                                        pvecback,
                                                        pvecback_size,
                                                        pba->error_message),
               pba->error_message,
               pba->error_message);
  }

  return _SUCCESS_;
}

/**
 * Background quantities at given conformal time tau.
 */

int background_at_tau(
                      struct background *pba,
                      double tau,
                      enum vecback_format return_format,
                      enum interpolation_method inter_mode,
                      int * last_index,
                      double * pvecback
                      ) {

  double z;

  class_call(background_z_of_tau(pba,tau,&z),
             pba->error_message,
             pba->error_message);

  class_call(background_at_z(pba,z,return_format,inter_mode,last_index,pvecback),
             pba->error_message,
             pba->error_message);

  return _SUCCESS_;
}

/**
 * Conformal time at given redshift.
 */

int background_tau_of_z(
                        struct background *pba,
                        double z,
                        double * tau
                        ) {

  int last_index;

  class_test(z < pba->z_table[pba->bt_size-1],
             pba->error_message,
             "out of range: z=%e < z_min=%e\n",z,pba->z_table[pba->bt_size-1]);

  class_test(z > pba->z_table[0],
             pba->error_message,
             "out of range: z=%e > z_max=%e\n",z,pba->z_table[0]);

  class_call(array_interpolate_spline(
                                      pba->z_table,
                                      pba->bt_size,
                                      pba->tau_table,
                                      pba->d2tau_dz2_table,
                                      1,
                                      z,
                                      &last_index,
                                      tau,
                                      1,
                                      pba->error_message),
             pba->error_message,
             pba->error_message);

  return _SUCCESS_;
}

/**
 * Redshift at given conformal time.
 */

int background_z_of_tau(
                        struct background *pba,
                        double tau,
                        double * z
                        ) {

  int last_index;

  class_test(tau < pba->tau_table[0],
             pba->error_message,
             "out of range: tau=%e < tau_min=%e\n",tau,pba->tau_table[0]);

  class_test(tau > pba->tau_table[pba->bt_size-1],
             pba->error_message,
             "out of range: tau=%e > tau_max=%e\n",tau,pba->tau_table[pba->bt_size-1]);

  class_call(array_interpolate_spline(
                                      pba->tau_table,
                                      pba->bt_size,
                                      pba->z_table,
                                      pba->d2z_dtau2_table,
                                      1,
                                      tau,
                                      &last_index,
                                      z,
                                      1,
                                      pba->error_message),
             pba->error_message,
             pba->error_message);

  return _SUCCESS_;
}

/**
 * Function evaluating all background quantities which can be computed
 * analytically as a function of a and of {B} quantities.
 *
 * HQIV MODIFICATION: The Friedmann equation is modified to:
 *   3H² - γH = 8πG_eff ρ
 * 
 * Solution: H = [γ + √(γ² + 24πG_eff ρ)] / 6
 */

int background_functions(
                         struct background * pba,
                         double a,
                         double * pvecback_B,
                         enum vecback_format return_format,
                         double * pvecback
                         ) {

  double rho_tot;
  double rho_crit;
  double p_tot;
  double rho_r;
  double rho_m;
  double rho_ncdm,p_ncdm,pseudo_p_ncdm;
  int n_ncdm;
  double w_fld, dw_over_da, integral_fld;
  double phi, phi_prime;
  double dp_dloga;
  
  /* HQIV-specific variables */
  double gamma_hqiv, alpha_hqiv;
  double G_eff_ratio;
  double H_standard;
  double gamma_eff;  /* epoch-dependent γ: 0 pre-recombination, γ_theory post-recombination */

  rho_tot = 0.;
  p_tot = 0.;
  dp_dloga = 0.;
  rho_r=0.;
  rho_m=0.;

  /** HQIV emergent lattice: densities come from pre-filled table only; interpolate from table */
  if (pba->hqiv_emergent == _TRUE_) {
    int last_index = 0;
    double loga = log(a);
    int pvecback_size = (return_format == normal_info) ? pba->bg_size_normal :
                       ((return_format == short_info) ? pba->bg_size_short : pba->bg_size);
    class_call(array_interpolate_spline(pba->loga_table,
                                        pba->bt_size,
                                        pba->background_table,
                                        pba->d2background_dloga2_table,
                                        pba->bg_size,
                                        loga,
                                        &last_index,
                                        pvecback,
                                        pvecback_size,
                                        pba->error_message),
               pba->error_message,
               pba->error_message);
    return _SUCCESS_;
  }

  class_test(a <= 0.,
             pba->error_message,
             "a = %e instead of strictly positive",a);

  pvecback[pba->index_bg_a] = a;

  /* photons */
  pvecback[pba->index_bg_rho_g] = pba->Omega0_g * pow(pba->H0,2) / pow(a,4);
  rho_tot += pvecback[pba->index_bg_rho_g];
  p_tot += (1./3.) * pvecback[pba->index_bg_rho_g];
  dp_dloga += -(4./3.) * pvecback[pba->index_bg_rho_g];
  rho_r += pvecback[pba->index_bg_rho_g];

  /* baryons */
  pvecback[pba->index_bg_rho_b] = pba->Omega0_b * pow(pba->H0,2) / pow(a,3);
  rho_tot += pvecback[pba->index_bg_rho_b];
  p_tot += 0;
  rho_m += pvecback[pba->index_bg_rho_b];

  /* cdm */
  if (pba->has_cdm == _TRUE_) {
    pvecback[pba->index_bg_rho_cdm] = pba->Omega0_cdm * pow(pba->H0,2) / pow(a,3);
    rho_tot += pvecback[pba->index_bg_rho_cdm];
    p_tot += 0.;
    rho_m += pvecback[pba->index_bg_rho_cdm];
  }

  /* idm */
  if (pba->has_idm == _TRUE_) {
    pvecback[pba->index_bg_rho_idm] = pba->Omega0_idm * pow(pba->H0,2) / pow(a,3);
    rho_tot += pvecback[pba->index_bg_rho_idm];
    p_tot += 0.;
    rho_m += pvecback[pba->index_bg_rho_idm];
  }

  /* dcdm */
  if (pba->has_dcdm == _TRUE_) {
    pvecback[pba->index_bg_rho_dcdm] = pvecback_B[pba->index_bi_rho_dcdm];
    rho_tot += pvecback[pba->index_bg_rho_dcdm];
    p_tot += 0.;
    rho_m += pvecback[pba->index_bg_rho_dcdm];
  }

  /* dr */
  if (pba->has_dr == _TRUE_) {
    pvecback[pba->index_bg_rho_dr] = pvecback_B[pba->index_bi_rho_dr];
    rho_tot += pvecback[pba->index_bg_rho_dr];
    p_tot += (1./3.)*pvecback[pba->index_bg_rho_dr];
    dp_dloga += -(4./3.) * pvecback[pba->index_bg_rho_dr];
    rho_r += pvecback[pba->index_bg_rho_dr];
  }

  /* Scalar field */
  if (pba->has_scf == _TRUE_) {
    phi = pvecback_B[pba->index_bi_phi_scf];
    phi_prime = pvecback_B[pba->index_bi_phi_prime_scf];
    pvecback[pba->index_bg_phi_scf] = phi;
    pvecback[pba->index_bg_phi_prime_scf] = phi_prime;
    pvecback[pba->index_bg_V_scf] = V_scf(pba,phi);
    pvecback[pba->index_bg_dV_scf] = dV_scf(pba,phi);
    pvecback[pba->index_bg_ddV_scf] = ddV_scf(pba,phi);
    pvecback[pba->index_bg_rho_scf] = (phi_prime*phi_prime/(2*a*a) + V_scf(pba,phi))/3.;
    pvecback[pba->index_bg_p_scf] =(phi_prime*phi_prime/(2*a*a) - V_scf(pba,phi))/3.;
    rho_tot += pvecback[pba->index_bg_rho_scf];
    p_tot += pvecback[pba->index_bg_p_scf];
    dp_dloga += 0.0;
    rho_r += 3.*pvecback[pba->index_bg_p_scf];
    rho_m += pvecback[pba->index_bg_rho_scf] - 3.* pvecback[pba->index_bg_p_scf];
  }

  /* ncdm */
  if (pba->has_ncdm == _TRUE_) {
    for (n_ncdm=0; n_ncdm<pba->N_ncdm; n_ncdm++) {
      class_call(background_ncdm_momenta(
                                         pba->q_ncdm_bg[n_ncdm],
                                         pba->w_ncdm_bg[n_ncdm],
                                         pba->q_size_ncdm_bg[n_ncdm],
                                         pba->M_ncdm[n_ncdm],
                                         pba->factor_ncdm[n_ncdm],
                                         1./a-1.,
                                         NULL,
                                         &rho_ncdm,
                                         &p_ncdm,
                                         NULL,
                                         &pseudo_p_ncdm),
                 pba->error_message,
                 pba->error_message);

      pvecback[pba->index_bg_rho_ncdm1+n_ncdm] = rho_ncdm;
      rho_tot += rho_ncdm;
      pvecback[pba->index_bg_p_ncdm1+n_ncdm] = p_ncdm;
      p_tot += p_ncdm;
      pvecback[pba->index_bg_pseudo_p_ncdm1+n_ncdm] = pseudo_p_ncdm;
      dp_dloga += (pseudo_p_ncdm - 5*p_ncdm);
      rho_r += 3.* p_ncdm;
      rho_m += rho_ncdm - 3.* p_ncdm;
    }
  }

  /* Lambda */
  if (pba->has_lambda == _TRUE_) {
    pvecback[pba->index_bg_rho_lambda] = pba->Omega0_lambda * pow(pba->H0,2);
    rho_tot += pvecback[pba->index_bg_rho_lambda];
    p_tot -= pvecback[pba->index_bg_rho_lambda];
  }

  /* fluid with w(a) and constant cs2 */
  if (pba->has_fld == _TRUE_) {
    pvecback[pba->index_bg_rho_fld] = pvecback_B[pba->index_bi_rho_fld];
    class_call(background_w_fld(pba,a,&w_fld,&dw_over_da,&integral_fld), pba->error_message, pba->error_message);
    pvecback[pba->index_bg_w_fld] = w_fld;
    rho_tot += pvecback[pba->index_bg_rho_fld];
    p_tot += w_fld * pvecback[pba->index_bg_rho_fld];
    dp_dloga += (a*dw_over_da-3*(1+w_fld)*w_fld)*pvecback[pba->index_bg_rho_fld];
  }

  /* relativistic neutrinos (and all relativistic relics) */
  if (pba->has_ur == _TRUE_) {
    pvecback[pba->index_bg_rho_ur] = pba->Omega0_ur * pow(pba->H0,2) / pow(a,4);
    rho_tot += pvecback[pba->index_bg_rho_ur];
    p_tot += (1./3.) * pvecback[pba->index_bg_rho_ur];
    dp_dloga += -(4./3.) * pvecback[pba->index_bg_rho_ur];
    rho_r += pvecback[pba->index_bg_rho_ur];
  }

  /* interacting dark radiation */
  if (pba->has_idr == _TRUE_) {
    pvecback[pba->index_bg_rho_idr] = pba->Omega0_idr * pow(pba->H0,2) / pow(a,4);
    rho_tot += pvecback[pba->index_bg_rho_idr];
    p_tot += (1./3.) * pvecback[pba->index_bg_rho_idr];
    rho_r += pvecback[pba->index_bg_rho_idr];
  }

  /**
   * HQIV MODIFICATION: Modified Friedmann equation with epoch-dependent γ_eff(a).
   *
   * Physical principle: Before recombination the baryon-photon plasma is tightly
   * coupled; radiation pressure + Thomson scattering wash out local horizon
   * anisotropies → γ_eff = 0.  After recombination baryons decouple and each
   * matter parcel feels its own causal horizon → γ_eff = γ_theory = 0.40 (fixed
   * by Brodie's backward-hemisphere overlap integral / discrete hockey-stick
   * combinatorics).
   *
   * Pre-recombination:  3H² = 8πG ρ                (standard Friedmann)
   * Post-recombination: 3H² - γ_eff H = 8πG ρ     (HQIV modified)
   *
   * The γH term acts as an effective horizon energy density that drives
   * late-time acceleration without a cosmological constant.
   *
   * G_eff is derived exactly from the Friedmann equation (no α parameter):
   *   G_eff/G_0 = H(3H - γ_eff_class) / [H_0(3H_0 - γ_eff0_class)]
   *
   * CLASS units: H in 1/Mpc, rho_tot = (8πG/3) ρ_physical.
   */

  H_standard = sqrt(rho_tot - pba->K/a/a);

  if (pba->hqiv.hqiv_on == _TRUE_) {
    double gamma_theory = pba->hqiv.gamma_hqiv;   /* 0.40 — fixed, not a free parameter */
    double a_rec = 1.0 / 1090.0;                  /* recombination scale factor */
    double da_width = 0.02 * a_rec;               /* smooth transition width (~few % of a_rec) */
    double rho_eff = rho_tot - pba->K/a/a;

    /* Smooth step: γ_eff(a) = γ_theory × ½[1 + tanh((a - a_rec)/Δa)] */
    gamma_eff = gamma_theory * 0.5 * (1.0 + tanh((a - a_rec) / da_width));

    {
      double gamma_class = gamma_eff * pba->H0;   /* dimensional: [Mpc⁻¹] */
      double b = gamma_class / 3.0;
      pvecback[pba->index_bg_H] = (b + sqrt(b*b + 4.0 * rho_eff)) / 2.0;
    }

    /* G_eff/G_0 from the exact Friedmann relation with curvature (paper §E.1) */
    class_call(hqiv_G_eff(pvecback[pba->index_bg_H], pba->H0, pba->hqiv.alpha_hqiv,
                          pvecback[pba->index_bg_H], gamma_eff,
                          pba->hqiv.gamma_hqiv,      /* γ_theory at z=0 */
                          pba->K / (a * a),           /* K/a² */
                          pba->K,                     /* K (for z=0 normalization) */
                          &(pba->hqiv.G_eff_ratio)),
               pba->error_message, pba->error_message);
    pba->hqiv.phi_horizon = pvecback[pba->index_bg_H] / pba->H0;

    /* α fully dynamic: α_eff = χφ/6 so f = χ/(χ+1); else use input alpha */
    {
      double alpha_eff = (pba->hqiv.alpha_dynamic == _TRUE_)
                        ? (pba->hqiv.chi_hqiv * pvecback[pba->index_bg_H] / 6.0)
                        : pba->hqiv.alpha_hqiv;
      hqiv_inertia_factor(alpha_eff, pvecback[pba->index_bg_H],
                          pba->hqiv.chi_hqiv, pba->hqiv.fmin_hqiv, &(pba->hqiv.inertia_factor));
    }
  }
  else {
    gamma_eff = 0.0;
    pvecback[pba->index_bg_H] = H_standard;
  }

  /** - compute derivative of H with respect to conformal time.
   * Standard GR: H' = -(3/2)(ρ+p)a + K/a  (from Raychaudhuri + Friedmann).
   * HQIV: differentiating 3H²-γ_eff H = 3(ρ-K/a²) w.r.t. conformal time τ gives
   *   (6H-γ_eff_class) H' = -9H(ρ+p)a + 6HK/a,  γ_eff_class = γ_eff × H0.
   * This reduces to the standard GR result when γ_eff → 0. */
  {
    double gamma_eff_class = gamma_eff * pba->H0;
    double H_now           = pvecback[pba->index_bg_H];
    double denom           = 6.0 * H_now - gamma_eff_class;
    if (denom > 1e-30 * pba->H0) {
      pvecback[pba->index_bg_H_prime] = (-9.0 * H_now * (rho_tot + p_tot) * a
                                         + 6.0 * H_now * pba->K / a) / denom;
    } else {
      pvecback[pba->index_bg_H_prime] = - (3./2.) * (rho_tot + p_tot) * a + pba->K/a;
    }
  }

  /* Total energy density*/
  pvecback[pba->index_bg_rho_tot] = rho_tot;

  /* Total pressure */
  pvecback[pba->index_bg_p_tot] = p_tot;

  /* Derivative of total pressure w.r.t. conformal time */
  pvecback[pba->index_bg_p_tot_prime] = a*pvecback[pba->index_bg_H]*dp_dloga;
  if (pba->has_scf == _TRUE_) {
    pvecback[pba->index_bg_p_prime_scf] = pvecback[pba->index_bg_phi_prime_scf]*
      (-pvecback[pba->index_bg_phi_prime_scf]*pvecback[pba->index_bg_H]/a-2./3.*pvecback[pba->index_bg_dV_scf]);
    pvecback[pba->index_bg_p_tot_prime] += pvecback[pba->index_bg_p_prime_scf];
  }

  /** - compute critical density.
   * Standard GR: rho_crit = H².
   * HQIV: From 3H² - γ_eff H = 8πG ρ, we have H² = ρ_tot + (γ_eff H_0 / 3) H.
   * The critical density is still H² so that Ω_m + Ω_r + Ω_horizon = 1. */
  if (pba->hqiv.hqiv_on == _TRUE_) {
    double H_here = pvecback[pba->index_bg_H];
    rho_crit = H_here * H_here;
  }
  else {
    rho_crit = rho_tot - pba->K/a/a;
  }
  class_test(rho_crit <= 0.,
             pba->error_message,
             "rho_crit = %e instead of strictly positive",rho_crit);

  /** - compute relativistic density to total density ratio */
  pvecback[pba->index_bg_Omega_r] = rho_r / rho_crit;

  /** - compute other quantities in the exhaustive, redundant format */
  if (return_format == long_info) {

    pvecback[pba->index_bg_rho_crit] = rho_crit;
    pvecback[pba->index_bg_Omega_m] = rho_m / rho_crit;
    pvecback[pba->index_bg_time] = pvecback_B[pba->index_bi_time];
    pvecback[pba->index_bg_rs] = pvecback_B[pba->index_bi_rs];
    pvecback[pba->index_bg_D] = pvecback_B[pba->index_bi_D];
    pvecback[pba->index_bg_f] = pvecback_B[pba->index_bi_D_prime]/( pvecback_B[pba->index_bi_D]*a*pvecback[pba->index_bg_H]);

    if (pba->has_varconst == _TRUE_) {
      class_call(background_varconst_of_z(pba,
                                          1./a-1.,
                                          &(pvecback[pba->index_bg_varc_alpha]),
                                          &(pvecback[pba->index_bg_varc_me])
                                          ),
                 pba->error_message,
                 pba->error_message);
    }
  }

  return _SUCCESS_;
}

/**
 * Single place where the fluid equation of state is defined.
 */

int background_w_fld(
                     struct background * pba,
                     double a,
                     double * w_fld,
                     double * dw_over_da_fld,
                     double * integral_fld
                     ) {

  double Omega_ede = 0.;
  double dOmega_ede_over_da = 0.;
  double d2Omega_ede_over_da2 = 0.;
  double a_eq, Omega_r, Omega_m;

  switch (pba->fluid_equation_of_state) {
  case CLP:
    *w_fld = pba->w0_fld + pba->wa_fld * (1. - a);
    break;
  case EDE:
    Omega_ede = (pba->Omega0_fld - pba->Omega_EDE*(1.-pow(a,-3.*pba->w0_fld)))
      /(pba->Omega0_fld+(1.-pba->Omega0_fld)*pow(a,3.*pba->w0_fld))
      + pba->Omega_EDE*(1.-pow(a,-3.*pba->w0_fld));
    dOmega_ede_over_da = - pba->Omega_EDE* 3.*pba->w0_fld*pow(a,-3.*pba->w0_fld-1.)/(pba->Omega0_fld+(1.-pba->Omega0_fld)*pow(a,3.*pba->w0_fld))
      - (pba->Omega0_fld - pba->Omega_EDE*(1.-pow(a,-3.*pba->w0_fld)))*(1.-pba->Omega0_fld)*3.*pba->w0_fld*pow(a,3.*pba->w0_fld-1.)/pow(pba->Omega0_fld+(1.-pba->Omega0_fld)*pow(a,3.*pba->w0_fld),2)
      + pba->Omega_EDE*3.*pba->w0_fld*pow(a,-3.*pba->w0_fld-1.);
    Omega_r = pba->Omega0_g * (1. + 3.044 * 7./8.*pow(4./11.,4./3.));
    Omega_m = pba->Omega0_b;
    if (pba->has_cdm == _TRUE_) Omega_m += pba->Omega0_cdm;
    if (pba->has_idm == _TRUE_) Omega_m += pba->Omega0_idm;
    a_eq = Omega_r/Omega_m;
    *w_fld = - dOmega_ede_over_da*a/Omega_ede/3./(1.-Omega_ede)+a_eq/3./(a+a_eq);
    break;
  }

  switch (pba->fluid_equation_of_state) {
  case CLP:
    *dw_over_da_fld = - pba->wa_fld;
    break;
  case EDE:
    d2Omega_ede_over_da2 = 0.;
    *dw_over_da_fld = - d2Omega_ede_over_da2*a/3./(1.-Omega_ede)/Omega_ede
      - dOmega_ede_over_da/3./(1.-Omega_ede)/Omega_ede
      + dOmega_ede_over_da*dOmega_ede_over_da*a/3./(1.-Omega_ede)/(1.-Omega_ede)/Omega_ede
      + a_eq/3./(a+a_eq)/(a+a_eq);
    break;
  }

  switch (pba->fluid_equation_of_state) {
  case CLP:
    *integral_fld = 3.*((1.+pba->w0_fld+pba->wa_fld)*log(1./a) + pba->wa_fld*(a-1.));
    break;
  case EDE:
    class_stop(pba->error_message,"EDE implementation not finished\n");
    break;
  }

  return _SUCCESS_;
}

/**
 * Single place where the variation of fundamental constants is defined.
 */

int background_varconst_of_z(
                             struct background* pba,
                             double z,
                             double* alpha,
                             double* me
                             ){

  switch(pba->varconst_dep){

  case varconst_none:
    *alpha = 1.;
    *me = 1.;
    break;

  case varconst_instant:
    if (z>pba->varconst_transition_redshift){
      *alpha = pba->varconst_alpha;
      *me = pba->varconst_me;
    }
    else{
      *alpha = 1.;
      *me = 1.;
    }
    break;
  }
  return _SUCCESS_;
}

/**
 * Initialize the background structure.
 */

int background_init(
                    struct precision * ppr,
                    struct background * pba
                    ) {

  if (pba->background_verbose > 0) {
    printf("Running CLASS version %s\n",_VERSION_);
    printf("Computing background\n");
  }

  class_test(pba->shooting_failed == _TRUE_,
             pba->error_message,
             "Shooting failed, try optimising input_get_guess(). Error message:\n\n%s",
             pba->shooting_error);

  class_call(background_indices(pba),
             pba->error_message,
             pba->error_message);

  class_call(background_checks(ppr,pba),
             pba->error_message,
             pba->error_message);

  class_call(background_solve(ppr,pba),
             pba->error_message,
             pba->error_message);

  class_call(background_find_equality(ppr,pba),
             pba->error_message,
             pba->error_message);

  class_call(background_output_budget(pba),
             pba->error_message,
             pba->error_message);

  pba->is_allocated = _TRUE_;

  return _SUCCESS_;

}

/**
 * Free all memory space allocated by background_init().
 */

int background_free(
                    struct background *pba
                    ) {

  class_call(background_free_noinput(pba),
             pba->error_message,
             pba->error_message);

  class_call(background_free_input(pba),
             pba->error_message,
             pba->error_message);

  pba->is_allocated = _FALSE_;

  return _SUCCESS_;
}

/**
 * Free only the memory space NOT allocated through input_read_parameters().
 */

int background_free_noinput(
                            struct background *pba
                            ) {

  free(pba->tau_table);
  free(pba->z_table);
  free(pba->loga_table);
  free(pba->d2tau_dz2_table);
  free(pba->d2z_dtau2_table);
  free(pba->background_table);
  free(pba->d2background_dloga2_table);

  return _SUCCESS_;
}

/**
 * Free pointers inside background structure which were allocated in input_read_parameters().
 */

int background_free_input(
                          struct background *pba
                          ) {

  int k;

  if (pba->Omega0_ncdm_tot != 0.) {
    for (k=0; k<pba->N_ncdm; k++) {
      free(pba->q_ncdm[k]);
      free(pba->w_ncdm[k]);
      free(pba->q_ncdm_bg[k]);
      free(pba->w_ncdm_bg[k]);
      free(pba->dlnf0_dlnq_ncdm[k]);
    }
    free(pba->ncdm_quadrature_strategy);
    free(pba->ncdm_input_q_size);
    free(pba->ncdm_qmax);
    free(pba->q_ncdm);
    free(pba->w_ncdm);
    free(pba->q_ncdm_bg);
    free(pba->w_ncdm_bg);
    free(pba->dlnf0_dlnq_ncdm);
    free(pba->q_size_ncdm);
    free(pba->q_size_ncdm_bg);
    free(pba->M_ncdm);
    free(pba->T_ncdm);
    free(pba->ksi_ncdm);
    free(pba->deg_ncdm);
    free(pba->Omega0_ncdm);
    free(pba->m_ncdm_in_eV);
    free(pba->factor_ncdm);
    if (pba->got_files!=NULL)
      free(pba->got_files);
    if (pba->ncdm_psd_files!=NULL)
      free(pba->ncdm_psd_files);
    if (pba->ncdm_psd_parameters!=NULL)
      free(pba->ncdm_psd_parameters);
  }

  if (pba->Omega0_scf != 0.) {
    if (pba->scf_parameters != NULL)
      free(pba->scf_parameters);
  }
  return _SUCCESS_;
}

/**
 * Assign value to each relevant index in vectors of background quantities.
 */

int background_indices(
                       struct background *pba
                       ) {

  int index_bg;
  int index_bi;

  pba->has_cdm = _FALSE_;
  pba->has_idm = _FALSE_;
  pba->has_ncdm = _FALSE_;
  pba->has_dcdm = _FALSE_;
  pba->has_dr = _FALSE_;
  pba->has_scf = _FALSE_;
  pba->has_lambda = _FALSE_;
  pba->has_fld = _FALSE_;
  pba->has_ur = _FALSE_;
  pba->has_idr = _FALSE_;
  pba->has_curvature = _FALSE_;
  pba->has_varconst  = _FALSE_;

  if (pba->Omega0_cdm != 0.)
    pba->has_cdm = _TRUE_;

  if (pba->Omega0_idm != 0.)
    pba->has_idm = _TRUE_;

  if (pba->Omega0_ncdm_tot != 0.)
    pba->has_ncdm = _TRUE_;

  if (pba->Omega0_dcdmdr != 0.) {
    pba->has_dcdm = _TRUE_;
    if (pba->Gamma_dcdm != 0.)
      pba->has_dr = _TRUE_;
  }

  if (pba->Omega0_scf != 0.)
    pba->has_scf = _TRUE_;

  if (pba->Omega0_lambda != 0.)
    pba->has_lambda = _TRUE_;

  if (pba->Omega0_fld != 0.)
    pba->has_fld = _TRUE_;

  if (pba->Omega0_ur != 0.)
    pba->has_ur = _TRUE_;

  if (pba->Omega0_idr != 0.)
    pba->has_idr = _TRUE_;

  if (pba->sgnK != 0)
    pba->has_curvature = _TRUE_;

  if (pba->varconst_dep != varconst_none)
    pba->has_varconst = _TRUE_;

  index_bg=0;

  class_define_index(pba->index_bg_a,_TRUE_,index_bg,1);
  class_define_index(pba->index_bg_H,_TRUE_,index_bg,1);
  class_define_index(pba->index_bg_H_prime,_TRUE_,index_bg,1);

  pba->bg_size_short = index_bg;

  class_define_index(pba->index_bg_rho_g,_TRUE_,index_bg,1);
  class_define_index(pba->index_bg_rho_b,_TRUE_,index_bg,1);
  class_define_index(pba->index_bg_rho_cdm,pba->has_cdm,index_bg,1);
  class_define_index(pba->index_bg_rho_idm,pba->has_idm,index_bg,1);
  class_define_index(pba->index_bg_rho_ncdm1,pba->has_ncdm,index_bg,pba->N_ncdm);
  class_define_index(pba->index_bg_p_ncdm1,pba->has_ncdm,index_bg,pba->N_ncdm);
  class_define_index(pba->index_bg_pseudo_p_ncdm1,pba->has_ncdm,index_bg,pba->N_ncdm);
  class_define_index(pba->index_bg_rho_dcdm,pba->has_dcdm,index_bg,1);
  class_define_index(pba->index_bg_rho_dr,pba->has_dr,index_bg,1);
  class_define_index(pba->index_bg_phi_scf,pba->has_scf,index_bg,1);
  class_define_index(pba->index_bg_phi_prime_scf,pba->has_scf,index_bg,1);
  class_define_index(pba->index_bg_V_scf,pba->has_scf,index_bg,1);
  class_define_index(pba->index_bg_dV_scf,pba->has_scf,index_bg,1);
  class_define_index(pba->index_bg_ddV_scf,pba->has_scf,index_bg,1);
  class_define_index(pba->index_bg_rho_scf,pba->has_scf,index_bg,1);
  class_define_index(pba->index_bg_p_scf,pba->has_scf,index_bg,1);
  class_define_index(pba->index_bg_p_prime_scf,pba->has_scf,index_bg,1);
  class_define_index(pba->index_bg_rho_lambda,pba->has_lambda,index_bg,1);
  class_define_index(pba->index_bg_rho_fld,pba->has_fld,index_bg,1);
  class_define_index(pba->index_bg_w_fld,pba->has_fld,index_bg,1);
  class_define_index(pba->index_bg_rho_ur,pba->has_ur,index_bg,1);
  class_define_index(pba->index_bg_rho_tot,_TRUE_,index_bg,1);
  class_define_index(pba->index_bg_p_tot,_TRUE_,index_bg,1);
  class_define_index(pba->index_bg_p_tot_prime,_TRUE_,index_bg,1);
  class_define_index(pba->index_bg_Omega_r,_TRUE_,index_bg,1);
  class_define_index(pba->index_bg_rho_idr,pba->has_idr,index_bg,1);

  pba->bg_size_normal = index_bg;

  class_define_index(pba->index_bg_rho_crit,_TRUE_,index_bg,1);
  class_define_index(pba->index_bg_Omega_m,_TRUE_,index_bg,1);
  class_define_index(pba->index_bg_conf_distance,_TRUE_,index_bg,1);
  class_define_index(pba->index_bg_ang_distance,_TRUE_,index_bg,1);
  class_define_index(pba->index_bg_lum_distance,_TRUE_,index_bg,1);
  class_define_index(pba->index_bg_time,_TRUE_,index_bg,1);
  class_define_index(pba->index_bg_rs,_TRUE_,index_bg,1);
  class_define_index(pba->index_bg_D,_TRUE_,index_bg,1);
  class_define_index(pba->index_bg_f,_TRUE_,index_bg,1);
  class_define_index(pba->index_bg_varc_alpha,pba->has_varconst,index_bg,1);
  class_define_index(pba->index_bg_varc_me,pba->has_varconst,index_bg,1);

  pba->bg_size = index_bg;

  index_bi=0;

  class_define_index(pba->index_bi_tau,_TRUE_,index_bi,1);
  class_define_index(pba->index_bi_rho_dcdm,pba->has_dcdm,index_bi,1);
  class_define_index(pba->index_bi_rho_dr,pba->has_dr,index_bi,1);
  class_define_index(pba->index_bi_rho_fld,pba->has_fld,index_bi,1);
  class_define_index(pba->index_bi_phi_scf,pba->has_scf,index_bi,1);
  class_define_index(pba->index_bi_phi_prime_scf,pba->has_scf,index_bi,1);

  pba->bi_B_size = index_bi;

  class_define_index(pba->index_bi_time,_TRUE_,index_bi,1);
  class_define_index(pba->index_bi_rs,_TRUE_,index_bi,1);
  class_define_index(pba->index_bi_D,_TRUE_,index_bi,1);
  class_define_index(pba->index_bi_D_prime,_TRUE_,index_bi,1);

  pba->bi_size = index_bi;

  return _SUCCESS_;

}

/**
 * This is the routine where the distribution function f0(q) of each
 * ncdm species is specified.
 */

int background_ncdm_distribution(
                                 void * pbadist,
                                 double q,
                                 double * f0
                                 ) {
  struct background * pba;
  struct background_parameters_for_distributions * pbadist_local;
  int n_ncdm,lastidx;
  double ksi;
  double qlast,dqlast,f0last,df0last;
  double *param;

  pbadist_local = pbadist;
  pba = pbadist_local->pba;
  param = pba->ncdm_psd_parameters;
  n_ncdm = pbadist_local->n_ncdm;
  ksi = pba->ksi_ncdm[n_ncdm];

  if (pba->got_files[n_ncdm]==_TRUE_) {

    lastidx = pbadist_local->tablesize-1;
    if (q<pbadist_local->q[0]) {
      *f0 = pbadist_local->f0[0];
    }
    else if (q>pbadist_local->q[lastidx]) {
      qlast=pbadist_local->q[lastidx];
      f0last=pbadist_local->f0[lastidx];
      dqlast=qlast - pbadist_local->q[lastidx-1];
      df0last=f0last - pbadist_local->f0[lastidx-1];

      *f0 = f0last*exp(-(qlast-q)*df0last/f0last/dqlast);
    }
    else{
      class_call(array_interpolate_spline(
                                          pbadist_local->q,
                                          pbadist_local->tablesize,
                                          pbadist_local->f0,
                                          pbadist_local->d2f0,
                                          1,
                                          q,
                                          &pbadist_local->last_index,
                                          f0,
                                          1,
                                          pba->error_message),
                 pba->error_message,     pba->error_message);
    }
  }
  else{
    *f0 = 1.0/pow(2*_PI_,3)*(1./(exp(q-ksi)+1.) +1./(exp(q+ksi)+1.));
  }

  return _SUCCESS_;
}

/**
 * This function is only used for the purpose of finding optimal
 * quadrature weights.
 */

int background_ncdm_test_function(
                                  void * pbadist,
                                  double q,
                                  double * test
                                  ) {

  double c = 2.0/(3.0*_zeta3_);
  double d = 120.0/(7.0*pow(_PI_,4));
  double e = 2.0/(45.0*_zeta5_);

  *test = pow(2.0*_PI_,3)/6.0*(c*q*q-d*q*q*q-e*q*q*q*q);

  return _SUCCESS_;
}

/**
 * This function finds optimal quadrature weights for each ncdm species.
 */

int background_ncdm_init(
                         struct precision *ppr,
                         struct background *pba
                         ) {

  int index_q, k,tolexp,row,status,filenum;
  double f0m2,f0m1,f0,f0p1,f0p2,dq,q,df0dq,tmp1,tmp2;
  struct background_parameters_for_distributions pbadist;
  FILE *psdfile;

  pbadist.pba = pba;

  class_alloc(pba->q_ncdm, sizeof(double*)*pba->N_ncdm,pba->error_message);
  class_alloc(pba->w_ncdm, sizeof(double*)*pba->N_ncdm,pba->error_message);
  class_alloc(pba->q_ncdm_bg, sizeof(double*)*pba->N_ncdm,pba->error_message);
  class_alloc(pba->w_ncdm_bg, sizeof(double*)*pba->N_ncdm,pba->error_message);
  class_alloc(pba->dlnf0_dlnq_ncdm, sizeof(double*)*pba->N_ncdm,pba->error_message);

  class_alloc(pba->q_size_ncdm,sizeof(int)*pba->N_ncdm,pba->error_message);
  class_alloc(pba->q_size_ncdm_bg,sizeof(int)*pba->N_ncdm,pba->error_message);
  class_alloc(pba->factor_ncdm,sizeof(double)*pba->N_ncdm,pba->error_message);

  for (k=0, filenum=0; k<pba->N_ncdm; k++) {
    pbadist.n_ncdm = k;
    pbadist.q = NULL;
    pbadist.tablesize = 0;
    if ((pba->got_files!=NULL)&&(pba->got_files[k]==_TRUE_)) {
      psdfile = fopen(pba->ncdm_psd_files+filenum*_ARGUMENT_LENGTH_MAX_,"r");
      class_test(psdfile == NULL,pba->error_message,
                 "Could not open file %s!",pba->ncdm_psd_files+filenum*_ARGUMENT_LENGTH_MAX_);
      for (row=0,status=2; status==2; row++) {
        status = fscanf(psdfile,"%lf %lf",&tmp1,&tmp2);
      }
      rewind(psdfile);
      pbadist.tablesize = row-1;

      class_alloc(pbadist.q,sizeof(double)*pbadist.tablesize,pba->error_message);
      class_alloc(pbadist.f0,sizeof(double)*pbadist.tablesize,pba->error_message);
      class_alloc(pbadist.d2f0,sizeof(double)*pbadist.tablesize,pba->error_message);
      for (row=0; row<pbadist.tablesize; row++) {
        status = fscanf(psdfile,"%lf %lf",
                        &pbadist.q[row],&pbadist.f0[row]);
      }
      fclose(psdfile);
      class_call(array_spline_table_lines(pbadist.q,
                                          pbadist.tablesize,
                                          pbadist.f0,
                                          1,
                                          pbadist.d2f0,
                                          _SPLINE_EST_DERIV_,
                                          pba->error_message),
                 pba->error_message,
                 pba->error_message);
      filenum++;
    }

    if (pba->ncdm_quadrature_strategy[k]==qm_auto) {
      class_alloc(pba->q_ncdm[k],_QUADRATURE_MAX_*sizeof(double),pba->error_message);
      class_alloc(pba->w_ncdm[k],_QUADRATURE_MAX_*sizeof(double),pba->error_message);

      class_call(get_qsampling(pba->q_ncdm[k],
                               pba->w_ncdm[k],
                               &(pba->q_size_ncdm[k]),
                               _QUADRATURE_MAX_,
                               ppr->tol_ncdm,
                               pbadist.q,
                               pbadist.tablesize,
                               background_ncdm_test_function,
                               background_ncdm_distribution,
                               &pbadist,
                               pba->error_message),
                 pba->error_message,
                 pba->error_message);
      class_realloc(pba->q_ncdm[k],pba->q_size_ncdm[k]*sizeof(double), pba->error_message);
      class_realloc(pba->w_ncdm[k],pba->q_size_ncdm[k]*sizeof(double), pba->error_message);

      if (pba->background_verbose > 0) {
        printf("ncdm species i=%d sampled with %d points for purpose of perturbation integration\n",
               k+1,
               pba->q_size_ncdm[k]);
      }

      class_alloc(pba->q_ncdm_bg[k],_QUADRATURE_MAX_BG_*sizeof(double),pba->error_message);
      class_alloc(pba->w_ncdm_bg[k],_QUADRATURE_MAX_BG_*sizeof(double),pba->error_message);

      class_call(get_qsampling(pba->q_ncdm_bg[k],
                               pba->w_ncdm_bg[k],
                               &(pba->q_size_ncdm_bg[k]),
                               _QUADRATURE_MAX_BG_,
                               ppr->tol_ncdm_bg,
                               pbadist.q,
                               pbadist.tablesize,
                               background_ncdm_test_function,
                               background_ncdm_distribution,
                               &pbadist,
                               pba->error_message),
                 pba->error_message,
                 pba->error_message);

      class_realloc(pba->q_ncdm_bg[k],pba->q_size_ncdm_bg[k]*sizeof(double), pba->error_message);
      class_realloc(pba->w_ncdm_bg[k],pba->q_size_ncdm_bg[k]*sizeof(double), pba->error_message);

      if (pba->background_verbose > 0) {
        printf("ncdm species i=%d sampled with %d points for purpose of background integration\n",
               k+1,
               pba->q_size_ncdm_bg[k]);
      }
    }
    else{
      pba->q_size_ncdm_bg[k] = pba->ncdm_input_q_size[k];
      pba->q_size_ncdm[k] = pba->ncdm_input_q_size[k];
      class_alloc(pba->q_ncdm_bg[k],pba->q_size_ncdm_bg[k]*sizeof(double),pba->error_message);
      class_alloc(pba->w_ncdm_bg[k],pba->q_size_ncdm_bg[k]*sizeof(double),pba->error_message);
      class_alloc(pba->q_ncdm[k],pba->q_size_ncdm[k]*sizeof(double),pba->error_message);
      class_alloc(pba->w_ncdm[k],pba->q_size_ncdm[k]*sizeof(double),pba->error_message);
      class_call(get_qsampling_manual(pba->q_ncdm[k],
                                      pba->w_ncdm[k],
                                      pba->q_size_ncdm[k],
                                      pba->ncdm_qmax[k],
                                      pba->ncdm_quadrature_strategy[k],
                                      pbadist.q,
                                      pbadist.tablesize,
                                      background_ncdm_distribution,
                                      &pbadist,
                                      pba->error_message),
                 pba->error_message,
                 pba->error_message);
      for (index_q=0; index_q<pba->q_size_ncdm[k]; index_q++) {
        pba->q_ncdm_bg[k][index_q] = pba->q_ncdm[k][index_q];
        pba->w_ncdm_bg[k][index_q] = pba->w_ncdm[k][index_q];
      }
      if (pba->background_verbose > 0) {
        printf("ncdm species i=%d sampled with %d points for purpose of background and perturbation integration using the manual method\n",
               k+1,
               pba->q_size_ncdm[k]);
      }
    }

    class_alloc(pba->dlnf0_dlnq_ncdm[k],
                pba->q_size_ncdm[k]*sizeof(double),
                pba->error_message);


    for (index_q=0; index_q<pba->q_size_ncdm[k]; index_q++) {
      q = pba->q_ncdm[k][index_q];
      class_call(background_ncdm_distribution(&pbadist,q,&f0),
                 pba->error_message,pba->error_message);

      for (tolexp=_PSD_DERIVATIVE_EXP_MIN_; tolexp<_PSD_DERIVATIVE_EXP_MAX_; tolexp++) {

        if (index_q == 0) {
          dq = MIN((0.5-ppr->smallest_allowed_variation)*q,2*exp(tolexp)*(pba->q_ncdm[k][index_q+1]-q));
        }
        else if (index_q == pba->q_size_ncdm[k]-1) {
          dq = exp(tolexp)*2.0*(pba->q_ncdm[k][index_q]-pba->q_ncdm[k][index_q-1]);
        }
        else{
          dq = exp(tolexp)*(pba->q_ncdm[k][index_q+1]-pba->q_ncdm[k][index_q-1]);
        }

        class_call(background_ncdm_distribution(&pbadist,q-2*dq,&f0m2),
                   pba->error_message,pba->error_message);
        class_call(background_ncdm_distribution(&pbadist,q+2*dq,&f0p2),
                   pba->error_message,pba->error_message);

        if (fabs((f0p2-f0m2)/f0)>sqrt(ppr->smallest_allowed_variation)) break;
      }

      class_call(background_ncdm_distribution(&pbadist,q-dq,&f0m1),
                 pba->error_message,pba->error_message);
      class_call(background_ncdm_distribution(&pbadist,q+dq,&f0p1),
                 pba->error_message,pba->error_message);
      df0dq = (+f0m2-8*f0m1+8*f0p1-f0p2)/12.0/dq;
      if (fabs(f0)==0.)
        pba->dlnf0_dlnq_ncdm[k][index_q] = -q;
      else
        pba->dlnf0_dlnq_ncdm[k][index_q] = q/f0*df0dq;
    }

    pba->factor_ncdm[k]=pba->deg_ncdm[k]*4*_PI_*pow(pba->T_cmb*pba->T_ncdm[k]*_k_B_,4)*8*_PI_*_G_
      /3./pow(_h_P_/2./_PI_,3)/pow(_c_,7)*_Mpc_over_m_*_Mpc_over_m_;

    if ((pba->got_files!=NULL)&&(pba->got_files[k]==_TRUE_)) {
      free(pbadist.q);
      free(pbadist.f0);
      free(pbadist.d2f0);
    }
  }


  return _SUCCESS_;
}

/**
 * For a given ncdm species: given the quadrature weights, the mass
 * and the redshift, find background quantities by a quick weighted sum.
 */

int background_ncdm_momenta(
                            double * qvec,
                            double * wvec,
                            int qsize,
                            double M,
                            double factor,
                            double z,
                            double * n,
                            double * rho,
                            double * p,
                            double * drho_dM,
                            double * pseudo_p
                            ) {

  int index_q;
  double epsilon;
  double q2;
  double factor2;

  factor2 = factor*pow(1+z,4);

  if (n!=NULL) *n = 0.;
  if (rho!=NULL) *rho = 0.;
  if (p!=NULL) *p = 0.;
  if (drho_dM!=NULL) *drho_dM = 0.;
  if (pseudo_p!=NULL) *pseudo_p = 0.;

  for (index_q=0; index_q<qsize; index_q++) {

    q2 = qvec[index_q]*qvec[index_q];
    epsilon = sqrt(q2+M*M/(1.+z)/(1.+z));

    if (n!=NULL) *n += q2*wvec[index_q];
    if (rho!=NULL) *rho += q2*epsilon*wvec[index_q];
    if (p!=NULL) *p += q2*q2/3./epsilon*wvec[index_q];
    if (drho_dM!=NULL) *drho_dM += q2*M/(1.+z)/(1.+z)/epsilon*wvec[index_q];
    if (pseudo_p!=NULL) *pseudo_p += pow(q2/epsilon,3)/3.0*wvec[index_q];
  }

  if (n!=NULL) *n *= factor2/(1.+z);
  if (rho!=NULL) *rho *= factor2;
  if (p!=NULL) *p *= factor2;
  if (drho_dM!=NULL) *drho_dM *= factor2;
  if (pseudo_p!=NULL) *pseudo_p *=factor2;

  return _SUCCESS_;
}

/**
 * When the user passed the density fraction Omega_ncdm or
 * omega_ncdm in input but not the mass, infer the mass with Newton iteration method.
 */

int background_ncdm_M_from_Omega(
                                 struct precision *ppr,
                                 struct background *pba,
                                 int n_ncdm
                                 ) {
  double rho0,rho,n,M,deltaM,drhodM;
  int iter,maxiter=50;

  rho0 = pba->H0*pba->H0*pba->Omega0_ncdm[n_ncdm];
  M = 0.0;

  background_ncdm_momenta(pba->q_ncdm_bg[n_ncdm],
                          pba->w_ncdm_bg[n_ncdm],
                          pba->q_size_ncdm_bg[n_ncdm],
                          M,
                          pba->factor_ncdm[n_ncdm],
                          0.,
                          &n,
                          &rho,
                          NULL,
                          NULL,
                          NULL);

  class_test(rho0<rho,pba->error_message,
             "The value of Omega for the %dth species, %g, is less than for a massless species! It should be atleast %g. Check your input.",
             n_ncdm,pba->Omega0_ncdm[n_ncdm],pba->Omega0_ncdm[n_ncdm]*rho/rho0);

  M = rho0/n;
  for (iter=1; iter<=maxiter; iter++) {

    background_ncdm_momenta(pba->q_ncdm_bg[n_ncdm],
                            pba->w_ncdm_bg[n_ncdm],
                            pba->q_size_ncdm_bg[n_ncdm],
                            M,
                            pba->factor_ncdm[n_ncdm],
                            0.,
                            NULL,
                            &rho,
                            NULL,
                            &drhodM,
                            NULL);

    deltaM = (rho0-rho)/drhodM;
    if ((M+deltaM)<0.0) deltaM = -M/2.0;
    M += deltaM;
    if (fabs(deltaM/M)<ppr->tol_M_ncdm) {
      pba->M_ncdm[n_ncdm] = M;
      break;
    }
  }
  class_test(iter>=maxiter,pba->error_message,
             "Newton iteration could not converge on a mass for some reason.");
  return _SUCCESS_;
}

/**
 * Perform some check on the input background quantities.
 */

int background_checks(
                      struct precision* ppr,
                      struct background* pba
                      ) {

  int n_ncdm;
  double rho_ncdm_rel,rho_nu_rel;
  double N_dark;
  double w_fld, dw_over_da, integral_fld;
  int filenum=0;

  class_test((pba->Omega0_g<=0) || (pba->Omega0_b<=0),
             pba->error_message,
             "CLASS is conceived to work in a universe containing at least two species: photons and baryons. You could work in the limit where Omega_g or Omega_b are very small, but not zero");

  class_test(fabs(pba->h * 1.e5 / _c_  / pba->H0 -1.)>ppr->smallest_allowed_variation,
             pba->error_message,
             "inconsistency between Hubble and reduced Hubble parameters: you have H0=%f/Mpc=%fkm/s/Mpc, but h=%f",pba->H0,pba->H0/1.e5* _c_,pba->h);

  if (pba->has_fld == _TRUE_) {

    class_call(background_w_fld(pba,0.,&w_fld,&dw_over_da,&integral_fld), pba->error_message, pba->error_message);

    class_test(w_fld >= 1./3.,
               pba->error_message,
               "Your choice for w(a--->0)=%g is suspicious, since it is bigger than 1/3 there cannot be radiation domination at early times\n",
               w_fld);
  }

  if (pba->has_varconst == _TRUE_) {
    class_test(pba->varconst_alpha <= 0,
               pba->error_message,
               "incorrect fine structure constant before transition");
    class_test(pba->varconst_me <= 0,
               pba->error_message,
               "incorrect effective electron mass before transition");
    class_test(pba->varconst_transition_redshift < 0,
               pba->error_message,
               "incorrect transition redshift");
  }

  if (pba->background_verbose > 0) {

    if (pba->has_ncdm == _TRUE_) {

      for (n_ncdm=0;n_ncdm<pba->N_ncdm; n_ncdm++) {

        if (pba->got_files[n_ncdm] == _TRUE_) {
          printf(" -> ncdm species i=%d read from file %s\n",n_ncdm+1,pba->ncdm_psd_files+filenum*_ARGUMENT_LENGTH_MAX_);
          filenum++;
        }

        printf(" -> non-cold dark matter species with i=%d has m_i = %e eV (so m_i / omega_i =%e eV)\n",
               n_ncdm+1,
               pba->m_ncdm_in_eV[n_ncdm],
               pba->m_ncdm_in_eV[n_ncdm]*pba->deg_ncdm[n_ncdm]/pba->Omega0_ncdm[n_ncdm]/pba->h/pba->h);

        background_ncdm_momenta(pba->q_ncdm_bg[n_ncdm],
                                pba->w_ncdm_bg[n_ncdm],
                                pba->q_size_ncdm_bg[n_ncdm],
                                0.,
                                pba->factor_ncdm[n_ncdm],
                                0.,
                                NULL,
                                &rho_ncdm_rel,
                                NULL,
                                NULL,
                                NULL);

        rho_nu_rel = 56.0/45.0*pow(_PI_,6)*pow(4.0/11.0,4.0/3.0)*_G_/pow(_h_P_,3)/pow(_c_,7)*
          pow(_Mpc_over_m_,2)*pow(pba->T_cmb*_k_B_,4);

        printf(" -> ncdm species i=%d sampled with %d (resp. %d) points for purpose of background (resp. perturbation) integration. In the relativistic limit it gives Delta N_eff = %g\n",
               n_ncdm+1,
               pba->q_size_ncdm_bg[n_ncdm],
               pba->q_size_ncdm[n_ncdm],
               rho_ncdm_rel/rho_nu_rel);
      }
    }

    if (pba->has_idr == _TRUE_) {
      N_dark = pba->Omega0_idr/7.*8./pow(4./11.,4./3.)/pba->Omega0_g;
      printf(" -> dark radiation Delta Neff %e\n",N_dark);
    }
    
    /* HQIV-specific output */
    if (pba->hqiv.hqiv_on == _TRUE_) {
      printf(" -> HQIV cosmology enabled:\n");
      printf("    gamma_theory = %.4f (fixed from overlap integral; no free parameter)\n", pba->hqiv.gamma_hqiv);
      printf("    gamma_eff(a): 0 pre-recombination, smoothly -> %.2f post-recombination\n", pba->hqiv.gamma_hqiv);
      printf("    chi = %.4f, f_min = %.4f\n", pba->hqiv.chi_hqiv, pba->hqiv.fmin_hqiv);
      printf("    G_eff/G_0 = H(3H - gamma_eff*H0) / [H0(3H0 - gamma*H0)] (exact, no alpha)\n");
    }
  }

  return _SUCCESS_;
}

/**
 * Integrate the background over time, allocate and fill the background table.
 */

int background_solve(
                     struct precision *ppr,
                     struct background *pba
                     ) {

  /** HQIV emergent lattice mode: fully parameter-free 4D object from pre-computed table */
  if (pba->hqiv_emergent == _TRUE_) {
    class_call(hqiv_emergent_lattice_solve(ppr, pba),
               pba->error_message, pba->error_message);
    return _SUCCESS_;
  }

  struct background_parameters_and_workspace bpaw;
  double * pvecback_integration;
  double * pvecback;
  double comoving_radius=0.;
  double conformal_distance;

  extern int evolver_rk(EVOLVER_PROTOTYPE);
  extern int evolver_ndf15(EVOLVER_PROTOTYPE);
  int (*generic_evolver)(EVOLVER_PROTOTYPE) = evolver_ndf15;

  double loga_ini, loga_final;
  double D_today;
  int index_loga, index_scf;
  int * used_in_output;
  int h0_iter;
  double H_today_closure, omega_b_closure, omega_g_closure, omega_ur_closure, h_new_closure;

  int n_ncdm;

  bpaw.pba = pba;
  class_alloc(pvecback,pba->bg_size*sizeof(double),pba->error_message);
  bpaw.pvecback = pvecback;

  class_alloc(pvecback_integration,pba->bi_size*sizeof(double),pba->error_message);

  for (h0_iter = 0; h0_iter < 20; h0_iter++) {

  class_call(background_initial_conditions(ppr,pba,pvecback,pvecback_integration,&(loga_ini)),
             pba->error_message,
             pba->error_message);

  loga_final = 0.;
  pba->bt_size = ppr->background_Nloga;

  class_alloc(pba->tau_table,pba->bt_size * sizeof(double),pba->error_message);
  class_alloc(pba->z_table,pba->bt_size * sizeof(double),pba->error_message);
  class_alloc(pba->loga_table,pba->bt_size * sizeof(double),pba->error_message);

  class_alloc(pba->d2tau_dz2_table,pba->bt_size * sizeof(double),pba->error_message);
  class_alloc(pba->d2z_dtau2_table,pba->bt_size * sizeof(double),pba->error_message);

  class_alloc(pba->background_table,pba->bt_size * pba->bg_size * sizeof(double),pba->error_message);
  class_alloc(pba->d2background_dloga2_table,pba->bt_size * pba->bg_size * sizeof(double),pba->error_message);

  class_alloc(used_in_output, pba->bt_size*sizeof(int), pba->error_message);

  for (index_loga=0; index_loga<pba->bt_size; index_loga++) {
    pba->loga_table[index_loga] = loga_ini + index_loga*(loga_final-loga_ini)/(pba->bt_size-1);
    used_in_output[index_loga] = 1;
  }

  switch (ppr->background_evolver) {

  case rk:
    generic_evolver = evolver_rk;
    if (pba->background_verbose > 1) {
      printf("%s\n", "Chose rk as generic_evolver");
    }
    break;

  case ndf15:
    generic_evolver = evolver_ndf15;
    if (pba->background_verbose > 1) {
      printf("%s\n", "Chose ndf15 as generic_evolver");
    }
    break;
  }

  class_call(generic_evolver(background_derivs,
                             loga_ini,
                             loga_final,
                             pvecback_integration,
                             used_in_output,
                             pba->bi_size,
                             &bpaw,
                             ppr->tol_background_integration,
                             ppr->smallest_allowed_variation,
                             background_timescale,
                             ppr->background_integration_stepsize,
                             pba->loga_table,
                             pba->bt_size,
                             background_sources,
                             NULL,
                             pba->error_message),
             pba->error_message,
             pba->error_message);

  pba->age = pvecback_integration[pba->index_bi_time]/_Gyr_over_Mpc_;
  pba->conformal_age = pvecback_integration[pba->index_bi_tau];
  if (pba->has_dcdm == _TRUE_) {
    pba->Omega0_dcdm = pvecback_integration[pba->index_bi_rho_dcdm]/pba->H0/pba->H0;
  }
  if (pba->has_dr == _TRUE_) {
    pba->Omega0_dr = pvecback_integration[pba->index_bi_rho_dr]/pba->H0/pba->H0;
  }
  D_today = pvecback_integration[pba->index_bi_D];

  for (index_loga=0; index_loga < pba->bt_size; index_loga++) {

    pba->background_table[index_loga*pba->bg_size+pba->index_bg_D]*= 1./D_today;

    conformal_distance = pba->conformal_age - pba->tau_table[index_loga];
    pba->background_table[index_loga*pba->bg_size+pba->index_bg_conf_distance] = conformal_distance;

    if (pba->sgnK == 0) { comoving_radius = conformal_distance; }
    else if (pba->sgnK == 1) { comoving_radius = sin(sqrt(pba->K)*conformal_distance)/sqrt(pba->K); }
    else if (pba->sgnK == -1) { comoving_radius = sinh(sqrt(-pba->K)*conformal_distance)/sqrt(-pba->K); }

    pba->background_table[index_loga*pba->bg_size+pba->index_bg_ang_distance] = comoving_radius/(1.+pba->z_table[index_loga]);
    pba->background_table[index_loga*pba->bg_size+pba->index_bg_lum_distance] = comoving_radius*(1.+pba->z_table[index_loga]);
  }

  class_call(array_spline_table_lines(pba->z_table,
                                      pba->bt_size,
                                      pba->tau_table,
                                      1,
                                      pba->d2tau_dz2_table,
                                      _SPLINE_EST_DERIV_,
                                      pba->error_message),
             pba->error_message,
             pba->error_message);

  class_call(array_spline_table_lines(pba->tau_table,
                                      pba->bt_size,
                                      pba->z_table,
                                      1,
                                      pba->d2z_dtau2_table,
                                      _SPLINE_EST_DERIV_,
                                      pba->error_message),
             pba->error_message,
             pba->error_message);

  class_call(array_spline_table_lines(pba->loga_table,
                                      pba->bt_size,
                                      pba->background_table,
                                      pba->bg_size,
                                      pba->d2background_dloga2_table,
                                      _SPLINE_EST_DERIV_,
                                      pba->error_message),
             pba->error_message,
             pba->error_message);

  pba->Neff = (pba->background_table[pba->index_bg_Omega_r]
               *pba->background_table[pba->index_bg_rho_crit]
               -pba->background_table[pba->index_bg_rho_g])
    /(7./8.*pow(4./11.,4./3.)*pba->background_table[pba->index_bg_rho_g]);

  if (pba->background_verbose > 0) {
    printf(" -> age = %f Gyr\n",pba->age);
    printf(" -> conformal age = %f Mpc\n",pba->conformal_age);
    printf(" -> N_eff = %g (summed over all species that are non-relativistic at early times) \n",pba->Neff);
  }

  pba->Omega0_m = pba->background_table[(pba->bt_size-1)*pba->bg_size+pba->index_bg_Omega_m];
  pba->Omega0_r = pba->background_table[(pba->bt_size-1)*pba->bg_size+pba->index_bg_Omega_r];

  if (pba->hqiv.hqiv_on == _TRUE_) {
    /* In HQIV, the "missing energy" is the horizon term γH, not Λ.
     * From H² = ρ + (γ_class/3)H → Ω_m + Ω_r + Ω_horizon = 1
     * where Ω_horizon = γ_class/(3 H_today). */
    double H_today = pba->background_table[(pba->bt_size-1)*pba->bg_size+pba->index_bg_H];
    double g_class = pba->hqiv.gamma_hqiv * pba->H0;
    double Omega_hz = (H_today > 0.0) ? g_class / (3.0 * H_today) : 0.0;
    pba->Omega0_de = Omega_hz;
    if (pba->background_verbose > 0) {

      /* --- HQIV energy budget --- */
      double h_ratio = H_today / pba->H0;
      double h_actual_kms = H_today * 299792.458;  /* H in km/s/Mpc */

      /* Omega_k physical: -K/H_actual² (Omega0_k uses H_input normalization) */
      double Omega_k_phys  = (H_today > 0.0) ? (-pba->K / (H_today * H_today)) : 0.0;

      printf(" -> HQIV energy budget at z=0:\n");
      printf("    H_input (h=%.3f)     = %.6e Mpc^-1  (%.1f km/s/Mpc)\n",
             pba->h, pba->H0, pba->h * 100.);
      printf("    H_actual(z=0)        = %.6e Mpc^-1  (%.1f km/s/Mpc, ratio = %.4f)\n",
             H_today, h_actual_kms, h_ratio);
      printf("    Omega_m              = %.6f\n", pba->Omega0_m);
      printf("    Omega_r              = %.6e\n", pba->Omega0_r);
      printf("    Omega_k (phys)       = %.6f  (input Omega_k = %.4f)\n", Omega_k_phys, pba->Omega0_k);
      printf("    Omega_horizon (γH)   = %.6f  (replaces dark energy)\n", Omega_hz);
      printf("    Sum (closure)        = %.6f\n",
             pba->Omega0_m + pba->Omega0_r + Omega_k_phys + Omega_hz);

      /* --- Energy conservation: CMB blackbody vs horizon + post-recomb energy ---
       *
       * At recombination (γ_eff = 0), all energy is in matter + radiation.
       * Today, the same energy budget is: matter + radiation + horizon energy.
       *
       * In a comoving volume:
       *   E_matter is conserved (ρ_m a³ = const)
       *   E_radiation loses energy to expansion (ρ_r a³ ∝ a⁻¹)
       *   E_horizon = (γ_class H / 3) a³ has GROWN from 0 at recombination
       *
       * The radiation energy "lost" between recombination and today equals:
       *   ΔE_rad = ρ_r(a_rec) a_rec³ - ρ_r(a_0) a_0³
       *          = ρ_r(a_0) a_0³ × (a_0/a_rec - 1)
       *          = ρ_r(a_0) × z_rec
       *
       * The horizon energy today (per comoving volume, a_0=1):
       *   E_horizon = γ_class × H_today / 3
       *
       * These should be related: the horizon "collected" the energy that
       * radiation lost to expansion.  Ratio = E_horizon / ΔE_rad.
       */
      {
        double rho_g_0 = pba->background_table[(pba->bt_size-1)*pba->bg_size+pba->index_bg_rho_g];
        double rho_ur_0 = pba->has_ur ? (pba->Omega0_ur * pba->H0 * pba->H0) : 0.0;
        double rho_r_0 = rho_g_0 + rho_ur_0;
        double z_rec = 1090.0;
        double a_rec = 1.0 / (1.0 + z_rec);

        /* Radiation energy lost since recombination (per comoving vol, a_0=1) */
        double rho_r_rec_comov = rho_r_0 * (1.0 + z_rec);   /* ρ_r(a_rec) × a_rec³ / a_0³ = ρ_r_0 × (1+z) */
        double delta_E_rad = rho_r_rec_comov - rho_r_0;      /* = ρ_r_0 × z_rec */

        /* Horizon energy density today */
        double rho_horizon_0 = g_class * H_today / 3.0;

        /* Total energy at recombination (comoving, no horizon term since γ_eff=0) */
        double rho_b_0 = pba->background_table[(pba->bt_size-1)*pba->bg_size+pba->index_bg_rho_b];
        double rho_m_rec_comov = rho_b_0;  /* matter conserved: ρ_m a³ = const */
        double E_rec = rho_m_rec_comov + rho_r_rec_comov;

        /* Total energy today = matter + radiation + horizon */
        double E_today = rho_b_0 + rho_r_0 + rho_horizon_0;

        printf("\n -> HQIV energy conservation (recombination → today, comoving volume):\n");
        printf("    ρ_CMB(z=0)           = %.6e  (T₀ = 2.725 K blackbody)\n", rho_g_0);
        printf("    ρ_radiation(z=0)     = %.6e  (photons + neutrinos)\n", rho_r_0);
        printf("    ρ_horizon(z=0)       = %.6e  (γH/3 accumulated since z_rec)\n", rho_horizon_0);
        printf("    ΔE_rad lost(z_rec→0) = %.6e  (radiation redshift loss)\n", delta_E_rad);
        printf("    ρ_horizon / ΔE_rad   = %.4f  (horizon captured fraction)\n",
               delta_E_rad > 0 ? rho_horizon_0 / delta_E_rad : 0.0);
        printf("    E(z_rec, comov)      = %.6e  (matter + radiation, no horizon)\n", E_rec);
        printf("    E(z=0, comov)        = %.6e  (matter + radiation + horizon)\n", E_today);
        printf("    E(z=0) / E(z_rec)    = %.6f\n", E_rec > 0 ? E_today / E_rec : 0.0);
      }

      /* --- Time dilation: wall-clock vs observed ---
       *
       * Wall-clock age: from integrating dt = da/(aH) with the actual H(a).
       * Observers measure H_actual(z=0), not H_input.
       *
       * Hubble time:  t_H = 1/H_actual  (upper bound on age)
       * ΛCDM inferred age: t_LCDM ≈ (2/3)/H_0_observed × correction ≈ 13.8 Gyr
       *   (for H_0 = 73.4 km/s/Mpc, Planck ΛCDM gives ~13.8 Gyr)
       *
       * The time dilation factor shows how much older the universe really is
       * compared to what a ΛCDM observer would infer from the same local H.
       */
      {
        double t_wall = pba->age;  /* Gyr */
        double t_Hubble_actual = 1.0 / H_today / _Gyr_over_Mpc_;  /* Gyr */
        double t_Hubble_input = 1.0 / pba->H0 / _Gyr_over_Mpc_;   /* Gyr */

        /* ΛCDM age for Ω_m=0.3, Ω_Λ=0.7:
         * t = (1/H_0) × (2/3√Ω_Λ) × arcsinh(√(Ω_Λ/Ω_m))
         *   = 0.9643 / H_0  (for standard Planck-like Ω_m=0.3) */
        double t_LCDM = 0.9643 / pba->H0 / _Gyr_over_Mpc_;

        /* What a ΛCDM observer measuring H_actual would infer */
        double t_LCDM_from_Hactual = 0.9643 / H_today / _Gyr_over_Mpc_;

        printf("\n -> HQIV time dilation:\n");
        printf("    Wall-clock age                = %.2f Gyr  (actual, integrated)\n", t_wall);
        printf("    1/H_actual(z=0)               = %.2f Gyr  (Hubble time)\n", t_Hubble_actual);
        printf("    1/H_input                     = %.2f Gyr  (Hubble time, input h=%.3f)\n",
               t_Hubble_input, pba->h);
        printf("    ΛCDM age (H_input, h=%.3f)    = %.2f Gyr  (ΛCDM observer, Ω_m=0.3)\n",
               pba->h, t_LCDM);
        printf("    ΛCDM age (H_actual=%.1f km/s)  = %.2f Gyr  (ΛCDM at actual H)\n",
               h_actual_kms, t_LCDM_from_Hactual);
        printf("    Dilation: wall / ΛCDM(H_input)= %.2f×  (universe %.1f× older than ΛCDM says)\n",
               t_LCDM > 0 ? t_wall / t_LCDM : 0.0, t_LCDM > 0 ? t_wall / t_LCDM : 0.0);
        printf("    Geodesic compression: H_actual/H_input = %.4f\n", h_ratio);
      }
    }
    /* H0 closure: shoot for h so H(a=1) = H0 (fully dynamic H) */
    if (pba->hqiv.H0_closure == _TRUE_) {
      H_today_closure = pba->background_table[(pba->bt_size-1)*pba->bg_size+pba->index_bg_H];
      if (fabs(H_today_closure - pba->H0) < 1e-6 * pba->H0)
        break;
      omega_b_closure = pba->Omega0_b * pba->h * pba->h;
      omega_g_closure = pba->Omega0_g * pba->h * pba->h;
      omega_ur_closure = pba->Omega0_ur * pba->h * pba->h;
      h_new_closure = pba->h * (H_today_closure / pba->H0);
      pba->h = h_new_closure;
      pba->H0 = h_new_closure * 1.e5/_c_/_Mpc_over_m_;
      pba->Omega0_b = omega_b_closure / (h_new_closure * h_new_closure);
      pba->Omega0_g = omega_g_closure / (h_new_closure * h_new_closure);
      pba->Omega0_ur = omega_ur_closure / (h_new_closure * h_new_closure);
      pba->K = -pba->Omega0_k * pba->H0 * pba->H0;
    } else {
      break;
    }
  } else {
    break;
  }
  }  /* end for h0_iter (H0_closure) */

  if (pba->hqiv.hqiv_on != _TRUE_) {
    pba->Omega0_de = 1. - (pba->Omega0_m + pba->Omega0_r + pba->Omega0_k);
  }

  pba->Omega0_nfsm =  pba->Omega0_b;
  if (pba->has_cdm == _TRUE_)
    pba->Omega0_nfsm += pba->Omega0_cdm;
  if (pba->has_idm == _TRUE_)
    pba->Omega0_nfsm += pba->Omega0_idm;
  if (pba->has_dcdm == _TRUE_)
    pba->Omega0_nfsm += pba->Omega0_dcdm;
  for (n_ncdm=0;n_ncdm<pba->N_ncdm; n_ncdm++) {
    if (pba->M_ncdm[n_ncdm] > ppr->M_nfsm_threshold) {
      pba->Omega0_nfsm += pba->Omega0_ncdm[n_ncdm];
    }
  }

  free(pvecback);
  free(pvecback_integration);
  free(used_in_output);

  return _SUCCESS_;

}

/**
 * Assign initial values to background integrated variables.
 */

int background_initial_conditions(
                                  struct precision *ppr,
                                  struct background *pba,
                                  double * pvecback,
                                  double * pvecback_integration,
                                  double * loga_ini
                                  ) {

  double a;

  double rho_ncdm, p_ncdm, rho_ncdm_rel_tot=0.;
  double f,Omega_rad, rho_rad;
  int counter,is_early_enough,n_ncdm;
  double scf_lambda;
  double rho_fld_today;
  double w_fld,dw_over_da_fld,integral_fld;

  a = ppr->a_ini_over_a_today_default;

  if (pba->has_ncdm == _TRUE_) {

    for (counter=0; counter < _MAX_IT_; counter++) {

      is_early_enough = _TRUE_;
      rho_ncdm_rel_tot = 0.;

      for (n_ncdm=0; n_ncdm<pba->N_ncdm; n_ncdm++) {

        class_call(background_ncdm_momenta(pba->q_ncdm_bg[n_ncdm],
                                           pba->w_ncdm_bg[n_ncdm],
                                           pba->q_size_ncdm_bg[n_ncdm],
                                           pba->M_ncdm[n_ncdm],
                                           pba->factor_ncdm[n_ncdm],
                                           1./a-1.0,
                                           NULL,
                                           &rho_ncdm,
                                           &p_ncdm,
                                           NULL,
                                           NULL),
                   pba->error_message,
                   pba->error_message);
        rho_ncdm_rel_tot += 3.*p_ncdm;
        if (fabs(p_ncdm/rho_ncdm-1./3.)>ppr->tol_ncdm_initial_w) {
          is_early_enough = _FALSE_;
        }
      }
      if (is_early_enough == _TRUE_) {
        break;
      }
      else {
        a *= _SCALE_BACK_;
      }
    }
    class_test(counter == _MAX_IT_,
               pba->error_message,
               "Search for initial scale factor a such that all ncdm species are relativistic failed.");
  }

  Omega_rad = pba->Omega0_g;
  if (pba->has_ur == _TRUE_) {
    Omega_rad += pba->Omega0_ur;
  }
  if (pba->has_idr == _TRUE_) {
    Omega_rad += pba->Omega0_idr;
  }
  rho_rad = Omega_rad*pow(pba->H0,2)/pow(a,4);
  if (pba->has_ncdm == _TRUE_) {
    rho_rad += rho_ncdm_rel_tot;
  }
  if (pba->has_dcdm == _TRUE_) {
    pvecback_integration[pba->index_bi_rho_dcdm] =
      pba->Omega_ini_dcdm*pba->H0*pba->H0*pow(a,-3);
    if (pba->background_verbose > 3)
      printf("Density is %g. Omega_ini=%g\n",pvecback_integration[pba->index_bi_rho_dcdm],pba->Omega_ini_dcdm);
  }

  if (pba->has_dr == _TRUE_) {
    if (pba->has_dcdm == _TRUE_) {
      f = 1./3.*pow(a,6)*pvecback_integration[pba->index_bi_rho_dcdm]*pba->Gamma_dcdm/pow(pba->H0,3)/sqrt(Omega_rad);
      pvecback_integration[pba->index_bi_rho_dr] = f*pba->H0*pba->H0/pow(a,4);
    }
    else{
      pvecback_integration[pba->index_bi_rho_dr] = 0.0;
    }
  }

  if (pba->has_fld == _TRUE_) {

    rho_fld_today = pba->Omega0_fld * pow(pba->H0,2);

    class_call(background_w_fld(pba,a,&w_fld,&dw_over_da_fld,&integral_fld), pba->error_message, pba->error_message);

    pvecback_integration[pba->index_bi_rho_fld] = rho_fld_today * exp(integral_fld);

  }

  if (pba->has_scf == _TRUE_) {
    scf_lambda = pba->scf_parameters[0];
    if (pba->attractor_ic_scf == _TRUE_) {
      pvecback_integration[pba->index_bi_phi_scf] = -1/scf_lambda*
        log(rho_rad*4./(3*pow(scf_lambda,2)-12))*pba->phi_ini_scf;
      if (3.*pow(scf_lambda,2)-12. < 0) {
        pvecback_integration[pba->index_bi_phi_scf] = 1./scf_lambda;
        if (pba->background_verbose > 0) {
          printf(" No attractor IC for lambda = %.3e ! \n ",scf_lambda);
        }
      }
      pvecback_integration[pba->index_bi_phi_prime_scf] = 2.*a*sqrt(V_scf(pba,pvecback_integration[pba->index_bi_phi_scf]))*pba->phi_prime_ini_scf;
    }
    else {
      printf("Not using attractor initial conditions\n");
      pvecback_integration[pba->index_bi_phi_scf] = pba->phi_ini_scf;
      pvecback_integration[pba->index_bi_phi_prime_scf] = pba->phi_prime_ini_scf;
    }
    class_test(!isfinite(pvecback_integration[pba->index_bi_phi_scf]) ||
               !isfinite(pvecback_integration[pba->index_bi_phi_scf]),
               pba->error_message,
               "initial phi = %e phi_prime = %e -> check initial conditions",
               pvecback_integration[pba->index_bi_phi_scf],
               pvecback_integration[pba->index_bi_phi_scf]);
  }

  class_call(background_functions(pba, a, pvecback_integration, normal_info, pvecback),
             pba->error_message,
             pba->error_message);

  class_test(fabs(pvecback[pba->index_bg_Omega_r]-1.) > ppr->tol_initial_Omega_r,
             pba->error_message,
             "Omega_r = %e, not close enough to 1. Decrease a_ini_over_a_today_default in order to start from radiation domination.",
             pvecback[pba->index_bg_Omega_r]);

  class_test(pvecback[pba->index_bg_H] <= 0.,
             pba->error_message,
             "H = %e instead of strictly positive",pvecback[pba->index_bg_H]);

  pvecback_integration[pba->index_bi_time] = 1./(2.* pvecback[pba->index_bg_H]);

  pvecback_integration[pba->index_bi_tau] = 1./(a * pvecback[pba->index_bg_H]);

  pvecback_integration[pba->index_bi_rs] = pvecback_integration[pba->index_bi_tau]/sqrt(3.);

  pvecback_integration[pba->index_bi_D] = 1.;
  pvecback_integration[pba->index_bi_D_prime] = 2.*a*pvecback[pba->index_bg_H];

  *loga_ini = log(a);

  return _SUCCESS_;

}

/**
 * Find the time of radiation/matter equality.
 */

int background_find_equality(
                             struct precision *ppr,
                             struct background *pba
                             ) {

  double Omega_m_over_Omega_r=0.;
  int index_tau_minus = 0;
  int index_tau_plus = pba->bt_size-1;
  int index_tau_mid = 0;
  double tau_minus,tau_plus,tau_mid=0.;
  double * pvecback;

  while ((index_tau_plus - index_tau_minus) > 1) {

    index_tau_mid = (int)(0.5*(index_tau_plus+index_tau_minus));

    Omega_m_over_Omega_r = pba->background_table[index_tau_mid*pba->bg_size+pba->index_bg_Omega_m]
      /pba->background_table[index_tau_mid*pba->bg_size+pba->index_bg_Omega_r];

    if (Omega_m_over_Omega_r > 1)
      index_tau_plus = index_tau_mid;
    else
      index_tau_minus = index_tau_mid;

  }

  tau_minus = pba->tau_table[index_tau_minus];
  tau_plus =  pba->tau_table[index_tau_plus];

  class_alloc(pvecback,pba->bg_size*sizeof(double),pba->error_message);

  while ((tau_plus - tau_minus) > ppr->tol_tau_eq) {

    tau_mid = 0.5*(tau_plus+tau_minus);

    class_call(background_at_tau(pba,tau_mid,long_info,inter_closeby,&index_tau_minus,pvecback),
               pba->error_message,
               pba->error_message);

    Omega_m_over_Omega_r = pvecback[pba->index_bg_Omega_m]/pvecback[pba->index_bg_Omega_r];

    if (Omega_m_over_Omega_r > 1)
      tau_plus = tau_mid;
    else
      tau_minus = tau_mid;

  }

  pba->a_eq = pvecback[pba->index_bg_a];
  pba->H_eq = pvecback[pba->index_bg_H];
  pba->z_eq = 1./pba->a_eq -1.;
  pba->tau_eq = tau_mid;

  if (pba->background_verbose > 0) {
    printf(" -> radiation/matter equality at z = %f\n",pba->z_eq);
    printf("    corresponding to conformal time = %f Mpc\n",pba->tau_eq);
  }

  free(pvecback);

  return _SUCCESS_;

}


/**
 * Subroutine for formatting background output
 */

int background_output_titles(
                             struct background * pba,
                             char titles[_MAXTITLESTRINGLENGTH_]
                             ) {

  int n;
  char tmp[40];

  class_store_columntitle(titles,"z",_TRUE_);
  class_store_columntitle(titles,"proper time [Gyr]",_TRUE_);
  class_store_columntitle(titles,"conf. time [Mpc]",_TRUE_);
  class_store_columntitle(titles,"H [1/Mpc]",_TRUE_);
  class_store_columntitle(titles,"comov. dist.",_TRUE_);
  class_store_columntitle(titles,"ang.diam.dist.",_TRUE_);
  class_store_columntitle(titles,"lum. dist.",_TRUE_);
  class_store_columntitle(titles,"comov.snd.hrz.",_TRUE_);
  class_store_columntitle(titles,"(.)rho_g",_TRUE_);
  class_store_columntitle(titles,"(.)rho_b",_TRUE_);
  class_store_columntitle(titles,"(.)rho_cdm",pba->has_cdm);
  class_store_columntitle(titles,"(.)rho_idm",pba->has_idm);
  if (pba->has_ncdm == _TRUE_) {
    for (n=0; n<pba->N_ncdm; n++) {
      class_sprintf(tmp,"(.)rho_ncdm[%d]",n);
      class_store_columntitle(titles,tmp,_TRUE_);
      class_sprintf(tmp,"(.)p_ncdm[%d]",n);
      class_store_columntitle(titles,tmp,_TRUE_);
    }
  }
  class_store_columntitle(titles,"(.)rho_lambda",pba->has_lambda);
  class_store_columntitle(titles,"(.)rho_fld",pba->has_fld);
  class_store_columntitle(titles,"(.)w_fld",pba->has_fld);
  class_store_columntitle(titles,"(.)rho_ur",pba->has_ur);
  class_store_columntitle(titles,"(.)rho_idr",pba->has_idr);
  class_store_columntitle(titles,"(.)rho_crit",_TRUE_);
  class_store_columntitle(titles,"(.)rho_dcdm",pba->has_dcdm);
  class_store_columntitle(titles,"(.)rho_dr",pba->has_dr);

  class_store_columntitle(titles,"(.)rho_scf",pba->has_scf);
  class_store_columntitle(titles,"(.)p_scf",pba->has_scf);
  class_store_columntitle(titles,"(.)p_prime_scf",pba->has_scf);
  class_store_columntitle(titles,"phi_scf",pba->has_scf);
  class_store_columntitle(titles,"phi'_scf",pba->has_scf);
  class_store_columntitle(titles,"V_scf",pba->has_scf);
  class_store_columntitle(titles,"V'_scf",pba->has_scf);
  class_store_columntitle(titles,"V''_scf",pba->has_scf);

  class_store_columntitle(titles,"(.)rho_tot",_TRUE_);
  class_store_columntitle(titles,"(.)p_tot",_TRUE_);
  class_store_columntitle(titles,"(.)p_tot_prime",_TRUE_);

  class_store_columntitle(titles,"Omega_r(z)",_TRUE_);
  class_store_columntitle(titles,"Omega_m(z)",_TRUE_);

  class_store_columntitle(titles,"gr.fac. D",_TRUE_);
  class_store_columntitle(titles,"gr.fac. f",_TRUE_);

  class_store_columntitle(titles,"rel. alpha",pba->has_varconst);
  class_store_columntitle(titles,"rel. m_e",pba->has_varconst);

  /* HQIV: horizon energy density and derived quantities */
  class_store_columntitle(titles,"gamma_eff",pba->hqiv.hqiv_on);
  class_store_columntitle(titles,"Omega_horizon",pba->hqiv.hqiv_on);
  class_store_columntitle(titles,"G_eff/G0",pba->hqiv.hqiv_on);

  return _SUCCESS_;
}

/**
 * Subroutine for writing the background output
 */

int background_output_data(
                           struct background *pba,
                           int number_of_titles,
                           double *data
                           ) {

  int index_tau, storeidx, n;
  double *dataptr, *pvecback;

  for (index_tau=0; index_tau<pba->bt_size; index_tau++) {
    dataptr = data + index_tau*number_of_titles;
    pvecback = pba->background_table + index_tau*pba->bg_size;
    storeidx = 0;

    class_store_double(dataptr,1./pvecback[pba->index_bg_a]-1.,_TRUE_,storeidx);
    class_store_double(dataptr,pvecback[pba->index_bg_time]/_Gyr_over_Mpc_,_TRUE_,storeidx);
    class_store_double(dataptr,pba->conformal_age-pvecback[pba->index_bg_conf_distance],_TRUE_,storeidx);
    class_store_double(dataptr,pvecback[pba->index_bg_H],_TRUE_,storeidx);
    class_store_double(dataptr,pvecback[pba->index_bg_conf_distance],_TRUE_,storeidx);
    class_store_double(dataptr,pvecback[pba->index_bg_ang_distance],_TRUE_,storeidx);
    class_store_double(dataptr,pvecback[pba->index_bg_lum_distance],_TRUE_,storeidx);
    class_store_double(dataptr,pvecback[pba->index_bg_rs],_TRUE_,storeidx);
    class_store_double(dataptr,pvecback[pba->index_bg_rho_g],_TRUE_,storeidx);
    class_store_double(dataptr,pvecback[pba->index_bg_rho_b],_TRUE_,storeidx);
    class_store_double(dataptr,pvecback[pba->index_bg_rho_cdm],pba->has_cdm,storeidx);
    class_store_double(dataptr,pvecback[pba->index_bg_rho_idm],pba->has_idm,storeidx);
    if (pba->has_ncdm == _TRUE_) {
      for (n=0; n<pba->N_ncdm; n++) {
        class_store_double(dataptr,pvecback[pba->index_bg_rho_ncdm1+n],_TRUE_,storeidx);
        class_store_double(dataptr,pvecback[pba->index_bg_p_ncdm1+n],_TRUE_,storeidx);
      }
    }
    class_store_double(dataptr,pvecback[pba->index_bg_rho_lambda],pba->has_lambda,storeidx);
    class_store_double(dataptr,pvecback[pba->index_bg_rho_fld],pba->has_fld,storeidx);
    class_store_double(dataptr,pvecback[pba->index_bg_w_fld],pba->has_fld,storeidx);
    class_store_double(dataptr,pvecback[pba->index_bg_rho_ur],pba->has_ur,storeidx);
    class_store_double(dataptr,pvecback[pba->index_bg_rho_idr],pba->has_idr,storeidx);
    class_store_double(dataptr,pvecback[pba->index_bg_rho_crit],_TRUE_,storeidx);
    class_store_double(dataptr,pvecback[pba->index_bg_rho_dcdm],pba->has_dcdm,storeidx);
    class_store_double(dataptr,pvecback[pba->index_bg_rho_dr],pba->has_dr,storeidx);

    class_store_double(dataptr,pvecback[pba->index_bg_rho_scf],pba->has_scf,storeidx);
    class_store_double(dataptr,pvecback[pba->index_bg_p_scf],pba->has_scf,storeidx);
    class_store_double(dataptr,pvecback[pba->index_bg_p_prime_scf],pba->has_scf,storeidx);
    class_store_double(dataptr,pvecback[pba->index_bg_phi_scf],pba->has_scf,storeidx);
    class_store_double(dataptr,pvecback[pba->index_bg_phi_prime_scf],pba->has_scf,storeidx);
    class_store_double(dataptr,pvecback[pba->index_bg_V_scf],pba->has_scf,storeidx);
    class_store_double(dataptr,pvecback[pba->index_bg_dV_scf],pba->has_scf,storeidx);
    class_store_double(dataptr,pvecback[pba->index_bg_ddV_scf],pba->has_scf,storeidx);

    class_store_double(dataptr,pvecback[pba->index_bg_rho_tot],_TRUE_,storeidx);
    class_store_double(dataptr,pvecback[pba->index_bg_p_tot],_TRUE_,storeidx);
    class_store_double(dataptr,pvecback[pba->index_bg_p_tot_prime],_TRUE_,storeidx);

    class_store_double(dataptr,pvecback[pba->index_bg_Omega_r],_TRUE_,storeidx);
    class_store_double(dataptr,pvecback[pba->index_bg_Omega_m],_TRUE_,storeidx);

    class_store_double(dataptr,pvecback[pba->index_bg_D],_TRUE_,storeidx);
    class_store_double(dataptr,pvecback[pba->index_bg_f],_TRUE_,storeidx);

    class_store_double(dataptr,pvecback[pba->index_bg_varc_alpha],pba->has_varconst,storeidx);
    class_store_double(dataptr,pvecback[pba->index_bg_varc_me],pba->has_varconst,storeidx);

    /* HQIV: horizon energy density columns */
    if (pba->hqiv.hqiv_on == _TRUE_) {
      double a_here = pvecback[pba->index_bg_a];
      double H_here = pvecback[pba->index_bg_H];
      double a_rec = 1.0 / 1090.0;
      double da_w  = 0.02 * a_rec;
      double g_eff = pba->hqiv.gamma_hqiv * 0.5 * (1.0 + tanh((a_here - a_rec) / da_w));
      double g_class = g_eff * pba->H0;
      /* Ω_horizon = γ_eff_class H / (3 H²) = γ_eff_class / (3 H)
       * Fractional energy from horizon term: 3H² - γH = 3ρ  → Ω_hz = γH/(3H²) */
      double Omega_hz = (H_here > 0.0) ? g_class / (3.0 * H_here) : 0.0;
      /* G_eff/G_0 from exact Friedmann relation with TODAY's γ_eff(z=0) = γ_theory
       * (denominator always uses the z=0 value for consistent normalization) */
      double g_class_today = pba->hqiv.gamma_hqiv * pba->H0;
      double denom_G = pba->H0 * (3.0 * pba->H0 - g_class_today);
      double G_ratio = 1.0;
      if (denom_G > 0.0) {
        double num_G = H_here * (3.0 * H_here - g_class);
        G_ratio = (num_G > 0.0) ? num_G / denom_G : 0.0;
      }
      class_store_double(dataptr, g_eff, _TRUE_, storeidx);
      class_store_double(dataptr, Omega_hz, _TRUE_, storeidx);
      class_store_double(dataptr, G_ratio, _TRUE_, storeidx);
    }
  }

  return _SUCCESS_;
}


/**
 * Subroutine evaluating the derivative with respect to loga
 * of quantities which are integrated (tau, t, etc).
 */

int background_derivs(
                      double loga,
                      double* y,
                      double* dy,
                      void * parameters_and_workspace,
                      ErrorMsg error_message
                      ) {

  struct background_parameters_and_workspace * pbpaw;
  struct background * pba;
  double * pvecback, a, H, rho_M;

  pbpaw = parameters_and_workspace;
  pba =  pbpaw->pba;
  pvecback = pbpaw->pvecback;

  a = exp(loga);

  class_call(background_functions(pba, a, y, normal_info, pvecback),
             pba->error_message,
             error_message);

  H = pvecback[pba->index_bg_H];

  dy[pba->index_bi_time] = 1./H;

  dy[pba->index_bi_tau] = 1./a/H;

  class_test(pvecback[pba->index_bg_rho_g] <= 0.,
             error_message,
             "rho_g = %e instead of strictly positive",pvecback[pba->index_bg_rho_g]);

  dy[pba->index_bi_rs] = 1./a/H/sqrt(3.*(1.+3.*pvecback[pba->index_bg_rho_b]/4./pvecback[pba->index_bg_rho_g]))*sqrt(1.-pba->K*y[pba->index_bi_rs]*y[pba->index_bi_rs]);

  rho_M = pvecback[pba->index_bg_rho_b];
  if (pba->has_cdm == _TRUE_) {
    rho_M += pvecback[pba->index_bg_rho_cdm];
  }
  if (pba->has_idm == _TRUE_){
    rho_M += pvecback[pba->index_bg_rho_idm];
  }

  dy[pba->index_bi_D] = y[pba->index_bi_D_prime]/a/H;
  {
    double G_eff_growth = 1.0;
    if (pba->hqiv.hqiv_on == _TRUE_) {
      /* G_eff/G_0 from the HQIV Friedmann equation, epoch-dependent via γ_eff(a).
       * Pre-recomb: γ_eff ≈ 0, so G_eff/G_0 ≈ 1 (standard growth).
       * Post-recomb: γ_eff = γ_theory, G_eff follows the exact formula.
       * Capped to [0.2, 5.0] for numerical stability in the ODE integrator. */
      G_eff_growth = pba->hqiv.G_eff_ratio;
      if (G_eff_growth > 5.0) G_eff_growth = 5.0;
      if (G_eff_growth < 0.2) G_eff_growth = 0.2;
    }
    dy[pba->index_bi_D_prime] = -y[pba->index_bi_D_prime] + 1.5*a*G_eff_growth*rho_M*y[pba->index_bi_D]/H;
  }

  if (pba->has_dcdm == _TRUE_) {
    dy[pba->index_bi_rho_dcdm] = -3.*y[pba->index_bi_rho_dcdm] - pba->Gamma_dcdm/H*y[pba->index_bi_rho_dcdm];
  }

  if ((pba->has_dcdm == _TRUE_) && (pba->has_dr == _TRUE_)) {
    dy[pba->index_bi_rho_dr] = -4.*y[pba->index_bi_rho_dr]+pba->Gamma_dcdm/H*y[pba->index_bi_rho_dcdm];
  }

  if (pba->has_fld == _TRUE_) {
    dy[pba->index_bi_rho_fld] = -3.*(1.+pvecback[pba->index_bg_w_fld])*y[pba->index_bi_rho_fld];
  }

  if (pba->has_scf == _TRUE_) {
    dy[pba->index_bi_phi_scf] = y[pba->index_bi_phi_prime_scf]/a/H;
    dy[pba->index_bi_phi_prime_scf] = - 2*y[pba->index_bi_phi_prime_scf] - a*dV_scf(pba,y[pba->index_bi_phi_scf])/H ;
  }

  return _SUCCESS_;

}

/**
 * At some step during the integraton of the background equations,
 * this function extracts the qantities that we want to keep memory
 * of, and stores them in a row of the background table.
 */

int background_sources(
                       double loga,
                       double * y,
                       double * dy,
                       int index_loga,
                       void * parameters_and_workspace,
                       ErrorMsg error_message
                       ) {

  struct background_parameters_and_workspace * pbpaw;
  struct background * pba;
  double a;
  double * bg_table_row;

  pbpaw = parameters_and_workspace;
  pba =  pbpaw->pba;

  bg_table_row = pba->background_table + index_loga*pba->bg_size;

  a = exp(loga);

  pba->z_table[index_loga] = MAX(0.,1./a-1.);

  pba->tau_table[index_loga] = y[pba->index_bi_tau];

  class_call(background_functions(pba, a, y, long_info, bg_table_row),
             pba->error_message,
             error_message);

  return _SUCCESS_;

}

/**
 * Evalute the typical timescale for the integration of he background
 * over loga=log(a/a_0).
 */

int background_timescale(
                         double loga,
                         void * parameters_and_workspace,
                         double * timescale,
                         ErrorMsg error_message
                         ) {

  *timescale = 1.;
  return _SUCCESS_;
}

/**
 * Function outputting the fractions Omega of the total critical density
 * today, and also the reduced fractions omega=Omega*h*h
 */

int background_output_budget(
                             struct background* pba
                             ) {

  double budget_matter, budget_radiation, budget_other,budget_neutrino;
  int index_ncdm;

  budget_matter = 0;
  budget_radiation = 0;
  budget_other = 0;
  budget_neutrino = 0;

  if (pba->background_verbose > 1) {

    printf(" ---------------------------- Budget equation ----------------------- \n");

    printf(" ---> Nonrelativistic Species \n");
    class_print_species("Bayrons",b);
    budget_matter+=pba->Omega0_b;
    if (pba->has_cdm == _TRUE_) {
      class_print_species("Cold Dark Matter",cdm);
      budget_matter+=pba->Omega0_cdm;
    }
    if (pba->has_idm == _TRUE_){
      class_print_species("Interacting DM - idr,b,g",idm);
      budget_matter+=pba->Omega0_idm;
    }
    if (pba->has_dcdm == _TRUE_) {
      class_print_species("Decaying Cold Dark Matter",dcdm);
      budget_matter+=pba->Omega0_dcdm;
    }

    if (pba->N_ncdm > 0) {
      printf(" ---> Non-Cold Dark Matter Species (incl. massive neutrinos)\n");
    }
    if (pba->N_ncdm > 0) {
      for (index_ncdm=0;index_ncdm<pba->N_ncdm;++index_ncdm) {
        printf("-> %-26s%-4d Omega = %-15g , omega = %-15g\n","Non-Cold Species Nr.",index_ncdm+1,pba->Omega0_ncdm[index_ncdm],pba->Omega0_ncdm[index_ncdm]*pba->h*pba->h);
        budget_neutrino+=pba->Omega0_ncdm[index_ncdm];
        budget_matter+=pba->Omega0_ncdm[index_ncdm];
      }
    }

    printf(" ---> Relativistic Species \n");
    class_print_species("Photons",g);
    budget_radiation+=pba->Omega0_g;
    if (pba->has_ur == _TRUE_) {
      class_print_species("Ultra-relativistic relics",ur);
      budget_radiation+=pba->Omega0_ur;
    }
    if (pba->has_dr == _TRUE_) {
      class_print_species("Dark Radiation (from decay)",dr);
      budget_radiation+=pba->Omega0_dr;
    }
    if (pba->has_idr == _TRUE_) {
      class_print_species("Interacting Dark Radiation",idr);
      budget_radiation+=pba->Omega0_idr;
    }

    if ((pba->has_lambda == _TRUE_) || (pba->has_fld == _TRUE_) || (pba->has_scf == _TRUE_) || (pba->has_curvature == _TRUE_)) {
      printf(" ---> Other Content \n");
    }
    if (pba->has_lambda == _TRUE_) {
      class_print_species("Cosmological Constant",lambda);
      budget_other+=pba->Omega0_lambda;
    }
    if (pba->has_fld == _TRUE_) {
      class_print_species("Dark Energy Fluid",fld);
      budget_other+=pba->Omega0_fld;
    }
    if (pba->has_scf == _TRUE_) {
      class_print_species("Scalar Field",scf);
      budget_other+=pba->Omega0_scf;
    }
    if (pba->has_curvature == _TRUE_) {
      class_print_species("Spatial Curvature",k);
      budget_other+=pba->Omega0_k;
    }

    printf(" ---> Total budgets \n");
    printf(" Radiation                        Omega = %-15g , omega = %-15g \n",budget_radiation,budget_radiation*pba->h*pba->h);
    printf(" Non-relativistic                 Omega = %-15g , omega = %-15g \n",budget_matter,budget_matter*pba->h*pba->h);
    if (pba->N_ncdm > 0) {
      printf(" - Non-Free-Streaming Matter      Omega = %-15g , omega = %-15g \n",pba->Omega0_nfsm,pba->Omega0_nfsm*pba->h*pba->h);
      printf(" - Non-Cold Dark Matter           Omega = %-15g , omega = %-15g \n",budget_neutrino,budget_neutrino*pba->h*pba->h);
    }
    if ((pba->has_lambda == _TRUE_) || (pba->has_fld == _TRUE_) || (pba->has_scf == _TRUE_) || (pba->has_curvature == _TRUE_)) {
      printf(" Other Content                    Omega = %-15g , omega = %-15g \n",budget_other,budget_other*pba->h*pba->h);
    }
    printf(" TOTAL                            Omega = %-15g , omega = %-15g \n",budget_radiation+budget_matter+budget_other,(budget_radiation+budget_matter+budget_other)*pba->h*pba->h);
    printf(" -------------------------------------------------------------------- \n");
  }

  return _SUCCESS_;
}

/**
 * HQIV emergent lattice solve — fully parameter-free 4D object from Planck lattice.
 *
 * Reads a pre-computed table from Python forward_4d_evolution (columns: loga, rho_r_comov, rho_b_comov, T).
 * Finds the slice where T == T0_cmb, rescales so a = 1 there ("today"), sets emergent H0, Omega_m, age, etc.
 * Lattice table densities must be in CLASS units (same normalization as background_table).
 */
int hqiv_emergent_lattice_solve(struct precision *ppr, struct background *pba) {

  FILE *f;
  char line[512];
  int i, n;
  double a_final;
  double *loga_table, *rho_r_table, *rho_b_table, *T_table;
  double gamma_theory, a_rec, da_width;
  double rho_today, H0_emergent;
  double gamma_eff, gamma_class, rho_phys_r, rho_phys_b, rho_tot, a_here, H_here;
  int idx_today;
  double frac, loga_today, scale;

  class_test(pba->hqiv_lattice_table[0] == '\0',
             pba->error_message,
             "hqiv_emergent requires hqiv_lattice_table file path");

  f = fopen(pba->hqiv_lattice_table, "r");
  class_test(f == NULL, pba->error_message, "Cannot open lattice table %s", pba->hqiv_lattice_table);

  /* Count data lines (skip lines starting with #) */
  n = 0;
  while (fgets(line, sizeof(line), f) != NULL) {
    if (line[0] == '#' || line[0] == '\n') continue;
    n++;
  }
  rewind(f);

  class_alloc(loga_table, n * sizeof(double), pba->error_message);
  class_alloc(rho_r_table, n * sizeof(double), pba->error_message);
  class_alloc(rho_b_table, n * sizeof(double), pba->error_message);
  class_alloc(T_table, n * sizeof(double), pba->error_message);

  i = 0;
  while (fgets(line, sizeof(line), f) != NULL && i < n) {
    if (line[0] == '#' || line[0] == '\n') continue;
    if (sscanf(line, "%lf %lf %lf %lf",
               &loga_table[i], &rho_r_table[i], &rho_b_table[i], &T_table[i]) == 4)
      i++;
  }
  fclose(f);
  n = i;

  class_test(n < 2, pba->error_message, "Lattice table must have at least 2 rows");

  /* Find slice where T == T0_cmb (linear interpolation) */
  idx_today = 0;
  for (i = 0; i < n - 1; i++) {
    if (T_table[i] >= pba->T_cmb && T_table[i + 1] < pba->T_cmb) {
      idx_today = i;
      break;
    }
  }
  if (T_table[n - 1] >= pba->T_cmb) idx_today = n - 1;

  frac = (T_table[idx_today + 1] != T_table[idx_today])
         ? (pba->T_cmb - T_table[idx_today]) / (T_table[idx_today + 1] - T_table[idx_today])
         : 0.0;
  if (idx_today >= n - 1) frac = 0.0;
  loga_today = loga_table[idx_today] + frac * (loga_table[idx_today + 1] - loga_table[idx_today]);

  a_final = exp(loga_today);
  scale = 1.0 / a_final;
  for (i = 0; i < n; i++)
    loga_table[i] += log(scale);

  /* Physical densities today (a = 1): rho_r_comov/a^4 = rho_r_comov, rho_b_comov/a^3 = rho_b_comov */
  rho_today = rho_r_table[idx_today] + rho_b_table[idx_today];

  /* Emergent H0 from modified Friedmann at a = 1: 3*H^2 - gamma_eff*H0*H = 3*rho, H = H0 => H0^2*(3 - gamma) = 3*rho */
  gamma_theory = pba->hqiv.gamma_hqiv;
  class_test(gamma_theory >= 3.0, pba->error_message, "hqiv_emergent requires gamma_hqiv < 3");
  H0_emergent = sqrt(3.0 * rho_today / (3.0 - gamma_theory));
  pba->H0 = H0_emergent;
  pba->h = pba->H0 * _Mpc_over_m_ / (1.e5 / _c_);  /* h = H0 in 100 km/s/Mpc units */

  /* Allocate CLASS tables */
  pba->bt_size = n;
  class_alloc(pba->tau_table, pba->bt_size * sizeof(double), pba->error_message);
  class_alloc(pba->z_table, pba->bt_size * sizeof(double), pba->error_message);
  class_alloc(pba->loga_table, pba->bt_size * sizeof(double), pba->error_message);
  class_alloc(pba->d2tau_dz2_table, pba->bt_size * sizeof(double), pba->error_message);
  class_alloc(pba->d2z_dtau2_table, pba->bt_size * sizeof(double), pba->error_message);
  class_alloc(pba->background_table, pba->bt_size * pba->bg_size * sizeof(double), pba->error_message);
  class_alloc(pba->d2background_dloga2_table, pba->bt_size * pba->bg_size * sizeof(double), pba->error_message);

  a_rec = 1.0 / 1090.0;
  da_width = 0.02 * a_rec;

  /* Fill tables: loga, z, tau (integrated), and background_table */
  pba->tau_table[0] = 0.0;
  for (i = 0; i < n; i++) {
    pba->loga_table[i] = loga_table[i];
    a_here = exp(loga_table[i]);
    pba->z_table[i] = exp(-loga_table[i]) - 1.0;

    rho_phys_r = rho_r_table[i] / (a_here * a_here * a_here * a_here);
    rho_phys_b = rho_b_table[i] / (a_here * a_here * a_here);
    rho_tot = rho_phys_r + rho_phys_b;

    gamma_eff = gamma_theory * 0.5 * (1.0 + tanh((a_here - a_rec) / da_width));
    gamma_class = gamma_eff * pba->H0;
    H_here = (gamma_class / 3.0 + sqrt((gamma_class / 3.0) * (gamma_class / 3.0) + 4.0 * rho_tot)) / 2.0;

    if (i > 0) {
      double a_prev = exp(loga_table[i - 1]);
      double H_prev = pba->background_table[(i - 1) * pba->bg_size + pba->index_bg_H];
      pba->tau_table[i] = pba->tau_table[i - 1]
                         + 0.5 * (1.0 / (a_prev * H_prev) + 1.0 / (a_here * H_here))
                           * (loga_table[i] - loga_table[i - 1]);
    }

    pba->background_table[i * pba->bg_size + pba->index_bg_a] = a_here;
    pba->background_table[i * pba->bg_size + pba->index_bg_rho_g] = rho_phys_r;
    pba->background_table[i * pba->bg_size + pba->index_bg_rho_b] = rho_phys_b;
    pba->background_table[i * pba->bg_size + pba->index_bg_H] = H_here;
    pba->background_table[i * pba->bg_size + pba->index_bg_H_prime] = 0.0;  /* spline will approximate */
    if (pba->has_cdm == _TRUE_)
      pba->background_table[i * pba->bg_size + pba->index_bg_rho_cdm] = 0.0;
    if (pba->has_lambda == _TRUE_)
      pba->background_table[i * pba->bg_size + pba->index_bg_rho_lambda] = 0.0;
    pba->background_table[i * pba->bg_size + pba->index_bg_rho_tot] = rho_tot;
    pba->background_table[i * pba->bg_size + pba->index_bg_p_tot] = rho_phys_r / 3.0;
    pba->background_table[i * pba->bg_size + pba->index_bg_p_tot_prime] = -(4.0 / 3.0) * rho_phys_r * a_here * H_here;
    pba->background_table[i * pba->bg_size + pba->index_bg_rho_crit] = H_here * H_here;
    pba->background_table[i * pba->bg_size + pba->index_bg_Omega_r] = rho_phys_r / (H_here * H_here);
    pba->background_table[i * pba->bg_size + pba->index_bg_Omega_m] = rho_phys_b / (H_here * H_here);
    pba->background_table[i * pba->bg_size + pba->index_bg_D] = a_here;
    pba->background_table[i * pba->bg_size + pba->index_bg_f] = 1.0;
    pba->background_table[i * pba->bg_size + pba->index_bg_time] = 0.0;  /* filled below */
    pba->background_table[i * pba->bg_size + pba->index_bg_rs] = 0.0;
  }

  /* Proper time integration */
  for (i = 1; i < n; i++) {
    double H_prev = pba->background_table[(i - 1) * pba->bg_size + pba->index_bg_H];
    double H_curr = pba->background_table[i * pba->bg_size + pba->index_bg_H];
    pba->background_table[i * pba->bg_size + pba->index_bg_time] =
      pba->background_table[(i - 1) * pba->bg_size + pba->index_bg_time]
      + 0.5 * (1.0 / H_prev + 1.0 / H_curr) * (pba->loga_table[i] - pba->loga_table[i - 1]);
  }

  pba->age = pba->background_table[(n - 1) * pba->bg_size + pba->index_bg_time] / _Gyr_over_Mpc_;
  pba->conformal_age = pba->tau_table[n - 1];

  /* Conformal/angular/luminosity distances at each row */
  for (i = 0; i < n; i++) {
    double conformal_distance = pba->conformal_age - pba->tau_table[i];
    double comoving_radius = conformal_distance;  /* K = 0 */
    pba->background_table[i * pba->bg_size + pba->index_bg_conf_distance] = conformal_distance;
    pba->background_table[i * pba->bg_size + pba->index_bg_ang_distance] = comoving_radius / (1.0 + pba->z_table[i]);
    pba->background_table[i * pba->bg_size + pba->index_bg_lum_distance] = comoving_radius * (1.0 + pba->z_table[i]);
  }

  class_call(array_spline_table_lines(pba->z_table, pba->bt_size, pba->tau_table, 1,
                                      pba->d2tau_dz2_table, _SPLINE_EST_DERIV_, pba->error_message),
             pba->error_message, pba->error_message);
  class_call(array_spline_table_lines(pba->tau_table, pba->bt_size, pba->z_table, 1,
                                      pba->d2z_dtau2_table, _SPLINE_EST_DERIV_, pba->error_message),
             pba->error_message, pba->error_message);
  class_call(array_spline_table_lines(pba->loga_table, pba->bt_size, pba->background_table, pba->bg_size,
                                      pba->d2background_dloga2_table, _SPLINE_EST_DERIV_, pba->error_message),
             pba->error_message, pba->error_message);

  pba->Omega0_m = pba->background_table[(n - 1) * pba->bg_size + pba->index_bg_Omega_m];
  pba->Omega0_r = pba->background_table[(n - 1) * pba->bg_size + pba->index_bg_Omega_r];
  pba->Omega0_b = pba->Omega0_m;
  pba->Omega0_g = pba->Omega0_r;
  pba->Omega0_ur = 0.0;
  pba->Omega0_de = gamma_theory / 3.0;  /* horizon term replaces dark energy */
  pba->Omega0_k = 0.0;
  pba->K = 0.0;
  pba->sgnK = 0;

  free(loga_table);
  free(rho_r_table);
  free(rho_b_table);
  free(T_table);

  if (pba->background_verbose > 0) {
    printf(" -> HQIV emergent lattice mode activated\n");
    printf(" -> T0 = %.4f K defines today (a=1)\n", pba->T_cmb);
    printf(" -> Emergent H0 = %.4f km/s/Mpc\n", pba->H0 * 299792.458);
    printf(" -> Emergent Omega_m (baryon-only) = %.6f\n", pba->Omega0_m);
    printf(" -> Emergent global age = %.2f Gyr\n", pba->age);
  }

  return _SUCCESS_;
}

/**
 * Scalar field potential and its derivatives with respect to the field _scf
 */

double V_e_scf(struct background *pba,
               double phi
               ) {
  double scf_lambda = pba->scf_parameters[0];

  return  exp(-scf_lambda*phi);
}

double dV_e_scf(struct background *pba,
                double phi
                ) {
  double scf_lambda = pba->scf_parameters[0];

  return -scf_lambda*V_e_scf(pba,phi);
}

double ddV_e_scf(struct background *pba,
                 double phi
                 ) {
  double scf_lambda = pba->scf_parameters[0];

  return pow(-scf_lambda,2)*V_e_scf(pba,phi);
}

double V_p_scf(
               struct background *pba,
               double phi) {
  double scf_alpha  = pba->scf_parameters[1];
  double scf_A      = pba->scf_parameters[2];
  double scf_B      = pba->scf_parameters[3];

  return  pow(phi - scf_B,  scf_alpha) +  scf_A;
}

double dV_p_scf(
                struct background *pba,
                double phi) {

  double scf_alpha  = pba->scf_parameters[1];
  double scf_B      = pba->scf_parameters[3];

  return   scf_alpha*pow(phi -  scf_B,  scf_alpha - 1);
}

double ddV_p_scf(
                 struct background *pba,
                 double phi) {
  double scf_alpha  = pba->scf_parameters[1];
  double scf_B      = pba->scf_parameters[3];

  return  scf_alpha*(scf_alpha - 1.)*pow(phi -  scf_B,  scf_alpha - 2);
}

double V_scf(
             struct background *pba,
             double phi) {
  return  V_e_scf(pba,phi)*V_p_scf(pba,phi);
}

double dV_scf(
              struct background *pba,
              double phi) {
  return dV_e_scf(pba,phi)*V_p_scf(pba,phi) + V_e_scf(pba,phi)*dV_p_scf(pba,phi);
}

double ddV_scf(
               struct background *pba,
               double phi) {
  return ddV_e_scf(pba,phi)*V_p_scf(pba,phi) + 2*dV_e_scf(pba,phi)*dV_p_scf(pba,phi) + V_e_scf(pba,phi)*ddV_p_scf(pba,phi);
}