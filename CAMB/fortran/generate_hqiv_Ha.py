#!/usr/bin/env python3
"""
Generate correct HQIV H(a)/H0 table for CAMB.

Background: H(a)/H0 = sqrt(Omega_r/a^4 + Omega_b/a^3 + Omega_eff)
where Omega_eff enforces flatness (acts as effective horizon energy density).
Varying G(a) and inertia factor act only at the perturbation level.

Parameters match params_hqiv_covariant_test.ini:
  H0 = 70 km/s/Mpc, ombh2 = 0.048, omch2 = 0, omnuh2 = 0, omk = 0
  T_CMB = 2.7255 K, N_eff = 3.044
"""
import numpy as np

H0 = 70.0
h = H0 / 100.0
h2 = h**2

ombh2 = 0.048
Omega_b = ombh2 / h2

T_CMB = 2.7255
N_eff = 3.044
Omega_gamma_h2 = 2.4728e-5 * (T_CMB / 2.7255)**4
Omega_gamma = Omega_gamma_h2 / h2
Omega_nu = N_eff * (7.0/8.0) * (4.0/11.0)**(4.0/3.0) * Omega_gamma
Omega_r = Omega_gamma + Omega_nu

Omega_eff = 1.0 - Omega_r - Omega_b

print(f"Omega_b   = {Omega_b:.8f}")
print(f"Omega_r   = {Omega_r:.6e}")
print(f"Omega_eff = {Omega_eff:.8f}")

n_pts = 5001
a_min = 1e-7
a_max = 1.0

ln_a = np.linspace(np.log(a_min), np.log(a_max), n_pts)
a_arr = np.exp(ln_a)
a_arr[-1] = 1.0

H_over_H0 = np.sqrt(Omega_r / a_arr**4 + Omega_b / a_arr**3 + Omega_eff)

print(f"\nTable: {n_pts} points, a=[{a_arr[0]:.2e}, {a_arr[-1]:.2e}]")
print(f"H/H0 range: [{H_over_H0[-1]:.6f}, {H_over_H0[0]:.6e}]")
print(f"H/H0(a=1) = {H_over_H0[-1]:.10f}")

age_integrand = 1.0 / (a_arr * H_over_H0 * H0 * 1e3 / 3.0856775814913673e22)
age_s = np.trapz(age_integrand, a_arr)
age_gyr = age_s / (365.25e9 * 86400)
print(f"Approx universe age = {age_gyr:.2f} Gyr")

out = np.column_stack([a_arr, H_over_H0])
np.savetxt("hqiv_Ha.txt", out, header="a  H_over_H0", comments="", fmt="%.10e")
print("Written hqiv_Ha.txt")
