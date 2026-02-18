"""
HQIV Cosmology — Background integration.

A_eff = A_std + horizon term. With varying G(t) = G0 (Θ0/Θ(t))^0.6.
Radiation ∝ a^{-4}, matter ∝ a^{-3}. Integrates in ln(a).
Dimensionally: ä/a has units 1/s^2; we use horizon term β H^2 (so A_eff in 1/s^2).
Produces: H(a), age, proper time at z=14, d_A(z=14), and H(a) table for CLASS.
"""

import numpy as np

try:
    from scipy.integrate import solve_ivp
    _has_scipy = True
except ImportError:
    _has_scipy = False

# Physical constants (SI)
c = 2.99792458e8   # m/s
G0 = 6.67430e-11   # m^3 kg^{-1} s^{-2}
H0_SI = 2.2e-18    # ~70 km/s/Mpc in 1/s
Mpc = 3.0856775814913673e22  # m

# HQIV parameters (from paper)
beta = 1.02
Omega_m_vis = 0.048
Omega_r0 = 9.0e-5    # radiation today
h = 0.70
H0 = h * 100 * 1e3 / Mpc  # 1/s
H0_km = h * 100       # km/s/Mpc

# Today's horizon Θ0 = 2c/H0
Theta0 = 2 * c / H0


def G_ratio(a, H_a):
    """G(t)/G0 = (Θ0/Θ(t))^0.6, Θ = 2c/H."""
    Theta = 2 * c / (H_a + 1e-100)
    return (Theta0 / Theta) ** 0.6


def hqiv_rhs(ln_a, y):
    """
    State y = [H]. We evolve H.
    A_eff = A_std + β H^2 (dimensionally 1/s^2; horizon acts like effective Λ).
    dH/d(ln a) = (A_eff - H^2)/H = A_eff/H - H.
    """
    a = np.exp(ln_a)
    H = max(y[0], 1e-4 * H0)  # avoid RHS blow-up when H -> 0
    if H <= 0:
        return [0.0]

    rho_m = 3 * H0**2 * Omega_m_vis / (8 * np.pi * G0) * (1 / a**3)
    rho_r = 3 * H0**2 * Omega_r0 / (8 * np.pi * G0) * (1 / a**4)
    G_rat = G_ratio(a, H)
    A_std = -(4 * np.pi * G_rat * G0 / 3) * (rho_m + 2 * rho_r)
    # Horizon term: β H^2 (1/s^2); a_min = β c H → effective A term β H^2
    A_eff = A_std + beta * H**2
    dH_dln_a = (A_eff - H**2) / H
    return [dH_dln_a]


def hqiv_rhs_log(ln_a, y):
    """State y = [ln(H)]. du/d(ln a) = (A_eff - H^2)/H^2 = A_eff/H^2 - 1. Better scaled."""
    ln_H = y[0]
    H = np.exp(ln_H)
    H = max(H, 1e-4 * H0)
    a = np.exp(ln_a)
    rho_m = 3 * H0**2 * Omega_m_vis / (8 * np.pi * G0) * (1 / a**3)
    rho_r = 3 * H0**2 * Omega_r0 / (8 * np.pi * G0) * (1 / a**4)
    G_rat = G_ratio(a, H)
    A_std = -(4 * np.pi * G_rat * G0 / 3) * (rho_m + 2 * rho_r)
    A_eff = A_std + beta * H**2
    du_dln_a = (A_eff - H**2) / (H**2)
    return [du_dln_a]


def _rk4_step(ln_a, y, dln_a, rhs):
    """Single RK4 step: y_new = y + dln_a * (k1+2*k2+2*k3+k4)/6."""
    k1 = np.array(rhs(ln_a, y))
    k2 = np.array(rhs(ln_a + 0.5 * dln_a, y + 0.5 * dln_a * k1))
    k3 = np.array(rhs(ln_a + 0.5 * dln_a, y + 0.5 * dln_a * k2))
    k4 = np.array(rhs(ln_a + dln_a, y + dln_a * k3))
    return y + (dln_a / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)


def _run_background_numpy(ln_a_vals, H_init, rhs):
    """Integrate with RK4 over ln_a_vals (numpy-only)."""
    y = np.array([H_init], dtype=float)
    H_arr = np.zeros(len(ln_a_vals))
    H_arr[0] = H_init
    for i in range(len(ln_a_vals) - 1):
        dln_a = ln_a_vals[i + 1] - ln_a_vals[i]
        y = _rk4_step(ln_a_vals[i], y, dln_a, rhs)
        # Keep H positive
        y[0] = max(y[0], 1e-30)
        H_arr[i + 1] = y[0]
    return H_arr


def run_background(a_min=1e-6, a_max=1.0, n_pts=4000):
    """
    Integrate from a_min to a_max. Initial H set from matter+radiation only.
    Uses scipy.solve_ivp (LSODA + ln(H) state) if available, else numpy RK4.
    a_min=1e-6 avoids stiff early phase; early age from 1e-12 to 1e-6 is negligible.
    """
    ln_a_min = np.log(a_min)
    ln_a_max = np.log(a_max)
    H2_init = (8 * np.pi * G0 / 3) * (
        3 * H0**2 * Omega_r0 / (8 * np.pi * G0) / a_min**4
        + 3 * H0**2 * Omega_m_vis / (8 * np.pi * G0) / a_min**3
    )
    H_init = np.sqrt(max(H2_init, 1e-100))
    ln_a_vals = np.linspace(ln_a_min, ln_a_max, n_pts)

    if _has_scipy:
        # Stiff ODE: evolve ln(H) for better scaling; Radau
        y0_log = [np.log(H_init)]
        sol = solve_ivp(
            hqiv_rhs_log,
            (ln_a_min, ln_a_max),
            y0_log,
            method="LSODA",  # auto stiff/non-stiff
            rtol=1e-5,
            atol=1e-9,
        )
        if not sol.success:
            raise RuntimeError("HQIV background integration failed: " + sol.message)
        H_sol = np.exp(sol.y[0])
        # Interpolate onto uniform ln(a) grid; floor H so age/d_A integrals stay finite
        H_arr = np.interp(ln_a_vals, sol.t, H_sol)
        H_arr = np.maximum(H_arr, 1e-4 * H0)
    else:
        H_arr = _run_background_numpy(ln_a_vals, H_init, lambda ln_a, y: hqiv_rhs(ln_a, list(y)))

    a_arr = np.exp(ln_a_vals)
    return a_arr, H_arr


def age_gyr(a_arr, H_arr):
    """Proper time from a_min to a: t = int da/(a H). dt = da/(a H), so t = int 1/(a H) da."""
    a_arr = np.asarray(a_arr)
    H_arr = np.asarray(H_arr)
    # t(a) = int_{a_min}^{a} da'/(a' H(a'))
    da = np.diff(a_arr)
    mid_a = (a_arr[:-1] + a_arr[1:]) / 2
    mid_H = (H_arr[:-1] + H_arr[1:]) / 2
    integrand = 1.0 / (mid_a * mid_H)
    t_s = np.zeros_like(a_arr)
    t_s[1:] = np.cumsum(integrand * da)
    return t_s / (3.1536e16)  # in Gyr


def proper_time_at_z(a_arr, H_arr, z_target):
    """Proper time (Gyr) from big bang to scale factor a = 1/(1+z)."""
    a_target = 1.0 / (1 + z_target)
    t_gyr = age_gyr(a_arr, H_arr)
    # Interpolate t at a_target
    idx = np.searchsorted(a_arr, a_target)
    if idx >= len(a_arr):
        return float(t_gyr[-1])
    if idx == 0:
        return 0.0
    t_val = np.interp(a_target, a_arr, t_gyr)
    return float(t_val)


def angular_diameter_distance(a_arr, H_arr, z_target):
    """d_A(z) = (1/(1+z)) * c * int_{a(z)}^{1} da / (a^2 H(a)) in Mpc."""
    a_target = 1.0 / (1 + z_target)
    mask = (a_arr >= a_target) & (a_arr <= 1.0)
    if not np.any(mask):
        return 0.0
    a_slice = np.concatenate([[a_target], a_arr[mask]])
    a_slice = np.unique(a_slice)
    H_slice = np.interp(a_slice, a_arr, H_arr)
    da = np.diff(a_slice)
    mid_a = (a_slice[:-1] + a_slice[1:]) / 2
    mid_H = np.interp(mid_a, a_arr, H_arr)
    integrand = 1.0 / (mid_a**2 * mid_H)
    comoving_m = c * np.sum(integrand * da)
    d_A_Mpc = (comoving_m / (1 + z_target)) / Mpc
    return float(d_A_Mpc)


def main():
    a_arr, H_arr = run_background()
    t_gyr = age_gyr(a_arr, H_arr)
    age_today_gyr = float(np.interp(1.0, a_arr, t_gyr))
    t_z14_gyr = proper_time_at_z(a_arr, H_arr, 14.0)
    dA_z14 = angular_diameter_distance(a_arr, H_arr, 14.0)

    print("HQIV background (sandbox)")
    print("  Universe age today: {:.2f} Gyr".format(age_today_gyr))
    print("  Proper time at z=14: {:.0f} Myr".format(t_z14_gyr * 1000))
    print("  d_A(z=14): {:.2f} Gpc".format(dA_z14 / 1000))
    print("  H0: {:.2f} km/s/Mpc".format(H0_km))

    # Save H(a) table for CLASS: columns a, H/H0
    H_over_H0 = H_arr / H0
    out = np.column_stack([a_arr, H_over_H0])
    np.savetxt(
        "hqiv_Ha.txt",
        out,
        header="a  H_over_H0",
        comments="",
        fmt="%.6e",
    )
    print("  Written hqiv_Ha.txt (a, H/H0) for CLASS.")

    return a_arr, H_arr, age_today_gyr, t_z14_gyr, dA_z14


if __name__ == "__main__":
    main()
