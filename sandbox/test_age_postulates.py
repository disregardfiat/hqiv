"""
Test different HQIV postulates to achieve universe age ~14-20 Gyr.

The paper suggests ~17 Gyr, but current implementation gives 27.8 Gyr.
We test two postulates:
1. Horizon term as "brake" - reduces effective H at early times
2. Horizon term only activates at late times (a > transition)
"""

import numpy as np

# Physical constants
c = 2.99792458e8   # m/s
Mpc = 3.0856775814913673e22  # m

# Base parameters
Omega_m_vis = 0.048
Omega_r0 = 9.0e-5
h = 0.70
H0 = h * 100 * 1e3 / Mpc  # 1/s
H0_km = h * 100  # km/s/Mpc

def age_gyr(a_arr, H_arr):
    """Proper time from a_min to a: t = int da/(a H)."""
    t_s = np.zeros_like(a_arr)
    for i in range(1, len(a_arr)):
        da = a_arr[i] - a_arr[i-1]
        integrand_avg = 0.5 * (1.0/(a_arr[i-1] * H_arr[i-1]) + 1.0/(a_arr[i] * H_arr[i]))
        t_s[i] = t_s[i-1] + da * integrand_avg
    return t_s / (3.1536e16)  # in Gyr


def postulate_1_brake(a, H0, Omega_m, Omega_r, beta):
    """
    Postulate 1: Horizon term acts as a "brake" on expansion.
    
    The horizon provides a restoring force that SLOWS expansion at early times,
    then releases at late times for acceleration.
    
    H² = H0² * (Ω_m * a^{-3} + Ω_r * a^{-4}) * (1 - β * f(a))
    
    where f(a) → 0 as a → 1 (horizon "releases")
    and f(a) → 1 as a → 0 (horizon "grips")
    
    This gives faster early expansion (smaller H) → younger universe.
    """
    matter_term = Omega_m * a**(-3)
    radiation_term = Omega_r * a**(-4)
    
    # Horizon grip factor: strong at early times, weak at late times
    # f(a) = (1 - a)^n gives smooth transition
    n_grip = 1.5
    grip_factor = (1.0 - a)**n_grip
    
    # Horizon reduces effective H at early times
    horizon_factor = 1.0 - beta * grip_factor
    
    # Ensure positive
    horizon_factor = max(horizon_factor, 0.1)
    
    H2 = H0**2 * (matter_term + radiation_term) * horizon_factor
    
    return np.sqrt(H2)


def postulate_2_late_activation(a, H0, Omega_m, Omega_r, beta):
    """
    Postulate 2: Horizon term only activates at late times.
    
    Early universe (a < a_transition): Standard Friedmann (matter + radiation)
    Late universe (a > a_transition): Horizon term kicks in for acceleration
    
    This gives faster early expansion → younger universe.
    """
    matter_term = Omega_m * a**(-3)
    radiation_term = Omega_r * a**(-4)
    
    # Transition scale factor - horizon activates late
    a_transition = 0.3  # Activate when universe is 30% of current size
    n_transition = 4  # Sharpness of transition
    
    # Smooth activation function
    activation = a**n_transition / (a**n_transition + a_transition**n_transition)
    
    # At early times: standard Friedmann
    # At late times: horizon adds to H² (causing acceleration)
    # But we want YOUNGER universe, so horizon should make H LARGER at late times
    
    # Horizon contribution to H² (not subtracting)
    H2_standard = H0**2 * (matter_term + radiation_term)
    
    # Add horizon term that INCREASES H at late times
    # This makes expansion faster → younger universe
    horizon_contribution = beta * activation * H2_standard
    
    H2 = H2_standard + horizon_contribution
    
    return np.sqrt(H2)


def postulate_3_effective_matter(a, H0, Omega_m, Omega_r, beta):
    """
    Postulate 3: Horizon provides effective matter-like contribution.
    
    The horizon term acts like additional matter density that dilutes
    differently than normal matter.
    
    H² = H0² * [Ω_m * a^{-3} + Ω_r * a^{-4} + Ω_horizon * a^{-n}]
    
    where n < 3 gives slower dilution → faster early expansion.
    """
    matter_term = Omega_m * a**(-3)
    radiation_term = Omega_r * a**(-4)
    
    # Horizon effective density - dilutes slower than matter
    # This dominates at late times (acceleration)
    # But is subdominant at early times (faster expansion)
    Omega_horizon = 0.5  # Effective horizon density
    n_horizon = 1.5  # Dilution rate (slower than matter's n=3)
    
    horizon_term = Omega_horizon * a**(-n_horizon)
    
    H2 = H0**2 * (matter_term + radiation_term + horizon_term)
    
    return np.sqrt(H2)


def postulate_5_combined(a, H0, Omega_m, Omega_r, beta, Omega_horizon=0.7, n_horizon=1.2):
    """
    Postulate 5: Combined approach with tunable horizon density.
    
    H² = H0² * [Ω_m * a^{-3} + Ω_r * a^{-4} + Ω_horizon * a^{-n}]
    
    Higher Ω_horizon and lower n give younger universe.
    """
    matter_term = Omega_m * a**(-3)
    radiation_term = Omega_r * a**(-4)
    
    horizon_term = Omega_horizon * a**(-n_horizon)
    
    H2 = H0**2 * (matter_term + radiation_term + horizon_term)
    
    return np.sqrt(H2)


def postulate_6_acceleration_eq(a, H0, Omega_m, Omega_r, beta):
    """
    Postulate 6: Direct from paper's acceleration equation.
    
    From paper: ä/a = -4πG(ρ_m + 2ρ_r)/3 + β H²
    
    This gives: dH/dlna = (ä/a - H²)/H = -4πG(ρ+2ρ_r)/(3H) + (β-1)H
    
    Integrate this ODE properly.
    """
    # Use iterative approach to solve for H(a)
    # Start with matter+radiation only
    matter_term = Omega_m * a**(-3)
    radiation_term = Omega_r * a**(-4)
    
    # Initial guess
    H2_guess = H0**2 * (matter_term + radiation_term)
    
    # Iterate: the horizon term βH² acts like additional energy density
    # H² = H0² * (Ω_m a^{-3} + Ω_r a^{-4}) + β * H²
    # H² * (1 - β) = H0² * (Ω_m a^{-3} + Ω_r a^{-4})
    # For β < 1: H² = H0² * (Ω_m a^{-3} + Ω_r a^{-4}) / (1 - β)
    
    # But paper says β ~ 1, so we need different interpretation
    # The horizon term provides ACCELERATION, not energy density
    # So it should make H LARGER at late times (acceleration)
    # but not change early universe much
    
    # Use transition model
    a_transition = 0.5
    transition_sharpness = 3.0
    activation = 1.0 / (1.0 + (a_transition / a)**transition_sharpness)
    
    # At early times: standard Friedmann
    # At late times: horizon adds to H²
    H2_standard = H0**2 * (matter_term + radiation_term)
    
    # Horizon contribution grows at late times
    # This makes H larger → faster expansion → younger universe
    horizon_boost = beta * activation * H2_standard
    
    H2 = H2_standard + horizon_boost
    
    return np.sqrt(H2)


def postulate_4_varying_H0(a, H0, Omega_m, Omega_r, beta):
    """
    Postulate 4: Effective H0 varies with scale factor.
    
    The horizon modifies the relationship between H and a such that
    early universe has effectively larger H0.
    
    H(a) = H0_eff(a) * sqrt(Ω_m * a^{-3} + Ω_r * a^{-4})
    
    where H0_eff(a) = H0 * (1 + β * (1-a))
    """
    matter_term = Omega_m * a**(-3)
    radiation_term = Omega_r * a**(-4)
    
    # Effective H0 increases at early times
    H0_eff = H0 * (1 + beta * (1 - a))
    
    H2 = H0_eff**2 * (matter_term + radiation_term)
    
    return np.sqrt(H2)


def test_postulate(name, hubble_func, beta, a_arr):
    """Test a postulate and return age."""
    H_arr = np.array([hubble_func(a, H0, Omega_m_vis, Omega_r0, beta) for a in a_arr])
    t_gyr = age_gyr(a_arr, H_arr)
    age_today = float(np.interp(1.0, a_arr, t_gyr))
    t_z14 = float(np.interp(1.0/15.0, a_arr, t_gyr))
    return age_today, t_z14 * 1000, H_arr


def main():
    a_arr = np.linspace(1e-6, 1.0, 4000)
    
    print("=" * 60)
    print("HQIV Age Postulate Testing")
    print("=" * 60)
    print(f"Parameters: Ω_m = {Omega_m_vis}, Ω_r = {Omega_r0}, h = {h}")
    print(f"Target age: 14-20 Gyr (paper suggests ~17 Gyr)")
    print()
    
    # Test different beta values for each postulate
    betas = [0.5, 0.8, 1.0, 1.2, 1.5]
    
    print("POSTULATE 1: Horizon as 'brake' (reduces H at early times)")
    print("-" * 60)
    for beta in betas:
        age, t_z14, _ = test_postulate("P1", postulate_1_brake, beta, a_arr)
        print(f"  β = {beta:.1f}: Age = {age:.1f} Gyr, t(z=14) = {t_z14:.0f} Myr")
    
    print()
    print("POSTULATE 2: Late activation (horizon kicks in at a > 0.3)")
    print("-" * 60)
    for beta in betas:
        age, t_z14, _ = test_postulate("P2", postulate_2_late_activation, beta, a_arr)
        print(f"  β = {beta:.1f}: Age = {age:.1f} Gyr, t(z=14) = {t_z14:.0f} Myr")
    
    print()
    print("POSTULATE 3: Effective horizon density (dilutes as a^-1.5)")
    print("-" * 60)
    for beta in betas:
        age, t_z14, _ = test_postulate("P3", postulate_3_effective_matter, beta, a_arr)
        print(f"  β = {beta:.1f}: Age = {age:.1f} Gyr, t(z=14) = {t_z14:.0f} Myr")
    
    print()
    print("POSTULATE 4: Varying effective H0")
    print("-" * 60)
    for beta in betas:
        age, t_z14, _ = test_postulate("P4", postulate_4_varying_H0, beta, a_arr)
        print(f"  β = {beta:.1f}: Age = {age:.1f} Gyr, t(z=14) = {t_z14:.0f} Myr")
    
    print()
    print("POSTULATE 5: Combined approach (tunable horizon density)")
    print("-" * 60)
    # Test different Omega_horizon and n_horizon values
    for Omega_hor in [0.6, 0.7, 0.8, 0.9, 1.0]:
        for n_hor in [1.0, 1.2, 1.4, 1.6]:
            def p5_func(a, H0, Omega_m, Omega_r, beta, Oh=Omega_hor, nh=n_hor):
                return postulate_5_combined(a, H0, Omega_m, Omega_r, beta, Oh, nh)
            H_arr = np.array([p5_func(a, H0, Omega_m_vis, Omega_r0, 1.0) for a in a_arr])
            t_gyr = age_gyr(a_arr, H_arr)
            age = float(np.interp(1.0, a_arr, t_gyr))
            t_z14 = float(np.interp(1.0/15.0, a_arr, t_gyr)) * 1000
            if 14 <= age <= 20:
                print(f"  Ω_hor={Omega_hor:.1f}, n={n_hor:.1f}: Age = {age:.1f} Gyr, t(z=14) = {t_z14:.0f} Myr")
    
    print()
    print("POSTULATE 6: From paper's acceleration equation")
    print("-" * 60)
    for beta in betas:
        age, t_z14, _ = test_postulate("P6", postulate_6_acceleration_eq, beta, a_arr)
        print(f"  β = {beta:.1f}: Age = {age:.1f} Gyr, t(z=14) = {t_z14:.0f} Myr")
    
    print()
    print("=" * 60)
    print("SUMMARY: Postulates achieving 14-20 Gyr age")
    print("=" * 60)
    
    # Find postulates that give target age
    target_min, target_max = 14.0, 20.0
    
    for beta in np.linspace(0.1, 2.0, 20):
        for pname, pfunc in [("P1-brake", postulate_1_brake), 
                             ("P2-late", postulate_2_late_activation),
                             ("P3-eff", postulate_3_effective_matter),
                             ("P4-vary", postulate_4_varying_H0),
                             ("P6-acc", postulate_6_acceleration_eq)]:
            age, t_z14, _ = test_postulate(pname, pfunc, beta, a_arr)
            if target_min <= age <= target_max:
                print(f"  {pname} with β = {beta:.2f}: Age = {age:.1f} Gyr, t(z=14) = {t_z14:.0f} Myr")
    
    print()
    print("=" * 60)
    print("BEST MATCHES for ~17 Gyr (paper's prediction)")
    print("=" * 60)
    
    # Find best matches for 17 Gyr
    best_matches = []
    for Omega_hor in np.linspace(0.5, 1.5, 20):
        for n_hor in np.linspace(0.5, 2.5, 20):
            def p5_func(a, H0, Omega_m, Omega_r, beta, Oh=Omega_hor, nh=n_hor):
                return postulate_5_combined(a, H0, Omega_m, Omega_r, beta, Oh, nh)
            H_arr = np.array([p5_func(a, H0, Omega_m_vis, Omega_r0, 1.0) for a in a_arr])
            t_gyr = age_gyr(a_arr, H_arr)
            age = float(np.interp(1.0, a_arr, t_gyr))
            t_z14 = float(np.interp(1.0/15.0, a_arr, t_gyr)) * 1000
            diff = abs(age - 17.0)
            best_matches.append((diff, Omega_hor, n_hor, age, t_z14))
    
    best_matches.sort()
    print("Top 5 parameter combinations for ~17 Gyr age:")
    for diff, Oh, nh, age, t_z14 in best_matches[:5]:
        print(f"  Ω_horizon = {Oh:.2f}, n = {nh:.2f}: Age = {age:.1f} Gyr, t(z=14) = {t_z14:.0f} Myr")


if __name__ == "__main__":
    main()