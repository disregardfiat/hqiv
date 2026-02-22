"""
HQIV Cosmology — Scalar Perturbation Module.

This module implements the linear scalar perturbation equations for both
standard ΛCDM and the Horizon-Quantized Informational Vacuum (HQIV) framework.

Paper Reference: Section 7 "Linear Perturbation Equations"
- Eq. 9: Continuity equation (δ' + θ = -3Φ')
- Eq. 10: Modified Euler equation with inertia reduction
- Eq. 11: Poisson equation with horizon correction

All equations are in Newtonian gauge and Fourier space.

Author: HQIV Team
"""

import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import interp1d

# Import background module
from .background import c, Mpc_m, G0_SI

# =============================================================================
# Scalar Perturbation Equations
# =============================================================================

def inertia_reduction_factor(alpha, phi, form='thermo'):
    """
    Compute the inertia reduction factor f(α, φ) with thermodynamic floor.
    
    From the action (β-free formulation):
        f(α, φ) = max( α / (α + c φ / 6), f_min )
    
    where:
        - a_min = χ × c × φ (thermodynamic floor)
        - χ ≈ 0.172 from full light-cone average (Brodie's overlap integral)
        - f_min ≈ 0.01 from informational cutoff
    
    This is equivalent to:
        f = α / (α + a_min / 6)
    
    Parameters
    ----------
    alpha : float or array
        Local acceleration magnitude |a_local|
    phi : float or array
        Local inverse horizon scale φ ≈ c H(a)
    form : str
        'thermo' for thermodynamic floor (default)
        'sqrt' for HQIV square-root form
        'brodie' for Brodie's linear form without floor
        
    Returns
    -------
    f : float or array
        Inertia reduction factor (f_min ≤ f ≤ 1)
    
    Notes
    -----
    When α >> a_min (high acceleration), f → 1 (standard inertia).
    When α ~ a_min (low acceleration), f < 1 (reduced inertia).
    The floor f_min ensures the model doesn't break in deep MOND regime.
    """
    if form == 'thermo':
        # Thermodynamic floor from Brodie's thermodynamics
        # χ from full light-cone average
        chi = 0.172  # Thermodynamic efficiency factor
        
        # Minimum acceleration from horizon information
        a_min = chi * c * phi
        
        # Brodie's linear form with thermodynamic floor
        f = alpha / (alpha + a_min / 6.0)
        
        # Floor from informational cutoff
        f_min = 0.01
        return np.maximum(f, f_min)
    
    elif form == 'sqrt':
        # HQIV square-root form (paper after eq. 8)
        ratio = c * phi / alpha
        # Clamp to avoid negative values
        ratio = np.minimum(ratio, 0.99)
        return np.sqrt(1.0 - ratio)
    
    elif form == 'brodie':
        # Brodie's form without floor
        return alpha / (alpha + c * phi / 6.0)
    
    else:
        raise ValueError(f"Unknown form: {form}")


def horizon_correction_term(k, a, background, phi=None):
    """
    Compute the horizon correction term in the Poisson equation.
    
    Paper eq. 11:
        ∇²Φ = 4πG(a) ρ_m δ + horizon_correction
    
    The horizon correction arises from the modified Einstein equation:
        G_μν + (2c²/Θ) g_μν = (8πG/c⁴) T_μν
    
    Parameters
    ----------
    k : float
        Comoving wavenumber [1/Mpc]
    a : float
        Scale factor
    background : CosmologicalBackground
        Background cosmology object
    phi : float, optional
        Local inverse horizon scale. If None, computed from background.
        
    Returns
    -------
    correction : float
        Horizon correction term [dimensionless, in Fourier space]
    """
    if not background.hqiv_on:
        return 0.0
    
    if phi is None:
        phi = background.phi_horizon(a)
    
    # The horizon term modifies the Poisson equation on large scales
    # (small k). A simple model:
    # correction ∝ (k_horizon / k)² where k_horizon ~ H(a)/c
    
    k_horizon = background.H(a) / c * Mpc_m  # [1/Mpc]
    
    # Smooth transition: no correction on small scales (large k)
    # Full correction on large scales (small k)
    scale_factor = 1.0 / (1.0 + (k / k_horizon)**2)
    
    # Magnitude from horizon term
    beta = background.beta_a(a)
    H = background.get_H(a)
    
    # The horizon term contributes an effective density
    correction = beta * (H / background.H0)**2 * scale_factor
    
    return correction


class ScalarPerturbations:
    """
    Solver for linear scalar perturbations in Fourier space.
    
    Implements the continuity, Euler, and Poisson equations in
    Newtonian gauge for both ΛCDM and HQIV frameworks.
    
    Parameters
    ----------
    background : CosmologicalBackground
        Background cosmology object
    k : float
        Comoving wavenumber [1/Mpc]
    hqiv_on : bool
        If True, include HQIV modifications
    inertia_form : str
        Form of inertia reduction: 'sqrt' or 'brodie'
    
    State Vector
    ------------
    y = [δ, θ, Φ]
    where:
        δ = matter density contrast
        θ = velocity divergence (∇·v in Fourier space)
        Φ = Newtonian potential
    
    Notes
    -----
    The equations (paper eqs. 9-11) in conformal time η:
    
    Continuity: δ' + θ = -3 Φ'
    
    Euler (modified): θ' + ℋ θ = -k² Φ / f(α, φ)
    
    Poisson: k² Φ = -(3/2) ℋ² Ω_m(a) δ a² + horizon_term × a²
    
    where ℋ = aH is the conformal Hubble parameter.
    """
    
    def __init__(self, background, k, hqiv_on=True, inertia_form='sqrt'):
        self.bg = background
        self.k = k  # [1/Mpc]
        self.hqiv_on = hqiv_on and background.hqiv_on
        self.inertia_form = inertia_form
        
        # Storage for solution
        self.a_arr = None
        self.delta_arr = None
        self.theta_arr = None
        self.Phi_arr = None
        
    def Omega_m_a(self, a):
        """
        Matter density parameter at scale factor a.
        
        Ω_m(a) = Ω_m a⁻³ / (H(a)/H0)²
        """
        E = self.bg.H_over_H0(a)
        return self.bg.Om_m * a**(-3) / E**2
    
    def conformal_H(self, a):
        """
        Conformal Hubble parameter ℋ = aH in units of 1/Mpc.
        
        Note: H is in 1/s, so ℋ = aH/c gives units of 1/m,
        then multiply by Mpc_m to get 1/Mpc.
        """
        return a * self.bg.get_H(a) / c * Mpc_m
    
    def estimate_alpha(self, a, delta, theta, Phi):
        """
        Estimate local acceleration magnitude α from perturbations.
        
        Paper discussion: α ≈ |a_local| ≈ (k/a) |v| or from potential gradient
        
        In Fourier space:
            v ~ θ / k (velocity from divergence)
            ∇Φ ~ k Φ (potential gradient)
        
        The local acceleration is approximately:
            α ~ k |Φ| (from potential gradient)
            or α ~ (k/a) |v| ~ |θ| / a (from velocity)
        
        We use a combination for stability.
        """
        # Velocity-based estimate
        v_est = np.abs(theta) / (self.k + 1e-10)  # Avoid division by zero
        
        # Potential gradient estimate
        grad_Phi = self.k * np.abs(Phi)
        
        # Combined estimate (use maximum for stability)
        alpha = np.maximum(v_est / (a + 1e-10), grad_Phi)
        
        # Add a floor to avoid numerical issues
        alpha_min = 1e-20  # Small but non-zero
        alpha = np.maximum(alpha, alpha_min)
        
        return alpha
    
    def derivatives(self, a, y):
        """
        Compute derivatives d(y)/da for the state vector.
        
        State: y = [δ, θ, Φ]
        
        We need to convert from d/dη to d/da:
            d/dη = aH d/da = ℋ d/da
            d/da = (1/ℋ) d/dη
        
        Equations in d/da form:
            dδ/da = -(1/ℋ)(θ + 3 dΦ/dη)
            dθ/da = -(1/ℋ)(ℋ θ + k² Φ / f)
            dΦ/da = -(1/ℋ) Φ'
        
        But we need Φ' from the Poisson equation:
            k² Φ = -(3/2) ℋ² Ω_m(a) δ a² + horizon_term
        
        So Φ' involves δ' and θ'.
        """
        delta, theta, Phi = y
        
        # Conformal Hubble parameter
        Hc = self.conformal_H(a)  # [1/Mpc]
        
        # Matter density parameter
        Om_a = self.Omega_m_a(a)
        
        # Poisson equation: solve for Φ given δ
        # k² Φ = -(3/2) ℋ² Ω_m(a) δ a² + horizon_correction
        horizon_corr = horizon_correction_term(self.k, a, self.bg) if self.hqiv_on else 0.0
        Phi_from_poisson = -(1.5 * Hc**2 * Om_a * delta * a**2 + horizon_corr * a**2) / (self.k**2 + 1e-10)
        
        # Use Poisson-computed Φ for consistency
        Phi = Phi_from_poisson
        
        # Compute Φ' (time derivative)
        # From time derivative of Poisson equation
        # This requires δ' and θ', so we iterate
        
        # Inertia reduction factor
        if self.hqiv_on:
            alpha = self.estimate_alpha(a, delta, theta, Phi)
            phi = self.bg.phi_horizon(a)
            f = inertia_reduction_factor(alpha, phi, form=self.inertia_form)
            # Clamp f to avoid numerical issues
            f = np.clip(f, 0.01, 1.0)
        else:
            f = 1.0
        
        # Derivatives (in d/da form)
        # dδ/da = -(θ + 3 Φ') / ℋ
        # dθ/da = -(ℋ θ + k² Φ / f) / ℋ
        # dΦ/da = -Φ' / ℋ
        
        # For Φ', we use the adiabatic approximation on large scales
        # and the full Poisson equation on small scales
        # Simplified: assume Φ' ≈ 0 on superhorizon scales, evolve from Poisson otherwise
        
        # dΦ/da from differentiating Poisson equation
        # This is complex, so we use a simplified approach:
        # Φ' ≈ -(1/2) ℋ Φ in matter domination (from metric evolution)
        
        # For now, use a simpler approach: evolve Φ from the Poisson constraint
        # and only evolve δ and θ dynamically
        
        # dδ/da
        # In matter domination: δ' + θ = 0 (neglecting Φ' for now)
        # More generally: δ' + θ = -3 Φ'
        
        # Estimate Φ' from the change in the Poisson source
        # Φ' ~ d/dη[-(3/2) ℋ² Ω_m δ a² / k²]
        # This is complex, so we use a simpler approximation
        
        # For the initial implementation, we use the standard approach:
        # dδ/da = -θ / ℋ (neglecting Φ' term for subhorizon modes)
        # This is valid for k >> ℋ (subhorizon)
        
        # For superhorizon modes, we need the full treatment
        # k in 1/Mpc, Hc in 1/Mpc
        k_over_Hc = self.k / (Hc + 1e-10)
        
        if k_over_Hc > 1.0:  # Subhorizon
            dPhi_da = 0.0  # Φ approximately constant in matter domination
            ddelta_da = -theta / Hc
        else:  # Superhorizon
            # On superhorizon scales, Φ is constant in standard cosmology
            # but can evolve in HQIV
            dPhi_da = 0.0
            ddelta_da = -theta / Hc - 3.0 * dPhi_da
        
        # dθ/da
        # θ' + ℋ θ = -k² Φ / f
        # dθ/da = -(ℋ θ + k² Φ / f) / ℋ = -θ - k² Φ / (f ℋ)
        dtheta_da = -theta - self.k**2 * Phi / (f * Hc)
        
        return np.array([ddelta_da, dtheta_da, dPhi_da])
    
    def initial_conditions(self, a_start):
        """
        Set initial conditions for adiabatic mode in radiation era.
        
        In the radiation era (a << a_eq), for adiabatic perturbations:
        - Superhorizon (kη << 1): Φ = const, δ = -2Φ, θ = -k²η Φ/2
        - Subhorizon: acoustic oscillations
        
        For simplicity, we start with superhorizon initial conditions
        and let the evolution handle the transition.
        
        Parameters
        ----------
        a_start : float
            Initial scale factor
            
        Returns
        -------
        y0 : array
            Initial state [δ, θ, Φ]
        """
        # Initial potential amplitude (arbitrary, will normalize later)
        Phi_0 = 1e-5
        
        # Conformal time at a_start
        eta = self.bg.get_eta(a_start) / Mpc_m  # [Mpc]
        Hc = self.conformal_H(a_start)  # [1/Mpc]
        
        # kη determines if mode is superhorizon or subhorizon
        k_eta = self.k * eta
        
        if k_eta < 0.1:  # Superhorizon
            # Standard adiabatic initial conditions
            # In radiation era: Φ = const, δ = -2Φ (for CDM)
            # But we have only baryons, so use matter-like ICs
            
            # For matter in radiation era:
            # δ grows logarithmically, but we start with constant Φ
            delta_0 = -2.0 * Phi_0  # Simplified
            theta_0 = -0.5 * self.k**2 * eta * Phi_0  # Small initial velocity
        else:
            # Subhorizon: start with small perturbations
            delta_0 = Phi_0 * np.cos(k_eta)
            theta_0 = Phi_0 * self.k * np.sin(k_eta)
            Phi_0 = Phi_0 * np.cos(k_eta) / (1.0 + 0.1 * k_eta)
        
        return np.array([delta_0, theta_0, Phi_0])
    
    def solve(self, a_start=1e-9, a_end=1.0, n_pts=500):
        """
        Solve the perturbation equations from a_start to a_end.
        
        Parameters
        ----------
        a_start : float
            Initial scale factor
        a_end : float
            Final scale factor
        n_pts : int
            Number of output points
            
        Returns
        -------
        a_arr : array
            Scale factor array
        delta_arr : array
            Density contrast evolution
        theta_arr : array
            Velocity divergence evolution
        Phi_arr : array
            Newtonian potential evolution
        """
        # Ensure background is computed
        if self.bg._a_arr is None:
            self.bg.compute_background()
        
        # Initial conditions
        y0 = self.initial_conditions(a_start)
        
        # Time points for output
        a_eval = np.logspace(np.log10(a_start), np.log10(a_end), n_pts)
        
        # Solve ODE
        # Use d(log a)/dt = H, so we integrate in log a for stability
        def deriv_log_a(log_a, y):
            a = 10.0**log_a
            dyda = self.derivatives(a, y)
            return dyda * a  # d(y)/d(log a) = a d(y)/da
        
        log_a_start = np.log10(a_start)
        log_a_end = np.log10(a_end)
        log_a_eval = np.log10(a_eval)
        
        sol = solve_ivp(
            deriv_log_a,
            [log_a_start, log_a_end],
            y0,
            method='RK45',
            t_eval=log_a_eval,
            rtol=1e-8,
            atol=1e-10
        )
        
        if not sol.success:
            print(f"Warning: ODE solver did not converge: {sol.message}")
        
        # Extract solution
        self.a_arr = 10.0**sol.t
        self.delta_arr = sol.y[0]
        self.theta_arr = sol.y[1]
        self.Phi_arr = sol.y[2]
        
        # Recompute Φ from Poisson equation for consistency
        for i, a in enumerate(self.a_arr):
            Hc = self.conformal_H(a)
            Om_a = self.Omega_m_a(a)
            horizon_corr = horizon_correction_term(self.k, a, self.bg) if self.hqiv_on else 0.0
            self.Phi_arr[i] = -(1.5 * Hc**2 * Om_a * self.delta_arr[i] * a**2 + horizon_corr * a**2) / (self.k**2 + 1e-10)
        
        return self.a_arr, self.delta_arr, self.theta_arr, self.Phi_arr
    
    def growth_factor(self):
        """
        Compute the linear growth factor D(a).
        
        D(a) = δ(a) / δ(a=1), normalized to 1 today.
        
        Returns
        -------
        a_arr : array
            Scale factor array
        D_arr : array
            Growth factor array
        """
        if self.delta_arr is None:
            raise RuntimeError("Run solve() first")
        
        D_arr = self.delta_arr / self.delta_arr[-1]
        return self.a_arr, D_arr


def compute_transfer_function(background, k_array, a_start=1e-9, a_end=1.0, 
                               hqiv_on=True, n_pts=500):
    """
    Compute transfer function T(k) for a range of wavenumbers.
    
    The transfer function relates the initial potential to the late-time
    density field:
        δ(k, a) = T(k) × Φ_initial(k)
    
    Parameters
    ----------
    background : CosmologicalBackground
        Background cosmology
    k_array : array
        Wavenumbers [1/Mpc]
    a_start : float
        Initial scale factor
    a_end : float
        Final scale factor
    hqiv_on : bool
        Include HQIV modifications
    n_pts : int
        Number of time steps
        
    Returns
    -------
    T_k : array
        Transfer function values
    delta_k : array
        Final density contrast for each k
    """
    T_k = np.zeros(len(k_array))
    delta_k = np.zeros(len(k_array))
    
    for i, k in enumerate(k_array):
        pert = ScalarPerturbations(background, k, hqiv_on=hqiv_on)
        a_arr, delta_arr, theta_arr, Phi_arr = pert.solve(a_start, a_end, n_pts)
        
        # Transfer function: ratio of final to initial amplitude
        T_k[i] = delta_arr[-1] / delta_arr[0] if delta_arr[0] != 0 else 0.0
        delta_k[i] = delta_arr[-1]
        
        if (i + 1) % 10 == 0:
            print(f"  Computed k = {k:.3f} Mpc⁻¹ ({i+1}/{len(k_array)})")
    
    return T_k, delta_k


def compute_growth_suppression(background, k_test=0.1, a_start=1e-9, a_end=1.0):
    """
    Compute the growth suppression factor in HQIV vs ΛCDM.
    
    Paper prediction: D_HQIV / D_ΛCDM ≈ 0.36 at a = 1
    
    Parameters
    ----------
    background : CosmologicalBackground
        Background cosmology (should be HQIV)
    k_test : float
        Test wavenumber [1/Mpc]
    a_start : float
        Initial scale factor
    a_end : float
        Final scale factor
        
    Returns
    -------
    suppression : float
        Ratio D_HQIV / D_ΛCDM at a = 1
    """
    # Create ΛCDM background for comparison
    lcdm_bg = CosmologicalBackground(
        H0=background.H0_km,
        Om_m=background.Om_m,
        hqiv_on=False
    )
    lcdm_bg.compute_background()
    
    # Solve for both
    pert_hqiv = ScalarPerturbations(background, k_test, hqiv_on=True)
    a_h, d_h, _, _ = pert_hqiv.solve(a_start, a_end)
    
    pert_lcdm = ScalarPerturbations(lcdm_bg, k_test, hqiv_on=False)
    a_l, d_l, _, _ = pert_lcdm.solve(a_start, a_end)
    
    # Growth factors
    D_h = d_h / d_h[-1]  # Normalized to 1 today
    D_l = d_l / d_l[-1]
    
    # Suppression at intermediate redshift (e.g., a = 0.5)
    a_test = 0.5
    D_h_test = np.interp(a_test, a_h, D_h)
    D_l_test = np.interp(a_test, a_l, D_l)
    
    suppression = D_h_test / D_l_test
    
    return suppression, a_h, D_h, a_l, D_l


# =============================================================================
# Testing and visualization
# =============================================================================

def test_scalar_evolution():
    """Test scalar perturbation evolution and create plots."""
    import matplotlib.pyplot as plt
    
    print("Testing scalar perturbation evolution...")
    
    # Create background cosmologies
    hqiv_bg = CosmologicalBackground(H0=73.2, Om_m=0.031, hqiv_on=True,
                                       beta=1.02, Om_h=1.00, n_h=1.04)
    lcdm_bg = CosmologicalBackground(H0=67.4, Om_m=0.315, hqiv_on=False)
    
    hqiv_bg.compute_background()
    lcdm_bg.compute_background()
    
    # Test wavenumber (k ~ 0.1 Mpc⁻¹, typical for large-scale structure)
    k_test = 0.1  # Mpc⁻¹
    
    # Solve perturbations
    print(f"\nSolving for k = {k_test} Mpc⁻¹...")
    
    pert_hqiv = ScalarPerturbations(hqiv_bg, k_test, hqiv_on=True)
    a_h, delta_h, theta_h, Phi_h = pert_hqiv.solve(a_start=1e-9, a_end=1.0)
    
    pert_lcdm = ScalarPerturbations(lcdm_bg, k_test, hqiv_on=False)
    a_l, delta_l, theta_l, Phi_l = pert_lcdm.solve(a_start=1e-9, a_end=1.0)
    
    # Growth factors
    D_h = delta_h / delta_h[-1]
    D_l = delta_l / delta_l[-1]
    
    # Create plots
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    
    # Density contrast evolution
    ax = axes[0, 0]
    ax.loglog(a_h, np.abs(delta_h), 'b-', label='HQIV', linewidth=2)
    ax.loglog(a_l, np.abs(delta_l), 'r--', label='ΛCDM', linewidth=2)
    ax.set_xlabel('Scale factor a')
    ax.set_ylabel('|δ|')
    ax.set_title(f'Density Contrast Evolution (k = {k_test} Mpc⁻¹)')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # Growth factor comparison
    ax = axes[0, 1]
    ax.semilogx(a_h, D_h, 'b-', label='HQIV', linewidth=2)
    ax.semilogx(a_l, D_l, 'r--', label='ΛCDM', linewidth=2)
    ax.set_xlabel('Scale factor a')
    ax.set_ylabel('D(a)')
    ax.set_title('Linear Growth Factor')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # Potential evolution
    ax = axes[1, 0]
    ax.loglog(a_h, np.abs(Phi_h), 'b-', label='HQIV', linewidth=2)
    ax.loglog(a_l, np.abs(Phi_l), 'r--', label='ΛCDM', linewidth=2)
    ax.set_xlabel('Scale factor a')
    ax.set_ylabel('|Φ|')
    ax.set_title('Newtonian Potential Evolution')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # Growth suppression
    ax = axes[1, 1]
    suppression = D_h / np.interp(a_h, a_l, D_l)
    ax.semilogx(a_h, suppression, 'b-', linewidth=2)
    ax.axhline(1.0, color='k', linestyle='--', alpha=0.5)
    ax.axhline(0.36, color='g', linestyle=':', alpha=0.7, label='Paper prediction (0.36)')
    ax.set_xlabel('Scale factor a')
    ax.set_ylabel('D_HQIV / D_ΛCDM')
    ax.set_title('Growth Suppression Factor')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('hqiv_solver/scalar_evolution.png', dpi=150)
    plt.close()
    print("\nSaved scalar_evolution.png")
    
    # Print summary
    print("\n" + "=" * 60)
    print("Scalar Perturbation Summary")
    print("=" * 60)
    print(f"Final δ (HQIV): {delta_h[-1]:.4f}")
    print(f"Final δ (ΛCDM): {delta_l[-1]:.4f}")
    print(f"Growth suppression at a=0.5: {np.interp(0.5, a_h, suppression):.4f}")
    print("=" * 60)
    
    return a_h, delta_h, a_l, delta_l


if __name__ == "__main__":
    test_scalar_evolution()