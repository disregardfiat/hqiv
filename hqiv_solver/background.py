"""
HQIV Cosmology — Background Evolution Module.

This module implements the background cosmology for both standard ΛCDM
and the Horizon-Quantized Informational Vacuum (HQIV) framework.

Paper Reference: Section 5 "Background Dynamics"
- Eq. 13: Modified acceleration equation with horizon term βH²
- Eq. 14: Modified Friedmann constraint with horizon density

Author: HQIV Team
"""

import numpy as np
from scipy.integrate import solve_ivp, quad
from scipy.interpolate import interp1d

# =============================================================================
# Physical Constants (SI units)
# =============================================================================
c = 2.99792458e8       # Speed of light [m/s]
G0_SI = 6.67430e-11    # Gravitational constant [m³ kg⁻¹ s⁻²]
Mpc_m = 3.0856775814913673e22  # Megaparsec in meters
Gyr_s = 3.1536e16      # Gigayear in seconds


def hqiv_age(gamma=0.40, Omega_m0=0.06, h=0.732, Neff=3.046):
    """
    Universe age in Gyr for HQIV background.

    Modified Friedmann in dimensionless form:
        h² - (γ/3)h = Ω(a)
    where h = H/H0.  At a = 1 the closure condition h(1) = 1 requires
        Ω_total(1) = 1 − γ/3,
    so an effective constant-density term
        Ω_eff = 1 − γ/3 − Ω_m − Ω_r
    appears alongside matter and radiation.  This is not "dark energy"
    but the geometric consequence of the γH horizon term; omitting it
    gives H(a=1) ≈ 0.32 H0 instead of H0.
    """
    Omega_r0 = 2.47e-5 / h**2 * (1 + 0.2271 * Neff)  # photons + neutrinos
    omega_eff = 1.0 - gamma / 3.0 - Omega_m0 - Omega_r0

    def E(a):
        """Dimensionless Hubble E(a) = H(a)/H0 from 3h² - γh = 3Ω(a)."""
        rho_tot = Omega_m0 / a**3 + Omega_r0 / a**4 + omega_eff
        disc = gamma**2 + 36.0 * rho_tot
        return (gamma + np.sqrt(disc)) / 6.0

    def integrand(a):
        return 1.0 / (a * E(a))

    integral, _ = quad(integrand, 1e-10, 1.0, epsabs=1e-10, epsrel=1e-10, limit=1000)
    H0_km_s_Mpc = h * 100.0
    age_Gyr = integral * (977.8 / H0_km_s_Mpc)
    return age_Gyr


# Planck units
hbar = 1.054571817e-34  # Reduced Planck constant [J·s]
k_B = 1.380649e-23      # Boltzmann constant [J/K]
T_CMB = 2.7255          # CMB temperature today [K]


class CosmologicalBackground:
    """
    Background cosmology solver for ΛCDM and HQIV frameworks.
    
    This class computes the background expansion history H(a), conformal time η(a),
    and related quantities needed for perturbation evolution.
    
    Parameters
    ----------
    H0 : float
        Hubble constant today [km/s/Mpc]
    Om_m : float
        Matter density parameter today (Ω_m)
    Om_r : float, optional
        Radiation density parameter today (Ω_r). If None, computed from T_CMB.
    Om_Lambda : float, optional
        Dark energy density parameter (ΛCDM only). If None, set to 1 - Om_m - Om_r.
    hqiv_on : bool
        If True, use HQIV modified expansion. If False, use standard ΛCDM.
    beta : float
        HQIV horizon-smoothing parameter (paper eq. 1)
    Om_h : float
        HQIV effective horizon density parameter (Ω_horizon)
    n_h : float
        HQIV horizon dilution exponent (paper eq. 15)
    alpha_G : float
        Exponent for varying G(a) (paper eq. 6)
    
    Notes
    -----
    The HQIV modified Friedmann equation (paper eq. 15):
        H² = H0² [Ω_m a⁻³ + Ω_r a⁻⁴ + Ω_h a⁻ⁿ]
    
    The horizon term provides late-time acceleration without dark energy.
    """
    
    def __init__(self, H0=73.2, Om_m=0.031, Om_r=None, Om_Lambda=None,
                 hqiv_on=True, beta=1.02, Om_h=1.00, n_h=1.04, alpha_G=0.6, gamma=0.40):
        # Store parameters
        self.H0_km = H0  # km/s/Mpc
        self.H0 = H0 * 1000.0 / Mpc_m  # Convert to 1/s
        self.Om_m = Om_m
        self.hqiv_on = hqiv_on
        self.beta = beta
        self.Om_h = Om_h
        self.n_h = n_h
        self.alpha_G = alpha_G
        self.gamma = gamma  # Thermodynamic coefficient from Brodie
        
        # Compute radiation density from CMB temperature if not provided
        if Om_r is None:
            # ρ_γ = π²/15 × (k_B T)⁴ / (ℏ c)³
            # This gives energy density, convert to mass density: ρ = ε/c²
            epsilon_gamma = (np.pi**2 / 15.0) * (k_B * T_CMB)**4 / (hbar * c)**3
            rho_gamma = epsilon_gamma / c**2  # Convert energy density to mass density
            
            # Critical density today: ρ_c = 3H0² / (8πG)
            rho_c = 3.0 * self.H0**2 / (8.0 * np.pi * G0_SI)
            
            self.Om_r = rho_gamma / rho_c
            # Add neutrino contribution (N_eff ≈ 3.046)
            # Neutrinos contribute ~0.2271 times photon density per species
            self.Om_r *= (1.0 + 3.046 * 7.0/8.0 * (4.0/11.0)**(4.0/3.0))
        else:
            self.Om_r = Om_r
        
        # ΛCDM dark energy (only used when hqiv_on=False)
        if Om_Lambda is None:
            self.Om_Lambda = 1.0 - self.Om_m - self.Om_r
        else:
            self.Om_Lambda = Om_Lambda
        
        # Storage for computed quantities
        self._a_arr = None
        self._H_arr = None
        self._eta_arr = None
        self._t_arr = None
        self._H_interp = None
        self._eta_interp = None
        
    def H(self, a):
        """
        Hubble parameter H(a) in units of 1/s.
        
        From Brodie's thermodynamics with γ parameter (paper: φ/c² for curvature):
            G_μν + γ(φ/c²) g_μν = (8π G_eff / c⁴) T_μν

        Numerics use c = ℏ = 1 (φ = H). Modified Friedmann:
            3H² - γ H - 8πG ρ = 0
        
        For normalization H(a=1) = H0, we need:
            ρ(a=1) = (3H0² - γH0) / (8πG) = ρ_crit(1 - γ/3)
        
        In dimensionless form (h = H/H₀):
            3h² - γh - 3Ω(a) = 0
        
        Solving:
            h(a) = [γ + √(γ² + 36Ω(a))] / 6
        
        where Ω(a) includes an effective Ω_Λ term to satisfy normalization.
        
        Parameters
        ----------
        a : float or array
            Scale factor
            
        Returns
        -------
        H : float or array
            Hubble parameter in 1/s
        """
        a = np.asarray(a)
        
        if self.hqiv_on:
            # γ from Brodie's thermodynamics (default 0.55)
            gamma = getattr(self, 'gamma', 0.40)
            
            # To satisfy H(a=1) = H0:
            # At a=1: 3(1)² - γ(1) - 3Ω(1) = 0 → Ω(1) = 1 - γ/3
            # With Ω_m + Ω_r = ~0.031 + 8e-5 ≈ 0.031
            # We need Ω_Λ_eff = (1 - γ/3) - Ω_m - Ω_r
            
            omega_total = self.Om_m + self.Om_r
            omega_eff = 1.0 - gamma/3.0 - omega_total
            
            # Effective density: Ω_m a⁻³ + Ω_r a⁻⁴ + Ω_eff
            rho_norm = self.Om_m * a**(-3) + self.Om_r * a**(-4) + omega_eff
            
            # Solve: h = [γ + √(γ² + 36ρ)] / 6
            disc = gamma**2 + 36.0 * rho_norm
            h = (gamma + np.sqrt(disc)) / 6.0
            
            H = h * self.H0
        else:
            # Standard ΛCDM
            matter_term = self.Om_m * a**(-3)
            radiation_term = self.Om_r * a**(-4)
            lambda_term = self.Om_Lambda
            H2 = self.H0**2 * (matter_term + radiation_term + lambda_term)
            H = np.sqrt(H2)
        
        return H
    
    def H_km(self, a):
        """Hubble parameter in km/s/Mpc."""
        return self.H(a) * Mpc_m / 1000.0
    
    def H_over_H0(self, a):
        """Hubble parameter in units of H0 (E(a) = H(a)/H0)."""
        return self.H(a) / self.H0
    
    def Hc(self, a):
        """
        Conformal Hubble parameter ℋ = aH in units of 1/s.
        
        This is used extensively in the perturbation equations.
        """
        return a * self.H(a)
    
    def G_eff(self, a):
        """
        Effective gravitational coupling G_eff(a).
        
        Paper eq. 6:
            G(a) = G0 × (H(a)/H0)^α
        
        This varying G is a key prediction of HQIV, derived from
        horizon-scale modifications to the gravitational interaction.
        
        Parameters
        ----------
        a : float or array
            Scale factor
            
        Returns
        -------
        G_eff : float or array
            Effective gravitational constant in units of G0
        """
        if self.hqiv_on:
            return (self.H(a) / self.H0)**self.alpha_G
        else:
            return np.ones_like(a) if hasattr(a, '__len__') else 1.0
    
    def phi_horizon(self, a):
        """
        Geometric horizon field φ(a).
        
        From the action (β-free formulation):
            φ ≡ 2c² / Θ_local
        
        In homogeneous FLRW background (using horizon radius):
            Θ_local = 2c / H  →  φ = cH
        
        Note: This is different from θ = 3H (expansion scalar).
        The thermodynamic/informational term uses horizon radius, not expansion.
        
        Parameters
        ----------
        a : float or array
            Scale factor
            
        Returns
        -------
        phi : float or array
            Horizon field in units of 1/s
        """
        if self.hqiv_on:
            # Paper-correct: φ = cH (not 3cH)
            return c * self.H(a)
        else:
            # Standard cosmology: no horizon field
            return np.zeros_like(a) if hasattr(a, '__len__') else 0.0
    
    def phi_conf(self, a):
        """
        Conformal horizon field φ_conf = a × φ = c × aH = c × ℋ.
        
        This is the conformal-time version used in perturbation equations.
        """
        if self.hqiv_on:
            return c * a * self.H(a)
        else:
            return np.zeros_like(a) if hasattr(a, '__len__') else 0.0
    
    def beta_a(self, a):
        """
        Time-dependent horizon-smoothing parameter β(a).
        
        For β-free formulation, we use a simple evolution that gives
        β → 1 as a → 1 (today).
        
        The evolution follows from structure formation:
            β(a) = β_early + (β_late - β_early) × [1 - exp(-(a/a_eq)^γ)]
        
        where:
            β_early ≈ 0.75 (early universe, high anisotropy)
            β_late = self.beta (default 1.02)
            a_eq ~ 0.3 (structure formation transition)
        """
        if self.hqiv_on:
            # Early universe: more horizon anisotropy (β ~ 0.7-0.8)
            beta_early = 0.75
            beta_late = self.beta  # Late-time value (default 1.02)
            
            # Transition scale: around matter-radiation equality to structure formation
            a_transition = 0.3
            gamma = 1.0
            
            # Evolution function
            evolution = 1.0 - np.exp(-(a / a_transition)**gamma)
            
            return beta_early + (beta_late - beta_early) * evolution
        else:
            return np.ones_like(a) if hasattr(a, '__len__') else 1.0
    
    def compute_background(self, a_start=1e-9, a_end=1.0, n_pts=2000):
        """
        Compute and store background quantities on a grid.
        
        This precomputes H(a), η(a), t(a) for use in perturbation evolution.
        
        Parameters
        ----------
        a_start : float
            Initial scale factor (default: a ≈ 10⁻⁹, deep radiation era)
        a_end : float
            Final scale factor (default: a = 1, today)
        n_pts : int
            Number of grid points
        """
        # Use logarithmic spacing for better resolution at early times
        self._a_arr = np.logspace(np.log10(a_start), np.log10(a_end), n_pts)
        
        # Compute H(a)
        self._H_arr = self.H(self._a_arr)
        
        # Compute conformal time η = ∫ c da / (a² H)
        # In log space: d(ln a) = da/a, so da = a d(ln a)
        # η = ∫ c d(ln a) / (a H)
        log_a = np.log(self._a_arr)
        integrand_eta = c / (self._a_arr * self._H_arr)  # c / (a H)
        self._eta_arr = np.zeros_like(self._a_arr)
        for i in range(1, len(self._a_arr)):
            dlog_a = log_a[i] - log_a[i-1]
            avg_integrand = 0.5 * (integrand_eta[i] + integrand_eta[i-1])
            self._eta_arr[i] = self._eta_arr[i-1] + avg_integrand * dlog_a
        
        # Compute proper time t = ∫ da / (a H)
        # In log space: t = ∫ d(ln a) / H
        # But we need t = ∫ da/(aH) = ∫ d(ln a) / H
        integrand_t = 1.0 / self._H_arr  # 1/H for d(ln a) integration
        self._t_arr = np.zeros_like(self._a_arr)
        for i in range(1, len(self._a_arr)):
            dlog_a = log_a[i] - log_a[i-1]
            avg_integrand = 0.5 * (integrand_t[i] + integrand_t[i-1])
            self._t_arr[i] = self._t_arr[i-1] + avg_integrand * dlog_a
        
        # Create interpolation functions
        self._H_interp = interp1d(self._a_arr, self._H_arr, kind='cubic', fill_value='extrapolate')
        self._eta_interp = interp1d(self._a_arr, self._eta_arr, kind='cubic', fill_value='extrapolate')
        self._t_interp = interp1d(self._a_arr, self._t_arr, kind='cubic', fill_value='extrapolate')
        
        return self._a_arr, self._H_arr, self._eta_arr
    
    def get_H(self, a):
        """Interpolated H(a) in 1/s."""
        if self._H_interp is None:
            self.compute_background()
        return float(self._H_interp(a))
    
    def get_eta(self, a):
        """Interpolated conformal time η(a) in meters."""
        if self._eta_interp is None:
            self.compute_background()
        return float(self._eta_interp(a))
    
    def get_t(self, a):
        """Interpolated proper time t(a) in seconds."""
        if self._t_interp is None:
            self.compute_background()
        return float(self._t_interp(a))
    
    def age_today(self):
        """Universe age today in Gyr. Uses CLASS-consistent hqiv_age when HQIV is on."""
        if self.hqiv_on:
            h = self.H0_km / 100.0
            return hqiv_age(gamma=self.gamma, Omega_m0=self.Om_m, h=h)
        if self._t_arr is None:
            self.compute_background()
        return self._t_arr[-1] / Gyr_s
    
    def proper_time_at_z(self, z):
        """Proper time from big bang to redshift z in Gyr."""
        a = 1.0 / (1.0 + z)
        return self.get_t(a) / Gyr_s
    
    def conformal_time_today(self):
        """Conformal time today η₀ in Mpc."""
        if self._eta_arr is None:
            self.compute_background()
        return self._eta_arr[-1] / Mpc_m
    
    def sound_horizon(self, a_rec=1/1100):
        """
        Sound horizon at recombination in Mpc.
        
        r_s = ∫ c_s da / (a² H) from a=0 to a_rec
        
        where c_s ≈ c/√3 is the sound speed in the photon-baryon fluid.
        """
        if self._a_arr is None:
            self.compute_background()
        
        mask = self._a_arr <= a_rec
        a_slice = self._a_arr[mask]
        H_slice = self._H_arr[mask]
        
        # Sound speed in photon-baryon fluid (simplified)
        c_s = c / np.sqrt(3.0)
        
        integrand = c_s / (a_slice**2 * H_slice)
        r_s = 0.0
        for i in range(1, len(a_slice)):
            da = a_slice[i] - a_slice[i-1]
            r_s += 0.5 * (integrand[i] + integrand[i-1]) * da
        
        return r_s / Mpc_m
    
    def summary(self):
        """Print summary of background cosmology."""
        if self._a_arr is None:
            self.compute_background()
        
        print("=" * 60)
        print("HQIV Background Cosmology Summary" if self.hqiv_on else "ΛCDM Background Cosmology Summary")
        print("=" * 60)
        print(f"H0 = {self.H0_km:.1f} km/s/Mpc")
        print(f"Ω_m = {self.Om_m:.4f}")
        print(f"Ω_r = {self.Om_r:.2e}")
        if self.hqiv_on:
            print(f"Ω_horizon = {self.Om_h:.4f}")
            print(f"n_horizon = {self.n_h:.4f}")
            print(f"β = {self.beta:.4f}")
            print(f"α_G = {self.alpha_G:.4f}")
        else:
            print(f"Ω_Λ = {self.Om_Lambda:.4f}")
        print("-" * 60)
        print(f"Universe age: {self.age_today():.2f} Gyr")
        print(f"Conformal time today: {self.conformal_time_today():.2f} Mpc")
        print(f"Sound horizon at recombination: {self.sound_horizon():.2f} Mpc")
        print(f"Proper time at z=10: {self.proper_time_at_z(10)*1000:.0f} Myr")
        print(f"Proper time at z=14: {self.proper_time_at_z(14)*1000:.0f} Myr")
        print("=" * 60)


# =============================================================================
# Convenience functions for quick testing
# =============================================================================

def compare_cosmologies():
    """Compare HQIV and ΛCDM background evolution."""
    import matplotlib.pyplot as plt
    
    # Create both cosmologies
    hqiv = CosmologicalBackground(H0=73.2, Om_m=0.031, hqiv_on=True, 
                                   beta=1.02, Om_h=1.00, n_h=1.04)
    lcdm = CosmologicalBackground(H0=67.4, Om_m=0.315, hqiv_on=False)
    
    # Compute backgrounds
    hqiv.compute_background()
    lcdm.compute_background()
    
    # Print summaries
    hqiv.summary()
    print()
    lcdm.summary()
    
    # Create comparison plots
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    
    a = hqiv._a_arr
    
    # H(a) comparison
    ax = axes[0, 0]
    ax.loglog(a, hqiv.H_km(a), 'b-', label='HQIV', linewidth=2)
    ax.loglog(a, lcdm.H_km(a), 'r--', label='ΛCDM', linewidth=2)
    ax.set_xlabel('Scale factor a')
    ax.set_ylabel('H(a) [km/s/Mpc]')
    ax.set_title('Hubble Parameter Evolution')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # G_eff(a) for HQIV
    ax = axes[0, 1]
    ax.semilogx(a, hqiv.G_eff(a), 'b-', linewidth=2)
    ax.axhline(1.0, color='k', linestyle='--', alpha=0.5)
    ax.set_xlabel('Scale factor a')
    ax.set_ylabel('G_eff / G0')
    ax.set_title('Effective Gravitational Coupling (HQIV)')
    ax.grid(True, alpha=0.3)
    
    # β(a) evolution
    ax = axes[1, 0]
    ax.semilogx(a, hqiv.beta_a(a), 'b-', linewidth=2)
    ax.axhline(1.0, color='k', linestyle='--', alpha=0.5)
    ax.set_xlabel('Scale factor a')
    ax.set_ylabel('β(a)')
    ax.set_title('Horizon-Smoothing Parameter Evolution')
    ax.grid(True, alpha=0.3)
    
    # Proper time comparison
    ax = axes[1, 1]
    ax.semilogx(a, hqiv._t_arr / Gyr_s, 'b-', label='HQIV', linewidth=2)
    ax.semilogx(a, lcdm._t_arr / Gyr_s, 'r--', label='ΛCDM', linewidth=2)
    ax.set_xlabel('Scale factor a')
    ax.set_ylabel('Proper time t [Gyr]')
    ax.set_title('Cosmic Time Evolution')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('hqiv_solver/background_comparison.png', dpi=150)
    plt.close()
    print("\nSaved background_comparison.png")
    
    return hqiv, lcdm


if __name__ == "__main__":
    compare_cosmologies()