"""
HQIV Background Cosmology Wrapper
=================================

Provides a clean interface to the HQIV background solver for use in N-body simulations.
Uses the exact H(z) from the full HQIV modified Friedmann equation.

Paper Reference: paper/main.tex, Section 5
    3H² - γH = 8πG_eff ρ(a)

Author: HQIV Team
"""

import sys
from pathlib import Path
import numpy as np
from scipy.interpolate import interp1d

# Add hqiv_solver to path
hqiv_solver_path = Path(__file__).parent.parent.parent / 'hqiv_solver'
if str(hqiv_solver_path) not in sys.path:
    sys.path.insert(0, str(hqiv_solver_path))

try:
    from background import CosmologicalBackground
    _HAS_HQIV_SOLVER = True
except ImportError:
    _HAS_HQIV_SOLVER = False
    print("Warning: hqiv_solver.background not found. Using approximate H(a).")


class HQIVBackground:
    """
    Wrapper for HQIV background cosmology.
    
    This class provides a clean interface to the exact HQIV background solver,
    with caching and interpolation for efficient use in N-body simulations.
    
    Parameters
    ----------
    H0 : float
        Hubble constant today [km/s/Mpc]
    Omega_m : float
        Matter density parameter (baryons only in HQIV)
    gamma : float
        HQIV thermodynamic coefficient from Brodie's overlap integral
    alpha_G : float
        Exponent for varying G(a)
    Neff : float
        Effective neutrino species
    """
    
    def __init__(self, H0=73.2, Omega_m=0.06, gamma=0.40, alpha_G=0.6, Neff=3.046):
        self.H0 = H0
        self.Omega_m = Omega_m
        self.gamma = gamma
        self.alpha_G = alpha_G
        self.Neff = Neff
        
        # Storage for interpolation
        self._a_table = None
        self._H_table = None
        self._H_interp = None
        self._G_eff_table = None
        self._G_eff_interp = None
        
        if _HAS_HQIV_SOLVER:
            self._init_solver()
        else:
            self._init_approximate()
    
    def _init_solver(self):
        """Initialize using the full HQIV solver."""
        self._solver = CosmologicalBackground(
            H0=self.H0,
            Om_m=self.Omega_m,
            hqiv_on=True,
            gamma=self.gamma,
            alpha_G=self.alpha_G
        )
        self._solver.compute_background()
        
        self._a_table = self._solver._a_arr
        self._H_table = self._solver._H_arr
        self._G_eff_table = self._solver.G_eff(self._a_table)
        
        # Create interpolators
        self._H_interp = interp1d(
            self._a_table, self._H_table, 
            kind='cubic', fill_value='extrapolate'
        )
        self._G_eff_interp = interp1d(
            self._a_table, self._G_eff_table,
            kind='cubic', fill_value='extrapolate'
        )
        
        self._has_exact_solver = True
    
    def _init_approximate(self):
        """Initialize using approximate formulas (fallback)."""
        # Create a table using the approximate formula
        self._a_table = np.logspace(-4, 0, 1000)
        
        # Compute radiation density
        h = self.H0 / 100.0
        Omega_r = 2.47e-5 / h**2 * (1 + 0.2271 * self.Neff)
        
        # Approximate H(a) from modified Friedmann: 3H² - γH = 8πG ρ
        # In dimensionless form: h = (γ + sqrt(γ² + 36Ω(a))) / 6
        # where Ω(a) = Ω_m a⁻³ + Ω_r a⁻⁴ + Ω_Λ_eff
        # and Ω_Λ_eff = 1 - γ/3 - Ω_m - Ω_r (to satisfy h(a=1) = 1)
        
        omega_eff = 1.0 - self.gamma/3.0 - self.Omega_m - Omega_r
        
        rho_norm = (self.Omega_m * self._a_table**(-3) + 
                   Omega_r * self._a_table**(-4) + 
                   omega_eff)
        
        disc = self.gamma**2 + 36.0 * rho_norm
        h_dimless = (self.gamma + np.sqrt(disc)) / 6.0
        
        # Convert to physical H
        H0_SI = self.H0 * 1000.0 / 3.0856775814913673e22  # 1/s
        self._H_table = h_dimless * H0_SI
        
        # G_eff(a) = (H(a)/H0)^alpha_G
        self._G_eff_table = h_dimless ** self.alpha_G
        
        # Create interpolators
        self._H_interp = interp1d(
            self._a_table, self._H_table,
            kind='cubic', fill_value='extrapolate'
        )
        self._G_eff_interp = interp1d(
            self._a_table, self._G_eff_table,
            kind='cubic', fill_value='extrapolate'
        )
        
        self._has_exact_solver = False
    
    def H_of_a(self, a):
        """
        Get Hubble parameter at scale factor a.
        
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
        result = self._H_interp(a)
        if result.ndim == 0:
            return float(result)
        return result
    
    def H_of_z(self, z):
        """
        Get Hubble parameter at redshift z.
        
        Parameters
        ----------
        z : float or array
            Redshift
            
        Returns
        -------
        H : float or array
            Hubble parameter in 1/s
        """
        z = np.asarray(z)
        a = 1.0 / (1.0 + z)
        return self.H_of_a(a)
    
    def H_km_of_a(self, a):
        """
        Get Hubble parameter in km/s/Mpc.
        
        Parameters
        ----------
        a : float or array
            Scale factor
            
        Returns
        -------
        H : float or array
            Hubble parameter in km/s/Mpc
        """
        Mpc_m = 3.0856775814913673e22
        H = self.H_of_a(a)
        return H * Mpc_m / 1000.0
    
    def H_km_of_z(self, z):
        """Get Hubble parameter at redshift z in km/s/Mpc."""
        z = np.asarray(z)
        a = 1.0 / (1.0 + z)
        return self.H_km_of_a(a)
    
    def G_eff_of_a(self, a):
        """
        Get effective gravitational constant at scale factor a.
        
        Parameters
        ----------
        a : float or array
            Scale factor
            
        Returns
        -------
        G_eff : float or array
            Effective G in units of G0
        """
        a = np.asarray(a)
        result = self._G_eff_interp(a)
        if result.ndim == 0:
            return float(result)
        return result
    
    def G_eff_of_z(self, z):
        """Get effective gravitational constant at redshift z."""
        z = np.asarray(z)
        a = 1.0 / (1.0 + z)
        return self.G_eff_of_a(a)
    
    def age_today(self):
        """Get universe age today in Gyr from t = ∫ da/(a H(a)), consistent with our H(a)."""
        Gyr_s = 3600 * 24 * 365.25 * 1e9
        if self._a_table is None or self._H_table is None:
            self._init_solver() if _HAS_HQIV_SOLVER else self._init_approximate()
        # _H_table is in 1/s; age = ∫ da/(a H) = ∫ d(ln a)/H
        log_a = np.log(self._a_table)
        integrand = 1.0 / self._H_table
        age_s = np.trapz(integrand, log_a)
        return age_s / Gyr_s
    
    def summary(self):
        """Print summary of background cosmology."""
        print("=" * 60)
        print("HQIV Background Cosmology")
        print("=" * 60)
        print(f"H0 = {self.H0:.1f} km/s/Mpc")
        print(f"Ω_m = {self.Omega_m:.4f}")
        print(f"γ = {self.gamma:.4f}")
        print(f"α_G = {self.alpha_G:.4f}")
        print(f"Using exact solver: {self._has_exact_solver}")
        print("-" * 60)
        print(f"H(z=0) = {self.H_km_of_z(0):.2f} km/s/Mpc")
        print(f"H(z=1) = {self.H_km_of_z(1):.2f} km/s/Mpc")
        print(f"H(z=10) = {self.H_km_of_z(10):.2f} km/s/Mpc")
        print(f"G_eff(z=0) = {self.G_eff_of_z(0):.4f}")
        print(f"G_eff(z=10) = {self.G_eff_of_z(10):.4f}")
        print(f"Universe age: {self.age_today():.2f} Gyr")
        print("=" * 60)


# =============================================================================
# Convenience functions for direct use
# =============================================================================

def get_H_at_redshift(z, gamma=0.40, Omega_m=0.06, H0=73.2, **kwargs):
    """
    Returns the exact H(z) from the full HQIV modified Friedmann equation.
    
    This is the main convenience function for getting H at any redshift.
    
    Parameters
    ----------
    z : float or array
        Redshift
    gamma : float
        HQIV thermodynamic coefficient (default 0.40)
    Omega_m : float
        Matter density parameter (default 0.06)
    H0 : float
        Hubble constant today in km/s/Mpc (default 73.2)
    **kwargs : dict
        Additional parameters (alpha_G, Neff, etc.)
        
    Returns
    -------
    H : float or array
        Hubble parameter in km/s/Mpc
        
    Examples
    --------
    >>> H = get_H_at_redshift(0.3)
    >>> print(f"H(z=0.3) = {H:.2f} km/s/Mpc")
    
    >>> # Array of redshifts
    >>> z = np.array([0, 0.5, 1, 2, 5])
    >>> H = get_H_at_redshift(z)
    """
    bg = HQIVBackground(H0=H0, Omega_m=Omega_m, gamma=gamma, **kwargs)
    return bg.H_km_of_z(z)


def get_H_at_scale_factor(a, gamma=0.40, Omega_m=0.06, H0=73.2, **kwargs):
    """
    Returns the exact H(a) from the full HQIV modified Friedmann equation.
    
    Parameters
    ----------
    a : float or array
        Scale factor
    gamma : float
        HQIV thermodynamic coefficient (default 0.40)
    Omega_m : float
        Matter density parameter (default 0.06)
    H0 : float
        Hubble constant today in km/s/Mpc (default 73.2)
    **kwargs : dict
        Additional parameters
        
    Returns
    -------
    H : float or array
        Hubble parameter in km/s/Mpc
    """
    bg = HQIVBackground(H0=H0, Omega_m=Omega_m, gamma=gamma, **kwargs)
    return bg.H_km_of_a(a)


def get_G_eff_at_redshift(z, gamma=0.40, Omega_m=0.06, H0=73.2, alpha_G=0.6, **kwargs):
    """
    Returns G_eff(z) from the HQIV varying gravitational coupling.
    
    Paper eq. 6: G(a) = G0 × (H(a)/H0)^α
    
    Parameters
    ----------
    z : float or array
        Redshift
    gamma : float
        HQIV thermodynamic coefficient
    Omega_m : float
        Matter density parameter
    H0 : float
        Hubble constant today in km/s/Mpc
    alpha_G : float
        Exponent for varying G (default 0.6)
    **kwargs : dict
        Additional parameters
        
    Returns
    -------
    G_eff : float or array
        Effective gravitational constant in units of G0
    """
    bg = HQIVBackground(H0=H0, Omega_m=Omega_m, gamma=gamma, alpha_G=alpha_G, **kwargs)
    return bg.G_eff_of_z(z)


# =============================================================================
# Module-level default background (for efficiency)
# =============================================================================

_default_background = None

def get_default_background(H0=73.2, Omega_m=0.06, gamma=0.40, alpha_G=0.6):
    """
    Get or create the default background cosmology.
    
    This caches the background solver for efficiency.
    
    Parameters
    ----------
    H0 : float
        Hubble constant [km/s/Mpc]
    Omega_m : float
        Matter density parameter
    gamma : float
        HQIV thermodynamic coefficient
    alpha_G : float
        Exponent for varying G
        
    Returns
    -------
    bg : HQIVBackground
        Background cosmology object
    """
    global _default_background
    
    if _default_background is None:
        _default_background = HQIVBackground(
            H0=H0, Omega_m=Omega_m, gamma=gamma, alpha_G=alpha_G
        )
    
    return _default_background


# =============================================================================
# Testing
# =============================================================================

def test_background_wrapper():
    """Test the background wrapper."""
    print("Testing HQIV Background Wrapper...")
    
    # Create background
    bg = HQIVBackground(H0=73.2, Omega_m=0.06, gamma=0.40, alpha_G=0.6)
    
    # Print summary
    bg.summary()
    
    # Test convenience functions
    print("\nConvenience function tests:")
    print(f"  H(z=0) = {get_H_at_redshift(0):.2f} km/s/Mpc")
    print(f"  H(z=0.3) = {get_H_at_redshift(0.3):.2f} km/s/Mpc")
    print(f"  H(z=1) = {get_H_at_redshift(1):.2f} km/s/Mpc")
    print(f"  G_eff(z=0) = {get_G_eff_at_redshift(0):.4f}")
    print(f"  G_eff(z=10) = {get_G_eff_at_redshift(10):.4f}")
    
    return bg


if __name__ == "__main__":
    test_background_wrapper()