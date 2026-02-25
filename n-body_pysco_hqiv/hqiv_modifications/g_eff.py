"""
HQIV Effective Gravitational Constant Module
=============================================

Implements the varying gravitational coupling G_eff(φ) and horizon correction terms.

Paper Reference: paper/main.tex, Section 3
    G(a) = G0 × (H(a)/H0)^α

where α ≈ 0.6 is derived from QI galaxy rotation curve matching.

Also implements the horizon correction term in the Poisson equation:
    ∇²Φ = 4πG_eff(φ) ρ_m δ + horizon_correction

Author: HQIV Team
"""

import numpy as np
from numpy.fft import fftn, ifftn, fftfreq
import configparser
from pathlib import Path

# Physical constants (SI units)
c = 2.99792458e8       # Speed of light [m/s]
G0_SI = 6.67430e-11    # Gravitational constant [m³ kg⁻¹ s⁻²]
Mpc_m = 3.0856775814913673e22  # Megaparsec in meters


def load_config():
    """Load configuration from config_hqiv.ini."""
    config_path = Path(__file__).parent / 'config_hqiv.ini'
    config = configparser.ConfigParser()
    config.read(config_path)
    return config


# Load default parameters
_config = load_config()
_alpha_G = float(_config.get('hqiv_parameters', 'alpha_G', fallback='0.6'))
_gamma = float(_config.get('hqiv_parameters', 'gamma', fallback='0.40'))
_Omega_m = float(_config.get('cosmology', 'Omega_m', fallback='0.06'))
_H0 = float(_config.get('cosmology', 'H0', fallback='73.2'))


def effective_gravitational_constant(H_a, H0, alpha=_alpha_G):
    """
    Compute the effective gravitational coupling G_eff(a).
    
    Paper Reference: paper/main.tex, Section 3
        G(a) = G0 × (H(a)/H0)^α
    
    This varying G is a key prediction of HQIV, derived from
    horizon-scale modifications to the gravitational interaction.
    
    Parameters
    ----------
    H_a : float or array
        Hubble parameter at scale factor a [1/s]
    H0 : float
        Hubble constant today [1/s]
    alpha : float
        Exponent for varying G (default 0.6)
        
    Returns
    -------
    G_eff : float or array
        Effective gravitational constant in units of G0
        
    Notes
    -----
    At early times (high z), H(a) > H0, so G_eff > G0.
    At late times (low z), H(a) < H0, so G_eff < G0.
    
    This affects:
    - Structure growth rate
    - CMB acoustic peak positions
    - Stellar evolution timescales
    """
    return (H_a / H0) ** alpha


def horizon_correction_term(k, a, H_a, H0, Omega_m, gamma=_gamma, box_size=100.0, N_grid=256):
    """
    Compute the horizon correction term in the Poisson equation.
    
    Paper Reference: paper/main.tex, Eq. 11
        ∇²Φ = 4πG_eff(φ) ρ_m δ + horizon_correction
    
    The horizon correction arises from the modified Einstein equation:
        G_μν + γ(φ/c²) g_μν = (8πG_eff/c⁴) T_μν
    
    Parameters
    ----------
    k : array
        Wavenumber magnitude [1/Mpc]
    a : float
        Scale factor
    H_a : float
        Hubble parameter [1/s]
    H0 : float
        Hubble constant today [1/s]
    Omega_m : float
        Matter density parameter
    gamma : float
        HQIV thermodynamic coefficient
    box_size : float
        Box size [Mpc]
    N_grid : int
        Grid size
        
    Returns
    -------
    correction : array
        Horizon correction term (additive to Poisson RHS)
        
    Notes
    -----
    The horizon term modifies the Poisson equation on large scales
    (small k). It provides an effective density contribution that
    mimics some effects of dark energy.
    """
    # Horizon wavenumber
    k_horizon = H_a / c * Mpc_m  # [1/Mpc]
    
    # Smooth transition: no correction on small scales (large k)
    # Full correction on large scales (small k)
    scale_factor = 1.0 / (1.0 + (k / k_horizon)**2)
    
    # Magnitude from horizon term
    # The horizon term contributes an effective density proportional to γ H²
    H_ratio = H_a / H0
    
    # Effective density from horizon
    rho_horizon_eff = gamma * H_ratio**2
    
    # Correction term (in Fourier space, additive to Poisson source)
    correction = rho_horizon_eff * scale_factor * a**2
    
    return correction


def modified_poisson_rhs(density, a, H_a, H0, Omega_m, gamma=_gamma, 
                          alpha_G=_alpha_G, dx=1.0, include_horizon=True):
    """
    Compute the modified right-hand side of the Poisson equation.
    
    Paper Reference: paper/main.tex, Eq. 11
        ∇²Φ = 4πG_eff(φ) ρ_m δ + horizon_correction
    
    In code units (normalized to critical density):
        ∇²Φ = (3/2) Ω_m G_eff a δ + horizon_term
    
    Parameters
    ----------
    density : array
        Density contrast δ [N, N, N]
    a : float
        Scale factor
    H_a : float
        Hubble parameter [1/s]
    H0 : float
        Hubble constant today [1/s]
    Omega_m : float
        Matter density parameter
    gamma : float
        HQIV thermodynamic coefficient
    alpha_G : float
        Exponent for varying G
    dx : float
        Grid spacing [Mpc]
    include_horizon : bool
        Whether to include horizon correction
        
    Returns
    -------
    rhs : array
        Modified RHS of Poisson equation
    G_eff : float
        Effective gravitational constant (in units of G0)
    """
    N = density.shape[0]
    
    # Effective gravitational constant
    G_eff = effective_gravitational_constant(H_a, H0, alpha_G)
    
    # Standard Poisson RHS with varying G
    # In our units: ∇²Φ = (3/2) Ω_m G_eff a δ
    rhs = 1.5 * Omega_m * G_eff * a * density
    
    if include_horizon:
        # Add horizon correction in Fourier space
        rhs_k = fftn(rhs)
        
        # k-vectors
        kx = fftfreq(N, d=dx) * 2 * np.pi
        ky = fftfreq(N, d=dx) * 2 * np.pi
        kz = fftfreq(N, d=dx) * 2 * np.pi
        KX, KY, KZ = np.meshgrid(kx, ky, kz, indexing='ij')
        k_mag = np.sqrt(KX**2 + KY**2 + KZ**2)
        
        # Horizon correction
        horizon_corr_k = horizon_correction_term(
            k_mag, a, H_a, H0, Omega_m, gamma
        )
        
        # Add to RHS (the correction is already in Fourier space form)
        # We add it as an effective density contribution
        rhs_k += horizon_corr_k * np.prod(rhs.shape)  # Normalize
        
        rhs = np.real(ifftn(rhs_k))
    
    return rhs, G_eff


class EffectiveGravitySolver:
    """
    Class to manage effective gravity computation throughout the simulation.
    
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
    N_grid : int
        Grid size
    box_size : float
        Box size [Mpc]
    """
    
    def __init__(self, H0=_H0, Omega_m=_Omega_m, gamma=_gamma, alpha_G=_alpha_G,
                 N_grid=256, box_size=100.0):
        self.H0 = H0 * 1000.0 / Mpc_m  # Convert to 1/s
        self.H0_km = H0  # Keep in km/s/Mpc for reference
        self.Omega_m = Omega_m
        self.gamma = gamma
        self.alpha_G = alpha_G
        self.N_grid = N_grid
        self.box_size = box_size
        self.dx = box_size / N_grid
        
        # Precompute k-grid
        self._init_kgrid()
        
        # Storage
        self.G_eff = None
        self.horizon_correction = None
        
    def _init_kgrid(self):
        """Initialize k-space grid."""
        kx = fftfreq(self.N_grid, d=self.dx) * 2 * np.pi
        ky = fftfreq(self.N_grid, d=self.dx) * 2 * np.pi
        kz = fftfreq(self.N_grid, d=self.dx) * 2 * np.pi
        
        self.KX, self.KY, self.KZ = np.meshgrid(kx, ky, kz, indexing='ij')
        self.k_mag = np.sqrt(self.KX**2 + self.KY**2 + self.KZ**2)
        self.k_mag[0, 0, 0] = 1.0  # Avoid division by zero
        
    def compute_G_eff(self, H_a):
        """Compute effective gravitational constant."""
        self.G_eff = effective_gravitational_constant(H_a, self.H0, self.alpha_G)
        return self.G_eff
    
    def compute_rhs(self, density, a, H_a, include_horizon=True):
        """
        Compute modified Poisson RHS.
        
        Parameters
        ----------
        density : array
            Density contrast
        a : float
            Scale factor
        H_a : float
            Hubble parameter
        include_horizon : bool
            Include horizon correction
            
        Returns
        -------
        rhs : array
            Modified RHS for Poisson equation
        """
        # Compute G_eff
        self.compute_G_eff(H_a)
        
        # Standard Poisson RHS
        rhs = 1.5 * self.Omega_m * self.G_eff * a * density
        
        if include_horizon:
            # Horizon correction in Fourier space
            k_horizon = H_a / c * Mpc_m
            scale_factor = 1.0 / (1.0 + (self.k_mag / k_horizon)**2)
            
            H_ratio = H_a / self.H0
            horizon_eff = self.gamma * H_ratio**2 * scale_factor * a**2
            
            # Add horizon contribution
            rhs_k = fftn(rhs)
            rhs_k += horizon_eff * self.N_grid**3
            rhs = np.real(ifftn(rhs_k))
            
            self.horizon_correction = horizon_eff
        
        return rhs
    
    def solve_poisson_fft(self, rhs):
        """
        Solve Poisson equation using FFT.
        
        ∇²Φ = rhs
        
        In Fourier space:
        Φ_k = -rhs_k / k²
        """
        rhs_k = fftn(rhs)
        
        # Avoid division by zero at k=0
        phi_k = -rhs_k / (self.k_mag**2)
        phi_k[0, 0, 0] = 0
        
        phi = np.real(ifftn(phi_k))
        
        return phi
    
    def compute_acceleration(self, phi):
        """
        Compute gravitational acceleration from potential.
        
        g = -∇Φ
        
        In Fourier space:
        g_k = -i k Φ_k
        """
        phi_k = fftn(phi)
        
        # Acceleration components
        gx_k = -1j * self.KX * phi_k
        gy_k = -1j * self.KY * phi_k
        gz_k = -1j * self.KZ * phi_k
        
        gx = np.real(ifftn(gx_k))
        gy = np.real(ifftn(gy_k))
        gz = np.real(ifftn(gz_k))
        
        return np.stack([gx, gy, gz])
    
    def full_poisson_solve(self, density, a, H_a, include_horizon=True):
        """
        Full Poisson solve: RHS → potential → acceleration.
        
        Parameters
        ----------
        density : array
            Density contrast
        a : float
            Scale factor
        H_a : float
            Hubble parameter
        include_horizon : bool
            Include horizon correction
            
        Returns
        -------
        phi : array
            Gravitational potential
        acceleration : array
            Gravitational acceleration [3, N, N, N]
        G_eff : float
            Effective G
        """
        # Compute modified RHS
        rhs = self.compute_rhs(density, a, H_a, include_horizon)
        
        # Solve Poisson equation
        phi = self.solve_poisson_fft(rhs)
        
        # Compute acceleration
        acceleration = self.compute_acceleration(phi)
        
        return phi, acceleration, self.G_eff


# =============================================================================
# Testing
# =============================================================================

def test_g_eff():
    """Test effective gravitational constant computation."""
    import matplotlib.pyplot as plt
    
    print("Testing G_eff computation...")
    
    # Test G_eff vs redshift
    H0 = 73.2 * 1000.0 / Mpc_m  # 1/s
    
    z = np.linspace(0, 10, 1000)
    a = 1.0 / (1.0 + z)
    
    # Simple H(a) model (matter-dominated for illustration)
    # In reality, use full HQIV H(a)
    H_a = H0 * np.sqrt(0.06 * a**(-3) + 0.94)  # Approximate
    
    # Compute G_eff
    G_eff = effective_gravitational_constant(H_a, H0, alpha=0.6)
    
    print(f"G_eff today (z=0): {G_eff[0]:.4f} G0")
    print(f"G_eff at z=5: {G_eff[np.argmin(np.abs(z-5))]:.4f} G0")
    print(f"G_eff at z=10: {G_eff[-1]:.4f} G0")
    
    # Test horizon correction
    N = 64
    box_size = 100.0
    dx = box_size / N
    
    solver = EffectiveGravitySolver(H0=73.2, Omega_m=0.06, gamma=0.40,
                                     N_grid=N, box_size=box_size)
    
    # Test density field
    density = np.zeros((N, N, N))
    density[N//2, N//2, N//2] = 1.0  # Point mass
    
    a_test = 0.5
    H_test = H0 * np.sqrt(0.06 * a_test**(-3) + 0.94)
    
    phi, acc, G_eff_val = solver.full_poisson_solve(density, a_test, H_test)
    
    print(f"\nPoisson solve test:")
    print(f"  G_eff = {G_eff_val:.4f} G0")
    print(f"  Phi max = {phi.max():.4e}")
    print(f"  Acc max = {np.max(np.sqrt(np.sum(acc**2, axis=0))):.4e}")
    
    # Plot
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    
    ax = axes[0]
    ax.plot(z, G_eff, 'b-', linewidth=2)
    ax.axhline(1.0, color='k', linestyle='--', alpha=0.5)
    ax.set_xlabel('Redshift z')
    ax.set_ylabel('G_eff / G0')
    ax.set_title('Varying Gravitational Coupling')
    ax.grid(True, alpha=0.3)
    
    ax = axes[1]
    slice_idx = N//2
    im = ax.imshow(phi[:, :, slice_idx], origin='lower', 
                   extent=[0, box_size, 0, box_size])
    ax.set_title('Potential Φ (slice)')
    ax.set_xlabel('x [Mpc]')
    ax.set_ylabel('y [Mpc]')
    plt.colorbar(im, ax=ax)
    
    plt.tight_layout()
    plt.savefig('g_eff_test.png', dpi=150)
    plt.close()
    print("Saved g_eff_test.png")
    
    return G_eff


if __name__ == "__main__":
    test_g_eff()