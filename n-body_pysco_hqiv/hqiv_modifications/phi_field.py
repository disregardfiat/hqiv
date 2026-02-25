"""
HQIV Phi Field Module
=====================

Computes the geometric horizon field φ(x) = 2c²/Θ_local(x) from the action principle.

Paper Reference: paper/main.tex, Section 4
    φ(x) ≡ 2c² / Θ_local(x)

where Θ_local(x) is the local causal-horizon radius measured by fundamental observers.
In homogeneous FLRW background: Θ = 2c/H, so φ = cH.

The field can also be computed from the expansion scalar θ = ∇·v of the timelike
congruence of fundamental observers.

Author: HQIV Team
"""

import numpy as np
from numpy.fft import fftn, ifftn, fftfreq
from scipy.interpolate import interp1d
import configparser
from pathlib import Path

# Physical constants (SI units)
c = 2.99792458e8       # Speed of light [m/s]
Mpc_m = 3.0856775814913673e22  # Megaparsec in meters


def load_config():
    """Load configuration from config_hqiv.ini."""
    config_path = Path(__file__).parent / 'config_hqiv.ini'
    config = configparser.ConfigParser()
    config.read(config_path)
    return config


def compute_theta_local(velocity_field, dx, nthreads=4):
    """
    Compute the local expansion scalar θ = ∇·v on the PM grid.
    
    This is the divergence of the peculiar velocity field in comoving coordinates.
    
    Parameters
    ----------
    velocity_field : array
        Velocity field [3, N, N, N] in comoving coordinates [km/s or code units]
    dx : float
        Grid spacing in comoving coordinates [Mpc]
    nthreads : int
        Number of threads for FFT
        
    Returns
    -------
    theta : array
        Expansion scalar θ = ∇·v [1/Mpc or matching velocity units]
        
    Notes
    -----
    In Fourier space: θ_k = i k · v_k
    """
    # Get grid dimensions
    N = velocity_field.shape[1]
    
    # Compute k-vectors
    kx = fftfreq(N, d=dx) * 2 * np.pi
    ky = fftfreq(N, d=dx) * 2 * np.pi
    kz = fftfreq(N, d=dx) * 2 * np.pi
    
    KX, KY, KZ = np.meshgrid(kx, ky, kz, indexing='ij')
    
    # Transform velocities to Fourier space
    vx_k = fftn(velocity_field[0])
    vy_k = fftn(velocity_field[1])
    vz_k = fftn(velocity_field[2])
    
    # Divergence in Fourier space: θ_k = i k · v_k
    theta_k = 1j * (KX * vx_k + KY * vy_k + KZ * vz_k)
    
    # Transform back to real space
    theta = np.real(ifftn(theta_k))
    
    return theta


def compute_Theta_local(theta, H_a, a, method='expansion_scalar'):
    """
    Compute the local causal-horizon radius Θ_local.
    
    Paper Reference: paper/main.tex, Section 2
        Θ_local(x) = proper distance along past light cone to nearest null surface
    
    In practice, we use approximations based on available quantities.
    
    Parameters
    ----------
    theta : array
        Local expansion scalar θ = ∇·v
    H_a : float
        Hubble parameter at scale factor a [1/s]
    a : float
        Scale factor
    method : str
        Method to compute Θ_local:
        - 'expansion_scalar': Use θ directly (Θ ~ 3c/θ for local expansion)
        - 'background': Use FLRW background value Θ = 2c/H
        - 'hybrid': Combine background with local perturbations
        
    Returns
    -------
    Theta_local : array
        Local horizon radius [meters]
        
    Notes
    -----
    In homogeneous FLRW: Θ = 2c/H
    With perturbations: Θ_local varies with local expansion rate
    """
    # Background horizon radius
    Theta_background = 2 * c / H_a
    
    if method == 'background':
        # Use background value everywhere
        return np.full_like(theta, Theta_background, dtype=np.float64)
    
    elif method == 'expansion_scalar':
        # Use local expansion scalar
        # θ = 3H in background, so Θ_local ~ 3c/|θ|
        # Add small regularization to avoid division by zero
        theta_abs = np.abs(theta) + 1e-100
        
        # Convert theta to physical units if needed
        # theta is in comoving units, need to convert
        # θ_physical = θ_comoving / a (approximately)
        theta_physical = theta_abs / a
        
        # Horizon from expansion: Θ ~ 3c/θ_physical
        Theta_local = 3 * c / theta_physical
        
        # Clamp to reasonable range
        Theta_min = 0.1 * Theta_background
        Theta_max = 10.0 * Theta_background
        Theta_local = np.clip(Theta_local, Theta_min, Theta_max)
        
        return Theta_local
    
    elif method == 'hybrid':
        # Combine background with local perturbations
        # Θ_local = Θ_background * (1 + δΘ/Θ)
        # where δΘ/Θ ~ -δθ/θ (perturbation in expansion)
        
        theta_background = 3 * H_a  # Background expansion scalar
        
        # Relative perturbation in expansion
        delta_theta_rel = (theta - theta_background) / (theta_background + 1e-100)
        
        # Horizon perturbation (opposite sign)
        delta_Theta_rel = -delta_theta_rel
        
        # Clamp perturbation amplitude
        delta_Theta_rel = np.clip(delta_Theta_rel, -0.5, 0.5)
        
        # Local horizon
        Theta_local = Theta_background * (1 + delta_Theta_rel)
        
        return Theta_local
    
    else:
        raise ValueError(f"Unknown method: {method}")


def phi_from_expansion_scalar(theta, H_a, a, method='hybrid'):
    """
    Compute φ directly from expansion scalar.
    
    Paper Reference: paper/main.tex, Section 4
        φ(x) ≡ 2c² / Θ_local(x)
        
    In FLRW background: φ = cH
    
    Parameters
    ----------
    theta : array
        Local expansion scalar θ = ∇·v
    H_a : float
        Hubble parameter at scale factor a [1/s]
    a : float
        Scale factor
    method : str
        Method to compute Θ_local
        
    Returns
    -------
    phi : array
        Horizon field φ [m/s²] (units of acceleration)
        
    Notes
    -----
    The field φ has units of acceleration (m/s²) because:
    - Θ has units of length (m)
    - φ = 2c²/Θ has units of m²/s² / m = m/s²
    """
    Theta_local = compute_Theta_local(theta, H_a, a, method=method)
    phi = 2 * c**2 / Theta_local
    
    return phi


def compute_phi_field(density, velocity_field, H_a, a, dx, 
                      method='hybrid', use_density_fallback=True):
    """
    Compute the full φ field on the PM grid.
    
    This is the main function to compute the horizon field for the simulation.
    
    Parameters
    ----------
    density : array
        Density contrast δ on the grid
    velocity_field : array
        Velocity field [3, N, N, N]
    H_a : float
        Hubble parameter [1/s]
    a : float
        Scale factor
    dx : float
        Grid spacing [Mpc]
    method : str
        Method for computing Θ_local
    use_density_fallback : bool
        If True, use density contrast as fallback when velocity field is unavailable
        
    Returns
    -------
    phi : array
        Horizon field φ [m/s²]
    theta : array
        Expansion scalar θ [1/Mpc]
        
    Notes
    -----
    Paper equation: φ(x) = 2c²/Θ_local(x)
    
    In the low-acceleration regime (large Θ), φ is small.
    In the high-acceleration regime (small Θ), φ is large.
    """
    # Compute expansion scalar from velocity field
    if velocity_field is not None:
        theta = compute_theta_local(velocity_field, dx)
    elif use_density_fallback:
        # Fallback: estimate θ from density using linear theory
        # In linear theory: θ ~ -f H δ, where f is growth rate
        # For HQIV, f is modified, but we use approximation
        f_growth = 1.0  # Approximate growth rate d ln D / d ln a
        H_Mpc = H_a * Mpc_m  # Convert to 1/Mpc
        theta = -f_growth * H_Mpc * density
    else:
        raise ValueError("No velocity field and density fallback disabled")
    
    # Compute φ from expansion scalar
    phi = phi_from_expansion_scalar(theta, H_a, a, method=method)
    
    return phi, theta


def compute_phi_gradient(phi, dx):
    """
    Compute the spatial gradient ∇φ of the horizon field.
    
    This gradient is crucial for the vorticity source term.
    
    Parameters
    ----------
    phi : array
        Horizon field φ [N, N, N]
    dx : float
        Grid spacing [Mpc]
        
    Returns
    -------
    grad_phi : array
        Gradient of φ [3, N, N, N] in units of [φ_units / Mpc]
        
    Notes
    -----
    Paper Reference: paper/main.tex, Section 6
        The vorticity source is proportional to (k × ∇φ)
    """
    N = phi.shape[0]
    
    # k-vectors for Fourier differentiation
    kx = fftfreq(N, d=dx) * 2 * np.pi
    ky = fftfreq(N, d=dx) * 2 * np.pi
    kz = fftfreq(N, d=dx) * 2 * np.pi
    
    KX, KY, KZ = np.meshgrid(kx, ky, kz, indexing='ij')
    
    # Transform to Fourier space
    phi_k = fftn(phi)
    
    # Gradient in Fourier space: ∇φ_k = i k φ_k
    grad_phi_x = np.real(ifftn(1j * KX * phi_k))
    grad_phi_y = np.real(ifftn(1j * KY * phi_k))
    grad_phi_z = np.real(ifftn(1j * KZ * phi_k))
    
    return np.stack([grad_phi_x, grad_phi_y, grad_phi_z])


class PhiFieldSolver:
    """
    Class to manage the φ field computation throughout the simulation.
    
    This class caches the background H(a) and provides efficient methods
    to compute φ at each timestep.
    
    Parameters
    ----------
    H0 : float
        Hubble constant today [km/s/Mpc]
    Omega_m : float
        Matter density parameter
    gamma : float
        HQIV thermodynamic coefficient
    N_grid : int
        Grid size
    box_size : float
        Box size [Mpc]
    """
    
    def __init__(self, H0=73.2, Omega_m=0.06, gamma=0.40, 
                 N_grid=256, box_size=100.0):
        self.H0 = H0 * 1000.0 / Mpc_m  # Convert to 1/s
        self.Omega_m = Omega_m
        self.gamma = gamma
        self.N_grid = N_grid
        self.box_size = box_size
        self.dx = box_size / N_grid
        
        # Precompute background H(a) table
        self._precompute_background()
        
        # Storage for current phi field
        self.phi = None
        self.theta = None
        self.grad_phi = None
        
    def _precompute_background(self, a_min=1e-4, a_max=1.0, n_pts=1000):
        """Precompute background H(a) for interpolation."""
        from scipy.integrate import quad
        
        self.a_bg = np.logspace(np.log10(a_min), np.log10(a_max), n_pts)
        self.H_bg = np.zeros(n_pts)
        
        # Compute radiation density
        h = self.H0 * Mpc_m / 1000.0 / 100.0  # h = H0 / 100 km/s/Mpc
        Omega_r = 2.47e-5 / h**2 * (1 + 0.2271 * 3.046)
        
        for i, a in enumerate(self.a_bg):
            # HQIV modified Friedmann: 3H² - γH = 8πG ρ
            # Solve for H: H = (γ + sqrt(γ² + 24 Ω(a) H0²)) / 6
            rho_norm = self.Omega_m * a**(-3) + Omega_r * a**(-4)
            disc = self.gamma**2 + 24.0 * rho_norm
            self.H_bg[i] = (self.gamma + np.sqrt(disc)) / 6.0 * self.H0
        
        # Create interpolation function
        self.H_interp = interp1d(self.a_bg, self.H_bg, kind='cubic', 
                                  fill_value='extrapolate')
    
    def H_of_a(self, a):
        """Get Hubble parameter at scale factor a."""
        return float(self.H_interp(a))
    
    def update(self, density, velocity_field, a, method='hybrid'):
        """
        Update the φ field for the current timestep.
        
        Parameters
        ----------
        density : array
            Density contrast δ
        velocity_field : array
            Velocity field [3, N, N, N]
        a : float
            Scale factor
        method : str
            Method for computing Θ_local
        """
        H_a = self.H_of_a(a)
        
        self.phi, self.theta = compute_phi_field(
            density, velocity_field, H_a, a, self.dx, method=method
        )
        
        self.grad_phi = compute_phi_gradient(self.phi, self.dx)
        
        return self.phi
    
    def get_phi_at_positions(self, positions):
        """
        Interpolate φ to particle positions using CIC.
        
        Parameters
        ----------
        positions : array
            Particle positions [N_part, 3] in comoving coordinates [Mpc]
            
        Returns
        -------
        phi_particles : array
            φ values at particle positions
        """
        if self.phi is None:
            raise RuntimeError("Call update() first")
        
        # CIC interpolation
        phi_p = np.zeros(len(positions))
        
        for i, pos in enumerate(positions):
            ix = int(np.floor(pos[0] / self.dx)) % self.N_grid
            iy = int(np.floor(pos[1] / self.dx)) % self.N_grid
            iz = int(np.floor(pos[2] / self.dx)) % self.N_grid
            
            ux = pos[0] / self.dx - np.floor(pos[0] / self.dx)
            uy = pos[1] / self.dx - np.floor(pos[1] / self.dx)
            uz = pos[2] / self.dx - np.floor(pos[2] / self.dx)
            
            for dx in [0, 1]:
                for dy in [0, 1]:
                    for dz in [0, 1]:
                        w = ((1-ux) if dx==0 else ux) * \
                            ((1-uy) if dy==0 else uy) * \
                            ((1-uz) if dz==0 else uz)
                        phi_p[i] += w * self.phi[
                            (ix+dx) % self.N_grid,
                            (iy+dy) % self.N_grid,
                            (iz+dz) % self.N_grid
                        ]
        
        return phi_p


# =============================================================================
# Convenience functions for testing
# =============================================================================

def test_phi_field():
    """Test the φ field computation with a simple density field."""
    import matplotlib.pyplot as plt
    
    print("Testing φ field computation...")
    
    # Create a simple test case
    N = 64
    box_size = 100.0  # Mpc
    dx = box_size / N
    
    # Create a density field with a single overdensity
    density = np.zeros((N, N, N))
    x = np.linspace(0, box_size, N, endpoint=False)
    X, Y, Z = np.meshgrid(x, x, x, indexing='ij')
    
    # Gaussian overdensity at center
    center = box_size / 2
    r = np.sqrt((X - center)**2 + (Y - center)**2 + (Z - center)**2)
    density = 5.0 * np.exp(-r**2 / (10.0)**2)
    
    # Create a simple velocity field (infall toward center)
    velocity = np.zeros((3, N, N, N))
    velocity[0] = -(X - center) / (r + 1e-10) * density * 0.1
    velocity[1] = -(Y - center) / (r + 1e-10) * density * 0.1
    velocity[2] = -(Z - center) / (r + 1e-10) * density * 0.1
    
    # Compute φ field
    H_a = 73.2 * 1000.0 / Mpc_m  # H0 in 1/s
    a = 0.5
    
    phi, theta = compute_phi_field(density, velocity, H_a, a, dx)
    
    print(f"φ range: {phi.min():.2e} to {phi.max():.2e} m/s²")
    print(f"θ range: {theta.min():.2e} to {theta.max():.2e} 1/Mpc")
    
    # Plot slice
    fig, axes = plt.subplots(1, 3, figsize=(15, 4))
    
    ax = axes[0]
    im = ax.imshow(density[:, :, N//2], origin='lower', extent=[0, box_size, 0, box_size])
    ax.set_title('Density δ')
    ax.set_xlabel('x [Mpc]')
    ax.set_ylabel('y [Mpc]')
    plt.colorbar(im, ax=ax)
    
    ax = axes[1]
    im = ax.imshow(theta[:, :, N//2], origin='lower', extent=[0, box_size, 0, box_size])
    ax.set_title('Expansion scalar θ')
    ax.set_xlabel('x [Mpc]')
    ax.set_ylabel('y [Mpc]')
    plt.colorbar(im, ax=ax)
    
    ax = axes[2]
    im = ax.imshow(phi[:, :, N//2], origin='lower', extent=[0, box_size, 0, box_size])
    ax.set_title('Horizon field φ [m/s²]')
    ax.set_xlabel('x [Mpc]')
    ax.set_ylabel('y [Mpc]')
    plt.colorbar(im, ax=ax)
    
    plt.tight_layout()
    plt.savefig('phi_field_test.png', dpi=150)
    plt.close()
    print("Saved phi_field_test.png")
    
    return phi, theta


if __name__ == "__main__":
    test_phi_field()